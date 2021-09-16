module saf_files

  use utilities

  implicit none
  !
  ! This module provides routines for accessing SAF file headers and data
  !
  ! File format tips: The header line appear to be terminated by 0D 0A ... in the sample
  ! file which is a CR LF (13, 10) 

  integer :: fileunit

  integer :: current_frame
  integer :: sign_type,var_type,var_len
  integer :: header_lines
  integer :: numdps

  character*256 :: saf_filename

    integer(kind=2),allocatable :: saf_data_int2(:,:,:)
    integer(kind=4),allocatable :: saf_data_int4(:,:,:)
    real, allocatable :: saf_data_flt(:,:,:)

    save

contains


  subroutine init_saf_file(filename,ierr)
    implicit none
    character*(*) filename
    integer :: form,ierr

    saf_filename = filename

    ! test that file can be opened

    ! use formatted for test
    form = 1
    
    call open_saf_file(form,ierr)

    ! close fileunit after test

    call close_saf_file(ierr)

  end subroutine init_saf_file




  subroutine open_saf_file(form,ierr)
    implicit none
    integer :: ierr
    integer :: form

    call find_free_unit_number(fileunit)

    if (form.eq.1) then 
       ! formatted
       open(fileunit,file=trim(saf_filename),form='formatted',iostat=ierr)
    elseif (form.eq.2) then 
       ! unformatted
       open(fileunit,file=trim(saf_filename),form='binary',iostat=ierr)
    endif

    if (ierr.ne.0) then 
       write(0,'(A,a,a,i6)') 'Error opening file:',trim(saf_filename),':Error =',ierr
    endif

    rewind (fileunit)

    return

  end subroutine open_saf_file

  subroutine close_saf_file(ierr)
    implicit none
    integer :: ierr
    integer :: form


    close(fileunit,iostat=ierr)

    if (ierr.ne.0) then 
       write(0,'(A,a,a,i6)') 'Error closing file:',trim(saf_filename),':Error =',ierr
    endif

    return

  end subroutine close_saf_file



  subroutine read_saf_file_header(xpixel,ypixel,nframes,data_type,data_len,reference_string,ierr)
    implicit none

    integer ::  xpixel,ypixel, nframes, data_type,data_len,ierr
    character*(*) :: reference_string

    logical :: data
    character*256 :: input_line
    integer :: ioerr,line_cnt
    integer :: xpixls,ypixls
    character*10 :: datype
    character*10 :: bytord
    character*10 :: keywrd
    integer :: imsize
    integer :: st_index

    ! Initialization
    imsize = 0
    xpixls = 0
    ypixls = 0
    numdps = 0

    current_frame = 0


    datype = ' '
    bytord = ' '
    keywrd = ' '

    ! Init Output values

    xpixel = 0 
    ypixel = 0
    nframes = 0
    data_type = 0
    data_len = 0
    reference_string = ' '

    ! start read

    ! open the file for formatted access 

    call open_saf_file(1,ierr)

    ierr = 0
    data = .false.
    sign_type = 0
    var_type  = 0
    var_len   = 0


    line_cnt = 0

    do while (.not.data)

       read(fileunit,'(a256)',iostat=ioerr) input_line

       line_cnt = line_cnt + 1

       !write(6,*) 'IN:'//trim(input_line)//':'

       if (ioerr.eq.0) then 

          if (line_cnt.eq.1.and.upcase(input_line(1:6)).ne.'HDSIZE') then 
             ! File is not in SAF format
             ierr = -1
             return
          elseif (upcase(input_line(1:6)).eq.'XPIXLS') then
             ! X pixels

             read(input_line(7:),*) xpixls

          elseif (upcase(input_line(1:6)).eq.'YPIXLS') then
             ! Y pixels
             read(input_line(7:),*) ypixls

          elseif (upcase(input_line(1:6)).eq.'DATYPE') then
             ! Type of data stored 
             ! Options are:
             ! Int8
             ! Int16
             ! Int32
             ! Int64
             ! Flt32
             ! Flt64
             ! RGB24  (3 bytes/pixel RGB)
             ! Can also have a U prepended for Unsigned 
             ! e.g. UInt16 = Unsigned 16 bit integer
             ! Number of bytes/pixel can also be determined by IMSIZE/(XPIXLS*YPIXLS)
             read(input_line(7:),*) datype

             ! Determine signed or unsigned
             if (datype(1:1).eq.'U'.or.datype(1:1).eq.'u') then 
                ! unsigned
                sign_type = 1
                st_index  = 2
             else 
                ! signed 
                sign_type = 2
                st_index  = 1
             endif
             ! Determine variable type - int, flt, rgb

             if (upcase(datype(st_index:st_index+2)).eq.'INT') then
                var_type = 1
             elseif (upcase(datype(st_index:st_index+2)).eq.'FLT') then
                var_type = 2
             elseif (upcase(datype(st_index:st_index+2)).eq.'RGB') then
                var_type = 3
             endif

             ! Determine variable length
             read(datype(st_index+3:),*) var_len


          elseif (upcase(input_line(1:6)).eq.'BYTORD') then
             ! Byte ordering ... usually Little endian appears to be used
             ! LH = Low/High = Little endian (typical PC) 
             ! HL = High/Low = Big endian 
             ! VX = Vax = byte order used on a vax
             read(input_line(7:),*) bytord

          elseif (upcase(input_line(1:6)).eq.'IMSIZE') then
             ! Number of bytes in the image
             read(input_line(7:),*) imsize


          elseif (upcase(input_line(1:6)).eq.'KEYWRD') then
             ! Type if data in the safe file 
             ! IMG = image file
             ! CMAP= colour map indexed image file 
             read(input_line(7:),*) keywrd


          elseif (upcase(input_line(1:6)).eq.'NUMDPS') then
             ! Number of data points ... in this case I think images or pictures in the file
             read(input_line(7:),*) numdps

          elseif (upcase(input_line(1:4)).eq.'DATA') then
             data = .true.
             header_lines = line_cnt
          endif

       else
          ! Io error on input
          ierr = ioerr
          return
       endif


    end do


    !
    ! Verify data read in and assign information to input arguments
    !

    write(6,'(a)') 'READ_SAF_HEADER: DATA READ IN:'
    write(6,'(a,i8)') 'XPIXELS     : ',xpixls
    write(6,'(a,i8)') 'YPIXELS     : ',ypixls
    write(6,'(a,i8)') 'NUMDPS      : ',numdps
    write(6,'(a,i8)') 'IMSIZE      : ',imsize
    write(6,'(a,g12.5)') 'BYTES/PIXEL : ',imsize/(xpixls*ypixls)
    write(6,'(a,a)')  'IMAGE TYPE:',trim(keywrd)
    write(6,'(a,a)')  'BYTE ORDER:',trim(bytord)
    write(6,'(a,a)')  'DATA TYPE :',trim(datype)
    write(6,'(a,i8)') ' -SIGNED  :',sign_type
    write(6,'(a,i8)') ' -TYPE    :',var_type
    write(6,'(a,i8)') ' -LEN-BITS:',var_len
    write(6,'(a,i8)') ' -LEN-BYTE:',var_len/8

    if (xpixls.le.0.or.ypixls.le.0.or.numdps.le.0.or.imsize.le.0) then 
       write(6,'(a)') 'ERROR IM IMAGE DATA: XPIXLS, YPIXLS, NUMDPS OR IMSIZE INCONSISTENT'
       ierr = -2
       return
    endif


    ! Note that although the data in the file may be UNSIGNED - fortran doesn't have an unsigned data type so 
    ! some effort must be made when reading in the data to translate it properly. 

    xpixel = xpixls
    ypixel = ypixls
    nframes = numdps
    data_type = var_type
    data_len = var_len/8

    call close_saf_file(ierr)

  end subroutine read_saf_file_header


  subroutine read_saf_data(xpixel,ypixel)
    implicit none
    integer :: xpixel,ypixel
    integer :: in,ix,iy
    integer :: ierr

    ! Decide which temporary storage to allocate
    if (var_type.eq.1) then 
       ! integer
       if (var_len.eq.16) then

           allocate(saf_data_int2(xpixel,ypixel,numdps))

       elseif (var_len.eq.32) then 

           allocate(saf_data_int4(xpixel,ypixel,numdps))

       endif
    elseif (var_type.eq.2) then 
       ! float
       allocate(saf_data_flt(xpixel,ypixel,numdps))

    elseif (var_type.eq.3) then 
       ! rgb not supported yet
       write(0,*) 'RGB data type not yet supported - STOP'
       stop
    endif

    call open_saf_file(2,ierr)

    ! scan for data string

    call remove_header

    ! load data

    do in = 1,numdps

       !write(0,*) 'reading frame :',in
       
       if (var_type.eq.1) then 

          if (var_len.eq.16) then 

             read(fileunit) ((saf_data_int2(ix,iy,in),ix=1,xpixel),iy=1,ypixel)

          elseif (var_len.eq.32) then 

             read(fileunit) ((saf_data_int4(ix,iy,in),ix=1,xpixel),iy=1,ypixel)

          endif

!          read(fileunit,inform) ((temp_int(ix,iy),ix=1,xpixel),iy=1,ypixel)

       elseif (var_type.eq.2) then 

          read(fileunit) ((saf_data_flt(ix,iy,in),ix=1,xpixel),iy=1,ypixel)

       endif


    end do

    call close_saf_file(ierr)

  end subroutine read_saf_data

  subroutine empty_saf_data
    implicit none

    if (allocated(saf_data_int2)) deallocate(saf_data_int2)
    if (allocated(saf_data_int4)) deallocate(saf_data_int4)
    if (allocated(saf_data_flt)) deallocate(saf_data_flt)

  end subroutine empty_saf_data


  
  subroutine get_saf_frame(frame,xpixel,ypixel,nframe)
    implicit none

    integer :: xpixel,ypixel,nframe
    real :: frame(xpixel,ypixel)

    character*4 :: form_2byte,form_4byte,inform
    integer :: ix,iy,in

    integer :: maxvalue_i,minvalue_i
    integer :: ierr

    ierr = 0
    
       !write(0,*) 'reading frame :',in
       
       if (var_type.eq.1) then 

          if (var_len.eq.16) then 

             frame = saf_data_int2(:,:,nframe)

          elseif (var_len.eq.32) then 

             frame = saf_data_int4(:,:,nframe)

          endif

       elseif (var_type.eq.2) then 

          frame = saf_data_flt(:,:,nframe)

       endif


  end subroutine get_saf_frame




  subroutine read_frame(frame,xpixel,ypixel,nframe)
    implicit none

    integer :: xpixel,ypixel,nframe
    real :: frame(xpixel,ypixel)

    integer(kind=2),allocatable :: temp_int2(:,:)
    integer(kind=4),allocatable :: temp_int4(:,:)
    real, allocatable :: temp_flt(:,:)
    
    character*4 :: form_2byte,form_4byte,inform
    integer :: ix,iy,in

    integer :: maxvalue_i,minvalue_i
    integer :: ierr

    ierr = 0

    call open_saf_file(2,ierr)


    ! Decide which temporary storage to allocate
    if (var_type.eq.1) then 
       ! integer
       if (var_len.eq.16) then

           allocate(temp_int2(xpixel,ypixel))

       elseif (var_len.eq.32) then 

           allocate(temp_int4(xpixel,ypixel))

       endif
    elseif (var_type.eq.2) then 
       ! float
       allocate(temp_flt(xpixel,ypixel))

    elseif (var_type.eq.3) then 
       ! rgb not supported yet
       write(0,*) 'RGB data type not yet supported - STOP'
       stop
    endif

    ! Set up binary format string
    !if (var_len.eq.16) then 
       ! 2 bytes
    !   inform = form_2byte
    !elseif (var_len.eq.32) then 
       ! 4 bytes
    !   inform = form_4byte
    !else
       ! data length not yet supported
    !   write(0,*) 'Data length:',var_len,' not yet supported - STOP'
    !   stop
    !endif


    ! scan for data string

    call remove_header


    ! scan for desired frame - it is the last one read - note that the following technique essentially assumes unsigned data
    
    
    do in = 1,nframe

       !write(0,*) 'reading frame :',in
       
       if (var_type.eq.1) then 

          if (var_len.eq.16) then 

             read(fileunit) ((temp_int2(ix,iy),ix=1,xpixel),iy=1,ypixel)

          elseif (var_len.eq.32) then 

             read(fileunit) ((temp_int4(ix,iy),ix=1,xpixel),iy=1,ypixel)

          endif

!          read(fileunit,inform) ((temp_int(ix,iy),ix=1,xpixel),iy=1,ypixel)

       elseif (var_type.eq.2) then 

          read(fileunit) ((temp_flt(ix,iy),ix=1,xpixel),iy=1,ypixel)

       endif


    end do


    ! Assign result to input frame

    if (var_type.eq.1) then 

       if (var_len.eq.16) then 
          
          maxvalue_i = maxval(temp_int2)
          minvalue_i = minval(temp_int2)

          if (minvalue_i.lt.0) then
             write(6,'(a,i10)') 'WARNING: MIN Value less than zero ... signed data needs to be rescaled:',minvalue_i
          endif

          !write(0,'(a,2i12)') 'MAX/MIN :', maxvalue_i, minvalue_i

          frame = temp_int2


       elseif (var_len.eq.32) then 
          
          frame = temp_int4

       endif

    elseif (var_type.eq.2) then
       
       frame = temp_flt

    endif


    if (allocated(temp_int2)) deallocate(temp_int2)
    if (allocated(temp_int4)) deallocate(temp_int4)
    if (allocated(temp_flt)) deallocate(temp_flt)

    call close_saf_file(ierr)


  end subroutine read_frame



    subroutine find_data
      implicit none
      logical :: data

      character*256 :: input_line
      integer :: ioerr,ierr
      

      data = .false.

      current_frame = 0

    do while (.not.data)

       read(fileunit,'(a256)',iostat=ioerr) input_line

       write(6,*) 'FD:',trim(input_line)

       if (ioerr.eq.0) then
          if (upcase(input_line(1:4)).eq.'DATA') then
             data = .true.
          endif

       else
          ! Io error on input
          ierr = ioerr

          write(0,*) 'ERROR encountered finding data line: STOP'
          stop
          !return
       endif


    end do



    end subroutine find_data

    subroutine remove_header
      implicit none
      integer :: in
      character*1 :: input_line
      integer :: line_cnt


      !write(6,*) 'RH:', header_lines

      line_cnt = 0

      do while (line_cnt.lt.header_lines)
         read(fileunit) input_line
         ! count lines by the line feed character - ascii 10 - lines should be terminated by either CR-LF for DOS files
         ! or just LF in Unix
         if (iachar(input_line).eq.10) then 
            line_cnt = line_cnt + 1
            !write(6,*) 'RH:CNT:',line_cnt
         endif
      end do


         
    end subroutine remove_header





end module saf_files
