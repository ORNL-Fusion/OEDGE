module analyse_lp
  use utilities

  private

  real,allocatable :: filedata(:,:,:)


  integer,parameter :: timebins=44
  real,parameter :: tstart = 2300
  real,parameter :: deltat = 50


  ! 5 bits of data for averaging/analysis ... jsat, Te, ne, Psin, deltaR, Power 
  integer,parameter :: ndata=6
  character*5 :: data_desc(ndata)

  integer :: nlp
  integer,allocatable :: nlp_count(:)
  integer,allocatable :: nlp_id(:)


  real :: times(timebins,2)

  !real :: base_profiles(timebins,nprobes,ndata)
  !real :: base_counts(timebins,nprobes,ndata)
  !real :: base_int_axis(timebins)

  real,allocatable :: base_profiles(:,:,:)
  real,allocatable :: base_counts(:,:,:)
  real :: base_int_axis(timebins)

  real :: deltabin(ndata)

  real,allocatable :: intensity_axis(:,:)
  real,allocatable :: profiles(:,:,:,:)
  real,allocatable :: counts(:,:,:,:)


  real :: maxrange(ndata),minrange(ndata)
  integer :: data_bins
  character*(100) :: fschan

  public :: setup_lp,accumulate_data, read_lp, analyse_print_data

  real :: mindata,maxdata


contains


  subroutine setup_lp(nbins)
    implicit none

    integer :: in,nbins

    integer :: ierr,it,ib




    if (allocated(intensity_axis)) deallocate(intensity_axis)
    if (allocated(profiles)) deallocate(profiles)
    if (allocated(counts)) deallocate(counts)


    ! allocate space for data accumulation
    allocate(intensity_axis(nbins,ndata),stat=ierr)
    allocate(profiles(nbins,timebins,ndata,nlp),stat=ierr)
    allocate(counts(nbins,timebins,ndata,nlp),stat=ierr)

    allocate(base_profiles(timebins,ndata,nlp),stat=ierr)
    allocate(base_counts(timebins,ndata,nlp),stat=ierr)



    data_bins = nbins

    ! set data descriptors

    data_desc(1) = 'Jsat'
    data_desc(2) = 'Te'
    data_desc(3) = 'ne'
    data_desc(4) = 'Psin'
    data_desc(5) = 'DeltaR'


    ! setup time sampling windows

    do it = 1,timebins
       if (it.eq.1) then 
          times(it,1) = tstart
          times(it,2) = times(it,1) + deltat
       else
          times(it,1) = times(it-1,2) 
          times(it,2) = times(it,1) + deltat 
       endif

       base_int_axis(it) = times(it,2) - deltat/2.0

    enddo


    ! set data ranges for binning

    do in = 1,ndata

       minrange(in) = 0.0
       maxrange(in) = 100.0

    end do

    ! binning is not needed for psi or deltar
    ! adjust maxrange for ne 
    maxrange(3) = maxrange(3) * 5e11


    ! set up binning information
    do in = 1,ndata
       deltabin(in) = (maxrange(in)-minrange(in))/nbins
       do ib = 1,data_bins
          intensity_axis(ib,in) = minrange(in) + (ib-0.5) * deltabin(in)
          !write(0,'(a,2i6,5(2x,g12.5))') in,ib,intensity_axis(ib,in),minrange(in),maxrange(in),deltabin(in)
       end do
    end do


    ! zero out data accumulation arrays

    profiles = 0.0
    counts = 0.0

    base_profiles = 0.0
    base_counts = 0.0

  end subroutine setup_lp

  subroutine accumulate_data
    implicit none

    integer :: nlines

    integer :: in,it,pn,id
    integer :: bin

    !write(0,*) 'Accumulate:',nlines

    do in = 1,nlp

       do id = 1,nlp_count(in)

          ! check chisq - set cut off arbitrarily at 0.05
          if (filedata(id,6,in) .lt. 0.05) then 

             do it = 1,timebins

                ! is data point in time window
                if (filedata(id,1,in).ge.times(it,1).and.filedata(id,1,in).le.times(it,2)) then 
                   ! record all 5 data points needed
                   !
                   ! jsat  = 2   ndata = 1
                   !

                   call record_data(filedata(id,2,in),1,it,in)
                   ! te    = 3  ndata = 2
                   call record_data(filedata(id,3,in),2,it,in)
                   ! ne    = 4  ndata = 3
                   call record_data(filedata(id,4,in),3,it,in)


                   ! psin and deltar need to be averaged but binning isn't required ... however, averaging should be easy still. 
                   ! psin  = 5  ndata = 4
                   call record_data(filedata(id,5,in),4,it,in)
                   ! deltar= 8  ndata = 5
                   call record_data(filedata(id,8,in),5,it,in)

                   !write(0,*) 'data:',in,it,filedata(in,4),filedata(in,5)
                   !write(0,*) filedata(in,8),it,bin

                end if
             end do
          end if

       end do

    end do

  end subroutine accumulate_data


  subroutine record_data(data_value,data_index,it,in)
    implicit none
    real :: data_value
    integer :: data_index,it,in
    integer :: bin

    bin = min(int(max((min(data_value,maxrange(data_index)) - minrange(data_index)),0.0) / deltabin(data_index) + 1),data_bins)
    profiles(bin,it,data_index,in) =  profiles(bin,it,data_index,in) + data_value
    counts(bin,it,data_index,in)   =  counts(bin,it,data_index,in)   + 1.0

  end subroutine record_data



  subroutine analyse_print_data(filenameout,printopt)
    implicit none
    character*(*) :: filenameout
    integer :: in,it,printopt,id,ip
    integer :: baserange,peaknum
    integer :: basepos(1)
    real :: base_count, base_intensity
    real, allocatable :: tmp_count(:)

    integer :: outunit,ios

    call find_free_unit_number(outunit)

    open(unit=outunit,file=trim(filenameout),iostat=ios)

    write(0,'(a,2x,i10)') 'Output file:'//trim(filenameout),data_bins

    ! calculate the average profiles

    do ip = 1,nlp
       do id = 1,ndata
          do in = 1,data_bins
             do it = 1,timebins

                if (counts(in,it,id,ip).ne.0.0) then 
                   profiles(in,it,id,ip) = profiles(in,it,id,ip) / counts(in,it,id,ip)
                else
                   profiles(in,it,id,ip) = 0.0
                endif

                !write(0,'(a,4i6,5(2x,g12.5))') 'Profiles:',in,it,id,ip,profiles(in,it,id,ip),counts(in,it,id,ip)
             end do
          end do

       end do
    end do

    ! write out the results


    !    if (printopt.ne.0) then 
    !       write(outunit,'(a,2(2x,g12.5))') 'Maximum and minimum range imposed for bins:',maxrange,minrange
    !       write(outunit,'(a,2(2x,g12.5))') 'Maximum and minimum range found in data   :',maxdata,mindata
    !    endif


    do it = 1,timebins


       do ip = 1,nlp

          do id = 1,ndata

             if (allocated(tmp_count)) deallocate(tmp_count)

             allocate(tmp_count(data_bins))
             tmp_count = counts(:,it,id,ip)

             basepos = maxloc(tmp_count)
             peaknum = basepos(1)

             !write(0,*) 'basepos=',peaknum

             baserange = 2

             base_intensity = 0.0
             base_count = 0.0

             if (id.eq.1.or.id.eq.2.or.id.eq.3) then 
                do in = max(peaknum-baserange,1),min(peaknum+baserange,data_bins)
                   base_intensity = base_intensity + profiles(in,it,id,ip) * counts(in,it,id,ip)
                   base_count = base_count + counts(in,it,id,ip)
                end do
             elseif (id.eq.4.or.id.eq.5) then 
                base_intensity = sum(profiles(:,it,id,ip)*counts(:,it,id,ip))
                base_count = sum(counts(:,it,id,ip))
             endif


             if (base_count.ne.0.0) then 
                base_intensity = base_intensity / base_count
             else
                base_intensity = 0.0
             endif

             base_profiles(it,id,ip) = base_intensity
             base_counts(it,id,ip) = base_count


          end do
       end do

       if (printopt.ne.0) then 

          write(outunit,'(a,i6,a,g12.5,a,g12.5,a)') 'Time bins (Set # ',it,' ) = ', times(it,1) ,':',times(it,2),' ms'



          do ip = 1,nlp
             do id = 1,ndata

                write(outunit,'(a,i6,a)')  'DATA for probe n=',ip,' Data type :'//data_desc(id)
                write(outunit,'(2(a,2x,g18.8))')  'Base Intensity Average = ',base_profiles(it,id,ip), ' Counts contributing                   = ', base_counts(it,id,ip)
                write(outunit,'(a)')    '  Intensity Axis    Counts      Average Intensity '
                do in = 1,data_bins
                   write(outunit,'(7(2x,g12.5))') intensity_axis(in,id),counts(in,it,id,ip),profiles(in,it,id,ip)
                end do

             end do
          end do

       endif

    end do


    ! Summary of baseline data
    write(outunit,'(a)') 'Summary of baseline LP data over time:'
    write(outunit,'(a,2x,g12.5)') 'Width of time windows (ms): ', deltat
    write(outunit,'(a)')    '   Time(ms)       Jsat           Te          ne             Psin             R              ProbeN'
    do ip = 1, nlp
       do it = 1,timebins
          ! only write out non-zero entries
          if (base_counts(it,1,ip).gt.0) then 
             write(outunit,'(7(2x,g12.5))') base_int_axis(it),base_profiles(it,1,ip), base_profiles(it,2,ip),base_profiles(it,3,ip),base_profiles(it,4,ip),base_profiles(it,5,ip),nlp_id(ip)
          endif 
       end do
       write(6,*)

    end do

    if (allocated(filedata)) deallocate(filedata)


  end subroutine analyse_print_data




  subroutine read_lp(filename,ios)
    implicit none

    integer :: fileunit
    integer :: ios,ierr
    integer :: nprobes,max_line_count

    character*(*) :: filename

    integer :: line_cnt,in


    write(0,'(a,i7,4a)') 'Read file:'//trim(filename)

    call find_free_unit_number(fileunit)

    open(unit=fileunit,file=filename,status='old',iostat=ios)

    if (ios.ne.0) then 
       write(0,*) 'FILE not found:',ios
       return
    endif


    if (allocated(filedata)) deallocate(filedata)
    if (allocated(nlp_count)) deallocate(nlp_count)
    if (allocated(nlp_id)) deallocate(nlp_id)

    ! prescan for number of lines in the file
    ! find data block of interest

    call read_lp_info(nprobes,max_line_count,fileunit)
    
    ! 8 data items/line
    ! time(msec)	   jsat(A/cm2)	      temp(eV)	  dens(cm-3)	        psin	   csq	              probeID     delrsepout2
    !  2401.36	      11.3374	      8.00457	  2.55826e+13	      1.04254	  0.000320742	          11	 -7.24792e-05

    allocate(filedata(max_line_count,8,nprobes),stat=ierr)

    ! read in the data file

    call read_lp_data(fileunit)

    ! exit after file has been read in

    close(fileunit)

  end subroutine read_lp

  subroutine read_lp_info(nprobes,max_line_count,fileunit)
    implicit none
    integer :: max_line_count,nprobes
    character*512 :: line

    integer :: ios
    integer :: fileunit
    integer :: in,ierr,pn

    real :: line_data(8)
    integer,parameter :: maxprobes = 10
    integer :: line_count(maxprobes)
    integer :: probeid(maxprobes)

    rewind(fileunit)

    ios = 0

    pn = 0
    line_count = 0
    max_line_count = 0
    nprobes = 0

    do while (ios.eq.0) 

       read(fileunit,'(a512)',iostat=ios) line

       if (ios.eq.0.and.line(1:1).eq.' ') then
          ! read in data line
          read(line,*) (line_data(in),in=1,8)
          if (pn.ne.line_data(7)) then 
             pn = line_data(7)
             nprobes = nprobes+1
             probeid(nprobes) = pn

             if (nprobes.gt.maxprobes) then 
                write(0,*) 'ERROR reading LP data: Too many probes > ', maxprobes
                stop 'Too many probes'
             endif
             line_count(nprobes) = line_count(nprobes) +1
          else
             line_count(nprobes) = line_count(nprobes) + 1
          endif

       endif

    end do

    max_line_count = maxval(line_count(1:nprobes))

    if (allocated(nlp_count)) deallocate(nlp_count)
    if (allocated(nlp_id)) deallocate(nlp_id)

    nlp = nprobes

    allocate(nlp_count(nprobes),stat=ierr)
    allocate(nlp_id(nprobes),stat=ierr)

    nlp_count = line_count(1:nprobes)
    nlp_id = probeid(1:nprobes)

  end subroutine read_lp_info


  subroutine read_lp_data(fileunit)
    implicit none

    character*512 :: line

    integer :: ios
    integer :: fileunit
    integer :: in

    real :: line_data(8)

    integer,parameter :: maxprobes = 10
    integer :: line_count(maxprobes)
    integer :: nprobes,pn

    rewind(fileunit)

    ios = 0

    pn = 0
    line_count = 0
    nprobes = 0

    do while (ios.eq.0) 

       read(fileunit,'(a512)',iostat=ios) line

       !write(0,*) 'Line:',trim(line)

       if (ios.eq.0.and.line(1:1).eq.' ') then
          ! read in data line
          read(line,*) (line_data(in),in=1,8)

          !write(0,'(a,8(2x,g12.5))') 'Line:',(line_data(in),in=1,8)

          if (pn.ne.line_data(7)) then 

             pn = line_data(7)
             nprobes = nprobes+1

             if (nprobes.gt.maxprobes) then 
                write(0,*) 'ERROR reading LP data: Too many probes > ', maxprobes
                stop 'Too many probes'
             endif

             line_count(nprobes) = line_count(nprobes) +1
          else
             line_count(nprobes) = line_count(nprobes) + 1
          endif

          filedata(line_count(nprobes),:,nprobes) = line_data

          !write(0,'(a,2i6,8(2x,g12.5))') 'Data:',nprobes,line_count(nprobes),(filedata(line_count(nprobes),in,nprobes),in=1,8)

       endif

    end do


  end subroutine read_lp_data



end module analyse_lp
