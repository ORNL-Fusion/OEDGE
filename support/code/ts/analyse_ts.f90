module analyse_ts
  use utilities
  use error_handling

  private

  real,allocatable :: filedata(:,:)
  integer :: nlines_filedata, ncols_filedata


  !real,parameter :: start_time1=1500.0,end_time1=2500.0
  !real,parameter :: start_time2=2500.0,end_time2=4500.0

  !real,parameter :: psimin=1.0, psimax=1.15
  !integer,parameter :: psibins=40
  !integer,parameter :: timebins=2

  integer,parameter :: ndata=3
  integer,parameter :: nelm=3

  real :: start_time, end_time,psi_min,psi_max,npsi_bins


  !real :: times(timebins,2)
  real :: delta_psi

  real,allocatable :: psiaxis(:)
  real,allocatable :: profiles(:,:,:)
  real,allocatable :: counts(:,:)

  public :: setup_ts_bins,read_ts,accumulate_ts_data,analyse_print_ts_data,deallocate_ts_data


contains


  subroutine setup_ts_bins(stime,etime,psimin,psimax,npsibins)
    implicit none
    real :: stime,etime,psimin,psimax
    integer :: npsibins

    integer :: in,ierr

    write(0,'(a,4(1x,f8.3),i10)') 'TS:Start bin setup',stime,etime,psimin,psimax,npsibins

    ! Set input binning parameters

    start_time = stime
    end_time = etime
    psi_min = psimin
    psi_max = psimax
    npsi_bins = npsibins


    !times(1,1) = start_time1
    !times(1,2) = end_time1

    !times(2,1) = start_time2
    !times(2,2) = end_time2

    ! Allocate arrays for storage

    if (allocated(psiaxis)) deallocate(psiaxis)
    allocate(psiaxis(npsi_bins),stat=ierr)    
    if (ierr.ne.0) call errmsg('SETUP_TS_BINS:','Problem allocating PSIAXIS')

    delta_psi = (psi_max-psi_min)/npsi_bins

    do in = 1,npsi_bins
       psiaxis(in) = psimin + (in-0.5) * delta_psi
    end do


    if (allocated(profiles)) deallocate(profiles)
    allocate(profiles(npsi_bins,ndata,nelm),stat=ierr)    
    if (ierr.ne.0) call errmsg('SETUP_TS_BINS:','Problem allocating PROFILES')

    if (allocated(counts)) deallocate(counts)
    allocate(counts(npsi_bins,nelm),stat=ierr)    
    if (ierr.ne.0) call errmsg('SETUP_TS_BINS:','Problem allocating COUNTS')


    profiles = 0.0
    counts = 0.0


    write(0,*) 'TS:End bin setup',ierr


  end subroutine setup_ts_bins

  subroutine accumulate_ts_data(nlines,elm_filt)
    implicit none

    integer :: nlines
    logical :: elm_filt

    integer :: in,it
    integer :: psibin

    write(0,*) 'Accumulate TS Data:',nlines

    do in = 1,nlines


       !          write(0,*) 'data:',in,it,filedata(in,8)

       if (filedata(in,8).ge.start_time.and.filedata(in,8).le.end_time) then 

          ! Accumulate total
          ! Determine psibin
          psibin = min(max((min(filedata(in,1),psi_max) - psi_min),0.0) / delta_psi + 1,npsi_bins)


          ! ne
          profiles(psibin,1,1) = profiles(psibin,1,1) + filedata(in,5)
          ! Te
          profiles(psibin,2,1) = profiles(psibin,2,1) + filedata(in,6) 
          ! Rsep
          profiles(psibin,3,1) = profiles(psibin,3,1) + filedata(in,11)

          counts(psibin,1) = counts(psibin,1) + 1.0

          !write(0,*) filedata(in,8),it,psibin

          ! Accumulate INELM - select small subset close to ELM to see if there is a difference
          if (elm_filt) then 

          if (filedata(in,13).ne.0.0) then 

             ! look at +/- 1ms time window
             if (filedata(in,12).ge.-1.0.and.filedata(in,12).le.1.0) then 
                ! ne
                profiles(psibin,1,2) = profiles(psibin,1,2) + filedata(in,5)
                ! Te
                profiles(psibin,2,2) = profiles(psibin,2,2) + filedata(in,6)
                ! Rsep
                profiles(psibin,3,2) = profiles(psibin,3,2) + filedata(in,11)

                counts(psibin,2) = counts(psibin,2) + 1.0
             endif


          else

             ! Accumulate Out of ELM - outside range -5ms to 5ms 
             ! ne
             profiles(psibin,1,3) = profiles(psibin,1,3) + filedata(in,5)

             ! Te
             profiles(psibin,2,3) = profiles(psibin,2,3) + filedata(in,6)

             ! Rsep
             profiles(psibin,3,3) = profiles(psibin,3,3) + filedata(in,11)

             counts(psibin,3) = counts(psibin,3) + 1.0

          endif

          endif 

       end if
    end do

    write(0,*) 'End Accumulate TS Data:',nlines
    


  end subroutine accumulate_ts_data


  subroutine analyse_print_ts_data(ofilename,elm_filt)
    implicit none
    character*(*) :: ofilename
    logical :: elm_filt

    integer :: in,it,is
    integer :: ounit,ounit2

    character*512 :: ofile1,ofile2,ofile3
    character*512 :: tmp_file
    integer :: ierr
    logical :: res1,res2

    !  integer :: exp_tmin,exp_tmax

    !  character*5 :: tmins,tmaxs


    ! calculate the average profiles

    write(0,*) 'PRINT:Calculating Averages:'

    do in = 1,npsi_bins
       do is = 1,nelm
          if (counts(in,is).ne.0.0) then 
             profiles(in,:,is) = profiles(in,:,is) / counts(in,is)
          else
             profiles(in,:,is) = 0.0
          endif
       end do
    end do

    write(0,*) 'PRINT: Averaging complete:'


    ! write out the results

    ! Two sets are required 
    ! a) revised TS data file unchanged except for the addition of ELM time data
    ! b) file containing the bin averaged data ... all 3 sets 

    ! The input file name is used as the stem 

    ! first output file

    if (elm_filt) then 
       ! print modified TS file with elm filter data added

       tmp_file = trim(ofilename)
       ofile1 = trim(tmp_file)//'.elm'

       call find_free_unit_number(ounit) 

       write(0,*) 'First output:'//trim(ofile1),':',ounit,':'

       open(unit=ounit,file=trim(ofile1),access='sequential',iostat=ierr)

       write(ounit,'(1x,13(2x,a,2x))') 'Psin','Lpol(m)','DeltaR','DeltaZ','Ne','Te','Ne*Te(torr)','Time(ms)','R(m)','Z(m)','Rsep(m)','Elm_offset(ms)','N_Elm'

       do in = 1,nlines_filedata
          write(ounit,'(1x,13(1x,g13.6))') (filedata(in,is),is=1,ncols_filedata)
       end do

       close(ounit)

       write(0,*) 'First output complete:'

    endif


    !  exp_tmin = int(log10(start_time))+1

    !  if (exp_tmin.eq.3) then 
    !     write(tmins(1:3),'(i3)') int(start_time)
    !  elseif (exp_tmin.eq.4) then
    !     write(tmins(1:4),'(i4)') int(start_time)
    !  endif


    !  exp_tmax = int(log10(end_time))+1

    !  if (exp_tmax.eq.3) then 
    !     write(tmaxs(1:3),'(i3)') int(end_time)
    !  elseif (exp_tmax.eq.4) then
    !     write(tmaxs(1:4),'(i4)') int(end_time)
    !  endif

    !  write(0,*) 'Times:',trim(tmins),':',trim(tmaxs),':'


    ! second output file


    write(ofile2,'(a,a,i4,a,i4,a)') trim(ofilename),'_',int(start_time),'-',int(end_time),'.average'

    call find_free_unit_number(ounit2) 

    write(0,*) 'Second output:',trim(ofile2),':',ounit2,':'

    open(unit=ounit2,file=trim(ofile2),access='sequential',iostat=ierr)
    inquire(unit=ounit2,opened=res1)
    inquire(file=trim(ofile2),opened=res2)
    if (res1.ne.res2) call errmsg('PRINT TS:','CAUTION: Second output file name is not correct')


    write(ounit2,'(a)') ' Thomson Average Analysis'
    write(ounit2,'(a,f12.3,a,f12.3,a)') ' Time range = ', start_time ,':',end_time,' ms'
    write(ounit2,'(a)') ' All Data binned:'
    write(ounit2,'(a)') '    PSIN    Rsep       ne        Te     Counts'


    do in = 1,npsi_bins

       if (counts(in,1).gt.0) then 
          write(ounit2,'(7(2x,g12.5))') psiaxis(in),profiles(in,3,1),profiles(in,1,1),profiles(in,2,1),counts(in,1)
       endif

    end do

    if (elm_filt) then 
       ! print out averaged data related to ELMs if elm filtering performed

       write(ounit2,*)
       write(ounit2,'(a)') ' ELM Data binned +/-1ms:'
       write(ounit2,'(a)') '    PSIN     Rsep       ne        Te     Counts'

       do in = 1,npsi_bins

          if (counts(in,2).gt.0) then 
             write(ounit2,'(7(2x,g12.5))') psiaxis(in),profiles(in,3,2),profiles(in,1,2),profiles(in,2,2),counts(in,2)
          endif

       end do

       write(ounit2,*)
       write(ounit2,'(a)') ' Non-ELM Data binned >5ms:'
       write(ounit2,'(a)') '    PSIN     Rsep      ne        Te     Counts'

       do in = 1,npsi_bins

          if (counts(in,3).gt.0) then 
             write(ounit2,'(7(2x,g12.5))') psiaxis(in),profiles(in,3,3),profiles(in,1,3),profiles(in,2,3),counts(in,3)
          endif

       end do


    endif


    close(ounit2)
    write(0,*) 'Second output complete:'



  end subroutine analyse_print_ts_data




  subroutine read_ts(filename,nlines,elm_filt,elm_filename)
    implicit none

    integer :: nlines
    logical :: elm_filt
    character*(*) :: elm_filename

    integer :: fileunit
    integer :: ios,ierr

    character*(*) :: filename
    character*512 :: line
    integer :: line_cnt,in

    integer :: ncols, nextra, fileformat,headersize
    real,allocatable :: line_data(:)


    write(0,*) 'READ_TS: Reading TS data'


    call find_free_unit_number(fileunit)

    if (allocated(filedata)) deallocate(filedata)


    open(unit=fileunit,file=filename,status='old',iostat=ios)

    if (ios.ne.0) then 
       write(0,*) 'FILE not found:',ios
       stop 'STOP: Data file not found:'//trim(filename)
    endif

    ! Assume standard format
    ! Read off header

    ! Read header and decide which format of file is being loaded ... 
    
    read(fileunit,'(a512)') line
    rewind(fileunit)
    if (line(1:4) .eq. 'time') then 
    ! file format = 2
    ! 11 data items/line but different
       ! 1       2      3      4      5     6     7   8    9     10    11
       !time   lpol    psin   rsep    r     z    den  te  press inner core


       fileformat = 2
       headersize = 1

       ncols = 11
       nextra = 2
       ncols_filedata = ncols+nextra

    else
    ! fileformat = 1
    !
    ! up to 11 data items/line
    ! Psin	Lpol(m)	DeltaR	DeltaZ	Ne()	Te()	Ne*Te(torr)	Time_ms	R(m)	Z(m)	N_lpol/Rsep(m)
    !  1          2        3       4     5       6         7               8     9       10       11


       fileformat = 1
       headersize = 2

       ncols = 11
       nextra = 2
       ncols_filedata = ncols+nextra

    endif

    ! prescan for number of lines in the file


    line_cnt = 0

    do while (ios.eq.0) 

       read(fileunit,*,iostat=ios) line
       if (ios.eq.0) line_cnt = line_cnt+1

    end do

    write(0,'(a,i7)') 'Read file:'//trim(filename),line_cnt

    ! Reduce count by headersize

    line_cnt = line_cnt -headersize

    ! save number of data lines
    nlines = line_cnt
    nlines_filedata = line_cnt


    ! allocate space for elm data as well 
    !allocate(filedata(line_cnt,11),stat=ierr)
    allocate(filedata(line_cnt,ncols_filedata),stat=ierr)

    allocate(line_data(ncols_filedata),stat=ierr)

    ! initialize
    filedata = 0.0

    ! read in the data file

    rewind(fileunit)
    line_cnt =0 
    ios = 0

    do while (ios.eq.0) 

       read(fileunit,'(a512)',iostat=ios) line
       ! write(0,'(a,a)') ':',trim(line)
       if (ios.eq.0) then 
          line_cnt = line_cnt +1
          ! read in line of data - 11 data elements
          if (line_cnt.gt.headersize) then 
             !read(line,*) (filedata(line_cnt-2,in),in=1,11)
             read(line,*) (line_data(in),in=1,ncols)

             if (fileformat.eq.1) then 
                filedata(line_cnt-headersize,:) = line_data
             elseif (fileformat.eq.2) then 
                ! psin
                filedata(line_cnt-headersize,1) = line_data(3)
                ! lpol
                filedata(line_cnt-headersize,2) = line_data(2)
                ! deltar -> inner
                filedata(line_cnt-headersize,3) = line_data(10)
                ! deltaz -> core
                filedata(line_cnt-headersize,4) = line_data(11)
                ! ne
                filedata(line_cnt-headersize,5) = line_data(7)
                ! te
                filedata(line_cnt-headersize,6) = line_data(8)
                ! press (ne*Te)
                filedata(line_cnt-headersize,7) = line_data(9)
                ! time (ms)
                filedata(line_cnt-headersize,8) = line_data(1)
                ! r
                filedata(line_cnt-headersize,9) = line_data(5)
                ! z
                filedata(line_cnt-headersize,10) = line_data(6)
                ! n_lpol -> rsep
                filedata(line_cnt-headersize,11) = line_data(4)
              endif

          endif

       end if

    end do

    ! exit after file has been read in

    close(fileunit)

    ! filter elms - assign elm index and temporal offset to nearest ELM

    

    if (elm_filt) then 
       write(0,*) 'READ_TS: Filtering ELMs:'
       call filter_ts_data(elm_filename,filedata,nlines,ncols,nextra)
    endif


    write(0,*) 'READ_TS: Reading complete:'

    deallocate(line_data)


  end subroutine read_ts



  subroutine filter_ts_data(elm_filename,ts_data,nlines,ncols,nextra)
    use filter_elms
    implicit none
    integer :: nlines,ncols,nextra
    character*(*) elm_filename
    real,allocatable :: ts_data(:,:)
    !real lp_data(nlines,ncols+nextra)
    integer :: ierr
    integer :: in
    real :: elm_time_offset
    integer :: elmref
    logical :: inelm

    !sd8
    ! Load ELM time data from file and then add ELM information to the TS data 
    !

    real :: elm_start_input,elm_end_input, elm_effect_start_input, elm_effect_end_input

    
    write(0,*) 'FILTER_TS_DATA: Start:',nlines,ncols,nextra

    ! try elm time window of -5 to +5
    elm_start_input = -5.0
    elm_end_input   =  5.0
    elm_effect_start_input = -5.0
    elm_effect_end_input   =  5.0

    call load_elm_times(elm_filename,ierr)

    if (ierr.ne.0) return


    call set_elm_criteria(elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)

    ! An elmref value of 0 means that it is outside time window that is associated with an ELM
    ! A negative elmref value is inside a possible ELM effect window
    ! A positive elmref value is inside the most likely ELM effect window
    do in = 1, nlines

       call get_elm_time(ts_data(in,8),elm_time_offset,elmref,inelm)

       ts_data(in,ncols+1) = elm_time_offset
       ts_data(in,ncols+2) = 0
       if (inelm) then 
          ts_data(in,ncols+2) = elmref
       else
          ts_data(in,ncols+2) = -elmref
       endif

    end do

    write(0,*) 'FILTER_TS_DATA: End'

  end subroutine filter_ts_data



  subroutine deallocate_ts_data

    if (allocated(filedata)) deallocate(filedata)
    if (allocated(psiaxis)) deallocate(psiaxis)
    if (allocated(profiles)) deallocate(profiles)
    if (allocated(counts)) deallocate(counts)

  end subroutine deallocate_ts_data


end module analyse_ts
