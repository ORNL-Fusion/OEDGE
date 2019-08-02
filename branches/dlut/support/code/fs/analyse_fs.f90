module analyse_fs
  use utilities

  private

  real,allocatable :: filedata(:,:)


  integer,parameter :: timebins=44
  real,parameter :: tstart = 2300
  real,parameter :: deltat = 50
  
  !integer,parameter :: ndata=2 

  real :: times(timebins,2)
  real :: base_profiles(timebins)
  real :: base_counts(timebins)
  real :: base_int_axis(timebins)

  real :: delta_intensity

  real,allocatable :: intensity_axis(:)
  real,allocatable :: profiles(:,:)
  real,allocatable :: counts(:,:)
  !real,allocatable :: averages(:,:)

  real :: maxrange,minrange
  integer :: data_bins
  character*(100) :: fschan

  public :: setup_fs,accumulate_data, read_fs, analyse_print_data

  real :: mindata,maxdata


contains


  subroutine setup_fs(minval,maxval,nbins,channel)
    implicit none
    real :: minval,maxval
    integer :: in,nbins
    character*(*) :: channel
    integer :: ierr,it

    if (allocated(intensity_axis)) deallocate(intensity_axis)
    if (allocated(profiles)) deallocate(profiles)
    if (allocated(counts)) deallocate(counts)
    !if (allocated(averages)) deallocate(averages)


    ! allocate space for data accumulation

    allocate(intensity_axis(nbins),stat=ierr)
    allocate(profiles(nbins,timebins),stat=ierr)
    allocate(counts(nbins,timebins),stat=ierr)
    !allocate(averages(nbins,timebins),stat=ierr)


    write(0,'(a,i6,4(2x,g12.5))') '1 Setup:Channel:'//trim(channel),nbins,minval,maxval

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

    base_profiles = 0.0
    base_counts = 0.0


    !times(1,1) = 2300
    !times(1,2) = 2500

    !times(2,1) = 2700
    !times(2,2) = 3000

    !times(3,1) = 3000
    !times(3,2) = 3500

    !times(4,1) = 3500
    !times(4,2) = 4000

    !times(5,1) = 4000
    !times(5,2) = 4500


    mindata =  1e30
    maxdata = -1e30

    maxrange = maxval
    minrange = minval

    data_bins = nbins

    fschan = channel

    !times(2,1) = start_time2
    !times(2,2) = end_time2

    ! Set up binning information

    delta_intensity = (maxrange-minrange)/nbins

    do in = 1,data_bins
       intensity_axis(in) = minval + (in-0.5) * delta_intensity
    end do

    ! zero out data accumulation arrays

    profiles = 0.0
    counts = 0.0
    !averages = 0.0

    !write(0,'(a)') '2 Setup:Channel:'//trim(channel)

  end subroutine setup_fs

  subroutine accumulate_data(nlines)
    implicit none

    integer :: nlines

    integer :: in,it
    integer :: bin

    !write(0,*) 'Accumulate:',nlines

    do in = 1,nlines

       do it = 1,timebins
          
          !write(0,*) 'data:',in,it,filedata(in,4),filedata(in,5)

          if (filedata(in,4).ge.times(it,1).and.filedata(in,4).le.times(it,2)) then 
             bin = min(int(max((min(filedata(in,5),maxrange) - minrange),0.0) / delta_intensity + 1),data_bins)
             ! binned intensity
             profiles(bin,it) = profiles(bin,it) + filedata(in,5)
             counts(bin,it) = counts(bin,it) + 1.0

             !write(0,*) 'data:',in,it,filedata(in,4),filedata(in,5)
             !write(0,*) filedata(in,8),it,bin

          end if
       end do
    end do

  end subroutine accumulate_data


  subroutine analyse_print_data(filenameout,printopt)
    implicit none
    character*(*) :: filenameout
    integer :: in,it,printopt
    integer :: baserange,peaknum
    integer :: basepos(1)
    real :: base_count, base_intensity
    real, allocatable :: tmp_count(:)

    integer :: outunit,ios
    integer :: elm_bin
    real :: total_count,elm_count


    call find_free_unit_number(outunit)

    open(unit=outunit,file=trim(filenameout),iostat=ios)

    write(0,'(a,2x,i10)') 'Output file:'//trim(filenameout),data_bins

    elm_count = 0.0
    total_count = 0.0


    ! calculate the average profiles

    do in = 1,data_bins
       do it = 1,timebins
          !do is = 1,ndata
          if (counts(in,it).ne.0.0) then 
             profiles(in,it) = profiles(in,it) / counts(in,it)
          else
             profiles(in,it) = 0.0
          endif
          !end do

          !write(0,'(a,i6,5(2x,g12.5))') 'Profiles:',in,it,profiles(in,it),counts(in,it)
       end do
    end do

    ! write out the results


    write(outunit,'(a)') 'Filter Scope Analysis: Channel = ',trim(fschan)

    if (printopt.ne.0) then 
       write(outunit,'(a,2(2x,g12.5))') 'Maximum and minimum range imposed for bins:',maxrange,minrange
       write(outunit,'(a,2(2x,g12.5))') 'Maximum and minimum range found in data   :',maxdata,mindata
    endif


    do it = 1,timebins

       if (printopt.ne.0) then 
          write(outunit,'(a,i6,a,g12.5,a,g12.5,a)') 'Time bins (Set # ',it,' ) = ', times(it,1) ,':',times(it,2),' ms'
       endif

       if (allocated(tmp_count)) deallocate(tmp_count)

       allocate(tmp_count(data_bins))
       tmp_count = counts(:,it)

       basepos = maxloc(tmp_count)
       peaknum = basepos(1)

       !write(0,*) 'basepos=',peaknum

       baserange = 1

       base_intensity = 0.0
       base_count = 0.0


       do in = max(peaknum-baserange,1),min(peaknum+baserange,data_bins)
          base_intensity = base_intensity + profiles(in,it) * counts(in,it)
          base_count = base_count + counts(in,it)
       end do

       if (base_count.ne.0.0) then 
          base_intensity = base_intensity / base_count
       else
          base_intensity = 0.0
       endif

       base_profiles(it) = base_intensity
       base_counts(it) = base_count


       
       ! Estimate total ELM time by summing counts greater than 1.2 X base level intensity 
       elm_bin =  min(int(max((min(1.2 * base_intensity,maxrange) - minrange),0.0) / delta_intensity + 1),data_bins)
       elm_count = elm_count + sum(tmp_count(elm_bin:data_bins))
       total_count = total_count + sum(tmp_count)


       if (printopt.ne.0) then 
          write(outunit,'(2(a,2x,g18.8))')  'Location of max counts = ',basepos,        ' Width used in Intensity Calculation +/- ', baserange
          write(outunit,'(2(a,2x,g18.8))')  'Base Intensity Average = ',base_intensity, ' Counts contributing                   = ', base_count


          write(outunit,'(a)')    '    Counts      Intensity_Axis       Average Intensity '


          do in = 1,data_bins
             write(outunit,'(7(2x,g12.5))') counts(in,it),intensity_axis(in),profiles(in,it)
          end do
       endif

    end do


    ! Summary of baseline emission data
    write(outunit,'(a)') 'Summary of baseline emission data over time:'
    write(outunit,'(a,2x,g12.5)') 'Width of time windows (ms): ', deltat
    write(outunit,'(a,3(2x,g12.5))') 'Fraction of ELM time - ELM counts - Total Counts - % ELM : ', elm_count, total_count, elm_count/total_count*100.0
    write(outunit,'(a)')    '    Counts           Time(ms)      Average Intensity '
    do it = 1,timebins
       write(outunit,'(7(2x,g12.5))') base_counts(it),base_int_axis(it),base_profiles(it)
    end do

    if (allocated(filedata)) deallocate(filedata)


  end subroutine analyse_print_data




  subroutine read_fs(filename,channel,nlines)
    implicit none

    integer :: nlines

    integer :: fileunit
    integer :: ios,ierr

    character*(*) :: filename
    character*(*) :: channel
    integer :: line_cnt,in

    call find_free_unit_number(fileunit)

    if (allocated(filedata)) deallocate(filedata)

    write(0,*) trim(filename)

    open(unit=fileunit,file=filename,status='old',iostat=ios)

    if (ios.ne.0) then 
       write(0,*) 'FILE not found:',ios
    endif

    ! prescan for number of lines in the file
    ! find data block of interest

    call read_fs_chan(line_cnt,fileunit,channel)

    ! save number of data lines
    nlines = line_cnt

    ! 5 data items/line
    !  Index1 Index2  Index3         Time(ms)       Intensity(ph/m2/s) 
    !   0       0       1           1500.1801       1.4626619e+17

    write(0,*) 'line_cnt:',line_cnt

    allocate(filedata(line_cnt,5),stat=ierr)

    ! read in the data file

    call read_fs_chan(line_cnt,fileunit,channel)

    ! exit after file has been read in

    close(fileunit)

    write(0,'(a,i7,4a)') 'Read file:'//trim(filename),line_cnt,'  Channel:',trim(fschan),':',trim(channel)

  end subroutine read_fs

  subroutine read_fs_chan(line_cnt,fileunit,channel)
    implicit none
    integer :: line_cnt
    character*(*) :: channel
    character*512 :: line
    logical :: channel_found
    integer :: ios
    integer :: fileunit
    integer :: a,in
    logical :: load_data

    load_data = .false.

    if (allocated(filedata)) load_data = .true.


    rewind(fileunit)

    ios = 0
    line_cnt = 0
    channel_found = 0

    a=len_trim(channel)

    write(0,*) 'Chan:',trim(channel)


    do while (ios.eq.0) 

       read(fileunit,'(a512)',iostat=ios) line

       if (ios.eq.0.and.channel_found) then
          !write(0,*) 'found:',load_data,trim(line)

          line_cnt = line_cnt + 1
       endif


       if (line(1:11).eq.'FS CHANNEL:') then 
          ! found a filter scope channel - check to see if it is the correct one
          if (channel_found) then 
             ! if reached next set of FS data then exit
             ios = 1
             channel_found = .false.
             line_cnt = line_cnt -1 
          elseif (line(13:min(13+a-1,len(line))).eq.trim(channel)) then 
             channel_found = .true.
          endif

       endif

       if (ios.eq.0.and.channel_found.and.line_cnt.gt.0) then 

          ! if the data space has been allocated then fill it
          if (load_data) then 

             read(line,*) (filedata(line_cnt,in),in=1,5)
             !write(0,'(a,5(2x,g12.5))') 'FS:'//trim(fschan),(filedata(line_cnt,in),in=1,5)
             mindata = min(mindata,filedata(in,5))
             maxdata = max(maxdata,filedata(in,5))

          endif

       endif






    end do


  end subroutine read_fs_chan





end module analyse_fs
