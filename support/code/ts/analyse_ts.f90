module analyse_ts
  use utilities



  real,allocatable :: filedata(:,:)


  real,parameter :: start_time1=1500.0,end_time1=2500.0
  real,parameter :: start_time2=2500.0,end_time2=4500.0

  real,parameter :: psimin=1.0, psimax=1.15
  integer,parameter :: psibins=40
  integer,parameter :: timebins=2
  integer,parameter :: ndata=2 


  real :: times(timebins,2)
  real :: deltapsi

  real :: psiaxis(psibins)
  real :: profiles(psibins,timebins,ndata)
  real :: counts(psibins,timebins,ndata)



contains


  subroutine setup_bins
    implicit none
    integer :: in

    times(1,1) = start_time1
    times(1,2) = end_time1

    times(2,1) = start_time2
    times(2,2) = end_time2

    deltapsi = (psimax-psimin)/psibins

    do in = 1,psibins
       psiaxis(in) = psimin + (in-0.5) * deltapsi
    end do

    profiles = 0.0
    counts = 0.0


  end subroutine setup_bins

  subroutine accumulate_data(nlines)
    implicit none

    integer :: nlines

    integer :: in,it
    integer :: psibin

    !write(0,*) 'Accumulate:',nlines

    do in = 1,nlines

       do it = 1,timebins
          
!          write(0,*) 'data:',in,it,filedata(in,8)

          if (filedata(in,8).ge.times(it,1).and.filedata(in,8).le.times(it,2)) then 
             psibin = max((min(filedata(in,1),psimax) - psimin),0.0) / deltapsi + 1
             ! ne
             profiles(psibin,it,1) = profiles(psibin,it,1) + filedata(in,5)
             counts(psibin,it,1) = counts(psibin,it,1) + 1.0

             ! Te
             profiles(psibin,it,2) = profiles(psibin,it,2) + filedata(in,6)
             counts(psibin,it,2) = counts(psibin,it,2) + 1.0

             !write(0,*) filedata(in,8),it,psibin

          end if
       end do
    end do



  end subroutine accumulate_data


  subroutine analyse_print_data
    implicit none
    integer :: in,it,is


    ! calculate the average profiles

    do in = 1,psibins
       do it = 1,timebins
          do is = 1,ndata
             if (counts(in,it,is).ne.0.0) then 
                profiles(in,it,is) = profiles(in,it,is) / counts(in,it,is)
             else
                profiles(in,it,is) = 0.0
             endif
          end do
       end do
    end do

    ! write out the results

    write(6,'(a)') 'Divertor Thomson Average Analysis'
    write(6,'(a,g12.5,a,g12.5,a)') 'FIRST  Time bins (Col 1,2 : ne,Te) = ', times(1,1) ,':',times(1,2),' ms'
    write(6,'(a,g12.5,a,g12.5,a)') 'SECOND Time bins (Col 3,4 : ne,Te) = ', times(2,1) ,':',times(2,2),' ms'
    write(6,'(a)')    '    PSIN           ne(1)        Te(1)     Counts(1)        ne(2)         Te(2)    Counts(2)'


    do in = 1,psibins

       write(6,'(7(2x,g12.5))') psiaxis(in),profiles(in,1,1),profiles(in,1,2),counts(in,1,1),profiles(in,2,1),profiles(in,2,2),counts(in,2,1)
    end do

    if (allocated(filedata)) deallocate(filedata)


  end subroutine analyse_print_data




  subroutine read_ts(filename,nlines)
    implicit none

    integer :: nlines

    integer :: fileunit
    integer :: ios,ierr

    character*(*) :: filename
    character*512 :: line
    integer :: line_cnt,in

    call find_free_unit_number(fileunit)

    if (allocated(filedata)) deallocate(filedata)


    open(unit=fileunit,file=filename,status='old',iostat=ios)

    if (ios.ne.0) then 
       write(0,*) 'FILE not found:',ios
    endif

    ! Assume standard format
    ! Read off header


    ! prescan for number of lines in the file


    line_cnt = 0

    do while (ios.eq.0) 

       read(fileunit,*,iostat=ios) line
       if (ios.eq.0) line_cnt = line_cnt+1

    end do

    write(0,'(a,i7)') 'Read file:'//trim(filename),line_cnt

    ! Reduce count by 2 for header lines

    line_cnt = line_cnt -2

    ! save number of data lines
    nlines = line_cnt

    ! 10 data items/line
    ! Psin	Lpol(m)	DeltaR	DeltaZ	Ne()	Te()	Ne*Te(torr)	Time_ms	R(m)	Z(m)	N_lpol
    !  1          2        3       4     5       6         7               8     9       10       11

    allocate(filedata(line_cnt,11),stat=ierr)

    ! read in the data file

    rewind(fileunit)
    line_cnt =0 
    ios = 0

    do while (ios.eq.0) 

       read(fileunit,'(a512)',iostat=ios) line
       ! write(0,'(a,a)') ':',trim(line)
       if (ios.eq.0) then 
          line_cnt = line_cnt +1
          ! read in line of data
          if (line_cnt.gt.2) then 
             read(line,*) (filedata(line_cnt-2,in),in=1,11)
          endif

       end if


    end do

    ! exit after file has been read in

    close(fileunit)


  end subroutine read_ts






end module analyse_ts
