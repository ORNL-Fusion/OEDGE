module langmuir_probes
  use error_handling
  implicit none


contains

  subroutine read_lp_data_file(iunit,lp_data,nlines,ncols)
    implicit none
    integer :: iunit,nlines,ncols
    real, allocatable :: lp_data(:,:)
    real, allocatable :: lp_axis(:),lp_proc_data(:,:)
    character*512 :: line
    integer :: ios,line_cnt,ierr,in,it


    ! determine number of lines of data - 2 header lines at beginning assumed

    !rewind(iunit)
    read(iunit,'(a)') line
    read(iunit,'(a)') line

    line_cnt = 0
    do while (ios.eq.0) 

       read(iunit,*,iostat=ios)
       if (ios.eq.0) then 
          line_cnt = line_cnt + 1
       endif

    end do


    ! Allocate storage
    if (allocated(lp_data)) deallocate(lp_data)
    allocate(lp_data(line_cnt,ncols),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('READ_LP_DATA_FILE:','ERROR ALLOCATING LP_DATA')
       stop 
    endif

    ! rewind input

    rewind(iunit)

    ! read input into array - remove header lines

    read(iunit,*)
    read(iunit,*)

    ! data format is expected to be:
    !   time(msec)	jsat(A/cm2)	temp(eV)	dens(cm-3)  		csq		probeID		delrsepin	delzsepin	delrsepout

    do in = 1,line_cnt
       read(iunit,*,iostat=ios) (lp_data(in,it),it=1,ncols)
    end do

    nlines = line_cnt


  end subroutine read_lp_data_file


  subroutine bin_lp_data_r(lp_axis,lp_proc_data,npts,ndata,lp_data,nlines,ncols,deltar,tmin,tmax,chisq_lim)
    implicit none
    real,allocatable :: lp_axis(:), lp_proc_data(:,:) 
    real,allocatable :: lp_data(:,:)
    integer :: nlines, ncols,ndata,npts
    real :: deltar
    real :: tmin,tmax,chisq_lim
    integer :: in,nrbins, ibin,ierr
    real :: rmin,rmax


    ! For data in the range Tmin to Tmax ... bin and average over range of deltar. 
    ! Averaging Te and Jsat
    
    if (.not.allocated(lp_data)) then 
       call errmsg('BIN_LP_DATA_R','LP_DATA NOT ALLOCATED')
       stop
    endif

    ! determine minR and maxR values in the desired range. 

    ! set max and min to out of range values
    rmin = 1e25
    rmax = -1e25

    do in = 1,nlines
       ! correct time window
       if (lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax.and.lp_data(in,5).le.chisq_lim) then 
          rmin = min(lp_data(in,9),rmin)
          rmax = max(lp_data(in,9),rmax)
       endif
    end do

    
    ! determine number of R bins
    nrbins = (rmax-rmin)/deltar + 1

    ! allocate storage for processed data
    
    ndata = 2
    npts = nrbins
    if (allocated(lp_axis)) deallocate(lp_axis)
    allocate(lp_axis(nrbins),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA_R','ERROR ALLOCATING LP_AXIS')
       stop
    endif
    lp_axis = 0.0


    if (allocated(lp_proc_data)) deallocate(lp_proc_data)
    allocate(lp_proc_data(nrbins,ndata+1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA_R','ERROR ALLOCATING LP_AXIS')
       stop
    endif
    lp_proc_data = 0.0


    ! loop through data and record results - use lp_axis to record counts on first run through

    do in = 1,nlines
       if (lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax.and.lp_data(in,5).le.chisq_lim) then 
          ibin = int((lp_data(in,9)-rmin)/deltar) + 1
          ! record data count
          lp_proc_data(ibin,ndata+1) = lp_proc_data(ibin,ndata+1) + 1.0
          ! record data sum
          ! jsat
          lp_proc_data(ibin,1) = lp_proc_data(ibin,1) + lp_data(in,2)
          ! te
          lp_proc_data(ibin,2) = lp_proc_data(ibin,2) + lp_data(in,3)

       endif
    end do

    ! loop through - calculate averages and assign axis values

    do in = 1,nrbins
       if (lp_proc_data(in,ndata+1).gt.0.0) then 
          lp_proc_data(in,1) = lp_proc_data(in,1)/lp_proc_data(in,ndata+1)
          lp_proc_data(in,2) = lp_proc_data(in,2)/lp_proc_data(in,ndata+1)
       else
          lp_proc_data(in,:) = 0.0
          !lp_proc_data(in,2) = 0.0
       endif

       ! assign axis value
       lp_axis(in) = (real(in) - 0.5)* deltar + rmin 

    end do 


  end subroutine bin_lp_data_r


  subroutine print_lp_bin_data(ounit,lp_axis,lp_proc_data,npts,ndata,ident)
    implicit none
    
    real,allocatable :: lp_axis(:),lp_proc_data(:,:)
    integer :: npts,ndata,ounit,in
    character*(*) :: ident

    if (.not.allocated(lp_axis)) then 
       call errmsg('PRINT_LP_BIN_DATA','LP_AXIS NOT ALLOCATED') 
       stop
    endif

    ! print out the processed/binned data
    ! at the present time this data is  LP_AXIS   LP_PROC_DATA(,1) = Jsat   LP_PROC_DATA(,2) = Te

    write(ounit,'(a,a)') 'ID : ',trim(ident)
    write(ounit,'(a)')   'R-Rsep          Jsat(A/cm2)          Te(eV)         Count'

    do in = 1,npts
       if (lp_proc_data(in,1).gt.0.0) then 
          write(ounit,'(10(1x,g18.8))') lp_axis(in),lp_proc_data(in,1),lp_proc_data(in,2),lp_proc_data(in,ndata+1)
       endif
    end do 


  end subroutine print_lp_bin_data




end module langmuir_probes
