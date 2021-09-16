module langmuir_probes
  use error_handling
  implicit none

  integer,parameter :: RBINOPT=0
  integer,parameter :: PSIBINOPT=1



  ! Initialization for current LP data file format
  ! Add three columns for including additional ELM information
  ! Added floating potential
  ! ncols  = 9 
  integer, parameter :: ncols  = 10 
  ! nextra - extra colums added to lp_data - 3 for ELM data - 1 for inner/outer identifier - 1 for flagging as outlier removed
  integer, parameter :: nextra = 6

  ! PSI is in lp_data(in,9) and R-Rsep_out is in LP_DATA(in,8)
  ! R-Rsep_in is in LP_DATA(in,7)
  ! vf is in lp_data(in,10)
  integer,parameter,private ::  timebin  = 1
  integer,parameter,private ::  jsatbin  = 2
  integer,parameter,private ::  tebin    = 3
  integer,parameter,private ::  nebin    = 4
  integer,parameter,private ::  chibin   = 5
  integer,parameter,private ::  probebin = 6
  integer,parameter,private ::  irbin    = 7  
  integer,parameter,private ::  rbin     = 8
  integer,parameter,private ::  psibin   = 9
  integer,parameter,private ::  vfbin    = 10

  integer,parameter,private ::  nzones = 2
  integer,parameter ::  inner   = 1
  integer,parameter ::  outer   = 2

  integer,parameter,private :: elm_time_col = 1
  integer,parameter,private :: elm_ref_col = 2
  integer,parameter,private :: elm_frac_col = 3
  integer,parameter,private :: zone_col = 4
  integer,parameter,private :: outlier_index = 5
  integer,parameter,private :: bin_col = 6

  ! allocate storage for processed data
  ! ndata = 3 ... jsat,te,vf
  integer,parameter :: ndata = 3 

  character*3,parameter :: nan = 'NaN'
  character*7,parameter :: line_form = '(a1024)'

  real,parameter :: psimin_limit = -0.5
  real,parameter :: psimax_limit =  2.0
  logical,parameter :: test_rlimits = .true.
  real,parameter :: rdiff_limit =  0.02   ! if the data point is within 0.02m of both strike points then exclude it
  real,parameter :: r_limit = 0.8         ! if the data point is greater than 0.8 meters from both strike points then exclude it

  ! Time filter variables
  logical:: filter_times
  integer :: time_window_cnt 
  real,allocatable :: time_windows(:,:)



contains

  subroutine read_lp_data_file(iunit,lp_data,tmin,tmax,nlines)
    implicit none
    integer :: iunit,nlines
    real :: tmin,tmax
    !
    ! Locals
    !
    real, allocatable :: lp_data(:,:)
    character*1024 :: line
    integer :: ios,tot_line_cnt,act_line_cnt,cur_line_cnt,ierr,in,it,nan_loc
    real,allocatable :: tmp_lp_data(:)


    write(0,'(a)') 'Reading LP Data file: Lines containing NaN are discarded'

    ! allocate temporary storage
    if (allocated(tmp_lp_data)) deallocate(tmp_lp_data)
    allocate(tmp_lp_data(ncols+nextra),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('READ_LP_DATA_FILE:','ERROR ALLOCATING TMP_LP_DATA')
       stop 
    endif

    ! determine number of lines of data - 2 header lines at beginning assumed

    !rewind(iunit)
    read(iunit,'(a)') line
    read(iunit,'(a)') line

    tot_line_cnt = 0
    act_line_cnt = 0
    ios = 0

    do while (ios.eq.0) 

       read(iunit,line_form,iostat=ios) line
       if (ios.eq.0) then 
          tot_line_cnt = tot_line_cnt + 1
          nan_loc=index(line,nan)
          if (nan_loc.eq.0) then 
             read(line,*) (tmp_lp_data(it),it=1,ncols)
             ! filter out invalid PSIn values and take only the desired time window
             if (check_lp_valid(tmp_lp_data,ncols,tmin,tmax)) then 
                act_line_cnt =  act_line_cnt + 1
             endif
          endif
       endif

    end do

    write(0,*) 'Total Lines = ',tot_line_cnt,' Data Lines = ',act_line_cnt

    ! Allocate storage
    if (allocated(lp_data)) deallocate(lp_data)
    allocate(lp_data(act_line_cnt,ncols+nextra),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('READ_LP_DATA_FILE:','ERROR ALLOCATING LP_DATA')
       stop 
    endif

    ! initialize the allocated storage
    lp_data = 0.0

    ! rewind input

    rewind(iunit)

    ! read input into array - remove header lines

    read(iunit,*)
    read(iunit,*)

    ! data format is expected to be:  
    ! old version 1 : 9 columns
    !   time(msec)	jsat(A/cm2)	temp(eV)	dens(cm-3)  		csq		probeID		delrsepin	delzsepin	delrsepout
    ! new version : 9 columns
    !   time(msec)	jsat(A/acm2)	temp(eV)	dens(cm-3)              csq	        probeID	        delrsepin	delrsepout	psin
    ! new version 2 : 10 columns
    !   time(msec)	jsat(A/acm2)	temp(eV)	dens(cm-3)              csq	        probeID	        delrsepin	delrsepout	psin    vf(V)

    cur_line_cnt = 0

    do in = 1,tot_line_cnt
       ! need to filter out and reject data lines that contain NaN quantities
       ! need to filter out and reject data lines with out of range PSIn values (accept -1, 2)
       read(iunit,line_form,iostat=ios) line
       
       nan_loc=index(line,nan)

       if (nan_loc.eq.0) then 

          tmp_lp_data = 0.0

          read(line,*,iostat=ios) (tmp_lp_data(it),it=1,ncols)

          if (check_lp_valid(tmp_lp_data,ncols,tmin,tmax)) then 

              cur_line_cnt = cur_line_cnt + 1
              ! If all the valid lines have already been read then exit
              if (cur_line_cnt.gt.act_line_cnt) exit

              ! copy over the data
              lp_data(cur_line_cnt,:) = tmp_lp_data

          endif
       endif

    end do

    nlines = act_line_cnt

    write(0,'(a,i8)') 'Finished: Reading LP Data file: nlines=',nlines

  end subroutine read_lp_data_file

  logical function check_lp_valid(tmp_lp_data,ncols,tmin,tmax)
    implicit none
    integer :: ncols
    real :: tmp_lp_data(ncols)
    real :: tmin,tmax

    check_lp_valid = .false.
    
    if (check_time_filter(tmp_lp_data(timebin)).and.&    ! check data point is in requested specific time windows
        (tmp_lp_data(timebin).ge.tmin.and.tmp_lp_data(timebin).le.tmax).and.& ! check data point in overall time frame
        (tmp_lp_data(psibin).ge.psimin_limit.and.tmp_lp_data(psibin).le.psimax_limit).and.& ! check data point is within valid range of PSI
        (test_rlimits.and.&
         ((abs(tmp_lp_data(irbin)-tmp_lp_data(rbin)).gt.rdiff_limit).and.&   ! check that the del_rsepin and del_rsepout are not the same
          (abs(tmp_lp_data(irbin)).lt.r_limit).and.(abs(tmp_lp_data(rbin)).lt.r_limit)))&
       ) then 
       check_lp_valid = .true.
    endif

  end function check_lp_valid


  subroutine flag_lp_data(lp_data,nlines)
    implicit none
    integer :: nlines
    real, allocatable :: lp_data(:,:)
    real*8 :: r_av,r_cnt,r_av2,r_cnt2
    real*8 :: max_drsep, min_drsep
    real*8 :: max_drsep_pfz, min_drsep_pfz
    real*8 :: max_drsep_sol, min_drsep_sol

    integer :: in,psiminind

    real :: psimin
    real :: rsplit
    real*8,parameter:: num_minval = -1e25 
    real*8,parameter:: num_maxval = 1e25 


    write(0,'(a)') 'Flagging Inner/Outer LP Data:'


    ! NOTE: This algorithm fails with sufficient sweeping ...also may fail for vertical inner target
    !       Lets try finding the average R-Rsep for all PSI<1 elements and then divide by 2. 
    !
    ! Loop through lp data and find the R-Rsep value associated with the minimum PSI value. This should be in the PFZ and should be one of the 
    ! two points defining the location between inner and outer



    psimin = 1e25
    psiminind = 0
    r_cnt = 0.0
    r_av = 0.0
    r_cnt2 = 0.0
    r_av2 = 0.0
    max_drsep = -num_minval
    max_drsep_pfz = -num_minval
    max_drsep_sol = -num_minval

    min_drsep = num_maxval
    min_drsep_pfz = num_maxval
    min_drsep_sol = num_maxval




    if (allocated(lp_data)) then 

       do in = 1,nlines

!          if (lp_data(in,psibin).lt.psimin) then
!             psiminind = in
!             psimin = lp_data(in,psibin)
!             rsplit = lp_data(in,rbin)
!          endif

          ! Record minimum R value in data in case there is no data inside the PFZ          
          max_drsep = max(max_drsep,dble(lp_data(in,rbin)))
          min_drsep = min(min_drsep,dble(lp_data(in,rbin)))

          if (lp_data(in,psibin).lt.1.0) then 
             min_drsep_pfz = min(min_drsep_pfz,dble(lp_data(in,rbin)))
             max_drsep_pfz = max(max_drsep_pfz,dble(lp_data(in,rbin)))

!             r_av = r_av + lp_data(in,rbin)
!             r_cnt = r_cnt + 1.0
          elseif (lp_data(in,psibin).ge.1.0) then 
             min_drsep_sol = min(min_drsep_sol,dble(lp_data(in,rbin)))
             max_drsep_sol = max(max_drsep_sol,dble(lp_data(in,rbin)))

!             r_av2 = r_av2 + lp_data(in,rbin)
!             r_cnt2 = r_cnt2 + 1.0
          endif

       end do


!       if (r_cnt.ne.0) then 
!          rsplit = (r_av/r_cnt)
!       elseif (min_drsep.ge.-0.01.and.max_drsep.gt.0) then 
!          rsplit = min_drsep
!       elseif (min_drsep.lt.max_drsep.and.max_drsep.lt.-0.01) then 
!          rsplit = max_drsep
!       else 
!          rsplit = (r_av2/r_cnt2)
!       endif


       write(0,'(a,4g18.8)') 'SEPS:',min_drsep_sol,max_drsep_sol,min_drsep_pfz,max_drsep_pfz

       
       ! Inner target data is included ...
       if (min_drsep_sol.lt.min_drsep_pfz) then 

          ! Outer target data is included    
          if (max_drsep_sol.gt.max_drsep_pfz) then 
             rsplit = (min_drsep_pfz+max_drsep_pfz)/2.0

          ! Inner only    
          else
             rsplit = max_drsep_pfz
          endif
             
       else

          ! Outer target only
          if (max_drsep_sol.gt.max_drsep_pfz) then 
             rsplit = min_drsep_pfz

          endif

       endif

       rsplit = -0.04


       ! rsplit now defines the splitting point between inner and outer - loop through the data and apply flag
       
       write(0,'(a,g12.5)') 'FLAGGING INNER/OUTER: SPLIT R-RSEP VALUE = ',rsplit

       do in = 1,nlines

          if (lp_data(in,rbin).lt.rsplit) then
             lp_data(in,ncols+zone_col) = INNER
          else
             lp_data(in,ncols+zone_col) = OUTER
          endif

          ! report inconsistent psi and r values
          !if ((lp_data(in,rbin)*lp_data(in,psibin)).lt.0.0.and.lp_data(in,ncols+4).eq.OUTER) then 
          !   write (6,'(a,12(1x,g18.8))') 'INCONSISTENT OUTER:RSEPOUT,PSIN:', (lp_data(cur_line_cnt,it),it=1,ncols)
          !endif

       end do

    endif

    write(0,'(a,g18.8)') 'Finished Flagging Inner/Outer LP Data:',rsplit


  end subroutine flag_lp_data


  subroutine bin_lp_data(lp_axis,lp_axis_psi,lp_axis_r,lp_proc_data,npts,lp_data,nlines,binopt,deltabin,tmin,tmax,chisq_lim,elm_filt,remove_outlier,outlier_mult,n_avs,n_elm_fractions,elm_fractions)
    implicit none
    real,allocatable :: lp_axis_psi(:,:,:), lp_axis_r(:,:,:), lp_proc_data(:,:,:,:) 
    real,allocatable :: lp_axis(:,:)
    real,allocatable :: lp_data(:,:)
    integer n_elm_fractions,n_avs
    real,allocatable :: elm_fractions(:,:)

    integer :: nlines,npts
    integer :: binopt,izone
    real :: deltabin
    real :: tmin,tmax,chisq_lim
    integer :: in,if,nbins, ibin,ierr,testbin,iz
    real :: rmin,rmax,psimin,psimax,testmin
    logical :: elm_filt,remove_outlier

    ! local variables

    integer :: n_av
    real,allocatable :: pre_average(:,:,:,:)

    real :: outlier_mult
    real :: outlier_gt_mult,outlier_lt_mult
    integer :: outlier_cnt

    ! PSI is in lp_data(in,9) and R-Rsep is in LP_DATA(in,8)
    !rbin = 8
    !psibin = 9
    !vfbin = 10
    !
    ! Adding code to calculate the median as well as the bin average
    ! - need median of binned data .. inner/outer .. outliers removed and not removed
    ! 



    write(0,'(a,g12.5)') 'Start Binning LP Data:',outlier_mult


    !outlier_mult = 3.0
    outlier_cnt = 0

    if (outlier_mult.lt.0.0) then 
       outlier_lt_mult = 1.0/abs(outlier_mult)
       outlier_gt_mult = abs(outlier_mult)
    else
       outlier_lt_mult = 0.0
       outlier_gt_mult = outlier_mult
    endif


    ! For data in the range Tmin to Tmax ... bin and average over range of deltar. 
    ! Averaging Te and Jsat
    ! Also bin and average data by inter-ELM fraction


    if (.not.allocated(lp_data)) then 
       call errmsg('BIN_LP_DATA','LP_DATA NOT ALLOCATED')
       stop
    endif

    ! determine minR and maxR values in the desired range. 

    ! set max and min to out of range values
    rmin = 1e25
    rmax = -1e25

    psimin = 1e25
    psimax = -1e25

    !
    ! PSI is in 9, R-Rsep is in 8
    !
    do in = 1,nlines
       ! correct time window
       if (lp_data(in,timebin).ge.tmin.and.lp_data(in,timebin).le.tmax.and.lp_data(in,chibin).le.chisq_lim) then 
          rmin = min(lp_data(in,rbin),rmin)
          rmax = max(lp_data(in,rbin),rmax)
          psimin = min(lp_data(in,psibin),psimin)
          psimax = max(lp_data(in,psibin),psimax)
       endif
    end do

    write(0,*) 'Bins:',binopt,deltabin,rbin,psibin,psimax,psimin

    ! determine number of R bins
    if (binopt.eq.RBINOPT) then 
       nbins = (rmax-rmin)/deltabin + 1
       testmin = rmin
       testbin = rbin
    else
       nbins = (psimax-psimin)/deltabin + 1
       testmin = psimin
       testbin = psibin
    endif


    ! allocate storage for processed data
    ! ndata = 3 ... jsat,te,vf
    !ndata = 3
    npts = nbins

    write(0,*) 'Sizes:',nbins,nzones

    if (allocated(lp_axis)) deallocate(lp_axis)
    allocate(lp_axis(nbins,nzones),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA','ERROR ALLOCATING LP_AXIS')
       stop
    endif
    lp_axis = 0.0

    if (allocated(lp_axis_psi)) deallocate(lp_axis_psi)
    allocate(lp_axis_psi(nbins,n_avs,nzones),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA','ERROR ALLOCATING LP_AXIS_PSI')
       stop
    endif
    lp_axis_psi = 0.0

    if (allocated(lp_axis_r)) deallocate(lp_axis_r)
    allocate(lp_axis_r(nbins,n_avs,nzones),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA','ERROR ALLOCATING LP_AXIS_R')
       stop
    endif
    lp_axis_r = 0.0


    if (allocated(lp_proc_data)) deallocate(lp_proc_data)
    allocate(lp_proc_data(nbins,ndata+1,n_avs,nzones),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA','ERROR ALLOCATING LP_PROC_DATA')
       stop
    endif
    lp_proc_data = 0.0


    if (allocated(pre_average)) deallocate(pre_average)
    allocate(pre_average(nbins,ndata+1,n_avs,nzones),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA','ERROR ALLOCATING PRE_AVERAGE')
       stop
    endif
    pre_average = 0.0

    write(0,*) 'Starting pre-average:'

    ! loop through data and record averages to give a base line for outlier removal
    do in = 1,nlines

        izone = int(lp_data(in,ncols+zone_col))
        if (izone.ne.INNER.and.izone.ne.OUTER) then
          write(0,'(i6,1x,3(1x,g18.8),1x,i8,8(1x,g18.8))') in,lp_data(in,1),lp_data(in,8),lp_data(in,9),int(lp_data(in,6)),&
               &lp_data(in,5),lp_data(in,10),lp_data(in,11),lp_data(in,12),lp_data(in,2),lp_data(in,3),lp_data(in,ncols+4)
        endif 
          
       ! Preaverage each of the possible categories ...
       ! 1. All data
       if (lp_data(in,timebin).ge.tmin.and.lp_data(in,timebin).le.tmax.and.lp_data(in,chibin).le.chisq_lim) then 

          n_av = 1
          ibin = int((lp_data(in,testbin)-testmin)/deltabin) + 1

          ! record bin for data point
          lp_data(in,ncols+bin_col) = ibin

          ! record data count
          pre_average(ibin,ndata+1,n_av,izone) = pre_average(ibin,ndata+1,n_av,izone) + 1.0
          ! record data sum
          ! jsat
          pre_average(ibin,1,n_av,izone) = pre_average(ibin,1,n_av,izone) + lp_data(in,jsatbin)
          ! te
          pre_average(ibin,2,n_av,izone) = pre_average(ibin,2,n_av,izone) + lp_data(in,tebin)

          ! bin over ELM filtered categories
          if (elm_filt) then 
             ! 2,3 ... ELM vs. no-ELM
             ! if ELM reference is greater than 0 then the data point is associated with an ELM
             ! a value less than 0 indicates a point that is in the peripheral region of an ELM 
             if (lp_data(in,ncols+elm_ref_col).gt.0) then 
                ! in - ELM
                n_av = 2
                ! record data count
                pre_average(ibin,ndata+1,n_av,izone) = pre_average(ibin,ndata+1,n_av,izone) + 1.0
                ! record data sum
                ! jsat
                pre_average(ibin,1,n_av,izone) = pre_average(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                ! te
                pre_average(ibin,2,n_av,izone) = pre_average(ibin,2,n_av,izone) + lp_data(in,tebin)
             elseif (lp_data(in,ncols+elm_ref_col).eq.0) then
                ! not in-ELM
                n_av = 3
                ! record data count
                pre_average(ibin,ndata+1,n_av,izone) = pre_average(ibin,ndata+1,n_av,izone) + 1.0
                ! record data sum
                ! jsat
                pre_average(ibin,1,n_av,izone) = pre_average(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                ! te
                pre_average(ibin,2,n_av,izone) = pre_average(ibin,2,n_av,izone) + lp_data(in,tebin)
             endif


             ! 4 .. 4+n_elm_fractions-1

             do if = 1,n_elm_fractions
                n_av = 3 + if

                if (lp_data(in,ncols+elm_frac_col).ge.elm_fractions(if,1).and.lp_data(in,ncols+elm_frac_col).le.elm_fractions(if,2).and.lp_data(in,ncols+elm_ref_col).eq.0) then 
                   ! record data count
                   pre_average(ibin,ndata+1,n_av,izone) = pre_average(ibin,ndata+1,n_av,izone) + 1.0
                   ! record data sum
                   ! jsat
                   pre_average(ibin,1,n_av,izone) = pre_average(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                   ! te
                   pre_average(ibin,2,n_av,izone) = pre_average(ibin,2,n_av,izone) + lp_data(in,tebin)
                endif

             end do
          endif
       endif


    end do

    write(0,*) 'Calculate pre-average:'
    ! calculate pre_average

    do iz = 1,nzones
       do in = 1,nbins
          do if = 1,n_avs
             if (pre_average(in,ndata+1,if,iz).gt.0.0) then 
                pre_average(in,1,if,iz) = pre_average(in,1,if,iz)/pre_average(in,ndata+1,if,iz)
                pre_average(in,2,if,iz) = pre_average(in,2,if,iz)/pre_average(in,ndata+1,if,iz)
             else
                pre_average(in,:,if,iz) = 0.0
             endif
          end do
       end do
    end do


    !write(6,'(a)') 'Pre-averages:'

    !do if = 1,n_avs
    !   write(6,'(a)')
    !   if (if.eq.1) then
    !      write(6,'(a)') ' OVERALL-AVERAGE'
    !   elseif (if.eq.2) then 
    !      write(6,'(a)') ' ELM-PEAK-AVERAGE'
    !   elseif (if.eq.3) then 
    !      write(6,'(a)') ' OUT-OF-ELM-AVERAGE'
    !   else
    !      n_av = if - 3
    !      write(6,'(a,f8.3,a,f8.3)') ' ELM-FRACTION-AVERAGE:',elm_fractions(n_av,1),' TO',elm_fractions(n_av,2)
    !   endif
    !   write(6,'(a)')   ' IBIN      R-Rsep(Bin)     R-Rsep(Av)         PSIN            Jsat(A/cm2)          Te(eV)         Count  '
    !   do in = 1,nbins
    !      if (pre_average(in,ndata+1,if).gt.0.0) then 
    !         write(6,'(i8,20(1x,g18.8))') in,lp_axis(in),lp_axis_r(in,if),lp_axis_psi(in,if),&
    !              pre_average(in,1,if),pre_average(in,2,if),pre_average(in,ndata+1,if)
    !      endif
    !   end do
    !end do


    ! loop through data and record results 
    write(0,*) 'Bin Average LP data:'


    do in = 1,nlines
       ! define zone for averaging
       izone = int(lp_data(in,ncols+zone_col))

       if (lp_data(in,timebin).ge.tmin.and.lp_data(in,timebin).le.tmax.and.lp_data(in,chibin).le.chisq_lim) then 

          !ibin = int((lp_data(in,9)-rmin)/deltar) + 1

          ibin = int((lp_data(in,testbin)-testmin)/deltabin) + 1

          ! record bin for data point
          !lp_data(in,ncols+bin_col) = ibin



          n_av = 1
          if (.not.remove_outlier.or.&
               &(remove_outlier.and.&
               &(((lp_data(in,jsatbin).le.outlier_gt_mult*pre_average(ibin,1,n_av,izone)).and. &
               & (lp_data(in,tebin).le.outlier_gt_mult*pre_average(ibin,2,n_av,izone)))      &
               &.and.((lp_data(in,jsatbin).gt.outlier_lt_mult*pre_average(ibin,1,n_av,izone)).and.&
               & (lp_data(in,tebin).gt.outlier_lt_mult*pre_average(ibin,2,n_av,izone)))&
               &  ))) then 


             ! Bin average for total ... between ELMs and during ELMs ("during" may not be very useful)
             ! record data count
             lp_proc_data(ibin,ndata+1,n_av,izone) = lp_proc_data(ibin,ndata+1,n_av,izone) + 1.0
             ! record data sum
             ! jsat
             lp_proc_data(ibin,1,n_av,izone) = lp_proc_data(ibin,1,n_av,izone) + lp_data(in,jsatbin)
             ! te
             lp_proc_data(ibin,2,n_av,izone) = lp_proc_data(ibin,2,n_av,izone) + lp_data(in,tebin)
             ! te
             lp_proc_data(ibin,3,n_av,izone) = lp_proc_data(ibin,3,n_av,izone) + lp_data(in,vfbin)
             ! Average of axis values psi and r
             lp_axis_psi(ibin,n_av,izone) = lp_axis_psi(ibin,n_av,izone) + lp_data(in,psibin)
             lp_axis_r(ibin,n_av,izone)   = lp_axis_r(ibin,n_av,izone)   + lp_data(in,rbin)



          else
             write(6,'(a,3i8,15(1x,g18.8))') 'OUTLIER removed:',in,ibin,n_av,lp_data(in,1),lp_data(in,2),lp_data(in,3),pre_average(ibin,1,n_av,izone),pre_average(ibin,2,n_av,izone),lp_data(in,12)
             lp_data(in,ncols+outlier_index) = 1.0
             outlier_cnt = outlier_cnt + 1
          endif

          if (elm_filt) then 

             ! 2,3 ... ELM vs. no-ELM
             ! if ELM reference is greater than 0 then the data point is associated with an ELM
             ! a value less than 0 indicates a point that is in the peripheral region of an ELM 
             if (lp_data(in,ncols+elm_ref_col).gt.0) then 
                ! in - ELM
                n_av = 2
                ! Do not apply outlier removal to in-ELM data since wide scatter is to be expected here
                !if (.not.remove_outlier.or.&
                !     &(remove_outlier.and.&
                !     &(((lp_data(in,2).le.outlier_mult*pre_average(ibin,1,n_av)).and.&
                !     & (lp_data(in,3).le.outlier_mult*pre_average(ibin,2,n_av)))&
                !&.and.((lp_data(in,2).ge.pre_average(ibin,1,n_av)/outlier_mult).and.&
                !& (lp_data(in,3).ge.pre_average(ibin,2,n_av)/outlier_mult))&
                !    &  ))) then 


                ! Bin average for total ... between ELMs and during ELMs ("during" may not be very useful)
                ! record data count
                lp_proc_data(ibin,ndata+1,n_av,izone) = lp_proc_data(ibin,ndata+1,n_av,izone) + 1.0
                ! record data sum
                ! jsat
                lp_proc_data(ibin,1,n_av,izone) = lp_proc_data(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                ! te
                lp_proc_data(ibin,2,n_av,izone) = lp_proc_data(ibin,2,n_av,izone) + lp_data(in,tebin)
                ! vf
                lp_proc_data(ibin,3,n_av,izone) = lp_proc_data(ibin,3,n_av,izone) + lp_data(in,vfbin)
                ! Average of axis values psi and r
                lp_axis_psi(ibin,n_av,izone) = lp_axis_psi(ibin,n_av,izone) + lp_data(in,psibin)
                lp_axis_r(ibin,n_av,izone)   = lp_axis_r(ibin,n_av,izone)   + lp_data(in,rbin)

                !else
                !   write(6,'(a,3i8,15(1x,g18.8))') 'OUTLIER removed:',in,ibin,n_av,lp_data(in,1),lp_data(in,2),lp_data(in,3),pre_average(ibin,1,n_av),pre_average(ibin,2,n_av),lp_data(in,12)
                !   outlier_cnt = outlier_cnt + 1
                !endif



             elseif (lp_data(in,ncols+elm_ref_col).eq.0) then 
                ! not associated with an ELM
                n_av = 3
                if (.not.remove_outlier.or.&
                     &(remove_outlier.and.&
                     &(((lp_data(in,jsatbin).le.outlier_gt_mult*pre_average(ibin,1,n_av,izone)).and.&
                     & (lp_data(in,tebin).le.outlier_gt_mult*pre_average(ibin,2,n_av,izone)))&
                     &.and.((lp_data(in,jsatbin).gt.outlier_lt_mult*pre_average(ibin,1,n_av,izone)).and.&
                     & (lp_data(in,tebin).gt.outlier_lt_mult*pre_average(ibin,2,n_av,izone)))&
                     &  ))) then 
                   ! Bin average for total ... between ELMs and during ELMs ("during" may not be very useful)
                   ! record data count
                   lp_proc_data(ibin,ndata+1,n_av,izone) = lp_proc_data(ibin,ndata+1,n_av,izone) + 1.0
                   ! record data sum
                   ! jsat
                   lp_proc_data(ibin,1,n_av,izone) = lp_proc_data(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                   ! te
                   lp_proc_data(ibin,2,n_av,izone) = lp_proc_data(ibin,2,n_av,izone) + lp_data(in,tebin)
                   ! vf
                   lp_proc_data(ibin,3,n_av,izone) = lp_proc_data(ibin,3,n_av,izone) + lp_data(in,vfbin)
                   ! Average of axis values psi and r
                   lp_axis_psi(ibin,n_av,izone) = lp_axis_psi(ibin,n_av,izone) + lp_data(in,psibin)
                   lp_axis_r(ibin,n_av,izone)   = lp_axis_r(ibin,n_av,izone)   + lp_data(in,rbin)

                else
                   write(6,'(a,3i8,15(1x,g18.8))') 'OUTLIER removed:',in,ibin,n_av,lp_data(in,1),lp_data(in,2),lp_data(in,3),pre_average(ibin,1,n_av,izone),pre_average(ibin,2,n_av,izone),lp_data(in,12)
                   lp_data(in,ncols+outlier_index) = 1.0
                   outlier_cnt = outlier_cnt + 1
                endif

             endif

             ! 4 .. 4+n_elm_fractions-1
             do if = 1,n_elm_fractions
                n_av = 3 + if

                if (lp_data(in,ncols+elm_frac_col).ge.elm_fractions(if,1).and.lp_data(in,ncols+elm_frac_col).le.elm_fractions(if,2).and.lp_data(in,ncols+elm_ref_col).eq.0) then 


                   if (.not.remove_outlier.or.&
                        &(remove_outlier.and.&
                        &(((lp_data(in,jsatbin).le.outlier_gt_mult*pre_average(ibin,1,n_av,izone)).and.&
                        & (lp_data(in,tebin).le.outlier_gt_mult*pre_average(ibin,2,n_av,izone)))&
                        &.and.((lp_data(in,jsatbin).gt.outlier_lt_mult*pre_average(ibin,1,n_av,izone)).and.&
                        & (lp_data(in,tebin).gt.outlier_lt_mult*pre_average(ibin,2,n_av,izone)))&
                        &  ))) then 


                      ! Bin average for total ... between ELMs and during ELMs ("during" may not be very useful)
                      ! record data count
                      lp_proc_data(ibin,ndata+1,n_av,izone) = lp_proc_data(ibin,ndata+1,n_av,izone) + 1.0
                      ! record data sum
                      ! jsat
                      lp_proc_data(ibin,1,n_av,izone) = lp_proc_data(ibin,1,n_av,izone) + lp_data(in,jsatbin)
                      ! te
                      lp_proc_data(ibin,2,n_av,izone) = lp_proc_data(ibin,2,n_av,izone) + lp_data(in,tebin)
                      ! vf
                      lp_proc_data(ibin,3,n_av,izone) = lp_proc_data(ibin,3,n_av,izone) + lp_data(in,vfbin)
                      ! Average of axis values psi and r
                      lp_axis_psi(ibin,n_av,izone) = lp_axis_psi(ibin,n_av,izone) + lp_data(in,psibin)
                      lp_axis_r(ibin,n_av,izone)   = lp_axis_r(ibin,n_av,izone)   + lp_data(in,rbin)

                   else
                      write(6,'(a,3i8,15(1x,g18.8))') 'OUTLIER removed:',in,ibin,n_av,lp_data(in,1),lp_data(in,2),lp_data(in,3),pre_average(ibin,1,n_av,izone),pre_average(ibin,2,n_av,izone),lp_data(in,12)
                      lp_data(in,ncols+outlier_index) = 1.0
                      outlier_cnt = outlier_cnt + 1
                   endif


                endif

             end do

          endif


       endif

    end do

    write(0,*) 'Finished ELM filtering:'

    ! loop through - calculate averages and assign axis values
    do iz = 1,nzones
       do in = 1,nbins
          do if = 1,n_avs
             if (lp_proc_data(in,ndata+1,if,iz).gt.0.0) then 
                lp_proc_data(in,1,if,iz) = lp_proc_data(in,1,if,iz)/lp_proc_data(in,ndata+1,if,iz)
                lp_proc_data(in,2,if,iz) = lp_proc_data(in,2,if,iz)/lp_proc_data(in,ndata+1,if,iz)
                lp_proc_data(in,3,if,iz) = lp_proc_data(in,3,if,iz)/lp_proc_data(in,ndata+1,if,iz)
                lp_axis_psi(in,if,iz) = lp_axis_psi(in,if,iz)/lp_proc_data(in,ndata+1,if,iz)
                lp_axis_r(in,if,iz) = lp_axis_r(in,if,iz)/lp_proc_data(in,ndata+1,if,iz)
             else
                lp_proc_data(in,:,if,iz) = 0.0
             endif
          end do

          ! assign axis value
          lp_axis(in,iz) = (real(in) - 0.5)* deltabin + testmin 

       end do

    end do

    write(0,'(a)') 'Finished Binning LP Data:'


  end subroutine bin_lp_data


  subroutine print_lp_bin_data(ounit,lp_axis,lp_axis_r,lp_axis_psi,lp_proc_data,npts,ident,n_avs,elm_filt,n_elm_fractions,elm_fractions,targ_flag)
    implicit none

    real,allocatable :: lp_axis(:,:),lp_axis_r(:,:,:),lp_axis_psi(:,:,:),lp_proc_data(:,:,:,:),elm_fractions(:,:)
    integer :: n_avs,n_elm_fractions,targ_flag
    integer :: npts,ounit,in,if
    character*(*) :: ident
    logical elm_filt

    integer :: irange
    integer :: iz
    integer :: start_zone,end_zone

    if (.not.allocated(lp_axis)) then 
       call errmsg('PRINT_LP_BIN_DATA','LP_AXIS NOT ALLOCATED') 
       stop
    endif

    if (targ_flag.eq.inner) then 
       start_zone=inner
       end_zone=inner
    elseif(targ_flag.eq.outer) then 
       start_zone=outer
       end_zone=outer
    else
       start_zone=1
       end_zone = nzones
    endif


    ! print out the processed/binned data
    ! at the present time this data is  LP_AXIS   LP_PROC_DATA(,1) = Jsat   LP_PROC_DATA(,2) = Te    LP_PROC_DATA(,3) = floating potential

    write(ounit,'(a,a)') 'ID : ',trim(ident)

    do if = 1,n_avs
       write(ounit,'(a)')
       if (if.eq.1) then
          write(ounit,'(a)') ' OVERALL-AVERAGE'
       elseif (if.eq.2) then 
          write(ounit,'(a)') ' ELM-PEAK-AVERAGE'
       elseif (if.eq.3) then 
          write(ounit,'(a)') ' OUT-OF-ELM-AVERAGE'
       else
          irange = if - 3
          write(ounit,'(a,f8.3,a,f8.3)') ' ELM-FRACTION-AVERAGE:',elm_fractions(irange,1),' TO',elm_fractions(irange,2)
       endif


       do iz = start_zone,end_zone
          if (iz.eq.inner) then 
             write(ounit,'(a)')   'Inner   BIN_coord     R-Rsep(Av)         Psin(Av)            Jsat(A/cm2)          Te(eV)      Vf(V)      Count  '
          elseif (iz.eq.outer) then 
             write(ounit,'(a)')   'Outer   BIN_coord     R-Rsep(Av)         Psin(Av)            Jsat(A/cm2)          Te(eV)      Vf(V)      Count  '
          else
             write(ounit,'(a)')   '        BIN_coord     R-Rsep(Av)         Psin(Av)            Jsat(A/cm2)          Te(eV)      Vf(V)      Count  '
          endif

          do in = 1,npts
             if (lp_proc_data(in,ndata+1,if,iz).gt.0.0) then 
                write(ounit,'(20(1x,g18.8))') lp_axis(in,iz),lp_axis_r(in,if,iz),lp_axis_psi(in,if,iz),&
                     lp_proc_data(in,1,if,iz),lp_proc_data(in,2,if,iz),lp_proc_data(in,3,if,iz),lp_proc_data(in,ndata+1,if,iz)
             endif
          end do
       end do
    end do

  end subroutine print_lp_bin_data


  subroutine print_lp_data(ounit,lp_data,nlines,ident,tmin,tmax,chisq_lim)
    implicit none

    real,allocatable :: lp_data(:,:)
    real :: tmin,tmax,chisq_lim
    integer :: nlines
    integer :: ounit,in
    character*(*) :: ident

    if (.not.allocated(lp_data)) then 
       call errmsg('PRINT_LP_DATA','LP_DATA NOT ALLOCATED') 
       stop
    endif

    ! print out the revised data with ELM offset added

    write(ounit,'(a,a)') ' ID : ',trim(ident)
    write(ounit,'(a)')   ' IN            Time(s)              R-Rsep             PSIN           Probe'//&
         &'       CHISQ     ELM-time              ELM-Index        ELM-frac        Jsat(A/cm2)              Te(eV)    Vf(V)    TARG_FLAG    OUT_FLAG'

    do in = 1,nlines
       if (lp_data(in,chibin).lt.chisq_lim.and.(lp_data(in,timebin).ge.tmin.and.lp_data(in,timebin).le.tmax).and.(lp_data(in,ncols+outlier_index).eq.0.0)) then 

          write(ounit,'(i6,1x,3(1x,g18.8),1x,i8,7(1x,g18.8),2i10)') in,lp_data(in,timebin),lp_data(in,rbin),lp_data(in,psibin),int(lp_data(in,probebin)),&
               &lp_data(in,chibin),lp_data(in,ncols+elm_time_col),lp_data(in,ncols+elm_ref_col),lp_data(in,ncols+elm_frac_col),&
               &lp_data(in,jsatbin),lp_data(in,tebin),lp_data(in,vfbin),int(lp_data(in,ncols+zone_col)),int(lp_data(in,ncols+outlier_index))
       endif
    end do


    write(ounit,'(a,a)') 
    write(ounit,'(a,a,a)') ' ID : ',trim(ident),': OUTLIERS :'
    write(ounit,'(a)')   ' IN            Time(s)              R-Rsep             PSIN           Probe'//&
         &'       CHISQ     ELM-time              ELM-Index        ELM-frac        Jsat(A/cm2)              Te(eV)     Vf(V)    TARG_FLAG   OUT_FLAG'
    do in = 1,nlines

       if (lp_data(in,chibin).lt.chisq_lim.and.(lp_data(in,timebin).ge.tmin.and.lp_data(in,timebin).le.tmax).and.(lp_data(in,ncols+outlier_index).ne.0.0)) then 
          write(ounit,'(i6,1x,3(1x,g18.8),1x,i8,7(1x,g18.8),2i10)') in,lp_data(in,timebin),lp_data(in,rbin),lp_data(in,psibin),int(lp_data(in,probebin)),&
               &lp_data(in,chibin),lp_data(in,ncols+elm_time_col),lp_data(in,ncols+elm_ref_col),lp_data(in,ncols+elm_frac_col),&
               &lp_data(in,jsatbin),lp_data(in,tebin),lp_data(in,vfbin),int(lp_data(in,ncols+zone_col)),int(lp_data(in,ncols+outlier_index))
       endif
    end do



  end subroutine print_lp_data




  subroutine filter_lp_data(elm_filename,lp_data,nlines,elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)
    use filter_elms
    implicit none
    integer :: nlines
    character*(*) elm_filename
    real,allocatable :: lp_data(:,:)
    !real lp_data(nlines,ncols+nextra)
    integer :: ierr
    integer :: in
    real :: elm_time_offset
    integer :: elmref
    logical :: inelm

    !
    ! Load ELM time data from file and then add ELM information to the LP data 
    !

    real :: elm_start_input,elm_end_input, elm_effect_start_input, elm_effect_end_input
    real :: elmfrac


    call load_elm_times(elm_filename,ierr)

    if (ierr.ne.0) return

    call set_elm_criteria(elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)

    ! An elmref value of 0 means that it is outside time window that is associated with an ELM
    ! A negative elmref value is inside a possible ELM effect window
    ! A positive elmref value is inside the most likely ELM effect window
    do in = 1, nlines

       call get_elm_time(lp_data(in,timebin),elm_time_offset,elmref,inelm,elmfrac,0)

       lp_data(in,ncols+elm_time_col) = elm_time_offset
       lp_data(in,ncols+elm_ref_col) = 0

       if (inelm) then 
          lp_data(in,ncols+elm_ref_col) = elmref
       else
          lp_data(in,ncols+elm_ref_col) = -elmref
       endif

       lp_data(in,ncols+elm_frac_col) = elmfrac

    end do


  end subroutine filter_lp_data


  subroutine setup_time_filter(time_filt,time_filename)
    use utilities
    implicit none
    logical :: time_filt
    character*(*) :: time_filename
    character*1024 :: line
    integer :: iunit,ierr,line_cnt,in,ios
    real :: t1,t2


    filter_times = .true.

    ios = 0
    ierr = 0

    call find_free_unit_number(iunit)

    open(unit=iunit,file=trim(time_filename),status='old',iostat=ierr)
    if (ierr.ne.0) then 
       call errmsg('setup_time_filter: problem opening time data file:'//trim(time_filename),ierr)
       filter_times = .false.
       return
    endif

    ! 
    ! Read in time data
    !
    ! first line is shot number followed by any number of time slices
    !


    line_cnt = 0

    do while (ios.eq.0) 

       read(iunit,line_form,iostat=ios) line
       !write(0,*) 'tw:',trim(line)
       if (ios.eq.0) then 
          read(line,*,iostat=ierr) t1,t2
          if (ierr.eq.0) then 
             line_cnt = line_cnt + 1
             !write(0,*) 'twa:',line_cnt,t1,t2
           endif
       endif
          
    end do
    
    if (allocated(time_windows)) deallocate(time_windows)
    allocate(time_windows(line_cnt,2),stat=ierr)
    time_window_cnt = line_cnt

    rewind(iunit)

    line_cnt = 0
    ios = 0
    ierr = 0

    do while (ios.eq.0) 

       read(iunit,line_form,iostat=ios) line
       if (ios.eq.0) then 
          read(line,*,iostat=ierr) t1,t2
          if (ierr.eq.0) then 
             line_cnt = line_cnt + 1
             !write(0,*) 'twb:',line_cnt,t1,t2
             time_windows(line_cnt,1) = t1
             time_windows(line_cnt,2) = t2
           endif
       endif
    end do

    close(iunit)


    write(0,*) 'Time Window filtering active:'
    do in=1,time_window_cnt
       write(0,'(a,i8,2(1x,g12.5))') 'Window: ', in,time_windows(in,1),time_windows(in,2)
    end do


  end subroutine setup_time_filter



  logical function check_time_filter(time)
    implicit none
    real :: time
    integer :: in
    
    if (filter_times.eq. .false.) then 
       check_time_filter = .true. 
       return
    endif

    check_time_filter = .false. 

    do in = 1,time_window_cnt
       if (time.ge.time_windows(in,1).and.time.le.time_windows(in,2)) then 
          check_time_filter = .true.
          return
       endif
    end do 

  end function check_time_filter
    



end module langmuir_probes
