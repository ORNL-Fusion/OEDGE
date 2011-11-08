module langmuir_probes
  use error_handling
  implicit none


contains

  subroutine read_lp_data_file(iunit,lp_data,nlines,ncols,nextra)
    implicit none
    integer :: iunit,nlines,ncols,nextra
    real, allocatable :: lp_data(:,:)
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
    allocate(lp_data(line_cnt,ncols+nextra),stat=ierr)
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
    ! old version 1 : 9 columns
    !   time(msec)	jsat(A/cm2)	temp(eV)	dens(cm-3)  		csq		probeID		delrsepin	delzsepin	delrsepout
    ! new version : 9 columns
    !   time(msec)	jsat(A/acm2)	temp(eV)	dens(cm-3)              csq	        probeID	        delrsepin	delrsepout	psin

    do in = 1,line_cnt
       read(iunit,*,iostat=ios) (lp_data(in,it),it=1,ncols)


       if (lp_data(in,8).gt.0.0.and.lp_data(in,9).lt.1.0) then 
          write (6,'(10(1x,g18.8))') (lp_data(in,it),it=1,ncols)
       endif


    end do

    nlines = line_cnt




  end subroutine read_lp_data_file


  subroutine bin_lp_data_r(lp_axis,lp_axis_psi,lp_proc_data,npts,ndata,lp_data,nlines,ncols,nextra,deltar,tmin,tmax,chisq_lim,elm_filt,remove_outlier)
    implicit none
    real,allocatable :: lp_axis_psi(:,:), lp_proc_data(:,:,:) 
    real,allocatable :: lp_axis(:)
    real,allocatable :: lp_data(:,:)

    integer :: nlines, ncols,ndata,npts,nextra
    real :: deltar
    real :: tmin,tmax,chisq_lim
    integer :: in,nrbins, ibin,ierr
    real :: rmin,rmax
    logical :: elm_filt,remove_outlier

    ! local variables

    integer :: rbin,psibin
    real,allocatable :: pre_average(:,:)
    real :: outlier_mult
    integer :: outlier_cnt

    ! PSI is in lp_data(in,9) and R-Rsep is in LP_DATA(in,8)
    rbin = 8
    psibin = 9

    outlier_mult = 5.0
    outlier_cnt = 0

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

    !
    ! PSI is in 9, R-Rsep is in 8
    !
    do in = 1,nlines
       ! correct time window
       if (lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax.and.lp_data(in,5).le.chisq_lim) then 
          !rmin = min(lp_data(in,9),rmin)
          !rmax = max(lp_data(in,9),rmax)
          rmin = min(lp_data(in,rbin),rmin)
          rmax = max(lp_data(in,rbin),rmax)
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

    if (allocated(lp_axis_psi)) deallocate(lp_axis_psi)
    allocate(lp_axis_psi(nrbins,3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA_R','ERROR ALLOCATING LP_AXIS_PSI')
       stop
    endif
    lp_axis_psi = 0.0


    if (allocated(lp_proc_data)) deallocate(lp_proc_data)
    allocate(lp_proc_data(nrbins,ndata+1,3),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA_R','ERROR ALLOCATING LP_PROC_DATA')
       stop
    endif
    lp_proc_data = 0.0


    if (allocated(pre_average)) deallocate(pre_average)
    allocate(pre_average(nrbins,ndata+1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('BIN_LP_DATA_R','ERROR ALLOCATING PRE_AVERAGE')
       stop
    endif
    pre_average = 0.0




    ! loop through data and record averages to give a base line for outlier removal

    do in = 1,nlines
       if (lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax.and.lp_data(in,5).le.chisq_lim) then 
          !ibin = int((lp_data(in,9)-rmin)/deltar) + 1
          ibin = int((lp_data(in,rbin)-rmin)/deltar) + 1

          ! record data count
          pre_average(ibin,ndata+1) = pre_average(ibin,ndata+1) + 1.0
          ! record data sum
          ! jsat
          pre_average(ibin,1) = pre_average(ibin,1) + lp_data(in,2)
          ! te
          pre_average(ibin,2) = pre_average(ibin,2) + lp_data(in,3)

       endif
    end do

    ! calculate pre_average

    do in = 1,nrbins
       if (pre_average(in,ndata+1).gt.0.0) then 
          pre_average(in,1) = pre_average(in,1)/pre_average(in,ndata+1)
          pre_average(in,2) = pre_average(in,2)/pre_average(in,ndata+1)
       else
          pre_average(in,:) = 0.0
       endif
    end do



    ! loop through data and record results - use lp_axis to record counts on first run through

    do in = 1,nlines
       if (lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax.and.lp_data(in,5).le.chisq_lim) then 

          !ibin = int((lp_data(in,9)-rmin)/deltar) + 1

          ibin = int((lp_data(in,rbin)-rmin)/deltar) + 1


          ! Bin average for total ... between ELMs and during ELMs ("during" is not very useful)

          ! record data count
          lp_proc_data(ibin,ndata+1,1) = lp_proc_data(ibin,ndata+1,1) + 1.0
          ! record data sum
          ! jsat
          lp_proc_data(ibin,1,1) = lp_proc_data(ibin,1,1) + lp_data(in,2)
          ! te
          lp_proc_data(ibin,2,1) = lp_proc_data(ibin,2,1) + lp_data(in,3)

          lp_axis_psi(ibin,1) = lp_axis_psi(ibin,1) + lp_data(in,psibin)

          if (elm_filt) then 
             ! add data around 1ms offset to ELM average
             if (lp_data(in,10).ge.0.5.and.lp_data(in,10).le.1.5) then 
                ! Associated with ELM 
                ! record data count
                lp_proc_data(ibin,ndata+1,2) = lp_proc_data(ibin,ndata+1,2) + 1.0
                ! record data sum
                ! jsat
                lp_proc_data(ibin,1,2) = lp_proc_data(ibin,1,2) + lp_data(in,2)
                ! te
                lp_proc_data(ibin,2,2) = lp_proc_data(ibin,2,2) + lp_data(in,3)

                lp_axis_psi(ibin,2) = lp_axis_psi(ibin,2) + lp_data(in,psibin)

             elseif (lp_data(in,11).eq.0.0) then 
                ! Between ELM

                if (.not.remove_outlier.or.&
                     &(remove_outlier.and.&
                     &(((lp_data(in,2).lt.outlier_mult*pre_average(ibin,1)).and.&
                     & (lp_data(in,3).lt.outlier_mult*pre_average(ibin,2)))&
                     !&.and.((lp_data(in,2).gt.pre_average(ibin,1)/outlier_mult).and.&
                     !& (lp_data(in,3).gt.pre_average(ibin,2)/outlier_mult))&
                     &  ))) then 



                   ! record data count
                   lp_proc_data(ibin,ndata+1,3) = lp_proc_data(ibin,ndata+1,3) + 1.0
                   ! record data sum
                   ! jsat
                   lp_proc_data(ibin,1,3) = lp_proc_data(ibin,1,3) + lp_data(in,2)
                   ! te
                   lp_proc_data(ibin,2,3) = lp_proc_data(ibin,2,3) + lp_data(in,3)

                   lp_axis_psi(ibin,3) = lp_axis_psi(ibin,3) + lp_data(in,psibin)

                else
                   write(0,'(a,2i8,10(1x,g18.8))') 'OUTLIER removed:',in,ibin,lp_data(in,1),lp_data(in,2),lp_data(in,3),pre_average(ibin,1),pre_average(ibin,2)
                   outlier_cnt = outlier_cnt + 1
                endif


             endif

          endif
       endif
    end do

    ! loop through - calculate averages and assign axis values

    do in = 1,nrbins
       if (lp_proc_data(in,ndata+1,1).gt.0.0) then 
          lp_proc_data(in,1,1) = lp_proc_data(in,1,1)/lp_proc_data(in,ndata+1,1)
          lp_proc_data(in,2,1) = lp_proc_data(in,2,1)/lp_proc_data(in,ndata+1,1)
          lp_axis_psi(in,1) = lp_axis_psi(in,1)/lp_proc_data(in,ndata+1,1)
       else
          lp_proc_data(in,:,1) = 0.0
       endif

       if (elm_filt) then 

          ! during ELM
          if (lp_proc_data(in,ndata+1,2).gt.0.0) then 
             lp_proc_data(in,1,2) = lp_proc_data(in,1,2)/lp_proc_data(in,ndata+1,2)
             lp_proc_data(in,2,2) = lp_proc_data(in,2,2)/lp_proc_data(in,ndata+1,2)
             lp_axis_psi(in,2) = lp_axis_psi(in,2)/lp_proc_data(in,ndata+1,2)
          else
             lp_proc_data(in,:,2) = 0.0
          endif

          ! between ELM
          if (lp_proc_data(in,ndata+1,3).gt.0.0) then 
             lp_proc_data(in,1,3) = lp_proc_data(in,1,3)/lp_proc_data(in,ndata+1,3)
             lp_proc_data(in,2,3) = lp_proc_data(in,2,3)/lp_proc_data(in,ndata+1,3)
             lp_axis_psi(in,3) = lp_axis_psi(in,3)/lp_proc_data(in,ndata+1,3)
          else
             lp_proc_data(in,:,3) = 0.0
          endif

       endif

       ! assign axis value
       lp_axis(in) = (real(in) - 0.5)* deltar + rmin 

    end do


  end subroutine bin_lp_data_r


  subroutine print_lp_bin_data(ounit,lp_axis,lp_axis_psi,lp_proc_data,npts,ndata,ident)
    implicit none
    
    real,allocatable :: lp_axis(:),lp_axis_psi(:,:),lp_proc_data(:,:,:)
    integer :: npts,ndata,ounit,in
    character*(*) :: ident

    if (.not.allocated(lp_axis)) then 
       call errmsg('PRINT_LP_BIN_DATA','LP_AXIS NOT ALLOCATED') 
       stop
    endif

    ! print out the processed/binned data
    ! at the present time this data is  LP_AXIS   LP_PROC_DATA(,1) = Jsat   LP_PROC_DATA(,2) = Te

    write(ounit,'(a,a)') 'ID : ',trim(ident)
    write(ounit,'(a)')   '       R-Rsep        PSIN            Jsat(A/cm2)          Te(eV)         Count  '//&
                              '       PSI(ELM)        Jsat(ELM)            Te(ELM)        Count(ELM)  '//&
                              '       PSI(NOE)        Jsat(NOE)            Te(NOE)        Count(NOE)'

    do in = 1,npts
       if (lp_proc_data(in,1,1).gt.0.0) then 
          write(ounit,'(20(1x,g18.8))') lp_axis(in),lp_axis_psi(in,1),&
                                                    lp_proc_data(in,1,1),lp_proc_data(in,2,1),lp_proc_data(in,ndata+1,1),&
                                                    lp_axis_psi(in,2),&  
                                                    lp_proc_data(in,1,2),lp_proc_data(in,2,2),lp_proc_data(in,ndata+1,2),&
                                                    lp_axis_psi(in,3),&
                                                    lp_proc_data(in,1,3),lp_proc_data(in,2,3),lp_proc_data(in,ndata+1,3)
       endif
    end do 


  end subroutine print_lp_bin_data


  subroutine print_lp_data(ounit,lp_data,nlines,ncols,nextra,ident,tmin,tmax,chisq_lim)
    implicit none
    
    real,allocatable :: lp_data(:,:)
    real :: tmin,tmax,chisq_lim
    integer :: nlines,ncols,nextra
    integer :: ounit,in
    character*(*) :: ident

    if (.not.allocated(lp_data)) then 
       call errmsg('PRINT_LP_DATA','LP_DATA NOT ALLOCATED') 
       stop
    endif

    ! print out the revised data with ELM offset added

    write(ounit,'(a,a)') 'ID : ',trim(ident)
    write(ounit,'(a)')   'IN            Time(s)              R-Rsep             PSIN           Probe'//&
                        &'       CHISQ     ELM-time              ELM-Index              Jsat(A/cm2)              Te(eV)'

    do in = 1,nlines
       if (lp_data(in,5).lt.chisq_lim.and.(lp_data(in,1).ge.tmin.and.lp_data(in,1).le.tmax)) then 
          write(ounit,'(i6,1x,3(1x,g18.8),1x,i8,10(1x,g18.8))') in,lp_data(in,1),lp_data(in,8),lp_data(in,9),int(lp_data(in,6)),&
                &lp_data(in,5),lp_data(in,10),lp_data(in,11),lp_data(in,2),lp_data(in,3)
       endif
    end do 


  end subroutine print_lp_data




  subroutine filter_lp_data(elm_filename,lp_data,nlines,ncols,nextra)
    use filter_elms
    implicit none
    integer :: nlines,ncols,nextra
    character*(*) elm_filename
    real,allocatable :: lp_data(:,:)
    !real lp_data(nlines,ncols+nextra)
    integer :: ierr
    integer :: in
    real :: elm_time_offset
    integer :: elmref
    logical :: inelm

    !sd8
    ! Load ELM time data from file and then add ELM information to the LP data 
    !
    
    real :: elm_start_input,elm_end_input, elm_effect_start_input, elm_effect_end_input


    ! try elm time window of -2 to +4
    elm_start_input = -2.0
    elm_end_input   =  4.0
    elm_effect_start_input = -2.0
    elm_effect_end_input   =  4.0
    
    call load_elm_times(elm_filename,ierr)

    if (ierr.ne.0) return


    call set_elm_criteria(elm_start_input,elm_end_input,elm_effect_start_input,elm_effect_end_input)

    ! An elmref value of 0 means that it is outside time window that is associated with an ELM
    ! A negative elmref value is inside a possible ELM effect window
    ! A positive elmref value is inside the most likely ELM effect window
    do in = 1, nlines

       call get_elm_time(lp_data(in,1),elm_time_offset,elmref,inelm)

       lp_data(in,ncols+1) = elm_time_offset
       lp_data(in,ncols+2) = 0
       if (inelm) then 
          lp_data(in,ncols+2) = elmref
       else
          lp_data(in,ncols+2) = -elmref
       endif

    end do


  end subroutine filter_lp_data



end module langmuir_probes
