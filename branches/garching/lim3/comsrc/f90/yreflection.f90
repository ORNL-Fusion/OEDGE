module yreflection
  use error_handling
  implicit none

  integer :: yreflection_opt
  real :: cmir_refl_lower, cmir_refl_upper
  real*8 :: cmir_refl_lower_dp, cmir_refl_upper_dp
  real*8 :: yreflection_event_count = 0
  real*8 :: relocation_count = 0

  real*8 :: lim_sep,deltay
  real*8 :: minrefl

  real*8 :: yref_upper_cnt = 0.0
  real*8 :: yref_lower_cnt = 0.0

  private :: lim_sep,cmir_refl_lower_dp,cmir_refl_upper_dp,minrefl
  save

contains

  subroutine init_reflection(ctwol)
    implicit none
    real :: ctwol

    ! initialize data
    yreflection_event_count = 0.0
    relocation_count = 0.0

    yref_upper_cnt = 0.0
    yref_lower_cnt = 0.0

    lim_sep = ctwol
    
    cmir_refl_lower_dp = cmir_refl_lower
    cmir_refl_upper_dp = cmir_refl_upper

    minrefl = 1.0d-7

  end subroutine init_reflection

  logical function check_reflected_region(y)
    implicit none
    real  ::  y,absy
    real*8 :: ynew

    ! check to see if the input y value lies in the region with mirrors on either side 
    ! If it does then the particle can't be launched here or does not belong here and 
    ! has leaked somehow.

    ynew = y
    
    check_reflected_region=.false.

    if ((ynew.ge.-lim_sep+cmir_refl_upper_dp.and.ynew.le.cmir_refl_lower_dp).or.&
        (ynew.ge.cmir_refl_upper_dp.and.ynew.le.lim_sep+cmir_refl_lower_dp)) then
       check_reflected_region=.true.
    endif


  end function check_reflected_region


  subroutine check_reflection(x,y,oldy,svy,debugl,ierr)
    implicit none
    real :: x,y,svy
    real,intent(in) :: oldy
    logical :: debugl
    integer :: ierr

    logical :: reflected
    real :: y_org,oldy_org,svy_org
    real :: y_tmp
    real*8 :: ynew,yprev


    real :: ran
    real,external :: getranf

    if (yreflection_opt.ne.0.and.(cmir_refl_lower.ge.0.0.or.cmir_refl_upper.le.0.0)) then 
       call errmsg('YREFLECTION:TEST_REFLECTION','Y-REFLECTION OPTION IS ON BUT REFLECTION&
                   & LOCATIONS HAVE NOT BEEN PROPERLY SPECIFIED : Y-REFLECTION DISABLED')
       yreflection_opt=0
       return
    endif

    ierr =0
    reflected = .false.

    if (debugl) then 
       y_org = y
       oldy_org = oldy
       svy_org = svy
    endif

    !
    ! jdemod - perform intermediate calculations in double precision
    !

    yprev = dble(oldy)
    ynew  = dble(y)

    if (ynew.lt.0) then 
       if ((yprev.gt.cmir_refl_lower_dp).and.(ynew.le.cmir_refl_lower_dp)) then 
          !
          !                   Reflection at lower (<0) mirror
          !

          deltay = (cmir_refl_lower_dp-ynew)
          ynew = cmir_refl_lower_dp + max(deltay,minrefl)

          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.
          
          yref_lower_cnt = yref_lower_cnt + 1.0

       elseif ((yprev.lt.(-lim_sep+cmir_refl_upper_dp)).and.(ynew.ge.(-lim_sep+cmir_refl_upper_dp))) then 
          !
          !                   Reflection at second upper (< 0) mirror
          !

          deltay = (-lim_sep+cmir_refl_upper_dp-ynew)
          ynew = -lim_sep+cmir_refl_upper_dp +  min(deltay,-minrefl)

          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

          yref_upper_cnt = yref_upper_cnt + 1.0

       endif

    elseif (ynew.ge.0) then 
       if ((yprev.lt.cmir_refl_upper_dp).and.(ynew.ge.cmir_refl_upper_dp)) then 
          !
          !                   Reflection at upper (>0) mirror
          !
          deltay = (cmir_refl_upper_dp - ynew)
          ynew = cmir_refl_upper_dp + min(deltay,-minrefl)

          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

          yref_upper_cnt = yref_upper_cnt + 1.0


       elseif ((yprev.gt.(lim_sep+cmir_refl_lower_dp)).and.(ynew.le.(lim_sep+cmir_refl_lower_dp))) then 
          !
          !                   Reflection at second lower (>0) mirror
          !
          deltay = (lim_sep +cmir_refl_lower_dp-ynew)
          ynew = lim_sep + cmir_refl_lower_dp + max(deltay,minrefl)

          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

          yref_lower_cnt = yref_lower_cnt + 1.0

       endif
    endif

    y = sngl(ynew)

    if (debugl.and.reflected) then
       write(error_message_data,'(a,i10,6(1x,g18.10))') 'REFLECTION OCCURRED:',int(yreflection_event_count),y,oldy,svy,svy_org
       call dbgmsg('CHECK_REFLECTION',error_message_data)
    endif

    if (check_reflected_region(y)) then 
       
       write(error_message_data,'(a,4(1x,g18.10))') 'REFLECTED PARTICLE HAS ENTERED MIRROR REGION - '//&
                                                  & 'TRY REDUCING SIMULATION TIMESTEPS AND MULTIPLIERS : DATA:', x,y,oldy,svy
       call errmsg('CHECK REFLECTION:WARNING:',error_message_data)
       
       y_tmp = y
       ran = getranf()
       y = cmir_refl_lower + ran * (cmir_refl_upper-cmir_refl_lower)

       write(error_message_data,'(a,i10,l4,10(1x,g20.12))') 'REFLECTION:',int(yreflection_event_count),reflected,x,deltay,y_org,ynew,oldy_org,oldy,ynew,yprev,cmir_refl_lower_dp,cmir_refl_upper_dp
       call dbgmsg('CHECK_REFLECTION',error_message_data)

       relocation_count = relocation_count +1.0

       write(error_message_data,'(a,2i12,10(1x,g20.12))') 'REVISED Y COORDINATES:',int(relocation_count),int(yreflection_event_count), y_org,oldy,y_tmp,y,svy_org,svy
       call errmsg('CHECK_REFLECTION PARTICLE RELOCATED:',error_message_data)


       ierr =1
    endif


  end subroutine check_reflection



  subroutine pr_yref_stats
    implicit none

!
! Print some statistics on yreflection events
! 
    if (yreflection_opt.ne.0) then 

        call prb
        call prc('  Summary of Y-Reflection events:')
        call prr('  Total number of reflection events:',sngl(yreflection_event_count))
        call prr('  Lower Mirror Location: ',cmir_refl_lower)
        call prr('     - Total Reflections from lower mirror:',sngl(yref_lower_cnt))
        call prr('  Upper Mirror Location: ',cmir_refl_upper)
        call prr('     - Total Reflections from lower mirror:',sngl(yref_upper_cnt))
        call prr('  Total number of relocated particles: ',sngl(relocation_count))

    endif

  end subroutine pr_yref_stats




end module yreflection
