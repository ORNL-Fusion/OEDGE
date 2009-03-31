module yreflection
  use error_handling
  implicit none

  integer :: yreflection_opt
  real :: cmir_refl_lower, cmir_refl_upper
  real :: yreflection_event_count = 0
  real :: relocation_count = 0
  real :: lim_sep
  save

contains

  subroutine init_reflection(ctwol)
    implicit none
    real :: ctwol

    ! initialize data
    yreflection_event_count = 0.0
    relocation_count = 0.0
    lim_sep = ctwol

  end subroutine init_reflection

  logical function check_reflected_region(y)
    implicit none
    real  ::  y,absy

    ! check to see if the input y value lies in the region with mirrors on either side 
    ! If it does then the particle can't be launched here or does not belong here and 
    ! has leaked somehow.
    
    absy = abs(y)
    check_reflected_region=.false.

    if (y.ge.cmir_refl_upper.and.y.le.-cmir_refl_lower) then
       check_reflected_region=.true.
    endif

  end function check_reflected_region


  subroutine check_reflection(y,oldy,svy,debugl,ierr)
    implicit none
    real :: y,oldy,svy
    logical :: debugl
    integer :: ierr

    logical :: reflected
    real :: y_org,oldy_org,svy_org
    real :: y_tmp

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
!    if (debugl) then 
       y_org = y
       oldy_org = oldy
       svy_org = svy
!    endif

    if (y.lt.0) then 
       if ((oldy.gt.cmir_refl_lower).and.(y.le.cmir_refl_lower)) then 
          !
          !                   Reflection at lower (<0) mirror
          !
          ! cmir_refl_lower is <0

          y = cmir_refl_lower + abs(abs(y)+cmir_refl_lower)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

       elseif ((oldy.lt.(-lim_sep+cmir_refl_upper)).and.(y.ge.(-lim_sep+cmir_refl_upper))) then 
          !
          !                   Reflection at second upper (< 0) mirror
          !
          ! cmir_refl_upper is >0
          
          y = -lim_sep+cmir_refl_upper - abs(abs(y)-cmir_refl_upper)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

       endif

    elseif (y.ge.0) then 
       if ((oldy.lt.cmir_refl_upper).and.(y.ge.cmir_refl_upper)) then 
          !
          !                   Reflection at upper (>0) mirror
          !

          y = cmir_refl_upper - abs(abs(y)-cmir_refl_upper)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

       elseif ((oldy.gt.(lim_sep+cmir_refl_lower)).and.(y.le.(lim_sep+cmir_refl_lower))) then 
          !
          !                   Reflection at second lower (>0) mirror
          !

          y = lim_sep + cmir_refl_lower + abs(abs(y)+cmir_refl_lower)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0
          reflected = .true.

       endif
    endif

    if (debugl.and.reflected) then
       write(error_message_data,'(a,i10,6(1x,g18.10))') 'REFLECTION:',int(yreflection_event_count),y_org,oldy_org,y,oldy,svy_org,svy
       call dbgmsg('CHECK_REFLECTION',error_message_data)
    endif

    if (check_reflected_region(y)) then 
       
       !write(error_message_data,'(a,3(1x,g18.10))') 'REFLECTED PARTICLE HAS ENTERED MIRROR REGION - '//&
       !                                           & 'TRY REDUCING SIMULATION TIMESTEPS AND MULTIPLIERS : DATA:', y,oldy,svy
       !call errmsg('CHECK REFLECTION:WARNING:',error_message_data)
       
       y_tmp = y

       !ran = getranf()
       !y = cmir_refl_lower + ran * (cmir_refl_upper-cmir_refl_lower)

       relocation_count = relocation_count +1.0
       write(error_message_data,'(a,i12,4(1x,g18.10))') 'REVISED Y COORDINATES:',int(relocation_count), y_org,oldy,y_tmp,y
       call errmsg('CHECK_REFLECTION PARTICLE RELOCATED:',error_message_data)


       ierr =1
    endif


  end subroutine check_reflection




end module yreflection
