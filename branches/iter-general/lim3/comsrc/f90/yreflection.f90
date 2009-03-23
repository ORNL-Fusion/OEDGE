module yreflection
  use error_handling
  implicit none

  integer :: yreflection_opt
  real :: cmir_refl_lower, cmir_refl_upper
  real :: yreflection_event_count = 0
  save

contains

  subroutine test_reflection(y,oldy,svy)
    implicit none
    real :: y,oldy,svy

    if (yreflection_opt.ne.0.and.(cmir_refl_lower.eq.0.0.or.cmir_refl_upper.eq.0.0)) then 
       call errmsg('YREFLECTION:TEST_REFLECTION','Y-REFLECTION OPTION IS ON BUT REFLECTION&
                   & LOCATIONS HAVE NOT BEEN PROPERLY SPECIFIED : Y-REFLECTION DISABLED')
       yreflection_opt=0
       return
    endif


    if (y.lt.0) then 
       if ((oldy.gt.cmir_refl_lower).and.(y.le.cmir_refl_lower)) then 
          !
          !                   Reflection at lower (<0) mirror
          !

          y = cmir_refl_lower + abs(y-cmir_refl_lower)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0

       elseif ((oldy.lt.-cmir_refl_upper).and.(y.ge.-cmir_refl_upper)) then 
          !
          !                   Reflection at second upper (< 0) mirror
          !
          
          y = -cmir_refl_upper - abs(y-cmir_refl_upper)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0

       endif

    elseif (y.ge.0) then 
       if ((oldy.lt.cmir_refl_upper).and.(y.ge.cmir_refl_upper)) then 
          !
          !                   Reflection at upper (>0) mirror
          !

          y = cmir_refl_upper - abs(y-cmir_refl_upper)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0


       elseif ((oldy.gt.-cmir_refl_lower).and.(y.le.-cmir_refl_lower)) then 
          !
          !                   Reflection at second lower (>0) mirror
          !

          y = -cmir_refl_lower + abs(y-cmir_refl_lower)
          svy = -svy
          yreflection_event_count = yreflection_event_count + 1.0


       endif
    endif


  end subroutine test_reflection




end module yreflection
