module iter_bm
use error_handling
implicit none

  ! ITER BM shaping inputs
  real    :: rtor_setback,rslot_setback,bm_tor_wid,slot_tor_wid,lambda_design
  ! ITER BM parameters
  real :: c_lim, y_re

  save

contains


  subroutine calc_iter_limiter_parameters
    implicit none
    !
    ! See "ITER Limiter Notes. Ideal Shape of the face of the outer wall BMs for startup/rampdown. 
    !      Set of 18 BMs in a toroidal ring at the outer wall" by P.C. Stangeby for details on the 
    !      formulae and calculations used here. Also see notes and comments in edge.f
    !


    y_re = ( bm_tor_wid * (1-exp(-rslot_setback/lambda_design)) +slot_tor_wid*(1-exp(-rtor_setback/lambda_design)))/ &
         ( (1-exp(-rslot_setback/lambda_design)) +(1-exp(-rtor_setback/lambda_design)))

    c_lim = (lambda_design * (1-exp(-rtor_setback/lambda_design))) / (bm_tor_wid - y_re)

  end subroutine calc_iter_limiter_parameters


end module iter_bm

