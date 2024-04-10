subroutine smooth_phi()
  use all_variables, only : global_parameters, zones
  use Msmoothing_vars
  implicit none
  integer*4 :: k

  do k=1,global_parameters%N_Zones
     call init_RHS_smoothing(zones(k))
  end do

  call compute_smoothing_rhs()
  call solve_implicit_pastix(CSC_smoothing)
  call redistribute_smoothing_sol(CSC_smoothing)

end subroutine smooth_phi
