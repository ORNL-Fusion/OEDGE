subroutine solve_neutrals()
  use all_variables, only : global_parameters, zones
  use Mneutral_vars
  implicit none
  integer*4 :: k

  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k)
  !$OMP DO
  do k=1,global_parameters%N_Zones
     call init_RHS_neutrals(zones(k))
     call compute_recycling_neutral_sources(zones(k))
     call add_RHS_neutrals_sources(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER
  !$OMP END PARALLEL

  call compute_neutral_RHS()
  call solve_implicit_pastix(CSC_neutral)
  call redistribute_neutral_sol(CSC_neutral)

end subroutine solve_neutrals
