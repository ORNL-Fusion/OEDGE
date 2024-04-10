subroutine compute_test_advection_sources()
  use test_var
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: k
  do k=1,global_parameters%N_zones
     call compute_test_advection_sourceN(k)
     call compute_test_advection_sourceG(k)
     call compute_test_advection_sourceT(k)
  end do
end subroutine compute_test_advection_sources
