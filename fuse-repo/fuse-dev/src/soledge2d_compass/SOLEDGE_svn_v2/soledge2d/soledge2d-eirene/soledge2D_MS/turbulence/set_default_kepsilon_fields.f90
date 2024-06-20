subroutine set_default_kepsilon_fields()
  use all_variables , only : global_parameters, zones
  use Mdefinitions
  implicit none
  integer*4 :: Nx, Nz, k
  do k=1,global_parameters%N_zones
     zones(k)%kepsilon(STEP_OLD)%k=1.
     zones(k)%kepsilon(STEP_OLD)%epsilon=1.
  end do
end subroutine set_default_kepsilon_fields
