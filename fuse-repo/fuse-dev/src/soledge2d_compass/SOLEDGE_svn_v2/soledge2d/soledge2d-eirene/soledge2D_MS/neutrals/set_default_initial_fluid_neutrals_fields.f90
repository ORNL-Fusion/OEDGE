subroutine set_default_initial_fluid_neutrals_fields()
  use all_variables , only : global_parameters, zones
  use Mdefinitions
  implicit none
  integer*4 :: Nx, Nz, k
  do k=1,global_parameters%N_zones
     zones(k)%neutrals%density=1.d-1
  end do
end subroutine set_default_initial_fluid_neutrals_fields
