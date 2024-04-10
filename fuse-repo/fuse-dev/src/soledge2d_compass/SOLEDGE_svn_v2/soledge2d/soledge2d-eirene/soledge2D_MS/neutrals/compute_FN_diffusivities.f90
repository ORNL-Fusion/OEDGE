subroutine compute_FN_diffusivities()
  use all_variables, only : reference_parameters, zones, global_parameters
  use Mneutral_vars
  implicit none
  integer*4 :: k
  do k=1,global_parameters%N_zones
     zones(k)%neutrals%Dn=FN_diffusivity&
          /(reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
  end do
end subroutine compute_FN_diffusivities
