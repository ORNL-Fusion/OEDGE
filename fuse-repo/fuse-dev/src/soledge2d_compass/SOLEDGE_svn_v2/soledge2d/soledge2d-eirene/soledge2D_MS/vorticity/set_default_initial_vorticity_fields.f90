subroutine set_default_initial_vorticity_fields()
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 :: k
  do k=1,global_parameters%N_Zones
     zones(k)%electric_fields(1)%phi=3.d0*zones(k)%species(0)%var(1)%temperature
     zones(k)%electric_fields(1)%vorticity=0.D0
  end do
end subroutine set_default_initial_vorticity_fields
