subroutine compute_coupling_terms(zone)
  use all_variables, only : global_parameters
  use MZone
  implicit none
  type(TZone),intent(inout) :: zone
  call compute_collision_times(zone)
  call compute_friction_force(zone)
  call compute_thermal_force(zone)
  call compute_energy_exchange(zone)
  call compute_collisional_velocity_heat_flux(zone)
end subroutine compute_coupling_terms
