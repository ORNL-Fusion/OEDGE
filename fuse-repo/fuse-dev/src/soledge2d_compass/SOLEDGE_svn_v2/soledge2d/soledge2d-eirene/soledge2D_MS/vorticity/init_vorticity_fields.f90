subroutine init_vorticity_fields()
  use all_variables, only : zones, global_parameters, flags, drift_flags
  implicit none
  if(drift_flags%vorticity_restart) then
     call load_vorticity_fields()
  else
     call set_default_initial_vorticity_fields()
  end if
end subroutine init_vorticity_fields
