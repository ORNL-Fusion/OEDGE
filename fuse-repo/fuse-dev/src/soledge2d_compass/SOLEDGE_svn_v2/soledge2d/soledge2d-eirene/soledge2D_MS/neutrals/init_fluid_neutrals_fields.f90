subroutine init_fluid_neutrals_fields()
  use all_variables, only : zones, global_parameters, flags
  implicit none
  call read_fluid_neutrals_input_file()
  call compute_FN_diffusivities()
  if(flags%restart) then
     call load_fluid_neutrals_fields()
  else
     call set_default_initial_fluid_neutrals_fields()
  end if
end subroutine init_fluid_neutrals_fields
