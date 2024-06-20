subroutine init_plasma()
  use all_variables, only : zones, global_parameters, flags
  implicit none
  if(flags%restart) then
     call load_plasma()
  else
     call set_default_initial_plasma()
  end if
end subroutine init_plasma
