subroutine init_kepsilon_fields()
  use all_variables, only : global_parameters, flags
  implicit none
  if(flags%restart) then
     call load_kepsilon_fields()
  else
     call set_default_kepsilon_fields()
  end if
end subroutine init_kepsilon_fields
