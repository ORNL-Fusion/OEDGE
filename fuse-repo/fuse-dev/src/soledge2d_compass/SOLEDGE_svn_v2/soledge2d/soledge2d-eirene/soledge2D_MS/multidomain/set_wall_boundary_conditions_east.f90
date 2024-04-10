subroutine set_wall_boundary_conditions_east(zone,STEP)
  use all_variables, only : global_parameters
  use Mzone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
end subroutine set_wall_boundary_conditions_east
