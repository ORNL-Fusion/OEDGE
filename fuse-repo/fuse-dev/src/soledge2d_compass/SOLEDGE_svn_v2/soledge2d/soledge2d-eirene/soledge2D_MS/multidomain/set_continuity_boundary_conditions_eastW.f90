subroutine set_continuity_boundary_conditions_eastW(zone,STEP)
  use all_variables, only : global_parameters
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  zone%electric_fields(STEP)%vorticity(1:Nx,Nz+1) = zone%electric_fields(STEP)%vorticity(1:Nx,Nz)
end subroutine set_continuity_boundary_conditions_eastW
