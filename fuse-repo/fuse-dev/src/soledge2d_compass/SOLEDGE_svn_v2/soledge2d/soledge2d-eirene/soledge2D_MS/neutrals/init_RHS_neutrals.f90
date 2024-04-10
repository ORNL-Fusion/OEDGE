subroutine init_RHS_neutrals(zone)
  use all_variables, only : global_variables
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 ::  Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Zone%neutrals%RHS=zone%neutrals%density(1:Nx,1:Nz)/global_variables%dt
end subroutine init_RHS_neutrals
