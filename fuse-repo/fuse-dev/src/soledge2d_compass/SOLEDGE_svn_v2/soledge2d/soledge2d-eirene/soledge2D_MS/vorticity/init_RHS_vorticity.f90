subroutine init_RHS_vorticity(zone)
  use all_variables, only : global_variables
  use MZone
  use MOperator
  use MDiffusion_perp
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 ::  Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Zone%electric_fields(1)%RHS=zone%electric_fields(1)%vorticity(1:Nx,1:Nz)
end subroutine init_RHS_vorticity
