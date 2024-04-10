subroutine set_continuity_boundary_conditions_south(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : global_parameters
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4, intent(in) :: STEP
  integer*4 :: n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     zone%species(n)%var(STEP)%density(0,1:Nz) = zone%species(n)%var(STEP)%density(1,1:Nz) 
     zone%species(n)%var(STEP)%Gamma(0,1:Nz) = zone%species(n)%var(STEP)%Gamma(1,1:Nz) 
     zone%species(n)%var(STEP)%temperature(0,1:Nz) = zone%species(n)%var(STEP)%temperature(1,1:Nz)
  end do
end subroutine set_continuity_boundary_conditions_south
