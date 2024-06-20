subroutine set_wall_boundary_conditions_north(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : global_parameters, global_variables
  use MZone
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     zone%species(n)%var(STEP)%density(Nx+1,1:Nz) = zone%species(n)%var(STEP)%density(Nx,1:Nz) !global_variables%min_density
     zone%species(n)%var(STEP)%Gamma(Nx+1,1:Nz) = 0.
     zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz) = zone%species(n)%var(STEP)%temperature(Nx,1:Nz)!  global_variables%min_temperature
  end do
end subroutine set_wall_boundary_conditions_north
