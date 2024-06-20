subroutine set_continuity_boundary_conditions_south_ke(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : boundary_conditions, global_parameters&
       ,reference_parameters, element_variables, ballooning_parameters
  use Meirene_vars
  use Mfeedback_control
  use MZone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: Nx,Nz
  real*8 :: rs0,c0,n0,T0eV
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  n0=reference_parameters%fields%n0
  T0eV=reference_parameters%fields%T0eV
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  ! set density BC to the completely ionized ion for the element
  ! for the other ions, copy density (grad approx 0)
  zone%kepsilon(STEP)%k(0,1:Nz) = zone%kepsilon(STEP)%k(1,1:Nz)
  zone%kepsilon(STEP)%epsilon(0,1:Nz) = zone%kepsilon(STEP)%epsilon(1,1:Nz)
end subroutine set_continuity_boundary_conditions_south_ke
