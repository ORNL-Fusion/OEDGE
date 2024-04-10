subroutine gather_neighboring_values_ke(zone,SIDE,STEP)
#include "compile_opt.inc"
  use all_variables, only : zones, global_parameters
  use Mzone
  implicit none
  integer*4 :: n
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: SIDE
  integer*4,intent(in) :: STEP
  integer*4 :: Neighbor
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  select case(SIDE)
  case(N_NORTH)
     Neighbor=zone%Neighbors(1)
     zone%kepsilon(STEP)%k(Nx+1,1:Nz)=zones(Neighbor)%kepsilon(STEP)%k(1,1:Nz)
     zone%kepsilon(STEP)%epsilon(Nx+1,1:Nz)=zones(Neighbor)%kepsilon(STEP)%epsilon(1,1:Nz)
  case(N_SOUTH)
     Neighbor=zone%Neighbors(2)
     Nx_N=zones(Neighbor)%mesh%Nx
     zone%kepsilon(STEP)%k(0,1:Nz)=zones(Neighbor)%kepsilon(STEP)%k(Nx_N,1:Nz)
     zone%kepsilon(STEP)%epsilon(0,1:Nz)=zones(Neighbor)%kepsilon(STEP)%epsilon(Nx_N,1:Nz)
  case(N_EAST)
     Neighbor=zone%Neighbors(3)
     zone%kepsilon(STEP)%k(1:Nx,Nz+1)=zones(Neighbor)%kepsilon(STEP)%k(1:Nx,1)
     zone%kepsilon(STEP)%epsilon(1:Nx,Nz+1)=zones(Neighbor)%kepsilon(STEP)%epsilon(1:Nx,1)
  case(N_WEST)
     Neighbor=zone%Neighbors(4)
     Nz_N=zones(Neighbor)%mesh%Nz
     zone%kepsilon(STEP)%k(1:Nx,0)=zones(Neighbor)%kepsilon(STEP)%k(1:Nx,Nz_N)
     zone%kepsilon(STEP)%epsilon(1:Nx,0)=zones(Neighbor)%kepsilon(STEP)%epsilon(1:Nx,Nz_N)
  end select
end subroutine gather_neighboring_values_ke


