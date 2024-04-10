subroutine gather_neighboring_values(zone,SIDE,STEP)
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
     do n=0,global_parameters%N_ions
        zone%species(n)%var(STEP)%density(Nx+1,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%density(1,1:Nz)
        zone%species(n)%var(STEP)%Gamma(Nx+1,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%Gamma(1,1:Nz)
        zone%species(n)%var(STEP)%temperature(Nx+1,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%temperature(1,1:Nz)
     end do
#if VORTICITY_PASTIX == 1
     zone%electric_fields(STEP)%phi(Nx+1,1:Nz)=zones(Neighbor)%electric_fields(STEP)%phi(1,1:Nz)
#endif
  case(N_SOUTH)
     Neighbor=zone%Neighbors(2)
     Nx_N=zones(Neighbor)%mesh%Nx
     do n=0,global_parameters%N_ions
        zone%species(n)%var(STEP)%density(0,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%density(Nx_N,1:Nz)
        zone%species(n)%var(STEP)%Gamma(0,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%Gamma(Nx_N,1:Nz)
        zone%species(n)%var(STEP)%temperature(0,1:Nz)=zones(Neighbor)%species(n)%var(STEP)%temperature(Nx_N,1:Nz)
     end do
#if VORTICITY_PASTIX == 1
     zone%electric_fields(STEP)%phi(0,1:Nz)=zones(Neighbor)%electric_fields(STEP)%phi(Nx_N,1:Nz)
#endif
  case(N_EAST)
     Neighbor=zone%Neighbors(3)
     do n=0,global_parameters%N_ions
        zone%species(n)%var(STEP)%density(1:Nx,Nz+1)=zones(Neighbor)%species(n)%var(STEP)%density(1:Nx,1)
        zone%species(n)%var(STEP)%Gamma(1:Nx,Nz+1)=zones(Neighbor)%species(n)%var(STEP)%Gamma(1:Nx,1)
        zone%species(n)%var(STEP)%temperature(1:Nx,Nz+1)=zones(Neighbor)%species(n)%var(STEP)%temperature(1:Nx,1)
     end do
#if VORTICITY_PASTIX == 1
     zone%electric_fields(STEP)%phi(1:Nx,Nz+1)=zones(Neighbor)%electric_fields(STEP)%phi(1:Nx,1)
#endif
  case(N_WEST)
     Neighbor=zone%Neighbors(4)
     Nz_N=zones(Neighbor)%mesh%Nz
     do n=0,global_parameters%N_ions
        zone%species(n)%var(STEP)%density(1:Nx,0)=zones(Neighbor)%species(n)%var(STEP)%density(1:Nx,Nz_N)
        zone%species(n)%var(STEP)%Gamma(1:Nx,0)=zones(Neighbor)%species(n)%var(STEP)%Gamma(1:Nx,Nz_N)
        zone%species(n)%var(STEP)%temperature(1:Nx,0)=zones(Neighbor)%species(n)%var(STEP)%temperature(1:Nx,Nz_N)
     end do
#if VORTICITY_PASTIX == 1
     zone%electric_fields(STEP)%phi(1:Nx,0)=zones(Neighbor)%electric_fields(STEP)%phi(1:Nx,Nz_N)
#endif
  end select
end subroutine gather_neighboring_values


