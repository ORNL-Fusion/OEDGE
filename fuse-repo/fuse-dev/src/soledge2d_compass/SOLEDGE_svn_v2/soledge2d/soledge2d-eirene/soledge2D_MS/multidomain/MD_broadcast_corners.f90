subroutine MD_broadcast_corners(zone,STEP)
#include "compile_opt.inc"
  use all_variables, only : zones,global_parameters
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4 :: n
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  do n=0,global_parameters%N_ions
     North=zone%Neighbors(1)
     if(North.gt.0) then
        Nz_N=zones(North)%mesh%Nz
        zone%species(n)%corners%density(1,1,1)=zones(North)%species(n)%var(1)%density(1,0)
        zone%species(n)%corners%density(1,2,1)=zones(North)%species(n)%var(1)%density(1,Nz_N+1)
        zone%species(n)%corners%Gamma(1,1,1)=zones(North)%species(n)%var(1)%Gamma(1,0)
        zone%species(n)%corners%Gamma(1,2,1)=zones(North)%species(n)%var(1)%Gamma(1,Nz_N+1)
        zone%species(n)%corners%velocity(1,1,1)=zones(North)%species(n)%var(1)%velocity(1,0)
        zone%species(n)%corners%velocity(1,2,1)=zones(North)%species(n)%var(1)%velocity(1,Nz_N+1)
        zone%species(n)%corners%temperature(1,1,1)=zones(North)%species(n)%var(1)%temperature(1,0)
        zone%species(n)%corners%temperature(1,2,1)=zones(North)%species(n)%var(1)%temperature(1,Nz_N+1)
     end if
     South=zone%Neighbors(2)
     if(South.gt.0) then
        Nz_N=zones(South)%mesh%Nz
        Nx_N=zones(South)%mesh%Nx
        zone%species(n)%corners%density(2,1,1)=zones(South)%species(n)%var(1)%density(Nx_N,0)
        zone%species(n)%corners%density(2,2,1)=zones(South)%species(n)%var(1)%density(Nx_N,Nz_N+1)
        zone%species(n)%corners%Gamma(2,1,1)=zones(South)%species(n)%var(1)%Gamma(Nx_N,0)
        zone%species(n)%corners%Gamma(2,2,1)=zones(South)%species(n)%var(1)%Gamma(Nx_N,Nz_N+1)
        zone%species(n)%corners%velocity(2,1,1)=zones(South)%species(n)%var(1)%velocity(Nx_N,0)
        zone%species(n)%corners%velocity(2,2,1)=zones(South)%species(n)%var(1)%velocity(Nx_N,Nz_N+1)
        zone%species(n)%corners%temperature(2,1,1)=zones(South)%species(n)%var(1)%temperature(Nx_N,0)
        zone%species(n)%corners%temperature(2,2,1)=zones(South)%species(n)%var(1)%temperature(Nx_N,Nz_N+1)
     end if
     East=zone%Neighbors(3)
     if(East.gt.0) then
        Nx_N=zones(East)%mesh%Nx
        Nz_N=zones(East)%mesh%Nz
        zone%species(n)%corners%density(1,2,2)=zones(East)%species(n)%var(1)%density(Nx_N+1,1)
        zone%species(n)%corners%density(2,2,2)=zones(East)%species(n)%var(1)%density(0,1)
        zone%species(n)%corners%Gamma(1,2,2)=zones(East)%species(n)%var(1)%Gamma(Nx_N+1,1)
        zone%species(n)%corners%Gamma(2,2,2)=zones(East)%species(n)%var(1)%Gamma(0,1)
        zone%species(n)%corners%velocity(1,2,2)=zones(East)%species(n)%var(1)%velocity(Nx_N+1,1)
        zone%species(n)%corners%velocity(2,2,2)=zones(East)%species(n)%var(1)%velocity(0,1)
        zone%species(n)%corners%temperature(1,2,2)=zones(East)%species(n)%var(1)%temperature(Nx_N+1,1)
        zone%species(n)%corners%temperature(2,2,2)=zones(East)%species(n)%var(1)%temperature(0,1)
     end if
     West=zone%Neighbors(4)
     if(West.gt.0) then
        Nx_N=zones(West)%mesh%Nx
        Nz_N=zones(West)%mesh%Nz
        zone%species(n)%corners%density(1,1,2)=zones(West)%species(n)%var(1)%density(Nx_N+1,Nz_N)
        zone%species(n)%corners%density(2,1,2)=zones(West)%species(n)%var(1)%density(0,Nz_N)
        zone%species(n)%corners%Gamma(1,1,2)=zones(West)%species(n)%var(1)%Gamma(Nx_N+1,Nz_N)
        zone%species(n)%corners%Gamma(2,1,2)=zones(West)%species(n)%var(1)%Gamma(0,Nz_N)
        zone%species(n)%corners%velocity(1,1,2)=zones(West)%species(n)%var(1)%velocity(Nx_N+1,Nz_N)
        zone%species(n)%corners%velocity(2,1,2)=zones(West)%species(n)%var(1)%velocity(0,Nz_N)
        zone%species(n)%corners%temperature(1,1,2)=zones(West)%species(n)%var(1)%temperature(Nx_N+1,Nz_N)
        zone%species(n)%corners%temperature(2,1,2)=zones(West)%species(n)%var(1)%temperature(0,Nz_N)
     end if
     if(North.lt.0) then
        zone%species(n)%corners%density(1,1,1)=zone%species(n)%corners%density(1,1,2)
        zone%species(n)%corners%density(1,2,1)=zone%species(n)%corners%density(1,2,2)
        zone%species(n)%corners%Gamma(1,1,1)=zone%species(n)%corners%Gamma(1,1,2)
        zone%species(n)%corners%Gamma(1,2,1)=zone%species(n)%corners%Gamma(1,2,2)
        zone%species(n)%corners%velocity(1,1,1)=zone%species(n)%corners%velocity(1,1,2)
        zone%species(n)%corners%velocity(1,2,1)=zone%species(n)%corners%velocity(1,2,2)
        zone%species(n)%corners%temperature(1,1,1)=zone%species(n)%corners%temperature(1,1,2)
        zone%species(n)%corners%temperature(1,2,1)=zone%species(n)%corners%temperature(1,2,2)
     end if
     if(South.lt.0) then
        zone%species(n)%corners%density(2,1,1)=zone%species(n)%corners%density(2,1,2)
        zone%species(n)%corners%density(2,2,1)=zone%species(n)%corners%density(2,2,2)
        zone%species(n)%corners%Gamma(2,1,1)=zone%species(n)%corners%Gamma(2,1,2)
        zone%species(n)%corners%Gamma(2,2,1)=zone%species(n)%corners%Gamma(2,2,2)
        zone%species(n)%corners%velocity(2,1,1)=zone%species(n)%corners%velocity(2,1,2)
        zone%species(n)%corners%velocity(2,2,1)=zone%species(n)%corners%velocity(2,2,2)
        zone%species(n)%corners%temperature(2,1,1)=zone%species(n)%corners%temperature(2,1,2)
        zone%species(n)%corners%temperature(2,2,1)=zone%species(n)%corners%temperature(2,2,2)
     end if
     if(East.lt.0) then
        zone%species(n)%corners%density(1,2,2)=zone%species(n)%corners%density(1,2,1)
        zone%species(n)%corners%density(2,2,2)=zone%species(n)%corners%density(2,2,1)
        zone%species(n)%corners%Gamma(1,2,2)=zone%species(n)%corners%Gamma(1,2,1)
        zone%species(n)%corners%Gamma(2,2,2)=zone%species(n)%corners%Gamma(2,2,1)
        zone%species(n)%corners%velocity(1,2,2)=zone%species(n)%corners%velocity(1,2,1)
        zone%species(n)%corners%velocity(2,2,2)=zone%species(n)%corners%velocity(2,2,1)
        zone%species(n)%corners%temperature(1,2,2)=zone%species(n)%corners%temperature(1,2,1)
        zone%species(n)%corners%temperature(2,2,2)=zone%species(n)%corners%temperature(2,2,1)
     end if
     if(West.lt.0) then
        zone%species(n)%corners%density(1,1,2)=zone%species(n)%corners%density(1,1,1)
        zone%species(n)%corners%density(2,1,2)=zone%species(n)%corners%density(2,1,1)
        zone%species(n)%corners%Gamma(1,1,2)=zone%species(n)%corners%Gamma(1,1,1)
        zone%species(n)%corners%Gamma(2,1,2)=zone%species(n)%corners%Gamma(2,1,1)
        zone%species(n)%corners%velocity(1,1,2)=zone%species(n)%corners%velocity(1,1,1)
        zone%species(n)%corners%velocity(2,1,2)=zone%species(n)%corners%velocity(2,1,1)
        zone%species(n)%corners%temperature(1,1,2)=zone%species(n)%corners%temperature(1,1,1)
        zone%species(n)%corners%temperature(2,1,2)=zone%species(n)%corners%temperature(2,1,1)
     end if
  end do
end subroutine MD_broadcast_corners
