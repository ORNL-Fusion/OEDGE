subroutine MD_broadcast_corners_mesh
  use all_variables, only : zones,global_parameters
  use MZone
  implicit none
  integer*4 :: k
  integer*4 :: North,South,East,West
  integer*4 :: Nx_N,Nz_N
  do k=1,global_parameters%N_zones
     North=zones(k)%Neighbors(1)
     if(North.gt.0) then
        Nz_N=zones(North)%mesh%Nz
        zones(k)%mesh%cornerx(1,1,1)=zones(North)%mesh%x(1,0)
        zones(k)%mesh%cornerx(1,2,1)=zones(North)%mesh%x(1,Nz_N+1)
        zones(k)%mesh%cornerz(1,1,1)=zones(North)%mesh%z(1,0)
        zones(k)%mesh%cornerz(1,2,1)=zones(North)%mesh%z(1,Nz_N+1)
        zones(k)%mesh%cornerRg(1,1,1)=zones(North)%mesh%Rgeom(1,0)
        zones(k)%mesh%cornerRg(1,2,1)=zones(North)%mesh%Rgeom(1,Nz_N+1)
        zones(k)%mesh%cornerZg(1,1,1)=zones(North)%mesh%Zgeom(1,0)
        zones(k)%mesh%cornerZg(1,2,1)=zones(North)%mesh%Zgeom(1,Nz_N+1)
     end if
     South=zones(k)%Neighbors(2)
     if(South.gt.0) then
        Nz_N=zones(South)%mesh%Nz
        Nx_N=zones(South)%mesh%Nx
        zones(k)%mesh%cornerx(2,1,1)=zones(South)%mesh%x(Nx_N,0)
        zones(k)%mesh%cornerx(2,2,1)=zones(South)%mesh%x(Nx_N,Nz_N+1)
        zones(k)%mesh%cornerz(2,1,1)=zones(South)%mesh%z(Nx_N,0)
        zones(k)%mesh%cornerz(2,2,1)=zones(South)%mesh%z(Nx_N,Nz_N+1)
        zones(k)%mesh%cornerRg(2,1,1)=zones(South)%mesh%Rgeom(Nx_N,0)
        zones(k)%mesh%cornerRg(2,2,1)=zones(South)%mesh%Rgeom(Nx_N,Nz_N+1)
        zones(k)%mesh%cornerZg(2,1,1)=zones(South)%mesh%Zgeom(Nx_N,0)
        zones(k)%mesh%cornerZg(2,2,1)=zones(South)%mesh%Zgeom(Nx_N,Nz_N+1)
     end if
     East=zones(k)%Neighbors(3)
     if(East.gt.0) then
        Nx_N=zones(East)%mesh%Nx
        Nz_N=zones(East)%mesh%Nz
        zones(k)%mesh%cornerx(1,2,2)=zones(East)%mesh%x(Nx_N+1,1)
        zones(k)%mesh%cornerx(2,2,2)=zones(East)%mesh%x(0,1)
        zones(k)%mesh%cornerz(1,2,2)=zones(East)%mesh%z(Nx_N+1,1)
        zones(k)%mesh%cornerz(2,2,2)=zones(East)%mesh%z(0,1)
        zones(k)%mesh%cornerRg(1,2,2)=zones(East)%mesh%Rgeom(Nx_N+1,1)
        zones(k)%mesh%cornerRg(2,2,2)=zones(East)%mesh%Rgeom(0,1)
        zones(k)%mesh%cornerZg(1,2,2)=zones(East)%mesh%Zgeom(Nx_N+1,1)
        zones(k)%mesh%cornerZg(2,2,2)=zones(East)%mesh%Zgeom(0,1)
     end if
     West=zones(k)%Neighbors(4)
     if(West.gt.0) then
        Nx_N=zones(West)%mesh%Nx
        Nz_N=zones(West)%mesh%Nz
        zones(k)%mesh%cornerx(1,1,2)=zones(West)%mesh%x(Nx_N+1,Nz_N)
        zones(k)%mesh%cornerx(2,1,2)=zones(West)%mesh%x(0,Nz_N)
        zones(k)%mesh%cornerz(1,1,2)=zones(West)%mesh%z(Nx_N+1,Nz_N)
        zones(k)%mesh%cornerz(2,1,2)=zones(West)%mesh%z(0,Nz_N)
        zones(k)%mesh%cornerRg(1,1,2)=zones(West)%mesh%Rgeom(Nx_N+1,Nz_N)
        zones(k)%mesh%cornerRg(2,1,2)=zones(West)%mesh%Rgeom(0,Nz_N)
        zones(k)%mesh%cornerZg(1,1,2)=zones(West)%mesh%Zgeom(Nx_N+1,Nz_N)
        zones(k)%mesh%cornerZg(2,1,2)=zones(West)%mesh%Zgeom(0,Nz_N)
     end if     
 end do
end subroutine MD_broadcast_corners_mesh
