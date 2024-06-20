subroutine MD_broadcast_geometry()
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 :: k,n
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  integer*4 :: North,South,East,West
  integer*4 :: mNorth,mSouth,mEast,mWest
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     !North boundary
     North=zones(k)%Neighbors(1)
     mNorth=zones(k)%MagNeighbors(1)
     if(mNorth.eq.0) then
        zones(k)%mesh%Rgeom(Nx+1,1:Nz)=zones(North)%mesh%Rgeom(1,1:Nz)
        zones(k)%mesh%Zgeom(Nx+1,1:Nz)=zones(North)%mesh%Zgeom(1,1:Nz)
        zones(k)%mesh%Br(Nx+1,1:Nz)=zones(North)%mesh%Br(1,1:Nz)
        zones(k)%mesh%Bz(Nx+1,1:Nz)=zones(North)%mesh%Bz(1,1:Nz)
        zones(k)%mesh%Bphi(Nx+1,1:Nz)=zones(North)%mesh%Bphi(1,1:Nz)
        zones(k)%mesh%B(Nx+1,1:Nz)=zones(North)%mesh%B(1,1:Nz)
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     mSouth=zones(k)%MagNeighbors(2)
     if(mSouth.eq.0) then
        Nx_N=zones(South)%mesh%Nx
        zones(k)%mesh%Rgeom(0,1:Nz)=zones(South)%mesh%Rgeom(Nx_N,1:Nz)
        zones(k)%mesh%Zgeom(0,1:Nz)=zones(South)%mesh%Zgeom(Nx_N,1:Nz)
        zones(k)%mesh%Br(0,1:Nz)=zones(South)%mesh%Br(Nx_N,1:Nz)
        zones(k)%mesh%Bz(0,1:Nz)=zones(South)%mesh%Bz(Nx_N,1:Nz)
        zones(k)%mesh%Bphi(0,1:Nz)=zones(South)%mesh%Bphi(Nx_N,1:Nz)
        zones(k)%mesh%B(0,1:Nz)=zones(South)%mesh%B(Nx_N,1:Nz)
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     mEast=zones(k)%MagNeighbors(3)
     if(mEast.eq.0) then
        zones(k)%mesh%Rgeom(1:Nx,Nz+1)=zones(East)%mesh%Rgeom(1:Nx,1)
        zones(k)%mesh%Zgeom(1:Nx,Nz+1)=zones(East)%mesh%Zgeom(1:Nx,1)
        zones(k)%mesh%Br(1:Nx,Nz+1)=zones(East)%mesh%Br(1:Nx,1)
        zones(k)%mesh%Bz(1:Nx,Nz+1)=zones(East)%mesh%Bz(1:Nx,1)
        zones(k)%mesh%Bphi(1:Nx,Nz+1)=zones(East)%mesh%Bphi(1:Nx,1)
        zones(k)%mesh%B(1:Nx,Nz+1)=zones(East)%mesh%B(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     mWest=zones(k)%MagNeighbors(4)
     if(mWest.eq.0) then
        Nz_N=zones(West)%mesh%Nz
        zones(k)%mesh%Rgeom(1:Nx,0)=zones(West)%mesh%Rgeom(1:Nx,Nz_N)
        zones(k)%mesh%Zgeom(1:Nx,0)=zones(West)%mesh%Zgeom(1:Nx,Nz_N)
        zones(k)%mesh%Br(1:Nx,0)=zones(West)%mesh%Br(1:Nx,Nz_N)
        zones(k)%mesh%Bz(1:Nx,0)=zones(West)%mesh%Bz(1:Nx,Nz_N)
        zones(k)%mesh%Bphi(1:Nx,0)=zones(West)%mesh%Bphi(1:Nx,Nz_N)
        zones(k)%mesh%B(1:Nx,0)=zones(West)%mesh%B(1:Nx,Nz_N)
     end if
  end do
end subroutine MD_broadcast_geometry
