subroutine MD_broadcast_jacobian()
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
        zones(k)%jacobian%dTdR(Nx+1,1:Nz)=zones(North)%jacobian%dTdR(1,1:Nz)
        zones(k)%jacobian%dTdZ(Nx+1,1:Nz)=zones(North)%jacobian%dTdZ(1,1:Nz)
        zones(k)%jacobian%dPdR(Nx+1,1:Nz)=zones(North)%jacobian%dPdR(1,1:Nz)
        zones(k)%jacobian%dPdZ(Nx+1,1:Nz)=zones(North)%jacobian%dPdZ(1,1:Nz)
        zones(k)%jacobian%dRdT(Nx+1,1:Nz)=zones(North)%jacobian%dRdT(1,1:Nz)
        zones(k)%jacobian%dRdP(Nx+1,1:Nz)=zones(North)%jacobian%dRdP(1,1:Nz)
        zones(k)%jacobian%dZdT(Nx+1,1:Nz)=zones(North)%jacobian%dZdT(1,1:Nz)
        zones(k)%jacobian%dZdP(Nx+1,1:Nz)=zones(North)%jacobian%dZdP(1,1:Nz)
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     mSouth=zones(k)%MagNeighbors(2)
     if(mSouth.eq.0) then
        Nx_N=zones(South)%mesh%Nx
        zones(k)%jacobian%dTdR(0,1:Nz)=zones(South)%jacobian%dTdR(Nx_N,1:Nz)
        zones(k)%jacobian%dTdZ(0,1:Nz)=zones(South)%jacobian%dTdZ(Nx_N,1:Nz)
        zones(k)%jacobian%dPdR(0,1:Nz)=zones(South)%jacobian%dPdR(Nx_N,1:Nz)
        zones(k)%jacobian%dPdZ(0,1:Nz)=zones(South)%jacobian%dPdZ(Nx_N,1:Nz)
        zones(k)%jacobian%dRdT(0,1:Nz)=zones(South)%jacobian%dRdT(Nx_N,1:Nz)
        zones(k)%jacobian%dRdP(0,1:Nz)=zones(South)%jacobian%dRdP(Nx_N,1:Nz)
        zones(k)%jacobian%dZdT(0,1:Nz)=zones(South)%jacobian%dZdT(Nx_N,1:Nz)
        zones(k)%jacobian%dZdP(0,1:Nz)=zones(South)%jacobian%dZdP(Nx_N,1:Nz)
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     mEast=zones(k)%MagNeighbors(3)
     if(mEast.eq.0) then
        zones(k)%jacobian%dTdR(1:Nx,Nz+1)=zones(East)%jacobian%dTdR(1:Nx,1)
        zones(k)%jacobian%dTdZ(1:Nx,Nz+1)=zones(East)%jacobian%dTdZ(1:Nx,1)
        zones(k)%jacobian%dPdR(1:Nx,Nz+1)=zones(East)%jacobian%dPdR(1:Nx,1)
        zones(k)%jacobian%dPdZ(1:Nx,Nz+1)=zones(East)%jacobian%dPdZ(1:Nx,1)
        zones(k)%jacobian%dRdT(1:Nx,Nz+1)=zones(East)%jacobian%dRdT(1:Nx,1)
        zones(k)%jacobian%dRdP(1:Nx,Nz+1)=zones(East)%jacobian%dRdP(1:Nx,1)
        zones(k)%jacobian%dZdT(1:Nx,Nz+1)=zones(East)%jacobian%dZdT(1:Nx,1)
        zones(k)%jacobian%dZdP(1:Nx,Nz+1)=zones(East)%jacobian%dZdP(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     mWest=zones(k)%MagNeighbors(4)
     if(mWest.eq.0) then
        Nz_N=zones(West)%mesh%Nz
        zones(k)%jacobian%dTdR(1:Nx,0)=zones(West)%jacobian%dTdR(1:Nx,Nz_N)
        zones(k)%jacobian%dTdZ(1:Nx,0)=zones(West)%jacobian%dTdZ(1:Nx,Nz_N)
        zones(k)%jacobian%dPdR(1:Nx,0)=zones(West)%jacobian%dPdR(1:Nx,Nz_N)
        zones(k)%jacobian%dPdZ(1:Nx,0)=zones(West)%jacobian%dPdZ(1:Nx,Nz_N)
        zones(k)%jacobian%dRdT(1:Nx,0)=zones(West)%jacobian%dRdT(1:Nx,Nz_N)
        zones(k)%jacobian%dRdP(1:Nx,0)=zones(West)%jacobian%dRdP(1:Nx,Nz_N)
        zones(k)%jacobian%dZdT(1:Nx,0)=zones(West)%jacobian%dZdT(1:Nx,Nz_N)
        zones(k)%jacobian%dZdP(1:Nx,0)=zones(West)%jacobian%dZdP(1:Nx,Nz_N)
     end if
  end do
end subroutine MD_broadcast_jacobian
