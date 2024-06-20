subroutine MD_broadcast_metric()
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
        zones(k)%metric_coefficients%cpp(Nx+1,1:Nz)=zones(North)%metric_coefficients%cpp(1,1:Nz)
        zones(k)%metric_coefficients%cpt(Nx+1,1:Nz)=zones(North)%metric_coefficients%cpt(1,1:Nz)
        zones(k)%metric_coefficients%ctt(Nx+1,1:Nz)=zones(North)%metric_coefficients%ctt(1,1:Nz)
        zones(k)%metric_coefficients%Jacobian(Nx+1,1:Nz)=zones(North)%metric_coefficients%Jacobian(1,1:Nz)
        zones(k)%metric_coefficients%G(Nx+1,1:Nz)=zones(North)%metric_coefficients%G(1,1:Nz)
        zones(k)%metric_coefficients%c_pp(Nx+1,1:Nz)=zones(North)%metric_coefficients%c_pp(1,1:Nz)
        zones(k)%metric_coefficients%c_pt(Nx+1,1:Nz)=zones(North)%metric_coefficients%c_pt(1,1:Nz)
        zones(k)%metric_coefficients%c_tt(Nx+1,1:Nz)=zones(North)%metric_coefficients%c_tt(1,1:Nz)
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     mSouth=zones(k)%MagNeighbors(2)
     if(mSouth.eq.0) then
        Nx_N=zones(South)%mesh%Nx
        zones(k)%metric_coefficients%cpp(0,1:Nz)=zones(South)%metric_coefficients%cpp(Nx_N,1:Nz)
        zones(k)%metric_coefficients%cpt(0,1:Nz)=zones(South)%metric_coefficients%cpt(Nx_N,1:Nz)
        zones(k)%metric_coefficients%ctt(0,1:Nz)=zones(South)%metric_coefficients%ctt(Nx_N,1:Nz)
        zones(k)%metric_coefficients%Jacobian(0,1:Nz)=zones(South)%metric_coefficients%Jacobian(Nx_N,1:Nz)
        zones(k)%metric_coefficients%G(0,1:Nz)=zones(South)%metric_coefficients%G(Nx_N,1:Nz)
        zones(k)%metric_coefficients%c_pp(0,1:Nz)=zones(South)%metric_coefficients%c_pp(Nx_N,1:Nz)
        zones(k)%metric_coefficients%c_pt(0,1:Nz)=zones(South)%metric_coefficients%c_pt(Nx_N,1:Nz)
        zones(k)%metric_coefficients%c_tt(0,1:Nz)=zones(South)%metric_coefficients%c_tt(Nx_N,1:Nz)
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     mEast=zones(k)%MagNeighbors(3)
     if(mEast.eq.0) then
        zones(k)%metric_coefficients%cpp(1:Nx,Nz+1)=zones(East)%metric_coefficients%cpp(1:Nx,1)
        zones(k)%metric_coefficients%cpt(1:Nx,Nz+1)=zones(East)%metric_coefficients%cpt(1:Nx,1)
        zones(k)%metric_coefficients%ctt(1:Nx,Nz+1)=zones(East)%metric_coefficients%ctt(1:Nx,1)
        zones(k)%metric_coefficients%Jacobian(1:Nx,Nz+1)=zones(East)%metric_coefficients%Jacobian(1:Nx,1)
        zones(k)%metric_coefficients%G(1:Nx,Nz+1)=zones(East)%metric_coefficients%G(1:Nx,1)
        zones(k)%metric_coefficients%c_pp(1:Nx,Nz+1)=zones(East)%metric_coefficients%c_pp(1:Nx,1)
        zones(k)%metric_coefficients%c_pt(1:Nx,Nz+1)=zones(East)%metric_coefficients%c_pt(1:Nx,1)
        zones(k)%metric_coefficients%c_tt(1:Nx,Nz+1)=zones(East)%metric_coefficients%c_tt(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     mWest=zones(k)%MagNeighbors(4)
     if(mWest.eq.0) then
        Nz_N=zones(West)%mesh%Nz
        zones(k)%metric_coefficients%cpp(1:Nx,0)=zones(West)%metric_coefficients%cpp(1:Nx,Nz_N)
        zones(k)%metric_coefficients%cpt(1:Nx,0)=zones(West)%metric_coefficients%cpt(1:Nx,Nz_N)
        zones(k)%metric_coefficients%ctt(1:Nx,0)=zones(West)%metric_coefficients%ctt(1:Nx,Nz_N)
        zones(k)%metric_coefficients%Jacobian(1:Nx,0)=zones(West)%metric_coefficients%Jacobian(1:Nx,Nz_N)
        zones(k)%metric_coefficients%G(1:Nx,0)=zones(West)%metric_coefficients%G(1:Nx,Nz_N)
        zones(k)%metric_coefficients%c_pp(1:Nx,0)=zones(West)%metric_coefficients%c_pp(1:Nx,Nz_N)
        zones(k)%metric_coefficients%c_pt(1:Nx,0)=zones(West)%metric_coefficients%c_pt(1:Nx,Nz_N)
        zones(k)%metric_coefficients%c_tt(1:Nx,0)=zones(West)%metric_coefficients%c_tt(1:Nx,Nz_N)
     end if
  end do
end subroutine MD_broadcast_metric
