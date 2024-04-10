subroutine MD_broadcast_masks()
  use all_variables, only : zones, global_parameters
  implicit none
  integer*4 :: k
  integer*4 :: Nx,Nz,Nx_N,Nz_N
  integer*4 :: North,South,East,West
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     !North boundary
     North=zones(k)%Neighbors(1)
     if(North.lt.0) then
        zones(k)%masks%chi1(Nx+1,:)=zones(k)%masks%chi1(Nx,:)
        zones(k)%masks%chi2(Nx+1,:)=zones(k)%masks%chi2(Nx,:)
        zones(k)%masks%chi3(Nx+1,:)=zones(k)%masks%chi3(Nx,:)
        zones(k)%masks%chi4(Nx+1,:)=zones(k)%masks%chi4(Nx,:)
     else
        zones(k)%masks%chi1(Nx+1,1:Nz)=zones(North)%masks%chi1(1,1:Nz)
        zones(k)%masks%chi2(Nx+1,1:Nz)=zones(North)%masks%chi2(1,1:Nz)
        zones(k)%masks%chi3(Nx+1,1:Nz)=zones(North)%masks%chi3(1,1:Nz)
        zones(k)%masks%chi4(Nx+1,1:Nz)=zones(North)%masks%chi4(1,1:Nz)
     end if
     !South boundary
     South=zones(k)%Neighbors(2)
     if(South.lt.0) then
        zones(k)%masks%chi1(0,:)=zones(k)%masks%chi1(1,:)
        zones(k)%masks%chi2(0,:)=zones(k)%masks%chi2(1,:)
        zones(k)%masks%chi3(0,:)=zones(k)%masks%chi3(1,:)
        zones(k)%masks%chi4(0,:)=zones(k)%masks%chi4(1,:)
     else
        Nx_N=zones(South)%mesh%Nx
        zones(k)%masks%chi1(0,1:Nz)=zones(South)%masks%chi1(Nx_N,1:Nz)
        zones(k)%masks%chi2(0,1:Nz)=zones(South)%masks%chi2(Nx_N,1:Nz)
        zones(k)%masks%chi3(0,1:Nz)=zones(South)%masks%chi3(Nx_N,1:Nz)
        zones(k)%masks%chi4(0,1:Nz)=zones(South)%masks%chi4(Nx_N,1:Nz)
     end if
     !East boundary
     East=zones(k)%Neighbors(3)
     if(East.lt.0) then
        zones(k)%masks%chi1(1:Nx,Nz+1)=zones(k)%masks%chi1(1:Nx,Nz)
        zones(k)%masks%chi2(1:Nx,Nz+1)=zones(k)%masks%chi2(1:Nx,Nz)
        zones(k)%masks%chi3(1:Nx,Nz+1)=zones(k)%masks%chi3(1:Nx,Nz)
        zones(k)%masks%chi4(1:Nx,Nz+1)=zones(k)%masks%chi4(1:Nx,Nz)
     else
        zones(k)%masks%chi1(1:Nx,Nz+1)=zones(East)%masks%chi1(1:Nx,1)
        zones(k)%masks%chi2(1:Nx,Nz+1)=zones(East)%masks%chi2(1:Nx,1)
        zones(k)%masks%chi3(1:Nx,Nz+1)=zones(East)%masks%chi3(1:Nx,1)
        zones(k)%masks%chi4(1:Nx,Nz+1)=zones(East)%masks%chi4(1:Nx,1)
     end if
     !West boundary
     West=zones(k)%Neighbors(4)
     if(West.lt.0) then
        zones(k)%masks%chi1(1:Nx,0)=zones(k)%masks%chi1(1:Nx,1)
        zones(k)%masks%chi2(1:Nx,0)=zones(k)%masks%chi2(1:Nx,1)
        zones(k)%masks%chi3(1:Nx,0)=zones(k)%masks%chi3(1:Nx,1)
        zones(k)%masks%chi4(1:Nx,0)=zones(k)%masks%chi4(1:Nx,1)
     else
        Nz_N=zones(West)%mesh%Nz
        zones(k)%masks%chi1(1:Nx,0)=zones(West)%masks%chi1(1:Nx,Nz_N)
        zones(k)%masks%chi2(1:Nx,0)=zones(West)%masks%chi2(1:Nx,Nz_N)
        zones(k)%masks%chi3(1:Nx,0)=zones(West)%masks%chi3(1:Nx,Nz_N)
        zones(k)%masks%chi4(1:Nx,0)=zones(West)%masks%chi4(1:Nx,Nz_N)
     end if
  end do
end subroutine MD_broadcast_masks
