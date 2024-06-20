subroutine allocate_geometry()
  use all_variables, only : global_parameters,zones
  implicit none
  integer*4 :: k,Nx,Nz
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(zones(k)%mesh%Rgeom(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%Zgeom(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%Rcorner(1:Nx+1,1:Nz+1))
     allocate(zones(k)%mesh%Zcorner(1:Nx+1,1:Nz+1))
     allocate(zones(k)%mesh%Br(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%Bz(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%Bphi(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%B(0:Nx+1,0:Nz+1))
     allocate(zones(k)%mesh%index(0:Nx+1,0:Nz+1))
     !set to non zero to avoid NaN
     zones(k)%mesh%Rgeom(0:Nx+1,0:Nz+1)=1.d-15
     zones(k)%mesh%Zgeom(0:Nx+1,0:Nz+1)=1.d-15
     zones(k)%mesh%Rcorner(1:Nx+1,1:Nz+1)=1.d-15
     zones(k)%mesh%Zcorner(1:Nx+1,1:Nz+1)=1.d-15
     zones(k)%mesh%Br(0:Nx+1,0:Nz+1)=1.d-15
     zones(k)%mesh%Bz(0:Nx+1,0:Nz+1)=1.d-15
     zones(k)%mesh%Bphi(0:Nx+1,0:Nz+1)=1.d-15
     zones(k)%mesh%B(0:Nx+1,0:Nz+1)=1.d-15
  end do
end subroutine allocate_geometry
