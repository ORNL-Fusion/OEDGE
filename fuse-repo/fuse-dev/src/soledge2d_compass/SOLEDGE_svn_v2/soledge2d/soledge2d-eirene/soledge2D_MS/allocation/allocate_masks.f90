subroutine allocate_masks()
  use all_variables, only : global_parameters,zones
  implicit none
  integer*4 :: k,Nx,Nz
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(zones(k)%masks%chi(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi1(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi2(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi3(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi4(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi5(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%chi6(0:Nx+1,0:Nz+1))
     allocate(zones(k)%masks%npts_around_penwall(1:Nx,1:Nz))
     allocate(zones(k)%masks%pts_around_penwall(1:Nx,1:Nz,1:4,1:3))
     zones(k)%masks%chi=0.D0
     zones(k)%masks%chi1=0.D0
     zones(k)%masks%chi2=0.D0
     zones(k)%masks%chi3=0.D0
     zones(k)%masks%chi4=0.D0
     zones(k)%masks%chi5=0.D0
     zones(k)%masks%chi6=0.D0
  end do
end subroutine allocate_masks
