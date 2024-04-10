subroutine allocate_neutrals()
  use all_variables, only : global_parameters,zones
  implicit none
  integer*4 :: k,Nx,Nz,nstep
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(zones(k)%neutrals%density(0:Nx+1,0:Nz+1))
     allocate(zones(k)%neutrals%Dn(0:Nx+1,0:Nz+1))
     allocate(zones(k)%neutrals%RHS(1:Nx,1:Nz))
     allocate(zones(k)%neutrals%Sn_nn(1:Nx,1:Nz))
  end do
end subroutine allocate_neutrals
