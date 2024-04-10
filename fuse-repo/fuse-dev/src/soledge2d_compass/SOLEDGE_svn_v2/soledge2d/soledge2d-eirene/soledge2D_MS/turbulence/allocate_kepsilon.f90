subroutine allocate_kepsilon
  use all_variables , only : global_parameters, zones
  use Mturbulence
  use Mdefinitions
  implicit none
  integer*4 :: k, Nx, Nz
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(zones(k)%kepsilon(STEP_OLD)%k(0:Nx+1,0:Nz+1))
     allocate(zones(k)%kepsilon(STEP_OLD)%epsilon(0:Nx+1,0:Nz+1))
     allocate(zones(k)%kepsilon(STEP_OLD)%mu_t(0:Nx+1,0:Nz+1))
     allocate(zones(k)%kepsilon(STEP_OLD)%Sk(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%Sepsilon(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%implicit_coefs%west_k(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%implicit_coefs%east_k(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%implicit_coefs%west_epsilon(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%implicit_coefs%east_epsilon(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%interchange(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%kh(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_OLD)%UE_shear(1:Nx,1:Nz))
     allocate(zones(k)%kepsilon(STEP_NEW)%k(0:Nx+1,0:Nz+1))
     allocate(zones(k)%kepsilon(STEP_NEW)%epsilon(0:Nx+1,0:Nz+1))
     allocate(zones(k)%kepsilon(STEP_NEW)%mu_t(0:Nx+1,0:Nz+1))
  end do
end subroutine allocate_kepsilon
