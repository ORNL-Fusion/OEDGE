subroutine allocate_vorticity()
  use all_variables, only : global_parameters,zones
  implicit none
  integer*4 :: k,Nx,Nz,nstep
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do nstep=1,2
        allocate(zones(k)%electric_fields(nstep)%phi(0:Nx+1,0:Nz+1))
        allocate(zones(k)%electric_fields(nstep)%vorticity(0:Nx+1,0:Nz+1))
     end do
     allocate(zones(k)%electric_fields(1)%phi_smooth(0:Nx+1,0:Nz+1))
     allocate(zones(k)%electric_fields(1)%pi_old(0:Nx+1,0:Nz+1))
     allocate(zones(k)%electric_fields(1)%current(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%SW(1:Nx,1:Nz))
     allocate(zones(k)%electric_fields(1)%FluxW(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%j_parallel(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%j_para_adv_W(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%j_perp_adv_W(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%j_diff_W(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%j_perp(1:Nx,1:Nz,1:4))
     allocate(zones(k)%electric_fields(1)%RHS(1:Nx,1:Nz))
  end do
end subroutine allocate_vorticity
