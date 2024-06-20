subroutine compute_total_radiation()
  use all_variables, only : globals, global_parameters, reference_parameters, zones
  use Mphysics
  implicit none
  integer*4 :: n, Nx, Nz, k
  do n=0,global_parameters%N_ions
     globals%total_radiation(n)=0.D0
  end do
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do n=0,global_parameters%N_ions
        globals%total_radiation(n)=globals%total_radiation(n)&
             +sum(zones(k)%species(n)%sources%rad(1:Nx,1:Nz)*&
             zones(k)%metric_coefficients%dvol_pu(1:Nx,1:Nz)&
             *reference_parameters%fields%n0*reference_parameters%fields%T0eV*eV&
             /reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_total_radiation
