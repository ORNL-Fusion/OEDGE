subroutine init_test_plasma_drift()
  use test_var
  use all_variables, only : zones, global_parameters,reference_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  real*8 :: psi,theta
  integer*4 :: i,j,Nx,Nz
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           !main ions
           psi=zones(k)%mesh%x(i,j)
           theta=2.d0*pi*zones(k)%mesh%z(i,j)
           zones(k)%species(1)%var(1)%density(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*psi))*test_n0/reference_parameters%fields%n0
           zones(k)%species(1)%var(1)%Gamma(i,j)=cos(theta)*sin(2.D0*pi*psi)&
                *test_Gamma0/reference_parameters%fields%n0/reference_parameters%fields%c0
           zones(k)%species(1)%var(1)%temperature(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*psi))*test_T0eV/reference_parameters%fields%T0eV
           !electrons
           zones(k)%species(0)%var(1)%temperature(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*psi))*test_T0eV/reference_parameters%fields%T0eV
           !potential
           zones(k)%electric_fields(1)%phi(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*psi))*test_phi0/reference_parameters%fields%T0eV
        end do
     end do
  end do
end subroutine init_test_plasma_drift
