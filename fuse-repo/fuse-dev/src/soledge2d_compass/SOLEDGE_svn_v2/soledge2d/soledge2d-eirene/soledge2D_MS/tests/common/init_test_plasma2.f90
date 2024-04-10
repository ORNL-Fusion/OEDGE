subroutine init_test_plasma2()
  use test_var
  use all_variables, only : zones, global_parameters,reference_parameters
  use Mphysics
  implicit none
  integer*4 :: k,n
  real*8 :: r,theta,test_theta,test_r
  integer*4 :: i,j,Nx,Nz
  do k=1,global_parameters%N_Zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do i=1,Nx
        do j=1,Nz
           !ion number 1 (He+)
           r=test_r(zones(k)%mesh%x(i,j))
           theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
           zones(k)%species(1)%var(1)%density(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin)))*test_n0/reference_parameters%fields%n0
           zones(k)%species(1)%var(1)%Gamma(i,j)=cos(theta)*sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin))&
                *test_Gamma0/reference_parameters%fields%n0/reference_parameters%fields%c0
           zones(k)%species(1)%var(1)%temperature(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin)))*test_T0eV/reference_parameters%fields%T0eV
           !ion number 2 (He2+)
           zones(k)%species(2)%var(1)%density(i,j)=0.5d0*(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin)))*test_n0/reference_parameters%fields%n0
           zones(k)%species(2)%var(1)%Gamma(i,j)=cos(theta)*sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin))&
                *test_Gamma0/reference_parameters%fields%n0/reference_parameters%fields%c0
           zones(k)%species(2)%var(1)%temperature(i,j)=2.d0*(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin)))*test_T0eV/reference_parameters%fields%T0eV
           !electrons
           zones(k)%species(0)%var(1)%temperature(i,j)=(1.D0+0.1D0*cos(theta)&
                *sin(2.D0*pi*(r-test_rmin)/(test_rmax-test_rmin)))*test_T0eV/reference_parameters%fields%T0eV
        end do
     end do
  end do
end subroutine init_test_plasma2
