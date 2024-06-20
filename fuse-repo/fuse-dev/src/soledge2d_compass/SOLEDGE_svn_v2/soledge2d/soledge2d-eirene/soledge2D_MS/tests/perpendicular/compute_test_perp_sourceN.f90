subroutine compute_test_perp_sourceN(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,Dp,n0
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  B0=test_B0
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  n0=test_n0
  Dp=transport_parameters%Dn_p(1)
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        test_sources(k,1)%Sn(i,j)=1/r/(R0+r*cos(theta))*(0.2D0*(R0+r*cos(theta))*Dp*n0*cos(theta)&
             *cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.2D0*r*cos(theta)**2&
             *Dp*n0*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)-0.4D0*r*(R0+r*cos(theta))&
             *Dp*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*pi**2/(rmax-rmin)**2)
        test_sources(k,1)%Sn(i,j)=-test_sources(k,1)%Sn(i,j)&
             /(reference_parameters%fields%n0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_perp_sourceN
