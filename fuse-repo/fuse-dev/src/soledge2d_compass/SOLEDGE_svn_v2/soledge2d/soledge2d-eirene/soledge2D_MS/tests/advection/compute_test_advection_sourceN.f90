subroutine compute_test_advection_sourceN(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  B0=test_B0
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        test_sources(k,1)%Sn(i,j)=-psi0 / a / (R0 + r * cos(theta)) / r * (-G0 * sin(theta) * sin(0.2D1&
             * pi * (r - rmin) / (rmax - rmin)) * (B0 ** 2 * R0 ** 2 /&
             (R0 + r * cos(theta)) ** 2 + psi0 ** 2 / a ** 2 / (R0 + r * cos(theta))&
             ** 2 * sin(theta) ** 2 + psi0 ** 2 / a ** 2 / (R0 + r * cos(theta)) ** 2&
             * cos(theta) ** 2) ** (-0.1D1 / 0.2D1) - G0 * cos(theta) &
             * sin(0.2D1 * pi * (r - rmin) / (rmax - rmin)) * (B0 ** 2 * R0 ** 2&
             / (R0 + r * cos(theta)) ** 2 + psi0 ** 2 / a ** 2 / (R0 + r&
             * cos(theta)) ** 2 * sin(theta) ** 2 + psi0 ** 2 / a ** 2 / (R0 +&
             r * cos(theta)) ** 2 * cos(theta) ** 2) ** (-0.3D1 / 0.2D1) * (0.2D1&
             * B0 ** 2 * R0 ** 2 / (R0 + r * cos(theta)) ** 3 * r * sin(theta)&
             + 0.2D1 * psi0 ** 2 / a ** 2 / (R0 + r * cos(theta)) ** 3 * sin(theta) ** 3&
             * r + 0.2D1 * psi0 ** 2 / a ** 2 / (R0 + r * cos(theta)) ** 3&
             * cos(theta) ** 2 * r * sin(theta)) / 0.2D1)
        test_sources(k,1)%Sn(i,j)=test_sources(k,1)%Sn(i,j)&
             /(reference_parameters%fields%n0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_advection_sourceN
