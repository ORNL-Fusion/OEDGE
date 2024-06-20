subroutine compute_test_advection_sourceG(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,mass,n0,T0
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  n0=test_n0
  T0=test_T0eV*eV/m_u*2.d0 ! 2 because Te=Ti
  G0=test_Gamma0
  B0=test_B0
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  mass=zones(k)%species(1)%element%mass
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        s2 = -psi0/a
        s4 = 1/(R0+r*cos(theta))
        s6 = 1/r
        s8 = -2*G0**2*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/sqrt(B0**2*R0**2/&
             (R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*sin(theta)
        s9 = 0.1D0*G0**2*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**3/(&
             n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2/sqrt(B0**2&
             *R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*n0*sin(theta)&
             -G0**2*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0&
             +0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/sqrt(B0**2*R0**2&
             /(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+&
             psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)**3*(2*B0**2*&
             R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3&
             *sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2&
             *r*sin(theta))/2
        s7 = s8+s9
        s5 = s6*s7
        s3 = s4*s5
        s1 = s2*s3
        s2 = -psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(theta)&
             )**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/&
             (R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))/mass-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/mass)
        test_sources(k,1)%SG(i,j) = s1+s2
        test_sources(k,1)%SG(i,j)=test_sources(k,1)%SG(i,j)&
             /(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_advection_sourceG
