subroutine compute_test_parallel_diffusion_sourceT(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_r,test_theta
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,mass,n0,T0,charge
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  n0=test_n0
  T0=test_T0eV
  G0=test_Gamma0
  B0=test_B0
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  mass=zones(k)%species(1)%element%mass
  charge=zones(k)%species(1)%charge
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        !for electrons
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s8 = 3444400000.D0*psi0/a/(R0+r*cos(theta))**2/(B0**2*R0**2/(R0+r*&
             cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*eV**2/m_e*sqrt(1000.D0)&
             *sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**5*&
             T0*sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))
        s9 = -3444400000.D0*psi0/a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*&
             cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)**2*eV**2/m_e*sqrt(1000.D0)&
             *sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**5&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(2*B0**2*R0**2/(R0&
             +r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3&
             *sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2&
             *r*sin(theta))
        s7 = s8+s9
        s6 = s7-0.8611D9*psi0/a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*cos(theta))**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*eV**2/m_e*sqrt(1000.D0)*&
             sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*T0**2&
             *sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2+3444400000.D0*psi0&
             /a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *cos(theta)**2)*eV**2/m_e*sqrt(1000.D0)*sqrt(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**5*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,0)%ST(i,j) = s1*s2/12.D0
        test_sources(k,0)%ST(i,j)=-test_sources(k,0)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)
        !for ions
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s8 = 0.2574D12*psi0/a/(R0+r*cos(theta))**2/(B0**2*R0**2/(R0+r*cos(theta))**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)/sqrt(mass)*eV**2/m_u*sqrt(1000.D0)&
             /charge**2*sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**5&
             *T0*sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))
        s9 = -0.2574D12*psi0/a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*cos(theta))**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)**2/sqrt(mass)*eV**2/m_u*sqrt(1000.D0)&
             /charge**2*sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**5&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(2*B0**2*R0**2&
             /(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*&
             sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s7 = s8+s9
        s6 = s7-0.6435D11*psi0/a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*cos(theta))**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*cos(theta)**2)/sqrt(mass)*eV**2/m_u*sqrt(1000.D0)&
             /charge**2*sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3&
             *T0**2*sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2+0.2574D12&
             *psi0/a/(R0+r*cos(theta))/r/(B0**2*R0**2/(R0+r*cos(theta)&
             )**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/&
             (R0+r*cos(theta))**2*cos(theta)**2)/sqrt(mass)*eV**2/m_u*sqrt(1000.D0&
             )/charge**2*sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**5&
             *T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,1)%ST(i,j) = s1*s2/12.D0
        test_sources(k,1)%ST(i,j)=-test_sources(k,1)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_parallel_diffusion_sourceT
