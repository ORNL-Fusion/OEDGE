subroutine compute_test_parallel_diffusion_sourceT(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,mass,n0,T0
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
  mass=zones(k)%species(1)%element%mass*m_u
  do i=1,Nx
     do j=1,Nz
        r=test_a*(0.6D0+0.8D0*zones(k)%mesh%x(i,j))
        theta=2.D0*pi*zones(k)%mesh%z(i,j)+test_theta_shift*pi*zones(k)%mesh%x(i,j)
        !for electrons
        s1 = -0.5D0/(2+r*cos(theta))
        s3 = 1/r
        s6 = -0.3587916665D12*ev**2*sqrt(100+0.1D2*cos(theta)*sin(2*0.3141592653589793D1&
             *(0.15625D1*r-0.75D0)))**3/m_e*sqrt(1000.D0)/(36/(2+r&
             *cos(theta))**2+0.25D0/(2+r*cos(theta))**2*sin(theta)**2+0.25D0/(2&
             +r*cos(theta))**2*cos(theta)**2)/(2+r*cos(theta))/r*sin(theta)**2*&
             sin(2*0.3141592653589793D1*(0.15625D1*r-0.75D0))**2
        s7 = -0.1435166666D11*ev**2*sqrt(100+0.1D2*cos(theta)*sin(2*0.3141592653589793D1&
             *(0.15625D1*r-0.75D0)))**5/m_e*sqrt(1000.D0)/(36/(2+r&
             *cos(theta))**2+0.25D0/(2+r*cos(theta))**2*sin(theta)**2+0.25D0/(2&
             +r*cos(theta))**2*cos(theta)**2)**2/(2+r*cos(theta))/r*sin(theta)*&
             sin(2*0.3141592653589793D1*(0.15625D1*r-0.75D0))*(72/(2+r*cos(theta))**3&
             *r*sin(theta)+0.5D0/(2+r*cos(theta))**3*sin(theta)**3*r+0.5D0&
             /(2+r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s5 = s6+s7
        s4 = s5+0.1435166666D11*ev**2*sqrt(100+0.1D2*cos(theta)*sin(2*0.3141592653589793D1&
             *(0.15625D1*r-0.75D0)))**5/m_e*sqrt(1000.D0)/(36/(2&
             +r*cos(theta))**2+0.25D0/(2+r*cos(theta))**2*sin(theta)**2+0.25D0/&
             (2+r*cos(theta))**2*cos(theta)**2)/(2+r*cos(theta))**2*sin(theta)*&
             *2*sin(2*0.3141592653589793D1*(0.15625D1*r-0.75D0))+0.1435166666D11&
             *ev**2*sqrt(100+0.1D2*cos(theta)*sin(2*0.3141592653589793D1*(0.15625D1&
             *r-0.75D0)))**5/m_e*sqrt(1000.D0)/(36/(2+r*cos(theta))**2+0.25D0&
             /(2+r*cos(theta))**2*sin(theta)**2+0.25D0/(2+r*cos(theta))**2*cos(theta)**2)&
             /(2+r*cos(theta))/r*cos(theta)*sin(2*0.3141592653589793D1*(0.15625D1*r-0.75D0))
        s2 = s3*s4
        test_sources(k,0)%ST(i,j) = s1*s2
        test_sources(k,0)%ST(i,j)=-test_sources(k,0)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_parallel_diffusion_sourceT
