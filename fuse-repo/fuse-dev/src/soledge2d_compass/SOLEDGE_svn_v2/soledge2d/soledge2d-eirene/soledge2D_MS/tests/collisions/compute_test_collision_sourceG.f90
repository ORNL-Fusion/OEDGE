subroutine compute_test_collision_sourceG(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,Dp,n0,nup,T0,loglambda
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
  real*8 :: m0,m1,m2
  real*8 :: q0,q1,q2
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
  T0=test_T0ev
  Dp=transport_parameters%Dn_p(1)
  nup=transport_parameters%nu_p(1)
  m0=zones(1)%species(0)%element%mass2*m_u
  m1=zones(1)%species(1)%element%mass2*m_u
  m2=zones(1)%species(2)%element%mass2*m_u
  q0=zones(1)%species(0)%charge
  q1=zones(1)%species(1)%charge
  q2=zones(1)%species(2)%charge
  loglambda=12.d0
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))

        !first  ion He+

        s3 = q1*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s8 = -0.425D-1*(m0/m1+m1/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s9 = 1/m1*(0.15D1*n0+0.15D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*eV**4*q0**2*q1**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(T0+0.1D0*T0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*eV/m1)**3*(3*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)
        s7 = s8*s9
        s9 = -0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s10 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)&
             *(0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*&
             cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))-0.1D0*q1**2/m1*psi0/a/(R0+r*cos(theta))&
             /r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s8 = s9*s10
        s6 = s7+s8
        s8 = s6
        s10 = -0.425D-1*(m0/m2+m2/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s11 = 1/m2*(0.125D1*n0+0.125D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV**4*q0**2*q2**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(2*T0+0.2D0*T0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))*eV/m2)**3*(3*G0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1&
             *n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)
        s9 = s10*s11
        s7 = s8+s9
        s8 = s7
        s11 = -0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s12 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m2)*(0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2&
             /(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))-0.2D0*q2**2/m2*psi0/a/(R0+r*cos(&
             theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r&
             *cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s10 = s11*s12
        s11 = psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(theta)&
             )**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/&
             (R0+r*cos(theta))**2*cos(theta)**2)*(-0.2D0*n0*sin(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV-0.1D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*eV)
        s9 = s10+s11
        s5 = s8+s9
        s6 = 1/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s4 = s5*s6
        s2 = s3*s4
        s4 = 0.425D-1*(m0/m1+m1/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s5 = 1/m1*(0.15D1*n0+0.15D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*eV**4*q0**2*q1**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(T0+0.1D0*T0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin)))*eV/m1)**3*(3*G0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)
        s3 = s4*s5
        s1 = s2+s3
        s3 = s1
        s5 = 0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s6 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)*&
             (0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r&
             *cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))-0.1D0*q1**2/m1*psi0/a/(R0+r*cos(theta))/&
             r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s4 = s5*s6
        s2 = s3+s4
        s3 = s2
        s6 = -0.425D-1*(m1/m2+m2/m1)*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m1
        s7 = 1/m2*(0.75D0*n0+0.75D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV**4*q1**2*q2**2*loglambda/sqrt((T0+0.1D0*T0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))*eV/m1+(2*T0+0.2D0*T0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin)))*eV/m2)**3*(G0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0*n0+0.5D-1*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)
        s5 = s6*s7
        s7 = -0.71D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s8 = 1/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1&
             +(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)&
             *(0.1D0*q1**2/m1*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+&
             r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))-0.2D0*q2**2/m2*psi0/a/(R0+r*cos(theta))&
             /r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s6 = s7*s8
        s4 = s5+s6

        test_sources(k,1)%SG(i,j) = s3+s4
        test_sources(k,1)%SG(i,j)=-test_sources(k,1)%SG(i,j)&
             /(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)
        test_sources(k,1)%SG(i,j)=test_sources(k,1)%SG(i,j)/(m_u*zones(k)%species(1)%element%mass2)

        !second ion He2+

        s3 = q2*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s8 = -0.425D-1*(m0/m1+m1/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s9 = 1/m1*(0.15D1*n0+0.15D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*eV**4*q0**2*q1**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(T0+0.1D0*T0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*eV/m1)**3*(3*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)
        s7 = s8*s9
        s9 = -0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s10 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)&
             *(0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r&
             *cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin))-0.1D0*q1**2/m1*psi0/a/(R0+r*cos(theta))&
             /r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s8 = s9*s10
        s6 = s7+s8
        s8 = s6
        s10 = -0.425D-1*(m0/m2+m2/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s11 = 1/m2*(0.125D1*n0+0.125D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV**4*q0**2*q2**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(2*T0+0.2D0*T0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))*eV/m2)**3*(3*G0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1&
             *n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)
        s9 = s10*s11
        s7 = s8+s9
        s8 = s7
        s11 = -0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s12 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m2)*(0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2&
             /(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))-0.2D0*q2**2/m2*psi0/a/(R0+r*cos(theta)&
             )/r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s10 = s11*s12
        s11 = psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+&
             r*cos(theta))**2*cos(theta)**2)*(-0.2D0*n0*sin(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))*eV-0.1D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*eV)
        s9 = s10+s11
        s5 = s8+s9
        s6 = 1/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s4 = s5*s6
        s2 = s3*s4
        s4 = 0.425D-1*(m0/m2+m2/m0)*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m0
        s5 = 1/m2*(0.125D1*n0+0.125D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV**4*q0**2*q2**2*loglambda/sqrt(eV*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(2*T0+0.2D0*T0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))*eV/m2)**3*(3*G0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*&
             n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)
        s3 = s4*s5
        s1 = s2+s3
        s3 = s1
        s5 = 0.71D0*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s6 = 1/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))/m2)*(0.1D0*q0**2/m0*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2&
             /(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2&
             +psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))-0.2D0*q2**2/m2*psi0/a/(R0+r*cos(theta)&
             )/r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*&
             cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s4 = s5*s6
        s2 = s3+s4
        s3 = s2
        s6 = 0.425D-1*(m1/m2+m2/m1)*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*sqrt(2.D0)/sqrt(0.3141592653589793D1)**3/epsilon_0**2/m1
        s7 = 1/m2*(0.75D0*n0+0.75D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*eV**4*q1**2*q2**2*loglambda/sqrt((T0+0.1D0*T0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))*eV/m1+(2*T0+0.2D0*T0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))*eV/m2)**3*(G0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((n0+0.1D0*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0*n0+0.5D-1*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)
        s5 = s6*s7
        s7 = 0.71D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s8 = 1/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1&
             +(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2&
             )*(0.1D0*q1**2/m1*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+&
             r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*T0*sin(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))-0.2D0*q2**2/m2*psi0/a/(R0+r*cos(theta))&
             /r/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*eV
        s6 = s7*s8
        s4 = s5+s6

        test_sources(k,2)%SG(i,j) = s3+s4
        test_sources(k,2)%SG(i,j)=-test_sources(k,2)%SG(i,j)&
             /(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)
        test_sources(k,2)%SG(i,j)=test_sources(k,2)%SG(i,j)/(m_u*zones(k)%species(2)%element%mass2)

     end do
  end do
end subroutine compute_test_collision_sourceG
