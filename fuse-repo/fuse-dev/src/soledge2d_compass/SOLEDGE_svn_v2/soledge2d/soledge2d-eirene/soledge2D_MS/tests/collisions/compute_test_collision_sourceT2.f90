subroutine compute_test_collision_sourceT2(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,Dp,n0,nup,chip
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
  real*8 :: s14,s15,s16,s17,s18,s19,s20
  real*8 :: mass,T0
  real*8 :: m0,m1,m2
  real*8 :: q0,q1,q2,loglambda
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  B0=test_B0
  T0=test_T0ev
  mass=zones(k)%species(1)%element%mass*m_u
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  n0=test_n0
  Dp=transport_parameters%Dn_p(1)
  nup=transport_parameters%nu_p(1)
  chip=transport_parameters%chii_p(1)
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

        !ion number 1 He+
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s11 = -0.71D-1*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q1**2&
             *(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(&
             n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             ))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0&
             +(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m1
        s12 = -0.142D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q1**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/&
             (0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0&
             +0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(n0+&
             0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m1
        s10 = s11+s12
        s11 = s10-0.71D-1*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q1**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m1
        s12 = s11
        s15 = 0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*q1**2
        s17 = (0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.1D0*G0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))**2*n0*sin(theta)+3*G0*sin(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-0.6D0*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))**2&
             /(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2&
             *n0*sin(theta))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m1)/m1
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q1**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)**2&
             /m1*(-0.2D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m0-0.1D0&
             *n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m1)
        s15 = s16*s17
        s13 = s14+s15
        s9 = s12+s13
        s11 = s9
        s13 = -0.71D-1*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q1**2&
             *(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             ))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m1
        s14 = -0.71D-1*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q1**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             ))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m1
        s12 = s13+s14
        s10 = s11+s12
        s11 = s10-0.355D-1*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q1**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(G0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/&
             ((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m1
        s12 = s11
        s15 = 0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q1**2
        s17 = (n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.1D0*G0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))**2*n0*sin(theta)+G0*sin(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))-0.5D-1*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))**2&
             /(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2&
             *n0*sin(theta))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m2)/m1
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q1**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m2)**2/m1*(-0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/m1-0.5D-1*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m2)
        s15 = s16*s17
        s13 = s14+s15
        s8 = s12+s13
        s9 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)
        s7 = s8*s9
        s10 = -0.355D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q1**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*&
             (G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))/m1)/m1
        s11 = -0.355D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q1**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*&
             (G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m2)/m1
        s9 = s10+s11
        s10 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(&
             theta)**2)**3*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2&
             /a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0&
             +r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s8 = s9*s10
        s6 = s7+s8
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,1)%ST(i,j)=s1*s2
        test_sources(k,1)%ST(i,j)=test_sources(k,1)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)

        !ion number 2
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s11 = -0.142D0*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q2**2&
             *(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m2)/m2
        s12 = -0.142D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q2**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m2)/m2
        s10 = s11+s12
        s11 = s10-0.355D-1*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))*q2**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))/m2)/m2
        s12 = s11
        s15 = 0.71D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q2**2
        s17 = (0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.5D-1*G0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin))**2/(0.5D0*n0+0.5D-1*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*sin(theta)+3*G0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))-0.6D0*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))**2/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))**2*n0*sin(theta))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q2**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))/m2)**2/m2*(-0.2D0*n0*sin(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/m0-0.5D-1*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))/m2)
        s15 = s16*s17
        s13 = s14+s15
        s9 = s12+s13
        s11 = s9
        s13 = -0.142D0*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q2**2&
             *(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             ))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s14 = -0.71D-1*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q2**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(0.5D0&
             *n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             ))))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(&
             0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s12 = s13+s14
        s10 = s11+s12
        s11 = s10-0.355D-1*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))*q2**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))&
             )/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1+(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s12 = s11
        s15 = 0.71D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*q2**2
        s17 = (n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.5D-1*G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))**2/(0.5D0*n0+0.5D-1*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*sin(theta)+G0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))-0.1D0*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2&
             *n0*sin(theta))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m2)/m2
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q2**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))/m2)**2/m2*(-0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/m1-0.5D-1*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m2)
        s15 = s16*s17
        s13 = s14+s15
        s8 = s12+s13
        s9 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)
        s7 = s8*s9
        s10 = -0.355D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q2**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-3*G0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s11 = -0.355D0*eV*(2*T0+0.2D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q2**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))))/((n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m1+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m2)/m2
        s9 = s10+s11
        s10 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(&
             theta)**2)**3*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2&
             /a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+&
             r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s8 = s9*s10
        s6 = s7+s8
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,2)%ST(i,j)=s1*s2
        test_sources(k,2)%ST(i,j)=test_sources(k,2)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)

        !electrons
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s11 = -0.71D-1*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q0**2&
             *(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(&
             n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(3*G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             ))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/&
             m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m0
        s12 = -0.142D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q0**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(3*G0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))))/((0.2D1*&
             n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m0+(n0&
             +0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m0
        s10 = s11+s12
        s11 = s10-0.71D-1*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(3*G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m0
        s12 = s11
        s15 = 0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*q0**2
        s17 = (0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-3*G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.6D0*G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))**2/(0.2D1*n0+0.2D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*sin(theta)+G0*sin(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))-0.1D0*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             **2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2&
             *n0*sin(theta))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m0
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(3*G0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/&
             m1)**2/m0*(-0.2D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m0-&
             0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/m1)
        s15 = s16*s17
        s13 = s14+s15
        s9 = s12+s13
        s11 = s9
        s13 = -0.71D-1*eV*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*q0**2&
             *(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(&
             0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(3*G0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))/m2)/m0
        s14 = -0.142D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q0**2*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(0.5D0*&
             n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(3*G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))/m2)/m0
        s12 = s13+s14
        s10 = s11+s12
        s11 = s10-0.355D-1*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(3*G0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m2)/m0
        s12 = s11
        s15 = 0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*q0**2
        s17 = (0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))&
             *(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s18 = (-3*G0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.6D0*G0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))**2/(0.2D1*n0+0.2D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*sin(theta)+G0*sin(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))-0.5D-1*G0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))**2/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))**2*n0*sin(theta))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))/m2)/m0
        s16 = s17*s18
        s14 = s15*s16
        s16 = -0.71D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s17 = (0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))/m2)**2/m0*(-0.2D0*n0*sin(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/m0-0.5D-1*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))/m2)
        s15 = s16*s17
        s13 = s14+s15
        s8 = s12+s13
        s9 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)
        s7 = s8*s9
        s10 = -0.355D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*&
             (3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))/m0+(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/m1)/m0
        s11 = -0.355D0*eV*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*q0**2*(0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))*(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*(3*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))/(0.2D1*n0+0.2D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-G0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin))/(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2&
             *pi*(r-rmin)/(rmax-rmin))))/((0.2D1*n0+0.2D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))/m0+(0.5D0*n0+0.5D-1*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))/m2)/m0
        s9 = s10+s11
        s10 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(&
             theta)**2)**3*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2&
             /a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+&
             r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s8 = s9*s10
        s6 = s7+s8
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,0)%ST(i,j)=s1*s2
        test_sources(k,0)%ST(i,j)=test_sources(k,0)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_collision_sourceT2
