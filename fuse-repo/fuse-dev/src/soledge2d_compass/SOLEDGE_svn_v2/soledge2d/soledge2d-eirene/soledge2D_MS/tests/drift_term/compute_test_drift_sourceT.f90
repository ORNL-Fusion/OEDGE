subroutine compute_test_drift_sourceT(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: psi,theta,thetashift
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,n0,phi0,T0,mass
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11
  real*8 :: s12,s13,s14,s15,s16,s17,s18,s19,s20,s21,s22,s23,s24
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  B0=test_B0
  phi0 = test_phi0
  G0=test_Gamma0
  rmax=test_rmax
  rmin=test_rmin
  thetashift = test_theta_shift
  a=test_a
  n0=test_n0
  T0=test_T0ev
  mass=zones(k)%species(1)%element%mass*m_u
  do i=1,Nx
     do j=1,Nz
        psi=zones(k)%mesh%x(i,j)
        theta=2.d0*pi*zones(k)%mesh%z(i,j)

        !ions 
	! ExB drift
	s1 = 1/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))
      	s4 = a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
      	s5 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
      	s3 = s4+s5
      	s7 = 0.1D0*(0.5D0*T0*cos(theta)*cos(2*pi*psi)*pi*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+0.5D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*n0*cos(theta)*&
		cos(2*pi*psi)*pi+2*mass*G0**2*cos(theta)**2*sin(2*pi*psi)/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*cos(2*pi*psi)*pi-0.1D0*mass*G0**2*&
		cos(theta)**3*sin(2*pi*psi)**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))**2*n0*cos(2*pi*psi)*pi)*phi0*sin(theta)
      	s8 = sin(2*pi*psi)*B0*R0/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		cos(theta+psi*thetashift*pi)**2)
      	s6 = s7*s8
      	s8 = -0.1D0*(5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*phi0*sin(theta)
      	s10 = sin(2*pi*psi)*B0
      	s12 = R0
      	s14 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
      	s15 = -2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		sin(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-&
		2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
      	s13 = s14*s15
      	s11 = s12*s13
      	s9 = s10*s11
      	s7 = s8*s9
      	s5 = s6+s7
      	s6 = s5
      	s9 = 0.2D0*(-0.25D0*T0*sin(theta)*sin(2*pi*psi)*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))-0.25D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*n0*&
		sin(theta)*sin(2*pi*psi)-mass*G0**2*cos(theta)*sin(2*pi*psi)**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*sin(theta)+0.5D-1*mass*G0**2*&
		cos(theta)**2*sin(2*pi*psi)**3/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))**2*n0*sin(theta))*phi0*cos(theta)
      	s10 = cos(2*pi*psi)*pi*B0*R0/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)&
		*cos(theta+psi*thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		cos(theta+psi*thetashift*pi)**2)
      	s8 = s9*s10
      	s10 = -0.2D0*(5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*phi0*cos(theta)*cos(2*pi*psi)
      	s11 = pi*B0*R0/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2*&
		(2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))
      	s9 = s10*s11
      	s7 = s8+s9
      	s4 = s6+s7
      	s2 = s3*s4
	test_sources(k,1)%ST(i,j) =s1*s2

	! BxgradB drift
	s1 = 1/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))
	s4 = a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
	s5 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
	s3 = s4+s5
	s7 = -(0.5D0*T0*cos(theta)*cos(2*pi*psi)*pi*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+0.5D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*n0*cos(theta)*&
		cos(2*pi*psi)*pi+2*mass*G0**2*cos(theta)**2*sin(2*pi*psi)/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*cos(2*pi*psi)*pi-0.1D0*mass*G0**2*&
		cos(theta)**3*sin(2*pi*psi)**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))**2*n0*cos(2*pi*psi)*pi)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))
	s8 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)*&
		(2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))*B0*R0
	s6 = s7*s8
	s9 = -0.2D0*(5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*T0*cos(theta)*cos(2*pi*psi)
	s10 = pi/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2*&
		(2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))*B0*R0
	s8 = s9*s10
	s10 = -(5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))
	s13 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
	s16 = -6*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*&
		cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)+0.16D1*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*a*sin(theta+psi*thetashift*pi)
	s15 = s16+2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi-&
		6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)*(0.8D0*a*&
		cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
	s16 = s15+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi)*thetashift*pi+0.16D1*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		sin(theta+psi*thetashift*pi)**3
	s14 = s16-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)+0.16D1*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*sin(theta+psi*thetashift*pi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)*thetashift*pi
	s12 = s13*s14
	s13 = B0*R0
	s11 = s12*s13
	s9 = s10*s11
	s7 = s8+s9
	s5 = s6+s7
	s7 = s5
	s10 = -0.25D0*T0*sin(theta)*sin(2*pi*psi)*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))-0.25D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*n0*sin(theta)*sin(2*pi*psi)-&
		mass*G0**2*cos(theta)*sin(2*pi*psi)**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*sin(theta)+0.5D-1*mass*G0**2*cos(theta)**2*sin(2*pi*psi)**3/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))**2*n0*sin(theta)
	s11 = (T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		cos(theta+psi*thetashift*pi)**2)**2
	s9 = s10*s11
	s10 = (-2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi))*B0*R0
	s8 = s9*s10
	s6 = s7+s8
	s7 = s6
	s10 = -0.1D0*(5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*T0*sin(theta)
	s12 = sin(2*pi*psi)/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2
	s13 = (-2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi))*B0*R0
	s11 = s12*s13
	s9 = s10*s11
	s11 = (5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))+mass*G0**2*cos(theta)**2*sin(2*pi*psi)**2/&
		(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))/2)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+&
		psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2
	s14 = -6*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		(-0.8D0*a*sin(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**4*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)
	s15 = s14-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*(-0.8D0*a*sin(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
	s13 = s15-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(-0.8D0*a*sin(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi)*thetashift*pi)
	s14 = B0*R0
	s12 = s13*s14
	s10 = s11*s12
	s8 = s9+s10
	s4 = s7+s8
	s2 = s3*s4
	test_sources(k,1)%ST(i,j) = test_sources(k,1)%ST(i,j) + s1*s2 

	! Adim
        test_sources(k,1)%ST(i,j) = test_sources(k,1)%ST(i,j)/(reference_parameters%fields%n0*reference_parameters%fields%T0eV/reference_parameters%fields%tau0)

	!electron
	! ExB drift
  	s1 = 1/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))
  	s4 = a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
      	s5 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
      	s3 = s4+s5
      	s6 = -0.25D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*phi0*sin(theta)
      	s8 = sin(2*pi*psi)*B0
      	s10 = R0
      	s12 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
      	s13 = -2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		sin(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-&
		2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
      	s11 = s12*s13
      	s9 = s10*s11
      	s7 = s8*s9
      	s5 = s6*s7
      	s7 = -0.5D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))*phi0*cos(theta)
      	s9 = cos(2*pi*psi)*pi*B0
      	s10 = R0/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2*&
		(2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))
      	s8 = s9*s10
      	s6 = s7*s8
      	s4 = s5+s6
      	s2 = s3*s4
	test_sources(k,0)%ST(i,j) = s1*s2

	! BxgradB drift
	s1 = 1/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))
        s4 = a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
        s5 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift*pi)+&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
        s3 = s4+s5
        s8 = 0.1D1*T0+0.1D0*T0*cos(theta)*sin(2*pi*psi)
        s10 = n0+0.1D0*n0*cos(theta)*sin(2*pi*psi)
        s11 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)*&
		(2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))
        s9 = s10*s11
        s7 = s8*s9
        s8 = B0*R0*T0*cos(theta)*cos(2*pi*psi)*pi
        s6 = s7*s8
        s9 = 0.5D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))**2*n0*cos(theta)*cos(2*pi*psi)
        s10 = pi/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2*(2*B0**2*R0**2/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi))*B0*R0
        s8 = s9*s10
        s10 = 5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))
        s13 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
        s16 = -6*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)+0.16D1*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*&
		sin(theta+psi*thetashift*pi)
        s15 = s16+2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi-&
		6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*sin(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
        s16 = s15+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi)*thetashift*pi+0.16D1*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**3
        s14 = s16-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)+0.16D1*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*sin(theta+psi*thetashift*pi)+2*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**3*(0.6D0+0.8D0*psi)*thetashift*pi
        s12 = s13*s14
        s13 = B0*R0
        s11 = s12*s13
        s9 = s10*s11
        s7 = s8+s9
        s5 = s6+s7
        s7 = s5
        s10 = 0.5D0*T0+0.5D-1*T0*cos(theta)*sin(2*pi*psi)
        s12 = n0+0.1D0*n0*cos(theta)*sin(2*pi*psi)
        s14 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
        s15 = -2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*&
		sin(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-&
		2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
        s13 = s14*s15
        s11 = s12*s13
        s9 = s10*s11
        s10 = B0*R0*T0*sin(theta)*sin(2*pi*psi)
        s8 = s9*s10
        s6 = s7+s8
        s7 = s6
        s10 = 0.25D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))**2*n0*sin(theta)
        s12 = sin(2*pi*psi)/(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		cos(theta+psi*thetashift*pi)**2)**2
        s13 = (-2*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*&
		cos(theta+psi*thetashift*pi)**2*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi))*B0*R0
        s11 = s12*s13
        s9 = s10*s11
        s11 = -5.D0/2.D0*(T0+0.1D0*T0*cos(theta)*sin(2*pi*psi))**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*psi))
        s14 = 1/((B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*&
		sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi*thetashift*pi)**2)**2)
        s16 = -6*B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*B0**2*R0**2/&
		(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*(-0.8D0*a*sin(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi)*thetashift*pi)-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*sin(theta+psi*thetashift*pi)**3*&
		(0.6D0+0.8D0*psi)*(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
        s17 = s16-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**3*sin(theta+psi*thetashift*pi)**2*(-0.8D0*a*sin(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
        s15 = s17-6*psi0**2/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**4*cos(theta+psi*thetashift*pi)**2*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*&
		(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)-2*psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*&
		cos(theta+psi*thetashift*pi))**3*cos(theta+psi*thetashift*pi)**2*(-0.8D0*a*sin(theta+psi*thetashift*pi)-&
		a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
        s13 = s14*s15
        s14 = B0*R0
        s12 = s13*s14
        s10 = s11*s12
        s8 = s9+s10
        s4 = s7+s8
        s2 = s3*s4
      	test_sources(k,0)%ST(i,j) = test_sources(k,0)%ST(i,j) + s1*s2
        ! Adim
        test_sources(k,0)%ST(i,j)= test_sources(k,0)%ST(i,j)/(reference_parameters%fields%n0*reference_parameters%fields%T0eV/reference_parameters%fields%tau0)

     end do
  end do
end subroutine compute_test_drift_sourceT
