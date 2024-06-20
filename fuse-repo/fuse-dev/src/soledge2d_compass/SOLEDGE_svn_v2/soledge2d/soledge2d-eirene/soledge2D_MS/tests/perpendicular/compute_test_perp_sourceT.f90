subroutine compute_test_perp_sourceT(k)
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
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  B0=test_B0
  T0=test_T0ev*ev
  mass=zones(k)%species(1)%element%mass*m_u
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  n0=test_n0
  Dp=transport_parameters%Dn_p(1)
  nup=transport_parameters%nu_p(1)
  chip=transport_parameters%chii_p(1)
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))

        !ions 
        s2 = 1/r
        s4 = 1/(R0+r*cos(theta))
        s6 = 0.2D0*(R0+r*cos(theta))*(5.D0/2.D0*T0+0.25D0*T0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))+mass*G0**2*cos(theta)**2*sin(2*pi*(r-&
             rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))**2/2)*Dp*n0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*&
             pi/(rmax-rmin)+0.2D0*r*cos(theta)**2*(5.D0/2.D0*T0+0.25D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))+mass*G0**2*cos(theta)**2*sin(&
             2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin)))**2/2)*Dp*n0*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/&
             (rmax-rmin)
        s7 = s6
        s9 = 0.2D0*r*(R0+r*cos(theta))*(0.5D0*T0*cos(theta)*cos(2*pi*(r-rmin)&
             /(rmax-rmin))*pi/(rmax-rmin)+2*mass*G0**2*cos(theta)**2*sin(2*pi*&
             (r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))**2*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)-0.2D0&
             *mass*G0**2*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+&
             0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*n0*cos(2*pi&
             *(r-rmin)/(rmax-rmin))*pi/(rmax-rmin))*Dp*n0*cos(theta)*cos(2*pi*(&
             r-rmin)/(rmax-rmin))*pi/(rmax-rmin)
        s10 = -0.4D0*r*(R0+r*cos(theta))*(5.D0/2.D0*T0+0.25D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))+mass*G0**2*cos(theta)**2*sin(2*pi&
             *(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin)))**2/2)*Dp*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             *pi**2/(rmax-rmin)**2
        s8 = s9+s10
        s5 = s7+s8
        s3 = s4*s5
        s1 = s2*s3
        s3 = 1/r/(R0+r*cos(theta))*(0.2D0*(R0+r*cos(theta))*Chip*(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*cos(2&
             *pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.2D0*r*cos(theta)**2*Chip&
             *(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(2&
             *pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.4D-1*r*(R0+r*cos(theta)&
             )*Chip*n0*cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax&
             -rmin)**2*T0-0.4D0*r*(R0+r*cos(theta))*Chip*(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))*pi**2/(rmax-rmin)**2)
        s5 = 1/r
        s7 = 1/(R0+r*cos(theta))
        s10 = (R0+r*cos(theta))*mass*nup*(n0+0.1D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))*(4*G0**2*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax&
             -rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2&
             *cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)-0.4D0*G0**2*cos(&
             theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta&
             )*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*n0*cos(2*pi*(r-rmin)/(rmax-rmin))&
             *pi/(rmax-rmin))/2
        s11 = r*cos(theta)*mass*nup*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*(4*G0**2*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*cos(2*pi&
             *(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)-0.4D0*G0**2*cos(theta&
             )**3*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))**3*n0*cos(2*pi*(r-rmin)/(rmax-rmin))*&
             pi/(rmax-rmin))/2
        s9 = s10+s11
        s10 = s9
        s12 = 0.1D0*r*(R0+r*cos(theta))*mass*nup*n0*cos(theta)*cos(2*pi*(r&
             -rmin)/(rmax-rmin))*pi/(rmax-rmin)*(4*G0**2*cos(theta)**2*sin(2*pi&
             *(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))**2*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)-0.4D0&
             *G0**2*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*n0*cos(2*pi*(r-rmin)&
             /(rmax-rmin))*pi/(rmax-rmin))
        s14 = r*(R0+r*cos(theta))/2
        s16 = mass*nup
        s18 = n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s20 = 8*G0**2*cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2&
             /(rmax-rmin)**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))**2-0.32D1*G0**2*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))&
             /(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*cos(2&
             *pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2*n0
        s19 = s20-8*G0**2*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2/&
             (n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*pi**2/(&
             rmax-rmin)**2+0.24D0*G0**2*cos(theta)**4*sin(2*pi*(r-rmin)/(rmax-rmin)&
             )**2/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**4&
             *n0**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2+0.8D0&
             *G0**2*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))**3/(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*n0*pi**2/(rmax-rmin)**2
        s17 = s18*s19
        s15 = s16*s17
        s13 = s14*s15
        s11 = s12+s13
        s8 = s10+s11
        s6 = s7*s8
        s4 = s5*s6
        s2 = s3+s4
        test_sources(k,1)%ST(i,j)=s1+s2
        test_sources(k,1)%ST(i,j)=-test_sources(k,1)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)

        !electrons
        s1 = 1/r/(R0+r*cos(theta))*(0.5D0*(R0+r*cos(theta))*(T0+0.1D0*T0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*Dp*n0*cos(theta)*cos(2*pi&
             *(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.5D0*r*cos(theta)**2*(T0+0.1D0&
             *T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*Dp*n0*cos(2*pi*&
             (r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.1D0*r*(R0+r*cos(theta))*T0*&
             cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2*Dp&
             *n0-0.1D1*r*(R0+r*cos(theta))*(T0+0.1D0*T0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*Dp*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin))*pi**2/(rmax-rmin)**2)
        s2 = 1/r/(R0+r*cos(theta))*(0.2D0*(R0+r*cos(theta))*Chip*(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*cos(2&
             *pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.2D0*r*cos(theta)**2*Chip&
             *(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(2&
             *pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.4D-1*r*(R0+r*cos(theta)&
             )*Chip*n0*cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax&
             -rmin)**2*T0-0.4D0*r*(R0+r*cos(theta))*Chip*(n0+0.1D0*n0*cos(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))*pi**2/(rmax-rmin)**2)
        test_sources(k,0)%ST(i,j)=s1+s2
        test_sources(k,0)%ST(i,j)=-test_sources(k,0)%ST(i,j)&
             /(kb*reference_parameters%fields%n0*reference_parameters%fields%T0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_perp_sourceT
