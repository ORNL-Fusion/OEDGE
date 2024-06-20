subroutine compute_test_perp_sourceG(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters,transport_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,Dp,n0,nup
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13
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
  nup=transport_parameters%nu_p(1)
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        s2 = 1/r
        s4 = 1/(R0+r*cos(theta))
        s6 = 0.2D0*(R0+r*cos(theta))*Dp*G0*cos(theta)**2*sin(2*pi*(r-rmin)&
             /(rmax-rmin))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))*n0*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)+0.2D0*r*cos(theta)**3&
             *Dp*G0*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*n0*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)
        s5 = s6+0.4D0*r*(R0+r*cos(theta))*Dp*G0*cos(theta)**2*cos(2*pi*(r-&
             rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*n0-0.4D-1*r*(R0+r*cos(theta))*Dp*&
             G0*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0**2*cos(2*pi*(r-rmin)/(&
             rmax-rmin))**2*pi**2/(rmax-rmin)**2-0.4D0*r*(R0+r*cos(theta))*Dp*G0&
             *cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*n0*pi**2/(rmax-rmin)**2
        s3 = s4*s5
        s1 = s2*s3
        s3 = 1/r
        s5 = 1/(R0+r*cos(theta))
        s8 = (R0+r*cos(theta))*nup*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*(2*G0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*pi&
             /(rmax-rmin)/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))-0.2D0*G0*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*cos(2*pi*(r-&
             rmin)/(rmax-rmin))*pi/(rmax-rmin))
        s9 = r*cos(theta)*nup*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*(2*G0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-&
             rmin)/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-0.2D0&
             *G0*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0*n0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*cos(2*pi*(r-rmin)&
             /(rmax-rmin))*pi/(rmax-rmin))
        s7 = s8+s9
        s8 = s7
        s10 = 0.2D0*r*(R0+r*cos(theta))*nup*n0*cos(theta)*cos(2*pi*(r-rmin&
             )/(rmax-rmin))*pi/(rmax-rmin)*(2*G0*cos(theta)*cos(2*pi*(r-rmin)/(&
             rmax-rmin))*pi/(rmax-rmin)/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))-0.2D0*G0*cos(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))&
             /(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0&
             *cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin))
        s12 = r*(R0+r*cos(theta))
        s13 = nup*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*&
             (-4*G0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*pi**2/(rmax-rmin)**2&
             /(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-0.8D0*&
             G0*cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2&
             /(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0+&
             0.8D-1*G0*cos(theta)**3*sin(2*pi*(r-rmin)/(rmax-rmin))/(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**3*n0**2*cos(2*pi*(&
             r-rmin)/(rmax-rmin))**2*pi**2/(rmax-rmin)**2+0.4D0*G0*cos(theta)**2&
             *sin(2*pi*(r-rmin)/(rmax-rmin))**2/(n0+0.1D0*n0*cos(theta)*sin(2*&
             pi*(r-rmin)/(rmax-rmin)))**2*n0*pi**2/(rmax-rmin)**2)
        s11 = s12*s13
        s9 = s10+s11
        s6 = s8+s9
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,1)%SG(i,j)=s1+s2
        test_sources(k,1)%SG(i,j)=-test_sources(k,1)%SG(i,j)&
             /(reference_parameters%fields%n0*reference_parameters%fields%c0/reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_perp_sourceG
