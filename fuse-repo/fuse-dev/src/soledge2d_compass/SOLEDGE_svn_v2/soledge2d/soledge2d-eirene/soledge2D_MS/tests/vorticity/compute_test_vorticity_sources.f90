subroutine compute_test_vorticity_sources(k)
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,T0,n0
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16
  real*8 :: me, mi
  Nx=zones(k)%mesh%Nx
  Nz=zones(k)%mesh%Nz
  psi0=test_psi0
  R0=test_R0
  G0=test_Gamma0
  n0=test_n0
  B0=test_B0
  T0=test_T0eV
  rmax=test_rmax
  rmin=test_rmin
  a=test_a
  me=m_e
  mi=m_u*global_parameters%element_list(1)%mass
  do i=1,Nx
     do j=1,Nz
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
        s1 = -psi0/a
        s3 = 1/(R0+r*cos(theta))
        s5 = 1/r
        s10 = 0.625D2*psi0/a/(R0+r*cos(theta))**2/sqrt(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0&
             **2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*sin(&
             theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))
        s11 = -0.3125D2*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*&
             cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)**3*log(4*pi*me/mi)*T0*&
             sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(2*B0**2*R0**2/(R0+r*cos(&
             theta))**3*r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(&
             theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(&
             theta))+0.625D2*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0&
             +r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+&
             psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s9 = s10+s11
        s10 = s9-0.125D4*psi0/a/(R0+r*cos(theta))**2/sqrt(B0**2*R0**2/(R0+&
             r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))&
             /(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*sin(theta)
        s11 = s10+0.625D3*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+&
             r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)**3*(-0.1D0*n0*sin(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(&
             2*pi*(r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))/(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(2*B0**2&
             *R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(&
             theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(&
             theta)**2*r*sin(theta))
        s12 = s11
        s14 = -0.125D4*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))+0.2D-1*n0*sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-&
             rmin))**2*T0-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin)))*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/(n0+0.1D0&
             *n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s15 = -0.125D3*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/(n0+&
             0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))**2*n0*sin(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin))
        s13 = s14+s15
        s8 = s12+s13
        s9 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             *sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))**3
        s7 = s8*s9
        s11 = -0.3125D2*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*&
             cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0&
             **2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*sin(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s12 = 0.625D3*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*sin(theta)*&
             sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/(n0&
             +0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s10 = s11+s12
        s11 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             **3*sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-&
             rmin)))**3*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2&
             /a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+r*&
             cos(theta))**3*cos(theta)**2*r*sin(theta))
        s9 = s10*s11
        s12 = -0.1875D3
        s15 = 0.5D-1*psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*sin(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s16 = -psi0/a/(R0+r*cos(theta))/r/sqrt(B0**2*R0**2/(R0+r*cos(theta&
             ))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*cos(theta)**2)*(-0.1D0*n0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/&
             (rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/(n0+0.1D0*&
             n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))
        s14 = s15+s16
        s15 = 1/sqrt(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(&
             theta)**2)
        s13 = s14*s15
        s11 = s12*s13
        s12 = sqrt(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*&
             T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s10 = s11*s12
        s8 = s9+s10
        s6 = s7+s8
        s4 = s5*s6
        s2 = s3*s4
        test_sources(k,1)%SW(i,j)=s1*s2
        test_sources(k,1)%SW(i,j)=-test_sources(k,1)%SW(i,j)/(reference_parameters%fields%W0&
             /reference_parameters%fields%tau0)
     end do
  end do
end subroutine compute_test_vorticity_sources
