subroutine compute_analytic_vorticity(k)
#include "compile_opt.inc"
  use test_var
  use all_variables, only : global_parameters,zones,reference_parameters
  use Mphysics
  implicit none
  integer*4,intent(in) :: k
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: r,theta,test_theta,test_r
  real*8 :: psi0,R0,G0,B0,rmax,rmin,a,T0,n0
  real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17
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
  allocate(test_vort(k)%W(0:Nx+1,0:Nz+1))
  allocate(test_vort(k)%phi(0:Nx+1,0:Nz+1))
  allocate(test_vort(k)%r(0:Nx+1,0:Nz+1))
  allocate(test_vort(k)%theta(0:Nx+1,0:Nz+1))
  do i=0,Nx+1
     do j=0,Nz+1
        r=test_r(zones(k)%mesh%x(i,j))
        theta=test_theta(zones(k)%mesh%z(i,j),zones(k)%mesh%x(i,j))
#if VORTICITY_PI == 0
        s1 = 1/r
        s3 = 1/(R0+r*cos(theta))
        s7 = -0.5D-2*n0*sin(theta)**2*sin(2*pi*(r-rmin)/(rmax-rmin))**2*r*&
             (R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(&
             R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2/(B0**2&
             *R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*&
             sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2))*log(4*pi*me/mi)*T0
        s8 = -0.5D-1*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))&
             *r**2*sin(theta)**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2&
             /(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             )*log(4*pi*me/mi)*T0*sin(2*pi*(r-rmin)/(rmax-rmin))
        s6 = s7+s8
        s7 = s6-0.5D-1*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)**2*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2&
             /r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             )*log(4*pi*me/mi)*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*&
             (2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2&
             /(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta)&
             )**3*cos(theta)**2*r*sin(theta))
        s8 = s7
        s11 = 0.5D-1*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)
        s12 = (-2*psi0**2/a**2/(R0+r*cos(theta))**3/r/(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*sin(theta)+psi0**2/a**2&
             /(R0+r*cos(theta))**2/r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*cos(theta)**2)**2*(2*B0**2*R0**2/(R0+r*cos(theta))**3&
             *r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+&
             2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(theta)))*log(&
             4*pi*me/mi)*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s10 = s11*s12
        s11 = 0.5D-1*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2&
             /(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2&
             ))*log(4*pi*me/mi)*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s9 = s10+s11
        s5 = s8+s9
        s6 = s5-0.2D-1*n0*cos(theta)**2*cos(2*pi*(r-rmin)/(rmax-rmin))**2*&
             pi**2/(rmax-rmin)**2*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta)&
             )**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0-0.1D0*(&
             n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*(R0+r*cos(theta)&
             )/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta&
             )**2)*log(4*pi*me/mi)*T0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))&
             *pi/(rmax-rmin)
        s7 = s6-0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             )))*r*cos(theta)**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))&
             **2*cos(theta)**2)*log(4*pi*me/mi)*T0*cos(2*pi*(r-rmin)/(rmax-rmin&
             ))*pi/(rmax-rmin)
        s8 = s7
        s10 = 0.1D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)**2*log(4*pi*me/mi)*T0*cos(theta)*cos(2*pi*(r-&
             rmin)/(rmax-rmin))*pi/(rmax-rmin)*(-2*B0**2*R0**2/(R0+r*cos(theta)&
             )**3*cos(theta)-2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**2*&
             cos(theta)-2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**3)
        s11 = 0.2D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))&
             )*r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))*pi**2/(rmax-rmin)**2
        s9 = s10+s11
        s4 = s8+s9
        s2 = s3*s4
        test_vort(k)%W(i,j) = s1*s2*global_parameters%element_list(1)%mass*m_u
#endif
#if VORTICITY_PI == 1
        s2 = 1/r
        s4 = 1/(R0+r*cos(theta))
        s8 = -0.5D-1*r**2*sin(theta)**2*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2&
             /(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             )*log(4*pi*me/mi)*T0*sin(2*pi*(r-rmin)/(rmax-rmin))
        s9 = -0.5D-2*r*(R0+r*cos(theta))*n0*sin(theta)**2*sin(2*pi*(r-rmin&
             )/(rmax-rmin))**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(&
             R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2&
             *cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2/(B0**2&
             *R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*&
             sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2))*log(&
             4*pi*me/mi)*T0
        s7 = s8+s9
        s8 = s7-0.5D-1*r*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)**2*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2&
             /r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)&
             **2))*log(4*pi*me/mi)*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin&
             ))*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(theta)+2*psi0**2/a**2&
             /(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0**2/a**2/(R0+r*cos(theta)&
             )**3*cos(theta)**2*r*sin(theta))
        s9 = s8
        s12 = 0.5D-1*r*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)
        s13 = (-2*psi0**2/a**2/(R0+r*cos(theta))**3/r/(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*sin(theta)+psi0**2/a**2&
             /(R0+r*cos(theta))**2/r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*cos(theta)**2)**2*(2*B0**2*R0**2/(R0+r*cos(theta))**3&
             *r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+&
             2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(theta)))*log(&
             4*pi*me/mi)*T0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s11 = s12*s13
        s12 = 0.5D-1*r*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*&
             (r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2&
             /(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2&
             ))*log(4*pi*me/mi)*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))
        s10 = s11+s12
        s6 = s9+s10
        s7 = s6-0.1D0*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*cos(theta)*cos(2*pi*(r-rmin)&
             /(rmax-rmin))*pi/(rmax-rmin)-0.1D0*r*cos(theta)**2*(n0+0.1D0*n0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*cos(2&
             *pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin)
        s8 = s7-0.2D-1*r*(R0+r*cos(theta))*n0*cos(theta)**2*cos(2*pi*(r-rmin)&
             /(rmax-rmin))**2*pi**2/(rmax-rmin)**2/(B0**2*R0**2/(R0+r*cos(theta)&
             )**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*cos(theta)**2)*log(4*pi*me/mi)*T0
        s9 = s8
        s11 = 0.1D0*r*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)**2*log(4*pi*me/mi)*T0*cos(theta)*cos(2*pi*(r-&
             rmin)/(rmax-rmin))*pi/(rmax-rmin)*(-2*B0**2*R0**2/(R0+r*cos(theta)&
             )**3*cos(theta)-2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**2*&
             cos(theta)-2*psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**3)
        s12 = 0.2D0*r*(R0+r*cos(theta))*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)*log(4*pi*me/mi)*T0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))*pi**2/(rmax-rmin)**2
        s10 = s11+s12
        s5 = s9+s10
        s3 = s4*s5
        s1 = s2*s3
        s3 = 1/r
        s5 = 1/(R0+r*cos(theta))
        s9 = -r**2*sin(theta)/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta&
             ))**2*cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**&
             2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))&
             **2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)&
             )*(-0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*T0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*sin(theta)*sin(2*pi*(r&
             -rmin)/(rmax-rmin)))
        s11 = -r*(R0+r*cos(theta))
        s13 = 1/((B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta&
             )**2)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2/r**2/(B0**2*&
             R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)&
             **2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2))
        s14 = (-0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0&
             *T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0*&
             cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin)))*(2*B0**2*R0**2/(R0+r*cos(theta))**3*r*sin(&
             theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2*psi0&
             **2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s12 = s13*s14
        s10 = s11*s12
        s8 = s9+s10
        s9 = s8
        s12 = r*(R0+r*cos(theta))
        s14 = 1/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)**2)
        s16 = -2*psi0**2/a**2/(R0+r*cos(theta))**3/r/(B0**2*R0**2/(R0+r*cos(&
             theta))**2+psi0**2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*cos(theta)**2)*sin(theta)+psi0**2/a**&
             2/(R0+r*cos(theta))**2/r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0&
             **2/a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*cos(theta)**2)**2*(2*B0**2*R0**2/(R0+r*cos(theta))**3*&
             r*sin(theta)+2*psi0**2/a**2/(R0+r*cos(theta))**3*sin(theta)**3*r+2&
             *psi0**2/a**2/(R0+r*cos(theta))**3*cos(theta)**2*r*sin(theta))
        s17 = -0.1D0*n0*sin(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0*&
             T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))-0.1D0*(n0+0.1D0*n0&
             *cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*sin(theta)*sin(2*pi&
             *(r-rmin)/(rmax-rmin))
        s15 = s16*s17
        s13 = s14*s15
        s11 = s12*s13
        s12 = r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/&
             a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(1/r**2-psi0**2/a**2/(R0+r*cos(theta))**2&
             /r**2/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta))**2*cos(theta)&
             **2))*(-0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin))*(T0+0.1D0&
             *T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))+0.2D-1*n0*sin(theta)&
             **2*sin(2*pi*(r-rmin)/(rmax-rmin))**2*T0-0.1D0*(n0+0.1D0*n0*cos(&
             theta)*sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin)))
        s10 = s11+s12
        s7 = s9+s10
        s8 = s7+(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(0.2D0*n0*cos(theta)*cos(2*pi*(r-rmin)/(rmax&
             -rmin))*pi/(rmax-rmin)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin)))+0.2D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))*T0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin))
        s9 = s8+r*cos(theta)/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2/a**2&
             /(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(0.2D0*n0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin)&
             )*pi/(rmax-rmin)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-rmin)/(rmax&
             -rmin)))+0.2D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(rmax-rmin)&
             ))*T0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax-rmin))
        s10 = s9
        s12 = -r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(&
             theta))**2*cos(theta)**2)**2*(0.2D0*n0*cos(theta)*cos(2*pi*(r-rmin&
             )/(rmax-rmin))*pi/(rmax-rmin)*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin)))+0.2D0*(n0+0.1D0*n0*cos(theta)*sin(2*pi*(r-rmin&
             )/(rmax-rmin)))*T0*cos(theta)*cos(2*pi*(r-rmin)/(rmax-rmin))*pi/(rmax&
             -rmin))*(-2*B0**2*R0**2/(R0+r*cos(theta))**3*cos(theta)-2*psi0**2&
             /a**2/(R0+r*cos(theta))**3*sin(theta)**2*cos(theta)-2*psi0**2/a**2&
             /(R0+r*cos(theta))**3*cos(theta)**3)
        s13 = r*(R0+r*cos(theta))/(B0**2*R0**2/(R0+r*cos(theta))**2+psi0**2&
             /a**2/(R0+r*cos(theta))**2*sin(theta)**2+psi0**2/a**2/(R0+r*cos(theta)&
             )**2*cos(theta)**2)*(-0.4D0*n0*cos(theta)*sin(2*pi*(r-rmin)/(&
             rmax-rmin))*pi**2/(rmax-rmin)**2*(T0+0.1D0*T0*cos(theta)*sin(2*pi*(&
             r-rmin)/(rmax-rmin)))+0.8D-1*n0*cos(theta)**2*cos(2*pi*(r-rmin)/(&
             rmax-rmin))**2*pi**2/(rmax-rmin)**2*T0-0.4D0*(n0+0.1D0*n0*cos(theta)&
             *sin(2*pi*(r-rmin)/(rmax-rmin)))*T0*cos(theta)*sin(2*pi*(r-rmin)&
             /(rmax-rmin))*pi**2/(rmax-rmin)**2)
        s11 = s12+s13
        s6 = s10+s11
        s4 = s5*s6
        s2 = s3*s4
        test_vort(k)%W(i,j) = (s1+s2)*global_parameters%element_list(1)%mass*m_u
#endif
        test_vort(k)%phi(i,j) = test_T0eV*(1.D0+0.1D0*cos(theta)*sin(2*pi*(r-&
             rmin)/(rmax-rmin)))*(-0.5*log(4.D0*pi*me/mi))
        test_vort(k)%r(i,j) = r
        test_vort(k)%theta(i,j) = theta
     end do
  end do
end subroutine compute_analytic_vorticity
