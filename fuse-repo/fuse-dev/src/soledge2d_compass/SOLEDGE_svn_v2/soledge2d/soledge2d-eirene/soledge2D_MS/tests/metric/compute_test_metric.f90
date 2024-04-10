subroutine compute_test_metric()
  use test_var
  use all_variables, only : global_parameters, zones, reference_parameters
  use Mphysics
  implicit none
  integer*4 :: k
  integer*4 :: i,j,Nx,Nz
  real*8 :: R0,psi0,a,thetashift,B0
  real*8 :: s1,s2,s3,s4,s5
  real*8 :: psi,theta
  B0=test_B0
  R0=test_R0
  a=test_a
  psi0=test_psi0
  thetashift=test_theta_shift
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     allocate(test_metric(k)%cpp(1:Nx,1:Nz))
     allocate(test_metric(k)%cpt(1:Nx,1:Nz))
     allocate(test_metric(k)%ctt(1:Nx,1:Nz))
     allocate(test_metric(k)%Jac(1:Nx,1:Nz))
     allocate(test_metric(k)%G(1:Nx,1:Nz))
     do i=1,Nx
        do j=1,Nz
           psi=zones(k)%mesh%x(i,j)
           theta=zones(k)%mesh%z(i,j)*2.D0*pi
           !Jac
           s1 = R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)
           s4 = a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0&
                *psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift&
                *pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift&
                *pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta&
                +psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*&
                pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
           s5 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)&
                +a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift&
                *pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta&
                +psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift*&
                pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
           s3 = s4+s5
           s2 = 1/s3
           test_metric(k)%Jac(i,j) = s1*s2

           !cpp
           s1 = a**2*(0.6D0+0.8D0*psi)**2*cos(theta+psi*thetashift*pi)**2/(-a&
                *(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta&
                +psi*thetashift*pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)&
                *thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0&
                *psi)*cos(theta+psi*thetashift*pi))**2
           s2 = a**2*(0.6D0+0.8D0*psi)**2*sin(theta+psi*thetashift*pi)**2/(-a&
                *(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta&
                +psi*thetashift*pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)&
                *thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0&
                *psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*&
                psi)*cos(theta+psi*thetashift*pi))**2
           test_metric(k)%cpp(i,j) = s1+s2
           test_metric(k)%cpp(i,j)=test_metric(k)%cpp(i,j)*reference_parameters%geometry%rs0**2.D0

           !cpt
           s1 = -a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)/(-a*(0.6D0+&
                0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift&
                *pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)&
                -(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*&
                cos(theta+psi*thetashift*pi))**2*(0.8D0*a*sin(theta+psi*thetashift&
                *pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)
           s2 = a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)/(-a*(0.6D0+0.8D0&
                *psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift&
                *pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift&
                *pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)&
                *sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+&
                psi*thetashift*pi))**2*(0.8D0*a*cos(theta+psi*thetashift*&
                pi)-a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)
           test_metric(k)%cpt(i,j) = s1+s2
           test_metric(k)%cpt(i,j)=test_metric(k)%cpt(i,j)*reference_parameters%geometry%rs0**2.D0/(2.D0*pi)

           !ctt
           s1 = (0.8D0*a*sin(theta+psi*thetashift*pi)+a*(0.6D0+0.8D0*psi)*cos(theta&
                +psi*thetashift*pi)*thetashift*pi)**2/(-a*(0.6D0+0.8D0*psi)*&
                sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)&
                +a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(&
                0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta&
                +psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi&
                *thetashift*pi))**2
           s2 = (0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta&
                +psi*thetashift*pi)*thetashift*pi)**2/(-a*(0.6D0+0.8D0*psi)*&
                sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)&
                +a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)-(&
                0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0+0.8D0*psi)*sin(theta&
                +psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)*cos(theta+psi&
                *thetashift*pi))**2
           test_metric(k)%ctt(i,j) = s1+s2
           test_metric(k)%ctt(i,j)=test_metric(k)%ctt(i,j)*reference_parameters%geometry%rs0**2.D0/(2.D0*pi)**2.D0

           !G
           s2 = psi0/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))*&
                sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta+psi*thetashift*pi)&
                +a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi)*thetashift*pi)/(&
                -a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta&
                +psi*thetashift*pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift&
                *pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0&
                +0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0*psi)&
                *cos(theta+psi*thetashift*pi))
           s3 = psi0/a/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))*&
                cos(theta+psi*thetashift*pi)*(0.8D0*a*cos(theta+psi*thetashift*pi)&
                -a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)/(&
                -a*(0.6D0+0.8D0*psi)*sin(theta+psi*thetashift*pi)*(0.8D0*a*sin(theta&
                +psi*thetashift*pi)+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift&
                *pi)*thetashift*pi)-(0.8D0*a*cos(theta+psi*thetashift*pi)-a*(0.6D0&
                +0.8D0*psi)*sin(theta+psi*thetashift*pi)*thetashift*pi)*a*(0.6D0+0.8D0&
                *psi)*cos(theta+psi*thetashift*pi))
           s1 = s2+s3
           s2 = 1/sqrt(B0**2*R0**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift&
                *pi))**2+psi0**2/a**2/(R0+a*(0.6D0+0.8D0*psi)*cos(theta+psi*&
                thetashift*pi))**2*sin(theta+psi*thetashift*pi)**2+psi0**2/a**2/(R0+&
                a*(0.6D0+0.8D0*psi)*cos(theta+psi*thetashift*pi))**2*cos(theta+psi&
                *thetashift*pi)**2)
           test_metric(k)%G(i,j) = s1*s2
           test_metric(k)%G(i,j)=test_metric(k)%G(i,j)*reference_parameters%geometry%R0
        end do

100     format(512es15.7)
     end do
  end do
end subroutine compute_test_metric
