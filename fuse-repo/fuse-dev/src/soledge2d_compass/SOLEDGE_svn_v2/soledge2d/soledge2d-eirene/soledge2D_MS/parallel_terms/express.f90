subroutine express(Nx,Nz,density,Gamma,T,Flux,G,Jac,mass,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,T,G,Jac,vpol
  real*8,dimension(1:Nx,0:Nz,1:3),intent(out) :: Flux
  real*8 :: mass
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j,1) =  Gamma(i,j)*G(i,j)+density(i,j)*vpol(i,j)
        FLUX(i,j,2) =  (Gamma(i,j)**2.d0/(Density(i,j))&
             +T(i,j)*Density(i,j)/mass)*G(i,j)+Gamma(i,j)*vpol(i,j)
        FLUX(i,j,3) =  (5.D0/2.D0*Gamma(i,j)*T(i,j)&
             +0.5D0*mass*Gamma(i,j)**3.D0/(Density(i,j))**2.D0)*G(i,j)&
             +(3.d0/2.d0*density(i,j)*T(i,j)+0.5d0*mass*Gamma(i,j)**2.d0/density(i,j))*vpol(i,j)
     end do
  ENDDO
  DO i=1,Nx
     do j=0,Nz
        do n=1,3
           FLUX(i,j,n) =  FLUX(i,j,n)*Jac(i,j)
        end do
     end do
  end DO
end subroutine express
