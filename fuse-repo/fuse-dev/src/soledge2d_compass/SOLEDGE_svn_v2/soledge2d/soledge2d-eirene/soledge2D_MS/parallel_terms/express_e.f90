subroutine express_e(Nx,Nz,density,Gamma,T,Flux,G,Jac,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,T,G,Jac,vpol
  real*8,dimension(1:Nx,0:Nz),intent(out) :: Flux
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  5.D0/2.D0*Gamma(i,j)*T(i,j)*G(i,j)+&
             3.d0/2.d0*density(i,j)*T(i,j)*vpol(i,j)
     end do
  ENDDO
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  FLUX(i,j)*Jac(i,j)
     end do
  end DO
end subroutine express_e
