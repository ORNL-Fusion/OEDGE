subroutine express_passive_scalar(Nx,Nz,density,Gamma,SP,Flux,G,Jac,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,SP,G,Jac,vpol
  real*8,dimension(1:Nx,0:Nz),intent(out) :: Flux
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  (Gamma(i,j)/density(i,j)*G(i,j)+vpol(i,j))*SP(i,j)
     end do
  ENDDO
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  FLUX(i,j)*Jac(i,j)
     end do
  end DO
end subroutine express_passive_scalar
