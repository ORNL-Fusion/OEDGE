subroutine express_W(Nx,Nz,density,Gamma,W,Flux,G,Jac)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,W,G,Jac
  real*8,dimension(1:Nx,0:Nz),intent(out) :: Flux
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  Gamma(i,j)*W(i,j)
     end do
  ENDDO
  DO i=1,Nx
     do j=0,Nz
        FLUX(i,j) =  FLUX(i,j)*Jac(i,j)*G(i,j)
     end do
  end DO
end subroutine express_W
