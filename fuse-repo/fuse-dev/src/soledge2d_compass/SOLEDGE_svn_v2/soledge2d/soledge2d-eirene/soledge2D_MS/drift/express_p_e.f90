subroutine express_p_e(Nx,Nz,density,T,Flux,Jac,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx,1:Nz),intent(in)::density,T,Jac,vpol
  real*8,dimension(0:Nx,1:Nz),intent(out) :: Flux
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=0,Nx
     do j=1,Nz
        FLUX(i,j) = 3.d0/2.d0*density(i,j)*T(i,j)*vpol(i,j)
     end do
  ENDDO
  DO i=0,Nx
     do j=1,Nz
        FLUX(i,j) =  FLUX(i,j)*Jac(i,j)
     end do
  end DO
end subroutine express_p_e
