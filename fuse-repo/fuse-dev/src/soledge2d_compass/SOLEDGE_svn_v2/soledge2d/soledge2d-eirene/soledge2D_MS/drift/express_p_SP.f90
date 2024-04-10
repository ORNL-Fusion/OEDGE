subroutine express_p_SP(Nx,Nz,SP,Flux,Jac,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx,1:Nz),intent(in)::SP,Jac,vpol
  real*8,dimension(0:Nx,1:Nz),intent(out) :: Flux
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=0,Nx
     do j=1,Nz
        FLUX(i,j) = SP(i,j)*vpol(i,j)
     end do
  ENDDO
  DO i=0,Nx
     do j=1,Nz
        FLUX(i,j) =  FLUX(i,j)*Jac(i,j)
     end do
  end DO
end subroutine express_p_SP
