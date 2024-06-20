subroutine eigenvalues_e(Nx,Nz,density,Gamma,lambda,Gmed,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in):: density,Gamma,Gmed,vpol
  real*8,dimension(1:Nx,0:Nz),intent(out)::lambda
  integer*4 i,j
  real*8 u,cs
  real*8 gam
  gam=5./3.
  do i=1,Nx
     do j=0,Nz
        lambda(i,j) =  gam*((Gamma(i,j)/density(i,j)*Gmed(i,j))+vpol(i,j))
     end do
  end do
end subroutine eigenvalues_e
