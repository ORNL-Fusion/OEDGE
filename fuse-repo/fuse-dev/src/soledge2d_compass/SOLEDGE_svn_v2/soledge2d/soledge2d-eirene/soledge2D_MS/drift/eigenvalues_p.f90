subroutine eigenvalues_p(Nx,Nz,lambda,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx,1:Nz),intent(in)::vpol
  real*8,dimension(0:Nx,1:Nz,1:3),intent(out)::lambda
  integer*4 i,j
  do i=0,Nx
     do j=1,Nz
        lambda(i,j,1) =  vpol(i,j)
        lambda(i,j,2) =  vpol(i,j)
        lambda(i,j,3) =  vpol(i,j)
     end do
  end do
end subroutine eigenvalues_p
