subroutine eigenvalues(Nx,Nz,density,Gamma,T,lambda,mass,Gmed,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,T,Gmed,vpol
  real*8,dimension(1:Nx,0:Nz,1:3),intent(out)::lambda
  real*8,intent(in) :: mass
  integer*4 i,j
  real*8 u,cs
  real*8 gam
  gam=5.D0/3.D0
  do i=1,Nx
     do j=0,Nz
        u=Gamma(i,j)/density(i,j)
        cs=sqrt(gam*T(i,j)/mass)
        lambda(i,j,1) =  u*Gmed(i,j)+vpol(i,j)
        lambda(i,j,2) =  (u-cs)*Gmed(i,j)+vpol(i,j)
        lambda(i,j,3) =  (u+cs)*Gmed(i,j)+vpol(i,j)
     end do
  end do
end subroutine eigenvalues
