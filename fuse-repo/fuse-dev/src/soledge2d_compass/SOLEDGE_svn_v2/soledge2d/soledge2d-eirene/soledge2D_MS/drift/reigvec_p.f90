subroutine reigvec_p(Nx,Nz,density,Gamma,rev,mass)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx,1:Nz)::density,Gamma
  real*8,dimension(0:Nx,1:Nz,1:3,1:3)::rev
  real*8,intent(in) :: mass
  integer*4 i,j
  real*8 u,gam
  real*8 eps
  eps=1.d-10
  gam=5.D0/3.D0
  do i=0,Nx
     do j=1,Nz
        u=Gamma(i,j)/(density(i,j)+eps)
        rev(i,j,1,1)=0.d0
        rev(i,j,2,1)=1.D0
        rev(i,j,3,1)=mass*u
        rev(i,j,1,2)=1.D0
        rev(i,j,2,2)=0.5D0*u
        rev(i,j,3,2)=0.D0
        rev(i,j,1,3)=0.D0
        rev(i,j,2,3)=0.D0
        rev(i,j,3,3)=1.D0
     end do
  end do
end subroutine reigvec_p
