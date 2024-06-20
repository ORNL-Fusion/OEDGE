subroutine reigvec(Nx,Nz,density,Gamma,T,rev,mass)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz)::density,Gamma,T
  real*8,dimension(1:Nx,0:Nz,1:3,1:3)::rev
  real*8 :: mass
  integer*4 i,j
  real*8 u,cs,gam
  real*8 eps
  eps=1.d-10
  gam=5.D0/3.D0
  do i=1,Nx
     do j=0,Nz
        cs=SQRT(gam*T(i,j)/mass)
        u=Gamma(i,j)/(density(i,j)+eps)
        rev(i,j,1,1)=1.d0
        rev(i,j,2,1)=u
        rev(i,j,3,1)=0.5d0*u*u*mass
        rev(i,j,1,2)=1.D0
        rev(i,j,2,2)=u-cs
        rev(i,j,3,2)=(0.5D0*u*u-u*cs+1.5D0*cs*cs)*mass
        rev(i,j,1,3)=1.D0
        rev(i,j,2,3)=u+cs
        rev(i,j,3,3)=(0.5D0*u*u+u*cs+1.5D0*cs*cs)*mass
     end do
  end do
end subroutine reigvec
