subroutine leigvec(Nx,Nz,density,Gamma,T,lev,mass)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(1:Nx,0:Nz),intent(in)::density,Gamma,T
  real*8,dimension(1:Nx,0:Nz,1:3,1:3),intent(out)::lev
  real*8 :: mass
  integer*4 i,j,k
  real*8 u,cs,gam
  real*8 norm,eps
  eps=1.d-10
  gam=5.D0/3.D0
  do i=1,Nx
     do j=0,Nz
        cs =SQRT(gam*T(i,j)/mass)
        u=Gamma(i,j)/(density(i,j)+eps)
        lev(i,j,1,1) = 0.5d0*u*u-1.5d0*cs*cs 
        lev(i,j,1,2) =  -u
        lev(i,j,1,3) = 1.D0/mass
        lev(i,j,2,1) = 0.5D0*u*u+1.5D0*u*cs
        lev(i,j,2,2) =  -u-1.5D0*cs
        lev(i,j,2,3) =  1.D0/mass
        lev(i,j,3,1) = 0.5D0*u*u-1.5D0*u*cs
        lev(i,j,3,2) =  -u+1.5D0*cs
        lev(i,j,3,3) = 1.D0/mass
     end do
  end do
  do k=1,3
     do i=1,Nx
        do j=0,Nz
           cs =SQRT(gam*T(i,j)/mass)
           lev(i,j,1,k)=lev(i,j,1,k)*(-2.D0/(3.D0*cs*cs))
           lev(i,j,2,k)=lev(i,j,2,k)/(3.D0*cs*cs)
           lev(i,j,3,k)=lev(i,j,3,k)/(3.D0*cs*cs)
        end do
     end do
  end do
end subroutine leigvec
