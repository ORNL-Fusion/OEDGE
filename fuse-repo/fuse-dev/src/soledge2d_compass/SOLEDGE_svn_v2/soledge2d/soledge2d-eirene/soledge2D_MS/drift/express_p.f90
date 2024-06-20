subroutine express_p(Nx,Nz,density,Gamma,T,Flux,Jac,mass,vpol)
  implicit none
  integer*4 Nx,Nz
  real*8,dimension(0:Nx,1:Nz),intent(in)::density,Gamma,T,Jac,vpol
  real*8,dimension(0:Nx,1:Nz,1:3),intent(out) :: Flux
  real*8 :: mass
  integer*4 i,j,n
  real*8  eps
  eps=1.d-10
  DO i=0,Nx
     do j=1,Nz
        FLUX(i,j,1) = density(i,j)*vpol(i,j)
        FLUX(i,j,2) = Gamma(i,j)*vpol(i,j)
        FLUX(i,j,3) = (3.d0/2.d0*density(i,j)*T(i,j)+0.5d0*mass*Gamma(i,j)**2.d0/density(i,j))*vpol(i,j)
     end do
  ENDDO
  DO i=0,Nx
     do j=1,Nz
        do n=1,3
           FLUX(i,j,n) =  FLUX(i,j,n)*Jac(i,j)
        end do
     end do
  end DO
end subroutine express_p
