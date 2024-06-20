subroutine invert_psi_tridiag(flux_surface)
  use MFlux_surface
  implicit none
  Type(Tflux_surface) :: flux_surface
  real*8 alpha,beta,gamma,fact
  integer*4 j,INFO
  real*8,allocatable:: DL(:),D(:),DU(:),B(:,:),x(:),z(:)
  integer*4 :: Nz
  external DGTSV
  Nz=flux_surface%Nz
  allocate(DL(1:Nz-1))
  allocate(D(1:Nz))
  allocate(DU(1:Nz-1))
  allocate(B(1:Nz,1))
  allocate(x(1:Nz),z(1:Nz))
  if(.not.flux_surface%properties%is_periodic) then
     do j=1,Nz-1
        DL(j)=flux_surface%tridiag%a(j+1)
     end do
     do j=1,Nz
        D(j)=flux_surface%tridiag%b(j)
     end do
     do j=1,Nz-1
        DU(j)=flux_surface%tridiag%c(j)
     end do
     do j=1,Nz
        B(j,1)=flux_surface%tridiag%S(j)
     end do
     call dgtsv(Nz,1, DL, D, DU, B, Nz, INFO )
     do j=1,Nz
        flux_surface%tridiag%buffer(j)=B(j,1)
     end do
  else
     alpha=flux_surface%tridiag%c(Nz)
     beta=flux_surface%tridiag%a(1)
     gamma=-flux_surface%tridiag%b(1)
     do j=1,Nz-1
        DL(j)=flux_surface%tridiag%a(j+1)
     end do
     do j=1,Nz
        D(j)=flux_surface%tridiag%b(j)
     end do
     D(1)=D(1)-gamma
     D(Nz)=D(Nz)-alpha*beta/gamma
     do j=1,Nz-1
        DU(j)=flux_surface%tridiag%c(j)
     end do
     do j=1,Nz
        B(j,1)=flux_surface%tridiag%S(j)
     end do
     call DGTSV(Nz,1, DL, D, DU, B, Nz, INFO )
     do j=1,Nz
        x(j)=B(j,1)
     end do
     do j=1,Nz-1
        DL(j)=flux_surface%tridiag%a(j+1)
     end do
     do j=1,Nz
        D(j)=flux_surface%tridiag%b(j)
     end do
     D(1)=D(1)-gamma
     D(Nz)=D(Nz)-alpha*beta/gamma
     do j=1,Nz-1
        DU(j)=flux_surface%tridiag%c(j)
     end do
     B(1:Nz,1)=0.d0
     B(1,1)=gamma
     B(Nz,1)=alpha
     call DGTSV(Nz,1, DL, D, DU, B, Nz, INFO )
     do j=1,Nz
        z(j)=B(j,1)
     end do
     fact=(x(1)+beta*x(Nz)/gamma)/&
          (1.d0+z(1)+beta*z(Nz)/gamma)
     do j=1,Nz
        flux_surface%tridiag%buffer(j)=x(j)-fact*z(j)
     end do
  end if
  deallocate(D,DU,DL,B,x,z)
end subroutine invert_psi_tridiag
