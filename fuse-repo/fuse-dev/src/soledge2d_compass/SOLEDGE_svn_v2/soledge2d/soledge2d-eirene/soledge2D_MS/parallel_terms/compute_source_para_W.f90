subroutine compute_source_para_W(zone,Fluxes_W,Nx,Nz)
  use all_variables, only : reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(out) :: Fluxes_W(1:Nx,1:Nz,1:4)
  real*8, allocatable :: Flu(:,:)
  real*8, allocatable :: density(:,:),Gamma(:,:)
  real*8, allocatable :: W(:,:)
  real*8, allocatable :: densitym(:,:),densityp(:,:)
  real*8, allocatable :: Gammam(:,:),Gammap(:,:)
  real*8, allocatable :: Wm(:,:),Wp(:,:)
  real*8, allocatable :: Fleft(:,:),Fright(:,:)
  real*8, allocatable :: lambdam(:,:),lambdap(:,:)
  real*8, allocatable :: fluxm(:,:),fluxp(:,:)
  real*8, allocatable :: Varm(:,:),Varp(:,:)
  real*8, allocatable :: Jm(:,:),Jp(:,:),Gm(:,:),Gp(:,:)
  real*8, allocatable :: Jmed(:,:),Gmed(:,:)
  real*8 :: gam
  integer*4 :: i,j,l
  real*8 :: R0,rs0

  allocate(Flu(1:Nx,0:Nz))
  allocate(density(1:Nx,0:Nz+1),Gamma(1:Nx,0:Nz+1),W(1:Nx,0:Nz+1))
  allocate(densitym(1:Nx,0:Nz),densityp(1:Nx,0:Nz))
  allocate(Gammam(1:Nx,0:Nz),Gammap(1:Nx,0:Nz))
  allocate(Wm(1:Nx,0:Nz),Wp(1:Nx,0:Nz))
  allocate(Fleft(1:Nx,0:Nz),Fright(1:Nx,0:Nz))
  allocate(lambdam(1:Nx,0:Nz),lambdap(1:Nx,0:Nz))
  allocate(fluxm(1:Nx,0:Nz),fluxp(1:Nx,0:Nz))
  allocate(Varm(1:Nx,0:Nz),Varp(1:Nx,0:Nz))
  allocate(Jm(1:Nx,0:Nz),Jp(1:Nx,0:Nz),Gm(1:Nx,0:Nz),Gp(1:Nx,0:Nz))
  allocate(Jmed(1:Nx,0:Nz),Gmed(1:Nx,0:Nz))

  gam=5.d0/3.d0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  density=zone%species(0)%var(1)%density(1:Nx,0:Nz+1)
  Gamma=zone%species(0)%var(1)%Gamma(1:Nx,0:Nz+1)
  W=zone%electric_fields(1)%vorticity(1:Nx,0:Nz+1)

  densitym=density(:,0:Nz)
  densityp=density(:,1:Nz+1)
  Gammam=Gamma(:,0:Nz)
  Gammap=Gamma(:,1:Nz+1)
  Wm=W(:,0:Nz)
  Wp=W(:,1:Nz+1)

  Gmed=(zone%metric_coefficients%G(1:Nx,0:Nz)&
       +zone%metric_coefficients%G(1:Nx,1:Nz+1))*0.5d0
  Jmed=(zone%metric_coefficients%jacobian(1:Nx,0:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx,1:Nz+1))*0.5d0

  call eigenvalues_e(Nx,Nz,densitym,Gammam,lambdam) !same as for Te
  call eigenvalues_e(Nx,Nz,densityp,Gammap,lambdap) !same as for Te
  do i=1,Nx
     do j=0,Nz
        lambdam(i,j)=lambdam(i,j)*Gmed(i,j)
        lambdap(i,j)=lambdap(i,j)*Gmed(i,j)
     end do
  end do
  !computation of the raw flux
  call express_W(Nx,Nz,densitym,Gammam,Wm,fluxm,Gmed,Jmed)
  call express_W(Nx,Nz,densityp,Gammap,Wp,fluxp,Gmed,Jmed)
  do i=1,Nx
     do j=0,Nz
        Varm(i,j)=Wm(i,j)*Jmed(i,j)
        Varp(i,j)=Wp(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_e(Nx,Nz,lambdam,lambdap,fluxm,fluxp,Varm,Varp,Fleft,Fright) !same as for Te

  DO i=1,Nx
     do j=0,Nz
        FLU(i,j)=(Fleft(i,j)+Fright(i,j))/(Jmed(i,j)*Gmed(i,j))
     end do
  ENDDO

  Fluxes_W=0.D0
  do i=1,Nx
     do j=1,Nz
        Fluxes_W(i,j,3)=FLU(i,j)*zone%metric_coefficients%sinepitch_east(i,j)
        Fluxes_W(i,j,4)=FLU(i,j-1)*zone%metric_coefficients%sinepitch_west(i,j)
     end do
  end do

  Fluxes_W=Fluxes_W*(2.d0*pi*R0/rs0)


end subroutine compute_source_para_W
