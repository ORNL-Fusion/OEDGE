subroutine compute_source_para_passive_scalar(zone,SP,Fluxes_SP,Nx,Nz)
  use all_variables, only : reference_parameters, global_variables
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  real*8,intent(in) :: SP(0:Nx+1,0:Nz+1)
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(out) :: Fluxes_SP(1:Nx,1:Nz,1:4)
  real*8, allocatable :: Flu(:,:)
  real*8, allocatable :: density(:,:),Gamma(:,:)
  real*8, allocatable :: SPred(:,:)
  real*8, allocatable :: densitym(:,:),densityp(:,:)
  real*8, allocatable :: Gammam(:,:),Gammap(:,:)
  real*8, allocatable :: SPm(:,:),SPp(:,:)
  real*8, allocatable :: Fleft(:,:),Fright(:,:)
  real*8, allocatable :: lambdam(:,:),lambdap(:,:)
  real*8, allocatable :: fluxm(:,:),fluxp(:,:)
  real*8, allocatable :: Varm(:,:),Varp(:,:)
  real*8, allocatable :: Jm(:,:),Jp(:,:),Gm(:,:),Gp(:,:)
  real*8, allocatable :: Jmed(:,:),Gmed(:,:)
  real*8,allocatable :: vpolm(:,:),vpolp(:,:)
  real*8 :: gam
  real*8 :: norm_ethetae, norm_ethetaw
  integer*4 :: i,j,l
  real*8 :: R0,rs0,c0

  allocate(Flu(1:Nx,0:Nz))
  allocate(density(1:Nx,0:Nz+1),Gamma(1:Nx,0:Nz+1),SPred(1:Nx,0:Nz+1))
  allocate(densitym(1:Nx,0:Nz),densityp(1:Nx,0:Nz))
  allocate(Gammam(1:Nx,0:Nz),Gammap(1:Nx,0:Nz))
  allocate(SPm(1:Nx,0:Nz),SPp(1:Nx,0:Nz))
  allocate(Fleft(1:Nx,0:Nz),Fright(1:Nx,0:Nz))
  allocate(lambdam(1:Nx,0:Nz),lambdap(1:Nx,0:Nz))
  allocate(fluxm(1:Nx,0:Nz),fluxp(1:Nx,0:Nz))
  allocate(Varm(1:Nx,0:Nz),Varp(1:Nx,0:Nz))
  allocate(Jm(1:Nx,0:Nz),Jp(1:Nx,0:Nz),Gm(1:Nx,0:Nz),Gp(1:Nx,0:Nz))
  allocate(Jmed(1:Nx,0:Nz),Gmed(1:Nx,0:Nz))
  allocate(vpolm(1:Nx,0:Nz),vpolp(1:Nx,0:Nz))

  gam=5.d0/3.d0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  density=zone%species(0)%var(1)%density(1:Nx,0:Nz+1)
  Gamma=zone%species(0)%var(1)%Gamma(1:Nx,0:Nz+1)
  SPred=SP(1:Nx,0:Nz+1)

  densitym=density(:,0:Nz)
  densityp=density(:,1:Nz+1)
  Gammam=Gamma(:,0:Nz)
  Gammap=Gamma(:,1:Nz+1)
  SPm=SPred(:,0:Nz)
  SPp=SPred(:,1:Nz+1)

  vpolm=(zone%species(1)%drifts%uEt(1:Nx,0:Nz)&
       +zone%species(1)%drifts%uBt(1:Nx,0:Nz))
  vpolp=(zone%species(1)%drifts%uEt(1:Nx,1:Nz+1)&
       +zone%species(1)%drifts%uBt(1:Nx,1:Nz+1))

  vpolm=(vpolm+vpolp)*0.5D0
  vpolp=vpolm

  Gmed=(zone%metric_coefficients%G(1:Nx,0:Nz)&
       +zone%metric_coefficients%G(1:Nx,1:Nz+1))*0.5d0
  Jmed=(zone%metric_coefficients%jacobian(1:Nx,0:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx,1:Nz+1))*0.5d0

  call eigenvalues_SP(Nx,Nz,densitym,Gammam,lambdam,Gmed,vpolm) !same as for Te
  call eigenvalues_SP(Nx,Nz,densityp,Gammap,lambdap,Gmed,vpolp) !same as for Te
  !computation of the raw flux
  call express_passive_scalar(Nx,Nz,densitym,Gammam,SPm,fluxm,Gmed,Jmed,vpolm)
  call express_passive_scalar(Nx,Nz,densityp,Gammap,SPp,fluxp,Gmed,Jmed,vpolp)
  do i=1,Nx
     do j=0,Nz
        Varm(i,j)=SPm(i,j)*Jmed(i,j)
        Varp(i,j)=SPp(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_e(Nx,Nz,lambdam,lambdap,fluxm,fluxp,Varm,Varp,Fleft,Fright) !same as for Te

  DO i=1,Nx
     do j=0,Nz
        FLU(i,j)=(Fleft(i,j)+Fright(i,j))/Jmed(i,j)
     end do
  ENDDO

  Fluxes_SP=0.D0
  do i=1,Nx
     do j=1,Nz
        norm_ethetae=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j+1)))
        norm_ethetaw=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j-1)))
        Fluxes_SP(i,j,3)=FLU(i,j)/norm_ethetae
        Fluxes_SP(i,j,4)=FLU(i,j-1)/norm_ethetaw
     end do
  end do

end subroutine compute_source_para_passive_scalar
