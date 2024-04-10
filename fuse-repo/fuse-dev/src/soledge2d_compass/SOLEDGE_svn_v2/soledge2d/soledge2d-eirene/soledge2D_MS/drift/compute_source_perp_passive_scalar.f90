subroutine compute_source_perp_passive_scalar(zone,SP,Fluxes_SP,Nx,Nz)
  use all_variables, only : reference_parameters, global_variables, flags, drift_flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(in) :: SP(0:Nx+1,0:Nz+1)
  real*8,intent(out) :: Fluxes_SP(1:Nx,1:Nz,1:4)
  real*8,allocatable :: Flu(:,:)
  real*8,allocatable :: SPred(:,:)
  real*8,allocatable :: SPm(:,:),SPp(:,:)
  real*8,allocatable :: vpolm(:,:),vpolp(:,:)
  real*8,allocatable :: Fbot(:,:),Ftop(:,:)
  real*8,allocatable :: lambdam(:,:),lambdap(:,:)
  real*8,allocatable :: fluxm(:,:),fluxp(:,:)
  real*8,allocatable :: Varm(:,:),Varp(:,:)
  real*8,allocatable :: Jmed(:,:),scppmed(:,:)
  real*8 :: gam
  integer*4 :: i,j,l,k
  real*8 :: R0,rs0,c0
  real*8 :: norm_epsin, norm_epsis

  allocate(Flu(0:Nx,1:Nz))
  allocate(SPm(0:Nx,1:Nz),SPp(0:Nx,1:Nz),SPred(0:Nx+1,1:Nz))
  allocate(vpolm(0:Nx,1:Nz),vpolp(0:Nx,1:Nz))
  allocate(Fbot(0:Nx,1:Nz),Ftop(0:Nx,1:Nz))
  allocate(lambdam(0:Nx,1:Nz),lambdap(0:Nx,1:Nz))
  allocate(fluxm(0:Nx,1:Nz),fluxp(0:Nx,1:Nz))
  allocate(Varm(0:Nx,1:Nz),Varp(0:Nx,1:Nz))
  allocate(Jmed(0:Nx,1:Nz),scppmed(0:Nx,1:Nz))

  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  SPred=SP(0:Nx+1,1:Nz)

  SPm=SPred(0:Nx,:)
  SPp=SPred(1:Nx+1,:)

  vpolm=0.D0
  vpolp=0.D0
  if(drift_flags%use_ExB_radial) then
     vpolm=vpolm+zone%species(1)%drifts%uEp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(1)%drifts%uEp(1:Nx+1,1:Nz)
  end if
  if(drift_flags%use_gradB_radial) then
     vpolm=vpolm+zone%species(1)%drifts%uBp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(1)%drifts%uBp(1:Nx+1,1:Nz)
  end if

  vpolm=(vpolm+vpolp)*0.5D0
  vpolp=vpolm


  Jmed=(zone%metric_coefficients%jacobian(0:Nx,1:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx+1,1:Nz))*0.5d0
  scppmed=(sqrt(zone%metric_coefficients%cpp(0:Nx,1:Nz))&
       +sqrt(zone%metric_coefficients%cpp(1:Nx+1,1:Nz)))*0.5d0

  vpolm=vpolm*scppmed*(2.d0*pi*R0/rs0)
  vpolp=vpolp*scppmed*(2.d0*pi*R0/rs0)

  call eigenvalues_p_SP(Nx,Nz,lambdam,vpolm)
  call eigenvalues_p_SP(Nx,Nz,lambdap,vpolp)
  !computation of the raw flux
  call express_p_SP(Nx,Nz,SPm,fluxm,Jmed,vpolm)
  call express_p_SP(Nx,Nz,SPp,fluxp,Jmed,vpolp)
  do i=0,Nx
     do j=1,Nz
        Varm(i,j)=SPm(i,j)*Jmed(i,j)
        Varp(i,j)=SPp(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_p_e(Nx,Nz,lambdam,lambdap,fluxm,fluxp,Varm,Varp,Fbot,Ftop)

  DO i=0,Nx
     do j=1,Nz
           FLU(i,j)=Fbot(i,j)/Jmed(i,j)+Ftop(i,j)/Jmed(i,j)
     end do
  ENDDO

  Fluxes_SP=0.D0

  do i=1,Nx
     do j=1,Nz
        norm_epsin=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i+1,j)))
        norm_epsis=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i-1,j)))
        Fluxes_SP(i,j,1)=FLU(i,j)/norm_epsin
        Fluxes_SP(i,j,2)=FLU(i-1,j)/norm_epsis
     end do
  end do

end subroutine compute_source_perp_passive_scalar
