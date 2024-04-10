subroutine compute_source_perp_e(zone,Fluxes_E,Nx,Nz)
  use all_variables, only : reference_parameters, global_variables, flags, drift_flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(out) :: Fluxes_E(1:Nx,1:Nz,1:4)
  real*8,allocatable :: Flu(:,:)
  real*8,allocatable :: density(:,:),Gamma(:,:)
  real*8,allocatable :: T(:,:)
  real*8,allocatable :: densitym(:,:),densityp(:,:)
  real*8,allocatable :: Gammam(:,:),Gammap(:,:)
  real*8,allocatable :: Tm(:,:),Tp(:,:)
  real*8,allocatable :: Em(:,:),Ep(:,:)
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
  allocate(density(0:Nx+1,1:Nz),Gamma(0:Nx+1,1:Nz),T(0:Nx+1,1:Nz))
  allocate(densitym(0:Nx,1:Nz),densityp(0:Nx,1:Nz))
  allocate(Gammam(0:Nx,1:Nz),Gammap(0:Nx,1:Nz))
  allocate(Tm(0:Nx,1:Nz),Tp(0:Nx,1:Nz))
  allocate(Em(0:Nx,1:Nz),Ep(0:Nx,1:Nz))
  allocate(vpolm(0:Nx,1:Nz),vpolp(0:Nx,1:Nz))
  allocate(Fbot(0:Nx,1:Nz),Ftop(0:Nx,1:Nz))
  allocate(lambdam(0:Nx,1:Nz),lambdap(0:Nx,1:Nz))
  allocate(fluxm(0:Nx,1:Nz),fluxp(0:Nx,1:Nz))
  allocate(Varm(0:Nx,1:Nz),Varp(0:Nx,1:Nz))
  allocate(Jmed(0:Nx,1:Nz),scppmed(0:Nx,1:Nz))

  gam=5.d0/3.d0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  density=zone%species(0)%var(1)%density(0:Nx+1,1:Nz)
  Gamma=zone%species(0)%var(1)%Gamma(0:Nx+1,1:Nz)
  T=zone%species(0)%var(1)%temperature(0:Nx+1,1:Nz)

  densitym=density(0:Nx,:)
  densityp=density(1:Nx+1,:)
  Gammam=Gamma(0:Nx,:)
  Gammap=Gamma(1:Nx+1,:)
  Tm=T(0:Nx,:)
  Tp=T(1:Nx+1,:)

  vpolm=0.D0
  vpolp=0.D0
  if(drift_flags%use_ExB_radial) then
     vpolm=vpolm+zone%species(0)%drifts%uEp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(0)%drifts%uEp(1:Nx+1,1:Nz)
  end if
  if(drift_flags%use_gradB_radial) then
     vpolm=vpolm+zone%species(0)%drifts%uBp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(0)%drifts%uBp(1:Nx+1,1:Nz)
  end if

  vpolm=(vpolm+vpolp)*0.5d0
  vpolp=vpolm

  Em=3.d0/2.d0*densitym*Tm
  Ep=3.d0/2.d0*densityp*Tp

  Jmed=(zone%metric_coefficients%jacobian(0:Nx,1:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx+1,1:Nz))*0.5d0
  scppmed=(sqrt(zone%metric_coefficients%cpp(0:Nx,1:Nz))&
       +sqrt(zone%metric_coefficients%cpp(1:Nx+1,1:Nz)))*0.5d0

  vpolm=vpolm*scppmed*(2.d0*pi*R0/rs0)
  vpolp=vpolp*scppmed*(2.d0*pi*R0/rs0)

  call eigenvalues_p_e(Nx,Nz,lambdam,vpolm)
  call eigenvalues_p_e(Nx,Nz,lambdap,vpolp)
  !computation of the raw flux
  call express_p_e(Nx,Nz,densitym,Tm,fluxm,Jmed,vpolm)
  call express_p_e(Nx,Nz,densityp,Tp,fluxp,Jmed,vpolp)
  do i=0,Nx
     do j=1,Nz
        Varm(i,j)=Em(i,j)*Jmed(i,j)
        Varp(i,j)=Ep(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_p_e(Nx,Nz,lambdam,lambdap,fluxm,fluxp,varm,varp,fbot,ftop)

  DO i=0,Nx
     do j=1,Nz
           FLU(i,j)=Fbot(i,j)/Jmed(i,j)+Ftop(i,j)/Jmed(i,j)
     end do
  ENDDO

  Fluxes_E=0.D0

  do i=1,Nx
     do j=1,Nz
        norm_epsin=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i+1,j)))
        norm_epsis=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i-1,j)))
        Fluxes_E(i,j,1)=FLU(i,j)/norm_epsin
        Fluxes_E(i,j,2)=FLU(i-1,j)/norm_epsis
     end do
  end do

end subroutine compute_source_perp_e
