subroutine compute_source_perp(zone,Fluxes_n,Fluxes_G,Fluxes_E,Nx,Nz,n_ion)
  use all_variables, only : reference_parameters, global_variables, flags, drift_flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz,n_ion
  real*8,intent(out) :: Fluxes_n(1:Nx,1:Nz,1:4)
  real*8,intent(out) :: Fluxes_G(1:Nx,1:Nz,1:4)
  real*8,intent(out) :: Fluxes_E(1:Nx,1:Nz,1:4)
  real*8,allocatable :: Flu(:,:,:)
  real*8,allocatable :: density(:,:),Gamma(:,:)
  real*8,allocatable :: T(:,:)
  real*8,allocatable :: densitym(:,:),densityp(:,:)
  real*8,allocatable :: Gammam(:,:),Gammap(:,:)
  real*8,allocatable :: Tm(:,:),Tp(:,:)
  real*8,allocatable :: Em(:,:),Ep(:,:)
  real*8,allocatable :: vpolm(:,:),vpolp(:,:)
  real*8,allocatable :: Fbot(:,:,:),Ftop(:,:,:)
  real*8,allocatable :: lambdam(:,:,:),lambdap(:,:,:)
  real*8,allocatable :: fluxm(:,:,:),fluxp(:,:,:)
  real*8,allocatable :: Varm(:,:,:),Varp(:,:,:)
  real*8,allocatable :: Jmed(:,:),scppmed(:,:)
  real*8 :: gam
  real*8 :: mass
  integer*4 :: i,j,l,k
  real*8 :: R0,rs0,c0
  real*8 :: norm_epsin, norm_epsis

  allocate(Flu(0:Nx,1:Nz,1:3))
  allocate(density(0:Nx+1,1:Nz),Gamma(0:Nx+1,1:Nz),T(0:Nx+1,1:Nz))
  allocate(densitym(0:Nx,1:Nz),densityp(0:Nx,1:Nz))
  allocate(Gammam(0:Nx,1:Nz),Gammap(0:Nx,1:Nz))
  allocate(Tm(0:Nx,1:Nz),Tp(0:Nx,1:Nz))
  allocate(Em(0:Nx,1:Nz),Ep(0:Nx,1:Nz))
  allocate(vpolm(0:Nx,1:Nz),vpolp(0:Nx,1:Nz))
  allocate(Fbot(0:Nx,1:Nz,1:3),Ftop(0:Nx,1:Nz,1:3))
  allocate(lambdam(0:Nx,1:Nz,1:3),lambdap(0:Nx,1:Nz,1:3))
  allocate(fluxm(0:Nx,1:Nz,1:3),fluxp(0:Nx,1:Nz,1:3))
  allocate(Varm(0:Nx,1:Nz,1:3),Varp(0:Nx,1:Nz,1:3))
  allocate(Jmed(0:Nx,1:Nz),scppmed(0:Nx,1:Nz))

  gam=5.d0/3.d0
  mass=zone%species(n_ion)%element%mass
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  density=zone%species(n_ion)%var(1)%density(0:Nx+1,1:Nz)
  Gamma=zone%species(n_ion)%var(1)%Gamma(0:Nx+1,1:Nz)
  T=zone%species(n_ion)%var(1)%temperature(0:Nx+1,1:Nz)

  densitym=density(0:Nx,:)
  densityp=density(1:Nx+1,:)
  Gammam=Gamma(0:Nx,:)
  Gammap=Gamma(1:Nx+1,:)
  Tm=T(0:Nx,:)
  Tp=T(1:Nx+1,:)
  vpolm=0.D0
  vpolp=0.D0
  if(drift_flags%use_ExB_radial) then
     vpolm=vpolm+zone%species(n_ion)%drifts%uEp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(n_ion)%drifts%uEp(1:Nx+1,1:Nz)
  end if
  if(drift_flags%use_gradB_radial) then
     vpolm=vpolm+zone%species(n_ion)%drifts%uBp(0:Nx,1:Nz)
     vpolp=vpolp+zone%species(n_ion)%drifts%uBp(1:Nx+1,1:Nz)
  end if

  vpolm=(vpolm+vpolp)*0.5d0
  vpolp=vpolm

  Em=3.d0/2.d0*densitym*Tm&
       +1.d0/2.d0*Gammam*Gammam/densitym*mass
  Ep=3.d0/2.d0*densityp*Tp&
       +1.d0/2.d0*Gammap*Gammap/densityp*mass

  Jmed=(zone%metric_coefficients%jacobian(0:Nx,1:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx+1,1:Nz))*0.5d0
  scppmed=(sqrt(zone%metric_coefficients%cpp(0:Nx,1:Nz))&
       +sqrt(zone%metric_coefficients%cpp(1:Nx+1,1:Nz)))*0.5d0

  vpolm=vpolm*scppmed*(2.d0*pi*R0/rs0)
  vpolp=vpolp*scppmed*(2.d0*pi*R0/rs0)

  call eigenvalues_p(Nx,Nz,lambdam,vpolm)
  call eigenvalues_p(Nx,Nz,lambdap,vpolp)
  !computation of the raw flux
  call express_p(Nx,Nz,densitym,Gammam,Tm,fluxm,Jmed,mass,vpolm)
  call express_p(Nx,Nz,densityp,Gammap,Tp,fluxp,Jmed,mass,vpolp)
  do i=0,Nx
     do j=1,Nz
        Varm(i,j,1)=densitym(i,j)*Jmed(i,j)
        Varp(i,j,1)=densityp(i,j)*Jmed(i,j)
        Varm(i,j,2)=Gammam(i,j)*Jmed(i,j)
        Varp(i,j,2)=Gammap(i,j)*Jmed(i,j)
        Varm(i,j,3)=Em(i,j)*Jmed(i,j)
        Varp(i,j,3)=Ep(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_p(Nx,Nz,lambdam,lambdap,fluxm,fluxp,varm,varp,Fbot,Ftop)

  DO i=0,Nx
     do j=1,Nz
        do l=1,3
           FLU(i,j,l)=Fbot(i,j,l)/Jmed(i,j)+Ftop(i,j,l)/Jmed(i,j)
        end do
     end do
  ENDDO

  Fluxes_n=0.D0
  Fluxes_G=0.D0
  Fluxes_E=0.D0

  do i=1,Nx
     do j=1,Nz
        norm_epsin=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i+1,j)))
        norm_epsis=0.5d0*(sqrt(zone%metric_coefficients%cpp(i,j))+sqrt(zone%metric_coefficients%cpp(i-1,j)))
        Fluxes_n(i,j,1)=FLU(i,j,1)/norm_epsin
        Fluxes_n(i,j,2)=FLU(i-1,j,1)/norm_epsis
        Fluxes_G(i,j,1)=FLU(i,j,2)/norm_epsin
        Fluxes_G(i,j,2)=FLU(i-1,j,2)/norm_epsis
        Fluxes_E(i,j,1)=FLU(i,j,3)/norm_epsin
        Fluxes_E(i,j,2)=FLU(i-1,j,3)/norm_epsis
     end do
  end do


end subroutine compute_source_perp
