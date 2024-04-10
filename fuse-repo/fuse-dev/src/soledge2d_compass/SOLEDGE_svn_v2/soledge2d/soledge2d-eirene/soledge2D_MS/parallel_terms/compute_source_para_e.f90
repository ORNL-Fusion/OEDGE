subroutine compute_source_para_e(zone,Fluxes_E,Nx,Nz)
  use all_variables, only : reference_parameters, global_variables, flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: Nx,Nz
  real*8,intent(out) :: Fluxes_E(1:Nx,1:Nz,1:4)
  real*8, allocatable :: Flu(:,:)
  real*8, allocatable :: density(:,:),Gamma(:,:)
  real*8, allocatable :: T(:,:)
  real*8, allocatable :: densitym(:,:),densityp(:,:)
  real*8, allocatable :: Gammam(:,:),Gammap(:,:)
  real*8, allocatable :: Tm(:,:),Tp(:,:)
  real*8, allocatable :: Em(:,:),Ep(:,:)
  real*8,allocatable :: vpolm(:,:),vpolp(:,:)
  real*8, allocatable :: uEtm(:,:),uEtp(:,:)
  real*8, allocatable :: uBtm(:,:),uBtp(:,:) 
  real*8, allocatable :: Fleft(:,:),Fright(:,:)
  real*8, allocatable :: lambdam(:,:),lambdap(:,:)
  real*8, allocatable :: fluxm(:,:),fluxp(:,:)
  real*8, allocatable :: Varm(:,:),Varp(:,:)
  real*8, allocatable :: Jm(:,:),Jp(:,:),Gm(:,:),Gp(:,:)
  real*8, allocatable :: Jmed(:,:),Gmed(:,:),scttmed(:,:)
  real*8, allocatable :: density_extend(:,:),Gamma_extend(:,:),T_extend(:,:)
  real*8, allocatable :: uEt_extend(:,:),uBt_extend(:,:)  
  real*8 :: gam
  integer*4 :: i,j,l
  real*8 :: R0,rs0,c0
  real*8 :: norm_ethetae, norm_ethetaw

  allocate(Flu(1:Nx,0:Nz))
  allocate(density(1:Nx,0:Nz+1),Gamma(1:Nx,0:Nz+1),T(1:Nx,0:Nz+1))
  allocate(densitym(1:Nx,0:Nz),densityp(1:Nx,0:Nz))
  allocate(Gammam(1:Nx,0:Nz),Gammap(1:Nx,0:Nz))
  allocate(Tm(1:Nx,0:Nz),Tp(1:Nx,0:Nz))
  allocate(Em(1:Nx,0:Nz),Ep(1:Nx,0:Nz))
  allocate(vpolm(1:Nx,0:Nz),vpolp(1:Nx,0:Nz))
  allocate(uEtm(1:Nx,0:Nz),uEtp(1:Nx,0:Nz))
  allocate(uBtm(1:Nx,0:Nz),uBtp(1:Nx,0:Nz))  
  allocate(Fleft(1:Nx,0:Nz),Fright(1:Nx,0:Nz))
  allocate(lambdam(1:Nx,0:Nz),lambdap(1:Nx,0:Nz))
  allocate(fluxm(1:Nx,0:Nz),fluxp(1:Nx,0:Nz))
  allocate(Varm(1:Nx,0:Nz),Varp(1:Nx,0:Nz))
  allocate(Jm(1:Nx,0:Nz),Jp(1:Nx,0:Nz),Gm(1:Nx,0:Nz),Gp(1:Nx,0:Nz))
  allocate(Jmed(1:Nx,0:Nz),Gmed(1:Nx,0:Nz),scttmed(1:Nx,0:Nz))

  gam=5.d0/3.d0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  density=zone%species(0)%var(1)%density(1:Nx,0:Nz+1)
  Gamma=zone%species(0)%var(1)%Gamma(1:Nx,0:Nz+1)
  T=zone%species(0)%var(1)%temperature(1:Nx,0:Nz+1)

  if(flags%use_weno) then
     allocate(density_extend(1:Nx,-1:Nz+2))
     allocate(Gamma_extend(1:Nx,-1:Nz+2))
     allocate(T_extend(1:Nx,-1:Nz+2))
     allocate(uEt_extend(1:Nx,-1:Nz+2))
     allocate(uBt_extend(1:Nx,-1:Nz+2))
     call extend_fields_before_weno(zone%number,0,Nx,Nz,density_extend,Gamma_extend,T_extend,uEt_extend,uBt_extend)
     call weno2_unif(Nx,Nz,density_extend,densitym,densityp)
     call weno2_unif(Nx,Nz,Gamma_extend,Gammam,Gammap)
     call weno2_unif(Nx,Nz,T_extend,Tm,Tp)
     call weno2_unif(Nx,Nz,uEt_extend,uEtm,uEtp)
     call weno2_unif(Nx,Nz,uBt_extend,uBtm,uBtp)
     vpolm=(uEtm+uBtm)
     vpolp=(uEtp+uBtp)
     deallocate(density_extend,Gamma_extend,T_extend,uEt_extend,uBt_extend)
  else
     densitym=density(:,0:Nz)
     densityp=density(:,1:Nz+1)
     Gammam=Gamma(:,0:Nz)
     Gammap=Gamma(:,1:Nz+1)
     Tm=T(:,0:Nz)
     Tp=T(:,1:Nz+1)
     vpolm=(zone%species(0)%drifts%uEt(1:Nx,0:Nz)&
          +zone%species(0)%drifts%uBt(1:Nx,0:Nz))
     vpolp=(zone%species(0)%drifts%uEt(1:Nx,1:Nz+1)&
          +zone%species(0)%drifts%uBt(1:Nx,1:Nz+1))
  end if

  vpolm=(vpolm+vpolp)*0.5D0
  vpolp=vpolm

  Em=3.d0/2.d0*densitym*Tm
  Ep=3.d0/2.d0*densityp*Tp

  Gmed=(zone%metric_coefficients%G(1:Nx,0:Nz)&
       +zone%metric_coefficients%G(1:Nx,1:Nz+1))*0.5d0
  Jmed=(zone%metric_coefficients%jacobian(1:Nx,0:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx,1:Nz+1))*0.5d0

  scttmed=(sqrt(zone%metric_coefficients%ctt(1:Nx,0:Nz))&
       +sqrt(zone%metric_coefficients%ctt(1:Nx,1:Nz+1)))*0.5d0

  vpolm=vpolm*scttmed*(2.d0*pi*R0/rs0)
  vpolp=vpolp*scttmed*(2.d0*pi*R0/rs0)


  call eigenvalues_e(Nx,Nz,densitym,Gammam,lambdam,Gmed,vpolm)
  call eigenvalues_e(Nx,Nz,densityp,Gammap,lambdap,Gmed,vpolp)
  !computation of the raw flux
  call express_e(Nx,Nz,densitym,Gammam,Tm,fluxm,Gmed,Jmed,vpolm)
  call express_e(Nx,Nz,densityp,Gammap,Tp,fluxp,Gmed,Jmed,vpolp)
  do i=1,Nx
     do j=0,Nz
        Varm(i,j)=Em(i,j)*Jmed(i,j)
        Varp(i,j)=Ep(i,j)*Jmed(i,j)
     end do
  end do
  call ROE_e(Nx,Nz,lambdam,lambdap,fluxm,fluxp,Varm,Varp,Fleft,Fright)

  DO i=1,Nx
     do j=0,Nz
        FLU(i,j)=(Fleft(i,j)+Fright(i,j))/Jmed(i,j)
     end do
  ENDDO

  Fluxes_E=0.D0
  do i=1,Nx
     do j=1,Nz
        norm_ethetae=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j+1)))
        norm_ethetaw=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j-1)))
        Fluxes_E(i,j,3)=FLU(i,j)/norm_ethetae
        Fluxes_E(i,j,4)=FLU(i,j-1)/norm_ethetaw
     end do
  end do

end subroutine compute_source_para_e
