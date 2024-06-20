subroutine compute_source_para(zone,Fluxes_n,Fluxes_G,Fluxes_E,Nx,Nz,n_ion)
  use all_variables, only : reference_parameters, global_variables, flags
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
  real*8,allocatable :: uEtm(:,:),uEtp(:,:)
  real*8,allocatable :: uBtm(:,:),uBtp(:,:) 
  real*8,allocatable :: revm(:,:,:,:),revp(:,:,:,:)
  real*8,allocatable :: levm(:,:,:,:),levp(:,:,:,:)
  real*8,allocatable :: Eigenm(:,:,:),Eigenp(:,:,:)
  real*8,allocatable :: Fleft(:,:,:),Fright(:,:,:)
  real*8,allocatable :: lambdam(:,:,:),lambdap(:,:,:)
  real*8,allocatable :: fluxm(:,:,:),fluxp(:,:,:)
  real*8,allocatable :: Varm(:,:,:),Varp(:,:,:)
  real*8,allocatable :: Omm(:,:,:),Omp(:,:,:)
  real*8,allocatable :: Zetam(:,:,:),Zetap(:,:,:)
  real*8,allocatable :: Jm(:,:),Jp(:,:),Gm(:,:),Gp(:,:)
  real*8,allocatable :: Jmed(:,:),Gmed(:,:),scttmed(:,:)
  real*8 :: gam
  real*8 :: mass
  integer*4 :: i,j,l,k
  real*8 :: R0,rs0,c0
  real*8 :: norm_ethetae, norm_ethetaw
  real*8,allocatable :: density_extend(:,:),Gamma_extend(:,:),T_extend(:,:)
  real*8,allocatable :: uEt_extend(:,:),uBt_extend(:,:)

  allocate(Flu(1:Nx,0:Nz,1:3))
  allocate(density(1:Nx,0:Nz+1),Gamma(1:Nx,0:Nz+1),T(1:Nx,0:Nz+1))
  allocate(densitym(1:Nx,0:Nz),densityp(1:Nx,0:Nz))
  allocate(Gammam(1:Nx,0:Nz),Gammap(1:Nx,0:Nz))
  allocate(Tm(1:Nx,0:Nz),Tp(1:Nx,0:Nz))
  allocate(Em(1:Nx,0:Nz),Ep(1:Nx,0:Nz))
  allocate(vpolm(1:Nx,0:Nz),vpolp(1:Nx,0:Nz))
  allocate(uEtm(1:Nx,0:Nz),uEtp(1:Nx,0:Nz))
  allocate(uBtm(1:Nx,0:Nz),uBtp(1:Nx,0:Nz))
  allocate(revm(1:Nx,0:Nz,1:3,1:3),revp(1:Nx,0:Nz,1:3,1:3))
  allocate(levm(1:Nx,0:Nz,1:3,1:3),levp(1:Nx,0:Nz,1:3,1:3))
  allocate(Eigenm(1:Nx,0:Nz,1:3),Eigenp(1:Nx,0:Nz,1:3))
  allocate(Fleft(1:Nx,0:Nz,1:3),Fright(1:Nx,0:Nz,1:3))
  allocate(lambdam(1:Nx,0:Nz,1:3),lambdap(1:Nx,0:Nz,1:3))
  allocate(fluxm(1:Nx,0:Nz,1:3),fluxp(1:Nx,0:Nz,1:3))
  allocate(Varm(1:Nx,0:Nz,1:3),Varp(1:Nx,0:Nz,1:3))
  allocate(Omm(1:Nx,0:Nz,1:3),Omp(1:Nx,0:Nz,1:3))
  allocate(Zetam(1:Nx,0:Nz,1:3),Zetap(1:Nx,0:Nz,1:3))
  allocate(Jm(1:Nx,0:Nz),Jp(1:Nx,0:Nz),Gm(1:Nx,0:Nz),Gp(1:Nx,0:Nz))
  allocate(Jmed(1:Nx,0:Nz),Gmed(1:Nx,0:Nz),scttmed(1:Nx,0:Nz))

  gam=5.d0/3.d0
  mass=zone%species(n_ion)%element%mass
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  c0=reference_parameters%fields%c0
  density=zone%species(n_ion)%var(1)%density(1:Nx,0:Nz+1)
  Gamma=zone%species(n_ion)%var(1)%Gamma(1:Nx,0:Nz+1)
  T=zone%species(n_ion)%var(1)%temperature(1:Nx,0:Nz+1)

  if(flags%use_weno) then
     allocate(density_extend(1:Nx,-1:Nz+2))
     allocate(Gamma_extend(1:Nx,-1:Nz+2))
     allocate(T_extend(1:Nx,-1:Nz+2))
     allocate(uEt_extend(1:Nx,-1:Nz+2))
     allocate(uBt_extend(1:Nx,-1:Nz+2))
     call extend_fields_before_weno(zone%number,n_ion,Nx,Nz,density_extend,Gamma_extend,T_extend,uEt_extend,uBt_extend)
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
     vpolm=(zone%species(n_ion)%drifts%uEt(1:Nx,0:Nz)&
          +zone%species(n_ion)%drifts%uBt(1:Nx,0:Nz))
     vpolp=(zone%species(n_ion)%drifts%uEt(1:Nx,1:Nz+1)&
          +zone%species(n_ion)%drifts%uBt(1:Nx,1:Nz+1))

  end if

  vpolm=(vpolm+vpolp)*0.5D0
  vpolp=vpolm

  Em=3.d0/2.d0*densitym*Tm&
       +1.d0/2.d0*Gammam*Gammam/densitym*mass
  Ep=3.d0/2.d0*densityp*Tp&
       +1.d0/2.d0*Gammap*Gammap/densityp*mass


  Gmed=(zone%metric_coefficients%G(1:Nx,0:Nz)&
       +zone%metric_coefficients%G(1:Nx,1:Nz+1))*0.5d0
  Jmed=(zone%metric_coefficients%jacobian(1:Nx,0:Nz)&
       +zone%metric_coefficients%jacobian(1:Nx,1:Nz+1))*0.5d0
  scttmed=(sqrt(zone%metric_coefficients%ctt(1:Nx,0:Nz))&
       +sqrt(zone%metric_coefficients%ctt(1:Nx,1:Nz+1)))*0.5d0

  vpolm=vpolm*scttmed*(2.d0*pi*R0/rs0)
  vpolp=vpolp*scttmed*(2.d0*pi*R0/rs0)

  call eigenvalues(Nx,Nz,densitym,Gammam,Tm,lambdam,mass,Gmed,vpolm)
  call eigenvalues(Nx,Nz,densityp,Gammap,Tp,lambdap,mass,Gmed,vpolp)
  !computation of the right eigen vectors
  call reigvec(Nx,Nz,densitym,Gammam,Tm,revm,mass)
  call reigvec(Nx,Nz,densityp,Gammap,Tp,revp,mass)
  !computation of the left eigen vectors
  call leigvec(Nx,Nz,densitym,Gammam,Tm,levm,mass)
  call leigvec(Nx,Nz,densityp,Gammap,Tp,levp,mass)
  !computation of the raw flux
  call express(Nx,Nz,densitym,Gammam,Tm,fluxm,Gmed,Jmed,mass,vpolm)
  call express(Nx,Nz,densityp,Gammap,Tp,fluxp,Gmed,Jmed,mass,vpolp)
  do i=1,Nx
     do j=0,Nz
        Varm(i,j,1)=densitym(i,j)*Jmed(i,j)
        Varp(i,j,1)=densityp(i,j)*Jmed(i,j)
        Varm(i,j,2)=Gammam(i,j)*Jmed(i,j)
        Varp(i,j,2)=Gammap(i,j)*Jmed(i,j)
        Varm(i,j,3)=Em(i,j)*Jmed(i,j)
        Varp(i,j,3)=Ep(i,j)*Jmed(i,j)
     end do
  end do
  call projection(Nx,Nz,levm,fluxm,ZETAM)
  call projection(Nx,Nz,levp,fluxp,ZETAP)
  call projection(Nx,Nz,levm,Varm,OMM)
  call projection(Nx,Nz,levp,Varp,OMP)
  call ROE(Nx,Nz,lambdam,lambdap,ZETAM,ZETAP,OMM,OMP,EIGENM,EIGENP)
  CALL PROJECTION(Nx,Nz,revm,EIGENM,Fleft)
  CALL PROJECTION(Nx,Nz,revp,EIGENP,Fright)

  DO i=1,Nx
     do j=0,Nz
        do l=1,3
           FLU(i,j,l)=Fleft(i,j,l)/Jmed(i,j)+Fright(i,j,l)/Jmed(i,j)
        end do
     end do
  ENDDO

  Fluxes_n=0.D0
  Fluxes_G=0.D0
  Fluxes_E=0.D0

  do i=1,Nx
     do j=1,Nz
        norm_ethetae=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j+1)))
        norm_ethetaw=0.5d0*(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j-1)))
        Fluxes_n(i,j,3)=FLU(i,j,1)/norm_ethetae
        Fluxes_n(i,j,4)=FLU(i,j-1,1)/norm_ethetaw
        Fluxes_G(i,j,3)=FLU(i,j,2)/norm_ethetae
        Fluxes_G(i,j,4)=FLU(i,j-1,2)/norm_ethetaw
        Fluxes_E(i,j,3)=FLU(i,j,3)/norm_ethetae
        Fluxes_E(i,j,4)=FLU(i,j-1,3)/norm_ethetaw
     end do
  end do

  call correct_BC_sources(zone,Fluxes_n,Fluxes_G,Fluxes_E,Nx,Nz,n_ion,fluxm,fluxp)

end subroutine compute_source_para
