subroutine radial_profiles_diffusivity_feedback2()
  use all_variables, only : global_parameters, reference_parameters, zones
  use MradialFeedback
  implicit none
  integer*4 :: n,k,i,j
  integer*4 :: Nx,Nxshift,Nz
  real*8 :: flux,grad,fluxE
  real*8,allocatable :: Dloc(:), xv(:)
  real*8,allocatable :: ChiEloc(:)
  real*8,allocatable :: ChiIloc(:)
  real*8,allocatable :: ne(:), Te(:), Ti(:)
  real*8 :: Dold, chieold, chiiold
  real*8 :: error,errorG
  real*8 :: Gain, keep, GainG
  real*8 :: eps
  logical :: chi_appeared
  eps=1.d-6
  allocate(Dloc(1:radialFeedbackData%Nxtot))
  allocate(Chieloc(1:radialFeedbackData%Nxtot))
  allocate(ChiIloc(1:radialFeedbackData%Nxtot))
  allocate(xv(1:radialFeedbackData%Nxtot))
  allocate(ne(1:radialFeedbackData%Nxtot))
  allocate(Te(1:radialFeedbackData%Nxtot))
  allocate(Ti(1:radialFeedbackData%Nxtot))
  Nxshift=0
  chi_appeared=.false.
  Gain=radialFeedbackData%Gain
  GainG=radialFeedbackData%GainG
  keep=radialFeedbackData%keep
  do n=1,radialFeedbackData%Nzones
     k=radialFeedbackData%zone_profile(n)
     Nx=zones(k)%mesh%Nx
     !recompute new local diffusivity from flux and desired gradient
     Nz=radialFeedbackData%Nz
     do i=1,Nx
        xv(i+Nxshift)=zones(k)%mesh%x(i,Nz)
        if(i+Nxshift.ne.1) then
           ! density
           errorG=(radialFeedbackData%input_n(i+Nxshift)-radialFeedbackData%input_n(i+Nxshift-1))/reference_parameters%fields%n0&
                -(zones(k)%species(1)%var(1)%density(i,Nz)-zones(k)%species(1)%var(1)%density(i-1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(i,Nz)-zones(k)%mesh%Rgeom(i-1,Nz))
           error=radialFeedbackData%input_n(i+Nxshift)/reference_parameters%fields%n0-zones(k)%species(1)%var(1)%density(i,Nz)
           Dloc(i+Nxshift)=2.D0*((radialFeedbackData%D(i+Nxshift)+radialFeedbackData%D(i+Nxshift-1))/2.D0+GainG*+errorG+Gain*error)-Dloc(i+Nxshift-1)
           Dloc(i+Nxshift)=min(Dloc(i+Nxshift),radialFeedbackData%Dmax)
           Dloc(i+Nxshift)=max(Dloc(i+Nxshift),radialFeedbackData%Dmin)
           ! Te
           error=radialFeedbackData%input_Te(i+Nxshift)/reference_parameters%fields%T0eV-zones(k)%species(0)%var(1)%temperature(i,Nz)
           errorG=(radialFeedbackData%input_Te(i+Nxshift)-radialFeedbackData%input_Te(i+Nxshift-1))/reference_parameters%fields%T0eV&
                -(zones(k)%species(0)%var(1)%temperature(i,Nz)-zones(k)%species(0)%var(1)%temperature(i-1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(i,Nz)-zones(k)%mesh%Rgeom(i-1,Nz))
           ChiEloc(i+Nxshift)=2.D0*((radialFeedbackData%chie(i+Nxshift)+radialFeedbackData%chie(i+Nxshift-1))/2.D0+GainG*errorG+Gain*error)-ChiEloc(i+Nxshift-1)
           ChiEloc(i+Nxshift)=min(ChiEloc(i+Nxshift),radialFeedbackData%Dmax)
           ChiEloc(i+Nxshift)=max(ChiEloc(i+Nxshift),radialFeedbackData%Dmin)
           ! Ti
           error=radialFeedbackData%input_Ti(i+Nxshift)/reference_parameters%fields%T0eV-zones(k)%species(1)%var(1)%temperature(i,Nz)
           errorG=(radialFeedbackData%input_Ti(i+Nxshift)-radialFeedbackData%input_Ti(i+Nxshift-1))/reference_parameters%fields%T0eV&
                -(zones(k)%species(1)%var(1)%temperature(i,Nz)-zones(k)%species(1)%var(1)%temperature(i-1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(i,Nz)-zones(k)%mesh%Rgeom(i-1,Nz)) 
           ChiIloc(i+Nxshift)=2.D0*((radialFeedbackData%chii(i+Nxshift)+radialFeedbackData%chii(i+Nxshift-1))/2.D0+GainG*errorG+Gain*error)-ChiIloc(i+Nxshift-1)
           ChiIloc(i+Nxshift)=min(ChiIloc(i+Nxshift),radialFeedbackData%Dmax)
           ChiIloc(i+Nxshift)=max(ChiIloc(i+Nxshift),radialFeedbackData%Dmin)
        else
           ! density
           error=radialFeedbackData%input_n(1)/reference_parameters%fields%n0-zones(k)%species(1)%var(1)%density(1,Nz)
           errorG=(radialFeedbackData%input_n(2)-radialFeedbackData%input_n(1))/reference_parameters%fields%n0&
                -(zones(k)%species(1)%var(1)%density(2,Nz)-zones(k)%species(1)%var(1)%density(1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(2,Nz)-zones(k)%mesh%Rgeom(1,Nz))
           Dloc(1)=radialFeedbackData%D(1)+GainG*errorG+Gain*error
           Dloc(1)=min(Dloc(1),radialFeedbackData%Dmax)
           Dloc(1)=max(Dloc(1),radialFeedbackData%Dmin)
           ! Te
           error=radialFeedbackData%input_Te(1)/reference_parameters%fields%T0eV-zones(k)%species(0)%var(1)%temperature(1,Nz)
           errorG=(radialFeedbackData%input_Te(2)-radialFeedbackData%input_Te(1))/reference_parameters%fields%T0eV&
                -(zones(k)%species(0)%var(1)%temperature(2,Nz)-zones(k)%species(0)%var(1)%temperature(1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(2,Nz)-zones(k)%mesh%Rgeom(1,Nz))
           ChiEloc(1)=radialFeedbackData%chie(1)+GainG*errorG+Gain*error
           ChiEloc(1)=min(ChiEloc(1),radialFeedbackData%Dmax)
           ChiEloc(1)=max(ChiEloc(1),radialFeedbackData%Dmin)
           ! Ti
           error=radialFeedbackData%input_Ti(1)/reference_parameters%fields%T0eV-zones(k)%species(1)%var(1)%temperature(1,Nz)
           errorG=(radialFeedbackData%input_Ti(2)-radialFeedbackData%input_Ti(1))/reference_parameters%fields%T0eV&
                -(zones(k)%species(1)%var(1)%temperature(2,Nz)-zones(k)%species(1)%var(1)%temperature(1,Nz))
           errorG=errorG/(zones(k)%mesh%Rgeom(2,Nz)-zones(k)%mesh%Rgeom(1,Nz))
           ChiIloc(1)=radialFeedbackData%chii(1)+Gain*error+GainG*errorG
           ChiIloc(1)=min(ChiIloc(1),radialFeedbackData%Dmax)
           ChiIloc(1)=max(ChiIloc(1),radialFeedbackData%Dmin)
        end if
        ne(i+Nxshift)=zones(k)%species(1)%var(1)%density(i,Nz)*reference_parameters%fields%n0
        Te(i+Nxshift)=zones(k)%species(0)%var(1)%temperature(i,Nz)*reference_parameters%fields%T0eV
        Ti(i+Nxshift)=zones(k)%species(1)%var(1)%temperature(i,Nz)*reference_parameters%fields%T0eV
        if(zones(k)%masks%chi(i,Nz).eq.1) then
           chi_appeared=.true.
        end if
        if(chi_appeared) then
           Dloc(i+Nxshift)=max(Dloc(i+Nxshift-1),1.D0)
           ChiEloc(i+Nxshift)=max(ChiEloc(i+Nxshift-1),1.D0)
           ChiIloc(i+Nxshift)=max(ChiIloc(i+Nxshift-1),1.D0)
        end if

     end do
     Nxshift=Nxshift+Nx
  end do

  radialFeedbackData%D=Dloc
  radialFeedbackData%chie=chieloc
  radialFeedbackData%chii=chiiloc


  open(unit=10,file='diffusion',status='unknown')
  do i=1,radialFeedbackData%Nxtot
     write(10,1) xv(i), radialFeedbackData%D(i), radialFeedbackData%chie(i), radialFeedbackData%chii(i)
  end do
  close(10)
  open(unit=11,file='feedbackProfs',status='unknown')
  do i=1,radialFeedbackData%Nxtot
     write(11,1) xv(i), ne(i), Te(i), Ti(i), radialFeedbackData%input_n(i), radialFeedbackData%input_Te(i), radialFeedbackData%input_Ti(i)
  end do
  close(11)
  radialFeedbackData%x=xv

  do k=1,global_parameters%N_zones
     call update_radial_diffusivity_feedback(zones(k))
  end do

  call MD_broadcast_transport_coefficients()

  deallocate(Dloc,xv,chieloc,chiiloc)
1 format(512es15.7)
end subroutine radial_profiles_diffusivity_feedback2
