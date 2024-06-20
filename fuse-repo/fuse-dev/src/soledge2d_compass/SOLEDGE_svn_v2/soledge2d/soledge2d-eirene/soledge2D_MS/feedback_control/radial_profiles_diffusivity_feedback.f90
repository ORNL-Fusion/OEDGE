subroutine radial_profiles_diffusivity_feedback()
  use all_variables, only : global_parameters, reference_parameters, zones
  use MradialFeedback
  implicit none
  integer*4 :: n,k,i,j
  integer*4 :: Nx,Nxshift,Nz
  real*8 :: flux,grad,fluxE
  real*8,allocatable :: Dloc(:), Dold(:), fluxv(:), fluxEv(:), gradv(:), xv(:), gradEv(:), gradIv(:), fluxIv(:)
  real*8,allocatable :: ChiEloc(:), ChiEold(:) 
  real*8,allocatable :: ChiIloc(:), ChiIold(:)
  real*8,allocatable :: ne(:), Te(:), Ti(:)
  real*8 :: eps, D0
  logical :: chi_appeared
  eps=1.d-6
  allocate(Dloc(1:radialFeedbackData%Nxtot))
  allocate(Chieloc(1:radialFeedbackData%Nxtot))
  allocate(ChiIloc(1:radialFeedbackData%Nxtot))
  allocate(Dold(1:radialFeedbackData%Nxtot))
  allocate(ChiEold(1:radialFeedbackData%Nxtot))
  allocate(ChiIold(1:radialFeedbackData%Nxtot))
  allocate(fluxv(1:radialFeedbackData%Nxtot))
  allocate(gradv(1:radialFeedbackData%Nxtot))
  allocate(fluxEv(1:radialFeedbackData%Nxtot))
  allocate(gradEv(1:radialFeedbackData%Nxtot))
  allocate(fluxIv(1:radialFeedbackData%Nxtot))
  allocate(gradIv(1:radialFeedbackData%Nxtot))
  allocate(xv(1:radialFeedbackData%Nxtot))
  allocate(ne(1:radialFeedbackData%Nxtot))
  allocate(Te(1:radialFeedbackData%Nxtot))
  allocate(Ti(1:radialFeedbackData%Nxtot))
  Nxshift=0
  D0=(reference_parameters%geometry%rs0**2.)/reference_parameters%fields%tau0
  chi_appeared=.false.
  do n=1,radialFeedbackData%Nzones
     k=radialFeedbackData%zone_profile(n)
     Nx=zones(k)%mesh%Nx
     !recompute new local diffusivity from flux and desired gradient
     Nz=radialFeedbackData%Nz
     do i=1,Nx
        xv(i+Nxshift)=zones(k)%mesh%x(i,Nz)
        ! density
        flux=radialFeedbackData%set(n)%fluxN(i)
        fluxv(i+Nxshift)=abs(flux)
        grad=radialFeedbackData%input_gradN(i+Nxshift)*reference_parameters%geometry%rs0/reference_parameters%fields%n0
        gradv(i+Nxshift)=grad
        if(grad.ne.0.D0) then
           Dloc(i+Nxshift)=(flux)/(-grad)
        else
           Dloc(i+Nxshift)=radialFeedbackData%Dmax/D0
        end if
        Dloc(i+Nxshift)=min(Dloc(i+Nxshift),radialFeedbackData%Dmax/D0)
        Dloc(i+Nxshift)=max(Dloc(i+Nxshift),radialFeedbackData%Dmin/D0)
        Dold(i+Nxshift)=zones(k)%species(1)%transport_perp%D_p(i,Nz)
        ! Te
        flux=radialFeedbackData%set(n)%fluxTe(i)
        fluxEv(i+Nxshift)=flux
        grad=radialFeedbackData%input_gradTe(i+Nxshift)*reference_parameters%geometry%rs0/reference_parameters%fields%T0eV
        gradEv(i+Nxshift)=grad
        if(grad.ne.0.D0) then
           ChiEloc(i+Nxshift)=(flux)&
                /(-(zones(k)%species(0)%var(1)%density(i,Nz)*grad))
        else
           ChiEloc(i+Nxshift)=radialFeedbackData%Dmax/D0
        end if
        ChiEloc(i+Nxshift)=min(ChiEloc(i+Nxshift),radialFeedbackData%Dmax/D0)
        ChiEloc(i+Nxshift)=max(ChiEloc(i+Nxshift),radialFeedbackData%Dmin/D0)
        ChiEold(i+Nxshift)=zones(k)%species(0)%transport_perp%chi_p(i,Nz)
        ! Ti
        flux=radialFeedbackData%set(n)%fluxTi(i)
        fluxIv(i+Nxshift)=flux
        grad=radialFeedbackData%input_gradTi(i+Nxshift)*reference_parameters%geometry%rs0/reference_parameters%fields%T0eV
        gradIv(i+Nxshift)=grad
        if(grad.ne.0.D0) then
           ChiIloc(i+Nxshift)=(flux)&
                /(-(zones(k)%species(1)%var(1)%density(i,Nz)*grad))
        else
           ChiIloc(i+Nxshift)=radialFeedbackData%Dmax/D0
        end if
        ChiIloc(i+Nxshift)=min(ChiIloc(i+Nxshift),radialFeedbackData%Dmax/D0)
        ChiIloc(i+Nxshift)=max(ChiIloc(i+Nxshift),radialFeedbackData%Dmin/D0)
        ChiIold(i+Nxshift)=zones(k)%species(1)%transport_perp%chi_p(i,Nz)
        ne(i+Nxshift)=zones(k)%species(1)%var(1)%density(i,Nz)*reference_parameters%fields%n0
        Te(i+Nxshift)=zones(k)%species(0)%var(1)%temperature(i,Nz)*reference_parameters%fields%T0eV
        Ti(i+Nxshift)=zones(k)%species(1)%var(1)%temperature(i,Nz)*reference_parameters%fields%T0eV
        if(zones(k)%masks%chi(i,Nz).eq.1) then
           chi_appeared=.true.
        end if
        if(chi_appeared) then
           Dloc(i+Nxshift)=Dloc(i+Nxshift-1)
           ChiEloc(i+Nxshift)=ChiEloc(i+Nxshift-1)
           ChiIloc(i+Nxshift)=ChiIloc(i+Nxshift-1)
        end if
     end do
     Nxshift=Nxshift+Nx
  end do
  open(unit=10,file='diffusion',status='unknown')
  do i=1,radialFeedbackData%Nxtot
     write(10,1) xv(i), Dold(i)*D0, Dloc(i)*D0, fluxv(i), gradv(i), ChiEold(i)*D0, fluxEv(i), gradEv(i), ChiIold(i)*D0, fluxIv(i), gradIv(i)
  end do
  close(10)
  open(unit=11,file='feedbackProfs',status='unknown')
  do i=1,radialFeedbackData%Nxtot
     write(11,1) xv(i), ne(i), Te(i), Ti(i), radialFeedbackData%input_n(i), radialFeedbackData%input_Te(i), radialFeedbackData%input_Ti(i)
  end do
  close(11)
  radialFeedbackData%x=xv
  Dloc(1)=0.5/D0
  chieloc(1)=1.5/D0
  chiiloc(1)=1.5/D0
  radialFeedbackData%D=Dloc
  radialFeedbackData%chie=chieloc
  radialFeedbackData%chii=chiiloc
  deallocate(Dloc,Dold,fluxv,gradv,xv,chiEloc,chieold,chiiloc,chiiold,fluxEv,gradEv,fluxIv,gradIv)
1 format(512es15.7)
end subroutine radial_profiles_diffusivity_feedback
