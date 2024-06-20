subroutine compute_explicit_diffusive_perp_source_terms(zone)
#include "compile_opt.inc"
  use all_variables, only : zones, global_parameters, flags, kepsilon_param, reference_parameters, drift_flags
  use MZone
  use MOperator
  use MDiffusion_perp
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  real*8,allocatable :: Diffusivity_p(:,:)
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8 :: Corners(2,2,2)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(Diffusivity_p(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
  do n=1,global_parameters%N_ions
     !compute source term for density equation
     Diffusivity_p=zone%species(n)%transport_perp%D_p
     Field=zone%species(n)%var(1)%density
     Corners=zone%species(n)%corners%density
     Fluxes=explicit_perp_diffusion(Field,Corners,&
          Diffusivity_p,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxn=zone%species(n)%fluxes%fluxn-Fluxes
     zone%species(n)%sources%Sn=zone%species(n)%sources%Sn+Source
     !compute source term for velocity equation
     Diffusivity_p=zone%species(n)%transport_perp%nu_p&
          *zone%species(n)%var(1)%density
     Field=zone%species(n)%var(1)%velocity
     Corners=zone%species(n)%corners%velocity
     Fluxes=explicit_perp_diffusion(Field,Corners,&
          Diffusivity_p,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxG=zone%species(n)%fluxes%fluxG-Fluxes
     zone%species(n)%sources%SG=zone%species(n)%sources%SG+Source
  end do
  do n=0,global_parameters%N_ions
     !compute source term for temperature equation
     Diffusivity_p=zone%species(n)%transport_perp%chi_p&
          *zone%species(n)%var(1)%density
     Field=zone%species(n)%var(1)%temperature
     Corners=zone%species(n)%corners%temperature
     Fluxes=explicit_perp_diffusion(Field,Corners,&
          Diffusivity_p,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxE=zone%species(n)%fluxes%fluxE-Fluxes
     zone%species(n)%sources%SE=zone%species(n)%sources%SE+Source
  end do
#if VORTICITY_PASTIX == 1
  !compute source term for vorticity equation
  Diffusivity_p=zone%species(0)%transport_perp%zeta_p
  Field=zone%electric_fields(1)%vorticity
  Corners=zone%electric_fields(1)%cornersW
  Fluxes=explicit_perp_diffusion(Field,Corners,&
       Diffusivity_p,zone,Nx,Nz)
  call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!  call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
  if(drift_flags%jdiffW) then
     zone%electric_fields(1)%j_diff_W=Fluxes
  else
     zone%electric_fields(1)%j_diff_W=0.d0
  end if
#endif
  if(flags%turbulence_model.eq.1) then
     !k 
!!$     Diffusivity_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_k&
!!$          *reference_parameters%fields%tau0*reference_parameters%fields%k0&
!!$          /(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$          *zone%species(0)%var(1)%density
     Diffusivity_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_k&
          *zone%species(0)%var(1)%density
     Field=zone%kepsilon(1)%k
     Corners=zone%kepsilon(1)%cornersK
     Fluxes=explicit_perp_diffusion(Field,Corners,&
          Diffusivity_p,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%kepsilon(1)%Sk=zone%kepsilon(1)%Sk+Source
     !epsilon
!!$     Diffusivity_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_epsilon&
!!$          *reference_parameters%fields%tau0*reference_parameters%fields%k0&
!!$          /(reference_parameters%geometry%rs0*reference_parameters%fields%c0)&
!!$          *zone%species(0)%var(1)%density
     Diffusivity_p=zone%kepsilon(1)%mu_t/kepsilon_param%sigma_epsilon&
          *zone%species(0)%var(1)%density
     Field=zone%kepsilon(1)%epsilon
     Corners=zone%kepsilon(1)%cornersEpsilon
     Fluxes=explicit_perp_diffusion(Field,Corners,&
          Diffusivity_p,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%kepsilon(1)%Sepsilon=zone%kepsilon(1)%Sepsilon+Source
  end if
  deallocate(Fluxes,Source,Diffusivity_p,Field)
end subroutine compute_explicit_diffusive_perp_source_terms
