subroutine compute_vorticity(zone,STEP,mode_BC)
#include "compile_opt.inc"
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use MZone
  use MOperator
  use MDiffusion_perp
  use MDefinitions
  use Mvorticity_vars
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: STEP
  integer*4,intent(in) :: mode_BC
  real*8,allocatable :: Diffusivity_p(:,:)
  real*8,allocatable :: Diffusivity_t(:,:)
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Fluxes2(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8 :: Corners(2,2,2)
  real*8,parameter :: eps=1.d-6
  integer*4 :: Nx,Nz
  integer*4 :: i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Fluxes2(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(Diffusivity_p(0:Nx+1,0:Nz+1))
  allocate(Diffusivity_t(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
  Fluxes=0.d0
#if VORTICITY_PI == 1
  !Pi part
  Diffusivity_p=global_parameters%element_list(1)%mass/(zone%mesh%B**2+eps)
  Diffusivity_t=Diffusivity_p
  if(STEP.eq.STEP_OLD) then
     Corners=zone%electric_fields(1)%cornersPi_old
     Field=zone%electric_fields(1)%pi_old
  else
     Corners=zone%species(1)%corners2%density*zone%species(1)%corners2%temperature
     Field=zone%species(1)%var(2)%density*zone%species(1)%var(2)%temperature
  end if
  Fluxes=explicit_perp_diffusion(Field,Corners,&
       Diffusivity_p,zone,Nx,Nz)
  Fluxes2=in_flux_surface_perp_diffusion(Field,Corners,&
       Diffusivity_p,Diffusivity_t,zone,Nx,Nz)
  Fluxes=Fluxes+Fluxes2
  if(mode_BC.eq.ZERO_FLUX) then
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!     call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
  end if
#endif
  !phi part
  Diffusivity_p=global_parameters%element_list(1)%mass*zone%species(0)%var(STEP)%density/(zone%mesh%B**2+eps)
  Diffusivity_t=Diffusivity_p
  Field=zone%electric_fields(STEP)%phi
  Corners=zone%electric_fields(STEP)%cornersPhi
  Fluxes=Fluxes+explicit_perp_diffusion(Field,Corners,&
       Diffusivity_p,zone,Nx,Nz)
  Fluxes2=in_flux_surface_perp_diffusion(Field,Corners,&
       Diffusivity_p,Diffusivity_t,zone,Nx,Nz)
  Fluxes=Fluxes+Fluxes2
  if(mode_BC.eq.ZERO_FLUX) then
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
!     call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
  end if
  if(STEP.eq.STEP_OLD) then
     zone%electric_fields(1)%j_perp=-Fluxes
  else
     zone%electric_fields(1)%j_perp=(zone%electric_fields(1)%j_perp+Fluxes)/global_variables%dt_vort
  end if
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%electric_fields(STEP)%vorticity(1:Nx,1:Nz)=Source
  deallocate(Fluxes,Source,Fluxes2,Field,Diffusivity_p,Diffusivity_t)
end subroutine compute_vorticity
