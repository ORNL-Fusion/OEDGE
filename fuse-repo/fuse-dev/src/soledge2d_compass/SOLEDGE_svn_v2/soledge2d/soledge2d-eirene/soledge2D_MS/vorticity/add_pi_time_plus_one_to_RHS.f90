subroutine add_pi_time_plus_one_to_RHS(zone,mode_BC)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables
  use MZone
  use MOperator
  use MDiffusion_perp
  use Mvorticity_vars
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4,intent(in) :: mode_BC
  real*8,allocatable :: Diffusivity_p(:,:)
  real*8,allocatable :: Diffusivity_t(:,:)
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Fluxes2(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8 :: Corners(2,2,2)
  real*8,parameter :: eps=1.d-6
  integer*4 :: Nx,Nz,i,j
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Fluxes2(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(Diffusivity_p(0:Nx+1,0:Nz+1))
  allocate(Diffusivity_t(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
  Diffusivity_p=global_parameters%element_list(1)%mass/(zone%mesh%B**2+eps)
  Diffusivity_t=Diffusivity_p
  Field=zone%species(1)%var(2)%density*zone%species(1)%var(2)%temperature
  Corners=zone%species(1)%corners2%density*zone%species(1)%corners2%temperature
  Fluxes=explicit_perp_diffusion(Field,Corners,&
       Diffusivity_p,zone,Nx,Nz)
  Fluxes2=in_flux_surface_perp_diffusion(Field,Corners,&
       Diffusivity_p,Diffusivity_t,zone,Nx,Nz)
  Fluxes=Fluxes+Fluxes2
  if(mode_BC.eq.ZERO_FLUX) then
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz) ! or switch on and set grad phi as BC perp
!     call set_perp_flux_to_zero_on_core(zone,Fluxes,Nx,Nz)
  end if
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%electric_fields(1)%RHS=zone%electric_fields(1)%RHS-Source
  deallocate(Fluxes,Source,Fluxes2,Field,Diffusivity_p,Diffusivity_t)
end subroutine add_pi_time_plus_one_to_RHS
