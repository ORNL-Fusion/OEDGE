subroutine compute_explicit_vpinch_perp_source_terms(zone)
  use all_variables, only : zones, global_parameters
  use MZone
  use MOperator
  use Mpinch_perp
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: Nx,Nz
  real*8,allocatable :: v_pinch(:,:)
  real*8,allocatable :: Field(:,:)
  real*8,allocatable :: Fluxes(:,:,:)
  real*8,allocatable :: Source(:,:)
  real*8 :: Corners(2,2,2)
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  allocate(v_pinch(0:Nx+1,0:Nz+1))
  allocate(Field(0:Nx+1,0:Nz+1))
  do n=1,global_parameters%N_ions
     !compute source term for density equation
     v_pinch=zone%species(n)%transport_perp%v_pinch
     Field=zone%species(n)%var(1)%density
     Fluxes=explicit_perp_pinch(Field,v_pinch,zone,Nx,Nz)
     call set_perp_flux_to_zero_on_the_wall(zone,Fluxes,Nx,Nz)
     Source=divergence(zone,Fluxes,Nx,Nz)
     zone%species(n)%fluxes%fluxn=zone%species(n)%fluxes%fluxn-Fluxes
     zone%species(n)%sources%Sn=zone%species(n)%sources%Sn+Source
  end do
  deallocate(Fluxes,Source,v_pinch,Field)
end subroutine compute_explicit_vpinch_perp_source_terms
