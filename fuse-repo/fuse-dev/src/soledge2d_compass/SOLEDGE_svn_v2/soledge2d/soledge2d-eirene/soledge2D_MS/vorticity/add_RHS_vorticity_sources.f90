subroutine add_RHS_vorticity_sources(zone)
  use all_variables,only : global_variables
  use MZone
  use Moperator
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  real*8, allocatable :: Fluxes(:,:,:)
  real*8, allocatable :: Source(:,:)
  Nx = zone%mesh%Nx
  Nz = zone%mesh%Nz
  allocate(Fluxes(1:Nx,1:Nz,1:4))
  allocate(Source(1:Nx,1:Nz))
  Fluxes=zone%electric_fields(1)%j_para_adv_W(1:Nx,1:Nz,1:4)&
       +zone%electric_fields(1)%j_diff_W(1:Nx,1:Nz,1:4)&
       +zone%electric_fields(1)%j_perp_adv_W(1:Nx,1:Nz,1:4)
  Source=divergence(zone,Fluxes,Nx,Nz)
  zone%electric_fields(1)%RHS=zone%electric_fields(1)%RHS+Source*global_variables%dt
  deallocate(Fluxes,Source)
end subroutine add_RHS_vorticity_sources
