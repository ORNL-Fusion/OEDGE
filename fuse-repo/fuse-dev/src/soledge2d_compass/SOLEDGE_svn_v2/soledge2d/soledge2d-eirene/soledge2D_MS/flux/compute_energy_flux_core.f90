subroutine compute_energy_flux_core(zone,nion,Flux)
  use MZone
  use all_variables, only : reference_parameters
  use Mphysics
  implicit none
  Type(TZone),intent(in) :: zone
  integer*4,intent(in) :: nion
  real*8,intent(out) :: Flux
  integer*4 :: i,j,Nx,Nz
  real*8 :: rs0,R0,n0,c0,T0
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Flux=0.d0
  if(zone%Neighbors(1).eq.-1) then
     do j=1,Nz
        Flux=Flux-zone%species(nion)%fluxes%fluxE(Nx,j,1)*zone%metric_coefficients%ds_north_PU(Nx,j)&
             *n0*c0*kB*T0*rs0/(2.d0*pi*R0)
     end do
  end if
  if(zone%Neighbors(2).eq.-1) then
     do j=1,Nz
        Flux=Flux+zone%species(nion)%fluxes%fluxE(1,j,2)*zone%metric_coefficients%ds_south_PU(1,j)&
             *n0*c0*kB*T0*rs0/(2.d0*pi*R0)
     end do
  end if
end subroutine compute_energy_flux_core
