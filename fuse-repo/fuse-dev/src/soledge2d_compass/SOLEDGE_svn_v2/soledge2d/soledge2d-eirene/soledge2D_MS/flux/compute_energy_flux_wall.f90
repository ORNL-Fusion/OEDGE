subroutine compute_energy_flux_wall(zone,nion,Flux)
  use MZone
  use all_variables, only : reference_parameters, flags
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
  if(zone%Neighbors(1).eq.-2) then
     do j=1,Nz
        Flux=Flux+zone%species(nion)%fluxes%fluxE(Nx,j,1)*zone%metric_coefficients%ds_north_PU(Nx,j)&
             *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
     end do
  end if
  if(zone%Neighbors(2).eq.-2) then
     do j=1,Nz
        Flux=Flux-zone%species(nion)%fluxes%fluxE(1,j,2)*zone%metric_coefficients%ds_south_PU(1,j)&
             *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
     end do
  end if
  if(zone%Neighbors(3).eq.-3) then
     do i=1,Nx
        Flux=Flux+zone%species(nion)%fluxes%fluxE(i,Nz,3)*zone%metric_coefficients%ds_east_PU(i,Nz)&
             *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
     end do
  end if
  if(zone%Neighbors(4).eq.-3) then
     do i=1,Nx
        Flux=Flux-zone%species(nion)%fluxes%fluxE(i,1,4)*zone%metric_coefficients%ds_west_PU(i,1)&
             *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
     end do
  end if
  if(flags%is_Pen) then
     do i=1,Nx
        do j=1,Nz
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i+1,j).eq.1)) then
              Flux=Flux+zone%species(nion)%fluxes%fluxE(i,j,1)*zone%metric_coefficients%ds_north_PU(i,j)&
                   *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i-1,j).eq.1)) then
              Flux=Flux-zone%species(nion)%fluxes%fluxE(i,j,2)*zone%metric_coefficients%ds_south_PU(i,j)&
                   *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j+1).eq.1)) then
              Flux=Flux+zone%species(nion)%fluxes%fluxE(i,j,3)*zone%metric_coefficients%ds_east_PU(i,j)&
                   *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
           end if
           if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j-1).eq.1)) then
              Flux=Flux-zone%species(nion)%fluxes%fluxE(i,j,4)*zone%metric_coefficients%ds_west_PU(i,j)&
                   *n0*c0*kb*T0*rs0/(2.d0*pi*R0)
           end if
        end do
     end do
  end if
end subroutine compute_energy_flux_wall
