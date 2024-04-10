subroutine compute_recycling_neutral_sources(zone)
  use all_variables,only : global_variables,reference_parameters
  use MZone
  use Mneutral_vars
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: i,j,k
  integer*4 :: Nx,Nz
  real*8 :: R0,n0,c0,rs0,tau0,DVOL
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  zone%neutrals%Sn_nn=0.d0
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  tau0=reference_parameters%fields%tau0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  do i=1,Nx
     do j=1,Nz
        DVOL=zone%metric_coefficients%dvol_DD(i,j)
        !north
        if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i+1,j).eq.1)) then
           zone%neutrals%Sn_nn(i,j)=zone%neutrals%Sn_nn(i,j)+&
                zone%species(1)%fluxes%fluxn(i,j,1)&
                *zone%metric_coefficients%ds_north_DD(i,j)&
                /DVOL
        end if
        if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i-1,j).eq.1)) then
           zone%neutrals%Sn_nn(i,j)=zone%neutrals%Sn_nn(i,j)-&
                zone%species(1)%fluxes%fluxn(i,j,2)&
                *zone%metric_coefficients%ds_south_DD(i,j)&
                /DVOL
        end if
        !east
        if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j+1).eq.1)) then
           zone%neutrals%Sn_nn(i,j)=zone%neutrals%Sn_nn(i,j)+&
                zone%species(1)%fluxes%fluxn(i,j,3)&
                *zone%metric_coefficients%ds_east_DD(i,j)&
                /DVOL
        end if
        !west
        if((zone%masks%chi2(i,j).eq.0).and.(zone%masks%chi2(i,j-1).eq.1)) then
           zone%neutrals%Sn_nn(i,j)=zone%neutrals%Sn_nn(i,j)-&
                zone%species(1)%fluxes%fluxn(i,j,4)&
                *zone%metric_coefficients%ds_west_DD(i,j)&
                /DVOL
        end if
     end do
     if(zone%neighbors(1).lt.-1) then
        DVOL=zone%metric_coefficients%dvol_DD(i,Nz)
        zone%neutrals%Sn_nn(i,Nz)=zone%neutrals%Sn_nn(i,Nz)+&
             zone%species(1)%fluxes%fluxn(i,Nz,1)&
             *zone%metric_coefficients%ds_north_DD(i,Nz)&
             /DVOL
     end if
     if(zone%neighbors(2).lt.-1) then
        DVOL=zone%metric_coefficients%dvol_DD(i,1)
        zone%neutrals%Sn_nn(i,1)=zone%neutrals%Sn_nn(i,1)-&
             zone%species(1)%fluxes%fluxn(i,1,2)&
             *zone%metric_coefficients%ds_south_DD(i,1)&
             /DVOL
     end if
     if(zone%neighbors(3).eq.-3) then
        DVOL=zone%metric_coefficients%dvol_DD(i,Nz)
        zone%neutrals%Sn_nn(i,Nz)=zone%neutrals%Sn_nn(i,Nz)+&
             zone%species(1)%fluxes%fluxn(i,Nz,3)&
             *zone%metric_coefficients%ds_east_DD(i,Nz)&
             /DVOL
     end if
     if(zone%neighbors(4).eq.-3) then
        DVOL=zone%metric_coefficients%dvol_DD(i,1)
        zone%neutrals%Sn_nn(i,1)=zone%neutrals%Sn_nn(i,1)-&
             zone%species(1)%fluxes%fluxn(i,1,4)&
             *zone%metric_coefficients%ds_west_DD(i,1)&
             /DVOL
     end if
  end do
  zone%neutrals%Sn_nn=zone%neutrals%Sn_nn*FN_recycling_coefficient
end subroutine compute_recycling_neutral_sources
