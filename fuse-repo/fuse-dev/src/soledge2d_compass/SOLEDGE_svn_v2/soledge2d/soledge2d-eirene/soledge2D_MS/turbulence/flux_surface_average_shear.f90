subroutine flux_surface_average_shear()
  use all_variables, only : global_variables, zones, global_parameters, flux_surfaces
  implicit none
  integer*4 :: Nz,i_psi,j
  integer*4 :: n,izone
  real*8 :: npt
  real*8 :: UEshear
  integer*4 :: npsi
  do npsi=1,global_parameters%N_psi_surfaces
     Nz=flux_surfaces(npsi)%Nz
     i_psi=flux_surfaces(npsi)%properties%i_psi
     !concat
     UEshear=0.D0
     npt=0.D0
     do n=1,flux_surfaces(npsi)%properties%n_zones
        izone=flux_surfaces(npsi)%properties%zones(n)
        do j=1,zones(izone)%mesh%Nz
           if(Zones(izone)%masks%chi2(i_psi,j).eq.0) then
              UEshear=UEshear+Zones(izone)%kepsilon(1)%UE_shear(i_psi,j)*&
                   (zones(izone)%mesh%z_plus_1half(i_psi,j)-zones(izone)%mesh%z_minus_1half(i_psi,j))/&
                   sqrt(zones(izone)%metric_coefficients%ctt(i_psi,j))
              npt=npt+&
                   (zones(izone)%mesh%z_plus_1half(i_psi,j)-zones(izone)%mesh%z_minus_1half(i_psi,j))/&
                   sqrt(zones(izone)%metric_coefficients%ctt(i_psi,j))
           end if
        end do
     end do
     if(npt.gt.0.D0) then
        UEshear=UEshear/npt
     else
        UEshear=0.D0
     end if
     !unconcat
     do n=1,flux_surfaces(npsi)%properties%n_zones
        izone=flux_surfaces(npsi)%properties%zones(n)
        do j=1,zones(izone)%mesh%Nz
           zones(izone)%kepsilon(1)%UE_shear(i_psi,j)=UEshear
        end do
     end do
  end do
end subroutine flux_surface_average_shear
