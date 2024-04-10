subroutine compute_vtheta_shear(zone)
  use all_variables, only : global_parameters, reference_parameters
  use Moperator
  use MZone
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: i,j
  real*8 :: uebot,uetop
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  ! dd : S_star = S /(c0^2/rs0^2)
  do i=1,Nx
     do j=1,Nz
        !top
        uetop=(zone%species(1)%drifts%uet(i+1,j)-zone%species(1)%drifts%uet(i,j))&
             /(zone%mesh%x(i+1,j)-zone%mesh%x(i,j))*&
             (zone%metric_coefficients%cpp(i,j)+zone%metric_coefficients%cpp(i+1,j))*0.5d0&
             +(zone%species(1)%drifts%uet(i,j+1)-zone%species(1)%drifts%uet(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*&
             zone%metric_coefficients%cpt(i,j)
        !bot
        uebot=(zone%species(1)%drifts%uet(i,j)-zone%species(1)%drifts%uet(i-1,j))&
             /(zone%mesh%x(i,j)-zone%mesh%x(i-1,j))*&
             (zone%metric_coefficients%cpp(i,j)+zone%metric_coefficients%cpp(i-1,j))*0.5d0&
             +(zone%species(1)%drifts%uet(i,j+1)-zone%species(1)%drifts%uet(i,j-1))&
             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*&
             zone%metric_coefficients%cpt(i,j)
!!$        zone%kepsilon(1)%UE_shear(i,j)=(zone%species(1)%drifts%uet(i+1,j)-zone%species(1)%drifts%uet(i-1,j))&
!!$             /(zone%mesh%x(i+1,j)-zone%mesh%x(i-1,j))*&
!!$             zone%metric_coefficients%cpp(i,j)&
!!$             +(zone%species(1)%drifts%uet(i,j+1)-zone%species(1)%drifts%uet(i,j-1))&
!!$             /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*&
!!$             zone%metric_coefficients%cpt(i,j)
        zone%kepsilon(1)%UE_shear(i,j)=(uebot*uebot+uetop*uetop)*0.5d0
        zone%kepsilon(1)%UE_shear(i,j)=zone%kepsilon(1)%UE_shear(i,j)/&
             zone%metric_coefficients%cpp(i,j)
        zone%kepsilon(1)%UE_shear(i,j)=zone%kepsilon(1)%UE_shear(i,j)&
             /reference_parameters%fields%c0**2
        zone%kepsilon(1)%UE_shear(i,j)=zone%kepsilon(1)%UE_shear(i,j)*(1.D0-zone%masks%chi2(i,j))&
             *(1.D0-zone%masks%chi2(i,j+1))*(1.D0-zone%masks%chi2(i,j-1))&
             *(1.D0-zone%masks%chi2(i+1,j))*(1.D0-zone%masks%chi2(i-1,j))
     end do
  end do
end subroutine compute_vtheta_shear
