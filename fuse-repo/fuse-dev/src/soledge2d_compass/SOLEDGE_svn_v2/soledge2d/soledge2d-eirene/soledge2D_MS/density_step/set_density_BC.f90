subroutine set_density_BC(flux_surface,n_ion)
  use all_variables, only : zones
  use Mflux_surface
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,izone,i_psi,Nz_N
  real*8 :: coef
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  if(flux_surface%properties%is_periodic) then
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(1)
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(Nz)
  else
     !BC
     !east
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     if(zones(izone)%Neighbors(3).gt.-6) then
        Nz_N=zones(izone)%mesh%Nz
        coef=(zones(izone)%metric_coefficients%G(i_psi,Nz_N)+zones(izone)%metric_coefficients%G(i_psi,Nz_N-1))&
             /(zones(izone)%metric_coefficients%G(i_psi,Nz_N+1)+zones(izone)%metric_coefficients%G(i_psi,Nz_N))&
             *(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))&
             /(zones(izone)%mesh%z(i_psi,Nz_N)-zones(izone)%mesh%z(i_psi,Nz_N-1))
        flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)*(1.d0+coef)&
             -flux_surface%tridiag%buffer(Nz-1)*coef
!        flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
     else
        Nz_N=zones(izone)%mesh%Nz
        flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)
     end if
     !west
     izone=flux_surface%properties%zones(1)
     if(zones(izone)%Neighbors(4).gt.-6) then
        izone=flux_surface%properties%zones(1)
        coef=(zones(izone)%metric_coefficients%G(i_psi,2)+zones(izone)%metric_coefficients%G(i_psi,1))&
             /(zones(izone)%metric_coefficients%G(i_psi,1)+zones(izone)%metric_coefficients%G(i_psi,0))&
             *(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))&
             /(zones(izone)%mesh%z(i_psi,2)-zones(izone)%mesh%z(i_psi,1))
        flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)*(1.d0+coef)&
             -flux_surface%tridiag%buffer(2)*coef
!        flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)
    else
        flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%density(i_psi,0)
     end if
  end if
end subroutine set_density_BC
