subroutine apply_density_boundary_conditions(flux_surface,n_ion)
  use all_variables, only : global_variables, zones
  use Mflux_surface
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi,Nz_N
  integer*4 :: n,izone
  real*8 :: coef
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  if(.not.flux_surface%properties%is_periodic) then
     !west BC : n(0)=n(1)
     izone=flux_surface%properties%zones(1)
     if(zones(izone)%Neighbors(4).gt.-6) then
        izone=flux_surface%properties%zones(1)
        coef=(zones(izone)%metric_coefficients%G(i_psi,2)+zones(izone)%metric_coefficients%G(i_psi,1))&
             /(zones(izone)%metric_coefficients%G(i_psi,1)+zones(izone)%metric_coefficients%G(i_psi,0))&
             *(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))&
             /(zones(izone)%mesh%z(i_psi,2)-zones(izone)%mesh%z(i_psi,1))
        Flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+&
             flux_surface%tridiag%a(1)*(1.d0+coef)
        Flux_surface%tridiag%c(1)=flux_surface%tridiag%c(1)-&
             flux_surface%tridiag%a(1)*coef
        flux_surface%tridiag%a(1)=0.d0
!        Flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+&
!             flux_surface%tridiag%a(1)
!        flux_surface%tridiag%a(1)=0.d0
     else
        flux_surface%tridiag%s(1)=flux_surface%tridiag%s(1)-flux_surface%tridiag%a(1)*&
             zones(izone)%species(n_ion)%var(1)%density(i_psi,0)
        flux_surface%tridiag%a(1)=0.d0
     end if
     !east BC : n(Nz+1)=n(Nz)
     if(zones(izone)%Neighbors(3).gt.-6) then
        izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
        Nz_N=zones(izone)%mesh%Nz
        coef=(zones(izone)%metric_coefficients%G(i_psi,Nz_N)+zones(izone)%metric_coefficients%G(i_psi,Nz_N-1))&
             /(zones(izone)%metric_coefficients%G(i_psi,Nz_N+1)+zones(izone)%metric_coefficients%G(i_psi,Nz_N))&
             *(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))&
             /(zones(izone)%mesh%z(i_psi,Nz_N)-zones(izone)%mesh%z(i_psi,Nz_N-1))
        flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+&
            flux_surface%tridiag%c(Nz)*(1.d0+coef)
        flux_surface%tridiag%a(Nz)=flux_surface%tridiag%a(Nz)-&
             flux_surface%tridiag%c(Nz)*coef
        flux_surface%tridiag%c(Nz)=0.d0
!        flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+&
!             flux_surface%tridiag%c(Nz)
!        flux_surface%tridiag%c(Nz)=0.d0
     else
        izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
        Nz_N=zones(izone)%mesh%Nz
        flux_surface%tridiag%s(Nz)=flux_surface%tridiag%s(Nz)-flux_surface%tridiag%c(Nz)*&
             zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)
        flux_surface%tridiag%c(Nz)=0.d0
     end if
  end if
end subroutine apply_density_boundary_conditions
