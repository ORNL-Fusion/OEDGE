subroutine set_temperature_drift_BC(flux_surface,n_ion)
  use Mflux_surface
  use Mphysics
  use all_variables, only : global_variables, zones, reference_parameters, global_parameters
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz
  integer*4 :: izone,i_psi
  integer*4 :: Nz_N
  real*8 :: Teps,gammai,Ge
  real*8 :: R0,rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Nz=flux_surface%Nz
  Teps=global_variables%Teps
  i_psi=flux_surface%properties%i_psi
  if(flux_surface%properties%is_periodic) then
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(1)
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(Nz)
  else
     !BC
     !west
     izone=flux_surface%properties%zones(1)
     gammai=zones(izone)%species(n_ion)%transport_para%gamma
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)&
          *(1.d0+(gammai-2.5d0)*(zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,1)+zones(izone)%species(n_ion)%var(1)%density(i_psi,1)*&
          (zones(izone)%species(n_ion)%drifts%uEt(i_psi,1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,1))&
          *sqrt(zones(izone)%metric_coefficients%ctt(i_psi,1)*(2.d0*pi*R0/rs0))/zones(izone)%metric_coefficients%G(i_psi,1))/& !think about implicit
          zones(izone)%metric_coefficients%G(i_psi,1)*(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))/&
          ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1),Teps))**(2.5d0)&
          *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,1)&
          /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,1)))
     !east
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)*(1.d0-(gammai-2.5d0)&
          *(zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N)+zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N)*&
          (zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N))&
          *sqrt(zones(izone)%metric_coefficients%ctt(i_psi,Nz_N)*(2.d0*pi*R0/rs0))/zones(izone)%metric_coefficients%G(i_psi,Nz_N))/&
          zones(izone)%metric_coefficients%G(i_psi,Nz_N)*(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))/&
          ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N),Teps))**(2.5d0)&
          *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,Nz_N)&
          /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,Nz_N)))
  end if
end subroutine set_temperature_drift_BC
