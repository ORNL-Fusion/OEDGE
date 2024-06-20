subroutine apply_temperature_drift_boundary_conditions(flux_surface,n_ion)
  use all_variables, only : global_variables, zones, reference_parameters, global_parameters
  use Mflux_surface
  use Mphysics
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: izone
  integer*4 :: Nz_N
  real*8 :: gammai,Teps,Ge
  real*8 :: R0,rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  Teps=global_variables%Teps
  if(.not.flux_surface%properties%is_periodic) then
     !west BC : 
     izone=flux_surface%properties%zones(1)
     gammai=zones(izone)%species(n_ion)%transport_para%gamma
     flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+flux_surface%tridiag%a(1)*&
          (1.d0+(gammai-2.5d0)*(zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,1)+zones(izone)%species(n_ion)%var(1)%density(i_psi,1)*&
          (zones(izone)%species(n_ion)%drifts%uEt(i_psi,1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,1))&
          *sqrt(zones(izone)%metric_coefficients%ctt(i_psi,1)*(2.d0*pi*R0/rs0))/zones(izone)%metric_coefficients%G(i_psi,1))/&
          zones(izone)%metric_coefficients%G(i_psi,1)*(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))/&
          ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1),Teps))**(2.5d0)&
          *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,1)&
          /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,1)))
     flux_surface%tridiag%a(1)=0.d0
     !east BC : 
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+flux_surface%tridiag%c(Nz)*&
          (1.d0-(gammai-2.5d0)*(zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N)+zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N)*&
          (zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N))&
          *sqrt(zones(izone)%metric_coefficients%ctt(i_psi,Nz_N)*(2.d0*pi*R0/rs0))/zones(izone)%metric_coefficients%G(i_psi,Nz_N))/&
          zones(izone)%metric_coefficients%G(i_psi,Nz_N)*(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))/&
          ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N),Teps))**(2.5d0)&
          *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,Nz_N)&
          /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,Nz_N)))
     flux_surface%tridiag%c(Nz)=0.d0
  end if
end subroutine apply_temperature_drift_boundary_conditions
