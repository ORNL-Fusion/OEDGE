subroutine apply_temperature_boundary_conditions(flux_surface,n_ion)
  use all_variables, only : global_variables, zones, global_parameters, transport_parameters
  use Mflux_surface
  use Mphysics
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: izone
  integer*4 :: Nz_N
  real*8 :: gammai,Teps, delta_e
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  Teps=global_variables%Teps
  if(.not.flux_surface%properties%is_periodic) then
     !west BC : 
     izone=flux_surface%properties%zones(1)
     if(zones(izone)%Neighbors(4).eq.-5) then
        flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+&
             flux_surface%tridiag%a(1)
        flux_surface%tridiag%a(1)=0.d0
     elseif (zones(izone)%Neighbors(4) .eq. -6) then
        flux_surface%tridiag%s(1)=flux_surface%tridiag%s(1)-&
             flux_surface%tridiag%a(1)*zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)
     else
        if(n_ion.eq.0) then
           delta_e=transport_parameters%delta_e
           gammai=2.d0/(1.D0-delta_e)-0.5d0*log(&
                (2.D0*pi*m_e/(m_u*global_parameters%element_list(1)%mass))*&
                (1+zones(izone)%species(1)%var(1)%temperature(i_psi,1)/zones(izone)%species(0)%var(1)%temperature(i_psi,1))&
                /(1.D0-delta_e)**2.D0)
           gammai=min(max(gammai,2.5d0),40.d0)
        else
           gammai=zones(izone)%species(n_ion)%transport_para%gamma
        end if
        flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+flux_surface%tridiag%a(1)*&
             (1.d0+(gammai-2.5d0)*zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,1)/&
             zones(izone)%metric_coefficients%G(i_psi,1)*(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))/&
             ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1),Teps))**(2.5d0)&
             *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,1)&
             /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,1)))
        flux_surface%tridiag%a(1)=0.d0
     end if
     !east BC : 
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     if(zones(izone)%Neighbors(3).eq.-5) then
        flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+&
             flux_surface%tridiag%c(Nz)
        flux_surface%tridiag%c(Nz)=0.d0
     elseif (zones(izone)%Neighbors(3).eq.-6) then
        flux_surface%tridiag%s(Nz)=flux_surface%tridiag%s(Nz)-&
             flux_surface%tridiag%c(Nz)*zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)
     else
        if(n_ion.eq.0) then
           delta_e=transport_parameters%delta_e
           gammai=2.d0/(1.D0-delta_e)-0.5d0*log(&
                (2.D0*pi*m_e/(m_u*global_parameters%element_list(1)%mass))*&
                (1+zones(izone)%species(1)%var(1)%temperature(i_psi,Nz_N)/zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N))&
                /(1.D0-delta_e)**2.D0)
           gammai=min(max(gammai,2.5d0),40.d0)
        else
           gammai=zones(izone)%species(n_ion)%transport_para%gamma
        end if
        flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+flux_surface%tridiag%c(Nz)*&
             (1.d0-(gammai-2.5d0)*zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N)/&
             zones(izone)%metric_coefficients%G(i_psi,Nz_N)*(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))/&
             ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N),Teps))**(2.5d0)&
             *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,Nz_N)&
             /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,Nz_N)))
        flux_surface%tridiag%c(Nz)=0.d0
     end if
  end if
end subroutine apply_temperature_boundary_conditions
