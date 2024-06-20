subroutine set_temperature_BC(flux_surface,n_ion)
  use Mflux_surface
  use all_variables, only : global_variables, zones, global_parameters, transport_parameters
  use Mphysics
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz
  integer*4 :: izone,i_psi
  integer*4 :: Nz_N
  real*8 :: Teps,gammai,delta_e
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
     if(zones(izone)%Neighbors(4).eq.-5) then
        flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)
     elseif (zones(izone)%Neighbors(4) .eq. -6) then
        flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)
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
        flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)&
             *(1.d0+(gammai-2.5d0)*zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,1)/& !think about implicit
             zones(izone)%metric_coefficients%G(i_psi,1)*(zones(izone)%mesh%z(i_psi,1)-zones(izone)%mesh%z(i_psi,0))/&
             ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1),Teps))**(2.5d0)&
             *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,1)&
             /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,1)))
     end if
     !east
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     if(zones(izone)%Neighbors(3).eq.-5) then
        flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
     elseif (zones(izone)%Neighbors(3).eq.-6) then
        flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)
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
        flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)*(1.d0-(gammai-2.5d0)&
             *zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N)/&
             zones(izone)%metric_coefficients%G(i_psi,Nz_N)*(zones(izone)%mesh%z(i_psi,Nz_N+1)-zones(izone)%mesh%z(i_psi,Nz_N))/&
             ((max(zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N),Teps))**(2.5d0)&
             *zones(izone)%species(n_ion)%transport_para%kappa(i_psi,Nz_N)&
             /zones(izone)%species(n_ion)%var(1)%log_Lambda(i_psi,Nz_N)))
     end if
  end if
end subroutine set_temperature_BC
