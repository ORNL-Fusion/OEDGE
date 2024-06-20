subroutine set_velocity_drift_BC(flux_surface,n_ion)
  use Mphysics
  use Mflux_surface
  use all_variables, only : global_variables, zones, reference_parameters, global_parameters
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz
  integer*4 :: izone,i_psi
  integer*4 :: Nz_N
  real*8 :: cs,uE,G,cs1,uE1,G1,upara1,upara,Mtheta,ctt,ctt1,Mach1
  real*8 :: R0,rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  if(flux_surface%properties%is_periodic) then
     flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(1)
     flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(Nz)
  else
     !BC
     !west
     izone=flux_surface%properties%zones(1)
     uE1 = zones(izone)%species(n_ion)%drifts%uEt(i_psi,1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,1)
     G1=zones(izone)%metric_coefficients%G(i_psi,1)
     ctt1=zones(izone)%metric_coefficients%ctt(i_psi,1)
     cs1 = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1)+&
          zones(izone)%species(0)%var(1)%temperature(i_psi,1))/zones(izone)%species(n_ion)%element%mass)
     upara1=zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,1)/zones(izone)%species(n_ion)%var(1)%density(i_psi,1)
     ue = zones(izone)%species(n_ion)%drifts%uEt(i_psi,0)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,0)
     G=zones(izone)%metric_coefficients%G(i_psi,0)
     ctt=zones(izone)%metric_coefficients%ctt(i_psi,0)
     cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)+&
          zones(izone)%species(0)%var(1)%temperature(i_psi,0))/zones(izone)%species(n_ion)%element%mass)
     Mach1=upara1/cs1
     Mtheta=(upara1*G1+uE1*sqrt(ctt1)*(2.d0*pi*R0/rs0))/(cs1*G1)
     if(global_variables%sign_metric.eq.1.d0) then
        if(Mtheta.ge.-1.d0) then
           !force Mach to -1 in Ghost cell
           upara=-cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
           flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
        else
           upara=cs*Mach1
           flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
        end if
     else
        if(Mtheta.le.1.d0) then
           !force Mach to 1 in Ghost cell
           upara=cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
           flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
        else
           upara=cs*Mach1
           flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
        end if
     end if
     !east
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N+1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N+1)
     G=zones(izone)%metric_coefficients%G(i_psi,Nz_N+1)
     ctt=zones(izone)%metric_coefficients%ctt(i_psi,Nz_N+1)
     cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)+&
          zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1))/zones(izone)%species(n_ion)%element%mass)
     uE1 = zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N)
     G1=zones(izone)%metric_coefficients%G(i_psi,Nz_N)
     ctt1=zones(izone)%metric_coefficients%ctt(i_psi,Nz_N)
     cs1 = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N)+&
          zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N))/zones(izone)%species(n_ion)%element%mass)
     upara1=zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N)/zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N)
     Mach1=upara1/cs1
     Mtheta=(upara1*G1+uE1*sqrt(ctt1)*(2.d0*pi*R0/rs0))/(cs1*G1)
     if(global_variables%sign_metric.eq.1.d0) then
        if(Mtheta.le.1.d0) then
           !force Mach to 1 in Ghost cell
           upara=cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
           flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
        else
           upara=cs*Mach1
           flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
        end if
     else
        if(Mtheta.ge.-1.d0) then
           !force Mach to -1 in Ghost cell
           upara=-cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
           flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
        else
           upara=cs*Mach1
           flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
        end if
     end if
  end if
end subroutine set_velocity_drift_BC
