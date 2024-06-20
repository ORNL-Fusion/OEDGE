subroutine apply_velocity_drift_boundary_conditions(flux_surface,n_ion)
  use all_variables, only : global_variables, zones, reference_parameters
  use Mflux_surface
  use Mphysics
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz,offset,i_psi
  integer*4 :: n,izone
  real*8 :: cs,Mach
  real*8 :: Mtheta, uE, ctt, G, upara
  integer*4 :: Nz_N
  real*8 :: R0, rs0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0 
  Nz=flux_surface%Nz
  i_psi=flux_surface%properties%i_psi
  if(.not.flux_surface%properties%is_periodic) then
     !west BC : 
     izone=flux_surface%properties%zones(1)
     if(zones(izone)%Neighbors(4).eq.-5) then
        !dG/ds = 0
        flux_surface%tridiag%b(1)=flux_surface%tridiag%b(1)+&
             flux_surface%tridiag%a(1)
        flux_surface%tridiag%a(1)=0.d0
     else
        Mach=zones(izone)%species(n_ion)%var(1)%Mach(i_psi,1)
        uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,1)
        G=zones(izone)%metric_coefficients%G(i_psi,1)
        ctt=zones(izone)%metric_coefficients%ctt(i_psi,1)
        cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,1)+&
             zones(izone)%species(0)%var(1)%temperature(i_psi,1))/zones(izone)%species(n_ion)%element%mass)
        Mtheta=(Mach+uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G/cs)
        if(global_variables%sign_metric.eq.1.d0) then
           if(Mtheta.ge.-1.d0) then !if M<-1 then M=-1
              uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,0)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,0)
              G=zones(izone)%metric_coefficients%G(i_psi,0)
              ctt=zones(izone)%metric_coefficients%ctt(i_psi,0)
              cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,0))/zones(izone)%species(n_ion)%element%mass)
              upara=-cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
              flux_surface%tridiag%S(1)=flux_surface%tridiag%S(1)-flux_surface%tridiag%a(1)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
              flux_surface%tridiag%a(1)=0.d0
           else !dM/ds = 0
              flux_surface%tridiag%S(1)=flux_surface%tridiag%S(1)+flux_surface%tridiag%a(1)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,0)&
                   *sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,0))/zones(izone)%species(n_ion)%element%mass)*abs(Mach)
              flux_surface%tridiag%a(1)=0.d0
           end if
        else
           if(Mtheta.le.1.d0) then !if M<-1 then M=-1
              uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,0)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,0)
              G=zones(izone)%metric_coefficients%G(i_psi,0)
              ctt=zones(izone)%metric_coefficients%ctt(i_psi,0)
              cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,0))/zones(izone)%species(n_ion)%element%mass)
              upara=cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
              flux_surface%tridiag%S(1)=flux_surface%tridiag%S(1)-flux_surface%tridiag%a(1)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,0)*upara
              flux_surface%tridiag%a(1)=0.d0
           else !dM/ds = 0
              flux_surface%tridiag%S(1)=flux_surface%tridiag%S(1)-flux_surface%tridiag%a(1)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,0)&
                   *sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,0))/zones(izone)%species(n_ion)%element%mass)*abs(Mach)
              flux_surface%tridiag%a(1)=0.d0
           end if
        end if
     end if
     !east BC : n(Nz+1)=n(Nz)
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     if(zones(izone)%Neighbors(3).eq.-5) then
        !dG/dS=0
        flux_surface%tridiag%b(Nz)=flux_surface%tridiag%b(Nz)+&
             flux_surface%tridiag%c(Nz)
        flux_surface%tridiag%c(Nz)=0.d0
     else
        Mach=zones(izone)%species(n_ion)%var(1)%Mach(i_psi,Nz_N)
        uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N)
        G=zones(izone)%metric_coefficients%G(i_psi,Nz_N)
        ctt=zones(izone)%metric_coefficients%ctt(i_psi,Nz_N)
        cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N)+&
             zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N))/zones(izone)%species(n_ion)%element%mass)
        Mtheta=(Mach+uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G/cs)
        if(global_variables%sign_metric.eq.1.d0) then
           if(Mtheta.le.1.d0) then
              uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N+1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N+1)
              G=zones(izone)%metric_coefficients%G(i_psi,Nz_N+1)
              ctt=zones(izone)%metric_coefficients%ctt(i_psi,Nz_N+1)
              cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1))/zones(izone)%species(n_ion)%element%mass)
              upara=cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
              flux_surface%tridiag%S(Nz)=flux_surface%tridiag%S(Nz)-flux_surface%tridiag%c(Nz)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
              flux_surface%tridiag%c(Nz)=0.d0
           else
              flux_surface%tridiag%S(Nz)=flux_surface%tridiag%S(Nz)-flux_surface%tridiag%c(Nz)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)&
                   *sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1))/zones(izone)%species(n_ion)%element%mass)*abs(Mach)
              flux_surface%tridiag%c(Nz)=0.d0
           end if
        else
           if(Mtheta.ge.-1.d0) then
              uE = zones(izone)%species(n_ion)%drifts%uEt(i_psi,Nz_N+1)+zones(izone)%species(n_ion)%drifts%uBt(i_psi,Nz_N+1)
              G=zones(izone)%metric_coefficients%G(i_psi,Nz_N+1)
              ctt=zones(izone)%metric_coefficients%ctt(i_psi,Nz_N+1)
              cs = sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1))/zones(izone)%species(n_ion)%element%mass)
              upara=-cs-uE*sqrt(ctt)*(2.d0*pi*R0/rs0)/G
              flux_surface%tridiag%S(Nz)=flux_surface%tridiag%S(Nz)-flux_surface%tridiag%c(Nz)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)*upara
              flux_surface%tridiag%c(Nz)=0.d0
           else
              flux_surface%tridiag%S(Nz)=flux_surface%tridiag%S(Nz)+flux_surface%tridiag%c(Nz)*&
                   zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)&
                   *sqrt((zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1)+&
                   zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1))/zones(izone)%species(n_ion)%element%mass)*abs(Mach)
              flux_surface%tridiag%c(Nz)=0.d0
           end if
        end if
     end if
  end if
end subroutine apply_velocity_drift_boundary_conditions
