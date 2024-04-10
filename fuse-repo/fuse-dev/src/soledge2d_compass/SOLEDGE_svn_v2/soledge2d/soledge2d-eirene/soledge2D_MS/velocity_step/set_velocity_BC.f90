subroutine set_velocity_BC(flux_surface,n_ion)
  use Mflux_surface
  use all_variables, only : global_variables, zones
  implicit none
  Type(Tflux_surface),intent(inout) :: flux_surface
  integer*4,intent(in) :: n_ion
  integer*4 :: Nz
  integer*4 :: izone,i_psi
  integer*4 :: Nz_N
  real*8 :: cs
  Nz=flux_surface%Nz
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
        flux_surface%tridiag%buffer(0)=zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,0)
     else
        cs=sqrt((zones(izone)%species(0)%var(1)%temperature(i_psi,0)&
             +zones(izone)%species(n_ion)%var(1)%temperature(i_psi,0))&
             /zones(izone)%species(n_ion)%element%mass)
        if(global_variables%sign_metric.eq.1.d0) then
           if(zones(izone)%species(n_ion)%var(1)%Mach(i_psi,1).ge.-1.d0) then
              flux_surface%tridiag%buffer(0)=-cs*zones(izone)%species(n_ion)%var(1)%density(i_psi,0)
           else
              flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)
           end if
        else
           if(zones(izone)%species(n_ion)%var(1)%Mach(i_psi,1).le.1.d0) then
              flux_surface%tridiag%buffer(0)=cs*zones(izone)%species(n_ion)%var(1)%density(i_psi,0)
           else
              flux_surface%tridiag%buffer(0)=flux_surface%tridiag%buffer(1)
           end if
        end if
     end if
     !east
     izone=flux_surface%properties%zones(flux_surface%properties%n_zones)
     Nz_N=zones(izone)%mesh%Nz
     if(zones(izone)%Neighbors(3).eq.-5) then
        flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
     elseif (zones(izone)%Neighbors(3).eq.-6) then
        flux_surface%tridiag%buffer(Nz+1)=zones(izone)%species(n_ion)%var(1)%Gamma(i_psi,Nz_N+1)
     else
        cs=sqrt((zones(izone)%species(0)%var(1)%temperature(i_psi,Nz_N+1)&
             +zones(izone)%species(n_ion)%var(1)%temperature(i_psi,Nz_N+1))&
             /zones(izone)%species(n_ion)%element%mass)
        if(global_variables%sign_metric.eq.1.d0) then
           if(zones(izone)%species(n_ion)%var(1)%Mach(i_psi,Nz_N).le.1.d0) then
              flux_surface%tridiag%buffer(Nz+1)=cs*zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)
           else
              flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
           end if
        else
           if(zones(izone)%species(n_ion)%var(1)%Mach(i_psi,Nz_N).ge.-1.d0) then
              flux_surface%tridiag%buffer(Nz+1)=-cs*zones(izone)%species(n_ion)%var(1)%density(i_psi,Nz_N+1)
           else
              flux_surface%tridiag%buffer(Nz+1)=flux_surface%tridiag%buffer(Nz)
           end if
        end if
     end if
  end if
end subroutine set_velocity_BC
