subroutine limit_drift_velocities(zone)
  use all_variables, only : global_parameters, zones, reference_parameters, global_variables, drift_flags
  use Mzone
  use Mphysics
  implicit none
  Type(Tzone),intent(inout) :: zone
  real*8 :: Mach_lim, Mach_lim_rad
  real*8 :: R0,rs0,c0
  integer*4 :: n
  real*8,parameter :: eps = 1.d-10
  Mach_lim=drift_flags%Mach_lim
  Mach_lim_rad=drift_flags%Mach_lim_rad
  c0=reference_parameters%fields%c0
  R0=reference_parameters%geometry%R0
  rs0=reference_parameters%geometry%rs0
  do n=0,global_parameters%N_ions
     !limit on drift velocity
!!$     zone%species(n)%drifts%uEt=sign(min(abs(zone%species(n)%drifts%uEt),sqrt((zone%species(1)%var(1)%temperature&
!!$          +zone%species(0)%var(1)%temperature)/zone%species(1)%element%mass)*&
!!$          abs(zone%metric_coefficients%G)/sqrt(zone%metric_coefficients%ctt+eps)*rs0/(2.d0*pi*R0)*Mach_lim),zone%species(n)%drifts%uEt)
!!$     zone%species(n)%drifts%uBt=sign(min(abs(zone%species(n)%drifts%uBt),sqrt((zone%species(1)%var(1)%temperature&
!!$          +zone%species(0)%var(1)%temperature)/zone%species(1)%element%mass)*&
!!$          abs(zone%metric_coefficients%G)/sqrt(zone%metric_coefficients%ctt+eps)*rs0/(2.d0*pi*R0)*Mach_lim),zone%species(n)%drifts%uBt)
!!$     zone%species(n)%drifts%uEp=sign(min(abs(zone%species(n)%drifts%uEp),sqrt((zone%species(1)%var(1)%temperature&
!!$          +zone%species(0)%var(1)%temperature)/zone%species(1)%element%mass)*Mach_lim_rad),zone%species(n)%drifts%uEp)
     zone%species(n)%drifts%uEt=zone%species(n)%drifts%uEt*Mach_lim
     zone%species(n)%drifts%uEp=zone%species(n)%drifts%uEp*Mach_lim_rad
  end do
end subroutine limit_drift_velocities
