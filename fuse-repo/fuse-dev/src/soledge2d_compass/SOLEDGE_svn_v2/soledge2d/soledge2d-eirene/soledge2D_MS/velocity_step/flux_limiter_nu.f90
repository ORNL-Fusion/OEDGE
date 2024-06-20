subroutine flux_limiter_nu(zone)
  use all_variables, only : global_parameters, reference_parameters, global_variables, flags
  use MZone
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 i,j,k,Nx,Nz,n
  real*8 :: gradT
  real*8 :: c0,T0,T0eV,R0,n0,tau0
  real*8 :: Teps,fSH,fmax
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  tau0=reference_parameters%fields%tau0
  T0eV=reference_parameters%fields%T0eV
  R0=reference_parameters%geometry%R0
  Nx=Zone%mesh%Nx
  Nz=Zone%mesh%Nz
  Teps=global_variables%Teps
  !for electrons
  Zone%species(0)%transport_para%nu=0.D0
  !for ions
  do n=1,global_parameters%N_ions
     if(flags%non_zero_parallel_viscosity) then
        if(Zone%species(n)%transport_para%flux_limiter_nu.ne.0.d0) then
           do i=1,Nx
              do j=1,Nz
                 fSH=zone%species(n)%transport_para%nu0(i,j)/&
                      zone%species(n)%var(1)%log_Lambda(i,j)&
                      *max(zone%species(n)%var(1)%temperature(i,j),Teps)**2.5d0*&
                      (zone%species(n)%var(1)%velocity(i,j+1)-zone%species(n)%var(1)%velocity(i,j-1))&
                      /(zone%mesh%z(i,j+1)-zone%mesh%z(i,j-1))*zone%Metric_coefficients%G(i,j)&
                      *(2.d0*pi*R0*m_u*zone%species(n)%element%mass*n0*c0)/tau0
                 fmax=zone%species(n)%transport_para%flux_limiter_nu*zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)&
                      *n0*T0eV*eV
                 zone%species(n)%transport_para%nu(i,j)=zone%species(n)%transport_para%nu0(i,j)/&
                      (1.D0+abs(fSH/fmax))
              end do
           end do
        else
           Zone%species(n)%transport_para%nu=0.D0
        end if
     else
        Zone%species(n)%transport_para%nu=0.D0
     end if
  end do
end subroutine flux_limiter_nu
