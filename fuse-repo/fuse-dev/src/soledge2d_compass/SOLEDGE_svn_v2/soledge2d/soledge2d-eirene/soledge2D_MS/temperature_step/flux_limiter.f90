subroutine flux_limiter(zone)
  use all_variables, only : global_parameters, reference_parameters
  use MZone
  use Mphysics
  implicit none
  Type(TZone),intent(inout) :: zone
  integer*4 i,j,k,Nx,Nz,n
  real*8 :: gradT
  real*8 :: c0,T0,T0eV
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  T0eV=reference_parameters%fields%T0eV
  Nx=Zone%mesh%Nx
  Nz=Zone%mesh%Nz
  do n=0,global_parameters%N_ions
     if(Zone%species(n)%transport_para%flux_limiter.ne.0.d0) then
        do i=1,Nx
           do j=1,Nz
              if(Zone%masks%chi2(i,j).eq.0) then
                 gradT=zone%metric_coefficients%G(i,j)&
                      *(Zone%species(n)%var(1)%temperature(i,j+1)-Zone%species(n)%var(1)%temperature(i,j-1))&
                      /(Zone%mesh%z(i,j+1)-Zone%mesh%z(i,j-1))
                 Zone%species(n)%transport_para%kappa(i,j)=Zone%species(n)%transport_para%kappa0(i,j)/&
                      (1.D0+1.d0/Zone%species(n)%transport_para%flux_limiter*abs(Zone%species(n)%var(1)%temperature(i,j)&
                      *gradT/Zone%species(n)%var(1)%density(i,j))&
                      *Zone%species(n)%transport_para%kappa0(i,j)/zone%species(n)%var(1)%log_Lambda(i,j)&
                      *T0eV*sqrt(zone%species(n)%element%mass2*m_u)*c0*eV/(T0**1.5d0*kB**1.5d0))
              else
                 Zone%species(n)%transport_para%kappa(i,j)=Zone%species(n)%transport_para%kappa0(i,j)
              end if
           end do
        end do
     else
        Zone%species(n)%transport_para%kappa=Zone%species(n)%transport_para%kappa0
     end if
  end do
end subroutine flux_limiter
