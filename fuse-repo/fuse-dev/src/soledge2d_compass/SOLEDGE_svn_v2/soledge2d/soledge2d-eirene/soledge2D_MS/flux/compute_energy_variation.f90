subroutine compute_energy_variation(zone,nion,Source)
  use MZone
  use all_variables, only : reference_parameters, flags, global_parameters, global_variables
  implicit none
  Type(TZone),intent(in) :: zone
  integer*4,intent(in) :: nion
  real*8,intent(out) :: Source
  integer*4 :: i,j,Nx,Nz,m
  real*8 :: rs0,R0,n0,c0,T0,tau0
  real*8 :: m_i
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  tau0=reference_parameters%fields%tau0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Source=0.d0
  m_i=zone%species(nion)%element%mass
  do i=1,Nx
     do j=1,Nz
        Source=Source+(1.5d0*zone%species(nion)%var(2)%temperature(i,j)*zone%species(nion)%var(2)%density(i,j)&
             +0.5D0*zone%species(nion)%var(2)%Gamma(i,j)**2.D0/zone%species(nion)%var(2)%density(i,j)*m_i&
             -1.5d0*zone%species(nion)%var(1)%temperature(i,j)*zone%species(nion)%var(1)%density(i,j)&
             -0.5D0*zone%species(nion)%var(1)%Gamma(i,j)**2.D0/zone%species(nion)%var(1)%density(i,j)*m_i)&
             /global_variables%dt&
             *zone%metric_coefficients%dvol_pu(i,j)&
             *kb*T0*n0/tau0&
             *(1.d0-zone%masks%chi2(i,j))
     end do
  end do
end subroutine compute_energy_variation
