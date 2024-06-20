subroutine compute_particle_source(zone,nion,Source)
  use MZone
  use all_variables, only : reference_parameters, flags
  implicit none
  Type(TZone),intent(in) :: zone
  integer*4,intent(in) :: nion
  real*8,intent(out) :: Source
  integer*4 :: i,j,Nx,Nz
  real*8 :: rs0,R0,n0,c0,T0,tau0
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  n0=reference_parameters%fields%n0
  c0=reference_parameters%fields%c0
  T0=reference_parameters%fields%T0
  tau0=reference_parameters%fields%tau0
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  Source=0.d0
  do i=1,Nx
     do j=1,Nz
        Source=Source+zone%species(nion)%sources%volumic_sources_n(i,j)&
             *zone%metric_coefficients%dvol_pu(i,j)&
             *n0/tau0&
             *(1.d0-zone%masks%chi2(i,j))
     end do
  end do
end subroutine compute_particle_source
