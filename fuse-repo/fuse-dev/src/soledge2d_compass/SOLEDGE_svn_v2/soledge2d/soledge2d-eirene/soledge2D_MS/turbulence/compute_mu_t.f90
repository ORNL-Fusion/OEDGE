subroutine compute_mu_t(zone)
  use all_variables, only : kepsilon_param, reference_parameters
  use MZone
  use Mphysics
  implicit none
  type(TZone), intent(inout) :: zone
  real*8 :: eps
  real*8,allocatable :: rho_L(:,:), c_s(:,:)
!!$  zone%kepsilon(1)%mu_t=max(min(kepsilon_param%Cmu&
!!$       *zone%kepsilon(1)%k*zone%kepsilon(1)%k/(zone%kepsilon(1)%epsilon)*reference_parameters%fields%k0**2/&
!!$       reference_parameters%fields%epsilon0,kepsilon_param%mu_max),kepsilon_param%mu_min)*reference_parameters%fields%epsilon0/&
!!$       reference_parameters%fields%k0**2

!!$  allocate(rho_L(0:zone%mesh%Nx+1,0:zone%mesh%Nz+1))
!!$  ! Larmor radius in m
!!$  rho_L=sqrt((zone%species(0)%var(1)%temperature+zone%species(1)%var(1)%temperature)*reference_parameters%fields%T0eV&
!!$       *zone%species(1)%element%mass*m_u/(eV*zone%species(1)%charge*zone%mesh%B*zone%mesh%B))       
!!$
!!$  zone%kepsilon(1)%mu_t=max(min(kepsilon_param%Cmu*sqrt(zone%kepsilon(1)%k)*&
!!$       rho_L&
!!$       *sqrt(reference_parameters%fields%k0),kepsilon_param%mu_max),kepsilon_param%mu_min)/&
!!$       (reference_parameters%geometry%rs0*reference_parameters%fields%c0)
!!$  deallocate(rho_L)

  allocate(c_s(0:zone%mesh%Nx+1,0:zone%mesh%Nz+1))
  c_s=sqrt((zone%species(0)%var(1)%temperature+zone%species(1)%var(1)%temperature)&
       /zone%species(1)%element%mass)*reference_parameters%fields%c0
  zone%kepsilon(1)%mu_t=kepsilon_param%Cmu*zone%kepsilon(1)%k*&
       zone%mesh%Rgeom/c_s*reference_parameters%fields%k0
  zone%kepsilon(1)%mu_t=min(zone%kepsilon(1)%mu_t,kepsilon_param%mu_max)
!  zone%kepsilon(1)%mu_t=max(zone%kepsilon(1)%mu_t,kepsilon_param%mu_min&
!       *sqrt(zone%species(0)%var(1)%temperature*reference_parameters%fields%T0eV/50.d0))
  zone%kepsilon(1)%mu_t=max(zone%kepsilon(1)%mu_t,kepsilon_param%mu_min)
  zone%kepsilon(1)%mu_t=zone%kepsilon(1)%mu_t/&
       (reference_parameters%geometry%rs0**2/reference_parameters%fields%tau0)
  deallocate(c_s)

end subroutine compute_mu_t
