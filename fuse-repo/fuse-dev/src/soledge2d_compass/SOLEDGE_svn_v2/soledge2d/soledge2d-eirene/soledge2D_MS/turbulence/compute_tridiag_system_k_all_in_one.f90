subroutine compute_tridiag_system_k_all_in_one(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters, kepsilon_param
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  real*8 :: dvol,ds_east,ds_west
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: rs0,dt
  real*8 :: etaT
  real*8 :: grad1, grad2
  rs0=reference_parameters%geometry%rs0
  dt=global_variables%dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  etaT=penalisation_parameters%eta**2
  do i=1,Nx
     do j=1,Nz
        dvol=zone%metric_coefficients%dvol_dd(i,j)
        ds_west=zone%metric_coefficients%ds_west_dd(i,j)
        ds_east=zone%metric_coefficients%ds_east_dd(i,j)
        zone%species(0)%tridiag%a(i,j)=0.D0
        zone%species(0)%tridiag%b(i,j)=zone%species(0)%var(2)%density(i,j)/dt
        zone%species(0)%tridiag%c(i,j)=0.D0
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%var(1)%density(i,j)*zone%kepsilon(1)%k(i,j)/dt
        !perp
        zone%species(0)%tridiag%a(i,j)=zone%species(0)%tridiag%a(i,j)-zone%kepsilon(1)%implicit_coefs%west_k(i,j)&
             *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
             *(zone%species(0)%var(1)%density(i,j-1)+zone%species(0)%var(1)%density(i,j))*0.5d0
        zone%species(0)%tridiag%b(i,j)=zone%species(0)%tridiag%b(i,j)+zone%kepsilon(1)%implicit_coefs%west_k(i,j)&
             *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
             *(zone%species(0)%var(1)%density(i,j-1)+zone%species(0)%var(1)%density(i,j))*0.5d0&
             +zone%kepsilon(1)%implicit_coefs%east_k(i,j)&
             *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
             *(zone%species(0)%var(1)%density(i,j+1)+zone%species(0)%var(1)%density(i,j))*0.5d0
        zone%species(0)%tridiag%c(i,j)=zone%species(0)%tridiag%c(i,j)-zone%kepsilon(1)%implicit_coefs%east_k(i,j)&
             *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
             *(zone%species(0)%var(1)%density(i,j+1)+zone%species(0)%var(1)%density(i,j))*0.5d0
        !explicit source
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
             +zone%kepsilon(1)%Sk(i,j)
        !penalisation
        zone%species(0)%tridiag%b(i,j)=zone%species(0)%tridiag%b(i,j)&
             +zone%masks%chi2(i,j)/penalisation_parameters%eta
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
             +zone%masks%chi2(i,j)/penalisation_parameters%eta*1.d-10        
        !k-epsilon source
!        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
!             -zone%species(0)%var(1)%density(i,j)*zone%kepsilon(1)%epsilon(i,j)
        !interchange
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
             +zone%kepsilon(1)%interchange(i,j)*max(zone%kepsilon(1)%k(i,j),1.d-5)&
             *zone%species(0)%var(1)%density(i,j)
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
             +zone%species(1)%sources%Sn_n(i,j)*zone%kepsilon(1)%k(i,j)
        !non linear saturation
        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
             -kepsilon_param%deltaOmega*reference_parameters%fields%tau0*reference_parameters%fields%k0&
             *zone%kepsilon(1)%k(i,j)*zone%kepsilon(1)%k(i,j)&
             *zone%species(0)%var(1)%density(i,j)
!!$        ! background source
!!$        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
!!$             +zone%species(0)%var(1)%density(i,j)*2.d0*zone%species(0)%var(1)%temperature(i,j)
        !shear --> moved to epsilon eqn
!!$        zone%species(0)%tridiag%S(i,j)=zone%species(0)%tridiag%S(i,j)&
!!$             -kepsilon_param%tauV*zone%kepsilon(1)%UE_shear(i,j)&
!!$             *zone%kepsilon(1)%k(i,j)*zone%species(0)%var(1)%density(i,j)&
!!$             *reference_parameters%fields%tau0*reference_parameters%fields%c0**2&
!!$             /reference_parameters%geometry%rs0**2
     end do
  end do
end subroutine compute_tridiag_system_k_all_in_one
