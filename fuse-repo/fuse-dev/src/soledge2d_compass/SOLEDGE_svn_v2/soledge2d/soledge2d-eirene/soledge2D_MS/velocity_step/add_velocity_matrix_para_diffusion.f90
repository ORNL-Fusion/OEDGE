subroutine add_velocity_matrix_para_diffusion(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n
  real*8 :: dvol,ds_east,ds_west
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: rs0,R0,Teps
  real*8,allocatable :: nu(:,:)
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Teps=global_variables%Teps
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  allocate(nu(0:Nx+1,0:Nz+1))
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call add_velocity_matrix_para_diffusion_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_velocity_matrix_para_diffusion


subroutine add_velocity_matrix_para_diffusion_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dvol,ds_east,ds_west
  real*8 :: rs0,R0,Teps
  real*8 :: nu, nu_m1, nu_p1
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Teps=global_variables%Teps
  nu=zone%species(n)%transport_para%nu(i,j)/&
       zone%species(n)%var(1)%log_Lambda(i,j)
  nu_m1=zone%species(n)%transport_para%nu(i,j-1)/&
       zone%species(n)%var(1)%log_Lambda(i,j-1)
  nu_p1=zone%species(n)%transport_para%nu(i,j+1)/&
       zone%species(n)%var(1)%log_Lambda(i,j+1)
  dvol=zone%metric_coefficients%dvol_dd(i,j)
  ds_west=zone%metric_coefficients%ds_west_dd(i,j)
  ds_east=zone%metric_coefficients%ds_east_dd(i,j)
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
       -(nu_m1*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0&
       +nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       /zone%species(n)%var(1)%density(i,j-1)*(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +(nu_m1*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0&
       +nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       /zone%species(n)%var(1)%density(i,j)*(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))&
       +(nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +nu_p1*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       /zone%species(n)%var(1)%density(i,j)*(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)&
       -(nu*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +nu_p1*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       /zone%species(n)%var(1)%density(i,j+1)*(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
end subroutine add_velocity_matrix_para_diffusion_point
