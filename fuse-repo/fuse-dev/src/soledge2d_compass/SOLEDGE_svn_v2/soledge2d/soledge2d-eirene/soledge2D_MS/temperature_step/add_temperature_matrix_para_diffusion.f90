subroutine add_temperature_matrix_para_diffusion(zone)
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
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call add_temperature_matrix_para_diffusion_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_matrix_para_diffusion


subroutine add_temperature_matrix_para_diffusion_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dvol,ds_east,ds_west
  real*8 :: rs0,R0,Teps
  real*8 :: kappa, kappa_m1, kappa_p1
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Teps=global_variables%Teps
  kappa=zone%species(n)%transport_para%kappa(i,j)/&
       zone%species(n)%var(1)%log_Lambda(i,j)
  kappa_p1=zone%species(n)%transport_para%kappa(i,j+1)/&
       zone%species(n)%var(1)%log_Lambda(i,j+1)
  kappa_m1=zone%species(n)%transport_para%kappa(i,j-1)/&
       zone%species(n)%var(1)%log_Lambda(i,j-1)
  dvol=zone%metric_coefficients%dvol_dd(i,j)
  ds_west=zone%metric_coefficients%ds_west_dd(i,j)
  ds_east=zone%metric_coefficients%ds_east_dd(i,j)
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
       -(kappa_m1*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0&
       +kappa*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +(kappa_m1*(max(zone%species(n)%var(1)%temperature(i,j-1),Teps))**2.5d0&
       +kappa*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       +(kappa*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +kappa_p1*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)&
       -(kappa*(max(zone%species(n)%var(1)%temperature(i,j),Teps))**2.5d0&
       +kappa_p1*(max(zone%species(n)%var(1)%temperature(i,j+1),Teps))**2.5d0)*0.5d0&
       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
end subroutine add_temperature_matrix_para_diffusion_point
