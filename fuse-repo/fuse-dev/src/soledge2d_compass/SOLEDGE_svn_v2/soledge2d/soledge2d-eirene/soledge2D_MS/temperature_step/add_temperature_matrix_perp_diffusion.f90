subroutine add_temperature_matrix_perp_diffusion(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call add_temperature_matrix_perp_diffusion_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_matrix_perp_diffusion


subroutine add_temperature_matrix_perp_diffusion_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dvol,ds_east,ds_west
  real*8 :: rs0,R0,Teps
  rs0=reference_parameters%geometry%rs0
  R0=reference_parameters%geometry%R0
  Teps=global_variables%Teps
  dvol=zone%metric_coefficients%dvol_dd(i,j)
  ds_west=zone%metric_coefficients%ds_west_dd(i,j)
  ds_east=zone%metric_coefficients%ds_east_dd(i,j)
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
       -zone%species(n)%implicit_coefs%west_temperature(i,j)&
       *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       *(zone%species(n)%var(1)%density(i,j-1)+zone%species(n)%var(1)%density(i,j))*0.5d0
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +zone%species(n)%implicit_coefs%west_temperature(i,j)&
       *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       *(zone%species(n)%var(1)%density(i,j-1)+zone%species(n)%var(1)%density(i,j))*0.5d0&
       +zone%species(n)%implicit_coefs%east_temperature(i,j)&
       *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       *(zone%species(n)%var(1)%density(i,j+1)+zone%species(n)%var(1)%density(i,j))*0.5d0
  zone%species(n)%tridiag%c(i,j)= zone%species(n)%tridiag%c(i,j)&
       -zone%species(n)%implicit_coefs%east_temperature(i,j)&
       *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       *(zone%species(n)%var(1)%density(i,j+1)+zone%species(n)%var(1)%density(i,j))*0.5d0
end subroutine add_temperature_matrix_perp_diffusion_point
