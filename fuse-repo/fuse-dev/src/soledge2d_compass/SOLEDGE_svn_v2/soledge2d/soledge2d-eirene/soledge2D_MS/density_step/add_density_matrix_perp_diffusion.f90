subroutine add_density_matrix_perp_diffusion(zone)
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
  real*8 :: rs0,dt
  rs0=reference_parameters%geometry%rs0
  dt=global_variables%dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           call add_density_matrix_perp_diffusion_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_density_matrix_perp_diffusion


subroutine add_density_matrix_perp_diffusion_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dvol,ds_east,ds_west
  integer*4 :: Nx,Nz
  real*8 :: rs0,dt
  rs0=reference_parameters%geometry%rs0
  dt=global_variables%dt
  dvol=zone%metric_coefficients%dvol_dd(i,j)
  ds_west=zone%metric_coefficients%ds_west_dd(i,j)
  ds_east=zone%metric_coefficients%ds_east_dd(i,j)
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)-zone%species(n)%implicit_coefs%west_density(i,j)&
       *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)+zone%species(n)%implicit_coefs%west_density(i,j)&
       *ds_west/dvol/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       +zone%species(n)%implicit_coefs%east_density(i,j)&
       *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)-zone%species(n)%implicit_coefs%east_density(i,j)&
       *ds_east/dvol/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))
end subroutine add_density_matrix_perp_diffusion_point
