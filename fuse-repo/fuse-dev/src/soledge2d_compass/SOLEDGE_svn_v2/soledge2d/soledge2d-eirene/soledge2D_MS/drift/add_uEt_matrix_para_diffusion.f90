subroutine add_uEt_matrix_para_diffusion(zone)
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
           call add_uEt_matrix_para_diffusion_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_uEt_matrix_para_diffusion


subroutine add_uEt_matrix_para_diffusion_point(zone,i,j,n)
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
  nu=global_variables%eta_para_smooth
  nu_m1=nu
  nu_p1=nu
  dvol=zone%metric_coefficients%dvol_dd(i,j)
  ds_west=zone%metric_coefficients%ds_west_dd(i,j)
  ds_east=zone%metric_coefficients%ds_east_dd(i,j)
  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
       -(nu_m1+nu)*0.5d0&
       *(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j-1)))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       *(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +(nu_m1+nu)*0.5d0&
       *(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j-1)))*0.5d0&
       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
       *(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))&
       +(nu+nu_p1)*0.5d0&
       *(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j+1)))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       *(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)&
       -(nu+nu_p1)*0.5d0&
       *(sqrt(zone%metric_coefficients%ctt(i,j))+sqrt(zone%metric_coefficients%ctt(i,j+1)))*0.5d0&
       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
       *(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
!!$  zone%species(n)%tridiag%a(i,j)=zone%species(n)%tridiag%a(i,j)&
!!$       -(nu_m1+nu)*0.5d0&
!!$       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
!!$       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
!!$       *(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))
!!$  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
!!$       +(nu_m1+nu)*0.5d0&
!!$       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j-1))*0.5d0&
!!$       *ds_west/dvol*zone%metric_coefficients%sinepitch_west(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j)-zone%mesh%z(i,j-1))&
!!$       *(1.d0-zone%masks%chi4(i,j-1))*(1.d0-zone%masks%chi4(i,j))&
!!$       +(nu+nu_p1)*0.5d0&
!!$       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
!!$       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
!!$       *(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
!!$  zone%species(n)%tridiag%c(i,j)=zone%species(n)%tridiag%c(i,j)&
!!$       -(nu+nu_p1)*0.5d0&
!!$       *(zone%metric_coefficients%G(i,j)+zone%metric_coefficients%G(i,j+1))*0.5d0&
!!$       *ds_east/dvol*zone%metric_coefficients%sinepitch_east(i,j)*(2.D0*pi*R0/rs0)/(zone%mesh%z(i,j+1)-zone%mesh%z(i,j))&
!!$       *(1.d0-zone%masks%chi4(i,j+1))*(1.d0-zone%masks%chi4(i,j))
end subroutine add_uEt_matrix_para_diffusion_point
