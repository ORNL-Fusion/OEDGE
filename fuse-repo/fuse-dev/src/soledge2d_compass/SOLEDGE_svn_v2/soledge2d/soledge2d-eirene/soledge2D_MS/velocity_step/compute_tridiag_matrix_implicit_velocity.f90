subroutine compute_tridiag_matrix_implicit_velocity(zone)
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
           call compute_tridiag_matrix_implicit_velocity_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine compute_tridiag_matrix_implicit_velocity


subroutine compute_tridiag_matrix_implicit_velocity_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: rs0,dt
  rs0=reference_parameters%geometry%rs0
  dt=global_variables%dt
  zone%species(n)%tridiag%a(i,j)=0.D0
  zone%species(n)%tridiag%b(i,j)=1.D0/dt
  zone%species(n)%tridiag%c(i,j)=0.d0
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%var(1)%Gamma(i,j)/dt
end subroutine compute_tridiag_matrix_implicit_velocity_point
