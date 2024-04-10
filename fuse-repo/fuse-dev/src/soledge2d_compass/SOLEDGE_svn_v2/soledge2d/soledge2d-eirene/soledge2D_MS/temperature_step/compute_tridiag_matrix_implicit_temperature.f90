subroutine compute_tridiag_matrix_implicit_temperature(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4 :: n
  integer*4 :: i,j
  integer*4 :: Nx,Nz
  real*8 :: dt
  dt=global_variables%dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     do i=1,Nx
        do j=1,Nz
           call compute_tridiag_matrix_implicit_temperature_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine compute_tridiag_matrix_implicit_temperature


subroutine compute_tridiag_matrix_implicit_temperature_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone),intent(inout) :: zone
  integer*4,intent(in) :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt, m_i
  m_i=zone%species(n)%element%mass
  dt=global_variables%dt
  zone%species(n)%tridiag%a(i,j)=0.d0
  zone%species(n)%tridiag%b(i,j)=1.5d0*zone%species(n)%var(2)%density(i,j)/dt
  zone%species(n)%tridiag%c(i,j)=0.d0
  zone%species(n)%tridiag%S(i,j)=1.5d0*zone%species(n)%var(1)%density(i,j)*zone%species(n)%var(1)%temperature(i,j)/dt&
       -0.5d0*m_i*(zone%species(n)%var(2)%Gamma(i,j)**2.D0/zone%species(n)%var(2)%density(i,j)&
       -zone%species(n)%var(1)%Gamma(i,j)**2.D0/zone%species(n)%var(1)%density(i,j))/dt
end subroutine compute_tridiag_matrix_implicit_temperature_point
