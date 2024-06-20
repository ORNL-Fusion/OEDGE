subroutine compute_velocity_source(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  real*8 :: alpham,alphap
  real*8 :: alpham1,alphap1
  real*8 :: cs,m_i
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  dt=global_variables%dt
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           zone%species(n)%tridiag%S(i,j)=zone%species(n)%var(1)%Gamma(i,j)/dt
        end do
     end do
  end do
end subroutine compute_velocity_source
