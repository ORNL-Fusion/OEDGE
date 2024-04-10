subroutine add_velocity_Eparallel_term(zone)
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
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           call add_velocity_Eparallel_term_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_velocity_Eparallel_term


subroutine add_velocity_Eparallel_term_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  real*8 :: alpham,alphap
  real*8 :: alpham1,alphap1
  real*8 :: cs,m_i
  m_i=zone%species(n)%element%mass
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%species(n)%charge*zone%species(n)%var(1)%density(i,j)&
       *zone%electric_fields(1)%E(i,j)/m_i
  zone%species(n)%sources%volumic_sources_G(i,j)=zone%species(n)%sources%volumic_sources_G(i,j)&
       +zone%species(n)%charge*zone%species(n)%var(1)%density(i,j)&
       *zone%electric_fields(1)%E(i,j)/m_i
end subroutine add_velocity_Eparallel_term_point
