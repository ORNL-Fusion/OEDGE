subroutine add_velocity_am_source(zone)
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
           call add_velocity_am_source_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_velocity_am_source


subroutine add_velocity_am_source_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dt
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%species(n)%sources%SG_am(i,j)
  zone%species(n)%sources%volumic_sources_G(i,j)=zone%species(n)%sources%volumic_sources_G(i,j)&
       +zone%species(n)%sources%SG_am(i,j)
end subroutine add_velocity_am_source_point
