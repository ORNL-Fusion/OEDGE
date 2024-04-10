subroutine add_temperature_coupling_terms(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n,m
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           call add_temperature_coupling_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_coupling_terms


subroutine add_temperature_coupling_terms_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  integer*4 :: m
  do m=0,global_parameters%N_ions
     zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
          +zone%species(n)%coupling_terms%Q(i,j,m)&
          +zone%species(n)%coupling_terms%R(i,j,m)*zone%species(n)%var(1)%velocity(i,j)
     zone%species(n)%sources%volumic_sources_E(i,j)=zone%species(n)%sources%volumic_sources_E(i,j)&
          +zone%species(n)%coupling_terms%Q(i,j,m)&
          +zone%species(n)%coupling_terms%R(i,j,m)*zone%species(n)%var(1)%velocity(i,j)
  end do
end subroutine add_temperature_coupling_terms_point
