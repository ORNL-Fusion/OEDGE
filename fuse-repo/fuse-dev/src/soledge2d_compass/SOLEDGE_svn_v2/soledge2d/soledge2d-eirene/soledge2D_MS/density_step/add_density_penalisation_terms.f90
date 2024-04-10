subroutine add_density_penalisation_terms(zone)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  real*8 :: dt
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  dt=global_variables%dt
  do i=1,Nx
     do j=1,Nz
        do n=1,global_parameters%N_ions
           call add_density_penalisation_terms_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_density_penalisation_terms


subroutine add_density_penalisation_terms_point(zone,i,j,n)
  use all_variables, only : global_parameters, reference_parameters,&
       global_variables, penalisation_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  real*8 :: dt
  dt=global_variables%dt
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%masks%chi4(i,j)/penalisation_parameters%eta2&
       *global_variables%min_density*10.d0
  zone%species(n)%tridiag%b(i,j)=zone%species(n)%tridiag%b(i,j)&
       +zone%masks%chi4(i,j)/penalisation_parameters%eta2
end subroutine add_density_penalisation_terms_point
