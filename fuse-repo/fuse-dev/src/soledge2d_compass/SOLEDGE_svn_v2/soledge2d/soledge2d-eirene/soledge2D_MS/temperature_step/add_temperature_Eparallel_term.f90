subroutine add_temperature_Eparallel_term(zone)
  use all_variables, only : global_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4 :: i,j,n
  integer*4 :: Nx,Nz
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do i=1,Nx
     do j=1,Nz
        do n=0,global_parameters%N_ions
           call add_temperature_Eparallel_term_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_Eparallel_term


subroutine add_temperature_Eparallel_term_point(zone,i,j,n)
  use all_variables, only : global_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%species(n)%charge*zone%species(n)%var(1)%Gamma(i,j)&
       *zone%electric_fields(1)%E(i,j)
  zone%species(n)%sources%volumic_sources_E(i,j)=zone%species(n)%sources%volumic_sources_E(i,j)&
       +zone%species(n)%charge*zone%species(n)%var(1)%Gamma(i,j)&
       *zone%electric_fields(1)%E(i,j)
end subroutine add_temperature_Eparallel_term_point
