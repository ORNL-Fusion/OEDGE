subroutine add_temperature_explicit_source(zone)
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
           call add_temperature_explicit_source_point(zone,i,j,n)
        end do
     end do
  end do
end subroutine add_temperature_explicit_source


subroutine add_temperature_explicit_source_point(zone,i,j,n)
  use all_variables, only : global_parameters
  use Mzone
  use Mphysics
  implicit none
  type(Tzone) :: zone
  integer*4,intent(in) :: i,j,n
  zone%species(n)%tridiag%S(i,j)=zone%species(n)%tridiag%S(i,j)&
       +zone%species(n)%sources%SE(i,j)
end subroutine add_temperature_explicit_source_point
