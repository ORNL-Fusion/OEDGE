subroutine allocate_tridiagonal_systems(zone)
  use all_variables, only : global_parameters
  use MZone 
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: n
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     allocate(zone%species(n)%tridiag%a(1:Nx,1:Nz))
     allocate(zone%species(n)%tridiag%b(1:Nx,1:Nz))
     allocate(zone%species(n)%tridiag%c(1:Nx,1:Nz))
     allocate(zone%species(n)%tridiag%S(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%east_density(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%west_density(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%east_velocity(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%west_velocity(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%east_temperature(1:Nx,1:Nz))
     allocate(zone%species(n)%implicit_coefs%west_temperature(1:Nx,1:Nz))
  end do
end subroutine allocate_tridiagonal_systems
