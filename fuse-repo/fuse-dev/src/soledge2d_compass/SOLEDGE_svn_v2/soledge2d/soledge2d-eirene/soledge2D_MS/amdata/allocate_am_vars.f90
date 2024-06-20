subroutine allocate_am_vars()
  use all_variables, only : global_parameters, zones
  implicit none
  integer*4 :: n,k
  integer*4 :: Nx,Nz
  do k=1,global_parameters%N_zones
     Nx=zones(k)%mesh%Nx
     Nz=zones(k)%mesh%Nz
     do n=1,global_parameters%N_ions
        allocate(zones(k)%species(n)%am_vars%ionization_rate_coefficient(1:Nx,1:Nz))
        allocate(zones(k)%species(n)%am_vars%recombination_rate_coefficient(1:Nx,1:Nz))
        allocate(zones(k)%species(n)%am_vars%radiation_function_excitation(1:Nx,1:Nz))
        allocate(zones(k)%species(n)%am_vars%radiation_function_recombination_bremsstrahlung(1:Nx,1:Nz))
     end do
  end do
end subroutine allocate_am_vars
