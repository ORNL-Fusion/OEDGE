subroutine allocate_species_arrays(zone)
  use all_variables, only : global_parameters, global_variables
  use MZone 
  implicit none
  type(TZone),intent(inout) :: zone
  integer*4 :: Nx,Nz
  integer*4 :: n
  Nx=zone%mesh%Nx
  Nz=zone%mesh%Nz
  do n=0,global_parameters%N_ions
     allocate(zone%species(n)%transport_perp%D_p(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%nu_p(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%chi_p(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%zeta_p(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%D_t(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%nu_t(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%chi_t(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%zeta_t(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_perp%v_pinch(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_para%kappa(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_para%kappa0(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_para%nu(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%transport_para%nu0(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%density(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%velocity(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%temperature(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%Gamma(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%Mach(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(2)%density(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(2)%Gamma(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(2)%velocity(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(2)%Mach(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(2)%temperature(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%var(1)%log_Lambda(0:Nx+1,0:Nz+1))
     allocate(zone%species(n)%penalisation_memories%alpham(1:Nx,1:Nz))
     allocate(zone%species(n)%penalisation_memories%alphap(1:Nx,1:Nz))
     !init to avoid NaN
     zone%species(n)%var(1)%density=1.d-15
     zone%species(n)%var(1)%Gamma=1.d-15
     zone%species(n)%var(1)%temperature=1.d-15
     zone%species(n)%var(2)%density=1.d-15
     zone%species(n)%var(2)%Gamma=1.d-15
     zone%species(n)%var(2)%temperature=1.d-15
     zone%species(n)%var(1)%log_Lambda=1.d-15
     zone%species(n)%transport_perp%v_pinch=0.d0
     allocate(zone%species(n)%coupling_terms%tau(1:Nx,1:Nz,0:global_parameters%N_ions))
     allocate(zone%species(n)%coupling_terms%R(1:Nx,1:Nz,0:global_parameters%N_ions))
     allocate(zone%species(n)%coupling_terms%Q(1:Nx,1:Nz,0:global_parameters%N_ions))
     allocate(zone%species(n)%coupling_terms%qu(1:Nx,0:Nz+1,0:global_parameters%N_ions))
  end do
  allocate(zone%electric_fields(1)%E(1:Nx,1:Nz))
  allocate(zone%electric_fields(1)%Etheta(1:Nx,1:Nz))
  allocate(zone%electric_fields(1)%Epsi(1:Nx,1:Nz))
  !shared fields
  allocate(zone%shared_fields%log_Lambda(0:Nx+1,0:Nz+1))
  allocate(zone%shared_fields%Zeff(0:Nx+1,0:Nz+1))
end subroutine allocate_species_arrays
