subroutine allocate_species_properties(N_species)
 use all_variables, only : global_parameters, boundary_conditions,&
      transport_parameters, element_variables
 implicit none
 integer*4, intent(in) :: N_species

 allocate(global_parameters%element_list(0:N_species))
 allocate(boundary_conditions%BCn_model(1:N_species))
 allocate(boundary_conditions%BCT_model(1:N_species))
 allocate(boundary_conditions%BCn(1:N_species))
 allocate(boundary_conditions%BCTi(1:N_species))
 allocate(transport_parameters%Dn_t(1:N_species)) 
 allocate(transport_parameters%Dn_p(1:N_species)) 
 allocate(transport_parameters%nu_t(1:N_species)) 
 allocate(transport_parameters%nu_p(1:N_species)) 
 allocate(transport_parameters%chii_t(1:N_species)) 
 allocate(transport_parameters%chii_p(1:N_species)) 
 allocate(transport_parameters%v_pinch(1:N_species)) 
 allocate(transport_parameters%flux_limiter(0:N_species)) 
 allocate(transport_parameters%flux_limiter_nu(1:N_species)) 
 allocate(element_variables(1:N_species))

end subroutine allocate_species_properties
