subroutine allocate_integrals()
  use all_variables, only : global_parameters, globals
  implicit none
  allocate(globals%flux_tot_out_ac(0:global_parameters%N_ions))
  allocate(globals%flux_tot_in_ac(0:global_parameters%N_ions))
  allocate(globals%flux_totE_out_ac(0:global_parameters%N_ions))
  allocate(globals%flux_totE_in_ac(0:global_parameters%N_ions))
  allocate(globals%source_n(0:global_parameters%N_ions))
  allocate(globals%source_E(0:global_parameters%N_ions))
  allocate(globals%variation_E(0:global_parameters%N_ions))
  allocate(globals%stored_E(0:global_parameters%N_ions))
  allocate(globals%variation_N(0:global_parameters%N_ions))
  allocate(globals%stored_N(0:global_parameters%N_ions))
  allocate(globals%source_ionz_tot(1:global_parameters%N_species))
  allocate(globals%source_E_ionz(0:global_parameters%N_species))
  allocate(globals%total_radiation(0:global_parameters%N_ions))
  allocate(globals%tot_N(0:global_parameters%N_ions))
  allocate(globals%tot_E(0:global_parameters%N_ions))
end subroutine allocate_integrals
