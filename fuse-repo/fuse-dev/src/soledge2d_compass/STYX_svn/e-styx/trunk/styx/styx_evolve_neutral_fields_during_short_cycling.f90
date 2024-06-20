subroutine styx_evolve_neutral_fields_during_short_cycling
  use all_variables, only : global_parameters, global_variables, reference_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_comusr
  use styx2eirene
  implicit none
  integer :: isp,ipls
  real(dp) :: Sat(natmi,Neir_cells)
  real(dp) :: Smo(nmoli,Neir_cells)
  real(dp) :: Sti(nioni,Neir_cells)
  real(dp) :: Dt_refresh

  ! time step between refreshing sources

  Dt_refresh = ns_refresh*global_variables%dt*reference_parameters%fields%tau0

  ! sources resolved by atom species
  Sat=0._dp 
  do ipls=1,global_parameters%n_ions
    isp=global_parameters%ions_list(ipls,1)
    Sat(isp,:)=Sat(isp,:)+Sn_at(:,ipls)
  enddo

  Smo=0._dp
  Sti=0._dp
  do ipls=1,global_parameters%n_ions
    Smo(1,:)=Smo(1,:)+Sn_mol(:,ipls)
    Sti(1,:)=Sti(1,:)+Sn_tion(:,ipls)
  enddo

  atom_density = atom_density-Sat*Dt_refresh
  mol_density  = mol_density-Smo*Dt_refresh
  tion_density = tion_density-Sti*Dt_refresh
  
end subroutine styx_evolve_neutral_fields_during_short_cycling
