subroutine step_in_soledge_advection()
  use all_variables
  use test_var
  use Mdefinitions
  implicit none
  integer*4 :: k,npsi,n_ion
  integer*4 :: j,i

  do k=1,global_parameters%N_zones
     call MD_broadcast_plasma(zones(k),STEP_OLD)
  end do

  do k=1,global_parameters%N_zones
     call MD_broadcast_corners(zones(k))
  end do

  do k=1,global_parameters%N_zones
     call compute_velocity_and_Mach(zones(k),STEP_OLD)
     call compute_electron_density_and_velocity(zones(k),STEP_OLD)
     call compute_coulomb_logarithm(zones(k),STEP_OLD)
     call reset_source_terms(zones(k))
     call compute_explicit_parallel_advection(zones(k))
  end do

  global_variables%dt=1.d-5
  penalisation_parameters%eta=1.d-6
  penalisation_parameters%eta2=penalisation_parameters%eta*1.d-2
  global_variables%min_density=1.d-8
  global_variables%min_temperature=1.d-5
  penalisation_parameters%dump_pen=500.
  penalisation_parameters%gain=0.5
  global_variables%Teps=1.d0/reference_parameters%fields%T0eV

  !################################################################################
  !#######################    density    ##########################################
  do k=1,global_parameters%N_zones
     call compute_tridiag_matrix_implicit_density(zones(k))
     call add_density_explicit_source(zones(k))
     call add_density_penalisation_terms(zones(k))
     call add_test_source(zones(k),DENSITY_FIELD)
  end do

  do npsi=1,global_parameters%N_psi_surfaces
     do n_ion = 1,global_parameters%n_ions
        call concat_tridiag_system(flux_surfaces(npsi),n_ion)
        call apply_density_boundary_conditions(flux_surfaces(npsi),n_ion)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_density_BC(flux_surfaces(npsi))
        call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,DENSITY_FIELD)
     end do
  end do
  !################################################################################
  !################################################################################

  !################################################################################
  !#######################    velocity    #########################################
  do k=1,global_parameters%N_zones
     call compute_electric_field(zones(k))
     call compute_tridiag_matrix_implicit_velocity(zones(k))
     call compute_velocity_source_implicit_correction(zones(k))
     call add_velocity_explicit_source(zones(k))
     call add_velocity_divb_term(zones(k))
     call add_velocity_Eparallel_term(zones(k))
     call add_velocity_penalisation_terms(zones(k))
     call add_test_source(zones(k),VELOCITY_FIELD)
  end do

  do npsi=1,global_parameters%N_psi_surfaces
     do n_ion = 1,global_parameters%n_ions
        call concat_tridiag_system(flux_surfaces(npsi),n_ion)
        call apply_velocity_boundary_conditions(flux_surfaces(npsi),n_ion)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_velocity_BC(flux_surfaces(npsi),n_ion)
        call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,VELOCITY_FIELD)
     end do
  end do
  !################################################################################
  !################################################################################

  do k=1,global_parameters%N_zones
     call clean_solution1(zones(k))
     call compute_velocity_and_Mach(zones(k),STEP_NEW)
     call compute_electron_density_and_velocity(zones(k),STEP_NEW)
  end do

  !################################################################################
  !#######################    temperature   #######################################
  do k=1,global_parameters%N_zones
     call compute_tridiag_matrix_implicit_temperature(zones(k))
     call compute_temperature_source_implicit_correction(zones(k))
     call add_temperature_penalisation_terms(zones(k))
     call add_temperature_explicit_source(zones(k))
     call add_temperature_Eparallel_term(zones(k))
     call add_test_source(zones(k),TEMPERATURE_FIELD)
  end do

  do npsi=1,global_parameters%N_psi_surfaces
     do n_ion = 0,global_parameters%n_ions
        call concat_tridiag_system(flux_surfaces(npsi),n_ion)
        call apply_temperature_boundary_conditions(flux_surfaces(npsi),n_ion)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_temperature_BC(flux_surfaces(npsi),n_ion)
        call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,TEMPERATURE_FIELD)
     end do
  end do
  !################################################################################
  !################################################################################

end subroutine step_in_soledge_advection
