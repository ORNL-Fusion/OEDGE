subroutine step_in_soledge_am()
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
  end do

  penalisation_parameters%eta=1.d-6
  penalisation_parameters%eta2=penalisation_parameters%eta*1.d-2
  global_variables%min_density=1.d-8
  global_variables%min_temperature=1.d-5
  penalisation_parameters%dump_pen=500.
  penalisation_parameters%gain=0.5
  global_variables%Teps=1.d0/reference_parameters%fields%T0eV

  global_variables%dt=1.d-2

  do k=1,global_parameters%N_zones
     call compute_am_vars(zones(k))
     call compute_am_sources(zones(k))
  end do

  !################################################################################
  !#######################    density    ##########################################
  do k=1,global_parameters%N_zones
     call compute_tridiag_matrix_implicit_density(zones(k))
     call add_density_am_source(zones(k))
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
  !#######################   temperature    #######################################
  do k=1,global_parameters%N_zones
     do n_ion = 0,global_parameters%n_ions
        zones(k)%species(n_ion)%var(2)%temperature=1.D0
     end do
  end do
  !################################################################################
  !################################################################################


  do k=1,global_parameters%N_zones
     call check_fields(zones(k))
     call compute_local_residuals(zones(k))
  end do

  call compute_residuals()

  global_variables%tempus=global_variables%tempus&
       +global_variables%dt*reference_parameters%fields%tau0

end subroutine step_in_soledge_am
