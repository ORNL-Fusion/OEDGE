subroutine step_in_soledge(n_ite)
#include "compile_opt.inc"
  use all_variables
  use Mdefinitions
  implicit none
  integer*4,intent(in) :: n_ite
  integer*4 :: k,npsi,n_ion
  integer*4 :: j,i
  real*8 :: new_dt, dt_var

  call compute_ions_core_flux()
  call compute_ballooning_weighted_Score()
  !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(k,npsi,n_ion)
  !$OMP DO
  do k=1,global_parameters%N_zones
     call MD_broadcast_plasma(zones(k),STEP_OLD)
     call compute_velocity_and_Mach(zones(k),STEP_OLD)
     if(flags%turbulence_model.eq.1) then
        call MD_broadcast_kepsilon(zones(k),STEP_OLD)
     end if
  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP DO
  do k=1,global_parameters%N_zones
     call MD_broadcast_corners(zones(k),STEP_OLD)
  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP MASTER
  if(flags%turbulence_model.eq.1) then
     call compute_implicit_tridiag_coefficients_ke()
  end if
  !$OMP END MASTER
  !$OMP BARRIER


#if VORTICITY_PASTIX == 1
  if(n_ite.eq.1) then ! init old pi for vorticity
     !$OMP DO	
     do k=1,global_parameters%N_zones
        call store_pi_old(zones(k))
     end do
     !$OMP END DO
  end if
#endif

  !$OMP MASTER
  If(drift_flags%solve_phi.and.(drift_flags%use_ExB_radial.or.drift_flags%use_ExB_poloidal)) then
     call smooth_phi()
  end If
  !$OMP END MASTER
  !$OMP BARRIER

  !$OMP DO	
  do k=1,global_parameters%N_zones
!     call compute_pinch_velocity(zones(k))
#if VORTICITY_PASTIX == 1
     If(drift_flags%solve_phi) then
        call MD_broadcast_corners_phi2(zones(k),STEP_OLD)
     end If
     If(drift_flags%solve_phi) then
        call compute_perp_velocity(zones(k))					!### Leybros modif ###
     end If
#endif
  end do
  !$OMP END DO
  !$OMP BARRIER


  !$OMP MASTER
  call MD_broadcast_vpinch()
#if VORTICITY_PASTIX == 1
  if(drift_flags%solve_drift) then
     call Broadcast_drift2()
  end if
#endif
  !$OMP END MASTER
  !$OMP BARRIER 

  !################################################################################
  !################################################################################

  !$OMP DO
  do k=1,global_parameters%N_zones
     if(drift_flags%solve_drift) then
        call correct_drifts_in_pen_cells(zones(k))
     end if
     call limit_drift_velocities(zones(k))
     call compute_electron_density_and_velocity(zones(k),STEP_OLD)
     if(flags%turbulence_model.eq.1) then
        call MD_broadcast_corners_ke(zones(k),STEP_OLD)
        call compute_mu_t(zones(k))
        call compute_vtheta_shear(zones(k))
        call compute_interchange(zones(k))
        call recompute_diffusion_coefficients(zones(k))
     end if
     call compute_coulomb_logarithm(zones(k),STEP_OLD)
     call compute_coupling_terms(zones(k))
     call flux_limiter_nu(zones(k))
     call flux_limiter(zones(k))
#if VORTICITY_PASTIX == 1
     call MD_broadcast_vorticity(zones(k),STEP_OLD)
     call MD_broadcast_corners_W(zones(k),STEP_OLD)
     call vorticity_continuity_penalisation(zones(k),STEP_OLD)
#endif
     call reset_source_terms(zones(k))
     call compute_explicit_diffusive_perp_source_terms(zones(k))
     call compute_explicit_vpinch_perp_source_terms(zones(k))
     call compute_electron_perp_particle_flux(zones(k))
     call compute_explicit_pseudo_convective_perp_source_termsI(zones(k)) !for ions
     call compute_explicit_pseudo_convective_perp_source_termsE(zones(k)) !for electrons
     if(flags%turbulence_model.eq.1) then
        call compute_explicit_pseudo_convective_perp_source_terms_ke(zones(k)) !for k-epsilon
     end if
     if(drift_flags%use_ExB_radial.or.drift_flags%use_gradB_radial) then 
        call compute_explicit_perpendicular_advection(zones(k))
     end if
     call compute_explicit_parallel_advection(zones(k))
     call compute_pdivv_term(zones(k))
     call compute_collisional_velocity_heat_flux_divergence(zones(k))
     if(flags%neutral_model.eq.2) then
        call compute_neutral_sources(zones(k))
     end if
  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP DO
  do k=1,global_parameters%N_zones
     call compute_dt(zones(k),global_variables%dt_table(k))
  end do
  !$OMP END DO
  !$OMP BARRIER

  !$OMP MASTER
  if(flags%turbulence_model.eq.1) then
     call flux_surface_average_shear()
  end if
  if(flags%neutral_model.eq.2) then
     !fixed dt to avoid refilling the neutrals matrix
     new_dt = minval(global_variables%dt_table)*global_parameters%CFL
     dt_var = global_variables%dt/new_dt
     if((dt_var.lt.0.8d0).or.(dt_var.gt.1.25d0)) then
        ! new dt (else nothing changes, dt is not updated)
        global_variables%dt=new_dt
        call recompute_neutral_matrix()
     end if
  else
     global_variables%dt=minval(global_variables%dt_table)*global_parameters%CFL
  end if
  call MD_broadcast_flux_limiters()
  !$OMP END MASTER
  !$OMP BARRIER 

  !$OMP DO
  do k=1,global_parameters%N_zones
     call compute_am_vars(zones(k))
     call compute_am_sources(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER


  !################################################################################
  !#######################    density    ##########################################
  !$OMP DO
  do k=1,global_parameters%N_zones
     call compute_tridiag_system_density(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER 

  !$OMP DO
  do npsi=1,global_parameters%N_psi_surfaces
     do n_ion = 1,global_parameters%n_ions
        call concat_tridiag_system(flux_surfaces(npsi),n_ion)
        call apply_density_boundary_conditions(flux_surfaces(npsi),n_ion)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_density_BC(flux_surfaces(npsi),n_ion)
        call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,DENSITY_FIELD)
     end do
  end do
  !$OMP END DO
  !$OMP BARRIER
  !################################################################################
  !################################################################################

  !################################################################################
  !#######################    velocity    #########################################
  !$OMP DO
  do k=1,global_parameters%N_zones
     call compute_electric_field(zones(k))
     call compute_velocity_source_implicit_correction(zones(k))
     call compute_tridiag_system_velocity(zones(k))
  end do
  !$OMP END DO
  !$OMP BARRIER 
  !$OMP DO
  do npsi=1,global_parameters%N_psi_surfaces
     do n_ion = 1,global_parameters%n_ions
        call concat_tridiag_system(flux_surfaces(npsi),n_ion)
        if(drift_flags%solve_drift) then				!### Leybros modif ###
           call apply_velocity_drift_boundary_conditions(flux_surfaces(npsi),n_ion)
	else
           call apply_velocity_boundary_conditions(flux_surfaces(npsi),n_ion)
	end if
        call invert_psi_tridiag(flux_surfaces(npsi))
        if(drift_flags%solve_drift) then				!### Leybros modif ###
           call set_velocity_drift_BC(flux_surfaces(npsi),n_ion)
	else
           call set_velocity_BC(flux_surfaces(npsi),n_ion)
	end if
        call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,VELOCITY_FIELD)
     end do
  end do
  !$OMP END DO
  !$OMP BARRIER
  !################################################################################
  !################################################################################

  !$OMP DO
  do k=1,global_parameters%N_zones
     call clean_solution1(zones(k))
     call compute_velocity_and_Mach(zones(k),STEP_NEW)
     call compute_electron_density_and_velocity(zones(k),STEP_NEW)
  end do
  !$OMP END DO
  !$OMP BARRIER 

  if(flags%solve_temperature) then
     !################################################################################
     !#######################    temperature   #######################################
     !$OMP DO
     do k=1,global_parameters%N_zones
        call compute_temperature_source_implicit_correction(zones(k)) 
        call add_temperature_parallel_viscosity_implicit_term(zones(k))
        call compute_tridiag_system_temperature(zones(k))
     end do
     !$OMP END DO
     !$OMP BARRIER 
     !$OMP DO
     do npsi=1,global_parameters%N_psi_surfaces
        do n_ion = 0,global_parameters%n_ions
           call concat_tridiag_system(flux_surfaces(npsi),n_ion)
	   if(drift_flags%solve_drift) then				!### Leybros modif ###
              call apply_temperature_drift_boundary_conditions(flux_surfaces(npsi),n_ion)
	   else
              call apply_temperature_boundary_conditions(flux_surfaces(npsi),n_ion)
	   end if
           call invert_psi_tridiag(flux_surfaces(npsi))
	   if(drift_flags%solve_drift) then				!### Leybros modif ###
              call set_temperature_drift_BC(flux_surfaces(npsi),n_ion)
	   else
              call set_temperature_BC(flux_surfaces(npsi),n_ion)
	   end if
           call unconcat_tridiag_system(flux_surfaces(npsi),n_ion,TEMPERATURE_FIELD)
        end do
     end do
     !$OMP END DO
     !$OMP BARRIER
     !################################################################################
     !################################################################################

  end if

  if(flags%turbulence_model.eq.1) then
     !################################################################################
     !#######################       k       ##########################################
     !$OMP DO
     do k=1,global_parameters%N_zones
        call compute_kepsilon_source_implicit_correction(zones(k))
        call compute_tridiag_system_k_all_in_one(zones(k))
     end do
     !$OMP END DO
     !$OMP BARRIER 

     !$OMP DO
     do npsi=1,global_parameters%N_psi_surfaces
        call concat_tridiag_system(flux_surfaces(npsi),0)
        call apply_k_boundary_conditions(flux_surfaces(npsi),0)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_k_BC(flux_surfaces(npsi))
        call unconcat_tridiag_system(flux_surfaces(npsi),0,K_FIELD)
     end do
     !$OMP END DO
     !$OMP BARRIER
     !################################################################################
     !################################################################################

     !################################################################################
     !#######################   epsilon     ##########################################
     !$OMP DO
     do k=1,global_parameters%N_zones
        call compute_tridiag_system_epsilon_all_in_one(zones(k))
     end do
     !$OMP END DO
     !$OMP BARRIER 

     !$OMP DO
     do npsi=1,global_parameters%N_psi_surfaces
        call concat_tridiag_system(flux_surfaces(npsi),0)
        call apply_epsilon_boundary_conditions(flux_surfaces(npsi),0)
        call invert_psi_tridiag(flux_surfaces(npsi))
        call set_epsilon_BC(flux_surfaces(npsi))
        call unconcat_tridiag_system(flux_surfaces(npsi),0,EPSILON_FIELD)
     end do
     !$OMP END DO
     !$OMP BARRIER
     !################################################################################
     !################################################################################
  end if

  !$OMP DO
  do k=1,global_parameters%N_zones
     call compute_implicit_fluxes(zones(k))
     call clean_solution2(zones(k))
     if(flags%turbulence_model.eq.1) then
        call clean_solution_ke(zones(k))
     end if
     call check_fields(zones(k))
     call compute_local_residuals(zones(k))
  end do
  !$OMP END DO

  !$OMP BARRIER
  !$OMP MASTER
  call compute_residuals()

  global_variables%dt_vort=global_variables%dt_vort&
       +global_variables%dt

  global_variables%tempus=global_variables%tempus&
       +global_variables%dt*reference_parameters%fields%tau0
  !$OMP END MASTER 

  !$OMP END PARALLEL


end subroutine step_in_soledge
