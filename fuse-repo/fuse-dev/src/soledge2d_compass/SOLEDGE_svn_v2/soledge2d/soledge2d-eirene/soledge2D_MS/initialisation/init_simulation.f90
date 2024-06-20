subroutine init_simulation(my_pe)

#include "compile_opt.inc"

  use all_variables
  use Mlog_message
  use MradialFeedback
  use Mdefinitions
#if USE_EIRENE == 1
  use Meirene_vars
  use Mfeedback_control
#endif
  implicit none
  integer,intent(in) :: my_pe
  integer*4 :: ier
  integer*4 :: k

#if USE_EIRENE == 1
  include 'mpif.h'
#endif

  if (my_pe == 0) then
     Interp_Data2%n_tor=1
     Interp_Data2%ang_max=360.d0
     call write_log_message(msg_soledge_started)
     call read_input_file()
     call read_drift_flags()
     if(flags%is_SLAB) then
        call read_slab_file()
     end if
     call read_mesh()

     call init_super_tridiag()

     allocate(global_variables%dt_table(1:global_parameters%N_zones))
     call allocate_arrays_memory()

     call prepare_penalisation_masks()

     call compute_reference_geometry()
     if(flags%is_slab) then
        call compute_metric_coefficients_slab()
     else
        call compute_jacobian()
        call compute_metric_coefficients()
     end if
     call compute_sign_metric()
     call make_Xpoint_cells_orthogonal()

     call compute_reference_parameters()  
     call compute_transport_coefficients()

     call init_plasma()

     if(flags%radialFeedback) then
        call read_feedback_data()
        if(flags%restart) then
           call load_transport_coefficients() 
        else
           radialFeedbackData%D=radialFeedbackData%input_D
           radialFeedbackData%chie=radialFeedbackData%input_Chi
           radialFeedbackData%chii=radialFeedbackData%input_Chi
           call init_feedback_plasma()
        end if
        do k=1,global_parameters%N_zones
           call update_radial_diffusivity_feedback(zones(k))
        end do
     end if

     call compute_implicit_tridiag_coefficients()

     call load_amdata()
     call allocate_am_vars()

     call allocate_integrals()

     inquire(File='triangles.h5',exist=flags%use_triangles)
     if(flags%use_triangles) then
        call init_triangles()
        call interpolate_magnetic()
     end if

     if(.not.flags%restart) then
        global_variables%tempus=0.d0
     else
        call load_simulation_global_parameters()
     end if

     allocate(global_variables%residuals(0:global_parameters%N_ions))

     call load_penalisation_parameters()
     global_variables%Teps=1.d0/reference_parameters%fields%T0eV

#if VORTICITY_PASTIX == 1
     call init_vorticity()
     if(drift_flags%solve_phi) then
        call init_drifts()
     end if
#endif

     if(flags%turbulence_model.eq.1) then
        call init_kepsilon()
     end if
     if(flags%neutral_model.eq.2) then
        global_variables%dt = 1.d10
        ! initialisation of dt (the value should not be used / just here to avoid division by zero)
        call init_fluid_neutrals()
     end if

     call init_custom_plots()

     call write_neutral_model()

     call write_log_message(msg_loop_started)

  endif

#if USE_EIRENE == 1
  call mpi_bcast(global_parameters%n_iterations,1,mpi_integer,0,mpi_comm_world,ier)
  call mpi_bcast(flags%neutral_model,1,mpi_logical,0,mpi_comm_world,ier)
  if(flags%neutral_model.eq.1) then
     call init_eirene_coupling()
     if (my_pe == 0) then
        if(eirene_vars%feedback.ne.0) then
           call init_control(Control_data,error_data)
        end if
     endif
     call styx_compute_pump_surface()
  end if
#endif

end subroutine init_simulation
