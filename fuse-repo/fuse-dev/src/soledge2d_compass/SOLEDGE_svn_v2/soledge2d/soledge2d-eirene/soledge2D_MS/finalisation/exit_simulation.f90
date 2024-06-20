subroutine exit_simulation(my_pe)
#include "compile_opt.inc"
  use Mlog_message
  use all_variables
  use Mdefinitions
  use Mfeedback_control
  use Mneutral_vars
  implicit none
  integer,intent(in) :: my_pe
  integer*4 :: k,i,j,n

  if (my_pe == 0) then
     do k=1,global_parameters%N_zones
        call MD_broadcast_plasma(zones(k),STEP_OLD)
     end do
     if(flags%use_triangles) then
        call interpolate_plasma()
        call interpolate_wall_fluxes()
        call interpolate_extra_fields()
        call compute_Zeff()
     end if
     call save_plasma()
     call save_fluxes()
     call save_simulation_global_parameters()
     call write_log_message(msg_soledge_ended)
#if VORTICITY_PASTIX == 1
     call save_vorticity_fields()
     call save_currents()
!!$     do k=1,global_parameters%N_zones
!!$        call current_balance(zones(k))
!!$     end do
#endif
     call save_magnetics()
     call save_control(Control_data,error_data)

     if(drift_flags%solve_phi) then
	call save_drift()
     endif

     if(flags%turbulence_model.eq.1) then
        call save_kepsilon()
     end if

     if(flags%neutral_model.eq.2) then
        call clear_pastix(CSC_neutral)
     end if

     call save_transport_coefficients()

  endif

#if USE_EIRENE == 1
  call exit_eirene_coupling()
#endif

end subroutine exit_simulation
