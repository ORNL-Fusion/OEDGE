program soledge2d

#include "compile_opt.inc"

  use all_variables
  use ifport
  use Mdefinitions

#if USE_EIRENE == 1
  use styx2eirene, only : Sn_intg,dt_eirene
  use eirmod_cpes
#endif

  implicit none
  integer*4 :: n_ite,nspe,k
  logical :: is_eirene_iteration,shall_we_run_eirene
  double precision ::  omp_get_wtime
  double precision :: t_bg_ite,t_end_ite,t_ite

!!$#if USE_PASTIX == 1 && USE_EIRENE == 0
!!$  integer :: my_pe, nprs
!!$#endif
#if USE_EIRENE == 1 
  include 'mpif.h'
  integer*4 :: ier
  call mpi_init(IER)
  call mpi_comm_size (mpi_comm_world,nprs,ier)
  call mpi_comm_rank (mpi_comm_world,my_pe,ier)
#endif
#if USE_EIRENE == 0 
  integer :: my_pe
  my_pe=0 ! default value if no mpi
#endif

  call init_simulation(my_pe)
  if (my_pe==0) then
     t_bg_ite=0.d0
     t_end_ite=0.d0
     t_ite=0.d0
     open(unit=6, carriagecontrol='fortran')
  endif

  !Beginning of the temporal loop
  do n_ite = 1, global_parameters%n_iterations
     if(my_pe.eq.0) then
        t_bg_ite = omp_get_wtime()
        call step_in_soledge(n_ite)
        do nspe=2,global_parameters%N_ions
           do k=1,global_parameters%N_zones
              zones(k)%species(nspe)%var(1)%temperature=zones(k)%species(1)%var(1)%temperature
              zones(k)%species(nspe)%var(2)%temperature=zones(k)%species(1)%var(2)%temperature
           end do
        end do
        if(flags%neutral_model.eq.2) then
           call solve_neutrals()
        end if
        if(flags%radialFeedback) then
           call radial_profiles_diffusivity_feedback2()
        end if
     end if
     ! the following call to eirene is executed by every process
#if USE_EIRENE == 1
     is_eirene_iteration=shall_we_run_eirene()
     if(is_eirene_iteration .or. (my_pe > 0)) then
        call update_neutrals_with_eirene(n_ite)
     end if
#endif
#if VORTICITY_PASTIX == 1
     if(my_pe.eq.0) then
        if(drift_flags%solve_phi) then
           if(modulo(dble(n_ite),dble(drift_flags%N_solve_phi)).eq.0) then
	      global_variables%dt_vort_old=global_variables%dt_vort !save old time step
              call solve_vorticity()
              global_variables%dt_vort=0.d0 !reset
           end if
        end if
     end if
#endif
     if(my_pe.eq.0) then
        call monitor_fluxes_and_sources(n_ite)
        call custom_plot_save(n_ite)
        call update_all_fields()
        t_end_ite = omp_get_wtime()
        t_ite=t_end_ite-t_bg_ite
        call write_progress_bar(n_ite,40,t_ite)
        call save_progress(n_ite)
        call save_residuals(n_ite)
     end if
  end do
  !end of the temporal loop

  call exit_simulation(my_pe)

end program soledge2d
