subroutine styx_init_eirene(ical)
  use all_variables, only : interp_data2, global_parameters, reference_parameters
  use Mphysics
  use styx2eirene
  use eirmod_precision ! was forgotten
  use eirmod_comprt, only : iunout
  use eirmod_parmmod, only : NSTS,NLIM,NSPZ
  ! sheath
  use eirmod_ctrig, only : PTRIX,PTRIY
  use eirmod_ccona, only : eps60

  use eirmod_cstep
  use eirmod_comusr
  use eirmod_cpes
  implicit none
  
  integer*4,intent(in) :: ical 
 
  include 'mpif.h'

  if (ical == 1) then

    call styx_setup_general_parameters()
    call styx_init_surface_data(1)
!#ifdef S2D
!    call map_triangles_to_soledge_mesh()
!#endif

    call styx_set_time_dependent_mode(timedep)
    call styx_set_eirene_magnetic_field()
    call styx_get_singly_charged_ions()
    call styx_initialize_pfc_models()
    call styx_read_cx_reactions()

    call set_eirene_input_file()

    call styx_complete_first_initialization_phase()

  elseif (ical == 2) then
  
    if (my_pe == 0) then 
      call styx_check_input_parameters()
      call styx_init_surface_data(2)
    endif

    call styx_init_sheath_model()
  
    if (my_pe == 0) then
      call eirene_alloc_cstep()
      call styx_get_amd_tweaks()
     ! get cell volumes to rescale styx cells which have been cut
      call styx_get_eirene_volumes
      call styx_initialize_short_cycling_procedure()
      call styx_establish_species_index_mappings()
    endif

    call styx_allocate_eirene_interface_structures()
  endif ! ical == 2


end subroutine styx_init_eirene
