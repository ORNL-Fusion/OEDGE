subroutine run_monte_carlo()
  use all_variables, only : interp_data2, global_parameters, reference_parameters
  use Meirene_vars
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_cpes
  use eirmod_cestim
  use eirmod_comusr
  use eirmod_cstep
  use eirmod_clgin
  use eirmod_comnnl
  use eirmod_cspei
  use eirmod_csdvi
  use eirmod_czt1
  use eirmod_comsou
  use eirmod_coutau
  use eirmod_cupd

  use eirmod_ccona
  use eirmod_comxs

  use styx2eirene
  implicit none


  call styx_reinit_eirene
 
  if (my_pe == 0) then
     call styx_feed_plasma_to_eirene()     
     call styx_feed_derived_data_to_eirene()
     call eirene_get_fluxes()
     call styx_rescale_puff_rate()   
     call styx_feed_wall_parameters_to_eirene()
     if (timedep) call styx_load_neutrals_from_the_past()
  endif

  ! broadcast new data to be used in monte carlo calculation
  ! test TK3X
 ! if (nprs > 1) then
     call styx_broadcast_eirene_variables
 ! endif

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!                MONTE CARLO calculation                   !!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  call eirene_mcarlo

  ! get the data out (per strata, from group leaders)

  if (nprs > 1) then
    call styx_mpi_reduce_eirene_variables
  endif
 
  

  if (my_pe == 0) then
    ! get data out (sum over strata)
    call eirene_if4cop()
    ! store census array
    call eirene_tmstep()
  endif


end subroutine run_monte_carlo
