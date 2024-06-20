subroutine styx_short_cycle()
  use all_variables, only : interp_data2, global_parameters, reference_parameters
  use eirmod_precision
  use eirmod_parmmod
  use eirmod_ccona
  use styx2eirene
  use eirmod_comusr
  use eirmod_comprt, only : iunout
  use eirmod_ccona
  implicit none

  integer*4 :: itri,itor,icell
  
  if (sc_level >= 3) then
    ! recalculate rate coefficients, energy losses
    call styx_feed_plasma_to_eirene
    call styx_feed_derived_data_to_eirene
  endif

  if (sc_level == 4) then
    call styx_evolve_neutral_fields_during_short_cycling
  endif

  if (sc_level >= 2) then
    ! rescale neutral fields and update derived data
    call styx_rescale_neutral_fields
    call styx_recalculate_neutrals_velocity_fields
  endif
      
  if (sc_level >= 2) then       
    call styx_preav_sources
  endif

!!!!!!!!!!!!!!!!!!   units at this point   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! atom density in m-3
! Sn in m-3.s-1
! Sm in (m.s-1.m-3).s-1
! SE in Watt.m-3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do itor=1,Ntor_cells
    do itri=1,Ntri_styx
      icell = itri + (itor-1)*(Ntri_styx+1)  
      Interp_data2%tri_Sn(itri,:,itor) = Sn_tot(icell,:)
      Interp_data2%tri_SG(itri,1:,itor) = Sm_tot(icell,:)
      Interp_data2%tri_SE(itri,:,itor) = SE_tot(icell,:)
    enddo
  enddo
end subroutine styx_short_cycle
