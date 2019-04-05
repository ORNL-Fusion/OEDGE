module allocate_storage_div


contains


  subroutine allocate_dynamic_storage

    ! comsrc
       ! use mod_adas_data_spec
       use mod_cadas
       use mod_cadas2
       use mod_cedge2d
       use mod_cgeom
       use mod_cioniz
       use mod_clocal
       use mod_cneut2
       use mod_cnoco
       use mod_colours
       ! use mod_comgra
       use mod_commv
       ! use mod_comprt
       use mod_comsol
       use mod_comtor
       use mod_div1
       use mod_div2
       use mod_divxy
       ! use mod_comxs
       use mod_cyield
       ! use mod_dperpz
       use mod_driftvel
       !use mod_grmdata
       ! use mod_grminfo
       use mod_line_profile
       ! use mod_hc_global_opts
       ! use mod_outbuffer
       use mod_outcom
       ! use mod_out_unstruc
       use mod_parmmod
       ! use mod_particle_specs
       use mod_pindata
       ! use mod_reader
       use mod_reiser_com
       ! use mod_signal_com
       use mod_slcom
       use mod_slout
       ! use mod_sol22pcx
       ! use mod_sol22pei
       ! use mod_sol22phelpi
       ! use mod_sol22pmom
       ! use mod_sol23_input
       use mod_solcommon
       ! use mod_solparams
       use mod_solswitch
       use hc_storage_setup
       
    ! div6

       use mod_adpak_com
       use mod_baffles
       use mod_cfd_osm_com
       use mod_pin_cfd
       use mod_sol23_com
       use mod_cneut
       use mod_cneut2
       use mod_comhr
       use mod_crand
       use mod_diagvel
       use mod_div3
       ! use mod_div4
       use mod_div5
       use mod_div6
       use mod_div7
       use mod_divbra
       use mod_dynam1
       use mod_dynam2
       use mod_dynam3
       use mod_dynam4
       use mod_dynam5
       use mod_fperiph_com
       !use mod_hc_com
       use mod_grbound
       use mod_local_baffles
       use mod_inel
       use mod_promptdep
       use mod_slcom_sol28
       use mod_solrk
       use mod_walls_com
       use mod_temp
       use mod_transcoef


       
    implicit none

    !
    ! Replacement for DIVIMP common blocks by dynamic allocation of arrays
    ! - call to this routine may be moved if/when the array size parameters
    !   are moved to inputs
    !

    ! comsrc
       ! call allocate_mod_adas_data_spec
       call allocate_mod_cadas
       call allocate_mod_cadas2
       call allocate_mod_cedge2d
       call allocate_mod_cgeom
       call allocate_mod_cioniz
       call allocate_mod_clocal
       call allocate_mod_cneut2
       call allocate_mod_cnoco
       call allocate_mod_colours
       ! call allocate_mod_comgra
       call allocate_mod_commv
       ! call allocate_mod_comprt
       call allocate_mod_comsol
       call allocate_mod_comtor
       call allocate_mod_div1
       call allocate_mod_div2
       call allocate_mod_divxy
       ! call allocate_mod_comxs
       call allocate_mod_cyield
       ! call allocate_mod_dperpz
       call allocate_mod_driftvel
       !call allocate_mod_grmdata
       ! call allocate_mod_grminfo
       call allocate_mod_line_profile
       ! call allocate_mod_hc_global_opts
       ! call allocate_mod_outbuffer
       call allocate_mod_outcom
       ! call allocate_mod_out_unstruc
       call allocate_mod_parmmod
       ! call allocate_mod_particle_specs
       call allocate_mod_pindata
       ! call allocate_mod_reader
       call allocate_mod_reiser_com
       ! call allocate_mod_signal_com
       call allocate_mod_slcom
       call allocate_mod_slout
       ! call allocate_mod_sol22pcx
       ! call allocate_mod_sol22pei
       ! call allocate_mod_sol22phelpi
       ! call allocate_mod_sol22pmom
       ! call allocate_mod_sol23_input
       call allocate_mod_solcommon
       ! call allocate_mod_solparams
       call allocate_mod_solswitch
       call allocate_hc_storage

    ! div6

       call allocate_mod_adpak_com
       call allocate_mod_baffles
       call allocate_mod_cfd_osm_com
       call allocate_mod_pin_cfd
       call allocate_mod_sol23_com
       call allocate_mod_cneut
       call allocate_mod_cneut2
       call allocate_mod_comhr
       call allocate_mod_crand
       call allocate_mod_diagvel
       call allocate_mod_div3
       ! call allocate_mod_div4
       call allocate_mod_div5
       call allocate_mod_div6
       call allocate_mod_div7
       call allocate_mod_divbra
       call allocate_mod_dynam1
       call allocate_mod_dynam2
       call allocate_mod_dynam3
       call allocate_mod_dynam4
       call allocate_mod_dynam5
       call allocate_mod_fperiph_com
       !call allocate_mod_hc_com
       call allocate_mod_grbound
       call allocate_mod_local_baffles
       call allocate_mod_inel
       call allocate_mod_promptdep
       call allocate_mod_slcom_sol28
       call allocate_mod_solrk
       call allocate_mod_walls_com
       call allocate_mod_temp
       call allocate_mod_transcoef




  end subroutine allocate_dynamic_storage



  subroutine deallocate_dynamic_storage
    
    ! comsrc
       ! use mod_adas_data_spec
       use mod_cadas
       use mod_cadas2
       use mod_cedge2d
       use mod_cgeom
       use mod_cioniz
       use mod_clocal
       use mod_cneut2
       use mod_cnoco
       use mod_colours
       ! use mod_comgra
       use mod_commv
       ! use mod_comprt
       use mod_comsol
       use mod_comtor
       use mod_div1
       use mod_div2
       use mod_divxy
       ! use mod_comxs
       use mod_cyield
       ! use mod_dperpz
       use mod_driftvel
       !use mod_grmdata
       ! use mod_grminfo
       use mod_line_profile
       ! use mod_hc_global_opts
       ! use mod_outbuffer
       use mod_outcom
       ! use mod_out_unstruc
       use mod_parmmod
       ! use mod_particle_specs
       use mod_pindata
       ! use mod_reader
       use mod_reiser_com
       ! use mod_signal_com
       use mod_slcom
       use mod_slout
       ! use mod_sol22pcx
       ! use mod_sol22pei
       ! use mod_sol22phelpi
       ! use mod_sol22pmom
       ! use mod_sol23_input
       use mod_solcommon
       ! use mod_solparams
       use mod_solswitch
       use hc_storage_setup
       
    ! div6

       use mod_adpak_com
       use mod_baffles
       use mod_cfd_osm_com
       use mod_pin_cfd
       use mod_sol23_com
       use mod_cneut
       use mod_cneut2
       use mod_comhr
       use mod_crand
       use mod_diagvel
       use mod_div3
       ! use mod_div4
       use mod_div5
       use mod_div6
       use mod_div7
       use mod_divbra
       use mod_dynam1
       use mod_dynam2
       use mod_dynam3
       use mod_dynam4
       use mod_dynam5
       use mod_fperiph_com
       !use mod_hc_com
       use mod_grbound
       use mod_local_baffles
       use mod_inel
       use mod_promptdep
       use mod_slcom_sol28
       use mod_solrk
       use mod_walls_com
       use mod_temp
       use mod_transcoef




    implicit none

    ! comsrc
       ! call deallocate_mod_adas_data_spec
       call deallocate_mod_cadas
       call deallocate_mod_cadas2
       call deallocate_mod_cedge2d
       call deallocate_mod_cgeom
       call deallocate_mod_cioniz
       call deallocate_mod_clocal
       call deallocate_mod_cneut2
       call deallocate_mod_cnoco
       call deallocate_mod_colours
       ! call deallocate_mod_comgra
       call deallocate_mod_commv
       ! call deallocate_mod_comprt
       call deallocate_mod_comsol
       call deallocate_mod_comtor
       call deallocate_mod_div1
       call deallocate_mod_div2
       call deallocate_mod_divxy
       ! call deallocate_mod_comxs
       call deallocate_mod_cyield
       ! call deallocate_mod_dperpz
       call deallocate_mod_driftvel
       !call deallocate_mod_grmdata
       ! call deallocate_mod_grminfo
       call deallocate_mod_line_profile
       ! call deallocate_mod_hc_global_opts
       ! call deallocate_mod_outbuffer
       call deallocate_mod_outcom
       ! call deallocate_mod_out_unstruc
       call deallocate_mod_parmmod
       ! call deallocate_mod_particle_specs
       call deallocate_mod_pindata
       ! call deallocate_mod_reader
       call deallocate_mod_reiser_com
       ! call deallocate_mod_signal_com
       call deallocate_mod_slcom
       call deallocate_mod_slout
       ! call deallocate_mod_sol22pcx
       ! call deallocate_mod_sol22pei
       ! call deallocate_mod_sol22phelpi
       ! call deallocate_mod_sol22pmom
       ! call deallocate_mod_sol23_input
       call deallocate_mod_solcommon
       ! call deallocate_mod_solparams
       call deallocate_mod_solswitch
       call deallocate_hc_storage

    ! div6

       call deallocate_mod_adpak_com
       call deallocate_mod_baffles
       call deallocate_mod_cfd_osm_com
       call deallocate_mod_pin_cfd
       call deallocate_mod_sol23_com
       call deallocate_mod_cneut
       call deallocate_mod_cneut2
       call deallocate_mod_comhr
       call deallocate_mod_crand
       call deallocate_mod_diagvel
       call deallocate_mod_div3
       ! call deallocate_mod_div4
       call deallocate_mod_div5
       call deallocate_mod_div6
       call deallocate_mod_div7
       call deallocate_mod_divbra
       call deallocate_mod_dynam1
       call deallocate_mod_dynam2
       call deallocate_mod_dynam3
       call deallocate_mod_dynam4
       call deallocate_mod_dynam5
       call deallocate_mod_fperiph_com
       !call deallocate_mod_hc_com
       call deallocate_mod_grbound
       call deallocate_mod_local_baffles
       call deallocate_mod_inel
       call deallocate_mod_promptdep
       call deallocate_mod_slcom_sol28
       call deallocate_mod_solrk
       call deallocate_mod_walls_com
       call deallocate_mod_temp
       call deallocate_mod_transcoef



  end subroutine deallocate_dynamic_storage


end module allocate_storage_div
