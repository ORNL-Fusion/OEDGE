module allocate_storage_div


contains


  subroutine allocate_dynamic_input
    use mod_comtor
    use mod_dynam4
    use mod_cgeom
    use mod_driftvel
    use mod_walls_com
    use mod_slcom
    use mod_solswitch
    use mod_parmmod
    
    implicit none

    ! jdemod - some code is executed as part of the input file processing - mostly written by Steve - so ALL of the variables in these
    !          modules will be allocated at the input step for now since they do not depend on maximp or maxizs

    ! arrays that are specifically used in the input file
    call allocate_mod_comtor_input
    call allocate_mod_dynam4_input
    call allocate_mod_cgeom_input
    call allocate_mod_driftvel_input
    call allocate_mod_walls_com_input
    call allocate_mod_solswitch_input
    call allocate_mod_slcom_input


    ! However, Steve has a lot of initialization code that is executed in the READIN routine. Ideally these would
    ! have come after read in but he has a lot of code mixed in. Trying to untangle the parts he initializes and also determine
    ! which initializations are actually required before the input file is READ is too much work at the present time.
    !
    ! Instead - everything except the variables dependent on maxizs and maximp will get allocated before the input
    ! file is read. The rest get allocated after. Careful checking to make sure that this works will be required. If it is ok
    ! we can migrate more variables to dynamic allocation afterward. At the present time, only maximp and maxizs are being used
    ! for true dynamic allocation. 


    
  end subroutine allocate_dynamic_input



  subroutine allocate_dynamic_before_input
  
       use mod_cadas2
       use mod_cedge2d
       use mod_cgeom
       use mod_cneut2
       use mod_cnoco
       use mod_colours
       use mod_comsol
       use mod_comtor
       use mod_div1
       use mod_div4
       use mod_divxy
       use mod_cyield
       use mod_driftvel
       use mod_line_profile
       use mod_parmmod
       use mod_pindata
       use mod_slout

       !use mod_solcommon
       !use mod_solswitch
       !use mod_solrk

       use mod_allocate_sol22_storage

    ! div6

       use mod_adpak_com
       use mod_baffles
       use mod_cfd_osm_com
       use mod_pin_cfd
       use mod_sol23_com
       use mod_cneut2
       use mod_comhr
       use mod_crand
       use mod_dynam5
       use mod_fperiph_com
       use mod_grbound
       use mod_local_baffles
       use mod_promptdep
       use mod_slcom_sol28
       use mod_walls_com
       use mod_temp
       use mod_transcoef
       use mod_slcom


    implicit none

    ! jdemod - some code is executed as part of the input file processing - mostly written by Steve - so ALL of the variables in these
    !          modules will be allocated at the input step for now since they do not depend on maximp or maxizs

       call allocate_mod_cadas2(maxnks)
       call allocate_mod_cedge2d
       call allocate_mod_cgeom
       call allocate_mod_cneut2
       call allocate_mod_cnoco
       call allocate_mod_colours
       call allocate_mod_comsol
       call allocate_mod_comtor
       call allocate_mod_div1
       call allocate_mod_div4
       call allocate_mod_divxy
       call allocate_mod_cyield
       call allocate_mod_driftvel
       call allocate_mod_line_profile
       call allocate_mod_parmmod
       call allocate_mod_pindata
       call allocate_mod_slout

       !call allocate_mod_solcommon
       !call allocate_mod_solswitch
       !call allocate_mod_solrk

       call allocate_sol22_storage
       
       call allocate_mod_adpak_com
       call allocate_mod_baffles
       call allocate_mod_cfd_osm_com
       call allocate_mod_pin_cfd
       call allocate_mod_sol23_com
       call allocate_mod_comhr
       call allocate_mod_crand
       call allocate_mod_dynam5
       call allocate_mod_fperiph_com
       call allocate_mod_grbound
       call allocate_mod_local_baffles
       call allocate_mod_promptdep
       call allocate_mod_slcom_sol28
       call allocate_mod_walls_com
       call allocate_mod_temp
       call allocate_mod_transcoef

       call allocate_mod_slcom
    
  end subroutine allocate_dynamic_before_input



  subroutine allocate_dynamic_storage(nimps,nimps2,nizs)

    ! jdemod - all of the modules included in this allocation routine contain variables with dependencies on
    !          quantities that are specified in the input file. So far these include
    !          - maxizs - maximum ionization state of the impurity particles
    !          - maximp - maximum number of impurity particles to be followed
    !
    !          by allocating these variables after input the memory usage of OEDGE can be constrained for each specific
    !          case
    !
       use mod_cadas
       use mod_cioniz
       use mod_clocal
       use mod_commv
       use mod_div2
       use mod_reiser_com
       use hc_storage_setup
       use comhc
      
   ! div6

       use mod_cneut
       use mod_diagvel
       use mod_div3
       use mod_div5
       use mod_div6
       use mod_div7
       use mod_divbra
       use mod_dynam1
       use mod_dynam2
       use mod_dynam3
       use mod_dynam4
       use mod_inel
       use mod_comtor
       
       use mod_sl_oldplasma
       use mod_sl_input

       use mod_sl_eircom
       use mod_sl_easmesh

       use bfield
       
    implicit none
    integer :: nimps,nimps2,nizs

    ! set parameters based on values that have been read in
    if (nimps.gt.0) then
       ! set maximum impurity count to be one greater than input just because :)
       maximp = nimps+nimps2+1
    endif

    max_impurities = maximp

    ! jdemod
    ! Set maximum ionization state to maximum specified in the input file
    ! Maxizs needs to be set to the maximum of cion and nizs
    ! Some of the ionization rate code in iztau calculates for all the states of the impurity
    ! while other parts stop at nizs. As a result, maxizs needs to be set to the maximum of
    ! cion and nizs
    
    maxvizs  = maxizs
    maxrtnsd = maxizs+1
    
    !
    ! Replacement for DIVIMP common blocks by dynamic allocation of arrays
    ! - call to this routine may be moved if/when the array size parameters
    !   are moved to inputs
    !

       call allocate_mod_cadas(maxnks,maxizs)
       call allocate_mod_cioniz
       call allocate_mod_clocal
       call allocate_mod_commv
       call allocate_mod_div2
       call allocate_mod_reiser_com
       call allocate_hc_storage
       call allocate_comhc
       
    ! div6

       call allocate_mod_cneut
       call allocate_mod_diagvel
       call allocate_mod_div3
       call allocate_mod_div5
       call allocate_mod_div6
       call allocate_mod_div7
       call allocate_mod_divbra
       call allocate_mod_dynam1
       call allocate_mod_dynam2
       call allocate_mod_dynam3
       call allocate_mod_dynam4
       call allocate_mod_inel

       call allocate_mod_sl_oldplasma
       call allocate_mod_sl_input

       call allocate_bfield

       call allocate_mod_sl_eircom
       call allocate_mod_sl_easmesh
       
  end subroutine allocate_dynamic_storage

  subroutine deallocate_dynamic_storage

    ! jdemod - deallocate all of the allocated variables in the code
    !          whenever an array is added to a module it should be added to the appropriate allocation
    !          routine (either normal or input if it is needed in the input block) and also added to the deallocation routine
    !
       use mod_cadas
       use mod_cadas2
       use mod_cedge2d
       use mod_cgeom
       use mod_cioniz
       use mod_clocal
       use mod_cneut2
       use mod_cnoco
       use mod_colours
       use mod_commv
       use mod_comsol
       use mod_comtor
       use mod_div1
       use mod_div2
       use mod_divxy
       use mod_cyield
       use mod_driftvel
       use mod_line_profile
       use mod_parmmod
       use mod_pindata
       use mod_reiser_com
       use mod_slcom
       use mod_slout

       !use mod_solcommon
       !use mod_solswitch
       !use mod_solrk
       use mod_allocate_sol22_storage

       use hc_storage_setup
       use comhc
       
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
       use mod_grbound
       use mod_local_baffles
       use mod_inel
       use mod_promptdep
       use mod_slcom_sol28
       use mod_walls_com
       use mod_temp
       use mod_transcoef

       use mod_sl_oldplasma
       use mod_sl_input

       use bfield

       use mod_sl_eircom
       use mod_sl_easmesh
       
       implicit none

       call deallocate_mod_cadas
       call deallocate_mod_cadas2
       call deallocate_mod_cedge2d
       call deallocate_mod_cgeom
       call deallocate_mod_cioniz
       call deallocate_mod_clocal
       call deallocate_mod_cneut2
       call deallocate_mod_cnoco
       call deallocate_mod_colours
       call deallocate_mod_commv
       call deallocate_mod_comsol
       call deallocate_mod_comtor
       call deallocate_mod_div1
       call deallocate_mod_div2
       call deallocate_mod_divxy
       call deallocate_mod_cyield
       call deallocate_mod_driftvel
       call deallocate_mod_line_profile
       call deallocate_mod_parmmod
       call deallocate_mod_pindata
       call deallocate_mod_reiser_com
       call deallocate_mod_slcom
       call deallocate_mod_slout

       !call deallocate_mod_solcommon
       !call deallocate_mod_solswitch
       !call deallocate_mod_solrk
       call deallocate_sol22_storage

       call deallocate_hc_storage
       call deallocate_comhc
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
       call deallocate_mod_grbound
       call deallocate_mod_local_baffles
       call deallocate_mod_inel
       call deallocate_mod_promptdep
       call deallocate_mod_slcom_sol28
       call deallocate_mod_walls_com
       call deallocate_mod_temp
       call deallocate_mod_transcoef

       call deallocate_mod_sl_oldplasma
       call deallocate_mod_sl_input
       call deallocate_bfield
       call deallocate_mod_sl_eircom

       call deallocate_mod_sl_easmesh

     end subroutine deallocate_dynamic_storage


end module allocate_storage_div
