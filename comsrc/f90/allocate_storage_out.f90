module allocate_storage_out


contains


  subroutine allocate_dynamic_storage
    !use mod_dynam1
    use mod_dynam2
    use mod_dynam3
    use mod_div6
    use mod_cneut
    use mod_comtor
    use mod_dynam2
    !use mod_clocal
    use mod_cioniz
    use hc_storage_setup

    implicit none

    !
    ! Replacement for DIVIMP common blocks by dynamic allocation of arrays
    ! - call to this routine may be moved if/when the array size parameters
    !   are moved to inputs
    !

    !call allocate_dynam1
    call allocate_dynam2
    call allocate_dynam3
    call allocate_div6
    call allocate_cneut
    call allocate_comtor
    !call allocate_clocal
    call allocate_cioniz
    call allocate_hc_storage


  end subroutine allocate_dynamic_storage



  subroutine deallocate_dynamic_storage
    !use mod_dynam1
    use mod_dynam2
    use mod_dynam3
    use mod_div6
    use mod_cneut
    use mod_comtor
    use mod_dynam2
    !use mod_clocal
    use mod_cioniz
    use hc_storage_setup

    implicit none

    !call deallocate_dynam1
    call deallocate_dynam2
    call deallocate_dynam3
    call deallocate_div6
    call deallocate_cneut
    call deallocate_comtor
    !call deallocate_clocal
    call deallocate_cioniz
    call deallocate_hc_storage

  end subroutine deallocate_dynamic_storage


end module allocate_storage_out
