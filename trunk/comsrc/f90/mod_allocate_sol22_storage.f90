module mod_allocate_sol22_storage

  implicit none

  public:: allocate_sol22_storage, deallocate_sol22_storage

  
contains
  
  subroutine allocate_sol22_storage
    use debug_options
    use mod_solswitch
  use mod_solcommon
  use mod_solrk
    implicit none
    
    call pr_trace('ALLOCATE_SOL22_STORAGE','ALLOCATE')
       call allocate_mod_solcommon
       call allocate_mod_solswitch
       call allocate_mod_solrk


  end subroutine allocate_sol22_storage



  subroutine deallocate_sol22_storage
    use debug_options
  use mod_solswitch
  use mod_solcommon
  use mod_solrk
    implicit none
    call pr_trace('ALLOCATE_SOL22_STORAGE','DEALLOCATE')
    
       call deallocate_mod_solcommon
       call deallocate_mod_solswitch
       call deallocate_mod_solrk

  end subroutine deallocate_sol22_storage

  !subroutine allocate_solswitch_dynamic_input
  !  use mod_solswitch
  !  implicit none
  !
  !  call allocate_mod_solswitch_input
  !
  !end subroutine allocate_solswitch_dynamic_input


end module mod_allocate_sol22_storage
