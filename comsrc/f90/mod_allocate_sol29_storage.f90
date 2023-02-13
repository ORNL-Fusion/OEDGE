module mod_allocate_sol29_storage
  implicit none
  
  public:: allocate_sol29_storage, deallocate_sol29_storage
  
contains

  ! Allocate storage wrapper.
  subroutine allocate_sol29_storage
    use mod_solcommon
    implicit none
    
    call pr_trace('ALLOCATE_SOL29_STORAGE','ALLOCATE')
    call allocate_mod_solcommon29
    
  end subroutine allocate_sol29_storage
  
  ! Deallocate storage wrapper.
  subroutine deallocate_sol29_storage
    use mod_solcommon
    implicit none
    
    call pr_trace('ALLOCATE_SOL29_STORAGE','DEALLOCATE')
    call deallocate_mod_solcommon29
  
  end subroutine deallocate_sol29_storage

end module mod_allocate_sol29_storage
