module mod_save

  use mod_params


!C                                                                               
!      COMMON /SAVE/ R,I,L                                                       
!      REAL          R(9,MAXPUT)                                                 
!      INTEGER       I(9,MAXPUT)                                                 
!      LOGICAL       L(3,MAXPUT)                                                 

  implicit none
  private

      REAL,public::          R(9,MAXPUT)                                                 
      INTEGER,public::       I(9,MAXPUT)                                                 
      LOGICAL,public::       L(3,MAXPUT)                                                 
 
  public :: allocate_mod_save, deallocate_mod_save


contains

  subroutine allocate_mod_save
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_save


  subroutine deallocate_mod_save
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_save



end module mod_save
