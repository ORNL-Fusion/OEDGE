module mod_gcom1

  use mod_params

!C
!C     THIS COMMON BLOCK SUPPLIES THE FEW VARIABLES REQUIRED FOR
!C     THE SIMPLE ROUTINES WRITTEN TO SIMULATE THE GHOST PACKAGE
!C     ALL THEY DO IS SAVE THE DATA TO A FILE. CURRENTLY UNIT 12
!C
!      COMMON /GCOM1/ GNOUT,GPICS,GLIM
!      INTEGER GNOUT,GPICS,GLIM


  implicit none
  private

      INTEGER,public:: GNOUT,GPICS,GLIM
  
  public :: allocate_mod_gcom1, deallocate_mod_gcom1


contains

  subroutine allocate_mod_gcom1
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_gcom1


  subroutine deallocate_mod_gcom1
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_gcom1



end module mod_gcom1
