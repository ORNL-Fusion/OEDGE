module mod_comgra

  use mod_params

!c
!c     Common block for graphics in OUT 
!c
!      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,NPLOTS,ISPOT
!      REAL            CXMIN,CXMAX,CYMIN,CYMAX
!      INTEGER         IPLOTS,NPLOTS,ISPOT
!c

  implicit none
  private

      REAL,public::            CXMIN,CXMAX,CYMIN,CYMAX
      INTEGER,public::         IPLOTS,NPLOTS,ISPOT
  
  public :: allocate_mod_comgra, deallocate_mod_comgra


contains

  subroutine allocate_mod_comgra
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

  end subroutine allocate_mod_comgra


  subroutine deallocate_mod_comgra
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_comgra



end module mod_comgra
