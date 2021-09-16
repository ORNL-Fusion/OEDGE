module mod_adas_data_spec

  use mod_params


!c        
!c     This common block exists to allow the load_divdata_array 
!c     rotuine to be reused without requiring the ADAS data specification
!c     to be reread from the input file. 
!c
!c     CADAS_SWITCH is initialized to zero which deactivates this option 
!c
!      common /adas_spec/ cisele,ciselr,ciselx,ciseld,cadasyr,
!     >                   cadas_switch,cadasid,cadasex
!      integer cisele,ciselr,ciselx,ciseld,cadasyr,cadas_switch
!      character cadasid*80,cadasex*3

  implicit none
  private


      integer,public::  cisele,ciselr,ciselx,ciseld,cadasyr,cadas_switch
      character,public::  cadasid*80,cadasex*3


  
  public :: allocate_mod_adas_data_spec, deallocate_mod_adas_data_spec


contains

  subroutine allocate_mod_adas_data_spec
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

  end subroutine allocate_mod_adas_data_spec


  subroutine deallocate_mod_adas_data_spec
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_adas_data_spec



end module mod_adas_data_spec
