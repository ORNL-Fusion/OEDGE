module mod_grminfo

  use mod_params



!c
!c     Common block information for GRM plots 
!c
!      common /grminfo/ pageplots,textsize,axistextsize
!c
!      integer pageplots,textsize,axistextsize
!c



  implicit none
  private

      integer,public:: pageplots,textsize,axistextsize

  
  public :: allocate_mod_grminfo, deallocate_mod_grminfo


contains

  subroutine allocate_mod_grminfo
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)


  end subroutine allocate_mod_grminfo


  subroutine deallocate_mod_grminfo
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_grminfo



end module mod_grminfo
