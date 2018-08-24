module mod_colours

  use mod_params

!c     -*-Fortran-*-
!c
!c     Colour definitions for use in OUT
!c
!c
!      integer maxncols,maxmarkers
!      parameter(maxncols=255,maxmarkers=8)
!      
!c
!      common /colours/ colour,ncols,defcol,icol,start_col
!      save /colours/
!c
!      integer colour(maxncols),ncols,defcol,icol,col,start_col
!c 
!c     Set up selected markers for plotting
!c
!      integer marker(maxmarkers) 
!c
!      data marker /232,250,224,225,227,248,228,229/
!c


  implicit none
  private


!     Colour definitions for use in OUT

      integer,public::  maxncols,maxmarkers
      parameter(maxncols=255,maxmarkers=8)
!
      integer,public:: colour(maxncols),ncols,defcol,icol,col,start_col
! 
!     Set up selected markers for plotting
!
      integer,public:: marker(maxmarkers) 
!
      data marker /232,250,224,225,227,248,228,229/
!
  
  public :: allocate_mod_colours, deallocate_mod_colours


contains

  subroutine allocate_mod_colours
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

  end subroutine allocate_mod_colours


  subroutine deallocate_mod_colours
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

  end subroutine deallocate_mod_colours



end module mod_colours
