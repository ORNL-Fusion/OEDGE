module mod_sl_eircom

implicit none


! jdemod
! This module contains eirene support arrays that used to be in common blocks
! which now require dynamic allocation due to conversion of relevant parameters
! to variables

!      COMMON  /EIRWALCOM/ walln,wallr,wallz,wallw,ebgki
!      INTEGER walln,wallw(MAXPTS),ebgki
!      REAL    wallr(MAXPTS,2),wallz(MAXPTS,2)

      INTEGER :: walln,ebgki
      INTEGER,allocatable :: wallw(:)
      REAL,allocatable  :: wallr(:,:),wallz(:,:)

      
      ! NOTE: This grid common block is not related to EIRENE but a module was needed to
      !       put it in that was least likely to conflict with others since mod_slcom
      !       did not work and it was undesirable to create another module.
      !       At some point in time it would be nice to clean up the code but it would be a
      !       significant effort for mostly aesthetic gain. 
  ! jdemod - moving more common blocks into mod_slcom
  ! COMMON /GRID/ iktop,irout,irin
  ! INTEGER       iktop(MAXNRS),irout(MAXNRS),irin(MAXNRS)
  INTEGER,public,allocatable ::  iktop(:),irout(:),irin(:)

      
contains


  subroutine allocate_mod_sl_eircom
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr
    call allocate_array(wallw,maxpts,'wallw',ierr)
    call allocate_array(wallr,maxpts,2,'wallr',ierr)
    call allocate_array(wallz,maxpts,2,'wallz',ierr)

    call allocate_array(iktop,maxnrs,'iktop',ierr)
    call allocate_array(irout,maxnrs,'irout',ierr)
    call allocate_array(irin,maxnrs,'irin',ierr)

  end subroutine allocate_mod_sl_eircom


  subroutine deallocate_mod_sl_eircom
    implicit none

    if (allocated(wallw)) deallocate(wallw)
    if (allocated(wallr)) deallocate(wallr)
    if (allocated(wallz)) deallocate(wallz)

    if (allocated(iktop)) deallocate(iktop)
    if (allocated(irout)) deallocate(irout)
    if (allocated(irin)) deallocate(irin)

  end subroutine deallocate_mod_sl_eircom


end module mod_sl_eircom
