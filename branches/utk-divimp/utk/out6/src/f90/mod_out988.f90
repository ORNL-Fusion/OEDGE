!     -*-Mode:f90-*-
MODULE MOD_OUT988

  IMPLICIT none
  PRIVATE

! Module routine declarations:
  PUBLIC :: ALLOC_DUM, DEALLOC_DUM

! Derived types:
  TYPE, PUBLIC :: type_options988
     INTEGER :: idum1
  ENDTYPE type_options988


! Allocatables:
  TYPE(type_options988), PUBLIC, ALLOCATABLE, SAVE :: dum1(:)

!  Public variables:
!  INTEGER, PUBLIC, SAVE :: NOBJ, NPIXEL, NCHORD

  ! Parameters:
  INTEGER    P_DUM1
  PARAMETER (P_DUM1=1)


CONTAINS

  SUBROUTINE ALLOC_DUM(ndum1)
    INTEGER, INTENT(IN) :: ndum1
    IF (ALLOCATED(dum1)) THEN  
       WRITE(0,*) 'DEALLOCATING!'
       DEALLOCATE (dum1)
    ENDIF
    WRITE(0,*) 'ALLOCATING!',ndum1
    ALLOCATE (dum1(ndum1))
    RETURN
  END SUBROUTINE ALLOC_DUM

  SUBROUTINE DEALLOC_DUM
    IF (.NOT.ALLOCATED(dum1)) RETURN
    DEALLOCATE (dum1)
    RETURN
  END SUBROUTINE DEALLOC_DUM

END MODULE MOD_OUT988
