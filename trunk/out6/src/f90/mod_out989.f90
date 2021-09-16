!     -*-Mode:f90-*-
MODULE MOD_OUT989

  IMPLICIT none
  PRIVATE

! Module routine declarations:
  PUBLIC :: ALLOC_DUM, DEALLOC_DUM

!...
  INTEGER, PUBLIC, PARAMETER :: MAXNIMAGE = 10

! Derived types:
  TYPE, PUBLIC :: type_options989
     INTEGER   :: imagetype
     CHARACTER :: fimage*512
     INTEGER   :: nxbin,nybin
     INTEGER   :: nimg
     INTEGER   :: img_type (MAXNIMAGE)
     REAL      :: img_scale(MAXNIMAGE)
     CHARACTER :: img_fname(MAXNIMAGE)*512
     INTEGER   :: img_nxbin(MAXNIMAGE)
     INTEGER   :: img_nybin(MAXNIMAGE)
     INTEGER   :: maptype
     CHARACTER :: fmap*512
     INTEGER   :: out_opt
     CHARACTER :: out_suffix*512
  ENDTYPE type_options989

  TYPE, PUBLIC :: type_header
!... Standard:
     INTEGER   :: shot
     INTEGER   :: frame
     REAL      :: time     
     INTEGER   :: channel
!... Old:
     CHARACTER :: id*8             
     INTEGER   :: size
     CHARACTER :: codec*8 
     CHARACTER :: date_time*20
     REAL      :: trigger
     CHARACTER :: lens*24
     CHARACTER :: filter*24
     CHARACTER :: view*64
     INTEGER   :: numFrames
     CHARACTER :: camera*64
     INTEGER*2 :: width      
     INTEGER*2 :: height
     INTEGER*2 :: depth     
     INTEGER   :: orient
     INTEGER*2 :: taps 
     INTEGER*2 :: color
     INTEGER*2 :: hBin 
     INTEGER*2 :: left
     INTEGER*2 :: right
     INTEGER*2 :: vBin
     INTEGER*2 :: top
     INTEGER*2 :: bottom
     INTEGER*2 :: offset(2)
     REAL      :: gain(2)
     INTEGER   :: preExp
     INTEGER   :: shutter
     INTEGER   :: strobe
     REAL      :: temperature
  ENDTYPE type_header

!  Allocatables:
  TYPE(type_header), PUBLIC, ALLOCATABLE, SAVE :: dum1(:)

! Public variables:
!  INTEGER, PUBLIC, SAVE :: m, n
!  REAL*8, PUBLIC, ALLOCATABLE, SAVE :: A(:,:), x(:), b(:)

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

END MODULE MOD_OUT989
