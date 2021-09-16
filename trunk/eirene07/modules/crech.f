!  16.01.07  subroutine for reinitialisation of point list introduced


      MODULE CRECH

      USE PRECISION

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: PPOINT,
     .          FIRST_POINT, LAST_POINT,
     .          INSTOR, INUMS, INOSF,
     .          LZR, CRECH_REINIT 

      TYPE :: PPOINT
        REAL(DP) :: XPL2D, YPL2D
        INTEGER :: NPL2D, NUMSUR
        TYPE(PPOINT), POINTER :: NXTPNT
      END TYPE PPOINT

      TYPE(PPOINT), POINTER :: FIRST_POINT, LAST_POINT
      INTEGER, SAVE :: INSTOR=0, INUMS, INOSF
      LOGICAL :: LZR

      contains

C     The following subroutine is for reinitialization of EIRENE (DMH)

      SUBROUTINE CRECH_REINIT
      IMPLICIT NONE
      TYPE(PPOINT), POINTER :: CUR

      CUR => FIRST_POINT
      DO WHILE (ASSOCIATED(CUR))
         FIRST_POINT => CUR
         CUR => CUR%NXTPNT
         DEALLOCATE(FIRST_POINT)
      END DO

      NULLIFY(FIRST_POINT)
      NULLIFY(LAST_POINT)

      INSTOR=0
      return
      end subroutine CRECH_REINIT
      

      END MODULE CRECH
