C
      SUBROUTINE PLTIN (RIB,FCN,ANF,END,INN,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CLMSUR
      USE CTRCEI

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: RIB, ANF, END
      REAL(DP) :: FCN
      INTEGER, INTENT(OUT) :: INN
      LOGICAL :: LBOX
      EXTERNAL FCN

      IF (RIB.EQ.1..OR.RIB.EQ.0.)
     .                CALL PLTKI (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.EQ.1.5) CALL PLTKA (FCN,ANF,END,INN,TRCPLT,LBOX)
      IF (RIB.LT.0.)  CALL PLTKU (FCN,ANF,END,INN,TRCPLT,LBOX)
      RETURN
      END