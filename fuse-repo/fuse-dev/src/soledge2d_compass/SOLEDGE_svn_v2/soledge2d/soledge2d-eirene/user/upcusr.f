C
C
      SUBROUTINE EIRENE_UPCUSR(WS,IND)
C
C  USER SUPPLIED COLLISION ESTIMATOR, VOLUME AVERAGED
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_COMXS
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: WS
      INTEGER, INTENT(IN) :: IND
C
C     WS=WEIGHT/SIGTOT=WEIGHT/(VEL*ZMFPI)=WEIGHT/(VEL*SIGMA,MACR.)
C
C  FOR PARTICLE DENSITY IN CELL NO. NCELL
C     COLV(1,NCELL)=COLV(1,NCELL)+WS
C
      RETURN
      END
