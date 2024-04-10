!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
 
 
      FUNCTION EIRENE_FEHVDS1 (IREI,K)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: EIRENE_FEHVDS1, EHVDS, EIRENE_FTABEI1,
     .            EIRENE_RATE_COEFF, ERATE
      INTEGER :: II, KK
 
      EIRENE_FEHVDS1=0.D0
      KK=NREAHV(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            EIRENE_FEHVDS1=EHVDS1(IREI,1)
        CASE (-2)
            EIRENE_FEHVDS1=6.
        CASE (-3)
            EIRENE_FEHVDS1=10.0
        CASE (-4)
            EIRENE_FEHVDS1=8.6
        CASE (-5)
            EIRENE_FEHVDS1=0.5
        CASE (-6)
            EIRENE_FEHVDS1=0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        EHVDS = EIRENE_RATE_COEFF(KK,TEINL(K),0._DP,.FALSE.,0,ERATE)
        EHVDS=EXP(MAX(-100._DP,EHVDS+FACREA(KK,2)))
        EIRENE_FEHVDS1=EHVDS*DEIN(K)/(EIRENE_FTABEI1(IREI,K)+EPS60)
      END IF
 
      RETURN
      END
