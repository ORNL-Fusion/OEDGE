

      FUNCTION FEHVDS1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: FEHVDS1, EHVDS, FTABEI1
      INTEGER :: II, KK

      FEHVDS1=0.D0
      KK=NREAHV(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEHVDS1=EHVDS1(IREI,1)
        CASE (-2)
            FEHVDS1=6.
        CASE (-3)
            FEHVDS1=10.0
        CASE (-4)
            FEHVDS1=8.6
        CASE (-5)
            FEHVDS1=0.5
        CASE (-6)
            FEHVDS1=0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        EHVDS = CREAC(9,1,KK)
        DO II=8,1,-1
          EHVDS = EHVDS*TEINL(K) + CREAC(II,1,KK)
        END DO
        EHVDS=EXP(MAX(-100._DP,EHVDS+FACREA(KK)))
        FEHVDS1=EHVDS*DEIN(K)/(FTABEI1(IREI,K)+EPS60)
      END IF

      RETURN
      END
