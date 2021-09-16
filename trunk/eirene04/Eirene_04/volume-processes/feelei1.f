

      FUNCTION FEELEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) ELEIC(9), FEELEI1, PLS, DEIMIN, EE, FTABEI1, DSUB, ELEI
      INTEGER :: J, I, KK, II

      FEELEI1=0.D0
      KK=NELREI(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEELEI1=-EIONHE
        CASE (-2)
            FEELEI1=EELDS1(IREI,1)
        CASE (-3)
            FEELEI1=-1.5*TEIN(K)
        CASE (-4)
            FEELEI1=-EIONH
        CASE (-5)
            FEELEI1=-10.5
        CASE (-6)
            FEELEI1=-25.0
        CASE (-7)
            FEELEI1=EELDS1(IREI,1)
        CASE (-8)
            FEELEI1=-10.5
        CASE (-9)
            FEELEI1=-15.5
        CASE (-10)
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            FEELEI1=-0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN
          ELEI = CREAC(9,1,KK)
          DO II=8,1,-1
            ELEI = ELEI*TEINL(K) + CREAC(II,1,KK)
          END DO
          ELEI=EXP(MAX(-100._DP,ELEI))
          FEELEI1=-ELEI*DEIN(K)*FACREA(KK)/(FTABEI1(IREI,K)+EPS60)
        ELSE
          ELEIC(1:JELREI(IREI)) = CREAC(9,1:JELREI(IREI),KK)
          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          DO J=1,JELREI(IREI)
            DO II=8,1,-1
              ELEIC(J)=ELEIC(J)*TEINL(K)+CREAC(II,J,KK)
            END DO
          END DO
          ELEI = ELEIC(9)
          DO I=8,1,-1
            ELEI=ELEI*PLS+ELEIC(I)
          END DO
          EE=MAX(-100._DP,ELEI+FACREA(KK)+DEINL(K))
          FEELEI1=-EXP(EE)/(FTABEI1(IREI,K)+EPS60)
        END IF
      END IF

      RETURN
      END
