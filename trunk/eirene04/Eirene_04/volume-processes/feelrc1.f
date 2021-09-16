

      FUNCTION FEELRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: ELRC1(9), PLS, DELE, EE, FEELRC1, FTABRC1, DSUB, 
     .            DEIMIN, ELRC
      INTEGER :: J, I, KK, II

      FEELRC1=0.D0
      KK=NELRRC(IRRC)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        CASE (-2)
            FEELRC1=EELRC1(IRRC,1)*FTABRC1(IRRC,K)
        CASE (-3)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        END SELECT
      ELSE IF (KK > 0) THEN
        IF (JELRRC(IRRC) == 1) THEN
          ELRC = CREAC(9,1,KK)
          DO II=8,1,-1
            ELRC = ELRC*TEINL(K) + CREAC(II,1,KK)
          END DO
          ELRC=ELRC+FACREA(KK)
          ELRC=EXP(MAX(-100._DP,ELRC))
          FEELRC1=-ELRC*DEIN(K)
        ELSE
          ELRC1(1:JELRRC(IRRC)) = CREAC(9,1:JELRRC(IRRC),KK)
          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          DO J=1,JELRRC(IRRC)
            DO II=8,1,-1
              ELRC1(J)=ELRC1(J)*TEINL(K)+CREAC(II,J,KK)
            END DO
          END DO
          ELRC = ELRC1(9)
          DO I=8,1,-1
            ELRC=ELRC*PLS+ELRC1(I)
          END DO
          EE=MAX(-100._DP,ELRC+DEINL(K)+FACREA(KK))
          FEELRC1=-EXP(EE)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          FEELRC1=FEELRC1+DELE*FTABRC1(IRRC,K)
        END IF
      END IF

      RETURN
      END
