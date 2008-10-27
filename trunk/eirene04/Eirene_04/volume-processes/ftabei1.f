

      FUNCTION FTABEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), FTABEI1, DEIMIN, DSUB, PLS, TBEI
      INTEGER :: J, I, II, KK

      TBEI=0.D0
      KK = NREAEI(IREI)
      IF (JEREAEI(IREI) == 1) THEN
        IF (MOD(IFTFLG(KK,2),100) == 10) THEN
          TBEI=CREAC(1,1,KK)
        ELSE
          TBEI = CREAC(9,1,KK)
          DO II=8,1,-1
            TBEI = TBEI*TEINL(K) + CREAC(II,1,KK)
          END DO
          TBEI = TBEI + FACREA(KK)
          TBEI=EXP(MAX(-100._DP,TBEI))
        ENDIF
        IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)
      ELSE
        TBEIC(1:JEREAEI(IREI)) = CREAC(9,1:JEREAEI(IREI),KK)
        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
        PLS=MAX(DEIMIN,DEINL(K))-DSUB
        DO J=1,JEREAEI(IREI)
          DO II=8,1,-1
            TBEIC(J)=TBEIC(J)*TEINL(K)+CREAC(II,J,KK)
          END DO
        END DO
        TBEI = TBEIC(9)
        DO I=8,1,-1
          TBEI=TBEI*PLS+TBEIC(I)
        END DO
        TBEI=TBEI+FACREA(KK)
        IF (IFTFLG(KK,2) < 100) TBEI=TBEI+DEINL(K)
        TBEI=MAX(-100._DP,TBEI)
        TBEI=EXP(TBEI)
      ENDIF

      FTABEI1 = TBEI

      RETURN
      END

