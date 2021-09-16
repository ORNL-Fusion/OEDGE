

      FUNCTION FTABRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: TBRCC(9), PLS(1), COUN(0:9,1), CF(9,0:9)
      REAL(DP) :: DEIMIN, DSUB, FTABRC1, ZX, TBRC
      INTEGER :: II, KK

      TBRC=0.D0
      KK = NREARC(IRRC)

      IF (KK == 0) THEN
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)
      ELSEIF (JEREARC(IRRC) == 1) THEN
        IF (MOD(IFTFLG(KK,2),100) == 10) THEN
          TBRC=CREAC(1,1,kk)
        ELSE
          TBRC = CREAC(9,1,KK)
          DO II=8,1,-1
            TBRC = TBRC*TEINL(K) + CREAC(II,1,KK)
          END DO
          TBRC=TBRC+FACREA(KK)
          TBRC=EXP(MAX(-100._DP,TBRC))
        END IF
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)
      ELSE
        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
        PLS(1)=MAX(DEIMIN,DEINL(K))-DSUB
        CALL CDEFN(TEINL(K),PLS,KK,COUN,1,CF,.TRUE.,.FALSE.,
     ,             .TRUE.)
        TBRC=COUN(1,1)*FACREA(KK)
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)
      ENDIF

      FTABRC1 = TBRC

      RETURN
      END

