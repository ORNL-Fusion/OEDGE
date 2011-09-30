!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FTABRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: DEIMIN, DSUB, FTABRC1, ZX, TBRC, PL, RATE_COEFF, ERATE
      INTEGER :: II, KK

      TBRC=0.D0
      KK = NREARC(IRRC)

      IF (KK == 0) THEN
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)

      ELSE

!pb        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
!pb        PL=MAX(DEIMIN,DEINL(K))-DSUB
        PL=MAX(DEIMIN,DEINL(K))

        TBRC = RATE_COEFF(KK,TEINL(K),PL,.TRUE.,1,ERATE)*FACREA(KK,1)
        IF (IFTFLG(KK,2) < 100) TBRC=TBRC*DEIN(K)
        
      ENDIF

      FTABRC1 = TBRC

      RETURN
      END

