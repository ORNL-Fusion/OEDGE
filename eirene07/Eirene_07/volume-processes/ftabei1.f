!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced


      FUNCTION FTABEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), FTABEI1, DEIMIN, DSUB, PLS, TBEI, RATE_COEFF
      INTEGER :: J, I, II, KK

      TBEI=0.D0
      KK = NREAEI(IREI)

!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
!pb      PLS=MAX(DEIMIN,DEINL(K))-DSUB
      PLS=MAX(DEIMIN,DEINL(K))

      TBEI = RATE_COEFF(KK,TEINL(K),PLS,.TRUE.,1)*FACREA(KK)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)

      FTABEI1 = TBEI

      RETURN
      END

