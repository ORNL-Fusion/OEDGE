!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
 
 
      FUNCTION EIRENE_FTABEI1 (IREI,K)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), EIRENE_FTABEI1, DEIMIN, DSUB, PLS, TBEI,
     .            EIRENE_RATE_COEFF
      REAL(DP) :: ERATE
      INTEGER :: J, I, II, KK
 
      TBEI=0.D0
      KK = NREAEI(IREI)
 
!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
!pb      PLS=MAX(DEIMIN,DEINL(K))-DSUB
      PLS=MAX(DEIMIN,DEINL(K))
 
      TBEI = EIRENE_RATE_COEFF(KK,TEINL(K),PLS,.TRUE.,1,ERATE)*
     .       FACREA(KK,1)
      IF (IFTFLG(KK,2) < 100) TBEI=TBEI*DEIN(K)
 
      EIRENE_FTABEI1 = TBEI
 
      RETURN
      END
 
