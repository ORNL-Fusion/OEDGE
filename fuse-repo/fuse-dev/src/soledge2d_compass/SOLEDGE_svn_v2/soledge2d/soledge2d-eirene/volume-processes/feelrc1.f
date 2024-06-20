!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  19.12.06: bremsstrahlung added
!dr  31.07.07: bug fix: tein(j) --> tein(k)
 
      FUNCTION EIRENE_FEELRC1 (IRRC,K)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT
      USE EIRMOD_CCONA
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: ELRC1(9), PLS, DELE, EE, EIRENE_FEELRC1,
     .            EIRENE_FTABRC1, DSUB,
     .            DEIMIN, ELRC, EIRENE_ENERGY_RATE_COEFF, BREMS, Z,
     .            eirene_ngffmh
      INTEGER :: J, I, KK, II
      LOGICAL :: LADAS
 
      EIRENE_FEELRC1=0.D0
      KK=NELRRC(IRRC)
 
      IF (KK < 0) THEN
 
        SELECT CASE (KK)
        CASE (-1)
            EIRENE_FEELRC1=-1.5*TEIN(K)*EIRENE_FTABRC1(IRRC,K)
        CASE (-2)
            EIRENE_FEELRC1=EELRC1(IRRC,1)*EIRENE_FTABRC1(IRRC,K)
        CASE (-3)
            EIRENE_FEELRC1=-1.5*TEIN(K)*EIRENE_FTABRC1(IRRC,K)
        END SELECT
 
      ELSE IF (KK > 0) THEN
 
        IF (JELRRC(IRRC) == 1) THEN
          ELRC = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.FALSE.,0)
          ELRC=ELRC+FACREA(KK,2)
          ELRC=EXP(MAX(-100._DP,ELRC))
          EIRENE_FEELRC1=-ELRC*DEIN(K)
        ELSE
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELRC = EIRENE_ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELRC+DEINL(K)+FACREA(KK,2))
          EIRENE_FEELRC1=-EXP(EE)
        END IF
 
        LADAS = EIRENE_IS_RTCEW_ADAS(KK)
        IF (LADAS.AND.(NCHRGP(IPLS) /= 0)) THEN
          Z = NCHRGP(IPLS)
          BREMS = 1.54E-32_DP * TEIN(K)**0.5 * Z**2 *
     .            eirene_ngffmh(Z**2 * 13.6_DP/TEIN(K))*
     .            DEIN(K)*FACREA(KK,1)/ELCHA
          EIRENE_FEELRC1=EIRENE_FEELRC1 + BREMS
        END IF
 
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          EIRENE_FEELRC1=EIRENE_FEELRC1+DELE*EIRENE_FTABRC1(IRRC,K)
        END IF
      END IF
 
      RETURN
      END
