c  0707: new, for PI processes, copied and adapted from feelei1.f

      FUNCTION FEELPI1 (IRPI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: ELPIC(9), FEELPI1, PLS, DEIMIN, EE, FTABPI3,
     .            ELPI, ENERGY_RATE_COEFF, DELE
      INTEGER :: J, I, KK, II

      FEELPI1=0.D0
      KK=NELRPI(IRPI)
      IF (KK < 0) THEN
        IF (KK == -1) THEN
          FEELPI1 = EPLPI3(IRPI,1,1) 
        ELSE
          WRITE (IUNOUT,*) 'ERROR IN FEELPI1: KK IS LT 0. '
          WRITE (IUNOUT,*) 
     .       'BUT THERE SHOULD BE NO DEFAULT PI PROCESSES '
          CALL EXIT_OWN(1)
        END IF
      ELSE IF (KK > 0) THEN
        IF (JELRPI(IRPI) == 1) THEN
          ELPI = ENERGY_RATE_COEFF(KK,TEINL(K),0._DP,.TRUE.,0)
          FEELPI1=-ELPI*DEIN(K)*FACREA(KK,1)/(FTABPI3(IRPI,K)+EPS60)
        ELSEIF(JELRPI(IRPI) == 9) THEN
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))
          ELPI = ENERGY_RATE_COEFF(KK,TEINL(K),PLS,.FALSE.,1)
          EE=MAX(-100._DP,ELPI+FACREA(KK,2)+DEINL(K))
          FEELPI1=-EXP(EE)/(FTABPI3(IRPI,K)+EPS60)
        ELSE
          WRITE (IUNOUT,* ) 'ERROR IN FEELPI1, INVALID JELRPI '
          CALL EXIT_OWN(1)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          FEELPI1=FEELPI1+DELE
        END IF
      END IF

      RETURN
      END
