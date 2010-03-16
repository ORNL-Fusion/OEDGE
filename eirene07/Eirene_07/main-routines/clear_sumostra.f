c  19.12.05:  stvw, stvws = 0   included
      SUBROUTINE CLEAR_SUMOSTRA

      USE PRECISION
      USE PARMMOD
      USE CSPEI
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMSOU
      USE CESTIM

      IMPLICIT NONE

      INTEGER :: ISPC

C
C**** CLEAR WORK AREA FOR SUM OVER STRATA ****************************
C
      SMESTV = 0._DP
      SMESTS = 0._DP
      STV    = 0._DP
      STVS   = 0._DP
      STVW   = 0._DP
      STVWS  = 0._DP
      STVC   = 0._DP
      STVCS  = 0._DP
C  ARRAYS: EE,FF,....
      EE     = 0._DP
      EES    = 0._DP
      FF     = 0._DP
      FFS    = 0._DP
C  BGK-ARRAYS:
      IF (NSIGI_BGK.GT.0) THEN
        STVS_BGK=0._DP
        EES_BGK=0._DP
        STV_BGK=0._DP
        EE_BGK=0._DP
      ENDIF
C  COP-ARRAYS:
      IF (NSIGI_COP.GT.0) THEN
        STVS_COP=0._DP
        EES_COP=0._DP
        STV_COP=0._DP
        EE_COP=0._DP
      ENDIF

C  SPECTRA
      IF ((NSTRAI > 1) .AND. (NSMSTRA > 0)) THEN
        DO ISPC=1,NADSPC
          SMESTL(ISPC)%PSPC%SPC = 0._DP
          IF (NSIGI_SPC > 0) THEN
            SMESTL(ISPC)%PSPC%SGM = 0._DP
            SMESTL(ISPC)%PSPC%SDV = 0._DP
          END IF
          SMESTL(ISPC)%PSPC%SPCINT = 0._DP
          SMESTL(ISPC)%PSPC%SGMS = 0._DP
          SMESTL(ISPC)%PSPC%STVS = 0._DP
          SMESTL(ISPC)%PSPC%EES = 0._DP
        END DO
      END IF

      RETURN

      END SUBROUTINE CLEAR_SUMOSTRA
