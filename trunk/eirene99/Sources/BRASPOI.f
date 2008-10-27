      MODULE BRASPOI

      PRIVATE
      PUBLIC :: CELLSIM, CELLMUL, SIMARR, MULARR,
     .          EAELS, EMELS, EIELS,
     .          EAPLS, EMPLS, EIPLS, PAPLS, PMPLS, PIPLS,
     .          PDENAS, PDENMS, PDENIS, EDENAS, COPVS

      INCLUDE 'PARMUSR'
      TYPE :: CELLSIM
        REAL*8 :: VALUES
        INTEGER :: ICS
        TYPE(CELLSIM), POINTER :: NXTSIM
      END TYPE CELLSIM

      TYPE :: CELLMUL
        REAL*8 :: VALUEM
        INTEGER :: IART,ICM
        TYPE(CELLMUL), POINTER :: NXTMUL
      END TYPE CELLMUL

      TYPE :: SIMARR
        TYPE(CELLSIM),POINTER :: PSIM
      END TYPE SIMARR

      TYPE :: MULARR
        TYPE(CELLMUL),POINTER :: PMUL
      END TYPE MULARR

      TYPE(SIMARR), DIMENSION(NSTRA),SAVE :: EAELS, EMELS, EIELS,
     .                                       EAPLS, EMPLS, EIPLS
      TYPE(MULARR), DIMENSION(NSTRA),SAVE :: PAPLS, PMPLS, PIPLS,
     .                                       PDENAS, PDENMS, PDENIS,
     .                                       EDENAS, COPVS

      END MODULE BRASPOI
