      MODULE EIRMOD_CRAND
 
      USE EIRMOD_PRECISION
 
      IMPLICIT NONE
 
      PUBLIC
 
      REAL(DP), DIMENSION(64), PUBLIC, SAVE ::
     R FG1, FG2, FG3,
     R FM1, FM2, FM3,
     R FC1, FC2, FC3,
     R FI1, FI2, FI3
 
      INTEGER, PUBLIC, SAVE ::
     I INIV1,  INIV2,  INIV3,  INIV4,  IRNDVC, IRNDVH,
     I INTMAX, ISEEDR
 
      END MODULE EIRMOD_CRAND
