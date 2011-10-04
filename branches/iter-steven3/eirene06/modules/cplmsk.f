      MODULE CPLMSK

      USE PRECISION

      IMPLICIT NONE

      PUBLIC

      REAL(DP), PUBLIC, SAVE ::
     R        X0PL, Y0PL

      REAL(SP), PUBLIC, SAVE ::
     R        LENX,   LENY,   MINX,   MINY,   MAXX,   MAXY,
     R        STPSZX, STPSZY

      INTEGER, PUBLIC, SAVE ::
     I         INTNRX, INTNRY, MINLX,  MINLY,  MAXLX,  MAXLY

      LOGICAL, PUBLIC, SAVE ::
     L         LOGX,   LOGY,   GRIDX,  GRIDY,  FITX,   FITY

      END MODULE CPLMSK
