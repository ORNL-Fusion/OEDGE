      MODULE CLMSUR

      USE PRECISION

      IMPLICIT NONE

      PUBLIC

!  data for one particular additional surface MSURF, but without label
!  only used in PLOT.F

      REAL(DP), PUBLIC, SAVE ::
     R A,       B,       C,       D,       E,       F,
     R XL1,     XL2,     YL1,     YL2,
     R AHALB,   BHALB,   XM,      YM,
     R SINA,    COSA,    A0,      A1,      ZL1,     ZL2,     ZPLT,
     R XC1,     XC2,     YC1,     YC2,
     R ALIN(9), XLIN(9), YLIN(9), ZLIN(9),
     R A0S(9),  A1S(9),  A2S(9),  A3S(9),  A4S(9),  A5S(9),
     R A6S(9),  A7S(9),  A8S(9),  A9S(9)

      INTEGER, PUBLIC, SAVE ::
     I MLIN,    MSCN

      END MODULE CLMSUR
