C
C
C
      SUBROUTINE EIRENE_DOPEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,EPS)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 1 ===> Es liegt eine Doppelebenen zur u2,-u3 Ebene vor  *
*     Die Normalform lautet: lambda(1) * u1 ** 2       = 0           *
*     d.h.: nue=0.
*     Aus dieser Gleichung ergibt sich die Loesung:                  *
*     u1 = 0.                                                        *
*     Diese Gleichung wird ruecktransformiert und mit                *
*     Hilfe eines Punktes, der die Normalform erfuellt, kann die     *
*     Ebenengleichung     e1: bo + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE EIRMOD_PRECISION
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3), EPS,
     >                       M(3), NUE
      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3
      REAL(DP) :: NORM, P(3)
C
C     DATA              EPS  / 5.D-10 /
C
*     Loesung der Gleichung:
*     U(1) = 0.0
*     U(2) = 0.0
*     U(3) = 0.0
*
*     Aufstellen des Normalenvektors:
*
      B1 = M(1)
      B2 = M(2)
      B3 = M(3)
*
*     Normieren:
*
      NORM = B1**2 + B2**2 + B3**2
      NORM = SQRT(NORM)
*
      B1 = B1 / NORM
      B2 = B2 / NORM
      B3 = B3 / NORM
*
*     Ein Punkt der Ebene:   P = (0.0,1.,1.)
*     Transformation:
      P(1) = EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,2) + EV(3,3) + M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
      END
