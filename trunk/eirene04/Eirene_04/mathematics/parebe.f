C
C
C
      SUBROUTINE PAREBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3,EPS)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 2 ===> Es liegen 2 parallele Ebenen zur u2,-u3 Ebene vor*
*     Die Normalform lautet: lambda(1) * u1 ** 2 + nue = 0           *
*     Aus dieser Gleichung ergeben sich die Loesungen:               *
*     u1 = sqrt ( -nue/lambda(1))   und  u1 = -sqrt(-nue/lambda(1)   *
*     Diese beiden Gleichungen werden ruecktransformiert und mit     *
*     Hilfe eines Punktes, der die Normalform erfuellt koennen die   *
*     Ebenengleichungen   e1: bo + b1*x + b2 * y + b3 * z            *
*                         e2: co + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3, C0, C1, C2, C3
      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3), EPS,
     >                      M(3), NUE
      REAL(DP) :: NORM, P(3), U(3)

C
C     DATA              EPS  / 5.D-10 /
C
*     1. Loesung der Gleichung:
*
      U(1) = SQRT(-NUE/LAMBDA(1))
      U(2) = 0.0
      U(3) = 0.0
*
*     Aufstellen des Normalenvektors:'
*
      B1 = EV(1,1) * U(1)
      B2 = EV(2,1) * U(1)
      B3 = EV(3,1) * U(1)
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
*     Ein Punkt der Ebene:   P = (U(1),1.,1.)
*     Transformation:
      P(1) = EV(1,1) * U(1) + EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,1) * U(1) + EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,1) * U(1) + EV(3,2) + EV(3,3) + M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
*
*     2. Loesung der Gleichung:
*
      U(1) = -SQRT(-NUE/LAMBDA(1))
      U(2) = 0.0
      U(3) = 0.0
*
*     Aufstellen des Normalenvektors:'
*
      C1 = EV(1,1) * U(1)
      C2 = EV(2,1) * U(1)
      C3 = EV(3,1) * U(1)
*
*     Normieren:
*
      NORM = C1**2 + C2**2 + C3**2
      NORM = SQRT(NORM)
*
      C1 = C1 / NORM
      C2 = C2 / NORM
      C3 = C3 / NORM
*
*     Ein Punkt der Ebene:   P = (U(1),1.,1.)
*     Transformation:
      P(1) = EV(1,1) * U(1) + EV(1,2) + EV(1,3) + M(1)
      P(2) = EV(2,1) * U(1) + EV(2,2) + EV(2,3) + M(2)
      P(3) = EV(3,1) * U(1) + EV(3,2) + EV(3,3) + M(3)
*
      C0 = -(C1 * P(1) + C2 * P(2) + C3 * P(3))
*
      END
