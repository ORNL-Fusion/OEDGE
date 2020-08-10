C
C
C
      SUBROUTINE SCHEBE(EV,LAMBDA,NUE,M,B0,B1,B2,B3,C0,C1,C2,C3)
**********************************************************************
*                                                     1. JUNI 1988   *
*     Inde = 3 ===> Es liegen 2 sich schneidende Ebenen vor.         *
*     Schnittgerade: u3 - Achse                                      *
*     Die Normalform lautet:                                         *
*      lambda(1) * u1 ** 2 + lambda(2) * u2**2  = 0                  *
*     d.h. nue=0.
*     Aus dieser Gleichung ergeben sich die Loesungen:               *
*     u1 = sqrt ( -lambda(2)/lambda(1)) * u2                         *
*     u1 = - sqrt ( -lambda(2)/lambda(1)) * u2                       *
*     Fuer jede dieser Gleichung werden Normalenvektor und ein Punkt *
*     der Ebene ruecktransformiert, so dass die Ebenengleichungen    *
*                         e1: bo + b1*x + b2 * y + b3 * z            *
*                         e2: co + b1*x + b2 * y + b3 * z            *
*     aufgestellt werden.                                            *
**********************************************************************
*
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: LAMBDA( 3 ), EV(3,3),
     >                      M(3), NUE
      REAL(DP), INTENT(OUT) :: B0, B1, B2, B3, C0, C1, C2, C3
      REAL(DP) :: NORM, P(3), U(3)

*
*
*     1. Loesung der Gleichung:
*     1. Normalenvektor
      U(1) = 1.0
      U(2) = -SQRT(-LAMBDA(2)/LAMBDA(1))
      U(3) = 0.0
*
      NORM = U(1)**2 + U(2)**2
      NORM = SQRT(NORM)
*
      U(1) = U(1) / NORM
      U(2) = U(2) / NORM
*
*     Ruecktransformation des Normalenvektors:'
*
      B1 = EV(1,1) * U(1) + EV(1,2) * U(2)
      B2 = EV(2,1) * U(1) + EV(2,2) * U(2)
      B3 = EV(3,1) * U(1) + EV(2,3) * U(2)
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
*     Ein Punkt der Ebene:   P = (sqrt(-lambda(2)/lambda(1),1.,0.)
*     Ruecktransformation:
      P(1) = EV(1,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(1,2)+M(1)
      P(2) = EV(2,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(2,2)+M(2)
      P(3) = EV(3,1)*SQRT(-LAMBDA(2)/LAMBDA(1)) + EV(3,2)+M(3)
*
      B0 = -(B1 * P(1) + B2 * P(2) + B3 * P(3))
*
*
*     2. Loesung der Gleichung:
*     Normalenvektor:
      U(1) = 1.0
      U(2) = SQRT(-LAMBDA(2)/LAMBDA(1))
      U(3) = 0.0

      NORM = U(1)**2 + U(2)**2
      NORM = SQRT(NORM)
*
      U(1) = U(1) / NORM
      U(2) = U(2) / NORM
*
*     Ruecktransformation des Normalenvektors:'
*
      C1 = EV(1,1) * U(1) + EV(1,2) * U(2)
      C2 = EV(2,1) * U(1) + EV(2,2) * U(2)
      C3 = EV(3,1) * U(1) + EV(3,2) * U(2)
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
*     Ein Punkt der Ebene:   P = (-SQRT(-LAMBDA(2)/LAMBDA(1)),1.,0.)
*     Transformation:
      P(1) = -EV(1,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(1,2) + M(1)
      P(2) = -EV(2,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(2,2) + M(2)
      P(3) = -EV(3,1)*SQRT(-LAMBDA(2)/LAMBDA(1))+EV(3,2) + M(3)
*
      C0 = -(C1 * P(1) + C2 * P(2) + C3 * P(3))
      END
