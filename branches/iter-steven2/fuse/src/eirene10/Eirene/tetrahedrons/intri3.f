

      LOGICAL FUNCTION INTRI3 (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X,Y,Z)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3,
     .                        X, Y, Z
      REAL(DP) :: A1, A2, A3, A, ARTRI3
      LOGICAL :: L1, L2

      A1=ARTRI3(X1,Y1,Z1,X2,Y2,Z2,X,Y,Z)
      A2=ARTRI3(X2,Y2,Z2,X3,Y3,Z3,X,Y,Z)
      A3=ARTRI3(X3,Y3,Z3,X1,Y1,Z1,X,Y,Z)
      A=ARTRI3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

      L1 = MIN(A1,A2,A3).GE.0.D0
      L2 = ABS(A-A1-A2-A3) < 1.D-3
      INTRI3 = L1 .AND. L2

      RETURN
      END
