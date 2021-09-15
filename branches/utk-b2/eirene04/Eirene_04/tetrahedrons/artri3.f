

      FUNCTION ARTRI3 (X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3
      REAL(DP) :: A1, A2, A3, ARTRI3

      A1 = Y1*Z2 + Y3*Z1 + Y2*Z3 - Y3*Z2 - Y1*Z3 - Y2*Z1 
      A2 = Z1*X2 + Z3*X1 + Z2*X3 - Z3*X2 - Z1*X3 - Z2*X1
      A3 = X1*Y2 + X3*Y1 + X2*Y3 - X3*Y2 - X1*Y3 - X2*Y1 

      ARTRI3=0.5_DP * SQRT(A1**2+A2**2+A3**2)
      RETURN
      END
