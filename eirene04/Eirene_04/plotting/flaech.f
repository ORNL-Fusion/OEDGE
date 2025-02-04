

      FUNCTION FLAECH (X1,Y1,X2,Y2,X3,Y3,X4,Y4)

      USE PRECISION

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1, Y1, X2, Y2, X3, Y3, X4, Y4
      REAL(DP) :: FLAECH1, FLAECH2, FLAECH

      FLAECH1 = 0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
      FLAECH2 = 0.5*(X1*(Y3-Y4)+X3*(Y4-Y1)+X4*(Y1-Y3))
      FLAECH = ABS(FLAECH1) + ABS(FLAECH2)
      END
