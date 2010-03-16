C
C
      FUNCTION PARA1 (X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: PARA1
      PARA1=-A/E*X*X-D/E*X-F/E
      RETURN
      END
