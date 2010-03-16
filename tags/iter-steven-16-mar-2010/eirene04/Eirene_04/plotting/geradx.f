C
C
      FUNCTION GERADX(Y)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(INOUT) :: Y
      REAL(DP) :: GERADX
      GERADX=Y
      Y=A0*Y+A1
      RETURN
      END
