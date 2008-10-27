C
C
      FUNCTION HYPM(X)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X
      REAL(DP) :: HYPM, DH, RAD, AH
      AH=A*X*X
      DH=D*X
      RAD=E*E/(4.*C*C)-(F+AH+DH)/C
      IF (RAD.GE.0.) GOTO 1
      HYPM=1.D50
      RETURN
1     HYPM=-E/(2.*C)-SQRT(RAD)
      RETURN
      END
