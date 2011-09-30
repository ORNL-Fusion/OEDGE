C
C
      FUNCTION PARA2U (XIN)
      USE CLMSUR
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: XIN
      REAL(DP) :: RAD, Y, DH, PARA2U
      DH=-E/(2.*C)
      RAD=(E*E)/(4.*C*C)-(D*XIN+F)/C
      IF (RAD.GE.0.) GOTO 1
      PARA2U=1.D50
      RETURN
1     Y=-SQRT(RAD)
      PARA2U=Y-DH
      RETURN
      END
