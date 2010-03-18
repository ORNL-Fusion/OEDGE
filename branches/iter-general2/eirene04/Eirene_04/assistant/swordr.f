C
C
      SUBROUTINE SWORDR (X,NX)
C
C  INVERTIEREN DER REIHENFOLGE DER WERTE AUF EINEM ARRAY:
C  X(1)-->X(NX)
C  X(2)-->X(NX-1)
C       .
C       .
C       .
C
      USE PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NX
      REAL(DP), INTENT(INOUT) :: X(NX)
      REAL(DP) :: H
      INTEGER :: IS, I
      
      IS=NX/2
      DO 1 I=1,IS
      H=X(I)
      X(I)=X(NX-I+1)
1     X(NX-I+1)=H
      RETURN
      END
