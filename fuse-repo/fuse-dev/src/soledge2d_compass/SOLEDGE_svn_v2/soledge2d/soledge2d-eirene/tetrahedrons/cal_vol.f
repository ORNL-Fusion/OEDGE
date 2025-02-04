 
 
      FUNCTION EIRENE_CAL_VOL (P1,P2,P3,P4)
C  RETURNS VOLUME OF TETRAHEDRON DEFINED BY THE 4 POINTS P1,P2,P3,P4
      USE EIRMOD_PRECISION
      USE EIRMOD_CTETRA
      IMPLICIT NONE
 
      REAL(DP) :: EIRENE_CAL_VOL
      REAL(DP), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3)
      REAL(DP) :: X1, X2, X3, X4, Y1, Y2, Y3, Y4, Z1, Z2, Z3, Z4
      REAL(DP) :: X1X2, X1X3, X1X4, Y1Y2, Y1Y3, Y1Y4, Z1Z2, Z1Z3, Z1Z4
 
      X1=P1(1)
      X2=P2(1)
      X3=P3(1)
      X4=P4(1)
      Y1=P1(2)
      Y2=P2(2)
      Y3=P3(2)
      Y4=P4(2)
      Z1=P1(3)
      Z2=P2(3)
      Z3=P3(3)
      Z4=P4(3)
      X1X2=X1-X2
      Y1Y2=Y1-Y2
      Z1Z2=Z1-Z2
      X1X3=X1-X3
      Y1Y3=Y1-Y3
      Z1Z3=Z1-Z3
      X1X4=X1-X4
      Y1Y4=Y1-Y4
      Z1Z4=Z1-Z4
      EIRENE_CAL_VOL = ITETHAND *
     .          (  X1X2 * Y1Y3 * Z1Z4
     .           + X1X4 * Y1Y2 * Z1Z3
     .           + X1X3 * Y1Y4 * Z1Z2
     .           - X1X4 * Y1Y3 * Z1Z2
     .           - X1X2 * Y1Y4 * Z1Z3
     .           - X1X3 * Y1Y2 * Z1Z4 ) / 6.D0
 
      RETURN
      END
