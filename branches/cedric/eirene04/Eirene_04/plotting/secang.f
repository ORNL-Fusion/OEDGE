C
C
      SUBROUTINE SECANG (B0,B1,B2,RAD,ANG1,ANG2)

      USE PRECISION
      USE CCONA
      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: B0, B1, B2, RAD
      REAL(DP), INTENT(OUT) :: ANG1, ANG2
      REAL(DP) :: B0Q, B1Q, B2Q, X2, Y2, B0B1, RD, XNEN, X1, Y1
C
      ANG1=0.
      ANG2=360.
C
      IF (ABS(B1).LT.EPS10.AND.ABS(B2).LT.EPS10) RETURN
C
      IF (ABS(B1).LT.EPS10) THEN
        Y1=-B0/B2
        IF (ABS(Y1).GE.RAD) RETURN
        X1=SQRT(RAD*RAD-Y1*Y1)
        X2=-X1
        Y2=Y1
      ELSEIF (ABS(B2).LT.EPS10) THEN
        X1=-B0/B1
        IF (ABS(X1).GE.RAD) RETURN
        Y1=SQRT(RAD*RAD-X1*X1)
        X2=X1
        Y2=-Y1
      ELSE
        B1Q=B1*B1
        B2Q=B2*B2
        B0Q=B0*B0
        XNEN=B1Q+B2Q
        B0B1=B0*B1
        RD=B0B1*B0B1/XNEN/XNEN-(B0Q-B2Q*RAD*RAD)/XNEN
        IF (RD.LT.0.) RETURN
        X1=-B0B1/XNEN+SQRT(RD)
        X2=-B0B1/XNEN-SQRT(RD)
        Y1=-B0/B2-B1/B2*X1
        Y2=-B0/B2-B1/B2*X2
      ENDIF
C
      ANG1=ATAN2(Y1,X1)*RADDEG
      IF (ANG1.LT.0.) ANG1=ANG1+360.
      ANG2=ATAN2(Y2,X2)*RADDEG
      IF (ANG2.LT.0.) ANG2=ANG2+360.
C
      RETURN
      END