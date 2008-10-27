C
C
      SUBROUTINE SPOINT (A0,A1,A2,A3,X,Y,Z,XR,YR,ZR,P,T)
C                                                     X+ LAMBDA * XR
C  SCHNITTPUNKT DER EBENE A0,A1,A2,A3 MIT DER GERADEN Y+ LAMBDA * YR
C                                                     Z+ LAMBDA * ZR
C  OUTPUT:  P(3): SCHNITTPUNKT
C           T:    LAMBDA
C
      USE PRECISION
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A0, A1, A2, A3, X, Y, Z, XR, YR, ZR
      REAL(DP), INTENT(OUT) :: P(3), T
      REAL(DP) :: XNEN

      XNEN=A1*XR+A2*YR+A3*ZR
      IF (ABS(XNEN).LT.1.D-30) THEN
        P(1)=1.D60
        P(2)=1.D60
        P(3)=1.D60
        T=1.D60
      ELSE
        T=-(A0+A1*X+A2*Y+A3*Z)/XNEN
        P(1)=X+T*XR
        P(2)=Y+T*YR
        P(3)=Z+T*ZR
      ENDIF
      RETURN
      END
