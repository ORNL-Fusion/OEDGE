C
C
      SUBROUTINE SCCONE(X0,Y0,Z0,VX,VY,VZ,ALF,T1,T2,BX,BY,BZ,CX,CY,CZ,
     .                  DANG,A,I,XP,YP,JA,JE,IXS)
C
C  CALLED FROM CONE
C
      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(*)
      REAL(DP), INTENT(OUT) :: XP(*),YP(*)
      REAL(DP), INTENT(IN) :: X0, Y0, Z0, VX, VY, VZ, ALF, T1, T2,
     .                      BX, BY, BZ, CX, CY, CZ, DANG
      INTEGER, INTENT(IN) :: I, JA, JE, IXS
      REAL(DP) :: EPS12, XK2, YK2, PXX1, PYY1, PZZ1, PXX2, PYY2, PZZ2,
     .          PX2, PY2, PZ2, XK1, YK1, PHI, DET, ALAM1, ALAM2, ROOT,
     .          XX, YY, ZZ, CC, PX21, PY21, PZ21, ALAMDA, AA, BB, XN,
     .          RAD1, RAD2, PX1, PY1, PZ1
      INTEGER :: J, IX
      LOGICAL LERR

      DATA EPS12 /1.E-12/

      LERR=.FALSE.
      IX=IXS-1
      RAD1=T1*TAN(ALF)
      RAD2=T2*TAN(ALF)
      PX1=X0+T1*VX
      PY1=Y0+T1*VY
      PZ1=Z0+T1*VZ
      PX2=X0+T2*VX
      PY2=Y0+T2*VY
      PZ2=Z0+T2*VZ
      DO 100 J=JA,JE
        PHI=(J-1)*DANG
        XK1=RAD1*COS(PHI)
        YK1=RAD1*SIN(PHI)
        PXX1=XK1*BX+YK1*CX
        PYY1=XK1*BY+YK1*CY
        PZZ1=XK1*BZ+YK1*CZ
        XK2=RAD2*COS(PHI)
        YK2=RAD2*SIN(PHI)
        PXX2=XK2*BX+YK2*CX
        PYY2=XK2*BY+YK2*CY
        PZZ2=XK2*BZ+YK2*CZ
        PX21=PXX2-PXX1
        PY21=PYY2-PYY1
        PZ21=PZZ2-PZZ1
        IF (LERR) THEN
          ALAMDA=0.
          GOTO 101
        ENDIF
        IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
          XN=PX21*A(2)+PY21*A(3)+PZ21*A(4)
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
          IF (ABS(XN).LT.EPS12) THEN
            WRITE (6,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
            WRITE (6,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*PXX1+A(3)*PYY1+A(4)*PZZ1))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE. SET ALAMDA=0.'
              WRITE (6,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (6,*) 'WITH PLANE '
              ALAMDA=0.
              LERR=.TRUE.
            ENDIF
          ENDIF
C
        ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
          AA=(A(5)*PX21+A(8)*PY21+A(9)*PZ21)*PX21+
     .       (A(6)*PY21+A(10)*PZ21)*PY21+A(7)*PZ21*PZ21
          BB=(A(2)+2.*A(5)*PXX1+A(8)*PYY1+A(9)*PZZ1)*PX21+
     .       (A(3)+2.*A(6)*PYY1+A(8)*PXX1+A(10)*PZZ1)*PY21+
     .       (A(4)+2.*A(7)*PZZ1+A(9)*PXX1+A(10)*PYY1)*PZ21
          CC=A(1)+(A(2)+A(5)*PXX1+A(8)*PYY1+A(9)*PZZ1)*PXX1+
     .       (A(3)+A(6)*PYY1+A(10)*PZZ1)*PYY1+(A(4)+A(7)*PZZ1)*PZZ1
C
          IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
            LERR=.TRUE.
            ALAMDA=0.
            WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'AA=BB=0.'
            GOTO 101
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 101
            ENDIF
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ELSE IF (ABS(CC).LT.EPS12) THEN
            ALAM1=0.
            ALAM2=-BB/AA
          ELSE
            DET=BB*BB/(4*AA*AA)-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (6,*) 'ERROR IN SUBR. SCCONE, SET ALAMDA=0. '
              WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (6,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 101
            ENDIF
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
          IF (LERR) GOTO 101
C  DECIDE, WHICH ONE OF THE 2 SOLUTIONS TO TAKE
          IF (ALAM1*ALAM2.LT.0.) THEN
            ALAMDA=MAX(ALAM1,ALAM2)
          ELSE IF (ABS(ALAM1).GT.ABS(ALAM2)) THEN
            ALAMDA=ALAM2
          ELSE
            ALAMDA=ALAM1
          ENDIF
          IF (ALAMDA.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (6,*) 'INTERSECTION IN WRONG DIRECTION'
            ALAMDA=0.
            LERR=.TRUE.
            GOTO 101
          ENDIF
        ENDIF
C
101     XX=PXX1+ALAMDA*PX21
        YY=PYY1+ALAMDA*PY21
        ZZ=PZZ1+ALAMDA*PZ21
        IX=IX+1
        CALL PL3D(XX,YY,ZZ,XP(IX),YP(IX))
100   CONTINUE
      RETURN
      END
