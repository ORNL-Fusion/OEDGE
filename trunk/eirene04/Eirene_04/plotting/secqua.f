C
C
      SUBROUTINE SECQUA(XX,YY,ZZ,VX,VY,VZ,A,I,ALAM1,ALAM2,LERR)
C
C MORE GENERAL THAN SHNITT OR CSCONE, CALLED FROM PL3D ITSELF
C
      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: A(*)
      REAL(DP), INTENT(IN) :: XX, YY, ZZ, VX, VY, VZ
      REAL(DP), INTENT(OUT) :: ALAM1, ALAM2
      INTEGER, INTENT(IN) :: I
      LOGICAL, INTENT(OUT) :: LERR
      REAL(DP) :: EPS12, XN, AA, DET, ROOT, BB, CC

      DATA EPS12 /1.E-12/

      LERR=.FALSE.
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
        IF (ABS(XN).LT.EPS12) THEN
          WRITE (6,*) 'ERROR IN SUBR. SECQUA. SET ALAM1=ALAM2=0.'
          WRITE (6,*) 'NO INTERSECTION FOUND WITH PLANE'
          ALAM1=0.
          ALAM2=0.
          LERR=.TRUE.
        ELSE
          ALAM1=(-A(1)-(A(2)*XX+A(3)*YY+A(4)*ZZ))/XN
          ALAM2=ALAM1
        ENDIF
C
      ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
        AA=(A(5)*VX+A(8)*VY+A(9)*VZ)*VX+(A(6)*VY+A(10)*VZ)*VY+
     .      A(7)*VZ*VZ
        BB=(A(2)+2.*A(5)*XX+A(8)*YY+A(9)*ZZ)*VX+
     .     (A(3)+2.*A(6)*YY+A(8)*XX+A(10)*ZZ)*VY+
     .     (A(4)+2.*A(7)*ZZ+A(9)*XX+A(10)*YY)*VZ
        CC=A(1)+(A(2)+A(5)*XX+A(8)*YY+A(9)*ZZ)*XX+
     .     (A(3)+A(6)*YY+A(10)*ZZ)*YY+(A(4)+A(7)*ZZ)*ZZ
C
        IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
          LERR=.TRUE.
          ALAM1=0.
          ALAM2=0.
          WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
          WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
          WRITE (6,*) 'AA=BB=0.'
        ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
          ALAM1=-CC/BB
          ALAM2=-CC/BB
        ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
          DET=-CC/AA
          IF (DET.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'DETERMINANT DET= ',DET
            ALAM1=0.
            ALAM2=0.
            LERR=.TRUE.
          ELSE
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ENDIF
        ELSE IF (ABS(CC).LT.EPS12) THEN
          ALAM1=0.
          ALAM2=-BB/AA
        ELSE
          DET=BB*BB/(4*AA*AA)-CC/AA
          IF (DET.LT.0.) THEN
            WRITE (6,*) 'ERROR IN SUBR. SECQUA, SET ALAM1=ALAM2=0. '
            WRITE (6,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (6,*) 'DETERMINANT DET= ',DET
            ALAM1=0.
            ALAM2=0.
            LERR=.TRUE.
          ELSE
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
        ENDIF
      ENDIF
      RETURN
      END
