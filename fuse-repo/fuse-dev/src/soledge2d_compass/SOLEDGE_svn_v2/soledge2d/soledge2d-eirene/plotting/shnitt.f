C
C
      SUBROUTINE EIRENE_SHNITT(P,PX,PY,PZ,VX,VY,VZ,A,I,XP,YP,JA,JE,IXS)
C
C  CALLED FROM ZYLIND
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: P(3,*), A(*)
      REAL(DP), INTENT(OUT) :: XP(*), YP(*)
      REAL(DP), INTENT(IN) :: PX, PY, PZ, VX, VY, VZ
      INTEGER, INTENT(IN) :: I, JA, JE, IXS
 
      REAL(DP) :: EPS12, XX, YY, ZZ, XN, ALAMDA, ALAM2, DET, ROOT,
     .            ALAM1, AA, BB, CC
      INTEGER :: J, IX
      LOGICAL :: LERR
 
      DATA EPS12 /1.E-12_DP/
 
      LERR=.FALSE.
      IX=IXS-1
C     WRITE (iunout,*) ' SHNITT  I = ',I
      IF (I.LE.4) THEN
C  SCHNITTKURVE MIT EBENE A1+A2X+A3Y+A4Z=0
C       WRITE (iunout,*) ' VX,VY,VZ ',VX,VY,VZ
C       WRITE (iunout,*) ' A ',A(1),A(2),A(3),A(4)
        XN=VX*A(2)+VY*A(3)+VZ*A(4)
C       WRITE (iunout,*) ' XN ',XN
        DO 100 J=JA,JE
          XX=P(1,J)+PX
          YY=P(2,J)+PY
          ZZ=P(3,J)+PZ
          IF (LERR) THEN
            ALAMDA=0.
            GOTO 101
          ENDIF
C  BERECHNE SCHNITTPUNKT (ALAMDA) VON XX+ALAMDA*VX, YY+...,ZZ+...
C  MIT DER EBENE. ALAMDA MUSS POSITIV SEIN, SONST FALSCHE EINGABE
          IF (ABS(XN).LT.EPS12) THEN
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
            WRITE (iunout,*) 'NO INTERSECTION FOUND WITH PLANE'
            ALAMDA=0.
            LERR=.TRUE.
          ELSE
            ALAMDA=(-A(1)-(A(2)*XX+A(3)*YY+A(4)*ZZ))/XN
            IF (ALAMDA.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT. SET ALAMDA=0.'
              WRITE (iunout,*) 'NO INTERSECTION IN POSITIV DIRECTION'
              WRITE (iunout,*) 'WITH PLANE '
              ALAMDA=0.
              LERR=.TRUE.
            ENDIF
          ENDIF
C
101       XX=XX+ALAMDA*VX
          YY=YY+ALAMDA*VY
          ZZ=ZZ+ALAMDA*VZ
          IX=IX+1
          CALL EIRENE_PL3D(XX,YY,ZZ,XP(IX),YP(IX))
100     CONTINUE
      ELSEIF (I.GT.4) THEN
C  SCHNITTKURVE MIT VOLLER GLEICHUNG 2TER ORDNUNG
        AA=(A(5)*VX+A(8)*VY+A(9)*VZ)*VX+(A(6)*VY+A(10)*VZ)*VY+
     .      A(7)*VZ*VZ
        DO 200 J=JA,JE
          XX=P(1,J)+PX
          YY=P(2,J)+PY
          ZZ=P(3,J)+PZ
          IF (LERR) THEN
            ALAMDA=0.
            GOTO 201
          ENDIF
          BB=(A(2)+2.*A(5)*XX+A(8)*YY+A(9)*ZZ)*VX+
     .       (A(3)+2.*A(6)*YY+A(8)*XX+A(10)*ZZ)*VY+
     .       (A(4)+2.*A(7)*ZZ+A(9)*XX+A(10)*YY)*VZ
          CC=A(1)+(A(2)+A(5)*XX+A(8)*YY+A(9)*ZZ)*XX+
     .       (A(3)+A(6)*YY+A(10)*ZZ)*YY+(A(4)+A(7)*ZZ)*ZZ
C
          IF (ABS(AA).LT.EPS12.AND.ABS(BB).LT.EPS12) THEN
            LERR=.TRUE.
            ALAMDA=0.
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
            WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
            WRITE (iunout,*) 'AA=BB=0.'
            GOTO 201
          ELSEIF (ABS(AA).LT.EPS12.AND.ABS(BB).GT.EPS12) THEN
            ALAM1=-CC/BB
            ALAM2=-CC/BB
          ELSEIF (ABS(BB).LT.EPS12.AND.ABS(AA).GT.EPS12) THEN
            DET=-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 201
            ENDIF
            ALAM1=SQRT(DET)
            ALAM2=-SQRT(DET)
          ELSE IF (ABS(CC).LT.EPS12) THEN
            ALAM1=0.
            ALAM2=-BB/AA
          ELSE
            DET=BB*BB/(4*AA*AA)-CC/AA
            IF (DET.LT.0.) THEN
              WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0. '
              WRITE (iunout,*) 'NO INTERSECTION WITH SURFACE FOUND '
              WRITE (iunout,*) 'DETERMINANT DET= ',DET
              ALAMDA=0.
              LERR=.TRUE.
              GOTO 201
            ENDIF
            ROOT=SQRT(DET)
            ALAM1=-BB/(2.*AA)+ROOT
            ALAM2=-BB/(2.*AA)-ROOT
          ENDIF
          IF (LERR) GOTO 201
C  DECIDE, WHICH ONE OF THE 2 SOLUTIONS TO TAKE
          IF (ALAM1*ALAM2.LT.0.) THEN
            ALAMDA=MAX(ALAM1,ALAM2)
          ELSE IF (ABS(ALAM1).GT.ABS(ALAM2)) THEN
            ALAMDA=ALAM2
          ELSE
            ALAMDA=ALAM1
          ENDIF
          IF (ALAMDA.LT.0.) THEN
            WRITE (iunout,*) 'ERROR IN SUBR. SHNITT, SET ALAMDA=0.'
            WRITE (iunout,*) 'INTERSECTION IN WRONG DIRECTION'
            ALAMDA=0.
            LERR=.TRUE.
            GOTO 201
          ENDIF
201       XX=XX+ALAMDA*VX
          YY=YY+ALAMDA*VY
          ZZ=ZZ+ALAMDA*VZ
          IX=IX+1
          CALL EIRENE_PL3D(XX,YY,ZZ,XP(IX),YP(IX))
200     CONTINUE
      ENDIF
      RETURN
      END
