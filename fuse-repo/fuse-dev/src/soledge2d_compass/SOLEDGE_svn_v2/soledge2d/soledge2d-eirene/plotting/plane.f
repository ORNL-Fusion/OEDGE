C
C
      SUBROUTINE EIRENE_PLANE (A0,A1,A2,A3,RL,N1,EPS,
     .  AL,XL,YL,ZL,AL0,XL1,YL1,ZL1,XL2,YL2,ZL2,XL3,YL3,ZL3,IO,NF,NUM)
C
C  PLOT PLANE, SECTION INSIDE A BOX
C
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: A0, A1, A2, A3, RL, EPS
      INTEGER, INTENT(IN) :: N1, IO, NUM
      REAL(DP), INTENT(IN) :: AL(N1,*),  XL(N1,*),  YL(N1,*),  ZL(N1,*),
     .                      AL0(N1,*), XL1(N1,*), YL1(N1,*), ZL1(N1,*),
     .                      XL2(N1,*), YL2(N1,*), ZL2(N1,*),
     .                      XL3(N1,*), YL3(N1,*), ZL3(N1,*)
      LOGICAL NF
      REAL(DP) :: P(3,36), XYZG(3,36), ANGLE(36), CORD(108), PS(3)
      REAL(DP) :: B1, B2, B3, DET, EIRENE_DETER, TEST, T, DX, DY, DZ, 
     .            HELP, XMIT, YMIT, ZMIT, PI, ANG, X1, X2, Y1, Y2, 
     .            Z1, Z2
      INTEGER :: K, J, IPOINT, IP, ILN, I, II, ISORT, ICOUNT, ICHECK,
     .           IS
C
      IS=0
      IF (RL.GT.0) THEN
        X1=XL1(1,NUM)
        X2=XL2(1,NUM)
        Y1=YL1(1,NUM)
        Y2=YL2(1,NUM)
        Z1=ZL1(1,NUM)
        Z2=ZL2(1,NUM)
        DX=XL2(1,NUM)-XL1(1,NUM)
        DY=YL2(1,NUM)-YL1(1,NUM)
        DZ=ZL2(1,NUM)-ZL1(1,NUM)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y1,Z1,DX,0._DP,0._DP,P(1,1),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y1,Z1,0._DP,0._DP,DZ,P(1,2),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y1,Z2,DX,0._DP,0._DP,P(1,3),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X2,Y1,Z2,0._DP,0._DP,DZ,P(1,4),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X2,Y1,Z1,0._DP,DY,0._DP,P(1,5),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y2,Z1,DX,0._DP,0._DP,P(1,6),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y2,Z1,0._DP,0._DP,DZ,P(1,7),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y2,Z2,DX,0._DP,0._DP,P(1,8),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X2,Y2,Z2,0._DP,0._DP,DZ,P(1,9),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X2,Y1,Z2,0._DP,DY,0._DP,P(1,10),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y1,Z1,0._DP,DY,0._DP,P(1,11),T)
        CALL EIRENE_SPOINT
     .  (A0,A1,A2,A3,X1,Y1,Z2,0._DP,DY,0._DP,P(1,12),T)
        IPOINT=12
        DO 10 I=1,IPOINT
          IF ((P(1,I)-X1.GE.-EPS.AND.X2-P(1,I).GE.-EPS).AND.
     .       (P(2,I)-Y1.GE.-EPS.AND.Y2-P(2,I).GE.-EPS).AND.
     .       (P(3,I)-Z1.GE.-EPS.AND.Z2-P(3,I).GE.-EPS)) THEN
            XYZG(1,IS+1)=P(1,I)
            XYZG(2,IS+1)=P(2,I)
            XYZG(3,IS+1)=P(3,I)
            IS=IS+1
          ENDIF
10      CONTINUE
      ELSEIF (RL.GT.-10.) THEN
        ILN=-RL
        DO 20 I=1,ILN
          IP=I+1
          DO 21 J=IP,ILN
C  SCHNITTPUNKT EBENE I, J, UND A0,A1,A2,A3
            DET=EIRENE_DETER(A1,XL(I,NUM),XL(J,NUM),A2,YL(I,NUM),
     .                YL(J,NUM),A3,ZL(I,NUM),ZL(J,NUM))
            IF (ABS(DET).LE.EPS) GOTO 21
            B1=-A0
            B2=-AL(I,NUM)
            B3=-AL(J,NUM)
            PS(1)=EIRENE_DETER(B1,B2,B3,A2,YL(I,NUM),YL(J,NUM),
     .                  A3,ZL(I,NUM),ZL(J,NUM))/DET
            PS(2)=EIRENE_DETER(A1,XL(I,NUM),XL(J,NUM),B1,B2,B3,
     .                  A3,ZL(I,NUM),ZL(J,NUM))/DET
            PS(3)=EIRENE_DETER(A1,XL(I,NUM),XL(J,NUM),A2,YL(I,NUM),
     .                  YL(J,NUM),B1,B2,B3)/DET
C  CHECKE, OB ALLE ANDEREN LINEAREN UNGLEICHUNGEN ERFUELLT SIND
            DO 40 K=1,ILN
              IF (K.EQ.I.OR.K.EQ.J) GOTO 40
              TEST=PS(1)*XL(K,NUM)+PS(2)*YL(K,NUM)+
     .             PS(3)*ZL(K,NUM)+AL(K,NUM)
              IF (TEST.GT.EPS) GOTO 21
40          CONTINUE
            IS=IS+1
            XYZG(1,IS)=PS(1)
            XYZG(2,IS)=PS(2)
            XYZG(3,IS)=PS(3)
21        CONTINUE
20      CONTINUE
      ENDIF
C
      IF (IS.LE.2) THEN
        WRITE (iunout,*) 'WENIGER ALS 3 SCHNITTPUNKTE ',
     .              'DER EBENE MIT DEM QUADER GEFUNDEN'
        WRITE (iunout,*) 'NUMMER DER FLAECHE: ',NUM
        RETURN
      ENDIF
C
      XMIT=0.
      YMIT=0.
      ZMIT=0.
      DO 200 I=1,IS
         XMIT=XMIT+XYZG(1,I)
         YMIT=YMIT+XYZG(2,I)
200      ZMIT=ZMIT+XYZG(3,I)
      XMIT=XMIT/DBLE(IS)
      YMIT=YMIT/DBLE(IS)
      ZMIT=ZMIT/DBLE(IS)
C
C     BESTIMME DIE WINKEL
      PI=4.*ATAN(1.)
      ICHECK=0
      ICOUNT=0
150   IF ((ABS(A2).LT.1.E-6.AND.ABS(A3).LT.1.E-6).OR.ICHECK.GT.1) THEN
        DO 300 I=1,IS
          ANGLE(I)=ATAN2((XYZG(3,I)-ZMIT),(XYZG(2,I)-YMIT))/PI*180.
300     CONTINUE
        ICHECK=1
      ELSEIF ((ABS(A1).LT.1.E-6.AND.ABS(A3).LT.1.E-6).OR.
     .        MOD(ICHECK,2).EQ.1) THEN
        DO 400 I=1,IS
          ANGLE(I)=ATAN2((XYZG(3,I)-ZMIT),(XYZG(1,I)-XMIT))/PI*180.
400     CONTINUE
        ICHECK=2
      ELSE
        DO 500 I=1,IS
          ANGLE(I)=ATAN2((XYZG(2,I)-YMIT),(XYZG(1,I)-XMIT))/PI*180.
500     CONTINUE
        ICHECK=3
      ENDIF
C
C     SORTIERE NACH WINKELN
      DO 600 I=1,IS-1
        ANG=ANGLE(I)
        ISORT=I
        DO 610 J=I+1,IS
          IF (ANGLE(J).LT.ANG) THEN
            ISORT=J
            ANG=ANGLE(J)
          ENDIF
610     CONTINUE
        DO 620 J=1,3
          HELP=XYZG(J,I)
          XYZG(J,I)=XYZG(J,ISORT)
620       XYZG(J,ISORT)=HELP
        HELP=ANGLE(I)
        ANGLE(I)=ANGLE(ISORT)
        ANGLE(ISORT)=HELP
600   CONTINUE
C
      ICOUNT=ICOUNT+1
      DO 700 I=1,IS-1
        DO 700 J=I+1,IS
           IF (ABS(ANGLE(J)-ANGLE(I)).LT.1.E-6.AND.
     .        (XYZG(1,I)-XYZG(1,J))**2+(XYZG(2,I)-XYZG(2,J))**2+
     .        (XYZG(3,I)-XYZG(3,J))**2.GT.1.E-6.AND.
     .         ICOUNT.LT.3) GOTO 150
700   CONTINUE
C
      II=0
      DO 1000 I=1,IS
        DO 1100 J=1,3
          II=II+1
1100      CORD(II)=XYZG(J,I)
1000  CONTINUE
      CALL EIRENE_PL3Q(CORD,IS,IO,NF)
C
      RETURN
      END
