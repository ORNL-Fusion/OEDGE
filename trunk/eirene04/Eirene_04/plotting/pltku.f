C
C
      SUBROUTINE PLTKU (FCN,XANF,XEND,INN,TRCPLT,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CLMSUR
      USE CPLOT

      IMPLICIT NONE

      REAL(DP) :: FCN
      REAL(DP), INTENT(IN) :: XANF, XEND
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) :: TRCPLT, LBOX

      REAL(DP) :: XTRAN, YTRAN, XI, ETA, DX, XINC, XOUT, YOUT, XS, YS,
     .          XOK, YOK, XXO, YYO, XP, YP, XTT, YTT, XT, YT, GL, GS,
     .          X, Y, XX, YY, DXX
      INTEGER :: INC, IFLAG, I, J
      LOGICAL :: LIN, LINO, LTST, LPLA, LPLAO, LPL, LGL, L1, L2, L3, L4

      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      IF (TRCPLT) WRITE (6,*) 'PLTKU'
      IF (LBOX) THEN
        CALL GRJMP(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
      ENDIF
C
      INC=100
100   CONTINUE
      XINC=DBLE(INC)
      DX=(XEND-XANF)/XINC
      INN=0
      IFLAG=0
C
      DO 1 I=1,MLIN
1        ALIN(I)=ALIN(I)+ZLIN(I)*ZPLT
      DO 2 I=1,MSCN
         A0S(I)=A0S(I)+(A3S(I)+A6S(I)*ZPLT)*ZPLT
         A1S(I)=A1S(I)+A8S(I)*ZPLT
         A2S(I)=A2S(I)+A9S(I)*ZPLT
2     CONTINUE
      DO 10 I=1,INC+1
         X=XANF+(I-1)*DX
         Y=FCN(X)
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         LPLA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
         LIN=.TRUE.
         DO 11 J=1,MLIN
            GL=ALIN(J)+XLIN(J)*XX+YLIN(J)*YY
            LIN=LIN.AND.GL.LE.0.
11       CONTINUE
         DO 12 J=1,MSCN
            GS=A0S(J)+(A1S(J)+A4S(J)*XX)*XX+(A2S(J)+
     .         A5S(J)*YY)*YY+A7S(J)*XX*YY
            LIN=LIN.AND.GS.LE.0.
12       CONTINUE
C
         IF (I.EQ.1) THEN
C*****FIRST POINT OF LINE
            IF (LIN.AND.LPLA) THEN
               IF (LZR) CALL GRJMP (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
               INN=24
               IFLAG=24
            ENDIF
         ELSE
C*****ALL OTHER POINTS
C*****POINT OUT OF PLOT AREA
            IF (.NOT.LPLA.AND..NOT.LPLAO) GOTO 16
C*****POINT OUT OF CONFIGURATION
            IF (.NOT.LIN.AND..NOT.LINO) GOTO 16
C*****POINT INSIDE PLOTAREA AND INSIDE THE CONFIGURATION
            IF (LIN.AND.LINO.AND.LPLA.AND.LPLAO) THEN
               IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,1)
            ELSE
C*****LINE CROSSES BOUNDARY
C***** ... OF CONFIGURATION
               LGL=LIN.AND.LINO
C***** ... OF PLOT AREA
               LPL=LPLA.AND.LPLAO
               IF (.NOT.LPL) THEN
                  L1=XX.LE.XL1.OR.XXO.LE.XL1
                  L2=XX.GE.XL2.OR.XXO.GE.XL2
                  L3=YY.LE.YL1.OR.YYO.LE.YL1
                  L4=YY.GE.YL2.OR.YYO.GE.YL2
                  IF (LPLA) THEN
                     XOK=XX
                     YOK=YY
                     XOUT=XXO
                     YOUT=YYO
                  ELSE
                     XOK=XXO
                     YOK=YYO
                     XOUT=XX
                     YOUT=YY
                  ENDIF
               ENDIF
C*****INTERPOLATE POINT ON CONFIGURATION BOUNDARY
               IF (.NOT.LGL) THEN
                  IF (LIN) THEN
                     XS=XANF+(I-1)*DX
                     DXX=-DX
                  ELSE
                    XS=XANF+(I-2)*DX
                    DXX=DX
                  ENDIF
C
13                DXX=DXX*0.5
                  XT=XS+DXX
                  YT=FCN(XT)
                  XTT=XTRAN(XT+XM,YT+YM)
                  YTT=YTRAN(XT+XM,YT+YM)
                  LTST=.TRUE.
                  DO 14 J=1,MLIN
                     GL=ALIN(J)+XLIN(J)*XTT+YLIN(J)*YTT
                     LTST=LTST.AND.GL.LE.0.
14                CONTINUE
                  DO 15 J=1,MSCN
                     GS=A0S(J)+(A1S(J)+A4S(J)*XTT)*XTT+
     .                  (A2S(J)+A5S(J)*YTT)*YTT+A7S(J)*XTT*YTT
                     LTST=LTST.AND.GS.LE.0.
15                CONTINUE
                  IF (LTST) XS=XS+DXX
                  IF (ABS(DXX).GT.1.E-3) GOTO 13
C
                  YS=FCN(XS)
                  XP=XTRAN(XS+XM,YS+YM)
                  YP=YTRAN(XS+XM,YS+YM)
                  IF (XP.GE.XL1.AND.XP.LE.XL2.AND.
     .                YP.GE.YL1.AND.YP.LE.YL2) THEN
C*****INTERPOLATED POINT INSIDE PLOT AREA, CAN BE PLOTTED
                     IF (LIN) THEN
                        IF (LZR) THEN
                           CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                                 REAL(YP,KIND(1.E0)))
                           CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                                 REAL(YY,KIND(1.E0)))
                        ENDIF
                        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                           CALL STCOOR (XP,YP,0)
                           CALL STCOOR (XX,YY,1)
                        ENDIF
                        INN=24
                        IFLAG=24
                     ELSE
                        IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                                       REAL(YP,KIND(1.E0)))
                        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                           CALL STCOOR (XP,YP,1)
                        ENDIF
                        IFLAG=0
                     ENDIF
                  ELSE
C*****INTERPOLATED POINT OUTSIDE PLOTAREA
                     L1=XP.LE.XL1
                     L2=XP.GE.XL2
                     L3=YP.LE.YL1
                     L4=YP.GE.YL2
                     XOUT=XP
                     YOUT=YP
                     IF (LPL) THEN
                        IF (LIN) THEN
                           XOK=XX
                           YOK=YY
                        ELSE
                           XOK=XXO
                           YOK=YYO
                        ENDIF
                     ENDIF
                     LPL=.FALSE.
                  ENDIF
               ENDIF
               IF (.NOT.LPL) THEN
C*****INTERPOLATE POINT ON BOUNDARY OF PLOTAREA
                  IF (L1) XP=XL1
                  IF (L2) XP=XL2
                  IF (L3) YP=YL1
                  IF (L4) YP=YL2
                  IF (L1.OR.L2) YP=YOK-(XOK-XP)*(YOK-YOUT)/(XOK-XOUT)
                  IF (L3.OR.L4) XP=XOK-(YOK-YP)*(XOK-XOUT)/(YOK-YOUT)
                  IF (IFLAG.EQ.0) THEN
                     IF (LZR) THEN
                        CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                              REAL(YP,KIND(1.E0)))
                        CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                              REAL(YY,KIND(1.E0)))
                     ENDIF
                     IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                        CALL STCOOR (XP,YP,0)
                        CALL STCOOR (XX,YY,1)
                     ENDIF
                     INN=24
                     IFLAG=24
                  ELSE
                     IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                                    REAL(YP,KIND(1.E0)))
                     IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                       CALL STCOOR (XP,YP,1)
                     ENDIF
                     IFLAG=0
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
16       XXO=XX
         YYO=YY
         LINO=LIN
         LPLAO=LPLA
10    CONTINUE
      IF (INN.GT.0.OR.INC.GE.400) RETURN
      INC=INC*2
      GOTO 100
      END
