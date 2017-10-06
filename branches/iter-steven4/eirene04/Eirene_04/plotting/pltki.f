C
C
      SUBROUTINE PLTKI (FCN,XANF,XEND,INN,TRCPLT,LBOX)

      USE PRECISION
      USE PARMMOD
      USE CRECH
      USE CCONA
      USE CLMSUR
      USE CPLOT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XANF, XEND
      REAL(DP) :: FCN
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) :: TRCPLT,LBOX

      REAL(DP) :: XTRAN, YTRAN, XI, ETA, GERAX, GERAY, XG, YG, X, Y, DX,
     .          XP, YP, XX, YY, XXO, YYO, XINC
      INTEGER :: I, ISIDE, IFLAG, IJUMP, INC

      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (6,*) 'PLTKI'
      IF (LBOX) THEN
        CALL GRJMP(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
        CALL GRDRW(REAL(XL2,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL2,KIND(1.E0)))
        CALL GRDRW(REAL(XL1,KIND(1.E0)),REAL(YL1,KIND(1.E0)))
      ENDIF
      INC=100
      XINC=DBLE(INC)
      DX=(XEND-XANF)/XINC
      IFLAG=0
      ISIDE=0
      INN=0
      XX=0.
      YY=0.
      DO 431 I=1,INC+1
         X=XANF+(I-1)*DX
         Y=FCN(X)
         XXO=XX
         YYO=YY
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         IJUMP=ISIDE
         ISIDE=0
         IF (XX.LE.XL1+EPS10) ISIDE=3
         IF (XX.GE.XL2-EPS10) ISIDE=2
         IF (YY.LE.YL1+EPS10) ISIDE=4
         IF (YY.GE.YL2-EPS10) ISIDE=1
         IF (ISIDE.NE.0) THEN
             IF (IFLAG.EQ.0) GOTO 431
             IJUMP=ISIDE
             IF (IFLAG.NE.0) GOTO 432
         ENDIF
C        ISIDE=0
         IF (IFLAG.EQ.0.AND.I.NE.1) GOTO 432
439      IF (I.NE.1) IFLAG=24
         IF (IFLAG.EQ.0) THEN
            IF (LZR) CALL GRJMP (REAL(XX,KIND(1.E0)),
     .                           REAL(YY,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
         IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                        REAL(YY,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XX,YY,1)
         ENDIF
         INN=MAX0(INN,IFLAG)
         GOTO 431
432      GOTO (451,452,453,454),IJUMP
         GOTO 439
451      XP=GERAX(YL2)
         YP=YL2
         GOTO 429
452      XP=XL2
CPB      YP=YTRAN(XP,FCN(XP-XM)+YM)
         YP=GERAY(XL2)
         GOTO 429
453      XP=XL1
CPB      YP=YTRAN(XP,FCN(XP-XM)+YM)
         YP=GERAY(XL1)
         GOTO 429
454      XP=GERAX(YL1)
         YP=YL1
         GOTO 429
C
429      CONTINUE
         IF (IFLAG.EQ.0) THEN
            IF (LZR) CALL GRJMP (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,0)
         ENDIF
         IF (IFLAG.NE.0) THEN
            IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                           REAL(YP,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL STCOOR (XP,YP,1)
         ENDIF
         INN=MAX0(INN,IFLAG)
         IF (IFLAG.EQ.0) GOTO 439
         IFLAG=0
C
431   CONTINUE
      RETURN
      END
