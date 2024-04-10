C
C
      SUBROUTINE EIRENE_PLTKA (FCN,XANF,XEND,INN,TRCPLT,LBOX)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CRECH
      USE EIRMOD_CLMSUR
      USE EIRMOD_CPLOT
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
 
      REAL(DP) :: FCN
      REAL(DP), INTENT(IN) :: XANF, XEND
      INTEGER, INTENT(OUT) :: INN
      LOGICAL, INTENT(IN) ::  TRCPLT, LBOX
 
      REAL(DP) :: XTRAN, YTRAN, XI, ETA, GERAX, XG, YG, XINC, DX,
     .          X, Y, XB1, YB1, XP, YP, XB2, YB2, XX, YY, XXO, YYO
      INTEGER :: I, IFLAG, INC
      LOGICAL :: LBX, LBXO, LPA, L1, L2, L3, L4
 
      XTRAN(XI,ETA)=XI*COSA-ETA*SINA
      YTRAN(XI,ETA)=XI*SINA+ETA*COSA
!pb      GERAY(XG)=(YY-YYO)/(XX-XXO)*XG+YY-XX*(YY-YYO)/(XX-XXO)
      GERAX(YG)=(YG-YY)*(XX-XXO)/(YY-YYO)+XX
      IF (TRCPLT) WRITE (iunout,*) 'PLTKA'
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
      INN=0
      X=XANF
      Y=FCN(X)
      XX=XTRAN(X+XM,Y+YM)
      YY=YTRAN(X+XM,Y+YM)
C*****LPA=.T. IF POINT INSIDE PLOTAREA
      LPA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
C*****LBX=.T. IF POINT INSIDE BOX
      LBX=XX.GE.XC1.AND.XX.LE.XC2.AND.YY.GE.YC1.AND.YY.LE.YC2
      IF (LPA.AND..NOT.LBX) THEN
        IF (LZR) CALL GRJMP
     .  (REAL(XX,KIND(1.E0)),REAL(YY,KIND(1.E0)))
        IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR (XX,YY,0)
        IFLAG=24
        INN=24
      ENDIF
      DO 1 I=1,INC
         X=XANF+I*DX
         Y=FCN(X)
         XXO=XX
         YYO=YY
         LBXO=LBX
         XX=XTRAN(X+XM,Y+YM)
         YY=YTRAN(X+XM,Y+YM)
         LPA=XX.GE.XL1.AND.XX.LE.XL2.AND.YY.GE.YL1.AND.YY.LE.YL2
         LBX=XX.GE.XC1.AND.XX.LE.XC2.AND.YY.GE.YC1.AND.YY.LE.YC2
         IF (IFLAG.EQ.0.AND.(LBX.OR..NOT.LPA)) GOTO 1
         IF (IFLAG.NE.0.AND.(LPA.AND..NOT.LBX)) THEN
            IF (LZR) CALL GRDRW (REAL(XX,KIND(1.E0)),
     .                           REAL(YY,KIND(1.E0)))
            IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR
     .  (XX,YY,1)
         ELSE
            IF (LBX.OR.LBXO) THEN
               XB1=XC1
               XB2=XC2
               YB1=YC1
               YB2=YC2
            ELSE
               XB1=XL1
               XB2=XL2
               YB1=YL1
               YB2=YL2
            ENDIF
            L1=XX.LE.XB1.OR.XXO.LE.XB1
            L2=XX.GE.XB2.OR.XXO.GE.XB2
            L3=YY.LE.YB1.OR.YYO.LE.YB1
            L4=YY.GE.YB2.OR.YYO.GE.YB2
            IF (L1) XP=XB1
            IF (L2) XP=XB2
            IF (L3) YP=YB1
            IF (L4) YP=YB2
CPB   IM FALLE EINER "SCHIEFEN" ELLIPSE FEHLER MOEGLICH
            IF (L1.OR.L2) YP=FCN(XP-XM)+YM
            IF (L3.OR.L4) XP=GERAX(YP)
            IF (IFLAG.EQ.0) THEN
               IF (LZR) THEN
                  CALL GRJMP
     .  (REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
                  CALL GRDRW
     .  (REAL(XX,KIND(1.E0)),REAL(YY,KIND(1.E0)))
               ENDIF
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) THEN
                  CALL EIRENE_STCOOR (XP,YP,0)
                  CALL EIRENE_STCOOR (XX,YY,1)
               ENDIF
               IFLAG=24
               INN=24
            ELSE
               IF (LZR) CALL GRDRW (REAL(XP,KIND(1.E0)),
     .                              REAL(YP,KIND(1.E0)))
               IF (.NOT.LZR.OR.PLSTOR.OR.PLNUMS) CALL EIRENE_STCOOR
     .  (XP,YP,1)
               IFLAG=0
            ENDIF
         ENDIF
1     CONTINUE
      RETURN
      END
