C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C  UPDATE:  25. 9.1990 GROTEN
      SUBROUTINE GRDRHS(PA,NP,PT,XX,YY)
      DIMENSION PA(15),PT(3,NP+1),XX(*),YY(*)

      PARAMETER(NROW=128,NROWD2=NROW/2,NCOLD2=NROWD2)
      DIMENSION XD(NROW),YD(NROW),TAB(NROW,NROW),IY(64)
      INTEGER HBX,HBY,HBXM1,HBYM1,HBXA,HBYA

      NP1=NP+1
      LBX=PA(12)
      LBY=PA(14)
      HBX=PA(13)
      HBY=PA(15)
      XMIN=PA(6)
      XMAX=PA(7)
      YMIN=PA(8)
      YMAX=PA(9)
      ZMIN=PA(10)
      IYDIF=HBY-LBY
      IXDIF=HBX-LBX
      IF (LBX.GE.1) GOTO 10
         PA(5)=10.
         GOTO 9999
   10 IF (HBX-LBX+1 .LE. NROWD2) GOTO 20
         PA(5)=11.
         GOTO 9999
   20 IF (HBX.GT.LBX) GOTO 30
         PA(5)=12.
         GOTO 9999
   30 HBXM1=HBX-1
      LBXP1=LBX+1
      DO 31 I=LBX,HBXM1
         IF (XX(I).LE.XX(I+1)) GOTO 31
            PA(5)=13.
            GOTO 9999
   31 CONTINUE
      IF (XMIN.LE.XX(LBX)) GOTO 35
         PA(5)=14.
         GOTO 9999
   35 IF (XMAX.GE.XX(HBX)) GOTO 40
         PA(5)=15.
         GOTO 9999
   40 IF (XMIN.NE.XMAX) GOTO 50
         PA(5)=16.
         GOTO 9999
   50 IF (LBY.GE.1) GOTO 60
         PA(5)=20.
         GOTO 9999
   60 IF (IYDIF+1 .LE. NCOLD2) GOTO 70
         PA(5)=21.
         GOTO 9999
   70 IF (HBY.GT.LBY) GOTO 80
         PA(5)=22.
         GOTO 9999
   80 HBYM1=HBY-1
      DO 81 I=LBY,HBYM1
         IF (YY(I).LE.YY(I+1)) GOTO 81
            PA(5)=23.
            GOTO 9999
   81 CONTINUE
      IF (YMIN.LE.YY(LBY)) GOTO 85
         PA(5)=24.
         GOTO 9999
   85 IF (YMAX.GE.YY(HBY)) GOTO 90
         PA(5)=25.
         GOTO 9999
   90 IF (YMIN.NE.YMAX) GOTO 100
         PA(5)=26.
         GOTO 9999
C
C  SORTIERE DIE PUNKTE NACH Y AUFSTEIGEND
C
  100 X=NP+2
      II=ALOG(X)/ALOG(2.)
      K=2**II-1
      II=II-1
      DO 181 LL=1,II
         K=K/2
         DO 180 J=1,K
            IA=J+K
            DO 170 I=IA,NP,K
               I1=I+1
               L=I-K+1
               IF ( PT(2,L).LE.PT(2,I1) ) GOTO 170
               PTX=PT(1,I1)
               PTY=PT(2,I1)
               PTZ=PT(3,I1)
               PT(1,I1)=PT(1,L)
               PT(2,I1)=PT(2,L)
               PT(3,I1)=PT(3,L)
               L=L-K
               LI=L-1
               DO 160 IN=J,LI,K
                  IF ( PT(2,L).LE.PTY ) GOTO 165
                  IM=L+K
                  PT(1,IM)=PT(1,L)
                  PT(2,IM)=PT(2,L)
                  PT(3,IM)=PT(3,L)
                  L=L-K
  160          CONTINUE
  165          IM=L+K
               PT(1,IM)=PTX
               PT(2,IM)=PTY
               PT(3,IM)=PTZ
  170       CONTINUE
  180 CONTINUE
 181  CONTINUE
C
C   GRENZEN DER Y-INTERVALLE
C
  200 IR=0
      DO 201 IK=2,NP1
         IF (PT(2,IK).GE.YY(LBY)) GOTO 202
         IR=IR-1
  201 CONTINUE
  202 IY(1)=IK
      DO 205 I=LBY+1,HBY
         IKM=IK
         J=I-LBY+1
         DO 203 IK=IKM,NP1
            IF (PT(2,IK).GT.YY(I)) GOTO 204
  203    CONTINUE
         IK=NP1+1
  204    IY(J)=IK
  205 CONTINUE
      IR=IR-(NP1-IK+1)
C
C  SORTIERE DIE PUNKTE INNERHALB DER Y-INTERVALLE NACH X AUFSTEIGEND
C
      DO 399 KY=1,IYDIF
         IYU=IY(KY)-1
         NPY=IY(KY+1)-IY(KY)
         X=NPY+2
         II=ALOG(X)/ALOG(2.)
         K=2**II-1
         II=II-1
         DO 381 LL=1,II
            K=K/2
            DO 380 J=1,K
               IA=J+K
               DO 370 I=IA,NPY,K
                  I1=I+IYU
                  L=I-K+IYU
                  IF ( PT(1,L).LE.PT(1,I1) ) GOTO 370
                  PTX=PT(1,I1)
                  PTY=PT(2,I1)
                  PTZ=PT(3,I1)
                  PT(1,I1)=PT(1,L)
                  PT(2,I1)=PT(2,L)
                  PT(3,I1)=PT(3,L)
                  L=L-K
                  LI=L-IYU
                  IF (LI.LT.J) GOTO 365
                  DO 360 IN=J,LI,K
                     IF ( PT(1,L).LE.PTX ) GOTO 365
                     IM=L+K
                     PT(1,IM)=PT(1,L)
                     PT(2,IM)=PT(2,L)
                     PT(3,IM)=PT(3,L)
                     L=L-K
  360             CONTINUE
  365             IM=L+K
                  PT(1,IM)=PTX
                  PT(2,IM)=PTY
                  PT(3,IM)=PTZ
  370          CONTINUE

  380    CONTINUE
 381    CONTINUE
  399 CONTINUE
C
C  ERZEUGE DIE DOPPELTEN X/Y- KOORDINATEN
C
      LBXA=LBX
      HBXA=HBX
      LBYA=LBY
      HBYA=HBY
      LBX=1
      PA(12)=LBX
      HBX=2*(HBXA-LBXA)
      PA(13)=HBX
      LBY=1
      PA(14)=LBY
      HBY=2*(HBYA-LBYA)
      PA(15)=HBY
      XD(LBX)=XX(LBXA)
      LX1=LBX+1
      LX2=HBX-2
      J=LBXA+1
      DO 210 I=LX1,LX2,2
         XD(I)=XX(J)
         XD(I+1)=XX(J)
         J=J+1
  210 CONTINUE
      XD(HBX)=XX(HBXA)
      YD(LBY)=YY(LBYA)
      LY1=LBY+1
      LY2=HBY-2
      J=LBYA+1
      DO 220 I=LY1,LY2,2
         YD(I)=YY(J)
         YD(I+1)=YY(J)
         J=J+1
 220  CONTINUE
      YD(HBY)=YY(HBYA)
C
C      EINORDNEN DER PUNKTE IN DIE MATRIX
C
      XMI=XX(LBXA)
      XMA=XX(HBXA)
      DO 490 I=1,IYDIF
      IKMI=IY(I)
      IKMA=IY(I+1)-1
      DO 410 IK=IKMI,IKMA
         IF (PT(1,IK).GE.XMI) GOTO 411
         IR=IR-1
  410 CONTINUE
  411 IKMI=IK
      IK=IKMA
      IF (IKMI.GT.IKMA) GOTO 421
      DO 420 J=IKMI,IKMA
         IF ( PT(1,IK).LE.XMA) GOTO 421
         IR=IR-1
         IK=IK-1
  420 CONTINUE
  421 IKMA=IK
      I1=2*I-1
      I2=I1+1
      IK=IKMI
      DO 480 J=1,IXDIF
      J1=2*J-1
      J2=J1+1
      IF ( PT(1,IK).LE.XD(J2) .AND. IK.LE.IKMA ) GOTO 440
      TAB(J1,I1)=ZMIN
      TAB(J1,I2)=ZMIN
      TAB(J2,I1)=ZMIN
      TAB(J2,I2)=ZMIN
      GOTO 480
  440 SU=0.
      IP=0
      IKMI=IK
      DO 450 IK=IKMI,IKMA
         IF ( PT(1,IK).GT.XD(J2) ) GOTO 451
         SU=SU+PT(3,IK)
         IP=IP+1
  450 CONTINUE
  451 SU=SU/IP
      TAB(J1,I1)=SU
      TAB(J1,I2)=SU
      TAB(J2,I1)=SU
      TAB(J2,I2)=SU
  480 CONTINUE
  490 CONTINUE
      CALL GRDRDM(PA,NROW,TAB,XD,YD)
      PA(12)=LBXA
      PA(13)=HBXA
      PA(14)=LBYA
      PA(15)=HBYA
      PA(5)=IR
      IF (IR.EQ.0) GOTO 777
 9999 IR=PA(5)
  777 END
