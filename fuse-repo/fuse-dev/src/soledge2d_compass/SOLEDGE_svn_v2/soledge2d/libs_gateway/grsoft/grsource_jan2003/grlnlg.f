C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 16. 8.1990 GROTEN
      SUBROUTINE GRLNLG(ANF,ACM,BCM,AP,CM1,CM2,ALO,P1,P2,P3,P4,AXZA)
C---- Groten 16.8.90 fuer neues GKS und CRAY
      LOGICAL LIN
      REAL   ALO(8)
      INTEGER   AXZA(32)

      LIN=AXZA(1).GE.0
      ZW12=ANF+CM1+AP
      ZW=ACM*P2
      ZW1=BCM*P1
      CCM=ABS(CM2)
      GRENZ=ACM*P1+BCM*P2+.1
      Z=ZW12
      L=1
    1 IF ((Z.LT.GRENZ .OR. ANF.NE.0.) .AND. Z.GT.-.1) GOTO 2
      AP=GRENZ+CM2-Z-.1
      GOTO 6
    2 ZZ=Z*P1+ZW
      ZZZ=Z*P2+ZW1
      IF (.NOT.LIN) GOTO 21
      L=L+1
      IP=AXZA(1)-L+2
      IF (CM2.LE.0) GOTO 20
      AXZA(1)=L
      IP=L
   20 PIII=P3
      IF (AXZA(IP).LE.0) PIII=P4
      GOTO 22
   21 PIII=P3
   22 ZT=ZZ+PIII*P2
      ZZT=ZZZ+PIII*P1
      CALL GRDRW(ZZ,ZZZ)
      CALL GRDRW(ZT,ZZT)
      CALL GRJMP(ZZ,ZZZ)
      IF (LIN) GOTO 5
      IF (Z.GE.GRENZ-.2 .AND. ANF.EQ.0. .OR. Z.LE..1 .AND. ANF.NE.0.)
     1   GOTO 7
      K=8
      IF (CCM.LT.12) K=2
      IF (CM2.GE.0.) GOTO 3
      ZZ=ZZ+CM2*P1
      ZZZ=ZZZ+CM2*P2
    3 DO 4 I=1,K
      IP=I
      IF (CM2.LT.0.) IP=K-I+1
      ZC=ZZ+ALO(IP)*P1
      ZZD=ZZZ+ALO(IP)*P2
      ZT=ZC+P4*P2
      ZZT=ZZD+P4*P1
      CALL GRDRW(ZC,ZZD)
      CALL GRDRW(ZT,ZZT)
      CALL GRJMP(ZC,ZZD)
    4 CONTINUE
    5 Z=Z+CM2
      GOTO 1
    6 CALL GRDRW(ACM,BCM)
    7 END
