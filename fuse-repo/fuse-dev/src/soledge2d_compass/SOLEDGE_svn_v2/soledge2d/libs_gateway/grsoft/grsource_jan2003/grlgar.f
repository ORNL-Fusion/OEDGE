C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE   6. 11. 90 M. BUSCH    SAVE
C UPDATE  10. 10. 95 M. BUSCH    ISA wegen Fehler bei Power WS
C-----------------------------------------------------------------------
      SUBROUTINE GRLGAR(SA,ISA,SAA,ALO,FAKT,IDECA)
      REAL        ALO(8),SA(3),SAA(4)
      INTEGER     ISA(2)
      CHARACTER (len=8) PICC

      SAVE

      X=ALOG10(SA(3))+1E-4
      I1=X
      IF(X.LT.0) I1=I1-1
      X=ALOG10(SA(2))-1E-4
      I2=X
      IF(X.GT.0.) I2=I2+1
      IF(I1.EQ.I2) I1=I1-1
      IDECA=0
   11 IDECA=IDECA+1
      II=I2-I1
      II=MOD(II,IDECA)
      IF(II.GT.0) I2=I2+IDECA-II
      FAKT=SA(1)/(I2-I1)
      SAA(2)=IDECA*FAKT
      IF( SAA(2).LT.2. ) GOTO 11
      SA(3)=10E0**I1*.9995
      SA(2)=10E0**I2*1.0005
      SAA(3)=I1
      SAA(4)=I2
      ALO(1)=SAA(2)*.3010
      ALO(2)=SAA(2)*.6990
      IF(SAA(2).LT.12.) GOTO 99
      ALO(4)=ALO(2)
      ALO(2)=SAA(2)*.4771
      ALO(3)=SAA(2)*.6021
      ALO(5)=SAA(2)*.7782
      ALO(6)=SAA(2)*.8451
      ALO(7)=SAA(2)*.9031
      ALO(8)=SAA(2)*.9542
      GOTO 99
C
      ENTRY GRLOG(SA,ISA,SAA,ALO,P1,P2,P3,IDECA)
      ISK=ISA(1)
      HOBR=.3
      IF (MOD(ISK,100).NE.-1) HOBR=.5
      HOBRK=HOBR*.666667
      PIII=P3+HOBR-.3
      PIV=2.*HOBR
      PV=HOBRK*.5
      ZW=90.*P2
      IZZ=SAA(3)
      SA(1)=SA(1)-.01
      Z=-HOBR
   13 ZZ=Z*P1-P2*PIII
      ZZZ=Z*P2-P1*P3
      CALL GRCHRC(HOBR,ZW,18)
      CALL GRTXT(ZZ,ZZZ,2,'10')
      ZZ=ZZ+PIV*P1-PV*P2
      ZZZ=ZZZ+PV*P1+PIV*P2
      CALL GRPCTR(IZZ,LPIC,PICC)
      CALL GRCHRC(HOBRK,ZW,18)
      CALL GRTXT(ZZ,ZZZ,LPIC,PICC)
      IF (IZZ.EQ.SAA(4)) GOTO 99
      IF (SAA(2).LE.4.) GOTO 15
      IP=2
      ZZ=ZZ-(HOBR+PV)*P1
      ZZZ=ZZZ-(HOBR+PV)*P2
      K=8
      IF(SAA(2).LT.12.) K=2
      DO 14 I=1,K
         ZT=ZZ+ALO(I)*P1
         ZZT=ZZZ+ALO(I)*P2
         II=MOD(IP,10)
         CALL GRPCTR(II,LPIC,PICC)
         CALL GRTXT(ZT,ZZT,LPIC,PICC)
         IP=IP+1
         IF(K.EQ.2)IP=5
   14    CONTINUE
   15 IZZ=IZZ+IDECA
      Z=Z+SAA(2)
      GOTO 13
   99 END
