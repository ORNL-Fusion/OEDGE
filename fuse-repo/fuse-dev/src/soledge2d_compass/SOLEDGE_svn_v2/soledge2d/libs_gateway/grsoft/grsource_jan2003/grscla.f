C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE   6. 11. 90 M. BUSCH    SAVE
C UPDATE  16. 12. 90 G. Groten   EPS
C UPDATE   8.  5. 91 G. GROTEN
C UPDATE   8.  1. 92 G. Groten   EPS neu wegen PC's
C UPDATE  10. 10. 95 G. Groten   ISA wegen eines Fehlers bei den Power WS
C------------------------------------------------------------------------
      SUBROUTINE GRSCLA(SA,ISA,SAA,EINT,FAKT)
      CHARACTER (len=8) PICC
      REAL   SA(3),SAA(4),RZAHL(15)
C---- Maschinenabhaengiges Epsilon fuer REAL
      REAL EPS,EPSI
      INTEGER   EINT,ISA(2)
      LOGICAL   LGROZA

      EQUIVALENCE (IGROZA,LGROZA)

      SAVE

      DATA RZAHL
     $/1.,2.,5.,10.,20.,50.,100.,200.,500.,1E3,2E3,5E3,1E4,2E4,5E4/

C---- EPS  : kleinste positive reelle Zahl X fuer die gilt '1.+X > 1.'
      EPS  = .000015258791
      EPSI = EPS*.5
   2  IF ( 1.+EPSI .GT. 1. ) THEN
         EPS = EPSI
         EPSI = EPS*.5
         GOTO 2
      ENDIF
      ISK=ISA(1)
      IF (SA(2).EQ.0. .AND. SA(3).EQ.0.) THEN
         SA(2)= EPS**5
         SA(3)=-EPS**5
      ENDIF

    1 ZW=ALOG10(MAX(ABS(SA(2)),ABS(SA(3))))
      ZW=ZW-1E-5*ABS(ZW)
      EINT=ZW
      IF (ZW.NE.EINT) EINT=EINT+1
      EINT=EINT-5
      ZW=10E0**(-EINT)
      IF (3.5*(SA(2)-SA(3))*ZW.LT.SA(1)) THEN
         FAKT=5.*ZW
         SA(2)=SA(2)+SA(1)/FAKT
         SA(3)=SA(3)-SA(1)/FAKT
         GOTO 1
      ENDIF

      ZW12=(SA(2)-SA(3))*ZW/SA(1)*1.4
      LGROZA=.FALSE.
      IF (ZW12.LT.12E0 .AND. MOD(ISK,100).EQ.1) THEN
         ISK=ISK+1
         LGROZA=.TRUE.
      ENDIF

      IF (SA(3).GT.-.01*SA(2)) THEN
         ZW13=.75*ZW12
      ELSE
         ZW13=ZW12
      ENDIF

      DO 5 I=1,14
         IF (RZAHL(I).GT.ZW13) GOTO 6
    5 CONTINUE

    6 SAA(4)=RZAHL(I)
      ISA(2)=IGROZA
      ISA(1)=ISK
      FAKT=1.4/ZW12
      SAA(2)=RZAHL(I)*FAKT
      ZW1=SA(3)*ZW
      ZW12=ZW1-1E-5
      DINT=ZW12/SAA(4)
      IDINT=DINT
      IF (DINT.NE.IDINT.AND.DINT.GE.0) IDINT=IDINT+1
      SAA(3)=IDINT*SAA(4)
      SAA(1)=(SAA(3)-ZW1)*FAKT
      FAKT=FAKT*ZW

    9 IF (SAA(4).GT.5. .AND. EINT.NE.0) THEN
         SAA(4)=IFIX(SAA(4))/10
         SAA(3)=IFIX(SAA(3))/10
         EINT=EINT+1
         GOTO 9
      ENDIF

      IF (SAA(3).LT.-100000.) THEN
         SAA(3)=SAA(3)+SAA(4)
         SAA(1)=SAA(1)+SAA(2)
      ENDIF

      GOTO 9999
C----------------------------------------------------------------------
      ENTRY GRLIN(SA,ISA,SAA,P1,P2,P3)
      ISK=ISA(1)
      IGROZA=ISA(2)
      HOBR=.26
      IF (MOD(ISK,100).NE.1 .AND. .NOT.LGROZA) HOBR=.4
      PIII=P3+HOBR-.26
      ZW=90.*P2
      CALL GRCHRC(HOBR,ZW,18)
      ZW1=SAA(1)-HOBR
      ZW12=SAA(3)
      Z=ZW1

   11 IF (Z.LE.SA(1)) THEN
         ZZ=(Z+.12)*P1-P2*PIII
         ZZZ=(Z+.12)*P2-P1*PIII
         KP=ZW12
         CALL GRPCTR(KP,LPIC,PICC)
         CALL GRTXT(ZZ,ZZZ,LPIC,PICC)
         ZW12=ZW12+SAA(4)*MOD(ISK,100)
         Z=Z+SAA(2)*MOD(ISK,100)
         GOTO 11
      ENDIF

 9999 END
