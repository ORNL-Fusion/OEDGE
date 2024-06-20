C@PROCESS NOGOSTMT OPT(3) IL(DIM) NOSDUMP FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 16.8.90 G.Groten:   Strichdicke durch 2000*i neues GKS
C UPDATE 09.6.89 G.GROTEN:   Strichdicke durch +2000
C UPDATE 19.9.86 G.GROTEN:   WENN XXCM UND YYCM BEIDE NEGATIV, KEINE
C                            VERSCHIEBUNG UM 3. BZW. 3.5 CM. KEINE
C                            STAUCHUNG AM RECHTEN RAND ODER OBEN.
C UPDATE 7.5.85  H.GERLACH   NBIL<0      KEIN VORSCHUB!
C UPDATE 27.9.90 M. BUSCH    Daten werden von den Routinen KURVEF,
C                            SKURF, PSKURF formatfrei geschrieben
C                            hier formatfrei lesen
C                            frueher (20A4) d.h. CRAY unvertraeglich
C UPDATE  6.11.90 M. BUSCH   SAVE
C UPDATE 17.12.90 M. BUSCH   ANPASSUNG AN GTS GRAL MARKER
C UPDATE 10. 1.91 M. BUSCH   FESTE CM ANGABEN BISHER NUR FUER DIN A3
C                            COMMON GRPIC EINGEFUEGT - BELIEBIG BILD-
C                            GROESSE MOEGLICH
C UPDATE 10. 4.91 M. BUSCH   I/O AUF UNIT 91,92 STAAT 4,36
C UPDATE 21. 5.91 M. BUSCH   Startpunkt mit Marker
C UPDATE 18. 7.91 M. BUSCH   Logik :
C                            Bild wird erst beim 2.Aufruf von
C                            GRBLD bzw. explizit durch ein GRNXTF gemalt
C                            alle Bilder bis auf das letzte werden von
C                            GRBLD ausgegeben - hat man nur 1 Bild so
c                            wird dies wie auch das letzte Bild
c                            durch GREND beendet.
C                            Vorgehen ist notwendig, damit nach dem CALL
C                            von GRBLD anschliessend Legende, Text usw.
C                            ausgegeben werden kann
C                            damit alles bzgl. Skalierung richtig laeuft
c                            muss der GRNXTF Aufruf als ERSTES erfolgen
C UPDATE 25. 9.91 G. Groten  Startpunkt mit Marker korrigiert
C UPDATE 22.10.91 G. Groten  logarithmisch ohne Achsen: IISK,JJSK<=-999
      SUBROUTINE GRBLD(XXCM,YYCM,IISK,JJSK,XMI,XMA,YMI,YMA,NBIL)
      IMPLICIT NONE
      INTEGER IP,I1,L,K,J,LL,IST1,IST,ISY,JSY,IFA
      REAL XL,YL,XQ,YQ,PPP3,QQQ3,CMX2,CMY2,qq3,cmx1,cmy1
      REAL xcm, xmax, rxh, ryh, xmin, ymax, ymin, ycm, xmxcm1, ymxcm1
      REAL ya, xa, xdcpic, ydcpic, pp, xmi, xxcm, yycm, cilber, yma, xma
      REAL ymi,hobr,ymaxf,fakty,pp3,yp,xp,ycms,a1,xcms,yus,xls,xmaxf
      REAL faktx,a4,a2,a3
      INTEGER kp, intx1, jsk, intx2, inty1, inty2, isk, nsclp, nsclc
      INTEGER nsclv, nkurv, jjsk, iisk, nbil, nbild, i, idecay, mo, lpic
      INTEGER idecax
      COMMON /GRCIL/ CILBER(8)
CDEC$ PSECT /GRCIL/ NOSHR
      INTEGER   AXZAX(32),AXZAY(32),IO
      LOGICAL   IGROZA,JGROZA,BLDFRA
C---- Grisy used in GRLGND
      INTEGER PISY
      COMMON /GRISY/ PISY(14)
CDEC$ PSECT /GRISY/ NOSHR
      LOGICAL   ANFANG,VORAUS,EXCH,SPREIZ,FLBIT
      LOGICAL   MARKER
      INTEGER   ISYS   ,POTX,POTY,ARTLAB,RETLAB,BLERI,BTOBO,BSEPO,
     * BFIPO,BAUS,BHEPO,AAUS
      PARAMETER(IO=92,ISYS=91)
      INTEGER  SYMVEC(13), RAHMEN, FLPIC, KSX(2), KSY(2)
      REAL   SX(3),SXX(4),SY(3),SYY(4),A(100),B(100)
      REAL   ALO1(8),ALO2(8)
      REAL   SH(3,5)
      CHARACTER(len=8) PICC
CCC   COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
CCC
CCC
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMXCM1,YMXCM1, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      INTEGER INTLIN,INTSYM,SIZMRK

      EQUIVALENCE (PP(13),INTLIN),(PP(16),INTSYM),(PP(17),SIZMRK)
CCC
CCC
      EQUIVALENCE (SX(1),XCM), (SX(2),XMAX), (SX(3),XMIN), (KSX(1),ISK)
      EQUIVALENCE (SY(1),YCM), (SY(2),YMAX), (SY(3),YMIN), (KSY(1),JSK)
      EQUIVALENCE (KSX(2),IGROZA),(KSY(2),JGROZA)
      EQUIVALENCE (SXX(1),CMX1), (SXX(2),CMX2), (SXX(3),INTX1),
     1            (SXX(4),INTX2)
      EQUIVALENCE (SYY(1),CMY1), (SYY(2),CMY2), (SYY(3),INTY1),
     1            (SYY(4),INTY2)

      SAVE /GRPP/, /GRPIC/,/GRCIL/,/GRISY/

      DATA   SH     / .25   , .25   , .25   ,
     1                .25   , .0625 , .03125,
     2                .03125, .09375, .03125,
     3                .5    , .25   , .25   ,
     4                1.     , .0625 ,1.      /
C     ANPSSUNG AN GTS GRAL MARKER
C     NORMIERT 1-5, GTSGRAL MARKER -101,..-110
C     MARKER 1 NICHT, DA NUR PUNKT

      DATA SYMVEC/2,3,4,5,-101,-102,-103,-108,-109,-110,-111,-112,-113/


      FLBIT=.FALSE.
      BLDFRA=XXCM.GE.0. .OR. YYCM.GE.0.
      IF (BLDFRA) THEN
         XA=PP(1)
         YA=PP(2)
      ENDIF
      GOTO 444
C
C GRFL      FOLIE: DRAWS KURVES INTO AN OLD WINDOW FROM PREVIOUS GRBLD
C
      ENTRY GRFL
      FLBIT=.TRUE.

  444 REWIND ISYS
      NKURV=CILBER(3)
  777 IF (CILBER(8).GT.0.AND.NBIL.GT.0)  CALL GRNXTF
      NBILD=NKURV
      CILBER(3)=0.
      IF (FLBIT) GOTO 26
C
C  AUFSTELLEN DER STANDARDPARAMETERLISTE
C  UEBERNAHME DER ARGUMENTE, FALLS SIE VON 666 VERSCHIEDEN
C
C BUSCH (VON HIER
C YUS Y UNTERE ECKE
C XLS XLINKE ECKE
      YUS = 3.5 * YMXCM1/28.7
      XLS = 3. * XMXCM1/39.5
      XCMS = 33. * XMXCM1 / 39.5
      YCMS = 24. * YMXCM1 / 28.7
      XCM=ABS(XXCM)
      IF (XCM.EQ.666.) XCM= XCMS
      YCM=ABS(YYCM)
      IF (YCM.EQ.666.) YCM= YCMS
      ISK=IISK
      IF (IISK.EQ.666) ISK=1
      JSK=JJSK
      IF (JJSK.EQ.666) JSK=1
      XMIN=CILBER(4)
      IF (XMI.NE.666.) XMIN=XMI
      XMAX=CILBER(5)
      IF (XMA.NE.666.) XMAX=XMA
      YMIN=CILBER(6)
      IF (YMI.NE.666.) YMIN=YMI
      YMAX=CILBER(7)
      IF (YMA.NE.666.) YMAX=YMA
      IF (NBIL.NE.666.AND.NBIL.NE.-666) NBILD=ABS(NBIL)
      IF(ISK.GT.-999) ISK=MOD(ISK,300)
      IF (IABS(ISK).EQ.200 .OR. IABS(ISK).EQ.100) ISK=0
      IF (JSK.GT.-999) JSK=MOD(JSK,300)
      IF (IABS(JSK).EQ.200 .OR. IABS(JSK).EQ.100) JSK=0
      IF (    NKURV  .NE. 0
     1    .AND. XCM  .GE. 1.
     2    .AND. YCM  .GE. 1.
     3    .AND. XMAX .GE. XMIN
     4    .AND. YMAX .GE. YMIN
     5    .AND. (ISK.GE.0 .OR. XMIN.GT.0)
     6    .AND. (JSK.GE.0 .OR. YMIN.GT.0) )
     7    GOTO 3
      WRITE(IO,1)
      WRITE(IO,2)XCM,YCM,ISK,JSK,XMIN,XMAX,YMIN,YMAX,NBILD,NKURV
      WRITE(6,1)
      WRITE(6,2)XCM,YCM,ISK,JSK,XMIN,XMAX,YMIN,YMAX,NBILD,NKURV
    1 FORMAT(5X,'GRBLD: Fehler in Bildparameter')
    2 FORMAT(1X,2F7.2,2I5,4E10.2,2I4)
      STOP
C
C  BERECHNUNG DER ACHSENEINTEILUNGEN
C
3     IF (CILBER(8).EQ.0 .AND. BLDFRA) THEN
C BUSCH
         A1=XLS+PP(1)
         A2=YUS+PP(2)
         A3=PP(1)+XLS+XCM
         A4=PP(2)+YUS+YCM
         CALL GRSCLC(A1,A2,A3,A4)
      ELSE
         A1=PP(1)
         A2=PP(2)
         A3=PP(1)+XCM
         A4=PP(2)+YCM
         CALL GRSCLC(A1,A2,A3,A4)
      ENDIF
      CALL GRSCLV(0.,0.,XCM,YCM)
      IF (BLDFRA) THEN
C BUSCH
         A2=XLS + XCMS -XA
         XCM=AMIN1(XCM,A2)
         A2=YUS + YCMS -YA
         YCM=AMIN1(YCM,A2)
      ENDIF
      SPREIZ=.FALSE.
      XMAXF=XMAX
      IF (ISK.GE.0) GOTO 8
      CALL GRLGAR(SX,KSX,SXX,ALO1,FAKTX,IDECAX)
      GOTO 9
    8 CALL GRSCLA(SX,KSX,SXX,POTX,FAKTX)
      IF(XMAXF.NE.XMAX) SPREIZ=.TRUE.
    9 IF(.NOT.SPREIZ) GOTO 92
      WRITE(IO,91)
   91 FORMAT(//,' BILD : X-INTERVALL GESPREIZT')
   92 SPREIZ=.FALSE.
      YMAXF=YMAX
      IF (JSK.GE.0) GOTO 10
      CALL GRLGAR(SY,KSY,SYY,ALO2,FAKTY,IDECAY)
      GOTO 11
   10 CALL GRSCLA(SY,KSY,SYY,POTY,FAKTY)
      IF(YMAX.NE.YMAXF) SPREIZ=.TRUE.
   11 IF(.NOT.SPREIZ) GOTO 112
      WRITE(IO,111)
  111 FORMAT(//,' BILD : Y-INTERVALL GESPREIZT')
  112 CONTINUE
      DO 113 I=1,14
         PISY(I)=14
  113 CONTINUE
      WRITE(IO,200) XMIN,XMAX,YMIN,YMAX
  200 FORMAT(2X,'BILD: XMIN',12X,'XMAX',12X,'YMIN',12X,'YMAX',/,5X,
     1       4(2X,E14.7))
C
C  ZEICHNEN DER BESCHRIFTUNGEN
C
      IF((ISK.EQ.0.OR.ISK.LE.-999).AND.(JSK.EQ.0.OR.JSK.LE.-999))GOTO 26
      IF (ISK.GE.0) GOTO 12
      IF (ISK.LE.-999) GOTO 13
      CALL GRLOG(SX,KSX,SXX,ALO1,1.,0.,.77,IDECAX)
      GOTO 13
   12 IF (ISK.EQ.0) GOTO 13
      CALL GRLIN(SX,KSX,SXX,1.,0.,.6)
      IF (POTX.EQ.0) GOTO 13
      HOBR=.4
      IF (MOD(ISK,100).EQ.1) HOBR=.3
      CALL GRCHRC(HOBR,0.,18)
      A1=2.3-7*HOBR
      CALL GRTXT(A1,-1.5,5,'*10**')
      CALL GRPCTR(POTX,LPIC,PICC)
      CALL GRTXTC(LPIC,PICC)
   13 IF (JSK.GE.0) GOTO 14
      IF (JSK.LE.-999) GOTO 15
      CALL GRLOG(SY,KSY,SYY,ALO2,0.,1.,.4,IDECAY)
      GOTO 15
   14 IF (JSK.EQ.0) GOTO 15
      CALL GRLIN(SY,KSY,SYY,0.,1.,.3)
      IF (POTY.EQ.0) GOTO 15
      HOBR=.4
      IF (MOD(JSK,100).EQ.1) HOBR=.3
      CALL GRCHRC(HOBR,90.,18)
      A2=2.3-7*HOBR
      CALL GRTXT(-1.1,A2,5,'*10**')
      CALL GRPCTR(POTY,LPIC,PICC)
      CALL GRTXTC(LPIC,PICC)
   15 CALL GRJMP(0.,0.)
C
C  ZEICHNEN DER ACHSEN
C
      XP=0.
      YP=0.
      IF (ISK.NE.0 .AND. ISK.GT.-999 ) GOTO 16
      CALL GRJMP(XCM,0.)
      GOTO 19
   16 IF (ISK.LT.0) GOTO 17
      DO 161 I=1,32
         AXZAX(I)=0
  161 CONTINUE
      PP3=-.2
      MO=MOD(ISK,100)
      DO 162 I=2,32,MO
         AXZAX(I)=1
  162 CONTINUE
      GOTO 18
   17 PP3=-.3
      CMX1=0.
      AXZAX(1)=-1
   18 CALL GRLNLG(0.,XCM,0.,XP,CMX1,CMX2,ALO1,1.,0.,PP3,-.2,AXZAX)
   19 IF (JSK.NE.0 .AND. JSK.GT.-999) GOTO 20
      CALL GRJMP(XCM,YCM)
      GOTO 23
   20 IF (IABS(ISK).GT.100 .OR. IABS(JSK).GT.100) CALL GRSPTS(15)
      IF (JSK.LT.0) GOTO 21
      DO 201 J=1,32
         AXZAY(J)=0
  201 CONTINUE
      QQ3=-XCM
      IF (JSK.LE.200) QQ3=.2
      PPP3=.2
      QQQ3=.2
      IF (JSK.GT.100) QQQ3=-XCM
      MO=MOD(JSK,100)
      DO 202 J=2,32,MO
         AXZAY(J)=1
  202 CONTINUE
      GOTO 22
   21 PPP3=.3
      QQQ3=.3
      QQ3=.2
      IF (JSK.LT.-200) QQ3=-XCM
      IF (JSK.LT.-100) QQQ3=-XCM
      CMY1=0.
      AXZAY(1)=-1
   22 CALL GRLNLG(0.,XCM,YCM,YP,CMY1,CMY2,ALO2,0.,1.,QQQ3,QQ3,AXZAY)
   23 IF (ISK.NE.0 .AND. ISK.GT.-999) GOTO 24
      CALL GRJMP(0.,YCM)
      GOTO 25
   24 IF (IABS(ISK).GT.100) PP3=YCM
      QQ3=.2
      IF (IABS(ISK).GT.200) QQ3=-YCM
      CALL GRLNLG(XCM,0.,YCM,-XP,0.,-CMX2,ALO1,1.,0.,-PP3,QQ3,AXZAX)
   25 IF (IABS(ISK).GT.100 .OR. IABS(JSK).GT.100) CALL GRSPTS(18)
      IF (JSK.NE.0 .AND. JSK.GT.-999)
     1CALL GRLNLG(YCM,0.,0.,-YP,0.,-CMY2,ALO2,0.,1.,-PPP3,-.2,AXZAY)
C
C  ZEICHNEN DER KURVEN IN EINEM BILD
C  SYMBOL UND KRAFT(FARBE) DER KURVEN
C
   26 IF (NBILD.GE.NKURV) NBILD=NKURV
      LL=0
c---- loop mit i von 1 bis NBILD
      I=1
   42 IF (I.GT.NBILD) GOTO 421
      ANFANG=.TRUE.
      VORAUS=.TRUE.
      IST1=0
   27 READ(ISYS) IST,ISY,A,B
      IF(IST1.NE.0) GOTO 36
      IST1=IST
      IF (I.LE.14) PISY(I)=ISY
      JSY=MOD(ISY,2000)
      I1=JSY/200+1
      CALL GRNWPN(I1)
      IF (ISY.LT.2000) GOTO 300
      I1=18+ISY/2000*2
      CALL GRSPTS(I1)
      CALL GRCHRC(.5,0.,I1)
      GOTO 31
  300 CALL GRSPTS(18)
      CALL GRCHRC(.3,0.,18)
   31 ISY=MOD(ISY,200)
      IF(ISY.GT.13) GOTO 33
      IF(ISY.GE.0) GOTO 32
      ISY=-ISY
      ARTLAB=1
      IFA=ISY-10
      IF (IFA.LT.4) IFA=4
      IF (IFA.GT.19) IFA=19
      IFA=IFA/4
      CALL GRNWPN(IFA)
      CALL GRSPTS(ISY)
      GOTO 35
   32 ARTLAB=2
      GOTO 35
   33 CONTINUE
      IF(.NOT.(ISY.LT.100.OR.(ISY.GT.113.AND.ISY.LE.118))) GOTO 34
      ARTLAB=3
      IF(ISY.GT.113)
     1   CALL GRDSH(SH(1,ISY-113),SH(2,ISY-113),SH(3,ISY-113))
      GOTO 35
   34 ARTLAB=4
   35 IP=ISY
      MARKER=.FALSE.
      IF (ISY.GE.100) THEN
         IP=ISY-100
         IF ( IP.LE.13) THEN
C           erster Punkt als Marker
            MARKER=.TRUE.
         ENDIF
      ENDIF
      IF ( IP.LE.13) THEN
C        ANPSSUNG AN GTS GRAL MARKER
         IP = SYMVEC(IP)
      ENDIF
C
      GOTO 37

   36 XP=XL
      YP=YL
      XQ=A(1)
      YQ=B(1)
      RETLAB = 1
      GOTO 500

   37 XL=A(100)
      YL=B(100)
      L=IST1
      IF(IST1.GT.100) L=100
      IF(L.EQ.1 .AND. ISY.LE.13) GOTO 371
      K=2
      KP=1
      GOTO 38

  371 K=1
      KP=2
      A(2)=XMAX+(XMAX-XMIN)*.5
      B(2)=YMAX+(YMAX-YMIN)*.5
      BAUS=21
      ANFANG=.FALSE.
   38 IF(K.GT.L) GOTO 40
      XP=A(KP)
      YP=B(KP)
      XQ=A(K)
      YQ=B(K)
      RETLAB = 2
      GOTO 500

   39 ANFANG=.FALSE.
      KP=K
      K=K+1
      GOTO 38

   40 IF (L.GE.IST1) GOTO 41
      IST1=IST1-100
      GOTO 27

   41 CALL GRSPTS(18)
      CALL GRNWPN(1)
      CALL GRDSH(1.,0.,1.)
      I=I+1
      GOTO 42
C---- End loop mit i -----------
  421 CONTINUE
      IF(LL.LE.0) GOTO 44
      WRITE(IO,43) LL
   43 FORMAT(1X,I10,' PUNKTE LIEGEN AUSSERHALB DER GRENZEN')
   44 NKURV=NKURV-NBILD
      CILBER(8)=XCM
      IF (NKURV.NE.0) THEN
         GOTO 777
      ENDIF
      CILBER(1)=XCM
      CILBER(2)=YCM
      REWIND ISYS
      RETURN
  500 IF(XQ.GE.XMIN) GOTO 501
      BLERI=2
      GOTO 503
  501 IF(XQ.LE.XMAX) GOTO 502
      BLERI=3
      GOTO 503
  502 BLERI=1
  503 IF(YQ.GE.YMIN) GOTO 504
      BTOBO=5
      GOTO 506
  504 IF(YQ.LE.YMAX) GOTO 505
      BTOBO=7
      GOTO 506
  505 BTOBO=1
  506 BSEPO=BLERI*BTOBO
      IF(BSEPO.NE.1) LL=LL+1
      IF(.NOT.ANFANG) GOTO 520
  510 IF(XP.GE.XMIN) GOTO 511
      BLERI=2
      GOTO 513
  511 IF(XP.LE.XMAX) GOTO 512
      BLERI=3
      GOTO 513
  512 BLERI=1
  513 IF(YP.GE.YMIN) GOTO 514
      BTOBO=5
      GOTO 516
  514 IF(YP.LE.YMAX) GOTO 515
      BTOBO=7
      GOTO 516
  515 BTOBO=1
  516 BFIPO=BLERI*BTOBO
      IF(BFIPO.NE.1) LL=LL+1
      GOTO 521
  520 BFIPO=BAUS
  521 BAUS=BSEPO
      AAUS=BFIPO
      EXCH=.FALSE.
  522 IF(BFIPO.EQ.1 .AND. BSEPO.EQ.1) GOTO 540
      IF(BFIPO/2*2.EQ.BFIPO .AND. BSEPO/2*2.EQ.BSEPO) GOTO 5000
      IF(BFIPO/3*3.EQ.BFIPO .AND. BSEPO/3*3.EQ.BSEPO) GOTO 5000
      IF(BFIPO/5*5.EQ.BFIPO .AND. BSEPO/5*5.EQ.BSEPO) GOTO 5000
      IF(BFIPO/7*7.EQ.BFIPO .AND. BSEPO/7*7.EQ.BSEPO) GOTO 5000
      IF(BFIPO.NE.1) GOTO 523
      BHEPO=BFIPO
      RXH=XP
      RYH=YP
      BFIPO=BSEPO
      XP=XQ
      YP=YQ
      BSEPO=BHEPO
      XQ=RXH
      YQ=RYH
      EXCH=.NOT.EXCH
  523 IF(BFIPO/2*2.NE.BFIPO) GOTO 524
      YP=YP+(YQ-YP)*(XMIN-XP)/(XQ-XP)
      XP=XMIN
      GOTO 530
  524 IF(BFIPO/3*3.NE.BFIPO) GOTO 525
      YP=YP+(YQ-YP)*(XMAX-XP)/(XQ-XP)
      XP=XMAX
      GOTO 530
  525 IF(BFIPO/5*5.NE.BFIPO) GOTO 526
      XP=XP+(XQ-XP)*(YMIN-YP)/(YQ-YP)
      YP=YMIN
      GOTO 530
  526 IF(BFIPO/7*7.NE.BFIPO) GOTO 530
      XP=XP+(XQ-XP)*(YMAX-YP)/(YQ-YP)
      YP=YMAX
  530 IF(XP.GE.XMIN) GOTO 531
      BLERI=2
      GOTO 533
  531 IF(XP.LE.XMAX) GOTO 532
      BLERI=3
      GOTO 533
  532 BLERI=1
  533 IF(YP.GE.YMIN) GOTO 534
      BTOBO=5
      GOTO 536
  534 IF(YP.LE.YMAX) GOTO 535
      BTOBO=7
      GOTO 536
  535 BTOBO=1
  536 BFIPO=BLERI*BTOBO
      GOTO 522
  540 IF(.NOT.EXCH) GOTO 541
      RXH=XP
      RYH=YP
      XP=XQ
      YP=YQ
      XQ=RXH
      YQ=RYH
  541 IF(.NOT.VORAUS) GOTO 550
      IF(ISK.GE.0) GOTO 542
      XP=ALOG10(XP/XMIN)*FAKTX
      GOTO 543
  542 XP=(XP-XMIN)*FAKTX
  543 IF(JSK.GE.0) GOTO 544
      YP=ALOG10(YP/YMIN)*FAKTY
      GOTO 545
  544 YP=(YP-YMIN)*FAKTY
c 545 CALL GRJMP(XP,YP)
c     Busch 21.5.91
c     Startpunkt mit/ohne Marker
  545 IF ( MARKER) THEN
          CALL GRJMPS(XP,YP,IP)
       ELSE
          CALL GRJMP(XP,YP)
      ENDIF
  550 IF(ISK.GE.0) GOTO 551
      XQ=ALOG10(XQ/XMIN)*FAKTX
      GOTO 552
  551 XQ=(XQ-XMIN)*FAKTX
  552 IF(JSK.GE.0) GOTO 553
      YQ=ALOG10(YQ/YMIN)*FAKTY
      GOTO 560
  553 YQ=(YQ-YMIN)*FAKTY
  
 560  select case(ARTLAB)

      case(1)

      IF(ANFANG .AND. AAUS.EQ.1) CALL GRDRW(XP+.02,YP+.02)
      IF(BAUS.NE.1) GOTO 5000
      CALL GRJMP(XQ,YQ)
      CALL GRDRW(XQ+.02,YQ+.02)

      case(2)

      IF(ANFANG .AND. AAUS.EQ.1) CALL GRJMPS(XP,YP,IP)
      IF(BAUS.EQ.1) CALL GRJMPS(XQ,YQ,IP)

      case(3)

      CALL GRDRW(XQ,YQ)

      case(4)

      CONTINUE
      CALL GRDRWS(XQ,YQ,IP)

      end select

 5000 VORAUS=.FALSE.
      IF(BAUS.NE.1) VORAUS=.TRUE.
      
      if ( RETLAB == 1 ) then
          goto 37
      elseif ( RETLAB == 2 ) then
          goto 39
      endif
 

      END
