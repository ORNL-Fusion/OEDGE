CC@PROCESS OPT(3) NOSDUMP NOGOSTMT FIPS(F) IL(DIM) MAP
C******************************************************* GRHHNL ********
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C NDIMX =   DIE WIRKLICHE GESAMT-DIMENSION IN X-RICHTUNG
C NX    =   WIEVIELE RASTERSTELLEN IN X-RICHTUNG GENOMMEN WERDEN SOLLEN
C NY    =   WIEVIELE RASTERSTELLEN IN Y-RICHTUNG GENOMMEN WERDEN SOLLEN
C OPTION =  MAX 8 ZEICHEN; ENTHAELT EIN K, WENN ISTAX USW. BENUTZT!
C ISTAX,INCX,ISTAY,INCY WIE BEI GRVFLD UND GRGFLD
C----------------------------------------------------------------------
C     CMS VERSION ALS KOMPRIMIERENDES ZWISCHENPROGRAMM
C
C 09.08.90 bei negativem INTEG kann man mit LT2 Intensitaeten angeben.
C 25.04.88 POLARKOORDINATEN COMMON /GRPOLAR/ IPOLAR. >0 GRAD. <0 BOGEN.
C 29.02.88 OPTION BIS ZU 8 CHARACTER, INCREMENTS ISTAX USW., KOMPRIMIERE
C          FUER ZWISCHENPROGRAMM. GROTEN
C 17.04.89 XA, Dynamic Common GRDYN, Dimensionen vergr.(10.5.90 weg)
C 10.12.90 BUSCH, PP(9) DURCH COMMON GRPIC ERSETZT XMXCMQ
C 10. 1.91 BUSCH, 39.5 DURCH XMXCMQ ERSETZT
C 04.07.91 Groten : Fehler bei negativem INT(N) korrigiert.
C 20.07.91 Groten : Fehler bei der Farbreihenfolge korrigiert
C 19.08.91 BUSCH  : NUMARG eliminiert - es muessen 14 Argumente
C                   kodiert werden
C 21.11.91 Groten : Farben von GRHHNL und von GRGFLD einander angepasst
C UPDATE 20. 2.1991 Busch , GKS Treiber erlaubt max. 300 Werte pro Kurve
C                   GPL Aufruf durch GRLN Aufruf ersetzt
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C 17. 1.92 Groten : siehe GRHHFL; IFC
C 11. 5.92 GROTEN : Krumml.transf. Koord.; IPOLAR=987654321-987654323
C  5. 8. 93 Busch   COMMON GRGKS und GRWK wegen Identifikation der
C                   Unixversion, dort keine GCRSG, da die Graphik
C                   nur in einem Segment (war notwendig wegen
C                   GRHPAN im VM
C 3.12.93 Busch Farben 3+4 , 5+6 fuer UNIX falsch
C 3.03.98 Groten Flaechenfuellung Abfrage geaendert (SEL), wegen IEEE,EX1,EX2
C----------------------------------------------------------------------------
      SUBROUTINE GRHHNL(NDIMX,TAB,NX,XX,NY,YY,L1,WERT,INTEG,OPTION,
     +                  ISTAX,INCX,ISTAY,INCY)
      IMPLICIT NONE
C     .. Previously untyped names ..
      INTEGER IFC,MAL,MAL1,MAN,MATMAX,MAXINX,MAXINY,INCX,INCY,ISTAX,
     +        ISTAY,L1,NDIMX,NX,NY
      REAL TAB,WERT,XX,YY
      INTEGER INTEG
      REAL DU1,DU2,DU3,AUX,AUY,AX,AY,HV,HX,HY
      INTEGER I,IND,J,K,L7,LENG
C     ..
C---- PC: IFC=72; CMS: IFC=512  in GRHHFL gleichsetzen!
      PARAMETER (IFC=512,MAN=IFC,MATMAX=MAN*MAN,MAXINX=MAN,MAXINY=MAN,
     +          MAL=3*MAXINX+3*MAXINY,MAL1=(MAL*6)/5)
      CHARACTER (len=*) OPTION
      CHARACTER (len=7)OPT1
      LOGICAL KOMPR
      DIMENSION TAB(NDIMX,*),XX(*),YY(*),WERT(*),INTEG(*)
      CHARACTER (len=1) MRK(MATMAX)
      INTEGER IGRSHW,ERRUN,IGR3,SHOWPR,IGR3PL
C---- COMMON 'GRDYN' WIRD AUCH IN GRHHFL, CRMENU UND CRPLOT VERWENDET.
      COMMON /GRDYN/HV(MATMAX),MRK,AX(MAL),AY(MAL),AUX(MAL1),AUY(MAL1),
     +       HX(MAXINX),HY(MAXINY),DU1,DU2,DU3
CDEC$ PSECT /GRDYN/ NOSHR
      COMMON /GRCOM1/IGRSHW,ERRUN,IGR3,SHOWPR,IGR3PL
CDEC$ PSECT /GRCOM1/ NOSHR
      SAVE  /GRDYN/,/GRCOM1/

      IF (NX.GT.MAXINX .OR. NY.GT.MAXINY .OR. NX*NY.GT.MATMAX) THEN
          WRITE (6,FMT=*) 'GRHHNL :'
          WRITE (6,FMT=*)
     +      ' Indexmenge in X- oder Y-Richtung ist zu gross!'
          WRITE (6,FMT=*) ' Maximale Groessen fuer NX , NY : ',MAXINX,
     +      MAXINY
          WRITE (6,FMT=*) ' Maximale Groesse  fuer NX * NY : ',MATMAX
          RETURN
      END IF
C---- NUMARG NOT ON PC'S, NOT STANDARD FORTRAN. ONLY TO REDUCE ERRORS.
      LENG = MIN(8,LEN(OPTION))
      K = INDEX(OPTION(1:LENG),'K')
C     KOMPR=NUMARG().EQ.14 .AND. K.GT.0
      KOMPR = K .GT. 0
      IF (KOMPR) THEN
          OPT1 = ' '
          L7 = MIN(LENG,7)
          IF (K.EQ.0) THEN
              OPT1(1:L7) = OPTION(1:LENG)
          ELSE IF (K.EQ.1) THEN
              OPT1(1:L7) = OPTION(2:LENG)
          ELSE IF (K.EQ.LENG) THEN
              OPT1(1:L7) = OPTION(1:K-1)
          ELSE
              OPT1(1:L7) = OPTION(1:K-1)//OPTION(K+1:LENG)
          END IF
          IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
              CALL GRHL(NDIMX,TAB(ISTAX,ISTAY),NX,XX(ISTAX),NY,
     +                  YY(ISTAY),L1,WERT,INTEG,OPT1,MRK,AX,AY,AUX,AUY,
     +                  MAL)
          ELSE
              IND = 0
              DO 17 J = ISTAY,ISTAY + (NY-1)*INCY,INCY
                  DO 16 I = ISTAX,ISTAX + (NX-1)*INCX,INCX
                      IND = IND + 1
                      HV(IND) = TAB(I,J)
   16             CONTINUE
   17         CONTINUE
              IND = 0
              DO 18 I = ISTAX,ISTAX + (NX-1)*INCX,INCX
                  IND = IND + 1
                  HX(IND) = XX(I)
   18         CONTINUE
              IND = 0
              DO 19 I = ISTAY,ISTAY + (NY-1)*INCY,INCY
                  IND = IND + 1
                  HY(IND) = YY(I)
   19         CONTINUE
              CALL GRHL(NX,HV,NX,HX,NY,HY,L1,WERT,INTEG,OPT1,MRK,AX,AY,
     +                  AUX,AUY,MAL)
          END IF
      ELSE
          CALL GRHL(NDIMX,TAB,NX,XX,NY,YY,L1,WERT,INTEG,OPTION,MRK,AX,
     +              AY,AUX,AUY,MAL)
      END IF
      END
CC@PROCESS OPT(3) NOSDUMP NOGOSTMT FIPS(F)
C     IM WESENTLICHEN DAS ALTE GRHHNL; 'MRK' JETZT PARAMETER
      SUBROUTINE GRHL(N1,TAB,M1,XX,M2,YY,L1,WERT,INT,OPTION,MRK,X,Y,XAU,
     +                YAU,MAL)
      IMPLICIT NONE
C     .. Previously untyped names ..
      REAL EX1,EX2,PIFAC,ZWOPI
      INTEGER MAXDEF,MAXHHL,L1,M1,M2,MAL,N1
      REAL XCURR,XDCPIC,XMAXDC,XMXCMQ,XUNITS,YCURR,YDCPIC,YMXCMQ,YUNITS
      INTEGER IOB,IOBFA,IPOLAR,LAGE,MAXPKT,NCHAR,NRLN,NRMR,NSCLC,NSCLP,
     +        NSCLV,NWERT
      REAL PP,WER
      INTEGER IFA
      REAL ABS1A,ABS2A,ABS3A,CP,D1,D2,DX,DY,F,FF5,FF6,FF7,FF8,H,HX,HY,
     +     OO3,OO4,PIF,PP5,PP6,PP7,PP8,QQ3,QQ4,RMU,SP,U,U1,U2,V,V1,V2,
     +     W11,W12,W21,W22,WEMI,X1,X2,XA,XE,XH,XK,XLI,XMAXCM,XRE,XS,XUP,
     +     Y1,Y2,YA,YE,YK,YLI,YMAXCM,YRE,YS,YUP,Z0,Z11,Z12,Z21,Z22,ZMAX,
     +     ZMIN
      INTEGER I,IALH,IALV,IDASH,IDUM,IER,IERR,II,III,IIMAX,IIMIN,IM,
     +        IMAX,IMAXM1,IN,INCOL,INCOT,INCX,INCY,INTA,INTLIN,IPO,IRX,
     +        IRY,J,JFA,JJ,JJJ,JJMAX,JJMIN,JM,JMAX,JMAXM1,K,KIND,LAB,
     +        LANG,LP,LTYPE,MAL1,MAXX,MAXY,MIM,N,NA,NAME,NAME1,NEND,
     +        NMAX,NOP,NRF,NSG,NX,NY
C     ..
      PARAMETER (MAXDEF=28,EX1=-75.7501E20,EX2=-75.7499E20)
      PARAMETER (ZWOPI=3.14159265358979*2,PIFAC=ZWOPI/360.)
      PARAMETER (MAXHHL=100)
      COMMON /GRPP/PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRWK/IXGKS,IWK,IWKWIS
CDEC$ PSECT /GRWK/ NOSHR
      COMMON /GRGKS/GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      COMMON /GRHHPN/IOB,NCHAR,LAGE,SKAL,DREHF
CDEC$ PSECT /GRHHPN/ NOSHR
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /SCALE/XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
      COMMON /GRWEFA/IOBFA,NWERT,WER(MAXHHL),IFA(MAXHHL)
CDEC$ PSECT /GRWEFA/ NOSHR
      COMMON /GRPOLAR/IPOLAR
CDEC$ PSECT /GRPOLAR/ NOSHR
      COMMON /GRPIC/FLPIC,NSCLC,NSCLV,NSCLP,RAHMEN,XMXCMQ,YMXCMQ,XDCPIC,
     +       YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /GRCOM1/IGRSHW,ERRUN,IGR3,SHOWPR,IGR3PL
CDEC$ PSECT /GRCOM1/ NOSHR
      SAVE /GRPP/,/GRWK/,/GRGKS/,/GRHHPN/,/GRREST/,/SCALE/,/GRWEFA/,
     $      /GRPOLAR/, /GRPIC/,/GRCOM1/
      EQUIVALENCE (PP(13),INTLIN)

C     REAL TAB(N1, * ),WERT( * ),XX( * ),YY( * ),RAX(5),RAY(5)
      REAL TAB(N1,M2),WERT(L1),XX(M1),YY(M2),RAX(5),RAY(5)
      CHARACTER (len=1) MRK(M1,M2)
      REAL X(MAL),Y(MAL),XAU(MAL+MAL/5),YAU(MAL+MAL/5),TRAN(2,3)
      CHARACTER (len=22) KETTE1,KETTE3
      CHARACTER (len=29) KETTE2
      CHARACTER (len=10) KETTE4
      CHARACTER (len=*) OPTION
      CHARACTER (len=7) MARK
      CHARACTER GKSTYP*7
      INTEGER RAHMEN,FLPIC,IXGKS,IWK,IWKWIS
      INTEGER INT(L1),FARBE,SKAL,DREHF,INFA(31),ICOL(16)
      INTEGER IGRSHW,ERRUN,IGR3,SHOWPR,IGR3PL,IPOLAB
      REAL XLN(300),YLN(300),XMR(300),YMR(300)
      CHARACTER (len=1) MR
      LOGICAL BIT(7),NOLINE,OBZA,LOGI
  
     
      DATA KETTE2/'LAGE =(          ,          )'/
      DATA KETTE1/'MIN: WERT =           '/
      DATA KETTE3/'MAX: WERT =           '/
      DATA MARK/'HRNEIGV'/
      DATA INFA/1,2,3,4,5,6,7,8,1,1,1,1,1,1,1,1,1,2,2,2,2,4,4,4,4,3,3,3,
     +     3,3,3/
      DATA ICOL/1,2,4,3,6,5,7,8,9,10,11,12,13,14,15,16/
C MB 3.12.93
      IF (GKSTYP.EQ.'GLIGKS' .OR. GKSTYP.EQ.'XGKS') THEN
          ICOL(3) = 3
          ICOL(4) = 4
          ICOL(5) = 5
          ICOL(6) = 6
      END IF
C
      MAL1 = MAL*1.2
      PIF = PIFAC
      IPOLAB = ABS(IPOLAR)
      IF (IPOLAB.EQ.1234567890) THEN
          IF (IPOLAR.LT.0) THEN
              PIF = 1.
          END IF
      END IF
      XMAXCM = XMXCMQ
      YMAXCM = YMXCMQ

      IF (NRLN.GT.1) THEN
          CALL GPL(NRLN,XLN,YLN)
          NRLN = 0
      END IF
      IF (NRMR.GT.0) THEN
          CALL GPM(NRMR,XMR,YMR)
          NRMR = 0
      END IF
      CALL GQPLCI(IERR,INCOL)
      CALL GQTXCI(IERR,INCOT)
      CALL GQTXAL(IERR,IALH,IALV)
      CALL GSTXAL(0,4)
      H = PP(14)
      DY = 1.7*H
      DX = 0.
      HY = 0.
      HX = H
      IMAX = M1
      JMAX = M2
      NMAX = L1
C---- die naechsten Zuweisungen nur um den Compiler zu beruhigen.
      JJ = 0
      NAME = 0
      NAME1 = 0
      III = 0
      JJJ = 0
      JMAXM1 = 0
      U1 = 0.
      U2 = 0.
      V2 = 0.
      W12 = 0.
      W22 = 0.
      FF5 = 0.
      FF6 = 0.
      FF7 = 0.
      FF8 = 0.
C---- Ende nur Beruhigung
      PP5 = PP(5)
      PP6 = PP(6)
      PP7 = PP(7)
      PP8 = PP(8)
      ABS1A = PP(10)
      ABS2A = PP(11)
      ABS3A = PP(12)
      INTA = INTLIN
      IF (IPOLAB.NE.1234567890) THEN
          MAXY = MAXDEF
          MAXX = MAXY* (PP(3)-PP(1))/ (PP(4)-PP(2))
      ELSE
          MAXY = 3*MAXDEF/4
          MAXX = 44*MAXY/7* (XX(IMAX)-XX(1))/ (ZWOPI/PIF)
      END IF
CCC   'SOLID' LINE
      PP(10) = 1.
      PP(11) = 0.
      PP(12) = 1.
      CALL GSLN(1)
      DO 10 I = 1,7
          BIT(I) = .FALSE.
   10 CONTINUE
      DO 12 I = 1,LEN(OPTION)
          DO 11 J = 1,7
              IF (OPTION(I:I).EQ.MARK(J:J)) BIT(J) = .TRUE.
   11     CONTINUE
   12 CONTINUE

C     INTERAKTIV UEBER PANEL

      IF (IGR3.EQ.1) THEN
          CALL GQOPSG(IERR,NAME)
          IF (IERR.EQ.0) CALL GCLSG
C        QUERY NUMBER OF OPEN WK'S
C        NUR EINE IST OFFEN
          CALL GQOPWK(1,IERR,NOP,IWK)
C        QUERY NUMBER OF SEGMENTS ON WK
          CALL GQSGWK(IWK,1,IERR,NSG,NAME)
          NAME = 1
          DO 444 I = 1,NSG
              CALL GQSGWK(IWK,I,IERR,IDUM,NAME1)
C           WRITE(93,*) '0 NAME1=', NAME1
              IF (NAME1.GE.NAME) NAME = NAME + 1
  444     CONTINUE
C     MIT GR3PAN
          IF (NAME1.EQ.NAME) NAME = NAME + 1
C        WRITE(93,*) '1 NAME=', NAME
          IF ((IXGKS.NE.1234567890) .AND.
     +        (GKSTYP.NE.'GLIGKS')) CALL GCRSG(NAME)
c            CALL GCRSG(NAME)
      END IF
C-----------------------------------------------------------------------
C---------- BESCHRIFTUNG -----------------------------------------------
C-----------------------------------------------------------------------
      IF (BIT(4) .OR. BIT(5)) THEN
          CALL GSCHXP(1.)
          IF (IOB.EQ.1234567890) THEN
              IF (LAGE.LE.3) THEN
                  DX = 0.
                  DY = 1.7*H
                  HX = H
                  HY = 0.
                  XUP = 0.
                  YUP = 1.
                  IF (LAGE.EQ.1) THEN
                      XS = PP(3) + NCHAR*H
                      YS = PP(4) - H
                  ELSE IF (LAGE.EQ.2) THEN
                      XS = PP(1) - (NCHAR+21)*H
                      YS = PP(4) - H
                  ELSE
                      XS = (PP(1)+PP(3)-21.*H)/2.
                      YS = PP(2) - (NCHAR+1)*H
                  END IF
              ELSE
                  DX = 1.7*H
                  DY = 0.
                  HX = 0.
                  HY = H
                  XUP = -1.
                  YUP = 0.
                  IF (LAGE.EQ.4) THEN
                      XS = PP(1) + H
                      YS = PP(4) + NCHAR*H
                  ELSE IF (LAGE.EQ.5) THEN
                      XS = PP(1) + H
                      YS = PP(2) - (NCHAR+21)*H
                  ELSE
                      XS = PP(3) + (NCHAR+1)*H
                      YS = (PP(2)+PP(4)-21.*H)/2.
                  END IF
              END IF
          ELSE
              XUP = 0.
              YUP = 1.
              XS = PP(3) + 3.
              YS = PP(4) - H
          END IF
          PP(5) = PP(1)
          PP(6) = PP(2)
          PP(7) = PP(3)
          PP(8) = PP(4)
          XUNITS = 1.
          YUNITS = 1.
          XRE = XMAXCM
          YRE = YMAXCM
          CALL GSWN(1,0.,XRE,0.,YRE)
          CALL GSCHUP(XUP,YUP)
          CALL GSCHH(PP(14))
      END IF
C-----------------------------------------------------------------------
C---------- LAGE EINES MINIMUMS UND EINES MAXIMUMS (EXTREMA) -----------
C-----------------------------------------------------------------------
      IF (BIT(4)) THEN
          ZMAX = TAB(1,1)
          ZMIN = TAB(1,1)
          IIMIN = 1
          IIMAX = 1
          JJMAX = JMAX
          JJMIN = JMAX
          DO 523 I = 1,IMAX
              DO 522 J = 1,JMAX
                  IF (TAB(I,J).LT.EX1 .OR. TAB(I,J).GT.EX2) THEN
                      IF (ZMAX.LT.TAB(I,J)) THEN
                          ZMAX = TAB(I,J)
                          IIMAX = I
                          JJMAX = J
                      END IF
                      IF (ZMIN.GT.TAB(I,J)) THEN
                          ZMIN = TAB(I,J)
                          IIMIN = I
                          JJMIN = J
                      END IF
                  END IF
  522         CONTINUE
  523     CONTINUE
          IIMIN = IIMIN - 1
          IIMAX = IIMAX - 1
          JJMIN = JJMIN - 1
          JJMAX = JJMAX - 1
          CALL GRFTOC(ZMIN,KETTE1(13:),LANG)
          LANG = LANG + 12
          CALL GTX(XS,YS,KETTE1(:LANG))
          YS = YS - DY
          XS = XS + DX
          IF (IOB.EQ.1234567890 .AND. SKAL.EQ.3) THEN
              CALL GRFTOC(XX(IIMIN+1),KETTE2(8:),LP)
          ELSE
              CALL GRPCTR(IIMIN,LP,KETTE2(8:))
          END IF
          IPO = 8 + LP
          KETTE2(IPO:IPO) = ','
          IPO = IPO + 1
          IF (IOB.EQ.1234567890 .AND. SKAL.EQ.3) THEN
              CALL GRFTOC(YY(JJMIN+1),KETTE2(IPO:),LP)
          ELSE
              CALL GRPCTR(JJMIN,LP,KETTE2(IPO:))
          END IF
          IPO = IPO + LP
          KETTE2(IPO:IPO) = ')'
          CALL GTX(XS,YS,KETTE2(:IPO))
          YS = YS - 2.*DY
          XS = XS + 2.*DX
          CALL GRFTOC(ZMAX,KETTE3(13:),LANG)
          LANG = LANG + 12
          CALL GTX(XS,YS,KETTE3(:LANG))
          YS = YS - DY
          XS = XS + DX
          IF (IOB.EQ.1234567890 .AND. SKAL.EQ.3) THEN
              CALL GRFTOC(XX(IIMAX+1),KETTE2(8:),LP)
          ELSE
              CALL GRPCTR(IIMAX,LP,KETTE2(8:))
          END IF
          IPO = 8 + LP
          KETTE2(IPO:IPO) = ','
          IPO = IPO + 1
          IF (IOB.EQ.1234567890 .AND. SKAL.EQ.3) THEN
              CALL GRFTOC(YY(JJMAX+1),KETTE2(IPO:),LP)
          ELSE
              CALL GRPCTR(JJMAX,LP,KETTE2(IPO:))
          END IF
          IPO = IPO + LP
          KETTE2(IPO:IPO) = ')'
          CALL GTX(XS,YS,KETTE2(:IPO))
          YS = YS - 3.*DY
          XS = XS + 3.*DX
      END IF
C----------------------------------------------------------------------
C-------- HOEHEN (NIVEAUS) UND ZUGEHOERIGE LINIEN (FARBE, INTENSITAETEN)
C-----------------------------------------------------------------------
      IF (BIT(5)) THEN
          CALL GTX(XS,YS,' WERT')
          CALL GTX(XS+11.*HX,YS+11.*HY,'STRICH')
          YS = YS - 2.*DY
          XS = XS + 2.*DX
          X1 = NMAX
          IF (IOB.EQ.1234567890) THEN
              IF (LAGE.LE.3) THEN
                  X2 = YS/DY
              ELSE
                  X2 = (XMAXCM-XS)/DX
C              X2=(39.5-XS)/DX
              END IF
          ELSE
              X2 = YS/DY
          END IF
          NEND = MIN1(X1,X2)
          DO 531 N = 1,NEND
              CALL GRFTOC(WERT(N),KETTE4,LANG)
              CALL GTX(XS,YS,KETTE4(:LANG))
              X1 = XS + 11.*HX
              X2 = X1 + 6.*HX
              Y1 = YS + 11.*HY
              Y2 = Y1 + 6.*HY
              XLN(1) = X1 - .3*DX
              YLN(1) = Y1 + .3*DY
              IF (INT(N).GT.0) THEN
                  IN = MAX(MOD(INT(N),32),1)
                  FARBE = ICOL(INFA(IN))
                  IDASH = MAX(MOD(INT(N)/32,64),1)
              ELSE
                  FARBE = ICOL(INFA(MAX(MOD(ABS(INT(N)),32),1)))
                  IDASH = MAX(MOD(ABS(INT(N))/32,8),1)
                  IN = 18 + MOD(ABS(INT(N))/256,8)*2
              END IF
              INTLIN = IN
              CALL GRSPTS(INTLIN)
              CALL GSPLCI(FARBE)
              XLN(2) = X2 - .3*DX
              YLN(2) = Y2 + .3*DY
              CALL GRSCDL(2,XLN,YLN,IDASH,' ',MAL1,XAU,YAU)
              YS = YS - DY
              XS = XS + DX
              INTLIN = INTA
              CALL GRSPTS(INTLIN)
              CALL GSPLCI(1)
  531     CONTINUE

C     NUR FUER AUFRUF UEBER INTERAKTIVE PANEL
          IF (IGR3.EQ.1) THEN
              CALL GQOPSG(IERR,NA)
              IF (IERR.EQ.0) CALL GCLSG
              NAME = NAME + 1
C           WRITE(93,*) '2 NAME=', NAME
c           CALL GCRSG(NAME)
              IF (IXGKS.NE.1234567890 .AND.
     +            GKSTYP.NE.'GLIGKS') CALL GCRSG(NAME)
          END IF
      END IF
C----------------------------------------------------------------------
C------------------------- HOEHENLINIEN ODER RAHMEN -------------------
C----------------------------------------------------------------------
      IF (BIT(1) .OR. BIT(2) .OR. BIT(3) .OR. BIT(6) .OR. BIT(7)) THEN
          D1 = (PP(4)-PP(2))/ (PP(3)-PP(1))
          IF (IPOLAB.EQ.1234567890) THEN
C---------------------------------- POLARKOORDINATEN -------------------
              U = XX(1)*PIF
              XA = YY(1)*COS(U)
              XE = XA
              YA = YY(1)*SIN(U)
              YE = YA
              DO 1717 I = 1,IMAX
                  U = XX(I)*PIF
                  CP = COS(U)
                  SP = SIN(U)
                  DO 1716 J = 1,JMAX,JMAX - 1
                      XK = YY(J)*CP
                      YK = YY(J)*SP
                      XA = MIN(XA,XK)
                      XE = MAX(XE,XK)
                      YA = MIN(YA,YK)
                      YE = MAX(YE,YK)
 1716             CONTINUE
 1717         CONTINUE
              IF (D1.GT. (YE-YA)/ (XE-XA)) THEN
                  U = YA
                  YA = (U+YE-D1* (XE-XA))*.5
                  YE = (U+YE+D1* (XE-XA))*.5
              ELSE
                  U = XA
                  XA = (U+XE- (YE-YA)/D1)*.5
                  XE = (U+XE+ (YE-YA)/D1)*.5
              END IF
          ELSE IF (IPOLAB.GE.987654321 .AND. IPOLAB.LE.987654323) THEN
C---------------------------------- beliebiebige KOORDINATEN ----------
              CALL GRHHTR(XX(1),YY(1),XA,YA)
              XE = XA
              YE = YA
              DO 1719 I = 1,IMAX
                  DO 1718 J = 1,JMAX,JMAX - 1
                      CALL GRHHTR(XX(I),YY(J),XK,YK)
                      XA = MIN(XA,XK)
                      XE = MAX(XE,XK)
                      YA = MIN(YA,YK)
                      YE = MAX(YE,YK)
 1718             CONTINUE
 1719         CONTINUE
              DO 1713 J = 1,JMAX
                  DO 1712 I = 1,IMAX,IMAX - 1
                      CALL GRHHTR(XX(I),YY(J),XK,YK)
                      XA = MIN(XA,XK)
                      XE = MAX(XE,XK)
                      YA = MIN(YA,YK)
                      YE = MAX(YE,YK)
 1712             CONTINUE
 1713         CONTINUE
              IF (IPOLAB.EQ.987654322) THEN
                  IF (D1.GT. (YE-YA)/ (XE-XA)) THEN
                      U = YA
                      YA = (U+YE-D1* (XE-XA))*.5
                      YE = (U+YE+D1* (XE-XA))*.5
                  ELSE
                      U = XA
                      XA = (U+XE- (YE-YA)/D1)*.5
                      XE = (U+XE+ (YE-YA)/D1)*.5
                  END IF
              END IF
          ELSE
C---------- cartesische Koordinaten
              XA = XX(1)
              XE = XX(IMAX)
              YA = YY(1)
              YE = YY(JMAX)
          END IF
          IF (IPOLAB.NE.987654323) THEN
              PP(5) = XA
              PP(6) = YA
              PP(7) = XE
              PP(8) = YE
          ELSE
              XA = PP(5)
              YA = PP(6)
              XE = PP(7)
              YE = PP(8)
          END IF
          IF (XA.GT.XE) THEN
              XH = XA
              XA = XE
              XE = XH
          END IF
          IF (YA.GT.YE) THEN
              XH = YA
              YA = YE
              YE = XH
          END IF
          XUNITS = (XE-XA)/ (PP(3)-PP(1))
          YUNITS = (YE-YA)/ (PP(4)-PP(2))
          XLI = XA - PP(1)*XUNITS
          YLI = YA - PP(2)*YUNITS
          XRE = XE + (XMAXCM-PP(3))*XUNITS
          YRE = YE + (YMAXCM-PP(4))*YUNITS
          CALL GSWN(1,XLI,XRE,YLI,YRE)
          CALL GSCHH(PP(14)*YUNITS)
          CALL GSCHXP(XUNITS/YUNITS)
          IF (IOB.EQ.1234567890 .AND. DREHF.NE.0) THEN
              IF (DREHF.EQ.1) THEN
                  TRAN(1,1) = 0.
                  TRAN(1,2) = - (PP(3)-PP(1))/ (PP(4)-PP(2))
                  TRAN(1,3) = (PP(3)*PP(4)-PP(1)*PP(2))/ (PP(4)-PP(2))/
     +                        XMAXCM
                  TRAN(2,1) = (PP(4)-PP(2))/ (PP(3)-PP(1))
                  TRAN(2,2) = 0.
                  TRAN(2,3) = (PP(3)*PP(2)-PP(1)*PP(4))/ (PP(3)-PP(1))/
     +                        XMAXCM
              ELSE IF (DREHF.EQ.2) THEN
                  TRAN(1,1) = -1.
                  TRAN(1,2) = 0.
                  TRAN(1,3) = (PP(1)+PP(3))/XMAXCM
                  TRAN(2,1) = 0.
                  TRAN(2,2) = -1.
                  TRAN(2,3) = (PP(2)+PP(4))/XMAXCM
              ELSE
                  TRAN(1,1) = 0.
                  TRAN(1,2) = (PP(3)-PP(1))/ (PP(4)-PP(2))
                  TRAN(1,3) = (PP(1)*PP(4)-PP(2)*PP(3))/ (PP(4)-PP(2))/
     +                        XMAXCM
                  TRAN(2,1) = - (PP(4)-PP(2))/ (PP(3)-PP(1))
                  TRAN(2,2) = 0.
C              TRAN(2,3)=(PP(3)*PP(4)-PP(1)*PP(2))/(PP(3)-PP(1))/YMAXCM
                  TRAN(2,3) = (PP(3)*PP(4)-PP(1)*PP(2))/ (PP(3)-PP(1))/
     +                        XMAXCM
              END IF
              CALL GSSGT(NAME,TRAN)
          END IF
      ELSE
          GO TO 600
      END IF
C----------------------------------------------------------------------
C------------------------- HOEHENLINIEN -------------------------------
C----------------------------------------------------------------------
      IF (.NOT.BIT(1)) GO TO 311
      N = 0
   31 N = N + 1
      IF (N.GT.NMAX) GO TO 311
      Z0 = WERT(N)
      DO 14 J = 1,JMAX
          DO 13 I = 1,IMAX
              MRK(I,J) = '8'
   13     CONTINUE
   14 CONTINUE
      OBZA = .FALSE.
      IDASH = ABS(INT(N))/32
      IF (IDASH.GE.64) THEN
          OBZA = .TRUE.
      END IF
      IF (OBZA) THEN
          CALL GRFTOC(Z0,KETTE4,LANG)
      ELSE
          KETTE4 = ' '
          LANG = 1
      END IF
      IF (INT(N).GT.0) THEN
          IN = MAX(MOD(INT(N),32),1)
          FARBE = ICOL(INFA(IN))
          IDASH = MAX(MOD(INT(N)/32,64),1)
      ELSE
          FARBE = ICOL(INFA(MAX(MOD(ABS(INT(N)),32),1)))
          IDASH = MOD(ABS(INT(N))/32,8)
          IN = 18 + MOD(ABS(INT(N))/256,8)*2
      END IF
      INTLIN = IN
      CALL GRSPTS(INTLIN)
      CALL GSPLCI(FARBE)
      CALL GSTXCI(FARBE)

C------- UNTERER RAND

      IMAXM1 = IMAX - 1
      II = 0
   17 II = II + 1
      IF (II.GT.IMAXM1) GO TO 171
      J = 1
      Y2 = YY(J)
      I = II
      X1 = XX(I)
      X2 = XX(I+1)
      Z12 = TAB(I,J) - Z0
      Z22 = TAB(I+1,J) - Z0
      IF (SEL(Z12,Z22) .OR. MRK(I,J).NE.'8') GO TO 17
      IF (Z12.NE.Z22) THEN
          X(1) = (X1*Z22-X2*Z12)/ (Z22-Z12)
          IF (Z12.LT.Z22) THEN
              NRF = 2
          ELSE
              NRF = 1
          END IF
      ELSE
          X(1) = X1
          NRF = 9
      END IF
      Y(1) = Y2
      LAB = 1
      GO TO 1000

C---------- RECHTER RAND

  171 JMAXM1 = JMAX - 1
      JJ = 0
   18 JJ = JJ + 1
      IF (JJ.GT.JMAXM1) GO TO 181
      I = IMAXM1
      X1 = XX(I+1)
      J = JJ
      Y1 = YY(J)
      Y2 = YY(J+1)
      Z11 = TAB(I+1,J) - Z0
      Z12 = TAB(I+1,J+1) - Z0
      MR = MRK(I,J)
      IF (SEL(Z11,Z12).OR.MR.EQ.'9' .OR. MR.EQ.'0' .OR.
     +    MR.EQ.'4') GO TO 18
      IF (Z11.NE.Z12) THEN
          Y(1) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
          IF (Z11.LT.Z12) THEN
              NRF = 4
          ELSE
              NRF = 3
          END IF
      ELSE
          Y(1) = Y1
          NRF = 10
      END IF
      X(1) = X1
      LAB = 2
      GO TO 1000

C---------- OBERER RAND
  181 III = 0
   19 III = III + 1
      IF (III.GT.IMAXM1) GO TO 191
      II = IMAXM1 + 1 - III
      J = JMAX - 1
      Y1 = YY(J+1)
      I = II
      X1 = XX(I)
      X2 = XX(I+1)
      Z11 = TAB(I,J+1) - Z0
      Z21 = TAB(I+1,J+1) - Z0
      MR = MRK(I,J)
      IF (SEL(Z11,Z21).OR. MR.EQ.'9' .OR. MR.EQ.'2' .OR.
     +    MR.EQ.'6') GO TO 19
      IF (Z11.NE.Z21) THEN
          X(1) = (X1*Z21-X2*Z11)/ (Z21-Z11)
          IF (Z11.LT.Z21) THEN
              NRF = 6
          ELSE
              NRF = 5
          END IF
      ELSE
          X(1) = X2
          NRF = 11
      END IF
      Y(1) = Y1
      LAB = 3
      GO TO 1000

C---------- LINKER RAND

  191 JJJ = 0
   20 JJJ = JJJ + 1
      IF (JJJ.GT.JMAXM1) GO TO 201
      JJ = JMAXM1 + 1 - JJJ
      I = 1
      X2 = XX(I)
      J = JJ
      Y1 = YY(J)
      Y2 = YY(J+1)
      Z21 = TAB(I,J) - Z0
      Z22 = TAB(I,J+1) - Z0
      MR = MRK(I,J)
      IF (SEL(Z21,Z22).OR. MR.EQ.'9' .OR. MR.EQ.'1' .OR.
     +    MR.EQ.'5') GO TO 20
      IF (Z21.NE.Z22) THEN
          Y(1) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
          IF (Z21.GT.Z22) THEN
              NRF = 7
          ELSE
              NRF = 8
          END IF
      ELSE
          NRF = 12
          Y(1) = Y2
      END IF
      X(1) = X2
      LAB = 4
      GO TO 1000

C---------- INNERES GEBIET

  201 U2 = XX(1)
      II = 0
   30 II = II + 1
      IF (II.GT.IMAXM1) GO TO 31
      U1 = U2
      U2 = XX(II+1)
      W12 = TAB(II,1) - Z0
      W22 = TAB(II+1,1) - Z0
      V2 = YY(1)
      JJ = 0
   29 JJ = JJ + 1
      IF (JJ.GT.JMAXM1) GO TO 30
      W11 = W12
      W21 = W22
      W12 = TAB(II,JJ+1) - Z0
      W22 = TAB(II+1,JJ+1) - Z0
      V1 = V2
      V2 = YY(JJ+1)
      IF (SEL(W22,W21).OR.SEL(W22,W12)
     + .OR. MRK(II,JJ).EQ.'9' .OR.
     + (SEL(W11,W22).AND.MRK(II,JJ).NE.'1'.AND.
     + MRK(II,JJ).NE.'7'))  GO TO 29
      I = II
      J = JJ
      X1 = U1
      X2 = U2
      Z11 = W12
      Z21 = W22
      Y(1) = V2
      Y1 = V2
      IF (Z11.NE.Z21) THEN
          X(1) = (X1*Z21-X2*Z11)/ (Z21-Z11)
          IF (Z11.LT.Z21) THEN
              NRF = 6
          ELSE
              NRF = 5
          END IF
      ELSE
          X(1) = X2
          NRF = 11
      END IF
      LAB = 5
      GO TO 1000
C-----------------------------------------------------------------------
C------------------------- ABSCHLUSS HOEHENLINIEN ----------------------
C-----------------------------------------------------------------------
C------  SET POLYLINE COLOUR INDEX
  311 CALL GSPLCI(1)
      INTLIN = INTA
      CALL GRSPTS(INTLIN)
C-----------------------------------------------------------------------
C-------- GRADIENTENFELD -----------------------------------------------
C-----------------------------------------------------------------------
      IF (BIT(6) .OR. BIT(7)) THEN
          IF (BIT(6)) THEN
              I = INDEX(OPTION,'G') + 1
              KIND = 2
          ELSE
              I = INDEX(OPTION,'V') + 1
              KIND = -2
          END IF
          IF (I.LE.LEN(OPTION) .AND. OPTION(I:I).GE.'1' .AND.
     +        OPTION(I:I).LE.'8') THEN
              READ (OPTION(I:I),FMT='(I1)') FARBE
CCC         SET POLYLINE COLOUR INDEX
              CALL GSPLCI(FARBE)
              IOBFA = 0
          ELSE
C ------ UEBERNIMM UND SORTIERE DIE WERTE AUFSTEIGEND
              IOBFA = 1234567890
              NWERT = MIN(MAXHHL,NMAX)
              DO 7007 I = 1,NWERT
                  WER(I) = WERT(I)
C------------- Farbe ohne ICOL ! das wird in GRGFLD gesetzt !
                  IF (INT(I).GT.0) THEN
                      FARBE = INFA(MAX(MOD(INT(I),32),1))
                  ELSE
                      FARBE = INFA(MAX(MOD(ABS(INT(I)),32),1))
                  END IF
                  IFA(I) = FARBE
 7007         CONTINUE
              DO 7009 I = 1,NWERT - 1
                  WEMI = WER(I)
                  K = I
                  DO 7008 J = I + 1,NWERT
                      IF (WER(J).LT.WEMI) THEN
                          WEMI = WER(J)
                          K = J
                      END IF
 7008             CONTINUE
                  JFA = IFA(K)
                  IFA(K) = IFA(I)
                  IFA(I) = JFA
                  WER(K) = WER(I)
                  WER(I) = WEMI
 7009         CONTINUE
          END IF
          INCX = (IMAX-1)/MAXX + 1
          INCY = (JMAX-1)/MAXY + 1
          IF (IPOLAB.EQ.1234567890 .OR.
     +        (IPOLAB.GE.987654321.AND.IPOLAB.LE.987654323)) THEN
              INCX = 1
              INCY = 1
          END IF
          IRX = MOD(IMAX-1,INCX)
          IRY = MOD(JMAX-1,INCY)
          NX = (IMAX-1)/INCX + 1
          NY = (JMAX-1)/INCY + 1
          OO3 = PP(3)
          OO4 = PP(4)
          QQ3 = PP(1) + (IMAX-1-IRX)* (OO3-PP(1))/ (IMAX-1)
          QQ4 = PP(2) + (JMAX-1-IRY)* (OO4-PP(2))/ (JMAX-1)
          PP(3) = QQ3
          PP(4) = QQ4
          XUNITS = (PP(7)-PP(5))/ (PP(3)-PP(1))
          YUNITS = (PP(8)-PP(6))/ (PP(4)-PP(2))
          XLI = PP(5) - PP(1)*XUNITS
          YLI = PP(6) - PP(2)*YUNITS
          XRE = PP(7) + (XMAXCM-PP(3))*XUNITS
          YRE = PP(8) + (YMAXCM-PP(4))*YUNITS
          CALL GSWN(1,XLI,XRE,YLI,YRE)

          IF (IPOLAB.EQ.1234567890 .OR.
     +        (IPOLAB.GE.987654321.AND.IPOLAB.LE.987654323)) THEN
              FF5 = PP(5)
              FF6 = PP(6)
              FF7 = PP(7)
              FF8 = PP(8)
              PP(5) = XX(1)
              PP(6) = YY(1)
              PP(7) = XX(IMAX)
              PP(8) = YY(JMAX)
          END IF
          CALL GRGFLD(N1,TAB,1,INCX,NX,1,INCY,NY,KIND,.4,.15,25.,.05,
     +                IER)
          IF (IPOLAB.EQ.1234567890 .OR.
     +        (IPOLAB.GE.987654321.AND.IPOLAB.LE.987654323)) THEN
              PP(5) = FF5
              PP(6) = FF6
              PP(7) = FF7
              PP(8) = FF8
          END IF

          PP(3) = OO3
          PP(4) = OO4
          XUNITS = (PP(7)-PP(5))/ (PP(3)-PP(1))
          YUNITS = (PP(8)-PP(6))/ (PP(4)-PP(2))
          XLI = PP(5) - PP(1)*XUNITS
          YLI = PP(6) - PP(2)*YUNITS
          XRE = PP(7) + (XMAXCM-PP(3))*XUNITS
          YRE = PP(8) + (YMAXCM-PP(4))*YUNITS
          CALL GSWN(1,XLI,XRE,YLI,YRE)
      END IF
C-----------------------------------------------------------------------
C     WER NICHT OPTIONEN='R','N' ANGIBT, KANN MIT GRAXS WEITERMACHEN
C                                GRSCLV IN USERKOORDINATEN
C-----------------------------------------------------------------------
      IF (.NOT. (BIT(2).OR.BIT(3))) GO TO 600
C-----------------------------------------------------------------------
C------------------------- RAND ----------------------------------------
C-----------------------------------------------------------------------
CCC   SET POLYLINE COLOUR INDEX
      CALL GSPLCI(1)
      CALL GSLN(1)
      IF (BIT(2)) THEN
          IF (IPOLAB.EQ.1234567890) THEN
C---------------------------------- POLARKOORDINATEN -------------------
              U = (XX(IMAX)-XX(1))*PIF/ZWOPI
              IF (U.LT.0.999 .OR. U.GT.1.001) THEN
                  DO 8855 IPO = 1,IMAX,IMAX - 1
                      U = XX(IPO)*PIF
                      CP = COS(U)
                      SP = SIN(U)
                      RAX(1) = YY(1)*CP
                      RAY(1) = YY(1)*SP
                      RAX(2) = YY(JMAX)*CP
                      RAY(2) = YY(JMAX)*SP
                      CALL GPL(2,RAX,RAY)
 8855             CONTINUE
              END IF
              DO 8857 IPO = 1,JMAX,JMAX - 1
                  DO 8856 K = 1,IMAX
                      U = XX(K)*PIF
                      X(K) = YY(IPO)*COS(U)
                      Y(K) = YY(IPO)*SIN(U)
 8856             CONTINUE
                  CALL GRLN(X,Y,IMAX)
 8857         CONTINUE
          ELSE IF (IPOLAB.GE.987654321 .AND. IPOLAB.LE.987654323) THEN
C------------------------- BELIEBIGE     KOORDINATEN -------------------
              DO 8867 IPO = 1,JMAX,JMAX - 1
                  DO 8866 K = 1,IMAX
                      CALL GRHHTR(XX(K),YY(IPO),X(K),Y(K))
 8866             CONTINUE
                  CALL GRLN(X,Y,IMAX)
 8867         CONTINUE
              DO 8877 IPO = 1,IMAX,IMAX - 1
                  DO 8876 K = 1,JMAX
                      CALL GRHHTR(XX(IPO),YY(K),X(K),Y(K))
 8876             CONTINUE
                  CALL GRLN(X,Y,JMAX)
 8877         CONTINUE
          ELSE
              RAX(1) = XA
              RAY(1) = YA
              RAX(2) = XE
              RAY(2) = YA
              RAX(3) = XE
              RAY(3) = YE
              RAX(4) = XA
              RAY(4) = YE
              RAX(5) = XA
              RAY(5) = YA
              CALL GPL(5,RAX,RAY)
          END IF
      END IF
C-----------------------------------------------------------------------
C--------------- NETZSTRICHELUNG AM RAND -------------------------------
C-----------------------------------------------------------------------
      IF (BIT(3)) THEN
          D1 = .0075* (YE-YA)
          IF (IPOLAB.EQ.1234567890) THEN
C--------------------------------- POLARKOORDINATEN -------------------
              DO 7788 IPO = 1,IMAX
                  U = XX(IPO)*PIF
                  CP = COS(U)
                  SP = SIN(U)
                  RAX(1) = YY(JMAX)*CP
                  RAY(1) = YY(JMAX)*SP
                  RAX(2) = 1.015*YY(JMAX)*CP
                  RAY(2) = 1.015*YY(JMAX)*SP
                  CALL GPL(2,RAX,RAY)
 7788         CONTINUE
              U = (XX(IMAX)-XX(1))*PIF/ZWOPI
              IF (U.LT.0.999 .OR. U.GT.1.001) THEN
                  U = XX(1)*PIF
                  CP = COS(U)
                  SP = SIN(U)
                  DO 7789 IPO = 1,JMAX
                      RAX(1) = YY(IPO)*CP
                      RAY(1) = YY(IPO)*SP
                      RAX(2) = RAX(1) + D1*SP
                      RAY(2) = RAY(1) - D1*CP
                      CALL GPL(2,RAX,RAY)
 7789             CONTINUE
              END IF
          ELSE IF (IPOLAB.GE.987654321 .AND. IPOLAB.LE.987654323) THEN
C------------------------ beliebige     KOORDINATEN -------------------
              RMU = -.0075
              DO 7798 JM = 1,JMAX,JMAX - 1
                  U = YY(JM) + (YY(JMAX)-YY(1))*RMU
                  RMU = -RMU
                  DO 7797 IPO = 1,IMAX
                      CALL GRHHTR(XX(IPO),YY(JM),RAX(1),RAY(1))
                      CALL GRHHTR(XX(IPO),U,RAX(2),RAY(2))
                      CALL GPL(2,RAX,RAY)
 7797             CONTINUE
 7798         CONTINUE
              RMU = -.0075
              DO 7800 IM = 1,IMAX,IMAX - 1
                  U = XX(IM) + (XX(IMAX)-XX(1))*RMU
                  RMU = -RMU
                  DO 7799 IPO = 1,JMAX
                      CALL GRHHTR(XX(IM),YY(IPO),RAX(1),RAY(1))
                      CALL GRHHTR(U,YY(IPO),RAX(2),RAY(2))
                      CALL GPL(2,RAX,RAY)
 7799             CONTINUE
 7800         CONTINUE
          ELSE
              D2 = D1 + D1
              YS = YA
              XLN(1) = XS
              XLN(2) = XS
              DO 511 I = 1,IMAX
                  XS = XX(I)
                  XLN(1) = XS
                  YLN(1) = YS
                  V = D1
                  IF ((I-1)/10*10.EQ. (I-1)) V = D2
                  Y2 = YS + V
                  XLN(2) = XS
                  YLN(2) = Y2
                  CALL GPL(2,XLN,YLN)
  511         CONTINUE
              YS = YE
              DO 512 I = 1,IMAX
                  XS = XX(I)
                  XLN(1) = XS
                  YLN(1) = YS
                  V = D1
                  IF ((I-1)/10*10.EQ. (I-1)) V = D2
                  Y2 = YS - V
                  XLN(2) = XS
                  YLN(2) = Y2
                  CALL GPL(2,XLN,YLN)
  512         CONTINUE
              D1 = .0075* (XE-XA)
              D2 = D1 + D1
              XS = XA
              DO 513 J = 1,JMAX
                  YS = YY(J)
                  XLN(1) = XS
                  YLN(1) = YS
                  V = D1
                  IF ((J-1)/10*10.EQ. (J-1)) V = D2
                  X2 = XS + V
                  XLN(2) = X2
                  YLN(2) = YS
                  CALL GPL(2,XLN,YLN)
  513         CONTINUE
              XS = XE
              DO 514 J = 1,JMAX
                  YS = YY(J)
                  XLN(1) = XS
                  YLN(1) = YS
                  V = D1
                  IF ((J-1)/10*10.EQ. (J-1)) V = D2
                  X2 = XS - V
                  XLN(2) = X2
                  YLN(2) = YS
                  CALL GPL(2,XLN,YLN)
  514         CONTINUE
          END IF
      END IF
C---- ZURUECKSETZEN ALTE KOORDINATEN ('N' UND/ODER 'R') ---------------

C     NUR FUER AUFRUF UEBER INTERAKTIVE PANEL
      IF (IGR3.EQ.1) THEN
          CALL GQOPSG(IERR,NA)
          IF (IERR.EQ.0) CALL GCLSG
          NAME = NAME + 1
C        WRITE(93,*) '3 NAME=', NAME
c        CALL GCRSG(NAME)
          IF (IXGKS.NE.1234567890 .AND.
     +        GKSTYP.NE.'GLIGKS') CALL GCRSG(NAME)
      END IF

Cwar  CALL  GRSCLV(PP5,PP6,PP7,PP8)
      PP(5) = PP5
      PP(6) = PP6
      PP(7) = PP7
      PP(8) = PP8
      XUNITS = (PP(7)-PP(5))/ (PP(3)-PP(1))
      YUNITS = (PP(8)-PP(6))/ (PP(4)-PP(2))
      XLI = PP(5) - PP(1)*XUNITS
      YLI = PP(6) - PP(2)*YUNITS
      XRE = PP(7) + (XMAXCM-PP(3))*XUNITS
      YRE = PP(8) + (YMAXCM-PP(4))*YUNITS
      CALL GSWN(1,XLI,XRE,YLI,YRE)
C-----------------------------------------------------------------------
C-------- ZURUECKSETZEN DASH, RETURN
C-----------------------------------------------------------------------
  600 PP(10) = ABS1A
      PP(11) = ABS2A
      PP(12) = ABS3A
      IF (PP(11).LT.0.0035) THEN
CCC   'SOLID' LINE
          LTYPE = 1
      ELSE IF (PP(12).LT.PP(10)) THEN
CCC   'DASHED-DOTTED' LINE
          LTYPE = 4
      ELSE IF (ABS(PP(10)-PP(12)).LT.0.1E-2 .AND. PP(10).LT.PP(11)) THEN
CCC   'DOTTED' LINE
          LTYPE = 3
      ELSE
CCC   'DASHED' LINE
          LTYPE = 2
      END IF
CCC   SET LINETYPE
      CALL GSLN(LTYPE)
      CALL GSPLCI(INCOL)
      CALL GSTXCI(INCOT)
      CALL GSTXAL(IALH,IALV)
      CALL GRJMP(0.,0.)
      GO TO 9999
C-----------------------------------------------------------------------
C     EINE ART SUBROUTINE: VERFOLGEN EINER HOEHENLINIE
C-----------------------------------------------------------------------
 1000 K = 1
      NOLINE = .FALSE.
 1001 K = K + 1
      IF (K.GT.MAL) THEN
          F = MAL
          WRITE (*,FMT=*) ' GRHHNL: a contour line is shortened'
          NOLINE = .TRUE.
      END IF
      LOGI = I .LT. 1 .OR. I .GE. IMAX .OR. J .LT. 1 .OR. J .GE. JMAX
      IF (.NOT.LOGI) LOGI = MRK(I,J) .EQ. '9'
      IF (LOGI) THEN
          IF (K.GT.2) THEN
              IF (IPOLAB.EQ.1234567890) THEN
C----------------------------------------------- POLARKOORDINATEN ------
                  DO 7777 IPO = 1,K - 1
                      U = X(IPO)*PIF
                      X(IPO) = Y(IPO)*COS(U)
                      Y(IPO) = Y(IPO)*SIN(U)
 7777             CONTINUE
              ELSE IF (IPOLAB.GE.987654321 .AND.
     +                 IPOLAB.LE.987654323) THEN
C--------------------------------------- beliebige    KOORDINATEN ------
                  DO 7771 IPO = 1,K - 1
                      CALL GRHHTR(X(IPO),Y(IPO),U,V)
                      X(IPO) = U
                      Y(IPO) = V
 7771             CONTINUE
              END IF
              MIM = MIN(K-1,MAL1)
              CALL GRSCDL(MIM,X,Y,IDASH,KETTE4(:LANG),MAL1,XAU,YAU)
          END IF
          GO TO (17,18,19,20,29) LAB
      ELSE
          IF (NOLINE) THEN
              IF (K.GT.3) THEN
                  IF (IPOLAB.EQ.1234567890) THEN
C----------------------------------------------- POLARKOORDINATEN ------
                      DO 7778 IPO = 1,K - 2
                          U = X(IPO)*PIF
                          X(IPO) = Y(IPO)*COS(U)
                          Y(IPO) = Y(IPO)*SIN(U)
 7778                 CONTINUE
                  ELSE IF (IPOLAB.GE.987654321 .AND.
     +                     IPOLAB.LE.987654323) THEN
C--------------------------------------- beliebige    KOORDINATEN ------
                      DO 7772 IPO = 1,K - 2
                          CALL GRHHTR(X(IPO),Y(IPO),U,V)
                          X(IPO) = U
                          Y(IPO) = V
 7772                 CONTINUE
                  END IF
                  MIM = MIN(K-2,MAL1)
                  CALL GRSCDL(MIM,X,Y,IDASH,KETTE4(:LANG),MAL1,XAU,YAU)
              END IF
              K = 1
          END IF
          NOLINE = .FALSE.
          GO TO (1100,1200,1300,1400,1500,1600,1700,1800,3100,3200,3300,
     +           3400) NRF
      END IF
C----------------------- 1100,1200 : VON UNTEN
 1100 Y1 = Y2
      Y2 = YY(J+1)
      Z11 = Z12
      Z21 = Z22
      Z12 = TAB(I,J+1) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
      IF (Z22.GT.0.) THEN
          X(K) = X2
          Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
          NRF = 8
          IF (Z12.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '0'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          I = I + 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z12.GT.0.) THEN
              NOLINE = TAB(I+1,J+1) .GT. EX1 .AND. TAB(I+1,J+1) .LT. EX2
              Y(K) = Y2
              X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
              J = J + 1
              GO TO 1001
          ELSE
              X(K) = X1
              IF (Z11.EQ.Z12) THEN
                  Y(K) = Y2
                  IF (Z12.EQ.Z22) THEN
                      K = K + 1
                      X(K) = X2
                      Y(K) = Y2
                      I = I + 1
                      NRF = 8
                  ELSE
                      J = J + 1
                      NRF = 1
                  END IF
              ELSE
                  NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
                  Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
                  NRF = 3
                  I = I - 1
              END IF
              GO TO 1001
          END IF
      END IF
 1200 Y1 = Y2
      Y2 = YY(J+1)
      Z11 = Z12
      Z21 = Z22
      Z12 = TAB(I,J+1) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
      IF (Z12.GT.0.) THEN
          X(K) = X1
          Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
          NRF = 4
          IF (Z22.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '1'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          I = I - 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z22.GT.0.) THEN
              Y(K) = Y2
              X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
              NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
              J = J + 1
              GO TO 1001
          ELSE
              X(K) = X2
              IF (Z22.NE.Z21) THEN
                  Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
                  NOLINE = TAB(I+1,J+1) .GT. EX1 .AND.
     +                     TAB(I+1,J+1) .LT. EX2
                  NRF = 7
                  I = I + 1
              ELSE
                  Y(K) = Y2
                  IF (Z12.EQ.Z22) THEN
                      K = K + 1
                      X(K) = X1
                      Y(K) = Y2
                      NRF = 4
                      I = I - 1
                  ELSE
                      J = J + 1
                      NRF = 2
                  END IF
              END IF
              GO TO 1001
          END IF
      END IF
C----------------------- 1300,1400 : VON RECHTS
 1300 X2 = X1
      X1 = XX(I)
      Z21 = Z11
      Z22 = Z12
      Z11 = TAB(I,J) - Z0
      Z12 = TAB(I,J+1) - Z0
      NOLINE = TAB(I+1,J+1) .GT. EX1 .AND. TAB(I+1,J+1) .LT. EX2
      IF (Z12.GT.0.) THEN
          Y(K) = Y2
          X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
          NRF = 1
          IF (Z11.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '2'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          J = J + 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z11.GT.0.) THEN
              X(K) = X1
              Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
              NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
              I = I - 1
              GO TO 1001
          ELSE
              Y(K) = Y1
              IF (Z11.EQ.Z21) THEN
                  X(K) = X1
                  IF (Z11.EQ.Z12) THEN
                      K = K + 1
                      X(K) = X1
                      Y(K) = Y2
                      J = J + 1
                      NRF = 1
                  ELSE
                      I = I - 1
                      NRF = 3
                  END IF
              ELSE
                  X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
                  NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
                  NRF = 6
                  J = J - 1
              END IF
              GO TO 1001
          END IF
      END IF
 1400 X2 = X1
      X1 = XX(I)
      Z21 = Z11
      Z22 = Z12
      Z11 = TAB(I,J) - Z0
      Z12 = TAB(I,J+1) - Z0
      NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
      IF (Z11.GT.0.) THEN
          Y(K) = Y1
          X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
          NRF = 5
          IF (Z12.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '3'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          J = J - 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z12.GT.0.) THEN
              X(K) = X1
              Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
              NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
              I = I - 1
              GO TO 1001
          ELSE
              Y(K) = Y2
              IF (Z22.EQ.Z12) THEN
                  X(K) = X1
                  IF (Z12.EQ.Z11) THEN
                      K = K + 1
                      X(K) = X1
                      Y(K) = Y1
                      NRF = 5
                      J = J - 1
                  ELSE
                      I = I - 1
                      NRF = 4
                  END IF
              ELSE
                  X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
                  NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
                  NRF = 2
                  J = J + 1
              END IF
              GO TO 1001
          END IF
      END IF
C----------------------- 1500,1600 : VON OBEN
 1500 Y2 = Y1
      Y1 = YY(J)
      Z12 = Z11
      Z22 = Z21
      Z11 = TAB(I,J) - Z0
      Z21 = TAB(I+1,J) - Z0
      NOLINE = TAB(I+1,J+1) .GT. EX1 .AND. TAB(I+1,J+1) .LT. EX2
      IF (Z21.GT.0.) THEN
          X(K) = X2
          Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
          NRF = 7
          IF (Z11.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '4'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          I = I + 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z11.GT.0.) THEN
              Y(K) = Y1
              X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
              NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
              J = J - 1
              GO TO 1001
          ELSE
              X(K) = X1
              IF (Z11.EQ.Z12) THEN
                  Y(K) = Y1
                  IF (Z11.EQ.Z21) THEN
                      K = K + 1
                      X(K) = X2
                      Y(K) = Y1
                      NRF = 7
                      I = I + 1
                  ELSE
                      J = J - 1
                      NRF = 5
                  END IF
              ELSE
                  Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
                  NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
                  NRF = 4
                  I = I - 1
              END IF
              GO TO 1001
          END IF
      END IF
 1600 Y2 = Y1
      Y1 = YY(J)
      Z12 = Z11
      Z22 = Z21
      Z11 = TAB(I,J) - Z0
      Z21 = TAB(I+1,J) - Z0
      NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
      IF (Z11.GT.0.) THEN
          X(K) = X1
          Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
          NRF = 3
          IF (Z21.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '5'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          I = I - 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z21.GT.0.) THEN
              Y(K) = Y1
              X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
              NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
              J = J - 1
              GO TO 1001
          ELSE
              X(K) = X2
              IF (Z22.EQ.Z21) THEN
                  Y(K) = Y1
                  IF (Z11.EQ.Z21) THEN
                      K = K + 1
                      X(K) = X1
                      Y(K) = Y1
                      I = I - 1
                      NRF = 3
                  ELSE
                      J = J - 1
                      NRF = 6
                  END IF
              ELSE
                  Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
                  NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
                  I = I + 1
                  NRF = 8
              END IF
              GO TO 1001
          END IF
      END IF
C----------------------- 1700,1800 : VON LINKS
 1700 X1 = X2
      X2 = XX(I+1)
      Z11 = Z21
      Z12 = Z22
      Z21 = TAB(I+1,J) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      NOLINE = TAB(I,J+1) .GT. EX1 .AND. TAB(I,J+1) .LT. EX2
      IF (Z22.GT.0.) THEN
          Y(K) = Y2
          X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
          NRF = 2
          IF (Z21.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '6'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          J = J + 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z21.GT.0.) THEN
              X(K) = X2
              Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
              NOLINE = TAB(I+1,J+1) .GT. EX1 .AND. TAB(I+1,J+1) .LT. EX2
              I = I + 1
              GO TO 1001
          ELSE
              Y(K) = Y1
              IF (Z11.EQ.Z21) THEN
                  X(K) = X2
                  IF (Z21.EQ.Z22) THEN
                      K = K + 1
                      X(K) = X2
                      Y(K) = Y2
                      NRF = 2
                      J = J + 1
                  ELSE
                      NRF = 7
                      I = I + 1
                  END IF
              ELSE
                  X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
                  NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
                  NRF = 5
                  J = J - 1
              END IF
              GO TO 1001
          END IF
      END IF
 1800 X1 = X2
      X2 = XX(I+1)
      Z11 = Z21
      Z12 = Z22
      Z21 = TAB(I+1,J) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      NOLINE = TAB(I,J) .GT. EX1 .AND. TAB(I,J) .LT. EX2
      IF (Z21.GT.0.) THEN
          Y(K) = Y1
          X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
          NRF = 6
          IF (Z22.GT.0.) THEN
              MRK(I,J) = '9'
          ELSE
              IF (MRK(I,J).EQ.'8') THEN
                  MRK(I,J) = '7'
              ELSE
                  MRK(I,J) = '9'
              END IF
          END IF
          J = J - 1
          GO TO 1001
      ELSE
          MRK(I,J) = '9'
          IF (Z22.GT.0.) THEN
              X(K) = X2
              Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
              NOLINE = TAB(I+1,J) .GT. EX1 .AND. TAB(I+1,J) .LT. EX2
              I = I + 1
              GO TO 1001
          ELSE
              Y(K) = Y2
              IF (Z22.EQ.Z12) THEN
                  X(K) = X2
                  IF (Z22.EQ.Z21) THEN
                      K = K + 1
                      X(K) = X2
                      Y(K) = Y1
                      NRF = 6
                      J = J - 1
                  ELSE
                      I = I + 1
                  END IF
              ELSE
                  X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
                  NOLINE = TAB(I+1,J+1) .GT. EX1 .AND.
     +                     TAB(I+1,J+1) .LT. EX2
                  J = J + 1
                  NRF = 1
              END IF
              GO TO 1001
          END IF
      END IF
C------------------------------------ VON UNTEN ------------------------
 3100 Y1 = Y2
      Y2 = YY(J+1)
      Z11 = Z12
      Z21 = Z22
      Z12 = TAB(I,J+1) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      MRK(I,J) = '9'
      IF (Z22.EQ.Z21) THEN
          X(K) = X2
          Y(K) = Y1
          K = K + 1
          X(K) = X2
          Y(K) = Y2
          IF (Z12.GT.0E0) THEN
              J = J + 1
              NRF = 1
          ELSE
              IF (Z12.LT.0E0) THEN
                  J = J + 1
                  NRF = 2
              ELSE
                  K = K + 1
                  X(K) = X1
                  Y(K) = Y2
                  K = K + 1
                  X(K) = X1
                  Y(K) = Y1
                  J = J - 1
                  NRF = 11
              END IF
          END IF
      ELSE
          IF(.NOT.SEL(Z12,Z22) ) THEN
              Y(K) = Y2
              X(K) = (X1*Z22-X2*Z12)/ (Z22-Z12)
              K = K + 1
              Y(K) = Y1
              X(K) = X1
              K = K + 1
          END IF
          X(K) = X2
          Y(K) = Y1
          I = I + 1
          IF (Z22.GT.0E0) THEN
              NRF = 8
          ELSE
              NRF = 7
          END IF
      END IF
      GO TO 1001
C------------------------------------ VON RECHTS -----------------------
 3200 X2 = X1
      X1 = XX(I)
      Z21 = Z11
      Z22 = Z12
      Z11 = TAB(I,J) - Z0
      Z12 = TAB(I,J+1) - Z0
      MRK(I,J) = '9'
      IF (Z12.EQ.Z22) THEN
          X(K) = X2
          Y(K) = Y2
          K = K + 1
          X(K) = X1
          Y(K) = Y2
          IF (Z11.GT.0E0) THEN
              I = I - 1
              NRF = 3
          ELSE
              IF (Z12.LT.0E0) THEN
                  I = I - 1
                  NRF = 2
              ELSE
                  K = K + 1
                  X(K) = X1
                  Y(K) = Y1
                  K = K + 1
                  X(K) = X2
                  Y(K) = Y1
                  I = I + 1
                  NRF = 12
              END IF
          END IF
      ELSE
          IF (.NOT.SEL(Z11,Z12)) THEN
              X(K) = X1
              Y(K) = (Y1*Z12-Y2*Z11)/ (Z12-Z11)
              K = K + 1
              X(K) = X2
              Y(K) = Y1
              K = K + 1
          END IF
          X(K) = X2
          Y(K) = Y2
          J = J + 1
          IF (Z12.GT.0E0) THEN
              NRF = 1
          ELSE
              NRF = 2
          END IF
      END IF
      GO TO 1001
C------------------------------------ VON OBEN -------------------------
 3300 Y2 = Y1
      Y1 = YY(J)
      Z12 = Z11
      Z22 = Z21
      Z11 = TAB(I,J) - Z0
      Z21 = TAB(I+1,J) - Z0
      MRK(I,J) = '9'
      IF (Z11.EQ.Z12) THEN
          X(K) = X1
          Y(K) = Y2
          K = K + 1
          X(K) = X1
          Y(K) = Y1
          IF (Z21.GT.0E0) THEN
              J = J - 1
              NRF = 6
          ELSE
              IF (Z21.LT.0E0) THEN
                  J = J - 1
                  NRF = 5
              ELSE
                  K = K + 1
                  X(K) = X2
                  Y(K) = Y1
                  K = K + 1
                  X(K) = X2
                  Y(K) = Y2
                  J = J + 1
                  NRF = 9
              END IF
          END IF
      ELSE
          IF (.NOT.SEL(Z11,Z21) ) THEN
              Y(K) = Y1
              X(K) = (X1*Z21-X2*Z11)/ (Z21-Z11)
              K = K + 1
              X(K) = X2
              Y(K) = Y2
              K = K + 1
          END IF
          X(K) = X1
          Y(K) = Y2
          I = I + 1
          IF (Z11.GT.0E0) THEN
              NRF = 3
          ELSE
              NRF = 4
          END IF
      END IF
      GO TO 1001
C--------------------------------- VON LINKS --------------------------
 3400 X1 = X2
      X2 = XX(I+1)
      Z11 = Z21
      Z12 = Z22
      Z21 = TAB(I+1,J) - Z0
      Z22 = TAB(I+1,J+1) - Z0
      MRK(I,J) = '9'
      IF (Z11.EQ.Z21) THEN
          X(K) = X1
          Y(K) = Y1
          K = K + 1
          X(K) = X2
          Y(K) = Y1
          IF (Z22.GT.0E0) THEN
              I = I + 1
              NRF = 8
          ELSE
              IF (Z22.LT.0E0) THEN
                  I = I + 1
                  NRF = 7
              ELSE
                  K = K + 1
                  X(K) = X2
                  Y(K) = Y2
                  K = K + 1
                  X(K) = X1
                  Y(K) = Y2
                  I = I - 1
                  NRF = 10
              END IF
          END IF
      ELSE
          IF (.NOT.SEL(Z21,Z22) ) THEN
              X(K) = X2
              Y(K) = (Y1*Z22-Y2*Z21)/ (Z22-Z21)
              K = K + 1
              X(K) = X1
              Y(K) = Y2
              K = K + 1
          END IF
          X(K) = X1
          Y(K) = Y1
          J = J - 1
          IF (Z21.GT.0E0) THEN
              NRF = 6
          ELSE
              NRF = 5
          END IF
      END IF
      GO TO 1001

      contains 
      
      function SEL(UU,VV) result(erg)
      REAL,intent(in) ::UU,VV
      LOGICAL erg
      
      erg =  UU > 0. .AND. VV > 0. .OR. UU < 0. .AND. VV < 0.
      END function SEL
    
 9999 END subroutine grhl
