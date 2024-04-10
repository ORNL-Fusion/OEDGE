C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C----------------------------------------------------------------------
C     UPDATE 15.6.87    BERECHNUNG VON FU NACH NEWTON-COTES     GROTEN
C     UPDATE 04.05.88   POLARKOORDINATEN : IPOLAR, KORREKTUR 28.9.88
C     UPDATE 19.07.91   GROTEN
C     UPDATE 21.11.91   GROTEN; EPS anders berechnet wegen RS 6000
C----------------------------------------------------------------------
      SUBROUTINE GRGFLD (IMAX, F, ISTAX, INCX, NX, ISTAY, INCY, NY,
     >                   KIND, HOEHE, BREITE, WINKEL, DICKE, IER)


C     DAS UNTERPROGRAMM ZEICHNET EIN FELD VON PFEILEN IN EIN RASTER


C***********************************************************************
C     PARAMETER
C
C     IMAX   : INTEGER
C              1. DIMENSION DER MATRIX   F
C     F    : REAL, (IMAX,*)-MATRIX
C              TABELLE DER FUNKTIONSWERTE, DEREN GRADIENTENFELD
C              GEZEICHNET WERDEN SOLL.
C     ISTAX,ISTAY: INTEGER
C              STARTINDEX FUER X BZW. Y ZUM NEHMEN DER WERTE AUS DER
C              MATRIX F      (1. KOORDINATE: X)
C     INCX,INCY: INTEGER
C              INKREMENT ZUM NEHMEN DER WERTE AUS DER MATRIX  F
C     NX,NY  : INTEGER
C              ANZAHL RASTERPUNKTE (WERTE) IN X- BZW. Y-RICHTUNG
C     KIND   : INTEGER  ^=0
C                KIND POSITIV : GRADIENTEN MIT LAENGE UND RICHTUNG.
C                KIND NEGATIV : ALLE PFEILE GLEICH LANG (RICHTUNGEN)
C                KIND=+/-1       PFEILSPITZE RECHTWINKLIG
C                     +/-2       PFEILSPITZEN SPITZ (MIT "WINKEL")
C                     +/-3       PFEILSCHAFT HAT "DICKE", SPITZE WIE 1
C                     +/-4       PFEILSCHAFT HAT "DICKE", SPITZE WIE 2
C                     +/-5       PFEIL IST OHNE SPITZE, NUR STRICH
C                 KIND PUR       PFEILMITTE IM ZENTRUM
C                 KIND +/- 1000  PFEILFUSS IM ZENTRUM
C     IER    : INTEGER
C              FEHLERPARAMETER
C              EINGABE: =1>    NACHRICHT UND RUECKSPRUNG
C                       =2>    RUECKSPRUNG
C                       SONST> NACHRICHT UND STOP
C              AUSGABE: =10>  KIND HAT EINEN NICHT VORGESEHENEN WERT
C                       =20>  ISTAX+(NX-1)*INCX > IMAX
C     HOEHE  : REAL
C              HOEHE JEDER PFEILSPITZE (WIRD RELATIV ZUR PFEILLAENGE
C              GESETZT; IST JEDOCH HOEHE, BREITE ODER DICKE NEGATIV,
C              WIRD DIE RASTERGROESSE ALS MASSSTAB GENOMMEN).
C     BREITE : REAL
C              HALBE BREITE JEDER PFEILSPITZE, BEZUGSGROESSE SIEHE HOEHE
C     WINKEL : REAL
C              WINKEL ZWISCHEN PFEILSCHAFT UND UNTERER ECKE DER PFEIL-
C              SPITZE; ANGABE IN GRAD; VERNUENFTIG: WERTE UNTER 90 GRAD.
C              ANGABE NUR WIRKSAM, WENN KIND= 2, 4
C     DICKE  : HALBE DICKE DES PFEILS, BEZUGSGROESSE WIE BEI HOEHE
C              ANGABE NUR WIRKSAM, WENN KIND= 3, 4
C     COMMON /GRWEFA/ IOBFA,NWERT,WERT(MAXHHL),IFA(MAXHHL)
C            IOBFA    =1234567890, MACHT DIE SCHICHTFAERBUNG AKTIV
C            NWERT    ANZAHL DER PAARE WERT-IFA (FUNKTIONSWERT-FARBE)
C                     NWERT<=16; NICHT MEHR ALS 17 HOEHENSCHICHTEN
C            WERT( )  AUFSTEIGEND SORTIERTE FUNKTIONSWERTE
C            IFA( )   ZUGEHOERIGE FARBEN.
C                     DER FUNKTIONSWERT F(I,J) BEKOMMT EINEN PFEIL IN
C                     DER FARBE IFA(K), WENN F(I,J)<=WERT(K) FUER DAS
C                     KLEINSTE K. WENN EIN FUNKTIONSWERT GROESSER IST
C                     ALS WERT(NWERT), BEKOMMT SEIN PFEIL DIE VOR
C                     AUFRUF VON GRGFLD GESETZTE FARBE.
C     COMMON /GRPOLAR/ IPOLAR       ABS(IPOLAR)=1234567890  POLARKOORDI.
C                                   IPOLAR<0 IN BOGENMASS
C                                   IPOLAR>0 IN GRAD
C     VORHER GRSCLV(PHIMIN,RMIN,PHIMAX,RMAX) AUFRUFEN!
C***********************************************************************
      INTEGER IMAX, ISTAX, INCX, NX, ISTAY, INCY, NY, KIND, IER
      REAL HOEHE, BREITE, WINKEL, DICKE
      REAL F(IMAX,*)

C***********************************************************************
C     PLOTPARAMETER

      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPOLAR/ IPOLAR
CDEC$ PSECT /GRPOLAR/ NOSHR

      SAVE /GRPP/, /GRPOLAR/

C***********************************************************************
C     LOKALE KONSTANTEN
C
C     EPS    : REAL
C              ZUR UEBERPRUEFUNG ZWEIER REAL-VARIABLEN AUF GLEICHHEIT
C     MAXST  : INTEGER, MAX. ANZAHL STUETZPUNKTE FUER EINEN PFEIL
C***********************************************************************

      REAL EPS, EPS1, PIFACT
      PARAMETER(MAXST=7, PIFACT=6.2831853072/360.)

      PARAMETER (GITZA=1./(720.*720.),EX1=-75.7501E20,EX2=-75.7499E20)

      PARAMETER (MAXHHL=100)

C***********************************************************************
C     LOKALE VARIABLEN
C
C     ICX    : INTEGER
C              HILFSVARIABLE; ZAEHLT DIE PFEILE, DIE AUF EINER HOEHE
C              LIEGEN
C     I,J    : INTEGER
C              LAUFINDIZES IN DO-LOOPS
C     DURCHM : REAL
C              MAX. MOEGLICHE LAENGE DER PFEILE
C     LMAX   : REAL
C              MAX. VORHANDENE PFEILLAENGE
C     XWO,YWO: INTEGER
C              KOORDINATEN DES RASTERPUNKTES, AN DEM EIN PFEIL GEZEICH-
C              NET WERDEN SOLL
C     XAB,YAB: INTEGER
C              SCHRITTWEITE IN X-/Y-RICHTUNG, UM ZUM NAECHSTEN RASTER-
C              PUNKT ZU GELANGEN
C     EPSREL : REAL
C              = EPS * LMAX
C     HREL,BREL,DREL: REAL
C              HOEHE UND BREITE DER PFEILSPITZE UND DICKE DES PFEILS IN
C              ABHAENGIGKEIT VON DER RASTERGROESSE
C     VL     : REAL
C              VEKTORLAENGE
C     PP5-PP8: REAL
C              PLOTPARAMETER PP(5)-PP(8)
C     ISTUE  : INTEGER, AKT. ANZAHL DER PUNKTE FUER EINEN PFEIL
C     XX(MAXST),YY(MAXST): REAL, STUETZPUNKTE ZUR ZEICHNUNG EINES PFEILS
C                     IM FOLGENDEN: IN KLAMMERN DIE BEDEUTUNG, FALLS
C                     DER PFEIL DICKER GEZEICHNET WERDEN SOLL
C        XX(1),YY(1) : BASIS DES PFEILS (EINE BASIS)
C        XX(2),YY(2) : MITTELPUNKT DER UNTEREN KANTE DER PFEILSPITZE
C                      (EIN MITTELPUNKT)
C        XX(3),YY(3) : EINE UNTERE ECKE DER PFEILSPITZE
C        XX(4),YY(4) : PFEILSPITZE, ZIELPUNKT
C        XX(5),YY(5) : EINE UNTERE ECKE DER PFEILSPITZE
C        XX(6),YY(6) : ENTSPRICHT XX(2,J),YY(2,J)
C                      (DER ZWEITE MITTELPUNKT)
C        XX(7),YY(7) : ---   (ZWEITE BASIS)
C        PIFACT      : REAL, =PI/180.
C        LENGTH      : REAL, LAENGE EINES PFEILS
C        XM,YM       : REAL, MITTELPUNKT DER UNTEREN PFEILSPITZENKANTE,
C                      WENN DIE BEIDEN UNTEREN 'ECKEN' DURCH EINE GERADE
C                      LINIE VERBUNDEN WERDEN
C        PFEKRO      : REAL, ABSTAND ZWISCHEN DEN BEIDEN PFEILSPITZEN
C        VERSCH      : REAL, HILFSVARIABLE BEI DER BERECHNUNG DER PUNKTE
C                      FUER PFEILE DES TYPS ABS(KIND)=4.
C        ABSK        : REAL, = ABS(KIND)
C        BWINK       : REAL, 'WINKEL' IN BOGENMASS
C***********************************************************************

      INTEGER ICX, I,J,ISTUE, KIN, ICOL(16)
      REAL LENGTH, XM, YM, VERSCH, PFEKRO, ABSK, BWINK
      REAL DURCHM, LMAX, XWO, YWO, XAB, YAB, EPSREL, HREL,BREL, DREL
      REAL VL,  PP5, PP6, PP7, PP8
      REAL XX(MAXST), YY(MAXST)
      REAL BCA,BSA,CA,SA,DSA,DCA,FACLAN,HMASS,TANBWI,UVA,UVB,WINK
      LOGICAL FEWX,FEWY,ZEN,FEST
      COMMON /GRWEFA/ IOBFA,NWERT,WERT(MAXHHL),IFA(MAXHHL)
CDEC$ PSECT /GRWEFA/ NOSHR
      SAVE /GRWEFA/
      REAL NEWCOT(11,10)
      DATA NEWCOT /360.,360.,  0.,  0.,  0.,  0.,  0., 0., 0., 0., 0.,
     >             120.,480.,120.,  0.,  0.,  0.,  0., 0., 0., 0., 0.,
     >              90.,270.,270., 90.,  0.,  0.,  0., 0., 0., 0., 0.,
     >              60.,240.,120.,240., 60.,  0.,  0., 0., 0., 0., 0.,
     >              60.,240.,105.,135.,135., 45.,  0., 0., 0., 0., 0.,
     >              40.,160., 80.,160., 80.,160., 40., 0., 0., 0., 0.,
     >              40.,160., 70., 90., 90., 70.,160.,40., 0., 0., 0.,
     >              30.,120., 60.,120., 60.,120., 60.,120.,30.,0., 0.,
     >              30., 90., 90., 60., 90., 90., 60., 90.,90.,30.,0.,
     >              24., 96., 48., 96., 48., 96., 48., 96.,48.,96.,24./
      DATA FR/0/,BREL/0/,RAD/0/,HREL/0/,PIF/0/,DREL/0/,PFEKRO/0/,IR/0/,
     $     XWO/0/,YWO/0/,DELR/0/,DELPHI/0/,VERSCH/0/,SP/0/,CP/0/,
     $     YMIT/0/,XMIT/0/
      DATA ICOL/1,2,4,3,6,5,7,8,9,10,11,12,13,14,15,16/

C---- EPS : kleinste positive reelle Zahl X fuer die gilt '1.+X > 1.'

C21.11.91  EPS=(2.-(4./3.+2./3.))**2/(4.-(8./3.+4./3.))
      EPS = .00390625
      EPS1 = EPS*.5
    1 IF ( 1.+EPS1 .NE. 1. ) THEN
         EPS = EPS1
         EPS1 = EPS*.5
         GOTO 1
      ENDIF

C     SICHERN VON PP(5) - PP(8) UND NEUE SKALIERUNG DES ZEICHENFELDES

      PP5 = PP(5)
      PP6 = PP(6)
      PP7 = PP(7)
      PP8 = PP(8)
      CALL GRSCLV (PP(1), PP(2), PP(3), PP(4))
      ZEN=ABS(KIND).GT.1000
      KIN=SIGN(MOD(ABS(KIND),1000),KIND)
      ABSK = ABS(KIN)

C     FEHLERABFRAGEN:

      IF (KIN.EQ.0 .OR. ABSK.GT.5 )THEN

C     A) KIND KORREKT?

        IF (IER .NE. 2)
     >    WRITE(*,*) '*** GRVFLD: THE VALUE OF ''KIND'' IS INCORRECT'
        IF ((IER.NE.1) .AND. (IER.NE.2))  STOP 'GRVFLD 10'
        IER = 10

      ELSEIF (ISTAX + (NX-1)*INCX .GT. IMAX) THEN

C       B) WIRD DIE 1. DIMENSION DER MATRIX F UEBERSCHRITTEN?

        IF (IER .NE. 2)
     >    WRITE(*,*) '*** GRVFLD: (ISTAX + (NX-1)*INCX) > NDIM1'
        IF ((IER.NE.1) .AND. (IER.NE.2))  STOP 'GRVFLD  20'
        IER = 20

      ELSE

C       C) PARAMETER SIND IN ORDNUNG

        IER = 0

C       MAX. EINHEITSLAENGE BESTIMMEN

        IF (ABS(IPOLAR).NE.1234567890) THEN
           XAB = (PP(3)-PP(1)) / (NX-1)
           YAB = (PP(4)-PP(2)) / (NY-1)
        ELSE
           IF (IPOLAR.LT.0.) THEN
              PIF=1.
           ELSE
              PIF=PIFACT
           ENDIF
           XMI=PP6*COS(PIF*PP5)
           XMA=XMI
           YMI=PP6*SIN(PIF*PP5)
           YMA=YMI
           RAD=PP6
           DELPHI=(PP7-PP5)/(NX-1)
           DO 76 J=1,2
              DO 75 I=0,NX-1
                 U=(PP5+DELPHI*I)*PIF
                 X=RAD*COS(U)
                 Y=RAD*SIN(U)
                 XMI=MIN(XMI,X)
                 XMA=MAX(XMA,X)
                 YMI=MIN(YMI,Y)
                 YMA=MAX(YMA,Y)
 75           CONTINUE
              RAD=PP8
 76        CONTINUE
           IF ((YMA-YMI)*(PP(3)-PP(1)).LT. (XMA-XMI)*(PP(4)-PP(2))) THEN
              FR=(PP(3)-PP(1))/(XMA-XMI)
              XMIT=PP(3)-FR*XMA
              YMIT=(PP(4)+PP(2)-FR*(YMA+YMI))/2.
           ELSE
              FR=(PP(4)-PP(2))/(YMA-YMI)
              YMIT=PP(4)-FR*YMA
              XMIT=(PP(3)+PP(1)-FR*(XMA+XMI))/2.
           ENDIF
           DELR=(PP8-PP6)/(NY-1)
           XAB=DELR*FR
           YAB=DELPHI*PIF*(PP6+PP8)/2*FR
           IR=0
        ENDIF
        IF (ZEN) THEN
          DURCHM = MIN(XAB,YAB)*.5
        ELSE
          DURCHM = MIN(XAB,YAB)*.85
        ENDIF

C       MASSE DER PFEILSPITZE UND DICKE DES PFEILS IN ABHAENGIGKEIT VON
C       DER RASTERGROESSE

        WINK=WINKEL
        IF (ABS(MOD(WINKEL-90.,180.)).LE.EPS) WINK=WINKEL-1.
        IF (ABS(WINK).LT.EPS) WINK=1.
        BWINK = WINK*PIFACT
        TANBWI=TAN(BWINK)
        FEST=HOEHE.LT.0. .OR. BREITE.LT.0. .OR. DICKE.LT.0.
        IF (FEST) THEN
           IF (ABS(HOEHE) .GT. 1) THEN
             HREL = DURCHM
           ELSE
             HREL = DURCHM*ABS(HOEHE)
           ENDIF

           IF (ABS(BREITE) .GT. 1) THEN
             BREL = DURCHM
           ELSE
             BREL = DURCHM*ABS(BREITE)
           ENDIF

           IF (ABS(DICKE).GT.1 .OR. ABS(DICKE).GT.ABS(BREITE)) THEN
             DREL = BREL
           ELSE
             DREL = DURCHM*ABS(DICKE)
           ENDIF

           PFEKRO = HREL - BREL/ABS(TANBWI)
           IF (PFEKRO .LT. 0.) THEN
             VERSCH = DREL/BREL*HREL
             PFEKRO = 0.
           ELSE
             VERSCH = PFEKRO + DREL/TANBWI
           ENDIF

        ENDIF

        IF ( ABSK.GE.3 ) THEN
          ISTUE = 7
        ELSE
          ISTUE = 6
        ENDIF


C       MAX. LAENGE BESTIMMEN

        FEWY=INCY.LE.10
        FEWX=INCX.LE.10
        IF (KIN.LT.0) THEN
           LMAX=1.
        ELSE
          LMAX = 0.
          DO 20, J= ISTAY, ISTAY+(NY-2)*INCY, INCY
            IF (ABS(IPOLAR).EQ.1234567890) THEN
               RAD=FR*(PP6+DELR*(.5+IR))
            ENDIF
            DO 10, I= ISTAX, ISTAX+(NX-2)*INCX, INCX
              UVA=0.
              DO 5 K=J,J+INCY,1
                FU1=F(I+INCX,K)
                FU2=F(I,K)
                IF (FU1.GT.EX1 .AND. FU1.LT.EX2) GOTO 10
                IF (FU2.GT.EX1 .AND. FU2.LT.EX2) GOTO 10
                IF (FEWY) THEN
                   UVA=UVA+NEWCOT(K-J+1,INCY)*(FU1-FU2)
                ELSE
                   TRA=720./INCY
                   IF (K.EQ.J .OR. K.EQ.J+INCY) TRA=360./INCY
                   UVA=UVA+TRA*(FU1-FU2)
                ENDIF
    5         CONTINUE
              UVB=0.
              DO 6 K=I,I+INCX,1
                FU1=F(K,J+INCY)
                FU2=F(K,J)
                IF (FU1.GT.EX1 .AND. FU1.LT.EX2) GOTO 10
                IF (FU2.GT.EX1 .AND. FU2.LT.EX2) GOTO 10
                IF (FEWX) THEN
                  UVB=UVB+NEWCOT(K-I+1,INCX)*(FU1-FU2)
                ELSE
                  TRA=720./INCX
                  IF (K.EQ.I .OR. K.EQ.I+INCX) TRA=360./INCX
                  UVB=UVB+TRA*(FU1-FU2)
                ENDIF
    6         CONTINUE
              IF (ABS(IPOLAR).EQ.1234567890) THEN
                 UVA=UVA/RAD*DELR
                 UVB=UVB*DELPHI*PIF/FR
              ENDIF
              VL=SQRT(UVA**2+UVB**2)
              LMAX = MAX(LMAX,VL)
   10       CONTINUE
          IR=IR+1
   20     CONTINUE
        ENDIF
         EPSREL = EPS*LMAX
         FACLAN = DURCHM/LMAX
C        BASIS BERECHNEN UND ENTWEDER LAENGE 'NORMALISIEREN' ODER ZIEL-
C        PUNKT BERECHNEN

         IF (ABS(IPOLAR).NE.1234567890) THEN
            XWO = PP(1) - XAB*.5
            YWO = PP(2) + YAB*.5
         ELSE
            U=(PP5+DELPHI*.5)*PIF
            RAD = FR*(PP6+DELR/2)
            IR=0
         ENDIF
         HMASS = FACLAN*.5
         CALL GQPLCI(J,NORFAR)
         NWE=MIN(NWERT,MAXHHL)

         IR=0
         DO 40, J= ISTAY, ISTAY+(NY-2)*INCY, INCY
            ICX=0
            DO 30, I= ISTAX, ISTAX+(NX-2)*INCX, INCX
              ICX=ICX+1

              IF (ABS(IPOLAR).EQ.1234567890) THEN
                 PHI=(PP5+DELPHI*(ICX-.5))*PIF
                 CP=COS(PHI)
                 SP=SIN(PHI)
                 XWO=XMIT+RAD*CP
                 YWO=YMIT+RAD*SP
              ELSE
                 XWO=XWO+XAB
              ENDIF

              UVA=0.
              DO 35 K=J,J+INCY,1
                FU1=F(I+INCX,K)
                FU2=F(I,K)
                IF (FU1.GT.EX1 .AND. FU1.LT.EX2) GOTO 30
                IF (FU2.GT.EX1 .AND. FU2.LT.EX2) GOTO 30
                IF (FEWY) THEN
                   UVA=UVA+NEWCOT(K-J+1,INCY)*(FU1-FU2)
                ELSE
                   TRA=720./INCY
                   IF (K.EQ.J .OR. K.EQ.J+INCY) TRA=360./INCY
                   UVA=UVA+TRA*(FU1-FU2)
                ENDIF
   35         CONTINUE
              UVB=0.
              DO 36 K=I,I+INCX,1
                FU1=F(K,J+INCY)
                FU2=F(K,J)
                IF (FU1.GT.EX1 .AND. FU1.LT.EX2) GOTO 30
                IF (FU2.GT.EX1 .AND. FU2.LT.EX2) GOTO 30
                IF (FEWX) THEN
                  UVB=UVB+NEWCOT(K-I+1,INCX)*(FU1-FU2)
                ELSE
                  TRA=720./INCX
                  IF (K.EQ.I .OR. K.EQ.I+INCX) TRA=360./INCX
                  UVB=UVB+TRA*(FU1-FU2)
                ENDIF
   36         CONTINUE
              IF (ABS(IPOLAR).EQ.1234567890) THEN
                 U=UVA/RAD*DELR
                 UVB=UVB*DELPHI*PIF/FR
                 UVA=UVB*CP-U*SP
                 UVB=UVB*SP+U*CP
              ENDIF
              IF (IOBFA.EQ.1234567890) THEN
                FU=0.
                DO 38 K=J,J+INCY,1
                  FUX=0.
                  DO 37 L=I,I+INCX,1
                    FU1=F(L,K)
                    IF (FU1.GT.EX1 .AND. FU1.LT.EX2) GOTO 30
                    IF (FEWX) THEN
                       FUX=FUX+NEWCOT(L-I+1,INCX)*FU1
                    ELSE
                       TRA=720./INCX
                       IF (L.EQ.I .OR. L.EQ.I+INCX) TRA=360./INCX
                       FUX=FUX+TRA*FU1
                    ENDIF
   37             CONTINUE
                  IF (FEWY) THEN
                     FU=FU+NEWCOT(K-J+1,INCY)*FUX
                  ELSE
                     TRA=720./INCY
                     IF (K.EQ.J .OR. K.EQ.J+INCY) TRA=360./INCY
                     FU=FU+TRA*FUX
                  ENDIF
   38           CONTINUE
                FU=FU*GITZA
                DO 39 K=1,NWE
                  IF (FU.LE.WERT(K)) GOTO 391
   39           CONTINUE
                KFA=NORFAR
                GOTO 392
  391           KFA=ICOL(IFA(K))
  392           CALL GSPLCI(KFA)
              ENDIF
              VL=SQRT(UVA**2+UVB**2)
              IF (KIN.LT.0) THEN
                IF (VL.NE.0.) THEN
                  UVA=UVA/VL
                  UVB=UVB/VL
                  VL=1.
                ELSE
                  UVA=0.
                  UVB=0.
                  VL=0.
                ENDIF
              ENDIF
              IF (ZEN) THEN
                XX(1) = XWO
                YY(1) = YWO
                XX(4) = XWO + UVA*FACLAN
                YY(4) = YWO + UVB*FACLAN
              ELSE
                XX(1) = XWO - UVA*HMASS
                YY(1) = YWO - UVB*HMASS
                XX(4) = XWO + UVA*HMASS
                YY(4) = YWO + UVB*HMASS
              ENDIF

              LENGTH = FACLAN*VL
              IF (LENGTH.GT.EPSREL) THEN
                 CA = (XX(4)-XX(1))/LENGTH
                 SA = (YY(4)-YY(1))/LENGTH
              ELSE
                 CA=1.
                 SA=0.
              ENDIF

C             EINE PFEILSPITZE WIRD NUR GEZEICHNET, WENN DIESE NICHT
C             LAENGER ALS DER ZUGEHOERIGE PFEIL IST

              IF ((FEST .AND. LENGTH .LT. HREL) .OR. ABSK.EQ.5) THEN

                XX(2) = XX(4)
                YY(2) = YY(4)
                CALL GPL (2,XX,YY)

              ELSE

C               BERECHNUNG DES PUNKTES, AN DEM DIE PFEILSPITZE ANLIEGT

                IF (.NOT.FEST) THEN
                   IF (ABS(HOEHE) .GT. 1) THEN
                     HREL = LENGTH
                   ELSE
                     HREL = LENGTH*ABS(HOEHE)
                   ENDIF
                   IF (ABS(BREITE) .GT. 1) THEN
                     BREL = LENGTH
                   ELSE
                     BREL = LENGTH*ABS(BREITE)
                   ENDIF

                   IF(ABS(DICKE).GT.1.OR.ABS(DICKE).GT.ABS(BREITE)) THEN
                     DREL = BREL
                   ELSE
                     DREL = LENGTH*ABS(DICKE)
                   ENDIF

                   PFEKRO = HREL - BREL/ABS(TANBWI)
                   IF (PFEKRO .LT. 0.) THEN
                     VERSCH = DREL/BREL*HREL
                     PFEKRO = 0.
                   ELSE
                     VERSCH = PFEKRO + DREL/TANBWI
                   ENDIF
                ENDIF
                XX(2) = XX(1) + (LENGTH-HREL) * CA
                YY(2) = YY(1) + (LENGTH-HREL) * SA
                XM = XX(2)
                YM = YY(2)

                IF ((ABSK .EQ. 2) .OR. (ABSK .EQ. 4)) THEN
                   XX(2) = XX(1) + (LENGTH-PFEKRO) * CA
                   YY(2) = YY(1) + (LENGTH-PFEKRO) * SA
                ENDIF

C               BERECHNUNG DER BEIDEN UNTEREN ECKEN DER PFEILSPITZE

                BSA = BREL*SA
                BCA = BREL*CA
                XX(3) = XM-BSA
                YY(3) = YM+BCA
                XX(5) = XM+BSA
                YY(5) = YM-BCA
                XX(6) = XX(2)
                YY(6) = YY(2)

C               BERECHNUNG DER 2 MITTELPUNKTE UND DER 2 BASISPUNKTE,
C               FALLS DER PFEIL DICKER GEZEICHNET WERDEN SOLL

                IF ((ABSK .EQ. 3) .OR. (ABSK .EQ. 4)) THEN
                  IF (ABSK .EQ. 3) THEN
                    XM = XX(2)
                    YM = YY(2)
                  ELSE
                    XM = XX(1) + (LENGTH-VERSCH) * CA
                    YM = YY(1) + (LENGTH-VERSCH) * SA
                  ENDIF
                  DSA = DREL*SA
                  DCA = DREL*CA
                  XX(2) = XM-DSA
                  YY(2) = YM+DCA
                  XX(6) = XM+DSA
                  YY(6) = YM-DCA
                  XM = XX(1)
                  YM = YY(1)
                  XX(1) = XM-DSA
                  YY(1) = YM+DCA
                  XX(7) = XM+DSA
                  YY(7) = YM-DCA

                ENDIF

C               AUFRUF DER GKS-ROUTINE    **********

                CALL GPL (ISTUE,XX,YY)

              ENDIF
   30      CONTINUE
           IF (ABS(IPOLAR).NE.1234567890) THEN
              XWO = PP(1) - XAB*.5
              YWO = YWO + YAB
           ELSE
              IR=IR+1
              RAD=(PP6+DELR*(.5+IR))*FR
           ENDIF
   40    CONTINUE
      ENDIF

C     ZURUECKSETZEN DER SKALIERUNG

      CALL GRSCLV (PP5,PP6,PP7,PP8)

      END
