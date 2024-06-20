C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE:  8.  5. 1991 GROTEN
C UPDATE: 21. 11. 1991 GROTEN; EPS neu berechnet wegen RS 6000
C UPDATE:  8.  1. 1992 GROTEN; NORFAR
C UPDATE: 21. 12. 1992 GROTEN; GRNWPN aufrufen statt GSPLCI
C UPDATE: 22. 12. 1992 GROTEN; KIND=6 und 7 wie 3 und 4 aber gefuellt.
      SUBROUTINE GRVFLD (NDIM1, VA,VB,IFA, ISTAX,INCX,NX, ISTAY,INCY,NY,
     >                   KIND, HOEHE,BREITE,WINKEL, DICKE, IER)


C     DAS UNTERPROGRAMM ZEICHNET EIN FELD VON PFEILEN IN EIN RASTER


C***********************************************************************
C     PARAMETER
C
C     NDIM1  : INTEGER
C              1. DIMENSION DER MATRIZEN VA, VB
C     VA     : REAL, (NDIM1,*)-MATRIX
C              VA(I,J): RADIUS ODER ABWEICHUNG IN X-RICHTUNG DES PFEILS
C                       AM RASTERPUNKT (I,J)
C     VB     : REAL, (NDIM1,*)-MATRIX
C              VB(I,J): WINKEL ODER ABWEICHUNG IN Y-RICHTUNG DES PFEILS
C                       AM RASTERPUNKT (I,J)
C                       WINKELANGABE BEI KIND= -15 BIS -11 IN BOGENMASS
C                                    BEI KIND=  -5 BIS  -1 IN GRAD
C     IFA    : INTEGER, (NDIM1,*)-MATRIX
C              MATRIX MIT FARBWERTEN FUER DIE PFEILE.
C              WENN ALLE PFEILE GLEICH SEIN SOLLEN (FARBE VORHER MIT
C              GRNWPN SETZEN!), KANN MAN 0 ANGEBEN.
C     ISTAX,ISTAY: INTEGER
C              STARTINDEX FUER X BZW. Y ZUM NEHMEN DER WERTE AUS DEN
C              MATRIZEN VA, VB (1. KOORDINATE: X)
C     INCX,INCY: INTEGER
C              INKREMENT ZUM NEHMEN DER WERTE AUS DEN MATRIZEN VA, VB
C     NX,NY  : INTEGER
C              ANZAHL RASTERPUNKTE (WERTE) IN X- BZW. Y-RICHTUNG
C     KIND   : INTEGER
C              <0:   PFEILE ANGEGEBEN MIT RADIUS UND WINKEL (GRAD)
C              <-10:   "        "      "    "     "  WINKEL (BOGENMASS)
C              >0: PFEILE ANGEGEBEN MIT ABWEICHUNG IN X- UND Y-RICHTUNG
C                KIND=1|-1|-11: PFEILSPITZE RECHTWINKLIG
C                     2,-2,-12: PFEILSPITZEN SPITZ (MIT "WINKEL")
C                     3,-3,-13: PFEILSCHAFT HAT "DICKE", SPITZE WIE 1
C                     4,-4,-14: PFEILSCHAFT HAT "DICKE", SPITZE WIE 2
C                     5,-5,-15: PFEIL IST OHNE SPITZE, NUR STRICH
C                KIND  PUR    : PFEILMITTE IM RASTERRECHTECKSMITTELPUNKT
C                KIND +/- 1000: PFEILFUSSPUNKT IM """""""""""""""""""
C                KIND +/- 2000: "" "". PFEILLAENGE NICHT NORMIERT IN CM.
C                KIND +/- 3000: "" "". PFEILLAENGE NICHT NORMIERT USER.
C                     (3000 UND KIND<0 NUR FUER GLEICHEN MASSTAB X-Y)
C     HOEHE  : REAL
C              HOEHE JEDER PFEILSPITZE (WIRD NORMALERWEISE RELATIV
C              ZUR PFEILLAENGE GERECHNET; WENN ABER HOEHE, BREITE ODER
C              DICKE NEGATIV IST, WIRD ES AUF DIE RASTERGROESSE
C              BEZOGEN).
C     BREITE : REAL
C              HALBE BREITE JEDER PFEILSPITZE. ZUM MASSSTAB: SIEHE HOEHE
C     WINKEL : REAL
C              WINKEL ZWISCHEN PFEILSCHAFT UND UNTERER ECKE DER PFEIL-
C              SPITZE; ANGABE IN GRAD; VERNUENFTIG: WERTE UNTER 90 GRAD.
C              ANGABE NUR WIRKSAM, WENN KIND= 2, 4, -2, -4, -12, -14
C                                            +/- 1000  +/- 2000 +/- 3000
C     DICKE  : HALBE DICKE DES PFEILS. ZUM MASSSTAB SIEHE HOEHE.
C              ANGABE NUR WIRKSAM, WENN KIND= 3, 4, -3, -4, -13, -14
C                                            +/- 1000  +/- 2000 +/- 3000
C     IER    : INTEGER
C              FEHLERPARAMETER
C              EINGABE: =1>    NACHRICHT UND RUECKSPRUNG
C                       =2>    RUECKSPRUNG
C                       SONST> NACHRICHT UND STOP
C              AUSGABE: =10>  KIND HAT EINEN NICHT VORGESEHENEN WERT
C                       =20>  ISTAX+(NX-1)*INCX > NDIM1
C***********************************************************************

      INTEGER NDIM1, ISTAX, INCX, NX, ISTAY, INCY, NY, KIND, IER
      REAL HOEHE, BREITE, WINKEL, DICKE
      REAL VA(NDIM1,*), VB(NDIM1,*)
      INTEGER IFA(NDIM1,*)

C***********************************************************************
C     PLOTPARAMETER

      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/


C***********************************************************************
C     LOKALE KONSTANTEN
C
C     EPS    : REAL
C              ZUR UEBERPRUEFUNG ZWEIER REAL-VARIABLEN AUF GLEICHHEIT
C     PI     : REAL
C              4*ATAN(1)
C     MAXST  : INTEGER, MAX. ANZAHL STUETZPUNKTE FUER EINEN PFEIL
C***********************************************************************

      REAL EPS, PI, EPS1
      PARAMETER( PI=3.1415926536,MAXST=7 )

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
C     XWO,YWO: REAL
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
C     HILF   : REAL
C              HILFSVARIABLE
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
C                      FUER PFEILE DES TYPS ABS(KIND)=4 UND KIND=-14
C        ABSK        : INTEGER, = MOD(ABS(KIND),10)
C        BWINK       : REAL, 'WINKEL' IN BOGENMASS
C***********************************************************************

      INTEGER ICX, I, J, ISTUE, KIN, ABSK
      REAL LENGTH, XM, YM, VERSCH, PFEKRO, BWINK
      REAL DURCHM, LMAX, XWO, YWO, XAB, YAB, EPSREL, HREL,BREL, DREL
      REAL HILF,  PP5, PP6, PP7, PP8, PIFACT
      REAL XX(MAXST), YY(MAXST)
      REAL BCA,BSA,CA,SA,DSA,DCA,FACLA1,FACLA2,TANBWI,UVA,UVB,WINK
      LOGICAL ZEN,FEST,NONORM
      COMMON /GRPOLAR/ IPOLAR
CDEC$ PSECT /GRPOLAR/ NOSHR
      SAVE /GRPOLAR/
      DATA RAD/0./, BREL/0./, FR/0./, HREL/0./, DREL/0./, PFEKRO/0./,
     $     DELPHI/0./, VERSCH/0./, DELR/0./, YMIT/0./, XMIT/0./, IR/0/

C---- EPS  : kleinste positive reelle Zahl X fuer die gilt '1.+X > 1.'
C21.11.91 EPS = (2.-(4./3.+2./3.))**2/(4.-(8./3.+4./3.))
      EPS =.00390625
      EPS1=EPS*.5
    1 IF (1.+EPS1 .NE. 1.) THEN
         EPS=EPS1
         EPS1=EPS*.5
         GOTO 1
      ENDIF

C     SICHERN VON PP(5) - PP(8) UND NEUE SKALIERUNG DES ZEICHENFELDES

      PP5 = PP(5)
      PP6 = PP(6)
      PP7 = PP(7)
      PP8 = PP(8)
      CALL GRSCLV (PP(1), PP(2), PP(3), PP(4))
      ZEN=ABS(KIND).GT.1000
      NONORM=ABS(KIND).GT.2000
      KIN=SIGN(MOD(ABS(KIND),1000),KIND)
      ABSK = MOD(ABS(KIN),10)
      FEST=HOEHE.LT.0. OR. BREITE.LT.0. .OR. DICKE.LT.0.

C     FEHLERABFRAGEN:

      IF (KIN.EQ.0 .OR. KIN.GT.7 .OR. ABSK.GT.7 .OR. KIN.LT.-17)THEN

C     A) KIND KORREKT?

        IF (IER .NE. 2)
     >    WRITE(*,*) '*** GRVFLD: THE VALUE OF ''KIND'' IS INCORRECT'
        IF ((IER.NE.1) .AND. (IER.NE.2))  STOP 'GRVFLD 10'
        IER = 10

      ELSEIF (ISTAX + (NX-1)*INCX .GT. NDIM1) THEN

C       B) WIRD DIE 1. DIMENSION DER MATRIZEN VA, VB UEBERSCHRITTEN?

        IF (IER .NE. 2)
     >    WRITE(*,*) '*** GRVFLD: (ISTAX + (NX-1)*INCX) > NDIM1'
        IF ((IER.NE.1) .AND. (IER.NE.2))  STOP 'GRVFLD  20'
        IER = 20

      ELSE
C---------------------------------- ARGUMENTE SIND OHNE FEHLER ---------
        IER = 0

C       MAX. EINHEITSLAENGE BESTIMMEN

        IF (ABS(IPOLAR).NE.1234567890) THEN
           XAB = (PP(3)-PP(1)) / NX
           YAB = (PP(4)-PP(2)) / NY
           pif=1.
        ELSE
           IF (IPOLAR.LT.0.) THEN
              PIF=1.
           ELSE
              PIF=PI/180.
           ENDIF
           XMI=PP6*COS(PIF*PP5)
           XMA=XMI
           YMI=PP6*SIN(PIF*PP5)
           YMA=YMI
           RAD=PP6
           DELPHI=(PP7-PP5)/NX
           DO 76 J=1,2
              DO 75 I=0,NX
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
           DELR=(PP8-PP6)/NY
           RAD = FR*(PP6+DELR/2)
           XAB=DELR
           YAB=DELPHI*PIF*(PP6+PP8)/2
           IR=0
        ENDIF
        IF (ZEN) THEN
          DURCHM=MIN(XAB,YAB)*.5
        ELSE
          DURCHM = MIN (XAB,YAB)*.85
        ENDIF

C       MASSE DER PFEILSPITZE UND DICKE DES PFEILS IN ABHAENGIGKEIT VON
C       DER RASTERGROESSE

        WINK=WINKEL
        IF (ABS(MOD(WINKEL-90.,180.)).LE.EPS) WINK=WINKEL-1.
        IF (ABS(WINK).LT.EPS) WINK=1.
        PIFACT = PI / 180.
        BWINK = WINK*PIFACT
        TANBWI=TAN(BWINK)
        IF (KIN.LT.-10) PIFACT=1.
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

        LMAX = 0.
        DO 20, J= ISTAY, ISTAY+(NY-1)*INCY, INCY
          DO 10, I= ISTAX, ISTAX+(NX-1)*INCX, INCX

            IF (KIN .LT. 0) THEN
              HILF = ABS(VA(I,J))
            ELSE
              HILF = SQRT (VA(I,J)**2 + VB(I,J)**2)
            ENDIF

            LMAX = MAX(LMAX,HILF)

   10     CONTINUE
   20   CONTINUE

        EPSREL = EPS*LMAX
        IF (NONORM) THEN
           IF (ABS(KIND).GT.3000) THEN
              IF (ABS(IPOLAR).NE.1234567890) THEN
                 FACLA1=(PP(3)-PP(1))/(PP7-PP5)
                 FACLA2=(PP(4)-PP(2))/(PP8-PP6)
              ELSE
                 FACLA1=FR
                 FACLA2=FR
              ENDIF
           ELSE
              FACLA1=1.
              FACLA2=1.
           ENDIF
        ELSE
           FACLA1 = DURCHM/LMAX
           FACLA2 = FACLA1
        ENDIF

C       BASIS BERECHNEN UND ENTWEDER LAENGE 'NORMALISIEREN' ODER ZIEL-
C       PUNKT BERECHNEN

         IF (ABS(IPOLAR).NE.1234567890) THEN
            XWO = PP(1) + XAB*.5
            YWO = PP(2) + YAB*.5
         ELSE
            U=(PP5+DELPHI*.5)*PIF
            RAD = FR*(PP6+DELR/2)
            XWO = XMIT+RAD*COS(U)
            YWO = YMIT+RAD*SIN(U)
         ENDIF
         ICX = 1
         HMASS1 = FACLA1*.5
         HMASS2 = FACLA2*.5

         DO 40, J= ISTAY, ISTAY+(NY-1)*INCY, INCY
            DO 30, I= ISTAX, ISTAX+(NX-1)*INCX, INCX

               IF (IFA(1,1).GT.0) THEN
                  CALL GRNWPN(MAX(IFA(I,J),1))
               ENDIF
               UVA=VA(I,J)
               UVB=VB(I,J)
               IF (KIN .LT. 0) THEN

C       A)        BASIS UND 'NORMALISIERTE' LAENGE

                  UVA = FACLA1*ABS(UVA)
                  UVB = UVB * PIFACT
                  CA = COS(UVB)
                  SA = SIN(UVB)
                  IF (ZEN) THEN
                    XX(1)=XWO
                    YY(1)=YWO
                  ELSE
                    XX(1) = XWO-UVA*.5 * CA
                    YY(1) = YWO-UVA*.5 * SA
                  ENDIF
                  LENGTH = UVA
                  XX(4) = XX(1) + UVA*CA
                  YY(4) = YY(1) + UVA*SA
               ELSE

C       B)        BASIS UND ZIELPUNKT

                  IF (ZEN) THEN
                     XX(1) = XWO
                     YY(1) = YWO
                     XX(4) = XWO + UVA*FACLA1
                     YY(4) = YWO + UVB*FACLA2
                  ELSE
                     XX(1) = XWO - UVA*HMASS1
                     YY(1) = YWO - UVB*HMASS2
                     XX(4) = XWO + UVA*HMASS1
                     YY(4) = YWO + UVB*HMASS2
                  ENDIF
                  LENGTH = SQRT( (YY(1)-YY(4))**2 + (XX(1)-XX(4))**2 )
                  IF (LENGTH.GT.EPSREL) THEN
                     CA = (XX(4)-XX(1))/LENGTH
                     SA = (YY(4)-YY(1))/LENGTH
                  ELSE
                     CA=1.
                     SA=0.
                  ENDIF

               ENDIF

C              EINE PFEILSPITZE WIRD NUR GEZEICHNET, WENN DIESE NICHT
C              LAENGER ALS DER ZUGEHOERIGE PFEIL IST

               IF ((FEST .AND. LENGTH.LT.HREL) .OR. ABSK.EQ.5) THEN

                  XX(2) = XX(4)
                  YY(2) = YY(4)
                  CALL GPL (2,XX,YY)

               ELSE

C              BERECHNUNG DES PUNKTES, AN DEM DIE PFEILSPITZE ANLIEGT

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

                     IF (ABS(DICKE).GT.1.OR.ABS(DICKE).GT.ABS(BREITE))
     >               THEN
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

                  IF(ABSK.EQ.2.OR.ABSK.EQ.4.OR.ABSK.EQ.7) THEN
                     XX(2) = XX(1) + (LENGTH-PFEKRO) * CA
                     YY(2) = YY(1) + (LENGTH-PFEKRO) * SA
                  ENDIF

C              BERECHNUNG DER BEIDEN UNTEREN ECKEN DER PFEILSPITZE

                  BSA = BREL*SA
                  BCA = BREL*CA
                  XX(3) = XM-BSA
                  YY(3) = YM+BCA
                  XX(5) = XM+BSA
                  YY(5) = YM-BCA
                  XX(6) = XX(2)
                  YY(6) = YY(2)

C              BERECHNUNG DER 2 MITTELPUNKTE UND DER 2 BASISPUNKTE,
C              FALLS DER PFEIL DICKER GEZEICHNET WERDEN SOLL

                  IF(ABSK.EQ.3.OR.ABSK.EQ.4.OR.ABSK.EQ.6.OR.ABSK.EQ.7)
     $            THEN
                     IF (ABSK .EQ. 3 .OR. ABSK.EQ.6) THEN
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

C                 AUFRUF DER GKS-ROUTINE    **********

                  IF (ABSK.LE.5) THEN
                     CALL GPL (ISTUE,XX,YY)
                  ELSE
                     CALL GRFILL(ISTUE,XX,YY,1,1)
                  ENDIF

               ENDIF
              ICX = ICX + 1
              IF (ICX .GT. NX) THEN
                 ICX = 1
                 IF (ABS(IPOLAR).NE.1234567890) THEN
                    XWO = PP(1) + XAB*.5
                    YWO = YWO + YAB
                 ELSE
                    IR=IR+1
                    RAD=(PP6+DELR*(.5+IR))*FR
                    U=(PP5+DELPHI*.5)*PIF
                    XWO = XMIT+RAD*COS(U)
                    YWO = YMIT+RAD*SIN(U)
                 ENDIF
              ELSE
                 IF (ABS(IPOLAR).NE.1234567890) THEN
                    XWO = XWO + XAB
                 ELSE
                    U=(PP5+DELPHI*(ICX-.5))*PIF
                    XWO = XMIT+RAD*COS(U)
                    YWO = YMIT+RAD*SIN(U)
                 ENDIF
              ENDIF

   30      CONTINUE
   40    CONTINUE
      ENDIF

C     ZURUECKSETZEN DER SKALIERUNG

      CALL GRSCLV (PP5,PP6,PP7,PP8)

      END
