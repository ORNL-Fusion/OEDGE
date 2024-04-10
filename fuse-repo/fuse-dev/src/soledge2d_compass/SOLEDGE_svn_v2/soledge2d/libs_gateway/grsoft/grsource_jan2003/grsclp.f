C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C-----------------------------------------------------------------------
C     GRSCLP
C     GR-Software scale paper size
C     AUTHOR: Marlene Busch
C     Date:   14. 12. 1990
c
C     muss direkt nach GRSTRT bzw. GRNXTF gerufen werden
c     ein zweiter Aufruf innerhalb eines Bildes wird ignoriert
C
C     legt Groesse des Zeichenfeldes in cm fest (Breite*Hoehe)
C     max. 200 * 130 cm
C     Groesse kann jedoch beliebig sein - Versatec schneidet ab
C     skaliert Zeichenfeld  mit den cm-Angaben (ueberschreibt  damit
C     GRSCLC (0,0),(39.5,28.7)
C     skaliert Wertebereich  mit den cm-Angaben (ueberschreibt  damit
C     GRSCLV (0,0),(39.5,28.7)
C     Bild wird in jedem Falle verzerrungsfrei auf dem Bildschirm
C     dargestellt
C     Drucker: es gelten die gemachten cm Angaben
C              130cm ist Limit in der Hoehe
C              200cm ist vorlaeufig (Testzweck!) Limit in der Breite{
C     IBOX>1 --> RAHMEN wird gezeichent ( Abmessen Versatec)
C                =1 durchgezogen,=2 gestrichelt, =3 gepunktet,=4 strichp
C     IBOX=0     kein Rahmen
C     Update: 22.7.91 Busch, Fehler mit Rahmen beim 2.bis n-ten Bild
C                     XMAXCM,YMAXCM wurde genommen - wenn Benutzer
C                     jedoch mit GRSCLV Werte-Bereich umgesetzt hat
C                     muss man dafuer sorgen dass unabh. davon Rahmen
C                     immer um das gesamte Bild gezeichnet wird
C                     --> XMAXXM,YMAXCM auf gesamten Bereich  umrechnen
C     Update: 08.10.91 Busch siehe Note im Text
C                      betrifft Geltungbereich von GRSCLP, wenn nicht
C                      fuer jedes Bild aufgerufen wird
C     Update: 17.10.91 Busch siehe Note im Text
C                      GRSCLC, GRSCLV einmal zu Anfang-> Fehler ab
C                      zweiten Bild wird zu klein gezeichnet
C                      NEU: von GRNXTF wird ENTRY GRRAHM gerufen
C                           Steuerung dadurch verbessert
C                           FLAG FLPIC wird zurueckgesezt in ENTRY
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C Update 17.2.94 Busch Unix Version im immediate mode muss Fenster-
C                groesse (ueber GRSCLP gesetzt) mitkriegen
C                Im immediate mode wird standardmaessig nur 39.5/28.7
C                eroeffnet
C                dagegen im normalen Arbeitsmodus ( interaktive Abfrage
C                nach der Ausgabe) ->
C                Beim Arbeiten mit WISS wird gswkwn und gswkvp in
C                grgks.c gesetzt
C Update 30.8.96 GSWKVP masstabgetreues Bild
C                nur mit Vorwahl der Ausgabe moeglich
C                GSWKVP wird genau metrisch gesetzt 
C
C                sonst wird entweder das Standardfenster oder bei
C                Setzen von GRSCLP ein Fenser in dem gewuenschten
C                Format erstellt
C Update: 19.9.96 Plotter erlaubt max 130cm in der Hoehe, die Breite
C                 wurde mit 200 cm beibehalten

C-----------------------------------------------------------------------
      SUBROUTINE GRSCLP(XCM, YCM, IBOX )

      INTEGER FLPIC, RAHMEN, ICHECK
      REAL    XMAXCM ,YMAXCM ,XMAX,YMAX

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      REAL    X(5),Y(5)
      LOGICAL FLGROT
      CHARACTER GKSTYP*7

      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
C---- COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(17) , FLGROT
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRREST/,/SCALE/ ,/GRPP/,/GRPIC/ , /GRGKS/
      EQUIVALENCE (PP(16),INTSYM)

      ICHECK=0

c     WRITE(3,*) '--GRSCLP--','Xcm ,Ycm ',Xcm ,Ycm
c     WRITE(3,*) '--GRSCLP--','nsclp=', nsclp

C     ZAEHLER AUFRUFE GRSCLP
      NSCLP=NSCLP+1

C     GRSCLP wird in jedem Falle in GRNXTF aufgerufen
C     (und zwar ueber den ENTRY GRRAHM seit 17.10.91)
C     damit auch dann wenn Benutzer fuer weitere Bilder GRSCLP
C     unveraendert gelten lassen moechte ein Rahmen gezeichnet wird
C
C     falls Benutzer selbst dann GRSCLP aufruft, wird das bis dahin er-
C     zeugte Segment weggeworfen (das ist der Fall wenn NSCLP=2 ist)
C
c     WRITE(*,*) 'grsclp -NSCLP,NSCLC,NSCLV', NSCLP,NSCLC  ,NSCLV
c     WRITE(*,*) 'grsclp-FLPIC=', FLPIC
c     read(*,*) ijkl


C 8.10.91
c     IF ( NSCLP.EQ.2) THEN
      IF ( NSCLP.EQ.1) THEN
         NSCLP=NSCLP-1
         CALL GQOPSG(IERR,NAME)
         IF (IERR.EQ.0) THEN
            CALL GCLSG
            CALL GDSG(NAME)
            CALL GCRSG(NAME)
         ENDIF
      ENDIF



C     Aufrufe von GRSCLP innerhalb eines Bildes werden ignoriert-
C     nur der erste Aufruf zaehlt- dieser muss vor GRSCLC bzw.GRSCLV
c     stehen
C     FLPIC=1234567890 - GRSCLP wurde fuer ein Bild bereits aufgerufen
C     NSCLC=0 noch kein Aufruf von GRSCLC fur das Bild
C     in GRSTRT,GRNXTF wird FLPIC=0 gesetzt

      IF ( FLPIC.EQ.1234567890 ) THEN
C        GRSCLP AUFGERUFEN, ABER GRSCLC VORHER --> GRSCLP INGNORIEREN
         IF ( NSCLC.GT.0)  RETURN
C        WENN NACH GRNXTF GRSCLP ALS ERSTES AUFGERUFEN WIRD - NEUE
C        PAPER SKALIERUNG
         IF (NSCLC.EQ.0.AND.NSCLV.EQ.0.AND.NSCLP.EQ.1) GOTO 555
C        GRSCLP AUFGERUFEN, ABER GRSCLC VORHER --> GRSCLP INGNORIEREN
      ENDIF

 555  FLPIC=1234567890
C     WRITE(3,*) 'grsclp-FLPIC=', FLPIC



C     RAHMEN WENN >=1
      IF (IBOX.GE.1.AND.IBOX.LE.4) THEN
         RAHMEN=IBOX
           ELSEIF ( IBOX.EQ.0) THEN
           RAHMEN=0
             ELSE
             RAHMEN=1
      ENDIF

      IF ( YCM.GT.130) THEN
         YMAXCM=130.
         WRITE(*,*) 'Bild wird in der Hoehe abgeschnitten max. 130 cm'
      ELSE
         YMAXCM=YCM
      ENDIF

      IF ( XCM.GT.200) THEN
         XMAXCM=200.
         WRITE(*,*) 'Bild wird in der Breite abgeschnitten max. 200 cm'
         WRITE(*,*) 'Falls Sie eine laengeres Bild erstellen wollen   '
         WRITE(*,*) 'melden Sie sich bitte in der Programmberatung ZAM'
      ELSE
         XMAXCM=XCM
      ENDIF

c17.10.91 abfrage neu
cmb   IF ( FLPIC.EQ.1234567890.and.ICHECK.EQ.0) THEN
      IF ( FLPIC.EQ.1234567890) THEN
         PP(1)=0.
         PP(2)=0.
         PP(3)=XMAXCM
         PP(4)=YMAXCM
         PP(5)=0.
         PP(6)=0.
         PP(7)=XMAXCM
         PP(8)=YMAXCM
      ENDIF

CCC   AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

      IF (NRLN.GT.1) THEN
         CALL GPL(NRLN,XLN,YLN)
         XCURR=XLN(NRLN)
         YCURR=YLN(NRLN)
         NRLN=0
      ENDIF
      IF (NRMR.GT.0) THEN
         CALL GPM(NRMR,XMR,YMR)
         NRMR=0
      ENDIF

c Festlegung des Breite/ Hoehe Verhaeltnisses des Bildes

      IF ( XMAXCM.GE.YMAXCM) THEN
         XDCPIC=1
         YDCPIC=YMAXCM/XMAXCM
      ELSE
         YDCPIC=1
         XDCPIC=XMAXCM/YMAXCM
      ENDIF

C Viewport
      CALL GSVP(1,0.,XDCPIC,0.,YDCPIC)
c     write(3,*) 'grsclp-viewport', xdcpic,ydcpic
c     write(3,*) 'grsclp-xmaxcm,ymaxcm', xmaxcm,ymaxcm


C---- IN DEN ROUTINEN  GRSCLC, GRSCLV MUSS DAS WINDOW
C---- NEU BERECHNET WERDEN.
C---- EINHEITEN/CM : UNITS=VALUES/CM
C     Hier noch richtig, sich auf  PP(1..8) zu beziehen, da diese
C     gleich den cm in X-Richtung und Y-Richtung
c     gesetzt wurden


calt  XUNITS=  (PP(7)-PP(5))/XMAXCM
calt  YUNITS=  (PP(8)-PP(6))/YMAXCM
      XUNITS=  (PP(7)-PP(5))/(PP(3)-PP(1))
      YUNITS=  (PP(8)-PP(6))/(PP(4)-PP(2))

c     write(3,*) 'grsclp- xunits,yunits', xunits,yunits

C---- UMRECHNUNG DES WINDOW AUF DEN GESAMTEN VIEWPORTBEREICH

      XLI=PP(5)-PP(1)*XUNITS
      YLI=PP(6)-PP(2)*YUNITS

      XRE=PP(7)+(XMAXCM-PP(3))*XUNITS
      YRE=PP(8)+(YMAXCM-PP(4))*YUNITS

C---- WINDOW AUF NDC Viewport  ABBILDEN
      CALL GSWN(1,XLI,XRE,YLI,YRE)
c     write(3,*)'grsclp-gswn', XLI,XRE,YLI,YRE

C Abfrage nur fuer UNIX Version von Beseutung, da COMMON GRGKS
C in CMS Version mit GKSENV FORTRAN nicht gesetzt wird
C 17. 2.94
      IF ( GKSTYP.EQ.'GLIGKS'.OR.GKSTYP.EQ.'GRALGKS') THEN

C--- offene wk's
         CALL GQOPWK(1,ierr,nopn,iwk)
         do 1 i=1,nopn
            CALL GQOPWK(i,ierr,nopn,iwk)
            CALL GQWKC ( IWK, IERR,ICON, IWT)
            CALL GQWKCA(IWT,IERR,ICAT)
            call gqdsp(IWT,ierr,idun,rx,ry,lx,ly)
            ratio = ymaxcm / xmaxcm
            if (iwt.eq.217) then
                if (ratio .lt.1) then
                   rx = rx * (640.0 / lx / ratio)
                   ry = rx *ratio
                else
                   ry = ry * (640.0 / ly * ratio)
                   rx = ry / ratio
                endif
            elseif ( iwt.ge.210.and.iwt.le.216) then
                 if (rx .gt. ry ) then
                    rx=rx*0.95
                 else
                    ry=ry*0.95
                 endif
            endif
C fit viewport
            if ( rx.gt.xmaxcm.and.ry.gt.ymaxcm) then
               rx=xmaxcm
               ry=ymaxcm
            else
               xratio = xcm / rx
               yratio = ycm / ry
               if (yratio .gt.xratio) then
                  rx = xcm / yratio
                  ry = ycm / yratio
               else
                  rx = xcm / xratio
                  ry = ycm / xratio
               endif
            endif

            xmin=0.
            ymin=0.
            xmax=rx
            ymax=ry
C           CALL GSWKVP (iwk,xmin,xmax,ymin,ymax)
            CALL GSWKVP (iwk,xmin,xmaxcm/100.,ymin,ymaxcm/100.)
            if ( ratio.lt.1) then
               xmax = 1.0
               ymax = ratio
            else
               xmax = 1.0 /ratio
               ymax = 1.0
            endif
            CALL GSWKWN (iwk,xmin,xmax,ymin,ymax)
1        continue

      endif
C Ende (GLIGKS)

      CALL GRCHRC(PP(14),PP(15),INTSYM)

C  BILD MARKIEREN - RAHMEN
      IF ( RAHMEN.GT.0) THEN
C QUERY POLYINE COLOR
         CALL GQPLCI(IERR,ICOL)
C QUERY POLYINE TYP
         CALL GQLN(IERR,IPLTYP)
         CALL GSPLCI(1)
         CALL GSLN  ( RAHMEN )
C RAHMEN
         X(1)=XLI
         X(2)=XRE
         X(3)=XRE
         X(4)=XLI
         X(5)=XLI
         Y(1)=YLI
         Y(2)=Yli
         Y(3)=YRE
         Y(4)=YRE
         Y(5)=YLI
         CALL GPL (5, X,Y)

C RESTORE POLYLINE COLOR
         CALL GSPLCI(ICOL)
C RESTORE POLYLINE TYP
         CALL GSLN (IPLTYP)
      ENDIF

999   CONTINUE
      END


      SUBROUTINE GRQSCL (XCM, YCM)

      REAL XCM, YCM

      INTEGER FLPIC, RAHMEN

      COMMON /GRPIC/ FLPIC, NSCLC, NSCLV, NSCLP, RAHMEN,
     *               XMAXCM, YMAXCM, XDCPIC, YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR

      XCM = XMAXCM
      YCM = YMAXCM
      END
