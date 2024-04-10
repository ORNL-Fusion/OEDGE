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
C     max. 200 * 83 cm
C     Groesse kann jedoch beliebig sein - Versatec schneidet ab
C     skaliert Zeichenfeld  mit den cm-Angaben (ueberschreibt  damit
C     GRSCLC (0,0),(39.5,28.7)
C     skaliert Wertebereich  mit den cm-Angaben (ueberschreibt  damit
C     GRSCLV (0,0),(39.5,28.7)
C     Bild wird in jedem Falle verzerrungsfrei auf dem Bildschirm
C     dargestellt
C     VERSATEC Farbelektrostat: es gelten die gemachten cm Angaben
C              83cm ist Limit in der Hoehe
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
C-----------------------------------------------------------------------

      Subroutine  GRRAHM

      INTEGER FLPIC, RAHMEN, ICHECK
      REAL    XMAXCM ,YMAXCM

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      REAL    X(5),Y(5)
      LOGICAL FLGROT

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
      SAVE /GRREST/, /SCALE/, /GRPP/, /GRPIC/

      EQUIVALENCE (PP(16),INTSYM)

C     falls GRSCLP fuer mehrere Bilder gilt, d.h. nach GRNXTF nicht
c     explizit aufgerufen wird, wird GRRAHM von GRENV gerufen

ctest FLPIC=0
      IBOX=RAHMEN
      XCM=XMAXCM
      YCM=YMAXCM
c     NSCLP=NSCLP+1
      ICHECK=1


C     RAHMEN WENN >=1
      IF (IBOX.GE.1.AND.IBOX.LE.4) THEN
         RAHMEN=IBOX
           ELSEIF ( IBOX.EQ.0) THEN
           RAHMEN=0
             ELSE
             RAHMEN=1
      ENDIF

      IF ( YCM.GT.83) THEN
         YMAXCM=83.
         WRITE(*,*) 'Bild wird in der Hoehe abgeschnitten max. 83 cm'
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
c     IF ( FLPIC.EQ.1234567890.and.ICHECK.EQ.0) THEN
c        PP(1)=0.
c        PP(2)=0.
c        PP(3)=XMAXCM
c        PP(4)=YMAXCM
c        PP(5)=0.
c        PP(6)=0.
c        PP(7)=XMAXCM
c        PP(8)=YMAXCM
c     ENDIF

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
