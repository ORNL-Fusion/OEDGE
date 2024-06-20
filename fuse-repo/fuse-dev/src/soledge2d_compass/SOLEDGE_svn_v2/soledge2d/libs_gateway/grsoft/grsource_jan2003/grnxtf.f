CCCCCC@PROCESS OPT(3) NOSDUMP NOGOSTMT
***********************************************************************
** *** COPYRIGHT 1990 BY FORSCHUNGSZENTRUM (KFA) JUELICH GMBH    ***  *
**                                                                    *
** NAME:    GRNXTF                                                    *
**                                                                    *
** ZWECK:   NEXT FRAME ROUTINE FUER FUER CRAY VERSION DER GR-SOFTWARE *
**                                                                    *
** AUTOR :   M. BUSCH,        (ZDV131@DJUKFA11)                       *
**                                                                    *
** DATUM :   25. 9. 1990                                              *
** UPDATE:   20.12. 1990 GRSCLP BELIEBIG BILDGROESSE                  *
* UPDATE :22.07.91 M. Busch
C     Default fuer GRWIN ist: fuer jedes Bild muss neuer Bildausschnitt
c     definiert werden (historisch)
* UPDATE : 8.10.91 M. BUSCH  GRSCLP wird in jedem falle in GRNXTF
*                           aufgerufen - siehe Note und auch Aenderung
*                           in GRSCLP zu dem Datum - hier faellt Aufruf
*                           von GRRAHM weg
* UPDATE :21.10.91 M. BUSCH  statt grsclp wird entry grrahm gerufen
* Update : 13.11.91 Busch request choice nur fuer Terminal
* Update : 15.11.91 Busch query segmtents on wk nicht fuer wk=mo
* Update: 26. 2.92 Busch COMMON GRREST eingefuegt
* Update:  8. 2.93 Busch UNIX  COMMON GRCOM1 neu
* Update: 13. 8. 93 Busch Anpassung an verschidene GKS Implementationen
*                         unter UMIX
* Update 26.1.94 Busch evtl. vorhandene Marker ausgeben
***********************************************************************
C
      SUBROUTINE GRNXTF
      REAL XLN(300), YLN(300), YMR(300), XMR(300)
CCC   COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
CCC
CCC
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      INTEGER INTLIN,  ERRUN
      INTEGER FLPIC, RAHMEN, ERR
      INTEGER IGRSHOW, IGR3, SHOWPR, IGR3PL
      LOGICAL FLGROT
      CHARACTER GKSTYP*7, STR*1
      REAL MOUT(2,3)
      EQUIVALENCE (PP(13),INTLIN), (PP(16),ICOLOR), (PP(17),SIZMRK)
      EQUIVALENCE (PP(9),IFONT), (PP(18),FLGROT)
CCC
CCC
      COMMON /GRWK/ IXGKS, IWK, IWKWIS
CDEC$ PSECT /GRWK/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      COMMON /GRREST/ MAXPKT, NRLN, XCURR, YCURR, XLN, YLN, NRMR, XMR,
     1       YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRPIC/ FLPIC, NSCLC, NSCLV, NSCLP, RAHMEN, XMAXCM, YMAXCM,
     1       XDCPIC, YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /GRWIND/ X1, X2, Y1, Y2, IWIN
CDEC$ PSECT /GRWIND/ NOSHR
      COMMON /GRCOM1/ IGRSHOW, ERRUN, IGR3, SHOWPR, IGR3PL
CDEC$ PSECT /GRCOM1/ NOSHR
      SAVE /GRPP/,/GRWK/ ,/GRGKS/,/GRREST/,/GRPIC/ ,/GRWIND/,/GRCOM1/


CCC   AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN
CCC
      IF (NRLN .GT. 1) THEN
        CALL GPL(NRLN, XLN, YLN)
        NRLN = 0
      END IF

CCC   AUSGABE EVT. VORHANDENER GRJMPS-DATEN
CCC
      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF

      XCURR = 0.
      YCURR = 0.

      SHOWPR = 1

C     in jedem Falle bei grnxtf einen optionalen Rahmen zeichen, da man
C     nicht weiss, ob Benutzer erneut GRSCLP,GRSCLC,GRSCLV aufruft
c     und der Viewport, wenn nichts anderes gesagt, erhalten bleibt


c 8.10.91
c     IF ( NSCLP.EQ.0.AND. RAHMEN .GT. 0 ) CALL GRRAHM


CCC     WENN FLGROT='Y', DANN MUSS DAS BILD GEDREHT WERDEN:
      IF (FLGROT) CALL GEVTM(0.5, 0.5, 0., 0., 3.141592*0.5, 1., 1., 1,
     1   MOUT)

      IF (GKSTYP .EQ. 'GLIGKS') THEN
        CALL GSWKWN(IWK, 0., XDCPIC, 0., YDCPIC)
C       Setzen der Segmenttransformation zum Drehen
        IF (FLGROT) CALL GSSGT(1, MOUT)
        CALL GRPAN(0)

      ELSE

        CALL GCLSG

C     query wk connection and type
        CALL GQWKC(IWK, IERR, ICON, IWT)
C     qyery wk cat.
        CALL GQWKCA(IWT, IERR, ICAT)


        IF (GKSTYP .EQ. 'GRALGKS') THEN
          IF (FLGROT) CALL GSSGT(1, MOUT)
          IF (ICAT .EQ. 2) CALL GRQST(IWK, 1, IDUM, 1, STR)
          CALL GCLRWK(IWK, 1)
          CALL GCRSG(1)

        ELSE


          CALL GACWK(IWK)
C         nicht fuer WISS
          IF (ICAT.NE.3) CALL GSWKWN(IWK, 0., XDCPIC, 0., YDCPIC)

c 15.11.91
C ICAT=(output,input,inout,wiss,mi,mo)
c     IF (ICAT.EQ.4) goto 100

          CALL GQSGWK(IWKWIS, 1, ERR, NSG, NAME)
          DO 3, I = 1, NSG
            CALL GQSGWK(IWKWIS, I, ERR, IDUM, NAME)
C        Setzen der Segmenttransformation zum Drehen
            IF (FLGROT) CALL GSSGT(NAME, MOUT)
            CALL GASGWK(IWK, NAME)
    3     CONTINUE

  100     CONTINUE

          CALL GUWK(IWK, 1)

          IF (ICAT .EQ. 2) CALL GRQCH(IWK, 1, IDUM, IDUM)

          CALL GDAWK(IWK)

          CALL GCLRWK(IWK, 1)
          CALL GCLRWK(IWKWIS, 1)


C****** wegen XGKS Fehler Loeschen aller Segmente
C****** beim GCLRWK auf CGMO werden Segmente nicht geloescht
C****** -> explizit loeschen
c     IF (ICAT.EQ.4) THEN


c     CALL GQSGWK ( IWKWIS, 1, ERR, NSG, NAME )
c      WRITE(*,*) 'NSG=', NSG
c     IF (NSG.EQ.0) goto 5
c     do 4 i=1,NSG
c        CALL GQSGWK ( IWKWIS, I, ERR, IDUM, NAME )
c        IF ( ERR.EQ.0) CALL GDSG(NAME)
c4    CONTINUE
c     ENDIF

    5     NAME = NAME + 1

          CALL GCRSG(NAME)

        END IF

      END IF

CCC   DEFAULTWERTE WINDOW UND VIEWPORT FEHLEN NOCH


C     GRSCLP gilt nach dem Aufruf unveranedert fuer das gesamte
C     Programm ( nach neu gerufen werden !)
C     GRSCLV,GRSCLC nach GRNXTF als erster Aufruf dann
C     wird, falls GRSCLP Falg FLPIC gesetzt GRSCLP mit den
C     alten Werten gerufen und gegebenfalls Rahmen gezeichent

C     Default fuer GRWIN ist: fuer jedes Bild muss neuer Bildausschnitt
c     definiert werden (historisch)

      IWIN = 0

C     in jedem Falle bei grnxtf einen optionalen Rahmen zeichen, da man
C     nicht weiss, ob Benutzer erneut GRSCLP,GRSCLC,GRSCLV aufruft
c     und der Viewport, wenn nichts anderes gesagt, erhalten bleibt

c 18.10.91 entry zu grsclp
      IF (FLPIC .EQ. 1234567890) CALL GRRAHM


C     ZAEHER GRSCLC AUFRUFE (BILD)
      NSCLC = 0
      NSCLP = 0
      NSCLV = 0

      SHOWPR = 0

      IGRSHOW = 0
      IGR3 = 0
      IGR3PL = 0
      RETURN

      END
