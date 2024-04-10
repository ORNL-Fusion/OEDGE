CCCC@PROCESS OPT(3) NOSDUMP NOGOSTMT
***********************************************************************
** *** COPYRIGHT 1990 BY FORSCHUNGSZENTRUM (KFA) JUELICH GMBH    ***  *
**                                                                    *
** NAME:    GREND                                                     *
**                                                                    *
** ZWECK:   ENDROUTINE FUER DIE CRAY VERSION DER GR-SOFTWARE          *
**                                                                    *
** AUTOR :   M. BUSCH,        (ZDV131@DJUKFA11)                       *
**                                                                    *
** DATUM :   25. 9. 1990                                              *
** UPDATE:    6.11. 1990   SAVE                                       *
* UPDATE : 8.10.91 M. BUSCH  GRSCLP wird in jedem falle in GRNXTF
*                           aufgerufen - siehe Note und auch Aenderung
*                           in GRSCLP zu dem Datum - hier faellt Aufruf
*                           von GRRAHM weg
* Update: 13.11.91 Busch request  choice nur fuer Terminal
* Update : 15.11.91 Busch query segmtents on wk nicht fuer wk=mo
* Update: 26. 2.92 Busch COMMON GRREST eingefuegt
* Update: 13. 8.93 Busch Anpassuing an verschiedene GKS Implementatione
* Update: 26.1.93 Busch Ausgabe evtl. vorhandener Marker
***********************************************************************
***********************************************************************
C
      SUBROUTINE GREND
      IMPLICIT NONE
      REAL XLN(300), YLN(300), YMR(300), XMR(300), pp
CCC   COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
CCC
CCC
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      INTEGER INTLIN, ERR, FLPIC, RAHMEN,IO,ISYS,nsclv,nsclp,nrln,nrmr
      INTEGER idum,iwkwis,nsg,name,i,icolor,ifont,ixgks,maxpkt,nsclc
      INTEGER ierr,isegnr,iwk,icon,iwt,icat
      PARAMETER (ISYS=91,IO=92)
      REAL MOUT(2,3),xcurr,ycurr,xdcpic,ydcpic,sizmrk,xmaxcm,ymaxcm
      CHARACTER GKSTYP*7, STR*1
      LOGICAL FLGROT
      EQUIVALENCE (PP(13),INTLIN), (PP(16),ICOLOR), (PP(17),SIZMRK)
      EQUIVALENCE (PP(9),IFONT), (PP(18),FLGROT)
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

      SAVE /GRPP/, /GRWK/,/GRGKS/ ,/GRREST/, /GRPIC/

CCC   Entfernen von temporaeren Dateien

      CLOSE(ISYS,STATUS='DELETE')
      CLOSE(IO,STATUS='DELETE')

CCC   AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN
CCC
      IF (NRLN .GT. 1) THEN
        CALL GPL(NRLN, XLN, YLN)
        XCURR = XLN(NRLN)
        YCURR = YLN(NRLN)
        NRLN = 0
      END IF
CCC   AUSGABE EVT. VORHANDENER GRJMPS-DATEN
CCC
      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF

C     in jedem Falle bei grnxtf einen optionalen Rahmen zeichen, da man
C     nicht weiss, ob Benutzer erneut GRSCLP,GRSCLC,GRSCLV aufruft
c     und der Viewport, wenn nichts anderes gesagt, erhalten bleibt

c 8.10.91
c     IF ( NSCLP.EQ.0.AND. RAHMEN .GT. 0 ) CALL GRRAHM


CCC     WENN FLGROT='Y', DANN MUSS DAS BILD GEDREHT WERDEN:
      IF (FLGROT) CALL GEVTM(0.5, 0.5, 0., 0., 3.141592*0.5, 1., 1., 1,
     1   MOUT)



      IF (GKSTYP .EQ. 'GLIGKS') THEN
C        Setzen der Segmenttransformation zum Drehen
        IF (FLGROT) CALL GSSGT(1, MOUT)

        CALL GRPAN(1)
        CALL GECLKS

      ELSE

CCC     INQUIRE NAME OF OPEN SEGMENT
        CALL GQOPSG(IERR, ISEGNR)

        CALL GCLSG

c     query wk connection and type
        CALL GQWKC(IWK, IERR, ICON, IWT)
c     qeury wk cat
        CALL GQWKCA(IWT, IERR, ICAT)

        IF (GKSTYP .EQ. 'GRALGKS') THEN
          IF (FLGROT) CALL GSSGT(1, MOUT)
          IF (ICAT .EQ. 2) CALL GRQST(IWK, 1, IDUM, 1, STR)

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

          IF (ICAT .EQ. 2) CALL GRQCH(IWK, 1, IDUM, IDUM)

          CALL GDAWK(IWKWIS)

          CALL GCLWK(IWKWIS)
        END IF

        CALL GDAWK(IWK)
        CALL GCLWK(IWK)

        CALL GCLKS

      END IF


      GO TO 4711
CCC   ENTRY NUR FUER DEN SPEZIALFALL, DASS KEIN BILD MEHR GEPLOTTET
CCC   WIRD, SONDERN NUR ALLE WORKSTATIONS ABGESCHLOSSEN WERDEN.
CCC   Z.B. WIRD ENTRY BEI CRAY-DATEN (CRAYPL) BENUTZT.
CCC
      ENTRY GRENDE

CCC   Entfernen von temporaeren Dateien

      CLOSE(ISYS,STATUS='DELETE')
      CLOSE(IO,STATUS='DELETE')

      CALL GECLKS
 4711 END
