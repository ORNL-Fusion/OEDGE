C@PROCESS NOSDUMP NOGOSTMT OPT(3)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C     changed: 29. 5. 90 Marlene Busch
C     GRTXTSCN eingefuegt - Umlaute und sz auf GTSGRAL GKS umgestellt
C     changed: 18.10.90  Marlene Busch
C     changed:  5.09.91 Gerd Groten Transformationen
C     changed:  6.10.91 Gerd Groten :90 grad und sehr verschiedener
C                                    Achsenmassstab korrigiert.
C              query text extent - im CMS ist dies die Dummy Workstation
C                                  4711 mit dem WK Indentifier 2
C                                - sonst angenommen, das die Workstation
C                                  die Workstation ID 1 hat
C Update 14.11.91 Busch GQTXX, nicht fuer CGM WK und NICHT GTSGRAL GKS
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
* Update: 18. 3.92 Busch COMMON GRTXPO eingefuegt, current position
*                        fuer Text neu (frueher: nur eine curr. pos.
*                        fuer Text und Linien
* Update 28.1.94 Busch beruecksichtigen, dass zuvor gstxal gesetzt
C ----------------------------------------------------------------------
      SUBROUTINE GRTXT(X,Y,LTEXT,TEXT)

      REAL X,Y
      INTEGER   LTEXT
      CHARACTER(len=*) TEXT

      INTEGER IWKID,LTE, LTE1
      PARAMETER ( PIFAC=6.2831853072/360)
      REAL    XEND,YEND,XRECT(4),YRECT(4),SCLV(4),SCLC(4)
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      CHARACTER  (len=256) TEXT1
      CHARACTER (len=7) GKSTYP
      LOGICAL IF1,IF2

C---- COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /SCALE/ XMAXDC,XUNITS,YUNITS
CDEC$ PSECT /SCALE/ NOSHR
      COMMON /GRTXPO/ XTCURR,YTCURR
CDEC$ PSECT /GRTXPO/ NOSHR
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
C 2.9.91
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRPIC/ FLPIC,NSCLC,NSCLV, NSCLP, RAHMEN,
     $               XMAXCM,YMAXCM, XDCPIC,YDCPIC
CDEC$ PSECT /GRPIC/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRPP/ ,/SCALE/,/GRTXPO/ , /GRREST/, /GRPIC/ ,/GRGKS/

C     query number of Workstations
C     falls der Dummy Treiber mit WK Typ 4711 existiert ist man im
C     CMS und IWK=2 wird gesetzt, da query text extent auf die dummy
C     wk gemacht wird
C     sonst wird angenommen, das GKS auf der CRAY oder einem anderen
C     GKS laeuft und nur eine WK mit WKID 1 offen ist
C
      CALL GQEWK(1,IERR,NUMBER,IWKT)
      DO 40 I=1,NUMBER
      CALL GQEWK(I,IERR,NUMBER,IWKT)
c     testweise , da Aerger mit Versatec Routinen
      IF (IWKT.GT.16000) GOTO 40
      IF ( IWKT.EQ.4711 .and. ierr.eq.0 ) THEN
         IWKID=2
         goto  50
      ENDIF
 40   CONTINUE
      IWKID=1
 50   CONTINUE

c 14.11.91
C      query wk connection and type
      CALL GQWKC(IWKID,IERR,ICon,IWT)
      CALL GQWKCA (IWT,IERR, ICAT)

C---- Nimm bei unmoeglichem LTEXT die Stringlaenge
      IF ( LTEXT.LT.0 .OR. LTEXT.GT.LEN(TEXT) ) THEN
         LTE=LEN(TEXT)
      ELSE IF (LTEXT.EQ.0) THEN
         GOTO 99
      ELSE
         LTE=LTEXT
      ENDIF

C---- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

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

C---- WINKEL FUER DIE TEXTAUSGABE
      SCFCTR=XUNITS/YUNITS
      IWIN=PP(15)+0.5
      IF1=ABS(SCFCTR-1.).LE.0.1
C 6.10.91
      IF2=IWIN/90*90.EQ.IWIN .AND. SCFCTR.LT.1E5 .AND. SCFCTR.GT.1E-5
C 2.9.91
      IF ( IF1.OR.IF2 ) THEN
C------- EINFACHER TEXTAUFRUF IM ORIGINAL-WINDOW
         CALL GRTXSCN ( TEXT, LTE, TEXT1 ,LTE1 )
         CALL GTX( X,Y,TEXT1(:LTE1))

C------- CURRENT POSITION MERKEN
c 14.11.91
c        nicht gqtxx , wenn NICHT GTSGRAL-GKS und MO- wktyp
         IF (iwkid.eq.1.and.ICAT.EQ.4.and.gkstyp.ne.'GLIGKS') goto 100
         CALL GQTXX (IWKID,X,Y,TEXT1(:LTE1),IERR,XEND,YEND,XRECT,YRECT)

C        IERR=39 ,wenn GQTXX fuer GKSM aufgerufen wird (GKS Version 3.2)
C        IERR=0  ,wenn GQTXX fuer CGM  aufgerufen wird (GKS Version 3.2)
         IF ( IERR .EQ. 39 )  THEN
            XEND=X
            YEND=Y
C 2.9.91
         endif

         XTCURR=XEND
         YTCURR=YEND
 100     continue
      ELSE
C------- TEXTAUSGABE NUR IM GEAENDERTEN WINDOW MOEGLICH
C------- WINDOW WIRD AUF CM ZURUECKGESETZT
         CALL GQCNTN(IERR,NUMTR)
         CALL GQNT(NUMTR,IERR,SCLV,SCLC)
         CALL GQCHH(IERR,CHH)
         CALL GQCHUP(IERR,CHUX,CHUY)
         CALL GQCHXP(IERR,CHXP)
C------- NEUES WINDOW IN CM , NUMTR=2 Clipping
C 5.9.91
         IF (NUMTR.EQ.1) THEN
            CALL GSWN(NUMTR,0.,XMAXCM,0.,YMAXCM)
         ELSE
            CALL GSWN(NUMTR,PP(1),PP(3),PP(2),PP(4))
         ENDIF
C------  CHARACTERGROESSE IN CM
         CHCM=PP(14)
         CALL GSCHH(CHCM)
C------  WINKEL
         WINKEL=PP(15)
         CHANG=WINKEL*PIFAC
         SC=SIN(CHANG)
         CC=COS(CHANG)
         CALL GSCHUP(-SC,CC)
C------- EXPANSIONFAKTOR AUF 1 SETZEN
         CALL GSCHXP(1.)
C 2.9.91
         XCOORD=PP(1)+(X-PP(5))/XUNITS
         YCOORD=PP(2)+(Y-PP(6))/YUNITS

         CALL GRTXSCN ( TEXT, LTE, TEXT1 ,LTE1 )
         CALL GTX(XCOORD,YCOORD,TEXT1(:LTE1))

C------- CURRENT POSITION MERKEN
C 14.11. 91
         IF (iwkid.eq.1.and.ICAT.EQ.4.and.gkstyp.ne.'GLIGKS') goto 200
         CALL GQTXX(IWKID,XCOORD,YCOORD,TEXT1(:LTE1),IERR,XEND,YEND,
     >              XRECT,YRECT)
C        IERR=39 ,wenn GQTXX fuer GKSM aufgerufen wird (GKS Version 3.2)
C        IERR=0  ,wenn GQTXX fuer CGM  aufgerufen wird (GKS Version 3.2)
         IF ( IERR .EQ. 39 )  THEN
            XEND=XCOORD
            YEND=YCOORD
         endif
 200     continue
C------- ALLE ORIGINALWERTE MUESSEN WIEDER EINGESETZT WERDEN
         CALL GSWN(NUMTR,SCLV(1),SCLV(2),SCLV(3),SCLV(4))
         CALL GSCHH(CHH)
         CALL GSCHXP(CHXP)
         CALL GSCHUP(CHUX,CHUY)
C 2.9.91
         XTCURR=PP(5)+(XEND-PP(1))*XUNITS
         YTCURR=PP(6)+(YEND-PP(2))*YUNITS
      ENDIF
  99  continue
      end
