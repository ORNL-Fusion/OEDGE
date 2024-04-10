C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRNWPN(IPEN)

C---- Farben sind bei GTSGRAL anders definiert und zwar:
C---- 1=weiss,2=rot,3=gruen, 4=blau, 5=gelb, 6=pink, 7=tuerkis, 8=schwar
C---- Datum: 31.5 90. Marlene Busch
C---- Update 26.7 91. Marlene Busch IFONT--> PP(16)
C---- Update 27.11.91 Marlene Busch Farben umsetzen nur fuer GTSGRAL GKS
C                     wenn Dummy WK 4711 dan GTSGRAL
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C UPDATE 12. 5.1992 Busch , LOGICAL GTSGRAL (auf CRAY wurde die GR-Farbe
C                           nicht umgesetzt, jetzt Kriterium
C                           auf CRAY wenn CGM WK 321 existiert)

C---- IPEN=1 WEISS
C---- IPEN=2 ROT
C---- IPEN=3 BLAU
C---- IPEN=4 GRUEN
C---- IPEN=5 PINK/VIOLETT
C---- IPEN=6 GELB
C---- IPEN=7 TUERKIS

      INTEGER FARBEN(8),INTLIN,ICOLOR,IFONT
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      LOGICAL FLGROT, GTSGRAL
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRREST/, /GRPP/
      EQUIVALENCE (PP(13),INTLIN),(PP(16),ICOLOR),(PP(17),SIZMRK)
      EQUIVALENCE (PP(9),IFONT),(PP(18),FLGROT)

      DATA FARBEN /1,2,4,3,6,5,7,8/

C---- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

      IF(NRLN.GT.1) THEN
        CALL GPL(NRLN,XLN,YLN)
        XCURR=XLN(NRLN)
        YCURR=YLN(NRLN)
        NRLN=0
      ENDIF
      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF

      GTSGRAL=.TRUE.
C     GKS laeuft und nur eine WK mit WKID 1 offen ist
C
      CALL GQEWK(1,IERR,NUMBER,IWKT)
      DO 40 I=1,NUMBER
      CALL GQEWK(I,IERR,NUMBER,IWKT)
c     testweise , da Aerger mit Versatec Routinen
      IF (IWKT.GT.16000) GOTO 40
      IF ( IWKT.EQ.4711 .and. ierr.eq.0 ) THEN
         goto  50
      ENDIF
c damit auch auf CRAY laeuft (CGM)
      IF ( IWKT.EQ.321.and. ierr.eq.0 ) THEN
         goto  50
      ENDIF
 40   CONTINUE
      GTSGRAL=.FALSE.
 50   CONTINUE


C---- Umsetzen der Farbe
      IF ( IPEN.GE.1.AND.IPEN.LE.8 ) THEN
C     Umsetzen nur fuer GRAL GKS
         IF (GTSGRAL) THEN
            ICOL=FARBEN(IPEN)
         ELSE
            if ( ipen.ne.8) then
               ICOL=IPEN
            else
               icol=1
            endif
         ENDIF
      ELSE
         ICOL = IPEN
      ENDIF

C---- SET TEXT COLOUR INDEX
      CALL GSTXCI(ICOL)

C---- SET POLYLINE COLOUR INDEX
      CALL GSPLCI(ICOL)

C---- SET POLYMARKER COLOUR INDEX
      CALL GSPMCI(ICOL)

C---- SET FILL AREA COLOUR INDEX
      CALL GSFACI(ICOL)
cc ---> COMMON  GRPP
      ICOLOR=IPEN
      END
