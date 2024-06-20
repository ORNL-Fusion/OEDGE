C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRSPTS(ISPOTS)
      INTEGER ISPOTS
C-----------------------------------------------------------------------
C     Update:  3. 8. 90  M. Busch
C              ISPOTS <= 18     Linewidth 1
C          >18 ISPOTS <= 20     Linewidth 2
C          >20 ISPOTS <= 22     Linewidth 3
C          >22 ISPOTS <= 24     Linewidth 4
C          >24 ISPOTS <= 26     Linewidth 5
C          >26 ISPOTS <= 28     Linewidth 6
C          >28 ISPOTS <= 30     Linewidth 7
C          >30 ISPOTS <= 32     Linewidth 8
C          >32 ISPOTS           Linewidth 8
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C-----------------------------------------------------------------------
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
C---- COMMONBLOCK DER STANDARDWERTE BZW. DER GEAENDERTEN TABELLENWERTE
      COMMON /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      INTEGER INTLIN
      EQUIVALENCE (PP(13),INTLIN)
      SAVE /GRREST/, /GRPP/

C---- AUSGABE EVT. VORHANDENER GRJMP- + GRDRW-DATEN

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

      INTLIN=ISPOTS

      IF (ISPOTS.LE.18) THEN
         WIDTH=1.
         ELSEIF (ISPOTS.LE.20) THEN
            WIDTH=2.
            ELSEIF (ISPOTS.LE.22) THEN
               WIDTH=3.
               ELSEIF (ISPOTS.LE.24) THEN
                  WIDTH=4.
                  ELSEIF (ISPOTS.LE.26) THEN
                     WIDTH=5.
                     ELSEIF (ISPOTS.LE.28) THEN
                        WIDTH=6.
                        ELSEIF (ISPOTS.LE.30) THEN
                           WIDTH=7.
                           ELSEIF (ISPOTS.LE.32) THEN
                           WIDTH=8.
                              ELSE
                              WIDTH=8.
      ENDIF
      CALL GSLWSC(WIDTH)

      END
