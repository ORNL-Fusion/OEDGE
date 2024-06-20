C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE:   2. 4. 1991 BUSCH
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C BUSCH 8.6.93 FILL AREA STYLE HOLLOW
      SUBROUTINE GRFILL(N,XX,YY,ISTYLE,ITYPE)
      INTEGER N,ISTYLE,ITYPE
      REAL    XX(N),YY(N)

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR

      SAVE /GRREST/
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

      CALL GQLN ( IERR, IPLTYP)
      ISTYL=ISTYLE
      IF (ISTYL.EQ.2) THEN
         ISTYL=0
         WRITE(6,*) 'Pattern fuer Fill Area nicht implementiert.'
      ELSE IF(ISTYL.GT.3.OR.ISTYL.LT.0) THEN
         WRITE(6,*) 'FALSCHER 3. PARAMETER IN GRFILL.'
         WRITE(6,*) 'Moeglich ist: ISTYLE=0,1,3'
         ISTYL=0
      ENDIF

C---- SET FILL AREA INTERIOR STYLE
      CALL GSFAIS(ISTYL)

C---- SET FILL AREA STYLE INDEX
      IF(ISTYL.EQ.3) CALL GSFASI(ITYPE)

C---- FILL AREA
      CALL GFA(N,XX,YY)

      XCURR=XX(N)
      YCURR=YY(N)
C BUSCH 8.6.93 FILL AREA STYLE HOLLOW
      CALL GSFAIS(0)

      END
