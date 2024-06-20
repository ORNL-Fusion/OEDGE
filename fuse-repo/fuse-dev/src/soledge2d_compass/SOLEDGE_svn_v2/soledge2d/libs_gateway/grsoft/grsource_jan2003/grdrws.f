C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRDRWS(X,Y,NR)
C----------------------------------------------------------------------
C     AUTHOR: GROTEN
C     UPDATE: M.BUSCH 11.3.91
C     Anpassung an GTSGRAL GKS
C     es gibt 5 Standardmarker 1,2,3,4,5
c     und 14 implementationsabh. Marker -101, bis -114
C     Das GKS der TU Berlin hatte die Marker 1 bis 53, wobei 1-5 die
C     Standardmarker sind
C     ---> 1,2,3,4,5 wie frueher
C     ---> 6-53 auf -101 bis -114 verteilen
C               -MOD(N-5,14) - 100
C     UPDATE: G.Groten 18.11.91 : filled marker GRFLLS
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C UPDATE 11. 3.1991 Busch , Marker auf Vektor XMR, YMR sammeln
C UPDATE  1. 6.1992 Busch , s. Notiz zum Datum
C UPDATE 14. 1.1993 Busch , GRJMPS und nachfolgend GRDRWS mit zwei
C                   verschiedenen Symbolen lieferte stattdessen das 2.
C                   Symbol doppelt
C ---------------------------------------------------------------------
      REAL      X,Y
      INTEGER   NR,NS

      REAL   XH(2),YH(2),XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR

      SAVE /GRREST/, /MKRTP/
C---- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

      IF (NRLN.EQ.MAXPKT ) THEN
         CALL GPL(NRLN,XLN,YLN)
         XCURR=XLN(NRLN)
         YCURR=YLN(NRLN)
         NRLN=0
      ENDIF
      IF (NRMR.EQ.MAXPKT) THEN
         CALL GPM(NRMR,XMR,YMR)
         NRMR=0
      ENDIF

C---- LINIE ZIEHEN VON CURR. POSITION BIS X#Y
c 1. 6.92 jetzt wieder Polyline fuer 2 Werte (hier kein Sammeln!)
      XH(1)=XCURR
      XH(2)=X
      YH(1)=YCURR
      YH(2)=Y
      CALL GPL(2,XH,YH)
c.    NRLN=NRLN+1
c.    XLN(NRLN)=X
c.    YLN(NRLN)=Y

      NS=NR
      IF ( NS.LT.100) THEN
         IF (NS.GT.5) THEN
            NS = - MOD(NS-5,14) -100
            IF (NS .EQ.-100) NS=-114
         ENDIF

         IF (NS.NE.MKTYP) THEN
c           IF (NRMR.GT.1)  CALL GRPTS (XMR,YMR,NRMR,MKTYP)
C 14.1.93
            IF (NRMR.GT.0)  CALL GRPTS (XMR,YMR,NRMR,MKTYP)
            MKTYP=NS
C---------- SET MARKER TYPE
            CALL GSMK(NS)
         ENDIF
c        CALL GPM(1,X,Y)
         NRMR=NRMR+1
         XMR(NRMR)=X
         YMR(NRMR)=Y
      ELSE
         CALL GRFLLS(1,X,Y,NS)
      ENDIF
      XCURR=X
      YCURR=Y
      END
