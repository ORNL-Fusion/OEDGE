C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRJMPS(X,Y,NR)
C----------------------------------------------------------------------
c     M.Busch 11.3.91
C     Anpassung an GTSGRAL GKS
C     es gibt 5 Standardmarker 1,2,3,4,5
c     und 14 implementationsabh. Marker -101, bis -114
C     Das GKS der TU Berlin hatte die Marker 1 bis 53, wobei 1-5 die
C     Standardmarker sind
C     ---> 1,2,3,4,5 wie frueher
C     ---> 6-53 auf -101 bis -114 verteilen
C               -MOD(N-5,14) - 100
C Update : 12.9.91 Groten
C Update : 18.11.91 Groten : Fillfonts (GRFLLS)
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C----------------------------------------------------------------------
      REAL      X,Y
      INTEGER   NR,NS

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR

      SAVE /GRREST/, /MKRTP/
C---- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

      IF (NRLN.GT.1) THEN
         CALL GPL(NRLN,XLN,YLN)
         XCURR=XLN(NRLN)
         YCURR=YLN(NRLN)
         NRLN=0
      ENDIF
      IF (NRMR.EQ.MAXPKT) THEN
         CALL GRPTS(XMR,YMR,NRMR,NS)
         NRMR=0
      ENDIF

C---- Update Groten 12.9.91 : NR darf nicht ueberschrieben werden
      NS = NR
      IF (NS.LT.100) THEN
         IF (NS.GT.5) THEN
            NS = - MOD(NS-5,14) -100
            IF (NS .EQ.-100) NS=-114
         ENDIF

         IF (NS.NE.MKTYP) THEN
c           wenn sich der Markertyp aendert, muessen gesammelte
C           vorher Marker geplottet werden
            IF ( NRMR.GT.0) CALL GRPTS(XMR,YMR,NRMR,NS)
            NRMR=0

            MKTYP=NS
C---------- SET MARKER TYPE
            CALL GSMK(NS)

         ENDIF

         NRMR=NRMR+1
         XMR(NRMR)=X
         YMR(NRMR)=Y


      ELSE
         CALL GRFLLS(1,X,Y,NS)
      ENDIF
      XCURR=X
      YCURR=Y

      END
