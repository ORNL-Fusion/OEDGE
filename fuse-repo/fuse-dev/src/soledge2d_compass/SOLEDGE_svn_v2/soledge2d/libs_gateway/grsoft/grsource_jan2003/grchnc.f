C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
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
C Update: 12. 9. 91 Groten
C UPDATE 20. 2.1991 Busch , GKS Treiber erlaubt max. 300 Werte pro Kurve
C                   GPL Aufruf durch GRLN Aufruf ersetzt
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C UPDATE 11. 3.1992 Busch , Abfrage KFA Symbole > 100
C----------------------------------------------------------------------
      SUBROUTINE GRCHNC(XX,YY,M,NR)
      REAL      XX(*),YY(*)
      INTEGER   M,NR,NS

      REAL      XXH(2),YYH(2)
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR

      SAVE /GRREST/, /MKRTP/
C---- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

      IF (M.LT.1) THEN
Ctest-   WRITE(*,*) 'GRCHNC: less than 1 point'
      ELSE
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

C------- VERBINDUNG ZWISCHEN GRCHN-ENDE | CURR. POS.
C------- MIT GRCHNC-ANFANG

         XXH(1)=XCURR
         XXH(2)=XX(1)
         YYH(1)=YCURR
         YYH(2)=YY(1)
         CALL GPL(2,XXH,YYH)
         IF (M.GT.1) CALL GRLN(XX,YY,M)

C------- Update 12.9.91 Groten:NR nicht ueberschreiben; oft Konstante!
         NS = NR
         IF (NS.LT.100) THEN
            IF (NS.GT.5) THEN
               NS = - MOD(NS-5,14) -100
               IF (NS .EQ.-100) NS=-114
            ENDIF

            IF (NS.NE.MKTYP) THEN
               MKTYP=NS
C---------- SET MARKER TYPE
               CALL GSMK(NS)
            ENDIF
            CALL GRPTS(XX,YY,M,NS)
         ELSE
            CALL GRFLLS(M,XX,YY,NS)
         ENDIF

         XCURR=XX(M)
         YCURR=YY(M)
      ENDIF
      END
