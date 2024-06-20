C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
      SUBROUTINE GRPTS(X,Y,N,NR)
C----------------------------------------------------------------------
C     G.Groten 12.9.91
C     Neuaufnahme, um zu GRLN und GRCHN das Analogon zu haben.
C
C     X,Y: Vektoren mit N Punkten, an denen Marker mit der Nummer NR
C          gezeichnet werden sollen.
C
C     Es gibt 5 Standardmarker 1,2,3,4,5
c     und 14 implementationsabh. Marker -101, bis -114
C     Das GKS der TU Berlin hatte die Marker 1 bis 53, wobei 1-5 die
C     Standardmarker sind
C     ---> 1,2,3,4,5 wie frueher
C     ---> 6-53 auf -101 bis -114 verteilen
C               -MOD(N-5,14) - 100
C----------------------------------------------------------------------
C UPDATE 20. 2.1991 Busch , GKS Treiber erlaubt max. 300 Marker pro
C                   Aufruf
C----------------------------------------------------------------------
      INTEGER   N,NR,NS
      REAL      X(*),Y(*)

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /MKRTP/ MKTYP
CDEC$ PSECT /MKRTP/ NOSHR
      SAVE /GRREST/,/MKRTP/
      IF (N.GT.0) THEN
C------- AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN

         IF (NRLN.GT.1) THEN
            CALL GPL(NRLN,XLN,YLN)
            XCURR=XLN(NRLN)
            YCURR=YLN(NRLN)
            NRLN=0
         ENDIF

C------- AUSGABE EVT. VORHANDENER GRJMPS,.. -DATEN

         IF (NRMR.GT.0) THEN
            CALL GPM(NRMR,XMR,YMR)
            NRMR=0
         ENDIF

         NS = NR
         IF (NS.LT.100) THEN
            IF (NS.GT.5 .AND. NS.LT.100 ) THEN
               NS = - MOD(NS-5,14) -100
               IF (NS .EQ.-100) NS=-114
            ENDIF

            IF (NS.NE.MKTYP) THEN
               MKTYP=NS
C------------- SET MARKER TYPE
               CALL GSMK(NS)
            ENDIF

           DO 1 I=1,N,MAXPKT
              CALL GPM(MIN(MAXPKT,N-I+1),X(I),Y(I))
  1        CONTINUE
         ELSE
            CALL GRFLLS(N,X,Y,NS)
         ENDIF
         XCURR=X(N)
         YCURR=Y(N)
      ENDIF

      END
