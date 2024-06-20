C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM) FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C  UPDATE: 25. 9.1990 BUSCH
C UPDATE 20. 2.1991 Busch , GKS Treiber erlaubt max. 300 Werte pro Kurve
C                   GPL Aufruf durch GRLN Aufruf ersetzt
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
      SUBROUTINE GRLNCN(XX,YY,M)
      INTEGER   M
      REAL      XXH(2),YYH(2),XX(1),YY(1)
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      SAVE /GRREST/
CCC   AUSGABE EVT. VORHANDENER GRJMP- UND GRDRW-DATEN
CCC
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
CCC   LINESTYLE, FARBE UND STRICHDICKE NOCH NICHT
CCC   IMPLEMENTIERT
CCC
      XXH(1)=XCURR
      XXH(2)=XX(1)
      YYH(1)=YCURR
      YYH(2)=YY(1)
      CALL GPL(2,XXH,YYH)
      CALL GRLN(XX,YY,M)
      XCURR=XX(M)
      YCURR=YY(M)
      END
