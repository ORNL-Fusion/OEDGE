C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM) FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C  UPDATE 25. 9.1990 BUSCH
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
      SUBROUTINE GRJMP(X,Y)
      REAL    X,Y
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
CCC
      NRLN=1
      XLN(NRLN)=X
      YLN(NRLN)=Y
      XCURR=X
      YCURR=Y
      END
