C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C UPDATE 8. 8. 1990 GROTEN
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
      SUBROUTINE GRDRW(X,Y)
      REAL   X,Y

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      SAVE  /GRREST/

C---- AUSGABE VORHANDENER GRJMP- UND GRDRW-DATEN,
C---- WENN 300 ZUSAMMEN SIND

      IF (NRLN.EQ.maxpkt) THEN
         CALL GPL(NRLN,XLN,YLN)
         XLN(1)=XLN(maxpkt)
         YLN(1)=YLN(maxpkt)
         NRLN=1
      ENDIF

      IF (NRMR.GT.0) THEN
         CALL GPM(NRMR,XMR,YMR)
         NRMR=0
      ENDIF


      NRLN=NRLN+1
      XLN(NRLN)=X
      YLN(NRLN)=Y
      XCURR=X
      YCURR=Y
      END
