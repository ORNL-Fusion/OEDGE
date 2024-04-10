C@process opt(3) nosdump nogostmt
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C  UPDATE    19.12.1990 BUSCH
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
      SUBROUTINE GRDSH(A1,A2,A3)
      REAL   A1,A2,A3

      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      INTEGER  LTYPE
      COMMON   /GRPP/ PP(18)
CDEC$ PSECT /GRPP/ NOSHR
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      SAVE /GRPP/, /GRREST/

      PP(10)=A1
      PP(11)=A2
      PP(12)=A3

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

      IF (PP(11).LT.0.0035) THEN
C------- 'SOLID' LINE
         LTYPE=1
      ELSE IF(PP(12).LT.PP(10)) THEN
C------ 'DASHED-DOTTED' LINE
         LTYPE=4
      ELSE IF(ABS(PP(10)-PP(12)).LT.0.1E-2.AND.PP(10).LT.PP(11)) THEN
C------- 'DOTTED' LINE
          LTYPE=3
      ELSE
C------- 'DASHED' LINE
          LTYPE=2
      ENDIF

C---- SET LINETYPE

      CALL GSLN(LTYPE)

      END
