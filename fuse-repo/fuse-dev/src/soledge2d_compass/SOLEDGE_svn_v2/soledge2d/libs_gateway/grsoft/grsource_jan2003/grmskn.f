C@PROCESS OPT(3) NOSDUMP NOGOSTMT IL(DIM) FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C     AUTHOR:  GROTEN
C     UPDATE:  28. 2.1991 BUSCH
C UPDATE 26. 2.1991 Busch , COMMON GRREST neu
C                   GPM durch GRPTS ersetzt
C UPDATE 14. 6.1993 Busch , GCLSG nur wenn Segment offen
C UPDATE 21.03.1994 BUsch, mehrere SEgmente nicht moeglich bei GLIGKS
      SUBROUTINE GRMSKN
C
      REAL    XLN(300),YLN(300),XMR(300),YMR(300)
      CHARACTER GKSTYP*7
      COMMON /GRREST/MAXPKT,NRLN,XCURR,YCURR,XLN,YLN,
     $               NRMR,XMR,YMR
CDEC$ PSECT /GRREST/ NOSHR
      COMMON /GRMSK/ MOB,NAME1,NAME2
CDEC$ PSECT /GRMSK/ NOSHR
      COMMON /GRGKS/ GKSTYP
CDEC$ PSECT /GRGKS/ NOSHR
      SAVE /GRREST/,/GRMSK/ ,/GRGKS/

      IF ( GKSTYP.EQ.'GLIGKS') GOTO 999
C
      IF(NRLN.GT.1) THEN
        CALL GPL(NRLN,XLN,YLN)
        NRLN=0
      ENDIF
      IF(NRMR.GT.0) THEN
        CALL GPM(NRMR,XMR,YMR)
        NRMR=0
      ENDIF
      CALL GQOPSG(IERR,NAME)
      if (ierr.eq.0) CALL GCLSG
      NAME=NAME+1
      CALL GCRSG(NAME)
      NAME1=NAME
999   END
