C@PROCESS OPT(3) NOSDUMP NOGOSTMT
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
C     SCHLIESSEN DES AKTUELLEN SEGMENTES UND EROEFFNEN EINES NEUEN
C     SEGMENTES
C     AUTHOR: M. BUSCH (ZDV131@ZAM001.KFA-JUELICH.DE)
C     DATUM:  12. 2.1992
      SUBROUTINE GRNWSG
      CALL GQOPSG(IERR,NAME)
      IF (IERR.EQ.0) THEN
         CALL GCLSG
         NAME=NAME+1
         CALL GCRSG(NAME)
      ENDIF
      END
