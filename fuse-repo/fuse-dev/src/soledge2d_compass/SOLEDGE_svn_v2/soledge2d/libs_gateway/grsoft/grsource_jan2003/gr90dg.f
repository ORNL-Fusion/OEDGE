CC@PROCESS OPT(3) NOSDUMP NOGOSTMT FIPS(F)
C***********************************************************************
C
C     Copyright 1991 by Forschungszentrum (KFA) Juelich GMBH
C
C***********************************************************************
c-----------------------------------------------------------------------
C     UPDATE:10.10. 90    M.BUSCH
C                         CRAY Version setzt nur das Flag
c-----------------------------------------------------------------------
      SUBROUTINE GR90DG

      LOGICAL FLGROT
      COMMON /GRPP/ PP(17),FLGROT
CDEC$ PSECT /GRPP/ NOSHR
      SAVE /GRPP/

      REAL MOUT(2,3)

      FLGROT=.TRUE.
CJHeinen
      CALL GEVTM(0.5,0.5,0.,0.,3.141592*0.5,1.,1.,1,MOUT)
      CALL GSSGT(1,MOUT)
CJHeinen

      GOTO 99

CCC   LIEGENDES RECHTECK
CCC
      ENTRY GR00DG

      FLGROT=.FALSE.
CJHeinen
      CALL GEVTM(0.5,0.5,0.,0.,0.0,1.,1.,1,MOUT)
      CALL GSSGT(1,MOUT)
CJHeinen

 99   END
