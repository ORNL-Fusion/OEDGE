C
C
      SUBROUTINE SIGTST(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST)
C
C  ONLY FOR TESTING LINE INTEGRAL ROUTINES
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IFIRST, JJJ
      REAL(DP), INTENT(IN) :: ZDS, DUMMY1, DUMMY2
      REAL(DP), INTENT(OUT) :: PSIG(0:NSPZ+2), ARGST(0:NSPZ+2,NRAD)
      INTEGER :: ISP, ICELL
      SAVE
C
      DO 100 ISP=0,NSPZ+2
        PSIG(ISP)=0.
        DO 100 ICELL=1,NRAD
          ARGST(ISP,ICELL)=0.
100   CONTINUE
C
      WRITE (6,*) 'FROM SIGTST ',JJJ,ZDS,NRCELL,NPCELL,NTCELL,NACELL
      RETURN
      END