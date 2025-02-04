

      SUBROUTINE SIGUSR(IFIRST,JJJ,ZDS,DUMMY1,PSIG,DUMMY2,ARGST,
     .                  XD0,YD0,ZD0,XD1,YD1,ZD1)
C
C  INPUT:
C          IFIRST: FLAG FOR INITIALISATION
C          NCELL:   INDEX IN TALLY ARRAYS FOR CURRENT ZONE
C          JJJ:    INDEX OF SEGMENT ALONG CHORD
C          ZDS:    LENGTH OF SEGMENT NO. JJJ
C  OUTPUT: CONTRIB. FROM CELL NCELL AND CHORD SEGMENT JJJ TO:
C          THE H ALPHA FLUX PSIG(I),I=0,4 CONTRIBUTIONS
C          FROM ATOMS, MOLECULES, TEST IONS AND BULK IONS
C          THE INTEGRANT ARGST SUCH THAT INTEGR.(ARGST*DL) = PSIG
C
      USE PRECISION
      USE PARMMOD
      USE CESTIM
      USE COMUSR
      USE COMPRT
      USE CSPEI
      USE COMSOU
      USE CLOGAU
      USE COMXS
      USE COMSIG
      USE CUPD
      USE CZT1
      USE CTRCEI
      USE CGRID
      USE CCONA
      USE CADGEO
      USE CLGIN
      USE CSDVI
      USE COUTAU
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: JJJ
      INTEGER, INTENT(INOUT) :: IFIRST
      REAL(DP), INTENT(INOUT) :: PSIG(0:NSPZ),ARGST(0:NSPZ,NRAD)
      REAL(DP), INTENT(IN) :: ZDS,DUMMY1,DUMMY2,XD0,YD0,ZD0,XD1,YD1,ZD1
      RETURN
      END
