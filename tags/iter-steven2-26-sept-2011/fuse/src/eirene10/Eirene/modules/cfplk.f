      MODULE CFPLK
C   parameters for fokker planck collision operator
      USE PRECISION

      IMPLICIT NONE

      PUBLIC 

      REAL(DP), PUBLIC, SAVE ::
     R E0PAR, VELPAR, VELPER, VLXPAR, VLYPAR, VLZPAR, SIGPAR, TAUE,
     R BVEC(3), BBX, BBY, BBZ
      LOGICAL, public, save :: LCART,COLFLAG

      CONTAINS


      SUBROUTINE DEALLOC_CFPLK

      RETURN
      END SUBROUTINE DEALLOC_CFPLK


      END MODULE CFPLK
