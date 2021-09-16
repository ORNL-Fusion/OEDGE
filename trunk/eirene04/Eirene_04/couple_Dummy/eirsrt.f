
C
      SUBROUTINE EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: fluxes(*)
      REAL(DP), INTENT(IN) :: DELTAT, B2BRM, B2RD, B2Q, B2VP
      logical, INTENT(IN) :: lstop, ltime
C
      return
      end
