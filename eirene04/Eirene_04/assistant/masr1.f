C
C
C*DK MASR1
      SUBROUTINE MASR1 (A,B)
      USE PRECISION
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B
      WRITE (6,60) A,B
60    FORMAT (1X,A8,3X,1PE12.4)
      RETURN
      END
