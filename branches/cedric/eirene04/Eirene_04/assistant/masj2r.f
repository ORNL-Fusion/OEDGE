C
C*DK MASJ2R
      SUBROUTINE MASJ2R (A,I,J,B)
      USE PRECISION
      IMPLICIT NONE
      CHARACTER(24), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B
      INTEGER, INTENT(IN) :: I, J
      WRITE (6,60) A,I,J,B
60    FORMAT (1X,A24,3X,I6,3X,I6,3X,1PE12.4)
      RETURN
      END