C
C
C*DK MASRR4
      SUBROUTINE MASRR4 (A,B,C,D,E,N)
      USE PRECISION
      IMPLICIT NONE
      CHARACTER(22), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      REAL(DP), INTENT(IN) :: B(N), C(N), D(N), E(N)
      INTEGER :: J
      WRITE (6,60) A
60    FORMAT (1X,A22)
      DO 5 J=1,N
5        WRITE (6,61) J,B(J),C(J),D(J),E(J)
61    FORMAT (1X,I4,1X,4(1PE12.4,3X))
      RETURN
      END