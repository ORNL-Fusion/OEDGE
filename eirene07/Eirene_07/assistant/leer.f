C
C
C*DK LEER
      SUBROUTINE LEER (N)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      INTEGER :: I
      DO 1 I=1,N
         WRITE (iunout,60)
1        CONTINUE
60    FORMAT ('      ')
      RETURN
      END
