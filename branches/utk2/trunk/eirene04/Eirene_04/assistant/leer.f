C
C
C*DK LEER
      SUBROUTINE LEER (N)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      INTEGER :: I
      DO 1 I=1,N
         WRITE (6,60)
1        CONTINUE
60    FORMAT ('      ')
      RETURN
      END
