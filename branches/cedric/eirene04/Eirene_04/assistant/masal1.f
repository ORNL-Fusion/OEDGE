C
C
C*DK MASAL1
      SUBROUTINE MASAL1 (A,NL,N)
      IMPLICIT NONE
      CHARACTER(6), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: NL(N)
      INTEGER :: I
      WRITE (6,60) A,(NL(I),I=1,N)
60    FORMAT (1X,A6/(20L3))
      RETURN
      END