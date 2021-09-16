C
C
C*DK HEADNG
      SUBROUTINE HEADNG (A,N)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: N
      INTEGER :: I
      CHARACTER(1) :: U(72)
      DATA U/72*'='/
      WRITE (6,'(1X,A)') A
      WRITE (6,60) (U(I),I=1,N)
60    FORMAT (1X,72A1)
      RETURN
      END
