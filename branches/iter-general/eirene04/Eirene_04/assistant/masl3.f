C
C
C*DK MASL3
      SUBROUTINE MASL3 (A,B,C,D)
      IMPLICIT NONE
      CHARACTER(*), INTENT(IN) :: A
      LOGICAL, INTENT(IN) :: B, C, D 
      WRITE (6,60) A,B,C,D
60    FORMAT (1X,A,5X,3L3)
      RETURN
      END