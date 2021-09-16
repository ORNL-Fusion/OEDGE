C
C
      SUBROUTINE PRLLO (P1,P2,P3,P4,P5,IO,NF)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3), P5(3)
      INTEGER, INTENT(IN) :: IO
      LOGICAL, INTENT(IN) :: NF
      REAL(DP) :: CORD(15)
      INTEGER :: I

      DO 1 I=1,3
      CORD(I)=P1(I)
      CORD(3+I)=P2(I)
      CORD(6+I)=P4(I)
      CORD(9+I)=P5(I)
1     CORD(12+I)=P3(I)
      CALL PL3Q (CORD,5,IO,NF)
      RETURN
      END
