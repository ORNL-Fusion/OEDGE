C
C
      SUBROUTINE EIRENE_PRLLO (P1,P2,P3,P4,P5,IO,NF)
 
      USE EIRMOD_PRECISION
 
      IMPLICIT NONE
c Yannick 
c      REAL(DP), INTENT(IN) :: P1(3), P2(3), P3(3), P4(3), P5(3)
      REAL(DP), INTENT(IN) :: P1, P2, P3, P4, P5
      INTEGER, INTENT(IN) :: IO
      LOGICAL, INTENT(IN) :: NF
      REAL(DP) :: CORD(15)
      INTEGER :: I
c Yannick
c      DO 1 I=1,3
      I=1
      CORD(1)=P1
      CORD(4)=P2
      CORD(7)=P4
      CORD(10)=P5
1     CORD(13)=P3
      CALL EIRENE_PL3Q (CORD,5,IO,NF)
      RETURN
      END
