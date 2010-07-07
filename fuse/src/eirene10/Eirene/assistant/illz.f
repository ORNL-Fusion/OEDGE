C
c
      FUNCTION ILLZ(N,V,IV)
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: V(*)
      INTEGER, INTENT(IN) :: N, IV
      INTEGER :: J, ILLZ
      J=1
      IF (IV.LT.0) J=1-(N-1)*IV
      DO 17,ILLZ=1,N
CODER    IF (V(J).LT.0) GOTO 18    BEI REAL ODER INTEGER
CPB      IF (.NOT.V(J)) GOTO 18
         IF (V(J)) GOTO 18
   17    J=J+IV
   18 ILLZ=ILLZ-1
      END
