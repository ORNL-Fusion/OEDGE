

      SUBROUTINE DEVCSF(N,A,LDA,EVAL,EVEC,LDEVEC)
      USE PRECISION
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N, LDA, LDEVEC
      REAL(DP), INTENT(IN) :: A(LDA,N)
      REAL(DP), INTENT(OUT) :: EVAL(N), EVEC(LDEVEC,N)
      INTEGER :: IERR, J, I, IAA
      REAL(DP), ALLOCATABLE :: EVEC1(:,:),EVAL1(:),
     .                         FV1(:),FV2(:)

      entry evcsf (n,a,lda,eval,evec,ldevec)

      ALLOCATE (EVEC1(LDA,LDA))
      ALLOCATE (EVAL1(LDA))
      ALLOCATE (FV1(LDA))
      ALLOCATE (FV2(LDA))

      CALL RS(LDA,N,A,EVAL1,1,EVEC1,FV1,FV2,IERR)
      DO 10,I=1,N
         EVAL(I) = EVAL1(N-I+1)
         EVEC(1:LDEVEC,I) = EVEC1(1:LDEVEC,N-I+1)
10    CONTINUE

      DEALLOCATE (EVEC1)
      DEALLOCATE (EVAL1)
      DEALLOCATE (FV1)
      DEALLOCATE (FV2)
      END
