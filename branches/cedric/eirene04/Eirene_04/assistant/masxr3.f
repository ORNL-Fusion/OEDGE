C
C
C*DK MASXR3
      SUBROUTINE MASXR3 (A,B,C,D,KI,KE,KST,N,M,TX)
      USE PRECISION
      IMPLICIT NONE
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(14), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: KI, KE, KST, N, M
      REAL(DP), INTENT(IN) :: B(M),C(N,M), D(N,M)
      INTEGER :: J, K
      IF (KE.LT.KI) RETURN
      WRITE (6,60) A
60    FORMAT (1X,A14)
      DO 5 J=1,M
5        WRITE (6,61) B(J),(TX(K),C(K,J),D(K,J),K=KI,KE,KST)
61    FORMAT (1X,1PE12.4,2X,6(A8,1X,2(1PE12.4,2X)))
      RETURN
      END