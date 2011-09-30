C
C
C*DK MASYR1
      SUBROUTINE MASYR1 (A,B,NL,I1,N0,N1,M0,M1,TX)
      USE PRECISION
      USE COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(9), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I1, N0, N1, M0, M1
      REAL(DP), INTENT(IN) :: B(N0:N1,M0:M1)
      LOGICAL, INTENT(IN) :: NL(N0:N1,M0:M1)
      CHARACTER(8), INTENT(IN) :: TX(*)
      CHARACTER(8) :: KK(100)
      INTEGER :: KKK(100)
      INTEGER :: J, KE, K
      KE=0
      DO 10 J=1,N1
        IF (.NOT.NL(J,I1)) GOTO 10
        KE=KE+1
        KK(KE)=TX(J)
        KKK(KE)=J
10    CONTINUE
      IF (KE.EQ.0) RETURN
      WRITE (iunout,60) A
60    FORMAT (1X,A9)
5        WRITE (iunout,61) (KK(K),B(KKK(K),I1),K=1,KE)
61    FORMAT (3(1X,1A8,1PE12.4,2X))
      RETURN
      END
