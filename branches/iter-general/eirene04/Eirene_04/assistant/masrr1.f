C
C*DK MASRR1
      SUBROUTINE MASRR1 (A,B,N,IS)
      USE PRECISION
      IMPLICIT NONE
C  IS IST DIE ANZAHL DER SPALTEN (MAXIMAL: 6)
      CHARACTER(11), INTENT(IN) :: A
      REAL(DP), INTENT(IN) :: B(*)
      INTEGER, INTENT(IN) :: N, IS
      INTEGER :: IE, I, J, IA
      WRITE (6,60) A
60    FORMAT (1X,A11)
      IA=1
      IE=MIN0(N,IA-1+IS)
      DO 62 J=1,N
      IF (IA.GT.IE) GOTO 63
      WRITE (6,61) (I,B(I),I=IA,IE)
      IA=IE+1
      IE=MIN0(N,IE+IS)
62    CONTINUE
63    CONTINUE
61    FORMAT(6(1X,I4,1X,1PE12.4))
      RETURN
      END