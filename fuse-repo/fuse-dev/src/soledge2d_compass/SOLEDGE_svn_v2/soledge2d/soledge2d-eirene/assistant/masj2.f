C
C*DK MASJ2
      SUBROUTINE EIRENE_MASJ2 (A,I,J)
      USE EIRMOD_PRECISION
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      IMPLICIT NONE
      CHARACTER(16), INTENT(IN) :: A
      INTEGER, INTENT(IN) :: I, J
      REAL(DP) :: FJ
      IF (J.LT.1000000) THEN
        WRITE (iunout,60) A,I,J
      ELSEIF (J.LT.100000000) THEN
        WRITE (iunout,70) A,I,J
      ELSE
        fj=j
        WRITE (iunout,80) A,I,fj
      ENDIF
60    FORMAT (1X,A16,2(3X,I6))
70    FORMAT (1X,A16,(3X,I6,3X,I8))
80    FORMAT (1X,A16,(3X,I6,3X,1pe12.4))
      RETURN
      END
