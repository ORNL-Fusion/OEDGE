 

      SUBROUTINE COFACT (A4,A3,I,J)
      USE PRECISION
      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: A4(4,4)
      REAL(DP), INTENT(OUT) :: A3(3,3)
      INTEGER, INTENT(IN) :: I, J
      INTEGER :: II, IS, JJ, JS

      IS = 0
      DO II=1,4
        IF (II /= I) IS = IS + 1
        JS = 0
        DO JJ=1,4
          IF (JJ /= J) JS = JS + 1
          IF ((II /= I) .AND. (JJ /= J)) A3(IS,JS) = A4(II,JJ)
        END DO
      END DO

      RETURN
      END
