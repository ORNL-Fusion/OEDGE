

      SUBROUTINE BUBBLE (IAR,N)

      USE PRECISION
      IMPLICIT NONE

      INTEGER, INTENT(INOUT) :: IAR(N)
      INTEGER, INTENT(IN) :: N

      INTEGER IH, K, J, IBOUND

      IBOUND = N

      DO WHILE (IBOUND > 1)
        K = 1
        J = 1
        DO WHILE (IBOUND > J)
          IF (IAR(J) > IAR(J+1)) THEN
            IH = IAR(J)
            IAR(J) = IAR(J+1)
            IAR(J+1) = IH
            K = J
          END IF
          J = J + 1
        END DO
        IBOUND = K
      END DO

      RETURN
      END 
