      SUBROUTINE Clock(s)
      IMPLICIT none

      DOUBLE PRECISION s

      REAL rdum1

      CALL CLOCK@(rdum1)

      s = DBLE(rdum1)

      RETURN
      END

      SUBROUTINE Clock_1970(s)
      IMPLICIT none

      DOUBLE PRECISION s

      INTEGER time

      CALL SECONDS_SINCE_1980@(s)

      RETURN
      END

