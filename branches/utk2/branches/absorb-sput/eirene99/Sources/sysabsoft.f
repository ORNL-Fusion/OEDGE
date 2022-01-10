      SUBROUTINE Clock(s)
      IMPLICIT none

      DOUBLE PRECISION s

      REAL etime

      REAL val(2)

      s = DBLE(etime(val))

      RETURN
      END



      SUBROUTINE Clock_1970(s)
      IMPLICIT none

      DOUBLE PRECISION s

      INTEGER time

      s = DBLE(time())

      RETURN
      END



