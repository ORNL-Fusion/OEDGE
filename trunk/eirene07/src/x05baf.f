c
c     ==================================================================
c     SYSTEM CALLS
c     ==================================================================
c
      SUBROUTINE Clock(s)
      IMPLICIT none
      DOUBLE PRECISION s
      REAL etime
      REAL rdum1,val(2)
      rdum1 = etime(val)
      s = DBLE(val(1))
      RETURN
      END

      SUBROUTINE Clock_1970(s)
      IMPLICIT none
      DOUBLE PRECISION s
      INTEGER time
      s = DBLE(time())
      RETURN
      END
c
c     ==================================================================
c     SUBSITUTE SUBROUTINES
c     ==================================================================
c
      DOUBLE PRECISION FUNCTION CPU_TIME()
c      DOUBLE PRECISION FUNCTION X05BAF()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION s
      CALL Clock(s)
      CPU_TIME = s
      RETURN
      END

