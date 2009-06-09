c
c     ==================================================================
c     UPDTAED ROUTINES
c     ==================================================================
c
      REAL*8 FUNCTION X05BAF()
C     .. Scalar Arguments ..
      IMPLICIT NONE
      REAL*8 s

      CALL Clock(s)
      x05baf = s

      RETURN
      END
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
      SUBROUTINE TABPRC(i1,i2,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      WRITE(0,*) 'WARNING: Calling TABPRC routine substitute'
      RETURN
      END

      SUBROUTINE EXIT
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CALL GREND
      STOP
      END

      SUBROUTINE F04ARF(d1,i1,d2,i2,d3,d4,i3)
      IMPLICIT none
      INTEGER i1,i2,i3
      DOUBLE PRECISION d1(i1,*),d2(*),d3(*),d4(*)
      WRITE(0,*) 'WARNING: Calling F04ARF routine substitute'
      RETURN
      END

      SUBROUTINE F04ATF(d1,i1,d2,i2,d3,d4,i3,d5,d6,i5)
      IMPLICIT none
      INTEGER          IAA,i1,i2,i3,i4,i5
      PARAMETER       (IAA=50)
      DOUBLE PRECISION d1(i1,*),d2(*),d3(*),d4(IAA,IAA),d5(IAA),d6(IAA)
      WRITE(0,*) 'WARNING: Calling F04ATF routine substitute'
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02AJF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling X02AJF routine substitute'
      X02AJF = 1.0D+14
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02AKF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling X02AKF routine substitute'
      X02AKF = 1.0D+14
      RETURN
      END

      DOUBLE PRECISION FUNCTION X02ALF()
      IMPLICIT NONE
      WRITE(0,*) 'WARNING: Calling F04ALF routine substitute'
      X02ALF = 1.0D+14
      RETURN
      END



