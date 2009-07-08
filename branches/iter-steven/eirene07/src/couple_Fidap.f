C EIRENE07 COMPILATION
C ===== SOURCE: bubble.f


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
C ===== SOURCE: df2triqua.f
      subroutine df2triqua (vector, r, s, t)

      USE PRECISION

      implicit none
      real(dp), intent(out) :: vector(3,3,27)
      real(dp), intent(in) :: r, s, t
      real(dp) :: halb, viertel, achtel, emr, ems, emt, epr, eps, ept,
     .            emrq, emsq, emtq, em2r, em2s, em2t, ep2r, ep2s, ep2t,
     .            rmrq, smsq, tmtq, rprq, spsq, tptq

      halb = 0.5_dp
      viertel = 0.25_dp
      achtel = 0.125_dp

      emr = 1._dp - r
      ems = 1._dp - s
      emt = 1._dp - t

      em2r = 1._dp - 2._dp * r
      em2s = 1._dp - 2._dp * s
      em2t = 1._dp - 2._dp * t

      epr = 1._dp + r
      eps = 1._dp + s
      ept = 1._dp + t

      ep2r = 1._dp + 2._dp * r
      ep2s = 1._dp + 2._dp * s
      ep2t = 1._dp + 2._dp * t

      emrq = 1._dp - r*r
      emsq = 1._dp - s*s
      emtq = 1._dp - t*t

      rmrq = r - r*r
      smsq = s - s*s
      tmtq = t - t*t

      rprq = r + r*r
      spsq = s + s*s
      tptq = t + t*t


!     vector(k,j,i) = df(Ni)/dkdj, 1<i<27, k,j aus {r,s,t}


!      vector(1,1) = -achtel * s * t * em2r * ems * emt
!                  = -achtel * em2r * smsq * tmtq
      vector(1,1,1) = -achtel * (-2._dp) * smsq* tmtq
      vector(1,2,1) = -achtel * em2r * em2s * tmtq
      vector(1,3,1) = -achtel * em2r * smsq * em2t

!      vector(2,1) = -achtel * r * t * emr * em2s * emt
!                  = -achtel * rmrq * em2s * tmtq
      vector(2,1,1) = -achtel * em2r * em2s * tmtq
      vector(2,2,1) = -achtel * rmrq * (-2._dp) * tmtq
      vector(2,3,1) = -achtel * rmrq * em2s * em2t

!      vector(3,1) = -achtel * r * s * emr * ems * em2t
!                  = -achtel * rmrq * smsq * em2t
      vector(3,1,1) = -achtel * em2r * smsq * em2t
      vector(3,2,1) = -achtel * rmrq * em2s * em2t
      vector(3,3,1) = -achtel * rmrq * smsq * (-2._dp)




!      vector(1,2) = viertel * s * t * (-2._dp) * r * ems * emt
!                  = viertel * (-2._dp) * r * smsq * tmtq
      vector(1,1,2) = viertel * (-2._dp) * smsq * tmtq
      vector(1,2,2) = viertel * (-2._dp) * r * em2s * tmtq
      vector(1,3,2) = viertel * (-2._dp) * r * smsq * em2t

!      vector(2,2) = viertel * t * emrq * em2s * emt
!                  = viertel * emrq * em2s * tmtq
      vector(2,1,2) = viertel * (-2._dp) * r * em2s * tmtq
      vector(2,2,2) = viertel * emrq * (-2._dp) * tmtq
      vector(2,3,2) = viertel * emrq * em2s * em2t

!      vector(3,2) = viertel * s * emrq * ems * em2t
!                  = viertel * emrq * smsq * em2t
      vector(3,1,2) = viertel * (-2._dp) * r * smsq * em2t
      vector(3,2,2) = viertel * emrq * em2s * em2t
      vector(3,3,2) = viertel * emrq * smsq * (-2._dp)




!      vector(1,3) = achtel * s * t * ep2r * ems * emt
!                  = achtel * ep2r * smsq * tmtq
      vector(1,1,3) = achtel * 2._dp * smsq * tmtq
      vector(1,2,3) = achtel * ep2r * em2s * tmtq
      vector(1,3,3) = achtel * ep2r * smsq * em2t

!      vector(2,3) = achtel * r * t * epr * em2s * emt
!                  = achtel * rprq * em2s * tmtq
      vector(2,1,3) = achtel * ep2r * em2s * tmtq
      vector(2,2,3) = achtel * rprq * (-2._dp) * tmtq
      vector(2,3,3) = achtel * rprq * em2s * em2t

!      vector(3,3) = achtel * r * s * epr * ems * em2t
!                  = achtel * rprq * smsq * em2t
      vector(3,1,3) = achtel * ep2r * smsq * em2t
      vector(3,2,3) = achtel * rprq * em2s * em2t
      vector(3,3,3) = achtel * rprq * smsq * (-2._dp)




!      vector(1,4) = viertel * t * em2r * emsq * emt
!                  = viertel * em2r * emsq * tmtq
      vector(1,1,4) = viertel * (-2._dp) * emsq * tmtq
      vector(1,2,4) = viertel * em2r * (-2._dp)* s * tmtq
      vector(1,3,4) = viertel * em2r * emsq * em2t

!      vector(2,4) = viertel * r * t * emr * (-2._dp) * s * emt
!                  = viertel * rmrq * (-2._dp) * s * tmtq
      vector(2,1,4) = viertel * em2r * (-2._dp) * s * tmtq
      vector(2,2,4) = viertel * rmrq * (-2._dp) * tmtq
      vector(2,3,4) = viertel * rmrq * (-2._dp) * s * em2t

!      vector(3,4) = viertel * r * emr * emsq * em2t
!                  = viertel * rmrq * emsq * em2t
      vector(3,1,4) = viertel * em2r * emsq * em2t
      vector(3,2,4) = viertel * rmrq * (-2._dp) * s * em2t
      vector(3,3,4) = viertel * rmrq * emsq * (-2._dp)




!      vector(1,5) = -halb * t * (-2._dp) * r * emsq * emt
!                  = -halb * (-2._dp) * r * emsq * tmtq
      vector(1,1,5) = -halb * (-2._dp) * emsq * tmtq
      vector(1,2,5) = -halb * (-2._dp) * r * (-2._dp) * s * tmtq
      vector(1,3,5) = -halb * (-2._dp) * r * emsq * em2t

!      vector(2,5) = -halb * t * emrq * (-2._dp) * s * emt
!                  = -halb * emrq * (-2._dp) * s * tmtq
      vector(2,1,5) = -halb * (-2._dp) * r * (-2._dp) * s * tmtq
      vector(2,2,5) = -halb * emrq * (-2._dp) * tmtq
      vector(2,3,5) = -halb * emrq * (-2._dp) * s * em2t

!      vector(3,5) = -halb * emrq * emsq * em2t
      vector(3,1,5) = -halb * (-2._dp) * r * emsq * em2t
      vector(3,2,5) = -halb * emrq * (-2._dp) * s * em2t
      vector(3,3,5) = -halb * emrq * emsq * (-2._dp)




!      vector(1,6) = -viertel * t * ep2r * emsq * emt
!                  = -viertel * ep2r * emsq * tmtq
      vector(1,1,6) = -viertel * 2._dp * emsq * tmtq
      vector(1,2,6) = -viertel * ep2r * (-2._dp) * s * tmtq
      vector(1,3,6) = -viertel * ep2r * emsq * em2t


!      vector(2,6) = -viertel * r * t * epr * (-2._dp) * s * emt
!                  = -viertel * rprq * (-2._dp) * s * tmtq
      vector(2,1,6) = -viertel * ep2r * (-2._dp) * s * tmtq
      vector(2,2,6) = -viertel * rprq * (-2._dp) * tmtq
      vector(2,3,6) = -viertel * rprq * (-2._dp) * s * em2t

!      vector(3,6) = -viertel * r * epr * emsq * em2t
!                  = -viertel * rprq * emsq * em2t
      vector(3,1,6) = -viertel * ep2r * emsq * em2t
      vector(3,2,6) = -viertel * rprq * (-2._dp) * s * em2t
      vector(3,3,6) = -viertel * rprq * emsq * (-2._dp)




!      vector(1,7) = achtel * s * t * em2r * eps * emt
!                  = achtel * em2r * spsq * tmtq
      vector(1,1,7) = achtel * (-2._dp) * spsq * tmtq
      vector(1,2,7) = achtel * em2r * ep2s * tmtq
      vector(1,3,7) = achtel * em2r * spsq * em2t

!      vector(2,7) = achtel * r * t * emr * ep2s * emt
!                  = achtel * rmrq * ep2s * tmtq
      vector(2,1,7) = achtel * em2r * ep2s * tmtq
      vector(2,2,7) = achtel * rmrq * 2._dp * tmtq
      vector(2,3,7) = achtel * rmrq * ep2s * em2t

!      vector(3,7) = achtel * r * s * emr * eps * em2t
!                  = achtel * rmrq * spsq * em2t
      vector(3,1,7) = achtel * em2r * spsq * em2t
      vector(3,2,7) = achtel * rmrq * ep2s * em2t
      vector(3,3,7) = achtel * rmrq * spsq * (-2._dp)




!      vector(1,8) = -viertel * s * t * (-2._dp) * r * eps * emt
!                  = -viertel * (-2._dp) * r * spsq * tmtq
      vector(1,1,8) = -viertel * (-2._dp) * spsq * tmtq
      vector(1,2,8) = -viertel * (-2._dp) * r * ep2s * tmtq
      vector(1,3,8) = -viertel * (-2._dp) * r * spsq * em2t

!      vector(2,8) = -viertel * t * emrq * ep2s * emt
!                  = -viertel * emrq * ep2s * tmtq
      vector(2,1,8) = -viertel * (-2._dp) * r * ep2s * tmtq
      vector(2,2,8) = -viertel * emrq * 2._dp * tmtq
      vector(2,3,8) = -viertel * emrq * ep2s * em2t

!      vector(3,8) = -viertel * s * emrq * eps * em2t
!                  = -viertel * emrq * spsq * em2t
      vector(3,1,8) = -viertel * (-2._dp) * r * spsq * em2t
      vector(3,2,8) = -viertel * emrq * ep2s * em2t
      vector(3,3,8) = -viertel * emrq * spsq * (-2._dp)




!      vector(1,9) = -achtel * s * t * ep2r * eps * emt
!                  = -achtel * ep2r * spsq * tmtq
      vector(1,1,9) = -achtel * 2._dp * spsq * tmtq
      vector(1,2,9) = -achtel * ep2r * ep2s * tmtq
      vector(1,3,9) = -achtel * ep2r * spsq * em2t

!      vector(2,9) = -achtel * r * t * epr * ep2s * emt
!                  = -achtel * rprq * ep2s * tmtq
      vector(2,1,9) = -achtel * ep2r * ep2s * tmtq
      vector(2,2,9) = -achtel * rprq * 2._dp * tmtq
      vector(2,3,9) = -achtel * rprq * ep2s * em2t

!      vector(3,9) = -achtel * r * s * epr * eps * em2t
!                  = -achtel * rprq * spsq * em2t
      vector(3,1,9) = -achtel * ep2r * spsq * em2t
      vector(3,2,9) = -achtel * rprq * ep2s * em2t
      vector(3,3,9) = -achtel * rprq * spsq * (-2._dp)




!      vector(1,10) = viertel * s * em2r * ems * emtq
!                   = viertel * em2r * smsq * emtq
      vector(1,1,10) = viertel * (-2._dp) * smsq * emtq
      vector(1,2,10) = viertel * em2r * em2s * emtq
      vector(1,3,10) = viertel * em2r * smsq * (-2._dp) * t

!      vector(2,10) = viertel * r * emr * em2s * emtq
!                   = viertel * rmrq * em2s * emtq
      vector(2,1,10) = viertel * em2r * em2s * emtq
      vector(2,2,10) = viertel * rmrq * (-2._dp) * emtq
      vector(2,3,10) = viertel * rmrq * em2s * (-2._dp) * t

!      vector(3,10) = viertel * r * s * emr * ems * (-2._dp) * t
!                   = viertel * rmrq * smsq * (-2._dp) * t
      vector(3,1,10) = viertel * em2r * smsq * (-2._dp) * t
      vector(3,2,10) = viertel * rmrq * em2s * (-2._dp) * t
      vector(3,3,10) = viertel * rmrq * smsq * (-2._dp)




!      vector(1,11) = -halb * s * (-2._dp) * r * ems * emtq
!                   = -halb * (-2._dp) * r * smsq * emtq
      vector(1,1,11) = -halb * (-2._dp) * smsq * emtq
      vector(1,2,11) = -halb * (-2._dp) * r * em2s * emtq
      vector(1,3,11) = -halb * (-2._dp) * r * smsq * (-2._dp) * t

!      vector(2,11) = -halb * emrq * em2s * emtq
      vector(2,1,11) = -halb * (-2._dp) * r * em2s * emtq
      vector(2,2,11) = -halb * emrq * (-2._dp) * emtq
      vector(2,3,11) = -halb * emrq * em2s * (-2._dp) * t

!      vector(3,11) = -halb * s * emrq * ems * (-2._dp) * t
!                   = -halb * emrq * smsq * (-2._dp) * t
      vector(3,1,11) = -halb * (-2._dp) * r * smsq * (-2._dp) * t
      vector(3,2,11) = -halb * emrq * em2s * (-2._dp) * t
      vector(3,3,11) = -halb * emrq * smsq * (-2._dp)




!      vector(1,12) = -viertel * s * ep2r * ems * emtq
!                   = -viertel * ep2r * smsq * emtq
      vector(1,1,12) = -viertel * 2._dp * smsq * emtq
      vector(1,2,12) = -viertel * ep2r * em2s * emtq
      vector(1,3,12) = -viertel * ep2r * smsq * (-2._dp) * t

!      vector(2,12) = -viertel * r * epr * em2s * emtq
!                   = -viertel * rprq * em2s * emtq
      vector(2,1,12) = -viertel * ep2r * em2s * emtq
      vector(2,2,12) = -viertel * rprq * (-2._dp) * emtq
      vector(2,3,12) = -viertel * rprq * em2s * (-2._dp) * t

!      vector(3,12) = -viertel * r * s * epr * ems * (-2._dp) * t
!                   = -viertel * rprq * smsq * (-2._dp) * t
      vector(3,1,12) = -viertel * ep2r * smsq * (-2._dp) * t
      vector(3,2,12) = -viertel * rprq * em2s * (-2._dp) * t
      vector(3,3,12) = -viertel * rprq * smsq * (-2._dp)




!     vector(1,13) = -halb * em2r * emsq * emtq
      vector(1,1,13) = -halb * (-2._dp) * emsq * emtq
      vector(1,2,13) = -halb * em2r * (-2._dp) * s * emtq
      vector(1,3,13) = -halb * em2r * emsq * (-2._dp) * t

!      vector(2,13) = -halb * r * emr * (-2._dp) * s * emtq
!                   = -halb * rmrq * (-2._dp) * s * emtq
      vector(2,1,13) = -halb * em2r * (-2._dp) * s * emtq
      vector(2,2,13) = -halb * rmrq * (-2._dp) * emtq
      vector(2,3,13) = -halb * rmrq * (-2._dp) * s * (-2._dp) * t

!      vector(3,13) = -halb * r * emr * emsq * (-2._dp) * t
!                   = -halb * rmrq * emsq * (-2._dp) * t
      vector(3,1,13) = -halb * em2r * emsq * (-2._dp) * t
      vector(3,2,13) = -halb * rmrq * (-2._dp) * s * (-2._dp) * t
      vector(3,3,13) = -halb * rmrq * emsq * (-2._dp)




!      vector(1,14) = (-2._dp) * r * emsq * emtq
      vector(1,1,14) = (-2._dp) * emsq * emtq
      vector(1,2,14) = (-2._dp) * r * (-2._dp) * s * emtq
      vector(1,3,14) = (-2._dp) * r * emsq * (-2._dp) * t

!      vector(2,14) = emrq * (-2._dp) * s * emtq
      vector(2,1,14) = (-2._dp) * r * (-2._dp) * s * emtq
      vector(2,2,14) = emrq * (-2._dp) * emtq
      vector(2,3,14) = emrq * (-2._dp) * s * (-2._dp) * t

!      vector(3,14) = emrq * emsq * (-2._dp) * t
      vector(3,1,14) = (-2._dp) * r * emsq * (-2._dp) * t
      vector(3,2,14) = emrq * (-2._dp) * s * (-2._dp) * t
      vector(3,3,14) = emrq * emsq * (-2._dp)




!      vector(1,15) = halb * ep2r * emsq * emtq
      vector(1,1,15) = halb * 2._dp * emsq * emtq
      vector(1,2,15) = halb * ep2r * (-2._dp) * s * emtq
      vector(1,3,15) = halb * ep2r * emsq * (-2._dp) * t

!      vector(2,15) = halb * r * epr * (-2._dp) * s * emtq
!                   = halb * rprq * (-2._dp) * s * emtq
      vector(2,1,15) = halb * ep2r * (-2._dp) * s * emtq
      vector(2,2,15) = halb * rprq * (-2._dp) * emtq
      vector(2,3,15) = halb * rprq * (-2._dp) * s * (-2._dp) * t

!      vector(3,15) = halb * r * epr * emsq * (-2._dp) * t
!                   = halb * rprq * emsq * (-2._dp) * t
      vector(3,1,15) = halb * ep2r * emsq * (-2._dp) * t
      vector(3,2,15) = halb * rprq * (-2._dp) * s * (-2._dp) * t
      vector(3,3,15) = halb * rprq * emsq * (-2._dp)




!      vector(1,16) = -viertel * s * em2r * eps * emtq
!                   = -viertel * em2r * spsq * emtq
      vector(1,1,16) = -viertel * (-2._dp) * spsq * emtq
      vector(1,2,16) = -viertel * em2r * ep2s * emtq
      vector(1,3,16) = -viertel * em2r * spsq * (-2._dp) * t

!      vector(2,16) = -viertel * r * emr * ep2s * emtq
!                   = -viertel * rmrq * ep2s * emtq
      vector(2,1,16) = -viertel * em2r * ep2s * emtq
      vector(2,2,16) = -viertel * rmrq * 2._dp * emtq
      vector(2,3,16) = -viertel * rmrq * ep2s * (-2._dp ) * t

!      vector(3,16) = -viertel * r * s * emr * eps * (-2._dp) * t
!                   = -viertel * rmrq * spsq * (-2._dp) * t
      vector(3,1,16) = -viertel * em2r * spsq * (-2._dp) * t
      vector(3,2,16) = -viertel * rmrq * ep2s * (-2._dp) * t
      vector(3,3,16) = -viertel * rmrq * spsq * (-2._dp)




!      vector(1,17) = halb * s * (-2._dp) * r * eps * emtq
!                   = halb * (-2._dp) * r * spsq * emtq
      vector(1,1,17) = halb * (-2._dp) * spsq * emtq
      vector(1,2,17) = halb * (-2._dp) * r * ep2s * emtq
      vector(1,3,17) = halb * (-2._dp) * r * spsq * (-2._dp) * t

!      vector(2,17) = halb * emrq * ep2s * emtq
      vector(2,1,17) = halb * (-2._dp) * r * ep2s * emtq
      vector(2,2,17) = halb * emrq * 2._dp * emtq
      vector(2,3,17) = halb * emrq * ep2s * (-2._dp) * t

!      vector(3,17) = halb * s * emrq * eps * (-2._dp) * t
!                   = halb * emrq * spsq * (-2._dp) * t
      vector(3,1,17) = halb * (-2._dp) * r * spsq * (-2._dp) * t
      vector(3,2,17) = halb * emrq * ep2s * (-2._dp) * t
      vector(3,3,17) = halb * emrq * spsq * (-2._dp)




!      vector(1,18) = viertel * s * ep2r * eps * emtq
!                   = viertel * ep2r * spsq * emtq
      vector(1,1,18) = viertel * 2._dp * spsq * emtq
      vector(1,2,18) = viertel * ep2r * ep2s * emtq
      vector(1,3,18) = viertel * ep2r * spsq * (-2._dp) * t

!      vector(2,18) = viertel * r * epr * ep2s * emtq
!                   = viertel * rprq * ep2s * emtq
      vector(2,1,18) = viertel * ep2r * ep2s * emtq
      vector(2,2,18) = viertel * rprq * 2._dp * emtq
      vector(2,3,18) = viertel * rprq * ep2s * (-2._dp) * t

!      vector(3,18) = viertel * r * s * epr * eps * (-2._dp) * t
!                   = viertel * rprq * spsq * (-2._dp) * t
      vector(3,1,18) = viertel * ep2r * spsq * (-2._dp) * t
      vector(3,2,18) = viertel * rprq * ep2s * (-2._dp) * t
      vector(3,3,18) = viertel * rprq * spsq * (-2._dp)




!      vector(1,19) = achtel * s * t * em2r * ems * ept
!                   = achtel * em2r * smsq * tptq
      vector(1,1,19) = achtel * (-2._dp) * smsq * tptq
      vector(1,2,19) = achtel * em2r * em2s * tptq
      vector(1,3,19) = achtel * em2r * smsq * ep2t

!      vector(2,19) = achtel * r * t * emr * em2s * ept
!                   = achtel * rmrq * em2s * tptq
      vector(2,1,19) = achtel * em2r * em2s * tptq
      vector(2,2,19) = achtel * rmrq * (-2._dp) * tptq
      vector(2,3,19) = achtel * rmrq * em2s * ep2t

!      vector(3,19) = achtel * r * s * emr * ems * ep2t
!                   = achtel * rmrq * smsq * ep2t
      vector(3,1,19) = achtel * em2r * smsq * ep2t
      vector(3,2,19) = achtel * rmrq * em2s * ep2t
      vector(3,3,19) = achtel * rmrq * smsq * 2._dp




!      vector(1,20) = -viertel * s * t * (-2._dp) * r * ems * ept
!                   = -viertel * (-2._dp) * r * smsq * tptq
      vector(1,1,20) = -viertel * (-2._dp) * smsq * tptq
      vector(1,2,20) = -viertel * (-2._dp) * r * em2s * tptq
      vector(1,3,20) = -viertel * (-2._dp) * r * smsq * ep2t

!      vector(2,20) = -viertel * t * emrq * em2s * ept
!                   = -viertel * emrq * em2s * tptq
      vector(2,1,20) = -viertel * (-2._dp) * r * em2s * tptq
      vector(2,2,20) = -viertel * emrq * (-2._dp) * tptq
      vector(2,3,20) = -viertel * emrq * em2s * ep2t

!      vector(3,20) = -viertel * s * emrq * ems * ep2t
!                   = -viertel * emrq * smsq * ep2t
      vector(3,1,20) = -viertel * (-2._dp) * r * smsq * ep2t
      vector(3,2,20) = -viertel * emrq * em2s * ep2t
      vector(3,3,20) = -viertel * emrq * smsq * 2._dp




!      vector(1,21) = -achtel * s * t * ep2r * ems * ept
!                   = -achtel * ep2r * smsq * tptq
      vector(1,1,21) = -achtel * 2._dp * smsq * tptq
      vector(1,2,21) = -achtel * ep2r * em2s * tptq
      vector(1,3,21) = -achtel * ep2r * smsq * ep2t

!      vector(2,21) = -achtel * r * t * epr * em2s * ept
!                   = -achtel * rprq * em2s * tptq
      vector(2,1,21) = -achtel * ep2r * em2s * tptq
      vector(2,2,21) = -achtel * rprq * (-2._dp) * tptq
      vector(2,3,21) = -achtel * rprq * em2s * ep2t

!      vector(3,21) = -achtel * r * s * epr * ems * ep2t
!                   = -achtel * rprq * smsq * ep2t
      vector(3,1,21) = -achtel * ep2r * smsq * ep2t
      vector(3,2,21) = -achtel * rprq * em2s * ep2t
      vector(3,3,21) = -achtel * rprq * smsq * 2._dp




!      vector(1,22) = -viertel * t * em2r * emsq * ept
!                   = -viertel * em2r * emsq * tptq
      vector(1,1,22) = -viertel * (-2._dp) * emsq * tptq
      vector(1,2,22) = -viertel * em2r * (-2._dp) * s * tptq
      vector(1,3,22) = -viertel * em2r * emsq * ep2t

!      vector(2,22) = -viertel * r * t * emr * (-2._dp) * s * ept
!                   = -viertel * rmrq * (-2._dp) * s * tptq
      vector(2,1,22) = -viertel * em2r * (-2._dp) * s * tptq
      vector(2,2,22) = -viertel * rmrq * (-2._dp) * tptq
      vector(2,3,22) = -viertel * rmrq * (-2._dp) * s * ep2t

!      vector(3,22) = -viertel * r * emr * emsq * ep2t
!                   = -viertel * rmrq * emsq * ep2t
      vector(3,1,22) = -viertel * em2r * emsq * ep2t
      vector(3,2,22) = -viertel * rmrq * (-2._dp) * s * ep2t
      vector(3,3,22) = -viertel * rmrq * emsq * 2._dp




!      vector(1,23) = halb * t * (-2._dp) * r * emsq * ept
!                   = halb * (-2._dp) * r * emsq * tptq
      vector(1,1,23) = halb * (-2._dp) * emsq * tptq
      vector(1,2,23) = halb * (-2._dp) * r * (-2._dp) * s * tptq
      vector(1,3,23) = halb * (-2._dp) * r * emsq * ep2t

!      vector(2,23) = halb * t * emrq * (-2._dp) * s * ept
!                   = halb * emrq * (-2._dp) * s * tptq
      vector(2,1,23) = halb * (-2._dp) * r * (-2._dp) * s * tptq
      vector(2,2,23) = halb * emrq * (-2._dp) * tptq
      vector(2,3,23) = halb * emrq * (-2._dp) * s * ep2t

!      vector(3,23) = halb * emrq * emsq * ep2t
      vector(3,1,23) = halb * (-2._dp) * r * emsq * ep2t
      vector(3,2,23) = halb * emrq * (-2._dp) * s * ep2t
      vector(3,3,23) = halb * emrq * emsq * 2._dp




!      vector(1,24) = viertel * t * ep2r * emsq * ept
!                   = viertel * ep2r * emsq * tptq
      vector(1,1,24) = viertel * 2._dp * emsq * tptq
      vector(1,2,24) = viertel * ep2r * (-2._dp) * s * tptq
      vector(1,3,24) = viertel * ep2r * emsq * ep2t

!      vector(2,24) = viertel * r * t * epr * (-2._dp) * s * ept
!                   = viertel * rprq * (-2._dp) * s * tptq
      vector(2,1,24) = viertel * ep2r * (-2._dp) * s * tptq
      vector(2,2,24) = viertel * rprq * (-2._dp) * tptq
      vector(2,3,24) = viertel * rprq * (-2._dp) * s * ep2t

!      vector(3,24) = viertel * r * epr * emsq * ep2t
!                   = viertel * rprq * emsq * ep2t
      vector(3,1,24) = viertel * ep2r * emsq * ep2t
      vector(3,2,24) = viertel * rprq * (-2._dp) * s * ep2t
      vector(3,3,24) = viertel * rprq * emsq * 2._dp




!      vector(1,25) = -achtel * s * t * em2r * eps * ept
!                   = -achtel * em2r * spsq * tptq
      vector(1,1,25) = -achtel * (-2._dp) * spsq * tptq
      vector(1,2,25) = -achtel * em2r * ep2s * tptq
      vector(1,3,25) = -achtel * em2r * spsq * ep2t

!      vector(2,25) = -achtel * r * t * emr * ep2s * ept
!                   = -achtel * rmrq * ep2s * tptq
      vector(2,1,25) = -achtel * em2r * ep2s * tptq
      vector(2,2,25) = -achtel * rmrq * 2._dp * tptq
      vector(2,3,25) = -achtel * rmrq * ep2s * ep2t

!      vector(3,25) = -achtel * r * s * emr * eps * ep2t
!                   = -achtel * rmrq * spsq * ep2t
      vector(3,1,25) = -achtel * em2r * spsq * ep2t
      vector(3,2,25) = -achtel * rmrq * ep2s * ep2t
      vector(3,3,25) = -achtel * rmrq * spsq * 2._dp




!      vector(1,26) = viertel * s * t * (-2._dp) * r * eps * ept
!                   = viertel * (-2._dp) * r * spsq * tptq
      vector(1,1,26) = viertel * (-2._dp) * spsq * tptq
      vector(1,2,26) = viertel * (-2._dp) * r * ep2s * tptq
      vector(1,3,26) = viertel * (-2._dp) * r * spsq * ep2t

!      vector(2,26) = viertel * t * emrq * ep2s * ept
!                  = viertel * emrq * ep2s * tptq
      vector(2,1,26) = viertel * (-2._dp) * r * ep2s * tptq
      vector(2,2,26) = viertel * emrq * 2._dp * tptq
      vector(2,3,26) = viertel * emrq * ep2s * ep2t

!      vector(3,26) = viertel * s * emrq * eps * ep2t
!                   = viertel * emrq * spsq * ep2t
      vector(3,1,26) = viertel * (-2._dp) * r * spsq * ep2t
      vector(3,2,26) = viertel * emrq * ep2s * ep2t
      vector(3,3,26) = viertel * emrq * spsq * 2._dp




!      vector(1,27) = achtel * s * t * ep2r * eps * ept
!                   = achtel * ep2r * spsq * tptq
      vector(1,1,27) = achtel * 2._dp * spsq * tptq
      vector(1,2,27) = achtel * ep2r * ep2s * tptq
      vector(1,3,27) = achtel * ep2r * spsq * ep2t

!      vector(2,27) = achtel * r * t * epr * ep2s * ept
!                   = achtel * rprq * ep2s * tptq
      vector(2,1,27) = achtel * ep2r * ep2s * tptq
      vector(2,2,27) = achtel * rprq * 2._dp * tptq
      vector(2,3,27) = achtel * rprq * ep2s * ep2t

!      vector(3,27) = achtel * r * s * epr * eps * ep2t
!                   = achtel * rprq * spsq * ep2t
      vector(3,1,27) = achtel * ep2r * spsq * ep2t
      vector(3,2,27) = achtel * rprq * ep2s * ep2t
      vector(3,3,27) = achtel * rprq * spsq * 2._dp


      return
      end subroutine df2triqua
C ===== SOURCE: dftriqua.f
      subroutine dftriqua (vector, r, s, t)

      USE PRECISION

      implicit none
      real(dp), intent(out) :: vector(3,27)
      real(dp), intent(in) :: r, s, t
      real(dp) :: halb, viertel, achtel, emr, ems, emt, epr, eps, ept,
     .            emrq, emsq, emtq, em2r, em2s, em2t, ep2r, ep2s, ep2t

      halb = 0.5_dp
      viertel = 0.25_dp
      achtel = 0.125_dp

      emr = 1._dp - r
      ems = 1._dp - s
      emt = 1._dp - t

      em2r = 1._dp - 2._dp * r
      em2s = 1._dp - 2._dp * s
      em2t = 1._dp - 2._dp * t

      epr = 1._dp + r
      eps = 1._dp + s
      ept = 1._dp + t

      ep2r = 1._dp + 2._dp * r
      ep2s = 1._dp + 2._dp * s
      ep2t = 1._dp + 2._dp * t

      emrq = 1._dp - r*r
      emsq = 1._dp - s*s
      emtq = 1._dp - t*t

      vector(1,1) = -achtel * s * t * em2r * ems * emt
      vector(2,1) = -achtel * r * t * emr * em2s * emt
      vector(3,1) = -achtel * r * s * emr * ems * em2t

      vector(1,2) = viertel * s * t * (-2._dp) * r * ems * emt
      vector(2,2) = viertel * t * emrq * em2s * emt
      vector(3,2) = viertel * s * emrq * ems * em2t

      vector(1,3) = achtel * s * t * ep2r * ems * emt
      vector(2,3) = achtel * r * t * epr * em2s * emt
      vector(3,3) = achtel * r * s * epr * ems * em2t

      vector(1,4) = viertel * t * em2r * emsq * emt
      vector(2,4) = viertel * r * t * emr * (-2._dp) * s * emt
      vector(3,4) = viertel * r * emr * emsq * em2t

      vector(1,5) = -halb * t * (-2._dp) * r * emsq * emt
      vector(2,5) = -halb * t * emrq * (-2._dp) * s * emt
      vector(3,5) = -halb * emrq * emsq * em2t

      vector(1,6) = -viertel * t * ep2r * emsq * emt
      vector(2,6) = -viertel * r * t * epr * (-2._dp) * s * emt
      vector(3,6) = -viertel * r * epr * emsq * em2t

      vector(1,7) = achtel * s * t * em2r * eps * emt
      vector(2,7) = achtel * r * t * emr * ep2s * emt
      vector(3,7) = achtel * r * s * emr * eps * em2t

      vector(1,8) = -viertel * s * t * (-2._dp) * r * eps * emt
      vector(2,8) = -viertel * t * emrq * ep2s * emt
      vector(3,8) = -viertel * s * emrq * eps * em2t

      vector(1,9) = -achtel * s * t * ep2r * eps * emt
      vector(2,9) = -achtel * r * t * epr * ep2s * emt
      vector(3,9) = -achtel * r * s * epr * eps * em2t


      vector(1,10) = viertel * s * em2r * ems * emtq
      vector(2,10) = viertel * r * emr * em2s * emtq
      vector(3,10) = viertel * r * s * emr * ems * (-2._dp) * t

      vector(1,11) = -halb * s * (-2._dp) * r * ems * emtq
      vector(2,11) = -halb * emrq * em2s * emtq
      vector(3,11) = -halb * s * emrq * ems * (-2._dp) * t

      vector(1,12) = -viertel * s * ep2r * ems * emtq
      vector(2,12) = -viertel * r * epr * em2s * emtq
      vector(3,12) = -viertel * r * s * epr * ems * (-2._dp) * t

      vector(1,13) = -halb * em2r * emsq * emtq
      vector(2,13) = -halb * r * emr * (-2._dp) * s * emtq
      vector(3,13) = -halb * r * emr * emsq * (-2._dp) * t

      vector(1,14) = (-2._dp) * r * emsq * emtq
      vector(2,14) = emrq * (-2._dp) * s * emtq
      vector(3,14) = emrq * emsq * (-2._dp) * t

      vector(1,15) = halb * ep2r * emsq * emtq
      vector(2,15) = halb * r * epr * (-2._dp) * s * emtq
      vector(3,15) = halb * r * epr * emsq * (-2._dp) * t

      vector(1,16) = -viertel * s * em2r * eps * emtq
      vector(2,16) = -viertel * r * emr * ep2s * emtq
      vector(3,16) = -viertel * r * s * emr * eps * (-2._dp) * t

      vector(1,17) = halb * s * (-2._dp) * r * eps * emtq
      vector(2,17) = halb * emrq * ep2s * emtq
      vector(3,17) = halb * s * emrq * eps * (-2._dp) * t

      vector(1,18) = viertel * s * ep2r * eps * emtq
      vector(2,18) = viertel * r * epr * ep2s * emtq
      vector(3,18) = viertel * r * s * epr * eps * (-2._dp) * t


      vector(1,19) = achtel * s * t * em2r * ems * ept
      vector(2,19) = achtel * r * t * emr * em2s * ept
      vector(3,19) = achtel * r * s * emr * ems * ep2t

      vector(1,20) = -viertel * s * t * (-2._dp) * r * ems * ept
      vector(2,20) = -viertel * t * emrq * em2s * ept
      vector(3,20) = -viertel * s * emrq * ems * ep2t

      vector(1,21) = -achtel * s * t * ep2r * ems * ept
      vector(2,21) = -achtel * r * t * epr * em2s * ept
      vector(3,21) = -achtel * r * s * epr * ems * ep2t

      vector(1,22) = -viertel * t * em2r * emsq * ept
      vector(2,22) = -viertel * r * t * emr * (-2._dp) * s * ept
      vector(3,22) = -viertel * r * emr * emsq * ep2t

      vector(1,23) = halb * t * (-2._dp) * r * emsq * ept
      vector(2,23) = halb * t * emrq * (-2._dp) * s * ept
      vector(3,23) = halb * emrq * emsq * ep2t

      vector(1,24) = viertel * t * ep2r * emsq * ept
      vector(2,24) = viertel * r * t * epr * (-2._dp) * s * ept
      vector(3,24) = viertel * r * epr * emsq * ep2t

      vector(1,25) = -achtel * s * t * em2r * eps * ept
      vector(2,25) = -achtel * r * t * emr * ep2s * ept
      vector(3,25) = -achtel * r * s * emr * eps * ep2t

      vector(1,26) = viertel * s * t * (-2._dp) * r * eps * ept
      vector(2,26) = viertel * t * emrq * ep2s * ept
      vector(3,26) = viertel * s * emrq * eps * ep2t

      vector(1,27) = achtel * s * t * ep2r * eps * ept
      vector(2,27) = achtel * r * t * epr * ep2s * ept
      vector(3,27) = achtel * r * s * epr * eps * ep2t

      return
      end subroutine dftriqua
C ===== SOURCE: eirsrt.f
C

      SUBROUTINE EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)
      USE PRECISION
      USE PARMMOD
      USE BRASPOI
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE CTRCEI
      USE CCOUPL
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: LSTOP, LTIME
      REAL(DP), INTENT(IN) :: DELTAT, FLUXES(*),
     .                      B2BRM, B2RD, B2Q, B2VP
C
C
C
      RETURN

      END
C ===== SOURCE: flaplace.f
      function flaplace (r,s,t,xvec,yvec,zvec,tevec)

      use precision
      use ccona
      
      implicit none

      real(dp), intent(in) :: r, s, t
      real(dp), intent(in), dimension(27) :: xvec, yvec, zvec, tevec

      real(dp) :: u1(3), u2(3), u3(3), eu1(3),eu2(3),eu3(3),guik(3,3),
     .            gradphi(3), eo1(3), eo2(3), eo3(3), gradphi2(3),
     .            goik(3,3), dfguik(3,3,3)
      real(dp) :: dxdr, dxds, dxdt, dydr, dyds, dydt, dzdr, dzds, dzdt,
     .            e123, gikd2te, sumj, dtedxj, flaplace
      integer :: ii, j, k, l

      REAL(DP), save :: DFVEC(3,27), DF2VEC(3,3,27)
      real(dp), save :: rold=1.E20_DP, sold=1.E20_DP, told=1.E20_DP

      
      if (abs(r-rold) + abs(s-sold) + abs(t-told) > eps30) then

! berechne 1. Ableitungen der Shapefunktionen
         CALL DFTRIQUA(DFVEC, r, s, t)

! berechne 2. Ableitungen der Shapefunktionen
         CALL DF2TRIQUA(DF2VEC, r, s, t)

         rold = r
         sold = s
         told = t

      end if

      dxdr=sum(xvec*dfvec(1,:))
      dxds=sum(xvec*dfvec(2,:))
      dxdt=sum(xvec*dfvec(3,:))
      dydr=sum(yvec*dfvec(1,:))
      dyds=sum(yvec*dfvec(2,:))
      dydt=sum(yvec*dfvec(3,:))
      dzdr=sum(zvec*dfvec(1,:))
      dzds=sum(zvec*dfvec(2,:))
      dzdt=sum(zvec*dfvec(3,:))
      
      eu1 = (/ dxdr, dydr, dzdr /)
      eu2 = (/ dxds, dyds, dzds /)
      eu3 = (/ dxdt, dydt, dzdt /)
      
      e123 = eu1(1)*eu2(2)*eu3(3) + eu2(1)*eu3(2)*eu1(3) +
     .       eu3(1)*eu1(2)*eu2(3) - eu1(3)*eu2(2)*eu3(1) -
     .       eu2(3)*eu3(2)*eu1(1) - eu3(3)*eu1(2)*eu2(1)
      
      call kreuzprod(eu2,eu3,e123,eo1)
      call kreuzprod(eu3,eu1,e123,eo2)
      call kreuzprod(eu1,eu2,e123,eo3)
      
      guik(1,1)=dxdr*dxdr + dydr*dydr + dzdr*dzdr
      guik(1,2)=dxdr*dxds + dydr*dyds + dzdr*dzds
      guik(1,3)=dxdr*dxdt + dydr*dydt + dzdr*dzdt
      
      guik(2,1)=dxds*dxdr + dyds*dydr + dzds*dzdr
      guik(2,2)=dxds*dxds + dyds*dyds + dzds*dzds
      guik(2,3)=dxds*dxdt + dyds*dydt + dzds*dzdt
      
      guik(3,1)=dxdt*dxdr + dydt*dydr + dzdt*dzdr
      guik(3,2)=dxdt*dxds + dydt*dyds + dzdt*dzds
      guik(3,3)=dxdt*dxdt + dydt*dydt + dzdt*dzdt
      
      goik(1,1) = sum(eo1*eo1)
      goik(1,2) = sum(eo1*eo2)
      goik(1,3) = sum(eo1*eo3)
      
      goik(2,1) = sum(eo2*eo1)
      goik(2,2) = sum(eo2*eo2)
      goik(2,3) = sum(eo2*eo3)
      
      goik(3,1) = sum(eo3*eo1)
      goik(3,2) = sum(eo3*eo2)
      goik(3,3) = sum(eo3*eo3)
      
      u1 = 1._dp/sqrt(guik(1,1))*eu1
      u2 = 1._dp/sqrt(guik(2,2))*eu2
      u3 = 1._dp/sqrt(guik(3,3))*eu3
      
!     berechne Ableitungen der guik
      do ii=1,3
         do k=1,3
            do l=1,3
!     dguik(i,k) / dx_l
               dfguik(ii,k,l)=sum(xvec*df2vec(ii,l,:))*
     .              sum(xvec*dfvec(k,:))
     .              + sum(xvec*dfvec(ii,:))*
     .              sum(xvec*df2vec(k,l,:))
     .              + sum(yvec*df2vec(ii,l,:))*
     .              sum(yvec*dfvec(k,:)) 
     .              + sum(yvec*dfvec(ii,:))*
     .              sum(yvec*df2vec(k,l,:))
     .              + sum(zvec*df2vec(ii,l,:))*
     .              sum(zvec*dfvec(k,:)) 
     .              + sum(zvec*dfvec(ii,:))*
     .              sum(zvec*df2vec(k,l,:))
            end do
         end do
      end do
      
      gikd2te=0._dp
      do ii=1,3
         do k=1,3
            gikd2te = gikd2te +
     .           goik(ii,k) * sum(tevec*df2vec(ii,k,:))
         end do
      end do
      
      sumj = 0._dp
      do j = 1, 3
         dtedxj = sum(tevec*dfvec(j,:))
         do l = 1, 3
            do ii = 1,3
               do k = 1, 3
                  sumj = sumj + goik(ii,k) * goik(j,l) * dtedxj *
     .                 (dfguik(ii,l,k) + dfguik(k,l,ii) -
     .                 dfguik(ii,k,l))
               end do
            end do
         end do
      end do
      
      flaplace = gikd2te - 0.5_dp * sumj
      
      return
      end function flaplace
C ===== SOURCE: gradf.f
      subroutine gradf (r,s,t,xvec,yvec,zvec,tevec,dtedx,dtedy,dtedz)

      use precision
      use ccona

      implicit none

      real(dp), intent(in) :: r, s, t
      real(dp), intent(in), dimension(27) :: xvec, yvec, zvec, tevec
      real(dp), intent(out) :: dtedx,dtedy,dtedz

      REAL(DP), save :: DFVEC(3,27)
      REAL(DP) ::AJ(3,3), AJM1(3,3), AJM1T(3,3)
      real(dp), save :: rold=1.E20_DP, sold=1.E20_DP, told=1.E20_DP

      real(dp) :: deti

      
      if (abs(r-rold) + abs(s-sold) + abs(t-told) > eps30) then

! berechne 1. Ableitungen der Shapefunktionen
         CALL DFTRIQUA(DFVEC, r, s, t)

         rold = r
         sold = s
         told = t

      end if

      AJ(1,1) = SUM(XVEC*DFVEC(1,:))
      AJ(1,2) = SUM(YVEC*DFVEC(1,:))
      AJ(1,3) = SUM(ZVEC*DFVEC(1,:))
      AJ(2,1) = SUM(XVEC*DFVEC(2,:))
      AJ(2,2) = SUM(YVEC*DFVEC(2,:))
      AJ(2,3) = SUM(ZVEC*DFVEC(2,:))
      AJ(3,1) = SUM(XVEC*DFVEC(3,:))
      AJ(3,2) = SUM(YVEC*DFVEC(3,:))
      AJ(3,3) = SUM(ZVEC*DFVEC(3,:))

      AJM1T(1,1) = AJ(2,2)*AJ(3,3) - AJ(2,3)*AJ(3,2)
      AJM1T(1,2) = AJ(2,3)*AJ(3,1) - AJ(2,1)*AJ(3,3)
      AJM1T(1,3) = AJ(2,1)*AJ(3,2) - AJ(3,1)*AJ(2,2)
      AJM1T(2,1) = AJ(3,2)*AJ(1,3) - AJ(1,2)*AJ(3,3)
      AJM1T(2,2) = AJ(3,3)*AJ(1,1) - AJ(3,1)*AJ(1,3)
      AJM1T(2,3) = AJ(3,1)*AJ(1,2) - AJ(3,2)*AJ(1,1)
      AJM1T(3,1) = AJ(1,2)*AJ(2,3) - AJ(1,3)*AJ(2,2)
      AJM1T(3,2) = AJ(1,3)*AJ(2,1) - AJ(2,3)*AJ(1,1)
      AJM1T(3,3) = AJ(1,1)*AJ(2,2) - AJ(1,2)*AJ(2,1)
      
      DETI = 1._DP / (AJ(1,1)*AJ(2,2)*AJ(3,3)
     .              + AJ(1,2)*AJ(2,3)*AJ(3,1)
     .              + AJ(1,3)*AJ(2,1)*AJ(3,2)
     .              - AJ(3,1)*AJ(2,2)*AJ(1,3)
     .              - AJ(3,2)*AJ(2,3)*AJ(1,1)
     .              - AJ(3,3)*AJ(2,1)*AJ(1,2))

      AJM1 = TRANSPOSE(AJM1T) * DETI
        
      DTEDX = SUM(TEVEC*DFVEC(1,:)) * AJM1(1,1) +
     .        SUM(TEVEC*DFVEC(2,:)) * AJM1(1,2) +
     .        SUM(TEVEC*DFVEC(3,:)) * AJM1(1,3)
      DTEDY = SUM(TEVEC*DFVEC(1,:)) * AJM1(2,1) +
     .        SUM(TEVEC*DFVEC(2,:)) * AJM1(2,2) +
     .        SUM(TEVEC*DFVEC(3,:)) * AJM1(2,3)
      DTEDZ = SUM(TEVEC*DFVEC(1,:)) * AJM1(3,1) +
     .        SUM(TEVEC*DFVEC(2,:)) * AJM1(3,2) +
     .        SUM(TEVEC*DFVEC(3,:)) * AJM1(3,3)

      return
      end subroutine gradf
C ===== SOURCE: if0prm.f


      SUBROUTINE IF0PRM(IUNIN)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CINIT
      USE CTETRA

      IMPLICIT NONE

      INTEGER :: NBRICK, NQUAD, NLINE
      COMMON /CINTFA/ NBRICK, NQUAD, NLINE

      INTEGER, INTENT(IN) :: IUNIN
      INTEGER :: NFLA, NCUTB, NCUTL, NDXA, NDYA, IPL, NTARGI, IT, IPRT,
     .           NAINB, IAIN, NAOTB, IAOT, NTARGSTP,
     .           NCOPI, NOD, NTT, I1, I2, I3
      INTEGER, ALLOCATABLE :: NTGPRT(:)
      CHARACTER(72) :: ZEILE
C
      READ (IUNIN,*)
      DO IPL=1,NPLS
        READ (IUNIN,*)
      END DO
C  NUMBER OF TARGET SOURCES ON B2 SURFACES: NTARGI
      READ (IUNIN,'(I6)') NTARGI
!pb      NSTEP = MAX(NSTEP,NTARGI)
C  NUMBER OF PARTS PER TARGET SOURCE
      NTARGSTP=0
      IF (NTARGI.GT.0) THEN
        DO IT=1,NTARGI
          READ (IUNIN,*) I1, I2, I3
          IF (I3 /= 0) NTARGSTP = NTARGSTP + 1
        END DO
      END IF
      NSTEP = MAX(NSTEP,NTARGSTP)
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM B2 INTO EIRENE
C  HERE: B2 VOLUME TALLIES
      READ (IUNIN,'(I6)') NAINB
      NAIN = MAX(NAIN,NAINB)
      DO IAIN=1,NAINB
        READ (IUNIN,*)
        READ (IUNIN,*)
        READ (IUNIN,*)
      END DO
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO B2
C  HERE: EIRENE SURFACE TALLIES
      READ (IUNIN,'(I6)') NAOTB
      DO IAOT=1,NAOTB
        READ (IUNIN,*)
      END DO
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      NCOPI=4
      NCPVI=max(NCOPI*NPLS,13)
      NCOP = NCOPI
      NCPV = MAX(NCPVI,13)
C
      OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 30

      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'NO. OF NODES') == 0)
         READ (30,'(A72)') ZEILE
      END DO
      READ (30,*) NCOORD

      NBRICK = 0
      NQUAD = 0
      NLINE = 0
      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'ELEMENT GROUPS') == 0)
         READ (30,'(A72)') ZEILE
      END DO

      DO
         READ (30,'(A72)',END=9) ZEILE
         I1 = INDEX(ZEILE,'ELEMENTS:')
         DO WHILE (I1 == 0)
            READ (30,'(A72)',END=9) ZEILE
            I1 = INDEX(ZEILE,'ELEMENTS:')
         END DO
         I2 = INDEX(ZEILE,'NODES:')
         I3 = INDEX(ZEILE,'GEOMETRY:')
         READ (ZEILE(I2+7:I3-1),*) NOD
         READ (ZEILE(I1+10:I2),*) NTT
         IF (NOD == 27) THEN
           NBRICK = NBRICK + NTT
         ELSE IF (NOD == 9) THEN
           NQUAD = NQUAD + NTT
         ELSE IF (NOD == 3) THEN
           NLINE = NLINE + NTT
         ELSE
           EXIT
         END IF
      END DO

 9    CONTINUE
      NTETRA = NBRICK*48
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      NDX = 1
      NDY = 1
      NFL = 1
      NDXP = NDX+1
      NDYP = NDY+1

      RETURN
      END
C ===== SOURCE: indmap.f
C
C
      SUBROUTINE INDMAP(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG)
C
C     INDEX MAPPING FOR BRAAMS DATA FIELDS. DATA IN DUMMY ZONES
C     (CUTS OR BOUNDARY ZONES) MAY BE NEEDED AND THUS ARE KEPT
C     AND DUBLICATED IN CASE NCUTL GT NCUTB
C
C     NCUTB= NUMBER OF CELLS IN IX DIRECTION PER CUT IN BRAAMS
C     NCUTL= NUMBER OF CELLS IN IX DIRECTION PER CUT IN LINDA (AND
C            THUS ALSO IN EIRENE) GEOMETRY
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX,NDY,NFL,NDXA,NDYA,NFLA,NCUTB,NCUTL,
     .                       NPPLG
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL)
      REAL(DP), INTENT(INOUT) :: DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IY, IF, IENDD, INB, IINID, IINIV, IENDV, IPART
C
C  LOOP FOR THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=FIELD(IX,IY,IF)
C
C
C      NDX DIRECTION: IX=0: NOT MODIFIED
C                     IX=I(CUT): USE CUT VALUE
C                     IX=I(LAST X ZONE): MOVE TO NDXA+1
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(IX,IY)=FIELD(INB,IY,IF)
213         CONTINUE
            DUMMY(IINID,IY) = FIELD(INB+1,IY,IF)
            IF (IENDD.NE.IINID) DUMMY(IENDD,IY) = FIELD(INB+NCUTB,IY,IF)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMAP: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMAP: INCONSISTENCY IN NUMBER OF '
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED.  '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .   IENDD-IINID+1
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: indmpi.f
C
C
      SUBROUTINE INDMPI(FIELD,DUMMY,NDX,NDY,NFL,NDXA,NDYA,NFLA,
     .                  NCUTB,NCUTL,NPOINT,NPPLG,NSTR,ISTR)
C
C     INDEX MAPPING: INVERS TO SUBR. INDMAP
C

      USE PRECISION
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NPOINT(2,*)
      INTEGER, INTENT(IN) :: NDX,NDY,NFL,NDXA,NDYA,NFLA,NCUTB,NCUTL,
     .                       NPPLG,NSTR,ISTR
      REAL(DP), INTENT(INOUT) :: FIELD(0:NDX+1,0:NDY+1,NFL,NSTR)
      REAL(DP), INTENT(INOUT) :: DUMMY(0:NDX+1,0:NDY+1)
      INTEGER :: IX, IY, IF, IENDD, INB, IINID, IINIV, IENDV, IPART
C
C  LOOP OVER THE SPECIES
C
      DO 500 IF=1,NFLA
C
C  INITIALIZE DUMMY
C
        DO 10 IY=0,NDY+1
          DO 10 IX=0,NDX+1
10          DUMMY(IX,IY)=0.
C
C
C      NDX DIRECTION
C
C  NPOINT(1,1)=1
C  NPOINT(2,NPPLG)=NDXA+1
C
        IF (NCUTB.LT.0) GOTO 990
        DO 211 IPART = 1,NPPLG
C  "VALID REGION"
          IINIV= NPOINT(1,IPART)
          IENDV= NPOINT(2,IPART)-1
C  "CUT REGION" AND LAST X ZONE IX = NDXA+1
          IF (IPART.LT.NPPLG) THEN
            IINID= NPOINT(2,IPART)
            IENDD= NPOINT(1,IPART+1)-1
            IF (IENDD-IINID+1.NE.NCUTL) GOTO 991
          ELSE
            IINID= NDXA+1
            IENDD= NDXA+1
          ENDIF
          DO 212 IY=0,NDYA+1
            DO 213 IX = IINIV,IENDV
              INB=IX-(IPART-1)*(NCUTL-NCUTB)
              DUMMY(INB,IY)=FIELD(IX,IY,IF,ISTR)
213         CONTINUE
            DUMMY(INB+1,IY)=FIELD(IINID,IY,IF,ISTR)
            IF (IENDD.NE.IINID)
     .          DUMMY(INB+NCUTB,IY)=FIELD(IENDD,IY,IF,ISTR)
212       CONTINUE
211     CONTINUE
        DO 220 IY=0,NDYA+1
          DO 220 IX=0,NDXA+1
            FIELD(IX,IY,IF,ISTR)=DUMMY(IX,IY)
220     CONTINUE
C
500   CONTINUE
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. INDMPI: THIS SUBR. IS VALID ONLY'
      WRITE (iunout,*) 'NCUTB>=0 BUT NCUTB = ',NCUTB
      CALL EXIT_OWN(1)
991   WRITE (iunout,*) 
     .  'ERROR IN SUBR. INDMPI: INCONSISTENCY IN NUMBER OF'
      WRITE (iunout,*) 'ZONES PER CUT FROM LINDA GEOMETRY DETECTED. '
      WRITE (iunout,*) 'NCUTL = ',NCUTL, ' IENDD-IINID+1 = ',
     .   IENDD-IINID+1
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: infcop.f
cc  coupling to fidap: ein brick = 48 tetraeder, so dass
cc  mindestens jeder fidap knoten eine ecke eines tetraeder ist.
cc  durch die zerlegung in 48 (nicht nur 5 oder 6) tetraeder koennen
cc  die bricks (quader) beliebig orientiert sein, es geht immer so auf,
cc  dass eine brick-seitenflaeche, zerlegt in tetraeder-dreieckseiten,
cc  auf die zerlegung der nachbarflaeche passt.
cc
C   THIS CODE SEGMENT CONTAINES VARIOUS SUBROUTINES NEEDED FOR
C   INTERFACING THE EIRENE CODE TO PLASMA FLUID CODES.
C   IT READS GEOMETRICAL DATA (MESHES) FROM FILE FT30
C   AND PRODUCES THE EIRENE INPUT DATA (BLOCK 2).
C   IT READS PLASMA BACKGROUND DATA FROM FILE (FT31) OR COMMON BLOCKS,
C   IT THEN PRODUCES INPUT DATA FOR EIRENE
C   INPUT BLOCK 5 (PLASMA DATA) AND BLOCK 7 (SURFACE RECYCLING SOURCES)
C
C
*DK COUPLE
C
      SUBROUTINE INFCOP
C
C
C     THIS SUBROUTINE DEFINES THE PLASMA MODEL IN CASE OF A COUPLED
C     NEUTRAL-PLASMA CALCULATION
C
C     THE ENTRY "IF0COP" RECEIVES GEOMETRICAL INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF1COP" RECEIVES PLASMA INPUT DATA FROM AN
C     EXTERNAL FILE (E.G. OTHER PLASMA CODES)
C     AND PREPARES THEM FOR AN EIRENE RUN
C
C     THE ENTRY "IF2COP" PREPARES THE SOURCE SAMPLING DISTRIBUTION
C     FROM THE EXTERNAL DATA, AND MAY OVERWRITE OTHER INPUT
C     DATA FROM BLOCKS 1 TO 13 AS WELL
C
C     THE ENTRIES "IF3COP, IF4COP" RETURN  RESULTS TO AN EXTERNAL FILE
C
C
C
      USE PRECISION
      USE PARMMOD
      USE BRASPOI
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CPLOT
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CCOUPL
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE CTETRA
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE CSTEP
      USE CTEXT
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI
      USE CTRIG
      USE module_avltree

      IMPLICIT NONE

      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL
C
      REAL(DP) :: SCALN(0:NFL)
      REAL(DP) :: DI(NPLS), VP(NPLS), CSP(NPLS)
      REAL(DP) :: EPEL(NRTAL), HELPAR(NRTAL), DIVI(NRTAL)
      REAL(DP) :: EMSAMP(NCOOR), EM1SMP(NCOOR), EM2SMP(NCOOR)
      REAL(DP) :: ABSAMP(NCOOR), AB1SMP(NCOOR), AB2SMP(NCOOR)
      REAL(DP) :: RIX(NCOOR), RIY(NCOOR), RIZ(NCOOR)
      REAL(DP) :: XKAP(NCOOR), YKAP(NCOOR), ZKAP(NCOOR)
      REAL(DP) :: EMIS(NCOOR), ABSORB(NCOOR), VOLSUM(NCOOR)
C
      REAL(DP) :: PUX(NRAD),PUY(NRAD),PVX(NRAD),PVY(NRAD),
     R          TORL(NSTRA,NGITT),EFLX(NSTRA),
     R          DUMMY(0:NDXP,0:NDYP),
     R          ESHT(NSTEP,NGITT),ORI(NSTEP,NGITT),
     R          SFNIT(0:NSTEP),SFEIT(0:NSTEP),
     R          SFEET(0:NSTEP),SHEAE(0:NSTEP),SHEAI(0:NSTEP)
      INTEGER, ALLOCATABLE :: NRWL(:)
C
      LOGICAL :: LSHORT=.FALSE., LSTOP=.TRUE., LSTP, LSHIFT

      INTEGER :: NBRICK, NQUAD, NLINE
      COMMON /CINTFA/ NBRICK, NQUAD, NLINE

      INTEGER, ALLOCATABLE, SAVE :: NODBRICK(:,:), NODQUAD(:,:)

      INTEGER, ALLOCATABLE :: IHELP(:), NOSTS(:), NQDSTS(:), NOSPEC(:),
     .                        NOSTEP(:), NQDSTP(:)
      REAL(DP) :: TIS(NPLS),DIS(NPLS),VXS(NPLS),VYS(NPLS),VZS(NPLS)
      REAL(SP) :: PARM(16)
      INTEGER :: ITSIDE(3,4)

      REAL(DP) :: EMAXW, ESUM, ETOT, EADD, PERWI, PARWI, FLX, FLXI, 
     .          VPX, VPY, GAMMA, SHEATH, VT, PN1, VTEST, VR, DRR,
     .          PERW, PARW, CS, CUR, X1, Y1, Z1, X2, Y2, Z2, PHI1, PHI2,
     .          STEP, THMAX, BXS, BYS, BZS, PM1, TE, VPZ, ZMFPTH,
     .          EEMAX, RP1, OR, EESHT, V, T, TES, VL, XCO, YCO, ZCO,
     .          DIST, DIST1, DIST2, DIST3, DIST4, BVAC, RECADD, EEADD,
     .          SUMN, SUMEE, RECTOT, FTABRC1, FEELRC1, XGAU, YGAU, ZGAU,
     .          DIXDX, DIYDY, DIZDZ, DETI, FAC, SIGMX
      REAL(DP) :: TQVEC(27), XVEC(27), YVEC(27), ZVEC(27), GK(3), 
     .          DFVEC(3,27), AJ(3,3), AJM1(3,3), RIXVEC(27), RIYVEC(27), 
     .          RIZVEC(27), E(3,3), AJM1T(3,3), TEVEC(27), 
     .          DF2VEC(3,3,27), chrsym(3,3,3)
      real(dp) :: dtdxyz(3,27)
      real(dp) :: u1(3), u2(3), u3(3), eu1(3),eu2(3),eu3(3),guik(3,3),
     .            gradphi(3), eo1(3), eo2(3), eo3(3), gradphi2(3), 
     .            goik(3,3), dfguik(3,3,3)
      real(dp) :: dxdr, dxds, dxdt, dydr, dyds, dydt, dzdr, dzds, dzdt,
     .            e123, dphidx, dphidy, dphidz, weifakt, sumik,
     .            xlaplace, test, testte, dtedxj,flaplace
      REAL(DP) :: TEG, DTEDXG, DTEDYG, DTEDZG, gikd2te, sumj, suml
      INTEGER :: IFIRST, NREC11, NDXY, NR1STQ, ISTRAI, IY, IRC, NFALSE,
     .           JC, J, NDXYM, NDYAM, NDXAM, N1, N2, N3IN, N3EN, ITFAL,
     .           IAOT, IAIN, IR, IP, IR1, IP1, K, MPER, N3, NT, IERROR,
     .           IPL, IMODE, IFRST, ITARG, INEW, IOLD, ISP, ISTEP,
     .           N1EN, N1IN, N2EN, N2IN, IPRT, NDZA, I, NTGPRI, IT, 
     .           IEPLS, I4, IIPLS, IG, NPES, ITRI, NTACT, IN, NTOLD, IS,
     .           ISTS, IC, IDU, I1, I2, I3, NOD, NTT, IEL, NBR, M,
     .           NQU, IQ, IMATCH, LANF, LEND, ICO, ITET, ICOLUMN, IANF,
     .           NREAD, ISRFSI, ISR, IIRC, IRRC, INC, IPO, L, II, JJ,
     .           IPLSTI, IPLSV, IPLV
      INTEGER :: INDCO(27),INDF1(9),INDF2(9),INDF3(9),INDF4(9),INDF5(9),
     .           INDF6(9),NSORQUAD(9)
      INTEGER, INTENT(IN) :: ISTRAA, ISTRAE, NEW_ITER
      INTEGER, EXTERNAL :: IDEZ
      REAL(DP), ALLOCATABLE :: TEF(:),DENF(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: DRSTDXYZ(:,:,:),
     .                         DTEDX(:),DTEDY(:),DTEDZ(:),XLAPLA(:)
C
      CHARACTER(72) :: ZEILE, NAMENT, SPECNAME
      CHARACTER(10) :: CHR, FORM
      CHARACTER(6) :: CITARG
      CHARACTER(72), ALLOCATABLE :: ENTITY(:)
      CHARACTER(1000) :: LINE
      LOGICAL, ALLOCATABLE :: LUSED(:)
C
C
      DATA ITSIDE /1,2,3,
     .             1,4,2,
     .             2,4,3,
     .             3,4,1/
C
C
      ENTRY IF0COP
C
      LSHORT=.FALSE.
C
      GOTO 99990
C
C  TO INITALIZE THE SHORT CYCLING, THE GEOMETRY HAS TO BE
C  DEFINED ONCE (ENTRY: INTER0)
C
      ENTRY INTER0
      LSHORT=.TRUE.
99990 CONTINUE
C
      IERROR=0
C
      IMODE=IABS(NMODE)
C
      IF (.NOT.LSHORT.AND.ITIMV.LE.1) THEN
        WRITE (iunout,*) '        SUBROUTINE INFCOP IS CALLED  '
C  READ INPUT DATA OF BLOCK 14
C  SAVE INPUT DATA OF BLOCK 14 FOR SHORT CYCLE ON COMMON CCOUPL
        CALL LEER(1)
        CALL ALLOC_CCOUPL(1)
        READ (IUNIN,'(5L1)') LSYMET,LBALAN
        IF (TRCINT)
     .  WRITE (iunout,*) ' LSYMET,LBALAN = ',LSYMET,LBALAN
        DO 20 IPL=1,NPLSI
          READ (IUNIN,'(2I6,2E12.4)') I,IFLB(IPL),FCTE(IPL),BMASS(IPL)
          IF (TRCINT)
     .    WRITE (iunout,*) IPL,IFLB(IPL),FCTE(IPL),BMASS(IPL)
20      CONTINUE
C  NUMBER OF DIFFERENT ENTITIES: NTARGI
        READ (IUNIN,'(I6)') NTARGI
        WRITE (iunout,*) '        NTARGI= ',NTARGI
        CALL LEER(1)
        ALLOCATE (NOSTS(NTARGI))
        ALLOCATE (NOSTEP(NTARGI))
        ALLOCATE (ENTITY(NTARGI))
        DO 30 IT=1,NTARGI
          nament=repeat(' ',72)
          READ (IUNIN,'(3I6,1X,A72)') I,NOSTS(IT),NOSTEP(IT),NAMENT
          LANF=verify(nament,' ')
          LEND=verify(nament,' ',.true.)
          ENTITY(IT)=repeat(' ',72)
          ENTITY(IT)(1:lend-lanf+1) = nament(lanf:lend)
          WRITE (iunout,*) ' ENTITY ',ENTITY(IT)
          CALL UPPERCASE (ENTITY(IT))
          WRITE (iunout,*) ' ENTITY ',ENTITY(IT)
          IF ((NOSTEP(IT) < 0).OR.(NOSTEP(IT) > NTARGI)) THEN
            WRITE (iunout,*) ' NOSTEP MUST BE >= 0 AND <= NTARGI '
            WRITE (iunout,*) ' NOSTEP, NTARGI ', NOSTEP(IT), NTARGI 
            CALL EXIT(1)
          END IF
          IF (TRCINT)
     .      WRITE (iunout,'(3I6,1X,A72)') 
     .            IT,NOSTS(IT),NOSTEP(IT),ENTITY(IT)
          IF (TRCINT) CALL LEER(1)
30      CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM FIDAP INTO EIRENE
C  HERE: FIDAP VOLUME TALLIES
        READ (IUNIN,'(I6)') NAINB
        NAIN = MAX(NAIN,NAINB)
        CALL ALLOC_CCOUPL(2)
        WRITE (iunout,*) '        NAINI = ',NAINB
        IF (NAINB.GT.NAIN) THEN
          CALL MASPRM ('NAIN',4,NAIN,'NAINB',5,NAINB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EXIT(1)
        ENDIF
        IF (TRCINT.AND.NAINB.GT.0)
     .      WRITE (iunout,*) 'I,NAINS(IAIN),NAINT(IAIN)'
        DO 40 IAIN=1,NAINB
          READ (IUNIN,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
          READ (IUNIN,'(A72)') TXTPLS(IAIN,12)
          READ (IUNIN,'(2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAINS(IAIN),NAINT(IAIN)
            WRITE (iunout,'(1X,A72)') TXTPLS(IAIN,12)
            WRITE (iunout,'(1X,2A24)') TXTPSP(IAIN,12),TXTPUN(IAIN,12)
          ENDIF
40      CONTINUE
C  READ ADDITIONAL DATA TO BE TRANSFERRED FROM EIRENE INTO FIDAP
C  HERE: EIRENE SURFACE TALLIES
        READ (IUNIN,'(I6)') NAOTB
        WRITE (iunout,*) '        NAOTI = ',NAOTB
        IF (NAOTB.GT.NLIMPS) THEN
          CALL MASPRM ('NLIMPS',6,NLIMPS,'NAOTB',5,NAOTB,IERROR)
          WRITE (iunout,*) 'EXIT CALLED FROM SUBR. INFCOP '
          CALL EXIT(1)
        ENDIF
        IF (TRCINT.AND.NAOTB.GT.0)
     .      WRITE (iunout,*) 'I,NAOTS(IAOT),NAOTT(IAOT)'
        DO 50 IAOT=1,NAOTB
          READ (IUNIN,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          IF (TRCINT) THEN
            WRITE (iunout,'(6I6)') I,NAOTS(IAOT),NAOTT(IAOT)
          ENDIF
50      CONTINUE

        READ (IUNIN,'(6E12.4)') ZMFPTHI, SIGMX, TDGTEMX
      ENDIF
C
C READING BLOCK 14 FROM FORMATTED INPUT FILE (IUNIN) FINISHED
C
C SAVE SOME MORE INPUT DATA FOR SHORT CYCLE ON COMMON CCOUPL
      LNLPLG=NLPLG
      LNLDRF=NLDRFT
      LTRCFL=TRCFLE
      NSTRI=NSTRAI
      DO 60 ISTRA=1,NSTRAI
        LNLVOL(ISTRA)=NLVOL(ISTRA)
60    CONTINUE
      NMODEI=NMODE
      NFILNN=NFILEN
C
C  DEFINE ADDITIONAL TALLIES FOR COUPLING (UPDATED IN SUBR. UPTCOP
C                                              AND IN SUBR. COLLIDE)
      NCPVI=max(NPLSI,13)
      IF (NCPVI.GT.NCPV) THEN
        WRITE (iunout,*) 'FROM INTERFACING SUBROUTINE INFCOP: '
        CALL MASPRM('NCPV',4,NCPV,'NCPVI',5,NCPVI,IERROR)
        CALL EXIT(1)
      ENDIF

      ICPVE(1)=3
      ICPRC(1)=1
      TXTTAL(1,NTALM)=
     .  'TOTAL SAMPLED EMISSION (copv(1)=eppht ?)                    '
      TXTSPC(1,NTALM)='                          '
      TXTUNT(1,NTALM)='W/CM**3                   '

      ICPVE(2)=3
      ICPRC(2)=1
      TXTTAL(2,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION (copv(2))                          '
      TXTSPC(2,NTALM)='                          '
      TXTUNT(2,NTALM)='W/CM**3                   '

      ICPVE(3)=3
      ICPRC(3)=1
      TXTTAL(3,NTALM)=
     .  'EMISSION, THICK FRACTION (copv(3))                          '
      TXTSPC(3,NTALM)='                          '
      TXTUNT(3,NTALM)='W/CM**3                   '

      ICPVE(4)=3
      ICPRC(4)=1
      TXTTAL(4,NTALM)=
     .  'ABSORPTION THICK FRACTION (copv(4))                         '
      TXTSPC(4,NTALM)='                          '
      TXTUNT(4,NTALM)='W/CM**3                   '

      ICPVE(5)=3
      ICPRC(5)=1
      TXTTAL(5,NTALM)=
     .  'EMISSION THIN FRACTION (copv(5))                            '
      TXTSPC(5,NTALM)='                          '
      TXTUNT(5,NTALM)='W/CM**3                   '

      ICPVE(6)=3
      ICPRC(6)=1
      TXTTAL(6,NTALM)=
     .  'ABSORPTION THIN FRACTION (copv(6))                          '
      TXTSPC(6,NTALM)='                          '
      TXTUNT(6,NTALM)='W/CM**3                   '

      ICPVE(7)=1
      ICPRC(7)=1
      TXTTAL(7,NTALM)=
     .  'I_x (copv(7))                                               '
      TXTSPC(7,NTALM)='                          '
      TXTUNT(7,NTALM)='AMP                       '

      ICPVE(8)=1
      ICPRC(8)=1
      TXTTAL(8,NTALM)=
     .  'I_y (copv(8))                                               '
      TXTSPC(8,NTALM)='                          '
      TXTUNT(8,NTALM)='AMP                       '

      ICPVE(9)=1
      ICPRC(9)=1
      TXTTAL(9,NTALM)=
     .  'I_z (copv(9))                                               '
      TXTSPC(9,NTALM)='                          '
      TXTUNT(9,NTALM)='AMP                       '

      ICPVE(10)=0
      ICPRC(10)=0
      TXTTAL(10,NTALM)=
     .  'Kappa   (copv(10))                                          '
      TXTSPC(10,NTALM)='                          '
      TXTUNT(10,NTALM)='AMP                       '

      ICPVE(11)=0
      ICPRC(11)=0
      TXTTAL(11,NTALM)=
     .  'Laplace Operator (copv(11))                                 '
      TXTSPC(11,NTALM)='                          '
      TXTUNT(11,NTALM)='AMP                       '

      ICPVE(12)=3
      ICPRC(12)=1
      TXTTAL(12,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION BY COLLISION ESTIMATOR (copv(12))  '
      TXTSPC(12,NTALM)='                          '
      TXTUNT(12,NTALM)='W/CM**3                   '

      ICPVE(13)=3
      ICPRC(13)=1
      TXTTAL(13,NTALM)=
     .  'TOTAL SAMPLED ABSORPTION BY TRACKLENGTH ESTIMATOR (copv(13))'
      TXTSPC(13,NTALM)='                          '
      TXTUNT(13,NTALM)='W/CM**3                   '
C
      OPEN (UNIT=29,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 29
C
      OPEN (UNIT=30,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 30
C
C  READ IN DATA TO SET UP GEOMETRY FOR NEUTRAL GAS TRANSPORT CODE
C  STATEMENT NUMBER 1000 ---> 1999
C
C  AT PRESENT THE DATA COME FROM THE FILE FT30
C  THIS PART WILL HAVE TO BE MODIFIED AS SOON AS BRAAMS PROVIDES
C  CELL VERTICES AND CUT DESCRIBTION
C
1000  CONTINUE
C
      INMP1I=0
      INMP2I=0
      INMP3I=0
      NTET = 0
      NCOOR = 0
      NTBAR = 0
      NTSEITE = 0
      ncltet = 0
      itethand=1
      ALLOCATE (COORTET(NCOORD))
      DO I=1,NCOORD
        NULLIFY(COORTET(I)%PTET)
      END DO

      nr1p2 = nr1st
      np2t3 = np2nd
      nr1ori = nr1st
      np2ori = np2nd
      nt3ori = nt3rd
      IF (SUM(NOSTEP(1:NTARGI)) > 0) THEN
        NGITT = NQUAD*8
      ELSE
        NGITT = 1
      END IF

      CALL ALLOC_CSTEP

      RRSTEP(1:NSTEP,1) = 0.D0
      ALLOCATE(NRWL(NTARGI))
      NRWL=0

      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'NODAL COORDINATES') == 0)
         READ (30,'(A72)') ZEILE
      END DO

      DO IC=1,NCOORD
         READ (30,*) IDU,XTETRA(IC),YTETRA(IC),ZTETRA(IC)
      END DO
      NCOOR=NCOORD

!  CONVERT TO CM 
      XTETRA(1:NCOORD)=XTETRA(1:NCOORD)*0.1_DP
      YTETRA(1:NCOORD)=YTETRA(1:NCOORD)*0.1_DP
      ZTETRA(1:NCOORD)=ZTETRA(1:NCOORD)*0.1_DP

      write (iunout,*) ' xmin = ',minval(xtetra(1:ncoor))
      write (iunout,*) ' xmax = ',maxval(xtetra(1:ncoor))
      write (iunout,*) ' ymin = ',minval(ytetra(1:ncoor))
      write (iunout,*) ' ymax = ',maxval(ytetra(1:ncoor))
      write (iunout,*) ' zmin = ',minval(ztetra(1:ncoor))
      write (iunout,*) ' zmax = ',maxval(ztetra(1:ncoor))

      ALLOCATE (NODBRICK(27,NBRICK))
      ALLOCATE (NODQUAD(9,NQUAD))
      ALLOCATE (NQDSTS(NQUAD))
      ALLOCATE (NQDSTP(NQUAD))
      ALLOCATE (LUSED(NQUAD))
      NODBRICK = 0
      NODQUAD = 0
      NQDSTS = 0
      NQDSTP = 0
      LUSED = .FALSE.

      NBR = 0
      NQU = 0

      READ (30,'(A72)') ZEILE
      DO WHILE (INDEX(ZEILE,'ELEMENT GROUPS') == 0)
         READ (30,'(A72)') ZEILE
      END DO

      DO
        READ (30,'(A72)',END=9) ZEILE
        I1 = INDEX(ZEILE,'ELEMENTS:')
        DO WHILE (I1 == 0)
          READ (30,'(A72)',END=9) ZEILE
          I1 = INDEX(ZEILE,'ELEMENTS:')
        END DO
        I2 = INDEX(ZEILE,'NODES:')
        I3 = INDEX(ZEILE,'GEOMETRY:')
        READ (ZEILE(I2+7:I3-1),*) NOD

        READ (ZEILE(I1+10:I2),*) NTT

        READ (30,'(A72)') ZEILE
        I1 = INDEX(ZEILE,'ENTITY NAME:')
        DO WHILE (I1 == 0)
          READ (30,'(A72)') ZEILE
          I1 = INDEX(ZEILE,'ENTITY NAME:')
        END DO
        
        LANF=i1+12+verify(zeile(i1+12:),' ')-1
        LEND=verify(zeile,' ',.true.)
        NAMENT = ZEILE(lanf:lend)
        write (iunout,*) ' nament: ',nament
        CALL UPPERCASE(NAMENT)
        ISTS=0
        ISTEP=0
        DO IT=1,NTARGI
          IF (ENTITY(IT) == NAMENT) THEN
            ISTS=NOSTS(IT)
            ISTEP=NOSTEP(IT)
          END IF
        END DO
        
        IF (NOD == 27) THEN
! READ COORDINATE NUMBERS OF BRICKS
          DO IEL=1,NTT
            NBR = NBR + 1
            READ (30,*) IDU,(NODBRICK(J,NBR),J=1,NOD)
          END DO
        ELSE IF (NOD == 9) THEN
! READ COORDINATE NUMBERS OF QUADRILATERALS
          DO IEL=1,NTT
            NQU = NQU + 1
            READ (30,*) IDU,(NODQUAD(J,NQU),J=1,NOD)
!pb            CALL BUBBLE (NODQUAD(1:NOD,NQU),NOD)
            NQDSTS(NQU) = NLIM+ISTS
            NQDSTP(NQU) = ISTEP
          END DO
        ELSE
          EXIT
        END IF
         
      END DO

 9    CONTINUE

      NTACT=0
      DO I=1,NBRICK
        INDCO = NODBRICK(:,I)
C  SET UP TETRAHEDRONS
        CALL MAKE_TETRA_48 (INDCO)

        ncell=NCELL+1
        if (nrtal.ne.nrad) ncltAL(ntact+1:ntet) = ncell

        IMATCH=0
        DO IQ=1,NQUAD

!  check side 1 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(11) == NODQUAD(9,IQ)) THEN
            INDF1 = (/ INDCO(21), INDCO(20), INDCO(19), 
     .                 INDCO(10), INDCO(1),  INDCO(2),
     .                 INDCO(3),  INDCO(12), INDCO(11) /)
            CALL BUBBLE (INDF1,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF1 == NSORQUAD)) THEN
              INMTIT(1,NTACT+1) = NQDSTS(IQ)
              INMTIT(1,NTACT+2) = NQDSTS(IQ)
              INMTIT(1,NTACT+3) = NQDSTS(IQ)
              INMTIT(1,NTACT+4) = NQDSTS(IQ)
              INMTIT(1,NTACT+5) = NQDSTS(IQ)
              INMTIT(1,NTACT+6) = NQDSTS(IQ)
              INMTIT(1,NTACT+7) = NQDSTS(IQ)
              INMTIT(1,NTACT+8) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+1,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+2,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+3,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+4,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+5,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+6,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+7,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+8,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 2 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(15) == NODQUAD(9,IQ)) THEN
            INDF2 = (/ INDCO(3),  INDCO(6),  INDCO(9), 
     .                 INDCO(18), INDCO(27), INDCO(24),
     .                 INDCO(21), INDCO(12), INDCO(15) /)
            CALL BUBBLE (INDF2,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF2 == NSORQUAD)) THEN
              INMTIT(1,NTACT+9) = NQDSTS(IQ)  
              INMTIT(1,NTACT+10) = NQDSTS(IQ) 
              INMTIT(1,NTACT+11) = NQDSTS(IQ)
              INMTIT(1,NTACT+12) = NQDSTS(IQ)
              INMTIT(1,NTACT+13) = NQDSTS(IQ)
              INMTIT(1,NTACT+14) = NQDSTS(IQ)
              INMTIT(1,NTACT+15) = NQDSTS(IQ)
              INMTIT(1,NTACT+16) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+9,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+10,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+11,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+12,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+13,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+14,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+15,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+16,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 3 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(17) == NODQUAD(9,IQ)) THEN
            INDF3 = (/ INDCO(9),  INDCO(8),  INDCO(7), 
     .                 INDCO(16), INDCO(25), INDCO(26),
     .                 INDCO(27), INDCO(18), INDCO(17) /)
            CALL BUBBLE (INDF3,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF3 == NSORQUAD)) THEN
              INMTIT(1,NTACT+17) = NQDSTS(IQ)
              INMTIT(1,NTACT+18) = NQDSTS(IQ)
              INMTIT(1,NTACT+19) = NQDSTS(IQ)
              INMTIT(1,NTACT+20) = NQDSTS(IQ)
              INMTIT(1,NTACT+21) = NQDSTS(IQ)
              INMTIT(1,NTACT+22) = NQDSTS(IQ)
              INMTIT(1,NTACT+23) = NQDSTS(IQ)
              INMTIT(1,NTACT+24) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+17,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+18,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+19,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+20,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+21,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+22,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+23,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+24,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 4 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(13) == NODQUAD(9,IQ)) THEN
            INDF4 = (/ INDCO(7),  INDCO(4),  INDCO(1), 
     .                 INDCO(10), INDCO(19), INDCO(22),
     .                 INDCO(25), INDCO(16), INDCO(13) /)
            CALL BUBBLE (INDF4,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF4 == NSORQUAD)) THEN
              INMTIT(1,NTACT+25) = NQDSTS(IQ) ! side 5
              INMTIT(1,NTACT+26) = NQDSTS(IQ) ! side 6
              INMTIT(1,NTACT+27) = NQDSTS(IQ)
              INMTIT(1,NTACT+28) = NQDSTS(IQ)
              INMTIT(1,NTACT+29) = NQDSTS(IQ)
              INMTIT(1,NTACT+30) = NQDSTS(IQ)
              INMTIT(1,NTACT+31) = NQDSTS(IQ)
              INMTIT(1,NTACT+32) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+25,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+26,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+27,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+28,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+29,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+30,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+31,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+32,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 5 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(5) == NODQUAD(9,IQ)) THEN
            INDF5 = (/ INDCO(3),  INDCO(2),  INDCO(1), 
     .                 INDCO(4),  INDCO(7),  INDCO(8),
     .                 INDCO(9),  INDCO(6),  INDCO(5) /)
            CALL BUBBLE (INDF5,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF5 == NSORQUAD)) THEN
              INMTIT(1,NTACT+33) = NQDSTS(IQ) 
              INMTIT(1,NTACT+34) = NQDSTS(IQ) 
              INMTIT(1,NTACT+35) = NQDSTS(IQ)
              INMTIT(1,NTACT+36) = NQDSTS(IQ)
              INMTIT(1,NTACT+37) = NQDSTS(IQ)
              INMTIT(1,NTACT+38) = NQDSTS(IQ)
              INMTIT(1,NTACT+39) = NQDSTS(IQ)
              INMTIT(1,NTACT+40) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+33,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+34,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+35,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+36,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+37,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+38,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+39,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+40,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  check side 6 of brick for match with covering quadrangle iq
!  center of brick side and center of quadrangle should match 
!  independent of sequence of points 
          IF (INDCO(23) == NODQUAD(9,IQ)) THEN
            INDF6 = (/ INDCO(19), INDCO(20), INDCO(21), 
     .                 INDCO(24), INDCO(27), INDCO(26),
     .                 INDCO(25), INDCO(22), INDCO(23) /)
            CALL BUBBLE (INDF6,9)
            NSORQUAD = NODQUAD(1:9,IQ)
            CALL BUBBLE (NSORQUAD,9)
            IF (ALL(INDF6 == NSORQUAD)) THEN
              INMTIT(1,NTACT+41) = NQDSTS(IQ)
              INMTIT(1,NTACT+42) = NQDSTS(IQ)
              INMTIT(1,NTACT+43) = NQDSTS(IQ)
              INMTIT(1,NTACT+44) = NQDSTS(IQ)
              INMTIT(1,NTACT+45) = NQDSTS(IQ)
              INMTIT(1,NTACT+46) = NQDSTS(IQ)
              INMTIT(1,NTACT+47) = NQDSTS(IQ)
              INMTIT(1,NTACT+48) = NQDSTS(IQ)
              ISTEP = NQDSTP(IQ)
              IF (ISTEP /= 0) THEN
                CALL TET_STEP (ISTEP,NTACT+41,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+42,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+43,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+44,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+45,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+46,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+47,1,NRWL(ISTEP))
                CALL TET_STEP (ISTEP,NTACT+48,1,NRWL(ISTEP))
              END IF
              LUSED(IQ) = .TRUE.
              IMATCH = IMATCH + 1
            END IF
          END IF

!  more than 6 matches are impossible for a brick
!  so move to next brick
          IF (IMATCH == 6) EXIT
        END DO 

        ntact=ntet

      END DO

      WRITE (iunout,*) ' NUMBER OF SIDES OF TETRAHEDRON ON ',
     .            'NON-DEFAULT SURFACES ',COUNT(INMTIT>0)
      WRITE (iunout,*) ' NQUAD*8 ',NQUAD*8

      IF (.NOT.ALL(LUSED)) THEN
        WRITE (iunout,*) 
     .    ' SURFACES WITHOUT MATCHING BRICK FOUND IN IF0COP '
        DO IQ=1,NQUAD
          IF (.NOT.LUSED(IQ)) WRITE (iunout,*) IQ
        END DO
        CALL EXIT(1)
      END IF

!PB [1-sqrt(3/5)]-1,  0,  1- [....] 
      GK(1) = (1._DP - SQRT(3._DP/5._DP)) -1._DP
      GK(2) = 0._DP
      GK(3) = 1._DP - (1._DP - SQRT(3._DP/5._DP))

      DO I=1,1
        J = 0
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,I)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,I)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,I)), K=1,27) /)
        DO IR=1,3
           DO IS=1,3
              DO IT=1,3
                J=J+1
                CALL TRIQUAINT(TQVEC,GK(IR),GK(IS),GK(IT))
                XGAU=SUM(TQVEC*XVEC)*10._dp
                YGAU=SUM(TQVEC*YVEC)*10._dp
                ZGAU=SUM(TQVEC*ZVEC)*10._dp
                WRITE (iunout,'(2I6,3i3,3E20.10)') I,J,ir,is,it,
     .                XGAU,YGAU,ZGAU
              END DO
           END DO
        END DO
      END DO

!  CALCULATE PARTIAL DERIVATIVES DR/DX, DR/DY, DR/DZ, DS/DX, DS/DX, DS/DZ,
!  DT/DX, DT/DY AND DT/DZ FOR EACH BRICK (HEXAHEDRON) AT THE CENTER POINT

      ALLOCATE (DRSTDXYZ(3,3,NBRICK))
      CALL DFTRIQUA(DFVEC, 0._DP, 0._DP, 0._DP)

      DO I=1,NBRICK
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,I)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,I)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,I)), K=1,27) /)
        AJ(1,1) = SUM(XVEC*DFVEC(1,:))
        AJ(1,2) = SUM(YVEC*DFVEC(1,:))
        AJ(1,3) = SUM(ZVEC*DFVEC(1,:))
        AJ(2,1) = SUM(XVEC*DFVEC(2,:))
        AJ(2,2) = SUM(YVEC*DFVEC(2,:))
        AJ(2,3) = SUM(ZVEC*DFVEC(2,:))
        AJ(3,1) = SUM(XVEC*DFVEC(3,:))
        AJ(3,2) = SUM(YVEC*DFVEC(3,:))
        AJ(3,3) = SUM(ZVEC*DFVEC(3,:))

        AJM1T(1,1) = AJ(2,2)*AJ(3,3) - AJ(2,3)*AJ(3,2)
        AJM1T(1,2) = AJ(2,3)*AJ(3,1) - AJ(2,1)*AJ(3,3)
        AJM1T(1,3) = AJ(2,1)*AJ(3,2) - AJ(3,1)*AJ(2,2)
        AJM1T(2,1) = AJ(3,2)*AJ(1,3) - AJ(1,2)*AJ(3,3)
        AJM1T(2,2) = AJ(3,3)*AJ(1,1) - AJ(3,1)*AJ(1,3)
        AJM1T(2,3) = AJ(3,1)*AJ(1,2) - AJ(3,2)*AJ(1,1)
        AJM1T(3,1) = AJ(1,2)*AJ(2,3) - AJ(1,3)*AJ(2,2)
        AJM1T(3,2) = AJ(1,3)*AJ(2,1) - AJ(2,3)*AJ(1,1)
        AJM1T(3,3) = AJ(1,1)*AJ(2,2) - AJ(1,2)*AJ(2,1)

        DETI = 1._DP / (AJ(1,1)*AJ(2,2)*AJ(3,3) 
     .                + AJ(1,2)*AJ(2,3)*AJ(3,1)   
     .                + AJ(1,3)*AJ(2,1)*AJ(3,2)   
     .                - AJ(3,1)*AJ(2,2)*AJ(1,3)   
     .                - AJ(3,2)*AJ(2,3)*AJ(1,1)   
     .                - AJ(3,3)*AJ(2,1)*AJ(1,2))  

        AJM1 = TRANSPOSE(AJM1T) * DETI

        E = matmul(aj,ajm1)

        DRSTDXYZ(:,:,I) = AJM1
        if ((i==2333).or.(i==2340)) then
           write (iunout,*) ' i = ',i
           write (iunout,*) ' koordinaten '
           do j=1,27
              write (iunout,'(i6,4es12.4)') j,xvec(j),yvec(j),zvec(j)
           end do
           write (iunout,*) ' aj '
           write (iunout,'(3es12.4)') aj(1,1),aj(1,2),aj(1,3)
           write (iunout,'(3es12.4)') aj(2,1),aj(2,2),aj(2,3)
           write (iunout,'(3es12.4)') aj(3,1),aj(3,2),aj(3,3)
           write (iunout,*) ' ajm1 '
           write (iunout,'(3es12.4)') ajm1(1,1),ajm1(1,2),ajm1(1,3)
           write (iunout,'(3es12.4)') ajm1(2,1),ajm1(2,2),ajm1(2,3)
           write (iunout,'(3es12.4)') ajm1(3,1),ajm1(3,2),ajm1(3,3)
           write (iunout,*) ' determinate(aj) ',deti
!pb           if (i==2340) stop ' infcop '
        end if
      END DO


      DEALLOCATE (NODQUAD)
      DEALLOCATE (NQDSTS)
      DEALLOCATE (NQDSTP)
      DEALLOCATE (LUSED)

C
C
C  TRANSFER FLAGS
C
      NAINI=NAINB
C
C
      NLPLG=.FALSE.
      NLTET=.TRUE.
      LEVGEO=5
      NR1ST=NTET+1
      NLPOL=.FALSE.
      NP2ND=1
      NLTOR=.FALSE.
      NT3RD=1
      NR1TAL=NCLTAL(NTET)+1
      NP2TAL=NP2ND
      NT3TAL=NT3RD

      DO IN=NR1ST,NR1ST+NRADD
        NCLTAL(IN)=NCLTAL(NTET)+IN-NTET
      ENDDO

      CALL LEER(2)
      CALL HEADNG(' CASE REDEFINED IN COUPLE_TETRA: ',32)
      WRITE (iunout,*) 'NLPLG,NLFEM ',NLPLG,NLFEM
      WRITE (iunout,*) 'NLPOL       ',NLPOL
      WRITE (iunout,*) 'NR1ST,NP2ND ',NR1ST,NP2ND
      CALL LEER(2)
CTRIG E
C
      RETURN
C
C   GEOMETRY DEFINITION PART FINISHED
C
      ENTRY IF1COP
C
C   NOW READ THE PLASMA STATE GIVEN BY BRAAMS
C   AT PRESENT THE DATA COME FROM THE FILE FT31
C   FURTHERMORE: SCALING TO EIRENE UNITS AND INDEX MAPPING
C   STATEMENT NO. 2000 ---> 2999
C
C  IN CASE OF "SHORT CYCLE" THE PLASMA STATE IS TRANSFERRED VIA COMMON
C
      LSHORT=.FALSE.
      CALL LEER(1)
      WRITE (iunout,*) 'IF1COP CALLED '
      IF (NLPLAS) WRITE (iunout,*) 'PLASMA DATA EXPECTED ON BRAEIR'
      IF (.NOT.NLPLAS) 
     .  WRITE (iunout,*) 'PLASMA DATA EXPECTED ON FORT.31'
C  SKIP READING PLASMA, IF NLPLAS
      IF(NLPLAS) GOTO 2100
C
      GOTO 99991
C
C  IN CASE OF "SHORT CYCLE" OR TIME DEP. MODE
C  THE PLASMA STATE IS TRANSFERRED VIA COMMON
C  ONLY SCALING TO EIRENE UNITS AND INDEX MAPPING NEEDS TO BE DONE HERE
C
      ENTRY INTER1
      LSHORT=.TRUE.
      GOTO 2100
C
99991 CONTINUE
C
C
C  TRANSFER PROFILES
C
      IF (.NOT.(INDPRO(1).EQ.9.OR.INDPRO(2).EQ.9.OR.INDPRO(3).EQ.9.OR.
     .          INDPRO(4).EQ.9)) RETURN
C
C
      IF (NFLA.GT.NFL) THEN
        WRITE (iunout,*) ' PARAMETER ERROR DETECTED IN INFCOP '
        WRITE (iunout,*) ' NFLA MUST BE <= NFL'
        WRITE (iunout,*) ' NFLA,NFL = ',NFLA,NFL
        CALL EXIT(1)
      ENDIF

2100  CONTINUE
C
C  RESET 2D ARRAYS ONTO 1D EIRENE ARRAYS, RESCALE TO EIRENE UNITS
C  AND CONVERT BRAAMS VECTORS INTO CARTHESIAN EIRENE VECTORS
C
C  UNITS CONVERSION FACTORS
      T=1._DP/11600._DP
      V=1.
      VL=1.
CTRIG A
C  VACCUM DATA NEEDED FOR REGION OUTSIDE B2-MESH
      TVAC=0.02
      DVAC=1.E2_DP
      VVAC=0.
      BVAC=1.
CTRIG E
      DO 2105 IPLS=1,NPLSI
        D(IPLS)=FCTE(IPLS)*1.E-6_DP
        FL(IPLS)=FCTE(IPLS)
2105  CONTINUE
C
      ALLOCATE (NOSPEC(NPLSI))
      ALLOCATE (LUSED(NPLSI))

      ALLOCATE (DTEDX(NRTAL))
      ALLOCATE (DTEDY(NRTAL))
      ALLOCATE (DTEDZ(NRTAL))
      ALLOCATE (xlapla(NRTAL))

      OPEN (UNIT=31,ACCESS='SEQUENTIAL',FORM='FORMATTED')
      REWIND 31
      READ (31,'(A1000)') LINE

      ICOLUMN = 0
      IANF = 1
      NOSPEC = 0
      FORM='(A       )'
      LUSED=.FALSE.
      DO IPLS=1,NPLSI
        IF (LEN_TRIM(CDENMODEL(IPLS)) > 0) LUSED(IPLS)=.TRUE.
      END DO
      DO 
!  POSITION OF NEXT NON-BLANK CHARACTER
        I1 = VERIFY(LINE(IANF:),' ')
        IF (I1 == 0) EXIT
        I2 = INDEX(LINE(IANF+I1:),' ')
        ICOLUMN=ICOLUMN+1
        IF (ICOLUMN > 6) THEN
          IF (I2 < 10) THEN
            WRITE (FORM(3:3),'(I1)') I2
          ELSE
            WRITE (FORM(3:4),'(I2)') I2
          END IF
          SPECNAME=REPEAT(' ',72)
          READ(LINE(IANF+I1-1:IANF+I1+I2-2),FORM) SPECNAME(1:i2)
          SPECNAME=TRIM(SPECNAME)
          CALL UPPERCASE (SPECNAME)
          DO IPLS=1,NPLSI
            IF (SPECNAME == TEXTS(NSPAMI+IPLS)) THEN
              IF (NOSPEC(IPLS) == 0) THEN
                NOSPEC(IPLS) = ICOLUMN-6
                LUSED(IPLS)=.TRUE.
              ELSE
                WRITE (iunout,*) 
     .           ' ERROR! SPECIES ',SPECNAME,' FOUND TWICE '
              END IF
              EXIT
            END IF
          END DO
        END IF
        IANF=IANF+I1+I2-1
      END DO
      
      IF (.NOT. ALL(LUSED)) THEN
        WRITE (iunout,*) ' ERROR DETECTED IN INFCOP! '
        WRITE (iunout,*) 
     .    ' BULK PARTICLE DEFINITION IN EIRENE INPUT DOES NOT'
        WRITE (iunout,*) ' MATCH SPECIFICATIONS IN FIDAP PLASMA FILE '
        WRITE (iunout,*) ' CALCULATION ABANDONNED '
        CALL EXIT(1)
      END IF

      NREAD=maxval(NOSPEC(1:nplsi))
      ALLOCATE (TEF(NCOOR))
      ALLOCATE (DENF(0:MAX(NPLSI,NREAD),NCOOR))

      DO IC=1,NCOOR
        READ (31,'(A1000)') LINE
        READ (LINE,*) ICO,XCO,YCO,ZCO,TEF(IC),
     .                (DENF(IPLS,IC),IPLS=0,NREAD)
      END DO

      TEIN=TVAC
      TIIN=TVAC
      DIIN=DVAC
      DO ITET=1,NTET
        DIST1=1._DP/SQRT((XTETRA(NTECK(1,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(1,ITET))-YTCEN(ITET))**2)
        DIST2=1._DP/SQRT((XTETRA(NTECK(2,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(2,ITET))-YTCEN(ITET))**2)
        DIST3=1._DP/SQRT((XTETRA(NTECK(3,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(3,ITET))-YTCEN(ITET))**2)
        DIST4=1._DP/SQRT((XTETRA(NTECK(4,ITET))-XTCEN(ITET))**2+
     .                   (YTETRA(NTECK(4,ITET))-YTCEN(ITET))**2)
        DIST=DIST1+DIST2+DIST3+DIST4
        TEIN(ITET) = T*(TEF(NTECK(1,ITET))*DIST1 + 
     .                  TEF(NTECK(2,ITET))*DIST2 +
     .                  TEF(NTECK(3,ITET))*DIST3 +
     .                  TEF(NTECK(4,ITET))*DIST4)/DIST
        DO IPLS=1,NPLSI
          IPL=NOSPEC(IPLS)
          IF (IPL == 0) CYCLE
          IPLSTI=MPLSTI(IPLS)
          TIIN(IPLSTI,ITET) = TEIN(ITET)
          DIIN(IPLS,ITET) = D(IPLS)*
     .                     (DENF(IPL,NTECK(1,ITET))*DIST1 + 
     .                      DENF(IPL,NTECK(2,ITET))*DIST2 +
     .                      DENF(IPL,NTECK(3,ITET))*DIST3 +
     .                      DENF(IPL,NTECK(4,ITET))*DIST4)/DIST
        END DO
      END DO

      write (iunout,*) ' temax = ',maxval(tein(1:NTET))
      write (iunout,*) ' temin = ',minval(tein(1:NTET))
!  koordinaten der gaussknoten in r,s,t
!PB [1-sqrt(3/5)]-1,  0,  1- [....] 
      GK(1) = (1._DP - SQRT(3._DP/5._DP)) -1._DP
      GK(2) = 0._DP
      GK(3) = 1._DP - (1._DP - SQRT(3._DP/5._DP))

        

      dtedx = 0._dp
      dtedy = 0._dp
      dtedz = 0._dp
      xlapla = 0._dp
!   eckpunkte des bricks (27 pro brick)
      do ir=1,nr1tal-1
        XVEC(1:27) = (/ (XTETRA(NODBRICK(K,IR)), K=1,27) /)
        YVEC(1:27) = (/ (YTETRA(NODBRICK(K,IR)), K=1,27) /)
        ZVEC(1:27) = (/ (ZTETRA(NODBRICK(K,IR)), K=1,27) /)
        TEVEC(1:27) = T*(/ (TEF(NODBRICK(K,IR)), K=1,27) /)
!   calculate gradient T in each cell (centre of brick)
        call gradf (0._dp, 0._dp, 0._dp, xvec, yvec, zvec, tevec,
     .              dtedx(ir), dtedy(ir), dtedz(ir))
!   gradient length in each cell, for comparison with photon mean free path
        TEDTEDX(IR) =  TEVEC(14) / DTEDX(IR)
        TEDTEDY(IR) =  TEVEC(14) / DTEDY(IR)
        TEDTEDZ(IR) =  TEVEC(14) / DTEDZ(IR)

!   calculate laplace T in each cell (centre of brick)
        xlapla(ir) = flaplace(0._dp, 0._dp, 0._dp, 
     .                         xvec, yvec, zvec, tevec)
      end do

      DEALLOCATE (TEF)
      DEALLOCATE (DENF)
      DEALLOCATE (NOSPEC)
      DEALLOCATE (LUSED)

C
C
      RETURN
C
2999  CONTINUE
C
C  PLASMA PROFILES ARE NOW READ IN
C
      ENTRY IF2COP(ITARG)
      IF (ITARG.GT.NTARGI) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'SOURCE DATA FOR STRATUM ISTRA= ',ITARG
        WRITE (iunout,*) 
     .    'CANNOT BE DEFINED IN IF2COP. CHANGE INDSRC(ISTRA)'
        CALL LEER(1)
        RETURN
      ENDIF
C
C  NEXT DEFINE FLUXES, TEMPERATURES AND VELOCITIES AT THE TARGETS
C  (FLUXES IN AMP/(CM ALONG TARGET), TEMPERATURES IN EV, VELOCITIES IN CM/SEC)
C   FNIXB*FL (FNIYB*FL) ARE GIVEN IN AMP
C  STATEMENT NO 3000 ---> 3999
C
3000  CONTINUE
C
      
      TESTEP(ITARG,:)=0.D0
      TISTEP(:,ITARG,:)=0.D0
      DISTEP(:,ITARG,:)=0.D0
      VXSTEP(:,ITARG,:)=0.D0
      VYSTEP(:,ITARG,:)=0.D0
      VZSTEP(:,ITARG,:)=0.D0
      FLSTEP(:,ITARG,:)=0.D0
      DO IG=1,NRWL(ITARG)
        IN=IRSTEP(ITARG,IG)
        TESTEP(ITARG,IG)=TEIN(IN)
        TISTEP(1:NPLSTI,ITARG,IG)=TIIN(1:NPLSTI,IN)
        DISTEP(1:NPLSI,ITARG,IG)=DIIN(1:NPLSI,IN)
        VXSTEP(1:NPLSV,ITARG,IG)=VXIN(1:NPLSV,IN)
        VYSTEP(1:NPLSV,ITARG,IG)=VYIN(1:NPLSV,IN)
        VZSTEP(1:NPLSV,ITARG,IG)=VZIN(1:NPLSV,IN)
        DO IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          CSP(IPLS)=CVEL2A*SQRT((TIIN(IPLSTI,IN)+TEIN(IN))/
     .              RMASSP(IPLS))
          FLSTEP(IPLS,ITARG,IG)=ELCHA*0.5*DIIN(IPLS,IN)*CSP(IPLS)
        END DO
      END DO
      nrwl(itarg)=nrwl(itarg)+1
C
      IF (TRCSOU) CALL LEER(2)
C
C  INITALIZE FUNCTION STEP (FOR RANDOM SAMPLING ALONG TARGET)
C  SET SOME SOURCE PARAMETERS EXPLICITLY TO ENFORCE INPUT CONSISTENCY
C
      IIPLS=1
      IEPLS=NPLSI
      FLUX(ITARG)=STEP(IIPLS,IEPLS,NRWL(ITARG),ITARG)
C
      NLPLS(ITARG)=.TRUE.
      NLATM(ITARG)=.FALSE.
      NLMOL(ITARG)=.FALSE.
      NLION(ITARG)=.FALSE.
C
      NLSRF(ITARG)=.TRUE.
      NLPNT(ITARG)=.FALSE.
      NLLNE(ITARG)=.FALSE.
      NLVOL(ITARG)=.FALSE.
      NLCNS(ITARG)=.FALSE.
C
      NSRFSI(ITARG)=1
      INDIM(1,ITARG)=4
      I4=IDEZ(INT(SORLIM(1,ITARG)),4,4)
      SORLIM(1,ITARG)=I4*1000+104
      SORIND(1,ITARG)=ITARG
C  IN CASE INDIM=4: INSOR,INDGRD... ARE REDUNDANT
      NRSOR(1,ITARG)=-1
      NPSOR(1,ITARG)=-1
      IF (INDSRC(ITARG).LT.6) THEN
        WRITE (iunout,*) 'MESSAGE FROM IF2COP: '
        WRITE (iunout,*) 'SOURCE STRENGTH AND SPATIAL DISTRIBUTION FOR '
        WRITE (iunout,*) 'STRATUM ',ISTRA,' MODIFIED.'
        CALL MASR1('FLUX=   ',FLUX(ISTRA))
        WRITE (iunout,*) 'USE STEP FUNCTION ISTEP= ',ITARG,
     .    ' FROM BLOCK 14'
      ENDIF
C
      IF (INDSRC(ITARG).EQ.6) THEN
C  DEFINE SOURCE FOR TARGET RECYCLING STRATUM ITARG
C  ASSUME NOW: ITARG=ISTRA
C  DEFAULTS ARE ALREADY SET IN SUBR. INPUT.
C
        CALL FTCRI(ITARG,CITARG)
        TXTSOU(ITARG)= 'SURFACE RECYCLING SOURCE NO.'//CITARG
        NPTS(ITARG)=NPTC(ITARG,1)
        NINITL(ITARG)=ITARG*1001
        NSPEZ(ITARG)=-1
        SORIFL(1,ITARG)=NIFLG(ITARG,1)
        SORWGT(1,ITARG)=1.
        IF (NIXY(ITARG,1).EQ.1) THEN
C TARGET RECYCLING SOURCE AT POLOIDAL SURFACE NPES
          NEMODS(ITARG)=3
          NAMODS(ITARG)=1
          SORENI(ITARG)=3.
          SORENE(ITARG)=0.5
        ELSEIF (NIXY(ITARG,1).EQ.2) THEN
C WALL RECYCLING SOURCE AT RADIAL SURFACE NPES
          NEMODS(ITARG)=2
          NAMODS(ITARG)=1
          SORENI(ITARG)=2.
          SORENE(ITARG)=0.
        ENDIF
C
C  SORAD1,...: USE POLYGON MESH, IE. SORAD1,... ARE REDUNDANT.
C
C  VELOCITY SPACE DISTRIBUTION
        SORCOS(ITARG)=1.
        SORMAX(ITARG)=0.
C
C
C  DO 2028 LOOP FROM SUBR. INPUT
        THMAX=MAX(0._DP,MIN(PIHA,SORMAX(ITARG)*DEGRAD))
        IF (NAMODS(ITARG).EQ.1) THEN
          RP1=SORCOS(ITARG)+1.
          SORCOS(ITARG)=1./RP1
          IF (ABS(COS(THMAX)).LE.EPS10) THEN
            SORMAX(ITARG)=1.
          ELSE
            SORMAX(ITARG)=1.-COS(THMAX)**RP1
          ENDIF
        ELSEIF (NAMODS(ITARG).EQ.2) THEN
          SORCOS(ITARG)=SORCOS(ITARG)*DEGRAD
          SORMAX(ITARG)=THMAX
        ENDIF
        NLSYMT(0)=NLSYMT(0).AND.NLSYMT(ITARG)
        NLSYMP(0)=NLSYMP(0).AND.NLSYMP(ITARG)
C
      ENDIF
C
C  SOURCE DEFINITION FOR TARGET RECYCLING STRATUM ITARG COMPLETED
C
3999  CONTINUE
C
C  TARGET DATA ARE DEFINED NOW
C
      DO 5000 IG=1,NGITT
        ELSTEP(:,ITARG,IG)=0.
5000  CONTINUE
C
C  COMPUTE EXACT SURFACE ENERGY FLUXES FOR COMPARISON WITH SAMPLED
C  E-FLUX "ETOTP". THIS IS ONLY FOR DIAGNOSTICS PURPOSES
C  E.G. TO CHECK CONSISTENCY OF BOUNDARY CONDITIONS
C  STATEMENT NO. 6000 ---> 6500
C
      IF (.NOT.TRCSOU) GOTO 6500
C
      write (iunout,*) ' exact surface energy flux computation '
      write (iunout,*) ' not yet available of tetrahedrons '
      GOTO 6500

 4812 continue
      EEMAX=0.
      EESHT=0.
C
      DO 6011 IG=1,NRWL(ITARG)-1
        OR=ORI(ITARG,IG)
C  COMPUTE SHEATH POTENTIAL ESHT(ITARG,IG)
C  USE ALL NPLSI SPECIES, NOT JUST IFL=NSPZI,NSPZE
        ESHT(ITARG,IG)=0.D0
        IF (NEMODS(ITARG).EQ.3.OR.NEMODS(ITARG).EQ.5.OR.
     .      NEMODS(ITARG).EQ.7) THEN
          IF (IGSTEP(ITARG,IG).GT.200000) THEN
            ITRI=IRSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-200000
            DO 6005 IPL=1,NPLSI
              IPLV=MPLSV(IPL)
              PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLV,ITARG,IG)+
     .             PTRIY(NPES,ITRI)*VYSTEP(IPLV,ITARG,IG))*OR
              VPZ=VZSTEP(IPLV,ITARG,IG)
              VP(IPL)=SQRT(PM1**2+VPZ**2)
              DI(IPL)=DISTEP(IPL,ITARG,IG)
6005        CONTINUE
            TE=TESTEP(ITARG,IG)
            CUR=0.
            GAMMA=0.
            ESHT(ITARG,IG)=SHEATH(TE,DI,VP,NCHRGP,GAMMA,CUR,NPLSI,
     .                           -ITARG)
          ELSEIF (IGSTEP(ITARG,IG).LT.200000) THEN
            ITRI=IRSTEP(ITARG,IG)
            NPES=IGSTEP(ITARG,IG)-100000
            ESHT(ITARG,IG)= 0
          ENDIF
        ENDIF
C
        IF (IGSTEP(ITARG,IG).LT.200000) GOTO 6010
        ITRI=IRSTEP(ITARG,IG)
        NPES=IGSTEP(ITARG,IG)-200000
        DO 6009 IPLS=1,NPLSI
          IF (FLSTEP(IPLS,ITARG,IG).EQ.0.D0) GOTO 6009
          IPLSTI=MPLSTI(IPLS)
          IPLSV=MPLSV(IPLS)
          VT=SQRT(2.*TISTEP(IPLSTI,ITARG,IG)/BMASS(IPLS))*CVEL2A
C  VELOCITY COMPONENT NORMAL TO TARGET SURFACE
C  I.E., POLOIDAL COMPONENT V-POL
C  ASSUMING ORTHOGONAL TARGET
          PM1=(PTRIX(NPES,ITRI)*VXSTEP(IPLSV,ITARG,IG)+
     .         PTRIY(NPES,ITRI)*VYSTEP(IPLSV,ITARG,IG))*OR
C  VELOCITY COMPONENT PARALLEL TO TARGET SURFACE
C  I.E., RADIAL PLUS TOROIDAL COMPONENT, V-RAD + V-TOR
C  AGAIN: ASSUMING ORTHOGONAL TARGET
          VPX=VXSTEP(IPLSV,ITARG,IG)-PM1*PPLNX(IY,NPES)*OR
          VPY=VYSTEP(IPLSV,ITARG,IG)-PM1*PPLNY(IY,NPES)*OR
          VPZ=VZSTEP(IPLSV,ITARG,IG)-0.
          PN1=SQRT(VPX**2+VPY**2+VPZ**2)
          PERW=0.
          PARW=0.
          IF (VT.GT.0.) THEN
            PERW=PM1/VT
            PARW=PN1/VT
          ENDIF
C
          CS=SQRT((1.*TISTEP(IPLSTI,ITARG,IG)+
     .                TESTEP(ITARG,IG))/BMASS(IPLS))*CVEL2A
C THE MACH NUMBER BOUNDARY CONDITION ONLY AFFECTS THE PARALLEL TO B
C MOMENTUM, I.E., NOT THE RADIAL VELOCITY
          VTEST=SQRT(PM1**2+VPZ**2)
          VTEST=VTEST/(CS+EPS60)
          VR=SQRT(VPX**2+VPY**2)
          WRITE (iunout,*) 'IPLS,ITARG,IG,MACH ',IPLS,ITARG,IG,VTEST
C         WRITE (iunout,*) 'POL., TOR., RAD. ',PM1,VPZ,VR
          CALL LEER(1)
C TARGET ENERGY FLUXES
          DRR=RRSTEP(ITARG,IG+1)-RRSTEP(ITARG,IG)
          IF (NEMODS(ITARG).EQ.1) THEN
            EADD=SORENI(ITARG)
          ELSEIF (NEMODS(ITARG).EQ.2.OR.NEMODS(ITARG).EQ.3) THEN
            EADD=SORENI(ITARG)*TISTEP(IPLSTI,ITARG,IG)+SORENE(ITARG)*
     .             TESTEP(ITARG,IG)
          ELSEIF (NEMODS(ITARG).GE.4) THEN
            PERWI=PERW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            PARWI=PARW/SQRT(BMASS(IPLS)/RMASSP(IPLS))
            EADD=EMAXW(TISTEP(IPLSTI,ITARG,IG),PERWI,PARWI)
          ENDIF
          ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
          ELSTEP(IPLS,ITARG,IG)=ELSTEP(IPLS,ITARG,IG)+ESUM
          EEMAX=EEMAX+ESUM*DRR
C  ADD SHEATH ACCELERATION
          EADD=NCHRGP(IPLS)*ESHT(ITARG,IG)
          ESUM=EADD*FLSTEP(IPLS,ITARG,IG)
          EESHT=EESHT+ESUM*DRR
          ELSTEP(IPLS,ITARG,IG)=ELSTEP(IPLS,ITARG,IG)+ESUM
6009    CONTINUE
        GOTO 6011
6010    CONTINUE
C  TO BE WRITTEN
6011  CONTINUE
C
      CALL LEER(1)
      WRITE (iunout,*) 'TARGET DATA: TARGET NO. ITARG=ISTRA= ',ITARG
      WRITE (iunout,*) 'IG, ARC, P-FLUX, E-FLUX, TE, TI, SHEATH/TE'
      DO 6100 IG=1,NRWL(ITARG)-1
        WRITE (iunout,'(1X,I4,1P,6E11.3)')
     .             IG,RRSTEP(ITARG,IG),FLSTEP(0,ITARG,IG),
     .             ELSTEP(0,ITARG,IG),
     .             TESTEP(ITARG,IG),TISTEP(1,ITARG,IG),
     .             ESHT(ITARG,IG)/(TESTEP(ITARG,IG)+EPS60)
6100  CONTINUE
      WRITE (iunout,'(1X,I4,1P,1E11.3)') NRWL(ITARG),
     .                                 RRSTEP(ITARG,NRWL(ITARG))
      CALL MASR1 ('EEMAX    ',EEMAX)
      CALL MASR1 ('EESHT    ',EESHT)
C
      ETOT=EEMAX+EESHT
      EFLX(ITARG)=EEMAX+EESHT
      WRITE (iunout,*) 'PARTICLE FLUX(IPLS), IPLS=1,NPLSI '
      WRITE (iunout,'(1X,1P,6E12.4)') (FLTOT(ISP,ITARG),ISP=1,NPLSI)
      CALL LEER(1)
      WRITE (iunout,*) 'ENERGY FLUX '
      WRITE (iunout,'(1X,1P,1E12.4)') EFLX(ITARG)
      CALL LEER(2)
C
6300  CONTINUE
C
C
C  SET SOME OTHER DATA SPECIFIC FOR EIRENE CODE REQUIREMENTS
C  STATEMENT NO. 6500 ---> 6999
C
6500  CONTINUE
C
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN IF2COP: NGITT TOO SMALL '
      CALL EXIT(1)
      RETURN
C
C
      ENTRY IF3COP(ISTRAA,ISTRAE,NEW_ITER)
C
C
      WRITE (iunout,*) ' IF3COP IS CALLED, ISTRAA,ISTRAE '
      WRITE (iunout,*) ISTRAA,ISTRAE
      LSHORT=.FALSE.
      LSTOP=.TRUE.
      IFIRST=0
      GOTO 99992
C
      ENTRY INTER3(LSTP,IFRST,ISTRAA,ISTRAE,NEW_ITER)
C
C  ENTRY FOR SHORT CYCLE FROM SUBR. EIRSRT
C
C  IFIRST=0: RESTORE DATA FROM A PREVIOUS EIRENE RUN, SET REFERENCE
C            DATA FOR "STOP-CRITERION" SNIS,SEES,SEIS
C  IFIRST>0: MODIFY SOURCE TERMS ACCORDING TO NEW PLASMA CONDITIONS,
C            COMPARE INTEGRALS WITH SNIS,...., AND DECIDE TO STOP OR
C            CONTINUE SHORT CYCLE (LSTOP)
C
      LSHORT=.TRUE.
      LSTOP=LSTP
      IFIRST=IFRST
C
99992 CONTINUE
C
      DO 10000 ISTRAI=ISTRAA,ISTRAE
C
        IF (XMCP(ISTRAI).LE.1.) GOTO 10000
C
        IF (LSHORT) GOTO 7000
C
C  READ DATA FROM STRATUM NO. ISTRAI BACK INTO WORKING SPACE
C  IF REQUIRED
C
        IF (ISTRAI.EQ.IESTR) THEN
C  NOTHING TO BE DONE
        ELSEIF ((NFILEN.EQ.1.OR.NFILEN.EQ.2).AND.IESTR.NE.ISTRAI) THEN
          IESTR=ISTRAI
          CALL RSTRT(ISTRAI,NSTRAI,NESTM1,NESTM2,NADSPC,
     .               ESTIMV,ESTIMS,ESTIML,
     .               NSDVI1,SDVI1,NSDVI2,SDVI2,
     .               NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .               NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .               NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .               NSIGI_SPC,TRCFLE)
        ELSE
          WRITE (iunout,*) 'ERROR IN INFCOP: STRATUM ISTRAI= ',ISTRAI
          WRITE (iunout,*) 'IS NOT AVAILABLE. EXIT CALLED'
          CALL EXIT(1)
        ENDIF
C
C  DATA TRANSFER BACK FROM EIRENE TO EXTERNAL CODE
C  STATEMENT NO 7000 ---> 7999
C
7000    CONTINUE
C
C  SCALE SURFACE SOURCES PER UNIT FLUX, FOR OTHER SOURCES USE
C  EIRENE SCALINGS
        IF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).NE.0.) THEN
C  FLUX FROM EIRENE TO PLASMA CODE: NEGATIVE
          FLX=-WTOTP(0,ISTRAI)
          FLXI=1./FLX
        ELSEIF (ISTRAI.LE.NTARGI.AND.WTOTP(0,ISTRAI).EQ.0.) THEN
          WRITE (iunout,*) 'NO FLUX FROM STRATUM NO. ISTRAI= ',ISTRAI
          WRITE (iunout,*) 'NO DATA RETURNED FOR THIS STRATUM'
          GOTO 7999
        ELSEIF (ISTRAI.GT.NTARGI) THEN
          FLXI=1.
        ENDIF
C
        IF (.NOT.LSHORT) GOTO 7400
C  SHORT LOOP CORRECTION FINISHED
C
7400    CONTINUE
C
C
C  ADD CONTRIBUTIONS FROM VOLUME RECOMBINATION SOURCE
C
C
        EPEL = 0.D0
        EMIS = 0._DP
        ABSORB = 0._DP
        VOLSUM = 0._DP
        ABSAMP = 0._DP
        AB1SMP = 0._DP
        AB2SMP = 0._DP
        EMSAMP = 0._DP
        EM1SMP = 0._DP
        EM2SMP = 0._DP
        RIX = 0._DP
        RIY = 0._DP
        RIZ = 0._DP
        XKAP = 0._DP
        YKAP = 0._DP
        ZKAP = 0._DP

        IF (NLVOL(ISTRAI)) THEN
C
C  EPEL: CALCULATED EMISSION (NOT SAMPLED)
C
          RECTOT = 0._DP
          DO ISRFSI=1,NSRFSI(ISTRAI)
            ISR=ISRFSI
            ISTEP=SORIND(ISR,ISTRAI)
            DO IPLS=1,NPLSI
              DO IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
                IF ((ISTEP.EQ.0).OR.(ISTEP.EQ.IRRC)) THEN
                  SUMN=0.0
                  SUMEE=0.0
                  DO ITET=1,NTET
                    INC=NCLTAL(ITET)
                    IF (NSTORDR >= NRAD) THEN
                      RECADD=-TABRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                      EEADD=  EELRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                    ELSE
                      RECADD=-FTABRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                      EEADD=  FEELRC1(IRRC,ITET)*DIIN(IPLS,ITET)*ELCHA
                    END IF
                    SUMN=SUMN+RECADD*VOL(ITET)
                    EPEL(INC)=EPEL(INC)+EEADD*VOL(ITET)
                    SUMEE=SUMEE+EEADD*VOL(ITET)
                  END DO         ! ITET
                  RECTOT = RECTOT + SUMN
                  WRITE (iunout,*) 'IPLS,IRRC ',IPLS,IRRC
                  CALL MASR2('SUMN, SUMEE     ',
     .                        SUMN,SUMEE)
                END IF
              END DO             ! IIRC
            END DO               ! IPLS 
          END DO                 ! ISRFSI
          do i=1,nsbox_tal
             epel(i) = epel(i)/voltal(i)
          END DO                 
        END IF

C  NEW VERSION FOR KAPPA_RAD
C  KAPPA FROM  NET SOURCE, THICK FRACTION, AND LAPLACE T 
C  (ANY T WILL DO, ALL ARE THE SAME)

        do i=1,nsbox_tal
           COPV(10,I)=-(COPV(3,I)+copv(4,I))/(xlapla(I)+EPS60)
           COPV(11,I)=xlapla(i)

! if kappa is negative put sources for this cell into thin fraction 
! note:  the variance for the total (sigma_cop(nspvi+1)) is not affected by this
!        but the individual variances sigma_cop(3...6) are not correct anymore.
            LSHIFT=COPV(10,I) < 0._DP
            IF (LSHIFT) THEN
              COPV(5,I)  = COPV(5,I) + COPV(3,I)
              COPV(3,I)  = 0._DP
              COPV(6,I)  = COPV(6,I) + COPV(4,I)
              COPV(4,I)  = 0._DP
!   set radiative vector flux I  of thick fraction to zero in this cell
              COPV(7,I)  = 0._DP
              COPV(8,I)  = 0._DP
              COPV(9,I)  = 0._DP
!   set kappa to zero in this cell
              COPV(10,I) = 0._DP
            END IF
        end do


        do i=3,13
          HELPAR(1:NSBOX_TAL) = COPV(I,1:NSBOX_TAL)
          CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .                 COPVI(I,ISTRAI),
     .                 NR1TAL,NP2TAL,NT3TAL,NBMLT)
          COPV(I,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
        end do


! NEXT:  FILL ADDV, FOR PLOTTING, VARIANCES AND DIAGNOSTICS ARRAYS

!  addv(1):  calculated emission
!  addv(2):  calculated emission - sampled absorption
!  addv(3):  standard deviation (%) for addv(2)
!  addv(4):  sampled emission - sampled absorption
!  addv(5):  standard deviation (%) for sampled emission - sampled absorption
!! NEXT: INTERMEDIATE DIAGNOSTICS,  6-8: CURRENTLY OUT
!!  addv(6):  dT/dx,  and later: div(vec(I)
!!  addv(7):  dT/dy
!!  addv(8):  dT/dz

        ADDV(1,:) = EPEL
!       ADDV(2,:) = EPEL + EPHPHT = EPEL + COPV(2)
        ADDV(2,:) = EPEL + COPV(2,:)

!  total emission is calculated. Standard deviation in absolute units
!  rescale sigma(,..) 1. )to absolute values, then
!                     2.) from total absorption to net emission.
!                     since total emission is an additive constant,
!                     step 2.) is not necessary
        WHERE (ADDV(2,:) /= 0._DP)
          ADDV(3,:) = SIGMA_COP(2,:)/100._DP*ABS(COPV(2,:)) 
        ELSE WHERE
          ADDV(3,:) = 0._DP
        END WHERE
!

!  total emission is sampled. 
        ADDV(4,:) = COPV(3,:) + COPV(4,:) + COPV(5,:) + COPV(6,:)
!  rescale sigma(,..) to absolute values 
        ADDV(5,:) = SIGMA_COP(NCPVI+1,:)/100._DP*ABS(ADDV(4,:))
!
        HELPAR(1:NSBOX_TAL) = ADDV(1,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(1,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(1,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
        
        HELPAR(1:NSBOX_TAL) = ADDV(2,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(2,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(2,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

!  STD DEVIATION OF CELL AVERAGE, ABSOLUTE VALUES
        ADDVI(3,ISTRAI) = SGMS_COP(2)/100._DP*ABS(COPVI(2,ISTRAI))
        
        HELPAR(1:NSBOX_TAL) = ADDV(4,1:NSBOX_TAL)
        CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
     .               ADDVI(4,ISTRAI),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        ADDV(4,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

!  STD DEVIATION OF CELL AVERAGE, ABSOLUTE VALUES
        ADDVI(5,ISTRAI) = SGMS_COP(NCPVI+1)/100._DP*
     .                    ABS(ADDVI(4,ISTRAI))

!       ADDV(6,1:NSBOX_TAL) = DTEDX(1:NSBOX_TAL)
!       ADDV(7,1:NSBOX_TAL) = DTEDY(1:NSBOX_TAL)
!       ADDV(8,1:NSBOX_TAL) = DTEDZ(1:NSBOX_TAL)

!       HELPAR(1:NSBOX_TAL) = ADDV(6,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(6,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(6,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
!       HELPAR(1:NSBOX_TAL) = ADDV(7,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(7,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(7,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)
!       HELPAR(1:NSBOX_TAL) = ADDV(8,1:NSBOX_TAL)
!       CALL INTTAL (HELPAR,VOLTAL,1,1,NSBOX_TAL,
!    .               ADDVI(8,ISTRAI),
!    .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
!       ADDV(8,1:NSBOX_TAL) = HELPAR(1:NSBOX_TAL)

C
        IF (.NOT.LSYMET) GOTO 7500
C
C  SECONDLY SYMMETRISE EIRENE ARRAYS ACCORDING TO SYMMETRY IN MODEL
C
C
C   THIRDLY WRITE EIRENE ARRAYS (1D) ONTO BRAAMS ARRAYS (2D)
C   AND RESCALE TO PROPER UNITS: #/CELL/STRATUM FLUX
C   # STANDS FOR PARTICLES (SNI), MOMENTUM (SMO)
C   AND ENERGY (SEE,SEI)
C
7500    CONTINUE
C
C   NEXT:
C   IF LSHORT: CRITERION TO STOP SHORT CYCLE,
C   IF NOT LSHORT: RESCALE SURFACE SOURCE STRATA
C                  UNITS: # PER UNIT TARGET PLATE FLUX
C
C
C
C   THIRDLY:
C   INDEX MAPPING BACK TO BRAAMS IMPLEMENTATION OF LINDA GEOMETRY
C
C
7700    CONTINUE
C
7999    CONTINUE

c  extrapolate tallies from cell averages to fidap-vertices (27 for each brick)

        DO ITET=1,NTET
          INC=NCLTAL(ITET)
          DO J=1,4
            ipo=nteck(j,itet)
            dist1=1._dp/sqrt((xtetra(ipo)-xtcen(itet))**2 +
     .                       (ytetra(ipo)-ytcen(itet))**2 +
     .                       (ztetra(ipo)-ztcen(itet))**2) 
            volsum(ipo)=volsum(ipo)+dist1
            EMIS(ipo) = EMIS(ipo) + EPEL(NCLTAL(ITET))*dist1
            ABSORB(ipo) = ABSORB(ipo)+
     .                    EPHPHT(NCLTAL(ITET))*dist1
            EMSAMP(ipo) = EMSAMP(IPO) + COPV(1,NCLTAL(ITET))*DIST1
            ABSAMP(ipo) = ABSAMP(IPO) + COPV(2,NCLTAL(ITET))*DIST1
            EM1SMP(ipo) = EM1SMP(IPO) + COPV(3,NCLTAL(ITET))*DIST1
            AB1SMP(ipo) = AB1SMP(IPO) + COPV(4,NCLTAL(ITET))*DIST1
            EM2SMP(ipo) = EM2SMP(IPO) + COPV(5,NCLTAL(ITET))*DIST1
            AB2SMP(ipo) = AB2SMP(IPO) + COPV(6,NCLTAL(ITET))*DIST1
            RIX(ipo)    = RIX(IPO) + COPV(7,NCLTAL(ITET))*DIST1
            RIY(ipo)    = RIY(IPO) + COPV(8,NCLTAL(ITET))*DIST1
            RIZ(ipo)    = RIZ(IPO) + COPV(9,NCLTAL(ITET))*DIST1
!   kappa_rad from grad Te:  out
!   copv(11,11,12) have now another meaning
!pb         xkap(ipo)   = xkap(IPO) + COPV(10,NCLTAL(ITET))*DIST1
!pb         ykap(ipo)   = ykap(IPO) + COPV(11,NCLTAL(ITET))*DIST1
!pb         zkap(ipo)   = zkap(IPO) + COPV(12,NCLTAL(ITET))*DIST1
!   kappa_rad from laplace Te
            xkap(ipo)   = xkap(IPO) + COPV(10,NCLTAL(ITET))*DIST1
          END DO                 ! J
        END DO                   ! ITET

        where (abs(volsum(1:ncoor)) > 1.D-10) 
          EMIS(1:ncoor) = EMIS(1:ncoor)/volsum(1:ncoor)
          ABSORB(1:ncoor) = ABSORB(1:ncoor)/volsum(1:ncoor)
          EMSAMP(1:ncoor) = EMSAMP(1:ncoor)/volsum(1:ncoor)
          ABSAMP(1:ncoor) = ABSAMP(1:ncoor)/volsum(1:ncoor)
          EM1SMP(1:ncoor) = EM1SMP(1:ncoor)/volsum(1:ncoor)
          EM2SMP(1:ncoor) = EM2SMP(1:ncoor)/volsum(1:ncoor)
          AB1SMP(1:ncoor) = AB1SMP(1:ncoor)/volsum(1:ncoor)
          AB2SMP(1:ncoor) = AB2SMP(1:ncoor)/volsum(1:ncoor)
          RIX(1:ncoor) = RIX(1:ncoor)/volsum(1:ncoor)*elcha
          RIY(1:ncoor) = RIY(1:ncoor)/volsum(1:ncoor)*elcha
          RIZ(1:ncoor) = RIZ(1:ncoor)/volsum(1:ncoor)*elcha
!pb       XKAP(1:ncoor) = XKAP(1:ncoor)/volsum(1:ncoor)
!pb       YKAP(1:ncoor) = YKAP(1:ncoor)/volsum(1:ncoor)
!pb       ZKAP(1:ncoor) = ZKAP(1:ncoor)/volsum(1:ncoor)
          XKAP(1:ncoor) = XKAP(1:ncoor)/volsum(1:ncoor)
        end where

! TALLIES FOR FIDAP INTERFACE: DONE
!  write output file for FIDAP, stream 35


        IF (ISTRAI == 1) THEN
          WRITE (35,*) NSTRAI
          WRITE (35,*) NCOORD
        END IF
        WRITE (35,*) TXTSOU(ISTRAI)
        WRITE (35,'(A9,1X,10A20)') 'NO. POINT','EMIS_1','ABSORPTION_1',
     .                            'EMIS_2','ABSORBTION_2',
     .                            'Kappa','I_x','I_y','I_z'
        DO IR=1,NCOORD
          WRITE (35,'(I9,1X,10ES20.10)') IR, EM1SMP(IR),AB1SMP(IR),
     .                                  EM2SMP(IR),AB2SMP(IR),
     .                                  xkap(ir),
     .                                  RIX(IR),RIY(IR),RIZ(IR)
        END DO
        WRITE (35,*)
C
C  DATA TRANSFER BACK TO PLASMA CODE FINISHED FOR STRATUM NO. ISTRAI
C
10000 CONTINUE
C
      RETURN
C
      ENTRY IF4COP
C
      NREC11=NOUTAU
      OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
      IRC=3
      WRITE (11,REC=IRC) RCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ALLOCATE (IHELP(NOUTAU))
      IHELP=0
      JC=0
      DO K=1,NPTRGT
        DO J=1,10*NSTEP
          JC=JC+1
          IHELP(JC)=ICCPL1(J,K)
          IF (JC == NOUTAU) THEN
            IRC=IRC+1
            WRITE (11,REC=IRC) IHELP
            IF (TRCINT.OR.TRCFLE) 
     .        WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
            JC=0
          END IF
        END DO
      END DO
      IF (JC > 0) THEN
        IRC=IRC+1
        WRITE (11,REC=IRC) IHELP
        IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      END IF
      DEALLOCATE (IHELP)
      IRC=IRC+1
      WRITE (11,REC=IRC) ICCPL2
      IRC=IRC+1
      WRITE (11,REC=IRC) LCCPL
      IF (TRCINT.OR.TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
C
      IF (LSHORT) LSTOP=LSTP
C
      IF (.NOT.LSTOP) RETURN
C
      IF (.NOT.(LBALAN)) GOTO 11000
C
C  BALANCES, SHOULD BE DONE ONLY AT THE END OF BRAAMS RUN
C  AT THE END OF AN EIRENE RUN THE BALANCES MAY BE OFF AT LEAST AT
C  THE BEGINNING OF THE CYCLING PROCEDURE, BECAUSE THE PLASMA STILL
C  HAS TO ADJUST TO THE NEW SOURCES
C
C
C
11000 CONTINUE
C
      DEALLOCATE (NODBRICK)
      RETURN
C
      END


      subroutine kreuzprod (a,b,vnorm,c)
      use precision
      implicit none
      real(dp), intent(in) :: a(3), b(3), vnorm
      real(dp), intent(out) :: c(3)

      c(1) = (a(2)*b(3) - a(3)*b(2)) / vnorm
      c(2) = (a(3)*b(1) - a(1)*b(3)) / vnorm
      c(3) = (a(1)*b(2) - a(2)*b(1)) / vnorm

      return
      end subroutine kreuzprod
C ===== SOURCE: mshproj.f
C
C
      SUBROUTINE MSHPROJ(X1,Y1,X2,Y2,X3,Y3,X4,Y4,PUX,PUY,PVX,PVY,
     .                   NDXA,NR1ST,IY)

      USE PRECISION

      IMPLICIT NONE
      REAL(DP), INTENT(IN) :: X1(*),Y1(*),X2(*),Y2(*),X3(*),Y3(*),
     .                      X4(*),Y4(*)
      INTEGER, INTENT(IN) :: NDXA,NR1ST,IY
      REAL(DP) :: PUX(*),PUY(*),PVX(*),PVY(*)
      REAL(DP) :: EPS60, D12, D34, D13, D24, DUX, DUY, DVX, DVY,
     .          PUPV, PVPV
      INTEGER :: IX, IN
      EPS60 = 1.E-60_DP
C
C
      DO 1 IX=1,NDXA
C
C  CALCULATE THE NORM OF THE VECTORS (POINT2-POINT1),....
C
        D12 = SQRT((X2(IX)-X1(IX))*(X2(IX)-X1(IX))+(Y2(IX)-Y1(IX))*
     .        (Y2(IX)-Y1(IX)))+EPS60
        D34 = SQRT((X4(IX)-X3(IX))*(X4(IX)-X3(IX))+(Y4(IX)-Y3(IX))*
     .        (Y4(IX)-Y3(IX)))+EPS60
        D13 = SQRT((X3(IX)-X1(IX))*(X3(IX)-X1(IX))+(Y3(IX)-Y1(IX))*
     .        (Y3(IX)-Y1(IX)))+EPS60
        D24 = SQRT((X4(IX)-X2(IX))*(X4(IX)-X2(IX))+(Y4(IX)-Y2(IX))*
     .        (Y4(IX)-Y2(IX)))+EPS60
C
C  CALCULATE THE BISSECTING VECTORS, BUT NOT NORMALISED YET
C
        DUX = (X2(IX)-X1(IX))/D12 + (X4(IX)-X3(IX))/D34
        DUY = (Y2(IX)-Y1(IX))/D12 + (Y4(IX)-Y3(IX))/D34
        DVX = (X3(IX)-X1(IX))/D13 + (X4(IX)-X2(IX))/D24
        DVY = (Y3(IX)-Y1(IX))/D13 + (Y4(IX)-Y2(IX))/D24
C
C  CALCULATE THE COMPONENTS OF THE TWO UNIT VECTOR (= PROJECTION RATE)
C
        IN=IY+(IX-1)*NR1ST
        PUX(IN) = DUX/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PUY(IN) = DUY/(SQRT(DUX*DUX+DUY*DUY)+EPS60)
        PVX(IN) = DVX/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
        PVY(IN) = DVY/(SQRT(DVX*DVX+DVY*DVY)+EPS60)
C
C  ORTHOGONORMALIZE, CONSERVE ORIENTATION (E.SCHMIDT)
C
        PUPV=PUX(IN)*PVX(IN)+PUY(IN)*PVY(IN)
        PVX(IN)=PVX(IN)-PUPV*PUX(IN)
        PVY(IN)=PVY(IN)-PUPV*PUY(IN)
        PVPV=SQRT(PVX(IN)*PVX(IN)+PVY(IN)*PVY(IN))+EPS60
        PVX(IN)=PVX(IN)/PVPV
        PVY(IN)=PVY(IN)/PVPV
C
1     CONTINUE
      RETURN
      END
C ===== SOURCE: neutr.f
C
C
*//NEUTR//
C=======================================================================
C          S U B R O U T I N E   N E U T R
C=======================================================================
      SUBROUTINE NEUTR(KARD,NDIMX,NDIMY,NDIMF,DUMMY,LDMX,LDMY,LDMF,
     .                 LDNS,IS)

      USE PRECISION

      IMPLICIT NONE
      INTEGER, INTENT(IN):: KARD,NDIMX,NDIMY,NDIMF,LDMX,LDMY,LDMF,
     .                      LDNS,IS
      INTEGER :: ND1,LIM,IX,IY,III,IF
      REAL(DP), INTENT(IN) :: DUMMY(0:LDMX+1,0:LDMY+1,LDMF,LDNS)
C
      ND1 = NDIMX
      LIM = (ND1/5)*5 - 4
      DO  500  IF = 1,NDIMF
        DO  110  IY = 1,NDIMY
          DO  100  IX = 1,LIM,5
  100     WRITE(KARD,910) (DUMMY(IX-1+III,IY,IF,IS),III = 1,5)
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX,IY,IF,IS),IX = LIM+5,ND1)
  110   CONTINUE
  500 CONTINUE
      RETURN
  910 FORMAT(5(E16.8))
*//END NEUTR//
      END
C ===== SOURCE: plasm.f


*//PLASM//
C=======================================================================
C          S U B R O U T I N E   P L A S M
C=======================================================================
      SUBROUTINE PLASM(KARD,NDIMX,NDIMY,NDIMF,N,M,NF,DUMMY)

      USE PRECISION

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: KARD,NDIMX,NDIMY,NDIMF,N,M,NF
      REAL(DP), INTENT(OUT) :: DUMMY(0:N+1,0:M+1,NF)
      INTEGER :: ND1, LIM, IX, IY, IF, III
      ND1 = NDIMX + 2
      LIM = (ND1/5)*5 - 4
      DO    110  IF = 1,NDIMF
      DO    110  IY = 0,NDIMY+1
      DO    100  IX = 1,LIM,5
100     READ(KARD,910) (DUMMY(-1+IX-1+III,IY,IF),III = 1,5)
        IF( (LIM+4).EQ.ND1 )     GOTO 110
        READ(KARD,910) (DUMMY(-1+IX,IY,IF),IX = LIM+5,ND1)
110   CONTINUE
      RETURN
910   FORMAT(5(E16.8))
*//END PLASM//
      END
C ===== SOURCE: statis_cop.f
C  STANDARD DEVIATIONS FOR COPV TALLIES, AND FOR SOME SPECIAL
C                      ALGEBRAIC FUNCTIONS (SUMS, ....) AMONGST THEM.
C
C  THIS VERSION:  COUPLING TO FIDAP FOR RADIATION TRANSFER
C                 COPV(1:6), SEE SUBR. UPTCOP

C                 ONE EXTRA COPV-TALLY: COPV(NCPVI+1), WHICH IS AN
C                 ALGEBRAIC TALLY (SUM): COPV(3)+COPV(4)+COPV(5)+COPV(6)
C
      SUBROUTINE STATIS_COP
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CSDVI
      USE CSDVI_COP
      USE COUTAU

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP, LT
      REAL(DP) :: XNM, SD2, DS, ZFLUXQ, SD2S, SDS, SDI, SDE,
     .          D2S, SG, DSA, DD, D, SG2, DA, FSIG, XN, ZFLUX,
     .          SD1, SD1S , SAV
      INTEGER :: IPLS, NSB, IR, ICO, IPL, NR1, NP2, NT3, ICPV,
     .           I,J,IIN,NRW,IRU,ITS,ITL,ISCO,IGI,IGE,NSYM,NSYH,
     .           IGS,IG,IT,IP,J1,J2
      REAL(DP), ALLOCATABLE, SAVE :: VECTOR(:), 
     .          SD(:)   
      INTEGER, ALLOCATABLE, SAVE ::                               
     .                              IND(:,:),   IIND(:),    INDSS(:,:)
C
      SAVE
C
      ENTRY STATS0_COP

      IMETCL = 0
      NCLMT = 0
      LMETSP = .FALSE.

      IF (NSIGI_COP.EQ.0) RETURN
C
      IF (.NOT.ALLOCATED(IND)) THEN
        AllOCATE (IND(NRTAL,8))  
        AllOCATE (IIND(NRTAL)) 
        AllOCATE (INDSS(NRTAL,8))
        AllOCATE (VECTOR(MAX(NRTAL,NLMPGS)))
        AllOCATE (SD(0:(MAX(NRTAL,NLMPGS))))
        SD=0._DP
      END IF

C  FILL IIND, INDSS ARRAYS FOR THOSE "AVERAGE" CELLS, TO WHICH "REAL" 
C  CELL IR ALSO CONTRIBUTES
C  IIND: HOW MANY CELLS
C  INDSS: WHICH CELLS
      CALL INDTAL(IND,NRTAL,NR1TAL,NP2TAL,NT3TAL,NBMLT)
      DO IR=1,NSBOX_TAL
        IIND(IR)=0
        IIN=0
        DO J=1,8
          IF (IND(IR,J).NE.0) THEN
            IIND(IR)=IIND(IR)+1
            IIN=IIN+1
            INDSS(IR,IIN)=J
          ENDIF
        ENDDO
      ENDDO
C
      RETURN

C
      ENTRY STATS1_COP(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
C
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
      NRW=NSIN

      IF (NSIGI_COP == 0) RETURN

C  HISTORY HAS TOUCHED NCLMT CELLS.
C  IT CONTRIBUTES TO NCLMTS CELLS (AVERAGES)
C  THESE CELL NUMBERS ARE STORED HERE ON ICLMT ARRAY
        NCLMTS = MAX(NCLMT,NCLMTS)
        DO I=1,NCLMT
          IR = ICLMT(I)
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            IF (IMETCL(IRU) == 0) THEN
              NCLMTS = NCLMTS+1
              IMETCL(IRU) = NCLMTS
              ICLMT(NCLMTS) = IRU
            END IF
          END DO
        END DO
C
C  STANDARD DEVIATION FOR ALL TALLIES COPV(1:NCPVI)
C
      DO 1012 ICPV=1,NCPVI
C  ARE THERE CONTRIBUTIONS TO THE REQUESTED TALLY FROM THIS HISTORY?
        IGS=ICPV
        ITL=NTALM
        ISCO = 0
        IF (LMETSP(NSPAN(ITL)+IGS-1)) ISCO = 1
        IF (ISCO == 0) GOTO 1012
C
        IF (.NOT.LP.AND..NOT.LT) GOTO 1005
C  USE SYMMETRY IN POLOIDAL/Y AND/OR TOROIDAL/Z COORDINATE
        IGI=IGS
        IGE=IGS
        IF (LP) THEN
          NSYM=NP2
          NSYH=(NSYM-1)/2
          DO 1003 IG=IGI,IGE
          DO 1003 IR=1,NR1
          DO 1003 IT=1,NT3
          DO 1003 IP=1,NSYH
                J1=IR+((IT-1)*NP2+IP-1)*NR1
                J2=IR+((IT-1)*NP2+NSYM-IP-1)*NR1
                SAV=(COPV(ICPV,J1)+COPV(ICPV,J2))*0.5
                COPV(ICPV,J1)=SAV
                COPV(ICPV,J2)=SAV
1003      CONTINUE
        ENDIF
        IF (LT) THEN
          NSYM=NT3
          NSYH=(NSYM-1)/2
          DO 1004 IG=IGI,IGE
          DO 1004 IR=1,NR1
          DO 1004 IP=1,NP2
          DO 1004 IT=1,NSYH
                J1=IR+((IT-1)*NP2+IP-1)*NR1
                J2=IR+((NSYM-IT-1)*NP2+IP-1)*NR1
                SAV=(COPV(ICPV,J1)+COPV(ICPV,J2))*0.5
                COPV(ICPV,J1)=SAV
                COPV(ICPV,J2)=SAV
1004      CONTINUE
        ENDIF
1005    CONTINUE
C
C  FILL ARRAY VECTOR, FOR INDIVIDUAL COPV-TALLIES
C  VECTOR IS FILLED ONLY FOR THOSE CELLS, 
C  WHICH HAVE BEEN TOUCHED BY THIS HISTROY
C  VECTOR CONTAINS THE SUM FROM ALL HISTORIES IN THESE CELLS UP TO THE
C  PRESENT HISTORY
       DO ICO = 1,NCLMT
         IR = ICLMT(ICO)
         VECTOR(ICO)=COPV(ICPV,IR)
       END DO

        SD1S = 0.D0
C  FILL ARRAY SD WITH THE INDIVIDUAL CONTRIBUTION FROM THIS HISTROY,
C  IN EACH CELL THAT HAS BEEN TOUCHED BY THIS HISTORY
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTOR(ICO)-SDVIA_COP(ICPV,IR)
          SD1S=SD1S+SD1
          SDVIA_COP(ICPV,IR)=VECTOR(ICO)
          SD(IR) = SD1
C  FILL SD ALSO FOR OTHER "CELL", TO WHICH CELL IR CONTRIBUTES
C  I.E., AVERAGES OVER COORDINATES OR OVER THE ENTIRE COMPUTATIONAL DOMAIN
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SD(IR)
          SIGMA_COP(ICPV,IR)=SIGMA_COP(ICPV,IR)+SD1*SD1
          SD(IR)=0._DP
        END DO

        SGMS_COP(ICPV)=SGMS_COP(ICPV)+SD1S*SD1S
1012  CONTINUE
C

C  STATISTICS FOR ALGEBRAIC FUNCTION OF TALLIES, AS NEEDED FOR COUPLING
C  THIS PART IS SPECIFIC FOR THIS VERSION OF COUPLE_....
C  HERE:
C  STATISTICS FOR TOTAL SAMPLED EMISSION- TOTAL ABSORPTION
C   
      ICPV=NCPVI+1
      IF (LMETSP(NSPAN(NTALM)+2) .OR.
     .    LMETSP(NSPAN(NTALM)+3) .OR.
     .    LMETSP(NSPAN(NTALM)+4) .OR.
     .    LMETSP(NSPAN(NTALM)+5) ) THEN

C  FILL ARRAY VECTOR, FOR INDIVIDUAL COPV-TALLIES
C  VECTOR IS FILLED ONLY FOR THOSE CELLS, 
C  WHICH HAVE BEEN TOUCHED BY THIS HISTROY
C  VECTOR CONTAINS THE SUM FROM ALL HISTORIES IN THESE CELLS UP TO THE
C  PRESENT HISTORY
       DO ICO = 1,NCLMT
         IR = ICLMT(ICO)
         VECTOR(ICO)=COPV(3,IR)+COPV(4,IR)+COPV(5,IR)+COPV(6,IR)
       END DO

        SD1S = 0.D0
C  FILL ARRAY SD WITH THE INDIVIDUAL CONTRIBUTION FROM THIS HISTROY,
C  IN EACH CELL THAT HAS BEEN TOUCHED BY THIS HISTORY
        DO ICO = 1,NCLMT
          IR = ICLMT(ICO)
          SD1 = VECTOR(ICO)-SDVIA_COP(ICPV,IR)
          SD1S=SD1S+SD1
          SDVIA_COP(ICPV,IR)=VECTOR(ICO)
          SD(IR) = SD1
C  FILL SD ALSO FOR OTHER "CELL", TO WHICH CELL IR CONTRIBUTES
C  I.E., AVERAGES OVER COORDINATES OR OVER THE ENTIRE COMPUTATIONAL DOMAIN
          DO IIN=2,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO ICO = 1,NCLMTS
          IR = ICLMT(ICO)
          SD1=SD(IR)
          SIGMA_COP(ICPV,IR)=SIGMA_COP(ICPV,IR)+SD1*SD1
          SD(IR)=0._DP
        END DO

        SGMS_COP(ICPV)=SGMS_COP(ICPV)+SD1S*SD1S
      END IF

C
1020  CONTINUE
      RETURN
C
      ENTRY STATS2_COP(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C           (=1, FALLS NUR EIN BEITRAG UNGLEICH 0, ODER (KUENSTLICH
C            ERZWUNGEN) FALLS GAR KEIN BEITRAG UNGLEICH NULL)
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.D0) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NCPVI.EQ.0) GOTO 2200
C
      DO 2112 ICPV=1,NCPVI
C
        DS=0.
        DO 2011 IR=1,NSB
          SD1=COPV(ICPV,IR)
          DS=DS+SD1
          DO 2016 IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
2016      CONTINUE
2011    CONTINUE

        DO 2111 IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_COP(ICPV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_COP(ICPV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_COP(ICPV,IR)=STV_COP(ICPV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_COP(ICPV,IR)=EE_COP(ICPV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
2111    CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_COP(ICPV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_COP(ICPV)=SG*FSIG
C
        STVS_COP(ICPV)=STVS_COP(ICPV)+SG2*ZFLUXQ/XNM/XN
        EES_COP(ICPV)=EES_COP(ICPV)+DS*ZFLUX/XN
2112  CONTINUE
C
C  STATISTICS FOR ALGEBRAIC FUNCTION OF TALLIES, AS NEEDED FOR COUPLING
C  HERE:
C  STATISTICS FOR TOTAL SAMPLED EMISSION- TOTAL ABSORPTION
C
      ICPV = NCPVI+1
      DS=0.
      DO IR=1,NSB
        SD1=SUM(COPV(3:6,IR))
        DS=DS+SD1
        DO IIN=1,IIND(IR)
          J=INDSS(IR,IIN)
          IRU=IND(IR,J)
          SD(IRU)=SD(IRU)+SD1
        ENDDO
      ENDDO
      DO IR=1,NSB
        D=SD(IR)
        DD=D*D
        DA=ABS(D)
        SG2=MAX(0._DP,SIGMA_COP(ICPV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
        SG=SQRT(SG2)/(DA+EPS60)
        SIGMA_COP(ICPV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
        STV_COP(ICPV,IR)=STV_COP(ICPV,IR)+SG2*ZFLUXQ/XNM/XN
        EE_COP(ICPV,IR)=EE_COP(ICPV,IR)+D*ZFLUX/XN
      END DO
      D2S=DS*DS
      DSA=ABS(DS)
      SG2=MAX(0._DP,SGMS_COP(ICPV)-D2S/XN)
      SG=SQRT(SG2)/(DSA+EPS60)
      SGMS_COP(ICPV)=SG*FSIG
C
      STVS_COP(ICPV)=STVS_COP(ICPV)+SG2*ZFLUXQ/XNM/XN
      EES_COP(ICPV)=EES_COP(ICPV)+DS*ZFLUX/XN
C
2200  CONTINUE
      RETURN
      END
C ===== SOURCE: triquaint.f
      subroutine triquaint (vector, r, s, t)

      USE PRECISION

      implicit none
      real(dp), intent(out) :: vector(27)
      real(dp), intent(in) :: r, s, t
      real(dp) :: halb, viertel, achtel, emr, ems, emt, epr, eps, ept,
     .            emrq, emsq, emtq

      halb = 0.5_dp
      viertel = 0.25_dp
      achtel = 0.125_dp

      emr = 1._dp - r
      ems = 1._dp - s
      emt = 1._dp - t

      epr = 1._dp + r
      eps = 1._dp + s
      ept = 1._dp + t

      emrq = 1._dp - r*r
      emsq = 1._dp - s*s
      emtq = 1._dp - t*t

      vector(1) = -achtel * r * s * t * emr * ems * emt
      vector(2) = viertel * s * t * emrq * ems * emt
      vector(3) = achtel * r * s * t * epr * ems * emt
      vector(4) = viertel * r * t * emr * emsq * emt
      vector(5) = -halb * t * emrq * emsq * emt
      vector(6) = -viertel * r * t * epr * emsq * emt
      vector(7) = achtel * r * s * t * emr * eps * emt
      vector(8) = -viertel * s * t * emrq * eps * emt
      vector(9) = -achtel * r * s * t * epr * eps * emt

      vector(10) = viertel * r * s * emr * ems * emtq
      vector(11) = -halb * s * emrq * ems * emtq
      vector(12) = -viertel * r * s * epr * ems * emtq
      vector(13) = -halb * r * emr * emsq * emtq
      vector(14) = emrq * emsq * emtq
      vector(15) = halb * r * epr * emsq * emtq
      vector(16) = -viertel * r * s * emr * eps * emtq
      vector(17) = halb * s * emrq * eps * emtq
      vector(18) = viertel * r * s * epr * eps * emtq

      vector(19) = achtel * r * s * t * emr * ems * ept
      vector(20) = -viertel * s * t * emrq * ems * ept
      vector(21) = -achtel * r * s * t * epr * ems * ept
      vector(22) = -viertel * r * t * emr * emsq * ept
      vector(23) = halb * t * emrq * emsq * ept
      vector(24) = viertel * r * t * epr * emsq * ept
      vector(25) = -achtel * r * s * t * emr * eps * ept
      vector(26) = viertel * s * t * emrq * eps * ept
      vector(27) = achtel * r * s * t * epr * eps * ept

      return
      end subroutine triquaint
C ===== SOURCE: upscop.f
C
C
      SUBROUTINE UPSCOP
      RETURN
      END
C ===== SOURCE: uptcop.f
C
      SUBROUTINE UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
C
C  UPDATE TALLIES COPV FOR COUPLING TO OTHER CODES
C
C  THIS VERSION: COUPLING TO FIDAP, FOR RADIATION TRANSFER
C
C  COPV(1): TOTAL SAMPLED EMISSION, SAME AS EPPHT IN LOCATE.F
C                                   COLLISION ESTIMATOR
C  COPV(2): TOTAL ABSORPTION, SAME AS EPHPHT IN UPDATE.F
C                             COLLISION ESTIMATOR IN CELL OF EMISSION
C                             TRACKLENGTH ESTIMATOR ELSE
c  copv(1), copv(2) are redundant, in principle, but they (and their
c                   variances), are used in infcop, entry if3cop.
C
C  COPV(3): SAMPLED EMISSION, THICK PART, ESTIMATOR AS COPV(1)
C  COPV(4): ABSORPTION OF THICK PART, ESTIMATORS AS WITH COPV(2)
C  COPV(5): SAMPLED EMISSION, THIN PART, ESTIMATOR AS COPV(2)
C  COPV(6): ABSORPTION OF THIN PART, ESTIMATORS AS WITH COPV(2)
C
C  COPV(7), COPV(8), COPV(9):  VEC-I
C           NOW: SUM COPV(3)+....COPV(6), FOR VARIANCE (STATIS_COP) ON ADDV(5)
C                IN STATIS_COP: THIS SUM IS ON SIGMA_COP(NCPVI+1)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE CGEOM
      USE COMPRT
      USE CSDVI
      USE COMXS
      USE PHOTON

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                        XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: P, WTRSIG, EION, V0_PARB, PARMOM_0, DIST, WTR, XLCRIT
      INTEGER :: IFIRST, ISP, ICOU, IRDO, IRD, IAOT, IROT, UPDF
      INTEGER, SAVE :: NMTSP
      LOGICAL, SAVE :: LTHICK
CDR
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
C
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI
C
      ENDIF

C
C  WV=WEIGHT/VEL
C
C  PHOTONS
      IF (ITYP.EQ.0) THEN

!  EMISSION

        IF (NSTCLL > 0) THEN
          IRD=NCLTAL(NSTCLL)
!  TOTAL SAMPLED EMISSION
          COPV(1,IRD) = COPV(1,IRD) + STEMIS
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+1)=.TRUE.
C
C  COPV(3): THICK FRACTION OF SAMPLED EMISSION
C  COPV(5): THIN  FRACTION OF SAMPLED EMISSION
C

!  ZMFP < ZMFPTHI*XLCRIT  => THICK FRACTION
!  XLCRIT = TEMPERATURE GRADIENT LENGTH (CM) IN FLIGHT DIRECTION

          XLCRIT = TEDTEDX(IRD)*VELX+TEDTEDY(IRD)*VELY+TEDTEDZ(IRD)*VELZ

          IF (ABS(XLCRIT)*ZMFPI > 1._dp/ZMFPTHI) THEN
! XLCRIT/zmfp > ZMFPHT  => THICK FRACTION 
!  EMISSION PART 1
            COPV(3,IRD) = COPV(3,IRD) + STEMIS  
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+3)=.TRUE.
            ISP = 4
            LTHICK = .TRUE.
          ELSE
!  XLCRIT/zmfp < ZMFPHT  => THIN FRACTION 
!  EMISSION PART 2
            COPV(5,IRD) = COPV(5,IRD) + STEMIS  
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+5)=.TRUE.
            ISP = 6
            LTHICK = .FALSE.
          END IF
          IF (IMETCL(IRD) == 0) THEN
            NCLMT = NCLMT+1
            ICLMT(NCLMT) = IRD
            IMETCL(IRD) = NCLMT
          END IF
          NSTCLL = -1
        END IF

!  EMISSION DONE, NOW: ABSORPTION

C  COPV(2): TOTAL ABSORPTION
C  COPV(4): THICK FRACTION OF ABSORPTION
C  COPV(6): THIN  FRACTION OF ABSORPTION
C
        DO ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)

!  xlcrit = temperature gradient length (cm) in flight direction
          XLCRIT = TEDTEDX(IRD)*VELX+TEDTEDY(IRD)*VELY+TEDTEDZ(IRD)*VELZ
          IF (ABS(XLCRIT)*ZMFPI > 1._dp/ZMFPTHI) THEN
            ISP = 4
            LTHICK = .TRUE.
          ELSE
            ISP = 6
            LTHICK = .FALSE.
          END IF

          IF (IMETCL(IRD) == 0) THEN
            NCLMT = NCLMT+1
            ICLMT(NCLMT) = IRD
            IMETCL(IRD) = NCLMT
          END IF
C
          IF (LGVAC(IRDO,0)) CYCLE
C
          if (ncou.gt.1) then
          XSTOR(1:mstor1,1:mstor2) = XSTOR2(1:mstor1,1:mstor2,ICOU)
          XSTORV(1:nstorv) = XSTORV2(1:nstorv,ICOU)
          endif

          WTRSIG=WTR*(SIGTOT-SIGBGK)

C  COLLISION ESTIMATOR FOR COPV(2), COPV(4), COPV(6)

          IF ((LAST_EVENT%IFLAG == 1) .AND.
     .        (LAST_EVENT%NCELL == IRD)) THEN

!  use collision estimator for first cell ("brick") along track
!  in case of a collision sample 1 (the full weight)
!  in case of no collision sample 0 
            IF ((IFLAG == 4).OR.(IFLAG == 5)) THEN
!  TOTAL
              COPV(2,IRD) = COPV(2,IRD) - WEIGHT*E0
              COPV(12,IRD) = COPV(12,IRD) - WEIGHT*E0
!  FRACTION ACCORDING TO ISP
              COPV(ISP,IRD) = COPV(ISP,IRD) - WEIGHT*E0
!pb              if (nltrc) write(iunout,*) ' start cell '
!pb              if (nltrc) write(iunout,*) ' ird ',ird,LAST_EVENT%NCELL
!pb              if (nltrc) write(iunout,*) ' weight*e0 ',weight*e0
            ELSE
!             add 0 ==> nothing to be done
!pb              if (nltrc) write(iunout,*) ' 0 added, ird ', ird
            END IF

          ELSE

C TRACKLENGTH ESTIMATOR FOR COPV(2), COPV(4), COPV(6)

!  TOTAL, TRACKLENGTH
            COPV(2,IRD) = COPV(2,IRD) - WTRSIG*E0
            COPV(13,IRD) = COPV(13,IRD) - WTRSIG*E0
!  FRACTION ACCORDING TO ISP, TRACKLENGTH
            COPV(ISP,IRD) = COPV(ISP,IRD) - WTRSIG*E0
C
!pb            if (nltrc) write (iunout,*) ' normal case '
!pb            if (nltrc) write (iunout,*) ' ird, wtrsig*e0 ',ird,WTRSIG*E0

          END IF
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+2)=.TRUE.
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+ISP)=.TRUE.

!  I_x, I_y, I_z , radiative  vector flux, to be related to  kappa* grad(Te)
          IF (LTHICK) THEN
            COPV(7,IRD) = COPV(7,IRD) + WTR*E0*VEL*VELX
            COPV(8,IRD) = COPV(8,IRD) + WTR*E0*VEL*VELY
            COPV(9,IRD) = COPV(9,IRD) + WTR*E0*VEL*VELZ
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+7)=.TRUE.
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+8)=.TRUE.
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+9)=.TRUE.
          END IF
C
C
        END DO
C
      ENDIF
C
      RETURN
      END
