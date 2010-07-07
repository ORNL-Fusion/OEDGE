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
