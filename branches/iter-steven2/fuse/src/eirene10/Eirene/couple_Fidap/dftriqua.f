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
