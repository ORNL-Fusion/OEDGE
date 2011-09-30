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
