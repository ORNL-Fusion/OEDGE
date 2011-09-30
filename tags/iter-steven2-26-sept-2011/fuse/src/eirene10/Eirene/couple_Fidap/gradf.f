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
