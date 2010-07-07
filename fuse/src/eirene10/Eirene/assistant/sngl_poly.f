      function sngl_poly (cf, al, rcmin, rcmax, fpp, ifexmn, ifexmx) 
     .                   result(cou)

      use precision
      
      implicit none

      real(dp), intent(in) :: cf(9), fpp(6)
      real(dp), intent(in) :: al, rcmin, rcmax
      integer, intent(in) :: ifexmn, ifexmx
      real(dp) :: cou, fp(6), s01, s02, ds12, expo1, expo2, ccxm1, 
     .            ccxm2, extrap 
      integer :: ii, if8, ifex

      if ((ifexmn .ne. 0) .and. (al < rcmin)) then

C  ELAB BELOW MINIMUM ENERGY FOR FIT:

        FP = FPP
        ifex = ifexmn

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN
        IF (IFEXMN.LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(SIGMA)
          S01=RCMIN
          S02=LOG(2._DP)+RCMIN
          DS12=S02-S01
          EXPO1=CF(9)
          EXPO2=CF(9)
          DO 1 II=1,8
            IF8=9-II
            EXPO1=EXPO1*S01+CF(IF8)
            EXPO2=EXPO2*S02+CF(IF8)
 1        CONTINUE
          CCXM1=EXPO1
          CCXM2=EXPO2
          FP(1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
          FP(2)=      (CCXM2-CCXM1)/DS12
          FP(3)=0.D0
C
          IFEX=5
        ENDIF

        COU=EXTRAP(AL,IFEX,FP(1),FP(2),FP(3))
        cou = log(cou)

      elseif ((ifexmx .ne. 0) .and. (al > rcmax)) then

C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:

C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
        FP = FPP
        COU=EXTRAP(AL,IFEXMX,FP(4),FP(5),FP(6))
        cou = log(cou)

      else

C  ELAB WITHIN ENERGY RANGE OF FIT:

        cou = cf(9)

        do ii = 8, 1, -1
          cou = cou * al + cf(ii)
        end do

      end if

      return
      end function sngl_poly
