      subroutine EIRENE_dbl_poly (cf, al, pl, cou, dum, ji, je,
     .                     rcmin, rcmax, fp, ifexmn, ifexmx )
 
      use EIRMOD_precision
 
      implicit none
 
      real(dp), intent(in) :: cf(9,9), fp(6)
      real(dp), intent(in) :: al, pl, rcmin, rcmax
      real(dp), intent(out) :: cou, dum(9)
      integer, intent(in) :: ji, je, ifexmn, ifexmx
      real(dp) :: s01, s02, ds12, expo1, expo2, ccxm1, ccxm2,
     .            fpar1, fpar2, fpar3, eirene_extrap
      integer :: kk, jj, j, i, ii, ifex
 
      dum = 0._dp
 
      if ((ifexmn .ne. 0) .and. (al < rcmin)) then
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(<S*V>)
        S01=RCMIN
        S02=LOG(2._DP)+RCMIN
        DS12=S02-S01
        EXPO1=0.
        EXPO2=0.
        DO 1 J=1,9
          JJ=J-1
          DO 1 I=1,9
            II=I-1
            EXPO1=EXPO1+S01**II*PL**JJ*CF(I,J)
            EXPO2=EXPO2+S02**II*PL**JJ*CF(I,J)
 1      CONTINUE
        CCXM1=EXPO1
        CCXM2=EXPO2
        FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
        FPAR2=      (CCXM2-CCXM1)/DS12
        FPAR3=0.D0
C
        IFEX=5
        COU=EIRENE_EXTRAP(AL,IFEX,FPAR1,FPAR2,FPAR3)
        cou=log(cou)
 
      else
 
        do jj = je, ji, -1
          dum(jj) = cf(9,jj)
          do kk = 8, 1, -1
            dum(jj) = dum(jj) * al + cf(kk,jj)
          end do
        end do
 
        cou = dum(9)
 
        do jj = 8, 1, -1
          cou = cou * pl + dum(jj)
        end do
 
      end if
 
      return
      end subroutine EIRENE_dbl_poly
