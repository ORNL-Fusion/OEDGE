!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  24.11.06: get extrapolation parameters for polynomial fit only
!pb  07.12.06: double declaration of dsub removed
 
      function EIRENE_rate_coeff (ir, p1, p2, lexp, iprshft, erate)
     .                     result (rate)
 
!   lexp:     rate coefficient in cm**3/sec
!  .not.lexp: ln(rate coefficient in cm**3/sec)
 
      use EIRMOD_precision
      use EIRMOD_parmmod
      use EIRMOD_comxs
      use EIRMOD_comprt, only: iunout
 
      implicit none
 
      integer, intent(in) :: ir, iprshft
      real(dp), intent(in) :: p1, p2
      logical, intent(in) :: lexp
      real(dp), intent(out) :: erate
      real(dp) :: rate, eirene_sngl_poly, dum(9), rcmin, rcmax,
     .            fp(6), q1, q2,
     .            ALPCR, SCR, SCRRAD, E_ALPCR, E_SCR, E_SCRRAD,
     .            E_ALPCR_T, E_SCR_T, E_SCRRAD_T
      real(dp), save :: xlog10e, xln10, dsub
      real(dp), allocatable, save :: pop0(:), pop1(:), pop2(:), qcol2(:)
      integer :: jfexmn, jfexmx
      integer, save :: ifirst=0, ifsub=0
 
      interface
        function EIRENE_intp_adas (ad,p1,p2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: adas_data
          type(adas_data), pointer :: ad
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function EIRENE_intp_adas
 
        function EIRENE_intp_table (tb,p1,p2) result(res)
          use EIRMOD_precision
          use EIRMOD_comxs, only: hydkin_data
          type(hydkin_data), pointer :: tb
          real(dp), intent(in) :: p1, p2
          real(dp) :: res
        end function EIRENE_intp_table
      end interface
 
 
      if (.not.reacdat(ir)%lrtc) then
        write (iunout,*) ' no data for rate coefficient available',
     .                    ' for reaction ',ir
        call EIRENE_exit_own(1)
      end if
 
      rate = 0._dp
      erate = 0._dp
 
      if ((reacdat(ir)%rtc%ifit == 1) .or.
     .    (reacdat(ir)%rtc%ifit == 2)) then
        rcmin  = reacdat(ir)%rtc%poly%rcmn
        rcmax  = reacdat(ir)%rtc%poly%rcmx
        fp     = reacdat(ir)%rtc%poly%fparm
        jfexmn = reacdat(ir)%rtc%poly%ifexmn
        jfexmx = reacdat(ir)%rtc%poly%ifexmx
      end if
 
      if (mod(iftflg(ir,2),100) == 10) then
 
        rate = reacdat(ir)%rtc%poly%dblpol(1,1)
 
      elseif (reacdat(ir)%rtc%ifit == 1) then
 
        rate = eirene_sngl_poly(reacdat(ir)%rtc%poly%dblpol(1:9,1),p1,
     .                   rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtc%ifit == 2) then
 
        if (ifsub == 0) then
          ifsub = 1
          dsub = log(1.e8_dp)
        end if
 
        q2 = p2
        if (iprshft > 0) q2 = q2 - dsub
 
        call EIRENE_dbl_poly
     .       (reacdat(ir)%rtc%poly%dblpol,p1,q2,rate,dum,1,9,
     .                rcmin, rcmax, fp, jfexmn, jfexmx)
        if (lexp) rate = exp(max(-100._dp,rate))
 
      else if (reacdat(ir)%rtc%ifit == 3) then
 
! ADAS
        if (ifirst == 0) then
          ifirst = 1
          xln10 = log(10._dp)
          xlog10e = 1._dp/xln10
        end if
 
        q1 = xlog10e*p1
        q2 = xlog10e*p2
        rate = eirene_intp_adas(reacdat(ir)%rtc%adas,q1,q2)
 
        if (lexp) then
          rate=10._dp**rate
        else
          rate = xln10*rate
        end if
 
      else if (reacdat(ir)%rtc%ifit == 4) then
 
        q1 = exp(p1)
        rate = eirene_intp_table(reacdat(ir)%rtc%hyd,q1,p2)
 
      else if (reacdat(ir)%rtc%ifit == 5) then
 
! H-colrad
 
        if (.not.allocated(pop0)) then
          allocate(pop0(40))
          allocate(pop1(40))
          allocate(pop2(40))
          allocate(qcol2(40))
        end if
 
        QCOL2 = 0._DP
        Q1 = EXP(P1)
        Q2 = EXP(P2)
        CALL EIRENE_H_COLRAD(Q1, Q2 ,QCOL2,POP0,POP1,POP2,
     .                ALPCR,    SCR,    SCRRAD,
     .                E_ALPCR,  E_SCR,  E_SCRRAD,
     .                E_ALPCR_T,E_SCR_T,E_SCRRAD_T)
 
        rate = log(scr)
        erate = log(-e_scr)
 
      end if
 
      return
 
      end function EIRENE_rate_coeff
 
