module mod_sol22_sources



  implicit none
  !
  ! This module contains all the code related to setting up and returning the source values along the field line.
  ! - particles, power, momentum, pressure ... both actual and estimates
  !
  !




contains


  real*8 function press(s,tecur,ticur)
    use mod_solparams
    use mod_solcommon
    use mod_solswitch
    implicit none
    real*8 s,tecur,ticur
    !
    !     This function returns the value of the pressure at a position
    !     s along the field line, at the moment the only contribution
    !     other than the Pinf is the Momentum loss to neutrals term.
    !
    real*8 rpos,pmomloss_tmp

    !real*8 pmomloss,estpint,majrpos,estppress
    !external pmomloss,estpint,majrpos,estppress

    pmomloss_tmp = pmomloss(s,0,1.0d0,tecur,ticur)

    !
    !     jdemod - pmomloss returns 0.0 if the option is off so this 
    !              extra check code is not required. 
    !
    !      if (actswnmom.eq.0.0) then
    !         if (actswmajr.eq.4.0) then
    !            rpos = majrpos(s)
    !            press = (pinf + estpint(s)+ padd)/rpos
    !         else
    !            press = pinf + padd
    !         endif
    !      else
    !
    if (actswmajr.eq.4.0) then
       rpos = majrpos(s)
       press = (pinf + estpint(s) + pmomloss_tmp + padd  + estppress(s))/rpos
    else
       press = pinf + pmomloss_tmp + padd +estppress(s)
    endif
    !
    !      endif
    !
    return
  end function press



  real*8 function pintupdt(s,n,te,ti)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s,n,te,ti
    !
    !     PINTUPDT: This routine updates the accumulated integral
    !               of the major radius pressure correction term.
    !
    !     Note: Initialization done by call from INITVAL
    !
    !
    !
    common /pint/ lastpint,lasts,rlasts,plasts,nlasts
    real*8 lastpint,lasts,rlasts,plasts,nlasts

    real*8 rpos,pcur
    !real*8 majrpos
    !external majrpos
    !
    !     Only if major radius correction option is ON
    !
    if (actswmajr.ne.4.0) then
       pintupdt = 0.0
       return
    endif

    if (s.eq.0.0.or.s.lt.lasts) then
       lasts = 0.0
       lastpint = 0.0
       pintupdt = 0.0
       plasts = pstatic0
       rlasts = r0init
       return
    elseif (s.eq.lasts) then
       pintupdt = lastpint
       return
    endif

    rpos = majrpos(s)
    pcur = n * (te+ti) * econv

    pintupdt = lastpint + 0.5 * (pcur+plasts) * (rpos-rlasts)

    rlasts = rpos
    lasts = s
    lastpint = pintupdt
    plasts = pcur

    return
  end function pintupdt



  real*8 function estpint(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     ESTPINT: This routine returns an estimate of the major radius
    !              integrated pressure correction term at a value
    !              of S - since the last update.
    !

    common /pint/ lastpint,lasts,rlasts,plasts,nlasts
    real*8 lastpint,lasts,rlasts,plasts,nlasts

    real*8 rpos
    !real*8 majrpos
    !external majrpos
    !
    !     Only if major radius correction option is ON
    !
    if (actswmajr.ne.4.0) then
       estpint = 0.0
       return
    endif

    if (s.eq.0.0.or.s.lt.lasts) then
       estpint = 0.0
       return
    elseif (s.eq.lasts) then
       estpint = lastpint
       return
    endif

    rpos = majrpos(s)

    estpint = lastpint + plasts * (rpos-rlasts)

    return
  end function estpint



  real*8 function pmomloss(s,opt,vcur,tecur,ticur)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use mod_sol22pmom
    implicit none
    real*8 s,vcur,tecur,ticur
    integer opt
    !
    !     This function will provide the momentum loss
    !     integrated to a point s. The options that depend
    !     on data transferred from Nimbus have been partially
    !     implemented.
    !
    !     The OPT argument is used when one requires the
    !     UPDATE/ESTIMATOR method of calculating an integral
    !     of the momentum source term. At present it is only
    !     of use for momentum option 4. However the mechanism
    !     may be employed for other momentum options as required.
    !     It is implemented in the same fashion as the ancillary
    !     power integration terms.
    !
    !     A value of 0 is a request for an estimate and a value of
    !     1 is an indicator to update the integral.
    !
    !      real*8 :: lasts,lastsmom,lastsrc,lastv,lastte

    real*8 src,teav,vav

    !real*8 srci,rcxmult,estscx
    !external srci,rcxmult,estscx
    real*8 :: vtmp,te_base

    logical temp_opt

    integer top,bot,mid,in
    !
    !     Initialize and set to zero for S = 0
    !
    pmomloss = 0.0
    !
    !     Exit for S = 0
    !
    if (s.eq.0) then 
       lastv = vcur
       lasts = s
       lastte = tecur
       return
    endif

    if (actswnmom.eq.0.0.or.(actswnmom.eq.6.0.and.(.not.pinavail))) then
       pmomloss = 0.0
    elseif (actswnmom.eq.1.0.or.(actswnmom.eq.7.0.and.(.not.pinavail))) then
       if (s.le.actlenmom*ringlen) then
          pmomloss = s * smom0
       else
          pmomloss = actlenmom * ringlen * smom0
       endif
    elseif (actswnmom.eq.2.0.or.(actswnmom.eq.8.0.and.(.not.pinavail))) then
       if (s.le.actlenmom*ringlen) then
          pmomloss = smom0 * actlammom * ringlen *(1.0d0 - exp (-s/(actlammom*ringlen)))
       else
          pmomloss = smom0 * actlammom * ringlen *(1.0d0 - exp (-actlenmom/actlammom))
          !            (1.0d0 - exp (-(actlenmom*ringlen)/(lammom*ringlen)))
       endif
    elseif (actswnmom.eq.3.0) then
       pmomloss = smom0 * srci(s)
    elseif (actswnmom.eq.4.0) then
       !
       !        The momentum loss is calculated using the term
       !        Smom = -m * vb * rcxmom * Siz(s)
       !
       !
       !        Return an estimate of Smom in the next interval
       !
       temp_opt=.true.
       te_base = 10.0

       if (temp_opt) then 
          vtmp = -sqrt(2.0*te_base*econv/(mb*mconv))
       else
          ! set to absolute so that flow reversal doesn't cause a pressure drop  (?)
          vtmp = -abs(vcur)
       endif


       if (opt.eq.0) then

          if (s.eq.soffset.or.s.lt.lasts) then
             pmomloss = 0.0
             return
          elseif (s.eq.lasts) then
             pmomloss = lastsmom
             return
          endif

          src = srci(s)

          pmomloss = lastsmom  &
               - mb * mconv * (vtmp+lastv)/2.0  &
               * rcxmom * rcxmult(lastte) &
               * (src-lastsrc) * smom_mult 

          !
          !           This now needs its own return statement to avoid
          !           double multiplication by smom_mult
          !

          return

       elseif (opt.eq.1) then

          vav = (vtmp + lastv)/2.0
          teav = (tecur+lastte)/2.0

          if (s.eq.soffset.or.s.lt.lasts) then

             lasts = soffset
             lastte = tecur
             lastv = vtmp
             lastsmom = 0.0
             lastsrc = 0.0
             pmomloss = 0.0

             !              write (6,'(a,10g12.5)') 'Mom lasts:',s,lasts,lastsrc,lastsmom,lastv,vtmp

             return
          elseif (s.eq.lasts) then
             pmomloss = lastsmom
             return
          endif

          src = srci(s)

          pmomloss = lastsmom - mb * mconv * vav * rcxmom * rcxmult(teav) * (src-lastsrc) * smom_mult

          lastv = vtmp
          lastte = tecur
          lastsmom = pmomloss
          lasts = s
          lastsrc = src
          !
          !           This now nees its own return statement to avoid
          !           double multiplication by smom_mult
          !
          return

       endif

    elseif (actswnmom.eq.5.0) then
       pmomloss = 0.0
    elseif (pinavail.and.(actswnmom.eq.6.0.or.actswnmom.eq.7.or.actswnmom.eq.8.0)) then
       !
       !        Options based on PIN when pin data is
       !        available.
       !
       !        Search for right cell
       !
       call binsearch(s,in)
       !
       if (in.eq.1) then
          pmomloss = intmomsrc(in) * s / sptscopy(1)
       else
          pmomloss = ( intmomsrc(in-1) + ( (intmomsrc(in)-intmomsrc(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1))))
       endif

    elseif (pinavail.and.(actswnmom.eq.9.0.or.actswnmom.eq.10.0)) then
       !
       !     WF'96: estimates pmomloss from NIMBUS nH, EH, and CXSIG()
       !
       pmomloss = estscx(s,nlast,ticur)

    endif
    !
    !     Apply overall multiplier to MOST options - some options exit before
    !     reaching this point.
    !
    pmomloss = pmomloss * smom_mult


    return
  end function pmomloss



  real*8 function rcxmult(t)
    use mod_solparams
    use mod_solcommon
    implicit none
    real*8 t
    !
    !     This function calculates a Cx/IZ multiplier for use
    !     in the neutral momentum loss term that contributes
    !     to the total pressure.
    !
    integer firstx
    data firstx /0/

    real*8 coeffa,coeffb
    !
    !     Initialize coefficients on first iteration
    !
    if (firstx.eq.0) then
       coeffb = 6.907/(tcxmom-1.0)
       coeffa = 1000.0 * exp(coeffb)
       firstx = 1
    endif

    if (t.le.tcxcut) then
       rcxmult = 1.0d0
    else

       rcxmult = coeffa * exp(-coeffb*t)

       if (rcxmult.lt.1.0) rcxmult = 1.0d0
       if (rcxmult.gt.1500.0) rcxmult = 1500.0d0

    endif

    return
  end function rcxmult



  real*8 function cond(s,t)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     Returns the first convective energy component
    !
    real*8 s,t

    !real*8 gamam
    !external gamma

    if (actswcond.eq.0.0) then
       cond = 0.0d0
       return
    endif

    cond = 5.0d0/2.0d0 * gamma(s) * econv * t

    return
  end function cond



  real*8 function conv(s,n,t)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     Calculates kinetic convective term
    !
    real*8 s,n,t

    real*8 tmpgam
    !real*8 gamma
    !external gamma

    if (actswconv.eq.0.0) then
       conv = 0.0d0
       return
    endif

    tmpgam = gamma(s)

    conv = 0.5d0 * mb * mconv * (tmpgam/n)**2 * tmpgam
    !
    !      write(6,'(a6,6g14.6)') 'conv:',mb,mconv,
    !     >                       tmpgam,tmpgam/n,n,conv
    !
    return
  end function conv


  real*8 function pradupdt(s,n,nold,te,teold)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     PRADUPDT: This returns the integrated value of the
    !               radiation losses up to the point s.
    !
    !     This function updates the lasts and lastprad values after
    !     calculating a Prad contribution based on the average temperatures
    !     over the interval.
    !
    real*8 s,n,te,teold,nold

    common /praddata/ lastprad,lasts
    real*8 lastprad,lasts

    real*8 teav,nav

    integer in

    if (actswprad.eq.0.0) then
       pradupdt = 0.0
    elseif (actswprad.eq.1.0) then

       if (s.lt.lenr) then
          pradupdt = lamr*prad0*(1-exp(-s/lamr))
       else
          pradupdt = lamr*prad0*(1-exp(-lenr/lamr))
       endif
    elseif (actswprad.eq.2.0) then

       if (s.eq.soffset.or.s.lt.lasts) then
          lasts = soffset
          lastprad = 0.0
          pradupdt = 0.0
          return
       elseif (s.eq.lasts) then
          pradupdt = lastprad
          return
       endif

       teav = (te+teold)/2.0 / talimp
       teav = max(1.0d-6,teav)

       nav  = (n+nold)/2.0

       pradupdt = lastprad +  alfimp *  nav**2 * 2.0d-31 / (teav**ex1imp+teav**ex2imp) * (s-lasts)
       !
       !        write(6,*) 'pradupdt1:',pradupdt,lastprad,s,lasts,nav,teav
       !
       lastprad = pradupdt
       lasts = s
       !
       !     PRAD = RADSRC_MULT * PINQE
       !
    elseif (actswprad.eq.3.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          pradupdt = (intqe(in) * s / sptscopy(1))
       else
          pradupdt = (( intqe(in-1) + ((intqe(in)-intqe(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif
       !
       !     PRAD = RADSRC_MULT * (EXTERNAL RADIATION SOURCE)
       !
    elseif (actswprad.eq.4.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          pradupdt = (intrad(in) * s / sptscopy(1))
       else
          pradupdt = ( ( intrad(in-1) + ( (intrad(in)-intrad(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif
    elseif (actswprad.eq.6.0) then
       !
       !        re-use lamr and lenr for a rectangular radiation source
       !
       if (s.lt.lamr) then
          pradupdt = 0.0
       elseif (s.lt.lenr) then 
          pradupdt = (s-lamr)/(lenr-lamr)*prad0
       else
          pradupdt = prad0
       endif
    endif

    pradupdt = radsrc_mult * pradupdt

    return
  end function pradupdt



  real*8 function estprad(s,n,te)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     ESTPRAD: This returns the radiative energy loss -
    !     approximately integrated to the current point s.
    !
    !     This function returns an ESTIMATE of the Prad contribution at
    !     the given S-position ... the lastprad and lasts values are
    !     updated after every R-K iteration.
    !
    real*8 s,n,te

    common /praddata/ lastprad,lasts
    real*8 lastprad,lasts

    real*8 teav

    integer in


    if (actswprad.eq.0.0) then
       estprad = 0.0
    elseif (actswprad.eq.1.0) then

       if (s.lt.lenr) then
          estprad = lamr*prad0*(1-exp(-s/lamr))
       else
          estprad = lamr*prad0*(1-exp(-lenr/lamr))
       endif
    elseif (actswprad.eq.2.0) then

       if (s.eq.soffset.or.s.lt.lasts) then
          estprad = 0.0
          return
       elseif (s.eq.lasts) then
          estprad = lastprad
          return
       endif

       teav = te / talimp
       teav = max(1.0d-6,teav)

       estprad = lastprad +  alfimp *  n**2 * 2.0d-31 / (teav**ex1imp+teav**ex2imp) * (s-lasts)
       !
       !     PRAD = PINQE * MULT
       !
    elseif (actswprad.eq.3.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          estprad = (intqe(in) * s / sptscopy(1))
       else
          estprad = ( ( intqe(in-1) + ( (intqe(in)-intqe(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    elseif (actswprad.eq.4.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          estprad = (intrad(in) * s / sptscopy(1))
       else
          estprad = ( ( intrad(in-1) + ( (intrad(in)-intrad(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    elseif (actswprad.eq.6.0) then
       !
       !        re-use lamr and lenr for a rectangular radiation source
       !
       if (s.lt.lamr) then
          estprad = 0.0
       elseif (s.lt.lenr) then 
          estprad = (s-lamr)/(lenr-lamr)*prad0
       else
          estprad = prad0
       endif
    endif
    !
    !     apply rad src multiplier
    !
    estprad = radsrc_mult * estprad

    !if (debug_s22) write(6,'(a,10(1x,g12.5))') 'ESTPRAD:',actswprad,s,n,te,estprad,lastprad,lasts,prad0
    !
    !      write(6,*) 'estprad:',estprad,lastprad,s,lasts,n,te
    !
    return
  end function estprad


  real*8 function peiupdt(s,n,nold,te,ti,teold,tiold,fval)
    use mod_solparams
    use mod_solswitch
    use mod_sol22pei
    use mod_solcommon
    implicit none
    !
    !     This returns the electron-ion energy exchange component -
    !     approximately integrated to the current point s.
    !
    !     This function updates the lasts and lastpei values after
    !     calculating a Pei contribution based on the average temperatures
    !     over the interval.
    !
    real*8 s,n,te,ti,teold,tiold,nold,fval
    !real*8 lnlam
    !external lnlam
    !
    !      common /pei/ lastpei,lasts
    !      real*8 lasts,lastpei
    !
    real*8 tiav,teav,nav

    if (actswpei.eq.0.0) then
       peiupdt = 0.0
       return
    endif

    if (s.eq.soffset.or.s.lt.lasts) then
       lasts = soffset
       lastpei = 0.0
       peiupdt = 0.0
       fval = 0.0
       return
    elseif (s.eq.lasts) then
       peiupdt = lastpei
       if (actswpei.eq.3.0) then
          fval = peiupdt
          peiupdt = 0.0
       endif
       return
    endif

    teav = (te+teold)/2.0
    tiav = (ti+tiold)/2.0
    nav  = (n+nold)/2.0

    peiupdt = lastpei +  peicf * ((1.14e-32 * nav**2 * ( teav - tiav )) / (mb* teav**(1.5))) *(s-lasts) * lnlam(nav,teav) / 15.0
    !
    !      write(6,*) 'peiupdt1:',peiupdt,lastpei,s,lasts,nav,teav,tiav
    !
    lastpei = peiupdt
    lasts = s

    if (actswpei.eq.3.0) then
       fval = peiupdt
       peiupdt = 0.0
    endif
    !
    return
  end function peiupdt



  real*8 function estpei(s,n,te,ti)
    use mod_solparams
    use mod_solswitch
    use mod_sol22pei
    use mod_solcommon
    implicit none
    !
    !     This returns the electron-ion energy exchange component -
    !     approximately integrated to the current point s.
    !
    !     This function returns an ESTIMATE of the Pei contribution at
    !     the given S-position ... the lastpei and lasts values are
    !     updated after every R-K iteration.
    !
    real*8 s,n,te,ti

    !
    !      common /pei/ lastpei,lasts
    !      real*8  lasts,lastpei
    !
    !real*8 lnlam
    !external lnlam

    if (actswpei.eq.0.0.or.actswpei.eq.3.0) then
       estpei = 0.0
       return
    endif
    if (s.eq.soffset.or.s.lt.lasts) then
       estpei = 0.0
       return
    elseif (s.eq.lasts) then
       estpei = lastpei
       return
    endif

    estpei = lastpei + peicf * ((1.14e-32 * n**2 * (te-ti)) / (mb * te**(1.5))) * (s-lasts) * lnlam(n,te) / 15.0
    !
    !      write(6,*) 'estpei:',estpei,lastpei,s,lasts,n,te,ti
    !
    return
  end function estpei



  real*8 function lnlam(n,te)
    use mod_lambda
    implicit none
    real*8 n,te
    !
    !     LNLAM: This function returns the value of Lambda used in the
    !            Pei energy transfer formulae.
    !
    ! The definition of lambda used in the code has been centralized into
    ! the mod_lambda.f90 code in comsrc/f90. Controlled by the optional input
    ! lambda_opt and lambda_val. lambda_val defaults to 15.0 if not input. 

    lnlam = coulomb_lambda(real(n),real(te))


    !if (n.le.0.0.or.te.le.0.0.or.((1.5e13 * te**(1.5) / sqrt(n)).le.1.0)) then
    !   lnlam = 15.0
    !else
    !   lnlam = log(1.5e13 * te**(1.5) / sqrt(n))
    !endif

    return
  end function lnlam



  real*8 function phelpiupdt(s,n,nold,t,told)
    use mod_solparams
    use mod_solswitch
    use mod_sol22phelpi
    use mod_solcommon
    implicit none
    !
    !     Calculates losses to electrons due to hydrogenic cooling
    !
    real*8 s,n,t,told,nold

    !
    !      common /phelpi/ lasts,lastphelp,lastsrc
    !      real*8 lasts,lastphelp,lastsrc
    !
    real*8 tav,nav

    real*8  helpi,src
    !real*8 srci
    !external srci
    integer in

    if (actswphelp.eq.0.0.or.(actswphelp.eq.3.0.and.(.not.pinavail))) then
       phelpiupdt = 0.0
       return
    elseif(actswphelp.eq.1.0.or.(actswphelp.eq.2.0.and.(.not.pinavail))) then

       if (s.eq.soffset.or.s.lt.lasts) then
          lasts = soffset
          lastphelp = 0.0
          lastsrc = 0.0
          phelpiupdt = 0.0
          return
       elseif (s.eq.lasts) then
          phelpiupdt = lastphelp
          return
       elseif (t.lt.tcutqe) then
          phelpiupdt = lastphelp
          return
       endif

       tav = (t+told)/2.0
       nav = (n+nold)/2.0

       helpi = 17.5 + (5.0+37.5/tav)*(1.0+ 0.25/tav)*log10(1e21/nav)
       src = srci(s)

       phelpiupdt = lastphelp + helpi * econv * (src-lastsrc)

       lastphelp = phelpiupdt
       lasts = s
       lastsrc = src

    elseif (actswphelp.eq.2.0.or.actswphelp.eq.3.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          phelpiupdt = (intqe(in) * s / sptscopy(1))
       else
          phelpiupdt = (( intqe(in-1) +( (intqe(in)-intqe(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    endif

    return
  end function phelpiupdt



  real*8 function estphelpi(s,n,t)
    use mod_solparams
    use mod_solswitch
    use mod_sol22phelpi
    use mod_solcommon
    implicit none
    !
    !     Calculates losses to electrons due to hydrogenic cooling
    !     Returns an ESTIMATE at the given s,n,t but does not
    !     update the values for the interval.
    !
    real*8 s,n,t
    !
    !      common /phelpi/ lasts,lastphelp,lastsrc
    !      real*8 lasts,lastphelp,lastsrc
    !
    real*8  helpi,src
    !real*8 srci
    !external srci

    integer in

    if (actswphelp.eq.0.0.or.(actswphelp.eq.3.0.and.(.not.pinavail))) then
       estphelpi = 0.0
       return
    elseif(actswphelp.eq.1.0.or.(actswphelp.eq.2.0.and.(.not.pinavail))) then

       if (s.eq.soffset.or.s.lt.lasts) then
          estphelpi = 0.0
          return
       elseif (s.eq.lasts) then
          estphelpi = lastphelp
          return
       elseif (t.lt.tcutqe) then
          estphelpi = lastphelp
          return
       endif

       helpi = 17.5 + (5.0 + 37.5/t) * (1.0+ 0.25/t) * log10(1e21/n)
       src = srci(s)

       estphelpi = lastphelp + helpi * econv * (src-lastsrc)

    elseif (actswphelp.eq.2.0.or.actswphelp.eq.3.0) then

       call binsearch(s,in)

       if (in.eq.1) then
          estphelpi = (intqe(in) * s / sptscopy(1))
       else
          estphelpi = (( intqe(in-1) +( (intqe(in)-intqe(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    endif

    return
  end function estphelpi



  real*8 function pcxupdt(s,t,told)
    use mod_solparams
    use mod_solswitch
    use mod_sol22pcx
    use mod_solcommon
    implicit none
    !
    !     Calculates the charge exchange loss term for each step. Stores
    !     the approximated integral.
    !
    !     Updates the integral after new Te, Ti are found for the interval.
    !     Uses the average.
    !
    real*8 s,t,told
    !
    !      common /pcx/ lasts,lastpcx,lastsrc
    !      real*8 lasts,lastpcx,lastsrc
    !
    real*8 tav
    real*8   src
    !real*8 srci,pinqid
    !external srci,pinqid

    integer in

    tav = (t+told)/2.0

    if (actswpcx.eq.0.0.or.(actswpcx.eq.3.0.and.(.not.pinavail))) then
       pcxupdt = 0.0
       return
    elseif(actswpcx.eq.1.0.or.((actswpcx.eq.2.0.or.actswpcx.eq.5.0).and.(.not.pinavail))) then

       if (s.eq.soffset.or.s.lt.lasts) then
          lasts = soffset
          lastpcx = 0.0
          lastsrc = 0.0
          pcxupdt = 0.0
          return
       elseif (s.eq.lasts) then
          pcxupdt = lastpcx
          return
       elseif (t.lt.tcutcx) then
          !
          !        Approximation - if new T is below the cutoff then set addition
          !                        for interval to zero.
          !            
          pcxupdt = lastpcx
          return
       endif

       src = srci(s)

       pcxupdt = lastpcx+(1.5*tav) * econv *ceicf * (src-lastsrc)

       lastpcx = pcxupdt
       lasts = s
       lastsrc = src
       return

    elseif (actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.actswpcx.eq.5) then

       call binsearch(s,in)

       if (in.eq.1) then
          pcxupdt = (intqi(in) * s / sptscopy(1))
       else
          pcxupdt = (( intqi(in-1) + ( (intqi(in)-intqi(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    elseif (actswpcx.eq.4.0) then

       pcxupdt = pinqid(s,1)

    endif

    return
  end function pcxupdt



  real*8 function estpcx(s,t)
    use mod_solparams
    use mod_solswitch
    use mod_sol22pcx
    use mod_solcommon
    implicit none
    !
    !     Calculates the charge exchange loss term for each step. Stores
    !     the approximated integral.
    !
    !     Returns an estimate of PCX over an interval. Does not update
    !     overall integral.
    !
    real*8 s,t
    !
    !      common /pcx/ lasts,lastpcx,lastsrc
    !      real*8 lasts,lastpcx,lastsrc
    !
    real*8   src
    !real*8 srci,pinqid
    !external srci,pinqid

    integer in

    if (actswpcx.eq.0.0.or.(actswpcx.eq.3.0.and.(.not.pinavail))) then
       estpcx = 0.0
       return
    elseif(actswpcx.eq.1.0.or.((actswpcx.eq.2.0.or.actswpcx.eq.5).and.(.not.pinavail))) then

       if (s.eq.soffset.or.s.lt.lasts) then
          estpcx = 0.0
          return
       elseif (s.eq.lasts) then
          estpcx = lastpcx
          return
       elseif (t.lt.tcutcx) then
          estpcx = lastpcx
          return
       endif

       src = srci(s)

       estpcx = lastpcx + (1.5*t) * econv *ceicf * (src-lastsrc)

    elseif (actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.actswpcx.eq.5) then
       !
       !        Find cell which particle is in
       !
       call binsearch(s,in)

       if (in.eq.1) then
          estpcx = (intqi(in) * s / sptscopy(1))
       else
          estpcx = ( ( intqi(in-1) + ( (intqi(in)-intqi(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif

    elseif (actswpcx.eq.4.0) then

       estpcx = pinqid(s,0)

    endif

    return
  end function estpcx



  real*8 function gamma(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     Sets the value of the flux = nv = n0v0- int[0 to s] (S(s))
    !     Goes to zero at the end of the ionization source. If using
    !     a normalized ionization source.
    !
    real*8 s

    real*8 rpos
    !real*8 srcf,majrpos,srcrec
    !external srcf,majrpos,srcrec
    !
    !      if (s.gt.ssrcfi) then
    !         gamma = 0.0
    !         return
    !      endif
    !
    if (actswmajr.eq.4.0) then
       rpos = majrpos(s)
       gamma = ( gamma0 * r0init + srcf(s) - srcrec(s)) / rpos
    else
       gamma = gamma0 + srcf(s) - srcrec(s)
    endif

    !if (debug_s22) write(6,'(a,10(1x,g12.5))') 'GAMMA:',s,gamma,gamma0,srcf(s),srcrec(s)
    !write(6,'(a,10(1x,g12.5))') 'GAMMA:',s,gamma,gamma0,srcf(s),srcrec(s)


    return
  end function gamma



  real*8 function srcrec(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     SRCREC: This rotuine returns the value of the integrated
    !             recombination particle source to the point s.
    !
    !             At this time it is only supported for PIN
    !             iterated plasma calculations. It will only work
    !             correctly in conjunction with ionization
    !             options 2 or 8 and with perpendicular flux
    !             option 2.
    !
    integer in

    if (actswrecom.eq.0.0.or.(actswrecom.eq.1.and.(.not.pinavail))) then
       srcrec = 0.0
       return
    elseif ((pinavail.and.actswrecom.eq.1.0).or.actswrecom.eq.2) then
       !
       !        Search for right cell
       !
       call binsearch(s,in)

       if (in.eq.1) then
          srcrec =  intrecsrc(in) * s / sptscopy(1)
       else
          srcrec = intrecsrc(in-1) +( (intrecsrc(in)-intrecsrc(in-1) ) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))
       endif

    endif

    return
  end function srcrec



  real*8 function srcf(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    ! NOTE: WARNING: srcf(s) uses the integration array intionsrc while srci(s) intioniz
    !                MOST of the time these should be the same. However, they are calculated
    !                at different points in the code and careful consideration is needed
    !                in trying to simplify the code. The best approach to this would likely be
    !                to completely eliminate the intioniz array and replace it with intionsrc.
    !
    !
    !     This returns the integral over the source flux function
    !     from 0 to s. Includes Ionization and cross-field losses.
    !     Returns the values from the array intionsrc.
    !
    real*8 s

    integer i,in

    srcf = 0.0
    !
    !     For S=0 - integration=0 so return
    !
    if (s.eq.0.0) return

    if (actswmajr.eq.4.0) then
       !
       !        Major Radius corrected ion source
       !
       if (s.gt.ssrcfi) then
          srcf = pnormfact * r0init
          return
       endif
       !
       !       ALL ionization options for major radius correction
       !       have been pre-integrated because of the R(s) factor
       !       in the integral - which depends on the grid. The
       !       Cross-field or RCONST/GPERPCOR term has been included during
       !       the pre-integration.
       !
       !        Search for right cell
       !
       call binsearch(s,in)

       if (in.eq.1) then
          srcf = fnorm * (intionsrc(in) * s / sptscopy(1))
       else
          srcf = fnorm*(( intionsrc(in-1)+( (intionsrc(in)-intionsrc(in-1))*(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif
       !
       !     Regular treatment
       !
    else

       if (actswion.eq.0.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.0.0 .and.(.not.pinavail))) then

          srcf = expsrc(s) + gperpf(s)

       elseif (actswion.eq.3.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.3.0.and.(.not.pinavail))) then
          !
          !          Add in triangular ionization source
          !
          srcf = trisrc(s)  +  gperpf(s)

       elseif (actswion.eq.4.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.4.0.and.(.not.pinavail))) then
          !
          !          Add in rectangular ionization source
          !
          srcf = rectsrc(s) +  gperpf(s)

       elseif (actswion.eq.6.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.6.0.and.(.not.pinavail))) then
          !
          !          Add in s5gauss ionization source
          !
          srcf = s5gauss(s) +  gperpf(s)

       elseif (actswion.eq.9.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.9.0.and.(.not.pinavail))) then
          !
          !          Add in s5gauss2 ionization source
          !
          srcf = s5gauss2(s) +  gperpf(s)

       elseif (((actswion.eq.1.0.or.actswion.eq.2.0).and.pinavail) &
            .or.((actswioni.eq.11.or.actswioni.eq.15).and.(.not.pinavail))) then
          !
          !          Search for right cell
          !
          call binsearch(s,in)

          if (in.eq.1) then
             srcf = fnorm * (intionsrc(in) * s / sptscopy(1)) + gperpf(s)
          else
             srcf = fnorm *( ( intionsrc(in-1) + ( (intionsrc(in)-intionsrc(in-1)) *  &
                  (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))) + gperpf(s)
          endif
       elseif (actswion.eq.16.or.(actswioni.eq.16.and.(.not.pinavail))) then 
          !
          ! Double source + gperp
          !
          srcf = dblsrc(s) + gperpf(s)

       else
          !
          !          Code has reached an error condition and should stop.
          !
          write (6,'(a,2(1x,g12.5),l6)') 'ERROR in SOLASCV:'//&
               ' Invalid Ionization Source Options: ' ,actswion,&
               actswioni,pinavail
          write (0,'(a,2(1x,g12.5),l6)') 'ERROR in SOLASCV:'//&
               ' Invalid Ionization Source Options: ' ,actswion,&
               actswioni,pinavail
          stop 'SOL22: Invalid Ionizaition Source Option'
       endif
       !
       !     Endif for swmajr
       !
    endif

    !if (debug_s22)  write(6,'(a,3(1x,g12.5))') 'SRCF:',s,srcf

    return
  end function srcf



!  real*8 function srcf(s)
!    use mod_solparams
!    use mod_solswitch
!    use mod_solcommon
!    implicit none
!    !
!    !     This returns the integral over the source flux function
!    !     from 0 to s. Includes Ionization and cross-field losses.
!    !     Returns the values from the array intionsrc.
!    !
!    real*8 s
!
!    integer i,in
!
!    srcf = 0.0
!    !
!    !     For S=0 - integration=0 so return
!    !
!    if (s.eq.0.0) return
!
!
!    srcf = srci(s) + gperpf(s)
!
!    write(6,'(a,10(1x,g12.5))') 'SRCF:',srci(s),gperpf(s)
!
!    !if (debug_s22)  write(6,'(a,3(1x,g12.5))') 'SRCF:',s,srcf
!
!    return
!  end function srcf

  


  real*8 function srci(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     This returns the integral over ONLY the ionization source
    !     from 0 to s. Includes just Ionization.
    !     Returns values from the array intioniz.
    !
    real*8 s

    integer i,in
    !real*8 expsrc,trisrc,rectsrc,s5gauss,s5gauss2
    !external expsrc,trisrc,rectsrc,s5gauss,s5gauss2

    srci = 0.0
    !
    !     For S=0 - integration=0 so return
    !
    if (s.eq.0.0) return

    if (actswmajr.eq.4.0) then
       !
       !       Major Radius corrected ion source
       !
       if (s.gt.ssrcfi) then
          srci = pnormfact * r0init
          return
       endif
       !
       !       ALL ionization options for major radius correction
       !       have been pre-integrated because of the R(s) factor
       !       in the integral - which depends on the grid. The
       !       Cross-field or RCONST/GPERPCOR term has been included during
       !       the pre-integration.
       !
       !        Search for right cell
       !
       call binsearch(s,in)

       if (in.eq.1) then
          srci = fnorm2 * (intioniz(in) * s / sptscopy(1))
       else
          srci = fnorm2 *(( intioniz(in-1) + ( (intioniz(in)-intioniz(in-1)) * &
               (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif


    else

       !
       !     Regular treatment
       !
   
       if (actswion.eq.0.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.0.0.and.(.not.pinavail))) then

          srci = expsrc(s)

       elseif (actswion.eq.3.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.3.0.and.(.not.pinavail))) then
          !
          !          Add in triangular ionization source
          !
          srci = trisrc(s)

       elseif (actswion.eq.4.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)  &
            .and.actswioni.eq.4.0.and.(.not.pinavail))) then
          !
          !          Add in rectangular ionization source
          !
          srci = rectsrc(s)

       elseif (actswion.eq.6.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.6.0.and.(.not.pinavail))) then
          !
          !          Add in gaussian ionization source
          !
          srci = s5gauss(s)

       elseif (actswion.eq.9.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0) &
            .and.actswioni.eq.9.0.and.(.not.pinavail))) then
          !
          !          Add in gaussian squared ionization source
          !
          srci = s5gauss2(s)

       elseif (((actswion.eq.1.0.or.actswion.eq.2.0).and.pinavail) &
            .or.((actswioni.eq.11.or.actswioni.eq.15).and.(.not.pinavail))) then
          !
          !          Search for right cell
          !
          call binsearch(s,in)

          if (in.eq.1) then
             srci = fnorm2 * (intioniz(in) * s / sptscopy(1))
          else
             srci = fnorm2*((intioniz(in-1) + ( (intioniz(in)-intioniz(in-1))*(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
          endif

       elseif (actswion.eq.16.or.(actswioni.eq.16.and.(.not.pinavail))) then 
          srci = dblsrc(s)
       else
          !
          !         Code has reached an error condition and should stop.
          !
          write (6,'(a,2(1x,g12.5),l6)') 'ERROR in SOL22 Input:'//&
               ' Invalid Ionization Source Options: ' ,actswion,&
               actswioni,pinavail
          write (0,'(a,2(1x,g12.5),l6)') 'ERROR in SOL22 Input:'//&
               ' Invalid Ionization Source Options: ' ,actswion,&
               actswioni,pinavail
          stop 'SOL22: Invalid Ionizaition Source Option'
       endif
       !
       !     Endif for swmajr
       !
    endif

    return
  end function srci



  real*8 function expsrc(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     EXPSRC: This function returns the integral of
    !     the exponential ionization source from soffset to s.
    !
    if (s.gt.ssrcfi) then
       expsrc = pnormfact
    else
       expsrc = s0 * (ssrcdecay*(1.0-exp(-(s-soffset)/ssrcdecay)))
    endif

    return
  end function expsrc



  real*8 function trisrc(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     TRISRC: This function returns the integral of
    !     the triangular ionization source from soffset to s.
    !
    if (s.lt.ssrcst) then
       trisrc = 0.0
    elseif (s.gt.ssrcfi) then
       trisrc = pnormfact
    elseif (s.lt.ssrcmid) then
       trisrc = s0 * ( 0.5 * (s-ssrcst)**2)
    elseif (s.le.ssrcfi) then
       trisrc = pnormfact - s0 * (0.5 * (ssrcfi-s)**2)
    endif

    return
  end function trisrc



  real*8 function rectsrc(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     RECTSRC: This function returns the integral of
    !     the rectangular ionization source from soffset to s.
    !
    if (s.lt.ssrcst) then
       rectsrc = 0.0
    elseif (s.gt.ssrcfi) then
       rectsrc = pnormfact
    elseif (s.le.ssrcfi) then
       rectsrc = s0 * (s-ssrcst)
    endif

    return
  end function rectsrc

  real*8 function rectsrc2(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     RECTSRC2: This function returns the integral of
    !     the rectangular ionization source from soffset to s.
    !
    if (s.lt.ssrcst2) then
       rectsrc2 = 0.0
    elseif (s.gt.ssrcfi2) then
       rectsrc2 = pnormfact
    elseif (s.le.ssrcfi2) then
       rectsrc2 = s02 * (s-ssrcst2)
    endif

    return
  end function rectsrc2


  
  real*8 function dblsrc(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     DBLSRC: This function returns the integral of
    !     an ionization source consisting of two separate sources
    !     weighted by dblsrc_frac and 1-dbsrc_frac
    !
    if (dblsrc_opt.eq.0) then ! exp+rect
       dblsrc = dblsrc_frac * expsrc(s) + (1.0-dblsrc_frac) * rectsrc2(s)
    elseif (dblsrc_opt.eq.1) then ! rect+rect
       dblsrc = dblsrc_frac * rectsrc(s) + (1.0-dblsrc_frac) * rectsrc2(s)
    endif

    return
  end function dblsrc


  real*8 function s5gauss(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     S5GAUSS:  This function returns the integral of
    !     an s**5 * exp(-s) ionization source from soffset to s.
    !
    real*8 eas2,stmp

    stmp = s - soffset

    if (stmp.lt.ssrcst) then
       s5gauss = 0.0
    elseif (stmp.gt.ssrcfi) then
       s5gauss = pnormfact
    elseif (stmp.le.ssrcfi) then
       eas2 = exp(-s5alph * stmp**2)
       s5gauss = s0 * ( 1.0 / s5alph3 -0.5 * stmp**4/s5alph * eas2 - stmp**2/s5alph2 * eas2 - eas2/s5alph3)
    endif

    return
  end function s5gauss



  real*8 function s5gauss2(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     S5GAUSS2:  This function returns the integral of
    !     an s**5 * exp(-s) ionization source from soffset to s shifted
    !     from the origin by a specified value (s5offset).
    !
    real*8 eas2,stmp

    stmp = s - soffset + s5offset

    if (stmp.lt.ssrcst) then
       s5gauss2 = 0.0
    elseif (stmp.gt.ssrcfi) then
       s5gauss2 = pnormfact
    elseif (stmp.le.ssrcfi) then
       eas2 = exp(-s5alph * stmp**2)
       s5gauss2 = s0 * ( s5startval -0.5 * stmp**4/s5alph * eas2 - stmp**2/s5alph2 * eas2 - eas2/s5alph3)
    endif

    return
  end function s5gauss2



  subroutine preint(startn,npts,spts,src,intsrc,srcsum,ringlen,flage2d,flagmajr,sbnd,rbnd,gperpn)
    use mod_solparams
    implicit none

    integer npts,startn
    real*8 spts(mxspts),src(mxspts),intsrc(mxspts),srcsum
    real*8 ringlen
    real*8 sbnd(0:mxspts),rbnd(0:mxspts),gperpn
    real*8 asum
    real    flage2d,flagmajr
    !
    !     PREINT: This subroutine pre-integrates a numerical source
    !             spread along the ring.
    !
    !             intionsrc(S) = INTEGRAL(0 to S) of ionsrc(S) + gperp(s)
    !
    !             The values are supplied at each grid point represented
    !             by sptscopy(i) - and are linearly interpolated.
    !
    !             This routine can be used to integrate any source
    !             array.
    !
    !             Flage2d = 0.0  - do full integration
    !                     = >=1.0  - Ignore first and last half cell
    !
    !
    integer i,ik
    real*8 rmean1,rmean2
    !
    !     Do integral including MAJOR RADIUS
    !
    !     Change to integration routine to support Major Radius
    !     Not quite the same as previously since it uses actual
    !     cell boundaries.
    !

    srcsum = 0.0
    intsrc = 0.0
    
    do ik = startn,npts

       if (flagmajr.eq.4.0) then
          rmean1 = (rbnd(ik-1)+(rbnd(ik)+rbnd(ik-1))/2.0)/2.0
          rmean2 = (rbnd(ik)+(rbnd(ik)+rbnd(ik-1))/2.0)/2.0
       else
          rmean1 = 1.0
          rmean2 = 1.0
       endif

       if (.not.(ik.eq.startn.and.flage2d.ne.0.0)) then
          srcsum = srcsum + (src(ik)+gperpn) * (spts(ik)-sbnd(ik-1))*rmean1
       endif
       !
       !           Do in two parts to get integral at cell centre.
       !

       intsrc(ik) = srcsum

       !
       !           Does not need startn adjustment because integral is
       !           only over 1/2 of the ring.
       !
       if(.not.(ik.eq.npts.and.flage2d.gt.0.0)) then
          srcsum = srcsum + (src(ik)+gperpn) * (sbnd(ik)-spts(ik))*rmean2
       endif
    end do

    intsrc(npts+1) = srcsum

    return
  end subroutine preint


  real*8 function estppress(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     Distribute the PP target pressure loss over the ring
    !
    real*8 s
    !
    !     ESTPPRESS: This routine calculates the amount of pressure transferred
    !                to the private plasma and then distributes this 
    !                as a function of S from the
    !                current main SOL ring.
    !
    !                Since everything is analytic at this point and does not
    !                depend on the changing plasma conditions for any of
    !                it's options - it does not have an update method as the
    !                other options may require.
    !
    !
    estppress = 0.0
    !
    !     Option is OFF
    !
    if (actswppress.eq.0.0) then

       estppress = 0.0
       !
       !     Power added below X-point
       !
    elseif (actswppress.eq.1.0) then

       if (s.lt.sxp) then
          estppress  = s/sxp * pp_press
       else
          estppress  = pp_press
       endif
       !
       !     Power added to specified distance along field line
       !
    elseif (actswppress.eq.2.0) then

       if (s.lt.pp_pow_dist*ringlen) then
          estppress  = s/(pp_pow_dist*ringlen) * pp_press
       else
          estppress  = pp_press
       endif

    endif
    !
    !     End of routine
    !
    return
  end function estppress





  !
  !     WF'95: DISTRIBUTE POWER ALONG S
  !
  real*8 function paes(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none

    real*8 s
    !real*8 majrpos,areaint
    !external majrpos,areaint

    real*8 rfact,rpos


    paes = pae

    if (actswmajr.eq.4.0) then
       rpos = majrpos(s)
       rfact = r0init/rpos
    else
       rfact = 1.0
    endif

    if (actswpow.eq.0.0) then
       ! total power being carried includes radiated power as well as target power - however - only works with some prad options
       !paes =  rfact * (pae + totprad)
       paes =  rfact * pae
    elseif (actswpow.eq.1.0) then
       ! total power being carried includes radiated power as well as target power - however - only works with some prad options
       paes = pae * (1.0 - s/(halfringlen))
       !paes = (pae+totprad) * (1.0 - s/(halfringlen))
    elseif (actswpow.eq.2.0) then
       ! total power being carried includes radiated power as well as target power - however - only works with some prad options
       if (s.le.sxp) then
          paes = (pae+totprad)
       else
          paes = (pae+totprad) * (1.0 - (s-sxp)/(halfringlen-sxp))
       endif
    elseif (actswpow.eq.3.0) then
       paes =  rfact * pae - pinpute * areaint(s) / rpos
    elseif (actswpow.eq.4.0) then
       paes = rfact* (pae - (pae + pae_end)*s/ringlen)
    elseif (actswpow.eq.5.0) then
       paes = rfact*(pae * (1.0 - s/(halfringlen)) - qesum * s / (halfringlen))
    elseif (actswpow.eq.6.0) then
       paes = rfact* ((pae - (pae + pae_end)*s/ringlen) - qesum * s / (halfringlen))
    elseif (actswpow.eq.7.0) then
       if (s.le.spow) then
          paes = pae
       else
          paes = pae * (1.0 - (s-spow)/(halfringlen-spow))
       endif
    elseif (actswpow.eq.8.0) then
       if (s.le.spow) then
          paes = pae
       else
          paes = rfact*(pae * (1.0 - (s-spow)/(halfringlen-spow)) - qesum * (s-spow) / (halfringlen-spow))
       endif
    elseif (actswpow.eq.9.0) then
       if (s.le.spow) then
          paes = pae
       elseif (s.le.spow2) then
          paes = rfact*(pae * (1.0 - (s-spow)/(spow2-spow)))
       else
          paes = 0.0
       endif
    elseif (actswpow.eq.10.0) then
       if (s.le.spow) then
          paes = pae
       elseif (s.le.spow2) then
          paes = rfact*(pae * (1.0 - (s-spow)/(spow2-spow)) - qesum * (s-spow) / (spow2-spow))
       else
          paes = 0.0
       endif

    elseif (actswpow.eq.11.0) then

       if (s.le.spow) then
          paes = pae
       else
          paes = rfact* (pae - (pae + pae_end)*(s-spow)/(ringlen-2.0*spow))
       endif

    endif

    return
  end function paes

  real*8 function qperpe(s,n,te,ti)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 :: s,n,te,ti
    
    if (actswqperpe.eq.0.0) then
       qperpe = 0.0
    elseif (actswqperpe.eq.1.0) then 
       qperpe = -(s/halfringlen)* (estprad(s,n,te) + estphelpi(s,n,te) + estpei(s,n,te,ti) + estppelec(s))
    endif
  end function qperpe
  
  real*8 function qperpi(s,n,te,ti)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 :: s,n,te,ti
    
    if (actswqperpi.eq.0.0) then
       qperpi = 0.0
    elseif (actswqperpi.eq.1.0) then 
       qperpi = -(s/halfringlen)* (estpcx(s,ti) - estpei(s,n,te,ti) + estppion(s))
    endif
    
  end function qperpi

  real*8 function pais(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !real*8 majrpos,areaint
    !external majrpos,areaint

    real*8 rfact,rpos

    pais = pai

    if (actswmajr.eq.4.0) then
       rfact = r0init/majrpos(s)
    else
       rfact = 1.0
    endif

    if (actswpow.eq.0.0) then
       pais =  rfact * pai
    elseif (actswpow.eq.1.0) then
       pais = pai * (1.0 - s/(halfringlen))
    elseif (actswpow.eq.2.0) then
       if (s.le.sxp) then
          pais = pai
       else
          pais = pai * (1.0 - (s-sxp)/(halfringlen-sxp))
       endif
    elseif (actswpow.eq.3.0) then
       pais =  rfact * pai -pinputi*areaint(s)/majrpos(s)
    elseif (actswpow.eq.4.0) then
       pais = rfact* (pai - (pai + pai_end)*s/ringlen)
    elseif (actswpow.eq.5.0) then
       pais = rfact*(pai * (1.0 - s/(halfringlen)) - qisum * s / (halfringlen))
    elseif (actswpow.eq.6.0) then
       pais = rfact* ((pai - (pai + pai_end)*s/ringlen) - qisum * s / (halfringlen))
    elseif (actswpow.eq.7.0) then
       if (s.le.spow) then
          pais = pai
       else
          pais = pai * (1.0 - (s-spow)/(halfringlen-spow))
       endif
    elseif (actswpow.eq.8.0) then
       if (s.le.spow) then
          pais = pai
       else
          pais = rfact*(pai * (1.0 - (s-spow)/(halfringlen-spow)) - qisum * (s-spow) / (halfringlen-spow))
       endif
    elseif (actswpow.eq.9.0) then
       if (s.le.spow) then
          pais = pai
       elseif (s.le.spow2) then

          pais = rfact*(pai * (1.0 - (s-spow)/(spow2-spow)))
       else
          pais = 0.0
       endif
    elseif (actswpow.eq.10.0) then
       if (s.le.spow) then
          pais = pai
       elseif (s.le.spow2) then
          pais = rfact*(pai * (1.0 - (s-spow)/(spow2-spow)) - qisum * (s-spow) / (spow2-spow))
       else
          pais = 0.0
       endif
    elseif (actswpow.eq.11.0) then

       if (s.le.spow) then
          pais = pai
       else
          pais = rfact* (pai - (pai + pai_end)*(s-spow)/(ringlen-2.0*spow))
       endif

    endif

    return
  end function pais
  !
  !     Distribute the PP electron target power loss over the ring
  !
  real*8 function estppelec(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     ESTPPELEC: This routine calculates the amount of power transferred
    !                to the private plasma as a function of S from the
    !                current main SOL ring.
    !
    !                Since everything is analytic at this point and does not
    !                depend on the changing plasma conditions for any of
    !                it's options - it does not have an update method as the
    !                other options may require.
    !
    !
    estppelec = 0.0
    !
    !     Option is OFF
    !
    if (actswppelec.eq.0.0) then

       estppelec = 0.0
       !
       !     Power added below X-point
       !
    elseif (actswppelec.eq.1.0) then

       if (s.lt.sxp) then
          estppelec  = s/sxp * ppelecpow
       else
          estppelec  = ppelecpow
       endif
       !
       !     Power added to specified distance along field line
       !
    elseif (actswppelec.eq.2.0) then

       if (s.lt.pp_pow_dist*ringlen) then
          estppelec  = s/(pp_pow_dist*ringlen) * ppelecpow
       else
          estppelec  = ppelecpow
       endif

    endif
    !
    !     End of routine
    !
    return
  end function estppelec

  !
  !     Distribute the PP ion target power loss over the ring
  !
  real*8 function estppion(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !
    !     ESTPPION:  This routine calculates the amount of power transferred
    !                to the private plasma as a function of S from the
    !                current main SOL ring.
    !
    !                Since everything is analytic at this point and does not
    !                depend on the changing plasma conditions for any of
    !                it's options - it does not have an update method as the
    !                other options may require.
    !
    !
    estppion = 0.0
    !
    !     Option is OFF
    !
    if (actswppion.eq.0.0) then

       estppion = 0.0
       !
       !     Power added below X-point
       !
    elseif (actswppion.eq.1.0) then

       if (s.lt.sxp) then
          estppion  = s/sxp * ppionpow
       else
          estppion  = ppionpow
       endif
       !
       !     Power added to specified distance along field line
       !
    elseif (actswppion.eq.2.0) then

       if (s.lt.pp_pow_dist*ringlen) then
          estppion  = s/(pp_pow_dist*ringlen) * ppionpow
       else
          estppion  = ppionpow
       endif

    endif
    !
    !     End of routine
    !
    return
  end function estppion






  subroutine initlen
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    !     INITLEN: This subroutine initializes the source lengths for
    !             the ionization and radiation sources.


    implicit none

    !        Absolute source distances

    if (lensind.eq.0) then

       if (lensst.ge.lensfi) lensst = 0.0
       ssrcst = min(lensst,halfringlen)

       ssrcfi = min(lensfi,halfringlen)

       if (ssrcst.ge.ssrcfi) ssrcst = 0.0
       ssrclen = (ssrcfi -ssrcst)

       ssrcmid = (ssrcfi + ssrcst) / 2.0

       ssrcdecay = min(lams,halfringlen)

       !        Relative Source distances


       !        Set up the ionization source limits for options where it is
       !        required.

    elseif (lensind.eq.1) then
       lensst = min(lensst,0.5d0)

       lensfi = min(lensfi,0.5d0)

       if (lensst.ge.lensfi) lensst = 0.0
       ssrcst = lensst * ringlen

       ssrcfi = lensfi * ringlen

       if (ssrcst.ge.ssrcfi) ssrcst = 0.0
       ssrclen = (ssrcfi -ssrcst)

       ssrcmid = (ssrcfi + ssrcst) / 2.0

       ssrcdecay = min(lams*ringlen,halfringlen)

       !     If using PIN data need to set the length of the ionization
       !     source to be the entire 1/2 ring - since the ionization
       !     will fall to zero at whatever point the PIN data
       !     stipulates - not an arbitrary point.

    endif
    if (pinavail.and.(actswion.eq.1.0.or.actswion.eq.2.0)) then
       ssrcst = 0.0d0
       ssrcfi = halfringlen
       ssrclen = (ssrcfi -ssrcst)
       ssrcmid = (ssrcfi + ssrcst) / 2.0

    endif

    !
    ! Set alternate values if swion=16 - analytic dual profile particle source
    !
    if (actswion.eq.16) then

       ! first source
       if (dblsrc_opt.eq.0) then ! exp+rect
          ssrcst = 0.0
          ssrcdecay = min(dble(dblsrc1_p1) * ringlen,halfringlen)
          ssrcfi    = min(dble(dblsrc1_p2) * ringlen, halfringlen)
       else ! rect + rect
          ssrcst = min(dble(dblsrc1_p1) * ringlen, halfringlen)
          ssrcfi = min(dble(dblsrc1_p2) * ringlen, halfringlen)
          if (ssrcst.ge.ssrcfi) ssrcst = 0.0
       endif
       ssrclen = ssrcfi-ssrcst
       ssrcmid = (ssrcfi+ssrcst)/2.0
       
       ! Second source
       ssrcst2 = min(dble(dblsrc2_p1) * ringlen, halfringlen)
       ssrcfi2 = min(dble(dblsrc2_p2) * ringlen, halfringlen)
       if (ssrcst2.ge.ssrcfi2) ssrcst2 = 0.0
       ssrclen2 = ssrcfi2-ssrcst2
       ssrcmid2 = (ssrcfi2+ssrcst2)/2.0
   
    endif

    
    !     Moved to assign_radiation_parameters
    !
    !     Adjust Radiation Source length
    !      lenr = min(lenri,halfringlen)

    !write (6,*) 'Lensrc:',lensst,lensfi,lenri,lenr
    !write (6,*) 'Ssrc  :',ssrcst,ssrcfi,ssrcmid,ssrclen,ringlen

    return



  end subroutine initlen



  real*8 function gperpf(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    !     GPERPF: This function returns the integrated perpendicular
    !             component of the flux from zero to S.

    real*8 s
    integer in
    real*8 gperp_add


    in = 0

    
    if (actswgperp.eq.0.0) then
       gperpf = 0.0
    elseif (actswgperp.eq.1.0.or.actswgperp.eq.2.0) then
       gperpf = gperpcor * (s-soffset)

       !        The overall correction factor is non-zero if there
       !        is a perpendicular flux component.

    elseif (actswgperp.eq.3.0.or.actswgperp.eq.4.0.or.actswgperp.eq.7.or.actswgperp.eq.8.0) then

       !           Search for right cell

       if (gperpcor.ne.0.0) then

          call binsearch(s,in)
          if (in.eq.1) then
             gperpf =  (intgperp(in) * s / sptscopy(1))
          else
             gperpf = ( intgperp(in-1) +( (intgperp(in)-intgperp(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))

             !        Set to zero in alternate case.

          endif

       else

          gperpf = 0.0

       endif

       !        Options 5 and 6 are the sum of two sources -
       !        One square and one uniform
       !        Option 5 - 1/2 ring balance
       !        Option 6 - Whole ring balance

    elseif (actswgperp.eq.5.0.or.actswgperp.eq.6.0) then
       if (s.lt.sgperpbeg) then
          gperpf = gperpcor * (s-soffset)
       elseif (s.ge.sgperpbeg.and.s.le.sgperpend) then
          gperpf = gperpcor * (s-soffset) +gperpcor2 * (s-sgperpbeg)
       elseif (s.gt.sgperpend) then
          gperpf = gperpcor * (s-soffset) +gperpcor2 * (sgperpend-sgperpbeg)
       endif
    endif

    !     ADD in additional Gperp Source/sink terms


    if (switch(swextra).gt.0.0) then

       !        Additional source

       gperp_add = 0.0

       if (s.gt.start_gextra_src.and.s.lt.stop_gextra_src) then

          gperp_add = gperp_add +gextra_src * ( s-start_gextra_src)

       elseif (s.ge.stop_gextra_src) then

          if (start_gextra_src.ne.stop_gextra_src) then

             gperp_add = gperp_add + gextra_src *(stop_gextra_src-start_gextra_src)

          else

             gperp_add = gperp_add + gextra_src

          endif

       endif


       !        Additional Sink
       if (s.gt.start_gextra_sink.and.s.lt.stop_gextra_sink) then

          gperp_add = gperp_add +gextra_sink * ( s-start_gextra_sink)

       elseif (s.ge.stop_gextra_sink) then

          if (start_gextra_src.ne.stop_gextra_src) then

             gperp_add = gperp_add + gextra_sink *(stop_gextra_sink-start_gextra_sink)

          else

             gperp_add = gperp_add + gextra_sink

          endif

       endif


       !        Modify net cross-field flux


       !     Endif for extra source

       gperpf = gperpf + gperp_add

       !      if (debug_s22) then
       !         write(6,'(a,g12.5,g12.5)') 'GPERPF:',s,gperpf
       !      endif

    endif

    return
  end function gperpf



  real*8 function nhs_s(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    !use mod_params
    !use mod_cadas

    !     WF'96: This returns the value of neutral density interpolated,
    !     for s from
    !     ( swnmom = 9 ) =>  nhs(ik), total neutral density
    !     ( swnmom = 10 ) =>  nhs0(ik), ie. the stagnant (primary)
    !                         neutral density.
    !     NB: the first cell density is taken as constant

    implicit none

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'
    !     include 'params'
    !     include 'cadas'

    real*8 s

    integer i,in,bot,top,mid
    if (.not.pinavail) then
       write(6,*) 'nh = ?, PIN not available'
       stop

    endif

    !        Search for right cell

    if (pinavail) then

       call binsearch(s,in)
       if (switch(swnmom).eq.9.0) then
          if (in.eq.1) then
             nhs_s = nhs(1)
          else
             nhs_s = 1.0 *( nhs(in-1) +( (nhs(in) - nhs(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))

             !        write(6,*) 's/spt,nhs/nhs',s/sptscopy(in),nhs_s/nhs(in)

          endif
       elseif (switch(swnmom).eq.10.0) then
          if (in.eq.1) then
             nhs_s = nhs0(1)
          else
             nhs_s = 1.0 *( nhs0(in-1) +( (nhs0(in) - nhs0(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))

             !        write(6,*) 's/spt,nhs/nhs0',s/sptscopy(in),nhs_s/nhs0(in)

          endif

       endif

       !       Code has reached an error condition and should stop.

    else
       write (6,*) 'ERROR in SOLASCV:'//' Invalid neutral density (s<0)'
       stop

    endif
    return


  end function nhs_s

  real*8 function ths_s(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    !use mod_params
    !use mod_cadas

    !     WF'96: This returns the value of neutral temperature interpolated,
    !     for s from ths(ik).

    !     NB: the first cell value is taken as constant

    implicit none

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'
    !     include 'params'
    !     include 'cadas'

    real*8 s

    integer i,in,bot,top,mid
    if (.not.pinavail) then
       write(6,*) 'th = ?, PIN not avail'
       stop 'function ths_s'

    endif

    !        Search for right cell

    if (pinavail) then

       call binsearch(s,in)
       if (in.eq.1) then
          ths_s = ths(in)
       else
          ths_s = 1.0 *( ths(in-1) +( (ths(in) - ths(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))

          !        write(6,*) 's/spt,ths/ths',s/sptscopy(in),ths_s/ths(in)

       endif

       !       Code has reached an error condition and should stop.

    else
       write (6,*) 'ERROR in SOLASCV:'//' Invalid neutral density (s<0)'
       stop 'mod_sol22_sources:ths_s:invalid neutral density'

    endif
    return



  end function ths_s

  real*8 function scx_s(s,n,Ti)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use mod_cx

    !     WF'96: returns the charge exchange mom. loss term (Pa/m)

    implicit none
    real*8 s,n,Ti

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    real*8  sigmavcx,sigma,vp_tmp,Eh_tmp,nh_tmp,ga_tmp,wn_iso,mass_tmp

    !real*8 nhs_s,Ths_s,gamma
    !external nhs_s,Ths_s,gamma

    !real get_bg_mass
    !external get_bg_mass
    ! within sol22 - the background plasma mass is assigned to the
    ! variable mb in the interface

    !     Initialize mass to mass of bg plasma ions

    !mass_tmp = get_bg_mass()
    mass_tmp = mb

    if (.not.pinavail) then
       scx_s = 0.0d0
       return

       !      wn_iso = 1.0/sqrt(3.0d0)

    endif
    wn_iso = 0.57735d0
    Eh_tmp = 1.5 * Ths_s(s)
    nh_tmp = nhs_s(s)

    ga_tmp = - gamma(s)
    if (n.gt.0.0) then
       vp_tmp = abs(ga_tmp) / n
    else
       vp_tmp = abs(ga_tmp) / n

    endif
    if (ga_tmp.eq.0.0.or.nh_tmp.eq.0.0) then
       scx_s = 0.0d0
       return

    endif

    call cxsig(1.0d0, mass_tmp, Eh_tmp, wn_iso, wn_iso,wn_iso, Ti, vp_tmp, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1,sigma, sigmavcx)

    !      if (s.ge.39.0) then
    !      write(6,*) 'nh,th:', nh_tmp,Eh_tmp
    !      write(6,*) 'ga,vp:', ga_tmp,vp_tmp
    !      write(6,*) 'n,Ti:',n,Ti
    !      write(6,*) 'scx_s:',scx_s
    !      write(6,*) '---'
    !      end if

    !      if (s.ge.39.0) then
    !        write(6,*) 's>39 in scx_s: stop'
    !        stop
    !      endif

    scx_s =  mb * mconv * ga_tmp * nh_tmp * sigmavcx
    return

    !

  end function scx_s



  real*8 function scxupdt(s,n,nold,Ti,Tiold)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     WF'96: This returns the charge exchange momentum loss -
    !     approximately integrated to the current s.

    !     This function updates lasts and lastscx values after
    !     calculating a Scx contribution based on the average values
    !     over the interval.

    implicit none

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'
    real*8 s,n,nold,Ti,Tiold

    common /scx/ lastscx,lasts
    real*8 lasts,lastscx,avescx,newscx,oldscx

    !real*8 scx_s
    !external scx_s

    if (.not.pinavail) then
       scxupdt = 0.0d0
       return

    endif
    if (s.eq.soffset.or.s.lt.lasts) then
       lasts = soffset
       lastscx = 0.0d0
       scxupdt = 0.0d0
       return
    elseif (s.eq.lasts) then
       scxupdt = lastscx
       return
    endif
    newscx = scx_s(s,n,Ti)
    oldscx = scx_s(lasts,nold,Tiold)
    avescx = 0.5d0 * (newscx + oldscx)

    !      if (s.ge.39.0) then
    !      write(6,*) '-----'
    !      write(6,*) 's:new,old',s,lasts
    !      write(6,*) 'Ti:new,old',Ti,Tiold
    !      write(6,*) 'n:new,old',n,nold
    !      write(6,*) 'scx:new,old,ave',newscx,oldscx,avescx
    !      write(6,*) 'scxupdt,lastscx:',scxupdt,lastscx
    !      end if

    scxupdt = lastscx + rcxmom * avescx * (s-lasts)
    lastscx = scxupdt

    lasts = s
    return


  end function scxupdt


  real*8 function estscx(s,n,Ti)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     WF'96: This returns the charge exchange momentum loss -
    !     approximately integrated to the current point s.

    !     This function returns an ESTIMATE of the Scx contribution at
    !     the given S-position ... the lastpei and lasts values are
    !     updated after every R-K iteration.

    implicit none

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'
    real*8 s,n,Ti

    common /scx/ lastscx,lasts
    real*8 lasts,lastscx,scx_tmp

    !real*8 scx_s
    !external scx_s

    if (.not.pinavail) then
       estscx = 0.0d0
       return

    endif
    if (s.eq.soffset.or.s.lt.lasts) then
       estscx = 0.0d0
       return
    elseif (s.eq.lasts) then
       estscx = lastscx
       return

    endif
    scx_tmp = scx_s(s,n,Ti)

    !      if (s.ge.40.0) then
    !      write(6,*) '----'
    !      write(6,*) 's-ls,scxtmp:',(s-lasts),scx_tmp
    !      write(6,*) 'n,Ti',n,Ti
    !      write(6,*) 'est,last',estscx,lastscx
    !      end if

    estscx = lastscx + rcxmom * scx_tmp * (s-lasts)
    return

  end function estscx




  ! external power term - electron  
  
!  real*8 function epowupdt(s)
!    use mod_solparams
!    use mod_solswitch
!    use mod_solcommon
!    implicit none
!    !
!    !     EPOWUPDT: This returns the integrated value of the
!    !               radiation losses up to the point s.
!    !
!    !
!    !     EPOW only involves integrated external terms and as a result
!    !     does not require an update for evolving plasma conditions during a
!    !     solver step. This code structure is maintained here to provide
!    !     future functionality if needed but at the moment epowupdt and estepow
!    !     do the same thing. 
!    !
!    real*8 s
!
!    !common /epowdata/ lastepow,lasts
!    !real*8 lastepow,lasts
!
!    !real*8 teav,nav
!
!    integer in
!
!    if (actswepow.eq.0.0) then
!       epowupdt = 0.0
!    else
!
!       call binsearch(s,in)
!
!       if (in.eq.1) then
!          epowupdt = (intepow(in) * s / sptscopy(1))
!       else
!          epowupdt = ( ( intepow(in-1) + ( (intepow(in)-intepow(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
!       endif
!    endif
!
!    return
!  end function epowupdt



  real*8 function estepow(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     ESTEPOW: This returns the external electron energy term
    !     approximately integrated to the current point s.
    !

    !
    !     EPOW only involves integrated external terms and as a result
    !     does not require an update for evolving plasma conditions during a
    !     solver step. This code structure is maintained here to provide
    !     future functionality if needed but at the moment epowupdt and estepow
    !     do the same thing. 
   
    !                Since everything is analytic at this point and does not
    !                depend on the changing plasma conditions for any of
    !                it's options - it does not have an update method as the
    !                other options may require.

    !
    real*8 s

    !common /praddata/ lastprad,lasts
    !real*8 lastprad,lasts

    !real*8 teav

    integer in

    if (actswepow.eq.0.0) then
       estepow = 0.0
    else
       call binsearch(s,in)
       if (in.eq.1) then
          estepow = (intepow(in) * s / sptscopy(1))
       else
          estepow = ( ( intepow(in-1) + ( (intepow(in)-intepow(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif
    endif
    return
  end function estepow
  

  ! external power term - ion
  
!  real*8 function ipowupdt(s)
!    use mod_solparams
!    use mod_solswitch
!    use mod_solcommon
!    implicit none
!    !
!    !     IPOWUPDT: This returns the integrated value of the
!    !               external ion power losses up to the point s.
!    !
!    !
!    !     IPOW only involves integrated external terms and as a result
!    !     does not require an update for evolving plasma conditions during a
!    !     solver step. This code structure is maintained here to provide
!    !     future functionality if needed but at the moment epowupdt and estepow
!    !     do the same thing. 
!    !
!    real*8 s
!
!    !common /epowdata/ lastepow,lasts
!    !real*8 lastepow,lasts
!
!    !real*8 teav,nav
!
!    integer in
!
!    if (actswipow.eq.0.0) then
!       ipowupdt = 0.0
!    else
!
!       call binsearch(s,in)
!
!       if (in.eq.1) then
!          ipowupdt = (intipow(in) * s / sptscopy(1))
!       else
!          ipowupdt = ( ( intipow(in-1) + ( (intipow(in)-intipow(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
!       endif
!    endif
!
!    return
!  end function ipowupdt



  real*8 function estipow(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     ESTIPOW: This returns the external ion energy term
    !     approximately integrated to the current point s.
    !

    !
    !     IPOW only involves integrated external terms and as a result
    !     does not require an update for evolving plasma conditions during a
    !     solver step. This code structure is maintained here to provide
    !     future functionality if needed but at the moment epowupdt and estepow
    !     do the same thing. 
    !
    !                Since everything is analytic at this point and does not
    !                depend on the changing plasma conditions for any of
    !                it's options - it does not have an update method as the
    !                other options may require.
    real*8 s

    !common /praddata/ lastprad,lasts
    !real*8 lastprad,lasts

    !real*8 teav

    integer in

    if (actswipow.eq.0.0) then
       estipow = 0.0
    else
       call binsearch(s,in)
       if (in.eq.1) then
          estipow = (intipow(in) * s / sptscopy(1))
       else
          estipow = ( ( intipow(in-1) + ( (intipow(in)-intipow(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
       endif
    endif
    return
  end function estipow
  

  
  real*8 function pinqid(s,opt)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    !use mod_params
    use mod_cadas
    implicit none
    real*8 s

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    !     include 'params'
    !     include 'cadas'

    !     PINQID: This function/subroutine implements the
    !             DIVIMP calculated version of PINQI.
    !             The function serves three purposes - selected
    !             by the "opt" switch.

    !             opt = -1 ... Initialize the source calculations
    !                          Pre-integrate any values required
    !                          Setup quantities for later use.

    !             opt = 0  ... Provide an estimate of the integrated
    !                          value of PINQID to the current S
    !                          position.

    !             opt = 1  ... Update the value of the PINQID integral
    !                          Depending on the sub-options
    !                          chosen for the calculation of PINQID
    !                          This function may not be required.


    !     Local variables

    integer opt
    integer in,ik,ir
    real*8 qidatiz(mxspts),qidmliz(mxspts),qidcx(mxspts),qidrec(mxspts)

    real*8 eav,sigvcx,qidsum


    !real*8 rcxmult
    ! external rcxmult

    !     If option is turned OFF - return zero.


    if (actswpcx.ne.4.0) then
       pinqid = 0.0
       return

       !     IF PIN is not available - then also return 0.0 - since
       !     PINQID can not be calculated.

    endif
    if (.not.pinavail) then
       pinqid = 0.0
       return

       !     INITIALIZATION

       !     Setup and initialization - the value passed back is the
       !     net integral of all of the PINQI components.

    endif

    !        Initialize all arrays to zero - this is equivalent to
    !        option zero in all cases - which turns the contribution
    !        to QID due to that effect - OFF.

    if (opt.eq.-1) then
       call qzero(qidatiz,mxspts)
       call qzero(qidmliz,mxspts)
       call qzero(qidrec,mxspts)

       call qzero(qidcx,mxspts)
       DO IK = 1, nptscopy
          PTESA(IK) = oldte(ik)
          PNESA(IK) = oldne(ik) * zb
          PNBS(IK) = oldne(ik)
          PNHS(IK) = nhs(ik)

          !         write (6,*) 'Values OVERWRITTEN due to Tcutoff:'

          !        Atomic ionization

       ENDDO

       if (actswqidatiz.eq.1) then
          do ik = 1,nptscopy
             if (oldti(ik).lt.tcutatiz) then
                !                  write (6,100) 'AT:',ik,oldti(ik),tcutatiz
                qidatiz(ik) = 0.0
             elseif (nhs(ik).gt.0.0) then
                qidatiz(ik) = 1.5 * ths(ik) * econv *nhs(ik)/(nhs(ik)+nh2s(ik)) * ionsrc(ik)
             endif

          end do

       elseif (actswqidatiz.eq.2) then
          write(year,'(i2.2)') iyearh

          call xxuid(useridh)
          ICLASS = 2

          CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
          DO IK = 1, nptscopy
             if (oldti(ik).lt.tcutatiz) then
                !                  write (6,100) 'AT:',ik,oldti(ik),tcutatiz
                qidatiz(ik) = 0.0
             else
                QIdATIZ(IK) = 1.5 * ths(ik) * econv* oldne(ik) * zb * nhs(ik) * PCOEF(IK,1)
             endif

          ENDDO

          !        Molecular Ionization

       endif

       if (actswqidmliz.eq.1) then
          do ik = 1,nptscopy
             if (oldti(ik).lt.tcutmliz) then
                !                  write (6,100) 'ML:',ik,oldti(ik),tcutmliz
                qidmliz(ik) = 0.0
             elseif (nh2s(ik).gt.0.0) then
                qidmliz(ik) = 3.0 * econv *nh2s(ik)/(nhs(ik)+nh2s(ik)) * ionsrc(ik)
             endif

          end do

       elseif (actswqidmliz.eq.2) then
          write(year,'(i2.2)') iyearh

          call xxuid(useridh)
          ICLASS = 2

          CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
          DO IK = 1, nptscopy
             if (oldti(ik).lt.tcutmliz) then
                !                  write (6,100) 'ML:',ik,oldti(ik),tcutmliz
                qidmliz(ik) = 0.0
             else
                QIdmlIZ(IK) = 3.0 * econv* oldne(ik) * zb * nh2s(ik) * PCOEF(IK,1)
             endif

          ENDDO


          !        Recombination

       endif

       if (actswqidrec.eq.1) then
          do ik = 1,nptscopy
             if (oldti(ik).lt.tcutrec) then
                !                  write (6,100) 'RC:',ik,oldti(ik),tcutrec
                qidrec(ik) = 0.0
             elseif (nh2s(ik).gt.0.0) then
                qidrec(ik) = -1.5 * oldti(ik) * econv *recsrc(ik)
             endif

          end do

          !        Charge Exchange

       endif

       !           Prescription from TN1396

       if (actswqidcx.eq.1) then

          do ik = 1,nptscopy
             if (oldti(ik).lt.tcutcx) then
                !                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
                qidcx(ik) = 0.0
             else
                qidcx(ik) = 1.5 * trefcx*rcxmult(oldti(ik))*ceicf* ionsrc(ik)

             endif
          end do

       elseif (actswqidcx.eq.2) then

          do ik = 1,nptscopy

             eav = 1.5 * (ths(ik)+oldti(ik)) / 2.0
             if (eav.gt.1000.0) then
                sigvcx = 1.0d-13
             else
                sigvcx = 1.0d-14 * eav**(1.0d0/3.0d0)

             endif
             if (oldti(ik).lt.tcutcx) then
                !                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
                qidcx(ik) = 0.0
             else
                QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))* econv * oldne(ik) * zb * nhs(ik) * sigvcx

             endif

             if (m0.eq.initm0) then

                if (ik.eq.1) write(6,*) 'QID Source terms:'
                write (6,'(i4,8(1x,g11.4))') ik,ths(ik),oldti(ik),oldne(ik),nhs(ik),nh2s(ik),sigvcx,qidcx(ik),ionsrc(ik)

             endif

          end do

       elseif (actswqidcx.eq.3) then
          write(year,'(i2.2)') iyearh

          call xxuid(useridh)
          ICLASS = 3
          CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
          DO IK = 1, nptscopy
             if (oldti(ik).lt.tcutcx) then
                !                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
                qidcx(ik) = 0.0
             else
                QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))* econv * oldne(ik) * zb * nhs(ik) * PCOEF(IK,1)
             endif

          ENDDO

       elseif (actswqidcx.eq.4) then

          do ik = 1,nptscopy

             eav = 1.5 * (ths(ik)+oldti(ik)) / 2.0
             if (eav.gt.1000.0) then
                sigvcx = 1.0d-13
             else
                sigvcx = 1.0d-14 * eav**(1.0d0/3.0d0)

             endif
             if (oldti(ik).lt.tcutcx) then
                !                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
                qidcx(ik) = 0.0
             else
                QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))* econv * oldne(ik) * zb * nhs(ik) * sigvcx

             endif

             if (qidcx(ik).gt.0.0) qidcx(ik) = 0.0

             if (m0.eq.initm0) then

                if (ik.eq.1) write(6,*) 'QID Source terms:'
                write (6,'(i4,8(1x,g11.4))') ik,ths(ik),oldti(ik),oldne(ik),nhs(ik),nh2s(ik),sigvcx,qidcx(ik),ionsrc(ik)

             endif
          end do

          !        Now that all four components have been calculated - sum them
          !        together and then integrate the total over S.

       endif
       do ik = 1,nptscopy
          qid(ik) = -1.0 *(qidatiz(ik) + qidmliz(ik) + qidrec(ik) + qidcx(ik))

          !        Perform pre-integration into array intqid

          !        This should give the equivalent of the PINQI array calculated
          !        locally - the other two options to this routine access this
          !        integrated array and return a linearly interpolated estimate
          !        of the integral as a function of S.

       end do

       call preint(startn,nptscopy,sptscopy,qid,intqid,qidsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       !        Print out a summary of results


       pinqid = qidsum
       if (m0.eq.initm0) then

          write(6,'(a,g13.6,i4)') 'Sol option 22: QIDsrcint :',qidsum,ringnum
          do ik = startn,nptscopy
             write(6,'(i3,10(1x,g9.3))') ik,sptscopy(ik),qidatiz(ik),qidmliz(ik),qidrec(ik),qidcx(ik),-qid(ik),&
                  -intqid(ik),intqid(ik),qisrc(ik),qesrc(ik)
          end do



          !     ESTIMATE or UPDATE  (at a later time the code for ESTIMATE and UPDATE
          !                          may differ if these quantities are calculated
          !                          for the current iteration - rather than the
          !                          previous one.)

       end if

       !        Find which cell the particle is in

    elseif (opt.eq.0.or.opt.eq.1) then

       call binsearch(s,in)
       if (in.eq.1) then
          pinqid = (intqid(in) * s / sptscopy(1))
       else
          pinqid = (( intqid(in-1) +( (intqid(in)-intqid(in-1)) *(s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))

       endif

       !     Formatted OUTPUT

    endif

100 format(a5,i4,2(1x,g13.5))
    return



  end function pinqid




  subroutine initval
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use error_handling
    use mod_io_units
    use mod_sol22_utils
    
    !     This subroutine calculates several of the quantities
    !     in the solcommon common block ... which are used elsewhere
    !     in the program.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    real*8 srcsum,momsum,gtmp,pinttmp
    real*8 recsum
    
    !real*8 srci,gamma,pintupdt
    !real*8 srcrec,srcf,gperpf
    !external srci,gamma,pintupdt,pinqid,srcrec,srcf,gperpf

    real*8 areasrc(mxspts),asum

    !     Temporary local variables

    integer i,ik

    real*8 tmp1,tmp2,tmp3,tmp4

    !     Calculate velocity at target - including specifed mach number

    gperpcor = 0.0

    !v0 = - m0 * sqrt ( (te0+ti0)/mb * econv/mconv)
    v0 = -m0 * getcs_sol22_dbl(te0,ti0)

    vpe2d = vpe2d * m1/e2dm0
    v1e2d = v1e2d * m1/e2dm0

    !     Calculate pressure at target and include mach number

    vpg   = vpg   * m1/e2dm0
    if (actswe2d.eq.0.0) then
       pstatic0 = n0 * (te0+ti0) * econv
       pinf0 = n0 * (te0+ti0) * (1+m0**2) * econv
       pinf = pinf0
       lastvel = v0
    elseif (actswe2d.eq.1.0) then
       pstatic0 = n1 * (te1+ti1) * econv
       pinf0 = n1 * (te1+ti1) * econv + n1 * vpe2d**2 * mb * mconv
       pinf = pinf0
       lastvel = vpe2d
    elseif (actswe2d.eq.2.0) then
       pstatic0 = n1 * (te1+ti1) * econv
       pinf0 = n1 * (te1+ti1) * econv + n1 * v1e2d**2 * mb * mconv
       pinf = pinf0
       lastvel = v1e2d
    elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.actswe2d.eq.9.0) then
       pstatic0 = n1 * (te1+ti1) * econv
       pinf0 = n1 * (te1+ti1) * econv + n1 * vpg**2 * mb * mconv
       pinf = pinf0
       lastvel = vpg

       !     Additional pressure used if Velocity correction option 2 is in
       !     use.

    endif

    padd = 0.0

    write (6,'(a,i4,6(1x,g20.12))') 'Pressure:',ringnum,pstatic0,pinf0,pinf

    !     Specify the R-value of the target -> r0init

    pinttmp = pintupdt(0.0d0,n0,te0,ti0)


    !     Calculate powers onto the target due to electrons and ions
    !     modified for kinetic transfer and corrective terms.

    r0init = rbnd(0)

    !      if (actswe2d.eq.0.0) then

    gammae = 5.0 + gamecor
    gammai = 2.5 +  0.5* m0**2 * (1.0 + te0/ti0) + gamcor

    !      else
    !         gammai = 2.5 +  0.5* m0**2 * (1.0 + te1/ti1) + gamcor
    !      endif

    !      if (actswe2d.eq.0.0) then

    !         pae = gammae * te0 * econv * n0 * abs(v0)
    !         pai = gammai * ti0 * econv * n0 * abs(v0)

    ! **************************************************************

    !      WARNING !!!!

    ! **************************************************************

    !      There is a possibility for the calculation of pae to
    !      become inconsistent if the calculation of pae_start and
    !      pai_start in the calcfluxes routine performs the
    !      calculations with  a mach number that is NOT equal
    !      to the value used for INITM0.

    !      This is not a problem at present - but if options
    !      are created that assign initial target velocities that
    !      are not equal to the sound speed then this could become
    !      a problem.

    !      D. Elder    March 6, 1998

    if (actswe2d.eq.1.0) then
       pae = gammae * te1 * econv * n1 * abs(vpe2d)
       pai = gammai * ti1 * econv * n1 * abs(vpe2d)
    elseif (actswe2d.eq.2.0) then
       pae = gammae * te1 * econv * n1 * abs(v1e2d)
       pai = gammai * ti1 * econv * n1 * abs(v1e2d)
    elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0) then
       pae = gammae * te1 * econv * n1 * abs(vpg)
       pai = gammai * ti1 * econv * n1 * abs(vpg)
    elseif (actswe2d.eq.9) then
       if (m0.eq.initm0) then
          pae = pae_start
          pai = pai_start
          tmp1 = gammae * te1 * econv * n1 * abs(vpg)
          tmp2 = gammai * ti1 * econv * n1 * abs(vpg)
       else
          pae = gammae * te1 * econv * n1 * abs(vpg)
          pai = gammai * ti1 * econv * n1 * abs(vpg)
          tmp1 = pae
          tmp2 = pai
       endif

       !write (6,'(a,2i4,6(1x,g13.6))')'Target PAE:', ringnum,ike2d_start,pae,tmp1,gammae,gammae*pae/tmp1
       !write (6,'(a,2i4,6(1x,g13.6))')'Target PAI:', ringnum,ike2d_start,pai,tmp2,gammai,gammai*pai/tmp2

    else
       if (m0.eq.initm0) then
          pae = pae_start
          pai = pai_start
          tmp1 = gammae * te0 * econv * n0 * abs(v0)
          tmp2 = gammai * ti0 * econv * n0 * abs(v0)
       else
          pae = gammae * te0 * econv * n0 * abs(v0)
          pai = gammai * ti0 * econv * n0 * abs(v0)
          tmp1 = pae
          tmp2 = pai
       endif

       ! It is possible that pae_start and pai_start are set to zero depending on where this code is being called from
       ! In these cases pae and pai should be set to tmp1 and tmp2 respectively. 
       if (pae.eq.0.0) then
          pae = tmp1
       endif

       if (pai.eq.0.0) then
          pai = tmp2
       endif

       if (pae.eq.0.0.or.pai.eq.0.0) then
          call errmsg('MOD_SOL22_SOURCES: INITVAL:','ERROR:PAE OR PAI ARE ZERO:STOPPING')
          stop 'MOD_SOL22_SOURCES:INITVAL:PAE OR PAI ARE ZERO'
       endif

       ! if pae and pai are assigned from pae_start,pai_start ... they should be equivalent to the values of tmp1,tmp2 unless MACH != 1
       ! but this should be checked and warning issued if not true. 
       if (pae.ne.tmp1.or.pai.ne.tmp2) then
          write(stddbg,'(a,20(1x,g18.8))') 'PAE,PAI WARNING (CHECK VALUES):',pae,tmp1,pai,tmp2
       endif

       !write (6,'(a,2i4,6(1x,g13.6))')'Target PAE:', ringnum,ike2d_start,pae,tmp1,gammae,gammae*pae/tmp1
       !write (6,'(a,2i4,6(1x,g13.6))')'Target PAI:', ringnum,ike2d_start,pai,tmp2,gammai,gammai*pai/tmp2

    endif

    !     Initialize the Area integral if swpow is 3.0

    !write (6,'(a,i4,6(1x,g13.6))')'Target Power:', ringnum,pae,pai

    !        Load up area array

    if (actswpow.eq.3.0) then
       do ik = startn,nptscopy
          areasrc(ik) = 1.0

       end do

       call preint(startn,nptscopy,sptscopy,areasrc,intarea,asum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       intarea(nptscopy+1) = asum
       pinpute = pae * r0init / asum
       pinputi = pai * r0init / asum

       !         write(6,*) 'Ring: ',ringnum,asum

       !         do ik = startn,nptscopy
       !            write(6,*) 'IntArea:',ik,intarea(ik)
       !         end do

    endif

    !     Calculated the target flux

    if (actswe2d.eq.1.0) then
       gamma0 = n1 * vpe2d * targfact
    elseif (actswe2d.eq.2.0) then
       gamma0 = n1 * v1e2d * targfact
    elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.actswe2d.eq.9.0) then
       gamma0 = n1 * vpg * targfact
    else
       gamma0 = n0 * v0 * targfact
    endif


    !     Initial conditions at target

    !     Recombination source setup - pre-calculate the integral.

    if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
       write (6,'(a,12g12.4)') 'Init Targ (Bound) Conditions:',te0,ti0,n0,v0,pinf0,initm0,n0*v0
    endif
    
    recsum = 0.0

    if ((pinavail.and.(actswrecom.eq.1.0)).or.actswrecom.eq.2) then
       !        Integrate over Recombination source term



       call preint(startn,nptscopy,sptscopy,recsrc,intrecsrc,recsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       if (m0.eq.initm0) then

          write(6,'(a,g13.6,i4)') 'Sol option 22: recsrcint :',recsum,ringnum
          do ik = startn,nptscopy
             write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),recsrc(ik),intrecsrc(ik)

          end do

       endif
    endif


    !     Calculate base value for ioniation source function

    !     Momentum source setup - pre-calculate the integral.

    call initioniz

    if (pinavail.and.(actswnmom.eq.6.or.actswnmom.eq.7.or.actswnmom.eq.8)) then

       !        Integrate over Momentum source term


       call preint(startn,nptscopy,sptscopy,momsrc,intmomsrc,momsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)
       if (m0.eq.initm0) then
          write(6,*) 'Sol option 22: momsrcint :',momsum
          do ik = startn,nptscopy
             write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),momsrc(ik),intmomsrc(ik)

          end do

       endif

       !     Calculate base value for momentum loss function.

    endif
    if (actswnmom.eq.0.0 ) then
       smom0 = 0.0
    elseif (actswnmom.eq.1.0 ) then
       smom0 = (pinf/(actlenmom * ringlen))*(1.0d0/actffric-1.0d0)
    elseif (actswnmom.eq.2.0 ) then
       smom0 = (pinf/(actlammom * ringlen))*(1.0d0/actffric-1.0d0)/ ( 1.0d0 -exp( - (actlenmom * ringlen) / (actlammom * ringlen)))
    elseif (actswnmom.eq.3.0 ) then
       !        jdemod - removed length division ... it is a bug
       !         smom0 = (pinf/(actlenmom * ringlen))*(1.0d0/actffric-1.0d0)
       !     >             * (1.0/srci(halfringlen))
       smom0 = pinf * (1.0d0/actffric-1.0d0)* (1.0/srci(halfringlen))
    elseif (actswnmom.eq.4.0 ) then
       smom0 = 0.0
    elseif (actswnmom.eq.5.0 ) then
       smom0 = 0.0
    endif


    !     Calculate base value for radiation function for PRAD option 1.
    if (actswprad.eq.6) then 
       prad0 = frr * (pae + pai) / (lenr-lamr)
    else
       prad0 = frr * (pae + pai) / (lamr * (1.0-exp(-lenr/lamr)))
    endif

    if (actswprad.eq.1.0.or.actswprad.eq.6) then 
       totprad = frr * (pae + pai)
    else
       totprad = 0.0
    endif

    if (debug_s22) write(6,'(a,10(1x,g12.5))') 'PRAD0:',frr,pae,pai,lamr,lenr


    pradsum = 0.0

    call qzero(intrad,mxspts)

    !     Caculate integrated radiation loss for Prad option 4
    !        Integrate over Radiation source term

    if (actswprad.eq.4.0) then 

       !        Print out radiation term for only first iteration on ring

       call preint(startn,nptscopy,sptscopy,radsrc,intrad,pradsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)
       if (m0.eq.initm0) then

          write(6,'(a,1x,g13.6,i4)')'Sol option 22: radsrcint :',pradsum,ringnum
          do ik = startn,nptscopy
             write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),radsrc(ik),intrad(ik)
          end do

       endif



    endif


    ! Add initialization for external electron and ion power terms

    if (actswepow.ne.0.0) then 

       !        Pre-integrate the numerical term to save time during execution

       call preint(startn,nptscopy,sptscopy,epowsrc,intepow,epowsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       !if (m0.eq.initm0) then
       if (debug_s22) then 
           write(6,'(a,1x,g13.6,i4)')'Sol option 22: epowsrcint :',epowsum,ringnum
           do ik = startn,nptscopy
              write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),epowsrc(ik),intepow(ik)
           end do
        endif
    endif

    if (actswipow.ne.0.0) then 

       !        Pre-integrate the numerical term to save time during execution

       call preint(startn,nptscopy,sptscopy,ipowsrc,intipow,ipowsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       !if (m0.eq.initm0) then
       if (debug_s22) then 
           write(6,'(a,1x,g13.6,i4)')'Sol option 22: ipowsrcint :',ipowsum,ringnum
           do ik = startn,nptscopy
              write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),ipowsrc(ik),intipow(ik)
           end do
       endif
    endif
     
    !     PIN Power Source Term SETUP - pre-calculate the integrals.
    call qzero(intqi,mxspts)

    !        Modify QISRC to remove any heating component

    if (pinavail.and.(actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.actswpcx.eq.5)) then
       if (actswpcx.eq.5.0) then
          do ik = startn,nptscopy
             if (qisrc(ik).lt.0.0) qisrc(ik) = 0.0
          end do

       endif


       !        Integrate over Ion energy source term

       call preint(startn,nptscopy,sptscopy,qisrc,intqi,qisum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)
       if (m0.eq.initm0) then

          write(6,'(a,1x,g13.6,i4)')'Sol option 22: qisrcint :',qisum,ringnum
          do ik = startn,nptscopy
             write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),qisrc(ik),intqi(ik)

          end do

       endif

    elseif (actswpcx.eq.4.0) then

       qisum = pinqid(0.0d0,-1)
    else
       qisum = 0.0

    endif


    !     Electron Source term.

    call qzero(intqe,mxspts)

    !        Integrate over Electron energy source term

    if (pinavail.and.(actswphelp.eq.2.0.or.actswphelp.eq.3.0)) then


       call preint(startn,nptscopy,sptscopy,qesrc,intqe,qesum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       if (m0.eq.initm0) then

          write(6,'(a,1x,g13.6,i4)')'Sol option 22: qesrcint :',qesum,ringnum
          do ik = startn,nptscopy
             write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),qesrc(ik),intqe(ik)

          end do

       endif
    else
       qesum = 0.0

    endif

    if (m0.eq.initm0) then

       !        Print out the NETFlux (Gamma) function.

       if (debug_s22) then 
          write (6,'(a,2(1x,g14.6))') 'Qesum, Qisum:',qesum,qisum
          write(6,'(a,g13.6,i4)') 'Sol option 22: GAMMA=nv :',gamma0,ringnum

          write(6,*) '  IK      S        GAMMA-ACT    GAMMA-CALC '//'   GAMMA0   SRCF-GPERPF     GPERPF   '//'     REC          SRCF'
          do ik = startn,nptscopy
             gtmp = gamma(sptscopy(ik))
             tmp1 = srcf(sptscopy(ik))
             tmp2 = gperpf(sptscopy(ik))
             tmp3 = srcrec(sptscopy(ik))
             tmp4 = gamma0
             write(6,'(i4,8(1x,g12.5))') ik,sptscopy(ik),gtmp, tmp4 + tmp1 - tmp3,tmp4,tmp1-tmp2,tmp2,tmp3,tmp1

          end do

       endif
       ! slmod begin - new


    endif
    return



  end subroutine initval


  subroutine initioniz
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     INITIONIZ: This subroutine initializes the
    !                ionization source depending on
    !                specified options and switches.

    !                Among other things - all sources
    !                that involve major radius
    !                correction option 4 are pre-integrated
    !                since it would be too costly to
    !                recalculate the integrals at
    !                each step. These are then linearly
    !                interpolated in the srci function.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    real*8 srcsum,momsum,gtmp,pinttmp
    real*8 tmp1,tmp2,tmp3,tmp4 

    !real*8 srci,gamma,pintupdt
    !real*8 srcf,gperpf,srcrec
    !external srci,gamma,pintupdt,srcf,gperpf,srcrec

    integer i,ik

    !real*8 majrpos
    !external majrpos

    real*8 tmpgamma,tmpsrci,tmprpos,tmpint(mxspts)

    !     IPP/08 Krieger - gtmp and srcsum should be initialized to avoid
    !     run time errors because they are included in write statements but
    !     not used in certain cases

    real*8 tmpsrcsum,ends,eas1,eas2

    gtmp=0.0
    srcsum=0.0

    !     Set the zero offset for Edge2D compatibility cases
    if (actswe2d.ne.0.0) then
       soffset = sptscopy(ike2d_start)
    else
       soffset = 0.0d0
    endif

    !     Check the ionization source position against the offset.

    if (ssrcst.lt.soffset) then
       ssrcst = soffset
       if (ssrcfi.lt.ssrcst) then
          write (7,*) 'ERROR ERROR ERROR - Invalid'//' Ionization Source Specification'
          ssrcfi = ssrcst + 1.0

       endif
       ssrcmid = (ssrcfi+ssrcst) / 2.0

       ssrclen = (ssrcfi-ssrcst)
       !         s5gausslen = s5gausslen - soffset


    endif

    if (actswion.eq.16) then 
       if (ssrcst2.lt.soffset) then
          ssrcst2 = soffset
          if (ssrcfi2.lt.ssrcst2) then
             write (7,*) 'ERROR: ssrcfi2 modified due to soffset (mod_sol22_sources) - Invalid Ionization Source Specification'
             ssrcfi2 = ssrcst2 + 1.0
             
          endif
          !ssrcmid = (ssrcfi+ssrcst) / 2.0
          !ssrclen = (ssrcfi-ssrcst)
          !         s5gausslen = s5gausslen - soffset

       endif
    endif
    
       !     Check and calculate the integral of the PIN ionization source
       !     without ANY cross-field terms - for use later in normalizing
       !     the analytic options.

    
    if (pinnorm.eq.1) then

       call preint(startn,nptscopy,sptscopy,ionsrc,tmpint,tmpsrcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       pnormfact = tmpsrcsum
    else
       pnormfact = -gamma0
    endif


       !     GPERP correction is only done when the
       !     Major radius corrector is off - the major radius
       !     corrector factors in GPERPCOR options.


       !     Calculate the gamma perp correction if the option is turned
       !     ON - otherwise set it to zero. Note that Gamma Perp correction
       !     will only apply to non-normalised ionization - otherwise it will
       !     (by definition) be equal to zero.


    
    if (actswgperp.eq.0.0) then

       gperpcor = 0.0

    elseif (actswgperp.eq.1.0) then

       !        Need to pre-calculate the ionization source
       !        using gperpcor = 0 - then find the difference
       !        and finally recalculate the correct gperpcor and
       !        then re-do the source pre-integration. This involves
       !        repeating the exact code from below twice - though
       !        only for this option. See below for the
       !        commented code.

       !        Only needs to be done for un-normalized sources.
       !        i.e. at present (pinavail.and.actswion.eq.2)

       !        Also need to include this calculation
       !        if RECOMBINATION is ON - and a normalized option has been
       !        specified.

       gperpcor = 0.0

       !          Perform source preintegration - include the R(s) factor
       !          if required.
       !          Also include perpendicular loss corrections.

       !          Set flag passed to preint routine.

       if ((pinavail.and.(actswion.eq.2.0.or.actswrecom.eq.1.0.or.actswrecom.eq.2)).or.(.not.pinavail.and.&
            (actswioni.eq.11.or.actswioni.eq.15))) then

          !          Set the normalization factor

          call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

          !          Use the above to calculate the gperpcor for the 1/2 ring

          fnorm = 1.0
          !fnorm2 = 1.0
          !gperpcor = 0.0
          
          gtmp = gamma(halfringlen)
          if (actswe2d.ne.0.0) then
             gperpcor = -gtmp/(halfringlen-soffset)
          else
             gperpcor = -gtmp/halfringlen

          endif

       endif

       !        GNET is calculated in the routine calcsoliz called from
       !        the beginning of CALCSOL_INTERFACE. It uses the fluxes
       !        at each target and the ionization over the whole ring
       !        to calculate a correction that is distributed over the
       !        whole ring.

       !        Again this is only non-zero for non-normalized sources.

       !write (6,'(a,5g18.7)') 'Gperpcor:1a:',gperpcor,gtmp,srcsum,0.5*ringlen,halfringlen

    elseif (actswgperp.eq.2.0) then
       if ((pinavail.or.(.not.pinavail.and.(actswioni.eq.11.or.actswioni.eq.15))).and.(actswion.eq.2.0.or.pinnorm.eq.1)) then
          gperpcor = gnet
          !write (6,'(a,2g18.7)') 'Gperpcor:2a:',gperpcor,gnet
       else
          gperpcor = 0.0

       endif

       !write (6,'(a,2g18.7)') 'Gperpcor:2b:',gperpcor,gnet

       !        Need to preint the gperp array so as to get the integrated
       !        gperp contribution.

    elseif (actswgperp.eq.3.0) then

       call preint(startn,nptscopy,sptscopy,gperp,intgperp,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       gperpcor = srcsum

       if (m0.eq.initm0) then

          if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then
             write(6,'(a,g13.6,i4)') 'Sol option 22: GPERPsrcint :',srcsum,ringnum
             do ik = startn,nptscopy
                if (ik.lt.100.or.ik.eq.(int(ik/(nptscopy/100)) * int(nptscopy/100))) then 
                   write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),intgperp(ik)
                endif
             end do
          endif
             
       endif

       !        This is the same as option 1 if the 1/2 flux tube is under
       !        ionized and will distribute the X-field flux proportional
       !        to density (as in option 3) if the 1/2 flux tube is
       !        over-ionized.

    elseif (actswgperp.eq.4.0) then

       !        Need to pre-calculate the ionization source
       !        using gperpcor = 0 - then find the difference
       !        and finally recalculate the correct gperpcor and
       !        then re-do the source pre-integration. This involves
       !        repeating the exact code from below twice - though
       !        only for this option. See below for the
       !        commented code.

       !        Only needs to be done for un-normalized sources.
       !        i.e. at present (pinavail.and.actswion.eq.2)

       !        Also need to include this calculation
       !        if RECOMBINATION is ON

       gperpcor = 0.0

       !          Perform source preintegration - include the R(s) factor
       !          if required.
       !          Also include perpendicular loss corrections.

       !          Set flag passed to preint routine.

       if ((pinavail.and.(actswion.eq.2.0.or.actswrecom.eq.1.0.or.actswrecom.eq.2)).or.(.not.pinavail.and.&
            (actswioni.eq.11.or.actswioni.eq.15))) then

          !          Set the normalization factor

          call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

          !          Use the above to calculate the gperpcor for the 1/2 ring

          fnorm = 1.0

          gtmp = gamma(halfringlen)
          if (actswe2d.ne.0.0) then
             gperpcor = -gtmp/(halfringlen-soffset)
          else
             gperpcor = -gtmp/(halfringlen)

             !          Need to re-distribute the flux - for the case where
             !          the flux-tube is over-ionized - i.e. gtmp > 0

          endif

          !             Reset the Gperp option to be the same as 1 - since
          !             it is functionally identical from here on.

          if (gperpcor.ge.0.0) then

             actswgperp = 1.0

             !             Redistribute the integrated flux in proportion to the
             !             density in each cell.

             !             Change back to total required from amount/unit length

          elseif (gperpcor.lt.0.0.and.m0.eq.initm0) then

             gperpcor = gperpcor * (halfringlen)

             !             Copy the density integral and renormalize it so that
             !             the total is equal to the required gperpcor.

             call preint(startn,nptscopy,sptscopy,oldne,tmpint,tmpsrcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)
             do ik = startn,nptscopy
                gperp(ik) = gperpcor/tmpsrcsum * oldne(ik)
                intgperp(ik) = gperpcor/tmpsrcsum * tmpint(ik)
             end do

             if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
                write(6,'(a,g13.6,i4,f7.3,2g13.6)')'Sol option 22: GPERPsrcint :',tmpsrcsum,ringnum,actswgperp,gperpcor,gtmp
                do ik = startn,nptscopy
                   if (ik.lt.100.or.ik.eq.(int(ik/(nptscopy/100)) * int(nptscopy/100))) then 
                       write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),intgperp(ik)
                   endif
                end do
             endif
          endif
       endif

    elseif (actswgperp.eq.5.0) then

       !        Need to pre-calculate the ionization source
       !        using gperpcor = 0 - then find the difference
       !        and finally recalculate the correct gperpcor and
       !        then re-do the source pre-integration. This involves
       !        repeating the exact code from below twice - though
       !        only for this option. See below for the
       !        commented code.

       !        Only needs to be done for un-normalized sources.
       !        i.e. at present (pinavail.and.actswion.eq.2)

       !        Also need to include this calculation
       !        if RECOMBINATION is ON - and a normalized option has been
       !        specified.

       gperpcor = 0.0

       !          Perform source preintegration - include the R(s) factor
       !          if required.
       !          Also include perpendicular loss corrections.

       !          Set flag passed to preint routine.

       if ((pinavail.and.(actswion.eq.2.0.or.actswrecom.eq.1.0.or.actswrecom.eq.2)).or.(.not.pinavail.and.&
            (actswioni.eq.11.or.actswioni.eq.15))) then

          !          Set the normalization factor

          call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

          !          Use the above to calculate the gperpcor for the 1/2 ring

          fnorm = 1.0

          gtmp = gamma(halfringlen)
          if (actswe2d.ne.0.0) then
             gperpcor = -gtmp/(halfringlen-soffset)
          else
             gperpcor = -gtmp/(halfringlen)
          endif

             !          Now that the total is available - distribute it
             !          between the two sources - uniform + rectangular.

             !          Square


          gperpcor2 = gperpcor * gperpfrac* halfringlen /(sgperpend-sgperpbeg)

          !          Reassign gperpcor

          !          Uniform

          gperpcor = gperpcor * (1.0-gperpfrac)

       endif


    elseif (actswgperp.eq.6.0) then
       !        GNET is calculated in the routine calcsoliz called from
       !        the beginning of CALCSOL_INTERFACE. It uses the fluxes
       !        at each target and the ionization over the whole ring
       !        to calculate a correction that is distributed over the
       !        whole ring.

       !        Again this is only non-zero for non-normalized sources.

       !        Cross-field source is split into two components

       if ((pinavail.or.(.not.pinavail.and.(actswioni.eq.11.or.actswioni.eq.15))).and.(actswion.eq.2.0.or.pinnorm.eq.1)) then
          gperpcor = gnet
       else
          gperpcor = 0.0
       endif


          !        Now that the total is available - distribute it
          !        between the two sources - uniform + rectangular.

          !        Square

       gperpcor2 = gperpcor * gperpfrac* halfringlen /(sgperpend-sgperpbeg)

       !        Reassign gperpcor

       !        Uniform

       gperpcor = gperpcor * (1.0-gperpfrac)

       !write(6,'(a,6(1x,g13.6))') 'gperpcor:',gnet,gperpcor,gperpcor2,ringlen,sgperpend,sgperpbeg

    elseif (actswgperp.eq.7.0.or.actswgperp.eq.8.0) then

       !        Need to preint the gperp array so as to get the integrated
       !        gperp contribution.


       call preint(startn,nptscopy,sptscopy,gperp,intgperp,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       gperpcor = srcsum

       if (m0.eq.initm0) then

          if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then
             write(6,'(a,3g13.6,i4)') 'Sol option 22: GPERPsrcint :',srcsum,gnet*ringlen,srcsum/(gnet*ringlen),ringnum
             do ik = startn,nptscopy
                if (ik.lt.100.or.ik.eq.(int(ik/(nptscopy/100)) * int(nptscopy/100))) then 
                   write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),intgperp(ik)
                endif
             end do
          endif
             
       endif
    endif

       !     Do code for major radius corrected ionization source.


    !        For an analytic source - set up the ionsrc array
    !        with the appropriate values and treat numerically
    !        from this point forward.

    !        For swioni = 11 - this is already done.

    if (actswmajr.eq.4.0) then

       !          Watch OUT for S-offsets when assigning this - check all
       !          options and fix bugs.

       !          Approximate s0 for now - it is not accurate since the
       !          integral over S is not the simple exponential
       !          given when the R-correction is applied. However, the
       !          incorrect numbers will be fixed by the fnorm
       !          factor which requires that the integral over
       !          the numerical source be equal to the target flux
       !          for the initial prescription cases.


       !           if (actswe2d.ne.0.0) then
       !             soffset = sptscopy(1)
       !           else
       !             soffset = 0.0d0
       !           endif

       if (actswion.eq.0.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.0.0.and.&
            (.not.pinavail))) then

          s0 = pnormfact *r0init /(ssrcdecay*(1.0-exp(-(ssrcfi-soffset)/ssrcdecay)))

          do ik = startn,nptscopy
             if (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0*exp(-(sptscopy(ik)-soffset)/ssrcdecay)
             else
                ionsrc(ik) = 0.0

             endif

          end do
       elseif (actswion.eq.3.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.3.0.and.&
            (.not.pinavail))) then

          !          Add in triangular source option


          s0 = pnormfact * r0init / ssrclen**2

          do ik = startn,nptscopy
             if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).lt.ssrcmid) then
                ionsrc(ik) = s0 * (sptscopy(ik) - ssrcst)
             elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * (ssrcfi - sptscopy(ik))

             endif

          end do
       elseif (actswion.eq.4.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.4.0.and.&
            (.not.pinavail))) then

          !          Add in triangular source option


          s0 = pnormfact * r0init / ssrclen

          do ik = startn,nptscopy
             if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0

             endif

          end do
       elseif (actswion.eq.6.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.6.0.and.&
            (.not.pinavail))) then

          !          Add in s5gauss source option

          s5alph = 2.5 / s5gausslen ** 2

          s0 = pnormfact * s5alph**3

          do ik = startn,nptscopy
             if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * sptscopy(ik)**5* exp(-s5alph*sptscopy(ik)**2)

             endif

          end do
       elseif (actswion.eq.9.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.9.0.and.&
            (.not.pinavail))) then

          !          Add in shifted S5*gaussian source option

          s5alph = 2.5 / s5gausslen ** 2

          s0 = pnormfact * s5alph**3

          do ik = startn,nptscopy
             if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
             elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * sptscopy(ik)**5* exp(-s5alph*(sptscopy(ik)+s5gausslen/2.0)**2)

             endif

          end do

       endif

          !        Perform source preintegration - include the R(s) factor
          !        Set flag passed to preint routine. Include Perpendicular
          !        flux correction factor.

       call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       if (m0.eq.initm0) then

          if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
             write(6,'(a,g13.6,i4)') 'Sol option 22: FLUXsrcint :',srcsum,ringnum
             do ik = startn,nptscopy
                if (ik.lt.100.or.ik.eq.(int(ik/(nptscopy/100)) * int(nptscopy/100))) then 
                   write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),ionsrc(ik),intionsrc(ik),gperpf(sptscopy(ik))
                endif
             end do
          endif
          
          !        Set the normalization factor

       endif

       if (pinavail.and.actswion.eq.2.0) then
          fnorm = 1.0
       elseif (pinnorm.eq.1) then
          fnorm = 1.0
       else
          fnorm = -gamma0 * r0init / srcsum
       endif

       
          !        Calculate the pre-integral of JUST the ionization
          !        by itself - for use in various source terms.
       call preint(startn,nptscopy,sptscopy,ionsrc,intioniz,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

       if (m0.eq.initm0) then

          if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
             write(6,'(a,g13.6,i4)') 'Sol option 22: IONsrcint :',srcsum,ringnum
             do ik = startn,nptscopy
                if (ik.lt.100.or.ik.eq.(int(ik/(nptscopy/100)) * int(nptscopy/100))) then 
                   write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),ionsrc(ik),intioniz(ik)
                endif
             end do
          endif
             
       endif

       if (pinavail.and.actswion.eq.2.0) then
          fnorm2 = 1.0
       elseif (pinnorm.eq.1) then
          fnorm2 = pnormfact/srcsum
       else
          fnorm2 = -gamma0 * r0init / srcsum
       endif

       !write(6,'(a,l8,20(1x,g12.5))') 'FNORM2A:',pinavail,actswion,pinnorm,fnorm2,fnorm,gamma0,srcsum

    else

          !      Normal options for when the major radius option is not
          !      turned on.

       !        Exponential decay options


       if (actswion.eq.0.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.0.0.and.&
            (.not.pinavail))) then

          s0 = pnormfact / (ssrcdecay*(1.0-exp(-(ssrcfi-soffset)/ssrcdecay)))

          !           Add code for triangular source option.

       elseif (actswion.eq.3.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.3.0.and.&
            (.not.pinavail))) then

          s0 = 4.0 * pnormfact / ssrclen**2

          !           Add code for rectangular source option.

       elseif (actswion.eq.4.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.4.0.and.&
            (.not.pinavail))) then

          s0 = pnormfact / ssrclen

          !           Add code for s5 * gaussian source option.

       elseif (actswion.eq.6.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.6.0.and.&
            (.not.pinavail))) then
          s5alph = 2.5 / s5gausslen **2
          s5alph2 = s5alph**2
          s5alph3 = s5alph**3

          s0 = pnormfact * s5alph3

          !           Add code for s5 * gaussian source option.

       elseif (actswion.eq.9.0.or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0).and.actswioni.eq.9.0.and.&
            (.not.pinavail))) then
          s5alph = 2.5 / s5gausslen **2
          s5offset = s5gausslen/2.0
          s5alph2 = s5alph**2

          s5alph3 = s5alph**3
          ends = ssrcfi + s5offset
          eas2 = exp(-s5alph*ends**2)

          eas1 = exp(-s5alph*s5offset**2)

          s0 = pnormfact /( (-ends**4/(2.0*s5alph)*eas2-ends**2/s5alph2*eas2-eas2/s5alph3)- &
               (-s5offset**4/(2.0*s5alph)*eas1-s5offset**2/s5alph2*eas1- eas1/s5alph3))

          s5startval = -1.0 *(-s5offset**4/(2.0*s5alph)*eas1-s5offset**2/s5alph2*eas1-eas1/s5alph3)

          !write(6,'(a,7g13.5)') 'Init9:',s0,s5offset,s5startval,s5alph,gamma(halfringlen),&
          !     (-ends**4/(2.0*s5alph)*eas2-ends**2/s5alph2*eas2-eas2/s5alph3),&
          !     - (-s5offset**4/(2.0*s5alph)*eas1-s5offset**2/s5alph2*eas1- eas1/s5alph3)

       elseif (actswion.eq.16) then 

          if (dblsrc_opt.eq.0) then ! exp + rect
             s0 = pnormfact / (ssrcdecay*(1.0-exp(-(ssrcfi-soffset)/ssrcdecay)))
             s02 = pnormfact / ssrclen2
          elseif (dblsrc_opt.eq.1) then ! rect+rect
             s0 =  pnormfact / ssrclen
             s02 = pnormfact / ssrclen2
          endif
          
       endif

       !write (6,'(a,9g12.4)') 'InitI:',actswion,actswioni,actswgperp,s0,ssrclen,gperpcor,pnormfact,gamma0
       !write (6,'(a,6g12.4)') 'LensI:',ssrcst,ssrcfi,ssrclen,ssrcmid,soffset,s5gausslen



       if ((pinavail.and.(actswion.eq.1.0.or.actswion.eq.2.0)).or.((actswioni.eq.11.or.actswioni.eq.15).and.&
            (.not.pinavail))) then

       !        Pre-calculate the integral of the source if necessary.

          call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

          !           Set up ionization source normalization if it is required.

          if (actswion.eq.1.0) then
             if (srcsum.ne.0.0) then
                fnorm = -gamma0/srcsum
             else
                fnorm = 1.0
             endif
          elseif (actswion.eq.2.0) then
             fnorm = 1.0
          endif

             !           Print out the source

          if (m0.eq.initm0) then
             if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
                write(6,'(a,g13.6,2i4)') 'Sol option 22: FLUXsrcint :',srcsum,ringnum,nptscopy
                do ik = startn,nptscopy
                   if (ik.lt.100.or.ik.eq.(int(ik/real(real(nptscopy)/100.0)) * int(nptscopy/100))) then 
                      write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),ionsrc(ik),intionsrc(ik),srcf(sptscopy(ik))
                   endif
                end do
             endif
          endif
          
             !           Calculate the pre-integral of JUST the ionization
             !           by itself - for use in various source terms.

          !           Set source normalization

          call preint(startn,nptscopy,sptscopy,ionsrc,intioniz,srcsum,ringlen,actswe2d,actswmajr,sbnd,rbnd,0.0d0)

          if (actswion.eq.1.0) then
             if (srcsum.ne.0.0) then
                fnorm2 = -gamma0/srcsum
             else
                fnorm2 = 1.0
             endif
          elseif (actswion.eq.2) then
             fnorm2 = 1.0
          endif

          !write(6,'(a,l8,20(1x,g12.5))') 'FNORM2B:',pinavail,actswion,pinnorm,fnorm2,fnorm,gamma0,srcsum

             !           Print out the source



          if (m0.eq.initm0) then

             if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then
                write(6,'(a,g13.6,i4)') 'Sol option 22: IONsrcint :',srcsum,ringnum
                do ik = startn,nptscopy
                   if (ik.lt.100.or.ik.eq.(int(ik/real(real(nptscopy)/100.0)) * int(nptscopy/100))) then 
                      write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),ionsrc(ik),intioniz(ik),srci(sptscopy(ik))
                   endif
                end do
             endif
          endif


          !     ENDIF of swmajr

       endif

    endif
    return
  end subroutine initioniz




  real*8 function majrpos(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    !     MAJRPOS: This routine returns the major radius
    !              position corresponding to the S-value
    !              passed to it. The sbnd and rbnd arrays
    !              that are used here are initialized before
    !              the solver is called and are adjusted
    !              from teh DIVIMP values to work away from
    !              each of the targets.

    !              Rbnd and Sbnd give the R and S boundary
    !              positions for each cell ik. The array
    !              starts with a zero index so that the
    !              the boundaries of each cell can be represented
    !              by the values at IK and IK-1. These boundary
    !              values are calculated from the underlying
    !              polygonal grid and thus REQUIRE that such a grid
    !              is in use. Errors will result if this information
    !              is not available. The values for cell boundaries
    !              are calculated in the TAU.D4A module at the
    !              same time as the KSS2 values.

    !              David Elder, Dec 12, 1995


    !     Local variables


    !     Since this routine will be called alot and searching will be
    !     inefficient - it will store the last cell called and check that
    !     first - after that it will perform a binary search on the sbnd array
    !     in order to find the correct cell. The lastik value is reset whenever
    !     S = 0 is encountered. It is initialized the first time to 1 - by the
    !     data statement.

    real*8 s
    integer lastik
    data lastik /1/


    !     Check if S in same cell as last time

    integer in,bot,top,mid

    if (s.ge.sbnd(lastik-1).and.s.le.sbnd(lastik)) then

       !     Not in last cell - Perform Binary search

       majrpos = (rbnd(lastik) - rbnd(lastik-1)) *((s-sbnd(lastik-1))/(sbnd(lastik)-sbnd(lastik-1)))+ rbnd(lastik-1)

       !        BOTTOM is the zeroth element in the SBND array and it should
       !        always be S=0

    else
       in  = 1
       bot = 1

       top = nptscopy

100    continue

       mid = (bot+top)/2
       if (s.le.sbnd(mid)) then
          top = mid
       else
          bot = mid + 1

       endif
       if (bot.eq.top) then
          in = top
          goto 200

       endif

       in = in +1
       if (in.gt.2000) then
          write (6,*) 'SOL22:MAJRPOS: Error in searching for cell: :',in,bot,top,mid,s
          stop 'SOL22:MAJRPOS: Error searching for cell'
       endif

       !        Found cell

       goto 100

200    majrpos = (rbnd(in) - rbnd(in-1)) *((s-sbnd(in-1))/(sbnd(in)-sbnd(in-1)))+ rbnd(in-1)

       !        Debug code to check for errors in search

       lastik = in
       if ( (abs(majrpos-rbnd(in))+abs(majrpos-rbnd(in-1))).ne.abs(rbnd(in)-rbnd(in-1))) then
          write (6,*) 'ERROR: MAJOR RADIUS INCORRECT'
          write (6,'(4i4)') nptscopy,bot,top,mid
          write (6,'(2i4,3g12.5)') in,lastik,majrpos,rbnd(in-1),rbnd(in)
          write (6,'(3g12.4)') s,sbnd(in-1),sbnd(in)
          stop 'SOL22:MAJRPOS: ERROR: MAJOR RADIUS INCORRECT'

       endif

    endif
    return



  end function majrpos

  subroutine binsearch(s,in)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    !     Search for right cell - record last in cell
    !     and check that first since most routines will be
    !     in the same cell more often than not.


    !     Local Variables

    integer in
    integer bot, top, mid, lastin

    !     Need to reset lastin if it is too large from
    !     prevous rings.

    data lastin /1/

    if (lastin.gt.nptscopy) lastin = 1
    if (lastin.eq.1.and.s.lt.sptscopy(lastin)) then
       ! slmod begin - new
       !...ARRAY BOUNDS:
       in = 1

       !      elseif ((lastin.gt.1).and.(s.ge.sptscopy(lastin-1))
       !     >       .and.(s.lt.sptscopy(lastin))) then
       ! slmod end
    elseif ((lastin.gt.1).and.(s.ge.sptscopy(MAX(1,lastin-1))).and.(s.lt.sptscopy(MAX(1,lastin)))) then
       in = lastin
    elseif (s.gt.sptscopy(nptscopy)) then
       in = nptscopy + 1
    else
       in  = 1
       bot = 1

       top = nptscopy

100    continue

       mid = (bot+top)/2
       if (s.lt.sptscopy(mid)) then
          top = mid
       else
          bot = mid + 1

       endif
       if (bot.eq.top) then
          in = top
          goto 200

       endif
       in = in +1
       !               write (6,*) 'Error in search:',in,bot,top,mid,s
       if (in.gt.2000) then
          stop 'SOL22:BINSEARCH: ERROR FINDING CELL'
       endif

       goto 100
200    continue

       lastin = in

    endif
    if (debug_s22) then 
       write(6,'(a,2i4,3(1x,g12.5))') 'BIN:',in,nptscopy,s,sptscopy(in-1),sptscopy(in)

    endif
    return

  end subroutine binsearch


  real*8 function areaint(s)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    real*8 s
    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'
    !
    !     Returns the estimates R(s)*ds  integral at position S.
    !
    integer top, bot,mid, in
    !
    !
    !      if (s.gt.sptscopy(nptscopy)) then
    !
    !         in = nptscopy +1
    !
    !         areaint = ( intarea(in-1) +
    !     >             ( (intarea(in)-intarea(in-1)) *
    !     >             (s-sptscopy(in-1))/((ringlen/2.0)-sptscopy(in-1))))
    !
    !      else
    !
    !      BINSEARCH returns in+1 for S > spts(npts) 
    !
    call binsearch(s,in)

    if (in.eq.1) then
       areaint = intarea(in) * s / sptscopy(1)
    else
       areaint = ( intarea(in-1) + ( (intarea(in)-intarea(in-1)) * (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1))))
    endif
    !
    !      endif
    !

    return
  end function areaint


      SUBROUTINE SOL22Headers
! ======================================================================

! subroutine: SOL22Headers


      !use mod_params
      !use mod_slcom
      use mod_solparams
      use mod_solcommon
!     INCLUDE 'params'
!     INCLUDE 'slcom'
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
      IMPLICIT none
      IF (sol22_osm_mode.LE.1) RETURN
      IF (miter.EQ.1) THEN
        WRITE(75,*)
        WRITE(75,'(A,A3,A8,1X,2A7,A10,1X,A10,A5,2(1X,A10))')'`','in','s','Ti','Te','ne','Vb','M','Ga','P'
        IF (sol22_outmode.GE.3) THEN
          WRITE(71,*)
          WRITE(71,'(A,A3,4(1X,A10,10X))')'`','in','pais','paes','peis','srcf'
        ENDIF
        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
          WRITE(72,*)
          WRITE(72,'(A,A3,1X,A7,3A10,10X,A9,2X,3A10)')'`','in','Ti','Pcf','Pcx','Pei','Pu/Pt','Conv','Cond','Total'
        ENDIF
        WRITE(73,*)
!        WRITE(73,*)
!        WRITE(73,'(A,A3,1X,A7,3A10,A9,2X,3A10)')
!     .    '`','in','Te','Pcf','PHi','Pei','Pu/Pt',
!     .    'Conv','Cond','Total'
        WRITE(73,'(A,A3,1X,A7,4A10,A9,2X,3A10)')'`','in','Te','Pcf','PHi','Pei','Prad','Pu/Pt','Conv','Cond','Total'
      ENDIF
      RETURN
99    STOP



    END SUBROUTINE SOL22Headers




    SUBROUTINE SOL22Output(loopstart,spts,npts,conde,condi,conve,convi,pcxv,peiv,phelpiv,pradv,te,ti,ne,vb,ga,act_press,pmloss,note)
! ======================================================================

! subroutine: SOL22Output

      !use mod_params
      !use mod_comtor
      !use mod_slcom
      use mod_solparams
      use mod_solcommon
      use mod_sol22_utils
!     INCLUDE 'params'
!     INCLUDE 'comtor'
!     INCLUDE 'slcom'
!     INCLUDE 'solparams'
!     INCLUDE 'solcommon'
      IMPLICIT none


      ! jdemod - this common block does not exist in any of the other SOL22 source code modules
      !          in addition - mxspts is no longer a constant due to the shift to dynamic allocation
      !          and arrays with variable size are not allowed in common blocks.
      !          These quantities are also not assigned a value in this routine so I have commented
      !          out dp4 and dp6 and will remove references to them in this routine.
      !
      !COMMON /OUTPUTJUNK/ dp4        ,dp6
      !REAL                dp4(MXSPTS),dp6(MXSPTS)
      !REAL     GetCs_sol22

      !REAL*8   cond,conv,paes,pais,pmomloss,gperpf
      !EXTERNAL cond,conv,paes,pais,pmomloss,gperpf

      INTEGER i,loopstart,npts
      !REAL*8  srcf,powi,powe,mach,te(MXSPTS),spts (MXSPTS),pmloss(MXSPTS),exp_press(MXSPTS),ne(MXSPTS),&
      REAL*8  powi,powe,mach,te(MXSPTS),spts (MXSPTS),pmloss(MXSPTS),exp_press(MXSPTS),ne(MXSPTS),&
           prad (MXSPTS),pcxv  (MXSPTS),act_press(MXSPTS),ga(MXSPTS),peiv (MXSPTS),pradv (MXSPTS),&
           phelpiv  (MXSPTS),ti(MXSPTS),condi(MXSPTS),conde (MXSPTS),vb(MXSPTS),convi(MXSPTS),conve (MXSPTS)
      CHARACTER*2 note(MXSPTS)
      ! jdemod - dumpai1 is printed below, declared in slcom but never assigned any value so I am just replacing
      ! it with a local variable assigned a value of zero - I am commenting them out in slcom as well
      real*8 :: dumpai1,dumpae1,dumpei1
      dumpai1 = 0.0
      dumpae1 = 0.0
      dumpei1 = 0.0
      
      !     Pcx > 0 => cooling
!     PHi > 0 => cooling
!     Pei > 0 => electron cooling and ion heating
      IF (sol22_osm_mode.LE.1) RETURN
      DO i = loopstart, npts
        powi = (pai - pais(spts(i))) - pcxv   (i) + peiv(i)
        powe = (pae - paes(spts(i))) - phelpiv(i) - peiv(i)
        mach = GetCs_sol22_dbl(te(i),ti(i))
!     .      srcf(spts(i)),gperpf(spts(i)),
        WRITE(70,'(1X,I3,F8.3,1X,2F7.2,1P,E10.2,1X,E10.2,0P,F5.2,1P,2(1X,E10.2),1X,2E10.2,0P,F10.4,A)')i,spts(i),ti(i),&
             te(i),ne(i),vb(i),DABS(vb(i)/mach),ga(i),act_press(i)/ECONV,ionsrc(i),gperpf(spts(i)),pmloss(i),note(i)
        IF (forcet.EQ.0.OR.forcet.EQ.2.OR.forcet.EQ.3) THEN
           WRITE(72,'(1X,I3,1X,F7.2,1P,3E10.2,0P,10X,F9.3,2X,1P,3E10.2,0P,A)')i,ti(i),-(pais(spts(i))-pai),-pcxv (i), &
                peiv(i),powi/pai,convi(i),condi(i),convi(i)+condi(i),note(i)
!...Prad!
        ENDIF                            
!        WRITE(73,'(1X,I3,1X,F7.2,1P,3E10.2,0P,F9.3,
!     .             2X,1P,3E10.2,0P,1X,2F8.2,A)')
!     .      i,te(i),-(paes(spts(i))-pae),-phelpiv(i),-peiv(i),powe/pae,
!     .      conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),note(i)
!
        ! jdemod - replace dp4 and dp6 with zeroes in the following output - since they are undefined
        !
        !WRITE(73,'(1X,I3,1X,F7.2,1P,4E10.2,0P,F9.3,'//' 2X,1P,3E10.2,0P,1X,2F8.2,A)')i,te(i),-(paes(spts(i))-pae),&
        !     -phelpiv(i),-peiv(i),-pradv(i),powe/pae,conve(i),conde(i),conve(i)+conde(i),dp4(i),dp6(i),note(i)
        WRITE(73,'(1X,I3,1X,F7.2,1P,4E10.2,0P,F9.3,'//' 2X,1P,3E10.2,0P,1X,2F8.2,A)')i,te(i),-(paes(spts(i))-pae),&
             -phelpiv(i),-peiv(i),-pradv(i),powe/pae,conve(i),conde(i),conve(i)+conde(i),0.0,0.0,note(i)
        IF (sol22_outmode.GE.3)WRITE(71,'(1X,I3,1P,4(1X,2E10.2),3X,E10.2,A)')i,pais(spts(i)),dumpai1,paes(spts(i)),&
             dumpae1,peiv(i),dumpei1,srcf(spts(i)),-1.0,pmomloss(spts(i),1,vb(i),te(i),ti(i)),note(i)
      ENDDO
      RETURN
99    STOP

! ======================================================================


! ======================================================================

    END SUBROUTINE SOL22Output


end module mod_sol22_sources


