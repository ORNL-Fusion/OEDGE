module mod_sol22_solver

  use mod_sol22_sources

  implicit none


contains



  subroutine solvstep(sinit,send,t1e,t1i,n,exitcond,imflag,negerrflag,vcount)
    !use sol22_input
    use sol22_debug
    use mod_solparams
    use mod_solcommon
    use mod_solswitch
    implicit none
    real*8 sinit,send,t1i,t1e,n
    integer exitcond,imflag,vcount,negerrflag
    !
    !     This routine uses an RK driver routine to solve from
    !     sinit to send - if it finds that this can't be done it
    !     returns a non-zero exitcondition. A zero exit condition
    !     means that it was successful in stepping from sinit to send.
    !
    !     Local variables
    !
    integer flag,imflag2
    real*8 s,h,newt1e,newt1i,errte,errti
    integer ierr
    logical errneg,imag

    real*8 fval

    !real*8 newn,pcxupdt,peiupdt,phelpiupdt,gamma,pintupdt
    !real*8 pmomloss,pradupdt,scxupdt
    !external newn,pcxupdt,peiupdt,phelpiupdt,gamma
    !external pintupdt,pmomloss,pradupdt,scxupdt

    !real*8,external :: srci,srcf

    real*8 ptmp
    real*8 tmpgam1,tmpgam2,vsep,vtmp
    real*8 v1,v2,vsub,vsup
    !real*8 smomtmp,scxtmp
    real*8 errmax,nimag,hnew,stmp,tmpgam

    real*8 hneg

    real*8 ntmp,tmpnimag,gtmp
    real*8 srtn
    integer tmpflag,iter

    !real*8 fegrad figrad
    !external fegrad,figrad
    !
    !     Set up the separation factor - if s exceeds this
    !     value * ssepmin - then one assumes that there
    !     is an acceptable solution for the entire ring at
    !     this Mach number..
    !
    real sepfact
    parameter (sepfact= 1.2)

    !
    ! Note: Finding imaginary results is normal for the solver when the step size is too large.
    !       This is particularly common on the very first step when the step size is the largest.
    !       This occurs more frequently due to a code update.
    !       This status is still reported but imaginary numbers occuring for the first solver
    !       step (simag=0.0) are no longer reported as errors in the sol22 driver routine. 
    !
    
    iter = 1
    imflag = 0
    negerrflag = 0
    imflag2 = 0
    exitcond = 0
    s = sinit
    slast = sinit

    nlast = n
    !
    !     Set the value of hlim - the minimum allowable
    !     step size in the equation solver. Based on the
    !     separation of the two points and the initial number
    !     of steps - keeping in mind that the default step
    !     size for default imaginary cases is 10 times
    !     smaller than the initial step-size.
    !
    hlim = (send - sinit) / dble(ndiv) / 10000.0
    hneg = (send - sinit) / dble(ndiv) / 500.0
    himag = (send - sinit) / dble(ndiv) / 50.0
    !
    !
    !     starting h-value
    !
    h = (send-sinit)/ dble(ndiv)

    if (debug_s22)  write (6,*) 'Hlim:',hlim,himag,sinit,send,h


2000 continue

    imflag2 = 0
    !
    !
    !------------------------------------------------------------------------
    !
    !     Try to perform an RK step of size h
    !
    !
    !
    call rkstep(s,t1e,t1i,newt1e,newt1i,errte,errti,h,m0,ierr)
    !
    !     Save new S value returned by rkstep - used to decide 
    !     whether to save debug information. 
    !     
    srtn = s + h
    !
    !write(6,'(a,10(1x,g20.12))') 'RKSTEP:',s,t1e,t1i,newt1i,newt1e,h,m0,ierr
    !
    !------------------------------------------------------------------------
    !
    !     Check error conditions from stepper
    !
    if (ierr.eq.1) then
       !
       !         Negative temperature - adjust h smaller and iterate
       !
       negerrflag = 1
       !
       !          write (6,'(a,5(1x,g12.5))')
       !     >         'Error!: Negative TEMP',h,hlim,himag,s,m0
       !          write (6,*) 'T:',t1e,newt1e,t1i,newt1i
       !
       !         ntmp = newn(s,t1e,t1i,tmpnimag,tmpflag)
       !          write (6,*) 'G:',tmpflag,gamma(s),ntmp,
       !               fegrad (s,t1e,t1i,ntmp),
       !               figrad (s,t1e,t1i,ntmp)
       !
       if (h.le.0.01*hneg) then
          if (lastiter) then
             write (6,*) 'ERROR: Negative temperature: h too small on last iteration'
             write (6,*) 'Iteration Stopping in SOLASCV:SOLVSTEP',iter
             write (6,*) 'Possibly caused by Physically Inconsistent Target Conditions'

             write (6,'(a,5(1x,g13.5))') 'Error!: Negative TEMP',h,hlim,himag,s,m0
             write (6,*) 'T:',t1e,newt1e,t1i,newt1i
             ntmp = newn(s,t1e,t1i,tmpnimag,tmpflag)
             gtmp = gamma(s)
             write (6,'(a,i4,4(1x,g13.5))') 'G:',tmpflag,gtmp,ntmp,fegrad (s,t1e,t1i,ntmp),figrad (s,t1e,t1i,ntmp)

             exitcond = 3
             snegerr = s
             return
          else
             write (6,*) 'ERROR: Negative temperature: h too small on last iteration'
             exitcond = 3
             snegerr = s
             return
          endif
       endif
       h = h / 10.0
       goto 2000

    elseif (ierr.eq.2) then
       !
       !
       !     Check for imaginary solutions.
       !

       simag = s

       !if (simag.eq.0.0) then
       !   write (6,'(a,8(1x, g12.5))') 'imag at 0.0 - reducing h:',h,m0,hlim,s,send,sinit,newt1e,newt1i
       !endif

       if (stopimag.and.(.not.lastiter)) then

          if (h.le.hlim) then
             exitcond = 1
             return
          elseif (h.gt.hlim) then
             h = h / 10.0
             goto 2000
          endif

       else

          if (h.le.hlim.and.(s+h).ne.send) then
             imflag = 1
             imflag2 = 1
             write (6,*) 'ERROR: Imaginary encountered: h too small:',ringnum,h,m0,lastiter,lastiters, &
                  stopimag,ierr,hlim,s,send,sinit
             write (6,*) 'Code will NOT stop'
             !               stop
             !
             exitcond = 4
             return

          elseif (h.gt.himag) then
             h = h /10.0
             goto 2000
          else
             imflag = 1
             imflag2 = 1
          endif
       endif
    elseif (ierr.eq.3) then!
       !
       !     Check for negative N
       !

       !        Negative N found - iterate smaller h OR set condition and exit
       !
       if (h.le.hlim) then
          exitcond = 6
          snegerr  = s
          return
       elseif (h.gt.hlim) then
          h = h / 10.0
          goto 2000
       endif

    elseif (ierr.eq.4) then
       !
       !     Check for NaNQ
       !
       !
       !        NaNQ found - iterate smaller h OR set condition and exit
       !
       if (h.le.hlim) then
          exitcond = 7
          snegerr  = s
          return
       elseif (h.gt.hlim) then
          h = h / 10.0
          goto 2000
       endif

    endif
    !
    !
    !------------------------------------------------------------------------
    !
    !     Need to test if h is too small here because if the
    !     code is not numerically stable then a condition could
    !     conceivably be generated where the h value just gets
    !     smaller and smaller without encountering imaginary
    !     numbers.
    !
    if (h.le.hlim) then

       if (lastiter) then
          exitcond = 4
          snegerr = s
          write (6,*) 'ERROR: delta h less than minimum :',h,hlim
          write (6,*) 'Possibly caused by Negative T: ',ierr
          return
       elseif (actswmach.ne.0.0.and.m0.lt.10.0*origm0) then

          if (s.gt.ssepmin.and.vsep.gt.vsepmin.and.m0.ne.origm0) then
             !
             !              Try iterating at smaller mach number
             !
             exitcond = 2

             return

          else
             !
             !              Try iterating at higher mach number
             !
             exitcond = 1
             snegerr  = s
             !
             !            write (6,*)'H:',h,hlim,snegerr
             !
             return

          endif

       else
          write (6,*) 'ERROR: delta h less than minimum :',h,hlim
          write (6,*) 'Program NOT stopping : ',ierr

          exitcond = 4
          snegerr  = s
          return
          !
          !            stop
          !
       endif

       !
       !         if (s.gt.lastiters) then
       !            exitcond = 4
       !            snegerr = s
       !            write (6,*) 'ERROR: delta h less than minimum :',
       !     >                  h,hlim
       !            write (6,*) 'Possibly caused by Negative T: ',ierr
       !            return
       !         else
       !            write (6,*) 'ERROR: delta h less than minimum :',
       !     >                  h,hlim
       !            write (6,*) 'Program NOT stopping : ',ierr
       !
       !            exitcond = 4
       !            snegerr  = s
       !            return
       !
       !            stop
       !
       !         endif
       !
    endif
    !
    !------------------------------------------------------------------------
    !
    !     Check to see if it was a successful step - if it was go on
    !     the next step using the revised h-step - if it wasn't retry
    !     with an appropriately revised h step.
    !
    errmax = max (abs(errte),abs(errti))

    errmax = errmax / eps
    !
    !      write (6,*) 'errs:',errmax,errte,errti
    !
    if (errmax.gt.1.0d0.and.imflag2.eq.0) then
       !
       !        Error is greater than the desired amount - retry with
       !        smaller h.
       !
       hnew = 0.9 * h * (errmax ** (-0.25))
       if (hnew.lt.0.1*h) then
          hnew = 0.1 * h
       endif
       s = slast

    else

       !
       !        Error is less - can afford to increase h a bit
       !
       stmp = slast + h
       if ((abs(stmp-send).le.hlim)) stmp = send

       s = stmp

       if (imflag2.eq.0) then
          if (errmax.gt.((5.0/0.96)**(-5.0))) then
             hnew = 0.96 * h * (errmax ** (-0.2))
          else
             hnew = 5.0 * h
          endif
       elseif (imflag2.eq.1) then
          hnew = himag
       endif
       !
       !         write (6,*) 'h1:',h,hnew,imflag2,s
       !         write (6,*) 'h2:',send,sinit,ndiv,dble(ndiv)
       !
       if ((s+hnew).ge.send) then
          hnew = send - stmp
       endif

       tmpgam = gamma(s)
       n = newn(s,newt1e,newt1i,nimag,flag)
       !
       !        Deal with NEWN return codes
       !
       if (flag.eq.0) then
          vsub=0.0
          vsup=0.0
          if (founds) then
             if (n.ne.0.0) vsub = tmpgam/n
             if (nimag.ne.0.0) vsup = tmpgam/nimag
          else
             if (nimag.ne.0.0) vsub = tmpgam/nimag
             if (n.ne.0.0) vsup = tmpgam/n
          endif
          vsep = abs(vsub-vsup)

       elseif (flag.eq.2) then
          !
          !           Negative N encountered - errcode/exitcond=6
          !
          exitcond = 6
          snegerr  = s
          return

       endif
       !
       !
       !        write(6,*) 'VSEPA:',stopimag,flag,founds,lastiters,vsep,vsepmin
       !
       if (stopimag.and.flag.eq.0) then
          if (vsep.le.vsepmin.and.(.not.founds)) then
             vsepmin = vsep
             ssepmin = s
             if (vsepmin.lt.5.0) then
                vcount = vcount + 1

                if (vcount.gt.1) then
                   founds = .true.
                   lastiters = s
                endif
             endif

             !               write(6,*) 'VSEP:',vsepmin,ssepmin,vcount,lastiters

             !
             !            elseif (s.gt.(sepfact*ssepmin).and.s.gt.(0.1*ringlen))
             !
          elseif (s.gt.(sepfact*ssepmin).and.ssepmin.ne.0.0.and.s.gt.(0.1*halfringlen)) then
             !
             !             elseif (vsep.gt.sepfact*vsepmin) then
             !
             !              Need to go back to last mach number that didn't
             !              work and start up from there.
             !
             !               if (ringnum.eq.6) then
             !                 write(6,*) 'lastiter:',lastiter,founds,s,ssepmin,
             !                           lastiters,m0
             !               endif
             !
             if (lastiter) then
                stopimag = .false.
             elseif (m0.ne.origm0) then
                exitcond = 2
                !
                !                 If founds is false - still set lastiters to
                !                 ssepmin - just in case the soltion with the
                !                 specified mach number resolution will not
                !                 result in a valid value of lastiters being found.
                !
                if (.not.founds) then
                   lastiters = ssepmin
                endif

                return
             endif

          elseif (halfringlen.lt.(sepfact*ssepmin)) then
             !
             !             In this case S can never be greater than
             !             sepfact * ssepmin and so some other calculation
             !             is needed - in this case  1/4 of the remaining
             !             distance S to the end of the ring from the
             !             point of minimum separation is used.
             !
             if (s.gt.(0.25*(halfringlen-ssepmin)+ssepmin)) then
                !
                !                Need to go back to last mach number that didn't
                !                work and start up from there.
                !
                !                 write(6,*) 'M iter:',lastiter,s,ssepmin,lastiters,m0
                !
                if (lastiter) then
                   stopimag = .false.
                elseif (m0.ne.origm0) then
                   exitcond = 2
                   return
                endif
             endif

          else
             vcount = 0

          endif
       endif
       !
       !     End of errmax if statement
       !
    endif

    !
    !------------------------------------------------------------------------
    !
    if ((s.ne.send).and.(abs(send-s).le.hlim)) then
       s = send
    endif

    if ((s.eq.send).or.(s.ne.slast)) then

       iter = iter + 1
       !
       !        Update the integrals at the end of
       !        each R-K step.
       !
       n = newn(s,newt1e,newt1i,nimag,flag)

       if (flag.eq.2) then
          !
          !           Negative N encountered - errcode/exitcond=6
          !
          exitcond = 6
          snegerr  = s
          return

       endif

       ptmp = peiupdt(s,n,nlast,newt1e,newt1i,t1e,t1i,fval)
       ptmp = pradupdt(s,n,nlast,newt1e,t1e)
       ptmp = phelpiupdt(s,n,nlast,newt1e,t1e)
       ptmp = pcxupdt(s,newt1i,t1i)
       ptmp = pintupdt(s,n,newt1e,newt1i)

       ptmp = pradupdt(s,n,nlast,newt1e,t1e)

       
       if (actswnmom.eq.9.or.actswnmom.eq.10) then
          ptmp = scxupdt(s,n,nlast,newt1i,t1i)
       endif

       tmpgam1 = gamma(s)

       vtmp = 0.0
       if (n.ne.0.0) vtmp = tmpgam1 / n

       lastvel = vtmp

       !if (actswnmom.eq.4) then
       !   smomtmp = pmomloss(s,1,vtmp,newt1e,newt1i)
       !endif

       !
       !        Calculate an estimate of the velocity gradient
       !        VGRAD for possible use in future viscosity options.
       !
       !        if (imflag2.eq.1) then
       !           vgrad = 0.0
       !        else
       tmpgam1 = gamma(s-h)
       tmpgam2 = gamma(s)

       v1=0.0
       if (nlast.ne.0.0) v1 = tmpgam1/nlast

       v2=0.0
       if (n.ne.0.0) v2 = tmpgam2/n

       vgrad = (v2-v1) / h
       !
       !        endif
       !
       nlast  = n
       slast  = s

       t1i = newt1i
       t1e = newt1e

       timax = max(timax,t1i)
       temax = max(temax,t1e)
       !
       !        Exit if error correction is ON and the
       !        temperature drops too much.
       !
       if (actswerror.ne.0.0.and. ((t1i.lt.dropfrac*timax).or. (t1e.lt.dropfrac*temax).or. &
            (t1i.le.tfloor).or. (t1e.le.tfloor))) then

          exitcond = 5
          snegerr = s
          return

       endif
       !
       !         sinit = s
       !
       if (s.eq.send) return
    endif

    if (hnew.le.0.0) then
       write (6,*) 'Hnew = ',hnew,' last H = ',h
       stop
    endif

    h = hnew

    !
    !     If debugging is on - record the data for the current step
    !
    if (debug_sol22_on) then 
       !write(0,*) 'mod_sol22_solver:',debug_sol22_on
       !if (debug_sol22_on.and.srtn.eq.slast) then 
       call save_s22_data(h,s,n,t1e,t1i,lastvel,gamma(s),srci(s),srcf(s),press(s,t1e,t1i))
    endif
    !
    goto 2000
    !
    !     End of Solvstep
    !
  end subroutine solvstep




  subroutine  rkstep(s,t1e,t1i,newt1e,newt1i,errte,errti,h,m0arg,ierr)
    use mod_solparams
    use mod_solcommon
    use mod_solrk
    use mod_solswitch
    implicit none
    !
    real*8 s,t1i,t1e,newt1i,newt1e,errte,errti,h,m0arg
    integer ierr
    !
    !     Local variables
    !
    integer i,j,flag
    real*8 ke(6),ki(6),nimag,n
    real*8 t2ep,t2ip
    real*8 si,ti,te


    !real*8 fegrad,figrad,newn,fgrad
    !external fegrad,figrad,newn,fgrad
    !
    !     Initialization
    !
    ierr = 0
    newt1e = t1e
    newt1i = t1i
    !
    !     Loop through calculating each step for te and ti
    !
    do i = 1,6
       !
       !        set up values at start of RK step
       !
       si = s + ai(i) * h
       ti = t1i
       te = t1e
       !
       !        recalculate Ti,Te for the current step
       !
       do j = 1,i-1

          ti = ti + bij(i,j) * ki(j)
          te = te + bij(i,j) * ke(j)

       end do
       !
       !        Test for valid te and ti values (i.e. > 0 ) - if this is
       !        not the case - exit - issue an error message and reduce
       !        the step-size by a factor of 10.
       !
       if ((te.lt.0.0).or.(ti.lt.0.0)) then
          ierr = 1
          return
       endif
       !
       !        Calculate n for this step
       !
       n = newn (si,te,ti,nimag,flag)
       !
       !        Check for imaginary n's
       !
       if (flag.eq.1) then
          !
          ierr = 2
          !
          if (((actswmach.eq.1.0.or.actswmach.eq.2.0).and.(.not.lastiter))&
               &.or.(lastiter.and.h.gt.himag))  return
          !
       elseif (flag.eq.2) then
          !
          !           Negative N encountered - ierr=3 - errcode/exitcond=6
          !
          ierr = 3
          return
          !
       endif
       !
       !        Check for NaNQ values in Te, Ti and N
       !
       if ((.not.(te.le.0.0.or.te.gt.0.0)).or.(.not.(ti.le.0.0.or.ti.gt.0.0)).or.(.not.(n.le.0.0.or.n.gt.0.0))) then
          write (6,'(a,3(1x,g13.5))') 'NaNQ Error:',te,ti,n
          ierr = 4
          return
       endif
       !
       !        Calculate the k-values
       !
       if (forcet.eq.0) then
          ke(i) = h * fegrad (si,te,ti,n)
          ki(i) = h * figrad (si,te,ti,n)
       elseif (forcet.eq.1) then
          ke(i) = h * fgrad (si,te,ti,n)
          ki(i) = ke(i)
       endif

    end do
    !
    !     If it has been successful - calculate the errors in the Te,Ti
    !     estimates
    !
    newt1e = t1e
    newt1i = t1i
    t2ep = t1e
    t2ip = t1i
    !
    do i = 1,6
       newt1e = newt1e + ci(i) * ke(i)
       newt1i = newt1i + ci(i) * ki(i)
       t2ep = t2ep + cip(i) * ke(i)
       t2ip = t2ip + cip(i) * ki(i)
    end do
    !
    !        Test for valid te and ti EXIT values (i.e. > 0 ) - if this
    !        is not the case - return - issue an error message and reduce
    !        the step-size by a factor of 10 in the calling routine.
    !
    if ((newt1e.lt.0.0).or.(newt1i.lt.0.0).or.(t2ep.lt.0.0).or.(t2ip.lt.0.0)) then
       ierr =1
       !         write(6,*) 'errneg2:',m0,h,s,newt1e,newt1i,
       !     >                            t2ep,t2ip
       return
    endif

    errte = newt1e - t2ep
    errti = newt1i - t2ip
    !
    !     NOTE: Due to constraints in supplementary databases and
    !           elsewhere - we do not want to return a value of Te
    !           or Ti less than an imposed minimum. If the values
    !           calculated for a step are less than that value then
    !           they are set equal to that value.
    !
    if (newt1e.lt.tfloor.and.newt1e.gt.0.0) newt1e = tfloor
    if (newt1i.lt.tfloor.and.newt1i.gt.0.0) newt1i = tfloor

    return
  end subroutine rkstep



  real*8 function fegrad(s,te,ti,n)
    use mod_solparams
    use mod_solcommon
    implicit none
    !
    !     This function returns the approximate derivate dte/ds and is
    !     used in the Runge-Kutta approximation method in attempting to
    !     estimate the functional value at the next step.
    !
    real*8 s,te,ti,n
    !
    !real*8 cond,estprad,estphelpi,estpei,estppelec,paes
    !external cond,estprad,estphelpi,estpei,estppelec,paes

    fegrad = (1.0/(k0e*te**2.5)) * (cond(s,te)+estprad(s,n,te) + estphelpi(s,n,te) + estpei(s,n,te,ti) + estepow(s) + estppelec(s) + paes(s) + qperpe(s,n,te,ti))

    !if (debug_s22) write(6,'(a,20(1x,g20.12))') 'Fegrad:',s,n,te,ti,gamma(s)/n,fegrad,&
    !     cond(s,te), estprad(s,n,te), estphelpi(s,n,te), estpei(s,n,te,ti),estppelec(s),paes(s),qperpe(s,n,te,ti),&
    !     (estprad(s,n,te) + estphelpi(s,n,te) + estpei(s,n,te,ti) + estppelec(s)),&
    !     (cond(s,te)+estprad(s,n,te) +  estphelpi(s,n,te) + estpei(s,n,te,ti) + estppelec(s) + paes(s) + qperpe(s,n,te,ti)),&
    !     (1.0/(k0e*te**2.5)), (1.0/(k0e*te**2.5))*&
    !     (cond(s,te)+estprad(s,n,te) +  estphelpi(s,n,te) + estpei(s,n,te,ti) + estppelec(s) + paes(s) + qperpe(s,n,te,ti)),&
    !     press(s,te,ti)

    return
  end function fegrad


  real*8 function figrad(s,te,ti,n)
    use mod_solparams
    use mod_solcommon
    implicit none
    !
    !     This function returns the approximate derivate dte/ds and is
    !     used in the Runge-Kutta approximation method in attempting to
    !     estimate the functional value at the next step.
    !
    real*8 s,te,ti,n

    !real*8 cond,conv,estpcx,estpei,estppion,pais
    !external cond,conv,estpcx,estpei,estppion,pais

    figrad = (1.0/(k0i*ti**2.5)) * (cond(s,ti)+ conv(s,n,ti) + estpcx(s,ti) - estpei(s,n,te,ti) + estipow(s) + estppion(s) + pais(s) + qperpi(s,n,te,ti))

    !if (debug_s22) write(6,'(a,20(1x,g20.12))') 'Figrad:',s,n,te,ti,gamma(s)/n,figrad,&
    !     cond(s,ti),conv(s,n,ti),estpcx(s,ti), -estpei(s,n,te,ti),estppion(s) ,pais(s), qperpi(s,n,te,ti), &
    !     estpcx(s,ti)-estpei(s,n,te,ti)+estppion(s),&
    !    (cond(s,ti) + conv(s,n,ti) + estpcx(s,ti) - estpei(s,n,te,ti) + estppion(s) + pais(s) + qperpi(s,n,te,ti)),&
    !    (1.0/(k0i*ti**2.5)), (1.0/(k0i*ti**2.5))*&
    !    (cond(s,ti) + conv(s,n,ti) + estpcx(s,ti) - estpei(s,n,te,ti) + estppion(s) + pais(s) + qperpi(s,n,te,ti)),&
    !    press(s,te,ti)

    return 
  end function figrad



  real*8 function fgrad(s,te,ti,n)
    use mod_solparams
    use mod_solcommon
    implicit none
    !
    !     This function returns the approximate derivate dt/ds for both species
    !     combined and is used in the Runge-Kutta approximation
    !     method in attempting to
    !     estimate the functional value at the next step.
    !
    real*8 s,te,ti,n

    !real*8 cond,estprad,estphelpi,estpcx,conv,paes,pais
    !real*8 estppelec,estppion
    !external cond,estprad,estphelpi,estpcx,conv,paes,pais
    !external estppelec,estppion

    fgrad = (1.0/(k0e*te**2.5+k0i*ti**2.5)) * (cond(s,te)+ cond(s,ti) &
         + conv(s,n,ti) + estprad(s,n,te) + estpcx(s,ti) &
         + estphelpi(s,n,te) + estppelec(s) + estppion(s) &
         + estepow(s) + estipow(s)&    ! Epow and Ipow need to be checked for sign conventions - both are added to their respective derivatives so this should be correct
         + paes(s) + pais(s))

    return
  end function fgrad

  real*8 function newn(s,te,ti,nimag,flag)
    use sol22_debug
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    !
    !     Calculates the density value from the given parameters by solving
    !     the quadratic equation for N - issues an error when part is
    !     imaginary and depending on options - may take corrective action.
    !
    real*8 s,te,ti,nimag
    integer flag

    real*8 imag,rest,postfact,tmpgam,tmpv,ptmp,tmppress,imag1,imag2

    !real*8 gamma,majrpos,ptmp,press
    !external gamma,press,majrpos

    flag = 0
    
    !if (s.eq.0) then
    !   newn = n0
    !   nimag = m0**2 * n0
    !   return
    !endif

    tmpgam = gamma(s)
    tmppress = press(s,te,ti)


    ! try scaling these by n0 so that the numbers aren't so large and numerical issues don't come into play
    imag1 = ((tmppress/(te+ti))/econv/n0) *  ((tmppress/(te+ti))/econv/n0)
    imag2 = - 4.0* ((((tmpgam/n0) * mconv) * (tmpgam/n0)) / econv) * (mb / (te+ti))
    imag = imag1 + imag2
    !imag = ((tmppress/(te+ti))/econv) **2/n0/n0 - 4.0* (((tmpgam * mconv) * tmpgam) / econv) * (mb / (te+ti))/n0/n0
    rest = (tmppress /(te+ti)) / econv/n0

    !if (s.eq.0) then 
    !   ! debug
    !   write(6,'(a,20(1x,g20.12))') 'NEWN:S=0:',s,imag,imag1,imag2,rest,n0,rest/n0, &
    !              tmppress,tmpgam,te,ti, &
    !              rest**2,-4.0*((tmpgam * mconv) / econv) &
    !              *(mb * tmpgam)/(te+ti),rest**2/n0/n0,&
    !              - 4.0* ((((tmpgam/n0) * mconv) * (tmpgam/n0)) / econv) * (mb / (te+ti)),&
    !              pinf,press(s,te,ti),pinf-press(s,te,ti),gamma0,gamma(s),gamma0-gamma(s)
    !   
    !
    !endif
    !
    if (imag.lt.0.0) then
       !
       !     set an error recovery condition ... if n becomes imaginary
       !     then set v = cs or v = lastvel and calculate n accordingly.
       !     OR adjust the pressure by adding enough so that the solution
       !     is not imaginary.
       !
       !
       if (velsw.eq.0) then
          tmpv =   sqrt( (te+ti)/mb * econv/mconv)
          newn =  abs( tmpgam / tmpv)
          nimag = newn
       elseif (velsw.eq.1) then
          tmpv =  lastvel
          newn =  abs( tmpgam / tmpv)
          nimag = newn
       elseif (velsw.eq.2) then

          ptmp = 4.0 * (econv*(te+ti)) * (mb * tmpgam) *  (tmpgam *mconv)
          padd = ptmp - tmppress + padd

          newn = ptmp / (econv*(te+ti))
          nimag= newn
          tmpv = tmpgam/newn

       elseif (velsw.eq.3) then

          ptmp = 4.0 * (econv*(te+ti)) * (mb * tmpgam) *  (tmpgam *mconv)
          !
          !            padd = ptmp - tmppress
          !
          newn = ptmp / (econv*(te+ti))
          nimag= newn
          tmpv = tmpgam/newn

       endif

       if (newn.lt.0.0) then
          write (6,*) 'NEWN<0 for Imaginary:Flow reversal?:',newn,tmpgam,tmpv,te,ti
          newn = -newn
          tmpv = -tmpv
       endif

       flag = 1

       !
       !        Debug
       !
       if (debug_sol22_on.and.pinavail) then
          write (6,'(a,10(1x,g14.6))') 'Newn:I',s,imag,rest, &
               tmpgam,te,ti, rest**2, &
               - 4.0*((tmpgam * mconv) / econv)*(mb * tmpgam)/(te+ti)
       endif


    else
       !
       !        Modify so that newn < 0 returns a flag / exit condition
       !
       if (founds.or.(lastiter.and.s.gt.lastiters)) then

          ! add scaling of n0 back in to obtain density
          newn = (rest + sqrt(imag))/2.0 * n0
          nimag = (rest - sqrt(imag))/2.0 * n0

          tmpv = tmpgam/newn

          if (newn.lt.0.0) then

             write (6,'(a,10(1x,g11.5))') 'Newn<0:',s,imag,rest, &
                  tmppress,tmpgam,te,ti, &
                  rest**2,-4.0*((tmpgam * mconv) / econv) &
                  *(mb * tmpgam)/(te+ti),newn 
             write(6,*) 'founds:',founds,lastiter,lastiters

             newn = abs(newn)
             flag = 2
             !return

          endif
       else

          ! add scaling of n0 back in to obtain density
          newn = (rest - sqrt(imag))/2.0 * n0
          nimag = (rest + sqrt(imag))/2.0 * n0

          tmpv = tmpgam/newn

          if (newn.lt.0.0) then

             write (6,*) 'Error: Newn < 0 :'
             write (6,'(a,20(1x,g13.5))') 'Newn:I',s,rest,imag,imag1,imag2, &
                  tmppress,tmpgam,te,ti, &
                  rest**2,-4.0*((tmpgam * mconv) / econv) &
                  *(mb * tmpgam)/(te+ti)
             write(6,*) 'founds:',founds,lastiter,lastiters

             newn = abs(newn)

             flag = 2
             !return

          endif
       endif
    endif

    if (debug_s22) then
       write(6,'(a,20(1x,g20.12))') 'NEWN:',s,rest,imag,imag1,imag2, &
            tmppress,tmpgam,te,ti, &
            rest**2/n0/n0,-4.0*((tmpgam * mconv) / econv) &
            *(mb * tmpgam)/(te+ti)/n0/n0,newn,tmpv,&
            flag
    endif


    return
  end function newn



  subroutine loadparms
    use mod_solrk

    !     include 'solrk'

    !     Loads Cash-Karp parameters

    implicit none

    !     Initialize

    integer i,j

    ai  = 0.0
    bij = 0.0
    ci  = 0.0
    cip = 0.0
    
    !do i = 1,6
    !   ai(i) = 0.0
    !   ci(i) = 0.0
    !   cip(i) = 0.0
    !   do j = 1,5
    !      bij(i,j) = 0.0
    !   end do
    !
    !end do

    ai(2) = 1.0/5.0
    ai(3) = 3.0/10.0
    ai(4) = 3.0/5.0
    ai(5) = 1.0
    ai(6) = 7.0/8.0

    bij(2,1) = 1.0/5.0
    bij(3,1) = 3.0/40.0
    bij(3,2) = 9.0/40.0
    bij(4,1) = 3.0/10.0
    bij(4,2) = - 9.0/10.0
    bij(4,3) = 6.0/5.0
    bij(5,1) = -11.0/54.0
    bij(5,2) = 5.0/2.0
    bij(5,3) = -70.0/27.0
    bij(5,4) = 35.0/27.0
    bij(6,1) = 1631.0/55296.0
    bij(6,2) = 175.0/512.0
    bij(6,3) = 575.0/13824.0
    bij(6,4) = 44275.0/110592.0
    bij(6,5) = 253.0/4096.0

    ci(1) = 37.0/378.0
    ci(3) = 250.0/621.0
    ci(4) = 125.0/594.0
    ci(6) = 512.0/1771.0

    cip(1) = 2825.0/27648.0
    cip(3) = 18575.0/48384.0
    cip(4) = 13525.0/55296.0
    cip(5) = 277.0/14336.0
    cip(6) = 1.0/4.0

    return



  end subroutine loadparms





  subroutine setsw(setind,pplasma,new_errlevel)
    use mod_solparams
    use mod_solswitch
    implicit none
    !     include 'solparams'
    !     include 'solswitch'

    !     SETSW: The purpose of this routine is to set the ACTIVE
    !            MODE switches that will be used for the specific
    !            run of SOL option 22. This allows Sol option 22
    !            to be run with different settings for each 1/2
    !            ring - thus allowing failed solutions to be replaced
    !            by simpler ones that will be more likely to work
    !            correctly and allowing the user to specify
    !            that problem rings be treated differently.

    !     Option - setind = 0 - set the switches to the input values
    !              setind = 1 - perform iterative error correction turning
    !                           off various options depending on the actswerror
    !                           switch and number of iterations.
    !              setind = 2 - Turn off all options and solve using conduction
    !                           ONLY


    !              errlevel=10- turn off equipartition if it was activated
    !              errlevel=9 - use uniform particles instead of d2n/dr2
    !              errlevel=8 - use only PINQI cooling contributions
    !              errlevel=7 - use 1/2 ring uniform power instead of whole
    !              errlevel=6 - use 1/2 ring uniform power + 1/2 ring particles
    !              errlevel=5 - 1/2 ring uniform particles + power at top
    !              errlevel=4 - 5 + turn off v^2 convection term
    !              errlevel=3 - 4 + no power terms
    !              errlevel=2 - 3 + no convective terms
    !              errlevel=1 - Conduction ONLY

    !     Default- setind = 0 - switches set to input values

    integer setind,pplasma,new_errlevel
    integer maxerrs
    parameter (maxerrs = 10)

    integer errlevel,errlevels(maxerrs)


    character*100 :: errlvltext(10)

    errlvltext(10)='Turn off equipartition if it was activated'
    errlvltext(9)='Use uniform particles instead of d2n/dr2'
    errlvltext(8)='Use only PINQI cooling contributions'
    errlvltext(7)='Use 1/2 ring uniform power instead of whole'
    errlvltext(6)='Use 1/2 ring uniform power + 1/2 ring particles'
    errlvltext(5)='1/2 ring uniform particles + power at top'
    errlvltext(4)='5 + turn off v^2 convection term'
    errlvltext(3)='4 + no power terms'
    errlvltext(2)='3 + no convective terms'
    errlvltext(1)='Conduction ONLY'
    
    new_errlevel = -1
    errlevel = 0

    !write(6,*) 'ERR:',pplasma,setind
    errlevels(10)= 9
    errlevels(9) = 8
    errlevels(8) = 7
    errlevels(7) = 6
    errlevels(6) = 5
    errlevels(5) = 4
    errlevels(4) = 3
    errlevels(3) = 2
    errlevels(2) = 1
    errlevels(1) = 0

    !        Load input values

    if (setind.eq.-1) then

       call setallsw(pplasma,0)

       !        Set error level

    elseif (setind.eq.-2.or.setind.gt.0) then
       if (setind.eq.-2) then
          errlevel = actswerror
       else
          errlevel = setind

       endif

       !        Initially set all quantities to INPUT values

       write(6,*) 'ERR:BEG:',pplasma,setind,errlevel,actswerror

       !        ERRLEVEL = 10 = Switch OFF equipartition if it was
       !                       turned ON

       call setallsw(pplasma,0)

       if (errlevel.eq.10) then
          if (actswpei.eq.1.0) then
             actswpei = 0.0
          else
             errlevel = errlevels(errlevel)

          endif

          !        ERRLEVEL = 9 = Switch to uniform gperp from d2n/dr2

       endif
       !write(6,*) 'ERR:I1 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp

       !           Turn OFF equipartition (from level 10)

       if (errlevel.eq.9) then

          !           Replace gradient proportional options with whole ring uniform

          actswpei = 0.0

          if (actswgperp.eq.7.0.or.actswgperp.eq.8.0) then

             actswgperp = 2.0

          else

             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 8 = Switch OFF any ION heating power

       !write(6,*) 'ERR:I2 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp,actswpcx

       !           Turn OFF equipartition

       if (errlevel.eq.8) then

          !           Replace whole ring uniform with 1/2 ring uniform

          actswpei = 0.0

          !           Turn OFF PINQI heating if active

          if (actswgperp.eq.7.0.or.actswgperp.eq.8.0)actswgperp = 2.0
          if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) then
             actswpcx= 5.0
          else
             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 7 = Eliminate whole ring uniform power options -

       !write(6,*) 'ERR:I3 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp,actswpcx

       !           Turn OFF equipartition

       if (errlevel.eq.7) then

          !           Replace whole ring gradient proportional particles with
          !           whole ring uniform option.

          actswpei = 0.0

          !           Eliminate PINQI heating

          if (actswgperp.eq.7.0.or.actswgperp.eq.8) actswgperp = 2.0

          !           Eliminate whole ring uniform power options - replace with
          !           1/2 ring uniform equivalents.

          if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
          if (actswpow.eq.4.0) then
             actswpow = 1.0
          elseif (actswpow.eq.6.0) then
             actswpow = 5.0
          else
             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 6 - 1/2 ring uniform target power + particles

       write(6,*) 'ERR:I4 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp

       !           Turn OFF equipartition

       if (errlevel.eq.6) then

          !           Eliminate PINQI heating

          actswpei = 0.0

          !           Replace power options with uniform 1/2 ring unless
          !           already set

          if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0

          if (actswpow.ne.1.0.or.actswgperp.ne.1.0) then
             if (actswpow.ne.0.0) actswpow = 1.0

             actswgperp = 1.0

          else

             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 5 - 1/2 ring uniform particles - all power in at top

       !write(6,*) 'ERR:I5 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp


       if (errlevel.eq.5) then

       !           Turn OFF equipartition
          actswpei = 0.0


          !           Eliminate PINQI heating

          if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0

          !           Replace power options with uniform 1/2 ring unless
          !           already set

          if (actswpow.ne.0.0.or.actswgperp.ne.1.0) then
             actswpow   = 0.0

             actswgperp = 1.0

          else

             errlevel = errlevels(errlevel)

          endif

          ! Turn off perpendicular power compensation terms
          actswqperpe = 0.0
          actswqperpi = 0.0
          
       endif


       !        ERRLEVEL = 4

       !        Turn off second convective term as well as level 5+

       !write(6,*) 'ERR:I6 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp


       !           Turn OFF equipartition

       if (errlevel.eq.4) then

          !           Eliminate PINQI heating

          actswpei = 0.0

          !           Replace power option - in at top
          !                   particle option - 1/2 ring uniform

          if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
          actswpow   = 0.0

          ! Turn off perpendicular power compensation terms
          actswqperpe = 0.0
          actswqperpi = 0.0

          actswgperp = 1.0

          if (actswconv.ne.0.0) then

             actswconv = 0.0

          else

             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 3 - All of above + power terms turned off

       !write(6,*) 'ERR:I7 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp

       !           Reset switches for other error levels

       if (errlevel.eq.3) then
          actswpei   = 0.0
          actswgperp = 1.0
          actswpow   = 0.0
          ! Turn off perpendicular power compensation terms
          actswqperpe = 0.0
          actswqperpi = 0.0


          !           Switch off power terms if necessary

          actswconv  = 0.0

          if (actswphelp.ne.0.0.or.actswpcx.ne.0.0.or.actswprad.ne.0.0.or.actswppion.ne.0.0.or.actswppelec.ne.0.0) then
             actswphelp = 0.0
             actswpei   = 0.0
             actswpcx   = 0.0
             actswprad  = 0.0
             actswppelec= 0.0

             actswppion = 0.0

          else

             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 2 - All of above + no convection terms

       !write(6,*) 'ERR:I8 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp
       if (errlevel.eq.2) then
          actswgperp = 1.0
          actswphelp = 0.0
          actswpei   = 0.0
          actswpcx   = 0.0
          actswppelec= 0.0
          actswppion = 0.0
          actswprad  = 0.0
          actswpow   = 0.0
          ! Turn off perpendicular power compensation terms
          actswqperpe = 0.0
          actswqperpi = 0.0


          actswconv  = 0.0

          if (actswcond.ne.0.0) then

             actswcond = 0.0

          else

             errlevel = errlevels(errlevel)

          endif

       endif

       !        ERRLEVEL = 1 - ALL OFF

       !write(6,*) 'ERR:I9 :',pplasma,setind,errlevel,actswerror,actswpei,actswgperp

       !           SET switches to conduction only

       if (errlevel.eq.1) then

          call setallsw(pplasma,1)

          !        SET the error switch to be the next in sequence of the error
          !        conditions to be used - the last is zero - meaning no more error
          !        provisions are available.

       endif
       new_errlevel = errlevel

       actswerror = errlevels(errlevel)

       if (errlevel.lt.10.and.errlevel.gt.0) then 
          write(6,*) 'ERR:END:',pplasma,setind,errlevel,actswerror,trim(errlvltext(errlevel+1))
       else
          write(6,*) 'ERR:END:',pplasma,setind,errlevel,actswerror
       endif
       !        Load maximum error settings

    elseif (setind.eq.0) then

       call setallsw(pplasma,1)

       !     Exit

    endif
    return



  end subroutine setsw


  subroutine setallsw(pplasma,ind)
    use mod_solparams
    use mod_solswitch
    implicit none
    !     include 'solparams'
    !     include 'solswitch'

    !     SETALLSW: This routine sets all of the switches to either their
    !               input values or to all OFF except conduction. The
    !               other options in the SETSW routine will tweak these
    !               to obtain other error solution conditions betweeen
    !               these two extremes.


    integer ind,pplasma


    if (ind.eq.0) then

       !        Set the switches differently for main SOL and private plasma

       !        Main SOL


       if (pplasma.eq.0) then
          !        Set all values to INPUT switch settings
          actswion  = switch(swion)
          actswioni = switch(swioni)
          actswcond = switch(swcond)
          actswconv = switch(swconv)
          actswprad = switch(swprad)
          actswphelp= switch(swphelp)
          actswpei  = switch(swpei)
          actswpcx  = switch(swpcx)
          actswppelec = switch(swppelec)
          actswppion  = switch(swppion)

          ! set external power switches
          actswepow = switch(swepow)
          actswipow = switch(swipow)
          
          actswppress = switch(swppress)

          !           PCX - DIVIMP PINQI - sub-option switches

          actswqidatiz = switch(swqidatiz)
          actswqidmliz = switch(swqidmliz)
          actswqidcx = switch(swqidcx)

          actswqidrec= switch(swqidrec)
          actswvisc1= switch(swvisc1)
          actswmach = switch(swmach)
          actswnmom = switch(swnmom)
          actswe2d  = switch(swe2d)
          actswpow  = switch(swpow)
          actswgperp= switch(swgperp)
          actswmajr = switch(swmajr)
          actswcore = switch(swcore)
          actswrecom= switch(swrecom)
          actswsmooth= switch(swsmooth)

          actswqperpe = switch(swqperpe)
          actswqperpi = switch(swqperpi)

          
          actswerror = switch(swerror)

       elseif (pplasma.eq.1) then
          !        Private plasma

          if (switch(swion).eq.1.0.or.switch(swion).eq.2.0.or.switch(swion).eq.8.0) then
             actswion  = switch(swion)
             actswioni = switch(swionp)
          else
             actswion  = switch(swionp)
             actswioni = switch(swioni)

          endif
          actswcond = switch(swcond)
          actswconv = switch(swconv)
          actswprad = switch(swprad)
          actswphelp= switch(swphelp)
          actswpei  = switch(swpei)

          !           Both of these options are always off for now in the PP - at
          !           some time they may be active but only (possibly) to distribute
          !           the power transferred from the main SOL - however - this is
          !           usually modelled using the already existing private plasma
          !           power distribution methods.

          actswpcx  = switch(swpcx)
          actswppelec = 0.0
          actswppion  = 0.0

          ! set external power switches - leave on in PP
          actswepow = switch(swepow)
          actswipow = switch(swipow)
          
          !            actswppelec  = switch(swppelec)
          !            actswppion   = switch(swppion)
          !            actswppress = switch(swppress)

          !           PCX - DIVIMP PINQI - sub-option switches

          actswppress = 0.0
          actswqidatiz = switch(swqidatiz)
          actswqidmliz = switch(swqidmliz)
          actswqidcx = switch(swqidcx)

          actswqidrec= switch(swqidrec)
          actswvisc1= switch(swvisc1)
          actswmach = switch(swmach)
          actswnmom = switch(swnmom)

          !           Set power to private plasma value

          actswe2d  = switch(swe2d)
          actswpow  = switch(swpowp)
          actswgperp= switch(swgperpp)
          actswmajr = switch(swmajr)
          actswcore = switch(swcore)
          actswrecom= switch(swrecom)
          actswsmooth= switch(swsmooth)

          actswqperpe = switch(swqperpe)
          actswqperpi = switch(swqperpi)


          actswerror = switch(swerror)

       endif

    elseif (ind.eq.1) then   ! turns all switches off

       !        Ionization Source is EXPONENTIAL
       actswion= 0.0


       !        Initial Ionization is set to exponential
       actswioni = 0.0


       !        First convection term is OFF
       actswcond = 0.0


       !        Kinetic convection term is OFF
       actswconv = 0.0


       !        Radiation is OFF
       actswprad = 0.0


       !        Phelpi is OFF
       actswphelp= 0.0

       ! External power terms are off
       actswepow = 0.0
       actswipow = 0.0
       
       !        Pei is OFF
       actswpei = 0.0


       !        Pcx is OFF
       actswpcx = 0.0

       !           PCX - DIVIMP PINQI - sub-option switches - all OFF
       actswqidatiz = 0.0
       actswqidmliz = 0.0
       actswqidcx = 0.0
       actswqidrec= 0.0

       !        Private plasma loss compensation options are OFF
       actswppelec = 0.0
       actswppion = 0.0
       actswppress = 0.0


       !        Viscosity is OFF
       actswvisc1 = 0.0


       !        Mach solver is OFF
       actswmach = 0.0


       !        Neutral Momentum Loss is OFF
       actswnmom = 0.0


       !        Edge2D Compatibility is set to INPUT value
       actswe2d = switch(swe2d)


       !        Distributed Power is OFF
       actswpow = 0.0


       !        Cross-field correction is OFF
       actswgperp = 0.0


       !        Major Radius Correction is OFF
       actswmajr = 0.0


       !        Core Option is OFF
       actswcore = 0.0

       !        Recombination option is OFF
       actswrecom = 0.0

       !        Smoothing is set to INPUT Value
       actswsmooth = switch(swsmooth)

       ! qperp power flux terms are off
       actswqperpe = 0.0
       actswqperpi = 0.0


       !        Error handling switch is now OFF
             actswerror = 0.0

    endif
    return



  end subroutine setallsw







end module mod_sol22_solver
