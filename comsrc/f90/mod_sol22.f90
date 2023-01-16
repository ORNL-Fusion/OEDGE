module mod_sol22



  !
  ! This module implements and streamlines the interface to SOL22 - 1D fluid
  ! equation solver that works from the specified target conditions.

  ! Required inputs:
  !
  ! target conditions - ne, Te, Ti - Cs calculated from Te,Ti - E field derived from plasma solution
  ! coordinates - set of coordinate distances in meters from the target identifying the locations where the plasma conditions
  !               need to be returned
  ! options/switches/parameters - input that defines the behaviour of each component of the solver
  !         - particle sources and sinks
  !         - energy/power sources and sinks
  !         - momentum sources/sinks - usually implemented through changes in the total pressure
  !
  ! automatic error correction
  !
  ! The code is based on SOL22 in OEDGE which is extremely large (the goal with this is to segment it from OEDGE
  ! and allow for standalone
  ! functionality. If Neutral code output is not available, the arrays for these will not be allocated and options
  ! will not require these.
  !
  ! Switch/parameter inputs will be read in from named files which will allow the solver options to change depending
  ! on selections made
  ! in the OEDGE/LIM input files. This will allow for different ionization options/decay lengths etc for different parts
  ! of the LIM mesh
  ! or OEDGE grid for example.
  !
  ! This code is initially intended for implementation in LIM. However, the goal is to clean it up and streamline it so
  ! that it can be put
  ! back into OEDGE as a cleaner and easier to understand implementation. 
  !

  use mod_sol22_solver
  use mod_sol22_plots
  use mod_sol22_output
  !use mod_sol22_support
  use mod_sol22_utils
  use mod_sol22_sources
  
  implicit none



contains



  !
  !     program calcsol
  !
  subroutine calcsol (spts,npts,errcode,serr,te,ti,ne,vb,exp_press,act_press,prad,ir,irsep,int_powrat,cprint,&
       cve,cvi,cde,cdi) 
    !use sol22_input
    use sol22_debug
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    use error_handling
    !use mod_solrk
    !use mod_params
    !use mod_slcom
    implicit none

    !COMMON /POWERFLOW/ cve        ,cvi        ,cde        ,cdi
    REAL               cve(*),cvi(*),cde(*),cdi(*)

    INTEGER     eflag    
    REAL*8      pmloss(MXSPTS)
    CHARACTER*2 note(MXSPTS)

    integer errcode,npts,ir,irsep,cprint
    real*8 serr
    real*8 spts(mxspts)
    real*8 te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts), exp_press(mxspts),act_press(mxspts),prad(mxspts)
    real*8 int_powrat(3)
    !
    !     CALCSOL:
    !
    !     This program employs a simple Runge-Kutta technique and several
    !     approximations
    !     to calculate the temperature, density and velocity distribution
    !     along the scrape off layer.
    !
    !     Modified to use adaptive step-size control and Cash-Karp
    !     parameters for embedded Runge-Kutte Method
    !
    !     David Elder, Aug 3, 1994
    !
    !
    real*8 sinit,send
    !
    !      real*8 h,s,ntmp
    !
    integer    i,k,flag,imflag,lastflag,vcount
    integer    loopstart
    integer    ik
    real*8    n,cs
    real*8    t1e,t1i,v1
    real*8    fval

    !real*8    fegrad,figrad,newn,gamma,peiupdt,pintupdt
    !real*8    phelpiupdt,pcxupdt,smomtmp
    !real*8    nimag,pmomloss,ptmp,pradupdt,scx_s,scxtmp
    !real*8    dm0,iters,vgradval,press,scxupdt
    !real*8    estppion,estppelec

    real*8    nimag,ptmp,dm0,iters,smomtmp,scxtmp

    !external   fegrad,figrad,newn,gamma,peiupdt,phelpiupdt,pcxupdt
    !external   vgradval,pintupdt,pmomloss,press,pradupdt
    !external   scx_s,scxupdt
    !external  estppion,estppelec
    !real*8,external :: srci,srcf

    real*8    ga(mxspts),vb2(mxspts)
    real*8    ne2(mxspts),vsupers(mxspts),vsubs(mxspts)
    real*8    vsound(mxspts),scxv(mxspts)
    real*8    peiv(mxspts),pcxv(mxspts),phelpiv(mxspts)
    real*8    peiv2(mxspts),pradv(mxspts)
    real*8    epowv(mxspts),ipowv(mxspts)
    
    real*8    pir(mxspts),pii(mxspts),vgradn(mxspts)

    real*8    condi(mxspts),conde(mxspts)
    real*8    convi(mxspts),conve(mxspts)
    real*8    int_powe,int_powi,int_conde,int_condi
    !
    real*8 :: tmpgam,tmppress,tmpn

    !
    !     Function declarations
    !
    !real*8 cond,conv,paes,pais
    !external cond,conv,paes,pais
    !
    !      real*8    tpress(mxspts)
    !      real*8    errmax, errti,errte
    !      real*8    newv1,errvb,vrat
    real*8    tauii
    !      real*8    hnew
    !      real*8    vsep,vsub,vsup,lastvsup,v2
    !
    integer      exitcond,negerrflag
    character*80 comment
    !
    !      logical      imag
    !
    !     logical      errneg
    !
    !     define arrays for k-values
    !
    !      real*8 ke(6),ki(6)
    !
    !     Initialization
    !
    ! slmod begin - new
    DATA t1i,t1e /0.0,0.0/
    !

    !
    ! Check for npts = 0 - should not happen but if it does issue an error message and return to caller
    !
    if (npts.eq.0) then
       call errmsg('MOD_SOL22:CALCSOL: NPTS is ZERO on entry - RETURNING - NPTS = ',npts)
       return
    endif
    
    ! set the Sol22 print option to match the option passed to calcsol - sol22_cprint is being replaced by a dedicated sol22_print
    ! option (TAG 299) (Keep sol22_cprint until all references have been removed)
    sol22_cprint = cprint
    ! over-write the value of sol22_print if cprint = 9 (all output) has been specified
    if (cprint.eq.9) then
       sol22_print = 2
    endif
    
    !
    !     jdemod - setting simag1 to the ring length should be big enough
    !
    !      simag1 = HI
    simag1 = ringlen
    !
    simag2 = halfringlen
    !ierror = MAXNKS
    eflag  = 0

    CALL DZero(conve ,MXSPTS)
    CALL DZero(convi ,MXSPTS)
    CALL DZero(conde ,MXSPTS)
    CALL DZero(condi ,MXSPTS)
    CALL DZero(pmloss,MXSPTS)
    ! slmod end
    !

    !
    !     Set up the value of startn for this ring ...
    !
    if (actswe2d.ne.0.0) then
       startn = ike2d_start
    else
       startn = 1
    endif

    !write(6,*) 'CALCSOL: Errlevel = ',actswerror
    !
    !     Set the minimum allowed temperature for the solver.
    !
    if (tmin.gt.0.0) then
       tfloor = tmin
    elseif (tmin.lt.0.0) then
       tfloor = min(abs(tmin),min(te0,ti0))
    else
       tfloor = 0.0
    endif

    timax = 0.0
    temax = 0.0

    call loadparms
    !
    !     Record the target density.
    !
    netarg = n0
    nefinal= n0

    n1center= n1
    n1final = n1

    errcode = 0
    serr = 0.0

    flag = 0
    imflag = 0
    lastflag = 0
    !      imag = .false.
    lastiter = .true.
    lastiters = 0.0
    founds = .true.
    !
    !     Recalculate the initm0 value (the value the solver starts at)
    !     for the case of Edge2D compatibility options where the solver
    !     starts at the middle of the first cell and the velocity
    !     there and thus the Mach number are specified from the
    !     Edge2D data.
    !
    if (actswmach.eq.3.0) then

       cs = -sqrt((te0+ti0)/mb * econv/mconv)

       if (cs.ne.0.0) then
          initm0 = abs(v0/cs)
       else
          initm0 = 1.0
       endif

       write (6,'(a,6g12.4)') 'Initm0:', te0,ti0,mb,v0,cs,initm0

    else

       initm0 = 1.0

    endif

    e2dm0 = 1.0

    if (actswe2d.ne.0.0) then

       cs = -sqrt((te1+ti1)/mb * econv/mconv)

       if (actswe2d.eq.1.0) then
          e2dm0 = vpe2d / cs
       elseif (actswe2d.eq.2.0) then
          e2dm0 = v1e2d / cs
       elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.actswe2d.eq.9.0) then
          e2dm0 = vpg / cs
       endif

    endif
    !
    !     Copy Spts for use in some routines - since it involves a lot
    !     code rewriting to move spts to a common block or to pass it
    !     to the routine srci - where it is needed.
    !
    !     Need to set values in ALL arrays for the npts+1 position
    !     Just copy over the values from npts.
    !     So ...
    !     Leave out calculating all the way to the mid-point.
    !
    !     spts(npts+1) is set to S-value at the boundary of the cell 
    !     in the calcsol_interface routine - this is the same as the 
    !     halfringlen value but is NOT the same as ringlen/2.0
    !
    nptscopy = npts
    do i = 1,nptscopy+1
       sptscopy(i) = spts(i)
    end do
    !
    !     Read in the parameters for the ring from a file
    !
    !      call readsol (spts,npts)
    !
    !     When operating in the mode that iteratively finds the lowest Mach
    !     number at the target that does not result in imaginary values
    !     of density or velocity when solving the equations ... the code
    !     loops back to this point ... and recalculates from here.
    !
    stopimag = .false.
    miter = 0
    !
    m0 = initm0
    m1 = e2dm0
    !
    if (actswmach.eq.1.0.or.actswmach.eq.2.0) then
       dm0 = deltam0
       origm0 = m0
       lastm0 = m0
       stopimag = .true.
       lastiter = .false.
    endif
    !

1000 continue

    miter = miter + 1

    if ((actswmach.eq.1.0.or.actswmach.eq.2.0).and.miter.gt.1) then
       lastm0 = m0
       m0 = m0 + dm0
       !
       !        Change initial target density if Mach solver option
       !        is 2.0

       if (actswmach.eq.2.0) then

          n0 = netarg * initm0 / m0
          nefinal = n0

          if (actswe2d.ne.0.0) then
             n1 = n1center * initm0 / m0
             n1final = n1
          endif

       endif

       if (miter.gt.maxiter) then
          write(6,*) 'Error: Program exceeds maximum iterations calculating m0'
          write(6,*) 'Iterations: ',miter-1,' m0 = ',m0
          stop
       elseif (m0.le.0.0) then
          write(6,*) 'Error: m0 less than zero :'
          write(6,*) 'Iterations: ',miter-1,' m0 = ',m0
          stop
       endif
       write(6,*) 'Svals:', miter,lastiters,m0,lastm0,exitcond
       lastflag = imflag
       imflag = 0
    endif

    if ((actswmach.eq.1.0.or.actswmach.eq.2.0).and.(m0.ne.origm0)) then
       founds = .false.
    else
       founds = .true.
    endif
    !      
    !     Set various values depending on switches ... for example the gamma
    !     factor for power to the plates is different depending on options.
    !     Pre-calculate the integral of the numerical ionization source - if
    !     there is one.
    !
    call initval
    ! slmod end 
    ! jdemod - moved call to sol22headers from the end of initval - moved to the interface routine
    !CALL SOL22Headers

1500 continue
    !
    !     Set up the initial values
    !
    if (spts(1).eq.0.0) then
       te(1) = te0
       ti(1) = ti0
       ne(1) = n0
       ne2(1) = m0**2 * n0
       vb(1) = v0
       ga(1) = n0 * v0
       pir(1) = -0.4 * n0 * econv * ti0
       pii(1) = -0.4 * n0 * econv * ti0

       exp_press(1) = pinf
       act_press(1) = pinf

       vgradn(1) = vgradval(spts(1),0)
       loopstart = startn + 1

       !write(6,'(a,20(1x,g12.5))') 'sol22 calcsol:', ne(1),te(1),ti(1),vb(1),ga(1),act_press(1)

    elseif (actswe2d.ne.0.0) then

       te(ike2d_start) = te1
       ti(ike2d_start) = ti1
       ne(ike2d_start) = n1
       ne2(ike2d_start) = m1**2 * n1
       !
       !        Different E2D options
       !
       if (actswe2d.eq.1.0) then
          vb(ike2d_start) = vpe2d
          ga(ike2d_start) = n1 * vpe2d
          v1 = vpe2d
       elseif (actswe2d.eq.2.0) then
          vb(ike2d_start) = v1e2d
          ga(ike2d_start) = n1 * v1e2d
          v1 = v1e2d
       elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.actswe2d.eq.9.0) then
          vb(ike2d_start) = vpg
          ga(ike2d_start) = n1 * vpg
          v1 = vpg
       endif

       pir(ike2d_start) = -0.4 * n1 * econv * ti1
       pii(ike2d_start) = -0.4 * n1 * econv * ti1

       exp_press(ike2d_start) = pinf
       act_press(ike2d_start) = pinf

       vgradn(ike2d_start) = vgradval(spts(ike2d_start),0)

       loopstart = startn+1

    else
       loopstart = startn
    endif
    !
    !     More initialization
    !
    vcount = 0

    vsepmin = abs(v0)
    ssepmin = 0.0
    nlast = n0
    slast = 0.0
    !
    !     If the options are set - initialize the UPDATE/ESTIMATE
    !     integral routines for S=Soffset..
    !
    if (actswnmom.eq.4) then
       if (actswe2d.eq.0) then
          smomtmp = pmomloss(soffset,1,v0,te0,ti0)
       else
          smomtmp = pmomloss(soffset,1,v1,te1,ti1)
       endif
       ptmp = press(soffset,te0,ti0)

       write (6,*) 'Smomtmp1:',smomtmp,ptmp

    elseif (actswnmom.eq.9.or.actswnmom.eq.10) then
       !
       !     WF'96: Initialize CX momentum loss routine
       !
       if (actswe2d.eq.0) then
          scxtmp = scxupdt(soffset,n0,n0,ti0,ti0)
       else
          scxtmp = scxupdt(soffset,n1,n1,ti1,ti1)
       endif
       ptmp = press(soffset,te0,ti0)
       !
       !         write (6,*) 'Smomtmp1:',scxtmp,ptmp
       !
    endif
    !
    !     The actual values of the arguments (except soffset) are
    !     irrelevant. These lines simply garantee the initialization
    !     of the integration routines for this 1/2 ring.
    !
    ptmp =  peiupdt(soffset,ne(1),ne(1),t1e, t1i,t1e,t1i,fval)
    ptmp =  pradupdt(soffset,ne(1),ne(1),t1e,t1e)
    ptmp = pcxupdt(soffset,t1i,t1i)
    ptmp = phelpiupdt(soffset,ne(1),ne(1),t1e,t1e)

    if (debug_s22) write(6,*) 'Init_power:', peiupdt(soffset,ne(1),ne(1),t1e, t1i,t1e,t1i,fval),&
         pradupdt(soffset,ne(1),ne(1),t1e,t1e),pcxupdt(soffset,t1i,t1i),phelpiupdt(soffset,ne(1),ne(1),t1e,t1e)


    ! slmod begin - new
    simag   = 0.0D0
    snegerr = 0.0D0
    ! slmod end
    !
    !     Loop through the S values
    !
    do i = loopstart,npts

       if (i.eq.1) then
          sinit = 0.0
          send = spts(i)
          t1e = te0
          t1i = ti0
          n = n0
          v1 = v0
       else
          sinit = spts (i-1)
          send = spts(i)
          t1e = te(i-1)
          t1i = ti(i-1)
          n = ne(i-1)
          v1 = vb(i-1)
       endif

       if (debug_sol22_on) then 
          !write(0,*) 'mod_sol22:',debug_sol22_on
          call save_s22_data(dble(i),sinit,n,t1e,t1i,v1,gamma(sinit),srci(sinit),srcf(sinit),press(sinit,t1e,t1i))
       endif
       !
       !        Inner loop for the steps between each S-value.
       !
       timax = max(timax,t1i)
       temax = max(temax,t1e)

       if (debug_s22) then
          tmppress = press(sinit,t1e,t1i)
          tmpgam = gamma(sinit)
          tmpn = newn(sinit,t1e,t1i,nimag,flag)
          write(6,'(a,20(1x,g20.12))') 'S1:',sinit,send,t1e,t1i,n,v1,n*v1,tmpgam,tmppress,tmpn,tmpn-n,nimag,flag
       endif

       call solvstep(sinit,send,t1e,t1i,n,exitcond,imflag,negerrflag,vcount)

       ! slmod begin - new
       !         CALL NoName01(imflag, i,  loopstart,spts,npts,
       !     .                 exitcond,  note,te,ti,ne,ne2,ga,vb,vb2,
       !     .                 pir,pii,vgradn,act_press,exp_press,
       !     .                 vsound,vsubs,vsupers,peiv,peiv2,pradv,pcxv,
       !     .                 phelpiv,scxv,convi,conve,condi,conde,errcode,
       !     .                 *2000)
       !         
       !         WRITE(0,*) 'Continuing nicely'
       ! slmod end
       !

       if (imflag.eq.1) then
          !
          !           An imaginary value was encountered - this value of
          !           the error code will be overwritten if a
          !           more serious errcode comes along.
          !
          ! The solver can often encounter an imaginary result for S=0 as the step size is adjusted.
          ! This is now not reported as an error. 
          ! 

          if (simag.eq.0.0) then
             errcode = 0
          else 
             errcode = 1
          endif

          serr = simag
          
          !     
          !           jdemod - copy the code setting simag1 and simag2 from NoName
          !
          IF (simag1.EQ.ringlen) THEN
             simag1 = simag - soffset
             simag2 = simag - soffset
          ELSE
             simag2 = simag - soffset
          ENDIF

       endif
       !
       !        Exit conditions 1 and 2 do not generate error codes
       !        because they are a normal part of the mach number
       !        solver - if imaginary results lead to the situation
       !        but an eventual solution is found - then errcode will
       !        have been set to 1 in the above statement.
       !
       if (exitcond.eq.1) then
          !
          !           Iterate at larger mach number
          !

          goto 1000

       elseif (exitcond.eq.2) then
          !
          !           Iterate at smaller mach number
          !
          if (dm0.le.m0res) then
             !
             !              If one has hit the mach number resolution
             !              limit then run through the entire
             !              ring solving one last time.
             !
             lastiter = .true.
             founds = .false.
             !
             !              Set M0 to the previous smaller value.
             !
             !               m0 = lastm0
             !
             write(6,*) 'Last iteration:',lastiters
             goto 1500
          else
             !
             !              Rest to last non-working m0 and reduce the
             !              delta m0 value by a factor of 10. - then
             !              iterate.
             !
             dm0 = dm0/10.0
             if (dm0.lt.m0res) dm0 = m0res
             m0 = lastm0
             goto 1000
          endif
       elseif (exitcond.eq.3) then
          !
          !           Negative temperature encountered - unexpected exit.
          !
          errcode = exitcond
          serr    = snegerr
          write(6,*) 'Exitcond = 3:',snegerr

          return

       elseif (exitcond.eq.4) then
          !
          !           Step size too small - unexpected exit.
          !           Possibly caused in response to T -> 0
          !
          !           Step size has hit lower limit and no other
          !           actions are available to try to push
          !           solution forward.
          !
          errcode = exitcond
          serr    = snegerr
          write(6,*) 'Exitcond = 4:',snegerr
          !
          ! slmod begin - new
          GOTO 2000
          !
          !            return
          ! slmod end
          !
       elseif (exitcond.eq.5) then
          !
          !           Temperature along ring dropped to less
          !           than dropfrac*tmax OR temperature hit
          !           minimum level limit for ring.
          !
          errcode = exitcond
          serr    = snegerr
          write(6,*) 'Exitcond = 5:',snegerr
          !
          ! slmod begin - new
          GOTO 2000
          !
          !            return
          ! slmod end
          !
       elseif (exitcond.eq.6) then
          !
          !           Negative density encountered in solver
          !           Code will exit - if the error corrector is OFF -
          !           it will be turned ON for this half-ring ONLY -
          !           an error message will be recorded to thsi effect.
          !           If the error corrector is ON - it will proceed
          !           normally.
          !
          errcode = exitcond
          serr    = snegerr
          write(6,*) 'Exitcond = 6:',snegerr
          !
          ! slmod begin - new
          GOTO 2000
          !
          !            return
          ! slmod end
          !
       elseif (exitcond.eq.7) then
          !
          !           NaNQ encountered in solver - usually indicative of
          !           a negative temperature condition.
          !
          !           Code will exit - if the error corrector is OFF -
          !           it will be turned ON for this half-ring ONLY -
          !           an error message will be recorded to this effect.
          !           If the error corrector is ON - it will proceed
          !           normally.
          !
          errcode = exitcond
          serr    = snegerr
          write(6,*) 'Exitcond = 7:',snegerr
          !
          ! slmod begin - new
          GOTO 2000
          !
          !            return
          ! slmod end
          !
       endif
       !
       !       Step was successfully solved - run through post processing
       !       and record all the quantities for the cell.
       !
       !       Tabulate the n,te,ti values for the grid point.
       !

       !if (i.gt.3) stop 'debug to this point - mod_sol22.f90'

       te(i) = t1e
       ti(i) = t1i

       ne(i) = newn(spts(i),t1e,t1i,nimag,flag)
       ne2(i) = nimag

       ga(i) = gamma(spts(i))
       vb(i) = ga(i) / ne(i)
       if (ga(i).eq.0.0) then
          vb2(i) = 0.0
       else
          vb2(i) = ga(i) / ne2(i)
       endif

       if (sol22_cprint.eq.3.or.sol22_cprint.eq.9.or.debug_s22) then
          ! The same info gets written for each condition - int(i/(npts/100)) needs to be evaluated as real just in case
          ! int(npts/100) = 0 
          if (i.lt.100.or.i.eq.int(i/real((real(npts)/100.0)))*int(npts/100)) then
             write(6,'(a,i8,9(1x,g12.5),2i2)') 'Step:',i,m0,send,te(i),ti(i),ne(i),vb(i),ga(i),ne(i)*vb(i),ionsrc(i),flag,negerrflag
          endif
       endif
       !
       !        Calculate viscosity estimate and limit terms
       !
       pir(i) = -0.4 * ne(i) * econv * ti(i)
       tauii = 2.5d0 * 2.09d13*(ti(i))**1.5* sqrt(mb) / (ne(i)  * 15.0d0)

       pii(i) = - 4.0d0/9.0d0 * ne(i) * econv * ti(i) * tauii * vgradval(spts(i),0)
       vgradn(i) = vgradval(spts(i),0)
       !
       !        Calculate current pressure
       !
       act_press(i) = ne(i) * econv * (te(i) + ti(i)) + mb * mconv * ne(i) * vb(i)**2
       !
       !        Calculate Expected pressure
       !
       !
       !         make sure pmomloss set up correctly for call from pressure
       !
       if (actswnmom.eq.4) then
          if (actswe2d.eq.0) then
             smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
          else
             smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
          endif
       endif

       exp_press(i) = press(spts(i),te(i),ti(i))

       if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then
          if (i.lt.100.or.i.eq.int(i/real((real(npts)/100.0)))*int(npts/100)) then
             write (6,'(a,8x,6(1x,g12.5))') 'Pres:',pinf,exp_press(i)-pinf-padd,exp_press(i)-pinf,exp_press(i),act_press(i),padd
          endif
       endif
       !
       !        Calculate velocities on last iteration
       !

       if (lastiter.or.m0.eq.1.0) then
          vsound(i) = -sqrt((te(i)+ti(i))/mb * econv/mconv)
          if (spts(i).ge.lastiters) then
             vsubs(i) = vb(i)
             vsupers(i) = vb2(i)
          else
             vsubs(i) = vb2(i)
             vsupers(i) = vb(i)
          endif
       endif
       !
       !        These functions are being called for the
       !        same value of s = send as is done in the
       !        solvstep subroutine - so they simply return the
       !        values of the function at the current s
       !        and ignore the other parameters.
       !
       !        These are used in power balance calculations ...
       !
       peiv(i) =  peiupdt(spts(i),ne(i),ne(i),t1e, t1i,t1e,t1i,fval)
       peiv2(i) = fval
       pradv(i) = pradupdt(spts(i),ne(i),ne(i),t1e,t1e)
       pcxv(i) = pcxupdt(spts(i),t1i,t1i)
       phelpiv(i) = phelpiupdt(spts(i),ne(i),ne(i),t1e,t1e)
       epowv(i) = estepow(spts(i))
       ipowv(i) = estipow(spts(i))
       
       if (actswnmom.eq.4) then
          smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
          ptmp = press(spts(i),te(i),ti(i))
          !            write (6,*) 'Smomtmp2:',smomtmp,ptmp
       elseif (actswnmom.eq.9.or.actswnmom.eq.10) then
          scxv(i) = scxupdt(spts(i),ne(i),ne(i),t1i,t1i)
       endif
       !
       !        Record electron and ion conduction components as well
       !        as electron and ion total power components.
       !
       convi(i) = cond(spts(i),t1i) + conv(spts(i),ne(i),t1i)
       conve(i) = cond(spts(i),t1e)

       condi(i) =  - (convi(i) + pcxv(i)-peiv(i) + ipowv(i) + estppion(spts(i)) +pais(spts(i)))

       conde(i) = -(conve(i)+pradv(i)+ phelpiv(i)+peiv(i) + epowv(i) + estppelec(spts(i)) +paes(spts(i)))
       !
       !
       ! slmod begin - new
       IF (imflag.EQ.1) THEN
          note(i) = ' i'
       ELSE
          note(i) = '  '
       ENDIF
       ! slmod end
       !
       !     This is the end do for the spts array
       !
    end do
    !
    !     Calculate Integrated Power Ratios - E,I,Total
    !
    !     Zero accumulators and results
    !
    int_powe = 0.0
    int_powi = 0.0
    int_conde= 0.0
    int_condi= 0.0

    do i = 1, npts

       int_powe = int_powe + (sbnd(i)-sbnd(i-1)) * (abs(conde(i)+ conve(i)))
       int_powi = int_powi + (sbnd(i)-sbnd(i-1)) * (abs(condi(i)+ convi(i)))

       int_conde= int_conde+ (sbnd(i)-sbnd(i-1)) * abs(conde(i))
       int_condi= int_condi+ (sbnd(i)-sbnd(i-1)) * abs(condi(i))
       !
       !         write (6,'(a,2i4,1x,1p,8e12.4)') 'PRAT:',ringnum,i, &
       !          abs(conve(i))+abs(conde(i)),conde(i),conve(i),&
       !          abs(convi(i))+abs(condi(i)),condi(i),convi(i),&
       !          sbnd(i),sbnd(i-1)
       !
    end do
    !
    !     Assign power ratios
    !
    int_powrat(1) = int_conde/int_powe
    int_powrat(2) = int_condi/int_powi
    int_powrat(3) = (int_conde+int_condi)/(int_powe+int_powi)
    !
    !      write(6,'(a,3(1x,g12.5))') 'PR:',int_powrat(1),int_powrat(2),&
    !                         int_powrat(3)
    !
    !
    !     Wrap up processing and print/plot the results
    !

    if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then

       write(6,*) 'Power Terms (QI,QE) (after): ',ringnum,nptscopy
       do ik = startn,nptscopy
          if (ik.lt.100.or.ik.eq.int(ik/real((real(nptscopy)/100.0)))*int(nptscopy/100)) then
             write(6,'(i4,20(1x,g13.6))') ik,sptscopy(ik),pcxv(ik),phelpiv(ik),conde(ik),conve(ik),&
                  pradv(ik),peiv(ik),estppelec(sptscopy(ik)),epowv(ik),ipowv(ik),paes(sptscopy(ik))
          endif
       end do

       write (6,*) '------'
    endif
    !
    !     Print out the model parameters
    !
    ! jdemod - only print this information if the print option is specified - could change this to a separate
    !          switch later.
    if (sol22_cprint.eq.3.or.sol22_cprint.eq.9) then 
       if (float((ir-irsep)/3).eq.(float(ir-irsep)/3.0)) then

          call echosolorg (spts,npts)

          if (imflag.eq.1.or.lastflag.eq.1) then
             call prbs
             call prs ('Caution: Imaginary roots encountered in n and v solutions')
             if (actswmach.eq.1.0.or.actswmach.eq.2.0) then
                call prs ('Target mach number was increased until imaginary roots vanished')
             endif
             call prbs
          endif

          call prbs
          call prs('Table of calculated SOL characteristics')
          call prbs
          write(comment,200)
          call prs(comment)
          do k = startn, npts
             write(comment,100) spts(k),te(k),ti(k),ne(k),vb(k)
             call prs (comment)
          end do
          !      
          !     Print out table of Viscosity estimates values
          !
          call prbs
          call prs('Table of calculated Viscosity values with Pressure')
          call prbs
          write(comment,300)
          call prs(comment)
          do k = startn, npts
             write(comment,100) spts(k),pir(k),pii(k),vgradn(k),act_press(k),exp_press(k)
             call prs (comment)
          end do
          !
          !     Print out tables of the velocity values
          !
          call prbs
          call prs('Tables of calculated SOL Velocity values')
          call prbs
          write(comment,400)
          call prs(comment)
          do i = startn, npts
             write(comment,100) spts(i),vsubs(i),vsupers(i),vsound(i),vb(i)
             call prs(comment)
          end do
          !
          !     End of IR=IRSEP if block
          !
       endif
    endif
    !
    !     Test graph variable ... if set then plot the results.
    !
    if (graph.eq.1.or.graphaux.eq.1.or.graphvel.eq.1) then
       CALL GPSTOP (100)
       CALL PAPER  (1)
       CALL HRDLIN(1)
       CALL HRDCHR(1)
    endif

    if (graph.eq.1) then
       call mkplot(spts,npts,te,ti,ne,vb,0)
       if (graphran.gt.0.0) then
          call mkplot(spts,npts,te,ti,ne,vb,1)
       endif
    endif

    if (graphaux.eq.1.and.(float((ir-irsep)/3).eq.(float(ir-irsep)/3.0))) then
       call mkauxplot(spts,npts,te,ti,ne,vb,phelpiv,peiv,peiv2,pcxv,pradv,0)
       if (graphran.gt.0.0) then
          call mkauxplot(spts,npts,te,ti,ne,vb,phelpiv,peiv,peiv2,pcxv,pradv,1)
       endif
    endif

    if (graphvel.eq.1) then
       call mkvelplot(spts,npts,te,ti,ne,vb,ne2,ga,0)
       if (graphran.gt.0.0) then
          call mkvelplot(spts,npts,te,ti,ne,vb,ne2,ga,1)
       endif
    endif

    if (graph.eq.1.or.graphaux.eq.1.or.graphvel.eq.1) then
       call grend
    endif
    !
    !     Assign final velocity to n0
    !
    n0 = nefinal
    !
    !     Post-process the contents of pradv (the integrated contribution
    !     of the radiation term) to obtain an estimate of the
    !     distributed radiation source.
    !
    !     jdemod - prad here likely needs pradv(i)-pradv(i-1)/(spts(i)-spts(i-1)) to get W/m3 out of it
    !            - need to check
    !
    do i = 1,npts

       if (i.lt.startn) then
          prad(i) = 0.0
       elseif (i.eq.1) then
          prad(i) = pradv(i)
       else
          prad(i) = pradv(i) - pradv(i-1)
       endif

    enddo
    ! slmod begin - new
2000 CONTINUE

    !...  Assign power channel data:       
    DO i = loopstart, npts
       cde(i) = conde(i)
       cdi(i) = condi(i)
       cve(i) = conve(i)
       cvi(i) = convi(i)
    ENDDO

    !...      
    CALL SOL22Output(loopstart,spts,npts,conde,condi,conve,convi,pcxv,peiv,phelpiv,pradv,te,ti,ne,vb,ga,act_press,pmloss,note)
    ! slmod end
    !
    !     Exit routine
    !
    return


100 format (6g13.5)
200 format (8x,'S',14x,'Te',13x,'Ti',13x,'Ne',13x,'Vb')
300 format (12x,'S',10x,'PiR',10x,'Pii',8x,'Vgrad',3x,'Act.Press.',3x,'Exp.Press.')
400 format (6x,'S',13x,'Vsub',12x,'Vsuper',11x,'Cs',14x,'Vused')

  end subroutine calcsol




  real*8 function vgradval(s,ind)
    use mod_solparams
    use mod_solcommon
    use mod_solswitch
    implicit none
    real*8 s
    integer ind
    !
    !     This function returns the current best estimate
    !     of the velocity gradient. This is used in conjunction
    !     with the vicosity solution options.
    !
    vgradval = vgrad

    return
  end function vgradval




  real*8 function find_ffric(ir,targid,actmomlen,actmomlam)
    !use mod_params
    use mod_solswitch
    use mod_solparams
    use mod_solcommon
    implicit none
    !     include 'params'
    !     include 'solparams'
    !     include 'solcommon'
    integer ir,targid

    !     Assign the actual FFRIC values for each ring from the data contained in
    !     the inputs FFRIC and EXTFFRIC. This option allows for a different
    !     amount of momentum loss to be specified for every half flux tube.

    real*8 :: actmomlen,actmomlam


    ! The parameter read in is the length for source option 1 and the exponential decay for option 2

    !     Negative input values use the default values
    !
    integer in

    find_ffric = ffric

    actmomlen = lenmom
    actmomlam = lammom

    !write(0,*) 'switch(swnmom):',switch(swnmom),actswnmom,ir,targid,n_extffric,ir

    if (n_extffric.eq.0) return        
    !
    do in = 1,n_extffric
       !
       if (ir.eq.int(extffric(in,1))) then 
          if (targid.eq.2) then 
             if (extffric(in,2).gt.0.0) then 
                find_ffric = extffric(in,2)
             endif
             if (extffric(in,3).gt.0.0) then 
                if (switch(swnmom).eq.2.0) then 
                   actmomlam = extffric(in,3)
                else
                   actmomlen = extffric(in,3)
                endif
             endif
          elseif(targid.eq.1) then 
             if (extffric(in,4).gt.0.0) then 
                find_ffric = extffric(in,4)
             endif
             if (extffric(in,5).gt.0.0) then 
                if (switch(swnmom).eq.2.0) then 
                   actmomlam = extffric(in,5)
                else
                   actmomlen = extffric(in,5)
                endif
             endif
          endif
          !write(0,*) 'switch(swnmom):',switch(swnmom),actswnmom,ir,targid,find_ffric,actmomlen,actmomlam

          return

       endif

    end do

    return

  end function find_ffric


  subroutine assign_radiation_parameters(ir,targid)
    use mod_solcommon
    implicit none
    !
    !     jdemod
    !
    !     Base radiation parameters are lenr, lamr, frr which
    !     define the basic radiation source. However, the amount
    !     of radiation needed can vary by flux tube to match the
    !     experimental temperature profiles along the field lines

    !     This code allows for ring by ring customization of
    !     the radiation.
    !
    !     Negative extended source values use the default values
    !     If data is not specified for a ring it also uses the default values
    !
    integer :: ir,targid


    !      lenr = lenri

    !     Adjust Radiation Source length

    integer :: in

    lenr = min(lenri,halfringlen)
    lamr = lamri
    frr  =  frri


    if (n_extradsrc.eq.0) return
    !
    do in = 1,n_extradsrc
       !
       if (ir.eq.int(extradsrc(in,1))) then 
          if (targid.eq.2) then 
             if (extradsrc(in,2).gt.0.0) then 
                lenr = extradsrc(in,2)
             endif
             if (extradsrc(in,3).gt.0.0) then 
                lamr = extradsrc(in,3)
             endif
             if (extradsrc(in,4).gt.0.0) then 
                frr = extradsrc(in,4)
             endif
          elseif(targid.eq.1) then 
             if (extradsrc(in,5).gt.0.0) then 
                lenr = extradsrc(in,5)
             endif
             if (extradsrc(in,6).gt.0.0) then 
                lamr = extradsrc(in,6)
             endif
             if (extradsrc(in,7).gt.0.0) then 
                frr = extradsrc(in,7)
             endif
          endif
          lenr = min(lenr,halfringlen)

          return

          !
       endif
       !
    end do
    return
  end subroutine assign_radiation_parameters


end module mod_sol22
