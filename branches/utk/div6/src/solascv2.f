c
c
c
c     program calcsol
c
      subroutine calcsol (spts,npts,errcode,serr,
     >                    te,ti,ne,vb,exp_press,act_press,
     >                    prad,ir,irsep,
     >                    int_powrat,cprint)
      !use sol22_input
      use sol22_debug
      implicit none
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      include 'solrk'
c slmod begin - new
c...Try to eliminate:      
      INCLUDE 'params'
      INCLUDE 'slcom'

      COMMON /POWERFLOW/ cve        ,cvi        ,cde        ,cdi
      REAL               cve(MXSPTS),cvi(MXSPTS),cde(MXSPTS),cdi(MXSPTS)
      
      INTEGER     eflag    
      REAL*8      pmloss(MXSPTS)
      CHARACTER*2 note(MXSPTS)
c slmod end   
c
      integer errcode,npts,ir,irsep,cprint
      real*8 serr
      real*8 spts(mxspts)
      real*8 te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts),
     >        exp_press(mxspts),act_press(mxspts),
     >        prad(mxspts)
      real*8 int_powrat(3)
c
c     CALCSOL:
c
c     This program employs a simple Runge-Kutta technique and several
c     approximations
c     to calculate the temperature, density and velocity distribution
c     along the scrape off layer.
c
c     Modified to use adaptive step-size control and Cash-Karp
c     parameters for embedded Runge-Kutte Method
c
c     David Elder, Aug 3, 1994
c
c
      real*8 sinit,send
c
c      real*8 h,s,ntmp
c
      integer    i,k,flag,imflag,lastflag,vcount
      integer    loopstart
      integer    ik
      real*8    n,cs
      real*8    t1e,t1i,v1
      real*8    fval
      real*8    fegrad,figrad,newn,gamma,peiupdt,pintupdt
      real*8    phelpiupdt,pcxupdt,smomtmp
      real*8    nimag,pmomloss,ptmp,pradupdt,scx_s,scxtmp
      real*8    dm0,iters,vgradval,press,scxupdt
      real*8    estppion,estppelec
      external   fegrad,figrad,newn,gamma,peiupdt,phelpiupdt,pcxupdt
      external   vgradval,pintupdt,pmomloss,press,pradupdt
      external   scx_s,scxupdt
      external  estppion,estppelec
      real*8,external :: srci,srcf

      real*8    ga(mxspts),vb2(mxspts)
      real*8    ne2(mxspts),vsupers(mxspts),vsubs(mxspts)
      real*8    vsound(mxspts),scxv(mxspts)
      real*8    peiv(mxspts),pcxv(mxspts),phelpiv(mxspts)
      real*8    peiv2(mxspts),pradv(mxspts)
      real*8    pir(mxspts),pii(mxspts),vgradn(mxspts)
c
      real*8    condi(mxspts),conde(mxspts)
      real*8    convi(mxspts),conve(mxspts)
      real*8    int_powe,int_powi,int_conde,int_condi
c
c     Function declarations
c
      real*8 cond,conv,paes,pais
      external cond,conv,paes,pais
c
c      real*8    tpress(mxspts)
c      real*8    tmpgam ,errmax, errti,errte
c      real*8    newv1,errvb,vrat
      real*8    tauii
c      real*8    hnew
c      real*8    vsep,vsub,vsup,lastvsup,v2
c
      integer      exitcond,negerrflag
      character*80 comment
c
c      logical      imag
c
c     logical      errneg
c
c     define arrays for k-values
c
c      real*8 ke(6),ki(6)
c
c     Initialization
c
c slmod begin - new
      DATA t1i,t1e /0.0,0.0/
c
c
c     jdemod - setting simag1 to the ring length should be big enough
c
c      simag1 = HI
      simag1 = ringlen
c
c     jdemod - halflen should probably be halfringlen since halflen is 
c              never calculated or assigned
c
c      simag2 = halflen
      simag2 = halfringlen
      ierror = MAXNKS
      eflag  = 0

      CALL DZero(conve ,MXSPTS)
      CALL DZero(convi ,MXSPTS)
      CALL DZero(conde ,MXSPTS)
      CALL DZero(condi ,MXSPTS)
      CALL DZero(pmloss,MXSPTS)
c slmod end
c
c     jdemod - supplemented by sol22_debug module functionality
      debug_s22 = .false.

c
c     Set up the value of startn for this ring ...
c
      if (actswe2d.ne.0.0) then
c
         startn = ike2d_start
c
      else
c
         startn = 1
c
      endif
c
      write(6,*) 'CALCSOL: Errlevel = ',actswerror
c
c     Set the minimum allowed temperature for the solver.
c
      if (tmin.gt.0.0) then
         tfloor = tmin
      elseif (tmin.lt.0.0) then
         tfloor = min(abs(tmin),min(te0,ti0))
      else
         tfloor = 0.0
      endif
c
      timax = 0.0
      temax = 0.0
c
      call loadparms
c
c     Record the target density.
c
      netarg = n0
      nefinal= n0
c
      n1center= n1
      n1final = n1
c
      errcode = 0
      serr = 0.0
c
      flag = 0
      imflag = 0
      lastflag = 0
c      imag = .false.
      lastiter = .true.
      lastiters = 0.0
      founds = .true.
c
c     Recalculate the initm0 value (the value the solver starts at)
c     for the case of Edge2D compatibility options where the solver
c     starts at the middle of the first cell and the velocity
c     there and thus the Mach number are specified from the
c     Edge2D data.
c
c
      if (actswmach.eq.3.0) then
c
         cs = -sqrt((te0+ti0)/mb * econv/mconv)
c
         if (cs.ne.0.0) then
            initm0 = abs(v0/cs)
         else
            initm0 = 1.0
         endif
c
         write (6,'(a,6g12.4)') 'Initm0:', te0,ti0,mb,v0,cs,initm0
c
      else
c
         initm0 = 1.0
c
      endif
c
      e2dm0 = 1.0
c
      if (actswe2d.ne.0.0) then
c
         cs = -sqrt((te1+ti1)/mb * econv/mconv)
c
         if (actswe2d.eq.1.0) then
            e2dm0 = vpe2d / cs
         elseif (actswe2d.eq.2.0) then
            e2dm0 = v1e2d / cs
         elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.
     >           actswe2d.eq.9.0) then
            e2dm0 = vpg / cs
         endif
c
      endif
c
c     Copy Spts for use in some routines - since it involves a lot
c     code rewriting to move spts to a common block or to pass it
c     to the routine srci - where it is needed.
c
c     Need to set values in ALL arrays for the npts+1 position
c     Just copy over the values from npts.
c     So ...
c     Leave out calculating all the way to the mid-point.
c
c     spts(npts+1) is set to S-value at the boundary of the cell 
c     in the calcsol_interface routine - this is the same as the 
c     halfringlen value but is NOT the same as ringlen/2.0
c
      nptscopy = npts
      do i = 1,nptscopy+1
         sptscopy(i) = spts(i)
      end do
c
c     Read in the parameters for the ring from a file
c
c      call readsol (spts,npts)
c
c     When operating in the mode that iteratively finds the lowest Mach
c     number at the target that does not result in imaginary values
c     of density or velocity when solving the equations ... the code
c     loops back to this point ... and recalculates from here.
c
c
c
      stopimag = .false.
      miter = 0
c
      m0 = initm0
      m1 = e2dm0
c
      if (actswmach.eq.1.0.or.actswmach.eq.2.0) then
         dm0 = deltam0
         origm0 = m0
         lastm0 = m0
         stopimag = .true.
         lastiter = .false.
      endif
c

1000  continue

c
      miter = miter + 1
c
      if ((actswmach.eq.1.0.or.actswmach.eq.2.0)
     >    .and.miter.gt.1) then
         lastm0 = m0
         m0 = m0 + dm0
c
c        Change initial target density if Mach solver option
c        is 2.0
c
         if (actswmach.eq.2.0) then
c
            n0 = netarg * initm0 / m0
            nefinal = n0
c
            if (actswe2d.ne.0.0) then
              n1 = n1center * initm0 / m0
              n1final = n1
            endif
c
         endif
c
         if (miter.gt.maxiter) then
             write(6,*) 'Error: Program exceeds maximum'//
     >                 ' iterations calculating m0'
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
c
      if ((actswmach.eq.1.0.or.actswmach.eq.2.0)
     >      .and.(m0.ne.origm0)) then
         founds = .false.
      else
         founds = .true.
      endif
c
c     Set various values depending on switches ... for example the gamma
c     factor for power to the plates is different depending on options.
c     Pre-calculate the integral of the numerical ionization source - if
c     there is one.
c
      call initval
c
 1500 continue
c
c     Set up the initial values
c
      if (spts(1).eq.0.0) then
         te(1) = te0
         ti(1) = ti0
         ne(1) = n0
         ne2(1) = m0**2 * n0
         vb(1) = v0
         ga(1) = n0 * v0
         pir(1) = -0.4 * n0 * econv * ti0
         pii(1) = -0.4 * n0 * econv * ti0
c
         exp_press(1) = pinf
         act_press(1) = pinf
c
         vgradn(1) = vgradval(spts(1),0)
         loopstart = startn + 1
      elseif (actswe2d.ne.0.0) then

         te(ike2d_start) = te1
         ti(ike2d_start) = ti1
         ne(ike2d_start) = n1
         ne2(ike2d_start) = m1**2 * n1
c
c        Different E2D options
c
         if (actswe2d.eq.1.0) then
            vb(ike2d_start) = vpe2d
            ga(ike2d_start) = n1 * vpe2d
            v1 = vpe2d
         elseif (actswe2d.eq.2.0) then
            vb(ike2d_start) = v1e2d
            ga(ike2d_start) = n1 * v1e2d
            v1 = v1e2d
         elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.
     >           actswe2d.eq.9.0) then
            vb(ike2d_start) = vpg
            ga(ike2d_start) = n1 * vpg
            v1 = vpg
         endif
c
         pir(ike2d_start) = -0.4 * n1 * econv * ti1
         pii(ike2d_start) = -0.4 * n1 * econv * ti1
c
         exp_press(ike2d_start) = pinf
         act_press(ike2d_start) = pinf
c
         vgradn(ike2d_start) = vgradval(spts(ike2d_start),0)
c
         loopstart = startn+1
c
      else
         loopstart = startn
      endif
c
c     More initialization
c
      vcount = 0
c
      vsepmin = abs(v0)
      ssepmin = 0.0
      nlast = n0
      slast = 0.0
c
c     If the options are set - initialize the UPDATE/ESTIMATE
c     integral routines for S=Soffset..
c
      if (actswnmom.eq.4) then
         if (actswe2d.eq.0) then
            smomtmp = pmomloss(soffset,1,v0,te0,ti0)
         else
            smomtmp = pmomloss(soffset,1,v1,te1,ti1)
         endif
         ptmp = press(soffset,te0,ti0)
c
         write (6,*) 'Smomtmp1:',smomtmp,ptmp
c
      elseif (actswnmom.eq.9.or.actswnmom.eq.10) then
c
c     WF'96: Initialize CX momentum loss routine
c
         if (actswe2d.eq.0) then
            scxtmp = scxupdt(soffset,n0,n0,ti0,ti0)
         else
            scxtmp = scxupdt(soffset,n1,n1,ti1,ti1)
         endif
         ptmp = press(soffset,te0,ti0)
c
c         write (6,*) 'Smomtmp1:',scxtmp,ptmp
c
      endif
c
c     The actual values of the arguments (except soffset) are
c     irrelevant. These lines simply garantee the initialization
c     of the integration routines for this 1/2 ring.
c
      ptmp =  peiupdt(soffset,ne(1),ne(1),
     >              t1e, t1i,
     >              t1e,t1i,fval)
      ptmp =  pradupdt(soffset,ne(1),ne(1),
     >              t1e,t1e)
      ptmp = pcxupdt(soffset,t1i,t1i)
      ptmp = phelpiupdt(soffset,
     >             ne(1),ne(1),t1e,t1e)
c slmod begin - new
      simag   = 0.0D0
      snegerr = 0.0D0
c slmod end
c
c     Loop through the S values
c
      do i = loopstart,npts
c
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
            call save_s22_data(dble(i),sinit,n,t1e,t1i,v1,gamma(sinit),
     >                         srci(sinit),srcf(sinit))
         endif
c
c        Inner loop for the steps between each S-value.
c
         timax = max(timax,t1i)
         temax = max(temax,t1e)
c
         if (debug_s22)
     >      write(6,'(a,10(1x,g12.5))') 'S1:',sinit,send,
     >            t1e,t1i,n


         call solvstep(sinit,send,t1e,t1i,n,exitcond,imflag,
     >                 negerrflag,vcount)
c slmod begin - new
c         CALL NoName01(imflag, i,  loopstart,spts,npts,
c     .                 exitcond,  note,te,ti,ne,ne2,ga,vb,vb2,
c     .                 pir,pii,vgradn,act_press,exp_press,
c     .                 vsound,vsubs,vsupers,peiv,peiv2,pradv,pcxv,
c     .                 phelpiv,scxv,convi,conve,condi,conde,errcode,
c     .                 *2000)
c         
c         WRITE(0,*) 'Continuing nicely'
c slmod end
c
         if (imflag.eq.1) then
c
c           An imaginary value was encountered - this value of
c           the error code will be overwritten if a
c           more serious errcode comes along.
c
            errcode = 1
            serr = simag

c     
c           jdemod - copy the code setting simag1 and simag2 from NoName
c
            IF (simag1.EQ.ringlen) THEN
               simag1 = simag - soffset
               simag2 = simag - soffset
            ELSE
               simag2 = simag - soffset
            ENDIF

         endif
c
c        Exit conditions 1 and 2 do not generate error codes
c        because they are a normal part of the mach number
c        solver - if imaginary results lead to the situation
c        but an eventual solution is found - then errcode will
c        have been set to 1 in the above statement.
c
         if (exitcond.eq.1) then
c
c           Iterate at larger mach number
c

            goto 1000
c
         elseif (exitcond.eq.2) then
c
c           Iterate at smaller mach number
c
            if (dm0.le.m0res) then
c
c              If one has hit the mach number resolution
c              limit then run through the entire
c              ring solving one last time.
c
               lastiter = .true.
               founds = .false.
c
c              Set M0 to the previous smaller value.
c
c               m0 = lastm0
c
               write(6,*) 'Last iteration:',lastiters
               goto 1500
            else
c
c              Rest to last non-working m0 and reduce the
c              delta m0 value by a factor of 10. - then
c              iterate.
c
               dm0 = dm0/10.0
               if (dm0.lt.m0res) dm0 = m0res
               m0 = lastm0
               goto 1000
            endif
         elseif (exitcond.eq.3) then
c
c           Negative temperature encountered - unexpected exit.
c
            errcode = exitcond
            serr    = snegerr
            write(6,*) 'Exitcond = 3:',snegerr
c
            return
c
         elseif (exitcond.eq.4) then
c
c           Step size too small - unexpected exit.
c           Possibly caused in response to T -> 0
c
c           Step size has hit lower limit and no other
c           actions are available to try to push
c           solution forward.
c
            errcode = exitcond
            serr    = snegerr
            write(6,*) 'Exitcond = 4:',snegerr
c
c slmod begin - new
            GOTO 2000
c
c            return
c slmod end
c
         elseif (exitcond.eq.5) then
c
c           Temperature along ring dropped to less
c           than dropfrac*tmax OR temperature hit
c           minimum level limit for ring.
c
            errcode = exitcond
            serr    = snegerr
            write(6,*) 'Exitcond = 5:',snegerr
c
c slmod begin - new
            GOTO 2000
c
c            return
c slmod end
c
         elseif (exitcond.eq.6) then
c
c           Negative density encountered in solver
c           Code will exit - if the error corrector is OFF -
c           it will be turned ON for this half-ring ONLY -
c           an error message will be recorded to thsi effect.
c           If the error corrector is ON - it will proceed
c           normally.
c
            errcode = exitcond
            serr    = snegerr
            write(6,*) 'Exitcond = 6:',snegerr
c
c slmod begin - new
            GOTO 2000
c
c            return
c slmod end
c
         elseif (exitcond.eq.7) then
c
c           NaNQ encountered in solver - usually indicative of
c           a negative temperature condition.
c
c           Code will exit - if the error corrector is OFF -
c           it will be turned ON for this half-ring ONLY -
c           an error message will be recorded to this effect.
c           If the error corrector is ON - it will proceed
c           normally.
c
            errcode = exitcond
            serr    = snegerr
            write(6,*) 'Exitcond = 7:',snegerr
c
c slmod begin - new
            GOTO 2000
c
c            return
c slmod end
c
         endif
c
c       Step was successfully solved - run through post processing
c       and record all the quantities for the cell.
c
c       Tabulate the n,te,ti values for the grid point.
c
         te(i) = t1e
         ti(i) = t1i
c
         ne(i) = newn(spts(i),t1e,t1i,nimag,flag)
         ne2(i) = nimag
c
         ga(i) = gamma(spts(i))
         vb(i) = ga(i) / ne(i)
         if (ga(i).eq.0.0) then
            vb2(i) = 0.0
         else
            vb2(i) = ga(i) / ne2(i)
         endif
c
c
         if (cprint.eq.3.or.cprint.eq.9) then
            write(6,'(a,i3,9(1x,g12.5),2i2)') 'Step:',
     >              i,m0,send,te(i),ti(i),
     >              ne(i),vb(i),ga(i),ne(i)*vb(i),ionsrc(i),flag,
     >              negerrflag
         endif
c
c
c        Calculate viscosity estimate and limit terms
c
         pir(i) = -0.4 * ne(i) * econv * ti(i)
         tauii = 2.5d0 * 2.09d13*(ti(i))**1.5* sqrt(mb)
     >                 / (ne(i)  * 15.0d0)
c
         pii(i) = - 4.0d0/9.0d0 * ne(i) * econv * ti(i) *
     >             tauii * vgradval(spts(i),0)
         vgradn(i) = vgradval(spts(i),0)
c
c        Calculate current pressure
c
         act_press(i) = ne(i) * econv * (te(i) + ti(i))
     >               + mb * mconv * ne(i) * vb(i)**2
c
c        Calculate Expected pressure
c
c
c         make sure pmomloss set up correctly for call from pressure
c
          if (actswnmom.eq.4) then
             if (actswe2d.eq.0) then
                smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
             else
                smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
             endif
         endif

         exp_press(i) = press(spts(i),te(i),ti(i))
c
         if (cprint.eq.3.or.cprint.eq.9) then

            write (6,'(a,6g13.5)') 'Pres:',pinf,
     >          exp_press(i)-pinf-padd,
     >          exp_press(i)-pinf,exp_press(i),act_press(i),padd
         endif
c
c        Calculate velocities on last iteration
c

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
c
c        These functions are being called for the
c        same value of s = send as is done in the
c        solvstep subroutine - so they simply return the
c        values of the function at the current s
c        and ignore the other parameters.
c
c        These are used in power balance calculations ...
c
         peiv(i) =  peiupdt(spts(i),ne(i),ne(i),
     >              t1e, t1i,t1e,t1i,fval)
         peiv2(i) = fval
         pradv(i) = pradupdt(spts(i),ne(i),ne(i),t1e,t1e)
         pcxv(i) = pcxupdt(spts(i),t1i,t1i)
         phelpiv(i) = phelpiupdt(spts(i),
     >             ne(i),ne(i),t1e,t1e)
c
         if (actswnmom.eq.4) then
            smomtmp = pmomloss(spts(i),1,vb(i),te(i),ti(i))
            ptmp = press(spts(i),te(i),ti(i))
c            write (6,*) 'Smomtmp2:',smomtmp,ptmp
         elseif (actswnmom.eq.9.or.actswnmom.eq.10) then
            scxv(i) = scxupdt(spts(i),ne(i),ne(i),
     >                  t1i,t1i)
         endif
c
c        Record electron and ion conduction components as well
c        as electron and ion total power components.
c
         convi(i) = cond(spts(i),t1i) + conv(spts(i),ne(i),t1i)
         conve(i) = cond(spts(i),t1e)
c
         condi(i) =  - (convi(i)
     >               + pcxv(i)-peiv(i)
     >               + estppion(spts(i)) +pais(spts(i)))
c
         conde(i) = -(conve(i)+pradv(i)+
     >             phelpiv(i)+peiv(i)
     >            + estppelec(spts(i)) +paes(spts(i)))
c
c
c slmod begin - new
        IF (imflag.EQ.1) THEN
          note(i) = ' i'
        ELSE
          note(i) = '  '
        ENDIF
c slmod end
c
c     This is the end do for the spts array
c
      end do
c
c     Calculate Integrated Power Ratios - E,I,Total
c
c     Zero accumulators and results
c
      int_powe = 0.0
      int_powi = 0.0
      int_conde= 0.0
      int_condi= 0.0
c

      do i = 1, npts
c
         int_powe = int_powe + (sbnd(i)-sbnd(i-1)) *
     >                         (abs(conde(i)+ conve(i)))
         int_powi = int_powi + (sbnd(i)-sbnd(i-1)) *
     >                         (abs(condi(i)+ convi(i)))
c
         int_conde= int_conde+ (sbnd(i)-sbnd(i-1)) * abs(conde(i))
         int_condi= int_condi+ (sbnd(i)-sbnd(i-1)) * abs(condi(i))
c
c         write (6,'(a,2i4,1x,1p,8e12.4)') 'PRAT:',ringnum,i,
c     >     abs(conve(i))+abs(conde(i)),conde(i),conve(i),
c     >     abs(convi(i))+abs(condi(i)),condi(i),convi(i),
c     >     sbnd(i),sbnd(i-1)
c
      end do
c
c     Assign power ratios
c
      int_powrat(1) = int_conde/int_powe
      int_powrat(2) = int_condi/int_powi
      int_powrat(3) = (int_conde+int_condi)/(int_powe+int_powi)
c
c      write(6,'(a,3(1x,g12.5))') 'PR:',int_powrat(1),int_powrat(2),
c     >                    int_powrat(3)
c
c
c     Wrap up processing and print/plot the results
c

      write(6,*) 'Power Terms (QI,QE) (after): ',ringnum,nptscopy
      do ik = startn,nptscopy
         write(6,'(i4,20(1x,g13.6))') ik,sptscopy(ik),
     >        pcxv(ik),phelpiv(ik),conde(ik),conve(ik),
     >        pradv(ik),peiv(ik),
     >        estppelec(sptscopy(ik)),paes(sptscopy(ik))
      end do
c
      write (6,*) '------'
c
c     Print out the model parameters
c
      if (float((ir-irsep)/3).eq.(float(ir-irsep)/3.0)) then
c
         call echosolorg (spts,npts)
c
      if (imflag.eq.1.or.lastflag.eq.1) then
         call prbs
         call prs ('Caution: Imaginary roots encountered'//
     >             ' in n and v solutions')
         if (actswmach.eq.1.0.or.actswmach.eq.2.0) then

            call prs ('Target mach number was increased until'//
     >             ' imaginary roots vanished')
         endif
         call prbs
      endif
c
      call prbs
      call prs('Table of calculated SOL characteristics')
      call prbs
      write(comment,200)
      call prs(comment)
      do k = startn, npts
         write(comment,100) spts(k),te(k),ti(k),ne(k),vb(k)
         call prs (comment)
      end do
c
c     Print out table of Viscosity estimates values
c
      call prbs
      call prs('Table of calculated Viscosity values with'//
     >         ' Pressure')
      call prbs
      write(comment,300)
      call prs(comment)
      do k = startn, npts
         write(comment,100) spts(k),pir(k),pii(k),vgradn(k),
     >                      act_press(k),exp_press(k)
         call prs (comment)
      end do
c
c     Print out tables of the velocity values
c
      call prbs
      call prs('Tables of calculated SOL Velocity values')
      call prbs
      write(comment,400)
      call prs(comment)
      do i = startn, npts
         write(comment,100) spts(i),vsubs(i),vsupers(i),
     >                      vsound(i),vb(i)
         call prs(comment)
      end do
c
c     End of IR=IRSEP if block
c
      endif
c
c     Test graph variable ... if set then plot the results.
c
      if (graph.eq.1.or.graphaux.eq.1.or.graphvel.eq.1) then
         CALL GPSTOP (100)
         CALL PAPER  (1)
         CALL HRDLIN(1)
         CALL HRDCHR(1)
      endif
c
      if (graph.eq.1) then
         call mkplot(spts,npts,te,ti,ne,vb,0)
         if (graphran.gt.0.0) then
            call mkplot(spts,npts,te,ti,ne,vb,1)
         endif
      endif
c
      if (graphaux.eq.1.and.
     >    (float((ir-irsep)/3).eq.(float(ir-irsep)/3.0))
     >    ) then
         call mkauxplot(spts,npts,te,ti,ne,vb,
     >                  phelpiv,peiv,peiv2,pcxv,pradv,0)
         if (graphran.gt.0.0) then
            call mkauxplot(spts,npts,te,ti,ne,vb,
     >                  phelpiv,peiv,peiv2,pcxv,pradv,1)
         endif
      endif
c
      if (graphvel.eq.1) then
         call mkvelplot(spts,npts,te,ti,ne,vb,ne2,ga,0)
         if (graphran.gt.0.0) then
            call mkvelplot(spts,npts,te,ti,ne,vb,ne2,ga,1)
         endif
      endif
c
      if (graph.eq.1.or.graphaux.eq.1.or.graphvel.eq.1) then
         call grend
      endif
c
c     Assign final velocity to n0
c
      n0 = nefinal
c
c     Post-process the contents of pradv (the integrated contribution
c     of the radiation term) to obtain an estimate of the
c     distributed radiation source.
c
      do i = 1,npts

         if (i.lt.startn) then
            prad(i) = 0.0
         elseif (i.eq.1) then
            prad(i) = pradv(i)
         else
            prad(i) = pradv(i) - pradv(i-1)
         endif

      enddo
c slmod begin - new
2000  CONTINUE

c...  Assign power channel data:       
      DO i = loopstart, npts
        cde(i) = conde(i)
        cdi(i) = condi(i)
        cve(i) = conve(i)
        cvi(i) = convi(i)
      ENDDO

c...      
      CALL SOL22Output(loopstart,spts,npts,conde,condi,conve,convi,
     .                 pcxv,peiv,phelpiv,pradv,te,ti,ne,vb,ga,
     .                 act_press,pmloss,note)
c slmod end
c
c     Exit routine
c
      return
c

 100  format (6g13.5)
 200  format (8x,'S',14x,'Te',13x,'Ti',13x,'Ne',13x,
     >        'Vb')
 300  format (12x,'S',10x,'PiR',10x,'Pii',8x,'Vgrad',
     >        3x,'Act.Press.',3x,'Exp.Press.')
 400  format (6x,'S',13x,'Vsub',12x,'Vsuper',11x,'Cs',14x,
     >        'Vused')
c

      end
c
c
c
      subroutine solvstep(sinit,send,t1e,t1i,n,exitcond,imflag,
     >                    negerrflag,vcount)
      !use sol22_input
      use sol22_debug
      implicit none
      real*8 sinit,send,t1i,t1e,n
      integer exitcond,imflag,vcount,negerrflag
      include 'solparams'
      include 'solcommon'
      include 'solswitch'
c
c     This routine uses an RK driver routine to solve from
c     sinit to send - if it finds that this can't be done it
c     returns a non-zero exitcondition. A zero exit condition
c     means that it was successful in stepping from sinit to send.
c
c     Local variables
c
      integer flag,imflag2
      real*8 s,h,newt1e,newt1i,errte,errti
      integer ierr
      logical errneg,imag
      real*8 newn,pcxupdt,peiupdt,phelpiupdt,gamma,pintupdt
      real*8 pmomloss,fval,pradupdt,scxupdt
      external newn,pcxupdt,peiupdt,phelpiupdt,gamma
      external pintupdt,pmomloss,pradupdt,scxupdt
      real*8,external :: srci,srcf

      real*8 ptmp
      real*8 tmpgam1,tmpgam2,vsep,vtmp
      real*8 v1,v2,vsub,vsup
      real*8 smomtmp,scxtmp
      real*8 errmax,nimag,hnew,stmp,tmpgam
c
      real*8 hneg
c
      real*8 ntmp,tmpnimag,fegrad,figrad,gtmp
      real*8 srtn
      integer tmpflag,iter
      external fegrad,figrad
c
c     Set up the separation factor - if s exceeds this
c     value * ssepmin - then one assumes that there
c     is an acceptable solution for the entire ring at
c     this Mach number..
c
      real sepfact
      parameter (sepfact= 1.2)
c
c
      iter = 1
      imflag = 0
      negerrflag = 0
      imflag2 = 0
      exitcond = 0
      s = sinit
      slast = sinit
c
      nlast = n
c
c     Set the value of hlim - the minimum allowable
c     step size in the equation solver. Based on the
c     separation of the two points and the initial number
c     of steps - keeping in mind that the default step
c     size for default imaginary cases is 10 times
c     smaller than the initial step-size.
c
      hlim = (send - sinit) / dble(ndiv) / 10000.0
      hneg = (send - sinit) / dble(ndiv) / 500.0
      himag = (send - sinit) / dble(ndiv) / 50.0
c
c      write (6,*) 'Hlim:',hlim,himag,sinit,send,h
c
c     starting h-value
c
      h = (send-sinit)/ dble(ndiv)
c

 2000 continue
c
      imflag2 = 0
c
c
c------------------------------------------------------------------------
c
c     Try to perform an RK step of size h
c
c
c
      call rkstep(s,t1e,t1i,newt1e,newt1i,errte,errti,h,m0,
     >            ierr)
c
c     Save new S value returned by rkstep - used to decide 
c     whether to save debug information. 
c     
      srtn = s + h
c
c      write(6,*) 'RKSTEP:',s,newt1i,newt1e,h,m0,ierr
c
c------------------------------------------------------------------------
c
c     Check error conditions from stepper
c
      if (ierr.eq.1) then
c
c         Negative temperature - adjust h smaller and iterate
c
          negerrflag = 1
c
c          write (6,'(a,5(1x,g12.5))')
c     >         'Error!: Negative TEMP',h,hlim,himag,s,m0
c          write (6,*) 'T:',t1e,newt1e,t1i,newt1i
c
c         ntmp = newn(s,t1e,t1i,tmpnimag,tmpflag)
c          write (6,*) 'G:',tmpflag,gamma(s),ntmp,
c               fegrad (s,t1e,t1i,ntmp),
c               figrad (s,t1e,t1i,ntmp)
c
         if (h.le.0.01*hneg) then
            if (lastiter) then
               write (6,*) 'ERROR: Negative temperature:'//
     >                     'h too small on last iteration'
               write (6,*) 'Iteration Stopping in SOLASCV:SOLVSTEP',
     >                      iter
               write (6,*) 'Possibly caused by'
     >                     //' Physically Inconsistent'
     >                     //' Target Conditions'
c
              write (6,'(a,5(1x,g13.5))')
     >            'Error!: Negative TEMP',h,hlim,himag,s,m0
              write (6,*) 'T:',t1e,newt1e,t1i,newt1i
              ntmp = newn(s,t1e,t1i,tmpnimag,tmpflag)
              gtmp = gamma(s)
              write (6,'(a,i4,4(1x,g13.5))') 'G:',tmpflag,gtmp,
     >            ntmp,
     >            fegrad (s,t1e,t1i,ntmp),
     >            figrad (s,t1e,t1i,ntmp)
c
               exitcond = 3
               snegerr = s
               return
            else
               write (6,*) 'ERROR: Negative temperature:'//
     >                     'h too small on last iteration'
               exitcond = 3
               snegerr = s
               return
            endif
         endif
         h = h / 10.0
         goto 2000
c
c
c     Check for imaginary solutions.
c
      elseif (ierr.eq.2) then
c
         simag = s
c
         if (simag.eq.0.0) then
              write (6,'(a,8(1x, g12.5))')
     >                         'imag:',
     >                         h,m0,
     >                         hlim,s,send,
     >                         sinit,newt1e,newt1i
         endif
c
         if (stopimag.and.(.not.lastiter)) then
c
            if (h.le.hlim) then
               exitcond = 1
               return
            elseif (h.gt.hlim) then
               h = h / 10.0
               goto 2000
            endif
c
         else
c
            if (h.le.hlim.and.(s+h).ne.send) then
               imflag = 1
               imflag2 = 1
               write (6,*) 'ERROR:'//
     >                     ' Imaginary encountered:'//
     >                     ' h too small:',ringnum,
     >                     h,m0,lastiter,lastiters,
     >                     stopimag,ierr,hlim,s,send,
     >                     sinit
               write (6,*) 'Program will NOT stop'
c               stop
c
               exitcond = 4
               return
c
            elseif (h.gt.himag) then
               h = h /10.0
               goto 2000
            else
               imflag = 1
               imflag2 = 1
            endif
         endif
c
c     Check for negative N
c
      elseif (ierr.eq.3) then
c
c        Negative N found - iterate smaller h OR set condition and exit
c
         if (h.le.hlim) then
            exitcond = 6
            snegerr  = s
            return
         elseif (h.gt.hlim) then
            h = h / 10.0
            goto 2000
         endif
c
c     Check for NaNQ
c
      elseif (ierr.eq.4) then
c
c        NaNQ found - iterate smaller h OR set condition and exit
c
         if (h.le.hlim) then
            exitcond = 7
            snegerr  = s
            return
         elseif (h.gt.hlim) then
            h = h / 10.0
            goto 2000
         endif
c
      endif
c
c
c------------------------------------------------------------------------
c
c     Need to test if h is too small here because if the
c     code is not numerically stable then a condition could
c     conceivably be generated where the h value just gets
c     smaller and smaller without encountering imaginary
c     numbers.
c
      if (h.le.hlim) then
c
         if (lastiter) then
            exitcond = 4
            snegerr = s
            write (6,*) 'ERROR: delta h less than minimum :',
     >                  h,hlim
            write (6,*) 'Possibly caused by Negative T: ',ierr
            return
         elseif (actswmach.ne.0.0.and.m0.lt.10.0*origm0) then
c
c
c
            if (s.gt.ssepmin.and.vsep.gt.vsepmin.and.m0.ne.origm0)
     >         then
c
c              Try iterating at smaller mach number
c
               exitcond = 2

               return
c
            else
c
c              Try iterating at higher mach number
c
               exitcond = 1
               snegerr  = s
c
c            write (6,*)'H:',h,hlim,snegerr
c
               return
c
            endif

c
         else
            write (6,*) 'ERROR: delta h less than minimum :',
     >                  h,hlim
            write (6,*) 'Program NOT stopping : ',ierr
c
            exitcond = 4
            snegerr  = s
            return
c
c            stop
c
         endif

c
c         if (s.gt.lastiters) then
c            exitcond = 4
c            snegerr = s
c            write (6,*) 'ERROR: delta h less than minimum :',
c     >                  h,hlim
c            write (6,*) 'Possibly caused by Negative T: ',ierr
c            return
c         else
c            write (6,*) 'ERROR: delta h less than minimum :',
c     >                  h,hlim
c            write (6,*) 'Program NOT stopping : ',ierr
c
c            exitcond = 4
c            snegerr  = s
c            return
c
c            stop
c
c         endif
c
      endif
c
c------------------------------------------------------------------------
c
c
c     Check to see if it was a successful step - if it was go on
c     the next step using the revised h-step - if it wasn't retry
c     with an appropriately revised h step.
c
      errmax = max (abs(errte),abs(errti))
c
      errmax = errmax / eps
c
c      write (6,*) 'errs:',errmax,errte,errti
c
      if (errmax.gt.1.0d0.and.imflag2.eq.0) then
c
c        Error is greater than the desired amount - retry with
c        smaller h.
c
         hnew = 0.9 * h * (errmax ** (-0.25))
         if (hnew.lt.0.1*h) then
            hnew = 0.1 * h
         endif
         s = slast
c
c         if (ringnum.eq.6)
c     >      write (6,'(a,4g12.5)') 'errs:',errmax,errte,errti,hnew
c
      else
c
c        Error is less - can afford to increase h a bit
c
         stmp = slast + h
         if ((abs(stmp-send).le.hlim))
     >      stmp = send
c
         s = stmp
c
         if (imflag2.eq.0) then
            if (errmax.gt.((5.0/0.96)**(-5.0))) then
               hnew = 0.96 * h * (errmax ** (-0.2))
            else
               hnew = 5.0 * h
            endif
         elseif (imflag2.eq.1) then
            hnew = himag
         endif
c
c         write (6,*) 'h1:',h,hnew,imflag2,s
c         write (6,*) 'h2:',send,sinit,ndiv,dble(ndiv)
c
         if ((s+hnew).ge.send) then
            hnew = send - stmp
         endif
c
         tmpgam = gamma(s)
         n = newn(s,newt1e,newt1i,nimag,flag)
c
c        Deal with NEWN return codes
c
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
c
         elseif (flag.eq.2) then
c
c           Negative N encountered - errcode/exitcond=6
c
            exitcond = 6
            snegerr  = s
            return
c
         endif
c
c
c        write(6,*) 'VSEPA:',stopimag,flag,founds,lastiters,
c     >             vsep,vsepmin
c
         if (stopimag.and.flag.eq.0) then
            if (vsep.le.vsepmin.and.(.not.founds)) then
               vsepmin = vsep
               ssepmin = s
               if (vsepmin.lt.5.0) then
                  vcount = vcount + 1
c
c                  if (ringnum.eq.6) then
c                     write(6,*) 'vcount:',s,vsep,vcount
c                  endif
c
                  if (vcount.gt.1) then
                     founds = .true.
                     lastiters = s
c
c                     if (ringnum.eq.6) then
c                        write(6,*) 'Founds = TRUE: S = ',lastiters
c                     endif
c
                  endif
               endif
c
               write(6,*) 'VSEP:',vsepmin,ssepmin,vcount,lastiters

c
c            elseif (s.gt.(sepfact*ssepmin).and.s.gt.(0.1*ringlen))
c
            elseif (s.gt.(sepfact*ssepmin).and.ssepmin.ne.0.0
     >              .and.s.gt.(0.1*halfringlen))
     >              then
c
c             elseif (vsep.gt.sepfact*vsepmin) then
c
c              Need to go back to last mach number that didn't
c              work and start up from there.
c
c               if (ringnum.eq.6) then
c                 write(6,*) 'lastiter:',lastiter,founds,s,ssepmin,
c                           lastiters,m0
c               endif
c
               if (lastiter) then
                  stopimag = .false.
               elseif (m0.ne.origm0) then
                  exitcond = 2
c
c                 If founds is false - still set lastiters to
c                 ssepmin - just in case the soltion with the
c                 specified mach number resolution will not
c                 result in a valid value of lastiters being found.
c
                  if (.not.founds) then
                     lastiters = ssepmin
                  endif
c
                  return
               endif
c
            elseif (halfringlen.lt.(sepfact*ssepmin)) then
c
c             In this case S can never be greater than
c             sepfact * ssepmin and so some other calculation
c             is needed - in this case  1/4 of the remaining
c             distance S to the end of the ring from the
c             point of minimum separation is used.
c
              if (s.gt.(0.25*(halfringlen-ssepmin)+ssepmin)) then
c
c                Need to go back to last mach number that didn't
c                work and start up from there.
c
c                 write(6,*) 'M iter:',lastiter,s,ssepmin,
c                           lastiters,m0
c
                 if (lastiter) then
                    stopimag = .false.
                 elseif (m0.ne.origm0) then
                    exitcond = 2
                    return
                 endif
              endif
c
            else
               vcount = 0
c
            endif
         endif
c
c     End of errmax if statement
c
      endif
c
c
c------------------------------------------------------------------------
c
c
      if ((s.ne.send).and.(abs(send-s).le.hlim)) then
         s = send
      endif
c
      if ((s.eq.send).or.(s.ne.slast)) then
c
         iter = iter + 1
c
c        Update the integrals at the end of
c        each R-K step.
c
         n = newn (s,newt1e,newt1i,nimag,flag)
c
         if (flag.eq.2) then
c
c           Negative N encountered - errcode/exitcond=6
c
            exitcond = 6
            snegerr  = s
            return
c
         endif
c
         ptmp = peiupdt(s,n,nlast,newt1e,newt1i,t1e,t1i,fval)
         ptmp = pradupdt(s,n,nlast,newt1e,t1e)
         ptmp = phelpiupdt(s,n,nlast,newt1e,t1e)
         ptmp = pcxupdt(s,newt1i,t1i)
         ptmp = pintupdt(s,n,newt1e,newt1i)
c
         if (actswnmom.eq.9.or.actswnmom.eq.10) then
            ptmp = scxupdt(s,n,nlast,newt1i,t1i)
         endif
c
         tmpgam1 = gamma(s)

         vtmp = 0.0
         if (n.ne.0.0) vtmp = tmpgam1 / n

         lastvel = vtmp
c
         if (actswnmom.eq.4) then
            smomtmp = pmomloss(s,1,vtmp,newt1e,newt1i)
         endif

c
c        Calculate an estimate of the velocity gradient
c        VGRAD for possible use in future viscosity options.
c
c        if (imflag2.eq.1) then
c           vgrad = 0.0
c        else
            tmpgam1 = gamma(s-h)
            tmpgam2 = gamma(s)
c
            v1=0.0
            if (nlast.ne.0.0) v1 = tmpgam1/nlast

            v2=0.0
            if (n.ne.0.0) v2 = tmpgam2/n
c
            vgrad = (v2-v1) / h
c
c        endif
c
         nlast  = n
         slast  = s
c
         t1i = newt1i
         t1e = newt1e
c
         timax = max(timax,t1i)
         temax = max(temax,t1e)
c
c        Exit if error correction is ON and the
c        temperature drops too much.
c
         if (actswerror.ne.0.0.and.
     >     ((t1i.lt.dropfrac*timax).or.
     >      (t1e.lt.dropfrac*temax).or.
     >      (t1i.le.tfloor).or.
     >      (t1e.le.tfloor))) then
c
            exitcond = 5
            snegerr = s
            return
c
          endif
c
c         sinit = s
c
         if (s.eq.send) return
      endif
c
      if (hnew.le.0.0) then
         write (6,*) 'Hnew = ',hnew,' last H = ',h
         stop
      endif
c
      h = hnew

c
c     If debugging is on - record the data for the current step
c
      if (debug_sol22_on.and.srtn.eq.slast) then 
         call save_s22_data(h,s,n,t1e,t1i,lastvel,
     >                      gamma(s),srci(s),srcf(s))
      endif
c
      goto 2000
c
c     End of Solvstep
c
      end
c
c
c
      subroutine  rkstep(s,t1e,t1i,newt1e,newt1i,errte,errti,h,m0arg,
     >                   ierr)
      implicit none
c
      real*8 s,t1i,t1e,newt1i,newt1e,errte,errti,h,m0arg
      integer ierr
c
      include 'solparams'
      include 'solcommon'
c
      include 'solrk'
      include 'solswitch'
c
c
c     Local variables
c
      integer i,j,flag
      real*8 ke(6),ki(6),nimag,n
      real*8 t2ep,t2ip
      real*8 si,ti,te,fegrad,figrad,newn,fgrad
      external fegrad,figrad,newn,fgrad
c
c     Initialization
c
      ierr = 0
c
c     Loop through calculating each step for te and ti
c
      do i = 1,6
c
c        set up values at start of RK step
c
         si = s + ai(i) * h
         ti = t1i
         te = t1e
c
c        recalculate Ti,Te for the current step
c
         do j = 1,i-1

            ti = ti + bij(i,j) * ki(j)
            te = te + bij(i,j) * ke(j)

         end do
c
c        Test for valid te and ti values (i.e. > 0 ) - if this is
c        not the case - exit - issue an error message and reduce
c        the step-size by a factor of 10.
c
         if ((te.lt.0.0).or.(ti.lt.0.0)) then
            ierr = 1
c            if (ringnum.eq.6) then
c               write(6,*) 'errnega:',m0,h,te,ti,s,si
c     >                    ,lastiter,lastiters
c            endif
            return
         endif
c
c        Calculate n for this step
c
         n = newn (si,te,ti,nimag,flag)
c
c        Check for imaginary n's
c
         if (flag.eq.1) then
c
            ierr = 2
c
            if (((actswmach.eq.1.0.or.actswmach.eq.2.0)
     >        .and.(.not.lastiter))
     >         .or.(lastiter.and.h.gt.himag)
     >         )  return
c
         elseif (flag.eq.2) then
c
c           Negative N encountered - ierr=3 - errcode/exitcond=6
c
            ierr = 3
            return
c
         endif
c
c        Check for NaNQ values in Te, Ti and N
c
         if ((.not.(te.le.0.0.or.te.gt.0.0)).or.
     >       (.not.(ti.le.0.0.or.ti.gt.0.0)).or.
     >       (.not.(n.le.0.0.or.n.gt.0.0))) then
c
             write (6,'(a,3(1x,g13.5))') 'NaNQ Error:',te,ti,n
             ierr = 4
             return
         endif
c
c        Calculate the k-values
c
         if (forcet.eq.0) then
            ke(i) = h * fegrad (si,te,ti,n)
            ki(i) = h * figrad (si,te,ti,n)
         elseif (forcet.eq.1) then
            ke(i) = h * fgrad (si,te,ti,n)
            ki(i) = ke(i)
         endif
c
c         write (6,*) 'rk:',i,s,te,ti
c         write (6,*) 'rk2:',h,n,ke(i),ki(i)

      end do
c
c     If it has been successful - calculate the errors in the Te,Ti
c     estimates
c
      newt1e = t1e
      newt1i = t1i
      t2ep = t1e
      t2ip = t1i
c
      do i = 1,6
         newt1e = newt1e + ci(i) * ke(i)
         newt1i = newt1i + ci(i) * ki(i)
         t2ep = t2ep + cip(i) * ke(i)
         t2ip = t2ip + cip(i) * ki(i)
      end do
c
c
c        Test for valid te and ti EXIT values (i.e. > 0 ) - if this
c        is not the case - return - issue an error message and reduce
c        the step-size by a factor of 10 in the calling routine.
c
      if ((newt1e.lt.0.0).or.(newt1i.lt.0.0).or.
     >    (t2ep.lt.0.0).or.(t2ip.lt.0.0)) then
         ierr =1
c         write(6,*) 'errneg2:',m0,h,s,newt1e,newt1i,
c     >                            t2ep,t2ip
         return
      endif
c
c
      errte = newt1e - t2ep
      errti = newt1i - t2ip
c
c     NOTE: Due to constraints in supplementary databases and
c           elsewhere - we do not want to return a value of Te
c           or Ti less than an imposed minimum. If the values
c           calculated for a step are less than that value then
c           they are set equal to that value.
c
      if (newt1e.lt.tfloor.and.newt1e.gt.0.0) newt1e = tfloor
      if (newt1i.lt.tfloor.and.newt1i.gt.0.0) newt1i = tfloor
c
c     exit
c
c 1000 continue
c
      return
      end
c
c
c
      real*8 function vgradval(s,ind)
      implicit none
      real*8 s
      integer ind
      include 'solparams'
      include 'solcommon'
      include 'solswitch'
c
c     This function returns the current best estimate
c     of the velocity gradient. This is used in conjunction
c     with the vicosity solution options.
c
      vgradval = vgrad
c
      return
      end
c
c
c
      real*8 function fegrad(s,te,ti,n)
      implicit none
c
c     This function returns the approximate derivate dte/ds and is
c     used in the Runge-Kutta approximation method in attempting to
c     estimate the functional value at the next step.
c
      include 'solparams'
      include 'solcommon'
c
      real*8 s,te,ti
c
      real*8 n,cond,estprad,estphelpi,estpei,estppelec,paes
      external cond,estprad,estphelpi,estpei,estppelec,paes
c
      fegrad = (1.0/(k0e*te**2.5)) * (cond(s,te)+estprad(s,n,te)
     >         + estphelpi(s,n,te) + estpei(s,n,te,ti)
     >         + estppelec(s) + paes(s))
      return
      end
c
c
c
      real*8 function figrad(s,te,ti,n)
      implicit none
c
c     This function returns the approximate derivate dte/ds and is
c     used in the Runge-Kutta approximation method in attempting to
c     estimate the functional value at the next step.
c
      include 'solparams'
      include 'solcommon'
c
      real*8 s,te,ti
c
      real*8 n,cond,conv,estpcx,estpei,estppion,pais
      external cond,conv,estpcx,estpei,estppion,pais
c
      figrad = (1.0/(k0i*ti**2.5)) * (cond(s,ti)+ conv(s,n,ti) +
     >         estpcx(s,ti) - estpei(s,n,te,ti)
     >         + estppion(s) + pais(s))
c
      return
      end
c
c
c
      real*8 function fgrad(s,te,ti,n)
      implicit none
c
c     This function returns the approximate derivate dt/ds for both species
c     combined and is used in the Runge-Kutta approximation
c     method in attempting to
c     estimate the functional value at the next step.
c
      include 'solparams'
      include 'solcommon'
c
      real*8 s,te,ti
c
      real*8 n,cond,estprad,estphelpi,estpcx,conv,paes,pais
      real*8 estppelec,estppion
      external cond,estprad,estphelpi,estpcx,conv,paes,pais
      external estppelec,estppion
c
      fgrad = (1.0/(k0e*te**2.5+k0i*ti**2.5)) * (cond(s,te)+ cond(s,ti)
     >         + conv(s,n,ti) + estprad(s,n,te) + estpcx(s,ti)
     >         + estphelpi(s,n,te) + estppelec(s) + estppion(s)
     >         + paes(s) + pais(s))
c
      return
      end
c
c
c     WF'95: DISTRIBUTE POWER ALONG S
c
      real*8 function paes(s)
      implicit none
      real*8 s,majrpos,areaint
      external majrpos,areaint
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 rfact,rpos
c
      paes = pae
c
      if (actswmajr.eq.4.0) then
         rpos = majrpos(s)
         rfact = r0init/rpos
      else
         rfact = 1.0
      endif
c
      if (actswpow.eq.0.0) then
         paes =  rfact * pae
      elseif (actswpow.eq.1.0) then
         paes = pae * (1.0 - s/(halfringlen))
      elseif (actswpow.eq.2.0) then
         if (s.le.sxp) then
            paes = pae
         else
            paes = pae * (1.0 - (s-sxp)/(halfringlen-sxp))
         endif
      elseif (actswpow.eq.3.0) then
         paes =  rfact * pae - pinpute * areaint(s) / rpos
      elseif (actswpow.eq.4.0) then
         paes = rfact* (pae
     >          - (pae + pae_end)*s/ringlen)
      elseif (actswpow.eq.5.0) then
         paes = rfact*(pae * (1.0 - s/(halfringlen))
     >          - qesum * s / (halfringlen))
      elseif (actswpow.eq.6.0) then
         paes = rfact* ((pae
     >          - (pae + pae_end)*s/ringlen)
     >          - qesum * s / (halfringlen))
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
            paes = rfact*(pae * (1.0 - (s-spow)/(halfringlen-spow))
     >          - qesum * (s-spow) / (halfringlen-spow))
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
            paes = rfact*(pae * (1.0 - (s-spow)/(spow2-spow))
     >          - qesum * (s-spow) / (spow2-spow))
         else
            paes = 0.0
         endif
c
      elseif (actswpow.eq.11.0) then
c
         if (s.le.spow) then
            paes = pae
         else
            paes = rfact* (pae
     >         - (pae + pae_end)*(s-spow)/(ringlen-2.0*spow))
         endif
c
      endif
c
      return
      end
c
c
c
      real*8 function pais(s)
      implicit none
      real*8 s,majrpos,areaint
      external majrpos,areaint
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 rfact,rpos
c
c
      pais = pai
c
      if (actswmajr.eq.4.0) then
         rfact = r0init/majrpos(s)
      else
         rfact = 1.0
      endif
c
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
         pais = rfact* (pai
     >          - (pai + pai_end)*s/ringlen)
      elseif (actswpow.eq.5.0) then
         pais = rfact*(pai * (1.0 - s/(halfringlen))
     >          - qisum * s / (halfringlen))
      elseif (actswpow.eq.6.0) then
         pais = rfact* ((pai
     >          - (pai + pai_end)*s/ringlen)
     >          - qisum * s / (halfringlen))
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
            pais = rfact*(pai * (1.0 - (s-spow)/(halfringlen-spow))
     >          - qisum * (s-spow) / (halfringlen-spow))
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
            pais = rfact*(pai * (1.0 - (s-spow)/(spow2-spow))
     >          - qisum * (s-spow) / (spow2-spow))
         else
            pais = 0.0
         endif
      elseif (actswpow.eq.11.0) then
c
         if (s.le.spow) then
            pais = pai
         else
            pais = rfact* (pai
     >         - (pai + pai_end)*(s-spow)/(ringlen-2.0*spow))
         endif
c
      endif
c
      return
      end
c
c     Distribute the PP electron target power loss over the ring
c
      real*8 function estppelec(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     ESTPPELEC: This routine calculates the amount of power transferred
c                to the private plasma as a function of S from the
c                current main SOL ring.
c
c                Since everything is analytic at this point and does not
c                depend on the changing plasma conditions for any of
c                it's options - it does not have an update method as the
c                other options may require.
c
c
      estppelec = 0.0
c
c     Option is OFF
c
      if (actswppelec.eq.0.0) then
c
         estppelec = 0.0
c
c     Power added below X-point
c
      elseif (actswppelec.eq.1.0) then
c
         if (s.lt.sxp) then
            estppelec  = s/sxp * ppelecpow
         else
            estppelec  = ppelecpow
         endif
c
c     Power added to specified distance along field line
c
      elseif (actswppelec.eq.2.0) then
c
         if (s.lt.pp_pow_dist*ringlen) then
            estppelec  = s/(pp_pow_dist*ringlen) * ppelecpow
         else
            estppelec  = ppelecpow
         endif
c
      endif
c
c     End of routine
c
      return
      end

c
c     Distribute the PP ion target power loss over the ring
c
      real*8 function estppion(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     ESTPPION:  This routine calculates the amount of power transferred
c                to the private plasma as a function of S from the
c                current main SOL ring.
c
c                Since everything is analytic at this point and does not
c                depend on the changing plasma conditions for any of
c                it's options - it does not have an update method as the
c                other options may require.
c
c
      estppion = 0.0
c
c     Option is OFF
c
      if (actswppion.eq.0.0) then
c
         estppion = 0.0
c
c     Power added below X-point
c
      elseif (actswppion.eq.1.0) then
c
         if (s.lt.sxp) then
            estppion  = s/sxp * ppionpow
         else
            estppion  = ppionpow
         endif
c
c     Power added to specified distance along field line
c
      elseif (actswppion.eq.2.0) then
c
         if (s.lt.pp_pow_dist*ringlen) then
            estppion  = s/(pp_pow_dist*ringlen) * ppionpow
         else
            estppion  = ppionpow
         endif
c
      endif
c
c     End of routine
c
      return
      end
c
c     Distribute the PP target pressure loss over the ring
c
      real*8 function estppress(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     ESTPPRESS: This routine calculates the amount of pressure transferred
c                to the private plasma and then distributes this 
c                as a function of S from the
c                current main SOL ring.
c
c                Since everything is analytic at this point and does not
c                depend on the changing plasma conditions for any of
c                it's options - it does not have an update method as the
c                other options may require.
c
c
      estppress = 0.0
c
c     Option is OFF
c
      if (actswppress.eq.0.0) then
c
         estppress = 0.0
c
c     Power added below X-point
c
      elseif (actswppress.eq.1.0) then
c
         if (s.lt.sxp) then
            estppress  = s/sxp * pp_press
         else
            estppress  = pp_press
         endif
c
c     Power added to specified distance along field line
c
      elseif (actswppress.eq.2.0) then
c
         if (s.lt.pp_pow_dist*ringlen) then
            estppress  = s/(pp_pow_dist*ringlen) * pp_press
         else
            estppress  = pp_press
         endif
c
      endif
c
c     End of routine
c
      return
      end
c
c
c
      real*8 function areaint(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     Returns the estimates R(s)*ds  integral at position S.
c
      integer top, bot,mid, in
c
c
c      if (s.gt.sptscopy(nptscopy)) then
c
c         in = nptscopy +1
c
c         areaint = ( intarea(in-1) +
c     >             ( (intarea(in)-intarea(in-1)) *
c     >             (s-sptscopy(in-1))/((ringlen/2.0)-sptscopy(in-1))))
c
c      else
c
c      BINSEARCH returns in+1 for S > spts(npts) 
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            areaint = intarea(in) * s / sptscopy(1)
         else
            areaint = ( intarea(in-1) +
     >             ( (intarea(in)-intarea(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
         endif
c
c      endif
c



      return
      end

c
c
c
      real*8 function newn(s,te,ti,nimag,flag)
      use sol22_debug
      implicit none
c
c     Calculates the density value from the given parameters by solving
c     the quadratic equation for N - issues an error when part is
c     imaginary and depending on options - may take corrective action.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 s,te,ti,nimag
      integer flag
c
      real*8 imag, rest,gamma,tmpgam,tmpv,press,tmppress
      real*8 posfact,majrpos,ptmp
      external gamma,press,majrpos
c
      flag = 0
c
      if (s.eq.0) then
         newn = n0
         nimag = m0**2 * n0
         return
      endif
c
      tmpgam = gamma(s)
      tmppress = press(s,te,ti)
c
      imag = (tmppress/ (econv*(te+ti)) )**2
     >       - 4.0* ((tmpgam * mconv) / econv) * (mb * tmpgam) / (te+ti)
      rest = tmppress /(econv*(te+ti))
c
      if (imag.lt.0.0) then
c
c     set an error recovery condition ... if n becomes imaginary
c     then set v = cs or v = lastvel and calculate n accordingly.
c     OR adjust the pressure by adding enough so that the solution
c     is not imaginary.
c
c
         if (velsw.eq.0) then
            tmpv =   sqrt( (te+ti)/mb * econv/mconv)
            newn =  abs( tmpgam / tmpv)
            nimag = newn
         elseif (velsw.eq.1) then
            tmpv =  lastvel
            newn =  abs( tmpgam / tmpv)
            nimag = newn
         elseif (velsw.eq.2) then
c
            ptmp = 4.0 * (econv*(te+ti)) * (mb * tmpgam)
     >                             *  (tmpgam *mconv)
            padd = ptmp - tmppress + padd
c
            newn = ptmp / (econv*(te+ti))
            nimag= newn
            tmpv = tmpgam/newn
c
         elseif (velsw.eq.3) then
c
            ptmp = 4.0 * (econv*(te+ti)) * (mb * tmpgam)
     >                             *  (tmpgam *mconv)
c
c            padd = ptmp - tmppress
c
            newn = ptmp / (econv*(te+ti))
            nimag= newn
            tmpv = tmpgam/newn
c
         endif
c
         if (newn.lt.0.0) then
            write (6,*) 'Flow reversal?:',newn,tmpgam,tmpv,te,ti
            newn = -newn
            tmpv = -tmpv
         endif
c
         flag = 1
c         write (6,1000) 'newn0:',te,ti,
c     >                imag,tmpv,newn,s
c
c
c        Debug
c
c      if (ringnum.eq.8.and.pinavail) then
c       write (6,'(a,8(1x,g13.5))') 'Newn:I',s,imag,rest,tmpgam,te,ti,
c     >       rest**2,
c     >     - 4.0*((tmpgam * mconv) / econv)*(mb * tmpgam)/(te+ti)
c      endif
c
      else
c
c        Modify so that newn < 0 returns a flag / exit condition
c
         if (founds.or.(lastiter.and.s.gt.lastiters)) then
c
            newn = (rest + sqrt(imag))/2.0
            nimag = (rest - sqrt(imag))/2.0
c
            if (newn.lt.0.0) then
c
               write (6,'(a,10(1x,g11.5))') 'Newn<0:',s,imag,rest,
     >              tmppress,tmpgam,te,ti,
     >               rest**2,-4.0*((tmpgam * mconv) / econv)
     >               *(mb * tmpgam)/(te+ti),newn
               write(6,*) 'founds:',founds,lastiter,lastiters
c
               newn = abs(newn)
               flag = 2
               return
c
            endif
c
c        Debug
c
           if (debug_sol22_on.and.pinavail) then
             write (6,'(a,10(1x,g14.6))') 'Newn:I',s,imag,rest,
     >            tmpgam,te,ti, rest**2,
     >            - 4.0*((tmpgam * mconv) / econv)*(mb * tmpgam)/(te+ti)
           endif
         else
c
            newn = (rest - sqrt(imag))/2.0
            nimag = (rest + sqrt(imag))/2.0
c
            if (newn.lt.0.0) then
c
               write (6,*) 'Error: Newn < 0 :'
               write (6,'(a,9(1x,g13.5))') 'Newn:I',s,imag,rest,
     >              tmppress,tmpgam,te,ti,
     >               rest**2,-4.0*((tmpgam * mconv) / econv)
     >               *(mb * tmpgam)/(te+ti)
               write(6,*) 'founds:',founds,lastiter,lastiters
c
               newn = abs(newn)
c
               flag = 2
               return
c
            endif
         endif
      endif
c
      return
      end
c
c
c
      real*8 function press(s,tecur,ticur)
      implicit none
      real*8 s,tecur,ticur
      include 'solparams'
      include 'solcommon'
      include 'solswitch'
c
c     This function returns the value of the pressure at a position
c     s along the field line, at the moment the only contribution
c     other than the Pinf is the Momentum loss to neutrals term.
c
      real*8 pmomloss,estpint,rpos,majrpos,pmomloss_tmp,estppress
      external pmomloss,estpint,majrpos,estppress
c
      pmomloss_tmp = pmomloss(s,0,1.0d0,tecur,ticur)
c
c
c     jdemod - pmomloss returns 0.0 if the option is off so this 
c              extra check code is not required. 
c
c      if (actswnmom.eq.0.0) then
c         if (actswmajr.eq.4.0) then
c            rpos = majrpos(s)
c            press = (pinf + estpint(s)+ padd)/rpos
c         else
c            press = pinf + padd
c         endif
c      else
c
         if (actswmajr.eq.4.0) then
            rpos = majrpos(s)
            press = (pinf + estpint(s) + pmomloss_tmp + padd 
     >               + estppress(s))
     >              /rpos
         else
            press = pinf + pmomloss_tmp + padd +estppress(s)
         endif
c
c      endif
c
      return
      end
c
c
c
      real*8 function pintupdt(s,n,te,ti)
      implicit none
      real*8 s,n,te,ti
c
c     PINTUPDT: This routine updates the accumulated integral
c               of the major radius pressure correction term.
c
c     Note: Initialization done by call from INITVAL
c
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      common /pint/ lastpint,lasts,rlasts,plasts,nlasts
      real*8 lastpint,lasts,rlasts,plasts,nlasts
c
      real*8 rpos,majrpos,pcur
      external majrpos
c
c     Only if major radius correction option is ON
c
      if (actswmajr.ne.4.0) then
         pintupdt = 0.0
         return
      endif
c
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
c
      rpos = majrpos(s)
      pcur = n * (te+ti) * econv
c
      pintupdt = lastpint + 0.5 * (pcur+plasts) * (rpos-rlasts)
c
      rlasts = rpos
      lasts = s
      lastpint = pintupdt
      plasts = pcur
c
      return
      end
c
c
c
      real*8 function estpint(s)
      implicit none
      real*8 s
c
c     ESTPINT: This routine returns an estimate of the major radius
c              integrated pressure correction term at a value
c              of S - since the last update.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      common /pint/ lastpint,lasts,rlasts,plasts,nlasts
      real*8 lastpint,lasts,rlasts,plasts,nlasts
c
      real*8 rpos,majrpos
      external majrpos
c
c     Only if major radius correction option is ON
c
      if (actswmajr.ne.4.0) then
         estpint = 0.0
         return
      endif
c
      if (s.eq.0.0.or.s.lt.lasts) then
         estpint = 0.0
         return
      elseif (s.eq.lasts) then
         estpint = lastpint
         return
      endif
c
      rpos = majrpos(s)
c
      estpint = lastpint + plasts * (rpos-rlasts)
c
      return
      end
c
c
c
      real*8 function pmomloss(s,opt,vcur,tecur,ticur)
      implicit none
      real*8 s,vcur,tecur,ticur
      integer opt
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      include 'sol22pmom'
c
c     This function will provide the momentum loss
c     integrated to a point s. The options that depend
c     on data transferred from Nimbus have been partially
c     implemented.
c
c     The OPT argument is used when one requires the
c     UPDATE/ESTIMATOR method of calculating an integral
c     of the momentum source term. At present it is only
c     of use for momentum option 4. However the mechanism
c     may be employed for other momentum options as required.
c     It is implemented in the same fashion as the ancillary
c     power integration terms.
c
c     A value of 0 is a request for an estimate and a value of
c     1 is an indicator to update the integral.
c
c      real*8 :: lasts,lastsmom,lastsrc,lastv,lastte
      real*8 src,srci,teav,rcxmult,estscx,vav
      external srci,rcxmult,estscx
      real*8 :: vtmp,te_base
c
      logical temp_opt

c
      integer top,bot,mid,in
c
c     Initialize and set to zero for S = 0
c
      pmomloss = 0.0
c
c     Exit for S = 0
c
      if (s.eq.0) then 
         lastv = vcur
         lasts = s
         lastte = tecur
         return
      endif
c
      if (actswnmom.eq.0.0.or.
     >   (actswnmom.eq.6.0.and.(.not.pinavail))) then
         pmomloss = 0.0
      elseif (actswnmom.eq.1.0.or.
     >   (actswnmom.eq.7.0.and.(.not.pinavail))) then
         if (s.le.actlenmom*ringlen) then
            pmomloss = s * smom0
         else
            pmomloss = actlenmom * ringlen * smom0
         endif
      elseif (actswnmom.eq.2.0.or.
     >   (actswnmom.eq.8.0.and.(.not.pinavail))) then
         if (s.le.actlenmom*ringlen) then
            pmomloss = smom0 * lammom * ringlen *
     >              (1.0d0 - exp (-s/(lammom*ringlen)))
         else
            pmomloss = smom0 * lammom * ringlen *
     >       (1.0d0 - exp (-actlenmom/lammom))
c    >       (1.0d0 - exp (-(actlenmom*ringlen)/(lammom*ringlen)))
         endif
      elseif (actswnmom.eq.3.0) then
         pmomloss = smom0 * srci(s)
      elseif (actswnmom.eq.4.0) then
c
c        The momentum loss is calculated using the term
c        Smom = -m * vb * rcxmom * Siz(s)
c
c
c        Return an estimate of Smom in the next interval
c
         temp_opt=.true.
         te_base = 10.0
         
         if (temp_opt) then 
            vtmp = -sqrt(2.0*te_base*econv/(mb*mconv))
         else
            ! set to absolute so that flow reversal doesn't cause a pressure drop
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
c
            src = srci(s)
c
            pmomloss = lastsmom
     >               - mb * mconv * (vtmp+lastv)/2.0
     >               * rcxmom * rcxmult(lastte)
     >               * (src-lastsrc) * smom_mult

c
c           This now nees its own return statement to avoid
c           double multiplication by smom_mult
c
c            if (ringnum.eq.17)
c     >         write(6,'(a,16(1x,g12.5))') ' Mom0:',
c     >         s,lasts,lastsmom,lastv,vtmp,(vtmp+lastv)/2.0,src,lastsrc,
c     >          (src-lastsrc),
c     >          lastte,rcxmom,rcxmult(lastte),
c     >          pmomloss,
c     >            - mb * mconv * lastv * rcxmom * (src-lastsrc) 
c     >             *smom_mult *rcxmult(lastte)
c

            return
c
         elseif (opt.eq.1) then
c
            vav = (vtmp + lastv)/2.0
            teav = (tecur+lastte)/2.0
c
            if (s.eq.soffset.or.s.lt.lasts) then
c
               lasts = soffset
               lastte = tecur
               lastv = vtmp
               lastsmom = 0.0
               lastsrc = 0.0
               pmomloss = 0.0

c               write (6,'(a,10g12.5)') 'Mom lasts:',s,lasts,
c     >                         lastsrc,lastsmom,lastv,vtmp
c
               return
            elseif (s.eq.lasts) then
               pmomloss = lastsmom
               return
            endif
c
            src = srci(s)
c
            pmomloss = lastsmom
     >               - mb * mconv * vav
     >               * rcxmom * rcxmult(teav)
     >               * (src-lastsrc) * smom_mult
c
c            if (ringnum.eq.17)
c     >         write(6,'(a,16(1x,g12.5))') ' Mom1:',
c     >          s,lasts,lastsmom,lastv,vtmp,vav,src,lastsrc,
c     >          (src-lastsrc),
c     >          teav,rcxmom,rcxmult(teav),
c     >          pmomloss,
c     >            - mb * mconv * vav * rcxmom * rcxmult(teav) 
c     >            * (src-lastsrc)* smom_mult
c
            lastv = vtmp
            lastte = tecur
            lastsmom = pmomloss
            lasts = s
            lastsrc = src
c
c           This now nees its own return statement to avoid
c           double multiplication by smom_mult
c
            return
c
         endif
c
      elseif (actswnmom.eq.5.0) then
         pmomloss = 0.0
      elseif (pinavail.and.(actswnmom.eq.6.0.or.
     >        actswnmom.eq.7.or.
     >        actswnmom.eq.8.0)) then
c
c        Options based on PIN when pin data is
c        available.
c
c        Search for right cell
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            pmomloss = intmomsrc(in) * s / sptscopy(1)
         else
            pmomloss = ( intmomsrc(in-1) +
     >             ( (intmomsrc(in)-intmomsrc(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
         endif
c
      elseif (pinavail.and.(actswnmom.eq.9.0
     >        .or.actswnmom.eq.10.0)) then
c
c     WF'96: estimates pmomloss from NIMBUS nH, EH, and CXSIG()
c
         pmomloss = estscx(s,nlast,ticur)
c
      endif
c
c     Apply overall multiplier to MOST options - some options exit before
c     reaching this point.
c
      pmomloss = pmomloss * smom_mult
c
c      if (ringnum.eq.17) 
c     >     write(6,'(a,16(1x,g12.5))') ' Mom2:',
c     >        s,lasts,lastsmom,smom0,srci(s),smom_mult,pmomloss
c
c
      return
      end
c
c
c
      real*8 function rcxmult(t)
      implicit none
      real*8 t
      include 'solparams'
      include 'solcommon'
c
c     This function calculates a Cx/IZ multiplier for use
c     in the neutral momentum loss term that contributes
c     to the total pressure.
c
      integer firstx
      data firstx /0/
c
      real*8 coeffa,coeffb
c
c     Initialize coefficients on first iteration
c
      if (firstx.eq.0) then
         coeffb = 6.907/(tcxmom-1.0)
         coeffa = 1000.0 * exp(coeffb)
         firstx = 1
      endif
c
      if (t.le.tcxcut) then
         rcxmult = 1.0d0
      else
c
         rcxmult = coeffa * exp(-coeffb*t)
c
         if (rcxmult.lt.1.0) rcxmult = 1.0d0
         if (rcxmult.gt.1500.0) rcxmult = 1500.0d0
c
      endif
c
      return
      end
c
c
c
      real*8 function cond(s,t)
      implicit none
c
c     Returns the first convective energy component
c
      real*8 s,t
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 gamma
      external gamma
c
      if (actswcond.eq.0.0) then
         cond = 0.0d0
         return
      endif
c
      cond = 5.0d0/2.0d0 * gamma(s) * econv * t
c
      return
      end
c
c
c
      real*8 function conv(s,n,t)
      implicit none
c
c     Calculates kinetic convective term
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 s,n,t
c
      real*8 tmpgam,gamma
      external gamma
c
      if (actswconv.eq.0.0) then
         conv = 0.0d0
         return
      endif
c
      tmpgam = gamma(s)
c
      conv = 0.5d0 * mb * mconv * (tmpgam/n)**2 * tmpgam
c
c      write(6,'(a6,6g14.6)') 'conv:',mb,mconv,
c     >                       tmpgam,tmpgam/n,n,conv
c
      return
      end
c
c
c
c      real*8 function prad(s)
c      implicit none
c
c     Returns the integrated radiative power loss using only a simple
c     exponential
c     for now.
c
c      real*8 s
c
c      include 'solparams'
c      include 'solswitch'
c      include 'solcommon'
c
c      if (actswprad.eq.0.0) then
c         prad = 0.0
c         return
c      elseif (actswprad.eq.1.0) then
c
c         if (s.lt.lenr) then
c            prad = lamr*prad0*(1-exp(-s/lamr))
c         else
c            prad = lamr*prad0*(1-exp(-lenr/lamr))
c         endif
c
c      endif
c
c      return
c      end
c
      real*8 function pradupdt(s,n,nold,te,teold)
      implicit none
c
c     PRADUPDT: This returns the integrated value of the
c               radiation losses up to the point s.
c
c     This function updates the lasts and lastprad values after
c     calculating a Prad contribution based on the average temperatures
c     over the interval.
c
      real*8 s,n,te,teold,nold
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      common /praddata/ lastprad,lasts
      real*8 lastprad,lasts
c
      real*8 teav,nav
c
      integer in
c
      if (actswprad.eq.0.0) then
         pradupdt = 0.0
c
c         return
c
      elseif (actswprad.eq.1.0) then
c
         if (s.lt.lenr) then
            pradupdt = lamr*prad0*(1-exp(-s/lamr))
         else
            pradupdt = lamr*prad0*(1-exp(-lenr/lamr))
         endif
c
c         return
c
      elseif (actswprad.eq.2.0) then
c
         if (s.eq.soffset.or.s.lt.lasts) then
            lasts = soffset
            lastprad = 0.0
            pradupdt = 0.0
            return
         elseif (s.eq.lasts) then
            pradupdt = lastprad
            return
         endif
c
         teav = (te+teold)/2.0 / talimp
         teav = max(1.0d-6,teav)
c
         nav  = (n+nold)/2.0
c
         pradupdt = lastprad +  alfimp *  nav**2
     >            * 2.0d-31 / (teav**ex1imp+teav**ex2imp)
     >            * (s-lasts)
c
c        write(6,*) 'pradupdt1:',pradupdt,lastprad,s,lasts,nav,teav
c
         lastprad = pradupdt
         lasts = s
c
c     PRAD = RADSRC_MULT * PINQE
c
      elseif (actswprad.eq.3.0) then

         call binsearch(s,in)
c
         if (in.eq.1) then
            pradupdt = (intqe(in) * s / sptscopy(1))
         else
            pradupdt = (
     >             ( intqe(in-1) +
     >             ( (intqe(in)-intqe(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
c         pradupdt = radsrc_mult * pradupdt
c
c
c     PRAD = RADSRC_MULT * (EXTERNAL RADIATION SOURCE)
c
      elseif (actswprad.eq.4.0) then

         call binsearch(s,in)
c
         if (in.eq.1) then
            pradupdt = (intrad(in) * s / sptscopy(1))
         else
            pradupdt = (
     >             ( intrad(in-1) +
     >             ( (intrad(in)-intrad(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
c         pradupdt = radsrc_mult * pradupdt
c
      elseif (actswprad.eq.6.0) then
c
c        re-use lamr and lenr for a rectangular radiation source
c
         if (s.lt.lamr) then
            pradupdt = 0.0
         elseif (s.lt.lenr) then 
            pradupdt = (s-lamr)/(lenr-lamr)*prad0
         else
            pradupdt = prad0
         endif
      endif

      pradupdt = radsrc_mult * pradupdt

c
      return
      end
c
c
c
      real*8 function estprad(s,n,te)
      implicit none
c
c     ESTPRAD: This returns the radiative energy loss -
c     approximately integrated to the current point s.
c
c     This function returns an ESTIMATE of the Prad contribution at
c     the given S-position ... the lastprad and lasts values are
c     updated after every R-K iteration.
c
      real*8 s,n,te
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      common /praddata/ lastprad,lasts
      real*8 lastprad,lasts
c
      real*8 teav
c
      integer in
c
      if (actswprad.eq.0.0) then
         estprad = 0.0
c
c         return
c
      elseif (actswprad.eq.1.0) then
c
         if (s.lt.lenr) then
            estprad = lamr*prad0*(1-exp(-s/lamr))
         else
            estprad = lamr*prad0*(1-exp(-lenr/lamr))
         endif
c
c         return
c
      elseif (actswprad.eq.2.0) then
c
         if (s.eq.soffset.or.s.lt.lasts) then
            estprad = 0.0
            return
         elseif (s.eq.lasts) then
            estprad = lastprad
            return
         endif
c
         teav = te / talimp
         teav = max(1.0d-6,teav)
c
         estprad = lastprad +  alfimp *  n**2
     >            * 2.0d-31 / (teav**ex1imp+teav**ex2imp)
     >            * (s-lasts)
c
c     PRAD = PINQE * MULT
c
      elseif (actswprad.eq.3.0) then
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            estprad = (intqe(in) * s / sptscopy(1))
         else
            estprad = (
     >             ( intqe(in-1) +
     >             ( (intqe(in)-intqe(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
      elseif (actswprad.eq.4.0) then

         call binsearch(s,in)
c
         if (in.eq.1) then
            estprad = (intrad(in) * s / sptscopy(1))
         else
            estprad = (
     >             ( intrad(in-1) +
     >             ( (intrad(in)-intrad(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
      elseif (actswprad.eq.6.0) then
c
c        re-use lamr and lenr for a rectangular radiation source
c
         if (s.lt.lamr) then
            estprad = 0.0
         elseif (s.lt.lenr) then 
            estprad = (s-lamr)/(lenr-lamr)*prad0
         else
            estprad = prad0
         endif
      endif
c
c     apply rad src multiplier
c
      estprad = radsrc_mult * estprad
c
c      write(6,*) 'estprad:',estprad,lastprad,s,lasts,n,te
c
      return
      end
c
c
c
      real*8 function peiupdt(s,n,nold,te,ti,
     >                      teold,tiold,fval)
      implicit none
c
c     This returns the electron-ion energy exchange component -
c     approximately integrated to the current point s.
c
c     This function updates the lasts and lastpei values after
c     calculating a Pei contribution based on the average temperatures
c     over the interval.
c
      real*8 s,n,te,ti,teold,tiold,nold,lnlam,fval
      external lnlam
c
      include 'solparams'
      include 'solswitch'
      include 'sol22pei'
      include 'solcommon'
c
c
c      common /pei/ lastpei,lasts
c      real*8 lasts,lastpei
c
      real*8 tiav,teav,nav
c
      if (actswpei.eq.0.0) then
         peiupdt = 0.0
         return
      endif
c
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
c
      teav = (te+teold)/2.0
      tiav = (ti+tiold)/2.0
      nav  = (n+nold)/2.0
c
      peiupdt = lastpei +  peicf * ((1.14e-32 * nav**2 *
     >               ( teav - tiav ))
     >              / (mb* teav**(1.5))) *(s-lasts)
     >                 * lnlam(nav,teav) / 15.0
c
c      write(6,*) 'peiupdt1:',peiupdt,lastpei,s,lasts,nav,teav,tiav
c
      lastpei = peiupdt
      lasts = s
c
      if (actswpei.eq.3.0) then
         fval = peiupdt
         peiupdt = 0.0
      endif
c
      return
      end
c
c
c
      real*8 function estpei(s,n,te,ti)
      implicit none
c
c     This returns the electron-ion energy exchange component -
c     approximately integrated to the current point s.
c
c     This function returns an ESTIMATE of the Pei contribution at
c     the given S-position ... the lastpei and lasts values are
c     updated after every R-K iteration.
c
      real*8 s,n,te,ti
c
      include 'solparams'
      include 'solswitch'
      include 'sol22pei'
      include 'solcommon'
c
c
c      common /pei/ lastpei,lasts
c      real*8  lasts,lastpei
c
      real*8 lnlam
      external lnlam
c
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
c
      estpei = lastpei + peicf * ((1.14e-32 * n**2 * (te-ti))
     >                 / (mb * te**(1.5))) * (s-lasts)
     >                 * lnlam(n,te) / 15.0
c
c      write(6,*) 'estpei:',estpei,lastpei,s,lasts,n,te,ti
c
      return
      end
c
c
c
      real*8 function lnlam(n,te)
      implicit none
      real*8 n,te
c
c     LNLAM: This function returns the value of Lambda used in the
c            Pei energy transfer formulae.
c
      if (n.le.0.0.or.te.le.0.0.or.
     >   ((1.5e13 * te**(1.5) / sqrt(n)).le.1.0)) then
         lnlam = 15.0
      else
         lnlam = log(1.5e13 * te**(1.5) / sqrt(n))
      endif
c
      return
      end
c
c
c
      real*8 function phelpiupdt(s,n,nold,t,told)
      implicit none
c
c     Calculates losses to electrons due to hydrogenic cooling
c
      real*8 s,n,t,told,nold
c
      include 'solparams'
      include 'solswitch'
      include 'sol22phelpi'
      include 'solcommon'
c
c
c      common /phelpi/ lasts,lastphelp,lastsrc
c      real*8 lasts,lastphelp,lastsrc
c
      real*8 tav,nav
c
      real*8  helpi,src,srci
      external srci
      integer in
c
      if (actswphelp.eq.0.0
     >    .or.(actswphelp.eq.3.0.and.(.not.pinavail))) then
         phelpiupdt = 0.0
         return
      elseif(actswphelp.eq.1.0.or.
     >      (actswphelp.eq.2.0.and.(.not.pinavail))) then

c
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
c
         tav = (t+told)/2.0
         nav = (n+nold)/2.0
c
         helpi = 17.5 + (5.0+37.5/tav)*(1.0+ 0.25/tav)*log10(1e21/nav)
         src = srci(s)
c
         phelpiupdt = lastphelp + helpi * econv * (src-lastsrc)
c
         lastphelp = phelpiupdt
         lasts = s
         lastsrc = src
c
      elseif (actswphelp.eq.2.0.or.actswphelp.eq.3.0) then
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            phelpiupdt = (intqe(in) * s / sptscopy(1))
         else
            phelpiupdt = (
     >             ( intqe(in-1) +
     >             ( (intqe(in)-intqe(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
c         if (ringnum.eq.7) then
c            write (6,'(a,i4,6(1x,g13.6))')
c     >       'E:',in,s,sptscopy(in),sptscopy(in-1),intqe(in),
c     >        intqe(in-1),phelpiupdt
c         endif
c

c
      endif

      return
      end
c
c
c
      real*8 function estphelpi(s,n,t)
      implicit none
c
c     Calculates losses to electrons due to hydrogenic cooling
c     Returns an ESTIMATE at the given s,n,t but does not
c     update the values for the interval.
c
      real*8 s,n,t
c
      include 'solparams'
      include 'solswitch'
      include 'sol22phelpi'
      include 'solcommon'
c
c
c      common /phelpi/ lasts,lastphelp,lastsrc
c      real*8 lasts,lastphelp,lastsrc
c
c
      real*8  helpi,src,srci
      external srci
c
      integer in
c
      if (actswphelp.eq.0.0.or.
     >    (actswphelp.eq.3.0.and.(.not.pinavail))) then
         estphelpi = 0.0
         return
      elseif(actswphelp.eq.1.0.or.
     >      (actswphelp.eq.2.0.and.(.not.pinavail))) then
c
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
c
         helpi = 17.5 + (5.0 + 37.5/t) * (1.0+ 0.25/t) * log10(1e21/n)
         src = srci(s)
c
         estphelpi = lastphelp + helpi * econv * (src-lastsrc)
c
      elseif (actswphelp.eq.2.0.or.actswphelp.eq.3.0) then
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            estphelpi = (intqe(in) * s / sptscopy(1))
         else
            estphelpi = (
     >             ( intqe(in-1) +
     >             ( (intqe(in)-intqe(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
      endif
c
      return
      end
c
c
c
      real*8 function pcxupdt(s,t,told)
      implicit none
c
c     Calculates the charge exchange loss term for each step. Stores
c     the approximated integral.
c
c     Updates the integral after new Te, Ti are found for the interval.
c     Uses the average.
c
      real*8 s,t,told
c
      include 'solparams'
      include 'solswitch'
      include 'sol22pcx'
      include 'solcommon'
c
c
c      common /pcx/ lasts,lastpcx,lastsrc
c      real*8 lasts,lastpcx,lastsrc
c
      real*8 tav
      real*8   src,srci,pinqid
      external srci,pinqid
c
      integer in
c
      tav = (t+told)/2.0
c
      if (actswpcx.eq.0.0.or.
     >   (actswpcx.eq.3.0.and.(.not.pinavail))) then
         pcxupdt = 0.0
         return
      elseif(actswpcx.eq.1.0.or.
     >      ((actswpcx.eq.2.0.or.actswpcx.eq.5.0)
     >        .and.(.not.pinavail))) then
c
         if (s.eq.soffset.or.s.lt.lasts) then
            lasts = soffset
            lastpcx = 0.0
            lastsrc = 0.0
            pcxupdt = 0.0
            return
         elseif (s.eq.lasts) then
            pcxupdt = lastpcx
            return
c
c        Approximation - if new T is below the cutoff then set addition
c                        for interval to zero.
c
         elseif (t.lt.tcutcx) then
            pcxupdt = lastpcx
            return
         endif
c
         src = srci(s)
c
         pcxupdt = lastpcx+(1.5*tav) * econv *ceicf * (src-lastsrc)
c
         lastpcx = pcxupdt
         lasts = s
         lastsrc = src
         return
c
      elseif (actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.
     >        actswpcx.eq.5) then
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            pcxupdt = (intqi(in) * s / sptscopy(1))
         else
            pcxupdt = (
     >             ( intqi(in-1) +
     >             ( (intqi(in)-intqi(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
c         if (ringnum.eq.7) then
c            write (6,'(a,i4,6(1x,g13.6))')
c     >       'I:',in,s,sptscopy(in),sptscopy(in-1),intqi(in),
c     >       intqi(in-1),pcxupdt
c         endif
c
      elseif (actswpcx.eq.4.0) then
c
         pcxupdt = pinqid(s,1)
c
      endif

      return
      end
c
c
c
      real*8 function estpcx(s,t)
      implicit none
c
c     Calculates the charge exchange loss term for each step. Stores
c     the approximated integral.
c
c     Returns an estimate of PCX over an interval. Does not update
c     overall integral.
c
      real*8 s,t
c
      include 'solparams'
      include 'solswitch'
      include 'sol22pcx'
      include 'solcommon'
c
c
c      common /pcx/ lasts,lastpcx,lastsrc
c      real*8 lasts,lastpcx,lastsrc
c
      real*8   src,srci,pinqid
      external srci,pinqid
c
      integer in
c
      if (actswpcx.eq.0.0.or.
     >   (actswpcx.eq.3.0.and.(.not.pinavail))) then
         estpcx = 0.0
         return
      elseif(actswpcx.eq.1.0.or.
     >      ((actswpcx.eq.2.0.or.actswpcx.eq.5)
     >          .and.(.not.pinavail))) then
c
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
c
         src = srci(s)
c
         estpcx = lastpcx + (1.5*t) * econv *ceicf * (src-lastsrc)
c
      elseif (actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.
     >        actswpcx.eq.5) then
c
c        Find cell which particle is in
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            estpcx = (intqi(in) * s / sptscopy(1))
         else
            estpcx = (
     >             ( intqi(in-1) +
     >             ( (intqi(in)-intqi(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
      elseif (actswpcx.eq.4.0) then
c
         estpcx = pinqid(s,0)
c
      endif
c
      return
      end
c
c
c
      real*8 function gamma(s)
      implicit none
c
c     Sets the value of the flux = nv = n0v0- int[0 to s] (S(s))
c     Goes to zero at the end of the ionization source. If using
c     a normalized ionization source.
c
      real*8 s
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 srcf,majrpos,rpos,srcrec
      external srcf,majrpos,srcrec
c
c      if (s.gt.ssrcfi) then
c         gamma = 0.0
c         return
c      endif
c
      if (actswmajr.eq.4.0) then
         rpos = majrpos(s)
         gamma = ( gamma0 * r0init + srcf(s) - srcrec(s)) / rpos
      else
         gamma = gamma0 + srcf(s) - srcrec(s)
      endif
c 
      if (debug_s22)
     >   write(6,'(a,4(1x,g12.5))') 'GAMMA:',s,gamma,gamma0

c
      return
      end
c
c
c
      real*8 function srcrec(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     SRCREC: This rotuine returns the value of the integrated
c             recombination particle source to the point s.
c
c             At this time it is only supported for PIN
c             iterated plasma calculations. It will only work
c             correctly in conjunction with ionization
c             options 2 or 8 and with perpendicular flux
c             option 2.
c
      integer in
c
      if (actswrecom.eq.0.0
     >   .or.(actswrecom.eq.1.and.(.not.pinavail))) then
         srcrec = 0.0
         return
      elseif ((pinavail.and.actswrecom.eq.1.0).or.actswrecom.eq.2) then
c
c        Search for right cell
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            srcrec =  intrecsrc(in) * s / sptscopy(1)
         else
            srcrec = intrecsrc(in-1) +
     >             ( (intrecsrc(in)-intrecsrc(in-1) ) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))
         endif
c
      endif
c
      return
      end
c
c
c
      real*8 function srcf(s)
      implicit none
c
c     This returns the integral over the source flux function
c     from 0 to s. Includes Ionization and cross-field losses.
c     Returns the values from the array intionsrc.
c
      real*8 s
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer i,in
      real*8 expsrc,trisrc,rectsrc,s5gauss,s5gauss2,gperpf
      external expsrc,trisrc,rectsrc,s5gauss,s5gauss2,gperpf
c
      srcf = 0.0
c
c     For S=0 - integration=0 so return
c
      if (s.eq.0.0) return
c
      if (actswmajr.eq.4.0) then
c
c        Major Radius corrected ion source
c
        if (s.gt.ssrcfi) then
          srcf = pnormfact * r0init
          return
        endif
c
c       ALL ionization options for major radius correction
c       have been pre-integrated because of the R(s) factor
c       in the integral - which depends on the grid. The
c       Cross-field or RCONST/GPERPCOR term has been included during
c       the pre-integration.
c
c        Search for right cell
c
         call binsearch(s,in)
c
          if (in.eq.1) then
            srcf = fnorm * (intionsrc(in) * s / sptscopy(1))
          else
            srcf = fnorm *(
     >             ( intionsrc(in-1) +
     >             ( (intionsrc(in)-intionsrc(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
          endif
c
c     Regular treatment
c
      else
c
        if (actswion.eq.0.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.0.0
     >      .and.(.not.pinavail))) then
c
           srcf = expsrc(s) + gperpf(s)
c
        elseif (actswion.eq.3.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.3.0
     >      .and.(.not.pinavail))) then
c
c          Add in triangular ionization source
c
           srcf = trisrc(s)  +  gperpf(s)
c
        elseif (actswion.eq.4.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.4.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srcf = rectsrc(s) +  gperpf(s)
c
        elseif (actswion.eq.6.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.6.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srcf = s5gauss(s) +  gperpf(s)
c
        elseif (actswion.eq.9.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.9.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srcf = s5gauss2(s) +  gperpf(s)
c
        elseif (((actswion.eq.1.0.or.actswion.eq.2.0)
     >        .and.pinavail)
     >        .or.((actswioni.eq.11.or.actswioni.eq.15)
     >        .and.(.not.pinavail))) then
c
c          Search for right cell
c
           call binsearch(s,in)
c
c
           if (in.eq.1) then
              srcf = fnorm * (intionsrc(in) * s / sptscopy(1))
     >                 + gperpf(s)
           else
              srcf = fnorm *(
     >             ( intionsrc(in-1) +
     >             ( (intionsrc(in)-intionsrc(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) )))
     >             + gperpf(s)
           endif
        else
c
c          Code has reached an error condition and should stop.
c
           write (6,'(a,2(1x,g12.5),l6)') 'ERROR in SOLASCV:'//
     >              ' Invalid Ionization Source Options: ' ,actswion,
     >                actswioni,pinavail
           stop 'SOL22: Invalid Ionizaition Source Option'
        endif
c
c     Endif for swmajr
c
      endif
c
      if (debug_s22)  
     >          write(6,'(a,3(1x,g12.5))') 'SRCF:',s,srcf
c
      return
      end
c
c
c
      real*8 function srci(s)
      implicit none
c
c     This returns the integral over ONLY the ionization source
c     from 0 to s. Includes just Ionization.
c     Returns values from the array intioniz.
c
      real*8 s
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer i,in
      real*8 expsrc,trisrc,rectsrc,s5gauss,s5gauss2
      external expsrc,trisrc,rectsrc,s5gauss,s5gauss2
c
      srci = 0.0
c
c     For S=0 - integration=0 so return
c
      if (s.eq.0.0) return
c
      if (actswmajr.eq.4.0) then
c
c       Major Radius corrected ion source
c
        if (s.gt.ssrcfi) then
          srci = pnormfact * r0init
          return
        endif
c
c       ALL ionization options for major radius correction
c       have been pre-integrated because of the R(s) factor
c       in the integral - which depends on the grid. The
c       Cross-field or RCONST/GPERPCOR term has been included during
c       the pre-integration.
c
c        Search for right cell
c
         call binsearch(s,in)
c
          if (in.eq.1) then
            srci = fnorm2 * (intioniz(in) * s / sptscopy(1))
          else
            srci = fnorm2 *(
     >             ( intioniz(in-1) +
     >             ( (intioniz(in)-intioniz(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
          endif
c
c     Regular treatment
c
      else
c
        if (actswion.eq.0.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.0.0
     >      .and.(.not.pinavail))) then
c
           srci = expsrc(s)
c
        elseif (actswion.eq.3.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.3.0
     >      .and.(.not.pinavail))) then
c
c          Add in triangular ionization source
c
           srci = trisrc(s)
c
        elseif (actswion.eq.4.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.4.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srci = rectsrc(s)
c
        elseif (actswion.eq.6.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.6.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srci = s5gauss(s)
c
        elseif (actswion.eq.9.0.or.
     >     ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >      .and.actswioni.eq.9.0
     >      .and.(.not.pinavail))) then
c
c          Add in rectangular ionization source
c
           srci = s5gauss2(s)
c
        elseif (((actswion.eq.1.0.or.actswion.eq.2.0)
     >        .and.pinavail)
     >        .or.((actswioni.eq.11.or.actswioni.eq.15)
     >        .and.(.not.pinavail))) then
c
c          Search for right cell
c
           call binsearch(s,in)
c
           if (in.eq.1) then
              srci = fnorm2 * (intioniz(in) * s / sptscopy(1))
           else
              srci = fnorm2 *(
     >             ( intioniz(in-1) +
     >             ( (intioniz(in)-intioniz(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) )))
           endif
        else
c
c         Code has reached an error condition and should stop.
c
           write (6,'(a,2(1x,g12.5),l6)') 'ERROR in SOLASCV:'//
     >              ' Invalid Ionization Source Options: ' ,actswion,
     >                actswioni,pinavail
           stop 'SOL22: Invalid Ionizaition Source Option'
        endif
c
c
c     Endif for swmajr
c
      endif
c
      return
      end
c
c
c
      real*8 function expsrc(s)
      implicit none
      real*8 s
c
c     EXPSRC: This function returns the integral of
c     the exponential ionization source from soffset to s.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c
      if (s.gt.ssrcfi) then
         expsrc = pnormfact
      else
         expsrc = s0 * (ssrcdecay*(1.0-exp(-(s-soffset)
     >                 /ssrcdecay)))
      endif
c
      return
      end
c
c
c
      real*8 function trisrc(s)
      implicit none
      real*8 s
c
c     TRISRC: This function returns the integral of
c     the triangular ionization source from soffset to s.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c
      if (s.lt.ssrcst) then
         trisrc = 0.0
      elseif (s.gt.ssrcfi) then
         trisrc = pnormfact
      elseif (s.lt.ssrcmid) then
         trisrc = s0 * ( 0.5 * (s-ssrcst)**2)
      elseif (s.le.ssrcfi) then
         trisrc = pnormfact - s0 * (0.5 * (ssrcfi-s)**2)
      endif
c
      return
      end
c
c
c
      real*8 function rectsrc(s)
      implicit none
      real*8 s
c
c     TRISRC: This function returns the integral of
c     the rectangular ionization source from soffset to s.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c
      if (s.lt.ssrcst) then
         rectsrc = 0.0
      elseif (s.gt.ssrcfi) then
         rectsrc = pnormfact
      elseif (s.le.ssrcfi) then
         rectsrc = s0 * (s-ssrcst)
      endif
c
      return
      end
c
c
c
      real*8 function s5gauss(s)
      implicit none
      real*8 s
c
c     S5GAUSS:  This function returns the integral of
c     an s**5 * exp(-s) ionization source from soffset to s.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c
      real*8 eas2,stmp
c
      stmp = s - soffset
c
      if (stmp.lt.ssrcst) then
         s5gauss = 0.0
      elseif (stmp.gt.ssrcfi) then
         s5gauss = pnormfact
      elseif (stmp.le.ssrcfi) then
         eas2 = exp(-s5alph * stmp**2)
         s5gauss = s0 * ( 1.0 / s5alph3
     >            -0.5 * stmp**4/s5alph * eas2
     >            - stmp**2/s5alph2 * eas2
     >            - eas2/s5alph3)
      endif
c
      return
      end
c
c
c
      real*8 function s5gauss2(s)
      implicit none
      real*8 s
c
c     S5GAUSS2:  This function returns the integral of
c     an s**5 * exp(-s) ionization source from soffset to s shifted
c     from the origin by a specified value (s5offset).
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c
      real*8 eas2,stmp
c
      stmp = s - soffset + s5offset
c
      if (stmp.lt.ssrcst) then
         s5gauss2 = 0.0
      elseif (stmp.gt.ssrcfi) then
         s5gauss2 = pnormfact
      elseif (stmp.le.ssrcfi) then
         eas2 = exp(-s5alph * stmp**2)
         s5gauss2 = s0 * ( s5startval
     >            -0.5 * stmp**4/s5alph * eas2
     >            - stmp**2/s5alph2 * eas2
     >            - eas2/s5alph3)
      endif
c
      return
      end
c
c
c
      subroutine preint(startn,npts,spts,src,intsrc,srcsum,ringlen,
     >                  flage2d,flagmajr,sbnd,rbnd,gperpn)
      implicit none
      include 'solparams'
      integer npts,startn
      real*8 spts(mxspts),src(mxspts),intsrc(mxspts),srcsum
      real*8 ringlen
      real*8 sbnd(0:mxspts),rbnd(0:mxspts),gperpn
      real*8 asum
      real    flage2d,flagmajr
c
c     PREINT: This subroutine pre-integrates a numerical source
c             spread along the ring.
c
c             intionsrc(S) = INTEGRAL(0 to S) of ionsrc(S) + gperp(s)
c
c             The values are supplied at each grid point represented
c             by sptscopy(i) - and are linearly interpolated.
c
c             This routine can be used to integrate any source
c             array.
c
c             Flage2d = 0.0  - do full integration
c                     = >=1.0  - Ignore first and last half cell
c
c
      integer i,ik
      real*8 rmean1,rmean2
c
c     Do integral including MAJOR RADIUS
c
c     Change to integration routine to support Major Radius
c     Not quite the same as previously since it uses actual
c     cell boundaries.
c
         srcsum = 0.0
c
         do ik = startn,npts
c
            if (flagmajr.eq.4.0) then
               rmean1 = (rbnd(ik-1)+(rbnd(ik)+rbnd(ik-1))/2.0)/2.0
               rmean2 = (rbnd(ik)+(rbnd(ik)+rbnd(ik-1))/2.0)/2.0
            else
               rmean1 = 1.0
               rmean2 = 1.0
            endif
c
            if (.not.(ik.eq.startn.and.flage2d.ne.0.0)) then
              srcsum = srcsum
     >             + (src(ik)+gperpn)
     >             * (spts(ik)-sbnd(ik-1))*rmean1
            endif
c
c           Do in two parts to get integral at cell centre.
c
            intsrc(ik) = srcsum
c
c           Does not need startn adjustment because integral is
c           only over 1/2 of the ring.
c
            if(.not.(ik.eq.npts.and.flage2d.gt.0.0)) then
               srcsum = srcsum
     >             + (src(ik)+gperpn)
     >             * (sbnd(ik)-spts(ik))*rmean2
            endif
c
         end do
c
         intsrc(npts+1) = srcsum
c
      return
      end
c
c
c
      subroutine initval
      implicit none
c
c     This subroutine calculates several of the quantities
c     in the solcommon common block ... which are used elsewhere
c     in the program.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 srcsum,momsum,srci,gamma,gtmp,pinttmp,pintupdt
      real*8 recsum,pinqid
      real*8 srcrec,srcf,gperpf
      external srci,gamma,pintupdt,pinqid,srcrec,srcf,gperpf
      real*8 areasrc(mxspts),asum
      integer i,ik
c
c     Temporary local variables
c
      real*8 tmp1,tmp2,tmp3,tmp4
c
      gperpcor = 0.0
c
c     Calculate velocity at target - including specifed mach number
c
      v0 = - m0 * sqrt ( (te0+ti0)/mb * econv/mconv)
c
      vpe2d = vpe2d * m1/e2dm0
      v1e2d = v1e2d * m1/e2dm0
      vpg   = vpg   * m1/e2dm0
c
c     Calculate pressure at target and include mach number
c
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
      elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.
     >        actswe2d.eq.9.0) then
        pstatic0 = n1 * (te1+ti1) * econv
        pinf0 = n1 * (te1+ti1) * econv + n1 * vpg**2 * mb * mconv
        pinf = pinf0
        lastvel = vpg
      endif
c
c     Additional pressure used if Velocity correction option 2 is in
c     use.
c
      padd = 0.0
c
      write (6,'(a,i4,6(1x,g13.6))') 'Pressure:',ringnum,
     >                               pstatic0,pinf0
c
      pinttmp = pintupdt(0.0d0,n0,te0,ti0)
c
c     Specify the R-value of the target -> r0init
c
      r0init = rbnd(0)
c
c
c     Calculate powers onto the target due to electrons and ions
c     modified for kinetic transfer and corrective terms.
c
      gammae = 5.0 + gamecor
c
c      if (actswe2d.eq.0.0) then
c
      gammai = 2.5 +  0.5* m0**2 * (1.0 + te0/ti0) + gamcor
c
c      else
c         gammai = 2.5 +  0.5* m0**2 * (1.0 + te1/ti1) + gamcor
c      endif
c
c      if (actswe2d.eq.0.0) then
c
c         pae = gammae * te0 * econv * n0 * abs(v0)
c         pai = gammai * ti0 * econv * n0 * abs(v0)
c
c **************************************************************
c
c      WARNING !!!!
c
c **************************************************************
c
c      There is a possibility for the calculation of pae to
c      become inconsistent if the calculation of pae_start and
c      pai_start in the calcfluxes routine performs the
c      calculations with  a mach number that is NOT equal
c      to the value used for INITM0.
c
c      This is not a problem at present - but if options
c      are created that assign initial target velocities that
c      are not equal to the sound speed then this could become
c      a problem.
c
c      D. Elder    March 6, 1998
c

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
c
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
c
         write (6,'(a,2i4,6(1x,g13.6))')
     >         'Target PAE:', ringnum,ike2d_start,pae,tmp1,gammae,
     >                     gammae*pae/tmp1
         write (6,'(a,2i4,6(1x,g13.6))')
     >         'Target PAI:', ringnum,ike2d_start,pai,tmp2,gammai,
     >                     gammai*pai/tmp2
c
      else
c
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
c
         write (6,'(a,2i4,6(1x,g13.6))')
     >         'Target PAE:', ringnum,ike2d_start,pae,tmp1,gammae,
     >                     gammae*pae/tmp1
         write (6,'(a,2i4,6(1x,g13.6))')
     >         'Target PAI:', ringnum,ike2d_start,pai,tmp2,gammai,
     >                     gammai*pai/tmp2
c
      endif
c
      write (6,'(a,i4,6(1x,g13.6))')
     >         'Target Power:', ringnum,pae,pai
c
c     Initialize the Area integral if swpow is 3.0
c
      if (actswpow.eq.3.0) then
c
c        Load up area array
c
         do ik = startn,nptscopy
            areasrc(ik) = 1.0
         end do
c
         call preint(startn,nptscopy,sptscopy,areasrc,intarea,asum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         intarea(nptscopy+1) = asum
c
         pinpute = pae * r0init / asum
         pinputi = pai * r0init / asum
c
c         write(6,*) 'Ring: ',ringnum,asum
c
c         do ik = startn,nptscopy
c            write(6,*) 'IntArea:',ik,intarea(ik)
c         end do
c
      endif
c
c     Calculated the target flux
c
      if (actswe2d.eq.1.0) then
         gamma0 = n1 * vpe2d * targfact
      elseif (actswe2d.eq.2.0) then
         gamma0 = n1 * v1e2d * targfact
      elseif (actswe2d.eq.3.0.or.actswe2d.eq.8.0.or.
     >        actswe2d.eq.9.0) then
         gamma0 = n1 * vpg * targfact
      else
         gamma0 = n0 * v0 * targfact
      endif
c
c     Initial conditions at target
c
      write (6,'(a,12g12.4)') 'Init Targ (Bound) Conditions:', 
     >                     te0,ti0,n0,v0,pinf0,initm0,n0*v0

c
c     Recombination source setup - pre-calculate the integral.
c
      recsum = 0.0
c
      if ((pinavail.and.
     >   (actswrecom.eq.1.0)).or.actswrecom.eq.2) then
c
c        Integrate over Recombination source term
c
         call preint(startn,nptscopy,sptscopy,recsrc,intrecsrc,
     >               recsum,
     >               ringlen,actswe2d,
     >               actswmajr,sbnd,rbnd,0.0d0)
c
c
         if (m0.eq.initm0) then
c
            write(6,'(a,g13.6,i4)') 'Sol option 22: recsrcint :',
     >                recsum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),recsrc(ik),
     >               intrecsrc(ik)
            end do
c
         endif
c
      endif
c
c     Calculate base value for ioniation source function
c
      call initioniz
c
c     Momentum source setup - pre-calculate the integral.
c
      if (pinavail.and.
     >   (actswnmom.eq.6.or.actswnmom.eq.7.or.
     >    actswnmom.eq.8)) then
c
c        Integrate over Momentum source term
c
         call preint(startn,nptscopy,sptscopy,momsrc,intmomsrc,
     >               momsum,
     >               ringlen,actswe2d,
     >               actswmajr,sbnd,rbnd,0.0d0)
c
c
         if (m0.eq.initm0) then

            write(6,*) 'Sol option 22: momsrcint :',momsum

            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),
     >               momsrc(ik),intmomsrc(ik)
            end do
c
         endif
c
      endif
c
c     Calculate base value for momentum loss function.
c
      if (actswnmom.eq.0.0 ) then
         smom0 = 0.0
      elseif (actswnmom.eq.1.0 ) then
         smom0 = (pinf/(actlenmom * ringlen))*(1.0d0/actffric-1.0d0)
      elseif (actswnmom.eq.2.0 ) then
         smom0 = (pinf/(lammom * ringlen))*(1.0d0/actffric-1.0d0)
     >           / ( 1.0d0 -
     >           exp( - (actlenmom * ringlen) / (lammom * ringlen)))
      elseif (actswnmom.eq.3.0 ) then
         smom0 = pinf * (1.0d0/actffric-1.0d0)
     >             * (1.0/srci(halfringlen))
c        jdemod - removed length division ... it is a bug
c         smom0 = (pinf/(actlenmom * ringlen))*(1.0d0/actffric-1.0d0)
c     >             * (1.0/srci(halfringlen))
      elseif (actswnmom.eq.4.0 ) then
         smom0 = 0.0
      elseif (actswnmom.eq.5.0 ) then
         smom0 = 0.0
      endif
c
c     Calculate base value for radiation function for PRAD option 1.
c
      if (actswprad.eq.5) then 
         prad0 = frr * (pae + pai) / (lenr-lamr)
      else
         prad0 = frr * (pae + pai) / (lamr * (1.0-exp(-lenr/lamr)))
      endif
c
c     Caculate integrated radiation loss for Prad option 4
c
      pradsum = 0.0
      call qzero(intrad,mxspts)
c
      if (actswprad.eq.4.0) then 
c
c        Integrate over Radiation source term
c
         call preint(startn,nptscopy,sptscopy,radsrc,intrad,pradsum,
     >               ringlen,actswe2d,
     >               actswmajr,sbnd,rbnd,0.0d0)
c
c        Print out radiation term for only first iteration on ring
c
         if (m0.eq.initm0) then

            write(6,'(a,1x,g13.6,i4)')
     >          'Sol option 22: radsrcint :',pradsum,ringnum
c
            do ik = startn,nptscopy
                write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),
     >               radsrc(ik),intrad(ik)
            end do
c
         endif
c
      endif
c
c
c     PIN Power Source Term SETUP - pre-calculate the integrals.
c
      call qzero(intqi,mxspts)
c
      if (pinavail.and.
     >   (actswpcx.eq.2.0.or.actswpcx.eq.3.0.or.actswpcx.eq.5)) then
c
c        Modify QISRC to remove any heating component
c
         if (actswpcx.eq.5.0) then
            do ik = startn,nptscopy
               if (qisrc(ik).lt.0.0) qisrc(ik) = 0.0
            end do
         endif
c
c        Integrate over Ion energy source term
c
         call preint(startn,nptscopy,sptscopy,qisrc,intqi,qisum,
     >               ringlen,actswe2d,
     >               actswmajr,sbnd,rbnd,0.0d0)
c
c
         if (m0.eq.initm0) then

            write(6,'(a,1x,g13.6,i4)')
     >          'Sol option 22: qisrcint :',qisum,ringnum
c
            do ik = startn,nptscopy
                write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),
     >               qisrc(ik),intqi(ik)
            end do
c
         endif
c
      elseif (actswpcx.eq.4.0) then
c
         qisum = pinqid(0.0d0,-1)
c
      else
         qisum = 0.0
      endif
c
c     Electron Source term.
c
      call qzero(intqe,mxspts)
c
      if (pinavail.and.
     >   (actswphelp.eq.2.0.or.actswphelp.eq.3.0)) then
c
c        Integrate over Electron energy source term
c
         call preint(startn,nptscopy,sptscopy,qesrc,intqe,
     >               qesum,
     >               ringlen,actswe2d,
     >               actswmajr,sbnd,rbnd,0.0d0)
c
c
         if (m0.eq.initm0) then
c
            write(6,'(a,1x,g13.6,i4)')
     >          'Sol option 22: qesrcint :',qesum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),
     >               qesrc(ik),intqe(ik)

            end do
c
         endif
c
      else
         qesum = 0.0
      endif
c
      if (m0.eq.initm0) then
c
         write (6,'(a,2(1x,g14.6))') 'Qesum, Qisum:',
     >                                     qesum,qisum
c
c        Print out the NETFlux (Gamma) function.
c
         write(6,'(a,g13.6,i4)') 'Sol option 22: GAMMA=nv :',
     >                gamma0,ringnum
         write(6,*) '  IK      S        GAMMA-ACT    GAMMA-CALC '
     >        //'   GAMMA0   SRCF-GPERPF     GPERPF   '
     >        //'     REC          SRCF'
c
         do ik = startn,nptscopy
             gtmp = gamma(sptscopy(ik))
             tmp1 = srcf(sptscopy(ik))
             tmp2 = gperpf(sptscopy(ik))
             tmp3 = srcrec(sptscopy(ik))
             tmp4 = gamma0
             write(6,'(i4,8(1x,g12.5))') ik,sptscopy(ik),
     >               gtmp, tmp4 + tmp1 - tmp3,
     >               tmp4,tmp1-tmp2,tmp2,tmp3,tmp1
         end do
c
      endif

c
c slmod begin - new
      CALL SOL22Headers
c slmod end
      return
      end
c
c
c
      subroutine initioniz
      implicit none
c
c     INITIONIZ: This subroutine initializes the
c                ionization source depending on
c                specified options and switches.
c
c                Among other things - all sources
c                that involve major radius
c                correction option 4 are pre-integrated
c                since it would be too costly to
c                recalculate the integrals at
c                each step. These are then linearly
c                interpolated in the srci function.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 srcsum,momsum,srci,gamma,gtmp,pinttmp,pintupdt
      real*8 srcf,gperpf,srcrec
      real*8 tmp1,tmp2,tmp3,tmp4 
      external srci,gamma,pintupdt,srcf,gperpf,srcrec
      integer i,ik
c
      real*8 majrpos
      external majrpos
c
      real*8 tmpgamma,tmpsrci,tmprpos,tmpint(mxspts)
      real*8 tmpsrcsum,ends,eas1,eas2
c
c     IPP/08 Krieger - gtmp and srcsum should be initialized to avoid
c     run time errors because they are included in write statements but
c     not used in certain cases
c
      gtmp=0.0
      srcsum=0.0
c
c     Set the zero offset for Edge2D compatibility cases

      if (actswe2d.ne.0.0) then
         soffset = sptscopy(ike2d_start)
      else
         soffset = 0.0d0
      endif
c
c     Check the ionization source position against the offset.
c
      if (ssrcst.lt.soffset) then
         ssrcst = soffset
         if (ssrcfi.lt.ssrcst) then
            write (7,*) 'ERROR ERROR ERROR - Invalid'//
     >                  ' Ionization Source Specification'
            ssrcfi = ssrcst + 1.0
         endif
c
         ssrcmid = (ssrcfi+ssrcst) / 2.0
         ssrclen = (ssrcfi-ssrcst)
c
c         s5gausslen = s5gausslen - soffset
c
      endif
c
c     Check and calculate the integral of the PIN ionization source
c     without ANY cross-field terms - for use later in normalizing
c     the analytic options.
c
      if (pinnorm.eq.1) then
c
         call preint(startn,nptscopy,sptscopy,ionsrc,tmpint,
     >               tmpsrcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         pnormfact = tmpsrcsum
c
      else
         pnormfact = -gamma0
      endif
c
c     GPERP correction is only done when the
c     Major radius corrector is off - the major radius
c     corrector factors in GPERPCOR options.
c
c
c     Calculate the gamma perp correction if the option is turned
c     ON - otherwise set it to zero. Note that Gamma Perp correction
c     will only apply to non-normalised ionization - otherwise it will
c     (by definition) be equal to zero.
c
      if (actswgperp.eq.0.0) then
c
         gperpcor = 0.0
c
      elseif (actswgperp.eq.1.0) then
c
         gperpcor = 0.0
c
c        Need to pre-calculate the ionization source
c        using gperpcor = 0 - then find the difference
c        and finally recalculate the correct gperpcor and
c        then re-do the source pre-integration. This involves
c        repeating the exact code from below twice - though
c        only for this option. See below for the
c        commented code.
c
c        Only needs to be done for un-normalized sources.
c        i.e. at present (pinavail.and.actswion.eq.2)
c
c        Also need to include this calculation
c        if RECOMBINATION is ON - and a normalized option has been
c        specified.
c
         if ((pinavail.and.(actswion.eq.2.0.or.
     >       actswrecom.eq.1.0.or.actswrecom.eq.2)).or.
     >       (.not.pinavail
     >        .and.(actswioni.eq.11.or.actswioni.eq.15))
     >       ) then
c
c          Perform source preintegration - include the R(s) factor
c          if required.
c          Also include perpendicular loss corrections.
c
c          Set flag passed to preint routine.
c
           call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)

c
c          Set the normalization factor
c
           fnorm = 1.0
c
c          Use the above to calculate the gperpcor for the 1/2 ring
c
           gtmp = gamma(halfringlen)
c
           if (actswe2d.ne.0.0) then
               gperpcor = -gtmp/(halfringlen-soffset)
           else
               gperpcor = -gtmp/(halfringlen)
           endif
c
         endif
c
         write (6,'(a,5g18.7)') 'Gperpcor:1a:',gperpcor,gtmp,
     >                           srcsum,0.5*ringlen,halfringlen
c
      elseif (actswgperp.eq.2.0) then
c
c        GNET is calculated in the routine calcsoliz called from
c        the beginning of CALCSOL_INTERFACE. It uses the fluxes
c        at each target and the ionization over the whole ring
c        to calculate a correction that is distributed over the
c        whole ring.
c
c        Again this is only non-zero for non-normalized sources.
c
c
         if ((pinavail.or.(.not.pinavail
     >        .and.(actswioni.eq.11.or.actswioni.eq.15)))
     >        .and.
     >      (actswion.eq.2.0.or.pinnorm.eq.1)) then
            gperpcor = gnet
            write (6,'(a,2g18.7)') 'Gperpcor:2a:',gperpcor,gnet
         else
            gperpcor = 0.0
         endif
c
         write (6,'(a,2g18.7)') 'Gperpcor:2b:',gperpcor,gnet
c
      elseif (actswgperp.eq.3.0) then
c
c        Need to preint the gperp array so as to get the integrated
c        gperp contribution.
c
         call preint(startn,nptscopy,sptscopy,gperp,intgperp,srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         gperpcor = srcsum
c
         if (m0.eq.initm0) then
c
            write(6,'(a,g13.6,i4)') 'Sol option 22: GPERPsrcint :',
     >                srcsum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),
     >               intgperp(ik)
            end do
         endif
c
      elseif (actswgperp.eq.4.0) then
c
c        This is the same as option 1 if the 1/2 flux tube is under
c        ionized and will distribute the X-field flux proportional
c        to density (as in option 3) if the 1/2 flux tube is
c        over-ionized.
c
         gperpcor = 0.0
c
c        Need to pre-calculate the ionization source
c        using gperpcor = 0 - then find the difference
c        and finally recalculate the correct gperpcor and
c        then re-do the source pre-integration. This involves
c        repeating the exact code from below twice - though
c        only for this option. See below for the
c        commented code.
c
c        Only needs to be done for un-normalized sources.
c        i.e. at present (pinavail.and.actswion.eq.2)
c
c        Also need to include this calculation
c        if RECOMBINATION is ON
c
         if ((pinavail.and.(actswion.eq.2.0.or.
     >        actswrecom.eq.1.0.or.actswrecom.eq.2)).or.
     >       (.not.pinavail
     >        .and.(actswioni.eq.11.or.actswioni.eq.15))
     >       ) then
c
c          Perform source preintegration - include the R(s) factor
c          if required.
c          Also include perpendicular loss corrections.
c
c          Set flag passed to preint routine.
c
           call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
c          Set the normalization factor
c
           fnorm = 1.0
c
c          Use the above to calculate the gperpcor for the 1/2 ring
c
           gtmp = gamma(halfringlen)
c
           if (actswe2d.ne.0.0) then
               gperpcor = -gtmp/(halfringlen-soffset)
           else
               gperpcor = -gtmp/(halfringlen)
           endif
c
c          Need to re-distribute the flux - for the case where
c          the flux-tube is over-ionized - i.e. gtmp > 0
c
           if (gperpcor.ge.0.0) then
c
c             Reset the Gperp option to be the same as 1 - since
c             it is functionally identical from here on.
c
              actswgperp = 1.0
c
           elseif (gperpcor.lt.0.0.and.m0.eq.initm0) then
c
c             Redistribute the integrated flux in proportion to the
c             density in each cell.
c
c             Change back to total required from amount/unit length
c
              gperpcor = gperpcor * (halfringlen)
c
              call preint(startn,nptscopy,sptscopy,oldne,tmpint,
     >               tmpsrcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
c             Copy the density integral and renormalize it so that
c             the total is equal to the required gperpcor.
c
              do ik = startn,nptscopy
                 gperp(ik) = gperpcor/tmpsrcsum * oldne(ik)
                 intgperp(ik) = gperpcor/tmpsrcsum * tmpint(ik)
              end do

              write(6,'(a,g13.6,i4,f7.3,2g13.6)')
     >                'Sol option 22: GPERPsrcint :',
     >                tmpsrcsum,ringnum,actswgperp,gperpcor,gtmp
c
              do ik = startn,nptscopy
                 write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),
     >               intgperp(ik)
              end do
c
           endif
c
         endif
c
      elseif (actswgperp.eq.5.0) then
c
         gperpcor = 0.0
c
c        Need to pre-calculate the ionization source
c        using gperpcor = 0 - then find the difference
c        and finally recalculate the correct gperpcor and
c        then re-do the source pre-integration. This involves
c        repeating the exact code from below twice - though
c        only for this option. See below for the
c        commented code.
c
c        Only needs to be done for un-normalized sources.
c        i.e. at present (pinavail.and.actswion.eq.2)
c
c        Also need to include this calculation
c        if RECOMBINATION is ON - and a normalized option has been
c        specified.
c
         if ((pinavail.and.(actswion.eq.2.0.or.
     >       actswrecom.eq.1.0.or.actswrecom.eq.2)).or.
     >       (.not.pinavail
     >        .and.(actswioni.eq.11.or.actswioni.eq.15))
     >       ) then
c
c          Perform source preintegration - include the R(s) factor
c          if required.
c          Also include perpendicular loss corrections.
c
c          Set flag passed to preint routine.
c
           call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
c          Set the normalization factor
c
           fnorm = 1.0
c
c          Use the above to calculate the gperpcor for the 1/2 ring
c
           gtmp = gamma(halfringlen)
c
           if (actswe2d.ne.0.0) then
               gperpcor = -gtmp/(halfringlen-soffset)
           else
               gperpcor = -gtmp/(halfringlen)
           endif
c
c          Now that the total is available - distribute it
c          between the two sources - uniform + rectangular.
c
c          Square
c
           gperpcor2 = gperpcor * gperpfrac
     >                 * halfringlen /(sgperpend-sgperpbeg)
c
c          Reassign gperpcor
c
c          Uniform
c
           gperpcor = gperpcor * (1.0-gperpfrac)
c
         endif
c

      elseif (actswgperp.eq.6.0) then
c
c        GNET is calculated in the routine calcsoliz called from
c        the beginning of CALCSOL_INTERFACE. It uses the fluxes
c        at each target and the ionization over the whole ring
c        to calculate a correction that is distributed over the
c        whole ring.
c
c        Again this is only non-zero for non-normalized sources.
c
c        Cross-field source is split into two components
c
c
         if ((pinavail.or.(.not.pinavail
     >        .and.(actswioni.eq.11.or.actswioni.eq.15)))
     >        .and.
     >      (actswion.eq.2.0.or.pinnorm.eq.1)) then
            gperpcor = gnet
         else
            gperpcor = 0.0
         endif
c
c        Now that the total is available - distribute it
c        between the two sources - uniform + rectangular.
c
c        Square
c
         gperpcor2 = gperpcor * gperpfrac
     >               * halfringlen /(sgperpend-sgperpbeg)
c
c        Reassign gperpcor
c
c        Uniform
c
         gperpcor = gperpcor * (1.0-gperpfrac)
c
         write(6,'(a,6(1x,g13.6))') 'gperpcor:',
     >                          gnet,gperpcor,gperpcor2,
     >                          ringlen,sgperpend,sgperpbeg

c
      elseif (actswgperp.eq.7.0.or.actswgperp.eq.8.0) then
c
c
c        Need to preint the gperp array so as to get the integrated
c        gperp contribution.
c
         call preint(startn,nptscopy,sptscopy,gperp,intgperp,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         gperpcor = srcsum
c
         if (m0.eq.initm0) then
c
            write(6,'(a,3g13.6,i4)') 'Sol option 22: GPERPsrcint :',
     >                srcsum,gnet*ringlen,srcsum/(gnet*ringlen),ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),gperp(ik),
     >               intgperp(ik)
            end do
         endif
c
      endif
c
c     Do code for major radius corrected ionization source.
c
      if (actswmajr.eq.4.0) then
c
c        For an analytic source - set up the ionsrc array
c        with the appropriate values and treat numerically
c        from this point forward.
c
c        For swioni = 11 - this is already done.
c
         if (actswion.eq.0.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.0.0
     >         .and.(.not.pinavail))) then
c
c          Watch OUT for S-offsets when assigning this - check all
c          options and fix bugs.
c
c          Approximate s0 for now - it is not accurate since the
c          integral over S is not the simple exponential
c          given when the R-correction is applied. However, the
c          incorrect numbers will be fixed by the fnorm
c          factor which requires that the integral over
c          the numerical source be equal to the target flux
c          for the initial prescription cases.
c
c
c           if (actswe2d.ne.0.0) then
c             soffset = sptscopy(1)
c           else
c             soffset = 0.0d0
c           endif
c
           s0 = pnormfact *r0init /(ssrcdecay*
     >                    (1.0-exp(-(ssrcfi-soffset)/ssrcdecay)))
c
           do ik = startn,nptscopy
c
              if (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0*exp(-(sptscopy(ik)-soffset)/ssrcdecay)
              else
                ionsrc(ik) = 0.0
              endif
c
           end do
c
         elseif (actswion.eq.3.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.3.0
     >         .and.(.not.pinavail))) then
c
c          Add in triangular source option
c
           s0 = pnormfact * r0init / ssrclen**2
c
           do ik = startn,nptscopy
c
              if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).lt.ssrcmid) then
                ionsrc(ik) = s0 * (sptscopy(ik) - ssrcst)
              elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * (ssrcfi - sptscopy(ik))
              endif
c
           end do
c
         elseif (actswion.eq.4.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.4.0
     >         .and.(.not.pinavail))) then
c
c          Add in triangular source option
c
           s0 = pnormfact * r0init / ssrclen
c
           do ik = startn,nptscopy
c
              if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0
              endif
c
           end do
c
         elseif (actswion.eq.6.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.6.0
     >         .and.(.not.pinavail))) then
c
c          Add in s5gauss source option
c
           s5alph = 2.5 / s5gausslen ** 2
           s0 = pnormfact * s5alph**3
c
           do ik = startn,nptscopy
c
              if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * sptscopy(ik)**5
     >                      * exp(-s5alph*sptscopy(ik)**2)
              endif
c
           end do
c
         elseif (actswion.eq.9.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.9.0
     >         .and.(.not.pinavail))) then
c
c          Add in shifted S5*gaussian source option
c
           s5alph = 2.5 / s5gausslen ** 2
           s0 = pnormfact * s5alph**3
c
           do ik = startn,nptscopy
c
              if (sptscopy(ik).lt.ssrcst) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).gt.ssrcfi) then
                ionsrc(ik) = 0.0
              elseif (sptscopy(ik).le.ssrcfi) then
                ionsrc(ik) = s0 * sptscopy(ik)**5
     >                  * exp(-s5alph*(sptscopy(ik)+s5gausslen/2.0)**2)
              endif
c
           end do
c
         endif
c
c        Perform source preintegration - include the R(s) factor
c        Set flag passed to preint routine. Include Perpendicular
c        flux correction factor.
c
         call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         if (m0.eq.initm0) then
c
            write(6,'(a,g13.6,i4)') 'Sol option 22: FLUXsrcint :',
     >                srcsum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),
     >               ionsrc(ik),intionsrc(ik),gperpf(sptscopy(ik))
            end do
         endif
c
c        Set the normalization factor
c
         if (pinavail.and.actswion.eq.2.0) then
            fnorm = 1.0
         elseif (pinnorm.eq.1) then
            fnorm = 1.0
         else
            fnorm = -gamma0 * r0init / srcsum
         endif
c
c
c
c        Calculate the pre-integral of JUST the ionization
c        by itself - for use in various source terms.
c
         call preint(startn,nptscopy,sptscopy,ionsrc,intioniz,srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         if (m0.eq.initm0) then
c
            write(6,'(a,g13.6,i4)') 'Sol option 22: IONsrcint :',
     >                srcsum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i4,3(1x,g13.6))') ik,sptscopy(ik),ionsrc(ik),
     >               intioniz(ik)
            end do
         endif
c
         if (pinavail.and.actswion.eq.2.0) then
            fnorm2 = 1.0
         elseif (pinnorm.eq.1) then
            fnorm2 = pnormfact/srcsum
         else
            fnorm2 = -gamma0 * r0init / srcsum
         endif
c
c      Normal options for when the major radius option is not
c      turned on.
c
       else
c
c
c
c        Exponential decay options
c
         if (actswion.eq.0.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.0.0
     >         .and.(.not.pinavail))) then
c
c
            s0 = pnormfact / (ssrcdecay*(1.0-exp(-(ssrcfi-soffset)
     >                                       /ssrcdecay)))
c
         elseif (actswion.eq.3.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.3.0
     >         .and.(.not.pinavail))) then
c
c           Add code for triangular source option.
c
            s0 = 4.0 * pnormfact / ssrclen**2
c
         elseif (actswion.eq.4.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.4.0
     >         .and.(.not.pinavail))) then
c
c           Add code for triangular source option.
c
            s0 = pnormfact / ssrclen
c
         elseif (actswion.eq.6.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.6.0
     >         .and.(.not.pinavail))) then
c
c           Add code for s5 * gaussian source option.
c
            s5alph = 2.5 / s5gausslen **2
            s5alph2 = s5alph**2
            s5alph3 = s5alph**3
            s0 = pnormfact * s5alph3
c
         elseif (actswion.eq.9.0.or.
     >        ((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >         .and.actswioni.eq.9.0
     >         .and.(.not.pinavail))) then
c
c           Add code for s5 * gaussian source option.
c
            s5alph = 2.5 / s5gausslen **2
            s5offset = s5gausslen/2.0
            s5alph2 = s5alph**2
            s5alph3 = s5alph**3
c
            ends = ssrcfi + s5offset
            eas2 = exp(-s5alph*ends**2)
            eas1 = exp(-s5alph*s5offset**2)
c
            s0 = pnormfact /
     >           ( (-ends**4/(2.0*s5alph)*eas2
     >              -ends**2/s5alph2*eas2
     >              -eas2/s5alph3)
     >           - (-s5offset**4/(2.0*s5alph)*eas1
     >              -s5offset**2/s5alph2*eas1
     >              - eas1/s5alph3))
c
            s5startval = -1.0 *
     >           (-s5offset**4/(2.0*s5alph)*eas1
     >            -s5offset**2/s5alph2*eas1
     >            -eas1/s5alph3)
c
            write(6,'(a,7g13.5)') 'Init9:',s0,s5offset,s5startval,
     >                             s5alph,gamma(halfringlen),
     >            (-ends**4/(2.0*s5alph)*eas2
     >              -ends**2/s5alph2*eas2
     >              -eas2/s5alph3),
     >           - (-s5offset**4/(2.0*s5alph)*eas1
     >              -s5offset**2/s5alph2*eas1
     >              - eas1/s5alph3)
c
         endif
c
         write (6,'(a,9g12.4)') 'InitI:',actswion,actswioni,
     >              actswgperp,s0,ssrclen,gperpcor,pnormfact,
     >                        gamma0
         write (6,'(a,6g12.4)') 'LensI:',ssrcst,ssrcfi,ssrclen,
     >                           ssrcmid,soffset,s5gausslen
c
c        Pre-calculate the integral of the source if necessary.
c
         if ((pinavail.and.
     >      (actswion.eq.1.0.or.actswion.eq.2.0)).or.
     >      ((actswioni.eq.11.or.actswioni.eq.15)
     >      .and.(.not.pinavail))) then
c
            call preint(startn,nptscopy,sptscopy,ionsrc,intionsrc,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
c           Set up ionization source normalization if it is required.
c
c
            if (actswion.eq.1.0) then
              if (srcsum.ne.0.0) then
                 fnorm = -gamma0/srcsum
              else
                 fnorm = 1.0
              endif
            elseif (actswion.eq.2.0) then
               fnorm = 1.0
            endif
c
c           Print out the source
c
            if (m0.eq.initm0) then
c
               write(6,'(a,g13.6,2i4)') 'Sol option 22: FLUXsrcint :',
     >                srcsum,ringnum,nptscopy
c
               do ik = startn,nptscopy
                  write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),
     >               ionsrc(ik),
     >               intionsrc(ik),srcf(sptscopy(ik))
               end do
c
            endif
c
c
c           Calculate the pre-integral of JUST the ionization
c           by itself - for use in various source terms.
c
            call preint(startn,nptscopy,sptscopy,ionsrc,intioniz,
     >               srcsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
c           Set source normalization
c
            if (actswion.eq.1.0) then
              if (srcsum.ne.0.0) then
                 fnorm2 = -gamma0/srcsum
              else
                 fnorm2 = 1.0
              endif
            elseif (actswion.eq.2) then
               fnorm2 = 1.0
            endif
c
c           Print out the source
c
            if (m0.eq.initm0) then
               write(6,'(a,g13.6,i4)') 'Sol option 22: IONsrcint :',
     >                srcsum,ringnum
c
               do ik = startn,nptscopy
                  write(6,'(i4,4(1x,g13.6))') ik,sptscopy(ik),
     >               ionsrc(ik),
     >               intioniz(ik),srci(sptscopy(ik))
               end do
            endif
c
         endif
c
c
c     ENDIF of swmajr
c
      endif
c
      return
      end
