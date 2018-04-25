c
c ======================================================================
c
c subroutine: ThomPP
c
      SUBROUTINE ThomPP(irlim1,irlim2,ikopt,targ_con_opt)
      IMPLICIT none
c
c     irlim1 - start ring
c     irlim2 - end ring
c     ikopt  - apply to whole or partial ring
c     targ_con_opt   =3 - do not apply conditions to target
c                    =4 - apply average values to target
c
c
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
c
c      INCLUDE 'solparams'
c      INCLUDE 'solswitch'
c      INCLUDE 'solcommon'
c

      INTEGER irlim1,irlim2,ikopt,targ_con_opt
      integer irstart,irend,ikstart,ikend

      REAL GetCs

      INTEGER    MAXTDAT     ,MAXCOLS
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)

      INTEGER thnnraw
c slmod begin - new
      REAL    thnraw(MAXTDAT,MAXCOLS)
c
c      REAL    dummy(MAXNKS,MAXNRS),thnraw(MAXTDAT,MAXCOLS)
c slmod end

      INTEGER ik,ir,ir1,ir2,in,tag(MAXNRS),num
      LOGICAL status,output,loaded
      REAL    te(2,MAXNRS),ne(2,MAXNRS),d1,d2,f
c slmod begin - new
      DATA loaded /.FALSE./

      SAVE
c slmod end

      output = .TRUE.
      status = .FALSE.
c slmod begin - new
      loaded = .false.
c slmod end

      CALL IZero(tag,MAXNRS)
      CALL RZero(te ,MAXNRS*2)
      CALL RZero(ne ,MAXNRS*2)
c
      if (.not.loaded) then
c slmod begin - new
         CALL LoadThomsonData(thnnraw,thnraw,MAXTDAT,MAXCOLS,NIL,NIL,2)
c
c         CALL LoadThomsonData(dummy,thnnraw,thnraw,MAXTDAT,MAXCOLS,2,2)
c slmod end
         loaded = .true.

         IF (output) THEN
            DO in = 1, thnnraw
               WRITE(PINOUT,'(A,I6,4F10.4,1P,E10.2,0P,F10.4)')
     .               'THOMSON : ', in,(thnraw(in,ir),ir=1,6)
            ENDDO
         ENDIF

c
c...     Find average temperature and density for rings with Thomson data:
c
         DO ir = 1, nrs

            num = 0
            DO in = 1, thnnraw
               IF (INT(thnraw(in,4)).EQ.ir) THEN
                  num = num + 1

                  f = 1.0 / REAL(num)

                  te(IKLO,ir) =  f * thnraw(in,6)
     .                        + (1.0 - f) * te(IKLO,ir)
                  ne(IKLO,ir) =  f * thnraw(in,5)
     .                        + (1.0 - f) * ne(IKLO,ir)
                  te(IKHI,ir) = te(IKLO,ir)
                  ne(IKHI,ir) = ne(IKLO,ir)

                  if (ir.ge.irtrap.and.ir.le.nrs) status = .TRUE.

                  IF (output)
     .               WRITE(PINOUT,'(A,2I6,F6.3,2F10.4,1P,2E10.2,0P)')
     .                    ' A : ',
     .                    ir,num,f,thnraw(in,6),te(IKLO,ir),
     .                             thnraw(in,5),ne(IKLO,ir)
               ENDIF
            ENDDO

            IF (num.GT.0) tag(ir) = 1

            IF (output)
     .              WRITE(PINOUT,'(A,3I6,F10.4,1P,E10.2,0P)') ' B : ',
     .                ir,tag(ir),num,te(IKLO,ir),ne(IKLO,ir)
         ENDDO
c


c...     Make sure there is some Thomson data in the PP:
         IF (ir.ge.irtrap.and.ir.le.nrs.and.
     .      (.NOT.status)) CALL ER('ThomPP','No data in PP',*99)

c...     Assign values to rings without Thomson data:
         DO ir = 1, nrs

           IF (ir.ge.irtrap.and.ir.le.nrs) then
              irstart = irtrap
              irend = nrs
           elseif (ir.ge.1.and.ir.le.irwall) then
              irstart = 1
              irend = irwall
           endif
c
           IF (tag(ir).EQ.0) THEN
             ir1 = ir
c slmod begin - new
c...fix
             DO WHILE (ir1.GT.irstart.AND.tag(ir1).EQ.0)
               ir1 = ir1 - 1
             ENDDO
             ir2 = ir
             DO WHILE (ir2.LT.irend+1.AND.tag(ir2).EQ.0)
               ir2 = ir2 + 1
             ENDDO
c
c             DO WHILE (ir.GT.irstart.AND.tag(ir1).EQ.0)
c               ir1 = ir1 - 1
c             ENDDO
c             ir2 = ir
c             DO WHILE (ir.LT.irend+1.AND.tag(ir2).EQ.0)
c               ir2 = ir2 + 1
c             ENDDO
c slmod end
c
             if (ir1.eq.irstart.and.ir2.eq.irend+1) then
               te(IKLO,ir) = 1.0
               ne(IKLO,ir) = 1.0e10
               te(IKHI,ir) = te(IKLO,ir)
               ne(IKHI,ir) = ne(IKLO,ir)
             elseif (ir1.EQ.irstart) THEN
               te(IKLO,ir) = te(IKLO,ir2)
               ne(IKLO,ir) = ne(IKLO,ir2)
               te(IKHI,ir) = te(IKLO,ir)
               ne(IKHI,ir) = ne(IKLO,ir)
             ELSEIF (ir2.EQ.irend+1) THEN
               te(IKLO,ir) = te(IKLO,ir1)
               ne(IKLO,ir) = ne(IKLO,ir1)
               te(IKHI,ir) = te(IKLO,ir)
               ne(IKHI,ir) = ne(IKLO,ir)
             ELSE
               d1 = REAL(ir  - ir1)
               d2 = REAL(ir2 - ir )
               te(IKLO,ir) = (d2*te(IKLO,ir1)
     .                     + d1*te(IKLO,ir2)) / (d1+d2)
               ne(IKLO,ir) = (d2*ne(IKLO,ir1)
     .                     + d1*ne(IKLO,ir2)) / (d1+d2)
               te(IKHI,ir) = te(IKLO,ir)
               ne(IKHI,ir) = ne(IKLO,ir)
             ENDIF
c
             IF (output)
     .           WRITE(PINOUT,'(A,3I6,F10.4,1P,E10.2,0P)') ' C : ',
     .                  ir,ir1,ir2,te(IKLO,ir),ne(IKLO,ir)
           ENDIF
c
         ENDDO

c
      endif
c

c
c...  Assign plasma to rings with Thomson data:
c
      DO ir = irlim1, irlim2
c
c        IF (ir.LT.irtrap.OR.ir.GT.nrs) CYCLE
c
        call set_ikvals(ir,ikstart,ikend,ikopt)
c
        DO ik = ikstart, ikend
          ktebs(ik,ir) = te(IKLO,ir)
          ktibs(ik,ir) = te(IKLO,ir)
          knbs (ik,ir) = ne(IKLO,ir)
          kvhs (ik,ir) = 0.0
        ENDDO

c...    Assign target values using Thomson data:
        IF (targ_con_opt.eq.4) THEN
c
           if (ikopt.eq.2.or.ikopt.eq.3) then

              kteds(idds(ir,1)) = ktebs(nks(ir),ir)
              ktids(idds(ir,1)) = ktibs(nks(ir),ir)
              knds (idds(ir,1)) = knbs (nks(ir),ir)
              kvds (idds(ir,1)) =
     .                GetCs(kteds(idds(ir,1)),ktids(idds(ir,1)))

           endif

           if (ikopt.eq.1.or.ikopt.eq.3) then

              kteds(idds(ir,2)) = ktebs(1      ,ir)
              ktids(idds(ir,2)) = ktibs(1      ,ir)
              knds (idds(ir,2)) = knbs (1      ,ir)
              kvds (idds(ir,2)) =
     .                -GetCs(kteds(idds(ir,2)),ktids(idds(ir,2)))
           endif


          ENDIF
c
      ENDDO

      RETURN
99    STOP
      END
c
c
c
      subroutine calcsol_interface (irlim1,irlim2,ikopt)
      use error_handling
      use debug_options
      use sol22_input
      use sol22_debug
      implicit none
      integer irlim1, irlim2,ikopt
c
c     The subroutine calcsol uses numerical methods (Runge-Kutta) to
c     solve the fluid equations along the field lines. This routine
c     acts as an interface between DIVIMP and the calcsol subroutine.
c     This approach was chosen so that the variables and components
c     of the solascv module could remain as independent of DIVIMP as
c     possible - thus facilitating changes and corrections to the code
c     based on developments in the stand-alone program module solascv.f
c     ... in particular ... at this time the calcsol series of
c     subroutines calculate all values in extended precision for
c     accuracy - which is not generally wanted or needed inside DIVIMP.
c     In addition, the calcsol routines have their own variable names
c     for many quantities already included in DIVIMP - in order to avoid
c     re-writing the code at this time - the interface routine assures
c     that all the quantities required in the calcsol subroutine series
c     are properly loaded.
c
c     David Elder,  Jan 24, 1995
c
c
c
      include 'params'
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
c slmod begin - new
      INCLUDE 'slcom'
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c

c
c
      COMMON /POWERFLOW/ cve        ,cvi        ,cde        ,cdi
      REAL               cve(MXSPTS),cvi(MXSPTS),cde(MXSPTS),cdi(MXSPTS)
      
      REAL Clock2

      REAL deltat
c slmod end
c
      integer ir,ircor,midnks,ik,ikstart,ikend
c
      real gperpa(maxnks,maxnrs)
      real oldknbs(maxnks,maxnrs)
      real grad_oldknbs(maxnks,maxnrs)
      real oldktibs(maxnks,maxnrs)
      real oldktebs(maxnks,maxnrs)
      real oldkvds(maxnds)
c
c     Radiation loss term from previous DIVIMP run
c     Impurity Ionization/recombination potential energy loss to e- 
c
      real div_tpowls(maxnks,maxnrs),div_tcooliz(maxnks,maxnrs)
      real div_cool(maxnks,maxnrs) 
c
      real gradi,grado,flux,len1,len2
      real quant2grad,polysidelen
      external quant2grad,polysidelen
      integer rc
c
      integer errcode,seterror,new_errlevel
      real*8 serr
c
c     For calling calcsol
c
      integer npts,id
      real*8 spts(mxspts),fact,slope,gam
      real*8 te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts),
     >        exp_press(mxspts),act_press(mxspts),
     >        prad(mxspts)
      real*8 ttarg
      real*8 int_powrat(3)
c
c     Assigning ACTFFRIC
c
      real*8 find_ffric 
      external find_ffric
c
c      integer ierr
c
c     Local variables
c
      real ds1,dp1,dt1,nb1,ds2,dp2,dt2,nb2
      integer applydef(maxnrs,2),in,pplasma,ikfirst,iklast
      integer sol22_iter
      data sol22_iter /0/
c
c     Temporary local variables until these are loaded in input file
c
      integer:: pfz_dist_opt
      real*8 :: pfz_dist_param(2)

c
c     TEMPORARY MULTIPLIER
c
c      real qesrc_mult
c      parameter (qesrc_mult=1.0)
c
c     Locals for smoothing
c
      real temid,timid,nmid,smax,asmexp,tmp
      
      call pr_trace('SOLASCV1','START CALCSOL_INTERFACE')

c
c     Make a call to initialize debugging if it is on
c
      if (debug_sol22.ne.0) then 
         ! debug_sol22 is manually set in the sol22_debug module
         ! parameters to the call are the ring number and ikopt for the half ring for which high res
         ! debugging data is required. 
         call init_sol22_debug(debug_sol22_ir,debug_sol22_ikopt) 
      endif

c
c     Set a flag internal to the SOLASCV code that indicates that
c     PIN data is available if the ionization option for it has been
c     specified.
c
      pinavail = piniter
c
c
c     If Switch(swe2d).lt.0 - use this to indicate using the
c     Edge2D background plasma as the seed to PIN for the
c     SOL option 22 iteration - flip the switch "+" for
c     subsequent iterations.
c
      if (.not.pinavail) then
c
         e2dstart = 0
c
         if (switch(swe2d).lt.0.0.or.
     >       switch(swioni).eq.12.or.
     >       switch(swioni).eq.13.or.
     >       switch(swioni).eq.14) then
c
            sol22_iter = sol22_iter + 1
c
            switch(swe2d) = abs(switch(swe2d))
            e2dstart = 1
c
            if (cre2d.eq.0) then
               call redge2d(0)
c               cre2d = 1
            endif
c
            if (switch(swioni).eq.14) then
               do ir = 1,nrs
                  do ik = 1,nks(ir)
                     oldknbs(ik,ir) = e2dnbs(ik,ir)
                     oldktibs(ik,ir) = e2dtibs(ik,ir)
                     oldktebs(ik,ir) = e2dtebs(ik,ir)
                     knbs(ik,ir) = e2dnbs(ik,ir)
                     ktebs(ik,ir)= e2dtebs(ik,ir)
                     ktibs(ik,ir)= e2dtibs(ik,ir)
                     kvhs(ik,ir) = e2dvhs(ik,ir)
                     kes(ik,ir)  = e2des(ik,ir)
                  end do
c
                  if (ir.lt.irsep) then
                     oldknbs(nks(ir),ir) = e2dnbs(1,ir)
                     oldktibs(nks(ir),ir) = e2dtibs(1,ir)
                     oldktebs(nks(ir),ir) = e2dtebs(1,ir)
                     knbs(nks(ir),ir) = e2dnbs(1,ir)
                     ktebs(nks(ir),ir)= e2dtebs(1,ir)
                     ktibs(nks(ir),ir)= e2dtibs(1,ir)
                     kvhs(nks(ir),ir) = e2dvhs(1,ir)
                     kes(nks(ir),ir)  = e2des(1,ir)
                  endif
c
               end do
c
            else
c
               do ir = irlim1,irlim2
c
                  call set_ikvals(ir,ikstart,ikend,ikopt)
c
                  do ik = ikstart,ikend
                     oldknbs(ik,ir) = e2dnbs(ik,ir)
                     oldktibs(ik,ir) = e2dtibs(ik,ir)
                     oldktebs(ik,ir) = e2dtebs(ik,ir)
                     knbs(ik,ir) = e2dnbs(ik,ir)
                     ktebs(ik,ir)= e2dtebs(ik,ir)
                     ktibs(ik,ir)= e2dtibs(ik,ir)
                     kvhs(ik,ir) = e2dvhs(ik,ir)
                     kes(ik,ir)  = e2des(ik,ir)
                  end do
               end do
c
            endif
c
            return
c
         elseif (switch(swioni).eq.11) then
c
c           Set-up for directly read Edge2D ionization
c
            if (cre2d.eq.0) then
               call redge2d(0)
c               cre2d = 1
            endif
c
            do ir = irlim1,irlim2
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) =  e2dnbs(ik,ir)
                  oldktibs(ik,ir) = e2dtibs(ik,ir)
                  oldktebs(ik,ir) = e2dtebs(ik,ir)
               end do
            end do
c
         elseif (switch(swioni).eq.15) then
c
c           Set-up for directly read Edge2D ionization
c
            if (cre2d.eq.0) then
               call redge2d(0)
c               cre2d = 1
            endif
c
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) =  e2dnbs(ik,ir)
                  oldktibs(ik,ir) = e2dtibs(ik,ir)
                  oldktebs(ik,ir) = e2dtebs(ik,ir)
                  knbs(ik,ir) = e2dnbs(ik,ir)
                  ktebs(ik,ir)= e2dtebs(ik,ir)
                  ktibs(ik,ir)= e2dtibs(ik,ir)
                  kvhs(ik,ir) = e2dvhs(ik,ir)
                  kes(ik,ir)  = e2des(ik,ir)
               end do
c
               if (ir.lt.irsep) then
c
                  oldknbs(nks(ir),ir) = e2dnbs(1,ir)
                  oldktibs(nks(ir),ir) = e2dtibs(1,ir)
                  oldktebs(nks(ir),ir) = e2dtebs(1,ir)
                  knbs(nks(ir),ir) = e2dnbs(1,ir)
                  ktebs(nks(ir),ir)= e2dtebs(1,ir)
                  ktibs(nks(ir),ir)= e2dtibs(1,ir)
                  kvhs(nks(ir),ir) = e2dvhs(1,ir)
                  kes(nks(ir),ir)  = e2des(1,ir)
c
               endif
c
            end do
c
c        Load previously calculated DIVIMP background for the first
c        iteration of SOL 22.
c
         elseif (switch(swioni).eq.16.or.switch(swioni).eq.17.or.
     >           switch(swioni).eq.18) then
c
c           Set-up for directly read DIV background plasma
c
            call readdivbg
c
c           Increment SOL 22 iteration count
c
            sol22_iter = sol22_iter + 1
c
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) =  knbs(ik,ir)
                  oldktibs(ik,ir) = ktibs(ik,ir)
                  oldktebs(ik,ir) = ktebs(ik,ir)
               end do
            end do
c
            do id = 1,nds
c
               oldkvds(id) = kvds(id)
c
            end do
c
c           Exit
c
            return
c
         endif
c
c
c
      elseif ((switch(swioni).eq.13.or.
     >         switch(swioni).eq.14)
     >         .and.sol22_iter.eq.1) then
c
         sol22_iter = sol22_iter + 1
c
         if (switch(swioni).eq.14) then
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) = e2dnbs(ik,ir)
                  oldktibs(ik,ir) = e2dtibs(ik,ir)
                  oldktebs(ik,ir) = e2dtebs(ik,ir)
                  knbs(ik,ir) = e2dnbs(ik,ir)
                  ktebs(ik,ir)= e2dtebs(ik,ir)
                  ktibs(ik,ir)= e2dtibs(ik,ir)
                  kvhs(ik,ir) = e2dvhs(ik,ir)
                  kes(ik,ir)  = e2des(ik,ir)
               end do
c
               if (ir.lt.irsep) then
                  knbs(nks(ir),ir) = e2dnbs(1,ir)
                  ktebs(nks(ir),ir)= e2dtebs(1,ir)
                  ktibs(nks(ir),ir)= e2dtibs(1,ir)
                  kvhs(nks(ir),ir) = e2dvhs(1,ir)
                  kes(nks(ir),ir)  = e2des(1,ir)
               endif
c
           end do
c
         else
c
            do ir = irlim1,irlim2
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) = e2dnbs(ik,ir)
                  oldktibs(ik,ir) = e2dtibs(ik,ir)
                  oldktebs(ik,ir) = e2dtebs(ik,ir)
                  knbs(ik,ir) = e2dnbs(ik,ir)
                  ktebs(ik,ir)= e2dtebs(ik,ir)
                  ktibs(ik,ir)= e2dtibs(ik,ir)
                  kvhs(ik,ir) = e2dvhs(ik,ir)
                  kes(ik,ir)  = e2des(ik,ir)
               end do
            end do
c
         endif
c
         return
c
c     ENDIF for .not.pinavail
c
      endif
c
c     For initial ionization option 17 - re-load the core from the
c     original plasma solution on every iteration. If the Tgrad option
c     is set to 0 then also load the previous private plasma
c     solution.
c
      if (switch(swioni).eq.16.or.switch(swioni).eq.17) then
c
c        Copy over target velocities so that mach 1 is not forced.
c
         do ir = irsep,nrs
            kvds(idds(ir,1)) = oldkvds(idds(ir,1))
            kvds(idds(ir,2)) = oldkvds(idds(ir,2))
         end do
c
      endif
c
      if (switch(swioni).eq.17) then
c
c        Reload core from original plasma solution - uses e2d variables
c        but will be loading a previous DIVIMP solution.
c
         do ir = 1,irsep-1
            do ik = 1,nks(ir)
               knbs(ik,ir) = e2dnbs(ik,ir)
               ktebs(ik,ir)= e2dtebs(ik,ir)
               ktibs(ik,ir)= e2dtibs(ik,ir)
               kvhs(ik,ir) = e2dvhs(ik,ir)
               kes(ik,ir)  = e2des(ik,ir)
            end do
         end do
c
c        If private plasma solution is also off - then reload original
c        for private plasma as well.
c
         if (ciopto.eq.0) then
c
            do ir = irtrap,nrs
               do ik = 1,nks(ir)
                  knbs(ik,ir) = e2dnbs(ik,ir)
                  ktebs(ik,ir)= e2dtebs(ik,ir)
                  ktibs(ik,ir)= e2dtibs(ik,ir)
                  kvhs(ik,ir) = e2dvhs(ik,ir)
                  kes(ik,ir)  = e2des(ik,ir)
               end do
            end do
c
         endif
c
      endif
c
c
c
c     FOR Radiation option 4 - read in the radiation data from 
c     a previous DIVIMP run. 
c
      call rzero(div_tpowls,maxnks*maxnrs) 
      call rzero(div_tcooliz,maxnks*maxnrs) 
      call rzero(div_cool,maxnks*maxnrs) 
c
      if (switch(swprad).eq.4) then 
c         
         call readdivaux('TPOWLS:',div_tpowls,
     >                             maxnks,maxnrs,1,1)
         call readdivaux('TCOOLIZ:',div_tcooliz,
     >                             maxnks,maxnrs,1,1)
c
c        Sum components of impurity electron cooling terms together  
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               div_cool(ik,ir) = div_tpowls(ik,ir)+
     >                           div_tcooliz(ik,ir) 
            end do
         end do
c
      endif  
c
c     For radiation option 5 calculate the distribution of radiation over the grid
c     given total radiation from each region. 
c     1) Inner divertor
c     2) Outer divertor
c     3) Inner and outer PFZ
c     4) Inner SOL to top
c     5) Outer SOL to top
c
c     Distribute Pin_region proportional to targ_flux(ir) * kvols(ik,ir) and 
c     integrate over region to get the appropriate total prad
c     Note: Pinqe and Pinqi options should be off when this option is used. 
c

      call pr_trace('CALCSOL_INTERFACE','BEFORE INIT')

c
c
c     Increase iteration count
c
      sol22_iter = sol22_iter + 1
c
c     Initialization
c
      call izero(cerr,maxnrs*2)
      call izero(cdeferr,maxnrs*2)
      call rzero(cserr,maxnrs*2)
      call rzero(cdefserr,maxnrs*2)
      call rzero(cdeferropt,maxnrs*2)
      call rzero(gperpa,maxnks*maxnrs)
c
c     Zero the correction factor array.
c
      call qzero(rconst,mxspts)
c
c     Zero the pressure array
c
c      call rzero(kpress,maxnks*maxnrs*2)
c
      call iinit(applydef,maxnrs*2,-1)
c
      if (ndef.gt.0) then
         do in = 1,ndef
            applydef(INT(deflist(in,1)),INT(deflist(in,2))) =
     .        deflist(in,3)
         end do
      end if
c
c     Read in EDGE 2D target data if Edge 2D compatibility is ON
c
      if (switch(swe2d).ne.0.0)then
         call redge2d(1)
      endif
c
c
c     Set up fixed input quantities
c
      title = 'DIVIMP'
      mb = crmb
      zb = rizb
      k0e = ck0
      k0i = ck0i
c
c     The following graph option will generate a series of plots based
c     on the SOL Option 22 data.They are currently turned off and in the
c     original code were read in as options.These are left for debugging
c     purposes and can be turned on by setting the options equal to 1.
c     This does mean however that the GHOST libraries would need to be
c     bound to DIVIMP for these plots. The grahing code can be easily
c     exised by commenting out the graph specific calls in each of the
c     subroutines. This will preserve the ability to print out the
c     relevant tables of information without needing to bind the GHOST
c     libraries.
c
c     David Elder                1995, May 16
c
      graph    = 1
      graphaux = 1
      graphvel = 1
      graphran = 0.0
c
c     Adjust Ionization switch settings if appropriate
c
c     Check to see if option set to use PIN data to
c     normalize the analytic sources.
c
      pinnorm = 0
c
      if (switch(swion).eq.8.0.and.pinavail) then
         switch(swion) = switch(swioni)
         switch(swioni) = 8.0
         pinnorm = 1
      endif
c
c     Calculate the second gradient taken from the OLD or
c     previous solution - uses values stored in KNBS -
c
      call rzero(grad_oldknbs,maxnks*maxnrs)
c
      if (switch(swgperp).eq.7.or.switch(swgperpp).eq.7.or.
     >    switch(swgperp).eq.8.or.switch(swgperpp).eq.8
     >    ) then
c
         do ir = irlim1,irlim2
            do ik = 1,nks(ir)
               call quantgrad(ir,kss(ik,ir),oldknbs,
     >                        knds,gradi,grado,2)
c
               len1 = polysidelen(ik,ir,INWARD41,rc)
               len2 = polysidelen(ik,ir,OUTWARD23,rc)
               flux = gradi*len1 -  grado*len2
c
               if (kareas(ik,ir).ne.0.0) then
                  grad_oldknbs(ik,ir)= flux/kareas(ik,ir)
               else
                  grad_oldknbs(ik,ir)= 1.0
               endif
c
c
c               write (6,'(a,2i4,15(1x,g12.5))') 'GP:',ik,ir,
c     >                 cdperp*grad_oldknbs(ik,ir),
c     >                 gradi,grado,flux,karea2(ik,ir),
c     >                 len1,len2,len1/len2
c
            end do
         end  do
c
      endif
c
c
c      write (6,*) 'Swion:',switch(swion),switch(swioni)
c
C
c     Loop through rings
c
c     All quantities that are used by the calcsol routine must be
c     copied into appropriate variables. There are two reasons for
c     this - first to allow the calcsol routines to use a different
c     level of precision than the rest of the code. Second to make the
c     contents of calcsol and its subroutines as independent of DIVIMP
c     as possible - so that changes and additions may be easily made.
c
      if (ctestsol.gt.0.0) then
         irlim1 = int(ctestsol)
         irlim2 = int(ctestsol)
      endif
c
c     Set debugging ring
c
      irdebug = irlim1
c
c      write (6,*) 'IR debug:',irdebug
c
c     Set the starting cell for all rings - if the E2D compatibility option
c     is set to 9 - set the starting cell not equal to 1 - may be an
c     input parameter in later revisions.
c
      if (switch(swe2d).eq.9) then
c
         ike2d_start = ike2d
c
      else
c
         ike2d_start = 1
c
      endif

c
c     Calculate the target particle and power outfluxes.
c
c     Temporarily set pfz_dist_opt and pfz_dist_param until these
c     are added to the input file
c
      pfz_dist_opt = 1
      pfz_dist_param(1) = sepdist2(idds(irsep+4,1))
      pfz_dist_param(2) = sepdist2(idds(irsep+4,2))
c
      call pr_trace('CALCSOL_INTERFACE','BEFORE CALCFLUXES')

c
      call calcfluxes(gtarg,ionptarg,elecptarg,e2dgtarg,
     >                presstarg,gamcor,gamecor,ike2d_start,
     >                g_pfzsol,pe_pfzsol,pi_pfzsol,pr_pfzsol,
     >                pfz_dist_opt,pfz_dist_param)
c
c     If pin is available
c     Call routine to calculate GPERP CORection factors
c
      call pr_trace('CALCSOL_INTERFACE','BEFORE CALCSOLIZ')

      if (pinavail) then
c
c        Call routine to print comparison of sources.
c
         if (cprint.eq.9) call ioniz_comp
c
         call calcsoliz(rconst,recfrac,gtarg,areasum,
     >                  gperpa,oldknbs,grad_oldknbs,pinion,pinrec,
     >                  gperprat,ike2d_start)
c
      elseif ((switch(swioni).eq.11.or.switch(swioni).eq.15)
     >         .and.(.not.pinavail)) then
         if (switch(swrecom).eq.2.0) then
            call calcsoliz(rconst,recfrac,gtarg,areasum,
     >                  gperpa,oldknbs,grad_oldknbs,e2dion,e2dhrec,
     >                  gperprat,ike2d_start)
         else
            call calcsoliz(rconst,recfrac,gtarg,areasum,
     >                  gperpa,oldknbs,grad_oldknbs,e2dion,pinrec,
     >                  gperprat,ike2d_start)
         endif
      endif
c
c     Use the EDGE2D background through ALL iterations
c     - essentially bypasses SOL 22 COMPLETELY -
c     Used as a method of iterative PIN testing while still
c     obtaining the SOL 22 flux debugging information from calcsoliz.
c
      call pr_trace('CALCSOL_INTERFACE','BEFORE EDGE2D SAVE')

      if (switch(swe2d).eq.4.or.switch(swe2d).eq.5.or.
     >    switch(swe2d).eq.6.or.switch(swe2d).eq.7) then
            do ir = irlim1,irlim2
               do ik = 1,nks(ir)
                  oldknbs(ik,ir) = e2dnbs(ik,ir)
                  oldktibs(ik,ir) = e2dtibs(ik,ir)
                  oldktebs(ik,ir) = e2dtebs(ik,ir)
                  knbs(ik,ir) = e2dnbs(ik,ir)
                  ktebs(ik,ir)= e2dtebs(ik,ir)
                  ktibs(ik,ir)= e2dtibs(ik,ir)
                  kvhs(ik,ir) = e2dvhs(ik,ir)
                  kes(ik,ir)  = e2des(ik,ir)
               end do
            end do
            return
      endif
c

      call pr_trace('CALCSOL_INTERFACE','BEFORE RING LOOP')

      do ir = irlim1, irlim2
c
         call set_ikvals(ir,ikstart,ikend,ikopt)
c
c
c        Bypass ALL of the code for ir = irwall and ir = irtrap
c        since these are rings with only virtual boundary cells.
c        They should not be used in SOL option 22 - the solutions
c        for rings irwall-1 and irtrap + 1 will eb copied in to
c        ensure valid values for ne,te,ti just in case they are
c        required.
c
         if (ir.eq.irwall.or.ir.eq.irtrap) goto 1000
c
c        Set the values of lengths to a maximum of 1/2 of the field
c        line length. (Same for lenr)
c
         ringlen = ksmaxs2(ir)
         ringnum = ir
c
c        Set private plasma ring indicator
c
         if (ringnum.gt.irwall) then
            pplasma = 1
         else
            pplasma = 0
         endif
c
c        If the private plasma ionization option specifies the
c        ANALYTIC option - then call this for Private Plasma
c        rings.
c
         if (switch(swionp).eq.-2.0.and.pplasma.eq.1) then
c
            call specplas(ir,ir,ikopt)
c
            do ik = 1,nks(ir)
               oldknbs(ik,ir) = knbs(ik,ir)
               oldktibs(ik,ir) = ktibs(ik,ir)
               oldktebs(ik,ir) = ktebs(ik,ir)
            end do
c
            goto 1000
c
c            cycle
c
c...  Assign PP plasma using Thomson data:
c
         elseif ((switch(swionp).eq.-3.0
     >        .or.switch(swionp).eq.-4.0)
     >       .and.pplasma.eq.1) then
c
c slmod begin - new
c...LITTLE FIX:
            call thompp(ir,ir,ikopt,int(abs(switch(swionp))))
c
c            call thompp(ir,ir,ikopt,int(abs(swionp)))
c slmod end
c
            do ik = 1, nks(ir)
               oldknbs (ik,ir) = knbs (ik,ir)
               oldktibs(ik,ir) = ktibs(ik,ir)
               oldktebs(ik,ir) = ktebs(ik,ir)
            end do
c
            goto 1000
c slmod begin - new
         elseif ((switch(swionp).eq.-5.0
     >        .or.switch(swionp).eq.-6.0)
     >       .and.pplasma.eq.1) then
c...        Assign private plasma from a listing of T and n in the input
c           file:
            call assignpp(ir,ir,ikopt,int(abs(switch(swionp))))
            do ik = 1, nks(ir)
               oldknbs (ik,ir) = knbs (ik,ir)
               oldktibs(ik,ir) = ktibs(ik,ir)
               oldktebs(ik,ir) = ktebs(ik,ir)
            end do
            goto 1000
c slmod end
c
c            cycle
c
c
c
         end if
c
c        Set up the effective GPERP/core correction term.
c
         gnet = rconst(ir)
c
c
c         if ((pinavail
c     >       .and.((pinnorm.eq.1.or.actswion.eq.2.0))
c     >       .and.(actswgperp.eq.2.0.or.actswgperp.eq.6.0
c     >             .or.actswgperp.eq.7.or.actswgperp.eq.8.0))
c     >       .or.
c     >       ((actswioni.eq.11.or.actswioni.eq.15).and.
c     >        (actswgperp.eq.2.0.or.actswgperp.eq.6.0
c     >         .or.actswgperp.eq.7.or.actswgperp.eq.8.0)
c     >       .and.(.not.pinavail))
c     >       ) then
c            gnet = rconst(ir)
c
c
c         else
c            gnet = 0.0
c         endif
c
c
c
c        Set up the power distribution length for the
c        power input options 7 and 8
c
         spow = spowbeg * ringlen
         spow2= spowlen * ringlen
c
c         write (6,*) 'SPOW:',spow,spow2
c
c        Set up gperp lengths for particle source compensation
c
         sgperpbeg = gperpbegf * ringlen
         sgperpend = gperpendf * ringlen
c
         midnks = ikmids(ir)
c         write(0,*) 'solascv:midnks:',midnks,ir,sol22_halfringlen_opt
         
c
c-----------------------------------------------------------------
c        OUTER TARGET
c-----------------------------------------------------------------
c
c
c        Only execute for the selected half-rings - ikopt
c
         if (ikopt.eq.1.or.ikopt.eq.3) then
            call pr_trace('CALCSOL_INTERFACE','START OUTER TARGET')
c
c        Check for sol22 debug start
c
c           The parameters are current ring snd local IKOPT
            if (debug_sol22.ne.0) call check_init_record_data(ir,1)

c
c        Set the values for a call - do 1/2 of the ring at a time.
c
c
c        Outer target first - i.e. elements 1 to 1/2 nks(ir)
c                                  target designated as 2  
c
c
c        Outer target
c
c
c        Reset switches if necessary
c
         seterror = 0
c
c        Set active switches
c
         if (applydef(ir,2).eq.-1) then
            call setsw(-1,pplasma,new_errlevel)
c
c           Reset Gperp switch for option 3 in case of under ionization
c
            if (( (pinavail.or.((.not.pinavail
     >          .and.(actswioni.eq.11.or.actswioni.eq.15))))
     >          .and.actswgperp.eq.3.0
     >          .and.rconst(ir).gt.0.0)) then
                actswgperp = 2.0
                write(6,*) 'Swgperp:',rconst(ir),actswgperp
            elseif ((.not.pinavail)
     >             .and.(actswioni.ne.11.and.actswioni.ne.15)
     >             .and.(actswgperp.eq.2.0.or.
     >                   actswgperp.eq.3.0.or.actswgperp.eq.7.0.or.
     >                   actswgperp.eq.8.0)) then
                actswgperp = 0.0
            elseif (actswgperp.eq.7.0.and.gperprat(ir).gt.5.0) then
                actswgperp = 2.0
            endif
c
         elseif (applydef(ir,2).ge.0) then
            call setsw(applydef(ir,2),pplasma,new_errlevel)
            actswerror = applydef(ir,2)
c
            if (actswerror.eq.0.0) seterror = 2
c
            cdeferr(ir,2) =  8
            cdefserr(ir,2) = 0.0
            cdeferropt(ir,2) = applydef(ir,2)
         endif
c
c        Initialize HALF ring length to be used
c
c        halfringlen = ksb(midnks,ir)
c
         if (sol22_halfringlen_opt.eq.0) then 
            halfringlen  = 0.5d0*ringlen
         elseif (sol22_halfringlen_opt.eq.1) then 
            ! ring midpoint is upper boundary of middle cell just below midpoint
            halfringlen = ksb(midnks+1,ir)
c            write(6,'(a,2i8,10(1x,g12.5))') 'Halfringlen1:',
c     >           midnks,ir,halfringlen,ringlen/2.0
c            write(0,'(a,2i8,10(1x,g12.5))') 'Halfringlen1:',
c     >           midnks,ir,halfringlen,ringlen/2.0
         endif
c
c        Initialize the ionization and radiation source lengths
c
         call initlen
c
c        Initialize the friction parameter for momentum loss options
c
         actffric = find_ffric(ir,2,actlenmom)
c
c        Initialize major radius and target condition data 
c
         targfact = 1.0
c
         if (actswmajr.eq.1.0) then
            if (actswe2d.eq.0.0) then
               id = idds(ir,2)
               targfact = rp(id) / r0
            elseif (actswe2d.ne.0.0) then
               targfact = rs(1,ir) / r0
            endif
         endif
c
         te0 = kteds(idds(ir,2))
         ti0 = ktids(idds(ir,2))
         n0 = knds(idds(ir,2))
         v0 = kvds(idds(ir,2))
c
         if (forcet.eq.1) then
            ttarg = (((5.0* te0 + 3.5 * ti0)*sqrt(te0+ti0))
     >                      / (8.5*sqrt(2.0)))**(2.0/3.0)
            te0 = ttarg
            ti0 = ttarg
            kteds(idds(ir,2)) = te0
            ktids(idds(ir,2)) = ti0
         endif
c
         if (actswe2d.eq.9.0) then
            n1 = e2dnbs(ike2d_start,ir)
            te1= e2dtebs(ike2d_start,ir)
            ti1= e2dtibs(ike2d_start,ir)
            v1e2d = e2dvhs(ike2d_start,ir)
            vpe2d = v1e2d
c
            vpg = (e2dgpara(ike2d_start,ir)+e2dgpara(ike2d_start+1,ir))
     >                 / (2.0 * n1)
c
            if (forcet.eq.1) then
               ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))
     >                      / (8.5 * sqrt(2.0)))**(2.0/3.0)
               te1 = ttarg
               ti1 = ttarg
            endif
c

c
            write(6,'(a,i4,g16.8,2f9.4,3g14.6)')
     >                 'E2d-9:SOL 22:First:',ir,n1,te1,ti1,v1e2d
     >                  ,vpe2d,vpg
c
         elseif (actswe2d.ne.0.0) then
            n1 = cellvals(ir,1,2)
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
            v1e2d = cellvals(ir,4,2)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*mb))
            vpg = gtarg(ir,2) / n1
c
            if (forcet.eq.1) then
               ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))
     >                      / (8.5 * sqrt(2.0)))**(2.0/3.0)
               te1 = ttarg
               ti1 = ttarg
            endif
c
            write(6,'(a,i4,g16.8,2f9.4,3g14.6)')
     >                 'E2d:SOL 22:First:',ir,n1,te1,ti1,v1e2d
     >                  ,vpe2d,vpg
c
         endif
c
         errcode = 0
         serr = 0.0
c
         npts = midnks
c
c        This loop would also be used to load other
c        quantities - like ionization source or background
c        neutral densities if these were to be used in the
c        calcsol subroutine.
c
         rbnd(0) = krb(0,ir)
         sbnd(0) = ksb(0,ir)
c
         do ik = 1,midnks
c
            spts(ik) = kss2(ik,ir)
c
c            write(6,*) 'in1:',ir,ik,spts(ik)
c
            fact = recfrac
c
            if (actswmajr.eq.2.0) then
               fact = rs(ik,ir)/r0 * fact
            elseif (actswmajr.eq.3.0) then
               fact = r0/rs(ik,ir) * fact
            endif
c
            if ((actswioni.eq.11.or.actswioni.eq.15)
     >           .and.(.not.pinavail)) then
               ionsrc(ik) = e2dion(ik,ir) * fact
            else
               ionsrc(ik) = pinion(ik,ir) * fact
            endif
c
            if (actswrecom.eq.2.and.(.not.pinavail)) then
               recsrc(ik) = e2dhrec(ik,ir)
            else
               recsrc(ik) = pinrec(ik,ir)
            endif
c
            gperp(ik)  = gperpa(ik,ir)
c
            if (oldktebs(ik,ir).lt.tcutqe) then
               qesrc(ik)  = 0.0
               write(6,'(a5,i4,2(1x,g13.5))') 'QE:',ik,
     >                         oldktebs(ik,ir),tcutqe
            else
               qesrc(ik)  = -pinqe(ik,ir) * fact
     >                                    * qesrc_mult
            endif
c
            if (oldktibs(ik,ir).lt.tcutcx) then
               qisrc(ik)  = 0.0
               write(6,'(a5,i4,2(1x,g13.5))') 'QI:',ik,
     >                         oldktibs(ik,ir),tcutcx
            else
               qisrc(ik)  = -pinqi(ik,ir) * fact
            endif
c
            radsrc(ik) = div_cool(ik,ir)
c
            nhs(ik) = PINATOM(IK,IR)
            nh2s(ik) = PINMOL(ik,ir)
            ths(ik) = pinena(ik,ir) * 0.6666
            nhs0(ik) = pinvdist(1,1,ik,ir) + pinvdist(2,1,ik,ir)
     >                + pinvdist(3,1,ik,ir)
c
c
c           Load values from last iteration for PINQID calculations
c
            oldne(ik) = oldknbs(ik,ir)
            oldte(ik) = oldktebs(ik,ir)
            oldti(ik) = oldktibs(ik,ir)
c
            rbnd(ik) = krb(ik,ir)
            sbnd(ik) = ksb(ik,ir)
c
            momsrc(ik) = pinmp(ik,ir) * fact
c
         end do
c
c        The end of the Spts array can either be the cell boundary
c        or the mid-point of the ring. The cell boundary will result
c        in more correct integrations since the contribution of the 
c        cells near the mid-point will be calculated more accurately.
c        (Though the difference would be expected to be small.)
c        However, the ring midpoint (0.5 * ringlen) makes it easier to 
c        scale or extend sources over the same spatial region by allowing
c        input to be specified proportional to the ring length. 
c
c         spts(npts+1) = sbnd(npts)
c
         spts(npts+1) = halfringlen
c
c        Find an estimated S-value of the X-point
c
         sxp = 0.0
c
         do ik = 1,midnks
            if (zs(1,ir).lt.zxp) then
               if (zs(ik,ir).ge.zxp) then
                  sxp = (kss2(ik,ir) - kss2(ik-1,ir)) *
     >                  (zxp-zs(ik-1,ir))
     >                  /(zs(ik,ir)-zs(ik-1,ir))
     >                  + kss2(ik-1,ir)
                  goto 100
               endif
            elseif (zs(1,ir).gt.zxp) then
               if (zs(ik,ir).lt.zxp) then
                  sxp = (kss2(ik,ir) - kss2(ik-1,ir)) *
     >                  (zxp-zs(ik-1,ir))
     >                  /(zs(ik,ir)-zs(ik-1,ir))
     >                  + kss2(ik-1,ir)
                  goto 100
               endif
            endif
         end do
100      continue
c
         write (6,*) 'Sxp:',ir,nks(ir),kss2(midnks,ir),
     >               sxp,ringlen
c
c        Set both target power values -
c
c        Note:  pae_end and pai_end have
c        a "-" ve because a "-" sign is used in the embedded power
c        equations.
c
         pae_start = -elecptarg(ir,2)
         pai_start = -ionptarg(ir,2)
c
         pae_end = -elecptarg(ir,1)
         pai_end = -ionptarg(ir,1)
c
c        Private Plasma target power loss re-distribution to main SOL
c
c        Turned off in private plasma. Electron then ion.
c
         if (actswppelec.eq.0.0.or.pplasma.eq.1) then

            ppelecpow = 0.0

         elseif (actswppelec.eq.1.0.or.actswppelec.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppelec = 0.0
               ppelecpow   = 0.0
c
            else
c
               ppelecpow = -elecptarg(ircor,2)
c
            endif

         elseif (actswppelec.eq.3.0) then
c
c           see comment above regarding "-" sign
c
            ppelecpow = -pe_pfzsol(ir,2)
c
         endif
c
c        PP Ion power loss term
c
         if (actswppion.eq.0.0.or.pplasma.eq.1) then

            ppionpow = 0.0

         elseif (actswppion.eq.1.0.or.actswppion.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppion = 0.0
               ppionpow   = 0.0
c
            else
c
               ppionpow = -ionptarg(ircor,2)
c
            endif
c
         elseif (actswppelec.eq.3.0) then
c
c           see comment above regarding "-" sign
c
            ppionpow = -pi_pfzsol(ir,2)
c
         endif
c
c        PP pressure loss term
c
         if (actswppress.eq.0.0.or.pplasma.eq.1) then

            pp_press = 0.0

         elseif (actswppress.eq.1.0.or.actswppress.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppress = 0.0
               pp_press   = 0.0
c
            else
c
               pp_press = presstarg(ircor,2)
c
            endif
c
            write(0,*) 'Press2:',actswppress,ir,ircor,pp_press
c
         elseif (actswppress.eq.3.0) then
c
c           see comment above regarding "-" sign
c
            pp_press = pr_pfzsol(ir,2)
c
         endif

c
c        Set up the ionization source - if choice option is selected.
c
         if (actswion.eq.5.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.5.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1369
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.5) then
                   actswion = 3.0
                else
                   actswioni = 3.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 13.0 - 10.0 * te0
                   ssrcfi = ssrcst + 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = ssrcst + 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.5) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion5O:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
c
         elseif ((actswion.eq.6.0.or.actswion.eq.9.0)
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.(actswioni.eq.6.0.or.actswioni.eq.9)
     >          .and.(.not.pinavail))) then
c
c           Adjust parameters for S**5 Gaussian ionization option.
c
c           This ionization source always starts at zero.
c           Thus the quantity in the input for the start position
c           is used instead as the width factor for the distribution.
c           If the Start point is listed as 0.0 then this is adjusted
c           to 1.0 for the width factor for the gaussian.
c
c           This code only needs to be executed once on each ring.
c
            if (ssrcdecay.eq.0.0) then
               s5gausslen = 1.0
            else
               s5gausslen = ssrcdecay
            endif
c
            ssrcst = 0.0
c
            ssrclen = ssrcfi - ssrcst
c
         elseif (actswion.eq.7.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.7.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1372
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.7) then
                   actswion = 6.0
                else
                   actswioni= 6.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 14.0 - 10.0 * te0
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.0
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.7) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion7O:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
c
         elseif (actswion.eq.10.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.10.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1372
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.10.0) then
                   actswion = 9.0
                else
                   actswioni= 9.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5*(14.0 - 10.0 * te0)
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.10.0) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion10O:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
c
         endif
c
c        Set up extra perpenicular source parameters
c
         if (switch(swextra).gt.0.0) then
c
c           Calculate source start and stop positions - for
c           first 1/2 ring use the parameters as specified.
c
c           Note: These stop and start positions may be
c           on either half of the flux tube or even overlap
c           the top. The integrations procedd outward from
c           each target so the values of the start and
c           stop positions will have to be adjusted but
c           the rest will remain the same.
c
            start_gextra_src = gextra_src_start * ringlen
            stop_gextra_src  = gextra_src_stop  * ringlen
            start_gextra_sink= gextra_sink_start * ringlen
            stop_gextra_sink = gextra_sink_stop * ringlen
c
c           Calculate extra source strength
c

c
            if (start_gextra_src.ne.stop_gextra_src) then
               gextra_src = abs(gtarg(ir,3)) * gextra_mult /
     >                   abs(start_gextra_src-stop_gextra_src)
            else
               gextra_src = abs(gtarg(ir,3)) * gextra_mult
            endif
c
c           Calculate extra sink strength
c
            if (start_gextra_sink.ne.stop_gextra_sink) then
               gextra_sink = -abs(gtarg(ir,3)) * gextra_mult /
     >                   abs(start_gextra_sink-stop_gextra_sink)
            else
               gextra_sink = -abs(gtarg(ir,3)) * gextra_mult
            endif
c
            write (6,*) 'Gextra 1:',ir,gextra_src,gextra_sink,
     >              gtarg(ir,3),ringlen,
     >              start_gextra_src,stop_gextra_src,
     >              start_gextra_sink,stop_gextra_sink
c

         endif
c
c        Error correcting branch point
c
 150     continue
c
c
c        Initialize output arrays to default error values
c
         call dinit(te,mxspts,10.0d0)
         call dinit(ti,mxspts,10.0d0)
         call dinit(ne,mxspts,1.0d19)
         call dinit(vb,mxspts,0.0d0)
c

         write (6,'(a,i4,a)') '************   START FIRST HALF ',
     >               ir,'  ***********'
         write (6,*) 'First:',ir,te0,ti0,n0,midnks
c
c        Zero power ratio information in case calcsol is not called.
c
         call qzero (int_powrat,3)
c
c        Override with detached plasma solution for the outer
c
         if (switch(swdetach).eq.1.0) then
c
            call detached_plasma(spts,npts,errcode,serr,
     >                           te,ti,ne,vb,ir,switch(swdetach),
     >                           te0,ti0,n0,v0,act_press,mb)
c
c        Solve using SOL22
c
         else
c
c slmod begin - new
c...Need a system routine for the other system files.  Rename 
c   CLOCK2:
            IF (osm_mode.GE.2) THEN
              deltat = Clock2()
              CALL OpenStorageFiles(ir,IKLO,'tmp.dat')
            ENDIF
c slmod end

            call calcsol (spts,npts,errcode,serr,
     >                 te,ti,ne,vb,exp_press,act_press,
     >                 prad,ir,irlim1,
     >                 int_powrat,cprint)
c slmod begin - new
            IF (osm_mode.GE.2) THEN
              CALL SOL22Status(IKLO,ir,deltat,serr,spts,npts,errcode)
            ENDIF
c slmod end
c

         endif
c
c        Handle Error Condition if Error Switch is set
c


         if ((   actswerror.ge.1.0
     >            .and.(errcode.eq.3.or.errcode.eq.4.or.
     >                  errcode.eq.5.or.errcode.eq.6.or.
     >                  errcode.eq.7))
     >       .or.(actswerror.eq.0.0.and.
     >           (errcode.eq.6.or.errcode.eq.7))) then
c
c           If a negative N error or NaNQ occurs without error correction -
c           Turn error correction ON.
c
            if (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7)
     >          .and.seterror.eq.0) then
               actswerror = 10
               seterror = 1
            elseif (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7)
     >              .and.seterror.eq.1) then

                seterror = 2
            elseif (seterror.eq.2) then
c
c               ERROR - negative N has been found EVEN with highest level
c                       of error correction - issue error messages and stop.
c
                call errmsg('SOLASCV:SOL22',
     >              ' Unsolvable Negative N error encountered.'//
     >              ' Program Stopping')
c
c                write (7,*) 'SOLASCV: SOL22:'//
c     >               ' Unsolvable Negative N error encountered'
c                write (7,*) 'Program will STOP'
c
                stop 'SOL22:NEG N'
c
            endif
c
c           Record error
c
            cdeferr(ir,2) = errcode
            cdefserr(ir,2) = serr
c
c           Set switches
c
            call setsw(-2,pplasma,new_errlevel)
c
c           Record error level of next attempt
c
c            cdeferropt(ir,2) = actswerror
c
            cdeferropt(ir,2) = new_errlevel
c
            call initlen
c
            write (6,*) 'Error Handler: SET OUTER:',cdeferr(ir,2),
     >                   cdefserr(ir,2),cdeferropt(ir,2)
c
c           Re-do Outer 1/2 ring
c
            goto 150
c
         endif
c
c        Save background density and temperature
c        Target first in case it has changed
c
         knds(idds(ir,2)) = n0
         kvds(idds(ir,2)) = v0
c
c        If E2D option 9 was in use - fill in the missing
c        plasma with a linear fit.
c
         if (actswe2d.eq.9.0.and.switch(swdetach).ne.1.0) then
c
c           Fill in missing values - depending on option specified.
c
            if (fillopt.eq.0) then
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0 + (te(ike2d_start)-te0) *
     >                      spts(ik)/spts(ike2d_start)
                  ti(ik) = ti0 + (ti(ike2d_start)-ti0) *
     >                      spts(ik)/spts(ike2d_start)
                  ne(ik) = n0 + (ne(ike2d_start)-n0) *
     >                      spts(ik)/spts(ike2d_start)
                  vb(ik) = v0 + (vb(ike2d_start)-v0) *
     >                      spts(ik)/spts(ike2d_start)
               end do
c
            elseif (fillopt.eq.1) then
c
               gam = n0 * v0
c
c              Extrapolate target electron temperature
c
               slope = (te(ike2d_start+1) - te(ike2d_start)) /
     >                 (spts(ike2d_start+1) - spts(ike2d_start))
c
               te0 = - spts(ike2d_start) * slope + te(ike2d_start)
c
               if (te0.le.0.0) te0 = te(ike2d_start)
c
c              Extrapolate target ion temperature
c
               slope = (ti(ike2d_start+1) - ti(ike2d_start)) /
     >                 (spts(ike2d_start+1) - spts(ike2d_start))
c
               ti0 = - spts(ike2d_start) * slope + ti(ike2d_start)
c
               if (ti0.le.0.0) ti0 = ti(ike2d_start)
c
               v0 = -sqrt((te0+ti0)/crmb * econv/mconv)
c
               n0 = gam/v0
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0 + (te(ike2d_start)-te0) *
     >                      spts(ik)/spts(ike2d_start)
                  ti(ik) = ti0 + (ti(ike2d_start)-ti0) *
     >                      spts(ik)/spts(ike2d_start)
                  ne(ik) = n0 + (ne(ike2d_start)-n0) *
     >                      spts(ik)/spts(ike2d_start)
                  vb(ik) = v0 + (vb(ike2d_start)-v0) *
     >                      spts(ik)/spts(ike2d_start)
               end do
c
c              Reset target values
c
               knds(idds(ir,2)) = n0
               kvds(idds(ir,2)) = v0
               kteds(idds(ir,2)) = te0
               ktids(idds(ir,2)) = ti0
c
            elseif (fillopt.eq.2) then
c
c              Fill with constant values from last point
c
c              Target values are NOT affected
c
               do ik = 1,ike2d_start-1
                  te(ik) = te(ike2d_start)
                  ti(ik) = ti(ike2d_start)
                  ne(ik) = ne(ike2d_start)
                  vb(ik) = vb(ike2d_start)
               end do
c
            elseif (fillopt.eq.3) then
c
c              Fill with values held constant at target conditions
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0
                  ti(ik) = ti0
                  ne(ik) = n0
                  vb(ik) = v0
               end do
c
            endif
c
         endif
c
c        Assign background ...
c
         do ik = 1,midnks
c
            ktebs(ik,ir) = te(ik)
            ktibs(ik,ir) = ti(ik)
            knbs(ik,ir) = ne(ik)
            oldknbs(ik,ir) = ne(ik)
            oldktibs(ik,ir) = ti(ik)
            oldktebs(ik,ir) = te(ik)
            kvhs(ik,ir) = vb(ik)
c slmod begin - new
c            oldkvhs(ik,ir) = vb(ik)

            osmcde(ik,ir) = cde(ik)
            osmcdi(ik,ir) = cdi(ik)
            osmcve(ik,ir) = cve(ik)
            osmcvi(ik,ir) = cvi(ik)
c slmod end
c
            kpress(ik,ir,1) = exp_press(ik)
            kpress(ik,ir,2) = act_press(ik)
            kprad(ik,ir)    = prad(ik)
c
         end do
c slmod begin - new
c... .TRUE. temporary:
         IF (.NOT.pinavail.AND.(rel_opt.EQ.1.OR.rel_opt.EQ.3))
     .     CALL CalcInitSrc(IKLO,ir)
c slmod end
c
c        Save mach numbers and error codes
c
         cmachno(ir,2) = m0
         cerr(ir,2)  = errcode
         cserr(ir,2) = serr
c
c        Save Power Ratios
c
         sol22_power_ratio(ir,2,1) = int_powrat(1)
         sol22_power_ratio(ir,2,2) = int_powrat(2)
         sol22_power_ratio(ir,2,3) = int_powrat(3)
c
         write(6,'(a,i4,3(1x,g12.5))') 'PR2:',ir,int_powrat(1),
     >                    int_powrat(2),int_powrat(3)
c
c        End of 1/2 ring selection
c

c           Check to print debugging data
c           The parameters are current ring snd local IKOPT
c           data will be printed if it has been collected
            if (debug_sol22.ne.0) call check_print_data(ir,1)

         endif
c
c
c-----------------------------------------------------------------
c        INNER TARGET
c-----------------------------------------------------------------
c
c
c        Only execute for the selected half-rings - ikopt
c

         if (ikopt.eq.2.or.ikopt.eq.3) then
            call pr_trace('CALCSOL_INTERFACE','START INNER TARGET')

c
c
c        Check for sol22 debug start
c
c           The parameters are current ring snd local IKOPT
            if (debug_sol22.ne.0) call check_init_record_data(ir,2)
c
c        Set the values for a call - do 1/2 of the ring at a time.
c
c        Do second half of the ring.
c
c        Inner target second - i.e. elements 1/2 nks(ir) to nks(ir)
c                                   target designated as 1  
c
c        Inner target
c
c        Reset switches
c
c
         seterror = 0
c
         if (applydef(ir,1).eq.-1) then
            call setsw(-1,pplasma,new_errlevel)
c
c           Reset Gperp switch for option 3 in case of under ionization
c
            if (( (pinavail.or.((.not.pinavail
     >          .and.(actswioni.eq.11.or.actswioni.eq.15))))
     >          .and.actswgperp.eq.3.0
     >          .and.rconst(ir).gt.0.0)) then
                actswgperp = 2.0
                write(6,*) 'Swgperp:',rconst(ir),actswgperp
            elseif ((.not.pinavail)
     >             .and.(actswioni.ne.11.and.actswioni.ne.15)
     >             .and. (actswgperp.eq.2.0.or.
     >                    actswgperp.eq.3.0.or.actswgperp.eq.7.or.
     >                    actswgperp.eq.8.0)) then
                actswgperp = 0.0
            elseif (actswgperp.eq.7.0.and.gperprat(ir).gt.5.0) then
                actswgperp = 2.0
            endif
c
         elseif (applydef(ir,1).ge.0) then
            call setsw(applydef(ir,1),pplasma,new_errlevel)
            actswerror = applydef(ir,1)
c
            if (actswerror.eq.0.0) seterror = 2
c
            cdeferr(ir,1) =  8
            cdefserr(ir,1) = 0.0
            cdeferropt(ir,1) = applydef(ir,1)
         endif
c
c        Initialize HALF ring length to be used
c        - keep in mind that the second half of the ring works down
c          from KSMAXS2(IR).    
c     
c        halfringlen = ksmaxs2(ir) - ksb(midnks,ir)
c
         if (sol22_halfringlen_opt.eq.0) then 
            halfringlen  = 0.5d0*ringlen
         elseif (sol22_halfringlen_opt.eq.1) then 
            ! ring midpoint is upper boundary of middle cell just below midpoint
            halfringlen = ksmaxs2(ir) - ksb(midnks+1,ir)
c            write(6,'(a,2i8,10(1x,g12.5))') 'Halfringlen2:',
c     >           midnks,ir,halfringlen,ringlen/2.0
c            write(0,'(a,2i8,10(1x,g12.5))') 'Halfringlen2:',
c     >           midnks,ir,halfringlen,ringlen/2.0
         endif
c
c        Initialize the ionization and radiation source lengths
c
         call initlen
c
c        Initialize the friction parameter for momentum loss options
c
         actffric = find_ffric(ir,1,actlenmom)
c
c        Initialize major radius and target condition data 
c
         targfact = 1.0
c
         if (actswmajr.eq.1.0) then
            if (actswe2d.eq.0.0) then
               id = idds(ir,1)
               targfact = rp(id) / r0
            elseif (actswe2d.ne.0.0) then
               targfact = rs(nks(ir),ir)/r0
            endif
         endif
c
c
c        The solvers always deal with 1/2 a ring at a time and need a
c        target velocity < 0 for all targets. 
c
c
         te0 = kteds(idds(ir,1))
         ti0 = ktids(idds(ir,1))
         n0 = knds(idds(ir,1))
         v0 = -kvds(idds(ir,1))
c
         if (forcet.eq.1) then
            ttarg = (((5.0* te0 + 3.5 * ti0) * sqrt(te0+ti0))
     >                      / (8.5 * sqrt(2.0)))**(2.0/3.0)
            te0 = ttarg
            ti0 = ttarg
            kteds(idds(ir,1)) = te0
            ktids(idds(ir,1)) = ti0
         endif
c
         if (actswe2d.eq.9.0) then
            n1 = e2dnbs(nks(ir)-ike2d_start+1,ir)
            te1= e2dtebs(nks(ir)-ike2d_start+1,ir)
            ti1= e2dtibs(nks(ir)-ike2d_start+1,ir)
            v1e2d = -e2dvhs(nks(ir)-ike2d_start+1,ir)
            vpe2d = v1e2d
c
            vpg = -(e2dgpara(nks(ir)-ike2d_start+1,ir)
     >             +e2dgpara(nks(ir)-ike2d_start+1+1,ir))
     >                 / (2.0 * n1)
c
            if (forcet.eq.1) then
               ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))
     >                      / (8.5 * sqrt(2.0)))**(2.0/3.0)
               te1 = ttarg
               ti1 = ttarg
            endif
c
            write(6,'(a,i4,g16.8,2f9.4,3g14.6)')
     >                 'E2d-9:SOL 22:First:',ir,n1,te1,ti1,v1e2d
     >                  ,vpe2d,vpg
c
         elseif (actswe2d.ne.0.0) then
            n1 = cellvals(ir,1,1)
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
            v1e2d = -cellvals(ir,4,1)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*mb))
            vpg = gtarg(ir,1) / n1
c
            if (forcet.eq.1) then
               ttarg = (((5.0* te1 + 3.5 * ti1) * sqrt(te1+ti1))
     >                      / (8.5 * sqrt(2.0)))**(2.0/3.0)
               te1 = ttarg
               ti1 = ttarg
            endif
c
            write(6,'(a,i4,g16.8,2f9.4,3g14.6)')
     >                 'E2d:SOL 22:First:',ir,n1,te1,ti1,v1e2d
     >                  ,vpe2d,vpg
c
         endif
c
         errcode = 0
         serr    = 0.0
c
         npts = nks(ir) - midnks
c
         do ik = nks(ir), midnks + 1 , -1
c
            spts(nks(ir) - ik + 1) =  ksmaxs2(ir) - kss2(ik,ir)
c
c           write(6,*) 'in2:',ir,ik,nks(ir)-ik+1,
c     >                   spts(nks(ir)-ik+1)
c
            fact = recfrac
c
            if (actswmajr.eq.2.0) then
               fact = rs(ik,ir)/r0 * fact
            elseif (actswmajr.eq.3.0) then
               fact = r0/rs(ik,ir) * fact
            endif
c
            if ((actswioni.eq.11.or.actswioni.eq.15)
     >          .and.(.not.pinavail)) then
               ionsrc(nks(ir)-ik+1) = e2dion(ik,ir) * fact
            else
               ionsrc(nks(ir)-ik+1) = pinion(ik,ir) * fact
            endif
c
            if (actswrecom.eq.2.and.(.not.pinavail)) then
               recsrc(nks(ir)-ik+1) = e2dhrec(ik,ir)
            else
               recsrc(nks(ir)-ik+1) = pinrec(ik,ir)
            endif
c
            gperp(nks(ir)-ik+1)  = gperpa(ik,ir)
c
            if (oldktebs(ik,ir).lt.tcutqe) then
               qesrc(nks(ir)-ik+1)  = 0.0
               write(6,'(a5,i4,2(1x,g13.5))') 'QE:',ik,
     >                         oldktebs(ik,ir),tcutqe
            else
               qesrc(nks(ir)-ik+1)  = -pinqe(ik,ir) * fact
     >                                 * qesrc_mult
            endif
c
            if (oldktibs(ik,ir).lt.tcutcx) then
               qisrc(nks(ir)-ik+1)  = 0.0
               write(6,'(a5,i4,2(1x,g13.5))') 'QI:',ik,
     >                         oldktibs(ik,ir),tcutcx
            else
               qisrc(nks(ir)-ik+1)  = -pinqi(ik,ir) * fact
            endif
c
            radsrc(ik) = div_cool(ik,ir)
c
            nhs(nks(ir)-ik+1) = PINATOM(IK,IR)
            nh2s(nks(ir)-ik+1) = PINMOL(ik,ir)
            ths(nks(ir)-ik+1) = pinena(ik,ir) * 0.6666
            nhs0(nks(ir)-ik+1) = pinvdist(1,1,ik,ir)
     >                  + pinvdist(2,1,ik,ir) + pinvdist(3,1,ik,ir)
c
c           Load values from last iteration for PINQID calculations
c
            oldne(nks(ir)-ik+1) = oldknbs(ik,ir)
            oldte(nks(ir)-ik+1) = oldktebs(ik,ir)
            oldti(nks(ir)-ik+1) = oldktibs(ik,ir)
c
c
c           Note that the rbnd array is shifted by one - it starts
c           with a 0 index.
c
            rbnd(nks(ir)-ik) = krb(ik,ir)
            sbnd(nks(ir)-ik) = ksmaxs2(ir) - ksb(ik,ir)
c
            momsrc(nks(ir)-ik+1) = - pinmp(ik,ir) * fact
c
         end do
c
c        Add last rbnd and sbnd elements.
c
c
         rbnd(nks(ir)-midnks) = krb(midnks,ir)
         sbnd(nks(ir)-midnks) = ksmaxs2(ir) - ksb(midnks,ir)
c
c        The end of the Spts array can either be the cell boundary
c        or the mid-point of the ring. The cell boundary will result
c        in more correct integrations since the contribution of the 
c        cells near the mid-point will be calculated more accurately.
c        (Though the difference would be expected to be small.)
c        However, the ring midpoint (0.5 * ringlen) makes it easier to 
c        scale or extend sources over the same spatial region by allowing
c        input to be specified proportional to the ring length. 
c
         spts(npts+1) = halfringlen
c
c        Find an estimated S-value of the X-point
c
         sxp = 0.0
c
         do ik = nks(ir), midnks + 1 , -1
            if (zs(nks(ir),ir).lt.zxp) then
               if (zs(ik,ir).ge.zxp) then
                  sxp = (kss2(ik+1,ir) - kss2(ik,ir)) *
     >                  (zxp-zs(ik+1,ir))
     >                  /(zs(ik,ir)-zs(ik+1,ir))
     >                  + (ksmaxs2(ir) - kss2(ik+1,ir))
                  goto 200
               endif
            elseif (zs(nks(ir),ir).gt.zxp) then
               if (zs(ik,ir).lt.zxp) then
                  sxp = (kss2(ik+1,ir) - kss2(ik,ir)) *
     >                  (zxp-zs(ik+1,ir))
     >                  /(zs(ik,ir)-zs(ik+1,ir))
     >                  + (ksmaxs2(ir)-kss2(ik+1,ir))
                  goto 200
               endif
            endif
         end do
200      continue
c
         write (6,*) 'Sxp:',ir,nks(ir),kss2(midnks,ir),
     >               sxp,ringlen
c
c
c        Set both target power values
c
c        Note:  pae_end and pai_end have
c        a "-" ve because a "-" sign is used in the embedded power
c        equations.
c
         pae_start = -elecptarg(ir,1)
         pai_start = -ionptarg(ir,1)
c
         pae_end = -elecptarg(ir,2)
         pai_end = -ionptarg(ir,2)
c
c
c        Private Plasma target power loss re-distribution to main SOL
c
c        Turned off in private plasma. Electron then ion.
c
         if (actswppelec.eq.0.0.or.pplasma.eq.1) then

            ppelecpow = 0.0

         elseif (actswppelec.eq.1.0.or.actswppelec.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppelec = 0.0
               ppelecpow   = 0.0
c
            else
c
               ppelecpow = -elecptarg(ircor,1)
c
            endif
c
         endif
c
c        PP Ion power loss term
c
         if (actswppion.eq.0.0.or.pplasma.eq.1) then

            ppionpow = 0.0

         elseif (actswppion.eq.1.0.or.actswppion.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppion = 0.0
               ppionpow   = 0.0
c
            else
c
               ppionpow = -ionptarg(ircor,1)
c
            endif
c
         endif
c
c        PP pressure loss term
c
         if (actswppress.eq.0.0.or.pplasma.eq.1) then

            pp_press = 0.0

         elseif (actswppress.eq.1.0.or.actswppress.eq.2.0) then
c
c           Determine corresponding PP ring to current ring.
c
            ircor = nrs - (ir-irsep)
c
            if (ircor.le.irtrap) then
c
c              Turn option off and set power to zero
c
               actswppress = 0.0
               pp_press   = 0.0
c
            else
c
               pp_press = presstarg(ircor,1)
c
            endif

            write(0,*) 'Press3:',actswppress,ir,ircor,pp_press

c
         endif
c
c
c        Set up the ionization source - if choice option is selected.
c
         if (actswion.eq.5.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.5.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1369
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.5) then
                   actswion = 3.0
                else
                   actswioni= 3.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 13.0 - 10.0 * te0
                   ssrcfi = ssrcst + 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = ssrcst + 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.5) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion5I:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
c
         elseif ((actswion.eq.6.0.or.actswion.eq.9.0)
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.(actswioni.eq.6.0.or.actswioni.eq.9)
     >          .and.(.not.pinavail))) then
c
c           Adjust parameters for S**5 Gaussian ionization option.
c
c           This ionization source always starts at zero.
c           Thus the quantity in the input for the start position
c           is used instead as the width factor for the distribution.
c           If the Start point is listed as 0.0 then this is adjusted
c           to 1.0 for the width factor for the gaussian.
c
c           This code only needs to be executed once on each ring.
c
            if (ssrcdecay.eq.0.0) then
               s5gausslen = 1.0
            else
               s5gausslen = ssrcdecay
            endif
c
            ssrcst = 0.0
c
            ssrclen = ssrcfi - ssrcst
c
         elseif (actswion.eq.7.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.7.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1372
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.7) then
                   actswion = 6.0
                else
                   actswioni= 6.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 14.0 - 10.0 * te0
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.0
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.7) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion7I:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
         elseif (actswion.eq.10.0
     >      .or.((actswion.eq.1.0.or.actswion.eq.2.0.or.actswion.eq.8.0)
     >          .and.actswioni.eq.10.0.and.(.not.pinavail))) then
c
c            Analyse target conditions and select appropriate analytic
c            ionization option
c
c            From TN1372
c
             if (n0.gt.1.0e19) then
c
                if (actswion.eq.10.0) then
                   actswion = 9.0
                else
                   actswioni= 9.0
                endif
c
                if (te0.le.1.3) then
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5* (14.0 - 10.0 * te0)
                else
                   ssrcst = 0.0
                   ssrcfi = halfringlen
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                   s5gausslen = 1.5
                endif
c
             elseif (n0.le.1.0e19) then
c
                if (actswion.eq.10.0) then
                   actswion = 4.0
                else
                   actswioni= 4.0
                endif
c
                if (te0.le.10.0) then
                   ssrcst = 0.0
                   ssrcfi = 13.0 - te0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                else
                   ssrcst = 0.0
                   ssrcfi = 2.0
                   ssrclen = ssrcfi - ssrcst
                   ssrcmid = (ssrcfi + ssrcst) / 2.0
                endif
             endif
c
             write(6,'(a,6g13.5)') 'Ion10I:',actswion,
     >                         actswioni,
     >                         ssrcst,ssrcfi,ssrcmid,ssrclen
c
         endif
c
c        Set up extra perpenicular source parameters
c
         if (switch(swextra).gt.0.0) then
c
c           Calculate source start and stop positions - for
c           first 1/2 ring use the parameters as specified.
c
c           Note: These stop and start positions may be
c           on either half of the flux tube or even overlap
c           the top. The integrations procedd outward from
c           each target so the values of the start and
c           stop positions will have to be adjusted but
c           the rest will remain the same.
c
c           For the second half ring - looking from the other
c           target - the positions of "start" and "stop" are
c           reversed.
c
c           e.g.               SOURCE
c                  start  stop
c           0        |      |                          SMAX
c           |-------------------------------------------|
c                                   |         |
c                                 start     stop
c                               SINK
c
c           But for S measured from the SMAX end of the
c           ring  SINK_stop becomes the start of the SINK
c           region and so on for the rest of the positions.
c
c
            start_gextra_src = (1.0-gextra_src_stop) * ringlen
            stop_gextra_src  = (1.0-gextra_src_start)  * ringlen
            start_gextra_sink= (1.0-gextra_sink_stop)* ringlen
            stop_gextra_sink = (1.0-gextra_sink_start) * ringlen
c
c           Calculate extra source strength
c
            if (start_gextra_src.ne.stop_gextra_src) then
               gextra_src = abs(gtarg(ir,3)) * gextra_mult /
     >                   abs(start_gextra_src-stop_gextra_src)
            else
               gextra_src = abs(gtarg(ir,3)) * gextra_mult
            endif
c
c           Calculate extra sink strength
c
            if (start_gextra_sink.ne.stop_gextra_sink) then
               gextra_sink = -abs(gtarg(ir,3)) * gextra_mult /
     >                   abs(start_gextra_sink-stop_gextra_sink)
            else
               gextra_sink = -abs(gtarg(ir,3)) * gextra_mult
            endif
c
            write (6,*) 'Gextra 1:',ir,gextra_src,gextra_sink,
     >              gtarg(ir,3),ringlen,
     >              start_gextra_src,stop_gextra_src,
     >              start_gextra_sink,stop_gextra_sink
c
         endif
c
c        Error correcting branch point
c
 250     continue
c
c        Initialize output arrays to default error values
c
         call dinit(te,mxspts,10.0d0)
         call dinit(ti,mxspts,10.0d0)
         call dinit(ne,mxspts,1.0d19)
         call dinit(vb,mxspts,0.0d0)
c
         write (6,'(a,i4,a)') '************   START SECOND HALF ',
     >               ir,'  ***********'
         write (6,*) 'Next :',ir,te0,ti0,n0,npts,nks(ir)
c
c        Zero power ratio information in case calcsol is not called.
c
         call qzero (int_powrat,3)
c
c        Override with detached plasma solution for the inner
c
         if (switch(swdetach).eq.2.0) then
c
            call detached_plasma(spts,npts,errcode,serr,
     >                           te,ti,ne,vb,ir,switch(swdetach),
     >                           te0,ti0,n0,v0,act_press,mb)
c
c        Solve using SOL22
c
         else
c
c slmod begin - new
c...Need a system routine for the other system files.  Rename 
c   CLOCK2:
            IF (osm_mode.EQ.2) THEN
              CALL OpenStorageFiles(ir,IKLO,'tmp.dat')
              deltat = Clock2()
            ENDIF
c slmod end
            call calcsol (spts,npts,errcode,serr,
     >                 te,ti,ne,vb,exp_press,act_press,
     >                 prad,ir,irlim1,
     >                 int_powrat,cprint)
c slmod begin - new
            IF (osm_mode.EQ.2) THEN
              CALL SOL22Status(IKHI,ir,deltat,serr,spts,npts,errcode)
            ENDIF
c slmod end
c
         endif
c
c        Handle Error Condition if Error Switch is set
c
         if ((   actswerror.ge.1.0
     >            .and.(errcode.eq.3.or.errcode.eq.4.or.
     >                  errcode.eq.5.or.errcode.eq.6.or.
     >                  errcode.eq.7))
     >       .or.(actswerror.eq.0.0.and.
     >           (errcode.eq.6.or.errcode.eq.7))) then
c
c           If a negative N error occurs without error correction -
c           Turn error correction ON.
c
            if (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7)
     >          .and.seterror.eq.0) then
               actswerror = 5.0
               seterror = 1
            elseif (actswerror.eq.0.0.and.(errcode.eq.6.or.errcode.eq.7)
     >              .and.seterror.eq.1) then
                seterror = 2
            elseif (seterror.eq.2) then
c
c               ERROR - negative N has been found EVEN with highest level
c                       of error correction - issue error messages and stop.
c
                write (6,*) 'SOLASCV: SOL22:'//
     >               ' Unsolvable Negative N error encountered'
                write (6,*) 'Program will STOP'
c
                write (7,*) 'SOLASCV: SOL22:'//
     >               ' Unsolvable Negative N error encountered'
                write (7,*) 'Program will STOP'
c
                stop
c
            endif
c
c           Record error
c
            cdeferr(ir,1) = errcode
            cdefserr(ir,1) = serr
c
c           Set switches
c
            call setsw(-2,pplasma,new_errlevel)
c
c           Record error level of next attempt
c
c            cdeferropt(ir,1) = actswerror
c
            cdeferropt(ir,1) = new_errlevel
c
c
            call initlen
c
            write (6,*) 'Error Handler: SET INNER:',cdeferr(ir,1),
     >                   cdefserr(ir,1),cdeferropt(ir,1)

c
c           Re-do Inner 1/2 ring
c
            goto 250
c
         endif
c
c        Save background density and temperature
c
         knds(idds(ir,1)) =  n0
         kvds(idds(ir,1)) = -v0
c
c        If E2D option 9 was in use - fill in the missing
c        plasma with a linear fit.
c
         if (actswe2d.eq.9.0.and.switch(swdetach).ne.2.0) then
c
c           Fill in missing values - depending on option specified.
c
            if (fillopt.eq.0) then
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0 + (te(ike2d_start)-te0) *
     >                      spts(ik)/spts(ike2d_start)
                  ti(ik) = ti0 + (ti(ike2d_start)-ti0) *
     >                      spts(ik)/spts(ike2d_start)
                  ne(ik) = n0 + (ne(ike2d_start)-n0) *
     >                      spts(ik)/spts(ike2d_start)
                  vb(ik) = v0 + (vb(ike2d_start)-v0) *
     >                      spts(ik)/spts(ike2d_start)
               end do
c
            elseif (fillopt.eq.1) then
c
               gam = n0 * v0
c
c              Extrapolate target electron temperature
c
               slope = (te(ike2d_start+1) - te(ike2d_start)) /
     >                 (spts(ike2d_start+1) - spts(ike2d_start))
c
               te0 = - spts(ike2d_start) * slope + te(ike2d_start)
c
               if (te0.le.0.0) te0 = te(ike2d_start)
c
c              Extrapolate target ion temperature
c
               slope = (ti(ike2d_start+1) - ti(ike2d_start)) /
     >                 (spts(ike2d_start+1) - spts(ike2d_start))
c
               ti0 = - spts(ike2d_start) * slope + ti(ike2d_start)
c
               if (ti0.le.0.0) ti0 = ti(ike2d_start)
c
               v0 = -sqrt((te0+ti0)/crmb * econv/mconv)
c
               n0 = gam/v0
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0 + (te(ike2d_start)-te0) *
     >                      spts(ik)/spts(ike2d_start)
                  ti(ik) = ti0 + (ti(ike2d_start)-ti0) *
     >                      spts(ik)/spts(ike2d_start)
                  ne(ik) = n0 + (ne(ike2d_start)-n0) *
     >                      spts(ik)/spts(ike2d_start)
                  vb(ik) = v0 + (vb(ike2d_start)-v0) *
     >                      spts(ik)/spts(ike2d_start)
               end do
c
c              Reset target values
c
               knds(idds(ir,1)) = n0
               kvds(idds(ir,1)) = -v0
               kteds(idds(ir,1)) = te0
               ktids(idds(ir,1)) = ti0
c
            elseif (fillopt.eq.2) then
c
c              Fill with constant values from last point
c
c              Target values are NOT affected
c
               do ik = 1,ike2d_start-1
                  te(ik) = te(ike2d_start)
                  ti(ik) = ti(ike2d_start)
                  ne(ik) = ne(ike2d_start)
                  vb(ik) = vb(ike2d_start)
               end do
c
            elseif (fillopt.eq.3) then
c
c              Fill with values held constant at target conditions
c
               do ik = 1,ike2d_start-1
                  te(ik) = te0
                  ti(ik) = ti0
                  ne(ik) = n0
                  vb(ik) = v0
               end do
c
            endif
c
         endif



c
c        Assign Background
c
         do ik = nks(ir), midnks + 1 , -1
c
            ktebs(ik,ir) = te(nks(ir) - ik + 1)
            ktibs(ik,ir) = ti(nks(ir) - ik + 1)
            knbs(ik,ir) = ne(nks(ir) - ik + 1)
            oldknbs(ik,ir) = ne(nks(ir) - ik + 1)
            oldktibs(ik,ir) = ti(nks(ir) - ik + 1)
            oldktebs(ik,ir) = te(nks(ir) - ik + 1)
            kvhs(ik,ir) = -vb(nks(ir) - ik + 1)
c slmod begin - new
c            oldkvhs(ik,ir) = -vb(nks(ir) - ik + 1)

            osmcde(ik,ir) = -cde(nks(ir)-ik+1)
            osmcdi(ik,ir) = -cdi(nks(ir)-ik+1)
            osmcve(ik,ir) = -cve(nks(ir)-ik+1)
            osmcvi(ik,ir) = -cvi(nks(ir)-ik+1)
c slmod end
c
            kpress(ik,ir,1) = exp_press(nks(ir) - ik + 1)
            kpress(ik,ir,2) = act_press(nks(ir) - ik + 1)
            kprad(ik,ir)    = prad(nks(ir)-ik+1)
c
         end do
c slmod begin - new
c... .TRUE. temporary:
         IF (.NOT.pinavail.AND.(rel_opt.EQ.1.OR.rel_opt.EQ.3))
     .     CALL CalcInitSrc(IKHI,ir)
c slmod end
c
c        Save mach numbers and error codes
c
         cmachno(ir,1) = m0
         cerr(ir,1) = errcode
         cserr(ir,1) = serr
c
c        Save Power Ratios
c
         sol22_power_ratio(ir,1,1) = int_powrat(1)
         sol22_power_ratio(ir,1,2) = int_powrat(2)
         sol22_power_ratio(ir,1,3) = int_powrat(3)
c
         write(6,'(a,i4,3(1x,g12.5))') 'PR1:',ir,int_powrat(1),
     >                    int_powrat(2),int_powrat(3)
c
c        End of 1/2 ring selection
c

c           The parameters are current ring snd local IKOPT
c           data will be printed if it has been collected
            if (debug_sol22.ne.0) call check_print_data(ir,2)


         endif


         call pr_trace('CALCSOL_INTERFACE','CALCULATE EFIELD')

C
C       CALCULATE ELECTRIC FIELD
C
C
C       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
C       SAME FACTOR USED IN CONVERTING T IN EV TO KT.
C
c       If OFIELD is turned ON then set the electric field to
c       zero for this case. Note: the electric field is
c       initialized to zero in the plasma.d3a module.
c
c       For OFIELD ... 0=off   1=on
c
        if (ofield.eq.0) then
c
c
          if (ikopt.eq.1.or.ikopt.eq.3) then


          if (kss2(1,ir).eq.0.0) then
            DS1 = KSS2(2,IR) - KSS2(1,IR)
            DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT1 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
C
            KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,2))= kes(1,ir)
          else
            DS2 = KSS2(2,IR) - KSS2(1,IR)
            DP2 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT2 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB2 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
c
            DS1 = KSS2(1,IR)
            DP1 = KNBS(1,IR)*KTEBS(1,IR)-
     >              KNDS(idds(ir,2))*KTEDS(idds(ir,2))
            DT1 = KTEBS(1,IR)-KTEDS(idds(ir,2))
            NB1 = 0.5*(KNBS(1,IR)+KNDS(idds(ir,2)))
c
            KES(1,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,2)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
c
          endif
c
C
          if (ikopt.eq.2.or.ikopt.eq.3) then
c
          if (kss2(nks(ir),ir).eq.ksmaxs2(ir)) then
            DS1 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
            DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
C
            KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
c
            KEDS(idds(ir,1))= kes(nks(ir),ir)
          else
c
            DS2 = KSS2(NKS(IR),IR) - KSS2(NKS(IR)-1,IR)
            DP2 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT2 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB2 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
c
            DS1 = ksmaxs2(ir) - KSS2(nks(ir),IR)
            DP1 = KNDS(idds(ir,1))*KTEDS(idds(ir,1))-
     >             KNBS(nks(ir),IR)*KTEBS(nks(ir),IR)
            DT1 = KTEDS(idds(ir,1))-KTEBS(nks(ir),IR)
            NB1 = 0.5*(KNBS(nks(ir),IR)+KNDS(idds(ir,1)))
c
            KES(nks(ir),IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
c
            KEDS(idds(ir,1)) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
          endif
c
          endif
C
          if (ikopt.eq.1) then
             ikfirst = ikstart +1
             iklast = ikmids(ir)
          elseif (ikopt.eq.2) then
             ikfirst = ikmids(ir) + 1
             iklast = ikend - 1
          elseif (ikopt.eq.3) then
             ikfirst = ikstart +1
             iklast = ikend -1
          endif
c
          DO 500 IK = ikfirst,iklast

            DS1 = KSS2(IK,IR) - KSS2(IK-1,IR)
            DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
            DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
            NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
            DS2 = KSS2(IK+1,IR) - KSS2(IK,IR)
            DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
            DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
            NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
            KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
C
C            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)
C
 500      CONTINUE
c
c         End of test on OFIELD
c

          write(6,*) 'KES:',keds(idds(ir,2)),keds(idds(ir,1)),
     >                      idds(ir,2),idds(ir,1),ir
          do ik =1,nks(ir)
             write (6,*) ik,kes(ik,ir)
          end do

c
        endif


c
c       If smoothing is turned ON
c
        if (switch(swsmooth).eq.1.0) then
c
c        WF'95:    CORRECTING(?) FOR PLATE ASSYMETRY
c
c        profile scaling of ne,Te,Ti to match midpoint values,
c        exponent asmexp determines smoothing range
c
         TEMID = 0.5*(KTEBS(MIDNKS+1,IR) + KTEBS(MIDNKS,IR))
         TIMID = 0.5*(KTIBS(MIDNKS+1,IR) + KTIBS(MIDNKS,IR))
         NMID = 0.5*(KNBS(MIDNKS+1,IR) + KNBS(MIDNKS,IR))
c         VMID = 0.5*(KVHS(MIDNKS+1,IR) + KVHS(MIDNKS,IR))
         SMAX = KSMAXS2(IR)
         ASMEXP = 1.0
c         write(6,*) 'assymetry', ir,ktibs(midnks,ir),
c     >   ktibs(midnks+1,ir), timid
c
        DO IK = ikstart , ikend
         IF (IK.LE.MIDNKS) THEN
           TMP = (KSS2(IK,IR) / (SMAX / 2.0))**ASMEXP
           KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +
     >        TMP*(TEMID/KTEBS(MIDNKS,IR) - 1.0))
           KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +
     >        TMP*(TIMID/KTIBS(MIDNKS,IR) - 1.0))
           KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +
     >        TMP*(NMID/KNBS(MIDNKS,IR) - 1.0))
c           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
c     >        TMP*(VMID/KVHS(MIDNKS,IR) - 1.0))
         ELSE
          TMP = ((SMAX-KSS2(IK,IR)) / (SMAX / 2.0))**ASMEXP
           KTEBS(IK,IR) = KTEBS(IK,IR)*(1.0 +
     >        TMP*(TEMID/KTEBS(MIDNKS+1,IR) - 1.0))
           KTIBS(IK,IR) = KTIBS(IK,IR)*(1.0 +
     >        TMP*(TIMID/KTIBS(MIDNKS+1,IR) - 1.0))
           KNBS(IK,IR) = KNBS(IK,IR)*(1.0 +
     >        TMP*(NMID/KNBS(MIDNKS+1,IR) - 1.0))
c           KVHS(IK,IR) = KVHS(IK,IR)*(1.0 +
c     >        TMP*(VMID/KVHS(MIDNKS+1,IR) - 1.0))
         ENDIF
        ENDDO
c
c      End of Smoothing IF
c
       endif
c
c      Branch location for virtual rings
c
 1000  continue
c
c     End of DO loop for IR
c
      end do
c
c     Reset the ionization options for swion=8.0
c
      if (switch(swioni).eq.8.0.and.pinavail) then
         switch(swioni) = switch(swion)
         switch(swion) = 8.0
      endif
c
c     End interface routine
c
      call pr_trace('CALCSOL_INTERFACE','BEFORE EXCEL PRINT')
c
c     Call routine to print out table of values to file for Excel plotting
c     (req. by Wojciech Fundamenski)
c
      call print_sol_excel
c
      return
      end
c
c
c
      real*8 function find_ffric(ir,targid,actmomlen)
      implicit none
      integer ir,targid
      include 'params'
      include 'solparams' 
      include 'solcommon'
      real*8 :: actmomlen
c
c     Assign the actual FFRIC values for each ring from the data contained in 
c     the inputs FFRIC and EXTFFRIC. This option allows for a different
c     amount of momentum loss to be specified for every half flux tube.
c
      integer in
c
      find_ffric = ffric
      actmomlen = lenmom
c
      if (n_extffric.eq.0) return        
c
      do in = 1,n_extffric
c     
         if (ir.eq.int(extffric(in,1))) then 
c     
            if (targid.eq.2) then 
               if (extffric(in,2).gt.0.0) then 
                  find_ffric = extffric(in,2)
               endif
               if (extffric(in,3).gt.0.0) then 
                  actmomlen = extffric(in,3)
               endif
            elseif(targid.eq.1) then 
               if (extffric(in,4).gt.0.0) then 
                  find_ffric = extffric(in,4)
               endif
               if (extffric(in,5).gt.0.0) then 
                  actmomlen = extffric(in,5)
               endif
            endif  
             
            return
      
         endif
c     
      end do
c
      return
      end 

