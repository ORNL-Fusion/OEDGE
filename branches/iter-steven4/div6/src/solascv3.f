c
c
c
      subroutine readsol(ierr)
      implicit none
c
      integer ierr
c
c     This subroutine reads the input parameters for the case
c     from the standard input or redirected from a file.
c
      include 'params'
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8 spts(mxspts)
c
c      integer ierr
c
c     Model parameters
c
      call rdi(forcet,.TRUE.,0,.true.,1,         'force te=ti'    ,ierr)
      CALL RDQ(initm0,.TRUE.,0.0d0,.FALSE.,0.0d0,'target mach num',IERR)
      CALL RDQ(deltam0,.TRUE.,0.0d0,.FALSE.,0.0d0,'delta mach num',IERR)
      CALL RDQ(m0res,.TRUE.,0.0d0,.FALSE.,0.0d0,'Resolution in m0',IERR)
c
c     ------------------------------------------------------------------
c
c     Ionization Source
c
      CALL RDI(lensind,.TRUE.,0,.TRUE.,1,     'ion source abs/rel',IERR)
      CALL RDQ(lensst,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src start',  IERR)
      CALL RDQ(lensfi,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src finish', IERR)
      CALL RDQ(lams,.TRUE. ,0.0d0,.FALSE.,0.0d0,'ion decay len   ',IERR)
c
c     Radiation source
c
      CALL RDQ(lenri,.TRUE.,0.0d0,.FALSE.,0.0d0,'rad source len  ',IERR)
      CALL RDQ(lamr,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad decay len   ',IERR)
      CALL RDQ(frr ,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad power mult  ',IERR)
      CALL RDQ(alfimp,.TRUE.,0.0d0,.FALSE.,0.0d0,'nimp/ne ratio  ',IERR)
      CALL RDQ(talimp,.TRUE. ,0.0d0,.FALSE.,0.0d0,'base Temp     ',IERR)
      CALL RDQ(ex1imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'ecponent 1  ',IERR)
      CALL RDQ(ex2imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 2  ',IERR)
c
c     Miscellaneous
c
      CALL RDQ(gamcor,.false.,0.0d0,.FALSE.,0.0d0,'i power corr. ',IERR)
      CALL RDQ(gamecor,.false.,0.0d0,.FALSE.,0.0d0,'e power corr.',IERR)
c
      CALL RDQ(ceicf,.TRUE. ,0.0d0,.FALSE.,0.0d0, 'CX power frac ',IERR)
      CALL RDQ(recfrac,.TRUE.,0.0d0,.TRUE.,1.0d0,'Recycle frac ',IERR)
c
      call rdq(peicf,.true.  ,0.0d0,.false.,0.0d0,'Pei Correction',ierr)
      call rdi(velsw,.true.  ,0,.true.,3,       'Vel Error Switch',ierr)
c
c     Power distribution
c
      call rdq(spowbeg,.true.,0.0d0,.true.,0.5d0,'Power Dist Beg',ierr)
      call rdq(spowlen,.true.,0.0d0,.true.,0.5d0,'Power Dist Len',ierr)
c
c     Gperp Distribution function
c
      call rdq(gperpfrac,.true.,0.0d0,.true.,1.0d0,
     >                                            'Part Dist Frac',ierr)
      call rdq(gperpbegf,.true.,0.0d0,.true.,0.5d0,'Part Dist Beg',ierr)
      call rdq(gperpendf,.true.,0.0d0,.true.,0.5d0,'Part Dist Len',ierr)
c
c     Extra Gperp source/sink - start and end positions.
c
      call rdq(gextra_mult,.true.,0.0d0,.false.,0.0d0,
     >                                         'Gextra flux mult',ierr)
      call rdq2(gextra_src_start,gextra_src_stop,
     >         .true.,0.0d0,.true.,1.0d0, 'Gextra SRC start/stop',ierr)
      call rdq2(gextra_sink_start,gextra_sink_stop,
     >         .true.,0.0d0,.true.,1.0d0,'Gextra SINK start/stop',ierr)
c
c     Field line length fraction for distributing the Private plasma
c     electron and ion power loads.
c
      call rdq(pp_pow_dist,.true.,0.0d0,.false.,0.0d0,'PP Pow Dist',
     >         ierr)
c
c     IK index for edge2d compatibility option 9
c
      call rdi(ike2d,.true.,1,.false.,0,'IK start Index for E2D-9',ierr)
      call rdi(fillopt,.true.,0,.true.,3,'Gap fill option - E2D-9',ierr)
c
c     Cutoff temperature for PINQID term
c
      call rdq(tcutqe,.true.,0.0d0,.false.,0.0d0,'Cut T-PINQE',ierr)
      call rdq(tcutatiz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDATIZ',ierr)
      call rdq(tcutmliz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDMLIZ',ierr)
      call rdq(tcutrec,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDREC',ierr)
      call rdq(tcutcx,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDCX',ierr)
      call rdq(trefcx,.true.,0.0d0,.false.,0.0d0,'REF T- QIDCX 1',ierr)
      call rdq(tmin,.false.,0.0d0,.false.,0.0d0,'Min. Allowed T',ierr)
      call rdq(dropfrac,.true.,0.0d0,.true.,1.0d0,'Allowed T-drop',ierr)
c
c     Momentum Source
c
      call rdq(smom_mult,.false.,0.0d0,.false.,0.0d0,
     >                                       'Mom.Loss Multiplier',ierr)
      call rdq(ffric,.false.,0.0d0,.false.,0.0d0,'Mom.loss frac  ',ierr)
      CALL RDQ(lenmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom source len',IERR)
      CALL RDQ(lammom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom decay len ',IERR)
      CALL RDQ(rcxmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'cx/iz ratio   ',IERR)
      call rdq(tcxmom,.TRUE.,1.0001d0,.FALSE.,0.0d0,'T for CXmult',ierr)
      call rdq(tcxcut,.true.,0.0d0,.false.,0.0d0,'Cut T -  CXmult',ierr)
c
c     Source term multipliers
c
      call rdq(qesrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE mult',ierr)
      call rdq(radsrc_mult,.false.,0.0d0,.FALSE.,0.0d0,
     >                                     'PINQE based PRAD mult',ierr)
c
c     call rdq(qisrc_mult,.false.,0.0d0,.false.,0.0d0,'PINQI mult',ierr)
c
c     ------------------------------------------------------------------
c
      CALL RDI(ndiv,.TRUE., 1,.FALSE., 0,'NUMBER OF STEPS    ',IERR)
c
c     Read in switches 0.0 is off, 1.0 is on.
c
      call rdr(switch(swion),.true. ,0.0,.false.,0.0, 'alt ion',ierr)
      call rdr(switch(swioni),.true. ,0.0,.false.,0.0,'init ion',ierr)
      call rdr(switch(swionp),.false.,0.0,.false.,0.0,'pp ion',ierr)
c
      if (switch(swionp).eq.-1.0) then
         if (switch(swion).eq.1.0.or.switch(swion).eq.2.0
     >       .or.switch(swion).eq.8.0) then
            switch(swionp) = switch(swioni)
         else
            switch(swionp) = switch(swion)
         endif
      endif
c
      call rdr(switch(swcond),.true.,0.0,.false.,0.0, 'cond sw',ierr)
      call rdr(switch(swconv),.true.,0.0,.false.,0.0, 'conv sw',ierr)
      call rdr(switch(swprad),.true.,0.0,.false.,0.0, 'prad sw',ierr)
      call rdr(switch(swphelp),.true.,0.0,.false.,0.0,'phelp sw',ierr)
      call rdr(switch(swpei),.true.,0.0,.false.,0.0,  'pei sw ',ierr)
      call rdr(switch(swpcx),.true.,0.0,.false.,0.0,  'pcx sw ',ierr)
c
c     PINQID switches
c
      call rdr(switch(swqidatiz),.true.,0.0,.false.,0.0,'atiz sw',ierr)
      call rdr(switch(swqidmliz),.true.,0.0,.false.,0.0,'mliz sw',ierr)
      call rdr(switch(swqidrec),.true.,0.0,.false.,0.0,'rec sw',ierr)
      call rdr(switch(swqidcx),.true.,0.0,.false.,0.0,'cx sw',ierr)
c
      call rdr(switch(swppelec),.true.,0.0,.false.,0.0,
     >                                               'pp elec sw',ierr)
      call rdr(switch(swppion),.true.,0.0,.false.,0.0,'pp ion sw',ierr)
c
c
c     jdemod - 
c     switch(swppress) is read in using tag 283 of unstructured input 
c     default value is 0 or OFF
c
c      call rdr(switch(swppress),.true.,0.0,.false.,0.0,'pp power sw',ierr)
c
c
      call rdr(switch(swvisc1),.true.,0.0,.true.,0.0,'visc1 sw',ierr)
      call rdr(switch(swnmom),.true.,0.0,.false.,0.0, 'N mom sw',ierr)
      call rdr(switch(swmach),.true.,0.0,.false.,0.0, 'mach sw',ierr)
c
c     Read in Edge2d compatibility switch and the subsequent values of
c     ne, Te and Ti at the centre point of the first cell
c
      call rdr(switch(swe2d),.false.,0.0,.false.,0.0, 'e2d sw',ierr)
      call rdr(switch(swpow),.true.,0.0,.false.,0.0, 'power sw',ierr)
      call rdr(switch(swpowp),.false.,0.0,.false.,0.0, 'pp pow',ierr)
c
      if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)
c
      call rdr(switch(swgperp),.true.,0.0,.false.,0.0,'GamPerp',ierr)
      call rdr(switch(swgperpp),.false.,0.0,.false.,0.0,'GamPerpP',ierr)
c
      if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)
c
c     Extra Gperp source/sink term
c
      call rdr(switch(swextra),.true.,0.0,.false.,0.0,'GP Src/Sink',
     >                                                          ierr)
c
      call rdr(switch(swmajr),.true.,0.0,.false.,0.0,'MajorRad',ierr)
      call rdr(switch(swcore),.true.,0.0,.false.,0.0,'Core Src',ierr)
      call rdr(switch(swrecom),.true.,0.0,.false.,0.0,'Recomb.',ierr)
      call rdr(switch(swsmooth),.true.,0.0,.false.,0.0,'Smooth',ierr)
      call rdr(switch(swdetach),.true.,0.0,.false.,0.0,'Detach',ierr)
      call rdr(switch(swerror),.true.,0.0,.false.,0.0,'ERROR',ierr)
      CALL RDRARN(deflist,ndef,mxspts,0.0,real(maxnrs),.FALSE.,
     >            0.0,MACHHI,2,'DEFAULT SOLVER DATA',IERR)
c
c     Default plots to off - this can be changed in the calcsol_interface
c     routine in solascv1.f
c
      graph = 0
      graphaux = 0
      graphvel = 0
      graphran = 0.0
c
c      CALL RDI(graph,.TRUE., 0,.TRUE., 1,'GRAPH OPTION       ',IERR)
c      CALL RDI(graphaux,.TRUE.,0,.TRUE.,1,'AUX GRAPH OPTION  ',IERR)
c      CALL RDI(graphvel,.TRUE.,0,.TRUE.,1,'VEL GRAPH OPTION  ',IERR)
c      call rdr(graphran,.true.,0.0,.false.,0.0,'CXMAX VALUE  ',ierr)
c
c
      return
      end
c
c
c
      subroutine echosol
      implicit none
c
c     This subroutine echoes the input values to standard out - it also
c     prints the additional calculated values - and is followed by the
c     tabular output of s,te,ti,n and v at each point S ... which would
c     be suitable for plotting on a spreadsheet or may be plotted by
c     calling GHOST routines if the graph option is set.
c
      include 'params'
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      include 'printopt'
c
      include 'cgeom'
      include 'comtor'
c slmod begin - new
      INCLUDE 'slcom'

      INTEGER fp,i1,i2
c slmod end
c
      integer i,irlim,ir,cnt,in
      real*8 spts(mxspts)
      character*20 errstr(9)
c
      real swtmp
c
c      character*24 sp
c      character  s1*6,s2*12
c
c
c      CHARACTER  COMENT*80
c
c     Initialization
c
c      sp = '                        '
c      s1 = '      '
c      s2 = '            '
c
      errstr(1) = 'Imaginary N FOUND'
      errstr(2) = 'INVALID          '
      errstr(3) = 'Negative T  ERROR'
      errstr(4) = 'Step-size   ERROR'
      errstr(5) = 'Excessive  T Drop'
      errstr(6) = 'Negative N  ERROR'
      errstr(7) = 'NaNQ-Solver ERROR'
      errstr(8) = 'SPECIFIED ERR OPT'
c
      IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >    ciopto.eq.3.or.ciopto.eq.4) THEN
         IRLIM = IRWALL
      ELSEIF (CIOPTO.EQ.1) THEN
         IRLIM = NRS
      ENDIF
c
      call prb
      CALL PRC ('  SOL OPTION 22:  SUB-OPTIONS AND RESULTS SUMMARY')
      call prb
      call prc (s1//'SOL22: SUMMARY OF OPTIONS AND INPUT VALUES')
      call prb
c
c     Options
c
      CALL PRQ (S1//'GAMMA CORRECTION FACTOR IN GAMMAI  ', GAMCOR)
      CALL PRQ (S1//'GAMMA CORRECTION FACTOR IN GAMMAE  ', GAMECOR)
c
      call prb
c
      if (recfrac.eq.1.0) then
         CALL PRC(S1//'RECYCLING SOURCE FRACTION IS OFF')
         CALL PRC(S1//'- RECYCLING SOURCE FRACTION = 1.0')
      elseif (recfrac.ne.1.0) then
         CALL PRC(S1//'RECYCLING SOURCE FRACTION IS ON')
         CALL PRQ(S1//'- RECYCLING SOURCE FRACTION = ',RECFRAC)
      endif
c
      CALL PRB
      CALL PRI (S1//'INITIAL NUMBER OF RUNGE-KUTTA STEPS BETWEEN'//
     >          ' EACH GRID KNOT IN SOLVER:', NDIV)
      CALL PRB
c
c
c     Indicate Forced Te=Ti or NOT
c
      if (forcet.eq.0) then
         CALL PRC (S1//'T FORCE OPTION 0: TE AND TI ARE FOLLOWED'//
     >                               ' INDEPENDENTLY.')
      ELSEIF (FORCET.EQ.1) THEN
         CALL PRC (S1//'T FORCE OPTION 1: TE AND TI ARE LOCKED'//
     >                               ' TOGETHER.')
         CALL PRC (SP//'BASED ON A COMBINED'//
     >                           ' ENERGY EQUATION.')
      endif
c
      call prb
c
c     Velocity Error Correction Option:
c
      if (velsw.eq.0) then
         CALL PRC (S1//'VEL/COR OPT   0 : VELOCITY SET TO CS WHEN'//
     >              ' IMAGINARY RESULT FOUND.')
      elseif (velsw.eq.1) then
         CALL PRC (S1//'VEL/COR OPT   1 : VELOCITY HELD CONSTANT WHEN'//
     >              ' IMAGINARY RESULT FOUND.')
      elseif (velsw.eq.2) then
         CALL PRC (S1//'VEL/COR OPT   2 : PRESSURE SET TO MINIMUM'//
     >              ' NECESSARY TO AVOID IMAGINARY.')
         CALL PRC (SP//'DENSITY IS SET FOR THIS'//
     >              ' PRESSURE VALUE.')
         CALL PRC (SP//'VELOCITY IS SET USING FLUX'//
     >              ' CONSERVATION.')
         CALL PRC (SP//'ADDITIONAL PRESSURE REQUI'//
     >              'RED IS CARRIED FORWARD.')
      elseif (velsw.eq.3) then
         CALL PRC (S1//'VEL/COR OPT   3 : PRESSURE SET TO MINIMUM'//
     >              ' NECESSARY TO AVOID IMAGINARY.')
         CALL PRC (SP//'DENSITY IS SET FOR THIS'//
     >              ' PRESSURE VALUE.')
         CALL PRC (SP//'VELOCITY IS SET USING FLUX'//
     >              ' CONSERVATION.')
         CALL PRC (SP//'ADDITIONAL PRESSURE REQUI'//
     >              'RED IS NOT CARRIED FORWARD.')
      endif
c
      if (lensind.eq.0) then
         CALL PRC (S1//'LENGTH OPTION 0 : IONIZATION LENGTHS ARE IN'//
     >                      ' ABSOLUTE UNITS (M)'
     >                   //' UNLESS')
         call prc (sp//'INDICATED OTHERWISE.')
      elseif (lensind.eq.1) then
         CALL PRC (S1//'LENGTH OPTION 1 : IONIZATION LENGTHS ARE IN'//
     >                    ' RELATIVE UNITS * SMAX'
     >                   //' UNLESS')
         call prc (sp//'INDICATED OTHERWISE.')
      endif
c
c
c     Main Ionization options
c
      call prb
      CALL PRC(S1//'MAIN IONIZATION OPTION: ')
c
      call prb
c
      if (switch(swion).eq.0.0) then
c
        CALL PRB
        CALL PRC (S1//'IONIZATION OPT 0: EXPONENTIAL'//
     >                    ' DECAY IONIZATION SOURCE')
        CALL PRB
C
        CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',
     >                         LENSFI)
        CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',
     >                         LAMS)
c
      elseif (switch(swion).eq.1.0) then
c
        CALL PRC (S1//'IONIZATION OPT 1: PIN DATA READ FOR'//
     >                       ' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS NORMALIZED TO NOVO')
        CALL PRB
        CALL PRC (SP//'DATA FOR IONIZATION SOURCE IS RETURNED')
        CALL PRC (SP//'FROM A PIN/NIMBUS RUN AND IS LINEARLY')
        CALL PRC (SP//'INTERPOLATED FOR VALUES BETWEEN GRID')
        CALL PRC (SP//'POINTS. THE IONIZATION SOURCE IS')
        CALL PRC (SP//'INTEGRATED OVER THE HALF-RING AND SET')
        CALL PRC (SP//'EQUAL TO THE FLOW TO THE TARGET NOVO')
        CALL PRC (SP//'AS A NORMALIZATION FACTOR.')
c
      elseif (switch(swion).eq.2.0) then
c
        CALL PRC (S1//'IONIZATION OPT 2: PIN DATA READ FOR'//
     >                        ' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS UNNORMALIZED.')
        CALL PRB
        CALL PRC (SP//'DATA FOR IONIZATION SOURCE IS RETURNED')
        CALL PRC (SP//'FROM A PIN/NIMBUS RUN AND IS LINEARLY')
        CALL PRC (SP//'INTERPOLATED FOR VALUES BETWEEN GRID')
        CALL PRC (SP//'POINTS. THE DATA IS USED AS-IS AND IS')
        CALL PRC (SP//'NOT NORMALIZED TO THE TARGET FLUX FOR')
        CALL PRC (SP//'THE RING.')
c
      elseif (switch(swion).eq.8.0) then
c
        CALL PRC (S1//'IONIZATION OPT 8: PIN DATA READ FOR'//
     >                    ' IONIZATION SOURCE')
        CALL PRC (SP//'DATA IS UNNORMALIZED')
        CALL PRB
        CALL PRC (SP//'THE INTEGRATED STRENGTH OF THE PIN'//
     >                      ' IONIZATION')
        CALL PRC (SP//'IS USED TO NORMALIZE THE ANALYTIC OPTION')
        CALL PRC (SP//'SPECIFIED FOR THE INITIAL PLASMA ON'//
     >                   ' SUBSEQUENT')
        CALL PRC (SP//'ITERATIONS.')
c
      elseif (switch(swion).eq.3.0) then
c
        CALL PRC(S1//'IONIZATION OPT 3: IMPOSED TRIANGULAR'//
     >                        ' IONIZATION SOURCE')
        CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
        CALL PRQ(SP//'          TO     :',LENSFI)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                        //' TO RING TARGET FLUX.')
        call prb
c
      elseif (switch(swion).eq.4.0) then
c
        CALL PRC(S1//'IONIZATION OPT 4: IMPOSED RECTANGULAR'//
     >                        ' IONIZATION SOURCE')
        CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
        CALL PRQ(SP//'          TO     :',LENSFI)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                        //' TO RING TARGET FLUX.')
        call prb
c
      elseif (switch(swion).eq.5.0) then
c
        CALL PRC(S1//'IONIZATION OPT 5: ALGORITHMIC RECT/TRI'//
     >                       ' IONIZATION SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET <  1.3 EV)')
        CALL PRC(SP//'     OR L1= 0.0           M (TET >= 1.3 EV)')
        CALL PRC(SP//'     TO L2=L1+2 (M)')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                             ' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                                 ' (TET>10EV)')

        call prb
c
      elseif (switch(swion).eq.6.0) then
c
        CALL PRC(S1//'IONIZATION OPT 6: IMPOSED S**5 GAUSSIAN'//
     >                               ' IONIZATION SOURCE')
        CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA * S**2)')
        CALL PRC(SP//'EXTENDING FROM :  0.0 (M)')
        CALL PRQ(SP//'CUT OFF AT     :',LENSFI)
        CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
        CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                        //' to ring target flux.')
        call prb
c
      elseif (switch(swion).eq.7.0) then
c
        CALL PRC(S1//'IONIZATION OPT 7: ALGORITHMIC'//
     >                       ' RECT/S**5GAUSS IONIZATION'//
     >                    ' SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//
     >                     ' (M) (TET <  1.3EV)')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//
     >                     ' (M) (TET >= 1.3EV)')
        CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                                 ' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                                 ' (TET>10EV)')

        call prb
c
      elseif (switch(swion).eq.9.0) then
c
        CALL PRC(S1//'IONIZATION OPT 9: IMPOSED OFFSET S**5'//
     >                      ' GAUSSIAN IONIZATION SOURCE')
        CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//
     >                             ' EXP(-ALPHA * (S+L)**2)')
        CALL PRC(SP//'EXTENDING FROM :  0.0 (M)')
        CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
        CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
        CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
        CALL PRC(SP//'AND L = WF/2.0')
        CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                        //' TO RING TARGET FLUX.')
        CALL PRB
c
      elseif (switch(swion).eq.10.0) then
c
        CALL PRC(S1//'IONIZ. OPTION 10: ALGORITHMIC RECT/'
     >                    //'OFFSET S**5GAUSS IONIZATION SOURCE')
        CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'
     >                    //' GAUSSIAN SOURCE')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//
     >                     ' (M) (TET <  1.3EV)')
        CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//
     >                     ' (M) (TET >= 1.3EV)')
        CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
        CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                                 ' (TET<10EV)')
        CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                                 ' (TET>10EV)')

        call prb
c
      endif
c
c     Initial Ionization Option
c
c
      if (switch(swion).eq.1.or.switch(swion).eq.2.or.
     >        switch(swion).eq.8) then
c
c        call prb
c        call prr (s1//'PIN ITERATED IONIZATION SOURCE OPTION'//
c     >                ' USED:',
c     >           switch(swion))
c
        call prb
        call prc (s1//'INITIAL SEED PLASMA OPTION:')
        call prb
c
        if (e2dstart.eq.1) then
c
          call prc(s1//'Seed plasma is read in from corresponding'//
     >            ' Edge2D case')
          call prb
c
        elseif (e2dstart.eq.0) then
c
c         regular seed plasma options
c
          if (switch(swioni).eq.0.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 0: EXPONENTIAL DECAY'//
     >                     ' IONIZATION SOURCE')
            CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',
     >                         LENSFI)
            CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',
     >                         LAMS)
c
          elseif (switch(swioni).eq.3.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 3: IMPOSED TRIANGULAR'//
     >                         ' IONIZATION SOURCE')
            CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
            CALL PRQ(SP//'          TO     :',LENSFI)
            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
c
          elseif (switch(swioni).eq.4.0) then
C
            CALL PRC(S1//'INIT IONIZ OPT 4: IMPOSED'//
     >                 ' RECTANGULAR IONIZATION SOURCE')
            CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
            CALL PRQ(SP//'          TO     :',LENSFI)
            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
c
          elseif (switch(swioni).eq.5.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 5: ALGORITHMIC'//
     >                  ' RECT/TRI IONIZATION SOURCE')
            CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET <'
     >                             //'  1.3 EV)')
            CALL PRC(SP//'     OR L1= 0.0           M (TET >'
     >                             //'= 1.3 EV)')
            CALL PRC(SP//'   TO L2=L1+2 (M)')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                               ' (TET>10EV)')
c
          elseif (switch(swioni).eq.6.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 6: IMPOSED S**5'//
     >                    ' GAUSSIAN IONIZATION SOURCE')
            CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA'
     >                          //' * S**2)')
            CALL PRC(SP//'EXTENDING FROM :  0.0 ')
            CALL PRQ(SP//'    CUT OFF AT :',LENSFI)
            CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
            CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
c
          elseif (switch(swioni).eq.7.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 7: ALGORITHMIC'//
     >                   ' RECT/S**5GAUSS IONIZATION')
            CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//
     >                   ' (M) (TET <  1.3EV)')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//
     >                   ' (M) (TET >= 1.3EV)')
            CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                               ' (TET>10EV)')
c
          elseif (switch(swioni).eq.9.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 9: IMPOSED'//
     >                 ' OFFSET S**5 GAUSSIAN')
            CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//
     >                           ' EXP(-ALPHA * (S+L)**2)')
            CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
            CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
            CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
            CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
            CALL PRC(SP//'AND L = WF/2.0')
            CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
c
          elseif (switch(swioni).eq.10.0) then
c
            CALL PRC(S1//'INIT IONIZ OP 10: ALGORITHMIC'//
     >                        ' RECT/OFFSET S**5GAUSS'
     >                  //' IONIZATION SOURCE')
            CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'
     >                           //' GAUSSIAN SOURCE')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//
     >                   ' (M) (TET <  1.3EV)')
            CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//
     >                   ' (M) (TET >= 1.3EV)')
            CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
            CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
            CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                               ' (TET>10EV)')
c
          elseif (switch(swioni).eq.11.0) then
c
            CALL PRC(s1//'INIT IONIZ OP 11: IONIZATION'//
     >                     ' SOURCE DATA READ FROM'//
     >                     ' EDGE2D INPUT FOR CASE.')
C
          elseif (switch(swioni).eq.12.0) then
c
            CALL PRC(S1//'INIT IONIZ OP 12: PIN IS RUN'//
     >                    ' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//
     >                           ' NOT INCLUDE PUFFING.')
c
          elseif (switch(swioni).eq.13.0) then
c
            CALL PRC(S1//'INIT IONIZ OP 13: PIN IS RUN'//
     >                      ' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'IN THE SOL TWICE BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//
     >                           ' PROPER PUFFING.')
c
          elseif (switch(swioni).eq.14.0) then
c
            CALL PRC(S1//'INIT IONIZ OPT 14: PIN IS RUN'//
     >                     ' WITH AN EDGE2D BACKGROUND')
            call prc(sp//'EVERYWHERE TWICE BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//
     >                           ' PROPER PUFFING.')
            CALL PRC(SP//'THIS SHOULD BE USED'//
     >                          ' IN CONJUNCTION'//
     >                          ' WITH CORE OPTION -1.')
c
          elseif (switch(swioni).eq.15.0) then
c
            CALL PRC(s1//'INIT IONIZ OP 15: IONIZATION'//
     >                     ' SOURCE DATA READ FROM'//
     >                     ' EDGE2D INPUT FOR CASE.')
            call prc(sp//'EDGE2D PLASMA ASSIGNED AS OLD'//
     >                   ' FOR ENTIRE BACKGROUND')
c
          elseif (switch(swioni).eq.16.0) then
c
            CALL PRC(S1//'INIT IONIZ OP 16: PIN IS RUN'//
     >             ' WITH A PREVIOUSLY CALCULATED DIVIMP BACKGROUND')
            call prc(sp//'EVERYWHERE - BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//
     >                           ' NOT INCLUDE PUFFING.')
c
          elseif (switch(swioni).eq.17.0) then
c
            CALL PRC(S1//'INIT IONIZ OP 17: PIN IS RUN'//
     >             ' WITH A PREVIOUSLY CALCULATED DIVIMP BACKGROUND')
            call prc(sp//'EVERYWHERE - BEFORE SOL 22 IS INVOKED.')
            CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//
     >                           ' NOT INCLUDE PUFFING.')
            call prc(sp//'EACH SUBSEQUENT ITERATION WILL RE-LOAD THE')
            call prc(sp//'CORE PLASMA SOLUTION (AND THE PRIVATE'//
     >                ' PLASMA SOLUTION')
            call prc(sp//'IF SPECIFIED) FROM THE ORIGINAL INPUT. ONLY')
            call prc(sp//'THE SOL SOLUTION WILL BE ITERATED.')
c
c         End of Seed plasma options
c
          endif
c
c       End of E2dstart
c
        endif
c
c     End of Ionization Option 1,2,8
c
      endif
c
c     Private plasma ionization option - if applicable
c
c
      if (ciopto.eq.1) then
c
c     Only if the TGRAD option is turned ON.
c
        call prb
        CALL PRC(S1//'PRIVATE PLASMA IONIZATION OPTION:')
        call prc(s1//'- NOTE.1: If the main Ionization Option has been'
     >              //' specified as')
        call prc(s1//'  PIN iterated (e.g. Options 1,2,8). Then the'
     >               //' value')
        call prc(s1//'  specified here is ONLY used for the seed'//
     >          ' plasma solution.')
        call prc(s1//'  The PIN iterated ionization option'//
     >                                       ' will be used')
        call prc(s1//'  the rest of the time - except in the case'//
     >         ' of option -2')
        call prc(s1//'  which will use the prescribed private plasma.')
        call prc(s1//'- NOTE.2: The Edge2D seed plasma option takes')
        call prc(s1//'  prcedence over option -2.')
c
c
c slmod begin - new
        if     (switch(swionp).eq.-6.0) then
          CALL PRC(S1//'PP IONIZ OPT  -6: A UNIFORM PLASMA IS ASSIGNED')
          CALL PRC(SP//'TO EACH RING IN THE PRIVATE FLUX ZONE FROM A')
          CALL PRC(SP//'LISTING OF TEMPERATURE AND DENSITY IN THE')
          CALL PRC(SP//'DIVIMP INPUT FILE.  THE TARGET VALUES ARE SET')
          CALL PRC(SP//'FROM THE BULK PLASMA VALUES.')
          fp = 7
          WRITE(fp,1158) 'Ring','Region','Te (eV)','Ti (eV)','n (m-3)',
     .                   'v (m s-1)'
          DO i1 = 1, osmnppv
            WRITE(fp,1159) (INT(osmppv(i1,i2)),i2=1,2),
     .                         (osmppv(i1,i2) ,i2=3,6)
          ENDDO
        elseif (switch(swionp).eq.-5.0) then
          CALL PRC(S1//'PP IONIZ OPT  -6: A UNIFORM PLASMA IS ASSIGNED')
          CALL PRC(SP//'TO EACH RING IN THE PRIVATE FLUX ZONE FROM A')
          CALL PRC(SP//'LISTING OF TEMPERATURE AND DENSITY IN THE')
          CALL PRC(SP//'DIVIMP INPUT FILE.  THE TARGET VALUES ARE NOT')
          CALL PRC(SP//'MODIFIED.')
          fp = 7
          WRITE(fp,1158) 'Ring','Region','Te (eV)','Ti (eV)','n (m-3)',
     .                   'v (m s-1)'
          DO i1 = 1, osmnppv
            WRITE(fp,1159) (INT(osmppv(i1,i2)),i2=1,2),
     .                         (osmppv(i1,i2) ,i2=3,6)
          ENDDO
1158      FORMAT(5X,2A8,2A10,2A12)
1159      FORMAT(5X,2I8,2F10.2,1P,2E12.2,0P)
        elseif (switch(swionp).eq.-4.0) then
          CALL PRC(S1//'PP IONIZ OPT  -3: THOMSON DATA APPLIED TO PP')
          CALL PRC(SP//'AVERAGE VALUE OF TH DATA ON A RING IS ASSIGNED')
          CALL PRC(SP//'TO EVERY CELL ON THE RING.  RINGS WITHOUT ')
          CALL PRC(SP//'DATA ARE INTERPOLATED.  THE TARGET FLUX IS')
          CALL PRC(SP//'ASSIGNED USING THE THOMSON DATA.')
        elseif (switch(swionp).eq.-3.0) then
          CALL PRC(S1//'PP IONIZ OPT  -3: THOMSON DATA APPLIED TO PP')
          CALL PRC(SP//'AVERAGE VALUE OF TH DATA ON A RING IS ASSIGNED')
          CALL PRC(SP//'TO EVERY CELL ON THE RING.  RINGS WITHOUT ')
          CALL PRC(SP//'DATA ARE INTERPOLATED.  THE TARGET FLUX IS')
          CALL PRC(SP//'SPECIFIED IN THE DIVIMP INPUT FILE.')
        elseif (switch(swionp).eq.-2.0) then
c
c        if (switch(swionp).eq.-2) then
c slmod end
c
          CALL PRC(S1//'PP IONIZ OPT  -2: ALL PP CONDITIONS'//
     >                                      ' ARE PRESCRIBED')
          CALL PRC(SP//'BY THE FOLLOWING FACTORS')
C
          CALL PRC(SP//'ALL DISTANCES ARE PROPORTIONS OF SMAX'//
     >          ' FOR THE RING')
C
          CALL PRC(SP//'QUANTITIES ARE LINEARLY INTERPOLATED'
     >                //' BETWEEN POINTS')
          CALL PRC(SP//'AND CONSTANT OUTSIDE RANGE AT END-POINT VALUES')
C
          write(coment,1160)  ctes1,ctef1,ctes2,ctef2
          call prc(coment)
          write(coment,1170)  ctis1,ctif1,ctis2,ctif2
          call prc(coment)
          write(coment,1180)  cnes1,cnef1,cnes2,cnef2
          call prc(coment)
          write(coment,1190)  cvbs1,cvbf1,cvbs2,cvbf2
          call prc(coment)
c
 1160     format(6x,'  (S,Te) : (0.0,Tet) -> (',f6.3,',',f6.3,
     >              '*Tet) -> (',
     >            f6.3,',',f6.3,'*Tet)')
 1170     format(6x,'  (S,Ti) : (0.0,Tit) -> (',f6.3,',',f6.3,
     >             '*Tit) -> (',
     >            f6.3,',',f6.3,'*Tit)')
 1180     format(6x,'  (S,Ne) : (0.0,Net) -> (',f6.3,',',f6.3,
     >                 '*Net) -> (',
     >            f6.3,',',f6.3,'*Net)')
 1190     format(6x,'  (S,Vb) : (0.0,Vbt) -> (',f6.3,',',f6.3,
     >              '*Vbt) -> (',
     >            f6.3,',',f6.3,'*Vbt)')
c
        elseif (switch(swionp).eq.-1) then
c
          CALL PRC(S1//'PP IONIZ OPT  -1: PP IONIZATION OPTIONS'//
     >                             ' ARE THE SAME')
          CALL PRC(SP//'AS FOR THE MAIN SOL.')
c

        elseif (switch(swionp).eq.0.0) then
c
          call prb
          CALL PRC(S1//'PP IONIZ OPT   0: EXPONENTIAL'//
     >                   ' DECAY IONIZATION SOURCE')
          CALL PRB
C
          CALL PRQ (SP//'LENGTH OF IONIZATION SOURCE          ',
     >                         LENSFI)
          CALL PRQ (SP//'DECAY LENGTH OF IONIZATION SOURCE    ',
     >                         LAMS)
c
        elseif (switch(swionp).eq.3.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   3: IMPOSED'//
     >                ' TRIANGULAR IONIZATION SOURCE')
          CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
          CALL PRQ(SP//'          TO     :',LENSFI)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
          CALL PRB
c
        elseif (switch(swionp).eq.4.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   4: IMPOSED'//
     >              ' RECTANGULAR IONIZATION SOURCE')
          CALL PRQ(SP//'EXTENDING FROM   :',LENSST)
          CALL PRQ(SP//'          TO     :',LENSFI)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
          CALL PRB
c
        elseif (switch(swionp).eq.5.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   5: ALGORITHMIC'//
     >                      ' RECT/TRI IONIZATION'
     >                     //' SOURCE')
          CALL PRC(SP//'IF NT > 1.0E19  - TRIANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= (13 - 10*TET) M (TET < '
     >                               //' 1.3 EV)')
          CALL PRC(SP//'     OR L1= 0.0           M (TET >='
     >                               //' 1.3 EV)')
          CALL PRC(SP//'     TO L2=L1+2 (M)')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                               ' (TET>10EV)')

          call prb
c
        elseif (switch(swionp).eq.6.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   6: IMPOSED'//
     >                        ' S**5 GAUSSIAN IONIZATION'
     >                        //' SOURCE')
          CALL PRC(SP//'OF FORM:  A * S**5 * EXP(-ALPHA'
     >                        //' * S**2)')
          CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
          CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
          CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
          CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
          call prb
c
        elseif (switch(swionp).eq.7.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   7: ALGORITHMIC'//
     >                  ' RECT/S**5GAUSS IONIZATION'//
     >                  ' SOURCE - ON')
          CALL PRC(SP//'IF NT > 1.0E19  - S5 GAUSSIAN SOURCE')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = (14-10TET)'//
     >                   ' (M) (TET <  1.3EV)')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = 1.0       '//
     >                   ' (M) (TET >= 1.3EV)')
          CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                                ' (TET>10EV)')

          call prb
c
        elseif (switch(swionp).eq.9.0) then
c
          CALL PRC(S1//'PP IONIZ OPT   9: IMPOSED'//
     >                     ' OFFSET S**5 GAUSSIAN'//
     >                    ' IONIZATION SOURCE')
          CALL PRC(SP//'OF FORM:  A * (S+L)**5 *'//
     >                           ' EXP(-ALPHA * (S+L)**2)')
          CALL PRC(SP//'EXTENDING FROM : 0.0 (M)')
          CALL PRQ(SP//'   CUT OFF AT  :',LENSFI)
          CALL PRC(SP//'WHERE          : ALPHA = 2.5 / WF**2')
          CALL PRQ(SP//'WITH WIDTH FACTOR (WF):',LAMS)
          CALL PRC(SP//'AND L = WF/2.0')
          CALL PRC(SP//'INTEGRAL OF SOURCE NORMALIZED'
     >                      //' TO RING TARGET FLUX.')
          call prb
c
        elseif (switch(swionp).eq.10.0) then
c
          CALL PRC(S1//'PP IONIZ OPT  10: ALGORITHMIC'//
     >                   ' RECT/OFFSET S**5GAUSS'
     >                  //' IONIZATION SOURCE')
          CALL PRC(SP//'IF NT > 1.0E19  - OFFSET S5'
     >                  //' GAUSSIAN SOURCE')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = (28-20TET)'//
     >                   ' (M) (TET <  1.3EV)')
          CALL PRC(SP//'   WITH WIDTH FACTOR WF = 2.0       '//
     >                   ' (M) (TET >= 1.3EV)')
          CALL PRC(SP//'   FROM L1= 0.0 TO L2 = 1/2 RING LENGTH')
          CALL PRC(SP//'IF NT <= 1.0E19 - RECTANGULAR SOURCE')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 13-TET (M)'//
     >                               ' (TET<10EV)')
          CALL PRC(SP//'   FROM L1= 0.0 (M) TO L2= 2 (M)'//
     >                               ' (TET>10EV)')

          call prb
c
        elseif (switch(swionp).eq.11.0) then
c
          CALL PRC(s1//'PP IONIZ OPT  11: IONIZATION'//
     >                   ' SOURCE DATA READ FROM'//
     >                   ' EDGE2D INPUT FOR CASE.')
C
        elseif (switch(swionp).eq.12.0) then
c
          CALL PRC(S1//'PP IONIZ OPT  12: PIN IS RUN'//
     >                  ' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')
          CALL PRC(SP//'THIS GENERATES POWER TERMS BUT WILL'//
     >                         ' NOT INCLUDE PUFFING.')
c
        elseif (switch(swionp).eq.13.0) then
c
          CALL PRC(S1//'PP IONIZ OPT  13: PIN IS RUN'//
     >                  ' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')
          CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//
     >                         ' PROPER PUFFING.')
c
        elseif (switch(swionp).eq.14.0) then
c
          CALL PRC(S1//'PP IONIZ OPT  14: PIN IS RUN'//
     >                  ' WITH AN EDGE2D BACKGROUND')
          call prc(sp//'IN THE SOL ONCE BEFORE SOL 22 IS INVOKED.')
          CALL PRC(SP//'THIS GENERATES POWER TERMS AND'//
     >                         ' PROPER PUFFING.')
          CALL PRC(SP//'THIS SHOULD BE USED'//
     >                        ' IN CONJUNCTION'//
     >                        ' WITH CORE OPTION -1.')
c
        elseif (switch(swionp).eq.15.0) then
c
          CALL PRC(s1//'PP IONIZ OPT  15: IONIZATION'//
     >                   ' SOURCE DATA READ FROM'//
     >                   ' EDGE2D INPUT FOR CASE.')
          call prc(sp//'EDGE2D PLASMA ASSIGNED AS OLD'//
     >                 ' FOR ENTIRE BACKGROUND')
C
        endif
c
c     CIOPTO = 0
c
      elseif (ciopto.eq.0.or.ciopto.eq.2.or.
     >        ciopto.eq.3.or.ciopto.eq.4) then
c
         CALL PRC(S1//'PRIVATE PLASMA IS NOT CALCULATED USING SOL 22.')
         CALL PRC(S1//'PRIVATE PLASMA WILL BE CONSTANT'//
     >                                       ' AT TARGET VALUES')
         CALL PRC(S1//'IF NO OTHER PRIVATE PLASMA OPTION HAS BEEN'//
     >                   ' SPECIFIED.')

c
c     ENDIF for CIOPTO
c
      endif
c
      call prb
c
c     Convection Term number 1
c
      if (switch(swcond).eq.0.0) then
         CALL PRC(S1//'5/2 NV*kT OPT 0 : THIS TERM IF OFF')
      else
         CALL PRC(S1//'5/2 NV*kT OPT 1 : THIS TERM IF ON')
      endif
c
      call prb
c
c     Convection Term number 2
c
      if (switch(swconv).eq.0.0) then
         CALL PRC(S1//'1/2 mV^3N OPT 0 : THIS TERM IS OFF')
      else
         CALL PRC(S1//'1/2 mV^3N OPT 1 : THIS TERM IS ON')
      endif
c
c     Prad term
c
      call prb
c
      if (switch(swprad).eq.0.0) then
         CALL PRC(S1//'PRAD OPTION   0 : PRAD TERM TERM  IS OFF')
      elseif (switch(swprad).eq.1.0) then
         CALL PRC(S1//'PRAD OPTION   1 : EXPONENTIAL'//
     >                       ' DECAY RADIATION SOURCE')
         CALL PRQ(SP//'LENGTH OF RADIATION SOURCE         ', LENR)
         CALL PRQ(SP//'DECAY LENGTH OF RADIATION  SOURCE  ', LAMR)
         CALL PRQ(SP//'SOURCE STRENGTH FRACTION (FRR)     ', FRR)
      elseif (switch(swprad).eq.2.0) then
         CALL PRC(S1//'PRAD OPTION   2 : ANALYTIC'//
     >                     ' (GARCHING) RADIATION SOURCE')
         CALL PRQ(SP//'ALPHA = NIMP/NE      ',ALFIMP)
         CALL PRQ(SP//'TBASE FACTOR         ',TALIMP)
         CALL PRQ(SP//'FIRST EXPONENT       ',EX1IMP)
         CALL PRQ(SP//'SECOND EXPONENT      ',EX2IMP)
      elseif (switch(swprad).eq.3.0) then
         CALL PRC(S1//'PRAD OPTION   3 : PROPORTIONAL TO PINQE')
         CALL PRC(SP//'IMPURITY RADIATION IS A GIVEN FACTOR')
         CALL PRC(SP//'TIMES THE CALCULATED PINQE.')
         CALL PRQ(SP//'MULTIPLICATION FACTOR = ',radsrc_mult)
      elseif (switch(swprad).eq.4.0) then
         CALL PRC(S1//'PRAD OPTION   4 : EXTERNAL RADIATION SOURCE')
         CALL PRC(SP//'IMPURITY RADIATION IS READ IN FROM A FILE')
         CALL PRC(SP//'USUALLY RESULTING FROM A PREVIOUS DIVIMP')
         CALL PRC(SP//'CASE.')
         CALL PRQ(SP//'MULTIPLICATION FACTOR = ',radsrc_mult)
      elseif (switch(swprad).eq.5.0) then 
         CALL PRC(S1//'PRAD OPTION   5 : EXTERNAL RADIATION SOURCE')
         call prc(s1//' **** WARNING ****')
         call prc(s1//' IN DEVELOPMENT - NOT YET IMPLEMENTED')
         call prc(s1//' **** WARNING ****')
         CALL PRC(SP//'TOTAL RADIATION TO BE APPLIED IN SOL22 IS')
         CALL PRC(SP//'SPECIFIED IN THE INPUT ON A PER REGION BASIS')
         CALL PRC(SP//'BASED ON EXPERIMENTAL BOLOMETRY DATA')
         CALL PRC(SP//'THIS IS THEN BE DISTRIBUTED TO EACH RING AND')
         CALL PRC(SP//'CELL USING VARIOUS MECHANISMS')
c         call prq(sp//'REGION N ', prad_for_region)
      elseif (switch(swprad).eq.6.0) then
         CALL PRC(S1//'PRAD OPTION   6 : RECTANGULAR'//
     >                       ' RADIATION SOURCE')
         CALL PRQ(SP//'END LENGTH OF RADIATION SOURCE     ', LENR)
         CALL PRQ(SP//'START LENGTH OF RADIATION  SOURCE  ', LAMR)
         CALL PRQ(SP//'SOURCE STRENGTH FRACTION (FRR)     ', FRR)
      endif
c
      call prb
c
c     Phelpi term
c
      if (switch(swphelp).eq.0.0) then
         CALL PRC(S1//'PHELPI OPTION 0 : PHELPI TERM  IS OFF')
      elseif (switch(swphelp).eq.1.0) then
         CALL PRC(S1//'PHELPI OPTION 1 : INTERNAL'//
     >                                ' PHELPI TERM  IS ON')
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swphelp).eq.2.0) then
         CALL PRC(S1//'PHELPI OPTION 2 : PINQE TERM  IS ON')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swphelp).eq.3.0) then
         CALL PRC(S1//'PHELPI OPTION 3 : PINQE TERM  IS ON')
         CALL PRC(SP//'INTERNAL PHELPI TERM IS OFF FOR SEED PLASMA')
      endif
c
      if (tcutqe.gt.0.0) then
         CALL PRB
         CALL PRQ(S1//'QE TEMPERATURE CUTOFF (EV) = ',TCUTQE)
         CALL PRC(S1//'- ALL CONTRIBUTIONS TO PHELPI TERM')
         call prc(s1//'  FROM REGIONS WITH TE<TCUT HAVE BEEN SET')
         call prc(s1//'  TO ZERO.')
      endif
c
      call prb
c
c     Pei Term
c
      if (switch(swpei).eq.0.0) then
         CALL PRC(S1//'PEI OPTION    0 : PEI TERM  IS OFF')
      elseif (switch(swpei).eq.1.0) then
         CALL PRC(S1//'PEI OPTION    1 : PEI TERM  IS ON')
         CALL PRQ(SP//'PEI CORRECTION FACTOR =',PEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpei).eq.3.0) then
         CALL PRC(S1//'PEI OPTION    3 : PEI TERM  IS'//
     >                ' CALCULATED')
         CALL PRC(SP//'BUT NOT USED.')
      endif
c
      call prb
c
c     PCX term
c
      if (switch(swpcx).eq.0.0) then
         CALL PRC(S1//'PCX OPTION    0 : PCX TERM  IS OFF')
      elseif (switch(swpcx).eq.1.0) then
         CALL PRC(S1//'PCX OPTION    1 : INTERNAL PCX TERM  IS ON')
         CALL PRQ(SP//'CX POWER COEFF (CEICF) =',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpcx).eq.2.0) then
         CALL PRC(S1//'PCX OPTION    2 : PINQI TERM IS ON')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRQ(SP//'CEICF FOR SEED PLASMA',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      elseif (switch(swpcx).eq.3.0) then
         CALL PRC(S1//'PCX OPTION    3 : PINQI TERM IS ON')
         CALL PRC(SP//'INTERNAL PCX TERM IS OFF FOR SEED PLASMA')
      elseif (switch(swpcx).eq.4.0) then
         CALL PRC(S1//'PCX OPTION    4 : DIVIMP QI TERM  IS ON')
         CALL PRC(SP//'INTERNAL TERM IS OFF FOR SEED PLASMA')

         call prb
         CALL PRC(S2//'THE FOUR COMPONENTS OF DIVIMP QI'//
     >            ' ARE CALCULATED USING')
         CALL PRC(S2//'THE FOLLOWING SUB-OPTIONS')
         call prb
c
c        Add switch print-outs for swqidatiz,swqidmliz,swqidcx,swqidrec
c
         CALL PRC(S2//'PINQID SUB-OPTION LISTING:')
c
         call prb
c
         if (switch(swqidatiz).eq.0.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 0 : OFF')
         elseif (switch(swqidatiz).eq.1.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 1 : ON')
            CALL PRC(S2//' - 3/2 K TATOM * (NH/(NH+NH2)) * SIONIZ')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTATIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidatiz).eq.2.0) then
            CALL PRC(S2//'ATOMIC IONIZATION OPT 2 : ON')
            CALL PRC(S2//' - 3/2 K TATOM * NE * NH * <SV>ATIZ (ADAS)')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTATIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         endif
c
         if (switch(swqidmliz).eq.0.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 0 : OFF')
         elseif (switch(swqidmliz).eq.1.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 1 : ON')
            CALL PRC(S2//' - 3.0 * (NH2/(NH+NH2)) * SIONIZ')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTMLIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidmliz).eq.2.0) then
            CALL PRC(S2//'MOLECULAR IONIZATION OPT 2 : ON')
            CALL PRC(S2//' - 3.0 * NE * NH2 * <SV>ATIZ (ADAS)')
            CALL PRC(S2//' - USES ATOMIC IONIZATION AS AN')
            CALL PRC(S2//'   APPROXIMATION TO MOLECULAR IONZIATION.')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTMLIZ)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         endif
c
         if (switch(swqidrec).eq.0.0) then
            CALL PRC(S2//'RECOMBINATION OPTION 0 : OFF')
         elseif (switch(swqidrec).eq.1.0) then
            CALL PRC(S2//'RECOMBINATION OPTION 1 : ON')
            CALL PRC(S2//' - -3/2 * K TI * NI * NE * <SV>REC (ADAS)')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTREC)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         endif
c
         if (switch(swqidcx).eq.0.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 0 : OFF')
         elseif (switch(swqidcx).eq.1.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 1 : ON')
            CALL PRC(S2//' - 3/2 K TREF * RCXMULT(TI) * CEICF *'//
     >                       ' SIONIZ')
            CALL PRQ(S2//' - TREF (EV) = ', TREFCX)
            CALL PRQ(S2//' - CEICF     = ',CEICF)
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.2.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 2 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH * <SV>CX')
            CALL PRC(S2//' - DEFINE EAV = 3/2 K (TATOM+ TI) / 2.0')
            CALL PRC(S2//' - <SV>CX = 10E-13           '//
     >                        '  FOR EAV  > 1000 EV')
            CALL PRC(S2//' - <SV>CX = 10E-14 EAV**(1/3)'//
     >                        '  FOR EAV =< 1000 EV')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.3.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 3 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH *'//
     >                       ' <SV>CX (ADAS)')
            CALL PRC(S2//' - ADAS CX DATA IS CURRENTLY UNRELIABLE')
            CALL PRC(S2//'   CHECK BEFORE USING THIS OPTION.')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
         elseif (switch(swqidcx).eq.4.0) then
            CALL PRC(S2//'CHARGE EXCHANGE OPTION 4 : ON')
            CALL PRC(S2//' - 3/2 * (TATOM - TI) * NE * NH * <SV>CX')
            CALL PRC(S2//' - DEFINE EAV = 3/2 K (TATOM+ TI) / 2.0')
            CALL PRC(S2//' - <SV>CX = 10E-13           '//
     >                        '  FOR EAV  > 1000 EV')
            CALL PRC(S2//' - <SV>CX = 10E-14 EAV**(1/3)'//
     >                        '  FOR EAV =< 1000 EV')
            CALL PRQ(S2//' - CUTOFF TEMPERATURE (EV) = ',TCUTCX)
            CALL PRC(S2//' - CONTRIBUTIONS = 0 IF TI < TCUT')
            CALL PRC(S2//' - TERM IS LIMITED TO ION COOLING ONLY')
         endif
c
      elseif (switch(swpcx).eq.5.0) then
         CALL PRC(S1//'PCX OPTION    5 : PINQI TERM IS ON')
         call prc(sp//'PINQI TERM IS CLIPPED - IT IS NOT ALLOWED')
         call prc(sp//'TO HEAT THE PLASMA.')
         CALL PRC(SP//'INTERNAL TERM IS ON FOR SEED PLASMA')
         CALL PRQ(SP//'CEICF FOR SEED PLASMA',CEICF)
         CALL PRC(SP//'SEE MANUAL FOR DETAILS.')
      endif
c
      if (switch(swpcx).ne.4.0.and.tcutcx.gt.0.0) then
         call prb
         CALL PRQ(S1//'QI TEMPERATURE CUTOFF (EV) = ',TCUTCX)
         CALL PRC(S1//'- ALL CONTRIBUTIONS FROM REGIONS WITH TI<TCUT')
         CALL PRC(S1//'  HAVE BEEN SET TO ZERO.')
      endif
c
c     Private plasma electron TARGET power LOSS compensation term
c
      call prb
c
      if (switch(swppelec).eq.0.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 0: PP ELECTRON POWER'//
     >                ' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppelec).eq.1.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 1: PP ELECTRON POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ELECTRON POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppelec).eq.2.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 2: PP ELECTRON POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ELECTRON POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppelec).eq.3.0) then
         CALL PRC(S1//'PP ELEC POW OPTION 3: PP ELECTRON POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IT IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      ENDIF
c
c     Private plasma ION TARGET power LOSS compensation term
c
      call prb
c
      if (switch(swppion).eq.0.0) then
         CALL PRC(S1//'PP ION POW OPTION 0 : PP ION POWER'//
     >                ' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppion).eq.1.0) then
         CALL PRC(S1//'PP ION POW OPTION 1 : PP ION POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppion).eq.2.0) then
         CALL PRC(S1//'PP ION POW OPTION 2 : PP ION POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS REMOVED FROM MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS POWER LOSS')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppion).eq.3.0) then
         CALL PRC(S1//'PP ION POW OPTION 3 : PP ION POWER'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'ION POWER LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IT IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      ENDIF
c
c     Private plasma PRESSURE LOSS compensation term
c
      call prb
c
      if (switch(swppress).eq.0.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 0 : PP PRESSURE'//
     >                ' LOSS COMPENSATION TERM IS OFF')
      elseif (switch(swppress).eq.1.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 1 : PP PRESSURE'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS ADDED TO MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS PRESSURE')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      elseif (switch(swppress).eq.2.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 2 : PP PRESSURE'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO EACH ELEMENT OF PP')
         call prc(sp//'TARGET IS ADDED TO MAIN SOL RINGS IN')
         call prc(sp//'THE CORRESPONDING POSITION RELATIVE TO THE')
         CALL PRC(SP//'SEPARATRIX. THIS PRESSURE')
         call prq(sp//'IS DISTRTIBUTED EVENLY OVER SMAX *',pp_pow_dist)
         call prc(sp//'FROM EACH TARGET ON THE MAIN SOL RINGS.')
      elseif (switch(swppress).eq.3.0) then
         CALL PRC(S1//'PP PRESSURE OPTION 3 : PP PRESSURE'//
     >                ' LOSS COMPENSATION TERM IS ON')
         call prc(sp//'PRESSURE LOST TO LOCAL PFZ TARGET')
         call prc(sp//'IS DISTRIBUTED TO MAIN SOL RINGS USING')
         call prc(sp//'ONE OF 5 POSSIBLE RING DISTRIBUTION OPTIONS')
         call prc(sp//'IS DISTRTIBUTED EVENLY TO THE XPOINT ON')
         call prc(sp//'THE MAIN SOL RINGS')
      ENDIF
c
c     Viscosity term
c
      if (switch(swvisc1).ne.0.0) then
         call prb
         CALL PRC(S1//'DENSITY VISCOSITY CORRECTION TERM'//
     >            '  IS NOT IMPLEMENTED')
      endif
c
c     Momentum Loss term
c
      call prb
c
      if (switch(swnmom).eq.0.0) then
         CALL PRC(S1//'MOMENTUM OPT  0 : MOMENTUM LOSS TERM IS OFF')
      elseif (switch(swnmom).eq.1.0) then
         CALL PRC(S1//'MOMENTUM OPT  1 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = SMOM0    S < L * SMAX        ')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(L * SMAX) * (1/FFRIC -1)')
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
            CALL PRQ(SP//'L    = ',LENMOM)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            CALL PRQ(SP//'DEFAULT L    = ',LENMOM)
            call prc(sp//'OVER-RIDDEN BY THE'//
     >                   ' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'         FFRIC  LENMOM    FFRIC LENMOM')
            call prc(sp//'RING    '//OUTER//'      '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,4(1x,f8.3))') 
     >            int(extffric(in,1)),extffric(in,2),extffric(in,3),
     >            extffric(in,4),extffric(in,5)
               call prc(sp//coment)
            end do 
            call prc(sp//'A VALUE <= 0.0 = DEFAULT')
         end if

      elseif (switch(swnmom).eq.2.0) then
         CALL PRC(S1//'MOMENTUM OPT  2 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = SMOM0 * EXP ( -S / (LM * SMAX))  FOR'
     >            //' S < L * SMAX')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(LM * SMAX)*(1/FFRIC -1) / EXPF')
         CALL PRC(SP//'EXPF = (1 - EXP(-(L * SMAX)/(LM * SMAX)))')
         CALL PRQ(SP//'LM   = ',LAMMOM)
         CALL PRQ(SP//'L    = ',LENMOM)
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            call prc(sp//'OVER-RIDDEN BY THE'//
     >                   ' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'RING    '//OUTER//'   '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,2(1x,f8.3))') 
     >            int(extffric(in,1)),extffric(in,2),extffric(in,3)
               call prc(sp//coment)
            end do 
         end if
      elseif (switch(swnmom).eq.3.0) then
         CALL PRC(S1//'MOMENTUM OPT  3 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = FACT * SIZ(S) / INT (0 TO L) '//
     >            '[ SIZ(S'') DS'' ] ')
         CALL PRC(SP//'FACT = PT *  (1/FFRIC -1) ')
         CALL PRC(SP//'L    = SMAX / 2.0')
         if (n_extffric.eq.0) then 
            CALL PRQ(SP//'FFRIC= ',FFRIC)
         else
            CALL PRQ(SP//'DEFAULT FFRIC= ',FFRIC)
            call prc(sp//'OVER-RIDDEN BY THE'//
     >                   ' FOLLOWING VALUES ON SPECIFIED RINGS:')
            call prc(sp//'RING    '//OUTER//'   '//INNER)
            do in = 1, n_extffric
               write(coment,'(i4,2x,2(1x,f8.3))') 
     >            int(extffric(in,1)),extffric(in,2),extffric(in,3)
               call prc(sp//coment)
            end do 
         end if
      elseif (switch(swnmom).eq.4.0) then
         CALL PRC(S1//'MOMENTUM OPT  4 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM = -MB * VB(S) * RCXMULT * RCXIZ'//
     >                      ' *  SIZ(S) ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)
         CALL PRC(SP//'RCXMULT = A * EXP(-B * T) ')
         CALL PRC(SP//'B = 6.907/(TCX - 1) ')
         CALL PRC(SP//'A = 1000 * EXP (B)  ')
         CALL PRQ(SP//'TCX = ',TCXMOM)
         CALL PRC(SP//'WHERE: 1.0 < RCXMULT < 1500.0 IS FORCED')
      elseif (switch(swnmom).eq.5.0) then
         CALL PRC(S1//'MOMENTUM OPT  5 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'SMOM IS READ DIRECTLY FROM NIMBUS')
      elseif (switch(swnmom).eq.6.0) then
         CALL PRC(S1//'MOMENTUM OPT  6 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'MOMENTUM LOSS TERM IS OFF')
      elseif (switch(swnmom).eq.7.0) then
         CALL PRC(S1//'MOMENTUM OPT  7 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'INITIAL NEUTRAL BG MOMENTUM LOSS TERM IS:')
         CALL PRC(SP//'SMOM = SMOM0    S < L * SMAX        ')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(L * SMAX) * (1/FFRIC -1)')
         CALL PRQ(SP//'FFRIC= ',FFRIC)
         CALL PRQ(SP//'L    = ',LENMOM)
      elseif (switch(swnmom).eq.8.0) then
         CALL PRC(S1//'MOMENTUM OPT  8 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'NEUTRAL BG MOMENTUN LOSS TERM IS TAKEN')
         CALL PRC(SP//'FROM PIN RESULTS EXCEPT FOR THE FIRST')
         CALL PRC(SP//'ITERATION, IN WHICH CASE: ')
         CALL PRC(SP//'INITIAL NEUTRAL BG MOMENTUM LOSS TERM IS:')
         CALL PRC(SP//'SMOM = SMOM0 * EXP ( -S / (LM * SMAX))  FOR'
     >            //' S < L * SMAX')
         CALL PRC(SP//'     = 0        S > L * SMAX        ')
         CALL PRC(SP//'SMOM0= PT/(LM * SMAX)*(1/FFRIC -1) / EXPF')
         CALL PRC(SP//'EXPF = (1 - EXP(-(L * SMAX)/(LM * SMAX)))')
         CALL PRQ(SP//'FFRIC= ',FFRIC)
         CALL PRQ(SP//'LM   = ',LAMMOM)
         CALL PRQ(SP//'L    = ',LENMOM)
      elseif (switch(swnmom).eq.9.0) then
         CALL PRC(S1//'MOMENTUM OPT  9 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'BASED ON CX EXCHANGE MOMENTUM CROSS-'//
     >                     'SECTIONS')
         CALL PRC(SP//'TAKEN FROM THE EDGE2D/NIMBUS'//
     >                ' IMPLEMENTATION.')
         CALL PRC(SP//'REF: WOJCIECH')
         CALL PRC(SP//'MOMENTUM LOSS MULTIPLIED BY RCXIZ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)
      elseif (switch(swnmom).eq.10.0) then
         CALL PRC(S1//'MOMENTUM OPT 10 : MOMENTUM LOSS TERM IS ON')
         CALL PRC(SP//'BASED ON CX EXCHANGE MOMENTUM CROSS-'//
     >             'SECTIONS TAKEN FROM')
         CALL PRC(SP//'EDGE2D/NIMBUS IMPLEMENTATION.')
         CALL PRC(SP//'REF: WOJCIECH')
         CALL PRC(SP//'     SOME SORT OF ALTERNATIVE TO'//
     >                  'OPTION 9')
         CALL PRC(SP//'     FACTORING IN H0 VELOCITY'//
     >                  ' DISTRIBUTIONS')
         CALL PRC(SP//'MOMENTUM LOSS MULTIPLIED BY RCXIZ')
         CALL PRQ(SP//'RCXIZ= ',RCXMOM)
      endif
c
c     Mulitplication factor for Momentum loss option
c
      if (switch(swnmom).ne.0.0) then
         CALL PRQ(SP//'MOMENTUM SOURCE MULT= ',smom_mult)
      endif
c
c     Iterative Mach number
c
      call prb
c
      if (switch(swmach).eq.0.0) then
         CALL PRC (S1//'MACH OPTION   0 : ITERATIVE MACH'//
     >                          ' SOLUTION IS OFF')
         CALL PRQ (SP//'IMPOSED MACH NUMBER AT THE TARGET ', M0)

      elseif (switch(swmach).eq.1.0) then
         CALL PRC (S1//'MACH OPTION   1 : ITERATIVE MACH'//
     >                          ' SOLUTION IS ON')
         CALL PRC (SP//'TARGET DENSITY IS FIXED WITH MACH NUMBER')
         CALL PRC (SP//'AS A RESULT, THE TARGET FLUX WILL CHANGE.')
      elseif (switch(swmach).eq.2.0) then
         CALL PRC (S1//'MACH OPTION   2 : ITERATIVE MACH'//
     >                          ' SOLUTION IS ON')
         CALL PRC (SP//'TARGET DENSITY VARIES WITH MACH NUMBER TO'//
     >                    ' MAINTAIN A CONSTANT TARGET FLUX.')
      elseif (switch(swmach).eq.3.0) then
         CALL PRC (S1//'MACH OPTION   3 : FIXED MACH NUMBER'//
     >                          ' SOLUTION IS ON')
         CALL PRC (SP//'MACH NUMBER IS SET BASED ON INPUT TARGET'//
     >                  ' FLOW VELOCITY RELATIVE TO SOUND SPEED.')
      endif
c
c     Edge2D compatibility
c
      call prb
c
      if (e2dstart.eq.1) then
         CALL PRC (S1//'WARNING: ')
         CALL PRB
         CALL PRC (S1//'EDGE2D STARTER PLASMA OPTION IS ON')
         CALL PRC (S1//'THE INITIAL PLASMA SOLUTION FOR SOL OPTION 22')
         CALL PRC (S1//'HAS BEEN READ FROM AN EDGE2D INPUT FILE.')
         CALL PRC (S1//'NONE OF THE SPECIFIED SWITCHES APPLY TO THE')
         CALL PRC (S1//'STARTER PLASMA.')
         CALL PRC (S1//'NOTE: THIS IS SET BY USING -OPT FOR THE EDGE2D')
         CALL PRC (S1//'         COMPATIBILITY OPTION.')
         CALL PRB
      endif
c
      if (switch(swe2d).eq.0.0) then
         CALL PRC(S1//'EDGE2D BG OPT 0 : COMPATIBILITY IS OFF')
         CALL PRC(SP//'SOLVER RUNS FROM THE ACTUAL TARGET')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES WILL ALSO BE USED IN')
         CALL PRC(SP//'THE SOLVER')
      elseif (switch(swe2d).eq.1.0) then
         CALL PRC(S1//'EDGE2D BG OPT 1 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'EDGE2D PRESSURE IS MATCHED AT FIRST')
         CALL PRC(SP//'CELL CENTRE')
      elseif (switch(swe2d).eq.2.0) then
         CALL PRC(S1//'EDGE2D BG OPT 2 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS READ FROM THE EDGE2D OUTPUT')
         CALL PRC(SP//'AT THE FIRST CELL CENTRE')
      elseif (switch(swe2d).eq.3.0) then
         CALL PRC(S1//'EDGE2D BG OPT 3 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.4.0) then
         CALL PRC(S1//'EDGE2D BG OPT 4 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS READ FROM THE EDGE2D OUTPUT')
         CALL PRC(SP//'AT THE FIRST CELL CENTRE')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.5.0) then
         CALL PRC(S1//'EDGE2D BG OPT 5 : COMPATIBILITY IS OFF')
         call prc(sp//'FLUXES ARE CALCULATED FROM THE TARGET')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.6.0) then
         CALL PRC(S1//'EDGE2D BG OPT 6 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'EDGE2D PRESSURE IS MATCHED AT FIRST')
         CALL PRC(SP//'CELL CENTRE')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.7.0) then
         CALL PRC(S1//'EDGE2D BG OPT 7 : COMPATIBILITY IS ON')
         CALL PRC(SP//'EDGE2D DATA IS READ FOR FIRST CELL')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'** NOTE: SOL 22 IS NOT RUN. THE EDGE2D')
         CALL PRC(SP//'         SOLUTION IS USED FOR THE BACKGROUND')
         CALL PRC(SP//'         PLASMA. THIS OPTION IS USED ONLY TO')
         CALL PRC(SP//'         PRODUCE DETAILED FLUX ANALYSES FOR')
         CALL PRC(SP//'         DEBUGGING.')
      elseif (switch(swe2d).eq.8.0) then
         CALL PRC(S1//'EDGE2D BG OPT 8 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF THE FIRST CELL')
         CALL PRC(SP//'EDGE2D DATA IS REQUIRED FOR THE ENTIRE RING.')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE FIRST CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'THESE FLUXES ARE EXTRACTED FROM THE EDGE2D')
         call prc(sp//'DOWN FLUX LISTING.')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES AT THE CELL CENTRE')
         CALL PRC(SP//'WILL ALSO BE USED IN THE SOLVER')
      elseif (switch(swe2d).eq.9.0) then
         CALL PRC(S1//'EDGE2D BG OPT 9 : COMPATIBILITY IS ON')
         CALL PRC(SP//'SOLVER RUNS FROM THE MIDDLE OF A SPECIFIED CELL')
         CALL PRI(SP//'STARTING KNOT INDEX = ',IKE2D)
         CALL PRC(SP//'EDGE2D DATA IS REQUIRED FOR THE ENTIRE RING.')
         CALL PRC(SP//'VELOCITY IS SET BY EXTRACTING THE EDGE2D')
         CALL PRC(SP//'FLUXES ENTERING AND LEAVING THE CELL')
         CALL PRC(SP//'AND TAKING THEIR AVERAGE AND DIVIDING BY')
         CALL PRC(SP//'EDGE2D CELL DENSITY.')
         call prc(sp//'THESE FLUXES ARE EXTRACTED FROM THE EDGE2D')
         call prc(sp//'DOWN FLUX LISTING.')
         call prc(sp//'IF E2D TARGET OPTION 5 IS SELECTED THEN THE')
         CALL PRC(SP//'EDGE2D DOWN POWER FLUXES AT THE CELL CENTRE')
         CALL PRC(SP//'WILL ALSO BE USED IN THE SOLVER')
      endif
c
c     Power Distribution Option
c
c
      call prb
c
      if (switch(swpow).eq.0.0) then
         CALL PRC(S1//'POWER DIST OPT 0: DISTRIBUTED INFLUX IS OFF')
      elseif (switch(swpow).eq.1.0) then
         CALL PRC(S1//'POWER DIST OPT 1: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
      elseif (switch(swpow).eq.2.0) then
         CALL PRC(S1//'POWER DIST OPT 2: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM THE X-POINT REGION TO THE MIDPOINT')
      elseif (switch(swpow).eq.3.0) then
         CALL PRC(S1//'POWER DIST OPT 3: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT.')
         CALL PRC(SP//'CORRECTED FOR MAJOR RADIUS EFFECTS. THIS')
         CALL PRC(SP//'OPTION USEFUL ONLY WITH MAJOR RADIUS'//
     >                 ' OPTION 4')
      elseif (switch(swpow).eq.4.0) then
         CALL PRC(S1//'POWER DIST OPT 4: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
      elseif (switch(swpow).eq.5.0) then
         CALL PRC(S1//'POWER DIST OPT 5: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                      //' DISTRIBUTED.')
      elseif (switch(swpow).eq.6.0) then
         CALL PRC(S1//'POWER DIST OPT 6: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                  //' DISTRIBUTED.')
      elseif (switch(swpow).eq.7.0) then
         CALL PRC(S1//'POWER DIST OPT 7: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpow).eq.8.0) then
         CALL PRC(S1//'POWER DIST OPT 8: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                 //' DISTRIBUTED.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpow).eq.9.0) then
         CALL PRC(S1//'POWER DIST OPT 9: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.10.0) then
         CALL PRC(S1//'POWER DIST OP 10: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                 //' DISTRIBUTED.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.11.0) then
         CALL PRC(S1//'POWER DIST OP 11: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE RING FROM :')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRC(SP//'(FROM EACH TARGET) ')
      endif
c
      call prb
c
c     Private Plasma Power Distribution
c
      if (switch(swpowp).eq.0.0) then
         CALL PRC(S1//'PP POWER OPT  0 : DISTRIBUTED INFLUX IS OFF')
      elseif (switch(swpowp).eq.1.0) then
         CALL PRC(S1//'PP POWER OPT  1 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
      elseif (switch(swpowp).eq.2.0) then
         CALL PRC(S1//'PP POWER OPT  2 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM THE X-POINT REGION TO THE MIDPOINT')
      elseif (switch(swpowp).eq.3.0) then
         CALL PRC(S1//'PP POWER OPT  3 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT.')
         CALL PRC(SP//'CORRECTED FOR MAJOR RADIUS EFFECTS. THIS')
         CALL PRC(SP//'OPTION USEFUL ONLY WITH MAJOR RADIUS'//
     >                 ' OPTION 4')
      elseif (switch(swpowp).eq.4.0) then
         CALL PRC(S1//'PP POWER OPT  4 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
      elseif (switch(swpowp).eq.5.0) then
         CALL PRC(S1//'PP POWER OPT  5 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'ALONG THE HALF RING TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                      //' DISTRIBUTED.')
      elseif (switch(swpowp).eq.6.0) then
         CALL PRC(S1//'PP POWER OPT  6 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE ENTIRE RING.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                  //' DISTRIBUTED.')
      elseif (switch(swpowp).eq.7.0) then
         CALL PRC(S1//'PP POWER OPT  7 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpowp).eq.8.0) then
         CALL PRC(S1//'PP POWER OPT  8 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'FROM A SPECIFIED POINT TO THE MIDPOINT')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                 //' DISTRIBUTED.')
         CALL PRQ(SP//'SPECIFIED POINT = SMAX * ',SPOWBEG)
      elseif (switch(swpowp).eq.9.0) then
         CALL PRC(S1//'PP POWER OPT  9 : DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpowp).eq.10.0) then
         CALL PRC(S1//'PP POWER OPT  10: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS DISTRIBITED LINEARLY')
         CALL PRC(SP//'OVER A SPECIFIED RANGE.')
         CALL PRC(SP//'IF PIN BASED VOLUME POWER SOURCES ARE ON')
         CALL PRC(SP//'THEN THE INTEGRAL OF THESE IS ALSO'
     >                 //' DISTRIBUTED.')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRQ(SP//'ENDING   POINT = SMAX * ',SPOWLEN)
      elseif (switch(swpow).eq.11.0) then
         CALL PRC(S1//'PP POWER OPT  11: DISTRIBUTED INFLUX IS ON')
         CALL PRC(SP//'TARGET POWER FLUX IS SUMMED FOR BOTH')
         CALL PRC(SP//'TARGETS AND DISTRIBUTED LINEARLY')
         CALL PRC(SP//'OVER THE RING FROM :')
         CALL PRQ(SP//'STARTING POINT = SMAX * ',SPOWBEG)
         CALL PRC(SP//'(FROM EACH TARGET) ')
      endif
c
      call prb
c
c     Perpendicular flux option
c
      if (switch(swgperp).eq.0.0) then
         CALL PRC(S1//'PERP FLUX OPT  0: TERM IS OFF')
      elseif (switch(swgperp).eq.1.0) then
         CALL PRC(S1//'PERP FLUX OPT  1: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//
     >            ' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
      elseif (switch(swgperp).eq.2.0) then
         CALL PRC(S1//'PERP FLUX OPT  2: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
      elseif (switch(swgperp).eq.3.0) then
         CALL PRC(S1//'PERP FLUX OPT  3: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperp).eq.4.0) then
         CALL PRC(S1//'PERP FLUX OPT  4: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE HALF FIELD LINE ->0')
         CALL PRC(SP//'AT THE MIDPOINT -  ')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperp).eq.5.0) then
         CALL PRC(S1//'PERP FLUX OPT  5: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//
     >            ' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX OVER THE 1/2 FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE'
     >             //' FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperp).eq.6.0) then
         CALL PRC(S1//'PERP FLUX OPT  6: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE ENTIRE FILED LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >           ' PERPENDICULAR FLUX OVER THE WHOLE FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE'
     >             //' FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperp).eq.7.0) then
         CALL PRC(S1//'PERP FLUX OPT  7: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX PROPORTIONAL')
         CALL PRC(SP//'TO THE SECOND GRADIENT IN THE DENSITY')
         CALL PRC(SP//'FROM AN EDGE2D SOLUTION OR A PREVIOUS')
         CALL PRC(SP//'SOL22 SOLUTION. THIS OPTION IS CHANGED')
         CALL PRC(SP//'TO A UNIFORM DISTRIBUTION FOR ANY RINGS')
         CALL PRC(SP//'WHERE THE RATIO OF EITHER THE POSITIVE')
         CALL PRC(SP//'OR NEGATIVE CONTRIBUTION TO THE INTEGRATED')
         CALL PRC(SP//'VALUE FOR TEH RING EXCEEDS 5.0' )
      elseif (switch(swgperp).eq.8.0) then
         CALL PRC(S1//'PERP FLUX OPT  8: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX CALCULATED')
         CALL PRC(SP//'FROM THE SECOND GRADIENT IN THE DENSITY')
         call PRC(SP//'AND AN IMPOSED VALUE FOR THE DIFFUSION')
         CALL PRR(SP//'RATE:  DPERP [M2/S]  = ',CDPERP)
         CALL PRC(SP//'ANY REMAINING EXECESS OR SHORTFALL OF')
         CALL PRC(SP//'FLUX IS COMPENSATED FOR BY THE APPLICATION')
         CALL PRC(SP//'OF AN ADDITIONAL UNIFORM SOURCE OVER THE')
         CALL PRC(SP//'ENTIRE FLUX TUBE.')
      endif
c
      call prb
c
      if (switch(swgperpp).eq.0.0) then
         CALL PRC(S1//'PP PERPF  OPT  0: TERM IS OFF')
      elseif (switch(swgperpp).eq.1.0) then
         CALL PRC(S1//'PP PERPF  OPT  1: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//
     >            ' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
      elseif (switch(swgperpp).eq.2.0) then
         CALL PRC(S1//'PP PERPF  OPT  2: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
      elseif (switch(swgperpp).eq.3.0) then
         CALL PRC(S1//'PP PERPF  OPT  3: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperpp).eq.4.0) then
         CALL PRC(S1//'PP PERPF  OPT  4: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE HALF FIELD LINE ->0')
         CALL PRC(SP//'AT THE MIDPOINT -  ')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED'//
     >            ' PERPENDICULAR FLUX.')
         CALL PRC(SP//'FOR UNDER-IONIZED CASES AND A SOURCE')
         CALL PRC(SP//'DISTRIBUTED PROPORTIONAL TO NE FROM THE')
         CALL PRC(SP//'PREVIOUS ITERATION FOR AN OVER-IONIZED')
         CALL PRC(SP//'CASE.')
      elseif (switch(swgperpp).eq.5.0) then
         CALL PRC(S1//'PP PERPF  OPT  5: TERM IS ON')
         CALL PRC(SP//'NET FLUX ALONG FIELD LINE ->0 AT THE'//
     >            ' MIDPOINT')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED PERPENDICULAR')
         CALL PRC(SP//'FLUX OVER THE 1/2 FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE')
         call prc(sp//'FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperpp).eq.6.0) then
         CALL PRC(S1//'PP PERPF  OPT  6: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER THE ENTIRE FILED LINE ->0')
         CALL PRC(SP//'USING AN EVENLY DISTRIBUTED PERPENDICULAR')
         CALL PRC(SP//'FLUX OVER THE WHOLE FIELD LINE')
         CALL PRC(SP//'AND A RECTANGULAR SOURCE BETWEEN THE')
         CALL PRC(SP//'FOLLOWING POINTS FROM EACH TARGET')
         CALL PRQ(SP//'START = SMAX * ',GPERPBEGF)
         CALL PRQ(SP//'END   = SMAX * ',GPERPENDF)
         CALL PRQ(SP//'FRACTION IN RECTANGULAR = ',GPERPFRAC)
         CALL PRQ(SP//'FRACTION UNIFORM        = ',1.0D0-GPERPFRAC)
      elseif (switch(swgperpP).eq.7.0) then
         CALL PRC(S1//'PERP FLUX OPT  7: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX PROPORTIONAL')
         CALL PRC(SP//'TO THE SECOND GRADIENT IN THE DENSITY')
         CALL PRC(SP//'FROM AN EDGE2D SOLUTION OR A PREVIOUS')
         CALL PRC(SP//'SOL22 SOLUTION. THIS OPTION IS CHANGED')
         CALL PRC(SP//'TO A UNIFORM DISTRIBUTION FOR ANY RINGS')
         CALL PRC(SP//'WHERE THE RATIO OF EITHER THE POSITIVE')
         CALL PRC(SP//'OR NEGATIVE CONTRIBUTION TO THE INTEGRATED')
         CALL PRC(SP//'VALUE FOR TEH RING EXCEEDS 5.0' )
      elseif (switch(swgperpP).eq.8.0) then
         CALL PRC(S1//'PERP FLUX OPT  8: TERM IS ON')
         CALL PRC(SP//'NET FLUX OVER ENTIRE FIELD LINE ->0')
         CALL PRC(SP//'USING A PERPENDICULAR FLUX CALCULATED')
         CALL PRC(SP//'FROM THE SECOND GRADIENT IN THE DENSITY')
         call PRC(SP//'AND AN IMPOSED VALUE FOR THE DIFFUSION')
         CALL PRR(SP//'RATE:  DPERP [M2/S]  = ',CDPERP)
         CALL PRC(SP//'ANY REMAINING EXECESS OR SHORTFALL OF')
         CALL PRC(SP//'FLUX IS COMPENSATED FOR BY THE APPLICATION')
         CALL PRC(SP//'OF AN ADDITIONAL UNIFORM SOURCE OVER THE')
         CALL PRC(SP//'ENTIRE FLUX TUBE.')
      endif
c
c     Extra Gperp source/sink
c
      if (switch(swextra).eq.0.0) then
         CALL PRC(S1//'EXTRA PERP FLUX 0 : TERM IS OFF')
      elseif (switch(swextra).eq.1.0) then
         CALL PRC(S1//'EXTRA PERP FLUX 1 : TERM IS ON')
         CALL PRC(SP//'AN EXTRA SOURCE AND SINK ARE SUPERIMPOSED ON')
         CALL PRC(SP//'THE FLOW PATTERN FOR THE FLUX TUBE. THIS')
         CALL PRC(SP//'SOURCE AND SINK EXACTLY CANCEL BUT WILL')
         CALL PRC(SP//'AFFECT THE FLOW PATTERN ON THE FLUX TUBE')
         call prq(sp//'SOURCE SIZE = (TOTAL TARGET FLUX ON RING) * ',
     >                 gextra_mult)
         call prc(sp//'THE SOURCE IS IMPOSED OVER THE REGION:')
         call prq2(sp//'   - SMAX  * [F1,F2]  = ',gextra_src_start,
     >                                            gextra_src_stop)
         call prc(sp//'THE SINK IS IMPOSED OVER THE REGION:')
         call prq2(sp//'   - SMAX  * [F1,F2]  = ',gextra_sink_start,
     >                                            gextra_sink_stop)
      endif
c
      call prb
c
c     Major Radius Option
c
c
      call prb
c
      if (switch(swmajr).eq.0.0) then
         CALL PRC(S1//'MAJOR RADIUS   0: MAJOR RADIUS FACTOR IS OFF')
      elseif (switch(swmajr).eq.1.0) then
         CALL PRC(S1//'MAJOR RADIUS   1: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL TARGET FLUX ADJUSTED BY RTARG/R0')
      elseif (switch(swmajr).eq.2.0) then
         CALL PRC(S1//'MAJOR RADIUS   2: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL IONIZATION ADJUSTED BY RCELL/R0')
      elseif (switch(swmajr).eq.3.0) then
         CALL PRC(S1//'MAJOR RADIUS   3: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'ALL IONIZATION ADJUSTED BY R0/RCELL')
      elseif (switch(swmajr).eq.4.0) then
         CALL PRC(S1//'MAJOR RADIUS   4: MAJOR RADIUS FACTOR IS ON')
         CALL PRC(SP//'GENERALIZED R-CORRECTION APPLIED')
      endif
c
c     Core Flux option
c
      if (switch(swcore).ne.0.0) then
         call prb
         CALL PRC(S1//'CORE IONIZATION CORRECTION TERM'//
     >            ' IS NOT IMPLEMENTED')
      endif
c
c      if (switch(swcore).eq.0.0) then
c         call prc('       Cross-field Core Ionization Source is OFF')
c      elseif (switch(swcore).eq.1.0) then
c         call prc('       Cross-field Core Ionization Source is ON')
c         call prq('       - Core source fraction : ',corefrac)
c         call prc('       - Core source flux is distribited linearly')
c         call prc('         along the entire ring to the midpoint')
c      elseif (switch(swcore).eq.2.0) then
c         call prc('       Cross-field Core Ionization Source is ON')
c         call prq('       - Core source fraction : ',corefrac)
c         call prc('       - Core source flux is distribited linearly')
c         call prc('         from the X-point region to the midpoint')
c      endif
c
c
c     Recombination term
c
      call prb
c
      if (switch(swrecom).eq.0.0) then
         CALL PRC(S1//'RECOMB OPT    0 : TERM IS OFF')
         CALL PRC(SP//'NO  RECOMBINATION PARTICLE SOURCE')
      elseif (switch(swrecom).eq.1.0) then
         CALL PRC(S1//'RECOMB OPT    1 : TERM IS ON')
         CALL PRC(SP//'RECOMBINATION PARTICLE SOURCE ADDED TO FLUX')
         CALL PRR(SP//'TE CUTOFF FOR RECOMBINATION = ',TRECCUT)
         CALL PRC(SP//'RECOMBINATION RATES ARE LIMITED BY THE'//
     >                             ' T LISTED.')
      elseif (switch(swrecom).eq.2.0) then
         CALL PRC(S1//'RECOMB OPT    2 : TERM IS ON')
         CALL PRC(SP//'RECOMBINATION PARTICLE SOURCE ADDED TO FLUX')
         CALL PRC(SP//'RECOMBINATION CALCULATED FROM EDGE2D BG'
     >                     //' INPUT')
      endif
c
      call prb
c
c     Smoothing Option
c
      if (switch(swsmooth).eq.0.0) then
         CALL PRC(S1//'SMOOTHING OPT 0 : SMOOTHING IS'//
     >            ' TURNED OFF')
         CALL PRC(sp//'VALUES WILL NOT MATCH AT THE'//
     >            ' MID-POINT')
      elseif (switch(swsmooth).eq.1.0) then
         CALL PRC(S1//'SMOOTHING OPT 1 : SMOOTHING IS'//
     >            ' TURNED ON')
         CALL PRC(SP//'BG QUANTITIES ARE ADJUSTED'//
     >            ' TO MATCH AT THE MID-POINT')
      endif

      call prb
c
c     Detached Solution Option
c
      if (switch(swdetach).eq.0.0) then
         call prc(s1//'DETACHED  OPT 0 : DETACHED PLASMA'//
     >                   ' OVERRIDE OPTION IS OFF')
      elseif (switch(swdetach).eq.1.0) then
         call prc(s1//'DETACHED  OPT 1 : DETACHED PLASMA MODEL'//
     >                           ' OVERRIDES')
         call prc(sp//'SOL22 FOR '//outer//
     >                              ' HALF RING (IK=1)')
c
       if (s21refsw.eq.0) then
        call prc(sp//'NOTE: REF OPTION 0 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc(sp//'NOTE: REF OPTION 1 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of PMAX. (MF=PMAX)')
        call prc(sp//'Converted to S for each'//
     >                                 ' ring.')
       elseif (s21refsw.eq.2) then
        call prc(sp//'NOTE: REF OPTION 2 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc(sp//'NOTE: REF OPTION 3 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along P.(MF=1.0)')
        call prc(sp//'Converted to S for each'//
     >                                 ' ring.')
       endif
c
c      OUTER
c
       call prc ('                       '//outer//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rat)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rat,l2rat
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rat)
       call prr ('                       Te at end of A: Te0 * ',
     >            terat)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirat)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrat)
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalph)
       call prr ('                       Radiated power in B (Qrad): D0
     >* ',qrat)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases lineraly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(D0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = D0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrat)
c
      elseif (switch(swdetach).eq.2.0) then
c
         call prc(s1//'DETACHED  OPT 2 : DETACHED PLASMA MODEL'//
     >                           ' OVERRIDES')
         call prc(sp//'SOL22 FOR'//INNER//
     >                              ' HALF RING (IK=NKS(IR))')
c
       if (s21refsw.eq.0) then
        call prc(sp//'NOTE: REF OPTION 0 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc(sp//'NOTE: REF OPTION 1 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in units of PMAX. (MF=PMAX)')
        call prc(sp//'Converted to S for each'//
     >                                 ' ring.')
       elseif (s21refsw.eq.2) then
        call prc(sp//'NOTE: REF OPTION 2 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc(sp//'NOTE: REF OPTION 3 - ')
        call prc(sp//'All distances are specified')
        call prc(sp//'in meters along P.(MF=1.0)')
        call prc(sp//'Converted to S for each'//
     >                                 ' ring.')
       endif
c
c
c      INNER
c
       call prb
       call prc ('                       '//INNER//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rati)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rati,l2rati
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rati)
       call prr ('                       Te at end of A: Te0 * ',
     >            terati)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirati)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrati)
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalphi)
       call prr ('                       Radiated power in B (Qrad): Q0
     >* ',qrati)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases linearly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrati)
c
      endif
c
      call prb
c
c     Error Correction
c
      if (switch(swerror).eq.0.0) then
         CALL PRC(S1//'ERROR OPTION  0 : AUTOMATIC ERROR'//
     >                ' CORRECTION IS TURNED OFF')
      elseif (switch(swerror).gt.0.0) then
         CALL PRC(S1//'ERROR OPTION  N : AUTOMATIC ERROR'//
     >                ' CORRECTION IS TURNED ON')
         CALL PRR(SP//'ERROR CORRECTION STARTS AT LEVEL :',
     >                    switch(swerror))
      endif
c
c     End of OPTIONS section
c
      call prb
      call prc(s1//'SOL 22: REPORT OF RESULTS')
      call prb

c
c     Summary of over/under ionization
c
      if (switch(swgperp).eq.2.or.switch(swgperp).eq.3.or.
     >    switch(swgperp).eq.7.or.switch(swgperp).eq.8.or.
     >    switch(swgperpp).eq.2.or.switch(swgperpp).eq.3.or.
     >    switch(swgperpp).eq.7.or.switch(swgperpp).eq.8
     >    ) then


         call prc(s1//'SUMMARY OF OVER/UNDER IONIZATION BY RING:')
         call prb

         if (switch(swgperp).eq.7.0.or.switch(swgperpp).eq.7.0) then
            call prc('     IR        RCONST    '//
     >              '  Gperp RATIO    Gperp OPTION USED')

            do ir = irsep,nrs
c
               if (ir.le.irwall) then
                  swtmp = switch(swgperp)
               else
                  swtmp = switch(swgperpp)
               endif
c
               if (swtmp.eq.7.0.and.gperprat(ir).gt.5.0) swtmp = 2.0
c
               write(coment,
     >               '(3x,i5,4x,1p,g13.6,0p,2(4x,f7.2))')
     >               ir,rconst(ir),gperprat(ir),swtmp
               call prc(coment)
            end do
c
         else
c
            call prc('     IR        RCONST     Gperp OPTION USED')
c
            do ir = irsep,nrs
c
               if (ir.le.irwall) then
                  swtmp = switch(swgperp)
               else
                  swtmp = switch(swgperpp)
               endif
c
               write(coment,
     >               '(3x,i5,4x,1p,g13.6,0p,4x,f7.2)')
     >               ir,rconst(ir),swtmp
               call prc(coment)
            end do
c
         endif
c
         call prb
c
      endif

c
c     Rings which were fixed by default conditions.
c
      if (switch(swerror).eq.0.0) then
c
         cnt = 0
         do ir = irsep,irlim
            if (cerr(ir,1).ne.0.or.cerr(ir,2).ne.0) then
               cnt = cnt + 1
               if (cnt.eq.1) then
                  CALL PRC(S1//'SUMMARY OF ERRORS:')
                  CALL PRC(S1//
     >                    '   - SOLVER HAD PROBLEMS WITH THESE RINGS:')
                  CALL PRC(S1//'RING        CODE   DESCRIPTION'//
     >                     '       POSITION     ERROR OPTION')
               endif
c
               if (cerr(ir,2).ne.0) then
                  write(coment,1060) ir,outer,
     >                            cerr(ir,2),errstr(cerr(ir,2)),
     >                            cserr(ir,2),cdeferropt(ir,2)
                  call prc(coment)
               endif
c
               if (cerr(ir,1).ne.0) then
                  write(coment,1060) ir,inner,
     >                               cerr(ir,1),errstr(cerr(ir,1)),
     >                               cserr(ir,1),cdeferropt(ir,1)
                  call prc(coment)
               endif
c
c
            endif
         end do
c
         call prb
c
 1060    format(6x,i4,1x,a5,':',i4,3x,a,1x,g13.6,1x,f7.1)

      elseif (switch(swerror).ne.0.0) then
c
         CALL PRB
         CALL PRC (S1//'ERROR CORRECTION OPTION ACTIVATED:')
         CALL PRB

c
c        Summarize corrected rings.
c
         cnt = 0
         do ir = irsep,irlim
            if (cdeferr(ir,1).ne.0.or.cdeferr(ir,2).ne.0) then
               cnt = cnt + 1
               if (cnt.eq.1) then
c
c                 Print Error Correction description
c
                  call prerrdesc(switch(swerror))
c
                  call PRC(S1//'     ERROR SOLVER HAD A PROBLEM'//
     >                  ' WITH THESE RINGS:')
C
                  CALL PRC(S1//'RING        CODE   DESCRIPTION'//
     >                     '       POSITION     ERROR OPTION')

c
               endif
c
               if (cdeferr(ir,2).ne.0) then
                  write(coment,1060) ir,outer,cdeferr(ir,2),
     >                            errstr(cdeferr(ir,2)),
     >                            cdefserr(ir,2),cdeferropt(ir,2)
                  call prc(coment)
               endif
c
               if (cdeferr(ir,1).ne.0) then
                  write(coment,1060) ir,inner,cdeferr(ir,1),
     >                            errstr(cdeferr(ir,1)),
     >                            cdefserr(ir,1),cdeferropt(ir,1)
                  call prc(coment)
               endif
c
            endif
         end do
c
         if (cnt.eq.0) then
c
            CALL PRC (S1//'NO RINGS REQUIRED ERROR'//
     >                    ' CORRECTION METHODS')
c
         endif
c
         call prb
c
      end if
c
c     Summary of Mach number results
c
      if (switch(swmach).eq.1.0.or.switch(swmach).eq.2.0) then
c
         call prc(s1//'MACH SOLVER ACTIVATED')
         CALL PRC(s1//'SUMMARY OF TARGET MACH NUMBERS:')
c
         do ir = irsep, irlim
            write(coment,1000) ir,inner,cmachno(ir,1),
     >                  outer,cmachno(ir,2)
            call prc(coment)
         end do
c
         call prb
c
         CALL PRQ (S1//'DELTA M0 FOR FIRST MACH ITERATION  : ',
     >             DELTAM0)
         CALL PRQ (S1//'ULTIMATE RESOLUTION OF MACH NUMBER : ',M0RES)
c
c        Format for Mach number table
c

1000     format('       Ring : ',i4,'  ',a5,' M#: ',f12.6,
     >          '  ',a5,' M#: ',f12.6)
c
      elseif (switch(swmach).eq.3.0) then
c
         call prc(s1//'FIXED TARGET MACH NUMBERS:')
         CALL PRC(s1//'SUMMARY OF TARGET MACH NUMBERS:')
c
         do ir = irsep, irlim
            write(coment,1000) ir,inner,cmachno(ir,1),
     >                            outer,cmachno(ir,2)
            call prc(coment)
         end do
c
      endif
c
      call prb
      call prc(s1//'SUMMARY OF INTEGRATED CONDUCTIVE POWER RATIOS:')
      call prb
c
      CALL PRC(S1//'RING        CONDe/Qetot  CONDi/Qitot'//
     >                         '  COND/Qtot')
c
      do ir = irsep,irlim
c
c        Target 2 - outer for JET geometries
c
         write(coment,1061) ir,outer,sol22_power_ratio(ir,2,1),
     >             sol22_power_ratio(ir,2,2),sol22_power_ratio(ir,2,3)
         call prc(coment)
c
c        Target 1 - inner for JET geometries
c
         write(coment,1061) ir,inner,sol22_power_ratio(ir,1,1),
     >             sol22_power_ratio(ir,1,2),sol22_power_ratio(ir,1,3)
         call prc(coment)
c
      end do

 1061    format(6x,i4,' ',a5,':',3(1x,f11.7))

      call prb
c
c
c
c     Output controls
c
c      call prb
c      if (graph.eq.0) then
c         call prc ('Plotting is turned OFF')
c      elseif (graph.eq.1) then
c         call prc ('Plotting is turned ON')
c      endif
c
c      if (graphaux.eq.0) then
c         call prc ('Auxilliary plots and tables are turned OFF')
c      elseif (graphaux.eq.1) then
c         call prc ('Auxilliary plots and tables are turned ON')
c      endif
c
c      if (graphvel.eq.0) then
c         call prc ('Velocity plots and tables are turned OFF')
c      elseif (graphvel.eq.1) then
c         call prc ('Velocity plots and tables are turned ON')
c      endif
c
c      call prb
c      if (graphran.eq.0.0) then
c         call prc ('Close up plots are turned OFF')
c      elseif (graphvel.eq.1) then
c         call prc ('Close up plots are turned ON')
c         call prr ('Range for close up plots is 0.0  to ',graphran)
c      endif
c
c
      return
      end
c
      subroutine echosolorg(spts,npts)
      implicit none
c
c     This subroutine echoes the input values to standard out - it also
c     prints the additional calculated values - and is followed by the
c     tabular output of s,te,ti,n and v at each point S ... which would
c     be suitable for plotting on a spreadsheet or may be plotted by
c     calling GHOST routines if the graph option is set.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer npts,i
      real*8 spts(mxspts)
c
      CHARACTER  COMENT*77
c
c
      call prbs
      call prs(title)
      call pris('Ring number: ',ringnum)
      call prbs
c
      CALL PQS (' Denisty at the target:             ', n0)
      CALL PQS (' Electron temperature at the target:', te0)
      CALL PQS (' Ion temperature at the target:     ', ti0)
      CALL PQS (' Plasma ion mass (amu)              ', mb)
      CALL PQS (' Flow velocity at the target        ', v0)
      CALL PQS (' Electron heat coefficient (K0e)    ', k0e)
      CALL PQS (' Ion heat coefficient      (K0i)    ', k0i)
c
c     Mach number options
c
      call prbs
      if (switch(swmach).eq.0.0) then
         call prs ('Iterative Mach number solution - OFF')
         CALL PQS (' Imposed mach number at the target  ', m0)
      elseif (switch(swmach).eq.1.or.switch(swmach).eq.2) then
         call prs ('Iterative Mach number solution - ON ')
         CALL PQS (' Initial mach number at the target  ',origm0)
         CALL PQS (' Final mach number at the target    ',m0)
         call pqs (' Initial Target Density             ',netarg)
         call pqs (' Final  Target Density              ',
     >             nefinal)
         CALL PQS (' Delta m0 for first iteration       ',
     >             deltam0)
         CALL PQS (' Ultimate resolution of m0          ',m0res)
         call pris(' Number of iterations to obtain m0  ',miter)
         call pqs (' Supersonic transition point (m)    ',
     >               lastiters)
      elseif (switch(swmach).eq.3.0) then
         call prs ('Iterative Mach number solution - OFF')
         CALL PrS (' MACH numbers imposed at target vary.')
      endif
c
      call prbs
      call prs ('Other calculated quantities:        ')
      call prbs
c
      call pqs (' Gammae - electron power coefficient', gammae)
      call pqs (' Gammai - ion power coefficient     ', gammai)
      call pqs (' Gamma correction factor in gammai  ', gamcor)
      call pqs (' Electron power to target           ', pae)
      call pqs (' Ion power to target                ', pai)
      call pqs (' Pressure at infinity(w.mach effect)', pinf)
      call pqs (' Value of Smax for ring (m)         ',ringlen)
c
      call prbs
      call prs ('Set of activated switches           ')
      call prbs
c
      if (switch(swcond).eq.0.0) then
         call prs(' 5/2 nv * kT     term  is OFF        ')
      else
         call prs(' 5/2 nv * kT     term  is ON         ')
      endif
c
      if (switch(swconv).eq.0.0) then
         call prs(' 1/2 m v^3 * n   term  is OFF        ')
      else
         call prs(' 1/2 m v^3 * n   term  is ON         ')
      endif
c
      if (switch(swprad).eq.0.0) then
         call prs(' Prad            term  is OFF        ')
      else
         call prs(' Prad            term  is ON         ')
      endif
c
      if (switch(swphelp).eq.0.0) then
         call prs('       Phelpi          term  is OFF        ')
      elseif (switch(swphelp).eq.1.0) then
         call prs('       Internal Phelpi term  is ON         ')
      elseif (switch(swphelp).eq.2.0) then
         call prs('       PINQE           term  is ON         ')
      endif
c
      if (switch(swpei).eq.0.0) then
         call prs(' Pei             term  is OFF        ')
      elseif (switch(swpei).eq.1.0) then
         call prs(' Pei             term  is ON         ')
      elseif (switch(swpei).eq.3.0) then
         call prs(' Pei             term  is Calculated (not used)')
      endif
c
      if (switch(swpcx).eq.0.0) then
         call prs('       Pcx             term  is OFF        ')
      elseif (switch(swpcx).eq.1.0) then
         call prs('       Internal Pcx    term  is ON         ')
         call pqs('       CX power coefficeint CEICF          ',ceicf)
      elseif (switch(swpcx).eq.2.0) then
         call prs('       PINQI           term  is ON         ')
      endif
c
      if (switch(swppelec).eq.0.0) then
         call prs(' PPelec          term  is OFF        ')
      elseif (switch(swppelec).eq.1.0) then
         call prs(' PPelec          term  is ON to Xpoint')
      elseif (switch(swppelec).eq.2.0) then
         call pqs(' PPelec          term  is ON to SMAX *',pp_pow_dist)
      endif
c
      if (switch(swppion).eq.0.0) then
         call prs(' PPion           term  is OFF        ')
      elseif (switch(swppion).eq.1.0) then
         call prs(' PPion           term  is ON to Xpoint')
      elseif (switch(swppion).eq.2.0) then
         call pqs(' PPion           term  is ON to SMAX *',pp_pow_dist)
      endif
c
      if (switch(swppress).eq.0.0) then
         call prs(' PP Press        term  is OFF        ')
      elseif (switch(swppress).eq.1.0) then
         call prs(' PP Press        term  is ON to Xpoint')
      elseif (switch(swppress).eq.2.0) then
         call pqs(' PP Press        term  is ON to SMAX *',pp_pow_dist)
      endif
c
c
c
      if (switch(swvisc1).eq.0.0) then
         call prs('       Density Viscosity correction term'//
     >            '  is NOT IMPLEMENTED')
      endif
c
c
      if (switch(swnmom).eq.0.0) then
         call prs('       Neutral BG Momentum loss term      is OFF')
      elseif (switch(swnmom).eq.1.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Smom = Smom0    S < L * Smax        ')
         call prs('              = 0        S > L * Smax        ')
         call prs('         Smom0= Pt/(L * Smax) * (1/Ffric -1)')
         call pqs('         Ffric= ',ffric)
         call pqs('         L    = ',lenmom)
      elseif (switch(swnmom).eq.2.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Smom = Smom0 * exp ( -S / (Lm * Smax))  for'
     >            //' S < L * Smax')
         call prs('              = 0        S > L * Smax        ')
         call prs('         Smom0= Pt/(Lm * Smax)*(1/Ffric -1) / expf')
         call prs('         expf = (1 - exp  (-(L*Smax)/(Lm * Smax)))')
         call pqs('         Ffric= ',ffric)
         call pqs('         Lm   = ',lammom)
         call pqs('         L    = ',lenmom)
      elseif (switch(swnmom).eq.3.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Smom = fact * Siz(S) / INT (0 to L) '//
     >            '[ Siz(S'') dS'' ] ')
         call prs('         fact = Pt *  (1/Ffric -1) ')
         call pqs('         Ffric= ',ffric)
         call prs('         L    = Smax / 2.0')
      elseif (switch(swnmom).eq.4.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Smom = -m * Vb * Rcxiz *  Siz(S) ')
         call pqs('         Rcxiz= ',rcxmom)
      elseif (switch(swnmom).eq.5.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Smom is read directly from NIMBUS')

      elseif (switch(swnmom).eq.6.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Neutral BG Momentun loss term is taken')
         call prs('         from PIN results except for the first')
         call prs('         iteration, in which case: ')
         call prs('         Neutral BG Momentum loss term is OFF')
      elseif (switch(swnmom).eq.7.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Neutral BG Momentun loss term is taken')
         call prs('         from PIN results except for the first')
         call prs('         iteration, in which case: ')
         call prs('         Initial neutral BG Momentum loss term is:')
         call prs('         Smom = Smom0    S < L * Smax        ')
         call prs('              = 0        S > L * Smax        ')
         call prs('         Smom0= Pt/(L * Smax) * (1/Ffric -1)')
         call pqs('         Ffric= ',ffric)
         call pqs('         L    = ',lenmom)
      elseif (switch(swnmom).eq.8.0) then
         call prs('       Neutral BG Momentum loss term      is ON')
         call prs('         Neutral BG Momentun loss term is taken')
         call prs('         from PIN results except for the first')
         call prs('         iteration, in which case: ')
         call prs('         Initial neutral BG Momentum loss term is:')
         call prs('         Smom = Smom0 * exp ( -S / (Lm * Smax))  for'
     >            //' S < L * Smax')
         call prs('              = 0        S > L * Smax        ')
         call prs('         Smom0= Pt/(Lm * Smax)*(1/Ffric -1) / expf')
         call prs('         expf = (1 - exp(-(L * Smax)/(Lm * Smax)))')
         call pqs('         Ffric= ',ffric)
         call pqs('         Lm   = ',lammom)
         call pqs('         L    = ',lenmom)
      endif
c
c     More switches
c
      if (switch(swe2d).eq.0.0) then
         call prs('       Edge2D BG compatibility is OFF')
      elseif (switch(swe2d).eq.1.0) then
         call prs('       Edge2D BG compatibility is ON')
         call prs('       - Edge2D data is read for first cell')
      elseif (switch(swe2d).eq.2.0) then
         call prs('       Edge2D BG compatibility is ON')
         call prs('       - Edge2D data is read for first cell')
         call prs('       - Edge2D Flux is matched at first')
         call prs('         cell centre')
      endif
c
      if (switch(swpow).eq.0.0) then
         call prs('       Distributed Power Influx is OFF')
      elseif (switch(swpow).eq.1.0) then
         call prs('       Distributed Power Influx is ON')
         call prs('       - Target Power flux is distribited linearly')
         call prs('         along the entire ring to the midpoint')
      elseif (switch(swpow).eq.2.0) then
         call prs('       Distributed Power Influx is ON')
         call prs('       - Target Power flux is distribited linearly')
         call prs('         from the X-point region to the midpoint')
      endif
c
      if (switch(swgperp).eq.0.0) then
         call prs('       Perpendicular Flux Correction is OFF')
      elseif (switch(swgperp).eq.1.0) then
         call prs('       Perpendicular Flux is ON')
         call prs('       - Net flux along field line ->0 at the'//
     >            ' midpoint')
         call prs('         using an evenly distributed'//
     >            ' perpendicular flux.')
      endif
c
      if (switch(swgperpp).eq.0.0) then
         call prs('       PP Perpendicular Flux Correction is OFF')
      elseif (switch(swgperpp).eq.1.0) then
         call prs('       PP Perpendicular Flux is ON')
         call prs('       - Net flux along field line ->0 at the'//
     >            ' midpoint')
         call prs('         using an evenly distributed'//
     >            ' perpendicular flux.')
      endif
c
      if (recfrac.ne.1.0) then
         call prs('       Recycling Source Fraction is ON')
         call pqs('       - Recycling source fraction = ',recfrac)
      endif
c
c      if (switch(swcore).eq.0.0) then
c         call prs('       Cross-field Core Ionization Source is OFF')
c      elseif (switch(swcore).eq.1.0) then
c         call prs('       Cross-field Core Ionization Source is ON')
c         call pqs('       - Core source fraction : ',corefrac)
c         call prs('       - Core source flux is distribited linearly')
c         call prs('         along the entire ring to the midpoint')
c      elseif (switch(swcore).eq.2.0) then
c         call prs('       Cross-field Core Ionization Source is ON')
c         call pqs('       - Core source fraction : ',corefrac)
c         call prs('       - Core source flux is distribited linearly')
c         call prs('         from the X-point region to the midpoint')
c      endif
c
c
c
      if (switch(swmajr).eq.0.0) then
         call prs('       Major Radius Factor is OFF')
      elseif (switch(swmajr).eq.1.0) then
         call prs('       Major Radius Factor is ON')
         call prs('       - All Target Flux adjusted by Rtarg/R0')
      elseif (switch(swmajr).eq.2.0) then
         call prs('       Major Radius Factor is ON')
         call prs('       - All ionization adjusted by Rcell/R0')
      elseif (switch(swmajr).eq.3.0) then
         call prs('       Major Radius Factor is ON')
         call prs('       - All ionization adjusted by R0/Rcell')
      elseif (switch(swmajr).eq.4.0) then
         call prs('       Major Radius Factor is ON')
         call prs('       - Generalized R-correction applied')
      endif
c
c     Ionization Source characteristics
c
      if (switch(swion).eq.0.0) then
c
      call prbs
      call prs ('     Ionization Source Characteristics  S(s)')
      call prs ('       Exponential Decay : ')
      call prbs
c
      call pqs ('       Length of ionization source        ', lensfi)
      call pqs ('       Decay length of ionization source  ', lams)
c
c      call pqs (' Base strength of ionization source ', s0)
c
      elseif (switch(swion).eq.1) then
c
      call prbs
      call prs ('     Ionization Source Characteristics  S(s)')
      call prbs
      call prs ('       Data for Ionization Source is returned')
      call prs ('       from a PIN/NIMBUS run and is linearly')
      call prs ('       interpolated for values between grid')
      call prs ('       points. The ionization source is')
      call prs ('       integrated over the half-ring and set')
      call prs ('       equal to the flow to the target NoVo')
      call prs ('       as a Normalization factor.')
c
      call pqs ('       Normalization factor = ',fnorm)
c
      elseif (switch(swion).eq.2) then
c
      call prbs
      call prs ('     Ionization Source Characteristics  S(s)')
      call prbs
      call prs ('       Data for Ionization Source is returned')
      call prs ('       from a PIN/NIMBUS run and is linearly')
      call prs ('       interpolated for values between grid')
      call prs ('       points. The data is used AS-IS and is')
      call prs ('       NOT Normalized to the target flux for')
      call prs ('       the ring.')
c
      elseif (switch(swion).eq.3.0) then
c
      call prs('       Triangular Ionization Source:  S(s)')
      call pqs('          Extending from   :',lensst)
      call pqs('                    to     :',lensfi)
      call prs('          Integral of source normalized'
     >                      //' to ring target flux.')
c
      elseif (switch(swion).eq.3.0) then
c
      call prs('       Rectangular (Flat) Ionization Source: S(s)')
      call pqs('          Extending from   :',lensst)
      call pqs('                    to     :',lensfi)
      call prs('          Integral of source normalized'
     >                      //' to ring target flux.')
c
      endif
c
      if (switch(swprad).eq.1.0) then
      call prbs
      call prs ('Radiation Source Characteristics  Prad')
      call prbs
c
      call pqs (' Length of radiation source         ', lenr)
      call pqs (' Decay length of radiation  source  ', lamr)
      call pqs (' Source strength fraction (frr)     ', frr)
      call pqs (' Strength of radiation source       ', prad0)

      elseif (switch(swprad).eq.6.0) then
      call prbs
      call prs ('Radiation Source Characteristics  Prad')
      call prbs
c
      call pqs (' End Length of radiation source     ', lenr)
      call pqs (' Start length of radiation  source  ', lamr)
      call pqs (' Source strength fraction (frr)     ', frr)
      call pqs (' Strength of radiation source       ', prad0)
c
      endif
c
c     Output controls
c
      call prbs
      if (graph.eq.0) then
         call prs ('Plotting is turned OFF')
      elseif (graph.eq.1) then
         call prs ('Plotting is turned ON')
      endif
c
      if (graphaux.eq.0) then
         call prs ('Auxilliary plots and tables are turned OFF')
      elseif (graphaux.eq.1) then
         call prs ('Auxilliary plots and tables are turned ON')
      endif
c
      if (graphvel.eq.0) then
         call prs ('Velocity plots and tables are turned OFF')
      elseif (graphvel.eq.1) then
         call prs ('Velocity plots and tables are turned ON')
      endif
c
      call prbs
      if (graphran.eq.0.0) then
         call prs ('Close up plots are turned OFF')
      elseif (graphvel.eq.1) then
         call prs ('Close up plots are turned ON')
         call prr ('Range for close up plots is 0.0  to ',graphran)
      endif
c
      call prbs
      call pris ('Initial number of Runge-Kutta steps between'//
     >          ' each S:', ndiv)
      call prbs
c
      call prs ('Set of S-values for Background calculations')
      call prbs
      do i = startn, npts
         call pqs('    S = ', spts(i))
      end do
      call prbs
      call prs ('End of Specification Section')
      call prbs
      call prs ('Results : ')
c
      return
      end
C
C
C  *********************************************************************
C  *  PRS:  PRINTS A CHARACTER STRING                                  *
C  *********************************************************************
C
      SUBROUTINE PRS(STRING)
      CHARACTER STRING*(*)
      WRITE (21,'(1X,A)') STRING
      RETURN
      END
C
C  *********************************************************************
C  *  PRBS:  PRINTS A BLANK LINE                                       *
C  *********************************************************************
C
      SUBROUTINE PRBS
      WRITE (21,'(1X)')
      RETURN
      END
C
C  *********************************************************************
C  *  PRIS:  PRINTS AN INTEGER                                         *
C  *********************************************************************
C
      SUBROUTINE PRIS (NAME, I)
      CHARACTER NAME*(*)
      INTEGER   I
      WRITE (21,'(1X,A,I7)') NAME,I
      RETURN
      END
C
C  *********************************************************************
C  *  PQS:  PRINTS A REAL*8 NUMBER                                     *
C  *********************************************************************
C
      SUBROUTINE PQS (NAME, R)
      CHARACTER NAME*(*)
      real*8      R
      IF (ABS(R).LT.0.1.OR.ABS(R).GE.1000.0) THEN
        WRitE (21,'(1X,A,1P,G15.8)') NAME,R
      ELSEIF (ABS(R).LT.1.0) THEN
        WRITE (21,'(1X,A,F15.8)') NAME,R
      ELSE
        WRITE (21,'(1X,A,F15.8)') NAME,R
      ENDIF
      RETURN
      END
c
c
c
      subroutine mkplot(dspts,npts,dte,dti,dne,dvb,xflag)
      implicit none
c
c     Calls GHOST routines to plot the temperature and other results.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer npts,xflag
      real*8 dspts(mxspts),dte(mxspts),dti(mxspts),
     >                dne(mxspts),dvb(mxspts)
c
c     local variables
c
      real  spts(mxspts),te(mxspts),ti(mxspts),
     >                ne(mxspts),vb(mxspts)
      integer i,ibrok,l,ispot,ncases
      character*50 xlabel,ylabel
      character*36 table,ref(4)
      real plots(mxspts,maxcases),scales(maxcases)
      real cxmin,cxmax,cymin,cymax
      real fmax,fmin
      integer lenstr
      external lenstr,fmax,fmin
c
c     Initialization
c
      do i = 1, mxspts
         spts(i) = dspts(i)
         te(i) = dte(i)
         ti(i) = dti(i)
         ne(i) = dne(i)
         vb(i) = dvb(i)
      end do
c
      ncases = 4
      ISPOT  = 12
c
      scales(1) = fmax(te,npts,1)
      scales(2) = fmax(ti,npts,1)
      scales(3) = fmax(ne,npts,1)
      scales(4) = abs(fmin(vb,npts,1))
c
      xlabel = 'Distance along the field line (m)'
      ylabel = 'Normalized Quantities'
      table = 'Table of Values'
c
      write(ref(1),'(''Scale: Te = '',g13.6)') scales(1)
      write(ref(2),'(''Scale: Ti = '',g13.6)') scales(2)
      write(ref(3),'(''Scale: Ne = '',g13.6)') scales(3)
      write(ref(4),'(''Scale: Vb = '',g13.6)') -scales(4)
c
c     draw titles
c
      call drawtitles(title,table,xlabel,ylabel)
c
c     Draw Frame
c
      call drawframe
c
c     draw scales
c
      cymin = -1.0
      cymax =  1.0
      cxmin = 0.0
      if (xflag.eq.0) then
         cxmax = fmax(spts,npts,1)
      elseif (xflag.eq.1)  then
         cxmax = graphran
      endif
c
      call drawscales (cxmin,cxmax,cymin,cymax)
c
c
c     draw plots - load data into plot array
c
      do i = startn,npts
         plots(i,1) = te(i)
         plots(i,2) = ti(i)
         plots(i,3) = ne(i)
         plots(i,4) = vb(i)
      end do
c
      call plotdata(spts,plots,npts,ref,ncases,
     >              cxmin,cymin,cxmax,cymax,scales)
c
      return
      end
c
c
c
      subroutine mkauxplot(dspts,npts,dte,dti,dne,dvb,
     >                  dphelpiv,dpeiv,dpeiv2,dpcxv,dpradv,xflag)
      implicit none
c
c     This routine calculates the various quantities that go
c     into calculating Te, Ti, n,v and plots them as a function of
c     S.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer npts,xflag
      real*8 dspts(mxspts),dte(mxspts),dti(mxspts),
     >                 dne(mxspts),dvb(mxspts)
      real*8 dphelpiv(mxspts),dpeiv(mxspts),dpcxv(mxspts)
      real*8 dpeiv2(mxspts),dpradv(mxspts)
c
c     local variables
c
      real spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
      real phelpiv(mxspts),peiv(mxspts),pcxv(mxspts)
      real peiv2(mxspts)
      integer i,ibrok,l,ispot
      character*50 xlabel,ylabel
      character*36 table,ref1,ref2,ref3,ref4,ref5,ref6
      character*112 comment
      integer lenstr
      external lenstr
c
      real cxmin,cxmax,cymin,cymax
c
      real condi(mxspts),conde(mxspts),conv1i(mxspts)
      real conv2i(mxspts),convi(mxspts),conve(mxspts)
      real pradv(mxspts)
      real powi(mxspts),powe(mxspts)
      real sume(mxspts), sumi(mxspts)
      real totci(mxspts),totpi(mxspts)
      real totce(mxspts),totpe(mxspts)
      real*8 s,t1i,t1e,n
      real*8 fegrad, figrad,pcx,pei,phelpi,conv,cond
      external fegrad,figrad,pcx,pei,phelpi,conv,cond
      real*8 paes,pais
      external paes,pais
c
c     Copy everything into local variables
c
      do i = 1,mxspts
         spts(i)   = dspts(i)
         te(i)     =  dte(i)
         ti(i)     =  dti(i)
         ne(i)     =  dne(i)
         vb(i)     =  dvb(i)
         phelpiv(i)=  dphelpiv(i)
         peiv(i)   =  dpeiv(i)
         peiv2(i)  =  dpeiv2(i)
         pcxv(i)   =  dpcxv(i)
         pradv(i)  =  dpradv(i)
      end do
c
c     Calculate all the auxiliary values - at each S-point
c     And the scale for the plot.
c
      cymin =  1.0e10
      cymax = -1.0e10
      ISPOT  = 12
c
c
      do i = startn,npts
         s   = spts(i)
         t1e = te(i)
         t1i = ti(i)
         n   = ne(i)
c
c
         if (actswcond.eq.0.0) then
            conve(i) = 0.0
            conv1i(i) = 0.0
         else
            conve(i) = cond(s,t1e)
            conv1i(i) = cond(s,t1i)
         endif
c
         if (actswconv.eq.0.0) then
            conv2i(i) = 0.0
         else
            conv2i(i) = conv(s,n,t1i)
         endif
c
c         if (actswprad.eq.0.0) then
c            pradv(i)  = 0.0
c         else
c            pradv(i)  = prad(s)
c         endif
c
         powe(i) = paes(s)
         powi(i) = pais(s)
c
c
         convi(i) = conv1i(i) + conv2i(i)
c
         conde(i)= -(conve(i)+pradv(i)+
     >         phelpiv(i)+peiv(i)+powe(i))
c
         condi(i)  = - (convi(i)
     >               + pcxv(i)-peiv(i)+powi(i))
c
         sume(i) = conde(i) + conve(i) + pradv(i)
     >             + phelpiv(i) + peiv(i)
         sumi(i) = condi(i) + convi(i)
     >             + pcxv(i) - peiv(i)
c
         totci(i) = condi(i) + convi(i)
         totce(i) = conde(i) + conve(i)
c
         totpi(i) = pcxv(i) - peiv(i) + powi(i)
         totpe(i) = phelpiv(i) + peiv(i) + powe(i) + pradv(i)
c
      end do
c
c
c
      if (actswpei.eq.3) then
         do i = startn,npts
            peiv(i) = peiv2(i)
         end do
      endif
c
      if (xflag.eq.0) then
c
c     Print out tables of Auxiliary values
c
      if (actswpei.eq.3.0) then
         call prbs
         call prs('Pei values are estimates NOT included in Totals:')
      endif
c
      call prbs
      call prs('Tables of calculated SOL Auxiliary values - E')
      call prbs
      write(comment,200)
      call prs(comment)
      do i = startn, npts
         write(comment,100) spts(i),conde(i),conve(i),
     >                 pradv(i),phelpiv(i),peiv(i),sume(i),powe(i),
     >                 totce(i),totpe(i)
         call prs (comment)
      end do
c
      call prbs
      call prs('Tables of calculated SOL Auxiliary values - I')
      call prbs
      write(comment,300)
      call prs(comment)
      do i = startn, npts
         write(comment,100) spts(i),condi(i),conv1i(i),
     >                 conv2i(i), pcxv(i),-peiv(i),sumi(i),powi(i),
     >                 totci(i), totpi(i)
         call prs (comment)
      end do
c
      endif
c
c
c
 100  format (g10.3,9g11.4)
 200  format (5x,'S',9x,'Cond',7x,'Conv',7x,'Prad',
     >        7x,'Phelp',6x,'Pei',8x,'Sum E',6x,'Pow e',
     >        5x,'Cond+Conv',4x,'Srcs')
 300  format (5x,'S',9x,'Cond',7x,'Conv1',5x,'Kinetic',5x,'Pcx',
     >        8x,'Pei',8x,'Sum I',6x,'Pow I',
     >        5x,'Cond+Conv',4x,'Srcs')
c
c
c     Set up the graph. E-power balance.
c
c
      xlabel = 'Distance along the field line (m)'
      ylabel = 'E-power balance quantities'
      table = 'Table of Values'
c
      write(ref1,'(''Conduction term '')')
      write(ref2,'(''5/2nvkTe        '')')
      write(ref3,'(''Prad            '')')
      write(ref4,'(''Phelpi          '')')
      write(ref5,'(''Pei             '')')
      write(ref6,'(''Sum of all terms'')')
c
c     calculate cymin and cymax ... if they are the same - skip plotting
c
c
c     Find Y range for the E-power plot
c
      cymin =  1.0e30
      cymax = -1.0e30
c
      call fminmax(cymin,cymax,conde,npts)
c
      if (actswcond.ne.0.0)
     >   call fminmax(cymin,cymax,conve,npts)
c
      if (actswprad.ne.0.0)
     >   call fminmax(cymin,cymax,pradv,npts)
c
      if (actswpei.ne.0.0)
     >   call fminmax(cymin,cymax,peiv,npts)
c
      if (actswphelp.ne.0.0)
     >   call fminmax(cymin,cymax,phelpiv,npts)
c
      call fminmax(cymin,cymax,sume,npts)
c
c
      if (cymin.lt.cymax) then
C
C---- DRAW TITLES
C
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRMAG (20)
      CALL LINCOL (1)
      CALL THICK  (2)
      L = LENSTR(TITLE)
      CALL PCSCEN (0.8, 0.95, TITLE(:L))
      L = LENSTR(TABLE)
      CALL PCSCEN (1.18, 0.87, TABLE(:L))
      L = LENSTR (XLABEL)
      CALL PCSCEN (0.5,0.025,XLABEL(:L))
      CALL THICK  (1)
      CALL CTRMAG (12)
      L = LENSTR (REF1)
      CALL plotst (1.06, 0.77, REF1(:L))
      L = LENSTR (REF2)
      CALL plotst (1.06, 0.72, REF2(:L))
      L = LENSTR (REF3)
      CALL plotst (1.06, 0.67, REF3(:L))
      L = LENSTR (REF4)
      CALL plotst (1.06, 0.62, REF4(:L))
      L = LENSTR (REF5)
      CALL plotst (1.06, 0.57, REF5(:L))
      L = LENSTR (REF6)
      CALL plotst (1.06, 0.52, REF6(:L))
      CALL CTRMAG (ISPOT)
C
C---- DRAW FRAMES
C
      CALL LINCOL (3)
      CALL POSITN (0.1, 0.1)
      CALL   JOIN (0.1, 0.9)
      CALL   JOIN (0.9, 0.9)
      CALL   JOIN (0.9, 0.1)
      CALL   JOIN (0.1, 0.1)
      CALL POSITN (0.93, 0.1)
      CALL   JOIN (0.93, 0.9)
      CALL   JOIN (1.35, 0.9)
      CALL   JOIN (1.35, 0.1)
      CALL   JOIN (0.93, 0.1)
      CALL POSITN (0.93, 0.85)
      CALL   JOIN (1.35, 0.85)
      CALL POSITN (0.93, 0.15)
      CALL   JOIN (1.35, 0.15)
      CALL POSITN (0.93, 0.20)
      CALL   JOIN (1.35, 0.20)
      CALL POSITN (0.93, 0.25)
      CALL   JOIN (1.35, 0.25)
C
c     Draw scales - xscale
c
      CXMIN = spts(1)
      if (xflag.eq.0) then
         cxmax = spts(npts)
      elseif (xflag.eq.1)  then
         cxmax = graphran
      endif
c
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
      CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
      CALL CTRMAG (10)
      CALL XSCALE
C
c     yscale
c
c
c     Set Y-scales
c
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
      CALL LINCOL (3)
      CALL CTRMAG (10)
      CALL YSCALE
      CALL PSPACE (0.0, 1.35, 0.11, 0.89)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRORI (90.0)
      CALL LINCOL (1)
      CALL CTRMAG (14)
      CALL THICK  (2)
      L = LENSTR (YLABEL)
      CALL PCSEND (0.02,0.99,YLABEL(:L))
      CALL CTRORI (0.0)
      CALL THICK  (1)
C
      CALL FULL
C
C---- DRAW PLOTS - E-power balance
C
c     Cond
c
      call plotln(1,0.77,npts,spts,conde,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Conv
c
      call plotln(2,0.72,npts,spts,conve,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Prad
c
      call plotln(3,0.67,npts,spts,pradv,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Phelpi
c
      call plotln(4,0.62,npts,spts,phelpiv,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Pei
c
      call plotln(5,0.57,npts,spts,peiv,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Sum of E-balance terms
c
      call plotln(6,0.52,npts,spts,sume,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
      call frame
c
c     End of E-power balance plot
c
      endif
c
c
c
c     Set up the graph. I-power balance.
c
c
c     If cymin = cymax then skip the plot as there will be no range
c     to plot.
c
c     Find Y range for the plot
c
      cymin =  1.0e30
      cymax = -1.0e30
c
      call fminmax(cymin,cymax,condi,npts)
c
      if (actswcond.ne.0.0.or.actswconv.ne.0.0)
     >   call fminmax(cymin,cymax,convi,npts)
c
      if (actswpcx.ne.0.0)
     >   call fminmax(cymin,cymax,pcxv,npts)
c
      if (actswpei.ne.0.0)
     >   call fminmax(cymin,cymax,peiv,npts)
c
      call fminmax(cymin,cymax,sumi,npts)
c
c     Test cymin,cymax
c
      if (cymin.lt.cymax) then
c
c
      xlabel = 'Distance along the field line (m)'
      ylabel = 'I-power balance quantities'
      ISPOT  = 12
      table = 'Table of Values'
c
      write(ref1,'(''Conduction term '')')
      write(ref2,'(''Convection term '')')
      write(ref3,'(''Pcx             '')')
      write(ref4,'(''Pei             '')')
      write(ref5,'(''Sum of all terms'')')
C
C---- DRAW TITLES
C
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRMAG (20)
      CALL LINCOL (1)
      CALL THICK  (2)
      L = LENSTR(TITLE)
      CALL PCSCEN (0.8, 0.95, TITLE(:L))
      L = LENSTR(TABLE)
      CALL PCSCEN (1.18, 0.87, TABLE(:L))
      L = LENSTR (XLABEL)
      CALL PCSCEN (0.5,0.025,XLABEL(:L))
      CALL THICK  (1)
      CALL CTRMAG (12)
      L = LENSTR (REF1)
      CALL plotst (1.06, 0.77, REF1(:L))
      L = LENSTR (REF2)
      CALL plotst (1.06, 0.72, REF2(:L))
      L = LENSTR (REF3)
      CALL plotst (1.06, 0.67, REF3(:L))
      L = LENSTR (REF4)
      CALL plotst (1.06, 0.62, REF4(:L))
      L = LENSTR (REF5)
      CALL plotst (1.06, 0.57, REF5(:L))
      CALL CTRMAG (ISPOT)
C
C---- DRAW FRAMES
C
      CALL LINCOL (3)
      CALL POSITN (0.1, 0.1)
      CALL   JOIN (0.1, 0.9)
      CALL   JOIN (0.9, 0.9)
      CALL   JOIN (0.9, 0.1)
      CALL   JOIN (0.1, 0.1)
      CALL POSITN (0.93, 0.1)
      CALL   JOIN (0.93, 0.9)
      CALL   JOIN (1.35, 0.9)
      CALL   JOIN (1.35, 0.1)
      CALL   JOIN (0.93, 0.1)
      CALL POSITN (0.93, 0.85)
      CALL   JOIN (1.35, 0.85)
      CALL POSITN (0.93, 0.15)
      CALL   JOIN (1.35, 0.15)
      CALL POSITN (0.93, 0.20)
      CALL   JOIN (1.35, 0.20)
      CALL POSITN (0.93, 0.25)
      CALL   JOIN (1.35, 0.25)
C
c     Draw scales - xscale
c
      CXMIN = spts(1)
      if (xflag.eq.0) then
         cxmax = spts(npts)
      elseif (xflag.eq.1)  then
         cxmax = graphran
      endif
c
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
      CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
      CALL CTRMAG (10)
      CALL XSCALE
C
c     yscale
c
c
c     Set Y-scales
c
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
      CALL LINCOL (3)
      CALL CTRMAG (10)
      CALL YSCALE
      CALL PSPACE (0.0, 1.35, 0.11, 0.89)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRORI (90.0)
      CALL LINCOL (1)
      CALL CTRMAG (14)
      CALL THICK  (2)
      L = LENSTR (YLABEL)
      CALL PCSEND (0.02,0.99,YLABEL(:L))
      CALL CTRORI (0.0)
      CALL THICK  (1)
C
      CALL FULL
C
C---- DRAW PLOTS -  I-Power Balance
C
c     Cond
c
      call plotln(1,0.77,npts,spts,condi,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Conv
c
      call plotln(2,0.72,npts,spts,convi,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Pcx
c
      call plotln(3,0.67,npts,spts,pcxv,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Pei
c
      call plotln(4,0.62,npts,spts,peiv,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
c     Sum of I-balance terms
c
      call plotln(5,0.57,npts,spts,sumi,
     >            cxmin,cxmax,cymin,cymax,1.0)
c
      call frame
c
c     End of I-power plot
c
      endif
c
c
      return
      end
c
c
c
      subroutine mkvelplot(dspts,npts,dte,dti,dne,dvb,dne2,dga,xflag)
      implicit none
c
c     This routine calculates the various quantities that go
c     into calculating Te, Ti, n,v and plots them as a function of
c     S.
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      integer npts,xflag
      real*8 dspts(mxspts),dte(mxspts),dti(mxspts),
     >                 dne(mxspts),dvb(mxspts)
      real*8 dne2(mxspts),dga(mxspts)
c
c     local variables
c
      real spts(mxspts),te(mxspts),ti(mxspts),ne(mxspts),vb(mxspts)
      real ne2(mxspts),ga(mxspts)
      character*50 xlabel,ylabel
      character*36 table,ref(3)
      character*78 comment
c
      real fmin,fmax
      integer lenstr,i
      external lenstr,fmin,fmax
c
      real cxmin,cxmax,cymin,cymax
c
      real vplots(mxspts,maxcases),scales(maxcases)
      integer ncases
c
c     Initialize
c
c
c     Copy everything into local varaibles
c
      do i = 1,mxspts
         spts(i)   =  dspts(i)
         te(i)     =  dte(i)
         ti(i)     =  dti(i)
         ne(i)     =  dne(i)
         vb(i)     =  dvb(i)
         ne2(i)    =  dne2(i)
         ga(i)     =  dga(i)
      end do
c
      ncases = 3
      do i = 1,ncases
         scales(i) = 0.0
      end do
c
c     Calculate all the velocities (super and sub sonic at each S-point)
c     And the scale for the plot.
c
      cymin =  1.0e10
      cymax = -1.0e10
c
      write (6,*) 'vsuper:',lastiters
c
      do i = startn,npts
c         if ((spts(i).lt.lastiters).or.(spts(i).eq.0.0)) then
         if ((spts(i).lt.lastiters)) then
            if (ne2(i).gt.0.0) then
               vplots(i,1) = ga(i)/ne2(i)
            else
               vplots(i,1) = 0.0
            endif
            vplots(i,2) = vb(i)
         elseif (spts(i).ge.lastiters) then
            vplots(i,1) = vb(i)
            if (ne2(i).gt.0.0) then
               vplots(i,2) = ga(i)/ne2(i)
            else
               vplots(i,2) = 0.0
            endif
         endif
         vplots(i,3) = -sqrt((te(i)+ti(i))/mb * econv/mconv)
      end do
c
c     Print
c
      if (xflag.eq.0) then
c
c     Print out tables of the velocity values
c
      call prbs
      call pris('Tables of calculated SOL Velocity values ',ringnum)
      call prbs
      write(comment,300)
      call prs(comment)
      do i = startn, npts
         write(comment,100) spts(i),vplots(i,1),vplots(i,2),
     >                      vplots(i,3),ga(i),ne(i)
         call prs(comment)
      end do
c
      endif
c
c
c
 100  format (6g13.5)
 300  format (6x,'S',11x,'Vsub',10x,'Vsuper',9x,'Cs',8x,
     >        'Gamma',9x,'Ne')

c
c     Set up the graph. Velocities
c
c
      xlabel = 'Distance along the field line (m)'
      ylabel = 'Velocity in m/s'
      table = 'Table of Values'
c
c     draw titles
c
      call drawtitles(title,table,xlabel,ylabel)
c
c     Draw Frame
c
      call drawframe
c
c     Calculate Min and Max values
c
      cxmin = 0.0
      cymax = 0.0
c
c
c      cxmax = fmax(spts,npts,1)
c
      if (xflag.eq.0) then
          cxmax = 10.0 * lastiters
c
c         cxmax = spts(npts)
c
      elseif (xflag.eq.1)  then
         cxmax = graphran
      endif
c
c      cymin = fmin(vplots,npts,ncases)
c
c
c      Need to establish a more usable value of cymin ... perhaps
c      four times the minimum value of Vb
c
       cymin= 4.0 * fmin(vb,npts,1)
c
c
c     draw scales
c
      call drawscales (cxmin,cxmax,cymin,cymax)
c
      write(ref(1),'(''Sub-sonic branch'')')
      write(ref(2),'(''Sup-sonic branch'')')
      write(ref(3),'(''Sound Speed (Cs)'')')

c     draw plots
c
      call plotdata(spts,vplots,npts,ref,ncases,
     >              cxmin,cymin,cxmax,cymax,scales)
c
c
      return
      end
c
c
c
      subroutine fminmax(ymin,ymax,vals,npts)
      implicit none
c
c     Finds the maximum and minimum values in the
c     array vals ... if they are greater/less than the
c     values passed in in ymin,ymax.
c
      integer npts,i
      real ymin,ymax,vals(npts)
c
      do i = 1,npts
         ymin = amin1(ymin,vals(i))
         ymax = amax1(ymax,vals(i))
      enddo
      return
      end
c
c
c
      subroutine plotln(ibrok,spot,npts,x,y,xmin,xmax,
     >                  ymin,ymax,scale)
      implicit none
c
c     This subroutine draws a curve on the graph.
c
      integer ibrok,npts,i
      real spot,x(npts),y(npts),xmin,xmax,ymin,ymax,scale
c
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (XMIN,XMAX,YMIN,YMAX)
C
      CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
c
      CALL POSITN (x(1), y(1)/scale)
c
      DO  I = 2,npts
        CALL JOIN (x(i), y(i)/scale)
      end do
c
      CALL LINCOL (2)
      CALL CTRMAG (12)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL POSITN (0.95,SPOT)
      CALL JOIN   (1.03,SPOT)
      CALL FULL
      return
      end
c
c
c
      subroutine plotdata(sdata,pvalues,pts,cases,ncases,
     >                     cxmin,cymin,cxmax,cymax,scales)
      implicit none
c
c     This routine plots the data in the two dimensional array.
c
      include 'solparams'
      include 'solcommon'
c
      integer ncases,pts
      real cxmin,cxmax,cymin,cymax
      real sdata(mxspts)
      real pvalues(mxspts,maxcases),scales(maxcases)
      character*(*) cases(ncases)
c
      integer i,j
      real values(mxspts),startpos,diffpos,plotpos
c
      startpos = 0.77
      diffpos = 0.05
c
      do i = 1,ncases
         plotpos= startpos - (i-1) * diffpos
         do j = 1,pts
            values(j) = pvalues(j,i)
         enddo
         call plotvals(i,plotpos,pts,sdata,values,
     >                 cxmin,cxmax,cymin,cymax,cases(i),scales(i))
      enddo
c
      call frame
c
      return
      end
c
c
c
      subroutine plotvals(ibrok,spot,npts,x,y,xmin,xmax,
     >                  ymin,ymax,name,scale)
      implicit none
c
c     This subroutine draws a curve on the graph.
c
      integer ibrok,npts,i
      real spot,x(npts),y(npts),xmin,xmax,ymin,ymax,scale
      character*(*) name
      integer lenstr,l
      external lenstr
c
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (XMIN,XMAX,YMIN,YMAX)
C
      CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
c
      if (scale.eq.0.0) then
         CALL POSITN (x(1), y(1))

c
         DO  I = 2,npts
            CALL JOIN (x(i), y(i))
         end do
      elseif (scale.ne.0.0) then
         CALL POSITN (x(1), y(1)/scale)

c
         DO  I = 2,npts
            CALL JOIN (x(i), y(i)/scale)
         end do
      endif
c
      CALL LINCOL (2)
      CALL CTRMAG (12)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL POSITN (0.95,SPOT)
      CALL JOIN   (1.03,SPOT)
      CALL FULL
      CALL THICK  (1)
      CALL CTRMAG (12)
      L = LENSTR (name)
      CALL plotst (1.06, spot, name(:L))
      call full
c
      return
      end
c
c
c
      subroutine drawscales(cxmin,cxmax,cymin,cymax)
      implicit none
c
c     Draw scales
c
      real cxmin,cxmax,cymin,cymax
c
c     xscale
c
c
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
      CALL MAP    (CXMIN, CXMAX, 0.1, 0.9)
      CALL CTRMAG (10)
      CALL XSCALE
C
c     yscale
c
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (0.0, 1.0, CYMIN, CYMAX)
      CALL LINCOL (3)
      CALL CTRMAG (10)
      CALL YSCALE
C
      CALL FULL
      return
      end
c
c
c
      subroutine drawframe
      implicit none
c
c     Draws plot frame work
c
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      call thick(1)
      call ctrmag(12)
      CALL LINCOL (3)
      CALL POSITN (0.1, 0.1)
      CALL   JOIN (0.1, 0.9)
      CALL   JOIN (0.9, 0.9)
      CALL   JOIN (0.9, 0.1)
      CALL   JOIN (0.1, 0.1)
      CALL POSITN (0.93, 0.1)
      CALL   JOIN (0.93, 0.9)
      CALL   JOIN (1.35, 0.9)
      CALL   JOIN (1.35, 0.1)
      CALL   JOIN (0.93, 0.1)
      CALL POSITN (0.93, 0.85)
      CALL   JOIN (1.35, 0.85)
      CALL POSITN (0.93, 0.15)
      CALL   JOIN (1.35, 0.15)
      CALL POSITN (0.93, 0.20)
      CALL   JOIN (1.35, 0.20)
      CALL POSITN (0.93, 0.25)
      CALL   JOIN (1.35, 0.25)
c
c
      return
      end
c
c
c
      subroutine drawtitles(title,table,xlabel,ylabel)
      implicit none
c
c     Draw labels
c
      integer l,lenstr
      external lenstr
      character*(*) title,xlabel,ylabel,table
c

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRMAG (20)
      CALL LINCOL (1)
      CALL THICK  (2)
      L = LENSTR(TITLE)
      CALL PCSCEN (0.8, 0.95, TITLE(:L))
      L = LENSTR(TABLE)
      CALL PCSCEN (1.18, 0.87, TABLE(:L))
      L = LENSTR (XLABEL)
      CALL PCSCEN (0.5,0.025,XLABEL(:L))
      CALL THICK  (1)
      CALL CTRMAG (12)
c
c     Ylabel
c
      CALL PSPACE (0.0, 1.35, 0.11, 0.89)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRORI (90.0)
      CALL LINCOL (1)
      CALL CTRMAG (14)
      CALL THICK  (2)
      L = LENSTR (YLABEL)
      CALL PCSEND (0.02,0.99,YLABEL(:L))
      CALL CTRORI (0.0)
      CALL THICK  (1)
c
      return
      end
c
c
c
      real function fmax(values,pts,ncases)
      implicit none
c
c     Finds maximum in array
c
      include 'solparams'
c
      real values(mxspts,maxcases)
      integer pts,ncases
c
      integer i,j
c
      fmax = -1.0e30
      do i = 1,ncases
         do j = 1,pts
            fmax = amax1(fmax,values(j,i))
         end do
      end do
      return
      end
c
c
c
      real function fmin(values,pts,ncases)
      implicit none
c
c     Finds minimum in array
c
      include 'solparams'
c
      real values(mxspts,maxcases)
      integer pts,ncases
c
      integer i,j
c
      fmin = 1.0e30
      do i = 1,ncases
         do j = 1,pts
            fmin = amin1(fmin,values(j,i))
         end do
      end do
      return
      end
c
c
c
      subroutine loadparms
      implicit none
c
      include 'solrk'
c
c     Loads Cash-Karp parameters
c
      integer i,j
c
c     Initialize
c
      do i = 1,6
         ai(i) = 0.0
         ci(i) = 0.0
         cip(i) = 0.0
         do j = 1,5
            bij(i,j) = 0.0
         end do
      end do
c
      ai(2) = 1.0/5.0
      ai(3) = 3.0/10.0
      ai(4) = 3.0/5.0
      ai(5) = 1.0
      ai(6) = 7.0/8.0
c
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
c
      ci(1) = 37.0/378.0
      ci(3) = 250.0/621.0
      ci(4) = 125.0/594.0
      ci(6) = 512.0/1771.0
c
      cip(1) = 2825.0/27648.0
      cip(3) = 18575.0/48384.0
      cip(4) = 13525.0/55296.0
      cip(5) = 277.0/14336.0
      cip(6) = 1.0/4.0
c
      return
      end
c
c
c
      subroutine calcfluxes(gtarg,ionptarg,elecptarg,e2dgtarg,
     >                      presstarg,gamcor,gamecor,ike2d_start,
     >                      g_pfzsol,pe_pfzsol,pi_pfzsol,pr_pfzsol,
     >                      pfz_dist_opt,pfz_dist_param)
      implicit none
      include 'params'
      include 'solparams'
      include 'solswitch'
c
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
c
      real*8 gtarg(mxspts,3)
      real*8 e2dgtarg(mxspts,3)
      real*8 ionptarg(mxspts,3)
      real*8 elecptarg(mxspts,3)
      real*8 presstarg(mxspts,3)
      real*8 gamcor,gamecor
      integer ike2d_start

      ! record total actual flux to sol and pfz target regions - inner and outer
      real*8 :: solpfz_fluxes(2,2,4) ! note: this code only works properly for single null/non-extended geometries
      ! distribute extra sources to sol from pfz
      real*8 :: g_pfzsol(maxpts,3)
      real*8 :: pe_pfzsol(maxpts,3)
      real*8 :: pi_pfzsol(maxpts,3)
      real*8 :: pr_pfzsol(maxpts,3)
      integer :: pfz_dist_opt
      real*8 :: pfz_dist_param(2)


      

c
c     This subroutine calculates the target partcle and
c     power fluxes based on the target data and the
c     SOL 22 options selected.
c
c
c     Local Variables
c

      real*8 :: dist_fact(maxpts,2)
      real*8 :: maxpress(2)
      integer :: maxpress_ir(2)

      integer ik,ir,in,id,ierr,ringno
      external ringno
      real rfact,v0
      real n1,te1,ti1,vpe2d,v1e2d,te0
      real mach0o,mach0i,rmeano,rmeani
      real gae,gaio,gaii

c
c     Initialization
c
      solpfz_fluxes = 0.0
c
      do ir = irsep,nrs
c
c        Calculate mach numbers -
c
         if (switch(swmach).eq.3.0) then
c
c           Outer
c
            v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
     >             * econv/mconv)
c
            mach0o = kvds(idds(ir,2)) / v0
c
c           Inner
c
            v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
     >             * econv/mconv)
c
            mach0i = kvds(idds(ir,1)) / v0
c
         else
            mach0o = 1.0
            mach0i = 1.0
         endif

c
c        Calculate heat transmission
c
         gae = 5.0 + gamecor
         gaio = 2.5 +  0.5* mach0o**2
     >         * (1.0 + kteds(idds(ir,2))/ktids(idds(ir,2)))
     >         + gamcor
         gaii = 2.5 +  0.5* mach0i**2
     >         * (1.0 + kteds(idds(ir,1))/ktids(idds(ir,1)))
     >         + gamcor
c
c        Calculate target fluxes (both particles and heat)
c        Assume target MACH number = 1.0
c
c        Outer
c
         if (switch(swe2d).eq.0.0.or.switch(swe2d).eq.5) then
c
           if (switch(swmajr).eq.4.0) then
              rmeano = rp(idds(ir,2))
              rmeani = rp(idds(ir,1))
           else
              rmeano = 1.0
              rmeani = 1.0
           endif
c
c           v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
c     >             * econv/mconv)
c
           v0 = -abs(kvds(idds(ir,2)))
c
           gtarg(ir,2)= knds(idds(ir,2)) * v0 * rmeano
c
           if (e2dtargopt.eq.5) then

              elecptarg(ir,2)= e2dpepara(1,ir)
              ionptarg(ir,2) = e2dpipara(1,ir)

           else

              elecptarg(ir,2)=gae*kteds(idds(ir,2))*econv*gtarg(ir,2)
              ionptarg(ir,2)=gaio*ktids(idds(ir,2))*econv*gtarg(ir,2)

           endif

           ! calculate target pressrure
           presstarg(ir,2) = knds(idds(ir,2)) * 
     >                   (kteds(idds(ir,2))+ktids(idds(ir,2))) 
     >                   * (1.0 + mach0o**2) * econv


c
c          Inner
c
c           v0 = -sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb
c     >             * econv/mconv)
c
           v0 = -abs(kvds(idds(ir,1)))

           gtarg(ir,1)= knds(idds(ir,1)) * v0 * rmeani
c
           if (e2dtargopt.eq.5) then

              elecptarg(ir,1)= -e2dpepara(nks(ir)+1,ir)
              ionptarg(ir,1) = -e2dpipara(nks(ir)+1,ir)
c
           else
c
              elecptarg(ir,1)=gae*kteds(idds(ir,1))*econv*gtarg(ir,1)
              ionptarg(ir,1)=gaii*ktids(idds(ir,1))*econv*gtarg(ir,1)
c
           endif

           ! calculate target pressrure
           presstarg(ir,1) = knds(idds(ir,1)) * 
     >                   (kteds(idds(ir,1))+ktids(idds(ir,1))) 
     >                   * (1.0 + mach0i**2) * econv

c
c          Totals
c
           gtarg(ir,3) = gtarg(ir,1) + gtarg(ir,2)
           elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
           ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
           presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
           write (6,'(a,i4,3g18.7)') 'DIVIMP  FLUXES:',ir,
     >                              gtarg(ir,1),
     >                               gtarg(ir,2),gtarg(ir,3)
c
           write (6,'(a,i4,3g18.7)') 'DIVIMP  ELEC P:',ir,
     >                              elecptarg(ir,1),
     >                           elecptarg(ir,2),elecptarg(ir,3)
c
           write (6,'(a,i4,3g18.7)') 'DIVIMP  ION  P:',ir,
     >                              ionptarg(ir,1),
     >                              ionptarg(ir,2),ionptarg(ir,3)
           write (6,'(a,i4,3g18.7)') 'DIVIMP  PRESS :',ir,
     >                              presstarg(ir,1),
     >                              presstarg(ir,2),presstarg(ir,3)
c
           if (fluxpts.gt.0.0) then
c
c            Find ring for EDGE2D data
c
             in = ringno(ir,fluxinfo,fluxpts,maxins,4,ierr)
c
c            Inner target
c
             e2dgtarg(ir,1) = -abs(fluxinfo(in,2))  /
     >                (dds(idds(ir,1))* 2.0 * PI * rp(idds(ir,1)))
c
c            Outer target
c
             e2dgtarg(ir,2) = -abs(fluxinfo(in,3))  /
     >                (dds(idds(ir,2))* 2.0 * PI * rp(idds(ir,2)))
c
             e2dgtarg(ir,3) = e2dgtarg(ir,1)+e2dgtarg(ir,2)
c
             if (cprint.eq.9)
     >          write (6,'(a,i4,3g18.7)') 'E2D RAW FLUXES:',ir,
     >                              e2dgtarg(ir,1),
     >                               e2dgtarg(ir,2),e2dgtarg(ir,3)
c
c            Modify for non-orthogonality and magnetic field
c
             e2dgtarg(ir,1)=e2dgtarg(ir,1)*kbfst(ir,1)
     >                                  /costet(idds(ir,1))
c
             e2dgtarg(ir,2)=e2dgtarg(ir,2)*kbfst(ir,2)
     >                                   /costet(idds(ir,2))
c
             e2dgtarg(ir,3) = e2dgtarg(ir,1)+e2dgtarg(ir,2)
c
c            Print outs
c
             if (cprint.eq.9) then
                write (6,'(a,i4,3g18.7)') 'E2D COR FLUXES:',ir,
     >                              e2dgtarg(ir,1),
     >                               e2dgtarg(ir,2),e2dgtarg(ir,3)
                write (6,'(a,i4,3g18.7)') 'RATIOS        :',ir,
     >                         e2dgtarg(ir,1)/gtarg(ir,1),
     >          e2dgtarg(ir,2)/gtarg(ir,2),e2dgtarg(ir,3)/gtarg(ir,3)
                write (6,*)
c
              write (6,'(a,i4,3g18.7)') 'E2D TARG FLUXES-E2D   :',ir,
     >                              e2dtarg(ir,6,1),
     >                               e2dtarg(ir,6,2)
              write (6,'(a,i4,3g18.7)') 'E2D TARG FLUXES-OLDDIV:',ir,
     >                              e2dtarg(ir,5,1),
     >                               e2dtarg(ir,5,2)
              endif

           endif
c
c
         elseif (switch(swe2d).eq.1.0.or.switch(swe2d).eq.6.0) then
c
           if (switch(swmajr).eq.4.0) then
              rmeano = rs(1,ir)
              rmeani = rs(nks(ir),ir)
           else
              rmeano = 1.0
              rmeani = 1.0
           endif
c
            n1 = cellvals(ir,1,2)
            te0= kteds(idds(ir,2))
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
            v1e2d = cellvals(ir,4,2)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,2)= n1 * vpe2d * rmeano
            elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
            ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)
            presstarg(ir,2) = n1 * (te1+ti1) * econv 
     >                        + n1 * vpe2d**2 *crmb * mconv
c
            n1 = cellvals(ir,1,1)
            te0= kteds(idds(ir,1))
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
            v1e2d = -cellvals(ir,4,1)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,1)= n1 * vpe2d * rmeani
            elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
            ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)
            presstarg(ir,1) = n1 * (te1+ti1) * econv 
     >                        + n1 * vpe2d**2 *crmb * mconv
c
            gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
            elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
            ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
            presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
         elseif (switch(swe2d).eq.2.0.or.switch(swe2d).eq.4.0) then
c
           if (switch(swmajr).eq.4.0) then
              rmeano = rs(1,ir)
              rmeani = rs(nks(ir),ir)
           else
              rmeano = 1.0
              rmeani = 1.0
           endif
c
            n1 = cellvals(ir,1,2)
            te0= kteds(idds(ir,2))
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
            v1e2d = cellvals(ir,4,2)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,2)= n1 * v1e2d * rmeano
            elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
            ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)
            presstarg(ir,2) = n1 * (te1+ti1) * econv 
     >                        + n1 * v1e2d**2 *crmb * mconv
c
            n1 = cellvals(ir,1,1)
            te0= kteds(idds(ir,1))
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
            v1e2d = -cellvals(ir,4,1)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,1)= n1 * v1e2d * rmeani
            elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
            ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)
            presstarg(ir,1) = n1 * (te1+ti1) * econv 
     >                        + n1 * v1e2d**2 *crmb * mconv
c
            gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
            elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
            ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
            presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
         elseif (switch(swe2d).eq.3.0.or.switch(swe2d).eq.7.0) then
c
           if (switch(swmajr).eq.4.0) then
              rmeano = rs(1,ir)
              rmeani = rs(nks(ir),ir)
           else
              rmeano = 1.0
              rmeani = 1.0
           endif
c
            n1 = cellvals(ir,1,2)
            te0= kteds(idds(ir,2))
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
            v1e2d = cellvals(ir,4,2)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,2)= -0.5 * (abs(e2dflux(1,ir))
     >                   + abs(e2dflux(2,ir))) * rmeano
            elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
            ionptarg(ir,2)=gaio*ti1*econv*gtarg(ir,2)
            presstarg(ir,2) = n1 * (te1+ti1) * econv 
     >                        + n1 * v1e2d**2 *crmb * mconv
c
c            write (6,*) 'E2d1:',ir,nks(ir),e2dflux(1,ir),
c     >                            e2dflux(2,ir),
c     >                            gtarg(ir,2)
c
c
            n1 = cellvals(ir,1,1)
            te0= kteds(idds(ir,1))
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
            v1e2d = -cellvals(ir,4,1)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
            gtarg(ir,1)= -0.5 * (abs(e2dflux(nks(ir),ir))
     >                        + abs(e2dflux(nks(ir)+1,ir)))
     >                         * rmeani
            elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
            ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)
            presstarg(ir,1) = n1 * (te1+ti1) * econv 
     >                        + n1 * v1e2d**2 *crmb * mconv
c
c            write (6,*) 'E2d1:',ir,nks(ir),e2dflux(nks(ir),ir),
c     >                            e2dflux(nks(ir)+1,ir),
c     >                            gtarg(ir,1)
c
c
c
            gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
            elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
            ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)
            presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
         elseif (switch(swe2d).eq.8.0) then
c
            if (switch(swmajr).eq.4.0) then
               rmeano = rs(1,ir)
               rmeani = rs(nks(ir),ir)
            else
               rmeano = 1.0
               rmeani = 1.0
            endif
c
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
c
            gtarg(ir,2)= -abs((e2dgpara(1,ir)+e2dgpara(2,ir))/2.0)
c
            if (e2dtargopt.eq.5) then

              elecptarg(ir,2)= (e2dpepara(1,ir)+e2dpepara(2,ir))/2.0
              ionptarg(ir,2) = (e2dpipara(1,ir)+e2dpipara(2,ir))/2.0

            else

              elecptarg(ir,2)=gae*te1*econv*gtarg(ir,2)
              ionptarg(ir,2) =gaio*ti1*econv*gtarg(ir,2)

            endif
c
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
c
            gtarg(ir,1)= -abs((e2dgpara(nks(ir)+1,ir)
     >                   +e2dgpara(nks(ir),ir))/2.0)
c
            if (e2dtargopt.eq.5) then
c
              elecptarg(ir,2)= -(e2dpepara(nks(ir),ir)
     >                         +e2dpepara(nks(ir)+1,ir))/2.0
              ionptarg(ir,2) = -(e2dpipara(nks(ir),ir)
     >                         +e2dpipara(nks(ir)+1,ir))/2.0
c
            else
c
               elecptarg(ir,1)=gae*te1*econv*gtarg(ir,1)
               ionptarg(ir,1) =gaii*ti1*econv*gtarg(ir,1)
c
            endif
c
            gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
            elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
            ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

            ! jdemod - turn off target pressure compensation for these options
            presstarg(ir,1) = 0.0
            presstarg(ir,2) = 0.0
            presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
         elseif (switch(swe2d).eq.9.0) then
c
            if (switch(swmajr).eq.4.0) then
               rmeano = rs(1,ir)
               rmeani = rs(nks(ir),ir)
            else
               rmeano = 1.0
               rmeani = 1.0
            endif
c
            gtarg(ir,2)=      (e2dgpara(ike2d_start,ir)
     >                        +e2dgpara(ike2d_start+1,ir))/2.0
c
            elecptarg(ir,2)=  (e2dpepara(ike2d_start,ir)
     >                        +e2dpepara(ike2d_start+1,ir))/2.0
c
            ionptarg(ir,2) =  (e2dpipara(ike2d_start,ir)
     >                        +e2dpipara(ike2d_start+1,ir))/2.0
c
            gtarg(ir,1)= -(e2dgpara(nks(ir)-ike2d_start+1,ir)
     >               +e2dgpara(nks(ir)-ike2d_start+1+1,ir))/2.0
c
            elecptarg(ir,1)=-(e2dpepara(nks(ir)-ike2d_start+1,ir)
     >               +e2dpepara(nks(ir)-ike2d_start+1+1,ir))/2.0
c
            ionptarg(ir,1) =-(e2dpipara(nks(ir)-ike2d_start+1,ir)
     >               +e2dpipara(nks(ir)-ike2d_start+1+1,ir))/2.0
c
            gtarg(ir,3) = gtarg(ir,2) + gtarg(ir,1)
            elecptarg(ir,3) = elecptarg(ir,1) + elecptarg(ir,2)
            ionptarg(ir,3) = ionptarg(ir,1) + ionptarg(ir,2)

            ! jdemod - turn off target pressure compensation for these options
            presstarg(ir,1) = 0.0
            presstarg(ir,2) = 0.0
            presstarg(ir,3) = presstarg(ir,1)+presstarg(ir,2) 
c
         endif
c
         write(6,'(a,i8,10(1x,g12.5))') 'Calcfluxes:',ir,
     >      presstarg(ir,1),presstarg(ir,2),presstarg(ir,3),
     >      elecptarg(ir,1),elecptarg(ir,2),elecptarg(ir,3),
     >      ionptarg(ir,1),ionptarg(ir,2),ionptarg(ir,3)

      end do

      !
      ! Accumulate data for main SOL and PFZ
      !


      do id = 1,2

          do ir = irsep,nrs

            if (ir.eq.irwall.or.ir.eq.irtrap) cycle

            if (ir.le.irwall) then 
               in=1
            else
               in=2
            endif

            solpfz_fluxes(in,id,1)=solpfz_fluxes(in,id,1) 
     >              + gtarg(ir,id) * dds(idds(ir,id))
            solpfz_fluxes(in,id,2) = solpfz_fluxes(in,id,2) 
     >              + elecptarg(ir,id)  * dds(idds(ir,id))
            solpfz_fluxes(in,id,3) = solpfz_fluxes(in,id,3) 
     >              +ionptarg(ir,id)   * dds(idds(ir,id))
            solpfz_fluxes(in,id,4) = solpfz_fluxes(in,id,4) 
     >              +presstarg(ir,id)  * dds(idds(ir,id))

         end do

         in = 1
             write(6,'(a,2i8,6(1x,g12.5))') 'SOL FLUX:',in,id,
     >           solpfz_fluxes(in,id,1),solpfz_fluxes(in,id,2),
     >           solpfz_fluxes(in,id,3),solpfz_fluxes(in,id,4)

          in = 2
             write(6,'(a,2i8,6(1x,g12.5))') 'PFZ FLUX:',in,id,
     >           solpfz_fluxes(in,id,1),solpfz_fluxes(in,id,2),
     >           solpfz_fluxes(in,id,3),solpfz_fluxes(in,id,4)


      end do

         !
         ! Calculate assignment of pfz fluxes to rings
         ! - use some decay profile ... linear or exponential? or something else
         !
         ! May not use them all but need to be calculated where the target
         ! geometry information is still available
         ! 
         ! sepdist2 contains the distance along the target to the target element 
         !
         

      ! find maximum pressure inner and outer
      maxpress=0.0
      maxpress_ir = 0

      do ir = irsep,nrs
         do id = 1,2
            if (presstarg(ir,id).gt.maxpress(id)) then 
               maxpress(id) = presstarg(ir,id)
               maxpress_ir(id) = ir
            endif
         end do
      end do

      dist_fact = 0

      ! flat distribution over distance pfz_dist_param
      if (pfz_dist_opt.eq.1) then 

         do ir = irsep,irwall
            do id = 1,2
              if (sepdist2(idds(ir,id)).le.pfz_dist_param(id)) then 
                 dist_fact (ir,id) = dds(idds(ir,id))
                 dist_fact(irwall+1,id) = dist_fact(irwall+1,id) 
     >                           + dist_fact(ir,id)
              endif
           end do
        end do

      ! linear decay distribution over distance pfz_dist_param (largest at separatrix)
      elseif (pfz_dist_opt.eq.2) then
      

         do ir = irsep,irwall
            do id = 1,2
              if (sepdist2(idds(ir,id)).le.pfz_dist_param(id)) then 
                 dist_fact (ir,id) = dds(idds(ir,id)) * 
     >             (pfz_dist_param(id)-sepdist2(idds(ir,id)))
     >                           /pfz_dist_param(id)
                 dist_fact(irwall+1,id) = dist_fact(irwall+1,id) 
     >                           + dist_fact(ir,id)
              endif
           end do
        end do



      ! exponential decay distribution with lambda = pfz_dist_param (largest at separatrix)
      elseif (pfz_dist_opt.eq.3) then

         do ir = irsep,irwall
            do id = 1,2
                 dist_fact (ir,id) = dds(idds(ir,id)) 
     >                 * exp(-sepdist2(idds(ir,id))/pfz_dist_param(id))
                 dist_fact(irwall+1,id) = dist_fact(irwall+1,id) 
     >                           + dist_fact(ir,id)
           end do
        end do


         
      ! linear decay to peak of pressure profile on target
      elseif (pfz_dist_opt.eq.4) then

         do ir = irsep,irwall
            do id = 1,2
              if (ir.le.maxpress_ir(id)) then 
                 dist_fact (ir,id) = dds(idds(ir,id)) * 
     >             (sepdist2(idds(maxpress_ir(id),id))
     >                       -sepdist2(idds(ir,id)))
     >                       /sepdist2(idds(maxpress_ir(id),id))
                 dist_fact(irwall+1,id) = dist_fact(irwall+1,id) 
     >                           + dist_fact(ir,id)
              endif
           end do
        end do


      ! exponential decay to peak of pressure profile on target
      elseif (pfz_dist_opt.eq.5) then

         do ir = irsep,irwall
            do id = 1,2
              if (ir.le.maxpress_ir(id)) then 
                 dist_fact (ir,id) = dds(idds(ir,id)) 
     >                    * exp(-sepdist2(idds(ir,id))
     >                       /sepdist2(idds(maxpress_ir(id),id)))
                 dist_fact(irwall+1,id) = dist_fact(irwall+1,id) 
     >                           + dist_fact(ir,id)
              endif
           end do
        end do


      endif


      ! normalize distribution 
      do ir = irsep,irwall
         do id = 1,2

           if (dist_fact(irwall+1,id).gt.0.0) then   
              dist_fact(ir,id) = dist_fact(ir,id)/dist_fact(irwall+1,id)
           else
              dist_fact(ir,id) = 0.0
           endif

         end do 
      end do 

      ! distribute pfz fluxes
      do ir = irsep,irwall

         do id = 1,2

            if (dds(idds(ir,id)).gt.0.0) then
              g_pfzsol(ir,id) = solpfz_fluxes(2,id,1)
     >                  *dist_fact(ir,id)/dds(idds(ir,id))
              pe_pfzsol(ir,id) = solpfz_fluxes(2,id,2)
     >                  *dist_fact(ir,id)/dds(idds(ir,id))
              pi_pfzsol(ir,id) = solpfz_fluxes(2,id,3)
     >                  *dist_fact(ir,id)/dds(idds(ir,id))
              pr_pfzsol(ir,id) = solpfz_fluxes(2,id,3)
     >                  *dist_fact(ir,id)/dds(idds(ir,id))
           else
              g_pfzsol(ir,id) = 0.0
              pe_pfzsol(ir,id) = 0.0
              pi_pfzsol(ir,id) = 0.0
              pr_pfzsol(ir,id) = 0.0
          endif

       write(6,'(a,2i8,6(1x,g12.5))') 'ADDITIONAL FLUX:',ir,id,
     >     g_pfzsol(ir,id),pe_pfzsol(ir,id), 
     >     pi_pfzsol(ir,id),pr_pfzsol(ir,id) 

         end do 
         
         

       end do


!      lambda = 0.01
!      dist_opt = 1


      ! Calculate distribution factors

      !do ir = irsep,irwall

       !  if (dist_opt.eq.1) then     ! flat

!         elseif (dist_opt.eq.2) then   ! linear decay

 !    elseif (dist_opt.eq.3) then ! exponential decay

  !       endif

         ! normalize the distribution


   !      dist_fact(ir,1) 



c
      return
      end
c
c
c
      subroutine calcsoliz(rconst,recfrac,gtarg,areasum,gperpa,
     >                     oldknbs,grad_oldknbs,ioniz_src,rec_src,
     >                     gperprat,ike2d_start)
      implicit none
      include 'params'
      include 'solparams'
      include 'solswitch'
      include 'printopt' 
c      include 'solcommon'
      real    gperpa(maxnks,maxnrs),oldknbs(maxnks,maxnrs)
      real    grad_oldknbs(maxnks,maxnrs)
      real    ioniz_src(maxnks,maxnrs),rec_src(maxnks,maxnrs)
      real*8 rconst(mxspts),gperprat(mxspts)
      real*8 gtarg(mxspts,3),areasum(mxspts)
      real*8 gamcor,gamecor,recfrac
      integer ike2d_start
c
c
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
c
c     CALCSOLIZ:
c
c     This routine calculates the total
c     ionization on each ring in order
c     to evaluate the net total excess or deficit of flux
c     onto the ring. These values for each ring are then
c     returned and used in the solver if gperp option 2 has
c     been specified.
c
c
c
c     Local variables
c
      integer ik,ir,in,id,startik,endik
      real rfact,v0,rmean1,rmean2
      real srcinteg (maxnks+3,maxnrs)
      real recinteg (maxnks+3,maxnrs)
      real srcsum,recsum,nesum,gradsum
      real gradsumneg,gradsumpos
      real n1,te1,ti1,vpe2d,v1e2d,te0
c
c     New variables for summaries
c
      real intgrad,intgradpos,intgradneg,inte2dsrc,intdist
      real inte2dhrec,intgrad2,ds,dcell,mulfact
      integer ikin,irin,ikout,irout
      real    influx,outflux,netflux,intnflux,flux1,flux2,flux3
      real    flux4
      real fluxst,fluxend,ioniz,rec,gperpd,gperpn,gdiv,sider
      real intgperpd,intgnet,intgdiv
      real initflux,dp,dp1,dp2,brat1,brat2,startflux
      real flux_const,endflux
c
      logical debug
c
c     Recylcling flux factor
c
      rfact = recfrac
      debug = .false.
c
c     Calculate the net fluxes on each ring after summing over
c     the ionization source.
c
c
      do ir = irsep,nrs
c
c
c        Now to sum up ionization and recombination (if on).
c
c        For this specific case ONLY the integral will
c        be done simply on a cell by cell basis - and
c        so the R-value used can be the mean between the
c        two end-points.
c
c
         srcsum = 0.0
         recsum = 0.0
         areasum(ir) = 0.0
         nesum = 0.0
         gradsum = 0.0
         gradsumpos = 0.0
         gradsumneg = 0.0
c
c        Set IK loop limits based on E2D option
c
         if (switch(swe2d).eq.0.0) then
            startik = 1
            endik   = nks(ir)
         elseif (switch(swe2d).gt.0.0) then
            startik = ike2d_start
            endik   = nks(ir) - ike2d_start + 1
         endif
c
c         write (6,*) 'E2d switch:',switch(swe2d),ike2d_start,
c     >                startik,endik
c

c
         do ik = startik,endik
c
            if (switch(swmajr).eq.4.0) then
               rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
               rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
            else
               rmean1 = 1.0
               rmean2 = 1.0
            endif
c
            if (.not.(ik.eq.startik.and.
     >         (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >          switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >          switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >          switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
               srcsum = srcsum
     >            + ioniz_src(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1
               recsum = recsum
     >            + rec_src(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1
               areasum(ir) = areasum(ir)
     >                     + (kss2(ik,ir)-ksb(ik-1,ir)) * rmean1
c
               nesum = nesum
     >            + oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))*rmean1
               gradsum = gradsum
     >                   + grad_oldknbs(ik,ir)
     >                   * (kss2(ik,ir)-ksb(ik-1,ir))*rmean1
c
               if (grad_oldknbs(ik,ir).lt.0.0) then
                  gradsumneg = gradsumneg
     >                   + grad_oldknbs(ik,ir)
     >                   * (kss2(ik,ir)-ksb(ik-1,ir))*rmean1
               else
                  gradsumpos = gradsumpos
     >                   + grad_oldknbs(ik,ir)
     >                   * (kss2(ik,ir)-ksb(ik-1,ir))*rmean1
               endif
c
            endif
c
c
c           Do in two parts to get integral at cell centre.
c
            srcinteg(ik,ir) = srcsum
            recinteg(ik,ir) = recsum
c
            if (.not.(ik.eq.endik.and.
     >         (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >          switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >          switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >          switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
               srcsum = srcsum
     >            + ioniz_src(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2
               recsum = recsum
     >            + rec_src(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2
               areasum(ir) = areasum(ir)
     >                     + (ksb(ik,ir)-kss2(ik,ir)) * rmean2
c
               nesum = nesum
     >            + oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))*rmean2
c
               gradsum = gradsum
     >                   + grad_oldknbs(ik,ir)
     >                   * (ksb(ik,ir)-kss2(ik,ir))*rmean2
c
               if (grad_oldknbs(ik,ir).lt.0.0) then
                  gradsumneg = gradsumneg
     >                   + grad_oldknbs(ik,ir)
     >                   * (ksb(ik,ir)-kss2(ik,ir))*rmean2
               else
                  gradsumpos = gradsumpos
     >                   + grad_oldknbs(ik,ir)
     >                   * (ksb(ik,ir)-kss2(ik,ir))*rmean2
               endif
c
            endif
c
            if ( ( (ksmaxs2(ir)/2.0).ge.ksb(ik-1,ir) )
     >         .and.((ksmaxs2(ir)/2.0).lt.ksb(ik,ir))) then
c
c              Assign amount of outer ionization source
c
               srcinteg(nks(ir)+2,ir) = srcsum
               recinteg(nks(ir)+2,ir) = recsum
            endif
c
         end do
c
         srcinteg(nks(ir)+1,ir) = srcsum
         recinteg(nks(ir)+1,ir) = recsum
c
         srcinteg(nks(ir)+3,ir) = srcsum - srcinteg(nks(ir)+2,ir)
         recinteg(nks(ir)+3,ir) = recsum - recinteg(nks(ir)+2,ir)
c
c        Calculate Rconst
c
         if (switch(swrecom).eq.0.0) then
            rconst(ir)= - (gtarg(ir,3)+rfact*srcinteg(nks(ir)+1,ir))
     >                 /areasum(ir)
         elseif (switch(swrecom).eq.1.0.or.switch(swrecom).eq.2) then
            rconst(ir)= - (gtarg(ir,3)+rfact*srcinteg(nks(ir)+1,ir)
     >                      - recinteg(nks(ir)+1,ir))
     >                 /areasum(ir)
         endif
c
c        Calculate Gperp fraction ratios - for use in Gperp option 7.
c
         if (switch(swgperp).eq.7.0.or.switch(swgperpp).eq.7.0) then
c
            if (gradsum.eq.0.0) then
               gperprat(ir)= 100.0
            else
               gperprat(ir)= max(abs(gradsumpos/gradsum),
     >                           abs(gradsumneg/gradsum))
            endif
c
         endif
c
c         write (6,'(a,i4,1p,10(1x,g12.5))') 'Rconst:',ir,rconst(ir),
c     >                        areasum(ir),
c     >         gradsum,gradsumpos,gradsumneg,gradsumpos/gradsum,
c     >                    gradsumneg/gradsum
c
c        Calculate a distributed gperp source if that is
c        required.
c
c
c        Distributed with d2ne/dr2 component and constant.
c
         if ((switch(swgperp).eq.8.and.ir.le.irwall).or.
     >       (switch(swgperpp).eq.8.and.ir.ge.irtrap)) then
c
            if (cprint.eq.9)
     >         write (6,'(a,i4,5g13.5)') 'GPERPA:',ir,nesum,gradsum,
     >                  rconst(ir),areasum(ir),rconst(ir)*areasum(ir)
c
c           Modify RCONST for GPERP term at given DPERP
c
            rconst(ir) = - (-rconst(ir)*areasum(ir)+cdperp*gradsum)
     >                     /areasum(ir)
c
            intgrad = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               gperpa(ik,ir) = cdperp * grad_oldknbs(ik,ir)
     >                         + rconst(ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                  endif
c
               endif
c
               if (cprint.eq.9)
     >           write(6,'(2i4,7(1x,g13.5))') ir,ik,gperpa(ik,ir),
     >                 oldknbs(ik,ir),cdperp*grad_oldknbs(ik,ir),
     >                 intgrad,
     >                 intgradpos,intgradneg
c
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                  endif
c
               endif

            end do


            if (cprint.eq.9)
     >         write (6,'(a,i4,7g13.5)') 'GPERPA-END:',ir,gradsum,
     >                  rconst(ir),areasum(ir),rconst(ir)*areasum(ir),
     >                  intgrad,intgradpos,intgradneg
c
         endif

c
c

c
c        Distributed proportional to ne or d2ne/dr2
c
         if ( ((switch(swgperp).eq.3.0.or.switch(swgperp).eq.7.0)
     >         .and.ir.le.irwall).or.
     >        ((switch(swgperpp).eq.3.or.switch(swgperpp).eq.7.0)
     >          .and.ir.gt.irwall)) then
c
            if (cprint.eq.9)
     >         write (6,'(a,i4,5g13.5)') 'GPERPA:',ir,nesum,gradsum,
     >                  rconst(ir),areasum(ir),rconst(ir)*areasum(ir)
c
            intgrad = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               if (switch(swgperp).eq.3.0.or.
     >             switch(swgperpp).eq.3.0)    then
c
                  gperpa(ik,ir)=oldknbs(ik,ir)
     >                          * rconst(ir)*areasum(ir)
     >                          / nesum
c
               elseif (switch(swgperp).eq.7.0.or.
     >                 switch(swgperpp).eq.7.0)   then
c
                  gperpa(ik,ir)=grad_oldknbs(ik,ir)
     >                          * rconst(ir)*areasum(ir)
     >                          / gradsum
c
               else
c
                  gperpa(ik,ir)=rconst(ir)
c
               endif
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))

                  endif
c
               endif
c
               if (cprint.eq.9)
     >            write(6,'(2i4,7(1x,g13.5))') ir,ik,gperpa(ik,ir),
     >                 oldknbs(ik,ir),grad_oldknbs(ik,ir),intgrad,
     >                 intgradpos,intgradneg
c
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))

                  endif
c
               endif

            end do


            if (cprint.eq.9)
     >         write (6,'(a,i4,7g13.5)') 'GPERPA-END:',ir,gradsum,
     >                  rconst(ir),areasum(ir),rconst(ir)*areasum(ir),
     >                  intgrad,intgradpos,intgradneg


         endif
c
c     ------------------------------------------------------------------
c
c
c        The following are a series of DIAGNOSTIC print outs related
c        mostly to the EDGE2D sources that are read in from the GHOST file.
c
c        They are only printed for a complete print out situation.
c
c
         if (cprint.eq.10) then
c
            call print_edge2d_flux_analysis(rconst,gtarg,areasum,ir,
     >               gperpa,oldknbs,grad_oldknbs,srcinteg,recinteg,
     >               ike2d_start)
         endif

c     ------------------------------------------------------------------
c
c        Close IR loop
c
      end do
c
c
c
c
c     Print out the results
c
c
      write (6,*) 'Ring Ionization and Flux Integration:'
      write (6,*) 'Recycling factor (rfact) = ',rfact
      write (6,*)
      if (switch(swrecom).eq.0.0) then
           write (6,*) 'Recombination is OFF'
      elseif (switch(swrecom).eq.1.0.or.switch(swrecom).eq.2) then
           write (6,*) 'Recombination is ON'
      endif
      write (6,*)
      write (6,100) outer,inner,outer,inner,outer,inner
      write (6,*)
      do ir = irsep,nrs
         write (6,200) ir,gtarg(ir,3),srcinteg(nks(ir)+1,ir),
     >                 rfact*srcinteg(nks(ir)+1,ir),rconst(ir),
     >                 areasum(ir),srcinteg(nks(ir)+2,ir),
     >                 srcinteg(nks(ir)+3,ir),gtarg(ir,2),
     >                 gtarg(ir,1),recinteg(nks(ir)+1,ir),
     >                 recinteg(nks(ir)+2,ir),
     >                 recinteg(nks(ir)+3,ir)
      end do

c
c     If the debug option is compiled ON then calculate
c     all the potential fluxes - all of the integrals and
c     their respective balances and corresponding Rconst
c     values for BOTH the Edge2d ionization AND the
c     DIVIMP calculated Pinion.
c
      if (debug) then

         call debugsoliz(rfact)

      endif
c
c     Format statements
c

 100  format(10x,'Total Flux',4x,'Total Ioniz',
     >       4x,'Rfact*Ioniz',5x,'Net Rconst',4x,'R(s) Integral',
     >       3x,a5,' Ioniz',3x,a5,' Ioniz',4x,a5,' Flux',
     >       4x,a5,' Flux',4x,'Total Recom',4x, a5,' Recom',
     >       4x, a5,' Recom' )
 200  format(i4,2x,12(e14.6,1x))
c
c
c
      return
      end
c
c
c
      subroutine debugsoliz(rfact)
      implicit none
      include 'params'
      include 'solparams'
      include 'solswitch'
      real rfact
c
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
c
c     This routine analyses the Flux/ionization balance
c     for the BG plasma using various methods of calculating
c     the target fluxes and ionization integrals including
c     ds integration and major radius corrected integration.
c
c
c
c
c     Local variables
c
      integer ik,ir,in,id
      real v0,rmean
      real srcinteg (maxnks+3,maxnrs)
      real srcsum
      real n1,te1,ti1,vpe2d,v1e2d,te0
c
      real*8 rconst(mxspts)
      real*8 gtarg(mxspts,3),areasum(mxspts)
      real gammas(maxnrs,8,2,3)
      real intsrc(maxnks+1,maxnrs,16)
      real srccomp(maxnks,maxnrs,24)
      real intarea(maxnks+1,maxnrs,16)
      real newrconst(maxnrs,16)
c
c      real fluxes(maxnks+1,maxnrs,16)
c
         real srcsume,srcsumre,srcsumr,srchalfe,srchalfre,
     >        srchalf,srchalfr,rmean1,rmean2,ds1,ds2,ds,
     >        ds1h,ds2h,dsh,areasuma,areasumr,areahalf,
     >        areahalfr,srcsuma,srchalfa,da,da1,da2,dah1,
     >        dah2,srcsumra,srchalfra,asuma,asumra,ahalfa,
     >        ahalfra,dah,e2dsrc,fluxtot,fluxout,fluxin,
     >        fluxtotr,fluxinr,fluxoutr,e2dsrcr
c
c
      e2dsrc = 0.0
      fluxtot = 0.0
c
      do ir = 1,nrs
         if (ir.ge.irsep) then
           fluxout =
     >       -e2dflux(1,ir) * dds2(idds(ir,2))
     >       * costet(idds(ir,2)) / kbfst(ir,2)
           fluxin =
     >       e2dflux(nks(ir)+1,ir) * dds2(idds(ir,1))
     >       * costet(idds(ir,1)) / kbfst(ir,1)


           fluxtot = fluxtot + fluxout + fluxin
c
           fluxoutr = -e2dflux(1,ir)
     >       * rp(idds(ir,2)) * dds2(idds(ir,2))
     >       * costet(idds(ir,2)) / kbfst(ir,2)
           fluxinr = e2dflux(nks(ir)+1,ir)
     >       * rp(idds(ir,1)) * dds2(idds(ir,1))
     >       * costet(idds(ir,1)) / kbfst(ir,1)
c

           fluxtotr = fluxtotr + fluxoutr + fluxinr

            write(6,*) 'E2d Flux: ',ir,e2dflux(1,ir),
     >           e2dbvel(1,ir),e2dnbs(1,ir),
     >           e2dflux(nks(ir)+1,ir),
     >           e2dbvel(nks(ir)+1,ir),
     >           e2dnbs(nks(ir),ir),dds2(idds(ir,2)),
     >           dds2(idds(ir,1))

         endif
         do ik = 1,nks(ir)
            e2dsrc = e2dsrc + e2dion(ik,ir) * karea2(ik,ir)
            e2dsrcr = e2dsrcr + e2dion(ik,ir) * karea2(ik,ir)
     >                * rs(ik,ir)
         end do
      end do
c
c
c
      write (6,*) 'EDGE2D: Total Target Flux    = ',fluxtot
      write (6,*) 'EDGE2D: Total Ionization     = ',e2dsrc
      write (6,*) 'EDGE2D: Ratio Ioniz/Flux     = ',e2dsrc/fluxtot
c
      write (6,*) 'EDGE2D: Total Target Flux (R)= ',fluxtotr
      write (6,*) 'EDGE2D: Total Ionization  (R)= ',e2dsrcr
      write (6,*) 'EDGE2D: Ratio Ioniz/Flux  (R)= ',
     >              e2dsrcr/fluxtotr
c
      do ir = irsep,nrs
c
c        Calculate target fluxes
c
c        Outer
c
c
c        DIVIMP  G0
c

           v0 = -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb
     >             * econv/mconv)
           gammas(ir,1,1,2)= knds(idds(ir,2)) * v0
           gammas(ir,1,2,2)= knds(idds(ir,2)) * v0 *rp(idds(ir,2))
c
c          Inner
c
           v0 = -sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb
     >             * econv/mconv)
           gammas(ir,1,1,1)= knds(idds(ir,1)) * v0
           gammas(ir,1,2,1)= knds(idds(ir,1)) * v0 *rp(idds(ir,1))
c
           gammas(ir,1,1,3)= gammas(ir,1,1,1)+gammas(ir,1,1,2)
           gammas(ir,1,2,3)= gammas(ir,1,2,1)+gammas(ir,1,2,2)
c
c        EDGE2D G0
c
c          Outer
c
           gammas(ir,2,1,2)= e2dflux(1,ir)
           gammas(ir,2,2,2)= e2dflux(1,ir) * rp(idds(ir,2))
c
c          Inner
c
           gammas(ir,2,1,1)= -e2dflux(nks(ir)+1,ir)
           gammas(ir,2,2,1)= -e2dflux(nks(ir)+1,ir) * rp(idds(ir,1))
c
           gammas(ir,2,1,3)= gammas(ir,2,1,1)+gammas(ir,2,1,2)
           gammas(ir,2,2,3)= gammas(ir,2,2,1)+gammas(ir,2,2,2)
c
c        EDGE2D G1
c
c
c          Outer
c
           gammas(ir,3,1,2)= e2dflux(2,ir)
           gammas(ir,3,2,2)= e2dflux(2,ir) * krb(1,ir)
c
c          Inner
c
           gammas(ir,3,1,1)= -e2dflux(nks(ir),ir)
           gammas(ir,3,2,1)= -e2dflux(nks(ir),ir) * krb(nks(ir)-1,ir)
c
           gammas(ir,3,1,3)= gammas(ir,3,1,1)+gammas(ir,3,1,2)
           gammas(ir,3,2,3)= gammas(ir,3,2,1)+gammas(ir,3,2,2)
c
c        EGDE2D G1/2 - A
c
c          Outer
c
           gammas(ir,4,1,2)= 0.5*(gammas(ir,2,1,2)+gammas(ir,3,1,2))
           gammas(ir,4,2,2)= 0.5*(gammas(ir,2,2,2)+gammas(ir,3,2,2))
c
c          Inner
c
           gammas(ir,4,1,1)= 0.5*(gammas(ir,2,1,1)+gammas(ir,3,1,1))
           gammas(ir,4,2,1)= 0.5*(gammas(ir,2,2,1)+gammas(ir,3,2,1))
c
           gammas(ir,4,1,3)= gammas(ir,4,1,1)+gammas(ir,4,1,2)
           gammas(ir,4,2,3)= gammas(ir,4,2,1)+gammas(ir,4,2,2)
c
c        EDGE2D G1/2 - B
c
            n1 = cellvals(ir,1,2)
            te0= kteds(idds(ir,2))
            te1= cellvals(ir,2,2)
            ti1= cellvals(ir,3,2)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
           gammas(ir,5,1,2)= n1 * vpe2d
           gammas(ir,5,2,2)= n1 * vpe2d * rs(1,ir)
c
            n1 = cellvals(ir,1,1)
            te0= kteds(idds(ir,1))
            te1= cellvals(ir,2,1)
            ti1= cellvals(ir,3,1)
            vpe2d = -sqrt((2.0*te0-te1+ti1)*econv/(mconv*crmb))
c
           gammas(ir,5,1,1)= n1 * vpe2d
           gammas(ir,5,2,1)= n1 * vpe2d * rs(nks(ir),ir)
c
           gammas(ir,5,1,3)= gammas(ir,5,1,1)+gammas(ir,5,1,2)
           gammas(ir,5,2,3)= gammas(ir,5,2,1)+gammas(ir,5,2,2)
c
c
c        EDGE2D G1/2 - C
c
            n1 = cellvals(ir,1,2)
            v1e2d = cellvals(ir,4,2)
c
           gammas(ir,6,1,2)= n1 * v1e2d
           gammas(ir,6,2,2)= n1 * v1e2d * rs(1,ir)
c
            n1 = cellvals(ir,1,1)
            v1e2d = -cellvals(ir,4,1)
c
           gammas(ir,6,1,1)= n1 * v1e2d
           gammas(ir,6,2,1)= n1 * v1e2d * rs(nks(ir),ir)
c
           gammas(ir,6,1,3)= gammas(ir,6,1,1)+gammas(ir,6,1,2)
           gammas(ir,6,2,3)= gammas(ir,6,2,1)+gammas(ir,6,2,2)
c
c
c
c        EDGE2D G0 - KAREA
c
c          Outer
c
           id = idds(ir,2)
           gammas(ir,7,1,2)= e2dflux(1,ir) * dds2(id) * costet(id)
     >                       / kbfst(ir,2)
           gammas(ir,7,2,2)= e2dflux(1,ir) * rp(idds(ir,2))
     >                 * dds2(id) * costet(id) / kbfst(ir,2)
c
c          Inner
c
           id = idds(ir,1)
           gammas(ir,7,1,1)= -e2dflux(nks(ir)+1,ir)
     >                       * dds2(id) * costet(id)
     >                       / kbfst(ir,1)
           gammas(ir,7,2,1)= -e2dflux(nks(ir)+1,ir) * rp(idds(ir,1))
     >                 * dds2(id) * costet(id) / kbfst(ir,1)
c
           gammas(ir,7,1,3)= gammas(ir,7,1,1)+gammas(ir,7,1,2)
           gammas(ir,7,2,3)= gammas(ir,7,2,1)+gammas(ir,7,2,2)

c
c        EGDE2D G1/2 - A - KAREA
c
c          Outer
c
c          Outer
c
           id = idds(ir,2)
           gammas(ir,8,1,2)= 0.5*(e2dflux(1,ir)+e2dflux(2,ir))
     >                       * dds2(id) * costet(id)
     >                       / kbfs(1,ir)
           gammas(ir,8,2,2)= 0.5*(e2dflux(1,ir)+e2dflux(2,ir))
     >                       * rs(1,ir)
     >                 * dds2(id) * costet(id) / kbfs(1,ir)
c
c          Inner
c
           id = idds(ir,1)
           gammas(ir,8,1,1)= -0.5*(e2dflux(nks(ir)+1,ir)
     >                            +e2dflux(nks(ir),ir))
     >                       * dds2(id) * costet(id)
     >                       / kbfs(nks(ir),ir)
           gammas(ir,8,2,1)= -0.5*(e2dflux(nks(ir)+1,ir)
     >                            +e2dflux(nks(ir),ir))
     >                           * rs(nks(ir),ir)
     >                 * dds2(id) * costet(id) / kbfs(nks(ir),ir)
c
           gammas(ir,8,1,3)= gammas(ir,8,1,1)+gammas(ir,8,1,2)
           gammas(ir,8,2,3)= gammas(ir,8,2,1)+gammas(ir,8,2,2)
c
c          That should be all of the gammas
c
c
c          Calculate the ionization source terms and all of
c          its components.
c
c
c        Now to sum up ionization.
c
c        For this specific case ONLY the integral will
c        be done simply on a cell by cell basis - and
c        so the R-value used can be the mean between the
c        two end-points.
c
c
         srcsume  = 0.0
         srcsumre = 0.0
         srcsum  = 0.0
         srcsumr = 0.0
         srchalfe = 0.0
         srchalfre = 0.0
         srchalf = 0.0
         srchalfr = 0.0
c
         srcsuma = 0.0
         srchalfa= 0.0
         srcsumra = 0.0
         srchalfra= 0.0
c
         areasuma = 0.0
         areahalf = 0.0
         areasumr = 0.0
         areahalfr= 0.0
c
         asuma = 0.0
         asumra = 0.0
         ahalfa = 0.0
         ahalfra = 0.0

c
c
         do ik = 1,nks(ir)
c
            rmean = (krb(ik-1,ir) + krb(ik,ir)) /2.0
c
            rmean1 = (krb(ik-1,ir) + rs(ik,ir)) /2.0
c
            rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
c
            ds1 = kss2(ik,ir) - ksb(ik-1,ir)
            ds2 = ksb(ik,ir)  - kss2(ik,ir)
            ds = ds1 + ds2

            ds1h = ds1
            ds2h = ds2

            da = karea2(ik,ir)
            da1 = da /2.0
            da2 = da /2.0
            dah1 = da / 2.0
            dah2 = da / 2.0

            if (ik.eq.1) then
               ds1h = 0.0
               dah1 = 0.0
            endif
            if (ik.eq.nks(ir)) then
               ds2h = 0.0
               dah2 = 0.0
            endif

            dsh = ds1h + ds2h
            dah = dah1 + dah2

            srcsumr = srcsumr + pinion(ik,ir)*ds1*rmean1
            srcsumre = srcsumre + e2dion(ik,ir)*ds1*rmean1
            srcsum = srcsum  + pinion(ik,ir)*ds1
            srcsume = srcsume + e2dion(ik,ir)*ds1
c
            srcsuma = srcsuma + e2dion(ik,ir) * da1
            srcsumra = srcsumra + e2dion(ik,ir) * da1 * rmean1
            srchalfa = srchalfa + e2dion(ik,ir) * dah1
            srchalfra = srchalfra + e2dion(ik,ir) * dah1 * rmean1
c
            srchalfr = srchalfr + pinion(ik,ir)*ds1h*rmean1
            srchalfre = srchalfre + e2dion(ik,ir)*ds1h*rmean1
            srchalf = srchalf + pinion(ik,ir)*ds1h
            srchalfe = srchalfe + e2dion(ik,ir)*ds1h
c
            areasuma = areasuma + ds1
            areahalf = areahalf + ds1h
            areasumr = areasumr + rmean1 * ds1
            areahalfr = areahalfr + rmean1 * ds1h
c
            asuma = asuma + da1
            asumra = asumra + da1 * rmean1
            ahalfa = ahalfa + dah1
            ahalfra = ahalfra + dah1 * rmean1
c
c
c           Do in two parts to get integral at cell centre.
c
            intsrc(ik,ir,1) = srcsum
            intsrc(ik,ir,2) = srcsumr
            intsrc(ik,ir,3) = srcsume
            intsrc(ik,ir,4) = srcsumre
            intsrc(ik,ir,5) = srchalf
            intsrc(ik,ir,6) = srchalfr
            intsrc(ik,ir,7) = srchalfe
            intsrc(ik,ir,8) = srchalfre
c
            intsrc(ik,ir,9) =  srcsuma
            intsrc(ik,ir,10) = srchalfa
            intsrc(ik,ir,11) = srcsumra
            intsrc(ik,ir,12) = srchalfra
c
            intarea(ik,ir,1) = areasuma
            intarea(ik,ir,2) = areahalf
            intarea(ik,ir,3) = areasumr
            intarea(ik,ir,4) = areahalfr
c
            intarea(ik,ir,5) = asuma
            intarea(ik,ir,6) = ahalfa
            intarea(ik,ir,7) = asumra
            intarea(ik,ir,8) = ahalfra
c
c
            srcsumr = srcsumr + pinion(ik,ir)*ds2*rmean2
            srcsumre = srcsumre + e2dion(ik,ir)*ds2*rmean2
            srcsum = srcsum  + pinion(ik,ir)*ds2
            srcsume = srcsume + e2dion(ik,ir)*ds2
c
            srcsuma = srcsuma + e2dion(ik,ir) * da2
            srcsumra = srcsumra + e2dion(ik,ir) * da2 * rmean2
            srchalfa = srchalfa + e2dion(ik,ir) * dah2
            srchalfra = srchalfra + e2dion(ik,ir) * dah2 * rmean2
c
            srchalfr = srchalfr + pinion(ik,ir)*ds2h*rmean2
            srchalfre = srchalfre + e2dion(ik,ir)*ds2h*rmean2
            srchalf = srchalf + pinion(ik,ir)*ds2h
            srchalfe = srchalfe + e2dion(ik,ir)*ds2h
c
            areasuma = areasuma + ds2
            areahalf = areahalf + ds2h
            areasumr = areasumr + rmean2 * ds2
            areahalfr = areahalfr + rmean2 * ds2h
c
            asuma = asuma + da2
            asumra = asumra + da2 * rmean2
            ahalfa = ahalfa + dah2
            ahalfra = ahalfra + dah2 * rmean2
c
c           Calculate Areas and other components.
c
            srccomp(ik,ir,1) = pinion(ik,ir)
            srccomp(ik,ir,2) = e2dion(ik,ir)
            srccomp(ik,ir,3) = pinion(ik,ir) * ds
            srccomp(ik,ir,4) = e2dion(ik,ir) * ds
            srccomp(ik,ir,5) = pinion(ik,ir) * dsh
            srccomp(ik,ir,6) = e2dion(ik,ir) * dsh
            srccomp(ik,ir,7) = pinion(ik,ir) * rmean
            srccomp(ik,ir,8) = e2dion(ik,ir) * rmean
            srccomp(ik,ir,9) = pinion(ik,ir) * rmean * ds
            srccomp(ik,ir,10) = e2dion(ik,ir) * rmean * ds
            srccomp(ik,ir,11) = pinion(ik,ir) * rmean * dsh
            srccomp(ik,ir,12) = e2dion(ik,ir) * rmean * dsh
c
            srccomp(ik,ir,17) = e2dion(ik,ir) * da
            srccomp(ik,ir,18) = e2dion(ik,ir) * dah
            srccomp(ik,ir,19) = e2dion(ik,ir) * rmean * da
            srccomp(ik,ir,20) = e2dion(ik,ir) * rmean * dah
c
            srccomp(ik,ir,13) = rmean * ds
            srccomp(ik,ir,14) = rmean * dsh
            srccomp(ik,ir,15) = rmean * da
            srccomp(ik,ir,16) = rmean * dah
c
            srccomp(ik,ir,21) = ds
            srccomp(ik,ir,22) = dsh
c
c
c
            intsrc(ik,ir,13) = srcsuma
            intsrc(ik,ir,14) = srchalfa
            intsrc(ik,ir,15) = srcsumra
            intsrc(ik,ir,16) = srchalfra
c
            intarea(ik,ir,9) = asuma
            intarea(ik,ir,10) = ahalfa
            intarea(ik,ir,11) = asumra
            intarea(ik,ir,12) = ahalfra
c

         end do
c
         intsrc(nks(ir)+1,ir,1) = srcsum
         intsrc(nks(ir)+1,ir,2) = srcsumr
         intsrc(nks(ir)+1,ir,3) = srcsume
         intsrc(nks(ir)+1,ir,4) = srcsumre
         intsrc(nks(ir)+1,ir,5) = srchalf
         intsrc(nks(ir)+1,ir,6) = srchalfr
         intsrc(nks(ir)+1,ir,7) = srchalfe
         intsrc(nks(ir)+1,ir,8) = srchalfre
c
         intsrc(nks(ir)+1,ir,9) =  srcsuma
         intsrc(nks(ir)+1,ir,10) = srchalfa
         intsrc(nks(ir)+1,ir,11) = srcsumra
         intsrc(nks(ir)+1,ir,12) = srchalfra
c
         intsrc(nks(ir)+1,ir,13) = srcsuma
         intsrc(nks(ir)+1,ir,14) = srchalfa
         intsrc(nks(ir)+1,ir,15) = srcsumra
         intsrc(nks(ir)+1,ir,16) = srchalfra
c
         intarea(nks(ir)+1,ir,1) = areasuma
         intarea(nks(ir)+1,ir,2) = areahalf
         intarea(nks(ir)+1,ir,3) = areasumr
         intarea(nks(ir)+1,ir,4) = areahalfr
c
         intarea(nks(ir)+1,ir,5) = asuma
         intarea(nks(ir)+1,ir,6) = ahalfa
         intarea(nks(ir)+1,ir,7) = asumra
         intarea(nks(ir)+1,ir,8) = ahalfra
c
         intarea(nks(ir)+1,ir,9) = asuma
         intarea(nks(ir)+1,ir,10) = ahalfa
         intarea(nks(ir)+1,ir,11) = asumra
         intarea(nks(ir)+1,ir,12) = ahalfra
c
         newrconst(ir,1) = - (gammas(ir,1,1,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,1))
     >                     / intarea(nks(ir)+1,ir,1)
         newrconst(ir,2) = - (gammas(ir,1,2,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,2))
     >                     / intarea(nks(ir)+1,ir,3)
c
         newrconst(ir,3) = - (gammas(ir,2,1,3)
     >                     + intsrc(nks(ir)+1,ir,3))
     >                     / intarea(nks(ir)+1,ir,1)
         newrconst(ir,4) = - (gammas(ir,2,2,3)
     >                     + intsrc(nks(ir)+1,ir,4))
     >                     / intarea(nks(ir)+1,ir,3)
c
         newrconst(ir,5) = - (gammas(ir,4,1,3)
     >                     + intsrc(nks(ir)+1,ir,7))
     >                     / intarea(nks(ir)+1,ir,2)
         newrconst(ir,6) = - (gammas(ir,4,2,3)
     >                     + intsrc(nks(ir)+1,ir,8))
     >                     / intarea(nks(ir)+1,ir,4)
c
         newrconst(ir,7) = - (gammas(ir,4,1,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,5))
     >                     / intarea(nks(ir)+1,ir,2)
         newrconst(ir,8) = - (gammas(ir,4,2,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,6))
     >                     / intarea(nks(ir)+1,ir,4)
c
         newrconst(ir,9) = - (gammas(ir,5,1,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,7))
     >                     / intarea(nks(ir)+1,ir,2)
         newrconst(ir,10) = - (gammas(ir,5,2,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,8))
     >                     / intarea(nks(ir)+1,ir,4)
c
         newrconst(ir,11) = - (gammas(ir,6,1,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,7))
     >                     / intarea(nks(ir)+1,ir,2)
         newrconst(ir,12) = - (gammas(ir,6,2,3)
     >                     + rfact*intsrc(nks(ir)+1,ir,8))
     >                     / intarea(nks(ir)+1,ir,4)
c
         newrconst(ir,13) = - (gammas(ir,7,1,3)
     >                     + intsrc(nks(ir)+1,ir,9))
     >                     / intarea(nks(ir)+1,ir,5)
         newrconst(ir,14) = - (gammas(ir,7,2,3)
     >                     + intsrc(nks(ir)+1,ir,11))
     >                     / intarea(nks(ir)+1,ir,7)
c
         newrconst(ir,15) = - (gammas(ir,8,1,3)
     >                     + intsrc(nks(ir)+1,ir,10))
     >                     / intarea(nks(ir)+1,ir,6)
         newrconst(ir,16) = - (gammas(ir,8,2,3)
     >                     + intsrc(nks(ir)+1,ir,12))
     >                     / intarea(nks(ir)+1,ir,8)
c
c        Calculate the Fluxes at each grid point
c

         do ik = 1,nks(ir)

            fluxes(ik,ir,1) = (gammas(ir,7,1,2) + intsrc(ik,ir,9)
     >                +  newrconst(ir,13) * intarea(ik,ir,5))
     >                *  srccomp(ik,ir,21)/karea2(ik,ir)
            fluxes(ik,ir,2) = ((gammas(ir,7,2,2) + intsrc(ik,ir,11)
     >                +  newrconst(ir,14) * intarea(ik,ir,7))
     >                *  srccomp(ik,ir,21)/karea2(ik,ir))
     >                /  rs(ik,ir)
            fluxes(ik,ir,3) = (gammas(ir,8,1,2) + intsrc(ik,ir,10)
     >                +  newrconst(ir,15) * intarea(ik,ir,6))
     >                *  srccomp(ik,ir,22)/karea2(ik,ir)/rs(ik,ir)
            fluxes(ik,ir,4) = ((gammas(ir,8,2,2) + intsrc(ik,ir,12)
     >                +  newrconst(ir,16) * intarea(ik,ir,8))
     >                *  srccomp(ik,ir,22)/karea2(ik,ir)/rs(ik,ir))
     >                /  rs(ik,ir)
c
c           Fluxes for regular Edge2d full and 1/2 cell cases
c           not using kareas - normal and R-corrected. Cell centred.
c
            fluxes(ik,ir,9) = (gammas(ir,2,1,2) + intsrc(ik,ir,3)
     >                +  newrconst(ir,3) * intarea(ik,ir,1))
            fluxes(ik,ir,10) = (gammas(ir,2,2,2) + intsrc(ik,ir,4)
     >                +  newrconst(ir,4) * intarea(ik,ir,3))
     >                /  rs(ik,ir)
            fluxes(ik,ir,11) = (gammas(ir,4,1,2) + intsrc(ik,ir,7)
     >                +  newrconst(ir,5) * intarea(ik,ir,2))
            fluxes(ik,ir,12) = (gammas(ir,4,2,2) + intsrc(ik,ir,8)
     >                +  newrconst(ir,6) * intarea(ik,ir,4))
     >                /  rs(ik,ir)

c
c           Fluxes for regular DIVIMP/PIN full and 1/2 cell cases
c           not using kareas - normal and R-corrected. Cell centred.
c
            fluxes(ik,ir,13) = (gammas(ir,1,1,2) + intsrc(ik,ir,1)
     >                +  newrconst(ir,1) * intarea(ik,ir,1))
            fluxes(ik,ir,14) = (gammas(ir,1,2,2) + intsrc(ik,ir,2)
     >                +  newrconst(ir,2) * intarea(ik,ir,3))
     >                /  rs(ik,ir)
            fluxes(ik,ir,15) = (gammas(ir,4,1,2) + intsrc(ik,ir,5)
     >                +  newrconst(ir,7) * intarea(ik,ir,2))
            fluxes(ik,ir,16) = (gammas(ir,4,2,2) + intsrc(ik,ir,6)
     >                +  newrconst(ir,8) * intarea(ik,ir,4))
     >                /  rs(ik,ir)

         end do
c
         do ik = 1,nks(ir)+1

            if (ik.eq.1) then
              fluxes(ik,ir,5) = gammas(ir,7,1,2)
     >                      *srccomp(1,ir,21)/karea2(1,ir)
              fluxes(ik,ir,6) = gammas(ir,7,2,2)
     >                      *srccomp(1,ir,21)/karea2(1,ir)
     >                      / rp(idds(ir,2))
              fluxes(ik,ir,7) = gammas(ir,8,1,2)
              fluxes(ik,ir,8) = gammas(ir,8,2,2) / rs(ik,ir)
            else
              fluxes(ik,ir,5) = (gammas(ir,7,1,2) + intsrc(ik-1,ir,13)
     >                +  newrconst(ir,13) * intarea(ik-1,ir,9))
     >                *  srccomp(ik-1,ir,21)/karea2(ik-1,ir)
              fluxes(ik,ir,6) = (gammas(ir,7,2,2) + intsrc(ik-1,ir,15)
     >                +  newrconst(ir,14) * intarea(ik-1,ir,11))
     >                *  srccomp(ik-1,ir,21)/karea2(ik-1,ir)/rs(ik-1,ir)
              fluxes(ik,ir,7) = (gammas(ir,8,1,2) + intsrc(ik-1,ir,14)
     >                +  newrconst(ir,15) * intarea(ik-1,ir,10))
     >                *  srccomp(ik-1,ir,22)/karea2(ik-1,ir)/rs(ik-1,ir)
              fluxes(ik,ir,8) = (gammas(ir,8,2,2) + intsrc(ik-1,ir,16)
     >                +  newrconst(ir,16) * intarea(ik-1,ir,12))
     >                *  srccomp(ik-1,ir,22)/karea2(ik-1,ir)/rs(ik-1,ir)
            endif
c
         end do
c
c     And NOW to somehow print ALL of this in a usable fashion!!!
c
c
c     Print only for the separatrix ring for now.
c
      if (ir.eq.irsep.or.ir.eq.irsep+1.or.ir.eq.irsep+2.or.
     >    ir.eq.irsep+5.or.ir.eq.irsep+7.or.ir.eq.irsep+11) then
c
         write (6,*) 'ANALYSIS for ring:', ir
         write (6,*) ' The Recycling Fraction is applied to'//
     >       ' the ionization source calculated below - THEN '
         write (6,*) ' Used in the formulae for RCONST'
         write (6,*) ' RECYCLING FRACTION = ',rfact
c
c        Target Fluxes
c
         write (6,*) ' ALL Fluxes and Integrals are per meter'//
     >               ' toroidally and not for the whole torus.'
         write (6,*) ' The units of the target flux is: ions/m2/s'
c
c
         write (6,*) 'Target Fluxes - calculated various ways:'
         write (6,*) 'Flux name:           Inner           Outer'//
     >      '         Total         R-Inner         R-outer'//
     >      '         R-Total'

c
         write (6,300) 'DIVIMP 0',gammas(ir,1,1,1),gammas(ir,1,1,2),
     >               gammas(ir,1,1,3),gammas(ir,1,2,1),
     >               gammas(ir,1,2,2),gammas(ir,1,2,3)
c
         write (6,300) 'EDGE2D 0',gammas(ir,2,1,1),gammas(ir,2,1,2),
     >               gammas(ir,2,1,3),gammas(ir,2,2,1),
     >               gammas(ir,2,2,2),gammas(ir,2,2,3)
c
         write (6,300) 'EDGE2D G1/2',gammas(ir,4,1,1),gammas(ir,4,1,2),
     >               gammas(ir,4,1,3),gammas(ir,4,2,1),
     >               gammas(ir,4,2,2),gammas(ir,4,2,3)
c
         write (6,300) 'EDGE2D VP',gammas(ir,5,1,1),gammas(ir,5,1,2),
     >               gammas(ir,5,1,3),gammas(ir,5,2,1),
     >               gammas(ir,5,2,2),gammas(ir,5,2,3)
c
         write (6,300) 'EDGE2D VCAV',gammas(ir,6,1,1),gammas(ir,6,1,2),
     >               gammas(ir,6,1,3),gammas(ir,6,2,1),
     >               gammas(ir,6,2,2),gammas(ir,6,2,3)
c
         write (6,300) 'EDGE2D POL 0',gammas(ir,7,1,1),gammas(ir,7,1,2),
     >               gammas(ir,7,1,3),gammas(ir,7,2,1),
     >               gammas(ir,7,2,2),gammas(ir,7,2,3)
c
         write (6,300) 'EDGE2D POL.5',gammas(ir,8,1,1),gammas(ir,8,1,2),
     >               gammas(ir,8,1,3),gammas(ir,8,2,1),
     >               gammas(ir,8,2,2),gammas(ir,8,2,3)
c
c        Total Ionization values and integrated areas
c
         write (6,*) 'Ionization and AREA Totals for each type:'
         write (6,*) 'Name:            Total Src       Total R-Src'//
     >               '      Total Area      Total R-Area'
c
         write (6,300) 'DIVIMP FULL',intsrc(nks(ir)+1,ir,1),
     >              intsrc(nks(ir)+1,ir,2),intarea(nks(ir)+1,ir,1),
     >              intarea(nks(ir)+1,ir,3)
c
         write (6,300) 'EDGE2D FULL',intsrc(nks(ir)+1,ir,3),
     >              intsrc(nks(ir)+1,ir,4),intarea(nks(ir)+1,ir,1),
     >              intarea(nks(ir)+1,ir,3)
c
         write (6,300) 'DIVIMP 1/2',intsrc(nks(ir)+1,ir,5),
     >              intsrc(nks(ir)+1,ir,6),intarea(nks(ir)+1,ir,2),
     >              intarea(nks(ir)+1,ir,4)
c
         write (6,300) 'EDGE2D 1/2',intsrc(nks(ir)+1,ir,7),
     >              intsrc(nks(ir)+1,ir,8),intarea(nks(ir)+1,ir,2),
     >              intarea(nks(ir)+1,ir,4)
c
         write (6,300) 'EDGE2D POL F',intsrc(nks(ir)+1,ir,9),
     >              intsrc(nks(ir)+1,ir,11),intarea(nks(ir)+1,ir,5),
     >              intarea(nks(ir)+1,ir,7)
c
         write (6,300) 'EDGE2D POL.5',intsrc(nks(ir)+1,ir,10),
     >              intsrc(nks(ir)+1,ir,12),intarea(nks(ir)+1,ir,6),
     >              intarea(nks(ir)+1,ir,8)

c
c
c        NET RCONST values
c
         write (6,*)  ' RCONST values:'
         write (6,*)  ' DIVIMP ionization sources * recycling factor'
         write (6,*)  'Name:           Rconst          R-Rconst'
c
         write (6,300) 'DIVIMP FULL',newrconst(ir,1),newrconst(ir,2)
         write (6,300) 'EDGE2D FULL',newrconst(ir,3),newrconst(ir,4)
         write (6,300) 'EDGE2D G1/2',newrconst(ir,5),newrconst(ir,6)
         write (6,300) 'DIVIMP G1/2',newrconst(ir,7),newrconst(ir,8)
         write (6,300) 'EDGE2D VP',newrconst(ir,9),newrconst(ir,10)
         write (6,300) 'EDGE2D VCAV',newrconst(ir,11),newrconst(ir,12)
         write (6,300) 'EDGE2D POL F',newrconst(ir,13),newrconst(ir,14)
         write (6,300) 'EDGE2D POL.5',newrconst(ir,15),newrconst(ir,16)
c
c
c        Now for the Ionization Source details!
c
         write (6,*)
         write (6,*) ' And NOW for the Ionization Source Details:'
c
         write (6,*) ' EDGE2D DATA - FULL RING -'//
     >               ' REGULAR and R-corrected - POL - KAREA'
         write (6,*) ' Integrals are to cell centres:'
         write (6,700)
         do ik = 1,nks(ir)
            write(6,600) ik,karea2(ik,ir),srccomp(ik,ir,21),
     >           intsrc(ik,ir,9),intsrc(ik,ir,11),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,15),
     >           srccomp(ik,ir,19),srccomp(ik,ir,17),
     >           intarea(ik,ir,5),intarea(ik,ir,7),
     >           fluxes(ik,ir,1),fluxes(ik,ir,2)
         end do
c
         write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//
     >               ' -  REGULAR and R-corrected - POL - KAREA'
         write (6,*) ' Integrals are to cell centres:'
         write (6,700)
         do ik = 1,nks(ir)
            write(6,600) ik,karea2(ik,ir),srccomp(ik,ir,22),
     >           intsrc(ik,ir,10),intsrc(ik,ir,12),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,16),
     >           srccomp(ik,ir,20),srccomp(ik,ir,18),
     >           intarea(ik,ir,6),intarea(ik,ir,8),
     >           fluxes(ik,ir,3),fluxes(ik,ir,4)
         end do
c
         write (6,*) ' EDGE2D DATA - FULL RING -'//
     >               ' REGULAR and R-corrected - POL - KAREA -'//
     >               ' AT CELL BOUNDARIES'
         write (6,*) ' Integrals are to cell boundaries:'
         write (6,750)
         do ik = 1,nks(ir)+1
            write(6,650) ik,karea2(ik,ir),srccomp(ik,ir,21),
     >           intsrc(ik,ir,13),intsrc(ik,ir,15),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,15),
     >           srccomp(ik,ir,19),srccomp(ik,ir,17),
     >           intarea(ik,ir,9),intarea(ik,ir,11),
     >           fluxes(ik,ir,5),fluxes(ik,ir,6),e2dflux(ik,ir),
     >           ksb(ik,ir)
         end do
         write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//
     >               ' -  REGULAR and R-corrected - POL - KAREA'//
     >               ' AT CELL BOUNDARIES'
         write (6,*) ' Integrals are to cell boundaries:'
         write (6,750)
         do ik = 1,nks(ir)+1
            write(6,650) ik,karea2(ik,ir),srccomp(ik,ir,22),
     >           intsrc(ik,ir,14),intsrc(ik,ir,16),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,16),
     >           srccomp(ik,ir,20),srccomp(ik,ir,18),
     >           intarea(ik,ir,10),intarea(ik,ir,12),
     >           fluxes(ik,ir,7),fluxes(ik,ir,8),e2dflux(ik,ir),
     >           ksb(ik,ir)
         end do
c
c
c
         write (6,*) ' EDGE2D DATA - FULL RING -'//
     >               ' REGULAR and R-corrected'
         write (6,*) ' Integrals are to cell centres:'
         write (6,450)
         do ik = 1,nks(ir)
            write(6,550) ik,kss2(ik,ir),
     >           intsrc(ik,ir,3),intsrc(ik,ir,4),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,13),
     >           srccomp(ik,ir,10),srccomp(ik,ir,4),fluxes(ik,ir,9),
     >           fluxes(ik,ir,10)
         end do
c
         write (6,*) ' EDGE2D DATA - RING - 1/2 CELLS'//
     >               ' -  REGULAR and R-corrected'
         write (6,*) ' Integrals are to cell centres:'
         write (6,450)
         do ik = 1,nks(ir)
            write(6,550) ik,kss2(ik,ir),
     >           intsrc(ik,ir,7),intsrc(ik,ir,8),
     >           e2dion(ik,ir),srccomp(ik,ir,8),srccomp(ik,ir,14),
     >           srccomp(ik,ir,12),srccomp(ik,ir,6),fluxes(ik,ir,11),
     >           fluxes(ik,ir,12)
         end do
c
c
c
         write (6,*) ' DIVIMP/PINION DATA - FULL RING'//
     >               ' - REGULAR and R-corrected'
         write (6,*) ' Integrals are to cell centres:'
         write (6,450)
         do ik = 1,nks(ir)
            write(6,550) ik,kss2(ik,ir),
     >           intsrc(ik,ir,1),intsrc(ik,ir,2),
     >           srccomp(ik,ir,1),srccomp(ik,ir,7),srccomp(ik,ir,13),
     >           srccomp(ik,ir,9),srccomp(ik,ir,3),
     >           fluxes(ik,ir,15),fluxes(ik,ir,16)
         end do
c
         write (6,*) ' DIVIMP/PINION DATA - RING - 1/2 CELLS'//
     >               ' -  REGULAR and R-corrected'
         write (6,*) ' Integrals are to cell centres:'
         write (6,450)
         do ik = 1,nks(ir)
            write(6,550) ik,kss2(ik,ir),
     >           intsrc(ik,ir,5),intsrc(ik,ir,6),
     >           srccomp(ik,ir,1),srccomp(ik,ir,7),srccomp(ik,ir,14),
     >           srccomp(ik,ir,11),srccomp(ik,ir,5),
     >           fluxes(ik,ir,15),fluxes(ik,ir,16)
         end do
c
         write (6,*) ' Cell AREA calculations for ring:',ir
c
         write(6,800)
         do ik = 1,nks(ir)
           in = korpg(ik,ir)
           write(6,900) ik,karea2(ik,ir),areap(in),
     >                hro(ik,ir)*drho(ik,ir)*
     >                hteta(ik,ir)*dthetag(ik,ir),
     >                hro(ik,ir),
     >                drho(ik,ir),hro(ik,ir)*drho(ik,ir),
     >                hteta(ik,ir),dthetag(ik,ir),
     >                hteta(ik,ir)*dthetag(ik,ir),
     >                karea2(ik,ir)/srccomp(ik,ir,21),
     >                hro(ik,ir)*drho(ik,ir)/kbfs(ik,ir)
         end do
c
 800     format(2x,'IK',7x,'KAREA2',8x,'AREAP', 5x,'HrdrHtdt',
     >         9x,'Hrho',9x,'dRHO',4x,'Hrho*dRHO',
     >         8x,'Hteta',6x,'dTHETAG',3x,'Hteta*dTHE',4x,'KAREA2/dS',
     >         2x,'HroDro*Bt/B')
 900     format(i4,11(1x,e12.5))
c
 300     format(a12,6(2x,e14.6))
c
 400     format(2x,'IK',3x,'S (m)',5x,'Intsrc',4x,'Intsrc*R',
     >          9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,
     >          'SRC * R * ds',4x,'SRC * ds')
 500     format(i4,f8.4,7e12.5)
c
 450     format(2x,'IK',3x,'S (m)',5x,'Intsrc',4x,'Intsrc*R',
     >          9x,'SRC',5x,'SRC * R',6x,'R * ds',1x,
     >          'SRC * R * ds',4x,'SRC * ds',8x,'G(s)',4x,'G(s,R)/R')
 550     format(i4,f8.4,9e12.5)
c
 600     format(i4,e12.5,f8.4,11e12.5)
 700     format(2x,'IK',6x,'Area ',4x,' S (m) ',3x,
     >          'Intsrc',4x,'Intsrc*R',
     >          9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,
     >          'SRC * R * ds',4x,'SRC * ds',4x,'Int Area',4x,
     >          'Int A (R)',6x,'G(s)',4x,'G(s,R)/R')
c
 650     format(i4,e12.5,f8.4,13e12.5)
 750     format(2x,'IK',6x,'Area ',4x,' S (m) ',3x,
     >          'Intsrc',4x,'Intsrc*R',
     >          9x,'SRC',5x,'SRC * R',6x,'R * ds',2x,
     >          'SRC * R * ds',4x,'SRC * ds',4x,'Int Area',4x,
     >          'Int A (R)',6x,'G(s) ',3x,'G(s,R)/R',5x,'E2DFlux',
     >          5x,'S bound')
c
      endif
c
c        Close IR loop
c
      end do



      return
      end
c
c
c
      real*8 function majrpos(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     MAJRPOS: This routine returns the major radius
c              position corresponding to the S-value
c              passed to it. The sbnd and rbnd arrays
c              that are used here are initialized before
c              the solver is called and are adjusted
c              from teh DIVIMP values to work away from
c              each of the targets.
c
c              Rbnd and Sbnd give the R and S boundary
c              positions for each cell ik. The array
c              starts with a zero index so that the
c              the boundaries of each cell can be represented
c              by the values at IK and IK-1. These boundary
c              values are calculated from the underlying
c              polygonal grid and thus REQUIRE that such a grid
c              is in use. Errors will result if this information
c              is not available. The values for cell boundaries
c              are calculated in the TAU.D4A module at the
c              same time as the KSS2 values.
c
c              David Elder, Dec 12, 1995
c
c
c     Local variables
c
c
c     Since this routine will be called alot and searching will be
c     inefficient - it will store the last cell called and check that
c     first - after that it will perform a binary search on the sbnd array
c     in order to find the correct cell. The lastik value is reset whenever
c     S = 0 is encountered. It is initialized the first time to 1 - by the
c     data statement.
c
      integer lastik
      data lastik /1/
      integer in,bot,top,mid
c
c
c     Check if S in same cell as last time
c
      if (s.ge.sbnd(lastik-1).and.s.le.sbnd(lastik)) then
c
         majrpos = (rbnd(lastik) - rbnd(lastik-1)) *
     >             ((s-sbnd(lastik-1))/(sbnd(lastik)-sbnd(lastik-1)))
     >             + rbnd(lastik-1)
c
c     Not in last cell - Perform Binary search
c
      else
c
c        BOTTOM is the zeroth element in the SBND array and it should
c        always be S=0
c
         in  = 1
         bot = 1
         top = nptscopy
c
 100     continue
c
            mid = (bot+top)/2
c
            if (s.le.sbnd(mid)) then
               top = mid
            else
               bot = mid + 1
            endif
c
            if (bot.eq.top) then
               in = top
               goto 200
             endif
c
             in = in +1
c
             if (in.gt.2000) then
                write (6,*) 'Error in search:',in,bot,top,mid,s
                stop 100
             endif
         goto 100
c
c        Found cell
c

 200     majrpos = (rbnd(in) - rbnd(in-1)) *
     >             ((s-sbnd(in-1))/(sbnd(in)-sbnd(in-1)))
     >             + rbnd(in-1)
c
         lastik = in
c
c        Debug code to check for errors in search
c
         if ( (abs(majrpos-rbnd(in))+abs(majrpos-rbnd(in-1)))
     >         .ne.abs(rbnd(in)-rbnd(in-1))) then
            write (6,*) 'ERROR: MAJOR RADIUS INCORRECT'
            write (6,'(4i4)') nptscopy,bot,top,mid
            write (6,'(2i4,3g12.5)') in,lastik,majrpos,
     >                               rbnd(in-1),rbnd(in)
            write (6,'(3g12.4)') s,sbnd(in-1),sbnd(in)
            stop 101
         endif
c
      endif
c
      return
      end
c
c
c
      subroutine binsearch(s,in)
      implicit none
      real*8 s
      integer in
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     Search for right cell - record last in cell
c     and check that first since most routines will be
c     in the same cell more often than not.
c
c
c     Local Variables
c
      integer bot, top, mid, lastin
      data lastin /1/
c
c     Need to reset lastin if it is too large from
c     prevous rings.
c
      if (lastin.gt.nptscopy) lastin = 1
c
      if (lastin.eq.1.and.s.lt.sptscopy(lastin)) then
         in = 1
c slmod begin - new
c...ARRAY BOUNDS:
      elseif ((lastin.gt.1).and.(s.ge.sptscopy(MAX(1,lastin-1)))
     >       .and.(s.lt.sptscopy(MAX(1,lastin)))) then
c
c      elseif ((lastin.gt.1).and.(s.ge.sptscopy(lastin-1))
c     >       .and.(s.lt.sptscopy(lastin))) then
c slmod end
         in = lastin
      elseif (s.gt.sptscopy(nptscopy)) then
         in = nptscopy + 1
      else
         in  = 1
         bot = 1
         top = nptscopy
c
 100     continue
c
            mid = (bot+top)/2
c
            if (s.lt.sptscopy(mid)) then
               top = mid
            else
               bot = mid + 1
            endif
c
            if (bot.eq.top) then
               in = top
               goto 200
            endif
c
            in = in +1
            if (in.gt.2000) then
c               write (6,*) 'Error in search:',in,bot,top,mid,s
               stop 111
            endif
            goto 100
c
 200     continue

         lastin = in
c
      endif
c
      if (debug_s22) then 
         write(6,'(a,2i4,3(1x,g12.5))') 'BIN:',in,nptscopy,s,
     >                   sptscopy(in-1),sptscopy(in)
      endif
c
      return
      end
c
c
c
      subroutine setsw(setind,pplasma,new_errlevel)
      implicit none
      integer setind,pplasma,new_errlevel
      include 'solparams'
      include 'solswitch'
c
c     SETSW: The purpose of this routine is to set the ACTIVE
c            MODE switches that will be used for the specific
c            run of SOL option 22. This allows Sol option 22
c            to be run with different settings for each 1/2
c            ring - thus allowing failed solutions to be replaced
c            by simpler ones that will be more likely to work
c            correctly and allowing the user to specify
c            that problem rings be treated differently.
c
c     Option - setind = 0 - set the switches to the input values
c              setind = 1 - perform iterative error correction turning
c                           off various options depending on the actswerror
c                           switch and number of iterations.
c              setind = 2 - Turn off all options and solve using conduction
c                           ONLY
c
c
c              errlevel=10- turn off equipartition if it was activated
c              errlevel=9 - use uniform particles instead of d2n/dr2
c              errlevel=8 - use only PINQI cooling contributions
c              errlevel=7 - use 1/2 ring uniform power instead of whole
c              errlevel=6 - use 1/2 ring uniform power + 1/2 ring particles
c              errlevel=5 - 1/2 ring uniform particles + power at top
c              errlevel=4 - 4 + turn off v^2 convection term
c              errlevel=3 - 3 + no power terms
c              errlevel=2 - 2 + no convective terms
c              errlevel=1 - Conduction ONLY
c
c     Default- setind = 0 - switches set to input values
c
      integer maxerrs
      parameter (maxerrs = 10)
      integer errlevel,errlevels(maxerrs)
c
      new_errlevel = -1
c
c     Set up error cross-references
c
c     IPP/08 Krieger - errlevel was used without being initialized
c      write(6,*) 'ERR:',pplasma,setind,errlevel
c
      errlevel = 0
      write(6,*) 'ERR:',pplasma,setind
c
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
c
      if (setind.eq.-1) then
c
c        Load input values
c
         call setallsw(pplasma,0)
c
      elseif (setind.eq.-2.or.setind.gt.0) then
c
c        Set error level
c
         if (setind.eq.-2) then
            errlevel = actswerror
         else
            errlevel = setind
         endif
c
         write(6,*) 'ERR:BEG:',pplasma,setind,errlevel,
     >                        actswerror
c
c        Initially set all quantities to INPUT values
c
         call setallsw(pplasma,0)
c
c        ERRLEVEL = 10 = Switch OFF equipartition if it was
c                       turned ON
c
         if (errlevel.eq.10) then
c
            if (actswpei.eq.1.0) then
               actswpei = 0.0
            else
               errlevel = errlevels(errlevel)
            endif
c
         endif

c
c        ERRLEVEL = 9 = Switch to uniform gperp from d2n/dr2
c
         write(6,*) 'ERR:I1 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp


         if (errlevel.eq.9) then
c
c           Turn OFF equipartition (from level 10)
c
            actswpei = 0.0
c
c           Replace gradient proportional options with whole ring uniform
c
            if (actswgperp.eq.7.0.or.actswgperp.eq.8.0) then
c
               actswgperp = 2.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I2 :',pplasma,setind,errlevel,
     >               actswerror,actswpei,actswgperp,actswpcx
c
c        ERRLEVEL = 8 = Switch OFF any ION heating power
c
         if (errlevel.eq.8) then
c
c           Turn OFF equipartition
c
            actswpei = 0.0
c
c           Replace whole ring uniform with 1/2 ring uniform
c
            if (actswgperp.eq.7.0.or.actswgperp.eq.8.0)
     >              actswgperp = 2.0
c
c           Turn OFF PINQI heating if active
c
            if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) then
               actswpcx= 5.0
            else
               errlevel = errlevels(errlevel)
            endif
c
         endif
c
         write(6,*) 'ERR:I3 :',pplasma,setind,errlevel,
     >               actswerror,actswpei,actswgperp,actswpcx
c
c        ERRLEVEL = 7 = Eliminate whole ring uniform power options -
c
         if (errlevel.eq.7) then
c
c           Turn OFF equipartition
c
            actswpei = 0.0
c
c           Replace whole ring gradient proportional particles with
c           whole ring uniform option.
c
            if (actswgperp.eq.7.0.or.actswgperp.eq.8) actswgperp = 2.0
c
c           Eliminate PINQI heating
c
            if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
c
c           Eliminate whole ring uniform power options - replace with
c           1/2 ring uniform equivalents.
c
            if (actswpow.eq.4.0) then
                actswpow = 1.0
            elseif (actswpow.eq.6.0) then
                actswpow = 5.0
            else
               errlevel = errlevels(errlevel)
            endif
c
         endif
c
         write(6,*) 'ERR:I4 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp
c
c        ERRLEVEL = 6 - 1/2 ring uniform target power + particles
c
         if (errlevel.eq.6) then
c
c           Turn OFF equipartition
c
            actswpei = 0.0
c
c           Eliminate PINQI heating
c
            if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
c
c           Replace power options with uniform 1/2 ring unless
c           already set
c
            if (actswpow.ne.1.0.or.actswgperp.ne.1.0) then
c
               if (actswpow.ne.0.0) actswpow = 1.0
               actswgperp = 1.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I5 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp

c
c        ERRLEVEL = 5 - 1/2 ring uniform particles - all power in at top
c
         if (errlevel.eq.5) then
c
c           Turn OFF equipartition
c
            actswpei = 0.0
c
c           Eliminate PINQI heating
c
            if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
c
c           Replace power options with uniform 1/2 ring unless
c           already set
c
            if (actswpow.ne.0.0.or.actswgperp.ne.1.0) then
c
               actswpow   = 0.0
               actswgperp = 1.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I6 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp
c
c
c        ERRLEVEL = 4
c
c        Turn off second convective term as well as level 5+
c
         if (errlevel.eq.4) then
c
c
c           Turn OFF equipartition
c
            actswpei = 0.0
c
c           Eliminate PINQI heating
c
            if (actswpcx.eq.2.0.or.actswpcx.eq.3.0) actswpcx= 5.0
c
c           Replace power option - in at top
c                   particle option - 1/2 ring uniform
c
            actswpow   = 0.0
            actswgperp = 1.0
c
            if (actswconv.ne.0.0) then
c
               actswconv = 0.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I7 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp
c
c        ERRLEVEL = 3 - All of above + power terms turned off
c
         if (errlevel.eq.3) then
c
c           Reset switches for other error levels
c
            actswpei   = 0.0
            actswgperp = 1.0
            actswpow   = 0.0
            actswconv  = 0.0
c
c           Switch off power terms if necessary
c
            if (actswphelp.ne.0.0.or.actswpcx.ne.0.0.or.
     >          actswprad.ne.0.0.or.actswppion.ne.0.0.or.
     >          actswppelec.ne.0.0) then
c
               actswphelp = 0.0
               actswpei   = 0.0
               actswpcx   = 0.0
               actswprad  = 0.0
               actswppelec= 0.0
               actswppion = 0.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I8 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp

c
c        ERRLEVEL = 2 - All of above + no convection terms
c
         if (errlevel.eq.2) then

            actswgperp = 1.0
            actswphelp = 0.0
            actswpei   = 0.0
            actswpcx   = 0.0
            actswppelec= 0.0
            actswppion = 0.0
            actswprad  = 0.0
            actswpow   = 0.0
            actswconv  = 0.0
c
            if (actswcond.ne.0.0) then
c
               actswcond = 0.0
c
            else
c
               errlevel = errlevels(errlevel)
c
            endif
c
         endif
c
         write(6,*) 'ERR:I9 :',pplasma,setind,errlevel,
     >                   actswerror,actswpei,actswgperp

c
c        ERRLEVEL = 1 - ALL OFF
c
         if (errlevel.eq.1) then
c
c           SET switches to conduction only
c
            call setallsw(pplasma,1)
c
         endif
c
c        SET the error switch to be the next in sequence of the error
c        conditions to be used - the last is zero - meaning no more error
c        provisions are available.
c
         new_errlevel = errlevel
         actswerror = errlevels(errlevel)
c
         write(6,*) 'ERR:END:',pplasma,setind,errlevel,
     >                           actswerror
c
      elseif (setind.eq.0) then
c
c        Load maximum error settings
c
         call setallsw(pplasma,1)
c
      endif
c
c     Exit
c
      return
      end
c
c
c
      subroutine setallsw(pplasma,ind)
      implicit none
      integer ind,pplasma
      include 'solparams'
      include 'solswitch'
c
c     SETALLSW: This routine sets all of the switches to either their
c               input values or to all OFF except conduction. The
c               other options in the SETSW routine will tweak these
c               to obtain other error solution conditions betweeen
c               these two extremes.
c
c
      if (ind.eq.0) then
c
c
c        Set the switches differently for main SOL and private plasma
c
c        Main SOL
c
         if (pplasma.eq.0) then
c
c        Set all values to INPUT switch settings
c
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
            actswppress = switch(swppress)
c
c           PCX - DIVIMP PINQI - sub-option switches
c
            actswqidatiz = switch(swqidatiz)
            actswqidmliz = switch(swqidmliz)
            actswqidcx = switch(swqidcx)
            actswqidrec= switch(swqidrec)
c
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
            actswerror = switch(swerror)
c
c
c        Private plasma
c
         elseif (pplasma.eq.1) then
c
            if (switch(swion).eq.1.0.or.switch(swion).eq.2.0
     >          .or.switch(swion).eq.8.0) then
               actswion  = switch(swion)
               actswioni = switch(swionp)
            else
               actswion  = switch(swionp)
               actswioni = switch(swioni)
            endif
c
            actswcond = switch(swcond)
            actswconv = switch(swconv)
            actswprad = switch(swprad)
            actswphelp= switch(swphelp)
            actswpei  = switch(swpei)
            actswpcx  = switch(swpcx)
c
c           Both of these options are always off for now in the PP - at
c           some time they may be active but only (possibly) to distribute
c           the power transferred from the main SOL - however - this is
c           usually modelled using the already existing private plasma
c           power distribution methods.
c
            actswppelec = 0.0
            actswppion  = 0.0
            actswppress = 0.0
c
c            actswppelec  = switch(swppelec)
c            actswppion   = switch(swppion)
c            actswppress = switch(swppress)
c
c           PCX - DIVIMP PINQI - sub-option switches
c
            actswqidatiz = switch(swqidatiz)
            actswqidmliz = switch(swqidmliz)
            actswqidcx = switch(swqidcx)
            actswqidrec= switch(swqidrec)
c
            actswvisc1= switch(swvisc1)
            actswmach = switch(swmach)
            actswnmom = switch(swnmom)
            actswe2d  = switch(swe2d)
c
c           Set power to private plasma value
c
            actswpow  = switch(swpowp)
            actswgperp= switch(swgperpp)
            actswmajr = switch(swmajr)
            actswcore = switch(swcore)
            actswrecom= switch(swrecom)
            actswsmooth= switch(swsmooth)
            actswerror = switch(swerror)
c
         endif
c
      elseif (ind.eq.1) then

c
c        Ionization Source is EXPONENTIAL
c
         actswion= 0.0
c
c        Initial Ionization is set to exponential
c
         actswioni = 0.0
c
c        First convection term is OFF
c
         actswcond = 0.0
c
c        Kinetic convection term is OFF
c
         actswconv = 0.0
c
c        Radiation is OFF
c
         actswprad = 0.0
c
c        Phelpi is OFF
c
         actswphelp= 0.0
c
c        Pei is OFF
c
         actswpei = 0.0
c
c        Pcx is OFF
c
         actswpcx = 0.0
c
c           PCX - DIVIMP PINQI - sub-option switches - all OFF
c
         actswqidatiz = 0.0
         actswqidmliz = 0.0
         actswqidcx = 0.0
         actswqidrec= 0.0
c
c        Private plasma loss compensation options are OFF
c
         actswppelec = 0.0
         actswppion = 0.0
         actswppress = 0.0
c
c        Viscosity is OFF
c
         actswvisc1 = 0.0
c
c        Mach solver is OFF
c
         actswmach = 0.0
c
c        Neutral Momentum Loss is OFF
c
         actswnmom = 0.0
c
c        Edge2D Compatibility is set to INPUT value
c
         actswe2d = switch(swe2d)
c
c        Distributed Power is OFF
c
         actswpow = 0.0
c
c        Cross-field correction is OFF
c
         actswgperp = 0.0
c
c        Major Radius Correction is OFF
c
         actswmajr = 0.0
c
c        Core Option is OFF
c
         actswcore = 0.0
c
c        Recombination option is OFF
c
         actswrecom = 0.0
c
c        Smoothing is set to INPUT Value
c
         actswsmooth = switch(swsmooth)
c
c        Error handling switch is now OFF
c
         actswerror = 0.0
c

      endif
c
      return
      end
c
c
c
      subroutine initlen
      implicit none
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     INITLEN: This subroutine initializes the source lengths for
c             the ionization and radiation sources.
c
c
      if (lensind.eq.0) then
c
c        Absolute source distances
c
         if (lensst.ge.lensfi) lensst = 0.0
c
         ssrcst = min(lensst,halfringlen)
         ssrcfi = min(lensfi,halfringlen)
c
         if (ssrcst.ge.ssrcfi) ssrcst = 0.0
c
         ssrclen = (ssrcfi -ssrcst)
         ssrcmid = (ssrcfi + ssrcst) / 2.0
c
         ssrcdecay = min(lams,halfringlen)
c
      elseif (lensind.eq.1) then
c
c        Relative Source distances
c
c
c        Set up the ionization source limits for options where it is
c        required.
c
         lensst = min(lensst,0.5d0)
         lensfi = min(lensfi,0.5d0)
c
         if (lensst.ge.lensfi) lensst = 0.0
c
         ssrcst = lensst * ringlen
         ssrcfi = lensfi * ringlen
c
         if (ssrcst.ge.ssrcfi) ssrcst = 0.0
c
         ssrclen = (ssrcfi -ssrcst)
         ssrcmid = (ssrcfi + ssrcst) / 2.0
c
         ssrcdecay = min(lams*ringlen,halfringlen)
c
      endif
c
c     If using PIN data need to set the length of the ionization
c     source to be the entire 1/2 ring - since the ionization
c     will fall to zero at whatever point the PIN data
c     stipulates - not an arbitrary point.
c
      if (pinavail.and.
     >     (actswion.eq.1.0.or.actswion.eq.2.0)) then
         ssrcst = 0.0d0
         ssrcfi = halfringlen
         ssrclen = (ssrcfi -ssrcst)
         ssrcmid = (ssrcfi + ssrcst) / 2.0
      endif
c
      write (6,*) 'Lensrc:',lensst,lensfi,lenri,lenr
      write (6,*) 'Ssrc  :',ssrcst,ssrcfi,ssrcmid,ssrclen,ringlen
c
c     Adjust Radiation Source length
c
      lenr = min(lenri,halfringlen)
c
      return
      end
c
c
c
      real*8 function gperpf(s)
      implicit none
      real*8 s
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
c     GPERPF: This function returns the integrated perpendicular
c             component of the flux from zero to S.
c
c
c
      integer in
      real*8 gperp_add
c
      if (actswgperp.eq.0.0) then
         gperpf = 0.0
      elseif (actswgperp.eq.1.0.or.actswgperp.eq.2.0) then
         gperpf = gperpcor * (s-soffset)
      elseif (actswgperp.eq.3.0.or.actswgperp.eq.4.0.or.
     >        actswgperp.eq.7.or.actswgperp.eq.8.0) then
c
c        The overall correction factor is non-zero if there
c        is a perpendicular flux component.
c
         if (gperpcor.ne.0.0) then

c
c           Search for right cell
c
            call binsearch(s,in)
c
            if (in.eq.1) then
               gperpf =  (intgperp(in) * s / sptscopy(1))
            else
               gperpf = ( intgperp(in-1) +
     >             ( (intgperp(in)-intgperp(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
            endif
c
c        Set to zero in alternate case.
c
         else
c
            gperpf = 0.0
c
         endif
c
      elseif (actswgperp.eq.5.0.or.actswgperp.eq.6.0) then
c
c        Options 5 and 6 are the sum of two sources -
c        One square and one uniform
c        Option 5 - 1/2 ring balance
c        Option 6 - Whole ring balance
c
         if (s.lt.sgperpbeg) then
            gperpf = gperpcor * (s-soffset)
         elseif (s.ge.sgperpbeg.and.s.le.sgperpend) then
            gperpf = gperpcor * (s-soffset) +
     >               gperpcor2 * (s-sgperpbeg)
         elseif (s.gt.sgperpend) then
            gperpf = gperpcor * (s-soffset) +
     >               gperpcor2 * (sgperpend-sgperpbeg)
         endif
c
      endif
c
c
c
c     ADD in additional Gperp Source/sink terms
c
      if (switch(swextra).gt.0.0) then
c
         gperp_add = 0.0
c
c        Additional source
c
         if (s.gt.start_gextra_src.and.s.lt.stop_gextra_src) then
c
            gperp_add = gperp_add +
     >                  gextra_src * ( s-start_gextra_src)
c
         elseif (s.ge.stop_gextra_src) then
c
            if (start_gextra_src.ne.stop_gextra_src) then
c
               gperp_add = gperp_add + gextra_src *
     >                    (stop_gextra_src-start_gextra_src)
c
            else
c
               gperp_add = gperp_add + gextra_src
c
            endif
c
         endif
c
c        Additional Sink
c
         if (s.gt.start_gextra_sink.and.s.lt.stop_gextra_sink) then
c
            gperp_add = gperp_add +
     >                  gextra_sink * ( s-start_gextra_sink)
c
         elseif (s.ge.stop_gextra_sink) then
c
            if (start_gextra_src.ne.stop_gextra_src) then
c
               gperp_add = gperp_add + gextra_sink *
     >                    (stop_gextra_sink-start_gextra_sink)
c
            else
c
               gperp_add = gperp_add + gextra_sink
c
            endif
c
         endif
c
c        Modify net cross-field flux
c
         gperpf = gperpf + gperp_add
c
c
c     Endif for extra source
c
      endif
c
c      if (debug_s22) then
c         write(6,'(a,g12.5,g12.5)') 'GPERPF:',s,gperpf
c      endif
c
      return
      end
c
c
c
      real*8 function nhs_s(s)
      implicit none
c
c     WF'96: This returns the value of neutral density interpolated,
c     for s from
c     ( swnmom = 9 ) =>  nhs(ik), total neutral density
c     ( swnmom = 10 ) =>  nhs0(ik), ie. the stagnant (primary)
c                         neutral density.
c     NB: the first cell density is taken as constant
c
      real*8 s
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      include 'params'
      include 'cadas'
c
      integer i,in,bot,top,mid
c
      if (.not.pinavail) then
        write(6,*) 'nh = ?, PIN not avial'
        stop
      endif
c
      if (pinavail) then
c
c        Search for right cell
c
         call binsearch(s,in)
c
         if (switch(swnmom).eq.9.0) then
            if (in.eq.1) then
               nhs_s = nhs(1)
            else
               nhs_s = 1.0 *
     >             ( nhs(in-1) +
     >             ( (nhs(in) - nhs(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
            endif
c
c        write(6,*) 's/spt,nhs/nhs',s/sptscopy(in),nhs_s/nhs(in)
c
         elseif (switch(swnmom).eq.10.0) then
            if (in.eq.1) then
               nhs_s = nhs0(1)
            else
               nhs_s = 1.0 *
     >             ( nhs0(in-1) +
     >             ( (nhs0(in) - nhs0(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
            endif
c
c        write(6,*) 's/spt,nhs/nhs0',s/sptscopy(in),nhs_s/nhs0(in)
c
         endif
c
      else
c
c       Code has reached an error condition and should stop.
c
        write (6,*) 'ERROR in SOLASCV:'//
     >              ' Invalid neutral density (s<0)'
        stop
      endif
c
      return
      end
c
c
      real*8 function ths_s(s)
      implicit none
c
c     WF'96: This returns the value of neutral temperature interpolated,
c     for s from ths(ik).
c
c     NB: the first cell value is taken as constant
c
      real*8 s
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      include 'params'
      include 'cadas'
c
      integer i,in,bot,top,mid
c
      if (.not.pinavail) then
        write(6,*) 'th = ?, PIN not avail'
        stop 'function ths_s'
      endif
c
      if (pinavail) then
c
c        Search for right cell
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            ths_s = ths(in)
         else
            ths_s = 1.0 *
     >             ( ths(in-1) +
     >             ( (ths(in) - ths(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)) ))
         endif
c
c        write(6,*) 's/spt,ths/ths',s/sptscopy(in),ths_s/ths(in)
c
      else
c
c       Code has reached an error condition and should stop.
c
        write (6,*) 'ERROR in SOLASCV:'//
     >              ' Invalid neutral density (s<0)'
        stop
      endif
c
      return
      end
c
c
c
      real*8 function scx_s(s,n,Ti)
      implicit none
c
c     WF'96: returns the charge exchange mom. loss term (Pa/m)
c
      real*8 s,n,Ti
      external cxsig
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      real*8  sigmavcx,sigma,nhs_s,Ths_s,gamma,
     >    vp_tmp,Eh_tmp,nh_tmp,ga_tmp,wn_iso,mass_tmp
      real get_bg_mass
      external nhs_s,Ths_s,gamma,get_bg_mass
c
c     Initialize mass to mass of bg plasma ions
c
      mass_tmp = get_bg_mass()
c
      if (.not.pinavail) then
        scx_s = 0.0d0
        return
      endif
c
c      wn_iso = 1.0/sqrt(3.0d0)
c
      wn_iso = 0.57735d0
      Eh_tmp = 1.5 * Ths_s(s)
      nh_tmp = nhs_s(s)
      ga_tmp = - gamma(s)
c
      if (n.gt.0.0) then
         vp_tmp = abs(ga_tmp) / n
      else
         vp_tmp = abs(ga_tmp) / n
      endif
c
      if (ga_tmp.eq.0.0.or.nh_tmp.eq.0.0) then
         scx_s = 0.0d0
         return
      endif
c
      call cxsig(1.0d0, mass_tmp, Eh_tmp, wn_iso, wn_iso,
     >    wn_iso, Ti, vp_tmp, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 1,
     >    sigma, sigmavcx)
c
       scx_s =  mb * mconv * ga_tmp * nh_tmp * sigmavcx
c
c      if (s.ge.39.0) then
c      write(6,*) 'nh,th:', nh_tmp,Eh_tmp
c      write(6,*) 'ga,vp:', ga_tmp,vp_tmp
c      write(6,*) 'n,Ti:',n,Ti
c      write(6,*) 'scx_s:',scx_s
c      write(6,*) '---'
c      end if
c
c      if (s.ge.39.0) then
c        write(6,*) 's>39 in scx_s: stop'
c        stop
c      endif
c
      return
      end
c
c     
c
      real function get_bg_mass()
      implicit none
      include 'params'
      include 'comtor'  
c
c     GET_BG_MASS: This routine returns the mass of the
c                  background plasma ions in AMU.  
c
      get_bg_mass = crmb
c
      return
      end    
c
cc
      real*8 function scxupdt(s,n,nold,Ti,Tiold)
      implicit none
c
c     WF'96: This returns the charge exchange momentum loss -
c     approximately integrated to the current s.
c
c     This function updates lasts and lastscx values after
c     calculating a Scx contribution based on the average values
c     over the interval.
c
      real*8 s,n,nold,Ti,Tiold
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      common /scx/ lastscx,lasts
c
      real*8 lasts,lastscx,scx_s,avescx,newscx,oldscx
      external scx_s
c
      if (.not.pinavail) then
         scxupdt = 0.0d0
         return
      endif
c
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

      scxupdt = lastscx + rcxmom * avescx * (s-lasts)
c
c      if (s.ge.39.0) then
c      write(6,*) '-----'
c      write(6,*) 's:new,old',s,lasts
c      write(6,*) 'Ti:new,old',Ti,Tiold
c      write(6,*) 'n:new,old',n,nold
c      write(6,*) 'scx:new,old,ave',newscx,oldscx,avescx
c      write(6,*) 'scxupdt,lastscx:',scxupdt,lastscx
c      end if
c
      lastscx = scxupdt
      lasts = s
c
      return
      end
c
c
      real*8 function estscx(s,n,Ti)
      implicit none
c
c     WF'96: This returns the charge exchange momentum loss -
c     approximately integrated to the current point s.
c
c     This function returns an ESTIMATE of the Scx contribution at
c     the given S-position ... the lastpei and lasts values are
c     updated after every R-K iteration.
c
      real*8 s,n,Ti
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
      common /scx/ lastscx,lasts
c
      real*8 lasts,lastscx,scx_s,scx_tmp
      external scx_s
c
      if (.not.pinavail) then
         estscx = 0.0d0
         return
      endif
c
      if (s.eq.soffset.or.s.lt.lasts) then
         estscx = 0.0d0
         return
      elseif (s.eq.lasts) then
         estscx = lastscx
         return
      endif
c
      scx_tmp = scx_s(s,n,Ti)
      estscx = lastscx + rcxmom * scx_tmp * (s-lasts)
c
c      if (s.ge.40.0) then
c      write(6,*) '----'
c      write(6,*) 's-ls,scxtmp:',(s-lasts),scx_tmp
c      write(6,*) 'n,Ti',n,Ti
c      write(6,*) 'est,last',estscx,lastscx
c      end if
c
      return
      end
c
c     WF'96: THIS CLUSTER OF ROUTINES WAS ADAPTED FROM NIMBUS,
c     WITH THE AIM OF CALCULATING THE ION-NEUTRAL
c     CHARGE EXCHANGE COLLISION RATE. CALLED VIA
c
c     pmomloss() =>  estscx() => scx_s() => cxsig()
c
C     CHARGE EXCHANGE   H + H+ -> H+ + H
c
C I   IZ       = ATOMIC NUMBER
C I   AN       = MASS OF NEUTRAL (AMU)
C I   EN       = NEUTRAL ENERGY (EV)
C I   WNX/Y/Z  = NEUTRAL DIRECTION COSINES
C I   TI       = ION TEMPERATURE (EV)
C I   vp       = velocity of PLASMA FLOW
C I   AI       = ION MASS (AMU)
C I   WPX/Y/Z  = FLOW DIRECTION COSINES
C I   IXTYPE   = C.X. CROSS SECTION MODEL
C O   SIGCX    = C.X. CROSS SECTION (EFFECTIVE)
C O   SVCX     = C.X. REACTION RATE
C
C     RETURNS MICROSCOPIC CHARGE EXCHANGE X-SECTION  BETWEEN H-ISOTOPES
C     AVERAGED OVER ALL TARGET VELOCITIES (MAXWELL + FLOW) FOR NEUTRALS
C     IN A PLASMA REGION CONTAINING H+,D+,DT+,T+.
C
c     NB:  originally in cgs;  converted to mks / eV , both in/output.
c
c
      subroutine CXSIG (IZ,AN,EN,WNX,WNY,WNZ,
     >            TI,VP,AI,WPX,WPY,WPZ,IXTYPE,
     >            SIGCX,SVCX)
      implicit none
c
      integer ip,ixtype
      real*8 IZ,AN,EN,WNX,WNY,WNZ,TI,VP,AI,
     >    WPX,WPY,WPZ,SIGCX,SVCX
      real*8 vn,vnx,vny,vnz,vpx,vpy,vpz,vr2,
     >   ES,ENS,ESTS,C,RRATE,ER,RM
c
      COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      real*8 tis,c1,c2,c3,c4
      integer idqh
c
      real*8 xtfct
      external xtfct
c
      IF(IZ.ge.2.0) GOTO 300
C
C     VELOCITY OF THE NEUTRAL
c
      vn = SQRT(EN/AN) * 1.3841E+06 * 0.01
      vnx = vn * WNX
      vny = vn * WNY
      vnz = vn * WNZ
c
c      write(6,*) 'wnxy',wnx,wny,wnz
c      write(6,*) 'vn,xyz',vn,vnx,vny,vnz
C
C     VELOCITY OF THE PLASMA
c
C      vp = SQRT((TE+TI)/AI) * 0.978E+06 * MACH
C
      vp  = abs (vp)
      vpx = vp * WPX
      vpy = vp * WPY
      vpz = vp * WPZ
c
c      write(6,*) 'Wxyz',WPX,WPY,WPZ
c      write(6,*) 'vp,xyz',vp,vpx,vpy,vpz
C
C     MEAN RELATIVE VELOCITY BETWEEN PLASMA AND NEUTRAL
c
      IF(VP.LE.0.0) THEN
        vr2 = vn ** 2
      ELSE
        vr2 = (vpx-vnx)**2 + (vpy-vny)**2 + (vpz-vnz)**2
      ENDIF
c
C     SCALED SHIFTED NEUTRAL ENERGY (EV)
c
      ES = 0.52197E-12 * vr2 * 1.0E+4
c
C     SHIFTED NEUTRAL ENERGY
c
      ENS = ES * AN
c
C     SCALED ION TEMPERATURE
c
      TIS = TI / AI
C
C     RATIO OF SHIFTED NEUTRAL ENERGY TO SCALED ION TEMP.
C
      ESTS = ES / TIS
c
      IF (ES.GE.4.0E+4) GO TO 20
      IF (ES.GE.5.0E+3) GO TO 10
      IF (ESTS.GE.5.0) GO TO 80
      IP = 1
      GO TO 30
   10 IF(ESTS.GE.26.0) GO TO 80
      IP = 2
      GO TO 30
   20 IF(ESTS.GE.100.0) GO TO 80
      IP = 3
   30 C = SQRT(AI/AN)
      C1 = SQRT(ENS/TI)
      C2 = C*C1
      C3 = 0.0
      IF (C2**2.LT.174.0) C3 = EXP(-C2**2)
      C4 = 2.0*C2
      IDQH = 1
      GO TO ( 40 , 50 , 60 ),IP
   40 CALL DQH04P(XTFCT,RM,ixtype)
      GO TO 70
   50 CALL DQH12P(XTFCT,RM,ixtype)
      GO TO 70
   60 CALL DQH32P(XTFCT,RM,ixtype)
   70 RRATE = 0.780939E+6 * RM * TI/AI * SQRT(AN/ENS)
     >          * 1.0E-6
c
      if (vn.gt.0.0) then
         SIGCX = RRATE / vn
      else
         SIGCX = 0.0
c         write (6,*) 'ERROR: SOLASCV - Momentum Option 9 :'//
c     >               ' Vn =< 0.0 ',vn,sigcx
      endif
c
      GOTO 90
c
C     MEAN RELATIVE NEUTRAL/ION ENERGY
c
   80 ER = 1.5 * AN * TIS + ENS
C
      CALL HXHP(ER/AN,ixtype,SIGCX)
C
C     CONSERVE THE REACTION RATE
C
      SIGCX = SIGCX * SQRT(ER/EN) * 1.0E-4
   90 CONTINUE
      SVCX = SIGCX * vn
      RETURN
c
  300 WRITE(6,*) 'ERROR: NON-HYDROGENIC CALL FROM CXSIG'
      write(6,*) 'iz=', iz
      write(6,*) 'var...', ip,ixtype,IZ,AN,EN,WNX,WNY,WNZ
     >    ,TI,VP,AI,WPX,WPY,WPZ,SIGCX,SVCX,vn,vnx,vny,vnz
     >    ,vpx,vpy,vpz,vr2,ES,ENS,ESTS,C,RRATE,ER,RM
      STOP
      RETURN
c
      END
c
C
C     RETURNS MICROSCOPIC CROSS SECTION SIGMA FOR CHARGE-EXCHANGE
C     REACTION BETWEEN HYDROGEN ATOMS/IONS
C
C     ER=RELATIVE ENERGY
C
C                    MICROSCOPIC CROSS SECTION OF REACTION
C                    (H+)+H-->H+(H+) IS GIVEN
C
C     ICXTYP=1:JANEV - 'ELEMENTARY PROCESSES IN HYDROGEN-HELIUM PLASMA',
C                       SPRINGER (1987))
C
c
      SUBROUTINE HXHP(ER,IXTYP,SIGMA)
      implicit none
C
      real*8 ER,SIGMA,ESCALE
      integer ixtyp
      external PNFIT
C
      real*8 R318(9)
c        DIMENSION R318(9)
      DATA R318/-3.274123792568E+01,
     >          -8.916456579806E-02,
     >          -3.016990732025E-02,
     >           9.205482406462E-03,
     >           2.400266568315E-03,
     >          -1.927122311323E-03,
     >           3.654750340106E-04,
     >          -2.788866460622E-05,
     >           7.422296363524E-07/
C
      IF( IXTYP.EQ.1 ) THEN
C                                   JANEV (1987)
          ESCALE = ER
          IF(ESCALE.LT.0.1) ESCALE = 0.1
          CALL PNFIT(9,R318(1),ESCALE,SIGMA)
      ELSE
          ESCALE = ER
          IF(ESCALE.GT.200.0) GO TO 10
          IF(ESCALE.LT.0.001) ESCALE=0.001
C                                   GREENLAND,1984 (E.LE.200 EV)
c
c         Caution: Alog10 was used previouly
c
          SIGMA = 4.549E-15 - 1.033E-15 * LOG10(ESCALE)
          RETURN
C                                   RIVIERE  (E.GT.200 EV)
   10     SIGMA = 6.937E-15 * (1.0 - 0.155*LOG10(ESCALE))**2
          IF(ESCALE.GE. 15000.0) THEN
             SIGMA = SIGMA/(1.0 + 0.1112E-14 * (ESCALE**3.3))
          ENDIF
      END IF
C
      RETURN
      END
c
c
C
C     ..................................................................
C
C        SUBROUTINE DQH04P
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
C                               0 TO +INFINITY).
C
C        USAGE
C           CALL DQH04P (FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 4-POINT GAUSSIAN-HERMITE
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 15.
C           FOR REFERENCE, SEE
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
C           TR00.1100 (MARCH 1964), PP.213-214.
C
C     ..................................................................
C
      SUBROUTINE DQH04P(FCT,Y,ixtype)
      implicit none
C
      REAL*8  X,Y,FCT
      integer ixtype
c
      COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      real*8 tis,c1,c2,c3,c4
      integer idqh
C
      EXTERNAL FCT
C
      X=.29306374202572440D1
      Y=.19960407221136762D-3*FCT(X,ixtype)
      X=.19816567566958429D1
      Y=Y+.17077983007413475D-1*FCT(X,ixtype)
      X=.11571937124467802D1
      Y=Y+.20780232581489188D0*FCT(X,ixtype)
      X=.38118699020732212D0
      Y=Y+.66114701255824129D0*FCT(X,ixtype)
      RETURN
      END
c
c
C.......................................................................
C.......................................................................
C
C
C        SUBROUTINE DQH12P
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
C                               0 TO +INFINITY).
C
C        USAGE
C           CALL DQH12P (FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 12-POINT GAUSSIAN-HERMITE
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 47.
C           FOR REFERENCE, SEE
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
C           TR00.1100 (MARCH 1964), PP.213-214.
C
C     ..................................................................
C
      SUBROUTINE DQH12P(FCT,Y,ixtype)
      implicit none
C
      REAL*8 X,Y,FCT
      integer ixtype
c
      COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      real*8 tis,c1,c2,c3,c4
      integer idqh
C
      EXTERNAL FCT
C
      X=.60159255614257397D1
      Y=.16643684964891089D-15*FCT(X,IXTYPE)
      X=.52593829276680444D1
      Y=Y+.65846202430781701D-12*FCT(X,IXTYPE)
      X=.46256627564237873D1
      Y=Y+.30462542699875639D-9*FCT(X,IXTYPE)
      X=.40536644024481495D1
      Y=Y+.40189711749414297D-7*FCT(X,IXTYPE)
      X=.35200068130345247D1
      Y=Y+.21582457049023336D-5*FCT(X,IXTYPE)
      X=.30125461375655648D1
      Y=Y+.56886916364043798D-4*FCT(X,IXTYPE)
      X=.25238810170114270D1
      Y=Y+.8236924826884175D-3*FCT(X,IXTYPE)
      X=.20490035736616989D1
      Y=Y+.70483558100726710D-2*FCT(X,IXTYPE)
      X=.15842500109616941D1
      Y=Y+.37445470503230746D-1*FCT(X,IXTYPE)
      X=.11267608176112451D1
      Y=Y+.12773962178455916D0*FCT(X,IXTYPE)
      X=.67417110703721224D0
      Y=Y+.28617953534644302D0*FCT(X,IXTYPE)
      X=.22441454747251559D0
      Y=Y+.42693116386869925D0*FCT(X,IXTYPE)
      RETURN
      END
c
C.......................................................................
C.......................................................................
C
C
C        SUBROUTINE DQH32P
C
C        PURPOSE
C           TO COMPUTE INTEGRAL(EXP(-X*X)*FCT(X), SUMMED OVER X FROM
C                               0 TO +INFINITY).
C
C        USAGE
C           CALL DQH32P (FCT,Y)
C           PARAMETER FCT REQUIRES AN EXTERNAL STATEMENT
C
C        DESCRIPTION OF PARAMETERS
C           FCT    - THE NAME OF AN EXTERNAL DOUBLE PRECISION FUNCTION
C                    SUBPROGRAM USED.
C           Y      - THE RESULTING DOUBLE PRECISION INTEGRAL VALUE.
C
C        REMARKS
C           NONE
C
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED
C           THE EXTERNAL DOUBLE PRECISION FUNCTION SUBPROGRAM FCT(X)
C           MUST BE FURNISHED BY THE USER.
C
C        METHOD
C           EVALUATION IS DONE BY MEANS OF 32-POINT GAUSSIAN-HERMITE
C           QUADRATURE FORMULA, WHICH INTEGRATES EXACTLY WHENEVER
C           FCT(X) IS A POLYNOMIAL UP TO DEGREE 127.
C           FOR REFERENCE, SEE
C           SHAO/CHEN/FRANK, TABLES OF ZEROS AND GAUSSIAN WEIGHTS OF
C           CERTAIN ASSOCIATED LAGUERRE POLYNOMIALS AND THE RELATED
C           GENERALIZED HERMITE POLYNOMIALS, IBM TECHNICAL REPORT
C           TR00.1100 (MARCH 1964), PP.213-214.
C
C     ..................................................................
C
      SUBROUTINE DQH32P(FCT,Y,ixtype)
      implicit none
C
      REAL*8 X,Y,FCT
      integer ixtype
c
      COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      real*8 tis,c1,c2,c3,c4
      integer idqh
C
      EXTERNAL FCT
C
      X=.10526123167960546D2
      Y=.55357065358569428D-48*FCT(X,IXTYPE)
      X=.9895287586829539D1
      Y=Y+.16797479901081592D-42*FCT(X,IXTYPE)
      X=.9373159549646721D1
      Y=Y+.34211380112557405D-38*FCT(X,IXTYPE)
      X=.8907249099964770D1
      Y=Y+.15573906246297638D-34*FCT(X,IXTYPE)
      X=.8477529083379863D1
      Y=Y+.25496608991129993D-31*FCT(X,IXTYPE)
      X=.8073687285010225D1
      Y=Y+.19291035954649669D-28*FCT(X,IXTYPE)
      X=.7689540164040497D1
      Y=Y+.7861797788925910D-26*FCT(X,IXTYPE)
      X=.7321013032780949D1
      Y=Y+.19117068833006428D-23*FCT(X,IXTYPE)
      X=.69652411205511075D1
      Y=Y+.29828627842798512D-21*FCT(X,IXTYPE)
      X=.66201122626360274D1
      Y=Y+.31522545665037814D-19*FCT(X,IXTYPE)
      X=.62840112287748282D1
      Y=Y+.23518847106758191D-17*FCT(X,IXTYPE)
      X=.59556663267994860D1
      Y=Y+.12800933913224380D-15*FCT(X,IXTYPE)
      X=.56340521643499721D1
      Y=Y+.52186237265908475D-14*FCT(X,IXTYPE)
      X=.53183252246332709D1
      Y=Y+.16283407307097204D-12*FCT(X,IXTYPE)
      X=.50077796021987682D1
      Y=Y+.39591777669477239D-11*FCT(X,IXTYPE)
      X=.47018156474074998D1
      Y=Y+.7615217250145451D-10*FCT(X,IXTYPE)
      X=.43999171682281376D1
      Y=Y+.11736167423215493D-8*FCT(X,IXTYPE)
      X=.41016344745666567D1
      Y=Y+.14651253164761094D-7*FCT(X,IXTYPE)
      X=.38065715139453605D1
      Y=Y+.14955329367272471D-6*FCT(X,IXTYPE)
      X=.35143759357409062D1
      Y=Y+.12583402510311846D-5*FCT(X,IXTYPE)
      X=.32247312919920357D1
      Y=Y+.8788499230850359D-5*FCT(X,IXTYPE)
      X=.29373508230046218D1
      Y=Y+.51259291357862747D-4*FCT(X,IXTYPE)
      X=.26519724354306350D1
      Y=Y+.25098369851306249D-3*FCT(X,IXTYPE)
      X=.23683545886324014D1
      Y=Y+.10363290995075777D-2*FCT(X,IXTYPE)
      X=.20862728798817620D1
      Y=Y+.36225869785344588D-2*FCT(X,IXTYPE)
      X=.18055171714655449D1
      Y=Y+.10756040509879137D-1*FCT(X,IXTYPE)
      X=.15258891402098637D1
      Y=Y+.27203128953688918D-1*FCT(X,IXTYPE)
      X=.12472001569431179D1
      Y=Y+.58739981964099435D-1*FCT(X,IXTYPE)
      X=.9692694230711780D0
      Y=Y+.10849834930618684D0*FCT(X,IXTYPE)
      X=.69192230581004458D0
      Y=Y+.17168584234908370D0*FCT(X,IXTYPE)
      X=.41498882412107868D0
      Y=Y+.23299478606267805D0*FCT(X,IXTYPE)
      X=.13830224498700972D0
      Y=Y+.27137742494130398D0*FCT(X,IXTYPE)
      RETURN
      END
c
C
C     POLYNOMIAL FORM USED FOR MOLECULAR REACTIONS WITH ELECTRONS
C     AND JANEV C.X. X.S.
C
C     SVJ=RATE COEFFICIENT (CM**3/SEC)
C     T=TEMPERATURE (EV)
C     ALOG(SVJ)=A(1)+A(2)*ALOG(T)+...+A(NJ)*ALOG(T)**(NJ-1)
C              =A(1)+X*(A(2)+X*(A(3)+X*(...+X*(A(N-1)+A(NJ)*X)...)))
C
      SUBROUTINE PNFIT(NJ,A,T,SVJ)
      implicit none
c
c      DIMENSION A(NJ)
c      REAL*8 TT
C
      integer nj,i
      real*8 a(10),t,svj,tt,sv,x
c
      TT = T
      X = LOG(TT)
      SV = A(NJ)
      DO 10 I = NJ-1, 1, -1
        SV = SV*X + A(I)
   10 CONTINUE
      SVJ = EXP(SV)
      RETURN
      END
c
c
c
      REAL*8 FUNCTION XTFCT(X,IXTYPE)
      implicit none
C
      integer ixtype
      REAL*8 X,D,D1,D2,F,ER,SIGMA
C
c
      COMMON /DQHCOM/ TIS,C1,C2,C3,C4,idqh
      real*8 tis,c1,c2,c3,c4
      integer idqh
C
      ER = TIS * X * X
c
C      GO TO ( 10 , 20 , 30 , 40 ), IDQH
c
   10 CALL HXHP(ER,ixtype,SIGMA)
c
c      GO TO 50
C   20 CALL HXHEP(ER,SIGMA)
C      GO TO 50
C   30 CALL XIPP(ER,SIGMA)
C      GO TO 50
C   40 CALL ELXS(KES,ER*AR,SIGMA)
c
   50 IF(C3.EQ.0.0) GO TO 60
      D = C4*X
      IF(D.GT.174.0) GO TO 60
      D = EXP(D)
      F = (D - 1.0D0/D)*C3
      GO TO 70
   60 D1 = C2*(2.0D0*X-C2)
      D2 = -C2*(2.0D0*X+C2)
      F = EXP(D1) - EXP(D2)
   70 XTFCT = X*X*F*SIGMA
      RETURN
      END
c
c
c
      real*8 function pinqid(s,opt)
      implicit none
      real*8 s
      integer opt
c
      include 'solparams'
      include 'solswitch'
      include 'solcommon'
c
      include 'params'
      include 'cadas'
c
c     PINQID: This function/subroutine implements the
c             DIVIMP calculated version of PINQI.
c             The function serves three purposes - selected
c             by the "opt" switch.
c
c             opt = -1 ... Initialize the source calculations
c                          Pre-integrate any values required
c                          Setup quantities for later use.
c
c             opt = 0  ... Provide an estimate of the integrated
c                          value of PINQID to the current S
c                          position.
c
c             opt = 1  ... Update the value of the PINQID integral
c                          Depending on the sub-options
c                          chosen for the calculation of PINQID
c                          This function may not be required.
c
c
c     Local variables
c
      integer in,ik,ir
      real*8 qidatiz(mxspts),qidmliz(mxspts),qidcx(mxspts),
     >        qidrec(mxspts)
      real*8 eav,sigvcx,qidsum
c
      real*8 rcxmult
      external rcxmult
c
c     If option is turned OFF - return zero.
c
      if (actswpcx.ne.4.0) then
         pinqid = 0.0
         return
      endif
c
c     IF PIN is not available - then also return 0.0 - since
c     PINQID can not be calculated.
c
      if (.not.pinavail) then
         pinqid = 0.0
         return
      endif
c
c     INITIALIZATION
c
c     Setup and initialization - the value passed back is the
c     net integral of all of the PINQI components.
c
      if (opt.eq.-1) then
c
c        Initialize all arrays to zero - this is equivalent to
c        option zero in all cases - which turns the contribution
c        to QID due to that effect - OFF.
c
         call qzero(qidatiz,mxspts)
         call qzero(qidmliz,mxspts)
         call qzero(qidrec,mxspts)
         call qzero(qidcx,mxspts)
c
         DO IK = 1, nptscopy
            PTESA(IK) = oldte(ik)
            PNESA(IK) = oldne(ik) * zb
            PNBS(IK) = oldne(ik)
            PNHS(IK) = nhs(ik)
         ENDDO
c
c         write (6,*) 'Values OVERWRITTEN due to Tcutoff:'
c
c        Atomic ionization
c
         if (actswqidatiz.eq.1) then
c
            do ik = 1,nptscopy
               if (oldti(ik).lt.tcutatiz) then
                  qidatiz(ik) = 0.0
c                  write (6,100) 'AT:',ik,oldti(ik),tcutatiz
               elseif (nhs(ik).gt.0.0) then
                  qidatiz(ik) = 1.5 * ths(ik) * econv *
     >                      nhs(ik)/(nhs(ik)+nh2s(ik)) * ionsrc(ik)
               endif
            end do
c
         elseif (actswqidatiz.eq.2) then
c
            write(year,'(i2.2)') iyearh
            call xxuid(useridh)
c
            ICLASS = 2
            CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
c
            DO IK = 1, nptscopy
               if (oldti(ik).lt.tcutatiz) then
                  qidatiz(ik) = 0.0
c                  write (6,100) 'AT:',ik,oldti(ik),tcutatiz
               else
                  QIdATIZ(IK) = 1.5 * ths(ik) * econv
     >                 * oldne(ik) * zb * nhs(ik) * PCOEF(IK,1)
               endif
            ENDDO
c
         endif
c
c        Molecular Ionization
c
         if (actswqidmliz.eq.1) then
c
            do ik = 1,nptscopy
               if (oldti(ik).lt.tcutmliz) then
                  qidmliz(ik) = 0.0
c                  write (6,100) 'ML:',ik,oldti(ik),tcutmliz
               elseif (nh2s(ik).gt.0.0) then
                  qidmliz(ik) = 3.0 * econv *
     >                 nh2s(ik)/(nhs(ik)+nh2s(ik)) * ionsrc(ik)
               endif
            end do
c
         elseif (actswqidmliz.eq.2) then
c
            write(year,'(i2.2)') iyearh
            call xxuid(useridh)
c
            ICLASS = 2
            CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
c
            DO IK = 1, nptscopy
               if (oldti(ik).lt.tcutmliz) then
                  qidmliz(ik) = 0.0
c                  write (6,100) 'ML:',ik,oldti(ik),tcutmliz
               else
                  QIdmlIZ(IK) = 3.0 * econv
     >                * oldne(ik) * zb * nh2s(ik) * PCOEF(IK,1)
               endif
            ENDDO
c
         endif
c
c
c        Recombination
c
         if (actswqidrec.eq.1) then
c
            do ik = 1,nptscopy
               if (oldti(ik).lt.tcutrec) then
                  qidrec(ik) = 0.0
c                  write (6,100) 'RC:',ik,oldti(ik),tcutrec
               elseif (nh2s(ik).gt.0.0) then
                  qidrec(ik) = -1.5 * oldti(ik) * econv *
     >                            recsrc(ik)
               endif
            end do
c
         endif
c
c        Charge Exchange
c
         if (actswqidcx.eq.1) then
c
c           Prescription from TN1396
c
            do ik = 1,nptscopy
c
               if (oldti(ik).lt.tcutcx) then
                  qidcx(ik) = 0.0
c                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
               else
                  qidcx(ik) = 1.5 * trefcx*rcxmult(oldti(ik))*ceicf
     >                         * ionsrc(ik)
               endif
c
            end do

         elseif (actswqidcx.eq.2) then
c
            do ik = 1,nptscopy
c
               eav = 1.5 * (ths(ik)+oldti(ik)) / 2.0
c
               if (eav.gt.1000.0) then
                  sigvcx = 1.0d-13
               else
                  sigvcx = 1.0d-14 * eav**(1.0d0/3.0d0)
               endif
c
               if (oldti(ik).lt.tcutcx) then
                  qidcx(ik) = 0.0
c                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
               else
                  QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))
     >                   * econv * oldne(ik) * zb * nhs(ik) * sigvcx
               endif
c
               if (m0.eq.initm0) then
c
                  if (ik.eq.1) write(6,*) 'QID Source terms:'
c
                  write (6,'(i4,8(1x,g11.4))') ik,ths(ik),
     >                    oldti(ik),oldne(ik),nhs(ik),nh2s(ik),
     >                    sigvcx,qidcx(ik),ionsrc(ik)
               endif
c
            end do
c
         elseif (actswqidcx.eq.3) then
c
            write(year,'(i2.2)') iyearh
            call xxuid(useridh)
c
            ICLASS = 3
            CALL ADASRD(YEAR,1,1,ICLASS,nptscopy,ptesa,pnesa,PCOEF)
            DO IK = 1, nptscopy
               if (oldti(ik).lt.tcutcx) then
                  qidcx(ik) = 0.0
c                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
               else
                  QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))
     >               * econv * oldne(ik) * zb * nhs(ik) * PCOEF(IK,1)
               endif
            ENDDO
c
         elseif (actswqidcx.eq.4) then
c
            do ik = 1,nptscopy
c
               eav = 1.5 * (ths(ik)+oldti(ik)) / 2.0
c
               if (eav.gt.1000.0) then
                  sigvcx = 1.0d-13
               else
                  sigvcx = 1.0d-14 * eav**(1.0d0/3.0d0)
               endif
c
               if (oldti(ik).lt.tcutcx) then
                  qidcx(ik) = 0.0
c                  write (6,100) 'CX:',ik,oldti(ik),tcutcx
               else
                  QIDCX(IK) = 1.5 * (ths(ik) - oldti(ik))
     >                   * econv * oldne(ik) * zb * nhs(ik) * sigvcx
               endif
c
               if (qidcx(ik).gt.0.0) qidcx(ik) = 0.0
c
               if (m0.eq.initm0) then
c
                  if (ik.eq.1) write(6,*) 'QID Source terms:'
c
                  write (6,'(i4,8(1x,g11.4))') ik,ths(ik),
     >                    oldti(ik),oldne(ik),nhs(ik),nh2s(ik),
     >                    sigvcx,qidcx(ik),ionsrc(ik)
               endif
c
            end do


         endif
c
c        Now that all four components have been calculated - sum them
c        together and then integrate the total over S.
c
         do ik = 1,nptscopy
            qid(ik) = -1.0 *
     >           (qidatiz(ik) + qidmliz(ik) + qidrec(ik) + qidcx(ik))
         end do
c
c        Perform pre-integration into array intqid
c
c        This should give the equivalent of the PINQI array calculated
c        locally - the other two options to this routine access this
c        integrated array and return a linearly interpolated estimate
c        of the integral as a function of S.
c
         call preint(startn,nptscopy,sptscopy,qid,intqid,qidsum,
     >               ringlen,actswe2d,actswmajr,
     >               sbnd,rbnd,0.0d0)
c
         pinqid = qidsum
c
c        Print out a summary of results
c
c
         if (m0.eq.initm0) then
            write(6,'(a,g13.6,i4)') 'Sol option 22: QIDsrcint :',
     >                qidsum,ringnum
c
            do ik = startn,nptscopy
               write(6,'(i3,10(1x,g9.3))') ik,sptscopy(ik),
     >            qidatiz(ik),qidmliz(ik),qidrec(ik),qidcx(ik),
     >            -qid(ik),-intqid(ik),intqid(ik),qisrc(ik),qesrc(ik)
            end do
         end if
c
c
c
c     ESTIMATE or UPDATE  (at a later time the code for ESTIMATE and UPDATE
c                          may differ if these quantities are calculated
c                          for the current iteration - rather than the
c                          previous one.)
c
      elseif (opt.eq.0.or.opt.eq.1) then
c
c        Find which cell the particle is in
c
         call binsearch(s,in)
c
         if (in.eq.1) then
            pinqid = (intqid(in) * s / sptscopy(1))
         else
            pinqid = (
     >             ( intqid(in-1) +
     >             ( (intqid(in)-intqid(in-1)) *
     >             (s-sptscopy(in-1))/(sptscopy(in)-sptscopy(in-1)))))
         endif
c
      endif

c
c     Formatted OUTPUT
c
 100  format(a5,i4,2(1x,g13.5))
c
      return
      end
c
c
c
      subroutine detached_plasma(spts,npts,errcode,serr,
     >                           te,ti,ne,vb,ir,swtarg,
     >                           te0,ti0,n0,v0,act_press,mb)
      implicit none
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
      integer errcode,npts,ir
      real    swtarg
c
      include 'solparams'
c
      real*8 serr,te0,ti0,n0,v0,mb
      real*8 spts(npts)
      real*8 te(npts),ti(npts),ne(npts),vb(npts),act_press(npts)
c
c     Local variables
c
      real*8 smax,pmax,news,lppa,lprad
      real*8 t1,ti1,t2,n1,v1,s,sl1,sl2,slv,na
      real*8 sl1a,sl1b
c
      double precision l1r,l1ri,l2r,l2ri,lvr,lvri,
     >                 ter,teri,tir,tiri,nr,nri,vbm,vbmi,qr,qri,
     >                 n_exp,n_expi
c
c     Extra parameters - have no effect unless specified in 
c     for S21            additional per ring data sets.
c
c     l1a ... allow for more detail in the description of the 
c              density evolution in region A - three linear 
c              fitted regions. 
c
      double precision l1a,l1ai,l1b,l1bi,nr1a,nr1ai,nr1b,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi
c
      integer ik
c
      serr = 0.0
      errcode = 0
c
c     Assign lengths
c
      smax = ksmaxs(ir)
      pmax = kpmaxs(ir)
c
c     Set up scaling lengths
c
      call load_s21params(ir,l1r,l1ri,l2r,l2ri,
     >                    lvr,lvri,ter,teri,tir,tiri,nr,nri,
     >                    vbm,vbmi,qr,qri,n_exp,n_expi,
     >                         l1a,nr1a,l1ai,nr1ai,
     >                         l1b,nr1b,l1bi,nr1bi,
     >                         ter1a,ter1ai,tir1a,tir1ai,
     >                         ter1b,ter1bi,tir1b,tir1bi)
c
c     Note : OUTER target values are loaded by default - only overwritten
c            if INNER target is being calculated. 
c
      if (swtarg.eq.1.0) then
c
c        Outer target - ik = 1
c
c         l1r = l1r
c         l2r = l2r
c         lvr = lvr
c
c         tir = tir
c         ter = ter
c
c        Pressure matching not allowed for SOL 21 option included in
c        SOL 22.
c
c         nr  = abs(nr)
c
         na  = n_exp
c
c
c         qr  = qr
c
      elseif (swtarg.eq.2.0) then
c
c        Inner target - ik = nks(ir)
c
         l1a = l1ai
         l1b = l1bi

         l1r = l1ri
         l2r = l2ri
         lvr = lvri

         tir = tiri
         ter = teri
c
         ter1a = ter1ai 
         ter1b = ter1bi 
         tir1a = tir1ai 
         tir1b = tir1bi 
c
c
c        Pressure matching not allowed for SOL 21 option included in
c        SOL 22.
c
         nr  = abs(nri)
         na  = n_expi

         nr1a = abs(nr1ai)
         nr1b = abs(nr1bi)

         qr  = qri

      endif
c
c     Length reference switch for SOL 21 - options
c
c        0 = relative to SMAX
c        1 = relative to PMAX
c        2 = absolute SMAX units
c        3 = absolute PMAX units
c
      if (s21refsw.eq.0) then
c
         sl1a = l1a * smax
         sl1b = l1b * smax
         sl1  = l1r * smax
         sl2  = l2r * smax
         slv  = lvr * smax
c
      elseif (s21refsw.eq.1) then
c
         sl1a = l1a * pmax
         sl1b = l1b * pmax
         sl1  = l1r * pmax
         sl2  = l2r * pmax
         slv  = lvr * pmax
c
      elseif (s21refsw.eq.2.or.s21refsw.eq.3) then
c
         sl1a = l1a 
         sl1b = l1b 
         sl1  = l1r
         sl2  = l2r
         slv  = lvr
c
      endif
c
c     Convert P-coordinates to S-coordinates
c
      if (s21refsw.eq.1.or.s21refsw.eq.3) then
c
         call cnvrtptos(sl1a,news,ir)
         sl1a = news
c
         call cnvrtptos(sl1b,news,ir)
         sl1b = news
c
         call cnvrtptos(sl1,news,ir)
         sl1 = news
c
         call cnvrtptos(sl2,news,ir)
         sl2 = news
c
         call cnvrtptos(slv,news,ir)
         slv = news
c
      endif
c
      lppa = (5.0+2.0*ti0/te0+15.0/te0)
     >             *n0*abs(v0)*te0*1.602192e-19
c
      lprad = qrat * lppa / (sl2-sl1)
c
      t1    = ter * te0
c
      ti1   = tir * ti0
c
      n1    = nr * n0
c
      t2    = (t1**3.5+7.0/(2.0*ck0)*(lppa*(sl2-sl1)
     >                    +0.5*(sl2-sl1)**2*lprad))**(2.0/7.0)
c
      v1    = (n0 * v0) / n1
c
c      write (6,*) 'SOL22/21o:',sl1,sl2,slv,smax,lppa,lprad
c      write (6,*) 'SOL22/21b:',t1,ti1,n1,t2,v1
c
c     Calculate the background for the given S-values
c
      do ik = 1,npts
c
         s = spts(ik)
c
c        Test S-value and then calculate Te,Ti,N and V
c
         if (s.lt.sl1a) then
c
            te(ik) = te0 + (te0*ter1a-te0) * (s/sl1a)
c
            ti(ik) = ti0 + (ti0*tir1a-ti0) * (s/sl1a)
c
            ne(ik)  = n0 + (n0*nr1a-n0)  * (s/sl1a)**na
c
            vb(ik)  = (n0 * v0)/ ne(ik)
c
         elseif (s.lt.sl1b) then
c
            te(ik) = te0*ter1a + (te0*ter1b-te0*ter1a) 
     >                          * (s-sl1a)/(sl1b-sl1a)
c
            ti(ik) = ti0*tir1a + (ti0*tir1b-ti0*tir1a)
     >                          * (s-sl1a)/(sl1b-sl1a)
c
            ne(ik)  = n0*nr1a + (n0*nr1b-n0*nr1a)  
     >                          * ((s-sl1a)/(sl1b-sl1a))**na
c
            vb(ik)  = (n0 * v0)/ ne(ik)
c
         elseif (s.le.sl1) then
c
            te(ik) = te0*ter1b + (t1-te0*ter1b) 
     >                          * (s-sl1b)/(sl1-sl1b)
c
            ti(ik) = ti0*tir1b + (ti1-ti0*tir1b)
     >                          * (s-sl1b)/(sl1-sl1b)
c
c            ne(ik)  = n0 + (n1-n0)  * (s/sl1)**na
c
            ne(ik)  = n0*nr1b + (n1-n0*nr1b)  
     >                          * ((s-sl1b)/(sl1-sl1b))**na
c
            vb(ik)  = (n0 * v0)/ ne(ik)
c
         elseif (s.le.sl2) then
c
            te(ik) = (t1**3.5+7.0/(2.0*ck0)*(lppa*(s-sl1)
     >                    +0.5*(s-sl1)**2*lprad))**(2.0/7.0)
c
            ti(ik) = te(ik)
c
            ne(ik)  = (n1*t1) / te(ik)
c
            if (s.le.slv) then
               vb(ik)= v1 * (slv-s)/(slv-sl1)
            else
               vb(ik)= 0.0
            endif
c
         else
c
            te(ik) = (t2**3.5+7.0/(2.0*ck0)*
     >                    ((1.0+qr)*lppa*(s-sl2)))**(2.0/7.0)
c
            ti(ik) = te(ik)
c
            ne(ik)  = (n1*t1) / te(ik)
c
            if (s.le.slv) then
               vb(ik)= v1 * (slv-s)/(slv-sl1)
            else
               vb(ik)= 0.0
            endif
c
         endif
c
c        Calculate pressure
c
         act_press(ik) = ne(ik) * econv * (te(ik) + ti(ik))
     >               + mb * mconv * ne(ik) * vb(ik)**2
c
      end do
c
c     End of Detached Plasma calculation
c
      return
      end
c
c
c
      subroutine ioniz_comp
      implicit none
c
      include 'params'
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
      include 'printopt'
c
c     IONIZ_COMP: This routine prints a comparison of the
c                 ionization data read in from EDGE2D and
c                 the ionization calculated by NIMBUS.
c
c
      integer ir,ik,in
c
      real ioniz(maxnrs,3,21)
      real sumiz(9,21),sum1,sum2
c
      call rzero(ioniz,maxnrs*3*21)
      call rzero(sumiz,9*21)
      sum1 = 0.0
      sum2 = 0.0
c
c     Sum up over rings
c
      write(6,*)
      write(6,*)  'RAW DATA:'
      write(6,*)
      write(6,*)  ' IK  IR   PINION  E2DION   PIN-AREA  E2D-AREA'

c
      DO IR=1,nrs
c
         do ik = 1,nks(ir)
c

             write (6,'(2(i4,1x),5g12.4)') ik,ir,pinion(ik,ir),
     >                 e2dion(ik,ir),karea2(ik,ir),e2dareas(ik,ir)
c
             if (.not.(ir.lt.irsep.and.ik.eq.nks(ir))) then
c
             if (ik.le.ikmids(ir)) then
                ioniz(ir,1,1) = ioniz(ir,1,1) +
     >                          karea2(ik,ir) * pinion(ik,ir)
                ioniz(ir,1,2) = ioniz(ir,1,2) +
     >                          karea2(ik,ir) * pinion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,3) = ioniz(ir,1,3) +
     >                          karea2(ik,ir) * pinion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,4) = ioniz(ir,1,4) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
                ioniz(ir,1,5) = ioniz(ir,1,5) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,6) = ioniz(ir,1,6) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,7) = ioniz(ir,1,7) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
                ioniz(ir,1,8) = ioniz(ir,1,8) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,9) = ioniz(ir,1,9) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,10) = ioniz(ir,1,10) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
                ioniz(ir,1,11) = ioniz(ir,1,11) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,12) = ioniz(ir,1,12) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,13) = ioniz(ir,1,13) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
                ioniz(ir,1,14) = ioniz(ir,1,14) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,15) = ioniz(ir,1,15) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,16) = ioniz(ir,1,16) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
                ioniz(ir,1,17) = ioniz(ir,1,17) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,18) = ioniz(ir,1,18) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,1,19) = ioniz(ir,1,19) +
     >                  karea2(ik,ir) * pinion(ik,ir) / hcorr(ik,ir)
                ioniz(ir,1,20) = ioniz(ir,1,20) +
     >                  karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,1,21) = ioniz(ir,1,21) +
     >                  karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)
     >                          *rs(ik,ir)/r0
c
c
             else
c
                ioniz(ir,2,1) = ioniz(ir,2,1) +
     >                          karea2(ik,ir) * pinion(ik,ir)
                ioniz(ir,2,2) = ioniz(ir,2,2) +
     >                          karea2(ik,ir) * pinion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,3) = ioniz(ir,2,3) +
     >                          karea2(ik,ir) * pinion(ik,ir)
     >                          *rs(ik,ir)/r0

c
                ioniz(ir,2,4) = ioniz(ir,2,4) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
                ioniz(ir,2,5) = ioniz(ir,2,5) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,6) = ioniz(ir,2,6) +
     >                          karea2(ik,ir) * e2dion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,2,7) = ioniz(ir,2,7) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
                ioniz(ir,2,8) = ioniz(ir,2,8) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,9) = ioniz(ir,2,9) +
     >                          karea2(ik,ir) * pinrec(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,2,10) = ioniz(ir,2,10) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
                ioniz(ir,2,11) = ioniz(ir,2,11) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,12) = ioniz(ir,2,12) +
     >                          e2dareas(ik,ir) * pinion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,2,13) = ioniz(ir,2,13) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
                ioniz(ir,2,14) = ioniz(ir,2,14) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,15) = ioniz(ir,2,15) +
     >                          e2dareas(ik,ir) * e2dion(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,2,16) = ioniz(ir,2,16) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
                ioniz(ir,2,17) = ioniz(ir,2,17) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,18) = ioniz(ir,2,18) +
     >                          e2dareas(ik,ir) * pinrec(ik,ir)
     >                          *rs(ik,ir)/r0
c
                ioniz(ir,2,19) = ioniz(ir,2,19) +
     >                  karea2(ik,ir) * pinion(ik,ir) / hcorr(ik,ir)
                ioniz(ir,2,20) = ioniz(ir,2,20) +
     >                  karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)
     >                          *r0/rs(ik,ir)
                ioniz(ir,2,21) = ioniz(ir,2,21) +
     >                  karea2(ik,ir) * pinion(ik,ir) /hcorr(ik,ir)
     >                          *rs(ik,ir)/r0
c
             endif
c
             endif
c
          end do
      end do
c
c     Summations
c
      do ir = 1,nrs
c
         do in = 1,21
            ioniz(ir,3,in) = ioniz(ir,1,in) + ioniz(ir,2,in)
         end do
c
c        Core totals
c
         if (ir.lt.irsep) then
            if (ir.ge.ircent) then
               do in = 1,21
c
                  sumiz(1,in) = sumiz(1,in) + ioniz(ir,3,in)
c
               end do
            else
               do in = 1,21
c
                  sumiz(9,in) = sumiz(9,in) + ioniz(ir,3,in)
c
               end do
            endif
c
c        Main SOL totals
c
         elseif (ir.ge.irsep.and.ir.lt.irwall) then
c
            do in = 1,21
c
               sumiz(2,in) = sumiz(2,in) + ioniz(ir,3,in)
               sumiz(3,in) = sumiz(3,in) + ioniz(ir,1,in)
               sumiz(4,in) = sumiz(4,in) + ioniz(ir,2,in)
c
            end do
c
c        PP totals
c
         elseif (ir.gt.irtrap) then
c
            do in = 1,21
c
               sumiz(5,in) = sumiz(5,in) + ioniz(ir,3,in)
               sumiz(6,in) = sumiz(6,in) + ioniz(ir,1,in)
               sumiz(7,in) = sumiz(7,in) + ioniz(ir,2,in)
c
            end do
c
         endif

      end do
c
c     Grand totals
c
      do in = 1,21
         sumiz(8,in) = sumiz(1,in) + sumiz(2,in) + sumiz(5,in)
      end do
c
c     Print out summaries for ALL
c
      write (6,*)
      write (6,*) 'Ionization SUMMARY - PIN Areas:'
      write (6,*)
      write (6,100)
c
      write (6,110) 'PIN     :',(sumiz(in,1),in=1,8),1.0
      write (6,110) 'PIN R0/R:',(sumiz(in,2),in=1,8),
     >       sumiz(8,2)/sumiz(8,1)
      write (6,110) 'PIN R/R0:',(sumiz(in,3),in=1,8),
     >       sumiz(8,3)/sumiz(8,1)
      write (6,110) 'E2D     :',(sumiz(in,4),in=1,8),
     >        sumiz(8,4)/sumiz(8,1),sumiz(8,4)/sumiz(8,2),
     >        sumiz(8,4)/sumiz(8,3)
      write (6,110) 'E2D R0/R:',(sumiz(in,5),in=1,8),
     >        sumiz(8,5)/sumiz(8,1),sumiz(8,5)/sumiz(8,2),
     >        sumiz(8,5)/sumiz(8,3)
      write (6,110) 'E2D R/R0:',(sumiz(in,6),in=1,8),
     >        sumiz(8,6)/sumiz(8,1),sumiz(8,6)/sumiz(8,2),
     >        sumiz(8,6)/sumiz(8,3)
      write (6,110) 'REC     :',(sumiz(in,7),in=1,8),1.0
      write (6,110) 'REC R0/R:',(sumiz(in,8),in=1,8),
     >       sumiz(8,8)/sumiz(8,7)
      write (6,110) 'REC R/R0:',(sumiz(in,9),in=1,8),
     >       sumiz(8,9)/sumiz(8,7)
      write(6,*)
c
      write (6,*)
      write (6,*) 'Ionization SUMMARY - E2D Areas - Estimated:'
      write (6,*)
      write (6,100)
c
      write (6,110) 'PIN     :',(sumiz(in,10),in=1,8),1.0
      write (6,110) 'PIN R0/R:',(sumiz(in,11),in=1,8),
     >       sumiz(8,11)/sumiz(8,10)
      write (6,110) 'PIN R/R0:',(sumiz(in,12),in=1,8),
     >       sumiz(8,12)/sumiz(8,10)

      write (6,110) 'E2D     :',(sumiz(in,13),in=1,8),
     >        sumiz(8,13)/sumiz(8,10),sumiz(8,13)/sumiz(8,11),
     >        sumiz(8,13)/sumiz(8,12)
      write (6,110) 'E2D R0/R:',(sumiz(in,14),in=1,8),
     >        sumiz(8,14)/sumiz(8,10),sumiz(8,14)/sumiz(8,11),
     >        sumiz(8,14)/sumiz(8,12)
      write (6,110) 'E2D R/R0:',(sumiz(in,15),in=1,8),
     >        sumiz(8,15)/sumiz(8,10),sumiz(8,15)/sumiz(8,11),
     >        sumiz(8,15)/sumiz(8,12)
      write (6,110) 'REC     :',(sumiz(in,16),in=1,8),1.0
      write (6,110) 'REC R0/R:',(sumiz(in,17),in=1,8),
     >       sumiz(8,17)/sumiz(8,16)
      write (6,110) 'REC R/R0:',(sumiz(in,18),in=1,8),
     >       sumiz(8,18)/sumiz(8,16)
      write(6,*)
      write (6,*) 'Ionization SUMMARY - PIN Areas / HCORR :'
      write (6,*)
      write (6,100)

      write (6,110) 'PIN     :',(sumiz(in,19),in=1,8),1.0
      write (6,110) 'PIN R0/R:',(sumiz(in,20),in=1,8),
     >       sumiz(8,20)/sumiz(8,19)
      write (6,110) 'PIN R/R0:',(sumiz(in,21),in=1,8),
     >       sumiz(8,21)/sumiz(8,19)
      write (6,110) 'E2D     :',(sumiz(in,13),in=1,8),
     >        sumiz(8,13)/sumiz(8,19),sumiz(8,13)/sumiz(8,20),
     >        sumiz(8,13)/sumiz(8,21)
      write (6,110) 'E2D R0/R:',(sumiz(in,14),in=1,8),
     >        sumiz(8,14)/sumiz(8,19),sumiz(8,14)/sumiz(8,20),
     >        sumiz(8,14)/sumiz(8,21)
      write (6,110) 'E2D R/R0:',(sumiz(in,15),in=1,8),
     >        sumiz(8,15)/sumiz(8,19),sumiz(8,15)/sumiz(8,20),
     >        sumiz(8,15)/sumiz(8,21)
c
c
      write (6,*)
      write (6,'(a,1x,3g12.4)') 'PIN I-Core Ioniz (D-A):',sumiz(9,1),
     >                          sumiz(9,2),sumiz(9,3)
      write (6,'(a,1x,3g12.4)') 'PIN I-Core Recom (D-A):',sumiz(9,7),
     >                          sumiz(9,8),sumiz(9,9)
      write (6,'(a,1x,3g12.4)') 'PIN I-Core Ioniz (D-E):',sumiz(9,10),
     >                          sumiz(9,11),sumiz(9,12)
      write (6,'(a,1x,3g12.4)') 'PIN I-Core Recom (D-E):',sumiz(9,16),
     >                          sumiz(9,17),sumiz(9,18)
      write (6,'(a,1x,3g12.4)') 'E2D Inner core Ioniz:',sumiz(9,4)
c
      write (6,*)
      write (6,*) 'Grand Total Ionizations:'
      write (6,*)
      write (6,*) 'Polygon Areas:'
      write (6,*)
      write (6,110) 'PIN     :',sumiz(8,1) + sumiz(9,1)
      write (6,110) 'PIN R0/R:',sumiz(8,2) + sumiz(9,2)
      write (6,110) 'PIN R/R0:',sumiz(8,3) + sumiz(9,3)
      write (6,110) 'E2D     :',sumiz(8,4) + sumiz(9,4)
      write (6,110) 'E2D R0/R:',sumiz(8,5) + sumiz(9,5)
      write (6,110) 'E2D R/R0:',sumiz(8,6) + sumiz(9,6)
      write (6,*)
      write (6,*) 'E2D Areas:'
      write (6,*)
      write (6,110) 'PIN     :',sumiz(8,10) + sumiz(9,10)
      write (6,110) 'PIN R0/R:',sumiz(8,11) + sumiz(9,11)
      write (6,110) 'PIN R/R0:',sumiz(8,12) + sumiz(9,12)
      write (6,110) 'E2D     :',sumiz(8,13) + sumiz(9,13)
      write (6,110) 'E2D R0/R:',sumiz(8,14) + sumiz(9,14)
      write (6,110) 'E2D R/R0:',sumiz(8,15) + sumiz(9,15)
      write (6,*)


c
 100  format('TYPE',5x,3x,'O-Core',3x,2x,'SOL Total',1x,
     >       2x,'SOL Outer',1x,2x,'SOL Inner',1x,2x,'PP Total',2x,
     >       2x,'PP Outer',2x,2x,'PP Inner',2x,1x,'Grand Total',
     >       6x,'Ratios')

 110  format(a8,1x,8g12.4,3f9.5)

 120  format('TYPE',5x,4x,a5,3x,4x,a5,3x,4x,'Total',
     >        6x,'Ratio')


      write(6,*)

      write(6,*) 'Ionization Summary by ring'
      write(6,*)
c
      write(6,120) outer,inner
c
      sum1 = 1.0
      sum2 = 0.0
c
      do ir = 1,nrs
c
         write(6,*)
         write(6,*) 'RING = ',ir
         write(6,*)
         write (6,*)
         write (6,*) '     PIN/DIV - Areas'
         write (6,*)
         write (6,110) 'PIN     :',(ioniz(ir,in,1),in=1,3),1.0
         write (6,110) 'PIN R0/R:',(ioniz(ir,in,2),in=1,3),
     >     ioniz(ir,3,2)/ioniz(ir,3,1)
         write (6,110) 'PIN R/R0:',(ioniz(ir,in,3),in=1,3),
     >     ioniz(ir,3,3)/ioniz(ir,3,1)
         write (6,110) 'E2D     :',(ioniz(ir,in,4),in=1,3),
     >     ioniz(ir,3,4)/ioniz(ir,3,1),ioniz(ir,3,4)/ioniz(ir,3,2),
     >     ioniz(ir,3,4)/ioniz(ir,3,3)
         write (6,110) 'E2D R0/R:',(ioniz(ir,in,5),in=1,3),
     >     ioniz(ir,3,5)/ioniz(ir,3,1),ioniz(ir,3,5)/ioniz(ir,3,2),
     >     ioniz(ir,3,5)/ioniz(ir,3,3)
         write (6,110) 'E2D R/R0:',(ioniz(ir,in,6),in=1,3),
     >     ioniz(ir,3,6)/ioniz(ir,3,1),ioniz(ir,3,6)/ioniz(ir,3,2),
     >     ioniz(ir,3,6)/ioniz(ir,3,3)
         write (6,110) 'REC     :',(ioniz(ir,in,7),in=1,3),1.0
         write (6,110) 'REC R0/R:',(ioniz(ir,in,8),in=1,3),
     >     ioniz(ir,3,8)/ioniz(ir,3,7)
         write (6,110) 'REC R/R0:',(ioniz(ir,in,9),in=1,3),
     >     ioniz(ir,3,9)/ioniz(ir,3,7)
         write (6,*)
         write (6,*) '     E2D - Areas         IRCENT=',ircent
         write (6,*)
         if (ir.ge.ircent.and.ir.ne.irtrap.and.ir.ne.irwall) then
            sum1 = sum1 + ioniz(ir,3,10)
            sum2 = sum2 + ioniz(ir,3,13)
         endif

         write (6,110) 'PIN     :',(ioniz(ir,in,10),in=1,3),1.0
         write (6,110) 'PIN R0/R:',(ioniz(ir,in,11),in=1,3),
     >     ioniz(ir,3,11)/ioniz(ir,3,10)
         write (6,110) 'PIN R/R0:',(ioniz(ir,in,12),in=1,3),
     >     ioniz(ir,3,12)/ioniz(ir,3,10)
         write (6,110) 'E2D     :',(ioniz(ir,in,13),in=1,3),
     >     ioniz(ir,3,13)/ioniz(ir,3,10),ioniz(ir,3,13)/ioniz(ir,3,11),
     >     ioniz(ir,3,13)/ioniz(ir,3,12)
         write (6,110) 'E2D R0/R:',(ioniz(ir,in,14),in=1,3),
     >     ioniz(ir,3,14)/ioniz(ir,3,10),ioniz(ir,3,14)/ioniz(ir,3,11),
     >     ioniz(ir,3,14)/ioniz(ir,3,12)
         write (6,110) 'E2D R/R0:',(ioniz(ir,in,15),in=1,3),
     >     ioniz(ir,3,15)/ioniz(ir,3,10),ioniz(ir,3,15)/ioniz(ir,3,11),
     >     ioniz(ir,3,15)/ioniz(ir,3,12)
         write (6,110) 'REC     :',(ioniz(ir,in,16),in=1,3),1.0
         write (6,110) 'REC R0/R:',(ioniz(ir,in,17),in=1,3),
     >     ioniz(ir,3,17)/ioniz(ir,3,16)
         write (6,110) 'REC R/R0:',(ioniz(ir,in,18),in=1,3),
     >     ioniz(ir,3,18)/ioniz(ir,3,16)
         write(6,*)
         write (6,110) 'RUNTOT:',sum1,sum2,sum2/sum1
c
c         write (6,*)
c         write (6,*) '     PIN/HCORR - Areas        '
c         write (6,*)
c         write (6,110) 'PIN     :',(ioniz(ir,in,19),in=1,3),1.0
c         write (6,110) 'PIN R0/R:',(ioniz(ir,in,20),in=1,3),
c     >     ioniz(ir,3,20)/ioniz(ir,3,19)
c         write (6,110) 'PIN R/R0:',(ioniz(ir,in,21),in=1,3),
c     >     ioniz(ir,3,21)/ioniz(ir,3,19)
c         write (6,110) 'E2D     :',(ioniz(ir,in,13),in=1,3),
c     >     ioniz(ir,3,13)/ioniz(ir,3,19),ioniz(ir,3,13)/ioniz(ir,3,20),
c     >     ioniz(ir,3,13)/ioniz(ir,3,21)
c         write (6,110) 'E2D R0/R:',(ioniz(ir,in,14),in=1,3),
c     >     ioniz(ir,3,14)/ioniz(ir,3,19),ioniz(ir,3,14)/ioniz(ir,3,20),
c     >     ioniz(ir,3,14)/ioniz(ir,3,21)
c         write (6,110) 'E2D R/R0:',(ioniz(ir,in,15),in=1,3),
c     >     ioniz(ir,3,15)/ioniz(ir,3,19),ioniz(ir,3,15)/ioniz(ir,3,20),
c     >     ioniz(ir,3,15)/ioniz(ir,3,21)
c         write(6,*)
c
c
      end do
c
      write (6,*)
      write (6,*) 'EDGE 2D without boundary rings:',
     >            sumiz(8,4)-ioniz(1,3,4)-ioniz(irwall,3,4)
     >                                   -ioniz(irtrap,3,4)
      write (6,*)
      write (6,*) 'Total Ratios - PIN/DIV Areas:'
      write (6,*)
      write (6,'(a,g12.4)')
     >            'E2D  /DIVIMP = ', sumiz(8,4)/sumiz(8,1)
      write (6,'(a,g12.4)')
     >            'E2DR1/DIVIMP = ', sumiz(8,5)/sumiz(8,1)
      write (6,'(a,g12.4)')
     >            'E2DR2/DIVIMP = ', sumiz(8,6)/sumiz(8,1)
      write(6,*)
      write (6,*)
      write (6,*) 'Total Ratios - E2D Areas:'
      write (6,*)
      write (6,'(a,g12.4)')
     >            'E2D  /DIVIMP = ', sumiz(8,13)/sumiz(8,10)
      write (6,'(a,g12.4)')
     >            'E2DR1/DIVIMP = ', sumiz(8,14)/sumiz(8,10)
      write (6,'(a,g12.4)')
     >            'E2DR2/DIVIMP = ', sumiz(8,15)/sumiz(8,10)
      write(6,*)
c
c      write (6,*) 'E2D Areas and correction factors:'
c
c      write (6,*) ' IK   IR   KAREAS     E2D-AREAS   H(K)*DR*DT/R'//
c     >            '   H(K)'//
c     >            '     HCORR(K)     '
c      do ir = 1,nrs
c         do ik = 1,nks(ir)
c            write (6,'(i4,1x,i4,1x,8g12.4)') ik,ir,karea2(ik,ir),
c     >         e2dareas(ik,ir),karea2(ik,ir)/hcorr(ik,ir),
c     >         hval(ik,ir)*abs(drho(ik,ir)*dthetag(ik,ir))/rs(ik,ir),
c     >         hval(ik,ir),hcorr(ik,ir)
c         end do
c      end do
c
      write (6,*)
c
c
      return
      end
c
c
c
      subroutine prerrdesc(errlevel)
      implicit none
      real errlevel
c
c     PRERRDESC: Prints a description of the ERROR correction
c                protocol used by the solver.
c
c
         call prb
c
         CALL PRC('  ERROR CORRECTION IN SOLVER IS TURNED ON')
c
         CALL PRI('    - INITIAL ERROR CORRECTION LEVEL SET TO : ',
     >                             NINT(ERRLEVEL))
C
         CALL PRC('    - ERROR CORRECTION MOVES FROM HIGHEST TO'//
     >            ' LOWEST LEVEL')
         CALL PRC('      UNTIL SOLUTION IS OBTAINED.')
c
         call prc('    - ALL CHANGES ARE CUMULATIVE')
C
         call prb

c
c              errlevel=10- turn off equipartition if it was activated
c              errlevel=9 - use uniform particles instead of d2n/dr2
c              errlevel=8 - use only PINQI cooling contributions
c              errlevel=7 - use 1/2 ring uniform power instead of whole
c              errlevel=6 - use 1/2 ring uniform power + 1/2 ring particles
c              errlevel=5 - 1/2 ring uniform particles + power at top
c              errlevel=4 - 5 + turn off v^2 convection term
c              errlevel=3 - 4 + no power terms
c              errlevel=2 - 3 + no convective terms
c              errlevel=1 - Conduction ONLY
c

         CALL PRC('  LISTING OF ERROR CORRECTION LEVELS:')
c

         CALL PRC('  10 - TURN OFF EQUIPARTITION IF IT IS ON')
         CALL PRC('   9 - REPLACE DENSITY GRADIENT DEPENDENT'//
     >                  ' CROSS-FIELD TERM WITH UNIFORM')
         call prc('   8 - NO HEATING BY PINQI IS ALLOWED.')
         CALL PRC('   7 - REPLACE WHOLE RING UNIFORM POWER WITH'//
     >                  ' HALF RING UNIFORM.')
         CALL PRC('   6 - HALF RING UNIFORM POWER AND HALF RING UNIFORM'
     >                   //' PARTICLES')
         CALL PRC('   5 - HALF RING UNIFORM PARTICLES'//
     >            ' AND POWER IN AT TOP')
         call prc('   4 - 1/2 M V^3 CONVECTIVE TERM TURNED OFF')
         CALL PRC('   3 - ALL ADDITIONAL POWER TERMS TURNED OFF')
         call prc('   2 - ALL CONVECTIVE TERMS TURNED OFF')
         CALL PRC('   1 - CONDUCTION ONLY - ANALYTIC IONIZATION ONLY.')
         call prb

      return
      end

c
c
c
      subroutine print_edge2d_flux_analysis(rconst,gtarg,areasum,ir,
     >               gperpa,oldknbs,grad_oldknbs,srcinteg,recinteg,
     >               ike2d_start)
      implicit none
      include 'params'
      include 'solparams'
      include 'solswitch'
      integer ir,ike2d_start
      real    gperpa(maxnks,maxnrs),oldknbs(maxnks,maxnrs)
      real    grad_oldknbs(maxnks,maxnrs)
      real*8 rconst(mxspts),areasum(mxspts)
      real*8 gtarg(mxspts,3)
      real srcinteg (maxnks+3,maxnrs)
      real recinteg (maxnks+3,maxnrs)

c
      include 'comtor'
      include 'cgeom'
      include 'pindata'
      include 'cedge2d'
c
c     PRINT_EDGE2D_FLUX_ANALYSIS:
c
c
c     This routine prints a number of different analyses of the
c     EDGE2D fluxes in the context of SOL option 22.
c
c
c
c     Local variables
c
      integer ik,in,id,startik,endik
      real rmean1,rmean2
c
c     New variables for summaries
c
      real intgrad,intgradpos,intgradneg,inte2dsrc,intdist
      real inte2dhrec,intgrad2,ds,dcell,mulfact
      integer ikin,irin,ikout,irout
      real    influx,outflux,netflux,intnflux,flux1,flux2,flux3
      real    flux4
      real fluxst,fluxend,ioniz,rec,gperpd,gperpn,gdiv,sider
      real intgperpd,intgnet,intgdiv
      real initflux,dp,dp1,dp2,brat1,brat2,startflux
      real flux_const,endflux


c
         if (switch(swe2d).eq.0.0) then
            startik = 1
            endik   = nks(ir)
         elseif (switch(swe2d).gt.0.0) then
            startik = ike2d_start
            endik   = nks(ir) - ike2d_start + 1
         endif
c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
c     ------------------------------------------------------------------
c
            write(6,'(20a)') 'FLUX:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'F:',ir,ik,
     >              kss(ik,ir),ds,
     >              gtarg(ir,2),inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (ik.eq.endik) then
                  if (intgrad2.gt.0.0) then 
                     mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
                  else
                     mulfact = 0.0
                  endif
               endif
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do
c
            ik = nks(ir)
c
            flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
            flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
            flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
c
            write (6,'(a,5g12.5)') 'NUM:',gtarg(ir,2),
     >                0.15*intgrad2,flux1,flux2,mulfact
c
c     ------------------------------------------------------------------
c
            write(6,'(a,2i4,1p,8(1x,g13.5))') 'F:',ir,ik,
     >            gtarg(ir,2),inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >            gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec,
     >            e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),
     >            (gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec)/
     >            (e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),
     >            gtarg(ir,1)
c
            write (6,*)
c
            write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact
            write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact
c
            write(6,'(20a)') 'FLUX:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1 * mulfact
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'FLUX:',ir,ik,
     >              kss(ik,ir),ds,
     >              gtarg(ir,2),inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec,
     >              e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1,
     >              (gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec)/
     >              (e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1),
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2 * mulfact
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do

c
            write (6,*)
c




c     ------------------------------------------------------------------

c
c           Analysis of the actual EDGE2D boundary fluxes
c
c
            write (6,500)


500   format('BND:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,
     >        4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,
     >        1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,
     >        3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,
     >        3x,'D-RAT',4x,2x,'E2DV-RAT')

c
           intgperpd = 0.0
           intgnet   = 0.0
           intgdiv   = 0.0


            do ik = startik,endik
c
               fluxst  = e2dflux(ik,ir)
               fluxend = e2dflux(ik+1,ir)
c
               ioniz = e2dion(ik,ir)
c
               rec   = e2dhrec(ik,ir)
c
               ds    = ksb(ik,ir) - ksb(ik-1,ir)
c
               gperpn= -(fluxst - fluxend
     >               + e2dion(ik,ir)*ds - e2dhrec(ik,ir)*ds)
c
               gperpd= gperpn/ds
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
c
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               gdiv =  0.15*grad_oldknbs(ik,ir)
c
               intgperpd = intgperpd + ds * gperpd
               intgnet   = intgnet   + ds * netflux
               intgdiv   = intgdiv   + ds * gdiv
c
               write (6,'(a,2i4,1p,15(1x,g11.4))') 'BND:',ik,ir,
     >                fluxst,ioniz*ds,-rec*ds,gperpn,fluxend,gperpd,
     >                netflux,gdiv,intgperpd,intgnet,intgdiv,
     >                gdiv/gperpd,netflux/gperpd
c
            end do
c
c     ------------------------------------------------------------------
c

            write (6,*)
            write (6,*) 'DOWN FLUX ANALYSIS:  1/2 CELL'
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            write(6,'(20a)') 'DFH:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0

            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'DFH:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
               endif
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do
c
            ik = nks(ir)
c
            flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
            flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
            flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
            flux4 =  e2dgpara(ik+1,ir)
c
            write (6,'(a,5g12.5)') 'NUM:',startflux,
     >                0.15*intgrad2,flux1,flux2,mulfact
c
c     ------------------------------------------------------------------
c
            write(6,'(a,2i4,1p,8(1x,g13.5))') 'F:',ir,ik,
     >            startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >            startflux+inte2dsrc+0.15*intgrad2-inte2dhrec,
     >            e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),
     >            (startflux+inte2dsrc+0.15*intgrad2-inte2dhrec)/
     >            (e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),
     >            flux4
c
            write (6,*)
c
            write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact
            write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact
c
            write(6,'(20a)') 'MDFH:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1 * mulfact
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MDFH:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2 * mulfact
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do
c
            write (6,*)
c
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c
            write (6,*)
            write (6,*) 'DOWN FLUX ANALYSIS:  FULL CELL'
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            write(6,'(20a)') 'DFF:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            startflux = e2dgpara(1,ir)

            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'DFF:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
               endif
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
               endif
c
            end do
c
            ik = nks(ir)
c
            flux1 =  gtarg(ir,2)+inte2dsrc+0.15*intgrad2-inte2dhrec
            flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
            flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
            flux4 = e2dgpara(ik+1,ir)
c
            write (6,'(a,5g12.5)') 'NUM:',startflux,
     >                0.15*intgrad2,flux1,flux2,mulfact
c
c     ------------------------------------------------------------------
c
c
            write(6,'(a,2i4,1p,8(1x,g13.5))') 'DFF:',ir,ik,
     >            startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >            startflux+inte2dsrc+0.15*intgrad2-inte2dhrec,
     >            e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),
     >            (startflux+inte2dsrc+0.15*intgrad2-inte2dhrec)/
     >            (e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),
     >            flux4
c
            write (6,*)
c
            write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact
            write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact
c
            write(6,'(20a)') 'MDFF:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            startflux = e2dgpara(1,ir)
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1 * mulfact
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
c
               flux1 =  startflux+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MDFF:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2 * mulfact
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
            end do
c
            write (6,*)
c
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c
            write (6,*)
            write (6,*) 'FLUX ANALYSIS:  1/2 CELL - REGULAR DIST.'
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            write(6,'(20a)') 'GFH:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0
            endflux   = (e2dgpara(nks(ir),ir)
     >                     + e2dgpara(nks(ir)+1,ir))/2.0
c
            flux_const= (startflux - endflux + srcinteg(nks(ir)+1,ir)
     >                      - recinteg(nks(ir)+1,ir))
     >                 /areasum(ir)

c

            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + flux_const*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'GFH:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,intgrad,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              flux_const*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/intgrad) + 1.0
               endif
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + flux_const*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do
c
            ik = nks(ir)
c
            flux1 =  gtarg(ir,2)+inte2dsrc+intgrad-inte2dhrec
            flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
            flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
            flux4 =  e2dgpara(ik+1,ir)
c
            write (6,'(a,5g12.5)') 'NUM:',startflux,
     >                0.15*intgrad2,flux1,flux2,mulfact
c
            write (6,*)
c
c     ------------------------------------------------------------------
c
c
            write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact
            write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact
c
            write(6,'(20a)') 'MGFH:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            startflux = (e2dgpara(1,ir)+e2dgpara(2,ir))/2.0
            endflux   = (e2dgpara(nks(ir),ir)
     >                     + e2dgpara(nks(ir)+1,ir))/2.0
c
            flux_const= (startflux -endflux + srcinteg(nks(ir)+1,ir)
     >                      - recinteg(nks(ir)+1,ir))
     >                 /areasum(ir)

c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + flux_const*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1 * mulfact
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MGFH:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,intgrad,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              flux_const*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + flux_const*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2 * mulfact
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               endif

            end do
c
            write (6,*)
c

c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c


            write (6,*)
            write (6,*) 'FLUX ANALYSIS:  FULL CELL - REGULAR DIST.'
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            write(6,'(20a)') 'GFF:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            startflux = e2dgpara(1,ir)
            endflux   = e2dgpara(nks(ir)+1,ir)
c
            flux_const= (startflux  -endflux + srcinteg(nks(ir)+1,ir)
     >                      - recinteg(nks(ir)+1,ir))
     >                 /areasum(ir)


            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
                  intgrad = intgrad
     >               + flux_const*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'GFF:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,intgrad,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              flux_const*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/(intgrad)) + 1.0
               endif
c
                  intgrad = intgrad
     >               + flux_const*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
               if (ik.eq.endik) then
                  mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
               endif
c
            end do
c
            ik = nks(ir)
c
            flux1 =  gtarg(ir,2)+inte2dsrc+intgrad-inte2dhrec
            flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)
            flux3 =  gtarg(ir,2)+inte2dsrc+intnflux-inte2dhrec
            flux4 = e2dgpara(ik+1,ir)
c
            write (6,'(a,5g12.5)') 'NUM:',startflux,
     >                0.15*intgrad2,flux1,flux2,mulfact
c
c     ------------------------------------------------------------------
c
            write(6,'(a,2i4,1p,8(1x,g13.5))') 'GFF:',ir,ik,
     >            startflux,inte2dsrc,intgrad,-inte2dhrec,
     >            startflux+inte2dsrc+intgrad-inte2dhrec,
     >            e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir),
     >            (startflux+inte2dsrc+intgrad-inte2dhrec)/
     >            (e2dnbs(nks(ir),ir)*e2dvhs(nks(ir),ir)),
     >            flux4
c
            write (6,*)
c
            write (6,*) 'GPERP MULTIPLICATION FACTOR = ',mulfact
            write (6,*) 'EFFECTIVE DPERP             = ',0.15*mulfact
c
            write(6,'(20a)') 'MGFF:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            startflux = e2dgpara(1,ir)
            endflux   = e2dgpara(nks(ir)+1,ir)
c
            flux_const= (startflux -endflux + srcinteg(nks(ir)+1,ir)
     >                      - recinteg(nks(ir)+1,ir))
     >                 /areasum(ir)

c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
c
                  intgrad = intgrad
     >               + flux_const*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1 * mulfact
c
                  intnflux = intnflux
     >               + netflux*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1
c
                  intdist = intdist
     >               + (kss2(ik,ir)-ksb(ik-1,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(kss2(ik,ir)-ksb(ik-1,ir))
     >                * rmean1

                  endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
c
               flux1 =  startflux+inte2dsrc+intgrad-inte2dhrec
               flux2 =  (e2dgpara(ik,ir) + e2dgpara(ik+1,ir))/2.0
c               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1
               flux3 =  startflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,17(1x,g11.4))') 'MGFF:',ir,ik,
     >              kss(ik,ir),ds,
     >              startflux,inte2dsrc,intgrad,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              flux_const*ds,e2dion(ik,ir)*ds,
     >              e2dhrec(ik,ir)*ds,netflux*ds,
     >              flux3,flux3/flux2
c
                  intgrad = intgrad
     >               + flux_const*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2 * mulfact
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2
c
                  intdist = intdist
     >               + (ksb(ik,ir)-kss2(ik,ir))
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*(ksb(ik,ir)-kss2(ik,ir))
     >                * rmean2

                  endif
c
            end do
c
            write (6,*)
c
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c
c
c           Analysis of the actual EDGE2D boundary fluxes
c
c
            write (6,*) 'DOWN FLUX ANALYSIS:'
            write (6,510)


510   format('D-G:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,
     >        4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,
     >        1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,
     >        3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,
     >        3x,'D-RAT',4x,2x,'E2DV-RAT')

c
           intgperpd = 0.0
           intgnet   = 0.0
           intgdiv   = 0.0


            do ik = startik,endik
c
               fluxst  = e2dgpara(ik,ir)
               fluxend = e2dgpara(ik+1,ir)
c
               ioniz = e2dion(ik,ir)
c
               rec   = e2dhrec(ik,ir)
c
               ds    = ksb(ik,ir) - ksb(ik-1,ir)
c
               gperpn= -(fluxst - fluxend
     >               + e2dion(ik,ir)*ds - e2dhrec(ik,ir)*ds)
c
               gperpd= gperpn/ds
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
c
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               gdiv =  0.15*grad_oldknbs(ik,ir)
c
               intgperpd = intgperpd + ds * gperpd
               intgnet   = intgnet   + ds * netflux
               intgdiv   = intgdiv   + ds * gdiv
c
               write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,
     >                fluxst,ioniz*ds,-rec*ds,gperpn,fluxend,gperpd,
     >                netflux,gdiv,intgperpd,intgnet,intgdiv,
     >                gdiv/gperpd,netflux/gperpd
c
            end do
c
c

c
c           Analysis of the actual EDGE2D boundary fluxes
c
c
c     ------------------------------------------------------------------
c
            write (6,*) 'DOWN FLUX ACTUAL PARTICLE ANALYSIS: KAREAS'
            write (6,520)


520   format('DPGA:',1x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,
     >        4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,
     >        1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,
     >        3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,
     >        3x,'D-RAT',4x,2x,'E2DV-RAT')

c
           intgperpd = 0.0
           intgnet   = 0.0
           intgdiv   = 0.0


            do ik = startik,endik
c
               in = korpg(ik,ir)
c
               sider = (rvertp(1,in) + rvertp(2,in)) /2.0
               fluxst  = e2dgdown(ik,ir)/ (2.0* PI * sider)
c
               sider = (rvertp(3,in) + rvertp(4,in)) /2.0
               fluxend = e2dgdown(ik+1,ir)/ (2.0 * PI * sider)
c
               ioniz = e2dion(ik,ir) * kareas(ik,ir)
c
               rec   = e2dhrec(ik,ir) * kareas(ik,ir)
c
               ds    = ksb(ik,ir) - ksb(ik-1,ir)
c
               gperpn= -(fluxst - fluxend
     >               + ioniz - rec)
c
               gperpd= gperpn/kareas(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
c
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               gdiv =  0.15*grad_oldknbs(ik,ir)
c
               intgperpd = intgperpd + ds * gperpd
               intgnet   = intgnet   + ds * netflux
               intgdiv   = intgdiv   + ds * gdiv
c
               write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,
     >                fluxst,ioniz,-rec,gperpn,fluxend,gperpd,
     >                netflux,gdiv,intgperpd,intgnet,intgdiv,
     >                gdiv/gperpd,netflux/gperpd
c
            end do
c
c
c
c           Analysis of the actual EDGE2D boundary fluxes
c
c     ------------------------------------------------------------------
c
            write (6,*) 'DOWN FLUX ACTUAL PARTICLE ANALYSIS: E2DAREAS'
            write (6,530)


530   format('DGPE:',1x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,
     >        4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,
     >        1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,2x,'G-DIV-E',3x,
     >        3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,
     >        3x,'D-RAT',4x,2x,'E2DV-RAT',2x,3x,'GDIV-A')

c
           intgperpd = 0.0
           intgnet   = 0.0
           intgdiv   = 0.0


            do ik = startik,endik
c
               in = korpg(ik,ir)
c
               sider = (rvertp(1,in) + rvertp(2,in)) /2.0
               fluxst  = e2dgdown(ik,ir)/ (2.0* PI * sider)
c
               sider = (rvertp(3,in) + rvertp(4,in)) /2.0
               fluxend = e2dgdown(ik+1,ir)/ (2.0 * PI * sider)
c
               ioniz = e2dion(ik,ir) * e2dareas(ik,ir)
c
               rec   = e2dhrec(ik,ir) * e2dareas(ik,ir)
c
               ds    = ksb(ik,ir) - ksb(ik-1,ir)
c
               gperpn= -(fluxst - fluxend
     >               + ioniz - rec)
c
               gperpd= gperpn/e2dareas(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
c
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               gdiv =  0.15*grad_oldknbs(ik,ir)
     >                   *kareas(ik,ir)/e2dareas(ik,ir)
c
               intgperpd = intgperpd + ds * gperpd
               intgnet   = intgnet   + ds * netflux
               intgdiv   = intgdiv   + ds * gdiv
c
               write (6,'(a,2i4,1p,15(1x,g11.4))') 'D-G:',ik,ir,
     >                fluxst,ioniz,-rec,gperpn,fluxend,gperpd,
     >                netflux,gdiv,intgperpd,intgnet,intgdiv,
     >                gdiv/gperpd,netflux/gperpd,
     >                0.15*grad_oldknbs(ik,ir)
c
            end do
c
c
c     ------------------------------------------------------------------
c     ------------------------------------------------------------------
c

            intgrad = 0.0
            intgrad2 = 0.0
            intgradpos = 0.0
            intgradneg = 0.0
            inte2dsrc = 0.0
            inte2dhrec = 0.0
            intdist = 0.0
            intnflux = 0.0
c
            write(6,*) 'Poloidal plane analysis:'
            write(6,'(20a)') 'FLUX:','  IR ',' IK ',
     >              '      S    ','   DS   ',
     >              '   GTARG     ',' |E2DION   ','  |GPERP   ',
     >              '  |RECSRC     ','|NETFLUX ',
     >              '     SUM    ',
     >              '   E2D FLUX    ' , 'RATIO   ','    GPERPD  ',
     >              '    E2DION  ','    RECSRC   ','   NETFLUX   ',
     >              '   SUM2'

c
            initflux = gtarg(ir,2) / kbfst(ir,2)
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               dp = kpb(ik,ir)-kpb(ik-1,ir)
               dp1 = kps2(ik,ir) - kpb(ik-1,ir)
               dp2 = kpb(ik,ir) - kps(ik,ir)
c
               if (.not.(ik.eq.startik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >             switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*dp1
     >                * rmean1
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*dp1
     >                * rmean1
c
                  intnflux = intnflux
     >               + netflux*dp1
     >                * rmean1
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*dp1
     >                * rmean1
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*dp1
     >                * rmean1
c
                  intdist = intdist
     >               + dp1
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*dp1
     >                * rmean1

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*dp1
     >                * rmean1

                  endif
c
               endif
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = rs(ik,ir)
               else
                  rmean1 = 1.0
               endif
c
               flux1 =  initflux+inte2dsrc+0.15*intgrad2-inte2dhrec
               flux2 =  e2dnbs(ik,ir)*e2dvhs(ik,ir)*rmean1/kbfs(ik,ir)
               flux3 =  initflux+inte2dsrc+intnflux-inte2dhrec
c
               write(6,'(a,2i4,2f9.4,1p,15(1x,g11.4))') 'F:',ir,ik,
     >              kps(ik,ir),dp,
     >              initflux,inte2dsrc,0.15*intgrad2,-inte2dhrec,
     >              intnflux,
     >              flux1,
     >              flux2,
     >              flux1/flux2,
     >              0.15 * grad_oldknbs(ik,ir)*dp,e2dion(ik,ir)*dp,
     >              e2dhrec(ik,ir)*dp,netflux*dp,
     >              flux3,flux3/flux2
c
c               if (ik.eq.endik) then
c                 mulfact = ((flux2-flux1)/(0.15*intgrad2)) + 1.0
c               endif
c
               if (.not.(ik.eq.endik.and.
     >            (switch(swe2d).eq.1.0.or.switch(swe2d).eq.2.0.or.
     >              switch(swe2d).eq.3.or.switch(swe2d).eq.4.or.
     >             switch(swe2d).eq.6.or.switch(swe2d).eq.7.or.
     >             switch(swe2d).eq.8.or.switch(swe2d).eq.9))) then
c
c
                  intgrad = intgrad
     >               + gperpa(ik,ir)*dp2
     >                * rmean2
c
                  intgrad2 = intgrad2
     >               + grad_oldknbs(ik,ir)*dp2
     >                * rmean2
c
                  inte2dsrc = inte2dsrc
     >               + e2dion(ik,ir)*dp2
     >                * rmean2
c
                  inte2dhrec = inte2dhrec
     >               + e2dhrec(ik,ir)*dp2
     >                * rmean2
c
                  intnflux = intnflux
     >               + netflux*dp2
     >                * rmean2
c
                  intdist = intdist
     >               + dp2
c
                  if (gperpa(ik,ir).le.0.0) then
                     intgradneg = intgradneg
     >                 + gperpa(ik,ir)*dp2
     >                * rmean2

                  else
                     intgradpos = intgradpos
     >                 + gperpa(ik,ir)*dp2
     >                * rmean2

                  endif
c
               endif

            end do
c





c
c           Analysis of the actual EDGE2D boundary fluxes - poloidal plane
c
c
c     ------------------------------------------------------------------
c
            write (6,*) 'Poloidally calculated:'
            write (6,501)


501   format('BND:',2x,'IK',2x,'IR',1x,2x,'FLUX-ST',2x,3x,'IONIZ',4x,
     >        4x,'REC',5x,1x,'GPERP-REQ',2x,2x,'FLUX-END',2x,
     >        1x,'G-REQ-DEN',2x,2x,'G-VPERP',3x,3x,'G-DIV',4x,
     >        3x,'|G-REQ',3x,2x,'|G-VPERP',2x,3x,'|G-DIV',3x,
     >        3x,'D-RAT',4x,2x,'E2DV-RAT')

c

           intgperpd = 0.0
           intgnet   = 0.0
           intgdiv   = 0.0

c
            do ik = startik,endik
c
               if (ik.eq.startik) then
                  brat1 = kbfst(ir,2)
c slmod begin - bug fix, 04.03.2010
c... Correct correction, swapping KBFST for KBFS when IR is used to index the array?
                  brat2 = (kbfs(ik  ,ir)+kbfs(ik+1,ir))/2.0
               elseif (ik.eq.endik) then
                  brat1 = (kbfs(ik-1,ir)+kbfs(ik  ,ir))/2.0
                  brat2 = kbfst(ir,1)
               else
                  brat1 = (kbfs(ik-1,ir)+kbfs(ik  ,ir))/2.0
                  brat2 = (kbfs(ik  ,ir)+kbfs(ik+1,ir))/2.0
               endif
c
c                  brat2 = (kbfst(ik,ir)+kbfs(ik+1,ir))/2.0
c               elseif (ik.eq.endik) then
c                  brat1 = (kbfst(ik-1,ir)+kbfs(ik,ir))/2.0
c                  brat2 = kbfst(ir,1)
c               else
c                  brat1 = (kbfst(ik-1,ir)+kbfs(ik,ir))/2.0
c                  brat2 = (kbfst(ik,ir)+kbfs(ik+1,ir))/2.0
c               endif
c slmod end
               fluxst  = e2dflux(ik,ir) /brat1
               fluxend = e2dflux(ik+1,ir) /brat2
c
               ioniz = e2dion(ik,ir)
c
               rec   = e2dhrec(ik,ir)
c
               dp    = kpb(ik,ir) - kpb(ik-1,ir)
c
               gperpn= -(fluxst - fluxend
     >               + e2dion(ik,ir)*dp - e2dhrec(ik,ir)*dp)
c
               gperpd= gperpn/dp
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               dcell = distin(ik,ir) + distout(ik,ir)
c
               influx = e2dnbs(ik,ir) * e2dvro(ik,ir)
c
               outflux = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               netflux = (influx - outflux)/dcell
c
               gdiv =  0.15*grad_oldknbs(ik,ir)
c
               intgperpd = intgperpd + dp * gperpd
               intgnet   = intgnet   + dp * netflux
               intgdiv   = intgdiv   + dp * gdiv
c
               write (6,'(a,2i4,1p,15(1x,g11.4))') 'BND:',ik,ir,
     >                fluxst,ioniz*dp,-rec*dp,gperpn,fluxend,gperpd,
     >                netflux,gdiv,intgperpd,intgnet,intgdiv,
     >                gdiv/gperpd,netflux/gperpd
c
            end do
c
c
c     ------------------------------------------------------------------
c
            write (6,*)
c
c           Analyse the EDGE2D perpendicular fluxes and compare them
c           to the calculated values.
c
            write(6,*) '  FLUX:','  IR ',' IK ','      S     ',
     >               '      DS     ',
     >               '     FLUX1   ','     FLUX2   ','      FLUX3 ',
     >               '  (F1-F2)/dcell',
     >               ' (F2-F3)/dcell',
     >               '  DIVFLUXD   ', '   E2D NBS  ',
     >               '    VRO'
c
            do ik = startik,endik
c
               if (switch(swmajr).eq.4.0) then
                  rmean1 = (rs(ik,ir) + krb(ik-1,ir)) /2.0
                  rmean2 = (krb(ik,ir) + rs(ik,ir)) /2.0
               else
                  rmean1 = 1.0
                  rmean2 = 1.0
               endif
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               ds = ksb(ik,ir)-ksb(ik-1,ir)
               dcell = distin(ik,ir) + distout(ik,ir)
c
               flux1 = e2dnbs(ikin,irin) * e2dvro(ikin,irin)
               flux2 = e2dnbs(ik,ir) * e2dvro(ik,ir)
               flux3 = e2dnbs(ikout,irout) * e2dvro(ikout,irout)
c
               write(6,'(a,2i4,1p,15(1x,g12.4))') 'E2DFLUX:',ir,ik,
     >              kss(ik,ir),ds,flux1,flux2,flux3,(flux1-flux2)/dcell,
     >              (flux2-flux3)/dcell,
     >              0.15*grad_oldknbs(ik,ir),
     >              e2dnbs(ik,ir),e2dvro(ik,ir)
c
            end do
c
c     ------------------------------------------------------------------
c
c
      return
      end
