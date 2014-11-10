c     -*-Fortran-*-
c
      SUBROUTINE DIV (title,equil,NIZS,NIMPS,NIMPS2,CPULIM,IONTIM,
     >                NEUTIM,SEED,NYMFS,FACTA,FACTB,ITER,NRAND)
c
!
! ammod begin.
!      Use ComHC ! Contains DIVIMP input file hydrocarbon-related options.
!      Use HC_Start ! Contains HC begin and end routines.
!      Use HC_Batch ! Used to call HC_Launch for C+ self-sputtered particles.
!      Use HC_WBC_Comp ! Records ion death statistics for WBC comparison.
!      Use HC_Utilities ! Sheath E-field calc by Brooks.
c
      use error_handling
      use debug_options
      use eckstein_2007_yield_data
      use subgrid_options
      use subgrid
      use ero_interface
c slmod begin
      use mod_interface
      use mod_divimp
      use mod_divimp_tdep
c slmod end
c
      implicit none
c
      include 'params'

c
      character*(*) title,equil
      INTEGER  NIZS,NIMPS,NYMFS,ITER,NRAND,NIMPS2
      REAL     IONTIM,NEUTIM,CPULIM
      REAL     FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
      DOUBLE PRECISION SEED
C
C  *********************************************************************
C  *                                                                   *
C  *   DIV:  MAIN CONTROLLING ROUTINE                                  *
C  *   ------------------------------                                  *
C  *                                                                   *
C  *                                                                   *
C  *                        CHRIS FARRELL (HUNTERSKIL)  FEB 1989       *
C  *                                                                   *
C  *                        JAMES SPENCE  (TESSELLA)    NOV 1990       *
C  *                                                                   *
C  *                        DAVID ELDER   (UTIAS)       JAN 1992       *
C  *                                                                   *
C  *********************************************************************
C
      include    'dynam1'
      include    'dynam3'
      include    'dynam4'
      include    'comtor'
      include    'cgeom'
      include    'cioniz'
      include    'commv'
      include    'cneut'
      include    'cneut2'
      include    'cnoco'
      include    'cadas'
      include    'clocal'
      include    'crand'
      include    'pindata'
      include    'diagvel'
      include    'cedge2d'
      include    'promptdep'
      include    'reiser_com'
      include    'printopt'
c
      include    'fperiph_com'
c
c      include    'div_com'
c
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div4'
      include 'div5'
      include 'div6'
      include 'div7'
c
      include    'particle_specs'
      include    'driftvel'
      include    'hc_global_opts'
c slmod begin - temp
      include 'slcom'

      integer i,fp,load_i
c slmod end


      integer :: perc
c
c     Output velocity along the field line from launch_one
c
      real vout


      integer fperiph
      external fperiph

      integer   verify_id
      external  verify_id

      real za02as
      external za02as

      real      yield
      external  yield

      integer ipos
      external ipos

      real ndrand
      external ndrand

      character*77 comment

      logical :: ero_record_data

      !
      ! Spara is now local to the transport routines - it is not
      ! available for printing in the remaining debug statements at
      ! this level - the quick workaround is to replace it with
      ! a value that is always zero.
      !
      real zero_spara


      CHARACTER WHAT(51)*10,FATE(11)*14,FACTOR*9,C(10)*9
c
      DATA  FATE  /'REACHED WALL',        'REACHED TARGET',
     >             'REACHED TIME',        'STARTED > NIZS',
     >             'RECOMBINED: 0',       'IONISED > NIZS',
     >             'REACHED MAIN',        'REMOVED       ',
     >             'HIT F-P TARGET',      'PROMPT REDEP  ',
     >             'REFLECTED NEUT'/
C
      DATA  WHAT  /'  PRIMARY ',   ' SECONDARY',   ' TERTIARY ',
     >             'QUATERNARY',   '  QUINARY ',   '  SIXTH   ',
     >             '  SEVENTH ',   '  EIGHTH  ',   '  NINTH   ',
     >             '  TENTH   ',   ' ELEVENTH ',   ' TWELFTH  ',
     >             'THIRTEENTH',   'FOURTEENTH',   ' FIFTEENTH',
     >             ' SIXTEENTH',   ' SEVENTEEN',   'EIGHTEENTH',
     >             'NINETEENTH',   ' TWENTIETH',   'TWENTY-1ST',
     >             'TWENTY-2ND',   'TWENTY-3RD',   'TWENTY-4TH',
     >             'TWENTY-5TH',   'TWENTY-6TH',   'TWENTY-7TH',
     >             'TWENTY-8TH',   'TWENTY-9TH',   ' THIRTIETH',
     >             'THIRTY-1ST',   'THIRTY-2ND',   'THIRTY-3RD',
     >             'THIRTY-4TH',   'THIRTY-5TH',   'THIRTY-6TH',
     >             'THIRTY-7TH',   'THIRTY-8TH',   'THIRTY-9TH',
     >             ' FORTIETH ',   ' FORTY-1ST',   ' FORTY-2ND',
     >             ' FORTY-3RD',   ' FORTY-4TH',   ' FORTY-5TH',
     >             ' FORTY-6TH',   ' FORTY-7TH',   ' FORTY-8TH',
     >             ' FORTY-9TH',   ' FIFTIETH ',
     >             '    ALL   '/



C
C-----------------------------------------------------------------------
C                   INITIALISATION
C-----------------------------------------------------------------------
C
      !call toggle_trace

      call pr_trace('DIV','BEGIN DIV')
c
C psmod
c
      debug0 = .false.
      debug_all=.false.
c
c
      Ftotal = 0.0
      COPTION = 0
c
c psmod
c
c      cisterrcnt = 0
c
c     The following "leakage" variables are used to monitor leakage
c     along the field lines away from the targets.
c
      hasleaked = .false.
      cleakt    =  0.0
      cleakp    =  1
c
c     Initialize extra source total particle influx  [#/m/s]
c
      neut2d_fytot = 0.0
c
c     If the ERO option is ON and neut2d_opt=2 calling for an ERO launch of particles 
c
      if (neut2d_opt.eq.2) then 

         if (ero_particle_launch_opt.gt.0) then 
            ! load and analyse the ero launch data in preparation for 
            ! particle launch in NEUT
            call load_ero_launch_data

         else
            
            call errmsg('ERROR: ERO OPTIONS:',
     >         'NEUT2D_OPT=2 but ERO particle launch option is OFF'//
     >         ': Turning off NEUT2D')
            neut2d_opt = 0

         endif
      endif

c
c     ERO interface initialization - default particle not recorded
c      
      if (ero_part_output_opt.eq.0) then 
         ero_record_data = .false.
      elseif (ero_part_output_opt.eq.1) then 
         ero_record_data = .true.
         call init_ero_part_output
      endif

c
c     Initialize recombination, reflection, and MTC variables
c
      recstruk    = 0.0
      recloss     = 0.0
      mtcrecstruk = 0.0
      recwalln    = 0.0
      mtcrecwalln = 0.0
      recMAIN     = 0.0
      recEXIT     = 0.0
      recATIZ     = 0.0
      recNEUT     = 0.0
      recCENT     = 0.0
      recTMAX     = 0.0
      recFAIL     = 0.0
c
c     Ion Reflection
c
      refflag     = 0
      refstruk    = 0.0
      refloss     = 0.0
      mtcrefstruk = 0.0
      refwalln    = 0.0
      mtcrefwalln = 0.0
      refMAIN     = 0.0
      refEXIT     = 0.0
      refATIZ     = 0.0
      refNEUT     = 0.0
      refCENT     = 0.0
      refTMAX     = 0.0
      refFAIL     = 0.0
c
      call pr_trace('DIV','1:')
c
      call rzero (recinf,14*7)
      call rzero (rectotcnt,12)
      call rzero (mtcinf,7*3)
      call rzero (mtctotcnt,12*3)
c
c     The following are used to monitor leakage across the separatrix
c     into the "CORE" plasma region.
c
      hasleakedcore = .false.
      nleakcore = 0
      totleakcore = 0.0
c
      procterm = .false.
c
      lpinz0 = .false.
c
c     Initialise target mach number array to 1.0's
c
      call rinit(cmachno,maxnrs*2,1.0)
C
C     CALCULATE WSSF FACTOR USED IN THE CALCULATION OF ABSFAC
C
      IF (NIMPS2.GT.0) THEN
        IF (CNEUTH.EQ.2) THEN
          IF (CNEUTB.EQ.2) THEN
c
c       Change to 1.0 - even though this is not accurate - it allows
c       display of plots which are needed even in these cases.
c
            WSSF = 1.0
          ELSE
            WSSF = FLOAT(NIMPS2) / FLOAT(NIMPS) + 1.0
          ENDIF
        ELSEIF (CNEUTB.EQ.2) THEN
          WSSF = FLOAT(NIMPS)/FLOAT(NIMPS2) + 1.0
        ELSE
          WSSF = 1.0
        ENDIF
      ELSEIF (CNEUTB.EQ.2) THEN
c
c       Change to 1.0 - even though this is not accurate - it allows
c       display of plots which are needed even in these cases.
c
        WSSF = 1.0
      ELSE
        WSSF = 1.0
      ENDIF
C
C---- SET UP FACTORS IN COMMON COMTAU
C
      IF (DEBUGL) CALL TEST

      TAUTIM = ZA02AS (1)
c slmod begin
      LOAD_I = -1
      NYMFS_GLOBAL = NYMFS  ! lame, but don't want to pass NYMFS is local and I don't want to pass it around -SL, 21/11/2011
c slmod end

      call pr_trace('DIV','BEFORE TAU')


      IF (ITER.EQ.1) CALL TAUIN1 (title,equil,NIZS,VFLUID)
      TAUTIM = ZA02AS (1) - TAUTIM
      WRITE(6,*) 'TIME USED IN SETUP: TAU SUBROUTINE :',TAUTIM,' S'

      call pr_trace('DIV','AFTER TAU')

C
C---- ZERO ARRAYS ETC
C


      CALL DZERO (DDLIMS, MAXNKS*MAXNRS*(MAXIZS+2))
      CALL DZERO (DDTS,   MAXNKS*MAXNRS*(MAXIZS+2))
      call dzero (ddvoid, 3)
c
c     Zero arrays for recording initial locations of leaking
c     impurities.
c
      call rzero (ncore,  maxnks*maxnrs)
      call rzero (nedge,  maxnks*maxnrs)
      call rzero (ntrap,  maxnks*maxnrs)
      call rzero (ndivert,maxnks*maxnrs)
      call rzero (nmsol,  maxnks*maxnrs)
c
c     Zero arrays for recording average ion forces and velocity
c
c psmod
c
      CALL RZERO (Fcell, MAXNKS*MAXNRS*MAXIZS)
      CALL RZERO (Ffi, MAXNKS*MAXNRS*MAXIZS)
      CALL RZERO (Fthi, MAXNKS*MAXNRS*MAXIZS)
      CALL RZERO (Fvbg, MAXNKS*MAXNRS*MAXIZS)
      CALL RZERO (DIFF, MAXNKS*MAXNRS*MAXIZS)
      CALL RZERO (VELavg, MAXNKS*MAXNRS*MAXIZS)


c
c psmod
c
c     (RIV)
c
      debugv = .false.
      if (cstepv.ne.0.0) debugv = .true.
      if (debugv) then
         call rzero (sdvs,   maxnks*maxnrs*(maxizs+2))
         call rzero (sdvs2,   maxnks*maxnrs*(maxizs+2))
         call rzero (sdvs3,   maxnks*maxnrs*maxizs*2)
         call rzero (sdvb,maxnks*maxnrs)
         call rzero (velspace, (2*nvel+2)*maxvizs*maxvnks)
         call rzero (velweight,(2*nvel+2)*maxvizs*maxvnks)
c
c        Set velplate equal to the sound speed on ring 8 at the inner pl
c        and use this to scale the rest of the distributional analysis.
c
c        Adjust it by multiplying by qtim ... this will be accounted
c        for at the end.
c
c
         velplate=9.79E3* SQRT( ktibs(1,irsep)/CRMI) * qtim
c
         do ik = 1,nks(injir)
            velcell(ik) = 9.79E3* SQRT( ktibs(ik,injir)/CRMI) * qtim
         end do


c
c        Calculate average background plasma velocity at all points
c        - multiply by qtim to save testing at each iteration
c
         do ir = 1,nrs
            do ik = 1, nks(ir)
               sdvb(ik,ir)=9.79E3* SQRT(2.0*ktibs(ik,ir)/CRMB) * qtim
            end do
         end do
c
      endif


c
      CALL RZERO (LIMS,   MAXNKS*MAXNRS*(MAXIZS+2)*MAXNTS)
      CALL RZERO (ELIMS,  MAXNKS*3*(MAXIZS+2))
      CALL RZERO (DEPS,   MAXNDS*MAXIZS)
      CALL RZERO (NEROS,  MAXNDS*5)
      call rzero (promptdeps,maxnds*6)
      CALL RZERO (WALLS,  MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RZERO (wallsn, maxpts+1)
      CALL RZERO (wallse, maxpts+1)
      CALL RZERO (wallse_i, maxpts+1)
      CALL RZERO (wallsi, maxpts+1)
      CALL RZERO (wallsiz, (maxpts+1) * MAXIZS)
      call rzero (wallseiz,(maxpts+1) * MAXIZS)
      CALL RZERO (TNTOTS, (MAXIZS+2)*4)
      CALL RZERO (TIZS,   MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RZERO (ZEFFS,  MAXNKS*MAXNRS*3)
      CALL RZERO (SDTZS,  MAXIZS+2)
      CALL RZERO (SDTZS2, MAXIZS+2)
      CALL DZERO (DOUTS,  MAXIZS*9)
      CALL DZERO (coreOUTS,  MAXIZS*9)
      CALL DZERO (DPARAS, MAXIZS*6)
      CALL RZERO (RIONS,  MAXIZS)
      CALL RZERO (WALKS,  MAXNWS*2)
      call rzero (cleakn, maxpts*(maxizs+1))
c
      call dzero (dvmaxv,6*(maxizs+1))
      call dzero (dvminv,6*(maxizs+1))
c
c      call dinit (dvmaxv,6*(maxizs+1),-1.0d20)
c      call dinit (dvminv,6*(maxizs+1),1.0d20)
c
c     Launchdat contains ancillary data for both book-keeping and
c     later specifying how to launch each initial and sputtered
c     particle. This eventually will allow for every particle to
c     be launched with specific and differing conditions if such
c     is desired while still making use of the existing frame work.
c     Initially LAUNCHDAT is all zeroes - which will be used for
c     the DEFAULT conditions.
c
c     David Elder.                        July 12, 1995
c
c     - at this time only launchdat(imp,2) is used - if this is 0
c       it specifies a regular target sputter/relaunch. If it is 1
c       it specifies a target relaunch originating from the far
c       periphery.
c
c     - launchdat(imp,3) now specifies whether the launched neutral
c       underwent any wall reflections - and if so - how many. It
c       will be 0.0 for those with no wall reflections
c
c     - launchdat(imp,1) is now used to specify a specific launch option
c       to be used in LAUNCH for each particle - this will allow for a
c       mixture of target and wall recycled particles to be launched in
c       the same group.
c
c     - launchdat(imp,4) - specifies the starting charge state for the ion
c                          on return from neut - this allows for ions and
c                          neutrals to be injected from a common source in
c                          neut but leave the ion tracking to div 
c
c     - launchdat(imp,5) - ERO related - did particle already enter ERO sample volume
c
c

      call rzero (launchdat, maximp*4)
c
c     See comment in LAUNCH subroutine for information on
c     contents of ionizdat array and wtsource array.
c
      call rzero (ionizdat , 2*2*2*2*5)
      call rzero (wtsource,maxpts*maxnrs*4*6)
      call rzero (wtdep,maxpts*(maxpts+1)*3)
c
      IW = 1
c
c     General variable initialization
C
      DO 100  IZ = 1, NIZS
         CICIZS(IZ) = 0.0
         CIFIZS(IZ) = HI
         CILIZS(IZ) = 0.0
         CISIZS(IZ) = 0.0
c
         cieizs(iz) = 0.0
         citizs(iz) = 0.0
c
         CICABS(IZ) = 0.0
         CIFABS(IZ) = HI
         CILABS(IZ) = 0.0
         CISABS(IZ) = 0.0
         CRTABS(IZ) = 0.0
         CICUTS(IZ) = 0.0
         CRTRCS(IZ) = 0.0
         CRVABS(IZ) = 0.0
         CRAVAV(IZ) = 0.0
         CTBS  (IZ) = 0.0
         CICRCS(IZ) = 0.0
         CIFRCS(IZ) = HI
         CILRCS(IZ) = 0.0
         CISRCS(IZ) = 0.0
         CXXX  (IZ) = 0.0
         CSSS  (IZ) = 0.0
         CNNN  (IZ) = 0.0
         CNNNX (IZ) = 0.0
         CNNNS (IZ) = 0.0
         CNNNT (IZ) = 0.0
         CNNNKT(IZ) = 0.0
         CNNNK (IZ) = 0.0
         cnorgs(iz) = 0.0
         cnorgr(iz) = 0.0
         cnorgz(iz) = 0.0
         CLLL  (IZ) = 0.0
         CLLLX (IZ) = 0.0
         CLLLS (IZ) = 0.0
         CMMM  (IZ) = 0.0
         CMMMX (IZ) = 0.0
         CMMMS (IZ) = 0.0
         CICLOS(IZ) = 0.0
         CIFLOS(IZ) = HI
         CILLOS(IZ) = 0.0
         CISLOS(IZ) = 0.0
         CSSSS(IZ)  = 0.0
         FPTARG(IZ) = 0.0
         acttarg(iz) = 0.0
  100 CONTINUE
      FPTTOT = 0.0
      FPTART = 0.0
      FPENT  = 0.0
      FPEXIT = 0.0
      CALL RZERO (CTEXS, 10)
      CICCOL = 0.0
      CICRXA = 0.0
      CIFRXA = HI
      CISRXA = 0.0
      CITRXA = 0.0
      CICRIN = 0.0
      CIFRIN = HI
      CISRIN = 0.0
      CITRIN = 0.0
      CIKRIN = 0.0
      CKTRIN = 0.0
      CICRNO = 0.0
      CIKRNO = 0.0
      CIrRNO = 0.0
      CIzRNO = 0.0
      CIsRNO = 0.0
      CKTRNO = 0.0
      CICRNJ = 0.0
      CKKMAX =-HI
      CKKMIN = HI
      CVVXC  = 0.0
      CVVZM  = 0.0
      CVVRM  = 0.0
      CVVSM  = 0.0
      CVVKM  = 0.0
      CVVXE  = 0.0
      CVVXP  = 0.0
      CVVXS  = 0.0
c
      cvvrefm = 0.0
      cvvnrfm = 0.0
      cvvfpref = 0.0
      cvvfpnrf = 0.0
c
      CALL RZERO (DCROSS, 4)


C
C     INITIAL VALUES
C
      DIFFR = SQRT(2.0*CDPERPFP*QTIM)
c
c     Set up values of the drift velocity for each
c     flux tube.
c

      call setup_drftv


! jdemod - moved to before prdata since that routine needs the files
!          properly opened in order to print error messages
!     Initialize all data which needs to be set only once in the case.
      If (global_hc_follow_option .ne. 0) Then

         Call global_HC_Begin

      End If
! ammod end.

      call pr_trace('DIV','AFTER INIT')


C
C-----------------------------------------------------------------------
C     PRINT SELECTED PARAMETERS SET UP ON FIRST ITERATION ONLY
C-----------------------------------------------------------------------
C
      IF (ITER.EQ.1) CALL PRDATA (NIZS,NIMPS,NIMPS2,nymfs)

      call pr_trace('DIV','AFTER PRDATA')

c
c     Assign NIMPS and NIMPS2 to local copies in div5 common block
c
      nimps_local = nimps
      nimps2_local = nimps2


c
c     Calculate transport coefficients from OSM
c
      if (cpinopt.eq.1.or.cpinopt.eq.4) then
         call oskin
      endif
c
      if (cprint.eq.1.or.cprint.eq.9) then
         call probescan
      endif
c

c
c     Exit at this point if testing SOL options ONLY - this will get
c     a partial print out of case options and summary of SOL results.
c
      if (ctestsol.ne.0.0) then
c
c        Print out POWER summary even when impurities are not run.
c
         call pr_power_summary
c
         write(6,*) 'Testing SOL options only - DIVIMP exiting'
         write(6,*) 'No particles launched'
         call prc ('DIVIMP EXIT - SOL Option test ONLY')
         call prc ('            - No particles launched')
         return
      endif
c
c     Add header to define start of simulation information
c
      if (iter.eq.1) then
         call prb
         call prchtml('--- IMPURITY SIMULATION'//
     >        ' DIAGNOSTIC INFORMATION ---','pr_runtime','0','B')
         call prb
      endif
c
C
C-----------------------------------------------------------------------
C
c     If an impurity injection launch based on impurity
c     ionization data has been selected then it is necessary to generate
c     the  cumulative probability function and cross-index listing
c     to be used for this launch.
C
C-----------------------------------------------------------------------
C
c      write(6,*) 'ciopte:',ciopte
c
c IPP/01 geier added ciote.eq.8
      if (cneuta.eq.1.and.
     >    (ciopte.eq.4.or.ciopte.eq.7.or.ciopte.eq.8)) then
c
c        Perform additional processing to set up injection options 4 and 7
c
         if (ciopte.eq.4) then
c
c           Load PIN/NIMBUS neutral impurity density
c           This must be normalised to the number of
c           impurities ionized to be consistent with
c           DIVIMP's storage of impurity ion densities
c           Note that absfac must be calculated after
c           the ion following in order to properly
c           account for self-sputtering.
c
c           Also note that, for following self-sputtered
c           neutrals, cneuta can be set to 0 so a logical
c           saying that PIN neutral impurity densities
c           were used is required when deciding how
c           to calculate absfac.
c
            lpinz0 = .true.
            write(6,*) 'Normalise pinz0..., zioniz = ',zioniz
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  ddlims(ik,ir,0) = pinz0(ik,ir)/zioniz
               enddo
            enddo
c
c        Injection option 7 - other setup code - fill PIN arrays with
c        corresponding data loaded from the neutral code.
c
c IPP/01 geier - add ciopte=8
         elseif (ciopte.eq.7.or.ciopte.eq.8) then
c
c           Set up PIN arrays with appropriate values
c
            zioniz = 0.0
            zioniz_tor = 0.0
c
            do ir = 1,nrs
               do ik = 1,nks(ir)
c
c                 Convert e2diz0 to a density
c
                  if (karea2(ik,ir).gt.0.0) then
                     e2diz0(ik,ir) = e2diz0(ik,ir) / karea2(ik,ir)
                  else
                     e2diz0(ik,ir) = 0.0
                  endif
c
c                 Assign PIN arrays
c
                  pinionz(ik,ir) = e2diz0(ik,ir)
                  pinz0(ik,ir) = e2dz0(ik,ir)
                  pinenz(ik,ir) = ktibs(ik,ir)
c
c                 Calculate a value of zioniz for the ion source rate.
c
                  zioniz = zioniz + e2diz0(ik,ir) * karea2(ik,ir)
                  zioniz_tor = zioniz_tor + e2diz0(ik,ir)*karea2(ik,ir)
     >                                       * rs(ik,ir)
c
               end do
            end do
c
            lpinz0 = .true.
            write(6,*) 'Normalise e2dz0..., zioniz = ',zioniz,zioniz_tor
c
c           Copy neutral density to ddlims as for pinz0 for injection
c           option 4.
c
            do ir = 1,nrs
               do ik = 1,nks(ir)
                  ddlims(ik,ir,0) = e2dz0(ik,ir)/zioniz
               enddo
            enddo
c
         endif
c
C-----------------------------------------------------------------------
c        Calculate injection probabilities
C-----------------------------------------------------------------------
c

         injnum = 0
         do 2500 ir = 1,nrs
         do 2500 ik = 1,nks(ir)
               iprob = pinionz(ik,ir) * karea2(ik,ir)
c
               if (cprint.eq.7.or.cprint.eq.9) then
                  write(6,'(a,2i6,3(1x,g15.8))') 'injp:',
     >                  ik,ir,pinionz(ik,ir),karea2(ik,ir),iprob
               endif
c
               if (iprob.gt.0.0) then
                  injnum = injnum +1
                  if (injnum.eq.1) then
                     injprob(injnum) = iprob
                  else
                     injprob(injnum) = injprob(injnum-1) + iprob
                  endif
                  injrind(injnum) = ir
                  injkind(injnum) = ik
               endif
 2500    continue
         do 2600 in = 1,injnum
            injprob(in) = injprob(in)/injprob(injnum)
            write(6,*) 'inj:',injkind(in),injrind(in),injprob(in)
 2600    continue
c
      endif


c
c     The code has been modified to run PIN after any SOL calculations
c     and so this extra call is unnecessary.
c
c
c     If Injection option is 4 - i.e. from PIN result - then be sure to
c     run PIN on the finalized background plasma - one does not
c     want to just use the PIN result from a past iteration of the
c     background plasma solver.
c
c      if (ciopte.eq.4.and.cneuta.eq.1) then
c
c       Set up for PIN call
c
c        CALL WRTPIN(title,equil,17)
c
c        write(0,*) 'Calling PIN for C injection data'
c        CALL INVOKEPIN(ACTPIN,pintim)
c        write(0,*) 'Return from PIN after ',pintim,' (seconds)'
c
c        totpintim = totpintim + pintim
c
c       Load PIN results
c
c        CALL READPIN
c
c      endif
c
C
C-----------------------------------------------------------------------
C    LAUNCH PRIMARY NEUTRALS EN MASSE
C-----------------------------------------------------------------------
C
C---- TNEUT : TOTAL NUMBER OF NEUTRALS LAUNCHED
C---- TATIZ : TOTAL NO OF IONS CREATED
C---- TWALL : TOTAL NO OF IONS REACHING WALL
C---- TWALLN: TOTAL NO OF NEUTRALS REACHING WALL
C---- TDEP  : TOTAL NO OF IONS DEPOSITED ON TARGET
C---- TTMAX : TOTAL NO OF NEUTRALS EXISTING AT TMAX
C---- TCENT : TOTAL NO OF NEUTRALS REACHING CENTRAL MIRROR
C---- TBYOND: TOTAL NO OF IONS IONISED BEYOND LIMIT
C---- TBELOW: TOTAL NO OF IONS RECOMBINING TO NEUTRALS
C---- TCUT  : TOTAL NO OF IONS EXISTING AT TMAX
C---- TSTRUK: TOTAL NO OF NEUTRALS STRIKING TARGET
C---- TFAIL : TOTAL NO OF FAILED NEUTRAL LAUNCHES
C---- TMAIN : TOTAL NEUTRALS REACHING MAIN PLASMA
C
      TNEUT  = 0.0
      TATIZ  = 0.0
      TWALL  = 0.0
      TWALLN = 0.0
      walltotn = 0.0
      TDEP   = 0.0
      TTMAX  = 0.0
      TCENT  = 0.0
      TBYOND = 0.0
      TBELOW = 0.0
      TCUT   = 0.0
      TSTRUK = 0.0
      TFAIL  = 0.0
      TMAIN  = 0.0
C
C---- RNEUT1: NUMBER OF PRIMARY NEUTRALS LAUNCHED
C---- RNEUT : NUMBER OF NEUTRALS LAUNCHED IN CURRENT BATCH
C---- RATIZ : NUMBER OF IONS CREATED IN CURRENT BATCH
C---- RSTRUK: NUMBER OF NEUTRALS STRIKING TARGET
C---- STATUS: 1,2,3..9 FOR PRIMARY,2ND,3RD LAUNCHES ETC, 10 FOR TOTAL.
C---- RMAIN : NUMBER REACHING MAIN PLASMA
C
C
C     PLACE CALL TO INITIALIZE THE SPUTTERING YIELD DATA SO THAT
C     IT IS AVAILABLE FOR ION AS WELL AS NEUT LAUNCH CASES.
C
C
C     LOAD YIELD COMMON BLOCK WITH APPROPRIATE DATA
C
      IF (CSPUTOPT.EQ.1) THEN
        CALL SYIELD (MATTAR,MATP,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
      ELSE IF (CSPUTOPT.EQ.2) THEN
        CALL SYLD93 (MATTAR,MATP,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
      ELSE IF (CSPUTOPT.EQ.3.or.csputopt.eq.4.or.csputopt.eq.5.or.
     >         csputopt.eq.6)THEN
        CALL SYLD96 (MATTAR,MATP,CNEUTD,
     >               CBOMBF,CBOMBZ,CION,CIZB,CRMB,CEBD)
        call init_eckstein_2007(mattar,matp)
      ENDIF
c
c     jdemod - print out sputtering yield data for current case
c
      if (cprint.eq.9) then 
         call print_sputtering_yields(mattar,matp,crmb)
      endif


C
C     SET YIELD MULTIPLICATION VALUES SO THAT THEY ARE AVAILABLE
C     FOR BOTH NEUTRAL AND ION INJECTION CASES.
C


c
c     Set defaults for targets
c
      DO ID = 1, NDS
c
c       Target Yield multipliers
c
        if (nymfs.gt.0.and.cymfs(1,1).eq.0) then
          KMFPS(ID) = CYMFS(1,3)
          KMFSS(ID) = CYMFS(1,4)
          KMFCS(ID) = CYMFS(1,5)
        else
          KMFPS(ID) = 1.0
          KMFSS(ID) = 1.0
          KMFCS(ID) = 1.0
        endif
c
      end do
c
c     Set defaults for walls
c
      DO ID = 1, wallpts
c
c       Target and Wall yield multipliers
c
        if (nymfs.gt.0.and.cymfs(1,1).eq.0) then
           if (wallpt(id,16).eq.1.or.wallpt(id,16).eq.4) then
              KMFPWS(ID) = CYMFS(1,3)
              KMFCWS(ID) = CYMFS(1,5)
           else
              KMFPWS(ID) = CYMFS(1,6)
              KMFCWS(ID) = CYMFS(1,7)
           endif
c
c          default wall reflection value
c
           wallpt(id,25) = cymfs(1,8)
c
        else
           KMFPWS(ID) = 1.0
           KMFCWS(ID) = 1.0
c
c          default wall reflection value
c
           wallpt(id,25) = 0.0
c
        endif
c
      end do
c
c     Map YMF's and reflection values to individual segments
c
      do in=1,nymfs
c
c        Loop through specified segments - wall indices only
c
         do id = cymfs(in,1),cymfs(in,2)
c
            if (id.ge.1.and.id.le.wallpts) then
c
c              Set reflection value
c
               wallpt(id,25) = cymfs(in,8)
c
               if (wallpt(id,16).eq.1.or.wallpt(id,16).eq.4) then
                  KMFPWS(ID) = CYMFS(in,3)
                  KMFCWS(ID) = CYMFS(in,5)
c
c                 Assign corresponding regular target values.
c
                  if (wallpt(id,18).ge.1.and.wallpt(id,18).le.nds) then
                     KMFPS(int(wallpt(ID,18))) = CYMFS(in,3)
                     KMFSS(int(wallpt(id,18))) = CYMFS(in,4)
c
c                    Check for target ion reflection and set flag
c
                     if (kmfss(int(wallpt(id,18))).le.-99.0) then
c
                        refflag = 1
c
                     endif
c
                     KMFCS(INT(wallpt(id,18))) = CYMFS(in,5)
                  endif
c
               else
                  KMFPWS(ID) = CYMFS(in,6)
                  KMFCWS(ID) = CYMFS(in,7)
               endif
            endif
         enddo
c
      enddo
c
c     ERODIV - if the option is set to turn off sources inside the ERO volume then loop through 
c              the wall/target locations and zero the yields for sections "substantially"
c              inside the ERO sample volume
c
c
      if (ero_remove_src_opt.gt.0) then 
         call ero_remove_src(wallpts,wallpt,kmfps,kmfcs,kmfss,
     >             kmfpws,kmfcws,maxpts,nwall_data,maxnds)

      endif



c     Print out some wall chatracteristics
c
c      write (6,*) 'Some wall values:'
c      do in = 1,wallpts
c         write (6,'(i5,3(1x,g12.5))') in,wallpt(in,18),
c     >         kmfss(wallpt(in,18)),wallpt(in,25)
c      end do
c
c      do in = 1,nds
c         write (6,'(a,2i5,4(1x,g12.5))') 'Wall index:',
c     >          in,wallindex(in),
c     >          kteds(in),ktids(in),knds(in),kvds(in)
c      enddo
c

      call pr_trace('DIV','AFTER SPUTTERING YIELDS')


! ammod begin.
C---- TOGGLE KINDS,KOUTDS BETWEEN 0 AND INFINITY
C---- (AND WHEN CTARGOPT = 0 OR 4 KBACDS, KFORDS)
! Note, this was moved above the call to neut to
! facilitate ion transport in the HC routines.
c
      DO 111 IR = 1, NRS
        IF (CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4) THEN
          IF (IR.GE.IRSEP) KBACDS(1,IR) = HI
          IF (IR.GE.IRSEP) KFORDS(NKS(IR),IR) = HI
        ENDIF
        DO 111 IK = 1, NKS(IR)
          IF (IR.EQ.1.OR.IR.EQ.IRTRAP.or.ir.eq.irtrap2)
     >       KINDS(IK,IR) = HI
          IF (IR.EQ.IRWALL.or.ir.eq.irwall2) KOUTDS(IK,IR) = HI
  111 CONTINUE


c
c     Launch neutrals
c
      call pr_trace('DIV','LAUNCH NEUTRALS')

      IF (CNEUTA.EQ.0) THEN
        STATUS = 1
        WRITE (6,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'
        WRITE (7,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'

        CALL NEUT (NATIZ,MATP,
     >             MATTAR,NIMPS,NIMPS2,FTOT,FYTOT,neut2d_fytot,
     >             RSTRUK,MTCSTRUK,
     >             RMAIN,REXIT,
     >             RATIZ,RNEUT,RWALLN,MTCWALLN,RCENT,RTMAX,SEED,NRAND,
     >             NEUTIM,RFAIL,NYMFS,STATUS)
        IF (NATIZ.EQ.0) GOTO 806
        WRITE(0,*) 'NEUT:NATIZ: ',NATIZ

      ELSE
        STATUS = 1
        NATIZ  = NIMPS
        RATIZ  = REAL (NIMPS)
        RNEUT  = 0.0
        RWALLN = 0.0
        MTCWALLN = 0.0
        RCENT  = 0.0
        RTMAX  = 0.0
        RSTRUK = 0.0
        MTCSTRUK = 0.0
        RFAIL  = 0.0
        RMAIN  = 0.0
        REXIT  = 0.0
c
      ENDIF

c
      RNEUT1 = RNEUT
      CLLL(-1) = REXIT
      CMMM(-1) = 0.0
      CNNN(-1) = RMAIN
C
C-----------------------------------------------------------------------
C     SELF-SPUTTERING CONTINUATION POINT ... PREPARE FOR MAIN LOOP
C-----------------------------------------------------------------------
C

      call pr_trace('DIV','SELF-SPUTTER LOOP START')

  200 CONTINUE
      IF (NIZS.GT.0 .AND. ITER.EQ.1) CALL TAUIN2 (NIZS)
      NPROD  = 0
      DO 201 M = 1, 2
        AVXPOS(M) = 0.0
        AVYPOS(M) = 0.0
        AVATIZ(M) = LO
  201 CONTINUE
      AVKPOS = 0.0
      AVSPOS = 0.0
      AVSMAX = 0.0
      AVVPOS = 0.0
      RWALL  = 0.0
      RDEP   = 0.0
      rfptarg = 0.0
      YLDTOT = 0.0
      YTHTOT = 0.0
      YLDMAX = 0.0
      RDIFFT = 0.0
      IF (STATUS.LE.50) THEN
         WRITE (6,9012) '***  FOLLOWING ',WHAT(STATUS),'   IONS    ***'
         WRITE (7,9012) '***  FOLLOWING ',WHAT(STATUS),'   IONS    ***'
      ELSE
         WRITE (6,9013) '***  FOLLOWING ',STATUS,' GENERATION IONS ***'
         WRITE (7,9013) '***  FOLLOWING ',STATUS,' GENERATION IONS ***'
      ENDIF
      DEBUGL = .FALSE.
      IF (CSTEPL.GT.0.0) THEN
        WRITE (6,9004) NINT(CSTEPL),QTIM
        DEBUGL = .TRUE.
      ENDIF
      CISTOT = 0.0
      CISMAX = 0.0
      STATIM = ZA02AS (1)
      PORM   = -1.0
      KK     = 1000 * ISECT
      KKLIM  = KK - 10
      DQTIM  = DBLE(QTIM)
      IMPLIM = 4
c
      num_entered_core = 0.0

      call pr_trace('DIV','MAIN LOOP START')

C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     M A I N   L O O P   S T A R T
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c sltmp
      IF (grdnmod.NE.0) iw = MAX(1,stopopt)

      tdep_save_n = 0

      DO 800  IMP = 1, NATIZ
c
c     jdemod - Commented out this debug line - only useful for reporting 
c              that essentially every 10% of ions are complete. 
c            - not sure why it would have a dependence on grdnmod either
c
         if (mod(imp,natiz/10).eq.0) then 
            perc = int((imp*10)/(natiz/10))
            write(0,'(a,i3,a,i8)') 
     >         'Following Ions: ',perc,' % complete. Particle # =',imp
         endif
c
c slmod begin
c        IF (.TRUE..AND.grdnmod.NE.0.AND.MOD(imp,natiz/10).EQ.0)
c        IF (sloutput) THEN 
c          IF ((natiz.GT.10.AND.
c     .         grdnmod.NE.0.AND.MOD(imp,natiz/10).EQ.0).OR.
c    .        (natiz.LE.10)) 
c     .      WRITE(0,*) 'debug imp:',imp,natiz
c        ENDIF
c slmod end
c
c       Particle initialization
c
        hasleaked = .false.
        hasleakedcore = .false.
        cleakp = 1
c
c       Recombination counter
c
        reccnt = 0
c
c       Initialize the IFATE value to zero
c
        ifate = 0
c
c        if (imp.eq.87) then
c           cstepl = 1.0
c           debugl = .true.
c        elseif (imp.eq.88) then
c           write (6,*) 'IMP 88: STOP'
c           stop
c        else
c           cstepl = 0.0
c           debugl = .false.
c        endif
c
C     WRITE(6,*) 'IMP:',IMP,' OF ', NATIZ
C
C------ SET LOCAL COPIES OF CHARACTERISTIC TIMES DATA, THESE MAY BE
C------ CHANGED WHEN TEMI CHANGES IF OPTIONS OTHER THAN 0 ARE USED.
C
        DO 777 IZ = 1, NIZS
          XXX(IZ) = ZMAX
          SSS(IZ) = 0.0
          SSSS(IZ)= 0.0
          DO 776 IR = 1, NRS
            DO 775 IK = 1, NKS(IR)
              LTOLDS(IK,IR,IZ) = CTEMAV
              LFPS  (IK,IR,IZ) = MAX (LO,KFPS(IK,IR,IZ))
              LFTS  (IK,IR,IZ) = MAX (LO,KFTS(IK,IR,IZ))
              LFSS  (IK,IR,IZ) = MAX (LO,KFSS(IK,IR,IZ))
              LLLFPS(IK,IR,IZ) = KKKFPS(IK,IR,IZ)
  775       CONTINUE
  776     CONTINUE
  777   CONTINUE
C
C------ SELECT IONS LAUNCH POINT FROM NEUT
C
c
c        SET IK,IR TO 0,0 - so that if they are changed in
c        the following that will fall through to the gridpos
c        routine. (As a double check.)
c
        ik = 0
        ir = 0
        idstart = 0
        idtype = -1
c
        IF (CNEUTA.EQ.0) THEN
          R     = XATIZS(IMP)
          Z     = YATIZS(IMP)
          VEL   = VINS(IMP) * QTIM

          SPUTY = SPUTYS(IMP)
          TEMI  = TEMTIZS(IMP)
          idstart = idatizs(imp,1)
          idtype  = idatizs(imp,2)
c
c         Assign starting wall index of particle
c
          if (idtype.eq.1.or.idtype.eq.2.or.idtype.eq.3) then
             iwstart = wallindex(idstart)
          else
             iwstart = idstart
          endif
c
          cist  = cistizs(imp)
c

C
C------ SELECT INJECTION POINTS FROM INPUT OPTIONS
C
        ELSEIF (CNEUTA.EQ.1) THEN
c
          IF (CIOPTE.EQ.1) THEN
            R     = CXSC
            Z     = CYSC
            PORM  = -1.0 * PORM
c
c            VEL   = 1.56E4 * SQRT (CTEM1/CRMI) * PORM * QTIM
c
            VEL   = 9.79e3 * SQRT (CTEM1/CRMI) * PORM * QTIM
c
            SPUTY = 1.0
c slmod begin
          ELSEIF (CIOPTE.EQ.9) THEN
c...        Impurity injection at a random point along a specified
c           line (SL 25.6.2003):

            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            R = CXSCA + RAN * (CXSCB - CXSCA)

c            NRAND = NRAND + 1
c            CALL SURAND2 (SEED, 1, RAN)
            Z = CYSCA + RAN * (CYSCB - CYSCA)


c            R     =  1.300 + RAN * (1.710 - 1.300)
c            Z     = -1.1
c            R     =  1.350 + RAN * (1.684 - 1.350)
c            Z     = -1.2
c            R     =  1.391 + RAN * (1.673 - 1.391)
c            Z     = -1.3

            PORM  = -1.0 * PORM
            VEL   = 9.79E+03 * SQRT (CTEM1 / CRMI) * PORM * QTIM
            SPUTY = 1.0
c slmod end
c jdemod begin
          ELSEIF (CIOPTE.EQ.10) THEN
c...        Impurity injection at a random location within a rectangle specified
c           by the corner pionts:
c           R1,Z1 = (CXSCA,CYSCA)
c           R2,Z2 = (CXSCB,CYSCB)
c
c           line (SL 25.6.2003):

            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            R = CXSCA + RAN * (CXSCB - CXSCA)

            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            Z = CYSCA + RAN * (CYSCB - CYSCA)


c            R     =  1.300 + RAN * (1.710 - 1.300)
c            Z     = -1.1
c            R     =  1.350 + RAN * (1.684 - 1.350)
c            Z     = -1.2
c            R     =  1.391 + RAN * (1.673 - 1.391)
c            Z     = -1.3

            PORM  = -1.0 * PORM
            VEL   = 9.79E+03 * SQRT (CTEM1 / CRMI) * PORM * QTIM
            SPUTY = 1.0
c jdemod end
          ELSEIF (CIOPTE.EQ.2.or.ciopte.eq.3.or.ciopte.eq.5
     >            .or.ciopte.eq.6) THEN
C
C           FOR NOW SIMPLY FIND GRID POINT CLOSEST TO DESIRED
C           INJECTION POSITION. THIS IS ALL THAT IS CURRENTLY DONE.
C           THE PARTIAL DISPLACEMENTS FROM THE GRID POINTS ARE IGNORED.
C
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            STMP = RAN*(INJF2-INJF1)*KSMAXS(INJIR)
     >             + INJF1*KSMAXS(INJIR)
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            IF (RAN.GT.0.5.and.(ciopte.eq.2.or.ciopte.eq.5)) THEN
C
C             OUTER PLATE - OTHERWISE INNER - for option 2
C
              STMP = KSMAXS(INJIR) - STMP
            ENDIF
C
C           FIND R,Z COORDINATES OF NEAREST GRID POINT
C
C-------- FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR INJI
C
            IK = 1
  778       IF (IK.LT.NKS(INJIR).AND.STMP.GT.KSS(IK,INJIR)) THEN
              IK = IK + 1
              GOTO 778
            ENDIF
  779       IF (IK.GT.1.AND.STMP.LE.KSS(IK-1,INJIR)) THEN
              IK = IK - 1
              GOTO 779
            ENDIF
            IF (IK.GT.1.AND.
     >         (STMP-KSS(IK-1,INJIR).LT.KSS(IK,INJIR)-STMP))
     >          IK = IK - 1
C
            ir = injir
c
            R = RS(IK,ir)
            Z = ZS(IK,ir)
c
c            write (6,*) 'stmp:',stmp,ik,injir,kss(ik,injir),
c     >                  r,z,ksmaxs(injir)
C
C           SET UP VELOCITY - BASED ON ION TEMPERATURE
C
            PORM  = -1.0 * PORM
c
c            VEL   = 1.56E4 * SQRT (CTEM1/CRMI) * PORM * QTIM
c
            if (ciopte.eq.2.or.ciopte.eq.3) then
               VEL   = 9.79e3 * SQRT (CTEM1/CRMI) * PORM * QTIM
            elseif (ciopte.eq.5.or.ciopte.eq.6) then
 7701          nrand = nrand + 1
               call surand2(seed,1,ran1)
               if (ran1.eq.0.0) goto 7701
               nrand = nrand + 1
               call surand2(seed,1,ran2)
               if (ran1.eq.0.0) goto 7701
               rgauss = sqrt(-2.0* log(ran1))*cos(2.0*PI*ran2)
               VEL   = 9.79e3 * SQRT (CTEM1/CRMI) * QTIM
     >                 * rgauss
            endif
c
            SPUTY = 1.0
c
          elseif (ciopte.eq.4.or.ciopte.eq.7.or.ciopte.eq.8) then
c
c           Find the initial injection bin
c
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            in = ipos(ran,injprob,injnum)
            ik = injkind(in)
            ir = injrind(in)
c
            R     =   rs(ik,ir)
            Z     =   zs(ik,ir)
c
c           Set up velocity based on energy
c
            PORM  = -1.0 * PORM
c
c           Verify conversion of pinenz ...
c           This :
c            VEL   = 1.56E4 * SQRT(2.0*pinenz(ik,ir)/CRMI)
c     >              * 0.5  * PORM * QTIM
c           assumes pinenz is specified in eV.
c
c            VEL   = 1.56E4 * SQRT(2.0*pinenz(ik,ir)/CRMI)
c     >              * 0.5  * PORM * QTIM
c
            if (ciopte.eq.8) then
c
               sigma_vel=sqrt(9.648e7*pinenz(ik,ir)/crmi)
               vel= ndrand(sigma_vel,0.)*qtim
c
            else
c
               VEL   = 9.79E3 * SQRT(2.0*pinenz(ik,ir)/CRMI)
     >              * 0.5  * PORM * QTIM
c
            endif
c
c
            SPUTY = 1.0
c slmod begin - t-dep
          ELSEIF (CIOPTE.EQ.11) THEN
c...        Ion injection from both the current source and the source stored
c           from a previous run:
            LOAD_I = -1
            NRAND = NRAND + 1
            CALL SURAND2(SEED, 1, RAN)
            write(0,*) 'branch',ran,tdep_load_frac
            IF (RAN.GT.TDEP_LOAD_FRAC) THEN
              write(0,*) 'standard'
              R     = CXSC
              Z     = CYSC
              PORM  = -1.0 * PORM
              VEL   = 9.79E3 * SQRT(CTEM1 / CRMI) * PORM * QTIM
              SPUTY = 1.0
            ELSE
              NRAND = NRAND + 1
              CALL SURAND2 (SEED, 1, RAN)
              LOAD_I =MIN(MAX(1,INT(REAL(TDEP_LOAD_N)*RAN)),TDEP_LOAD_N)
c... left off: need to sort out setting the charge state, and also the strange initialization of 
c maxciz from cizsc :
c jdemod - maxciz is the maximum charge state reached by THIS particle ... so it starts at initial value and gets incremented later
c          DO 792 JZ = CIZSC, MAXCIZ
c why would this be done?
c a bug?
              R = TDEP_LOAD(LOAD_I)%R
              Z = TDEP_LOAD(LOAD_I)%Z
              PORM  = -1.0 * PORM
              VEL  = TDEP_LOAD(LOAD_I)%VEL
c should be adjusting sputy, abs
              SPUTY = TDEP_LOAD(LOAD_I)%WEIGHT

              write(0,*) '_load',load_i,r,z,vel,sputy

            ENDIF
c slmod end
          ENDIF
c
c         Set initial ion temperature.
c
c Geier IPP/01 added  .or.ciopte.eq.8
c slmod begin
          if (load_i.ne.-1) then
            temi = tdep_load(load_i)%temp
          elseif (ciopte.eq.4.or.ciopte.eq.7.or.ciopte.eq.8) then
c
c          if (ciopte.eq.4.or.ciopte.eq.7.or.ciopte.eq.8) then
c slmod end
             TEMI = pinenz(ik,ir)
          else
             TEMI = CTEM1
          endif
c
c
c         Set initial time step for injected ions
c

          cist = 1.0
c
        ENDIF
c
c       Temporarily record starting position of ion so that it
c       can be recorded later if necessary.
c
        rstart = r
        zstart = z
C
        IF (IW.LT.MAXNWS) THEN
          WALKS(IW,1) = 100.0 * RMAX
          WALKS(IW,2) = 100.0 * ZMAX
          IW = IW + 1
          IF (IW.LT.MAXNWS) THEN
            WALKS(IW,1) = R
            WALKS(IW,2) = Z
            WALKS(IW+1,1) = HI
            WALKS(IW+1,2) = HI
            IW = IW + 1
          ENDIF
        ENDIF


c
c     ERO record particle data flags
c     

      if (ero_part_output_opt.eq.0) then 
         ero_record_data = .false.
      elseif (ero_part_output_opt.eq.1) then 
         if (launchdat(imp,5).eq.1) then 
            ero_record_data = .false.   ! already recorded as a neutral
         else
            call ero_init_part_track(r,z,ero_record_data)
         endif
      endif


C
C------ SET ARRAY INDICES ETC
C
        TSTEPL = CSTEPL
c
c       Moved assignment of TEMI up to position appropriate to
c       ion injection vs. neutral launch
c
c        TEMI   = CTEM1
c
c slmod begin
        IF (LOAD_I.NE.-1) THEN
          IZ     = NINT(TDEP_LOAD(LOAD_I)%CHARGE)
          MAXCIZ = IZ         ! not sure about this...
c
c         jdemod - maxciz is the maximum charge state reached by this particle and it is incremented 
c                  as IZ changes ... so should start at IZ. 
c
c          MAXCIZ = CIZSC         ! not sure about this...
c
c       jdemod - add code for particles launched using ERO distribution 
c              - can have any initial charge state - neutrals have been handled in neut
c       
        elseif (ero_particle_launch_opt.eq.1.and.neut2d_opt.eq.2) then 
           iz = nint(launchdat(imp,4))
           maxciz = iz
        ELSE
          IZ     = CIZSC
          MAXCIZ = CIZSC
        ENDIF


        RIZ    = REAL(IZ)
        DSPUTY = DBLE(SPUTY)
c
c        IZ     = CIZSC
c        RIZ    = REAL(IZ)
c        MAXCIZ = CIZSC
c        DSPUTY = DBLE (SPUTY) 
c slmod end
c
c       Move setting of initial value of cist so that time spent as
c       neutrals can be included.
c
c        CIST   = 1.0
c
c        cisterr = .false.
c
c       Initialize time in ionization state.
c
        cistiz = 1.0
c
        IT     = 1
        OLDZ   = Z
C
C------ FIND NEAREST POINT IN (R,Z) CONTOUR MAP
C
        if (debugl)
     >   write (6,'(a,3i7,6(1x,f12.5))') 'ION:',
     >           imp,ik,ir,rstart,zstart,sstart

c
        call gridpos(ik,ir,r,z,.true.,griderr)
c
        if (griderr) then
c
c          A grid error at this point indicates a serious
c          mistake ...  since all ions must travel inside
c          the plasma. Also all the R,Z values being
c          returned by NEUT or gnerated by the
c          ion injection algorithm should be inside
c          the plasma ... This could possibly be invoked
c          if an R,Z injection position outside the grid
c          was specified.
c
           write (6,'(a,3i6,3(1x,g12.5))')
     >          'GRID ERROR: Particle not on grid',
     >          imp,ik,ir,cist,r,z
           write (0,'(a,3i6,3(1x,g12.5))')
     >          'GRID ERROR: Particle not on grid',
     >          imp,ik,ir,cist,r,z
           goto 800
c
c           STOP
c
        endif

c        write (0,'(a,4i7,5(1x,f12.5))') 'ION:',
c     >           imp,ik,ir,iz,rstart,zstart,sstart,sputy



c
c       SET Initial S and CROSS postion for particles.
c
C
c slmod begin - t-dep
        if (load_i.ne.-1) then
          cross = tdep_load(load_i)%cross
          k     = kks(ir)
          s     = tdep_load(load_i)%s
        elseif (init_pos_opt.eq.0) then
c
c
c        if (init_pos_opt.eq.0) then
c slmod end
c
           CROSS  = 0.0
           K      = KKS(IR)
           if (cneuta.eq.1.and.(ciopte.eq.2.or.ciopte.eq.3
     >         .or.ciopte.eq.5.or.ciopte.eq.6)) then
              S  = STMP
           else
              S  = KSS(IK,IR)
           endif
c
        elseif (init_pos_opt.eq.1) then
c
           CROSS  = 0.0
           K      = KKS(IR)
           if (cneuta.eq.1.and.(ciopte.eq.2.or.ciopte.eq.3
     >         .or.ciopte.eq.5.or.ciopte.eq.6)) then
              S  = STMP
           else
              call getscross_approx(r,z,s,cross,ik,ir)
           endif
        endif
c
c       Record approximate starting S-distance from nearest target.
c
        sstart = min(s,ksmaxs(ir)-s)
c
c       Krieger IPP/08 - added debug output line
c       write(0,'(a,i6,a,f7.1)') 'Progress indicator: ion=',imp,
c    >                           '  time=',za02as(1)-statim
        if ( status.le.10.and.
     >      ((natiz.lt.1000).or.
     >       (natiz.lt.10000.and.(imp/10)*10.0.eq.imp).or.
     >       ((imp/100)*100.0.eq.imp)) ) then

           write (6,'(a,i6,2i4,5(1x,g13.5))')
     >                   'ION-A:',imp,ik,ir,rstart,zstart,sstart,
     >                       cist
c     >                       ,ZA02AS (1) - STATIM
        endif
c
c
        M  = 2
        if (cgridopt.eq.2) then
           if ((ir.ge.irsep.and.ir.le.irwall2).or.
     >        ( ((ir.ge.irtrap2.and.ir.le.nrs2).or.
     >           (ir.ge.irtrap.and.ir.le.nrs).or.
     >           (ir.lt.irsep)).and.
     >            ik.le.(nks(ir)/2))) m=1
        elseIF (IK.gt.NKS(IR)/2) then
           M = 1
        endif

c
c       SET Theta value for non-orthogonal transport
c
        if (northopt.eq.1.or.northopt.eq.3) then
c
c          Calculate theta is in ion_parallel_transport.f
c
           call calculate_theta(ik,ir,s,theta)
c
c
c       nonorth
c
c        IF (S.GT.KSS(IK,IR)) THEN
c          IF (IK.LT.NKS(IR)) THEN
c            THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
c     >                              *(THETAG(IK+1,IR)-THETAG(IK,IR))
c          ELSE
c            THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
c     >                              *(THETAT(IDDS(IR,1))-THETAG(IK,IR))
c          ENDIF
c        ELSE
c          IF (IK.GT.1) THEN
c            THETA = THETAG(IK,IR) + (KSS(IK,IR)-S)/KBACDS(IK,IR)
c     >                              *(THETAG(IK,IR)-THETAG(IK-1,IR))
c          ELSE
c            THETA = THETAG(IK,IR) + (KSS(IK,IR)-S)/KBACDS(IK,IR)
c     >                              *(THETAG(IK,IR)-THETAT(IDDS(IR,2)))
c          ENDIF
c        ENDIF
c
c       nonorth
c
        endif
c
        SMAX   = KSMAXS(IR)
c
c       Set range for poloidal drift velocity
c
c
c       The drift range is now done in setup_drftvel on a ring by ring basis
c
c        if (cpdrft.eq.1.or.cpdrft.eq.2.or.cpdrft.eq.3) then
c           sdrft_start = cdrftv_start*smax
c           sdrft_end   = cdrftv_end*smax
c        endif
c
        AVXPOS(M) = AVXPOS(M) + R * SPUTY
        AVYPOS(M) = AVYPOS(M) + Z * SPUTY
        AVATIZ(M) = AVATIZ(M) + SPUTY
c        write(0,'(a,5i6,5(1x,g12.5))') 'AVATIZ:',imp,ik,ir,iz,m,sputy,
c     >              avatiz(m),ratiz
        AVKPOS = AVKPOS + K * SPUTY
        AVSPOS = AVSPOS + MIN (S, SMAX-S) * SPUTY
        AVSMAX = AVSMAX + SMAX * SPUTY
        AVVPOS = AVVPOS + ABS(VEL)/QTIM * SPUTY
        XXX(IZ)= Z
        SSS(IZ)= MIN (S, SMAX-S)
        KKK    = K
        CVVXC  = CVVXC  + Z * SPUTY
        XTRIPP = 0.0
        XTRIPS = 0.0
C
        IF (IMP.LE.4.and.(status.le.100.or.
     >      real(int(status/100)).eq.real(status)/100.0)) THEN
          C(1) = FACTOR (KES  (IK,IR),7)
          C(2) = FACTOR (KVHS (IK,IR),8)
          C(3) = FACTOR (LFPS (IK,IR,1),2)
          C(4) = FACTOR (LFSS (IK,IR,1),3)
          C(5) = FACTOR (LFTS (IK,IR,1),4)
          C(6) = FACTOR (KFEGS(IK,IR)*KALPHS(1),11)
          C(7) = FACTOR (KFIGS(IK,IR)*KBETAS(1),11)
          WRITE (6,9010) IMP,KTEBS(IK,IR),KNBS(IK,IR),C(1),C(2),
     >      KFIZS(IK,IR,1),KFRCS(IK,IR,1),KPCHS(IK,IR,1),
     >      KPRCS(IK,IR,1),
     >      KFCXS(IK,IR,1),KNHS(IK,IR),C(3),C(4),C(5),KKS(IR),
     >      KFIZS(IK,IR,2),KTIBS(IK,IR),LLLFPS(IK,IR,1),C(6),C(7)
        ENDIF
C
C------ CHECK IF ION HAS STARTED ABOVE MAX IONIZATION STATE
C------ IF SET TI=TB FOR STATE IZ APPLIES, BETTER DO IT
C
        IF ((IZ .GT. CION) .OR. (IZ .GT. NIZS))    THEN
          TBYOND = TBYOND + SPUTY
          IFATE = 4
          GOTO 790
        ENDIF
        IF (IZ.EQ.CIZSET) TEMI = MAX (TEMI,KTIBS(IK,IR))
C
C------ CALCULATE TIME AT WHICH DIFFUSION WILL BE FIRST APPLIED.
C
        RCONST = HI
        IF (CDIFOP.EQ.0) THEN
C
C       THE FOLLOWING LINE WAS COMMENTED OUT SINCE IT SEEMED MORE
C       REASONABLE TO ALLOW RCONST=0 FOR ALL COLLISION OPTIONS
C       WHEN INSTANTANEOUS DIFFUSION IS SPECIFIED.
C       IF THIS CAUSES PROBLEMS IT CAN BE RESTORED.
C
c
c       This was restored so that for collison option 1 -
c       which specifies no diffusion the RCONST value
c       would be large enough to prevent it from occurring
c       Even if instantaneous diffusion was specified.
c
        IF (CIOPTB.NE.1)  RCONST = 0.0
C
C          RCONST = 0.0
c
        ELSEIF (CDIFOP.EQ.1) THEN
          IF (LFPS(IK,IR,IZ).GT.0.0) THEN
            NRAND = NRAND + 1
            CALL SURAND2 (SEED, 1, RAN)
            RCONST = -TEMI * LOG (RAN) / LFPS(IK,IR,IZ)
          ENDIF
        ENDIF
C
C------ ZERO COUNTERS
C
        DIFFUS = .FALSE.
        ZERO_SPARA  = 0.0
        QUANT  = 0.0
c
c        DIST   = DBLE (CIST)
c
c        distacc = dble(cist)
c        cistacc = cist
C
        IF (DEBUGL) THEN
          WRITE (6,9005)
c
          WRITE (6,9003) IMP,0.0,IK,IR,IZ,R,Z,S,K,
     >      THETA,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,'ION APPEARED'
c
        ENDIF
        IF (100*(IMP/100).EQ.IMP)
     >    WRITE (6,'('' DIV: ION'',I6,'' STARTING'')') IMP
C
c
c       Calculate starting region and set appropriate logical
c       variable and update the appropriate storage array. Then
c       as the particle enters each of the other regions - set
c       the logical variable to true and accumulate the count
c       in the appropriate starting bin.
c
        if (global_hc_follow_option.ne.0.and.
     >      (idtype.eq.2.or.idtype.eq.5)) then
c
           irstart = idatizs (imp,3)
           ikstart = idatizs (imp,4)
c
c          Transferred source terms and localization data  
c          based on entire molecule 
c
           incore = Travel_Locations (imp,1)
           inedge = Travel_Locations (imp,2)
           inmsol = Travel_Locations (imp,3)
           indiv  = Travel_Locations (imp,4)
           intrap = Travel_Locations (imp,5)
c
        else
c
           irstart = ir
           ikstart = ik
c
           incore = .false.
           inedge = .false.
           inmsol = .false.
           indiv  = .false.
           intrap = .false.
        endif
c
        if (irstart.lt.irsep) then
c
c          Add checking to see if the flags have been previously set
c
           ! Updating num_entered_core is done by DIVIMP for each HC.
           num_entered_core = num_entered_core + sputy
c
           ! jdemod - wtsource data for each HC particle is accumulated
           !          in the HC code ... continued for ions here  
           !
           if (.not.incore) then  
              incore = .true.
c
              ncore(ikstart,irstart) = ncore(ikstart,irstart) + sputy
c
              if (idtype.gt.0) then
                  wtsource(idstart,irstart,3,idtype) =
     >                    wtsource(idstart,irstart,3,idtype) + sputy
                  wtsource(idstart,irstart,4,idtype) =
     >               wtsource(idstart,irstart,4,idtype)+sputy*sstart
              endif
c
           endif
c
        else
c
           if (.not.inedge) then
c
c             In the edge so set some sub-divisions
c
              inedge = .true.
              nedge(ikstart,irstart) = nedge(ikstart,irstart) + sputy
c
              if (ir.gt.irwall.and.ir.le.nrs) then
c
c                Particle starting in trap region
c
                 intrap = .true.
                 ntrap(ikstart,irstart) = ntrap(ikstart,irstart) + sputy
              elseif ((z.ge.zxp.and.refct.eq.1).or.
     >                (z.le.zxp.and.refct.eq.0)) then
c
c                Divertor region
c
                 indiv = .true.
                 ndivert(ikstart,irstart) = ndivert(ikstart,irstart)
     >                                      +sputy
              else
c
c                Main SOL Region
c
                 inmsol = .true.
                 nmsol(ikstart,irstart) = nmsol(ikstart,irstart)+sputy
              endif
           endif
        endif
c
        CFLRXA = .TRUE.
        CFLREX = .TRUE.
        IF     (IR.GE.IRSEP) THEN
          INMAIN = .FALSE.
          CFLRIN = .TRUE.
        ELSE
          INMAIN = .TRUE.
          CFLRIN = .FALSE.
          CICRNJ = CICRNJ + SPUTY
          IF (CSTOP.EQ.1) THEN
            stopped_follow = stopped_follow+sputy
            IFATE = 7
            GOTO 790
          ENDIF
        ENDIF
c
c       Check for PROMPT REDEPOSITION of the ORIGINAL ION
c
        if (prompt_depopt.eq.1.or.prompt_depopt.eq.2) then
c
c          Check for prompt deposition
c
           call promptdep(ik,ir,id,r,z,riz,sputy,crmi,temi,
     >                    sheath_fraction,rc)
c
c          A return code of 1 indicates that prompt redeposition
c          has occurred. The routine also returns the impact energy
c          of the depositing ion.
c

           if (rc.eq.1) then
c
c ---        REFLECT ---- IONS may need to be reflected here
c
c            Set default non reflection condition
c
             reflect_ion = .false.
c
             if (kmfss(id).le.-99.0) then
c
c               Probability of reflection starts at 0.0 for a value
c               of -99.0 and rises to 1.0 for a value of -100.0 or more
c
                refprob = min(abs(kmfss(id)+99.0),1.0)
c
c               Select random number
c
                NRAND = NRAND + 1
                CALL SURAND2 (SEED, 1, RAN)
c
                if (ran.le.refprob) reflect_ion = .true.
c
             endif
c
c            Test for reflection
c
             if (reflect_ion) then
c
c               Call Launch_one to launch a single reflected neutral -
c               then branch to the appropriate location depending
c               on the result.
c
c
c               Follow reflected impurities
c
                ikorg = ik
                irorg = ir
c
                write(6,*) ' debug: launch_one from div',id
                call LAUNCH_ONE (IMP,R,Z,RIZPOS,ZIZPOS,id,iwstart,
     >                   rc,ctem1,cist,sputy,
     >                   refSTRUK,mtcrefstruk,refMAIN,refEXIT,
     >                   refATIZ,refNEUT,refWALLN,mtcrefwalln,
     >                   refCENT,refTMAX,
     >                   SEED,NRAND,
     >                   NEUTIM,refFAIL,7,vout,vrec,refloss)
c
c
c               Set vel to appropriate value by scaling the along the field
c               line vout value returned by launch_one
c
                vel = vout * qtim
c
c
                write (6,'(a,2i8,6g12.5)') 'REFLECT IMP:',imp,rc,
     >                       r,z,rizpos,zizpos,
     >                       temi,cist
c
c               For all results other than reionization to state 1 - the code
c               will stop following the particle at this point and exit
c               as it would have for the old recombination implementation.
c
                if (rc.ne.5) goto 790
c
c               Deal with particle that was re-ionized
c
                r = rizpos
                z = zizpos
                iz= 1
c
c               Look near last recorded ik,ir
c
                call gridpos(ik,ir,r,z,.false.,griderr)
c
c               If a griderr then exit anyway and issue error message
c               since this shouldn't happen.
c
                if (griderr) then

                   write (6,*) 'PTR ION ERROR: NOT ON GRID:',r,z
                   goto 790
c
                else
c
c                  If on grid - re-assign S-value - to cell centre if
c                  it left the cell - otherwise leave as is.
c
                   if (ik.ne.ikorg.or.ir.ne.irorg) then
c
c                     Also need to reset the SMAX value for the new ring
c
                      smax = ksmaxs(ir)
c
                      if (init_pos_opt.eq.0) then
                         s = kss(ik,ir)
                         cross = 0.0
                      elseif (init_pos_opt.eq.1) then
c
                         call getscross_approx(r,z,s,cross,ik,ir)
c
                      endif
c
                   endif
c
                endif
c
c               Continue as if the particle had not struck the target
c
             else
c
c            ION is NOT reflected - continue as normal
c
c
c              Check for self-sputtering
c
c              Record all the regular loss statistics
c
               CICABS(IZ) = CICABS(IZ) + SPUTY
c
c               CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c               CILABS(IZ) = MAX (CILABS(IZ), CIST)
c
               CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
               CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
               CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c              Recording addition to total elapsed time in state IZ
c
               cieizs(iz) = cieizs(iz) + cistiz * sputy
               citizs(iz) = citizs(iz) + sputy
c
               CRTABS(IZ) = CRTABS(IZ) + TEMI * SPUTY
               CRVABS(IZ) = CRVABS(IZ) + VEL * SPUTY
               CRAVAV(IZ) = CRAVAV(IZ) + ABS(VEL) * SPUTY
               CTBS  (IZ) = CTBS  (IZ) + KTEBS(IK,IR) * SPUTY
               IM         = MIN (INT(TEMI/(0.2*CTEB0))+1, 10)
               CTEXS(IM)  = CTEXS(IM) + TEMI * SPUTY
c
               acttarg(iz) = acttarg(iz) + sputy
C
               RDEP   = RDEP + SPUTY
c
               if (id.lt.1.or.id.gt.nds) then
                  write (6,*) 'DEPS Error:',id,iz,sputy
               else
                  DEPS(ID,IZ) = DEPS(ID,IZ) + SPUTY
                  NEROS(ID,1) = NEROS(ID,1) + SPUTY
               endif
c
c              Prompt deposition energy is returned by the promptdep
c              subroutine - it is calculated as some fraction of the
c              3 kTe sheath/MPS drop - there are no contributions
c              for ion velocity. This is because the ion should
c              have had no time for energy transfer collisions and
c              should simply have its' own initial energy.
c
               ENERGY = sheath_fraction * RIZ * KTEBS(IK,IR) +
     >           5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c
c              energy = sheath_fraction*riz*ktebs(ik,ir) + 2.0*TEMI
c
c              Record average energy
c
               promptdeps(id,5) = promptdeps(id,5) + sputy * energy
c slmod begin
               if (allocated(wall_flx)) then
c                 in = nimindex(id)   ! Changed from NIMINDEX to WALLINDEX since the former is only assigned if running PIN. -SL, 26/03/2012
                 if (nimindex(id).ne.0) then
                   if (nimindex(id).ne.wallindex(id)) then
                     write(0,*) 'error: nimindex and wallindex are '//  ! temporary check
     .                          'not the same, investigate'
                     stop
                   endif
                 endif
                 in = wallindex(id)   
                 wall_flx(in)%prompt = wall_flx(in)%prompt + sputy
               endif
c slmod end
c
               if (kmfss(id).ge.0.0) then
                  RYIELD = YIELD (6, MATTAR, ENERGY,
     >                     ktebs(ik,ir),ktibs(ik,ir)) * KMFSS(ID)
               elseif (kmfss(id).lt.0.0.and.kmfss(id).ge.-50.0) then
                  RYIELD = abs(KMFSS(ID))
               elseif (kmfss(id).le.-99.0) then
                  RYIELD = YIELD (6, MATTAR, ENERGY,
     .                            ktebs(ik,ir),ktibs(ik,ir))
               endif
c
               SPUNEW = SPUTY * RYIELD
               YLDTOT = YLDTOT + SPUNEW
               YLDMAX = MAX (YLDMAX, SPUNEW)
c
c              WRITE(6,*) 'SPUTTERED:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID),
c     >                   ryield,kmfss(id),
c     >                    mattar,energy,sputy,spunew,yldtot
C
               IF (SPUNEW.GT.CTRESH) THEN
                  YTHTOT = YTHTOT + SPUNEW
                  NPROD  = NPROD + 1
                  SNEWS(NPROD) = SPUNEW
                  IF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.
     >                CNEUTD.EQ.4) THEN
                    EMAX = CEMAXF * ENERGY
                    RMAXS(NPROD) = 1.0 / (1.0+CEBD/EMAX)**2
                  ELSE
                    RMAXS(NPROD) = 1.0
                  ENDIF
c
c                 Find target impact coordinates for promptly redeposited particles.
c
                  IF (S.LE.smax/2.0) THEN
                     S  = 0.0
                     IK = 1
                     id = verify_id(ik,ir,2)
                  ELSE
                     S  = SMAX
                     IK = NKS(IR)
                     id = verify_id(ik,ir,1)
                  ENDIF

c
c                 No adjustments for initial position or position on target options
c                 since promptly redeposited particles are assumed to have not
c                 travelled too far.
c
                  R = RP(ID)
                  Z = ZP(ID)
c
                  XPRODS(NPROD) = R
                  YPRODS(NPROD) = Z
c
c                 For segments with a fixed sputtering yield - allow for
c                 the energy of the sputtered particle to be set to a
c                 specific value.
c
                  if (cselfs.eq.2.and.
     >               (kmfss(id).lt.0.0.and.kmfss(id).ge.-50.0))
     >               then
                     eprods(nprod) = ctem1
                  else
                     eprods(nprod) = 0.0
                  endif
c
                  IDPRODS(NPROD) = ID
                  launchdat(nprod,2) = 0.0
c
c                 Record self-sputtering event
c
                  promptdeps(id,6) = promptdeps(id,6) + spunew
c
               ENDIF
c
c
c              jdemod - added update of wall deposition in the case of prompt deposition
c
c              Update wall deposition
c
               call update_walldep(ik,ir,iz,id,0,iwstart,idtype,sputy,
     >                             energy)

c
c              Exit due to prompt deposition
c
               ifate = 10

! ammod begin.
               ! WBC comparison addition for ion prompt deposition.
               call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
! ammod end.
               goto 790

             endif

           endif

        endif
C
C-----------------------------------------------------------------------
C       CONTINUATION POINT IN MAIN LOOP   ...   RANDOM NUMBERS
C-----------------------------------------------------------------------
C
  500   CONTINUE
c
c       Check to see if more random numbers are required
c
        if (debug0) write(0,*) 'Before GRN',kk,kklim
c
        call get_random_numbers(kk,kklim,nrand,seed)
c
        if (debug0) write(0,*) 'After  GRN',kk,kklim
c
        DEBUG = DEBUGL.AND.(CIST.GE.TSTEPL)

c
        call change_local_values
c
c       Update ion temperature
c
        TEMI = MAX (LO, TEMI + (KTIBS(IK,IR)-TEMI) * LFTS(IK,IR,IZ))
c
c       Execute transport step
c
        if (debug_all)
     >     write(6,*) 'DEBUG A:',s,smax,slast,theta,cross
c
        if (debug0) write(0,*) 'Before EX',ifate,ik,ir,iz,s,r,z
c
        call execute_transport_step(seed,nrand,neutim,ero_record_data)
c
        if (debug_all)
     >     write(6,*) 'DEBUG B:',s,smax,slast,theta,cross

c
        if (debug0) write(0,*) 'After EX',ifate,ik,ir,iz,s,r,z

c
        if (ifate.ne.0) goto 790
c
C
C-----------------------------------------------------------------------
C       BOTH ROUTES CONTINUE HERE ...
C-----------------------------------------------------------------------
c
c
c        if (kk.gt.(kklim+10)) then
c           write(6,*) '3:',kk,kklim
c        endif
c
c-------------------------------------------------------------------c
c
c
C
C------ CHECK FOR COLLISION
C
        KK = KK + 1
        IF ((TEMI*RANV(KK)) .LE. LFPS(IK,IR,IZ)) THEN
          CICCOL = CICCOL + SPUTY
        ENDIF
C
C------ SCORE PARTICLE IN ARRAYS
C
        DDLIMS(IK,IR,IZ) = DDLIMS(IK,IR,IZ) + DSPUTY

        DDTS  (IK,IR,IZ) = DDTS  (IK,IR,IZ) + DSPUTY * DBLE(TEMI)
c
        if (subgrid_opt.gt.0) then
           call getrz(ik,ir,s,cross,r,z,rzopt)
           call update_subgrid(r,z,iz,sputy)
        endif

c
        if (ddlims(ik,ir,iz).lt.0.0) then
           write (6,'(a,3i4,2g16.8)') 'DDLIM Error:',ik,ir,iz,
     >                          ddlims(ik,ir,iz),dsputy

        endif
c
c       Determine the time bin for the particle
c
        if (nts.gt.0) then
c
           IT = IPOS (sngl(CIST), CTIMES(1,IZ), NTS+1)
c
c          IF (CIST.GE.CTIMES(IT,IZ)) THEN
c
c          Record time dependent contribution at each time step
c
           LIMS(IK,IR,IZ,IT) = LIMS(IK,IR,IZ,IT) + SPUTY
c
c          IT = IT + 1
c          ENDIF
C
        endif
c
c
c       (RIV)
c
        call debug_velocity
c
c
c
        XXX(IZ) = MIN (XXX(IZ), Z)
        SSS(IZ) = MAX (SSS(IZ), MIN (S,SMAX-S))
        KKK     = KKK + K
        IF (CFLRIN) THEN
          IF (K.GT.KATIZS(IMP)) THEN
            XTRIPP  = XTRIPP + (Z - OLDZ)
          ELSE
            XTRIPS  = XTRIPS + (Z - OLDZ)
          ENDIF
        ENDIF
        OLDZ = Z
C
        IF (DIFFUS) THEN
          IF (ZERO_SPARA.LE.0.0) THEN
            DPARAS(IZ,1) = DPARAS(IZ,1) + DSPUTY
            IF (IR.GE.IRSEP) DPARAS(IZ,5) = DPARAS(IZ,5) + DSPUTY
          ELSE
            DPARAS(IZ,2) = DPARAS(IZ,2) + DSPUTY
            IF (IR.GE.IRSEP) THEN
              DPARAS(IZ,3) = DPARAS(IZ,3) + DSPUTY
              DPARAS(IZ,4) = DPARAS(IZ,4) + DSPUTY * DBLE(ZERO_SPARA)
            ENDIF
          ENDIF
        ENDIF

c
c       Check for a change of state
c

        if (debug0) write(0,*) 'Before CS',ifate,iz

        call check_ion_change_state(seed,nrand,neutim,nizs)

        if (debug0) write(0,*) 'After CS',ifate,iz

c
        if (ifate.ne.0) goto 790
c
c       Check for ion removal
c
        if (debug0) write(0,*) 'Before IR',ifate,iz

        call check_ion_removal(ifate,kk,ik,ir,iz,cist,cistiz,
     >                             ssss,s,smax,sputy)
c
        if (debug0) write(0,*) 'After IR ',ifate,iz,cist,debug

        if (ifate.ne.0) goto 790
C
C-----------------------------------------------------------------------
C       LOOP BACK IF CUTOFF TIME NOT YET REACHED
C-----------------------------------------------------------------------
C
        IF (DEBUG) THEN

         if (debug0)  write(0,*) 'IN Debug ',tstepl,cstepl,cist


  495     TSTEPL = TSTEPL + CSTEPL
          IF (TSTEPL.LE.CIST) GOTO 495
c
c
c nonorth
          WRITE (6,9003) IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >      THETA,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,'SET K, NEW S'

          if (debug0)
     >     WRITE (0,9003) IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >      THETA,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,'SET K, NEW S'
c nonorth
c          WRITE (6,9003) IMP,CIST,IK,IR,IZ,
c     >      R,Z,S,K,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,'SET K, NEW S'
        ENDIF
C
c       Record walks
c
        if (debug0) write(0,*) 'Before RW',iw,maxnws,r,z
c

        IF (IW.LT.MAXNWS) THEN
          WALKS(IW,1) = R
          WALKS(IW,2) = Z
          WALKS(IW+1,1) = HI
          WALKS(IW+1,2) = HI
          IW = IW + 1
        ENDIF

        if (debug0) write(0,*) 'Before TI',cist,cstmax

C
c       Increment time - check to see if less than maximum
c
        IF (CIST.LT.CSTMAX) THEN
c
c         Increment timing and go to beginning of particle iteration
c
          cist = cist + 1.0d0
c
c         Update time in ionization state IZ
c
          cistiz = cistiz + 1.0d0
c
c         Continue following particle
c
          if (debug0) write(0,*) 'Before LB',cist,cstmax

          GOTO 500

        ENDIF
C
C------ REACHED CUTOFF TIME
C
 780    CONTINUE
c        CICUTS(IZ) = CICUTS(IZ) + SPUTY
c        CRTRCS(IZ) = CRTRCS(IZ) + TEMI * SPUTY
c        TCUT = TCUT + SPUTY
c slmod begin
        IF (sloutput) IFATE = 3
c
cc        IFATE = 3
c slmod end
C      (GOTO 790)
C
C-----------------------------------------------------------------------
C       CURRENT ION OR SUB-ION FINISHED WITH ... END OF DO-LOOP
C-----------------------------------------------------------------------
C
  790   CONTINUE
        if (ifate.eq.3) then
           CICUTS(IZ) = CICUTS(IZ) + SPUTY
           CRTRCS(IZ) = CRTRCS(IZ) + TEMI * SPUTY
           TCUT = TCUT + SPUTY
c          IFATE = 3
c slmod begin - t-dep
c...       Store the state of the ion so that it can be re-launched in a 
c          subsequent run:                     
           IF (OPT_DIV%PSTATE.EQ.1) THEN
             IF (.NOT.ALLOCATED(TDEP_SAVE)) ALLOCATE(TDEP_SAVE(NATIZ))

             tdep_save_n = tdep_save_n + 1
             tdep_save(tdep_save_n)%r      = r
             tdep_save(tdep_save_n)%z	   = z
             tdep_save(tdep_save_n)%phi    = 0.0
             tdep_save(tdep_save_n)%s	   = s
             tdep_save(tdep_save_n)%cross  = cross
             tdep_save(tdep_save_n)%diag   = 0.0
             tdep_save(tdep_save_n)%temp   = temi
             tdep_save(tdep_save_n)%vel    = vel
             tdep_save(tdep_save_n)%charge = riz
             tdep_save(tdep_save_n)%weight = sputy
           ENDIF
c slmod end
        endif

c
        CISTOT = CISTOT + CIST * SPUTY
        CISMAX = MAX (CISMAX, CIST)
c
c       Record recombination counts
c
        if (reccnt.gt.10) then
           rectotcnt(11) = rectotcnt(11) + sputy
        else
           rectotcnt(reccnt) = rectotcnt(reccnt) + sputy
        endif
c
c       Reset time in ionization state IZ
c
        cistiz = 1.0
c
        DO 792 JZ = CIZSC, MAXCIZ
          RIONS(JZ) = RIONS(JZ) + SPUTY
          CXXX (JZ) = CXXX (JZ) + XXX(JZ) * SPUTY
          CSSS (JZ) = CSSS (JZ) + SSS(JZ) * SPUTY
          CSSSS(JZ) = CSSSS(JZ) + SSSS(JZ) * SPUTY
  792   CONTINUE
C
        IF (IW.LT.MAXNWS) THEN
          WALKS(IW,1) = R
          WALKS(IW,2) = Z
          WALKS(IW+1,1) = HI
          WALKS(IW+1,2) = HI
          IW = IW + 1
        ENDIF
C
        IF (CFLRIN) THEN
          CICRNO = CICRNO + SPUTY
          CIKRNO = CIKRNO + KATIZS(IMP) * SPUTY
          cirrno = cirrno + xatizs(imp) * sputy
          cizrno = cizrno + yatizs(imp) * sputy
          cisrno = cisrno + satizs(imp) * sputy
          CKTRNO = CKTRNO + KKK / CIST * SPUTY
        ENDIF
C
c nonorth
        IF (DEBUGL) WRITE (6,9003) IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >    THETA,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,FATE(IFATE)
c slmod begin
        IF (DEBUGL) WRITE(6,*) '  ITERATION LIMIT',cstmax
c slmod end
c nonorth
c        IF (DEBUGL) WRITE (6,9003) IMP,CIST,IK,IR,IZ,
c     >    R,Z,S,K,SMAX,VEL,TEMI,ZERO_SPARA,CROSS,SPUTY,IT,FATE(IFATE)
C
C
C        THE FOLLOWING TIME LIMIT CODE IS TAKEN DIRECTLY FROM LIM WHERE
C        IT WORKS WITH REASONABLE SUCCESS. THE SAME IS EXPECTED HERE. IT
C        IS NECESSARY TO PREVENT RUN-AWAY CONDITIONS FROM DEVELOPING. ON
C        THE WORKSTATION LONG RUNS ARE ACCEPTABLE, SO A TIME LIMIT OF
C        36000.0 SECONDS MAY BE USED.
C
C-----------------------------------------------------------------------
C       SEE IF TEST OF CPU TIME USED IS DUE
C       TRAP CASE WHERE TIMUSD=0 OCCURS.
C-----------------------------------------------------------------------
C
        IF (IMP.GE.IMPLIM) THEN
         TIMUSD = ZA02AS(1) - STATIM
         IF (TIMUSD.GT.0.0) THEN
          PARTIM = (CPULIM-NEUTIM-IONTIM) / TIMUSD
         ELSE
           PARTIM = 10.0
         ENDIF
         IF (PARTIM.GE.1.05) THEN
           IMPLIM = INT (REAL(IMP) * MIN (4.0,0.25+0.75*PARTIM))
         ELSE
C
C--------- HAVE RUN OUT OF CPU TIME, STOP ITERATION
C--------- THERE ARE SO MANY COUNTERS IT IS VIRTUALLY IMPOSSIBLE TO
C--------- WIND UP THE ROUTINE CLEANLY.  JUST WORK OUT HOW MANY IONS
C--------- HAVE BEEN FOLLOWED IN THE TIME ALLOTTED AND CALL IT A DAY.
C
           IMPLIM = IMP
           DO 791 J = IMP+1, NATIZ
             RATIZ = RATIZ - SPUTYS(J)
  791      CONTINUE
           NATIZ = IMP
           WRITE (6,'('' ERROR:  CPU TIME LIMIT REACHED'')')
           WRITE (6,'('' NUMBER OF IONS REDUCED TO'',I5)') NINT(RATIZ)
           CALL PRB
           CALL PRC ('ERROR:  CPU TIME LIMIT REACHED')
           CALL PRI ('NUMBER OF IMPURITY IONS REDUCED TO ',NINT(RATIZ))
           CALL PRC ('INCREASE CPU TIME OR INCREASE QUANTUM TIMESTEP')
           CALL PRC ('RESULTS WILL BE PRINTED BUT THEY SHOULD BE TREATED
     > WITH CAUTION ...')
           GOTO 805
         ENDIF
        ENDIF
C
  800 CONTINUE
C
  805 CONTINUE

      call pr_trace('DIV','MAIN LOOP END')

C
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
C     M A I N   L O O P   E N D
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C 
C
C---- CALCULATE AVERAGE TB, NB VALUES
C
      DO 876 M = 1, 2
c
        ik = 0
        ir = 0
c
        if (avatiz(m).ne.0.0) then
           call gridpos(ik,ir,avxpos(m)/avatiz(m),avypos(m)/avatiz(m),
     >                  .true.,griderr)
        else
           griderr = .true.
        endif
c
        if (griderr) then
           AVTPOS(M) = -1.0
           AVNPOS(M) = -1.0
        else
           AVTPOS(M) = KTEBS(IK,IR)
           AVNPOS(M) = KNBS (IK,IR)
        endif
  876 CONTINUE
C
C---- PRINT DETAILS AND ADD TO TOTALS
C
      write(7,'(/40X,a,4x,a)')  INNER,OUTER
c
c      if (zxp.le.z0) then
c         WRITE (7,'(/40X,''OUTER    INNER'')')
c      else
c         WRITE (7,'(/40X,''INNER    OUTER'')')
c      endif
c
c     changed format to allow for larger numbers, Krieger IPP/97
      WRITE (7,'(1X,A,2(F10.3,1x))')
c
     >           'NUMBER OF IONS INJECTED            ',
     >           AVATIZ(1),AVATIZ(2)
      WRITE (7,'(1X,A,2(F10.5,1x))')
     >           'AVERAGE R INJECTION POSITION       ',
     >            AVXPOS(1)/AVATIZ(1),AVXPOS(2)/AVATIZ(2)
      WRITE (7,'(1X,A,2(F10.5,1x))')
     >            'AVERAGE Z INJECTION POSITION       ',
     >            AVYPOS(1)/AVATIZ(1),AVYPOS(2)/AVATIZ(2)
      CALL PRR2 ('PLASMA ELECTRON TEMP AT MEAN (R,Z) ',
     >                                       AVTPOS(1),AVTPOS(2))
      CALL PRR2 ('PLASMA DENSITY AT MEAN (R,Z)       ',
     >                                       AVNPOS(1),AVNPOS(2))
      WRITE (7,'('' AVERAGE K VALUE AT INJECTION'',F24.6)')AVKPOS/RATIZ
      CALL PRR('AVERAGE S OR SMAX-S VALUE AT INJECTION   ',AVSPOS/RATIZ)
      CALL PRR('AVERAGE SMAX VALUE AT INJECTION          ',AVSMAX/RATIZ)
      CALL PRR('AVERAGE ABS(VEL) AT INJECTION  M/S       ',AVVPOS/RATIZ)
      CALL PRR('AVERAGE TIME TO FIRST DIFFUSION          ',RDIFFT/RATIZ)
      CALL PRr('NUMBER OF IONS PLATING OUT ON WALLS      ',RWALL)
      CALL PRr('NUMBER OF IONS PLATING ON TARGET         ',RDEP)
      call prr('NUMBER OF IONS PLATING ON FP TARGET      ',rfptarg)
      call prr('NUMBER OF IONS ENTERING CORE PLASMA      ',
     >                                                 num_entered_core)
      CALL PRI('TOTAL SELF-SPUTTERING YIELD              ',NINT(YLDTOT))
      CALL PRI('TOTAL YIELDS ABOVE THRESHOLD YIELD       ',NINT(YTHTOT))
      CALL PRR('MAXIMUM YIELD ANY NEUTRAL FRAGMENT       ',YLDMAX)
      WRITE (6,'(1X,A,I9 )') 'AVERAGE NUMBER OF ITERATIONS PER ION   ',
     >                                         NINT(sngl(CISTOT)/RATIZ)
      WRITE (6,'(1X,A,I13)') 'MAXIMUM NUMBER OF ITERATIONS       ',
     >                                                     NINT(CISMAX)

      TATIZ  = TATIZ  + RATIZ
      TDEP   = TDEP   + RDEP
      TWALL  = TWALL  + RWALL
      TNEUT  = TNEUT  + RNEUT
      TWALLN = TWALLN + RWALLN
      TMTCWALLN = TMTCWALLN + MTCWALLN
      TCENT  = TCENT  + RCENT
      TTMAX  = TTMAX  + RTMAX
      TSTRUK = TSTRUK + RSTRUK
      TMTCSTRUK = TMTCSTRUK + MTCSTRUK
      TFAIL  = TFAIL  + RFAIL
      TMAIN  = TMAIN  + RMAIN
      TEXIT  = TEXIT  + REXIT
      IONTIM = IONTIM + ZA02AS (1) - STATIM
C
C---- LAUNCH SECONDARIES, TERTIARIES, ETC IF REQUIRED
C
      IF ((cselfs.eq.1.or.cselfs.eq.2).AND.
     >     NPROD.GT.0 .AND.
     >     YTHTOT.GT.2.0 .AND. (STATUS.LT.CMAXGENS)) THEN
        STATUS = STATUS + 1
        CNEUTA = 0
        CNEUTB = 0
        IF (CNEUTD.EQ.4) CNEUTC = 1
        IF (CNEUTE.EQ.1) CNEUTE = 0
        IF (STATUS.LE.50) THEN
           WRITE(6,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'
           WRITE(7,9012) '***  LAUNCHING ',WHAT(STATUS),' NEUTRALS  ***'
        ELSE
           WRITE (6,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'
           WRITE (7,9013) '***  LAUNCHING ',STATUS,
     >                    ' GENERATION NEUTRALS  ***'
        ENDIF
        DO 8700 IMP = 1, NPROD
          SPUTYS(IMP) = SNEWS(IMP)
 8700   CONTINUE


! ammod begin.
c jdemod - some code changes required to make this option work - likely create a global
c          version of hc_self_sputter
c
c
c        If (global_HC_Follow_Option .eq. 0 .or. HC_Self_Sputter .eq. 0) Then
c          ! Launch self-sputtered particles via the DIVIMP method.
c
        CALL LAUNCH (1,NPROD,1,NATIZ,RSTRUK,MTCSTRUK,RMAIN,REXIT,
     >               RATIZ,RNEUT,RWALLN,MTCWALLN,RCENT,RTMAX,SEED,NRAND,
     >               NEUTIM,RFAIL,STATUS,6,MATTAR,3)

c
c        Else
c          ! Launch self-sputtered particles via DIVIMP-HC.
c          ! Note:  These will only be as a result of C+ striking
c          ! the target region and self sputtering (all CHx self
c          ! sputtering has already been taken care of).  Start
c          ! HC_Launch with CNEUTB and CNEUTC as set above.
c           CALL HC_Launch (1,NPROD,1,NATIZ,RSTRUK,MTCSTRUK,RMAIN,REXIT,
c     >               RATIZ,RNEUT,RWALLN,MTCWALLN,RCENT,RTMAX,SEED,
c     >               NRAND,NEUTIM,RFAIL,STATUS,6,MATTAR,3,CNEUTB,CNEUTC)
c        End If
c
! ammod end.


C        WRITE(6,*) 'LAUNCHING: ',WHAT(STATUS),' NEUTRALS - ',NATIZ
        IF (NATIZ.GT.0) GOTO 200
      ENDIF
C
C-----------------------------------------------------------------------
C     ALL PARTICLES COMPLETED - PRINT SUMMARY
C-----------------------------------------------------------------------
C
      call pr_trace('DIV','PRINT SUMMARY START')

      goto 8701

c
c     Try to capture a SIGUSR1 signal as a process termination.
c     This is used to issue an error message and (hopefully) print out
c     whatever information is available.
c
      entry divkill
      procterm = .true.
      call prb
      call prb
      call prc('DIVIMP PROCESS TERMINATED BY USER SIGNAL!!!!')
      call prc('-------   RESULTS ARE NOT RELIABLE ---------')
      call prc('USE FOR DEBUGGING PURPOSES ONLY !!!!!!!!!!!!')
      call prb
      call prb
      write (6,*) 'DIVIMP PROCESS TERMINATED BY USER SIGNAL!!!!'
      write (6,*) '-------   RESULTS ARE NOT RELIABLE ---------'
      write (6,*) 'USE FOR DEBUGGING PURPOSES ONLY !!!!!!!!!!!!'
c
c     Continue with wrap up execution.
c
 8701 continue




      WRITE(6,*) 'FINISHED PARTICLES'
c
c     Print out the parallel diffusive step characterization.
c
c     changed by Krieger, IPP 12/94
c
      do 8702 in=1,6
c
         dvparastep(in,nizs+1) = 0.0
         vparastep(in,nizs+1) = 0.0
         dvparacnt(in,nizs+1) = 0.0
c
         do iz = 1,nizs
c
            if (abs(dvparacnt(in,iz)).gt.1.e-14) then
              dvparanorm(in,iz)=dvparastep(in,iz)/dvparacnt(in,iz)
              vparanorm(in,iz) = vparastep(in,iz)/dvparacnt(in,iz)
            else
              dvparanorm(in,iz)=0.
              vparanorm(in,iz)=0.
            endif
c
            dvparastep(in,nizs+1)=dvparastep(in,nizs+1)
     >                           +dvparastep(in,iz)
            vparastep(in,nizs+1) = vparastep(in,nizs+1)
     >                           + vparastep(in,iz)
            dvparacnt(in,nizs+1) =dvparacnt(in,nizs+1)
     >                           +dvparacnt(in,iz)
c
         end do
c
         if (abs(dvparacnt(in,nizs+1)).gt.1.e-14) then
            dvparanorm(in,nizs+1)=dvparastep(in,nizs+1)
     >                            /dvparacnt(in,nizs+1)
            vparanorm(in,nizs+1) = vparastep(in,nizs+1)
     >                            /dvparacnt(in,nizs+1)
         else
            dvparanorm(in,iz)=0.
            vparanorm(in,iz)=0.
         endif
c
         if (abs(dsparacnt(in)).gt.1.e-14) then
            dsparanorm(in)=dsparastep(in)/dsparacnt(in)
         else
            dsparanorm(in)=0.
         endif
c
 8702 continue

c
c     changed by Krieger, IPP 12/94
c
      do iz = 1,nizs+1
        write(6,*)   'Summary of Velocity Diffusive steps: State=',iz
        write(6,2151)
        write(6,2152) (dvparacnt(in,iz),in=1,6)
        write(6,2152) (dvparanorm(in,iz),in=1,6)
        write(6,2152) (vparanorm(in,iz)/qtim,in=1,6)
        if (iz.ne.nizs+1) then
           write(6,2152) (dvmaxv(in,iz)/qtim,in=1,6)
           write(6,2152) (dvminv(in,iz)/qtim,in=1,6)
        endif
      end do
c
      write(6,*)   'Summary of Spatial Diffusive steps'
      write(6,2151)
      write(6,2152) (dsparacnt(in),in=1,4)
      write(6,2152) (dsparanorm(in),in=1,4)
c
 2151 format(5x,'OUTSOL-',8x,'OUTSOL+',8x,'INSOL -',8x,'INSOL +',
     >8x,'CORE - ',8x,'CORE + ')
 2152 format(6g15.6)
c
c
c
      SSEF = 0.0
      YEFF = 0.0
      IF (RNEUT1.GT.0.0) SSEF = (TNEUT-RNEUT1) / RNEUT1
      IF (FTOT.GT.0.0)   YEFF = (1.0 + SSEF) * FYTOT / FTOT
C
  806 CONTINUE

! ammod begin.
      ! Complete launch, follow and recording of all particles.
      ! Print hydrocarbon diagnostics.
      If (global_hc_follow_option .ne. 0) Then
!        write (0,*) "Doing HC End a"
        Call global_HC_End (NIMPS,NIMPS2,REXIT,RMAIN,RNEUT,NEUTIM,NIZS)
      End If
! ammod end.


c
! ammod begin
c
c     Print out some HC simulation related data (?)
c
c      call global_hc_print_temp
c
c
c      open (unit=91,file="temp.out")
c      ! Write rate DDLIMS data for injir
c      Write (91,17) "Final CIST",CIST,"CISTIZ",CISTIZ,"CISTOT",CISTOT,
c     >  "INJIR",INJIR
c      Do M=1,nks(injir)
c         write (91,19) "cell",M,"counts",DDLIMS(M,injir,1),"length",
c     >     ksb(M,injir)-ksb(M-1,injir),"Center S",kss(M,injir),
c     >     "KBFS",kbfs(M,injir),"BRatio",
c     >     bratio(M,injir),"area",kareas(M,injir),
c     >     "ratio L/A",(ksb(M,injir)-ksb(M-1,injir))/kareas(M,injir),
c     >     "ratio C/A",DDLIMS(M,injir,1)/kareas(M,injir),
c     >     "ratio C/L",DDLIMS(M,injir,1)/(ksb(M,injir)-ksb(M-1,injir)),
c     >     "ne",knbs(M,injir),"te",ktebs(M,injir),"ti",ktibs(M,injir),
c     >     "vb",kvhs(M,injir),"e",kes(M,injir)
c      End Do
c      write (91,*) "By Ring:",irwall
c      Do N=1,irwall-1
c       ! Find cell at injir S.
c
C       FOR NOW SIMPLY FIND GRID POINT CLOSEST TO DESIRED
C       INJECTION POSITION. THIS IS ALL THAT IS CURRENTLY DONE.
C       THE PARTIAL DISPLACEMENTS FROM THE GRID POINTS ARE IGNORED.
C
c        NRAND = NRAND + 1
c        CALL SURAND2 (SEED, 1, RANVAL)
c        STMP = RANVAL*(INJF2-INJF1)*KSMAXS(N)
c     >         + INJF1*KSMAXS(N)
c        NRAND = NRAND + 1
c        CALL SURAND2 (SEED, 1, RANVAL)
c        IF (RANVAL.GT.0.5.and.(ciopte.eq.2.or.ciopte.eq.5)) THEN
C
C         OUTER PLATE - OTHERWISE INNER - for option 2
C
c          STMP = KSMAXS(N) - STMP
c        ENDIF
c
c        ! FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR INJI
c        IK = 1
c  758   IF (IK.LT.NKS(N).AND.STMP.GT.KSS(IK,N)) THEN
c          IK = IK + 1
c          GOTO 758
c        ENDIF
c  759   IF (IK.GT.1.AND.STMP.LE.KSS(IK-1,N)) THEN
c          IK = IK - 1
c          GOTO 759
c        ENDIF
c        IF (IK.GT.1.AND.
c     >     (STMP-KSS(IK-1,N).LT.KSS(IK,N)-STMP))
c     >      IK = IK - 1
c        M = IK
c
c        If (N .lt. irsep) Then
c          M = 14
c        Else
c          M = 34
c        EndIf
c
c         write (91,18) "Ring",N,"cell",M,"counts",DDLIMS(M,N,1),
c     >     "length",ksb(M,N)-ksb(M-1,N),"Center S",kss(M,N),
c
c     >     "r centre",rs(M,N),"z centre",zs(M,N),
c     >     "KBFS",kbfs(M,injir),"BRatio",
c     >     bratio(M,N),"area",kareas(M,N),
c     >     "ratio L/A",(ksb(M,N)-ksb(M-1,injir))/kareas(M,injir),
c     >     "ratio C/A",DDLIMS(M,N,1)/kareas(M,N),
c     >     "ratio C/L",DDLIMS(M,N,1)/(ksb(M,N)-ksb(M-1,N))
c      End Do
c   17 Format (a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,I3)
c   !18 Format (a,1x,I3,1x,a,1x,I3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x,
c   !  >  a,1x,e9.3,1x,a,1x,e9.3,1x,
c   !  >  a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,f6.4,1x,a,1x,e9.3,1x,
c   !  >  a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x)
c   !19 Format (a,1x,I3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x,1x,
c   !  >  a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,f6.4,1x,a,1x,e9.3,1x,
c   !  >  a,1x,e9.3,1x,a,1x,e9.3,1x,a,1x,e9.3,1x)
c   18 Format (2(A,1X,I3,1X),20(A,1X,E11.5,1X))
c   19 Format (1(A,1X,I3,1X),20(A,1X,E11.5,1X))
c      close (unit=91)
! ammod end.

C
C     Calculate total leakage - from information for each charge
c     state.
C
      if (checkleak) then
        do 3003 in = 1,cleaksn
          DO 3003 IZ = 1,NIZS
            cleakn(in,nizs+1) = cleakn(in,nizs+1) + cleakn(in,iz)
3003    CONTINUE
      endif
C
C
C     CALCULATE THE WALL AND TRAP DEPOSITION TOTALS FOR ALL STATES
C
      DO 3000 IZ = 0,MAXIZS
        DO 3000 IR = 1,NRS
          DO 3000 IK = 1,NKS(IR)
             WALLS(IK,IR,MAXIZS+1) = WALLS(IK,IR,MAXIZS+1) +
     >                            WALLS(IK,IR,IZ)
3000  CONTINUE
C
C     Neutrals on Walls - also sum up wall deposition and erosion arrays
C
      write(6,'(a,5(1x,f9.2))') 'Walls Start:',
     >      wallse(maxpts+1),wallse_i(maxpts+1),
     >      wallsi(maxpts+1),wallsn(maxpts+1)
c
      do 3005 in = 1,wallpts
c
c         write (6,'(a,i5,3(1x,g12.5))') 'Wall Dep:',in,
c     >        wallsn(in),wallsi(in),wallse(in)
c
         walltotn = walltotn + wallsn(in)
         wallse(maxpts+1) = wallse(maxpts+1) + wallse(in)
         wallse_i(maxpts+1) = wallse_i(maxpts+1) + wallse_i(in)
         wallsn(maxpts+1) = wallsn(maxpts+1) + wallsn(in)
         wallsi(maxpts+1) = wallsi(maxpts+1) + wallsi(in)
         
         do iz= 1,nizs
            if (wallsiz(in,iz).gt.0.0) then 
               wallseiz(in,iz) = wallseiz(in,iz)/wallsiz(in,iz)
            endif
         enddo

         if (iz.lt.80) then 
            write(6,'(a,i7,3(1x,f12.6),256(1x,f9.3))') 
     >          'Walls Data:',in,
     >         wallpt(in,1),wallpt(in,2),wallpt(in,7),
     >         wallse(in),wallse_i(in),wallsi(in),wallsn(in),
     >         ((real(iz),wallsiz(in,iz),wallseiz(in,iz)),iz=1,nizs)
         else
            write(6,'(a,i7,256(1x,f9.2))') 'Walls Data:',in,
     >         wallpt(in,1),wallpt(in,2),wallpt(in,7),
     >         wallse(in),wallse_i(in),wallsi(in),wallsn(in)
         endif

3005  continue
c
      write(6,'(a,5(1x,f9.2))') 'Walls End:',
     >      wallse(maxpts+1),wallse_i(maxpts+1),
     >      wallsi(maxpts+1),wallsn(maxpts+1)


C     K. Schmid 2008 output charge state resolved wall impact information
      write (6, *) 'CHARGE RESOLVED WALL IMPACT INFO START: ', NIZS,
     >               wallpts
c
c     jdemod - the <NIZS> repeat specification does not appear to be 
c              standardized and breaks some compilers so I changed the 
c              repeat value to 100 which should be large enough in most
c              cases. 
c
3006  Format(i5,' ',g12.5,g12.5,100(' ',g12.5))
c3006  Format(i5,' ',g12.5,g12.5,<NIZS>(' ',g12.5))
c
      if (nizs.le.100) then 
         do in = 1,wallpts
C        write (6,*) in,' ',wallsn(in),' ',wallsi(in),' ',
C     >       wallsiz(in, 1:NIZS)
             write (6,3006) in,wallsn(in),wallsi(in),wallsiz(in, 1:NIZS)
         end do
         write (6,3006) -1, wallsn(maxpts+1),wallsi(maxpts+1),
     >      wallsiz(maxpts+1, 1:NIZS)
         write (6, *) 'END OF CHARGE RESOLVED WALL IMPACT INFO'
      else
         call errmsg('ERROR PRINTING CHARGE'//
     >               ' RESOLVED WALL IMPACT INFO: NIZS > 100')
      endif
c
c     jdemod end
c

      call pr_trace('DIV','CALCULATE TOTALS')

C
C
C     CALCULATE TOTALS
C
      DO 3010 IZ = 0,MAXIZS+1
        DO 3010 IK = 1,NKS(IRWALL)
         TNTOTS(IZ,1) = TNTOTS(IZ,1) + WALLS(IK,IRWALL,IZ)
3010  CONTINUE
c
c     Sum up wall loss particles that may not have been in
c     wall ring at the time of the particle loss
c
      DO IZ = 0,MAXIZS+1
         do ir = 1,irwall-1
            DO IK = 1,NKS(IR)
               TNTOTS(IZ,3) = TNTOTS(IZ,3) + WALLS(IK,IR,IZ)
            end do
         end do
      end do
c
      if (cgridopt.eq.2) then
         DO 3011 IZ = 0,MAXIZS+1
           DO 3011 IK = 1,NKS(IRWALL2)
            TNTOTS(IZ,1) = TNTOTS(IZ,1) + WALLS(IK,IRWALL2,IZ)
3011     CONTINUE
      endif
C
      DO 3020 IZ = 0,MAXIZS+1
        DO 3020 IK = 1,NKS(IRTRAP)
         TNTOTS(IZ,2) = TNTOTS(IZ,2) + WALLS(IK,IRTRAP,IZ)
3020  CONTINUE
c
c     Sum up wall loss particles that may not have been in
c     wall ring at the time of the particle loss
c
      DO IZ = 0,MAXIZS+1
         do ir = irtrap+1,nrs
            DO IK = 1,NKS(IR)
               TNTOTS(IZ,4) = TNTOTS(IZ,4) + WALLS(IK,IR,IZ)
            end do
         end do
      end do
c
      if (cgridopt.eq.2) then
         DO 3021 IZ = 0,MAXIZS+1
           DO 3021 IK = 1,NKS(IRtrap2)
            TNTOTS(IZ,2) = TNTOTS(IZ,2) + WALLS(IK,IRtrap2,IZ)
3021     CONTINUE
      endif


      call pr_trace('DIV','PRINT TOTALS')


C
C     NOW PRINT THESE IN THE SUMMARIES
C
      CALL PRB
c      CALL PRC ('---  S U M M A R Y   D E T A I L S  ---')
      CALL PRChtml ('---  S U M M A R Y   D E T A I L S  ---',
     >            'pr_summary','0','B')
      CALL PRB
c
c     Convert to printing out real numbers so that exact values
c     for fractional weights are printed.
c
      call prb
c      call prc ('ORIGINAL NEUTRAL PARTICLE SUMMARY:')
      call prchtml('ORIGINAL NEUTRAL PARTICLE SUMMARY:',
     >             'pr_neut','0','B')
      call prb
c
      CALL PRr0('TOTAL NUMBER OF NEUTRALS LAUNCHED        ',TNEUT)
      CALL PRr0('TOTAL NO OF NEUTRALS PLATING ON WALL+TRAP',TWALLN)
C
      if (mtcopt.eq.1) then
         call prr0('TOTAL NO OF THESE NEUTRALS HAVING MTC    ',
     >                                                   TMTCWALLN)
      ENDIF
C
      CALL PRr2('TOTAL NO OF NEUTRALS PLATING ON WALLS    ',
     >              TNTOTS(0,1),tntots(0,3))
      CALL PRr2('TOTAL NO OF NEUTRALS PLATING ON TRAP     ',
     >              TNTOTS(0,2),tntots(0,4))
      CALL PRr0('TOTAL NO OF NEUTRALS STRIKING TARGET     ',TSTRUK)
C
      if (mtcopt.EQ.1) THEN
         CALL PRr0('TOTAL NO OF THESE NEUTRALS HAVING MTC    ',
     >                                                   TMTCSTRUK)
      ENDIF
C
      CALL PRr0('TOTAL NO OF NEUTRALS STRIKING ALL SURFACES',
     >         walltotn)
      CALL PRr0('TOTAL NO OF NEUTS REACHING CENTRAL MIRROR',TCENT)
      CALL PRr0('TOTAL NO OF NEUTRALS EXISTING AT TMAX    ',TTMAX)
      CALL PRr0('TOTAL NO OF FAILED NEUTRAL LAUNCHES      ',TFAIL)
      CALL PRr0('TOTAL NEUTRALS ENTERING MAIN PLASMA      ',TMAIN)
      CALL PRr0('TOTAL OF THESE NEUTRALS EXITING MAIN P.  ',TEXIT)
c
      if (mtcopt.eq.1) then
        call prb
        call prc ('NEUTRAL MOMENTUM TRANSFER COLLISION SUMMARY:')
        call prb
        call prc ('ORIGINAL NEUTRALS:')
        call prr ('NUMBER OF ORIGINAL NEUTRALS               ',
     >                                 tneut)
        call prr ('NUMBER OF MTC EVENTS FOR ORIG NEUTRALS    ',
     >                                 mtcinf(1,1))
        if (mtcinf(1,1).gt.0.0.and.tneut.gt.0) then  

         call prr ('AVERAGE NUMBER OF MTC EVENTS/NEUTRAL     ',
     >                   mtcinf(1,1)/tneut)
         call prr ('AVERAGE TIME BETWEEN MTC EVENTS          ',
     >                    mtcinf(3,1)/mtcinf(1,1)*fsrate)
         call prr ('AVERAGE TIME TO FIRST MTC EVENT          ',
     >                    mtcinf(2,1)/mtcinf(1,1)*fsrate)
         call prr ('AVERAGE DISTANCE BETWEEN MTC EVENTS      ',
     >                    mtcinf(4,1)/mtcinf(1,1))
         call prr ('AVERAGE DISATNCE TO FIRST MTC EVENT      ',
     >                    mtcinf(5,1)/mtcinf(1,1))
c
        endif
c
        if (cfolrec.eq.1) then
c
          call prb
          call prc ('RECOMBINED NEUTRALS:')
          call prr ('NUMBER OF RECOMBINED NEUTRALS              ',
     >                                 recneut)
          call prr ('NUMBER OF MTC EVENTS FOR REC.NEUTRALS      ',
     >                                 mtcinf(1,2))
          if (mtcinf(1,2).gt.0.0) then

            call prr ('AVERAGE NUMBER OF MTC EVENTS/REC.NEUTRAL ',
     >                   mtcinf(1,2)/recneut)
            call prr ('AVERAGE TIME BETWEEN MTC EVENTS          ',
     >                    mtcinf(3,2)/mtcinf(1,2)*fsrate)
            call prr ('AVERAGE TIME TO FIRST MTC EVENT          ',
     >                    mtcinf(2,2)/mtcinf(1,2)*fsrate)
            call prr ('AVERAGE DISTANCE BETWEEN MTC EVENTS      ',
     >                    mtcinf(4,2)/mtcinf(1,2))
            call prr ('AVERAGE DISATNCE TO FIRST MTC EVENT      ',
     >                    mtcinf(5,2)/mtcinf(1,2))
          endif
        endif
c
        call prb
        call prc('NUMBER OF NEUTRALS HAVING "N" MTC EVENTS')
        call prc('                      ORIGINAL / RECOMB.')
c
        do in = 0,11
           if (in.ne.11) then
             write(comment,'(i5,1x,''MTC EVENTS: '',f12.3,2x,
     >                          f12.3)')
     >                    in,mtctotcnt(in,1),mtctotcnt(in,2)
           else
             write(comment,'(''  >10 MTC EVENTS: '',f12.3,
     >                       2x,f12.3)')
     >                    mtctotcnt(in,1),mtctotcnt(in,2)
           endif
c
           if (mtctotcnt(in,1).gt.0.0.or.mtctotcnt(in,2).gt.0.0)
     >           then
              call prc(comment)
           endif
c
        end do
c
      endif
c
      call prb
c      call prc ('ION SUMMARY:')
      call prchtml('ION SUMMARY:','pr_ion','0','B')
      CALL PRB
      CALL PRr0('TOTAL NUMBER OF IONS CREATED             ',TATIZ)
      CALL PRr0('TOTAL NO OF IONS PLATING ON WALL+TRAP    ',TWALL)
      CALL PRr0('TOTAL NO OF IONS PLATING ON WALL         ',
     >          TNTOTS(MAXIZS+1,1)-TNTOTS(0,1))
      CALL PRr0('TOTAL NO OF IONS PLATING ON TRAP         ',
     >          TNTOTS(MAXIZS+1,2)-TNTOTS(0,2))
      CALL PRr0('TOTAL NO OF IONS PLATING ON TARGET       ',TDEP)
C
      IF (FPOPT.EQ.3.or.fpopt.eq.5) THEN
c
        CALL PRC('FAR PERIPHERY OPTION: STATISTICS')
        CALL PRr0('  NUMBER ENTERING FP                     ',
     >                                  FPENT)
        CALL PRr0('  NUMBER EXITING FP TO PLASMA            ',
     >                                  FPEXIT)
        CALL PRr0('  NUMBER EXISTING TO TMAX IN FP          ',
     >                                  FPTTOT)
        CALL PRr0('  NUMBER LOST TO FP TARGET               ',
     >                                  FPTART)
C
      ENDIF
c
      if (cstop.eq.1) then
         call prr0('TOTAL IONS NOT FOLLOWED AFTER CORE ENTRY ',
     >              stopped_follow)
      endif
C
      CALL PRr0('TOTAL NO OF IONS IONISED BEYOND LIMIT    ',TBYOND)
      CALL PRr0('TOTAL IONS RECOMBINED TO FORM NEUTRALS   ',TBELOW)

c
      if (cfolrec.eq.1) then
c
c        Sum up over recombined neutral array to get results for
c        all recombined ions.
c
         do in = 1,14
            do id = 1,6
               recinf(in,7) = recinf(in,7) + recinf(in,id)
            end do
         end do
c
c        Normalize array values
c
         do in = 2,3
            do id = 1,7
               if (recinf(14,id).eq.0.0) then
                  recinf(in,id) = 0.0
               else
                  recinf(in,id) = recinf(in,id)/recinf(14,id)
               endif
            end do
         end do
c
         do in = 4,13
            do id = 1,7
               if (recinf(1,id).eq.0.0) then
                  recinf(in,id) = 0.0
               else
                  recinf(in,id) = recinf(in,id)/recinf(1,id)
               endif
            end do
         end do

c
         call prb
         call prc('SUMMARY OF RECOMBINED IMPURITY'//
     >                ' NEUTRAL RESULTS:')
         call prr0('  TOTAL NUMBER FOLLOWED            ',recneut)
         call prr0('  - TOTAL NUMBER ENDING ON TARGET  ',recstruk)
c
         if (mtcopt.eq.1) then
            CALL PRR0('   - NUMBER ENDING ON TARG AFTER MTC',
     >                                                MTCRECSTRUK)
         endif
c
         call prr0('  - TOTAL NUMBER ENDING ON WALLS   ',recwalln)
c
         if (mtcopt.eq.1) then
            call prr0('   - NUMBER ENDING ON WALL AFTER MTC',
     >                                                mtcrecwalln)
         endif
c
         call prr0('  - NUMBER EXCEEDING MAX REFLECTION',recloss)
c
         call prr0('  - NUMBER RE-IONIZED              ',recatiz)
         call prr0('  - NUMBER OF FAILED LAUNCHES      ',recfail)
         call prr0('  - NUMBER ENDING AT MAX TIME      ',rectmax)
         call prr0('  - NUMBER ENTERING MAIN           ',recmain)
         call prr0('  - NUMBER EXITING MAIN            ',recexit)
         call prr0('  - NUMBER STRIKING CENTRAL MIRROR ',reccent)
         if (tatiz.gt.0) then 
             call prr0('  AVERAGE NO. OF RECOMB/ORIG ION   ',
     >                          recneut/tatiz)
         endif

c
         call prb
c
         if (recinf(1,5).gt.0.0) then
            call prr('  TOTAL NUMBER OF REIONIZATIONS       ',
     >                  recinf(1,5))
            call prr('  AVERAGE TIME FOR RE-IONIZATION      ',
     >                       recinf(4,5)*qtim)
            call prr('  AVERAGE TIME TO FIRST RE-IONIZATION ',
     >                recinf(2,5)*qtim)
            call prr('  AVERAGE DISTANCE TO RE-IONIZATION   ',
     >                recinf(5,5))
            call prr('  AVERAGE DIST TO FIRST RE-IONIZATION ',
     >                recinf(3,5))
            call prr('  AVERAGE TEMP OF RE-IONIZED NEUTRALS ',
     >                recinf(6,5))
            call prr('  AVERAGE VEL OF RE-IONIZED NEUTRALS  ',
     >                recinf(7,5))
            call prr('  AVERAGE NE AT RECOM.POSITION        ',
     >                recinf(8,5))
            call prr('  AVERAGE NH AT RECOM.POSITION        ',
     >                recinf(9,5))
            call prr('  AVERAGE TE AT RECOM.POSITION        ',
     >                recinf(10,5))
            call prr('  AVERAGE R  AT RECOM.POSITION        ',
     >                recinf(11,5))
            call prr('  AVERAGE Z  AT RECOM.POSITION        ',
     >                recinf(12,5))
            call prr('  AVERAGE S  OR (SMAX-S) AT RECOM.    ',
     >                recinf(13,5))
            call prr('  AVERAGE VEL OF RE-IONIZED NEUTRALS  ',
     >                recinf(7,5))
            call prr('  RECOMBINATION VELOCITY MULTIPLIER   ',
     >                       cvrmult)
         endif
c
         call prb
         call prc('NUMBER OF IONS HAVING "N" RECOMBINATIONS')
c
         do in = 0,11
c
            if (in.ne.11) then
              write(comment,'(i5,'' RECOMBINATIONS: '',f12.3)')
     >                    in,rectotcnt(in)
            else
             write(comment,'(''  >10 RECOMBINATIONS: '',f12.3)')
     >                      rectotcnt(in)
            endif
c
            if (rectotcnt(in).gt.0.0) then
               call prc(comment)
            endif
c
         end do
c
      endif
c
c     Back to original
c
c
c     Print out summary for ions reflected as neutrals
c
      if (refflag.eq.1) then

         call prb
         call prc('SUMMARY OF REFLECTED (ION) IMPURITY'//
     >                ' NEUTRAL RESULTS:')
         call prr0('   TOTAL NUMBER FOLLOWED            ',refneut)
         call prr0('   - TOTAL NUMBER ENDING ON TARGET  ',refstruk)
c
         if (mtcopt.eq.1) then
            CALL PRR0('   - NUMBER ENDING ON TARG AFTER MTC',
     >                                                MTCREfSTRUK)
         endif
c
         call prr0('   - TOTAL NUMBER ENDING ON WALLS   ',refwalln)
c
         if (mtcopt.eq.1) then
            call prr0('   - NUMBER ENDING ON WALL AFTER MTC',
     >                                                mtcrefwalln)
         endif
c
         call prr0('   - NUMBER EXCEEDING MAX REFLECTION',refloss)
c
         call prr0('   - NUMBER RE-IONIZED              ',refatiz)
         call prr0('   - NUMBER OF FAILED LAUNCHES      ',reffail)
         call prr0('   - NUMBER ENDING AT MAX TIME      ',reftmax)
         call prr0('   - NUMBER ENTERING MAIN           ',refmain)
         call prr0('   - NUMBER EXITING MAIN            ',refexit)
         call prr0('   - NUMBER STRIKING CENTRAL MIRROR ',refcent)
         if (tatiz.gt.0) then 
            call prr0('   AVERAGE NO. OF RECOMB/ORIG ION   ',
     >                          recneut/tatiz)
         endif

         call prb
      endif
c
c
c
      call prr0('TOTAL IONS CROSSING SEPARATRIX       ',
     >                                           totleakcore)
      CALL PRr0('TOTAL NO OF IONS EXISTING AT TMAX    ',TCUT)
c
      call prb
c
c     Calculate and print TAUP(Core)
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            tmpncore = tmpncore + ncore (ik,ir)
            tmpiz = tmpiz + tizs (ik,ir,0)
         end do
      end do
c
      if (tmpiz.gt.0) then 
         call prr('TAUP CORE CALCULATED: ',tmpncore/tmpiz)        
      endif
c
      call prb
      call prchtml('SUMMARY OF INITIAL IMPURITY NEUTRAL IONIZATION',
     >             'pr_ioniz','0','B')
      call prc('FOR NEUTRALS ORIGINATING FROM SIDE 1 and SIDE 2')
      call prb
c
      if (refct.eq.1) then
         call prc('  NOTE: The (R,Z) values listed are for grid')
         call prc('  coordinates prior to the grid being reflected')
         call prc('  in the R-axis.')
      endif


      call pr_trace('DIV','PRINT ADDITIONAL SUMMARIES')
c
C-----------------------------------------------------------------------
c
c
c      The subroutine PRIONIZ takes care of printing the
c      summary of the ionization data requested. It checks
c      to see if the index combination is valid before
c      printing - so if reflected data are chosen
c      for printing and reflection is off the routine
c      prints nothing. The default values of irflct and
c      ifp of 1 each are always printed.
c
c      If values greater than 2 are requested - the printing
c      routine will generate the various SUMS of the independent
c      components - note that the S data will not be printed
c      OR calculated for values with a component in the main
c      plasma.
c
       do isol = 1,3
         do irflct = 1,2
           do ifp = 1,2
             call prioniz(isol,irflct,ifp,ionizdat)
           end do
         end do
       end do
c
      call prb
c
c
c     Print out the summary of average S position and field line
c     of ionization for all neutrals from all target elements.
c
c
c      call prchtml('ANALYSIS OF CORE LEAKAGE','pr_leakage','0','B')
c
c      call prleakage
c
      CALL PRR('TOTAL FLUX*YIELD (ALL SOURCES)           ',FYTOT)
      CALL PRR('TOTAL FLUX       (ALL SOURCES)           ',FTOT)
      CALL PRR('SELF-SPUTTERING ENHANCEMENT FACTOR       ',SSEF)
      CALL PRR('EFFECTIVE YIELD                          ',YEFF)
C
      DO 222 IR = 1, NRS
        IF (CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4) THEN
          IF (IR.GE.IRSEP) KBACDS(1,IR) = LO
          IF (IR.GE.IRSEP) KFORDS(NKS(IR),IR) = LO
        endif
        DO 222 IK = 1, NKS(IR)
          IF (IR.EQ.1.OR.IR.EQ.IRTRAP.or.ir.eq.irtrap2)
     >       KINDS(IK,IR) = LO
          IF (IR.EQ.IRWALL.or.ir.eq.irwall2) KOUTDS(IK,IR) = LO
  222 CONTINUE
C
C-----------------------------------------------------------------------
C     SCALING OF ARRAYS TO BIN SIZES, NUMBERS LAUNCHED, ETC
C-----------------------------------------------------------------------
C
c
      call set_normalization_factors(facta,factb,tatiz,tneut,
     >                               fsrate,qtim,nizs,cneuta,0)
C
      CLLL(0) = TEXIT
      CMMM(0) = 0.0
      CNNN(0) = TMAIN


c
c     BOUNDARY RING DATA:
c
c     DEAL WITH ANY DENSITY DATA IN THE BOUNDARY RINGS:
c     THESE RINGS ARE VIRTUAL AND DO NOT CONTAIN ANY
c     MEANINGFUL INFORMATION. IN ORDER TO AVOID
c     PROBLEMS IN SUMMARIES OR OTHER PLACES WHERE THIS
c     DATA MIGHT BE ACCIDENTALLY INCLUDED - IT IS SET
c     TO ZERO HERE.
c
c     The rings affected are IR=1, IR=IRWALL and IR=IRTRAP
c     for grid types 0 and 3 - JET and SONNET
c
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.
c slmod begin - ribbon grid
     >    cgridopt.eq.4.or.cgridopt.eq.5.or.
     >    cgridopt.eq.LINEAR_GRID.or.
     >    cgridopt.eq.RIBBON_GRID) then
c
c     >    cgridopt.eq.4.or.cgridopt.eq.5) then
c slmod end
c
         ir =1
         do ik = 1,nks(ir)
            do iz = -1,nizs
              ddlims(ik,ir,iz) = 0.0
            end do
         end do
c
         ir = irwall
         do ik = 1,nks(ir)
            do iz = -1,nizs
               ddlims(ik,ir,iz) = 0.0
            end do
         end do
c
         if (irtrap.lt.nrs) then
            ir = irtrap
            do ik = 1,nks(ir)
               do iz = -1,nizs
                  ddlims(ik,ir,iz) = 0.0
               end do
            end do
         endif
c
      endif

C
C---- DEAL WITH ANOMALY FOR CONTINUOUS RINGS, WHERE TWO POINTS ON THE
C---- RING ARE COINCIDENT - COMBINE FIRST & LAST POINTS ON RING.
C
      DO 4100 IZ = -1, NIZS
       ELIMS(1,3,IZ) = ELIMS(1,3,IZ) + ELIMS(NKS(IRSEP-1),3,IZ)
       ELIMS(NKS(IRSEP-1),3,IZ) = 0.0
       DO 4100 IR = 1, IRSEP-1
        DDTS  (1,IR,IZ) = DDTS  (1,IR,IZ) + DDTS  (NKS(IR),IR,IZ)
        DDLIMS(1,IR,IZ) = DDLIMS(1,IR,IZ) + DDLIMS(NKS(IR),IR,IZ)
        TIZS  (1,IR,IZ) = TIZS  (1,IR,IZ) + TIZS  (NKS(IR),IR,IZ)
        DDTS  (NKS(IR),IR,IZ) = 0.0D0
        DDLIMS(NKS(IR),IR,IZ) = 0.0D0
        TIZS  (NKS(IR),IR,IZ) = 0.0
c
        DO IT = 1, NTS
           LIMS(1,IR,IZ,IT)= LIMS(1,IR,IZ,IT)+LIMS(NKS(IR),IR,IZ,IT)
           LIMS(NKS(IR),IR,IZ,IT)= 0.0
        end do
c
c       (RIV)
c
        sdvs(1,ir,iz)  = sdvs(1,ir,iz) + sdvs(nks(ir),ir,iz)
        sdvs(nks(ir),ir,iz) = 0.0
c
        if (debugv) then
c
           sdvs2(1,ir,iz)  = sdvs2(1,ir,iz) + sdvs2(nks(ir),ir,iz)
           sdvs2(nks(ir),ir,iz) = 0.0
c
           sdvs3(1,ir,iz,1)=sdvs3(1,ir,iz,1)+sdvs3(nks(ir),ir,iz,1)
           sdvs3(nks(ir),ir,iz,1) = 0.0
           sdvs3(1,ir,iz,2)=sdvs3(1,ir,iz,2)+sdvs3(nks(ir),ir,iz,2)
           sdvs3(nks(ir),ir,iz,2) = 0.0
c
        endif

 4100 CONTINUE
c
C
C====================== Recorded Ion Velocity (RIV) ===================
C
C
      call check_ddlim(nizs,1)
c
      if (debugv) then

        do iz = 1,nizs

          DO 4210 IR = 1, NRS
            DO 4220 IK = 1, NKS(IR)
              IF (DDLIMS(IK,IR,IZ).GT.0.0D0) THEN
c
c               Also calculate the mean ion velocity.
c
                sdvs(IK,IR,IZ) = sdvs(IK,IR,IZ) / DDLIMS(IK,IR,IZ) /qtim
c
                sdvs2(IK,IR,IZ) = sdvs2(IK,IR,IZ) / DDLIMS(IK,IR,IZ)
     >                          / qtim**2
                sdtimp(ik,ir,iz)=(sdvs2(ik,ir,iz)-sdvs(ik,ir,iz)**2)
     >                         / 9.58084e7 * crmi
c
                sdvs3(ik,ir,iz,1)=sdvs3(ik,ir,iz,1)
     >                           /sdvs3(ik,ir,iz,2)/qtim
c
                sdvs3(ik,ir,iz,2)=sdvs3(ik,ir,iz,2)/ddlims(ik,ir,iz)

c
              ENDIF
 4220       CONTINUE
 4210     CONTINUE
        end do
c
c       Renormalize background velocity
c
        do ir = 1,nrs
           do ik = 1,nks(ir)
              sdvb(ik,ir) = sdvb(ik,ir)/qtim
           end do
        end do
c
c       Print out debug information on cells and their contents
c       where there are particles with Vz > Vb (local)
c
        write (6,*) 'SUMMARY OF IMPURITY ION VELOCITIES'//
     >              ' EXCEEDING LOCAL SOUND SPEED:'
c
        write (6,*) ' IK   IR   IZ  VB=(2kT/m)^0.5   AVERAGE VZ  '//
     >              '    FRACTION       NUMBER'
c
        do ir = 1,nrs
           do ik = 1,nks(ir)
              do iz = 1,nizs
c
                 if (sdvs3(ik,ir,iz,2).gt.0.0) then
                    write (6,'(3i5,1x,g13.5,1x,g13.5,1x,2g16.8)')
     >                    ik,ir,iz,sdvb(ik,ir),sdvs3(ik,ir,iz,1),
     >                    sdvs3(ik,ir,iz,2),sdvs3(ik,ir,iz,2)
     <                    *ddlims(ik,ir,iz)
c
                 endif
c
              end do
           end do
        end do
c
c
c       Calculate the distribution of velocities
c
        if (maxvnks.gt.nks(injir).and.nks(injir).gt.0) then
           ikv = nks(injir)
        else
           ikv = 1
        endif
c
        do ik = 1,ikv
c
           do 4228 iz = 1,nizs
c
              veltot(iz)  = 0.0
c
              do 4225 in = -nvel,nvel+1
                 if (velspace(in,iz,ik).ne.0.0) then
                    velspace(in,iz,ik) = velspace(in,iz,ik)
     >                              /velweight(in,iz,ik)/qtim
                 endif
c
                 veltot(iz) = veltot(iz) + velweight(in,iz,ik)
c
 4225         continue
c
              do 4226 in = -nvel,nvel+1
                 if (veltot(iz).ne.0.0) then
                    velweight(in,iz,ik) = velweight(in,iz,ik)
     >                       / veltot(iz)
                 endif
 4226         continue


 4228      continue
c
        end do
c
        do ik = 1,ikv
c
           if (ikv.eq.1) then
              write (6,*) 'Velocity Weight Distribution: ',
     >                                velplate/qtim
              vel = velplate/qtim
           else
              write (6,*) 'Velocity Weight Distribution: Knot=',ik,
     >                                velcell(ik)/qtim
              vel = velcell(ik)/qtim
           endif
c
           do 4227 in = -nvel,nvel+1
              if (in.le.0) then
                 velcoord(in) = (float(in)-0.5)*velsep
              elseif (in.gt.0) then
                velcoord(in) = (float(in-1)+0.5)*velsep
              endif
              write(6,'(i4,1x,f8.3,1x,f14.6,12(1x,g16.6))') in,
     >                   velcoord(in),velcoord(in)*vel,
     >                   (velweight(in,iz,ik),iz=1,nizs),
     >                   (velspace(in,iz,ik),iz=1,nizs)
 4227      continue

        end do

c
c       Print out velocity array
c
        WRITE(6,'(//1X,''TABLE OF VELOCITY VALUES:'')')
        DO 4230 IR = 1, NRS
           WRITE (6,9031) 'PRIMARY','TOTNEUT','IONIZ 1',
     >                 'IONIZ 2','IONIZ 3','IONIZ 4',
     >                 'IONIZ 5','IONIZ 6'
           WRITE (6,9032)
           DO 4235 IK = 1, NKS(IR)
              WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),
     >            (sdvs(ik,ir,iz),IZ=-1,NIZS)
 4235      CONTINUE
 4230   CONTINUE

c
c       Print out impurity ion temperature array
c
        WRITE(6,'(//1X,''TABLE OF ION TEMPERATURE VALUES:'')')
        DO 4236 IR = 1, NRS
           WRITE (6,9031) 'PRIMARY','TOTNEUT','IONIZ 1',
     >                 'IONIZ 2','IONIZ 3','IONIZ 4',
     >                 'IONIZ 5','IONIZ 6'
           WRITE (6,9032)
           DO 4238 IK = 1, NKS(IR)
             WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),
     >            (sdtimp(ik,ir,iz),IZ=-1,NIZS)
 4238      CONTINUE
 4236   CONTINUE
c
c       End of debugv
c
      endif

 9031 FORMAT(/1X,'  IK  IR    R      Z  ',12(2X,A7))
 9032 FORMAT(1X,131('-'))
 9033 FORMAT(1X,2I4,2F7.3,1P,12E9.2)
 9034 FORMAT(39X , 1P , 12E9.2 )


      call pr_trace('DIV','FORCES')

C
C====================== FORCES ===================================
C
c
C---- DEAL WITH ANOMALY FOR CONTINUOUS RINGS, WHERE TWO POINTS ON THE
C---- RING ARE COINCIDENT - COMBINE FIRST & LAST POINTS ON RING.
      DO IZ = 1, NIZS
       DO IR = 1, IRSEP-1
        Fcell(1,IR,IZ) = Fcell(1,IR,IZ)  + Fcell(NKS(IR),IR,IZ)
        Ffi(1,IR,IZ)   = Ffi(1,IR,IZ)    + Ffi(NKS(IR),IR,IZ)
        Fthi(1,IR,IZ)  = Fthi(1,IR,IZ)   + Fthi(NKS(IR),IR,IZ)
        Fvbg(1,IR,IZ)  = Fvbg(1,IR,IZ)   + Fvbg(NKS(IR),IR,IZ)
        DIFF(1,IR,IZ)  = DIFF(1,IR,IZ)   + DIFF(NKS(IR),IR,IZ)
        Velavg(1,IR,IZ)= Velavg(1,IR,IZ) + Velavg(NKS(IR),IR,IZ)
        Fcell(NKS(IR),IR,IZ)  = 0.0e0
        Ffi(NKS(IR),IR,IZ)    = 0.0e0
        Fthi(NKS(IR),IR,IZ)   = 0.0e0
        Fvbg(NKS(IR),IR,IZ)   = 0.0e0
        DIFF(NKS(IR),IR,IZ)   = 0.0e0
        Velavg(NKS(IR),IR,IZ) = 0.0e0
       END DO
      END DO
c---------------------------------------------------------------------c
c
c    -Calculate the average force per grid cell per charge state
c    -Calculate the average ion velocity / grid cell / charge state
c    -Output force values to the screen or to a data file
c
c---------------------------------------------------------------------c
c
      DO iz = 1,nizs
         DO ir = 1,nrs
            DO ik = 1,nks(ir)
               IF ( DDLIMS(ik,ir,iz).NE.0.0 ) THEN

                  Fcell(ik,ir,iz) = Fcell(ik,ir,iz)/DDLIMS(ik,ir,iz)
                  Ffi(ik,ir,iz)   = Ffi(ik,ir,iz)/DDLIMS(ik,ir,iz)
                  Fthi(ik,ir,iz)  = Fthi(ik,ir,iz)/DDLIMS(ik,ir,iz)
                  Fvbg(ik,ir,iz)  = Fvbg(ik,ir,iz)/DDLIMS(ik,ir,iz)
                  DIFF(ik,ir,iz)  = DIFF(ik,ir,iz)/DDLIMS(ik,ir,iz)
                  Velavg(ik,ir,iz)= Velavg(ik,ir,iz)/DDLIMS(ik,ir,iz)
     >                              /QTIM

               ENDIF
            END DO
         END DO
      END DO

C
C====================== TEMPERATURES ===================================
C
c     Include neutrals
c
      DO 4290 IZ = -1, NIZS
        DSUM1 = 0.0D0
        DSUM3 = 0.0D0
        DSUM4=  0.0D0
        DO 4250 IR = 1, NRS
          DO 4240 IK = 1, NKS(IR)
            DSUM1 = DSUM1 + DDTS  (IK,IR,IZ)
            DSUM3 = DSUM3 + DDLIMS(IK,IR,IZ)
            DSUM4 = DSUM4 + ddlims(ik,ir,iz) * sdtimp(ik,ir,iz)
 4240     CONTINUE
 4250   CONTINUE
        IF (DSUM3.GT.0.0D0) SDTZS(IZ) = DSUM1 / DSUM3
        IF (DSUM3.GT.0.0D0) SDTZS2(IZ) = DSUM4 / DSUM3
C
        DO 4270 IR = 1, NRS
          DO 4260 IK = 1, NKS(IR)
            IF (DDLIMS(IK,IR,IZ).GT.0.0D0) THEN
c
              DDTS(IK,IR,IZ) = DDTS(IK,IR,IZ) / DDLIMS(IK,IR,IZ)
c
            ENDIF
 4260     CONTINUE
 4270   CONTINUE
C
 4290 CONTINUE
C
C================= IONISATION & TIME-DEPENDENT =========================
C
      DO 4430 IZ = -1, NIZS
       DO 4420 IR = 1, NRS
        DO 4410 IK = 1, NKS(IR)
         if (kareas(ik,ir).ne.0.0) then
            FACT = FACTA(IZ) / KAREAS(IK,IR)
         else
            FACT = 0.0
         endif
c
         TIZS(IK,IR,IZ) = FACT * TIZS(IK,IR,IZ)
c
c        Normalize ionization of chemiclaly sputtered component.
c
         if (iz.eq.0) then
            chemizs(ik,ir) = fact *chemizs(ik,ir)
         endif
c
c         DO 4400 IT = 1, NTS
c          LIMS(IK,IR,IZ,IT) = FACT * LIMS(IK,IR,IZ,IT)
c 4400    CONTINUE
c
 4410   CONTINUE
 4420  CONTINUE
 4430 CONTINUE
C
C===================== NUMBER DENSITY ==================================
C
      DO 4630 IZ = -1, NIZS
       DO 4620 IR = 1, NRS
        DO 4610 IK = 1, NKS(IR)
          if (kareas(ik,ir).ne.0.0) then
             DACT = DBLE (FACTB(IZ) / KAREAS(IK,IR))
          else
             DACT = 0.0
          endif
c
          DDLIMS(IK,IR,IZ) = DACT * DDLIMS(IK,IR,IZ)
c
c         Normalize the time dependent portion too
c
          DO IT = 1, NTS
             LIMS(IK,IR,IZ,IT) = DACT * LIMS(IK,IR,IZ,IT)
     >             /(ctimes(it,iz)-ctimes(it-1,iz))
          END DO
c
c         Normalize the neutral density due to chemical sputtering
c
          if (iz.eq.0) then
             chemden(ik,ir) = dact * chemden(ik,ir)
          endif
c
 4610   CONTINUE
 4620  CONTINUE
 4630 CONTINUE
c
      call check_ddlim(nizs,2)
c
c     Normalize data on the subgrid if it is in use.
c
      call norm_subgrid(nizs,cneuta,tneut,tatiz,fsrate,qtim)
c
c     Void region - number density - no areas involved.
c
      do in = 1,3
         ddvoid(in) = ddvoid(in) * factb(0)
      end do
C
C================= DEPOSITION, NET EROSION AND WALLS ===================
C
      IF (NIZS.GT.0) THEN
        DO 881 IZ = 1, NIZS
          DO 880 ID = 1, NDS
            IF (DDS(ID).NE.0.0) THEN
              DEPS(ID,IZ) = DEPS(ID,IZ) / DDS(ID) * FACTA(0)
            ELSE
              DEPS(ID,IZ) = 0.0
            ENDIF
  880     CONTINUE
  881   CONTINUE
      ENDIF
C
      FACT = 0.0
      IF (TDEP.GT.0.0) FACT = TNEUT / TDEP
      DO 883 ID = 1, NDS
        IF (DDS(ID).NE.0.0) THEN
          NEROS(ID,1) =-NEROS(ID,1) / DDS(ID) * FACTA(0)
          NEROS(ID,2) = NEROS(ID,2) / DDS(ID) * FACTA(0)
c          NEROS(ID,2) = NEROS(ID,2) / DDS(ID) * FACTA(-1)
          NEROS(ID,3) = NEROS(ID,3) / DDS(ID) * FACTA(0)
        ELSE
          NEROS(ID,1) = 0.0
          NEROS(ID,2) = 0.0
          NEROS(ID,3) = 0.0
        ENDIF
        NEROS(ID,4) = NEROS(ID,1) + NEROS(ID,3)
        NEROS(ID,5) = FACT * NEROS(ID,1) + NEROS(ID,3)
  883 CONTINUE
C
C       ONLY NORMALIZE THE INTEGRATED DATA FOR NOW
C
        DO 889 IR = 1, NRS
          DO 888 IK = 1, NKS(IR)
            if (kareas(ik,ir).ne.0.0) then
               WALLS(IK,IR,MAXIZS+1) = WALLS(IK,IR,MAXIZS+1)
     >                             / KAREAS(IK,IR) * FACTA(0)
            else
c               WALLS(IK,IR,MAXIZS+1) = WALLS(IK,IR,MAXIZS+1)
c     >                             / KAREAS(IK,IR) * FACTA(0)
               WALLS(IK,IR,MAXIZS+1) = 0.0
            endif
  888     CONTINUE
  889   CONTINUE
C
C============== RADIATIVE POWER LOSS AND LINE RADIATION ================
C
      CALL RZERO (POWLS, MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RZERO (LINES, MAXNKS*MAXNRS*(MAXIZS+2))
      CALL RZERO (HPOWLS, MAXNKS*MAXNRS*2)
      CALL RZERO (HLINES, MAXNKS*MAXNRS*2)
c
      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) then
         CALL RZERO (e2dPOWLS, MAXNKS*MAXNRS*MAXe2dizs)
         CALL RZERO (e2dLINES, MAXNKS*MAXNRS*maxe2dizs)
      endif

      call pr_trace('DIV','RADIATED POWER CALC')

c
c     Choose Nocorona or ADAS
c
c     Nocorona
c
      if (cdatopt.eq.0) then
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
      DO 990 IR = 1, NRS
C
C---- SET TEMPS AND DENSITIES FOR EACH X BIN.  CONVERT TO CM**3
C---- SET BIN VOLUMES IN CM**3  (ASSUME 1M IN THIRD DIMENSION)
C---- COPY DENSITIES (CM**-3) FOR NOCORONA
C
        DO 900 IK = 1, NKS(IR)
          PTES(IK) = KTEBS(IK,IR)
          PNES(IK) = KNBS(IK,IR) * 1.E-6 * RIZB
 900   CONTINUE
        DO 920 IK = 1, NKS(IR)
          PDVOLS(IK) = 1.0E6 * KAREAS(IK,IR)
          DO 910 IZ = 0, NIZS
            PNZS(IZ+1,1,IK) = 1.0E-6 * SNGL(DDLIMS(IK,IR,IZ))
 910     CONTINUE
 920   CONTINUE
C
C------ CALCULATE POWER LOSS (W.CM**-3) AND LINE RAD (W.CM**-3)
C------ RECONVERT TO W.M**-3
C------ VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY
C------ GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.
C
c        CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NKS(IR))
c
        call ncrdlong (nks(ir))
c
        DO 940 IK = 1, NKS(IR)
          DO 930 IZ = 0, NIZS
            POWLS(IK,IR,IZ) = MAX (0.0, PRADIS(1,IZ+1,1,IK)*1.0E6)
            LINES(IK,IR,IZ) = MAX (0.0, PRADIS(3,IZ+1,1,IK)*1.0E6)
 930     CONTINUE
 940   CONTINUE
C
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1) LOCATIONS
C
        DO 960 IK = 1, NKS(IR)
          IF (DDLIMS(IK,IR,0).LE.0.0D0) THEN
            POWLS(IK,IR,-1) = 0.0
            LINES(IK,IR,-1) = 0.0
          ELSE
            POWLS(IK,IR,-1) = POWLS(IK,IR,0) *
     >                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
            LINES(IK,IR,-1) = LINES(IK,IR,0) *
     >                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
          ENDIF
  960   CONTINUE
  990 CONTINUE


c
c     Calculate E2DLINES and E2DPOWLS if EDGE2D data including
c     impurity data have been loaded for reference.
c

      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) then

C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
      DO IR = 1, NRS
C
C---- SET TEMPS AND DENSITIES FOR EACH X BIN.  CONVERT TO CM**3
C---- SET BIN VOLUMES IN CM**3  (ASSUME 1M IN THIRD DIMENSION)
C---- COPY DENSITIES (CM**-3) FOR NOCORONA
C
        DO  IK = 1, NKS(IR)
          PTES(IK) = KTEBS(IK,IR)
          PNES(IK) = KNBS(IK,IR) * 1.E-6 * RIZB
        end do
        DO IK = 1, NKS(IR)
          PDVOLS(IK) = 1.0E6 * KAREAS(IK,IR)
          DO IZ = 0, cre2dIZS
             PNZS(IZ+1,1,IK) = 1.0E-6 * e2dnzs(ik,ir,iz)
          end do
        end do
C
C------ CALCULATE POWER LOSS (W.CM**-3) AND LINE RAD (W.CM**-3)
C------ RECONVERT TO W.M**-3
C------ VALUES OF -1 ARE RETURNED WHERE THE TEMPERATURE OR DENSITY
C------ GOES OUTSIDE THE ALLOWABLE RANGE - THESE ARE CHECKED FOR BELOW.
C
c        CALL RDLONG (PTES,PNES,PNZS,PDVOLS,PRADIS,NKS(IR))
c
        call ncrdlong (nks(ir))
c
        DO IK = 1, NKS(IR)
          DO IZ = 0, cre2dIZS
            e2dPOWLS(IK,IR,IZ) = MAX (0.0, PRADIS(1,IZ+1,1,IK)*1.0E6)
            e2dLINES(IK,IR,IZ) = MAX (0.0, PRADIS(3,IZ+1,1,IK)*1.0E6)
          end do
        end do
c
      end do
c
c     Corresponds to if (cre2d.eq.1...
c
      end if
c
c
c     ADAS - use ADAS data when ADPAK is not available
c
      elseif (cdatopt.eq.1.or.cdatopt.eq.2.or.cdatopt.eq.3) then
c
C
      DO 1190 IR = 1, NRS
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
        DO 1100 IK = 1, NKS(IR)
          PTESA(IK) = KTEBS(IK,IR)
          PNESA(IK) = KNBS(IK,IR) * RIZB
          PNBS(IK) = KNBS(IK,IR)
          PNHS(IK) = PINATOM(IK,IR)
 1100   CONTINUE
        DO 1120 IK = 1, NKS(IR)
          DO 1110 IZ = 0, NIZS
            PNZSA(IK,IZ) = SNGL(DDLIMS(IK,IR,IZ))
 1110     CONTINUE
 1120   CONTINUE
C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C
        write(year,'(i2.2)') iyearz
        call xxuid(useridz)
c        YEAR = '89'
c        YEARDF = '89'
        ICLASS = 5
        MIZS = MIN(CION-1,NIZS)
        DO 1130 IZ = 0, MIZS
          CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
     +                PCOEF(1,IZ+1))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ+1))
          DO 1135 IK = 1, NKS(IR)
            LINES(IK,IR,IZ) = PCOEF(IK,IZ+1)*PNESA(IK)*PNZSA(IK,IZ)
            POWLS(IK,IR,IZ) = LINES(IK,IR,IZ)
c            write (6,'(a,3i5,3g16.8)') 'Debug DIV:',ir,ik,iz,
c     >              pcoef(ik,iz+1),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      DIV:',
c     >            lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
 1135     CONTINUE
 1130   CONTINUE
        ICLASS = 4
        MIZS = MIN(CION,NIZS)
        DO 1140 IZ = 1, MIZS
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
     +                PCOEF(1,IZ))
c          CALL ADASRD(YEAR,YEARDF,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
c     +                PCOEF(1,IZ))
          DO 1145 IK = 1, NKS(IR)
            POWLS(IK,IR,IZ) = POWLS(IK,IR,IZ)
     +                        + PCOEF(IK,IZ)*PNESA(IK)*PNZSA(IK,IZ)
 1145     CONTINUE
 1140   CONTINUE
C
C------ DEAL WITH PRIMARY NEUTRALS STORED IN DDLIMS(,,-1)
C
        DO 1160 IK = 1, NKS(IR)
          IF (DDLIMS(IK,IR,0).LE.0.0) THEN
            POWLS(IK,IR,-1) = 0.0
            LINES(IK,IR,-1) = 0.0
          ELSE
            POWLS(IK,IR,-1) = POWLS(IK,IR,0) *
     +                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
            LINES(IK,IR,-1) = LINES(IK,IR,0) *
     +                        SNGL(DDLIMS(IK,IR,-1) / DDLIMS(IK,IR,0))
c
c            write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
c     >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      POW:',
c     >           lines(ik,ir,iz), powls(ik,ir,iz),ddlims(ik,ir,iz)
c
          ENDIF
 1160   CONTINUE
C
C------ REPEAT FOR HYDROGEN POWER
C
        write(year,'(i2.2)') iyearh
        call xxuid(useridh)
c        YEAR = '89'
c        YEARDF = '89'
        ICLASS = 5
        CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
c        CALL ADASRD(YEAR,YEARDF,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
        DO 1136 IK = 1, NKS(IR)
          HLINES(IK,IR,0) = PCOEF(IK,1)*PNESA(IK)*PNHS(IK)
          HPOWLS(IK,IR,0) = HLINES(IK,IR,0)
 1136   CONTINUE
        ICLASS = 4
        CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
c        CALL ADASRD(YEAR,YEARDF,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
        DO 1146 IK = 1, NKS(IR)
          HPOWLS(IK,IR,1) = PCOEF(IK,1)*PNESA(IK)*PNBS(IK)
 1146   CONTINUE
 1190 CONTINUE

c
c     Calculate E2DLINES and E2DPOWLS if EDGE2D data including
c     impurity data have been loaded for reference.
c

      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.0) then
c
c
c
      DO IR = 1, NRS
C
C---- LOAD POWER DATA ONE RING AT A TIME.
C
        DO IK = 1, NKS(IR)
          PTESA(IK) = KTEBS(IK,IR)
          PNESA(IK) = KNBS(IK,IR) * RIZB
          PNBS(IK) = KNBS(IK,IR)
c
          if (cpinopt.eq.1.or.cpinopt.eq.4) then
             PNHS(IK) = PINATOM(IK,IR)
          else
             PNHS(IK) = E2DATOM(IK,IR)
          endif
c
        end do
        DO IK = 1, NKS(IR)
          DO IZ = 0, cre2dIZS
             PNZSA(IK,IZ) = e2dnzs(ik,ir,iz)
          end do
        end do

C
C------ GET POWER LOSS FROM ADAS DATA FILES. LOAD TOTAL LINE RADIATION
C------ INTO LINES AND ADD RECOMBINATION AND BREMSSTRAHLUNG POWER TO
C------ GET TOTAL RADIATIVE LOSSES
C
        write(year,'(i2.2)') iyearz
        call xxuid(useridz)
c
        ICLASS = 5
        MIZS = MIN(CION-1,cre2dizs)
c
        DO IZ = 0, MIZS
          CALL ADASRD(YEAR,CION,IZ+1,ICLASS,NKS(IR),PTESA,PNESA,
     +                PCOEF(1,IZ+1))
          DO  IK = 1, NKS(IR)
            e2dLINES(IK,IR,IZ) = PCOEF(IK,IZ+1)*PNESA(IK)*PNZSA(IK,IZ)
            e2dPOWLS(IK,IR,IZ) = e2dLINES(IK,IR,IZ)
c
c            write (6,'(a,3i5,3g16.8)') 'Debug E2d:',ir,ik,iz,
c     >              pcoef(ik,iz+1),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      E2D:',
c     >          e2dlines(ik,ir,iz),e2dpowls(ik,ir,iz),e2dnzs(ik,ir,iz)
c
          end do
        end do
c
        ICLASS = 4
        MIZS = MIN(CION,cre2dizs)
c
        DO IZ = 1, MIZS
c
          CALL ADASRD(YEAR,CION,IZ,ICLASS,NKS(IR),PTESA,PNESA,
     +                PCOEF(1,IZ))
          DO IK = 1, NKS(IR)
            e2dPOWLS(IK,IR,IZ) = e2dPOWLS(IK,IR,IZ)
     +                        + PCOEF(IK,IZ)*PNESA(IK)*PNZSA(IK,IZ)
c
c            write (6,'(a,3i5,3g16.8)') 'Debug POW:',ir,ik,iz,
c     >              pcoef(ik,iz),pnesa(ik),pnzsa(ik,iz)
c            write (6,'(a,15x,3g16.8)') '      POW:',
c     >          e2dlines(ik,ir,iz),e2dpowls(ik,ir,iz),e2dnzs(ik,ir,iz)
c
          end do
        end do
c
      end do
c
c     End of EDGE2D Power calculations
c
      end if
c
c     End of selection for nocorona/ADAS
c
      endif
c

      call pr_trace('DIV','GLOBAL TOTALS')
C
C================ GLOBAL TOTALS OVER ENTIRE PLASMA =====================
C
c
c     Analyse the Radiation source. Determine mean ne and niz
c     from the radiating volume. (Assumed to be highest 2/3 of
c     radiating cells.
c
      call radproc(nizs,rions,nimps)
c
      CALL DZERO (DTOTS, 46)
      call dzero (ditots,(maxizs+2)*2)
      call dzero (ptots,10*(maxizs+3))
c
      call rzero (e2dtots,9)
      call rzero (e2dptots,6*(maxe2dizs+1))
c
      DO 4050 IR = 1, NRS
        DO 4040 IK = 1, NKS(IR)
          DACT = KAREAS(IK,IR)
c
c         Be sure to EXCLUDE the repeated cell that closes
c         the core rings - otherwise - it's area will be
c         counted twice.
c
          IF (IR.LT.IRSEP.and.ik.ne.nks(ir)) THEN
            DTOTS(1) = DTOTS(1) + DACT
            DTOTS(2) = DTOTS(2) + DACT * DBLE(KNBS(IK,IR))
            DTOTS(33) = DTOTS(33) + DACT * DBLE(HPOWLS(IK,IR,0))
     >                            + DACT * DBLE(HPOWLS(IK,IR,1))
            DTOTS(34) = DTOTS(34) + DACT * DBLE(HLINES(IK,IR,0))
     >                            + DACT * DBLE(HLINES(IK,IR,1))
            if ((kss(ik,ir).ge.0.25*kss(nks(ir)-1,ir)).and.
     >          (kss(ik,ir).le.0.75*kss(nks(ir)-1,ir))) then
                dtots(36) = dtots(36) + dact
            endif
c
          ELSE
            DTOTS(15) = DTOTS(15) + DACT
            DTOTS(16) = DTOTS(16) + DACT * DBLE(KNBS(IK,IR))
            DTOTS(31) = DTOTS(31) + DACT * DBLE(HPOWLS(IK,IR,0))
     >                            + DACT * DBLE(HPOWLS(IK,IR,1))
            DTOTS(32) = DTOTS(32) + DACT * DBLE(HLINES(IK,IR,0))
     >                            + DACT * DBLE(HLINES(IK,IR,1))
          ENDIF
c
c         EXCLUDE repeated core point.
c
          if (ir.eq.irsep-1.and.ik.ne.nks(ir)) then
             ditots(maxizs+2,1) = ditots(maxizs+2,1) + dact
          endif
c
          if ((ir.lt.irsep).and.
     >        (kss(ik,ir).ge.0.25*kss(nks(ir)-1,ir)).and.
     >        (kss(ik,ir).le.0.75*kss(nks(ir)-1,ir))) then
              dtots(36) = dtots(36) + dact
c
              if (ir.eq.irsep-1)
     >             ditots(maxizs+2,2)=ditots(maxizs+2,2)+dact
c
          endif
c
c         Add up area in the divertor region.
c
          if (cgridopt.eq.0) then
            IF (ZS(IK,IR).GT.ZXP) THEN
              dtots(44)= dtots(44) + dact
            ENDIF
          elseif (cgridopt.eq.1.or.
     >            cgridopt.eq.2.or.cgridopt.eq.3) then
            IF (ZS(IK,IR).LT.ZXP) THEN
              dtots(44)= dtots(44) + dact
            ENDIF
          endif
c
c         Sum over ionization states
c
          DO 4030 IZ = 1, NIZS
            IF (IR.LT.IRSEP.and.ik.ne.nks(ir)) THEN
              DTOTS(3) = DTOTS(3) + DACT * DDLIMS(IK,IR,IZ)
c
              if ((kss(ik,ir).ge.0.25*kss(nks(ir)-1,ir)).and.
     >            (kss(ik,ir).le.0.75*kss(nks(ir)-1,ir))) then
                  dtots(35) = dtots(35) + dact * ddlims(ik,ir,iz)
                  dtots(36) = dtots(36) + dact
              endif
c
              DTOTS(4) = DTOTS(4) + DACT * DBLE(POWLS(IK,IR,IZ))
              DTOTS(5) = DTOTS(5) + DACT * DBLE(LINES(IK,IR,IZ))
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1.and.
     >             iz.le.cre2dizs) then
                 e2dtots(1) = e2dtots(1) + dact * e2dnzs(ik,ir,iz)
                 e2dtots(2) = e2dtots(2) + dact * e2dpowls(ik,ir,iz)
                 e2dtots(3) = e2dtots(3) + dact * e2dlines(ik,ir,iz)
              endif
c
              if (ir.eq.irsep-1) then
                ditots(iz,1) = ditots(iz,1) + dact * ddlims(ik,ir,iz)
                ditots(maxizs+1,1) = ditots(maxizs+1,1) + dact *
     >                             ddlims(ik,ir,iz)

                if ((kss(ik,ir).ge.0.25*kss(nks(ir)-1,ir)).and.
     >             (kss(ik,ir).le.0.75*kss(nks(ir)-1,ir))) then
                    ditots(iz,2) = ditots(iz,2)
     >                               + dact * ddlims(ik,ir,iz)
                    ditots(maxizs+1,2) = ditots(maxizs+1,2) + dact *
     >                             ddlims(ik,ir,iz)
                endif
              endif
c
            ELSEif (ir.ge.irsep.and.ir.le.irwall) then
c
              DTOTS(6) = DTOTS(6) + DACT * DDLIMS(IK,IR,IZ)
              DTOTS(7) = DTOTS(7) + DACT * DBLE(POWLS(IK,IR,IZ))
              DTOTS(8) = DTOTS(8) + DACT * DBLE(LINES(IK,IR,IZ))
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1.and.
     >             iz.le.cre2dizs) then
                 e2dtots(4) = e2dtots(4) + dact * e2dnzs(ik,ir,iz)
                 e2dtots(5) = e2dtots(5) + dact * e2dpowls(ik,ir,iz)
                 e2dtots(6) = e2dtots(6) + dact * e2dlines(ik,ir,iz)
              endif
c
            ELSEif (ir.ge.irwall.and.ir.le.nrs) then
c
              DTOTS(37) = DTOTS(37) + DACT * DDLIMS(IK,IR,IZ)
              DTOTS(38) = DTOTS(38) + DACT * DBLE(POWLS(IK,IR,IZ))
              DTOTS(39) = DTOTS(39) + DACT * DBLE(LINES(IK,IR,IZ))
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1.and.
     >             iz.le.cre2dizs) then
                 e2dtots(7) = e2dtots(7) + dact * e2dnzs(ik,ir,iz)
                 e2dtots(8) = e2dtots(8) + dact * e2dpowls(ik,ir,iz)
                 e2dtots(9) = e2dtots(9) + dact * e2dlines(ik,ir,iz)
              endif
            ENDIF
c
            if (xpoint_up) then
c
              IF (ZS(IK,IR).GT.ZXP) THEN
                DTOTS(21)= DTOTS(21)+ DACT * DDLIMS(IK,IR,IZ)
              ENDIF
              IF (ZS(IK,IR).GT.CZD) THEN
                DTOTS(22)= DTOTS(22)+ DACT * DDLIMS(IK,IR,IZ)
              ENDIF
c
            elseif (.not.xpoint_up) then
c
              IF (ZS(IK,IR).LT.ZXP) THEN
                DTOTS(21)= DTOTS(21)+ DACT * DDLIMS(IK,IR,IZ)
              ENDIF
              IF (ZS(IK,IR).LT.CZD) THEN
                DTOTS(22)= DTOTS(22)+ DACT * DDLIMS(IK,IR,IZ)
              ENDIF
            endif
c
 4030     CONTINUE
C
C         FIND NEUTRAL CONTRIBUTION TO TOTAL POWER AND LINE RAD
C         Also find contributions for each ionization state
C
          do 4035 iz = 0,nizs
c
c
c           ptots(1,iz) - power in core
c                 2     - line rad in core
c                 3     - power in main SOL
c                 4     - line rad in main SOL
c                 5     - power in PP
c                 6     - line rad in PP
c
c                 7     - power in main SOL below XP (DIVERTOR)
c                 8     - line rad          below        "
c                 9     - power in main SOL above XP
c                10     - line rad          above XP
c
c               maxizs+1 = sum of ion power and line rad
c               maxizs+2 = total ion and neutral power and line rad
c
            if (ir.lt.irsep.and.ik.ne.nks(ir)) then
c
              ptots(1,iz) = ptots(1,iz) + dact * dble(powls(ik,ir,iz))
c
              ptots(2,iz) = ptots(2,iz) + dact * dble(lines(ik,ir,iz))
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.iz.le.cre2dizs) then
                 e2dptots(1,iz)=e2dptots(1,iz)+dact*e2dpowls(ik,ir,iz)
                 e2dptots(2,iz)=e2dptots(2,iz)+dact*e2dlines(ik,ir,iz)
              endif
c
            elseif (ir.ge.irsep.and.ir.le.irwall) then
c
              ptots(3,iz) = ptots(3,iz) + dact * dble(powls(ik,ir,iz))

              ptots(4,iz) = ptots(4,iz) + dact * dble(lines(ik,ir,iz))
c
              if (zs(ik,ir).gt.zxp) then
                 if (xpoint_up) then
                    ptots(7,iz) = ptots(7,iz) +
     >                                dact * dble(powls(ik,ir,iz))
                    ptots(8,iz) = ptots(8,iz) +
     >                                dact * dble(lines(ik,ir,iz))
                 else
                    ptots(9,iz) = ptots(9,iz) +
     >                                dact * dble(powls(ik,ir,iz))
                    ptots(10,iz) = ptots(10,iz) +
     >                                dact * dble(lines(ik,ir,iz))
                 endif
              else
                 if (xpoint_up) then
                    ptots(9,iz) = ptots(9,iz) +
     >                                dact * dble(powls(ik,ir,iz))
                    ptots(10,iz) = ptots(10,iz) +
     >                                dact * dble(lines(ik,ir,iz))
                 else
                    ptots(7,iz) = ptots(7,iz) +
     >                                dact * dble(powls(ik,ir,iz))
                    ptots(8,iz) = ptots(8,iz) +
     >                                dact * dble(lines(ik,ir,iz))
                 endif
              endif
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.iz.le.cre2dizs) then
                 e2dptots(3,iz)=e2dptots(3,iz)+dact*e2dpowls(ik,ir,iz)
                 e2dptots(4,iz)=e2dptots(4,iz)+dact*e2dlines(ik,ir,iz)
              endif
c
            elseif (ir.ge.irtrap.and.ir.le.nrs) then
c
              ptots(5,iz) = ptots(5,iz) + dact * dble(powls(ik,ir,iz))

              ptots(6,iz) = ptots(6,iz) + dact * dble(lines(ik,ir,iz))
c
              if ((cre2d.eq.1.or.cre2d.eq.2).and.iz.le.cre2dizs) then
                 e2dptots(5,iz)=e2dptots(5,iz)+dact*e2dpowls(ik,ir,iz)
                 e2dptots(6,iz)=e2dptots(6,iz)+dact*e2dlines(ik,ir,iz)
              endif
c
            endif
4035      continue
c
4040    CONTINUE
4050  CONTINUE
c
c     Sum up over ptots arrays -
c
      if (tatiz.ne.0) then
         dact = dble(tneut)/dble(tatiz)
      else
         dact = 1.0
      endif
c
      do in = 1,10
         do iz = 1,nizs
            ptots(in,maxizs+1) = ptots(in,maxizs+1) + ptots(in,iz)
         end do
c
         ptots(in,maxizs+2) = ptots(in,maxizs+1) + ptots(in,0) * dact
c
      end do
c
c     Reassign to dtots array for compatibility
c     Assigns neutral contributions to special positions
c
c     neutral power main
c
      DTOTS(23)= ptots(1,0)
c
c     neutral line main
c
      DTOTS(24)= ptots(2,0)
c
c     neutral power sol
c
      DTOTS(25)= ptots(3,0)
c
c     neutral power trap
c
      DTOTS(40)= ptots(5,0)
c
c     neutral line sol
c
      DTOTS(26)= ptots(4,0)
c
c     neutral line trap
c
      DTOTS(41)= ptots(6,0)
c
c
      if (tneut.ne.0.and.(tatiz.le.tneut)) then
c
c
c       Code changes because the absolute factor is calculated for
c       one ION entering the system/second - the previous formula
c       would have been valid for an absolute factor based on
c       one NEUTRAL entering/second.
c
c        dact = dble(tatiz)/dble(tneut)
c
        if (tatiz.ne.0) then
           dact = dble(tneut)/dble(tatiz)
        else
           dact = 1.0
        endif
c
c       total power main
c
c        dtots(27) = dtots(23) + dact * dtots(4)
c
        dtots(27) = dact * dtots(23) + dtots(4)
c
c
c       total line main
c
c        dtots(28) = dtots(24) + dact * dtots(5)
c
        dtots(28) = dact * dtots(24) +  dtots(5)
c
c
c       total power sol+trap
c
c        dtots(29) = dtots(25) + dtots(40) + dact * dtots(7)
c     >                        + dact * dtots(38)
c
        dtots(29) = dact*dtots(25) + dact*dtots(40) + dtots(7)
     >                        + dtots(38)
c
c       total line sol+trap
c
c        dtots(30) = dtots(26) + dtots(41) + dact * dtots(8)
c     >                        + dact * dtots(39)
c
        dtots(30) = dact*dtots(26) + dact*dtots(41) + dtots(8)
     >                        + dact * dtots(39)
c
      else
        dtots(27) = dtots(4)
        dtots(28) = dtots(5)
        dtots(29) = dtots(7) + dtots(38)
        dtots(30) = dtots(8) + dtots(39)
      endif


      call pr_trace('DIV','ZEFFS')

c
C
C======================== Z EFFECTIVES, ETC ============================
C
C
C
C    CALCULATE THE ABSOLUTE SCALING FACTOR. THIS VALUE WILL BE ZERO
C    FOR THOSE CASES WHERE THE ABSFAC IS A MEANINGLESS QUANTITY.
C    E.G. CASES CONSISTING SOLELY OF WALL LAUNCHES.
C    IT WILL BE MODIFIED BY A "WALL SOURCE STRENGTH FACTOR" FOR
C    THOSE CASES WHERE A COMBINATION OF PLATE AND WALL LAUNCHES IS
C    MODELLED. IT WILL BE ZERO FOR OTHER SUPPLEMENTARY LAUNCHES.
C
C    DAVID ELDER  FEB 5 / 1992
C
C    NOTE: IF IT IS DECIDED TO SET ABSFAC = 1.0 FOR THE OUT OF
C          CONDITION CASES - SIMPLY MOVE THE DEFINITION
C          OF ABSFAC TO THE POINTS WHERE THE DIFFERENT WSSF'S
C          ARE DEFINED.
C
C    NOTE2: THE FACTOR WSSF IS DEFINED AND CALCULATED BEFORE ANY
C           OF THE INPUT OPTIONS ARE MODIFIED DURING THE RUNNING
C           OF THE CODE. E.G. CNEUTB SET TO ZERO FOR SELF-SPUTTER
C           RE - LAUNCHES.
C
C    NOTE3: ABSFAC MODIFIED TO BE EQUAL TO 1.0 FOR INJECTION CASES
C           THIS ALLOWS TABLES OF DENSITIES TO BE OBTAINED IN OUT
c
c    NOTE4: ABSFAC has historically been scaled to 1 ion/sec entering
c           the system. However, this resulted in a bug in the treatment
c           of neutral particle densities - these densities were scaled by
c           1/TNEUT - the total number of neutrals entering the system.
c           However, ABSFAC represented the total number of ions entering the
c           system - thus the neutral densities were incorrect
c           by a factor of TATIZ/TNEUT. This could be corrected by scaling
c           the neutral data to one ion entering the system/sec or by scaling
c           all of the results to 1 particle/s where a "particle" is a neutral
c           for initial neutral launch cases and and ion for ion injection
c           cases. This means that absfac needs to be modified by removing the
c           TATIZ/TNEUT factor. In other words absfac becomes the same as absfac_neut
C
c
c     jdemod - I think Karl's code is appropriate for any case with a specified
c              absolute factor and self-sputtering assuming the absolute
c              scaling factor specified is not intended to scale the entire
c              source including self sputtering.
c
C     IPP/08 Krieger - added code for proper handling of self sputtering
c     cases with homogeneous wall launch and predefined ABSFAC
c     This was added by K. Schmid but I do not know if it is still
c     required in the new version
c      if(nabsfac.gt.0.0.and.cselfs.gt.0.and.cneutbSAV.eq.2) then
c         IF (TNEUT.GT.0.0.and.(tatiz.le.tneut)) then
c            ABSFAC = WSSF*(nabsfac + SSEF * nabsfac)*CSEF
c     >                          + neut2d_fytot
c            absfac_neut = absfac
c         endif
c         if (cneuta.eq.0.and.tneut.gt.0.0) then
c            absfac_ion = absfac * tatiz/tneut
c         else
c            absfac_ion = absfac
c         endif
c
c
c     jdemod - comment out the special case for SOL 21 - it is a bug waiting
c              happen when running cases
c
c      if (nabsfac.gt.0.0.and.cioptf.eq.21) then
c         absfac = qrat * nabsfac /(dtots(27)+dtots(29)) + neut2d_fytot
c         absfac_neut = absfac
c         if (cneuta.eq.0.and.tneut.gt.0.0) then
c            absfac_ion = absfac * tatiz/tneut
c         else
c            absfac_ion = absfac
c         endif
c      elseif (nabsfac.gt.0.0.and.cioptf.ne.21) then
c
c     SSEF is only non-zero if self-sputtering occurs - so the following code
c     can be generalized by just including SSEF.
c
      IF (sloutput) write(0,*) 'nabsfac=',nabsfac

      if (nabsfac.gt.0.0) then
         absfac = (nabsfac + SSEF*nabsfac) + neut2d_fytot
         absfac_neut = absfac
         if (cneuta.eq.0.and.tneut.gt.0.0) then
            absfac_ion = absfac * tatiz/tneut
         else
            absfac_ion = absfac
         endif
         write(6,'(a,10(1x,g12.5))') 'NABSFAC CALCULATION:',nabsfac,
     >            ssef,
     >            neut2d_fytot,absfac_ion,absfac_neut
c slmod begin
        if (nabsfac.eq.1.0) then
          write(6,*) 
          write(6,*) '-------------------------------------------------'
          write(6,*) ' WARNING: ABSFACV_neut & _ion forced = 1.0'
          write(6,*) '-------------------------------------------------'
          write(6,*) 

          write(0,*) 
          write(0,*) '-------------------------------------------------'
          write(0,*) ' WARNING: ABSFACV_neut & _ion forced = 1.0'
          write(0,*) '-------------------------------------------------'
          write(0,*) 

          absfac      = nabsfac
          absfac_ion  = absfac
          absfac_neut = absfac
        endif
c slmod end
      elseif (lpinz0) then
         absfac = zioniz*(tatiz/nimps) + neut2d_fytot
         absfac_neut = absfac
         if (cneuta.eq.0.and.tneut.gt.0.0) then
            absfac_ion = absfac * tatiz/tneut
         else
            absfac_ion = absfac
         endif
c
         DO IR = 1,NRS
           DO IK = 1,NKS(IR)
             DDLIMS(IK,IR,0) = DDLIMS(IK,IR,0)/(tatiz/nimps)
           ENDDO
         ENDDO
c
      else
c
c     Calculate absolute scaling factor
c
c     jdemod - hard code the first normalization option
c            - old code can be removed later
c
c
c
c         if (.true.) then
c
            ABSFAC = 1.0
            absfac_neut = absfac
            absfac_ion = absfac
c
c           Note: the calculated absolute factor only makes senfse for 
c                 neutral launch cases - the various quantities are 
c                 not defined for ion injection cases - so the ABSFAC 
c                 should then be set to 1.0
c
            if (cneuta.eq.0.and.tneut.gt.0.0) then
               ABSFAC = WSSF*FTOT*YEFF*CSEF + neut2d_fytot
               absfac_neut = absfac
               absfac_ion = absfac * tatiz/tneut
            endif
            write(6,'(a,10(1x,g12.5))') 'ABSFAC CALCULATION:',
     >               wssf,ftot,yeff,csef,
     >            neut2d_fytot,absfac_ion,absfac_neut
c
c         else
c
c            ABSFAC = 1.0
c            absfac_neut = absfac
c            IF (TNEUT.GT.0.0.and.(tatiz.le.tneut)) then
c               ABSFAC = WSSF*FTOT*YEFF*CSEF*TATIZ/TNEUT
c     >                          + neut2d_fytot
c               absfac_neut = tneut/tatiz * absfac
c               absfac_ion = absfac
c            endif
c         endif
c
      endif

c
c     jdemod - absolute factor calculated - write to ERO particle track file if needed
c
      if (ero_part_output_opt.eq.1) then 
         call close_ero_part_output(absfac,tneut)
      endif

c
      WRITE(6,*) 'ABSFAC: ',ABSFAC,' PARTS: WSSF ',WSSF,' FTOT ',
     >   FTOT,' YEFF ',YEFF,' CSEF ',CSEF,' TATIZ ',TATIZ,' TNEUT ',
     >   TNEUT,'NBASFAC:',nabsfac,'NEUT2D_FYTOT:',neut2d_fytot
C slmod begin - t-dep
c...  Dump the particle distribution to a file:
      IF (opt_div%pstate.EQ.1) THEN
        CALL find_free_unit_number(fp)
        OPEN (UNIT=fp,FILE='raw.divimp_tdep_dist',ACCESS='SEQUENTIAL',
     .        STATUS='REPLACE')
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '{VERSION}'
        WRITE(fp,*    ) 1.0
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '{DATA RUN}'
        WRITE(fp,'(1P,E15.7,0P,A)') qtim  ,'  qtim'
        WRITE(fp,'(1P,E15.7,0P,A)') cstmax,'  cstmax'
        WRITE(fp,'(1P,E15.7,0P,A)') absfac,'  absfac'
        WRITE(fp,'(A)') '*'
        WRITE(fp,'(A)') '{DATA PARTICLE}'
        WRITE(fp,*    ) tdep_save_n
        WRITE(fp,'(A)') '*'
        DO i = 1, tdep_save_n
          WRITE(fp,'(I9,7F13.7,1P,E15.7,0P,F5.1,1P,E15.7,0P)') i,
     .      tdep_save(i)%r      ,
     .      tdep_save(i)%z      ,	  
     .      tdep_save(i)%phi    ,     
     .      tdep_save(i)%s      ,  
     .      tdep_save(i)%cross  ,
     .      tdep_save(i)%diag   ,
     .      tdep_save(i)%temp   ,
     .      tdep_save(i)%vel    ,
     .      tdep_save(i)%charge ,
     .      tdep_save(i)%weight
        ENDDO
        CLOSE(fp)
c...    Clear memory:
        DEALLOCATE(tdep_save)
        IF (ALLOCATED(tdep_load)) DEALLOCATE(tdep_load)
      ENDIF
C slmod end
C
      DO 1020 IR = 1, NRS
        DO 1010 IK = 1, NKS(IR)
C
C-------- SUM CLOUD DENSITIES OVER IONISATION STATES
C-------- CALCULATE TOTAL IMPURITY INFLUX, MULTIPLY BY SPUTTERING
C-------- ENHANCEMENT FACTOR...
C
          DSUM1 = 0.0D0
          DSUM2 = 0.0D0
          DO 1000 IZ = 1, NIZS
            DIZ   = DBLE(IZ)
            DSUM1 = DSUM1 +  DIZ  * DDLIMS(IK,IR,IZ)
            DSUM2 = DSUM2 +DIZ*DIZ* DDLIMS(IK,IR,IZ)
1000      CONTINUE
          DSUM1 = DSUM1 * DBLE(ABSFAC)
          DSUM2 = DSUM2 * DBLE(ABSFAC)
C
C-------- CALCULATE NIE,   ZB.NBTRUE,   ZEFF
C
          ZEFFS(IK,IR,1) = SNGL (DSUM1)
          ZEFFS(IK,IR,2) = RIZB * KNBS(IK,IR) - ZEFFS(IK,IR,1)
          IF (ZEFFS(IK,IR,2).GT.0.0) THEN
            ZEFFS(IK,IR,3) = (RIZB * ZEFFS(IK,IR,2) + SNGL(DSUM2)) /
     >                       (RIZB * KNBS(IK,IR))
          ELSE
            ZEFFS(IK,IR,3) = SNGL (DSUM2/DSUM1)
          ENDIF
1010    CONTINUE
1020  CONTINUE
c
C=======================================================================
c
c     IMPURITY CONTENT
c
      call rzero (impurity_content,(maxizs+1)*4*2)
c

c
      do ir = 1,nrs
        do ik = 1,nks(ir)
c
          IF (ABS(RS(IK,IR)-RXP).LT.0.4.and.
     >       ((cgridopt.eq.0.and.ZS(IK,IR).GT.ZXP).or.
     >       ((cgridopt.eq.1.or.cgridopt.eq.2.or.cgridopt.eq.3).and.
     >       ZS(IK,IR).LT.ZXP))) then
            DTOTS(9) = DTOTS(9) + DACT
            DTOTS(10)= DTOTS(10)+ DACT * DBLE(ZEFFS(IK,IR,1))
            DTOTS(11)= DTOTS(11)+ DACT * DBLE(ZEFFS(IK,IR,2))
            DTOTS(12)= DTOTS(12)+ DACT * DBLE(ZEFFS(IK,IR,3))
          ENDIF
c
c         Calculate the impurity content of various regions
c         on the grid. (Include impurity content reported by
c         supplemental code - PINZ0 - can be used by PIN, EIRENE or
c         may be loaded with a UEDGE result.)
c
c         X-point up
c
          if (xpoint_up) then
c
c            Divertor
c
             if (zs(ik,ir).ge.zxp) then
c
                dtots(43) = dtots(43) + ddlims(ik,ir,0)
     >                                  * kareas(ik,ir)
c
                dtots(46) = dtots(46) + pinz0(ik,ir)
     >                                  * kareas(ik,ir)
c
c               Outer
c
                if (ik.lt.nks(ir)/2) then
c
                   do iz = 0,nizs
c
                      impurity_content(iz,2,1)=impurity_content(iz,2,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,2,2)=
     >                           impurity_content(iz,2,2)
     >                         + e2dnzs(ik,ir,iz) * kareas(ik,ir)

                      endif
                   end do
c
c               Inner
c
                else
c
                   do iz = 0,nizs
c
                      impurity_content(iz,1,1)=impurity_content(iz,1,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c

                      if (iz.le.cre2dizs) then
                         impurity_content(iz,1,2)=
     >                            impurity_content(iz,1,2)
     >                          + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
                   end do
c
                endif
c
c            Main
c
             else
c
                dtots(42) = dtots(42) + ddlims(ik,ir,0)
     >                                  * kareas(ik,ir)
c
                dtots(45) = dtots(45) + pinz0(ik,ir)
     >                                  * kareas(ik,ir)
c
                if (ir.lt.irsep) then
c
                   do iz = 0,nizs
c
                      impurity_content(iz,4,1)=
     >                       impurity_content(iz,4,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,4,2)=
     >                       impurity_content(iz,4,2)
     >                     + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
                   end do
c
                else
c
                   do iz = 0,nizs
c
                      impurity_content(iz,3,1)=
     >                       impurity_content(iz,3,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,3,2)=
     >                       impurity_content(iz,3,2)
     >                     + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
                   end do
c
                endif

c
             endif
c
          elseif (.not.xpoint_up) then
c
c            Divertor
c
             if (zs(ik,ir).le.zxp) then
c
                dtots(43) = dtots(43) + ddlims(ik,ir,0)
     >                                  * kareas(ik,ir)
c
                dtots(46) = dtots(46) + pinz0(ik,ir)
     >                                  * kareas(ik,ir)
c
c               Inner
c
                if (ik.lt.nks(ir)/2) then
c
                   do iz = 0,nizs
c
                      impurity_content(iz,1,1)=impurity_content(iz,1,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,1,2)=
     >                          impurity_content(iz,1,2)
     >                        + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
                   end do
c
c               Outer
c
                else
c
                   do iz = 0,nizs
c
                      impurity_content(iz,2,1)=impurity_content(iz,2,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,2,2)=
     >                          impurity_content(iz,2,2)
     >                        + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
                   end do
c
                endif
c
c            Main
c
             else
c
                dtots(42) = dtots(42) + ddlims(ik,ir,0)
     >                                  * kareas(ik,ir)
c
                dtots(45) = dtots(45) + pinz0(ik,ir)
     >                                  * kareas(ik,ir)
c
c               Sum up for each charge state
c
                if (ir.lt.irsep) then
c
                   do iz = 0,nizs
c
                      impurity_content(iz,4,1)=
     >                       impurity_content(iz,4,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,4,2)=
     >                       impurity_content(iz,4,2)
     >                     + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
c
                   end do
c
                else
c
                   do iz = 0,nizs
c
                      impurity_content(iz,3,1)=
     >                       impurity_content(iz,3,1)
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
                      if (iz.le.cre2dizs) then
                         impurity_content(iz,3,2)=
     >                       impurity_content(iz,3,2)
     >                     + e2dnzs(ik,ir,iz) * kareas(ik,ir)
                      endif
c
                   end do
c
                endif
c
             endif
c
          endif
c
        end do
      end do
c
c
      call pr_trace('DIV','DTOTS')


c
      IF (DTOTS(9).NE.0.0D0) THEN
        DTOTS(10) = DTOTS(10) / DTOTS(9)
        DTOTS(11) = DTOTS(11) / DTOTS(9)
        DTOTS(12) = DTOTS(12) / DTOTS(9)
      ENDIF
      IF (DTOTS(1).NE.0.0D0) THEN
        DTOTS(13) = DTOTS(2)  / DTOTS(1)
        DTOTS(17) = DTOTS(3)  / DTOTS(1)
      ENDIF
      IF (DTOTS(13)*DTOTS(3).NE.0.0D0)
     >  DTOTS(14) = DTOTS(4)  / (DTOTS(13)*DTOTS(3))
c
      DO 4090 J = 1, 46
        STOTS(J) = SNGL(DTOTS(J))
4090  CONTINUE
c
      STOTS(18) = KPMAXS(IRSEP-1)
c
      do 4095 J = 1,10
        do 4095 iz = 0,maxizs+2
          sptots(j,iz) = sngl(ptots(j,iz))
4095  continue
c
      do in = 1,2
        do 4096 iz = 1,nizs
           if (ditots(maxizs+2,in).ne.0.0) then 
              sitots(iz,in) = ditots(iz,in) / ditots(maxizs+2,in)
           endif  

4096    continue
c
        if (ditots(maxizs+2,in).ne.0.0) then 
          sitots(maxizs+1,in) = ditots(maxizs+1,in)/ditots(maxizs+2,in)
        endif
c
      end do
c
C-----------------------------------------------------------------------
c
c     Print louver deposition information
c
C-----------------------------------------------------------------------
c
      call prb
c
      tmpsrc = wallsn(maxpts+1) + wallsi(maxpts+1)
      tmpmult = absfac * 2.0 * PI * R0
c
      call prb
      call prchtml('SUMMARY OF PARTICLE DEPOSITION/EROSION',
     >             'pr_depero','0','B')
      call prb
c
c     Louvers for JET grids only
c
      if (cgridopt.eq.0) then
         call prc('- LOUVER AREA IS ASSUMED TO BE THE OUTERMOST'//
     >         ' TWO ELEMENTS OF')
         call prc('  EACH TARGET')
      endif
c
      call prr('TOTAL EROSION         : ',wallse(maxpts+1))
      call prr('TOTAL IONIZED         : ',wallse_i(maxpts+1))
      call prr('TOTAL DEPOSITION      : ',tmpsrc)
      call prr('- NEUTRALS            : ',wallsn(maxpts+1))
      call prr('- IONS                : ',wallsi(maxpts+1))
c
      if (cgridopt.eq.0) then

         tmpdep = wallsn(wallindex(2)) + wallsn(wallindex(3))
c
         call prr('NEUTRAL DEPOSITION ON '//inner//' "LOUVER"  : ',
     >             tmpdep)

         if (tmpsrc.ne.0.0) then 
            call prr('APPROXIMATE '//inner//' C DEPOSITION RATE   : ',
     >             tmpdep/tmpsrc*tmpmult)
         endif
c
         tmpdep = wallsn(wallindex(nds-1)) + wallsn(wallindex(nds-2))
c
         call prr('NEUTRAL DEPOSITION ON '//outer//' "LOUVER"  : ',
     >             tmpdep)

         if (tmpsrc.ne.0.0) then 
           call prr('APPROXIMATE '//outer//' C DEPOSITION RATE   : ',
     >             tmpdep/tmpsrc*tmpmult)
         endif
c
      endif
c
      call prb
c
      call pr_trace('DIV','DEPS')


c     Calculate and summarize region deposition probabilities.
c
c     Based on contents of wallse wallse_i and wtdep array.
c
c     Calculate totals over 4 regions:
c
c     - Outer target
c     - Inner target
c     - Main vessel Wall
c     - PFZ Wall
c
c     with particle sources from each of the 4 regions - this gives 16 nubmers
c     in the end showing the probability that a particle beginning in one of
c     the four regions will deposit in one of the other four regions.
c
c
c     Sum total deposition (ion+neutral) into third array element
c
      do iw = 1,wallpts
c
         do in = 1,wallpts
c
            wtdep(iw,in,3) = wtdep(iw,in,1) + wtdep(iw,in,2)
c
         end do
c
      end do

      
c
c     Print and calculate wall deposition
c
      call pr_calc_walldep
c
      call prb

      call pr_trace('FINISHED CALC_WALLDEP')
C
C-----------------------------------------------------------------------
C                     PRINT REISER INFORMATION TO .LIM FILE
C                     FOR POSSIBLE PLOTS.
C-----------------------------------------------------------------------
C
c psmod
c
c     Write imp. number density, Velavg, Fcell, Ffi, Fthi, and Fvbg from
c     SOL region to .lim file to process 3D plots via Excel.
c
      if (cioptr.gt.0.and.(cprint.eq.8.or.cprint.eq.9)) then
         WRITE(6,*)'Writing impurity force data to .lim file'
         CALL DATA3DII(1)
      endif
c
c psmod
c
c
c     Calculate the wall distribution of any hydrogenic or impurity
c     radiation in POWLS and HPOWLS.
c
c     jdemod - this is quite computationally intensive so only calculate
c              it if the full printing options are on
c
      if (cprint.eq.9.or.cprint.eq.10) then 
         call calc_wallprad(nizs)
      endif

      call pr_trace('FINISHED CALC_WALLPRAD')

C
C-----------------------------------------------------------------------
C                     PRINT CLOSING MESSAGES
C-----------------------------------------------------------------------
C
      call pr_trace('DIV','CLOSING MESSAGES')


      IF (NIZS.GT.0)
     >     CALL MONPRI (FACTA(1),VFLUID,NIZS,SDTZS ,sdtzs2,
     >           STOTS,DOUTS,RIONS,TDEP,TWALL,DPARAS,DCROSS,
     >           TNTOTS,FPTARG,acttarg,sptots,sitots,coreouts,
     >           e2dtots,e2dptots,cre2d,cre2dizs,impurity_content)
c
      if (nizs.gt.0) call calcnt(nizs)
c
c     Print out line profile information
c
      call pr_line_profile
c
c     Print out any extra DIVIMP data to be saved from this run
c
      call wrtdivaux(nizs)
c
      CALL PRB
      CALL PRI ('NUMBER OF NEUTRALS FOLLOWED   ',NINT(TNEUT))
      CALL PRI ('NUMBER OF IONS FOLLOWED       ',NINT(TATIZ))
      WRITE (6,'('' NUMBER OF NEUTRALS FOLLOWED   '',G11.4)') TNEUT
      WRITE (6,'('' NUMBER OF IONS FOLLOWED       '',G11.4)') TATIZ
c
c
c     Print the pinchvel debug information to channel 5
c
      write(6,*) 'PINCHVEL DEBUGGING INFORMATION:'
      do in = -max_d_pinch_v,max_d_pinch_v
         write(6,'(a,f12.5,1x,1p,g18.10)') 'VEL:',
     >       in * d_pinch_vel, d_pinch_v(in)
      end do
c slmod begin
      CALL inOpenInterface('idl.divimp_summary',ITF_WRITE)
      i = nimps
      IF (cneuth.NE.-1) i = i + nimps2  ! Check for a supplementary launch
      CALL inPutData(i               ,'IONS_REQUESTED'       ,'N/A')
      CALL inPutData(tneut           ,'NEUTRALS_LAUNCHED'    ,'N/A')
      CALL inPutData(tfail           ,'NEUTRALS_FAILED'      ,'N/A')
      CALL inPutData(tatiz           ,'IONS_CREATED'         ,'N/A')
      CALL inPutData(num_entered_core,'IONS_REACHING_CORE'   ,'N/A')
      CALL inPutData(twall           ,'IONS_LOST_WALL'       ,'N/A')
      CALL inPutData(tdep            ,'IONS_LOST_TARGET'     ,'N/A')
      CALL inPutData(tbyond          ,'IONS_LOST_STATE_LIMIT','N/A')
      CALL inPutData(tbelow          ,'IONS_LOST_RECOMBINED' ,'N/A')
      CALL inPutData(cion            ,'ION_ATOMIC_NUMBER'    ,'N/A')
      CALL inPutData(nizs            ,'MAX_CHARGE_STATE'     ,'N/A')
      CALL inCloseInterface

      CALL OutputData(87,'END OF DIV')
c slmod end

      call pr_trace('DIV','END OF DIV')

c
c      if (cisterrcnt.ne.0) then
c         call prb
c         call prc('WARNING: DIVIMP TIMING ERRORS ENCOUNTERED.')
c         call prc('         NUMBER OF PARTICLE TRACKS EXCEEDING')
c         call pri('         MAXIMUM TIME RESOLUTION = ',cisterrcnt)
c         call prc('         MAXIMUM TIME STEPS =  16777216.00')
c         call prr('         MAXIMUM TIME_TRACK = ',16777216.00*qtim)
c      endif
c
c     In case of terminated procedure ... STOP the program at
c     this point.
c
      if (procterm) stop
C
C---- FORMATS ...
C
c nonorth
 9003 FORMAT(1X,I8,1x,F12.1,1x,2(1x,I4),1x,I2,2(1x,F12.5),1x,
     >       F9.3,1x,F6.2,1x,F12.5,1x,F8.3,1P,1x,E15.8,
     >       0P,1x,1x,F9.3,1P,1x,E10.3,0P,1x,F10.5,1x,F6.2,
     >       1xI3,:,1X,A,:,1x,F8.5)
c nonorth
c 9003 FORMAT(1X,I5,F9.1,2I3,I2,2F9.5,F8.3,F6.2,F8.3,1P,E15.8,
c     >  0P,F7.1,1P,E8.1,0P,F8.5,F5.2,I2,:,1X,A,:,F8.5)
 9004 FORMAT(//1X,'DIV DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,
     >  ' TIMESTEPS  (DELTA T =',1P,G10.3,' SECONDS).',//)
c
c nonorth
 9005 FORMAT(1X,' ------ION- ----TIME-  -IK- -IR- IZ  ',
     >  '---------R- ----------Z- -------S- ----K-  -----THETA- ',
     >  '---SMAX- ------DRIFTVEL-  ',
     >  '----TEMI- -PARADIFF- ----CROSS- -FRAC- -IT ',
     >  14('-'))
c 9005 FORMAT(1X,'--ION-----TIME-IK-IR-IZ',
c     >  '----R--------Z-------S------K---THETA--SMAX---',
c     >  '---DRIFTVEL------TEMI-PARADIFF-CROSS--FRAC-IT',
c     >  14('-'))
c nonorth
c 9005 FORMAT(1X,'--ION-----TIME-IK-IR-IZ',
c     >  '----R--------Z-------S------K----SMAX---',
c     >  '---DRIFTVEL------TEMI-PARADIFF-CROSS--FRAC-IT',
c     >  14('-'))
 9010 FORMAT(/1X,'ION',I5,' STARTING CONDITIONS:',1P,
     >  /5X,'   TEB',E9.2,',   NB ',E9.2,',    E1',A9  ,',   VB1',A9  ,
     >  /5X,'TAUIZ1',E9.2,',TAUEI1',E9.2,',CPROB1',E9.2,',RFRAC1',E9.2,
     >  /5X,'TAUCX1',E9.2,',   NH ',E9.2,', TAUP1',A9  ,', TAUS1',A9  ,
     >  /5X,' TAUT1',A9  ,',    K ',E9.2,',TAUIZ2',E9.2,',   TIB',E9.2,
     >  /5X,'LLLFPS',E9.2,',  FEG1',A9  ,',  FIG1',A9)
 9012 FORMAT(/1X,A,A,A)
 9013 FORMAT(/1X,A,I7,A)
 9022 FORMAT(1X,'DIV: ZENTRY',F6.3,', ZCREAT',F6.3,', ZTRIPP',F6.3,
     >  ', ZTRIPS',F6.3,', %P',F7.1,', %S',F7.1,'  (ION',I5,
     >  '  WEIGHT',F5.2,')',' TIME:',f10.2)
      RETURN
      END
c
c
c
      subroutine calcnt (nizs)
      implicit none
      integer nizs
      include 'params'
      include 'comtor'
      include 'cgeom'
      include 'dynam1'
c
c     Calculate total content in each ionization state in the trapped
c     region - this is used for SFT comparisons and calculation of
c     leakage. (Calculated only on the injection ring using the value
c     specified for CSTGRAD times SMAX for the injection ring)
c
c     NOTE: CSTGRAD is also used to specify the region of the ring
c           where temperature gradient forces -> 0 - this may or
c           may not be compatible with its use in this function to
c           define the impurity integration region - care should
c           be taken when looking at these results when using the options
c           which need cstgrad to be specified to turn off the temperature
c           gradient forces.
c
c
      real    nt(maxnrs,maxizs),cp,lastcp,nextcp,lastmid,nextmid
      real    scut,scut2
      integer ik,ir,iz
      character*77 comment
c
c     If the integration point has been specified equal to or less than
c     zero then this option is turned off - since the integral over the
c     entire plot range is printed out in th eOUT program.
c
      if (cstgrad.le.0.0) return
c
      call rzero(nt,maxizs*maxnrs)
c
c     Add up content in rings over each ionization state
c
      do iz = 1,nizs
         do ir = irsep,irwall
c
c           Start at scut and stop at scut2
c
            scut = cstgrad * ksmaxs(ir)
            scut2= ksmaxs(ir) - scut
c            write (6,*) 'scuts:',scut,scut2
c
            do ik = 1,nks(ir)
c
c              centrepoint
c
               cp   = kss(ik,ir)
c
c              last centre and cell boundary
c
               if (ik.ne.1) then
                  lastcp = kss(ik-1,ir)
                  lastmid= (lastcp+cp) / 2.0
               else
                  lastcp = 0.0
                  lastmid= 0.0
               endif
c
c              Next center and cell boundary
c
               if (ik.ne.nks(ir)) then
                  nextcp = kss(ik+1,ir)
                  nextmid= (nextcp+cp) / 2.0
               else
                  nextcp = ksmaxs(ir)
                  nextmid= ksmaxs(ir)
               endif
c
c              Add contribution if required
c
               if (lastmid.le.scut.and.nextmid.ge.scut) then
                  nt(ir,iz) = nt(ir,iz) + kareas(ik,ir) *
     >                        (nextmid-scut)/(nextmid-lastmid) *
     >                        ddlims(ik,ir,iz)
c
c             write (6,*) 'Start add:',ik,ir,iz,scut,lastmid,nextmid
c
c
               elseif (lastmid.le.scut2.and.nextmid.ge.scut2) then
c
                  nt(ir,iz) = nt(ir,iz) + kareas(ik,ir) *
     >                        (scut2-lastmid)/(nextmid-lastmid) *
     >                        ddlims(ik,ir,iz)
c
c             write (6,*) 'End add:',ik,ir,iz,scut2,lastmid,nextmid
c
c
               elseif (lastmid.ge.scut.and.nextmid.le.scut2) then
c
                  nt(ir,iz) = nt(ir,iz) + kareas(ik,ir) *
     >                        ddlims(ik,ir,iz)
c
               endif
            end do
         end do
      end do
c
c     Print various information depending on the initial quantities
c     specified.
c
      if ((cprint.eq.1.or.cprint.eq.9).and.(cneuta.eq.1).and.
     >   (ciopte.eq.2.or.ciopte.eq.3.or.ciopte.eq.5.or.ciopte.eq.6)
     >    .and.(cdperp.eq.0.0).and.(cizsc.eq.nizs)) then
         call prb
         call prr (' Impurity content is          : ',
     >              nt(injir,nizs))
         call pri2(' For ring and ionization state: ',injir,nizs)
         call prr2(' Integrated from / to (m)      : ',
     >             cstgrad * ksmaxs(injir), ksmaxs(injir)*(1.0-cstgrad))
         call prb
c
c
      elseif ((cprint.eq.1.or.cprint.eq.9).and.(cneuta.eq.1).and.
     >    (ciopte.eq.2.or.ciopte.eq.3.or.ciopte.eq.5.or.ciopte.eq.6)
     >    .and.(cdperp.eq.0.0).and.(cizsc.lt.nizs)) then
c
         call prb
         call pri(' Table of impurity content of ring :',injir)
         call prr2(' Integrated from / to (m)           : ',
     >             cstgrad * ksmaxs(injir), ksmaxs(injir)*(1.0-cstgrad))
         do iz = cizsc,nizs
            write(comment,4700) iz,nt(injir,iz)
            call prc(comment)
         end do
         call prb
c
c
      elseif ((cprint.eq.1.or.cprint.eq.9).and.(cneuta.eq.1).and.
     >    (ciopte.eq.2.or.ciopte.eq.3.or.ciopte.eq.5.or.ciopte.eq.6)
     >    .and.(cdperp.ne.0.0)) then
c
         call prb
         call prc(' Table of impurity content of SOL rings')
         write (comment,4710) (iz,iz=cizsc,nizs)
         call prc (comment)
         do ir = irsep,irwall
            call pri(' Summary for ring : ',ir)
            call prr2(' Integrated from / to (m)      : ',
     >             cstgrad * ksmaxs(ir), ksmaxs(ir)*(1.0-cstgrad))
            if ((nizs-cizsc).le.5) then
               write(comment,4720) (nt(ir,iz),iz=cizsc,nizs)
               call prc(comment)
            elseif ((nizs-cizsc).le.11) then
               write(comment,4720) (nt(ir,iz),iz=cizsc,cizsc+5)
               call prc(comment)
               write(comment,4720) (nt(ir,iz),iz=cizsc+6,nizs)
               call prc(comment)
            else
               write(comment,4720) (nt(ir,iz),iz=cizsc,cizsc+5)
               call prc(comment)
               write(comment,4720) (nt(ir,iz),iz=cizsc+6,cizsc+11)
               call prc(comment)
               write(comment,4720) (nt(ir,iz),iz=cizsc+11,nizs)
               call prc(comment)
            endif
         end do
         call prb
      endif
c

4690  format (i4,6(1x,g12.7))
4700  format (2x,'Ionization state: ',i4,' content = ',g12.5)
4710  format (18(1x,' IZ= ',i4,1x))
4720  format (18(1x,g10.4))

      return
      end
c
c
c
      subroutine prioniz (isol,ifp,irflct,pionizdat)
      implicit none
      integer isol, irflct, ifp
      real pionizdat(2,2,2,2,5)
c
      include 'params'
      include 'comtor'
      include 'fperiph_com'
c
c     PRIONIZ:
c
c     This routine prints the ionization data selected
c     by the indicators isol, irflct, and ifp for both
c     the "INNER" and "OUTER" target regions.
c
c     The raw data itself is in the array pionizdat -
c     which is compiled in the NEUT subroutine. This
c     data refers to the initial ionization characteristics
c     of the neutrals and does not deal with where these
c     particles eventually end up.
c
c
c
c          Here is a summary of the array and its contents
c
c          PIONIZDAT(2    ,2  ,2    ,2      ,5)
c                   isol  m   ifp   irflct  quant
c
c          isol = 1 for SOL information
c               = 2 for MAIN plasma information
c
c          m    = 1 for target 1   ik > nks(ir)/2
c               = 2 for target 2   ik =< nks(i2)/2
c
c          ifp  = 1 for a neutral resulting from a regular
c                   launch
c               = 2 for a neutral resulting from a Far Periphery
c                   relaunch
c
c          irflct = 1 for a neutral that has NOT been reflected
c                 = 2 for a neutral that has been reflected
c
c          quant= 1 total weight of neutrals
c               = 2 weighted R coordinate
c               = 3 weighted Z coordinate
c               = 4 weighted K value
c               = 5 weighted S value
c
      character*12 soltext(3)
      character*15 fptext(3)
      character*19 rfttext(3)
c
      integer ifplb,irflctlb,isollb,i,j,k
      integer ifpub,irflctub,isolub
      real n1,n2,r1,r2,z1,z2,k1,k2,s1,s2
c
      soltext(1) = ' SOL        '
      soltext(2) = ' MAIN       '
      soltext(3) = ' SOL + MAIN '
c
      fptext(1)  = ' ORDINARY      '
      fptext(2)  = ' FP            '
      fptext(3)  = ' ORDINARY + FP '
c
      rfttext(1) = ' NON-REFLECTED     '
      rfttext(2) = ' REFLECTED         '
      rfttext(3) = ' ALL (REF+NON-REF) '
c
c---------------------------------------------------------------------
c
c
c
c     IF all of the indices are less than 3 - i.e. no summary data
c     requested.
c
      if (isol.lt.3.and.irflct.lt.3.and.ifp.lt.3) then
         n1 = pionizdat(isol,1,ifp,irflct,1)
         n2 = pionizdat(isol,2,ifp,irflct,1)
c
         r1 = pionizdat(isol,1,ifp,irflct,2)
         r2 = pionizdat(isol,2,ifp,irflct,2)
c
         z1 = pionizdat(isol,1,ifp,irflct,3)
         z2 = pionizdat(isol,2,ifp,irflct,3)
c
         k1 = pionizdat(isol,1,ifp,irflct,4)
         k2 = pionizdat(isol,2,ifp,irflct,4)
c
         s1 = pionizdat(isol,1,ifp,irflct,5)
         s2 = pionizdat(isol,2,ifp,irflct,5)
c
      else
c
         n1 = 0.0
         n2 = 0.0
         r1 = 0.0
         r2 = 0.0
         z1 = 0.0
         z2 = 0.0
         k1 = 0.0
         k2 = 0.0
         s1 = 0.0
         s2 = 0.0
c
         if (isol.eq.3) then
            isollb = 1
            isolub = 2
         else
            isollb = isol
            isolub = isol
         endif
c
         if (irflct.eq.3) then
            irflctlb = 1
            irflctub = 2
         else
            irflctlb = irflct
            irflctub = irflct
         endif
c
         if (ifp.eq.3) then
            ifplb = 1
            ifpub = 2
         else
            ifplb = ifp
            ifpub = ifp
         endif
c
         do i = isollb,isolub
            do j = ifplb,ifpub
               do k = irflctlb,irflctub
                  n1 = n1 + pionizdat(i,1,j,k,1)
                  n2 = n2 + pionizdat(i,2,j,k,1)
                  r1 = r1 + pionizdat(i,1,j,k,2)
                  r2 = r2 + pionizdat(i,2,j,k,2)
                  z1 = z1 + pionizdat(i,1,j,k,3)
                  z2 = z2 + pionizdat(i,2,j,k,3)
                  k1 = k1 + pionizdat(i,1,j,k,4)
                  k2 = k2 + pionizdat(i,2,j,k,4)
                  s1 = s1 + pionizdat(i,1,j,k,5)
                  s2 = s2 + pionizdat(i,2,j,k,5)
               end do
            end do
         end do
c
      endif
c
c
      if (n1.gt.0.0.or.n2.gt.0.0) then
c
         call prc ('SUMMARY OF NEUTRAL IONIZATION FOR'
     >               //soltext(isol))
c
         if (fpropt.eq.1) then
            call prc (' FOR: '//fptext(ifp)//' LAUNCHED NEUTRALS')
         endif
c
         if (nrfopt.eq.1.or.nrfopt.eq.2) then
            call prc (' FOR: '//rfttext(irflct)//' PARTICLES')
         endif
c
         call prr2 (' NUMBER OF PARTICLES IN DATA  ',n1,n2)
c
         CALL PRR2 (' AVERAGE R FOR IONISATIONS    ',
     >                     r1/MAX(LO,n1),r2/MAX(LO,n2))
         CALL PRR2 (' AVERAGE Z FOR IONISATIONS    ',
     >                     z1/MAX(LO,n1),z2/MAX(LO,n2))
         CALL PRR2 (' AVERAGE K FOR IONISATIONS    ',
     >                     k1/MAX(LO,n1),k2/MAX(LO,n2))
c
c     Don't print Average S for data involving MAIN ionizations
c
         if (isol.eq.1) then
            call prr2 (' AVERAGE <S> FOR IONISATIONS  ',
     >                     s1/MAX(LO,n1),s2/MAX(LO,n2))
         endif
c
      endif
c
      return
      end
c
c
c
      subroutine probescan
      implicit none
c
c     PROBESCAN: Calculate the background plasma values along
c                a specified vertical and horizontal line at
c                specified R and Z coordinates. The vertical
c                probe corresponds to a CMOD diagnostic. The
c                horizontal to a JET midplane probe.
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     Local variable
c
      real r,z,s,rscan,d1,d2,d3,d4,d5,isat,zscan
      real te,ti,ne
      real fact, factc
      integer ir,ik,isection,sectcnt
      character*256 comment
      logical first,found
c
c
      call prb
      call prchtml('--- RECIPROCATING PROBE RESULTS ---',
     >             'pr_rcp','0','B')
      call prb
c
      rscan = crploc
      zscan = czploc
c
c     Different implementation for polygon based grids.
c
      if (cgridopt.eq.0.or.cgridopt.eq.3.and.pdopt.eq.1) then
c
c        Do vertical probe and then horizontal probe. If the scanning
c        value is -99.0 then ignore that scan. This is just an
c        alternate method to turn the table off.
c
c
         if (rscan.ne.-99.0.and.rlocnum.gt.0) then
c
            first = .true.
            isection = rlocnum
c
            do ir = 2,irwall-1
c
               found = .false.
               sectcnt = 0
c
               do ik = 1,nks(ir)
c
c                 Various distances
c
                  d1 = abs(rscan-krb(ik-1,ir))
                  d2 = abs(rscan-krb(ik,ir))
c
                  d3 = abs(rscan-rs(ik,ir))
c
                  d4 = abs(rs(ik,ir)-krb(ik-1,ir))
c
                  d5 = abs(krb(ik,ir)-rs(ik,ir))
c
c                 Intersection in first half of cell
c
                  if (d1.le.d4.and.d3.le.d4) then

                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then

c
c                       Intersection in first half cell
c
                        if (ik.eq.1) then


                           factc = d1/d4
                           te = kteds(idds(ir,2)) + factc *
     >                         (ktebs(1,ir)-kteds(idds(ir,2)))
                           ti = ktids(idds(ir,2)) + factc *
     >                         (ktibs(1,ir)-ktids(idds(ir,2)))
                           ne = knds(idds(ir,2)) + factc *
     >                         (knbs(1,ir)-knds(idds(ir,2)))
                           s = ksb(ik-1,ir) + factc *
     >                        (kss(ik,ir)-ksb(ik-1,ir))
                           z = kzb(ik-1,ir) + factc *
     >                        (zs(ik,ir)-kzb(ik-1,ir))
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 10
c
c                       Intersection in any other FIRST half cell
c
                        else
c
                           factc = d1/d4
                           fact  = (d1/d4 * (kps(ik,ir)-kpb(ik-1,ir))
     >                              +  kpb(ik-1,ir)-kps(ik-1,ir)) /
     >                               (kps(ik,ir)-kps(ik-1,ir))
c
                           te = ktebs(ik-1,ir) + fact *
     >                         (ktebs(ik,ir)-ktebs(ik-1,ir))
                           ti = ktibs(ik-1,ir) + fact *
     >                         (ktibs(ik,ir)-ktibs(ik-1,ir))
                           ne = knbs(ik-1,ir) + fact *
     >                         (knbs(ik,ir)-knbs(ik-1,ir))
c
                           s  = ksb(ik-1,ir) + factc *
     >                         (kss(ik,ir)-ksb(ik-1,ir))
                           z  = kzb(ik-1,ir) + factc *
     >                         (zs(ik,ir)-kzb(ik-1,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 10

                        endif
c

                     endif

c
c                 Intersection in second half of cell.
c
                  elseif (d2.le.d5.and.d3.le.d5) then


                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then

c
c                       Intersection in first half cell
c
                        if (ik.eq.nks(ir)) then

                           factc = d3/d5

                           te = ktebs(ik,ir) + factc *
     >                         (kteds(idds(ir,1))-ktebs(ik,ir))
                           ti = ktibs(ik,ir) + factc *
     >                         (ktids(idds(ir,1))-ktibs(ik,ir))
                           ne = knbs(ik,ir) + factc *
     >                         (knds(idds(ir,1))- knbs(ik,ir))
c
                           s  = kss(ik,ir) + factc *
     >                         (ksb(ik,ir)-kss(ik,ir))
                           z  = zs(ik,ir) + factc *
     >                         (kzb(ik,ir)-zs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
                           found = .true.
                           goto 10
c
c                       Intersection in any other FIRST half cell
c
                        else
c
                           factc =  d3/d5
                           fact  = (d3/d5 * (kpb(ik,ir)-kps(ik,ir))
     >                              +  kps(ik+1,ir)-kpb(ik,ir)) /
     >                               (kps(ik+1,ir)-kps(ik,ir))
c
                           te = ktebs(ik,ir) + fact *
     >                         (ktebs(ik+1,ir)-ktebs(ik,ir))
                           ti = ktibs(ik,ir) + fact *
     >                         (ktibs(ik+1,ir)-ktibs(ik,ir))
                           ne = knbs(ik,ir) + fact *
     >                         (knbs(ik+1,ir)-knbs(ik,ir))
c
                           s  = kss(ik,ir) + factc *
     >                         (ksb(ik,ir)-kss(ik,ir))
                           z  = zs(ik,ir) + factc *
     >                         (kzb(ik,ir)-zs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 10
c
                        endif
c
                     endif
c
                  endif
c
c              END DO FOR IK
c
               end do
c
 10            if (found) then
c
                  if (first) then
                     write (6,*) 'Table of Reciprocating'//
     >                      ' Probe diagnostic data: R= ',rscan
                     write (6,*) '- For Intersection Number = ',
     >                             rlocnum
                     write (6,200)
c
                     if (cprint.ne.1.and.cprint.ne.9) first = .false.
c
                  endif
c
                  write(6,300) ik,ir,te,ti,ne,isat,s,z,psitarg(ir,1)
c
                  if (cprint.eq.1.or.cprint.eq.9) then
                     if (first) then
                        call prr('Table of Probe Diagnostic data: R= ',
     >                           rscan)
                        call pri ('- For Intersection Number = ',
     >                             rlocnum)
c
                        write (comment,200)
                        call prc(comment)
c
                        first = .false.
c
                     endif
c
                     write(comment,300) ik,ir,te,ti,ne,isat,s,z,
     >                             psitarg(ir,1)
                     call prc(comment)
c
                  endif
c
               endif
c
c           End DO for IR
c
            end do
c
c        End if for RSCAN.ne.-99.0
c
         end if

c
c        Now - do the horizontal probe scan for a fixed Z.
c

         if (zscan.ne.-99.0.and.zlocnum.gt.0) then
c
            first = .true.
            isection = zlocnum
c
            do ir = 2,irwall-1
c
               found = .false.
               sectcnt = 0
c
               do ik = 1,nks(ir)
c
c                 Various distances
c
                  d1 = abs(zscan-kzb(ik-1,ir))
c
                  d2 = abs(zscan-kzb(ik,ir))
c
                  d3 = abs(zscan-zs(ik,ir))
c
                  d4 = abs(zs(ik,ir)-kzb(ik-1,ir))
c
                  d5 = abs(kzb(ik,ir)-zs(ik,ir))
c
c                 Intersection in first half of cell
c
                  if (d1.le.d4.and.d3.le.d4) then

                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then
c
c                       Intersection in first half cell
c
                        if (ik.eq.1) then
c
                           factc = d1/d4
c
                           te = kteds(idds(ir,2)) + factc *
     >                         (ktebs(1,ir)-kteds(idds(ir,2)))
                           ti = ktids(idds(ir,2)) + factc *
     >                         (ktibs(1,ir)-ktids(idds(ir,2)))
                           ne = knds(idds(ir,2)) + factc *
     >                         (knbs(1,ir)-knds(idds(ir,2)))
                           s = ksb(ik-1,ir) + factc *
     >                        (kss(ik,ir)-ksb(ik-1,ir))
                           r = krb(ik-1,ir) + factc *
     >                        (rs(ik,ir)-krb(ik-1,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 20
c
c                       Intersection in any other FIRST half cell
c
                        else
c
                           factc = d1/d4
                           fact  = (d1/d4 * (kps(ik,ir)-kpb(ik-1,ir))
     >                              +  kpb(ik-1,ir)-kps(ik-1,ir)) /
     >                               (kps(ik,ir)-kps(ik-1,ir))
c
                           te = ktebs(ik-1,ir) + fact *
     >                         (ktebs(ik,ir)-ktebs(ik-1,ir))
                           ti = ktibs(ik-1,ir) + fact *
     >                         (ktibs(ik,ir)-ktibs(ik-1,ir))
                           ne = knbs(ik-1,ir) + fact *
     >                         (knbs(ik,ir)-knbs(ik-1,ir))
c
                           s  = ksb(ik-1,ir) + factc *
     >                         (kss(ik,ir)-ksb(ik-1,ir))
                           r  = krb(ik-1,ir) + factc *
     >                         (rs(ik,ir)-krb(ik-1,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 20

                        endif
c

                     endif

c
c                 Intersection in second half of cell.
c
                  elseif (d2.le.d5.and.d3.le.d5) then
c
                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then
c
c                       Intersection in first half cell
c
                        if (ik.eq.nks(ir)) then

                           factc = d3/d5

                           te = ktebs(ik,ir) + factc *
     >                         (kteds(idds(ir,1))-ktebs(ik,ir))
                           ti = ktibs(ik,ir) + factc *
     >                         (ktids(idds(ir,1))-ktibs(ik,ir))
                           ne = knbs(ik,ir) + factc *
     >                         (knds(idds(ir,1))- knbs(ik,ir))
c
                           s  = kss(ik,ir) + factc *
     >                         (ksb(ik,ir)-kss(ik,ir))
                           r  = rs(ik,ir) + factc *
     >                         (krb(ik,ir)-rs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
                           found = .true.
                           goto 20
c
c                       Intersection in any other FIRST half cell
c
                        else
c
                           factc =  d3/d5
                           fact  = (d3/d5 * (kpb(ik,ir)-kps(ik,ir))
     >                              +  kps(ik+1,ir)-kpb(ik,ir)) /
     >                               (kps(ik+1,ir)-kps(ik,ir))
c
                           te = ktebs(ik,ir) + fact *
     >                         (ktebs(ik+1,ir)-ktebs(ik,ir))
                           ti = ktibs(ik,ir) + fact *
     >                         (ktibs(ik+1,ir)-ktibs(ik,ir))
                           ne = knbs(ik,ir) + fact *
     >                         (knbs(ik+1,ir)-knbs(ik,ir))
c
                           s  = kss(ik,ir) + factc *
     >                         (ksb(ik,ir)-kss(ik,ir))
                           r  = rs(ik,ir) + factc *
     >                         (krb(ik,ir)-rs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                           isat = ech*ne*
     >                            9.79e3*SQRT(0.5*(te+ti)
     >                            *(1.0+RIZB)/CRMB)
c
                           found = .true.
                           goto 20
c
                        endif
c
                     endif
c
                  endif
c
c              END DO FOR IK
c
               end do
c
 20            if (found) then
c
                  if (first) then
                     write (6,*) 'Table of Reciprocating'//
     >                      ' Probe diagnostic data: Z= ',zscan
                     write (6,*) '- For Intersection Number = ',
     >                             zlocnum
                     write (6,201)
c
                     if (cprint.ne.1.and.cprint.ne.9) first = .false.
c
                  endif
c
                  write(6,300) ik,ir,te,ti,ne,isat,s,r,psitarg(ir,1)
c
                  if (cprint.eq.1.or.cprint.eq.9) then
                     if (first) then
                        call prr('Table of Probe Diagnostic data: Z= ',
     >                           zscan)
                        call pri ('- For Intersection Number = ',
     >                             zlocnum)
c
                        write (comment,201)
                        call prc(comment)
c
                        first = .false.
c
                     endif
c
                     write(comment,300) ik,ir,te,ti,ne,isat,s,r,
     >                                  psitarg(ir,1)
                     call prc(comment)
c
                  endif
c
               endif
c
c           End DO for IR
c
            end do
c
         end if
c
c     Else - for grids without polygon data and all the extra grid
c     quantities calculated
c
      else
c
c        Scanning for both probes - regular grid
c
         if (rscan.ne.-99.0.and.rlocnum.gt.0) then
c
            first = .true.
            isection = rlocnum
c
            do ir = irsep,irwall-1
c
               found = .false.
               sectcnt = 0
c
               do ik = 1,nks(ir)-1
c
c                 Various distances
c
                  d1 = abs(rscan-rs(ik,ir))
c
                  d2 = abs(rscan-rs(ik+1,ir))
c
                  d3 = abs(rs(ik+1,ir)-rs(ik,ir))
c
                  if (d1.le.d3.and.d2.le.d3) then

                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then
c
c                       Found point - calculate values and
c                       goto next
c
                        te = ktebs(ik,ir) + d1/d3 *
     >                      (ktebs(ik+1,ir)-ktebs(ik,ir))
                        ti = ktibs(ik,ir) + d1/d3 *
     >                      (ktibs(ik+1,ir)-ktibs(ik,ir))
                        ne = knbs(ik,ir) + d1/d3 *
     >                      (knbs(ik+1,ir)-knbs(ik,ir))
                        s  = kss(ik,ir) + d1/d3 *
     >                      (kss(ik+1,ir)-kss(ik,ir))
c
                        z  = zs(ik,ir) + d1/d3 *
     >                      (zs(ik+1,ir)-zs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                        isat = ech*ne*
     >                         9.79e3*SQRT(0.5*(te+ti)
     >                         *(1.0+RIZB)/CRMB)
c
                        found = .true.
                        goto 30
c
                     endif

                  endif
c
c              END DO FOR IK
c
               end do
c
 30            if (found) then
c
                  if (first) then
                     write (6,*) 'Table of Reciprocating'//
     >                      ' Probe diagnostic data: R= ',rscan
                     write (6,*) '- For Intersection Number = ',
     >                             rlocnum
                     write (6,200)
c
                     if (cprint.ne.1.and.cprint.ne.9) first = .false.
c
                  endif
c
                  write(6,300) ik,ir,te,ti,ne,isat,s,z,psitarg(ir,1)
c
                  if (cprint.eq.1.or.cprint.eq.9) then
                     if (first) then
                        call prr('Table of Probe Diagnostic data: R= ',
     >                           rscan)
                        call pri ('- For Intersection Number = ',
     >                             rlocnum)
c
                        write (comment,200)
                        call prc(comment)
c
                        first = .false.
c
                     endif
c
                     write(comment,300) ik,ir,te,ti,ne,isat,s,z,
     >                                  psitarg(ir,1)
                     call prc(comment)
c
                  endif
c
               endif
c
c           End DO for IR
c
            end do
c
c        End if for RSCAN.ne.-99.0
c
         end if
c
c
c        Now - do the horizontal probe scan for a fixed Z.
c
         if (zscan.ne.-99.0.and.zlocnum.gt.0) then
c
            first = .true.
            isection = zlocnum
c
            do ir = irsep,irwall-1
c
               found = .false.
               sectcnt = 0
c
               do ik = 1,nks(ir)
c
c                 Various distances
c
                  d1 = abs(zscan-zs(ik,ir))
c
                  d2 = abs(zscan-zs(ik+1,ir))
c
                  d3 = abs(zs(ik+1,ir)-zs(ik,ir))
c
                  if (d1.le.d3.and.d2.le.d3) then

                     sectcnt = sectcnt + 1
c
                     if (sectcnt.eq.isection) then
c
c                       Found point - calculate values and
c                       goto next
c
                        te = ktebs(ik,ir) + d1/d3 *
     >                      (ktebs(ik+1,ir)-ktebs(ik,ir))
                        ti = ktibs(ik,ir) + d1/d3 *
     >                      (ktibs(ik+1,ir)-ktibs(ik,ir))
                        ne = knbs(ik,ir) + d1/d3 *
     >                      (knbs(ik+1,ir)-knbs(ik,ir))
                        s  = kss(ik,ir) + d1/d3 *
     >                      (kss(ik+1,ir)-kss(ik,ir))
c
                        r  = rs(ik,ir) + d1/d3 *
     >                      (rs(ik+1,ir)-rs(ik,ir))
c
c
c                          jdemod - removed factor of 0.5 - bug
c     
c                           isat = 0.5 * ech*ne*
c
                        isat = ech*ne*
     >                         9.79e3*SQRT(0.5*(te+ti)
     >                         *(1.0+RIZB)/CRMB)
c
                        found = .true.
                        goto 40
c
                     endif

                  endif
c
c              END DO FOR IK
c
               end do
c
 40            if (found) then
c
                  if (first) then
                     write (6,*) 'Table of Reciprocating'//
     >                      ' Probe diagnostic data: Z= ',zscan
                     write (6,*) '- For Intersection Number = ',
     >                             zlocnum
                     write (6,201)
c
                     if (cprint.ne.1.and.cprint.ne.9) first = .false.
c
                  endif
c
                  write(6,300) ik,ir,te,ti,ne,isat,s,r,psitarg(ir,1)
c
                  if (cprint.eq.1.or.cprint.eq.9) then
                     if (first) then
                        call prr('Table of Probe Diagnostic data: Z= ',
     >                           zscan)
                        call pri ('- For Intersection Number = ',
     >                             zlocnum)
c
                        write (comment,201)
                        call prc(comment)
c
                        first = .false.
c
                     endif
c
                     write(comment,300) ik,ir,te,ti,ne,isat,s,r,
     >                                  psitarg(ir,1)
                     call prc(comment)
c
                  endif
c
               endif
c
c           End DO for IR
c
            end do
c
c
         end if
c
      end if
c
c
c     Format statements
c

200   format(3x,'IK',3x,'IR',7X,'Te',8X,'Ti',10x,
     >       'Ne',4x,'Probe_Isat',9x,'s',9x,'Z',9x,'PSIn')
201   format(3x,'IK',3x,'IR',7X,'Te',8X,'Ti',10x,
     >       'Ne',4x,'Probe_Isat',9x,'s',9x,'R',9x,'PSIn')
300   format(2x,i3,2x,i3,3x,f9.3,x,f9.3,x,e13.5,x,e13.5,
     >       x,f8.3,x,f9.5,1x,f9.5)
c
      return
      end
c
c
c
      subroutine radproc(nizs,rions,nimps)
      implicit none
      include    'params'
      include    'dynam1'
      include    'dynam3'
      include    'comtor'
      include    'cgeom'
      include    'commv'
c
      integer nizs,nimps
      real    rions(maxizs)
c
c     RADPROC: This routine processes the impurity ionization array
c              trying to determine specific quantities and
c              characteristics about where the radiation is
c              occuring. These quantities include the volume weighted
c              average of the local densities within the radiating
c              volume. The radiating volume is defined as the volume
c              that radiates the top 2/3's of the energy.
c
      integer ik,ir,in,count(0:maxizs+1),top(0:maxizs+1),num,iz,iz2
c
c     String lengths for print outs
c
      integer len,lenstr
      external lenstr
c
      character*131 coment
c
      real rad(maxnks,maxnrs,0:maxizs+1),radtot(0:maxizs+1)
      real neav(0:maxizs+1),nizav(0:maxizs+1),volav(0:maxizs+1)
      real teav(0:maxizs+1),lzav(0:maxizs+1)
      real pradt(0:maxizs+1),ddtmp
      real frad(0:maxizs+1)
      real sdtot(maxnks,maxnrs)
c
      real radord(maxnks*maxnrs,0:maxizs+1,5)
c
      call rzero(radtot,maxizs+2)
      call rzero(rad,maxnks*maxnrs*(maxizs+2))
      call rzero(sdtot,maxnks*maxnrs)
c
c     Get array containg all Radiation and the TOTAL.
c
      do iz = 0,nizs
         do ir = 1,nrs
            do ik = 1,nks(ir)
               rad(ik,ir,iz) = powls(ik,ir,iz) * kareas(ik,ir)
               radtot(iz) = radtot(iz) + rad(ik,ir,iz)
               rad(ik,ir,nizs+1) = rad(ik,ir,nizs+1) + rad(ik,ir,iz)
               radtot(nizs+1) = radtot(nizs+1) + rad(ik,ir,iz)
               sdtot(ik,ir) = sdtot(ik,ir) + ddlims(ik,ir,iz)
            end do
         end do
      end do
c
c     Sort/order each ionization state and the total
c
c
      do iz = 0,nizs+1
         count(iz) = 0
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (rad(ik,ir,iz).gt.0.0) then
                  count(iz) = count(iz) + 1
                  radord(count(iz),iz,1) = rad(ik,ir,iz)
                  radord(count(iz),iz,2) = -1
                  radord(count(iz),iz,3) = -1
                  radord(count(iz),iz,4) = ik
                  radord(count(iz),iz,5) = ir
               endif
            end do
         end do
      end do
c
c     Now have unsorted array - need to sort it.
c
      call sortrad(radord,count,top,nizs)
c
c     Now that we have sorted arrays - can take the top 2/3 of the
c     total radiated power and calculated avearged densities.
c
      do iz = 0,nizs+1
c
         num = top(iz)
c
         frad(iz)  = 0.0
         pradt(iz) = 0.0
         neav(iz)  = 0.0
         nizav(iz) = 0.0
         teav(iz)  = 0.0
         volav(iz) = 0.0
c
c         if (count(iz).eq.0) cycle
c
         if (count(iz).ne.0) then
c
c         do while(num.ne.-1
c     >               .and.pradt(iz).lt.0.66*radtot(iz))
c
c           do 20 while(num.ne.-1
c     >               .and.pradt(iz).lt.0.66*radtot(iz))
c
 30       if (num.eq.-1.or.pradt(iz).ge.0.66*radtot(iz)) goto 20
c
            pradt(iz) = pradt(iz) + radord(num,iz,1)
c
            ik = radord(num,iz,4)
            ir = radord(num,iz,5)
c
            neav(iz) = neav(iz) + radord(num,iz,1) * knbs(ik,ir)
            teav(iz) = teav(iz) + radord(num,iz,1) * ktebs(ik,ir)
c
            if (iz.eq.nizs+1) then
c
               nizav(iz)= nizav(iz)+radord(num,iz,1)*sdtot(ik,ir)
c
            else
               nizav(iz)= nizav(iz)+radord(num,iz,1)*ddlims(ik,ir,iz)
            endif
c
            volav(iz)= volav(iz)+ kareas(ik,ir)
c
            num = radord(num,iz,2)
c
c         end do
c
            goto 30
c
 20      continue
c
c        Calculate averages
c
         if (volav(iz).ne.0.0) then
            neav(iz) = neav(iz) / pradt(iz)
            nizav(iz)= nizav(iz) / pradt(iz)
            teav(iz) = teav(iz) / pradt(iz)
            lzav(iz) = pradt(iz)/(volav(iz)*neav(iz)*nizav(iz))
c
c           Add checks for conditions that could cause a division by zero
c
            if (iz.eq.0) then
               if (cieizs(iz).le.0.0.or.citizs(iz).le.0.0) then
                  frad(iz) = 0.0
               else
                  frad(iz) = radtot(iz) /
     >                    ( fsrate*cieizs(iz)/citizs(iz) *
     >                      neav(iz) * 1.0)
               endif
            elseif (iz.le.nizs) then
               if (cieizs(iz).le.0.0.or.citizs(iz).le.0.0.or.
     >             rions(iz).le.0.0.or.nimps.le.0) then
                  frad(iz) = 0.0
               else
                  frad(iz) = radtot(iz) /
     >                    ( qtim*cieizs(iz)/citizs(iz) *
     >                      neav(iz) *
     >                      rions(iz)/float(nimps))
               endif
            else
               frad(iz) = 0.0
            endif
c
         else
            neav(iz) = 0.0
            nizav(iz)= 0.0
            teav(iz) = 0.0
            lzav(iz) = 0.0
            frad(iz) = 0.0
         endif
c
c        ENDIF For count(iz) .ne. 0
c
         endif

      end do
c
c     Print out summary of data
c
      write (6,*) 'SUMMARY of Radiation Source: LOCAL Conditions'
      write (6,*)
      write (6,1020)
c
      do iz = 0,nizs+1
c
         if (iz.eq.nizs+1) then
            write(6,1010) radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz)
         else
            write(6,1000) iz,radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz),frad(iz)
         endif
c
      end do
c
      call prc('SUMMARY of Radiation Source: LOCAL Conditions')
      call prb
      write (coment,1021)
      len = lenstr(coment)
      call prc(coment(1:len))
c
      do iz = 0,nizs+1
c
         if (iz.eq.nizs+1) then
            write(coment,1011) radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz)
            len = lenstr(coment)
            call prc(coment(1:len))
         else
            write(coment,1001) iz,radtot(iz),pradt(iz),neav(iz),
     >         nizav(iz)*absfac,
     >         nizav(iz),teav(iz),volav(iz),lzav(iz),frad(iz)
            len = lenstr(coment)
            call prc(coment(1:len))
         endif
c
      end do
c
c     Write ALL ionization states data to output file
c
      call prb
      iz = nizs+1
      call prc('TOTAL Radiation Source: Local Conditions')
      call prr('   Total Radiation         : ',radtot(iz))
      call prr('   Radiation Cutoff        : ',pradt(iz))
      call prr('   Average Ne              : ',neav(iz))
      call prr('   Average Ni (All states) : ',nizav(iz)*absfac)
      call prr('   Raw     Ni (All states) : ',nizav(iz))
      call prr('   Average Te              : ',teav(iz))
      call prr('   Radiating Volume        : ',volav(iz))
      call prr('   Estimated Lz Average    : ',lzav(iz))
      call prb
c
 1020 format (5x,'Ion State',5x,'Total Rad',4x,'Cutoff Rad',
     >             9x,'Ne AV',8x,'Niz AV',7x,'Niz Raw',7x,'Te AV',
     >             7x,'Volume',7x,'Lz EST',6x,'Frad')
 1010 format (8x,'ALL',3x,8(g13.5))
 1000 format (7x,i4,3x,9(g13.5))
c
 1021 format ('Ion',1x,'Tot Rad',1x,'Cut Rad',
     >             3x,'Ne AV',2x,'Niz AV',1x,'Niz Raw',3x,'Te AV',
     >             2x,'Volume',2x,'Lz EST',4x,'Frad')
 1011 format ('ALL',9(1x,1p,g7.1))
 1001 format (i3,9(1x,1p,g7.1))
c
      return
      end
c
c
c
      subroutine prleakage
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'printopt'
c
c     PRLEAKAGE: The purpose of this subroutine is to print out
c                an analysis of the core leakage data. The code
c                now tracks the source type of the particle by
c                wall or target element index, ring in which
c                the particle is ionized, and whether the particle
c                was produced by physical or chemical sputtering.
c
c                This code sums up over this data - producing both
c                element by element output and a summary of
c                leakage by wall and target segment.
c
c
c        Local Variables
c
         real totsource(5)
         real totsleak(5)
c
c         real targsrc(3,4)
c         real targleak(3,4)
c         real wallsrc(5,3)
c         real wallleak(5,3)
c
         real totsrc,totleak
c
         integer id,ir,in,itmp
c
         character*80 comment
c

         call rzero(targsrc,3*4)
         call rzero(targleak,3*4)
         call rzero(wallsrc,5*3)
         call rzero(wallleak,5*3)
c
c        Print out target sources
c
         if (cprint.eq.2.or.cprint.eq.9) then
c
            call prc('SUMMARY OF RING OF IONIZATION FROM'//
     >            ' EACH TARGET ELEMENT:')
            call prc('- ELEMENTS WITH ZERO PRODUCTION ARE'//
     >                 ' NOT LISTED')
            call prb
c
         endif
c
         do id = 1,nds
c
            totsource(1) = 0.0
            totsource(2) = 0.0
            totsource(3) = 0.0
c
            totsleak(1) = 0.0
            totsleak(2) = 0.0
            totsleak(3) = 0.0
c
            do in = 1,3
c
               do ir = 1,nrs
c
c                 Total ionization and leakage
c
                  totsource(in)= totsource(in)+wtsource(id,ir,1,in)
                  totsleak(in) = totsleak(in) +wtsource(id,ir,3,in)
c
c                 Average S for ionized particles
c
                  if (wtsource(id,ir,1,in).gt.0.0) then
                     wtsource(id,ir,2,in) = wtsource(id,ir,2,in)
     >                             / wtsource(id,ir,1,in)
                  else
                     wtsource(id,ir,2,in) = 0.0
                  endif
c
c                 Average starting S for core leaked particles
c
                  if (wtsource(id,ir,3,in).gt.0.0) then
                     wtsource(id,ir,4,in) = wtsource(id,ir,4,in)
     >                             / wtsource(id,ir,3,in)
                  else
                     wtsource(id,ir,4,in) = 0.0
                  endif
c
c                 Collect data for summaries
c
                  if (id.le.ndsin) then
                     targsrc(1,in)=targsrc(1,in)+wtsource(id,ir,1,in)
                     targleak(1,in)=targleak(1,in)+wtsource(id,ir,3,in)
                  else
                     targsrc(2,in)=targsrc(2,in)+wtsource(id,ir,1,in)
                     targleak(2,in)=targleak(2,in)+wtsource(id,ir,3,in)
                  endif
c
               end do
c
            end do
c
c           Print Target souce element by element
c
            if (cprint.eq.2.or.cprint.eq.9) then
c
c              Total Source
c
               if (totsource(1).gt.0.0.or.totsource(2).gt.0.0
     >             .or.totsource(3).gt.0.0) then

                  call prb
                  write (comment,9050) id,irds(id),rp(id),zp(id)
                  call prc(comment)

                  write (comment,9052)
     >               totsource(1)+totsource(2)+totsource(3),
     >               totsleak(1)+totsleak(2)+totsleak(3),
     >               (totsleak(1)+totsleak(2)+totsleak(3))/
     >               (totsource(1)+totsource(2)+totsource(3))
                  call prc(comment)

               endif
c
c              Physical Sputter Source
c
               in = 1
c
               if (totsource(in).gt.0.0) then

                  write (comment,9053)
     >               totsource(in),totsleak(in),
     >               totsleak(in)/totsource(in)
                  call prc(comment)
c
                  do ir = 1,nrs
c
                     if (wtsource(id,ir,1,in).gt.0.0) then

                       write (comment,9057) ir,wtsource(id,ir,1,in),
     >                    wtsource(id,ir,2,in),wtsource(id,ir,3,in),
     >                    wtsource(id,ir,4,in)
                        call prc(comment)

                     endif
c
                  enddo
c
               endif
c
c              Chemical Sputter Source
c
               in = 2
c
               if (totsource(in).gt.0.0) then

                  write (comment,9054)
     >               totsource(in),totsleak(in),
     >               totsleak(in)/totsource(in)
                  call prc(comment)
c
                  do ir = 1,nrs
c
                     if (wtsource(id,ir,1,in).gt.0.0) then

                       write (comment,9057) ir,wtsource(id,ir,1,in),
     >                    wtsource(id,ir,2,in),wtsource(id,ir,3,in),
     >                    wtsource(id,ir,4,in)
                        call prc(comment)

                     endif
c
                  enddo
c
               endif
c
c
c              Self-Sputter Source
c
               in = 3
c
               if (totsource(in).gt.0.0) then

                  write (comment,9058)
     >               totsource(in),totsleak(in),
     >               totsleak(in)/totsource(in)
                  call prc(comment)
c
                  do ir = 1,nrs
c
                     if (wtsource(id,ir,1,in).gt.0.0) then

                       write (comment,9057) ir,wtsource(id,ir,1,in),
     >                    wtsource(id,ir,2,in),wtsource(id,ir,3,in),
     >                    wtsource(id,ir,4,in)
                       call prc(comment)

                     endif
c
                  enddo
c
               endif
c
c           End of Print option
c
            endif
c
         enddo


         if (cprint.eq.2.or.cprint.eq.9) then
c
c           Print OUT WALL sources
c
            call prb
            call prc('SUMMARY OF RING OF IONIZATION FROM'//
     >            ' EACH WALL ELEMENT:')
            call prc('- ELEMENTS WITH ZERO PRODUCTION ARE'//
     >                 ' NOT LISTED')
            call prb
c
         endif
c
         do id = 1,wallpts
c
            totsource(4) = 0.0
            totsource(5) = 0.0
c
            totsleak(4) = 0.0
            totsleak(5) = 0.0
c
            do in = 4,5
c
               do ir = 1,nrs
c
c                 Total ionization and leakage
c
                  totsource(in)= totsource(in)+wtsource(id,ir,1,in)
                  totsleak(in) = totsleak(in) +wtsource(id,ir,3,in)
c
c                 Average S for ionized particles
c
                  if (wtsource(id,ir,1,in).gt.0.0) then
                     wtsource(id,ir,2,in) = wtsource(id,ir,2,in)
     >                             / wtsource(id,ir,1,in)
                  else
                     wtsource(id,ir,2,in) = 0.0
                  endif
c
c                 Average starting S for core leaked particles
c
                  if (wtsource(id,ir,3,in).gt.0.0) then
                     wtsource(id,ir,4,in) = wtsource(id,ir,4,in)
     >                             / wtsource(id,ir,3,in)
                  else
                     wtsource(id,ir,4,in) = 0.0
                  endif
c
c                 Collect Summary information
c
c                 1 = Inner target
c                 4 = Outer target
c                 2,3,5,6,7 = Main wall
c                 8 = PP Wall
c
c
                  itmp = in -3
c
c                 Inner target
c
                  if (wallpt(id,16).eq.4) then
c
                     wallsrc(1,itmp)=wallsrc(1,itmp)
     >                              +wtsource(id,ir,1,in)
                     wallleak(1,itmp)=wallleak(1,itmp)
     >                               +wtsource(id,ir,3,in)
c
c                 Outer Target
c
                  elseif (wallpt(id,16).eq.1) then
c
                     wallsrc(2,itmp)=wallsrc(2,itmp)
     >                              +wtsource(id,ir,1,in)
                     wallleak(2,itmp)=wallleak(2,itmp)
     >                               +wtsource(id,ir,3,in)
c
c                 PP Wall
c
                  elseif (wallpt(id,16).eq.8) then
c
                     wallsrc(3,itmp)=wallsrc(3,itmp)
     >                              +wtsource(id,ir,1,in)
                     wallleak(3,itmp)=wallleak(3,itmp)
     >                               +wtsource(id,ir,3,in)
c
c                 Main Wall
c
                  elseif (wallpt(id,16).eq.2.or.
     >                    wallpt(id,16).eq.3.or.
     >                    wallpt(id,16).eq.5.or.
     >                    wallpt(id,16).eq.6.or.
     >                    wallpt(id,16).eq.7) then
c
                     wallsrc(4,itmp)=wallsrc(4,itmp)
     >                              +wtsource(id,ir,1,in)
                     wallleak(4,itmp)=wallleak(4,itmp)
     >                               +wtsource(id,ir,3,in)
c
                  endif
c
               end do
c
            end do
c

            if (cprint.eq.2.or.cprint.eq.9) then
c
c              Total Source
c
               if (totsource(4).gt.0.0.or.totsource(5).gt.0.0) then
c
                  call prb
c
                  if (wallpt(id,18).le.0.0) then
c
                     write (comment,9051) id,wallpt(id,1),wallpt(id,2)
                     call prc(comment)
c
                  elseif (wallpt(id,18).gt.0.0) then
c
                     write (comment,9059) id,wallpt(id,1),wallpt(id,2),
     >                                 int(wallpt(id,18))
                     call prc(comment)
c
                  endif
c
                  write (comment,9053)
     >               totsource(4)+totsource(5),
     >               totsleak(4)+totsleak(5),
     >               (totsleak(4)+totsleak(5))/
     >               (totsource(4)+totsource(5))
                  call prc(comment)

               endif
c
c              Physical Sputter Source
c
               in = 4
c
               if (totsource(in).gt.0.0) then

                  write (comment,9053)
     >               totsource(in),totsleak(in),
     >               totsleak(in)/totsource(in)
                  call prc(comment)
c
                  do ir = 1,nrs
c
                     if (wtsource(id,ir,1,in).gt.0.0) then

                       write (comment,9057) ir,wtsource(id,ir,1,in),
     >                    wtsource(id,ir,2,in),wtsource(id,ir,3,in),
     >                    wtsource(id,ir,4,in)
                       call prc(comment)

                     endif
c
                  enddo
c
               endif
c
c              Chemical Sputter Source
c
               in = 5
c
               if (totsource(in).gt.0.0) then

                  write (comment,9054)
     >               totsource(in),totsleak(in),
     >               totsleak(in)/totsource(in)
                  call prc(comment)
c
                  do ir = 1,nrs
c
                     if (wtsource(id,ir,1,in).gt.0.0) then

                       write (comment,9057) ir,wtsource(id,ir,1,in),
     >                    wtsource(id,ir,2,in),wtsource(id,ir,3,in),
     >                    wtsource(id,ir,4,in)
                       call prc(comment)

                     endif
c
                  enddo
c
               endif
c
c          End of Print Option
c
           endif
c
        end do
c
c
c       Calculate and Print Leakage Summaries -
c
c       targsrc(x,y)  ->  x = 1 = inner target
c                         x = 2 = outer target
c                         x = 3 = totals over x
c
c                         y = 1 = physical
c                         y = 2 = chemical
c                         y = 3 = self
c                         y = 4 = totals over y
c
c       wallsrc(x,y)  ->  x = 1 = inner target
c                         x = 2 = outer target
c                         x = 3 = pp wall
c                         x = 4 = main wall
c                         x = 5 = totals over x
c
c                         y = 1 = physical
c                         y = 2 = chemical
c                         y = 3 = totals over y
c
c
c       Sum over x
c
        do in = 1,3
c
           targsrc(3,in) = targsrc(1,in)+targsrc(2,in)
           targleak(3,in) = targleak(1,in)+targleak(2,in)
c
        end do
c
c       Sum over y
c
        do in = 1,2
c
           targsrc(in,4) = targsrc(in,1)+targsrc(in,2)
     >                    +targsrc(in,3)
           targleak(in,4) = targleak(in,1)+targleak(in,2)
     >                     +targleak(in,3)
c
           wallsrc(5,in) = wallsrc(1,in) + wallsrc(2,in)+
     >                     wallsrc(3,in) + wallsrc(4,in)
           wallleak(5,in) = wallleak(1,in) + wallleak(2,in)+
     >                      wallleak(3,in) + wallleak(4,in)
c
        end do
c
        targsrc(3,4) = targsrc(3,1)+targsrc(3,2)+targsrc(3,3)
        targleak(3,4)= targleak(3,1)+targleak(3,2)+targleak(3,3)
c
        do in = 1,4
c
           wallsrc(in,3) = wallsrc(in,1) + wallsrc(in,2)
           wallleak(in,3) = wallleak(in,1) + wallleak(in,2)
c
        end do
c
        wallsrc(5,3) = wallsrc(5,1) + wallsrc(5,2)
        wallleak(5,3) = wallleak(5,1) + wallleak(5,2)
c
        totsrc = targsrc(3,4) + wallsrc(5,3)
        totleak = targleak(3,4) + wallleak(5,3)
c
c
c       Only print out if there is actually some core leakage.
c
        if (totsrc.gt.0.0.and.totleak.gt.0.0) then
c
           call prb
           call prc('SUMMARY OF LEAKAGE DATA:')
           call prb
c
c           if (zxp.lt.z0) then
c              call prc('********************************')
c              call prc('        WARNING!!!              ')
c              call prc(' Meaning of Inner/Outer MUST be ')
c              call prc(' reversed for X-point down grids')
c              call prc('********************************')
c           endif
c
           call prb
           write(comment,'(a,e12.3)')
     >                 ' Total Source Strength :',totsrc
           call prc(comment)
c
           write(comment,'(a,e12.3)')
     >              ' Total Leakage         :',totleak
           call prc(comment)
c
           write(comment,'(a,e12.3)')
     >              ' Total Percent Leakage :',
     >                totleak/totsrc*100.0
           call prc(comment)
c
           call prb
c
           call prc(' Source and Leakage Breakdown:')
c
c          Target Source Leakage
c
           if (targsrc(3,4).gt.0.0) then
              call prb
              call prc (' Source due to target ION flux:')
              call prb
              write(comment,9070)
              call prc(comment)
           endif
c
c          Inner Target
c
           if (targsrc(1,4).gt.0.0) then
c
c             Inner Target - total
c
              write(comment,9060) inner,targsrc(1,4),targleak(1,4),
     >            targsrc(1,4)/totsrc,targleak(1,4)/totleak
              call prc(comment)
c
c             Inner Target - Physical
c
              write(comment,9064) targsrc(1,1),targleak(1,1),
     >            targsrc(1,1)/totsrc,targleak(1,1)/totleak
              call prc(comment)
c
c             Inner Target - Chemical
c
              write(comment,9065) targsrc(1,2),targleak(1,2),
     >            targsrc(1,2)/totsrc,targleak(1,2)/totleak
              call prc(comment)
c
c             Inner Target - Self
c
              write(comment,9068) targsrc(1,3),targleak(1,3),
     >            targsrc(1,3)/totsrc,targleak(1,3)/totleak
              call prc(comment)
c
           endif
c
c          Outer Target
c
           if (targsrc(2,4).gt.0.0) then
c
c             Outer Target - total
c
              write(comment,9060) outer,targsrc(2,4),targleak(2,4),
     >            targsrc(2,4)/totsrc,targleak(2,4)/totleak
              call prc(comment)
c
c             Outer Target - Physical
c
              write(comment,9064) targsrc(2,1),targleak(2,1),
     >            targsrc(2,1)/totsrc,targleak(2,1)/totleak
              call prc(comment)
c
c             Outer Target - Chemical
c
              write(comment,9065) targsrc(2,2),targleak(2,2),
     >            targsrc(2,2)/totsrc,targleak(2,2)/totleak
              call prc(comment)
c
c             Outer Target - Self
c
              write(comment,9068) targsrc(2,3),targleak(2,3),
     >            targsrc(2,3)/totsrc,targleak(2,3)/totleak
              call prc(comment)
c
           endif
c
c          Totals
c
           if (targsrc(3,4).gt.0.0) then
c
c             Target Total
c
              write(comment,9066) targsrc(3,4),targleak(3,4),
     >            targsrc(3,4)/totsrc,targleak(3,4)/totleak
              call prc(comment)
c
c             Target - Physical
c
              write(comment,9064) targsrc(3,1),targleak(3,1),
     >            targsrc(3,1)/totsrc,targleak(3,1)/totleak
              call prc(comment)
c
c             Target - Chemical
c
              write(comment,9065) targsrc(3,2),targleak(3,2),
     >            targsrc(3,2)/totsrc,targleak(3,2)/totleak
              call prc(comment)
c
c             Target - Self
c
              write(comment,9068) targsrc(3,3),targleak(3,3),
     >            targsrc(3,3)/totsrc,targleak(3,3)/totleak
              call prc(comment)
c
          endif
c
c--------------------------------------------------------------
c
c
c          Wall Source leakage
c
           if (wallsrc(5,3).gt.0.0) then
              call prb
              call prc (' Source due to wall and target ATOM flux:')
              call prb
              write(comment,9070)
              call prc(comment)
           endif
c
c          Inner Target
c
           if (wallsrc(1,3).gt.0.0) then
c
c             Inner Target - total
c
              write(comment,9060) inner,wallsrc(1,3),wallleak(1,3),
     >            wallsrc(1,3)/totsrc,wallleak(1,3)/totleak
              call prc(comment)
c
c             Inner Target - Physical
c
              write(comment,9064) wallsrc(1,1),wallleak(1,1),
     >            wallsrc(1,1)/totsrc,wallleak(1,1)/totleak
              call prc(comment)
c
c             Inner Target - Chemical
c
              write(comment,9065) wallsrc(1,2),wallleak(1,2),
     >            wallsrc(1,2)/totsrc,wallleak(1,2)/totleak
              call prc(comment)
c
           endif
c
c          Outer Target
c
           if (wallsrc(2,3).gt.0.0) then
c
c             Outer Target - total
c
              write(comment,9060) outer,wallsrc(2,3),wallleak(2,3),
     >            wallsrc(2,3)/totsrc,wallleak(2,3)/totleak
              call prc(comment)
c
c             Outer Target - Physical
c
              write(comment,9064) wallsrc(2,1),wallleak(2,1),
     >            wallsrc(2,1)/totsrc,wallleak(2,1)/totleak
              call prc(comment)
c
c             Outer Target - Chemical
c
              write(comment,9065) wallsrc(2,2),wallleak(2,2),
     >            wallsrc(2,2)/totsrc,wallleak(2,2)/totleak
              call prc(comment)
c
           endif
c
c          PP Wall
c
           if (wallsrc(3,3).gt.0.0) then
c
c             PP Wall - Total
c
              write(comment,9062) wallsrc(3,3),wallleak(3,3),
     >            wallsrc(3,3)/totsrc,wallleak(3,3)/totleak
              call prc(comment)
c
c             PP Wall - Physical
c
              write(comment,9064) wallsrc(3,1),wallleak(3,1),
     >            wallsrc(3,1)/totsrc,wallleak(3,1)/totleak
              call prc(comment)
c
c             PP Wall - Chemical
c
              write(comment,9065) wallsrc(3,2),wallleak(3,2),
     >            wallsrc(3,2)/totsrc,wallleak(3,2)/totleak
              call prc(comment)
c
           endif
c
c          Main Wall
c
           if (wallsrc(4,3).gt.0.0) then
c
c             Main Wall - Total
c
              write(comment,9063) wallsrc(4,3),wallleak(4,3),
     >            wallsrc(4,3)/totsrc,wallleak(4,3)/totleak
              call prc(comment)
c
c             Main Wall - Physical
c
              write(comment,9064) wallsrc(4,1),wallleak(4,1),
     >            wallsrc(4,1)/totsrc,wallleak(4,1)/totleak
              call prc(comment)
c
c             Main Wall - Chemical
c
              write(comment,9065) wallsrc(4,2),wallleak(4,2),
     >            wallsrc(4,2)/totsrc,wallleak(4,2)/totleak
              call prc(comment)
c
           endif
c
c          Totals
c
           if (wallsrc(5,3).gt.0.0) then
c
c             Wall Total
c
              write(comment,9067) wallsrc(5,3),wallleak(5,3),
     >            wallsrc(5,3)/totsrc,wallleak(5,3)/totleak
              call prc(comment)
c
c             Wall - Physical
c
              write(comment,9064) wallsrc(5,1),wallleak(5,1),
     >            wallsrc(5,1)/totsrc,wallleak(5,1)/totleak
              call prc(comment)
c
c             Wall - Chemical
c
              write(comment,9065) wallsrc(5,2),wallleak(5,2),
     >            wallsrc(5,2)/totsrc,wallleak(5,2)/totleak
              call prc(comment)
c
          endif
c
       else
c
          call prc('NO CORE LEAKAGE -'//
     >            ' SOURCE ANALYSIS NOT PRINTED')

       endif
c
       call prb
c
c
c     Format statements
c
 9050 format('Target-ID:',i4,' Ring:',i4,'  R:',f8.3,'  Z:',f8.3)
 9051 format('Wall-ID:',i4,'  R:',f8.3,'  Z:',f8.3)
 9059 format('Wall-ID:',i4,'  R:',f8.3,'  Z:',f8.3,
     >       4x,'Matching Target-ID:',i4)
c
 9052 format('Totals  -Launched:',f8.3,'  Leak:',f8.3,
     >       '  Prob:',f8.5)
 9053 format('  Phys  -Launched:',f8.3,'  Leak:',f8.3,
     >       '  Prob:',f8.5)
 9054 format('  Chem  -Launched:',f8.3,'  Leak:',f8.3,
     >       '  Prob:',f8.5)
 9058 format('  Self  -Launched:',f8.3,'  Leak:',f8.3,
     >       '  Prob:',f8.5)
c
 9055 format(4x,'Ring',7x,'Number',7x,'Average S',5x,
     >       'Leaked to Core',4x,'Av. S Leak')
 9056 format(4x,i4,3x,f13.5,4x,f9.3,4x,f9.3,4x,f9.3)
c
 9057 format('    Ring:',i4,'  Num:',f8.3,'  S:',f8.3,
     > '  Leak:',f8.3,'  Av S:',f8.3)
c
c IPP/08 Krieger - changed format statements
 9060 format(' ',a5,' TARGET:',2(1x,g12.3),2(1x,g12.3))
c 9061 format(' ',a5,' TARGET:',2(1x,f9.3),2(1x,f9.5))
 9062 format(' PP WALL     :',2(1x,g12.3),2(1x,g12.3))
 9063 format(' MAIN WALL   :',2(1x,g12.3),2(1x,g12.3))
 9064 format('   -PHYSICAL :',2(1x,g12.3),2(1x,g12.3))
 9065 format('   -CHEMICAL :',2(1x,g12.3),2(1x,g12.3))
 9068 format('   -SELF     :',2(1x,g12.3),2(1x,g12.3))
 9066 format(' TOTAL ION   :',2(1x,g12.3),2(1x,g12.3))
 9067 format(' TOTAL ATOM  :',2(1x,g12.3),2(1x,g12.3))
 9070 format(18x,'Source',6x,'Leak',2x,'% Source',4x,'% Leak')
c
c 9060 format(' ',a5,' TARGET:',2(1x,f9.3),2(1x,f9.5))
c 9061 format(' ',a5,' TARGET:',2(1x,f9.3),2(1x,f9.5))
c 9062 format(' PP WALL     :',2(1x,f9.3),2(1x,f9.5))
c 9063 format(' MAIN WALL   :',2(1x,f9.3),2(1x,f9.5))
c 9064 format('   -PHYSICAL :',2(1x,f9.3),2(1x,f9.5))
c 9065 format('   -CHEMICAL :',2(1x,f9.3),2(1x,f9.5))
c 9068 format('   -SELF     :',2(1x,f9.3),2(1x,f9.5))
c 9066 format(' TOTAL ION   :',2(1x,f9.3),2(1x,f9.5))
c 9067 format(' TOTAL ATOM  :',2(1x,f9.3),2(1x,f9.5))
c 9070 format(18x,'Source',6x,'Leak',2x,'% Source',4x,'% Leak')
c
c     End of PRLEAKAGE
c
      return
      end
c
c
c
      subroutine promptdep(ik,ir,id,r,z,riz,sputy,massi,temi,
     >                     sheath_drop,rc)
      implicit none
      integer ik,ir,rc,id
      real r,z,temi,sheath_drop,riz,massi,sputy
c
c     Common blocks
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'promptdep'
c
c     PROMPTDEP: This routine estimates if the position (R,Z)
c                of an ion in charge state RIZ is within one
c                Larmor Radius of the target. It then estimates
c                the expected impact energy of this particle on
c                the target.
c
c     David Elder,    November 1997
c
c     Assumptions: 1) The distance to the target is defined by
c                     the closest target element on the ring
c                     of ionization. This will break down if
c                     the particle is far from the target -
c                     however this is not expected to be the
c                     case for promptly redeposited ions.
c
c
c
c     Local Variables
c
      integer ik_local,ir_local,it
      logical griderr
      real    larmor_radius
      real    larmor,b_field
      real    targ_dist,dist_to_point
      external dist_to_point
      external larmor
c slmod begin - tmp
      LOGICAL getrz_error
c slmod end
c
c     First check the particle grid position
c
      rc = 0
      ik_local = ik
      ir_local = ir
c
      call gridpos(ik_local,ir_local,r,z,.false.,griderr)
c
c     Find nearest target element - assume for now it is closest target
c     on current ring.
c
      if (ik_local.lt.nks(ir_local)/2) then
         it = 2
      else
         it = 1
      endif
c
c     Define target index
c
      id = idds(ir_local,it)
c slmod begin
      getrz_error = .FALSE.
      IF (id.EQ.0) THEN
        getrz_error = .TRUE.
        WRITE(0,*) 'WARNING promptdep: getrz_confusion, prompt '//
     .             'redeposition check lost'
c        WRITE(0,*) 'WHOA! PROBLEM!'
c        WRITE(0,*) griderr
c        WRITE(0,*) r,z
c        WRITE(0,*) ik_local,ir_local
c        WRITE(0,*) ik,ir
c        WRITE(0,*) it
        id = idds(irsep,2)
      ENDIF
c slmod end
c
c     Calculate Larmor radius - use toroidal field for now.
c
      b_field = bts(ik_local,ir_local)
c
      larmor_radius = larmor(massi,temi,b_field,riz)
c
c     Find distance to target from ionization position to linear
c     extension of target element.
c
      targ_dist = dist_to_point(r,z,rp(id),zp(id),thetas(id))
c slmod begin
      IF (getrz_error) targ_dist = 1.0E+20
c slmod end
c
c     Does Prompt depostion occur?
c
c     Use simple initial assumptions -
c
c     1) IF dist_to_point <= MPS then prompt redeposition
c        at reduced impact energy.
c     2) IF dist_to_point <= larmor_radius then redeposition
c        at full energy.
c
c
c      write(6,*) 'DEBUG PD:',ik,ir,targ_dist,
c     >           mps_thickness(ir_local,it),larmor_radius,bfield
c
c
      if (targ_dist.le.mps_thickness(ir_local,it)) then
c
c        Redeposition due to MPS effect
c
         sheath_drop = targ_dist/mps_thickness(ir_local,it)
     >                     * mps_energy(ir,it)
     >                     + (3.0-mps_energy(ir_local,it))
c
c        set return code - prompt deposition has occurred
c
c
c        Record prompt deposition
c
         promptdeps(id,1) = promptdeps(id,1) + sputy
         promptdeps(id,2) = promptdeps(id,2) + targ_dist*sputy
         promptdeps(id,3) = promptdeps(id,3) + larmor_radius*sputy
         promptdeps(id,4) = promptdeps(id,4) + sheath_drop*sputy
c
         rc = 1
c
      elseif (targ_dist.le.larmor_radius) then
c
c        Redeposition due to Larmor Radius effect
c
c        Inside Larmor_radius but outside MPS
c
         sheath_drop = 3.0
c
c
c        Record prompt deposition
c
         promptdeps(id,1) = promptdeps(id,1) + sputy
         promptdeps(id,2) = promptdeps(id,2) + targ_dist*sputy
         promptdeps(id,3) = promptdeps(id,3) + larmor_radius*sputy
         promptdeps(id,4) = promptdeps(id,4) + sheath_drop*sputy
c
c        set return code - prompt deposition has occurred
c
         rc = 1
c
      endif
C
      return
      end
c
c
c
      real function dist_to_point(rp,zp,r1,z1,theta)
      implicit none
      real rp,zp,r1,z1,theta
c
c     DIST_TO_POINT:
c
c     This routine calculates the perpendicular distance from the
c     point R,Z to the linear extension of the element specified by
c     a point on the element and the normal angle THETA to the target surface.
c
c     Theta is defined as the angle normal to each target
c     segment measured from the positive R-axis.
c
c     The line of the target is defined by the equation ar+bz=c
c
c     The normal to this line is the vector (a,b) = (cos(theta),sin(theta))
c
c     The perpendicular distance D is equal to the magnitude of the vector
c     from (r1,z1) to (rp,zp) projected onto the normal to the line.
c
c     B = p1 -> pp       D = | B cos (alpha) |
c     alpha = angle between B and N (normal to the line)
c
c     But B .dot. N =  | B | | N| cos (alpha)
c
c     Therefore   D = | B .dot. N|  / |N|
c
c     B = <rp-r1, zp-z1>   N = <cos(theta),sin(theta)>
c
c     D = |cos(theta) * (rp-z1) + sin(theta) * (zp -z1)| / |N|
c
c     |N| = 1
c
c     This is the formula used here in this function.
c
c
      dist_to_point = abs(cos(theta)*(rp-r1) + sin(theta)*(zp-z1))
c
      return
      end
c
c
c
      subroutine update_walldep(ik,ir,iz,idt,idw,iwstart,idtype,sputy,
     >                          eimp)
      implicit none
c
      integer ik,ir,iz,iwstart,idtype
      integer,intent(in) ::  idt,idw
      real sputy,eimp
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'dynam3'
c
c     UPDATE_WALLDEP:
c
c     This routine records the ion impact with the wall segment
c     closest to where the particle left the grid. If ID is not
c     equal to zero this is treated as a specified target index
c     and is mapped to the equivalent wall index through the
c     wallindex array. If it is zero - the R,Z of the cell centre
c     specified by the ik,ir indices is used to find the shortest
c     distance to the nearest wall segment centre - this is then
c     used for recording the exit wall segment for the ion.
c
c     David Elder     Nov 5, 1998
c
c     Added a wallsiz array that records the wall impact information in a
c     charge resolved maner
c
c     K. Schmid Feb. 2008 and june 2009
      real best,dsq,r,z
      integer ind,id
c
c     If neither a target nor wall segment has been specified find the wall segment centre
c     closest to the cell centre of the particle exit.
c
      if (cprint.eq.9) then 

         write(6,'(a,7i6,5(1x,g12.5))') 'Update_walldep:',ik,ir,iz,
     >         idt,idw,
     >         iwstart,idtype,sputy,eimp

      endif
c
c     Left grid through periphery ... add Teb factor to energy when element is known
c
      if (idt.eq.0.and.idw.eq.0) then
c
          R = RS(ik,ir)
          Z = ZS(ik,ir)
c
          BEST = HI
          DSQ  = HI
          IND = 1
          DO 650 ID = 1,WALLPTS
             DSQ = (WALLPT(ID,1)-R) ** 2 + (WALLPT(ID,2)-Z) ** 2
             IF (DSQ.LT.BEST) THEN
               BEST = DSQ
               IND   = ID
             ENDIF
650       CONTINUE
C
c          WRITE(6,'(a,2i5,5(1x,g13.6))')
c     >              'DSQ:',ind,WALLPTS,DSQ,R,Z,ROLD,ZOLD
C

          if (ind.lt.1.or.ind.gt.wallpts) then
c
             if (cprint.eq.9) then 
                write(6,'(a,4i10)') 'Wallsi: No wall found:',
     >                             idt,ind,ik,ir
             endif
c
             wallsi(maxpts+1) = wallsi(maxpts+1) + sputy
c
             wallsiz(maxpts+1,iz) = wallsiz(maxpts+1,iz) + sputy
             wallseiz(maxpts+1,iz)= wallseiz(maxpts+1,iz) + eimp * sputy
c
             if (iwstart.ge.1.and.iwstart.le.wallpts) then

                wtdep(iwstart,maxpts+1,1) =
     >                    wtdep(iwstart,maxpts+1,1) + sputy
             endif

          else

             wallsi(ind) = wallsi(ind) + sputy
             wallsiz(ind,iz)  = wallsiz(ind,iz) + sputy

c
c            Add plasma temperature at deposition element to impact energy
c
             eimp = eimp + 3.0 * iz * wallpt(ind,29) 

             wallseiz(ind,iz) = wallseiz(ind,iz) + eimp * sputy

c
c             write(6,*) 'ind case ind:',ind,' iz: ',iz,' sputy: ', sputy
c

             if (iwstart.ge.1.and.iwstart.le.wallpts) then

                wtdep(iwstart,ind,1) =
     >                    wtdep(iwstart,ind,1) + sputy
             endif

          endif

c
c      Target segment specified
c
       elseif (idt.gt.0) then
c
          if (wallindex(idt).ne.0) then

             wallsi(wallindex(idt)) = wallsi(wallindex(idt))+ sputy
             wallsiz(wallindex(idt), iz) =
     >            wallsiz(wallindex(idt), iz) + sputy
             wallseiz(wallindex(idt), iz) =
     >            wallseiz(wallindex(idt), iz) + eimp * sputy


c             write(6,*) 'idt case ind:',wallindex(idt),' iz: ',iz,
c     >          ' sputy: ', sputy,
c     >          ' wallsiz:', wallsiz(wallindex(idt), iz)

             if (iwstart.ge.1.and.iwstart.le.wallpts) then

                wtdep(iwstart,wallindex(idt),1) =
     >                    wtdep(iwstart,wallindex(idt),1) + sputy
             endif

          else

c             write (6,'(a,3i5)') 'Wallsi: target?:',idt,wallindex(idt)
c                   wallsi(maxpts+1) = wallsi(maxpts+1) + sputy
c                   wallsiz(maxpts+1, iz) = wallsiz(maxpts+1, iz) + sputy
c
             if (iwstart.ge.1.and.iwstart.le.wallpts) then

                wtdep(iwstart,maxpts+1,1) =
     >                    wtdep(iwstart,maxpts+1,1) + sputy
             endif

          endif
c
c      Wall segment specified
c
       elseif (idw.ge.1.and.idw.le.wallpts) then
c
             wallsi(idw) = wallsi(idw)+ sputy
             wallsiz(idw, iz) = wallsiz(idw, iz) + sputy
             wallseiz(idw, iz) = wallseiz(idw, iz) + eimp * sputy

c
c         write(6,*) 'idw case ind:',idw,' iz: ',iz,' sputy: ', sputy
c

             if (iwstart.ge.1.and.iwstart.le.wallpts) then

                wtdep(iwstart,idw,1) =
     >                    wtdep(iwstart,idw,1) + sputy
             endif

       else

c          write (6,'(a,3i5)') 'Wallsi: wall?:',idw

          wallsi(maxpts+1) = wallsi(maxpts+1) + sputy
          wallsiz(maxpts+1, iz) = wallsiz(maxpts+1, iz) + sputy
          wallseiz(maxpts+1, iz) = wallseiz(maxpts+1, iz) + eimp * sputy
c
          if (iwstart.ge.1.and.iwstart.le.wallpts) then

             wtdep(iwstart,maxpts+1,1) =
     >                    wtdep(iwstart,maxpts+1,1) + sputy
          endif

       endif
c
c
c
       return
       end
c
c
c
       integer function verify_id(ik,ir,itarg)
       implicit none
       integer ik,ir,itarg
c
       include 'params'
       include 'cgeom'
c
c      VERIFY_ID: This routine looks at the target segment
c                 associated with the given ik,ir indices and
c                 then returns the closest valid target
c                 segment id. The reason for this is to catch
c                 cases where a particle strikes the wall and
c                 target at the same time - it must be assured of
c                 a relaunch from a valid target index.
c
c                 David Elder,     Nov 11, 1998
c
c
       integer id
c
       id = idds(ir,itarg)
c
       do while (dds(id).eq.0.0)
c
          if (ir.le.irwall) then
c
             ir = irins(ik,ir)
             id = idds(ir,itarg)

          else

             ir = irouts(ik,ir)
             id = idds(ir,itarg)

          endif
c
       end do
c
       verify_id = id
c
       return
       end
c
c
c
      subroutine get_random_numbers(kk,kklim,nrand,seed)
      implicit none
      integer kk,nrand,kklim
      real*8 seed
      include 'params'
      include 'crand'

      IF (KK.GT.KKLIM) THEN
         CALL SURAND (SEED, KK, RANV)
         NRAND = NRAND + KK
         KK = 0
      ENDIF


      return
      end
c
c
c
      subroutine debug_velocity
      implicit none
      include    'params'
      include    'dynam3'
      include    'comtor'
      include    'cgeom'
      include    'diagvel'
c
      include 'div1'
      include 'div2'
      include 'div5'
      include 'div6'
c
      include    'particle_specs'
c
        if (debugv) then
c
           sdvs (ik,ir,iz)  = sdvs(ik,ir,iz)  + sputy * fvel
           sdvs2(ik,ir,iz)  = sdvs2(ik,ir,iz) + sputy * fvel**2.0
c
           if (abs(fvel).gt.sdvb(ik,ir)) then
c
              sdvs3(ik,ir,iz,1) = sdvs3(ik,ir,iz,1) + sputy * abs(fvel)
              sdvs3(ik,ir,iz,2) = sdvs3(ik,ir,iz,2) + sputy
c
           endif
c
c          Split up the velocity distribution and assign it
c          to an appropriate bin.
c
           in = int(fvel/(velsep*velplate))
           if (iabs(in).gt.nvel) in = isign(nvel,in)
           if (in.ge.0.and.fvel.ge.0.0) then
              in = in + 1
           endif
c
c          Set knot currently occupied.
c
           if (maxvnks.gt.nks(injir).and.nks(injir).gt.0) then
              ikv = ik
           else
              ikv = 1
           endif
c
           velspace(in,iz,ikv) = velspace(in,iz,ikv) + fvel * sputy
           velweight(in,iz,ikv) = velweight(in,iz,ikv) + sputy
        endif
c
      return
      end
c
c
c
      subroutine check_ion_change_state(seed,nrand,neutim,nizs)
      implicit none
      real    neutim
      real*8  seed
      integer nrand,nizs
      include    'params'
      include    'dynam3'
      include    'dynam4'
      include    'comtor'
      include    'cgeom'
      include    'cioniz'
      include    'commv'
      include    'cneut'
      include    'crand'
c
      include 'div1'
      include 'div3'
      include 'div4'
      include 'div5'
      include 'div6'
c
      include    'particle_specs'

c
c     Output velocity along the field line from launch_one
c
      real vout


      integer ipos
      external ipos


C
C-----------------------------------------------------------------------
C       IONISATION AND RECOMBINATION
C-----------------------------------------------------------------------
C
        KK = KK + 1
        IF (RANV(KK).LT.KPCHS(IK,IR,IZ)) THEN
          KK = KK + 1
          IF (RANV(KK).LT.KPRCS(IK,IR,IZ)) THEN
            CICRCS(IZ) = CICRCS(IZ) + SPUTY
c            CIFRCS(IZ) = MIN (CIFRCS(IZ), CIST)
c            CILRCS(IZ) = MAX (CILRCS(IZ), CIST)
            CIFRCS(IZ) = MIN (CIFRCS(IZ), sngl(CIST))
            CILRCS(IZ) = MAX (CILRCS(IZ), sngl(CIST))
            CISRCS(IZ) = CISRCS(IZ) + CIST * SPUTY
c
c           Recording addition to total elapsed time in state IZ
c
            cieizs(iz) = cieizs(iz) + cistiz * sputy
            citizs(iz) = citizs(iz) + sputy
            cistiz = 1.0
c
            IZ = IZ - 1
            IT = IPOS (sngl(CIST), CTIMES(1,IZ), NTS+1)
            RIZ = REAL(IZ)
c
c           Recombination to neutral
c
            IF (IZ.LT.1) THEN
              TBELOW = TBELOW + SPUTY
c
c              IFATE = 5
c
c             Record recombination statistics
c
              if (cfolrec.eq.0) then
                 reccnt = reccnt + 1
c
c                Set IFATE to 5 for recombined impurity
c
                 ifate = 5
c
c                 GOTO 790
c
              elseif (cfolrec.eq.1) then
c
c                Follow recombined impurities
c
                 cist_recstart = cist
c
                 ikorg = ik
                 irorg = ir
c
                 call LAUNCH_ONE (IMP,R,Z,RIZPOS,ZIZPOS,0,iwstart,
     >                   rc,temi,cist,sputy,
     >                   recSTRUK,mtcrecstruk,recMAIN,recEXIT,
     >                   recATIZ,recNEUT,recWALLN,mtcrecwalln,
     >                   recCENT,recTMAX,
     >                   SEED,NRAND,
     >                   NEUTIM,recFAIL,6,vout,vrec,recloss)
c
c
c               Set vel to appropriate value by scaling the along the field
c               line vout value returned by launch_one
c
                vel = vout * qtim
c
                 cist_elapsed = cist - cist_recstart
                 dist_travelled = 2.0 * abs(vel) *
     >                            cist_elapsed * qtim
c
                 recinf(1,rc) = recinf(1,rc) + sputy
                 reccnt = reccnt + 1
c
                 if (reccnt.eq.1) then
                    recinf(14,rc) = recinf(14,rc) + sputy
                    recinf(2,rc) = recinf(2,rc) + cist_elapsed * sputy
                    recinf(3,rc) = recinf(3,rc) + dist_travelled* sputy
                 endif
c
                 recinf(4,rc) = recinf(4,rc) + cist_elapsed * sputy
                 recinf(5,rc) = recinf(5,rc) + dist_travelled * sputy
                 recinf(6,rc) = recinf(6,rc) + temi * sputy
                 recinf(7,rc) = recinf(7,rc) + sputy * abs(vrec)
c
c                Other quantities
c
                 recinf(8,rc) = recinf(8,rc) + knbs(ikorg,irorg) * sputy
                 recinf(9,rc) = recinf(9,rc) + knhs(ikorg,irorg) * sputy
                 recinf(10,rc) = recinf(10,rc) +
     >                           ktebs(ikorg,irorg) * sputy
                 recinf(11,rc) = recinf(11,rc) + r  * sputy
                 recinf(12,rc) = recinf(12,rc) + z  * sputy
                 recinf(13,rc) = recinf(13,rc) + min(s,smax-s) * sputy
c
c                 write (6,'(a,2i8,6g12.5)') 'RECOM IMP:',imp,rc,
c     >                       r,z,rizpos,zizpos,
c     >                       temi,cist
c
c
c                For all results other than reionization to state 1 - the code
c                will stop following the particle at this point and exit
c                as it would have for the old recombination implementation.
c
                 if (rc.ne.5) then
c
c                   Set IFATE to 5 - though it is not lost directly due to recombination
c                   the particle is removed as part of the recombination process.
c
                    ifate = 5
c
c                   Exit from this routine
c
c                    goto 790
c
                 else
c
c                   Deal with particle that was re-ionized
c
                    r = rizpos
                    z = zizpos
                    iz= 1
c
c                   Look near last recorded ik,ir
c
                    call gridpos(ik,ir,r,z,.false.,griderr)
c
c                   If a griderr then exit anyway and issue error message
c                   since this shouldn't happen.
c
                    if (griderr) then

                       write (6,*) 'RECOMBINATION ERROR: NOT ON GRID:',r,z
c
c                      Set IFATE to 5 - lost via recombination
c
                       ifate = 5
c
c                       goto 790
c
                    else
c
c                      If on grid - re-assign S-value - to cell centre
c                      if it left the cell - otherwise leave as is.
c
                       if (ik.ne.ikorg.or.ir.ne.irorg) then
c
c                         Also need to reset the SMAX value for the new ring
c
                          smax = ksmaxs(ir)
c
                          if (init_pos_opt.eq.0) then
                             s = kss(ik,ir)
                             cross = 0.0
                          elseif (init_pos_opt.eq.1) then
c
                             call getscross_approx(r,z,s,cross,ik,ir)
c
                          endif
c
                       endif
c
                    endif


                 endif
c
              endif
c
            ENDIF
C
          ELSE
            CICIZS(IZ) = CICIZS(IZ) + SPUTY
c            CIFIZS(IZ) = MIN (CIFIZS(IZ), CIST)
c            CILIZS(IZ) = MAX (CILIZS(IZ), CIST)
            CIFIZS(IZ) = MIN (CIFIZS(IZ), sngl(CIST))
            CILIZS(IZ) = MAX (CILIZS(IZ), sngl(CIST))
            CISIZS(IZ) = CISIZS(IZ) + CIST * SPUTY
c
c           Recording addition to total elapsed time in state IZ
c
            cieizs(iz) = cieizs(iz) + cistiz * sputy
            citizs(iz) = citizs(iz) + sputy
            cistiz = 1.0
c
            TIZS(IK,IR,IZ) = TIZS(IK,IR,IZ) + SPUTY
            IZ = IZ + 1
            IT = IPOS (sngl(CIST), CTIMES(1,IZ), NTS+1)
            RIZ = REAL(IZ)
            IF (IZ.EQ.CIZSET) TEMI = MAX (TEMI,KTIBS(IK,IR))
            IF (IZ.GT.NIZS) THEN
              TBYOND = TBYOND + SPUTY
              IFATE = 6
c
c             IFATE is now used to define the exit condition
c
c              GOTO 790
c
            ENDIF
            MAXCIZ = MAX (MAXCIZ, IZ)
          ENDIF
        ENDIF


      return

      end
c
c
c
      subroutine check_ion_removal(ifate,kk,ik,ir,iz,cist,cistiz,
     >                             ssss,s,smax,sputy)
      implicit none
c
      include    'params'
      include    'commv'
      include    'cioniz'
      include    'crand'
c
      integer ifate,kk,ik,ir,iz
      real ssss(maxizs),s,sputy,smax
      real*8 cist,cistiz
C
C-----------------------------------------------------------------------
C       ION REMOVAL
C-----------------------------------------------------------------------

        KK = KK + 1
        IF( RANV(KK).LT.KPLOS(IK,IR,IZ) ) THEN
            IFATE = 8
            CICLOS(IZ) = CICLOS(IZ) + SPUTY
c            CIFLOS(IZ) = MIN( CIFLOS(IZ) , CIST )
c            CILLOS(IZ) = MAX( CILLOS(IZ) , CIST )
            CIFLOS(IZ) = MIN( CIFLOS(IZ) , sngl(CIST) )
            CILLOS(IZ) = MAX( CILLOS(IZ) , sngl(CIST) )
            CISLOS(IZ) = CISLOS(IZ) + CIST*SPUTY
c
c           Recording addition to total elapsed time in state IZ
c
            cieizs(iz) = cieizs(iz) + cistiz * sputy
            citizs(iz) = citizs(iz) + sputy
c
            SSSS(IZ)   = MAX( SSSS(IZ) , MIN( S , SMAX-S ) )
c
c           goto 790
c
        END IF

      return
      end

c
c
c
      subroutine change_local_values
      implicit none
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'cioniz'
      include    'clocal'
c
      include 'div1'
      include 'div2'
c
      include    'particle_specs'

c
C
C-----------------------------------------------------------------------
C       PARALLEL DIFFUSION AND ION TEMPERATURE
C-----------------------------------------------------------------------
C
C------ QUICKLY ADJUST RELEVANT CHARACTERISTIC TIMES
C------ FOR NON-STANDARD PLASMA OPTIONS.  BECAUSE OF THE SQUARE ROOTS
C------ THESE CALCULATIONS ARE DONE SPARINGLY TO PREVENT DRAMATIC
C------ SLOWING OF PROGRAM EXECUTION - FOR EXAMPLE, WHENEVER THE
C------ TEMPERATURE CHANGES BY 20% OR MORE.
C
        IF( ((CIOPTB.eq.2.OR.cioptb.eq.3.or.cioptb.eq.4.or.
     >       cioptb.eq.7.or.cioptb.eq.9.or.
     >       CIOPTC.eq.2.or.cioptc.eq.3)
     >       .AND.IR.GE.IRSPEC)
     >       .OR.
     >       CIOPTD.GE.3.or.cioptb.eq.5.or.cioptb.eq.10) THEN
            TEMOLD = LTOLDS(IK,IR,IZ)
            IF (TEMI.GT.1.2*TEMOLD .OR. TEMI.LT.0.8*TEMOLD) THEN
              LTOLDS(IK,IR,IZ) = TEMI
C
              IF (CIOPTB.EQ.2) THEN
                RATIO1 = SQRT (TEMOLD / TEMI)
                RATIO2 = 1.0 / SQRT (RATIO1)
                LFPS(IK,IR,IZ) = RATIO1 * LFPS(IK,IR,IZ)
                LLLFPS(IK,IR,IZ) = RATIO2 * LLLFPS(IK,IR,IZ)
              ELSEIF (CIOPTB.EQ.3.or.cioptb.eq.7.or.cioptb.eq.9
     >                .or.cioptb.eq.5.or.cioptb.eq.10) THEN
                RATIO2 = SQRT (TEMI / TEMOLD)
                LLLFPS(IK,IR,IZ) = RATIO2 * LLLFPS(IK,IR,IZ)
              ELSEIF (CIOPTB.EQ.4) THEN
                IF (TEMI.LE.KTIBS(IK,IR)*CRMI/CRMB) THEN
                  LFPS(IK,IR,IZ) = KFPS(IK,IR,IZ)
                  LLLFPS(IK,IR,IZ) = KKKFPS(IK,IR,IZ)
                ELSE
                  LFPS(IK,IR,IZ) = C350B*REAL(IZ*IZ)*KTIBS(IK,IR) *
     >              KNBS(IK,IR) / (TEMI ** 1.5)
                  LLLFPS(IK,IR,IZ) = QTIM *
     >              SQRT(4.88E8/(LFPS(IK,IR,IZ)*CRMI))
                ENDIF
              ENDIF
C
              IF (CIOPTC.EQ.2) THEN
                TAU = LFPS(IK,IR,IZ) / (2.0*TEMI)
                IF (TAU.GT.1.E-3) THEN
                  LFSS(IK,IR,IZ) = 1.0 - EXP(-TAU)
                ELSE
                  LFSS(IK,IR,IZ) = TAU
                ENDIF
              ELSEIF (CIOPTC.EQ.3) THEN
                IF (TEMI.LE.KTIBS(IK,IR)*CRMI/CRMB) THEN
                  LFSS(IK,IR,IZ) = KFSS(IK,IR,IZ)
                ELSE
                  TAU = C350A * REAL(IZ*IZ) * KNBS(IK,IR) /
     >                  (TEMI ** 1.5)
                  IF (TAU.GT.1.E-3) THEN
                    LFSS(IK,IR,IZ) = 1.0 - EXP(-TAU)
                  ELSE
                    LFSS(IK,IR,IZ) = TAU
                  ENDIF
                ENDIF
              ENDIF
C
              IF (CIOPTD.EQ.3) THEN
                TAU = C215A * REAL(IZ*IZ) * KNBS(IK,IR) /
     >                ((CRMI * KTIBS(IK,IR) + CRMB * TEMI) ** 1.5)
                IF (TAU.GT.1.E-3) THEN
                  LFTS(IK,IR,IZ) = 1.0 - EXP (-TAU)
                ELSE
                  LFTS(IK,IR,IZ) = TAU
                ENDIF
              ENDIF
            ENDIF
        ENDIF

      return
      end
c
c
c
      subroutine setup_drftv
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'driftvel'
      include 'fperiph_com'
c
c     SETUP_DRFTV:
c
c     Setup the drift velocity for each ring
c
c     Local variables
c
      real startp,endp
      integer irstart,irend,ir,in,ikmid,ik
      real calc_cs
      external calc_cs
c
c     Iniitalize to zero in all cases - this is all that is
c     done for option 0
c
      call rzero(pol_drftv,maxnrs)
c
c     Initialize start and end regions for each ring such that
c     a drift is not possible unless the values are set correctly
c
      do ir = 1,nrs
         sdrft_start(ir) = ksmaxs(ir)
         sdrft_end(ir) = 0.0
      end do

c
c      write(0,*) 'CPDRFT=',cpdrft
c
c     Exit this routine if the poloidal drift option is not active
c

      if (cpdrft.eq.0) return

c
c
c     Assign the flow velocity if the option is ON
c
c     Note: The following IF can be eliminated though it allows
c     for different behaviours to be supported in future options.
c
      IF (CPDRFT.EQ.1.or.cpdrft.eq.2.or.cpdrft.eq.3) THEN
c
c     Assign/calculate the drift velocity for each flux tube
c
         call get_drftv_rings(irstart,irend)
c
c     If no per ring data is entered assign the default value - old behaviour
c     Initialize to base value
c
c     The base value may be interpreted as a MACH number if the mach option is
c     set.
c
c
         if (drftvel_machopt.eq.0) then
c
            do ir = irstart,irend
               pol_drftv(ir) = CDRFTV * QTIM
            end do
c
         elseif (drftvel_machopt.eq.1.or.drftvel_machopt.eq.2) then
c
            call rzero(ringcs,maxnrs)
c
c     Opt 1 : 2 x TE based sound speed
c     Opt 2 : TE + TI based sound speed
c
            do ir = irstart,irend
               ikmid = ikmids(ir)
               ringcs(ir) =
     >              (calc_cs(ktebs(ikmid,ir),ktibs(ikmid,ir),crmb,
     >              drftvel_machopt)+
     >              calc_cs(ktebs(ikmid+1,ir),ktibs(ikmid+1,ir),crmb,
     >              drftvel_machopt))
     >              /2.0
c
c     Assign default drift velocity
c
               pol_drftv(ir) = CDRFTV * QTIM * ringcs(ir)
c
            enddo
c
         endif
c
c     Overwrite the default value with per ring information
c
         if (ndrftvel.gt.0) then
c
c     If the mach number option is ON - calculate the sound speed on
c     each ring at the mid-point - averaged over the two cells on
c     either side of the midpoint.
c     For MACH option 1 use 2 x TE,
c     For MACH option 2 use TE+TI
c
c     Assign input velocity directly for each ring
c
            if (drftvel_machopt.eq.0) then
c
c     Now assign per ring velocities
c
c     Assign per ring data
c
               do in = 1,ndrftvel
                  ir = int(ringdrftvel(in,1))
                  pol_drftv(ir) = ringdrftvel(in,2) * qtim
               enddo

c
c     Opt 1 : 2 x TE based sound speed
c     Opt 2 : TE + TI based sound speed
c
            elseif (drftvel_machopt.eq.1.or.drftvel_machopt.eq.2) then
c
c
c     Now assign detiled per ring velocities
c
               do in = 1,ndrftvel

                  ir = int(ringdrftvel(in,1))
c
                  pol_drftv(ir) = ringdrftvel(in,2) * ringcs(ir) * qtim
c
               enddo

            endif

         endif

c
c         write(0,*) 'Drift Dist Opt=',drft_distopt,
c     >                 cdrftv_start,cdrftv_end
c
c     Calculate the drift regions for each affected ring in terms of S
c
         do ir = irstart,irend

c
c     Distance specified as fraction of S parallel
c
            if (drft_distopt.eq.0) then
               sdrft_start(ir) = cdrftv_start * ksmaxs(ir)
               sdrft_end(ir)   = cdrftv_end * ksmaxs(ir)
c
c     Distance specified as a fraction of the poloidal distance along the field lines
c
            elseif (drft_distopt.eq.1) then

               startp = kpmaxs(ir) * cdrftv_start
               endp   = kpmaxs(ir) * cdrftv_end

               do ik = 1,nks(ir)
                                ! found startp bin
                  if ((kpb(ik-1,ir).le.startp).and.
     >                 (kpb(ik,ir).ge.startp)) then
c
                     sdrft_start(ir) = ksb(ik-1,ir) +
     >                    (ksb(ik,ir)-ksb(ik-1,ir)) *
     >                    (startp-kpb(ik-1,ir))/
     >                    (kpb(ik,ir)-kpb(ik-1,ir))
c
                  elseif ((kpb(ik-1,ir).le.endp).and.
     >                    (kpb(ik,ir).ge.endp)) then
c
                     sdrft_end(ir) = ksb(ik-1,ir) +
     >                    (ksb(ik,ir)-ksb(ik-1,ir)) *
     >                    (endp-kpb(ik-1,ir))/(kpb(ik,ir)-kpb(ik-1,ir))
c
                  endif

               end do

c
c     Distance specified relative to a specified Z value
c
            elseif (drft_distopt.eq.2) then
c
c     This option must be run on two half rings since the Z
c     values repeat (meaning that Z by itself can not
c     be distinguished between inner and outer the way the
c     S and P coordinates can be.
c
c     In this case cdrftv_start and cdrftv_end will contain
c     separate Z coordinates to apply to the first and second
c     halves of the ring
c
c
c     Need to first check if the Z value specified is
c     beyond the end of the flux tube
c
               if (abs(kzb(0,ir)-z0).lt.
     >             abs(cdrftv_start-z0)) then

                  sdrft_start = 0.0

               else


                  do ik = 1,nks(ir)/2


                     if (((kzb(ik-1,ir).le.cdrftv_start).and.
     >                    (kzb(ik,ir).ge.cdrftv_start)).or.
     >                    ((kzb(ik-1,ir).ge.cdrftv_start).and.
     >                    (kzb(ik,ir).le.cdrftv_start)) ) then
c
                        sdrft_start(ir) = ksb(ik-1,ir) +
     >                       (ksb(ik,ir)-ksb(ik-1,ir)) *
     >                       (cdrftv_start-kzb(ik-1,ir))/
     >                       (kzb(ik,ir)-kzb(ik-1,ir))
c
                     endif
                  end do

               endif

c     Need to first check if the Z value specified is
c     beyond the end of the flux tube

               if (abs(kzb(nks(ir),ir)-z0).lt.
     >             abs(cdrftv_end-z0)) then

                  sdrft_end = ksmaxs(ir)

               else


                  do ik = nks(ir),nks(ir)/2, -1
                     if (((kzb(ik-1,ir).le.cdrftv_end).and.
     >                    (kzb(ik,ir).ge.cdrftv_end)).or.
     >                    ((kzb(ik-1,ir).ge.cdrftv_end).and.
     >                    (kzb(ik,ir).le.cdrftv_end)) ) then
c
                        sdrft_end(ir) = ksb(ik-1,ir) +
     >                       (ksb(ik,ir)-ksb(ik-1,ir)) *
     >                       (cdrftv_end-kzb(ik-1,ir))/
     >                       (kzb(ik,ir)-kzb(ik-1,ir))
c
                     endif
                  end do

               endif

            endif

            write(6,'(a,i6,1(1x,g12.5),4(1x,f16.5))')
     >           'Poloidal drift:',ir,pol_drftv(ir)/qtim,
     >           sdrft_start(ir),sdrft_end(ir),ksmaxs(ir),kpmaxs(ir)
c            write(6,'(a,4(1x,f16.5))') 'Input:', cdrftv_start,
c     >                 cdrftv_end,startp,endp

            do ik = 1,nks(ir)
               write(6,'(a,2i6,10(1x,f16.5))') 'Poloidal data:',
     >             ik,ir,ksb(ik-1,ir),kss(ik,ir),ksb(ik,ir),
     >               kpb(ik-1,ir),kps(ik,ir),kpb(ik,ir),
     >               kzb(ik-1,ir),zs(ik,ir),kzb(ik,ir)
            end do

c            write(0,'(a,i6,1(1x,g12.5),4(1x,f16.5))')
c     >           'Poloidal drift:',ir,pol_drftv(ir)/qtim,
c     >           sdrft_start(ir),sdrft_end(ir),ksmaxs(ir),kpmaxs(ir)

         end do




      ENDIF


c
c     If far periphery transport is active
c
      if (fpopt.eq.5) then
c
c     Set up far periphery flow velocity
c
c     NOTE!: fp_flow_velocity is an array with an element
c     for each peripheral region
c
c
c        Set up flow velocity in each peripheral region
c
         if (fp_flow_opt.eq.0) then
            fp_flow_velocity = 0.0
         elseif (fp_flow_opt.eq.1) then
            do in = 1,num_fp_regions
               fp_flow_velocity(in) = pol_drftv(fp_rings(in))
            end do
         elseif (fp_flow_opt.eq.2) then
            fp_flow_velocity = fp_flow_velocity_input * qtim
         elseif (fp_flow_opt.eq.3) then
            do in = 1,num_fp_regions
               fp_flow_velocity(in) = pol_drftv(fp_virt_rings(in))
            end do
         endif


      endif

      return
      end
c
c
c
      real function calc_cs(te,ti,mb,opt)
      implicit none
      real te,ti,mb
      integer opt
c
c     Local variables
c
      real temp
c
c     Calculate the sound speed
c
      if (opt.eq.1) then
         temp = 2.0 * te
      elseif (opt.eq.2) then
         temp = te * ti
      endif

      calc_cs = 9.79E3 * SQRT (temp/mb)

      return
      end
c
c
c
      subroutine get_drftv_rings(irstart,irend)
      implicit none
      integer irstart,irend
      include 'params'
      include 'cgeom'
      include 'driftvel'
c
c       Option 3 applies only to the PFZ
c
c       Added drift region option - preferred procedure to using the CPDRFT method
c       to limit flow regions.
c
c       Region 1 - SOL + PFZ
c
        if (drft_region.eq.1) then
           irstart = irsep
           irend   = nrs
           if (cpdrft.eq.3) irstart = irtrap
c
c       Region 2 - SOL only
c
        elseif (drft_region.eq.2) then
           irstart = irsep
           irend   = irwall
c
c       Region 3 - PFZ only
c
        elseif (drft_region.eq.3) then
           irstart = irtrap
           irend   = nrs
c
c       Region 4 - CORE only
c
        elseif (drft_region.eq.4) then
           irstart = 2
           irend   = irsep-1
        endif

      return
      end
c
c
c
      subroutine get_drftv(pol_vel,start_s,end_s,ir)
      implicit none
c
      integer ir
      real pol_vel,start_s,end_s
c
c     Return the drift velocity characteristics - value and range
c     for the specified ring
c
      include 'params'
      include 'driftvel'
c
      pol_vel = pol_drftv(ir)
      start_s = sdrft_start(ir)
      end_s   = sdrft_end(ir)
c
      return
      end
c
c
c
      subroutine set_normalization_factors(facta,factb,tatiz,tneut,
     >                               fsrate,qtim,nizs,cneuta,normopt)
      implicit none
      include 'params'
      REAL     FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
      real     tatiz,tneut
      real     fsrate,qtim
      integer  cneuta,nizs
      integer  normopt

C
C     SET THE DEFAULT VALUE FOR THESE SCALING CONSTANTS TO 0.0
C     THEN IF THE NUMBER OF NEUTRALS LAUNCHED IS ZERO - AS OCCURS
C     FOR SOME ION INJECTION CASES - ALLOW THE SCALING FACTOR TO
C     BE SET TO THE INVERSE OF THE NUMBER OF IONS INJECTED.
C
c
c     Changed this normalization - ddlims(,,-1) should also be
c     normalized to total neutrals entering/second so that the
c     density is properly comparable to the ddlims(,,0) entry
c     containing the total neutral density.
c
c     jdemod - there is a problem with normalization of the density
c     data for neutrals. It is scaled by 1/TNEUT while the ions are
c     scaled by 1/TATIZ - however, ABSFAC is multiplied by TATIZ/TNEUT.
c     This results in the neutral densities being incorrect when they
c     are used in OUT if there are multiplied by ABSFAC.
c
c     To correct this situation - either the neutral data could be
c     scaled to one ion entering the system/sec - which doesn't make
c     a lot of sense logically - but then ABSFAC would be left as is.
c     On the other hand, all of the data could be scaled by the number
c     of particles launched for the simulation - ions for ion injection
c     and neutrals for a neutral launch. This would then put consistent
c     densities in all of the arrays - and the ABSFAC value would be
c     the total absolute particle influx - though this is usually neutrals.
c
c      IF (RNEUT1.GT.0.0) THEN
c         FACTA(-1) = 1.0 / RNEUT1
c
c
c     The original normalization code is included as a hard coded option
c     in case I want to replicate old cases in the future
c

      real particle_norm_factor
      integer iz
c
c     Modified normalization scalings (default)
c
c     Set normalization to zero initially
c
      facta = 0.0
      factb = 0.0
c
      if (normopt.eq.0) then
c
c        Neutral launch - cneuta = 0
c
         IF (cneuta.eq.0.and.tneut.gt.0) THEN
            particle_norm_factor = 1.0/TNEUT
         ELSEIF (CNEUTA.eq.1.and.tatiz.gt.0) THEN
            particle_norm_factor = 1.0/TATIZ
         else
            particle_norm_factor = 1.0
         ENDIF
C
         IF (NIZS.GT.0) THEN
           DO IZ = -1, NIZS
             facta(iz) = particle_norm_factor
             if (iz.lt.1) then
                factb(iz) = facta(iz) * fsrate
             else
                factb(iz) = facta(iz) * qtim
             endif
           enddo
         ENDIF
c
      elseif (normopt.eq.1) then
c
c        Original normalization code
c
         IF (TNEUT.GT.0.0) THEN
            FACTA(-1) = 1.0 / TNEUT
         ELSEIF (TATIZ.GT.0.0) THEN
            FACTA(-1) = 1.0 / TATIZ
         ENDIF
         FACTB(-1) = FACTA(-1) * FSRATE
C
         FACTA(0) = 0.0
         IF (TNEUT.GT.0.0) THEN
            FACTA(0) = 1.0 / TNEUT
         ELSEIF (TATIZ.GT.0.0) THEN
            FACTA(0) = 1.0 / TATIZ
         ENDIF
         FACTB(0) = FACTA(0) * FSRATE
C
         IF (NIZS.GT.0) THEN
           DO IZ = 1, NIZS
             IF (TATIZ.GT.0.0) THEN
               FACTA(IZ) = 1.0 / TATIZ
             ELSE
               FACTA(IZ) = 0.0
             ENDIF
             FACTB(IZ) = FACTA(IZ) * QTIM
           end do
         ENDIF
c
      endif
      return
      end


      subroutine print_sputtering_yields(matt,matp,mb)
      implicit none
c
c     jdemod - this routine prints out physical and chemical sputtering 
c              yields used in the simulation as a function of impact energy
c              and temperature
c
      real :: mb
      integer :: matt,matp
c
      real e0,temp
      integer, parameter :: itmax = 8
      integer, parameter :: inmax = 40
      integer in,it
      real e0min,e0max
      real tmin,tmax
      real tmpflux 
      real,external :: yield, yldchem

      real temps(itmax),e0vals(inmax)
      real yieldc(inmax,itmax)
      real yieldp(inmax)
      real yields(inmax)

      ! print out physical chemical and self sputtering yields 

      e0min = 5.0
      e0max = 200.0
      tmin = 300.0
      tmax = 1000.0

      ! use constant flux for now ... higher fluxes reduce effective temperature
      ! for chemical sputtering

      tmpflux = 1.0e19

      
      do it = 1,itmax
         temps(it) = tmin + (it-1) * (tmax-tmin)/(itmax-1)
      end do

      do in = 1,inmax
         e0vals(in) =  e0min + (in-1) * (e0max-e0min)/(inmax-1)
      end do
         
      do in = 1,inmax
         ! last two parameters should be Te,Ti - used only for W yields
         ! physical sputtering
         yieldp(in) = yield(matp,matt,e0vals(in),10.0,10.0)
         ! self-sputtering
         yields(in) = yield(6,matt,e0vals(in),10.0,10.0)

         do it = 1,itmax

            ! chemical sputtering
            yieldc(in,it)= yldchem(e0vals(in),tmpflux,
     >                             matp,matt,temps(it))
           
         end do 


      end do

      ! print results


      write(6,*) 'SUMMARY OF SPUTTERING YIELDS:'
      write(6,'(a5,3(5x,a3,5x),20(1x,g12.5))') 'Temps:','E0','YP','YS',
     >               (temps(it),it=1,itmax)
      do in = 1,inmax
         write(6,'(a5,23(1x,g12.5))') 'Yield:',e0vals(in),yieldp(in),
     >          yields(in),(yieldc(in,it),it=1,itmax)

      end do


      end
