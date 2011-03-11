c     -*-Fortran-*-
c
      SUBROUTINE TAUIN1 (title,equil,NIZS,VFLUID)
c      SUBROUTINE TAUIN1 (NIZS,VFLUID)
c slmod begin
      USE mod_grid_divimp
      USE mod_solps
      use mtc
      use bfield
c slmod end
      implicit none
      character*(*) title,equil
      INTEGER NIZS
c slmod begin
      INTEGER GetModel
c slmod end
      REAL    VFLUID , XJI , XJ , XJF
C
C  *********************************************************************
C  *                                                                   *
C  *  TAUIN1:  INITIALISING GEOMETRY ETC                               *
C  *                                                                   *
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *
C  *                                                                   *
C  *            JAMES SPENCE   (TESSELLA)    JULY     1991             *
C  *                                                                   *
c  *            LORNE HORTON   (JET)         JANUARY  1993             *
C  *                                                                   *
C  *            GRAHAM RADFORD (JET)         JULY     1993             *
C  *                                                                   *
C  * UPDATES     : ??/07/90 --- ALLOW RINGS TO BE MISSING FROM WALL    *
C  *               08/07/91 --- ALLOW VIRTUAL RINGS (IRING=1,NXW,NXW+1)*
C  *                        --- ALLOW EDGE2D INPUT                     *
C  *               15/08/91 --- CHANGES FOR TEMP. GRAD OPT 2           *
C  *               28/01/93 --- Modify to read standard geometry file  *
C  *                            only and to only read EDGE file if     *
C  *                            the data had been requested            *
C  *                        --- BRS and BZS removed as the necessary   *
C  *                            information (KBFS) can be read         *
C  *                            directly from the equil. file (SH)     *
C  *               28/03/93 --- Read necessary information from        *
C  *                            standard geometry file so that true    *
C  *                            plasma polygon areas can be calculated *
C  *                            (in TAUVOL, stored in KVOLS for use    *
C  *                            with PIN results)                      *
C  *               15/04/93 --- Remove 1st and last points on open     *
C  *                            rings for target options which place   *
C  *                            these points outside the plasma        *
C  *               02/04/93 --- For nonorthogonal geometries read      *
C  *                            in ITAGDV from equil. file and store   *
C  *                            in TAGDV.                              *
C  *               24/05/93 --- For nonorthogonal cells calculate      *
C  *                            angle of deviation from orthogonality  *
C  *                            and store cosine of this in COSAL1/2.  *
C  *                                                                   *
C  *********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'cedge2d'
      include 'comtor'
      include 'cioniz'
      include 'reader'
      include 'dynam5'
      include 'pindata'
      include 'baffles'
      include 'printopt'
c
      include 'fperiph_com'
c
      include 'reiser_com'
c
c     Include SLCOM for optional input values
c
      include 'slcom'
c
      CHARACTER MESAGE*72,C(10)*9,FACTOR*9
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      integer in
      INTEGER IDUMMY(2),NJ(MAXNRS),NI(MAXNKS)
      INTEGER ITAG(MAXNKS*MAXNRS,5)
      INTEGER KORX(MAXNKS,MAXNRS),KORYE(MAXNRS,MAXNKS)
      INTEGER ISEP,NINOMP,NINOWA,IREF,NZ,NZMAX,K0,IK0,IKMID
      INTEGER IND,IW,IK1,IK2,IK3,IK4,IR1,IR2,IR3,IR4,IKK,IKKP1,IKMAX
      integer ikr
      real    distref,refstep,distcur
      real    vhfact,efact
      REAL    DUMMY(MAXNKS*MAXNRS,9),HMASS,ZMASS,S,DELTAL,P,FACT
      REAL    sPERP,BEST,DSQ,R,Z,THIN,THOUT,DELTAR,DELTAZ
      real    sperpt
      REAL    MU,SMAX,SMID,SDIFF
      real    qovrt
      REAL    ALPHA1,ALPHA2,ALPHA3,ALPHA4
      REAL    D0TO1,D0TO2,D0TO3,D0TO4,D1TO2,D1TO4,D3TO2,D3TO4,XCOS
      real    acos_test
      external acos_test
      REAL    DTHETA,DTHETA1
cfm
      REAL    TAUCH, TAURECTOT
cfm
      real    rf,rb,zf,zb,deltal1,deltal2
c
c     Counter for state change errors
c
      integer error_count , iwarn_count,nwarn_count,pwarn_count
      real miniontau,minneutau
      real*8 minprob
c
      character*80 comment
c
c     Ring indices 
c
      integer startir, endir
c
c     Temporary
c
      real tcar(maxizs,3),trec(3),trec_r(3),tcx(3),tcx_r(3)
c
      integer ringno
      external ringno
c
c      real    e2dvirt(maxnrs,4,2)
c
      integer  rc
      real     cenlen,outlen,inlen,polysidelen
      real     calcwav,getfracs,dist1,dist2
      external calcwav,getfracs,polysidelen
C
      integer    jplft,jprgt,tmpcpin
      integer    tmpcneur,tmpctrap
C
c      REAL    RHO(MAXNKS,MAXNRS), THETA(MAXNKS,MAXNRS)
c
c     Moved to common block
c
c      REAL    HRO(MAXNKS,MAXNRS), HTETA(MAXNKS,MAXNRS)
C
c      PINITER has been moved to the common block COMTOR so that
c      it can be accessed from within the SOL option code so that
c      they can decide whether to use default analytic ionization
c      sources or the PIN data if available.
c
      LOGICAL PINCHECK,litersol
c
      PINCHECK = .TRUE.
      PINITER  = .FALSE.
c
c     Initialize DIV BG READ flag - to off
c
      crdivbg = 0

c
c     Initialize the DperpZ Delta S transport option
c     - if this option is active, additional deltaS 
c       steps for the ions will be calculated - these 
c       deltaS steps are actually the result of dperp steps
c       in the orthogonal direction not present in DIVIMP
c       directly. 
c     This needs to be initialized prior to the neutral launch 
c     routines since it will affect HC ion transport
c
      call init_dperpz 

c
c     Initialize the background plasma arrays
c
c     Centralize the initialization in order to allow
c     for different SOL options on different rings and
c     make the code re-entrant.
c
c
c     - these two moved from SOL module
c
      CALL RZERO (KVHS, MAXNKS*MAXNRS)
      CALL RZERO (KES,  MAXNKS*MAXNRS)
      call rzero (Ktebs,maxnks*maxnrs)
      call rzero (Ktibs,maxnks*maxnrs)
      call rzero (Knbs,maxnks*maxnrs)
c
C-----------------------------------------------------------------------
c     LOAD GRID FILES
C-----------------------------------------------------------------------
c
c     Load the grid files ... depending on the option.
c
c     The routines that load the grid files are expected
c     to place the grid information into the DIVIMP
c     arrays in a standard fashion. However, the nature
c     of the ITER grid with a double null configuration
c     will require changes in most routines in order to
c     properly deal with the extra targets.
c
      if (cgridopt.eq.0) then
c
c        JET GRID FILES
c
         call rjet
      elseif (cgridopt.eq.1) then
c
c        ASDEX GRID FILES
c
         call rasdex
      elseif (cgridopt.eq.2) then
c
c        ITER GRID FILES
c
         call riter
      elseif (cgridopt.eq.3.or.cgridopt.eq.4.or.cgridopt.eq.5) then
c
c        ASDEX/TEXTOR/CMOD/DIIID - SONNET MESH FILES
c
         call raug
c slmod begin
      elseif (cgridopt.eq.LINEAR_GRID) then
         call BuildLinearGrid
      elseif (cgridopt.eq.OSM_GRID) then
         call ImportOSMGrid
         cgridopt = 3  ! switch to SONNET grid specifier
c slmod end
c
c     jdemod - support for new grid option to be added
c
c
c      elseif (cgridopt.eq.GEN_GRID) then 
c
c        jdemod - cgridopt=7 (GEN_GRID)
c
c         call readgeneralisedgrid
c
      endif
c
c     Check for a valid grid - a grid with no rings is not valid
c
      if (nrs.eq.0) then
         call prc(' INVALID GRID: NRS = 0')
         call prc(' Check that a valid grid file name was specifed')
         write (6,*) ' INVALID GRID: NRS = 0'
         write (6,*) ' Check that a valid grid file name was specifed'
         write (0,*) ' INVALID GRID: NRS = 0'
         write (0,*) ' Check that a valid grid file name was specifed'
         stop
      endif
c
C-----------------------------------------------------------------------
c
c     SET INNER/OUTER Variables for print purposes
c
C-----------------------------------------------------------------------
c
c     Set up the meaning of the words INNER and OUTER in the context
c     of this specific case. An X-point up JET case has
c     INNER = 'INNER' and OUTER = 'OUTER' while and X-point down
c     case has this orientation reversed - this will be easier than
c     adding "IF" statements throughout repeating the test to define
c     INNER and OUTER.
c
c     X-point down - rings start numbering from INNER target clockwise
c
      if (sonnet_grid_sub_type.eq.1) then  
c
         INNER='LOWER'
         OUTER='UPPER'
c
      else 
         if (zxp.le.z0) then
            INNER = 'OUTER'
            OUTER = 'INNER'
            xpoint_up = .false.
c
c        X-point UP - clockwise from OUTER target - original default
c  
         else
            INNER = 'INNER'
            OUTER = 'OUTER'
            xpoint_up = .true.
         endif
      endif
c
c
C-----------------------------------------------------------------------
c     GRID ADJUSTMENT CODE - DEPENDING ON TARGET OPTIONS
C-----------------------------------------------------------------------
c
c     Note: must be included here because it may apply to different
c           grids. (AUG vs. JET)
c
C
C     Some of the target options move the target to between the
C     last grid point and the next one in. Thus putting the first
c     and last grid points on each SOL ring outside the plasma.
c     Unfortunately, many quantities are calculated over the
c     ENTIRE grid from 1 to nks(ir) for every ring. This can cause
c     problems for the virtual points which are actually outside the
c     plasma. There are at least two solutions to the problem. First,
c     one could go through the code with a fine-toothed comb and
c     assure that the first and last grid locations and all
c     variables calculated, read-in or assigned do not reference these
c     locations. Second, the grid can be adjusted to eliminate these
c     points. Both accomplish the same objective in assuring that
c     only the relevant points are referenced ... however the
c     second does it in one centralized location that will be
c     easier to adjust at a later time if the implementation is
c     inadequate at the cost of eliminating some of the available
c     information. At this time (3/22/93) I have decided to implement
c     option 2 because it will be faster and easier to adjust at a
c     later point in time.
c
c     David Elder
c
c
c     The following is a logical variable indicating
c     the state of the grid.
c
c     virtgrid = .true. (virtual points at ends of rings
c                        are included.)
c     virtgrid = .false.(virtual points at ends of rings
c                        are excluded.)
c
      virtgrid = .true.
c
c     The following code moves all the grid positions and background
c     positions read in from the grid file for the SOL and private
c     plasma regions up one space in the array and reduces by two the
c     number of elements on the ring - thus effectively elimintaing the
c     virtual points at each end of the SOL/PP rings.
c
      if (ctargopt.eq.0.or.ctargopt.eq.1.or.ctargopt.eq.2
     >    .or.ctargopt.eq.3.or.ctargopt.eq.6) then

c
c        First calculate the target points for later use in the
c        target and wall routines for target option 1.
c
         if (ctargopt.eq.1) then
            nplat = 0
            do 290 ir = irsep , nrs
               nplat = nplat + 1
               platco(nplat,1) = ir
               platco(nplat,2) = 0.5 *(rs(nks(ir),ir)
     >                                +rs(nks(ir)-1,ir))
               platco(nplat,3) = 0.5 *(zs(nks(ir),ir)
     >                                +zs(nks(ir)-1,ir))
               platco(nplat,4) = 0.5 *(rs(1,ir)+rs(2,ir))
               platco(nplat,5) = 0.5 *(zs(1,ir)+zs(2,ir))
 290        continue

c nonorth
            if (northopt.eq.1) then
c
c             Save shifted THETA values of the target points
c
              id = 0
              do ir = irwall, irsep, -1
                 id = id + 1
                 thetat(id) = 0.5 *(thetag(nks(ir),ir)
     >                             + thetag(nks(ir)-1,ir))
              enddo
              do ir = nrs, irtrap, -1
                 id = id + 1
                 thetat(id) = 0.5 *(thetag(nks(ir),ir)
     >                             + thetag(nks(ir)-1,ir))
              enddo
              do ir = irtrap, nrs
                 id = id + 1
                 thetat(id) = 0.5 *(thetag(1,ir) + thetag(2,ir))
              enddo
              do ir = irsep, irwall
                 id = id + 1
                 thetat(id) = 0.5 *(thetag(1,ir) + thetag(2,ir))
              enddo
            endif
c nonorth


c
c        Target option 6 - using polygon edges - calculate target segment
c        centre points and evaluate THETAT values as the linear interpolation
c        of the THETAG values on the first and second grid points - keep in
c        mind that the virtual points at the end of the ring have not
c        yet been removed - and so the cells whose boundaries will form the
c        targets are for ik = 2 and ik = nks(ir)-1 .
c
         elseif (ctargopt.eq.6) then
c
c           Calculate platco values
c
c           To handle target corner for FRCv1 grids - calculate platco values for 
c           one ring less than the separatrix (would be equivalent to NRS if PFZ existed
c           for FRC grids. 
c
            if (sonnet_grid_sub_type.eq.1) then 
               startir = irsep -1 
               endir   = nrs
            else
               startir = irsep
               endir   = nrs
            endif
c
            nplat = 0
            do 291 ir = startir , endir
               nplat = nplat + 1
               platco(nplat,1) = ir
c
c              Platco 2 and 3 correspond to the first target segment
c              from ID=1 to NDSIN (Inner target in old terminology)
c              ID=1 to NDSIN corresponds to IR, NKS(IR)-1 in this
c              terminology. Platco 4 and 5 correspond to points at
c              ID = NSIN+1 to NDS for IR, IK=2.
c
c              For the first case IR,NKS(IR) or ID=1,NDSIN the elements
c              rvertp(3,in) and rvertp(4,in) are the R-coordinates of
c              the corners of the side of the polygon defining the
c              target. (same for Z). Elements 1 and 2 of the RVERTP and
c              Zvertp arrays define the target at the IR,IK=1 or
c              ID=NDSIN+1,NDS end of the ring.
c
c              INNER target (ik=nks(ir)-1, ID =1,NDSIN)
c
               in = korpg(nks(ir)-1,ir)
c              IPP/08 Krieger - ensure index of nvertp is not zero
               if (in.ne.0.and.nvertp(max(1,in)).ne.0) then
                  platco(nplat,2) = 0.5 *(rvertp(4,in)+rvertp(3,in))
                  platco(nplat,3) = 0.5 *(zvertp(4,in)+zvertp(3,in))
c
                  if (cprint.eq.3.or.cprint.eq.9)
     >              write(6,'(a,4i5,2g13.6)')
     >                  'platco: '//inner//':',in,nks(ir)-1,ir,nplat,
     >                    platco(nplat,2),platco(nplat,3)
c
               else
                  platco(nplat,2) = 0.5 *(rs(nks(ir),ir)
     >                                +rs(nks(ir)-1,ir))
                  platco(nplat,3) = 0.5 *(zs(nks(ir),ir)
     >                                +zs(nks(ir)-1,ir))
               endif
c
c              OUTER target (ik=2, ID = NDSIN+1,NDS)
c
               in = korpg(2,ir)
c              IPP/08 Krieger - ensure index of nvertp is not zero
               if (in.ne.0.and.nvertp(max(1,in)).ne.0) then
                  platco(nplat,4) = 0.5 *(rvertp(2,in)+rvertp(1,in))
                  platco(nplat,5) = 0.5 *(zvertp(2,in)+zvertp(1,in))
c
                  if (cprint.eq.3.or.cprint.eq.9)
     >               write(6,'(a,4i5,2g13.6)') 
     >                 'platco:'//outer//':',in,2,ir,nplat,
     >                     platco(nplat,4),platco(nplat,5)
c
               else
                  platco(nplat,4) = 0.5 *(rs(1,ir)+rs(2,ir))
                  platco(nplat,5) = 0.5 *(zs(1,ir)+zs(2,ir))
               endif
 291        continue


c nonorth
            if (northopt.eq.1) then
c
c           Calculate estimates of the THETA values for the targets
c           by lineraly interpolating the cell contents. This is based
c           on the assumption that the virtual points actually contain
c           valid THETA information.
c
c           Save shifted THETA values of the target points
c
              write (6,*) 'THATAT values for ctargopt = 6'
c
              id = 0
              do ir = irwall, irsep, -1
                 id = id + 1
                 thetat(id) = calcwav(nks(ir)-1,ir,
     >                        thetag(nks(ir)-1,ir),thetag(nks(ir),ir))
c
c                 write (6,*) id,ir,nks(ir)-1,thetat(id),
c     >                       thetag(nks(ir),ir),
c     >                       thetag(nks(ir)-1,ir)
              enddo
              do ir = nrs, irtrap, -1
                 id = id + 1
                 thetat(id) = calcwav(nks(ir)-1,ir,
     >                        thetag(nks(ir)-1,ir),thetag(nks(ir),ir))
c
c                 write (6,*) id,ir,nks(ir)-1,thetat(id),
c     >                       thetag(nks(ir),ir),
c     >                       thetag(nks(ir)-1,ir)
              enddo
              do ir = irtrap, nrs
                 id = id + 1
                 thetat(id) = calcwav(2,ir,thetag(2,ir),thetag(1,ir))
c
c                 write (6,*) id,ir,2,thetat(id),thetag(1,ir),
c     >                       thetag(2,ir)
              enddo
              do ir = irsep, irwall
                 id = id + 1
                 thetat(id) = calcwav(2,ir,thetag(2,ir),thetag(1,ir))
c
c                 write (6,*) id,ir,2,thetat(id),thetag(1,ir),
c     >                       thetag(2,ir)
              enddo
            endif
c nonorth

         endif

c
c        Now move the array contents of everything read in from the
c        grid file. Effectively deleting the virtual points at the
c        ends of the SOL and TRAP rings.
c

         do 292 ir = irsep,nrs
c
c           Retain an averaged KBFS value which is applicable at the
c           plates. Before any of the values are shifted.
c           Inner plate is 1 and outer plate is 2.
c
            if (ctargopt.eq.6) then
c
c              Outer plate
c
               kbfst(ir,2) = calcwav(2,ir,kbfs(2,ir),kbfs(1,ir))
c
c              Inner plate
c
               kbfst(ir,1) = calcwav(nks(ir)-1,ir,
     >                       kbfs(nks(ir)-1,ir),kbfs(nks(ir),ir))
c
            else
c
c              Outer plate
c
               kbfst(ir,2) = 0.5 * (kbfs(1,ir) + kbfs(2,ir))
c
c              Inner plate
c
               kbfst(ir,1) =0.5*(kbfs(nks(ir),ir)+kbfs(nks(ir)-1,ir))
c
            endif
c
c
            do 294 ik = 1,nks(ir)-1
               RS(IK,IR)    = rs(ik+1,ir)
               ZS(IK,IR)    = zs(ik+1,ir)
c               RHO(IK,IR)   = rho(ik+1,ir)
               RHOG(IK,IR)   = rhog(ik+1,ir)
c
               if (northopt.eq.1) then
c nonorth
                 THETAG(IK,IR)= thetag(ik+1,ir)
c nonorth
               endif
c
c               THETA(IK,IR) = theta(ik+1,ir)
               HRO(IK,IR)   = hro(ik+1,ir)
               HTETA(IK,IR) = hteta(ik+1,ir)
               BTS(IK,IR)   = bts(ik+1,ir)
               KBFS(IK,IR)  = kbfs(ik+1,ir)
c slmod begin - new
               BRATIO(IK,IR) = BRATIO(IK+1,IR)
               IF (ALLOCATED(divimp_ik)) THEN
                 divimp_ik(ik,ir) = divimp_ik(ik+1,ir)
                 divimp_ir(ik,ir) = divimp_ir(ik+1,ir)
               ENDIF
c slmod end
               PSIFL(IK,IR) = psifl(ik+1,ir)
               KORY(IR,IK)  = kory(ir,ik+1)
               KORPG(IK,IR) = korpg(ik+1,ir)
c
               if (northopt.eq.1) then
c nonorth
                 TAGDV(IK,IR) = tagdv(ik+1,ir)
c nonorth
               endif
c
 294        continue
            nks(ir) = nks(ir) - 2
 292     continue
c
c        Set virtgrid to false ... indicating no virtual points
c
         virtgrid = .false.
c
c        Adjust the IKREF and IKT value to reflect the revised grid
c
c        Need to revise ITER ikt values as well.
c
         IKREF  = NKS(IRSEP)/2 + 1
         if (ikto.gt.0.or.ikti.lt.nks(irsep)+1) then
            IKT = IKT -1
            IKTO = IKTO -1
            IKTI = IKTI -1
         else
            ikto = 0
            ikti = nks(irsep) +1
         endif
c
c        Recalculate MKS based on revised grid
c
         MKS = 0
         DO 295 IR = 1, NRS
           MKS = MAX ( MKS, NKS(IR) )
  295    CONTINUE
c
c
c        The following code adjusts any data read in from an EDGE case
c        if that is necessary.
c
         if (cre2d.eq.1.or.cre2d.eq.2.or.cioptf.eq.99.or.
     >       cioptg.eq.99.or.cioptg.eq.90) then

            do 296 ir = irsep,nrs
c
               do 298 ik = 1,nks(ir)
c
c                 Main background quantities
c
                  KNBS (IK,IR) = knbs(ik+1,ir)
                  KVHS (IK,IR) = kvhs(ik+1,ir)
                  KTIBS(IK,IR) = ktibs(ik+1,ir)
                  KTEBS(IK,IR) = ktebs(ik+1,ir)
                  KES  (IK,IR) = kes(ik+1,ir)
c
c                 Fluid code values ...
c
                  e2dNBS (IK,IR) = e2dnbs(ik+1,ir)
                  e2dnes (IK,IR) = e2dnes(ik+1,ir)
c
                  do iz = 0,cre2dizs

c
                     if (ik.eq.1) then
c
                        e2dvts(ir,2,iz) = e2dvzs(1,ir,iz)

                     elseif (ik.eq.nks(ir)) then

                        e2dvts(ir,1,iz) = e2dvzs(nks(ir)+2,ir,iz)

                     endif

                     e2dnzs(ik,ir,iz) = e2dnzs(ik+1,ir,iz)
                     e2dvzs(ik,ir,iz) = e2dvzs(ik+1,ir,iz)


                  end do
c
                  e2dz0(IK,IR)    = e2dz0(ik+1,ir)
                  e2diz0(IK,IR)   = e2diz0(ik+1,ir)
c
                  e2dVHS(IK,IR)   = e2dvhs(ik+1,ir)
                  e2dTIBS(IK,IR)  = e2dtibs(ik+1,ir)
                  e2dTEBS(IK,IR)  = e2dtebs(ik+1,ir)
                  e2dES(IK,IR)    = e2des(ik+1,ir)
c
                  e2dion(IK,IR)   = e2dion(ik+1,ir)
                  e2datom(IK,IR)  = e2datom(ik+1,ir)
                  e2dhrec(IK,IR)  = e2dhrec(ik+1,ir)
                  e2drec(IK,IR)   = e2drec(ik+1,ir)
                  e2dcxrec(ik,ir) = e2dcxrec(ik+1,ir)
                  e2dvro(ik,ir)   = e2dvro(ik+1,ir)
                  drho(ik,ir)     = drho(ik+1,ir)
                  dthetag(ik,ir)  = dthetag(ik+1,ir)
                  e2dareas(ik,ir) = e2dareas(ik+1,ir)
                  e2dgpara(ik,ir) = e2dgpara(ik+1,ir)
                  e2dgdown(ik,ir) = e2dgdown(ik+1,ir)
                  e2dpepara(ik,ir)= e2dpepara(ik+1,ir)
                  e2dpedown(ik,ir)= e2dpedown(ik+1,ir)
                  e2dpipara(ik,ir)= e2dpipara(ik+1,ir)
                  e2dpidown(ik,ir)= e2dpidown(ik+1,ir)
c
 298           continue
 296        continue
c
         endif
c
      endif
c
C
C-----------------------------------------------------------------------
C     CALCULATE ELEMENTAL VOLUMES AND AREAS
C-----------------------------------------------------------------------
C
c      CALL TAUVOL
C
C-----------------------------------------------------------------------
c     CALCULATE OTHER GEOMETRY INFORMATION
C-----------------------------------------------------------------------
c
c slmod begin - new
      CALL SetupGrid
c slmod end
C
C-----------------------------------------------------------------------
C     CALCULATE "NEAREST NEIGHBOURS" FOR EACH POINT
C-----------------------------------------------------------------------
C
C
c     JET and SONNET GRIDS
c
c slmod begin
      if (nbr.gt.0.or.eirgrid.eq.1.or.
     .    cgridopt.EQ.LINEAR_GRID) then
c...    Generalized grid:
        CALL BuildMap

      elseif (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3.or.
     >        cgridopt.eq.4.or.cgridopt.eq.5) then

c
c      if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3.or.
c     >    cgridopt.eq.4.or.cgridopt.eq.5) then
c slmod end
c
c     Orthogonal grid treatment
c
c      if (northopt.eq.0.or.northopt.eq.2) then
c
c      Newer version from JET - requires correct values for IKTI, IKTO
c
c
      DO 340 IR = 1, NRS
       DO 330 IK = 1, NKS(IR)
        IF     (IR.EQ.1) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.IRSEP-1) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IKOUTS(IK,IR) = IK + IKTO
C  1ST AND LAST POINTS ARE IDENTICAL
          IF (IK.EQ.NKS(IR)) IKOUTS(IK,IR) = IKOUTS(1,IR)
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.IRSEP) THEN
          IF (IK.LE.IKTO) THEN
            IKINS(IK,IR)= IK
            IRINS(IK,IR)= NRS
          ELSEIF (IK.GE.IKTI) THEN
            IKINS(IK,IR)= IK - (NKS(IR) - NKS(NRS))
            IRINS(IK,IR)= NRS
          ELSE
            IKINS(IK,IR)= IK - IKTO
            IRINS(IK,IR)= IR - 1
          ENDIF
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.IRWALL) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR
        ELSEIF (IR.EQ.IRTRAP) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.NRS) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IF (IK.LE.IKTO) THEN
           IKOUTS(IK,IR)= IK
           IROUTS(IK,IR)= IRSEP
          ELSE
           IKOUTS(IK,IR)= IK + (NKS(IRSEP) - NKS(IR))
           IROUTS(IK,IR)= IRSEP
          ENDIF
        ELSE
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ENDIF
C
        IKIN = IKINS(IK,IR)
        IRIN = IRINS(IK,IR)
        KINDS(IK,IR) = SQRT ((RS(IK,IR)-RS(IKIN,IRIN))**2 +
     >                       (ZS(IK,IR)-ZS(IKIN,IRIN))**2)
        IKOUT = IKOUTS(IK,IR)
        IROUT = IROUTS(IK,IR)
        KOUTDS(IK,IR) = SQRT ((RS(IK,IR)-RS(IKOUT,IROUT))**2 +
     >                        (ZS(IK,IR)-ZS(IKOUT,IROUT))**2)

       if (cprint.eq.3.or.cprint.eq.9)   
     >     write (6,'(''GRID IK/IR:'',10i5)') ik,ir,ikto,ikti,
     >           ikouts(ik,ir),irouts(ik,ir),ikins(ik,ir),irins(ik,ir)

  330  CONTINUE
  340 CONTINUE

c
c slmod begin
c...     Moved below because KPS needs to be defined for the updated metric
c        calculation code to work.
c
c        Calculate THETAG for non-orthogonal  options where it is required.
c
c slmod end
c
c     ITER GRIDS
c
      elseif (cgridopt.eq.2) then
c
      DO 1340 IR = 1, NRS
       DO 1330 IK = 1, NKS(IR)
        IF     (IR.EQ.1) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.IRSEP-1) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IF (IK .LT. IKT1) THEN
            IKOUTS(IK,IR) = IK + IKT - 1
            IROUTS(IK,IR) = IR + 1
          ELSE
            IKOUTS(IK,IR) = IK + IKT2 - IKT1 + 1
            IROUTS(IK,IR) = IRSEP2
          ENDIF
          IKOUTS(NKS(IR),IR) = IKOUTS(1,IR)
          IROUTS(NKS(IR),IR) = IROUTS(1,IR)
        ELSEIF (IR.EQ.IRSEP) THEN
          IF (IK.LT.IKT) THEN
            IKINS(IK,IR)= IK
            IRINS(IK,IR)= IRSEP2-1
          ELSEIF (IK.GE.NKS(IR)-IKT+2) THEN
            IKINS(IK,IR)= NKS(IR)-IK+1
            IRINS(IK,IR)= NRS
          ELSE
            IKINS(IK,IR)= IK - IKT + 1
            IRINS(IK,IR)= IR - 1
          ENDIF
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.IRSEP2) THEN
          IF (IK.LE.IKT2) THEN
            IKINS(IK,IR)= NKS(NRS) - IK + 1
            IRINS(IK,IR)= NRS
          ELSEIF (IK.GE.NKS(IR)-IKT2+1) THEN
            IKINS(IK,IR)= IK - (NKS(IR) - NKS(NRS))
            IRINS(IK,IR)= IRSEP2 - 1
          ELSE
            IKINS(IK,IR)= IKT1+IK-IKT2-1
            IRINS(IK,IR)= IRSEP-1
          ENDIF
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF ((IR.EQ.IRWALL) .OR. (IR .EQ. IRWALL2)) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR
        ELSEIF ((IR.EQ.IRTRAP) .OR. (IR .EQ. IRTRAP2)) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ELSEIF (IR.EQ.NRS) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IF (IK.LT.IKT) THEN
           IKOUTS(IK,IR)= NKS(IRSEP)-IK+1
           IROUTS(IK,IR)= IRSEP
          ELSE
           IKOUTS(IK,IR)= IKT2-(IK-IKT)
           IROUTS(IK,IR)= IRSEP2
          ENDIF
        ELSEIF (IR.EQ.IRSEP2-1) THEN
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IF (IK.LT.IKT) THEN
           IKOUTS(IK,IR)= IK
           IROUTS(IK,IR)= IRSEP
          ELSE
           IKOUTS(IK,IR)= IK + (NKS(IRSEP2) - NKS(IR))
           IROUTS(IK,IR)= IRSEP2
          ENDIF
        ELSE
          IKINS(IK,IR)  = IK
          IRINS(IK,IR)  = IR - 1
          IKOUTS(IK,IR) = IK
          IROUTS(IK,IR) = IR + 1
        ENDIF
C
        IKIN = IKINS(IK,IR)
        IRIN = IRINS(IK,IR)
        KINDS(IK,IR) = SQRT ((RS(IK,IR)-RS(IKIN,IRIN))**2 +
     >                       (ZS(IK,IR)-ZS(IKIN,IRIN))**2)
        IKOUT = IKOUTS(IK,IR)
        IROUT = IROUTS(IK,IR)
        KOUTDS(IK,IR) = SQRT ((RS(IK,IR)-RS(IKOUT,IROUT))**2 +
     >                        (ZS(IK,IR)-ZS(IKOUT,IROUT))**2)
 1330  CONTINUE
 1340 CONTINUE
c
c     End for different grids
c
      endif
c
c---------------------------------------------------------------
c
c     Calculate distances along the reference line. (IKREF)
c
c
      call rzero (refdist,maxnrs)
c
      distcur = 0.0
c
      do ir = irsep,irwall
         refdist(ir) = distcur
         distcur = distcur + koutds(ikref,ir)
      end do
c
c
C-----------------------------------------------------------------------
c
c     If - after the virtual points have been removed - ikto = 0
c     and ikti = nks(irsep) +1 then this indicates a situation
c     like a Textor grid with NO private plasma region. For these
c     cases there are no rings in the private plasma even though
c     they were used to remove the virtual points from the core
c     plasma rings. So set NRS = IRWALL. (At least as a temporary
c     fix).
c
      if (cgridopt.ne.2.and.(ikto.eq.0.and.ikti.eq.nks(irsep)+1)) then
         nrs = irwall
         nrs2 = nrs
         irtrap = nrs + 1
         irtrap2 = irtrap
      endif
c
C-----------------------------------------------------------------------
c
c      Initialize the Cosines even for grids being treated as Orthogonal
c
       DO IR = 1 , NRS
          DO IK = 1 , NKS(IR)
             COSALI(IK,IR) = 1.0
             COSALO(IK,IR) = 1.0
             iking(ik,ir) = ikins(ik,ir)
             ikoutg(ik,ir) = ikouts(ik,ir)
          ENDDO
       ENDDO
c
C
C-----------------------------------------------------------------------
C     FOR NONORTHOGONAL CELLS CALCULATE COSINE OF THE ANGLE DEVIATION
C     FROM ORTHOGONALITY
C-----------------------------------------------------------------------
C
c
       if (northopt.eq.1.or.northopt.eq.3) then

c nonorth
C.. Method : (1) Construct four triangles using
C                 IK,IR   = point of interest
C                 IK1,IR1 = IK+1,IR (except last knot, for CTARGOPT
C                           equals 1, 2, 3, or 5 use target point
C                           in PLATCO)
C                 IK2,IR2 = point on next inner ring but on same row
C                           as IK (except for first ring in core plasma
C                           IR=1 and first ring in trap region IRTRAP)
C                 IK3,IR3 = IK-1,IR (except first knot, for CTARGOPT
C                           equals 1, 2, 3, or 5 use target point
C                           in PLATCO)
C                 IK4,IR4 = point on next outer ring but on same row
C                           as IK (except for wall ring IRWALL)
C             (2) Use cosine rule
C                 If either current cell or inner cell are nonorthogonal
C                    COSALI(IK,IR) = cosine from the normal - INWARDS
C                 If either current cell or outer cell are nonorthogonal
C                    COSALO(IK,IR) = cosine from the normal - OUTWARDS
C                 Else COSALI = COSALO = 1.0
C
C..                                         Adapted from EDGE2D 24/05/93
C

      DO IR = 1 , NRS

         IK  = 1
         IK1 = IK + 1
         IR1 = IR
         IK2 = IKINS(IK,IR)
         IR2 = IRINS(IK,IR)
         IK4 = IKOUTS(IK,IR)
         IR4 = IROUTS(IK,IR)
C
         D0TO1 = SQRT ((RS(IK1,IR1)-RS(IK,IR))**2 +
     >                 (ZS(IK1,IR1)-ZS(IK,IR))**2)
         D0TO2 = KINDS(IK,IR)
         D0TO4 = KOUTDS(IK,IR)
         D1TO2 = SQRT ((RS(IK2,IR2)-RS(IK1,IR1))**2 +
     >                 (ZS(IK2,IR2)-ZS(IK1,IR1))**2)
         D1TO4 = SQRT ((RS(IK4,IR4)-RS(IK1,IR1))**2 +
     >                 (ZS(IK4,IR4)-ZS(IK1,IR1))**2)
C
c        Target options 0 and 4
c
         IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4.or.ir.lt.irsep)THEN
             IF( D0TO2.GT.0.0.and.d0to1.gt.0.0)THEN
                 XCOS   = (D0TO2**2 - D1TO2**2 + D0TO1**2)
     >                    / (2.0 * D0TO1 * D0TO2)
                 ALPHA1 = ACOS_TEST(XCOS,1)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK2,IR2).GT.0 )
c sl begin
     >               COSALI(IK,IR) = COS(PI/2.0-ALPHA1)
c     >               COSALI(IK,IR) = COS( ALPHA1 )
c sl end
             ENDIF
             IF( D0TO4.GT.0.0.and.d0to1.gt.0.0)THEN
                 XCOS   = (D0TO4**2 - D1TO4**2 + D0TO1**2)
     >                    / (2.0 * D0TO4 * D0TO1)
                 ALPHA4 = ACOS_TEST(XCOS,2)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK4,IR4).GT.0 )
c sl begin
     >               COSALO(IK,IR) = COS(PI/2.0-ALPHA4)
c     >               COSALO(IK,IR) = COS( ALPHA4 )
c sl end
             ENDIF
c
c        All other target options
c
         ELSE
             DO ICOUNT = 1,NPLAT
                IF( PLATCO(ICOUNT,1).EQ.IR )ID = ICOUNT
             ENDDO
             D0TO3 = SQRT ((PLATCO(ID,4)-RS(IK,IR))**2 +
     >                     (PLATCO(ID,5)-ZS(IK,IR))**2)
             D3TO2 = SQRT ((RS(IK2,IR2)-PLATCO(ID,4))**2 +
     >                     (ZS(IK2,IR2)-PLATCO(ID,5))**2)
             D3TO4 = SQRT ((RS(IK4,IR4)-PLATCO(ID,4))**2 +
     >                     (ZS(IK4,IR4)-PLATCO(ID,5))**2)
C
             IF( D0TO2.GT.0.0.and.d0to1.gt.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO2**2 - D1TO2**2 + D0TO1**2)
     >                    / (2.0 * D0TO1 * D0TO2)
                 ALPHA1 = ACOS_TEST(XCOS,3)
                 XCOS   = (D0TO2**2 - D3TO2**2 + D0TO3**2)
     >                    / (2.0 * D0TO3 * D0TO2)
                 ALPHA2 = ACOS_TEST(XCOS,4)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK2,IR2).GT.0 )
     >               COSALI(IK,IR) = COS( 0.5*(ALPHA1+ALPHA2)-ALPHA1 )
             ENDIF
             IF( D0TO4.GT.0.0.and.d0to1.gt.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO4**2 - D3TO4**2 + D0TO3**2)
     >                    / (2.0 * D0TO3 * D0TO4)
                 ALPHA3 = ACOS_TEST(XCOS,5)
                 XCOS   = (D0TO4**2 - D1TO4**2 + D0TO1**2)
     >                    / (2.0 * D0TO4 * D0TO1)
                 ALPHA4 = ACOS_TEST(XCOS,6)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK4,IR4).GT.0 )
     >               COSALO(IK,IR) = COS( 0.5*(ALPHA3+ALPHA4)-ALPHA4 )
             ENDIF
         ENDIF
C
         DO IK = 2     , NKS(IR)-1
C
            IK1 = IK + 1
            IR1 = IR
            IK2 = IKINS(IK,IR)
            IR2 = IRINS(IK,IR)
            IK3 = IK - 1
            IR3 = IR
            IK4 = IKOUTS(IK,IR)
            IR4 = IROUTS(IK,IR)
C
            D0TO1 = SQRT ((RS(IK1,IR1)-RS(IK ,IR ))**2 +
     >                    (ZS(IK1,IR1)-ZS(IK ,IR ))**2)
            D0TO2 = KINDS(IK,IR)
            D0TO3 = SQRT ((RS(IK3,IR3)-RS(IK ,IR ))**2 +
     >                    (ZS(IK3,IR3)-ZS(IK ,IR ))**2)
            D0TO4 = KOUTDS(IK,IR)
            D1TO2 = SQRT ((RS(IK2,IR2)-RS(IK1,IR1))**2 +
     >                    (ZS(IK2,IR2)-ZS(IK1,IR1))**2)
            D1TO4 = SQRT ((RS(IK4,IR4)-RS(IK1,IR1))**2 +
     >                    (ZS(IK4,IR4)-ZS(IK1,IR1))**2)
            D3TO2 = SQRT ((RS(IK2,IR2)-RS(IK3,IR3))**2 +
     >                    (ZS(IK2,IR2)-ZS(IK3,IR3))**2)
            D3TO4 = SQRT ((RS(IK4,IR4)-RS(IK3,IR3))**2 +
     >                    (ZS(IK4,IR4)-ZS(IK3,IR3))**2)
C
            IF( D0TO2.GT.0.0.and.d0to1.gt.0.0.and.d0to3.gt.0.0)THEN
                XCOS   = (D0TO1**2 + D0TO2**2 - D1TO2**2)
     >                   / (2.0 * D0TO1 * D0TO2)
                ALPHA1 = ACOS_TEST(XCOS,7)
                XCOS   = (D0TO3**2 + D0TO2**2 - D3TO2**2)
     >                   / (2.0 * D0TO3 * D0TO2)
                ALPHA2 = ACOS_TEST(XCOS,8)
                IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK2,IR2).GT.0)
     >              COSALI(IK,IR) = COS( 0.5*(ALPHA1+ALPHA2) - ALPHA1 )
            ENDIF

            IF( D0TO4.GT.0.0.and.d0to3.gt.0.0.and.d0to1.gt.0.0)THEN
                XCOS   = (D0TO3**2 + D0TO4**2 - D3TO4**2)
     >                   / (2.0 * D0TO3 * D0TO4)
                ALPHA3 = ACOS_TEST(XCOS,9)
                XCOS   = (D0TO4**2 + D0TO1**2 - D1TO4**2)
     >                   / (2.0 * D0TO4 * D0TO1)
                ALPHA4 = ACOS_TEST(XCOS,10)
                IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK4,IR4).GT.0 )
     >              COSALO(IK,IR) = COS( 0.5*(ALPHA3+ALPHA4) - ALPHA4 )
            ENDIF
         ENDDO

         IK  = NKS(IR)
         IK2 = IKINS(IK,IR)
         IR2 = IRINS(IK,IR)
         IK3 = IK - 1
         IR3 = IR
         IK4 = IKOUTS(IK,IR)
         IR4 = IROUTS(IK,IR)
         D0TO2 = KINDS(IK,IR)
         D0TO3 = SQRT ((RS(IK3,IR3)-RS(IK,IR))**2 +
     >                 (ZS(IK3,IR3)-ZS(IK,IR))**2)
         D0TO4 = KOUTDS(IK,IR)
         D3TO2 = SQRT ((RS(IK2,IR2)-RS(IK3,IR3))**2 +
     >                 (ZS(IK2,IR2)-ZS(IK3,IR3))**2)
         D3TO4 = SQRT ((RS(IK4,IR4)-RS(IK3,IR3))**2 +
     >                 (ZS(IK4,IR4)-ZS(IK3,IR3))**2)
C
c        Target options 0 and 4
c
         IF( CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4.or.ir.lt.irsep)THEN
             IF( D0TO2.GT.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO3**2 + D0TO2**2 - D3TO2**2)
     >                    / (2.0 * D0TO3 * D0TO2)
                 ALPHA2 = ACOS_TEST(XCOS,11)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK2,IR2).GT.0 )
c sl begin
     >               COSALI(IK,IR) = COS(PI/2.0-ALPHA2)
c     >               COSALI(IK,IR) = COS( ALPHA2 )
c sl end
             ENDIF
             IF( D0TO4.GT.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO3**2 + D0TO4**2 - D3TO4**2)
     >                    / (2.0 * D0TO3 * D0TO4)
                 ALPHA3 = ACOS_TEST(XCOS,12)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK4,IR4).GT.0 )
c sl begin
     >               COSALO(IK,IR) = COS(PI/2.0-ALPHA3)
c     >               COSALO(IK,IR) = COS( ALPHA3 )
c sl end
             ENDIF
c
c        All other target options
c
         ELSE
C
             D0TO1 = SQRT ((PLATCO(ID,2)-RS(IK,IR))**2 +
     >                     (PLATCO(ID,3)-ZS(IK,IR))**2)
             D1TO2 = SQRT ((RS(IK2,IR2)-PLATCO(ID,2))**2 +
     >                     (ZS(IK2,IR2)-PLATCO(ID,3))**2)
             D1TO4 = SQRT ((RS(IK4,IR4)-PLATCO(ID,2))**2 +
     >                     (ZS(IK4,IR4)-PLATCO(ID,3))**2)
C
             IF( D0TO2.GT.0.0.and.d0to1.gt.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO1**2 + D0TO2**2 - D1TO2**2)
     >                    / (2.0 * D0TO1 * D0TO2)
                 ALPHA1 = ACOS_TEST(XCOS,13)
                 XCOS   = (D0TO3**2 + D0TO2**2 - D3TO2**2)
     >                    / (2.0 * D0TO3 * D0TO2)
                 ALPHA2 = ACOS_TEST(XCOS,14)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK2,IR2).GT.0 )
     >               COSALI(IK,IR) = COS( 0.5*(ALPHA1+ALPHA2)-ALPHA1 )
             ENDIF
c
             IF( D0TO4.GT.0.0.and.d0to1.gt.0.0.and.d0to3.gt.0.0)THEN
                 XCOS   = (D0TO3**2 + D0TO4**2 - D3TO4**2)
     >                    / (2.0 * D0TO3 * D0TO4)
                 ALPHA3 = ACOS_TEST(XCOS,15)
                 XCOS   = (D0TO4**2 + D0TO1**2 - D1TO4**2)
     >                    / (2.0 * D0TO4 * D0TO1)
                 ALPHA4 = ACOS_TEST(XCOS,16)
                 IF( TAGDV(IK,IR).GT.0 .OR. TAGDV(IK4,IR4).GT.0 )
     >               COSALO(IK,IR) = COS( 0.5*(ALPHA3+ALPHA4)-ALPHA4 )
             ENDIF
         ENDIF


      ENDDO
c
      if (cprint.eq.3.or.cprint.eq.9) then
         WRITE(6,'(//1X,''COS ALPHA:'')')
         DO IR = 1, NRS
           WRITE (6,*)' IK  IR   COSALI COSALO  TAGDV'
           DO IK = 1, NKS(IR)
             WRITE (6,'(2I4,2X,2F7.4,I5)') IK,IR,
     >                     COSALI(IK,IR),COSALO(IK,IR),TAGDV(IK,IR)
           ENDDO
         ENDDO
      endif

C
C-----------------------------------------------------------------------
C     CALCULATE TRUE NEAREST KNOTS FOR EACH POINT FOR NON-ORTHOGONAL
C     GEOMETRIES
C-----------------------------------------------------------------------
C     NOTE: 1.Uses the THETA array read in from the geometry file.
C           2.Results stored in arrays IKING(IK,IR) and IKOUTG(IK,IR).
C           3.IKING and IKOUTG will only differ from IKINS and IKOUTS
C             for the G type geometry files generated by GRID2D.
C           4.There is a mismatch of THETA values between rings IRSEP
C             and NRS for knots greater than that corresponding to
C             JPLFT.  This is accounted for here and in subroutine DIV
C             when crossing between these two rings (DTHETG).
C             Thus, on the separatrix, it is assumed that the knots
C             corresponding to JPRGT (IKTI-1) and JPLFT (IKTO+1) lie on
C             an orthogonal grid in which THETA has integral values.
C             Any other points may be non-orthogonal.
C
      DTHETG = THETAG(IKTI-1,IRSEP) - THETAG(IKTO+1,IRSEP) + 1.0
      DO IR = 1, NRS
        IF( IR.EQ.IRSEP )THEN
          DO IK = 1, IKTI-1
            IKING(IK,IR) = 1
            IRIN = IRINS(IK,IR)
            IKMAX = NKS(IRIN)
            IF (IRIN.LT.IRSEP) IKMAX = NKS(IRIN) - 1
            DO IKK = 1, IKMAX
              IKKP1 = MIN0( IKK+1,NKS(IRIN) )
              DTHETA = THETAG(IKK,IRIN) - THETAG(IK,IR)
              DTHETA1 = THETAG(IKKP1,IRIN) - THETAG(IK,IR)
              IF( DTHETA.LE.0.0 )THEN
                IKING(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKING(IK,IR) = IKKP1
              ENDIF
              IKING(IK,IR) = MAX0( IKING(IK,IR),1 )
            ENDDO
            IKIN = IKING(IK,IR)
            IF (TAGDV(IK,IR).EQ.0 .AND. TAGDV(IKIN,IRIN).EQ.0)
     >        IKING(IK,IR) = IKINS(IK,IR)
          ENDDO
          DO IK = IKTI, NKS(IR)
            IKING(IK,IR) = 1
            IRIN = IRINS(IK,IR)
            DO IKK = 1, NKS(IRIN)
              IKKP1 = MIN0( IKK+1,NKS(IRIN) )
              DTHETA = THETAG(IKK,IRIN) - (THETAG(IK,IR) - DTHETG)
              DTHETA1 = THETAG(IKKP1,IRIN) - (THETAG(IK,IR) - DTHETG)
              IF( DTHETA.LE.0.0 )THEN
                IKING(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKING(IK,IR) = IKKP1
              ENDIF
              IKING(IK,IR) = MAX0( IKING(IK,IR),1 )
            ENDDO
          ENDDO
        ELSE
          DO IK = 1, NKS(IR)
            IKING(IK,IR) = 1
            IRIN = IRINS(IK,IR)
            IKMAX = NKS(IRIN)
            IF (IRIN.LT.IRSEP) IKMAX = NKS(IRIN) - 1
            DO IKK = 1, IKMAX
              IKKP1 = MIN0( IKK+1,NKS(IRIN) )
              DTHETA = THETAG(IKK,IRIN) - THETAG(IK,IR)
              DTHETA1 = THETAG(IKKP1,IRIN) - THETAG(IK,IR)
              IF( DTHETA.LE.0.0 )THEN
                IKING(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKING(IK,IR) = IKKP1
              ENDIF
              IKING(IK,IR) = MAX0( IKING(IK,IR),1 )
            ENDDO
            IKIN = IKING(IK,IR)
            IF (TAGDV(IK,IR).EQ.0 .AND. TAGDV(IKIN,IRIN).EQ.0)
     >        IKING(IK,IR) = IKINS(IK,IR)
          ENDDO
        ENDIF
        IF( IR.EQ.NRS )THEN
          DO IK = 1, NKS(IR)
            IROUT = IROUTS(IK,IR)
            DO IKK = 1, IKTO
              IKKP1 = MIN0( IKK+1,IKTO )
              DTHETA = THETAG(IKK,IROUT) - THETAG(IK,IR)
              DTHETA1 = THETAG(IKKP1,IROUT) - THETAG(IK,IR)
              IF( DTHETA.LE.0.0 )THEN
                IKOUTG(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKOUTG(IK,IR) = IKKP1
              ENDIF
              IKOUTG(IK,IR) = MAX0( IKOUTG(IK,IR),1 )
            ENDDO
            DO IKK = IKTI, NKS(IROUT)
              IKKP1 = MIN0( IKK+1,NKS(IROUT) )
              DTHETA = THETAG(IKK,IROUT) - DTHETG - THETAG(IK,IR)
              DTHETA1 = THETAG(IKKP1,IROUT) - DTHETG - THETAG(IK,IR)
              IF( DTHETA.LE.0.0 )THEN
                IKOUTG(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKOUTG(IK,IR) = IKKP1
              ENDIF
              IKOUTG(IK,IR) = MAX0( IKOUTG(IK,IR),1 )
            ENDDO
          ENDDO
        ELSE
          DO IK = 1, NKS(IR)
            IROUT = IROUTS(IK,IR)
            IKMAX = NKS(IROUT)
            IF (IRIN.LT.IRSEP) IKMAX = NKS(IROUT) - 1
            DO IKK = 1, IKMAX
              IKKP1 = MIN0( IKK+1,NKS(IROUT) )
              DTHETA = THETAG(IKK,IROUT) - THETAG(IK,IR)
              DTHETA1 = THETAG(IKKP1,IROUT) - THETAG(IK,IR)
              IF( DTHETA.LE.0.0 )THEN
                IKOUTG(IK,IR) = IKK
                IF( ABS(DTHETA1).LT.ABS(DTHETA) )IKOUTG(IK,IR) = IKKP1
              ENDIF
              IKOUTG(IK,IR) = MAX0( IKOUTG(IK,IR),1 )
            ENDDO
            IKOUT = IKOUTG(IK,IR)
            IF (TAGDV(IK,IR).EQ.0 .AND. TAGDV(IKOUT,IROUT).EQ.0)
     >        IKOUTG(IK,IR) = IKOUTS(IK,IR)
          ENDDO
        ENDIF
c
        if (cprint.eq.3.or.cprint.eq.9) then
           WRITE(6,*)' IR  IK   IKING   IKINS  IKOUTG  IKOUTS'
           DO IK = 1, NKS(IR)
             WRITE(6,'(2I4,4I8)')IR,IK,IKING(IK,IR),IKINS(IK,IR),
     >                           IKOUTG(IK,IR),IKOUTS(IK,IR)
           ENDDO
        endif
c
      ENDDO
c nonorth
      endif
c slmod begin 
c
c Find IKING and IKOUTG... check this out and make sure it should
c be here...
c
      IF (grdnmod.NE.0) THEN
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            iking (ik,ir) = ikins (ik,ir)
            ikoutg(ik,ir) = ikouts(ik,ir)
          ENDDO
        ENDDO
      ENDIF
c slmod end
c
c
      if (cprint.eq.3.or.cprint.eq.9) then
         write(6,*) 'Connection map:'
c
         do 3000 ir = 1,nrs
           do 3010 ik = 1,nks(ir)
             write(6,'(a,2i5,a,4i5,2(1x,g13.6))') 
     >        'ind:',ik,ir,':',ikins(ik,ir),irins(ik,ir),
     >        ikouts(ik,ir),irouts(ik,ir),kinds(ik,ir),koutds(ik,ir)
 3010      continue
 3000    continue
      endif
c
C-----------------------------------------------------------------------
c     Calculate the cross-field distance ratios for transport
c     from one cell to the next. Based on kinds and koutds.
C-----------------------------------------------------------------------
c
      if (cfdopt.eq.0) then
c
c        Original standard option - cross-field cell boundary is
c        midway between cell centres.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               finds(ik,ir) = 0.5
               foutds(ik,ir) = 0.5
            end do
         end do
c
      elseif (cfdopt.eq.1) then
c
c        Cross-field distances are based on the distance
c        from cell centre to the polygon boundary on the
c        line-segment joining the cell centres.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
c
               finds(ik,ir)=getfracs(ik,ir,ikin,irin,1)
c
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               foutds(ik,ir)=getfracs(ik,ir,ikout,irout,3)
c
            end do
         end do
c
c        Adjust finds and foutds for boundary rings -
c
c        Core boundary
c
         ir =1
         do ik = 1,nks(ir)
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
            foutds(ik,ir) = 1.0 - finds(ikout,irout)
         end do
c
c        PP boundary
c
c afmod begin
c         IF (ippchange) THEN
c        Krieger IPP/05 - next loop only if trap exists	 
c
c        jdemod - only execute the following code if the PFZ exists
c
	   IF (irtrap.LE.nrs) then
	     ir=irtrap
             do ik = 1,nks(ir)
              ikout = ikouts(ik,ir)
              irout = irouts(ik,ir)
              foutds(ik,ir) = 1.0 - finds(ikout,irout)
             end do
	   endif  
c
c
c	 else 
c           ir =irtrap
c           do ik = 1,nks(ir)
c              ikout = ikouts(ik,ir)
c              irout = irouts(ik,ir)
c              foutds(ik,ir) = 1.0 - finds(ikout,irout)
c           end do
c         end if
c
c         ir =irtrap
c         do ik = 1,nks(ir)
c            ikout = ikouts(ik,ir)
c            irout = irouts(ik,ir)
c            foutds(ik,ir) = 1.0 - finds(ikout,irout)
c         end do
c afmod end
c
c        Wall boundary
c
         ir =irwall
         do ik = 1,nks(ir)
            ikin = ikins(ik,ir)
            irin = irins(ik,ir)
            finds(ik,ir) = 1.0 - foutds(ikin,irin)
         end do
c
c        Print - out - debug information
c
         if (cprint.eq.3.or.cprint.eq.9) then
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               ikin = ikins(ik,ir)
               irin = irins(ik,ir)
               ikout = ikouts(ik,ir)
               irout = irouts(ik,ir)
c
               write (6,'(a,6i4,8g12.5)') 'FRAC:',ik,ir,
     >                    ikin,irin,
     >                      ikout,irout,finds(ik,ir),
     >                      foutds(ik,ir),
     >                   finds(ik,ir)+foutds(ikin,irin),
     >                   finds(ikout,irout)+foutds(ik,ir),
     >                   kinds(ik,ir),koutds(ik,ir),
     >                   cosali(ik,ir)*kinds(ik,ir),
     >                   cosalo(ik,ir)*koutds(ik,ir)
c
               if (abs(1.0-(finds(ik,ir)+foutds(ikin,irin)))
     >                 .gt.1.0e-3) then
c
                  in = korpg(ik,ir)
c
c                 IPP/08 Krieger - ensure index of nvertp is not zero
                  if (in.gt.0.and.nvertp(max(1,in)).gt.0) then
c
                     write (6,'(a,6i4,2g15.6)') 'IN  ERROR:',
     >                       ik,ir,ikin,irin,
     >                       in,nvertp(in),rs(ik,ir),zs(ik,ir)

                     write (6,'(a,4g15.6)') 'SIDE 1:',
     >                 rvertp(nvertp(in),in),zvertp(nvertp(in),in),
     >                 rvertp(1,in),zvertp(1,in)

                  else

                     write (6,'(a,5i4,2g15.6)') 'IN  ERROR:',
     >                       ik,ir,ikin,irin,
     >                       in,rs(ik,ir),zs(ik,ir)

                  endif

               endif

               if (abs(1.0-(foutds(ik,ir)+finds(ikout,irout)))
     >                      .gt.1.0e-3) then

                  in = korpg(ik,ir)
c

c                 IPP/08 Krieger - ensure index of nvertp is not zero
                  if (in.gt.0.and.nvertp(max(1,in)).gt.0) then

                     write (6,'(a,6i4,2g15.6)') 'OUT ERROR:',
     >                       ik,ir,ikout,irout,
     >                       in,nvertp(in),rs(ik,ir),zs(ik,ir)

                     write (6,'(a,4g15.6)') 'SIDE 3:',
     >                 rvertp(2,in),zvertp(2,in),
     >                 rvertp(3,in),zvertp(3,in)


                  else


                     write (6,'(a,5i4,2g15.6)') 'OUT ERROR:',
     >                       ik,ir,ikout,irout,
     >                       in,rs(ik,ir),zs(ik,ir)

                  endif

               endif
c
            end do
         end do
c
         endif
c
      endif
c
c
C-----------------------------------------------------------------------
c     Calculate the GRID distances and the actual perpendicular
c     to field line distances between each set of cell centres.
C-----------------------------------------------------------------------
c
c     Cosali and cosalo are set to 1.0 for orthogonal grid options.
c
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
            distin(ik,ir) = finds(ik,ir) * kinds(ik,ir) * cosali(ik,ir)
            distout(ik,ir)= foutds(ik,ir)* koutds(ik,ir)* cosalo(ik,ir)
c
            if (cprint.eq.1.or.cprint.eq.3.or.cprint.eq.9) then 
               write(6,'(a,2i6,10(1x,g12.5))') 'DISTS:', ik,ir,
     >                 distin(ik,ir),distout(ik,ir)
c
            endif

c
         end do
      end do
c
c     Calculate total distances.
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
            ikin = ikins(ik,ir)
            irin = irins(ik,ir)
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
c
            tdistin(ik,ir) = distin(ik,ir) + distout(ikin,irin)
            tdistout(ik,ir)= distout(ik,ir)+ distin(ikout,irout)
c
         end do
      end do
c
c
C-----------------------------------------------------------------------
c     Check that all grid polygons are ordered correctly.
C-----------------------------------------------------------------------
c
c slmod begin
c...  This does not work with generalized geometry:

      IF (grdnmod.EQ.0) call grid_check
c
c      call grid_check
c slmod end
c
C
C-----------------------------------------------------------------------
C     CALCULATE THE PLATES AND WALLS
C-----------------------------------------------------------------------
C
C     TARGETS
C
      DO ir = 1, nrs
        idds(ir,1) = MAX(MIN(idds(ir,1),999),-999)
        idds(ir,2) = MAX(MIN(idds(ir,2),999),-999)
      ENDDO
c
      CALL DOTARG

c      CALL BuildTargets
c      CALL OutputData(86,'After calling DOTARG')
c      CALL DumpGrid('After calling DOTARG')
     

C
C     NEUTRAL WALLS
C
c
c     If the NIMBUS wall has been specified for use - then
c     set the neutral wall option to the GRID wall
c     for the first iteration of DOWALL - the coordinates
c     of the actual wall will not be available until after PIN
c     is run.
c
c     Record the original option values for later use.
c
      tmpcneur = cneur
      tmpctrap = ctrap
c
      if (cneur.eq.5) then
         cneur = 4
      endif
c
      if (ctrap.eq.5) then
         ctrap = 4
      endif
c
c slmod begin - new

      IF (nbr.GT.0.OR.grdnmod.NE.0.OR.eirgrid.EQ.1) THEN
c...    Generalized grid:
        CALL BuildNeutralWall
      ELSE
        CALL DOWALL
      ENDIF

c...  This assigns the xVESM arrays from the WALLPT array.  It also
c     adds the wall segments that are specified in the DIVIMP
c     input file (if any, see unstructured input tag 077).  It may
c     be better to put this elsewhere, especially if the NIMBUS wall
c     is overwritten after a call to EIRENE in some circumstances (not 
c     sure on this at the moment).  This routine is called for all grids
c     except JET grids: (NOW called for JET grids if Eirene is run)

      IF (pincode.EQ.1.OR.pincode.EQ.2.OR.pincode.EQ.3.OR.
     .    pincode.EQ.4.OR.pincode.EQ.5) 
     .  CALL AssignNIMBUSWall

c...  Build vacuum grid:

      IF (stopopt .EQ.3.OR.stopopt .EQ.108.OR.stopopt3.EQ. 6.OR.
     .    stopopt3.EQ.7.OR.stopopt3.EQ. 12.OR.stopopt3.EQ.13.OR.
     .    stopopt3.EQ.9.OR.
     .    (eirbgk.GE.2.AND.eirbgk.LE.4.AND.vacnseg.GT.0.AND.
     .     pincode.NE.4.AND.pincode.NE.5))
     .  CALL BuildVacuumGrid

c...  Setup data for toroidal NBLOCK switching surfaces:
      IF (eirnasdat.GT.0) CALL SetupToroidalSurfaces

      IF (grdnmod.NE.0) THEN
c...    Assign core (R,ZCW) and grid boundary (R,ZIW) polygons for ion
c       transport routines in DIVIMP:
        CALL BuildGridPolygons
      ELSE
C
C       ION WALLS - JUST THE POLYGON POINTS FOR WALLEDGE ROUTINE.
C
        CALL IONWALL
      ENDIF
c
c      CALL DOWALL
c
c
C
C     ION WALLS - JUST THE POLYGON POINTS FOR WALLEDGE ROUTINE.
C
c      CALL IONWALL
c slmod end
C
c
c     This routine examines grid orthogonality. Calculates the orthogonal
c     angle for each cell and the difference from orthogonality for the
c     target segments.
c
      call calcorth
c
C
C-----------------------------------------------------------------------
C     CALCULATE DISTANCES ALONG CONTOURS S, MAGNETIC FIELD FACTORS ETC
C-----------------------------------------------------------------------
C
c
c     KINS is the probability of an inward step - typically it is set
c     to 0.5 - however - this can cause an unrealistic pinch effect
c     since in a curved geometry - the probability of stepping inward
c     is not the same as stepping outward. This is approximated by using
c     the length of the sides of the polygon which lie parallel to the
c     magnetic field as approximate estimators of the probability of a
c     particle making a step in that direction.
c
c
      DO 320 IR = 1, NRS
        DO 310 IK = 1, NKS(IR)
          IF (IR.EQ.1) THEN
            KINS(IK,IR) = 0.0
          ELSE
            KINS(IK,IR) = 0.5
          ENDIF
  310   CONTINUE
c
c       The following code calculates the distances along the field
c       lines - they mey be calculated two ways:
c
c       0) Calculated by joining the centre points of cells.
c          This is the older way and works for grids lacking polygon
c          information.
c
c       1) Calculated through the mid-points of the polygon sides that
c          are parallel to the field lines - these are slightly longer
c          but the added precision is needed for non-orthogonal
c          treatment option 3 and for SOL Option 22. IF PDOPT=1
c          has been selected the quantities calculated for KSS2, KPS2 ...
c          are substituted for the values used for KSS and KPS.
c
c          David Elder and Steven Lisgo, August 26, 1996
c
c          Calculate standard quantities KSS ...
c

        IF (CTARGOPT.EQ.0.or.ctargopt.eq.4.OR.IR.LT.IRSEP) THEN
           S = 0.0
           P = 0.0
           KSS(1,IR) = 0.0
           KPS(1,IR) = 0.0
        ELSEIF (CTARGOPT.EQ.1.or.ctargopt.eq.2.or.
     >          ctargopt.eq.3.or.ctargopt.eq.5.or.ctargopt.eq.6) THEN
           DELTAL = SQRT ((RP(IDDS(IR,2))-RS(1,IR)) ** 2 +
     >                   (ZP(IDDS(IR,2))-ZS(1,IR)) ** 2)
           if (ctargopt.eq.5) then
              S      = KBFS(1,IR) * DELTAL
           elseif (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >             ctargopt.eq.3.or.ctargopt.eq.6) then
              S      = KBFST(IR,2) * DELTAL
           endif
           P      = DELTAL
           KSS(1,IR) = S
           KPS(1,IR) = P
        ENDIF
C
        DO 315 IK = 2, nks(ir)
          DELTAL = SQRT ((RS(IK,IR)-RS(IK-1,IR)) ** 2 +
     >                   (ZS(IK,IR)-ZS(IK-1,IR)) ** 2)
          S      = S + 0.5*(KBFS(IK,IR)+KBFS(IK-1,IR)) * DELTAL
          P      = P + DELTAL
          KSS(IK,IR) = S
          KPS(IK,IR) = P
  315   CONTINUE
C
C
        IF (CTARGOPT.EQ.0.or.ctargopt.eq.4.or.ir.lt.irsep) then
           KSMAXS(IR) = S
           KPMAXS(IR) = P
        ELSEIF (CTARGOPT.EQ.1.or.ctargopt.eq.2.or.
     >          ctargopt.eq.3.or.ctargopt.eq.5.or.ctargopt.eq.6) THEN
           DELTAL = SQRT ((RP(IDDS(IR,1))-RS(NKS(IR),IR)) ** 2 +
     >                   (ZP(IDDS(IR,1))-ZS(NKS(IR),IR)) ** 2)
           if (ctargopt.eq.5) then
              S      = S + KBFS(NKS(IR),IR) * DELTAL
           elseif (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >             ctargopt.eq.3.or.ctargopt.eq.6) then
              S      = S + KBFST(IR,1) * DELTAL
           endif
           P      = P + DELTAL
           KSMAXS(IR) = S
           KPMAXS(IR) = P
        ENDIF
C
        IF (IR.LT.IRSEP) THEN
          KKS(IR)  = 1.0 + 0.05 * REAL(IR-IRSEP)
        ELSEIF (ir.ge.irsep.and.ir.le.irwall2) then
          KKS(IR)  = 1.0 + 0.05 * REAL(IR-IRSEP)
        ELSEIF (ir.ge.irtrap2.and.ir.le.nrs2) then
          KKS(IR)  = 1.0 + 0.05 * REAL(IR-NRS2-1)
        ELSEIF (ir.ge.irsep2.and.ir.le.irwall) then
          KKS(IR)  = 1.0 + 0.05 * REAL(IR-IRSEP2)
        ELSEIF (ir.ge.irtrap.and.ir.le.nrs) then
          KKS(IR)  = 1.0 + 0.05 * REAL(IR-NRS-1)
        ENDIF
C
  320 CONTINUE
c
c     Calculate the KSS2 and KPS2 arrays - based on the complete
c     grid - using the mid-points of the sides to define the S limits
c     of the cell. KSB contains the upper S-boundary of each cell.
c     KPB - contains the upper P boundary of each cell.
c
c     JET and SONNET grids
c
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.
     >    cgridopt.EQ.LINEAR_GRID) then
c
c        Only possible for grids with polygon information
c
      DO 322 IR = 1, NRS
c
        if (cprint.eq.3.or.cprint.eq.9)
     >     write (6,'(a,i6)') 'Values of KSS2 KPS2'//
     >                        ' KSS KSB for ring - ',ir
c
c       Calculate beginning of ring.
c
        in = korpg(1,ir)
c
c       If a polygon exists for the cell.
c
c       The code assumes that all cells will exist for a given
c       ring and that there won't be gaps in rings without defined
c       polygons. If the first cell doesn't have a polygon it just
c       uses the KSS values found previously for KSS2.
c
c
        if (in.ne.0) then
c
           rb = (rvertp(1,in) + rvertp(2,in)) /2.0
           rf = (rvertp(3,in) + rvertp(4,in)) /2.0
           zb = (zvertp(1,in) + zvertp(2,in)) /2.0
           zf = (zvertp(3,in) + zvertp(4,in)) /2.0
c
           ksb(0,ir) = 0.0
           kpb(0,ir) = 0.0
c
           krb(0,ir) = rb
           krb(1,ir) = rf
c
           kzb(0,ir) = zb
           kzb(1,ir) = zf
c
c          -------------------------------------
c 
c          IK = 1
c
c          Should be indented ...
c
           IF (CTARGOPT.EQ.0.or.ctargopt.eq.4.OR.IR.LT.IRSEP) THEN
              S = 0.0
              P = 0.0
              KSS2(1,IR) = 0.0
              KPS2(1,IR) = 0.0
              DELTAL1 = 0.0
              DELTAL2 = sqrt((rs(1,ir)-rf)**2+(zs(1,ir)-zf)**2)
              ksb(1,ir) = deltal2 * kbfs(1,ir)
              kpb(1,ir) = deltal2
           ELSEIF (CTARGOPT.EQ.1.or.ctargopt.eq.2.or.
     >          ctargopt.eq.3.or.ctargopt.eq.5.or.ctargopt.eq.6) THEN
c
              DELTAL1 = SQRT ((RP(IDDS(IR,2))-RS(1,IR)) ** 2 +
     >                   (ZP(IDDS(IR,2))-ZS(1,IR)) ** 2)
c
              DELTAL2 = sqrt((rs(1,ir)-rf)**2+(zs(1,ir)-zf)**2)
c
              if (ctargopt.eq.5) then
                 S      = KBFS(1,IR) * DELTAL1
              elseif (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >                ctargopt.eq.3.or.ctargopt.eq.6) then
                 S      = KBFST(IR,2) * DELTAL1
              endif
              P      = DELTAL1
c
              KSS2(1,IR) = S
              KPS2(1,IR) = P
              S = S + DELTAL2 *KBFS(1,IR)
              P = P + DELTAL2
              KSB(1,IR) = S
              kpb(1,ir) = P
           ENDIF
c
           ik =1
c
           if (cprint.eq.3.or.cprint.eq.9)
     >        write(6,'(2(i4),4(2x,f12.4))') ir,ik,kss2(ik,ir),
     >                   kps2(ik,ir),ksb(ik,ir),kss(ik,ir)
c
c          -------------------------------------
c
c          IK = 2, nks(ir)-1
c
c          Calculate bulk of ring
c
           DO 324 IK = 2, nks(ir)-1
              in = korpg(ik,ir)
              rb = (rvertp(1,in) + rvertp(2,in)) /2.0
              rf = (rvertp(3,in) + rvertp(4,in)) /2.0
              zb = (zvertp(1,in) + zvertp(2,in)) /2.0
              zf = (zvertp(3,in) + zvertp(4,in)) /2.0
c
              DELTAL1 = SQRT ((RS(IK,IR)-rb) ** 2 +
     >                   (ZS(IK,IR)-zb) ** 2)
              DELTAL2 = SQRT ((RS(IK,IR)-rf) ** 2 +
     >                   (ZS(IK,IR)-zf) ** 2)
c
              S      = S + KBFS(IK,IR) * DELTAL1
              P      = P + DELTAL1
              KSS2(IK,IR) = S
              KPS2(IK,IR) = P
              S = S + DELTAL2 *KBFS(ik,IR)
              P = P + DELTAL2
              KSB(ik,IR) = S
              kpb(ik,ir) = P
c
              krb(ik,ir) = rf
              kzb(ik,ir) = zf
c
           if (cprint.eq.3.or.cprint.eq.9)
     >         write(6,'(2(i4),4(2x,f12.4))') ir,ik,kss2(ik,ir),
     >                        kps2(ik,ir),kss(ik,ir),ksb(ik,ir)
c
  324      CONTINUE
C
c          -------------------------------------
c
c          IK = nks(ir)
c 
C          Calculate end of rings
c
           in = korpg(nks(ir),ir)
           rb = (rvertp(1,in) + rvertp(2,in)) /2.0
           rf = (rvertp(3,in) + rvertp(4,in)) /2.0
           zb = (zvertp(1,in) + zvertp(2,in)) /2.0
           zf = (zvertp(3,in) + zvertp(4,in)) /2.0
c
           krb(nks(ir),ir) = rf
           kzb(nks(ir),ir) = zf
c
           IF (CTARGOPT.EQ.0.or.ctargopt.eq.4.or.ir.lt.irsep) then
c
              DELTAL1 = sqrt((rs(nks(ir),ir)-rb)**2+
     >                    (zs(nks(ir),ir)-zb)**2)
              DELTAL2 = 0.0
c
              S = S + DELTAL1 *KBFS(nks(ir),IR)
              P = P + DELTAL1
c
              KSS2(nks(ir),IR) = S
              KPS2(nks(ir),IR) = P
              ksb(nks(ir),ir) = S
              kpb(nks(ir),ir) = P
c
              KSMAXS2(IR) = S
              KPMAXS2(IR) = P
c
           ELSEIF (CTARGOPT.EQ.1.or.ctargopt.eq.2.or.
     >          ctargopt.eq.3.or.ctargopt.eq.5.or.ctargopt.eq.6) THEN
c
             DELTAL1 = SQRT ((RS(nks(ir),IR)-rb) ** 2 +
     >                   (ZS(nks(ir),IR)-zb) ** 2)
             DELTAL2 = SQRT ((RP(IDDS(IR,1))-RS(NKS(IR),IR)) ** 2 +
     >                   (ZP(IDDS(IR,1))-ZS(NKS(IR),IR)) ** 2)
c
             S      = S + KBFS(nks(ir),IR) * DELTAL1
             P      = P + DELTAL1
             KSS2(IK,IR) = S
             KPS2(IK,IR) = P
c
             if (ctargopt.eq.5) then
                S      = S + KBFS(NKS(IR),IR) * DELTAL2
             elseif (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >                ctargopt.eq.3.or.ctargopt.eq.6) then
                S      = S + KBFST(IR,1) * DELTAL2
             endif
c
             P = P + DELTAL2
c
             KSB(nks(ir),IR) = S
             KPB(nks(ir),IR) = P
             KSMAXS2(IR) = S
             KPMAXS2(IR) = P
           ENDIF
c
           ik = nks(ir)
c
           if (cprint.eq.3.or.cprint.eq.9) then
c
              write(6,'(2(i4),4(2x,f12.4))') ir,ik,kss2(ik,ir),
     >                 kps2(ik,ir),kss(ik,ir),ksb(ik,ir)
c
              write(6,'(a,i4,4(2x,f12.4))') 'Total:',ir,
     >           ksmaxs2(ir),kpmaxs2(ir),ksmaxs(ir),kpmaxs(ir)
c
           endif
c
c
c       For rings without polygons defined - KSS2 -> KSS and the
c       KSB values -> (KSS(IK+1,ir)-KSS(ik,IR))/2.0
c
        elseif (in.eq.0) then
c
           ksb(0,ir) = 0.0
           kpb(0,ir) = 0.0
c
           if (ir.ge.irsep) then
              krb(0,ir) = rp(idds(ir,2))
              kzb(0,ir) = zp(idds(ir,2))
           else
              krb(0,ir) = (rs(1,ir) + rs(nks(ir)-1,ir))/2.0
              kzb(0,ir) = (zs(1,ir) + zs(nks(ir)-1,ir))/2.0
           endif
c
c          Loop through ring
c
           do ik = 1,nks(ir)
              kss2(ik,ir) = kss(ik,ir)
              kps2(ik,ir) = kps(ik,ir)
c
c             Set cell boundaries for the ring without polygons.
c
              if (ik.eq.nks(ir)) then
                 ksb(ik,ir) = ksmaxs(ir)
                 kpb(ik,ir) = kpmaxs(ir)
                 if (ir.ge.irsep) then
                    krb(ik,ir) = rp(idds(ir,1))
                    kzb(ik,ir) = zp(idds(ir,1))
                 else
                    krb(ik,ir) = krb(0,ir)
                    kzb(ik,ir) = kzb(0,ir)
                 endif
              else
                 ksb(ik,ir) = (kss(ik+1,ir) + kss(ik,ir)) / 2.0
                 kpb(ik,ir) = (kps(ik+1,ir) + kps(ik,ir)) / 2.0
                 krb(ik,ir) = (rs(ik+1,ir) +rs(ik,ir)) / 2.0
                 kzb(ik,ir) = (zs(ik+1,ir) +zs(ik,ir)) / 2.0
              endif
c
c             Print
c
           if (cprint.eq.3.or.cprint.eq.9)
     >        write(6,'(2(i4),4(2x,f12.4))') ir,ik,kss2(ik,ir),
     >                    kps2(ik,ir),kss(ik,ir),ksb(ik,ir)
           end do
c
           ksmaxs2(ir) = ksmaxs(ir)
           kpmaxs2(ir) = kpmaxs(ir)

c
           if (cprint.eq.3.or.cprint.eq.9)
     >         write(6,'(a,i4,4(2x,f12.4))') 'Total:',ir,ksmaxs2(ir),
     >                           kpmaxs2(ir),ksmaxs(ir),kpmaxs(ir)
c
        endif

C
  322 CONTINUE
c
c     Assign the KSS2 and KPS2 values to KSS and KPS if the PDOPT=1 has
c     been specified.
c
c     PDOPT=1 - PARALLEL DISTANCE OPTION -------------------------------
c
         if (pdopt.eq.1) then
c
            do ir = 1,nrs
c
               ksmaxs(ir) = ksmaxs2(ir)
               kpmaxs(ir) = kpmaxs2(ir)
c
               do ik = 1,nks(ir)
                  kss(ik,ir) = kss2(ik,ir)
                  kps(ik,ir) = kps2(ik,ir)
c
               end do
c
            end do
c
c
         endif
c
c     End of grid option IF for calculating KSS2 ...
c
      endif
      CALL OUTPUTDATA(87,'AFTER KKS')
c
C-----------------------------------------------------------------------
c
c     Calculate the S-reflection positions if the parallel ion
c     reflection option is turned on for a broken grid.
c
C-----------------------------------------------------------------------
c
      if (s_reflect_opt.ne.0) then 
         call calc_s_reflect
      endif 
c
C-----------------------------------------------------------------------
c
c     Enclose the following in an IF-statment so it is only executed if
c     the diffusive pinch correction is turned on.
c
C-----------------------------------------------------------------------
c
c      if (cdiffopt.ne.0) then
c
c     Calculate the individual segment size ratios for use in another
c     option.
c
c     For grids with polygon information
c
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.
     >    cgridopt.EQ.LINEAR_GRID) then
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               kpsiz(ik,ir) = kpb(ik,ir) - kpb(ik-1,ir)
c
               if (cprint.eq.3.or.cprint.eq.9) then
                  write (6,'(a,2i6,4(1x,g15.8))') 'kpsiz:',
     >              ik,ir,kpsiz(ik,ir),
     >              ksb(ik,ir)-ksb(ik-1,ir),
     >              (ksb(ik,ir)-ksb(ik-1,ir))/kpsiz(ik,ir)
               endif
c
            end do
         end do
c
c     For grids without polygon information
c
      else

         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (ik.eq.1) then
                  deltal = (kps(ik+1,ir) - kps(ik,ir))/2.0
     >              +    SQRT ((RS(ik,IR)-RS(nks(ir)-1,IR)) ** 2 +
     >                         (ZS(ik,IR)-ZS(nks(ir)-1,IR)) ** 2)/2.0
               elseif (ik.eq.nks(ir)) then
                  deltal = (kps(ik,ir) - kps(ik-1,ir))/2.0
     >              +    SQRT ((RS(2,IR)-RS(ik,IR)) ** 2 +
     >                         (ZS(2,IR)-ZS(ik,IR)) ** 2)/2.0
               else
                  deltal = (kps(ik+1,ir) - kps(ik,ir))/2.0 +
     >                  (kps(ik,ir) - kps(ik-1,ir))/2.0
               endif
c
               kpsiz (ik,ir) = deltal
c
            end do
         end do
c
      endif
c
c     Calculate required values for perpendicular step options 1 and 2
c
      if (cdiffopt.eq.1.or.cdiffopt.eq.2) then
c
c
c     Calculate the ratios of segment lengths for each cell.
c
c     These are done here in order to save calculating them at every
c     diffusive step and so save some CPU.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
            ikin  = ikins(ik,ir)
            irin  = irins(ik,ir)
c
            if (ir.eq.1) then
c
               kprat2(ik,ir,1) = -1.0
c
               if (kpsiz(ik,ir).eq.kpsiz(ikout,irout)) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = kpsiz(ik,ir)/
     >                        (kpsiz(ik,ir)-kpsiz(ikout,irout))
               endif
c
            elseif (ir.eq.irwall) then
c
c              Inner ratio
c
               if (kpsiz(ikin,irin).eq.kpsiz(ik,ir)) then
                  kprat2(ik,ir,1) = HI
               else
                  kprat2(ik,ir,1) = kpsiz(ik,ir)/
     >                       (kpsiz(ikin,irin)-kpsiz(ik,ir))
               endif
c
c              Outer ratio
c
               kprat2(ik,ir,2) = -1.0
c
c
            elseif (ir.eq.irtrap) then
c
c              Inward ratio
c
               kprat2(ik,ir,1) = -1.0
c
c              Outward ratio
c
               if (kpsiz(ik,ir).eq.kpsiz(ikout,irout)) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = kpsiz(ik,ir)/
     >                        (kpsiz(ik,ir)-kpsiz(ikout,irout))
               endif
c
            else
c
c              Inner ratio
c
               if (kpsiz(ikin,irin).eq.kpsiz(ik,ir)) then
                  kprat2(ik,ir,1) = HI
               else
                  kprat2(ik,ir,1) = kpsiz(ik,ir)/
     >                       (kpsiz(ikin,irin)-kpsiz(ik,ir))
               endif
c
c              Outer ratio
c
               if (kpsiz(ik,ir).eq.kpsiz(ikout,irout)) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = kpsiz(ik,ir)/
     >                        (kpsiz(ik,ir)-kpsiz(ikout,irout))
               endif
c
            endif

            if (cprint.eq.3.or.cprint.eq.9)
     >         write (6,'(a,2i4,2g16.8)') 'kprat2:',ir,ik,
     >                 kprat2(ik,ir,1),kprat2(ik,ir,2)
c
            end do
         end do
c
c     Calculate required values for perpendicular step option 3
c
      elseif (cdiffopt.eq.3) then
c
c        Calculate the kprat2 quantity for each 1/2 cell using cell
c        boundary lengths as well as the cell centre size.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
            ikout = ikouts(ik,ir)
            irout = irouts(ik,ir)
            ikin  = ikins(ik,ir)
            irin  = irins(ik,ir)
c
            if (ir.eq.1) then
c
               kprat2(ik,ir,1) = -1.0
c
               if (kpsiz(ik,ir).eq.kpsiz(ikout,irout)) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = kpsiz(ik,ir)/
     >                        (kpsiz(ik,ir)-kpsiz(ikout,irout))
               endif
c
            elseif (ir.eq.irwall) then
c
c              Inner ratio
c
               if (kpsiz(ikin,irin).eq.kpsiz(ik,ir)) then
                  kprat2(ik,ir,1) = HI
               else
                  kprat2(ik,ir,1) = kpsiz(ik,ir)/
     >                       (kpsiz(ikin,irin)-kpsiz(ik,ir))
               endif
c
c              Outer ratio
c
               kprat2(ik,ir,2) = -1.0
c
c
            elseif (ir.eq.irtrap) then
c
c              Inward ratio
c
               kprat2(ik,ir,1) = -1.0
c
c              Outward ratio
c
               if (kpsiz(ik,ir).eq.kpsiz(ikout,irout)) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = kpsiz(ik,ir)/
     >                        (kpsiz(ik,ir)-kpsiz(ikout,irout))
               endif
c
            else
c
c              Inner ratio
c
               inlen = polysidelen(ik,ir,INWARD41,rc)
c
               if (rc.eq.1)
     >             write (6,*) 'Error: POLYSIDELEN : ',ik,ir,1,rc
c
               cenlen= kpsiz(ik,ir)

               outlen = polysidelen(ik,ir,OUTWARD23,rc)
c
               if (rc.eq.1)
     >             write (6,*) 'Error: POLYSIDELEN : ',ik,ir,1,rc
c
               if (inlen.eq.cenlen) then
                  kprat2(ik,ir,1) = HI
               else
                  kprat2(ik,ir,1) = cenlen/
     >                       (inlen-cenlen)
               endif
c
c              Outer ratio
c
               if (cenlen.eq.outlen) then
                  kprat2(ik,ir,2) = HI
               else
                  kprat2(ik,ir,2) = cenlen/
     >                        (cenlen-outlen)
               endif
c
            endif

            if (cprint.eq.3.or.cprint.eq.9) then
               write (6,'(a,2i4,5g16.8)') 'kprat2:',ir,ik,
     >                 kprat2(ik,ir,1),kprat2(ik,ir,2)
c     >                 ,inlen,cenlen,outlen
            endif
c
c           End of ik and ir do loops
c
            end do
         end do

c
      endif
c
c
c     Endif for CDIFFOPT
c
c      endif
c
c
C-----------------------------------------------------------------------
C     CALCULATE DISTANCES FORWARD AND BACKWARD ALONG THE FIELD LINES
c     BETWEEN CELL CENTRES
C-----------------------------------------------------------------------
c
C
      DO 360 IR = 1, NRS
        IF (IR.GE.IRSEP) THEN
          IF (CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4) THEN
            KBACDS(1,IR) = 0.0
          ELSE
            KBACDS(1,IR) = KSS(1,IR)
          ENDIF
        ELSE
          KBACDS(1,IR) = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
        ENDIF
        KFORDS(1,IR) = KSS(2,IR) - KSS(1,IR)
        DO 350 IK = 2, NKS(IR)-1
          KBACDS(IK,IR) = KSS(IK,IR) - KSS(IK-1,IR)
          KFORDS(IK,IR) = KSS(IK+1,IR) - KSS(IK,IR)
  350   CONTINUE
        KBACDS(NKS(IR),IR) = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
        IF (IR.GE.IRSEP) THEN
          IF (CTARGOPT.EQ.0 .OR. CTARGOPT.EQ.4) THEN
            KFORDS(NKS(IR),IR) = 0.0
          ELSE
            KFORDS(NKS(IR),IR) = KSMAXS(IR) - KSS(NKS(IR),IR)
          ENDIF
        ELSE
          KFORDS(NKS(IR),IR) = KSS(2,IR) - KSS(1,IR)
        ENDIF
  360 CONTINUE
c
c
c---------------------------------------------------------------
c
c     Calculate the cell containing the mid-point of each ring
c     - based on the S-distance between the two targets.
c     - and based on the distance between cell centres.
c
c     CAUTION: This array contains the index of the cell center
c              immediately less than the S-centre value of the
c              ring. SOME DIVIMP code uses this as the
c              mid-point index, other DIVIMP code requires that
c              IKMID = IKMIDS(IR) + 1 - be aware of this when
c              modifying or adding code - check the definition
c              of IKMID being used in the loacl routine.
c
      call izero(ikmids,maxnrs)
c
      do ir = 1,nrs
         smid = ksmaxs(ir)/2.0
         do ik = 1 ,nks(ir)
            if (kss(ik,ir).gt.smid) then
               ikmids(ir) = ik-1
               goto 365
            endif
         end do
c
 365     continue
c
c         write (6,'(a,3i5,3g14.5)')
c     >            'IKMIDS:',ir,ik,ikmids(ir),smid,kss(ikmids(ir),ir),
c     >                          kss(ikmids(ir)-1,ir)
c
      end do

c
C-----------------------------------------------------------------------
c     Calculate the R,Z values of the inner and outer mid-planes
C-----------------------------------------------------------------------
c
      call find_midplane
C
C-----------------------------------------------------------------------
C     CALCULATE ELEMENTAL VOLUMES AND AREAS
C-----------------------------------------------------------------------
C
      CALL TAUVOL
c
C-----------------------------------------------------------------------
c     Write out the grid - formatted and including vertices.
C-----------------------------------------------------------------------
c
      if ((cprint.eq.6.or.cprint.eq.9)
     >    .and.(cgridopt.eq.0.or.cgridopt.eq.3.or.
     >          cgridopt.eq.LINEAR_GRID)) then
         call writegrd(cgridopt)
      endif
C
C-----------------------------------------------------------------------
C     INITIALIZE 3D B FIELD VECTORS IN THE bfield module
C-----------------------------------------------------------------------
c
      call setup_bvectors
c
C-----------------------------------------------------------------------
c     Calculate the separatrix area and effective area
C-----------------------------------------------------------------------
c
      call calc_asep_eff
c
      call print_average_sepdist 
c
C-----------------------------------------------------------------------
c     Calculate the flux based EDGE2D target conditions if the
c     information is available to do so ...
C-----------------------------------------------------------------------
c
c     Calculate the flux based target density IF flux data available.
c
      if ((cre2d.eq.1.or.cre2d.eq.2).and.fluxpts.gt.0.0) then

         do ir = irsep,nrs

c
c           Find ring for EDGE2D data
c
            in = ringno(ir,fluxinfo,fluxpts,maxins,4,ierr)
c
c           Inner target
c
c           Load actual flux
c
            if (dds(idds(ir,1)).ne.0.0.and.
     >          costet(idds(ir,1)).ne.0.0.and.
     >          rp(idds(ir,1)).ne.0.0) then
c
               e2dtarg(ir,6,1) =  fluxinfo(in,2)  /
     >           (dds(idds(ir,1))* 2.0 * PI * rp(idds(ir,1)))
     >            *kbfst(ir,1)/costet(idds(ir,1))
c
            else

               e2dtarg(ir,6,1) = 0.0

            endif
c
c           Calculate density - assuming E2D velocity
c
            if (e2dtarg(ir,4,1).ne.0.0) then
               e2dtarg(ir,7,1) = e2dtarg(ir,6,1) / e2dtarg(ir,4,1)
            else
               e2dtarg(ir,7,1) = 0.0
            endif
c
c           Calculate Velocity - assuming E2D density
c
            if (e2dtarg(ir,1,1).ne.0.0) then
               e2dtarg(ir,8,1) = e2dtarg(ir,6,1) / e2dtarg(ir,1,1)
            else
               e2dtarg(ir,8,1) = 0.0
            endif
c
c           Outer target
c
c           Load actual flux
c
            if (dds(idds(ir,2)).ne.0.0.and.
     >          costet(idds(ir,2)).ne.0.0.and.
     >          rp(idds(ir,2)).ne.0.0) then
c
               e2dtarg(ir,6,2) =   fluxinfo(in,3)  /
     >              (dds(idds(ir,2))* 2.0 * PI * rp(idds(ir,2)))
     >              *kbfst(ir,2)/costet(idds(ir,2))
c
            else

               e2dtarg(ir,6,2) = 0.0

            endif
c
c           Calculate density - assuming E2D velocity
c
            if (e2dtarg(ir,4,2).ne.0.0) then
               e2dtarg(ir,7,2) = e2dtarg(ir,6,2) / e2dtarg(ir,4,2)
            else
               e2dtarg(ir,7,2) = 0.0
            endif
c
c           Calculate Velocity - assuming E2D density
c
            if (e2dtarg(ir,1,2).ne.0.0) then
               e2dtarg(ir,8,2) = e2dtarg(ir,6,2) / e2dtarg(ir,1,2)
            else
               e2dtarg(ir,8,2) = 0.0
            endif
c
            write (6,'(a,2i4,8(1x,g12.5))') 'E2DTARG:FLUX ',ir,in,
     >               e2dtarg(ir,7,1),e2dtarg(ir,6,1),e2dtarg(ir,8,1),
     >               e2dtarg(ir,8,1)/e2dtarg(ir,4,1),
     >               e2dtarg(ir,7,2),
     >               e2dtarg(ir,6,2),e2dtarg(ir,8,2),
     >               e2dtarg(ir,8,2)/e2dtarg(ir,4,2)
c
            end do
c
      endif
c slmod begin
      if (northopt.eq.3) then
         IF (stopopt2.EQ.900) THEN
           CALL CalcMetricQuickandDirty
c         ELSEIF (cgridopt.EQ.LINEAR_GRID) THEN
c           WRITE(0,*) 'NOT CALCULATING THETAG METRIC'
c           thetag = 1.0
         ELSEIF (.TRUE..OR.grdnmod.NE.0) THEN
           IF (grdnmod.EQ.0) THEN
             WRITE(0,*) 'WARNING tauin1: Calling the new THETAG '
             WRITE(0,*) '  calculation code for all cases as of '
             WRITE(0,*) '  10/03/2010'
           ENDIF
           CALL CalcMetric
         ELSEIF (connected) THEN
           WRITE(0,*) 'NOT CALCULATING THETAG METRIC'
           thetag = 1.0
         ELSE
c... Get rid of this after validating that the new routines give the same
c    result:
           WRITE(0,*)
           WRITE(0,*) '*** CALLING THE OLD THETA CALCULATION '//
     .                'CODE ***'
           WRITE(0,*)
           call calctheta_old
         ENDIF
      endif

c...  Make a token call to SetupRelaxation, in case SOL22/24 is not
c     being used.  This should be moved somewhere more appropriate,
c     but not sure where just yet:
      IF (relmode.NE.0) CALL SetupRelaxation

c...  Load supplimental .RAW file(s) if requested:
      IF (osm_store.GE.0) CALL LoadPIN


c...  Read in the results from a SOLPS/B2 simulation (from Rozhansky at
c     the moment): 
      IF (solps_opt.GT.0) CALL LoadSOLPSData

c slmod end
c
C-----------------------------------------------------------------------
c     CALCULATE BACKGROUND PLASMA
C-----------------------------------------------------------------------
c
c     Call the routine to calculate the Background Plasma
c
c     The following routine will calculate the background plasma over
c     the entire grid based on the options selected. This includes
c     combining various different SOL and plasma options.
c
C-----------------------------------------------------------------------
c
      call bgplasma(title,equil)
c slmod begin - temp
      CALL DB('DONE CALCULATING BACKGROUND PLASMA')
c slmod end
C
C-----------------------------------------------------------------------
c
c     Write a DIVIMP plasma background file - intended for internal use
c
c     Use the DIVIMP Print option flag to print the background plasma
c     file. I would not expect other debugging and background plasma file
c     writing to be done at the same time. The purpose of this option
c     is to save some time when one wants to use the same identical
c     DIVIMP plasma for a series of cases. Especially when PIN execution
c     is involved, these calculations can take a long time. This allows
c     them to be done once and the data copied for later use.
c
c     Since this file isn't too large - switch to writing it all the
c     time
c
c      if (cprint.eq.10) then
c
         call wrtdivbg
c
c      endif
c
C-----------------------------------------------------------------------
C
c     If a PIN wall has been specified for NIMBUS - it is necessary
c     to recalculate all of the NEUTRAL Wall at this point
c     Prior to any Wall launch probability calculations.
c
c     Eirene always runs on the wall supplied by DIVIMP and so does
c     not have the problems associated with the NIMBUS wall definition.
c
c     However, if the wall option has been specified as 5 then run 
c     the remap code for both.
c
      if (cpinopt.gt.0.and.pincode.eq.0.and.
     >   (tmpcneur.eq.5.or.tmpctrap.eq.5)) then
c
c        If the wall has been redefined to INCLUDE the baffles then
c        call the routine that will apply this to the vesm related arrays
c        and reorganize and recalculate where necessary teh fluxes in the
c        PIN flux arrays.
c
         if (nbufle.gt.0.or.nbufmx.gt.0) then

            call redef_pinwalldata

         endif
c
c        Map the NIMBUS segments onto the WALLCO array.
c
         if (tmpcneur.eq.5) then
c
            call nimwall
c
         endif
c
c        Map the NIMBUS segments onto the WALLCO2 array.
c
         if (tmpctrap.eq.5) then
c
            call nimwall2
c
         endif
c
c        Redo the neutral walls
c
         if (tmpcneur.eq.5)   cneur = 2
         if (tmpctrap.eq.5)   ctrap = 3
c
         call dowall
c
c        Set up the Indices (pointers) into the NIMBUS/EIRNE
c        arrays - for all of the Wallpts and for the
c        target segments.
c
         call nimind
c
c        Reset the wall options
c
         cneur = tmpcneur
         ctrap = tmpctrap
c
c     For Eirene99 - assume the wall definitions match and run nimind
c     to generate the appropriate index relationships 
c
      elseif (cpinopt.gt.0.and.
     .        (pincode.eq.2.or.pincode.eq.3.or.
     .         pincode.eq.4.or.pincode.eq.5)) then 
c
c        Reset the wall definition variables to original values
c
         if (tmpcneur.eq.5) cneur = tmpcneur
         if (tmpctrap.eq.5) ctrap = tmpctrap
c
c        Set up the Indices (pointers) into the NIMBUS/EIRENE
c        arrays - for all of the Wallpts and for the
c        target segments.
c
         call nimind
c
      endif
c
C-----------------------------------------------------------------------
c
c     Deal with UEDGE background plasma issues - set up wall fluxes if
c     necessary.
c
      if (uedge_bg.eq.1.and.(cneur.eq.7.or.ctrap.eq.7)) then
         call setup_uedge_wall
      endif
c
c     Check fluid code carbon and hydrogen fluxes
c
      if (uedge_bg.eq.1) then
         call check_fluxes
      endif
c
C-----------------------------------------------------------------------
c
c
c     After the walls and targets have been finalized - assign
c     the wall and target temperature to the individual segments.
c
      call dotemp
c
C
C-----------------------------------------------------------------------
C     AFTER PIN HAS BEEN EXECUTED - IT IS NECESSARY TO FINISH THE WALL
C     DEFINITIONS BY CALCULATING THE DISTRIBUTED WALL LAUNCH
C     PROBABILITY.  THIS IS DELAYED UNTIL AFTER A POSSIBLE PIN EXECUTION
C     BECAUSE THE CALCULATION MAY REQUIRE THE IMPURITY SPUTTERING DATA
C     GENERATED BY PIN
C-----------------------------------------------------------------------
C
      CALL CALCWP
c
C-----------------------------------------------------------------------
c
c  This is done after the background plasma is calculated - ideally this
c  might be best at the very end of TAU - but this requires checking the
c  assumptions used in any kvhs references in the following code.
C
C     MULTIPLY BY ITERATION TIME FACTORS
C
c     Properly normalize the velocity and electric field.
c
      VHFACT = QTIM
      EFACT  = QTIM * QTIM * EMI / CRMI
      DO 9800 IR = 1,NRS
        DO 9800 IK = 1,NKS(IR)
          KVHS(IK,IR) = KVHS(IK,IR) * VHFACT
          KES(IK,IR)  = KES(IK,IR)  * EFACT
9800  CONTINUE
c
C-----------------------------------------------------------------------
c
      VFLUID = KVHS(1,IRSEP)
      IF (CTEB0.LE.0.0) CTEB0 = KTEBS(1,IRSEP)
      WRITE (6,'(/,'' FLUID VELOCITY    VFLUID'',1P,E12.4)') VFLUID
C
C-----------------------------------------------------------------------
C     TEMPERATURE GRADIENT FORCES IN THE SOL
C-----------------------------------------------------------------------
C
      DO 660 IZ = 1, NIZS
        IF     (CIOPTM.EQ.0) THEN
          KALPHS(IZ) = 0.0
        ELSEIF (CIOPTM.EQ.1.or.cioptm.eq.3) THEN
          KALPHS(IZ) = 0.71 * REAL(IZ*IZ)
        ELSEIF (CIOPTM.EQ.2) THEN
          KALPHS(IZ) = 1.5*(1.0-0.6934*(1.3167**(-IZ)))*REAL(IZ*IZ)
        ENDIF
C
        IF     (CIOPTN.EQ.0) THEN
          KBETAS(IZ) = 0.0
        ELSEIF (CIOPTN.EQ.1.or.cioptn.eq.3) THEN
          MU = CRMI / (CRMI+CRMB)
          KBETAS(IZ)=-3.0*(1.0-MU-5.0*ROOT2*(1.1*MU**2.5-0.35*MU**1.5)*
     >      REAL(IZ*IZ)) / (2.6 - 2.0*MU + 5.4*MU*MU)
        ELSEIF (CIOPTN.EQ.2) THEN
          KBETAS(IZ)= CHZO * REAL(IZ*IZ) /
     >                (REAL(CZO) + SQRT(0.5*(1.0+CRMB/CRMI)))
        ENDIF
  660 CONTINUE
C
      FACT = QTIM * QTIM * EMI / CRMI
      CALL RZERO (KFEGS, MAXNKS*MAXNRS)
      CALL RZERO (KFIGS, MAXNKS*MAXNRS)
c
c     Do SOL and private plasma
c
      DO 680 IR = IRSEP, NRS
        SMAX = KSMAXS(IR)
        KFEGS(1,IR) = FACT * (KTEBS(2,IR)-KTEBS(1,IR))/KFORDS(1,IR)
        KFIGS(1,IR) = FACT * (KTIBS(2,IR)-KTIBS(1,IR))/KFORDS(1,IR)
        if (kss(1,ir).eq.0) then
           kfeds(idds(ir,2)) = kfegs(1,ir)
           kfids(idds(ir,2)) = kfigs(1,ir)
        else
           kfeds(idds(ir,2)) = fact * (ktebs(1,ir)-kteds(idds(ir,2)))
     >                             /kss(1,ir)
           kfids(idds(ir,2)) = fact * (ktibs(1,ir)-ktids(idds(ir,2)))
     >                             /kss(1,ir)
        endif
        DO 670 IK = 2,nks(ir)-1
c
          KFEGS(IK,IR) = ((KTEBS(IK,IR)-KTEBS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KTEBS(IK+1,IR)-KTEBS(IK,IR))/KFORDS(IK,IR)) * FACT * 0.5
          KFIGS(IK,IR) = ((KTIBS(IK,IR)-KTIBS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KTIBS(IK+1,IR)-KTIBS(IK,IR))/KFORDS(IK,IR)) * FACT * 0.5
c
  670   CONTINUE
        KFEGS(NKS(IR),IR) = FACT *
     >    (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))/KBACDS(NKS(IR),IR)
        KFIGS(NKS(IR),IR) = FACT *
     >    (KTIBS(NKS(IR),IR)-KTIBS(NKS(IR)-1,IR))/KBACDS(NKS(IR),IR)
        if (kss(nks(ir),ir).eq.smax) then
           kfeds(idds(ir,1)) = kfegs(nks(ir),ir)
           kfids(idds(ir,1)) = kfigs(nks(ir),ir)
        else
           kfeds(idds(ir,1))=fact*(ktebs(nks(ir),ir)-kteds(idds(ir,1)))
     >                             /(smax-kss(nks(ir),ir))
           kfids(idds(ir,1))=fact*(ktibs(nks(ir),ir)-ktids(idds(ir,1)))
     >                             /(smax-kss(nks(ir),ir))
        endif
c slmod begin - 04/01/2020
c
c       Setting the KFEGS and KFIGS to 0.0 creats a startling artifact
c       for higher charge states when the other parallel force terms
c       are small.  This nulling out of the thermal force at the midpoint
c       of each ring is a good idea when there's a temperature 
c       discontinuty, but it's not necessary for SOL28 (and a problem for
c       some ITER cases I found):
c
        if (getmodel(3,ir).eq.28) cycle
c        write(0,*) '****** blanking mid-point grad terms! ******'
c slmod end
c
c       Fix the mid-point of the ring where the solutions join and force
c       KFEGS and KFIGS to be zero for ikmid and ikmid+1 
c
        kfegs(ikmids(ir),ir) = 0.0 
        kfegs(ikmids(ir)+1,ir) = 0.0 
c
        kfigs(ikmids(ir),ir) = 0.0 
        kfigs(ikmids(ir)+1,ir) = 0.0 
c
  680 CONTINUE
c
c
C-----------------------------------------------------------------------
c
c     Do FORCES in core plasma.
c
      do ir = 1,irsep-1
c
        KFEGS(1,IR) = ((KTEBS(nks(ir),IR)-KTEBS(nks(ir)-1,IR))
     >                 /KBACDS(nks(ir),IR) +
     >      (KTEBS(2,IR)-KTEBS(1,IR))/KFORDS(1,IR)) * FACT * 0.5
        KFIGS(1,IR) = ((KTIBS(nks(ir),IR)-KTIBS(nks(ir)-1,IR))
     >                 /KBACDS(nks(ir),IR) +
     >      (KTIBS(2,IR)-KTIBS(1,IR))/KFORDS(1,IR)) * FACT * 0.5
c
        DO  IK = 2,nks(ir)-1
          KFEGS(IK,IR) = ((KTEBS(IK,IR)-KTEBS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KTEBS(IK+1,IR)-KTEBS(IK,IR))/KFORDS(IK,IR)) * FACT * 0.5
          KFIGS(IK,IR) = ((KTIBS(IK,IR)-KTIBS(IK-1,IR))/KBACDS(IK,IR) +
     >      (KTIBS(IK+1,IR)-KTIBS(IK,IR))/KFORDS(IK,IR)) * FACT * 0.5
c
c        write (6,'(a,2i4,6g18.10)') 'KFIGS:',ik,ir,
c     >         ktibs(ik-1,ir),ktibs(ik,ir),
c     >         ktibs(ik+1,ir),kfords(ik,ir),kbacds(ik,ir),kfigs(ik,ir)
c
        end do
c
        KFEGS(NKS(IR),IR) = kfegs(1,ir)
        KFIGS(NKS(IR),IR) = kfigs(1,ir)
c
      end do
c
C-----------------------------------------------------------------------
c
c     Calculate the Kinetic correction to the Ion and electron
c     temperature gradient forces.
c
C-----------------------------------------------------------------------
c
      call calcapp_fgradmod
c
c psmod
c
c-----------------------------------------------------------------------
c
C     Calculate the Background Velocity Gradient for each grid cell.
c
c-----------------------------------------------------------------------
c
      CALL VBGRAD 
c
c-----------------------------------------------------------------------
c
c     Calculate the Coefficients used in the Coulomb Collision 
c     subroutine provided by Dirk Reiser
c
c-----------------------------------------------------------------------
c
      CALL COEFF(NIZS)
c
c-----------------------------------------------------------------------
c
c     Compare gradient scale lengths of the temperature and velocity 
c     gradients with the minimum values for which Dirk Reiser's Drift-
c     kinetic Model is valid. 
c
c-----------------------------------------------------------------------
c
      IF (CIOPTR.EQ.1.OR.CIOPTR.EQ.2) CALL GradScaleLengthCheck
c
c psmod
c
C
C-----------------------------------------------------------------------
C     CROSS FIELD DIFFUSION AND PINCH
C-----------------------------------------------------------------------
C
      SPERP = SQRT (2.0 * CDPERP * QTIM)
      sperpt = SQRT (2.0 * CDPERPT * QTIM)
c
c     Assign reference cell.
c
      ikrefsol = nks(irsep)/2 + 1
      ikrefcore= ikins(ikrefsol,irsep)
      ikrefpp  = 1
      sdperpref= sperp
      sdperppp =sdperpref * (kinds(ikrefpp,irsep)*cosali(ikrefpp,irsep))
     >                  /(kinds(ikrefsol,irsep)*cosali(ikrefsol,irsep))
c
      if (cioptj.eq.0.or.cioptj.eq.1.or.cioptj.eq.2) then
         DO IR = 1, NRS
            DO IK = 1, NKS(IR)
               IF     (CIOPTJ.EQ.0.or.cioptj.eq.2.OR.IR.LT.IRSEP) THEN
                  if (ir.le.irwall) then
                     KPERPS(IK,IR) = SPERP
                  elseif (ir.ge.irtrap) then
                     KPERPS(IK,IR) = SPERPT
                  endif
               ELSEIF (CIOPTJ.EQ.1) THEN
                  KPERPS(IK,IR) = SPERP * KNBS(IK,IRSEP) / KNBS(IK,IR)
               ENDIF
            end do
         end do
c
      elseif (cioptj.eq.3.or.cioptj.eq.4) then
c
c        Reevaluate Dperp spatially based on Constant number of
c        cross-field steps for each cell.
c
         do ir = 1,irwall
c
            if (ir.lt.irsep) then
               ikr = ikrefcore
            else
               ikr = ikrefsol
            endif
c
            distref = distout(ikr,ir) + distin(ikr,ir)
c
            refstep = distref/sperp
c
            do ik = 1,nks(ir)
c
               distcur = distout(ik,ir) + distin(ik,ir)
c
               kperps(ik,ir) = distcur/refstep
c
            end do
c
         end do


         if (cioptj.eq.3) then
c
c           Loop to match what is done in Edge2D at this time.
c
            do ir = irtrap,nrs
               do ik = 1,nks(ir)
c
                  if (ik.gt.ikti) then
                     ikr = (ik -ikti -1) + ikto
                  else
                     ikr = ik
                  endif
c
                  kperps(ik,ir) = kperps(ikr,irsep)
c
               end do
            end do
c
         elseif (cioptj.eq.4) then
c
c           Try to calculate a Dperp variation for the private plasma
c           Assume that the Dperp in the first cells across the
c           target at the IK=1 index plate - are a constant - equal
c           to the value on the separatrix in cell 1.
c

            do ir = irtrap,nrs
c
               ikr = 1
c
               distref = distout(ikr,ir) + distin(ikr,ir)
c
               refstep = distref/kperps(ikr,irsep)
c
               do ik = 1,nks(ir)
c
                  distcur = distout(ik,ir) + distin(ik,ir)
c
                  kperps(ik,ir) = distcur/refstep
c
               end do
c
            end do
c
         endif
c
      endif
c
c-------------------------------------------------------------
c
c     PINCH FLOWS
c
c     Call routine to set up pinch flows based on input options
c
c
      call calculate_pinch


c slmod end
c
C
C-----------------------------------------------------------------------
C     ION REMOVAL LOSS TIMES
C-----------------------------------------------------------------------
C
C     QOVRT IS CALCULATED OUTSIDE THE LOOP AND THEN ASSIGNED
C     BECAUSE AN OPTIMIZATION BUG AT JET WOULD GENERATE A
C     DIVISION BY ZERO ERROR
C
C
      IF (TLOSS.GT.0.0) THEN
         QOVRT = QTIM / TLOSS
      ELSE
         QOVRT = 0.0
      ENDIF
      DO 600 IR       = 1 , NRS
         DO 610 IK    = 1 , NKS(IR)
            DO 620 IZ = 1 , NIZS
               KPLOS(IK,IR,IZ) = QOVRT
  620       CONTINUE
  610    CONTINUE
  600 CONTINUE
C
C-----------------------------------------------------------------------
C     IONISATION AND RECOMBINATION TIMES
C-----------------------------------------------------------------------
C---- NOTE:
C---- THE NORMAL ROUTE IS TO DISABLE IONISATION BEYOND LEVEL NIZS.
C---- A SPECIAL CASE IS ALLOWED WHEN CIOPTA=0 SUCH THAT IONISATION FROM
C---- LEVEL NIZS TO NIZS+1 IS ALLOWED, BUT WHEN IT OCCURS THE ION IS
C---- NO LONGER FOLLOWED
C
      IF (CIOPTA.EQ.0 .OR. CIOPTA.EQ.3 .OR. CIOPTA.EQ.4) THEN
        MIZS = MIN (CION, NIZS+1)
      ELSE
        MIZS = MIN (CION, NIZS)
      ENDIF
C
C---- SET IONISATION / E-I RECOMBINATION TIME INTERVALS    CFIZS,CFRCS
C
      CALL IZTAU (CRMI,crmb,CION,RIZB,CIOPTA,cprint)
C
C---- SET C-X TIMES         KFCXS
C
      CALL CXREC (NIZS,CION,CIOPTI,RIZB,CRMB,CVCX,
     >            CNHC,CNHO,CLAMHX,CLAMHY,cprint,cpinopt)
C
C---- SET PROBABILITY OF EITHER AN IONISATION OR A RECOMBINATION
C---- SET PROPORTION OF THESE WHICH WILL BE RECOMBINATIONS
C---- PREVENT ANY IONISATION BEYOND MAXIMUM LIMIT SPECIFIED IF REQUIRED
C
      IF (MIZS.GT.1) THEN
        DO 770 IZ = 1, MIZS-1
         DO 760 IR = 1, NRS
          DO 750 IK = 1, NKS(IR)
C---- START WITH CHARACTERISTIC TIMES(**-1)
            TAUCH = 1.0/KFIZS(IK,IR,IZ)
            TAURECTOT = 0.0
            IF (KFRCS(IK,IR,IZ) .GT. 0.0) THEN
              TAUCH = TAUCH + 1.0/KFRCS(IK,IR,IZ)
              TAURECTOT = TAURECTOT + 1.0/KFRCS(IK,IR,IZ)
            ENDIF
            IF (KFCXS(IK,IR,IZ). GT. 0.0) THEN
              TAUCH = TAUCH + 1.0/KFCXS(IK,IR,IZ)
              TAURECTOT = TAURECTOT + 1.0/KFCXS(IK,IR,IZ)
            ENDIF
            KPCHS(IK,IR,IZ) = 1.0-exp(-QTIM*TAUCH)
            KPRCS(IK,IR,IZ) = TAURECTOT/TAUCH
            KPCHS(IK,IR,IZ) = MIN (1.0, KPCHS(IK,IR,IZ))
  750     CONTINUE
  760    CONTINUE
  770   CONTINUE
      ENDIF
C
      IF (MIZS .LE. MIN(CION,NIZS)) THEN
        DO 790 IR = 1, NRS
          DO 780 IK = 1, NKS(IR)
            TAUCH = 0.0
            IF (KFRCS(IK,IR,MIZS).GT.0.0) THEN
              TAUCH = TAUCH + 1/KFRCS(IK,IR,MIZS)
            ENDIF
            IF (KFCXS(IK,IR,MIZS).GT.0.0) THEN
              TAUCH = TAUCH + 1/KFCXS(IK,IR,MIZS)
            ENDIF
            KPCHS(IK,IR,MIZS) = 1.0-exp(-QTIM*TAUCH)
            KPRCS(IK,IR,MIZS) = 1.0
            KPCHS(IK,IR,MIZS) = MIN (1.0, KPCHS(IK,IR,MIZS))
  780     CONTINUE
  790   CONTINUE
      ENDIF
c
C---- SET IONISATION PROBABILITIES FOR NEUT ...
C---- (SAVES REPEATED CALCULATION EVERY ITERATION)
C---- THEY ARE ALL MULTIPLIED BY THE "IONISATION RATE FACTOR" IRF
C---- WHICH TYPICALLY MIGHT BE 0.2 TO GIVE DEEPER IONISATION.
C
      DO 795 IR = 1, NRS
        DO 795 IK = 1, NKS(IR)
          KPCHS(IK,IR,0) = MIN (1.0, CIRF * FSRATE / KFIZS(IK,IR,0))
          kpizs(ik,ir) = kpchs(ik,ir,0)
c
c         The change of state probability array for neutrals now has
c         three components -
c
c         1) Ionization (as before)
c         2) Charge Exchange Collision (ionization - currently
c                                       not implemented) (CXC)
c         3) Momentum transfer collision (MTC)
c
c          kpcxc(ik,ir)   = 0.0
c
          kpmtc = mtc_prob(crmb,crmi,ik,ir,fsrate)


c          if (mtcopt.eq.0) then
c             kpmtc(ik,ir) = 0.0
c          elseif (mtcopt.eq.1) then
c
c             if (cpinopt.eq.1) then
c
c                kpmtc(ik,ir)   = 16.0 * crmb / (3.0 * ( crmb + crmi))
c     >           * ( kelighi * knbs(ik,ir) + kelighg * pinatom(ik,ir))
c     >           * fsrate
c
c             else
c
c                kpmtc(ik,ir)   = 16.0 * crmb / (3.0 * ( crmb + crmi))
c     >           * ( kelighi * knbs(ik,ir) + kelighg * e2datom(ik,ir))
c     >           * fsrate
c
c             endif
c
c          endif
c
          if (kpmtc(ik,ir).ne.0.0) then
c
             lmfp_mtc(ik,ir) = 9.79e3 * SQRT (ktibs(ik,ir)/CRMI)
     >                             /(kpmtc(ik,ir)/fsrate)
c
          endif
c
c          kpchs(ik,ir,0) = kpchs(ik,ir,0) + kpcxc(ik,ir) + kpmtc(ik,ir)
c
          kpchs(ik,ir,0) = kpchs(ik,ir,0) + kpmtc(ik,ir)
          kpchs(ik,ir,0) = min(1.0,kpchs(ik,ir,0))
c
  795 CONTINUE
c
      if (cprint.eq.3.or.cprint.eq.9) then
         write (6,*) 'MTC Mean free paths from UEDGE:'
c
         do ir = irsep,nrs
            do ik = 1,nks(ir)
               write (6,'(2i4,10(1x,g15.6))') ik,ir,
     >                   kpmtc(ik,ir),lmfp_mtc(ik,ir)
            end do
         end do
c
      endif
c
C-----------------------------------------------------------------------
c
c     NEUT2D_OPT =1
c
c     Set up the raw data of neutral influxes for the neut2d option
c     code - for option 1 - take the UEDGE recombination rates as
c     the original basis.
c
C-----------------------------------------------------------------------
c
      if (neut2d_opt.eq.1.or.
     >   ((cneutb.eq.5.or.cneuth.eq.5).and.cgridopt.eq.3)) then
c
         call rzero (neut2d_raw,maxnks*maxnrs)
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               neut2d_raw(ik,ir) = e2drec(ik,ir) + e2dcxrec(ik,ir)
            end do
         end do
c
      endif
c
c
C-----------------------------------------------------------------------
c
c     CALL CALC_MPS
c
c     Calculate an estimate of  the Magnetic Pre-sheath (MPS) thickness
c     for each target segment - this must be done after the BG plasma has
c     been calculated since the temperatures are required
c
c     Also calculates the actual magnetic field angle on the
c     target segment.
c
C-----------------------------------------------------------------------
c
      call calc_mps
c
C-----------------------------------------------------------------------
c
c     CALC_TARGFLUXDATA 
c
c     Calculate the particle and heat fluxes onto the target surfaces  
c     due to direct hydrogen effects - ion and neutral - it does not
c     include
c     the affect of either radiation or impurities at this time.
c
C-----------------------------------------------------------------------
c
      call calc_targfluxdata
      call calc_wallfluxdata 
c
c----------------------------------------------------------------------- 
c
c     FP_WALLDIST 
c
c     Calculate the distance from the outermost real cells on the grid
c     to the wall - this distance is calculated perpendicular to the 
c     flux tube from the cell's center - this is used in fpopt 5 to 
c     define the width of the SOL for each element of the grid at the
c     edge of the plasma. It defines the width of the far peripheral 
c     region.   
c
c----------------------------------------------------------------------- 
c
      if (fpopt.eq.5) then 
         call setup_fp
      endif 
c
c
C-----------------------------------------------------------------------
c
c     Associate some plasma conditions with each element of the wall
c
C-----------------------------------------------------------------------
c
      call assign_wall_plasma 
c
C
C-----------------------------------------------------------------------
C     PRINT TABLES OF RESULTS
C-----------------------------------------------------------------------
C     THIS CAUSES SOME PROBLEMS ON THE JET IBM BECAUSE THE
C     DATA RECORDS ARE TOO LONG. SO LEAVE IT COMMENTED OUT FOR
C     NOW.
C
C     WRITE(6,*) 'TABLE OF S VALUES: RING NO. THEN BY IK INDEX'
C     WRITE(6,*) 'RING  MAX IK  S VALUES ....'
C     DO 925 IR = 1,NRS
C        WRITE(6,9050) IR,NKS(IR),(KSS(IK,IR),IK=1,NKS(IR))
C925  CONTINUE
C
      DO 920 IR = 1, NRS
        WRITE (6,9002)
        DO 910 IK = 1, NKS(IR)
          C(1) = FACTOR (KES(IK,IR),7)
          C(2) = FACTOR (KVHS(IK,IR),8)
          C(3) = FACTOR (KFEGS(IK,IR)*KALPHS(1),11)
          C(4) = FACTOR (KFIGS(IK,IR)*KBETAS(1),11)
          WRITE (6,9003) IK,IR,RS(IK,IR),ZS(IK,IR),
     >      BTS(IK,IR),KTEBS(IK,IR),KTIBS(IK,IR),
     >      KNBS(IK,IR),C(1),C(2),KSS(IK,IR),KBFS(IK,IR),C(3),C(4)
  910   CONTINUE
  920 CONTINUE
C
      DO 940 IR = 1, NRS
        WRITE (9,9022)
        DO 930 IK = 1, NKS(IR)
          WRITE (9,9023) IK,IR,RS(IK,IR),ZS(IK,IR),
     >      IKINS(IK,IR),IRINS(IK,IR),IKOUTS(IK,IR),IROUTS(IK,IR),
     >      KINDS(IK,IR),KOUTDS(IK,IR),KBACDS(IK,IR),KFORDS(IK,IR),
     >      KAREAS(IK,IR),KINS(IK,IR),KPS(IK,IR)
  930   CONTINUE
  940 CONTINUE
C
      write (6,*) 'REACTION RATE SUMMARY:'
c
      DO 960 IR = IRSEP-1, NRS
        WRITE (6,9032) 1,1,1,1,1,NIZS,NIZS,NIZS,NIZS,NIZS
        DO 950 IK = 1, NKS(IR)
          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),
     >      KNHS(IK,IR),KNBS(IK,IR),KTEBS(IK,IR),
     >      KFIZS(IK,IR,0),KPCHS(IK,IR,0),kpmtc(ik,ir),
     >     (KFIZS(IK,IR,IZ),KFRCS(IK,IR,IZ),KFCXS(IK,IR,IZ),
     >      KPCHS(IK,IR,IZ),KPRCS(IK,IR,IZ),IZ=1,NIZS,MAX(1,NIZS-1))
  950   CONTINUE
  960 CONTINUE
c
c      write (6,*) 'STATE CHANGE SUMMARY:'
c
c      DO IR = IRSEP-1, NRS
c
c        DO IK = 1, NKS(IR)
c          WRITE (6,'(2i5,14(1x,g12.5))') IK,IR,
c     >        (kpchs(ik,ir,iz),iz=0,nizs)
c          WRITE (6,'(10x,14(1x,g12.5))')
c     >        (kprcs(ik,ir,iz),iz=0,nizs)
c        end do
c      end do
c
c     CHECK VALIDITY OF CHANGE OF STATE PROBABILITIES
c
      error_count = 0
      iwarn_count = 0
      nwarn_count = 0
      pwarn_count = 0
      miniontau = hi
      minneutau = hi
      minprob = 1.0d0/2.147483647d9
c
c     Loop through cells and charge states.
c
      do ir = 1,nrs
c
c        Check only real rings
c
         if (ir.ne.1.and.ir.ne.irwall.and.ir.ne.irtrap) then
c
         do ik = 1,nks(ir)
c            if (kpchs(ik,ir,iz).ge.1.0) then
c               error_count = error_count + 1
c            endif
c
            do iz = 1,nizs

               if (kpchs(ik,ir,iz).ge.1.0) then
                  error_count = error_count +1
                  write(6,*) 'Error Count:',error_count,ik,ir,iz,
     >                        kpchs(ik,ir,iz)
               endif


               if ((kfizs(ik,ir,iz).lt.qtim.and.
     >              kfizs(ik,ir,iz).gt.0.0).or.
     >             (kfrcs(ik,ir,iz).lt.qtim.and.
     >              kfrcs(ik,ir,iz).gt.0.0).or.
     >             (kfcxs(ik,ir,iz).lt.qtim.and.
     >              kfcxs(ik,ir,iz).gt.0.0)) then
c
                   iwarn_count = iwarn_count + 1
                   miniontau = min(miniontau,kfizs(ik,ir,iz))
                   miniontau = min(miniontau,kfrcs(ik,ir,iz))
                   miniontau = min(miniontau,kfcxs(ik,ir,iz))
c
               endif
c
               if (kpchs(ik,ir,iz).gt.0.0.and.
     >             kpchs(ik,ir,iz).lt.minprob) then
                   pwarn_count = pwarn_count +1
c
                   write (6,'(a,3i4,3(1x,g14.6))') 'PWARN:',
     >                        ik,ir,iz,
     >                        kpchs(ik,ir,iz),minprob


               endif
c
            end do
c
c           Neutrals
c

            if (kpchs(ik,ir,0).gt.0.0.and.
     >          kpchs(ik,ir,0).lt.minprob) then
                pwarn_count = pwarn_count +1
c
                write (6,'(a,3i4,3(1x,g14.6))') 'PWARN:',
     >                        ik,ir,0,
     >                        kpchs(ik,ir,0),minprob

            endif
c
c
            if (kfizs(ik,ir,0).lt.fsrate.and.
     >          kfizs(ik,ir,0).gt.0.0) then

                nwarn_count = nwarn_count + 1
                minneutau = min(minneutau,kfizs(ik,ir,0))

            endif

         end do
c
c     End of Ring number check
c
      endif
c
c
      end do
c
      if (error_count.gt.0.0) then
c
         write (0,*) 'ERROR: STATE CHANGE PROBABILITY IS 1.0 IN ',
     >                       error_count, ' CELLS'
         write (0,*) '       REDUCE SIMULATION TIME STEPS'
c         write (0,*) '       PROGRAM STOPPING'
c
         write (6,*) 'ERROR: STATE CHANGE PROBABILITY IS 1.0 IN ',
     >                       error_count, ' CELLS'
         write (6,*) '       REDUCE SIMULATION TIME STEPS'
c         write (6,*) '       PROGRAM STOPPING'
c
         write (comment,'(a,i7,a)')
     >                   'ERROR: STATE CHANGE PROBABILITY'//
     >                     ' IS 1.0 IN ',
     >                       error_count, ' CELLS'
         call prc (comment)
         call prc ('       REDUCE SIMULATION TIME STEPS')
c         call prc ('       PROGRAM STOPPING')
         call prb

c         stop
      endif
c
      if (iwarn_count.gt.0.0) then
c
         write (0,'(a,i7,a)') 'WARNING : ION TIME STEP IS LARGER THAN'//
     >                         ' CHARACTERISTIC TIMES IN ',
     >                       iwarn_count, ' CELLS'
         write (0,'(a)') '          REDUCE SIMULATION ION TIME STEP'
c
         write (6,'(a,i7,a)') 'WARNING : ION TIME STEP IS LARGER THAN'//
     >                         ' CHARACTERISTIC TIMES IN ',
     >                       iwarn_count, ' CELLS'
         write (6,'(a)') '          REDUCE SIMULATION ION TIME STEP'
c
         write (comment,'(a,i7,a)')
     >                 'WARNING : ION TIME STEP IS LARGER THAN'//
     >                         ' CHAR. TIMES IN ',
     >                       iwarn_count, ' CELLS'
         call prc (comment)
         call prc('          REDUCE SIMULATION ION TIME STEP')
         call prb
c
      endif
c
      if (nwarn_count.gt.0.0) then
c
         write (0,'(a,i7,a)') 'WARNING : NEUTRAL TIME STEP'//
     >                        ' IS LARGER THAN'//
     >                        ' CHARACTERISTIC TIMES IN ',
     >                       nwarn_count, ' CELLS'
         write (0,'(a)') '          REDUCE SIMULATION NEUTRAL TIME STEP'
c
         write (6,'(a,i7,a)') 'WARNING : NEUTRAL TIME STEP'//
     >                        ' IS LARGER THAN'//
     >                        ' CHARACTERISTIC TIMES IN ',
     >                       nwarn_count, ' CELLS'
         write (6,'(a)') '          REDUCE SIMULATION NEUTRAL TIME STEP'
c
         write (comment,'(a,i7,a)')
     >                   'WARNING : NEUTRAL TIME STEP IS LARGER'//
     >                         ' THAN CHAR. TIMES IN ',
     >                       nwarn_count, ' CELLS'
         call prc (comment)
         call prc('          REDUCE SIMULATION NEUTRAL TIME STEP')
         call prb
c
      endif
c
      if (pwarn_count.gt.0.0) then
c
         write (0,'(a,i7,a)') 'WARNING : CHANGE OF STATE PROBABILITY'//
     >               ' IS LESS THAN MIN. IN ',
     >                       pwarn_count, ' CASES'
c
         write (6,'(a,i7,a)') 'WARNING : CHANGE OF STATE PROBABILITY'//
     >               ' IS LESS THAN MIN. IN ',
     >                       pwarn_count, ' CASES'
c
         write (comment,'(a,i7,a)')
     >              'WARNING : CHANGE OF STATE PROB.'//
     >               ' IS LESS THAN MIN. IN ',
     >                       pwarn_count, ' CASES'
         call prc (comment)
         call prb
c
      endif
c
c      DO 960 IR = IRSEP-1, NRS
c        WRITE (6,9032) 1,1,1,1,1,NIZS,NIZS,NIZS,NIZS,NIZS
c        DO 950 IK = 1, NKS(IR)
c          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),
c     >      KNHS(IK,IR),KNBS(IK,IR),KTEBS(IK,IR),
c     >      KFIZS(IK,IR,0),KPCHS(IK,IR,0),kpmtc(ik,ir),
c     >     (KFIZS(IK,IR,IZ),KFRCS(IK,IR,IZ),KFCXS(IK,IR,IZ),
c     >      KPCHS(IK,IR,IZ),KPRCS(IK,IR,IZ),IZ=1,NIZS,MAX(1,NIZS-1))
c  950   CONTINUE
c  960 CONTINUE
c
c     Print out data concerning the targets and target conditions.
c
      if (cprint.eq.3.or.cprint.eq.9) then
         write (6,*) 'Target data summary:'

         do 970 id = 1,nds
            write(6,9060) irds(id),kteds(id),ktids(id),knds(id),ikds(id)
  970    continue
c
      endif
c
9060  format(5x,i4,3(1x,g12.5),3x,i4)


c
c     Debug test ...
c
      ik=45
      ir=15

      do iz = 1,6
         write(6,'(a,i5,10(1x,g12.5))') 'DEBUG RATES:',iz,
     >           kpchs(ik,ir,iz),kprcs(ik,ir,iz), 
     >           kfizs(ik,ir,iz),kfrcs(ik,ir,iz),
     >           kfcxs(ik,ir,iz)
      end do

c
c     Print out rate information data
c
C
c      DO 960 IR = IRSEP-1, IRSEP
c        WRITE (6,9032) 1,1,1,1,1,NIZS,NIZS,NIZS,NIZS,NIZS
c        DO 950 IK = 1, NKS(IR)
c          WRITE (6,9033) IK,IR,RS(IK,IR),ZS(IK,IR),
c     >      KNHS(IK,IR),KFIZS(IK,IR,0),KPCHS(IK,IR,0),
c     >     (KFIZS(IK,IR,IZ),KFRCS(IK,IR,IZ),KFCXS(IK,IR,IZ),
c     >      KPCHS(IK,IR,IZ),KPRCS(IK,IR,IZ),IZ=1,NIZS,MAX(1,NIZS-1))
c  950   CONTINUE
c  960 CONTINUE
c
C
C
C       GENERATE AN ANALYSIS OF THE IONIZATION DATA FOR CONSISTENCY
C       CHECK PURPOSES.
C
        IF ((CPINOPT.EQ.1.OR.CPINOPT.EQ.4).AND.PINCHECK) THEN
           CALL PINPRN
           CALL PINCHK
        ENDIF
c slmod begin - new

c      WRITE(0,*) 'DUMPING TRIANGLES AT THE END OF TAUIN1'
c      CALL DumpTriangles
c      WRITE(0,*) 'DUMPING QUADRANGLES AT THE END OF TAUIN1'
c      CALL DumpQuadrangles

c...Add print option:
      CALL OutputData(87,'END OF TAUIN1')
      CALL OutputEIRENE(65,'END OF TAUIN1')

      IF (outmode.GE.3) CALL DB('END OF TAUIN1')
      IF (sloutput) WRITE(0,*) 'END OF TAUIN1'
c slmod end
c
      RETURN
C
C
c 9001 FORMAT(/1X,'GEOMETRY DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
c     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
c     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
c     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
c     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
c     >  /5X,'TOROIDAL FIELD     BPHI =  ',F8.3,
c     >  /5X,'NO OF POINTS NP     =',I6,5X,'MAX NO. ROWS MKS    =',I6,
c     >  /5X,'SEPARATRIX   IRSEP  =',I6,5X,'WALL         IRWALL =',I6,
c     >  /5X,'FIRST TRAP   IRTRAP =',I6,5X,'NO OF RINGS  NRS    =',I6,
c     >  /5X,'K SPLIT PT   IKT    =',I6,5X,'K REF POINT  IKREF  =',I6)
 9002 FORMAT(/1X,' IK IR    R           Z          BPH',
     >  'I     TEB     TIB     NB      E1      VB          S      ',
     >  'BTOT/BTHETA    FEG1    FIG1',/1X,131('-'))
 9003 FORMAT(1X,2I3,2F12.8,f7.3,2f8.1,1P,E8.1,0P,2A9,G14.8,F8.2,3X,2A9)
c 9006 FORMAT(/1X,'TAUIN1: AREA OF SOL+TRAP =',F9.6,',  MAIN P =',F9.6,/
c 9011 FORMAT(/1X,'EDGE PLASMA DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
c     >  /5X,'RUN                 =',A,
c     >  /5X,'PLASMA       HMASS  =',I6,5X,'IMPURITY     ZMASS  =',I6,
c     >  /5X,'FIRST RING   IRCENT =',I6,5X,'LOST INN RGS NINOMP =',I6,
c     >  /5X,'LOST OUT RGS NINOWA =',I6)
c 9012 FORMAT(1P,6E12.4)
 9022 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
     >  ' IROUT  INDIST OUTDIST BACDIST FORDIST VOLUME  INPROB',
     >  '  POLDIST',/1X,131('-'))
 9023 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,2X,4F8.3,F8.5,F8.4,F8.3)
 9032 FORMAT(/1X,' IK IR    R      Z       NH0      NE      TE   ',
     >    ' TAUIZ0  CHPROB0',
     >    ' MTCPROB ',
     >  2('TAUIZ',I1,'  TAUEI',I1,'  TAURC',I1,'  CHPROB',I1,
     >  ' RCFRAC',I1,' '),/1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,2X,1P,6E8.1,2(1P,4E8.1,0P,F8.4))
c 9040 FORMAT(
c     >   5X,'NO OF R PTS       NXS   =  ',I5,',   DELTAR =',F9.5,
c     >  /5X,'NO OF Z PTS       NYS   =  ',I5,',   DELTAZ =',F9.5,
c     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
c     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
c     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
c     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
c     >  /5X,'INNER TARGET PTS  NDSIN =  ',I5,
c     >  /5X,'TOTAL TARGET PTS  NDS   =  ',I5)
c 9042 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
c     >  ' IROUT  ID   INTHETA OUTTHETA THETA   DIST',
c     >  /1X,131('-'))
c 9043 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,I7,2X,3F8.2,F8.5)
c 9044 FORMAT(I7,E14.0,A20)
 9050 FORMAT(I5,I5,50F8.3)
      END
C
C
C
      SUBROUTINE TAUIN2 (NIZS)
      
      !
      ! The Taus module calculates the characteristic times as 
      ! used in DIVIMP. It requires a call to an init routine to 
      ! set locally stored data - after that it is stand-alone. 
      !
      ! jde
      !
      use taus

      IMPLICIT  NONE
      INTEGER   NIZS
C
C***********************************************************************
C
C       SETS UP VALUES IN COMMON BLOCKS COMTAU/COMT2
C       THESE VALUES MAY DEPEND ON CTEMAV, AVERAGE IONISATION TEMP
C       RETURNED FROM A CALL TO NEUT.
C
C***********************************************************************
C
C     INCLUDE   "PARAMS"
      include 'params'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C     INCLUDE   "CIONIZ"
      include 'cioniz'
C     INCLUDE   "COMTOR"
      include 'comtor'
C
      INTEGER   IZ,IK,IR,I
      !REAL      LAMBDA,ROOTMI,ROOTTT
      !PARAMETER (LAMBDA=15.0)
      !REAL      FTAU,FTAUP,FTAUS,FTAUT,RIZSQR,STAU,TAU
      CHARACTER C(10)*9,FACTOR*9
C
      !ROOTMI = SQRT (CRMI)
      !ROOTTT = SQRT (CTEMAV)

      call init_taus(crmb,crmi,rizb,cioptb,cioptc,
     >               cioptd,czenh,cizeff,ctemav,irspec,qtim)

C
C-----------------------------------------------------------------------
C                     SET UP CFPS, CFSS, CFTS
C  NOTE 215: EXTRA HEATING OPTION, SET UP CONSTANT C215A.
C  NOTE 350: EXTRA STOPPING AND COLLISION OPTIONS, SET UP CONSTANTS.
C-----------------------------------------------------------------------
C
C
      DO  IZ = 1, MIN (CION,NIZS)
        DO IR = 1, NRS
          DO IK = 1, NKS(IR)

            call eval_taus(ik,ir,iz,knbs(ik,ir),ktibs(ik,ir),
     >                  kfps(ik,ir,iz),
     >                  kkkfps(ik,ir,iz),kfss(ik,ir,iz),kfts(ik,ir,iz))
          enddo
        enddo
      enddo

c
C
C---- PRINT SOME RESULTS - for debugging
C
      if (cprint.eq.5.or.cprint.eq.9) then
c
      DO 960 IR = IRSEP, IRWALL-1
        WRITE (6,9002) 1,1,1,NIZS,NIZS,NIZS
        DO 950 IK = 1, NKS(IR)
          C(1) = FACTOR (KFPS(IK,IR,1),2)
          C(2) = FACTOR (KFSS(IK,IR,1),3)
          C(3) = FACTOR (KFTS(IK,IR,1),4)
          C(4) = FACTOR (KFPS(IK,IR,NIZS),2)
          C(5) = FACTOR (KFSS(IK,IR,NIZS),3)
          C(6) = FACTOR (KFTS(IK,IR,NIZS),4)
          WRITE (6,9003) IK,IR,RS(IK,IR),ZS(IK,IR),(C(I),I=1,6)
  950   CONTINUE
  960 CONTINUE
      endif
c
c
      RETURN
C
 9002 FORMAT(/1X,' IK IR    R      Z      ',
     >  2('TAUP',I1,'    TAUS',I1,'    TAUT',I1,'    '),
     >  /1X,131('-'))
 9003 FORMAT(1X,2I3,2F7.3,2X,6A9)
      END
C
C
C
      SUBROUTINE TAUCHK (IK,IR,IZ,SPARA,TI)
      IMPLICIT NONE
      INTEGER IK,IR,IZ
      REAL SPARA,TI
C
C  *********************************************************************
C  *                                                                   *
C  *  TAUCHK: THIS ROUTINE CURRENT VALUES OF TAU PARALLEL, ETC         *
C  *  AGAINST THEIR "EXPECTED" VALUES ACCORDING TO THE ORIGINAL        *
C  *  FORMULAE.  IT IS CALLED PERIODICALLY FROM THE DIV ROUTINE.       *
C  *                                                                   *
C  *  SEE THE FOLLOWING NOTES FOR THE SPECIFICATIONS :-                *
C  *                                                                   *
C  *             GENERAL  OPTION0  OPTION1  OPTION2  OPTION3  OPTION4  *
C  *  COLLISION   8,50     3,121              103      284      350    *
C  *  FRICTION    50       3,121              103      350       -     *
C  *  HEATING     3        3,121     103             215,350     -     *
C  *                                                                   *
C  *                                      C.M.FARRELL   MAY 1989       *
C  *                                                                   *
C  *********************************************************************
C
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CLOCAL"
      include 'clocal'
C
      REAL K,MI,TB,MB,NB,ZB,ZI,LAM,ZENH,ZEFF,TAUP,TAUS,TAUT
      REAL DPARA,DS,TEMP
      CHARACTER C(3)*9,FACTOR*9
C
      K      = KKS(IR)
      MI     = CRMI
      TB     = KTIBS(IK,IR)
      MB     = CRMB
      NB     = KNBS(IK,IR) * 1.E-18
      ZB     = RIZB
      ZI     = REAL(IZ)
      LAM    = 15.0
      ZENH   = CZENH
      ZEFF   = REAL(CIZEFF)
C
      IF (CIOPTB.EQ.0.OR.CIOPTB.EQ.3.OR.IR.LT.IRSPEC.OR.
     >    (CIOPTB.EQ.4.AND.TI.LE.TB*MI/MB)) THEN
        TAUP = MI*SQRT(TB)*TI/(6.8E4*SQRT(MB)*NB*ZB*ZB*ZI*ZI*LAM*ZENH)
      ELSEIF (CIOPTB.EQ.1) THEN
        TAUP = HI
      ELSEIF (CIOPTB.EQ.2) THEN
        TAUP = SQRT(MI)*TI**1.5/(6.8E4*NB*ZB*ZEFF*ZI*ZI*LAM*ZENH)
      ELSEIF (CIOPTB.EQ.4) THEN
        TAUP = TI**2.5*MB/(9.E4*SQRT(MI)*TB*NB*ZB*ZB*ZI*ZI*LAM*ZENH)
      ENDIF
C
      DPARA  = 1.22E8*TI*TAUP/MI
      IF (CIOPTB.EQ.3.AND.IR.GE.IRSPEC) THEN
        DS   = SQRT(2.*DPARA*TAUP)
      ELSE
        DS   = SQRT(2.*DPARA*QTIM)
      ENDIF
C
      IF ((CIOPTC.EQ.0.or.cioptc.eq.4).OR.IR.LT.IRSPEC.OR.
     >    (CIOPTC.EQ.3.AND.TI.LE.TB*MI/MB)) THEN
        TAUS = MI*TB**1.5/(6.8E4*SQRT(MB)*(1.+MB/MI)*NB*ZB*ZB*ZI*ZI*LAM*
     >         ZENH)
      ELSEIF (CIOPTC.EQ.1) THEN
        TAUS = HI
      ELSEIF (CIOPTC.EQ.2) THEN
        TAUS = TAUP
      ELSEIF (CIOPTC.EQ.3) THEN
        TAUS =TI**1.5*SQRT(MI)/(9.E4*(1.+MI/MB)*NB*ZB*ZB*ZI*ZI*LAM*ZENH)
      ENDIF
C
      IF     (CIOPTD.EQ.0.OR.(CIOPTD.NE.3.AND.IR.LT.IRSPEC)) THEN
        TAUT = MI*TB**1.5/(1.4E5*SQRT(MB)*NB*ZB*ZB*ZI*ZI*LAM*ZENH)
      ELSEIF (CIOPTD.EQ.1) THEN
        TAUT = HI
      ELSEIF (CIOPTD.EQ.2) THEN
        TAUT = 0.0
      ELSEIF (CIOPTD.EQ.3) THEN
        TAUT = (MI*TB+MB*TI)**1.5/(1.4E5*SQRT(MB*MI)*ZB*ZB*NB*ZI*ZI*
     >         LAM*ZENH)
      ENDIF
C
      TEMP  = CTEMAV
      CTEMAV = TI
      C(1)  = FACTOR (LFPS(IK,IR,IZ),2)
      CTEMAV = TEMP
      C(2)  = FACTOR (LFSS(IK,IR,IZ),3)
      C(3)  = FACTOR (LFTS(IK,IR,IZ),4)
      WRITE (6,9001) IK,IR,K,IZ,TB,NB,TI,LTOLDS(IK,IR,IZ),TAUP,TAUS,
     >  TAUT,DS,C(1),C(2),C(3),SPARA
      RETURN
C
 9001 FORMAT(/1X,'TAUCHK: IK',I3,' IR',I3,' K',F6.3,' ZI',I2,' TB',F9.2,
     >  ' NB',F7.3,' TI',F9.2,' TOLD',F9.2,/4X,1P,'EXPECTED VALS:- ',
     >  'TAUP',E9.2,' TAUS',E9.2,' TAUT',E9.2,' DELTAS',F11.5,
     >  ' (OR 0.0)',/4X,'ACTUAL VALS:-   TAUP',
     >  A9,' TAUS',A9,' TAUT',A9,' DELTAS',F11.5)
      END
C
C
C
C@PROCESS OPT(0),VECTOR(LEV(0))
C
C
      SUBROUTINE TAUIN3 (NIZS,ABSFAC)
      IMPLICIT  NONE
      INTEGER   NIZS
      REAL ABSFAC
      WRITE (6,*) ' TAUIN3: NOT IMPLEMENTED...'
      STOP
      END
C
C
      SUBROUTINE TAUVOL
      IMPLICIT NONE
C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "COMTOR"
      include 'comtor'
      include 'cedge2d'
c slmod begin
      INCLUDE 'slcom'
c slmod end
C
C  *********************************************************************
C  *                                                                   *
C  *  TAUVOL: THIS FUNCTION CALCUALTES THE AREA IN THE POLOIDAL PLANE  *
C  *  ASSOCIATED WITH EACH (IK,IR) GRID POINT, AND THE ELEMENTAL       *
C  *  VOLUMES.                                                         *
C  *                                                                   *
C  *  CHRIS FARRELL  (HUNTERSKIL)  JUNE 1989                           *
C  *                                                                   *
C  *  MODIFIED TO CALCULATE THE TRUE PLASMA CELL AREA USING THE        *
C  *  POSITIONS OF THE CELL VERTICES AS WELL AS THE DIVIMP             *
C  *  APPROXIMATION.  FOR COMPARISION THE AREA OBTAINED FROM THE XY    *
C  *  GRID NEAREST NEIGHBOURS IS ALSO CALCULATED.  THE VOLUME OF EACH  *
C  *  CELL IS NOW LOADED INTO KVOLS AND KVOL2, DEPENDING ON WHICH      *
C  *  AREA CALCULATION WAS USED.                                       *
C  *                                                                   *
C  *  LORNE HORTON (JET) MARCH 1993                                    *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IK,IR,JK(5,5),JR(5,5),I,J,IX,IY,K,KP,L,LP1,cnt
      REAL    S1,S2,S3,S4,RRS(5,5),ZZS(5,5),S,XX
      REAL    ACHK,RNGCHK,TOTCHK
      REAL    TOTA,TOTV,TOTA2,TOTV2,tmpval
c slmod begin
      REAL*8  :: AREA_SUM
      LOGICAL :: WARNING_MESSAGE = .TRUE.
c slmod end
      real    apara,apol,afact
C
      TOTA   = 0.0
      TOTV   = 0.0
      TOTA2  = 0.0
      TOTV2  = 0.0
      TOTCHK = 0.0
C
      DO 910 IR = 1, NRS
c slmod begin
      IF (IDRING(IR).EQ.BOUNDARY) CYCLE
c slmod end
      KTOTAS(IR) = 0.0
      KTOTVS(IR) = 0.0
      KTOTA2(IR) = 0.0
      KTOTV2(IR) = 0.0
      RNGCHK     = 0.0
      WRITE (9,9005)
      DO 900 IK = 1, NKS(IR)
C
C---- SET UP 5 BY 5 MATRIX OF (IK,IR) VALUES SURROUNDING POINT, USING
C---- POSITIONS 1,3,5 ONLY IN THE MATRIX.  THE INTERMEDIATE LOCATIONS
C---- WILL BE USED TO REPRESENT THE MIDPOINTS.
C
c slmod begin
c...  Sorry for the spagetti, but the (legacy) code between here and
c     205 is completely broken for linear device grids.  I have set
c     KAREAS to be assigned from KAREAS2 below:
      IF (cgridopt.EQ.LINEAR_GRID) GOTO 205
c slmod end
      JK(1,3) = IK - 1
      IF (IK.EQ.1.AND.IR.LT.IRSEP) JK(1,3) = NKS(IR) - 1
      IF (IK.EQ.1.AND.IR.GE.IRSEP) JK(1,3) = 1
      JK(3,3) = IK
      JK(5,3) = IK + 1
      IF (IK.EQ.NKS(IR).AND.IR.LT.IRSEP) JK(5,3) = 2
      IF (IK.EQ.NKS(IR).AND.IR.GE.IRSEP) JK(5,3) = NKS(IR)
      JR(1,3) = IR
      JR(3,3) = IR
      JR(5,3) = IR
      JK(1,1) = IKINS(JK(1,3),JR(1,3))
      JK(3,1) = IKINS(JK(3,3),JR(3,3))
      JK(5,1) = IKINS(JK(5,3),JR(5,3))
      JR(1,1) = IRINS(JK(1,3),JR(1,3))
      JR(3,1) = IRINS(JK(3,3),JR(3,3))
      JR(5,1) = IRINS(JK(5,3),JR(5,3))
      JK(1,5) = IKOUTS(JK(1,3),JR(1,3))
      JK(3,5) = IKOUTS(JK(3,3),JR(3,3))
      JK(5,5) = IKOUTS(JK(5,3),JR(5,3))
      JR(1,5) = IROUTS(JK(1,3),JR(1,3))
      JR(3,5) = IROUTS(JK(3,3),JR(3,3))
      JR(5,5) = IROUTS(JK(5,3),JR(5,3))
C
C---- SET (R,Z) VALUES AT THE 9 NEIGHBOURS AND AT THEIR MIDPOINTS.
C---- NOTE THAT THE FOUR "MID-MIDPOINTS" AT (2,2), (2,4), (4,2), (4,4)
C---- CAN ONLY BE CALCULATED ONCE THE OTHER MIDPOINTS ARE KNOWN, AND
C---- THAT THEY COULD BE CALCULATED IN TWO WAYS - USING I-1,I+1 OR
C---- J-1,J+1 - THE LATTER IS CHOSEN HERE.
C
      DO 100 I = 1, 5
       DO 100 J = 1, 5
        IF ((I.EQ.2.OR.I.EQ.4).AND.(J.EQ.2.OR.J.EQ.4)) THEN
         GOTO 100
        ELSEIF (I.EQ.2.OR.I.EQ.4) THEN
         RRS(I,J) = .5*(RS(JK(I-1,J),JR(I-1,J))+RS(JK(I+1,J),JR(I+1,J)))
         ZZS(I,J) = .5*(ZS(JK(I-1,J),JR(I-1,J))+ZS(JK(I+1,J),JR(I+1,J)))
        ELSEIF (J.EQ.2.OR.J.EQ.4) THEN
         RRS(I,J) = .5*(RS(JK(I,J-1),JR(I,J-1))+RS(JK(I,J+1),JR(I,J+1)))
         ZZS(I,J) = .5*(ZS(JK(I,J-1),JR(I,J-1))+ZS(JK(I,J+1),JR(I,J+1)))
        ELSE
         RRS(I,J) = RS(JK(I,J),JR(I,J))
         ZZS(I,J) = ZS(JK(I,J),JR(I,J))
        ENDIF
  100 CONTINUE
C
c
c     The values need to be adjusted for conditions
c     where the target does not coincide with the last
c     grid positions. In these cases the 2,4 array values need
c     to be set with the values of the target coordinates
c     for the ring in question. This would be either
c     row 4 for the ik=nks end of a ring greater than irsep
c     or row 2 for the ik=1 end of the ring.
c
      if (ik.eq.1.and.ir.ge.irsep.and.
     >    (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >     ctargopt.eq.3.or.ctargopt.eq.5.or.
     >     ctargopt.eq.6)
     >    ) then
         rrs(2,1) = rp(idds(irins(ik,ir),2))
         zzs(2,1) = zp(idds(irins(ik,ir),2))
         rrs(2,3) = rp(idds(ir,2))
         zzs(2,3) = zp(idds(ir,2))
         rrs(2,5) = rp(idds(irouts(ik,ir),2))
         zzs(2,5) = zp(idds(irouts(ik,ir),2))
      endif
c
      if (ik.eq.nks(ir).and.ir.ge.irsep.and.
     >    (ctargopt.eq.1.or.ctargopt.eq.2.or.
     >     ctargopt.eq.3.or.ctargopt.eq.5.or.
     >     ctargopt.eq.6)
     >    ) then
         rrs(4,1) = rp(idds(irins(ik,ir),1))
         zzs(4,1) = zp(idds(irins(ik,ir),1))
         rrs(4,3) = rp(idds(ir,1))
         zzs(4,3) = zp(idds(ir,1))
         rrs(4,5) = rp(idds(irouts(ik,ir),1))
         zzs(4,5) = zp(idds(irouts(ik,ir),1))
      endif
c
c
      DO 110 I = 2, 4, 2
       DO 110 J = 2, 4, 2
        RRS(I,J) = 0.5 * (RRS(I,J-1) + RRS(I,J+1))
        ZZS(I,J) = 0.5 * (ZZS(I,J-1) + ZZS(I,J+1))
  110 CONTINUE
C
C---- CALCULATE THE FOUR AREAS (SOME MAY BE 0) AND SUM THEM.
C---- AREA OF 4-SIDED IRREGULAR FIGURE IS THAT OF TWO TRIANGLES.
C
      KAREAS(IK,IR) = 0.0
      DO 200 I = 2, 3
       DO 200 J = 2, 3
        S1 = SQRT ((RRS(I,J+1)-RRS(I,J))**2 + (ZZS(I,J+1)-ZZS(I,J))**2)
        S2 = SQRT ((RRS(I+1,J)-RRS(I,J))**2 + (ZZS(I+1,J)-ZZS(I,J))**2)
        S3 = SQRT ((RRS(I,J+1)-RRS(I+1,J+1))**2 +
     >             (ZZS(I,J+1)-ZZS(I+1,J+1))**2)
        S4 = SQRT ((RRS(I+1,J)-RRS(I+1,J+1))**2 +
     >             (ZZS(I+1,J)-ZZS(I+1,J+1))**2)
        XX = SQRT ((RRS(I,J+1)-RRS(I+1,J))**2 +
     >             (ZZS(I,J+1)-ZZS(I+1,J))**2)
        IF (S1.LE.0.0.OR.S2.LE.0.0.OR.S3.LE.0.0.OR.S4.LE.0.0) GOTO 200
        S    = 0.5 * (S1 + S2 + XX)
        IF (S*(S-S1)*(S-S2)*(S-XX).GE.0.0) THEN
           KAREAS(IK,IR) = KAREAS(IK,IR)
     +                  + SQRT (S * (S-S1) * (S-S2) * (S-XX))
        ELSE
           WRITE(6,*) ' S EXP <0 (AREA):1,2', S,S1,S2,XX,I,J,
     +                KAREAS(IK,IR)
        ENDIF
        S    = 0.5 * (S3 + S4 + XX)
        IF (S*(S-S3)*(S-S4)*(S-XX).GE.0.0) THEN
           KAREAS(IK,IR) = KAREAS(IK,IR)
     +                  + SQRT (S * (S-S3) * (S-S4) * (S-XX))
        ELSE
           WRITE(6,*) ' S EXP <0 (AREA):3,4', S,S3,S4,XX,I,J,
     +                KAREAS(IK,IR)
        ENDIF
  200 CONTINUE
C
C---- CALCULATE AREA OF CELL AS DEDUCED FROM XY GRID
C
c     This has been commented out because for a large XY
c     grid - the computational time is ridiculous!
c
c     The iteration count is sigma(ir) 1500*1500*nks(ir)
c     IR is 77 for ITER grids and nks(ir) has a mean of
c     40 or so ... yielding a loop count of 6.9e9
c     iterations! For a non-essential chunk of code!
c
c
       ACHK = 0.0
c      DO IX = 1, NXS
c       DO IY = 1, NYS
c        IF (IRXYS(IX,IY).EQ.IR.AND.IKXYS(IX,IY).EQ.IK
c     +      .AND.IFXYS(IX,IY).EQ.1) ACHK = ACHK + DR*DZ
c       ENDDO
c      ENDDO
C
C---- CALCULATE TRUE CELL AREA FROM POSITION OF POLYGON VERTICES.
C---- THIS AREA WILL ONLY BE AVAILABLE FOR PLASMA POINTS.  AREAS
C---- OF VIRTUAL POINTS ARE SET TO ZERO.  NOTE THAT THE ALGORITHM
C---- FOR CALCULATING THE AREA OF A POLYGON FROM THE POSITION OF
C---- ITS VERTICES ASSUMES THAT THE VERTICES ARE ORDERED CLOCKWISE
C---- AROUND THE PERIMETER.
C
c slmod begin
 205  CONTINUE
      AREA_SUM = 0.0D0
c slmod end
      KP = KORPG(IK,IR)
      KAREA2(IK,IR) = 0.0
      IF (KP.GT.0) THEN
        DO L = 1, NVERTP(KP)
          LP1 = L + 1
          IF (L.EQ.NVERTP(KP)) LP1 = 1
          KAREA2(IK,IR) = KAREA2(IK,IR) + (RVERTP(LP1,KP)*ZVERTP(L,KP)
     +                                    - RVERTP(L,KP)*ZVERTP(LP1,KP))
        ENDDO
c
c       Ensure that the area is greater than zero.
c
        KAREA2(IK,IR) = 0.5 * abs(KAREA2(IK,IR))
c slmod begin
        IF     (KAREA2(IK,IR).LT.-1.0E-10) THEN 
          CALL ER('TAUVOL','Malformed cell detected, stopping',*99)       
        ELSEIF (KAREA2(IK,IR).LT.1.0E-06) THEN 
          IF (WARNING_MESSAGE) THEN
            WRITE(0,*) 'WARNING: Small area cell detected, using '//
     +                 'double precision calculation in TAUVOL'
            WARNING_MESSAGE = .FALSE.
          ENDIF
          KAREA2(IK,IR) = 0.0
          DO L = 1, NVERTP(KP)
            LP1 = L + 1
            IF (L.EQ.NVERTP(KP)) LP1 = 1
            AREA_SUM = AREA_SUM + 
     +                 DBLE(RVERTP(LP1,KP)) * DBLE(ZVERTP(L  ,KP)) -
     +                 DBLE(RVERTP(L  ,KP)) * DBLE(ZVERTP(LP1,KP))
c            AREA_SUM = AREA_SUM + DBLE(RVERTP(LP1,KP) * ZVERTP(L  ,KP))-
c     +                            DBLE(RVERTP(L  ,KP) * ZVERTP(LP1,KP))
            IF (.TRUE.) THEN
              KAREA2(IK,IR)=KAREA2(IK,IR)+(RVERTP(LP1,KP)*ZVERTP(L,KP)
     +                                   - RVERTP(L,KP)*ZVERTP(LP1,KP))
              WRITE(6,*) 'WARNING: SMALL AREA CELL DETECTED'
              WRITE(6,*) '  KAREA2:',ik,ir,l
              WRITE(6,*) '        :',karea2(ik,ir)
              WRITE(6,*) '        :',area_sum
              WRITE(6,*) '        :',lp1,kp,rvertp(lp1,kp)
              WRITE(6,*) '        :',l  ,kp,rvertp(l  ,kp)
              WRITE(6,*) '        :',RVERTP(LP1,KP)*ZVERTP(L,KP)
              WRITE(6,*) '        :',l  ,kp,rvertp(l  ,kp)
              WRITE(6,*) '        :',lp1,kp,zvertp(lp1,kp)
              WRITE(6,*) '        :',RVERTP(L,KP)*ZVERTP(LP1,KP)
            ENDIF
          ENDDO
          KAREA2(IK,IR) = SNGL(0.5D0*AREA_SUM)
        ENDIF
c slmod end
      ENDIF
c slmod begin
      IF (cgridopt.EQ.LINEAR_GRID) THEN
        kareas = karea2
      ENDIF
c slmod end
C
C---- STORE VOLUMES AS WELL AS AREAS
C
c     For ASDEX grids the volumes are loaded in the RASDEX
c     subroutine.
c
c     For ASDEX UPGRADE calculate the values here as for JET.
c
      if (cgridopt.eq.0.or.cgridopt.eq.3.or.
     .    cgridopt.eq.LINEAR_GRID) then
         KVOLS(IK,IR) = 2.0*PI*RS(IK,IR)*KAREAS(IK,IR)
      endif
      KVOL2(IK,IR) = 2.0*PI*RS(IK,IR)*KAREA2(IK,IR)

c slmod begin
c...  Scale KVOLS and KVOL2 if only a fraction of the torus is 
c     being modelled in EIRENE:
      IF (eirtorfrac.NE.1.0) THEN
        IF (ik.EQ.1.AND.ir.EQ.1) THEN
          WRITE(0     ,*) 'SCALING TOROIDAL VOLUMES'
          WRITE(PINOUT,*) 'SCALING TOROIDAL VOLUMES',eirtorfrac
        ENDIF
        kvols(ik,ir) = kvols(ik,ir) * eirtorfrac
        kvol2(ik,ir) = kvol2(ik,ir) * eirtorfrac
      ENDIF
c slmod end
c
c     After having calculated KVOLS and KAREAS using both methods -
c     store the one that will be used for calculations in DIVIMP
c     as selected by the CVOLOPT input parameter, into the
c     KAREAS and KVOLS arrays - the other will be saved in the
c     second array. Note - since polygons do not exist for all the
c     recorded grid points - the polygon method will not be able
c     to calculate areas for all cells - in these cases just use zeroes
c     since if one is using this method it should not be used in
c     conjunction with other methods that follow particles into
c     these cells. The KAREA2 array always contains areas based on
c     polygons because it is used explicitly in conjunction with many
c     PIN quantities - as such - it is NOT copied over.
c
c     Note: if CVOLOPT = 0 nothing needs to be done to the arrays
c           because this is the default at this time. If CVOLOPT=1
c           has been specified one needs to change the arrays.
c
      if (cvolopt.eq.1) then
         kareas(ik,ir) = karea2(ik,ir)
         kvols(ik,ir) = kvol2(ik,ir)
      endif
C
      WRITE (9,9006) IK,IR,RS(IK,IR),ZS(IK,IR),KAREAS(IK,IR),
     >  KVOLS(IK,IR),KAREA2(IK,IR),KVOL2(IK,IR),ACHK,CVMF(IK,IR),
     >  e2dareas(ik,ir)
C
C---- SUM OVER RING
C
      IF (IK.NE.1.OR.IR.GE.IRSEP) THEN
        KTOTAS(IR) = KTOTAS(IR) + KAREAS(IK,IR)
        KTOTVS(IR) = KTOTVS(IR) + KVOLS(IK,IR)
        KTOTA2(IR) = KTOTA2(IR) + KAREA2(IK,IR)
        KTOTV2(IR) = KTOTV2(IR) + KVOL2(IK,IR)
        RNGCHK     = RNGCHK     + ACHK
      ENDIF
  900 CONTINUE
      WRITE (9,9007) IR,KTOTAS(IR),KTOTVS(IR),KTOTA2(IR),KTOTV2(IR),
     +               RNGCHK
C
C---- TOTALS
C
      TOTA   = TOTA   + KTOTAS(IR)
      TOTV   = TOTV   + KTOTVS(IR)
      TOTA2  = TOTA2  + KTOTA2(IR)
      TOTV2  = TOTV2  + KTOTV2(IR)
      TOTCHK = TOTCHK + RNGCHK
  910 CONTINUE
      WRITE (9,9008) TOTA,TOTV,TOTA2,TOTV2,TOTCHK
c
c     Calculate some characteristics of Areas
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
      write(6,*) 'Area factors:'
c
      do ir = 1,nrs
        do ik = 1,nks(ir)
           if (ik.eq.1) then
             apol = karea2(ik,ir) / ksb(ik,ir)
           else
             apol = karea2(ik,ir) / (ksb(ik,ir)-ksb(ik-1,ir))
           endif
c
           apara = apol / kbfs(ik,ir)
c
           afact = apara / rs(ik,ir)
c
           write (6,'(''Areas:'',2i4,5(2x,e15.6))') ir,ik,
     >                karea2(ik,ir),kbfs(ik,ir),
     >                apol,apara,afact
c
        end do
      end do
c
      endif
c
c     Loop through and check for zero volume cells.
c     Print an informational message.
c
      cnt = 0
      do ir = 1,nrs
         do ik = 1,nks(ir)
            if (kareas(ik,ir).le.0.0) then
c
               if (cnt.eq.0) then
                  write (6,*)
                  write (6,*) 'Listing of ZERO Volume cells:'
                  write (6,*)
                  cnt = cnt + 1
               endif
c
               if (ir.eq.1.or.ir.eq.irwall.or.ir.eq.irtrap) then
                  write (6,'(a10,2(1x,i4),2(1x,g12.5))')
     >                  ' ',ik,ir,kareas(ik,ir),karea2(ik,ir)
               else
                  write (6,'(a10,2(1x,i4),2(1x,g12.5))')
     >               'WARNING!:',ik,ir,kareas(ik,ir),karea2(ik,ir)
                  if (cnt.eq.1) then
                     call prb
                     call prc('WARNING!: GRID CONTAINS ZERO'//
     >                  ' VOLUME CELLS IN NON-BOUNDARY RINGS')
                     call prb
                     cnt = cnt + 1
                  endif
               endif
            endif
         end do
      end do
c
c     End of TAUVOL
c
C
 9005 FORMAT(/1X,' IK IR    R      Z    KAREAS     KVOLS    KAREA2  ',
     >  '   KVOL2    ACHECK     VMF     E2D-AREA',/1X,131('-'))
 9006 FORMAT(1X,2I3,2F7.3,1P,7E10.3)
 9007 FORMAT(/1X,'RING: ',I3,11X,1P,5E10.3)
 9008 FORMAT(//1X,'TOTAL: ',13X,1P,5E10.3)
C
      RETURN
c slmod begin
 99   WRITE(0,*) '  IK,IR=',ik,ir
      STOP
c slmod end
      END
c
c
c
      SUBROUTINE SKORXYDIV(MP,NP,ITAG,
     *                  MR,NR,MX,KORX,NI,
     *                  MC,NC,MY,KORY,NJ)
c
      implicit none
c
c      IMPLICIT REAL*8(A-H,O-Z)
C
c
C     THIS ROUTINE IS IDENTICAL TO THE NIMBUS ROUTINE
C     SKORXY ... WHICH IS USED BY NIMBUS TO DECODE THE
C     ITAG ARRAY AND SET UP THE CROSS-REFERENCING
C     INDICES - SO THAT THE DATA IN VARIOUS ARRAYS CAN
C     BE PROPERLY ASSIGNED TO THE VARIOUS GRID COORDINATES.
C     SINCE DIVIMP READS THE EQUILIBRIUM AND EDGE2D FILES FOR
C     SOME OF THIS DATA, IT ALSO NEEDS TO BE ABLE TO ASSIGN
C     THE ARRAY VALUES TO THE PROPER BINS. THIS ROUTINE IS
C     REQUIRED WHETHER OR NOT PIN IS INVOKED. NOTE: AT SOME
C     TIME THE NIMBUS ROUTINE MAY CHANGE - IN WHICH CASE
C     IT MAY BE NECESSARY TO CHECK THIS ROUTINE TO MAKE SURE
C     THAT IT STILL BEHAVES IN A CONSISTENT FASHION.
C
C     LDH AND JDE, FEB 9/93
C
      integer  mp,np,mr,nr,mx,my,nc,mc
      integer  ITAG(MP,5),KORX(MR,MX),NI(MR),KORY(MC,MY),NJ(MC)
C
C I            NP = NUMBER OF POINTS
C I     ITAG(K,1) = INDEX OF FLUX LINE  ( <= NC )
C              2  =   "    " NORMAL TO FLUX LINES  ( <= NR )
C                   (=0 IF NO SWEEP IN RHO-DIRECTION)
C                   (-INDEX OF NORMAL  IF PERIODIC POINT AT THE TOP)
C              3  = TAG FOR X SWEEP
C              4  =  "   "  Y   "
C              5  =  "   "  M.C.
C I            MR = MAXIMUM NUMBER OF ROWS (= NORMAL TO FLUX LINES)
C O            NR = NUMBER OF ROWS
C I            MX = MAXIMUM NUMBER OF POINTS ALONG RHO
C O     KORX(J,I) = LOCATION IN VECTORS SUCH AS ITAG OF THE I-TH
C                   ELEMENT FOR SWEEPING ALONG X OF THE J-TH ROW
C O         NI(J) = NUMBER OF ELEMENTS IN ROW J
C I            MC = MAXIMUM NUMBER OF COLUMNS (= FLUX LINES)
C O            NC = NUMBER OF COLUMNS
C I            MY = MAXIMUM NUMBER OF POINTS ALONG THETA
C O     KORY(I,J) = LOCATION IN VECTORS SUCH AS ITAG OF THE J-TH
C                   ELEMENT FOR SWEEPING ALONG Y OF THE I-TH COLUMN
C O         NJ(I) = NUMBER OF ELEMENTS IN ROW I
C
C     COMPUTE NR,NC AND CHECK THAT THE WORM IS GIVEN BY COLUMNS
C
c     local variables
c
      integer i,j,k,i0,j0,i1,j1
c
      NC=0
      NR=0
      I0=ITAG(1,1)
      J0=ITAG(1,2)
      DO 240 K=2,NP
      I=ITAG(K,1)
      NC=MAX0(NC,I)
      J=ITAG(K,2)
      NR=MAX0(NR,IABS(J))
      IF(I.GT.MC .OR. J.GT.MR) GOTO 250
      IF(I-I0) 250,210,220
  210 IF(J.LE.0) GOTO 230
      IF(J.LE.J0) GOTO 250
  220 J0=J
  230 I0=I
  240 CONTINUE
      GOTO 280
  250 WRITE(6,260) K
  260 FORMAT(' ERRORE A K=',I3)
c      CALL EXIT(6)
      CALL EXIT
  280 CONTINUE
C
C               ACCESS VECTORS FOR SWEEPING ALONG X
      DO 330 J=1,NR
      I=0
      DO 320 K=1,NP
      J1=ITAG(K,2)
      IF(J1.NE.J) GOTO 320
      I=I+1
      IF(I.GT.MX) THEN
       WRITE(6,315) J
  315  FORMAT(' TOO MANY POINT AT ROW',I3)
c       CALL EXIT(6)
       CALL EXIT
       ENDIF
      KORX(J,I)=K
  320 CONTINUE
      NI(J)=I
  330 CONTINUE
C
C               ACCESS VECTORS FOR SWEEPING ALONG Y
      DO 380 I=1,NC
      J=0
      DO 350 K=1,NP
      I1=ITAG(K,1)
      IF(I1.NE.I) GOTO 350
      J=J+1
      IF(J.GT.MY) THEN
       WRITE(6,345) J
  345  FORMAT(' TOO MANY POINT AT COLUMN',I3)
c       CALL EXIT(6)
       CALL EXIT
       ENDIF
      KORY(I,J)=K
  350 CONTINUE
      NJ(I)=J
  380 CONTINUE
C
      RETURN
      END
c
c
c
      subroutine rjet
      implicit none
c
c     The purpose of this routine is to read the JET grids into the
c     standard DIVIMP ararys ... hopefully the majority of code to
c     interface to JET grids can be placed here and the interface
c     into DIVIMP standardized to make it easier to add new grid
c     configuartion to DIVIMP. Some of the data in these routines
c     will need to be customized and spread throughout the code
c     because DIVIMP was not written originally to handle double
c     null configurations. (This is for the ITER case.)
c
c     David Elder,    June 18 , 1993
c

C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CNEUT"
c      include 'cneut'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CIONIZ"
      include 'cioniz'
C     INCLUDE "READER"
      include 'reader'
C     INCLUDE "DYNAM5"
      include 'dynam5'
c
      include 'pindata'
      include 'baffles'
c
c     Edge2d values
c
      include 'cedge2d'
c slmod begin
      include 'slcom'
c slmod end
c
      CHARACTER MESAGE*72,C(10)*9,FACTOR*9
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IDUMMY(2),NJ(MAXNRS),NI(MAXNKS)
      INTEGER ITAG(MAXNKS*MAXNRS,5),ITAGDV(MAXNKS*MAXNRS)
      INTEGER KORX(MAXNKS,MAXNRS),KORYE(MAXNRS,MAXNKS)
      INTEGER ISEP,NINOMP,NINOWA,IREF,NZ,NZMAX,K0,IK0,IKMID
      INTEGER IND,IW,in
      REAL    DUMMY(MAXNKS*MAXNRS,11),HMASS,ZMASS,S,DELTAL,P,FACT
      REAL    sPERP,BEST,DSQ,R,Z,THIN,THOUT,DELTAR,DELTAZ
      REAL    MU,SMAX,SMID,SDIFF
      real    qovrt
      real    rtmp,ztmp
c
c     Hybrid wall data
c
      real rvmod(mves),zvmod(mves)
C
      integer    jplft,jprgt
c
      integer    ighost,ifluid,ieshot
      real       xji,xj,xjf,teslic
C
c      REAL    RHO(MAXNKS,MAXNRS)
C     REAL    THETA(MAXNKS,MAXNRS)
c
c     Moved to common block
c
c      REAL    HRO(MAXNKS,MAXNRS), HTETA(MAXNKS,MAXNRS)
C
C
c
C-----------------------------------------------------------------------
C     Extract geometry details from GRID2D output file
C-----------------------------------------------------------------------
C
      REWIND (4)
C
      MESAGE = 'Extracting GRID2D geometry details'
  100 READ (4,'(A72)',ERR=9999,END=110) BUFFER
C     WRITE (9,'(1X,A72,1X,A6)') BUFFER,'TAUIN1'
      IF (BUFFER(1:6).EQ.'    NP')
     +  READ (4,*) NP,IDUMMY(1),ISHOT,TSLICE,R0,Z0,CBPHI
      IF (NP.GT.MAXNKS*MAXNRS) THEN
        WRITE (6,*) ' NUMBER OF POINTS IN EQUILIBRIUM FILE EXCEEDS'//
     +              ' ALLOCATION'
        WRITE (6,*) ' NP = ',NP,' MAXNKS*MAXNRS = ',MAXNKS*MAXNRS
        STOP ' '
      ENDIF
      IF (BUFFER(1:6).EQ.'NYREFG')
     +  READ (4,*) IDUMMY(1),REFCT
      IF (BUFFER(1:6).EQ.' JPLFT')
     +  READ (4,*) JPLFT,JPRGT,RXP,ZXP,IDUMMY(1),IDUMMY(2),
     +              IRSEP,IRWALL
      IF (BUFFER(1:6).EQ.'ITAG  ')
     +  READ (4,'(12I6)') ((ITAG(I,J),J=1,5),I=1,NP)
      IF (BUFFER(1:6).EQ.'RMESH ') READ (4,*) (DUMMY(K,1),K=1,NP)
      IF (BUFFER(1:6).EQ.'ZMESH ') READ (4,*) (DUMMY(K,2),K=1,NP)
      IF (BUFFER(1:6).EQ.'RHO   ') READ (4,*) (DUMMY(K,3),K=1,NP)
      IF (BUFFER(1:6).EQ.'THETA ') READ (4,*) (DUMMY(K,4),K=1,NP)
      IF (BUFFER(1:6).EQ.'HRO   ') READ (4,*) (DUMMY(K,5),K=1,NP)
      IF (BUFFER(1:6).EQ.'HTETA ') READ (4,*) (DUMMY(K,6),K=1,NP)
c
c      IF (BUFFER(1:6).EQ.'BFI   ') READ (4,*) (DUMMY(K,7),K=1,NP)
c
      IF (BUFFER(1:6).EQ.'SH    ') READ (4,*) (DUMMY(K,8),K=1,NP)
      IF (BUFFER(1:6).EQ.'PSI   ') READ (4,*) (DUMMY(K,11),K=1,NP)
      IF (BUFFER(1:6).eq.'NVES  ') then
         read (4,*) nves
c
c        Check value of NVES
c
         if (nves.gt.maxpts) then
            write (6,*) 'ERROR: NUMBER OF VESSEL WALL ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE'
            write (6,*) 'NVES = ',nves,'  >  MAXPTS = ',maxpts
            call prc ('ERROR: NUMBER OF VESSEL WALL ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE')
            write (7,'(a7,i6,a13,i6)') 'NVES = ',nves,
     >                              '  >  MAXPTS = ',maxpts
            stop
         endif
      endif
c
      if (buffer(1:6).eq.'RVES  ') read (4,*) (rves(k),k=1,nves)
      if (buffer(1:6).eq.'ZVES  ') read (4,*) (zves(k),k=1,nves)
c
c     READ in information on baffles if any are present
c
c     Primary baffle -
c
      if (buffer(1:6).eq.'NBUFLE') then
         read (4,*) nbufle
c
c        Check value of NBUFLE
c
         if (nbufle.gt.mbufle) then
            write (6,*) 'ERROR: NUMBER OF PRIMARY BAFFLE'//
     >               ' ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE'
            write (6,*) 'NBUFLE = ',nbufle,'  >  MBUFLE = ',mbufle
            call prc('ERROR: NUMBER OF PRIMARY BAFFLE'//
     >               ' ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE')
            write (7,'(a7,i6,a13,i6)') 'NBUFLE = ',nbufle,
     >                              '  >  MBUFLE = ',mbufle
            stop
         endif
      endif
c
      if (buffer(1:6).eq.'RBUFLE') read (4,*) (rbufle(k),k=1,nbufle)
      if (buffer(1:6).eq.'ZBUFLE') read (4,*) (zbufle(k),k=1,nbufle)
c
c     Other defined baffles
c
      if (buffer(1:6).eq.'NBUFMX') then
c
         read (4,*) nbufmx
c
c        Check value of NBUFMX
c
         if (nbufmx.gt.mbufx) then
            write (6,*) 'ERROR: NUMBER OF SECONDARY BAFFLES'//
     >               ' EXCEEDS ARRAY SIZE'
            write (6,*) 'NBUFMX = ',nbufx,'  >  MBUFX = ',mbufx
            call prc('ERROR: NUMBER OF SECONDARY BAFFLES'//
     >               ' EXCEEDS ARRAY SIZE')
            write (7,'(a7,i6,a13,i6)') 'NBUFMX = ',nbufx,
     >                              '  >  MBUFX = ',mbufx
            stop
         endif
      endif
c
      if (buffer(1:6).eq.'NBUFX1') then
c
         read (4,*) (nbufx(k),k=1,nbufmx)
c
c        Check value of NBUFX
c
         do k = 1,nbufmx

            if (nbufx(k).gt.mbufle) then
               write (6,*) 'ERROR: NUMBER OF SECONDARY BAFFLE'//
     >               ' ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE'
               write (6,*) 'J = ',k,' NBUFX(J) = ',nbufx(k),
     >                     '  >  MBUFLE = ',mbufle
               call prc( 'ERROR: NUMBER OF SECONDARY BAFFLE'//
     >               ' ELEMENTS'//
     >               ' EXCEEDS ARRAY SIZE')
               write (7,'(a4,i4,a8,i6,a13,i6)') 'J = ',k,
     >                   ' NBUFX(J) = ',nbufx(k),
     >                   '  >  MBUFLE = ',mbufle
               stop
            endif
         end do
      endif
c
      if (buffer(1:6).eq.'RBUFX1') read (4,*)
     >              ((rbufx(i,j),j=1,nbufx(i)),i=1,nbufmx)
      if (buffer(1:6).eq.'ZBUFX1') read (4,*)
     >              ((zbufx(i,j),j=1,nbufx(i)),i=1,nbufmx)
c
      IF (BUFFER(2:7).EQ.'NPOLYP') READ (4,*) NPOLYP
c
c     The polygons are organized as follows:
c     The corner points are ordered 1 to 4 - they count clockwise for
c     JET geometries read straight from the GRID2D file.
c
c     1         2
c     4         3
c
c     With this arrangement of corner points for the cell - the elements
c     1 and 2 form the target when the index IK = 1 - the elements 3 and
c     4 form the target at the other end of the ring when the index
c     IK = NKS(IR). JET grids are calculated for an X-Point up
c     configuration - no matter what the actual arrangement for the
c     shot was. Thus 1 and 2 are associated with the "OUTER" target and
c     3 and 4 with the "INNER" - though this terminology no longer
c     necessarily coincides with the reactor configuration and is being
c     phased out.
c
      IF (BUFFER(2:7).EQ.'NVERTP') READ (4,*) (NVERTP(K),K=1,NPOLYP)
      IF (BUFFER(2:7).EQ.'RVERTP')
     +  READ (4,*) ((RVERTP(J,K),J=1,NVERTP(K)+1),K=1,NPOLYP)
      IF (BUFFER(2:7).EQ.'ZVERTP')
     +  READ (4,*) ((ZVERTP(J,K),J=1,NVERTP(K)+1),K=1,NPOLYP)
      IF (BUFFER(2:7).EQ.'KORPG ') READ (4,*) (DUMMY(K,9),K=1,NP)
      IF (BUFFER(1:6).EQ.'ITAGDV') READ (4,'(12I6)') (ITAGDV(K),K=1,NP)
      IF (BUFFER(1:6).eq.' AREAP') read (4,*) (areap(k),k=1,npolyp)
      GOTO 100
  110 CONTINUE
      WRITE (9,'('' K,ITAG'',/(1X,6I8))') (K,(ITAG(K,L),L=1,5),K=1,NP)
C
      IRTRAP = IRWALL + 1
c

      IKTO  = JPRGT - 1
      IKTI  = JPLFT + 1
c
c     ???? This used to be coded ... it is unclear why ... leave it
c          in but commented out ... just in case a reason becomes
c          apparent later.
c
c     IF (ZXP.LT.0.0) THEN
c       ZXP = -ZXP
c       Z0  = -Z0
c     ENDIF
C
C  DECODE RING AND KNOT NUMBER INFORMATION FROM ITAG USING SKORXYDIV
C
      CALL SKORXYDIV(MAXNKS*MAXNRS,NP,ITAG,
     >            MAXNKS,NR,MAXNRS,KORX,NI,
     >            MAXNRS,NC,MAXNKS,KORY,NJ)
      NRS = NC
      MKS = 0
      DO 120 IR = 1, NRS
        NKS(IR) = NJ(IR)
        MKS = MAX ( MKS, NKS(IR) )
  120 CONTINUE
      IKREF  = NKS(IRSEP)/2 + 1
C
C  LOAD GEOMETRY INFORMATION INTO MATRICES USING KORY,
C  CALCULATE MAXIMUM EXTENT OF GRID AND CORRECT THETAG ARRAY
C
      RMIN = HI
      RMAX = LO
      ZMIN = HI
      ZMAX = LO
      DO 190 IR = 1, NRS
        DO 180 IK = 1, NKS(IR)
          RS(IK,IR)    = DUMMY(KORY(IR,IK),1)
          ZS(IK,IR)    = DUMMY(KORY(IR,IK),2)
          RHOG(IK,IR)  = DUMMY(KORY(IR,IK),3)
          THETAG(IK,IR)= DUMMY(KORY(IR,IK),4)
          HRO(IK,IR)   = DUMMY(KORY(IR,IK),5)
          HTETA(IK,IR) = DUMMY(KORY(IR,IK),6)
c
c          BTS(IK,IR)   = DUMMY(KORY(IR,IK),7)
c
          BTS(IK,IR) = CBPHI * R0 / RS(IK,IR)
c
          KBFS(IK,IR)  = 1.0 / DUMMY(KORY(IR,IK),8)
c slmod begin
          bratio(ik,ir) = dummy(kory(ir,ik),8)
c slmod end
          PSIFL(IK,IR) = DUMMY(KORY(IR,IK),11)
          KORPG(IK,IR) = DUMMY(KORY(IR,IK),9)
          TAGDV(IK,IR) = ITAGDV(KORY(IR,IK))
          RMIN = MIN (RMIN, RS(IK,IR))
          RMAX = MAX (RMAX, RS(IK,IR))
          ZMIN = MIN (ZMIN, ZS(IK,IR))
          ZMAX = MAX (ZMAX, ZS(IK,IR))
  180   CONTINUE
c
c       Comment out this assignment - for cross-field transport
c       to work correctly in the last cell of core rings the
c       theta value has to be able to go to the maximum value
c       for the ring at the cell centre just as with
c       KSS(NKS(IR),IR)  
c
c        IF (IR.LT.IRSEP) THETAG(NKS(IR),IR) = THETAG(1,IR)
c
  190 CONTINUE
C
      WRITE (6,9001) ISHOT,TSLICE,R0,Z0,RXP,ZXP,RMIN,ZMIN,
     >  RMAX,ZMAX,CBPHI,
     >  NP,MKS,IRSEP,IRWALL,IRTRAP,NRS,IKTO,IKREF
c
c     Before processing any baffles - modify the wall to account
c     for any hybrid wall options.
c
      if (ihybrid.gt.0) then
           call rdhybd(ihybrid,nves,rvmod,zvmod,rves,zves)
           call polchg(0,.true.,nves,rvmod,zvmod,rves,zves)
      endif
c
c     BAFFLE processing -
c
c     Redefine the vessel wall to incorporate the baffles if the
c     baffles actually touch the wall. Link together all
c     contiguous baffles and the wall to create a redefined
c     vessel wall that incorporates the additional structures.
c     Change the baffle array contents to retain any
c     non-incorporated baffles.
c
      write (6,*) 'BAFFLES:',nbufle,nbufmx,redefopt
c
      wallredef = 0
c
      if ((nbufle.gt.0.or.nbufmx.gt.0).and.
     >    (cneur.eq.4.or.ctrap.eq.4).and.
     >     redefopt.eq.1) then
c
c        Make copy of original vessel coordinates
c
         nvesorg = nves
c
         do in = 1,nves
            rvesorg(in) = rves(in)
            zvesorg(in) = zves(in)
         end do
c
c        Redefine the vessel wall to include contiguous baffles
c
         call redefine_vessel(nves,
     >                     rves,zves,nbufle,
     >                     rbufle,zbufle,
     >                     nbufmx,nbufx,rbufx,zbufx,
     >                     node_origin,wallredef)
      endif
c
c
c     Correct for counter-clockwise oriented vessel wall
c     description if the option has been set. The vessel
c     wall option intended for use is #4 - however - option
c     #6 will first invert the ordering of the wall segments
c     so that they will be listed clockwise and will thus
c     be compatible with the code and later grid descriptions.
c
      if (cneur.eq.6) then
c
         cneur = 4
c
         do in = 1,nves/2
            rtmp = rves(in)
            ztmp = zves(in)
            rves(in) = rves(nves-in+1)
            zves(in) = zves(nves-in+1)
            rves(nves-in+1) = rtmp
            zves(nves-in+1) = ztmp
         end do
      endif
c
c     Write out the vessel wall coordinates
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         write(6,*) 'Vessel Wall Coordinates:',nves

         do in = 1,nves
c
            write(6,'(i5,4x,3p,bz,f8.2,4x,f8.2)') in,rves(in),zves(in)
c
         end do
c
      endif

C
C-----------------------------------------------------------------------
C     If EDGE plasma profiles have been requested:
C-----------------------------------------------------------------------
C
      CRUN = 'NO EDGE DATA READ'

c
c     Read Edge2D file into Edge2d arrays
c
      if (cre2d.eq.1.or.
     >    CIOPTF.EQ.99.OR.CIOPTG.EQ.99.or.cioptg.eq.90) then
c
          call redge2d(0)
c
      endif
c

c
c      Copy the EDGE2D data into DIVIMP background arrays if required.
c
c
      IF (CIOPTF.EQ.99 .OR. CIOPTG.EQ.99.or. cioptg.eq.90) THEN

c
         do ir = 1,nrs
            do ik = 1,nks(ir)

               KNBS (IK,IR) = e2dnbs(ik,ir)
               KVHS (IK,IR) = e2dvhs(ik,ir)
               KTIBS(IK,IR) = e2dtibs(ik,ir)
               KTEBS(IK,IR) = e2dtebs(ik,ir)
               KES  (IK,IR) = 0.0
C
C              OVERRIDE THE INPUT ELECTRIC FIELD IF THE OPTION IS SET.
C
               IF (OFIELD.EQ.1) then
                  KES(IK,IR) = 0.0
               else
                  KES  (IK,IR) = e2des(ik,ir)
               endif
C
            end do
         end do
c
      end if
c
c
c     Variables required in other parts of the code
c
      IRSEP2 = IRSEP
      IRWALL2 = IRWALL
      IRTRAP2 = IRTRAP
      nrs2    = nrs
      ikt     = ikti
c
c     Calculate the corners of the targets
c
      write (6,*) 'Target corners'
c
c     IK =1 target polygon corners
c
      ik = 2
c
      ir = irwall
      call wrpoly(ik,ir)
c
      ir = irwall-1
      call wrpoly(ik,ir)
c
      ir = irtrap+1
      call wrpoly(ik,ir)
c
      ir = irtrap
      call wrpoly(ik,ir)
c
c     IK = nks(ir) target polygon corners
c
      ir = irtrap
      ik = nks(ir)-1
      call wrpoly(ik,ir)
c
      ir = irtrap+1
      ik = nks(ir)-1
      call wrpoly(ik,ir)
c
      ir = irwall-1
      ik = nks(ir)-1
      call wrpoly(ik,ir)
c
      ir = irwall
      ik = nks(ir)-1
      call wrpoly(ik,ir)
c
c     Write out points across target
c
      write (6,*) 'Targets: Counter-clockwise'
c
c     Outside
c
      ik = 2
      do ir = irwall,irsep,-1
         in = korpg(ik,ir)
         if (in.eq.0) then
            write (6,*) 'No polygon ik,ir:',ik,ir
         else
            write (6,9065) rvertp(2,in),zvertp(2,in)
         endif
      end do
      do ir = nrs,irtrap,-1
         in = korpg(ik,ir)
         if (in.eq.0) then
            write (6,*) 'No polygon ik,ir:',ik,ir
         else
            write (6,9065) rvertp(2,in),zvertp(2,in)
         endif
         if (ir.eq.irtrap+1)
     >           write (6,9065) rvertp(1,in),zvertp(1,in)
      end do
c
c     Inside
c
      do ir = irtrap,nrs
         ik = nks(ir) -1
         in = korpg(ik,ir)
         if (in.eq.0) then
            write (6,*) 'No polygon ik,ir:',ik,ir
         else
            write (6,9065) rvertp(4,in),zvertp(4,in)
         endif
      end do
      do ir = irsep,irwall
         ik = nks(ir) -1
         in = korpg(ik,ir)
         if (in.eq.0) then
            write (6,*) 'No polygon ik,ir:',ik,ir
         else
            write (6,9065) rvertp(4,in),zvertp(4,in)
         endif
         if (ir.eq.irwall-1)
     >           write (6,9065) rvertp(3,in),zvertp(3,in)
      end do
c slmod begin
      in = 0
      DO ir = 1, nrs
        DO ik = 1, nks(ir)
          IF (korpg(ik,ir).LT.MAXNKS*MAXNRS) in = MAX(korpg(ik,ir),in)
        ENDDO
      ENDDO
      vpolmin = (MAXNKS*MAXNRS - in) / 2 + in
      vpolyp  = vpolmin

c...  Assign PSIn values for the targets:
      psitarg = 0.0
      DO ir = 1, nrs
        psitarg(ir,2) = psifl(1      ,ir)       
        psitarg(ir,1) = psifl(nks(ir),ir)       
      ENDDO      

      CALL OutputData(85,'End of RJET')

c      z0  = -z0
c      zxp = -zxp
c      DO id = 1, in
c        zvertp(1:4,id) = -zvertp(1:4,id)
c      ENDDO 
c      CALL OutputData(86,'End of RJET, really this time')
c      CALL DumpGrid('End of RJET')
c slmod end
c
c
c     END RJET
c
      return
C
 9999 IERR = 1
      WRITE (7,'(1X,A,3(/1X,A))')
     >  'TAUIN1: ERROR READING GEOMETRY DATA :-',MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
C
 9001 FORMAT(/1X,'GEOMETRY DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
     >  /5X,'TOROIDAL FIELD     BPHI =  ',F8.3,
     >  /5X,'NO OF POINTS NP     =',I6,5X,'MAX NO. ROWS MKS    =',I6,
     >  /5X,'SEPARATRIX   IRSEP  =',I6,5X,'WALL         IRWALL =',I6,
     >  /5X,'FIRST TRAP   IRTRAP =',I6,5X,'NO OF RINGS  NRS    =',I6,
     >  /5X,'K SPLIT PT   IKTO   =',I6,5X,'K REF POINT  IKREF  =',I6)
 9002 FORMAT(/1X,' IK IR    R      Z     BPH',
     >  'I     TEB     TIB     NB      E1      VB        S  ',
     >  'BTOT/BTHETA    FEG1    FIG1',/1X,131('-'))
 9003 FORMAT(1X,2I3,3F7.3,2f8.1,1P,E8.1,0P,2A9,2F8.2,3X,2A9)
 9006 FORMAT(/1X,'TAUIN1: AREA OF SOL+TRAP =',F9.6,',  MAIN P =',F9.6,/)
 9011 FORMAT(/1X,'EDGE PLASMA DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
     >  /5X,'RUN                 =',A,
     >  /5X,'PLASMA       HMASS  =',I6,5X,'IMPURITY     ZMASS  =',I6,
     >  /5X,'FIRST RING   IRCENT =',I6,5X,'LOST INN RGS NINOMP =',I6,
     >  /5X,'LOST OUT RGS NINOWA =',I6)
 9012 FORMAT(1P,6E12.4)
 9022 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
     >  ' IROUT  INDIST OUTDIST BACDIST FORDIST VOLUME  INPROB',
     >  '  POLDIST',/1X,131('-'))
 9023 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,2X,4F8.3,F8.5,F8.4,F8.3)
 9032 FORMAT(/1X,' IK IR    R      Z       NH0   TAUIZ0  CHPROB0 ',
     >  2('TAUIZ',I1,'  TAUEI',I1,'  TAURC',I1,'  CHPROB',I1,
     >  ' RCFRAC',I1,' '),/1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,2X,1P,3E8.1,2(1P,4E8.1,0P,F8.4))
 9040 FORMAT(
     >   5X,'NO OF R PTS       NXS   =  ',I5,',   DELTAR =',F9.5,
     >  /5X,'NO OF Z PTS       NYS   =  ',I5,',   DELTAZ =',F9.5,
     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
     >  /5X,'INNER TARGET PTS  NDSIN =  ',I5,
     >  /5X,'TOTAL TARGET PTS  NDS   =  ',I5)
 9042 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
     >  ' IROUT  ID   INTHETA OUTTHETA THETA   DIST',
     >  /1X,131('-'))
 9043 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,I7,2X,3F8.2,F8.5)
 9044 FORMAT(I7,E14.0,A40)
 9050 FORMAT(I5,I5,50F8.3)
 9065 format(1p,6(1x,e11.4))
c
      end
c
c
c
      subroutine wrpoly(ik,ir)
      implicit none
      integer ik,ir
c
      include 'params'
      include 'cgeom'
c
c     WRPOLY: Writes out the polygon at coordinate ik,ir - if
c             one exists.
c
      integer in,i
c
      in = korpg(ik,ir)
c
      if (in.eq.0) then
         write (6,*) 'No polygon for ik,ir:',ik,ir
      else
         write (6,*) 'Polygon:',ik,ir
         do i = 1,nvertp(in)
            write (6,*) 'Corner:',i,' R,Z=',
     >              rvertp(i,in),zvertp(i,in)
         end do
      endif
c
      return
      end
c
c
c
      subroutine rasdex
      IMPLICIT none
c
c     As with the above routine ... this is intended to read ASDEX
c     grids into the standard DIVIMP arrays. The original code
c     upon which this routine was based was written at Julich.
c     Unfrotunately it was hard-coded for a specific geometry file ...
c     this needs to be corrected, when more information about the
c     characteristics of ASDEX grids is available ... this upgrade
c     will be made.
c
c     David Elder, June 18, 1993
c
C
C  *********************************************************************
C  *                                                                   *
C  *  RASDEX  :  INITIALISING GEOMETRY ETC
C  *                                                                   *
C  *            Chris Farrell  (Hunterskil)  February 1989             *
C  *                                                                   *
C  *            James Spence   (Tessella)    July     1991             *
C  *                                                                   *
C  * updates     : ??/07/90 --- Allow rings to be missing from wall    *
C  *               08/07/91 --- Allow virtual rings (iring=1,nxw,nxw+1)*
C  *                        --- Allow EDGE2D input                     *
C  *               15/08/91 --- Changes for temp. grad opt 2           *
C  *                                                                   *
C  *********************************************************************
C
      REAL    VFLUID , XJI , XJ , XJF
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cioniz'
      include 'reader'
      include 'dynam5'
      CHARACTER MESAGE*72,C(10)*9,FACTOR*9,FORM*72
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IDUMMY((MAXNKS+2)*MAXNRS),NJ(MAXNRS),NI(MAXNKS+2)
      INTEGER ISEP,IWALL,ITRAP,NSOL,NRSOL,NINOWA,IREF,NZ,K0,IK0,IKMID
      LOGICAL FLUSH
      REAL    DUMMY((MAXNKS+2)*MAXNRS,5),HMASS,ZMASS,S,DELTAL,P,FACT
      REAL    BTHETA,SPERP,BEST,DSQ,R,Z,THIN,THOUT,DELTAR,DELTAZ
      REAL    MU,SMAX,SMID,SDIFF
c
      real brs(maxnks,maxnrs),bzs(maxnks,maxnrs)
C
c     Set variable Crun
c
      Crun = 'ASDEX GRID'

c
C
C-----------------------------------------------------------------------
C     Extract Geometry details from first file
C-----------------------------------------------------------------------
C
      REWIND (10)
C
      FLUSH  = .FALSE.
      MESAGE = 'Extracting FLUSH header details'
  100 READ (10,'(A72)',ERR=9999,END=9999) BUFFER
C     WRITE (9,'(1X,A72,1X,A6)') BUFFER,'TAUINAS'
      IF (FLUSH) THEN
        CALL FORMMR(BUFFER(35:72),FORM,1)
        READ (BUFFER(35:72),FORM) TSLICE
        FLUSH = .FALSE.
      ENDIF
      IF (BUFFER(3:12).EQ.' FLUSH *  ') THEN
C        CALL FORMI(BUFFER(35:72),FORM)
C        READ (BUFFER(35:72),FORM) ISHOT
        READ (BUFFER(35:72),*) ISHOT
        FLUSH = .TRUE.
      ENDIF
      IF (BUFFER(3:12).EQ.'      RW. ') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,2)
C          READ (BUFFER(20:72),FORM) RMIN,RMAX
          READ (BUFFER(20:72),*) RMIN,RMAX
      ENDIF
      IF (BUFFER(3:12).EQ.'      ZW. ') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,2)
C          READ (BUFFER(20:72),FORM) ZMIN,ZMAX
          READ (BUFFER(20:72),*) ZMIN,ZMAX
      ENDIF
      IF (BUFFER(3:12).EQ.'      XP0 ') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,1)
C          READ (BUFFER(20:72),FORM) R0
          READ (BUFFER(20:72),*) R0
      ENDIF
      IF (BUFFER(3:12).EQ.'      YP0 ') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,1)
C          READ (BUFFER(20:72),FORM) Z0
          READ (BUFFER(20:72),*) Z0
      ENDIF
      IF (BUFFER(3:12).EQ.' BT       ') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,1)
C          READ (BUFFER(20:72),FORM) CBPHI
          READ (BUFFER(20:72),*) CBPHI
      ENDIF
      IF (BUFFER(3:12).EQ.'      RSEP') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,1)
C          READ (BUFFER(20:30),FORM) RXP
          READ (BUFFER(20:30),*) RXP
      ENDIF
      IF (BUFFER(3:12).EQ.'      ZSEP') THEN
C          CALL FORMMR(BUFFER(20:72),FORM,1)
C          READ (BUFFER(20:30),FORM) ZXP
          READ (BUFFER(20:30),*) ZXP
      ENDIF
      IF (BUFFER(3:15).NE.'   K IRING   ') GOTO 100
C
      MESAGE = 'Expecting tables of IK,IR,R,Z,BR,BZ values'
      MKS = 0
  150 IK  = 0
  160 READ (10,'(A72)',ERR=9999,END=170) BUFFER
      WRITE (9,'(1X,A72,1X,A6,2i4)') BUFFER,'TAUINAS',ik,ir
      IF (BUFFER(3:15).EQ.'             ') GOTO 150
      IF (BUFFER(3:15).EQ.'   K IRING   ') GOTO 150
      IK = IK + 1
C      CALL FORM2I4R(BUFFER,FORM)
C      READ (BUFFER,FORM) K,IR,RS(IK,IR),ZS(IK,IR),BRS(IK,IR),BZS(IK,IR)
      READ (BUFFER,*) K,IR,RS(IK,IR),ZS(IK,IR),BRS(IK,IR),BZS(IK,IR)
      IF (IK.GT.MAXNKS) WRITE (6,*) ' Error!  ik=',IK
      IF (IR.GT.MAXNRS) WRITE (6,*) ' Error!  ir=',IR
      NKS(IR) = IK
      MKS = MAX (MKS, NKS(IR))
      NRS = IR
      GOTO 160
  170 CONTINUE
C
C---- FIND ACTUAL SEPARATRIX , WALL & 1ST PRIVATE REGION RINGS
C
      IRSEP  = 0
      IRWALL = NRS
C
      DO 175 IR = 2 , NRS
         IF( IRSEP.NE.0 .AND. NKS(IR).LT.NKS(IR-1) ) IRWALL = IR - 1
         IF( IRSEP.EQ.0 .AND. NKS(IR).GT.NKS(IR-1) ) IRSEP  = IR
  175 CONTINUE
      IRTRAP = IRWALL + 1
      IRSEP2 = IRSEP
      IRWALL2 = IRWALL
      IRTRAP2 = IRTRAP
      nrs2    = nrs
C
      NRSOL  = IRWALL - IRSEP + 1
C
C---- Convert cm to m, Gauss to Tesla units; calculate toroidal parts.
C
      R0   = 0.01 * R0
      Z0   = 0.01 * Z0
      RXP  = 0.01 * RXP
      ZXP  = 0.01 * ZXP
      RMIN = 0.01 * RMIN
      ZMIN = 0.01 * ZMIN
      RMAX = 0.01 * RMAX
      ZMAX = 0.01 * ZMAX
      CBPHI = 0.0001 * CBPHI
      WRITE (6,9001) ISHOT,TSLICE,R0,Z0,RXP,ZXP,RMIN,ZMIN,
     >  RMAX,ZMAX,CBPHI
C
      RMIN = HI
      RMAX = LO
      ZMIN = HI
      ZMAX = LO
      DO 190 IR = 1, NRS
        DO 180 IK = 1, NKS(IR)
          RS(IK,IR)  = 0.01 * RS(IK,IR)
          ZS(IK,IR)  = 0.01 * ZS(IK,IR)
          BRS(IK,IR) = 0.0001 * BRS(IK,IR)
          BZS(IK,IR) = 0.0001 * BZS(IK,IR)
          BTS(IK,IR) = CBPHI * R0 / RS(IK,IR)
          RMIN = MIN (RMIN, RS(IK,IR))
          RMAX = MAX (RMAX, RS(IK,IR))
          ZMIN = MIN (ZMIN, ZS(IK,IR))
          ZMAX = MAX (ZMAX, ZS(IK,IR))
  180   CONTINUE
  190 CONTINUE

C
C-----------------------------------------------------------------------
C     Extract Densities, Temperatures etc from second file
C-----------------------------------------------------------------------
C
      if (cioptf.eq.99.or.cioptg.eq.99.or.cioptg.eq.90) then
         REWIND (11)
         MESAGE = 'EXPECTING TABLES OF IK,IR,DNIB,TIB,TEB,UPB'
  151    IK  = 0
  161    READ (11,'(A72)',ERR=9999,END=171) BUFFER
C        WRITE (9,'(1X,A72,1X,A6)') BUFFER,'TAUINAS'
         IF (BUFFER(3:15).EQ.'             ') GOTO 151
         IF (BUFFER(3:15).EQ.'   K IRING   ') GOTO 151
         IK = IK + 1
C         CALL FORE2I4R(BUFFER,FORM)
         READ (BUFFER,*) K,IR,KNBS(IK,IR),
     >                   KTIBS(IK,IR),KTEBS(IK,IR),KVHS(IK,IR)
         KTEBS(IK,IR) = KTEBS(IK,IR) / 1.6022E-19
         KTIBS(IK,IR) = KTIBS(IK,IR) / 1.6022E-19
         KVHS(IK,IR) = KVHS(IK,IR) *100.
         KES(IK,IR) = 0.0
         IF (IK.GT.MAXNKS) WRITE (6,*) ' Error!  ik=',IK
         IF (IR.GT.MAXNRS) WRITE (6,*) ' Error!  ir=',IR
         NKS(IR) = IK
         MKS = MAX (MKS, NKS(IR))
         NRS = IR
         GOTO 161
171      CONTINUE
      endif
C
      IKT = 25
      IKTO = 25
      IKTI = nks(irsep) - ikt +1
      IKREF  = NKS(IRSEP)/2 + 1
      NP = K
      NINOWA = 0
      HMASS = 0.0
      ZMASS = 0.0
      WRITE (6,9011) ISHOT,TSLICE,CRUN,NP,MKS,IRSEP,IRWALL,
     >  IRTRAP,IRCENT,IKT,NRS,NINT(HMASS),NINT(ZMASS),IKREF,NINOWA
C
C-----------------------------------------------------------------------
C     Extract rho and theta data for later volume calculations ...
C     (Note that SOL and Trap rings have two extra points in file 12)
C-----------------------------------------------------------------------
C
      REWIND (14)
C
      NP1 = 0
      MESAGE = 'Extracting rho and theta data'
  152 IK  = 0
  162 READ (14,'(A72)',ERR=9999,END=172) BUFFER
      IF (BUFFER(3:15).EQ.'             ') GOTO 152
      IF (BUFFER(3:15).EQ.'   K IRING   ') GOTO 152
      IK = IK + 1
C      CALL FORE2I4R(BUFFER,FORM)
      READ (BUFFER,*) K,IR,KVOLS(IK,IR)
      GOTO 162
  172 CONTINUE
  295 CONTINUE

C
C-----------------------------------------------------------------------
C     Calculate distances along contours S, magnetic field factors etc
C-----------------------------------------------------------------------
C
      DO 320 IR = 1, NRS
        DO 310 IK = 1, NKS(IR)
          BTHETA = SQRT (BRS(IK,IR)**2 + BZS(IK,IR)**2)
          KBFS(IK,IR) = SQRT (BTHETA**2 + BTS(IK,IR)**2) / BTHETA
  310   CONTINUE
  320 CONTINUE
C
C
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,A,3(/1X,A))')
     >  'RASDEX: ERROR READING GEOMETRY DATA :-',MESAGE,
     >  'Last line read :-',BUFFER
      RETURN
C
 9001 FORMAT(/1X,'Geometry data for Shot',I6,',  Time',F8.3,' :-',
     >  /5X,'Centre      (Ro  ,Zo  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X point     (Rxp ,Zxp ) = (',F8.3,',',F8.3,')',
     >  /5X,'Lower left  (Rmin,Zmin) = (',F8.3,',',F8.3,')',
     >  /5X,'Upper right (Rmax,Zmax) = (',F8.3,',',F8.3,')',
     >  /5X,'Toroidal Field     Bphi =  ',F8.3)
 9002 FORMAT(/1X,' ik ir    R      Z     BR     BZ    Bph',
     >  'i     TeB     TiB     NB      E1      VB        S  ',
     >  'Btot/Btheta    FeG1    FiG1',/1X,131('-'))
 9003 FORMAT(1X,2I3,5F7.3,2F8.1,1P,E8.1,0P,2A9,2F8.2,3X,2A9)
 9006 FORMAT(/1X,'TAUINAS: AREA OF SOL+TRAP =',F9.6,',
     >       MAIN P =',F9.6,/)
 9011 FORMAT(/1X,'Equilibrium data for Shot',I6,',  Time',F8.3,' :-',
     >  /5X,'Run                 =',A,
     >  /5X,'No of points NP     =',I6,5X,'Max no. rows MKS    =',I6,
     >  /5X,'Separatrix   IRSEP  =',I6,5X,'Wall         IRWALL =',I6,
     >  /5X,'First trap   IRTRAP =',I6,5X,'First ring   IRCENT =',I6,
     >  /5X,'K split pt   IKT    =',I6,5X,'No of rings  NRS    =',I6,
     >  /5X,'Plasma       HMASS  =',I6,5X,'Impurity     ZMASS  =',I6,
     >  /5X,'K ref point  IKREF  =',I6,5x,'Lost rings   NINOWA =',I6)
 9012 FORMAT(1P,6E12.4)
 9013 FORMAT(/1X,'Volume data :-  Number of pts NP1 =',I5)
 9022 FORMAT(/1X,' ik ir    R      Z   ikin irin ikout',
     >  ' irout  InDist OutDist BacDist ForDist Volume  InProb',
     >  '  PolDist',/1X,131('-'))
 9023 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,2X,4F8.3,F8.5,F8.4,F8.3)
 9032 FORMAT(/1X,' ik ir    R      Z       NH0   tauiz0  chprob0 ',
     >  2('tauiz',I1,'  tauei',I1,'  taurc',I1,'  chprob',I1,
     >  ' rcfrac',I1,' '),/1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,2X,1P,3E8.1,2(1P,4E8.1,0P,F8.4))
 9040 FORMAT(
     >   5X,'No of R pts       NXS   =  ',I5,',   DeltaR =',F9.5,
     >  /5X,'No of Z pts       NYS   =  ',I5,',   DeltaZ =',F9.5,
     >  /5X,'Centre      (Ro  ,Zo  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X point     (Rxp ,Zxp ) = (',F8.3,',',F8.3,')',
     >  /5X,'Lower left  (Rmin,Zmin) = (',F8.3,',',F8.3,')',
     >  /5X,'Upper right (Rmax,Zmax) = (',F8.3,',',F8.3,')',
     >  /5X,'Inner Target Pts  NDSIN =  ',I5,
     >  /5X,'Total Target Pts  NDS   =  ',I5)
 9042 FORMAT(/1X,' ik ir    R      Z   ikin irin ikout',
     >  ' irout  id   InTheta OutTheta Theta   Dist',
     >  /1X,131('-'))
 9043 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,I7,2X,3F8.2,F8.5)
 9044 FORMAT(I7,E14.0,A40)
      END
c
c
c
      subroutine riter
c
c     This subroutine is intended to read and convert ITER grids.
c
      IMPLICIT none
C
C  *********************************************************************
C  *                                                                   *
C  *  RITER  :  INITIALISING GEOMETRY FOR ITER GRIDS                   *
C  *                                                                   *
C  *********************************************************************
C
      REAL    VFLUID , XJI , XJ , XJF
      REAL RICHTABX,RICHTABY,RICHTCDX,RICHTCDY,MUE
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cioniz'
      include 'reader'
      include 'dynam5'
      CHARACTER MESAGE*72,C(10)*9,FACTOR*9,FORM*72
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IDUMMY((MAXNKS+2)*MAXNRS),NJ(MAXNRS),NI(MAXNKS+2)
C     INTEGER ITAG((MAXNKS+2)*MAXNRS,5),KORY(MAXNRS,MAXNKS)
      INTEGER ISEP,IWALL,ITRAP,NSOL,NRSOL,NINOWA,IREF,NZ,K0,IK0,IKMID
      LOGICAL FLUSH
      REAL    DUMMY((MAXNKS+2)*MAXNRS,5),HMASS,ZMASS,S,DELTAL,P,FACT
      REAL    BTHETA,SPERP,BEST,DSQ,R,Z,THIN,THOUT,DELTAR,DELTAZ
      REAL    MU,SMAX,SMID,SDIFF
c
      real brs(maxnks,maxnrs), bzs(maxnks,maxnrs)
C
C-----------------------------------------------------------------------
C     Extract Geometry details from first file
C-----------------------------------------------------------------------
C
      REWIND (10)
C
c     Set variable Crun
c
      Crun = 'ITER GRID'
c
      FLUSH  = .FALSE.
      MESAGE = 'Extracting FLUSH header details'
  100 READ (10,'(A80)',ERR=9999,END=9999) BUFFER1
      IF (FLUSH) THEN
C        CALL FORMMR(BUFFER1(35:72),FORM,1)
C        READ (BUFFER1(35:72),FORM) TSLICE
        READ (BUFFER1(35:72),*) TSLICE
        FLUSH = .FALSE.
      ENDIF
      IF (BUFFER1(3:12).EQ.' FLUSH *  ') THEN
C        CALL FORMI(BUFFER1(35:72),FORM)
C        READ (BUFFER1(35:72),FORM) ISHOT
        READ (BUFFER1(35:72),*) ISHOT
        FLUSH = .TRUE.
      ENDIF
      IF (BUFFER1(3:12).EQ.'      RW. ') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,2)
C          READ (BUFFER1(20:72),FORM) RMIN,RMAX
          READ (BUFFER1(20:72),*) RMIN,RMAX
      ENDIF
      IF (BUFFER1(3:12).EQ.'      ZW. ') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,2)
C          READ (BUFFER1(20:72),FORM) ZMIN,ZMAX
          READ (BUFFER1(20:72),*) ZMIN,ZMAX
      ENDIF
      IF (BUFFER1(3:12).EQ.'      XP0 ') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,1)
C          READ (BUFFER1(20:72),FORM) R0
          READ (BUFFER1(20:72),*) R0
      ENDIF
      IF (BUFFER1(3:12).EQ.'      YP0 ') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,1)
C          READ (BUFFER1(20:72),FORM) Z0
          READ (BUFFER1(20:72),*) Z0
      ENDIF
      IF (BUFFER1(3:12).EQ.' BT       ') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,1)
C          READ (BUFFER1(20:72),FORM) CBPHI
          READ (BUFFER1(20:72),*) CBPHI
      ENDIF
      IF (BUFFER1(3:12).EQ.'      RSEP') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,1)
C          READ (BUFFER1(20:30),FORM) RXP
          READ (BUFFER1(20:30),*) RXP
      ENDIF
      IF (BUFFER1(3:12).EQ.'      ZSEP') THEN
C          CALL FORMMR(BUFFER1(20:72),FORM,1)
C          READ (BUFFER1(20:30),FORM) ZXP
          READ (BUFFER1(20:30),*) ZXP
      ENDIF
      IF (BUFFER1(3:15).NE.'   K IRING   ') GOTO 100
C
      MESAGE = 'Expecting tables of IK,IR,R,Z,BR,BZ values'
      MKS = 0
  150 IK  = 0
  160 READ (10,'(A80)',ERR=9999,END=170) BUFFER1
      IF (BUFFER1(3:15).EQ.'             ') GOTO 150
      IF (BUFFER1(3:15).EQ.'   K IRING   ') GOTO 150
      IK = IK + 1
C      CALL FORM2I4R(BUFFER1,FORM)
C      READ (BUFFER1,FORM) K,IR,RS(IK,IR),ZS(IK,IR),BRS(IK,IR),BZS(IK,IR
      READ (BUFFER1,*) K,IR,RS(IK,IR),ZS(IK,IR),BRS(IK,IR),BZS(IK,IR)
      IF (IK.GT.MAXNKS) WRITE (6,*) ' Error!  ik=',IK
      IF (IR.GT.MAXNRS) WRITE (6,*) ' Error!  ir=',IR
      NKS(IR) = IK
      MKS = MAX (MKS, NKS(IR))
      NRS = IR
      GOTO 160
  170 CONTINUE
C
C---- FIND ACTUAL SEPARATRIX , WALL & 1ST PRIVATE REGION RINGS
C
      IRSEP  = 0
      IRWALL = 62
      IRWALL2 = 31
      NRS2 = 46
C
      DO 175 IR = 2 , NRS
         IF( IRSEP.NE.0 .AND. NKS(IR).LT.NKS(IR-1) ) IRWALL = IR - 1
         IF( IRSEP.EQ.0 .AND. NKS(IR).GT.NKS(IR-1) ) IRSEP  = IR
  175 CONTINUE
      IRTRAP = IRWALL + 1
      IRTRAP2 = IRWALL2 + 1
      IRSEP2 = 47
C
      NRSOL  = IRWALL - IRSEP + 1
c
      write (6,*) 'Ring indices:', irsep,irwall2,irtrap2,nrs2,
     >             irsep2,irwall,irtrap,nrs
C
C---- Convert cm to m, Gauss to Tesla units; calculate toroidal parts.
C
      R0   = 0.01 * R0
      Z0   = 0.01 * Z0
      RXP  = 0.01 * RXP
      ZXP  = 0.01 * ZXP
      RMIN = 0.01 * RMIN
      ZMIN = 0.01 * ZMIN
      RMAX = 0.01 * RMAX
      ZMAX = 0.01 * ZMAX
      CBPHI = 0.0001 * CBPHI
      WRITE (6,9001) ISHOT,TSLICE,R0,Z0,RXP,ZXP,RMIN,ZMIN,
     >  RMAX,ZMAX,CBPHI
C
      RMIN = HI
      RMAX = LO
      ZMIN = HI
      ZMAX = LO
      DO 190 IR = 1, NRS
        DO 180 IK = 1, NKS(IR)
          RS(IK,IR)  = 0.01 * RS(IK,IR)
          ZS(IK,IR)  = 0.01 * ZS(IK,IR)
          BRS(IK,IR) = 0.0001 * BRS(IK,IR)
          BZS(IK,IR) = 0.0001 * BZS(IK,IR)
          BTS(IK,IR) = CBPHI * R0 / RS(IK,IR)
          RMIN = MIN (RMIN, RS(IK,IR))
          RMAX = MAX (RMAX, RS(IK,IR))
          ZMIN = MIN (ZMIN, ZS(IK,IR))
          ZMAX = MAX (ZMAX, ZS(IK,IR))
  180   CONTINUE
  190 CONTINUE
C
C-----------------------------------------------------------------------
C     Extract Densities, Temperatures etc from second file
C-----------------------------------------------------------------------
C
      if (cioptf.eq.99.or.cioptg.eq.99.or.cioptg.eq.90) then
        REWIND (11)
        MESAGE = 'EXPECTING TABLES OF IK,IR,DNIB,TIB,TEB,UPB'
  151   IK  = 0
  161   READ (11,'(A72)',ERR=9999,END=171) BUFFER1
        IF (BUFFER1(3:15).EQ.'             ') GOTO 151
        IF (BUFFER1(3:15).EQ.'   K IRING   ') GOTO 151
        IK = IK + 1
C        CALL FORE2I4R(BUFFER1,FORM)
C        READ (BUFFER1,FORM) K,IR,KNBS(IK,IR),
C     >                   KTIBS(IK,IR),KTEBS(IK,IR),KVHS(IK,IR)
        READ (BUFFER1,*) K,IR,KNBS(IK,IR),
     >                   KTIBS(IK,IR),KTEBS(IK,IR),KVHS(IK,IR)
        KTEBS(IK,IR) = KTEBS(IK,IR) / 1.6022E-19
        KTIBS(IK,IR) = KTIBS(IK,IR) / 1.6022E-19
c
c        KVHS(IK,IR) = KVHS(IK,IR) *100.
c
c       Do not change units of kvhs for ITER run
c
        KES(IK,IR) = 0.0
        IF (IK.GT.MAXNKS) WRITE (6,*) ' Error!  ik=',IK
        IF (IR.GT.MAXNRS) WRITE (6,*) ' Error!  ir=',IR
        NKS(IR) = IK
        MKS = MAX (MKS, NKS(IR))
        NRS = IR
        GOTO 161
171     CONTINUE
      endif
C
      IKT = 15
      IKT1 = 15
      IKT2 = 16
c      IKREF  = NKS(IRSEP)/2 + 1
      NP = K
      NINOWA = 0
      HMASS = 0.0
      ZMASS = 0.0
      WRITE (6,9011) ISHOT,TSLICE,CRUN,NP,MKS,IRSEP,IRWALL,
     >  IRTRAP,IRCENT,IKT,NRS,NINT(HMASS),NINT(ZMASS),IKREF,NINOWA
C
C-----------------------------------------------------------------------
C     Extract rho and theta data for later volume calculations ...
C     (Note that SOL and Trap rings have two extra points in file 12)
C-----------------------------------------------------------------------
C
      REWIND (14)
C
      NP1 = 0
      MESAGE = 'Extracting rho and theta data'
  152 IK  = 0
  162 READ (14,'(A72)',ERR=9999,END=172) BUFFER1
      IF (BUFFER1(3:15).EQ.'             ') GOTO 152
      IF (BUFFER1(3:15).EQ.'   K IRING   ') GOTO 152
      IK = IK + 1
C      CALL FORE2I4R(BUFFER1,FORM)
C      READ (BUFFER1,FORM) K,IR,KVOLS(IK,IR)
      READ (BUFFER1,*) K,IR,KVOLS(IK,IR)
      GOTO 162
  172 CONTINUE
  295 CONTINUE
C
C-----------------------------------------------------------------------
C     Calculate distances along contours S, magnetic field factors etc
C-----------------------------------------------------------------------
C
      DO 320 IR = 1, NRS
        DO 310 IK = 1, NKS(IR)
          BTHETA = SQRT (BRS(IK,IR)**2 + BZS(IK,IR)**2)
          KBFS(IK,IR) = SQRT (BTHETA**2 + BTS(IK,IR)**2) / BTHETA
  310   CONTINUE
  320 CONTINUE
c
c     End
c
      RETURN
C
 9999 IERR = 1
      WRITE (7,'(1X,A,3(/1X,A))')
     >  'RITER: ERROR READING GEOMETRY DATA :-',MESAGE,
     >  'Last line read :-',BUFFER1
      RETURN
C
 9001 FORMAT(/1X,'Geometry data for Shot',I6,',  Time',F8.3,' :-',
     >  /5X,'Centre      (Ro  ,Zo  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X point     (Rxp ,Zxp ) = (',F8.3,',',F8.3,')',
     >  /5X,'Lower left  (Rmin,Zmin) = (',F8.3,',',F8.3,')',
     >  /5X,'Upper right (Rmax,Zmax) = (',F8.3,',',F8.3,')',
     >  /5X,'Toroidal Field     Bphi =  ',F8.3)
 9002 FORMAT(/1X,' ik ir    R      Z     BR     BZ    Bph',
     >  'i     TeB     TiB     NB      E1      VB        S  ',
     >  'Btot/Btheta    FeG1    FiG1',/1X,131('-'))
 9003 FORMAT(1X,2I3,5F7.3,2F8.1,1P,E8.1,0P,2A9,2F8.2,3X,2A9)
 9006 FORMAT(/1X,'TAUINAS: AREA OF SOL+TRAP =',F9.6,',
     >       MAIN P =',F9.6,/)
 9011 FORMAT(/1X,'Equilibrium data for Shot',I6,',  Time',F8.3,' :-',
     >  /5X,'Run                 =',A,
     >  /5X,'No of points NP     =',I6,5X,'Max no. rows MKS    =',I6,
     >  /5X,'Separatrix   IRSEP  =',I6,5X,'Wall         IRWALL =',I6,
     >  /5X,'First trap   IRTRAP =',I6,5X,'First ring   IRCENT =',I6,
     >  /5X,'K split pt   IKT    =',I6,5X,'No of rings  NRS    =',I6,
     >  /5X,'Plasma       HMASS  =',I6,5X,'Impurity     ZMASS  =',I6,
     >  /5X,'K ref point  IKREF  =',I6,5x,'Lost rings   NINOWA =',I6)
 9012 FORMAT(1P,6E12.4)
 9013 FORMAT(/1X,'Volume data :-  Number of pts NP1 =',I5)
 9022 FORMAT(/1X,' ik ir    R      Z   ikin irin ikout',
     >  ' irout  InDist OutDist BacDist ForDist Volume  InProb',
     >  '  PolDist',/1X,131('-'))
 9023 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,2X,4F8.3,F8.5,F8.4,F8.3)
 9032 FORMAT(/1X,' ik ir    R      Z       NH0   tauiz0  chprob0 ',
     >  2('tauiz',I1,'  tauei',I1,'  taurc',I1,'  chprob',I1,
     >  ' rcfrac',I1,' '),/1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,2X,1P,3E8.1,2(1P,4E8.1,0P,F8.4))
 9040 FORMAT(
     >   5X,'No of R pts       NXS   =  ',I5,',   DeltaR =',F9.5,
     >  /5X,'No of Z pts       NYS   =  ',I5,',   DeltaZ =',F9.5,
     >  /5X,'Centre      (Ro  ,Zo  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X point     (Rxp ,Zxp ) = (',F8.3,',',F8.3,')',
     >  /5X,'Lower left  (Rmin,Zmin) = (',F8.3,',',F8.3,')',
     >  /5X,'Upper right (Rmax,Zmax) = (',F8.3,',',F8.3,')',
     >  /5X,'Inner Target Pts  NDSIN =  ',I5,
     >  /5X,'Total Target Pts  NDS   =  ',I5)
 9042 FORMAT(/1X,' ik ir    R      Z   ikin irin ikout',
     >  ' irout  id   InTheta OutTheta Theta   Dist',
     >  /1X,131('-'))
 9043 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,I7,2X,3F8.2,F8.5)
 9044 FORMAT(I7,E14.0,A40)
      END
c     
c     
c     
      subroutine raug
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cedge2d'
c     slmod begin
      include 'slcom'
c...  temp
      include 'pindata'

c...  temp
      CHARACTER title*174,desc*1024,job*72,equil*60
      REAL      facta(-1:MAXIZS),factb(-1:MAXIZS)

      INTEGER i1
c     slmod end
c     
c     The purpose of this code is to read in ASDEX UPGRADE
c     geometry files and process them for use in DIVIMP. The
c     geometry file has a fixed format, it may change at some
c     time, in which case the following code will need to be adapted.
c     
c     The information included in the geometry file is:
c     1) Corner points of the cell
c     2) Centre point of the cell
c     3) Magnetic field ratio ... it is uncertain (right now
c     2pm Aug 5/93) as to whether this is Bpol/Btor OR Bpol/Btot
c     For the time being I am assuming that it is the quantity
c     needed for target flux calculations.
c     
c     DIVIMP will use its own routines to estimate the cell areas and
c     volumes. (Until this data can be made available from the ASDEX
c     grid files.)
c     
c     This routine is implemented using fixed format statements which
c     should correspond to the current structure of the geometry file.
c     The examples of geometry files currently available share the
c     same format.
c     
c     Furthermore, the grid structure is identical in each case, even
c     though the coordinates vary. Thus the following "rules" are used
c     to decode the grid information.
c     
c     1) The zero and 25 rings are boundary rings for the fluid code
c     2) The zero and 97 elements of each ring are boundary points for
c     the targets.
c     3) Rings 0 to 7 contain the specifications for the trap and
c     core plasma rings. Elements 0 to 24 and 73 to 97 are in the
c     trapped plasma. Elements 25 to 72 are in the core ... starting
c     at the point just above the X-point. If 25 and 72 do not
c     lie on top of one another then an additional point matching
c     25 is generated and added to the end so that it will form
c     a closed ring.
c     4) The other rings represent the SOL region.
c     
c     
c     David Elder, Aug 5, 1993
c     
c----------------------------------------------------------------------
c     
c     Asdex grid - characteristic numbers
c     
c     Need to add one to the cut points and maximums because
c     the AUG grids are indexed from zero.
c     6+1  = 7
c     24+1 = 25
c     73+1 = 74
c     25+1 = 26
c     97+1 = 98
c     
c     
c     Moved into common block COMTOR so that they may be read from the
c     input file - this will generalize the reading of Asdex Upgrade
c     grids and those with a similar structure (Sonnet grid format)
c     
c     integer cutring,cutpt1,cutpt2,maxrings,maxkpts
c     data cutring /7/
c     data cutpt1 /25/
c     data cutpt2 /74/
c     data maxrings /26/
c     data maxkpts /98/
c     
      integer gridunit
      parameter (gridunit=4)
c     
      integer indexcnt,indexnadj,indexiradj,indexikadj
c     
      integer ik,ir,in,id,ikold,irold,loop_cnt
      integer max_ikold,irtmp
      integer ikn,irn
c     
      integer in1,in2,in3
c     
      character*100 buffer
      real rvert(4),zvert(4)
      real rcent,zcent
      real rshift,zshift
      real brat
c
c     jdemod - Add factors to scale grid if desired
c
      real rscale_grid,zscale_grid
c slmod begin
      real b_scale
c slmod end
c     
c     double precision rvert(4),zvert(4)
c     double precision rcent,zcent
c     double precision rshift,zshift
c     double precision brat
c     
c     Count of number of PSI data points
c     
c     
      integer npsi,ios
      real psi
c     
c     Add indicator for extra boundary cells offset into the grid.
c     
      integer ix_cell_offset  
c     
c     jdemod - nopriv moved to common cgeom
c     
c     logical nopriv
c     
c     Initialization
c     
      ios = 0
      indexcnt = 0
      indexnadj = 0
      indexiradj = 0
      indexikadj = 0
c
      nves = 0 
      npsi = 0
c
c     Init grid shift factors
c
      rshift   = 0.0
      zshift   = 0.0
c
c     jdemod - Init grid scaling factors
c
      rscale_grid = 1.0
      zscale_grid = 1.0
c     
      max_ikold = 1
      ikold = 1
      irold = 1
      call izero(nks,maxnrs)
c     
      ix_cell_offset = 0
c     
      crun = ' SONNET GRID (AUG,CMOD,DIIID,TDV...)'
      tslice = 0.0
c     slmod begin - tr
c...  Check if it is a quasi-double-null grid:
      READ(gridunit,'(A100)') buffer
      IF (sloutput) WRITE(0,*) 'BUFFER:'//buffer(1:20)//':'
      IF     (buffer(1:17).EQ.'QUASI-DOUBLE-NULL') THEN ! A couple of DIII-D grid still using this...
         IF (sloutput) WRITE(0,*) 'CALLING ReadQuasiDoubleNull'
         CALL ReadQuasiDoubleNull(gridunit,ik,ir,rshift,zshift,
     .        indexiradj)
         GOTO 300
      ELSEIF (buffer(1:19).EQ.'GENERALISED_GRID_SL') THEN
         WRITE(0,*) 'CALLING ReadGeneralisedGrid_SL'
        CALL ReadGeneralisedGrid_SL(gridunit,ik,ir,rshift,zshift,
     .                              indexiradj)
        GOTO 300
      ELSEIF (buffer(1:20).EQ.'GENERALISED_GRID_OSM') THEN
        IF (sloutput) WRITE(0,*) 'CALLING ReadGeneralisedGrid_OSM'
        CALL ReadGeneralisedGrid_OSM(gridunit,ik,ir,rshift,zshift,
     .                               indexiradj)
        GOTO 300
      ELSEIF (buffer(1:16).EQ.'GENERALISED_GRID') THEN
        IF (sloutput) WRITE(0,*) 'CALLING ReadGeneralisedGrid'
        CALL ReadGeneralisedGrid(gridunit,ik,ir,rshift,zshift,
     .       indexiradj)
        GOTO 300
      ELSE
         IF (sloutput) WRITE(0,*) 'Standard RAUG grid load'
         BACKSPACE(gridunit)
      ENDIF
c     slmod end
c     
c     Read in lines until '======'
c     
c     write (6,*) 'Working ... '
c     
 100  read(gridunit,'(a100)',end=1000) buffer
c     
c     Scan through any header information - the grid comes after the
c     row of '======= ... '
c     
c     If the characteristic grid numbers are in the grid header - extract
c     them for later use.
c     
      if (buffer(1:5).eq.'GEOM:'.or.
     >     buffer(1:5).eq.'Geom:'.or.
     >     buffer(1:5).eq.'geom:') then
c     
c     The agreed upon order for the GEOM: line if it is present
c     is the following ...
c     
c     GEOM:  #rings  #cutring    #knots   #cutpoint1  #cutpoint2
c     
         read (buffer(6:),*) maxrings,cutring,maxkpts,cutpt1,cutpt2
c     
c     Check for existence of private plasma for grid - this process will
c     INCLUDE the boundary cells in the core.
c     
c     Setting these values to 1 and maxkpts will also result in a grid
c     without a PFZ but in this case the boundary cells will be stripped 
c     from the core rings.
c     
         if ((cutpt1.eq.0.or.cutpt1.eq.1)) then 
            nopriv = .true.
         else
            nopriv = .false.
         endif

      endif
c     
c     
c     
      if (buffer(1:12).eq.'CELL_OFFSET:'.or.
     >     buffer(1:12).eq.'Cell_Offset:'.or.
     >     buffer(1:12).eq.'cell_offset:') then
c     
c     If this line is present it specifies the IX offset into the 
c     grid where the target virtual cells can be found. The boundary
c     cells at the edge of the grid are stripped and the entire 
c     grid shifted so that the embedded boundary cells are moved to 
c     the ends. This process is repeated with the background plasma
c     data before the call to maptodiv.  
c     
         read (buffer(13:),*) ix_cell_offset
c     
c     write(0,*) 'IX_CELL_OFFSET:',ix_cell_offset
c     
      endif
c     
c     Check also for grid shift/displacement information
c     
      if (buffer(1:6).eq.'SHIFT:'.or.
     >     buffer(1:6).eq.'Shift:'.or.
     >     buffer(1:6).eq.'shift:') then
c     
c     The data is RSHIFT ZSHIFT
c     
         read (buffer(7:),*) rshift,zshift
c     
      endif
c     
c     
c     jdemod - Check for grid scaling information
c            - these values are used to scale up the cell and wall dimensions
c            - the default values are 1.0
c     
      if (buffer(1:6).eq.'SCALE:'.or.
     >     buffer(1:6).eq.'Scale:'.or.
     >     buffer(1:6).eq.'scale:') then
c     
c     The data is RSCALE_GRID, ZSCALE_GRID
c     
         read (buffer(7:),*) rscale_grid,zscale_grid
c     
      endif
c slmod begin
c...  Some grids required rescaling of the magnetic field ratio:
      if (buffer(1:8).eq.'B-SCALE:'.or.
     >    buffer(1:8).eq.'B-Scale:'.or.
     >    buffer(1:8).eq.'B-scale:') then
         read (buffer(7:),*) b_scale 
      endif   
c slmod end
c     
      if (buffer(4:8).ne.'=====') goto 100
c     
      write(6,*) 'GEOM:',maxrings,cutring,maxkpts,cutpt1,cutpt2
c     
c     write(0,*) 'GEOM:',maxrings,cutring,maxkpts,cutpt1,cutpt2
c     
c     Start reading in elements - three lines - different formats
c     And a fourth separator line
c     
c     The orientation of the polygons in the SONNET grids must match
c     that for JET grids - the original read statements assigned the
c     corners off by 90 degrees.
c     
c     
c     200  read(12,9000,end=300)in,ik,ir,rvert(1),zvert(1),rvert(2),zvert(2)
c     read(12,9001,end=2000) brat,rcent,zcent
c     read(12,9002,end=2000) rvert(4),zvert(4),rvert(3),zvert(3)
c     read(12,'(a)',end=2000) buffer
c     
c     Results 180 degrees out of phase
c     
c     200  read(12,9000,end=300)in,ik,ir,rvert(4),zvert(4),rvert(1),zvert(1)
c     read(12,9001,end=2000) brat,rcent,zcent
c     read(12,9002,end=2000) rvert(3),zvert(3),rvert(2),zvert(2)
c     read(12,'(a)',end=2000) buffer
c     
c     200  read(12,9000,end=300)in,ik,ir,rvert(1),zvert(1),rvert(4),zvert(4)
c     read(12,9001,end=2000) brat,rcent,zcent
c     read(12,9002,end=2000) rvert(2),zvert(2),rvert(3),zvert(3)
c     read(12,'(a)',end=2000) buffer
c     
c     200  read(12,9000,end=300)in,ik,ir,rvert(4),zvert(4),rvert(3),zvert(3)
c     read(12,9001,end=2000) brat,rcent,zcent
c     read(12,9002,end=2000) rvert(1),zvert(1),rvert(2),zvert(2)
c     read(12,'(a)',end=2000) buffer
c     
 200  read(gridunit,9000,end=300) in,ik,ir,
     >     rvert(2),zvert(2),rvert(3),zvert(3)
c     
c     write(6,'(a,3i5)') 'READ  :',ik,ir,in
c     
      read(gridunit,9001,end=2000) brat,rcent,zcent
      read(gridunit,9002,end=2000) rvert(1),zvert(1),rvert(4),zvert(4)
      read(gridunit,'(a)',end=2000) buffer
c     
c     Adjust all values for any grid shifts or scalings 
c     
      rcent = rscale_grid * (rcent + rshift)
      zcent = zscale_grid * (zcent + zshift)
c     
      do loop_cnt = 1,4
c     
         rvert(loop_cnt) = rscale_grid * (rvert(loop_cnt) + rshift)
         zvert(loop_cnt) = zscale_grid * (zvert(loop_cnt) + zshift)
c     
      end do    
c     
c     
c     write(15,9000) in,ik,ir,rvert(1),zvert(1),rvert(2),zvert(2)
c     write(15,9001) brat,rcent,zcent
c     write(15,9002) rvert(4),zvert(4),rvert(3),zvert(3)
c     write(15,'(a)') buffer
c     
c     Assign all the values to the appropriate arrays
c     Korpg is a cross-reference to the cell vertices. nvertp
c     contains the number of vertices for a cell ... which is
c     fixed at 4 for the present.
c     
c     Increment all of the indices by an adjustment because DIVIMP does not
c     use zeroth element array addressing for these arrays.
c     
      indexcnt = indexcnt + 1
c     
      if (indexcnt.eq.1) then
         indexnadj  = 1 - in 
         indexiradj = 1 - ir 
         indexikadj = 1 - ik 
c     
c     if (in.eq.-1) then
c     indexadj = 2
c     elseif (in.eq.0) then
c     indexadj = 1
c     elseif (in.eq.1) then
c     indexadj = 0
c     endif
      endif
c     
      ik = ik + indexikadj
      ir = ir + indexiradj
      in = in + indexnadj
c     
c     If IX_CELL_OFFSET is non-zero then the normal 
c     boundary cells have to be ignored/removed before the 
c     grid is shifted to move the second set of boundary cells 
c     to the targets.    
c     
c     write(6,'(a,4i5)') 'BEFORE:',ik,ir,in,ix_cell_offset
c     
      if (ix_cell_offset.ne.0.and.(ik.eq.1.or.ik.eq.maxkpts)) goto 150
c     
c     Find corrected value of IK for the cell - keeping in mind that the 2
c     end boundary cells have been removed. 
c     
      if (ix_cell_offset.ne.0) then 
c     
         if (ik.le.(ix_cell_offset+1) ) then 
            ik = ik + (maxkpts-(ix_cell_offset+1)-1) -1 
         else
            ik = ik - (ix_cell_offset+1)
         endif     
c     
      endif
c     
c     write(6,'(a,3i5)') 'MID   :',ik,ir,in 
c     
c     
c     Less than or equal to the cutring is the PFZ/core 
c     IK,IR values for these need further adjustments.
c     
      if (ir.le.cutring) then
         if (ik.le.cutpt1) then
            ir = maxrings + ir
         elseif (ik.ge.cutpt2) then
c     
c     Moved code to add repeat cells in core to after 
c     entire grid is loaded - this allows a split 
c     grid to be properly handled. 
c     
c     
            ir = maxrings + ir
            ik = ik - cutpt2 + cutpt1 +1
c     
         else
            ik = ik - cutpt1
         endif
      endif
c     
c     write(6,'(a,3i5)') 'AFTER :',ik,ir,in 
c     write(6,*)  
c     
      korpg(ik,ir) = in
      nvertp(in) = 4
      do 210 id = 1,nvertp(in)
         rvertp(id,in) = rvert(id)
         zvertp(id,in) = zvert(id)
 210  continue
      rs(ik,ir) = rcent
      zs(ik,ir) = zcent
      bratio(ik,ir) = brat
c     
c     Check to see if a ring end has been passed.
c     
c     For the trap rings ... this will be filled out twice
c     first when it hits the cut and again at the end of the
c     ring ... i.e. the two times when ir changes.
c     
c     
      if (ir.ne.irold) then 
         nks(irold) = max(nks(irold),max_ikold)
         max_ikold = 1  
      endif 
c     
      ikold = ik
      irold = ir
      max_ikold = max(max_ikold,ikold)
c     
c     Exit condition ... can be supplemented with the END of file
c     exit in the first read statement at line 200.
c     

 150  continue

      if (ir.eq.maxrings.and.ik.eq.maxkpts) then
         nks(ir) = max(nks(ir),max_ikold)
         npolyp = in
         

c     
c     Add boundary cells for sonnet grids without virtual cells - 
c     Add rings and add cells at the end of rings 
c     Add cells by duplicating the cells from the adjacent ring with zero volume
c     
         if (sonnet_grid_sub_type.eq.2) then 
c     
c     
c     If needed add boundary cells/rings to the grid - this is for grids without virtual cells around the edges
c     
c     Variables affected - 
c     nvertp(ip)
c     rvertp(1:4,ip)
c     zvertp(1:4,ip)
c     
c     korpg(ik,ir)
c     bratio(ik,ir)
c     rs(ik,ir)
c     zs(ik,ir)
c     nrs
c     nks(ir)
c     
c     maxrings (total rings read in - sonnet format)
c     cutring  (last ring in core)
c     cutpt1,cutpt2
c     
c     Add virtual rings first 
c     Then add virtual target points at the ends of each SOL/PFZ ring
c     
c     Procedure only works with SN at the moment
c     irsep = cutring + 1
c     irwall = maxrings
c     nrs = maxrings + cutring
c     
c     
            if (maxrings+cutring+3 .gt.maxnrs) then 
               write(0,*) 'ERROR ADDING BOUNDARY RINGS:'//
     >              ' MAXNRS TOO SMALL = ',maxnrs
               write(6,*) 'ERROR ADDING BOUNDARY RINGS:'//
     >              ' MAXNRS TOO SMALL = ',maxnrs
               stop
            endif
c
c     Create space at ir = maxrings+1 for SOL and PFZ boundary
c     Add space for 2 rings in the array
c
            do irn = maxrings+cutring,maxrings+1,-1
               do ikn = 1,nks(irn)
                  korpg(ikn,irn+2) = korpg(ikn,irn)
                  rs(ikn,irn+2) = rs(ikn,irn)
                  zs(ikn,irn+2) = zs(ikn,irn)
                  bratio(ikn,irn+2) = bratio(ikn,irn)
               end do 
               nks(irn+2) = nks(irn)
            end do
c     
c     Create space at ir = 1 for PFZ boundary
c     
            do irn = maxrings+cutring+2,1,-1
               do ikn = 1,nks(irn)
                  korpg(ikn,irn+1) = korpg(ikn,irn)
                  rs(ikn,irn+1) = rs(ikn,irn)
                  zs(ikn,irn+1) = zs(ikn,irn)
                  bratio(ikn,irn+1) = bratio(ikn,irn)
               end do 
               nks(irn+1) = nks(irn)
            end do
c     
c     Fill in the core boundary
c     
            call add_boundary_ring(in,1,2,INWARD41)
c     
c     Fill in the PFZ boundary
c     
            call add_boundary_ring(in,maxrings+3,maxrings+4,INWARD41)
c     
c     Fill in SOL boundary 
c     
            call add_boundary_ring(in,maxrings+2,maxrings+1,OUTWARD23)
c     
c     Adjust grid ring parameters
c     
            maxrings = maxrings+2
            ir = ir + 2
            cutring = cutring+1
c     
c     Add poloidal boundary cells at the ends of all rings in SOL/PFZ
c     
            do irn = cutring+1,maxrings+cutring
               call add_boundary_targets(in,irn)
            end do
c
c     Adjust cell counts appropriately
c
            ik = ik + 2
            maxkpts  = maxkpts + 2
            cutpt1 = cutpt1 + 1
            cutpt2 = cutpt2 + 1
c
         endif
c
         goto 300
      endif
c     
c     Loop back
c     
      goto 200
c     
c     
c     Continuation point after grid has been read in
c     
c     
 300  continue
c     slmod begin

c      CALL OutputData(86,'300 of RAUG')    
c         STOP 'RAUG MID'

      CALL DB('RAUG: Done reading grid')

      vpolmin = (MAXNKS*MAXNRS - npolyp) / 2 + npolyp
      vpolyp  = vpolmin
c     slmod end
c     
c     Check for complete grid and issue error message
c     
      if (ir.ne.maxrings.or.ik.ne.maxkpts) then
c     
         write (6,*) 'Error reading SONNET geometry file!'
         write (6,*) 'Short read in data - grid may be incomplete'
         write (6,*) 'Last read: ring:',ir,' knot:',ik
         write (6,*) 'Number of polygons:',in
         write (6,*) 'GEOM:',maxrings,cutring,maxkpts,cutpt1,cutpt2
c     
         write (0,*) 'Error reading SONNET geometry file!'
         write (0,*) 'Short read in data - grid may be incomplete'
         write (0,*) 'Last read: ring:',ir,' knot:',ik
         write (0,*) 'Number of polygons:',in
         write (0,*) 'GEOM:',maxrings,cutring,maxkpts,cutpt1,cutpt2
c     
         npolyp = in
c     
         stop
c     
      endif
c     
c     Add repeat cells to core rings - those less than or equal to cutring
c     
c     For grids with NO PFZ - this needs to be modified since the "boundary" cells
c     still have a finite but vey small volume - if these are deleted from the core
c     rings (as they have been) then the core will not match up with the main SOL since
c     the boundary cells are removed in the main SOL by code found later. As a result, 
c     we need to manufacture a new cell by combining the two guard cells with one 
c     of the adjacent cells and then duplicating that result.
c     
c     This procedure applies for LIMITER machine grids generated by UEDGE. A different
c     procedure would be necessary for cells with zero volumes.  
c     
      if (sonnet_grid_sub_type.eq.0.or.sonnet_grid_sub_type.eq.2) then

         do ir = 1,cutring 
c     
            if (nopriv) then 
c     
c     Modify the along field line coordinates of the last cell of the core
c     ring so that they will match up with the first cell.  
c     
               in1 = korpg(1,ir) 
               in2 = korpg(nks(ir),ir) 
c     
               rvertp(4,in2) = rvertp(1,in1)
               zvertp(4,in2) = zvertp(1,in1)
               rvertp(3,in2) = rvertp(2,in1)
               zvertp(3,in2) = zvertp(2,in1)
c     
c     Also need to modify one corner of the first non-boundary main SOL cell. 
c     so that the target closes and the corners match up correctly. 
c     
c     NOTE: if the cells inserted actually had zero volume none of this 
c     additional processing would be needed. 
c     
c     This modification should not be done for the FRC grid.  
c     
               if (ir.eq.cutring) then 
c     
                  in3 = korpg(nks(ir+1)-1,ir+1)
                  rvertp(4,in3) = rvertp(3,in2)
                  zvertp(4,in3) = zvertp(3,in2)
c     
               endif 
c     
            endif
c     
            nks(ir) = nks(ir) + 1
c     
            ik = nks(ir)           
c     
c     Add a repeat point to the core
c     plasma rings - to close them.
c     
            rs(ik,ir) = rs(1,ir)
            zs(ik,ir) = zs(1,ir)
            korpg(ik,ir) = korpg(1,ir)
            bratio(ik,ir) = bratio(1,ir)
c     
         end do
c     
      endif 


c     
c     
c     READ in any additional information at the end of the grid file
c     
c     1) Target PSI values
c     2) NEUTRAL WALL coordinates (if any)
c     
      
 400  read(gridunit,'(a)',end=500) buffer

c     
c     Check for neutral wall 
c     
      if (buffer(1:13).eq.'NEUTRAL WALL:') then 

         read(buffer(14:),*,IOSTAT=ios) nves
c     
c     Support old format where number of elements was on next line
c     Trapped by error in reading the new format
c     
         if (ios.ne.0) then 
            read(gridunit,*) nves
         endif
c     
         do in = 1,nves
            read(gridunit,*,end=4000) rves(in),zves(in)
c     
c     Need to apply same shift and scale to neutral wall as to grid
c     
            rves(in) = rscale_grid * (rves(in) + rshift)
            zves(in) = zscale_grid * (zves(in) + zshift)
c     
         end do


      endif  
      IF (.TRUE.) THEN
         nvesm = nves
         DO i1 = 1, nves-1
            rvesm(i1,1) = rves(i1)
            zvesm(i1,1) = zves(i1)
            rvesm(i1,2) = rves(i1+1)
            zvesm(i1,2) = zves(i1+1)
         ENDDO
         rvesm(i1,1) = rves(i1)
         zvesm(i1,1) = zves(i1)
         rvesm(i1,2) = rves(1)
         zvesm(i1,2) = zves(1)
      ENDIF
c     
c     Check for PSI values
c     
      if (buffer(1:4).eq.'PSI:') then
c     
c        Initialize psitarg array to zero
c
         psitarg = 0.0
c
         read(buffer(5:),*) npsi
c     
c     The PSI values are listed with one on each line
c     indexed by knot and ring index based on the SONNET 
c     grid coordinates. 
c     
         do in = 1,npsi
c     
            read(gridunit,*,end=4000) ik,ir,psi
c
c     In the case of sonnet_grid_sub_type=2 - boundary rings have been added
c     and these will not have corresponding psi values - set psi values to 0.0?
c
c     
c     IK values equal to 1 correspond to the INSIDE target
c     for an X-point down DIIID grid - this would correspond
c     to target number "2" in the DIVIMP nomenclature
c     
c
c     jdemod - the table of PSI values at the end of the grid file has
c              never been indexed from 0 - it has always started at 1 for the 
c              first ring - however, this is due to the fact that UEDGE 
c              appears to have just skipped printing values for boundary 
c              rings - the correction is still needed to map properly to 
c              the grid. 
c
            ir = ir + indexiradj
c
c
c     In the case of sonnet_grid_sub_type=2 - boundary rings have been added
c     and these will not have corresponding psi values - set psi values to 0.0 as default
c     Adding one to ir should account for the indexing. 
c
            if (sonnet_grid_sub_type.eq.2) then
               ir = ir +1 
               write(6,'(a,4i8,g12.5)') 'PSITARG:',ir,indexiradj,
     >                        cutring,ik,psi
            endif
c     
            if (ir.le.cutring) ir = ir + maxrings 
c     
c     Check for first or second half of ring
c     
            if (ik.lt.maxkpts/2) then 
c     
               psitarg(ir,2) = psi
c     
            else
c     
               psitarg(ir,1) = psi
c     
            endif
c     
         end do
c     
      endif
c     

      goto 400


 500  continue
c     
c     
c     Check to see if vessel wall was READ IN - If required
c     
      if ((cneur.eq.4.or.ctrap.eq.4.or.cneur.eq.5.or.ctrap.eq.5)
     >     .and.nves.eq.0) then
c     
         write (6,*) 'Error reading SONNET geometry file!'
         write (6,*) 'The Neutral Wall data was not appended'
         write (6,*) '(but is requored for specified options.)'
         call prc('Error reading SONNET geometry file!')
         call prc('The Neutral Wall data was not appended')
         call prc('(but is requored for specified wall options)')
         stop
c     
      end if
c     
      write (6,'(a,9(1x,i6))') 'Index:',
     >     indexcnt,indexnadj,indexiradj,indexikadj,
     >     ik,ir,in,maxrings,
     >     maxkpts
c     
      do ir = 1,maxrings+cutring
         write (6,'(a,5i5)') 'RINGS:',ir,nks(ir),maxrings,cutring
      end do
c     
      do ir = 1,maxrings+cutring
         do ik = 1,nks(ir)
            write (55,'(2g18.10,2i5)') rs(ik,ir),zs(ik,ir), ik,ir
         end do
      end do
c     
c     
c     Set total number of rings and total number of polygons
c     
      refct = 0
c     
      if (nopriv) then
         nrs = maxrings
      else
         nrs = maxrings + cutring
      endif
c     
      irsep = cutring +1
      irwall = maxrings
c     
c     These may need adjusting for meshes without a private plasma
c     
      irtrap  = irwall +1
      irtrap2 = irtrap
      irwall2 = irwall
      nrs2 = nrs
c     
c     
c     For the FRC grid replace any boundary cells whose coordinates are all 
c     zeroes with cell vertex values taken from the adjacent real cell.
c     
c     Rings affected are IR=1 and IR=NRS/IRWALL and the first cells on each ring
c     
c     NOTE!: The polygons created here are in fact later deleted - however, the RS,ZS
c     values retain significance.
c     
      if (sonnet_grid_sub_type.eq.1) then 
c     
c     Ring 1 and Ring NRS first from knot 2 to knot nks(ir)-1
c     
         ir = 1
c     
c     Map outward for first ring - do all cells on ring since extra cells on core
c     rings are not later stripped while they are for the other rings. 
c     The "virtual" cells at the ends of this ring have already been removed to 
c     the PFZ where they are deleted. 
c     
c     In addition the core ring needs to be lined up with the ik+1 cell in the 
c     adjacent main SOL ring since they is a discontinuity in cell alignment. 
c     
         do ik = 1,nks(ir)
c     

            if ((rs(ik,ir).eq.rshift) .and.
     >           (zs(ik,ir).eq.zshift)) then 
c     
               in1 = korpg(ik,ir)
c     
c     Line up cell with appropriate one in next ring
c     
               if ((ir+1).eq.irsep) then 
                  in2 = korpg(ik+1,ir+1)
               else
                  in2 = korpg(ik,ir+1)
               endif
c     
c     Copy cell vertices from next ring
c     
               rvertp(3,in1) = rvertp(4,in2)
               zvertp(3,in1) = zvertp(4,in2)
               rvertp(2,in1) = rvertp(1,in2)
               zvertp(2,in1) = zvertp(1,in2)
c     
c     Duplicate these for other vertices to get a 
c     zero volume cell
c     
               rvertp(4,in1) = rvertp(3,in1)
               zvertp(4,in1) = zvertp(3,in1)
               rvertp(1,in1) = rvertp(2,in1)
               zvertp(1,in1) = zvertp(2,in1)
c     
               rs(ik,ir) = (rvertp(3,in1) + rvertp(2,in1))/2.0
               zs(ik,ir) = (zvertp(3,in1) + zvertp(2,in1))/2.0
               bratio(ik,ir) = bratio(ik,ir+1)
c     
            endif
c     
         end do 
c     
         ir = nrs
c     
c     Map inward for last ring
c     
         do ik = 2,nks(ir)-1

c     
            if ((rs(ik,ir).eq.rshift) .and.
     >           (zs(ik,ir).eq.zshift)) then 
c     
               in1 = korpg(ik,ir)
               in2 = korpg(ik,ir-1)
c     
c     Copy cell vertices from next ring
c     
               rvertp(4,in1) = rvertp(3,in2)
               zvertp(4,in1) = zvertp(3,in2)
               rvertp(1,in1) = rvertp(2,in2)
               zvertp(1,in1) = zvertp(2,in2)
c     
c     Duplicate these for other vertices to get a 
c     zero volume cell
c     
               rvertp(3,in1) = rvertp(4,in1)
               zvertp(3,in1) = zvertp(4,in1)
               rvertp(2,in1) = rvertp(1,in1)
               zvertp(2,in1) = zvertp(1,in1)
c     
               rs(ik,ir) = (rvertp(4,in1) + rvertp(1,in1))/2.0
               zs(ik,ir) = (zvertp(4,in1) + zvertp(1,in1))/2.0
               bratio(ik,ir) = bratio(ik,ir-1)
c     
            endif
c     
         end do 
c     
c     Assign vertices for the cells at the
c     ends of the rings for cells in the main SOL. 
c     
         do ir = irsep,nrs
c     
c     First end - ik = 1 
c     
            ik = 1
c     
            if ((rs(ik,ir).eq.rshift) .and.
     >           (zs(ik,ir).eq.zshift)) then 
c     
               in1 = korpg(ik,ir)
               in2 = korpg(ik+1,ir)
c     
c     Copy cell vertices from next cell
c     
               rvertp(4,in1) = rvertp(1,in2)
               zvertp(4,in1) = zvertp(1,in2)
               rvertp(3,in1) = rvertp(2,in2)
               zvertp(3,in1) = zvertp(2,in2)
c     
c     Duplicate these for other vertices to get a 
c     zero volume cell
c     
               rvertp(1,in1) = rvertp(4,in1)
               zvertp(1,in1) = zvertp(4,in1)
               rvertp(2,in1) = rvertp(3,in1)
               zvertp(2,in1) = zvertp(3,in1)
c     
               rs(ik,ir) = (rvertp(4,in1) + rvertp(3,in1))/2.0
               zs(ik,ir) = (zvertp(4,in1) + zvertp(3,in1))/2.0
               bratio(ik,ir) = bratio(ik+1,ir)
c     
            endif
c     
c     Second end - ik = nks(ir)
c     
            ik = nks(ir)
c     
            if ((rs(ik,ir).eq.rshift) .and.
     >           (zs(ik,ir).eq.zshift)) then 
c     
               in1 = korpg(ik,ir)
               in2 = korpg(ik-1,ir)
c     
c     Copy cell vertices from next cell
c     
               rvertp(1,in1) = rvertp(4,in2)
               zvertp(1,in1) = zvertp(4,in2)
               rvertp(2,in1) = rvertp(3,in2)
               zvertp(2,in1) = zvertp(3,in2)
c     
c     Duplicate these for other vertices to get a 
c     zero volume cell
c     
               rvertp(4,in1) = rvertp(1,in1)
               zvertp(4,in1) = zvertp(1,in1)
               rvertp(3,in1) = rvertp(2,in1)
               zvertp(3,in1) = zvertp(2,in1)
c     
               rs(ik,ir) = (rvertp(1,in1) + rvertp(2,in1))/2.0
               zs(ik,ir) = (zvertp(1,in1) + zvertp(2,in1))/2.0
               bratio(ik,ir) = bratio(ik-1,ir)
c     
            endif  
c     
         end do 
c     
      endif
c     
c     Approximate X-point values and centre point
c     
c     Use approximate limiter tip for Textor style limiter
c     grids - or those for which there is no private plasma
c     specified at all.
c     
      if (nopriv.or.(cutpt1.eq.1.and.cutpt2.eq.maxkpts)) then
         rxp = (rs(1,irsep)+ rs(nks(irsep),irsep)) / 2.0
         zxp = (zs(1,irsep)+ zs(nks(irsep),irsep)) / 2.0

      else
         rxp = 0.25 * (rs(cutpt1,nrs)+rs(cutpt1+1,nrs)
     >        +rs(1,cutring)+rs(nks(cutring)-1,cutring))
         zxp = 0.25 * (zs(cutpt1,nrs)+zs(cutpt1+1,nrs)
     >        +zs(1,cutring)+zs(nks(cutring)-1,cutring))
      endif
c     
c     Need to verify that the innermost core ring on the grid actually contains valid r,s data
c     especially when r0 and z0 are used elsewhere - use the r,z coordinates for the second 
c     ring of the core since it should be guaranted to contain real data.
c     
      r0 = 0.0
      z0 = 0.0
      irtmp = 2
      do 70 ik = 1,nks(irtmp)-1
         r0 = r0 + rs(ik,irtmp)
         z0 = z0 + zs(ik,irtmp)
 70   continue
      r0 = r0/(nks(irtmp)-1)
      z0 = z0/(nks(irtmp)-1)
c     
      write (6,*) 'Calculated: rxp,zxp,r0,z0 = ',rxp,zxp,r0,z0
c     
c     
c     Set up value of KBFS and BTS
c     
c     Also set up values of RMIN,RMAX,ZMIN,ZMAX
c     
      rmin = hi
      rmax = lo
      zmin = hi
      zmax = lo
c     
      do 50 ir = 1,nrs
         do 60 ik = 1,nks(ir)
            if (bratio(ik,ir).ne.0.0) then
               kbfs(ik,ir) = 1.0 / bratio(ik,ir)
            else
               write(6,'(a,2i6,4(1x,g12.5))') 
     >              'WARNING: RAUG: BRATIO=0:',ik,ir,bratio(ik,ir)
               kbfs(ik,ir) = 1.0
            endif
c     
c     Assign values for toroidal magnetic field
c     
            if (rs(ik,ir).ne.0.0) then
               BTS(IK,IR) = CBPHI * R0 / RS(IK,IR)
            else 
               write(6,'(a,2i6,2g14.6)') 
     >              'WARNING: RAUG: R=0 FOR CELL:',
     >              ik,ir,rs(ik,ir),zs(ik,ir)
               bts(ik,ir) = cbphi
            endif
c     
            RMIN = MIN(rmin,rs(ik,ir))
            rmax = max(rmax,rs(ik,ir))
            zMIN = MIN(zmin,zs(ik,ir))
            zmax = max(zmax,zs(ik,ir))
 60      continue
 50   continue
c     
      if (cprint.eq.3.or.cprint.eq.9) then 
         write (6,*) 'KBFS:'
c     
         do ir = 1,nrs
c     
            write (6,'(''Ring:'',2i4)') ir,nks(ir)
            do ik = 1,nks(ir)
               write (6,'(2g18.10,i5,2g18.10)') bratio(ik,ir),
     >              kbfs(ik,ir),
     >              korpg(ik,ir),rs(ik,ir),zs(ik,ir)
            end do
         end do
      endif   
c     
c     
c     Other values
c     
c     These settings for IKT will only work for a symmetric grid
c     where the number of points on the private plasma sections at
c     each end are the same.
c     

      IKT = cutpt1
      IKTO= IKT
c     
c     IKTI = nks(irsep) - ikt + 1
c     
c     This formula should work for assymmetric legs
c     
      if (ix_cell_offset.ne.0) then   
         ikti = nks(irsep) - ( (maxkpts-2) - cutpt2 ) 
      else
         ikti = nks(irsep) - ( maxkpts - cutpt2 )
      endif
c     
      write (6,*) 'Cut values: ikt...:',ikt,ikto,ikti,
     >     cutpt1,cutpt2,maxkpts,
     >     cutring, maxrings, irsep,nrs
c     
c     NINOWA = 0
c     HMASS = 0.0
c     ZMASS = 0.0
c     
c     
c     Calculate KORY for non-JET grids so that the TRAN file
c     can be written correctly. Note: at this point the 
c     boundary cells are still on the grid so no additional 
c     work has to be done. 
c     
c     KORY is used to index and store data in the appropriate order
c     for use at JET. 
c     
      in = 0
      do ir = 1,nrs
         do ik = 1,nks(ir)
            in = in+1
            kory(ir,ik) = in 
         end do
      end do
c     
c     If it has been specified to read the background plasma
c     solution from an external file. It is assumed that the
c     values are the result of a B2-Eirene run and the files
c     will be in the format of the B2 to Eirene coupling. The
c     code to read this was taken from this interface and should
c     be identical. Furthermore, care must be taken to ensure that
c     the grids and numbers of values match up. The code makes
c     the same assumptions as iterated at the beginning of the
c     RAUG subroutine ... in regards to mapping the plasma values
c     onto the DIVIMP grid.
c     
C     
C-----------------------------------------------------------------------
C     Extract Densities, Temperatures etc from second file
C-----------------------------------------------------------------------
C     
      if (cre2d.eq.2.or.cioptf.eq.99.or.
     >     cioptg.eq.99.or.cioptg.eq.90) then

         call b2repl(maxrings,maxkpts,cutring,cutpt1,cutpt2,readaux,
     >        rizb,crmb,cion,ix_cell_offset)

      endif
c     
c     Write diagnostics
c     
      write (6,'(a,4(1x,i8))') 
     >     'nrs, irsep, irwall, irtrap as initially read:',
     >     nrs,irsep,irwall,irtrap
      write (6,*) 'nks as initially read:',(nks(ir),ir=1,nrs)
      write (6,'(a,4(1x,i8))') 
     >     'maxrings, maxkpts, cutring, cutpt1, cutpt2:',
     >     maxrings, maxkpts, cutring, cutpt1, cutpt2
c     
      do 10 ir = 1,cutring
         do 20 ik = 1,nks(ir)
            write (diagunit,'(2i5,3(1x,g18.12))') 
     >           ik,ir,rs(ik,ir),zs(ik,ir),kbfs(ik,ir)
 20      continue
         write(diagunit,'(a)')
 10   continue
      do 30 ir = cutring+1,nrs
         do 40 ik = 1,nks(ir)
            write (diagunit,'(2i5,3(1x,g18.12))') 
     >           ik,ir,rs(ik,ir),zs(ik,ir),kbfs(ik,ir)
 40      continue
         write(diagunit,'(a)')
 30   continue
c     
c     jdemod - Output the grid before modifications are made
c     
      call OutputGrid2(67,'RAUG: before modifications')

c     slmod begin 
      IF (quasidn) CALL PrepQuasiDoubleNull

c...  Tailor/cut grid to wall:
      IF (grdnmod.GT.0) THEN
c...    Get rid of poloidal boundary cells (to be added again below
c       after grid manipulations are complete):
        DO ir = irsep, nrs
          CALL DeleteCell(nks(ir),ir)
          CALL DeleteCell(1      ,ir)
        ENDDO
c...    Modify the grid based on entries in the GRDMOD array assigned 
c       from the input file:
        CALL TailorGrid
c...    Add virtual boundary cells, which will be stripped off later:
        IF (CTARGOPT.EQ.0.OR.CTARGOPT.EQ.1.OR.CTARGOPT.EQ.2.OR.
     .      CTARGOPT.EQ.3.OR.CTARGOPT.EQ.6) 
     .      CALL AddPoloidalBoundaryCells
      ENDIF      
c     
c     jdemod - write out grid after modifications
c     
      call OutputGrid2(68,'RAUG after grid modifications')

c     WRITE(0,*) 'GRDNMOD set = -1'
c     grdnmod = -1 
c     slmod end
c     
c     
c     Set IKREF value
c     
      IKREF  = NKS(IRSEP)/2 + 1
c     
c     REMOVE THIS FOR NOW - Having no polygon or cell information for
c     the boundary rings can cause some issues with central mirror reflection
c     and far periphery entry - needs to be fixed up. 
c     
c     Eliminate any record of polygons for the core, wall and
c     trap wall rings since these are not proper rings in the
c     first place.
c     
c     if (sonnet_grid_sub_type.eq.0) then  
c     
c     write (6,*) 'Eliminating KORPG for :', 1,irwall,irtrap
c     
c     CALL ISet(korpg(1,1     ),MAXNKS,0)
c     CALL ISet(korpg(1,irwall),MAXNKS,0)
c     CALL ISet(korpg(1,irtrap),MAXNKS,0)
c     
c     endif 
c     
c     Core
c     
c     ir = 1
c     do ik = 1,nks(ir)
c     korpg(ik,ir) = 0
c     end do
c     
c     Wall
c     
c     ir = irwall
c     do ik = 1,nks(ir)
c     korpg(ik,ir) = 0
c     end do
c     
c     Trap Wall
c     
c     ir = irtrap
c     do ik = 1,nks(ir)
c     korpg(ik,ir) = 0
c     end do
c     
      return
c     
c     Error exit conditions
c     
 1000 write (6,*) 'Error reading SONNET geometry file!'
      write (6,*) 'No starting line ===== found.'
      stop

 2000 write (6,*) 'Error reading SONNET geometry file!'
      write (6,*) 'Short read in data - missing data'
      stop

 4000 write (6,*) 'Error reading SONNET geometry file!'
      write (6,*) 'INCOMPLETE READ ON ADDITIONAL DATA'
      call prc('Error reading SONNET geometry file!')
      call prc('INCOMPLETE READ ON ADDITIONAL DATA')
      stop
c     
c     Format statements
c     
c     The brain dead SUN compiler will not accept constant character
c     string in an input format statement even if the match the input 
c     text - changing to use x edit descriptor to space data for input -
c     use original format statements for output. 
c     

 9000 format(10x,i5,4x,i3,1x,i3,4x,
     >     e17.10e2,1x,e17.10e2,8x,e17.10e2,1x,e17.10e2)
 9001 format(18x,e17.10e2,14x,
     >     e17.10e2,1x,e17.10e2)
 9002 format(30x,e17.10e2,1x,e17.10e2,8x,
     >     e17.10e2,1x,e17.10e2)
c     
c     Original formats for writing 
c     
 9003 format(3x,'Element',i5,' = (',i3,',',i3,'): (',
     >     e17.10e2,',',e17.10e2,')',6x,'(',e17.10e2,',',e17.10e2,')')
 9004 format(3x,'Field ratio  = ',e17.10e2,13x,
     >     '(',e17.10e2,',',e17.10e2,')')
 9005 format(3x,26x,'(',e17.10e2,',',e17.10e2,')',6x,
     >     '(',e17.10e2,',',e17.10e2,')')
c     
c     
      end
c     
c     
c     

      subroutine add_boundary_ring(in,irn,irref,ref_side)
      implicit none
      include 'params'
      include 'cgeom'
      
      integer, intent(in) ::  irn,irref,ref_side
      integer in 

c     
c     Local variables
c     
      integer ikn,inref

      nks(irn) = nks(irref)

      do ikn = 1,nks(irn)
         inref = korpg(ikn,irref)

         call add_boundary_cell(in,ikn,irn,inref,ref_side)

         bratio(ikn,irn) = bratio(ikn,irref)

      end do

c     
      return
      end
c     
c     
c     
      subroutine add_boundary_targets(in,irn)
      implicit none
      include 'params'
      include 'cgeom'
c
      integer in,irn
c
      integer ikn,inref

      if (nks(irn)+2 .gt.maxnks) then 
         write(0,*) 'ERROR ADDING BOUNDARY AT TARGETS:'//
     >              ' MAXNKS TOO SMALL = ',maxnks
         write(6,*) 'ERROR ADDING BOUNDARY AT TARGETS:'//
     >              ' MAXNKS TOO SMALL = ',maxnks
         stop
      endif

c     
c     Make space for new cells
c     
      do ikn = nks(irn),1,-1
         korpg(ikn+1,irn) = korpg(ikn,irn)
         rs(ikn+1,irn) = rs(ikn,irn)
         zs(ikn+1,irn) = zs(ikn,irn)
         bratio(ikn+1,irn) = bratio(ikn,irn)
      end do
c     
c     Add first cell
c     
      ikn=1
      inref=korpg(ikn+1,irn)
      bratio(ikn,irn) = bratio(ikn+1,irn)
c     
      call add_boundary_cell(in,ikn,irn,inref,DOWN12)
c     
c     Add last cell
c     
      ikn = nks(irn)+2
      inref=korpg(ikn-1,irn)
      bratio(ikn,irn) = bratio(ikn-1,irn)
c     
      call add_boundary_cell(in,ikn,irn,inref,UP34)
c     
      nks(irn) = nks(irn) + 2

      return
      end
c     
c     
c     
      subroutine add_boundary_cell(in,ikn,irn,inref,ref_side)
      implicit none
      include 'params'
      include 'cgeom'
      
      integer in,inref,ref_side,ikn,irn

      in = in + 1
c     
c     Add new cell data
c     
      korpg(ikn,irn) = in
      nvertp(in) = 4
c     
c     For outer boundary need to set vertices 2,3 to adjacent real cell values - then set
c     vertex 2 = 1 and 3 = 4
c     
c     
c     Inner edge of grid
c     
      if (ref_side.eq.INWARD41) then 

         rvertp(3,in)=rvertp(4,inref) 
         rvertp(2,in)=rvertp(1,inref) 
         rvertp(4,in)=rvertp(3,in) 
         rvertp(1,in)=rvertp(2,in) 

         zvertp(3,in)=zvertp(4,inref) 
         zvertp(2,in)=zvertp(1,inref) 
         zvertp(4,in)=zvertp(3,in) 
         zvertp(1,in)=zvertp(2,in) 

         rs(ikn,irn) = (rvertp(3,in)+rvertp(2,in))/2.0
         zs(ikn,irn) = (zvertp(3,in)+zvertp(2,in))/2.0

c     
c     Outer edge of grid
c     
      elseif (ref_side.eq.OUTWARD23) then 

         rvertp(4,in)=rvertp(3,inref) 
         rvertp(1,in)=rvertp(2,inref) 
         rvertp(2,in)=rvertp(1,in) 
         rvertp(3,in)=rvertp(4,in) 

         zvertp(4,in)=zvertp(3,inref) 
         zvertp(1,in)=zvertp(2,inref) 
         zvertp(2,in)=zvertp(1,in) 
         zvertp(3,in)=zvertp(4,in) 

         rs(ikn,irn) = (rvertp(4,in)+rvertp(1,in))/2.0
         zs(ikn,irn) = (zvertp(4,in)+zvertp(1,in))/2.0

c     
c     End of ring
c     
      elseif (ref_side.eq.UP34) then 

         rvertp(2,in)=rvertp(3,inref) 
         rvertp(1,in)=rvertp(4,inref) 
         rvertp(4,in)=rvertp(1,in) 
         rvertp(3,in)=rvertp(2,in) 

         zvertp(2,in)=zvertp(3,inref) 
         zvertp(1,in)=zvertp(4,inref) 
         zvertp(4,in)=zvertp(1,in) 
         zvertp(3,in)=zvertp(2,in) 

         rs(ikn,irn) = (rvertp(2,in)+rvertp(1,in))/2.0
         zs(ikn,irn) = (zvertp(2,in)+zvertp(1,in))/2.0
c     
c     Start of ring
c     
      elseif (ref_side.eq.DOWN12) then 

         rvertp(4,in)=rvertp(1,inref) 
         rvertp(3,in)=rvertp(2,inref) 
         rvertp(1,in)=rvertp(4,in) 
         rvertp(2,in)=rvertp(3,in) 

         zvertp(4,in)=zvertp(1,inref) 
         zvertp(3,in)=zvertp(2,inref) 
         zvertp(1,in)=zvertp(4,in) 
         zvertp(2,in)=zvertp(3,in) 

         rs(ikn,irn) = (rvertp(4,in)+rvertp(3,in))/2.0
         zs(ikn,irn) = (zvertp(4,in)+zvertp(3,in))/2.0

      endif

      return
      end

c
c
c
      subroutine b2repl(mrings,mkpts,cutring,cutpt1,cutpt2,readaux,
     >                  rizb,crmb,cion,ix_cell_offset)
      implicit none
      integer mrings,mkpts,cutring,cutpt1,cutpt2,readaux,cion,
     >        ix_cell_offset
      real    rizb,crmb
      include 'params'
      include 'cgeom'
      include 'cedge2d'
c
c     B2REPL:
c
c     This subroutine is a duplicate of the B2WRPL routine except
c     that the subroutine called to manipulate the data READ's instead
c     of WRITE's.
c
c     The values are read into temporary variables which match the B2
c     declarations ... which are then mapped onto the DIVIMP meshed
c     variables. These values can then be used later by DIVIMP in the
c     normal fashion. The MAPTODIV routine handles this mapping. It
c     can be invoked to do either cell centre direct mapping or to
c     average cell edge values and assign them to the DIVIMP cell
c     centred quantity - as happens with the velocity.
c
c     Note: Where possible ... only one array is used for multiple
c     quantities in order to reduce the amount of storage required for
c     temporary variables.
c
c
c     nfla = number of ion species in the Braams data file ... this
c     can be hard-coded or read-in from the input datafile.
c
c
c     temporary variable file ?
c
c      include 'commgp.f'
c
c     NFLA -> cgeom
c     MAXNFLA -> params
c
c     MAXNFLA - limits the storage allocated to the temporary variable used
c     to hold the B2 plasma solution before being copied into the
c     DIVIMP arrays.
c
c     NFLA - indicates the number of fluids present in the B2 solution and
c     is entered in the input file NFLA <= MAXNFLA
c
c      integer nfla
c
c     Single fluid background
c
c      parameter (nfla=1)
c
c     Multi-fluid background - Carbon + Hydrogen
c
c     parameter (nfla=7)
c
      integer maxix,maxiy,maxis
      parameter (maxix=maxnks,maxiy=maxnrs,maxis=maxnfla)
c
      real tempvol(maxnks,maxnrs)
      integer kp,l,lp1

c
      real ndummy(0:maxix+1,0:maxiy+1,maxis)
      real tdummy(0:maxix+1,0:maxiy+1)
c
      integer ix,iy,is,nplasf,nx,ny,nxd,nyd,ir,nplasaux
      integer ik,iz,ios
      real    tmpne
      character*200 buffer
c
c     Unit number for Braams data - this could also be read from the
c     datafile
c
      parameter (nplasf=11,nplasaux=12)
c
c      external gfsub3r
c
c     Set array size parameters for data loading routines. 
c
c     MAXKPTS=MKPTS - this quantity specifies the number of knots actually 
c     present in every row of the grid file and the plasma solution. When 
c     ix_cell_offset is non-zero the boundary cells at the ends of the rows
c     are stripped and the data modified so that the boundary cells located 
c     at the locations  ix_cell_offset and ix_cell_offset+1 are at the 
c     appropriate ends of the grid. However, for the code to work properly it
c     needs to believe that the data that has been read in and reorganized 
c     actually applies to a grid that has 2 fewer cells in each row. To do this
c     the value of mkpts is modified here at the beginning of the read routine. 
c     A non-zero value of ix_cell_offset implies two things - first the 
c     normal boundary cells must be removed and second that the grid and data 
c     must be rearranged so that the boundary cells emedded in the grid are found
c     at the ends. This data manipulation is done in the gfsub3r routine. The 
c     value of mkpts is restored to its original value at the end of this
c     routine.
c



      do ir = 1,nrs
         do ik = 1,nks(ir)

            KP = KORPG(IK,IR)
            tempvol(IK,IR) = 0.0
            IF (KP.GT.0) THEN
            DO L = 1, NVERTP(KP)
                 LP1 = L + 1
                 IF (L.EQ.NVERTP(KP)) LP1 = 1
                 tempvol(IK,IR) = tempvol(IK,IR) 
     >                     + (RVERTP(LP1,KP)*ZVERTP(L,KP)
     >                       - RVERTP(L,KP)*ZVERTP(LP1,KP))
            ENDDO
c
c           Ensure that the area is greater than zero.
c
            tempvol(IK,IR) = 0.5 * abs(tempvol(IK,IR))
            ENDIF
         end do 
      end do


      if (ix_cell_offset.gt.0) then 
         mkpts = mkpts -2
      endif
c
      nx = mkpts-2
      ny = mrings-2
c
      nxd = maxix
      nyd = maxiy
c
c     Initialization
c
      cre2dizs = -1
c
      write (6,*) 'NFLA:',nfla
c
c     calls to read routine
c
c     species densities  (ni)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)
c
      call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >         nfla,maxnfla,ndummy(0,0,1),knbs,maxnks,maxnrs,1.0,0)
c
c     Load the impurity species data if available
c     NOTE: B2E has started running case with multiple impurities in the
c           solution - we are ONLY interested for now in the profiles
c           of the impurity that is being run in DIVIMP - to facilitate this
c           the following code has been modified to selectively load the
c           appropriate densities when possible. 
c
      if (nfla.gt.1) then
c
c        Set number of charge states recorded in E2DNZS array
c
c         cre2d    = 2
c
         cre2dizs = nfla-1
c
         do iz = 1,nfla-1
c
           call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >         nfla,maxnfla,ndummy(0,0,iz+1),e2dnzs(1,1,iz),
     >         maxnks,maxnrs,1.0,0)
c
         end do
c
      endif
c
c     poloidal velocity  (uu)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)
c
c     radial velocity    (vv)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)
c
c
c     electron temperature (te)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
      call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >       1,1,tdummy(0,0),ktebs,maxnks,maxnrs,1.0/1.6e-19,0)
c
c      CALL PRRMATDIV(KTEBS,MAXNKS,nks(irsep),NRS,6,'TE')
c
c     ion temperature      (ti)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
      call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >       1,1,tdummy(0,0),ktibs,maxnks,maxnrs,1.0/1.6e-19,0)
c
c      CALL PRRMATDIV(KTIBS,MAXNKS,nks(irsep),NRS,6,'TI')
c
c     unknown              (pr)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c
c     parallel velocity    (up)
c     This is supposed to be at the cell boundaries ... this
c     should mean that the array is 1 element larger on each ring.
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)

      if (fc_v_interp_opt.eq.0) then  
c
c        Map as cell boundary velocity to cell boundary - into e2dbvel
c
         call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >      nfla,maxnfla,ndummy(0,0,1),e2dbvel,maxnks+1,maxnrs,1.0,0)
c
c        Map as cell boundary to cell centre velocity - into kvhs 
c
         call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >         nfla,maxnfla,ndummy(0,0,1),kvhs,maxnks,maxnrs,1.0,1)
c
      elseif (fc_v_interp_opt.eq.1) then 
c
c        Map as cell centre to cell centre velocity - into kvhs 
c
         call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >         nfla,maxnfla,ndummy(0,0,1),kvhs,maxnks,maxnrs,1.0,0)
c
      endif

c
c
c     Load the impurity species velocity data if available
c
      if (nfla.gt.1) then
c
         do iz = 1,nfla-1

c
c           do ix = 0,nx
c              do iy = 0,iy
c                 write (6,'(a,3i4,1x,g12.4)') 'V-dummy:',ix,iy,iz,
c     >                           ndummy(ix,iy,iz+1)
c              end do
c           end do
c
c
           call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >         nfla,maxnfla,ndummy(0,0,iz+1),e2dvzs(1,1,iz),
     >         maxnks,maxnrs,1.0,1)

c
         end do
c
c
c        Extract the impurity ion velocities from the westout files.
c
c        Copy Hydrogen ion velocity to the neutral entry of the
c        e2dvzs array.
c
c         write (6,*) 'E2DVZS:'
c
c         do ir = 1,nrs
c            do ik = 1,nks(ir)+1
c               e2dvzs(ik,ir,0) = kvhs(ik,ir)
c
c               write (6,'(2i4,7(1x,g12.4))')
c     >             ik,ir,(e2dvzs(ik,ir,iz),iz=0,6)
c
c            end do
c         end do
c
      endif
c
c     Bthet/Btot ratio     (pit)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c     unknown              (fnix)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)
c 
c     Map fnix as cell boundary quantity  
c
      call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >      nfla,maxnfla,ndummy(0,0,1),e2dflux,maxnks+1,maxnrs,1.0,0)
c
c     unknown              (fniy)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,nfla,maxnfla,ndummy(0,0,1),
     >             ix_cell_offset)
c
c
c     unknown              (feix)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c     unkonown             (feiy)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c     unknown              (feex)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c     unknown              (feey)
c
      call gfsub3r(nplasf,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
c
c     The electric field is calculated later - after KSS
c     is calculated
c
c-----------------------------------------------------------------
c
c     Read in data from an auxiliary input file - if one is specified
c
c     The quantities in the auxiliary plasma data file are:
c
c     1) Neutral Hydrogen Density
c     2) Neutral Carbon/Impurity Density
c     3) Recombination source rate
c     4) CX-recombination source rate
c     5) Impurity ionization rate for C0->C1+
c
      if (readaux.eq.1) then
c
c        Convert this code to a keyed loop that continues through the 
c        file until an EOF is reached - also to support old files - check 
c        to see if the first line is an AUX: tag and if it is not - just 
c        read in using the older methods.
c

            read(nplasaux,'(a200)',end=200,err=200) buffer
c
         if (buffer(1:20).eq.'FLUID CODE AUX FILE:') then 
c
            ios = 0
c
            do while (ios.eq.0) 

               read(nplasaux,'(a200)',end=200,err=200,
     >              iostat=ios) buffer
c
c              Read in Neutral Hydrogen Density - save in E2DATOM -
c                    Copied to KNHS in CXREC.
c
               if (buffer(1:11 ).eq.'H0 DENSITY:') then 
                  call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
                  call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                 1,1,tdummy(0,0),e2datom,maxnks,maxnrs,1.0,0)
               endif
c
c              Read in neutral impurity density
c
               if (buffer(1:20).eq.'IMP NEUTRAL DENSITY:') then 

                  call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
                  call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                 1,1,tdummy(0,0),e2dz0,maxnks,maxnrs,1.0,0)
c
c                 Copy into e2dnzs(ik,ir,0)
c
                  do ir = 1,maxnrs
                     do ik = 1,maxnks
                        e2dnzs(ik,ir,0) = e2dz0(ik,ir)
                     end do
                  end do
c
               endif
c
c              Read in C+ regular recombination rate
c
               if (buffer(1:10).eq.'C+ EI REC:') then 
 
                  call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
                  call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                 1,1,tdummy(0,0),e2drec,maxnks,maxnrs,1.0,0)
               endif 

c
c              Read in C+ CX recombination rate
c
               if (buffer(1:10).eq.'C+ CX REC:') then 
                  call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
                  call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                 1,1,tdummy(0,0),e2dcxrec,maxnks,maxnrs,1.0,0)
               endif 

c
c              Read in C0->C1+ ionization rate
c              May be copied to PINIONZ for injection option 7.
c
               if (buffer(1:) .eq. 'C0->C+ IONIZATION:') then 

                  call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
                  call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                 1,1,tdummy(0,0),e2diz0,maxnks,maxnrs,1.0,0)
c
              endif 
c
c
            end do
c
         else
c
            backspace nplasaux
c
c           Read in Neutral Hydrogen Density - save in E2DATOM -
c                    Copied to KNHS in CXREC.
c
            call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
            call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >          1,1,tdummy(0,0),e2datom,maxnks,maxnrs,1.0,0)
c
c           Read in neutral impurity density
c
            call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
            call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >          1,1,tdummy(0,0),e2dz0,maxnks,maxnrs,1.0,0)
c
c           Copy into e2dnzs(ik,ir,0)
c
            do ir = 1,maxnrs
               do ik = 1,maxnks
                  e2dnzs(ik,ir,0) = e2dz0(ik,ir)
               end do
            end do
c
c           Read in C+ regular recombination rate
c
            call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
            call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >          1,1,tdummy(0,0),e2drec,maxnks,maxnrs,1.0,0)
c
c           Read in C+ CX recombination rate
c
            call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
            call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >          1,1,tdummy(0,0),e2dcxrec,maxnks,maxnrs,1.0,0)
c
c           Read in C0->C1+ ionization rate
c           May be copied to PINIONZ for injection option 7.
c
            call gfsub3r(nplasaux,nx,ny,nxd,nyd,1,1,tdummy(0,0),
     >             ix_cell_offset)
            call maptodiv(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >          1,1,tdummy(0,0),e2diz0,maxnks,maxnrs,1.0,0)
c
         endif
c
 200     continue     
c
c        Convert these values in the rec and cxrec arrays from
c        particles/s to particles/m-toroidally/s
c
c        Later convert e2diz0 to a density when cell areas are known
c
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               if (rs(ik,ir).ne.0.0) then
                  e2drec(ik,ir)  = e2drec(ik,ir)  /(2.0*PI*rs(ik,ir))
                  e2dcxrec(ik,ir)= e2dcxrec(ik,ir)/(2.0*PI*rs(ik,ir))
                  e2diz0(ik,ir)  = e2diz0(ik,ir)  /(2.0*PI*rs(ik,ir))
               endif
            end do
         end do
c
      endif
c

c
c     Now that all of the fluid code results have been loaded
c     - try to figure out which section of the array is needed
c       and move the rest of the data to compensate. 
c     - if (nfla-1 = cion) then assume that the correct 
c       species is the only one in the fluid code results file. 
c
c      if (nfla-1.lt.cion) then 
c
c         write(0,*) 'DATA FOR IMPURITY OF ATOMIC NUMBER:', cion,
c     >                 ' IS NOT IN THE FLUID CODE DATA FILE'
c
c      elseif (nfla-1.gt.cion.and.cion.eq.6.and.nfla.ne.17) then 
c
c        Note: 17 fluids is usually H+C+Ne and so does not need 
c              adjusting.
c
c        The fluid code results file contains additional fluids
c        beyond the one that is being looked at here. Need to select
c        which section of the array to use.   
c
c        Also - assume for now that we are only looking for Carbon
c        data from the fluid code file. In this case there is likely 
c        only helium before the Carbon in the file and possibly other
c        species after it - shift the array down by two charge states 
c        to remove the He and then set the upper bound to just C.
c
c         write (6,*) 'Eliminating Fluids:',cion
c
c         cre2dizs = cion             
c
c         do ir = 1,nrs
c            do ik = 1,nks(ir) 
c               do iz = 3,3+cion-1
c                  e2dnzs(ik,ir,iz-2) = e2dnzs(ik,ir,iz)
c                  e2dvzs(ik,ir,iz-2) = e2dvzs(ik,ir,iz)
c               end do
c            end do
c         end do
c
c      endif 
c
c
      if (nfla-1.lt.cion) then 
c
         write(0,*) 'DATA FOR IMPURITY OF ATOMIC NUMBER:', cion,
     >                 ' IS NOT IN THE FLUID CODE DATA FILE'
c
      elseif (nfla-1.gt.cion) then 
c
c        Limit number of charge states to only those of interest and assume
c        that they are first in the fluid file.  
c
         cre2dizs = cion
c
      endif
c
c     Continue processing 
c
      if (cre2dizs.gt.0) then

         do ir = 1,nrs
            do ik = 1,nks(ir)

               tmpne = knbs(ik,ir)

               do iz = 1,cre2dizs

                  tmpne = tmpne + iz * e2dnzs(ik,ir,iz)

               end do

               if (tmpne.gt.1.1*knbs(ik,ir)) then
                  write (6,*) '***NOTE***'

                  write(6,'(a,2i4,4(1x,g12.5))')
     >               'IONIZ:',ik,ir,knbs(ik,ir),tmpne,
     >               e2dnzs(ik,ir,1),e2dnzs(ik,ir,2)
               endif
c
c              Assign to fluid code array
c
               e2dnes(ik,ir) = tmpne
c
            end do

         end do

      endif
c
c     Copy background that has been loaded into the specific fluid code
c     arrays and assign appropriate values to the e2dtarg array. 
c
      call rzero(e2des,maxnks*maxnrs)
c
c
      write(6,*) 'Fluid code solution:'
c
      do ir = 1,nrs
c
         do ik = 1,nks(ir)
c
            e2dnbs(ik,ir) = knbs(ik,ir)
            e2dtebs(ik,ir)= ktebs(ik,ir)
            e2dtibs(ik,ir)= ktibs(ik,ir)
            e2dvhs(ik,ir) = kvhs(ik,ir)
c
c
            write (6,'(a,2i4,10g12.4)') 'FC:',ik,ir,e2dnbs(ik,ir),
     >                  e2dtebs(ik,ir),e2dtibs(ik,ir),e2dvhs(ik,ir),
     >                  e2des(ik,ir),kvhs(ik,ir),e2dbvel(ik,ir),
     >                  e2dflux(ik,ir),tempvol(ik,ir)
c
         end do
c
         write(6,*)
c
      end do         

c
c     Calculate the FLUID CODE target conditions
c     This must be done before the virtual points are
c     stripped away.
c     
c     The values in E2Dtarg are:
c     1 - Density
c     2 - Te
c     3 - Ti
c     4 - Vb
c     5 - Gamma = Ne * Vb
c     
c     Note: Vb is calculated as the sound speed for the given target
c           temperatures OR the cell boundary velocity
c     
c     NOTE: NKS(IR) has NOT yet been modified for the removal of
c           any virtual points (when flag=0) and thus nks(ir)
c           will point to the virtual point.
c     
c
c     NOTE: For UEDGE target conditions - pull the data from the 
c           virtual cells for Te and Ti.
c
c     
      write(6,*) 'FC Target Conditions:'
c     
      do ir = irsep,nrs

c
c        Target 1 - this is  the OUTER target for most sonnet grids
c                   or grids with the X-POINT at the bottom
c        Target 2 - this is  the INNER target for most sonnet grids
c                   or grids with the X-POINT at the bottom
c
c        The indexing of the e2dtarg array is:
c        e2dtarg(ring#, quantity index#, target#) 
c
c
c        Ne 
c
         if (fc_ne_calc_opt.eq.0) then   
            e2dtarg(ir,1,2) = e2dnbs(1,ir)
            e2dtarg(ir,1,1) = e2dnbs(nks(ir),ir)
         elseif (fc_ne_calc_opt.eq.1) then  
            e2dtarg(ir,1,2) = (e2dnbs(1,ir)+e2dnbs(2,ir))/2.0
            e2dtarg(ir,1,1) = (e2dnbs(nks(ir),ir)
     >                        +e2dnbs(nks(ir)-1,ir))/2.0
         elseif (fc_ne_calc_opt.eq.2) then  
            e2dtarg(ir,1,2) = e2dnbs(2,ir)
            e2dtarg(ir,1,1) = e2dnbs(nks(ir)-1,ir)
         endif    
c
c        Te
c
         if (fc_te_calc_opt.eq.0) then   
            e2dtarg(ir,2,2) = e2dtebs(1,ir)
            e2dtarg(ir,2,1) = e2dtebs(nks(ir),ir)
         elseif (fc_te_calc_opt.eq.1) then  
            e2dtarg(ir,2,2) = (e2dtebs(1,ir)+e2dtebs(2,ir))/2.0
            e2dtarg(ir,2,1) = (e2dtebs(nks(ir),ir)
     >                       +e2dtebs(nks(ir)-1,ir))/2.0
         elseif (fc_te_calc_opt.eq.2) then  
            e2dtarg(ir,2,2) = e2dtebs(2,ir)
            e2dtarg(ir,2,1) = e2dtebs(nks(ir)-1,ir)
         endif    
c
c        Ti
c
         if (fc_ti_calc_opt.eq.0) then   
            e2dtarg(ir,3,2) = e2dtibs(1,ir)
            e2dtarg(ir,3,1) = e2dtibs(nks(ir),ir)
         elseif (fc_ti_calc_opt.eq.1) then  
            e2dtarg(ir,3,2) = (e2dtibs(1,ir)+e2dtibs(2,ir))/2.0
            e2dtarg(ir,3,1) = (e2dtibs(nks(ir),ir)
     >                       +e2dtibs(nks(ir)-1,ir))/2.0
         elseif (fc_ti_calc_opt.eq.2) then  
            e2dtarg(ir,3,2) = e2dtibs(2,ir)
            e2dtarg(ir,3,1) = e2dtibs(nks(ir)-1,ir)
         endif    
c
c        Vb
c
         if (fc_v_calc_opt.eq.0) then  
c
            e2dtarg(ir,4,2) = kvhs(1,ir)
c
c           The code used to use nks(ir)-1 - however, I do not
c           believe that this corresponds to the virtual cell at the
c           second target where the actual velocity at the cell 
c           boundary was to be stored. Since at this point the cells
c           number from 1 to nks(ir)+1 since the velocity is a cell 
c           boundary quantity. 
c
c           e2dtarg(ir,4,1) = kvhs(nks(ir)-1,ir)
c
            e2dtarg(ir,4,1) = kvhs(nks(ir),ir)
c
            write(6,'(a,4(1x,g12.5))') 'FC:V:',
     >               kvhs(nks(ir)-1,ir),kvhs(nks(ir),ir),
     >               e2dbvel(nks(ir)-1,ir),e2dbvel(nks(ir),ir)
c
c        Vb = Cs 
c
         elseif (fc_v_calc_opt.eq.1) then  

             e2dtarg(ir,4,2) = 
     >          -9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
     >               *(1.0+rizb)/crmb)

             e2dtarg(ir,4,1) = 
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
     >               *(1.0+rizb)/crmb)

         endif
c
c        Base flux value from Ne * Vb
c
         e2dtarg(ir,5,2) = e2dtarg(ir,1,2) * e2dtarg(ir,4,2)
         e2dtarg(ir,5,1) = e2dtarg(ir,1,1) * e2dtarg(ir,4,1)
c
c        Write out
c     
  101    format(a,i4,1p,4(1x,e13.5))
  102    format(a,i4,4(1x,f13.5))
c     
c     
c         if (cprint.eq.4.or.cprint.eq.9) then
c           write(6,*)   'NKS:',nks(ir),nj(ir)
           write(6,101) 'Ne:',ir,e2dtarg(ir,1,1),e2dtarg(ir,1,2)
           write(6,102) 'Te:',ir,e2dtarg(ir,2,1),e2dtarg(ir,2,2)
           write(6,102) 'Ti:',ir,e2dtarg(ir,3,1),e2dtarg(ir,3,2)
           write(6,101) 'Vb:',ir,e2dtarg(ir,4,1),e2dtarg(ir,4,2)
           write(6,101) 'Cs:',ir,
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
     >                    ),
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
     >                    )
           write(6,102) 'Ra:',ir,e2dtarg(ir,4,1)/
     >           (9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
     >                    )),
     >                           -e2dtarg(ir,4,2)/
     >           (9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
     >                    ))
c
           write(6,101) 'Ga:',ir,e2dtarg(ir,5,1),e2dtarg(ir,5,2)
c         endif
      
      end do
c
c     Restore mkpts value
c
      if (ix_cell_offset.gt.0) then 
         mkpts = mkpts -2
      endif
c
      return
      end
c
c
c
      subroutine gfsub3r(kard,nx,ny,ndimx,ndimy,
     >                   ns,ndims,dummy,ix_cell_offset)
      implicit none
c
      include 'params'
c
      real,allocatable :: dummy_temp(:,:,:)
c
      integer kard,nx,ny,ndimx,ndimy,ndims,nd1,lim,is,ix,iy,iii
      integer ns,ix_cell_offset
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)
c
      integer ierr
c
c     Increase value of nx to account for the two end boundary cells 
c     that are to be deleted.
c
      if (ix_cell_offset.ne.0) then 
c 
         nx = nx+2
c
      endif
c
c     READ the data from the plasma file
c
      nd1 = nx+2
      lim = (nd1/5)*5 - 4
      do 110 is = 1,ns
        do 110 iy = 0,ny+1
          do 100 ix = 1,lim,5
             read(kard,910,end=999)
     >                 (dummy(ix-2+iii,iy,is),iii=1,5)
 100      continue
          if ((lim+4).eq.nd1) goto 110
          read(kard,910,end=999) (dummy(ix-1,iy,is),ix=lim+5,nd1)
 110  continue

c
c     If ix_cell_offset is non-zero then:
c     (A) Strip boundary cells from the ends of the rows
c     (B) Rearrange the data in the array so that the internal
c         boundary cells at locations ix_cell_offset and ix_cell_offset+1 
c         are moved to the ends - ix_cell_offset goes to location ny-1 while
c         ix_cell_offset+1 goes to 1
c
      if (ix_cell_offset.ne.0) then 
c
c        Delete end boundary data by moving the array down in ix
c                    
         do is = 1,ns
            do iy = 0,ny+1
               do ix = 0,nx+1
                  dummy(ix,iy,is) = dummy(ix+1,iy,is)
               end do
            end do
         end do
c
c        Restore value of nx to its input value (which did not include the boundary cells)
c
         nx = nx -2 
c
c        Move the data in the array so that it is shifted by the amount ix_cell_offset
c        which will move the embedded boundary cells to the ends of the rows.
c
c        Allocate space for temp dummy array.
c
         allocate(dummy_temp(0:ndimx+1,0:ndimy+1,1:ndims),STAT=ierr)
c
         if (ierr.ne.0) then 
            write(0,*) 'ERROR:GFSUB3R: DUMMY array for FC'//
     >                 ' data could not be allocated: ',
     >                  ndimx,ndimy,ndims
            write(6,*) 'ERROR:GFSUB3R: DUMMY array for FC'//
     >                  ' data could not be allocated: ',
     >                  ndimx,ndimy,ndims
            stop
         endif        
c
c        Copy/move the data
c
         do is = 1,ns
            do iy = 0,ny+1
               do ix = 0,nx+1

                  if (ix.le.(ix_cell_offset-1)) then

                     dummy_temp(ix-(ix_cell_offset-1)+nx+1,iy,is) = 
     >                                        dummy(ix,iy,is)

                  else

                     dummy_temp(ix-(ix_cell_offset-1)-1,iy,is) =
     >                            dummy(ix,iy,is)

                  endif

               end do
            end do
         end do

         do is = 1,ns
            do iy = 0,ny+1
               do ix = 0,nx+1
                  dummy(ix,iy,is) = dummy_temp(ix,iy,is)
               end do
            end do
         end do

         deallocate(dummy_temp)
  
      endif 

      return
 910  format(5(e16.8))

 999  continue
      write (6,*) 'ERROR: EOF found in B2 formatted data file.'
      write (6,*) '       Check Geometry and number of fluids'//
     >            ' in data file for consistency.'            
      write (0,*) 'ERROR: EOF found in B2 formatted data file.'
      write (0,*) '       Check Geometry and number of fluids'//
     >            ' in data file for consistency.'            
      return
      end
c
c
c
      subroutine maptodiv(cutring,cutpt1,cutpt2,nx,ny,ndimx,ndimy,
     >            ns,ndims,dummy,divarr,dim1,dim2,scalef,valtype)
      implicit none
      integer cutring,cutpt1,cutpt2
      integer nx,ny,ndimx,ndimy,ndims,dim1,dim2,valtype,ns
      real divarr(dim1,dim2),scalef
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)
      include 'params'
      include 'cgeom'
c slmod begin
      include 'slcom'
c slmod end
c
c     MAPTODIV:
c
c     The purpose of this subroutine is to take the data read into
c     the dummy array from a Braams output file in subroutine
c     B2REPL and map it correctly into the SONNET geometry
c     that has already been read in. The characteristic sizes and
c     cut points of the ASDEX geometry are currently read-in in the
c     input data file.
c
c     This routine uses the same assumptions as are outlined in
c     the subroutine RAUG to determine the assignment of values
c     from the Braams arrays to the DIVIMP ones. At the present
c     time, the B2 code writes values for the virtual points on
c     each ring.
c
c     The background plasma data is in the first element of the
c     arrays ... other dimensions will contain the values for
c     impurity species if B2 was run as a multi-fluid code.
c
c     Valtype defines whether the array to be mapped contains
c     cell-centered quantities (valtype=0) or cell-edge ones
c     (valtype=1). In the case of cell-edge quantities the
c     routine saves the values at the first and last cell edge to
c     the first and last elements of the array - except for
c     core rings - since these are virtual points that will be
c     stripped. The rest of the elements are filled with the
c     average of the two adjacent edge quantities. Valtype=1
c     assumes that the first element in the array is the first
c     edge value at the boundary of the guard cell and that
c     the last cell in the ring contains an undefined value.
c
c
      integer ir,ik,mrings,mkpts,ikold,irold,ikact,iract
c
c     The maximum number of knots/ring in the DIVIMP case is nx+2
c     The maximum number of rings in the DIVIMP case is ny+2
c
      mkpts = nx + 2
      mrings= ny + 2
c
c
      ikold = 1
      irold = 1
c
c     Valtype = 0 - cell-centered mapping
c
      if (valtype.eq.0) then
c
         do 100 ir = 1,mrings
            do 100 ik = 1,mkpts
c
               if (ir.le.cutring) then
                  if (ik.le.cutpt1) then
                     iract = mrings + ir
                     ikact = ik
                  elseif (ik.ge.cutpt2) then
c
                     if (ik.eq.cutpt2) then
c
c                       Add a repeat point to the core
c                       plasma rings - to close them.
c
                        ikold = ikold +1
                        divarr(ikold,irold) = divarr(1,irold)
                     endif
c
                     iract = mrings + ir
                     ikact = ik - cutpt2 + cutpt1 +1
c
                  else
                     iract = ir
                     ikact = ik - cutpt1
                  endif
               else
                  iract = ir
                  ikact = ik
               endif

c               IF (kvols(ikact,iract).EQ.0.0) THEN
c               write (SLOUT,'(A,4I6,1P,G11.3,0P)') 
c     .                  'maptodiv 0:',ik-1,ir-1,ikact,iract,
c     >                  dummy(ik-1,ir-1,1)
c               ELSE
c               write (SLOUT,'(A,4I6,1P,G11.3,0P)') 
c     .                  'maptodiv 0:',ik-1,ir-1,ikact,iract,
c     >                  dummy(ik-1,ir-1,1)
cc     >                  dummy(ik-1,ir-1,1)/kvols(ikact,iract)*
cc     .                  rs(ik,ir)/rxp*1.6E-19*1.0E-06,kvols(ikact,iract)
c               ENDIF

               divarr(ikact,iract) = dummy(ik-1,ir-1,1)*scalef
               ikold = ikact
               irold = iract
 100     continue
c
c     Valtype =1 - cell-edge based quantities
c
      elseif (valtype.eq.1) then
c
c slmod begin - new
c
c See MAPTOB2 routine.
c
        CALL WN('MaptoDiv','Need corrections for edge quantities')
c slmod end
         do 200 ir = 1,mrings
            do 200 ik = 1,mkpts
c
               if (ir.le.cutring) then
                  if (ik.le.cutpt1) then
                     iract = mrings + ir
                     ikact = ik
c
                     if (ik.eq.1) then
                        divarr(ikact,iract) = dummy(ik-1,ir-1,1)*scalef
                     else
                        divarr(ikact,iract) =
     >                      0.5 *(dummy(ik-2,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
                     endif
c
                  elseif (ik.ge.cutpt2) then
c
                     iract = mrings + ir
                     ikact = ik - cutpt2 + cutpt1 +1
c
                     if (ik.eq.cutpt2) then
c
c                       Add a repeat point to the core
c                       plasma rings - to close them.
c
                        ikold = ikold +1
                        divarr(ikold,irold) = divarr(1,irold)
c
c                       Also deal with next average
c
                        divarr(ikact,iract) =
     >                      0.5 *(dummy(cutpt1-1,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
c
                     elseif (ik.eq.mkpts) then
c
c                       Copy boundary value to last cell.
c
                        divarr(ikact,iract) = dummy(ik-2,ir-1,1)*scalef
c
                     else
c
                        divarr(ikact,iract) =
     >                      0.5 *(dummy(ik-2,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
c
                     endif
c
                  else
                     iract = ir
                     ikact = ik - cutpt1
c
                     if (ik.eq.cutpt1+1) then
c
                        divarr(ikact,iract) =
     >                      0.5 *(dummy(cutpt2-2,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
c
                     else
c
                        divarr(ikact,iract) =
     >                      0.5 *(dummy(ik-2,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
c
                     endif
                  endif
               else
c
                  iract = ir
                  ikact = ik
c
                  if (ik.eq.1) then
c
                     divarr(ikact,iract) = dummy(ik-1,ir-1,1)*scalef
c
                  elseif (ik.eq.mkpts) then
c
                     divarr(ikact,iract) = dummy(ik-2,ir-1,1)*scalef
c
                  else
c
                     divarr(ikact,iract) =
     >                      0.5 *(dummy(ik-2,ir-1,1)
     >                          + dummy(ik-1,ir-1,1))*scalef
c
                  endif
c
               endif
               ikold = ikact
               irold = iract
 200     continue
c
      endif
c
      return
      end
c
c
c
      subroutine calcefb2
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     CALCEFB2:
c
c     This subroutine uses a simple prescription to estimate
c     the electric field from the background density and
c     temperature and the distances between grid points.
c
c     The densities and temperatures will have been read in from
c     the converged results of a B2-Eirene run. The calculation is
c     not performed when the data is read in because the KSS values
c     of the grid points have not yet been calculated at that
c     time.
c
      real ds1,dp1,dt1,nb1,ds2,dp2,dt2,nb2,irlimit
      INTEGER ik,ir
C
C       CALCULATE ELECTRIC FIELD
C
C
C       IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
C       SAME FACTOR USED IN CONVERTING T IN EV TO KT.
C
      call rzero(kes,maxnks*maxnrs)
c
c     Calculate electric field only if the OFIELD option
c     is turned off.
c
c      if (ofield.eq.0) then
c
c
      if (ciopto.eq.0.or.ciopto.eq.2.or.
     >    ciopto.eq.3.or.ciopto.eq.4) then
         irlimit = irwall
c slmod begin
      elseif (ciopto.eq.1.or.ciopto.eq.5) then
c
c      elseif (ciopto.eq.1) then
c slmod end
         irlimit = nrs
      endif
c
c
c
      do 600 ir = irsep,irlimit
        DS1 = KSS(2,IR) - KSS(1,IR)
        DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
        DT1 = (KTEBS(2,IR)-KTEBS(1,IR))
        NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
C
        if ((ds1 .ne. 0.0) .and. (nb1.ne.0.0)) then 
           KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
        else
           kes(1,ir) = 0.0
           write(6,'(a,2i8,4(1x,g12.5))') 
     >             'KES CALCULATION: IK,IR,DS1,NB1:',
     >              1,ir,ds1,nb1
        endif

C
        DS1 = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
        DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >         -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
        DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
        NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))
C
c        KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
        if ((ds1 .ne. 0.0) .and. (nb1.ne.0.0)) then 
           KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
        else
           kes(1,ir) = 0.0
           write(6,'(a,2i8,4(1x,g12.5))') 
     >             'KES CALCULATION: IK,IR,DS1,NB1:',
     >              nks(ir),ir,ds1,nb1
        endif
C
C        WRITE(6,*) 'KES:IR:',KES(1,IR),KES(NKS(IR),IR)
C
        DO 500 IK = 2,NKS(IR)-1

            DS1 = KSS(IK,IR) - KSS(IK-1,IR)
            DP1 = KNBS(IK,IR)*KTEBS(IK,IR)-KNBS(IK-1,IR)*KTEBS(IK-1,IR)
            DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
            NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
            DS2 = KSS(IK+1,IR) - KSS(IK,IR)
            DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)-KNBS(IK,IR)*KTEBS(IK,IR)
            DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
            NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
           if ((ds1 .ne. 0.0) .and. (nb1.ne.0.0).and.
     >         (ds2 .ne. 0.0) .and. (nb2.ne.0.0)) then 
               KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
           else
              kes(1,ir) = 0.0
              write(6,'(a,2i8,4(1x,g12.5))') 
     >            'KES CALCULATION: IK,IR,DS1,NB1,ds2,nb2:',
     >              ik,ir,ds1,nb1,ds2,nb2
           endif
c            KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
c     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
C
C            WRITE(6,*) 'KES:',IK,IR,KES(IK,IR)
C
 500    CONTINUE
c
c       Set EFIELD to zero for the midpoint of the ring where inner 
c       and outer solutions join
c
        kes(ikmids(ir),ir) = 0.0
        kes(ikmids(ir)+1,ir) = 0.0
c

600   CONTINUE

c
c     End of OFIELD if
c
c      endif
c
c
      return
      end
c
c
c
      subroutine calcleq
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'comsol'
c slmod begin
      include 'slcom'
c slmod end
c
c     This routine calculates the equivalent source lengths
c     for the PIN ionization source fro each ring at each plate.
c     This data is used to characterize the source for comparison
c     to the DIVIMP models.
c
c     The Lequiv is calculated from the equation ...
c     Leq GammaTarget = INT(0 to smax/2) {Gamma(s)ds}
c
c     Which is equivalent to ...
c
c     Lequiv = Smax/2 - INT INT (0 to Smax/2) {S(s)ds} /
c                                INT (0 to smax/2) {S(s)ds}
c
c
c     This routine utilizes the same routines as are used in the
c     Soledge subroutine to calculate the BG plasma characteristics.
c
c
      integer          ir,i,sopt
      double precision srcion,dist,smax,ionint2,ionint2i,cis1,cis2
      external         srcion,cis1,cis2
c
c
      do 100 ir = irsep,nrs
c slmod begin
         IF (idring(ir).EQ.BOUNDARY) CYCLE

         IF (osm_model(IKLO,ir).EQ.28) CYCLE
c slmod end
c
c        Setup the values required for the ionization routine to
c        work correctly.
c
         smax = ksmaxs(ir)
         dist = smax/2.0
c
         LENSRC = SMAX /2.0
         LAMSRC = SMAX /2.0
         IRN = IR
         MAXIK = NKS(IR)
         TMAXS = smax
         DO 500 I = 1,MAXIK
           TSS(I) = KSS(I,IR)
500      CONTINUE
         FNORM = 1.0
         FNORMI = 1.0
c
c        The ionization source value S(s) can be provided by
c        the SRCION function if sopt = 2 is specified.
c
         sopt = 2
c
c        Set-up the values needed for source function
c        integration.
c
c        Ion source
c
         sionl = 0.0d0
         pionl = 0
         iionl = 0.0d0
c
c
         ionint  = cis1(dist,srcion,sopt,0)
         ioninti = cis1(dist,srcion,sopt,1)
c
         ionint2 = cis2(dist,srcion,sopt,0)
         ionint2i= cis2(dist,srcion,sopt,1)
c
         if (ionint.eq.0.0) then
            cleq(ir,2) = 0.0
         else
            cleq(ir,2) = smax/2.0
     >           - ionint2/ionint
         endif
c
         if (ioninti.eq.0.0) then
            cleq(ir,1) = 0.0
         else
            cleq(ir,1) = smax/2.0
     >           - ionint2i/ioninti
         endif
         if (cprint.eq.9)
     >      write(6,*) 'cleq:',ir,smax/2.0,ionint,ionint2,
     >          cleq(ir,2) ,ioninti,ionint2i,cleq(ir,1)
c
c
c
 100  continue
c
c
c
      return
      end
c
c
c
      subroutine calcorth
      implicit none
c
c     This subroutine calculates the orthogonality
c     characteristic of each grid cell with a complete
c     polygon defined. Ones which have odd shapes
c     are flagged for future examination. In addition,
c     a correction factor for cross-field transport
c     can be calculated that can be applied when moving
c     in non-orthogonal grid regions.
c
c     This routine also calculates the orthogonal angle for each
c     target segment.
c
c
c     David Elder,   1994 March 24
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     Local Variables
c
      real rp1,zp1,rp2,zp2,rp3,zp3,rp4,zp4
      real ri,zi,m1,m2,b1,b2
      real asqr,bsqr,csqr
      real cr1,cr2,cz1,cz2,ds1,ds2,ds3,ds4
      real rdum1
      integer ik,ir,in,kp,id
c
c     Externals
c
      real atan2c,acos_test
      external atan2c,acos_test
c
c
      do 10 ir = 1, nrs
        do 10 ik = 1,nks(ir)
          kp = korpg(ik,ir)
          if (kp.gt.0) then
            if (nvertp(kp).eq.4) then
c
              rp1 = (rvertp(1,kp) + rvertp(2,kp))/2.0
              rp2 = (rvertp(2,kp) + rvertp(3,kp))/2.0
              rp3 = (rvertp(3,kp) + rvertp(4,kp))/2.0
              rp4 = (rvertp(4,kp) + rvertp(1,kp))/2.0
c
              zp1 = (zvertp(1,kp) + zvertp(2,kp))/2.0
              zp2 = (zvertp(2,kp) + zvertp(3,kp))/2.0
              zp3 = (zvertp(3,kp) + zvertp(4,kp))/2.0
              zp4 = (zvertp(4,kp) + zvertp(1,kp))/2.0
c
              if (rp3-rp1.ne.0.0) then 
                 m1 = (zp3-zp1)/(rp3-rp1)
              else
                 m1 = 0.0
              endif
              b1 = zp1 - m1 * rp1
c
              if (rp4-rp2.ne.0.0) then 
                 m2 = (zp4-zp2)/(rp4-rp2)
              else
                 m2 = 0.0
              endif
              b2 = zp2 - m2 * rp2
c
c             Code modified to deal with degenerate or 
c             zero volume cells.  
c
              if ((rp3.eq.rp1).and.
     >            (rp4.eq.rp2)) then  
c                
c                Force bsqr to zero to get default solution
c
                 ri = rp3
                 zi = zp3   
c
              elseif (rp3.eq.rp1) then
c
c                Case of m1 = INF (vertical line)
c                ri = rp3 = rp1
c
                 ri = rp3
                 zi = m2 * ri + b2
c
              elseif (rp4.eq.rp2) then
c
c                Case of m2 = INF (vertical line)
c                ri = rp4 = rp2
c
                 ri = rp4
                 zi = m1 * ri + b1
c
              elseif (m1.eq.m2) then 
c                
c                Force bsqr to zero to get default solution
c
                 ri = rp3
                 zi = zp3   
c
              else
                 ri = (b2-b1) / (m1-m2)
                 zi = m1 * ri + b1
              endif
c
              asqr = (rp3-rp4)**2 + (zp3-zp4)**2
c
c             For zero or almost zero volume cells - eliminate the 
c             problem of numerical inaccuracy causing a false cell 
c             orthogonality by setting extremely small values of bsqr and
c             csqr to zero - this will result in the cell being reported 
c             as orthogonal. 
c
              bsqr = (rp3-ri )**2 + (zp3-zi )**2
              if (bsqr.lt.1.0e-12) bsqr = 0.0

              csqr = (rp4-ri )**2 + (zp4-zi )**2
              if (csqr.lt.1.0e-12) bsqr = 0.0

c
c              write(6,'(a,2i6,10(1x,g12.5))') 'CALC_ORTH:',
c     >           ik,ir,asqr,bsqr,csqr,rp3-rp1,rp4-rp2,m1,m2
c
              if (bsqr.le.0.0.or.csqr.le.0.0) then
                 cosalph(ik,ir) = 0.0
                 sinalph(ik,ir) = 1.0
              else
                 cosalph(ik,ir)=(asqr-bsqr-csqr)/
     >                          (-2.0*sqrt(bsqr*csqr))
                 sinalph(ik,ir)= 1.0 - cosalph(ik,ir)**2
              endif
c
c             Calculate the orthogonal angle for target elements.
c
              if (ir.ge.irsep.and.ik.eq.1) then
c
c                Target number 2 - outer for JET , inner for Sonnet
c
c                Corners 2,1 form the target
c
                 id = idds(ir,2)
c
                 target_orth_angle(id) =
     >                 atan2c(zs(ik,ir)-zp1,rs(ik,ir)-rp1)
c
              elseif (ir.ge.irsep.and.ik.eq.nks(ir)) then
c
c                Target number 1 - inner for JET , outer for Sonnet
c
c                Corners 3,4 form the target
c
                 id = idds(ir,1)
c
                 target_orth_angle(id) =
     >                 atan2c(zs(ik,ir)-zp3,rs(ik,ir)-rp3)
c
              endif
c
            else
              cosalph(ik,ir) = 0.0
              sinalph(ik,ir) = 1.0
            endif
c
          else
             cosalph(ik,ir) = 0.0
             sinalph(ik,ir) = 1.0
          endif
c sl begin
          ALPH(IK,IR) = ACOS_TEST(COSALPH(IK,IR),17) * 180.0 / PI
c sl end
c
c
 10   continue

c
c     Print out results for debugging
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         write (6,*) 'Calculate Cell Orthogonality'
         do 20 ir = 1,nrs
           write (6,*) 'Ring Number:',ir
           do 20 ik = 1,nks(ir)
             kp = korpg(ik,ir)
             if (kp.gt.0) then
              rdum1 = acos_test(cosalph(ik,ir),18)  ! Avoids recursive I/O error - SL, 04.03.2010
              write(6,'(2i4,3(1x,g13.6),2i4)') ir,ik,cosalph(ik,ir),
     >          sinalph(ik,ir),raddeg*rdum1,
     >          kp,nvertp(kp)
             endif
 20      continue
c
c        Print target orthogonal angles
c
         write (6,*)
         write (6,*) 'Calculate Target Orthogonal Angles',nds
         do id = 1,nds
           write(6,'(3i4,3(1x,g13.6))') id,irds(id),ikds(id),
     >        raddeg*thetas(id),raddeg*target_orth_angle(id),
     >        raddeg*abs(thetas(id)-target_orth_angle(id))
c
         end do
         write (6,*)

c
      endif

      return
      end
c
c
c
      subroutine writegrd(cgridopt)
      implicit none
      integer cgridopt
      include 'params'
      include 'cgeom'
c
c     This routine writes out information about the grid to a
c     separate file. Currently assigned to fort.25. It includes
c     vertices - cell centre information as well as magnetic field.
c
c     This information is only printed if the expanded print option
c     has been selected. It is printed in a format compatible with
c     the SONET grid generator output.
c
c     David Elder, March 8, 1995
c
      integer iounit,iounit2
      parameter (iounit= 25,iounit2=27)
c
      integer ik,ir,in,ic
c
      integer tmpkorpg(maxnks,maxnrs), tmpnvertp(maxnks*maxnrs)
      integer tmpnks(maxnrs)
c
      real tmprvertp(5,maxnks*maxnrs),tmpzvertp(5,maxnks*maxnrs)
c
      write(iounit,6001)
      write(iounit,6011)
      write(iounit,6005)
c
c     Write out grid data - map to the indexing used for SONNET grids.
c
c
c     Generate "fake" polygons for cells that do not have them.
c
c      call rzero(tmpkorpg,maxnrs*maxnks)
c      call rzero(tmpnvertp,maxnrs*maxnks)
c      call rzero(tmprvertp,5*maxnrs*maxnks)
c      call rzero(tmpzvertp,5*maxnrs*maxnks)
c
c     Copy equilibrium data
c
c      do ir = 1,nrs
c         tmpnks(ir) = nks(ir)
c         do ik = 1,nks(ir)
c            in = korpg(ik,ir)
c            tmpkorpg(ik,ir) = in
c            do ic = 1,5
c               tmprvertp(ic,in) = rvertp(ic,in)
c               tmpzvertp(ic,in) = rvertp(ic,in)
c            end do
c         end do
c      end do
c
c     Shift grid and generate extra cells.
c
c      do ir = irsep,nrs
c         do
c         tmpkorpg(ik,
c





c
c
c 200  read(gridunit,9000,end=300)in,ik,ir,
c     >                             rvert(2),zvert(2),rvert(3),zvert(3)
c      read(gridunit,9001,end=2000) brat,rcent,zcent
c      read(gridunit,9002,end=2000) rvert(1),zvert(1),rvert(4),zvert(4)
c      read(gridunit,'(a)',end=2000) buffer
c
c      write(15,9000) in,ik,ir,rvert(1),zvert(1),rvert(2),zvert(2)
c      write(15,9001) brat,rcent,zcent
c      write(15,9002) rvert(4),zvert(4),rvert(3),zvert(3)
c      write(15,'(a)') buffer
c
c     Assign all the values to the appropriate arrays
c     Korpg is a cross-reference to the cell vertices. nvertp
c     contains the number of vertices for a cell ... which is
c     fixed at 4 for the present.
c
c     Increment all of the indices by 1 because DIVIMP does not
c     use zeroth element array addressing for these arrays.
c
c      indexcnt = indexcnt + 1
c
c      if (indexcnt.eq.1) then
c         if (in.eq.0) then
c            indexadj = 1
c         elseif (in.eq.1) then
c            indexadj = 0
c         endif
c      endif
c
c      ik = ik + indexadj
c      ir = ir + indexadj
c      in = in + indexadj
c
c      if (ir.le.cutring) then
c         if (ik.le.cutpt1) then
c            ir = maxrings + ir
c         elseif (ik.ge.cutpt2) then
c
c            if (ik.eq.cutpt2) then
c
c              Add a repeat point to the core
c              plasma rings - to close them.
c
c               ikold = ikold +1
c               rs(ikold,irold) = rs(1,irold)
c               zs(ikold,irold) = zs(1,irold)
c               korpg(ikold,irold) = korpg(1,irold)
c               bratio(ikold,irold) = bratio(1,irold)
c            endif
c
c            ir = maxrings + ir
c            ik = ik - cutpt2 + cutpt1 +1
c
c         else
c            ik = ik - cutpt1
c         endif
c      endif
c
c      korpg(ik,ir) = in
c      nvertp(in) = 4
c      do 210 id = 1,nvertp(in)
c         rvertp(id,in) = rvert(id)
c         zvertp(id,in) = zvert(id)
c 210  continue
c      rs(ik,ir) = rcent
c      zs(ik,ir) = zcent
c      bratio(ik,ir) = brat
c
c     Check to see if a ring end has been passed.
c
c     For the trap rings ... this will be filled out twice
c     first when it hits the cut and again at the end of the
c     ring ... i.e. the two times when ir changes.
c
c
c      if (ir.ne.irold)  nks(irold) = ikold
c
c      ikold = ik
c      irold = ir
c
c     Exit condition ... can be supplemented with the END of file
c     exit in the first read statement at line 200.
c
c      if (ir.eq.maxrings.and.ik.eq.maxkpts) then
c         nks(ir) = maxkpts
c         npolyp = in
c         goto 300
c      endif
c
c     Loop back
c
c      goto 200
c


      do ir = 1,nrs
        do ik = 1,nks(ir)
          in = korpg(ik,ir)
c
c         Check for 4 vertex polygons
c
          if (nvertp(in).ne.4) then
             write(iounit,6009) in,ik,ir,nvertp(in)
          else
c
c           Print indexing info and first two vertices
c
            write(iounit,6002) in,ik,ir,rvertp(1,in),zvertp(1,in),
     >                       rvertp(2,in),zvertp(2,in)
            write (iounit2,7000)  rvertp(1,in),zvertp(1,in),
     >                            rs(ik,ir),zs(ik,ir),ik,ir
            write (iounit2,7001)  rvertp(2,in),zvertp(2,in)
            write (iounit2,7001)  rvertp(3,in),zvertp(3,in)
            write (iounit2,7001)  rvertp(4,in),zvertp(4,in)
            write (iounit2,7001)  rvertp(1,in),zvertp(1,in)

 7000       format(4(g17.10,2x),2(i6,2x))
 7001       format(2(g17.10,2x))
c
c           Print magnetic field and center coordinates
c
            write(iounit,6003) 1.0/kbfs(ik,ir),rs(ik,ir),zs(ik,ir)
c
c           Print second two vertices
c
            write(iounit,6004) rvertp(4,in),zvertp(4,in),
     >                       rvertp(3,in),zvertp(3,in)
c
c           Print cell area information - include rho,theta
c           data for JET grids
c
            if (cgridopt.eq.0) then
              write(iounit,6010) kareas(ik,ir),karea2(ik,ir)
c              write(iounit,6008) kareas(ik,ir),karea2(ik,ir),
c     >          rhog(ik,ir)*hro(ik,ir)*thetag(ik,ir)*hteta(ik,ir)
              write(iounit,6007) rhog(ik,ir),hro(ik,ir),
     >                            thetag(ik,ir),hteta(ik,ir)
            else
              write(iounit,6010) kareas(ik,ir),karea2(ik,ir)
            endif
          endif
          write(iounit,6006)
        end do
      end do

      return
c
C-----------------------------------------------------------------------
c
c     Format statements
c
C-----------------------------------------------------------------------
c
6001  format('Element Output:')
6002  format('Element:',i4,1x,'=',1x,'(',i2,',',i2,'):',1x,
     >       '(',f12.7,',',f12.8,')',3x,'(',f12.8,',',f12.8,')')
6003  format('Field Ratio  = ',e17.10,7x,'(',f12.8,',',f12.8,')')
6004  format(24x,'(',f12.8,',',f12.8,')',3x,'(',f12.8,',',f12.8,')')
6005  format('==========================================================
     >=======================================================')
6006  format('----------------------------------------------------------
     >-------------------------------------------------------')
6007  format('(Rho,Hrho)= ','(',e12.7,',',e12.7,')',
     >       '   (The,Hthe)= ','(',e12.7,',',e12.7,')')
6008  format('Cell Area Data:  KAREAS = ',g17.10,'  KAREA2 = ', g17.10,
     >       ' Rho*Theta =',e17.10)
6009  format('Element:',i4,2x,'=',1x,'(',i4,',',i4,'):',4x,
     >       'ERROR - Not a four-vertex polygon: N = ',i4)
6010  format('Cell Area Data:  KAREAS = ',g17.10,'  KAREA2 = ', g17.10)
6011  format(' ')
c
c     Sonnet Grid format.
c
 9000 format(3x,'Element',i5,' = (',i3,',',i3,'): (',
     >       e17.10e2,',',e17.10e2,')',6x,'(',e17.10e2,',',e17.10e2,')')
 9001 format(3x,'Field ratio  = ',e17.10e2,13x,
     >       '(',e17.10e2,',',e17.10e2,')')
 9002 format(3x,26x,'(',e17.10e2,',',e17.10e2,')',6x,
     >       '(',e17.10,',',e17.10,')')

      end
c
c
c
      real function calcwav(ik,ir,q1,q2)
      implicit none
      integer ik,ir
      real q1,q2
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     This function interpolates the value of theta between
c     two grid points for target option 6.
c
      real    d1,d2,drp1,dzp1,drp2,dzp2,d3
      integer in, findnum
      external findnum
c
      calcwav = 0.0
c
      in = findnum(platco,nplat,ir)
c
      if (ik.lt.nks(ir)/2) then
         drp2 = (rs(ik-1,ir) - platco(in,4))
         dzp2 = (zs(ik-1,ir) - platco(in,5))
         drp1 = (rs(ik,ir) - platco(in,4))
         dzp1 = (zs(ik,ir) - platco(in,5))
         d1 = drp1**2 + dzp1**2
         d2 = drp2**2 + dzp2**2
         d3 = ((drp1-drp2)**2+(dzp1-dzp2)**2)
c
      else
         drp2 = (rs(ik+1,ir) - platco(in,2))
         dzp2 = (zs(ik+1,ir) - platco(in,3))
         drp1 = (rs(ik,ir) - platco(in,2))
         dzp1 = (zs(ik,ir) - platco(in,3))
         d1 = drp1**2 + dzp1**2
         d2 = drp2**2 + dzp2**2
         d3 = ((drp1-drp2)**2+(dzp1-dzp2)**2)
c
      endif
c
      if (d1.eq.0.0.and.d2.eq.0.0) then
         calcwav = q1
      else
         calcwav = ( d2 * q1  + d1 * q2) / (d1+d2)
      endif
c
      return
      end
c
c
c
      subroutine redge2d(flag)
      implicit none
      integer flag
c
c     The purpose of this routine is to read the EDGE2D ghost file
c     and to extract specific quantities from it for use in DIVIMP.
c     In the present instance it is being used only by SOL option
c     22 to extract the EDGE2D values at the mid-point of the first
c     cell so that they can be used in a modified version of
c     SOL option 22 which should be based on the EDGE2D boundary
c     conditions. Initially this routine is based on teh code to read
c     the EDGE2D background forun in rjet - it will only work with
c     JET grids and assumes that teh rjet routine has already been
c     called so that the necessary values have been setup.
c
c     David Elder,  1995 Oct 26
c
c     P.S. The flag value will allow this routine to be used for
c          various purposes -
c          currently:
c          0 - read Edge2D BG and store into specified variables
c          1 - read Edge2D BG and record only first cell center
c              values for use in SOL option 22.
c
c     P.P.S. It is VERY IMPORTANT to make sure whether the virtual
c            points have been stripped from the DIVIMP grids before
c            this routine is called or not - this will affect the IK
c            indexing of the Edge2D arrays here.
c
c

C     INCLUDE "PARAMS"
      include 'params'
C     INCLUDE "CGEOM"
      include 'cgeom'
C     INCLUDE "CNEUT"
c      include 'cneut'
C     INCLUDE "COMTOR"
      include 'comtor'
C     INCLUDE "CIONIZ"
      include 'cioniz'
C     INCLUDE "READER"
      include 'reader'
C     INCLUDE "DYNAM5"
      include 'dynam5'
c
      include 'cadas'
c
c     Variables to hold edge2d values
c
      include 'cedge2d'
c
      CHARACTER MESAGE*72,C(10)*9,FACTOR*9
      INTEGER IK,IR,K,NP,L,J,I,NR,NC,NXW,IEXTRA,JK,JR,MIZS,IZ,IERR,ID
      integer in
      INTEGER IX,IY,IKIN,IKOUT,IRIN,IROUT,MKS,NP1,ICOUNT,INEXT,KNEXT
      INTEGER IDUMMY(2),NJ(MAXNRS),NI(MAXNKS)
      INTEGER ITAG(MAXNKS*MAXNRS,5),ITAGDV(MAXNKS*MAXNRS)
      INTEGER KORX(MAXNKS,MAXNRS),KORYE(MAXNRS,MAXNKS)
      INTEGER ISEP,NINOMP,NINOWA,IREF,NZ,NZMAX,K0,IK0,IKMID
      INTEGER IND,IW
      REAL    DUMMY(MAXNKS*MAXNRS,12+MAXE2DIZS),HMASS,ZMASS,S,DELTAL,
     >             P,FACT
      REAL    sPERP,BEST,DSQ,R,Z,THIN,THOUT,DELTAR,DELTAZ
      REAL    MU,SMAX,SMID,SDIFF
      real    qovrt
C
      integer    jplft,jprgt
      integer    ighost,ifluid,ieshot
      real       xji,xj,xjf,teslic
c
      real teav,tiav,cslocal,ldist,rmid,zmid,cosalpha
      real    tmpne
c
      real e2dpi(maxnks,maxnrs),e2dpe(maxnks,maxnrs)
C
c
C
C-----------------------------------------------------------------------
C     Extract Densities, Temperatures etc
C     NB:  EDGE is often run on a subset of the rings found in the
C          equilibrium file.  For this reason a separate KORY matrix
C          must be read from the EDGE output and used to assign
C          information to specific ring and knot pairs.  Points which
C          have no EDGE information (ie were ignored in the calculation)
C          are mapped from the nearest available neighbour.
C     NB:  Some (but not all) EDGEoutput files have been reformatted
C          to have identifying labels.  These labels are not yet
C          used as keys, thus making the order of the file unimportant,
C          so we will just skip them for now.
C-----------------------------------------------------------------------
C
c
c     Initialize cre2dizs
c
      cre2dizs = -1
c
      REWIND (11)
c
c     Check the IGHOST parameter to see which version of the GHOST file
c     is in use.
c
      ighost = 0
c
      read (11,'(a72)',iostat=ierr) buffer
c
      if (buffer(1:12).eq.'       HMASS')
     >    READ (buffer,'(71X,I1)') IGHOST
c
      write (6,*) 'REDGE2D:Ighost:',ighost
      write (6,*) 'Buff0 :',buffer
c
      if (ighost.ne.0) then
c
         write (6,*) 'Reading NEW Ghost file'
         rewind(11)
c
c        Read in an indexed GHOST file.
c
 201     read (11,'(a72)',end=203,err=202,iostat=ierr) buffer
c
         if (buffer(1:12).eq.'       HMASS') then
            READ (11,*) HMASS
         endif
c
         if (buffer(1:7).eq.'  ISHOT') then
            READ (11,9044) IESHOT,TESLIC,CRUN
            write(6,*) 'Buff1 :',buffer
            write (6,*) 'ishot:',IESHOT,TESLIC,CRUN
c
c
c           Read edge shot and time into separate variables
c           and check against
c           values from equilibrium file.  Should really
c           be stopping the code
c           if they don't match!
c
            IF (ISHOT.NE.IESHOT .OR. ABS(TSLICE-TESLIC).GT.1.0E-3) THEN
             WRITE (6,*) ' Mismatched EDGE and EQUILIBRIUM shot/times!'
             WRITE (6,*) ' EQUILIBRIUM: ISHOT = ',ISHOT, ' TIME = ',
     >                     TSLICE
             WRITE (6,*) ' EDGE:        ISHOT = ',IESHOT,' TIME = ',
     >                     TESLIC
            ENDIF
         endif
c
         if (buffer(1:7).eq.'    NP') then
            READ (11,*) NP1,IFLUID,NZ
            write (6,*) 'Buff2:',buffer
            write (6,*) 'np   :',NP1,IFLUID,NZ
c
            IF (NP1.GT.MAXNKS*MAXNRS)
     >         WRITE (6,*) ' Error! EDGE np=',NP1
            NZMAX = NZ
            IF (NZ.LE.0) NZMAX = 1
         endif
c
         if (buffer(1:5).eq.' ITAG') then
            READ (11,'(12I6)') ((ITAG(K,L),L=1,5),K=1,NP1)
         endif
c
         if (buffer(1:3).eq.' NR') then
            READ (11,*) NR
            write (6,*) 'Buff3:',buffer
            write (6,*) 'nr   :',NR
            IF (NR.GT.MAXNKS) WRITE (6,*) ' Error! EDGE nr=',NR
         endif
c
         if (buffer(1:3).eq.' NI') then
            READ (11,*) (NI(J),J=1,NR)
            DO 212 J = 1, NR
               IF (NI(J).LE.0) GOTO 212
               READ (11,'(A72)') BUFFER
               READ (11,*) (KORX(J,I),I=1,NI(J))
 212        CONTINUE
         endif
c
         if (buffer(1:7).eq.'    NC ') then
            READ (11,*) NC,NXW,IEXTRA,ZMASS
            write (6,*) 'Buff4:',buffer
            write (6,*) 'nc   :',NC,NXW,IEXTRA,ZMASS
            IF (NC+IEXTRA.GT.MAXNRS) WRITE (6,*) ' Error! EDGE nc=',NC
         endif
C
         if (buffer(1:3).eq.' NJ'.and.buffer(1:4).ne.' NJP') then
            READ (11,*) (NJ(I),I=1,NC+IEXTRA)
            DO 214 I = 1, NC+IEXTRA
               READ (11,'(A72)') BUFFER
               READ (11,*) (KORYE(I,J),J=1,NJ(I))
 214        CONTINUE
         endif
c
         if (buffer(1:6).eq.' DEN (')
     >      READ (11,9012) (DUMMY(K,1),K=1,NP1)
c
         if (buffer(1:6).eq.' VTE (')
     >      READ (11,9012) (DUMMY(K,2),K=1,NP1)
c
         if (buffer(1:6).eq.' TEV (')
     >      READ (11,9012) (DUMMY(K,3),K=1,NP1)
c
         if (buffer(1:6).eq.' DA (A')
     >      READ (11,9012) (DUMMY(K,10),K=1,NP1)
c
         if (buffer(1:6).eq.' VRO (')  then
c     >      READ (11,9012) (DUMMY(K,11),K=1,NP1)
             READ (11,9012) (DUMMY(K,11),K=1,NP1)
             write (6,*) 'READING VRO:'
         endif
c
         if (buffer(1:6).eq.' TEVE ')
     >      READ (11,9012) (DUMMY(K,4),K=1,NP1)
c
         if (buffer(1:6).eq.' SOUN ')
     >      READ (11,9012) (DUMMY(K,6),K=1,NP1)
c
         if (buffer(1:6).eq.' VPOT ')
     >      READ (11,9012) (DUMMY(K,5),K=1,NP1)
c
         if (buffer(1:6).eq.' PMACH')
     >      READ (11,9012) (DUMMY(K,7),K=1,NP1)
c
         if (buffer(1:6).eq.' DENPI')
     >      READ (11,9012) (DUMMY(K,8),K=1,NP1)
c
         if (buffer(1:6).eq.' DENPE')
     >      READ (11,9012) (DUMMY(K,9),K=1,NP1)
c
c        Neutral Impurity density
c
         if (buffer(1:6).eq.' DZ (N') then
            READ (11,9012) (DUMMY(K,12),K=1,NP1)
            cre2dizs = 0
            write (6,*) 'DZ:',cre2dizs
         endif
c
c
c        Impurity density from EDGE2D
c
         if (buffer(1:7).eq.' DENZ :') then
c
            read(buffer(12:15),'(i4)') iz
c
            write (6,*) 'denz:',iz,'|',buffer(12:15),'|',
     >                        cre2dizs
c
            if (iz.ge.1.and.iz.le.maxe2dizs) then
               READ (11,9012) (DUMMY(K,12+iz),K=1,NP1)
               cre2dizs = max(cre2dizs,iz)
               write (6,*) 'DENZ:',iz,'|',buffer(12:15),'|',
     >                             cre2dizs
            endif
c
         endif
c
         goto 201

202      write(6,*) 'INPUT Error Reading Edge2d Plasma:'
         write(6,*) 'Error Code = ',ierr

203      continue

         IF (IFLUID.ne.2) THEN
            DO K=1,NP1
               DUMMY(K,4) = DUMMY(K,3)
            ENDDO
         ENDIF
c
c
      else
c
c        Read in an old GHOST file without indices
c
c        This just uses the existing code so there are
c        a number of redundant if (ighost.ne.1) statements
c
C-----------------------------------------------------------------------
      write (6,*) 'Reading OLD Ghost file'
      REWIND (11)
      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) HMASS

c  Read edge shot and time into separate variables and check against
c  values from equilibrium file.  Should really be stopping the code
c  if they don't match!

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,9044) IESHOT,TESLIC,CRUN
      IF (ISHOT.NE.IESHOT .OR. ABS(TSLICE-TESLIC).GT.1.0E-3) THEN
        WRITE (6,*) ' Mismatched EDGE and EQUILIBRIUM shot/times!'
        WRITE (6,*) ' EQUILIBRIUM: ISHOT = ',ISHOT, ' TIME = ',TSLICE
        WRITE (6,*) ' EDGE:        ISHOT = ',IESHOT,' TIME = ',TESLIC
      ENDIF


      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) NP1,IFLUID,NZ
      IF (NP1.GT.MAXNKS*MAXNRS) WRITE (6,*) ' Error! EDGE np=',NP1
      NZMAX = NZ
      IF (NZ.LE.0) NZMAX = 1

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) (DUMMY(K,1),K=1,NP1)
      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) (DUMMY(K,1),K=1,NP1)

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,'(12I6)') ((ITAG(K,L),L=1,5),K=1,NP1)


      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) NR
      IF (NR.GT.MAXNKS) WRITE (6,*) ' Error! EDGE nr=',NR

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) (NI(J),J=1,NR)

      DO 210 J = 1, NR
        IF (NI(J).LE.0) GOTO 210
        IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
        READ (11,*) (KORX(J,I),I=1,NI(J))
  210 CONTINUE

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) NC,NXW,IEXTRA,ZMASS
      IF (NC+IEXTRA.GT.MAXNRS) WRITE (6,*) ' Error! EDGE nc=',NC

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,*) (NJ(I),I=1,NC+IEXTRA)
      DO 220 I = 1, NC+IEXTRA
        IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
        READ (11,*) (KORYE(I,J),J=1,NJ(I))
  220 CONTINUE
C
C     DEN(K)                 -----  ION DENSITY    (CM-3)
C     TEV(K)                 -----  ION TEMP       (EV)
C     TEVE(K)                -----  ELE TEMP       (EV)
C     VTE(K)                 -----  PARALLEL SPEED (CM/SEC)
C     EFIELD(K)              -----  ELECTRIC FIELD (V/CM) as only given
C                                   for average ring then assign this to
C                                   all SOL rings.
C
      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,9012) (DUMMY(K,1),K=1,NP1)

      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,9012) (DUMMY(K,2),K=1,NP1)
      IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
      READ (11,9012) (DUMMY(K,3),K=1,NP1)

      DO 230 I = 1, 8
        IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
        READ (11,9012) (DUMMY(K,4),K=1,NP1)
  230 CONTINUE
      IF (IFLUID.EQ.2) THEN
        IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
        READ (11,9012) (DUMMY(K,4),K=1,NP1)
        DO I = 1, 3
          IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
          READ (11,9012) (DUMMY(K,5),K=1,NP1)
        ENDDO
      ELSE
        DO K=1,NP1
          DUMMY(K,4) = DUMMY(K,3)
        ENDDO
      ENDIF
      DO 235 I = 1, (NZMAX*3) + 8
        IF (IGHOST.NE.0) READ (11,'(A72)') BUFFER
        READ (11,9012) (DUMMY(K,5),K=1,NP1)
  235 CONTINUE
c
c     End of if statement for old GHOST file compatibility
c
      endif

c
c
c
C
C  FIND SEPARATRIX RING IN EDGE GRID AND CALCULATE HOW MANY
C  OF THE GRID2D RINGS WERE NOT USED IN THE PLASMA CORE (NINOMP)
C  AND THE SOL (NINOWA)
C
      ISEP  = 0
      DO 200 K = 1, NP1
C  FIRST POINT NEXT TO EITHER TARGET IS ON THE FIRST SOL RING (ISEP)
        IF (ISEP.EQ.0  .AND. (ITAG(K,4).EQ.4 .OR. ITAG(K,4).EQ.5) )
     >                                        ISEP = ITAG(K,1)
        IF (ISEP.NE.0) GOTO 205
  200 CONTINUE
C
  205 NINOMP = IRSEP - ISEP
      NINOWA = (NRS - NC) - NINOMP
C
      IRCENT = NINOMP + 1
c
      write (6,*) 'info:',irsep,isep,ninomp,nrs,nc,ninowa,ircent
C
C  VIRTUAL POINTS:  THE INNERMOST AND OUTERMOST RINGS AND THE END
C                   POINTS OF THE OPEN RINGS ARE VIRTUAL IN EDGE,
C                   THAT IS THEY ARE USED ONLY IN THE BOUNDARY
C                   CONDITION SO THAT THE SOLUTION AT THESE POINTS
C                   IS NOT MEANINGFUL.  ** IN THE EDGE CODES THESE
C                   POINTS WERE OUTSIDE THE PLASMA **  FOR DIVIMP
C                   WE'LL MAP THE NEAREST PLASMA POINT ONTO THE
C                   VIRTUAL POINTS IN ORDER TO GENERATE A BIGGER
C                   PLASMA
C
C  START WITH THE THREE VIRTUAL RINGS (1,NXW,NXW+1)
C
      DO 260 L             = 1 , 3
C
         IF     ( L.EQ.1 ) THEN
                  I        = 1
         ELSE IF( L.EQ.2 ) THEN
                  I        = NXW
         ELSE
                  I        = NXW + 1
         END IF
C
         ICOUNT            = 0
         DO 240 J          = 1 , NJ(I)
            K              = KORYE(I,J)
            IF(I TAG(K,4).LT.0 ) ICOUNT = ICOUNT + 1
  240    CONTINUE
C
         IF( ICOUNT.EQ.NJ(I) ) THEN
             IF( I.EQ.1 .OR. I.EQ.NXW+1 ) THEN
                 INEXT     = I + 1
             ELSE
                 INEXT     = I - 1
             END IF
             DO 250 J      = 1 , NJ(I)
                K          = KORYE(I,J)
                KNEXT      = KORYE(INEXT,J)
                DUMMY(K,1) = DUMMY(KNEXT,1)
                DUMMY(K,2) = DUMMY(KNEXT,2)
                DUMMY(K,3) = DUMMY(KNEXT,3)
                DUMMY(K,4) = DUMMY(KNEXT,4)
                DUMMY(K,5) = DUMMY(KNEXT,5)
                DUMMY(K,6) = DUMMY(KNEXT,6)
                DUMMY(K,7) = DUMMY(KNEXT,7)
                DUMMY(K,8) = DUMMY(KNEXT,8)
                DUMMY(K,9) = DUMMY(KNEXT,9)
                DUMMY(K,10)= DUMMY(KNEXT,10)
                DUMMY(K,11)= DUMMY(KNEXT,11)
                do iz = 0,maxe2dizs
                   DUMMY(K,12+iz)= DUMMY(KNEXT,12+iz)
                end do
c
  250        CONTINUE
         END IF
C
  260 CONTINUE
C
c  Note: The following is not necessary since Edge2d already has plasma
c  mapped onto the virtual ring end-points - and these values have
c  some significance.
c
c  David Elder, Oct 18, 1995
c
C  NOW MAP PLASMA ONTO THE END POINTS OF OPEN RINGS
C  (PTS 1 AND NJ(I) OF RINGS I= ISEP, NC+IEXTRA)
C
c      DO 261 I = ISEP, NC+IEXTRA
c         DO 262 L = 1,6
c            DUMMY(KORYE(I,1),L)     = DUMMY(KORYE(I,2),L)
c            DUMMY(KORYE(I,NJ(I)),L) = DUMMY(KORYE(I,NJ(I)-1),L)
c  262    CONTINUE
c  261 CONTINUE
C
C.......................................................................
C
      DO 265 I = 1 , CNVMF
C
         IF     ( CIRNG0(I).EQ.-10 ) THEN
                  CIRNG0(I) = IRSEP - 1
         ELSE IF( CIRNG0(I).EQ.-20 ) THEN
                  CIRNG0(I) = IRSEP
         ELSE IF( CIRNG0(I).EQ.-30 ) THEN
                  CIRNG0(I) = IRWALL
         ELSE IF( CIRNG0(I).EQ.-40 ) THEN
                  CIRNG0(I) = IRWALL + 1
         ELSE IF( CIRNG0(I).EQ.-50 ) THEN
                  CIRNG0(I) = NRS
         END IF
         IF( CIRNG0(I).LT.1   ) CIRNG0(I) = 1
         IF( CIRNG0(I).GT.NRS ) CIRNG0(I) = NRS
C
         IF     ( CIRNG1(I).EQ.-10 ) THEN
                  CIRNG1(I) = IRSEP - 1
         ELSE IF( CIRNG1(I).EQ.-20 ) THEN
                  CIRNG1(I) = IRSEP
         ELSE IF( CIRNG1(I).EQ.-30 ) THEN
                  CIRNG1(I) = IRWALL
         ELSE IF( CIRNG1(I).EQ.-40 ) THEN
                  CIRNG1(I) = IRWALL + 1
         ELSE IF( CIRNG1(I).EQ.-50 ) THEN
                  CIRNG1(I) = NRS
         END IF
         IF( CIRNG1(I).LT.1   ) CIRNG1(I) = 1
         IF( CIRNG1(I).GT.NRS ) CIRNG1(I) = NRS
C
         IF( CIRNG0(I).GT.CIRNG1(I) ) THEN
             WRITE (7,'(/2(A,I2,A,I3))')
     >          ' *** ERROR *** CIRNG0(' , I , ') = ' , CIRNG0(I) ,
     >          '  .GT.  CIRNG1(' , I , ') = ' , CIRNG1(I)
             STOP
         END IF
C
  265 CONTINUE
c
c     FLAG = 0
c
c     Regular situation - read Edge2D quantites
c
c
      if (flag.eq.0) then

C
C  MAP EDGE PLASMA ONTO FULL GRID
C
      DO 280 IR = 1 , NRS
C
         IREF   = IR - IRCENT + 1
C  USE INNERMOST EDGE RING FOR MISSING INNER RINGS
         IF( IR.LT.IRCENT )   IREF = 1
         IF( IR.GT.IRWALL-NINOWA .AND. IR.LE.IRWALL ) THEN
             IF( CNIWA.EQ.0 ) IREF = 0
C  USE OUTERMOST EDGE SOL RING FOR MISSING SOL RINGS
             IF( CNIWA.EQ.1 ) IREF = NXW
         END IF
         IF( IR.GT.IRWALL )   IREF = IREF - NINOWA
c
c      write (0,*) 'IR:',ir,iref,ircent,irwall,ninowa,cniwa,nks(ir)
c     >                       ,nj(ir)
C
         DO 270 IK = 1 , nks(ir)
            CVMF(IK,IR) = 1.0E+00
            IF( IREF.NE.0 ) THEN
                DO 275 I = 1 , CNVMF
                   IF( IR.GE.CIRNG0(I) .AND. IR.LE.CIRNG1(I) ) THEN
                       IF     ( IK.LE.CJ0(I) ) THEN
                                XJI         = 1.00E+00
                                XJ          = IK
                                XJF         = CJ0(I)
                                CVMF(IK,IR) = CVMF0(I) +
     >                                        (CVMF1(I)-CVMF0(I)) *
     >                                        (XJ-XJI)/(XJF-XJI)
                       ELSE IF( IK.GE.nks(IR)-CJ1(I)+1 ) THEN
                                XJI         = nks(IR)-CJ1(i)+1
                                XJ          = IK
                                XJF         = nks(IR)
                                CVMF(IK,IR) = CVMF1(I) +
     >                                        (CVMF2(I)-CVMF1(i)) *
     >                                        (XJ-XJI)/(XJF-XJI)
                       ELSE
                                CVMF(IK,IR) = CVMF1(I)
                       END IF
                       GOTO 276
                   END IF
  275           CONTINUE
C
  276           K  = KORYE(IREF,IK)
c
c               if (cprint.eq.4.or.cprint.eq.9) then
c                  WRITE(0,1) K , IK , IR, iref
c   1              FORMAT('REDGE2D: K=',I5,' IK=',I5,
c     >                    ' IR=',I5,
c     >                    ' IREF=',i5)
c               endif
c
                E2DNBS (IK,IR) = MAX (DUMMY(K,1),1.0) * 1.0E6
c
                do iz = 0,cre2dizs
                   e2dnzs(ik,ir,iz) = dummy(k,12+iz) * 1.0e6
                end do
c
c               Copy neutral content to separate array for compatibility
c               elsewhere.
c
                e2dz0(ik,ir) = e2dnzs(ik,ir,0)
c
c               Normalization of Velocity by Quantum time step
c               is in the main TAU routine.
c
                E2DVHS (IK,IR) = DUMMY(K,2) * CVMF(IK,IR) * 0.01
                e2dvro(ik,ir)  = dummy(k,11) * 0.01
c
                E2DTIBS(IK,IR) = MAX (DUMMY(K,3), LO)
                E2DTEBS(IK,IR) = MAX (DUMMY(K,4), LO)
                e2dmach(ik-1,ir) = dummy(k,7)
c
c               Only assign these following quantities for rings with an
c               actual EDGE2D solution - otherwise zero.
c
                if (ir.lt.ircent.or.
     >             (ir.GT.IRWALL-NINOWA .AND. IR.LE.IRWALL )) THEN
                   e2dpi(ik,ir)   = 0.0
                   e2dpe(ik,ir)   = 0.0
                   E2DION(ik,ir)  = 0.0
                   E2DATOM(ik,ir) = 0.0
                else
                   e2dpi(ik,ir)   = dummy(k,8) * 1.0e10
                   e2dpe(ik,ir)   = dummy(k,9) * 1.0e10
                   E2DION(ik,ir)  = dummy(k,6) * 1.0e6
                   E2DATOM(ik,ir) = dummy(k,10) * 1.0e6
                endif
c
                IF( IR.GE.IRSEP .AND. IR.LE.IRWALL ) THEN
                    IF( IEXTRA.GT.0 ) THEN
                        K0   = KORYE(NC+IEXTRA,IK)
                    ELSE
                        K0   = K
                    END IF
                    E2DES(IK,IR) = DUMMY(K0,5) * 100.0 * QTIM*QTIM
     >                           * EMI / CRMI
                ELSE
                    E2DES(IK,IR) = 0.0
                END IF
            ELSE
                E2DNBS (IK,IR) = 1.0E6
                E2DVHS (IK,IR) = 0.0
                E2DTIBS(IK,IR) = LO
                E2DTEBS(IK,IR) = LO
                E2DES  (IK,IR) = 0.0
                e2dpi(ik,ir)   = 0.0
                e2dpe(ik,ir)   = 0.0
                E2DION (ik,ir) = 0.0
                e2dvro(ik,ir)  = 0.0
                e2datom (ik,ir) = 0.0
                e2dmach(ik-1,ir) = 0.0
            END IF
C

  270    CONTINUE
C
  280 CONTINUE

c
c     Calculate electron density
c
      if (cre2dizs.gt.1) then

         do ir = 1,nrs
            do ik = 1,nks(ir)

               tmpne = e2dnbs(ik,ir)

               do iz = 1,cre2dizs

                  tmpne = tmpne + iz * e2dnzs(ik,ir,iz)

               end do

               if (tmpne.gt.1.1*e2dnbs(ik,ir)) then
                  write (6,*) '***NOTE***'

                  write(6,'(a,2i4,4(1x,g12.5))')
     >               'IONIZ:',ik,ir,e2dnbs(ik,ir),tmpne,
     >               e2dnzs(ik,ir,1),e2dnzs(ik,ir,2)
               endif
c
c              Assign to fluid code array
c
               e2dnes(ik,ir) = tmpne
c
            end do

         end do

      endif

C
c
c       Calculate the Edge2D target conditions
c       This must be done before the virtual points are
c       stripped away.
c
c       The values in E2Dtarg are:
c       1 - Density
c       2 - Te
c       3 - Ti
c       4 - Vb
c       5 - Gamma = Ne * Vb
c
c       Note: Vb is calculated as the sound speed for the given target
c             temperatures times the E2D reported mach number.
c
c       NOTE: NKS(IR) has NOT yet been modified for the removal of
c             any virtual points (when flag=0) and thus nks(ir)
c             will point to the virtual point.
c
c
        write(6,*) 'E2Dtarg:'
c
        do ir = irsep,nrs
c
c        Target 1 - this is  the INNER target forJET grids
c                   or grids with the X-POINT at the top
c        Target 2 - this is  the OUTER target for JET grids
c                   or grids with the X-POINT at the top
c
c        The indexing of the e2dtarg array is:
c        e2dtarg(ring#, quantity index#, target#) 
c
c
c        Ne 
c
         if (fc_ne_calc_opt.eq.0) then   
            e2dtarg(ir,1,2) = e2dnbs(1,ir)
            e2dtarg(ir,1,1) = e2dnbs(nks(ir),ir)
         elseif (fc_ne_calc_opt.eq.1) then  
            e2dtarg(ir,1,2) = (e2dnbs(1,ir)+e2dnbs(2,ir))/2.0
            e2dtarg(ir,1,1) = (e2dnbs(nks(ir),ir)
     >                        +e2dnbs(nks(ir)-1,ir))/2.0
         elseif (fc_ne_calc_opt.eq.2) then  
            e2dtarg(ir,1,2) = e2dnbs(2,ir)
            e2dtarg(ir,1,1) = e2dnbs(nks(ir)-1,ir)
         endif    
c
c        Te
c
         if (fc_te_calc_opt.eq.0) then   
            e2dtarg(ir,2,2) = e2dtebs(1,ir)
            e2dtarg(ir,2,1) = e2dtebs(nks(ir),ir)
         elseif (fc_te_calc_opt.eq.1) then  
            e2dtarg(ir,2,2) = (e2dtebs(1,ir)+e2dtebs(2,ir))/2.0
            e2dtarg(ir,2,1) = (e2dtebs(nks(ir),ir)
     >                       +e2dtebs(nks(ir)-1,ir))/2.0
         elseif (fc_te_calc_opt.eq.2) then  
            e2dtarg(ir,2,2) = e2dtebs(2,ir)
            e2dtarg(ir,2,1) = e2dtebs(nks(ir)-1,ir)
         endif    
c
c        Ti
c
         if (fc_ti_calc_opt.eq.0) then   
            e2dtarg(ir,3,2) = e2dtibs(1,ir)
            e2dtarg(ir,3,1) = e2dtibs(nks(ir),ir)
         elseif (fc_ti_calc_opt.eq.1) then  
            e2dtarg(ir,3,2) = (e2dtibs(1,ir)+e2dtibs(2,ir))/2.0
            e2dtarg(ir,3,1) = (e2dtibs(nks(ir),ir)
     >                       +e2dtibs(nks(ir)-1,ir))/2.0
         elseif (fc_ti_calc_opt.eq.2) then  
            e2dtarg(ir,3,2) = e2dtibs(2,ir)
            e2dtarg(ir,3,1) = e2dtibs(nks(ir)-1,ir)
         endif    
c
c        Vb
c
         if (fc_v_calc_opt.eq.0) then  
c
            e2dtarg(ir,4,2) = e2dmach(1,ir) *
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
     >               *(1.0+rizb)/crmb)
c
            e2dtarg(ir,4,1) = e2dmach(nks(ir)-1,ir) *
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
     >                    *  (1.0+RIZB)/crmb)
c
c        Vb = Cs 
c
         elseif (fc_v_calc_opt.eq.1) then  

             e2dtarg(ir,4,2) = 
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
     >               *(1.0+rizb)/crmb)

             e2dtarg(ir,4,1) = 
     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
     >               *(1.0+rizb)/crmb)

         endif
c
c        Base flux value from Ne * Vb
c
         e2dtarg(ir,5,2) = e2dtarg(ir,1,2) * e2dtarg(ir,4,2)
         e2dtarg(ir,5,1) = e2dtarg(ir,1,1) * e2dtarg(ir,4,1)

c
c ** Leave old code for now - will be deleted later
c
c          Outer target
c
c           e2dtarg(ir,1,2) = e2dnbs(2,ir)
c           e2dtarg(ir,2,2) = e2dtebs(2,ir)
cc           e2dtarg(ir,2,2) = (e2dtebs(1,ir)+e2dtebs(2,ir))/2.0
c           e2dtarg(ir,3,2) = e2dtibs(2,ir)
c           e2dtarg(ir,4,2) = e2dmach(1,ir) *
c     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,2) + e2dtarg(ir,2,2))
c     >               *(1.0+rizb)/crmb)
c           e2dtarg(ir,5,2) = e2dtarg(ir,1,2) * e2dtarg(ir,4,2)
c
c          Inner target
c
c           e2dtarg(ir,1,1) = e2dnbs(nks(ir)-1,ir)
c           e2dtarg(ir,2,1) = e2dtebs(nks(ir)-1,ir)
cc           e2dtarg(ir,2,1) = (e2dtebs(nks(ir)-1,ir)
cc     >                       + e2dtebs(nks(ir),ir))/2.0
c           e2dtarg(ir,3,1) = e2dtibs(nks(ir)-1,ir)
c           e2dtarg(ir,4,1) = e2dmach(nks(ir)-1,ir) *
c     >           9.79E3 * SQRT (0.5*(e2dtarg(ir,3,1) + e2dtarg(ir,2,1))
c     >                    *  (1.0+RIZB)/crmb)
c           e2dtarg(ir,5,1) = e2dtarg(ir,1,1) * e2dtarg(ir,4,1)
c
c
c          Write out
c
  101      format(a4,i4,2e13.5)
c
c
           if (cprint.eq.4.or.cprint.eq.9) then
c             write(6,*)   'NKS:',nks(ir),nj(ir)
             write(6,101) 'Ne:',ir,e2dtarg(ir,1,1),e2dtarg(ir,1,2)
             write(6,101) 'Te:',ir,e2dtarg(ir,2,1),e2dtarg(ir,2,2)
             write(6,101) 'Ti:',ir,e2dtarg(ir,3,1),e2dtarg(ir,3,2)
             write(6,101) 'Vb:',ir,e2dtarg(ir,4,1),e2dtarg(ir,4,2)
             write(6,101) 'Ga:',ir,e2dtarg(ir,5,1),e2dtarg(ir,5,2)
           endif

        end do
c
c       Set up the e2dbvel based on the mach number array
c
c       The arrays e2dmach and e2dbvel are NOT adjusted for virtual
c       grids because they represent cell boundary values.
c
        call rzero (e2dbvel,(maxnks+1)*maxnrs)
        call rzero (e2dflux,(maxnks+1)*maxnrs)
c
        write(6,*) 'E2D Boundary velocities'
c
        do ir = 1,nrs

          if (ir.ne.irwall.and.ir.ne.irtrap) then

            do ik = 1,nks(ir)-1

              teav = 0.5*(e2dtebs(ik,ir)+e2dtebs(ik+1,ir))
c
c             Deal with ion temperature at targets slightly
c             differently
c
              if (ik.eq.1) then
                 tiav = e2dtibs(2,ir)
              elseif (ik.eq.nks(ir)-1) then
                 tiav = e2dtibs(nks(ir)-1,ir)
              else
                 tiav = 0.5*(e2dtibs(ik,ir)+e2dtibs(ik+1,ir))
              endif
c
c             Calculate sound speed
c
              cslocal = 9.79e3 * sqrt((teav+tiav)/crmb)

              e2dbvel(ik,ir) = e2dmach(ik,ir)*cslocal

              e2dflux(ik,ir) = e2dbvel(ik,ir)*e2dnbs(ik+1,ir)

              if (cprint.eq.4.or.cprint.eq.9) then
                write (6,'(''mach:'',i4,i4,6(1x,e14.6))') ir,ik,
     >             e2dmach(ik,ir),e2dbvel(ik,ir),e2dflux(ik,ir),
     >             e2dnbs(ik,ir),e2dvhs(ik,ir),cslocal
              endif
c
c
c             Calculate the separation L(Btheta/B) for each cell
c             boundary. This needs the length of the side AND
c             the angle of impact of the side with the line drawn to
c             the cell centre.
c
c             Since the grid has not yet been modified to remove the
c             virtual cells - this is OK.
c
c              in = korpg(ik,ir)
c
c              write (6,*) 'grid info:',ik,ir,in
c
c             Take polygon corners 3,4
c
c              ldist = sqrt( (rvertp(3,in)-rvertp(4,in))**2 +
c     >                      (zvertp(3,in)-zvertp(4,in))**2)
c
c              rmid = 0.5*(rvertp(3,in)+rvertp(4,in))
c              zmid = 0.5*(zvertp(3,in)+zvertp(4,in))
c
c              write (6,'(''siz1:'',i4,i4,7(1x,e14.6))') ir,ik,
c     >             hteta(ik,ir),thetag(ik,ir),
c     >             hteta(ik,ir)*thetag(ik,ir),
c     >             hro(ik,ir),rhog(ik,ir),hro(ik,ir)*rhog(ik,ir),
c     >             hteta(ik,ir)*thetag(ik,ir)
c     >             *hro(ik,ir)*rhog(ik,ir)
c
c              if (ik.ge.2) then
c                 dthetag(ik,ir) = 0.5*(thetag(ik+1,ir)
c     >                            -thetag(ik-1,ir))
c                 if (ir.eq.irsep) then
c                   drho(ik,ir) =  (rhog(ik,ir+1)-rhog(ik,ir))
c                 elseif (ir.eq.irwall-1) then
c                   drho(ik,ir) =  (rhog(ik,ir)-rhog(ik,ir-1))
c                 elseif (ir.eq.irtrap+1) then
c                   drho(ik,ir) =  (rhog(ik,ir)-rhog(ik,ir+1))
c                 elseif (ir.eq.nrs) then
c                   drho(ik,ir) =  (rhog(ik,ir-1)-rhog(ik,ir))
c                 else
c                   drho (ik,ir) = 0.5* (rhog(ik,ir+1)
c     >                            -rhog(ik,ir-1))
c                 endif
c              endif
c
c              write (6,'(''siz2:'',i4,i4,7(1x,e14.6))') ir,ik,
c     >             hteta(ik,ir),dthetag(ik,ir),
c     >             hteta(ik,ir)*dthetag(ik,ir),
c     >             hro(ik,ir),drho(ik,ir),hro(ik,ir)*drho(ik,ir),
c     >             hteta(ik,ir)*dthetag(ik,ir)
c     >             *hro(ik,ir)*drho(ik,ir)
c
            end do

          endif

        end do
c
c        do ir = irsep,nrs
c
c           if (ir.eq.irwall.or.ir.eq.irtrap) cycle
c
c           if (ir.ne.irwall.and.ir.ne.irtrap) then
c
c                Write out the various velocities
c
c                 write(6,'(''E2D:vels:Outer:'',i4,i4,
c     >                    4(1x,g13.6))')
c     >                 ir,ik, -sqrt((2.0*e2dtebs(1,ir)
c     >                         -e2dtebs(2,ir)
c     >                         +e2dtibs(2,ir))*emi/crmb),
c     >          (e2dbvel(1,ir) + e2dbvel(2,ir))/2.0,
c     >          (e2dflux(1,ir) + e2dflux(2,ir))/2.0/e2dnbs(2,ir),
c     >           e2dvhs(2,ir)
c
c                 write(6,'(''E2D:vels:Inner:'',i4,i4,
c     >                    4(1x,g13.6))')
c     >                 ir,ik, -sqrt((2.0*e2dtebs(nks(ir),ir)
c     >                         -e2dtebs(nks(ir)-1,ir)
c     >                         +e2dtibs(nks(ir)-1,ir))*emi/crmb),
c     >          (e2dbvel(nks(ir)-1,ir) + e2dbvel(nks(ir)-2,ir))/2.0,
c     >          (e2dflux(nks(ir)-1,ir) + e2dflux(nks(ir)-2,ir))/2.0
c     >                    /e2dnbs(nks(ir)-1,ir),
c     >           e2dvhs(nks(ir)-1,ir)
c
c           endif
c
c        end do
c
        write(6,*) 'Power Densities Targets: Edge2D'
        write(6,'(2x,''ir'',11x,''Ion'',6x,''Electron'')')
c
        do ir = irsep,nrs
c
c         write out the power densities at each target.
c
          write(6,'(i4,2(1x,g13.6),a)') ir,e2dpi(2,ir),e2dpe(2,ir),
     >         'Outer'

          write(6,'(i4,2(1x,g13.6),a)') ir,e2dpi(nks(ir)-1,ir),
     >                e2dpe(nks(ir)-1,ir),'Inner'
c
        end do
c
c       Calculate estimate of E2D cell areas
c
        do ir = 1,nrs
c slmod begin
c...      ir-1 below gives out of bounds, temp fix...
          IF (ir.EQ.1) CYCLE
c slmod end
c
           do ik = 1,nks(ir)
c
              if (ik.gt.1.and.ik.lt.nks(ir)) then
c
                 if (ik.eq.2) then
                    dthetag(ik,ir) = (thetag(ik+1,ir)
     >                            -thetag(ik,ir))
                 elseif (ik.eq.nks(ir)-1) then
                    dthetag(ik,ir) = (thetag(ik,ir)
     >                            -thetag(ik-1,ir))
                 else
                    dthetag(ik,ir) = 0.5* (thetag(ik+1,ir)
     >                            -thetag(ik-1,ir))
                 endif
c
                 if (ir.eq.irsep) then
                    drho(ik,ir) =  (rhog(ik,ir+1)-rhog(ik,ir))
                 elseif (ir.eq.irwall-1) then
                    drho(ik,ir) =  (rhog(ik,ir)-rhog(ik,ir-1))
                 elseif (ir.eq.irtrap+1) then
                    drho(ik,ir) =  (rhog(ik,ir)-rhog(ik,ir+1))
                 elseif (ir.eq.nrs) then
                    drho(ik,ir) =  (rhog(ik,ir-1)-rhog(ik,ir))
                 elseif (ir.le.irwall) then
                    drho (ik,ir) = 0.5* (rhog(ik,ir+1)
     >                            -rhog(ik,ir-1))
                 elseif (ir.gt.irwall) then
                    drho (ik,ir) = 0.5* (rhog(ik,ir-1)
     >                            -rhog(ik,ir+1))
                 endif
c
                 e2dareas(ik,ir) = abs(hro(ik,ir)*hteta(ik,ir)
     >                            *drho(ik,ir)*dthetag(ik,ir))

              endif
c

c
c              write (6,'(''siz2:'',i4,i4,7(1x,e14.6))') ir,ik,
c     >             hteta(ik,ir),dthetag(ik,ir),
c     >             hteta(ik,ir)*dthetag(ik,ir),
c     >             hro(ik,ir),drho(ik,ir),hro(ik,ir)*drho(ik,ir),
c     >             e2dareas(ik,ir)
c
           end do
        end do
c

c
c       Read Auxilliary input file if readaux is specified.
c
c
        if (cprint.eq.9) write(6,*) 'READ-e2d-AUX:',readaux
c
        if (readaux.eq.1) then
c
           call reade2daux
c
        endif
c
c        if (ctargopt.eq.0.or.ctargopt.eq.1.or.ctargopt.eq.2
c     >    .or.ctargopt.eq.3.or.ctargopt.eq.6) then
c
c
c         The following code adjusts any data read in from an EDGE2 case
c         if that is necessary. (adjusts for DIVIMP removal of
c         virtual cells)
c
c            do 296 ir = irsep,nrs
c
c               do 298 ik = 1,nks(ir)
c                  e2dNBS (IK,IR) = e2dnbs(ik+1,ir)
c
c                  do iz = 0,cre2dizs
c                     e2dnzs(ik,ir,iz) = e2dnzs(ik+1,ir,iz)
c                  end do
c
c                  e2dVHS (IK,IR) = e2dvhs(ik+1,ir)
c                  e2dTIBS(IK,IR) = e2dtibs(ik+1,ir)
c                  e2dTEBS(IK,IR) = e2dtebs(ik+1,ir)
c                  e2dES  (IK,IR) = e2des(ik+1,ir)
c                  e2dion (IK,IR) = e2dion(ik+1,ir)
c                  e2dvro (ik,ir) = e2dvro(ik+1,ir)
c                  drho(ik,ir)    = drho(ik+1,ir)
c                  dthetag(ik,ir) = dthetag(ik+1,ir)
c                  e2dareas(ik,ir)= e2dareas(ik+1,ir)
c                  e2dgpara(ik,ir)= e2dgpara(ik+1,ir)
c                  e2dgdown(ik,ir)= e2dgdown(ik+1,ir)
c
c 298           continue
c 296        continue
c        endif
c
c       If readaux has been invoked - assign the e2dgdown values
c       to the fluxinfo array - if it does not already contain
c       other data - this is for use in e2dtargopt 2,3,4
c
        if (readaux.eq.1.and.fluxpts.le.0) then
           fluxpts = nrs-irsep+1
           do in = 1,fluxpts
              ir = in + irsep -1
              fluxinfo(in,1) = ir
              fluxinfo(in,2) = e2dgdown(nks(ir),ir)
              fluxinfo(in,3) = e2dgdown(2,ir)
              fluxinfo(in,4) = 0.0
           end do
         endif

c
c       Calculate an estimate of the Edge2D recombination - based on
c       specified option.
c
C
C---- ONE RING AT A TIME.
C
        DO IR = 1, NRS
          DO IK = 1, NKS(IR)
c
c            Artificial Temperature cutoff.
c
             if (ktebs(ik,ir).lt.treccut) then
                PTESA(IK) = treccut
             else
                PTESA(IK) = e2dTEBS(IK,IR)
            endif
c
            PNESA(IK) = e2dNBS(IK,IR) * RIZB
            PNBS(IK) = e2dNBS(IK,IR)
            PNHS(IK) = e2dATOM(IK,IR)
          ENDDO
C
c         Recombination
c
          write(year,'(i2.2)') iyearh
          call xxuid(useridh)
          ICLASS = 1
c
          if (crecopt.eq.4) then
             CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
          else
             call otherrec(NKS(IR),PTESA,PNESA,PCOEF,crecopt)
          endif
c
          DO IK = 1, NKS(IR)
            e2dhREC(IK,IR) = PNESA(IK)*PNBS(IK)*PCOEF(IK,1)
          ENDDO
        ENDDO

c
c
c     Elseif - for flag = 1 - only load some values
c     Load the set of values in the first cells of the SOL
c     rings into the cellvals array for use in SOL option 22.
c
      elseif (flag.eq.1) then
c
c         Read just the first cell  element values.
c
      DO 290 IR = irsep , NRS
C
         IREF   = IR - IRCENT + 1
C  USE INNERMOST EDGE RING FOR MISSING INNER RINGS
         IF( IR.LT.IRCENT )   IREF = 1
         IF( IR.GT.IRWALL-NINOWA .AND. IR.LE.IRWALL ) THEN
             IF( CNIWA.EQ.0 ) IREF = 0
C  USE OUTERMOST EDGE SOL RING FOR MISSING SOL RINGS
             IF( CNIWA.EQ.1 ) IREF = NXW
         END IF
         IF( IR.GT.IRWALL )   IREF = IREF - NINOWA
C
c        Ik = 1 (Outer end of ring - index = 2)
c
c        For the array cellvalls(ir,4,2)
c        index 1->4 is for:
c        1 - density
c        2 - Te
c        3 - Ti
c        4 - velocity as calculated in E2D
c
c        index 1->2
c        1 - inner  - ik = nks(ir) +2 because virtual points have been
c                          stripped by - or nj(ir)
c        2 - outer  - ik = 2
c
c
         K  = KORYE(IREF,2)
c
         cellvals(ir,1,2) = MAX (DUMMY(K,1),1.0) * 1.0E6
         cellvals(ir,2,2)  = MAX (DUMMY(K,4), LO)
         cellvals(ir,3,2)  = MAX (DUMMY(K,3), LO)
         cellvals(ir,4,2)  = DUMMY(K,2) * 0.01
c
         K  = KORYE(IREF,nj(iref)-1)
c
         cellvals(ir,1,1) = MAX (DUMMY(K,1),1.0) * 1.0E6
         cellvals(ir,2,1)  = MAX (DUMMY(K,4), LO)
         cellvals(ir,3,1)  = MAX (DUMMY(K,3), LO)
         cellvals(ir,4,1)  = DUMMY(K,2) * 0.01
c
c        write out
c
         write (6,*) 'cellvals:',ir,iref,k,nj(ir),nj(iref)
         write (6,*) 'Inner:',cellvals(ir,1,1),cellvals(ir,2,1),
     >               cellvals(ir,3,1),cellvals(ir,4,1)
         write (6,*) 'Outer:',cellvals(ir,1,2),cellvals(ir,2,2),
     >               cellvals(ir,3,2),cellvals(ir,4,2)

C
C
  290 CONTINUE
c
      endif
c
c

      return
C
 9999 IERR = 1
      WRITE (7,'(1X,A,3(/1X,A))')
     >  'TAUIN1: ERROR READING GEOMETRY DATA :-',MESAGE,
     >  'LAST LINE READ :-',BUFFER
      RETURN
C
 9001 FORMAT(/1X,'GEOMETRY DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
     >  /5X,'TOROIDAL FIELD     BPHI =  ',F8.3,
     >  /5X,'NO OF POINTS NP     =',I6,5X,'MAX NO. ROWS MKS    =',I6,
     >  /5X,'SEPARATRIX   IRSEP  =',I6,5X,'WALL         IRWALL =',I6,
     >  /5X,'FIRST TRAP   IRTRAP =',I6,5X,'NO OF RINGS  NRS    =',I6,
     >  /5X,'K SPLIT PT   IKTO   =',I6,5X,'K REF POINT  IKREF  =',I6)
 9002 FORMAT(/1X,' IK IR    R      Z     BPH',
     >  'I     TEB     TIB     NB      E1      VB        S  ',
     >  'BTOT/BTHETA    FEG1    FIG1',/1X,131('-'))
 9003 FORMAT(1X,2I3,3F7.3,2f8.1,1P,E8.1,0P,2A9,2F8.2,3X,2A9)
 9006 FORMAT(/1X,'TAUIN1: AREA OF SOL+TRAP =',F9.6,',  MAIN P =',F9.6,/)
 9011 FORMAT(/1X,'EDGE PLASMA DATA FOR SHOT',I6,',  TIME',F8.3,' :-',
     >  /5X,'RUN                 =',A,
     >  /5X,'PLASMA       HMASS  =',I6,5X,'IMPURITY     ZMASS  =',I6,
     >  /5X,'FIRST RING   IRCENT =',I6,5X,'LOST INN RGS NINOMP =',I6,
     >  /5X,'LOST OUT RGS NINOWA =',I6)
 9012 FORMAT(1P,6E12.4)
 9022 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
     >  ' IROUT  INDIST OUTDIST BACDIST FORDIST VOLUME  INPROB',
     >  '  POLDIST',/1X,131('-'))
 9023 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,2X,4F8.3,F8.5,F8.4,F8.3)
 9032 FORMAT(/1X,' IK IR    R      Z       NH0   TAUIZ0  CHPROB0 ',
     >  2('TAUIZ',I1,'  TAUEI',I1,'  TAURC',I1,'  CHPROB',I1,
     >  ' RCFRAC',I1,' '),/1X,131('-'))
 9033 FORMAT(1X,2I3,2F7.3,2X,1P,3E8.1,2(1P,4E8.1,0P,F8.4))
 9040 FORMAT(
     >   5X,'NO OF R PTS       NXS   =  ',I5,',   DELTAR =',F9.5,
     >  /5X,'NO OF Z PTS       NYS   =  ',I5,',   DELTAZ =',F9.5,
     >  /5X,'CENTRE      (RO  ,ZO  ) = (',F8.3,',',F8.3,')',
     >  /5X,'X POINT     (RXP ,ZXP ) = (',F8.3,',',F8.3,')',
     >  /5X,'LOWER LEFT  (RMIN,ZMIN) = (',F8.3,',',F8.3,')',
     >  /5X,'UPPER RIGHT (RMAX,ZMAX) = (',F8.3,',',F8.3,')',
     >  /5X,'INNER TARGET PTS  NDSIN =  ',I5,
     >  /5X,'TOTAL TARGET PTS  NDS   =  ',I5)
 9042 FORMAT(/1X,' IK IR    R      Z   IKIN IRIN IKOUT',
     >  ' IROUT  ID   INTHETA OUTTHETA THETA   DIST',
     >  /1X,131('-'))
 9043 FORMAT(1X,2I3,2F7.3,I5,I3,I8,I3,I7,2X,3F8.2,F8.5)
 9044 FORMAT(I7,E14.0,A40)
 9050 FORMAT(I5,I5,50F8.3)
      end
c
c
c
      subroutine TARGFLUX
c slmod begin
      USE mod_eirene06
c slmod end 
      IMPLICIT NONE
C
C*********************************************************************
C
C     TARGFLUX: Calculates the DIVIMP estimate of particle flux
c               onto the target.
C
C*********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     include "cneut"
c      include 'cneut'
c     include "dynam2"
c      include 'dynam2'
c     include "dynam3"
c      include 'dynam3'
c     include "dynam4"
c      include 'dynam4'
c     include "pindata"
c
      include 'pindata'
      include 'cadas'
      include 'outbuffer'
c
c slmod begin - new
      INCLUDE 'slcom'

      INTEGER fp,i1,i2,i
      REAL    puffsrc,addion,addiont,rc
c slmod end

C
      INTEGER IK,IR,id
C
      real TOTFLX, TOTREC
c
      REAL FLUX(MAXNRS),OUFLUX(MAXNRS),INFLUX(MAXNRS),
     >     OTOTFLX,ITOTFLX,cs,
     >     OUHEATE(MAXNRS),INHEATE(MAXNRS),HEATE(MAXNRS),
     >     OUHEATI(MAXNRS),INHEATI(MAXNRS),HEATI(MAXNRS),
     >     coreiz,soltotflx,pptotflx,pinescpd
c
      real alliz,ringiz,ringarea,ringneut,ringizdist,ringizsepdist
c
      real totflxt, soltotflxt,pptotflxt,fact
      real pincoreizt,pinsolizt,pinppizt ,allizt,totrect
      real totrece,totrecte,eircor,pincortmp
      integer saveunit
C
C---- Calculate the Particle Fluxes into the Divertor
C
c---- Initialization
c
      write (6,*) 'Hioniz2:',hioniz
c
      pincor = 1.0
c
      call rzero(flux,maxnrs)
      call rzero(ouflux,maxnrs)
      call rzero(influx,maxnrs)
      call rzero(heate,maxnrs)
      call rzero(ouheate,maxnrs)
      call rzero(inheate,maxnrs)
      call rzero(heati,maxnrs)
      call rzero(ouheati,maxnrs)
      call rzero(inheati,maxnrs)
c
      TOTFLX = 0.0
      ITOTFLX = 0.0
      OTOTFLX = 0.0
c
      soltotflx=0.0
      pptotflx =0.0
c
      totflxt = 0.0
      soltotflxt=0.0
      pptotflxt =0.0
c
      DO 200 ID=1,NDS,1
        IK = IKDS(ID)
        IR = IRDS(ID)
c
c        if (cioptf.eq.22) then
c           if (ik.eq.1) then
c              CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,2)
c           elseif (ik.eq.nks(ir)) then
c              CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,1)
c           endif
c        else
c           CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB)
c        endif
c
        cs = abs(kvds(id) )
c
        IF (IK.EQ .NKS(IR)) THEN
          INFLUX(IR) = KNDS (ID) * CS / KBFST (IR,1) *
     >                 COSTET (ID) * DDS2 (ID)
c slmod begin
          IF (supflx(IKHI,ir).EQ.1) influx(ir) = influx(ir) * 1.0E-15
c slmod end
          INHEATI(IR) =  2 * ECH * ( KTIDS (ID) ) * INFLUX (IR)
          INHEATE(IR) =  5 * ECH * ( KTEDS (ID) ) * INFLUX (IR)
          ITOTFLX = ITOTFLX + INFLUX(IR)
        ELSEIF (IK.EQ.1) THEN
          OUFLUX(IR) = KNDS (ID) * CS / KBFST (IR,2) *
     >                 COSTET (ID) * DDS2 (ID)
c slmod begin
          IF (supflx(IKLO,ir).EQ.1) ouflux(ir) = ouflux(ir) * 1.0E-15
c slmod end
          OUHEATI(IR) =  2 * ECH * ( KTIDS (ID) ) * OUFLUX (IR)
          OUHEATE(IR) =  5 * ECH * ( KTEDS (ID) ) * OUFLUX (IR)
          OTOTFLX = OTOTFLX + OUFLUX(IR)
        ENDIF
c
c       Calculate totals for each ring
c
        write (6,'(a,3i4,4(1x,g12.5))') 'targflux:',id,ir,ik,
     >                             dds2(id),cs,costet(id)
c
 200  CONTINUE
c
      do ir = irsep,nrs
         FLUX(IR)=OUFLUX(IR)+INFLUX(IR)
         HEATE(IR)=OUHEATE(IR)+INHEATE(IR)
         HEATI(IR)=OUHEATI(IR)+INHEATI(IR)
c
         if (ir.le.irwall) then
            soltotflx = soltotflx + flux(ir)
            soltotflxt = soltotflxt + ouflux(ir) * rp(idds(ir,2))/r0 +
     >                 influx(ir) * rp(idds(ir,1))/r0
         else
            pptotflx = pptotflx + flux(ir)
            pptotflxt = pptotflxt + ouflux(ir) * rp(idds(ir,2))/r0 +
     >                 influx(ir) * rp(idds(ir,1))/r0
         endif
c
      end do
c
      totflxt = soltotflxt + pptotflxt
C
      TOTFLX = OTOTFLX + ITOTFLX
C
C---- Calculate the number of recombinations
C
      TOTREC = 0.0
      totrect = 0.0
c
      TOTRECe = 0.0
      totrecte = 0.0
c
C
C---- ONE RING AT A TIME.
C
      DO IR = 1, NRS
        DO IK = 1, NKS(IR)
          TOTREC = TOTREC + DIVREC(IK,IR)*KAREA2(IK,IR)
          TOTRECt = TOTRECt + DIVREC(IK,IR)*KAREA2(IK,IR)*rs(ik,ir)/r0
          TOTRECe = TOTRECe + PINREC(IK,IR)*KAREA2(IK,IR)
          TOTRECte = TOTRECte+ PINREC(IK,IR)*KAREA2(IK,IR)*rs(ik,ir)/r0
        ENDDO
      ENDDO
c
c     Calculate ionization in various regions - also various sums for each ring
c
      pincoreiz = 0.0
      pinsoliz = 0.0
      pinppiz = 0.0
      alliz = 0.0
c
      pincoreizt = 0.0
      pinsolizt = 0.0
      pinppizt = 0.0
      allizt = 0.0
c
      do ir = 1,nrs
c
         ringiz = 0.0
         ringarea = 0.0
         ringneut = 0.0
c
         do ik = 1,nks(ir)
c
            if (ir.lt.irsep) then
c
c               if (ik.eq.nks(ir)) then
c                  write (6,*) 'Coreiz:',ir,pinion(ik,ir),
c     >                                   pinion(1,ir)
c               else
c
                if (ik.ne.nks(ir)) then
c
                  pincoreiz = pincoreiz
     >                   + pinion(ik,ir)*karea2(ik,ir)
                  pincoreizt = pincoreizt
     >               + pinion(ik,ir)*karea2(ik,ir)*rs(ik,ir)/r0
c
               endif
c
            elseif (ir.ge.irsep.and.ir.le.irwall) then
               pinsoliz = pinsoliz
     >                  + pinion(ik,ir)*karea2(ik,ir)
               pinsolizt = pinsolizt
     >           + pinion(ik,ir)*karea2(ik,ir)*rs(ik,ir)/r0
            elseif (ir.gt.irwall) then
               pinppiz = pinppiz
     >                 + pinion(ik,ir)*karea2(ik,ir)
               pinppizt = pinppizt
     >           + pinion(ik,ir)*karea2(ik,ir)*rs(ik,ir)/r0
            endif
c
            if (.not.(ik.eq.nks(ir).and.ir.lt.irsep)) then
               alliz = alliz + pinion(ik,ir) * karea2(ik,ir)
               allizt = allizt + pinion(ik,ir) * karea2(ik,ir)
     >                                      *rs(ik,ir)/r0
               ringiz = ringiz + pinion(ik,ir) * karea2(ik,ir)
               ringarea = ringarea + karea2(ik,ir)
               ringneut = ringneut + pinatom(ik,ir) *karea2(ik,ir)
            endif
c
            if (cprint.eq.9)
     >         write (6,'(2i4,6(1x,g16.8))') ik,ir,alliz,pincoreiz,
     >             pinsoliz,pinppiz,pinion(ik,ir),karea2(ik,ir)
c
         end do
c
c        Calculate ionization and neutral content on each ring
c
         piniz_info(ir,1) = ringiz
         piniz_info(ir,3) = ringneut
c
         if (ringarea.ne.0.0) then
            piniz_info(ir,2) = ringiz/ringarea
            piniz_info(ir,4) = ringneut/ringarea
         endif
c
      end do
c
c      CALL PRRMATDIV(PINREC,MAXNKS,nks(irsep),NRS,6,'RECOMBINATION')
c      CALL PRRMATDIV(PINION,MAXNKS,nks(irsep),NRS,6,'IONIZATION')
c
c  Calculate correction factor to ensure particle balance.  The PIN
c  recycling source, SRECYC, will be somewhat different from TOTFLX
c  since PIN uses a toroidal correction --> Use the fraction of
c  particles escaping to the core and the pump from PIN and renormalise
c  to what DIVIMP thinks is the total outflux of ions.
c
c
c
c
c
c     NIMBUS
c
c
c     ** NOTE **
c
c     TEMPORARY FIX TO GET TEH PIN OUTPUT INTO ANOTHER FILE
c
      saveunit = datunit
      datunit  = pinunit
c
      write(6,*) 'PINUNIT:',pinunit,datunit,saveunit,tmpunit
c
      if (pincode.eq.0) then
c
c      Change definition of "pinescpd" - or in other words - the additional
c      PIN source that may be contributed/lost through puffing and escape
c      to the pump and albedo regions. This form is taken from Lorne's fax
c      of Sept 9/96.
c
c         pinescpd = ((1.0-acthpcpuf)*hescal + (hescpd-hescal))
c
         pinescpd = ((1.0-acthpcpuf)*phfgal+phfuga)*(srecyc+srecom)



c
         PINCORtmp = (TOTFLXt+TOTRECt)/(allizt+pinescpd)
c
c        Since the correction factor is near 1 - keep the
c        analysis for consistency checking - but set the
c        actual correction to 1.0
c
         if (abs(pincortmp-1.0).gt.1) then
            pincor = pincortmp
         else
            pincor = 1.0
         endif
c
c        Temporarily force PINCOR=1.0
c
         pincor = 1.0

c
c         PINCOR2 = (TOTFLX+TOTREC)/(HIONIZ+pinesctmp)
c
c
      call prc('NIMBUS/PIN Information on Hydrogenic Behaviour')
      WRITE (7,'(1X,''PIN RANDOM NUMBER SEED:'',10X,I15)') PINSEED
      call prc('  NOTE: See LDH Memo - 9 Sept. 1996')
      call prc('  The following are for toroidal geometry:')
      call prr('  SRECYC - recycling source        :',srecyc)
      call prr('  SRECOM - recombination source    :',srecom)
      call prr('  PINFLX - SRECYC - HPUFF          :',
c     >              srecyc-acthpcpuf*hescal_last)
     >        srecyc -  (acthpcpuf*phfgal_last*
     >             ((srecyc-acthpcpuf*hescal_last+srecom)/
     >             (1.0-acthpcpuf*phfgal_last))
     >           + ppcpuf*(srecyc-acthpcpuf*hescal_last)))
c
      call prr('  HESCPD - HESCAL + core losses    :',hescpd)
      call prr('  CORELS - core losses             :',
     >                                          hescpd-hescal)
      call prr('  HESCAL - albedo losses (pump)    :',hescal)
      call prr('  HESCLK - losses to bypass leak   :',hesclk)
      call prr('  HPUFF  - total puff              :',
     >            acthpcpuf*phfgal_last*
     >             ((srecyc-acthpcpuf*hescal_last+srecom)/
     >             (1.0-acthpcpuf*phfgal_last))
     >              +ppcpuf*(srecyc-acthpcpuf*hescal_last))
      call prr('  HPCPUF - repuffed albedo fraction:',acthpcpuf)
      call prr('  PPCPUF - puffed target fraction  :',ppcpuf)
      call prr('  TOTSOU - source (SRECYC+SRECOM)  :',srecyc
     >                                               +srecom)
      call prr('  PHFGAL - HESCAL/TOTSOU           :',phfgal)
      call prr('  PHFUGA - (TOT.LOSS-HESCAL)/TOTSOU:',phfuga)
      call prc('DIVIMP/PIN Information on Hydrogenic Behaviour')
      call prc('  The following are for cylindrical/toroidal'//
     >                             ' geometry:')
      call prr2('  PIN Ionization (HIONIZ)  :',hioniz,allizt)
      call prr2('    PIN Core Ionization    :',pincoreiz,
     >                                         pincoreizt)
      call prr2('    PIN SOL  Ionization    :',pinsoliz,
     >                                         pinsolizt)
      call prr2('    PIN PP   Ionization    :',pinppiz,pinppizt)
      call prr ('    PIN H2   Ionization    :',h2ioniz)
      call prr('  PIN Escaped (Toroidal)   :',pinescpd)
      call prc('    PIN Escaped = [(1-HPCPUF)*HESCAL+'//
     >                                      '(HESCPD-HESCAL)]  OR')
      call prc('                  [((1-HPCPUF)*PHFGAL+PHFUGA)'//
     >                                      '*(SRECYC+SRECOM)]')
      call prr2('  DIV flux to targets      :',totflx,totflxt)
      call prr2('    DIV flux to SOL targets:',soltotflx,
     >                                        soltotflxt)
      call prr2('    DIV flux to PP  targets:',pptotflx,pptotflxt)
      call prr2('  DIV total recombinations :',totrec,totrect)
      call prr2('DIV Total Source (REC+FLX) :',totrec+totflx,
     >             totrect+totflxt)
      call prc('PINCOR = (DIV Total Source) / (HIONIZ + PINESCPD)')
      call prr('PINCOR factor (DIV/PIN)(Toroidal Quantities)= ',
     >         pincortmp)
      call prr('Actual Scaling factor used for PIN results  :',
     >         pincor)
      call prb
      call prc ('PIN/NIMBUS calculated GAUGE pressures:')
      call prc ('PUMP GAUGE:')
      call prr2('     IN/OUT  FLUXES      :',gaugedat(1,1),
     >                                       gaugedat(2,1))
      call prr ('     DENSITY             :',gaugedat(3,1))
      call prr ('     PRESSURE (MB)       :',gaugedat(4,1))
      call prc ('IN-VESSEL GAUGE:')
      call prr2('     IN/OUT  FLUXES      :',gaugedat(1,2),
     >                                       gaugedat(2,2))
      call prr ('     DENSITY             :',gaugedat(3,2))
      call prr ('     PRESSURE (MB)       :',gaugedat(4,2))
      call prb
c
c
c         write (6,*) 'PINCOR:',pincor,totflx,totrec,hioniz
c         write (6,*) 'srec:',srecyc,srecom,hescpd,pinescpd
c
c
c     Sonnet grids
c
      elseif (pincode.eq.1.or.pincode.eq.2.or.pincode.eq.3.or.
     .        pincode.eq.4.or.pincode.eq.5) then
c
c
c slmod begin - new

c...    Calculate ionisation on additional cells:
        addion  = 0.0
        addiont = 0.0
        i1 = 0
        DO i = 1, asc_ncell*ascncut
c...      3D vacuum grid:
          i1 = i1 + 1
          IF (i1.EQ.asc_ncell+1) i1 = 1

c...      Estimate the r value of the cell centroid (lame method):
          rc = 0.0
          DO i2 = 1, ascnvertex(i1), 2
            rc = rc + ascvertex(i2,i1)
          ENDDO
          rc = rc / REAL(ascnvertex(i1))

c...      Need some check of the option used to calculate
c         the additional cell volume:
          addion  = addion  + ascdata(i,1) * asc_vol(i)
          addiont = addiont + ascdata(i,1) * asc_vol(i) / rxp * rc

          WRITE(SLOUT,*) 'ADDVOL:',i1,asc_vol(i)
        ENDDO

        IF (.TRUE.) THEN
          fact = 2.0 * PI * rxp * eirtorfrac

          alliz      = fact * alliz 
          pincoreiz  = fact * pincoreiz 
          pinsoliz   = fact * pinsoliz 
          pinppiz    = fact * pinppiz 
          totrece    = fact * totrece
          totflx     = fact * totflx * eirsrcmul
          pptotflx   = fact * pptotflx * eirsrcmul
          soltotflx  = fact * soltotflx  * eirsrcmul
          totrec     = fact * totrec 

          fact = 2.0 * PI * r0 * eirtorfrac

          allizt     = fact * allizt
          pincoreizt = fact * pincoreizt
          pinsolizt  = fact * pinsolizt
          pinppizt   = fact * pinppizt
          totrecte   = fact * totrecte
          totflxt    = fact * totflxt * eirsrcmul
          pptotflxt  = fact * pptotflxt * eirsrcmul
          soltotflxt = fact * soltotflxt * eirsrcmul
          totrect    = fact * totrect
        ENDIF

c...    Calculate the total puffed source:
        puffsrc = 0.0
        IF (pincode.EQ.5) THEN
          DO i1 = 1, nstrata
            IF (NINT(strata(i1)%type).EQ.3) 
     .        puffsrc = puffsrc + strata(i1)%flux / ECH
          ENDDO
        ELSE
          DO i1 = 1, eirnpuff
            puffsrc = puffsrc + eirpuff(i1,MAXPUFF) / ECH
          ENDDO
        ENDIF

        EIRCOR = (TOTFLX + TOTRECe + puffsrc) / (alliz+hescpd+addion)

c...    Only output to .dat file for the last PIN call (also,
c       find the appropriate replacement for rel_count):

c        write(0,*) 'REL_COUNT:', rel_count, citersol,nitersol
c
c       jdemod
c     
c       Note rel_count appears to start at zero ... so this output is never
c       printed if iteration is on since rel_count is always less than nitersol
c       on the last iteration
c

        IF (citersol.EQ.0.OR.nitersol.EQ.0.OR.
     .      (rel_count+1).EQ.nitersol) THEN
          fp = datunit
        ELSE
          fp = PINOUT
        ENDIF

        WRITE(fp,*)
        WRITE(fp,85) 'DIVIMP/EIRENE Information on Hydrogenic '//
     .               'Behaviour'
        WRITE(fp,85) 'The following are for cylindrical/toroidal '//
     .               'geometry. Note that it is the'
        WRITE(fp,85) 'toroidal target flux that is currently passed'//
     .               ' to EIRENE:'
        WRITE(fp,87) 'EIRENE ionization        :',alliz    ,allizt
        WRITE(fp,87) '  EIR Core ionization    :',pincoreiz,pincoreizt
        WRITE(fp,87) '  EIR SOL  ionization    :',pinsoliz ,pinsolizt
        WRITE(fp,87) '  EIR PP   ionization    :',pinppiz  ,pinppizt
        WRITE(fp,87) '  EIR ADD  ionization    :',addion   ,addiont
        WRITE(fp,87) 'EIRENE recombinations    :',totrece  ,totrecte
        WRITE(fp,87) 'EIRENE puffed source (D) :',puffsrc  ,puffsrc
        WRITE(fp,87) 'DIV flux to targets      :',totflx   ,totflxt
        WRITE(fp,87) '  DIV flux to SOL targets:',soltotflx,soltotflxt
        WRITE(fp,87) '  DIV flux to PP  targets:',pptotflx ,pptotflxt
        WRITE(fp,87) 'DIV recombinations       :',totrec   ,totrect
        WRITE(fp,87) 'DIV total source(REC+FLX):',totrec +totflx ,
     .                                            totrect+totflxt
c...    TOTFLX and TOTFLXT should really be passed back from EIRENE
c       and reported here -- fix:
        WRITE(fp,87) 'EIR total source(REC+FLX):',
     .    totrece +totflx +puffsrc,
     .    totrecte+totflxt+puffsrc

        WRITE(fp,89) 'Flux lost to core boundary',hescpd-hescal,
     .               ' (',(hescpd-hescal)/(totrecte+totflxt+puffsrc)*
     .               100.0,'% of toroidal EIR total source)'
        WRITE(fp,89) 'Flux lost to pump(s)     :',hescal,
     .               ' (',(hescal       )/(totrecte+totflxt+puffsrc)*
     .               100.0,'% of toroidal EIR total source)'
        WRITE(fp,85) 'EIRCOR = (EIR total source)/(EIR ioniz+'//
     .               'flux lost)'
        WRITE(fp,86) 'EIRCOR (CYL/TOR) =        ',
     .             eircor,
     .             (totflxt+totrecte+puffsrc)/(allizt+hescpd+addiont)
        WRITE(fp,85) 'Factor NOT applied to EIRENE data.  This '//
     .               'information is'
        WRITE(fp,85) 'only reported for the *final* call to EIRENE.'
        WRITE(fp,*)
        WRITE(fp,86) 'Toroidal fraction (EIRTORFRAC)=',eirtorfrac
        WRITE(fp,86) 'Global src multi   (EIRSRCMUL)=',eirsrcmul
 
        IF (pincode.EQ.5) THEN
          CALL HD(fp,'  EIRENE STRATA SUMMARY','EIRNUMPAR-HD',5,77)
          CALL PRB
          WRITE(fp,'(4X,A8,A6,3A12)') 
     .      'STRATUM','TYPE','NO. TRACKS','FLUXT','PTRASH(%)'
          DO i1 = 1, nstrata
            WRITE(fp,'(4X,I8,F6.1,I12,1P,E12.4,0P,F12.4)') i1,
     .        strata(i1)%type,
     .        strata(i1)%ipanu,
     .        strata(i1)%fluxt,
     .        strata(i1)%ptrash / strata(i1)%fluxt * 100.0
          ENDDO
        ELSE
          CALL HD(fp,'  NUMBER OF PARTICLE TRACKS','EIRNUMPAR-HD',5,67)
          CALL PRB
          WRITE(fp,90) 'STRATUM','TYPE','NO. TRACKS'
          DO i1 = 1, eirnstrata+eirnpuff
            WRITE(fp,91) i1,eirntracks(i1)
          ENDDO
        ENDIF

85      FORMAT(3X,A)
86      FORMAT(3X,A,   2(1X,F12.3))
87      FORMAT(3X,A,1P,2(1X,E12.4:))
88      FORMAT(3X,A,1P,1X,E12.4)
89      FORMAT(3X,A,1P,1X,E12.4,0P,A,F5.1,A)
90      FORMAT(4X,A8,A13)
91      FORMAT(4X,I8,I13) 

        CALL PRB
c
c         EIRCOR = (TOTFLX+TOTREC)/(alliz+totrece)
c
c
c      call prc('DIVIMP/EIRENE Information on Hydrogenic Behaviour')
c      call prc(' The following are for cylindrical/toroidal'//
c     >                             ' geometry:')
c      call prr2(' EIRENE Ionization (HIONIZ):',alliz,allizt)
c      call prr2('    EIR Core Ionization    :',pincoreiz,
c     >                                         pincoreizt)
c      call prr2('    EIR SOL  Ionization    :',pinsoliz,
c     >                                         pinsolizt)
c      call prr2('    EIR PP   Ionization    :',pinppiz,pinppizt)
c      call prr2(' EIRENE Tot. recombinations:',totrece,totrecte)
c      call prr2(' DIV flux to targets       :',totflx,totflxt)
c      call prr2('    DIV flux to SOL targets:',soltotflx,
c     >                                        soltotflxt)
c      call prr2('    DIV flux to PP  targets:',pptotflx,pptotflxt)
c      call prr2(' DIV total recombinations  :',totrec,totrect)
c      call prr2('DIV Total Source (REC+FLX) :',totrec+totflx,
c     >             totrect+totflxt)
c      call prr2('EIR Total Source (REC+FLX) :',totrece+alliz,
c     >             totrecte+allizt)
c
c      call prc('EIRCOR = (DIV REC+FLX) / (HIONIZ + EIRREC)')
c      call prc('- Does not include corrections for particles lost')
c      call prc('  from EIRENE due to pumping or to the core.')
c      call prr('EIRCOR factor (DIV/EIR) = ',eircor)
c      call prc('Factor NOT applied to EIRENE data.')
c slmod end
c
      endif
c
      if (pinprint.eq.1) then
c
         call prb
c
         call prc('RING SUMMARY OF PIN IONIZATION AND'//
     >         ' NEUTRAL CONTENT:')
c
         call prc('  RING   IONIZATION     IZ DENS.    NEUTRALS'//
     >            '    NEUT DENS.')
         do ir = 1,nrs
c
            if (ir.eq.irsep) then
               call prc('    ------ Main SOL       ----------')
            elseif (ir.eq.irtrap) then
               call prc('    ------ Private Plasma ----------')
            endif
c
            write(comment,'(2x,i4,1p,4(1x,g12.4))') ir,
     >                       (piniz_info(ir,id),id=1,4)
            call prc(comment)
c
         enddo
c
         ringizdist = 0.0
         ringizsepdist = 0.0
         ringiz = 0.0
c
         do ir = 1,irsep-1
c
            ringizdist = ringizdist + piniz_info(ir,1) *
     >                                abs(middist(ir,2))
            ringiz     = ringiz      + piniz_info(ir,1)
c
            do ik = 1,nks(ir)-1

               ringizsepdist = ringizsepdist +
     >                         pinion(ik,ir)*karea2(ik,ir)
     >                         *separatrix_dist(ik,ir)

            end do
c
         end do
c
         if (ringiz.ne.0.0) then
c
            call prb
            call prr('MID-PLANE   AVERAGED DEPTH OF CORE'//
     >            ' IONIZATION (M) : ',ringizdist/ringiz)
            call prr('CELL-CENTER AVERAGED DEPTH OF CORE'//
     >            ' IONIZATION (M) : ',ringizsepdist/ringiz)
         endif
c

         call prb
      endif
c
c     ** NOTE **
c
c     RESTORE VALUE OF DATUNIT
c
      datunit = saveunit
      write(6,*) 'PINUNIT:',pinunit,datunit,saveunit,tmpunit
c
      return
      END
c     
c     
c     
      subroutine OSKIN
      IMPLICIT NONE
C     
C*********************************************************************
C     
C     OSKIN :  Extracts Transport Coefficients from the Onion
C     Skin background plasma model generated in DIVIMP
C     
C     Ray Monk (RHBNC) September 1993
C     
C     Using the skeleton of Program 'OUT' by
C     Chris Farrell  (Hunterskil)  February 1989
C     
C     *******************************************************************
C     
c     include "params"
      include 'params'
c     include "cgeom"
      include 'cgeom'
c     
c     Edge2D data
c     
      include 'cedge2d'
c     
c     Transport Data
c     
      include 'transcoef'
c     
      include 'comtor'
      include 'pindata'
      include 'printopt'
c     slmod begin
c...  TMP
      include 'slcom'
c     slmod end
c     
C     
C***********************************************************************
C     
C     Useful predefined variables :-
C     
C     NRS                  Number of Onion Skins perp to s
C     NKS(NRS)             Number of Elements along s on given ring
C     RS,ZS(NKS,NRS)       Machine coordinates of element
C     KTEBS(NKS,NRS)       Electron Temperature (eV)
C     KTIBS(NKS,NRS)       Ion Temperature (eV)
C     KNBS(NKS,NRS)        Electron Density (m**-3)
C     KES(NKS,NRS)         Electric Field (V/m)
C     KVHS(NKS,NRS)        Flow Velocity (m/s)
C     IRSEP                First SOL ring outside core
C     IRWALL               Last open ring inside the wall
C     RXP,ZXP              Coordinates of X-point
C     R0,Z0                Coordinates of Plasma Centre
C     KSS(NKS,NRS)         Distance along s
C     KSMAXS(NRS)          Maximum Distance along s
C     CIZB (RIZB)          Plasma Ion Charge
C     CRMB                 Plasma Ion Mass
C     QTIM                 Timestep for Ions (s)
C     BPHI                 Toroidal Field Bt (T)
C     KBFS(NKS,NRS)        Incidence Angle (Btot/Bpol)
C     
C***********************************************************************
C     
      INTEGER IK,IR,IKMID,OUEND,INEND,
     >     ID,LOOP
C     
      REAL FLUX(MAXNRS),OUFLUX(MAXNRS),INFLUX(MAXNRS),
     >     TOTFLX,OTOTFLX,ITOTFLX,
     >     totflx_save,ototflx_save,itotflx_save,
     >     netflx_save,onetflx_save,inetflx_save,
     >     DELN,ODELN,IDELN,
     >     DELTE,ODELTE,IDELTE,
     >     DELTI,ODELTI,IDELTI,
     >     IONIS(MAXNRS),OIONIS(MAXNRS),IIONIS(MAXNRS),
     >     NETFLX(MAXNRS),ONETFLX(MAXNRS),INETFLX(MAXNRS),
     >     NETQE(MAXNRS),PYTHGI,PYTHGO,
     >     IAREA,OAREA,PYTHG,AREA,SEP,OSEP,ISEP,
     >     IAREAO,OAREAO,AREAO,
     >     TOTQI,OTOTQI,ITOTQI,
     >     TOTQE,OTOTQE,ITOTQE,
c     
     >     OUHEATE(MAXNRS),INHEATE(MAXNRS),HEATE(MAXNRS),
     >     OUHEATI(MAXNRS),INHEATI(MAXNRS),HEATI(MAXNRS),
     >     cs,coreiz,alliz,soliz,ppiz,
     >     isepf,isepb,osepf,osepb,sepb,sepf,
     >     ni,no,nav,gai,gae,
     >     distx,distxmin,dist0,dist0min,
     >     tionis,itionis,otionis,
     >     surf(maxnks),delqe(maxnks),delqi(maxnks),deldp(maxnks),
     >     qetoto,qetoti,qitoto,qitoti,dptoto,dptoti
      real qetot,qitot
c     
c     Power Loss terms
c     
c     integer dpploss
c     parameter (dpploss=1)
c     
      real ionploss(maxnrs),iionploss(maxnrs),oionploss(maxnrs)
      real elecploss(maxnrs),ielecploss(maxnrs),oelecploss(maxnrs)
      real tiploss, teploss,peifact
c     
      real qecond(maxnks,maxnrs),
     >     qicond(maxnks,maxnrs),dpflux(maxnks,maxnrs),
     >     qeconv(maxnks,maxnrs),qiconv(maxnks,maxnrs),
     >     surfa(maxnks,maxnrs),
     >     dgradn(maxnks,maxnrs)
      real cfgam(maxnrs),icfgam(maxnrs),ocfgam(maxnrs)
      real cfqe(maxnrs),ocfqe(maxnrs),icfqe(maxnrs)
      real cfqi(maxnrs),ocfqi(maxnrs),icfqi(maxnrs)
      real delouti,delouto,qeconvi,qeconvo,qiconvi,qiconvo
      real qeouto,qeouti,qiouto,qiouti,qeoutconvi,qeoutconvo,
     >     qioutconvi,qioutconvo
c     
      real dgradt,dperpt,aperp,apdg
c     
      integer avmod
c     
      integer ikin,irin,ikout,irout,ikstart,ikend,ird
      integer irskip,irtmp,in1,in2,in
      real    totgradn,totgradte,totgradti
      real    totcfflux,totcffluxe2d, gradnsep,totrec
      real    ototcfflux,itotcfflux
      logical i
      real    flxval(maxnrs,13,3)
      logical usedpav,usexpav
      real    dpav,xpavi,xpavio,xpavii,xpave,xpaveo,xpavei,dptmp,
     >     dpavi,dpavo,fact,xpavt,xpavti,xpavto,
     >     tmpnflx,tmpion,tmpflx
c     
c     For calculation of volumetric quantities.
c     
      real srcs(maxnks,maxnrs,8)
      real srcsint(maxnrs,12)
      real srcstot(12)
c     
      real Quant2Grad
      external quant2grad
      real grad2n(maxnks,maxnrs)
      real grad2ti(maxnks,maxnrs)
      real grad2te(maxnks,maxnrs)
      real totgrad2n(maxnrs)
      real totgrad2ti(maxnrs)
      real totgrad2te(maxnrs)
c     
c     real pei(maxnks,maxnrs)
c     
c     Spline variables
c     
      external tg01b
      integer n1
      real x1(maxnrs),f1(maxnrs),dist
      REAL FDASH1(maxnrs),WORK(3*maxnrs),TG01B
c     slmod begin
      IF (grdnmod.NE.0.OR.iflexopt(8).EQ.11) THEN
         WRITE(0,*)
         WRITE(0,*)'-------------------------------------------------'
         WRITE(0,*) '           NOT EXECUTING OSKIN ROUTINE'
         WRITE(0,*)'-------------------------------------------------'
         WRITE(0,*)
         RETURN
      ENDIF
c     slmod end
c     
C     >     QLOSS(MAXNRS)
c     
c     CHARACTER coment*72
C     
C***********************************************************************
C     
C     Local variables :-
C     
C     TITLE                Title of DIVIMP Run
C     JOB                  Code for the DIVIMP Job
C     IR                   Current Onion Skin Ring (across s)
C     IK                   Current Onion Skin Ring (along s)
C     IKMID                Mid Point between Targets
C     INEND                Inner Target Divertor Entrance
C     OUEND                Outer Target Divertor Entrance
C     
C***********************************************************************
C     
C     
C-----------------------------------------------------------------------
C---- START OF MAIN ROUTINE
C-----------------------------------------------------------------------
C     
c     
c     Set Gamma values - these should be set to match the SOL option
c     and eventually may be set automatically
c     
      gai = 3.5
      gae = 5.0
c     
c     
      dpav  = 0.0
      dpavi = 0.0
      dpavo = 0.0
c     
      xpave  = 0.0
      xpavei = 0.0
      xpaveo = 0.0
c     
      xpavi  = 0.0
      xpavii = 0.0
      xpavio = 0.0
c     
      xpavt  = 0.0
      xpavti = 0.0
      xpavto = 0.0
c     
      usedpav = .false.
      usexpav = .false.
c     
c     irskip = (irwall - irsep) / 2
c     
      irskip = dpavopt
c     
c     irskip = 7
c     
      write (6,*) 'irskip:',irskip,irwall,irsep,' R0 = ',r0
c     
C     
C---- Determine the Midpoint between the Targets
C     
c     IR=IRSEP
c     DO 80 IK=1,NKS(IR),1
c     IF (KSS(IK,IR).GE.(0.5*KSMAXS(IR))) THEN
c     IKMID=IK
c     GOTO 90
c     ENDIF
c     80   CONTINUE
c     90   continue
C     
C---- Locate Elements at Divertor Entrance (i.e. at X-point height)
C     
C     
C---- and Locate Elements at Plasma Midplane
c     
c     Do these for all cases since we want the R-coordinates for
c     the inner and outer mid-planes.
c     
c     slmod begin - new
c     ARRAY BOUNDS: IR used to be set to IRSEP above, but it is now commented
c     out.  IR needs to be set for the following code.
      ir = irsep
c     slmod end
      ikmid = ikmids(irsep)+1
C     
C     Outer Plates...
C     
      distxmin = hi
      dist0min = hi
      DO 100 IK=1,IKMID-1,1
         distx = ((ZS(IK,IR)-zxp)**2)
         IF (distx.lt.distxmin) THEN
            distxmin = distx
            OUEND=IK
         ENDIF
         dist0 = (ZS(IK,IR)-Z0)**2
         IF (dist0.lt.dist0min) THEN
            dist0min = dist0
            OUMID=IK
         ENDIF
 100  CONTINUE

C     
C     Inner Plates...
C     
      distxmin = hi
      dist0min = hi
      DO 120 IK=NKS(IR),IKMID,-1
         distx = (ZS(IK,IR)-zxp)**2
         IF (distx.lt.distxmin) THEN
            distxmin = distx
            INEND=IK
         ENDIF
         dist0 = (ZS(IK,IR)-Z0)**2
         IF (dist0.lt.dist0min) THEN
            dist0min = dist0
            INMID=IK
         ENDIF
 120  CONTINUE

c     
c     Note that rcinner will be descending order and rcouter will be
c     ascending - that is why rcouter will typically be used for
c     plotting.
c     
c     do ir = irsep,irwall
c     rcinner(ir) = rs(inmid,ir)
c     rcouter(ir) = rs(oumid,ir)
c     end do
c     
c     
C     
C---- Calculate the Particle Fluxes into the Divertor
C     
      DO 200 ID=1,NDS,1
         IK = IKDS(ID)
         IR = IRDS(ID)
c     
         if (dprcopt.eq.0) then
            fact = 1.0
         else
            fact = rp(id)/r0
         endif
c     
c     if (cioptf.eq.22) then
c     if (ik.lt.ikmid) then
c     CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,2)
c     elseif (ik.ge.ikmid) then
c     CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB) * cmachno(ir,1)
c     endif
c     else
c     CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
c     >                (1.0+RIZB)/CRMB)
c     endif
c     
         cs = abs(kvds(id))
c     
         IF (IK.EQ .NKS(IR)) THEN
            INFLUX(IR) = KNDS (ID) * CS / KBFST (IR,1) *
     >           COSTET (ID) * DDS2 (ID) * fact
            INHEATI(IR) =  GAI * ECH * ( KTIDS (ID) ) * INFLUX (IR)
            INHEATE(IR) =  GAE * ECH * ( KTEDS (ID) ) * INFLUX (IR)
         ELSEIF (IK.EQ.1) THEN
            OUFLUX(IR) = KNDS (ID) * CS / KBFST (IR,2) *
     >           COSTET (ID) * DDS2 (ID) * fact
            OUHEATI(IR) =  GAI * ECH * ( KTIDS (ID) ) * OUFLUX (IR)
            OUHEATE(IR) =  GAE * ECH * ( KTEDS (ID) ) * OUFLUX (IR)
         ENDIF
c     
C     
c     
 200  CONTINUE
c     
c     Copy target fluxes and save them for later
c     
      totflx_save = 0.0
      ototflx_save = 0.0
      itotflx_save = 0.0
c     
      do ir = irsep,irwall-1 
         totflx_save = totflx_save + ouflux(ir)+influx(ir)
         ototflx_save = ototflx_save + ouflux(ir)
         itotflx_save = itotflx_save + influx(ir)
      end do 
c     
c     Calculate a variety of source terms based on the
c     background plasma. This routine will likely be moved from
c     here at some point in time.
c     

      call rzero(srcs,maxnks*maxnrs*8)
      call rzero(srcsint,maxnrs*12)
      call rzero(srcstot,12)
c     
      do ir = irsep,nrs
c     
         if (ir.ne.irwall.and.ir.ne.irtrap) then
c     
c     if (ir.eq.irwall.or.ir.eq.irtrap) cycle
c     
            do ik = 1,nks(ir)
c     
c     Major Radius correction if requested
c     
               if (dprcopt.eq.0.or.dprcopt.eq.2) then
                  fact = 1.0
               elseif (dprcopt.eq.1) then
                  fact = rs(ik,ir)/r0
               endif
c     
c     Pei
c     
c     peifact = 0.5 * log(1.5e13 * ktebs(ik,ir)**1.5
c     >                   / sqrt(knbs(ik,ir)))  / 15.0
c     
               peifact = dppei * log(1.5e13 * ktebs(ik,ir)**1.5
     >              / sqrt(knbs(ik,ir)))  / 15.0


               srcs(ik,ir,1) = - 1.14e-32 * knbs(ik,ir)**2 *
     >              (ktebs(ik,ir)-ktibs(ik,ir)) /
     >              (crmb * ktebs(ik,ir)**1.5)
     >              * peifact
     >              * fact


c     
c     Phelpi
c     
               srcs(ik,ir,2) = (17.5 + (5.0+37.5/ktebs(ik,ir)) *
     >              (1.0+0.25/ktebs(ik,ir)) *
     >              log10(1.0e21/knbs(ik,ir))) *
     >              (pinion(ik,ir)*dprec)* ech
     >              * fact
c     
c     Pcx
c     
               srcs(ik,ir,3) = 1.5 * ktibs(ik,ir) * ech * 1.0 *
     >              (pinion(ik,ir) *dprec)
     >              * fact
c     
c     Smom
c     
               srcs(ik,ir,4) = - crmb * amu *kvhs(ik,ir)/qtim* 1.0*
     >              ( pinion(ik,ir) *dprec)
     >              * fact
c     
c     PINQI
c     
               srcs(ik,ir,5) =  pinqi(ik,ir)
     >              * fact * dprec
c     
c     PINQE
c     
               srcs(ik,ir,6) =  pinqe(ik,ir) * dp_pinqe_mult
     >              * fact * dprec
c     
c     Recombination
c     
               srcs(ik,ir,7) =  pinrec(ik,ir)
     >              * fact * dprec
c     
c     Ionization
c     
               srcs(ik,ir,8) =  pinion(ik,ir)
     >              * fact * dprec
c     
c     
c     Integrate over ring
c     
               srcsint(ir,1) = srcsint(ir,1)+srcs(ik,ir,1)*karea2(ik,ir)
               srcsint(ir,2) = srcsint(ir,2)+srcs(ik,ir,2)*karea2(ik,ir)
               srcsint(ir,3) = srcsint(ir,3)+srcs(ik,ir,3)*karea2(ik,ir)
               srcsint(ir,4) = srcsint(ir,4)+srcs(ik,ir,4)*karea2(ik,ir)
               srcsint(ir,5) = srcsint(ir,5)
     >              +dprec*pinion(ik,ir)*karea2(ik,ir)
               srcsint(ir,6) = srcsint(ir,6)
     >              +dprec*pinqi(ik,ir)*karea2(ik,ir)
               srcsint(ir,7) = srcsint(ir,7)
     >              +dprec*pinqe(ik,ir)*karea2(ik,ir)*dp_pinqe_mult
               srcsint(ir,8) = srcsint(ir,8)+abs(srcs(ik,ir,4))
     >              *karea2(ik,ir)
c     
               srcsint(ir,9) = srcsint(ir,5)+srcs(ik,ir,5)*karea2(ik,ir)
               srcsint(ir,10)= srcsint(ir,6)+srcs(ik,ir,6)*karea2(ik,ir)
               srcsint(ir,11)= srcsint(ir,7)+srcs(ik,ir,7)*karea2(ik,ir)
               srcsint(ir,12)= srcsint(ir,8)+srcs(ik,ir,8)*karea2(ik,ir)
c     
            end do
c     
c     Add up sources over grid
c     
            do in = 1,12
               srcstot(in) = srcstot(in) + srcsint(ir,in)
            end do
c     
c     Ring Test EndIF
c     
         end if
c     
      end do
c     
c     Print out the tallies of information - cell by cell and ring by ring
c     
 2000 format(8x,'IR',2x,'IK',9x,'Karea',11x, 'Pei',
     >     8x,'Phelpi',11x,'Pcx',10x,'Smom',11x,'Siz',
     >     8x,'PIN QI',8x,'PIN QE')
 2010 format(6x,2(i4),8(1x,g13.6))
 2020 format(a6,2(i4),14x,7(1x,g13.6))
 2025 format(a6,8x,14x,3(14x),1x,g13.6)
 2030 format(a14,14x,7(1x,g13.6))
 2035 format(a14,14x,3(14x),1x,g13.6)
c     
c     Only print if the Dperp extractor print out is ON.
c     
      if (cprint.eq.1.or.cprint.eq.9) then

c     
         write(6,*)
         write(6,*) 'Volumetric Source Term Summary:'
         write(6,*)
         write(6,*) 'Dperp Extractor Recycle Fraction = ',dprec
c     
         do ir = irsep,nrs
c     
c     if (ir.eq.irwall.or.ir.eq.irtrap) cycle
c     
            if (ir.ne.irwall.and.ir.ne.irtrap) then
c     
               write (6,2000)
c     
               do ik = 1,nks(ir)
c     
                  write(6,2010) ir,ik,karea2(ik,ir),
     >                 (srcs(ik,ir,in),in=1,4),dprec*pinion(ik,ir),
     >                 dprec*pinqi(ik,ir),dprec*pinqe(ik,ir)
     >                 * dp_pinqe_mult
               end do
c     
               write(6,2020) 'Total:',ir,nks(ir),
     >              (srcsint(ir,in),in=1,7)
               write(6,2025) 'Smom: ',srcsint(ir,8)
               write(6,*)

c     
c     Endif for IR test
c     
            endif

         end do
c     
         write (6,2030) 'Grand Totals:',(srcstot(in),in=1,7)
         write (6,2035) '        Smom:',srcstot(8)
         write(6,*)
c     
      endif
c     
c     
c     
c     Calculate ring totals
c     
      do ir = irsep,nrs
c     
         FLUX(IR)=OUFLUX(IR)+INFLUX(IR)
         HEATE(IR)=OUHEATE(IR)+INHEATE(IR)
         HEATI(IR)=OUHEATI(IR)+INHEATI(IR)
c     
      end do
c     
C     
C---- Determine the Elemental Volumes and Associated Pin
C---- Ionisation Levels (using Lorne's 'real' volume)
C     
c     CALL PRRMATDIV(E2DION,MAXNKS,nks(irsep),NRS,6,
c     >              'EDGE2D IONIZATION')
c     

      DO 220 IR=1,nrs,1
c     
         IONIS(IR)=0.0
         IIONIS(IR)=0.0
         OIONIS(IR)=0.0
c     
c     Power Loss
c     
         ionploss(ir) = 0.0
         iionploss(ir) = 0.0
         oionploss(ir) = 0.0
         elecploss(ir) = 0.0
         ielecploss(ir) = 0.0
         oelecploss(ir) = 0.0
c     
c     write (6,*) 'ir:',ir,nks(ir)
c     
         DO 230 IK=1,NKS(IR),1
c     
c     Using karea2
c     
            if (dprcopt.eq.0.or.dprcopt.eq.2) then
               fact = 1.0
            elseif (dprcopt.eq.1) then
               fact = rs(ik,ir)/r0
            endif
c     

            if (cioptg.eq.99.and.cre2d.eq.0) then
c     
               IONIS(IR)=IONIS(IR)+(E2DION(IK,IR)
     >              *Karea2(IK,IR))*fact
c     
               if (kss(ik,ir).lt.(ksmaxs(ir)/2.0)) then
                  OIONIS(IR)=OIONIS(IR)+(E2DION(IK,IR)
     >                 *Karea2(IK,IR))*fact
               else
                  IIONIS(IR)=IIONIS(IR)+(E2DION(IK,IR)
     >                 *Karea2(IK,IR))*fact
               endif
c     
            else
c     
               IONIS(IR)=IONIS(IR)+( (PINION(IK,IR)*dprec)
     >              *Karea2(IK,IR))*fact
c     
               if (dpploss.eq.1) then
                  ionploss(ir) = ionploss(ir) +
     >                 (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
                  elecploss(ir) = elecploss(ir) +
     >                 (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                 * dp_pinqe_mult
               elseif (dpploss.eq.2) then
                  ionploss(ir) = ionploss(ir)
     >                 + (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                 - srcs(ik,ir,1) * karea2(ik,ir) * fact
                  elecploss(ir) = elecploss(ir)
     >                 + (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                 * dp_pinqe_mult
     >                 + srcs(ik,ir,1) * karea2(ik,ir) * fact
               endif
c     
               if (kss(ik,ir).lt.(ksmaxs(ir)/2.0)) then
                  OIONIS(IR)=OIONIS(IR)+((PINION(IK,IR)*dprec)
     >                 *Karea2(IK,IR))*fact
c     
                  if (dpploss.eq.1) then
                     oionploss(ir) = oionploss(ir) +
     >                    (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
                     oelecploss(ir) = oelecploss(ir) +
     >                    (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    * dp_pinqe_mult
                  elseif (dpploss.eq.2) then
                     oionploss(ir) = oionploss(ir)
     >                    + (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    - srcs(ik,ir,1) * karea2(ik,ir) * fact
                     oelecploss(ir) = oelecploss(ir)
     >                    + (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    * dp_pinqe_mult
     >                    + srcs(ik,ir,1) * karea2(ik,ir) * fact
                  endif
c     
               else
                  IIONIS(IR)=IIONIS(IR)+((PINION(IK,IR)*dprec)
     >                 *Karea2(IK,IR))*fact
c     
                  if (dpploss.eq.1) then
                     iionploss(ir) = iionploss(ir) +
     >                    (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
                     ielecploss(ir) = ielecploss(ir) +
     >                    (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    * dp_pinqe_mult
                  elseif (dpploss.eq.2) then
                     iionploss(ir) = iionploss(ir)
     >                    + (pinqi(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    - srcs(ik,ir,1) * karea2(ik,ir) * fact
                     ielecploss(ir) = ielecploss(ir)
     >                    + (pinqe(ik,ir)*dprec)*Karea2(IK,IR)*fact
     >                    * dp_pinqe_mult
     >                    + srcs(ik,ir,1) * karea2(ik,ir) * fact
                  endif
c     
               endif
c     
c     
            endif


 230     CONTINUE
 220  CONTINUE

c     
c     Calculate ionization in core
c     
      coreiz = 0.0
      soliz = 0.0
      ppiz = 0.0
      alliz = 0.0
      do ir = 1,nrs
c     
         if (ir.lt.irsep) then
            coreiz = coreiz + ionis(ir)
         elseif (ir.ge.irsep.and.ir.le.irwall) then
            soliz = soliz + ionis(ir)
         elseif (ir.gt.irwall) then
            ppiz = ppiz + ionis(ir)
         endif
c     
         alliz = alliz + ionis(ir)
      end do
c     
      write(6,*) 'Core ionization from NIMBUS: ',coreiz
      write(6,*) 'All  ionization from NIMBUS: ',alliz
      write(6,*) 'Pin Correction factor: ',pincor
      write(6,*) 'Dperp Extractor Imposed Recycle Fraction = ',dprec
c     
c     
      write (6,1210)
c     
      DO 240 IR=IRSEP,nrs,1
c     
c     Net particle fluxes
c     
         NETFLX(IR) =FLUX(IR)  -IONIS(IR)
         INETFLX(IR)=INFLUX(IR)-IIONIS(IR)
         ONETFLX(IR)=OUFLUX(IR)-OIONIS(IR)
         write(6,1200) ir,INETFLX(IR),INFLUX(IR),IIONIS(IR),
     >        ONETFLX(IR),OUFLUX(IR),OIONIS(IR),
     >        NETFLX(IR),FLUX(IR),IONIS(IR)
c     
 240  CONTINUE

c     
      write(6,'(a)') 'Ion and Electron Energy Source Terms'
      write(6,1240) inner,outer
c     
      tiploss = 0.0
      teploss = 0.0
c     
      do ir = irsep,nrs
         tiploss = tiploss + ionploss(ir)
         teploss = teploss + elecploss(ir)

         if (cprint.eq.1.or.cprint.eq.9) then
            write(6,1230) 'Ion :',ir,ionploss(ir),iionploss(ir),
     >           oionploss(ir)
            write(6,1230) 'Elec:',ir,elecploss(ir),ielecploss(ir),
     >           oelecploss(ir)
         endif

      end do
      write (6,*) 'Total Ion      Ploss = ',tiploss
      write (6,*) 'Total Electron Ploss = ',teploss
c     
 1240 format(1x,' IR ',10x,'Total',10x,a,10x,a)
 1230 format(a,i4,3(1x,g14.6))

C     
C---- Calculate the NET Particle Flux Entering the Divertor Region
C     ~~~
      OTOTFLX=0.0
      ITOTFLX=0.0
      TOTFLX=0.0
c     
      OTOTQI=0.0
      ITOTQI=0.0
      TOTQI=0.0
c     
      OTOTQE=0.0
      ITOTQE=0.0
      TOTQE=0.0
c     
      Tionis  = 0.0
      iTionis = 0.0
      oTionis = 0.0
c     
      do ir = irsep,irwall-1
c     
c     Total Ionization

         Tionis = tionis + ionis (ir)
         iTionis = itionis + iionis (ir)
         oTionis = otionis + oionis (ir)
c     
c     Total particle fluxes
c     
         ITOTFLX=ITOTFLX + INETFLX(IR)
         OTOTFLX=OTOTFLX + ONETFLX(IR)
         TOTFLX =TOTFLX  + NETFLX(IR)
c     
c     Total heat fluxes
c     
         OTOTQI=OTOTQI+OUHEATI(IR)-oionploss(ir)
         OTOTQE=OTOTQE+OUHEATE(IR)-oelecploss(ir)
c     
         ITOTQI=ITOTQI+INHEATI(IR)-iionploss(ir)
         ITOTQE=ITOTQE+INHEATE(IR)-ielecploss(ir)
c     
         TOTQI=TOTQI+HEATI(IR)-ionploss(ir)
         TOTQE=TOTQE+HEATE(IR)-elecploss(ir)
c     
      end do
c     
c     Save netfluxes 
c     
      netflx_save = totflx
      onetflx_save = ototflx
      inetflx_save = itotflx
c     
c     Calculate the Net cross-field fluxes for
c     each ring.
c     
c     
      if (dpfluxopt.eq.0) then

         ir = irsep

         cfgam(ir) = totflx - netflx(ir)
         icfgam(ir) = itotflx - inetflx(ir)
         ocfgam(ir) = ototflx - onetflx(ir)
c     
         cfqe(ir) = totqe - heate(ir) + elecploss(ir)
         icfqe(ir) = itotqe - inheate(ir) + ielecploss(ir)
         ocfqe(ir) = ototqe - ouheate(ir) + oelecploss(ir)
c     
         cfqi(ir) = totqi - heati(ir) + ionploss(ir)
         icfqi(ir) = itotqi - inheati(ir) + iionploss(ir)
         ocfqi(ir) = ototqi - ouheati(ir) + oionploss(ir)
c     
         if (cprint.eq.1.or.cprint.eq.9) then
c     
            write(6,2001) 'cf :',ir,cfgam(ir),cfqe(ir),cfqi(ir)
            write(6,2001) 'cfi:',ir,icfgam(ir),icfqe(ir),icfqi(ir)
            write(6,2001) 'cfo:',ir,ocfgam(ir),ocfqe(ir),ocfqi(ir)
c     
         endif
c     
         do ir = irsep+1,irwall -1
            cfgam(ir)  = cfgam(ir-1)  - netflx(ir)
            icfgam(ir) = icfgam(ir-1) - inetflx(ir)
            ocfgam(ir) = ocfgam(ir-1) - onetflx(ir)
c     
            cfqe(ir)  = cfqe(ir-1)  - heate(ir) + elecploss(ir)
            icfqe(ir) = icfqe(ir-1) - inheate(ir)+ ielecploss(ir)
            ocfqe(ir) = ocfqe(ir-1) - ouheate(ir)+ oelecploss(ir)
c     
            cfqi(ir)  = cfqi(ir-1)  - heati(ir)+ ionploss(ir)
            icfqi(ir) = icfqi(ir-1) - inheati(ir)+ iionploss(ir)
            ocfqi(ir) = ocfqi(ir-1) - ouheati(ir)+ oionploss(ir)
c     
            if (cprint.eq.1.or.cprint.eq.9) then

               write(6,2001) 'cf :',ir,cfgam(ir),cfqe(ir),cfqi(ir)
               write(6,2001) 'cfi:',ir,icfgam(ir),icfqe(ir),icfqi(ir)
               write(6,2001) 'cfo:',ir,ocfgam(ir),ocfqe(ir),ocfqi(ir)
c     
            endif
c     
         end do
c     
 2001    format(a,1x,i4,3(1x,g14.6))

c     
c     DPFLUXOPT
c     
      elseif (dpfluxopt.eq.1) then

         ir = irsep

         cfgam(ir) = totflx - 0.5 * netflx(ir)
         icfgam(ir) = itotflx - 0.5 * inetflx(ir)
         ocfgam(ir) = ototflx - 0.5 * onetflx(ir)
c     
         cfqe(ir) = totqe - 0.5 * heate(ir)
         icfqe(ir) = itotqe - 0.5 * inheate(ir)
         ocfqe(ir) = ototqe - 0.5 * ouheate(ir)
c     
         cfqi(ir) = totqi - 0.5 * heati(ir)
         icfqi(ir) = itotqi - 0.5 * inheati(ir)
         ocfqi(ir) = ototqi - 0.5 * ouheati(ir)
c     
c     write(6,*) 'cf :',ir,cfgam(ir),cfqe(ir),cfqi(ir)
c     write(6,*) 'cfi:',ir,icfgam(ir),icfqe(ir),icfqi(ir)
c     write(6,*) 'cfo:',ir,ocfgam(ir),ocfqe(ir),ocfqi(ir)
c     
c     do ir = irsep+1,irwall -1
c     cfgam(ir)  = cfgam(ir-1)
c     >                -0.5*(netflx(ir)+netflx(ir-1))
c     icfgam(ir) = icfgam(ir-1)
c     >                -0.5*(inetflx(ir)+inetflx(ir-1))
c     ocfgam(ir) = ocfgam(ir-1)
c     >                -0.5*(onetflx(ir)+onetflx(ir-1))
c     
c     cfqe(ir)  = cfqe(ir-1)
c     >               - 0.5*(heate(ir)+heate(ir-1))
c     icfqe(ir) = icfqe(ir-1)
c     >               - 0.5*(inheate(ir)+inheate(ir-1))
c     ocfqe(ir) = ocfqe(ir-1)
c     >               - 0.5*(ouheate(ir)+ouheate(ir-1))
c     
c     cfqi(ir)  = cfqi(ir-1)
c     >               - 0.5*(heati(ir)+heati(ir-1))
c     icfqi(ir) = icfqi(ir-1)
c     >               - 0.5*(inheati(ir)+inheati(ir-1))
c     ocfqi(ir) = ocfqi(ir-1)
c     >               - 0.5*(ouheati(ir)+ouheati(ir-1))
c     
c     write(6,*) 'cf :',ir,cfgam(ir),cfqe(ir),cfqi(ir)
c     write(6,*) 'cfi:',ir,icfgam(ir),icfqe(ir),icfqi(ir)
c     write(6,*) 'cfo:',ir,ocfgam(ir),ocfqe(ir),ocfqi(ir)
c     
c     end do
c     


      endif
c     
c     
      write(6,*) 'Initial Values: (flux and ionization Main SOL)'
      write(6,*) 'ITOTFLX= ',itotflx,' OTOTFLX= ',ototflx,
     >     ' TOTFLX= ',totflx
      write(6,*) 'ITionis= ',itionis,' OTionis= ',otionis,
     >     ' Tionis= ',tionis
      write(6,*) ' ITFlux= ',itotflx+itionis,
     >     ' OTFlux= ',ototflx+otionis,
     >     '  Tflux= ',totflx+tionis
c     

c     
C---- Find the Total Surface Area of the Plasma and Mean Onion
C---- Skin Seperation...(Density Gradient!)
C     
c     
c     This is total plasma area adjacent to the main plasma.
c     
c     


      if (dpmethod.eq.0.or.dpmethod.eq.1) then
c     
c     
c     
         DO 280 IR=IRSEP,IRWALL-1,1
c     
            ikmid = ikmids(ir)+1
c     
            if (dpsuml.eq.1) then
               ouend =1
               inend =nks(ir)
            endif


c     
c     Implement alternate solution methods but leave in the code
c     Select by hard-wired switch
c     


            if (dpmethod.eq.0) then
c     
               AREA=0.0
               iarea= 0.0
               oarea= 0.0
c     
               DO 290 IK=OUEND,INEND,1
c     
                  if (dprcopt.eq.0) then
                     fact = 1.0
                  else
                     fact = rs(ik,ir)/r0
                  endif
c     
                  PYTHG=SQRT( (RS(IK,IR)-RS(IK+1,IR))**2 +
     >                 (ZS(IK,IR)-ZS(IK+1,IR))**2 ) * fact
c     
c     
                  AREA=AREA+ PYTHG
c     
c     
                  if (ik.lt.ikmid) then
                     oAREA=oAREA+ PYTHG
                  else
                     iAREA=iAREA+ PYTHG
                  endif
c     
 290           CONTINUE
c     
c     Calculate Gradients for both inner and outer legs.
c     
c     
               if (ir.eq.irsep) then

                  OSEP = SQRT((RS(OUMID,IR)-RS(OUMID,IR+1))**2 +
     >                 (ZS(OUMID,IR)-ZS(OUMID,IR+1))**2 )
                  ISEP = SQRT((RS(INMID,IR)-RS(INMID,IR+1))**2 +
     >                 (ZS(INMID,IR)-ZS(INMID,IR+1))**2 )

                  ODELN  = (KNBS(OUMID,IR)-KNBS(OUMID,IR+1))/OSEP
                  ODELTI = (KTIBS(OUMID,IR)-KTIBS(OUMID,IR+1))/OSEP
                  ODELTE = (KTEBS(OUMID,IR)-KTEBS(OUMID,IR+1))/OSEP
c     
                  IDELN  = (KNBS(inMID,IR)-KNBS(inMID,IR+1))/ISEP
                  IDELTI = (KTIBS(inMID,IR)-KTIBS(inMID,IR+1))/ISEP
                  IDELTE = (KTEBS(inMID,IR)-KTEBS(inMID,IR+1))/ISEP

               elseif (ir.eq.irwall-1) then

                  OSEP = SQRT((ABS(RS(OUMID,IR-1)-RS(OUMID,IR)))**2 +
     >                 (ABS(ZS(OUMID,IR-1)-ZS(OUMID,IR)))**2)
                  ISEP = SQRT((ABS(RS(INMID,IR-1)-RS(INMID,IR)))**2 +
     >                 (ABS(ZS(INMID,IR-1)-ZS(INMID,IR)))**2)

                  ODELN  = (KNBS(OUMID,IR-1)-KNBS(OUMID,IR))/OSEP
                  ODELTI = (KTIBS(OUMID,IR-1)-KTIBS(OUMID,IR))/OSEP
                  ODELTE = (KTEBS(OUMID,IR-1)-KTEBS(OUMID,IR))/OSEP
c     
                  IDELN  = (KNBS(inMID,IR-1)-KNBS(inMID,IR))/ISEP
                  IDELTI = (KTIBS(inMID,IR-1)-KTIBS(inMID,IR))/ISEP
                  IDELTE = (KTEBS(inMID,IR-1)-KTEBS(inMID,IR))/ISEP

               else

                  OSEPF = SQRT((RS(OUMID,IR)-RS(OUMID,IR+1))**2 +
     >                 (ZS(OUMID,IR)-ZS(OUMID,IR+1))**2 )
                  ISEPF = SQRT((RS(INMID,IR)-RS(INMID,IR+1))**2 +
     >                 (ZS(INMID,IR)-ZS(INMID,IR+1))**2 )

                  OSEPB = SQRT((RS(OUMID,IR-1)-RS(OUMID,IR))**2 +
     >                 (ZS(OUMID,IR-1)-ZS(OUMID,IR))**2 )
                  ISEPB = SQRT((RS(INMID,IR-1)-RS(INMID,IR))**2 +
     >                 (ZS(INMID,IR-1)-ZS(INMID,IR))**2 )

                  ODELN  = ((KNBS(OUMID,IR)-KNBS(OUMID,IR+1))/OSEPF+
     >                 (KNBS(OUMID,IR-1)-KNBS(OUMID,IR))/OSEPB)/2.0

                  ODELTI = ((KTIBS(OUMID,IR)-KTIBS(OUMID,IR+1))/OSEPF+
     >                 (KTIBS(OUMID,IR-1)-KTIBS(OUMID,IR))/OSEPB)/2.0

                  ODELTE = ((KTEBS(OUMID,IR)-KTEBS(OUMID,IR+1))/OSEPF+
     >                 (KTEBS(OUMID,IR-1)-KTEBS(OUMID,IR))/OSEPB)/2.0
c     
                  IDELN  = ( (KNBS(inMID,IR)-KNBS(inMID,IR+1)) /ISEPF+
     >                 (KNBS(inMID,IR-1)-KNBS(inMID,IR)) /ISEPB )/2.0

                  IDELTI = ((KTIBS(inMID,IR)-KTIBS(inMID,IR+1))/ISEPF+
     >                 (KTIBS(inMID,IR-1)-KTIBS(inMID,IR))/ISEPB)/2.0

                  IDELTE = ((KTEBS(inMID,IR)-KTEBS(inMID,IR+1))/ISEPF+
     >                 (KTEBS(inMID,IR-1)-KTEBS(inMID,IR))/ISEPB)/2.0

                  OSEP = (OSEPF+OSEPB )/2.0
                  ISEP = (ISEPF+ISEPB )/2.0

c     write (6,*) 'TESTI:',ir,inmid,ideln,isepf,isepb,iarea
c     write (6,*) 'TESTI:',ir,knbs(inmid,ir-1),knbs(inmid,ir),
c     >              knbs(inmid,ir+1)
c     write (6,*) 'TESTI:',rs(inmid,ir-1),zs(inmid,ir-1),
c     >               rs(inmid,ir),
c     >               zs(inmid,ir),rs(inmid,ir+1),zs(inmid,ir+1)
c     write (6,*) 'TESTI:',( (KNBS(inMID,IR)-KNBS(inMID,IR+1))
c     >             /ISEPF+
c     >             (KNBS(inMID,IR-1)-KNBS(inMID,IR)) /ISEPB )/2.0
c     write (6,*) 'TESTO:',ir,oumid,odeln,osepf,osepb,oarea
c     write (6,*) 'TESTO:',ir,knbs(oumid,ir-1),knbs(oumid,ir),
c     >              knbs(oumid,ir+1)
c     write (6,*) 'TESTO:',rs(oumid,ir-1),zs(oumid,ir-1),
c     >               rs(oumid,ir),
c     >               zs(oumid,ir),rs(oumid,ir+1),zs(oumid,ir+1)
c     write (6,*) 'TESTO:',((KNBS(OUMID,IR)-KNBS(OUMID,IR+1))
c     >              /OSEPF+
c     >           (KNBS(OUMID,IR-1)-KNBS(OUMID,IR))/OSEPB)/2.0
c     
c---------------------------------------------------------------------
c     
               endif
c     
c     Density values
c     
               NI = knbs(inmid,ir)
               NO = knbs(oumid,ir)
               NAV = (NI+NO)/2.0
c     
c     Overall gradient - take average of inner/outer
c     
               SEP = (OSEP+ISEP)/2.0
c     
               DELN   = (ODELN  + IDELN)/2.0
               DELTE  = (ODELTE + IDELTE)/2.0
               DELTI  = (ODELTI + IDELTI)/2.0
c     
               TOTFLX=TOTFLX-NETFLX(IR)
               TOTQI=TOTQI-HEATI(IR)
               TOTQE=TOTQE-HEATE(IR)
c     
               ITOTFLX= ITOTFLX -INETFLX(IR)
               ITOTQI = ITOTQI  -INHEATI(IR)
               ITOTQE = ITOTQE  -INHEATE(IR)
c     
               OTOTFLX= OTOTFLX -ONETFLX(IR)
               OTOTQI = OTOTQI  -OUHEATI(IR)
               OTOTQE = OTOTQE  -OUHEATE(IR)
c     
c     TOTflx
c     

C     
C---- Calculate the Inboard/Outboard Dperp Values
C     
c     
c     Total - using averaged mid-plane gradients
c     
               DPERP(IR)    = (TOTFLX/AREA)/DELN
               CHIPERPI(IR) = (TOTQI/AREA)/(DELTI*ech)/NAV
               CHIPERPE(IR) = (TOTQE/AREA)/(DELTE*ech)/nav
c     
c     Outer
c     
               ODPERP(IR)    = (OTOTFLX/OAREA)/ODELN
               OCHIPERPI(IR) = (OTOTQI/OAREA)/(ODELTI*ech)/no
               OCHIPERPE(IR) = (OTOTQE/OAREA)/(ODELTE*ech)/no
c     
c     Inner
c     
               IDPERP(IR)    = (ITOTFLX/IAREA)/IDELN
               ICHIPERPI(IR) = (ITOTQI/IAREA)/(IDELTI*ech)/ni
               ICHIPERPE(IR) = (ITOTQE/IAREA)/(IDELTE*ech)/ni
c     
c     DPMETHOD
c     

            elseif (dpmethod.eq.1) then
c     
c     Method #2
c     

               do ik = ouend,inend
c     
                  if (dprcopt.eq.0) then
                     fact = 1.0
                  else
                     fact = rs(ik,ir)/r0
                  endif
c     
                  if (ik.eq.1) then

                     surf(ik) = (SQRT( (RS(IK,IR)-RS(IK+1,IR))**2 +
     >                    (ZS(IK,IR)-ZS(IK+1,IR))**2 )/ 2.0
     >                    + SQRT( (RS(IK,IR)-RP(idds(ir,2)))**2 +
     >                    (ZS(IK,IR)-ZP(idds(ir,2)))**2 )
     >                    ) * fact

                  elseif (ik.eq.nks(ir)) then

                     surf(ik) = (SQRT( (RS(IK,IR)-RP(idds(ir,1)))**2 +
     >                    (ZS(IK,IR)-ZP(idds(ir,1)))**2 )
     >                    + SQRT( (RS(IK,IR)-RS(IK-1,IR))**2 +
     >                    (ZS(IK,IR)-ZS(IK-1,IR))**2 )/ 2.0
     >                    ) * fact

                  else

                     surf(ik) = (SQRT( (RS(IK,IR)-RS(IK+1,IR))**2 +
     >                    (ZS(IK,IR)-ZS(IK+1,IR))**2 )/ 2.0
     >                    + SQRT( (RS(IK,IR)-RS(IK-1,IR))**2 +
     >                    (ZS(IK,IR)-ZS(IK-1,IR))**2 )/ 2.0
     >                    ) * fact

                  endif


c     
                  if (ir.eq.irsep) then
c     
                     SEP = SQRT((RS(ik,IR)-RS(ik,IR+1))**2 +
     >                    (ZS(ik,IR)-ZS(ik,IR+1))**2 )

                     DELN  = (KNBS(ik,IR)-KNBS(ik,IR+1))/SEP

                     DELTI = (KTIBS(ik,IR)-KTIBS(ik,IR+1))/SEP
                     DELTE = (KTEBS(ik,IR)-KTEBS(ik,IR+1))/SEP
c     
                  elseif (ir.eq.irwall-1) then
c     
                     SEP = SQRT((RS(ik,IR)-RS(ik,IR-1))**2 +
     >                    (ZS(ik,IR)-ZS(ik,IR-1))**2 )

                     DELN  = (KNBS(ik,IR-1)-KNBS(ik,IR))/SEP

                     DELTI = (KTIBS(ik,IR-1)-KTIBS(ik,IR))/SEP
                     DELTE = (KTEBS(ik,IR-1)-KTEBS(ik,IR))/SEP
c     
                  else

                     SEPF = SQRT((RS(ik,IR)-RS(ik,IR+1))**2 +
     >                    (ZS(ik,IR)-ZS(ik,IR+1))**2 )

                     SEPB = SQRT((RS(ik,IR)-RS(ik,IR-1))**2 +
     >                    (ZS(ik,IR)-ZS(ik,IR-1))**2 )

                     DELN  = ((KNBS(ik,IR)-KNBS(ik,IR+1))/SEPF+
     >                    (KNBS(ik,IR-1)-KNBS(ik,IR))/SEPB)/2.0

                     DELTI = ((KTIBS(ik,IR)-KTIBS(ik,IR+1))/SEPF+
     >                    (KTIBS(ik,IR-1)-KTIBS(ik,IR))/SEPB)/2.0

                     DELTE = ((KTEBS(ik,IR)-KTEBS(ik,IR+1))/SEPF+
     >                    (KTEBS(ik,IR-1)-KTEBS(ik,IR))/SEPB)/2.0
c     
                  endif
c     
                  deldp(ik) = surf(ik) * deln
                  delqe(ik) = surf(ik) * knbs(ik,ir) * ech * delte
                  delqi(ik) = surf(ik) * knbs(ik,ir) * ech * delti
c     
c     End of do ik
c     
               end do
c     
c     Now we have all the single cell fluxes
c     in the del variables - need to sum over them
c     and equate them to the remaining outflux to
c     estimate the transport coefficients for the
c     ring.
c     
               dptoto = 0.0
               qetoto = 0.0
               qitoto = 0.0
c     
               do ik = ouend,ikmid -1
                  dptoto = dptoto + deldp(ik)
                  qetoto = qetoto + delqe(ik)
                  qitoto = qitoto + delqi(ik)
               end do
c     
               dptoti = 0.0
               qetoti = 0.0
               qitoti = 0.0
c     
               do ik = ikmid,inend
                  dptoti = dptoti + deldp(ik)
                  qetoti = qetoti + delqe(ik)
                  qitoti = qitoti + delqi(ik)
               end do
c     
c     Take ratio with totals to extract dperp, xperp ...
c     
               TOTFLX=TOTFLX-NETFLX(IR)
               TOTQI=TOTQI-HEATI(IR)
               TOTQE=TOTQE-HEATE(IR)
c     
               ITOTFLX= ITOTFLX -INETFLX(IR)
               ITOTQI = ITOTQI  -INHEATI(IR)
               ITOTQE = ITOTQE  -INHEATE(IR)
c     
               OTOTFLX= OTOTFLX -ONETFLX(IR)
               OTOTQI = OTOTQI  -OUHEATI(IR)
               OTOTQE = OTOTQE  -OUHEATE(IR)
c     
C     
C---- Calculate the Inboard/Outboard Dperp Values
C     
c     
c     Total - using averaged mid-plane gradients
c     
               DPERP(IR)    = TOTFLX/(dptoto+dptoti)
               CHIPERPI(IR) = TOTQI/(qitoto+qitoti)
               CHIPERPE(IR) = TOTQE/(qetoto+qetoti)
c     
c     Outer
c     
               ODPERP(IR)    = OTOTFLX/dptoto
               OCHIPERPI(IR) = OTOTQI/qitoto
               OCHIPERPE(IR) = OTOTQE/qetoto
c     
c     Inner
c     
               IDPERP(IR)    = ITOTFLX/dptoti
               ICHIPERPI(IR) = ITOTQI/qitoti
               ICHIPERPE(IR) = ITOTQE/qetoti
c     

            endif

c     
c     Dperp
c     
            write(6,*) 'Indices: IKMID=',ikmid,
     >           ' IKINMID=',inmid,' IKOUMID=',oumid
            write(6,*) 'Z0 = ',z0
            write(6,500)
            write(6,400) 'Dperp '//inner//':',ir,idperp(ir),itotflx,
     >           ideln,iarea,isep,
     >           influx(ir),iionis(ir),
     >           inetflx(ir),knbs(inmid,ir),isepb,isepf
            write(6,400) 'Dperp '//outer//':',ir,odperp(ir),ototflx,
     >           odeln,oarea,osep,
     >           ouflux(ir),oionis(ir),
     >           onetflx(ir),knbs(oumid,ir),osepb,osepf
            write(6,400) 'Dperp Total:',ir,dperp(ir),totflx,deln,area,
     >           sep,flux(ir),ionis(ir),netflx(ir),nav

c     
c     Xperpi
c     
            write(6,700)
            write(6,400) 'Xperpi '//Inner,ir,ichiperpi(ir),itotqi,
     >           idelti,iarea,isep,
     >           inheati(ir),ktibs(inmid,ir),isepb,isepf
            write(6,400) 'Xperpi '//Outer,ir,ochiperpi(ir),ototqi,
     >           odelti,oarea,osep,
     >           ouheati(ir),ktibs(oumid,ir),osepb,osepf
            write(6,400) 'Xperpi Total',ir,chiperpi(ir),totqi,delti,
     >           area,sep,heati(ir)
c     
c--------------------------------------------------------------------
c     
c     Xperpe
c     

            write(6,900)
            write(6,400) 'Xperpe '//Inner,ir,ichiperpe(ir),itotqe,
     >           idelte,iarea,isep,
     >           inheate(ir),ktebs(inmid,ir),isepb,isepf
            write(6,400) 'Xperpe '//Outer,ir,ochiperpe(ir),ototqe,
     >           odelte,oarea,osep,
     >           ouheate(ir),ktebs(oumid,ir),osepb,osepf
            write(6,400) 'Xperpe Total',ir,chiperpe(ir),totqe,delte,
     >           area,sep,heate(ir)




c     
 280     CONTINUE
c     
c     
c     
      elseif (dpmethod.eq.2) then
c     
c     Calculate the gradients everywhere in the SOL
c     Te, Ti, N
c     
c     
c     First Calculate Area elements - for whole grid
c     
         do ir = irsep,irwall-1
c     
            do ik = 1,nks(ir)
c     
c     Calculate Area elements
c     
               if (dprcopt.eq.0) then
                  fact = 1.0
               else
                  fact = rs(ik,ir)/r0
               endif
c     
c     Calculate the surface area based on cell-centres
c     
               if (dparea.eq.0) then
c     
                  if (ik.eq.1) then

                     surfa(ik,ir) = (kps(ik,ir) + kps(ik+1,ir)) / 2.0
     >                    * fact

                  elseif (ik.eq.nks(ir)) then

                     surfa(ik,ir) = (kpmaxs(ir)
     >                    -  0.5 *(kps(ik,ir)+kps(ik-1,ir)))
     >                    * fact
                  else

                     surfa(ik,ir) = 0.5 *(kps(ik+1,ir)- kps(ik-1,ir))
     >                    * fact
                  endif
c     
c     Calculate the surface area based on the outward cell boundary.
c     
               elseif (dparea.eq.1) then
c     
                  in = korpg(ik,ir)
c     
                  if (in.eq.0) then
                     surfa(ik,ir) = 0.0
                  else
                     surfa(ik,ir) = sqrt((rvertp(2,in)-rvertp(3,in))**2+
     >                    (zvertp(2,in)-zvertp(3,in))**2)
     >                    * fact
                  endif
c     
c     
c     Calculate the surface area based on polygon edges
c     
               elseif (dparea.eq.2) then
c     
                  surfa(ik,ir) = (kpb(ik,ir) - kpb(ik-1,ir))
     >                 * fact
c     
               endif
c     
            end do
         end do
c     
c     
c     Calculate actual gradients - for whole grid - no matter what
c     options are specified.
c     
c     
c     DPNAV Greater than zero - use normal derivative finding estimates
c     
         if (dpnav.ge.0) then

            do ir = irsep,irwall-1
c     
               do ik = 1,nks(ir)
c     
c     Calculate gradients - for E2d if data is also given.
c     
c     Calculate using Method 0 - equally weighted averaging
c     
                  if (dpnav.eq.0) then
c     

                     if (ir.eq.irsep) then
c     
                        SEP = koutds(ik,ir)
                        ikout = ikouts(ik,ir)
                        irout = irouts(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ik,IR)
     >                                     -KNBS(ikout,irout))/SEP
                           dgradn(ik,ir) = gradn(ik,ir)
                           gradTI(ik,ir) = (KTIBS(ik,IR)
     >                                     -KTIBS(ikout,IRout))/SEP
                           gradtE(ik,ir) = (KTEBS(ik,IR)
     >                                     -KTEBS(ikout,IRout))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ik,IR)
     >                             -e2dNBS(ikout,irout))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ik,IR)
     >                             -e2dTIBS(ikout,IRout))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ik,IR)
     >                             -e2dTEBS(ikout,IRout))/SEP
                           endif
c     
                        else
                           gradn(ik,ir)  = 0.0
                           dgradn(ik,ir) = gradn(ik,ir)
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
c     
                     elseif (ir.eq.irwall-1) then
c     
                        SEP = kinds(ik,ir)
                        ikin = ikins(ik,ir)
                        irin = irins(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                                     -KNBS(ik,IR))/SEP
                           dgradn(ik,ir) = gradn(ik,ir)
                           gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                                      -KTIBS(ik,IR))/SEP
                           gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                                     -KTEBS(ik,IR))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                             -e2dNBS(ik,ir))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                             -e2dTIBS(ik,IR))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                             -e2dTEBS(ik,IR))/SEP
                           endif
                        else
                           gradn(ik,ir)  = 0.0
                           dgradn(ik,ir) = gradn(ik,ir)
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
c     
                     else
c     
                        SEPF = koutds(ik,ir)
                        ikout = ikouts(ik,ir)
                        irout = irouts(ik,ir)
c     
                        SEPB = kinds(ik,ir)
                        ikin = ikins(ik,ir)
                        irin = irins(ik,ir)
c     
                        if (sepf.ne.0.0.and.sepb.ne.0.0) then
                           gradn(ik,ir)  = ((KNBS(ikin,IRin)
     >                                      -KNBS(ik,IR))/SEPB+
     >                            (KNBS(ik,IR)-KNBS(ikout,irout))
     >                              /SEPF)/2.0
                           dgradn(ik,ir)  = ((KNBS(ikin,IRin)
     >                                       -KNBS(ik,IR))/SEPB-
     >                          (KNBS(ik,IR)-KNBS(ikout,irout))/SEPF)
                           gradTI(ik,ir) = ((KTIBS(ikin,IRin)
     >                                      -KTIBS(ik,IR))/SEPB+
     >                          (KTIBS(ik,IR)-KTIBS(ikout,IRout))
     >                             /SEPF)/2.0
                           gradtE(ik,ir) = ((KTEBS(ikin,IRin)
     >                                      -KTEBS(ik,IR))/SEPB+
     >                          (KTEBS(ik,IR)-KTEBS(ikout,IRout))
     >                             /SEPF)/2.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = ((e2dNBS(ikin,IRin)
     >                             -e2dNBS(ik,IR))/SEPB+
     >                             (e2dNBS(ik,IR)
     >                             -e2dNBS(ikout,irout))/SEPF)/2.0
                              e2dgradTI(ik,ir) = ((e2dTIBS(ikin,IRin)
     >                             -e2dTIBS(ik,IR))/SEPB+
     >                             (e2dTIBS(ik,IR)
     >                             -e2dTIBS(ikout,IRout))/SEPF)/2.0
                              e2dgradtE(ik,ir) = ((e2dTEBS(ikin,IRin)
     >                             -e2dTEBS(ik,IR))/SEPB+
     >                             (e2dTEBS(ik,IR)
     >                             -e2dTEBS(ikout,IRout))/SEPF)/2.0
                           endif
                        elseif (sepf.eq.0.0) then
                           if (sepb.eq.0.0) then
                              gradn(ik,ir)  = 0.0
                              dgradn(ik,ir) = gradn(ik,ir)
                              gradTI(ik,ir) = 0.0
                              gradtE(ik,ir) = 0.0
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = 0.0
                                 e2dgradTI(ik,ir) = 0.0
                                 e2dgradtE(ik,ir) = 0.0
                              endif
                           else
                              gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                                        -KNBS(ik,IR))/SEPB
                              dgradn(ik,ir) = gradn(ik,ir)
                              gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                                        -KTIBS(ik,IR))/SEPB
                              gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                                        -KTEBS(ik,IR))/SEPB
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                                -e2dNBS(ik,ir))/SEPb
                                 e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                                -e2dTIBS(ik,IR))/SEPb
                                 e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                                -e2dTEBS(ik,IR))/SEPb
                              endif
                           endif
                        elseif (sepb.eq.0.0) then
                           if (sepf.eq.0.0) then
                              gradn(ik,ir)  = 0.0
                              dgradn(ik,ir) = gradn(ik,ir)
                              gradTI(ik,ir) = 0.0
                              gradtE(ik,ir) = 0.0
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = 0.0
                                 e2dgradTI(ik,ir) = 0.0
                                 e2dgradtE(ik,ir) = 0.0
                              endif
                           else
                              gradn(ik,ir) =(KNBS(ik,IR)
     >                                      -KNBS(ikout,irout))/SEPF
                              dgradn(ik,ir) = gradn(ik,ir)
                              gradTI(ik,ir)=(KTIBS(ik,IR)
     >                                      -KTIBS(ikout,IRout))/SEPF
                              gradtE(ik,ir)=(KTEBS(ik,IR)
     >                                      -KTEBS(ikout,IRout))/SEPF
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = (e2dNBS(ik,IR)
     >                                -e2dNBS(ikout,irout))/SEPf
                                 e2dgradTI(ik,ir) = (e2dTIBS(ik,IR)
     >                                -e2dTIBS(ikout,IRout))/SEPf
                                 e2dgradtE(ik,ir) = (e2dTEBS(ik,IR)
     >                                -e2dTEBS(ikout,IRout))/SEPf
                              endif
                           endif
                        endif
                     endif
c     
c     Calculate using Method 1 - weighted by ring separation
c     
                  elseif (dpnav.eq.1) then

                     if (ir.eq.irsep) then
c     
                        SEP = koutds(ik,ir)
                        ikout = ikouts(ik,ir)
                        irout = irouts(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ik,IR)
     >                                     -KNBS(ikout,irout))/SEP
                           gradTI(ik,ir) = (KTIBS(ik,IR)
     >                                     -KTIBS(ikout,IRout))/SEP
                           gradtE(ik,ir) = (KTEBS(ik,IR)
     >                                     -KTEBS(ikout,IRout))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ik,IR)
     >                             -e2dNBS(ikout,irout))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ik,IR)
     >                             -e2dTIBS(ikout,IRout))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ik,IR)
     >                             -e2dTEBS(ikout,IRout))/SEP
                           endif
c     
                        else
                           gradn(ik,ir)  = 0.0
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
c     
                     elseif (ir.eq.irwall-1) then
c     
                        SEP = kinds(ik,ir)
                        ikin = ikins(ik,ir)
                        irin = irins(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                                     -KNBS(ik,IR))/SEP
                           gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                                     -KTIBS(ik,IR))/SEP
                           gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                                     -KTEBS(ik,IR))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                             -e2dNBS(ik,ir))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                             -e2dTIBS(ik,IR))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                             -e2dTEBS(ik,IR))/SEP
                           endif
                        else
                           gradn(ik,ir)  = 0.0
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
c     
                     else
c     
                        SEPF = koutds(ik,ir)
                        ikout = ikouts(ik,ir)
                        irout = irouts(ik,ir)
c     
                        SEPB = kinds(ik,ir)
                        ikin = ikins(ik,ir)
                        irin = irins(ik,ir)
c     
                        if (sepf.ne.0.0.and.sepb.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                          -KNBS(ikout,irout))/(sepb+sepf)
                           gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                          -KTIBS(ikout,IRout))/(sepb+sepf)
                           gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                          -KTEBS(ikout,IRout))/(sepb+sepf)
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                             -e2dNBS(ikout,irout))/(sepb+sepf)
                              e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                             -e2dTIBS(ikout,IRout))/(sepb+sepf)
                              e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                             -e2dTEBS(ikout,IRout))/(sepb+sepf)
                           endif
                        elseif (sepf.eq.0.0) then
                           if (sepb.eq.0.0) then
                              gradn(ik,ir)  = 0.0
                              gradTI(ik,ir) = 0.0
                              gradtE(ik,ir) = 0.0
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = 0.0
                                 e2dgradTI(ik,ir) = 0.0
                                 e2dgradtE(ik,ir) = 0.0
                              endif
                           else
                              gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                                        -KNBS(ik,IR))/SEPB
                              gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                                        -KTIBS(ik,IR))/SEPB
                              gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                                        -KTEBS(ik,IR))/SEPB
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                                -e2dNBS(ik,ir))/SEPb
                                 e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                                -e2dTIBS(ik,IR))/SEPb
                                 e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                                -e2dTEBS(ik,IR))/SEPb
                              endif
                           endif
                        elseif (sepb.eq.0.0) then
                           if (sepf.eq.0.0) then
                              gradn(ik,ir)  = 0.0
                              gradTI(ik,ir) = 0.0
                              gradtE(ik,ir) = 0.0
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = 0.0
                                 e2dgradTI(ik,ir) = 0.0
                                 e2dgradtE(ik,ir) = 0.0
                              endif
                           else
                              gradn(ik,ir) =(KNBS(ik,IR)
     >                                      -KNBS(ikout,irout))/SEPF
                              gradTI(ik,ir)=(KTIBS(ik,IR)
     >                                      -KTIBS(ikout,IRout))/SEPF
                              gradtE(ik,ir)=(KTEBS(ik,IR)
     >                                      -KTEBS(ikout,IRout))/SEPF
                              if (cre2d.eq.1.or.cre2d.eq.2) then
                                 e2dgradn(ik,ir)  = (e2dNBS(ik,IR)
     >                                -e2dNBS(ikout,irout))/SEPf
                                 e2dgradTI(ik,ir) = (e2dTIBS(ik,IR)
     >                                -e2dTIBS(ikout,IRout))/SEPf
                                 e2dgradtE(ik,ir) = (e2dTEBS(ik,IR)
     >                                -e2dTEBS(ikout,IRout))/SEPf
                              endif
                           endif
                        endif
                     endif
c     
c     Elseif (dpnav = 2)
c     
                  elseif (dpnav.eq.2) then
c     
c     Calculate gradient going forward to next ring ONLY
c     
c     
                     if (ir.eq.irwall-1) then
c     
                        SEP = kinds(ik,ir)
                        ikin = ikins(ik,ir)
                        irin = irins(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ikin,IRin)
     >                                      -KNBS(ik,IR))/SEP
                           gradTI(ik,ir) = (KTIBS(ikin,IRin)
     >                                      -KTIBS(ik,IR))/SEP
                           gradtE(ik,ir) = (KTEBS(ikin,IRin)
     >                                      -KTEBS(ik,IR))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ikin,IRin)
     >                             -e2dNBS(ik,ir))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ikin,IRin)
     >                             -e2dTIBS(ik,IR))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ikin,IRin)
     >                             -e2dTEBS(ik,IR))/SEP
                           endif
                        else
                           gradn(ik,ir)  = 0.0
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
c     
                     else
c     
                        SEP = koutds(ik,ir)
                        ikout = ikouts(ik,ir)
                        irout = irouts(ik,ir)
c     
                        if (sep.ne.0.0) then
                           gradn(ik,ir)  = (KNBS(ik,IR)
     >                                      -KNBS(ikout,irout))/SEP
                           gradTI(ik,ir) = (KTIBS(ik,IR)
     >                                      -KTIBS(ikout,IRout))/SEP
                           gradtE(ik,ir) = (KTEBS(ik,IR)
     >                                      -KTEBS(ikout,IRout))/SEP
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = (e2dNBS(ik,IR)
     >                             -e2dNBS(ikout,irout))/SEP
                              e2dgradTI(ik,ir) = (e2dTIBS(ik,IR)
     >                             -e2dTIBS(ikout,IRout))/SEP
                              e2dgradtE(ik,ir) = (e2dTEBS(ik,IR)
     >                             -e2dTEBS(ikout,IRout))/SEP
                           endif
c     
                        else
                           gradn(ik,ir)  = 0.0
                           gradTI(ik,ir) = 0.0
                           gradtE(ik,ir) = 0.0
                           if (cre2d.eq.1.or.cre2d.eq.2) then
                              e2dgradn(ik,ir)  = 0.0
                              e2dgradTI(ik,ir) = 0.0
                              e2dgradtE(ik,ir) = 0.0
                           endif
                        endif
                     endif
c     
c     Endif for dpnav
c     
                  endif
c     
c     End of first ik,ir loops
c     
c     write (6,*) 'grad:',ir,ik,surfa(ik,ir),gradte(ik,ir),
c     >              gradti(ik,ir),gradn(ik,ir)


               end do
            end do
c     
c     End of DPNAV >= 0
c     
         elseif (dpnav.lt.0) then
c     
c     Need to take each "ik sheaf" of elements as one group
c     This requires an odd looping construct.
c     
c     The following assumes that all of the rings in the
c     main SOL have the same number of IK elements and as
c     such uses nks(irsep) as the base value.
c     
            do ik = 1,nks(irsep)
c     
c     Load arrays for passing to the spline routine
c     
               N1 = (irwall-1) -irsep + 1
               dist = 0.0
c     
c     Density
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  x1(in) = dist
                  f1(in) = knbs(ik,ir)
                  dist = dist + koutds(ik,ir)
               end do
c     
               CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  gradn(ik,ir) = -fdash1(in)
               end do
c     
c     Te
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  f1(in) = ktebs(ik,ir)
               end do
c     
               CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  gradte(ik,ir) = -fdash1(in)
               end do
c     
c     Ti
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  f1(in) = ktibs(ik,ir)
               end do
c     
               CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
               do ir = irsep,irwall-1
                  in = ir - irsep +1
                  gradti(ik,ir) = -fdash1(in)
               end do
c     
c     Calculate gradients for E2D data if available
c     
               if (cre2d.eq.1.or.cre2d.eq.2) then
c     
c     E2D Density
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     f1(in) = e2dnbs(ik,ir)
                  end do
c     
                  CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     e2dgradn(ik,ir) = -fdash1(in)
                  end do
c     
c     E2D Te
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     f1(in) = e2dtebs(ik,ir)
                  end do
c     
                  CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     e2dgradte(ik,ir) = -fdash1(in)
                  end do
c     
c     E2D Ti
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     f1(in) = e2dtibs(ik,ir)
                  end do
c     
                  CALL TB04A (N1,X1,F1,FDASH1,WORK)
c     
                  do ir = irsep,irwall-1
                     in = ir - irsep +1
                     e2dgradti(ik,ir) = -fdash1(in)
                  end do
c     
               endif

            end do
c     
c     Endif for DPNAV
c     
         endif
c     
c     Modify the Gradients if the non-orthogonal correction has
c     been turned on. (This is a first order correction based on
c     the individual cell non-orthogonality as an estimate of
c     the amount of change required in the calculated gradient.)
c     
         if (dporth.eq.1) then
            do ir = irsep,irwall-1
               do ik = 1,nks(ir)
                  gradn(ik,ir) = gradn(ik,ir) / sinalph(ik,ir)
                  gradte(ik,ir) = gradte(ik,ir) / sinalph(ik,ir)
                  gradti(ik,ir) = gradti(ik,ir) / sinalph(ik,ir)
c     
                  if (cre2d.eq.1.or.cre2d.eq.2) then
                     e2dgradn(ik,ir) = e2dgradn(ik,ir)/sinalph(ik,ir)
                     e2dgradte(ik,ir)= e2dgradte(ik,ir)/sinalph(ik,ir)
                     e2dgradti(ik,ir)= e2dgradti(ik,ir)/sinalph(ik,ir)
                  endif
c     
               end do
            end do
c     
c     End of DPORTH
c     
         endif
c     
c     Smooth the gradients radially by averaging.
c     
         if (dpsmooth.gt.0) then
c     
            do ir = irsep,irwall -1
c     
c     
c     Loop over all radial elements
c     
               do ik = 1,nks(ir)
c     
                  totgradn = 0.0
                  totgradte = 0.0
                  totgradti = 0.0
c     
c     Loop over range of averaging
c     
                  avmod = 0
c     
                  do in = ir-dpsmooth,ir+dpsmooth

                     if (in.lt.irsep) then

c     
c     totgradn  = totgradn + gradn(ik,irsep)
c     totgradte = totgradte+ gradte(ik,irsep)
c     totgradti = totgradti+ gradti(ik,irsep)
c     
                        avmod = avmod + 1

                     elseif (in.gt.irwall-1) then

c     
c     totgradn  = totgradn + gradn(ik,irwall-1)
c     totgradte = totgradte+ gradte(ik,irwall-1)
c     totgradti = totgradti+ gradti(ik,irwall-1)
c     
                        avmod = avmod + 1

                     else

                        totgradn  = totgradn + gradn(ik,ir)
                        totgradte = totgradte+ gradte(ik,ir)
                        totgradti = totgradti+ gradti(ik,ir)

                     endif
                  end do

                  gradn(ik,ir) = totgradn / (2.0*dpsmooth +1.0-avmod)
                  gradte(ik,ir) = totgradte / (2.0*dpsmooth +1.0-avmod)
                  gradti(ik,ir) = totgradti / (2.0*dpsmooth +1.0-avmod)

               end do

            end do
c     
c     End for DPSMOOTH
c     
         endif
c     
c     
c     Calculate the Dperp for each ring
c     
c     
c     
c     If OUTER ring losses are to be included - calculate this first
c     
c     
         delouti = 0.0
         delouto = 0.0
c     
         if (dpouter.eq.1) then
c     
            if (dpsuml.eq.0) then
               ikstart = ouend
               ikend = inend
            elseif (dpsuml.eq.1) then
               ikstart = 1
               ikend = nks(ir)
            endif
c     
            do ik = ikstart,ikend
               if (kss(ik,irwall-1).lt.(0.5*ksmaxs(irwall-1))) then
                  delouto = delouto + surfa(ik,irwall-1)*
     >                                      gradn(ik,irwall-1)
               else
                  delouti = delouti + surfa(ik,irwall-1)
     >                                      *gradn(ik,irwall-1)
               endif
            end do
         endif
c     
         write (6,*) 'delout:',delouto,delouti,delouto+delouti
c     
         do ir = irsep,irwall-1
c     
            dptoto = 0.0
            dptoti = 0.0
c     
            if (dpsuml.eq.0) then
               ikstart = ouend
               ikend = inend
            elseif (dpsuml.eq.1) then
               ikstart = 1
               ikend = nks(ir)
            endif
c     
            do ik = ikstart,ikend
c     
               deldp(ik) = surfa(ik,ir) * gradn(ik,ir)
c     
               if (ik.lt.ikmid) then
                  dptoto = dptoto + deldp(ik)
               else
                  dptoti = dptoti + deldp(ik)
               endif
c     
            end do
c     
            write (6,*) 'dptot:',ir,dptoto,dptoti,dptoto+dptoti
c     
c     This implicitly assumes that the Dperp on the outer ring is
c     equal to the Dperp on the current ring.
c     
            if (usedpav) then
               if (dpouter.eq.0) then
                  DPERP(IR)  = cfgam(ir)/(dptoto+dptoti)
                  ODPERP(IR) = ocfgam(ir)/dptoto
                  IDPERP(IR) = icfgam(ir)/dptoti
               elseif (dpouter.eq.1) then
                  DPERP(IR)  = (cfgam(ir)+dpav*(delouti+delouto))
     >                 /(dptoto+dptoti)
                  ODPERP(IR) = (ocfgam(ir)+dpavo*delouto)/(dptoto)
                  IDPERP(IR) = (icfgam(ir)+dpavi*delouti)/(dptoti)
               endif
            else
               if (dpouter.eq.0) then
                  DPERP(IR)  = cfgam(ir)/(dptoto+dptoti)
                  ODPERP(IR) = ocfgam(ir)/dptoto
                  IDPERP(IR) = icfgam(ir)/dptoti
               elseif (dpouter.eq.1) then
c
c                 jdemod - for ir=irwall-1 - the denominator in these
c                          calculations is exactly zero.
c                        - to avoid this - for ir=irwall-1 - assign
c                          dperp(irwall-1) = dperp(irwall-2)
c                        - this fix needs checking to make sure it is
c                          logically correct                         
c
                  if (ir.ne.irwall-1) then 
                     DPERP(IR)  = cfgam(ir)
     >                           /(dptoto+dptoti-(delouti+delouto))
                     ODPERP(IR) = ocfgam(ir)/(dptoto-delouto)
                     IDPERP(IR) = icfgam(ir)/(dptoti-delouti)
                  else
                     DPERP(IR)  = dperp(ir-1)
                     ODPERP(IR) = odperp(ir-1)
                     IDPERP(IR) = idperp(ir-1)
                  endif

               endif
            endif

c     
c     Calculate dpav if the option is set and there are enough rings
c     
            if (dpavopt.gt.0.and.ir.eq.(irsep+irskip+3)) then
               irtmp = irsep + irskip
               dpav = 0.25*(dperp(irtmp) + dperp(irtmp+1)
     >              +dperp(irtmp+2) + dperp(irtmp+3))
               dpavo = 0.25*(odperp(irtmp) + odperp(irtmp+1)
     >              +odperp(irtmp+2) + odperp(irtmp+3))
               dpavi = 0.25*(idperp(irtmp) + idperp(irtmp+1)
     >              +idperp(irtmp+2) + idperp(irtmp+3))
               usedpav = .true.

               if (cprint.eq.1.or.cprint.eq.9) then
                  write (6,*) 'DPERP AV: rings ',irtmp,' to ',
     >                 irtmp+3,' = ',dpav

               endif

            endif

c     
c     Record values for each ring
c     
            flxval(ir,1,1) = dptoto
            flxval(ir,1,2) = dptoti
            flxval(ir,1,3) = dptoto + dptoti
            flxval(ir,2,1) = delouto
            flxval(ir,2,2) = delouti
            flxval(ir,2,3) = delouti + delouto
            flxval(ir,3,1) = ocfgam(ir)
            flxval(ir,3,2) = icfgam(ir)
            flxval(ir,3,3) = cfgam(ir)

c     
         end do
c     
c     Now calculate the Xperps - add in cross-field convection if
c     turned on = 5/2 kT Gam(Perp) = 5/2 kT Dperp dn/dr
c     
c     
c     If the outer ring loss option is turned on - calculate the
c     OUTER ring conductive loss term assuming the Xperp for the
c     outer ring equals the Xperp for the ring being calculated.
c     Also - if the cross-field convective losses are turned
c     on calculate these for the last ring.
c     
c     
c     Zero out the values in case the option is not in use and
c     for initialization purposes in case they are.
c     
         qeouti= 0.0
         qeouto= 0.0
         qiouti= 0.0
         qiouto= 0.0
c     
         qeoutconvi= 0.0
         qeoutconvo= 0.0
         qioutconvi= 0.0
         qioutconvo= 0.0
c     
         if (dpouter.eq.1) then
c     
            if (dpsuml.eq.0) then
               ikstart = ouend
               ikend = inend
            elseif (dpsuml.eq.1) then
               ikstart = 1
               ikend = nks(ir)
            endif
c     
            ir = irwall - 1
            ird = irwall - 2
c     
            if (dpxpratio.gt.0.0) then
               dptmp = dpxpratio
            elseif (dpxpratio.eq.-1.0) then
               dptmp = cdperp
            elseif (usedpav) then
               dptmp = dpav
            else
               dptmp = dperp(ird)
            endif
c     
            do ik = ikstart,ikend
c     
               if (ik.lt.ikmid) then
                  qeouto = qeouto
     >                 + surfa(ik,ir)*knbs(ik,ir)*ech*gradte(ik,ir)
                  qiouto = qiouto
     >                 + surfa(ik,ir)*knbs(ik,ir)*ech*gradti(ik,ir)

               else
                  qeouti = qeouti
     >                 + surfa(ik,ir)*knbs(ik,ir)*ech*gradte(ik,ir)
                  qiouti = qiouti
     >                 + surfa(ik,ir)*knbs(ik,ir)*ech*gradti(ik,ir)
               endif
c     
               if (dpconv.eq.1) then
                  if (ik.lt.ikmid) then
                     qeoutconvo = qeoutconvo
     >                    +surfa(ik,ir)*2.5*dptmp
     >                    *ktebs(ik,ir)*ech*gradn(ik,ir)
                     qioutconvo = qioutconvo
     >                    +surfa(ik,ir)*2.5*dptmp
     >                    *ktibs(ik,ir)*ech*gradn(ik,ir)

                  else
                     qeoutconvi = qeoutconvi
     >                    +surfa(ik,ir)*2.5*dptmp*ktebs(ik,ir)
     >                    *ech*gradn(ik,ir)
                     qioutconvi = qioutconvi
     >                    +surfa(ik,ir)*2.5*dptmp*ktibs(ik,ir)
     >                    *ech*gradn(ik,ir)
                  endif
               endif
c     

            end do

         endif
c     
         if (cprint.eq.1.or.cprint.eq.9) then
            write (6,*) 'qeout:', qeouto,qeouti,qiouto,qiouti
            write (6,*) 'qeoutconv:', qeoutconvo,qeoutconvi,
     >           qioutconvo,qioutconvi
         endif
c     
c     
c     Calculate Xperp values
c     
         do ir = irsep, irwall -1
c     
c     Zero conductive loss terms
c     
            qetoti = 0.0
            qitoti = 0.0
            qetoto = 0.0
            qitoto = 0.0
c     
c     
c     Zero convective loss terms
c     
            qeconvi = 0.0
            qeconvo = 0.0
            qiconvi = 0.0
            qiconvo = 0.0
c     
            if (dpsuml.eq.0) then
               ikstart = ouend
               ikend = inend
            elseif (dpsuml.eq.1) then
               ikstart = 1
               ikend = nks(ir)
            endif
c     
            if (dpxpratio.gt.0.0) then
               dptmp = dpxpratio
            elseif (dpxpratio.eq.-1.0) then
               dptmp = cdperp
            elseif (usedpav) then
               dptmp = dpav
            else
               dptmp = dperp(ir)
            endif
c     
            do ik = ikstart,ikend
c     
               delqe(ik) = surfa(ik,ir) * knbs(ik,ir) 
     >                                      * ech * gradte(ik,ir)
               delqi(ik) = surfa(ik,ir) * knbs(ik,ir) 
     >                                      * ech * gradti(ik,ir)
               qecond(ik,ir) = delqe(ik)
               qicond(ik,ir) = delqi(ik)
c     
               if (ik.lt.ikmid) then
                  qetoto = qetoto + delqe(ik)
                  qitoto = qitoto + delqi(ik)
               else
                  qetoti = qetoti + delqe(ik)
                  qitoti = qitoti + delqi(ik)
               endif
c     
               if (dpconv.eq.1) then
c     
c     if (dpxpratio.gt.0.0) then
c     

                  qeconv(ik,ir) = surfa(ik,ir)*2.5*dptmp
     >                 *ktebs(ik,ir)*ech*gradn(ik,ir)
                  qiconv(ik,ir) = surfa(ik,ir)*2.5*dptmp
     >                 *ktibs(ik,ir)*ech*gradn(ik,ir)

c     
c     else
c     
c     qeconv(ik,ir) = surfa(ik,ir)*2.5*dperp(ir)
c     >                      *ktebs(ik,ir)*ech*gradn(ik,ir)
c     qiconv(ik,ir) = surfa(ik,ir)*2.5*dperp(ir)
c     >                      *ktibs(ik,ir)*ech*gradn(ik,ir)
c     endif
c     
c     
                  if (ik.lt.ikmid) then
                     qeconvo = qeconvo + qeconv(ik,ir)
                     qiconvo = qiconvo + qiconv(ik,ir)
                  else
                     qeconvi = qeconvi + qeconv(ik,ir)
                     qiconvi = qiconvi + qiconv(ik,ir)
                  endif
               endif

            end do
c     
            if (cprint.eq.1.or.cprint.eq.9) then
               write (6,*) 'q:',ir,qetoti,qitoti,qetoto,qitoto
               write (6,*) 'conv:',ir,qeconvi,qeconvo,
     >              qiconvi,qiconvo
            endif
c     
c     Calculate the Xperps
c     
            if (dpxpratio.gt.0.0) then


               qitot = (qitoto+qitoti-(qiouti+qiouto)
     >              +(qiconvi+qiconvo)
     >              -(qioutconvi+qioutconvo))

               if (qitot.ne.0.0) then
                  CHIPERPI(IR) = cfqi(ir)/qitot
c     >              / (qitoto+qitoti-(qiouti+qiouto)
c     >              +(qiconvi+qiconvo)
c     >              -(qioutconvi+qioutconvo))
               else
                  chiperpi(ir) = 0.0
               endif
c     
               qetot = (qetoto+qetoti-(qeouti+qeouto)
     >              +(qeconvi+qeconvo)
     >              -(qeoutconvi+qeoutconvo))

               
               if (qetot.ne.0.0) then 
                  CHIPERPE(IR) = cfqe(ir)/qetot
c     >              / (qetoto+qetoti-(qeouti+qeouto)
c     >              +(qeconvi+qeconvo)
c     >              -(qeoutconvi+qeoutconvo))
               else
                  chiperpe(ir)=0.0
               endif
c     


               if ((qetot+qitot).ne.0.0) then 
                  xperpt(ir) =  (cfqi(ir)+cfqe(ir))/(qetot+qitot)
c     >              / (qitoto+qitoti+qetoto+qetoti
c     >              -(qiouti+qiouto)-(qeouti+qeouto)
c     >              +(qiconvi+qiconvo)-(qioutconvi+qioutconvo)
c     >              +(qeconvi+qeconvo)-(qeoutconvi+qeoutconvo))
               else
                  xperpt(ir) = 0.0
               endif
c     
c     Outer
c     
               if ((qitoto-qiouto+qiconvo-qioutconvo).ne.0.0) then 
                  OCHIPERPI(IR) = ocfqi(ir)
     >                 / (qitoto-qiouto+qiconvo-qioutconvo)
               else
                  OCHIPERPI(IR) = 0.0
               endif

c     
               if ((qetoto-qeouto+qeconvo-qeoutconvo).ne.0.0) then 
                  OCHIPERPE(IR) = ocfqe(ir)
     >                 / (qetoto-qeouto+qeconvo-qeoutconvo)
               else
                  OCHIPERPE(IR) = 0.0
               endif

c     
               if ((qitoto+qetoto
     >              -(qiouto)-(qeouto)
     >              +(qiconvo)-(qioutconvo)
     >              +(qeconvo)-(qeoutconvo)).ne.0.0) then 
                  oxperpt(ir) =  (ocfqi(ir)+ocfqe(ir))
     >                 / (qitoto+qetoto
     >                 -(qiouto)-(qeouto)
     >                 +(qiconvo)-(qioutconvo)
     >                 +(qeconvo)-(qeoutconvo))
               else
                  oxperpt(ir) =  0.0
               endif

c     
c     Inner
c     
               if ((qitoti-qiouti+qiconvi-qioutconvi).ne.0.0) then 
                  ICHIPERPI(IR) = icfqi(ir)
     >                 / (qitoti-qiouti+qiconvi-qioutconvi)
               else 
                  ICHIPERPI(IR) = 0.0
               endif
c     
               if ((qetoti-qeouti+qeconvi-qeoutconvi).ne.0.0) then 
                  ICHIPERPE(IR) = icfqe(ir)
     >                 / (qetoti-qeouti+qeconvi-qeoutconvi)
               else
                  ICHIPERPE(IR) = 0.0
               endif

c     
               if ((qitoti+qetoti
     >              -(qiouti)-(qeouti)
     >              +(qiconvi)-(qioutconvi)
     >              +(qeconvi)-(qeoutconvi)).ne.0.0) then 
                  ixperpt(ir) =  (icfqi(ir)+icfqe(ir))
     >                 / (qitoti+qetoti
     >                 -(qiouti)-(qeouti)
     >                 +(qiconvi)-(qioutconvi)
     >                 +(qeconvi)-(qeoutconvi))
               else
                  ixperpt(ir) =  0.0
               endif

            elseif (usexpav) then

               if ((qitoto+qitoti).ne.0.0) then 
                  CHIPERPI(IR) = (cfqi(ir)-(qiconvi+qiconvo)
     >                 +(qioutconvi+qioutconvo)
     >                 +xpavi*(qiouti+qiouto))
     >                 / (qitoto+qitoti)
               else
                  CHIPERPI(IR) = 0.0
               endif
               

               if ((qetoto+qetoti).ne.0.0) then
                  CHIPERPE(IR) = (cfqe(ir)-(qeconvi+qeconvo)
     >                 +(qeoutconvi+qeoutconvo)
     >                 +xpave*(qeouti+qeouto))
     >                 / (qetoto+qetoti)
               else
                  CHIPERPE(IR) = 0.0
               endif

               if ((qitoto+qitoti+qetoto+qetoti).ne.0.0) then 
                  xperpt(ir) = ( (cfqi(ir)+cfqe(ir))
     >                 -(qiconvi+qiconvo)+(qioutconvi+qioutconvo)
     >                 +xpavt*(qiouti+qiouto)
     >                 -(qeconvi+qeconvo)+(qeoutconvi+qeoutconvo)
     >                 +xpavt*(qeouti+qeouto))
     >                 / (qitoto+qitoti+qetoto+qetoti)
               else
                  xperpt(ir) = 0.0
               endif

c     
c     Outer
c     
               if ((qitoto).ne.0.0) then
                  OCHIPERPI(IR) = (ocfqi(ir)
     >                 -qiconvo+qioutconvo+xpavio*qiouto)
     >                 / (qitoto)
               else
                  OCHIPERPI(IR) = 0.0
               endif
               
               if ((qetoto).ne.0.0) then 
                  OCHIPERPE(IR) = (ocfqe(ir)
     >                 -qeconvo+qeoutconvo+xpaveo*qeouto)
     >                 / (qetoto)
               else
                  OCHIPERPE(IR) = 0.0
               endif

               
               if ((qitoto+qetoto).ne.0.0) then 
                  oxperpt(ir) = ( (ocfqi(ir)+ocfqe(ir))
     >                 -(qiconvo)+(qioutconvo)
     >                 +xpavto*(qiouto)
     >                 -(qeconvo)+(qeoutconvo)
     >                 +xpavto*(qeouto))
     >                 / (qitoto+qetoto)
               else 
                  oxperpt(ir) = 0.0
               endif


c     
c     Inner
c     
               if ((qitoti)
     >              .ne.0.0) then                
                  ICHIPERPI(IR) = (icfqi(ir)
     >                 -qiconvi+qioutconvi+xpavii*qiouti)
     >                 / (qitoti)
               else
                  ICHIPERPI(IR) = 0.0
               endif

               if ((qetoti)
     >              .ne.0.0) then                
                  ICHIPERPE(IR) = (icfqe(ir)
     >                 -qeconvi+qeoutconvi+xpavei*qeouti)
     >                 / (qetoti)
               else
                  ICHIPERPE(IR) = 0.0
               endif


               if ((qitoti+qetoti)
     >              .ne.0.0) then                
                  ixperpt(ir) = ( (icfqi(ir)+icfqe(ir))
     >                 -(qiconvi)+(qioutconvi)
     >                 +xpavti*(qiouti)
     >                 -(qeconvi)+(qeoutconvi)
     >                 +xpavti*(qeouti))
     >                 / (qitoti+qetoti)
               else
                  ixperpt(ir) = 0.0
               endif


            else

               if ((qitoto+qitoti-(qiouti+qiouto))
     >              .ne.0.0) then                
                  CHIPERPI(IR) = (cfqi(ir)-(qiconvi+qiconvo)
     >                 +(qioutconvi+qioutconvo))
     >                 / (qitoto+qitoti-(qiouti+qiouto))
               else
                  CHIPERPI(IR) = 0.0
               endif


               if ((qetoto+qetoti-(qeouti+qeouto))
     >              .ne.0.0) then                
                  CHIPERPE(IR) = (cfqe(ir)-(qeconvi+qeconvo)
     >                 +(qeoutconvi+qeoutconvo))
     >                 / (qetoto+qetoti-(qeouti+qeouto))
               else
                  CHIPERPE(IR) = 0.0
               endif


               if ((qitoto+qitoti+qetoto+qetoti
     >              -(qiouti+qiouto)-(qeouti+qeouto))
     >              .ne.0.0) then                
                  xperpt(ir) = ( (cfqi(ir)+cfqe(ir))
     >                 -(qiconvi+qiconvo)+(qioutconvi+qioutconvo)
     >                 -(qeconvi+qeconvo)+(qeoutconvi+qeoutconvo))
     >                 / (qitoto+qitoti+qetoto+qetoti
     >                 -(qiouti+qiouto)-(qeouti+qeouto))
               else
                  xperpt(ir) = 0.0
               endif
c     
c     Outer
c     
               if ((qitoto-qiouto)
     >              .ne.0.0) then                
                  OCHIPERPI(IR) = (ocfqi(ir)-qiconvo+qioutconvo)
     >                 / (qitoto-qiouto)
               else
                  OCHIPERPI(IR) = 0.0
               endif


               if ((qetoto-qeouto)
     >              .ne.0.0) then                
                  OCHIPERPE(IR) = (ocfqe(ir)-qeconvo+qeoutconvo)
     >                 / (qetoto-qeouto)
               else
                  OCHIPERPE(IR) = 0.0
               endif


               if ((qitoto+qetoto
     >              -(qiouto)-(qeouto))
     >              .ne.0.0) then                
                  oxperpt(ir) = ( (ocfqi(ir)+ocfqe(ir))
     >                 -(qiconvo)+(qioutconvo)
     >                 -(qeconvo)+(qeoutconvo))
     >                 / (qitoto+qetoto
     >                 -(qiouto)-(qeouto))
               else
                  oxperpt(ir) = 0.0
               endif
c     
c     Inner
c     
               if ((qitoti-qiouti)
     >              .ne.0.0) then                
                  ICHIPERPI(IR) = (icfqi(ir)-qiconvi+qioutconvi)
     >                 / (qitoti-qiouti)
               else
                  ICHIPERPI(IR) = 0.0
               endif


               if ((qetoti-qeouti)
     >              .ne.0.0) then                
                  ICHIPERPE(IR) = (icfqe(ir)-qeconvi+qeoutconvi)
     >                 / (qetoti-qeouti)
               else
                  ICHIPERPE(IR) = 0.0
               endif


               if ((qitoti+qetoti
     >              -(qiouti)-(qeouti))
     >              .ne.0.0) then                
                  ixperpt(ir) = ( (icfqi(ir)+icfqe(ir))
     >                 -(qiconvi)+(qioutconvi)
     >                 -(qeconvi)+(qeoutconvi))
     >                 / (qitoti+qetoti
     >                 -(qiouti)-(qeouti))
               else
                  ixperpt(ir) = 0.0
               endif

            endif

c     
c     Calculate Averages if option is turned on
c     
            if (dpavopt.gt.0.and.ir.eq.(irsep+irskip+3)) then
               irtmp = irsep + irskip
c     
               xpave = 0.25*(chiperpe(irtmp) + chiperpe(irtmp+1)
     >              +chiperpe(irtmp+2) + chiperpe(irtmp+3))
               xpaveo= 0.25*(ochiperpe(irtmp) + ochiperpe(irtmp+1)
     >              +ochiperpe(irtmp+2) + ochiperpe(irtmp+3))
               xpavei= 0.25*(ichiperpe(irtmp) + ichiperpe(irtmp+1)
     >              +ichiperpe(irtmp+2) + ichiperpe(irtmp+3))
c     
               xpavi = 0.25*(chiperpi(irtmp) + chiperpi(irtmp+1)
     >              +chiperpi(irtmp+2) + chiperpi(irtmp+3))
               xpavio= 0.25*(ochiperpi(irtmp) + ochiperpi(irtmp+1)
     >              +ochiperpi(irtmp+2) + ochiperpi(irtmp+3))
               xpavii= 0.25*(ichiperpi(irtmp) + ichiperpi(irtmp+1)
     >              +ichiperpi(irtmp+2) + ichiperpi(irtmp+3))
c     
               xpavt = 0.25*(xperpt(irtmp) + xperpt(irtmp+1)
     >              +xperpt(irtmp+2) + xperpt(irtmp+3))
               xpavto = 0.25*(oxperpt(irtmp) + oxperpt(irtmp+1)
     >              +oxperpt(irtmp+2) + oxperpt(irtmp+3))
               xpavti = 0.25*(ixperpt(irtmp) + ixperpt(irtmp+1)
     >              +ixperpt(irtmp+2) + ixperpt(irtmp+3))
c     
               usexpav = .true.

               if (cprint.eq.1.or.cprint.eq.9) then
                  write (6,*) 'XPERP AV: rings ',irtmp,
     >                       ' to ',irtmp+3,' = ',
     >                        xpave,xpavi
               endif

            endif
c     
c     Record values for checking and printing
c     
c     
c     XperpE
c     
            flxval(ir,4,1) = qetoto
            flxval(ir,4,2) = qetoti
            flxval(ir,4,3) = qetoti + qetoto
c     
            flxval(ir,5,1) = qeconvo
            flxval(ir,5,2) = qeconvi
            flxval(ir,5,3) = qeconvo + qeconvi
c     
            flxval(ir,6,1) = qeouto
            flxval(ir,6,2) = qeouti
            flxval(ir,6,3) = qeouto + qeouti
c     
            flxval(ir,7,1) = qeoutconvo
            flxval(ir,7,2) = qeoutconvi
            flxval(ir,7,3) = qeoutconvo + qeoutconvi
c     
            flxval(ir,8,1) = ocfqe(ir)
            flxval(ir,8,2) = icfqe(ir)
            flxval(ir,8,3) = ocfqe(ir) + icfqe(ir)
c     
c     XperpI
c     
            flxval(ir,9,1) = qitoto
            flxval(ir,9,2) = qitoti
            flxval(ir,9,3) = qitoti + qitoto
c     
            flxval(ir,10,1) = qiconvo
            flxval(ir,10,2) = qiconvi
            flxval(ir,10,3) = qiconvo + qiconvi
c     
            flxval(ir,11,1) = qiouto
            flxval(ir,11,2) = qiouti
            flxval(ir,11,3) = qiouto + qiouti
c     
            flxval(ir,12,1) = qioutconvo
            flxval(ir,12,2) = qioutconvi
            flxval(ir,12,3) = qioutconvo + qioutconvi
c     
            flxval(ir,13,1) = ocfqi(ir)
            flxval(ir,13,2) = icfqi(ir)
            flxval(ir,13,3) = ocfqi(ir) + icfqi(ir)
c     
         end do


      endif

c     
c     Try calculating second derivative of N, Ti, Te and see what
c     these yield when applied to extracting teh value of Dperp.
c     
c     

c     
c     First calculate the second derivatives and their integrals over S
c     
      do ir = irsep,nrs
c     
c     Perform summation either over the entire ring or only
c     around the confined plasma
c     
         if (dpsuml.eq.0) then
            ikstart = ouend
            ikend = inend
         elseif (dpsuml.eq.1) then
            ikstart = 1
            ikend = nks(ir)
         endif
c     
         do ik = ikstart,ikend
            grad2n(ik,ir) = quant2grad(ir,kss(ik,ir),knbs,knds)
            grad2ti(ik,ir)= quant2grad(ir,kss(ik,ir),ktibs,ktids)
            grad2te(ik,ir)= quant2grad(ir,kss(ik,ir),ktebs,kteds)
c     
c     Include major radius correction if requested
c     
            if (dprcopt.eq.0) then
               fact = 1.0
            else
               fact = rs(ik,ir)/r0
            endif
c     

            totgrad2n(ir) =totgrad2n(ir)
     >           +grad2n(ik,ir)*surfa(ik,ir)*fact
            totgrad2ti(ir)=totgrad2ti(ir)
     >           +grad2ti(ik,ir)*surfa(ik,ir)*fact
            totgrad2te(ir)=totgrad2te(ir)
     >           +grad2te(ik,ir)*surfa(ik,ir)*fact
c     
         end do

      end do
c     
c     Calculate alternate values of Dperp/Xperp based on the second
c     derivative in density.
c     

c     
      if (cprint.eq.1.or.cprint.eq.9) then
c     
c     
c     Print out values for Dperp, Chiperpe,i
c     Both total for each ring as well as inner/outer values.
c     
         call prb
         call prchtml('  TRANSPORT COEFFICIENTS EXTRACTOR',
     >        'pr_transport_coeff','0','B')
         call prc('  Input Options:')
         call prb
c     
         call prc ('    METHOD OPTION:')
         if (dpmethod.eq.0) then
            call prc('      0 - Using Inner and Outer'//
     >               ' mid-plane gradients')
         elseif (dpmethod.eq.1) then
            call prc('      1 - Using Gradients for each cell'
     >           //' independently')
         elseif (dpmethod.eq.2) then
            call prc('      2 - Using Gradients for each cell'
     >           //' independently')
         endif
c     
         call prc('    SUMMATION OPTION:')
         if (dpsuml.eq.0) then
            call prc('      0 - For Transport above the X-point')
         elseif (dpsuml.eq.1) then
            call prc('      1 - For the entire field line')
         endif
c     
         call prc('    OUTER RING LOSS OPTION:')
         if (dpouter.eq.0) then
            call prc('      0 - Ignoring Perpendicular Fluxes across'//
     >           ' the last ring')
         elseif (dpsuml.eq.1) then
            call prc('      1 - Including Perpendicular Fluxes across'//
     >           ' the last ring')
         endif
c     
         call prc('    PERPENDICULAR CONVECTION OPTION:')
         if (dpconv.eq.0) then
            call prc('      0 - Ignoring Perpendicular Convection')
         elseif (dpconv.eq.1) then
            call prc('      1 - Including Perpendicular Convection')
         endif
c     
         call prc('    CELL CROSS-FIELD AREA OPTION:')
         if (dparea.eq.0) then
            call prc('      0 - Cross-field areas calculated at'//
     >           ' cell centres')
         elseif (dparea.eq.1) then
            call prc('      1 - Cross-field areas calculated at'//
     >           ' cell boundaries')
         endif
c     
         call prc('    CELL CENTRE CORRECTION OPTION: (Use Option 0)')
         if (dpfluxopt.eq.0) then
            call prc('      0 - Net Fluxes are calculated'//
     >               ' for whole ring')
         elseif (dpfluxopt.eq.1) then
            call prc('      1 - Net Fluxes are calculated'//
     >               ' for ring centre')
         endif
c     
         call prc('    USE OF AVERAGE DPERP OPTION:')
         if (dpavopt.eq.0) then
            call prc('      0 - Dperp/Xperp for Outer'//
     >               ' ring losses floats')
         elseif (dpavopt.gt.0) then
            call pri('      N - Dperp/Xperp Average'//
     >               ' Option set to',irskip)

            call pri2('          Averages over rings',
     >           irsep+irskip,irsep+irskip+3)

            call prr ('          Dperp  Average = ',dpav)
            call prr ('          XperpE Average = ',xpave)
            call prr ('          XperpI Average = ',xpavi)
         endif
c     
         call prc('    MAJOR RADIUS CORRECTION OPTION:')
         if (dprcopt.eq.0) then
            call prc('      0 - Not corrected for Major-Radius effects')
         elseif (dprcopt.eq.1) then
            call prc('      1 - Major Radius Correction'//
     >               ' Rcell/R0 applied')
            call prc('       to Areas and to Ionization Source.')
         elseif (dprcopt.eq.2) then
            call prc('      2 - Major Radius Correction'//
     >               ' Rcell/R0 applied')
            call prc('       to Areas ONLY - not to ionization source.')
         endif
c     
         call prc('    GRADIENT CALCULATION METHOD:')
         if (dpnav.eq.-1) then
            call prc('     -1 - Gradient calculated using derivatives'//
     >           ' taken from ')
            call prc('          cubic spline iterpolation procedure')
         elseif (dpnav.eq.0) then
            call prc('      0 - Gradient calculated taking average of'//
     >           ' slopes to ')
            call prc('          neighbouring cross-field cells')
         elseif (dpnav.eq.1) then
            call prc('      1 - Gradient calculated taking '//
     >           'slope between ')
            call prc('          neighbour cross-field cells.')
            call prc('          Value at current cell is not used.')
         elseif (dpnav.eq.2) then
            call prc('      2 - Gradient calculated taking '//
     >           'slope between current and next outward')
            call prc('          cross-field cell only.')
         endif
c     
         call prc('    GRADIENT SMOOTHING OPTION:')
         if (dpsmooth.eq.0) then
            call prc('      0 - Smoothing is turned OFF')
         elseif (dpsmooth.gt.0) then
            call pri('      >0- Smoothing is turned ON.  Opt =',
     >           dpsmooth)
            call pri('          Smoothing performed over N knots. N=',
     >           2*dpsmooth+1)
            call prc('          Spaced evenly radially around'
     >           //' the cell of interest')

         endif
c     
         call prc('    POWER LOSS OPTION:')
         if (dpploss.eq.0) then
            call prc('      0 - Volume Power Loss terms are OFF')
         elseif (dpploss.eq.1) then
            call prc('      1 - Volume Power Losses from PIN for'
     >           //' ions and electrons.')
         elseif (dpploss.eq.2) then
            call prc('      2 - Volume Power Losses from PIN for'
     >           //' ions and electrons.')
            call prc('          Electron-Ion heat transfer (Pei) also'
     >           //' included.')
            call prr('          Pei correction factor for OSKIN = ',
     >               dppei)
         endif
c     
         call prc('    NON-ORTHOGONAL CORRECTION OPTION:')
         if (dporth.eq.0) then
            call prc('      0 - Non-othogonal grid correction'
     >           // ' to gradients is  OFF')
         elseif (dporth.eq.1) then
            call prc('      1 - Non-othogonal grid correction'
     >           //' to gradients is  ON')
         endif
c     
         call prr('    RECYLCE FRACTION APPLIED TO EXTRACTOR: ',dprec)
c     
         if (dpxpratio.eq.-1.0) then
            call prc('    DPERP/XPERP FIXED RATIO OPTION IS OFF')
            call prr('    DPERP VALUE FOR XPERP EXTRACTION FIXED AT :',
     >           cdperp)
         elseif (dpxpratio.le.0.0) then
            call prc('    DPERP/XPERP FIXED RATIO OPTION IS OFF')
            if (dpavopt.eq.0) then
               call prc('    DPERP VALUE FOR XPERP EXTRACTION'//
     >              ' IS VALUE CALCULATED FOR EACH RING')
            elseif (dpavopt.gt.0) then
               call prr('    DPERP VALUE FOR XPERP EXTRACTION'//
     >              ' IS CALCULATED AVERAGE VALUE: ',dpav)
            endif
         elseif (dpxpratio.gt.0.0) then
            call prr('    DPERP/XPERP FIXED RATIO HAS BEEN SET TO: ',
     >           dpxpratio )
         endif
c     
c     
         call prb
         call prr('    Value of GammaI used for Xperp calculation:'
     >        , gai)
         call prr('    Value of GammaE used for Xperp calculation:'
     >        ,gae)
         call prb
         call prc('  Table of Dperp values extracted from OSM')
         call prc('  Ring         '//Inner//'          '//Outer//
     >        '          Total')
         do ir = irsep , irwall-2
            write(coment,300) ir,idperp(ir),odperp(ir),dperp(ir)
            call prc(coment)
         end do
c     
         call prc(' Table of Xperp ION values extracted from OSM')
         call prc('  Ring         '//Inner//'           '//Outer//
     >        '          Total')
         do ir = irsep , irwall-2
            write(coment,300) ir,ichiperpi(ir),ochiperpi(ir),
     >           chiperpi(ir)
            call prc(coment)
         end do
c     
         call prc(' Table of Xperp ELECTRON extracted from OSM')
         call prc('  Ring         '//Inner//'           '//Outer//
     >        '          Total')
         do ir = irsep , irwall-2
            write(coment,300) ir,ichiperpe(ir),ochiperpe(ir),
     >           chiperpe(ir)
            call prc(coment)
         end do
c     
         call prc(' Table of Xperp TOTAL extracted from OSM')
         call prc('  Ring         '//Inner//'           '//Outer//
     >        '          Total')
         do ir = irsep , irwall-2
            write(coment,300) ir,ixperpt(ir),oxperpt(ir),
     >           xperpt(ir)
            call prc(coment)
         end do
c     
c     Detailed Ring summaries
c     
         write (6,*) 'Dperp Components'
         do ir = irsep,irwall-1
            write(6,'(i4,9(e12.4,1x))') ir,
     >           ((flxval(ir,in1,in2),in1=1,3),in2=1,3)
         end do
c     
         write (6,*) 'Dperp Components - TOTAL'
         write (6,*) 'DPTOT(IR) = SIGMA(IK)'//
     >        ' (AREA(IK,IR) * GRADN(IK,IR))'
         write (6,*) 'DELTOT(IR)= SIGMA(IK)'//
     >        '(AREA(IK,IRWALL-1) * GRADN(IK,IRWALL-1))'
         write (6,*) 'CFGAM(IRSEP) = TOTFLX-NETFLX(IRSEP)'
         write (6,*) 'CFGAM(IR)    = CFGAM(IR-1) - NETFLX(IR)'
         write (6,*)
         write (6,*) ' IR     DPTOT        DELTOT'//
     >        '        CFGAM        DPERP        NETFLX'//
     >        '       IONIS         FLUX       TOTFLUX'
c     
         write (6,*) 'MAIN SOL:'
c     
         tmpnflx = 0.0
         tmpion  = 0.0
         tmpflx  = 0.0
c     
         do ir = irsep,irwall
            write(6,'(i4,8(e12.4,1x))') ir,
     >           flxval(ir,1,3),flxval(ir,2,3),flxval(ir,3,3),
     >           dperp(ir),netflx(ir),ionis(ir),flux(ir),totflx
            tmpflx = tmpflx + flux(ir)
            tmpion = tmpion + ionis(ir)
            tmpnflx = tmpnflx + netflx(ir)
         end do
         write(6,'(''TOT:'',4(13x),3(e12.4,1x))')
     >        tmpnflx,tmpion,tmpflx
c     
         write (6,*) 'PP:'

         tmpnflx = 0.0
         tmpion  = 0.0
         tmpflx  = 0.0
c     
         do ir = irtrap,nrs
            write(6,'(i4,4(13x),4(e12.4,1x))') ir,
     >           netflx(ir),ionis(ir),flux(ir)
            tmpflx = tmpflx + flux(ir)
            tmpion = tmpion + ionis(ir)
            tmpnflx = tmpnflx + netflx(ir)
         end do
c     
         write(6,'(''TOT:'',4(13x),3(e12.4,1x))')
     >        tmpnflx,tmpion,tmpflx
c     
         write (6,*) 'CORE:'

         tmpion  = 0.0
c     
         do ir = 1,irsep-1
            write(6,'(i4,5(13x),4(e12.4,1x))') ir,
     >           ionis(ir)
            tmpion = tmpion + ionis(ir)
         end do
c     
         write(6,'(''TOT:'',5(13x),3(e12.4,1x))')
     >        tmpion
         write(6,'(a,g12.6)') 'TOTAL: ALL IZ = ',alliz
c     
c     Dperp - Alternate calculation
c     
         write (6,*)
         write (6,*) 'Dperp Components - ALTERNATE'
         write (6,*)
         write (6,*) 'Dperp = [ Phi(ir) - INT(Siz(ir)) ] / '
         write (6,*) '        INT [Aperp * { (dn/dr)|in -'//
     >        ' (dn/dr)|out }]'
         write (6,*)
         write (6,*) ' IR     APERP   INT(DN/DR|IN-DN/DR|OUT)'//
     >        '        INT(AP DN/DR)        DPERP        NETFLX'//
     >        '       IONIS         FLUX   '
c     
         do ir = irsep+1,irwall-2
            aperp = 0.0
            dgradt = 0.0
            apdg = 0.0

            do ik = 1,nks(ir)

               aperp = aperp + surfa(ik,ir)
               dgradt = dgradt + dgradn(ik,ir)
               apdg = apdg + surfa(ik,ir) * dgradn(ik,ir)

            end do
c     
            dperpt = netflx(ir) / apdg
c     
            write(6,'(i4,8(e12.4,1x))') ir,
     >           aperp,dgradt,apdg,dperpt,
     >           netflx(ir),ionis(ir),flux(ir)
c     
         end do
c     
c     X perp components
c     
         write(6,*)
c     
         write (6,*) 'XperpE Components'
         do ir = irsep,irwall-1
            write(6,'(i4,15(e12.4,1x))') ir,
     >           ((flxval(ir,in1,in2),in1=4,8),in2=1,3)
         end do
c     
         write (6,*) 'Xperpi Components'
         do ir = irsep,irwall-1
            write(6,'(i4,15(e12.4,1x))') ir,
     >           ((flxval(ir,in1,in2),in1=9,13),in2=1,3)
         end do
c     
c     End of cprint IF
c     
      endif
c     
c     
c     Calculate the core outflux based on the density gradients in the
c     calculated background and the imposed Dperp given for the case.
c     
c     write (6,*) 'cioptg:',cioptg,cre2d
c     
c     if (cioptg.eq.99) then
c     
c     write (6,*) 'irsep:',irsep,nks(irsep)
c     
      totcfflux = 0.0
      itotcfflux = 0.0
      ototcfflux = 0.0
c     
      ir = irsep
c     
      do ik = 1,nks(irsep)
         ikin = ikins(ik,irsep)
         irin = irins(ik,irsep)
         if (irin.eq.irsep-1) then
            in = korpg(ik,irsep)
            area = sqrt( (rvertp(1,in)-rvertp(4,in))**2
     >           +(zvertp(1,in)-zvertp(4,in))**2)
            gradnsep = (knbs(ikin,irin)-knbs(ik,irsep))
     >           /kinds(ik,irsep)

            if (kss(ik,irsep).lt.(ksmaxs(irsep)/2.0)) then 
c     
               ototcfflux = ototcfflux +  cdperp * area * gradnsep
c     
            else
c     
               itotcfflux = itotcfflux +  cdperp * area * gradnsep
c     
            endif
c     
            totcfflux = totcfflux + cdperp * area * gradnsep
c     
c     IPP/08 Krieger - fixed debug output -> writing area,
c     gradnsep and totcfflux makes only sense if they have been
c     computed before

            if (cprint.eq.1.or.cprint.eq.9) then

               write (6,*) 'CF:',ik,ir,ikin,irin,in,area,gradnsep,
     >              totcfflux,cdperp
            endif
c     
         endif


      end do
c     
      if (cprint.eq.1.or.cprint.eq.9) then
         call prb
         call prc('  Calculation of Cross-field Flux from Core '//
     >        'for actual BG plasma')
         call prr('  Value of Dperp assumed:',cdperp)
         call prr('  Total CF Flux       = ',totcfflux)
         call prr2('  Total '//inner//'/'//outer//' CF Flux     = ',
     >        itotcfflux,ototcfflux)
         call prr('  Total Target Flux       = ',totflx_save)
         call prr2('  Total '//inner//'/'//outer//' Target Flux = ',
     >        itotflx_save,ototflx_save)
         call prr('  Total Ionization       = ',tionis)
         call prr2('  Total '//inner//'/'//outer//' Ionization  = ',
     >        itionis,otionis)
         call prr('  Total NET Flux       = ',netflx_save)
         call prr2('  Total '//inner//'/'//outer//' NET Flux = ',
     >        inetflx_save,onetflx_save)

         call prb
c     
c     call targflux(totflx,totrec)
c     
      endif
c     
c     endif
c     
      if (cre2d.eq.1.or.cre2d.eq.2) then
c     
         totcfflux = 0.0
         itotcfflux = 0.0
         ototcfflux = 0.0
         ir = irsep
c     
c     write (6,*) 'irsep:',irsep,nks(irsep)
c     
         do ik = 1,nks(irsep)
            ikin = ikins(ik,irsep)
            irin = irins(ik,irsep)
            if (irin.eq.irsep-1) then
               in = korpg(ik,irsep)
               area = sqrt((rvertp(1,in)-rvertp(4,in))**2
     >              +(zvertp(1,in)-zvertp(4,in))**2)
               gradnsep = (e2dnbs(ikin,irin)-e2dnbs(ik,irsep))
     >              /kinds(ik,irsep)
c     
               if (kss(ik,irsep).lt.(ksmaxs(irsep)/2.0)) then 
c     
                  ototcfflux = ototcfflux +  cdperp * area * gradnsep
c     
               else
c     
                  itotcfflux = itotcfflux +  cdperp * area * gradnsep
c     
               endif

               totcfflux = totcfflux + cdperp * area * gradnsep

            endif

            if (cprint.eq.1.or.cprint.eq.9) then
               write (6,*) 'CF:',ik,ir,ikin,irin,in,area,gradnsep,
     >              totcfflux,cdperp
            endif

         end do

         if (cprint.eq.1.or.cprint.eq.9) then
            call prb
            call prc('  Calculation of Cross-field Flux'//
     >           ' from Core for FLUID CODE Solution')
            call prr('  Value of Dperp assumed:',cdperp)
c     
            call prr('  Total CF Flux       = ',totcfflux)
            call prr2('  Total '//inner//'/'//outer//' CF Flux     = ',
     >           itotcfflux,ototcfflux)
c     
            call prb
c     
         endif
      endif


c     
c     
c     Format statements
c     

 300  format (2x,i4,3x,3(2x,g13.5))
 400  format (a13,1x,i3,1x,11g10.3)
 500  format ('Title:',7x,'Ring',5x, 'Dperp','TotOutFlux',
     >     5x,'Ngrad',
     >     6x,'Area',7x,'Sep',6x,'Flux',5x,'Ionis',
     >     3x,'Netflux',6x,'Nmid',
     >     3x,'Sep IN',2x,'Sep OUT')
 700  format ('Title:',7x,'Ring',4x, 'Xperpi',5x,'TotQi',
     >     4x,'Tigrad',
     >     6x,'Area',7x,'Sep',6x,'Heat',5x,'Timid',
     >     3x,'Sep IN',2x,'Sep OUT')
 900  format ('Title:',7x,'Ring',4x, 'Xperpe',5x,'TotQe',
     >     4x,'Tegrad',
     >     6x,'Area',7x,'Sep',6x,'Heat',5x,'Temid',
     >     4x,'Sep IN',3x,'Sep OUT')
 1200 format(i4,9(g11.3))
 1210 format('Ring',5x,'NfluxI',6x,'FluxI',5x,'Iionis',
     >     5x,'NfluxO',
     >     6x,'FluxO',5x,'Oionis',6x,'Nflux',7x,'Flux',6x,'ionis')
C     

      return
      END
c
c
c
      subroutine calc_divrec(totrec)
      IMPLICIT NONE
      real totrec
C
C*********************************************************************
C
C     CALC_DIVREC: Calculates the DIVIMP estimate of particle
c                  recombintaion source.
C
C*********************************************************************
C
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'cadas'
c
c      include 'cedge2d'
C
C
      INTEGER IK,IR,id
C
      real totrect
c
      write (6,*) 'Calculating DIVREC:'
C
C---- Calculate the number of recombinations
C---- Store the recombination array for later use
C
      TOTREC = 0.0
      totrect = 0.0
C
C---- ONE RING AT A TIME.
C
      DO IR = 1, NRS
        DO IK = 1, NKS(IR)
          if (ktebs(ik,ir).lt.treccut) then
             PTESA(IK) = treccut
          else
             PTESA(IK) = KTEBS(IK,IR)
          endif
          PNESA(IK) = KNBS(IK,IR) * RIZB
          PNBS(IK) = KNBS(IK,IR)
          PNHS(IK) = PINATOM(IK,IR)
        ENDDO
C
c       Recombination
c
        write(year,'(i2.2)') iyearh
        call xxuid(useridh)
        ICLASS = 1
c
        if (crecopt.eq.4) then
           CALL ADASRD(YEAR,1,1,ICLASS,NKS(IR),PTESA,PNESA,PCOEF)
        else
           CALL otherrec(NKS(IR),PTESA,PNESA,PCOEF,crecopt)
        endif
c
        DO IK = 1, NKS(IR)
          DIVREC(IK,IR) = (PNESA(IK)*PCOEF(IK,1))*PNBS(IK)
          TOTREC = TOTREC + DIVREC(IK,IR)*KAREA2(IK,IR)
          TOTRECt = TOTRECt + DIVREC(IK,IR)*KAREA2(IK,IR)*rs(ik,ir)/r0
        ENDDO
      ENDDO

c
c      CALL PRRMATDIV(DIVREC,MAXNKS,nks(irsep),NRS,6,'RECOMBINATION')
c
      return
      end
c
c
c
      subroutine nimwall
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
c
c     NIMWALL: This routine extracts the coordinates of
c              the NIMBUS MAIN WALL from the data read in
c              from PIN and assigns it to the array wallco.
c              This is then used in a call to DOWALL to
c              calculate the actual complete wall - linking
c              this specification with the DIVIMP targets.
c              THE NIMBUS and DIVIMP segmenst across the
c              target are the same IF target option 6 has
c              been selected - HOWEVER - the actual numbers
c              for the coordinates differ by a small
c              amounts due to rounding - so the actual
c              DIVIMP values are used.
c
c
c              NIMWALL2 performs the same function for points in
c              the private plasma wall region.
c
c
c     If the vessel has been redefined to include any baffles
c     then this code changes the contents of the NIMBUS vessel arrays
c     and adjusts the flux arrays to match.
c
c
      integer in,ik,ir,cnt
      logical add
c
c
c     Initialize
c
      add=.false.
      cnt = 0
c
c      in = 1
c
c      do while (in.le.nvesm)
c
c
      do in =1 ,nvesm
c
c         Need to load from the second end of the "Corner" point at
c         end to the first point of the corner segment at the other.
c
c         Take the FIRST corner segment to the LAST corner segment.
c
c         The "OUTER" corner segment is coded as 2
c              INNER            "                5
c
          if (cprint.eq.3.or.cprint.eq.9)
     >       write(6,*) 'Nim wall:',add,in,cnt,rvesm(in,1),
     >                   zvesm(in,1),rvesm(in,2),zvesm(in,2),
     >                   jvesm(in)
c
          if (add) then
             cnt = cnt + 1
             wallco(cnt,1) = rvesm(in,1)
             wallco(cnt,2) = zvesm(in,1)
c
          endif
c
c          if ((.not.add).and.jvesm(in).eq.2.or.jvesm(in+1).eq.3) then
c
          if ((.not.add).and.jvesm(in).eq.2.or.jvesm(in).eq.3.or.
     >                       jvesm(in).eq.5.or.jvesm(in).eq.6.or.
     >                       jvesm(in).eq.7) then
c
             add = .true.
c
          endif
c
c          if (add.and.jvesm(in).eq.5.or.jvesm(in-1).eq.4) then
c
          if (add.and.jvesm(in+1).eq.4) then
c
             add=.false.
             goto 100

c
c             in = nvesm
c
          endif
c
c          in = in + 1
c
      end do
c
 100  nwall = cnt
c
c
      write (6,*) 'NIMBUS wall:',nvesm,nwall
      do in = 1,nvesm
         write(6,'(i5,4(1x,f13.5),i5)') in,rvesm(in,1),
     >           zvesm(in,1),rvesm(in,2),zvesm(in,2),jvesm(in)
      end do


      return
c
      end
c
c
c
      subroutine nimwall2
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
c
c     NIMWALL2:This routine extracts the coordinates of
c              the NIMBUS TRAP WALL from the data read in
c              from PIN and assigns it to the array wallco2.
c              This is then used in a call to DOWALL to
c              calculate the actual completed wall - linking
c              this specification with the DIVIMP targets.
c              THE NIMBUS and DIVIMP segments across the
c              target are the same IF target option 6 has
c              been selected - HOWEVER - the actual numbers
c              for the coordinates differ by a small
c              amounts due to rounding - so the actual
c              DIVIMP values are used.
c
c
c              NIMWALL performs the same function for points in
c              the private plasma wall region.
c
      integer in,ik,ir,cnt
      logical add
c
c     Initialize
c
      add=.false.
      cnt = 0
c
c      in = 1
c
c      do while (in.le.nvesm)
c
      do in = 1,nvesm
c
c         Need to load from the second end of the "Corner" point at
c         end to the first point of the corner segment at the other.
c
c         Take the FIRST corner segment to the LAST corner segment.
c
c         The "OUTER" corner segment is coded as 2
c              INNER            "                5
c         The  TRAP   segments are      coded as 8
c
          if (cprint.eq.3.or.cprint.eq.9)
     >        write(6,*) 'Nim wall2:',add,in,cnt,rvesm(in,1),
     >                        zvesm(in,1),jvesm(in)
c
          if (add) then
             cnt = cnt + 1
             wallco2(cnt,1) = rvesm(in,1)
             wallco2(cnt,2) = zvesm(in,1)
          endif
c
          if ((.not.add).and.(jvesm(in).eq.8.or.
     >             jvesm(in).eq.9.or.jvesm(in).eq.10)) then
             add = .true.
          endif
c
c         The PP wall is actually usually the last segment
c         in the specification - the following test is
c         only useful in case that changes. Otherwise - the
c         values all the way to NVESM will be loaded into
c         wallco2.
c
          if (add.and.(jvesm(in+1).ne.8.and.
     >                 jvesm(in+1).ne.9.and.
     >                 jvesm(in+1).ne.10)) then
c
             add=.false.
c
c             in = nvesm
c
              goto 100
          endif
c
c         in = in + 1
c
      end do
c
 100  nwall2 = cnt
c
c
      if (cprint.eq.3.or.cprint.eq.9) then
c
         write (6,*) 'NIMBUS wall:',nvesm,nwall2
         do in = 1,nvesm
            write(6,'(i5,2(1x,f13.5),i5)') in,rvesm(in,1),
     >           zvesm(in,1),jvesm(in)
         end do
c
      endif
c
      return
c
      end
c
c
c
      subroutine nimind
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
c slmod begin
      include 'slcom'
c slmod end
c
c     NIMIND: This routine calculates and assigns indices (pointers)
c             into the arrays of NIMBUS/PIN flux data from the
c             re-mapped PIN/NIMBUS wall segments.
c
c             For the entire wall (including targets)- it stores the
c             numbers in the wallpt(ind,17) array element.
c
c             It also does an independent calculation for the targets
c             based on the NDS target index. The elements for
c             1,ndsin, ndsin+1 and nds - are set to null index
c             pointers because they do not refer to real target
c             segments for target option 6. The length of these
c             target segments dds(id) has been set to zero. This
c             index is stored in the nimindex(id) and wallindex arrays.
c
      integer in,id
      integer wlo,wli,wt,wl,ndso,ndsi
c
c     Variables to identify wall inconsistencies
c
      integer nmain,nit,not,npp
      integer ndmain,ndit,ndot,ndpp
      integer otarg,itarg, extra_index_in,extra_index_out
      real tol
      parameter(tol=0.00005)
c slmod begin
      IF (grdnmod.NE.0) THEN
c...    Skip this, as it is done in AssignNimbusWall (I think):
        WRITE(0,*) 'SKIPPING CALL TO NIMIND - TROUBLE?'
        RETURN
      ENDIF
c slmod end
c
c     Initialize
c
      call izero(nimindex,maxnds)
      call izero(wallindex,maxnds)
c
c     Cycle through the NIMBUS wall array - assigning the
c     appropriate indices - also check to make sure that they are
c     the correct corresponding indices.
c
      in = 1
      extra_index_in  = -1
      extra_index_out = -1
c
c     Targets
c
      ndso = ndsin+2
      ndsi = 2
c
c     Walls
c
      wlo = wltrap2+1
      wli = wlwall2+1
c
      wt = wltrap1
      wl = wlwall1
c
      write (6,*) 'Wallpt:',wlwall1,wlwall2,wltrap1,wltrap2,wallpts
c
      write (6,*) 'Num:',wallpts,nvesm
c
      if (wallpts.ne.nvesm) then

         call prc('WARNING ERROR: WALL DEFINITION'//
     >          ' FROM NIMBUS DOES NOT MATCH DIVIMP')
         call pri2('             : TOTAL WALL ELEMENTS (D/N) ',
     >             wallpts,nvesm)
c         call prc ('PROGRAM STOPPING')

         write(0,*) 'WARNING ERROR: WALL DEFINITION'//
     >              ' FROM NIMBUS DOES NOT MATCH DIVIMP'
         write(0,*) '             : TOTAL WALL ELEMENTS (D/N) ',
     >               wallpts,nvesm
c         write(0,*) 'PROGRAM STOPPING'

         write(6,*) 'WARNING ERROR: WALL DEFINITION'//
     >              ' FROM NIMBUS DOES NOT MATCH DIVIMP'
         write(6,*) '             : TOTAL WALL ELEMENTS (D/N) ',
     >               wallpts,nvesm
c
c        Only STOP if the case is running impurities - for a
c        background plasma solution allow the case to continue.
c
         if ((ctestsol.ne.-1).and.((nvesm-wallpts).gt.2)) then
c
            write(6,*) 'PROGRAM STOPPING'
            write(0,*) 'PROGRAM STOPPING'
            write(7,*) 'PROGRAM STOPPING'
            stop
c
         endif
c
c        Find the offending element/s - PIN/NIMBUS has now started
c        splitting some target elements into two pieces - need to
c        determine which piece and combine the PIN/NIMBUS results
c        for this segment so that it correctly maps to the one
c        segment used in DIVIMP.
c
c
c        Where is the extra segment ?
c
         nmain = 0
         not   = 0
         nit   = 0
         npp   = 0

c
         do in = 1,nvesm
c
c           Main wall
c

            if (jvesm(in).eq.2.or.
     >          jvesm(in).eq.3.or.
     >          jvesm(in).eq.5.or.
     >          jvesm(in).eq.6.or.
     >          jvesm(in).eq.7) then


                nmain = nmain + 1

c
c           Outer target
c
            elseif(jvesm(in).eq.1) then

                not = not +1

c
c           Inner target
c
            elseif(jvesm(in).eq.4) then

                nit = nit+1


c
c           Private plasma wall
c
            elseif(jvesm(in).eq.8.or.
     >             jvesm(in).eq.9.or.
     >             jvesm(in).eq.10) then

                npp = npp +1


            endif

         end do

         ndmain = wlwall2-wlwall1+1
         ndpp   = wltrap2-wltrap1+1
         ndit   = (wltrap1-1)-(wlwall2+1) +1
         ndot   = wallpts-(wltrap2+1) + 1


         write(6,*) 'Main:', nmain, ndmain
         write(6,*) 'IT  :', nit, ndit
         write(6,*) 'OT  :', not, ndot
         write(6,*) 'PP  :', npp, ndpp

c
c        Loop through NIMBUS wall and identify the element the
c        has been split.
c
c
         otarg  = wlo
         itarg  = wli
c
         do in = 1,nvesm

            if (jvesm(in).eq.4.and.(nit.ne.ndit).and.
     >          extra_index_in.eq.-1) then
c
               if (abs(wallpt(itarg,20)-rvesm(in,1)).gt.tol.or.
     >             abs(wallpt(itarg,21)-zvesm(in,1)).gt.tol.or.
     >             abs(wallpt(itarg,22)-rvesm(in,2)).gt.tol.or.
     >             abs(wallpt(itarg,23)-zvesm(in,2)).gt.tol) then

                   write(6,*) 'Found INNER Point:',in,itarg
                   write(6,'(8g13.6)')
     >               wallpt(itarg,20),rvesm(in,1),
     >               wallpt(itarg,21),zvesm(in,1),
     >               wallpt(itarg,22),rvesm(in,2),
     >               wallpt(itarg,23),zvesm(in,2)
                   write(6,'(8g13.6)')
     >               wallpt(itarg+1,20),rvesm(in+1,1),
     >               wallpt(itarg+1,21),zvesm(in+1,1),
     >               wallpt(itarg+1,22),rvesm(in+1,2),
     >               wallpt(itarg+1,23),zvesm(in+1,2)


                   extra_index_in = in

               endif

               itarg = itarg + 1

            endif
c
c           Outer
c
            if (jvesm(in).eq.1.and.(not.ne.ndot).and.
     >          extra_index_out.eq.-1) then
c
               if (abs(wallpt(otarg,20)-rvesm(in,1)).gt.tol.or.
     >             abs(wallpt(otarg,21)-zvesm(in,1)).gt.tol.or.
     >             abs(wallpt(otarg,22)-rvesm(in,2)).gt.tol.or.
     >             abs(wallpt(otarg,23)-zvesm(in,2)).gt.tol) then

                   write(6,*) 'Found OUTER Point:',in,otarg
                   write(6,'(8g13.6)')
     >               wallpt(otarg,20),rvesm(in,1),
     >               wallpt(otarg,21),zvesm(in,1),
     >               wallpt(otarg,22),rvesm(in,2),
     >               wallpt(otarg,23),zvesm(in,2)
                   write(6,'(8g13.6)')
     >               wallpt(otarg+1,20),rvesm(in+1,1),
     >               wallpt(otarg+1,21),zvesm(in+1,1),
     >               wallpt(otarg+1,22),rvesm(in+1,2),
     >               wallpt(otarg+1,23),zvesm(in+1,2)

                   extra_index_out = in

               endif

               otarg = otarg + 1


            endif


         end do
c
         write(6,*) 'Some extra_index quantities:'
c
         if (extra_index_in.ne.-1) then
            write (6,*) 'INNER:'
            write(6,'(10g10.3)')
     >         fluxhw(extra_index_in),fluxhw(extra_index_in+1),
     >         flxhw2(extra_index_in),flxhw2(extra_index_in+1),
     >         flxhw3(extra_index_in),flxhw3(extra_index_in+1),
     >         flxhw4(extra_index_in),flxhw4(extra_index_in+1),
     >         flxhw5(extra_index_in),flxhw5(extra_index_in+1)
         endif
c
         if (extra_index_in.ne.-1) then
            write (6,*) 'OUTER:'
            write(6,'(10g10.3)')
     >         fluxhw(extra_index_out),fluxhw(extra_index_out+1),
     >         flxhw2(extra_index_out),flxhw2(extra_index_out+1),
     >         flxhw3(extra_index_out),flxhw3(extra_index_out+1),
     >         flxhw4(extra_index_out),flxhw4(extra_index_out+1),
     >         flxhw5(extra_index_out),flxhw5(extra_index_out+1)
         endif
c
c     End of wall fix up routine
c
      endif
c
      do in = 1,nvesm
c
         if (in.eq.(extra_index_in+1)) cycle
         if (in.eq.(extra_index_out+1)) cycle
c
c        Main Wall including divertor and corners
c
         if (jvesm(in).eq.2.or.
c
c     >      (jvesm(in).eq.1.and.jvesm(in+1).ne.1
c     >                     .and.jvesm(in+1).ne.2).or.
c     >      (jvesm(in).eq.4.and.jvesm(in-1).ne.4
c     >                     .and.jvesm(in-1).ne.5).or.
c
     >      jvesm(in).eq.3.or.
     >      jvesm(in).eq.5.or.jvesm(in).eq.6.or.
     >      jvesm(in).eq.7) then
c
            wallpt(wl,17) = in

c
            if (cprint.eq.3.or.cprint.eq.9) then
c
               write (6,'(''N:'',4(1x,f18.10),2i4)')
     >            rvesm(in,1),zvesm(in,1),
     >            rvesm(in,2),zvesm(in,2),in,jvesm(in)
               write (6,'(''W:'',4(1x,f18.10))')
     >      wallpt(wl,1) + wallpt(wl,5) * cos(wallpt(wl,8)),
     >      wallpt(wl,2) + wallpt(wl,5) * sin(wallpt(wl,8)),
     >      wallpt(wl,1) + wallpt(wl,6) * cos(wallpt(wl,9)),
     >      wallpt(wl,2) + wallpt(wl,6) * sin(wallpt(wl,9))
               write (6,'(''T:'',4(1x,f18.10),2i4)')
     >            wallpt(wlo,1),
     >            wallpt(wlo,2),
     >            rp(ndso),zp(ndso),wl
               write (6,*)
c
            endif
c
            wl = wl + 1
c
c        Outer target
c
         elseif (jvesm(in).eq.1) then
c
            wallindex(ndso) = wlo
c
            wallpt(wlo,17) = in
            wallpt(wlo,18) = ndso
            nimindex(ndso) = in
c
c
c
            if (cprint.eq.3.or.cprint.eq.9) then
c
               write (6,'(''N:'',4(1x,f18.10),2i4)')
     >            rvesm(in,1),zvesm(in,1),
     >            rvesm(in,2),zvesm(in,2),in,jvesm(in)
               write (6,'(''W:'',4(1x,f18.10))')
     >      wallpt(wlo,1) + wallpt(wlo,5) * cos(wallpt(wlo,8)),
     >      wallpt(wlo,2) + wallpt(wlo,5) * sin(wallpt(wlo,8)),
     >      wallpt(wlo,1) + wallpt(wlo,6) * cos(wallpt(wlo,9)),
     >      wallpt(wlo,2) + wallpt(wlo,6) * sin(wallpt(wlo,9))
               write (6,'(''O:'',4(1x,f18.10),2i4)')
     >            wallpt(wlo,1),
     >            wallpt(wlo,2),
     >            rp(ndso),zp(ndso),wlo,ndso
               write (6,*)
c
            endif
c
            wlo = wlo + 1
            ndso = ndso + 1
c
c        Inner Target
c
         elseif (jvesm(in).eq.4) then
c
            wallindex(ndsi) = wli
c
            wallpt(wli,17) = in
            wallpt(wli,18) = ndsi
            nimindex(ndsi) = in
c
            if (cprint.eq.3.or.cprint.eq.9) then
c
               write (6,'(''N:'',4(1x,f18.10),2i4)')
     >            rvesm(in,1),zvesm(in,1),
     >            rvesm(in,2),zvesm(in,2),in,jvesm(in)
               write (6,'(''W:'',4(1x,f18.10))')
     >      wallpt(wli,1) + wallpt(wli,5) * cos(wallpt(wli,8)),
     >      wallpt(wli,2) + wallpt(wli,5) * sin(wallpt(wli,8)),
     >      wallpt(wli,1) + wallpt(wli,6) * cos(wallpt(wli,9)),
     >      wallpt(wli,2) + wallpt(wli,6) * sin(wallpt(wli,9))
               write (6,'(''I:'',4(1x,f18.10),2i4)')
     >            wallpt(wli,1),
     >            wallpt(wli,2),
     >            rp(ndsi),zp(ndsi),wli,ndsi

               write (6,*)
c
            endif
c
            wli = wli + 1
            ndsi = ndsi + 1
c
c        Private Plasma indices
c
         elseif (jvesm(in).eq.8.or.
     >           jvesm(in).eq.9.or.
     >           jvesm(in).eq.10) then
c
            wallpt(wt,17) = in
c
            if (cprint.eq.3.or.cprint.eq.9) then
c
               write (6,'(''N:'',4(1x,f18.10),2i4)')
     >            rvesm(in,1),zvesm(in,1),
     >            rvesm(in,2),zvesm(in,2),in,jvesm(in)
                  write (6,'(''T:'',4(1x,f18.10))')
     >      wallpt(wt,1) + wallpt(wt,5) * cos(wallpt(wt,8)),
     >      wallpt(wt,2) + wallpt(wt,5) * sin(wallpt(wt,8)),
     >      wallpt(wt,1) + wallpt(wt,6) * cos(wallpt(wt,9)),
     >      wallpt(wt,2) + wallpt(wt,6) * sin(wallpt(wt,9))
               write (6,'(''T:'',4(1x,f18.10),2i4)')
     >            wallpt(wlo,1),
     >            wallpt(wlo,2),
     >            rp(ndso),zp(ndso),wt

               write (6,*)
c
            endif
c
            wt = wt + 1
c
         endif
c
      end do
c
c     Set the indices of the NON-target parts of nimindex
c
      nimindex(1) = 0
      nimindex(ndsin) = 0
      nimindex(ndsin+1) = 0
      nimindex(nds) = 0
c
      if (cprint.eq.3.or.cprint.eq.9) then

         do in = 1,nds
            write(6,'(a,3i5)') 'Wall/Nim index:',in,
     >                     wallindex(in),nimindex(in)
         end do
      endif
c
      return
c
      end
c
c
c
      SUBROUTINE otherrec(NPTS,TE,NE,COEF,opt)
      implicit none
      integer npts,opt
      REAL TE(NPTS), NE(NPTS), COEF(NPTS)
C
c     OTHERREC - Calculate the Gordeev,Janev or NRL  recombination
c                rate coefficients for the array of temperatures
c                passed into the routine. NE is not currently used
c                for these options.
c
c
      integer ik
      real tekev,betan
c
      do ik = 1,npts
c
         if (te(ik).gt.0.0) then
c
            TEKEV = TE(ik)*1.0E-03
            BETAN = 13.6E+00 / TE(ik)
C
c           All these formulae are for coefficients in CGS - convert to
c           MKS by multiplying by 1.0e-6.
c
            if (opt.eq.0) then
c
               coef(ik) = 0.0
c
            elseIF ( opt.EQ.1 ) THEN
c
C              GORDEEV coefficients
c
               IF(TEKEV.LE.4.0) THEN
                  coef(ik)=1.48E-14/(SQRT(TEKEV)*(1.0+43.4*TEKEV))
               ELSE
                  coef(ik)=3.58E-16/TEKEV**1.388
               END IF
c
            ELSE IF ( opt.EQ.2 ) THEN
c
C              JANEV Coefficients (NL=1S ONLY) ?
c
               coef(ik) = 3.92E+00 * 1.00E-14
     >               * ( BETAN**1.5 / ( BETAN + 0.35E+00 ) )
c
            ELSE IF ( opt.EQ.3 ) THEN
c
C              NRL coefficients
c
               coef(ik) = 5.2E-14 * SQRT(BETAN)
     >              * ( 0.43E+00 + 0.500E+00 * ALOG(BETAN)
     >                        + 0.469E+00 / BETAN**0.3333333 )
            endif
c
c           Convert to MKS
c
            coef(ik) = coef(ik) * 1.0e-6
c
c
c            write(6,*) 'rec:',ik,te(ik),ne(ik),coef(ik),tekev,betan
c
         else
c
c           For zero Te <= 0.0  - set the cross-section to 0.0
c
            coef(ik) = 0.0
c
         endif
c
      end do
c
      return
      end
c
c
c
      real function getfracs(ik1,ir1,ik2,ir2,side)
      implicit none
      integer ik1,ir1,ik2,ir2,opt,side
c
c     Commons
c
      include 'params'
      include 'cgeom'
c
c     GETFRACS:
c
c     This routine calculates the distance from the cell centres
c     specified by (ik1,ir1) and (ik2,ir2) to the point of intersection
c     with the polygon walls of the cell (ik1,ir1). The ratio of these
c     distances is also returned. These calculations may be performed
c     with different options that may calculate the fraction in
c     alternative ways. The objective is to obtain an estimate of the
c     proper cross-field distance that would need to be travelled between
c     these two cells. This may be further modified elsewhere by
c     multiplication by KINDS/KOUTDS and COSALI/COSALO - which may
c     be factored into the calculations here if necessary.
c
      integer in,npoly,nsides,ind
      real rp1e,zp1e,rp2e,zp2e,rs1e,zs1e,rs2e,zs2e,rie,zie
      real d1e,d2e,d3e
c
c      real elim
c      parameter (elim=1.0e-4)
c
c     If the same cell is requested return a ratio of 0.5.
c     References to the same cell occur on the boundaries for
c     the wall, the core and the private plasma.
c
c
      if (ik1.eq.ik2.and.ir1.eq.ir2) then
         getfracs = 0.5
         return
      endif
c
c     Does a polygon exist for the primary cell - if not return a
c     ratio of 0.5.
c
      npoly = korpg(ik1,ir1)
c
      if (npoly.gt.0) then
         nsides = nvertp(npoly)
      else
         nsides = 0
      endif
c
c      write (6,*) 'getfracs:',ik1,ir1,ik2,ir2,npoly,nsides
c
      if (npoly.gt.0.and.nsides.ge.3) then
c
c        Polygon exists and has sufficient sides for an intersection.
c
c        Look for intersection point between the line joining the centres
c        of the two grid points and one of the polygon sides.
c
         rp1e = rs(ik1,ir1)
         zp1e = zs(ik1,ir1)
c
         rp2e = rs(ik2,ir2)
         zp2e = zs(ik2,ir2)
c
c        Look only at side passed in on call - unless side=0 -
c        then check for the side actually intersected.
c
         if (side.eq.0) then

            do in = 1,nsides
c
               if (in.eq.1) then

                  rs1e=rvertp(nsides,npoly)
                  zs1e=zvertp(nsides,npoly)
                  rs2e=rvertp(in,npoly)
                  zs2e=zvertp(in,npoly)

               else

                  rs1e=rvertp(in-1,npoly)
                  zs1e=zvertp(in-1,npoly)
                  rs2e=rvertp(in,npoly)
                  zs2e=zvertp(in,npoly)

               endif
c
               call calc_intersect(rp1e,zp1e,rp2e,zp2e,
     >                    rs1e,zs1e,rs2e,zs2e,rie,zie,ind)
c
c               write(6,'(a,4i4,4g15.6)') 'Int-sect:',ind,in,
c     >                            nsides,npoly,
c     >                            rp1e,zp1e,rp2e,zp2e
c               write(6,'(a,6g15.6)') '        :',rie,zie,
c     >                    rs1e,zs1e,rs2e,zs2e
c
c              The intersection routine returns an indicator
c              value of 0 if the intersection actually occurs
c              on the line segments conneting the 2 pairs of
c              points entered.
c
               if (ind.eq.0.or.ind.eq.1) then
                  d1e = sqrt( (rp1e-rie)**2 + (zp1e-zie)**2)
                  d2e = sqrt( (rp2e-rie)**2 + (zp2e-zie)**2)
                  d3e = sqrt( (rp2e-rp1e)**2 + (zp2e-zp1e)**2)

                  getfracs = d1e/(d1e+d2e)
               endif
c
            end do
c
         else
c
c           Calculate intersection only for specific side.
c
            in = side
c
            if (in.eq.1) then

               rs1e=rvertp(nsides,npoly)
               zs1e=zvertp(nsides,npoly)
               rs2e=rvertp(in,npoly)
               zs2e=zvertp(in,npoly)

            else

               rs1e=rvertp(in-1,npoly)
               zs1e=zvertp(in-1,npoly)
               rs2e=rvertp(in,npoly)
               zs2e=zvertp(in,npoly)

            endif
c
            call calc_intersect(rp1e,zp1e,rp2e,zp2e,
     >                    rs1e,zs1e,rs2e,zs2e,rie,zie,ind)
c
c            if (ir1.eq.1.or.ir1.eq.2) then
c
c               write(6,'(a,6i4,4g15.6)') 'Int-sect:',ind,in,
c     >                            ik1,ir1,
c     >                            nsides,npoly,
c     >                            rp1e,zp1e,rp2e,zp2e
c               write(6,'(a,6g15.6)') '        :',rie,zie,
c     >                    rs1e,zs1e,rs2e,zs2e
c            endif
c
c           The intersection routine returns an indicator
c           value of 0 if the intersection actually occurs
c           on the line segments conneting the 2 pairs of
c           points entered.
c
            if (ind.eq.0.or.ind.eq.1) then
               d1e = sqrt( (rp1e-rie)**2 + (zp1e-zie)**2)
               d2e = sqrt( (rp2e-rie)**2 + (zp2e-zie)**2)
               d3e = sqrt( (rp2e-rp1e)**2 + (zp2e-zp1e)**2)
c
c            if (ir1.eq.1.or.ir1.eq.2)
c     >          write(6,'(a,5g15.6,4i4)') 'inf:',d1e,d2e,d3e,
c     >                          d1e+d2e,
c     >                          d1e/(d1e+d2e),
c     >                          ik1,ir1,
c     >                          ik2,ir2
c
c
               getfracs = d1e/(d1e+d2e)

            endif
c
         endif
c
         if (ind.ne.0.and.ind.ne.1) then
c
c        It only comes here if no intersection was found - in
c        which case return a fraction of 0.5 - this could happen
c        due to errors in numerical precision for lines that
c        are almost vertical or horizontal - but not quite.
c
            getfracs = 0.5
c
         else
c
c           Check to see if the fraction should be 0.0 - does the
c           intersection point equal the cell centre of the this cell?
c
c           If the fraction is close to zero or 1.0 set it exactly
c           Cells on a grid should not be so disproportionate
c
            if (getfracs.lt.0.001) then
               getfracs = 0.0
            elseif (getfracs.gt.0.999) then
               getfracs = 1.0
            endif
c
          endif
c
c
c         write (6,*) 'GETFRACS:',in,ik1,ir1,ik2,ir2,getfracs
c         write (6,*) '       R:',rie,rp1e,rp2e,'Z:',zie,zp1e,zp2e
c
      else
         getfracs = 0.5
      endif
c
      return
      end
c
c
c
      subroutine calc_intersect(rp1,zp1,rp2,zp2,
     >                    rs1,zs1,rs2,zs2,ri,zi,ind)
      implicit none
      real rp1,zp1,rp2,zp2
      real rs1,zs1,rs2,zs2
      real ri, zi
      integer ind
c
c     CALC_INTERSECT:
c
c     This routine calculates the intersection point between the
c     two lines defined by [(rp1,zp1),(rp2,zp2)] and
c     [(rs1,zs1),(rs2,zs2)]. It also sets the indicator value
c     to zero if this point lies on both line segments defined
c     by the given end-points and to one if it does not. It also sets
c     ind to two if there is no intersection and to three if the
c     points represent exactly the same lines.
c
c     David Elder           Oct 1, 1996
c
      integer in
      real*8 m1,b1,m2,b2,zsdist,zpdist
      real*8 rp1e,zp1e,rp2e,zp2e,rs1e,zs1e,rs2e,zs2e,rie,zie
c
      logical vertlin1,vertlin2
c
c     Vertical line indicators - i.e. m -> INF
c
      vertlin1 = .false.
      vertlin2 = .false.
c
c     Set to default value for ind - and for interscetion point
c     so that the variables are defined even for those cases
c     where there is no intersection.
c
      ind = 0
      rie = 0.0
      zie = 0.0
c
c     Assign input to local higher precision variables
c
      rp1e = rp1
      zp1e = zp1
      rp2e = rp2
      zp2e = zp2
      rs1e = rs1
      zs1e = zs1
      rs2e = rs2
      zs2e = zs2
c
c     Calculate Line 1 -
c
      if (rp1e.eq.rp2e) then
         vertlin1 = .true.
         m1 = 0.0
         b1 = rp1e
      else
         m1 = (zp1e-zp2e)/(rp1e-rp2e)
         b1 = zp1e - rp1e * m1
      endif
c
c     Calculate Line 2 -
c
      if (rs1e.eq.rs2e) then
         vertlin2 = .true.
         m2 = 0.0
         b2 = rs1e
      else
         m2 = (zs1e-zs2e)/(rs1e-rs2e)
         b2 = zs1e - rs1e * m2
      endif
c
c     Deal with cases
c
c
c     Line 1 vertical
c
      if (vertlin1) then
c
c        with Line 2 vertical
c
         if (vertlin2) then
            ind =2
         else
            rie = b1
            zie = m2 * rie + b2
         endif
c
c        Line 2 vertical
c
      elseif (vertlin2) then
         rie = b2
         zie = m1 * rie + b1
c
c     Same slopes
c
      elseif (m1.eq.m2) then
c
c        Exact same line?
c
         if (b1.eq.b2) then
            ind = 3
         else
            ind = 2
         endif
      else
c
c        Calculate intersection for last case.
c
         rie = (b2-b1)/(m1-m2)
         zie = m1 * rie + b1
      endif
c
      if (ind.eq.0) then
         zpdist = abs(zp1e-zp2e)
         zsdist = abs(zs1e-zs2e)
c
         if ((abs(zie-zp1e).le.zpdist).and.
     >       (abs(zie-zp2e).le.zpdist).and.
     >       (abs(zie-zs1e).le.zsdist).and.
     >       (abs(zie-zs2e).le.zsdist)) then
            ind = 0
         else
            ind = 1
         endif
      endif
c
      ri = rie
      zi = zie
c
c      write(6,*) 'data1:',m1,b1,vertlin1,m2,b2,vertlin2
c      write(6,*) 'data2:',zpdist,zsdist,rie,zie,ind
c
      return
      end
c
c
c
      real function get_sidelen(ik,ir,side,rc)
      implicit none
      integer ik,ir,side,rc
c
c     Commons
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     GET_SIDELEN:
c
c     This routine returns the length of the side of the polygon
c     specified by "side". If the cell does not have that many sides
c     or if the polygon has no sides - this routine will return the
c     distance across the cell centre and will set the value of
c     rc to 1.
c
c     side is defined to be one of four values that define the side
c         of interest in the specific geometry. 
c
c     side = INWARD41 - side facing core or into PP region.
c                       Polygon corners 1,4
c          = OUTWARD23- side facing outer wall or towards separatrix
c                       from inside the PP region. Polygon corners 2,3
c          = DOWN34   - side facing along field lines toward IK=NKS(IR)
c                       target. Polygon corners 4,3
c          = UP12     - side facing along the field lines toward the
c                       IK=1 target.Polygon corners 1,2
c
c     Normal rc is 0
c
c
      integer in,npoly,nsides,ind
      real rp1e,zp1e,rp2e,zp2e,rs1e,zs1e,rs2e,zs2e,rie,zie
c
c     Set default rc
c
      rc = 0
c
c     Does a polygon exist for the primary cell - if not return a
c     length extracted from the cell centre information.
c
      npoly = korpg(ik,ir)
      nsides = nvertp(npoly)
      in = side
c
      if (npoly.ne.0.and.nsides.ge.side) then
c
c        Polygon exists and has sufficient sides
c
c        Look only at side passed in on call.
c
c
         if (side.eq.INWARD41) then

            rs1e=rvertp(4,npoly)
            zs1e=zvertp(4,npoly)
            rs2e=rvertp(1,npoly)
            zs2e=zvertp(1,npoly)

         else

            rs1e=rvertp(in,npoly)
            zs1e=zvertp(in,npoly)
            rs2e=rvertp(in+1,npoly)
            zs2e=zvertp(in+1,npoly)

         endif
c
c         if (side.eq.1) then
c
c            rs1e=rvertp(nsides,npoly)
c            zs1e=zvertp(nsides,npoly)
c            rs2e=rvertp(in,npoly)
c            zs2e=zvertp(in,npoly)
c
c         else
c
c            rs1e=rvertp(in-1,npoly)
c            zs1e=zvertp(in-1,npoly)
c            rs2e=rvertp(in,npoly)
c            zs2e=zvertp(in,npoly)
c
c         endif
c
c        Calculate length ...
c
         get_sidelen = sqrt((rs2e-rs1e)**2 + (zs2e-zs1e)**2)
         rc = 0
c
      else
c
c        Invalid polygon - calculate a length using alternate methods
c
         get_sidelen = kpsiz(ik,ir)
         rc = 1
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
      subroutine wrtdivbg
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     WRTDIVBG: The purpose of this routine is to write out the
c               DIVIMP background plasma in a DIVIMP specific
c               format. This is an expediency for reading and
c               writing just the information DIVIMP
c               needs for running a case.
c
c     Instead of using a unit number - assign a file name. 
c
c
      integer of,ierr
      integer ik,ir,id
c
c     Write to unit 98
c
      parameter(of=98)
c
c     Local file name assigned for output is divimp_plasma.out
c
      OPEN(UNIT=of,FILE='divimp_plasma.out',STATUS='NEW',
     >     ERR=2000,iostat=ierr)
c
c
c     Write Title line
c
      write (of,10) 'DIVIMP BACKGROUND PLASMA:'
      write (of,200) nrs,irsep,nds
      write (of,10) 'KNOTS:'
      write (of,400)  (nks(ir),ir=1,nrs)
c
c     Write out BG quantities - volume
c
c
c     Density - volume and target
c
      write (of,10)  'KNBS:'
      write (of,500) ((knbs(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'KNDS:'
      write (of,500) (knds(id),id=1,nds)
c
c     Te - volume and target
c
      write (of,10)  'KTEBS:'
      write (of,500) ((ktebs(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'KTEDS:'
      write (of,500) (kteds(id),id=1,nds)
c
c     Ti - volume and target
c
      write (of,10)  'KTIBS:'
      write (of,500) ((ktibs(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'KTIDS:'
      write (of,500) (ktids(id),id=1,nds)
c
c     Velocity - volume and target
c
      write (of,10)  'KVHS:'
      write (of,500) ((kvhs(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'KVDS:'
      write (of,500) (kvds(id),id=1,nds)
c
c     Electric Field - volume and target
c
      write (of,10)  'KES:'
      write (of,500) ((kes(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'KEDS:'
      write (of,500) (keds(id),id=1,nds)
c 
c     Blank line at end
c
      write(of,'(/)')
c
      return
c
 2000 continue
c
      write (6,*) 'ERROR WRITING DIVIMP PLASMA BACKGROUND:'
     >       //' FILE EXISTS',ierr
      call pri('ERROR WRITING DIVIMP PLASMA BACKGROUND FILE:'
     >          //' ERR NO = ',ierr)
c
      return
c
c     Formatting
c
  10  format(a,1x,3(1x,i5))
 100  format(a40)
 200  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5)
 300  format('KNOTS:')
 400  format(12i6)
 500  format(6e18.10)
c
      end
c
c
c
      subroutine readdivbg
      implicit none
      include 'params'
      include 'cgeom'
      include 'cedge2d'
c
c     READDIVBG:The purpose of this routine is to read in the
c               DIVIMP background plasma in a DIVIMP specific
c               format. This is an expediency for reading and
c               writing just the information DIVIMP
c               needs for running a case.
c
      integer infile
c
c     Read from unit 98 
c
c     The script that runs DIVIMP copies the specified 
c     background plasma calculated by a previous DIVIMP
c     run "results/<old case name>.bgp" to a file called 
c     divimp_plasma.dat in the DIVIMP run directory. This 
c     file is locally connected to unit 98 and is closed 
c     by the end of the read routine. 
c
c slmod begin - new
      parameter(infile=98)
c
c slmod end
c
      character*120 buffer
      integer ik,ir,id
      integer tmpnrs,tmpnds,tmpirsep
      integer tmpnks(maxnrs)
c slmod begin - new
c...  Open background plasma file:
      OPEN(UNIT=infile,FILE='divimp_plasma.dat',STATUS='OLD',ERR=2000)
c slmod end
c
c     Initialize the CRE2D flag to indicate that data has been
c     read in to the e2d arrays and should be passed to OUT
c
      cre2d = 5
c
c     Zero out background - volume
c
      call rzero(knbs,maxnks*maxnrs)
      call rzero(ktibs,maxnks*maxnrs)
      call rzero(ktebs,maxnks*maxnrs)
      call rzero(kvhs,maxnks*maxnrs)
      call rzero(kes,maxnks*maxnrs)
c
c     target
c
      call rzero(knds,maxnds)
      call rzero(ktids,maxnds)
      call rzero(kteds,maxnds)
      call rzero(kvds,maxnds)
      call rzero(keds,maxnds)
c
c     READ Title line
c
      read (infile,10) buffer
c
      if (buffer(1:6).ne.'DIVIMP') then
         write (6,*) 'ERROR ERROR: NOT A DIVIMP PLASMA FILE'
         write (6,*) buffer
         write (6,*) 'PROGRAM EXITING:'
         write (0,*) 'ERROR ERROR: NOT A DIVIMP PLASMA FILE'
         write (0,*) buffer
         write (0,*) 'PROGRAM EXITING:'
         call prc( 'ERROR ERROR: NOT A DIVIMP PLASMA FILE')
         call prc(buffer)
         call prc('PROGRAM EXITING')
         stop
      endif
c
c
c
 20   read(infile,10,end=1000,err=2000) buffer
c
      if (buffer(1:4).eq.'NRS:') then
         read (buffer,200) tmpnrs,tmpirsep,tmpnds
c
c        Check to see if this matches the current grid.
c
         if (nrs.ne.tmpnrs.or.irsep.ne.tmpirsep.or.nds.ne.tmpnds)
     >          then
c
c           Grid characteristic mismatch - exit program.
c
            write (6,*) 'DIVIMP PLASMA FILE DOES NOT MATCH GRID:'
            write (6,*) 'NRS  :',nrs,tmpnrs
            write (6,*) 'IRSEP:',irsep,tmpirsep
            write (6,*) 'NDS  :',nds,tmpnds
            write (6,*) 'PROGRAM EXITING'
c
            call prc('DIVIMP PLASMA FILE DOES NOT MATCH GRID:')
            call pri2('NRS  :',nrs,tmpnrs)
            call pri2('IRSEP:',irsep,tmpirsep)
            call pri2('NDS  :',nds,tmpnds)
            call prc('PROGRAM EXITING')
c
            stop
c
         end if
c
      elseif(buffer(1:6).eq.'KNOTS:') then

         read (infile,400)  (tmpnks(ir),ir=1,nrs)
c
c        Check to see if knots match
c
c slmod begin
c...bug: SL, 5.7.2004
         do ir = 1, nrs
c         do ik = 1,nrs
c slmod end
            if (nks(ir).ne.tmpnks(ir)) then

               write (6,*) 'DIVIMP PLASMA FILE DOES NOT MATCH GRID:'
               write (6,*) 'IR     :',ir
               write (6,*) 'NKS(IR):',nks(ir),tmpnks(ir)
               write (6,*) 'PROGRAM EXITING'
c
               call prc('DIVIMP PLASMA FILE DOES NOT MATCH GRID:')
               call pri('IR     :',ir)
               call pri2('NKS(IR):',nks(ir),tmpnks(ir))
               call prc('PROGRAM EXITING')
c
               stop
            end if
         end do
c
      elseif (buffer(1:5).eq.'KNBS:') then
c
c        Density - volume
c
         read (infile,500) ((knbs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      elseif (buffer(1:5).eq.'KNDS:') then
c
c        Density - target
c
         read (infile,500) (knds(id),id=1,nds)

      elseif (buffer(1:6).eq.'KTEBS:') then
c
c        Te - volume
c
         read (infile,500) ((ktebs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      elseif (buffer(1:6).eq.'KTEDS:') then
c
c        Te - target
c
         read (infile,500) (kteds(id),id=1,nds)

      elseif (buffer(1:6).eq.'KTIBS:') then
c
c        Ti - volume
c
         read (infile,500) ((ktibs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      elseif (buffer(1:6).eq.'KTIDS:') then
c
c        Ti - target
c
         read (infile,500) (ktids(id),id=1,nds)

      elseif (buffer(1:5).eq.'KVHS:') then
c
c        Velocity - volume
c
         read (infile,500) ((kvhs(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      elseif (buffer(1:5).eq.'KVDS:') then
c
c        Velocity - target
c
         read (infile,500) (kvds(id),id=1,nds)

      elseif (buffer(1:4).eq.'KES:') then
c
c        Electric Field - volume
c
         read (infile,500) ((kes(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      elseif (buffer(1:5).eq.'KEDS:') then
c
c        Electric Field - target
c
         read (infile,500) (keds(id),id=1,nds)

      endif
c
c     Loop back for continued reading
c
      goto 20
c
c     End of file - 1000 - assumed normal exit
c
 1000 continue
c
c     Copy DIVIMP plasma input into fluid code variables since these
c     are not used when a DIVIMP plasma is loaded. This will allow
c     plotting of an evolving DIVIMP background solution.
c
c slmod begin - new
c...  Only write to the fluid code arrays if a fluid code solution is
c     not being read in for reference:
      if (cre2d.eq.0) then
         do ir = 1,nrs
            do ik = 1,nks(ir)
               e2dnbs(ik,ir) = knbs(ik,ir)
               e2dtebs(ik,ir)= ktebs(ik,ir)
               e2dtibs(ik,ir)= ktibs(ik,ir)
               e2dvhs(ik,ir) = kvhs(ik,ir)
               e2des(ik,ir)  = kes(ik,ir)

c               write (6,'(a,2i4,5g12.4)') 'FC:',ik,ir,e2dnbs(ik,ir),
c     >                     e2dtebs(ik,ir),e2dtibs(ik,ir),e2dvhs(ik,ir),
c     >                     e2des(ik,ir)

            end do
         end do
      endif

      CLOSE (infile)

      CALL OutputData(85,'after loading plasma')
c
c      do ir = 1,nrs
c         do ik = 1,nks(ir)
c            e2dnbs(ik,ir) = knbs(ik,ir)
c            e2dtebs(ik,ir)= ktebs(ik,ir)
c            e2dtibs(ik,ir)= ktibs(ik,ir)
c            e2dvhs(ik,ir) = kvhs(ik,ir)
c            e2des(ik,ir)  = kes(ik,ir)
cc
cc            write (6,'(a,2i4,5g12.4)') 'FC:',ik,ir,e2dnbs(ik,ir),
cc     >                  e2dtebs(ik,ir),e2dtibs(ik,ir),e2dvhs(ik,ir),
cc     >                  e2des(ik,ir)
cc
c         end do
c      end do
c slmod end
c
      return
c
 2000 continue
c
      write (0,*) 'ERROR READING IN DIVIMP PLASMA BACKGROUND:'
      write (6,*) 'ERROR READING IN DIVIMP PLASMA BACKGROUND:'
      call prc('ERROR READING IN DIVIMP PLASMA BACKGROUND:')
      stop
c
c     Formatting
c
  10  format(a)
 100  format(a40)
 200  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5)
 400  format(12i6)
 500  format(6e18.10)
c
      end
c
c
c
      subroutine wrtdivaux(nizs)
      implicit none
c
      integer nizs
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cioniz'
      include 'dynam1'
      include 'dynam3'
c
c     WRTDIVAUX:The purpose of this routine is to write out 
c               additional data from a DIVIMP run. This data
c               includes any desired impurity or simulation
c               related result that may need to be read in as
c               input to a future DIVIMP run. 
c
c               The associated read routine can be used to
c               access specific data from this file as well 
c               as data from the matching divimp background plasma
c               file. 
c
c     Instead of using a unit number - assign a file name. 
c
c
c     Write to unit 98
c
      integer of
      parameter(of=98)
c
c     Ionization potential energies of Carbon
c
      real Iiz(0:5)
      data Iiz /11.0, 24.0, 48.0, 64.0, 392.0, 490.0/  
c
      integer ik,ir,iz,ierr
      real tpowls(maxnks,maxnrs)
      real tcooliz(maxnks,maxnrs)
c
c     Initialize array that will contain total impurity emitted power in 
c     each cell.
c
      call rzero(tpowls,maxnks*maxnrs)
c
c     Sum over charge states and multiply by absfac 
c     to get total emitted power.
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            do iz = 0,nizs
               tpowls(ik,ir) = tpowls(ik,ir)+powls(ik,ir,iz)*absfac
            end do
         end do
      end do
C
C     Calculate cooling term due to ionization of impurities
c     ONLY for Carbon at the present time. Base the ionization and
c     recombination rate data on the density and reaction rates since it 
c     will give better statistics than using the actual DIVIMP ionization
c     results.
C
      if (cion.eq.6) then       
c
         call rzero(tcooliz,maxnks*maxnrs)
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               do iz = 0,nizs-1
                  tcooliz(ik,ir) = tcooliz(ik,ir) +
     >                Iiz(iz)*(ddlims(ik,ir,iz) / kfizs(ik,ir,iz)
     >                      -ddlims(ik,ir,iz+1) / kfrcs(ik,ir,iz+1))
               end do 
            end do
         end do
c
c        Scale by the absolute factor and convert from eV/m3 to W/m3 
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               tcooliz(ik,ir) = tcooliz(ik,ir) * ( ech * absfac)
            end do
         end do
c
      endif
c 
c     Local file name assigned for output is divimp_aux_data.out
c
      OPEN(UNIT=of,FILE='divimp_aux_data.out',STATUS='NEW',
     >     ERR=2000,iostat=ierr)
c
c     Write Title line
c
c     Note - the grid information is included so that the data from this
c     file can be verified against the grid that is in use by the current
c     simulation.   
c
      write (of,10) 'DIVIMP RUN ADDITIONAL OUTPUT:'
      write (of,200) nrs,irsep,nds,nizs
      write (of,10) 'KNOTS:'
      write (of,400)  (nks(ir),ir=1,nrs)
c
c     Write out Additional data quantities 
c
c     ABSFAC  - absolute scaling factor
c
      write (of,10)  'ABSFAC:'
      write (of,500) absfac
c
c     TOTAL POWLS: Total impurity radiation density
c                  in each cell
c        
      write (of,10)  'TPOWLS:'
      write (of,500) ((tpowls(ik,ir),ik=1,nks(ir)),ir=1,nrs)
      write (of,10)  'TCOOLIZ:'
      write (of,500) ((tcooliz(ik,ir),ik=1,nks(ir)),ir=1,nrs)
c
c     POWLS: Impurity radiation density for each charge state 
c     in each cell - scaled to 1 part/sec
c        
c      write (of,10)  'POWLS:'
c      write (of,500) (((powls(ik,ir,iz),
c     >            ik=1,nks(ir)),ir=1,nrs),iz=-1,nizs)
c
c     Add other tagged input in the same fashion.
c       
c 
c     Blank line at end
c
      write(of,'(/)')
c
      return
c
 2000 continue
c
      write (6,*) 'ERROR WRITING DIVIMP AUXILIARY DATA FILE:'
     >       //' FILE EXISTS',ierr
      call pri('ERROR WRITING DIVIMP AUXILIARY DATA FILE:'
     >          //' ERR NO = ',ierr)
c
      return
c
c     Formatting
c
  10  format(a,1x,3(1x,i5))
 100  format(a40)
 200  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5,'NIZS:',i5)
 300  format('KNOTS:')
 400  format(12i6)
 500  format(6e18.10)
c
      end
c
c
c
      subroutine readdivaux(tag,data_array,maxk,maxr,minz,maxz)
      implicit none
      include 'params'
      include 'cgeom'
      include 'cedge2d'
c
      character*(*) tag
      integer maxk,maxr,minz,maxz
      real data_array(maxk,maxr,minz:maxz)  
c
c     READDIVAUX:The purpose of this routine is to read in
c                specific data sets from the divimp_aux_data.dat
c                file and the divimp_plasma.dat file. 
c       
c                The tag of the item to be looked for is passed
c                into this routine along with the name of the 
c                variable where this will be stored. This
c                variable can have as many as 3 dimensions.
c
c                This routine is designed to allow for the selective
c                reading of specific data from the results files
c                of a previous DIVIMP run. Unfortunately, although
c                the raw data file contains the information of 
c                interest - the file itself is not setup in a 
c                useful fashion for extracting individually tagged
c                data sets. 
c
      integer infile
c
c     Read from unit 98 
c
c     The script that runs DIVIMP copies the specified 
c     background plasma calculated by a previous DIVIMP
c     run "results/<old case name>.bgp" to a file called 
c     divimp_plasma.dat in the DIVIMP run directory. This 
c     file is locally connected to unit 98 and is closed 
c     by the end of the read routine. 
c     
c     The auxiliary results file from a previous DIVIMP run is 
c     copied to divimp_aux_data.dat
c
c
      parameter(infile=98)
c
      character*120 buffer,filename
      integer ik,ir,id,iz
      integer tmpnrs,tmpnds,tmpirsep,tmpnizs
      integer tmpnks(maxnrs)
      integer len,lent,lenf,lenstr
      integer fileid
      external lenstr
c
c     Based on the value of the tag - decide which file to open
c
      lent= lenstr(tag)
c
      if (tag(1:lent).eq.'ABSFAC:'.or.
c     >    tag(1:lent).eq.'POWLS:'.or.
     >    tag(1:lent).eq.'TPOWLS:') then 
         filename = 'divimp_aux_data.dat'
         fileid = 1
      elseif (tag(1:lent).eq.'KNBS:'.or.  
     >        tag(1:lent).eq.'KNDS:'.or. 
     >        tag(1:lent).eq.'KTEBS:'.or. 
     >        tag(1:lent).eq.'KTEDS:'.or. 
     >        tag(1:lent).eq.'KTIBS:'.or. 
     >        tag(1:lent).eq.'KTIDS:'.or. 
     >        tag(1:lent).eq.'KVHS:'.or. 
     >        tag(1:lent).eq.'KVDS:'.or. 
     >        tag(1:lent).eq.'KES:'.or. 
     >        tag(1:lent).eq.'KEDS:') then 
         filename = 'divimp_plasma.dat'
         fileid=2
      else
c
c        Unrecognized tag
c
         write(6,*) 'UNRECOGNIZED TAG IN READDIVAUX:'
     >                     //tag(1:lent)
         write(6,*) 'Routine Exiting'
         write(0,*) 'UNRECOGNIZED TAG IN READDIVAUX:'
     >                     //tag(1:lent)
         write(6,*) 'Routine Exiting'
         return  
      endif
c
      len = lenstr(filename)
      OPEN(UNIT=infile,FILE=filename(1:len),STATUS='OLD',ERR=2000)
c
c
c     READ Title line
c
      read (infile,10) buffer
c
      if (buffer(1:6).ne.'DIVIMP') then
         write (6,*) 'ERROR ERROR: NOT A DIVIMP FILE'
         write (6,*) buffer
         write (6,*) 'READ DIV AUX EXITING:'
c
         write (0,*) 'ERROR ERROR: NOT A DIVIMP FILE'
         write (0,*) buffer
         write (0,*) 'READ DIV AUX EXITING:'
         return
      endif
c
c     Verify grid characteristics
c
 20   read(infile,10,end=1000,err=2000) buffer
c
      if (buffer(1:4).eq.'NRS:') then
c
         if (fileid.eq.1) then 
            read (buffer,210) tmpnrs,tmpirsep,tmpnds,tmpnizs
         elseif (fileid.eq.2) then 
            read (buffer,200) tmpnrs,tmpirsep,tmpnds
            tmpnizs = minz
         endif
c
c        Check to see if this matches the current grid.
c
         if (nrs.ne.tmpnrs.or.irsep.ne.tmpirsep.or.nds.ne.tmpnds)
     >          then
c
c           Grid characteristic mismatch - exit routine.
c
            write (6,*) 'DIVIMP FILE DOES NOT MATCH GRID:'
            write (6,*) 'NRS  :',nrs,tmpnrs
            write (6,*) 'IRSEP:',irsep,tmpirsep
            write (6,*) 'NDS  :',nds,tmpnds
            write (6,*) 'READ DIV AUX EXITING'
c
            write (0,*) 'DIVIMP FILE DOES NOT MATCH GRID:'
            write (0,*) 'NRS  :',nrs,tmpnrs
            write (0,*) 'IRSEP:',irsep,tmpirsep
            write (0,*) 'NDS  :',nds,tmpnds
            write (0,*) 'READ DIV AUX EXITING'
c
            return
c
         end if
c
      elseif(buffer(1:6).eq.'KNOTS:') then

         read (infile,400)  (tmpnks(ir),ir=1,nrs)
c
c        Check to see if knots match
c
         do ik = 1,nrs
            if (nks(ir).ne.tmpnks(ir)) then

               write (6,*) 'DIVIMP FILE DOES NOT MATCH GRID:'
               write (6,*) 'IR     :',ir
               write (6,*) 'NKS(IR):',nks(ir),tmpnks(ir)
               write (6,*) 'READ DIV AUX EXITING'
c
               write (0,*) 'DIVIMP FILE DOES NOT MATCH GRID:'
               write (0,*) 'IR     :',ir
               write (0,*) 'NKS(IR):',nks(ir),tmpnks(ir)
               write (0,*) 'READ DIV AUX EXITING'
c
               return
            end if
         end do
c
      elseif (buffer(1:lent).eq.tag(1:lent)) then
c
c        Read in arbitrary tagged field
c
         if (minz.ne.maxz) then  
            read (infile,500) (((data_array(ik,ir,iz),
     >             ik=1,nks(ir)),ir=1,nrs),iz=minz,tmpnizs)
         else 
            read (infile,500) ((data_array(ik,ir,1),
     >             ik=1,nks(ir)),ir=1,nrs)
         endif 
c
      endif  
c
c     Loop back for continued reading
c
      goto 20
c
c     End of file - 1000 - assumed normal exit
c
 1000 continue
c
      CLOSE (infile)

c
      return
c
 2000 continue
c
      write(6,*) 'ERROR READING IN DIVIMP DATA FILE: READDIVAUX:'
      call prc('ERROR READING IN DIVIMP DATA FILE; READDIVAUX:')
      return
c
c     Formatting
c
  10  format(a)
 100  format(a40)
 200  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5)
 210  format('NRS:',i5,'IRSEP:',i5,'NDS:',i5,'NIZS:',i5)
 400  format(12i6)
 500  format(6e18.10)
c
      end
c
c
c
      subroutine calcapp_fgradmod
      implicit none
      include  'params'
      include  'cgeom'
      include  'comtor'
      include  'cioniz'
c
c     CALCAPP_FGRADMOD:
c
c     The purpose of this routine is to calculate and apply
c     a modification
c     or correction factor that is applied to the temperature
c     gradient forces in an attempt to account for kinetic
c     effects in the application of those foeces.
c
c     fgrad option 0: No correction is applied - all factors
c                     are set equal to 1.0
c
c     fgrad option 1: Uses an adhoc formula used in UEDGE to
c                     account for these effects - based
c                     on scale lengths and mean free paths
c                     with a specifiable parameter.
c
c     fgrad option 2: Uses an adhoc formula used in UEDGE to
c                     account for these effects - based
c                     on scale lengths and mean free paths
c                     with a specifiable parameter. When the
c                     ion mean free path is greater than the
c                     distance to the target then FIG is set
c                     to zero.
c
c     fgrad option 3: Uses an adhoc formula taken from Garching
c                     for limiting the size of the temperature
c                     gradient.
c
c     David Elder, Nov 6, 1996
c
c
c     Local variables
c
      integer ik,ir,id,in
c
      real fact
c
      real fgradmode(maxnks,maxnrs),fgradmodi(maxnks,maxnrs)
c
      real lgrad(maxnks,maxnrs)
c
      real kpbs(maxnks,maxnrs),kpds(maxnds)
c
      real gradnpara(maxnks,maxnrs),gradtepara(maxnks,maxnrs)
      real gradppara(maxnks,maxnrs),gradtipara(maxnks,maxnrs)
c
      real lgradn(maxnks,maxnrs),lgradte(maxnks,maxnrs)
      real lgradp(maxnks,maxnrs),lgradti(maxnks,maxnrs)
c
      real limfp(maxnks,maxnrs),lemfp(maxnks,maxnrs)
c
      real lmfp(maxnks,maxnrs)
c
      call rzero(limfp,maxnks*maxnrs)
      call rzero(lemfp,maxnks*maxnrs)
c
c     Option 0: Off - no changes
c
      if (fgradopt.eq.0) then
c
c        Do NOT change the temperature gradients.
c
         call rinit(fgradmode,maxnrs*maxnks,1.0)
         call rinit(fgradmodi,maxnrs*maxnks,1.0)
c
c     Option 1, 2 and 4: UEDGE variations
c
      elseif (fgradopt.eq.1.or.fgradopt.eq.2.or.fgradopt.eq.4) then
c
c        Calculate the modifier factors based on the UEDGE formula
c        Correction is the same for FEG and FIG.
c
c
c        Calculate gradients
c
         call calc_grad(gradnpara,knbs,knds)
         call calc_grad(gradtepara,ktebs,kteds)
         call calc_grad(gradtipara,ktibs,ktids)
c
c        Calculate the STATIC pressure in each cell and at the targets.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               kpbs(ik,ir) = knbs(ik,ir) *
     >                          ech * (ktebs(ik,ir)+ktibs(ik,ir))
            end do
         end do
c
         do id = 1,nds
            kpds(id) = knds(id) * ech * (kteds(id)+ktids(id))
         end do
c
         call calc_grad(gradppara,kpbs,kpds)
c
c        Calculate the scale lengths
c
         call calc_scale(gradnpara,lgradn,knbs,nrs,nks)
         call calc_scale(gradtepara,lgradte,ktebs,nrs,nks)
         call calc_scale(gradtipara,lgradti,ktibs,nrs,nks)
         call calc_scale(gradppara,lgradp,kpbs,nrs,nks)
c
c        Calculate the Lmfp value - sum of mean free paths for
c        electrons and ions
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
c              Calculate the electron and ion mean free paths
c

               if (knbs(ik,ir).le.0.0) then
                  limfp(ik,ir) = 0.0
                  lemfp(ik,ir) = 0.0
               else
c
                  limfp(ik,ir)=ktibs(ik,ir)**2 * (1.0e16/knbs(ik,ir))
                  lemfp(ik,ir)=ktebs(ik,ir)**2 * (1.0e16/knbs(ik,ir))
c
               endif
c
               lmfp(ik,ir) = limfp(ik,ir) + lemfp(ik,ir)
c
            end do
         end do
c
c        Calculate the correction modifier.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
c              Determine minimum scale length for cell.
c
               lgrad(ik,ir) = min(lgradn(ik,ir),lgradp(ik,ir),
     >                     lgradte(ik,ir),lgradti(ik,ir))
c
c              Calculate correction
c
               if (lgrad(ik,ir).eq.0.0) then
                  fgradmodi(ik,ir) = 1.0
               else
                  fgradmodi(ik,ir) = 1.0 /
     >               (1.0 + fgradfact*(lmfp(ik,ir)/lgrad(ik,ir))**2)

               endif
c
               fgradmode(ik,ir) = fgradmodi(ik,ir)
c
            end do
         end do
c
c     Option 3: Garching B2 formula
c
      elseif (fgradopt.eq.3) then
c
c        Calculate the Limfp and lemfp values
c
         fact = qtim*qtim*emi/crmi
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
c              Calculate the electron and ion mean free paths
c
               if (knbs(ik,ir).le.0.0) then
                  limfp(ik,ir) = 0.0
                  lemfp(ik,ir) = 0.0
               else
c
                  limfp(ik,ir)=ktibs(ik,ir)**2 *
     >                        (fgradfact*1.5e16/knbs(ik,ir))
                  lemfp(ik,ir)=ktebs(ik,ir)**2 *
     >                        (fgradfact*1.5e16/knbs(ik,ir))
c
               endif
            end do
         end do
c
c        Calculate the correction modifier.
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
c
c              Determine minimum scale length for cell.
c              Include kfigs in the calculation because this
c              is a correction factor that will multiply kfigs.
c              The abs() is necessary so that the sign of
c              kfigs,kfegs does not change with the value.
c
c
c              FIG
c
               if (limfp(ik,ir).ne.0.0.and.kfigs(ik,ir).ne.0.0) then
c
                  lgrad (ik,ir) = 0.3*ktibs(ik,ir)/limfp(ik,ir)
c
                  fgradmodi(ik,ir) =
     >                min(abs(kfigs(ik,ir)/fact),
     >                    0.3*ktibs(ik,ir)/limfp(ik,ir))
     >                         / abs(kfigs(ik,ir)/fact)
c
               else
c
                  lgrad (ik,ir) = 0.0
c
                  fgradmodi(ik,ir) = 1.0
c
               endif
c
c              FEG
c
               if (lemfp(ik,ir).ne.0.0.and.kfegs(ik,ir).ne.0.0) then
c
                  fgradmode(ik,ir) =
     >                min(abs(kfegs(ik,ir)/fact),
     >                    0.3*ktebs(ik,ir)/lemfp(ik,ir))
     >                         / abs(kfegs(ik,ir)/fact)
c
               else
c
                  fgradmode(ik,ir) = 1.0
c
               endif

            end do
         end do
c
      endif
c
c     Apply the modifier to kfegs and kfigs
c
c
c     Multiply the force arrays KFEGS and KFIGS by the kinetic
c     correction factor.
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           FeG
c
            if (fgradopt.eq.4.and.
     >         ( min(kss(ik,ir),ksmaxs(ir)-kss(ik,ir))
     >               .le.lemfp(ik,ir) )) then
               kfegs(ik,ir) = 0.0
            else
               kfegs(ik,ir) = kfegs(ik,ir)*fgradmode(ik,ir)
            endif
c
c           FiG
c
            if ((fgradopt.eq.2.or.fgradopt.eq.4).and.
     >         ( min(kss(ik,ir),ksmaxs(ir)-kss(ik,ir))
     >               .le.limfp(ik,ir) )) then
               kfigs(ik,ir) = 0.0
            else
               kfigs(ik,ir) = kfigs(ik,ir)*fgradmodi(ik,ir)
            endif
c
         end do
      end do
c
c     Calculate the force of friction modifier if friction option
c     4 is in use.
c
c     Initialize friction modifier to 1.0 everywhere
c
      call rinit (kfssmod,maxnks*maxnrs,1.0)
c
c     Modify it for kinetic corrections - this array is used external
c     to this routine.
c
      if (cioptc.eq.4) then
c
         do ir = 1,nrs
c
            do ik = 1,nks(ir)

               if (min(kss(ik,ir),ksmaxs(ir)-kss(ik,ir))
     >               .le.limfp(ik,ir) ) then
                  kfssmod(ik,ir) = 0.0
c
c               else
c                  kfssmod(ik,ir) = 1.0
c
               endif
c
            end do
c
         end do
c
      endif
c
c
c     Print out values for debugging purposes
c
      if ((cprint.eq.3..or.cprint.eq.9).and.
     >    (fgradopt.eq.1.or.fgradopt.eq.2.or.fgradopt.eq.4)) then
c
         CALL PRRMATDIV(gradnpara,MAXNKS,nks(irsep),NRS,6,'GRAD-N')
         CALL PRRMATDIV(gradtepara,MAXNKS,nks(irsep),NRS,6,'GRAD-TE')
         CALL PRRMATDIV(gradtipara,MAXNKS,nks(irsep),NRS,6,'GRAD-TI')
         CALL PRRMATDIV(gradppara,MAXNKS,nks(irsep),NRS,6,'GRAD-P')
c
         CALL PRRMATDIV(kpbs,MAXNKS,nks(irsep),NRS,6,'PRESSURE')
c
         CALL PRRMATDIV(lgradn,MAXNKS,nks(irsep),NRS,6,'LGRAD-N')
         CALL PRRMATDIV(lgradte,MAXNKS,nks(irsep),NRS,6,'LGRAD-TE')
         CALL PRRMATDIV(lgradti,MAXNKS,nks(irsep),NRS,6,'LGRAD-TI')
         CALL PRRMATDIV(lgradp,MAXNKS,nks(irsep),NRS,6,'LGRAD-P')
c
         CALL PRRMATDIV(lmfp,MAXNKS,nks(irsep),NRS,6,'LMFP-TOTAL')
         CALL PRRMATDIV(lgrad,MAXNKS,nks(irsep),NRS,6,'LGRAD')
         CALL PRRMATDIV(lemfp,MAXNKS,nks(irsep),NRS,6,'ELECTRON-MFP')
         CALL PRRMATDIV(limfp,MAXNKS,nks(irsep),NRS,6,'ION-MFP')
         CALL PRRMATDIV(ktibs,MAXNKS,nks(irsep),NRS,6,'ION-TI')
         CALL PRRMATDIV(kfigs,MAXNKS,nks(irsep),NRS,6,'ION-FIGS')
         CALL PRRMATDIV(fgradmode,MAXNKS,nks(irsep),NRS,6,'FGRAD-MOD E')
         CALL PRRMATDIV(fgradmodi,MAXNKS,nks(irsep),NRS,6,'FGRAD-MOD I')
c
      endif
c
c
      return
      end
c
c
c
      subroutine calc_grad(valgrad,val,valtarg)
      implicit none
      include 'params'
      include 'cgeom'
c
      real valgrad(maxnks,maxnrs)
      real val(maxnks,maxnrs)
      real valtarg(maxnds)
      real endval
c
c     CALC_GRAD:
c
c     This subroutine calculates the parallel to the field
c     line gradients of the parameter passed in. It uses the
c     values at the targets passed in the valtarg array to
c     help define the gradients in the first cell.
c
c
      integer ik,ir,id
      real dist1, dist2
c
      do ir = 1, nrs
         do ik = 1,nks(ir)

            if (ik.eq.1) then
c
c              Distinguish between core and SOL
c
               dist2 = (kss(ik+1,ir)-kss(ik,ir))
c
               if (ir.lt.irsep) then
c
                  dist1 = (kss(nks(ir),ir) - kss(nks(ir)-1,ir))
                  endval= val(nks(ir)-1,ir)
c
               else
c
c                 IK=1 target (second)
c
                  id = idds(ir,2)

                  dist1 = (kss(ik,ir)-ksb(ik-1,ir))
                  endval = valtarg(id)
c
               endif
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((val(ik+1,ir)-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-endval)
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-endval)
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (val(ik+1,ir)-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            elseif (ik.eq.nks(ir)) then
c
c              Distinguish between core and SOL
c
               dist1 = (kss(ik,ir)-kss(ik-1,ir))

               if (ir.lt.irsep) then
c
                  dist2 = (kss(2,ir) - kss(1,ir))
                  endval= val(2,ir)
c
               else
c
c                 IK=NKS(IR) target (first)
c
                  id = idds(ir,1)
c
                  dist2 = ksb(ik,ir)-kss(ik,ir)
                  endval = valtarg(id)
c
               endif
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((endval-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-val(ik-1,ir))
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-val(ik-1,ir))
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (endval-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            else
c
c              General case - give average of slope to next and last
c              cells.
c
               dist2 = (kss(ik+1,ir)-kss(ik,ir))
               dist1 = (kss(ik,ir)-kss(ik-1,ir))
c
               if (dist1.ne.0.0.and.dist2.ne.0.0) then
                  valgrad(ik,ir) = ((val(ik+1,ir)-val(ik,ir))
     >                             /dist2
     >                           + (val(ik,ir)-val(ik-1,ir))
     >                             /dist1  )/2.0
               elseif (dist1.ne.0.0) then
                  valgrad(ik,ir) = (val(ik,ir)-val(ik-1,ir))
     >                                /dist1
               elseif (dist2.ne.0.0) then
                  valgrad(ik,ir) = (val(ik+1,ir)-val(ik,ir))
     >                             /dist2
               else
                  valgrad(ik,ir) = 0.0
               endif
c
            endif
c
         end do
      end do
c
      return
      end
c
c
c
      subroutine calc_scale(valgrad,valscale,val,nrs,nks)
      implicit none
      include 'params'
      real valgrad(maxnks,maxnrs)
      real val(maxnks,maxnrs)
      real valscale(maxnks,maxnrs)
      integer nrs
      integer nks(maxnrs)
c
c     CALC_SCALE:
c
c     This subroutine calculates the cell by cell scale
c     lengths of the quantities passed in.
c
      integer ik,ir
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           Calculate scale lengths - an error condition will
c           result in a very large scale length being assigned.
c
            if (valgrad(ik,ir).eq.0.0.or.val(ik,ir).eq.0.0) then
               valscale(ik,ir) = hi
            else
               valscale(ik,ir) = abs(
     >                       1.0/((1.0/val(ik,ir))*valgrad(ik,ir)))
            endif
         end do
      end do
c
      return
      end
c
c
c
      integer function sfind(s,ir)
      implicit none
      real s
      integer ir
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     SFIND: This routine finds the IK index of the cell on ring IR
c            that contains the given S-position. Depending on the
c            value of PDOPT - it will make different assumptions on
c            the locations of the cell boundaries. PDOPT=1 is taken to
c            imply that the ksb arrays are properly filled and
c            available. PDOPT=0 will simply find the closest KSS value.
c
c            David Elder, March 5, 1997
c
      integer ik,ikfound
c
      if (pdopt.eq.0) then
c
c        Find closest KSS
c
         do ik = 1,nks(ir)-1

c
c           Before first cell
c
            if (ik.eq.1.and.s.lt.kss(ik,ir)) then
               sfind = 1
               return
c
c           After last cell
c
            elseif (ik.eq.nks(ir).and.s.gt.kss(ik,ir)) then
               sfind = nks(ir)
               return
c
c           In one of the intervening cells
c
            elseif (s-kss(ik,ir).lt.(kss(ik+1,ir)-kss(ik,ir))) then
c
c              Check to make sure that the kss array is properly setup
c
               if ((kss(ik+1,ir)-kss(ik,ir)).gt.0.0) then
                  if ((s-kss(ik,ir)).gt.
     >               (0.5*(kss(ik+1,ir)-kss(ik,ir)))) then
                     sfind = ik +1
                     return
                  else
                     sfind = ik
                     return
                  endif
               endif
            endif

         end do

      elseif (pdopt.eq.1) then

         do ik = 1,nks(ir)
c
c           Check to see if S is with the bounds of each bin
c
            if ((s-ksb(ik-1,ir)).lt.(ksb(ik,ir)-ksb(ik-1,ir))) then
               sfind = ik
               return
            endif

         end do

      endif

c
c     Error condition - the code should never reach here
c
      call prr ('ERROR IN FUNCTION "SFIND" : ',s)
      write (6,*) 'ERROR IN FUNCTION "SFIND" : ',s
      sfind = nks(ir)/2
c
      return
      end
c
c
c
      subroutine set_ikvals(ir,ikstart,ikend,ikopt)
      implicit none
      integer ir,ikstart,ikend,ikopt
      include 'params'
      include 'cgeom'
c
c     SET_IKVALS: The purpose of this routine is to
c                 set the IK values to be used in the
c                 looping construct in tje calling routine
c                 based on the oprion specified by IKOPT
c
c                 IKOPT = 1 - first half of ring (OUTER for JET)
c                 IKOPT = 2 - second half of ring (INNER for JET)
c                 IKOPT = 3 - whole ring
c
c                 David Elder,   July, 1997
c
      if (ikopt.eq.1) then
c
c        First half of ring
c
         ikstart = 1
         ikend   = ikmids(ir)
c
      elseif (ikopt.eq.2) then
c
c        Second half of ring
c
         ikstart = ikmids(ir) + 1
         ikend   = nks(ir)
c
      elseif (ikopt.eq.3) then
c
c        Whole ring
c
         ikstart = 1
         ikend   = nks(ir)
      endif
c
c     ADJUST ENDS FOR CORE RINGS - which are repeated points
c
      if (ir.lt.irsep.and.ikend.eq.nks(ir)) then
c
         ikend = ikend -1
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
      subroutine reade2daux
      implicit none
c
c     Read E2D auxilliary input file
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'reader'
c
c     Variables to hold edge2d values
c
      include 'cedge2d'
c
c
c     Local variables
c
      integer ik,ir,in,rc,ikd
      integer core_rings,sol_rings,pp_rings,nrings
c
      real flux,polysidelen,sidelen,sider,field_ratio
      real side_orthog,dummy(maxnrs),sidez,atan2c
      external polysidelen,atan2c
      real theta1,theta2
c
c     Initializiation
c
      call rzero (e2dgdown,maxnrs*maxnks)
      call rzero (e2dgpara,maxnrs*maxnks)
      call rzero (e2dpedown,maxnrs*maxnks)
      call rzero (e2dpepara,maxnrs*maxnks)
      call rzero (e2dpidown,maxnrs*maxnks)
      call rzero (e2dpipara,maxnrs*maxnks)
c
      core_rings = irsep -1
      sol_rings  = irwall - irsep + 1
      pp_rings   = nrs-irtrap + 1
c
      nrings = sol_rings + max(core_rings,pp_rings)
c
      write (6,*) 'READAUX:',cprint ,
     >             core_rings,sol_rings,pp_rings,nrings
c
c     Read an auxilliary file from unit 12 - this code
c     is currently customized to read one specific set of
c     EDGE2D data in a formatted table.
c
 550  read(12,'(a132)',end=590) buffer
c
c      write(0,*) 'INT:',buffer
c

      if (buffer(1:6).eq.'#PFLXD') then
c
         call readauxarray(e2dgdown,core_rings,sol_rings,pp_rings,
     >                     nrings,'E2D-GDOWN')
c
      elseif (buffer(1:7).eq.'#QIFLXD') then
c
         call readauxarray(e2dpidown,core_rings,sol_rings,pp_rings,
     >                     nrings,'E2D-PI-DOWN')
c
      elseif (buffer(1:7).eq.'#QEFLXD') then
c
         call readauxarray(e2dpedown,core_rings,sol_rings,pp_rings,
     >                     nrings,'E2D-PE-DOWN')
c
      endif

      goto 550

 590  continue

c
c
c     Convert to fluxes/m2 from particles ...
c
c     Keep in mind - these are cell boundary quantities and first cell
c     is zero -
c
c      write(6,*) 'IK=1 IS A VIRTUAL CELL'
c      write(6,100)
c100   format(6x,2x,'IK',1x,'IKD',2x,'IR',2x, 'GDOWN-PART',1x,
c     >       2x,'GDOWN-ORTH',1x,2x,'ORTH-CORR',2x,
c     >       1x,'GDOWN-FINAL',1x,3x,'GE2D-BND',2x,
c     >       2x,'GE2D-CENT',2x)
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
            ikd = ik +1
c
            in = korpg(ikd,ir)
c
c           IPP/08 Krieger - ensure index of nvertp is not zero
            if (in.ne.0.and.nvertp(max(1,in)).ne.0) then
c
               sidelen = polysidelen(ikd,ir,DOWN12,rc)
c
               sider = (rvertp(1,in) + rvertp(2,in)) /2.0
               sidez = (zvertp(1,in) + zvertp(2,in)) /2.0
c
               field_ratio = (kbfs(ik,ir) + kbfs(ik+1,ir)) /2.0
c
               side_orthog = 1.0
c
               theta1 = atan2c(zvertp(2,in)-zvertp(1,in),
     >                         rvertp(2,in)-rvertp(1,in))-PI/2.0
               theta2 = atan2c(zs(ikd,ir)-sidez,
     >                         rs(ikd,ir)-sider)
c
               side_orthog = cos(theta2-theta1)
c
c              Particle flux
c
               e2dgpara(ikd,ir)  = e2dgdown(ikd,ir)
     >                / (sidelen* 2.0 * PI * sider)
     >                * field_ratio/side_orthog
c
c              Ion power flux
c
               e2dpipara(ikd,ir) = e2dpidown(ikd,ir)
     >                / (sidelen* 2.0 * PI * sider)
     >                * field_ratio/side_orthog
c
c              Electron power flux
c
               e2dpepara(ikd,ir) = e2dpedown(ikd,ir)
     >                / (sidelen* 2.0 * PI * sider)
     >                * field_ratio/side_orthog
c
c
c               write (6,'(a,3i4,1p,10(1x,g12.5))') 'G-DOWN',ik,ikd,
c     >               ir,
c     >               e2dgdown(ikd,ir),e2dgpara(ikd,ir)
c     >               *side_orthog,
c     >               side_orthog,e2dgpara(ikd,ir),
c     >               e2dflux(ik,ir),e2dnbs(ik,ir)*e2dvhs(ik,ir),
c     >               e2dpipara(ik,ir),e2dpepara(ik,ir)
c
c
            elseif (ikd.eq.nks(ir).and.korpg(ik,ir).ne.0) then
c
c              For the last listing on the ring - need to use
c              the data for the 3,4 cell face of the previous cell.
c
               in = korpg(ik,ir)
c
               sidelen = polysidelen(ik,ir,UP34,rc)
c
               sider = (rvertp(4,in) + rvertp(3,in)) /2.0
               sidez = (zvertp(4,in) + zvertp(3,in)) /2.0
c
               field_ratio = (kbfs(ik,ir) + kbfs(ik+1,ir)) /2.0
c
               side_orthog = 1.0
c
               theta1 = atan2c(zvertp(3,in)-zvertp(4,in),
     >                         rvertp(3,in)-rvertp(4,in))-PI/2.0
               theta2 = atan2c(zs(ikd,ir)-sidez,
     >                         rs(ikd,ir)-sider)
c
               side_orthog = cos(theta2-theta1)
c
               e2dgpara(ikd,ir) = e2dgdown(ikd,ir)
     >                 / (sidelen* 2.0 * PI * sider)
     >                 * field_ratio/side_orthog
c
c              Ion power flux
c
               e2dpipara(ikd,ir) = e2dpidown(ikd,ir)
     >                / (sidelen* 2.0 * PI * sider)
     >                * field_ratio/side_orthog
c
c              Electron power flux
c
               e2dpepara(ikd,ir) = e2dpedown(ikd,ir)
     >                / (sidelen* 2.0 * PI * sider)
     >                * field_ratio/side_orthog
c
c
c
c               write (6,'(a,3i4,1p,10(1x,g12.5))') 'G*DOWN',ik,ikd,
c     >               ir,
c     >               e2dgdown(ikd,ir),e2dgpara(ikd,ir)
c     >               *side_orthog,
c     >               side_orthog,e2dgpara(ikd,ir),
c     >               e2dflux(ik,ir),e2dnbs(ik,ir)*e2dvhs(ik,ir),
c     >               e2dpipara(ik,ir),e2dpepara(ik,ir)
c
            else
c
c               write(6,'(a,3i4,1p,1x,g12.5,a13,26x,4(1x,g12.5))')
c     >              'G-DOWN',ik,ikd,ir,
c     >             e2dgdown(ikd,ir),'  NO POLYGON ',e2dflux(ik,ir),
c     >             e2dnbs(ik,ir)*e2dvhs(ik,ir),
c     >               e2dpipara(ik,ir),e2dpepara(ik,ir)
c
c
            end if
c
         end do
c
      end do
c
      return
      end
c
c
c
      subroutine readauxarray(e2dfluxdata,
     >                        core_rings,sol_rings,pp_rings,
     >                        nrings,title)
      implicit none
c
      include 'params'
      include 'cgeom'
      include 'reader'
      include 'comtor'
c
      real e2dfluxdata(maxnks,maxnrs)
      integer core_rings,pp_rings,sol_rings,nrings
      character*(*) title
c
c     READAUXARRAY - This subroutine reads in an auxiliary
c                    information array from an extra Edge2d
c                    input file that specifies additional
c                    background information in a particular
c                    format.
c
c
c     Local Variables
c
      integer ikpp1,ikpp2,ikend
      integer nblocks,blkcnt,top_index,ik,ir,i,ib
      integer other_rings,first_sol
      integer blksize
      parameter (blksize=9)
      real dummy(maxnks,maxnrs)
c
      nblocks = (nrings-1)/blksize + 1
c
c     Initialize the core/pp separator values of ik
c
      ikpp1 = 0
      ikpp2 = 0
      ikend = 0
c
      write (6,*) 'Blocks:',nblocks,nrings

c
      do ib = 1,nblocks
c
         top_index = min(nrings-(ib-1)*blksize,blksize)
c
c        Read block separators
c
 555     read(12,'(a132)') buffer

         if (buffer(1:6).eq.'      '.or.
     >       buffer(1:6).eq.'------') goto 555

c
c        Start reading in blocks of data - process first line then
c        loop back.
c
c
c        Check for separator lines
c

 560     if (buffer(7:12).eq.'......') then
c
c           Set the location of the pp/core breaks in ik on the
c           first time through.
c
            if (ikpp1.eq.0) then
               ikpp1 = ik
            elseif (ikpp2.eq.0) then
               ikpp2 = ik
            endif
c
            goto 565
c
         elseif (buffer(1:6).eq.'      '.or.buffer(1:1).eq.'*') then
c
            if (ikend.eq.0) then
               ikend = ik
            endif
c
            goto 590
c
         endif

         read(buffer,*) ik,(dummy(ik,blksize*(ib-1)+i),i=1,top_index)

 565     read(12,'(a132)',end=600) buffer

c
         write (6,*) 'AUX:',buffer(1:80),ik
c
         goto 560

 590     continue

      end do

 600  continue

c
c     Assign the values in the dummy array to the DIVIMP array
c
c     The table contains a number of rings equal to nsol + max(npp,ncore)
c     It assumes that the sum of the Kpoints of each core+pp ring is equal
c     to the number of Kpoints in each SOL ring.
c
      other_rings = max(core_rings,pp_rings)
c
      first_sol = nrings - sol_rings + 1
c
c     Assign data to core and pp rings
c
c     Assign the main SOL first - assign ir index so that it maps onto
c     the actual SOL rings no matter whether there are more core or PP
c     rings in the actual print-out
c
      do ir = first_sol,nrings
         do ik = 1,ikend
            e2dfluxdata(ik,irsep+ir-first_sol) = dummy(ik,ir)
         end do
      end do
c
c     Assign core data
c
      do ir = other_rings-core_rings+1,first_sol-1
c
         do ik = ikpp1+1,ikpp2
c
            e2dfluxdata(ik-ikpp1,1+ir-(other_rings-core_rings+1)) =
     >                        dummy(ik,ir)
c
         end do
c
      end do

c
c     Assign private plasma data
c
      do ir = other_rings-pp_rings+1,first_sol-1
c
c        First part
c
         do ik = 1,ikpp1
c
            e2dfluxdata(ik,irtrap+ir-(other_rings-pp_rings+1)) =
     >                        dummy(ik,ir)
c
         end do
c
c        Second part
c
         do ik = ikpp2+1,ikend
c
            e2dfluxdata(ik-(ikpp2-ikpp1),
     >                  irtrap+ir-(other_rings-pp_rings+1)) =
     >                        dummy(ik,ir)
c
         end do
c
      end do
c
c     Print out array and continue
c
      if (cprint.eq.9) then
         CALL PRRMATDIV(e2dfluxdata,MAXNKS,nks(irsep),NRS,6,title)
      endif
c
      return
      end
c
c
c
      subroutine calc_mps
      implicit none
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'promptdep'
c
C-----------------------------------------------------------------------
c
c     CALC_MPS:
c
c     This routine calculates an estimate of the Magnetic Pre-sheath
c     thickness for each target segment. This information is used
c     in the prompt deposition calculations to define one of the
c     methods of prompt loss of ions to the targets.
c
c     David Elder,             Jan. 5 1998
c
C-----------------------------------------------------------------------
c
      integer ik,ir,id,it
      real    b_field
      real    larmor
      external larmor
c
      do id = 1,nds
c
         ik = ikds(id)
         ir = irds(id)
c
         if (ik.eq.1) then
            it = 2
         else
            it = 1
         endif
c
c        Extract magnetic field - use target based magnetic field
c        ratios but use first cell toroidal B-field for now since I don't
c        know if BTS is defined for the virtual cells.
c
c        jdemod - bug here calculating MPS data
c                 should just test kbfst <=1 
c
c         if (kbfst(ir,it).le.0.0.or.kbfst(ir,it).ge.1.0.or.
         if (kbfst(ir,it).le.1.0.or.
     >      id.eq.1.or.id.eq.ndsin.or.id.eq.ndsin+1.or.id.eq.nds) then
            b_field = 0.0
            mps_thickness(ir,it) = 0.0
            mps_energy(ir,it) = 0.0
         else
c
            b_field=kbfst(ir,it)/sqrt(kbfst(ir,it)**2-1.0)*bts(ik,ir)
c
            mps_thickness(ir,it) = larmor(crmb,ktids(id),b_field,rizb)
c
            mps_energy(ir,it) = abs(max(log(1.0/kbfst(ir,it)*
     >                 cos(thetas(id)-target_orth_angle(id))),-3.0))
c
         endif
c
         if (cprint.eq.3.or.cprint.eq.9) then
c
            write(6,'(a,i5,8(1x,g12.5))') 'MPS:',id,b_field,
     >                   mps_thickness(ir,it),mps_energy(ir,it),
     >                   thetas(id),target_orth_angle(id),
     >                   kbfst(ir,it),thetas(id)-target_orth_angle(id),
     >                   cos(thetas(id)-target_orth_angle(id))
c
         endif
c
      end do
c
      return
      end
c
c
c
      real function larmor(mz,Ez,B,Z)
      implicit none
c
      real mz,Ez,B,Z
c
c     LARMOR:
c
c     This function returns the value for the Larmor
c     radius in meters. (m)
c
c     Input:
c
c     mz = Mass           [amu]
c     Ez = Energy         [eV]
c     B  = Magnetic Field [Tesla]
c     Z  = Charge State
c
c     This routine returns a value of 0.0 for invalid input.
c
c
      if (Z.eq.0.or.B.eq.0.0.or.Ez.lt.0.0.or.mz.lt.0.0) then
         larmor = 0.0
      else
         larmor = 1.02e-4 * sqrt(mz*Ez) / (B * Z)
      endif
c
      return
c
      end
c
c
c
      subroutine find_midplane
      implicit none
      include 'params'
      include 'cgeom'
c
      integer ikmidin,ikmidout,ik,ir,ikmid,oumid,inmid,nr
      integer ikmidouts(maxnrs),ikmidins(maxnrs)
c
      real    dist0min,dist0,middistin,middistout
c
      do ir = 1,irwall
c
         ikmid = ikmids(ir)
C
C        Outer
C
         dist0min = hi
c
         DO IK=1,IKMID,1
           dist0 = (ZS(IK,IR)-Z0)**2
           IF (dist0.lt.dist0min) THEN
             dist0min = dist0
             IKMIDOUT=IK
           ENDIF
        end do
C
C       Inner
C
        dist0min = hi
c
        DO IK=NKS(IR),IKMID+1,-1
           dist0 = (ZS(IK,IR)-Z0)**2
           IF (dist0.lt.dist0min) THEN
             dist0min = dist0
             IKMIDIN=IK
           ENDIF
        end do
c
        if (ir.eq.irsep) then 
           ikmidout_sep = ikmidout
           ikmidin_sep = ikmidin 
        endif 
c
        rcouter(ir) = rs(ikmidout,ir)
        zcouter(ir) = zs(ikmidout,ir)
        rcinner(ir) = rs(ikmidin,ir)
        zcinner(ir) = zs(ikmidin,ir)
c
        ikmidouts(ir) = ikmidout
        ikmidins(ir)  = ikmidin
c
      end do
c
c     Construct mid-plane distances from the separatrix going inward and
c     outward for the inner and outer midplanes.
c
c     Inward from separatrix
c
      middistout = 0.0
      middistin = 0.0
      call rzero(separatrix_dist,maxnks*maxnrs)
c
      do ir = irsep-1,1,-1
c
         middistout = middistout - distout(ikmidouts(ir),ir)
         middistin  = middistin  - distout(ikmidins(ir),ir)
c
         middist(ir,1) = middistin
         middist(ir,2) = middistout
c
         middistout = middistout - distin(ikmidouts(ir),ir)
         middistin  = middistin  - distin(ikmidins(ir),ir)
c
         do ik = 1,nks(ir)-1
c
c           Calculate cell by cell approximate distances from the separatrix
c
            if (ir.eq.irsep-1) then

               separatrix_dist(ik,ir) = distout(ik,ir)

            else

               separatrix_dist(ik,ir) = tdistout(ik,ir) +
     >                                separatrix_dist(ikouts(ik,ir),ir)

            endif
c
         end do
c
      end do
c
c     Outward from separatrix
c
      middistout = 0.0
      middistin = 0.0
c
      do ir = irsep,irwall
c
         middistout = middistout + distin(ikmidouts(ir),ir)
         middistin  = middistin  + distin(ikmidins(ir),ir)
c
         middist(ir,1) = middistin
         middist(ir,2) = middistout
c
         middistout = middistout + distout(ikmidouts(ir),ir)
         middistin  = middistin  + distout(ikmidins(ir),ir)
c
      end do
c     
c     Assign midplane distances for the PFZ from equivalent
c     core flux surfaces for now
c
      do ir = nrs,irtrap,-1
c
         nr = irsep - (nrs-ir) -1 
c
         if (nr.gt.1) then 

            middist(ir,1) = middist(nr,1)
            middist(ir,2) = middist(nr,2)

         endif
c
      end do
c
c
c
      return
      end
c
c
c
      subroutine check_fluxes
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cedge2d'
c
c
      integer ik,ir,in,iz,rc,id,itarg
c
      real vel_h,vel_imp,side_len,bvel_h,bvel_i
      real totfluxes(10,5)
      real polysidelen
      real rtor
      external polysidelen
c
c     Calculate target ion fluxes
c
      call rzero(totfluxes,10*5)
c
      do id = 1,nds
c
         rtor = 2.0 * PI * rp(id)
c
         ik = ikds(id)
         ir = irds(id)
c
         itarg = 2
         if (id.le.ndsin) itarg = 1
c
         bvel_h = kvds(id)
         bvel_i = e2dvts(ir,itarg,0)
c
         totfluxes(1,1) = totfluxes(1,1) +
     >       knds(id) * abs(bvel_i) * dds2(id) * rtor
     >         /kbfs(ik,ir) * costet(id)
c
         totfluxes(1,5) = totfluxes(1,5) +
     >       knds(id) * abs(bvel_h) * dds2(id) * rtor
     >         /kbfs(ik,ir) * costet(id)
c
         vel_h = 1.56e4 * sqrt (ktibs(ik,ir)/crmb)
         vel_imp = 1.56e4 * sqrt (ktibs(ik,ir)/crmi)
c
         totfluxes(2,1) = totfluxes(2,1) +
     >     0.25*e2datom(ik,ir) * vel_h * dds2(id) * rtor
c
         totfluxes(3,1) = totfluxes(3,1) +
     >     0.25*e2dnzs(ik,ir,0) * vel_imp * dds2(id)*rtor
c
         write (6,'(a,4i4,1p,5(1x,g12.4))') 'V:',
     >      ik,ir,itarg,id,
     >      bvel_h,
     >      bvel_i,bvel_i/bvel_h,rp(id),dds2(id)
c
         do iz = 1,6

            totfluxes(3+iz,1) = totfluxes(3+iz,1) +
     >          e2dnzs(ik,ir,iz) * abs(bvel_i)
     >             * dds2(id) * rtor
     >         /kbfs(ik,ir) * costet(id)

            totfluxes(3+iz,5) = totfluxes(3+iz,5) +
     >          e2dnzs(ik,ir,iz) * abs(bvel_h)
     >             * dds2(id) * rtor
     >         /kbfs(ik,ir) * costet(id)

         end do

         write (6,'(a,2i4,2(1x,g12.5))') 'UEDGE:TARG OUT FLUXES:',
     >     ik,ir,
     >     0.25*e2datom(ik,ir) * vel_h * dds2(id) * rtor,
     >     0.25*e2dnzs(ik,ir,0) * vel_imp * dds2(id)*rtor
c
      end do
c
c     Calculate main wall neutral fluxes
c
      ir = irwall -1
c
      do ik = 1,nks(ir)
c
         vel_h = 1.56e4 * sqrt (ktibs(ik,ir)/crmb)
         vel_imp = 1.56e4 * sqrt (ktibs(ik,ir)/crmi)
c
         side_len = polysidelen(ik,ir,OUTWARD23,rc)
c
         in = korpg(ik,ir)
c
         if (in.eq.0.or.nvertp(in).ne.4) then
            rtor = 0.0
         else
            rtor = 2.0 * PI *
     >            (rvertp(2,in)+rvertp(3,in))/2.0
         endif
c
         totfluxes(2,2) = totfluxes(2,2) +
     >     0.25*e2datom(ik,ir)*vel_h * side_len * rtor
c
         totfluxes(3,2) = totfluxes(3,2) +
     >     0.25*e2dnzs(ik,ir,0)*vel_imp * side_len * rtor
c
         write (6,'(a,2i4,2(1x,g12.5))') 'UEDGE:WALL OUT FLUXES:',
     >     ik,ir,
     >     0.25*e2datom(ik,ir)*vel_h * side_len * rtor,
     >     0.25*e2dnzs(ik,ir,0)*vel_imp * side_len * rtor


      end do
c
c     Private plasma wall - neutral fluxes
c
      ir = irtrap + 1
c
      do ik = 1,nks(ir)
c
         vel_h = 1.56e4 * sqrt (ktibs(ik,ir)/crmb)
         vel_imp = 1.56e4 * sqrt (ktibs(ik,ir)/crmi)
c
         side_len = polysidelen(ik,ir,INWARD41,rc)
c
         in = korpg(ik,ir)
c
         if (in.eq.0.or.nvertp(in).ne.4) then
            rtor = 0.0
         else
            rtor = 2.0 * PI *
     >            (rvertp(4,in)+rvertp(1,in))/2.0
         endif
c
         totfluxes(2,3) = totfluxes(2,3) +
     >     0.25*e2datom(ik,ir) * vel_h * side_len*rtor
c
         totfluxes(3,3) = totfluxes(3,3) +
     >     0.25*e2dnzs(ik,ir,0) * vel_imp * side_len*rtor
c
      end do
c
c     Sum up
c
      do in = 1,9
         do id = 1,3
            totfluxes(in,4) = totfluxes(in,4)+totfluxes(in,id)
         end do
      end do
c
c     Print out
c
      write(6,300) 'Hydrogen Ions - Target   :',totfluxes(1,1),
     >         totfluxes(1,5)
      write(6,100) 'Hydrogen Neut - Target   :',totfluxes(2,1)
      write(6,100) 'Hydrogen Neut - Main Wall:',totfluxes(2,2)
      write(6,100) 'Hydrogen Neut - PP Wall  :',totfluxes(2,3)

      write(6,100) 'Target/Main Wall Yield   = ',const_yield
      write(6,100) 'PP Wall Yield            = ',
     >                            const_yield/10.0

      write(6,300) 'Carbon InFlux - Target     :',totfluxes(1,1)
     >       * const_yield,totfluxes(1,5) * const_yield
      write(6,100) 'Carbon InFlux - Target(?)  :',totfluxes(2,1)
     >                     * const_yield
      write(6,100) 'Carbon InFlux - Main Wall  :',totfluxes(2,2)
     >                     * const_yield
      write(6,100) 'Carbon InFlux - PP Wall    :',totfluxes(2,3)
     >                     * const_yield/10.0
      write(6,100) 'Total Carbon InFlux      :',
     >            totfluxes(1,1) * const_yield +
     >            totfluxes(2,1) * const_yield +
     >            totfluxes(2,2) * const_yield +
     >            totfluxes(2,3) * const_yield/10.0

      write(6,*) 'Carbon OutFlux:'

      do iz = 0,6

         write(6,400) 'Carbon OutFlux:', iz ,totfluxes(3+iz,4),
     >                totfluxes(3+iz,5)
         totfluxes(10,4) = totfluxes(10,4) + totfluxes(3+iz,4)

         if (iz.gt.0)
     >     totfluxes(10,3) = totfluxes(10,3) + totfluxes(3+iz,4)

      end do
c
      write(6,100) 'Carbon Neutral OutFlux - target   :',
     >                    totfluxes(3,1)
      write(6,100) 'Carbon Neutral OutFlux - main wall:',
     >                    totfluxes(3,2)
      write(6,100) 'Carbon Neutral OutFlux - pp wall  : ',
     >                    totfluxes(3,3)
c
      write(6,100) 'Total Carbon Neut OutFlux:',
     >     totfluxes(3,1)+totfluxes(3,2)+totfluxes(3,3)


      write(6,100) 'Total Carbon Ion OutFlux :',totfluxes(10,3)
      write(6,100) 'Total Carbon OutFlux     :',totfluxes(10,4)
c
      write(6,100) 'Total Excess/Deficit OF-IF:',
     >             totfluxes(10,4) -
     >            (totfluxes(1,1) * const_yield +
     >            totfluxes(2,1) * const_yield +
     >            totfluxes(2,2) * const_yield +
     >            totfluxes(2,3) * const_yield/10.0)
c
c     Formats
c
 100  format(a,5x,1p,g15.6)
 200  format(a,1x,i4,1p,g15.6)
 300  format(a,5x,1p,2(1x,g15.6))
 400  format(a,1x,i4,1p,2(1x,g15.6))

c
      return
      end
c
c
c
      subroutine setup_uedge_wall
      implicit none
      include 'params'
      include 'pindata'
      include 'cgeom'
      include 'comtor'
      include 'cedge2d'
c
c     SETUP_UEDGE_WALL: This routine will calculate and assign
c     wall fluxes of hydrogen for each of the segments of the wall.
c     This only works for wall option 7 which selects the outermost
c     polygon boundary for use as the neutral wall. The wall
c     fluxes contain only neutral hydrogen calculated from the
c     gas kinetic formula  F = nv/4 at the edge of the grid.
c
c     The energy loaded is simply the local ion temperature in the
c     corresponding grid cell.
c
c     The ion+atom flux array is loaded with the target ion flux plus
c     the target neutral flux since it is used in the flux dependent
c     portion of yield calculations.
c
c     David Elder                   1998, Oct 26
c
      integer ik,ir,id,in,icount,rc
      real neuth_vel,cs,ionflux
c
c     Loop through the wall
c
      do in = 1,wallpts
c
c        Construct flux and energy arrays for each segment
c
c        Set index for PIN equivalent wall data
c
         wallpt(in,17) = in
c
c        Main wall
c        PP Wall
c
         if (wallpt(in,16).eq.7.or.wallpt(in,16).eq.8) then
c
            if (wallpt(in,16).eq.7) then
               ir = irwall -1
            else
               ir = irtrap + 1
            endif
c
            ik = wallpt(in,24)
c
c           Only neutral fluxes for walls
c
            neuth_vel = 1.56e4 * sqrt (ktibs(ik,ir)/crmb)
c
            fluxhw(in) = 0.25 * e2datom(ik,ir) * neuth_vel
            flxhw2(in) = fluxhw(in)
            flxhw5(in) = ktibs(ik,ir)
c
c        Inner Target (target #1 - outer for Sonnet grids)
c        Outer Target (target #2 - inner for Sonnet grids)
c
         elseif (wallpt(in,16).eq.4.or.wallpt(in,16).eq.1) then
c
            id = wallpt(in,18)
            ik = ikds(id)
            ir = irds(id)
c
c           Neutral fluxes for walls + target ion fluxes
c
            neuth_vel = 1.56e4 * sqrt (ktids(id)/crmb)
c
            fluxhw(in) = 0.25*e2datom(ik,ir)*neuth_vel
c
            cs = abs(kvds(id))
c
            if (northopt.eq.0) then
c
               ionflux = KNDS(ID) * CS / KBFS(IK,IR)
c
            elseif (northopt.eq.1.or.northopt.eq.2
     >                .or.northopt.eq.3) then
c
               ionflux = KNDS(ID) * CS/KBFS(IK,IR)*COSTET(ID)
c
            endif

            flxhw2(in) = fluxhw(in) + ionflux
            flxhw5(in) = ktids(id)
c
            nimindex(id) = in
c
         endif
c
         if (cprint.eq.3.or.cprint.eq.9) then

            write (6,'(a,4i5,5(1x,g12.4))') 'Flux:',in,ik,
     >          ir,id,wallpt(in,16),wallpt(in,18),fluxhw(in),
     >          flxhw2(in), flxhw5(in)
         endif

c
      end do
c
      return
      end
c
c
c
      subroutine redef_pinwalldata
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'baffles'
c
c     REDEF_PINDATA: This routine redefines the NIMBUS wall and
c                    reorganizes and reassigns the data in the
c                    flux arrays to match the redefined walls.
c                    This includes spreading out the fluxes onto
c                    baffles over all of the appropriate surfaces.
c
c     This routine also performs sanity checks and issues error
c     messages indicating when and if problems occur and whether
c     the wall segments actually match up between the NIMBUS and
c     DIVIMP redefined walls.
c
c     The redefined DIVIMP wall is contained in the nves,rves,zves
c     arrays.
c
      integer in,iw,oin
c
      integer nvesmt,jvesmt(maxseg),tmpind,tmpseg
      real    rvest(maxseg,2),zvest(maxseg,2)
c
      integer maxbrokseg,curseg
      parameter (maxbrokseg=10)
      integer broksegcnt,broksegind(maxbrokseg,2)
      integer broksegused(maxbrokseg)
      real    brokseglen(maxbrokseg),totseglen
c
      integer bafcnt,bafseg
      integer bafind(mbufx+1,mbufle+2)
      real    baflen(mbufx+1,mbufle+1),segratio
c
c     Local wall redefintion
c
      integer nvesl,bafstart
      real rvesl(mves),zvesl(mves)
c
      write(6,*) 'REDEFINED NIMBUS WALL: PROCESSING'
c
c     Initialization
c
      broksegcnt = 0
      call izero(bafind,(mbufx+1)*(mbufle+1))
      call izero(broksegused,maxbrokseg)
      call rzero(baflen,(mbufx+1)*(mbufle+1))
c
c     Create local version of the NIMBUS wall
c
c     Find starting point of baffles in listing and only copy the wall up to
c     that point.
c
      bafstart = 0
c
      do in = 1,nvesm
c
c        Print out contents of flux data and vessel
c
c         write (6,'(a,2i4,8(1x,g12.5))') 'Flux:',in,jvesm(in),
c     >         rvesm(in,1),zvesm(in,1),fluxhw(in),flxhw2(in),
c     >         flxhw3(in),flxhw4(in),flxhw5(in),flxhw6(in)
c
c
         if ((jvesm(in).eq.9.or.jvesm(in).eq.10).and.bafstart.eq.0) then
c
            bafstart = in
c
         endif
c
      end do
c
c     Check bafstart - if it is zero then NIMBUS did NOT include
c     the baffle in the current run - simply exit in this case.
c     NIMBUS could do this if LBUFLE=F is specified in the
c     NIMBUS input. Returning should treat the results as if
c     no baffle was specified or included by NIMBUS.
c
      if (bafstart.eq.0) return

c
c------------------------------------------------------------------
c
c     CONTINUE PROCESSING IF BAFFLE FOUND
c
c
c     Copy wall up to the start of the NIMBUS baffles
c
      nvesl = bafstart -1
c
      do in = 1,nvesl

         rvesl(in) = rvesm(in,1)
         zvesl(in) = zvesm(in,1)
c
      end do
c
c     Call redefine_vessel to include the baffles into an extended
c     vessel wall.
c
      wallredef = 0
c
      call redefine_vessel(nvesl,
     >                     rvesl,zvesl,nbufle,
     >                     rbufle,zbufle,
     >                     nbufmx,nbufx,rbufx,zbufx,
     >                     node_origin,wallredef)
c
c     Copy redefined wall into new temporary version of the NIMBUS wall.
c
      nvesmt = nvesl
c
      do in = 1,nvesl

         if (in.eq.nvesl) then
            rvest(in,1) = rvesl(in)
            zvest(in,1) = zvesl(in)
            rvest(in,2) = rvesl(1)
            zvest(in,2) = zvesl(1)
         else
            rvest(in,1) = rvesl(in)
            zvest(in,1) = zvesl(in)
            rvest(in,2) = rvesl(in+1)
            zvest(in,2) = zvesl(in+1)
         endif
c
      end do
c
c     Construct jvesmt array values based on the original
c     source wall segments.
c
      do in = 1,nvesmt
c
         tmpind = node_origin(in,1)
         tmpseg = node_origin(in,2)

c
         if (tmpind.eq.1) then
c
c           Main wall segment
c
c
c           Assign JVESM - check for possible errors - mismatches
c           should only occur for wall segments that are broken to
c           connect the baffle.
c
            jvesmt(in) = jvesm(tmpseg)
c
            if (rvest(in,1).ne.rvesm(tmpseg,1).or.
     >          rvest(in,2).ne.rvesm(tmpseg,2).or.
     >          zvest(in,1).ne.zvesm(tmpseg,1).or.
     >          zvest(in,2).ne.zvesm(tmpseg,2)) then
c
                write (6,'(a,2i4,8(1x,g12.5))') 'Segment mismatch:',
     >                in,tmpseg,
     >                rvest(in,1),zvest(in,1),
     >                rvesm(tmpseg,1),zvesm(tmpseg,1),
     >                rvest(in,2),zvest(in,2),
     >                rvesm(tmpseg,2),zvesm(tmpseg,2)
c
c               Collect information on broken wall segments
c
                if (broksegcnt.lt.maxbrokseg)  then
c
                   broksegcnt = broksegcnt + 1
                   broksegind(broksegcnt,1) = in
                   broksegind(broksegcnt,2) = tmpseg
                   brokseglen(broksegcnt) =
     >                       sqrt( (rvest(in,1)-rvest(in,2))**2
     >                            +(zvest(in,1)-zvest(in,2))**2)
c
                endif
c
            endif
c
         else
c
c           Baffle Segment
c
c           Assign a value of 10 to a baffle segment - accumulate
c           information on baffle segment index and length
c
            jvesmt(in) = 10
c
c           Record wall index cross-reference for baffle segment
c
            bafind(tmpind,tmpseg)=in
            bafind(tmpind,mbufle+1) = bafind(tmpind,mbufle+1) +1
c
c           Record lengths and cumulative length of baffle segments
c           for later distribution of NIMBUS fluxes.
c
            baflen(tmpind,tmpseg) =
     >                       sqrt( (rvest(in,1)-rvest(in,2))**2
     >                            +(zvest(in,1)-zvest(in,2))**2)
            baflen(tmpind,mbufle+1) = baflen(tmpind,mbufle+1) +
     >                                baflen(tmpind,tmpseg)
c
         endif

c
c         write(6,'(a,4i4,8(1x,g12.5))') 'DATA:',in,tmpind,tmpseg,
c     >         jvesmt(in),
c     >         rvest(in,1),zvest(in,1),rvesm(tmpseg,1),zvesm(tmpseg,1),
c     >         rvest(in,2),zvest(in,2),rvesm(tmpseg,2),zvesm(tmpseg,2)
c
      end do
c
c     Find indexes for baffles in original NIMBUS wall FLUX data - these
c     are different from the indices in the array showing the
c     location and coordinates of the baffles.
c
      bafcnt = 1
c
      do in = 1,nvesm
c
c        Compound baffle segment - should be only one for each baffle
c
c        Baffles should be listed in correct order
c
c        Note the indexing goes from 2+ for baffles since the code
c        that redefines the wall uses the actual vessel wall as baffle
c        one.
c
         if (jvesm(in).eq.9) then
c
            bafcnt = bafcnt +1
            bafind(bafcnt,mbufle+2) = bafstart + bafcnt - 2
c
c            write(6,'(a,4i4,4(1x,g12.5))') 'BAFFLES::',in,jvesm(in),
c     >         bafcnt,bafind(bafcnt,mbufle+2)
c
         endif
c
      end do
c
c     Start reassigning the flux arrays to match the new wall definition
c     This assumes that they actually match up - earlier diagnostics will c
c     show whether this is true or not - process baffle segments differently.
c
c     Take data from original PIN flux arrays
c
      do in = 1,nvesmt
c
         tmpind = node_origin(in,1)
         tmpseg = node_origin(in,2)
c
         if (jvesmt(in).eq.10) then
c
c           Treat baffle segments differently
c
            bafseg = bafind(tmpind,mbufle+2)
            segratio = baflen(tmpind,tmpseg)/baflen(tmpind,mbufle+1)
c
c            write(6,'(a,4i4,8(1x,g12.5))') 'PIN:BAFF:',in,jvesmt(in),
c     >         tmpind,bafseg,segratio,
c     >         baflen(tmpind,tmpseg),baflen(tmpind,mbufle+1),
c     >         fluxhw_pin(bafseg),flxhw2_pin(bafseg),flxhw3_pin(bafseg),
c     >         flxhw5_pin(bafseg)
c
c
            fluxhw(in) = fluxhw_pin(bafseg) * segratio
            flxhw2(in) = flxhw2_pin(bafseg) * segratio
            flxhw3(in) = flxhw3_pin(bafseg) * segratio
            flxhw4(in) = flxhw4_pin(bafseg) * segratio
            flxhw6(in) = flxhw6_pin(bafseg) * segratio
c
c           This entry is not split because it is an average
c           impact energy and not a flux.
c
            flxhw5(in) = flxhw5_pin(bafseg)
c
         else
c
c           Wall segments
c
c            write(6,'(a,4i4,4(1x,g12.5))') 'PIN:WALL:',in,jvesmt(in),
c     >         tmpind,tmpseg,
c     >         fluxhw_pin(tmpseg),flxhw2_pin(tmpseg),flxhw3_pin(tmpseg),
c     >         flxhw5_pin(tmpseg)
c
            fluxhw(in) = fluxhw_pin(tmpseg)
            flxhw2(in) = flxhw2_pin(tmpseg)
            flxhw3(in) = flxhw3_pin(tmpseg)
            flxhw4(in) = flxhw4_pin(tmpseg)
            flxhw5(in) = flxhw5_pin(tmpseg)
            flxhw6(in) = flxhw6_pin(tmpseg)

         endif

      end do
c
c     Deal with broken segments - adjust the fluxes onto these segments
c
c     There should only be two broken segments - they may be part of the
c     same vessel wall segment.
c
      do in = 1,broksegcnt
c
         if (broksegused(in).ne.1) then
c
            curseg = broksegind(in,2)
            totseglen = brokseglen(in)
c
            do iw = in+1,broksegcnt
c
               if (curseg.eq.broksegind(iw,2)) then
c
c                 Add other segments to the total for the
c                 current segment.
c
                  totseglen = totseglen + brokseglen(iw)
c
               endif
c
            end do
c
            segratio = brokseglen(in) / totseglen
c
c            write(6,'(a,3i4,4(1x,g12.5))') 'BROK1:',in,
c     >         curseg,broksegcnt,segratio,brokseglen(in),
c     >         totseglen
c
c
            if (segratio.lt.1.0) then
c
c              Adjust all but flxhw5 which is not a flux
c
               oin = broksegind(in,1)
c
               fluxhw(oin) = fluxhw_pin(curseg) * segratio
               flxhw2(oin) = flxhw2_pin(curseg) * segratio
               flxhw3(oin) = flxhw3_pin(curseg) * segratio
               flxhw4(oin) = flxhw4_pin(curseg) * segratio
               flxhw6(oin) = flxhw6_pin(curseg) * segratio
c
c
c            write(6,'(a,4i4,6(1x,g12.5))') 'BROK2:',in,oin,
c     >         curseg,broksegcnt,segratio,brokseglen(in),
c     >         totseglen,fluxhw(oin),flxhw2(oin),flxhw5(oin)
c
c
c              Loop adjusting the other segments on this one
c              segment
c
               do iw = in+1,broksegcnt
c
                  if (curseg.eq.broksegind(iw,2)) then
c
c                    Add other segments to the total for the
c                    current segment.
c
                     segratio = brokseglen(iw) / totseglen
c
c                    Adjust all but flxhw5 which is not a flux
c
                     oin = broksegind(iw,1)
c
                     fluxhw(oin) = fluxhw_pin(curseg) * segratio
                     flxhw2(oin) = flxhw2_pin(curseg) * segratio
                     flxhw3(oin) = flxhw3_pin(curseg) * segratio
                     flxhw4(oin) = flxhw4_pin(curseg) * segratio
                     flxhw6(oin) = flxhw6_pin(curseg) * segratio

c
c            write(6,'(a,3i4,6(1x,g12.5))') 'BROK3:',oin,
c     >         curseg,broksegcnt,segratio,brokseglen(iw),
c     >         totseglen,fluxhw(oin),flxhw2(oin),flxhw5(oin)
c
c
                     broksegused(iw) =1
c
                  endif
c
               end do
c
            endif
c
            broksegused(in) = 1

         endif
c
      end do
c
c     Copy over the nves,rves,zves and jves information to complete
c     the redefinition
c
      nvesm = nvesmt
c
      do in = 1, nvesmt
c
         jvesm(in) = jvesmt(in)
c
         do iw = 1,2
c
            rvesm(in,iw) = rvest(in,iw)
            zvesm(in,iw) = zvest(in,iw)
c
         end do
c
      end do
c
c     Copy over the pump information to the end of the modified arrays
c
c     At present no fluxes for the pump elements are passed back -
c     this could be changed in future and modifications would be
c     required for the code below.
c
c
      do in = 1, nvesp
c
         jvesm(nvesm+in) = jvesm_pin(nvesm_pin+in)
c
         do iw = 1,2
c
            rvesm(nvesm+in,iw) = rvesm_pin(nvesm_pin+in,iw)
            zvesm(nvesm+in,iw) = zvesm_pin(nvesm_pin+in,iw)
c
         end do
c
      end do
c
c     Finished wall redefinition ... dealing with NIMBUS wall
c
      return
      end
c
c
c
      subroutine grid_check
      use error_handling
      implicit none
      include 'params'
      include 'cgeom'
c
c
c     GRID_CHECK: Loop through all of the grid polygons and make sure
c                 that adjacent polygons share the appropriate sides.
c
      integer ik,ir,ikn,irn,nv, in, testin,err,id,cnt
      character*200 error_comment
c
c     Initialize error counter
c
      err = 0
      cnt = 0
c
      do ir = 1,nrs
c
c         if (ir.ne.1.and.ir.ne.irwall.and.ir.ne.irtrap) then
c
         write(6,'(a,2i10)') 'GRID_CHECK: Checking ring:',ir,nks(ir)
c
            do ik = 1,nks(ir)
c
               in = korpg(ik,ir)
c
c              IPP/08 Krieger - ensure index of nvertp is not zero
               if (in.ne.0.and.nvertp(max(1,in)).ne.0) then
c
                  cnt = cnt + 1
c
c                 Check all sides of polygon for adjacent matches
c
c                 Along field line - forward or up
c
                  if (ik.ne.nks(ir)) then
c
                     testin = korpg(ik+1,ir)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0)
     >                           then

                        if (rvertp(3,in).ne.rvertp(2,testin).or.
     >                      zvertp(3,in).ne.zvertp(2,testin).or.
     >                      rvertp(4,in).ne.rvertp(1,testin).or.
     >                      zvertp(4,in).ne.zvertp(1,testin)) then
c
                            err = err + 1
c
                            write(error_comment,'(a,6i5)') 
     >                          'CELL GEOMETRY ERROR: SIDE UP     34:'
     >                                //'(IK1,IR1,IK2,IR2,IN1,IN2):',
     >                                ik,ir,ik+1,ir,
     >                                in,testin
                            
                            call errmsg('ERROR:TAU MODULE:'//
     >                                  'ROUTINE GRID_CHECK',
     >                         error_comment(1:len_trim(error_comment)))

c
                  write(6,'(i6,8(1x,g12.5))') in,((rvertp(id,in),
     >                           zvertp(id,in)),
     >                           id = 1,4)
                  write(6,'(i6,8(1x,g12.5))')testin,((rvertp(id,testin),
     >                           zvertp(id,testin)),
     >                           id = 1,4)
c
                        endif
c
                     endif
c
                  endif 
c
c
c                 Along field line - backward or down
c
                  if (ik.ne.1) then
c
                     testin = korpg(ik-1,ir)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0)then

                        if (rvertp(1,in).ne.rvertp(4,testin).or.
     >                      zvertp(1,in).ne.zvertp(4,testin).or.
     >                      rvertp(2,in).ne.rvertp(3,testin).or.
     >                      zvertp(2,in).ne.zvertp(3,testin)) then
c
                            err = err + 1
c
                            write(error_comment,'(a,6i5)') 
     >                           'CELL GEOMETRY ERROR: SIDE DOWN   12:'
     >                                //'(IK1,IR1,IK2,IR2,IN1,IN2):',                            
     >                                ik,ir,ik-1,ir,
     >                                in,testin
                            
                            call errmsg('ERROR:TAU MODULE:'//
     >                                  'ROUTINE GRID_CHECK',
     >                         error_comment(1:len_trim(error_comment)))

                  write(6,'(i6,8(1x,g12.5))') in,((rvertp(id,in),
     >                           zvertp(id,in)),
     >                           id = 1,4)
                  write(6,'(i6,8(1x,g12.5))')testin,((rvertp(id,testin),
     >                           zvertp(id,testin)),
     >                           id = 1,4)

c
                        endif
c
                     endif
c
                  endif
c
c
c                 Cross field line - outward
c
c                  if (ir.ne.irwall-1) then
                  if (ir.ne.irwall) then
c
                     ikn = ikouts(ik,ir)
                     irn = irouts(ik,ir)

                     testin = korpg(ikn,irn)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0.and.
     >                   ik.ne.ikn.and.ir.ne.irn ) then
c
                        if (rvertp(2,in).ne.rvertp(1,testin).or.
     >                      zvertp(2,in).ne.zvertp(1,testin).or.
     >                      rvertp(3,in).ne.rvertp(4,testin).or.
     >                      zvertp(3,in).ne.zvertp(4,testin)) then
c
                            err = err + 1
c
                            write(error_comment,'(a,6i5)') 
     >                            'CELL GEOMETRY ERROR: SIDE OUTWARD23:'
     >                                //'(IK1,IR1,IK2,IR2,IN1,IN2):',
     >                                ik,ir,ikn,irn,
     >                                in,testin
                            
                            call errmsg('ERROR:TAU MODULE:'//
     >                                  'ROUTINE GRID_CHECK',
     >                         error_comment(1:len_trim(error_comment)))

                  write(6,'(i6,8(1x,g12.5))') in,((rvertp(id,in),
     >                           zvertp(id,in)),
     >                           id = 1,4)
                  write(6,'(i6,8(1x,g12.5))')testin,((rvertp(id,testin),
     >                           zvertp(id,testin)),
     >                           id = 1,4)
c
                        endif
c
                     endif
c
                  endif
c
c
c                 Cross field line - inward
c
c                  if (ir.ne.2.and.ir.ne.irtrap+1) then
                  if (ir.ne.1.and.ir.ne.irtrap) then
c
                     ikn = ikins(ik,ir)
                     irn = irins(ik,ir)

                     testin = korpg(ikn,irn)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0.and.
     >                   ik.ne.ikn.and.ir.ne.irn ) then
c
                        if (rvertp(4,in).ne.rvertp(3,testin).or.
     >                      zvertp(4,in).ne.zvertp(3,testin).or.
     >                      rvertp(1,in).ne.rvertp(2,testin).or.
     >                      zvertp(1,in).ne.zvertp(2,testin)) then
c
                            err = err + 1
c
                            write(error_comment,'(a,6i5)') 
     >                            'CELL GEOMETRY ERROR: SIDE INWARD 41:'
     >                                //'(IK1,IR1,IK2,IR2,IN1,IN2):',
     >                                ik,ir,ikn,irn,
     >                                in,testin
                            
                            call errmsg('ERROR:TAU MODULE:'//
     >                                  'ROUTINE GRID_CHECK',
     >                         error_comment(1:len_trim(error_comment)))
c
              write(6,'(i6,8(1x,g12.5))') in,((rvertp(id,in),
     >                         zvertp(id,in)),
     >                         id = 1,4)
              write(6,'(i6,8(1x,g12.5))') testin,((rvertp(id,testin),
     >                         zvertp(id,testin)),
     >                         id = 1,4)
c
                        endif
c
                     endif
c
                  endif
c
               endif
c
            end do
c
c         endif
c
      end do
c
      write (6,*) 'POLYGON GEOMETRY ERRORS:', err, ' OF ', cnt
c
      return
c
      end
c
c
c
      subroutine calc_s_reflect
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c
c     CALC_S_REFLECT
c
c     This routine uses parallel GEOMETRY ERRORS (non adjacent parallel 
c     cell boudaries on rings that should be adjacent) to find the S-locations
c     of the breaks in the rings where parallel ion reflections should occur. 
c     The array is initially zeroed which will turn off parallel ion reflection
c     even if the option is activated so it should not have an effect on a
c     contiguous grid without parallel geometry errors. 
c
      integer ik,ir,in,testin,id
c
c     Zero array
c
      call rzero(s_reflect,maxnrs*2)  
c
      do ir = 1,nrs
c
         if (ir.ne.1.and.ir.ne.irwall.and.ir.ne.irtrap) then
c
            do ik = 1,nks(ir)
c
               in = korpg(ik,ir)
c
c              IPP/08 Krieger - ensure index of nvertp is not zero
               if (in.ne.0.and.nvertp(max(1,in)).ne.0) then
c
c                 Check parallel sides of polygon for adjacent matches
c
c                 Along field line - forward
c
                  if (ik.ne.nks(ir)) then
c
                     testin = korpg(ik+1,ir)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0)then

                        if (rvertp(3,in).ne.rvertp(2,testin).or.
     >                      zvertp(3,in).ne.zvertp(2,testin).or.
     >                      rvertp(4,in).ne.rvertp(1,testin).or.
     >                      zvertp(4,in).ne.zvertp(1,testin)) then
c
c                           Forward geometry error found - assume this is 
c                           upper limit in S for first half ring.  
c
                            s_reflect(ir,1) = kss(ik,ir) 
c
                            write(6,'(a,6i5,g12.5)') 
     >                                'S_REFLECT FORWARD:',
     >                                ik,ir,ik+1,ir,
     >                                in,testin,s_reflect(ir,1)
c
                        endif
c
                     endif
c
                  endif
c
c                 Along field line - backward
c
                  if (ik.ne.1) then
c
                     testin = korpg(ik-1,ir)
c
c                    IPP/08 Krieger - ensure index of nvertp is not zero
                     if (testin.ne.0.and.nvertp(max(1,testin)).ne.0)then

                        if (rvertp(1,in).ne.rvertp(4,testin).or.
     >                      zvertp(1,in).ne.zvertp(4,testin).or.
     >                      rvertp(2,in).ne.rvertp(3,testin).or.
     >                      zvertp(2,in).ne.zvertp(3,testin)) then
c
c                           Backward geometry error found - assume this is 
c                           lower limit in S for second half ring.  
c
                            s_reflect(ir,2) = kss(ik,ir) 
c
                            write(6,'(a,6i5,g12.5)') 
     >                                'S_REFLECT BACK   :',
     >                                ik,ir,ik-1,ir,
     >                                in,testin,s_reflect(ir,2)
c
c
                        endif
c
                     endif
c
                  endif
c
               endif
c
            end do
c
         endif
c
      end do
c
      if (cprint.eq.3.or.cprint.eq.9) then 
c
         write(6,'(a,i5)') 'PARALLEL ION REFLECTION OPT: ',s_reflect_opt
c 
         write(6,'(a)') 'S-VALUES FOR S REFLECTION:'
c
         do ir = 1,nrs
c                 
            write(6,'(a,i5,2(1x,f10.3))') 'RING:',ir,
     >                 s_reflect(ir,1),s_reflect(ir,2)
c
         end do
c
      end if  
c
c
c
      return
c
      end
c
c
c
c     Hybrid wall processing routine
c
      subroutine rdhybd(ihybrid,nves,rvmod,zvmod,rves,zves)
      implicit none
c
      integer ihybrid, nves
c      real*8  rvmod(nves), zvmod(nves), rves(nves), zves(nves)
      real*4  rvmod(nves), zvmod(nves), rves(nves), zves(nves)
c
      integer lhybrd, igeom, ngeom, ncoord, i, j, k
c      real*8  coord(4)
      real*4  coord(4)
      character buffer*80, string*80
c
      real thresh
      parameter(thresh=0.0002)
c
      data lhybrd/28/
c
      rewind(lhybrd)
c
c  read number of hybrid geometries available in input file
c
    5 read(lhybrd,'(a80)') buffer
      if (buffer(1:1).eq.'C' .or. buffer(1:1).eq.'c') goto 5
      read(buffer,*) string, ngeom
c
c  check that ihybrid has a valid value
c
      if (ihybrid.lt.1 .or. ihybrid.gt.ngeom) then
        write(6,*) ' INVALID VALUE FOR IHYBRID:  ', ihybrid
        stop 'INVALID HYBRID WALL SELECTED'
      endif
c
      igeom = 1
c
c  read number of coordinate pairs for next geometry
c
   10 read(lhybrd,'(a80)') buffer
      if (buffer(1:1).eq.'C' .or. buffer(1:1).eq.'c') goto 10
      read(buffer,*) string, ncoord
c
c  check number of coordinates against nves.  There should be
c  one less since the unmodified vessel was closed by setting
c  the last (extra) point equal to the first.
c
c      if (igeom.eq.ihybrid .and. ncoord.ne.nves-1) then
c
c     IN NIMBUS the wall has been closed by setting the first
c     point equal to the last one - this is not in the grid
c     file.
c
      if (igeom.eq.ihybrid .and. ncoord.ne.nves) then
        write(6,*) ' INCONSISTENT VALUES FOR NCOORD AND NVES: ',
     >                ncoord, nves
        stop 'HYBRID WALL DOES NOT MATCH ACTUAL WALL'
      endif
c
c  skip title line
c
      read(lhybrd,'(a80)') buffer
c
c  read the coordinates one line at a time
c
      do i = 1, ncoord
        read(lhybrd,*) j, (coord(k),k=1,4)
c
c  if this is the required geometry, check that unmodified
c  values are OK and load modifications
c
        if (igeom.eq.ihybrid) then
          do k=1,4
            coord(k) = coord(k)/1000.0
          enddo

          write(6,*) 'HYB:',i,rves(i),coord(1),zves(i),coord(2)

          if ((abs(coord(1)-rves(i)).gt.thresh)
     >      .and.(coord(3).lt.1.0e20.or.coord(4).lt.1.0e20)) then
            write(6,*) ' INCONSISTENT R-COORDINATES: ',
     >                 i,coord(1),rves(i)
            stop 'HYBRID: INCONSISTENT R-COORDINATES'
          endif
          if ((abs(coord(2)-zves(i)).gt.thresh)
     >      .and.(coord(3).lt.1.0e20.or.coord(4).lt.1.0e20)) then
          write(6,*) ' INCONSISTENT Z-COORDINATES: ',
     >                 i,coord(2),zves(i)
            stop 'HYBRID: INCONSISTENT Z-COORDINATES'
          endif
          rvmod(i) = coord(3)
          zvmod(i) = coord(4)
        endif
      enddo
c
c  If we haven't got to the required geometry, loop back
c
      if (ihybrid.gt.igeom) then
        igeom = igeom + 1
        goto 10
      endif
c
      return
      end
c
c    POLCHG from -
c*** /u/sim/edge2d/source/defaultn/pf2ds/data.f ***
c
C
C=======================================================================
      SUBROUTINE POLCHG( IOUT , LOOP , NPOL , RMOD , ZMOD
     &                 , RPOL , ZPOL )
      IMPLICIT NONE
C
C+ .....................................................................
C
C ROUTINE : CHANGE POLYGON COORDINATES
C           --  -  ---
C VERSION : V1.R1.M0
C
C PURPOSE : TO OVERWRITE THE COORDINATES OF A POLYGON IS THE MODIFIER
C           IS LESS THAN 1.0E+20
C
C INPUT   : (I*4) IOUT         > 0 --- PRINT TABLE TO CHANNEL 'IOUT'
C           (LOG) LOOP         = .T. --- LET LAST POINT EQUAL THE 1ST
C           (I*4) NPOL         = NO. OF COORDINATES
C           (R*4) RMOD()       = R-COORDINATE OF NEW POLYGON IF <1.0E+20
C           (R*4) ZMOD()       = Z-COORDINATE OF NEW POLYGON IF <1.0E+20
C
C I/O     : (R*4) RPOL()       = OLD R-COORDINATES (ON INPUT)
C                              = NEW R-COORDINATES (ON OUTPUT)
C           (R*4) ZPOL()       = OLD Z-COORDINATES (ON INPUT)
C                              = NEW Z-COORDINATES (ON OUTPUT)
C
C AUTHOR  : JAMES SPENCE  (K1/0/80)  EXT. 4865
C           JET
C
C HISTORY : V1.R1.M0 --- 30/04/96 --- CREATION
C
C- .....................................................................
C
C..INPUT
      INTEGER*4  IOUT       , NPOL
c      REAL*8     RMOD(NPOL) , ZMOD(NPOL)
      REAL*4     RMOD(NPOL) , ZMOD(NPOL)
      LOGICAL    LOOP
C
C..I/O
c      REAL*8     RPOL(NPOL) , ZPOL(NPOL)
      REAL*4     RPOL(NPOL) , ZPOL(NPOL)
C
C..PROGRAM
      INTEGER*4  I          , NMAX
c      REAL*8     RTRAN      , ZTRAN
      REAL*4     RTRAN      , ZTRAN
C
C-----------------------------------------------------------------------
C
      IF( LOOP ) THEN
          NMAX       = NPOL - 1
      ELSE
          NMAX       = NPOL
      END IF
C
C-----------------------------------------------------------------------
C
      DO 100 I       = 1 , NMAX
         RTRAN       = RPOL(I)
         ZTRAN       = ZPOL(I)
         IF( RMOD(I).LT.1.00E+20 .AND. ZMOD(I).LT.1.00E+20 ) THEN
             RPOL(I) = RMOD(I)
             ZPOL(I) = ZMOD(I)
             IF( IOUT.GT.0 )
     >           WRITE(IOUT,1000) I,RTRAN,ZTRAN,RMOD(I),ZMOD(I)
         ELSE IF( RMOD(I).LT.1.00E+20 .AND. ZMOD(I).GE.1.00E+20 ) THEN
             RPOL(I) = RMOD(I)
             IF( IOUT.GT.0 )
     >           WRITE(IOUT,1010) I , RTRAN , ZTRAN , RMOD(I)
         ELSE IF( RMOD(I).GE.1.00E+20 .AND. ZMOD(I).LT.1.00E+20 ) THEN
             ZPOL(I) = ZMOD(I)
             IF( IOUT.GT.0 )
     >           WRITE(IOUT,1020) I , RTRAN , ZTRAN , ZMOD(I)
         ELSE
             IF( IOUT.GT.0 )
     >           WRITE(IOUT,1030) I , RTRAN , ZTRAN
         END IF
  100 CONTINUE
C
C-----------------------------------------------------------------------
C
      IF( LOOP ) THEN
          RPOL(NPOL)     = RPOL(1)
          ZPOL(NPOL)     = ZPOL(1)
      END IF
C
C-----------------------------------------------------------------------
C
 1000 FORMAT( I6 , 1P , 2E12.4 , E12.4        , E12.4        )
 1010 FORMAT( I6 , 1P , 2E12.4 , E12.4        , 1X , 11('.') )
 1020 FORMAT( I6 , 1P , 2E12.4 , 1X , 11('.') , E12.4        )
 1030 FORMAT( I6 , 1P , 2E12.4 , 1X , 11('.') , 1X , 11('.') )
C
C-----------------------------------------------------------------------
C
      RETURN
      END
c
c
c
      subroutine calc_asep_eff
      implicit none
      include 'params'
      include 'cgeom'
c
c     CALC_ASEP_EFF: This routine calculates a number of 
c                    area quantities related to the 
c                    separatrix for use in estimating
c                    Dperp and Xperp
c      
      integer ik,ir,ikmidref,ikout,irout,rc
      real polysidelen,sidelen,dist,get_refdist,refdist
      external polysidelen,get_refdist
c
c     Setup 
c
c     Assign irout equal to irsep - much more straightforward
c
      ir = irsep -1
      irout = irsep 
c
      asep = 0.0
      asep_eff = 0.0
c
      asep_tor = 0.0
      asep_tor_eff = 0.0 
c
c     Find the IK Index of the Outer midplane and use this as the 
c     flux expansion reference point.  
c
      if (xpoint_up) then       
C
         ikmidref = ikmidout_sep
c
      else
C
         ikmidref = ikmidin_sep
c
      endif 
c
      refdist = get_refdist(ikmidref,irout,INWARD41,rc) 
c
c     Loop around last core ring - not including last element
c
      do ik =  1,nks(ir)-1
c
         sidelen = polysidelen(ik,ir,OUTWARD23,rc)
c
         ikout = ikouts(ik,ir)
         irout = irouts(ik,ir)
c
         dist = get_refdist(ikout,irout,INWARD41,rc)
c
c        Note - toroidal areas are using cell center R values
c        as opposed to the mid-points of the polygon sides. This is 
c        a small difference and can be corrected for by using 
c        exact values later if necessary. 
c
         asep = asep + sidelen
         asep_tor = asep_tor + sidelen * 2.0 * PI * rs(ikout,irout)  
c                 
         if (dist.ne.0.0) then  
c
            asep_eff = asep_eff + sidelen * refdist / dist
c
            asep_tor_eff = asep_tor_eff 
     >                   + sidelen * refdist / dist 
     >                     * 2.0 * PI * rs(ikout,irout)
c
         endif 
c
      end do
c
c
      return
      end
c
c
c
      real function get_refdist(ik,ir,side,rc)
      implicit none
      integer ik,ir,side,rc
c
c     Commons
c
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     GET_REFDIST:
c
c     This routine returns the distance from the cell center to
c     the midpoint of the specified side. 
c
c     If a problem is encountered then rc is set to 1.
c
c     Normal rc is 0
c
c
      integer in,npoly,nsides,ik2,ir2
      real rp1e,zp1e,rp2e,zp2e,rs1e,zs1e,rs2e,zs2e,rie,zie
c
c     Set default rc
c
c slmod begin - temp
c      WRITE(0,*) '-->',ik,ir,side,rc
c slmod end
      rc = 0
c
c     Does a polygon exist for the primary cell - if not return a
c     length extracted from the cell centre information.
c
      npoly = korpg(ik,ir)
      nsides = nvertp(npoly)
      in = side
c
      if (npoly.ne.0.and.nsides.ge.side) then
c
c        Polygon exists and has sufficient sides
c
c        Look only at side passed in on call.
c
c
         if (side.eq.INWARD41) then 
c
            rs1e=rvertp(4,npoly)
            zs1e=zvertp(4,npoly)
            rs2e=rvertp(1,npoly)
            zs2e=zvertp(1,npoly)
c
         else
c
            rs1e=rvertp(in,npoly)
            zs1e=zvertp(in,npoly)
            rs2e=rvertp(in+1,npoly)
            zs2e=zvertp(in+1,npoly)
c
         endif
c
c         if (side.eq.1) then
c
c            rs1e=rvertp(nsides,npoly)
c            zs1e=zvertp(nsides,npoly)
c            rs2e=rvertp(in,npoly)
c            zs2e=zvertp(in,npoly)
c
c         else
c
c            rs1e=rvertp(in-1,npoly)
c            zs1e=zvertp(in-1,npoly)
c            rs2e=rvertp(in,npoly)
c            zs2e=zvertp(in,npoly)
c
c         endif
c
         rie = (rs1e+rs2e)/2.0
         zie = (zs1e+zs2e)/2.0
c
c        Calculate distance ...
c
         get_refdist = sqrt((rs(ik,ir)-rie)**2+(zs(ik,ir)-zie)**2)
         rc = 0
c
      else
c
c        Invalid polygon - calculate a length using alternate methods
c
         if (side.eq.DOWN12) then 
c
            ik2 = max(ik-1,1)
            ir2 = ir   
c
         elseif (side.eq.UP34) then 
c
            ik2 = min(ik+1,nks(ir))
            ir2 = ir   
c
         elseif (side.eq.OUTWARD23) then 
c
            ik2 = ikouts(ik,ir)
            ir2 = irouts(ik,ir)
c
         elseif (side.eq.INWARD41) then 
c
            ik2 = ikins(ik,ir)
            ir2 = irins(ik,ir)
c
         endif          

         get_refdist = sqrt((rs(ik,ir)-rs(ik2,ir2))**2 
     >                    + (zs(ik,ir)-zs(ik2,ir2))**2)/2.0
         rc = 1
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
      subroutine calc_targfluxdata
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'printopt'
c slmod begin
      include 'slcom'
c slmod end
c
c     CALC_TARGFLUXDATA:
c
c     This routine calculates the heat and particle flux components onto each 
c     element of the target - it saves all of the elements as well
c     as the total for later print out. At present, it only includes
c     the effects of hydrogen ion and neutral fluxes - it does not include
c     either radiation or impurity energy fluxes.  
c
c
      integer ik,ir,in,id,wallind
      real ts,cs
      real gai,gae
c
      call rzero(targfluxdata,(maxnds+3)*4*4)
c
c     Loop through target elements
c
      do id = 1,nds
c
c        Do not calculate for invalid target elements
c
         if (dds(id).eq.0.0) cycle
c
c        Setup indexing
c
         ik = ikds(id)
         ir = irds(id)
c 
c        Set target identifier
c
         if (id.le.ndsin) then 
            in = 1
         else
            in = 2
         endif
c
c        Extract actual target speed (TS) and calculate sound speed (CS)
c
         TS = abs(kvds(id))
c
         CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
     >                 (1.0+RIZB)/CRMB) * cmachno(ir,in)
c
c        Calculate particle flux on element 
c
c        Ion particle flux 
c
c         targfluxdata(id,1,1) = KNDS(ID) * TS/KBFS(IK,IR) * COSTET(ID)
c
         targfluxdata(id,1,1) = KNDS(ID) * TS/KBFST(IR,IN) * COSTET(ID)
c
c        Calculate ion heat coefficients
c
         gae = 5.0
         gai = 2.5
c
c         gai = 2.5 +  0.5 * (ts/cs)**2 * (1.0 + kteds(id)/ktids(id))
c
c        Ion heat flux 
c
         targfluxdata(id,1,2) = (gae * rizb * kteds(id) +
     >                           gai * ktids(id)) * ech *
     >                           targfluxdata(id,1,1)
c
c        Ion recombination and potential energy
c
         targfluxdata(id,1,3) = (13.6 + 2.2) * ech *
     >                           targfluxdata(id,1,1)
c
c        Total for ions
c         
         targfluxdata(id,1,4) = targfluxdata(id,1,2) +
     >                          targfluxdata(id,1,3)



        write(6,'(a,i6,10(1x,g12.5))') 'TARGFLUX:',id,
     >       dds(id),dds2(id),targfluxdata(id,1,1),
     >       targfluxdata(id,1,2),knds(id),ts,costet(id),
     >       kbfst(ir,in)  




c
c         write(6,*) 'ID:',id,nimindex(id)
c
c        Atom and molecular fluxes 
c
         if (nimindex(id).ne.0) then 
c
            wallind = nimindex(id)
c
c slmod begin
            if (pincode.eq.2.or.pincode.eq.3.or.pincode.eq.4.or.
     .          pincode.eq.5) then
              IF (id.EQ.1) THEN               
                write(0,*) ! Look below for the same message...
                write(0,*) '-------------------------------'
                write(0,*) ' NOT SURE FLXHW USED CORRECTLY'
                write(0,*) '-------------------------------'
                write(0,*) 
              ENDIF
c
c...           Using the ion flux returned by EIRENE (FLXHW8).  This
c              quantity is available for EIRENE99 only:
c
c              targfluxdata(id,2,1) = max(flxhw2(wallind)-
c     >                                   flxhw8(wallind),0.0)
c
c jdemod -     flxhw6 contains the atomic flux directly - seems better to 
c              use it that the derived quantity flxhw2.
c
               targfluxdata(id,2,1) = max(flxhw6(wallind),0.0)
c
            else
c
               targfluxdata(id,2,1) = max(flxhw2(wallind)
     >                           -targfluxdata(id,1,1),0.0)
c
            endif
c
c            targfluxdata(id,2,1) = max(flxhw2(wallind)
c     >                         -targfluxdata(id,1,1),0.0)
c slmod end
c
c           Calculate atom heat
c
            targfluxdata(id,2,2) =  flxhw5(wallind) * ech *
     >                              targfluxdata(id,2,1)
c
c           Calcuate atom to molecule recombination energy
c
            targfluxdata(id,2,3) =  2.2 * ech *
     >                              targfluxdata(id,2,1)
c
c           Total for atoms
c         
            targfluxdata(id,2,4) = targfluxdata(id,2,2) +
     >                             targfluxdata(id,2,3)
c
c           Molecule particle flux
c 
            targfluxdata(id,3,1) = max(fluxhw(wallind)
     >                         -targfluxdata(id,2,1),0.0)
c
c           Calculate molecule heat
c
            targfluxdata(id,3,2) =  flxhw7(wallind) * ech *
     >                              targfluxdata(id,3,1)
c
c           Total for molecules
c         
            targfluxdata(id,3,4) = targfluxdata(id,3,2) 
c
c
c           Comparing DIVIMP and Eirene fluxes: 
c
            write(6,'(a,2i5,1x,f7.3,10(1x,g12.5))') 'D/E WFLUX:',
     >           id,nimindex(id),
     >        (targfluxdata(id,1,1)-flxhw8(wallind))
     >                 /targfluxdata(id,1,1)*100.0,
     >         targfluxdata(id,1,1),flxhw8(wallind),
     >         targfluxdata(id,2,1),flxhw6(wallind),
     >         targfluxdata(id,3,1),fluxhw(wallind)-flxhw6(wallind)
c
         endif
c
      end do  
c
c     Sum up these arrays to get inner and outer target as well as total 
c     contributions.
c
c     Target 1
c
      do id = 1,ndsin
c
         do in = 1,3

            targfluxdata(id,4,1) = targfluxdata(id,4,1) + 
     >                             targfluxdata(id,in,1)
            targfluxdata(id,4,2) = targfluxdata(id,4,2) + 
     >                             targfluxdata(id,in,2)
            targfluxdata(id,4,3) = targfluxdata(id,4,3) + 
     >                             targfluxdata(id,in,3)
            targfluxdata(id,4,4) = targfluxdata(id,4,4) + 
     >                             targfluxdata(id,in,4)
         end do 
c
         do in = 1,4

            targfluxdata(maxnds+1,1,in)=
     >                              targfluxdata(maxnds+1,1,in)+
     >                              targfluxdata(id,1,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+1,2,in)=
     >                              targfluxdata(maxnds+1,2,in)+
     >                              targfluxdata(id,2,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+1,3,in)=
     >                              targfluxdata(maxnds+1,3,in)+
     >                              targfluxdata(id,3,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+1,4,in)=
     >                              targfluxdata(maxnds+1,4,in)+
     >                              targfluxdata(id,4,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
         end do
c
      end do 
c
c     Target 2 
c

      do id = ndsin+1,nds
c
         do in = 1,3

            targfluxdata(id,4,1) = targfluxdata(id,4,1) + 
     >                             targfluxdata(id,in,1)
            targfluxdata(id,4,2) = targfluxdata(id,4,2) + 
     >                             targfluxdata(id,in,2)
            targfluxdata(id,4,3) = targfluxdata(id,4,3) + 
     >                             targfluxdata(id,in,3)
            targfluxdata(id,4,4) = targfluxdata(id,4,4) + 
     >                             targfluxdata(id,in,4)
         end do 
c
         do in = 1,4

            targfluxdata(maxnds+2,1,in)=
     >                              targfluxdata(maxnds+2,1,in)+
     >                              targfluxdata(id,1,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+2,2,in)=
     >                              targfluxdata(maxnds+2,2,in)+
     >                              targfluxdata(id,2,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+2,3,in)=
     >                              targfluxdata(maxnds+2,3,in)+
     >                              targfluxdata(id,3,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
            targfluxdata(maxnds+2,4,in)=
     >                              targfluxdata(maxnds+2,4,in)+
     >                              targfluxdata(id,4,in)*dds(id)
     >                              * 2.0 * PI * rp(id)
         end do
c
      end do 
c
c     Sum over both target totals
c
      do in = 1,4

         targfluxdata(maxnds+3,1,in) = 
     >                             targfluxdata(maxnds+1,1,in)+
     >                             targfluxdata(maxnds+2,1,in)
         targfluxdata(maxnds+3,2,in) = targfluxdata(maxnds+1,2,in)+
     >                             targfluxdata(maxnds+2,2,in)
         targfluxdata(maxnds+3,3,in) = targfluxdata(maxnds+1,3,in)+
     >                             targfluxdata(maxnds+2,3,in)
         targfluxdata(maxnds+3,4,in) = targfluxdata(maxnds+1,4,in)+
     >                             targfluxdata(maxnds+2,4,in)
      end do 

c
c     Print summary to unit 6 
c
      write(6,'(a)') 'Target Flux Summary:'
      write(6,'(a)') 'TYPE   ID   FLUX          HEAT  '//
     >            '       REC         TOTAL'
c
      do id = 1,ndsin
         write(6,'(a6,i4,4(1x,g12.5))') 'ION :',id,
     >             (targfluxdata(id,1,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'ATOM:',id,
     >             (targfluxdata(id,2,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'MOLE:',id,
     >             (targfluxdata(id,3,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'TOT :',id,
     >             (targfluxdata(id,4,in),in=1,4)
      end do
c
      write(6,'(a10,4(1x,g12.5))') 'TOTI '//inner,
     >             (targfluxdata(maxnds+1,1,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTA '//inner,
     >             (targfluxdata(maxnds+1,2,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTM '//inner,
     >             (targfluxdata(maxnds+1,3,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOT  '//inner,
     >             (targfluxdata(maxnds+1,4,in),in=1,4)
c
      do id = ndsin+1,nds
         write(6,'(a6,i4,4(1x,g12.5))') 'ION :',id,
     >             (targfluxdata(id,1,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'ATOM:',id,
     >             (targfluxdata(id,2,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'MOLE:',id,
     >             (targfluxdata(id,3,in),in=1,4)
         write(6,'(a6,i4,4(1x,g12.5))') 'TOT :',id,
     >             (targfluxdata(id,4,in),in=1,4)
      end do

      write(6,'(a10,4(1x,g12.5))') 'TOTI '//outer,
     >             (targfluxdata(maxnds+2,1,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTA '//outer,
     >             (targfluxdata(maxnds+2,2,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTM '//outer,
     >             (targfluxdata(maxnds+2,3,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOT  '//outer,
     >             (targfluxdata(maxnds+2,4,in),in=1,4)
c
      write(6,'(a)') 'GRAND TOTALS:' 
c
      write(6,'(a10,4(1x,g12.5))') 'TOTI '//outer,
     >             (targfluxdata(maxnds+3,1,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTA '//outer,
     >             (targfluxdata(maxnds+3,2,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOTM '//outer,
     >             (targfluxdata(maxnds+3,3,in),in=1,4)
      write(6,'(a10,4(1x,g12.5))') 'TOT  '//outer,
     >             (targfluxdata(maxnds+3,4,in),in=1,4)
c
      return
      end 

c
c
c
      subroutine calc_wallfluxdata
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'pindata'
      include 'printopt'
c slmod begin
      include 'slcom'
c slmod end
c
c     CALC_WALLFLUXDATA:
c
c     This routine calculates the heat and particle flux components onto each 
c     element of the wall - it saves all of the elements as well
c     as the total for later print out. At present, it only includes
c     the effects of hydrogen ion and neutral fluxes - it does not include
c     either radiation or impurity energy fluxes.  
c
c     Contents of the wallfluxdata array
c     
c     wallfluxdata(in,is,it) 
c
c     in = wall element
c     is = species        1 = ions 
c                         2 = atoms
c                         3 = molecules
c                         4 = total
c     it = type of flux   1 = particle flux
c                         2 = direct heat flux
c                         3 = recombination and potential energy
c                         4 = total  
c
c     All data is in W/m^2 until grand totals are calculated then the 
c     units are in W (for the whole toroid)
c
c
      integer ik,ir,in,id,it,itarg,is,segtype
      integer wallid,targid
      real ts,cs
      real gai,gae
c
      call rzero(wallfluxdata,(maxpts+5)*4*4)
c
c     Loop through wall elements
c
      do in = 1,wallpts
c
c        Extract wall and target indices from the wall data array
c        Target indices are non-zero only for REAL target elements. 
c
         wallid = int(wallpt(in,17))
         targid = int(wallpt(in,18))
c
c        Calculate ION components if this is a target element.
c
         if (targid.ne.0) then 
c
c           Setup indexing
c
            ik = ikds(targid)
            ir = irds(targid)
c 
c           Set target identifier
c
            if (targid.le.ndsin) then 
               itarg = 1
            else
               itarg = 2
            endif
c
c           Extract actual target speed (TS) and calculate sound speed (CS)
c
            TS = abs(kvds(targid))
c
            CS = 9.79E3 * SQRT (0.5*(KTEDS(targID)+KTIDS(targID))*
     >                 (1.0+RIZB)/CRMB) * cmachno(ir,itarg)
c
c           Calculate particle flux on element 
c
c           Ion particle flux 
c
            wallfluxdata(in,1,1) = KNDS(targID) 
     >                           * TS/KBFST(IR,Itarg) * COSTET(targid)
c
c           Calculate ion heat coefficients
c
            gae = 5.0
            gai = 2.5
c
c           gai = 2.5 +  0.5 * (ts/cs)**2 * (1.0 + kteds(id)/ktids(id))
c
c           Ion heat flux 
c
            wallfluxdata(in,1,2) = (gae * rizb * kteds(targid) +
     >                           gai * ktids(targid)) * ech *
     >                           wallfluxdata(in,1,1)
c
c           Ion recombination and potential energy
c
            wallfluxdata(in,1,3) = (13.6 + 2.2) * ech *
     >                           wallfluxdata(in,1,1)
c
c           Total for ions
c         
            wallfluxdata(in,1,4) = wallfluxdata(in,1,2) +
     >                             wallfluxdata(in,1,3)
c
         endif 
c
c        Atom and molecular fluxes 
c
         if (wallid.ne.0) then 
c
c slmod begin
            if (pincode.eq.2.or.pincode.eq.3.or.pincode.eq.4.or.
     .          pincode.eq.5) then
c              IF (in.EQ.1) THEN
c                write(0,*) 
c                write(0,*) '-------------------------------'
c                write(0,*) ' NOT SURE FLXHW USED CORRECTLY'
c                write(0,*) '-------------------------------'
c                write(0,*) 
c              ENDIF
c
c...           Using the ion flux returned by EIRENE (FLXHW8).  This
c              quantity is available for EIRENE99 only:
c
c              wallfluxdata(id,2,1) = max(flxhw2(wallid)-
c     >                                   flxhw8(wallid),0.0)
c
c jdemod -     flxhw6 contains the atomic flux directly - seems better to 
c              use it that the derived quantity flxhw2.
c
               wallfluxdata(in,2,1) = max(flxhw6(wallid),0.0)
c
            else
c
               wallfluxdata(in,2,1) = max(flxhw2(wallid)
     >                           -wallfluxdata(in,1,1),0.0)
c
            endif
c
c slmod end
c
c           Calculate atom heat
c
            wallfluxdata(in,2,2) =  flxhw5(wallid) * ech *
     >                              wallfluxdata(in,2,1)
c
c           Calcuate atom to molecule recombination energy
c
            wallfluxdata(in,2,3) =  2.2 * ech *
     >                              wallfluxdata(in,2,1)
c
c           Total for atoms
c         
            wallfluxdata(in,2,4) = wallfluxdata(in,2,2) +
     >                             wallfluxdata(in,2,3)
c
c           Molecule particle flux
c 
            wallfluxdata(in,3,1) = max(fluxhw(wallid)
     >                         -wallfluxdata(in,2,1),0.0)
c
c           Calculate molecule heat
c
            wallfluxdata(in,3,2) =  flxhw7(wallid) * ech *
     >                              wallfluxdata(in,3,1)
c
c           Total for molecules
c         
            wallfluxdata(in,3,4) = wallfluxdata(in,3,2) 
c
         endif
c
      end do  
c
c     Sum up these arrays to get inner and outer target, pfz wall,
c     main wall and total contributions.
c
      do in = 1,wallpts
c
c        Sum total fluxes of each type for each element. 
c
         do is = 1,3
c
c           Sum up each type of flux
c
            do it = 1,3
c
               wallfluxdata(in,4,it) = wallfluxdata(in,4,it)
     >                               + wallfluxdata(in,is,it)
c
            end do
c
            wallfluxdata(in,4,4) = wallfluxdata(in,4,4)
     >                           + wallfluxdata(in,is,4)
c
         end do 
c
c
c        The following grand totals are calculated on a complete torus
c        basis factoring in the major radius. 
c
c
c        Check the tag for the wall element to see which type of 
c        segment it is:
c
c        1 = Target 1 (Outer target for X-point up grids)
c        4 = Target 2 (Inner target for X-point up grids)
c        7 = (and others 2,3) Main Vessel wall
c        8 = PFZ wall
c        9,10 = baffle segments - added to PFZ only applies 
c                                 for some JET grids 
c
         segtype = wallpt(in,16)
c
c        First Target (Inner for X-point up - Outer for down)
c  
         if (segtype.eq.4) then          
c
            do it = 1,4
               do is = 1,4
c
                  wallfluxdata(maxpts+1,is,it)=
     >                              wallfluxdata(maxpts+1,is,it)+
     >                              wallfluxdata(in,is,it)*wallpt(in,7)
     >                              * 2.0 * PI * wallpt(in,1)
               end do
            end do   
c
c        Second Target (Outer for X-point up - Inner for down)
c  
         elseif (segtype.eq.1) then          
c
            do it = 1,4
               do is = 1,4

                  wallfluxdata(maxpts+2,is,it)=
     >                              wallfluxdata(maxpts+2,is,it)+
     >                              wallfluxdata(in,is,it)*wallpt(in,7)
     >                              * 2.0 * PI * wallpt(in,1)
               end do
            end do   
c
c        Main Vessel wall elements
c
         elseif (segtype.eq.7.or.segtype.eq.2.or.segtype.eq.3) then 
c
            do it = 1,4
               do is = 1,4

                  wallfluxdata(maxpts+3,is,it)=
     >                              wallfluxdata(maxpts+3,is,it)+
     >                              wallfluxdata(in,is,it)*wallpt(in,7)
     >                              * 2.0 * PI * wallpt(in,1)
               end do
            end do   
c
c        Private Flux Zone wall elements   
c
         elseif (segtype.eq.8.or.segtype.eq.9.or.segtype.eq.10) then
c
            do it = 1,4
               do is = 1,4

                  wallfluxdata(maxpts+4,is,it)=
     >                              wallfluxdata(maxpts+4,is,it)+
     >                              wallfluxdata(in,is,it)*wallpt(in,7)
     >                              * 2.0 * PI * wallpt(in,1)
               end do
            end do   
c
         endif
c
      end do 
c
c     Create grand totals summing over all regions 
c
      do ir = 1,4
c
         do it = 1,4
c
            do is = 1,4 

               wallfluxdata(maxpts+5,is,it) = 
     >                             wallfluxdata(maxpts+5,is,it) +
     >                             wallfluxdata(maxpts+ir,is,it)
            end do
         end do
      end do

c
c     Print summary to unit 6 
c
      write(6,*)
      write(6,'(a)') 'Wall Flux Summary:'
      write(6,'(a)') 'TYPE   ID TARG#  FLUX  '//
     >            '        HEAT  '//
     >            '       REC         TOTAL'
c
c     Wall elements
c
      do in = 1,wallpts
         write(6,'(a6,2i4,1p,4(1x,e12.5))') 'ION :',in,
     >             int(wallpt(in,18)),
     >             (wallfluxdata(in,1,it),it=1,4)
         write(6,'(a6,2i4,1p,4(1x,e12.5))') 'ATOM:',in,
     >             int(wallpt(in,18)),
     >             (wallfluxdata(in,2,it),it=1,4)
         write(6,'(a6,2i4,1p,4(1x,e12.5))') 'MOLE:',in,
     >             int(wallpt(in,18)),
     >             (wallfluxdata(in,3,it),it=1,4)
         write(6,'(a6,2i4,1p,4(1x,e12.5))') 'TOT :',in,
     >             int(wallpt(in,18)),
     >             (wallfluxdata(in,4,it),it=1,4)
      end do
c
      write(6,*)   
c
c     First target
c
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTI '//inner,
     >             (wallfluxdata(maxpts+1,1,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTA '//inner,
     >             (wallfluxdata(maxpts+1,2,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTM '//inner,
     >             (wallfluxdata(maxpts+1,3,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOT  '//inner,
     >             (wallfluxdata(maxpts+1,4,it),it=1,4)
c
c     Second target
c
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTI '//outer,
     >             (wallfluxdata(maxpts+2,1,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTA '//outer,
     >             (wallfluxdata(maxpts+2,2,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTM '//outer,
     >             (wallfluxdata(maxpts+2,3,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOT  '//outer,
     >             (wallfluxdata(maxpts+2,4,it),it=1,4)
c
c     Main Wall
c
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTI  MAIN',
     >             (wallfluxdata(maxpts+3,1,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTA  MAIN',
     >             (wallfluxdata(maxpts+3,2,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTM  MAIN',
     >             (wallfluxdata(maxpts+3,3,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOT   MAIN',
     >             (wallfluxdata(maxpts+3,4,it),it=1,4)
c
c     PFZ Wall
c
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTI   PFZ',
     >             (wallfluxdata(maxpts+4,1,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTA   PFZ',
     >             (wallfluxdata(maxpts+4,2,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTM   PFZ',
     >             (wallfluxdata(maxpts+4,3,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOT    PFZ',
     >             (wallfluxdata(maxpts+4,4,it),it=1,4)
c
      write(6,'(a)') 'GRAND TOTALS:' 
c
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTI GRAND',
     >             (wallfluxdata(maxnds+5,1,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTA GRAND',
     >             (wallfluxdata(maxnds+5,2,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOTM GRAND',
     >             (wallfluxdata(maxnds+5,3,it),it=1,4)
      write(6,'(a10,1p,4(1x,e12.5))') 'TOT  GRAND',
     >             (wallfluxdata(maxnds+5,4,it),it=1,4)
c
      return
      end 
c
c 
c
      subroutine calc_wall_intersections(n_int,max_int,r_int,z_int,
     >                       rstart,zstart,rend,zend,ignore_end)
      implicit none
      integer n_int,max_int
      real*8 r_int(max_int),z_int(max_int)
      real*8 rstart,zstart,rend,zend
      logical ignore_end
c
      include 'params'
      include 'comtor'
c
c     CALC_WALL_INTERSECTIONS:
c
c     This routine loops through the entire wall an reports all
c     of the intersections of the given line (rstart,zstart) (rend,zend) 
c     with the wall. The ignore_end flag instructs the code to not
c     calculate any intersections for the end-point since it is known
c     to be on the wall. 
c
      integer in
      real*8 rnew,znew,tnew,tnorm
      logical sect
c
      real eps
      parameter (eps=1.0e-6)
c
      n_int = 0
      call rzero(r_int,max_int)
      call rzero(z_int,max_int)
c
      do in = 1,wallpts
c
c         Initialize for call - most not needed
c
          rnew = 0.0
          znew = 0.0         
          tnew = 0.0
          tnorm= 0.0
          sect = .false.
c
          CALL INTCALCDP(Rend,Zend,Rstart,Zstart,
     >                 dble(WALLPT(IN,1)),dble(WALLPT(IN,2)),
     >                 dble(WALLPT(IN,8)),dble(WALLPT(IN,9)),
     >                 dble(WALLPT(IN,5)),dble(WALLPT(IN,6)),
     >                 RNEW,ZNEW,TNEW,tnorm,
     >                 SECT,nrfopt)
c
c         Record any intersections found
c       
          if (sect) then 
c
c            If not (endpoint ignored and this is the endpoint)
c
             if (.not.(ignore_end.and.
     >             (abs(rnew-rend).lt.eps.and.
     >              abs(znew-zend).lt.eps))) then 
c
c               Record point
c
                if (n_int.lt.max_int) then 
                   n_int = n_int + 1
                   r_int(n_int) = rnew
                   z_int(n_int) = znew   
                endif 
c
             endif
c
          endif  
c
      end do
c  
      return
      end
c
c      
c
       subroutine calc_wallprad(nizs)
       implicit none
c
       integer nizs 
c
       include 'params'
       include 'cgeom'
       include 'comtor'
       include 'dynam3'
       include 'printopt'
c
c      CALC_WALLPRAD:
c
c      This routine calculates the total radiated power onto each wall segment.
c      It also calculates the total radiated H and IMP power. 
c      In order to calculate the amount on each wall segment it has to 
c      allocate the amount from each cell that would strike each wall elememt.
c      This routine uses the data in HPOWLS and POWLS to make these calculations.
C      The POWLS data is multiplied by the ABSFAC quantity to get absolute scalings
c      unless ABSFAC is zero. 
c
       integer in,ir,ik,iz,it
       real*8 imp_prad,h_prad,tot_prad
       real*8 rcent,zcent
c
       integer n_int,max_int,backward,forward,ivert
       parameter(max_int=10,backward=-1,forward=1)
       real*8 r_int(max_int),z_int(max_int)
       integer n_intersections(maxpts+1),next_vert
       external next_vert 
c
       real*8 datan2c
       real*8 wall_angle(maxpts+1),net_angle(maxpts),result_angle
       external datan2c
c
       integer nstart,nend,ntest,segtype
       real*8 rtest,ztest
c
c      Initialization 
c
       call rzero(wallprad,(maxpts+6)*3)
c
c      For every cell on the grid - this code must loop through every element of 
c      the wall.  
c       
       do ir = 1,nrs
          do ik = 1,nks(ir)
c
c            Sum up impurity and hydrogen radiation contributions
c              
c            Impurity
c
             imp_prad = 0.0
c 
             do iz = 0,nizs
c
                imp_prad = imp_prad + powls(ik,ir,iz)
c
             end do
c
             if (absfac.gt.0.0) imp_prad = imp_prad * absfac
c
c            Hydrogen
c
             h_prad = 0.0
c 
             do iz = 0,1
c
                h_prad = h_prad + hpowls(ik,ir,iz)
c
             end do
c
c            Scale for cell area and major radius - converting to watts
c
             imp_prad = imp_prad * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
             h_prad   = h_prad   * kareas(ik,ir) * 2.0 * PI * rs(ik,ir) 
             tot_prad = h_prad + imp_prad
c
c            Don't perform geometry calculations unless there is a radiation contribution 
c
             if (tot_prad.eq.0.0) cycle
c
c            Record totals of radiated source in highest array element - these
c            will be compared to totals integrated over all wall elements
c            and should be the same.
c
             wallprad(maxpts+6,1) = wallprad(maxpts+6,1) + imp_prad
             wallprad(maxpts+6,2) = wallprad(maxpts+6,2) + h_prad
             wallprad(maxpts+6,3) = wallprad(maxpts+6,3) + tot_prad
c
c            Calculate the geometry setup information - can each wall vertex be 
c            seen from the current cell and what are the angles to all of the 
c            vertices.
c 
             rcent = rs(ik,ir)
             zcent = zs(ik,ir)  
c
             do in = 1,wallpts
c
c               Calculate initial intersections and whether LOS to each
c               vertex is clear. Last point and first point of the wall
c               should be identical
c
                n_intersections(in) = 0
c
                rtest = wallpt(in,20)
                ztest = wallpt(in,21)
c
                call calc_wall_intersections(ntest,max_int,r_int,z_int,
     >                       rcent,zcent,rtest,ztest,.true.)
c             
                n_intersections(in) = ntest 
c
c               Calculate wall_angle array
c
c               Get angles to the start and end of the wall segment
c 
                wall_angle(in) = datan2c(ztest-zcent,rtest-rcent)
c
c               Convert angles to 0 to 2 PI range  from -PI to PI
c
                if (wall_angle(in).lt.0.0)
     >              wall_angle(in)=wall_angle(in)+2.0*PI
c
             end do
c
c            Assign values to complete the wall
c
             n_intersections(wallpts+1) = n_intersections(1)
             wall_angle(wallpts+1) = wall_angle(1)
c
c            Calculate the net angle for each segment of the wall
c
             do in = 1,wallpts
c
                net_angle(in) = wall_angle(in+1)- wall_angle(in)
c
                if (net_angle(in).gt.PI) 
     >              net_angle(in) = net_angle(in) - 2.0 * PI
                if (net_angle(in).lt.-PI) 
     >              net_angle(in) = net_angle(in) + 2.0 * PI 
c
             end do 
c
c  
c            Basic LOS rules:
c	     
c            1) Wall is listed in a basically clockwise order so assuming that 
c               all angles are positive in the range 0.0 to 2 PI then the angle
c               to the leading point of a wall segment will always be less 
c               than the angle to the second point for a forward oriented wall 
c               segment. 
c            2) All wall elements that have some exposure to radiation from the
c               cell being examined will be forward going assuming all radiation 
c               is from within the closed vessel. 
c            3) If the LOS from the cell to the ends of the wall segment cross the
c               the wall at any point - then the view to that wall segment is 
c               assumed to be blocked.
c            4) If the LOS from the cell to the ends of the wall segment does
c               not cross the wall then the LOS is considered unblocked and the
c               entire segment is visible to the source. 
c            5) If either end of a backward going wall segment is blocked then 
c               none of that line segment can be seen from the cell. 
c            6) If the wall segment is forward going but one vertex is blocked then
c               the wall element is partially obscured and the proportion of the
c               segment exposed to the source cell is approximately calculated.  
c            7) Only clockwise oriented segments can receive radiation (net_angle<0)
c            8) The total exposure of a partially blocked segment is limited by
c               the visible vertex and either the next or last visible vertex along 
c               the wall. 
c 
c

             do in = 1,wallpts
c
                if (net_angle(in).lt.0.0) then 
c
c                  If net_angle is less than zero then the wall element is 
c                  forward going or clockwise relative to the observation point.
c              
c                  Check number of intersections for the wall element vertices  
c
c                  Totals only need updating for unobstructed and partially 
c                  obstructed views.
c
c              
c                  Wall segment is unobstructed 
c
                   if (n_intersections(in).eq.0.and.
     >                 n_intersections(in+1).eq.0) then 
c
                      result_angle = -net_angle(in) 
c
c
c                  Clockwise wall segment with one blocked vertex.
c                  Need to calculate the required angle
c                  Need angle of last unobstructed vertex
c
                   elseif (n_intersections(in).gt.0.and.
     >                    n_intersections(in+1).eq.0) then 
c		    
                      ivert = next_vert(n_intersections,
     >                                  wallpts+1,in,backward)
c		    
                      result_angle = wall_angle(in+1)- wall_angle(ivert)
c		    
                      if (result_angle.gt.PI) 
     >                    result_angle = result_angle - 2.0 * PI
                      if (result_angle.lt.-PI) 
     >                    result_angle = result_angle + 2.0 * PI 

                      result_angle = -result_angle
c
c
c                  Need angle of next unobstructed vertex
c
                   elseif (n_intersections(in).eq.0.and.
     >                    n_intersections(in+1).gt.0) then
		  
                      ivert = next_vert(n_intersections,
     >                                  wallpts+1,in+1,forward)
c		  
                      result_angle = wall_angle(ivert)- wall_angle(in)
c		  
                      if (result_angle.gt.PI) 
     >                    result_angle = result_angle - 2.0 * PI
                      if (result_angle.lt.-PI) 
     >                    result_angle = result_angle + 2.0 * PI 

                      result_angle = -result_angle
c
c
c                  LOS obstructed
c
                   elseif (n_intersections(in).gt.0.and.
     >                     n_intersections(in+1).gt.0) then
c
                       result_angle = 0.0 
c
                   endif
c
c                  Calculate contributions to this wall segment
c
                   wallprad(in,1) = wallprad(in,1) 
     >                         + imp_prad * result_angle/(2.0*PI)
                   wallprad(in,2) = wallprad(in,2) 
     >                         + h_prad * result_angle/(2.0*PI)
                   wallprad(in,3) = wallprad(in,3) 
     >                         + tot_prad * result_angle/(2.0*PI)
                end if

             end do
c
         end do
c
      end do  
c
c     Print outs of results
c
      do in = 1,wallpts
c
c        The following grand totals are calculated on a complete torus
c        basis factoring in the major radius. 
c
c
c        Check the tag for the wall element to see which type of 
c        segment it is:
c
c        1 = Target 1 (Outer target for X-point up grids)
c        4 = Target 2 (Inner target for X-point up grids)
c        7 = (and others 2,3) Main Vessel wall
c        8 = PFZ wall
c        9,10 = baffle segments - added to PFZ only applies 
c                                 for some JET grids 
c
         segtype = wallpt(in,16)
c
c        First Target (Outer for X-point up - inner for down)
c  
         if (segtype.eq.1) then          
c
            wallprad(maxpts+1,1) = wallprad(maxpts+1,1) + wallprad(in,1)
            wallprad(maxpts+1,2) = wallprad(maxpts+1,2) + wallprad(in,2)
            wallprad(maxpts+1,3) = wallprad(maxpts+1,3) + wallprad(in,3)
c
c        Second Target (Inner for X-point up - Outer for down)
c  
         elseif (segtype.eq.4) then          
c
            wallprad(maxpts+2,1) = wallprad(maxpts+2,1) + wallprad(in,1)
            wallprad(maxpts+2,2) = wallprad(maxpts+2,2) + wallprad(in,2)
            wallprad(maxpts+2,3) = wallprad(maxpts+2,3) + wallprad(in,3)
c
c        Main Vessel wall elements
c
         elseif (segtype.eq.7.or.segtype.eq.2.or.segtype.eq.3) then
c
            wallprad(maxpts+3,1) = wallprad(maxpts+3,1) + wallprad(in,1)
            wallprad(maxpts+3,2) = wallprad(maxpts+3,2) + wallprad(in,2)
            wallprad(maxpts+3,3) = wallprad(maxpts+3,3) + wallprad(in,3)
c
c        Private Flux Zone wall elements   
c
         elseif (segtype.eq.8.or.segtype.eq.9.or.segtype.eq.10) then
c
            wallprad(maxpts+4,1) = wallprad(maxpts+4,1) + wallprad(in,1)
            wallprad(maxpts+4,2) = wallprad(maxpts+4,2) + wallprad(in,2)
            wallprad(maxpts+4,3) = wallprad(maxpts+4,3) + wallprad(in,3)
c
         endif
c
      end do 
c
c     Sum up grand totals
c
      do in = 1,4
c 
         do it = 1,3

            wallprad(maxpts+5,it) = wallprad(maxpts+5,it) 
     >                            + wallprad(maxpts+in,it)
c
         end do
c
      end do
c
c     Convert all of the wall element data to W/m2
c
      do in = 1,wallpts
         do it = 1,3
            if (wallpt(in,7).gt.0.0.and.wallpt(in,1).gt.0.0) then
               wallprad(in,it) = wallprad(in,it)
     >                           /(wallpt(in,7)*2.0*PI*wallpt(in,1))
            endif
         end do
      end do
c
c
c     Print summary to unit 6 
c
      write(6,*) 
      write(6,'(a)') 'WALL Radiation Flux Summary [W/m2]:'
      write(6,'(a)') ' ID   TYPE    IMPURITY     HYDROGEN     TOTAL'
c     Wall elements
c
      do in = 1,wallpts
         write(6,'(2(i4,2x),4(1x,g12.5))') in,int(wallpt(in,16)),
     >             (wallprad(in,it),it=1,3)
      end do
c
      write(6,*) 
c
c     First target
c
      write(6,*) 'Regional Radiation Totals in [W]'

      write(6,'(a10,4(1x,g12.5))') 'TOT  '//inner,
     >             (wallprad(maxpts+1,it),it=1,3)
      write(6,'(a10,4(1x,g12.5))') 'TOT  '//outer,
     >             (wallprad(maxpts+2,it),it=1,3)
      write(6,'(a10,4(1x,g12.5))') 'TOT   MAIN',
     >             (wallprad(maxpts+3,it),it=1,3)
      write(6,'(a10,4(1x,g12.5))') 'TOT    PFZ',
     >             (wallprad(maxpts+4,it),it=1,3)
      write(6,'(a10,4(1x,g12.5))') 'TOTAL  SEG',
     >             (wallprad(maxpts+5,it),it=1,3)
      write(6,'(a10,4(1x,g12.5))') 'TOTAL  SRC',
     >             (wallprad(maxpts+6,it),it=1,3)


      return
      end
c
c
c
      integer function next_vert(n_intersections,maxp,startin,step)
      implicit none
      integer startin,maxp,step
      integer n_intersections(maxp)
c
c     NEXT_VERT: This routine finds the next entry in n_intersections
c                which is equal to zero - moving in the direction
c                defined by step and starting at startin. 
c
c
      integer current
c
      current = startin
c  
      do while (n_intersections(current).ne.0)
c
         current = current + step
c          
         if (current.gt.maxp) current = 1
         if (current.lt.1) current = maxp
c
      end do
c 
      next_vert = current
c
      return
      end
c
c
c
      subroutine assign_wall_plasma
      implicit none
      include 'params'
      include 'walls_com'
      include 'cgeom'
c
c
c     ASSIGN_WALL_PLASMA:
c
c     This routine associates some plasma conditions with each wall
c     element. These conditions will only be used by certain specific 
c     options and do not imply the existence of such a plasma. 
c 
c     Option 0: default - simply associates the plasma conditions in 
c     the nearest cell with the wall. Across the target these data
c     are replaced with the target conditions.
c
c     Option 1: Linear decay over a specified decay length. Minimum values
c     are imposed.
c
c     Option 2: Exponential decay with specified factor. 
c
      integer in,id,ikp,irp
      real te_min,ti_min,ne_min
      parameter(te_min=0.5,ti_min=0.5,ne_min=1.0e15)
c
c     Loop around wall
c
      do in = 1,wallpts
c
c        Check for target elements
c
         if (wallpt(in,16).eq.1.or.wallpt(in,16).eq.4) then
c
            if (wallpt(in,18).ne.0) then 
c
               id = wallpt(in,18)
c
               wallpt(in,29) = kteds(id)
               wallpt(in,30) = ktids(id)
               wallpt(in,31) = knds(id)
c
            else
c
               wallpt(in,29) = te_min
               wallpt(in,30) = ti_min
               wallpt(in,31) = ne_min
c
            endif
c
         else
c
c           Default option
c 
            if (wall_plasma_opt.eq.0) then 
c
               ikp = wallpt(in,26)
               irp = wallpt(in,27)
c
               wallpt(in,29) = ktebs(ikp,irp)
               wallpt(in,30) = ktibs(ikp,irp)
               wallpt(in,31) = knbs(ikp,irp)
c
c           Linear or exponential from grid with interpolation
c
            elseif (wall_plasma_opt.eq.1.or.
     >              wall_plasma_opt.eq.2) then 
c
c               Put in subroutine 
c
                call calc_wall_plasma(in,te_min,ti_min,ne_min)
c
            endif  
c
         endif 
c
      end do
c
      return
      end 
c
c
c
      subroutine calc_wall_plasma(in,te_min,ti_min,ne_min)
      implicit none
      integer in
      real te_min,ti_min,ne_min 
      include 'params'
      include 'walls_com'
      include 'cgeom'  
c
c     CALC_WALL_PLASMA: interpolate and extend calculation of the 
c     wall plasma. 
c
c     First - calculate interpolated plasma conditions at ring intersection
c     point.     
c
      real ri,zi
      real r1,z1,r2,z2,s,st,sn
      real d1,d2,fact
      real ne,te,ti,nen,ten,tin  
      integer ikp,irp,id
c
c     Load wall intersection point calculated in the WALLS.f routine - entries will be overwritten with 
c     the wall plasma conditions.
c
      ri = wallpt(in,29)
      zi = wallpt(in,30)
c
c     Plasma cell
c
      ikp = wallpt(in,26)
      irp = wallpt(in,27) 
c
c     Base plasma conditions 
c
      ne = knbs(ikp,irp)
      te = ktebs(ikp,irp)
      ti = ktibs(ikp,irp)
      st = kss(ikp,irp) 
c
c     R,Z bounds of the plasma cell
c
      r1= krb(ikp-1,irp) 
      z1= kzb(ikp-1,irp) 
c
      r2= krb(ikp,irp) 
      z2= kzb(ikp,irp) 
c
c     jdemod - incorporated Geier IPP/02 which calculates the actual distances to 
c              the point closest to the wall. Note that this corresponds to points
c              at the cell edge if the normal distance does not intersect the cell
c              axis.
c
c afmod begin
c IPP change from 02, only just noticed(?) or deliberately omitted(?)
c Geier IPP/02 added sqrt     
c      IF (ippchange) THEN
c        WRITE(0,*) 'CHECK THE IPP MODE IN calc_wall_plasma'
c
        d1 = sqrt((ri-r1)**2 + (zi-z1)**2) 
        d2 = sqrt((ri-r2)**2 + (zi-z2)**2)  
c
c      ELSE  
c        d1 = (ri-r1)**2 + (zi-z1)**2
c        d2 = (ri-r2)**2 + (zi-z2)**2 
c      ENDIF  
c
c
c      d1 = (ri-r1)**2 + (zi-z1)**2 
c      d2 = (ri-r2)**2 + (zi-z2)**2 
c afmod end
c
c     S value - approximate the S value of the location closest to the wall element
c
      s  = ksb(ikp-1,irp) 
     >   + d1/(d1+d2) * (ksb(ikp,irp)-ksb(ikp-1,irp))
c
c     In second half of cell
c
      if (d2.lt.d1) then  
c
c        At end of ring  
c
         if (ikp.eq.nks(irp)) then 
c
            nen = knds(idds(irp,1))
            ten = kteds(idds(irp,1))
            tin = ktids(idds(irp,1))
            sn  = ksmaxs(irp)
c
         else
c
            nen = knbs(ikp+1,irp)
            ten = ktebs(ikp+1,irp)
            tin = ktibs(ikp+1,irp)
            sn  = kss(ikp+1,irp) 
c
         endif

c
c     In first half of cell
c
      else
c
c        At end of ring  
c
         if (ikp.eq.1) then 
c
            nen = knds(idds(irp,2))
            ten = kteds(idds(irp,2))
            tin = ktids(idds(irp,2))
            sn  = 0.0
c
         else
c
            nen = knbs(ikp-1,irp)
            ten = ktebs(ikp-1,irp)
            tin = ktibs(ikp-1,irp)
            sn  = kss(ikp-1,irp) 
c
         endif
c
      endif
c
c     Calculate base ne,te,ti
c 
c     This overwrites the ri,zi values that were in these arrays
c 
c afmod begin
c IPP change from 02, only just noticed(?) or deliberately omitted(?)
c Geier IPP/02 corrected signs and variable names for interpolation
c
c      IF (IPPCHANGE) THEN
c        WRITE(0,*) 'CHECK THIS IPP MOD AS WELL'
c        wallpt(in,29) = ten + (s-sn)/(st-sn) * (te-ten) 
c        wallpt(in,30) = tin + (s-sn)/(st-sn) * (ti-tin) 
c        wallpt(in,31) = nen + (s-sn)/(st-sn) * (ne-nen) 
c      ELSE
c
c      jdemod - original code was wrong - either Geier's code or the revision below fixes the problem
c
        wallpt(in,29) = te + (s-st)/(sn-st) * (ten-te) 
        wallpt(in,30) = ti + (s-st)/(sn-st) * (tin-ti) 
        wallpt(in,31) = ne + (s-st)/(sn-st) * (nen-ne) 
c
c      ENDIF
c
c      wallpt(in,29) = te + (s-sn)/(st-sn) * (ten-te) 
c      wallpt(in,30) = ti + (s-sn)/(st-sn) * (tin-ti) 
c      wallpt(in,31) = ne + (s-sn)/(st-sn) * (nen-ne) 
c afmod end
c
c     Apply linear or exponential decay to these values
c      
      if (wall_plasma_opt.eq.1) then 

         if (wall_plasma_fact.eq.0.0) then 
            fact = 1.0
         else
            fact = wallpt(in,28)/wall_plasma_fact
         endif
c
c        Check valid value of fact 
c           
         if (fact.ge.1.0) then 
c
            wallpt(in,29) = te_min
            wallpt(in,30) = ti_min
            wallpt(in,31) = ne_min
c
         else
c afmod begin
c Geier IPP/02 corrected indices -bug
c           IF (IPPCHANGE) THEN
c
c           jdemod - leave only the index corrected version
c
             wallpt(in,29) = wallpt(in,29) 
     >  		   + fact * (te_min-wallpt(in,29))
	     wallpt(in,30) = wallpt(in,30) 
     >  		   + fact * (ti_min-wallpt(in,30))
	     wallpt(in,31) = wallpt(in,31) 
     >  		   + fact * (ne_min-wallpt(in,31))
c
c           ELSE
c             wallpt(in,29) = wallpt(in,29) 
c     >                     + fact * (te_min-wallpt(in,29))
c             wallpt(in,29) = wallpt(in,30) 
c     >                     + fact * (ti_min-wallpt(in,30))
c             wallpt(in,29) = wallpt(in,31) 
c     >                     + fact * (ne_min-wallpt(in,31))
c           ENDIF
c
c            wallpt(in,29) = wallpt(in,29) 
c     >                    + fact * (te_min-wallpt(in,29))
c            wallpt(in,29) = wallpt(in,30) 
c     >                    + fact * (ti_min-wallpt(in,30))
c            wallpt(in,29) = wallpt(in,31) 
c     >                    + fact * (ne_min-wallpt(in,31))
c afmod end
c
         endif 
c
      elseif (wall_plasma_opt.eq.2) then 
c
         fact = exp(-wallpt(in,28)/wall_plasma_fact)
c
         wallpt(in,29) = max(wallpt(in,29) * fact,te_min)
         wallpt(in,29) = max(wallpt(in,30) * fact,ti_min)
         wallpt(in,29) = max(wallpt(in,31) * fact,ne_min)
c
      endif
c
      return 
      end  
c
c
c
      subroutine setup_fp
      use mod_fp_transport ! fp transport module
      implicit none
      include 'params'
      include 'cgeom'
      include 'fperiph_com'
      include 'driftvel'
      include 'comtor'
      include 'hc_global_opts'
c
c     SETUP_FP::
c
c     This routine sets up the far periphery plasma in specific 
c     far periphery plasma arrays - however - it also assigns these
c     background plasma values to the irwall and irtrap rings
c     respectively - depending on the fp_region that is being referenced.
c     The reason for this is so that the collisionality and other transport
c     coefficients will then be correctly calculated. 
c 
c
      integer ik,ir,in,iv
      real mindist
      integer fp_irmain,fp_irpfz
c
      integer id
c
c     Code default rings for fp geometry and possible plasma reference
c
      fp_irmain = irwall-1
      fp_irpfz   = irtrap+1
      
      fp_rings(fp_main) = fp_irmain
      fp_cells(fp_main) = nks(fp_irmain)
      fp_rings(fp_pfz)  = fp_irpfz
      fp_cells(fp_pfz)  = nks(fp_irpfz)
c
c     Loop over fp regions and set them up
c
c     Initialize plasma to 0.0 (F90 assigns the value to the entire array)
c
      fp_plasma = 0.0
c
      do in = 1,num_fp_regions
 
c        Assign plasma and geometry data for fp region
c     
         call assign_fp_data(fp_rings(in),in)
c
c        Assign wall data for fp region
c
         call assign_fp_wall(in)
c
      end do
c
c
c      write(6,*) 'WALL FOR FP:'
c      do in = 1,wallpts
c
c        Assign end points of wall segment
c
c         write(6,'(10g18.8)') 
c     >     wallpt(in,20),
c     >     wallpt(in,21)
c         write(6,'(10g18.8)')
c     >     wallpt(in,22),
c     >     wallpt(in,23)
c         write(6,'(10g18.8)')
c      end do
c
c
c     Write out the wall distances for debugging:
c
c     Main
c
      if (cprint.eq.3.or.cprint.ge.9) then

         do ik = 0,2*fp_cells(fp_main)
c
            write(6,'(a,3i6,3g18.8)') 'FP_WALLDIST: MAIN:',fp_irmain,ik,
     >         ik/2+1, fp_walldist(ik,fp_main),
     >         fp_wallcoords(ik,fp_main,1),fp_wallcoords(ik,fp_main,2)
         enddo
c
c        Pfz
c
         do ik = 1,fp_cells(fp_pfz)
c
            write(6,'(a,3i6,3g18.8)') 'FP_WALLDIST: PFZ:',fp_irpfz,ik,
     >         ik/2+1, fp_walldist(ik,fp_pfz),
     >         fp_wallcoords(ik,fp_pfz,1),fp_wallcoords(ik,fp_pfz,2)
c
         enddo

c
c        Print plasma conditions
c
         do in = 1,num_fp_regions
             do ik =1, fp_cells(in)
               write(6,'(a,2i6,7g12.5)') 'FP PLASMA:',in,ik,
     >             (fp_plasma(ik,in,id),id=1,7)
            end do
         end do

      endif
c
c     Set the global far periphery transport options
c     NOTE: CTEMAV does not have a real value at this point - it requires that the 
c           neutral launch routine be run - CTEMAV is the average temperature at
c           ionization - if CTEMAV is not properly assigned a value - a value of 2.0eV
c           is used.
c
c
c      
      call fp_set_options(cioptb,cioptc,cioptd,cioptn,cioptm,cpdrft,
     >                    pinchopt,crmb,crmi,rizb,
     >                    global_hc_follow_option,cdperpfp,
     >                    cstgrad,
     >                    czo,chzo,czenh,cizeff,ctemav,irspec,qtim,
     >                    debug_fp)
c
      return
      end
c
c
c
      subroutine assign_fp_wall(ireg)
      implicit none
      integer ireg
c
      include 'params'
      include 'cgeom'
      include 'walls_com'
      include 'fperiph_com'
c
c     ASSIGN_FP_WALL:
c
c
c     Loop over fp wall regions and determine the closest
c     distance to the walls along with the coordinates of the
c     wall intersections from the middle of each cell side as 
c     well as the cell corners when needed. 
c     Also record the R,Z coordinates of this intersection. 
c
c
      real*8 scale_len
      integer ir,side,ik,fp_vertex,iv
      real mindist
      integer ipolya,ipoly,ipolyb


c
c     Loop over the far periphery rings
c
      scale_len = (rmax-rmin) + (zmax-zmin)
c
      ir = fp_rings(ireg)
c
c     Define the polygon side for which distance data is sought
c     - this will need to be modified if additional fp regions are added
c
      if (ireg.eq.fp_main) then 
         side = OUTWARD23
      elseif (ireg.eq.fp_pfz) then 
         side = INWARD41
      endif
c
c     Need to loop through fp_vertex = 0 for all ik
c                          fp_vertex = 1 for ik = 1, nks(ir)-1
c     Set fp_vertex=-1,ik=1 and fp_vertex=1,ik=nks(ir)    
c
c
      do ik = 1,nks(ir)
c
c
c        Get polygon indices
c
         if (ik.eq.1) then
            ipolya = 0
            ipoly  = korpg(ik,ir)
            ipolyb = korpg(ik+1,ir)
         elseif (ik.eq.nks(ir)) then 
            ipolya = korpg(ik-1,ir)
            ipoly  = korpg(ik,ir)
            ipolyb = 0
         else
            ipolya = korpg(ik-1,ir)
            ipoly  = korpg(ik,ir)
            ipolyb = korpg(ik+1,ir)
         endif

c
c         write(0,'(a,6i6)') 'Processing Cell:',ik,ir,ireg,
c     >                           ipolya,ipoly,ipolyb
c
c        Loop over vertices - leave out -1 so as not to repeat calculations unnecessarily
c
         do fp_vertex = 0,1
c
c           Do not process cell nks(ir) with fp_vertex=1
c
            if (ik.eq.nks(ir).and.fp_vertex.eq.1) cycle
c
c           Check for valid polygon data 
c
            if (ipoly.eq.0.or.
     >         (ipolya.eq.0.and.fp_vertex.eq.-1).or.
     >         (ipolyb.eq.0.and.fp_vertex.eq.1)) then 
c
c              Assign default values - input for distance and ring coordinates for
c              intersections
c
               if ((xpoint_up.and.ik.le.nks(ir)/2).or.
     >            (.not.xpoint_up.and.ik.gt.nks(ir)/2)) then 
c
                  fp_walldist(2*ik-1+fp_vertex,ireg) = fpxmaxo
c
               else 
c
                  fp_walldist(2*ik-1+fp_vertex,ireg) = fpxmaxi   
c
               endif
c
               fp_wallcoords(2*ik-1+fp_vertex,ireg,1) = 
     >                                            rs(ik,fp_rings(ireg))
               fp_wallcoords(2*ik-1+fp_vertex,ireg,2) = 
     >                                            zs(ik,fp_rings(ireg))

            else
c
c               Polygon data exists - look to assign wall data to fp cell
c
                call calc_fp_wall_data(ik,ireg,fp_vertex,scale_len,
     >                               ipolya,ipoly,ipolyb,side)
c
            endif

         end do 

      end do
c
c     Assign data for the first and last indices on the ring
c
      if (xpoint_up) then 
c
         fp_walldist(0,ireg) = fpxmaxo
         fp_walldist(2*nks(ir),ireg) = fpxmaxi
c
      else
c
         fp_walldist(0,ireg) = fpxmaxi
         fp_walldist(2*nks(ir),ireg) = fpxmaxo
c
      endif
c
c     May need to change these coordinates since they are not ideal for cross
c     field losses in the first half of the first cells on the ring. 
c
      fp_wallcoords(0,ireg,1) = rp(idds(ir,2))
      fp_wallcoords(0,ireg,2) = zp(idds(ir,2))
c
      fp_wallcoords(2*nks(ir),ireg,1) = rp(idds(ir,1))
      fp_wallcoords(2*nks(ir),ireg,2) = zp(idds(ir,1))
c
c     Loop through and assign the min_fp_walldist for each cell in the FP
c     Checks for wall collisions will only occur when Cross exceeds this
c     value
c
      do ik = 1,nks(ir)
         mindist = HI
         do iv = -1,1
            mindist = min(mindist,fp_walldist(2*ik-1+iv,ireg))
         end do
         min_fp_walldist(ik,ireg) = mindist
      end do
c
      return
      end
c
c
c
      subroutine calc_fp_wall_data(ik,ireg,fp_vertex,scale_len,
     >                              ipolya,ipoly,ipolyb,side)
      implicit none
      integer ik,ireg,fp_vertex,ipolya,ipoly,ipolyb,side
      real*8 scale_len
c
      include 'params'
      include 'cgeom'
      include 'fperiph_com'
c
c     Local variables
c
      integer side2,ir
      real*8 theta
      real*8 ra,za,rb,zb,rsect,zsect,wdist
      logical sect_found
c
      real*8 dra,dza,rthet,datan2c
      external datan2c

c
      real*8 fp_theta,angle_average,theta1,theta2
      external fp_theta,angle_average
c
      ir = fp_rings(ireg)
      side2 = mod(side,4)+1 
c
      if (fp_vertex.eq.-1) then 

         theta1 = fp_theta(ipoly,side) 
         theta2 = fp_theta(ipolya,side)

         theta = angle_average(theta1,theta2)

         ra =  rvertp(side,ipoly)
         za =  zvertp(side,ipoly)

      elseif (fp_vertex.eq.0) then 
c
         theta = fp_theta(ipoly,side)
c
         ra = (rvertp(side2,ipoly) + rvertp(side,ipoly))/2.0
         za = (zvertp(side2,ipoly) + zvertp(side,ipoly))/2.0
c

      elseif (fp_vertex.eq.1) then 

         theta1 = fp_theta(ipoly,side) 
         theta2 = fp_theta(ipolyb,side)

         theta = angle_average(theta1,theta2)

         ra = rvertp(side2,ipoly)
         za = zvertp(side2,ipoly)

      endif
c
c     Assign other end of test line 
c
      rb = ra + scale_len * cos(theta)
      zb = za + scale_len * sin(theta)

c
c     Assign default values
c
      if ((xpoint_up.and.ik.le.nks(ir)/2).or.
     >    (.not.xpoint_up.and.ik.gt.nks(ir)/2)) then 
c
         fp_walldist(2*ik-1+fp_vertex,ireg) = fpxmaxo
c
      else 
c
         fp_walldist(2*ik-1+fp_vertex,ireg) = fpxmaxi   
c
      endif
c
      fp_wallcoords(2*ik-1+fp_vertex,ireg,1) = ra
      fp_wallcoords(2*ik-1+fp_vertex,ireg,2) = za
c
c     Find closest wall intersection outside the grid. 
c
      call find_wall_intsect(ra,za,rb,zb,rsect,zsect,wdist,
     >                       sect_found)
c
c     Assign results if intersection found
c 
      if (sect_found) then

         fp_walldist(2*ik-1+fp_vertex,ireg) = wdist
         fp_wallcoords(2*ik-1+fp_vertex,ireg,1) = rsect
         fp_wallcoords(2*ik-1+fp_vertex,ireg,2) = zsect

      endif

      dra = rsect-ra
      dza = zsect-za
      rthet = datan2c(dza,dra)

c
c      write(0,'(a,5i6,l4,10g18.8)') 'FP WALL:',2*ik-1+fp_vertex,
c     >   side2,side,
c     >   ir,ireg,sect_found,
c     >   theta*raddeg,ra,za,rb,zb,rsect,zsect,wdist
c
c      write(6,'(a,7i6,l4,12g18.8)') 'FP WALL:',2*ik-1+fp_vertex,
c     >   ik,fp_vertex,side2,side,
c     >   ir,ireg,sect_found,
c     >   theta*raddeg,rthet*raddeg,ra,za,rb,zb,rsect,zsect,wdist
c
c      write(64,'(6g18.8)') rsect,zsect
c
c      write(65,'(6g18.8)') ra,za
c      write(65,'(6g18.8)') rb,zb
c      write(65,'(6g18.8)') 
c
c      write(66,'(6g18.8)') rvertp(side,ipoly),zvertp(side,ipoly)
c      write(66,'(6g18.8)') rvertp(side2,ipoly),zvertp(side2,ipoly)
c      write(66,'(6g18.8)') 
c
c
c
      return
      end
c
c
c
      real*8 function angle_average(theta1,theta2)
      implicit none
      real*8 theta1,theta2
      include 'params'
c
c     Averages two angles to find the bisecting angle
c
      real*8 t1,t2,ta,tb
c
      if (theta1.lt.0.0) then 
         ta = theta1 + 2.0d0*PI
      else
         ta = theta1
      endif
c
      if (theta2.lt.0.0) then 
         tb = theta2 + 2.0d0*PI
      else
         tb = theta2
      endif
c
c     Both angles in 0 to 2PI
c
      t1 = max(ta,tb)
      t2 = min(ta,tb)
c
      if (t2.gt.t1-PI) then 
         angle_average=(t1+t2)/2.0d0
      else
         angle_average=(t1-2.0*PI+t2)/2.0d0
      endif
c
c      write(6,'(a,10g18.8)') 'ANGAV:',theta1*raddeg,theta2*raddeg,
c     >                   t1*raddeg,t2*raddeg,angle_average*raddeg
c

      return
      end
c
c
c
      real*8 function fp_theta(ipoly,side)
      implicit none
      integer ipoly,side
c
c     Returns the theta value between 0.0 and 2PI 
c     
      include 'params'
      include 'cgeom'
c
      integer side2
      real*8 datan2c,deltar,deltaz
      external datan2c
c
      side2 = mod(side,4)+1 
c  
      deltar = rvertp(side2,ipoly) - rvertp(side,ipoly)
      deltaz = zvertp(side2,ipoly) - zvertp(side,ipoly)
c
      fp_theta = datan2c(deltaz,deltar) + PI/2.0
c
c      if (fp_theta.lt.0.0) then 
c         fp_theta = fp_theta + 2.0*PI
c      endif
c
c      write(6,'(a,2i6,10g18.8)') 'ANGLE:',side2,side,fp_theta*raddeg,
c     >  deltaz,deltar,rvertp(side2,ipoly),zvertp(side2,ipoly),
c     >  rvertp(side,ipoly),zvertp(side,ipoly)
c
      return
      end
c
c
c
      subroutine find_wall_intsect(ra,za,rb,zb,rsect,zsect,
     >                             wdist,sect_found)
      implicit none
      real*8 ra,za,rb,zb,rsect,zsect,wdist
      logical sect_found
c
      include 'params'
      include 'walls_com'
c
c     Local variables
c
      real*8 r1,z1,r2,z2,rint,zint,dist
      integer in,sect
c
      sect_found = .false.
c
c     Loop through wall elements looking for the closest intersection
c
      wdist=HI
c
      do in = 1,wallpts
c
c        Assign end points of wall segment
c
         r1 = wallpt(in,20)
         z1 = wallpt(in,21)
         r2 = wallpt(in,22)
         z2 = wallpt(in,23)
c
c        call intersection routine
c
         call intsect2dp(ra,za,rb,zb,r1,z1,r2,z2,rint,zint,sect)
c
c        Check for intersection on both segments
c              
         if (sect.eq.1) then 
c
c           Set minimum distance from ra,za to rint,zint
c
            dist = sqrt((ra-rint)**2+(za-zint)**2)
c
            if (dist.lt.wdist) then 
               wdist = dist
               rsect = rint
               zsect = zint
c               write(6,'(a,3g18.8)') 'SECT :',rsect,zsect,dist

            endif
c
            sect_found = .true.
c
         endif
c
      enddo
c
      return
      end 
c
c
c
      subroutine assign_fp_data(ir,in)
      implicit none
      integer ir,in
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'fperiph_com'
c
c     ASSIGN_FP_PLASMA:
c
c     Assign plasma data to the peripheral region
c
      integer ik,id
      real dtebac,dtefor,dtibac,dtifor,dsbac,dsfor
      real fact
c     
c     Assign base plasma data  
c
c     The fp_plasma array contains the following data:
c
c     fp_plasma(ik,in,1) = Ne
c     fp_plasma(ik,in,2) = Te
c     fp_plasma(ik,in,3) = Ti
c     fp_plasma(ik,in,4) = Vb
c     fp_plasma(ik,in,5) = E
c     fp_plasma(ik,in,6) = Te gradient
c     fp_plasma(ik,in,7) = Ti gradient
c
c
      do ik = 1,nks(ir)
c
c        Peripheral plasma assigned from grid
c
         if (fp_plasma_opt.eq.0) then 
c
            fp_plasma(ik,in,1) = knbs(ik,ir)
            fp_plasma(ik,in,2) = ktebs(ik,ir)
            fp_plasma(ik,in,3) = ktibs(ik,ir)
            fp_plasma(ik,in,4) = kvhs(ik,ir)
            fp_plasma(ik,in,5) = kes(ik,ir)
c
         elseif (fp_plasma_opt.eq.1) then 
c
c           Te=Ti=INPUT - rest assigned from ring
c
            fp_plasma(ik,in,1) = knbs(ik,ir)
            fp_plasma(ik,in,2) = fp_te
            fp_plasma(ik,in,3) = fp_te
            fp_plasma(ik,in,4) = kvhs(ik,ir)
            fp_plasma(ik,in,5) = kes(ik,ir)
c
         elseif (fp_plasma_opt.eq.2) then 
c
            fp_plasma(ik,in,1) = fp_ne
            fp_plasma(ik,in,2) = fp_te
            fp_plasma(ik,in,3) = fp_te
            fp_plasma(ik,in,4) = kvhs(ik,ir)
            fp_plasma(ik,in,5) = kes(ik,ir)
c
         endif
c
c        Assign fp_plasma to associated outermost grid rings as well
c
         if (in.eq.fp_main) then 
            knbs(ik,irwall) = fp_plasma(ik,in,1)
            ktebs(ik,irwall)= fp_plasma(ik,in,2)
            ktibs(ik,irwall)= fp_plasma(ik,in,3)
            kvhs(ik,irwall) = fp_plasma(ik,in,4)
            kes(ik,irwall)  = fp_plasma(ik,in,5)
         elseif (in.eq.fp_pfz) then
            knbs(ik,irtrap) = fp_plasma(ik,in,1)
            ktebs(ik,irtrap)= fp_plasma(ik,in,2)
            ktibs(ik,irtrap)= fp_plasma(ik,in,3)
            kvhs(ik,irtrap) = fp_plasma(ik,in,4)
            kes(ik,irtrap)  = fp_plasma(ik,in,5)
         endif
c
      end do 
c
c     Assign cell geometry data (S-values) as well 
c
c      write(6,'(a,3i6,10(1x,g12.5))') 'FP_S: ',ir,in,nks(ir),ksmaxs(ir)

      do ik = 1,nks(ir)
c
         fp_s(2*ik-1-1,in) = ksb(ik-1,ir)
         fp_s(2*ik-1,in)   = kss(ik,ir)
         fp_s(2*ik-1+1,in) = ksb(ik,ir)
c
c         write(6,'(a,1i6,3(1x,g12.5)') 'FP_IK:',ik,
c     >       fp_s(2*ik-1-1,in),fp_s(2*ik-1,in),
c     >       fp_s(2*ik-1+1,in)
c
      end do
c
c     Calculate the FP temperature gradients
c
c     Scaling factor
c
      FACT = QTIM * QTIM * EMI / CRMI
c
c     Calculate gradients for first and last cells 
c
      ik = 1
      dtefor = fp_plasma(ik+1,in,2)-fp_plasma(ik,in,2)
      dtifor = fp_plasma(ik+1,in,3)-fp_plasma(ik,in,3)
      dsfor = fp_s(2*(ik+1)-1,in) - fp_s(2*ik-1,in)
      fp_plasma(ik,in,6) = (dtefor/dsfor) * fact * 0.5
      fp_plasma(ik,in,7) = (dtifor/dsfor) * fact * 0.5     
c
      ik = nks(ir)
      dtebac = fp_plasma(ik,in,2)- fp_plasma(ik-1,in,2)
      dtibac = fp_plasma(ik,in,3)- fp_plasma(ik-1,in,3)
      dsbac = fp_s(2*ik-1,in) - fp_s(2*(ik-1)-1,in)
      fp_plasma(ik,in,6) = (dtebac/dsbac) * fact * 0.5
      fp_plasma(ik,in,7) = (dtibac/dsbac) * fact * 0.5
c
      do ik = 2,nks(ir)-1
c
c        Zero out the entries where the plasma background is
c        discontinuous
c
         if (ik.eq.ikmids(ir).or.ik.eq.ikmids(ir)+1) then
            fp_plasma(ik,in,6) = 0.0
            fp_plasma(ik,in,7) = 0.0
         else
c
c           Calculate average gradients
c
            dtebac = fp_plasma(ik,in,2)- fp_plasma(ik-1,in,2)
            dtefor = fp_plasma(ik+1,in,2)-fp_plasma(ik,in,2)

            dtibac = fp_plasma(ik,in,3)- fp_plasma(ik-1,in,3)
            dtifor = fp_plasma(ik+1,in,3)-fp_plasma(ik,in,3)

            dsbac = fp_s(2*ik-1,in) - fp_s(2*(ik-1)-1,in)
            dsfor = fp_s(2*(ik+1)-1,in) - fp_s(2*ik-1,in)

            fp_plasma(ik,in,6)=(dtebac/dsbac + dtefor/dsfor)* fact * 0.5
            fp_plasma(ik,in,7)=(dtibac/dsbac + dtifor/dsfor)* fact * 0.5
         endif
c
      end do  

      return 
      end

c
c
c
      real function get_radial_sepdist(ik,ir)
      implicit none
      integer ik,ir
      include 'params'
      include 'cgeom'
c
      real get_refdist
      external get_refdist
c
c     
c
      integer rc
c
c     Approximate radial distance to the separatrix from cell ik,ir
c
c     returns 0.0 for cells not connected through to the LCFS
c
      if (irins(ik,irsep).eq.irsep-1) then 
c
c slmod begin
         get_radial_sepdist  = 
     >          sqrt( (rs(ik,ir)-rs(ik,irsep))**2 + 
     >                (zs(ik,ir)-zs(ik,irsep))**2 ) +
c
c         get_radial_sepdist  = 
c     >          sqrt( (rs(ik,ir)-rs(ik,irsep))**2 + 
c     >              + (zs(ik,ir)-zs(ik,irsep))**2 ) +
c slmod end
     >          get_refdist(ik,irsep,INWARD41,rc) 
c
      else
c 
         get_radial_sepdist = 0.0
c
      endif
c
      return
      end
c
c     
c
      real function weighted_sepdist(ik,ir,len)
      implicit none
      include 'params'
      include 'cgeom'
c
      integer ik,ir
      real len,dist
      real get_radial_sepdist
      external get_radial_sepdist 
c
      weighted_sepdist = 0.0
c
c     Poloidal length  
c
      len = kpb(ik,ir) - kpb(ik-1,ir)
c 
      dist = get_radial_sepdist(ik,ir)
c
      if (dist.eq.0) then 
         len = 0.0
      else
         weighted_sepdist = dist * len
      endif
c
      return 
      end          
c
c
c
      subroutine print_average_sepdist
      implicit none
      include 'params'
      include 'cgeom'
c
      integer ik,ir
      real totlen,totdist
      real weighted_sepdist
      external weighted_sepdist
      real len
c
c     Loop around sol rings
c
      do ir = irsep,irwall-1 
c
         totlen = 0.0
         totdist = 0.0
c
         do ik = 1,nks(ir)
c
            totdist = totdist 
     >              + weighted_sepdist(ik,ir,len)
            totlen  = totlen+len 
c
         end do
c
         write(6,'(a,3x,i6,3x,g12.5)') ' AVERAGE SEPDIST IR:',ir,
     >         totdist/totlen
c
      end do
c
      return
      end 
c
c
c
      subroutine set_bcomponents(br,bz,bt)
      implicit none
      include 'params'
      include 'cgeom'
c
      real br(maxnks,maxnrs)
      real bz(maxnks,maxnrs)
      real bt(maxnks,maxnrs)
c
c     This routine calculates the direction of the magnetic field vector
c     for each cell on the grid - assuming a cylindi=rical geometry. 
c
c     This code does not calculate the absolute magnitude of the B field only the 
c     normalized direction vector.
c
c     First need to calculate the cell axis vector for each cell on the 
c     grid. From these - caclulate the BT vector that scales with the 
c     other two and the bratio for the cell. Finally, renormalize the 
c     magnitude of the components so that the bfield vector has a magnitude
c     of 1. 
c
c     Start the bfield calculations at the ik=1 target for simplicity
c       
      integer ik,ir
      real btot
      real btmp
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            br(ik,ir) = krb(ik,ir)-krb(ik-1,ir)
            bz(ik,ir) = kzb(ik,ir)-kzb(ik-1,ir)
c
c           Temporarily set bt to contain Btot - which is 
c           Bpol * kbfs
c
            btot = kbfs(ik,ir) * sqrt(br(ik,ir)**2 + bz(ik,ir)**2)
c
c           Now use the bratio to calculate btotal and thus btoroidal
c
c           bratio = Bpol/B
c           kbfs   = 1/bratio = B/Bpol
c     
c           B = kbfs * Bpol
c
c           Bt = sqrt(B**2 - Br**2 - Bz**2) 
c     
            btmp = btot**2 - br(ik,ir)**2 -bz(ik,ir)**2
            if (btmp.ge.0.0) then 
               bt(ik,ir) = sqrt(btmp)
            else
               write(6,'(a,2i6,10(1x,g12.5))') 
     >              'BFIELD CALCULATION ERROR:',
     >              ik,ir,btot,br(ik,ir),bz(ik,ir),kbfs(ik,ir),btmp
               bt(ik,ir) = 0.0
            endif
c
c           Normalize the b-vector to magnitude 1
c
            if (btot.gt.0.0) then 
               br(ik,ir) = br(ik,ir)/btot
               bz(ik,ir) = bz(ik,ir)/btot
               bt(ik,ir) = bt(ik,ir)/btot
            endif
c
            write(6,'(a,2i6,10(1x,g12.5))') 'BFIELD:',ik,ir,
     >          kbfs(ik,ir),br(ik,ir),bz(ik,ir),bt(ik,ir),btot,
     >        sqrt(br(ik,ir)**2+bz(ik,ir)**2+bt(ik,ir)**2)

         end do
      end do
c
         

      return
      end
c
c
c
      subroutine calculate_pinch
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     local variables 
c
      integer :: startir,endir
      integer :: ik,ir
c
      real cell_dr,cell_dz,cell_norm
      logical abovexp,oninside,abovezpinch,abovexpsep

      real zpinch

      zpinch = -1.05


C-----------------------------------------------------------------------
c
c     Scale the Perpendicular Pinch Velocity
c
      vpinch = cvpinch * qtim
c
c     Initialization for spacially varying pinches
c
      kpinchs = 0.0
      kpinchs_para = 0.0
c
c slmod begin
c
c     Krieger IPP/97
c
c     added quadratic pinch term kpinchs; Krieger IPP/97
c     IMPORTANT: DIVIMP's sign convention is contrary to the usual
c     convention i.e. negative velocity means outward. To keep the
c     usual convention for the parameter file, the sign is changed
c     here.
c
c     Calculate pinch velocity as a function of grid cell
c
      if (pinchopt.eq.3) then  
c
        DO IR = 1, NRS
          DO IK = 1, NKS(IR)
            ! CVPINCH needs to be scaled by qtim as shown here
            if (ir.lt.irsep) then
               kpinchs(ik,ir) =  -cvpinch * qtim *
     >                      (kpmaxs(ir)/kpmaxs(irsep-1))**2
            else
              kpinchs(ik,ir) = 0.0
            endif
          end do
        end do
c
c     Setup values in the kpinchs array for option 6 - the sign
c     of the pinch changes depending on the R location of the cell
c     center - thus the pinch is generally all towards the center 
c     column or all towards the outer midplane depending on the sign
c     of cvpinch - the pinch is also applied inside the core plasma
c     for this option.
c
      elseif (pinchopt.eq.6.or.pinchopt.eq.7) then 
c
        if (pinchopt.eq.6) then 
           startir = 1
           endir = nrs
        elseif (pinchopt.eq.7) then
           startir = irsep
           endir = nrs
        endif
c
        DO IR = startir, endir
          DO IK = 1, NKS(IR)
            !
            ! vpinch is already scaled by QTIM (above) so multiplication is 
            ! not required. 
            ! Positive vpinch in the input file should be a net outward 
            ! pinch while negative will be net inward. 
            !
            if (rs(ik,ir).lt.rxp) then
               kpinchs(ik,ir) =   vpinch 
            else
               kpinchs(ik,ir) =   -vpinch
            endif
          end do
        end do
c
      elseif (pinchopt.eq.8.or.pinchopt.eq.9.or.pinchopt.eq.10.or.
     >        pinchopt.eq.11.or.pinchopt.eq.12.or.pinchopt.eq.13.or.
     >        pinchopt.eq.14.or.pinchopt.eq.15) then 
c
        if (pinchopt.eq.8) then 
           startir = 1
           endir = nrs
        elseif (pinchopt.eq.9.or.pinchopt.eq.10.or.pinchopt.eq.11.or.
     >          pinchopt.eq.12.or.pinchopt.eq.13.or.pinchopt.eq.14)then 
           startir = irsep
           endir = nrs
        elseif (pinchopt.eq.15) then 
           startir = irsep
           endir = irwall
        endif
c
        DO IR = startir, endir
          DO IK = 1, NKS(IR)
             !
             ! The pinch is assumed to be purely radial thus in 
             ! vector form is (Vr,0)
             !
             ! Projection of this parallel to S is Vr* norm_dr
             ! Projection of this perpendicular to the cell axis (cross)
             !   is Vr*norm_dz - the normal vector is (norm_dz,-norm_dr) - 
             ! Vr dot perp_vector gives only the first component since Vz=0
             !
             ! This should properly change the signs for the grid geometry
             !


             abovezpinch =(zs(ik,ir).gt.zpinch.and.(.not.xpoint_up)).or.
     >                 (zs(ik,ir).lt.zpinch.and.xpoint_up)  
             
             abovexp = (zs(ik,ir).gt.zxp.and.(.not.xpoint_up)).or.
     >                 (zs(ik,ir).lt.zxp.and.xpoint_up)  
             oninside = (rs(ik,ir).lt.r0)

             if (ir.ge.irsep.and.ir.le.irwall) then 
               abovexpsep=(zs(ik,irsep).gt.zxp.and.(.not.xpoint_up)).or.
     >                    (zs(ik,irsep).lt.zxp.and.xpoint_up)  
             else
                abovexpsep = .false.
             endif

             
             if (pinchopt.eq.9.or.
     >          (pinchopt.eq.10.and.abovexp).or.
     >          (pinchopt.eq.12.and.abovexp.and.oninside).or.
     >          (pinchopt.eq.14.and.abovezpinch).or.
     >          (pinchopt.eq.15.and.abovexpsep)) then 

            !
            ! vpinch is already scaled by QTIM (above) so multiplication is 
            ! not required. 
            ! Positive vpinch in the input file should be a net outward 
            ! pinch while negative will be net inward. 
            !

                cell_dr = krb(ik,ir) - krb(ik-1,ir)
                cell_dz = kzb(ik,ir) - kzb(ik-1,ir)

                cell_norm = sqrt(cell_dr**2 + cell_dz**2)

                cell_dr = cell_dr/cell_norm
                cell_dz = cell_dz/cell_norm

                kpinchs_para(ik,ir) = vpinch * cell_dr * kbfs(ik,ir)
                kpinchs(ik,ir)      = vpinch * cell_dz

                write(6,'(a,2i6,3(1x,g12.5))') 'RADIAL PINCH:',ik,ir,
     >              kpinchs(ik,ir),kpinchs_para(ik,ir)


             elseif ((pinchopt.eq.11.and.abovexp).or.
     >               (pinchopt.eq.13.and.abovexp.and.oninside)) then

                kpinchs(ik,ir) = vpinch
                kpinchs_para (ik,ir) = 0.0

                write(6,'(a,2i6,3(1x,g12.5))') 'CROSS PINCH:',ik,ir,
     >              kpinchs(ik,ir),kpinchs_para(ik,ir)
             else
                kpinchs_para(ik,ir) = 0.0
                kpinchs(ik,ir)      = 0.0
             endif

          end do
        end do
c
      endif
c
c      write(6,'(a)') 'Ring   Pinch velocity'
c      do ir = 1,irsep-1
c        write(6,'(i2.2,4x,1p,e10.3)') ir,-kpinchs(1,ir)/qtim
c      enddo
c
      return
      end

      real function acos_test(xcos,flag)
      use error_handling
      implicit none
      real xcos
      integer flag

c slmod begin
c...  Adding a bit of precision tolerance:
c      if (xcos.lt.-1.000005) then 
      if (xcos.lt.-0.9999995) then 
c
c      if (xcos.lt.-1.0) then 
c slmod end
c         write(0,*) 'ACOS ERROR: XCOS < -1.0: FLAG = ',
c     >                    flag,xcos,acos(-1.0)
c         write(6,*) 'ACOS ERROR: XCOS < -1.0: FLAG = ',
c     >                    flag,xcos,acos(-1.0)
c
         write(error_message_data,*) 'ACOS ERROR: XCOS < -1.0: FLAG = ',
     >                    flag,xcos,acos(-1.0)
         
         call dbgmsg('ACOS_TEST:',error_message_data)


         acos_test = acos(-1.0) 

c slmod begin
c...  Adding a bit of precision tolerance:
c      elseif (xcos.gt.1.000005) then 
      elseif (xcos.gt.0.9999995) then 
c
c      elseif (xcos.gt.1.0) then 
c slmod end
c
c          write(0,*) 'ACOS ERROR: XCOS > 1.0: FLAG = ',
c     >                      flag,xcos,acos(1.0)
c          write(6,*) 'ACOS ERROR: XCOS > 1.0: FLAG = ',
c     >                      flag,xcos,acos(1.0)
c
         write(error_message_data,*) 'ACOS ERROR: XCOS >  1.0: FLAG = ',
     >                    flag,xcos,acos(1.0)
         call dbgmsg('ACOS_TEST:',error_message_data)
c
         acos_test = acos(1.0)

      else
         acos_test=acos(xcos)
      endif

      return
      end
