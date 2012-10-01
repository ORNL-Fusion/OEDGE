c     -*-Fortran-*- 
c     @PROCESS NOOPT
      SUBROUTINE READIN (TITLE,desc,equil,NIZS,NIMPS,NIMPS2,CPULIM,
     >                   IERR,NYMFS,NITERS)
c      SUBROUTINE READIN (TITLE,NIZS,NIMPS,NIMPS2,CPULIM,
c     >                   IERR,NYMFS,NITERS)
      use error_handling
      implicit none
c
      INTEGER   IERR,NIZS,NIMPS,NYMFS,NITERS,NIMPS2
      REAL      CPULIM
      CHARACTER TITLE*(*),equil*(*),desc*(*)
C
C  *********************************************************************
C  *                                                                   *
C  *  READIN:   READS IN THE DATAFILE, PERFORMS VALIDITY CHECKS, ETC.  *
C  *                                                                   *
C  *                                                                   *
C  *  CHRIS FARRELL    FEBRUARY 1989                                   *
C  *                                                                   *
C  *  JAMES SPENCE     NOVEMBER 1990                                   *
C  *                                                                   *
C  *********************************************************************
C
      include    'params'
      include    'cgeom'
      include    'comtor'
      include    'cadas'
      include    'dynam4'
      include    'dynam5'
      include    'diagvel'
      include    'cedge2d'
      include    'pindata'
      include    'adpak_com'
      include    'promptdep'
      include    'reiser_com'
c
      include    'fperiph_com'
      include    'driftvel'
c
      include    'sol23_input'
c slmod begin - new
      INCLUDE    'slcom'
c slmod end
c
c
      INTEGER NQS,ISTEP,in
      REAL    ZO
c
c     jdemod - option for converting LP data input from particles/s to A/s
c     
      real lpdat_conv
c
c     Option indicating if the SOL23 parameter list is included in the data
c     file.
c
      integer readin_sol23_params
c
c slmod begin
      CALL InitializeVariables
c slmod end
      CALL RDC (TITLE, 'TITLE FOR RUN', IERR)
      CALL RDC (desc, 'DESCRIPTION OF RUN', IERR)
      call rdc (equil,'equilibrium file name',ierr)
c slmod begin
c...        Added cgridopt=6 for linear device grids:
c...  jde - added cgridopt=7 for generalized grids
c...  jde - added cgridopt=8 for ribbon grids
      call rdi (cgridopt,.true.,0,.true.,9,'GRID OPTION          ',ierr)
c
c      call rdi (cgridopt,.true.,0,.true.,3,'GRID OPTION          ',ierr)
c slmod end
      call rdi (northopt,.true.,0,.true.,3,'NON-ORTHO. GRID OPT  ',ierr)
      call rdi (pdopt,.true.,0,.true.,1,   'PARALLEL DIST (S) OPT',ierr)
      call rdi (cfdopt,.true.,0,.true.,1,'Proper Cross-field dist',ierr)
      call rdi (rzopt,.true.,0,.true.,3,   'CALCULATE ACTUAL R,Z ',ierr)
      call rdi (xygrid,  .true.,0,.true.,1,'XY GRID OPTION       ',ierr)
      call rdi (cvolopt,.true.,0,.true.,1, 'Cell Volumes from PGs',ierr)
      call rdr (cbphi,.true.,0.0,.false.,0.0,'ON-AXIS B-FIELD',    ierr)
      call rdi (cdatopt,.true.,0,.true.,3, 'Rad/ioniz data source',ierr)
      call rdc (useridh,'ADAS H userid',ierr)
      call rdi (iyearh,.true., 0,.true.,99,'ADAS H year          ',ierr)
      call rdc (useridz,'ADAS Z userid',ierr)
      call rdi (iyearz,.true., 0,.true.,99,'ADAS Z year          ',ierr)
      call rdc (mcfile,'Name of File for ADPAK/INEL atomic data',ierr)
      CALL RDI (CIOPTA,.TRUE., 0,.TRUE., 6,'IONIZATION OPT       ',IERR)
      CALL RDI (CIOPTB,.TRUE., 0,.TRUE.,14,'COLLISION OPT        ',IERR)
      CALL RDI (CIOPTR,.TRUE., 0,.TRUE., 2,'REISER COLLISION OPT ',IERR)
c
c     All other Reiser input is unstructured and optional 
c
c      CALL RDI (ASWITCH,.TRUE., 0,.TRUE., 1,'MAIN COEF SWITCH OPT',IERR)
c      CALL RDI (SK11,.TRUE., 0,.TRUE., 1,'DRIFT K11(v-vb)only OPT',IERR)
c      CALL RDI (SK12,.TRUE., 0,.TRUE., 1,'DRIFT K12(gradTb)   OPT',IERR)
c      CALL RDI (SK13,.TRUE., 0,.TRUE., 1,'DRIFT K13(gradvb)   OPT',IERR)
c      CALL RDI (SD11,.TRUE., 0,.TRUE., 1,'DIFF. D11(v-vb)only OPT',IERR)
c      CALL RDI (SD12,.TRUE., 0,.TRUE., 1,'DIFF. D12(gradTb) OPT  ',IERR)
c      CALL RDI (SD13,.TRUE., 0,.TRUE., 1,'DIFF. D13(gradvb) OPT  ',IERR)
c      Also - coulomb_log and linearpeak
c
      CALL RDI (CIOPTC,.TRUE., 0,.TRUE., 4,'FRICTION OPT         ',IERR)
      CALL RDI (CIOPTD,.TRUE., 0,.TRUE., 3,'HEATING OPT          ',IERR)
c slmod begin
      CALL RDI (CIOPTE,.TRUE., 0,.TRUE.,11,'INJECTION OPT        ',IERR)
c
c      CALL RDI (CIOPTE,.TRUE., 0,.TRUE.,10,'INJECTION OPT        ',IERR)
c slmod end
      CALL RDI (CIOPTF,.TRUE.,-1,.TRUE.,99,'SOL OPT              ',IERR)
c slmod begin
      call rdi (ccoreopt,.true.,-1,.true.,28,'Core Option        ',ierr)
      if (ccoreopt.gt.6.and.ccoreopt.ne.28) then
        write(0,*) 'ERROR READIN: INVALID CORE OPTION'
        ierr=1
        return
      endif
c
c      call rdi (ccoreopt,.true.,-1,.true.,6,'Core Option         ',ierr)
c slmod end
      CALL RDI (CIOPTG,.TRUE., 0,.TRUE.,99,'PLASMA DECAY OPT     ',IERR)
c slmod begin
c...  Need to store these because they are overwritten in BGPLASMA but 
c     are needed in GETMODEL:
      orgcioptf = cioptf
      orgcioptg = cioptg
c slmod end
      CALL RDRARN(BGPLASOPT,NBGPLAS,2*MAXNRS,-MACHHI,MACHHI,.FALSE.,
     >            0.0,MACHHI,7,'SET OF BG PLASMA OPTIONS BY RING',IERR)
      call rdi (cre2d ,.TRUE., 0,.true.,2 ,'READ E2D BG FOR REF  ',ierr)
      call rdi (e2dtargopt,.TRUE., 0,.true.,6 ,'EDGE2D TARG COND',ierr)
      CALL RDI (CIOPTI,.TRUE., 0,.TRUE., 9,'CX RECOMB OPT        ',IERR)
      CALL RDI (CDIFOP,.TRUE., 0,.TRUE., 2,'FIRST DIFFUSE OPT    ',IERR)
      CALL RDI (CIOPTJ,.TRUE., 0,.TRUE., 4,'DPERP OPTION         ',IERR)
      call rdi (cdiffopt,.true.,0,.true.,3,'PERP STEP PROB OPTION',ierr)
      call rdi (pinchopt,.true.,0,.true.,15,
     >                                     'PINCH VELOCITY OPTION',ierr)
      CALL RDI (CIOPTK,.TRUE., 0,.TRUE.,99,'TEB GRADIENT OPTION  ',IERR)
      CALL RDI (CIOPTL,.TRUE., 0,.TRUE.,99,'TIB GRADIENT OPTION  ',IERR)
      CALL RDI (CIOPTM,.TRUE., 0,.TRUE., 3,'TEB GRAD COEFF OPTION',IERR)
      CALL RDI (CIOPTN,.TRUE., 0,.TRUE., 3,'TIB GRAD COEFF OPTION',IERR)
      call rdi (cflatopt,.true.,0,.true.,2,'Flatten Te,i option',ierr)
      call rdr (ctegcut,.FALSE.,0.0,.TRUE.,0.5,'TEB GRAD. CUTOFF',ierr)
      call rdr (ctigcut,.FALSE.,0.0,.TRUE.,0.5,'TIB GRAD. CUTOFF',ierr)
      call rdi (fgradopt,.true.,0,.true.,4,'GRAD FORCE MOD OPTION',ierr)
      CALL RDI (CIOPTO,.TRUE., 0,.TRUE., 4,'TRAP TGRAD OPTION    ',IERR)
c      call rdi (ppforceopt,.TRUE.,0,.TRUE.,1,'TRAP FORCE SWITCH ',IERR)
      CALL RDI (CNEUTA,.TRUE., 0,.TRUE., 1,'CONTROL SWITCH       ',IERR)
      CALL RDI (CNEUTB,.TRUE., 0,.TRUE., 7,'LAUNCH OPTION        ',IERR)
      CALL RDI (CNEUTC,.TRUE., 0,.TRUE.,19,'VEL/ANGLE FLAG       ',IERR)
C
C     IF CNEUTH OR CNEUTI ARE -1 THEY NEED TO BE SET TO THE VALUES
C     READ IN FOR CNEUTB AND CNEUTC.
C
c
C     THESE VARIABLES WERE ADDED TO ALLOW LAUNCHES FROM THE WALLS,
C     IN ADDITION TO THE NORMAL PLATE LAUNCHES. HOWEVER, THE
C     CAPABILITY CAN BE USED TO SIMULATE TWO DIFFERING DISTRIBUTIONS
C     OF PARTICLES COMING FROM THE PLATES. PERHAPS DUE TO
C     PHYSICAL AND CHEMICAL SPUTTERING SIMULATNEOUSLY.
C
      CALL RDI (CNEUTH,.TRUE.,-1,.TRUE., 5,'SUP LAUNCH OPTION    ',IERR)
      CALL RDI (CNEUTI,.TRUE.,-1,.TRUE.,19,'SUP VEL/ANGLE FLAG   ',IERR)
      CALL RDI (NVAOPT,.TRUE.,-1,.TRUE.,19,'NEUT VEL/ANGLE FLAG  ',IERR)
      IF (NVAOPT.EQ.-1) NVAOPT = CNEUTC
      IF (CNEUTH.EQ.-1) CNEUTH = CNEUTB
      IF (CNEUTI.EQ.-1) CNEUTI = CNEUTC
      
c
c     jdemod - remove reference to IPPCHANGE here - not needed
c     
c afmod begin
c  found this IPP/02 change that wasn't in. Am putting it in to investigate.
c This causes non-IPP input files to crash, as this is in unstructured input
c in the current version and not here. Input files  must have diverged at 
c some point...
c So one code is synched, can probabyl ditch this?!
c      write(0,*)'IPPCHANGE:',ippchange
c      IF (ippchange) THEN
c Geier IPP/02 for wall plasma erosion/launch
c
c      CALL RDI (wall_plasma_opt,.TRUE.,-1,.TRUE., 2,
c     >                                      'wall plasma option  ',IERR)
c      CALL RDR (wall_plasma_fact ,.TRUE. ,0.0,.FALSE.,0.0,
c     >                                     'wall_plasma_fact     ',IERR)
c     
c      ENDIF
c afmod end      
c
c     These variables were added so that a 2-D grid source of
c     neutrals/cell could be added to the regular target sources.
c     A proportion of teh total particles to be launched in the
c     case will be diverted to this source proportional to
c     strength of this source mechansim.
c
      CALL RDI (neut2d_opt,.TRUE.,0,.TRUE., 1,
     >                                    'EXTRA 2D LAUNCH OPTION',IERR)
      CALL RDI (neut2d_vaopt,.TRUE.,-1,.TRUE.,19,
     >                                         'EXTRA 2D V/A FLAG',IERR)
c
      CALL RDI (CSPUTOPT,.TRUE., 1,.TRUE., 6,'SPUTTER SOURCE OPT ',IERR)
      CALL RDI (CCHEMOPT,.TRUE., 1,.TRUE.,11,'CHEMSPUT SOURCE OPT',IERR)
      CALL RDI (CNEUTD,.TRUE., 0,.TRUE., 8,'SPUTTER OPTION       ',IERR)
      CALL RDI (CNEUTD2,.TRUE.,-1,.TRUE.,8,'2ND SPUTTER OPTION   ',IERR)
      if (cneutd2.eq.-1) cneutd2 = cneutd
      CALL RDI (CSELFS,.TRUE., 0,.TRUE., 2,'SELF-SPUTTER OPTION  ',IERR)
      CALL RDI (CNEUTE,.TRUE., 0,.TRUE., 2,'NORMAL OPTION        ',IERR)
      CALL RDI (CNEUTF,.TRUE., 0,.TRUE., 1,'NEUT SPREADING       ',IERR)
      CALL RDI (CNEUTG,.TRUE., 0,.TRUE., 3,'INITIAL ION VELOCITY ',IERR)
      CALL RDI (CIONR ,.TRUE., 0,.TRUE., 2,'ION WALL OPTION      ',IERR)
      CALL RDI (CNEUR ,.TRUE., 0,.TRUE., 7,'NEUTRAL WALL OPTION  ',IERR)
      CALL RDI (CTRAP ,.TRUE., 0,.TRUE., 8,'TRAP WALL OPTION     ',IERR)
      call rdi (redefopt,.true.,0,.true.,1,'VESSEL REDEF. OPT    ',IERR)
      call rdi (cneutvel,.true.,0,.true.,2,'IMP.NEUT.VEL.TYPE.OPT',ierr)
c ammod - added options 3 and 4 to NRFOPT
      CALL RDI (NRFOPT,.TRUE., 0,.TRUE., 4,'NEUTRAL REFLECTION   ',IERR)
c
      CALL RDI (CFOLREC,.TRUE., 0,.TRUE.,1,'FOLLOW REC. NEUTRAL  ',IERR)
      call rdi (mtcopt,.true.,0,.true.,2,'NEUT.MOM.TRAN.COLL OPT ',ierr)
      call rdr (kelighi,.true.,0.0,.false.,0.0,'MTC COEFF 1 Imp-I',ierr)
      call rdr (kelighg,.true.,0.0,.false.,0.0,'MTC COEFF 2 Imp-N',ierr)
      call rdi (prompt_depopt,.true.,0,.true.,1,'ION PROMPT DEPOSITION'
     >                                                            ,ierr)
      CALL RDI (CTARGOPT,.TRUE.,0,.TRUE.,6,'TARGET POSITION OPT  ',IERR)
      call rdi (cmiropt, .true.,0,.true.,4,'TARGET MIRROR OPT    ',IERR)
      CALL RDI (CGEOOPT,.TRUE.,-1,.TRUE.,1,'GEOMETRY OPTION      ',IERR)
      CALL RDI (FPOPT,   .TRUE.,0,.TRUE.,5,'ION PERIPHERY OPT    ',IERR)
      call rdi (fpropt,.true.,  0,.true.,1,'FP RECYCLE OPT       ',ierr)
      CALL RDI (CPDRFT,.TRUE., 0,.TRUE., 3,'POL. DRIFT OPTION    ',IERR)
C
      CALL RDI (IRSPEC,.FALSE., 0 ,.FALSE., 0 ,'SPEC PLASMA RSPEC',IERR)
      IF (CIOPTB.EQ.0.AND.(CIOPTC.EQ.0.or.cioptc.eq.4).AND.
     >    (CIOPTD.EQ.0.OR.CIOPTD.EQ.3)) IRSPEC=2*MAXNRS
      call rdi (ircore,.true.,1,.false.,0,'CORE MIRROR RING SPEC' ,ierr)
c
      CALL RDR (CRMB, .TRUE. ,0.1 ,.FALSE.,0.0,'PLASMA ION MASS  ',IERR)
      CALL RDI (CIZB, .TRUE. , 1  ,.FALSE., 0 ,'PLASMA ION CHARGE',IERR)
      RIZB = REAL (CIZB)
C
      CALL RDR(CTEB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TEB0  ',IERR)
      CALL RDR(CTEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TEBP  ',IERR)
      CALL RDR(CTEBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBOUT',IERR)
      CALL RDR(CTEBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBIN ',IERR)
      CALL RDR(CTEBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TEBT  ',IERR)
      CALL RDR(CTEBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TEBOUTP',IERR)
      CALL RDR(CFEBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL1 ',IERR)
      CALL RDR(CFEBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL2 ',IERR)
      CALL RDR(CFEBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBT  ',IERR)
      CALL RDR(CFEB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEB2  ',IERR)
      CALL RDR(CTIB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TIB0  ',IERR)
      CALL RDR(CTIBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TIBP  ',IERR)
      CALL RDR(CTIBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBOUT',IERR)
      CALL RDR(CTIBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBIN ',IERR)
      CALL RDR(CTIBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TIBT  ',IERR)
      CALL RDR(CTIBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TIBOUTP',IERR)
      CALL RDR(CFIBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL1 ',IERR)
      CALL RDR(CFIBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL2 ',IERR)
      CALL RDR(CFIBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBT  ',IERR)
      CALL RDR(CFIB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIB2  ',IERR)
      CALL RDR(CNB0  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS AT 0   NB0   ',IERR)
      CALL RDR(CNEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES DENS NEBP  ',IERR)
      CALL RDR(CNBOUT,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBOUT ',IERR)
      CALL RDR(CNBIN ,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBIN  ',IERR)
      CALL RDR(CNBT  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS TRAP   NBT   ',IERR)
      CALL RDR(CNBOUP,.TRUE. ,0.0,.FALSE.,0.0,'DENS DECAY  CNBOUP',IERR)
C
C     ENTER LANGMUIR PROBE DATA ALONG PLATES FOR TEMPERATURE GRADIENT
C     OPTIONS 3 AND 4, AND PLASMA DECAY OPTIONS 2 AND 3. IN THE
C     CASE OF THE PLASMA DECAY OPTIONS THE DATA ENTERED HERE WILL
C     BE FILLED IN FOR ALL POINTS OF THE PLASMA. THUS GIVING A
C     UNIFORM CONDITION UNLESS A TEMPERATURE GRADIENT OPTION
C     MODIFIES THE TEMPERATURE.
C
C     FOR THOSE OPTIONS WITH ONLY ONE SET OF DATA THE INNER
C     SPECIFICATION IS USED.
C
C     NOTE LPDATI :- LANGMUIR PROBE DATA INNER ...
c
c     THE LANGMUIR PROBE DATA SWITCH IDENTIFIES THE THIRD QUANTITY
c     AS EITHER THE DENSITY OR THE Isat (PROBE SATURATION CURRENT)
c
C     THE DATA IS ORDERED AS  RING #,TEBP,TIBP,NBP (or ISAT)
C
      CALL RDI(lpdatsw,.TRUE. , 0 ,.TRUE., 2 ,'LP DATA SWITCH', IERR)
c
c     Inner
c
      call rdr3(te_mult_i,ti_mult_i,n_mult_i,.TRUE.,0.0,.FALSE.,0.0,
     >          'Inner Multipliers       ',IERR)
c
      CALL RDRARN(LPDATI,NLPDATI,MAXINS,-MACHHI,MACHHI,.FALSE.,
     >            0.0,MACHHI,3,'SET OF L. PROBE DATA INNER',IERR)
c
c     Outer
c
      call rdr3(te_mult_o,ti_mult_o,n_mult_o,.TRUE.,0.0,.FALSE.,0.0,
     >          'Outer Multipliers       ',IERR)
c
      CALL RDRARN(LPDATO,NLPDATO,MAXINS,-MACHHI,MACHHI,.FALSE.,
     >            0.0,MACHHI,3,'SET OF L. PROBE DATA OUTER',IERR)
c
c     Adjust the Target conditions - by multiplication factors
c
c     jdemod 
c
c     If the lpdat switch is set to 2 convert particles/s to Amperes for
c     compatibility elsewhere in the code. 
c      
      if (lpdatsw.eq.2) then 
c
c        Note: after conversion the data is then in A and the lpdatsw needs to be
c        set accordingly. 
c
         write(6,'(a)') 'INPUT: LPDATSW : Initial value 2 : '//
     >     'Isats converted from particles/s to Amp : LPDATSW set to 1'
c
         lpdat_conv = ech
         lpdatsw = 1
      else
         lpdat_conv = 1.0
      endif
c
      do in = 1,nlpdati
         lpdati(in,2) =  lpdati(in,2) * te_mult_i
         lpdati(in,3) =  lpdati(in,3) * ti_mult_i
         lpdati(in,4) =  lpdati(in,4) *  n_mult_i * lpdat_conv
      end do
c
      do in = 1,nlpdato
         lpdato(in,2) =  lpdato(in,2) * te_mult_o
         lpdato(in,3) =  lpdato(in,3) * ti_mult_o
         lpdato(in,4) =  lpdato(in,4) *  n_mult_o * lpdat_conv
      end do
C
c     Read in data for core rings - data varies depending on
c     core option selected.
c

      CALL RDRARN(coredat,Ncoredat,MAXINS,-MACHHI,MACHHI,.FALSE.,
     >            -machhi,MACHHI,4,'SET OF CORE DATA',IERR)
c
      CALL RDR(CDPERP,.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERP       ',IERR)
      CALL RDR(CDPERPT,.FALSE.,0.0,.FALSE.,0.0,'X DIFF DPERP-TRAP',IERR)
      if (cdperpt.le.0.0) cdperpt = cdperp
      CALL RDR(cVPINCH,.FALSE.,0.0,.FALSE.,0.0,'PERP PINCH VEL.  ',IERR)
      CALL RDR(CRMI,  .TRUE. ,0.1,.FALSE.,0.0,'IMPURITY ION MASS', IERR)
      CALL RDI(CION,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)
      CALL RDR(CEBD,  .TRUE. ,0.0,.FALSE.,0.0,'CHAR ENERGY EBD',   IERR)
      CALL RDI(CIZEFF,.TRUE. , 0 ,.FALSE., 0 ,'ZEFF(SELF)     ',   IERR)
      CALL RDR(CTEM1, .FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP   ',   IERR)
      CALL RDR(CTEM2 ,.FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP (2)',  IERR)
C
      CALL RDR(CXSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL R POSITION',IERR)
      CALL RDR(CYSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL Z POSITION',IERR)
      CALL RDI(CIZSC, .TRUE.,  1, .TRUE.,CION,'INITIAL IZ STATE',  IERR)
C
      CALL RDR(CNHC,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHC   ',IERR)
      CALL RDR(CNHO,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHO   ',IERR)
      CALL RDR(CLAMHX,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHX ',IERR)
      CALL RDR(CLAMHY,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHY ',IERR)
      CALL RDR(CVCX  ,.TRUE., 0.0,.FALSE.,0.0,'CXREC CONSTANT VCX',IERR)
      CALL RDR(CXNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  XNEAR ',IERR)
      CALL RDR(CYNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  YNEAR ',IERR)
      CALL RDR(CSNORM,.FALSE.,0.0,.FALSE.,0.0,'MEASURE THETA FROM',IERR)
      CSNORM = CSNORM / RADDEG
      CALL RDR(CVHIN, .FALSE.,0.0,.FALSE.,0.0,'INB. PLASMA FLOW V',IERR)
      CALL RDR(CEIN , .FALSE.,0.0,.FALSE.,0.0,'INBOARD ELEC FIELD',IERR)
      CALL RDR(CVHOUT,.FALSE.,0.0,.FALSE.,0.0,'OUT. PLASMA FLOW V',IERR)
      CALL RDR(CEOUT ,.FALSE.,0.0,.FALSE.,0.0,'OUTBRD ELEC FIELD ',IERR)
      CALL RDR(CZENH, .FALSE.,0.0,.FALSE.,0.0,'COLLIS ENHANC ZENH',IERR)
      CALL RDI(CIZSET,.FALSE., 0 ,.FALSE., 0 ,'SET TI=TB AT STATE',IERR)
      CALL RDR(CTRESH,.FALSE.,0.0,.FALSE.,0.0,'SELF-SPU THRESHOLD',IERR)
      CALL RDI(CBOMBZ,.TRUE.,  0 ,.FALSE., 0 ,'BOMBION CHARGE    ',IERR)
      CALL RDI(CBOMBF,.TRUE.,  0 ,.TRUE.,  7 ,'BOMBION FLAG 0:7  ',IERR)
      CALL RDR(CIRF  ,.TRUE., 0.0,.FALSE.,0.0,'IONISE RATE FACTOR',IERR)
      CALL RDR(CSEF  ,.TRUE., 0.0,.FALSE.,0.0,'SPUT ENHANCE FACT.',IERR)
      CALL RDR(CSOLEF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE E(S) ',IERR)
      CALL RDR(CSOLVF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE VH(S)',IERR)
      CALL RDR(CFL   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FL" ',IERR)
      CALL RDR(CFS   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FS" ',IERR)
      CALL RDR(CFRM  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "FRM"',IERR)
      CALL RDR(CKIN  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "KIN"',IERR)
      CALL RDR(CKOUT ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR"KOUT"',IERR)
      CALL RDR(CFRMIN,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMIN"',IERR)
      CALL RDR(CFRMAX,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMAX"',IERR)
c
      CALL RDR(CVBL1 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 1',IERR)
      CALL RDR(CVBM1 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 1  ',IERR)
      CALL RDR(CVBL2 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 2',IERR)
      CALL RDR(CVBM2 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 2  ',IERR)
C
      CALL RDI(IMODE, .TRUE.,  0, .TRUE.,   2  ,'OPERATION MODE',  IERR)
      CALL RDI(NIZS,  .TRUE.,  0, .TRUE.,MAXIZS,'MAX IZ STATE',    IERR)
      CALL RDI(NIMPS, .TRUE.,  1, .TRUE.,MAXIMP,'NO OF IONS',      IERR)
      CALL RDI(NIMPS2,.TRUE.,0,.TRUE.,MAXIMP-NIMPS,'NUM SUP IONS', IERR)
c
c     RESET NIMPS2 - if an Ion injection has been specified and
c                    not a neutral launch - set nimps2 = 0
c
      if (cneuta.eq.1) nimps2 = 0
c
      CALL RDR(FSRATE,.TRUE.,1.E-10,.TRUE.,1.0,'NEUT ITERATE TIME',IERR)
      CALL RDR(QTIM,  .TRUE.,1.E-10,.TRUE.,1.0,'DIV ITERATE TIME', IERR)
      CALL RDR(CPULIM,.TRUE.,0.0,.FALSE., 0.0 ,'CPU TIME LIMIT'  , IERR)
C
C---- READ IN TIME DEPENDENT DATA  (NOTE 128)
C
      CALL RDRAR(DWELTS(0),NQS,MAXIZS+1,0.0,machhi,.FALSE.,'DWELT',IERR)
      IF (IMODE.EQ.1.AND.NQS.LT.NIZS+1) GOTO 1002
      DWELTS(-1) = DWELTS(0)
      CALL RDRAR(DWELFS,NTS,MAXNTS,  0.0,machhi,.TRUE.,'T FACTORS',IERR)
      call rdr(ctimmax,.true.,0.0,.false.,10.0,'DWELL TIME LIMIT ',ierr)
C
C---- READ IN YIELD MODIFIER FUNCTIONS
C
      CALL RDRARN(CYMFS,NYMFS,MAXPTS+1,-machhi,machhi,.true.,-machhi,
     >                            machhi,
     >                            7,'SET YIELD(P,S,CT,PW,CW) VALS',IERR)
c
      CALL RDR(const_yield,.true.,0.0,.FALSE.,0.0,'Fixed Yield',IERR)
c
c---- Read in target plate temperature
c
      CALL RDR(ctargt ,.true.,0.0,.FALSE.,0.0,'Target Temp (K) ',IERR)
      CALL RDR(cwallt ,.true.,0.0,.FALSE.,0.0,'Wall Temp (K)   ',IERR)
      CALL RDR(cwallp ,.true.,0.0,.FALSE.,0.0,'PP Wall Temp (K)',IERR)
      CALL RDRARN(walltemp,NWLTEMP,MAXPTS,-MACHHI,MACHHI,.FALSE.,
     >                  0.0,MACHHI,2,'WALL SEGMENT TEMPERATURES',IERR)
C
C---- READ IN DEBUG OPTS 0:OFF >0 PRINT DIAGNOSTICS EVERY I'TH TIMESTEP
C---- READ IN RANDOM NUMBER GENERATOR SEED
C
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG NEUT ***',IERR)
      CSTEPN = REAL (ISTEP)
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG DIV  ***',IERR)
      CSTEPL = REAL (ISTEP)
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'**DEBUG VELOCITY**',IERR)
      CSTEPV = REAL (ISTEP)
      CALL RDI(CISEED,.FALSE., 0 ,.FALSE., 0 ,'RANDOM NUMBER SEED',IERR)
      CALL RDI(PINISEED,.FALSE.,0,.FALSE.,0,'PIN RANDOM NUM. SEED',IERR)
      CALL RDI(CPRINT,.true., 0 ,.true.,  19 ,'PRINT OPTION',IERR)
      CALL RDI(PINPRINT,.true., 0 ,.true.,1,'PIN DATA PRINT OPT',IERR)
C
      CALL RDI(NITERS,.TRUE. , 1 ,.TRUE. ,20 ,'NO. OF ITERATIONS ',IERR)
      IF (NITERS.GT.1 .AND. (CIOPTB.EQ.3.or.cioptb.eq.9)) GOTO 1003
      CALL RDI(CSTOP ,.TRUE. , 0 ,.TRUE. , 1 ,'STOP WHEN HIT MAIN',IERR)
      CALL RDI(CRECT ,.TRUE. , 0 ,.TRUE. , 99,'RECT GRID CALC.   ',IERR)
      CALL RDI(CZO   ,.TRUE. , 0 ,.FALSE., 0 ,'ZO PARAMETER      ',IERR)
      ZO = REAL(CZO)
      CHZO = 1.56*(1.+1.41*ZO)*(1.+0.52*ZO)/((1.+2.65*ZO)*
     >       (1.+0.285*ZO))
      CALL RDI(CNIWA ,.TRUE. , 0 ,.FALSE., 1 ,'CNIWA OPTION      ',IERR)
C
      CALL RDVMF('VEL. MULT. FACTOR',IERR)
C
      CALL RDR(TLOSS,.FALSE.,0.0,.FALSE.,0.0,'ION LOSS TIME', IERR)
C
      CALL RDR(CPA  ,.FALSE.,0.0,.FALSE.,0.0,'POWER DENSITY', IERR)
      CALL RDR(CK0  ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND', IERR)
      CALL RDR(CK0I ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND IONS',IERR)
C
c slmod begin
      CALL RDI(OFIELD,.TRUE. , 0 ,.TRUE. , 5 ,'EFIELD=0 FOR FILE',IERR)
c
c      CALL RDI(OFIELD,.TRUE. , 0 ,.TRUE. , 4 ,'EFIELD=0 FOR FILE',IERR)
c slmod end
      CALL RDR(CEFLEN,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD SRC LENGTH',IERR)
      CALL RDR(CEFfact,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD COLL-MULT',IERR)
      CALL RDR(CSOLLS,.TRUE.,0.0,.FALSE.,0.0, 'LS DECAY SOL12   ',IERR)
      CALL RDR(CSOLLT,.TRUE.,0.0,.FALSE.,0.0, 'SECOND DECAY DIST',IERR)
      CALL RDR(CFIZ  ,.TRUE.,0.0,.TRUE. ,1.0, 'SOURCE FRACTION  ',IERR)
      CALL RDR(CSOLLR,.TRUE.,0.0,.FALSE.,0.0, 'LR DECAY SOL12   ',IERR)
      CALL RDR(CSOLPR,.TRUE.,0.0,.FALSE.,0.0, 'PR CONST SOL12   ',IERR)
      CALL RDR(CSOLFR,.TRUE.,0.0,.FALSE.,0.0, 'FR SOURCE FRAC.  ',IERR)
      CALL RDI(CSOPT ,.TRUE., 0,.TRUE.,5, 'IONIZATION SOURCE OPT',IERR)
      CALL RDI(CPOPT ,.TRUE., 0,.TRUE.,3, 'RADIATION SOURCE OPT ',IERR)
      CALL RDI(SROOTOPT,.TRUE.,0,.TRUE.,1,'TREAT NEGATIVE ROOTS ',IERR)
      CALL RDI(FLUXROPT,.TRUE.,0,.TRUE.,1,'SRC FLUX RECIRC OPT',  IERR)
      CALL RDRARN(FLUXINFO,FLUXPTS,MAXINS,-MACHHI,MACHHI,.FALSE.,
     >            -machhi,MACHHI,3,'SET OF FLUX DATA',IERR)
      CALL RDR(CZD   ,.FALSE.,0.0,.FALSE.,0.0, 'Z DIVERTOR LIMIT',IERR)
C
C     INJECTION OPTION 2 - INPUT PARAMETERS
C
      CALL RDI(INJIR,.TRUE.,1,.TRUE.,MAXNRS,'RING NUMBER FOR INJ',IERR)
      CALL RDR(INJF1,.TRUE.,0.0,.TRUE.,1.0  ,'INJECTION AREA LB' ,IERR)
      CALL RDR(INJF2,.TRUE.,INJF1,.TRUE.,1.0,'INJECTION AREA UB' ,IERR)
c
      CALL RDR2(FPXMAXO,fpxmaxI,.TRUE.,0.0,.FALSE.,0.0,
     >                                       'FP DIFFUSION SIZE' ,IERR)
      CALL RDR2(FPTIMO,fptimI,.TRUE.,0.0,.FALSE.,0.0,
     >                                       'FP TARGET LOSS T ' ,IERR)
      CALL RDR(CDPERPFP,.FALSE.,0.0,.FALSE.,0.0,'FP DIFF RATE  ' ,IERR)
      IF (CDPERPFP.LT.0.0) CDPERPFP = CDPERP
C
      CALL RDRARN(WLPROB,NWLPROB,MAXPTS,-MACHHI,MACHHI,.FALSE.,
     >            0.0,MACHHI,2,'WALL LAUNCH PROB. MODIFIERS',IERR)
C
      CALL RDI(WLPABS,.TRUE. ,0  ,.TRUE. , 3 ,'ABS WALL PROB    ',IERR)
      CALL RDR(CDRFTV,.FALSE.,0.0,.FALSE.,0.0,'POL. DRIFT VEL.  ',IERR)
      call rdr2(cdrftv_start,cdrftv_end,.false.,0.0,.false.,1.0,
     >                                   'Start and End of DRFTV',ierr)
      CALL RDR(CEMAXF,.FALSE.,0.0,.FALSE.,0.0,'EMAX-FACTOR',      IERR)
      CALL RDR(CEIMP, .TRUE., 0.0,.FALSE.,0.0,'EIMP- WALL LAUNCH',IERR)
C
C     COMMANDS RELATING TO PIN EXECUTION
C
c slmod begin
      CALL RDI(CPINOPT,.TRUE.,0  ,.TRUE. ,4  ,'RUNPIN 0-NO 1-YES',IERR)
c 
c      CALL RDI(CPINOPT,.TRUE.,0  ,.TRUE. ,1  ,'RUNPIN 0-NO 1-YES',IERR)
c slmod end
      CALL RDC(CPINCOM,'COMMAND TO RUN PIN',IERR)
      READ(CPINCOM(11:80),'(A69)') ACTPIN
c
      CALL RDI(IHCORR,.TRUE.,0  ,.TRUE.,1 ,'PIN CELL AREA OPTION',IERR)
      CALL RDI(IHYBRID,.TRUE.,0  ,.TRUE. ,6 ,'HYBRID WALL IN PIN',IERR)
c
c     More PIN related values - PUFFING OPTIONS
c
      call rdi(pinpuff,.true.,0,.true.,2,    'PIN PUFFING OPTION',ierr)
      call rdi(swpvhpf,.true.,0,.true.,1,  'PUFF LOCATION OPTION',ierr)
      CALL RDR(hpcpuf, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF FRAC',IERR)
      CALL RDR(ppcpuf,.true.,0.0,.FALSE.,0.0, 'FLUX-> PUFF OPT 2',IERR)
c
c     CALL RDR(hextrl,.false.,0.0,.FALSE.,0.0,'EXTRA H FLUX-OPT2',IERR)
c     CALL RDR(phxtra,.TRUE.,0.0,.FALSE.,0.0,'PUFF FRAC OF EXTRA',IERR)
c
      CALL RDR(tpufh, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF TEMP ',IERR)
      call rdi2(jhpuf1(1),jhpuf1(2),.false.,0,.false.,1,
     >                                  'PUFF LOCATION INDICES 1',ierr)
      call rdi2(jhpuf2(1),jhpuf2(2),.false.,0,.false.,1,
     >                                  'PUFF LOCATION INDICES 2',ierr)
c
c slmod begin
      CALL RDI(CITERSOL,.TRUE.,0 ,.TRUE. ,2  ,'DO SOL AFTER PIN ',IERR)
c
c      CALL RDI(CITERSOL,.TRUE.,0 ,.TRUE. ,1  ,'DO SOL AFTER PIN ',IERR)
c slmod end
      CALL RDI(CSECSOL ,.TRUE.,-2,.TRUE. ,14 ,'SECONDARY SOL OPT',IERR)
      IF (CSECSOL.EQ.-2) CSECSOL = CIOPTF
      CALL RDI(CSOPT2,.TRUE.,0,.TRUE.,5,  'SECOND ION SOURCE OPT',IERR)
      CALL RDI(NITERSOL,.TRUE.,0 ,.FALSE.,maxpiniter,
     >                                        'NO. OF PIN ITER. ',IERR)
C
      CALL RDRARN(PLATCO,NPLAT,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     >           -MACHHI,MACHHI,4,'PLATE COORDINATES',IERR)
      CALL RDRARN(WALLCO,NWALL,MAXPTS,-MACHHI,MACHHI,.FALSE.,
     >           -MACHHI,MACHHI,1,'WALL COORDINATES',IERR)
      CALL RDRARN(WALLCO2,NWALL2,MAXPTS,-MACHHI,MACHHI,.FALSE.,
     >           -MACHHI,MACHHI,1,'WALL COORDINATES',IERR)
      CALL RDI(CMAXGENS,.TRUE., 0,.FALSE.,  0, 'MAX. GENERATIONS',IERR)
      CALL RDI(CIRHR,  .TRUE., 1,.TRUE.,MAXNRS,'HI-RES RING     ',IERR)
c
      CALL RDR(CNIN,.TRUE.,0.01,.FALSE.,0.0,'COSINE DIST. POWER',IERR)
      CALL RDR(CVAMULT,.TRUE.,LO,.FALSE.,0.0,'VELOCITY MULTIPLIER',IERR)
      CALL RDR(CVRMULT,.TRUE.,LO,.FALSE.,0.0,'REC. VEL. MULT',IERR)
c
      CALL RDRAR(cleaks,cleaksn,maxpts,0.0,machhi,.TRUE.,'LEAK- S',IERR)
c
c     Set up the flag indicating that leak checking code should be activ
c
      if (cleaksn.gt.0) then
         checkleak = .true.
      else
         checkleak = .false.
      endif
c
c     The following variables all apply to SOL option 21 which attempts
c     to model a detached plasma divertor configuration.
c
      call rdi(s21refsw,.true.,0,.true.,3,'S21 Length Ref Switch',ierr)
c
c     Load up Outer/Inner values
c
      call rdr2(terat,terati,.true.,0.0,.false.,0.0,
     >                                        'Te Ratio at L1 O/I',ierr)
      call rdr2(tirat,tirati,.true.,0.0,.false.,0.0,
     >                                        'Ti Ratio at L1 O/I',ierr)
      call rdr2(nrat,nrati,.false.,0.0,.false.,0.0,
     >                                        'Ne Ratio at L1 O/I',ierr)
      call rdr2(nalph,nalphi,.true.,0.0,.false.,0.0,
     >                                     'Ne Exponent at L1 O/I',ierr)
      call rdr2(qrat,qrati, .false.,0.0,.false.,0.0,
     >                                        'Emitted energy rat',ierr)
      call rdr2(l1rat,l1rati,.true.,0.0,.false.,0.0,
     >                                        'L1*SMAX boundary 1',ierr)
      call rdr2(l2rat,l2rati,.true.,0.0,.false.,0.0,
     >                                        'L2*SMAX boundary 2',ierr)
      call rdr2(lvrat,lvrati,.true.,0.0,.false.,0.0,
     >                                        'V->0 at lvrat*SMAX',ierr)
      call rdr2(vbmult,vbmulti,.true.,0.0,.false.,0.0,
     >                                        'Vmult for Region A',ierr)
c
c     Read in by ring parameters for SOL21 if any are specified.
c
      CALL RDRARN(S21PARMI,NS21I,MAXNRS,-machhi,machhi,.FALSE.,
     >           -MACHHI,MACHHI,9,'SOL21 PARAMS INNER',IERR)
      CALL RDRARN(S21PARMO,NS21O,MAXNRS,-MACHHI,MACHHI,.FALSE.,
     >           -MACHHI,MACHHI,9,'SOL21 PARAMS OUTER',IERR)
c
c     Impose an external ABSFAQ ... unless it is <= 0 then use the
c     DIVIMP calculated value.
c
      call rdr(nabsfac,.false.,0.0,.false.,0.0, 'ABSFAC modifier',ierr)
c
c     The following set of numbers specify the characteristic values of
c     the ASDEX U style grid file to be read in. The Number of rings
c     and elements per ring (a ring is a set of polygons in a row or
c     indexed by the same ring number), the cutring,(which specifies the
c     end of  rings for the core and trapped plasma), and the cutpts 1
c     and 2 which specify the points on the rings numbered less than the
c     cut ring where the splits for core and trapped or private plasma
c     occur.
c

      CALL RDI(maxrings,.TRUE.,0  ,.true.,maxnrs,'MAXRINGS in AUG',IERR)
      CALL RDI(maxkpts,.TRUE.,0  ,.true.,maxnks,'MAX PTS in AUG',IERR)
      CALL RDI(cutring,.TRUE.,0  ,.true.,maxrings,'CUTRING in AUG',IERR)
      CALL RDI(cutpt1,.TRUE.,0  ,.true.,maxkpts+1,'CUTPT1 in AUG',IERR)
      CALL RDI(cutpt2,.TRUE.,0  ,.true.,maxkpts+1 ,'CUTPT2 in AUG',IERR)
      CALL RDI(nfla,.TRUE.,1,.true.,maxnfla,'NUM.FLUIDS IN B2 BG',IERR)
      CALL RDI(readaux,.TRUE.,0,.true.,1,'READ AUX BG FILE',IERR)
c
      CALL RDR(Cstgrad,.TRUE.,0.0,.TRUE.,1.0, 'TGRAD FORCES -> 0',IERR)
c
c     Call the subroutine to read in the characteristic values required
c     for the AS-AD SOL model - SOL option 22
c
      call readsol(ierr)
c
      call rdr(ctestsol,.true.,-1.0,.false.,1.0,'TEST SOL OPT ONLY',
     >         ierr)
c
c     Overall Recombination Options
c
      call rdi(crecopt,.TRUE.,0  ,.true. ,4  ,'Recomb. Calc Opt',IERR)
      call rdr(treccut,.TRUE.,0.0,.false.,0.0,'Recomb. cutoff T',IERR)
c
c     Dperp extractor options
c
      CALL RDI(dpmethod,.TRUE.,0  ,.true.,2 ,'DPERP EXT METHOD',IERR)
      call rdi(dpsuml,.TRUE.,0  ,.true.,1 ,'DPERP EXT SUM LIMIT',IERR)
      call rdi(dpouter,.TRUE.,0  ,.true.,1 ,'DP EXT OUTER RING',IERR)
      call rdi(dpconv,.TRUE.,0  ,.true.,1 ,'DP EXT CONVECT LOSS',IERR)
      call rdi(dpfluxopt,.TRUE.,0  ,.true.,1 ,'CELL CENTRE FLUX',IERR)
      call rdi(dpavopt,.TRUE.,0  ,.false.,1 ,'AVERAGE DPERP OPT',IERR)
      call rdi(dprcopt,.TRUE.,0  ,.true.,2 ,'MAJOR RADIUS CORR',IERR)
      call rdi(dpsmooth,.true.,0,.false.,0 ,'GRADIENT SMOOTHING',IERR)
      call rdi(dpnav,.true.,-1,.true.,2 ,   'GRADIENT CALC METH',IERR)
      call rdi(dparea,.TRUE.,0  ,.true.,1 ,'CELL BOUND AREA',IERR)
      call rdi(dpploss,.TRUE.,0  ,.true.,2 ,'POWER LOSS TERMS',IERR)
      call rdi(dporth,.TRUE.,0 ,.true.,1,'DP NON-ORTH GRADIENT',IERR)
      call rdr(dppei,.TRUE.,0.0,.false.,0.0,'PEI Correction',IERR)
      call rdr(dprec,.TRUE.,0.0,.false.,0.0,'EXT Recycle Frac',ierr)
      call rdr(dpxpratio,.FALSE.,0.0,.false.,0.0,'Dp/Xp Ratio',ierr)
c
c     The following values apply to the specifiable SOL option that
c     allows two part linearly interpolated fits for each of
c     Te, Ti, ne and vb to be specified.
c
      call rdr(ctes1,.TRUE.,0.0,.false.,0.0,'Te distance  1',ierr)
      call rdr(ctef1,.TRUE.,0.0,.false.,0.0,'Te fact/mult 1',ierr)
      call rdr(ctes2,.TRUE.,0.0,.false.,0.0,'Te distance  2',ierr)
      call rdr(ctef2,.TRUE.,0.0,.false.,0.0,'Te fact/mult 2',ierr)
      call rdr(ctis1,.TRUE.,0.0,.false.,0.0,'Ti distance  1',ierr)
      call rdr(ctif1,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 1',ierr)
      call rdr(ctis2,.TRUE.,0.0,.false.,0.0,'Ti distance  2',ierr)
      call rdr(ctif2,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 2',ierr)
      call rdr(cnes1,.TRUE.,0.0,.false.,0.0,'Nb distance  1',ierr)
      call rdr(cnef1,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 1',ierr)
      call rdr(cnes2,.TRUE.,0.0,.false.,0.0,'Nb distance  2',ierr)
      call rdr(cnef2,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 2',ierr)
      call rdr(cvbs1,.TRUE.,0.0,.false.,0.0,'vb distance  1',ierr)
      call rdr(cvbf1,.TRUE.,0.0,.false.,0.0,'vb fact/mult 1',ierr)
      call rdr(cvbs2,.TRUE.,0.0,.false.,0.0,'vb distance  2',ierr)
      call rdr(cvbf2,.TRUE.,0.0,.false.,0.0,'vb fact/mult 2',ierr)
c
c     Reciprocating/Fast Scanning Probe - R location.
c
c     .... and horizontal probe location
c
      call rdi(rlocnum,.true.,1,.false.,0, 'Crossing number',ierr)
      call rdr(crploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)
      call rdi(zlocnum,.true.,1,.false.,0, 'Crossing number',ierr)
      call rdr(czploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)
c
c     Marfe Core options - for option 4,5,6
c
      call rdr(corefv,.true.,0.0,.true.,0.5,'CORE-VEL FRAC',ierr)
      call rdr(coreft,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC',ierr)
      call rdr(corefv2,.true.,0.0,.true.,0.5,'CORE-VEL FRAC2',ierr)
      call rdr(coreft2,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC2',ierr)
c
c     Factor for temperature gradient force modifier option
c
      call rdr(fgradfact,.false.,0.0,.false.,0.0,'FGRAD MOD FACT',ierr)
c
c     TEMPORARY - CORRECT Tgrad forces for Edge2D bug.
c
c      CALL RDI(fixtgrad,.TRUE.,0  ,.true.,1 ,'FIX E2D TGRAD',IERR)
c
c
c     Read in the INPUT options for SOL option 23
c
      call rdi(readin_sol23_params,.TRUE.,0 ,.true.,1,
     >                       'READ SOL23 PARMAS OPT',IERR)
c
      if (readin_sol23_params.eq.1) then

         call read_sol23_params(ierr)

      endif
c
c     Old sol23 input
c
c      call rdi(sol23_intopt,.TRUE.,0 ,.true.,1,
c     >                                  'SOL23-INTEGRATION OPT',IERR)
c      call rdi(sol23_adaptgrid,.TRUE.,0 ,.true.,5,
c     >                                  'SOL23-ADAPTIVE GRID',IERR)
c      call rdi(sol23_bndcond,.TRUE.,0 ,.true.,1000,
c     >                                  'SOL23-BND CONDITION OPT',IERR)
c      call rdi(sol23_seed,.TRUE.,1 ,.true.,2,
c     >                                  'SOL23-SEED PLASMA OPT',IERR)
c      call rdr(sol23_izlen,.true.,0.0,.false.,0.5,
c     >                                  'SOL23-IONIZATION LEN',ierr)
c      call rdr(sol23_izlam,.true.,0.0,.false.,0.5,
c     >                                  'SOL23-DECAY LENGTH',ierr)
c      call rdr(sol23_izoffset,.true.,0.0,.false.,0.5,
c     >                                  'SOL23-OFFSET LENGTH',ierr)
c      call rdr(sol23_momlen,.true.,0.0,.false.,0.0,
c     >                                  'MOMENTUM LENGTH',ierr)
c
c      call rdr(sol23_sm_s,.true.,0.0,.true.,1.0E8,
c     >                                  'SMOOTHING ALONG B',ierr)
c      call rdr(sol23_sm_r,.true.,0.0,.true.,1.0E8,
c     >                                  'SMOOTHING ACROSS B',ierr)
c      call rdr(sol23_relax,.true.,0.0,.true.,1.0,
c     >                                  'RELAXATION WITH PIN',ierr)
c      call rdr(sol23_maxtol,.true.,0.0,.true.,1.0,
c     >                             'MAX TOLERANCE IN CFD SOLN',ierr)
c      call rdr(sol23_rmstol,.true.,0.0,.true.,1.0,
c     >                             'RMS TOLERANCE IN CFD SOLN',ierr)
c      call rdi(sol23_perp,.true.,0 ,.true.,9 ,
c     >                                  'Q PERP OPTION',ierr)
c      call rdr(sol23_maxpow,.true.,-1.0,.true.,1.0E8,
c     >                                  'MAX. POWER LIMIT',ierr)
c
c
c
c     READ IN THE REST OF THE DATA FILE - MAXIMUM 50 LINES
c
c     The reason for this is so that the selected PIN input options
c     can be printed in the DIVIMP output file. At this time we do
c     not want to try interpreting them - simply echo then so that
c     it is clear which options were in use when a particular DIVIMP
c     case was run.
c
      CALL RDCAR(CNIMBIN,nlines,maxlines,'NIMBUS NAMELIST INPUT',IERR)
c
c--------------------------------------------------------------------------
c 
c     CHECKS AND VERIFICATION OF SOME INPUT VALUES:
c
c--------------------------------------------------------------------------   
c
c     SOME INPUT DATA VERIFICATION
c
C
C     If the target option specified requires the platco array to
C     store data (as in target option 1), then the geometry data
C     can not be loaded on top.
C
C     Note: if the target option is 3 or the wall option 4 or 5
C           and cgeoopt has not been specified then an error
C           has occurred.
C
      IF (CGEOOPT.EQ.-1.AND.(ctargopt.eq.3.or.
     >    cneur.eq.3)) then
         call prc('Invalid Input ... target or wall option that')
         call prc('requires pre-loaded data has been specified ')
         call prc('but specific pre-loaded grid data has not   ')
         call prc('selected.')
         stop
      endif
c
c     Check to make sure that Trap option 3 and a double-null ITER
c     grid are not specified at the same time - these options share
c     a common array with different meanings at this time and are thus
c     incompatible.
c
      if (cgridopt.eq.2.and.ctrap.eq.3) then
         call prc('ITER grid is incompatible with TRAP option 3.')
         call prc('Please respcify input - program ending ...')
         stop
      endif
c
c     Error - make sure PIN is running if NIMBUS wall specified.
c
      if (cpinopt.eq.0.and.(cneur.eq.5.or.ctrap.eq.5)) then
c
         call prc('ERROR: Incompatible NEUTRAL Wall/Trap'//
     >            ' Wall Options')
         call prc('Neutral Wall Option 5 or Trap Wall Option 5')
         call prc('specified without corresponding call to'//
     >            ' PIN/NIMBUS')
         call prc('Program stopping')
c
c        Write to screen for easy notification
c
         write(0,*) 'ERROR: Incompatible NEUTRAL Wall/Trap'//
     >            ' Wall Options'
         write(0,*) 'Neutral Wall Option 5 or Trap Wall Option 5'
         write(0,*) 'specified without corresponding call to'//
     >            ' PIN/NIMBUS'
         write(0,*) 'Program stopping'
c
         stop
      endif
c
c     Warning - if one of the SOL - Plasma Decay or Temperature Gradient
c               options is specified as 99 and at least one is NOT
c             - issue a Warning!
c
      if (cioptg.eq.99.or.cioptf.eq.99.or.
     >    cioptk.eq.99.or.cioptl.eq.99) then

         if (cioptg.ne.99.or.cioptf.ne.99.or.
     >       cioptk.ne.99.or.cioptl.ne.99) then
c
c           Issue Warning - since mixed file data specified.
c
            call prb
            call prc(' WARNING --- WARNING --- WARNING !!! ')
            call prb
            call prc(' POSSIBLY INCONSISTENT INPUT: ')
            call pri(' SOL OPTION          : ',cioptf)
            call pri(' PLASMA DECAY OPTION : ',cioptg)
            call pri(' Te GRADIENT OPTION  : ',cioptk)
            call pri(' Ti GRADIENT OPTION  : ',cioptl)
            call prb
            call prc(' AT LEAST ONE HAS BEEN SPECIFIED AS'//
     >               ' FILE INPUT (99)')
            call prc(' BUT NOT ALL - THIS MAY PRODUCE'//
     >               ' UNEXPECTED RESULTS!')
            call prb
c
         endif
c
      endif
c
c     e2dtargopt=2 requires flux data for each flux tube to be specified.
c     e2dtargopt=3,4
c
      if ((e2dtargopt.eq.2.or.e2dtargopt.eq.3.or.e2dtargopt.eq.4.or.
     >         e2dtargopt.eq.5)
     >        .and.(fluxpts.le.0.and.readaux.ne.1)) then
         call prc('ERROR: Incompatible Input.')
         call prc('- Fluid code target condition option 2, 3'//
     >            ' or 4 has been')
         call prc('  specified without corresponding fluxes for')
         call prc('  each target element - FLUXINFO - or using')
         call prc('  FC Target Option and Read Auxiliary Input to ')
         call prc('  extract down fluxes from an external file.')
         call prc('Program stopping')
         stop
      elseif (e2dtargopt.ne.0.and.cre2d.eq.0) then
         call prc('ERROR: Incompatible Input.')
         call prc('- A Fluid Code  target condition option has been')
         call prc('  specified without the corresponding option')
         call prc('  to also read the Fluid Code background'//
     >            ' for reference.')
         call prc('Program stopping')
         stop

      endif
c
c     Error - CIOPTI=8 is invalid except for Carbon impurity in a
c             deuterium plasma.
c
      if ((ciopti.eq.8.or.ciopti.eq.9).and.
     >    (crmb.ne.2.0.or.crmi.ne.12.0)) then
c
         call errmsg('INPUT ERROR:','Incompatible Input:'//
     >       ' CX Recombination option 8 is ONLY compatible '//
     >       '  with a carbon impurity in a deuterium plasma.')

         stop
c
       endif
c
c      IF Pinch option 4 is specified then a Probability Distribution
c      Function must be loaded - if not - the pinch option is turned 
c      OFF. 
c
       if (pinchopt.eq.4.and.pinch_npdf.le.0) then 
c
          call errmsg('INPUT ERROR:','PINCH OPTION 4'//
     >                ' SELECTED BUT PDF NOT SPECIFIED IN INPUT :'//
     >               ' PINCH OPTION TURNED'//
     >               ' OFF (SET TO 0)')
c
          pinchopt = 0
c
       endif
C
c      CSTMAX is set in the rundiv main program
c
c      CSTMAX = 10.0 / QTIM
c
c-------------  INITIALIZATION ROUTINES --------------------
c
c slmod begin - new 
c jdemod - note routine now used to assign some defaults to unstructured values
c          after the input file has been completely read in
      CALL ValidateUnstructuredInput

      CALL InitializeRelaxation
c slmod end
c
c     jdemod
c
c     Call init_modules to load some global variables into specific modules private storage
c
c     One example is the mtc module implementing momentum transfer collisions. 
c
      call init_modules(nizs)
c
c
c
      RETURN
C
 1002 CALL PRC ('READIN: NOT ENOUGH DWELL TIMES GIVEN')
      IERR = 1
      RETURN
 1003 CALL PRC ('READIN: ITERATIONS > 1 INCOMPATIBLE WITH COLL OPT 3')
      IERR = 1
      RETURN
      END
C
C
C
