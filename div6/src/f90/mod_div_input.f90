module mod_div_input
  use mod_params
  use mod_reader


  implicit none

  !
  ! input_file_type = 0 (untagged or structured), =1 (tagged or unstructured)
  ! Option 0 uses the legacy input code in subroutine readin
  !

  integer :: input_file_type = 1




contains


  subroutine read_params
    use mod_params
    use mod_io_units
    implicit none

    ! jdemod
    ! - access the input file
    ! - read in any updates to the parameter values from their default values
    ! - determine the type of input file - tagged or untagged 

    ! parameter tags start with the @ sign



    ! Also look for tags defining maxizs


!
!     jdemod - After NIZS and CION have been read in - a value can be
!              assigned to maxizs and used for storate allocation      
!
!     - moved from allocate_storage_div - move to after input file read
!
      maxizs   = max(nizs,cion)
!
!     Allocate the dwelts array which depends on maxizs and which is read
!     later in the input file - problem?? - allocated arrays have their size
!     in the number of elements entered in the input - just need to avoid referencing
!     maxizs when referring to them - maxnts = nts ?
!      
!     move to after read in?  


'+S05    Atomic number of impurity ions     Zi           '      6
      CALL RDI(CION,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)

'+I15    Maximum ionization state                        '      6
      CALL RDI(NIZS,  .TRUE.,  0, .FALSE.,MAXIZS,'MAX IZ STATE',   IERR)


    
  end subroutine read_params
  
    

  subroutine read_input_file

    implicit none





  end subroutine read_input_file







  






end module mod_div_input




  c     -*Fortran*- 
c     @PROCESS NOOPT
      SUBROUTINE READIN (TITLE,desc,equil,NIZS,NIMPS,NIMPS2,CPULIM,
     >                   IERR,NYMFS,NITERS)
c      SUBROUTINE READIN (TITLE,NIZS,NIMPS,NIMPS2,CPULIM,
c     >                   IERR,NYMFS,NITERS)
      use error_handling
      use debug_options
      use ero_interface
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_cadas
      use mod_dynam4
      use mod_dynam5
      use mod_diagvel
      use mod_cedge2d
      use mod_pindata
      use mod_adpak_com
      use mod_promptdep
      use mod_reiser_com
      use mod_fperiph_com
      use mod_driftvel
      use mod_sol23_input
      use mod_slcom
      use mod_sol22_input
      use comhc
      use mod_lambda
      implicit none

      INTEGER   IERR,NIZS,NIMPS,NYMFS,NITERS,NIMPS2
      REAL      CPULIM
      CHARACTER TITLE*(*),equil*(*),desc*(*)
!
!  *********************************************************************
!  *                                                                   *
!  *  READIN:   READS IN THE DATAFILE, PERFORMS VALIDITY CHECKS, ETC.  *
!  *                                                                   *
!  *                                                                   *
!  *  CHRIS FARRELL    FEBRUARY 1989                                   *
!  *                                                                   *
!  *  JAMES SPENCE     NOVEMBER 1990                                   *
!  *                                                                   *
!  *********************************************************************
!
!
      INTEGER NQS,ISTEP,in,id
      REAL    ZO
!
!     jdemod - option for converting LP data input from particles/s to A/s
!     
      real lpdat_conv
!
!     Option indicating if the SOL23 parameter list is included in the data
!     file.
!
      integer readin_sol23_params


      call pr_trace('READIN','START')

      CALL InitializeVariables


      
'+A01 Ref DIV Title' '86575 1650 ms  best case'
      CALL RDC (TITLE, 'TITLE FOR RUN', IERR)

'+A02 DIV Desc:' 'Text Describing the case'
      CALL RDC (desc, 'DESCRIPTION OF RUN', IERR)

'+A03 Equil File Name' 'sonnet_oskn_c08'
      call rdc (equil,'equilibrium file name',ierr)

'+G01    Grid Option    0-JET 1-ASDEX 2-ITER             '     3
      call rdi (cgridopt,.true.,0,.true.,9,'GRID OPTION          ',ierr)

'+G02    Non-Orthogonal Grid option 0-off 1-JET N.O.     '     3
      call rdi (northopt,.true.,0,.true.,3,'NON-ORTHO. GRID OPT  ',ierr)

'+G03    Parallel Distance Option 0-centers 1-edges      '     1
      call rdi (pdopt,.true.,0,.true.,1,   'PARALLEL DIST (S) OPT',ierr)

'+G04    Cross-field Distance Option 0-centres 1-edges   '     1
      call rdi (cfdopt,.true.,0,.true.,1,'Proper Cross-field dist',ierr)

'+G05    RZ calculation option 0-centers 1-Actual RZ     '     1
      call rdi (rzopt,.true.,0,.true.,3,   'CALCULATE ACTUAL R,Z ',ierr)

'+G06    XY Grid option 0-off 1-on                       '     0
      call rdi (xygrid,  .true.,0,.true.,1,'XY GRID OPTION       ',ierr)

'+G07    Cell Area Calculation Option 0-approx 1-polygon '     1
      call rdi (cvolopt,.true.,0,.true.,1, 'Cell Volumes from PGs',ierr)

'+S01    On-AXIS Toroidal B-field value                  '    2.0
      call rdr (cbphi,.true.,0.0,.false.,0.0,'ON-AXIS B-FIELD',    ierr)

'+D01    Source Data Option 0-Nocorona 1-Adas            '     2
      call rdi (cdatopt,.true.,0,.true.,3, 'Rad/ioniz data source',ierr)

'+D02    USERID for ADAS H data (*=use central database) '   '*'
      call rdc (useridh,'ADAS H userid',ierr)

'+D03    Year for ADAS H data                            '    96
      call rdi (iyearh,.true., 0,.true.,99,'ADAS H year          ',ierr)

'+D04    USERID for ADAS Z data (*=use central database) '   '*'
      call rdc (useridz,'ADAS Z userid',ierr)

'+D05    Year for ADAS Z data                            '    96
      call rdi (iyearz,.true., 0,.true.,99,'ADAS Z year          ',ierr)

'+D06 B2FRATES name:' '/u/progs/div6/adpak/C_rates.strahl'
      call rdc (mcfile,'Name of File for ADPAK/INEL atomic data',ierr)

'+T01    Ionization     0/1/2old 3ei 4no 5ei/dis 6no/dis '     3
      CALL RDI (CIOPTA,.TRUE., 0,.TRUE., 6,'IONIZATION OPT       ',IERR)

'+T02    Collision opt  0std 1inf 2zeff 3tpara           '    13
      CALL RDI (CIOPTB,.TRUE., 0,.TRUE.,14,'COLLISION OPT        ',IERR)

!     All other Reiser input is unstructured and optional 
'+T03    REISER                                          '     0
      CALL RDI (CIOPTR,.TRUE., 0,.TRUE., 2,'REISER COLLISION OPT ',IERR)



'+T04    Friction opt   0std 1inf 2tpara                 '     0
      CALL RDI (CIOPTC,.TRUE., 0,.TRUE., 4,'FRICTION OPT         ',IERR)

'+T05    Heating opt    0std 1inf 2zero 3Ti              '     0
      CALL RDI (CIOPTD,.TRUE., 0,.TRUE., 3,'HEATING OPT          ',IERR)

'+I01    Injection opt  1/2/3                            '     7
      CALL RDI (CIOPTE,.TRUE., 0,.TRUE.,14,'INJECTION OPT        ',IERR)

'+P01    SOL option     0,1,1a,2,3,4,5,9,10  99file      '     6
      CALL RDI (CIOPTF,.TRUE.,-1,.TRUE.,99,'SOL OPT              ',IERR)


'+P02    Core Option    0,1,2,3                          '     1
      call rdi (ccoreopt,.true.,-1,.true.,28,'Core Option        ',ierr)

      ! move to after input file read
      if (ccoreopt.gt.6.and.ccoreopt.ne.28) then
        write(0,*) 'ERROR READIN: INVALID CORE OPTION'
        ierr=1
        return
      endif

'+P03    Plasma decay   0std                 99file      '    91
      CALL RDI (CIOPTG,.TRUE., 0,.TRUE.,99,'PLASMA DECAY OPT     ',IERR)


      ! move to after input file read

      ! jdemod - code needs to move to start of bgplasma likely - or at least after input file read
      !...  Need to store these because they are overwritten in BGPLASMA but are needed in GETMODEL:
      orgcioptf = cioptf
      orgcioptg = cioptg


'+P04 ' 'BG PLASMA Options by Ring (PlasDec opts 90&91)  '
'    R1,R2,Sect, PlasDec, Sol, Teg, Tig, Core, Efield'     0
!
!     Example
!     10 23  3.0      4.0 22.0  0.0  0.0   0.0     3.0
!     26 33  3.0      4.0 22.0  0.0  0.0   0.0     3.0
      CALL RDRARN(BGPLASOPT,NBGPLAS,2*MAXNRS,-MACHHI,MACHHI,.FALSE.,
     >            0.0,MACHHI,11,'SET OF BG PLASMA OPTIONS BY RING',IERR)


'+F01    Read EDGE2D BG for reference  0=No  1=Yes       '     2
      call rdi (cre2d ,.TRUE., 0,.true.,2 ,'READ E2D BG FOR REF  ',ierr)

'+F02    Use EDGE2D Target Data Option 0=Off 1=Reg 2=Flux'     0
      call rdi (e2dtargopt,.TRUE., 0,.true.,6 ,'EDGE2D TARG COND',ierr)

'+T06    CX Recomb opt  0off 1on 2Vcx                    '     8
      CALL RDI (CIOPTI,.TRUE., 0,.TRUE., 9,'CX RECOMB OPT        ',IERR)

'+I02    First diffuse  0inst 1random 2tpara             '     0
      CALL RDI (CDIFOP,.TRUE., 0,.TRUE., 2,'FIRST DIFFUSE OPT    ',IERR)

'+T07    Dperp option   0const 1vary                     '     0
      CALL RDI (CIOPTJ,.TRUE., 0,.TRUE., 4,'DPERP OPTION         ',IERR)

'+T08 TN1272 Perpendicular step option 0-normal 1-core    '     3
      call rdi (cdiffopt,.true.,0,.true.,3,'PERP STEP PROB OPTION',ierr)

'+T09 TN14?? Pinch Velocity Option 0=off 1=all 2=main SOL '     0
      call rdi (pinchopt,.true.,0,.true.,15,'PINCH VELOCITY OPTION',ierr)

'+Q01    TeB Gradient   0lin 1lin/lin 2p/a   99file      '     0
      CALL RDI (CIOPTK,.TRUE., 0,.TRUE.,99,'TEB GRADIENT OPTION  ',IERR)

'+Q02    TiB Gradient   0lin 1lin/lin 2p/a   99file      '     0
      CALL RDI (CIOPTL,.TRUE., 0,.TRUE.,99,'TIB GRADIENT OPTION  ',IERR)

'+T10    TeB Grad Coeff 0off 1on                         '     1
      CALL RDI (CIOPTM,.TRUE., 0,.TRUE., 3,'TEB GRAD COEFF OPTION',IERR)

'+T11    TiB Grad Coeff 0off 1on                         '     1
      CALL RDI (CIOPTN,.TRUE., 0,.TRUE., 3,'TIB GRAD COEFF OPTION',IERR)

'+Q03 TN???? Te,Ti Flatten Option                         '     0
      call rdi (cflatopt,.true.,0,.true.,3,'Flatten Te,i option',ierr)

'+Q04 TN1447 Te flat upstream for S > SMAX * Teg cutoff = '    0.0
      call rdr (ctegcut,.FALSE.,0.0,.TRUE.,0.5,'TEB GRAD. CUTOFF',ierr)

'+Q05 TN1447 Ti flat upstream for S > SMAX * Tig cutoff = '    0.0
      call rdr (ctigcut,.FALSE.,0.0,.TRUE.,0.5,'TIB GRAD. CUTOFF',ierr)

'+T12    T-GRAD Force Modification Function 0-off 1-UEDGE'     0
      call rdi (fgradopt,.true.,0,.true.,4,'GRAD FORCE MOD OPTION',ierr)

'+P05    Trap Tgrad Opt 0off 1on                         '     1
      CALL RDI (CIOPTO,.TRUE., 0,.TRUE., 4,'TRAP TGRAD OPTION    ',IERR)

'+I03    Control switch 0atoms 1ions                     '     0
      CALL RDI (CNEUTA,.TRUE., 0,.TRUE., 1,'CONTROL SWITCH       ',IERR)

'+N01    Launch option  0distrib 1point 2asymp 3tip 4wall'     3
      CALL RDI (CNEUTB,.TRUE., 0,.TRUE., 7,'LAUNCH OPTION        ',IERR)

'+N02    Vel/angle flag 0-11                             '    14
      CALL RDI (CNEUTC,.TRUE., 0,.TRUE.,20,'VEL/ANGLE FLAG       ',IERR)
!
!     IF CNEUTH OR CNEUTI ARE -1 THEY NEED TO BE SET TO THE VALUES
!     READ IN FOR CNEUTB AND CNEUTC.
!
!
!     THESE VARIABLES WERE ADDED TO ALLOW LAUNCHES FROM THE WALLS,
!     IN ADDITION TO THE NORMAL PLATE LAUNCHES. HOWEVER, THE
!     CAPABILITY CAN BE USED TO SIMULATE TWO DIFFERING DISTRIBUTIONS
!     OF PARTICLES COMING FROM THE PLATES. PERHAPS DUE TO
!     PHYSICAL AND CHEMICAL SPUTTERING SIMULATNEOUSLY.
!
'+N03 TN487 Supplemental Launch Option (as above)         '     4
      CALL RDI (CNEUTH,.TRUE.,-1,.TRUE., 5,'SUP LAUNCH OPTION    ',IERR)

'+N04 TN487 Supplemental V/A flag      (as above)         '     3
      CALL RDI (CNEUTI,.TRUE.,-1,.TRUE.,19,'SUP VEL/ANGLE FLAG   ',IERR)

'+N05    Initial Neutral Vel/Ang flag (-1=above,0-13)    '    -1
      CALL RDI (NVAOPT,.TRUE.,-1,.TRUE.,19,'NEUT VEL/ANGLE FLAG  ',IERR)

      ! move to AFTER the input file has been completely read in
      IF (NVAOPT.EQ.-1) NVAOPT = CNEUTC
      IF (CNEUTH.EQ.-1) CNEUTH = CNEUTB
      IF (CNEUTI.EQ.-1) CNEUTI = CNEUTC
    

'+N06 TN1490 Supplemental 2D Neutral Launch 0=off 1=UEDGE '     0
      CALL RDI (neut2d_opt,.TRUE.,0,.TRUE., 2,
     >                                    'EXTRA 2D LAUNCH OPTION',IERR)
'+N07 TN1490 V/A Flag for 2D Neutral Launch (as regular)  '    15
      CALL RDI (neut2d_vaopt,.TRUE.,-1,.TRUE.,20,
     >                                         'EXTRA 2D V/A FLAG',IERR)

'+D07    Sputter data option 1-old 2-93                  '     4
      CALL RDI (CSPUTOPT,.TRUE., 1,.TRUE., 8,'SPUTTER SOURCE OPT ',IERR)

'+D08    Chemical Sputter Data Option                    '     9
      CALL RDI (CCHEMOPT,.TRUE., 1,.TRUE.,12,'CHEMSPUT SOURCE OPT',IERR)

'+N08    Sputter option 0std 1special 2mix 3self 4selfva1'     0
      CALL RDI (CNEUTD,.TRUE., 0,.TRUE., 8,'SPUTTER OPTION       ',IERR)

'+N09 TN1209 Secondary Sputter Option                     '     5
      CALL RDI (CNEUTD2,.TRUE.,-1,.TRUE.,8,'2ND SPUTTER OPTION   ',IERR)

      ! move to after input file read in
      if (cneutd2.eq.-1) cneutd2 = cneutd

'+I04    Self- Sputtering Option 0-off 1-on              '     1
      CALL RDI (CSELFS,.TRUE., 0,.TRUE., 2,'SELF-SPUTTER OPTION  ',IERR)

'+N10    Normal option  0std 1fromT                      '     0
      CALL RDI (CNEUTE,.TRUE., 0,.TRUE., 2,'NORMAL OPTION        ',IERR)

'+N11    NEUT spreading 0off 1on                         '     0
      CALL RDI (CNEUTF,.TRUE., 0,.TRUE., 1,'NEUT SPREADING       ',IERR)

'+I05    Init ion Vel   1                                '     1
      CALL RDI (CNEUTG,.TRUE., 0,.TRUE., 3,'INITIAL ION VELOCITY ',IERR)

'+G08 T   Ion Wall Option        0 to 2                   '     2
      CALL RDI (CIONR ,.TRUE., 0,.TRUE., 2,'ION WALL OPTION      ',IERR)

'+G09 T   Neutral Wall Option    0 to 4                   '     2 4
      CALL RDI (CNEUR ,.TRUE., 0,.TRUE., 7,'NEUTRAL WALL OPTION  ',IERR)

'+G10 T   Trap Wall Option       0 to 4                   '     3 4
      CALL RDI (CTRAP ,.TRUE., 0,.TRUE., 8,'TRAP WALL OPTION     ',IERR)

'+G11 T   Vessel Wall Redefinition Option (Baffle Removal)'     0
      call rdi (redefopt,.true.,0,.true.,1,'VESSEL REDEF. OPT    ',IERR)

'+N12 TN1488 Neutral Impurity Velocity Type Option        '     1
      call rdi (cneutvel,.true.,0,.true.,2,'IMP.NEUT.VEL.TYPE.OPT',ierr)

'+N13 T   Neutral Wall Reflection 0-off 1-on              '     0
      ! ammod - added options 3 and 4 to NRFOPT
      CALL RDI (NRFOPT,.TRUE., 0,.TRUE., 4,'NEUTRAL REFLECTION   ',IERR)

'+I06 TN1465 Follow Imp.Ions Recombined to Neutrals 0=off '     1
      CALL RDI (CFOLREC,.TRUE., 0,.TRUE.,1,'FOLLOW REC. NEUTRAL  ',IERR)

'+N14 TN1481 Imp Neutral Momentum.Tran.Coll Opt 0=off 1=on'     1
      call rdi (mtcopt,.true.,0,.true.,2,'NEUT.MOM.TRAN.COLL OPT ',ierr)

'+D09 TN1481 BG Ion Mom.Tran.Coll.Coef.     (kelighi)     '    5.0E-16
      call rdr (kelighi,.true.,0.0,.false.,0.0,'MTC COEFF 1 Imp-I',ierr)

'+D10 TN1481 BG Neutral Mom.Tran.Coll.Coef. (kelighg)     '    5.0E-16
      call rdr (kelighg,.true.,0.0,.false.,0.0,'MTC COEFF 2 Imp-N',ierr)

'+I07 TN1479 Ion Prompt Redeposition Option 0=off 1=on    '     0
      call rdi (prompt_depopt,.true.,0,.true.,1,'ION PROMPT DEPOSITION',ierr)

'+G12 T   Target Position Option 0 to 6                   '     6
      CALL RDI (CTARGOPT,.TRUE.,0,.TRUE.,6,'TARGET POSITION OPT  ',IERR)

'+I08 T   Target Mirror Option 0-off  1-on                '     0
      call rdi (cmiropt, .true.,0,.true.,4,'TARGET MIRROR OPT    ',IERR)

'+G13 T   Pre-defined geometry option -1-off 0-719 1-307  '    -1
      CALL RDI (CGEOOPT,.TRUE.,-1,.TRUE.,1,'GEOMETRY OPTION      ',IERR)

'+I09 T   Ion Periphery Option   0 to 3                   '     1
      CALL RDI (FPOPT,   .TRUE.,0,.TRUE.,6,'ION PERIPHERY OPT    ',IERR)

'+I10 TN996 Periphery Recycle Option       0-off 1-on     '     0
      call rdi (fpropt,.true.,  0,.true.,1,'FP RECYCLE OPT       ',ierr)

'+T13 TN505 Poloidal Velocity Drift Option 0-off 1-on     '     0
      CALL RDI (CPDRFT,.TRUE., 0,.TRUE., 3,'POL. DRIFT OPTION    ',IERR)

'+B01    Special plasma parameter           Rspec        '     1
      CALL RDI (IRSPEC,.FALSE., 0 ,.FALSE., 0 ,'SPEC PLASMA RSPEC',IERR)

      ! move to after input file read
      IF (CIOPTB.EQ.0.AND.(CIOPTC.EQ.0.or.cioptc.eq.4).AND.(CIOPTD.EQ.0.OR.CIOPTD.EQ.3)) IRSPEC=2*MAXNRS

'+G14    Central Mirror Ring Location       (IR)         '     1
      call rdi (ircore,.true.,1,.false.,0,'CORE MIRROR RING SPEC' ,ierr)

'+S02    Mass of plasma ions                Mb           '    2.0
      CALL RDR (CRMB, .TRUE. ,0.1 ,.FALSE.,0.0,'PLASMA ION MASS  ',IERR)

'+S03    Charge on plasma ions              Zb           '    1
      CALL RDI (CIZB, .TRUE. , 1  ,.FALSE., 0 ,'PLASMA ION CHARGE',IERR)
      RIZB = REAL (CIZB)

'+Q06    Temperature of electrons at 0      TeB0  (eV)   '   20.0
      CALL RDR(CTEB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TEB0  ',IERR)

'+Q07    Temperature of electrons at plates TeBP  (eV)   '   20.0
      CALL RDR(CTEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TEBP  ',IERR)

'+Q08    Temperature outer TeB step         TeBout(eV)   '    0.05
      CALL RDR(CTEBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBOUT',IERR)

'+Q09    Temperature inner TeB step         TeBin (eV)   '  250.0
      CALL RDR(CTEBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TEBIN ',IERR)

'+Q10    Temperature of trapped plasma      TeBt  (eV)   '   10.0
      CALL RDR(CTEBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TEBT  ',IERR)

'+Q11 TN1278 Te exp decay step in trap      TeBouP       '    0.01
      CALL RDR(CTEBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TEBOUTP',IERR)

'+Q12    Temperature gradient factor        feBL1        '    0.0
      CALL RDR(CFEBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL1 ',IERR)

'+Q13    Temperature gradient factor        feBL2        '    0.0
      CALL RDR(CFEBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBL2 ',IERR)

'+Q14    Temperature gradient factor        feBt         '    1.0
      CALL RDR(CFEBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEBT  ',IERR)

'+Q15    Temperature gradient factor        feB2         '    1.0
      CALL RDR(CFEB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FEB2  ',IERR)

'+Q16    Temperature of ions at 0           TiB0  (eV)   '   20.0
      CALL RDR(CTIB0 ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP AT 0   TIB0  ',IERR)

'+Q17    Temperature of ions at plates      TiBP  (eV)   '   20.0
      CALL RDR(CTIBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES TEMP TIBP  ',IERR)

'+Q18    Temperature outer TiB step         TiBout(eV)   '    0.05
      CALL RDR(CTIBOU,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBOUT',IERR)

'+Q19    Temperature inner TiB step         TiBin (eV)   '  250.0
      CALL RDR(CTIBIN,.TRUE. ,0.0,.FALSE.,0.0,'TEMP STEP   TIBIN ',IERR)

'+Q20    Temperature of trapped plasma      TiBt  (eV)   '   10.0
      CALL RDR(CTIBT ,.TRUE. ,0.0,.FALSE.,0.0,'TEMP TRAP   TIBT  ',IERR)

'+Q21 TN1278 Ti exp decay step in trap      TiBouP       '    0.01
      CALL RDR(CTIBOUP,.TRUE.,0.0,.FALSE.,0.0,'TEMP DECAY TIBOUTP',IERR)

'+Q22    Temperature gradient factor        fiBL1        '    0.0
      CALL RDR(CFIBL1,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL1 ',IERR)

'+Q23    Temperature gradient factor        fiBL2        '    0.0
      CALL RDR(CFIBL2,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBL2 ',IERR)

'+Q24    Temperature gradient factor        fiBt         '    1.0
      CALL RDR(CFIBT ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIBT  ',IERR)

'+Q25    Temperature gradient factor        fiB2         '    1.0
      CALL RDR(CFIB2 ,.TRUE. ,0.0,.FALSE.,0.0,'GRAD FACTOR FIB2  ',IERR)

'+Q26    Density at 0                       NB0   (m**-3)'    1.0E19
      CALL RDR(CNB0  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS AT 0   NB0   ',IERR)

'+Q27 TN1264 Density at plates              NEBP  (m**-3)'    4.3e19
      CALL RDR(CNEBP ,.TRUE. ,0.0,.FALSE.,0.0,'PLATES DENS NEBP  ',IERR)

'+Q28    Density outer NB step              NBout (m**-3)'    0.03
      CALL RDR(CNBOUT,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBOUT ',IERR)

'+Q29    Density inner NB step              NBin  (m**-3)'    5.0E18
      CALL RDR(CNBIN ,.TRUE. ,0.0,.FALSE.,0.0,'DENS STEP   NBIN  ',IERR)

'+Q30    Density of trapped plasma          NBt   (m**-3)'    1.0E19
      CALL RDR(CNBT  ,.TRUE. ,0.0,.FALSE.,0.0,'DENS TRAP   NBT   ',IERR)

'+Q31 TN1278 Ne exp decay step in trap      NBouP        '    0.01
      CALL RDR(CNBOUP,.TRUE. ,0.0,.FALSE.,0.0,'DENS DECAY  CNBOUP',IERR)


     

'+Q32 TN1347 Langmuir Probe Switch     0=Nb  1=Isat       '     1
!     ENTER LANGMUIR PROBE DATA ALONG PLATES FOR TEMPERATURE GRADIENT
!     OPTIONS 3 AND 4, AND PLASMA DECAY OPTIONS 2 AND 3. IN THE
!     CASE OF THE PLASMA DECAY OPTIONS THE DATA ENTERED HERE WILL
!     BE FILLED IN FOR ALL POINTS OF THE PLASMA. THUS GIVING A
!     UNIFORM CONDITION UNLESS A TEMPERATURE GRADIENT OPTION
!     MODIFIES THE TEMPERATURE.
!
!     FOR THOSE OPTIONS WITH ONLY ONE SET OF DATA THE INNER
!     SPECIFICATION IS USED.
!
!     NOTE LPDATI :- LANGMUIR PROBE DATA INNER ...
!
!     THE LANGMUIR PROBE DATA SWITCH IDENTIFIES THE THIRD QUANTITY
!     AS EITHER THE DENSITY OR THE Isat (PROBE SATURATION CURRENT)
!
!     THE DATA IS ORDERED AS  RING #,TEBP,TIBP,NBP (or ISAT)
!
      CALL RDI(lpdatsw,.TRUE. , 0 ,.TRUE., 2 ,'LP DATA SWITCH', IERR)

'+Q33    OUTER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
      call rdr3(te_mult_i,ti_mult_i,n_mult_i,.TRUE.,0.0,.FALSE.,0.0,'Inner Multipliers       ',IERR)

'+Q34 ' 'Probe data at outer plate (opt4) or total (opt3)'
'    Ring , TeBP , TiBP , NBP : Number of rows - '          0
      CALL RDRARN(LPDATI,NLPDATI,MAXINS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,3,'SET OF L. PROBE DATA INNER',IERR)
!
!     Outer
!

'+Q35    INNER Target Data Multipliers (Te,Ti,Nb):'  1.0   1.0   1.0
      call rdr3(te_mult_o,ti_mult_o,n_mult_o,.TRUE.,0.0,.FALSE.,0.0,'Outer Multipliers       ',IERR)

'+Q36 ' 'Probe data at inner plate          (T grad opt4)'
'    Ring , TeBP , TiBP , NBP : Number of rows - '          0
      CALL RDRARN(LPDATO,NLPDATO,MAXINS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,3,'SET OF L. PROBE DATA OUTER',IERR)


!     move to after input file read 
!      
!     Adjust the Target conditions - by multiplication factors
!
!     jdemod 
!
!     If the lpdat switch is set to 2 convert particles/s to Amperes for
!     compatibility elsewhere in the code. 
!      
      if (lpdatsw.eq.2) then 
!
!        Note: after conversion the data is then in A and the lpdatsw needs to be
!        set accordingly. 
!
         write(6,'(a)') 'INPUT: LPDATSW : Initial value 2 : '//
     >     'Isats converted from particles/s to Amp : LPDATSW set to 1'
!
         lpdat_conv = ech
         lpdatsw = 1
      else
         lpdat_conv = 1.0
      endif

      do in = 1,nlpdati
         lpdati(in,2) =  lpdati(in,2) * te_mult_i
         lpdati(in,3) =  lpdati(in,3) * ti_mult_i
         lpdati(in,4) =  lpdati(in,4) *  n_mult_i * lpdat_conv
      end do

      do in = 1,nlpdato
         lpdato(in,2) =  lpdato(in,2) * te_mult_o
         lpdato(in,3) =  lpdato(in,3) * ti_mult_o
         lpdato(in,4) =  lpdato(in,4) *  n_mult_o * lpdat_conv
      end do

'+Q37 ' 'CORE Plasma Data - for Core Options 1,2 and 3'
'    Ring , TeB , TiB , NB , Vb : Number of rows - '       0
!
!     Read in data for core rings - data varies depending on
!     core option selected.
!
      CALL RDRARN(coredat,Ncoredat,MAXINS,-MACHHI,MACHHI,.FALSE.,-machhi,MACHHI,4,'SET OF CORE DATA',IERR)

'+T14    Cross Field Diffusion factor       Dperp (m*m/s)'    0.3
      CALL RDR(CDPERP,.TRUE.,0.0,.FALSE.,0.0,'X DIFF DPERP       ',IERR)

'+T15    Trap Cross Field Diffusion factor  Dperpt(m*m/s)'   -1.0
      CALL RDR(CDPERPT,.FALSE.,0.0,.FALSE.,0.0,'X DIFF DPERP-TRAP',IERR)

      ! move after input file read
      if (cdperpt.le.0.0) cdperpt = cdperp

'+T16    Perpendicular Pinch Velocity        Vpinch (m/s)'    0.0
      CALL RDR(cVPINCH,.FALSE.,0.0,.FALSE.,0.0,'PERP PINCH VEL.  ',IERR)

'+S04    Mass of impurity ions              Mi           '   12.0
      CALL RDR(CRMI,  .TRUE. ,0.1,.FALSE.,0.0,'IMPURITY ION MASS', IERR)

'+S05    Atomic number of impurity ions     Zi           '      6
      CALL RDI(CION,  .TRUE. , 1 ,.FALSE., 0 ,'IMP ATOMIC NUMBER', IERR)

'+D11    Characteristic energy Ebd          Ebd   (eV)   '    0.0
      CALL RDR(CEBD,  .TRUE. ,0.0,.FALSE.,0.0,'CHAR ENERGY EBD',   IERR)

'+I11    Z effective (self)                 Zeff         '      1
      CALL RDI(CIZEFF,.TRUE. , 0 ,.FALSE., 0 ,'ZEFF(SELF)     ',   IERR)

'+S06    Initial temperature                Tem1  (eV)   '    0.5
      CALL RDR(CTEM1, .FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP   ',   IERR)

'+S07    Initial temperature (2)            Tem2  (eV)   '    0.0
      CALL RDR(CTEM2 ,.FALSE.,0.0,.FALSE.,0.0,'INITIAL TEMP (2)',  IERR)

'+S08    Initial R position of impurity     R0    (m)    '    2.65
      CALL RDR(CXSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL R POSITION',IERR)

'+S09    Initial Z position of impurity     Z0    (m)    '    1.6226
      CALL RDR(CYSC,  .FALSE.,0.0,.FALSE.,0.0,'INITIAL Z POSITION',IERR)

'+I12    Initial ionization state of impurity ions       '      1
      CALL RDI(CIZSC, .TRUE.,  1, .TRUE.,CION,'INITIAL IZ STATE',  IERR)

'+D12    Neutral hydrogen density parameter Nhc   (m**-3)'    1.0E15
      CALL RDR(CNHC,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHC   ',IERR)

'+D13                                       Nho   (m**-3)'    3.0E18
      CALL RDR(CNHO,  .FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  NHO   ',IERR)

'+D14                                       lamhx (m)    '    0.02
      CALL RDR(CLAMHX,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHX ',IERR)

'+D15                                       lamhy (m)    '    0.11
      CALL RDR(CLAMHY,.FALSE.,0.0,.FALSE.,0.0,'H DENSITY:  LAMHY ',IERR)

'+D16    Constant for CX Recomb option 2    Vcx   (m/s)  '    2.4E4
      CALL RDR(CVCX  ,.TRUE., 0.0,.FALSE.,0.0,'CXREC CONSTANT VCX',IERR)

'+B02    For average density "near" target  xnear (m)    '    0.6
      CALL RDR(CXNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  XNEAR ',IERR)

'+B03    For average density "near" target  ynear (m) +/-'    0.6
      CALL RDR(CYNEAR,.TRUE. ,0.0,.FALSE.,0.0,'NEAR PARAM  YNEAR ',IERR)

'+N15    Measure theta from T (degrees to X=0) for launch'  -90.0
      CALL RDR(CSNORM,.FALSE.,0.0,.FALSE.,0.0,'MEASURE THETA FROM',IERR)

      ! move to after input file read
      CSNORM = CSNORM / RADDEG

'+Q38    Inboard plasma flow velocity       Vhin  (m/s)  '    0.0
      CALL RDR(CVHIN, .FALSE.,0.0,.FALSE.,0.0,'INB. PLASMA FLOW V',IERR)

'+Q39    Inboard electric field             Ein   (V/m)  '    0.0
      CALL RDR(CEIN , .FALSE.,0.0,.FALSE.,0.0,'INBOARD ELEC FIELD',IERR)

'+Q40    Outboard plasma flow vel  (SOL5)   Vhout (m/s)  '    0.0
      CALL RDR(CVHOUT,.FALSE.,0.0,.FALSE.,0.0,'OUT. PLASMA FLOW V',IERR)

'+Q41    Outboard electric field   (SOL5)   Eout  (V/m)  '    0.0
      CALL RDR(CEOUT ,.FALSE.,0.0,.FALSE.,0.0,'OUTBRD ELEC FIELD ',IERR)

'+I13    Collision Enhancement Factor       Zenh         '    1.0
      CALL RDR(CZENH, .FALSE.,0.0,.FALSE.,0.0,'COLLIS ENHANC ZENH',IERR)

'+I14    Set Ti = max(Ti,Tb) when reaching state (0 off) '      0
      CALL RDI(CIZSET,.FALSE., 0 ,.FALSE., 0 ,'SET TI=TB AT STATE',IERR)

'+D17    Threshold yield for sput opt3            (eV)   '    0.1
      CALL RDR(CTRESH,.FALSE.,0.0,.FALSE.,0.0,'SELF-SPU THRESHOLD',IERR)

'+D18    Bombarding ion charge state        Zimp         '      0
      CALL RDI(CBOMBZ,.TRUE.,  0 ,.FALSE., 0 ,'BOMBION CHARGE    ',IERR)

'+D19    Bombion type 0Zb 1H 2D 3T 4He4 5C 6Zi 7O        '      0
      CALL RDI(CBOMBF,.TRUE.,  0 ,.TRUE.,  7 ,'BOMBION FLAG 0:7  ',IERR)

'+D20    Ionisation rate factor for neutrals          IRF'    1.0
      CALL RDR(CIRF  ,.TRUE., 0.0,.FALSE.,0.0,'IONISE RATE FACTOR',IERR)

'+D21    Sputtering Enhancement Factor                SEF'    1.0
      CALL RDR(CSEF  ,.TRUE., 0.0,.FALSE.,0.0,'SPUT ENHANCE FACT.',IERR)

'+P06    SOL Enhancement Factor - Electric Field    SOLEF'    1.0
      CALL RDR(CSOLEF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE E(S) ',IERR)

'+P07    SOL Enhancement Factor - Drift Velocity    SOLVF'    1.0
      CALL RDR(CSOLVF,.TRUE., 0.0,.FALSE.,0.0,'SOL  ENHANCE VH(S)',IERR)

'+P08    SOL1a Factor                               fl   '    0.01
      CALL RDR(CFL   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FL" ',IERR)

'+P09    SOL1a Factor                               fs   '    1.0
      CALL RDR(CFS   ,.FALSE.,0.0,.FALSE.,0.0,'SOL1A FACTOR "FS" ',IERR)

'+P10    SOL10 Reversal Mach Number                 fRM  '    1.0
      CALL RDR(CFRM  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "FRM"',IERR)

'+P11    SOL10 factor                               kin  '    1.0
      CALL RDR(CKIN  ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR "KIN"',IERR)

'+P12    SOL10 factor                               kout '    1.2
      CALL RDR(CKOUT ,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACTOR"KOUT"',IERR)

'+P13    SOL10 factor                               fRmin'    0.01
      CALL RDR(CFRMIN,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMIN"',IERR)

'+P14    SOL10 factor                               fRmax'    0.4
      CALL RDR(CFRMAX,.FALSE.,0.0,.FALSE.,0.0,'SOL10 FACT "FRMAX"',IERR)

'+P15    SOL6&7 Vb Length factor 1 (* SMAX)         VbL1 '    0.166
      CALL RDR(CVBL1 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 1',IERR)

'+P16    SOL6&7 Vb multiplication factor 1          VbM1 '    0.0
      CALL RDR(CVBM1 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 1  ',IERR)

'+P17    SOL6&7 Vb Length factor 2 (* SMAX)         VbL2 '    0.5
      CALL RDR(CVBL2 ,.TRUE.,0.0 ,.TRUE. ,1.0,'SOL11 LENGTH 2',IERR)

'+P18    SOL6&7 Vb multiplication factor 2          VbM2 '    0.0
      CALL RDR(CVBM2 ,.false.,0.0,.FALSE.,1.0,'SOL11 MULT 2  ',IERR)

'+S10    Operation Mode  1 Time-Dependent  2 Steady-State'      2
      CALL RDI(IMODE, .TRUE.,  0, .TRUE.,   2  ,'OPERATION MODE',  IERR)

      
!     jdemod - remove upper bounds checks for maxizs and maximp due
!              to dynamical allocation      

'+I15    Maximum ionization state                        '      6
      CALL RDI(NIZS,  .TRUE.,  0, .FALSE.,MAXIZS,'MAX IZ STATE',   IERR)

!
!     jdemod - After NIZS and CION have been read in - a value can be
!              assigned to maxizs and used for storate allocation      
!
!     - moved from allocate_storage_div - move to after input file read
!
      maxizs   = max(nizs,cion)
!
!     Allocate the dwelts array which depends on maxizs and which is read
!     later in the input file - problem?? - allocated arrays have their size
!     in the number of elements entered in the input - just need to avoid referencing
!     maxizs when referring to them - maxnts = nts ?
!      
!     move to after read in?  
      call allocate_mod_dynam4_input_special



      call pr_trace('READIN','AFTER SPECIAL ALLOCATION')

     
'+S11    Number of impurity ions to be followed          '    500
      CALL RDI(NIMPS, .TRUE.,  1, .FALSE.,MAXIMP,'NO OF IONS',     IERR)

'+S12 TN487 Number of Supplementary Neutrals to Launch   '      0
      CALL RDI(NIMPS2,.TRUE.,0,.FALSE.,MAXIMP-NIMPS,'NUM SUP IONS',IERR)
!
!     RESET NIMPS2 - if an Ion injection has been specified and
!                    not a neutral launch - set nimps2 = 0
!
      ! move to after input file read
      if (cneuta.eq.1) nimps2 = 0


'+S13    Quantum iteration time for atoms   fsrate (s)   '    1.0E-8
      CALL RDR(FSRATE,.TRUE.,1.E-10,.TRUE.,1.0,'NEUT ITERATE TIME',IERR)

'+S14    Quantum iteration time for ions    qtim   (s)   '    1.0E-7
      CALL RDR(QTIM,  .TRUE.,1.E-10,.TRUE.,1.0,'DIV ITERATE TIME', IERR)

'+S15 T   CPU time limit in seconds          cpulim (s)   '  80000.0
      CALL RDR(CPULIM,.TRUE.,0.0,.FALSE., 0.0 ,'CPU TIME LIMIT'  , IERR)



      
!---- READ IN TIME DEPENDENT DATA  (NOTE 128)
!
'+S16 ' 'Average Dwell Times (s) for each state 0,1,2..'
'        Number of dwell times given below :-'    0
      CALL RDRAR(DWELTS(0),NQS,MAXIZS+1,0.0,machhi,.FALSE.,'DWELT',IERR)
      IF (IMODE.EQ.1.AND.NQS.LT.NIZS+1) GOTO 1002
      DWELTS(-1) = DWELTS(0)


'+S17 ' 'Dwell Time Factors for time-dependent analysis'
'        Number of dwell time factors :-'  0
      CALL RDRAR(DWELFS,NTS,MAXNTS,  0.0,machhi,.TRUE.,'T FACTORS',IERR)


'+S18 T   Maximum dwell time for steady state (s)         '    0.5
      call rdr(ctimmax,.true.,0.0,.false.,10.0,'DWELL TIME LIMIT ',ierr)
!
!---- READ IN YIELD MODIFIER FUNCTIONS
!


! note default values should be set so that only modifications need to be read in      
'+D22 ' 'Set of Yield Modifiers for Primary, Secondary neutrals'
'      Number of rows of (X,Mpt,Mst,Mct,Mpw,Mcw,Refl) data-'  1
         0.0   0.0    1.0     1.0    1.0    1.0    1.0     1.0
      CALL RDRARN(CYMFS,NYMFS,MAXPTS+1,-machhi,machhi,.true.,-machhi,machhi,7,'SET YIELD(P,S,CT,PW,CW) VALS',IERR)

'+D23 TN? Fixed Yield Value for Sputter Data Option 4     '     0.001
      CALL RDR(const_yield,.true.,0.0,.FALSE.,0.0,'Fixed Yield',IERR)

      !
!---- Read in target plate temperature
!

'+D24 TN1209 Target Temperature (K) for Chem. Sputt. Opt. '   600.0
      CALL RDR(ctargt ,.true.,0.0,.FALSE.,0.0,'Target Temp (K) ',IERR)

'+D25       Main Wall Temperature (K) for Chem. Sputt.   '   453.0
      CALL RDR(cwallt ,.true.,0.0,.FALSE.,0.0,'Wall Temp (K)   ',IERR)

'+D26 TN1450 PP Wall   Temperature (K) for Chem. Sputt.   '   453.0
      CALL RDR(cwallp ,.true.,0.0,.FALSE.,0.0,'PP Wall Temp (K)',IERR)

'+D27 ' 'TN1450 Wall Temperatures (K) for specific segments'
'      Number of segment ranges (Index1 Index2 Temp):'      0
      CALL RDRARN(walltemp,NWLTEMP,MAXPTS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,2,'WALL SEGMENT TEMPERATURES',IERR)

      !
!---- READ IN DEBUG OPTS 0:OFF >0 PRINT DIAGNOSTICS EVERY I'TH TIMESTEP
!---- READ IN RANDOM NUMBER GENERATOR SEED
!
'+B04    Debug atoms (0 off, >0 print every nth timestep)'      0
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG NEUT ***',IERR)
      CSTEPN = REAL (ISTEP)

'+B05    Debug ions  (0 off, >0 print every nth timestep)'      0
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'*** DEBUG DIV  ***',IERR)
      CSTEPL = REAL (ISTEP)

'+B06    Debug ion velocity  (0 off, >1 on)              '      0
      CALL RDI(ISTEP ,.FALSE., 0 ,.FALSE., 0 ,'**DEBUG VELOCITY**',IERR)

      ! move to after input file read
      CSTEPV = REAL (ISTEP)
      ! jdemod - move setting of debug velocity flag to point where it
      ! is read in
      if (cstepv.ne.0.0) then
         debugv = .true.
      else
         debugv = .false.
      endif

'+S19    Random number seed  (0 generate new seed)       '      0
      CALL RDI(CISEED,.FALSE., 0 ,.FALSE., 0 ,'RANDOM NUMBER SEED',IERR)

'+H01    PIN Random number seed  (<0=1, 0 generate new)  '     -1
      CALL RDI(PINISEED,.FALSE.,0,.FALSE.,0,'PIN RANDOM NUM. SEED',IERR)

      ! move to after input file read
      !
!     PIN expects a value of -1 to generate a new seed while the DIVIMP
!     documentation indicates <0 = 1, 0=generate new ... so we need
!     to change the value of piniseed here to match the EIRENE
!     expectations.(Other values are passed along as is). 
!     
      if (piniseed.lt.0) then 
         piniseed = 1
      elseif (piniseed.eq.0) then 
         piniseed = -1
      endif
      

'+A04    Print option  (0 reduced, 1 full)               '     10 3
      CALL RDI(CPRINT,.true., 0 ,.true.,  19 ,'PRINT OPTION',IERR)

'+H02    PIN Data Print option  (0 reduced, 1 more)      '      0
      CALL RDI(PINPRINT,.true., 0 ,.true.,1,'PIN DATA PRINT OPT',IERR)

'+S20    Number of Iterations                            '      1
      CALL RDI(NITERS,.TRUE. , 1 ,.TRUE. ,20 ,'NO. OF ITERATIONS ',IERR)

      
      ! move to after input file read - checks
      IF (NITERS.GT.1 .AND. (CIOPTB.EQ.3.or.cioptb.eq.9)) GOTO 1003
1003  CALL PRC ('READIN: ITERATIONS > 1 INCOMPATIBLE WITH COLL OPT 3')
      IERR = 1

'+I16    Stop following ions reaching Main Plasm 0no 1yes'      0
      CALL RDI(CSTOP ,.TRUE. , 0 ,.TRUE. , 1 ,'STOP WHEN HIT MAIN',IERR)

'+G15    Rectangular grid for neutrals 0calculate 99file '      0
      CALL RDI(CRECT ,.TRUE. , 0 ,.TRUE. , 99,'RECT GRID CALC.   ',IERR)

'+D28    Temperature Gradient Coefficient parameter  ZO  '      4
      CALL RDI(CZO   ,.TRUE. , 0 ,.FALSE., 0 ,'ZO PARAMETER      ',IERR)

      ! move to after input file read
      ZO = REAL(CZO)
      CHZO = 1.56*(1.+1.41*ZO)*(1.+0.52*ZO)/((1.+2.65*ZO)*(1.+0.285*ZO))


'+F03    Plasma condition for missing SOL rings     CNIWA'      1
      CALL RDI(CNIWA ,.TRUE. , 0 ,.FALSE., 1 ,'CNIWA OPTION      ',IERR)

'+F04 ' 'EDGE1D/2D Deuterium drift vel. mult. factor VMF '
'            Number of VMF blocks                          '      0
'            Ring Range :-' -20     -30
'            J0 & J1    :-'   5       5
'+           VMF0,1,2   :-'   1.000   1.000   1.000
      CALL RDVMF('VEL. MULT. FACTOR',IERR)

'+I17    Ion removal loss time              Tloss  (s)   '    0.000
      CALL RDR(TLOSS,.FALSE.,0.0,.FALSE.,0.0,'ION LOSS TIME', IERR)

'+P19    Power density                      P/A    (W/m2)'    3.0E+07
      CALL RDR(CPA  ,.FALSE.,0.0,.FALSE.,0.0,'POWER DENSITY', IERR)

'+P20    Parallel heat conduction coeff     K0           '    2.0E+03
      CALL RDR(CK0  ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND', IERR)

'+P21    Parallel heat conduction -ions     K0I          '     58.9
      CALL RDR(CK0I ,.FALSE.,0.0,.FALSE.,0.0,'PAR HEAT COND IONS',IERR)

'+P22    Override input E-field from file E=0  0-off 1-on'      3
      CALL RDI(OFIELD,.TRUE. , 0 ,.TRUE. , 5 ,'EFIELD=0 FOR FILE',IERR)

'+P23 T1433 E-field Opt 4 - Length of E-field region *SMAX'    0.25
      CALL RDR(CEFLEN,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD SRC LENGTH',IERR)

'+P24 T1433 E-field Opt 4 - Te collisionality multiplier  '    2.0
      CALL RDR(CEFfact,.TRUE. ,0.0,.FALSE.,0.0,'EFIELD COLL-MULT',IERR)

'+P25 TN401 Decay length for ionization source Ls  SOL12  '    0.08
      CALL RDR(CSOLLS,.TRUE.,0.0,.FALSE.,0.0, 'LS DECAY SOL12   ',IERR)

'+P26 T     Second decay characteristic length            '    0.08
      CALL RDR(CSOLLT,.TRUE.,0.0,.FALSE.,0.0, 'SECOND DECAY DIST',IERR)

'+P27 T     Source fraction                               '    0.5
      CALL RDR(CFIZ  ,.TRUE.,0.0,.TRUE. ,1.0, 'SOURCE FRACTION  ',IERR)

'+P28 TN401 Decay length for radiative losses  Lr  SOL12  '    0.08
      CALL RDR(CSOLLR,.TRUE.,0.0,.FALSE.,0.0, 'LR DECAY SOL12   ',IERR)

'+P29 TN401 Coefficient for radiative losses   Pr  SOL12  '    0.0
      CALL RDR(CSOLPR,.TRUE.,0.0,.FALSE.,0.0, 'PR CONST SOL12   ',IERR)

'+P30 TN775 Radiation source strength fraction Frr        '    1.0
      CALL RDR(CSOLFR,.TRUE.,0.0,.FALSE.,0.0, 'FR SOURCE FRAC.  ',IERR)

'+P31 TN401 Source Ionization Option 0-lin 1-exp   SOL12  '    1
      CALL RDI(CSOPT ,.TRUE., 0,.TRUE.,5, 'IONIZATION SOURCE OPT',IERR)

'+P32 TN401 Source Radiation Option  0-lin 1-exp   SOL12  '    3
      CALL RDI(CPOPT ,.TRUE., 0,.TRUE.,3, 'RADIATION SOURCE OPT ',IERR)

'+P33 TN660 Imaginary Root Option                  SOL12+ '    1
      CALL RDI(SROOTOPT,.TRUE.,0,.TRUE.,1,'TREAT NEGATIVE ROOTS ',IERR)

'+P34 TN407 Flux Recirculation Option 0-off 1-on          '    0
      CALL RDI(FLUXROPT,.TRUE.,0,.TRUE.,1,'SRC FLUX RECIRC OPT',  IERR)

'+P35 ' 'TN????+407 Set of flux specifications  '
'    Ring , data1(I), data2(O), data3 : Rows - '          0
      CALL RDRARN(FLUXINFO,FLUXPTS,MAXINS,-MACHHI,MACHHI,.FALSE.,-machhi,MACHHI,3,'SET OF FLUX DATA',IERR)

'+B07 TN442 Z-value defining divertor region  (m)         '   1.7
      CALL RDR(CZD   ,.FALSE.,0.0,.FALSE.,0.0, 'Z DIVERTOR LIMIT',IERR)

'+I18 TN480 Ring for ion injection option 2        INJIR  '    1
      CALL RDI(INJIR,.TRUE.,1,.TRUE.,MAXNRS,'RING NUMBER FOR INJ',IERR)

'+I19 TN480 Factor for Region Lower Bound          INJF1  '   0.0
      CALL RDR(INJF1,.TRUE.,0.0,.TRUE.,1.0  ,'INJECTION AREA LB' ,IERR)

'+I20 TN480 Factor for Region Upper Bound          INJF2  '   0.0
      CALL RDR(INJF2,.TRUE.,INJF1,.TRUE.,1.0,'INJECTION AREA UB' ,IERR)

'+I21 TN443 X-max for Far-Periphery region (O/I) '       0.1  0.1
      CALL RDR2(FPXMAXO,fpxmaxI,.TRUE.,0.0,.FALSE.,0.0,'FP DIFFUSION SIZE' ,IERR)

'+I22 TN443 Far-Periphery Target Loss Time (O/I) '    1.0e-3  1.0e-3
      CALL RDR2(FPTIMO,fptimI,.TRUE.,0.0,.FALSE.,0.0,'FP TARGET LOSS T ' ,IERR)

'+I23 TN688 Far-Periphery Diffusion Rate ( < 0 = CDPERP ) '  -1.0
      CALL RDR(CDPERPFP,.FALSE.,0.0,.FALSE.,0.0,'FP DIFF RATE  ' ,IERR)
      IF (CDPERPFP.LT.0.0) CDPERPFP = CDPERP

'+N16 ' 'TN487 Launch probability modifiers for each     '
'    TN487 Wall segment range  #1  #2  mod.  :- '           0
      CALL RDRARN(WLPROB,NWLPROB,MAXPTS,-MACHHI,MACHHI,.FALSE.,0.0,MACHHI,2,'WALL LAUNCH PROB. MODIFIERS',IERR)

'+N17 TN721 Use Wall Probabilities as Absolute 0=No 1=Yes '    0
      CALL RDI(WLPABS,.TRUE. ,0  ,.TRUE. , 3 ,'ABS WALL PROB    ',IERR)

'+T17 TN505 Poloidal Drift Velocity (m/s)                 '   0.0
      CALL RDR(CDRFTV,.FALSE.,0.0,.FALSE.,0.0,'POL. DRIFT VEL.  ',IERR)

'+T18 TN    Poloidal Drift Range F1*SMAX < S < F2*SMAX'    0.02 0.98
      call rdr2(cdrftv_start,cdrftv_end,.false.,0.0,.false.,1.0,'Start and End of DRFTV',ierr)

'+D29      CEMAXF factor for sputtering velocities       '   1.0
      CALL RDR(CEMAXF,.FALSE.,0.0,.FALSE.,0.0,'EMAX-FACTOR',      IERR)

'+D30 TN521 Impact Energy for wall launch Vel. dist (eV)  '  100.0
      CALL RDR(CEIMP, .TRUE., 0.0,.FALSE.,0.0,'EIMP- WALL LAUNCH',IERR)

'+H03 TN408 Run PIN from inside DIVIMP  0-NO 1-YES        '    1
      CALL RDI(CPINOPT,.TRUE.,0  ,.TRUE. ,4  ,'RUNPIN 0-NO 1-YES',IERR)

'+H04 ' 'TN408 Pin: reire07                              '
      CALL RDC(CPINCOM,'COMMAND TO RUN PIN',IERR)
      READ(CPINCOM(11:80),'(A69)') ACTPIN

'+H05      PIN Cell Area Option (IHCORR)                 '    1
      CALL RDI(IHCORR,.TRUE.,0  ,.TRUE.,1 ,'PIN CELL AREA OPTION',IERR)

'+H06      PIN Hybrid Wall Option 0=off 1,2=selection    '    0
      CALL RDI(IHYBRID,.TRUE.,0  ,.TRUE. ,6 ,'HYBRID WALL IN PIN',IERR)

'+H07      PIN Puffing Option - 0=off 1=on               '    2
      call rdi(pinpuff,.true.,0,.true.,2,    'PIN PUFFING OPTION',ierr)

'+H08      PIN Puff Location switch - 0=main SOL 1=PP    '    0
      call rdi(swpvhpf,.true.,0,.true.,1,  'PUFF LOCATION OPTION',ierr)

'+H09      PIN Puff fraction (opt 1)                     '   0.0
      CALL RDR(hpcpuf, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF FRAC',IERR)

'+H10      PIN Recycle -> Puff fraction (puff opt 2)     '   0.16
      CALL RDR(ppcpuf,.true.,0.0,.FALSE.,0.0, 'FLUX-> PUFF OPT 2',IERR)

'+H11      PIN Puff Injection temperature (eV)           '   0.5
      CALL RDR(tpufh, .TRUE., 0.0,.FALSE.,0.0,'PIN RE-PUFF TEMP ',IERR)

'+H12      PIN Puff location indices JHPUF1(1 and 2)'  -17 -1000
      call rdi2(jhpuf1(1),jhpuf1(2),.false.,0,.false.,1,'PUFF LOCATION INDICES 1',ierr)

'+H13      PIN Puff location indices JHPUF2(1 and 2)'  -16 -1001
      call rdi2(jhpuf2(1),jhpuf2(2),.false.,0,.false.,1,'PUFF LOCATION INDICES 2',ierr)

'+P36 TN408 Calculate SOL iteratively? 0-NO 1-YES         '    1
      CALL RDI(CITERSOL,.TRUE.,0 ,.TRUE. ,2  ,'DO SOL AFTER PIN ',IERR)

'+P37 TN408 Secondary SOL option for iterative calculation'   -2
      CALL RDI(CSECSOL ,.TRUE.,-2,.TRUE. ,14 ,'SECONDARY SOL OPT',IERR)
      IF (CSECSOL.EQ.-2) CSECSOL = CIOPTF

'+P38 TN408 Ionization option for iterative SOL           '    3
      CALL RDI(CSOPT2,.TRUE.,0,.TRUE.,5,  'SECOND ION SOURCE OPT',IERR)

'+P39      Number of Pin/SOL iterations                  '    5 
      CALL RDI(NITERSOL,.TRUE.,0 ,.FALSE.,maxpiniter,'NO. OF PIN ITER. ',IERR)

'+G16 ' 'TN    Set of Plate coordinates                  '
'    TN    Ring #, Outer R,Z   , Inner R,Z :      '       0
      CALL RDRARN(PLATCO,NPLAT,MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,4,'PLATE COORDINATES',IERR)

'+G17 ' 'Wall coordinates                                '
'    TN    R,Z coordinates starting at outer plate   '   0
      CALL RDRARN(WALLCO,NWALL,MAXPTS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,1,'WALL COORDINATES',IERR)

'+G18 ' 'Wall coordinates - PFZ                          '
'    TN    R,Z coordinates                           '  0   
      CALL RDRARN(WALLCO2,NWALL2,MAXPTS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,1,'WALL COORDINATES',IERR)

'+D31 TN83? Maximum Number of sputtered generations       '   75
      CALL RDI(CMAXGENS,.TRUE., 0,.FALSE.,  0, 'MAX. GENERATIONS',IERR)

'+B08 TN    Ring number for detailed background data      '    8
      CALL RDI(CIRHR,  .TRUE., 1,.TRUE.,MAXNRS,'HI-RES RING     ',IERR)

'+N18 A18   Power of cosine release dist. (V/A 12,13)     '   1.0
      CALL RDR(CNIN,.TRUE.,0.01,.FALSE.,0.0,'COSINE DIST. POWER',IERR)

'+N19 TN???? Extra Velocity Multiplier for V/A 14         '   1.0
      CALL RDR(CVAMULT,.TRUE.,LO,.FALSE.,0.0,'VELOCITY MULTIPLIER',IERR)

'+N20 TN???? Velocity Multiplier for Recombined Neutrals  '   1.0
      CALL RDR(CVRMULT,.TRUE.,LO,.FALSE.,0.0,'REC. VEL. MULT',IERR)

'+C01 ' 'Set of S-distances for ion leakage diagnostic(m)'
'TN982     Number of S-values  :-'  0
! Example
!'TN982     Number of S-values  :-'  4
!   5.0
!   10.0
!   15.0
!   20.0
      CALL RDRAR(cleaks,cleaksn,maxpts,0.0,machhi,.TRUE.,'LEAK- S',IERR)

      ! move to after input file read
      !
!     Set up the flag indicating that leak checking code should be activ
!
      if (cleaksn.gt.0) then
         checkleak = .true.
      else
         checkleak = .false.
      endif
!
!     The following variables all apply to SOL option 21 which attempts
!     to model a detached plasma divertor configuration.
!
'+R01      DETACHED PLASMA: Length Scaling Switch 0-S 1-P'    0
      call rdi(s21refsw,.true.,0,.true.,3,'S21 Length Ref Switch',ierr)
!
!     Load up Outer/Inner values
!
'+R02 TN988 Detached Plasma Model: Te/Te0 at L1 O/I '  1.0    1.0
      call rdr2(terat,terati,.true.,0.0,.false.,0.0,'Te Ratio at L1 O/I',ierr)

'+R03 TN1439   Ti/Ti0 at Position L1            O/I '  1.0    1.0
      call rdr2(tirat,tirati,.true.,0.0,.false.,0.0,'Ti Ratio at L1 O/I',ierr)

'+R04 TN988    Ne/Ne0 at Position L1            O/I ' 10.0   10.0
      call rdr2(nrat,nrati,.false.,0.0,.false.,0.0,'Ne Ratio at L1 O/I',ierr)

'+R05 TN1496   Exponent for Ne Equation 1.0=lin O/I '  1.0    1.0
      call rdr2(nalph,nalphi,.true.,0.0,.false.,0.0,'Ne Exponent at L1 O/I',ierr)

'+R06 TN988    Qrad/Q0                          O/I '  5.0    5.0
      call rdr2(qrat,qrati, .false.,0.0,.false.,0.0,'Emitted energy rat',ierr)

'+R07 TN988    L1/SMAX ratio                    O/I '  0.1    0.1
      call rdr2(l1rat,l1rati,.true.,0.0,.false.,0.0,'L1*SMAX boundary 1',ierr)

'+R08 TN988    L2/SMAX ratio                    O/I '  0.2    0.2
      call rdr2(l2rat,l2rati,.true.,0.0,.false.,0.0,'L2*SMAX boundary 2',ierr)

'+R09 TN988    LV/SMAX ratio                    O/I '  0.2    0.2
      call rdr2(lvrat,lvrati,.true.,0.0,.false.,0.0,'V->0 at lvrat*SMAX',ierr)

'+R10 TN       Velocity multiplier at L1        O/I '  1.0    1.0
      call rdr2(vbmult,vbmulti,.true.,0.0,.false.,0.0,'Vmult for Region A',ierr)

'+R11  ' 'TN    DETACHED Plasma Specifications on a by ring basis'
'TN    INNER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
      CALL RDRARN(S21PARMI,NS21I,MAXNRS,-machhi,machhi,.FALSE.,-MACHHI,MACHHI,9,'SOL21 PARAMS INNER',IERR)

'+R12 ' 'TN    DETACHED Plasma Specifications on a by ring basis'
'TN    OUTER - IR TER TIR NR QR L1R L2R LVR VBM - N= '    0
      CALL RDRARN(S21PARMO,NS21O,MAXNRS,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,9,'SOL21 PARAMS OUTER',IERR)

'+D32 TN1007 ABSFAC or Power - Specified - Use if > 0.0   '   0.0
      call rdr(nabsfac,.false.,0.0,.false.,0.0, 'ABSFAC modifier',ierr)


'+G19 ASDEX U - GRID CHARACTERISTICS:  Number of Rings    '    26
!
!     NOTE: These items are now usually on the GEOM: line at the beginning
!           of the grid file or are calculated automatically from the grid 
!     The following set of numbers specify the characteristic values of
!     the ASDEX U style grid file to be read in. The Number of rings
!     and elements per ring (a ring is a set of polygons in a row or
!     indexed by the same ring number), the cutring,(which specifies the
!     end of  rings for the core and trapped plasma), and the cutpts 1
!     and 2 which specify the points on the rings numbered less than the
!     cut ring where the splits for core and trapped or private plasma
!     occur.
!
!     not sure these are used anymore at all
      CALL RDI(maxrings,.TRUE.,0  ,.true.,maxnrs,'MAXRINGS in AUG',IERR)

'+G20                                 Number of Knots    '    34
      CALL RDI(maxkpts,.TRUE.,0  ,.true.,maxnks,'MAX PTS in AUG',IERR)

'+G21                                 Cut ring           '     7
      CALL RDI(cutring,.TRUE.,0  ,.true.,maxrings,'CUTRING in AUG',IERR)

'+G22                                 Cut point 1        '     1
      CALL RDI(cutpt1,.TRUE.,0  ,.true.,maxkpts+1,'CUTPT1 in AUG',IERR)

'+G33                                 Cut point 2        '    34
      CALL RDI(cutpt2,.TRUE.,0  ,.true.,maxkpts+1 ,'CUTPT2 in AUG',IERR)

'+F05 SONNET-Number of Fluids in B2 Background Plasma File'     7
      CALL RDI(nfla,.TRUE.,1,.true.,maxnfla,'NUM.FLUIDS IN B2 BG',IERR)

'+F06       Read Aux. Background Input File    0=off 1=on'     0
      CALL RDI(readaux,.TRUE.,0,.true.,3,'READ AUX BG FILE',IERR)

'+D33  TN1200 Stgrad - Distance where Tgrad forces -> 0    '   1.0
      CALL RDR(Cstgrad,.TRUE.,0.0,.TRUE.,1.0, 'TGRAD FORCES -> 0',IERR)

'+S21    SOLTEST - 0.0 run normally -1.0 test SOL opt    '   -1.0
      call rdr(ctestsol,.true.,-1.0,.false.,1.0,'TEST SOL OPT ONLY',ierr)

'+D34    H Recombination Calculation Option 4-Adas 0+-oth'     4
      call rdi(crecopt,.TRUE.,0  ,.true. ,4  ,'Recomb. Calc Opt',IERR)

'+D35    Recombination Limit Cut-OFF Temperature (eV)    '    0.0
      call rdr(treccut,.TRUE.,0.0,.false.,0.0,'Recomb. cutoff T',IERR)

'+C02 TN1303 Dperp Extractor - methods used - 0 1 2       '     2
      CALL RDI(dpmethod,.TRUE.,0  ,.true.,2 ,'DPERP EXT METHOD',IERR)

'+C03 TN1310 Dperp Ext. 0=only to Xpoint 1=Full Field Line'     1
      call rdi(dpsuml,.TRUE.,0  ,.true.,1 ,'DPERP EXT SUM LIMIT',IERR)

'+C04 TN1310 Dperp Ext. Outer Ring Losses    0=off 1=on   '     1
      call rdi(dpouter,.TRUE.,0  ,.true.,1 ,'DP EXT OUTER RING',IERR)

'+C05 TN1309 Dperp Ext. Dperp  Convection    0=off 1=on   '     1
      call rdi(dpconv,.TRUE.,0  ,.true.,1 ,'DP EXT CONVECT LOSS',IERR)

'+C06 TN1311 Dperp Ext. 1/2 Cell Flux Corr.  0=off 1=on   '     0
      call rdi(dpfluxopt,.TRUE.,0  ,.true.,1 ,'CELL CENTRE FLUX',IERR)

'+C07 TN1311 Dperp Ext. Calc. Average Dperp. 0=off N=on   '     1
      call rdi(dpavopt,.TRUE.,0  ,.false.,1 ,'AVERAGE DPERP OPT',IERR)

'+C08 TN1311 Dperp Ext. Major Radius Corr.   0=off 1=on   '     0
      call rdi(dprcopt,.TRUE.,0  ,.true.,2 ,'MAJOR RADIUS CORR',IERR)

'+C09 TN1314 Dperp Ext. Gradient Smoothing   0=off N=on   '     0
      call rdi(dpsmooth,.true.,0,.false.,0 ,'GRADIENT SMOOTHING',IERR)

'+C10 TN     Dperp Ext. Gradient Calc Meth     -1,0,1     '     0
      call rdi(dpnav,.true.,-1,.true.,2 ,   'GRADIENT CALC METH',IERR)

'+C11 TN     Dperp Ext. Cross-field Area  0-centre 1-bound'     0
      call rdi(dparea,.TRUE.,0  ,.true.,1 ,'CELL BOUND AREA',IERR)

'+C12 TN     Dperp Ext. Power Loss Terms     0-off 1-on   '     1
      call rdi(dpploss,.TRUE.,0  ,.true.,2 ,'POWER LOSS TERMS',IERR)

'+C13 TN     Dperp Ext. Non-ortho Correction 0-off 1-on   '     1
      call rdi(dporth,.TRUE.,0 ,.true.,1,'DP NON-ORTH GRADIENT',IERR)

'+C14 TN     Dperp Ext. Pei Correction Factor             '    1.0
      call rdr(dppei,.TRUE.,0.0,.false.,0.0,'PEI Correction',IERR)

'+C14 TN1373 Dperp Ext. Source Recycle Correction Factor  '    1.0
      call rdr(dprec,.TRUE.,0.0,.false.,0.0,'EXT Recycle Frac',ierr)

'+C16 TN1445 Dperp Ext. Dperp/Xperp Fixed Ratio  0.0=off  '    0.4
      call rdr(dpxpratio,.FALSE.,0.0,.false.,0.0,'Dp/Xp Ratio',ierr)


'+P40 TN     Te S1 - First  S -value = ctes1 * SMAX       '    0.05
!
!     The following values apply to the specifiable SOL option that
!     allows two part linearly interpolated fits for each of
!     Te, Ti, ne and vb to be specified.
!
!
! The following lines specify the parameters for the linearly
! interpolated specified BG SOL option. This is purely empirical.
!                    S-value      Function Value
! The form is:         0.0            F0
!                       S1            F1
!                       S2            F2
! For S>S2 F=F2
!
      call rdr(ctes1,.TRUE.,0.0,.false.,0.0,'Te distance  1',ierr)

'+P41 TN     Te F1 - First  Te-value = ctef1 * te0        '    1.5
      call rdr(ctef1,.TRUE.,0.0,.false.,0.0,'Te fact/mult 1',ierr)

'+P42 TN     Te S2 - Second S -value = ctes2 * SMAX       '    0.3
      call rdr(ctes2,.TRUE.,0.0,.false.,0.0,'Te distance  2',ierr)

'+P43 TN     Te F2 - Second Te-value = ctef2 * te0        '    2.0
      call rdr(ctef2,.TRUE.,0.0,.false.,0.0,'Te fact/mult 2',ierr)

'+P44 TN     Ti S1 - First  S -value = ctis1 * SMAX       '    0.05
      call rdr(ctis1,.TRUE.,0.0,.false.,0.0,'Ti distance  1',ierr)

'+P45 TN     Ti F1 - First  Te-value = ctif1 * ti0        '    1.5
      call rdr(ctif1,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 1',ierr)

'+P46 TN     Ti S2 - Second S -value = ctis2 * SMAX       '    0.3
      call rdr(ctis2,.TRUE.,0.0,.false.,0.0,'Ti distance  2',ierr)

'+P47 TN     Ti F2 - Second Te-value = ctif2 * ti0        '    2.0
      call rdr(ctif2,.TRUE.,0.0,.false.,0.0,'Ti fact/mult 2',ierr)

'+P48 TN     Nb S1 - First  S -value = cnes1 * SMAX       '    0.05
      call rdr(cnes1,.TRUE.,0.0,.false.,0.0,'Nb distance  1',ierr)

'+P49 TN     Nb F1 - First  Te-value = cnef1 * ne0        '    2.0
      call rdr(cnef1,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 1',ierr)

'+P50 TN     Nb S2 - Second S -value = cnes2 * SMAX       '    0.35
      call rdr(cnes2,.TRUE.,0.0,.false.,0.0,'Nb distance  2',ierr)

'+P51 TN     Nb F2 - Second Te-value = cnef2 * ne0        '    2.0
      call rdr(cnef2,.TRUE.,0.0,.false.,0.0,'Nb fact/mult 2',ierr)

'+P52 TN     vb S1 - First  S -value = cvbs1 * SMAX       '    0.1
      call rdr(cvbs1,.TRUE.,0.0,.false.,0.0,'vb distance  1',ierr)

'+P53 TN     vb F1 - First  Te-value = cvbf1 * ti0        '    0.25
      call rdr(cvbf1,.TRUE.,0.0,.false.,0.0,'vb fact/mult 1',ierr)

'+P54 TN     vb S2 - Second S -value = cvbs2 * SMAX       '    0.4
      call rdr(cvbs2,.TRUE.,0.0,.false.,0.0,'vb distance  2',ierr)

'+P55 TN     vb F2 - Second Te-value = cvbf2 * ti0        '    0.0
      call rdr(cvbf2,.TRUE.,0.0,.false.,0.0,'vb fact/mult 2',ierr)

'+C17 Vertical   Reciprocating Probe - Intersection #  '     1
!
!     Reciprocating/Fast Scanning Probe - R location.
!
!     .... and horizontal probe location
!
      call rdi(rlocnum,.true.,1,.false.,0, 'Crossing number',ierr)

'+C18 Vertical   Reciprocating Probe - R-Value         '    3.25
      call rdr(crploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)

'+C19   Horizontal Reciprocating Probe - Intersection #  '     1
      call rdi(zlocnum,.true.,1,.false.,0, 'Crossing number',ierr)

'+C20   Horizontal Reciprocating Probe - Z-Value         '    0.0
      call rdr(czploc,.FALSE.,0.0,.false.,0.0,'Rec.probe loc',ierr)

      
'+P56 TN1424 Core Option4,5- Velocity decay fraction *SMAX'    0.05
!
!     Marfe Core options - for option 4,5,6
!
      call rdr(corefv,.true.,0.0,.true.,0.5,'CORE-VEL FRAC',ierr)

'+P57 TN1424 Core Option4,5- Te,Ti    decay fraction *SMAX'    0.05
      call rdr(coreft,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC',ierr)

'+P58 TN1427 Core Option4,5- Velocity decay frac 2   *SMAX'    0.5
      call rdr(corefv2,.true.,0.0,.true.,0.5,'CORE-VEL FRAC2',ierr)

'+P59 TN1427 Core Option4,5- Te,Ti    decay frac 2   *SMAX'    0.1
      call rdr(coreft2,.true.,0.0,.true.,0.5,'CORE-TEMP FRAC2',ierr)
!
!     Factor for temperature gradient force modifier option
!
'+D36 TN1429 T-GRAD Force Modification - Factor Applied   '    1.0
      call rdr(fgradfact,.false.,0.0,.false.,0.0,'FGRAD MOD FACT',ierr)
!
!     TEMPORARY - CORRECT Tgrad forces for Edge2D bug.
!
!      CALL RDI(fixtgrad,.TRUE.,0  ,.true.,1 ,'FIX E2D TGRAD',IERR)
!
!
!     Read in the INPUT options for SOL option 23
!
      call pr_trace('READIN','BEFORE SOL23 INPUT')
      
'+300 SOL23  Parameter read option  0-NO 1-YES            '     0
      call rdi(readin_sol23_params,.TRUE.,0 ,.true.,1,'READ SOL23 PARMAS OPT',IERR)
!
! Code not needed for unstructured input
!      
!      if (readin_sol23_params.eq.1) then
!
!         call read_sol23_params(ierr)
!
!      endif
!
!     READ IN THE REST OF THE DATA FILE - MAXIMUM 50 LINES
!
!     The reason for this is so that the selected PIN input options
!     can be printed in the DIVIMP output file. At this time we do
!     not want to try interpreting them - simply echo then so that
!     it is clear which options were in use when a particular DIVIMP
!     case was run.
!
      ! NIMBUS input namelist - not used - from &NIMBIN to &END if found in the input file
      ! NIMBUS reads the input file directory to get this info - so reading it in DIVIMP is
      ! purely for documentation purposes
!+H14 Nimbus input namelist: (see, for example, PINPGX for definitions)

      CALL RDCAR(CNIMBIN,nlines,maxlines,'NIMBUS NAMELIST INPUT',IERR)




! move to after input file read in
!
!--------------------------------------------------------------------------
! 
!     CHECKS AND VERIFICATION OF SOME INPUT VALUES:
!
!--------------------------------------------------------------------------   
!
!
!     SOME INPUT DATA VERIFICATION
!
!
!     If the target option specified requires the platco array to
!     store data (as in target option 1), then the geometry data
!     can not be loaded on top.
!
!    Note: if the target option is 3 or the wall option 4 or 5
!           and cgeoopt has not been specified then an error
!           has occurred.
!
      IF (CGEOOPT.EQ.-1.AND.(ctargopt.eq.3.or.
     >    cneur.eq.3)) then
         call prc('Invalid Input ... target or wall option that')
         call prc('requires pre-loaded data has been specified ')
         call prc('but specific pre-loaded grid data has not   ')
         call prc('selected.')
         stop
      endif
!
!     Check to make sure that Trap option 3 and a double-null ITER
!     grid are not specified at the same time - these options share
!     a common array with different meanings at this time and are thus
!     incompatible.
!
      if (cgridopt.eq.2.and.ctrap.eq.3) then
         call prc('ITER grid is incompatible with TRAP option 3.')
         call prc('Please respcify input - program ending ...')
         stop
      endif
!
!     Error - make sure PIN is running if NIMBUS wall specified.
!
      if (cpinopt.eq.0.and.(cneur.eq.5.or.ctrap.eq.5)) then
!
         call prc('ERROR: Incompatible NEUTRAL Wall/Trap Wall Options')
         call prc('Neutral Wall Option 5 or Trap Wall Option 5')
         call prc('specified without corresponding call to PIN/NIMBUS')
         call prc('Program stopping')
!
!        Write to screen for easy notification
!
         write(0,*) 'ERROR: Incompatible NEUTRAL Wall/Trap Wall Options'
         write(0,*) 'Neutral Wall Option 5 or Trap Wall Option 5'
         write(0,*) 'specified without corresponding call to PIN/NIMBUS'
         write(0,*) 'Program stopping'
         stop
      endif
!
!     Warning - if one of the SOL - Plasma Decay or Temperature Gradient
!               options is specified as 99 and at least one is NOT
!             - issue a Warning!
!
      if (cioptg.eq.99.or.cioptf.eq.99.or.cioptk.eq.99.or.cioptl.eq.99) then

         if (cioptg.ne.99.or.cioptf.ne.99.or.cioptk.ne.99.or.cioptl.ne.99) then
!
!           Issue Warning - since mixed file data specified.
!
            call prb
            call prc(' WARNING --- WARNING --- WARNING !!! ')
            call prb
            call prc(' POSSIBLY INCONSISTENT INPUT: ')
            call pri(' SOL OPTION          : ',cioptf)
            call pri(' PLASMA DECAY OPTION : ',cioptg)
            call pri(' Te GRADIENT OPTION  : ',cioptk)
            call pri(' Ti GRADIENT OPTION  : ',cioptl)
            call prb
            call prc(' AT LEAST ONE HAS BEEN SPECIFIED AS FILE INPUT (99)')
            call prc(' BUT NOT ALL - THIS MAY PRODUCE UNEXPECTED RESULTS!')
            call prb
         endif
      endif
!
!     e2dtargopt=2 requires flux data for each flux tube to be specified.
!     e2dtargopt=3,4
!
      if ((e2dtargopt.eq.2.or.e2dtargopt.eq.3.or.e2dtargopt.eq.4.or.e2dtargopt.eq.5)&
              .and.(fluxpts.le.0.and.readaux.ne.1)) then
         call prc('ERROR: Incompatible Input.')
         call prc('- Fluid code target condition option 2, 3 or 4 has been')
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
!
!     Error - CIOPTI=8 is invalid except for Carbon impurity in a
!             deuterium plasma.
!
      if ((ciopti.eq.8.or.ciopti.eq.9).and.(crmb.ne.2.0.or.crmi.ne.12.0)) then
         call errmsg('INPUT ERROR:','Incompatible Input: CX Recombination option 8 is ONLY compatible with a carbon impurity in a deuterium plasma.')
         stop 'Invalid CX option'
       endif
!
!      IF Pinch option 4 is specified then a Probability Distribution
!      Function must be loaded - if not - the pinch option is turned 
!      OFF. 
!
       if (pinchopt.eq.4.and.pinch_npdf.le.0) then 
!
          call errmsg('INPUT ERROR:','PINCH OPTION 4 SELECTED BUT PDF NOT SPECIFIED IN INPUT :'//&
                      ' PINCH OPTION TURNED OFF (SET TO 0)')
          pinchopt = 0
       endif

!     Issue an error message if fp_neut_opt is non-zero for
!     any periphery option except 3.        
!     
       if (fp_neut_opt.ne.0.and.fpopt.ne.3) then
           write(0,'(a,i8,a,i8)')  'WARNING (ERROR): PERIPHERY IONIZATION OPTION ',&
               fp_neut_opt,' HAS BEEN SPECIFIED WITH AN INCOMPATIBLE PERIPHERY OPTION',fpopt
           write(6,'(a,i8,a,i8)') 'WARNING (ERROR): PERIPHERY IONIZATION OPTION ',&
               fp_neut_opt,' HAS BEEN SPECIFIED WITH AN INCOMPATIBLE PERIPHERY OPTION',fpopt
       endif


!       
!      CSTMAX is set in the rundiv main program
!
!      CSTMAX = 10.0 / QTIM
!
!
!
       if (neut2d_opt.eq.2.and.ero_particle_launch_opt.eq.0) then 
         call prc('ERROR: Incompatible Input.')
         call prc('       Neut2d option = 2 specified without')
         call prc('       ERO particle launch being selected')
         call prc('       K40 - ERO particle launch option')
         call prc('       must be set and ERO particle data')
         call prc('       file must be available')
         neut2d_opt = 0

         write(0,'(a)') 'WARNING: NEUT2D option set to ERO'//&
           ' but ERO particle launch option (K40) turned'//&
           ' off. NEUT2D OPT set to 0'

       endif

        ! jdemod - hc_lambda_calc input option has been deprecated in
        ! favor of a global defintion of the coulomb logarithm in
        ! mod_lambda.f90
        ! 
       hc_lambda_calc = lambda_opt
        
      ! Check for cdperpc = -1.0 here.
      if (cdperpc.lt.0) then
        cdperpc = cdperp
      endif

      

       call pr_trace('READIN','END')

      RETURN

 1002 CALL PRC ('READIN: NOT ENOUGH DWELL TIMES GIVEN')
      IERR = 1
      RETURN
 1003 CALL PRC ('READIN: ITERATIONS > 1 INCOMPATIBLE WITH COLL OPT 3')
      IERR = 1
      RETURN
    END SUBROUTINE READIN




      subroutine read_sol23_params_unstruc

      use mod_params
      use mod_sol23_input
      implicit none
      integer ierr
!
!     READ_SOL23_PARAMS: This routine reads the SOL23 parameters
!     from the input file - it is invoked from IODIV. These have been
!     grouped together in the sol23 module in order to minimize
!     the impact of changes made to SOL 23 on teh rest of the code.
!
'+301 SOL23  ptipte ................................      '    1.0
      call rdr(sol23_par_ptipte,.true.,0.0,.true.,1.0E8,'par_ptipte       ',ierr)

'+302 SOL23  adaptnr                                      '     2
      call rdi(sol23_par_adaptnr,.true.,0  ,.true.,10  ,'par_adaptnr      ',ierr)

'+303 SOL23  debugflag                                    '     0
      call rdi(sol23_par_debugflag,.true.,0 ,.true.,10 ,'par_debugflag    ',ierr)

'+304 SOL23  debugnr                                      '     0
      call rdi(sol23_par_debugnr,.true.,0  ,.true.,10  ,'par_debugnr      ',ierr)

'+305 SOL23  refresh                                      '    100
      call rdi(sol23_par_refresh,.true.,0  ,.true.,10000,'par_refresh      ',ierr)

'+306 SOL23  artvisc2 ..............................      '    0.3
      call rdr(sol23_par_artvisc2,.true.,0.0,.true.,1.0E8,'par_artvisc2     ',ierr)

'+307 SOL23  artvisc4                                     '    0.01
      call rdr(sol23_par_artvisc4,.true.,0.0,.true.,1.0E8,'par_artvisc4     ',ierr)

'+308 SOL23  dtf                                          '    0.5
      call rdr(sol23_par_dtf,.true.,0.0,.true.,1.0E8,'par_dtf          ',ierr)

'+309 SOL23  dtg                                          '   1.0E24
      call rdr(sol23_par_dtg,.true.,0.0,.true.,1.0E36,'par_dtg          ',ierr)

'+310 SOL23  grid0                                        '    0.1
      call rdr(sol23_par_grid0,.true.,0.0,.true.,1.0E8,'par_grid0        ',ierr)

'+311 SOL23  gridexp ................................     '    1.25
      call rdr(sol23_par_gridexp,.true.,0.0,.true.,1.0E8,'par_gridexp      ',ierr)

'+312 SOL23  itermax                                      '    1000
      call rdi(sol23_par_itermax,.true.,0  ,.true.,10000,'par_itermax      ',ierr)

'+313 SOL23  ga1                                          '    10.0
      call rdr(sol23_par_ga1      ,.true.,0.0,.true.,1.0E8,'par_ga1          ',ierr)

'+314 SOL23  ga2                                          '    0.01
      call rdr(sol23_par_ga2    ,.true.,0.0,.true.,1.0E8,'par_ga2          ',ierr)

'+315 SOL23  ga3                                          '    0.2
      call rdr(sol23_par_ga3    ,.true.,0.0,.true.,1.0E8,'par_ga3          ',ierr)

'+316 SOL23  ga4   .................................      '    10.0
      call rdr(sol23_par_ga4     ,.true.,0.0,.true.,1.0E8,'par_ga4          ',ierr)

'+317 SOL23  updtqpit                                     '     0
      call rdi(sol23_par_updtqpit,.true.,0  ,.true.,10   ,'par_updtqpit     ',ierr)

'+318 SOL23  updtqp2                                      '    0.003
      call rdr(sol23_par_updtqp2 ,.true.,0.0,.true.,1.0E8,'par_updtqp2      ',ierr)

'+319 SOL23  updtqp3                                      '    0.03
      call rdr(sol23_par_updtqp3,.true.,0.0,.true.,1.0E8,'par_updtqp3      ',ierr)

'+320 SOL23  updtqpte                                     '    0.01
      call rdr(sol23_par_updtqpte,.true.,0.0,.true.,1.0E8,'par_updtqpte     ',ierr)

'+321 SOL23  garelax ................................     '    0.03
      call rdr(sol23_par_garelax,.true.,0.0,.true.,1.0E8,'par_garelax      ',ierr)

'+322 SOL23  gaiter                                       '     1
      call rdi(sol23_par_gaiter ,.true.,0  ,.true.,10000,'par_gaiter       ',ierr)

'+323 SOL23  limitte                                      '     1
      call rdi(sol23_par_limitte  ,.true.,0  ,.true.,10   ,'par_limitte      ',ierr)

'+324 SOL23  celldte                                      '    0.95
      call rdr(sol23_par_celldte,.true.,0.0,.true.,1.0E8,'par_celldte      ',ierr)

'+325 SOL23  updtbcte                                     '    30.0
      call rdr(sol23_par_updtbcte,.true.,0.0,.true.,1.0E8,'par_updtbcte     ',ierr)

'+326 SOL23  updtdel0 ..............................      '    1.0
      call rdr(sol23_par_updtdel0,.true.,0.0,.true.,1.0E8,'par_updtdel0     ',ierr)

'+327 SOL23  updtdel1                                     '    1.0
      call rdr(sol23_par_updtdel1,.true.,0.0,.true.,1.0E8,'par_updtdel1     ',ierr)

'+328 SOL23  updtdelm                                     '    0.1
      call rdr(sol23_par_updtdelm,.true.,0.0,.true.,1.0E8,'par_updtdelm     ',ierr)

'+329 SOL23  qbrelax                                      '    0.02
      call rdr(sol23_par_qbrelax,.true.,0.0,.true.,1.0E8,'par_qbrelax      ',ierr)

'+330 SOL23  gridg                                        '    0.03
      call rdr(sol23_par_gridg   ,.true.,0.0,.true.,1.0E8,'par_gridg        ',ierr)

'+331 SOL23  grid_dx0 ...............................     '    0.001
      call rdr(sol23_par_grid_dx0,.true.,0.0,.true.,1.0E8,'par_grid_dx0     ',ierr)

'+332 SOL23  gae                                          '    5.0
      call rdr(sol23_par_gae    ,.true.,0.0,.true.,1.0E8,'par_gae          ',ierr)

'+333 SOL23  gai                                          '    2.5
      call rdr(sol23_par_gai      ,.true.,0.0,.true.,1.0E8,'par_gai          ',ierr)

'+334 SOL23  tectrl                                       '     1
      call rdi(sol23_par_tectrl ,.true.,0  ,.true.,10   ,'par_tectrl       ',ierr)

'+335 SOL23  drflag                                       '     0
      call rdi(sol23_par_drflag ,.true.,0  ,.true.,10   ,'par_drflag       ',ierr)

'+336 SOL23  dsflag  ................................     '     0
      call rdi(sol23_par_dsflag  ,.true.,0  ,.true.,10   ,'par_dsflag       ',ierr)

'+337 SOL23  limrel1                                      '    0.03
      call rdr(sol23_par_limrel1 ,.true.,0.0,.true.,1.0E8,'par_limrel1      ',ierr)

'+338 SOL23  limrel2                                      '    0.1
      call rdr(sol23_par_limrel2 ,.true.,0.0,.true.,1.0E8,'par_limrel2      ',ierr)

'+339 SOL23  limrel3                                      '    0.03
      call rdr(sol23_par_limrel3 ,.true.,0.0,.true.,1.0E8,'par_limrel3      ',ierr)

'+340 SOL23  limrel4                                      '    0.03
      call rdr(sol23_par_limrel4 ,.true.,0.0,.true.,1.0E8,'par_limrel4      ',ierr)

'+341 SOL23  tulimit  ...............................     '    5.0
      call rdr(sol23_par_tulimit,.true.,0.0,.true.,1.0E8,'par_tulimit      ',ierr)

'+342 SOL23  g0relax                                      '    0.1
      call rdr(sol23_par_g0relax,.true.,0.0,.true.,1.0E8,'par_g0relax      ',ierr)

'+343 SOL23  p0relax                                      '    0.02
      call rdr(sol23_par_p0relax  ,.true.,0.0,.true.,1.0E8,'par_p0relax      ',ierr)

'+344 SOL23  dpdrflag                                     '     0
      call rdi(sol23_par_dpdrflag,.true.,0  ,.true.,10   ,'par_dpdrflag     ',ierr)

'+345 SOL23  dpdrtemin                                    '    1.0
      call rdr(sol23_par_dpdrtemin,.true.,0.0,.true.,1.0E8,'par_dpdrtemin    ',ierr)

'+346 SOL23  dpdrstep  ..............................     '    0.1
      call rdr(sol23_par_dpdrstep,.true.,0.0,.true.,1.0E8,'par_dpdrstep     ',ierr)

'+347 SOL23  nuflag                                       '     0
      call rdi(sol23_par_nuflag  ,.true.,0  ,.true.,10   ,'par_nuflag       ',ierr)

'+348 SOL23  nulimit                                      '   1.0E20
      call rdr(sol23_par_nulimit,.true.,0.0,.true.,1.0E24,'par_nulimit      ',ierr)

'+349 SOL23  vnmult                                       '    1.0
      call rdr(sol23_par_vnmult ,.true.,0.0,.true.,1.0E8,'par_vnmult       ',ierr)

'+350 SOL23  pinmom                                       '     0
      call rdi(sol23_par_pinmom  ,.true.,0  ,.true.,10   ,'par_pinmom       ',ierr)

'+351 SOL23  emolec                                       '    3.0
      call rdr(sol23_par_emolec ,.true.,0.0,.true.,1.0E8,'par_emolec       ',ierr)

'+352 SOL23  rec_heat  ...............................    '    13.6
      call rdr(sol23_par_rec_heat,.true.,0.0,.true.,1.0E8,'par_rec_heat     ',ierr)

'+353 SOL23  pinqimult                                    '    1.0
      call rdr(sol23_par_pinqimult,.true.,0.0,.true.,1.0E8,'par_pinqimult    ',ierr)

'+354 SOL23  pinqiflag                                    '     0
      call rdi(sol23_par_pinqiflag,.true.,0  ,.true.,10   ,'par_pinqiflag    ',ierr)

'+355 SOL23  prring0                                      '     7
      call rdi(sol23_par_prring0,.true.,0  ,.true.,10   ,'par_prring0      ',ierr)

'+356 SOL23  prring1                                      '     3
      call rdi(sol23_par_prring1 ,.true.,0  ,.true.,10   ,'par_prring1      ',ierr)

'+357 SOL23  qperp34 .................................    '     0
      call rdi(sol23_par_qperp34 ,.true.,0  ,.true.,10   ,'par_qperp34      ',ierr)

'+358 SOL23  qeiflag                                      '     1
      call rdi(sol23_par_qeiflag ,.true.,0  ,.true.,10   ,'par_qeiflag      ',ierr)

'+359 SOL23  chie                                         '     1
      call rdi(sol23_par_chie   ,.true.,0  ,.true.,10   ,'par_chie         ',ierr)

'+360 SOL23  joule                                        '     1
      call rdi(sol23_par_joule   ,.true.,0  ,.true.,10   ,'par_joule        ',ierr)

'+361 SOL23  fluxexp                                      '     1
      call rdi(sol23_par_fluxexp,.true.,0  ,.true.,10   ,'par_fluxexp      ',ierr)

'+362 SOL23  qrec    ..................................   '     1
      call rdi(sol23_par_qrec   ,.true.,0  ,.true.,10   ,'par_qrec         ',ierr)

'+363 SOL23  fzrad                                        '    1.0
      call rdr(sol23_par_fzrad,.true.,0.0,.true.,1.0E8,'par_fzrad        ',ierr)

'+364 SOL23  dvmrel                                       '     2
      call rdi(sol23_par_dvmrel,.true.,1,.true.,10,'par_dvmrel       ',ierr)

'+365 SOL23  intopt                                       '     1
      call rdi(sol23_intopt,.TRUE.,0 ,.true.,1,'intopt',IERR)

'+366 SOL23  adaptgrid                                    '     0
      call rdi(sol23_adaptgrid,.TRUE.,0 ,.true.,5,'adaptgrid',IERR)

'+367 SOL23  seed                                         '     1
      call rdi(sol23_seed,.TRUE.,1 ,.true.,2,'seed',IERR)

'+368 SOL23  izlen   ..................................   '    0.5
      call rdr(sol23_izlen,.true.,0.0,.false.,0.5,'izlen',ierr)

'+369 SOL23  izlam                                        '    0.03
      call rdr(sol23_izlam,.true.,0.0,.false.,0.5,'izlam',ierr)

'+370 SOL23  izoffset                                     '    0.03
      call rdr(sol23_izoffset,.true.,0.0,.false.,0.5,'izoffset',ierr)

'+371 SOL23  momlen                                       '    0.1
      call rdr(sol23_momlen,.true.,0.0,.false.,0.0,'momlen',ierr)

'+372 SOL23  relax                                        '    0.01
      call rdr(sol23_relax,.true.,0.0,.true.,1.0,'relax',ierr)

'+373 SOL23  maxtol                                       '    3.0E-3
      call rdr(sol23_maxtol,.true.,0.0,.true.,1.0,'maxtol',ierr)

'+374 SOL23  rmstol  .................................    '    3.0E-4
      call rdr(sol23_rmstol,.true.,0.0,.true.,1.0,'rmstol',ierr)

'+375 SOL23  Qperp                                        '    1
      call rdi(sol23_perp,.true.,0 ,.true.,10 ,'perp',ierr)

      return





        

        
      end subroutine read_sol23_params_unstruc




      subroutine read_sol22_params_unstruc

        implicit none


    !use mod_params
    use debug_options
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    implicit none

    !     This subroutine reads the input parameters for the case
    !     from the standard input or redirected from a file.

    !     include 'params'

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    integer ierr

    !      integer ierr

    !     Model parameters

    !real*8 spts(mxspts)

    call pr_trace('MOD_SOL22_INPUT','START OF READSOL')
    
'+201    Force Te=Ti through SOL 22  0=off  1=on         '     0 1
    call rdi(forcet,.TRUE.,0,.true.,1,         'force te=ti'    ,ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER FORCET')

'+202    Imposed mach number at the target               '    1.0
    CALL RDQ(initm0,.TRUE.,0.0d0,.FALSE.,0.0d0,'target mach num',IERR)

'+203    Delta mach number for initial iterative solution'    0.1
    CALL RDQ(deltam0,.TRUE.,0.0d0,.FALSE.,0.0d0,'delta mach num',IERR)

    !     ------------------------------------------------------------------

    !     Ionization Source

'+204    Maximum resolution in calculation of m0         ' 0.00001
    CALL RDQ(m0res,.TRUE.,0.0d0,.FALSE.,0.0d0,'Resolution in m0',IERR)

'+205    Ionization Source Lengths  0=Absolute 1=Relative'     1
    CALL RDI(lensind,.TRUE.,0,.TRUE.,1,     'ion source abs/rel',IERR)

'+206    Start of Ionization Source (for supported opts) '    0.0
    CALL RDQ(lensst,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src start',  IERR)

'+207    End or Length of Ionization Source              '    0.4
    CALL RDQ(lensfi,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src finish', IERR)

'+208    Decay length of ionization source               '    0.05
    CALL RDQ(lams,.TRUE. ,0.0d0,.FALSE.,0.0d0,'ion decay len   ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER IONIZATION SOURCE')
    !     Radiation source

'+209    Length of radiation source                      '    5.0
    CALL RDQ(lenri,.TRUE.,0.0d0,.FALSE.,0.0d0,'rad source len  ',IERR)

'+210    Decay length of radiation  source               '    0.5
    CALL RDQ(lamri,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad decay len  ',IERR)

'+211    Source strength fraction (frr)                  '    2.0
    CALL RDQ(frri ,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad power mult ',IERR)

'+212    Garching Model: Alpha = ni/ne ratio             '    1.0
    CALL RDQ(alfimp,.TRUE.,0.0d0,.FALSE.,0.0d0,'nimp/ne ratio  ',IERR)

'+213    Garching Model: Tbase = Tratio denominator      '   15.0
    CALL RDQ(talimp,.TRUE. ,0.0d0,.FALSE.,0.0d0,'base Temp     ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER RADIATION SOURCE')
    !     Miscellaneous

'+214    Garching Model: Exponent 1                      '    1.5
    CALL RDQ(ex1imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 1  ',IERR)

'+215    Garching Model: Exponent 2                      '   -3.0
    CALL RDQ(ex2imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 2  ',IERR)

'+216    Gamma correction factor in gammai               '    0.0
    CALL RDQ(gamcor,.false.,0.0d0,.FALSE.,0.0d0,'i power corr. ',IERR)

'+217    Gamma correction factor in gammae               '    0.0
    CALL RDQ(gamecor,.false.,0.0d0,.FALSE.,0.0d0,'e power corr.',IERR)

'+218    CX power coefficeint CEICF                      '    1.0
    CALL RDQ(ceicf,.TRUE. ,0.0d0,.FALSE.,0.0d0, 'CX power frac ',IERR)

'+219    Recycling  source fraction                      '    1.0
    CALL RDQ(recfrac,.TRUE.,0.0d0,.TRUE.,1.0d0,'Recycle frac ',IERR)

'+220    Pei Power Transfer Correction Factor            '    1.0
    call rdq(peicf,.true.  ,0.0d0,.false.,0.0d0,'Pei Correction',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MISC')
    !     Power distribution

'+221    Velocity Error Switch       0=Cs   1=const      '     1
    call rdi(velsw,.true.  ,0,.true.,3,       'Vel Error Switch',ierr)

'+222    Distributed Power Start position * SMAX         '    0.10
    call rdq(spowbeg,.true.,0.0d0,.true.,0.5d0,'Power Dist Beg',ierr)

    !     Gperp Distribution function

'+223    Distributed Power End   position * SMAX         '    0.50
    call rdq(spowlen,.true.,0.0d0,.true.,0.5d0,'Power Dist Len',ierr)

'+224    Distributed GPERP particle Fraction- non-uniform'    0.8
    call rdq(gperpfrac,.true.,0.0d0,.true.,1.0d0,'Part Dist Frac',ierr)

'+225    Distributed GPERP Start position * SMAX         '    0.0
    call rdq(gperpbegf,.true.,0.0d0,.true.,0.5d0,'Part Dist Beg',ierr)

    !     Extra Gperp source/sink - start and end positions.

'+226    Distributed GPERP End   position * SMAX         '    0.1
    call rdq(gperpendf,.true.,0.0d0,.true.,0.5d0,'Part Dist Len',ierr)

'+227    Gextra Source strength - Target flux multiplier '    0.1
    call rdq(gextra_mult,.true.,0.0d0,.false.,0.0d0,'Gextra flux mult',ierr)

'+228    Gextra Source Start/Stop * SMAX                 ' 0.2   0.35
    call rdq2(gextra_src_start,gextra_src_stop,.true.,0.0d0,.true.,1.0d0, 'Gextra SRC start/stop',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER GPERP')

    !     Field line length fraction for distributing the Private plasma
    !     electron and ion power loads.

'+229    Gextra Sink   Start/Stop * SMAX                 ' 0.65  0.8
    call rdq2(gextra_sink_start,gextra_sink_stop,.true.,0.0d0,.true.,1.0d0,'Gextra SINK start/stop',ierr)

'+230    PP target power loss redistribution range *SMAX '    0.1
    call rdq(pp_pow_dist,.true.,0.0d0,.false.,0.0d0,'PP Pow Dist',ierr)

'+231    Knot Start Index for E2D Option 9               '     8
    call rdi(ike2d,.true.,1,.false.,0,'IK start Index for E2D-9',ierr)

'+232    Plasma Fill option for missing knots - E2D Opt 9'     1
    call rdi(fillopt,.true.,0,.true.,3,'Gap fill option - E2D-9',ierr)

'+233    Qe Term - Temperature Cutoff (eV)               '    0.0
    call rdq(tcutqe,.true.,0.0d0,.false.,0.0d0,'Cut T-PINQE',ierr)

'+234    PINQID - Atomic Ionization    - T cutoff (eV)   '    0.0
    call rdq(tcutatiz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDATIZ',ierr)

'+235    PINQID - Molecular Ionization - T cutoff (eV)   '    0.0
    call rdq(tcutmliz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDMLIZ',ierr)

'+236    PINQID - Recombination        - T cutoff (eV)   '    0.0
    call rdq(tcutrec,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDREC',ierr)

'+237    Qi Term/PINQID-Charge Exchange- T cutoff (eV)   '    0.0
    call rdq(tcutcx,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDCX',ierr)

'+238    PINQID - CX option 1 - Reference T - (eV)       '    1.0
    call rdq(trefcx,.true.,0.0d0,.false.,0.0d0,'REF T- QIDCX 1',ierr)

'+239    Minimum Temperature allowed in Solver (spec<0)  '    0.1
    call rdq(tmin,.false.,0.0d0,.false.,0.0d0,'Min. Allowed T',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER PINQID')

'+240    Minimum T allowed as fraction of Tmax reached   '    0.5
    call rdq(dropfrac,.true.,0.0d0,.true.,1.0d0,'Allowed T-drop',ierr)

'+241    Momentum loss term multiplier   (Usually 1.0)   '    1.0
    call rdq(smom_mult,.false.,0.0d0,.false.,0.0d0,'Mom.Loss Multiplier',ierr)

'+242    Friction factor for Momentum loss formula       '    0.2
    call rdq(ffric,.false.,0.0d0,.false.,0.0d0,'Mom.loss frac  ',ierr)

'+243    Length of the Momentum loss region * Smax       '    0.1
    CALL RDQ(lenmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom source len',IERR)

'+244    Decay length of the Momentum loss region * Smax '    0.02
    CALL RDQ(lammom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom decay len ',IERR)

'+245    Ratio of CX to IZ events (fixed)                '    1.0
    CALL RDQ(rcxmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'cx/iz ratio   ',IERR)

'+246    Te cutoff for increased CX multiplier (eV)      '    5.0
    call rdq(tcxmom,.TRUE.,1.0001d0,.FALSE.,0.0d0,'T for CXmult',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MOMENTUM')

'+247    Te lower limit cutoff for CX multiplier (eV)    '    1.0
    call rdq(tcxcut,.true.,0.0d0,.false.,0.0d0,'Cut T -  CXmult',ierr)

'+248    PINQE multiplier                                '    1.0
    call rdq(qesrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE mult',ierr)

'+249    PRAD option 3 multiplier (x PINQE)              '    1.0
    call rdq(radsrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE based PRAD mult',ierr)

'+250    Initial number of stages for Runge Kutta steps  '    100
    CALL RDI(ndiv,.TRUE., 1,.FALSE., 0,'NUMBER OF STEPS    ',IERR)

    call pr_trace('MOD_SOL22_INPUT','BEFORE SWITCHES')

'+251    Switch: Ionization Opt : 0.0-exp 1.0+ - others  '    2.0
    call rdr(switch(swion),.true. ,0.0,.false.,0.0, 'alt ion',ierr)

'+252    Switch: Initial IonOpt : 0.0-exp 3.0+ - others  '    0.0
    call rdr(switch(swioni),.true. ,0.0,.false.,0.0,'init ion',ierr)

'+253    Switch: PPlasma IonOpt : 0.0-exp 3.0+ - others  '   -6.0 -5.0 -3.0 0.0
    call rdr(switch(swionp),.false.,0.0,.false.,0.0,'pp ion',ierr)


    ! move to after input file read - at the beginning of SOL22 processing?
    if (switch(swionp).eq.-1.0) then
       if (switch(swion).eq.1.0.or.switch(swion).eq.2.0.or.switch(swion).eq.8.0) then
          switch(swionp) = switch(swioni)
       else
          switch(swionp) = switch(swion)
       endif

    endif

'+254    Switch: 5/2 nv * kT    : 0.0-off 1.0-on         '    1.0
    call rdr(switch(swcond),.true.,0.0,.false.,0.0, 'cond sw',ierr)

'+255    Switch: 1/2 m v^3 * n  : 0.0-off 1.0-on         '    1.0
    call rdr(switch(swconv),.true.,0.0,.false.,0.0, 'conv sw',ierr)

'+256    Switch: Prad           : 0.0-off 1.0-on         '    3.0
    call rdr(switch(swprad),.true.,0.0,.false.,0.0, 'prad sw',ierr)

'+257    Switch: Phelpi         : 0.0-off 1.0-on         '    2.0
    call rdr(switch(swphelp),.true.,0.0,.false.,0.0,'phelp sw',ierr)

'+258    Switch: Pei            : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swpei),.true.,0.0,.false.,0.0,  'pei sw ',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER BASIC SWITCHES')

'+259    Switch: Pcx            : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swpcx),.true.,0.0,.false.,0.0,  'pcx sw ',ierr)

'+260    SUB-switch: Pcx Opt 4  : PINQID- Atomic Ioniz.  '    1.0
    call rdr(switch(swqidatiz),.true.,0.0,.false.,0.0,'atiz sw',ierr)

'+261    SUB-switch: Pcx Opt 4  : PINQID- Molecular Ioniz'    1.0
    call rdr(switch(swqidmliz),.true.,0.0,.false.,0.0,'mliz sw',ierr)

'+262    SUB-switch: Pcx Opt 4  : PINQID- Recombination  '    1.0
    call rdr(switch(swqidrec),.true.,0.0,.false.,0.0,'rec sw',ierr)

'+263    SUB-switch: Pcx Opt 4  : PINQID- Charge Exchange'    2.0
    call rdr(switch(swqidcx),.true.,0.0,.false.,0.0,'cx sw',ierr)

'+264    Switch: PP ElecLoss    : 0.0-off 1.0-XPT 2.0-DIS'    1.0
    call rdr(switch(swppelec),.true.,0.0,.false.,0.0,'pp elec sw',ierr)


    !     jdemod -
    !     switch(swppress) is read in using tag 283 of unstructured input
    !     default value is 0 or OFF

    !      call rdr(switch(swppress),.true.,0.0,.false.,0.0,'pp power sw',ierr)


'+265    Switch: PP IonLoss     : 0.0-off 1.0-XPT 2.0-DIS'    1.0
    call rdr(switch(swppion),.true.,0.0,.false.,0.0,'pp ion sw',ierr)

'+266    Switch: Visc 1 - N calc: 0.0-off 1.0-on         '    0.0
    call rdr(switch(swvisc1),.true.,0.0,.true.,0.0,'visc1 sw',ierr)

'+267    Switch: Momentum loss  : 0.0-off 1.0-on         '    9.0
    call rdr(switch(swnmom),.true.,0.0,.false.,0.0, 'N mom sw',ierr)

    call pr_trace('MOD_SOL22_INPUT','HALFWAY THROUGH SWITCHES')


'+268    Switch: Iterative Mach : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swmach),.true.,0.0,.false.,0.0, 'mach sw',ierr)

'+269    Switch: Edge 2D Data   : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swe2d),.false.,0.0,.false.,0.0, 'e2d sw',ierr)

'+270    Switch: Power Distrib. : 0.0-con 1.0-lin 2.0-xpt'    6.0
    call rdr(switch(swpow),.true.,0.0,.false.,0.0, 'power sw',ierr)

'+271    Switch: PPlasma PowDist: 0.0-con 1.0-lin 2.0-xpt'   10.0
    call rdr(switch(swpowp),.false.,0.0,.false.,0.0, 'pp pow',ierr)

    ! move to after input file read
    if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)

'+272    Switch: Gamma Perp     : 0.0-off 1.0-on         '    2.0
    call rdr(switch(swgperp),.true.,0.0,.false.,0.0,'GamPerp',ierr)

'+273    Switch: PP Gamma Perp  : 0.0-off 1.0-on         '    6.0
    call rdr(switch(swgperpp),.false.,0.0,.false.,0.0,'GamPerpP',ierr)

    !     Extra Gperp source/sink term

    ! move to after input file read
    if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)

'+274    Switch: GPero Src/Sink : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swextra),.true.,0.0,.false.,0.0,'GP Src/Sink',ierr)

'+275    Switch: Major Radius   : 0.0-off 1.0-nor 2.0-inv'    0.0
    call rdr(switch(swmajr),.true.,0.0,.false.,0.0,'MajorRad',ierr)

'+276    Switch: Core Gamma Src : 0.0-off 1.0-all 2.0-xpt'    0.0
    call rdr(switch(swcore),.true.,0.0,.false.,0.0,'Core Src',ierr)

'+277    Switch: Recomb. Src    : 0.0-off 1.0-PIN        '    1.0
    call rdr(switch(swrecom),.true.,0.0,.false.,0.0,'Recomb.',ierr)

'+278    Switch: Smooth Centres : 0.0-off 1.0-on         '    0.0
    call rdr(switch(swsmooth),.true.,0.0,.false.,0.0,'Smooth',ierr)

'+279    Switch: Detached Option: 0.0-off 1.0-out 2.0-in '    0.0
    call rdr(switch(swdetach),.true.,0.0,.false.,0.0,'Detach',ierr)

'+280    Switch: Error corrected: 0.0-off 1.0-cond       '   10.0
    call rdr(switch(swerror),.true.,0.0,.false.,0.0,'ERROR',ierr)

    call pr_trace('MOD_SOL22_INPUT','DONE SWITCHES')


    !CALL RDRARN(deflist,ndef,mxspts,0.0,real(maxnrs),.FALSE.,0.0,MACHHI,2,'DEFAULT SOLVER DATA',IERR)
    !
    ! jdemod - The max value should the maximum number of rings but to avoid dependency on params
    !          module I am removing this constraint. 

'+281 ' 'Automatic DEFAULT Solver condition switches     '
'    DEFAULT applied automatically to these rings'         0
    CALL RDRARN(deflist,ndef,mxspts,0.0,SOL22_HI,.FALSE.,0.0,SOL22_MACHHI,2,'DEFAULT SOLVER DATA',IERR)

    !     Default plots to off - this can be changed in the calcsol_interface
    !     routine in solascv1.f


    ! move to initialization

    !      CALL RDI(graph,.TRUE., 0,.TRUE., 1,'GRAPH OPTION       ',IERR)
    !      CALL RDI(graphaux,.TRUE.,0,.TRUE.,1,'AUX GRAPH OPTION  ',IERR)
    !      CALL RDI(graphvel,.TRUE.,0,.TRUE.,1,'VEL GRAPH OPTION  ',IERR)
    !      call rdr(graphran,.true.,0.0,.false.,0.0,'CXMAX VALUE  ',ierr)

    graph = 0
    graphaux = 0
    graphvel = 0
    graphran = 0.0


    call pr_trace('MOD_SOL22_INPUT','END OF READSOL')


    return




      end subroutine read_sol22_params_unstruc



        
