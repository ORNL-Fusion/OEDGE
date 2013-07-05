c
c     This module contains code related to the unstructured input options
c
c
c ======================================================================
c
c subroutine: ValidateUnstructuredInput
c
c
      SUBROUTINE ValidateUnstructuredInput
      IMPLICIT none
      include 'params'
      include 'comtor'

      INTEGER    MAXTAG
      PARAMETER (MAXTAG=1000) 
      COMMON /INPCHK/ ntaglist,taglist 
      INTEGER     ntaglist
      CHARACTER*3 taglist(MAXTAG)

      INTEGER    MAXCHK
      PARAMETER (MAXCHK=15)
      INTEGER      i1,i2,nchklist,mchklist(MAXCHK)
      CHARACTER*72 chklist(2*MAXCHK)
      CHARACTER*69 sp

c...  List of required unstructured input tags:
      DATA chklist
     .  /'026',' PIN selection  0-NIMBUS 1-EIRENE97 2-EIRENE99     (1)',
     .   '010',' Geometry data  0-standard 1-from DIVIMP              ',
     .   '021',' Input file     0-standard 1-from DIVIMP              ',
     .   '020','   Run time (CPU seconds)                             ',
     .   '022','   Material: target  1-Mo 2-C 3-W 4-Be                ',
     .   '024','             wall                                     ',
     .   '011','   Grid type     0-structured 1-generalized           ',
     .   '018','   Wall data     0-standard   1-seamless              ',
     .   '019','   Debug option  0-off                                ',
     .   'E11',' n-n collisions  0-off 1-standard mesh             (0)',

     .   'E12',' Lyman alpha opacity  0-off 1-rec 2-rec&ion        (0)',
     .   '058',' 1.1 Pressure gauge specification:                    ',
     .   '076',' 1.0 Surface properties:                              ',
     .   '077',' 1.0 Additional surfaces:                             ',
     .   '078','   Target data shift  JET-i/o CMOD,DIIID-o/i (m) (0.0)'/

c
c jdemod
c
c     For unstructured input items whose default value is intended to match a 
c     regular DIVIMP input - assign the DIVIMP values if the unstructured input
c     still contains its default value. At present this applies to some of the HC
c     code input data.
c
      call global_hc_assign_inputs
c
c jdemod
c
c

c
c jdemod 
c
c     For ion injection options 9 and 10 make sure that values have been 
c     specified for the endpoints/corners of the injection line/region      
c
      if (cneuta.eq.1.and.(ciopte.eq.9.or.ciopte.eq.10).and.
     >     (cxsca.eq.0.0.and.cysca.eq.0.0.and.
     >      cxscb.eq.0.0.and.cyscb.eq.0.0)) then 
c
c         Invalid inputs specified
c
          write(0,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte
          write(0,*) '       INJECTION REGION NOT SPECIFIED' 
          write(0,*) '       SEE INPUT ITEMS *I29 and *I30'
          write(6,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte
          write(6,*) '       INJECTION REGION NOT SPECIFIED' 
          write(6,*) '       SEE INPUT ITEMS *I29 and *I30'
          write(7,*) 'ERROR: Invalid Input - INJECTION OPTION =', ciopte
          write(7,*) '       INJECTION REGION NOT SPECIFIED' 
          write(7,*) '       SEE INPUT ITEMS *I29 and *I30'
c
       endif
c
      if (cneuta.eq.0.and.(cneutb.eq.6.or.cneutb.eq.7).and.
     >     (cxsca.eq.0.0.and.cysca.eq.0.0.and.
     >      cxscb.eq.0.0.and.cyscb.eq.0.0)) then 
c
c         Invalid inputs specified
c
          write(0,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb
          write(0,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
          write(0,*) '       SEE INPUT ITEMS *I29 and *I30'
c
          write(6,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb
          write(6,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
          write(6,*) '       SEE INPUT ITEMS *I29 and *I30'
c
          write(7,*) 'ERROR: Invalid Input- NEUT LAUNCH OPTION =',cneutb
          write(7,*) '       NEUTRAL LAUNCH REGION NOT SPECIFIED' 
          write(7,*) '       SEE INPUT ITEMS *I29 and *I30'
c
       endif
c
c jdemod
c
      RETURN
90    FORMAT(A,A3,A)
97    DO i1 = 1, MAXCHK
        IF (mchklist(i1).NE.1) THEN
          WRITE(6,90) '  ',chklist(2*i1-1),chklist(2*i1)
          WRITE(0,90) '  ',chklist(2*i1-1),chklist(2*i1)
        ENDIF
      ENDDO
      GOTO 99
98    WRITE(6,90) 'SAMPLE: ',chklist(2*i2-1),chklist(2*i2)
99    STOP
      END
c
c
c
c
c ======================================================================
c
c subroutine: InitializeUnstructuredInput
c
      subroutine InitializeUnstructuredInput
      use subgrid_options
      use ribbon_grid_options
      use sol22_input
      use allocatable_input_data
      implicit none
      
      INCLUDE 'params'

      INCLUDE 'slcom'
      INCLUDE 'cgeom'
      include 'cadas'
      INCLUDE 'comtor'
c      INCLUDE 'pindata'
      INCLUDE 'cedge2d'  

      INCLUDE 'solparams'
      INCLUDE 'solswitch'
      INCLUDE 'solcommon'

      include 'reiser_com'  
      include 'line_profile' 
c
      include 'driftvel'
      include 'fperiph_com'
      include 'dperpz'
c      include 'slcom_sol28'
c

c
c
c     Intializing Unstructured input data to default values. 
c
c -----------------------------------------------------------------------
c
c     David's options  
c
c -----------------------------------------------------------------------
c
c     TAG 282: SOL22
c
c     Initialization of Array input for tag 282 specifying ffric 
c     values on a ring by ring basis for both targets. 
c
      n_extffric = 0
      call qzero(extffric,maxnrs*3) 
c     
c     TAG 283: SOL22 - private plasma pressure loss option
c
c     Set the default for this value to OFF = 0
c
      switch(swppress) = 0.0

c     
c     TAG 284: SOL22 - debug SOL22 
c
c     Set the default for this value to OFF = 0
c
      debug_sol22 = 0
c     
c     TAG 284: SOL22 - debug SOL22 
c
c     Set the default for this value to OFF = 0
c
      debug_sol22_ir = 1
c     
c     TAG 284: SOL22 - debug SOL22 
c
c     Set the default for this value to OFF = 0
c
      debug_sol22_ikopt = 1
c
c -----------------------------------------------------------------------
c
c     TAG A05
c
c     ne_opt - selects how electron density is defined 
c     0 -> ne = nb    1 -> ne = nb + sigma nz (from fluid code)
c
      ne_opt   = 0
c
c -----------------------------------------------------------------------
c
c     TAG A06
c
c     Option to write a JET TRAN file for POST-PROCESSOR use
c     Written to Unit 41 - 0 = off - 1 = 0n .. default ON
c
      write_tran = 1
c
c -----------------------------------------------------------------------
c
c     TAG C21:
c      
c     Initialization of PINQE multiplier for the Dperp extractor - this 
c     value is usually 1.0 - however, by allowing for this multiplier 
c     it is possible to adjust for impurity radiated power when this 
c     information is not available. 
c
      dp_pinqe_mult = 1.0
c
c -----------------------------------------------------------------------
c 
c     TAG C22
c
c     line_profile_opt is set to 0 (off) by default. There are a number
c     of input values required when this option is active and these 
c     input values MUST immediately follow the unstructured input option.
c     
c     The first line contains a set of ADAS selectors in the standard 
c     format. 
c
c     'Text'   'ADASID'  ADASYR 'ADAS EXTENSION'  ISELE ISELR ISELX ISELD
c     - only ISELE is used at the present time. 
c
c     The second input line defines the LOS for the calculation and 
c     the instrument and bin characteristics.
c
c     'Text'  ROBS  ZOBS  THETA  DTHETA INSTRUMENT_WIDTH BIN_WIDTH
c
c     No default values are assigned to any sub-options when this 
c     option is off - all options MUST be specified when it is on. 
c
      line_profile_opt=0 
c
c -----------------------------------------------------------------------
c
c     TAG D37 and D38 
c     ADAS IONIZATION AND RECOMBINATION RATE MODIFIERS
c
c     These values can be used to modify the rates read in from ADAS
c     The default values should always be set to 1.0 and these
c     should be used with care if used at all.  
c
c     D37 - ADAS Ionization rate multiplier
c
      adas_iz_rate_mult = 1.0 
c
c     D38 - ADAS Recombination rate multiplier
c
      adas_rec_rate_mult = 1.0 
c
c
c -----------------------------------------------------------------------
c
c     TAG D39 : Alternate Sputter data specifier - used to select one of 
c               several custom sputter datasets - usually based
c               on different impact angles
c 
c     Set to normal incidence data as default
c
      extra_sputter_angle = 0.0
c
c
c -----------------------------------------------------------------------
c
c     TAG F11
c
c     This is set to 1 to indicate that a UEDGE/fluid code background has
c     been loaded and that the fluid code ionization data has been loaded
c     into PIN arrays. 
c
      uedge_bg = 0
c
c -----------------------------------------------------------------------
c
c     TAG F12:
c
c     fc_target_calc_option - this option affects the calculation  of the 
c     target conditions that are extracted from the background plasma 
c     of a fluid code solution. This affects the values that are assigned
c     to the e2dtarg array at the point when the plasma solution is read in. 
c 
c     Each of the sub-options can also be set explicitly to different 
c     values. 
c
c     The fc_target_calc_option's supported are: 
c
c     Option 0: EDGE2D standard
c               fc_v_calc_opt  = 0
c               fc_te_calc_opt = 1
c               fc_ti_calc_opt = 2
c               fc_ne_calc_opt = 2
c     Option 1: UEDGE standard  
c               fc_v_calc_opt  = 0
c               fc_te_calc_opt = 0
c               fc_ti_calc_opt = 0
c               fc_ne_calc_opt = 2
c     Option 2: Base B2/Eirene standard
c               fc_v_calc_opt  = 1
c               fc_te_calc_opt = 2
c               fc_ti_calc_opt = 2
c               fc_ne_calc_opt = 2
c     Option 3: Alternate B2/Eirene
c               fc_v_calc_opt  = 1
c               fc_te_calc_opt = 1
c               fc_ti_calc_opt = 1
c               fc_ne_calc_opt = 2
c
c     F13:
c
c     fc_ne_calc_opt: 0 : Ti is value in guard cell
c                   : 1 : Ti is arithmetic average of guard and first cell
c                   : 2 : Ti is value in real cell
c
c     F14:
c
c     fc_te_calc_opt: 0 : Te is value in guard cell
c                   : 1 : Te is arithmetic average of guard and first cell
c                   : 2 : Te is value in real cell
c
c     F15: 
c
c     fc_ti_calc_opt: 0 : Ti is value in guard cell
c                   : 1 : Ti is arithmetic average of guard and first cell
c                   : 2 : Ti is value in real cell
c
c     F16:
c
c     fc_v_calc_opt  : 0 : Velocity is set to value at boundary from fc      
c                    : 1 : Velocity is set to sound speed
c
c
c     Set defaults to match EDGE2D target condition interpretation option
c
      fc_target_calc_option = 0
      fc_v_calc_opt  = 0
      fc_te_calc_opt = 1
      fc_ti_calc_opt = 2
      fc_ne_calc_opt = 2
c
c -----------------------------------------------------------------------
c
c     F17:
c
c     B2 and B2.5 write a fort.31 file with an identical format - however
c     the parallel velocity arrays are different in the two cases - for 
c     B2 the velocity is defined on the cell edges while for B2.5 the
c     velocity array in this file is defined at the cell centers. This means
c     that the data need to be interpreted differently in the two cases. 
c     This option may also apply to the fluxes that are written to this
c     file.
c
c     option 0 = cell boundary value in file
c     option 1 = cell center value in file
c
      fc_v_interp_opt = 0
c
c -----------------------------------------------------------------------
c
c     TAG G23:
c
c     This option is used to tag a SONNET style grid as being an 
c     FRC (Field Reversed Configuration) custom grid. At the moment 
c     only one option is supported but further subtypes could be 
c     added for Sonnet grids requiring special processing. 
c
c     0 = stanard Sonnet grid = Default 
c     1 = FRC version 1 - type of Sonnet grid 
c         - used to set various FRC related options
c     2 = sonnet grid without boundary cells - boundary cells are added
c         - useful for carre grids
c
      sonnet_grid_sub_type = 0
c
c -----------------------------------------------------------------------
c
c     TAG G34:
c
c     This input value is used to define the machine type for 
c     the purpose of placing data in the TRAN file for JET 
c     post-processor use. The default value is for JET. Unfortunately,
c     the grid option itself is not sufficiently selective since it defines
c     formats of grids that may or may not be associated with specific
c     machines.
c
c       0 = jet
c       1 = diiid
c       2 = alcator (cmod)
c       3 = aug (asdex upgrade) 
c       4 = iter
c       5 = ignitor
c       6 = fire
c
      tmachine_opt = 0
c
c -----------------------------------------------------------------------
c
c     TAG G35:
c
c     This input value is used to define a shotid
c     for the case for cataloguing in the JET catalog
c     system. This must include the shot number. If the 
c     shot number is present in the case name then this 
c     value is not required - the case name may be used as
c     the shotid in that case. This input takes a default
c     value of ' '. Meaning that the shot number is expected to 
c     be a part of the case name. 
c
      divshotid = ' '
c
c
c -----------------------------------------------------------------------
c
c     TAG G36:
c
c     s_reflect_opt
c
c     This input value is used to activate parallel ion 
c     reflection. This will insert mirrors for parallel ion transport
c     at ONE parallel break in the grid along each field line. This is
c     used to model ion transport on double null half grids. The default
c     value is set to off or zero.
c
      s_reflect_opt = 0
c
c -----------------------------------------------------------------------
c
c     G37: Used by Steve for some sort of grid option - no default 
c          values specified
c
c -----------------------------------------------------------------------
c
c     The following tags are related to the subgrid option for recording 
c     more detailed data on a finer grid
c
c     G38: Base subgrid option ON/OFF
c     G39: R,Z dimensions of the grid region
c     G40: RMIN,RMAX of the subgrid region
c     G41: ZMIN,ZMAX of the subgrid region
c
c     Note: Although default values are assigned - these values
c           should be explicitly set when this option is used. 
c
c     G38: Base subgrid option - OFF
c
      subgrid_opt = 0
c
c     G39: Dimensions of subgrid
c
      sg_rdim=100
      sg_zdim=100
c
c     G40: RMIN and RMAX values of the gridded region
c
      sg_rmin = 1.0
      sg_rmax = 1.5
c
c     G41: ZMIN and ZMAX values of the gridded region
c
      sg_zmin = 0.0
      sg_zmax = 1.4
c
c -----------------------------------------------------------------------
c
c     Options related to ribbon grids
c     G42 - grid generation option - <i4>
c     G43 - intersection point averaging option - opt_block_av - <r4>
c     G44 - maximum R separation in grid generator - max_r_sep - <r4>
c     G45 - maximum S/Z separation in grid generator - max_s_sep - <r4>
c     G46 - min number of cells on ring - min_cells - <i4>
c     G47 - castem output identifier - <string>
c     G48 - min and max S for selecting intersection subset  2 x <r4>
c     G49 - min and max R for intersection subset grid generation 2 x <r4>
c     G50 - min and max S for intersection subset grid generation 2 x <r4>
c     G51 - length cutoff factor for ring generation <r4> default = 0.0
c     G52 - Cell spacing option
c     G53 - cell spacing factor
c     G54 - Input file option - RAY or CASTEM 
c
c------------------------------------------------------------------------
c
c     G42 - grid option
c           0 = unstructured
c           1 = structured
c           default = unstructured
c           
      rg_grid_opt = 0
c
c     G43 - block averaging option (removes blobs of intersection data)
c     
      rg_block_av = 0
c
c     G44 - maximum R separation between rows
c      
      rg_max_r_sep = 0.002
c
c     G45 - maximum S/Z separation between cells along row
c
      rg_max_s_sep = 0.5
c
c     G46 - minimum number of cells in a row
c
      rg_min_cells = 5
c
c     G47 - Castem data set to read in
c
      rg_castem_data = '100610'
c
c     G48 - min and max S for selecting intersection subset  2 x <r4>
c           if min=max then window option is not selected
c
      rg_int_win_mins = 0.0
      rg_int_win_maxs = 0.0
c
c     G49 - min and max R for intersection subset grid generation 2 x <r4>
c           These are only used for subset grid generation
c
      rg_minr=0.0
      rg_maxr=0.0
c
c     G50 - min and max S for intersection subset grid generation 2 x <r4>
c           These are only used for subset grid generation
c
      rg_mins=0.0
      rg_maxs=0.0
c
c     G51 - length cutoff to eliminate short rings far from the separatrix from 
c           ring generation ... some testing of this will be necessary
c           to obtain an optimal grid. default value is 0.0 which 
c           effectively turns this feature off.
c
      lcutoff = 0.0
c
c     G52 - Cell spacing option - default value is exponential with the exponent
c           specified by G53. A value of 1.0 for the cell spacing factor gives
c           a linear spacing. This option works better with structured grids at
c           the moment. 
c
      cell_spacing_option = 0

c     G53 - Cell spacing factor for determining the distribution of cells
c           between fixed points on rings. 
c           default = 1.0 which gives a linear spacing
c      
      cell_spacing_factor = 1.0
c
c     G54 - Intersection data input file format option 
c           0 = CASTEM
c           1 = RAY
c
      ribbon_input_format_opt = 1
c
c
c -----------------------------------------------------------------------
c
c     HC related variable initializations
c
c     TAG H15 to H64 and H90,H91:
c
c     Insert code to initialize the HC variables 
c
      call global_hc_init
c
c -----------------------------------------------------------------------
c
c     TAG I24 
c
c     init_pos_opt - this option affects the initial position of neutrals 
c     generated by ions which stike the target and subsequently recycle
c     as well as the initial position of ions that are formed from neutrals. 
c     Option 1 for neutrals uses the cross component at the time of target 
c     impact to estimate a more precise R,Z location for recycling - away 
c     from the center of the target. For neutral ionization this option will
c     invoke getscross_approx to get a guesstimate of appropriate S CROSS 
c     values for the initial position of the ion. 
c
      init_pos_opt = 1
c
c -----------------------------------------------------------------------
c
c     TAG I25
c
c     fp_neut_opt - this option turns on the possibility of neutral
c     ionization within a far peripheral region - it only applies to 
c     wall elements not target ones. 
c     = 0 = off
c     = 1 = on       
c
      fp_neut_opt = 0
c
c -----------------------------------------------------------------------
c
c     TAG I26
c
c     fp_plasma_opt - this option specifies how a plasma will be
c     determined for the far periphery if a plasma is required. This
c     is particularly related to option I25.
c     = 0 = plasma from nearest associated real grid cell
c     = 1 = Te from specified input, ne from grid
c     = 2 = Both Te and Ne from specified input
c
      fp_plasma_opt = 0
c
c -----------------------------------------------------------------------
c
c     TAG I27
c
c     fp_te - default far periphery temperature if option in use (in eV)
c
      fp_te = 10.0
c
c -----------------------------------------------------------------------
c
c     TAG I28
c
c     fp_ne - default far periphery density if option in use (m-3)
c
      fp_ne = 1.0e18
c
c -----------------------------------------------------------------------
c
c     TAG I29 and I30
c
c     These are endpoints or corner points for the line/box ion injection
c     options ciopte=9 or ciopte=10. Default values are not 
c     meaningful but are assigned to values of 0.0. Code for 
c     ciopte 9 and 10 will check if all values are zero and will 
c     issue and error and exit. 
c
      cxsca = 0.0 
      cysca = 0.0 
      cxscb = 0.0 
      cyscb = 0.0 
c
c -----------------------------------------------------------------------
c
c     TAG I31
c
c     Far periphery transport flow option
c     0 - no flow in far periphery
c     1 - flow in far periphery is the same as associated ring
c     2 - flow in far periphery is specified as input
c  
      fp_flow_opt = 0
c
c     TAG I32 
c
c     Far periphery variable to hold velocity input 
c
      fp_flow_velocity_input = 0.0
c
c
c -----------------------------------------------------------------------
c
c     TAG P60 
c
c     ngradopt - new density gradient option to allow for density 
c     variation within the 2PM specification. Works in conjunction 
c     with temperature gradient options.
c
      ngradopt = 0  
c
c -----------------------------------------------------------------------
c
c     TAG P61
c
c     override_bg_velocity_opt - Option to override the background 
c     velocity which is calculated by the other SOL options.
c     Option 0 : off
c     Option 1 : Prescribed flow - using data from osmns28 
c     Option 2 : Recalculate background flow using the density 
c                from the SOL option and source data from EIRENE
c
c     Default value is OFF
c
      override_bg_velocity_opt = 0  
c
c -----------------------------------------------------------------------
c
c     TAG Q42
c
c     Tags Q43 and Q43 specify a temperature on a ring by ring 
c     basis which is then used instead of the target temperature for
c     calculating the target heat flux and sputtering yields. At the 
c     present time these values are used directly in NEUT. 
c 
c     TAG Q42  
c
c     - Two parameters - IR TE - for INNER JET/OUTER SONNET
c
      nsheath_vali = 0
      call rzero(sheath_vali,maxnrs*2) 

c
c -----------------------------------------------------------------------
c
c     TAG Q43  
c
c     - Two parameters - IR TE - for OUTER JET/INNER SONNET
c
      nsheath_valo = 0
      call rzero(sheath_valo,maxnrs*2) 

c
c -----------------------------------------------------------------------
c
c     TAG Q44 - Core plasma profiles as a function of PSIN  
c
c     - Five parameters - PSIN TE TI NE VB
c
c     Data array is allocatable and does not need initialization
c     (note: use allocatable_data module
c
      ncoreprofile = 0

c
c -----------------------------------------------------------------------
c
c     TAG R13
c
c     The following R-tags are enhancements to the detached plasma 
c     prescription that allow for the shape of the density profile
c     in this region to be more precisely specified. This was done
c     to make it possible to specify profiles in this region which
c     correspond reasonably well to the measured Divertor Thomson
c     profiles. By specifying values for the following inputs - ONLY 
c     supported on a ring by ring basis - the default standard 
c     behaviour is used elsewhere - the density profiles in region A
c     of the detached plasma prescription can be more generally 
c     defined. 
c
c     R13 - array of input data for INNER JET/OUTER Sonnet
c
c     TAG R13
c
c     - 9 parameters - IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B
c
      aux_ns21i = 0
      call rzero(aux_s21parmi,maxnrs*9)
c
c
c -----------------------------------------------------------------------
c
c     TAG R14 
c
c     R14 - array of input data for OUTER JET/INNER Sonnet
c
c     - 9 parameters - IR L1A L1B NR1A NR1B TER1A TER1B TIR1A TIR1B
c
      aux_ns21o = 0
      call rzero(aux_s21parmo,maxnrs*9)
c
c -----------------------------------------------------------------------
c
c     TAG S22
c
c     S22 - Option to turn on neutral V/A flag debugging
c           OFF by default
c
      debug_neutv = 0
c
c -----------------------------------------------------------------------
c
c     TAG S23
c
c     S23 - maximum energy/velocity used in the neutral velocity 
c           debugging options.
c
      debug_neutv_einmax = 50.0
c
c -----------------------------------------------------------------------
c
c     TAG S24
c
c     S24 - Number of bins to divide the velocity distribution when 
c           debugging neutral velocity 
c
      debug_neutv_nbins = 500
c
c -----------------------------------------------------------------------
C 
C     TAG T19 TO T27
c
c     Initialization of Reiser options to default values
c
      aswitch = 0
      sk11 = 0
      sk12 = 0
      sk13 = 0
      sd11 = 0
      sd12 = 0
      sd13 = 0
      coulomb_log = 15.0
      linearpeak = 0
c
c -----------------------------------------------------------------------
C 
c     TAG T28
c
c     T28 - pinch_loc_opt - this option specifies the region of the 
c                           grid where the radial velocity in
c                           pinchopt 4 will be applied.
c
c         = 0 = Entire grid excluding PFZ (default)
c         = 1 = main SOL only
c         = 2 = Entire grid including PFZ 
c               - this requires a sign change to the value assigned to 
c                 Vr (or Vpinch depending on terminology)
c         = 3 = Main SOL only above Xpoint region 
c               (i.e. In cells adjacent to adjacent to Xpoint and above)
c         = 4 = Main SOL above Xpoint region + core
c                                 
c
      pinch_loc_opt = 0   
c
c
c -----------------------------------------------------------------------
C 
c     TAG T29 
c
c     T29 - pinch_npdf, pinch_pdf - this option loads the probability
c           distribution function to be used when randomly 
c           determining the value of the pinch/radial velocity at 
c           each time step. 
c
c           At the present time no default PDF is loaded - a check must
c           be added at the end of the read routine to make sure
c           that a PDF has been specified if pinchopt 4 has been selected.
c
      pinch_npdf = 0
c 
c -----------------------------------------------------------------------
c
c     TAG T30 
c
c     T30 - pinch correlation time. A new radial velocity value will be 
c           chosen periodically based on the value of this quantity.
c           The default value of 0.0 will result in a new velocity 
c           being selected every time step. This value is specified 
c           in second. On JET it is typically 5 to 20 microseconds.
c
      pinch_correlation_time = 0.0  
c 
c -----------------------------------------------------------------------
c
c     TAG T31
c
c     T31 - Drift region - specifies the region to which poloidal 
c           drifts should be applied. 
c           1 - SOL + PFZ
c           2 - SOL only
c           3 - PFZ only
c           4 - CORE only
c
c           Other options can easily be added as needed - the default is
c           option 1. 
c
      drft_region = 1 
c
c -----------------------------------------------------------------------
c
c     TAG T32
c
c     T32 - Drift Mach Option - Detailed drift velocity input on a ring
c           ring basis is specified as a mach number to be multiplied
c           by the sound speed at the top of the torus for each ring
c
c           Option 0 : OFF
c           Option 1 : CS calculated from 2*Te
c           Option 2 : CS calculated from Te+Ti
c         
c           Default value is 0 - OFF - data is specified in terms of 
c                                      velocity
c
      drftvel_machopt=0
c
c -----------------------------------------------------------------------
c
c     TAG T33
c
c     T33 - Detailed specifications of data ring by ring - this is 
c           an array listing 
c                      ring number       velocity/mach
c           Data does not need to be specified for each ring - the 
c           default value will be applied instead.
c
c           The number of array elements is initialized to zero  
c
      ndrftvel = 0 
c
c -----------------------------------------------------------------------
c
c     TAG T34
c
c     T34 - S displacement in 2D resulting from a perpendicular step in
c           paramagnetic "Z" direction in 3D - this actually moves the 
c           particle onto an adjacent flux tube - however, since 
c           DIVIMP is 2D - the effect is to actually move the particle
c           onto an adjacent identical flux tube at a different S location
c           - thus effectively giving a net S displacement.
c
c           The first approximation to this is to use
c
c           ds = cross_step * Btor/Bpol
c
c           A value of 0 for this option is OFF 
c                      1 is ON
c
c           Default is shown below 
c
      dperpz_opt = 0      
      base_dperpz_step = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG T35 - related to poloidal drift options - T31,T32,T33
c
c     T35 - Drift region calculation option
c
c     Option 0: Input values are specified in terms of S
c     Option 1: Input values are specified in terms of P (poloidal distance)
c     Option 2: Input values are specified in terms of Z (single null only)
c
c     Default value is S para
c
      drft_distopt = 0
c
c -----------------------------------------------------------------------
c
c     TAG W01
c
c     wall_plasma_opt - this option specifies the algorithm to be used
c     to define the plasma conditions (if any) associated with each 
c     element of the wall. The target elements are associated with the
c     target plasma conditions and are not affected by this option.
c
c     0 = use plasma conditions in associated cell 
c     1 = linear decay beyond last ring + interpolation along ring
c     2 = exponential decay beyond last ring + interpolation along ring
c 
      wall_plasma_opt = 0
c
c
c -----------------------------------------------------------------------
c
c     TAG W02
c
c     wall_plasma_fact - this is a scale factor to be used in the scaling
c     algorithms defined by wall_plasma_opt = 1,2 - it is set to 0.1 m for
c     now. There are minimum values for wall plasma conditions specified
c     in the code. These could be moved to optional input if required in
c     the future.
c 
      wall_plasma_fact = 0.1
c
c -----------------------------------------------------------------------
c
c slmod begin
c     TAG 077
c
c     Data for additional neutral wall surfaces
c
      eirnasdat = 1
      eirasdat  = 0.0
      eirasdat(1,1) = 998.0  ! This triggers the automated wall clipping for EIRENE
c slmod end
c
c     End of intialization
c      
      return
      end
