module unstructured_input

  implicit none




contains

  !     -*-Fortran-*- 
  !
  !
  ! ======================================================================
  !
  ! subroutine: ValidateUnstructuredInput
  !
  !
  SUBROUTINE ValidateUnstructuredInput
    use mod_params
    use mod_comtor
    IMPLICIT none
    ! 
    !      include 'params'
    !      include 'comtor'
    !
    !
    !     No validation at present for LIM
    !

    RETURN
  END SUBROUTINE ValidateUnstructuredInput
  !
  ! ======================================================================
  !
  SUBROUTINE InitializeOUTUnstructuredInput
    use mod_params
    use mod_comtor
    use mod_coords
    use mod_comxyt
    IMPLICIT none
    !
    !     This routine sets the OUT related Unstructured inputs to their 
    !     default values. 
    !
    !
    !      INCLUDE 'params'
    !      include 'comtor'
    !      include 'coords'
    !
    ! -----------------------------------------------------------------------
    !
    !     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
    !
    ! -----------------------------------------------------------------------
    !
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG O01:
    !
    !     This option allows an absolute scaling factor for the LIM
    !     run results to be specified in the OUT routine. It's default
    !     value is zero.
    !
    new_absfac = 0.0
    !         
    ! -----------------------------------------------------------------------
    !
    !     TAG O99:
    !
    !     Net erosion plot scaling option
    !     0 = normal (which is particles/m /particle entering the system)
    !     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
    !     2 = Not yet implemented
    !
    erosion_scaling_opt = 0
    !      


    return
  end SUBROUTINE InitializeOUTUnstructuredInput
  !
  ! ======================================================================
  !
  SUBROUTINE InitializeUnstructuredInput
    use mod_params
    use iter_bm
    use variable_wall
    use yreflection
    use mod_comtor
    use mod_coords
    use mod_comxyt
    use mod_soledge_input
    use mod_sol22_input
    use mod_sol22_input_lim
    use mod_vtig
    use mod_diagvel_unstruc
    use mod_comt2
    use mod_lambda
    use mod_assign_plasma_input
    use allocatable_input_data
    IMPLICIT none
    !
    !     This routine sets the Unstructured inputs to their 
    !     default values. 
    !
    !
    !      INCLUDE 'params'
    !      include 'comtor'
    !      include 'coords'
    !
    !     Initialize sol22 unstructured input
    !      
    !      call sol22_initialize_unstructured_input
    !      
    ! -----------------------------------------------------------------------
    !
    !     TAG D07:
    !
    !     Sputter data option - this option specifies which set of 
    !     data will be used for calculating sputtering yields. The 
    !     available options are:
    !     1 - original LIM - Bohdansky
    !     2 - Eckstein IPP9/82 (1993)
    !     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
    !     4 - specified constant yield
    !     5 - As 3 except a custom routine is used for W.  
    !
    !
    csputopt = 3

    !
    ! -----------------------------------------------------------------------
    !
    !     TAG D08:
    !
    !     Chemical Sputter data option - this option specifies which set of 
    !     data will be used for calculating chemical sputtering yields. The 
    !     available options are:
    !     1 - Garcia-Rosales/Roth 1994
    !     2 - Garcia-Rosales/Roth 1996
    !     3 - JET 1 - Garcia-Rosales Formula EPS94
    !     4 - JET 2 - Popiesczyk EPS95
    !     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
    !     6 - JET 4 - Haasz - submitted to JNM Dec 1995
    !     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
    !     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
    !     9 - Constant yield with value set to const_yield
    !    10 - Modified Haasz 1997 - mass dependence made for H 
    !                             - lower yield at low temperatures
    !    11 - Modified Haasz 1997 - mass dependence
    !                             - reduced yield at low plasma temps 
    !                             - modified surface temp dependence
    !
    !     Set default to unmodified Haasz data 
    !
    cchemopt = 8
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG D23:
    !
    !     Constant value for use with CPUTOPT option = 4 
    !
    const_yield = 0.01
    !     
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG D39 : Alternate Sputter data specifier - used to select one of 
    !               several custom sputter datasets - usually based
    !               on different impact angles
    !
    !     Set to normal incidence data as default
    ! 
    extra_sputter_angle = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG D98 : Initial sputtered particle Y-coordinate
    !               This quantity over rides the Y-coordinate generated
    !               from the neutral following routines and replaces it with
    !               the specified value. 
    !
    init_y_coord = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG D99 : Sputter impact energy option -
    !               0 = LIM standard
    !               1 = remove particle kinetic energy term 
    !                   leaving sheath and temperature
    !
    !     Start counting down to avoid collision with DIVIMP tags
    !
    !     Set to standard as default
    !
    impact_energy_opt = 0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG I04 : Self sputtering option (added to match feature in DIVIMP)
    !               This allows self-sputtering to be specified independently
    !               of the sputter option
    !               0 = off
    !               1 = on (default to match existing behaviour for most options)
    !               2 = on (fixed constant self-sputtering energy specified)
    !
    cselfs = 1
    !
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L01: The "L" series is designated for LIM unstructured input
    !
    !              This option is only in effect for limiters with a
    !              specified poloidal extent for 3D LIM. Basically, particles
    !              crossing the L_con/2 point which are on a flux tube not
    !              connecting to the limiter will have their P coordinate
    !              adjusted so they will be placed on a flux tube which
    !              connects to the limiter. 
    !
    !              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
    !
    shear_short_circuit_opt=0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
    !              0 = off ... Wall = CAW
    !              1..n = on      CAW  = function of Y (option specifies function)
    !
    !              Option 1 = linear wall from CAW at ywall_start to 
    !                         CAW_MIN at YHALF (half way point between limiters)
    !
    lim_wall_opt = 0
    !
    !     TAG L03: LIM Wall shape option - starting Y value for revised wall value
    !
    ywall_start = 0.0

    !
    !     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
    !
    caw_min = HI
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L05 to L9: Limiter shape parameters for EDGE option 11 - ITER
    !     L05: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
    !     L06: rslot_setback - radial setback from LCFS at slot half width
    !     L07: bm_tor_wid - toroidal half width of the BM (blanket module)
    !     L08: slot_tor_wid - toroidal half witdth of the center slot of BM 
    !     L09: lambda_design - design SOL decay length
    !
    !     Additional parameters* for EDGE option 12 - ITER limiter shape modified for pitch angle
    !
    !     *see below
    !     L23: Bth_Bphi_ratio ... magnetic field ratio 
    !     L24: p_0_value ... central p value for this limiter surface slice
    !     L25: Rho_p_pol ... constant for g(p) function
    !     L26: R_ow ... R value for calculating shadow line
    !
    rtor_setback  = 0.07
    rslot_setback = 0.01
    bm_tor_wid    = 0.5954 
    slot_tor_wid  = 0.03
    lambda_design = 0.015
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L10 to L12: Inputs related to Y-Reflection option
    !
    !     L10: yreflection_opt: 0=off  1+ on : default = 0.0 or off
    !     L11: cmir_refl_lower - less than 0.0 - indicates location of Y<0 mirror
    !     L12: cmir_refl_upper - greater than 0.0 - indicates location of Y>0 mirror
    !     
    !     Default locations are also 0.0 for off - these MUST be specified to turn the 
    !     option on. 
    !
    !     yreflection_event_count - global initialization of counter to 0.0
    !
    !     These variables are in the yreflection module
    !
    yreflection_opt = 0
    cmir_refl_lower = 0.0
    cmir_refl_upper = 0.0
    yreflection_event_count = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L13: Calculate 3D Power emissions 
    !
    !     calc_3d_power = 0 (off)
    !                   = 1 (on)
    !
    !     Calculating the 3D versions of powls and lines which are stored
    !     in lim5 and tiz3 in the dmpout routine is time consuming. The 
    !     default of this option is to allow for calculation but when these
    !     data aren't needed the calculation can be turned off.
    !
    calc_3d_power = 1
    !
    ! -----------------------------------------------------------------------
    ! 
    !     TAG L14 and L15: Specified Sputtering flux and energy function
    !     
    !     L14: External flux option: extfluxopt
    !          0 = off (limiter plasma conditions used for surface fluxes)
    !          1 = external flux data specified in X
    !          2 = external flux data specified in Y
    !          3 = external flux data specified in D (distance along limiter surface)
    !
    !          Note: -X,-Y and -D data apply to the Y<0 side of the limiter
    !                 X, Y and  D data apply to the Y>0 side of the limiter
    !
    !     L14:
    !
    extfluxopt = 0
    !
    !     L15: 
    !
    !     External flux and energy function in either X,Y or D space
    !        <coord>   <flux m-2s-1>    <Eimpact eV>
    !     extfluxdata
    !
    nextfluxdata = 0 
    !
    !     Array not allocated at this point so can't be initialized      
    !     extfluxdata = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L16 to L22: Inputs related to X and Y absorption surfaces
    !
    !     
    !     L16 : xabsorb_opt : 0 = off 1=on
    !     L17 : xabsorb     : -CAW < Xabs < CA ... if X > Xabs particle removed  
    !
    !     L18 : yabsorb_opt: 0=off >0 on  Note: frame options not currently supported
    !     L19 : yabsorb1a  : first absorption surface (Y>0)
    !     L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
    !     L21 : yabsorb2a  : second absorption surface (Y<0)
    !     L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
    !     
    !     Default locations are also 0.0 for off - these MUST be specified to turn the 
    !     option on. 
    !
    !     yreflection_event_count - global initialization of counter to 0.0
    !
    !     These variables are in the yreflection module
    !
    xabsorb_opt = 0
    xabsorb = 1e6
    !
    !     Default values of Y absorber surfaces are at frame 0 ... coordinate 0.0
    !     There are no valid defaults for absorber surfaces
    !
    yabsorb_opt = 0
    yabsorb1a = 0.0
    yabsorb1_frame = 0
    yabsorb2a = 0.0
    yabsorb2_frame = 0
    !
    ! -----------------------------------------------------------------------
    !
    !     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
    !
    !     L23: Bth_Bphi_ratio ... magnetic field ratio 
    !     L24: p_0_value ... central p value for this limiter surface slice
    !     L25: Rho_p_pol ... constant for g(p) function
    !     L26: R_ow ... R value for calculating shadow line (outer wall radius)
    !
    bth_bphi_ratio = 0.154
    p_0_value = 0.0
    rho_p_pol = 3.0
    r_ow = 4.0   ! need to get a good default value
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L27: Optional self-sputtering yield modifier input
    !              X  SS_YMF(Y<0)   SS_YMF(Y>0)
    !
    ss_nymfs = 0
    !
    !     can't initialize dynamically allocated input      
    !     ss_cymfs = 1.0
    !
    ! -----------------------------------------------------------------------
    !
    !     TAG L28: Background plasma flow beyond limiter edge.
    !              This is mono-directional in the entire simulation - does
    !              not change sign. 
    !
    vpflow_3D = 0.0
    !
    !-----------------------------------------------------------------------
    !
    !     TAG L29 and L30: Minimum and maximum P values for 3D volumetric 
    !                      injection of initial ions ... cneuta = 3
    !                      Set default values to 0.0.    
    !
    !     L29
    p0s = 0.0
    !     L30
    p0l = 0.0
    !
    !     These quantities only make sense for 3D simulations
    !
    !     L31: P reflection option - reflect ions at P boundaries
    !
    preflect_opt = 0
    !
    !     L32: P reflection boundary value ... +/- specified quantity
    !          A value of 0.0 will set the reflection boundary to 
    !          ABS(PS(-MAXNPS))+CPSUB
    !
    preflect_bound = 0.0
    !
    !     
    !     --------------------------------------------------------
    !     
    !     L33: Specify P bin boundaries which will supercede P bin
    !          widths in the input file
    !
    npbins = 0
    !
    !     can't initialize dynamically allocated inputs      
    !     pbin_bnds = 0.0
    !
    !     Options for alternate SOL specification
    !     - use SOL12/13 modified code imported into LIM
    !     - upper and lower bounds are the specified absorbing surfaces
    !     - all ionization sources are analytical
    !     - allow for uniform power/particle which will be similar to the
    !       base LIM options
    !     - mid point in profiles will not necessarily coincide with LIMITER tip
    !     - modify plasma for flux tubes where the limiter is present
    !     - break point is 1/2 between limiter and absorbing boundaries 
    !
    !     L34 
    !     - multiple poloidal plasma zones (with or without limiter surfaces)
    !       SET OF P1,P2,ZONE,SURF limiter poloidal bounds
    !
    !      
    nsurf = 0
    !
    !     can't initialize dynamically allocated inputs
    !     surf_bnds = 0.0
    !
    !      
    !     L35 - colprobe3d - 0=off, 1=on 
    !           this input needs to occur prior to dynamic allocation of storage
    !     
    colprobe3d=0
    !     
    ! -----------------------------------------------------------------------
    !
    !     SOLEDGE related options in LIM (SOL 12+) 
    !
    !
    !     L36 - soledge solver option
    !
    cioptf_soledge = 13
    !     
    !     L37 - ionization option 0,1,4,5
    !
    csopt = 0
    !
    !     L38 - Radiation option 0,1,2,3
    !
    cpopt = 0
    !
    !     L39 - CSOLLS - length of ionization source  (or decay length) [fraction of field line length]
    !
    csolls = 0.5
    !
    !     L40 - CSOLLT - length of source region for csopt 4,5
    !
    csollt = 0.5
    !
    !     L41 - csollr - length or decay length of radiation source
    !      
    csollr = 0.5
    !
    !     L42 - csolfr - fractional strength of radiation source (cpopt 2,3)
    !      
    csolfr = 0.0
    !
    !     L43 - csolpr - absolute strength of radiation source (cpopt 0,1)      
    !
    csolpr = 0.0
    !
    !     L44 - cfiz - fractional split between linear and exponential ionization sources
    !                  in csopt 4,5      
    !
    cfiz = 0.0
    !
    !     L45 - sol13_padd - additional pressure loss to be applied
    !
    sol13_padd = 0.0
    !     
    !     L46 - sol13_pdist - distance over which to add pressure - 0.0 = change pressure at target by factor
    !
    sol13_pdist = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     L47 - Te profile radial shift
    !     L48 - Ti profile radial shift
    !     L49 - ne profile radial shift
    !
    te_prof_shift = 0.0
    ti_prof_shift = 0.0
    ne_prof_shift = 0.0
    !     
    !
    !     L50 - Te profile multiplier      
    !     L51 - Ti profile multiplier      
    !     L52 - ne profile multiplier      
    !
    te_prof_mult = 1.0
    ti_prof_mult = 1.0
    ne_prof_mult = 1.0

    !     sazmod
    !     L53 - Variable absorbing boundary (that affects plasma solution) switch
    !     L54 - Y location of left step in wall
    !     L55 - Y location of right step in wall
    !     L56 - X location of left step
    !     L57 - X location of right step
    !
    !
    !                 |                _                   |
    !                 |               | |              ____|  <-- xabsorb2a
    !  xabsorb1a -->  |_____          | |             |    ^          _step
    !      _step      ^     |         | |             |    |
    !                 |     |         | |             |    |
    !                 |     |_________|_|_____________|    |
    !                 |     ^          ^              ^   yabsorb2a
    !             yabsorb1a |          |              |   
    !                yabsorb1a_step   probe         yabsorb2a_step
    !
    vary_absorb = 0.0
    yabsorb1a_step = 0.0
    yabsorb2a_step = 0.0
    xabsorb1a_step = 0.0
    xabsorb2a_step = 0.0

    !     sazmod
    !     L58 - Switch for dividing radial diffusion into regions
    !     L59 - Radial diffusion coefficient in region 1
    !     L60 - Radial diffusion coefficient in region 2
    !     L61 - Radial diffusion coefficient in region 3
    !     L62 - Radial diffusion coefficient in region 4
    !
    !     Regions defined as follows:
    !
    !                 |        1       _         2         |
    !                 |               | |              ____|
    !                 |_____          | |             |              
    !                       |         | |             |    
    !                       |    3    | |    4        |    
    !                       |_________|_|_____________|    
    !                                                                     

    dperp_reg_switch = 0.0
    dperp_reg1 = 0.0
    dperp_reg2 = 0.0
    dperp_reg3 = 0.0
    dperp_reg4 = 0.0

    !     sazmod
    !     L63 - Skip writing to raw file to save time.

    skip_raw = 0

    !     sazmod
    !     L64 - Modify the plasma velocity in front of the step (right half
    !     only for now). Unsure of the physics basis here, but surely a
    !     shortening of L would affect the velocity somehow.

    mod_v_fact = 1.0

    !     sazmod
    !     Allow to choose from an exponential distribution in the Y direction
    !     in 3D injection option. X and P are still chosen uniformly between
    !     X0S, X0L and P0S, P0L, but Y will be chosen between Y0S, Y0L
    !     according to exp(lambda*Y), i.e. a positive lambda gives the
    !     greatest probability of Y to be Y0L.
    !     L65 - Switch for exponential.
    choose_exp = 0

    !     L66 - Lambda for exponential.
    choose_exp_lambda = 0
    !     This doesn't belong here since it isn't an input option.
    choose_exp_fact = 0

    !     L67 - Overall scaling to apply to the background plasma.
    vel_mod = 1.0

    !
    !
    !-----------------------------------------------------------------------
    !      
    !     SOL Option 22 overlay switch
    !     
    !     L68 - number of SOL22 overlay sections. 0 turns off SOL22
    !
    !        SOL22 related options - read a set of custom inputs
    !        
    !        Integer - number of SOL22 overlays 
    !         
    !        'desc'   X1  X2  PZ1  PZ2 SOLVER_OPT  'sol22_parameters.txt'  
    !      
    !     Each line contains a description, an [X1,X2] and [PZ1,PZ2] range, solver_option, file name
    !     to apply the plasma solver overlay and a filename containing the
    !     parameters for the region.          
    !
    !     SOLVER_OPT = 0 - use soledge for specified region (2PMs)
    !     SOLVER_OPT = 1 - use sol22 for specified region (1D fluid code)         
    !              
    !
    nsol22_opt = 0
    !
    !     L69 - SOLEDGE_OPT - option turns on SOL 12/13 two point model SOL
    !
    soledge_opt = 0
    !      
    !-----------------------------------------------------------------------
    !
    !     PARAMETER SPECIFICATION
    !      
    !     Optional specifications for LIM input paramters - the default values
    !     are set in mod_params_lim.f90 but they can be over ridden
    !     by optional inputs at the start of the input file.
    !
    !     These inputs must be at the start of the LIM input file and will
    !     over-ride the default values.       
    !      
    !     L70 MAXNXS
    !     L71 MAXNYS
    !     L72 MAXNPS
    !     L73 MAXIZS 
    !     L74 MAXIMP
    !     L75 MAXQXS
    !     L76 MAXQYS
    !     L77 MAXY3D  
    !     L78 MAXNTS
    !     L79 MAXINS
    !     L80 maxpzone
    !     
    !     DEFAULTS are specified in mod_params_lim.f (mod_params)
    !            
    !
    !-----------------------------------------------------------------------
    !
    !     L90 - load vTi profiles for use in calculating Ti
    !     there are a specified number of profiles
    !     each profile is applied over a given range of X (radial)
    !     These values will be used for calculating CTEMBSI in this radial
    !     range.
    !     Different profiles can be specified for different radial ranges         
    !      
    !
    n_vtig_blocks = 0
    !     
    !-----------------------------------------------------------------------
    !     
    !     L91 - load vb profiles for assigning background plasma velocity
    !     there are a specified number of profiles
    !     each profile is applied over a given range of X (radial)
    !     Different profiles can be specified for different radial ranges         
    !      
    !
    n_vb_blocks = 0
    !
    !-----------------------------------------------------------------------
    !
    !     L92: Debug velocity option - 0 off 1 on
    !
    debug_v_opt = 0
    !
    !-----------------------------------------------------------------------
    !
    !     L93: Velocity switch for forces - use velplasma and efiled instead
    !     of cvhys and ceys for calculating the forces on the particles
    !     in each cell. This allows spatially varying Efield and plasma
    !     velocity            
    !
    vel_efield_opt = 0

    !
    ! -----------------------------------------------------------------------
    !
    !     L94: X reflection option - 0 off 1 on
    !
    xreflection_opt = 0
    !
    ! -----------------------------------------------------------------------
    !
    !     L95: X reflection boundary - only used if the option is ON
    !
    xreflect_bound = 0.0
    !
    ! -----------------------------------------------------------------------
    !
    !     L96: vTiG Option - 0 = off
    !                        1 =vTiG specified (change Ti)
    !                        2 =vTiG specified (constant Ti imposed after)
    !                        3 =dTi/ds directly specified (constant Ti)
    !
    vtig_opt = 0
    !
    ! -----------------------------------------------------------------------
    !
    !     L97: vb Option - 0 = off
    !                      1 = on  vb specified
    !
    vb_opt = 0      
    !     
    ! -----------------------------------------------------------------------
    !
    !     L98: SOL22 background overlay switch
    !          sol22_opt = 0 off
    !          sol22_opt = 1 on
    !      
    sol22_opt = 0
    !     
    ! -----------------------------------------------------------------------
    !
    !     L99: sf_tau - scaling factor for all the characteristic times
    !                 = 1.0 by default which does not change the calculations
    sf_tau = 1.0      
    !
    !     
    ! -----------------------------------------------------------------------
    !
    !     LA0: sf_vdiff - scaling factor for the velocity diffusive step size
    !                 = 1.0 by default which does not change the calculations
    sf_vdiff = 1.0      

    ! -----------------------------------------------------------------------      
    !
    !     LA1: ctimsc_win - time injection window [ctimsc,ctimsc_win]
    !     particles launched with a random start time in
    !     this window - 0.0 turns the option off      
    !     
    ctimsc_win = 0.0
    ! -----------------------------------------------------------------------      
    !
    !     LA2: cdwelt_sum - option to either record particle position 
    !                       AT the specified times dwelfs * dwelts
    !                       option 0 - record at time t
    !                       option 1 - SUM particle positions over each time
    !     window from [dwelfs(i)->dwelfs(i+1)] * dwelts(iz)
    !
    !     default is 0 - or instantaneous snapshot     
    !
    !     jdemod - remove cdwelt_sum option functionality because it isn't
    !              physically meaningful.                  
    !     
    !      cdwelt_sum = 0
    !
    !
    ! -----------------------------------------------------------------------      
    !
    !     LA3: Lambda_vary_opt
    !      
    !     Define below - grouped with other lambda related inputs      
    !
    !
    ! -----------------------------------------------------------------------      
    !
    !     TAG LA4: pzone_opt - defines the poloidal zones that will be set up
    !              this may also over-write the value of maxpzone        
    !
    !     0 = off (2D) unless colprobe3D=1 in which case 2 zones - and pzone_opt reset to 1
    !     1 = simple collector probe - maxpzone = 2 - probe is in zone 2
    !     2 = maxpzone = 2*MAXNPS+1 - each poloidal row in 3D has its own plasma
    !         calculation - array extends from 1..2*MAXNPS+1 with 1->-NPS 
    !     3 = user specified pzones - maxpzone = pzone_opt (The maximum zone identifier
    !         allowed in the surface input is limited to maxpzone)
    !
    !     This option must appear in the input file before dynamic allocation     
    !     after parameter specification. 
    !     
    !     
    pzone_opt = 0
    !      
    ! -----------------------------------------------------------------------
    !
    !     LA5: solver_axis_opt - this option affects how the points at which
    !                            the plasma solution in soledge and sol22 are
    !                            calculated.
    !
    !                            0 = axis is evenly spaced with sufficient points
    !                                for reasonable spatial resolution. The results
    !                                are interpolated onto the actual coordinates
    !                            1 = the actual cell center coordinates are loaded 
    !                                into the solver and these are used for the  
    !                                coordinates in the plasma solver.
    !     
    solver_axis_opt = 0

    ! -----------------------------------------------------------------------
    !
    !     LA6: absorb_plasma - this option specifies plasma conditions along
    !                          the Y<0 and Y>0 absorbing surfaces if they are
    !                          different from the conditions along the limiter
    !
    !               Y<0          Y>0  
    !     X1 , ne, Te, Ti    X2, ne, Te, Ti
    !     X1, X2 must be in ascending order - only X1 is checked         
    !
    nabsorb_plasma = 0

    ! -----------------------------------------------------------------------
    !
    !     LA7: absorb_surf_data - this option specifies the location of the 
    !                             absorbing surfaces as a function of plasma
    !                             zone and an x-coordinate range.
    !                             These values overwrite the values loaded
    !                             into yabsorb_surf by other options. 
    !                             A value of 0.0 for either absorbing surface
    !                             is ignored. 
    !
    !     PZ  X1   X2   YABSORB_SURF(Y<0)   YABSORB_SURF(Y>0)
    !
    nabsorb_surf = 0
    !
    !-----------------------------------------------------------------------
    !
    !     LA8: sol22_default_filename - this option specifies the default
    !                                   sol22 paramter file to be used 
    !                                   if specific ones are not specified
    !
    sol22_default_filename = 'sol22-default.txt'
    !     
    !      
    ! -----------------------------------------------------------------------
    ! -----------------------------------------------------------------------
    !
    !    T47 Coulomb logarithm calculation options
    ! 
    !     0  = default = constant (default value = 15.0)
    !     1  = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)  [HC code - ]
    !          Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in
    !          Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.
    !
    !     2  = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)  [LIM code]
    !     3  = log(1.5e13 * t**(1.5) / sqrt(n))   [SOL22 PEI term]
    !    
    !     Default option is 2 in LIM
    !
    lambda_opt = 2
    !
    !    T48 Coulomb logarithm calculation options
    !     Coulomb logarithm constant value - default value is 15.0 - this allows
    !     specification of alternate constant values for option 0.         
    !
    !     dafault value = 15.0
    !
    lambda_val = 15.0
    !      
    !     LA3: Lambda_vary_opt
    !
    !     This is a lambda option specific to LIM. LIM allows for the use
    !     of single values of lambda for the entire simulation calculated
    !     from the n,t value at the inboard edge of the limiter tip. This does
    !     not vary spatially while the DIVIMP option is either a specified
    !     constant value or spatially varying based on local plasma conditions
    !
    !     Option 0: Constant single value for the entire plasma calculated using
    !               the formula from lambda_opt (default)
    !     
    !     1: Spatially varying lambda values based on local plasma
    !        conditions      
    !
    lambda_vary_opt = 0
    !     
    ! -----------------------------------------------------------------------
    !
    !     TAG Q26:
    !     
    !     Specification of a density multiplier (gradient) to be applied
    !     to the outboard region. 
    !
    !     READ IN DENSITY GRADIENT INFORMATION, IF ANY
    !     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
    !     OF THE DENSITY
    !
    !     Turned off by default 
    !
    nnbg = 0 
    !
    ! -----------------------------------------------------------------------
    !
    !
    !
    !
    ! -----------------------------------------------------------------------
    !
    !
    !     TAG M01:
    !
    !     Set the initial velocity angle of neutrals as 0.0 unless otherwise
    !     specified when CNEUTC=17
    !
    CIANGN=0.0
    !
    !
    ! -----------------------------------------------------------------------
    !
    !
    !     End of initialization 
    !
    return
  end SUBROUTINE InitializeUnstructuredInput
  !
  ! ======================================================================
  !
  SUBROUTINE ReadUnstructuredInput(line2)
    use mod_params
    use iter_bm
    use variable_wall
    use yreflection
    use mod_comtor
    use mod_coords
    use mod_comxyt
    use mod_soledge_input
    use mod_cadas
    use mod_sol22_input
    use mod_sol22_input_lim
    use allocate_arrays
    use mod_vtig
    use mod_diagvel_unstruc
    use mod_comt2
    use mod_lambda
    use mod_assign_plasma_input
    !use allocatable_input
    use allocatable_input_data
    use debug_options
    use mod_io
    IMPLICIT none

    CHARACTER line2*(*),LINE*72,TAG*3,COMENT*72,cdum1*1024
    REAL      R,vol,z1
    INTEGER   I,ir,ierr,i1,i2
    !
    integer :: izone
    !     
    ! jdemod - added variable to hold line read when calling RDG1 to get 
    !          ADAS data.
    !
    character line3*512

    !
    !      INCLUDE 'params'
    !      include 'comtor'
    !      include 'coords'
    !
    !
    !      COMMON /INPUTCHK/ inputflag
    !      INTEGER           inputflag(100)
    !
    !      COMMON /MACHCOM/ machine
    !      CHARACTER*64     machine
    !
    !      INTEGER    MAXTAG
    !      PARAMETER (MAXTAG=1000)
    !      COMMON /INPCHK/ ntaglist,taglist
    !      INTEGER     ntaglist,idum1
    !
    !      CHARACTER*3 taglist(MAXTAG)
    !
    !     Function declaration for TAG T29 
    !
    !      real vtest,res,vr_pdf_int
    !      external vr_pdf_int
    !

    integer in,is
    !
    WRITE(line,'(A72)') line2

    WRITE(TAG,'(A3)') LINE(3:5)

    ierr = 0
    !
    !     jdemod - the code that independently reads sol22 input from a
    !     specified file INCLUDES the possibility of
    !     unstructured input values. These routines utilize the same
    !     code as is used for the LIM input file - so if a tag
    !     passed here starts with a '2' it needs to be redirected to
    !     the code to read the SOL22 input values.       
    !      
    if (tag(1:1).eq.'2') then       
       call sol22_unstructured_input(tag,line,ierr)
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D07 : Physical Spuuter Data option
       !
    elseIF (tag(1:3).EQ.'D07') THEN
       !
       !
       !     Physical Sputter data option - this option specifies which set of 
       !     data will be used for calculating physical sputtering yields. The 
       !     available options are:
       !     1 - original LIM - Bohdansky
       !     2 - Eckstein IPP9/82 (1993)
       !     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
       !     4 - specified constant yield
       !     5 - As 3 except a custom routine is used for W.  
       !     6 - 2007 Eckstein data where available - otherwise option 3
       !
       !
       CALL ReadI(line,csputopt,1,6,'Sputter Data option')
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D08 : Chemical Spuuter Data option
       !
    ELSEIF (tag(1:3).EQ.'D08') THEN
       !
       !
       !     Chemical Sputter data option - this option specifies which set of 
       !     data will be used for calculating chemical sputtering yields. The 
       !     available options are:
       !     1 - Garcia-Rosales/Roth 1994
       !     2 - Garcia-Rosales/Roth 1996
       !     3 - JET 1 - Garcia-Rosales Formula EPS94
       !     4 - JET 2 - Popiesczyk EPS95
       !     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
       !     6 - JET 4 - Haasz - submitted to JNM Dec 1995
       !     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
       !     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
       !     9 - Constant yield with value set to const_yield
       !    10 - Modified Haasz 1997 - mass dependence made for H 
       !                             - lower yield at low temperatures
       !    11 - Modified Haasz 1997 - mass dependence
       !                             - reduced yield at low plasma temps 
       !                             - modified surface temp dependence
       !
       !
       CALL ReadI(line,cchemopt,1,11,'Chemical Sputter Data option')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D23 : Yield for Sputter option 4 - constant
       ! 
    elseif (tag(1:3).EQ.'D23') THEN
       !
       !     Constant value for use with CPUTOPT option = 4 
       !
       CALL ReadR(line,const_yield,0.0,1.0,'Specified constant yield')

       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D39 : Alternate Sputter data specifier - used to select one of 
       !               several custom sputter datasets - usually based
       !               on different impact angles
       ! 
    elseif (tag(1:3).EQ.'D39') THEN
       !
       !     Secondary sputter data specifier
       !
       CALL ReadR(line,extra_sputter_angle,-10.0,90.0,'Extra Sputter Angle Opt')
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D98 : Initial sputtered particle Y-coordinate
       !               This quantity over rides the Y-coordinate generated
       !               from the neutral following routines and replaces it with
       !               the specified value. 
       !
    elseif (tag(1:3).EQ.'D98') THEN
       !
       !     Sputtered/launched particle initial Y coordinate
       !
       CALL ReadR(line,init_y_coord,-HI,HI,'Specified Initial Y coordinate for all Launched Ions')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG D99 : Sputter impact energy option -
       !               0 = LIM standard
       !               1 = remove particle kinetic energy term 
       !                   leaving sheath and temperature
       !
       !     Start label counting down to avoid collision with DIVIMP tags
       !     Set to standard as default
       !
    elseif (tag(1:3).EQ.'D99') THEN
       !
       !     Impurity particle impact energy option
       !
       CALL ReadI(line,impact_energy_opt,0,1,'Impurity impact Energy calculation option')
       !
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L01: The "L" series is designated for LIM unstructured input
       !
       !              This option is only in effect for limiters with a
       !              specified poloidal extent for 3D LIM. Basically, particles
       !              crossing the L_con/2 point which are on a flux tube not
       !              connecting to the limiter will have their P coordinate
       !              adjusted so they will be placed on a flux tube which
       !              connects to the limiter. 
       !
       !              Shear_short_circuit_opt 0: OFF
       !              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
       !
    elseif (tag(1:3).EQ.'L01') THEN
       CALL ReadI(line,shear_short_circuit_opt,0,1,'Shear Short Circuit Option')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
       !              0 = off ... Wall = CAW
       !              1..n = on      CAW  = function of Y (option specifies function)
       !
       !              Option 1 = linear wall from CAW at ywall_start to 
       !                         CAW_MIN at YHALF (half way point between limiters)
       !
    elseif (tag(1:3).EQ.'L02') THEN
       CALL ReadI(line,lim_wall_opt,0,1,'LIM wall option')
       !
       !     TAG L03: LIM Wall shape option - starting Y value for revised wall value
       !
    elseif (tag(1:3).EQ.'L03') THEN
       CALL ReadR(line,ywall_start,0.0,HI,'Starting Y value for alternate wall')
       !
       !     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
       !
    elseif (tag(1:3).EQ.'L04') THEN
       CALL ReadR(line,caw_min,-HI,0.0,'Distance to wall at Yhalf')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L05 to L09: Limiter shape parameters for EDGE option 11 - ITER
       !     L03: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
       !     L04: rslot_setback - radial setback from LCFS at slot half width
       !     L05: bm_tor_wid - toroidal half width of the BM (blanket module)
       !     L06: slot_tor_wid - toroidal half witdth of the center slot of BM 
       !     L07: lambda_design - design SOL decay length
       !
       !     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
       !
       !     *see below
       !
       !     L23: Bth/Bphi ... magnetic field ratio 
       !     L24: p_0 ... central p value for this limiter surface slice
       !     L25: Rho_p_pol
       !     L26: R_ow ... R value for calculating shadow line
       !
    elseif (tag(1:3).EQ.'L05') THEN
       CALL ReadR(line,rtor_setback,0.0,HI,'Radial setback at BM edge (M)')
    elseif (tag(1:3).EQ.'L06') THEN
       CALL ReadR(line,rslot_setback,0.0,HI,'Radial setback at inner slot edge (M)')
    elseif (tag(1:3).EQ.'L07') THEN
       CALL ReadR(line,bm_tor_wid,0.0,HI,'Toroidal Half width of BM (M)')
    elseif (tag(1:3).EQ.'L08') THEN
       CALL ReadR(line,slot_tor_wid,0.0,HI,'Toroidal half width of slot (M)')
    elseif (tag(1:3).EQ.'L09') THEN
       CALL ReadR(line,lambda_design,0.0,HI,'Design decay length (M)')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L10 to L12: Inputs related to Y-Reflection option
       !
       !     L10: yreflection_opt: 0=off  1+ on : default = 0.0 or off
       !     L11: cmir_refl_lower - less than 0.0 - indicates location of Y<0 mirror
       !     L12: cmir_refl_upper - greater than 0.0 - indicates location of Y>0 mirror
       !     
       !     Default locations are also 0.0 for off - these MUST be specified to turn the 
       !     option on. 
       !
       !     L10: Y reflection option flag 
       !
    elseif (tag(1:3).EQ.'L10') THEN
       CALL ReadI(line,yreflection_opt,0,2,'Y-Reflection Option')
       !
       !     TAG L11: Y < 0 Reflection location specification
       !
    elseif (tag(1:3).EQ.'L11') THEN
       CALL ReadR(line,cmir_refl_lower,-HI,0.0,'Y-reflection: Y<0 reflection boundary')
       !
       !     TAG L12: Y > 0 Reflection location specification
       !
    elseif (tag(1:3).EQ.'L12') THEN
       CALL ReadR(line,cmir_refl_upper,0.0,HI,'Y-reflection: Y>0 reflection boundary')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L13: Calculate 3D Power emissions 
       !
       !     calc_3d_power = 0 (off)
       !                   = 1 (on)
       !
       !     Calculating the 3D versions of powls and lines which are stored
       !     in lim5 and tiz3 in the dmpout routine is time consuming. The 
       !     default of this option is to allow for calculation but when these
       !     data aren't needed the calculation can be turned off.
       !
    elseif (tag(1:3).EQ.'L13') THEN
       CALL ReadI(line,calc_3d_power,0,1,'3D power calculation option')
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L14 and L15: Specified Sputtering flux and energy function
       !     
       !     L14: External flux option: extfluxopt
       !          0 = off (limiter plasma conditions used for surface fluxes)
       !          1 = external flux data specified in X
       !          2 = external flux data specified in Y
       !          3 = external flux data specified in D (distance along limiter surface)
       !
       !          Note: -X,-Y and -D data apply to the Y<0 side of the limiter
       !                 X, Y and  D data apply to the Y>0 side of the limiter
       !
       !     L14:
       !
    elseif (tag(1:3).EQ.'L14') THEN
       CALL ReadI(line,extfluxopt,0,3,'External sputtering flux option')
       !
       !     L15: 
       !
       !     External flux and energy function in either X,Y or D space
       !        <coord>   <flux m-2s-1>    <Eimpact eV>
       !     extfluxdata
       !
    elseif (tag(1:3).EQ.'L15') THEN
       !
       !         CALL RDRARN(extfluxdata,nextfluxdata,
       CALL divrda(extfluxdata,nextfluxdata,MAXINS,-MACHHI,MACHHI,.TRUE.,0.0,MACHHI,2,'External sputtering flux data',IERR)


       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L16 to L22: Inputs related to X and Y absorption surfaces
       !
       !     
       !     L16 : xabsorb_opt : 0 = off 1=on
       !     L17 : xabsorb     : -CAW < Xabs < CA ... if X > Xabs particle removed  
       !
       !     L18 : yabsorb_opt: 0=off 1+=on  Note: frame options not currently supported
       !     L19 : yabsorb1a  : first absorption surface 
       !     L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
       !     L21 : yabsorb2a  : 
       !     L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
       !     
       !     These variables are in the yreflection module
       !
    elseif (tag(1:3).EQ.'L16') THEN
       !       L16 : xabsorb_opt : 0 = off 1=on
       CALL ReadI(line,xabsorb_opt,0,1,'X-absorption Option')
       !
    elseif (tag(1:3).EQ.'L17') THEN
       !       L17 : xabsorb : -CAW < Xabs < CA ... if X > Xabs particle removed  
       CALL ReadR(line,xabsorb,-HI,HI,'X absorption surface - X > Xabs')

    elseif (tag(1:3).EQ.'L18') THEN
       !       L18 : yabsorb_opt: 0=off N=number of absorbers (1 or 2 supported)
       CALL ReadI(line,yabsorb_opt,0,2,'Y-Absorption Option')
       !
    elseif (tag(1:3).EQ.'L19') THEN
       !       L19 : yabsorb1a  : first absorption surface 
       CALL ReadR(line,yabsorb1a,0.0,HI,'Y location of first absorber')

    elseif (tag(1:3).EQ.'L20') THEN
       !       L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
       CALL ReadI(line,yabsorb1_frame,-1000,1000,'Frame for absorption surface')
       !
    elseif (tag(1:3).EQ.'L21') THEN
       !       L21 : yabsorb2a  : 
       CALL ReadR(line,yabsorb2a,-HI,0.0,'Y location of second absorber')

    elseif (tag(1:3).EQ.'L22') THEN
       !       L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
       CALL ReadI(line,yabsorb2_frame,-1000,1000,'Frame for absorption surface')
       !
       ! -----------------------------------------------------------------------
       !
       !     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
       !
       !     L23: Bth_Bphi_ratio ... magnetic field ratio 
       !     L24: p_0_value ... central p value for this limiter surface slice
       !     L25: Rho_p_pol ... constant for g(p) function
       !     L26: R_ow ... R value for calculating shadow line
       !
    elseif (tag(1:3).EQ.'L23') THEN
       CALL ReadR(line,bth_bphi_ratio,0.0,HI,'Magnetic field ratio at limiter')
    elseif (tag(1:3).EQ.'L24') THEN
       CALL ReadR(line,p_0_value,-HI,HI,'Base p value for slice across limiter')
    elseif (tag(1:3).EQ.'L25') THEN
       CALL ReadR(line,rho_p_pol,0.0,HI,'Rho value for calculating g(p) function')
    elseif (tag(1:3).EQ.'L26') THEN
       CALL ReadR(line,r_ow,0.0,HI,'R value for calculating location of shadowline')

       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L27: Optional self-sputtering yield modifier input
       !              X  SS_YMF(Y<0)   SS_YMF(Y>0)
       !
    elseif (tag(1:3).EQ.'L27') THEN

       !         CALL RDRARN(ss_cymfs,ss_nymfs,
       CALL divrda(ss_cymfs,ss_nymfs,MAXINS,-MACHHI,MACHLO,.TRUE.,0.0,MACHHI,2,'SET of SS YMF X,M(Y<0),M(Y>0)',IERR)
       if (ss_cymfs(1,1).gt.ss_cymfs(ss_nymfs,1)) then 
          call errmsg('READIN PARAMETER: ','CYMFS DATA MUST BE ENTERED IN ASCENDING ORDER IN X')
          stop
       endif
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG L28: Background plasma flow beyond limiter edge.
       !              This is mono-directional in the entire simulation - does
       !              not change sign. 
       !
    elseif (tag(1:3).EQ.'L28') THEN
       CALL ReadR(line,vpflow_3d,-HI,HI,'SOL flow outside 3D limiter region')
       !-----------------------------------------------------------------------
       !
       !     TAG L29 and L30: Minimum and maximum P values for 3D volumetric 
       !                      injection of initial ions ... cneuta = 3
       !                      Set default values to 0.0.    
       !
    elseif (tag(1:3).EQ.'L29') THEN
       CALL ReadR(line,p0s,-HI,HI,'Min P value of injection region')
    elseif (tag(1:3).EQ.'L30') THEN
       CALL ReadR(line,p0l,-HI,HI,'Max P value of injection region region')
       !
       !
       !-----------------------------------------------------------------------
       !
       !     L31: P reflection option - 0=off, 1 =on (reflect at P boundaries)
       !
    elseif (tag(1:3).EQ.'L31') THEN
       CALL ReadI(line,preflect_opt,0,1,'P-Reflection Option')
       !
       !     L32: P reflection boundary +/- P bound for reflection events
       !          0.0 sets the P boundary to ABS(PS(-MAXNPS))+CPSUB
       !
    elseif (tag(1:3).EQ.'L32') THEN
       CALL ReadR(line,preflect_bound,0.0,HI,'P-reflection: +/- reflection boundary')

       !
       !
       !     L33: Specify P bin boundaries which will supercede P bin
       !          widths in the input file
       !       
    elseif (tag(1:3).eq.'L33') then
       CALL divrda(pbin_bnds,npbins,2*MAxnps+1,-MACHHI,MACHHI,.TRUE.,'*L33:Set of Pbin boundaries',IERR)
       !         
       !     L34 
       !     - multiple limiter boundaries
       !
    elseif (tag(1:3).eq.'L34') then

       call divrda(surf_bnds,nsurf,-MACHHI,MACHHI,.TRUE.,-MACHHI,MACHHI,3,&
            '*L34:SET OF P1,P2,ZONE,SURF limiter poloidal bounds',IERR)

       !
       if (allocated(surf_bnds)) then 
          !        Verify surface bounds to make sure that they do not overlap
          do izone = 1,nsurf
             !            write(0,'(a,i8,5(1x,g12.5))') 'Surf bounds:',izone,
             !     >                surf_bnds(izone,1),surf_bnds(izone,2),
             !     >                surf_bnds(izone,3),surf_bnds(izone,4)

             if (surf_bnds(izone,1).gt.surf_bnds(izone,2)) then
                write(error_message_data,'(a,i8,2(1x,g12.5))') &
                     'P1 > P2 for zone specification', &
                     izone,surf_bnds(izone,1),surf_bnds(izone,2)
                call errmsg('ReadUnstructuredInput: *L34:',error_message_data)

                stop 'ERROR: IN POLOIDAL PLASMA ZONE SPECIFICATION:1'
             endif

             if (izone.lt.nsurf) then 
                if (surf_bnds(izone,2).gt.surf_bnds(izone+1,1)) then

                   write(0,*) 'Incorrect limiter poloidal extents',izone,izone+1,nsurf
                   write(error_message_data,'(a,i8,2(1x,g12.5))') 'P2 (zone) > P1(zone+1)',&
                        izone,surf_bnds(izone,2),surf_bnds(izone+1,1)
                   call errmsg('ReadUnstructuredInput: *L34:',error_message_data)

                   stop 'ERROR IN POLOIDAL PLASMA ZONE SPECIFICATION:2'
                endif
             endif

             ! perform valid zone checking at the end of iolim

          end do
       endif
       !     
       !       L35 colprobe3d option ... 0 off 1 on
       !
    elseif (tag(1:3).eq.'L35') then 
       CALL ReadI(line,colprobe3d,0,1,'Switch to activate collector probe 3d plasma')
       !        
       !     L36 - cioptf_soledge - solver option
       !
    elseif (tag(1:3).EQ.'L36') THEN
       CALL ReadI(line,cioptf_soledge,11,14,'SOLEDGE base solver option')

       !        
       !     L37 - csopt - ionization option 0,1,4,5
       !
    elseif (tag(1:3).EQ.'L37') THEN
       CALL ReadI(line,csopt,0,5,'SOLEDGE ionization option')
       !
       !     L38 - cpopt - Radiation option 0,1,2,3
       !
    elseif (tag(1:3).EQ.'L38') THEN
       CALL ReadI(line,cpopt,0,3,'SOLEDGE radiation option')

       cpopt = 0
       !
       !     L39 - CSOLLS - length of ionization source  (or decay length) [fraction of field line length]
       !
    elseif (tag(1:3).EQ.'L39') THEN
       CALL ReadR(line,csolls,0.0,1.0,'Ionization source length or lambda')
       !
       !     L40 - CSOLLT - length of source region for csopt 4,5
       !
    elseif (tag(1:3).EQ.'L40') THEN
       CALL ReadR(line,csollt,0.0,1.0,'Ionization source length for opt 4,5')
       !
       !     L41 - csollr - length or decay length of radiation source
       !     
    elseif (tag(1:3).EQ.'L41') THEN
       CALL ReadR(line,csollr,0.0,1.0,'Radiation source length or lambda')
       !
       !     L42 - csolfr - fractional strength of radiation source (cpopt 2,3)
       !      
    elseif (tag(1:3).EQ.'L42') THEN
       CALL ReadR(line,csolfr,0.0,HI,'Radiation source fraction')
       !
       !     L43 - csolpr - absolute strength of radiation source (cpopt 0,1)      
       !
    elseif (tag(1:3).EQ.'L43') THEN
       CALL ReadR(line,csolpr,0.0,HI,'Radiation source strength (absolute)')
       !
       !     L44 - cfiz - fractional split between linear and exponential ionization sources
       !                  in csopt 4,5      
       !
    elseif (tag(1:3).EQ.'L44') THEN
       CALL ReadR(line,cfiz,0.0,1.0,'Ionization source fraction split opt 4,5')
       !
       !     L45 - sol13_padd - additional pressure loss to be applied
       !
    elseif (tag(1:3).EQ.'L45') THEN
       CALL ReadR(line,sol13_padd,0.0,HI,'Additional pressure source')
       !     
       !     L46 - sol13_pdist - distance over which to add pressure - 0.0 = change pressure at target by factor
       !
    elseif (tag(1:3).EQ.'L46') THEN
       CALL ReadR(line,sol13_pdist,0.0,1.0,'Additional pressure distribution distance')
       ! -----------------------------------------------------------------------
       !
       !     L47 - Te profile radial shift
       !     L48 - Ti profile radial shift
       !     L49 - ne profile radial shift
       !
    elseif (tag(1:3).EQ.'L47') THEN
       CALL ReadR(line,te_prof_shift,-HI,HI,'Te profile shift')
    elseif (tag(1:3).EQ.'L48') THEN
       CALL ReadR(line,ti_prof_shift,-HI,HI,'Ti profile shift')
    elseif (tag(1:3).EQ.'L49') THEN
       CALL ReadR(line,ne_prof_shift,-HI,HI,'ne profile shift')
       !     
       !
       !     L50 - Te profile multiplier      
       !     L51 - Ti profile multiplier      
       !     L52 - ne profile multiplier      
       !
    elseif (tag(1:3).EQ.'L50') THEN
       CALL ReadR(line,te_prof_mult,-HI,HI,'Te profile mult')
    elseif (tag(1:3).EQ.'L51') THEN
       CALL ReadR(line,ti_prof_mult,-HI,HI,'Ti profile mult')
    elseif (tag(1:3).EQ.'L52') THEN
       CALL ReadR(line,ne_prof_mult,-HI,HI,'ne profile mult')

       !     sazmod
       !     L53 - Variable absorbing boundary (that affects plasma solution) switch
       !     L54 - Y location of left step in wall
       !     L55 - Y location of right step in wall
       !     L56 - X location of left step
       !     L57 - X location of right step 

    elseif (tag(1:3).eq.'L53') then
       call ReadI(line, vary_absorb, 0, 1, 'vary absorbtion boundary')
    elseif (tag(1:3).eq.'L54') then
       call ReadR(line, yabsorb1a_step, -HI, HI, 'Y loc of left step')
    elseif (tag(1:3).eq.'L55') then
       call ReadR(line, yabsorb2a_step, -HI, HI, 'Y loc of right step')
    elseif (tag(1:3).eq.'L56') then
       call ReadR(line, xabsorb1a_step, -HI, HI, 'X loc of left step')
    elseif (tag(1:3).eq.'L57') then
       call ReadR(line, xabsorb2a_step, -HI, HI, 'X loc of right step')

       !     sazmod
       !     L58 - Switch for dividing radial diffusion into regions
       !     L59 - Radial diffusion coefficient in region 1
       !     L60 - Radial diffusion coefficient in region 2
       !     L61 - Radial diffusion coefficient in region 3
       !     L62 - Radial diffusion coefficient in region 4

    elseif (tag(1:3).eq.'L58') then
       call ReadI(line, dperp_reg_switch, 0, 1, 'Dperp regions switch')
    elseif (tag(1:3).eq.'L59') then
       call ReadR(line, dperp_reg1, -HI, HI,'Radial Dperp in region 1')
    elseif (tag(1:3).eq.'L60') then
       call ReadR(line, dperp_reg2, -HI, HI,'Radial Dperp in region 2')
    elseif (tag(1:3).eq.'L61') then
       call ReadR(line, dperp_reg3, -HI, HI,'Radial Dperp in region 3')
    elseif (tag(1:3).eq.'L62') then
       call ReadR(line, dperp_reg4, -HI, HI,'Radial Dperp in region 4')


       !     sazmod
       !     L63 - Skip writing to raw file to save time if you never use that file.

    elseif (tag(1:3).eq.'L63') then
       call ReadI(line, skip_raw, 0, 1, 'Skip raw output')

       !     sazmod
       !     L64 - Modify the plasma velocity in front of the step (right only for now)

    elseif (tag(1:3).eq.'L64') then
       call ReadR(line,mod_v_fact,-HI,HI,'Modify velocity in step reg')

       !     sazmod 
       !     L65 - Siwtch to choose from an exponential distribution for Y in
       !     the 3D injection option. X and P still uniformly chosen.
       !     L66 - The lambda value for the exponential distribution. Positive
       !     would put the max probability at Y0L, negative Y0S.
    elseif (tag(1:3).eq.'L65') then
       call ReadI(line,choose_exp,0,1,'Choose Y exp. switch')
    elseif (tag(1:3).eq.'L66') then
       call ReadR(line,choose_exp_lambda,-HI,HI,'Lambda for exp.')

       !     sazmod
       !     Multiply the plasma velocity everywhere by vel_mod. This could be
       !     justified if say you think the experiment was in a regime of 
       !     intermediate collisionality and thus the velocities are too fast
       !     when using just a simple SOL prescription (so a vel_mod < 1).         
    elseif (tag(1:3).eq.'L67') then
       call ReadR(line,vel_mod,-HI,HI,'Plasma velocity mod factor')


    elseif (tag(1:3).eq.'L68') then

       !
       !        SOL22 related options - read a set of custom inputs
       !        
       !        Integer - number of SOL22 overlays 
       !         
       !        'desc'   X1  X2  PZ1  PZ2 SOLVER_OPT  'sol22_parameters.txt'  
       !      
       !     Each line contains a description, an [X1,X2] and [PZ1,PZ2] range, solver_option, file name
       !     to apply the plasma solver overlay and a filename containing the
       !     parameters for the region.          
       !
       !     SOLVER_OPT = 0 - use soledge for specified region (2PMs)
       !     SOLVER_OPT = 1 - use sol22 for specified region (1D fluid code)         
       !         
       !     
       call ReadI(line,nsol22_opt,0,100,'NSOL22_OPT: Number of solver specifications')
       !
       !        Allocate storage to hold the options
       !         
       !     
       !           Read in SOL22 specifications
       !         
       call read_sol22_input(nsol22_opt)
       !
       !               
    elseif (tag(1:3).eq.'L69') then
       !
       !     L69 - SOLEDGE_OPT - option turns on SOL 12/13 two point model SOL
       !
       call pr_trace('Unstructured Input *L69:','soledge_opt')
       call ReadI(line,soledge_opt,0,1,'SOLEDGE_OPT: Turns on use of SOL 12,13')
       !     
       !-----------------------------------------------------------------------
       !
       !     PARAMETER SPECIFICATION
       !      
       !     Optional specifications for LIM input paramters - the default values
       !     are set in mod_params_lim.f90 but they can be over ridden
       !     by optional inputs at the start of the input file.
       !
       !     These inputs must be at the start of the LIM input file and will
       !     over-ride the default values.       
       !      
       !     L70 MAXNXS
       !     L71 MAXNYS
       !     L72 MAXNPS
       !     L73 MAXIZS 
       !     L74 MAXIMP
       !     L75 MAXQXS
       !     L76 MAXQYS
       !     L77 MAXY3D  
       !     L78 MAXNTS
       !     L79 MAXINS
       !     L80 maxpzone
       !     
    elseif (tag(1:3).eq.'L70') then
       write(0,*) 'IYEARH:',iyearh,maxnxs
       if (iyearh.eq.-1) then 
          call ReadI(line,maxnxs,1,1000000,'MAXNXS: Max X cells in mesh')
       else
          call errmsg('Unstructured Input L70','Attempt to change MAXNXS after allocation')
       endif

    elseif (tag(1:3).eq.'L71') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxnys,1,1000000,'MAXNYS: Max Y cells in mesh')
       else
          call errmsg('Unstructured Input L71','Attempt to change MAXNYS after allocation')
       endif

    elseif (tag(1:3).eq.'L72') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxnps,1,100000,'MAXNPS: Max P cells in mesh')
       else
          call errmsg('Unstructured Input L72','Attempt to change MAXNPS after allocation')
       endif

    elseif (tag(1:3).eq.'L73') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxizs,1,100,'MAXIZS: Max Imp charge states')
       else
          call errmsg('Unstructured Input L73','Attempt to change MAXIZS after allocation')
       endif

    elseif (tag(1:3).eq.'L74') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maximp,1,100000000,'MAXIMP: Max Impurity to follow')
       else
          call errmsg('Unstructured Input L74','Attempt to change MAXIMP after allocation')
       endif

    elseif (tag(1:3).eq.'L75') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxqxs,1,10000000,'MAXQXS: Max X cells fine mesh')
       else
          call errmsg('Unstructured Input L75','Attempt to change MAXQXS after allocation')
       endif

    elseif (tag(1:3).eq.'L76') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxqys,1,10000000,'MAXQYS: Max Y cells fine mesh')
       else
          call errmsg('Unstructured Input L76','Attempt to change MAXQYS after allocation')
       endif

    elseif (tag(1:3).eq.'L77') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxy3d,1,10000000,'MAXQYS: Max Y 3D cells')
       else
          call errmsg('Unstructured Input L77','Attempt to change MAXY3D after allocation')
       endif

    elseif (tag(1:3).eq.'L78') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxnts,1,100,'MAXNTS: Max Time cells')
       else
          call errmsg('Unstructured Input L78','Attempt to change MAXNTS after allocation')
       endif

    elseif (tag(1:3).eq.'L79') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxins,1,10000,'MAXINS: Max size some input')
       else
          call errmsg('Unstructured Input L79','Attempt to change MAXINS after allocation')
       endif

    elseif (tag(1:3).eq.'L80') then
       if (iyearh.eq.-1) then 
          call ReadI(line,maxpzone,1,10000,'MAXPZONE: Max poloidal zones')        
       else
          call errmsg('Unstructured Input L80','Attempt to change MAXPZONE after allocation')
       endif

       !
       !-----------------------------------------------------------------------
       !
       !     L90 - load vTi profiles for use in calculating Ti
       !     there are a specified number of profiles
       !     each profile is applied over a given range of X (radial)
       !     These values will be used for calculating CTEMBSI in this radial
       !     range.
       !     Different profiles can be specified for different radial ranges         
       !     data format
       !     'title   ' nblocks
       !     Xmin Xmax   Number_of_lines_in_block   poloidal_zone   y_zone
       !     Poloidal zone corresponds to 3D poloidal zones defined in LIM  - 0 is all zones
       !     Yzone is +/- 1.0 to apply to one side of the limiter or the other - 0 is all y zones
       !     
       !
    elseif (tag(1:3).eq.'L90') then
       call ReadI(line,n_vtig_blocks,0,10,'Number of blocks of vTiG data')        

       if (n_vtig_blocks.gt.0) then 
          call read_v_data(n_vtig_blocks,vtig_range,vtig_ndata,vtig_data,vtig_zones)
       endif
       !     
       !-----------------------------------------------------------------------
       !     
       !     L91 - load vb profiles for assigning background plasma velocity
       !     there are a specified number of profiles
       !     each profile is applied over a given range of X (radial)
       !     Different profiles can be specified for different radial ranges         
       !     data format
       !     'title   ' nblocks
       !     Xmin Xmax   Number_of_lines_in_block   poloidal_zone   y_zone
       !     Poloidal zone corresponds to 3D poloidal zones defined in LIM  - 0 is all zones
       !     Yzone is +/- 1.0 to apply to one side of the limiter or the other - 0 is all y zones
       !      
       !
    elseif (tag(1:3).eq.'L91') then
       call ReadI(line,n_vb_blocks,0,10,'Number of blocks of vb data')        

       if (n_vtig_blocks.gt.0) then 
          call read_v_data(n_vb_blocks,vb_range,vb_ndata,vb_data,vb_zones)
       endif

       !
       !-----------------------------------------------------------------------
       !
       !     L92: Debug velocity option - 0 off 1 on
       !
       !
       !           
    elseif (tag(1:3).eq.'L92') then
       call ReadI(line,debug_v_opt,0,1,'Debug velocity switch')        
       !call allocate_mod_diagvel
       !
       !-----------------------------------------------------------------------
       !
       !     L93: Velocity switch for forces - use velplasma and efiled instead
       !     of cvhys and ceys for calculating the forces on the particles
       !     in each cell. This allows spatially varying Efield and plasma
       !     velocity            
       !
    elseif (tag(1:3).eq.'L93') then
       call ReadI(line,vel_efield_opt,0,1,'Velocity/efield data option')        
       !                  
       !
       !-----------------------------------------------------------------------
       !
       !     L94: Xreflection_opt - 0 off 1 on
       !     
    elseif (tag(1:3).eq.'L94') then
       call ReadI(line,xreflection_opt,0,1,'Xreflection option')        
       !                  
       !-----------------------------------------------------------------------
       !
       !     L95: Xreflect_bound - location of X mirror
       !     
    elseif (tag(1:3).eq.'L95') then
       call ReadR(line,xreflect_bound,-HI,HI,'Plasma velocity mod factor')
       !
       !-----------------------------------------------------------------------
       !
       !     L96: vTiG Option - 0 = off
       !                        1 =vTiG specified (change Ti)
       !                        2 =vTiG specified (constant Ti imposed after)
       !                        3 =dTi/ds directly specified (constant Ti)
       !     
    elseif (tag(1:3).eq.'L96') then
       call ReadI(line,vtig_opt,0,3,'vTiG option')        
       !
       !-----------------------------------------------------------------------
       !
       !     L97: vb Option - 0 = off
       !                      1 = on  vb specified
       !
    elseif (tag(1:3).eq.'L97') then
       call ReadI(line,vb_opt,0,1,'vb option')        
       !
       !-----------------------------------------------------------------------
       !
       !     L98: SOL22 Option - 0 = off
       !                         1 = on 
       !
    elseif (tag(1:3).eq.'L98') then
       call ReadI(line,sol22_opt,0,1,'SOL22 background overlay switch')        

       !     
       ! -----------------------------------------------------------------------
       !
       !     L99: sf_tau - scaling factor for all the characteristic times
       !                 = 1.0 by default which does not change the calculations
    elseif (tag(1:3).EQ.'L99') THEN
       CALL ReadR(line,sf_tau,0.0,HI,'Characteristic times scaling factor')
       !     
       ! -----------------------------------------------------------------------
       !
       !     LA0: sf_vdiff - scaling factor for velocity diffusion step size
       !                 = 1.0 by default which does not change the calculations
    elseif (tag(1:3).EQ.'LA0') THEN
       CALL ReadR(line,sf_vdiff,0.0,HI,'Velocity diffusion step size scaling factor')
       ! -----------------------------------------------------------------------      
       !
       !     LA1: ctimsc_win - time injection window [ctimsc,ctimsc_win]
       !     particles launched with a random start time in
       !     this window - 0.0 turns the option off      
       !     
    elseif (tag(1:3).EQ.'LA1') THEN
       CALL ReadR(line,ctimsc_win,-HI,HI,'End of particle injection time window')

       !
       ! -----------------------------------------------------------------------      
       !
       !     jdemod - remove cdwelt_sum option functionality because it isn't
       !              physically meaningful.                  
       !
       !     TAG LA2: cdwelt_sum - option to either record particle position 
       !                       AT the specified times dwelfs * dwelts
       !                       option 0 - record at time t
       !                       option 1 - SUM particle positions over each time
       !     window from [dwelfs(i)->dwelfs(i+1)] * dwelts(iz)
       !
       !     default is 0 - or instantaneous snapshot     
       !     
       !      elseif (tag(1:3).EQ.'LA2') THEN
       !        call ReadI(line,cdwelt_sum,0,1,
       !     >                 'Time dependent data collection option')        
       !        

       !
       ! -----------------------------------------------------------------------      
       !
       !     TAG LA3: lambda_vary_opt - grouped below with other lambda options
       !     
       ! -----------------------------------------------------------------------      
       !
       !     TAG LA4: pzone_opt - defines the poloidal zones that will be set up
       !              this may also over-write the value of maxpzone        
       !
       !     0 = off (2D) unless colprobe3D=1 in which case 2 zones
       !     1 = simple collector probe - maxpzone = 2 - probe is in zone 2
       !     2 = maxpzone = 2*NPS+1 - each poloidal row in 3D has its own plasma
       !         calculation
       !     3 = user specified pzones - maxpzone = nsurf (set when these options are
       !         read in. If nsurf = 0 and pzone_opt=3 after input - pzone_opt=0 is
       !         assigned.         
       !
       !     This option must appear in the input file before dynamic allocation     
       !     after parameter specification. 
       !     
       !     
    elseif (tag(1:3).EQ.'LA4') THEN
       call ReadI(line,pzone_opt,0,3,'Option defining poloidal plasma zones')        
       ! Ideally maxpzone should be set to match this but gets into
       ! issue with allocation unless I extract all the input
       ! arrays to a separate module and delay allocation until
       ! after the input file is read.
       ! Might be doable using rdrarn_alloc and/or alloc versions
       ! of all array reading routines. 
       ! But need to make sure these do not get re-allocated after being
       ! read in. 
       ! setting maxpzone is done at the beginning of iolim.f : readin
       !     
       !     -----------------------------------------------------------------------
       !
       !     LA5: solver_axis_opt - this option affects how the points at which
       !                            the plasma solution in soledge and sol22 are
       !                            calculated.
       !
       !                            0 = axis is evenly spaced with sufficient points
       !                                for reasonable spatial resolution. The results
       !                                are interpolated onto the actual coordinates
       !                            1 = the actual cell center coordinates are loaded 
       !                                into the solver and these are used for the  
       !                                coordinates in the plasma solver.
       !     
    elseif (tag(1:3).EQ.'LA5') THEN
       call ReadI(line,solver_axis_opt,0,1,'Defines the axis to be used in soledge and sol22 solvers')        
       ! -----------------------------------------------------------------------
       !
       !     LA6: absorb_plasma - this option specifies plasma conditions along
       !                          the Y<0 and Y>0 absorbing surfaces if they are
       !                          different from the conditions along the limiter
       !
       !               Y<0          Y>0  
       !     X1 , ne, Te, Ti    X2, ne, Te, Ti
       !     X1, X2 must be in ascending order - only X1 is checked         
       !     
    elseif (tag(1:3).EQ.'LA6') THEN
       call divrda(absorb_plasma,nabsorb_plasma,-MACHHI,MACHHI,.FALSE.,-MACHHI,MACHHI,7,&
            '*LA6:Plasma conditions along absorbing surfaces',IERR)
       if (nabsorb_plasma.gt.0) then 
          do in = 1,nabsorb_plasma
             write(0,*) 'ABS_PLASMA:',(absorb_plasma(in,is),is=1,8)
          end do
       endif

       !
       ! -----------------------------------------------------------------------
       !
       !     LA7: absorb_surf_data - this option specifies the location of the 
       !                             absorbing surfaces as a function of plasma
       !                             zone and an x-coordinate range.
       !                             These values overwrite the values loaded
       !                             into yabsorb_surf by other options. 
       !                             A value of 0.0 for either absorbing surface
       !                             is ignored. 
       !
       !     PZ1  PZ2  X1   X2   YABSORB_SURF(Y<0)   YABSORB_SURF(Y>0)
       !
    elseif (tag(1:3).EQ.'LA7') THEN
       call divrda(absorb_surf_data,nabsorb_surf,-MACHHI,MACHHI,.TRUE.,-MACHHI,MACHHI,5,&
            '*LA7:Absorbing surfaces by zone and X-range',IERR)

       !
       !-----------------------------------------------------------------------
       !
       !     LA8: sol22_default_filename - this option specifies the default
       !                                   sol22 paramter file to be used 
       !                                   if specific ones are not specified
       !
    elseif (tag(1:3).EQ.'LA8') THEN
       CALL RDC(sol22_default_filename,'SOL22 Default parameter file name', IERR)                                   
       !
       !        
       ! -----------------------------------------------------------------------
       !        
       ! -----------------------------------------------------------------------
       !
       !    T47 Coulomb logarithm calculation options
       ! 
       !     0  = default = constant (default value = 15.0)
       !     1  = 30.0 - 0.5 * LOG(ni) + 1.5 * LOG(ti)  [HC code - ]
       !          Originally in Sivukhin, D.V., Coulomb collisions in a fully ionized plasma in
       !          Review of Plasma Physics (Consultation Bureau, New York, 1966) Vol. 4, p.88.
       !
       !     2  = 17.3 - 0.5*LOG(n/1.0E20) + 1.5*LOG(t/1000.0)  [LIM code]
       !     3  = log(1.5e13 * t**(1.5) / sqrt(n))   [SOL22 PEI term]
       !    
    ELSEIF (tag(1:3).EQ.'T47') THEN
       CALL ReadI(line,lambda_opt,0,3,'Coulomb logarithm calc opt')
       !
       !    T48 Coulomb logarithm calculation options
       !     Coulomb logarithm constant value - default value is 15.0 - this allows
       !     specification of alternate constant values for option 0.         
       !
       !        
    ELSEIF (tag(1:3).EQ.'T48') THEN
       CALL ReadR(line,lambda_val,0.0,HI,'Coulomb logarithm const val')
       !        
       !     LA3: Lambda variation option
       !     0 = one value for the entire plasma (based on lambda_opt and plasma
       !         conditions from the inboard limiter tip         
       !     1 = varies depending on local conditions
       !
    ELSEIF (tag(1:3).EQ.'LA3') THEN
       CALL ReadI(line,lambda_vary_opt,0,1,'Lambda spatial variation option')
       !
       !-----------------------------------------------------------------------
       !
       !     TAG Q26:
       !
       !     Specification of a density multiplier (gradient) to be applied
       !     to the outboard region. 
       !
    elseif (tag(1:3).EQ.'Q26') THEN
       !
       !
       !     READ IN DENSITY GRADIENT INFORMATION, IF ANY
       !     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
       !     OF THE DENSITY
       !

       !         CALL RDRARN(MNBG,NNBG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
       CALL divrda(MNBG,NNBG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,&
            1,'SET OF Y,MNB VALUES',IERR)
       !
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG M01:
       !
       !     Set initial velocity angle of neutrals as 0.0 unless otherwise
       !     specified when CNEUTC=17
       !
    ELSEIF (TAG(1:3).EQ.'M01') THEN
       CALL ReadR(line,CIANGN,-180.0,180.0,'Initial neutral velocity angle')       
       ciangn = ciangn * degrad  
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG O01: 
       !
       !     Specify an alternate absolute factor in LIM
       !
    elseif (tag(1:3).eq.'O01') then 
       !
       !     This option allows an absolute scaling factor for the LIM
       !     run results to be specified in the OUT routine. It's default
       !     value is zero.
       !
       !        CALL ReadR(line,new_absfac,0.0,HI,
       !     >                   'Imposed ABSFAC in OUT')
       CALL ReadDP(line,new_absfac,0.0,HI,'Imposed ABSFAC in OUT')
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG O99:
       !
       !     Net erosion plot scaling option
       !     0 = normal (which is particles/m /particle entering the system)
       !     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
       !     2 = Not yet implemented
       !
    elseif (tag(1:3).eq.'O99') then 
       !      
       CALL ReadI(line,erosion_scaling_opt,0,2,'Erosion Scaling Option')
       !         
       !
       ! -----------------------------------------------------------------------
       !
       !     TAG does not match available input - signal ERROR and EXIT
       !

    ELSE
       CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
    endif
    !

    RETURN
    !
    ! There is an error:
    !
99  WRITE(stddbg,'(5X,3A)') 'LINE = "',trim(line),'"'
    WRITE(stddbg,'(5X,3A)') 'TAG  = "',trim(tag) ,'"'
    WRITE(stderr,'(5X,3A)') 'LINE = "',trim(line),'"'
    WRITE(stderr,'(5X,3A)') 'TAG  = "',trim(tag) ,'"'
    WRITE(stderr,*) '    DIVIMP HALTED'
    STOP
  END SUBROUTINE ReadUnstructuredInput



end module unstructured_input
