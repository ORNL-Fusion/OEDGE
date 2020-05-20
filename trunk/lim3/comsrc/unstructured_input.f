c     -*-Fortran-*- 
c
c
c ======================================================================
c
c subroutine: ValidateUnstructuredInput
c
c
      SUBROUTINE ValidateUnstructuredInput
      use mod_params
      use mod_comtor
      IMPLICIT none
c 
c      include 'params'
c      include 'comtor'
c
c
c     No validation at present for LIM
c

      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE InitializeOUTUnstructuredInput
      use mod_params
      use mod_comtor
      use mod_coords
      use mod_comxyt
      IMPLICIT none
c
c     This routine sets the OUT related Unstructured inputs to their 
c     default values. 
c
c
c      INCLUDE 'params'
c      include 'comtor'
c      include 'coords'
c
c -----------------------------------------------------------------------
c
c     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
c
c -----------------------------------------------------------------------
c
c
c -----------------------------------------------------------------------
c
c     TAG O01:
c
c     This option allows an absolute scaling factor for the LIM
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c
      new_absfac = 0.0
c         
c -----------------------------------------------------------------------
c
c     TAG O99:
c
c     Net erosion plot scaling option
c     0 = normal (which is particles/m /particle entering the system)
c     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
c     2 = Not yet implemented
c
      erosion_scaling_opt = 0
c      

      
      return
      end
c
c ======================================================================
c
      SUBROUTINE InitializeUnstructuredInput
      use mod_params
      use iter_bm
      use variable_wall
      use yreflection
      use mod_comtor
      use mod_coords
      use mod_comxyt
      use mod_soledge
      IMPLICIT none
c
c     This routine sets the Unstructured inputs to their 
c     default values. 
c
c
c      INCLUDE 'params'
c      include 'comtor'
c      include 'coords'
c
c -----------------------------------------------------------------------
c
c     TAG D07:
c
c     Sputter data option - this option specifies which set of 
c     data will be used for calculating sputtering yields. The 
c     available options are:
c     1 - original LIM - Bohdansky
c     2 - Eckstein IPP9/82 (1993)
c     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
c     4 - specified constant yield
c     5 - As 3 except a custom routine is used for W.  
c
c
      csputopt = 3

c
c -----------------------------------------------------------------------
c
c     TAG D08:
c
c     Chemical Sputter data option - this option specifies which set of 
c     data will be used for calculating chemical sputtering yields. The 
c     available options are:
c     1 - Garcia-Rosales/Roth 1994
c     2 - Garcia-Rosales/Roth 1996
c     3 - JET 1 - Garcia-Rosales Formula EPS94
c     4 - JET 2 - Popiesczyk EPS95
c     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
c     6 - JET 4 - Haasz - submitted to JNM Dec 1995
c     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
c     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
c     9 - Constant yield with value set to const_yield
c    10 - Modified Haasz 1997 - mass dependence made for H 
c                             - lower yield at low temperatures
c    11 - Modified Haasz 1997 - mass dependence
c                             - reduced yield at low plasma temps 
c                             - modified surface temp dependence
c
c     Set default to unmodified Haasz data 
c
      cchemopt = 8
c
c -----------------------------------------------------------------------
c
c     TAG D23:
c
c     Constant value for use with CPUTOPT option = 4 
c
      const_yield = 0.01
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
c -----------------------------------------------------------------------
c
c     TAG D98 : Initial sputtered particle Y-coordinate
c               This quantity over rides the Y-coordinate generated
c               from the neutral following routines and replaces it with
c               the specified value. 
c
      init_y_coord = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG D99 : Sputter impact energy option -
c               0 = LIM standard
c               1 = remove particle kinetic energy term 
c                   leaving sheath and temperature
c
c     Start counting down to avoid collision with DIVIMP tags
c
c     Set to standard as default
c
      impact_energy_opt = 0
c
c -----------------------------------------------------------------------
c
c     TAG I04 : Self sputtering option (added to match feature in DIVIMP)
c               This allows self-sputtering to be specified independently
c               of the sputter option
c               0 = off
c               1 = on (default to match existing behaviour for most options)
c               2 = on (fixed constant self-sputtering energy specified)
c
      cselfs = 1
c
c
c -----------------------------------------------------------------------
c
c     TAG L01: The "L" series is designated for LIM unstructured input
c
c              This option is only in effect for limiters with a
c              specified poloidal extent for 3D LIM. Basically, particles
c              crossing the L_con/2 point which are on a flux tube not
c              connecting to the limiter will have their P coordinate
c              adjusted so they will be placed on a flux tube which
c              connects to the limiter. 
c
c              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
c
      shear_short_circuit_opt=0
cg
c -----------------------------------------------------------------------
c
c     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
c              0 = off ... Wall = CAW
c              1..n = on      CAW  = function of Y (option specifies function)
c
c              Option 1 = linear wall from CAW at ywall_start to 
c                         CAW_MIN at YHALF (half way point between limiters)
c
      lim_wall_opt = 0
c
c     TAG L03: LIM Wall shape option - starting Y value for revised wall value
c
      ywall_start = 0.0

c
c     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
c
      caw_min = HI
c
c -----------------------------------------------------------------------
c
c     TAG L05 to L9: Limiter shape parameters for EDGE option 11 - ITER
c     L05: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
c     L06: rslot_setback - radial setback from LCFS at slot half width
c     L07: bm_tor_wid - toroidal half width of the BM (blanket module)
c     L08: slot_tor_wid - toroidal half witdth of the center slot of BM 
c     L09: lambda_design - design SOL decay length
c
c     Additional parameters* for EDGE option 12 - ITER limiter shape modified for pitch angle
c
c     *see below
c     L23: Bth_Bphi_ratio ... magnetic field ratio 
c     L24: p_0_value ... central p value for this limiter surface slice
c     L25: Rho_p_pol ... constant for g(p) function
c     L26: R_ow ... R value for calculating shadow line
c
      rtor_setback  = 0.07
      rslot_setback = 0.01
      bm_tor_wid    = 0.5954 
      slot_tor_wid  = 0.03
      lambda_design = 0.015
c
c -----------------------------------------------------------------------
c
c     TAG L10 to L12: Inputs related to Y-Reflection option
c
c     L10: yreflection_opt: 0=off  1+ on : default = 0.0 or off
c     L11: cmir_refl_lower - less than 0.0 - indicates location of Y<0 mirror
c     L12: cmir_refl_upper - greater than 0.0 - indicates location of Y>0 mirror
c     
c     Default locations are also 0.0 for off - these MUST be specified to turn the 
c     option on. 
c
c     yreflection_event_count - global initialization of counter to 0.0
c
c     These variables are in the yreflection module
c
      yreflection_opt = 0
      cmir_refl_lower = 0.0
      cmir_refl_upper = 0.0
      yreflection_event_count = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG L13: Calculate 3D Power emissions 
c
c     calc_3d_power = 0 (off)
c                   = 1 (on)
c
c     Calculating the 3D versions of powls and lines which are stored
c     in lim5 and tiz3 in the dmpout routine is time consuming. The 
c     default of this option is to allow for calculation but when these
c     data aren't needed the calculation can be turned off.
c
      calc_3d_power = 1
c
c -----------------------------------------------------------------------
c 
c     TAG L14 and L15: Specified Sputtering flux and energy function
c     
c     L14: External flux option: extfluxopt
c          0 = off (limiter plasma conditions used for surface fluxes)
c          1 = external flux data specified in X
c          2 = external flux data specified in Y
c          3 = external flux data specified in D (distance along limiter surface)
c
c          Note: -X,-Y and -D data apply to the Y<0 side of the limiter
c                 X, Y and  D data apply to the Y>0 side of the limiter
c
c     L14:
c
      extfluxopt = 0
c
c     L15: 
c
c     External flux and energy function in either X,Y or D space
c        <coord>   <flux m-2s-1>    <Eimpact eV>
c     extfluxdata
c
      nextfluxdata = 0
      extfluxdata = 0.0
c
c -----------------------------------------------------------------------
c
c     TAG L16 to L22: Inputs related to X and Y absorption surfaces
c
c
c     L16 : xabsorb_opt : 0 = off 1=on
c     L17 : xabsorb     : -CAW < Xabs < CA ... if X > Xabs particle removed  
c
c     L18 : yabsorb_opt: 0=off N=number of absorbers (1 or 2 supported)
c     L19 : yabsorb1a  : first absorption surface 
c     L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
c     L21 : yabsorb2a  : 
c     L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
c     
c     Default locations are also 0.0 for off - these MUST be specified to turn the 
c     option on. 
c
c     yreflection_event_count - global initialization of counter to 0.0
c
c     These variables are in the yreflection module
c
      xabsorb_opt = 0
      xabsorb = 1e6
c
c     Default values of Y absorber surfaces are at frame 0 ... coordinate 0.0
c     There are no valid defaults for absorber surfaces
c
      yabsorb_opt = 0
      yabsorb1a = 0.0
      yabsorb1_frame = 0
      yabsorb2a = 0.0
      yabsorb2_frame = 0
c
c -----------------------------------------------------------------------
c
c     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
c
c     L23: Bth_Bphi_ratio ... magnetic field ratio 
c     L24: p_0_value ... central p value for this limiter surface slice
c     L25: Rho_p_pol ... constant for g(p) function
c     L26: R_ow ... R value for calculating shadow line (outer wall radius)
c
      bth_bphi_ratio = 0.154
      p_0_value = 0.0
      rho_p_pol = 3.0
      r_ow = 4.0   ! need to get a good default value
c
c -----------------------------------------------------------------------
c
c     TAG L27: Optional self-sputtering yield modifier input
c              X  SS_YMF(Y<0)   SS_YMF(Y>0)
c
      ss_nymfs = 0
      ss_cymfs = 1.0
c
c -----------------------------------------------------------------------
c
c     TAG L28: Background plasma flow beyond limiter edge.
c              This is mono-directional in the entire simulation - does
c              not change sign. 
c
      vpflow_3D = 0.0
c
c-----------------------------------------------------------------------
c
c     TAG L29 and L30: Minimum and maximum P values for 3D volumetric 
c                      injection of initial ions ... cneuta = 3
c                      Set default values to 0.0.    
c
c     L29
      p0s = 0.0
c     L30
      p0l = 0.0
c
c     These quantities only make sense for 3D simulations
c
c     L31: P reflection option - reflect ions at P boundaries
c
      preflect_opt = 0
c
c     L32: P reflection boundary value ... +/- specified quantity
c          A value of 0.0 will set the reflection boundary to 
c          ABS(PS(-MAXNPS))+CPSUB
c
      preflect_bound = 0.0
c
c     
c     --------------------------------------------------------
c     
c     L33: Specify P bin boundaries which will supercede P bin
c          widths in the input file
c
      npbins = 0
      pbin_bnds = 0.0
c
c     Options for alternate SOL specification
c     - use SOL12/13 modified code imported into LIM
c     - upper and lower bounds are the specified absorbing surfaces
c     - all ionization sources are analytical
c     - allow for uniform power/particle which will be similar to the
c       base LIM options
c     - mid point in profiles will not necessarily coincide with LIMITER tip
c     - modify plasma for flux tubes where the limiter is present
c     - break point is 1/2 between limiter and absorbing boundaries 
c   
C
c     L34 
c     - multiple limiter boundaries
c
c      
      nsurf = 0
      surf_bnds = 0.0

c      
c     L35 - colprobe3d - 0=off, 1=on 
c        
      colprobe3d=0
c     
c -----------------------------------------------------------------------
c 
c     SOLEDGE related options in LIM (SOL 12+) 
c
c
c     L36 - soledge solver option
c
      cioptf_soledge = 13
c     
c     L37 - ionization option 0,1,4,5
c
      csopt = 0
c
c     L38 - Radiation option 0,1,2,3
c
      cpopt = 0
c
c     L39 - CSOLLS - length of ionization source  (or decay length) [fraction of field line length]
c
      csolls = 0.5
c
c     L40 - CSOLLT - length of source region for csopt 4,5
c
      csollt = 0.5
c
c     L41 - csollr - length or decay length of radiation source
c      
      csollr = 0.5
c
c     L42 - csolfr - fractional strength of radiation source (cpopt 2,3)
c      
      csolfr = 0.0
c
c     L43 - csolpr - absolute strength of radiation source (cpopt 0,1)      
c
      csolpr = 0.0
c
c     L44 - cfiz - fractional split between linear and exponential ionization sources
c                  in csopt 4,5      
c
      cfiz = 0.0
c
c     L45 - sol13_padd - additional pressure loss to be applied
c
      sol13_padd = 0.0
c     
c     L46 - sol13_pdist - distance over which to add pressure - 0.0 = change pressure at target by factor
c
      sol13_pdist = 0.0
c
c -----------------------------------------------------------------------
c
c     L47 - Te profile radial shift
c     L48 - Ti profile radial shift
c     L49 - ne profile radial shift
c
      te_prof_shift = 0.0
      ti_prof_shift = 0.0
      ne_prof_shift = 0.0
c     
c
c     L50 - Te profile multiplier      
c     L51 - Ti profile multiplier      
c     L52 - ne profile multiplier      
c
      te_prof_mult = 1.0
      ti_prof_mult = 1.0
      ne_prof_mult = 1.0
c      
c -----------------------------------------------------------------------
c
c     TAG Q26:
c
c     Specification of a density multiplier (gradient) to be applied
c     to the outboard region. 
c
C     READ IN DENSITY GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE DENSITY
c
c     Turned off by default 
C
      nnbg = 0 
c
c -----------------------------------------------------------------------
c
C
C
C
C -----------------------------------------------------------------------
C
C
C     TAG M01:
C
C     Set the initial velocity angle of neutrals as 0.0 unless otherwise
C     specified when CNEUTC=17
C
      CIANGN=0.0
C
c
c -----------------------------------------------------------------------
c
c
c     End of initialization 
c
      return
      end  
c
c ======================================================================
c
      SUBROUTINE ReadUnstructuredInput(line2)
      use mod_params
      use iter_bm
      use variable_wall
      use yreflection
      use mod_comtor
      use mod_coords
      use mod_comxyt
      use mod_soledge
      IMPLICIT none

      CHARACTER line2*(*),LINE*72,TAG*3,COMENT*72,cdum1*1024
      REAL      R,vol,z1
      INTEGER   I,ir,ierr,i1,i2
c
      integer :: izone
c     
c jdemod - added variable to hold line read when calling RDG1 to get 
c          ADAS data.
c
      character line3*512

c
c      INCLUDE 'params'
c      include 'comtor'
c      include 'coords'
c
c
c      COMMON /INPUTCHK/ inputflag
c      INTEGER           inputflag(100)
c
c      COMMON /MACHCOM/ machine
c      CHARACTER*64     machine
c
c      INTEGER    MAXTAG
c      PARAMETER (MAXTAG=1000)
c      COMMON /INPCHK/ ntaglist,taglist
c      INTEGER     ntaglist,idum1
c
c      CHARACTER*3 taglist(MAXTAG)
c
c     Function declaration for TAG T29 
c
c      real vtest,res,vr_pdf_int
c      external vr_pdf_int
c

      integer in
c
      WRITE(line,'(A72)') line2

      WRITE(TAG,'(A3)') LINE(3:5)

      ierr = 0

c
c -----------------------------------------------------------------------
c
c     TAG D07 : Physical Spuuter Data option
c
      IF (tag(1:3).EQ.'D07') THEN
c
c
c     Physical Sputter data option - this option specifies which set of 
c     data will be used for calculating physical sputtering yields. The 
c     available options are:
c     1 - original LIM - Bohdansky
c     2 - Eckstein IPP9/82 (1993)
c     3 - Eckstein IPP9/82 + Adjustments from Garcia/Rosales-Roth 1996
c     4 - specified constant yield
c     5 - As 3 except a custom routine is used for W.  
c     6 - 2007 Eckstein data where available - otherwise option 3
c
c
        CALL ReadI(line,csputopt,1,6,'Sputter Data option')
c
c
c -----------------------------------------------------------------------
c
c     TAG D08 : Chemical Spuuter Data option
c
      ELSEIF (tag(1:3).EQ.'D08') THEN
c
c
c     Chemical Sputter data option - this option specifies which set of 
c     data will be used for calculating chemical sputtering yields. The 
c     available options are:
c     1 - Garcia-Rosales/Roth 1994
c     2 - Garcia-Rosales/Roth 1996
c     3 - JET 1 - Garcia-Rosales Formula EPS94
c     4 - JET 2 - Popiesczyk EPS95
c     5 - JET 3 - Vietzke (from Phys.Processes.of.Interaction.Fus.Plas.with.Solids)
c     6 - JET 4 - Haasz - submitted to JNM Dec 1995
c     7 - JET 5 - Roth & Garcia-Rosales - submitted to JNM March 1996
c     8 - JET 6 - Haasz 1997 - Brian Mech's PhD thesis data
c     9 - Constant yield with value set to const_yield
c    10 - Modified Haasz 1997 - mass dependence made for H 
c                             - lower yield at low temperatures
c    11 - Modified Haasz 1997 - mass dependence
c                             - reduced yield at low plasma temps 
c                             - modified surface temp dependence
c
c
        CALL ReadI(line,cchemopt,1,11,'Chemical Sputter Data option')
c
c -----------------------------------------------------------------------
c
c     TAG D23 : Yield for Sputter option 4 - constant
c 
      elseif (tag(1:3).EQ.'D23') THEN
c
c     Constant value for use with CPUTOPT option = 4 
c
        CALL ReadR(line,const_yield,0.0,1.0,'Specified constant yield')

c
c -----------------------------------------------------------------------
c
c     TAG D39 : Alternate Sputter data specifier - used to select one of 
c               several custom sputter datasets - usually based
c               on different impact angles
c 
      elseif (tag(1:3).EQ.'D39') THEN
c
c     Secondary sputter data specifier
c
        CALL ReadR(line,extra_sputter_angle,-10.0,90.0,
     >             'Extra Sputter Angle Opt')
c
c
c -----------------------------------------------------------------------
c
c     TAG D98 : Initial sputtered particle Y-coordinate
c               This quantity over rides the Y-coordinate generated
c               from the neutral following routines and replaces it with
c               the specified value. 
c
      elseif (tag(1:3).EQ.'D98') THEN
c
c     Sputtered/launched particle initial Y coordinate
c
        CALL ReadR(line,init_y_coord,-HI,HI,
     >       'Specified Initial Y coordinate for all Launched Ions')
c
c -----------------------------------------------------------------------
c
c     TAG D99 : Sputter impact energy option -
c               0 = LIM standard
c               1 = remove particle kinetic energy term 
c                   leaving sheath and temperature
c
c     Start label counting down to avoid collision with DIVIMP tags
c     Set to standard as default
c
      elseif (tag(1:3).EQ.'D99') THEN
c
c     Impurity particle impact energy option
c
        CALL ReadI(line,impact_energy_opt,0,1,
     >                'Impurity impact Energy calculation option')
c
c
c
c -----------------------------------------------------------------------
c
c     TAG L01: The "L" series is designated for LIM unstructured input
c
c              This option is only in effect for limiters with a
c              specified poloidal extent for 3D LIM. Basically, particles
c              crossing the L_con/2 point which are on a flux tube not
c              connecting to the limiter will have their P coordinate
c              adjusted so they will be placed on a flux tube which
c              connects to the limiter. 
c
c              Shear_short_circuit_opt 0: OFF
c              Shear_short_circuit_opt 1: P = CPCO * ( 2*ran-1) = (-CPCO,+CPCO) 
c
      elseif (tag(1:3).EQ.'L01') THEN
        CALL ReadI(line,shear_short_circuit_opt,0,1,
     >                'Shear Short Circuit Option')
c
c -----------------------------------------------------------------------
c
c     TAG L02: LIM Wall shape option - allow CAW to vary as a function of Y
c              0 = off ... Wall = CAW
c              1..n = on      CAW  = function of Y (option specifies function)
c
c              Option 1 = linear wall from CAW at ywall_start to 
c                         CAW_MIN at YHALF (half way point between limiters)
c
      elseif (tag(1:3).EQ.'L02') THEN
        CALL ReadI(line,lim_wall_opt,0,1,'LIM wall option')
c
c     TAG L03: LIM Wall shape option - starting Y value for revised wall value
c
      elseif (tag(1:3).EQ.'L03') THEN
        CALL ReadR(line,ywall_start,0.0,HI,
     >               'Starting Y value for alternate wall')
c
c     TAG L04: LIM Wall shape option - value of CAW reached at midpoint (YHALF)
c
      elseif (tag(1:3).EQ.'L04') THEN
        CALL ReadR(line,caw_min,-HI,0.0,'Distance to wall at Yhalf')
c
c -----------------------------------------------------------------------
c
c     TAG L05 to L09: Limiter shape parameters for EDGE option 11 - ITER
c     L03: rtor_setback - radial setback from LCFS at the toroidal half width of the BM
c     L04: rslot_setback - radial setback from LCFS at slot half width
c     L05: bm_tor_wid - toroidal half width of the BM (blanket module)
c     L06: slot_tor_wid - toroidal half witdth of the center slot of BM 
c     L07: lambda_design - design SOL decay length
c
c     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
c
c     *see below
c
c     L23: Bth/Bphi ... magnetic field ratio 
c     L24: p_0 ... central p value for this limiter surface slice
c     L25: Rho_p_pol
c     L26: R_ow ... R value for calculating shadow line
c
      elseif (tag(1:3).EQ.'L05') THEN
        CALL ReadR(line,rtor_setback,0.0,HI,
     >          'Radial setback at BM edge (M)')
      elseif (tag(1:3).EQ.'L06') THEN
        CALL ReadR(line,rslot_setback,0.0,HI,
     >                    'Radial setback at inner slot edge (M)')
      elseif (tag(1:3).EQ.'L07') THEN
        CALL ReadR(line,bm_tor_wid,0.0,HI,
     >                    'Toroidal Half width of BM (M)')
      elseif (tag(1:3).EQ.'L08') THEN
        CALL ReadR(line,slot_tor_wid,0.0,HI,
     >                    'Toroidal half width of slot (M)')
      elseif (tag(1:3).EQ.'L09') THEN
        CALL ReadR(line,lambda_design,0.0,HI,
     >                    'Design decay length (M)')
c
c -----------------------------------------------------------------------
c
c     TAG L10 to L12: Inputs related to Y-Reflection option
c
c     L10: yreflection_opt: 0=off  1+ on : default = 0.0 or off
c     L11: cmir_refl_lower - less than 0.0 - indicates location of Y<0 mirror
c     L12: cmir_refl_upper - greater than 0.0 - indicates location of Y>0 mirror
c     
c     Default locations are also 0.0 for off - these MUST be specified to turn the 
c     option on. 
c
c     L10: Y reflection option flag 
c
      elseif (tag(1:3).EQ.'L10') THEN
        CALL ReadI(line,yreflection_opt,0,1,'Y-Reflection Option')
c
c     TAG L11: Y < 0 Reflection location specification
c
      elseif (tag(1:3).EQ.'L11') THEN
        CALL ReadR(line,cmir_refl_lower,-HI,0.0,
     >               'Y-reflection: Y<0 reflection boundary')
c
c     TAG L12: Y > 0 Reflection location specification
c
      elseif (tag(1:3).EQ.'L12') THEN
        CALL ReadR(line,cmir_refl_upper,0.0,HI,
     >               'Y-reflection: Y>0 reflection boundary')
c
c -----------------------------------------------------------------------
c
c     TAG L13: Calculate 3D Power emissions 
c
c     calc_3d_power = 0 (off)
c                   = 1 (on)
c
c     Calculating the 3D versions of powls and lines which are stored
c     in lim5 and tiz3 in the dmpout routine is time consuming. The 
c     default of this option is to allow for calculation but when these
c     data aren't needed the calculation can be turned off.
c
      elseif (tag(1:3).EQ.'L13') THEN
        CALL ReadI(line,calc_3d_power,0,1,'3D power calculation option')
c
c
c -----------------------------------------------------------------------
c 
c     TAG L14 and L15: Specified Sputtering flux and energy function
c     
c     L14: External flux option: extfluxopt
c          0 = off (limiter plasma conditions used for surface fluxes)
c          1 = external flux data specified in X
c          2 = external flux data specified in Y
c          3 = external flux data specified in D (distance along limiter surface)
c
c          Note: -X,-Y and -D data apply to the Y<0 side of the limiter
c                 X, Y and  D data apply to the Y>0 side of the limiter
c
c     L14:
c
      elseif (tag(1:3).EQ.'L14') THEN
        CALL ReadI(line,extfluxopt,0,3,'External sputtering flux'//
     >         ' option')
c
c     L15: 
c
c     External flux and energy function in either X,Y or D space
c        <coord>   <flux m-2s-1>    <Eimpact eV>
c     extfluxdata
c
      elseif (tag(1:3).EQ.'L15') THEN
c
         CALL RDRARN(extfluxdata,nextfluxdata,
     >               MAXINS,-MACHHI,MACHHI,.TRUE.,0.0,MACHHI,            
     >               2,'External sputtering flux data',IERR)


c
c -----------------------------------------------------------------------
c
c     TAG L16 to L22: Inputs related to X and Y absorption surfaces
c
c
c     L16 : xabsorb_opt : 0 = off 1=on
c     L17 : xabsorb     : -CAW < Xabs < CA ... if X > Xabs particle removed  
c
c     L18 : yabsorb_opt: 0=off N=number of absorbers (1 or 2 supported)
c     L19 : yabsorb1a  : first absorption surface 
c     L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
c     L21 : yabsorb2a  : 
c     L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
c     
c     These variables are in the yreflection module
c
      elseif (tag(1:3).EQ.'L16') THEN
c       L16 : xabsorb_opt : 0 = off 1=on
        CALL ReadI(line,xabsorb_opt,0,1,'X-absorption Option')
c
      elseif (tag(1:3).EQ.'L17') THEN
c       L17 : xabsorb : -CAW < Xabs < CA ... if X > Xabs particle removed  
        CALL ReadR(line,xabsorb,-HI,HI,
     >               'X absorption surface - X > Xabs')

      elseif (tag(1:3).EQ.'L18') THEN
c       L18 : yabsorb_opt: 0=off N=number of absorbers (1 or 2 supported)
        CALL ReadI(line,yabsorb_opt,0,2,'Y-Absorption Option')
c
      elseif (tag(1:3).EQ.'L19') THEN
c       L19 : yabsorb1a  : first absorption surface 
        CALL ReadR(line,yabsorb1a,-HI,HI,
     >               'Y location of first absorber')

      elseif (tag(1:3).EQ.'L20') THEN
c       L20 : yabsorb1_frame : frame reference for surface - 0 for no reflection cases
        CALL ReadI(line,yabsorb1_frame,-1000,1000,
     >                'Frame for absorption surface')
c
      elseif (tag(1:3).EQ.'L21') THEN
c       L21 : yabsorb2a  : 
        CALL ReadR(line,yabsorb2a,-HI,HI,
     >               'Y location of second absorber')

      elseif (tag(1:3).EQ.'L22') THEN
c       L22 : yabsorb2_frame : frame reference for surface - 0 for no reflection cases
        CALL ReadI(line,yabsorb2_frame,-1000,1000,
     >                'Frame for absorption surface')
c
c -----------------------------------------------------------------------
c
c     Additional parameters for EDGE option 12 - ITER limiter shape modified for pitch angle
c
c     L23: Bth_Bphi_ratio ... magnetic field ratio 
c     L24: p_0_value ... central p value for this limiter surface slice
c     L25: Rho_p_pol ... constant for g(p) function
c     L26: R_ow ... R value for calculating shadow line
c
      elseif (tag(1:3).EQ.'L23') THEN
        CALL ReadR(line,bth_bphi_ratio,0.0,HI,
     >          'Magnetic field ratio at limiter')
      elseif (tag(1:3).EQ.'L24') THEN
        CALL ReadR(line,p_0_value,-HI,HI,
     >           'Base p value for slice across limiter')
      elseif (tag(1:3).EQ.'L25') THEN
        CALL ReadR(line,rho_p_pol,0.0,HI,
     >            'Rho value for calculating g(p) function')
      elseif (tag(1:3).EQ.'L26') THEN
        CALL ReadR(line,r_ow,0.0,HI,
     >            'R value for calculating location of shadowline')

c
c -----------------------------------------------------------------------
c
c     TAG L27: Optional self-sputtering yield modifier input
c              X  SS_YMF(Y<0)   SS_YMF(Y>0)
c
      elseif (tag(1:3).EQ.'L27') THEN

         CALL RDRARN(ss_cymfs,ss_nymfs,
     >               MAXINS,-MACHHI,MACHLO,.TRUE.,0.0,MACHHI,            
     >               2,'SET of SS YMF X,M(Y<0),M(Y>0)',IERR)
         if (ss_cymfs(1,1).gt.ss_cymfs(ss_nymfs,1)) then 
            call errmsg('READIN PARAMETER: ','CYMFS DATA MUST'//
     >       ' BE ENTERED IN ASCENDING ORDER IN X')
            stop
         endif
c
c -----------------------------------------------------------------------
c
c     TAG L28: Background plasma flow beyond limiter edge.
c              This is mono-directional in the entire simulation - does
c              not change sign. 
c
      elseif (tag(1:3).EQ.'L28') THEN
        CALL ReadR(line,vpflow_3d,-HI,HI,
     >            'SOL flow outside 3D limiter region')
c-----------------------------------------------------------------------
c
c     TAG L29 and L30: Minimum and maximum P values for 3D volumetric 
c                      injection of initial ions ... cneuta = 3
c                      Set default values to 0.0.    
c
      elseif (tag(1:3).EQ.'L29') THEN
        CALL ReadR(line,p0s,-HI,HI,
     >            'Min P value of injection region')
      elseif (tag(1:3).EQ.'L30') THEN
        CALL ReadR(line,p0l,-HI,HI,
     >            'Max P value of injection region region')
c
c
c-----------------------------------------------------------------------
c
c     L31: P reflection option - 0=off, 1 =on (reflect at P boundaries)
c
      elseif (tag(1:3).EQ.'L31') THEN
        CALL ReadI(line,preflect_opt,0,1,'P-Reflection Option')
c
c     L32: P reflection boundary +/- P bound for reflection events
c          0.0 sets the P boundary to ABS(PS(-MAXNPS))+CPSUB
c
      elseif (tag(1:3).EQ.'L32') THEN
        CALL ReadR(line,preflect_bound,0.0,HI,
     >               'P-reflection: +/- reflection boundary')

c
c
c     L33: Specify P bin boundaries which will supercede P bin
c          widths in the input file
c       
      elseif (tag(1:3).eq.'L33') then
         CALL RDRAR(pbin_bnds,npbins,
     >        2*MAxnps+1,-MACHHI,MACHHI,.TRUE.,
     >       'Set of Pbin boundaries',IERR)

c         
c     L34 
c     - multiple limiter boundaries
c
      elseif (tag(1:3).eq.'L34') then 
         CALL RDRARN(surf_bnds,nsurf,max_nsurf,
     >               -MACHHI,MACHHI,.TRUE.,0.0,MACHHI,            
     >               1,'SET OF P1,P2 limiter poloidal bounds',IERR)

c        Verify surface bounds to make sure tha they do not overlap
         do izone = 1,nsurf
c            write(0,'(a,i8,2(1x,g12.5))') 'Surf bounds:',izone,
c     >                surf_bnds(izone,1),surf_bnds(izone,2)

            if (surf_bnds(izone,1).gt.surf_bnds(izone,2)) then
               write(0,*) 'Incorrect limiter poloidal extents',
     >                 izone,nsurf
               stop 'ERROR IN LIMITER POLOIDAL EXTENT SPECIFICATION'
            endif

            if (izone.lt.nsurf) then 
               if (surf_bnds(izone,2).gt.surf_bnds(izone+1,1)) then
                  write(0,*) 'Incorrect limiter poloidal extents',
     >                 izone,izone+1,nsurf
                  stop 'ERROR IN LIMITER POLOIDAL EXTENT SPECIFICATION'
               endif
            endif
        end do
c
c       L35 colprobe3d option ... 0 off 1 on
c
      elseif (tag(1:3).eq.'L35') then 
        CALL ReadI(line,colprobe3d,0,1,
     >                'Switch to activate collector probe 3d plasma')
c        
c     L36 - cioptf_soledge - solver option
c
      elseif (tag(1:3).EQ.'L36') THEN
        CALL ReadI(line,cioptf_soledge,12,14,
     >                'SOLEDGE base solver option')

c        
c     L37 - csopt - ionization option 0,1,4,5
c
      elseif (tag(1:3).EQ.'L37') THEN
        CALL ReadI(line,csopt,0,5,
     >                'SOLEDGE ionization option')
c
c     L38 - cpopt - Radiation option 0,1,2,3
c
      elseif (tag(1:3).EQ.'L38') THEN
        CALL ReadI(line,cpopt,0,3,
     >                'SOLEDGE radiation option')

      cpopt = 0
c
c     L39 - CSOLLS - length of ionization source  (or decay length) [fraction of field line length]
c
      elseif (tag(1:3).EQ.'L39') THEN
        CALL ReadR(line,csolls,0.0,1.0,
     >               'Ionization source length or lambda')
c
c     L40 - CSOLLT - length of source region for csopt 4,5
c
      elseif (tag(1:3).EQ.'L40') THEN
        CALL ReadR(line,csollt,0.0,1.0,
     >               'Ionization source length for opt 4,5')
c
c     L41 - csollr - length or decay length of radiation source
c     
      elseif (tag(1:3).EQ.'L41') THEN
        CALL ReadR(line,csollr,0.0,1.0,
     >               'Radiation source length or lambda')
c
c     L42 - csolfr - fractional strength of radiation source (cpopt 2,3)
c      
      elseif (tag(1:3).EQ.'L42') THEN
        CALL ReadR(line,csolfr,0.0,HI,
     >               'Radiation source fraction')
c
c     L43 - csolpr - absolute strength of radiation source (cpopt 0,1)      
c
      elseif (tag(1:3).EQ.'L43') THEN
        CALL ReadR(line,csolpr,0.0,HI,
     >               'Radiation source strength (absolute)')
c
c     L44 - cfiz - fractional split between linear and exponential ionization sources
c                  in csopt 4,5      
c
      elseif (tag(1:3).EQ.'L44') THEN
        CALL ReadR(line,cfiz,0.0,1.0,
     >               'Ionization source fraction split opt 4,5')
c
c     L45 - sol13_padd - additional pressure loss to be applied
c
      elseif (tag(1:3).EQ.'L45') THEN
        CALL ReadR(line,sol13_padd,0.0,HI,
     >               'Additional pressure source')
c     
c     L46 - sol13_pdist - distance over which to add pressure - 0.0 = change pressure at target by factor
c
      elseif (tag(1:3).EQ.'L46') THEN
        CALL ReadR(line,sol13_pdist,0.0,1.0,
     >               'Additional pressure distribution distance')
c -----------------------------------------------------------------------
c
c     L47 - Te profile radial shift
c     L48 - Ti profile radial shift
c     L49 - ne profile radial shift
c
      elseif (tag(1:3).EQ.'L47') THEN
        CALL ReadR(line,te_prof_shift,-HI,HI,'Te profile shift')
      elseif (tag(1:3).EQ.'L48') THEN
        CALL ReadR(line,ti_prof_shift,-HI,HI,'Ti profile shift')
      elseif (tag(1:3).EQ.'L49') THEN
        CALL ReadR(line,ne_prof_shift,-HI,HI,'ne profile shift')
c     
c
c     L50 - Te profile multiplier      
c     L51 - Ti profile multiplier      
c     L52 - ne profile multiplier      
c
      elseif (tag(1:3).EQ.'L50') THEN
        CALL ReadR(line,te_prof_mult,-HI,HI,'Te profile mult')
      elseif (tag(1:3).EQ.'L51') THEN
        CALL ReadR(line,ti_prof_mult,-HI,HI,'Ti profile mult')
      elseif (tag(1:3).EQ.'L52') THEN
        CALL ReadR(line,ne_prof_mult,-HI,HI,'ne profile mult')
c
c
c        
c -----------------------------------------------------------------------
c
c     TAG Q26:
c
c     Specification of a density multiplier (gradient) to be applied
c     to the outboard region. 
c
      elseif (tag(1:3).EQ.'Q26') THEN
c
c
C     READ IN DENSITY GRADIENT INFORMATION, IF ANY
C     FORM IS POSITION (AS PORTION OF L) AND VALUE AS A MULTIPLIER
C     OF THE DENSITY
C

         CALL RDRARN(MNBG,NNBG,MAXINS,-MACHLO,MACHHI,.TRUE.,0.0,MACHHI,            
     >                                     1,'SET OF Y,MNB VALUES',IERR)
C
C
C
C -----------------------------------------------------------------------
C
C     TAG M01:
C
C     Set initial velocity angle of neutrals as 0.0 unless otherwise
C     specified when CNEUTC=17
C
      ELSEIF (TAG(1:3).EQ.'M01') THEN
        CALL ReadR(line,CIANGN,-180.0,180.0,'Initial neutral velocity
     > angle')       
        ciangn = ciangn * degrad  
c
c
c -----------------------------------------------------------------------
c
c     ADD TAGS RELATED TO OUT - USING SERIES 'O' oooh :) ... for OUT
c
c -----------------------------------------------------------------------
c
c     TAG O01: 
c
c     Specify an alternate absolute factor in LIM
c
      elseif (tag(1:3).eq.'O01') then 
c
c     This option allows an absolute scaling factor for the LIM
c     run results to be specified in the OUT routine. It's default
c     value is zero.
c
c        CALL ReadR(line,new_absfac,0.0,HI,
c     >                   'Imposed ABSFAC in OUT')
        CALL ReadDP(line,new_absfac,0.0,HI,
     >                   'Imposed ABSFAC in OUT')
c
c -----------------------------------------------------------------------
c
c     TAG O99:
c
c     Net erosion plot scaling option
c     0 = normal (which is particles/m /particle entering the system)
c     1 = mm/hr   ABSFAC * 3600 / 1.22e26 for Beryllium!!
c     2 = Not yet implemented
c
      elseif (tag(1:3).eq.'O99') then 
c      
        CALL ReadI(line,erosion_scaling_opt,0,2,
     >                'Erosion Scaling Option')
c         
c
c -----------------------------------------------------------------------
c
c     TAG does not match available input - signal ERROR and EXIT
c

      ELSE
        CALL ER('ReadUnstructuredInput','Unrecognized tag',*99)
      endif
c

      RETURN
c
c There is an error:
c
99    WRITE(6,'(5X,3A)') 'LINE = "',line,'"'
      WRITE(6,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0    ,'(5X,3A)') 'LINE = "',line(1:LEN_TRIM(line)),'"'
      WRITE(0    ,'(5X,3A)') 'TAG  = "',tag ,'"'
      WRITE(0,*) '    DIVIMP HALTED'
      STOP
      END


c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadIR(line,ival,rval,imin,imax,tag)
      use mod_io_units
      !use mod_params
      use mod_slcom
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax
      REAL    rval

      INTEGER i
      REAL    r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i,r

      IF (i.LT.imin.OR.i.GT.imax)
     .  CALL ER('ReadI','Out of bounds: '//line,*99)

      ival = i
      rval = r

      WRITE(STDDBG,'(A)') line
      WRITE(STDDBG,'(5X,2A,I4,1P,E10.2)') tag,' = ',ival,rval

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadI(line,ival,imin,imax,tag)
      use mod_params
      use mod_slcom

      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i

      IF (i.LT.imin.OR.i.GT.imax) then 

        write (0,*)  'READI:ERROR:',i,imin,imax 
        CALL ER('ReadI','Out of bounds: '//line,*99)

      endif

      ival = i

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3I4)') 'I,IVAL,IMIN,IMAX = ',i,ival,imin,imax
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadC(line,cval,tag)
      use mod_params
      use mod_slcom

      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      CHARACTER line*(*),tag*(*),cval*(*)
      INTEGER fp,ival,imin,imax

      INTEGER i
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,cval

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,A)') tag,' = ',cval

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,2A)')    'CVAL = ''',cval,''''
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE Read2I(line,ival1,ival2,imin,imax,tag)
      use mod_params
      use mod_slcom

      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      CHARACTER line*72,tag*(*)
      INTEGER fp,ival1,ival2,imin,imax

      INTEGER i1,i2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,i1,i2

      IF (i1.LT.imin.OR.i1.GT.imax.OR.
     .    i2.LT.imin.OR.i2.GT.imax)
     .  CALL ER('Read2I','Out of bounds: '//line,*99)

      ival1 = i1
      ival2 = i2

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival1
      WRITE(STDDBG,'(5X,2A,I4)') tag,' = ',ival2

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      STOP
      END
c
c
c ======================================================================
c
c
c
      SUBROUTINE ReadR(line,rval,rmin,rmax,tag)
      use error_handling
      use mod_params
      use mod_slcom
      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval,rmin,rmax

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      REAL r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r

      IF (r.LT.rmin.OR.r.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval = r

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(2A,G10.3)') tag,' = ',rval

      RETURN

 98   call errmsg('READR','Problem reading unstructured input')
      WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r,rval,rmin,rmax
      STOP
      END
c
c ======================================================================
c
c
c
      SUBROUTINE ReadDP(line,dpval,rmin,rmax,tag)
      use error_handling
      use mod_params
      use mod_slcom
      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rmin,rmax
      real*8 dpval

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      REAL*8 r
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r

      IF (r.LT.rmin.OR.r.GT.rmax)
     .  CALL ER('ReadDP','Out of bounds: '//line,*99)

      dpval = r

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(2A,G10.3)') tag,' = ',dpval

      RETURN

 98   call errmsg('READDP','Problem reading unstructured input')
      WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,3G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r,dpval,rmin,rmax
      STOP
      END

c
c ======================================================================
c
c
c
      SUBROUTINE Read2R(line,rval1,rval2,rmin,rmax,tag)
      use mod_params
      use mod_slcom

      IMPLICIT none

      CHARACTER line*72,tag*(*)
      REAL rval1,rval2,rmin,rmax

c      INCLUDE 'params'
c      INCLUDE 'slcom'

      REAL r1,r2
      CHARACTER comment*72

      READ (line,*,ERR=98,END=98) comment,r1,r2

      IF (r1.LT.rmin.OR.r1.GT.rmax.OR.
     .    r2.LT.rmin.OR.r2.GT.rmax)
     .  CALL ER('ReadR','Out of bounds: '//line,*99)

      rval1 = r1
      rval2 = r2

      WRITE(STDDBG,'(A)')        line
      WRITE(STDDBG,'(2A,2G10.3)') tag,' = ',rval1,rval2

      RETURN
98    WRITE(DATUNIT,*) 'Problem reading unstructured input'
99    WRITE(DATUNIT,'(5X,2A)')    'LINE = ''',line,''''
      WRITE(DATUNIT,'(5X,2A)')    'TAG  = ''',tag,''''
      WRITE(DATUNIT,'(5X,A,6G10.3)')
     .  'R,RVAL,RMIN,RMAX = ',r1,r2,rval1,rval2,rmin,rmax
      STOP
      END







