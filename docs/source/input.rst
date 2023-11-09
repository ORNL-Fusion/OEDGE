Input Options
=============

OSM and DIVIMP (collectively called OEDGE here) share the same input file as they are two tightly integrated codes that share many of the same options. This inevitably means many options will only apply to OSM and not to DIVIMP, and vice-versa. OEDGE also has a long and rich history of development, thus many input options are considered defunct. 

The input options are separated into "Tags" consisting of a letter and a number, e.g., T13. There are also one set of Tags that start with the number 2 that relate to the SOL22 plasma solver within OSM. Generally, the Tags are supposed to represent a grouping of related input options, but this convention is not always followed. 

============ ============================
  `A Tags`_   
-----------------------------------------
  `A01`_      Title
  `A02`_      Case Description
  `A03`_      Equilibrium Grid File Name
  `A04`_      Print Option
  A05         To be documented. 
  A06         To be documented.
  A07         To be documented.
============ ============================

============ ==========================================
  `B Tags`_   
-------------------------------------------------------
  B01         Obsolete
  B02         Obsolete
  B03         Obsolete
  `B04`_      Debug NEUT
  `B05`_      Debug DIV
  `B06`_      Debug Ion Velocity
  `B07`_      Z-Limit for Divertor
  `B08`_      Ring Number for Detailed Background Data
============ ==========================================

============ ========================================================
  `C Tags`_   
---------------------------------------------------------------------
  `C01`_      Set of S-Values for Ion Leakage Diagnostic
  The Dperp/Chiperp Extractor
---------------------------------------------------------------------
  `C02`_      Extractor Methods
  `C03`_      Extractor Range
  `C04`_      Include Outer Ring Losses
  `C05`_      Dperp Convection
  `C06`_      1/2 Cell Flux Correction
  `C07`_      Calculate Average Coefficients
  `C08`_      Major Radius Correction
  `C09`_      Gradient Smoothing
  `C10`_      Gradient Calculation Method
  `C11`_      Cross-Field Area Option
  `C12`_      Power Loss Terms
  `C13`_      Non-Orthogonal Correction
  `C14`_      Pei Correction Factor
  `C15`_      Recycling Coefficient Correction
  `C16`_      Extractor Dperp/Xperp Ratio Specification
End Dperp/Chiperp Extractor
---------------------------------------------------------------------
  `C17`_      Vertical Reciprocating Probe - R Crossing Number
  `C18`_      Vertical Reciprocating Probe - R Location
  `C19`_      Horizontal Reciprocating Probe - Z Crossing Location
  `C20`_      Horizontal Reciprocating Probe - Z Location
  C21         To be documented.
  C22         To be documented.
============ ========================================================

============ =========================================================================================
  `D Tags`_   
------------------------------------------------------------------------------------------------------
  `D01`_      Ionization Data Source Option
  `D02`_      Source Data Option Specifications – Userid for H Database
  `D03`_      Source Data Option Specifications – H Data Year
  `D04`_      Source Data Option Specifications – Userid for Z (Impurity) Database
  `D05`_      Source Data Option Specifications – Z (Impurity) Data Year
  `D06`_      Atomic Data Source File Name (Adpak, Strahl …)
  `D07`_      Physical Sputter Data – Source Option
  `D08`_      Chemical Sputter Yield Option
  `D09`_      Momentum Transfer Collision – First Coefficient (Kelighi)
  `D10`_      Momentum Transfer Collision – Second Coefficient (Kelighg)
  `D11`_      Characteristic Energy Ebd (eV)
  `D12`_      Neutral Hydrogen Density Parameter – Nhc (m\ :sup:`-3`)
  `D13`_      Nho (m\ :sup:`-3`)
  `D14`_      Lamhx (m)
  `D15`_      Lamhy (m)
  `D16`_      Constant for CX Recomb Option 2 – Vcx (m/s)
  `D17`_      Threshold Yield for Self-Sputtering (eV)
  `D18`_      Bombarding Ion Charge State – Zimp
  `D19`_      Bombion (0-mb 1-H 2-D 3-T 4-He4 5-C 6-mi 7-O)
  `D20`_      Ionization Rate Factor for Neutrals – Irf
  `D21`_      Sputtering Enhancement Factor – Sef
  `D22`_      Set of Yield Modifiers for Primary, Secondary and Chemically Sputtered Neutrals
  `D23`_      Specified Fixed Yield Value for Sputter Data Option 4
  `D24`_      Target Temperature in K
  `D25`_      Wall Temperature in K
  `D26`_      Private Plasma Wall Temperature
  `D27`_      Specific Wall Segment Temperatures
  `D28`_      Temperature Gradient Coefficient Parameter – Z0
  `D29`_      Emax-Factor for Neut Launch Velocity – Emaxf
  `D30`_      Impurity Ion Impact Energy for Wall Launch Neutrals
  `D31`_      Maximum Number of Sputtered Generations
  `D32`_      Absfac or Power to Targets
  `D33`_      Stgrad: Gradient End-Point Specifier
  `D34`_      H Recombination Calculation Option
  `D35`_      H Recombination - Limiting Cutoff Temperature
  `D36`_      T-Grad Modification Factor
  D37         To be documented.
  D38         To be documented.
  D39         To be documented.
  D40         To be documented.
============ =========================================================================================

============ ========================================================================
  `F Tags`_   
-------------------------------------------------------------------------------------
  `F01`_       Read Fluid Code Background for Reference
  `F02`_       Fluid Code Target Data Usage Option
  `F03`_       Lost Sol Ring Option
  `F04`_       Velocity Multiplication Factors – For Data Read from Files
  `F05`_       Sonnet Grid: Number of Fluid Results in Background Plasma File (Nfla)
  `F06`_       Read Background Plasma Auxiliary Input File
  F07          To be documented.
  F08          To be documented.
  F09          To be documented.
  F10          To be documented.
  F11          To be documented.
  F12          To be documented.
  F13          To be documented.
  F14          To be documented.
  F15          To be documented.
  F16          To be documented.
  F17          To be documented.
  F18          To be documented.
  F19          To be documented.
  F20          To be documented.
============ ========================================================================

============ ===========================================================================================
  `G Tags`_   
--------------------------------------------------------------------------------------------------------
  `G01`_       Grid Option
  `G02`_       Non-orthogonal Option
  `G03`_       Parallel Distance Option
  `G04`_       Cross-Field Distance Option
  `G05`_       R, Z Calculation Option
  `G06`_       XY Grid Option
  `G07`_       Cell Area Calculation Option
  `G08`_       Ion Wall Option
  `G09`_       Neutral Wall Option
  `G10`_       Trap Wall Option
  `G11`_       Vessel Wall Redefinition Option (Baffle Inclusion)
  `G12`_       Target Position Option
  `G13`_       Pre-defined Geometry Selection Option
  `G14`_       Ring Location of Core Mirror – Ircore
  `G15`_       Rectangular Grid for Neutrals
  `G16`_       Set of Target Coordinates
  `G17`_       Set of Wall Coordinates
  `G18`_       Set of Trap Wall or ITER Second Wall Coordinates
  `G19`_       Sonnet Grid Characteristic Specifications – ASDEX U – Cmod – TEXTOR: Number of Rings
  `G20`_       Sonnet Grid: Number of Knots
  `G21`_       Sonnet Grid: Cut Ring Number
  `G22`_       Sonnet Grid: Cut Point 1
  `G23`_       Sonnet Grid: Cut Point 2
  G24-G56      To be documented.
============ ===========================================================================================

============ ===========================================================================================
  `H Tags`_   
--------------------------------------------------------------------------------------------------------
  `H01`_       PIN/Nimbus Random Number Seed 
  `H02`_       PIN/Nimbus Print Option
  `H03`_       Run PIN Option
  `H04`_       PIN Command Line
  `H05`_       PIN Cell Area Option (Ihcorr)
  `H06`_       PIN Hybrid Wall Option
  `H07`_       PIN Puffing Option
  `H08`_       PIN Puff Location Switch
  `H09`_       PIN Puff Fraction – Hpcpuf
  `H10`_       PIN Flux Puff Fraction – Ppcpuf
  `H11`_       PIN Puff Injection Temperature (Ev) – Tpufh
  `H12`_       PIN Puff Location Indices – JHPUF1
  `H13`_       PIN Puff Location Indices – JHPUF2
  `H14`_       Nimbus Namelist Input: Nimbin
Hydrocarbon Module Options
--------------------------------------------------------------------------------------------------------
  H15-H64      To be documented.
  H90-H91      To be documented.
============ ===========================================================================================

============ ===========================================================================================
  `I Tags`_   
--------------------------------------------------------------------------------------------------------
  `I01`_       Injection
  `I02`_       First Diffusion
  `I03`_       Control Switch
  `I04`_       Self Sputtering Option
  `I05`_       Initial Ion Velocity Option
  `I06`_       Follow Recombined Impurity Neutral Option
  `I07`_       Prompt Deposition Option
  `I08`_       Target Mirror Option
  `I09`_       Ion Periphery Option
  `I10`_       Periphery Recycle Option
  `I11`_       Z Effective (Self) – Zeff
  `I12`_       Initial Ionization State of Impurity Ions
  `I13`_       Collision Enhancement Factor – Zenh
  `I14`_       Set Ti = Max(ti, TB) When Reaching State
  `I15`_       Maximum Ionization State
  `I16`_       Stop Following Ions Reaching Main Plasma
  `I17`_       Ion Loss Time
  `I18`_       Ring Number for Ion Injection – Injection Option 2,3, 5,6 – Injir
  `I19`_       Injection Region -Lower Bound-Injection Option 2,3, 5,6 – INJ1
  `I20`_       Injection Region -Upper Bound-Injection Option 2,3, 5,6 – INJ2
  `I21`_       Far Periphery Width Definition
  `I22`_       Far Periphery Target Loss – Characteristic Time
  `I23`_       Far Periphery Diffusion Rate
  I24          To be documented.
  I25          To be documented.
  I26          To be documented.
  I27          To be documented.
  I28          To be documented.
  I29          To be documented.
  I30          To be documented.
  I31          To be documented.
  I32          To be documented.
  I33          To be documented.
  I34          To be documented.
  I35          To be documented.
  I36          To be documented.
  I37          To be documented.
  I38          To be documented.
============ ===========================================================================================

============ ===========================================================================================
  `K Tags`_  (ERO Interface?) 
--------------------------------------------------------------------------------------------------------
  K01-K??      To be documented.

============ ===========================================================================================

============ ===========================================================================================
  `N Tags`_   
--------------------------------------------------------------------------------------------------------
  `N01`_       Launch
  `N02`_       Vel/Angle Flag
  `N03`_       Supplementary Launch Option
  `N04`_       Supplementary Velocity/Angle Flag
  `N05`_       Initial Neutral Velocity/Angle Flag
  `N06`_       Extra 2D Neutral Launch Option
  `N07`_       2D Neutral Launch – Velocity/Angle Flag Option
  `N08`_       Sputter Option
  `N09`_       Secondary Sputter Option (TN1209)
  `N10`_       Normal
  `N11`_       Neut Spreading
  `N12`_       Impurity Neutral Velocity Type Option
  `N13`_       Neutral Reflection Option
  `N14`_       Impurity Neutral Momentum Transfer Collision Option
  `N15`_       Measure Theta from T Degrees for Launch
  `N16`_       Wall Launch Segment Probability Multipliers
  `N17`_       Absolute Wall Probabilities
  `N18`_       Power of Cosine Release Distribution (V/A Flag 12,13)
  `N19`_       Velocity Multiplier for Velocity/Angle Flag 14 and 15
  `N20`_       Velocity Multiplier for Recombined Ions
  `N21`_       External Sputtering Flux Data Source
============ ===========================================================================================

============ ============================
  `O Tags`_  
-----------------------------------------
  O01-O??      To be documented.
============ ============================

================ ===========================================================================================
  `P Tags`_   
------------------------------------------------------------------------------------------------------------
  `P01`_         SOL
  `P02`_         Core Plasma Options
  `P03`_         Plasma Decay
  `P04`_         Piece-Wise Background Plasma Option Inputs
  `P05`_         Trap Tgrad Option
  `P06`_         SOL Enhancement Factor – Electric Field – Solef
  `P07`_         SOL Enhancement Factor – Drift Velocity – Solvf
  `P08`_-`P14`_       SOL Parameters – FL, Fs, Frm, Kin, Kout, Frmin, Frmax
  `P15`_-`P18`_       SOL Parameters for Sol Options 6 and 7.
  `P19`_         Power Density – P/A
  `P20`_         Parallel Heat Conduction Coefficient – K0
  `P21`_         Parallel Ion Heat Conduction Coefficient – K0I
  `P22`_         Electric Field Option – Overrides Other E-field Options or Data
  `P23`_         Electric Field Option 4 – Source Length Specifier
  `P24`_         Electric Field Option 4 – Collisional Determination Factor
  `P25`_         Ionization Source-Characteristic Length- SOL12 to 15 – Ls
  `P26`_         Ionization Source- Second Characteristic Length – L2
  `P27`_         Ionization Source- Source Fraction – Fi
  `P28`_         Radiation Source-Characteristic Length- SOL12 to 15 – Lr
  `P29`_         Radiative Power Constant – SOL12 to 15 – PR/a (W/m2)
  `P30`_         Radiation Source Strength Fraction – Frr
  `P31`_         Ionization Source Option – SOL12 to 15
  `P32`_         Radiative Source Option – SOL12 to 15
  `P33`_         Imaginary Root Option
  `P34`_         Flux Recirculation Option – SOL 12 to 15
  `P35`_         Flux Recirculation – Source Specifications
  `P36`_         Iterate SOL Option
  `P37`_         Secondary SOL Option
  `P38`_         Ionization Option for Iterative SOL
  `P39`_         Number of PIN/SOL Iterations
  `P40-P55`_     Private Plasma (Trap) Specification Option Inputs
  `P56-P59`_     Input Parameters for Core Option 4 and 5 (Marfe Simulation)
  P60-P66        To be documented.
================ ===========================================================================================

================ ===========================================================================================
  `Q Tags`_   
------------------------------------------------------------------------------------------------------------
  `Q01`_          Teb Gradient Option
  `Q02`_          Tib Gradient Option
  `Q03`_          Forced Flat Temperature Gradient Option
  `Q04`_          Te Gradient Cut-Off for Flattening Option
  `Q05`_          Ti Gradient Cut-Off for Flattening Option
  `Q06-Q11`_      Electron Temperatures – TeB0, Tebp, Tebout, Tebin, Tebt, 
  `Q12-Q15`_      Electron Gradient Parameters – FEBL1, FEBL2, Febt, FEB2, 
  `Q16-Q21`_      Ion Temperatures - TiB0, Tibp, Tibout, Tibin, Tibt
  `Q22-Q25`_      Ion Gradients - FIBL1, FIBL2, Fibt, FIB2
  `Q26-Q31`_      Densities – NB0, Nebp, Nbout, Nbin, Nbt, Nboup
  `Q32`_          Langmuir Probe Data Switch
  `Q33`_          Inner/Both Target Data Multipliers
  `Q34`_          Langmuir Probe Data Input – Inner/Both Plate
  `Q35`_          Outer Target Data Multipliers
  `Q36`_          Langmuir Probe Data Input – Outer Plate
  `Q37`_          Core Plasma Input Data
  `Q38`_          Obsolete
  `Q39`_          Obsolete
  `Q40`_          Outboard Plasma Flow Vel (Sol 5,6,7) – Vhyout(m/s)
  `Q41`_          Outboard Electric Field (Sol 5,6, 7) – Eyout (V/m)
  Q42             To be documented.
  Q43             To be documented.
  Q44             To be documented.
  Q45             To be documented.
================ ===========================================================================================

============ ===========================================================================================
  `R Tags`_   
--------------------------------------------------------------------------------------------------------
  R01-R15      To be documented.

============ ===========================================================================================

================ ===========================================================================================
  `S Tags`_   
------------------------------------------------------------------------------------------------------------
  `S01`_          On Axis B-field Value
  `S02`_          Mass of Plasma Ions – Mb
  `S03`_          Charge on Plasma Ions – Zb
  `S04`_          Mass of Impurity Ions – Mi
  `S05`_          Atomic Number of Impurity Ions – Zi
  `S06`_          Initial Temperature – TEM1 (Ev)
  `S07`_          Initial Temperature (2) – TEM2 (Ev)
  `S08, S09`_     Initial R, Z Position of Impurity     
  `S10`_          DIVIMP Mode (1 Impulse, 2 Steady State, 0 Both)
  `S11`_          Number of Impurity Ions to Be Followed
  `S12`_          Number of Supplementary Particles to Be Followed
  `S13`_          Quantum Iteration Time in Neut – Fsrate (S)
  `S14`_          Quantum Iteration Time in Div – Qtim (S)
  `S15`_          CPU Time Limit (S)
  `S16`_          Average Dwell Times (S) for Each Charge State
  `S17`_          Dwell Time Factors for Time Dependent Analysis
  `S18`_          Maximum Dwell Time for Steady State
  `S19`_          Random Number Seed (0 Generate New Seed)
  `S20`_          Number of DIVIMP Iterations
  `S21`_          SOL Test Option
  S22             To be documented.
  S23             To be documented.
  S24             To be documented.
================ ===========================================================================================

============ ===========================================================================================
  `T Tags`_   
--------------------------------------------------------------------------------------------------------
  `T01`_       Ionization
  `T02`_       Collision
  `T03`_       Reiser Coulomb Collision Transport Option
  `T04`_       Friction
  `T05`_       Heating
  `T06`_       Cx Recomb
  `T07`_       Dperp Option
  `T08`_       Perpendicular Step Option
  `T09`_       Pinch Velocity Option
  `T10`_       Teb Grad Coeff Option
  `T11`_       Tib Grad Coeff Option
  `T12`_       Temperature Gradient Force – Modification Option
  `T13`_       Poloidal Drift Option
  `T14`_       Cross-Field Diffusion Rate – Dperp (m\ :sup:`2` /s)
  `T15`_       Cross-Field Diffusion Rate for Private Plasma Region – Dperpt (m\ :sup:`2` /s)
  `T16`_       Perpendicular Pinch Velocity – Cvpinch (m/s)
  `T17`_       Poloidal Drift Velocity – Vpol (m/s)
  `T18`_       Poloidal Drift Velocity – Range of Effect
  `T19`_       To be documented.
  `T20`_       To be documented.
  `T21`_       To be documented.
  `T22`_       To be documented.
  `T23`_       To be documented.
  `T24`_       To be documented.
  `T25`_       To be documented.
  `T26`_       To be documented.
  `T27`_       To be documented.
  `T28`_       To be documented.
  `T29`_       To be documented.
  `T30`_       To be documented.
  `T31`_       To be documented.
  `T32`_       To be documented.
  `T33`_       To be documented.
  `T34`_       To be documented.
  `T35`_       To be documented.
  `T36`_       To be documented.
  `T37`_       To be documented.
  `T38`_       To be documented.
  `T39`_       To be documented.
  `T40`_       To be documented.
  `T41`_       To be documented.
  `T42`_       To be documented.
  `T43`_       To be documented.
  `T44`_       To be documented.
  `T45`_       To be documented.
  `T46`_       To be documented.
  `T47`_       To be documented.
  `T48`_       To be documented.
  `T49`_       To be documented.
  `T50`_       To be documented.
  `T51`_       To be documented.
  `T52`_       To be documented.
  `T53`_       To be documented.
  `T54`_       To be documented.
  `T55`_       To be documented.
  `T56`_       To be documented.
  `T57`_       To be documented.
  `T58`_       To be documented.
  `T59`_       To be documented.
  `T60`_       To be documented.
  `T61`_       To be documented.
  `T62`_       To be documented.
============ ===========================================================================================

============ ============================
  `W Tags`_  
-----------------------------------------
  W01           To be documented.
  W02           To be documented.
============ ============================

============ ============================
  `Z Tags`_  (SOL 29) 
-----------------------------------------
  Z01-Z??      To be documented.
============ ============================

============ ========================================================================
  `200 Tags`_  (SOL 22)
-------------------------------------------------------------------------------------
`201`_         Force Te = Ti
`202`_         Initially Imposed Target Mach Number
`203`_         Initial Mach Number Step Size
`204`_         Ultimate Mach Number Resolution
`205`_         Ionization Source Length Switch
`206`_         Start of Ionization Source
`207`_         End or Length of Ionization Source
`208`_         Decay Factor or Width of Ionization Source
`209`_         Length of the Radiation Source
`210`_         Decay Length of Radiation Source
`211`_         Source Strength Fraction (Frr)
`212`_         Garching Radiation Model: Parameter Alpha
`213`_         Garching Radiation Model: Temperature Base
`214`_         Garching Radiation Model: First Exponent
`215`_         Garching Radiation Model: Second Exponent
`216`_         Correction Factor to Gamma (Ion)
`217`_         Correction Factor to Gamma (Electron)
`218`_         CX Power Coefficient
`219`_         Recycling Source Coefficient
`220`_         Pei (Equipartition) Correction Factor
`221`_         Velocity Error Switch
`222`_         Distributed Power Start Position
`223`_         Distributed Power End Position
`224`_         Compound Gperp – Fraction of Gperp in Rectangular Distribution
`225`_         Compound Gperp – Start of Region
`226`_         Compound Gperp – End of Region
`227`_         Extra Perpendicular Source/Sink Strength
`228`_         Range of S-values for Extra Source
`229`_         Range of S-values for Extra Sink
`230`_         Distance Factor for PP Power Loss Re-distribution
`231`_         Start Knot Index for EDGE2D Compatibility Option
`232`_         Fill Option for Skipped Cells in EDGE2D Compatibility Option 9
`233`_         Qe Term – Temperature Cutoff (eV)
`234`_         PINQID – Atomic Ionization – T Cutoff (eV)
`235`_         PINQID – Molecular Ionization – T Cutoff (eV)
`236`_         PINQID – Recombination – T Cutoff (eV)
`237`_         Qi Term / PINQID – Charge Exchange – T Cutoff (eV)
`238`_         PINQID – Charge Exchange Option 1 – Reference Temperature (eV)
`239`_         Minimum Temperature Allowed in Solver (SOL 22)
`240`_         Maximum Allowed Temperature Drop Fraction
`241`_         Momentum Loss Term Multiplier
`242`_         Friction Factor for Momentum Loss
`243`_         Length of Momentum Loss Region
`244`_         Decay Length of Momentum Loss
`245`_         Correction Ratio of CX to Ionization Events for Momentum Transfer
`246`_         Te Cut-Off for Increased CX Multiplier
`247`_         Te Lower Limit Cutoff for CX Multiplier
`248`_         PINQE (Electron Energy Loss) Term Multiplier
`249`_         Prad Option 3 – (Multiplier for PINQE)
`250`_         Initial Number of Runge-Kutta Steps Between Grid Points
SOL 22 Switches
-------------------------------------------------------------------------------------
`251`_         Ionization Option
`252`_         Initial Ionization Option:
`253`_         Private Plasma Ionization Option
`254`_         5/2 nv kT Term
`255`_         1/2 mv\ :sup:`3` n Term
`256`_         Prad Option
`257`_         Phelpi Option
`258`_         Pei Option
`259`_         Pcx Option
`260`_         PINQID – DIVIMP Calculated Qi – Atomic Ionization
`261`_         PINQID – DIVIMP Calculated Qi – Molecular Ionization
`262`_         PINQID – DIVIMP Calculated Qi – Recombination
`263`_         PINQID – DIVIMP Calculated Qi – Charge Exchange
`264`_         PP Target Electron Power Loss Redistribution Option
`265`_         PP Target Ion Power Loss Redistribution Option
`266`_         Viscosity Option
`267`_         Momentum Loss Option
`268`_         Iterative Mach Number Option
`269`_         Edge 2D Data Compatibility Option
`270`_         Power Distribution Option
`271`_         Private Plasma Power Distribution
`272`_         Gamma Perp Option
`273`_         Private Plasma Gamma Perp Option
`274`_         Extra Gperp Source/Sink Option
`275`_         Major Radius Option
`276`_         Core Flux Source
`277`_         Recombination Source Option
`278`_         Smoothing Option
`279`_         Detached Plasma Prescription Option
`280`_         Error Correction Level
`281`_         Automatic Default Error Correction

============ ========================================================================


A Tags
------

.. _A01:
A01 : Title
  This entry is composed of two string constants. The first is fixed and simply indicates that the first line of the input file is reserved for the title of the run. The second string, after leaving a space, (See the sample input file), contains a title for the case. Typical contents would be the name and series number and perhaps a reference to a note or other document that would describe why the case was run. These are useful in the future when cross-checking the cases that have been run. The title will appear in all print outs and on plots produced from the output of the case.

.. _A02:
A02 : Case Description
  This input is another DIVIMP character string value. The first string is an identifier tag and the second string should contain a description of the important features or options for the case. This description will be included in the header information of the .dat/.html case file and will be included in the description section of the posted case database. 

.. _A03:
A03 : Equilibrium Grid File Name
  This input line specifies the complete path-inclusive name for the equilibrium file that is being used to run the case. This information is then passed to PIN/NIMBUS, if it is being run, for its use in loading the equilibrium grid. At this time, the information is not directly used by DIVIMP itself. DIVIMP relies on the script file to pre-connect the equilibrium grid to the appropriate input unit number.

  e.g. 'Equil File Name' '/u/progs/div4/shots/g31627.v3'

.. _A04:
A04 : Print option (0 reduced, 1+)
  This option provides a means of selecting the types of printouts. The reduced printout (enter 0 here) is enough for most cases. The various other printouts give additional information on various aspects of the simulation. Option 1 includes extra print-outs on bin sizes, ionization rates, characteristic times, and other items and is occasionally required. Other values that are used in the code are 2, 3, 4, and 9. These supply different print-outs and may be used for debugging. Print option 10 requests the code to write out the calculated DIVIMP plasma background in a DIVIMP specific format. This plasma background may be read in by using plasma decay option 98.

  Print Option 0 : Standard DIVIMP print out. Adequate for most cases.

  Print Option 1 : Standard print out plus the following:
    - Dperp Extractor Print Out
    - Fast Scanning Probe Data
    - Private Plasma Impurity Content Data

  Print Option 2 : Standard Print out plus the following:
    - Detailed core leakage and source description information.

  Print Option 3 : Standard Print out plus the following:
    - Debug - Additional Geometric Data about grids, targets and walls

  Print Option 4 : Standard Print out plus the following:
    - Additional debug information about EDGE2D target conditions

  Print Option 5 : Standard Print out plus the following:
    - Additional information about the background plasma conditions
    - Some extra characteristic times data
    - Retention predictor values

  Print Option 6 : Standard Print out plus the following:
    - Writes the grid information to a separate file in a SONNET style format.

  Print Option 7 : Standard Print out plus the following:
    - Extra PIN related data for debugging.

  Print Option 9 : Standard Print option plus all other possible print-outs. This option is the same as turning on ALL print options from 1 to 8. It is not recommended for use unless necessary since it generates a great deal of output.

  Print Option 10: Standard Print out plus the following:
    - Writes the finalized background plasma to a DIVIMP specific format that can be read in using Plasma Decay option 98. 

B Tags
------
.. _B04:
B04 : Debug NEUT
  Generally set this value (and the next) to 0. Setting this value to 1 generates some debugging information in the neutrals part of the program. It generates a complete history of each neutral followed, consisting of one line of output indicating the launch parameters for the neutral, followed by one line of output after each neutral timestep, and ending with a line indicating the fate of the neutral (e.g. ionised, hit wall, etc.). Hence a great deal of output is generated if the timestep is small! Setting this value to (say) 100 generates the first and last lines plus an extra line after every 100th timestep, similarly values of 1000, 10000 etc. could be given here.

.. _B05:
B05 : Debug DIV
  Generates histories of ions tracked by DIV. Again, a value of 0 switches the option off, a value of 1 produces copious output, and values of 100, 1000, 10000 etc. are often more helpful. Note: when variable timesteps are in use, the printed lines of debug are unlikely to occur at exactly the specified intervals. They might occur at 101.7, 215.3, 306.4 instead of 100, 200, 300 for example. A variable timestep feature is not currently implemented in DIVIMP. The purpose of such a feature is to transport the ions more quickly, in terms of CPU time (i.e. ion iterations) in regions where events are occurring very slowly relative to the base timestep, as may be the case with the core plasma rings.

.. _B06:
B06 : Debug Ion Velocity
  This option will generate debugging information about the ion velocity distribution. The variables used for this are declared in the common block DIAGVEL. Among other things, this option will generate a distribution of the particle velocities as a multiple of the impurity ion sound speed at the inner target on the separatrix ring. This can be used for debugging collision options based on parallel velocity diffusion. This option will also generate average diffusive step size information for the spatial diffusion options. The value used to turn the option on is not significant at this time. Any value greater than zero will suffice.

.. _B07:
B07 : Z-Limit for divertor
  Defines an arbitrary Z-point above which is considered the divertor influenced region. This is useful in calculating some quantities that are later plotted. One example is the density decay in the outermost SOL rings which is affected by the behaviour near the plates. In essence this quantity is used to define (very roughly) a near divertor region.

.. _B08:
B08 : Ring number for Detailed Background Data
  Some of the higher numbered SOL options may allow calculations of the background on a much finer grid than the bin system that is used for particle accounting. In order to save storage, these high resolution background values are not stored for every ring. This parameter specifies the one ring of interest for which these high resolution background data will be stored. This is useful for debugging purposes and for looking at transitional effects which are not visible on the scale of the larger bin system. (Not available for all SOL options) 

C Tags
------
.. _C01:
C01 : Set of S-values for Ion Leakage Diagnostic
  The first line specifies the number of entries. This is limited by the value of MAXPTS, set in the PARAMS common block. The rest of the data, one number per line, specifies the S-bin values at which the leakage information will be collected. When an ion exceeds the value of S along the field line that is listed in this table then the count for that distance will be incremented by the weight of the particle. The particle will not be counted more than once for any of the distances. The distances are specified in meters.

  e.g.

  .. code-block::
    
    ' ' 'Set of S-distances for ion leakage diagnostic (m)'
    'TN982 Number of S-values :-' 5
      1.0
      2.0
      3.5
      6.0
     13.35

  The limitations on these values are that they be greater than zero and be recorded in ascending order.

.. _C02:
C02 : Extractor Methods
  These were various methods used initially to calculate the transport coefficients. Only Option 2 should be used now and in later releases this input option will be removed.

  Option 0,1: OBSOLETE

  Option 2: Calculates transport coefficients ring by ring over the entire grid. Gradients are calculated for each cell on the grid, based on the values in adjacent cells.

.. _C03:
C03 : Extractor Range
  This option specifies the section of the field-line over which the extractor will work. Either the whole ring from target to target or only the sub-section of the ring from X-point region around to the X-point region on the other side. (The X-point to X-point region is approximately defined by the first cell on the separatrix with it's center located above the X-point.) Generally, it seems best to use the whole field line in the calculations.

  Option 0: Xpoint

  Option 1: Whole Ring

.. _C04:
C04 : Include Outer Ring Losses
  Examining the gradients in density and temperature, it becomes clear that there must be a cross-field ion and heat flux across the outermost boundary of the grid. In order to calculate the transport coefficients correctly, particularly in the rings closer to the outside, it may become necessary to include this term in the transport coefficient calculations. This may be particularly true for cases with large cross-field gradients on the outermost field lines. For cases where there is little or no cross-field variation in the background plasma at the outermost ring, there will also be little in the way of outer ring losses and so this term may not play a role in these cases. The actual amount of outer ring losses requires knowledge of the transport coefficients. This difficulty is overcome by assuming that the transport coefficient for the outer ring is the same as the ring currently being analyzed. This allows the gradient summation for the outer ring to be combined with the one for the current ring and thus allow a value for the transport coefficient to be extracted. This method is used when no data about the transport coefficient is available. However, when the Average Calculation Option (described below) is turned on, the values of the transport coefficients used for the outer ring losses are those calculated by averaging the coefficients found over a set of rings closer to the separatrix - themselves calculated using the first method described. This average is then applied to the outermost ring losses in the transport coefficient calculations for those rings beyond the averaging region.

  Option 0: Off - outermost ring loss corrections are left out of the transport coefficient calculations. This assumes that there are no particle or heat flows across the outermost plasma ring.

  Option 1: On - The methods described above are applied to calculate an outermost ring perpendicular loss that is factored into the transport coefficient calculations.

.. _C05:
C05 : Dperp Convection
  The cross-field flow of particles carries heat and this can be included in the calculation of the Xperp coefficients. The Dperp values for each ring are calculated first, this then allows the convective heat contribution from the Dperps to be added to the Xperp calculations. Unfortunately, under certain circumstances, it can be difficult to extract a reliable value of Dperp. For example, when the particle balance for the ring is almost entirely due to the ionization and target loss (i.e. when the actual amount of cross-field flow required is small), the calculation of Dperp can be problematic. When this occurs, it is probably better to leave this option off and avoid adding the noise in the Dperp extraction to the Xperp values.

  Option 0: Off

  Option 1: On

.. _C06:
C06 : 1/2 Cell Flux Correction
  This is also an obsolete option that was used in examining alternative methods of calculating the transport coefficients. It should always be turned OFF.

  Option 0: Off

.. _C07:
C07 : Calculate Average Coefficients
  This option will calculate the average values of the transport coefficients over the rings IRSEP+N to irsep+N+3 and then apply this averaged value to the outer ring loss calculations for all rings greater than irsep+N+3. The reason for this is that as one moves out the outer ring losses become a bigger contribution to the remaining flux that must be accounted for by the transport coefficients. As a result, if the transport coefficients for the outer rings are allowed to float the calculations can become very sensitive and unstable leading to large variations of the transport coefficients on the outer-most rings. This option can help to stabilize these variations.

  Option 0: Off

  Option N: On - calculate average as noted above where this is the value of N indicated.

.. _C08:
C08 : Major Radius Correction
  This option is aimed at correcting the extractor for major radius effects (toroidal geometry). After testing it was found that this option did not cause significant variation in the extracted values of the transport coefficients. The code may no longer be up to date with the latest extractor options. It is recommended that this option be left turned OFF.

  Option 0: Off

.. _C09:
C09 : Gradient Smoothing
  One concern in the extractor was noise in the gradients calculated from the plasma background grid. This option allows for averaging/smoothing of the gradients in an attempt to even out large variations. It is usually left off, at least for initial evaluation purposes. If turned on, a small value of N is recommended.

  Option 0: Off

  Option N: On - Calculate the gradient in cell (i,j) by averaging the values of the gradients in cells (i-N,j), (i-N+1,j) ... (i,j) ... (i+N,j) and assigning the average to the (i,j) cell. Where the i'th index runs across the field lines.

.. _C10:
C10 : Gradient Calculation Method:
  Different methods used to calculate the cross-field gradients at each cell of the grid. Option 0 is usually used.

  Option -1: On - the gradients are taken from a routine that fits a cubic spline to the set of data along each set of knots perpendicular to the field line. The cubic spline interpolation itself is not used, only the gradients returned by this routine. This is actually functionally equivalent to option 0 once the details of the cubic spline interpolation routine were examined.

  Option 0: On - gradient is calculated by taking the average of the gradient outward to the next grid cell and inward to the last grid cell.

Option 1: On - gradient is calculated by taking the value of the function at the inward and outward grid cells and the total distance between the inward and outward grid cells, ignoring the values in the current cell.

  Option 2: On - gradient is calculated by taking the value of the function at the current cell and the next outward cell and the distance between them.

.. _C11:
C11 : Cross-field Area Option
  The "area" across which the cross-field flux is moving is required in order to estimate the transport coefficients. This could be calculated either on the cell centres or at the edge of the respective cells. This is typically a very small difference and this option was only implemented to test if there was any appreciable change in the calculated transport coefficients when this was combined with other options. The cell centre Areas are the recommended ones and the ones most often used

  Option 0: On - Use cell centre evaluated poloidal lengths

  Option 1: On - Use cell boundary evaluated poloidal lengths

.. _C12:
C12 : Power Loss Terms
  This option will include the PIN calculated power loss terms in the calculation of the energy transport coefficients (Xperp).

  Option 0: Off

  Option 1: On - the values of PINQE and PINQI are read from the PIN output and used in the calculation of the extracted transport coefficients.

  Option 2: On - the values of PINQE and PINQI are read from the PIN output and used in the calculation of the extracted transport coefficients. The value of PEI the equipartition energy transport term is calculated analytically and also added in the calculation of the transport coefficients. The PEI term is modified by the Pei correction factor described below.

.. _C13:
C13 : Non-orthogonal Correction:
  This option is a first approximation at correcting the calculated gradients for non-orthogonal grids. Since the calculation should use the perpendicular gradients, it is necessary to correct the calculated gradients based on cell centre positions and quantities for the actual perpendicular distance involved. It is recommended to use option 1. However, since the majority of the grid is usually orthogonal this is often just a small correction.

  Option 0: Off - No Non-orthogonal corrections are performed.

  Option 1: On - the calculated value of the gradient for the cell is adjusted by dividing by the SIN(cell orthogonal angle). Where the cell orthogonal angle is 90 degrees or Pi/2.0 for an orthogonal cell and tends towards zero or Pi for a degenerate cell.

.. _C14:
C14 : Pei Correction Factor
  This value multiplies the analytic expression for the Pei energy transfer term used in the extractor. (In option 2 of the Power Loss Term option). It would normally be set to 1.0, however, allowing a specifiable parameter one can gauge the effect of equipartition on the output.

.. _C15:
C15 : Recycling Coefficient Correction
  This value is also usually equal to 1.0. However, in some cases it is known that the normalization of the PIN terms is not correct due to recycling, pumping or other effects. In these cases, it is not possible to extract appropriate transport coefficients because the strength of the ionization source is typically too large. This quantity is used to multiply all of the source terms that are extracted from the PIN data in order to account for recycling/pumping loss. It is also used to match Edge2D data where a recycling source fraction less than one is sometimes used to induce plasma flow from the core in the fluid simulation.

.. _C16:
C16 : Extractor Dperp/Xperp Ratio Specification
  Often the value of Dperp obtained by the extractor can be unreliable because it relies on the small difference between two large quantities. (The net remaining SOL plasma outflux and the total integrated ionization source over a section of the grid.) This unfortunately affects the value of Xperp extracted from the background. The Xperp value itself does not have the same sensitivity as does the value of Dperp. For this reason it can be useful to run the extractor using a fixed ratio of Dperp/Xperp. This will allow an estimate of the Xperp value to be extracted while losing the error caused by the noisy nature of the extracted Dperp value. The cost is the additional constraint of the system of equations. The ratio between Dperp and Xperp is not necessarily well known and as a result the asigned value may significantly affect the extracted Xperp values. 

.. _C17:
C17 : Vertical Reciprocating Probe - R crossing Number
  DIVIMP will produce output for a vertical fast scanning probe. It requires two values. The first is the intersection count. The code starts at the first target (IK=1, OUTER for JET, INNER for SONNET) and searches for intersections at the R value specified for the probe. This quantity specifies the R-intersection for which the data will be output. A value of 0 will turn off the vertical fast scanning probe output.

.. _C18:
C18 : Vertical Reciprocating Probe - R location
  This specifies the R-location for the vertical fast scanning probe. A value of -99.0 will turn off the probe output.

.. _C19:
C19 : Horizontal Reciprocating Probe - Z crossing number
  DIVIMP will produce output for a horizontal fast scanning probe. It requires two values. The first is the intersection count. The code starts at the first target (IK=1, OUTER for JET, INNER for SONNET) and searches for intersections at the Z value specified for the probe. This input specifies the Z-intersection for which the data will be output. A value of 0 will turn off the horizontal fast scanning probe output.

.. _C20:
C20 : Horizontal Reciprocating Probe - Z location
  This specifies the Z-location for the horizontal fast scanning probe. A value of -99.0 will turn off the probe output.

D Tags
------
.. _D01:
D01 : Ionization Data Source Option
  Source Data Option 0: Ionization and radiation data are taken from the NOCORONA subroutine package.

  Source Data Option 1: Ionization and radiation data are taken from the ADAS subroutine package.

  Source Data Option 2: B2-FRATRES formatted atomic physics data is used. The specific file name must be input using the MC-Filename option. This option is used for both ADPAK and STRAHL databases.

  Source Data Option 3: INEL formatted atomic physics database. The file name must be specified using the MC-Filename option described below.

.. _D02:
D02 : Source Data Option Specifications - UserID for H Database
  This is a character string specifying the path leading to the ADAS database files for the hydrogen data to be used in the case. If a '*' is specified this instructs the code to use the ADAS central database whose location has been defined by specifying it in the environment variable ADASCENT.

.. _D03:
D03 : Source Data Option Specifications - H data year
  The number entered as input is the year of the data from the specified hydrogen database to be used in the DIVIMP calculations. The ADAS database may have multiple sets of data for each element from differing years. An example would be, 93, which would select the 1993 from the specified hydrogen database.

.. _D04:
D04 : Source Data Option Specifications - UserID for Z (impurity) Database
  This is a character string specifying the path leading to the ADAS database files for the impurity data to be used in the case. If a '*' is specified this instructs the code to use the ADAS central database whose location has been defined by specifying it in the environment variable ADASCENT.

.. _D05:
D05 : Source Data Option Specifications - Z (impurity) data year
  The number entered as input is the year of the data from the specified impurity database to be used in the DIVIMP calculations. The ADAS database may have multiple sets of data for each element from differing years. An example would be, 89, which would select the 1989 from the specified impurity database.

.. _D06:
D06 : Name of file containing ADPAK/INEL atomic database
  MC-Filename Option. This line takes a character string entry that gives the complete path for the atomic ionization data that is to be used for Atomic Data Options 2 and 3 (INEL, ADPAK and STRAHL).

.. _D07:
D07 : Physical Sputter Data - Source Option
  Sputter Data Option 1: Data is taken from earlier publications.

  Sputter Data Option 2: Data is taken from Eckstein (1993)

  Sputter Data Option 3: Data is based on Eckstein (1993) with small changes to H,D,T coefficients taken from Garcia-Rosales and Roth(1996)

  Sputter Data Option 4: Specified CONSTANT yield value.

  Sputter Data Option 5: Data is based on Eckstein (1993) with small changes to H,D,T coefficients taken from Garcia-Rosales and Roth(1996). A customized routine has been created for the W sputtering values and is used instead of the Eckstein values for this option.

.. _D08:
D08 : Chemical Sputter Yield Option
  Chemical Sputter Opt 1: DIVIMP implementation - yield formulae from Garcia-Rosales/Roth (1994)

  Chemical Sputter Opt 2: DIVIMP implementation - yield formulae from Garcia-Rosales/Roth (1996)

  Note: Options 3 to 7 correspond exactly to chemical sputtering options 1 through 5 as they are implemented in NIMBUS. The identical code for calculating these yields has been incorporated into DIVIMP.

  Chemical Sputter Opt 3: NIMBUS option 1 - Garcia-Rosales/Roth (EPS'94)

  Chemical Sputter Opt 4: NIMBUS option 2 - Pospieszczyk (EPS'95)

  Chemical Sputter Opt 5: NIMBUS option 3 - Vietzke

  Chemical Sputter Opt 6: NIMBUS option 4 - Haasz

  Chemical Sputter Opt 7: NIMBUS option 5 - Garcia-Rosales/Roth (1996)

  Chemical Sputter Opt 8: NIMBUS option 6 - Haasz (1997) - (Brian Mech PhD)

  Chemical Sputter Opt 9: Specified CONSTANT yield value.

  Chemical Sputter Opt 10: Based on Haasz (1997) - (Brian Mech PhD Thesis). Modified to reduce the yield to 1/5 of its value as the temperature drops from 10eV to >5eV. Constant at 1/5 below 5eV.

.. _D09:
D09 : Momentum Transfer Collision - First Coefficient (kelighi)
  Probability coefficient for a momentum transfer collision with a background ion. 

.. _D10:
D10 : Momentum Transfer Collision - Second Coefficient (kelighg)
  Probability coefficient for a momentum transfer collision with a background neutral particle. 

.. _D11:
D11 : Characteristic energy Ebd (eV)
  This item is used for NEUT cases where Vel/ang flag is 0,1,4 or 5. Typically one would use 8.3 (7.3*) for Carbon impurity. Reference: Note 41. This value may change as research defines the value more precisely. In the meantime, one can either enter a desired value from the literature or enter a value of 0.0. In the case of zero, the code reads the value from a hard-coded table taken from the literature.

.. _D12:
D12 : Neutral hydrogen density parameter - Nhc (m\ :sup:`-3`)
  The next few items refer to cases where Charge Exchange Recombination has been requested. It is not available for all impurities and should normally only be used in conjunction with a hydrogenic plasma. It can be used with other plasma types, but a warning message will be issued by the program. The parameters here are described in Note 89. The standard value for Nhc is 1.0E15.

.. _D13:
D13 : Nho (m\ :sup:`-3`)
  Also for CX-Recombination. Standard value 3.0E18.

.. _D14:
D14 : Lamhx (m)
  Also for CX-Recombination. Standard value 0.02.

.. _D15:
D15 : Lamhy (m)
  Also for CX-Recombination. Standard value 0.11.

.. _D16:
D16 : Constant for CX Recomb option 2 - Vcx (m/s)
  Also for CX-Recombination option 2 only. For the usual case (option1) Vcx is calculated as sqrt(2Tb/Mb) and varies with Tb(x). For this option we set Vcx to a constant value, typically 2.4E4. See Note 173.

.. _D17:
D17 :  Threshold yield for self-sputtering (eV)
  Required for cases using proper self-sputtering. The impact energies of ions returning to the targets are calculated, and these are then used to calculate Yields "Ys" using a specified set of sputtering data that includes self-sputtering yields. If Ys is greater than the threshold yield given here, then a fraction Ys of a neutral is sputtered off and followed until it too is eventually removed. If it ionizes and returns to the target then a new value of Ys is calculated and may result in a smaller fraction Ys1*Ys2 of a neutral being sputtered off again. A typical value for this parameter is 0.1. See Note 87.

.. _D18:
D18 : Bombarding ion charge state - Zimp
  This item is required for NEUT sputter options 1 and 2, which provide a simple method for modeling self-sputtering (not the same as the proper self-sputtering method, option 3). The ion species bombarding the limiter is defined using this item and the following item "bombion". For example, if "bombion"=5 and this item =4, then we are simulating the bombardment of the limiter with C4+ ions. See Notes 38 and 144.

.. _D19:
D19 : Bombion (0-Mb 1-H 2-D 3-T 4-He4 5-C 6-Mi 7-O)
  This item again required for NEUT sputter option 1. There are two special cases :- Setting this flag to 0 indicates that the bombarding ion type is the plasma background ion Mb, but of course the value of Zimp could be specified as different to Zb; setting this flag to 6 indicates that the bombarding ion type is the impurity ion Mi, ie the limiter material itself, with an appropriate value for Zimp. The remaining values for this flag allow a variety of other ion types to bombard the limiter. Note: the subscript "imp" is always used to mean "impact", while the subscript "i" is always used to indicate "impurity". See Note 144.

.. _D20:
D20 : Ionization rate factor for neutrals - IRF
  Adjustment factor applied to the ionization rate of the neutrals only. Is only required for NEUT cases, and normally is set to 1.0. Entering a value such as 0.2 would reduce the ionization rate of the neutrals. Reference: Note 146.

.. _D21:
D21 : Sputtering Enhancement Factor - SEF
  This is a correction factor applied to the calculation of Z effectives, etc. Normally it would be set to 1.0. Any other value causes the "total primary integrated flux*yield" to be adjusted by this factor prior to being used in the Z effectives formula. Reference: Note 152. 

.. _D22:
D22 : Set of Yield Modifiers for Primary, Secondary and Chemically sputtered neutrals
  One must specify here a set of yield multipliers that will be applied to physical and chemically sputtered particle yields originating from energetic ion impact on the targets and the vessel walls. In addition, one specifies a yield multiplier for self-sputtered particles that originate from the target - self-sputtering is not supported from other vessel surfaces. The last item on the line is a value for the reflection coefficient for neutral particles, with a value between 0.0 and 1.0.

  The "targets" are defined by the grid segments at the ends of the rings in the SOL and Private Plasma. The "vessel walls" are defined by a set of joined line segments which connect the outer corners of the INNER and OUTER targets together. This set of line-segments approximates the geometry of the actual vessel wall. Normally, set the "number of rows of data" to 1 and enter the following values on the following line:

  .. code-block::

    ' ' 'Set of Yield Modifiers for Primary, Secondary neutrals'  
    '      Number of rows of (X,Mpt,Mst,Mct,Mpw,Mcw,Refl) data-'  4 
             0.0   0.0    1.0     1.0    0.3    1.0    0.3    0.0
            36.0 130.0    1.0     1.0    0.3   15.0   15.0    0.0
           165.0 166.0    1.0     1.0    0.3    1.0    0.3    0.0
           167.0 177.0    1.0   -99.0    1.0    1.0    1.0   -1.0

  The leading two numbers are wall indices which specify the range of wall segments to which the specified yield and reflection modifiers will apply. If the first line contains an element labelled as 0.0 as it's first element then this set of values is taken to be the default set for the entire wall and is applied first, before later table entries change the values for specific ranges of wall sections. If there is an overlap in range of segments specified for two regions then the set of yield modifiers later in the list is the one used for any segments that overlap earlier specifications.

  There are six values specified on each line in addition to the segment indices. Each of these numbers is interpreted in different contexts depending on it's value. In order the six numbers represent the following quantities.

  Modifier for Physical Sputtering on Target Segments

    This option multiplies the calculated yield for physical sputtering on the specified wall element if it is a part of a target. A value of 1.0 is usually used unless there is a reason to suspect a change in the effectiveness of the physical sputtering process.

  Modifier for Self-Sputtering

    Self-sputtering multiplier. A positive value will multiply the calculated self-sputtering yield. A negative value in the range of [-50.0,0.0) will be used as an actual value for the yield of the fragment on the surface. Thus a value of -2.0 will result in a fixed yield value of 2.0 being used for this segment no matter what the nature of the impacting particle or it's energy. A value in the range of -99.0 to -100.0 will activate an ion reflection/prompt thermal re-emission mechanism by which the impacting ion will be neutralized and relaunched as a neutral particle with a given fixed energy specified in the input file and with it's trajectory selected from a cosine distribution. A value of -99.0 represents a probability of ion reflection of 1.0. A value of -100.0 is a ion reflection probability of 0.0.
    
  Modifier for Chemical Sputtering on Target Segments

    This option multiplies the calculated yield for chemical sputtering on the specified wall element if it is a part of a target. A value less than 1.0 (often in the range of 0.3 to 0.5) may be used as a method of modelling the quick break-up and prompt redeposition of methane fragments. This may be particularly true for target segments which may have significant plasma contact.
    
  Modifier for Physical Sputtering on Wall Segments

    This option multiplies the calculated yield for physical sputtering on the specified wall element if it does NOT form part of a target. A value of 1.0 is usually used unless there is a reason to suspect a change in the effectiveness of the physical sputtering process.
    
  Modifier for Chemical Sputtering on Wall Segments

    This option multiplies the calculated yield for chemical sputtering on the specified wall element if it is NOT a part of a target.
    
  Reflection Coefficient for Neutral Impact on Segments

    If the neutral reflection option is active then this quantity allows each element of the vessel wall/target to have a different reflection coefficient. A value of 0.0 will deactivate reflection for the given element, even if the refection option is turned on. A negative value for this quantity results in a prompt thermal re-emission mechanism being employed instead of normal reflection. In PTR the emission energy of the neutrals coming from the specific wall segment with a negative yield is specified in the input file by the input energy quantity. Otherwise the energy of the neutral is retained during a surface reflection

.. _D23:
D23 : Specified Fixed Yield Value for Sputter Data Option 4
  When Sputter Data Option 4 is in use the yield for all sputtering events is set to a fixed value. This input line specifies the value that will be used for that fixed value of the yield.

.. _D24:
D24 : Target Temperature in K
  This value is used in the chemical sputtering yield options for calculating the total chemical sputtering yield from each target segment. (Units of Kelvin)

.. _D25:
D25 : Wall Temperature in K
  The wall temperature is used in the same manner as the target temperature. It is a factor in the chemical sputtering yield formulae and is used in calculating the yield for each vessel wall segment. (Units of Kelvin)

.. _D26:
D26 : Private Plasma Wall Temperature
  The private plasma wall temperature is used in the same manner as the wall or target temperature. It is a factor in the chemical sputtering yield formulae and is used in calculating the yield for each vessel wall segment. (Units of Kelvin)

.. _D27:
D27 : Specific Wall Segment Temperatures
  This input contains a list of wall element index numbers and a temperature to be associated with each wall element. A different temperature can be specified for every element of target, wall and private plasma wall. Temperatures may be specified for any range of elements. If a temperature is not specified for a specific element - the temperature for that element is determined by the overall wall, target, and private plasma wall temperatures described previously.

  The format for this input is as follows:

  .. code-block::

    ' ' 'TN1450 Wall Temperatures in K for Specific segments'
    'Number of Segment Ranges (Index1 Index2 Temp):' 3
        35    40    800
       103   107   1000
       116   116   1200

  Each line contains a range of wall element indices followed by a temperature that should be applied to those indices. The actual index numbers that should be used can be obtained by running a case using the appropriate grid and wall options and then looking in the output ".dat" file for the wall element listing. This listing contains the location and the index numbers of each element of the vessel wall. 

.. _D28:
D28 : Temperature Gradient Coefficient parameter - ZO
  For use with TeB and TiB gradient coefficient (`Q01`_, `Q02`_) option 2. See note 412 for details.

.. _D29:
D29 : Emax-factor for NEUT launch velocity - Emaxf
  Required for NEUT cases with Vel/ang flag 4 or 5. This allows us to adjust the cut-off on the velocity distribution at launch. This factor is normally set to 1.0 but can be adjusted up or down (Note 93).

.. _D30:
D30 : Impurity ion impact energy for Wall Launch neutrals
  This is part of the formula for calculating the maximum launch energy available to wall launched neutrals. It represents the mean or average expected impact energy of the particles causing the neutral wall particles to sputter.

.. _D31:
D31 : Maximum Number of Sputtered Generations
  This parameter limits the number of iterations or generations for which a sputter fragment can exist. The purpose of this is to limit the CPU time for cases using recycling gases with high YMF values. Recycling gases and pumping rates are simulated through two mechanisms. First, they are defined to have a sputter yield of 1.0. This means that every particle striking the target would be relaunched. This runaway condition is limited by applying a yield modifier of, for example, 0.95, which gives an effective pumping rate at the targets of 5% of the incident recycling impurity flux. These particles will only be followed until their weight is less than the sputter threshold (usually 0.1) or their weight reaches (YMF)^n where n is the specified number of sputter generations. This is implemented by simply ceasing to follow the ion fragment after the specified number of target collisions. 

.. _D32:
D32 : ABSFAC or Power to Targets
  Normally this quantity will be calculated by the DIVIMP run. However, for purely ion injection cases it is not possible for DIVIMP to calculate an absolute calibration factor for the results. This quantity can be externally calculated or deduced for a given shot or other simulation and then entered into the code to check the absolute values of the results generated by DIVIMP against the experiment. This quantity may instead be interpreted as the Power onto the targets when SOL option 21 is used. This also yields an absolute calibration factor but has the additional effect of easing the definition of the region B power loss in terms of plate target power. The absolute factor would have values on the order of 10\ :sup:`19` for ABSFAC and 2.0 x 10\ :sup:`5` for the power onto the targets.

.. _D33:
D33 : STGrad - Gradient End-point Specifier
  This quantity is used for some simple options to arbitrarily turn off the temperature gradient force or velocity diffusion for values of S > stgrad * smax ... and yet still allow a temperature gradient for S > stgrad * smax and temperature gradient forces for the region less than Stgrad * smax (from both targets). This option was used in testing some specific modelling conditions involving various collision options and force balance effects and is not recommended for regular use. 

.. _D34:
D34 : H Recombination Calculation Option
  This option specifies the formula to be used to calculate the Hydrogenic Recombination. This term is calculated in DIVIMP because it is not passed back from NIMBUS. (It is available from EIRENE). The option specified here is passed to NIMBUS as an input to ensure that the same calculations are made in DIVIMP and NIMBUS for the hydrogenic recombination source. The following entry is the specification of a lower temperature limit that may be used in the calculation of the recombination. Some of the expressions or data used to calculate the recombination may not be reliable at very low temperatures. As a result, this mechanism was introduced to prevent this term from becoming exceedingly large at very low temperatures as seemed to be the trend in some cases.

  Option 0: OFF Hydrogenic Recombination is OFF

  Option 1: Gordeev Gordeev coefficients - as implemented in NIMBUS (May, 1996) - are used to calculate the recombination.

  Option 2: Janev Janev coefficients - as implemented in NIMBUS (May, 1996) - are used to calculate the recombination.

  Option 3: NRL NRL coefficients - as implemented in NIMBUS (May, 1996) - are used to calculate the recombination.

  Option 4: ADAS ADAS coefficients are used to calculate the recombination.

.. _D35:
D35 : H Recombination - Limiting Cutoff temperature (eV)
  This is a lower limit cut-off temperature, implemented in DIVIMP, which is applied in the calculation of the recombination in the above formulae. Usually, this quantity is set to 0.0 so that it does not interfere with the calculations. It should only be used when a specific situation where a specific cell with an exceptionally low temperature contains so much recombination that it prevents the solver from working altogether. 

.. _D36:
D36 : T-Grad Modification Factor
  This factor is used in the UEDGE and Garching correction formulae for the temperature gradient forces. These corrections are intended to compensate for kinetic effects which would be expected to weaken the temperature gradient force which is dependent on a good collisional coupling for energy transfer along the field lines. 

F Tags
------
.. _F01:
F01 : Read FLUID CODE (e.g. EDGE2D/UEDGE) Background for Reference
  Edge2D Reference Option 0: OFF- An Edge2D background will not be read.

  Edge2D Reference Option 1: ON- An Edge2D background - whose name was specified to the DIVIMP execution script - will be read in and stored for comparison to the calculated DIVIMP background plasma.

  Edge2D Reference Option 2: ON - SONNET based (i.e. B2 or UEDGE) background plasma is read in for reference.

  Note: The purpose of this option is to allow easy comparison between DIVIMP calculated background plasma results and equivalent Edge2D results. The values read in are passed to the OUT program and are then plotted on the same graph with the DIVIMP values in order to provide a direct method of comparison.

.. _F02:
F02 : Fluid Code Target Condition Usage Option
  This option will read the target conditions from an EDGE2D case using one of five different methods and make them available to the DIVIMP SOL background routines. It also allows these conditions to be extracted and passed to NIMBUS.

  FLUID CODE Target Opt 0: OFF. FLUID CODE data that has been read in not used to assign the initial target conditions.

  FLUID CODE Target Opt 1: ON.

  EDGE2D data is read using the standard method and is applied as the target conditions. The standard EDGE2D method uses the following formulae:

    N-target = N-centre-of-first-cell

    Ti-target= Ti-centre-of-first-cell

    Te-target= Te-centre-of-first-cell

    V-target = E2D-mach-number * Cs

  For UEDGE or other B2 like codes the data is read using the standard method and is applied as the target conditions. The standard UEDGE method uses the following formulae:

    N-target = N-centre-of-first-cell

    Ti-target= Ti-guard-cell

    Te-target= Te-guard-cell

    V-target = V-cell-edge (for each target)

  EDGE2D Target Opt 2: ON. EDGE2D data is read from the plasma file. In addition the EDGE2D down fluxes MUST be supplied so that the target conditions can be properly extracted.

    V-target,Te-target and Ti-target as option 1.

    N-target = G-target/V-target

  EDGE2D Target Opt 3: ON. EDGE2D target values are calculated as in option 1. However, the EDGE2D down fluxes are also loaded and are passed to NIMBUS as the target fluxes for each segment. The target conditions and fluxes will not necessarily agree.

  EDGE2D Target Opt 4: ON. EDGE2D data read for target. EDGE2D down fluxes are also read.

    N-target,Te-target and Ti-target as option 1.

    V-target = G-target/N-target

  EDGE2D Target Opt 5: ON. EDGE2D data read for target. EDGE2D down fluxes for both particles and power are also read. The target conditions are calculated as described in option 2.

  EDGE2D Target Opt 6: ON. EDGE2D data read for target. This is the same as option 1 except that the velocity is explicitly set to the sound speed based on the extracted target conditions rather than being based on the EDGE2D mach number. 

.. _F03:
F03 : Lost SOL ring option
  Lost SOL Rings 0 : Plasma Set to minimum values

  Lost SOL Rings 1 : Plasma Set to outer ring value

.. _F04:
F04 : Velocity Multiplication Factors - for data read from files
  This section specifies a block of data that defines factors to multiply the background velocity. The format of the input is as follows:

  .. code-block::
  
    ' ' ' Edge1D/2D Deuterium drift vel. mult. factor VMF '
    ' Number of VMF blocks ' 0 
    ' Ring range :- ' -20  -30
    ' J0 and J1 :-  '   5    5
    ' VMF0,1,2 :-   '  1.0  1.0  1.0

The meaning of the various quantities is...

    Number of VMF blocks. Defines the number of different sets of VMF data that will be entered. Each set of data consists of the three data lines listed above. Note that even when 0 sets of data are specified - at three lines of data must appear even though the information they contain is ignored. One reason for this is to keep the data entry format as part of the data file since the input format is somewhat cryptic. The ring range specifies the set of rings over which the data will be applied. These can be actual ring numbers, if such are known, or a set of symbolic negative numbers that specify specific rings. The negative numbers correspond to the following rings:

    Ring No. < 1 1 (First main plasma ring)

    -10 IRSEP -1 (Last main plasma ring)

    -20 IRSEP (Separatrix - or first SOL ring)

    -30 IRWALL (Wall ring - or last SOL ring)

    -40 irwall+1 (Trap wall - or first trap ring)

    -50 NRS (Last ring - innermost trap ring)

  The quantities J0 and J1 define the number of bins from each end of the ring over which the multipliers will be applied. So, values of 5 and 5 define the first 5 bins as region 0 and the last 5 as region 2. The three VMF factors specify the velocity multiplier that will be applied in each of these regions. In region 0 and 2 the value scales linearly up to the quantity in region 1 which will in turn apply over the central region.

.. _F05:
F05 : Sonnet Grid: Number of Fluid results in background plasma file (NFLA)
  This specifies the number of fluid solution that will be found in the B2/UEDGE background plasma solution that will be read by the case. Specifying this number enables the input file to be read correctly and will allow the code to load the impurity fluid results for comparison with the DIVIMP results.

.. _F06:
F06 : Read Background Plasma Auxiliary Input File
  This switch instructs DIVIMP to try to read in an auxiliary plasma input file. This file will typically contain additional information that is not included in the plasma file. For example, at this time, when this option is activated for a JET grid and EDGE2D background , DIVIMP will look for a file connected to unit #12 that contains the EDGE2D down flux information for the case. The expected name for this file is "shots/<gridfile>.<plasma descriptor>.aux". For an UEDGE case, the auxiliary file is expected to contain a variety of information in a B2 plasma file format. This includes such items as the hydrogenic neutral density, the carbon neutral density and the 2-D carbon neutral production rate on the grid among other items.

G Tags
------
.. _G01:
G01 : Grid Option
  Grid Option 0: Standard jet grid files

  Grid Option 1: Standard Asdex grid files. As implemented at Juelich (KFA)

  Grid Option 2: Standard ITER grid files. Works with one example of a double null ITER grid - processed at Juelich (KFA).

  Grid Option 3: Standard Sonnet grid files. Works with standard Asdex upgrade/B2 grid files. The number of points to the cut point and other parameters are entered as parameters in the input file. This is also the standard grid type usually used for DIIID, TdV and CMOD modeling.

.. _G02:
G02 : Non-Orthogonal Option
  Non-Orthogonal Opt 0: Standard. All transport and target angles are treated as orthogonal. No corrections for non-orthogonality are made.

  Non-Orthogonal Opt 1: JET. Targets and transport are treated using corrections for non-orthogonal grids. Non-orthogonal transport is implemented using ancillary information available with JET grids only.

  Non-Orthogonal Opt 2: Targets only. Target fluxes and other factors are corrected for target grid orthogonality in a manner identical with option 1. Ion transport is treated as orthogonal with no corrections made for grid orthogonality.

  Non-Orthogonal Opt 3: Generalized Non-orthogonal treatment. Both targets and ion transport are corrected for grid non-orthogonality. The non-orthogonal ion transport is implemented by the calculation of an additional orthogonal co-ordinate which is held constant when moving cross-field. This co-ordinate is calculated based on individual cell orthogonal characteristics. (e.g. center angles) and is functionally identical to the additional grid information that is available for JET grids.

.. _G03:
G03 : Parallel Distance Option
  Parallel Dist Opt 0: Cell Centers. This option affects particle accounting and ion transport and in addition should be selected in combination with cell area option 0. The boundary between cells for particle accounting purposes is half-way between the centres of the adjacent cells. The S-distances along the field lines are calculated by joining the centers of adjacent cells.

  Parallel Dist Opt 1: Polygon Boundaries. This option affects particle accounting and ion transport and in addition should be selected in combination with cell area option 1 and cross-field distance option 1. The S-distances along the field lines are calculated by joining the mid-points of the ends of the polygon that cross the field line to the center point of the cell. The S-coordinates of both the cell center and the polygon boundaries are recorded. An ion is in a specific bin if the S position of the ion lies between the S-boundaries for the cell.

.. _G04:
G04 : Cross-field Distance Option
  Cross-field Dist Opt 0: Cell centers. This option affects particle accounting and ion transport. It should be used in combination with area option 0 and parallel distance option 0. A particle is considered to have cross-field diffused into the next cell when it crosses the half-way point between the cell centers moving inward or outward.

  Cross-field Dist Opt 1: Approximate Polygon Boundaries. This option affects particle accounting and ion transport. It should be used in combination with Area Option 1 and parallel distance option 1. A particle is considered to have cross-field diffused into the next cell when it has stepped farther than distance of the intersection point of the polygon boundary with the line joining the two cell centers of the adjacent cells.

.. _G05:
G05 : R,Z Calculation Option
  R,Z Option 0 : Cell center R,Z values are used to estimate the particle position.

  R,Z Option 1 : The GETRZ subroutine is used to calculate an estimate of the actual R,Z position of the particle. At present this estimate only includes the parallel displacement and not the perpendicular component because of the difficulties in defining a perpendicular angle consistently through-out the cell.

.. _G06:
G06 : XY Grid Option
  XY Grid Option 0: Off- XY grids are NOT used to track impurity neutrals in DIVIMP. A bin-finding subroutine is used instead. The rectangular grid option described later in the code is NOT used.

  XY Grid Option 1: On- Use XY grids to track impurity neutrals. the rectangular grid option specified later is used to define whether the grid will be calculated or loaded.

.. _G07:
G07 : Cell Area Calculation Option
  Cell Area Option 0: Approximate. This is the original DIVIMP method that calculates cell areas based on the locations of the cell centres. It must be used with older grids that do not include the all the grid cell information. It is not recommended for use with current polygonal grids.

  Cell Area Option 1: Polygonal. This option uses the complete cell polygon information, specifically the co-ordinates of the corners, to calculate the proper area of each cell on the grid. It can not calculate areas for the virtual rings for which complete cell information is not available. 

.. _G08:
G08 : Ion Wall Option
  The Ion Wall options specify the boundaries for ion transport in the DIVIMP code. This option is combined with the target option to define the region of allowed ion transport. This option applies to both the main SOL outer wall and private plasma outer wall definitions for ion transport.

  **Ion Wall Option 0**: Ion walls mid-way between the center points of the outermost two rings on the grid. The outermost ring on both JET and Sonnet based grids is virtual and is used only to anchor the fluid solutions. As such, the plasma space mapped by the grid only extends to the midpoint between this virtual ring and the last real ring contained within it.

  **Ion Wall Option 1**: Ion walls are at the outermost ring of the SOL. The center points of the virtual ring are used to define the ion wall position.

  **Ion Wall Option 2**: Ion walls are located at the polygon edge of the outermost complete ring of polygons defining the plasma region. This is almost equivalent to option 0 except that it uses the actual polygon definition of the last real ring to construct the outer wall for ion transport. Option 0 was the original option for DIVIMP grids that contained only the centre points and not the compete polygonal description of the grid.

.. _G09:
G09 : Neutral Wall Option
  The Neutral Wall options specify the boundaries for neutral transport in the DIVIMP code. In all cases the different neutral wall options link to the target and private plasma wall options to define a vessel boundary for neutral transport. This definition applies to the neutral walls bordering the main SOL region. The private plasma wall region for neutrals is defined in the next option.

  **Neut Wall Option 0**: Neutral walls are half-way between the center points of the outermost two rings on the grid.

  **Neut Wall Option 1**: Neutral walls are created by joining the center points of the outermost (virtual) ring.

  **Neut Wall Option 2**: Neutral walls are specified by a set of coordinates entered in the input file.

  **Neut Wall Option 3**: Neutral walls are specified by a set of points that have been hard-coded for specific JET grids in the subroutine loadgeo - located in the pindiv.d4a module. This option is old and has been superseded by options 4 and 5 - however, it may prove useful if it is ever necessary to often use a specific grid for which pre-generated wall data is either unavailable or inaccurate.

  **Neut Wall Option 4**: Neutral walls are specified by the vessel coordinates that are specified in the GRID2D geometry file that is read-in by DIVIMP.

  **Neut Wall Option 5**: Neutral walls are read-in from the PIN/NIMBUS transfer file. Although, the position of the walls is the same as in the GRID2D file, PIN/NIMBUS generally sub-divides the vessel surface for its own purposes in calculating data - as such it is necessary to load the PIN/NIMBUS version of the wall specification when using PIN/NIMBUS wall data in conjunction with DIVIMP. Use of this option requires that PIN/NIMBUS be run from within DIVIMP.

  **Neut Wall Option 6**: This option deals with walls whose segments are ordered counter-clockwise. It will invert the ordering of these wall segments so that they conform to the clockwise standard required by DIVIMP. This option then sets the neutral wall option to be option 4. This option should ONLY be used for VERY old JET grids where this is known to be a problem.

  **Neut Wall Option 7**: The main wall for neutrals is located at the outermost polygon edge of the last real (non-boundary) ring of the grid.

.. _G10:
G10 : Trap Wall Option
  **Trap Wall Option 0**: The wall for neutrals in the trap region is midway between the outermost two rings.

  **Trap Wall Option 1**: The wall for neutrals in the trap region is at the outermost ring.

  **Trap Wall Option 2**: The wall for neutrals in the trap region is created by joining the private plasma end-points of the two plates with a straight line. The end-points are defined to be midway between the outermost two rings of the private plasma region in order to remain compatible with PIN.

  **Trap Wall Option 3**: The neutral wall in the trap region is specified by a series of line segments entered in the data file just after the entry for the outer wall specification. These points are joined to the end/corner points of the target to form a continuous outer boundary for neutrals.

  **Trap Wall Option 4**: The private plasma wall for neutrals is specified by a set of additional coordinates for the vessel wall taken from the GRID2D geometry file.

  **Trap Wall Option 5**: The private plasma wall for neutrals is read from the PIN/NIMBUS transfer file.

  **Trap Wall Option 7**: The private plasma wall for neutrals is located at the outermost polygon edge of the last real (non-boundary) ring of the grid.

.. _G11:
G11 : Vessel Wall Re-Definition Option
  The purpose of this option is to instruct DIVIMP to include any baffles that may be specified in the grid file as part of the vessel wall. The code follows the main vessel wall until a section where a baffle comes very close to the wall is found. The code then redefines the wall by following around the outside of the baffle until it rejoins the wall. This option is only meaningful in conjunction with wall option 4 and/or trap wall option 4. If trap or wall option 5 are specified, in which the wall definition is taken from the PIN/NIMBUS output, then the wall will automatically be redefined if baffles are present.

  **Vessel Redef Option 0**: OFF - Vessel will not be redefined.

  **Vessel Redef Option 1**: ON - Vessel will be redefined to include adjacent baffles. 

.. _G12:
G12 : Target Position Option
  **Target Option 0**: Target is located at second grid points on the SOL and trap rings. Virtual points are discarded.

  **Target Option 1**: Target is located mid-way between the centre points of the virtual cell and first real cell on the SOL and trap rings. Virtual cell centers are then discarded.

  **Target Option 2**: Target is specified by a set of points entered in the data file. One point for each end of each ring. Virtual points are discarded.

  **Target Option 3**: Target is specified by a set of points that are hard-coded. One point for each end of each ring. The set of points is selected by the geometry option. The virtual points are discarded.

  **Target Option 4**: Target is located at the center of the first (virtual) cell on the SOL and trap rings.

  **Target Option 5**: Target is specified by a set of points entered in the data file. One point for each end of each ring. Virtual points are not discarded.

  **Target Option 6**: Target is specified by the polygon boundaries of the last of the real plasma cells on each ring. The virtual cells are discarded. This option should be the one used with all grids for which complete polygon information is available since the target in these cases should correspond exactly with the polygon edges.

.. _G13:
G13 : Geometry Selection Option
  **Geometry Option -1**: Wall and Target position data are to be read from the input file for the target and wall options specified above if such action is appropriate.

  **Geometry Option 0**: Shot 24719 (one specific grid). Hard-coded target and wall data are available for the appropriate target and wall options.

  **Geometry Option 1**: Shot 26308 (one specific grid). Hard-coded target and wall data are available for the appropriate target and wall options. 

.. _G14:
G14 : Ring Location of Core Mirror - IRCORE
  This quantity allows the core ring, which will be used as the core mirror for impurity ion transport, to be specified. Any ions reaching this ring will be reflected. 

.. _G15:
G15 : Rectangular grid for neutrals 0calculate 99file
  This option is NOT needed anymore, unless the XY-grid option mentioned earlier in this document has been turned on. The XY grid was an older method for efficiently determining the grid element which an ion currently occupied. It was used originally because the cell polygon corner information was not available to DIVIMP and grid positioning was based on the closest cell center. Now that this information is generally available and used by DIVIMP, the XY grids are no longer required. This option remains in order to maintain compatibility with older grids.

  This option should be set to 0 whenever some new shot data is first used in order to create a set of index arrays relating the (ik,ir) grid to a straightforward rectangular (ix,iy) grid used for following neutrals. This calculation can be quite time consuming, especially if a fairly fine rectangular grid of say 201 by 200 elements is being used. Once the file has been created, set this option back to 99 for all subsequent runs using the same geometry data and connect the file to Unit number 13. In general, it has become preferable to use the grid position routine that uses the actual coordinates of the bin vertices to determine which bin a particle is in. This may require slightly more computational time while processing the particles, but it is more than compensated for by the greater accuracy in assigning the initial bin of the particles. Furthermore, for larger grids, it had become necessary to use XY grids with resolutions exceeding 1000 by 501 elements. Index arrays of this size require a significant amount of storage and CPU time to compute. Note, however, that the XY grids are still used in the OUT program for mapping the results onto an even XY grid for passing to the plot routines. As before, the indexing arrays are calculated only once and are then reused for subsequent runs. If one wishes to return to the previous system, it is necessary to adjust the values of parameters in the PARAMS, DIVXY and OUTXY common blocks.

.. _G16:
G16: Set of Target Coordinates
  This table specifies target co-ordinates for each ring. The first line is a comment. The second contains a comment and the number of lines in the table. The rest of the lines contain the co-ordinates of the targets for each ring. For this data to be used, the target option must be set appropriately and the geometry option set to -1. The contents of each line in the table are:

  Ring number R,Z Outer target R,Z Inner target
  The R,Z values are specified in meters. There are then 5 numbers on each line of the table. Target points must be specified for each ring.

.. _G17:
G17 : Set of Wall Coordinates
  This table specifies the set of wall segments that will make up the wall. The first line is a comment, the second is a comment plus the number of lines in the table, the rest are the R,Z coordinates of the end-points of the wall segments, one end-point per line. The end-point of the first section is the corner point of the outer plate. The sections then number clockwise from that point to the outer corner point of the inner plate. The wall in the private plasma cannot (at this time) be independently specified. Furthermore, the targets are assumed to be part of the wall.

.. _G18:
G18 : Set of Trap Wall or ITER second wall coordinates
  This table specifies the points that will be used either for trap wall option 3 - where a neutral wall is specified in the trap region or for the second wall in an ITER double null geometry. The first line in the input is a comment, the second is a comment plus the number of lines in the table, the rest are the R,Z coordinates of the end-points of the wall segments, one end-point per line. These points are joined by line segments and connected to the inside corners of the targets to give a continuous wall surface for neutrals inside the trap region. 

.. _G19:
G19 : Sonnet Grid Characteristic Specifications - Asdex U - CMOD - Textor - DIIID : Number of Rings
  These quantities are now usually included at the beginning of a SONNET style grid file. If there are no numbers present in the grid file then the values specified here will be used. The values in the grid file will always take precedence over numbers specified in the input file.

  The following set of values specify the characteristics of grids created by the Sonnet mesh generator. This generator creates grids compatible with the B2/Eirene code combination. Thus the grid is initially stored in an indexed XY array that needs to be divided into three segments for DIVIMP. The core plasma, the trapped or private plasma and the SOL. This first parameter specifies the number of rings in the Sonnet geometry file. (This would typically be the number of rows in the Sonnet mesh- since the Sonnet mesh is indexed from 0 to n - the number entered here is n+1)

.. _G20:
G20 : Sonnet Grid: Number of Knots
  This quantity specifies the number of knots along each ring. It is typically larger than the number of rings but there is no requirement that this be true. Again, since Sonnet numbers from 0 to m, the number entered here will be m+1.

.. _G21:
G21 : Sonnet Grid: Cut ring number
  This is the last ring composing the trapped and core plasma regions. (The last ring containing a B2 cut line.) Note: the terminology comes from B2 where the square geometry has two regions - the SOL is usually the upper half and the trapped plasma and core are the lower half. However there is no direct coupling between the trapped plasma and the core and so there are "insulating" lines (cut lines) inserted separating these regions on these rings. That is why they are termed "cut rings". DIVIMP needs to split these rings into core and trapped plasma rings which DIVIMP deals with separately and so needs to know the ring number to stop this procedure at - again since Sonnet grids are numbered from 0 - this number is the cut ring index number +1.

.. _G22:
G22 : Sonnet Grid: Cut point 1
  This index specifies the last element belonging to the trapped plasma region on the cut rings. All the elements greater than this and less than cut point 2 belong on the core plasma rings.

.. _G23:
G23 : Sonnet Grid: Cut point 2
  This specifies the first point after the second cut belonging to the trapped plasma rings. (see above.)

H Tags
------
.. _H01:
H01 : PIN/NIMBUS Random number seed (0 generate new seed)
  The random number generator is reset using a specific seed value before starting any PIN/NIMBUS. A value of 0 will generate a new seed for each PIN/NIMBUS run based onthe current date and time. A value of -1 will use the PIN/NIMBUS default seed value of 1 for all runs. A value greater than 0 will use the specified number for the PIN/NIMBUS seed. Occasionally it may be desirable to try and reproduce a case exactly (for debugging perhaps) when the random number seed can be read from the printout of the previous run and inserted here for the new run. In this case the given seed is used in place of any generated seed. Historically, PIN/NIMBUS has been run using a fixed seed value that does not vary between runs.

.. _H02:
H02 : PIN/NIMBUS Print option (0 regular, 1 more data)
  This option provides a means of selecting enhanced print outs for a PIN/NIMBUS run. Option 0 will result in a regular PIN/NIMBUS print out while option 1 will cause additional information on the PIN/NIMBUS run to be included in the output files. Option 0 is all that is required for most purposes. 

.. _H03:
H03: Run PIN Option
  Run PIN Option 0 : PIN is not executed from within DIVIMP.

  Run PIN Option 1 : Execute PIN from inside DIVIMP and import results.

  This on/off switch instructs DIVIMP, when running JET grids, to invoke PIN and pass it the background plasma characteristics. PIN then invokes NIMBUS, a hydrogenic neutral code that follows the hydrogen background given the initial DIVIMP plasma and generates the expected distribution of hydrogen ionization as well as carbon sputtering, charge exchange and several other quantities. This information can then be imported back into DIVIMP and used for several different purposes. (e.g. Recalculation of the background plasma using the PIN source ionization, calculation of charge exchange probabilities using the hydrogen neutral densities, plotting of Halpha radiation distributions, Carbon ion launch based on the distribution of primary Carbon ionization and several other possibilities.)

  Note: PIN can be run from DIVIMP in a number of ways. The whole PIN module can be bound to DIVIMP and invoked as a subroutine. Alternatively, PIN can be invoked as a stand-alone program that is called from inside DIVIMP and while running, the DIVIMP process is not active. Each procedure has different advantages and disadvantages and the option chosen is architecture and processor specific.

  If this option is invoked when DIVIMP is being run on a Sonnet grid then DIVIMP will write a background plasma file in the B2/EIRENE format and then invoke EIRENE to calculate the hydrogenic behaviour. EIRENE is an alternate Monte Carlo hydrogen code that is available and runs on Sonnet grids. It reads it's input from a data file in the EIRENE subdirectory. This information in this data file needs to be tailored to the specific machine specifications and contains all of the simulations parameters required by EIRENE. The name of this file is passed as a command-line parameter to the reirediv script. If a datafile name is not specified it uses a default file called "asdex.dat"

.. _H04:
H04 : PIN Command Line
  This is not a simple optional choice or specifiable value. It should contain the UNIX command to execute PIN with the appropriate options.

  e.g. rpindiv

  or

  reirediv

  The name (rpindiv) is the name of the command file that will execute PIN. The character string in UNIX is passed directly to the operating system by the SYSTEM subroutine and is then executed. The input for this option would be ignored if PIN has been bound to the DIVIMP executable and is being called as a subroutine. The script file obtains the name of the case and other information from environment variables that are set when the script that runs DIVIMP begins executing.

.. _H05:
H05 : PIN Cell Area Option (IHCORR)
  NIMBUS has two methods for calculating the cell areas. These two methods do NOT give the same results either for cell area or for final results from the neutral code. As a result we have added this as a DIVIMP input parameter that is passed directly to NIMBUS. IHCORR=0 will instruct NIMBUS to use the EDGE2D cell areas based on the formula Rho * Theta * Drho * Dtheta while IHCORR=1 instructs NIMBUS to use the cell polygon areas.

.. _H06:
H06 : PIN Hybrid Wall Option
  This option applies only to PIN runs on JET grids. It selects one of a limited set of hybrid wall specifications. The option seeks to simulate an "average" wall for JET by allowing for a wall with greater gaps between the plasma and the wall than is done when the grid file is created. The hybrid wall file contains wall descriptions for both MKII and MKIIa divertor configurations.

  Hybrid Wall Option 0: Off. This option specifies the use of the Standard Wall in the grid file.

  Hybrid Wall Option 1: On. Use the MKII hybrid wall specified in the ancillary input file "hybrid.dat".

  Hybrid Wall Option 2: On. Use the MKIIa hybrid wall specification found in the same input file as option 1.

  Hybrid Wall Option 3: On. Use the corrected MKIIa hybrid wall specification found in the same input file as option 1.

  Hybrid Wall Option 4: On. Use the modified corrected MKIIa hybrid wall specification found in the same input file as option 1.

.. _H07:
H07 : PIN Puffing Option
  This option specifies the form of the PIN/NIMBUS puffing that will be used if NIMBUS is run. Puffing is the process of adding additional sources to the hydrogenic Monte Carlo solution that is produced by NIMBUS. (See: TN1421, 10 Sept. 1996; Memos of L.Horton of 6 Aug. 1996, 16 Aug. 1996, 20 Aug.1996, 9 Sept.1996, 27 Oct.1996,)

  PIN/NIMBUS Puffing Option 0: Puffing is OFF.

  PIN/NIMBUS Puffing Option 1: Puffing is ON.

    Particles lost to the Albedo (PUMP) escape regions are re-injected with characteristics specified by the following quantities. The puffing will only occur on subsequent iterations of PIN because the neutral code does not know how many will be lost until it has finished and thus can not puff these at the beginning.

    Puff Fraction (proportion of losses reinjected) = HPCPUF (`H09`_)

    Puff Temperature (eV) = TPUFH (`H11`_)

    The location of the puff is specified by the location option in combination with the segment specifiers.

  PIN/NIMBUS Puffing Option 2: Puffing is ON.

    A portion of the target flux is injected as a puff instead.

    Fraction of target flux to be puffed = PPCPUF (`H10`_)

    Puff Temperature (eV) = TPUFH (`H11`_)

    The location of the puff is specified by the location option in combination with the segment specifiers.

.. _H08:
H08 : PIN Puff Location Switch
  This switch determines the general location of the puffed particles. The segment specifiers allow the region of the source to be fine-tuned.

  PIN/NIMBUS Puff Location option 0: From main SOL walls.

  PIN/NIMBUS Puff Location option 1: From private void walls.

.. _H09:
H09 : PIN Puff Fraction - HPCPUF
  This parameter applies to puff option 1 only. It specifies the fraction of all particles lost to the Albedo regions that will be re-injected by puffing. Thus a value of 1.0 will result in all Albedo losses from the system being re-injected. This will not work on the first PIN iteration - only on subsequent ones.

.. _H10:
H10 : PIN Flux Puff Fraction - PPCPUF
  This parameter applies to Puff option 2 only. It specifies a fraction of the total target flux that is redirected into a puffed hydrogen source. A value of 0.0 should result in no extra source even if puff option 2 is in use.

.. _H11:
H11 : PIN Puff Injection Temperature (eV) - TPUFH
  This is the initial temperature of the puffed hydrogen atoms. (Note: It is not clear at this time whether NIMBUS actually puffs these particles as atoms or molecules. It is believed that they are puffed as atoms but there is no confirmation of this at this time.)

.. _H12:
H12 : PIN Puff Location Indices - JHPUF1
  The exact meaning of these numbers is difficult to determine before running a case and becoming familiar with the particular grid. The following quantities - JHPUF1(1), JHPUF2(1), JHPUF1(2) and JHPUF2(2) are integers that define the indices of the wall segments from which puffing will occur. There are two relationships that need to be defined before these numbers can be interpreted. First, how are the indices of the outer grid segments defined inside the NIMBUS/EDGE2D code for an arbitrary grid? Second, how are the numbers specified here employed to determine valid puff segments?
  
  1) NIMBUS puffs neutrals in these cases from the outermost boundary of the plasma grid. The neutrals are then followed towards the outer wall - reflect from the wall - and then re-enter the plasma or are removed by other methods. (i.e. Albedo or core loss). The index numbers of these segments are defined in the following manner.
  
  All rings have the same number of cells. There is a certain number of cells along the separatrix ring on both the left and right divertor legs that are adjacent to private plasma cells. The rest are adjacent to cells in the core plasma. The number of the cell that is first adjacent to a core cell is designated as JPRGT. The number of the last cell on the separatrix adjacent to a core cell is called JPLFT. All of the cell indexing for main SOL wall launches is done relative to these two positions. In addition, although the private plasma rings in DIVIMP are treated as separate and smaller rings - in the fluid codes (EDGE2D or B2) and in the hydrogenic Monte Carlo codes (NIMBUS or EIRENE) - the private plasma rings are extensions of the corresponding core rings. Each of these combined private plasma + core rings has the same number of cells - in these other codes - as the main SOL rings. As a result, designation of the puffing segments in the private plasma is handled in a somewhat similar fashion to main SOL puffing.
  
  2) The numbers specified in these inputs are interpreted in the following way to determine at which boundary segments the puff should occur.
  
  In the case of main SOL puffing - Segments J - satisfying the following relationships are selected for puffing.
  
  A) JPRGT + JHPUF1(1) ( J ( JPRGT + JHPUF2(1)
  
  or JPLFT - JHPUF2(2) ( J ( JPLFT - JHPUF1(2)
  
  Examples: (1) To puff from the corner of the right target up to the X-point.
  
  JHPUF1(1) = -1000 JHPUF2(1) = -1
  
  ( Using -1000 will guarantee that all knots on the boundary of the grid from the right X-point down to the target will be selected. (Unless one is using an incredibly large grid.) A value of -1 will start the puffing at the first segment whose corresponding cell on the separatrix is adjacent to the private plasma. Unless the grid is very non-orthogonal JPRGT and JPLFT will roughly correspond to the cells just up from the Xpoint.
  
  ( To puff from the corner of the left target up to the X-point simply set
  
  JHPUF1(2) = -1000 JHPUF2(2) = -1
  
  ( This will do exactly the same as the above except for the left target
  
  (2) To puff from the entire main SOL wall.
  
  JHPUF1(1) = -1000 JHPUF2(1) = 1000
  
  ( This will select the entire wall.
  
  ( Wall segments will not be selected "twice" if overlapping regions from the right and left targets are specified. To turn OFF puffing relative to one of the conditions just specify values for JHPUF that can not be satisfied.
  
  JHPUF1(2) = 1 JHPUF2(2) = 0
  
  ( For the left target - this would require that JPLFT ( J ( JPLFT -1 which can not be satisfied and thus no segments would be selected for puffing through this condition.
  
  (3) To puff from the top of the torus. This can be specified relative to either target and knowledge of the values of JPRGT and JPLFT as well as what numbers correspond to what boundary elements is required before these can be specified with confidence. IF - (for example) - JPRGT = 10 and JPLFT = 30 then the middle segment between these two might be roughly 20. To select this segment for puffing one could choose to set the following values.
  
  JHPUF1(1) = 10 JHPUF2(1) = 10
  
  JHPUF1(2) = 1 JHPUF2(2) = 0
  
  ( These should select the segment numbered 20 relative to the right target. To puff from both the 20th and 21st - set the values in the first line to the following.
  
  JHPUF1(1) = 10 JHPUF2(1) = 11
  
  In the case of puffing in the private plasma the following relationship is used to define which segments to puff from. The values entered for JHPUF1(2) and JHPUF2(2) have no relevance for a private plasma boundary puff source.
  
  B) J < JPRGT or J > JPLFT
  
  and 1 + JHPUF1(1) ( J ( NR - JHPUF2(1)
  
  NR is the total number of segments on the ring.
  
  Examples: (1) To puff from the entire private plasma boundary.
  
  JHPUF1(1) = 0 JHPUF2(1) = 0
  
  (2) To puff from the right side of the private plasma boundary
  
  JHPUF1(1) = JPRGT JHPUF2(1) = 0
  
  (3) To puff from the left side of the private plasma boundary
  
  JHPUF1(1) = JPLFT JHPUF2(1) = 0
  
  (4) To puff from the middle segments of the private plasma boundary
  
  JHPUF1(1) = JPRGT-1 JHPUF2(1) = JPLFT+1
  
  These are indices that specify the bounds of the region from which the extra puffing will occur. The general region is specified by the location option. These parameters allow the specific segments or range of segments to be selected. The default values of JHPUF1 and JHPUF2 together should result in the entire main plasma wall or the entire private void wall being specified as the injection region.
  
  The default entry for this quantity is:

  .. code-block::

    ' PIN Puff Location Indices JHPUF1 (1 and 2) ' -1000 -1000

.. _H13:
H13 : PIN Puff Location Indices - JHPUF2
  See the previous entry for a detailed description of the meaning and interpretation of these parameters.
  
  These are indices that specify the bounds of the region from which the extra puffing will occur. The general region is specified by the location option. These parameters allow the specific segments or range of segments to be selected. The default values of JHPUF1 and JHPUF2 together should result in the entire main plasma wall or the entire private void wall being specified as the injection region.
  
  The default entry for this quantity is:

  .. code-block::
  
    ' PIN Puff Location Indices JHPUF2 (1 and 2) ' -1 -1

.. _H14:
H14 : Nimbus Namelist Input - NIMBIN
  The following block of entries are not used by DIVIMP - instead they are used by NIMBUS when a PIN run is requested. Since DIVIMP ignores these values and NIMBUS will ignore everything else except the namelist input - it is safer and more efficient to include the parameters for the NIMBUS run in the DIVIMP input file. Under these circumstances, it is always clear what options were used to run Nimbus/PIN with a specific DIVIMP case. In addition, the entries in the namelist are now analyzed by DIVIMP and printed in the data file with a description and their usual default values. These items are not printed or analyzed if NIMBUS is not run.
  
  The following is a list of the namelist parameters for PIN/NIMBUS which DIVIMP recognizes and for which it can produce a simple description. The contents of this list are neither complete nor extensive but simply serve as a quick reference for DIVIMP users using the PIN/NIMBUS code. The default value for an entry (the value that is assigned inside NIMBUS if this parameter is NOT specified) is included in [ ] if it is known.
  
  NHIST: Number of Neutral Histories [2000]
  
  IFCHAN: 0=No Channels 1=Yes [1]
  
  IFWALD: 0=off [0]
  
  1=request distributions of sputtering and power along walls in print file.
  
  IFPRIM: 0=do not follow impurity neutrals 1=do [1]')
  
  IZWALL: Atomic number of wall
  
  IAEMIS: 0=Mol. reemission [0]
  
  1=Atomic Reemission
  
  -1/-2=INUTPG @ EATMR + AT./MOL
  
  KINDPR: Print Switch 0=minimum [0]
  
  TWALL: Vessel Wall Temperature (C) [300]
  
  ZESCUT: Gap polygons above or equal are wall [INF]
  
  JXLM: Knot for projection beyond X-pt. (Use default) [0]
  
  JXRM: Knot for projection beyond X-pt. (Use default) [0]
  
  XC1: Point for projection beyond X-pt. (Use default) [RPX]
  
  YC1: Point for projection beyond X-pt. (Use default) [ZPX]
  
  IALB: Albedo condition 0=wall 1=albedo 2=void for P.P. void [0]
  
  LWALL: Use Actual Vessel As Wall (True/False) [T]
  
  LBUFLE: Use Baffle (Set False) (True/False) [T]
  
  LPWALL: Use Vessel Wall for private region (True/False) (forces LBUFLE=F) [F]
  
  LPSEG: Use Explicit Source segments around private void. (Set True!) (True/False) [F]
  
  IALBPG: Switch for turning specific wall segments into albedo regions - rely on pump files - [off]
  
  ALBPG: Switch for turning specific wall segments into albedo regions - rely on pump files - [off]
  
  ALBEPG: Switch for turning specific wall segments into albedo regions - rely on pump files - [off]
  
  ALATO: Switch for turning specific wall segments into albedo regions - rely on pump files - [off]
  
  LNWESC: Use New Escape Figure Method (Use T) (True/False) [T]
  
  EATMD: Energy (eV) of neutrals re-emitted as atoms. 0.0=Franck-Condon [0.025]
  
  MCX: Model for energy after CX - use default [0]
  
  NTSPUT: Turn on neutral sputtering of impurities? 0=off 1=on [1]
  
  IHOR: Switch for multi-group vel. distributions. 0=off 1=on [0]
  
  DECIMA: Decimation probability - leave as is - [0.7]
  
  MODATM: Model for atomic CX Losses [1]
  
  NCOLP: Max # of collisons before R.R. analog game (0 to 100,00) (use default) [0]
  
  ISEHHE: Model for elastic scattering. 0=none 1to7=diff. comb. of H,HZ,HE (use 0) [0]
  
  RNLITE: Reflected fraction of light impurity (?) [1]
  
  EWLITE: Energy (eV) of reflected light impurity (?) [0]
  
  INUTPG: Regions to be set as recyclers (use default) [All Wall]
  
  EATMR: Enrgy (eV) of forced reflected neutrals [EATMD]
  
  TDIV: Divertor wall temperature (C) [TWALL]
  
  LPUMP: Switch on pump - in pump structure - (T/F) [T]
  
  INPUMP: Channel for reading pumpfile - SET=18!! - [LDUMIO]
  
  FPUMP: Pump structure file name default = none ['' '']
  
  ALBPMP: Albedo for pump <0=use pump value [-1e30]
  
  PSEMPO: Transparency of Outer SOL DIV <0= use pump file [-1e30]
  
  PSEMPT: Transparency of TARGET DIV <0=use pump file [-1e30]
  
  PSEMPI: Transparency of Inner SOL DIV <0= use pump file [-1e30]
  
  PSEMPO: Transparency of CHEVRON <0= use pump file [-1e30]
  
  IPSEMP: Define wall regions to be semi-transparent - Rely on pump files - [none]
  
  PSEMP: Define wall regions to be semi-transparent- Rely on pump files - [0]
  
  PSEMPB: Transparency of baffle - 1e30 = do not use [1e30]')
  
  IPVOID: Flag for treatment of pump void walls - use -1! - [1]
  
  IVIEW: 0=std. Nimbus geom. map -1=user defined window [0]
  
  VIEW: (Rmin,Zmin,Rlen, Zlen) of user defined window GEOM from GRID2D - [GEOM]
  
  ITRIM: 0=no TRIM files 1=use TRIM files -use 1!- [0]
  
  FTRIM: TRIM file prefix - set in /NIMBIN/ - [CTRIMF]
  
  LFULL: Full setup for NIMBUS at every call (T/F) - use T - [T]
  
  ITARHZ: Switches (2) to determine when to use horizontal escape figure - use default - [2*MX]
  
  ICHKP: 0=Stop on polygon problems in NIMBUS 1=warn - use default - [1]
  
  LTIME: Time dependent Monte Carlo (T/F) - use F - [F]
  
  TWIDTH: Time slice width - only for time-dependent mode - [0]
  
  TWDMIN: Minimum time slice width - only for time dependent mode - [0]
  
  AYIZ: Enhanced Yield from H: Y''=AY+B [1]
  
  BYIZ: Enhanced Yield from H: Y''=AY+B [0]
  
  ICORRN: Random number correlation flag - use default - [UNDEF]
  
  ALBLK: Albedo of pump leaks <0=use pump file - [-1e30]
  
  DCUTCX: Maximum density for CX [1e30]
  
  TCUTCX: Minimum temperature for CX [-1e30]
  
  ITAU: Model for flux estimator <=1=estimated dist. >1=dist.- use default - [1]
  
  IYCHEM: Model for chemical sputtering - 0=off 1+=model selected - [0]
  
  EYCHEM: Energy (eV) of chemically sputtered C ''atom'' - [0.0]
  
  IDBHST: Number of histories to store trajectories - use default - [0]
  
  XGAUGE: R position - to override gauge location in pump file - [1e30]
  
  YGAUGE: Z position - to override gauge location in pump file - [1e30]
  
  RGAUGE: Radius - to override gauge location in pump file - [1e30]
  
  LGAUGE: Label - 'G' for pump gauge (leave as is) - 'K' for vessel - [''G'',''K'']
  
  MODEZR: Model for impurity ion recycling - use default - [1]
  
  ISPOFF: Switch off recycling in specifed macrozones - [0]
  
  MIMP: # imp switch - being developed - use 0 for now - [0]
  
  LRS: Override leak recycling segments in pump file : '<ID> X1 Y1 X2 Y2' - [60*'' '']
  
  GAP: Override Gap segments in pump file : '<ID> X1 Y1 X2 Y2 T' (transparency) - [60*'' '']

I Tags
------
.. _I01:
I01 : Injection
  **Injection Option * **: Disregarded when NEUT control switch not 0
  
  **Injection Option 1**: Inject ions at given (r,z) with given v0
  
  **Injection Option 2**: Inject uniformly on a given ring (`I18`_) between INJ1*SMAX (`I19`_) and INJ2*SMAX (`I20`_) relative to both plates. Given V0.
  
  **Injection Option 3**: Inject uniformly on a given ring (`I18`_) between INJ1*SMAX (`I19`_) and INJ2*SMAX (`I20`_). Given V0.
  
  **Injection Option 4**: Neutral impurity ionization profiles taken from a nimbus/pin run are used to generate a probability map for ion injection. The initial ion energy is taken from the nimbus/pin results.
  
  **Injection Option 5**: Inject uniformly on a given ring (`I18`_) between INJ1*SMAX (`I19`_) and INJ2*SMAX (`I20`_) relative to both plates.
  
    Initial velocity is calculated from:

      .. math::
    
        v_{init} = r_g * v_0'
  
    where 

      .. math::

        r_g = \sqrt{-2*ln(x_1)} \times cos(2 \pi x_2)
  
    x\ :sub:`1`,x\ :sub:`2` are uniform on [0, 1] and the value for v\ :sub:`0` is given.
  
  **Injection opt 6**: Inject uniformly on a given ring (`I18`_) between INJ1*SMAX (`I19`_) and INJ2*SMAX (`I20`_).
  
  Initial velocity is calculated from:
  
    .. math::
    
        v_{init} = r_g * v_0'
  
  where 

      .. math::

        r_g = \sqrt{-2*ln(x_1)} \times cos(2 \pi x_2)
  
    x\ :sub:`1`,x\ :sub:`2` are uniform on [0, 1] and the value for v\ :sub:`0` is given.
  
  **Injection Option 7**: Based on FLUID CODE results.
  
    Neutral impurity ionization source profiles taken from a FLUID CODE run are used to generate a probability map for ion injection. The initial ion energy is ??? also taken from the FLUID CODE data.

  **Injection Option 8**: Needs to be documented.

  **Injection Option 9**: Needs to be documented.

  **Injection Option 10**: Needs to be documented.

  **Injection Option 11**: Needs to be documented.

  **Injection Option 12**: Needs to be documented.

  **Injection Option 13**: Needs to be documented. 

  **Injection Option 14**: Needs to be documented.

.. _I02:
I02 : First diffusion
  **First diffusion 0**: Instant

  **First diffusion 1**: After randomly generated time interval

  **First diffusion 2**: After ion has existed for time = Tau parallel

.. _I03:
I03 : Control switch
  **Control switch 0**: NEUT on: follow atoms to ionization positions

  **Control switch 1**: NEUT off: inject ions as "initial state" option

.. _I04:
I04 : Self Sputtering Option
  **Self-Sputter Opt 0**: OFF. Self-sputtering cascade does not occur.

  **Self-Sputter Opt 1**: ON. Self-sputtering cascades are followed as usual based on the maximum number of generations of self-sputtering allowed to be followed, the minimum specified threshold yield allowed, and the calculated actual yield based on the particle impact energy for self-sputtering.

  **Self-Sputter Opt 2**: ON. Self-sputtering cascades are followed as usual based on the maximum number of generations of self-sputtering allowed to be followed and the minimum threshold yield. The yield for each segment is set at a separately specified fixed value. The energy of emitted sputter fragments is also fixed. This option was implemented to allow for modelling of prompt re-emission of impurity species ions (particularly carbon) as they strike a target surface through mechanisms other than physical sputtering

.. _I05:
I05 : Initial Ion Velocity Option
  **Initial Ion Vel 0**: v\ :sub:`n` = 0.0

  **Initial Ion Vel 1**: :math:`\pm 0.5 v_n = 0.0` along S

  **Initial Ion Vel 2**: v\ :sub:`n` = 0.0 along S away from Target

  **Initial Ion Vel 3**: :math:`\pm \sqrt{\$} \times v_n` along S, $ in (0,1)

.. _I06:
I06 : Follow Recombined Impurity Neutral Option
  **Follow Recombined Neutrals Opt 0**: OFF. Recombined Impurity neutrals are not-followed.

  **Follow Recombined Neutrals Opt 1**: ON. Recombined Impurity neutrals are followed in the neutral state until re-ionization or removal from the system by some other mechanism. The re-ionized neutral is then followed as an ion. The recombined neutral is given the poloidal plane component of a projected 3D isotropic velocity calculated based on the temperature of the recombining ion. Unless over-ridden by the neutral velocity type option.

.. _I07:
I07 : Prompt Deposition Option
  **Ion Prompt Dep. Opt 0**: OFF. Ion prompt deposition does not occur.

  **Ion Prompt Dep. Opt 1**: ON. Ion prompt deposition is allowed to occur. If an ion is created within a larmor radius of the target surface OR if it is created within the Magnetic Pre-Sheath thickness of the asociated target segment then it is assumed to promptly redeposit on the target surface. The ion's impact energy for self-sputtering purposes is calculated based on the ion's creation energy and it's precise location of ionization relative to the Magnetic Pre-Sheath. The impact energy is reduced appropriately for ionization occurring within the MPS. 

  **Ion Prompt Dep. Opt 2**: Needs to be documented.

  **Ion Prompt Dep. Opt 3**: Needs to be documented.

  **Ion Prompt Dep. Opt 4**: Needs to be documented.

.. _I08:
I08 : Target Mirror Option
  **Target Mirror Opt 0**: OFF. Target is treated normally.

    This option was implemented for testing purposes and to look at some specific transport physics problems. It should normally be turned OFF.

  **Target Mirror Opt 1**: ON. The target reflects all ions striking it. Both position and velocity are reflected. If the particle at position A makes a parallel step of length B that would carry it into the target located at the position X=0 (B>A and towards the target) then the ion's velocity is changed in sign and the position of the particle is adjusted to "|B-A|". Self-sputtering does not occur and particles can not be removed from the system at the target.

  **Target Mirror Opt 2**: ON. The target reflects all ions striking it. Only the position and not the velocity is reflected. If the particle at position A makes a parallel step of length B that would carry it into the target located at the position X=0 (B>A and towards the target) then the ion's velocity is unchanged and the position of the particle is adjusted to "|B-A|". Self-sputtering does not occur and particles can not be removed from the system at the target. 

.. _I09:
I09 : Ion Periphery Option
  **Ion Periphery Opt 0**: Hard wall - ions are immediately removed at wall impact.

  **Ion Periphery Opt 1**: Reflecting wall - ions are always reflected at wall impact. The only loss mechanisms are recombination and target plate impact.

  **Ion Periphery Opt 2**: No wall - ions are allowed to cross-field diffuse indefinitely. However they are always associated with the outermost ring of the grid system. Thus plots using this option could be invalid because this effectively increases the volume of the outermost ring - thus influencing density calculations.

  **Ion Periphery Opt 3**: Far periphery model at wall (`I21`_). The ion enters a region where there are three loss mechanisms invoked and governed by three parameters. First, the ion can diffuse (`I23`_) back to the plasma - it enters at the point where it left. Second, it can be lost to the far periphery target region with a specified characteristic loss time (`I22`_). Finally, it can diffuse to the "wall" which is specified to be a constant distance from the outermost ring of the grid. In addition, the diffusion coefficient in the far periphery can be specified independently of the value in the SOL and main plasma.

  **Ion Periphery Opt 4**: Special debugging wall option. The main wall is treated as a hard surface for particle loss. All ions striking the main wall are removed from the system. The private plasma wall is treated as a mirror. All particles striking the private plasma wall are reflected back into the plasma. The purpose of this option is to facilitate testing of the various Dperp options. If the ony sink is the outside wall then it should be possible to reach a condition where the density is constant inboard of the location of ion injection - if the Dperp and parallel transport are spatially uniform. This option is used in conjunction with the mirror target option.

.. _I10:
I10 : Periphery Recycle Option
  **FP Recycle Option 0**: Ions are lost and removed from the system upon far periphery target or wall impact.

  **FP Recycle Option 1**: Ions lost to the far periphery target and walls are relaunched from the edge of the nearest target.

.. _I11:
I11 : Z effective (self) - Zeff
  This item is only required for collision option 2, otherwise it is ignored. Typically 1.0, 2.0, 4.0, etc. See note 103.

.. _I12:
I12 : Initial ionization state of impurity ions
  This quantity must be specified for non-NEUT cases where ions are injected straight into the plasma. Normally it will be 1, so that singly ionized ions are injected, but it can be greater than 1 if desired. For NEUT cases, this value is set to 1 internally. 

.. _I13:
I13 : Collision Enhancement Factor - Zenh
  This item is used in every case and should normally be set to 1.0. It is a modifying factor used in calculating the collision times, heating times and slowing times and can be set to, for example, 2.0, 4.0, 8.0 to obtain modified plasma characteristics. Defined in Note 121.

.. _I14:
I14 : Set Ti = max(Ti,Tb) when reaching state (0 off)
  This item is rarely used and should normally be 0 for "off". Entering a value higher than 0 indicates an ionization state which, when reached, results in an ion being instantaneously heated (or cooled) to the local plasma background temperature at that point. For example, if this item is set to 2, then any ion ionizing to 2+ or recombining to 2+ has the plasma temperature at its current position assigned to it. Described in Note 121.

.. _I15:
I15 : Maximum ionization state
  This item gives the maximum charge state of interest. There are both time and storage advantages in restricting the number of charge states with this parameter. It should be used in conjunction with the "ionization option" above. For options 0,3 and 4, ions are allowed to ionize beyond the maximum charge state given here, but they are then recorded as having "ionized beyond limit" and ignored. For options 1,2,5 and 6, ions are just not allowed to ionize beyond the maximum state given here, so any ion reaching this state remains until it is removed. Note that if option 5 is used (which allows electron-ion recombination), ions which are "stuck" in the highest allowed ionization state are still able to recombine to lower states. Any singly-ionized ions which recombine are recorded and ignored. One special use of this item is in doing a neutral analysis only - if the maximum charge state is set to 0 then any ions created are immediately recorded as "ionized beyond limit", without being tracked. If the maximum state is given as the atomic number of the ion, then of course there are no restrictions on ionization or recombination.

.. _I16:
I16 : Stop following ions reaching Main Plasma
  **0**: No

  **1**: Yes

  This item has been used for testing SOL processes without using too much CPU time. Any ions reaching the main plasma are noted and their trajectories stopped.

.. _I17:
I17 : Ion Loss Time
  This quantity defines an arbitrary ion removal function with the characteristic time defined by this value. The probability of an ion being eliminated on any specific time step is QTIM/LOSSTIME. A value of 0.0 entered for this will turn this function off. (This would be the equivalent of an infinite loss time.)

  NOTE: It is very important to keep in mind the nature of random number generators when using this and other similar functions in a Monte Carlo code or any code that depends on random numbers. As an example, assume a characteristic loss time of one second and QTIM = 10\ :sup:`-7` s was specified, then the probability of elimination on any time step would be 10\ :sup:`-7`. Unfortunately, the spectrum of many simple random number generators yield a maximum of ~1/32500 or at best - a probability division of 10\ :sup:`-5`. If the condition test range is 0 <= ran < prob, then usually the result is obtained with a probability of ~ 10\ :sup:`-5 since the value 0.0 is usually a part of the spectrum of these simple random number generators. This can skew the likelihood of low probability events and in addition can be very difficult to discover.

.. _I18:
I18 : Ring number for ion injection - INJIR
  Specifies the ring number for ions to be injected on for injection options 2, 3, 5 and 6 (`I01`_).

.. _I19:
I19 : Injection region - Lower Bound - INJ1
  Lower bound injection multiplier for injection options 2, 3, 5 and 6 (`I01`_). The ions will be injected between INJ1*SMAX and INJ2*SMAX symmetrically from both plates. Thus values of .15 and .2 will result in ion injection in the two regions 0.15 * SMAX to 0.2 * SMAX and 0.80 * SMAX to 0.85 * SMAX. For options 3 and 5, the particles will be injected between INJ1*SMAX and INJ2*SMAX, wherever that may be, including values of INJ1, INJ2 > 0.5 and < 1.0.

.. _I20:
I20 : Injection region - Upper Bound - INJ2
  Upper bound as described above.

.. _I21:
I21 : Far Periphery Width definition
  This defines the "width" in meters of the far periphery region model (`I09`_). Ions must diffuse cross-field this far in order to be considered as striking the walls. Different values for the OUTER and INNER halves of the main plasma may be specified. The other loss mechanisms from the far periphery region include diffusing back into the main plasma and far periphery loss to the target plates with a characteristic time specified by the following quantity.

  e.g. 'TN443 X-max for Far Periphery Region (Outer/Inner)' 0.1 0.1

.. _I22:
I22 : Far Periphery Target Loss - characteristic time
  This defines the probability, in an arbitrary sense, that an ion entering the far periphery will experience sufficient parallel diffusion or other effects to carry it to the target plates outside the normal strike zones. The loss time is entered in units of seconds. Different values may be specified for the OUTER and INNER halves of the main plasma edge as described above for the Far Periphery width. These ions do not at this time result in self-sputtering - but the totals of such particles are recorded for later analysis. These particles may be re-launched from the corners of the targets closest to their point of entry into the Far Peripheral Region. (See the input option TN998 - Far Periphery Recycle Option described earlier.)

.. _I23:
I23 : Far Periphery Diffusion Rate
  This specifies the diffusion coefficient (m\ :sup:`2`/s) to be used for the far periphery region. This can be used to vary the ion sink strength of the periphery region by increasing the rate of cross-field diffusion. It can be set to any value. A value less than zero results in the diffusion coefficient in the far periphery region being the same as that used for the rest of the plasma.

K Tags
------

N Tags
------
.. _N01:
N01 : Launch
  **Launch Option 0**: Distributed launch along target

  **Launch Option 1**: At given (R,Z) (`S08`_, `S09`_)

  **Launch Option 2**: Homogeneously along walls - different wall definitions are possible.

  **Launch Option 3**: Distributed launch along Target due to hydrogenic ion impact utilizing PIN/NIMBUS data.

  **Launch Option 4**: Distributed launch along "Wall" surfaces due to hydrogen ATOM impact utilizing PIN/NIMBUS data for wall fluxes and wall definitions. Walls are defined as all vessel surfaces including the target segments. Wall coordinates are passed to DIVIMP from NIMBUS.

  **Launch Option 5**: 2D Distributed Launch (`N07`_). Neutral particles are launched from a selection of cell centres on the grid. The probability of launch from a given cell centre is taken from an externally supplied distribution function. The only currently supported distribution function is the Carbon neutral recombination rate that has been imported from some UEDGE modelling of DIIID.

.. _N02:
N02 : Velocity/Angle Flag
  **Velocity/Angle Flag 0**: 

    :math:`\theta = \pm \arcsin (\xi),\ \ \xi \in (0,1)`

    :math:`v_{in} = \sqrt{\frac{2E_{bd}}{m_i (\sqrt{\xi^{-1}} - 1)}},\ \ \xi \in (0,1)`

  **Velocity/Angle Flag 1**: 

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`v_{in} = v_{mult} \sqrt{\frac{2E_{bd}}{m_i (\sqrt{\xi_1^{-1}} - 1)}},\ \ \xi_1 \lt \xi_{1,max}`

    :math:`\xi_{1,max}` is a limit on the random number that is derived according to the restriction that the sputtered neutral must have enough energy to overcome the surface binding energy, E\ :sub:`bd`. If the maximum energy that could be transfered to a neutral in a collision is E\ :sub:`max`, then the corresponding max on the random number is (*derivation lost to the ether, but likely not too hard to redo*):

    :math:`\xi_{1,max} = (1+E-{bd}/E_{max})^{-2}`

    :math:`v_{mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}` 

    :math:`\beta = \arcsin(\sqrt{\xi_2}),\ \ \xi_2 \in (0,1)`

    :math:`\phi = 2 \pi \xi_3,\ \ \xi_3 \in (0,1)`

  **Velocity/Angle Flag 2**: 

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`v_{in} = \sqrt{2T_g/m_i} \sqrt{| ln(1-\xi) |},\ \ \xi \in (0,1)`

    Gas temperature, :math:`T_g`, given in (`S06`_).

  **Velocity/Angle Flag 3**: 
  
    :math:`\theta = \pm \arcsin(\sqrt{\xi}),\ \ \xi \in (0,1)` 

    :math:`v_{in} = \sqrt{2E_{in}/m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 4**: 

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`v_{in} = v_{mult} \sqrt{\frac{2E_{bd}}{m_i (\sqrt{\xi_1^{-1}} - 1)}},\ \ \xi_1 \lt \xi_{1,max}`

    See Velocity/Angle Flag 1 for a description of :math:`\xi_{1,max}`.

    :math:`v_{mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}` 

    :math:`\beta = \arcsin(\sqrt{\xi_2}),\ \ \xi_2 \in (0,1)`

    :math:`\phi = 2 \pi \xi_3,\ \ \xi_3 \in (0,1)`

    Velocities are limited by the Emax factor given in (`D29`_).

  **Velocity/Angle Flag 5**: 

    :math:`\theta = \pm \arcsin(\sqrt{\xi_1}),\ \ \xi_1 \in (0,1)` 

    :math:`v_{in} = \sqrt{\frac{2E_{bd}}{m_i (\sqrt{\xi_2^{-1}} - 1)}},\ \ \xi_2 \lt \xi_{2,max}`

    Velocities are limited by the Emax factor given in (`D29`_).

  **Velocity/Angle Flag 6**: 

    :math:`\theta = 0`

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 7**: 

    :math:`\theta = \pm \arccos((1-\xi)^{1/3}),\ \ \xi \in (0,1)`   "Free Jet" 

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 8**: 

    :math:`\theta = 2 \pi \xi,\ \ \xi \in (0,1)`   "Isotropic"

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 9**: 

    :math:`\theta = \pm \arcsin(\sqrt{\xi}),\ \ \xi \in (0,1)` 

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`
 
    Two given values are used alternately for :math:`E_{in}` (`S06`_, `S07`_).

  **Velocity/Angle Flag 10**: 

    :math:`\beta = \arccos((1-\xi_1)^{1/3}),\ \ \xi_1 \in (0,1)`   "3D Free Jet"

    :math:`\psi = 2 \pi \xi_2,\ \ \xi_2 \in (0,1)`

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 11**: 

    :math:`\beta = \arccos((1-\xi_1)^{1/3}),\ \ \xi_1 \in (0,1)`   "2.5D Free Jet"

    P0 assigned randomly in range (-Pmax, Pmax). *Unclear what P0 is.*

    :math:`v_{in} = \sqrt{2E_{in} / m_i}`

    Initial energy, :math:`E_{in}`, given in (`S06`_).

  **Velocity/Angle Flag 12**: 
    Emission at a constant energy specified by the input quantity :math:`E_{in}` (`S06`_, CTEM1) into a :math:`\cos^N` distribution where N is specified with `N18`_.

  **Velocity/Angle Flag 13**: 
    Emission at a temperature (:math:`T_g`) into a :math:`\cos^N` distribution. :math:`T_g` is specified by the same input quantity (`S06`_, CTEM1) as in Option 12. N is also a specified input as in V/A flag 12 (`N18`_).

  **Velocity/Angle Flag 14**: 

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`v_{in} = v_{1,mult} v_{2,mult} \sqrt{\frac{2E_{targ}}{m_i}}`

    :math:`v_{1,mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}` 

    :math:`\beta = \arcsin(\sqrt{\xi_1}),\ \ \xi_1 \in (0,1)`

    :math:`\phi = 2 \pi \xi_2,\ \ \xi_2 \in (0,1)`

    :math:`v_{2,mult}` is specified with `N19`_.

  **Velocity/Angle Flag 15**: 

    :math:`\theta = 2 \pi \xi,\ \ \xi \in (0,1)\ \ \ `   "Isotropic"

    :math:`v_{in} = v_{mult} \sqrt{T_i/m_i}`

    :math:`T_i` is the local ion temperature, and :math:`v_{mult}` is specified with `N19`_.

  **Velocity/Angle Flag 16**:

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`Y(E) = \frac{E}{(E+E_{bd})^3} (1 - \sqrt{\frac{E+E_{bd}}{G(1-G)E_{imp}}})`

    :math:`G = \frac{4m_i m_b}{(m_i + m_b)^2}`

    :math:`v_{in} = v_{mult} \sqrt{2*E / m_i}`

    :math:`v_{mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}` 

    :math:`\beta = \arcsin(\sqrt{\xi_1}),\ \ \xi_1 \in (0,1)`

    :math:`\phi = 2 \pi \xi_2,\ \ \xi_2 \in (0,1)`

    The energy, E, is selected randomly from Y(E). *Presumably this means using Y(E) as a probability distribution function.*

  **Velocity/Angle Flag 17**:

    :math:`\theta = 0\ \ \ `    Normal to surface

    :math:`v_{in} = \sqrt{\frac{2E_{bd}}{m_i (\sqrt{\xi^{-1}} - 1)}},\ \ \xi \in (0,1)`

  **Velocity/Angle Flag 18**:
    Same as Option 16 just :math:`\theta = 0`.

    :math:`\theta = 0\ \ \ `    Normal to surface

    :math:`Y(E) = \frac{E}{(E+E_{bd})^3} (1 - \sqrt{\frac{E+E_{bd}}{G(1-G)E_{imp}}})`

    :math:`G = \frac{4m_i m_b}{(m_i + m_b)^2}`

    :math:`v_{in} = v_{mult} \sqrt{2*E / m_i}`

    :math:`v_{mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}` 

    :math:`\beta = \arcsin(\sqrt{\xi_1}),\ \ \xi_1 \in (0,1)`

    :math:`\phi = 2 \pi \xi_2,\ \ \xi_2 \in (0,1)`

  **Velocity/Angle Flag 19**: 
    3D Isotropic

    :math:`\theta = \arctan(\tan(\beta)\cos(\phi))`

    :math:`v_{in} = v_{mult} \sqrt{2*E_{in} / m_i}`

    :math:`v_{mult} = \sqrt{| \cos(\beta)^2 + \sin(\beta)^2 cos(\phi)^2 |}`  

    :math:`\beta = 2 \pi \xi_1,\ \ \xi_1 \in (0,1)`

    :math:`\phi = 2 \pi \xi_2 - \pi,\ \ \xi_2 \in (0,1)`

    :math:`E_{in}` specified in `S06`_.

  **Velocity/Angle Flag 20**:
    Data from ERO for 2D particle source

.. _N03:
N03 : Supplementary Launch Option
  These are identical to the possible launch options and apply to any supplementary neutral launches.

  **Sup. Launch opt -1**: Set value to the same as the primary launch option.

.. _N04:
N04 : Supplementary Velocity/Angle Flag
  These are identical to the possible V/A flag options and apply to any supplementary neutral launches.

  **Sup. V/A flag opt -1**: Set value to the same as the primary V/A flag option.

.. _N05:
N05 : Initial Neutral Velocity/Angle Flag
  These are identical to the possible V/A flag options and apply to any the initial group of neutral launches only. This allows such for simulations of a free-space neutral pellet ablation followed by regular self-sputtering. Each of which has quite different characteristic Velocity/Angle source distributions.

  **Init. Neut. V/A opt -1**: Set value to the same as the primary V/A flag option.

.. _N06:
N06 : Extra 2D Neutral Launch Option
  **2D Neut. Launch Opt 0**: OFF. No 2D neutral source is used.

  **2D Neut. Launch Opt 1**: ON. A 2D source of impurity neutrals - equivalent to launch option 5 is launched in addition to all other specified impurity sources. Each source is weighted according to its relative production strength. This option is required in order to allow for both wall and target particle sources as well as 2D distributed impurity neutral sources. This could be either modelling recombined impurities or modelling a puff of some description in addition to regular target production.

.. _N07:
N07 : 2D Neutral Launch - Velocity/Angle flag option
  This takes the same options as those described in the Velocity/Angle Flag above. This Velocity/Angle flag is applied to any particles launched using the 2D Neutral Launch described in the previous option. 

.. _N08:
N08 : Sputter Option
  **Sputter Option 0**: Sputtering by background ions (Mb,Zb) only

    :math:`E_{imp} = T_b (2 + 3Z_b)`

    :math:`E_{max} = E_{imp} \Gamma (1 - \Gamma) - E_{bd}`

    :math:`\Gamma` is from...

  **Sputter Option 1**: Sputtering by specified ion type only

    :math:`E_{imp} = T_b (2 + 3Z_b)`

    :math`E_{max} = E_{imp},\ \ Z_{imp}\ `given

  **Sputter Option 3**: Initial sputtering by background ions only.
    Identical to Option 1. See comment at end of this entry.  

    :math:`E_{imp} = T_b (2 + 3Z_b)`

    :math:`E_{max} = E_{imp} \Gamma (1 - \Gamma) - E_{bd}`

  **Sputter Option 4**: Initial sputtering by background ions only.

    :math:`E_{imp} = T_b (2 + 3Z_b)`

    :math:`E_{max} = E_{imp} \Gamma (1 - \Gamma) - E_{bd}`

    The maximum energy of the distribution of self-sputtered particles resulting from the use of this option is multiplied by the EMAX-FACTOR described later (`D29`_, CEMAXF). Velocity/Angle flag 1 will be used for any self-sputtered particles.

  **Sputter Option 5**: Initial Chemical Sputtering Source. The formula used to calculate the chemical sputtering is specified using the Chemical Sputtering Source Option defined above. The surface temperatures of the walls and the target can be specified separately and are defined later in this document.

  **Sputter Option 6**: Initial COMBINED physical and chemical sputtering. Two groups of atoms are launched. The first is PHYSICALLY SPUTTERED using the Velocity/Angle flag specified by the Initial Neutral V/A flag option. The second group is CHEMICALLY SPUTTERED and uses a Velocity/Angle flag of 3 with the characteristic energy specified in the input. The ratio of particles launched through each mechanism is proportional to the total (BG FLUX) * (YIELD) for each sputter source.

  Note: In previous versions of DIVIMP - sputter options of 3 or more were used to indicate that self-sputtering would be active. This is no longer the case. A separate switch has been added to turn self-sputtering On or Off. Thus there is some duplication of definition in the various sputter options.

.. _N09:
N09 : Secondary Sputter Option (TN1209)
  These are identical to the possible Sputter options listed above and apply to any supplementary neutral launches.

  **Sup. Sputter opt -1**: Set value to the same as the primary sputter option.

.. _N10:
N10 : Normal
  **Normal option 0**: Measure :math:`\theta` from surface normal

  **Normal option 1**: Measure :math:`\theta` from T degrees to X=0, T given in `N15`_.

    Apply to primary neutrals only.

  **Normal option 2**: Measure :math:`\theta` from T degrees to X=0, T given in `N15`_.

    Apply to primary and self-sputtered neutrals.

.. _N11:
N11 : NEUT spreading
  **NEUT spreading 0**: OFF (Launch at meshpoints only)

  **NEUT spreading 1**: NOT SUPPORTED ... (ON -Launch with variation in Y0)

.. _N12:
N12 : Impurity Neutral Velocity Type Option
  **Velocity Type Opt 0**: OFF. Impurity neutrals will have a constant speed from creation until final removal.

  **Velocity Type Opt 1**: ON. Impurity neutral speed will change as the particle moves across the grid. The impurity neutral is assigned a speed based on the local plasma ion temperature. This speed changes as the neutral enters a new cell on the grid and is adjusted to match the local temperature.

  **Velocity Type Opt 2**: ON. Impurity neutral speed may change as the particle moves. The neutral will be assigned a new speed based on the local plasma ion temperature and the impurity mass whenever a Momentum Transfer Collision is calculated to take place.

.. _N13:
N13 : Neutral Reflection Option
  **Reflection Opt 0**: Off - Impurity neutrals striking the walls are recorded. They do not cause self-sputtering nor are they reflected.

  **Reflection Opt 1**: On - Impurity neutrals striking the wall are reflected specularly retaining the same energy as they had before impact. They do not cause sputtering at wall impact. This may be a poor approximation for most carbon-wall collisions. It might not be too bad for inert gas collisions.

  **Reflection Opt 2**: On - Impurity neutrals striking the wall are reflected with a cosine angular distribution retaining the same energy as they had before impact. They do not cause sputtering at wall impact. This may be poor approximation for most carbon-wall collisions. It might not be too bad for inert gas collisions.

.. _N14:
N14 : Impurity Neutral Momentum Transfer Collision Option
  **Momentum Transfer Collision Opt 0**: OFF. Neutrals travel in straight lines from creation until ionization.

  **Momentum Transfer Collision Opt 1**: ON. Neutrals will undergo 90 degree changes in flight path based on the probability calculated for the occurrence of a momentum transfer collision. The collision frequency is given by the following expression.

    :math:`\nu_{mfc} = \frac{16 m_b}{3(m_z + m_b)} (K_{elighi} n_b + K_{elighg} n_b0)`

    The variables :math:`K_{elighi}` and :math:`K_{elighg}` are given in `D09`_ and `D10`_, respectively. *This option could use a link or something explaining where it comes from.*

.. _N15:
N15 : Measure theta from T degrees for launch
  Note on definition of :math:`\theta`: in DIVIMP, :math:`\theta = 0.0` in the positive R direction, and is measured anti-clockwise. Hence an angle of -90.0 degrees might be used for launch uniformly about the -Z direction.

.. _N16:
N16 : Wall Launch Segment Probability Multipliers
  This applies to the Wall Launch option (`N01`_, option 2). The wall is composed (presently) of a set of line segments joining the points defining the outermost SOL ring, both plates, and the outermost trap ring. Depending on the geometry, there are approximately 60 of these segments. For the wall launch, all of the target segments automatically have a probability multiplier of zero - thus preventing wall launches from occurring at the targets. (The previous statement may no longer be true in some DIVIMP versions, i.e. the targets will be included as valid sources of wall-launched particles, so it is desirable to explicitly set the launch probability modifiers for these segments to zero.) In addition, this option can be used to further restrict and/or limit the probability of neutral launch from any particular wall segment or range of wall segments. Initially, the probability of launch from a specific segment is proportional to the length of that segment. This function will modify those probabilities.

    .. code-block::

      ' ' ' TN487 Launch Probability modifiers for each '
      ' TN487 wall segment range #1 #2 mod. :- ' 2
         1   28   0.0
        53   60   0.0

  The example above prevents neutral launch from the first half of the main wall region and from the walls of the trapped plasma region. The probability modifiers could as easily have been set to 0.1 to reduce but not eliminate the wall sources in these regions or to a value greater than 1.0 to enhance the probability of launch from these regions.

  This option is also used, in combination with the following option, to specify the proportional probability of launch from each wall segment. The lengths of the wall segments are ignored and the values in this table substituted and normalized to 1.0 to generate a wall launch distribution. Finally, if the wall launch probability data is being generated directly from the results of a PIN run then this table can be used to modify those numbers via a mechanism the same as that above.

.. _N17:
N17 : Absolute Wall Probabilities
  **Option 0**: Off - treat wall launch data as multipliers of initial wall segment launch probability which is calculated based on the relative size of the wall segments.

  **Option 1**: On - treat wall launch probabilities as actual segment probability weightings.

  When this is set to 0, the values in the table above are taken as relative multipliers of wall launch probabilities obtained through another method. If it is set to 1, then the values entered in the table above are used as the actual weight for launch from that specific wall segment with the sum total of all weights equal to one and the individual weights renormalized using that sum.

  **Option 2**: On - Wall Segment Launch probabilty is loaded from PIN and multiplied by the probability modifiers.

  **Option 3**: On - Wall Segment Launch probabilty is loaded from PIN and is NOT multiplied by the probability modifiers.

.. _N18:
N18 : Power of Cosine release distribution (V/A flag 12,13)
  This specifies the exponential power of the cosine distribution used to inject particles in V/A flag (`N02`_) options 12 and 13.

.. _N19:
N19 : Velocity Multiplier for Velocity/Angle Flag 14 and 15
  This number is used to multiply the velocity of all initial impurity neutral launches which use Velocity/Angle flags (`N02`_) 14 or 15. It does not apply to the velocity assigned to any recombined impurity neutrals.

.. _N20:
N20 : Velocity Multiplier for Recombined Ions
  This quantity multiplies the velocity assigned to recombined impurity neutrals (`I06`_).

.. _N21:
N21 : External Sputtering Flux Data Source
  **Option 0**: Geier file format for Argon

  **Option 1**: Import DIVIMP charge-resolved flux and energy data from a previous DIVIMP run

  Default: 1

O Tags
------

P Tags
------

.. _P01:
P01 : SOL
  **SOL option -1**: SOL1a (fl,fs) values given (N344)

    Background velocity and electric field are specified by formulae. See the sol.d6a source code module for a detailed description of the formulae.

  **SOL option 0**: SOL0

    Background velocity and electric field are set to zero everywhere.

  **SOL option 1**: SOL1

    Background velocity and electric field are specified by formulae. See the sol.d6a source code module for a detailed description of the formulae.

  **SOL option 2**: SOL2

    The velocity is set to Cs * (1-S/SMAX) while the electric field is set to -Te/SMAX.

  **SOL option 3**: SOL3

    Background velocity and electric field are specified by formulae. See the sol.d6a source code module for a detailed description of the formulae.

  **SOL option 4**: SOL4

    Background velocity and electric field are specified by formulae. See the sol.d6a source code module for a detailed description of the formulae.

  **SOL option 5**: SOL5

    Background velocity and electric field are held constant in the SOL with appropriate changes in sign for the different ends of the flux tube. The values assigned to the velocity and electric field are entered separately in the input file.

  **SOL option 6**: Efield = constant.

    Vb is set to three discrete constant values in three separate regions. The actual values assigned are controlled by the SOL option 6 and 7 parameters specified with `P15`_ - `P18`_.

  **SOL option 7**: Efield = constant.

    Vb is split into three linearly ramping regions with Vb going to zero at SMAX/2.0. Vb varies linearly across each region from the value specified at one boundary to the value specified for the next. The actual values assigned are controlled by the SOL option 6 and 7 parameters specified with `P15`_ - `P18`_.

  **SOL option 8**: Not used

  **SOL option 9**: SOL9 similar to SOL2 (N327)

    The velocity is set to the sound speed everywhere with appropriate sign for each half flux tube. The electric field is set to -Te/SMAX where Te is the temperature in the current cell and SMAX is the field line length from target to target.

  **SOL option 10**: SOL10 fRM, fRmin, fRmax, Kin, Kout and (fl,fs) values given in `P08`_ - `P14`_ (N353)

  **SOL option 12**: SOL12 - Pseudo self consistent

    Overrides Tgrad options. Electron and ion heat transport are equal. Te, Ti, Nb, Vh and E are calculated from the following equations. The value for K0 (`P20`_) is the same for both electrons and ions.

    *Need to type up equations here*

  **SOL option 13**: SOL13 - Pseudo self consistent Overrides Tgrad options. Electron and ion heat transport are independent. Te, Ti, Nb, Vh and E are calculated from the following equations. The value of K0 differs for electrons (`P21`_) and ions (`P21`_).

    *Need to type up equations here*

  **SOL option 14**: SOL14 - Pseudo self consistent Overrides Tgrad options in the SOL. Calculates Ti, Te, Nb, Vh and E from the following equations involving both heat conduction and convection.

    *Need to type up equations here*

  **SOL option 15**: SOL15 - Pseudo self consistent Overrides Tgrad options. Electron and ion heat transport are independent. Te, Ti, Nb, Vh and E are calculated from the following equations. The value of K0 differs for electrons (`P21`_) and ions (`P21`_).

    *Need to type up equations here*

  **SOL option 21**: Detached Plasma Prescription (see `R01`_ to `R12`_).

    The regions may be defined either in terms of SMAX or PMAX. This is controlled by the length switch. The following description uses SMAX as an example.

    Three regions:

      A: 0 < SMAX * L1

      B: L1 * SMAX < L2 * SMAX

      C: S > SMAX * L2

      A: Te at outer edge of A : Te0 * TeR1

      A: Ti at outer edge of A : Te0 * TiR1

      A: Ne at outer edge of A: Ne0 * NR1

      Te increases linearly in A to Te0*TR1 at SMAX*L1

      Ti increases linearly in A to Te0*TR1 at SMAX*L1

      Ne increases linearly in A to Ne0*NR1 at SMAX*L1

      Velocity in A: v(s)=N0V0/n(s)

      B: T=(T1**3.5+7/2K0*(Q0(s-L1) + 1/2(s-L1)^2*(Qrad/Lrad)))**(2/7)

      B: Radiated power in B (Qrad): Q0*QR

      C: T=(T2**3.5+7/2K0*(Qtot(s-L2)))**(2/7)

      Qtot = Q0 + Qrad

      B,C: N(s) = N1 * (T(s)/T1)**(-1)

      Velocity linearly->0 at SMAX*VR1

  **SOL Option 22**: SOL Option 22 is a multi-parameter OSM (Onion Skin Model) that uses a Runge-Kutta based solver to numerically evaluate the one-dimensional fluid equations. Starting from given target conditions and including a wide variety of effects, the option generates a complete background plasma solution for density, ion and electron temperature, parallel flow velocity and parallel electric field. The solutions represent what appears to be a reasonable approximation to the conditions found in the edge. As such, they should be a useful research tool in examining the behaviour of impurities as well as the hydrogenic species behaviour in the reactor and especially for comparing predicted observables to actual experimental results. It is important to verify that the solution generated, for a given set of target conditions, appear to be reasonable. The code can not definitively evaluate the legitimacy of a given background solution - this requires that the user look at the solution generated and evaluate it in the context of the physical situation being modelled. This SOL option has a large number of sub-options that are specified in a block at the bottom of the input file and are documented later in this manual. The equations that are solved are listed below - the solver can find solutions for Te=Ti or for Te and Ti evolving independently. Various switches and inputs controlling SOL 22 can be found in the `200 Tags`_ section.

    :math:`\frac{d}{ds} (\frac{5}{2} n(s) v(s) kT_e(s) - \kappa_{0e}T_e(s)^{5/2} \frac{dT_e(S)}{ds}) = -P_{rad}(s) - P_{helpi}(s) - P_{ei}(s)`

    :math:`\frac{d}{ds} (\frac{5}{2} n(s) v(s) kT_i(s) + \frac{1}{2} m n(s) v(s)^3  - \kappa_{0i}T_i(s)^{5/2} \frac{dT_i(S)}{ds}) = -P_{CX}(s) + P_{ei}(s)`

    :math:`\Gamma(s) = n(s) v(s) = n_0 v_0 + \int_0^s S(s')ds'`

    :math:`n(s) ((kT_e(s) + kT_i(s)) + mv(s)^2) = 2n_0(kT_e(s) + kT_i(s))`

    :math:`P_0 = (2kT_{i0} + 5kT_{e0}) n_0 v_0`

  **SOL option 23**: SOL option 23 solves the 1D fluid equations in a DIVIMP context. More deatil is provided in the SOL 23 specific documentation.

  **SOL Option 29**: This is an experimental SOL type where the plasma background within each cell is the average of the blobs counted within that cell. Blobs are launched from the separatrix as particles and are followed with characteristic decay times for the temperatures and densities within them. This is certainly too simple of a prescription for the SOL, so this solver is not recommended for use this time and is considered under development. Developed by Shawn Zamperini. 

  **SOL option 98**: Read data from DIVIMP generated background plasma file - this option must be used in combination with plasma decay option 98.

  **SOL option 99**: Read data from file - B2 or Edge2D depending on grid.

.. _P02:
P02 : Core Plasma Options
  **Core Option -1**: Ignore. This option will cause all of the core processing options to be bypassed. Any values set for the core plasma either initially or by routines other than the core plasma code section will be used for the core plasma. The values for the core will be either those specified in the INITPLASMA routine or some values applied through other options in the code.

  **Core Option 0**: Normal. This has been the standard DIVIMP option in the past. The temperature and density are constant on each core ring around it's length and increase by a specified step for each ring inside the separatrix ring. These quantities are specified by TebIn, TibIn and NbIn. The flow velocity is always zero in the core.

  **Core Option 1**: Core Conditions Specified (`Q37`_.). The quantities, Te, Ti and Nb are specified for each ring in the core separately. The values are constant along each ring.

  **Core Option 2**: Core Marfe Option. The values of Te, Ti, Ne, Vb at the X-point are specified for each ring (`Q37`_.). These ramp up linearly along the core ring to the standard DIVIMP values (specified as in option zero) at the inner and outer mid-planes.

  **Core Option 3**: Core Marfe Option. The values of Te, Ti,---, and Vb at the X-point are specified (`Q37`_, Note that the density is NOT specified in this option). These ramp up linearly from the X-point along the core ring to the standard DIVIMP values for that core ring (specified as in option zero) at the inner and outer mid-planes. The density along the ring is calculated by applying pressure conservation. A number must be supplied for the density in the input but it is ignored.

  **Core Option 4**: Core Marfe Option. The values of Te, Ti,---, and Vb near the X-point are specified. (`Q37`_, Note that the density is NOT specified in this option). In addition, location factors (`P56`_-`P59`_) are also specified for the temperature and velocity. The temperature and velocity ramp up linearly from the first location along the core ring to the standard DIVIMP values for that core ring (specified as in option zero) at the second location specified. The density along the ring is calculated by applying pressure conservation. A number must be supplied for the density in the input but it is ignored.

    e.g.

    V -> 0 for S < Vf1 * SMAX

    V -> Vb specified above for S = Vf1 * SMAX

    V ramps linearly from Vb down to zero

    for Vf1 * SMAX < S < Vf2 * SMAX

    V -> 0 for S > Vf2 * SMAX

    Te,Ti -> Tex,Tix for S =< Tf1 * SMAX

    Te,Ti ramp up linearly to the standard values

    for Tf1 * SMAX < S < Tf2 * SMAX

    Te,Ti -> Standard for S> Tf2 * SMAX

  **Core Option 5**: This option is exactly the same as option 4 except that instead of using the standard option for calculating the base conditions on each core ring - this information is expected to be entered on a ring by ring basis in the data input block used to specify the plasma background for the SOL (`Q34`_). See plasma decay option 7 (`P03`_). Use of this option does NOT require that plasma decay option 7 be specified - only that the data be entered in the appropriate portion of the input file. This allows alternative methods of specifying target conditions to be used while reading in the core conditions and simulating a Marfe. 

.. _P03:
P03 : Plasma Decay
  **Plasma decay 0**: Standard method (N309)

  **Plasma decay 1**: Exponential decay outboard using the distance along the target from the separatrix strike point as the distance for the decay. See `Q06`_-`Q11`_, `Q16`_-`Q21`_ for more information. 

  **Plasma decay 2**: Temperature and density taken from input data for rings in the SOL (`Q34`_, `Q36`_).

  **Plasma decay 3**: Temperature and density taken from input data for rings in the SOL. Inner (`Q34`_) and Outer (`Q36`_) plates may differ.

  **Plasma decay 4**: Temperature and density taken from input data for rings in the SOL and Trap. Inner (`Q34`_) and Outer (`Q36`_) plates may differ.

  **Plasma decay 5**: Exponential Decay Outboard using the straight-line distance of the target point from the separatrix strike point as the distance for the decay. See `Q06`_-`Q11`_, `Q16`_-`Q21`_ for more information. 

  **Plasma decay 6**: Exponential Decay Outboard as in option 5 except that different exponential decay factors are specified for the main SOL and Private Plasma target regions. See `Q06`_-`Q11`_, `Q16`_-`Q21`_ for more information. 

  **Plasma decay 7**: Temperature and density taken from input data for rings in the main SOL, Private Plasma and the Core. Inner (`Q34`_) and Outer (`Q36`_) plates may differ.

  **Plasma decay 90**: Compound Background Plasma Option. The background plasma is created by using different options for different half-rings throughout the plasma. The specific input lines are described below (`P04`_). The plasma solution can be iterated through several PIN iterations. Option 90 uses an EDGE2D or other fluid code solution as the basis and overlays the other specified pieces of the solution over the different segments of the grid.

  **Plasma decay 91**: Compound Background Plasma Option. This option also creates a piece wise assembled background plasma (`P04`_) by allowing a great deal of flexibility in what options are allowed for each half-ring. This option uses the given input values for the regular SOL options combined with an assumed plasma decay option of 4 to generate a base background plasma. The specified pieces are then overlaid on top of this.

  **Plasma decay 98**: Read data from a DIVIMP formatted plasma transfer file written by a previous DIVIMP run. A plasma file can be created on any DIVIMP run by setting the print option equal to 10.

  **Plasma decay 99**: From data file.

  Note: Most Plasma Decay options (e.g. 3 and 4) may be used to assign uniform values to the entire plasma region and not just at the targets. Other temperature gradient options may then take these values that have been assigned to the entire plasma and change them. For example, SOL option 13 will use the values entered through these options as target specifications. Temperature gradient option 1 on the other hand will interpret values entered through using plasma option 3 or 4 as the mid-plane temperatures and densities.

.. _P04:
P04 : Piece-Wise Background Plasma Option Inputs
  These input lines describe the options to be overlaid onto various pieces of the background plasma. For example, this can be used to allow for detachment of only a few rings at the inside target while solving for all the rest of the half-rings normally. The input consists of two values. I line indicating the number of lines of input - or pieces to be overlaid on the base background. This is followed by the specifications of options to be used for each piece including the rings to which those options should be applied.
  
  e.g.

  .. code-block::

    ' ' 'BG PLASMA Options by Ring (PlasDec Opts 90 & 91) '
    '  R1, R2, Sect, PlasDec, SOL, Teg, Tig, Core, Efield ' 2
        1  16     3        4    0    0    0     1       3
       17  28     2        4   21    0    0     0       3

  The input specifications are as follows. R1 to R2 represent the range of rings to be affected. In this case 1 to 16 and 17 to 28 respectively. The next integer represents the section of the ring to be affected. Section 1 = the first section of the ring (IK=1 to the midpoint) or the OUTER target for JET grids (INNER target for SONNET grids). Section 2 = the second section of the ring (IK = midpoint to NKS(IR)) or the INNER half of the ring for JET Grids (OUTER half for SONNET grids). Section 3 = the entire ring. PlasDec is the Plasma Decay option to be applied to the specified region. SOL is the SOL option for the specified region. Teg and Tig are the temperature gradient options to be applied. Finally, Core and E-field are the core and e-field options to be applied to the specified region. If the region is not a core ring, the core option will be ignored. Similarly, the SOL options will have no effect if a core ring region is being calculated. 

.. _P05:
P05 : Trap Tgrad option
  **Trap Tgrad 0**: Off. Tgrad options are not applied in the trapped plasma. Temperature and density are constant.

  **Trap Tgrad 1**: On. The specified temperature gradient or SOL option are applied to the trap region as if they were a standard SOL ring.

  **Trap Tgrad 2**: On. The Private Plasma Conditions are completely specified by two (x,f(x)) points for each of density, electron temperature, ion temperature and velocity. Each quantity is constant out to the mid-plane at the value of the second point. The specific input parameters defining these points are described later in this document (`P40-P55`_).

  **Trap Tgrad 3**: On. Plasma conditions in the private plasma region are calculated from experimental Thomson measurements. All data for each flux tube are averaged and this value is then assigned to all cells on the flux tube. The target conditions are taken from the input Langmuir probe values.

  **Trap Tgrad 4**: On. Plasma conditions in the private plasma region are calculated from experimental Thomson measurements. All data for each flux tube are averaged and this value is then assigned to all cells on the flux tube. The target conditions are then set to equal the flux tube values. 

.. _P06:
P06 : SOL Enhancement Factor - Electric Field - SOLEF
  To allow selective switching on/off of the electric force, for use with any SOL option. Normally set this value to 1.0, or set to 0.0 to switch off electric field.

.. _P07:
P07 : SOL Enhancement Factor - Drift Velocity - SOLVF
  Similarly, allows switching on/off of drift velocity - set to 0.0 or 1.0

.. _P08:
P08 : SOL 1a Factor - fl
  See SOL Option 1 (`P01`_).

.. _P09:
P09 : SOL 1a Factor - fs
  See SOL Option 1 (`P01`_).

.. _P10:
P10 : SOL 10 Reversal Mach Number - fRM
  See SOL Option 10 (`P01`_).

.. _P11:
P11 : SOL 10 Factor - kin
  See SOL Option 10 (`P01`_).

.. _P12:
P12 : SOL 10 Factor - kout
  See SOL Option 10 (`P01`_).

.. _P13:
P13 : SOL 10 Factor - fRmin
  See SOL Option 10 (`P01`_).

.. _P14:
P14 : SOL 10 Factor - fRmax
  See SOL Option 10 (`P01`_).

.. _P15:
P15 : SOL 6 & 7 Vb Length Factor 1 - VbL1
  See SOL 6 & 7 (`P01`_).

   These quantities VbL1, VbM1, VbL2, and VbM2 (`P15`_-`P18`_) control how the velocity is calculated in SOL options 6 and 7. In option 6, the Velocity is the Base Value entered previously out to VbL1 * SMAX then it takes the value VbM1 * Vhout until VbL2 * SMAX and then finally it is VbM2 * Vhout for S > VbL2 * SMAX and S < SMAX/2. The velocity is symmetric from each target except for a change in sign. SOL option 7 is similar except each of the values (0.0 , vhout) , (VbL1*SMAX, VbM1*vhout), (VbL2*SMAX, VbM2*vhout), (SMAX/2.0 , 0.0) are all linked by linear ramps from one point to the next, giving significant flexibility in assigning background flow velocities. (Including approximations to flow reversed situations.) 

.. _P16:
P16 : SOL 6 & 7 Vb Multiplication Factor 1 - VbM1
  See SOL 6 & 7 (`P01`_).

.. _P17:
P17 : SOL 6 & 7 Vb Length Factor 2 - VbL2
  See SOL 6 & 7 (`P01`_).

.. _P18:
P18 : SOL 6 & 7 Multiplication Factor 2 - VbM2
  See SOL 6 & 7 (`P01`_).

.. _P19:
P19 : Power Density - P/A
  This specifies the power flow along the field lines for Temperature gradient option 2 (`Q01`_).

.. _P20:
P20 : Parallel heat conduction coefficient - K0
  This is the conduction coefficient used for both electrons and ions in Temperature gradient options 2, 3, 4 and 6 (`Q01`_), and exclusively for electrons in options 5 and 7 as well as other SOL options.

.. _P21:
P21 : Parallel ion heat conduction coefficient - K0i
  This is the conduction coefficient used for ions in Temperature gradient options 5 and 7 (`Q01`_) as well as other SOL options.

.. _P22:
P22 : Electric field option - overrides other E-field options or data
  This option allows the behaviour of the electric field to be specified - including over-riding any electric field read in from a background plasma file.

  **Electric field 0**: Electric field as read from file is used.

  **Electric field 1**: Electric field is set to zero everywhere.

  **Electric field 2**: For use with simple SOL options. E -> 0 for S > a specified fraction of SMAX.

  **Electric field 3**: Electric field calculated using the standard formula based on pressure and temperature gradients is used everywhere. This replaces any other electric field that has been read in or specified. Note that in the following formula - T is in eV and p=nT where T is again in eV. Thus there are several instances where the factor e (the electronic charge) cancels out.

    :math:`E_{field}(s) = - \frac{1}{n} (\frac{dp}{ds}) - 0.71 \frac{dT}{ds}`

  **Electric field 4**: The formula from option 3 is used for half-rings that are deemed to be "collisional". For half-rings that are not "collisional", the following formula is used.

    :math:`F_{field}(s) = - \frac{\bar{T_e}}{2 L_{source}}\ \ [V/m]`

    "Collisionality" is determined by an ad hoc method. If the electron temperature at the midpoint of the ring is less than a certain multiple of the target Te then the half-ring is deemed to be non-collisional. The value of this multiplier is entered in the input data file as the second entry after this option.

    The source length used above is also arbitrary. The value for the initial source length is given as a multiplier of SMAX for a given ring. This is also described below. If PIN has been run then the value of Lequiv (the source equivalent length) that is calculated after PIN has been run, is used instead of this specified Lsource value.

.. _P23:
P23 : Electric Field Option 4 - Source Length Specifier
  This value specifies the length of the source to be used on the first iteration of the overriding E-field option. If PIN has been run then the calculated values of the ion source equivalent length at each end of each flux tube are used instead of this imposed value.

.. _P24:
P24 : Electric Field Option 4 - Collisional Determination Factor
  This value is used to decide whether a particular half flux tube is to be treated as if it was non-collisional. If the electron temperature at the midpoint is less than this value times the target electron temperature then this half flux tube will be treated as non-collisional for the purposes of this E-field option. The logical oddity arises from the fact that most of the background plasmas to which this option will be applied will have been calculated with the basic fluid assumption that the background is collisional throughout the range of solution.

.. _P25:
P25 : Ionization Source - Characteristic Length - SOL 12 to 15 - Ls
  This quantity defines the characteristic length for the ionization source options (`P31`_) used in the SOL models in SOL options 12 to 15 as a proportion of the total length of the SOL.

.. _P26:
P26 : Ionization Source - Second Characteristic Length - L2
  This defines the second characteristic length for ionization source options 4 and 5 (`P31`_) which involve a combination of linear and exponential functions superimposed.

.. _P27:
P27 : Ionization Source - Source Fraction - Fi
  This quantity specifies the ratio of source strengths between the two discrete ionization sources that are superimposed in ionization options 4 and 5 (`P31`_).

.. _P28:
P28 : Radiation Source - Characteristic Length - SOL 12 to 15 - Lr
  This quantity defines the characteristic length for the radiative source options (`P32`_) used in the SOL models in SOL options 12 to 15 as a proportion of the total length of the SOL.

.. _P29:
P29 : Radiative Power Constant - SOL12 to 15 - Pr/A (W/m2)
  This defines the radiative loss constant used in the power loss model for SOL options 12 to 15.

.. _P30:
P30 : Radiation Source strength fraction - Frr
  This quantity is used in radiation options 2 and 3 (`P32`_) to specify the radiation source strength relative to the total power at the plates. The exact description is in the equations under the respective radiation options.

.. _P31:
P31 : Ionization Source Option - SOL12 to 15
  **Ionization Option 0**: Ionization constant over 0 < s < Ls * SMAX

    With :math:`S_0 = - N_0 \times V_0`

  **Ionization Option 1**: Exponential decay described by:

    :math:`S_i(s) = S_0 e^{-s / (L_s SMAX)}`

    :math:`S_0 = \frac{-N_0 V_0}{L_s SMAX (1-e^{-L / (L_s SMAX)})}`

  Options 2 and 3 are equivalent to the Ionization source options that can be initially specified, except that they only apply to PIN ionization data and because of this can not be specified initially since PIN data is not immediately available.

  **Ionization Option 2**: The PIN ionization data for each ring is normalized so that the total ionization (to the midplane) is equal to the outflow at the plate.

  **Ionization Option 3**: The PIN ionization data is utilized in an unnormalized form so over/under ionization on each flux tube is allowed and other effects arising from the background ionization distribution can be realized. (e.g., Flow recirculation).

    In addition, the following normal ionization source options have been added, in conformance with TN740.

    *It seems there is missing documentation here.*

  **Ionization Option 4**: Combination of two constant sources of the type shown in option 0. The combination is described by the following equation.

    :math:`S(s) = -(1-F_i) \frac{N_0 V_0}{L_s} F_i \times \frac{N_0 V_0}{L_2}\ \ for\ 0 < s < L_s`

    :math:`S(s) = -F_i \frac{N_0 V_0}{L_2}\ \ for\ L_s < s < L_2`

  **Ionization Option 5**: Combination of a constant source (as in option 0) and an exponential source (as in option 1). The combination is described by the following equations.

    :math:`S(s) = S_0 e^{-s / L_s} - F_i \frac{N_0 V_0}{L_2}`

    :math:`S_0 = -(1-F_i) \frac{N_0 V_0}{L_s (1-e^{-L_2 / L_s})}`

.. _P32:
P32 : Radiative Source Option - SOL 12 to 15
  **Radiative Option 0**: Constant over 0 < s < Lr * SMAX

    With Pr/A given previously.

  **Radiative Option 1**: Exponential decay described by:

    :math:`\frac{P_r(s)}{A} = \frac{P_r0}{A} e^{-s / (L_r SMAX)}`

    With Pr0/A given previously.

  **Radiative Option 2**: Constant over 0 < s < Lr * SMAX with the constant value given by the following equation.

    :math:`\frac{P_r}{A} = F_{rr} \frac{P/A}{L_r}`

    Where F\ :sub:`rr` is the radiative source strength fraction that is defined previously.

  **Radiative Option 3**: Exponential decay of source described by the following equations.

    :math:`\frac{P_r(s)}{A} = \frac{P_{r0}}{A} e^{-s / L_r}`

    :math:`\frac{P_{r0}}{A} = F_{rr} \frac{P/A}{L_r (1-e^{-SMAX/2L_r})}`

    Where F\ :sub:`rr` is the radiative source strength fraction that is defined previously.

.. _P33:
P33 : Imaginary Root Option
  This option is another flag that is applicable only to SOL options 12 to 15. Some of the equations that are solved for the background plasma involve a square root. The quantities in the square root can become imaginary. This is representative of a transition from sub-sonic to super-sonic flow. There are two options available here. 

  **Option 0**: The flow is not restricted and any imaginary quantities are set to zero. 

  **Option 1**: The local velocity is capped at the sound speed and the equation that involves the calculation of the square root is not utilized if it will involve an imaginary result.

.. _P34:
P34 : Flux Recirculation Option - SOL 12 to 15
  This flag turns flux recirculation in SOL options 12 to 15 on and off. Flux recirculation is allowed for via over and under ionization on each flux tube. The source characteristics when the flux recirculation option is on are specified by the next entry. Note that if a multiplication factor of 1.0 is used, this is equivalent to no flux recirculation since the ionization on the ring is equal to the net influx. However, it does allow the following entry to be used to specify different source characteristics for each ring.

  **Flux Recirculation Option 0**: Flux Recirculation Off

  **Flux Recirculation Option 1**: Flux Recirculation ON

.. _P35:
P35 : Flux Recirculation - Source Specifications
  The first line contains the description of the entry. The second line has headings and then the number of lines to follow in the table of flux recirculation data. The contents of the table are:

    Ring number, Flux multiplication factor, Source Length, Source Decay Length

  These specify the characteristics of the ion source for each ring. The Flux multiplication factor represents the quantity of over/under ionization to be found on that ring/flux-tube. A value of 1.0 for this results in no flux recirculation. (i.e. outflux=ionization) A value of 0.5 means that there will be 1/2 as much ionization on the flux-tube as there is flux to the plates and a factor 2.0 will result in twice as much ionization. This results in flow reversal or increased flow to the plates occurring in the calculation of the background velocities in SOL 12 to 15. The Source Length and Source Decay Length are the same as are generally specified except that they can be customized for each flux tube.

.. _P36:
P36 : Iterate SOL Option
  **Iterate SOL Option 0**: Off - Do NOT calculate the SOL iteratively if PIN has been run. This allows the use of PIN data without requiring a recalculation of the background characteristics.

  **Iterate SOL Option 1**: On - Re-calculate the characteristics of the SOL a specified number of times (`P39`_), recalculating each time after PIN has been executed. If the correct options are selected, this will use the PIN generated data in the SOL recalculation.

.. _P37:
P37 : Secondary SOL option
  This specifies the SOL option (`P01`_) to be used in the plasma background recalculation. It is allowed to differ from the SOL option originally specified. Of course, if a SOL option is chosen that does not use the data generated by iteration then the calculated background plasma will be unchanged. This may be useful wen trying to stabilize the effect of PIN puffing. If this option is specified as "-2" then the iteration will automatically use the same SOL option as was specified for the first iteration.

.. _P38:
Ionization Option for Iterative SOL
  This is the same as the ionization source options specified previously for SOL 12 to 15 (`P31`_), except that only ionization options 2 and 3 are meaningful. Option 2 will normalize the ionization data on the flux tube so that the total ionization is equal to the total flux to the plates. This then ensures no flux recirculation and simply uses the PIN data to distribute the ionization along the flux tube. Option 3 takes the PIN ionization data directly (it is not adjusted or normalized except to convert to MKS units) and allows over and under ionization to occur and thus allows for the possibility of flow reversal and its effects on the background plasma characteristics.

.. _P39:
P39 : Number of Pin/SOL Iterations
  This quantity specifies the number of iterations (`P36`_) that will be made in a DIVIMP SOL option -> PIN -> DIVIMP SOL option iteration sequence in an attempt to generate a convergent background solution with the aid of the hydrogenic neutral code.

.. _P40-P55:
P40 - P55 : Private Plasma (Trap) Specification Option Inputs
  These inputs are used for Tgrad option 2 and SOL22 private plasma ionization option -2. The 16 parameters specified here provide a complete description of the density, temperatures and velocity of the background plasma in the private plasma. This completely prescriptive option was implemented because of mounting evidence that the models of the background plasma OSM that are appropriate for the main SOL may not apply in the private plasma. Furthermore, there are few diagnostics and a limited understanding of the physics that lead to the actual conditions in the private plasma. As such, it seemed best to allow for the option of prescribing the conditions in a manner appropriate for the specific case under consideration.

  This option uses the specified parameters to impose a two piece step-wise linear prescription for the density, temperatures and velocity. To make the velocity go to 0.0 at the midpoint, it is necessary to set VBF2 = 0.0.

  The sixteen parameters are:

    Electron temperature (Te): TES1 (P40), TEF1 (P41), TES2 (P42), TEF2 (P43)

    Ion temperature (Ti): TIS1 (P44), TIF1 (P45), TIS2 (P46), TIF2 (P47)

    Density (Nb): NES1 (P48), NEF1 (P49), NES2 (P50), NEF2 (P51)

    Velocity (Vb): VBS1 (P52), VBF1 (P53), VBS2 (P54), VBF2 (P55)

  They are employed in the following fashion where Q can be interpreted as Te, Ti, Nb or Vb.

    :math:`Q(s) = Q_0 + Q_0 (F_1 - 1) \frac{s}{S_1 \times SMAX}\ \ for\ S \le S_1 \times SMAX`

    :math:`Q(S) = F_1 Q_0 + Q_0 (F_2 - F_1) \frac{s - S_1 \ times SMAX}{S_2 \times SMAX - S_1 \times SMAX}\ \ for\ S_1 \ times SMAX \le S \le S_2 \times SMAX`

    :math:`Q(s) = F_2 Q_0\ \ for\ S \ge S_2 \times SMAX`

.. _P56-P59:
These four parameters (P56 - P59) are the quantities Vf1, Tf1, Vf2 and Tf2 described in the Core Option 4 and 5 entries much earlier in this document (`P02`_). They are multiplied by the core field line lengths to obtain the distances over which the core Marfe descriptions of options 4 and 5 are applied.

.. _P56:
P56 : Input Parameters for Core Option 4 and 5 (Marfe Simulation) - Vf1
  Velocity distance factor 1.

.. _P57:
P57 : Input Parameters for Core Option 4 and 5 (Marfe Simulation) - Tf1
  Temperature distance factor 1.

.. _P58:
P58 : Input Parameters for Core Option 4 and 5 (Marfe Simulation) - Vf2
  Velocity distance factor 2.

.. _P59:
P59 : Input Parameters for Core Option 4 and 5 (Marfe Simulation) - Tf2
  Temperature distance factor 2.

Q Tags
------
.. _Q01:
Q01 : TeB Gradient option
  Electron temperatures (TeB0, TeBP, TeBout, TeBin, TeBt) for the following options are specified at `Q06`_-Q11`_. The electron temprature gradient factors (feBL1, feBL2, feBt, feB2) are specified at `Q12`_-`Q15`_. 

  **TeB Gradient 0**: Linear, from TeB0 x feBt at Target to TeB0 at feBL x Smax, then constant

  **TeB Gradient 1**: Linear from Teb0 x Febt at target to TeB0 x Febt2 at Febl x Smax to TeB at Febl2 x Smax, then constant.

  **TeB Gradient 2**: P/A driven gradients. P/A (`P19`_) and K0 (`P20`_) are given

  **TeB Gradient 3**: P/A driven gradients. K0 (`P20`_) given. Based on input data/ring.

    :math:`\frac{P}{A} = (2kT_{iBP} + 5kT_{eBP}) \times N_{BP} \times V_{SBP}`

  **TeB Gradient 4**: P/A driven gradients. K0 (`P20`_) given. P/A calculated as above. The input data for inner and outer plates may differ. A factor of 7/4 is assumed in the heat transport equation.

  **TeB Gradient 5**: P/A driven gradients. K0 (`P20`_) given. The input data for inner and outer plates may differ. A factor of 7/4 is assumed in the heat transport equation.

    :math:`\frac{P}{A} = 5kT_{eBP} \times N_{BP} \times V_{SBP}`

  **TeB Gradient 6**: Identical to option 4 except a factor of 7/2 is assumed in the heat transport equation.

  **TeB Gradient 7**: Identical to option 5 except a factor of 7/2 is assumed in the heat transport equation.

  **TeB Gradient 98**: Read data from DIVIMP generated background plasma file - this option must be used in combination with plasma decay option 98.

  **TeB Gradient 99**: Read from file.

.. _Q02:
Q02 : TiB Gradient option
  Ion temperatures (TiB0, TiBP, TiBout, TiBin, TiBt) for the following options are specified at `Q16`_-Q21`_. The ion temprature gradient factors (fiBL1, fiBL2, fiBt, fiB2) are specified at `Q22`_-`Q25`_. 

  **TiB Gradient 0**: Linear, from TiB0 x fiBt at Target to TiB0 at fiBL x Smax, then constant

  **TiB Gradient 1**: Linear from Tib0 x Fibt at target to TiB0 x Fibt2 at Fibl x Smax to TiB at Fibl2 x Smax, then constant.

  **TiB Gradient 2**: P/A driven gradients. P/A (`P19`_) and K0 (`P20`_) are given

  **TiB Gradient 3**: P/A driven gradients. K0 (`P20`_) given. Based on input data/ring.

    :math:`\frac{P}{A} = (2kT_{iBP} + 5T_{eBP}) \times N_{BP} \times V_{SBP}`

  **TiB Gradient 4**: P/A driven gradients. K0 (`P20`_) given. P/A calculated as above. The input data for inner and outer plates may differ. A factor of 7/4 is assumed in the heat transport equation.

  **TiB Gradient 5**: P/A driven gradients. K0i (`P21`_) given. The input data for inner and outer plates may differ. A factor of 7/4 is assumed in the heat transport equation.

    :math:`\frac{P}{A} = 2kT_{iBP} \times N_{BP} \times V_{SBP}`

  **TiB Gradient 6**: Identical to option 4 except a factor of 7/2 is assumed in the heat transport equation.

  **TiB Gradient 7**: Identical to option 5 except a factor of 7/2 is assumed in the heat transport equation.

  **TiB Gradient 98**: Read data from DIVIMP generated background plasma file - this option must be used in combination with plasma decay option 98.

  **TiB Gradient 99**: Read from file

.. _Q03:
Q03 : Forced Flat Temperature Gradient Option
  Flatten Gradient Opt 0: OFF. Temperature profiles are not flattened.

  Flatten Gradient Opt 1: ON. Temperature profiles are flattened at their current values for S or (SMAX-S) > SMAX * Fcut

  Flatten Gradient Opt 2: ON. Temperature profiles for S > SMAX * Fcut are limited to a maximum of the value at the position SMAX * Fcut.

.. _Q04:
Q04 : Te Gradient Cut-off for Flattening Option
  This is the Fcut factor that will be used for the Electron Temperature profile flattening. If this value is specified as 0.0 then no flattening or cutting off of the temperature rise will occur.

.. _Q05:
Q05 : Ti Gradient Cut-off for Flattening Option
  This is the Fcut factor that will be used for the Ion Temperature profile flattening. If this value is specified as 0.0 then no flattening or cutting off of the temperature rise will occur. 

.. _`Q06-Q11`:
The "standard" plasma temperature and density profiles are described in Note 336. The subscript "e" refers to electrons, and "i" refers to ions. In plasma decay option 1, the quantities TeBout,TiBout and NBout are assumed to be characteristic decay lengths, and are used to determine an exponential decay along the Reference Line for the SOL. When plasma option 99 is used, these values are all ignored and the temperature and density profiles are read in from a file. These values are also ignored - except for the trap values - when various plasma decay and/or temperature gradient options are specified. Furthermore, these values can be used in conjunction with plasma decay options that read values from the input data file and temperature gradient options to generate a variety of background plasma profiles. 

.. _Q06:
Q06 : Temperatures - TeB0
  Temperature of electrons at the mid-point of the field line. For options which use this quantity. May be read using some Plasma Decay Options as an alternative (eV). For use with `Q01`_.

.. _Q07:
Q07 : Temperatures - TeBP
  Temperature of electrons at the targets. May be read using some Plasma Decay Options as an alternative (eV). For use with `Q01`_.

.. _Q08:
Q08 : Temperatures - TeBout
  Outboard step for scaling electron temperature. In some cases it will be a linear step in eV and in other cases it might represent a decay length or e-folding distance for electron temperature across the target into the main SOL (eV). For use with `Q01`_.

.. _Q09:
Q09 : Temperatures - TeBin
  Electron temperature step per ring moving inboard from the separatrix into the core plasma. This may be over-ridden by some of the plasma decay or core options (eV). For use with `Q01`_.

.. _Q10:
Q10 : Temperatures - TeBt
  Electron temperature of the trapped plasma - when a constant trap temperature option is in use (eV). For use with `Q01`_.

.. _Q11:
Q11 : Temperatures - TeBouP
  Step for scaling electron temperature at the target into the private plasma region. Used for exponential decay only. (at this time (eV). For use with `Q01`_.

.. _`Q12-Q15`:
Several temperature gradient options are available, as described in notes 351 and 388 and later. The gradients can be switched off by setting feBt,feB2,fiBt and fiB2 to 1.0. The temperature gradient options should be set to 99 when file data is being used. Similar to the temperature parameters - these will only be used when the appropriate gradient options are specified. The parameters labeled "L1" and "L2" specify distances along the field line in units of SMAX. The Bt and B2 quantities are multiplication factors that modify the "base" temperature for the ring. The plasma values assigned to each cell are linearly interpolated (depending on the temperature gradient options selected) between these points.

.. _Q12:
Q12 : Gradient Parameter - feBL1
   This is the first length factor. At a distance of (feBL1 * SMAX) from the target the temperature will rise to a value of (feB2 * Specified Upstream Temperature). 

.. _Q13:
Q13 : Gradient Parameter - feBL2
  This is the second length factor. At a distance of (feBL2 * SMAX), the temperature rises to the value specified for the Upstream Temperature and remains constant from this point until the mid-point of the ring - where it will meet the specification coming from the other target.

.. _Q14:
Q14 : Gradient Parameter - feBt
  Target multiplier factor. The temperature at the target is (feBt * Specified Upstream Temperature). The upstream temperature is either specified as TeB0 or may be entered using a plasma decay option that specifies a different value for each ring. 

.. _Q15:
Q15 : Gradient Parameter - feB2
  Multiplication factor for distance (feBL1 * SMAX)

.. _Q16: 
.. _Q17: 
.. _Q18:  
.. _Q19:
.. _Q20:
.. _Q21:
.. _Q22:
.. _Q23:
.. _Q24:
.. _Q25:
.. _Q16-Q21:
.. _Q22-Q25:
Q16 - Q25 : Ion Temperatures and Gradients
  Identical to the electron options above, just for ions instead. 

.. _Q26:
.. _Q27:
.. _Q28:
.. _Q29:
.. _Q30:
.. _Q31:
.. _Q26-Q31:
Q26 - Q31 : Densities
  Identicial to the electron temperature options above, just for the density instead.

.. _Q32:
Q32 : Langmuir Probe Data Switch
  This option specifies the form of the data entered in the Langmuir Probe Data input described below.

  **Probe Switch Option 0**: The third column is interpreted as target densities (m\ :sup:`-3`)

  **Probe Switch Option 1**: The third column is interpreted as target saturation currents. (Isat in A/m\ :sup:`2`)

.. _Q33:
Q33 : INNER/Both Target Data Multipliers
  This input item is a set of three numbers on the one line separated by spaces. These numbers are used to multiply the input data entered in the Langmuir Probe Data entry. This allows the input data to be quickly and easily modified. This is especially useful for cases where it may be suspected that Ti is not equal to Te. The order of the three multipliers is Te multiplier, Ti multiplier, Nb multiplier

  The following should be the default value of this input.

  .. code-block::

    'INNER target data multipliers (Ti,Te, Nb) : ' 1.0 1.0 1.0

.. _Q34:
Q34 : Langmuir Probe Data Input - Inner/Both plate
  This section is used to specify a set of data to be used as the base temperature and density for each specified ring. Both the ion and electron temperatures as well as the density can be independently specified for each ring. The data should be entered in ascending ring number order. The format of a line of data is as follows:

  .. code-block::

     ' ' 'Probe data at inner plate (opt 4) or both(opt 3)' 
     ' Ring   TeBP   TiBP       NBP        Number of rows:' 3 
          9   30.0   30.0   1.00e19
         10   15.0   20.0   8.10e18 
         11    5.0   14.5   1.93e18 

  For rings in the SOL for which a line of data is not specified, the values for the next inward - i.e. lower numbered ring are used. These values will be assigned to all bins on the specified ring - depending on the plasma decay option specified. Variations in temperature and density will be caused by imposed temperature gradient options ... or by various SOL options which apply a SOL model to calculate the values throughout the SOL based on the plate input data provided here. The values specified here are initially assigned to the respective halves of the grid up to the mid-plane for the inner and outer targets respectively.

.. _Q35:
Q35 : OUTER Target Data Multipliers
  This is equivalent to the INNER/Both target data multiplier described above. It applies to the OUTER target data that is entered after it.

.. _Q36:
Q36 : Langmuir Probe Data Input - Outer plate
  See above. Provides input data for the outer plate for cases where it differs.

.. _Q37:
Q37 : Core Plasma Input Data
  This section is used to specify a set of data for the temperature, density, and velocity used by the various core plasma options. These values will either specify the temperatures and densities to apply to the entire ring in the case of core option 1 - or they specify the values of temperature and velocity (and density) that apply at the X-point in the case of the other core options. Core Options 2 to 5 are used to simulate an X-point Marfe within the core plasma. The density is ignored for core option 3 to 5. It is calculated by applying conservation of pressure along the field line. When core options 2, 3 and 4 are in use, the base temperatures for the core rings are calculated using the standard core specification procedures. For core option 5, the base values of temperature and density for the core rings must be entered in the tables for Langmuir Probe Data described above; indexed by their core ring numbers.

  The format of a line of data is as follows - the data should be entered in ascending ring number order.

  .. code-block::

     ' ' 'CORE Plasma data - for Core Options 1 and 2      '
     ' Ring   TeBP   TiBP      NBP       Vb Number of rows:' 3
          4    5.0   10.0   1.0e19      0.0
          5    3.0    3.0   5.0e19    200.0
          6    2.0    2.0   1.0e20   2000.0

  In the case of core option 1 - for rings in the CORE for which a line of data is not specified, the values for the next inward - i.e. lower numbered ring are used. These values will initially be assigned to all bins on the specified ring. For the other core options, the data are applied ONLY to the specified rings.

.. _Q38:
Q38 : Inboard plasma flow velocity - Not supported - Vhyin (m/s)
  NOT SUPPORTED IN DIVIMP. For most cases we only have flow of the background plasma for the SOL, in which case this item should be set to 0.0. However, to allow for plasma flow inboard this item can be set to a constant value. Note that this option only applies to the inboard region and does not affect the SOL. This item described for LIM in Note 118.

.. _Q39:
Q39 : Inboard electric field - Not supported - Eyin (V/m)
  NOT SUPPORTED IN DIVIMP. Again, generally set to 0.0. A positive value gives an inboard electric field.

.. _Q40:
Q40 : Outboard plasma flow vel (SOL 5,6 & 7) - Vhyout(m/s)
  Required when SOL 5, 6 or 7 has been selected, otherwise set this value to 0.0. For details of SOL5 see Note 123. This also sets the base velocity for use in SOL options 6 and 7.

.. _Q41:
Q41 : Outboard electric field (SOL 5,6 & 7) - Eyout (V/m)
  Required when SOL 5, 6 or 7 has been selected, otherwise set this value to 0.0. This also sets the base velocity for use in SOL options 6 and 7.

R Tags
------

S Tags
------

T Tags
------

.. _T01:
T01 : Ionization
  **Ionization Option 0**: rates from S(z,Te) data. Ions not followed after state Z (given)
  
  **Ionization Option 1**: rates from S(z,Te) data. Ionization disabled after state Z (given)
  
  **Ionization Option 2**: rates taken as MAX (S(z,Te)). Ionization disabled after state Z (given)
  
  **Ionization Option 3**: rates from Abels van Maanen with E-I REC. Ions not followed after state Z (given)
  
  **Ionization Option 4**: rates from Abels van Maanen. Ions not followed after state Z (given)
  
  **Ionization Option 5**: rates from Abels van Maanen with E-I REC. Ionization disabled after state Z (given)
  
  **Ionization Option 6**: rates from Abels van Maanen. Ionization disabled after state Z (given)

.. _T02:
T02 : Collision
  *Some of the options here are missing complete documentation. Likely will need to investigate source code to fill out missing gaps.*

  **Collision option 0**:

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

  **Collision option 1**: :math:`\tau_{||} = \infty`, no diffusion outside of given Kspec only. Elsewhere reverts to collision option 0.
  
  **Collision option 2**:

    :math:`\tau_{||} = \frac{m_i}{T_i} \frac{T_i}{6.8 \times 10^4 n_b Z_b Z_{eff} Z_i^2 \lambda}`
  
    where Zeff given outside of given Kspec only.
  
    Elsewhere reverts to collision option 0.
  
  **Collision option 3**:

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    Time between K diff steps = tau para outside of given Kspec only. Elsewhere Time between Y diff steps = DeltaT
  
  **Collision option 4**:

    :math:`\tau_{||} = m_b \frac{T_i}{m_i} \frac{T_i}{9.0 \times 10^4 n_b T_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    When Ti > Tb.Mi/Mb and for rings >= given ring no. Elsewhere reverts to collision option 0. Time between S diff steps = tau para for rings >= given ring no. Elsewhere Time between Y diff steps = DeltaT
  
  **Collision Option 5**: Parallel velocity diffusion

    :math:`\Delta v = \sqrt{\frac{8kT_i}{\pi m_i}}`

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`
  
    Diffusion occurs if: 0 < (random) < dt/:math:`\tau_{||}`
  
  **Collision Option 6**: Parallel velocity diffusion

    :math:`\Delta v = \sqrt{\frac{8kT_i}{\pi m_i}} \frac{dt}{\tau_{||}}`
  
    at every time step

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`
  
  **Collision Option 7**:

    *Documentation fragmented*

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    Diffusive steps in the direction opposite of the particles velocity reverse the sign of that v. For rings greater than IRSPEC-1: (Unless special plasma parameter = 0)
  
    time between s diff steps = (t
  
    Elsewhere:
  
    time between s diff steps = taupara
  
  **Collision Option 8**:
  
    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    S diffusive steps are based on:

    :math:`\frac{2.0 kT_i}{m_i} \tau_{||}`

    Otherwise as Collision Option 0.
  
  **Collision Option 9**:

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    S diffusive steps are based on:

    :math:`\frac{2.0 kT_i}{m_i} \tau_{||}`
  
    For IR < IRSPEC:
  
    time between s diff steps = :math:`\tau_{||}`
  
    Elsewhere:
  
    time between s diff steps = (t
  
  **Collision Option 10**: Parallel velocity diffusion

    :math:`\Delta v = \sqrt{\frac{2kT_i}{m_i}}`

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    Diffusion occurs if: 0 < (random) < (t / :math:`\tau_{||}`)
  
  **Collision Option 11**: Parallel velocity diffusion

    :math:`\Delta v = \sqrt{\frac{2kT_i}{m_i}} \frac{\Delta t}{\tau_{||}}`

    At every timestep:

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`
  
  **Collision Option 12**: Parallel velocity diffusion

    :math:`\Delta v = R_G \sqrt{\frac{2kT_i}{m_i}} \sqrt{\frac{\Delta t}{\tau_{||}}}`

    at every time step, where

    :math:`R_G = \sqrt{-2ln(x_1)} cos(2 \pi x_2)`
  
    x1, x2 are uniform on [0,1]

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`
  
  **Collision Option 13**: Parallel velocity diffusion

    :math:`\Delta v = R_G \sqrt{\frac{2kT_i}{m_i}} \sqrt{\frac{\Delta t}{\tau_{||}}}`

    at every time step, where

    :math:`R_G = \sqrt{-2ln(x_1)} cos(2 \pi x_2)`

    x1, x2 are uniform on [0,1]

    :math:`\tau_{||} = m_i \frac{T_b}{m_b} \frac{T_i}{6.8 \times 10^4 n_b Z_b^2 Z_i^2 Z_{enh} \lambda} (1.0 + \frac{m_b}{m_i})`
  
  **Collision Option 14**: Parallel Velocity Diffusion
  
    Exactly the same as Option 13 EXCEPT that velocity diffusion is turned off (DELTAV=0.0) for S > FACTOR * SMAX for the ring, from each target. FACTOR is specified by the Stgrad parameter.

.. _T03:
T03 : Reiser Parallel Force Calculation Option
  **Reiser option 0**:
  
    Reiser coulomb collison option for calculating the parallel forces is turned off.
  
  **Reiser option 1**:
  
    Reiser coulomb collison option for calculating the parallel forces is turned on. The regular DIVIMP collison and friction options are not used. The Reiser coefficients are constant for each cell.
  
  **Reiser option 2**:
  
    Reiser coulomb collison option for calculating the parallel forces is turned on. The regular DIVIMP collison and friction options are not used. The Reiser transport coeffcients are recalculated at every time step based on the local conditions and impurity particle location along the field line. Background profiles are interpolated between cell centers. This may incur a significant computational cost.

.. _T04:
T04 : Friction
  **Friction option 0**:

    :math:`\tau_{stop} = m_i T_b \sqrt{\frac{T_b}{m_b}} (6.8 \times 10^4 (1.0 + \frac{m_b}{m_i} n_b Z_b^2 Z_i^2 Z_{enh} \lambda)^{-1}`

  **Friction option 1**: :math:`\tau_{stop} = \inf` outside of given Kspec only. Elsewhere reverts to friction option 0.
  
  **Friction option 2**: :math:`\tau_{stop} = \tau_{||}` outside of given Kspec only. Elsewhere reverts to friction option 0
  
  **Friction option 3**:

    :math:`\tau_{stop} = T_i \sqrt{\frac{T_i}{m_i}} (9.0 \times 10^4 (1.0 + \frac{m_i}{m_b} n_b Z_b^2 Z_i^2 Z_{enh} \lambda)^{-1}`

    when Ti > Tb.Mi/Mb and for rings >= given ring no. Elsewhere reverts to friction option 0.
  
  **Friction option 4**:

    :math:`\tau_{stop} = m_i T_b \sqrt{\frac{T_b}{m_b}} (6.8 \times 10^4 (1.0 + \frac{m_b}{m_i} n_b Z_b^2 Z_i^2 Z_{enh} \lambda)^{-1}`
  
    The friction goes to zero for a cell whose mean free path is less than the distance to the target.

.. _T05:
T05 : Heating
  **Heating option 0**:

    :math:`m_i T_b \sqrt{\frac{T_b}{m_b}} (1.4 \times 10^5 n_b Z_b^2 Z_i^2 Z_{enh} \lambda)^{-1}`
  
  **Heating option 1**: :math:`\tau_{heat} = \infty` outside of given Kspec only. Elsewhere reverts to heating option 0.
  
  **Heating option 2**: :math:`\tau_{heat} = 0` outside of given Kspec only. Elsewhere reverts to heating option 0.
  
  **Heating option 3**:

    :math:`\frac{(m_iT_b + m_bT_i)^{3/2}}{1.4 \times 10^5 \frac{m_i}{m_b} n_b Z_b^2 Z_i^2 Z_{enh} \lambda}`

    outside of given Kspec only. Elsewhere reverts to heating option 0

.. _T06:
T06 : CX Recomb
  **CX Recomb option 0**: No charge exchange recombination
  
  **CX Recomb option 1**: Nh = :ref:`Nho<D13>` (constant) and :ref:`Vcx<D16>` = sqrt(2Tb/Mb) where Nho given.
  
  **CX Recomb option 2**: Nh = :ref:`Nho<D13>` (constant) with constant :ref:`Vcx<D16>` (given) where Nho given.
  
  **CX Recomb option 3**: Nh = :ref:`Nhc<D12>`, Constant in core.
  
    Nh = :ref:`Nho<D13>` * exp( -S/:ref:`lamhx<D14>`), Exponential decay from the plates in the SOL.
  
    where :ref:`Nho<D13>`, :ref:`Nhc<D12>` and :ref:`lamhx<D14>` are given. S is the distance from the plates along the field lines.
  
  **CX Recomb option 4**: Nh from PIN with :ref:`Vcx<D16>` = sqrt(2 Tb / Mb)
  
  **CX Recomb option 5**: Nh from PIN with constant :ref:`Vcx<D16>`.
  
  **CX Recomb option 6**: Nh from PIN with Charge exchange Coefficient Data (CCD) taken from ADAS.
  
  **CX Recomb option 7**: ADPAK/INEL CX rates. CX rates are extracted from B2-FRATES or INEL formatted atomic information database. Nh is supplied either by PIN or by loading as an auxiliary quantity to the background plasma specification.
  
  **CX Recomb option 8**: Nh from PIN. <SIGMA V> cx rates have been taken from the PhD thesis of C.F.MAGGI and have been fitted to a three parameter exponential by Tom Rognlien (LLNL).
  
  **CX Recomb option 9**: Nh from PIN. <SIGMA V> cx rates have been taken from the PhD thesis of C.F.MAGGI and have been fitted to a modified three parameter exponential by Tom Rognlien (LLNL) The modified coefficients reduce the CX recombination rates extrapolated from the Maggi data for low temperature conditions.

.. _T07:
T07 : Dperp option
  **Dperp option 0**: constant
  
  **Dperp option 1**: Dperp = Dperp0.Nb0/Nb in SOL&Trap
  
    Constant Dperp0 in Main
  
  **Dperp option 2**: Dperp is held constant along the reference line at knot number NKS(IRSEP)/2 +1. Cross-field transport elsewhere is based on moving particles in proportion to their equivalent cross-field position on the reference line. Transport in the Private Plasma is mapped relative to the adjacent cell on the separatrix at the IK=1 index which is then mapped back to the reference line.
  
  **Dperp option 3**: UNTESTED. Spatially varying Dperp. The value of Dperp in each cell is allowed to changed so that the number of Dperp steps required to cross a cell remains consatnt along a field line. The Private Plasma Dperp is set to match the corresponding cells on the separatrix. The actual cross-field steps are done proportional to cell sizes so that a particle taking a step from one cell to the next can step back to its starting position despite the differing Dperp values in the adjacent cells.
  
  **Dperp option 4**: UNTESTED. As option 3 except that the private plasma Dperp values are varied to try to try to keep the number of steps consistent instead of matching the Dperp on the separatrix.

.. _T08:
T08 : Perpendicular Step Option
  **Perp. Step Opt 0**: Constant - The probability of inward and outward cross-field diffusive steps is constant everywhere at a value of 0.5
  
  **Perp. Step Opt 1**: Geometrically Varying for Core Only - The probability of making an inward or outward cross-field diffusive step within the core region is equal to the ratio of the lengths of the sides parallel to the field lines of a small cell located at the current particle position. In the SOL and private plasma the probability of a cross-field step is as calculated in option 0.
  
  **Perp. Step Opt 2**: Geometrically Varying for the entire grid - The probability of making an inward or outward cross-field diffusive step in all regions is equal to the ratio of the lengths of the sides parallel to the field lines of a small cell located at the current particle position.
  
  **Perp. Step Opt 3**: Geometrically Varying for the entire grid - The probability of making an inward or outward cross-field diffusive step in all regions is equal to the ratio of the lengths of the sides parallel to the field lines of a small cell located at the current particle position. This is calculated for each half-cell independently. 

.. _T09:
T09 : Pinch Velocity Option
  **Pinch Velocity Opt 0**: OFF. No Pinch Velocity Applied.
  
  **Pinch Velocity Opt 1**: ON. Pinch Velocity is applied :ref:`at the specified value<T16>` everywhere on the grid.
  
  **Pinch Velocity Opt 2**: ON. Pinch Velocity is applied ::ref:`at the specified value<T16>` only in the main SOL.
  
  **Pinch Velocity Opt 3**: ON. Pinch Velocity is applied with the :ref:`specified value<T16>` at the separatrix. The pinch is only applied in the core and is scaled proportional to the square of the poloidal field line length as a particle moves deeper into the core.

.. _T10:
T10 :  TeB Grad Coeff option
  **TeB Grad Coeff 0**: :math:`\alpha_e = 0`
  
  **TeB Grad Coeff 1**: :math:`\alpha_e = 0.71Z_i^2`
  
  **TeB Grad Coeff 2**: :math:`\alpha_e = 1.5 (1 - 0.6934(1.3167^{-Z_i}))Z_i^2`
    *An explanation of what this is would be nice.*
  
  **TeB Grad Coeff 3**: :math:`\alpha_e = 0.71Z_i^2`
    Feg is set to zero for S values as measured from either target that are greater than the :ref:`D33<specified value>`.

.. _T11:
T11 : TiB Grad Coeff option
  **TiB Grad Coeff 0**: :math:`\beta_i = 0`
  
  **TiB Grad Coeff 1**: :math:`\beta_i = \frac{-3 (1 - \mu - 5 Z_i^2 \sqrt{2 \mu} \mu (1.1 \mu - 0.35))}{(2.6 - 2 \mu + 5.4 \mu^2)},\ \  \mu = \frac{m_i}{m_i+m_b}`
  
  **TiB Grad Coeff 2**: :math:`\beta_i = \frac{H(Z_O) Z_i^2}{Z_O + \sqrt{0.5 (1 + m_b / m_i)}}`
  
    where :math:`H(Z_O) = 1.56 \frac{(1 + 1.41 Z_O) (1 + 0.52 Z_O)}{(1 + 2.65 Z_O)(1 + 0.285 Z_O)}`
  
    where :math:`Z_O` given at `D28`_.
  
  **TiB Grad Coeff 3**: :math:`\beta_i = \frac{-3 (1 - \mu - 5 Z_i^2 \sqrt{2 \mu} \mu (1.1 \mu - 0.35))}{(2.6 - 2 \mu + 5.4 \mu^2)},\ \  \mu = \frac{m_i}{m_i+m_b}`
  
    Fig is set to zero for S values as measured from either target that are greater than the :ref:`specified value<D33>`. 

.. _T12:
T12 : Temperature Gradient Force - Modification Option
  **T-Grad Mod Opt 0**: Off - the temperature gradient forces are not modified from their values calculated on the basis of fluid assumptions. (i.e. a sufficiently collisional plasma). The modification array is set to 1.0 everywhere.

    :math:`F_{T-grad-mod} = 1.0`

  **T-Grad Mod Opt 1**: On - the temperature gradient forces are modified by a multiplicative factor which is an attempt to include kinetic effects in the underlying force equations. The following formula is taken from the correction term used in the UEDGE code.

    :math:`F_{T-grad-mod} = \frac{1}{1 + \alpha_F (\frac{L_{mfp}}{L_{grad}})^2}`

    :math:`\alpha_F` is a specified parameter.

    :math:`L_{mfp} = \lambda_{ii} + \lambda_{ee}`

    :math:`\lambda_{ii} = 10^{16} \frac{T_i^2}{n_i}`

    :math:`\lambda_{ee} = 10^{16} \frac{T_e^2}{n_e}`

    :math:`L_{grad} = Min(L_n, L_{T_e}, L_{T_i}, L_{press}) = Min(scale\ length)`
  
  **T-Grad Mod Opt 2**: On - This is identical to option 1 except that the value for :math:`F_{T-grad-mod}` is equal to 0.0 for the ion temperature gradient force when within one ion mean free path length of the target.
  
  **T-grad Mod Opt 3**: On - the temperature gradient forces are modified by a multiplicative factor which is an attempt to include kinetic effects in the underlying temperature gradient calculations. This option is based on the methods used in the Garching-B2 code with an additional parameter to provide the ability to change the magnitude of the effect to allow for which specific mean-free path is to be used. This option limit the temperature gradient directly as opposed to generating a multiplicative factor that modifies the final calculated force. The temperature gradients are found using the following relation. (Applied to both electrons and ions.)

    :math:`\frac{dT}{ds} = min(\frac{dT}{ds}_{natural}, 0.3 \frac{T}{\alpha_F \lambda})`

    :math:`\alpha_F` is and input parameter, usually set to 1. Set :math:`\alpha_F` to 0.5 to match the present B2 at Garching (*this statement is a very old one, needs to be confirmed*). 

    :math`\lambda_{e,i} = 1.5 \times 10^{16} \frac{T_{e,i}^2}{n}`
  
  **T-grad Mod Opt 4**: On - This is identical to option 1 except that the value for :math:`F_{T-grad-mod}` is equal to 0.0 for both the ion and electron temperature gradient force when within one ion mean free path length of the target.

.. _T13:
T13 : Poloidal Drift option
  **Poloidal Drift 0**: Off - No poloidal drift velocity imposed.
  
  **Poloidal Drift 1**: On - :ref:`Additional poloidal drift velocity given<T17`. (This velocity is actually imposed along the field lines and NOT directly poloidally.) The poloidal velocity is only applied over :ref:`Additional "background plasma" poloidal drift velocity<T18>` is given. The velocity is added to the background plasma velocity and is coupled to the impurities through the force of friction.(This velocity is actually imposed along the field lines and NOT directly poloidally.) The poloidal velocity is only applied over Special plasma parameter - Rspec
  
  To define where any special collision, friction and heating options are to apply in DIVIMP enter a ring number > 0 here, and then the special options will apply for ring numbers greater than or equal to the value given. Typically, the ring no. of the separatrix will be entered, so that for example Collision option 3 applies just to the SOL and Trapped Plasma, reverting to Collision option 0 in the Main plasma. This parameter has no effect on the Velocity Diffusion collision options 12 and 13. If selected, these options apply to particle transport through-out the plasma. In general, this should be set to 0 so that the selected options will apply to the entire plasma. 

.. _T14:
T14 : Cross-field Diffusion Rate - Dperp (m*m/s)
  The basic perpendicular diffusion coefficient is specified for the whole model space. The values for Dperp can be 0.0, 0.5, 1.0, 3.0, etc. as desired.

.. _T15:
T15 : Cross-field Diffusion Rate for Private Plasma Region - Dperpt (m*m/s)
  This option allows a different cross-field transport rate for the private plasma to be specified. It was needed to mimic some effects in EDGE2D.

.. _T16:
T16 : Perpendicular Pinch Velocity - CVPINCH (m/s)
  Value of the perpendicular pinch velocity. At every time-step it will act on ions to move them cross-field. A positive velocity is used to designate inward motion toward the core or into the private plasma region. 

.. _T17:
T17 : Poloidal Drift Velocity - Vpol (m/s)
  This is the drift velocity used if :ref:`poloidal drift option<T13>` 1 is specified. This velocity is applied along the field lines in the SOL. It is not imposed directly on the motion in the poloidal plane.

.. _T18:
T18 : Poloidal Drift Velocity - Range of Effect
  This entry specifies two values. These are the start and stop points of the range of effect of the poloidal drift velocity. These numbers are expressed as a fraction of SMAX along each individual flux tube. If :ref:`poloidal drift option<T13>` 1 is specified then this velocity is applied along the field lines in the SOL for particles in the region from Factor1 * SMAX < S < Factor2 * SMAX. For example, if Factor1 is specified as 0.1 and Factor2 is 0.9 then the poloidal drift velocity will act on particles when their S-position is in the range 0.1 < S/SMAX < 0.9.

.. _T19:
T19 :


.. _T20:
T20 :


.. _T21:
T21 :


.. _T22:
T22 :


.. _T23:
T23 :


.. _T24:
T24 :


.. _T25:
T25 :


.. _T26:
T26 :


.. _T27:
T27 :


.. _T28:
T28 :


.. _T29:
T29 :


.. _T30:
T30 :


.. _T31:
T31 :


.. _T32:
T32 :


.. _T33:
T33 :


.. _T34:
T34 :


.. _T35:
T35 :


.. _T36:
T36 :


.. _T37:
T37 :


.. _T38:
T38 :


.. _T39:
T39 :


.. _T40:
T40 :


.. _T41:
T41 :


.. _T42:
T42 :


.. _T43:
T43 :


.. _T44:
T44 :


.. _T45:
T45 :


.. _T46:
T46 :


.. _T47:
T47 :


.. _T48:
T48 :


.. _T49:
T49 :


.. _T50:
T50 :


.. _T51:
T51 :


.. _T52:
T52 :


.. _T53:
T53 :


.. _T54:
T54 :


.. _T55:
T55 :


.. _T56:
T56 :


.. _T57:
T57 :


.. _T58:
T58 :


.. _T59:
T59 :


.. _T60:
T60 :


.. _T61:
T61 :


.. _T62:
T62 :

W Tags
------

Z Tags
------

200 Tags
--------

.. _201:
201 : Force Te = Ti
  **0**: Off - Te and Ti are calculated separately applying the source terms that are appropriate for each species in the independent heat transport equations.

  **1**: On - Te and Ti are forced to be equal each other at all points - source terms for the two are combined into one heat transport equation.

.. _202: 
202 : Initially Imposed Target Mach Number
  This is the value of the flow velocity initially imposed at the target as a multiple of the target sound speed. (A value of 1.0 is usually used initially for the target mach number). If the iterative mach solver option is turned on - then the value of the mach number at the target may move from this initial value as the solver searches for a smooth solution at the point of the super-sonic to sub-sonic transition. 

.. _203:
203 : Initial Mach Number Step Size
  When the iterative Mach solver is turned ON the values of the Mach number are initially stepped by this amount as the solver conducts its search. (Typically this value is set to 0.1 - so that mach number solver initially proceeds in increments of 0.1 - trying to bracket the critical target mach number.)

.. _204:
204 : Ultimate Mach Number Resolution
  The solver resolves the Mach number to this level of "accuracy". (Usually 0.00001 is used) . Due to instabilities encountered in the equations, the solution is found to bifurcate at the value of the critical mach number - even for exceptionally small changes in the Mach number of 10\ :sup:`-10` or less. As such, it has proven difficult to actually find a solution that smoothly traverses the transition region when examined on a small scale length. Usually, the solutions containing sonic transitions are adequate given the granularity of the grid on which the simulation is taking place. 

.. _205:
205 : Ionization Source Length Switch
  This option controls the interpretation of the length entries of the ionization source characteristics that are entered below.

  **0**: Source lengths are interpreted to be in absolute units (meters)

  **1**: Source lengths are expressed in relative units as a proportion of SMAX for each individual ring

.. _206: 
206 : Start of Ionization Source
  The interpretation of this number depends on the analytic ionization option selected. This number is interpreted as the starting S position (relative or absolute) of the ionization source. In the case of the triangular or rectangular ionization sources there will be no ionization for S less than the value listed here.

.. _207:
207 : End or Length of Ionization Source
  This specifies the end of the ionization source or effectively its length if the ionization source starts at 0.0. All of the source lengths are limited to a maximum of 1/2 of the field line. This is because the solver operates from each target out to the point mid-way from both targets along the field line and would ignore any source contributions outside this range.

.. _208:
208 : Decay Factor or Width of Ionization Source
  This is the characteristic decay length or the width factor of the ionization source depending on which ionization source option has been selected. In particular, this value is used for the decay length of the default Exponential Decay Ionization Source.

.. _209:
209 : Length of the Radiation Source
  Length of radiating source region in meters.

.. _210:
210 : Decay Length of Radiation Source
  Characteristic decay length for the exponential decay radiation source - Prad Option 1.

.. _211:
211 : Source Strength Fraction (Frr)
  This specifies the total power radiated by the radiation source term in terms of the power flux onto the target. A value of 3.0 means that the integrated strength of radiated losses will total 3 times the total target power flux for the specific ring.

.. _212:
212 : Garching Radiation Model: Parameter Alpha
  The Garching model for impurity radiation uses the following formula. The quantities Alpha, the base temperature and the two exponents are read from the input data file using these lines. The actual option is selected through the Radiative Source Switch (Prad Option 2) which is found below.

  .. math::
    P_{rad} = \alpha n_e^2 L_z(T_e)\ \ \ \ [W/m^2]

  .. math::
    L_z(T_e) = \frac{2.0 \times 10^{-31}}{T_n^{Exp1} + T_n^{Exp2}}

  .. math::
    T_n = Max(10^{-6}, T_e/T_{base})

.. _213:
213 : Garching Radiation Model: Temperature Base
  This is the base temperature in the above radiation formula. It is typically 15 eV.

.. _214:
214 : Garching Radiation Model: First Exponent
  This is the first exponent in the above radiation formula. It is typically 1.5.

.. _215:
215 : Garching Radiation Model: Second Exponent
  This is the second exponent in the above radiation formula. It is typically -3.

.. _216:
216 : Correction factor to Gamma (ion)
  Ion sheath heat transmission coefficient correction factor. This factor adjusts the total ion heat flux onto the targets by adjusting the Gamma factor for the ion heat sheath transmission. It is usually set to zero.

  The formula used to calculate Ion heat flux is (where M is the target Mach number):

  .. math::
    P_{ai} = \gamma_i \times kT_{i0} \times \Gamma_0

  .. math::
    \gamma_i = 2.5 * \frac{1}{2} M^2 (1.0 + T_{e0} / T_{i0}) + \gamma_{i-corr}

.. _217:
217 : Correction factor to Gamma (electron)
  Electron sheath heat transmission coefficient correction factor. This factor adjusts the total electron heat flux onto the targets by adjusting the Gamma factor for the electron heat sheath transmission.

  .. math::
    P_{ae} = \gamma_e \times kT_{e0} \times \Gamma_0

  .. math::
    \gamma_e = 5.0 \times \gamma_{e-corr}

.. _218:
218 : CX Power Coefficient
  Ion charge exchange is one of the power loss mechanisms that play a role in the ion heat equation. Roughly 3/2 kT is removed for each charge exchange event. This coefficient can be used to increase or decrease the number of charge exchange events occurring and thus adjust the significance of the associated terms. This factor is used in most of the Pcx Options. This is usually set to 1.0.

.. _219:
219 : Recycling Source Coefficient
  The recycling source coefficient should be set to 1.0 in most cases. It represents the proportion of target flux that is being recycled into the ionization source. This feature is particularly useful and necessary when dealing with PIN or EDGE2D ionization results where the recycling fraction was intentionally reduced to obtain a "flow" from the core. Effectively assuming that a certain proportion of the total ion source is effectively entering the system through the core and not through the ionization of particles recycling at the targets. This factor is used to multiply the ionization sources obtained from either PIN or EDGE2D and thus is a factor in determining the amount of cross-field flow required to equalize the sources and sinks on the field line. This cross-field flow can be added as a source through the use of the various Gperp options.

.. _220:
220 : Pei (Equipartition) Correction Factor
  This factor multiplies the Pei power transfer term in the Te and Ti heat equations, thus allowing some flexibility in the strength of equipartition. (Pei may be turned on and off completely through the use of the Pei Option Switch described below.) The value SHOULD always be set to 1.0 unless there is a good reason for changing it or the objective is to match results from other sources which may have used an alternative form of equipartition.

.. _221:
221 : Velocity Error Switch
  This switch affects the default action of the R-K solver when it calculates conditions that result in an imaginary quantity for the density. This condition likely arises in the system of equations because not all of the physical effects are necessarily included. Furthermore, the models of the included physical processes, although as good as possible, are approximations in many cases. Thus, over small scale lengths or in regions near or beyond a super-sonic/sub-sonic flow transition region the relationship between pressure, flux, Te and Ti can result in no solution for the density. This is unphysical since it is believed that the density is a smoothly varying physical quantity. However, the product of the density and the velocity, (Flux = :math:`\Gamma = nv`), is a well-defined quantity at all points. When the equations yield an imaginary n, the solver assigns a value of v, either equal to the local sound speed or equal to the last value of v that was correctly calculated and from this value for the velocity and the known parallel particle flux, calculates a local value of n. These methods seem to produce reasonable solutions of n and v for the cases that have been examined.

  **0**: C\ :sub:`s`

  **1**: Constant

  **2**: Pressure Adjustment - additional pressure required carried forward.

  **3**: Pressure Adjustment - additional pressure required is not carried forward.

.. _222:
222 : Distributed Power Start Position
  This parameter is used with power distribution (`270`_) options 7, 8 and 9. It specifies the starting position, as a function of SMAX for the ring, of the region where the power being carried by the electrons and ions should start reducing. This can be used for distributing the target power flux and any volume power terms if present. The simplest assumption in SOL 22 occurs when there are no heat sources or sinks along the ring and the power reaching the target is carried by the electrons and ions all the way from the mid-point ("top"). This is an adequate first order assumption - but in fact the power carried to the target or being radiated along the length of the ring would have likely entered the ring with some spatial distribution and not at the top. The power distribution options address the question of distributing this power loading.

.. _223:
223 : Distributed Power End Position
  This specifies the position (times SMAX for the ring) at which all the power has been input. (See previous entry and notes on power distribution options.)

.. _224:
224 : Compound Gperp - Fraction of Gperp in Rectangular Distribution
  Gperp is the amount of additional cross-field source or sink that is necessary to force particle conservation on an individual flux-tube. This additional source or sink can be distributed in a variety of ways - determined by the Gperp option described below. This parameter specifies the fraction of the required amount that will be added in a rectangular distribution in a region defined by the following two parameters. The remainder is added uniformly over the length of the flux tube. This parameter is used with Gperp options 5 and 6 (`272`_).

.. _225:
225 : Compound Gperp - Start of Region
  For Gperp options involving a specified region for adding the perpendicular flux - this option defines the start of this region. The starting point is calculated by multiplying this factor times the value of SMAX for the ring.

.. _226:
226 : Compound Gperp - End of Region
  For Gperp options involving a specified region for adding the perpendicular flux - this option defines the end point of this region. The ending point is calculated by multiplying this factor times the value of SMAX for the ring.

.. _227:
227 : Gextra Source strength - Target flux multiplier
  If the "extra" source/sink switch is turned on then an additional source is added to the ring over a specified region and a sink of exactly equal magnitude is removed over a second specified region. This artificial addition of a source and sink does not affect the global particle balance on a ring but will significantly affect the flow pattern.

  This first parameter specifies the magnitude of the source/sink to be applied to the ring. The magnitude is expressed as a fraction of the total target outflux for the entire flux tube. Both source and sink will exactly equal the particle influx specified by this parameter resulting in a net influx/outflux of zero. This option is designed to redistribute the influx/outflux on the flux tube in order to replicate experimentally observed flows at the mid-plane.

.. _228:
228 : Gextra Source Start/Stop * SMAX
  This entry consists of two numbers in the range [0.0, 1.0]. These are then multiplied by the Smax distance for the ring to obtain the region over which the additional source is uniformly added.

.. _229:
229 : Gextra Sink Start/Stop * SMAX
  This entry consists of two numbers in the range [0.0, 1.0]. These are then multiplied by the Smax distance for the ring to obtain the region over which the additional sink is uniformly removed.

.. _230:
230 : PP target power loss redistribution range * SMAX
  If the Start Knot Index for EDGE2D Compatibility Option (`231`)...

  This parameter specifies the starting knot number for EDGE2D compatibility option 9 (`269`_). SOL option 22 will start solving at this knot using the fluxes and background plasma conditions taken from EGE2D at this point. The following option describes what is done with the cells between the target and the starting point of the solver.

.. _231:
231 : Fill Option for Skipped Cells in EDGE2D Compatibility Option 9
  This option specifies how the plasma background is calculated for the cells between the target and the starting point of the SOL 22 solver when compatibility option 9 is invoked.

  **0**: Linear interpolation from the target conditions extracted from EDGE2D to the solver solution at the knot where the solver starts.

  **1**: Extrapolation. The cell values and target conditions are linearly extrapolated back to the target from the values found by the solver in the first two cells for which it has a solution. If this extrapolation would yield negative target temperatures then the values are held constant from the start cell to the target at the values for the solver start cell.

  **2**: Constant at solver start cell values. The cells from the target to the start cell are held constant at the start cell values. The target conditions themselves are extracted from the EDGE2D solution and are not changed.

  **3**: Constant at EDGE2D target values. The cells from the target to the start cell are held constant at the EDGE2D target values.

.. _232:
232 : Qe Term - Temperature Cutoff (eV)
  It has been found that some of the power terms calculated by various sources may be erroneous or very destabilizing for the solver at very low temperatures. As a result this cutoff temperature was introduced such that for any cell with a Te < Qe Term Tcutoff - the contributions from PINQE would be set to zero.

.. _233:
233 : PINQID - Atomic Ionization - T cutoff (eV)
  It has been found that some of the power terms calculated by various sources may be erroneous or very destabilizing for the solver at very low temperatures. As a result this cutoff temperature was introduced such that for any cell with a Ti < PINQID atomic Tcutoff - the contributions calculated for the atomic ionization portion of PINQID would be set to zero.

.. _234:
234 : PINQID - Molecular Ionization - T cutoff (eV)
  It has been found that some of the power terms calculated by various sources may be erroneous or very destabilizing for the solver at very low temperatures. As a result this cutoff temperature was introduced such that for any cell with a Ti < PINQID molecular Tcutoff - the contributions calculated for the molecular ionization portion of PINQID would be set to zero.

.. _235:
235 : PINQID - Recombination - T cutoff (eV)
  It has been found that some of the power terms calculated by various sources may be erroneous or very destabilizing for the solver at very low temperatures. As a result this cutoff temperature was introduced such that for any cell with a Ti < PINQID recombination Tcutoff - the contributions calculated for the recombination portion of PINQID would be set to zero.

.. _236:
236 : Qi Term / PINQID - Charge Exchange - T cutoff (eV)
  It has been found that some of the power terms calculated by various sources may be erroneous or very destabilizing for the solver at very low temperatures. As a result this cutoff temperature was introduced such that for any cell with a Ti < PINQID charge exchange Tcutoff - the contributions calculated for the charge exchange portion of PINQID would be set to zero.

  This cutoff term is also used for the PINQI quantity read in from NIMBUS. For Ti < Qi Term Tcutoff - the contributions from PINQI would be set to zero.

.. _237:
237 : PINQID - Charge Exchange Option 1 - Reference Temperature (eV)
  PINQID charge exchange option 1 (`262`_), described in the section on SOL 22 switches, requires a reference temperature in order to calculate the power loss through the ion charge exchange channel.

.. _238:
238 : Minimum Temperature Allowed in Solver (SOL 22)
  This is the minimum temperature allowed in the solver. If temperatures lower than the value entered here are encountered within the solver they are set equal to this value. Usually this is set to zero. However, in some extreme cases, usually when modelling detachment, recombination and charge exchange processes can result in cases where power and particle balances will result in a falling ion temperature - to such an extent that the temperature would be negative unless this is prevented.

.. _239:
239 : Maximum Allowed Temperature Drop Fraction
  This quantity specifies the maximum fractional temperature drop that is allowed on a ring. If the temperature drops to less than this fraction of it's maximum value then the solver generates an error condition and exits for this flux tube. At this point the error correcting code can adjust the solution options for this flux tube and re-run for this ring. This requires that the error correction option be activated. The solver records the maximum temperature reached along the ring and checks to see that the calculated temperature does not fall below this fraction of that value. Specifying a value of 0.0 will turn off this option. Specifying a value of 1.0 is the same as requiring only solutions with constantly increasing or at least constant temperatures along the field line.

.. _240:
240 : Momentum Loss Term Multiplier
  This is an arbitrary variable factor used to scale the momentum loss term if desired. It can be used to enhance or decrease the effect of calculated momentum loss terms for the purposes of exploring the effect, comparing with other codes utilizing a somewhat different model or to proportionally allow for other physical effects that are not a part of the current model. It should usually be set to 1.0. See the descriptions of Momentum Loss formulae in the DIVIMP guide and for the Momentum loss options below.

.. _241:
241 : Friction Factor for Momentum Loss
  This is a variable factor used in the simplified formulae for momentum loss. See the descriptions of Momentum Loss formulae in the DIVIMP guide and for the Momentum loss options below.

.. _242:
242 : Length of Momentum Loss region
  Length of momentum loss region for the simpler analytic momentum loss options. (0.25 is a typical value - distributing the momentum loss over the first quarter of the ring.)

.. _243:
243 : Decay Length of Momentum Loss
  The decay length is used as the scale length or exponential characteristic length in some of the analytic momentum source options. See the descriptions of Momentum Loss formulae in the DIVIMP guide and the descriptions of the Momentum Loss Option switch in the section below.

.. _244:
244 : Correction ratio of CX to Ionization events for Momentum transfer
  This is an extra correction multiplier for the number of CX momentum events (and thus pressure loss) to ionization events for those momentum sources that depend on the ionization source. It is typically set to 1.0 unless one wants to increase the effectiveness of the momentum loss term.

.. _245:
245 : Correction ratio of CX to Ionization events for Momentum transfer
  This is an extra correction multiplier for the number of CX momentum events (and thus pressure loss) to ionization events for those momentum sources that depend on the ionization source. It is typically set to 1.0 unless one wants to increase the effectiveness of the momentum loss term.

.. _246:
246 : Te Cut-off for increased CX multiplier
  At low Te the ratio of CX events to ionization events increases at a great rate - dependent on the local electron temperature. In order to replicate this and thus approximately estimate the momentum loss terms - it is necessary to increase the number of momentum transfer events occurring / ionization. The following formula is used to give the ratio of momentum loss to ionization rates.

  For temperatures above the Tcut value specified in this entry - the multiplier is equal to 1.0. Below this temperature - the multiplier increases rapidly. However, behaviour at very low temperatures is not well understood. So the following entry cuts off the multiplication factor for temperatures lower than the given cut-off.

.. _247:
247 : Te Lower Limit Cutoff for CX Multiplier
  The CX multiplier described in the previous entry returns a value of 1.0 for temperatures below this cut-off value.

.. _248:
248 : PIN Qe Term Multiplier
  This is a factor to allow for scaling of the electron energy loss term that is returned by the hydrogenic neutral code. It should usually be set to 1.0.

.. _249:
249 : PRAD option 3 multiplier (x PINQE)
  This option applies to PRAD option 3 (`256`_). Prad option 3 defines the radiative loss term to be a multiple of the PINQE (electron energy loss term that is calculated by PIN. This parameter specifies the value of the multiplier.

.. _250:
250 : Initial Number of Runge-Kutta Steps Between Grid Points
  This is the initial number of Runge-Kutta steps that the solver assumes between each consecutive grid point where a solution is necessary. The solver will increase or decrease the actual step-sizes taken depending on errors encountered and the solver's estimate of the error on any given step. The usual value entered here is 100.

SOL 22 Switches
The general options to the solver are invoked by a series of switches that can each take a number of values for each option. The following section describes the switches, their functions and the various acceptable input values. 

.. _251:
251 : Ionization Option
  There are, at present, eleven different ionization options that can be chosen. This number changes so if the code version is newer than this documentation then you might want to examine the code and the ECHOSOL subroutine for additional information.

  **0**: Exponential Decay. The ionization source is modelled as an exponential decay falling away from the target. This is the simplest and usually the default option. Lengths of the source region and the characteristic decay are as specified above in the ionization parameter section.

  .. math::
    S_{iz}(s) = S_0 e^{-s/\lambda}
  
  Where S\ :sub:`0` is the normalization factor - usually set so that the integral of the ionization source along the ring will equal the flux to the targets.

  **1**: PIN data. Normalized PIN data is read in and normalized to the target flux for each end of the flux tube individually.

  **2**: PIN data. Unnormalized PIN data is read in and used as is - there is a global normalization check performed to make sure that the integral of ionization over the entire grid is equal to the particle sources (usually target flux and recombination).

  **3**: Triangular Source. The ionization is distributed in a triangular shape between the Start position of the ionization Source to the End position (specified as parameters above). The integral over the triangle is normalized to the target flux for each end of the flux tube.

  **4**: Rectangular Source. The ionization is a rectangular shape from the Start of the ionization source to the End of the ionization source region. (Start and End are parameters specified above). The strength of the ionization is constant over the region. The integral of the source is normalized to the target flux for each end of the flux tube.

  **5**: Algorithmic Source 1. This option chooses between the Triangular and Rectangular sources and their characteristics applied to these sources based on the target conditions for the half-ring.

    The algorithm used is the following:

    If ntarget > 10\ :sup:`19` then a Triangular Source (option 3 is used).

    If Te < 1.3 eV

    Start of Triangular Source = 13 - 10 X Te

    End of Triangular Source = Start + 2.0

    If Te >= 1.3 eV

    Start of Triangular Source = 0.0

    End of Triangular Source = 2.0 meters

    If ntarget <= 10\ :sup:`19` then a Rectangular Source (option 4 is used).

    If Te < 10 eV

    Start of Rectangular Source = 0.0

    End of Rectangular Source = 13.0 - Te

    If Te >= 10 eV

    Start of Rectangular Source = 0.0

    End of Rectangular Source = 2.0 meters

  **6**: s\ :sup:`5` Gaussian. This option selects a source of the form s\ :sup:`5` exp\ :sup:`(-s\ :sup:`2`)` (i.e. s\ :sup:`5` times a Gaussian distribution). This may be a good analytic source for low temperature cases where the ionization is somewhat removed from the target. However, the almost complete lack of ionization immediately adjacent to the target can cause problems with the solver. The solver often needs to use a very small step-size in these regions and may run into conditions resulting in imaginary solutions. Either of these can significantly increase the computational time required by SOL 22 to provide a background plasma solution. The "Decay Length of Ionization Source" (`208`_) specified above is used as the width or decay factor for this distribution - it approximately specifies the location of the peak of this distribution. The source starts at S=0 and is cutoff at the End of ionization source specified in the parameter section. The total integrated source strength for each half ring is normalized to the target flux for each end of the flux tube. This form of a Gaussian was chosen for the source because of it's reasonable spatial distribution and because it is analytically integrable.

    .. math::
      S_{iz}(s) = As^5 e^{-\alpha s^2}

    .. math::
      \alpha = 2.5 / \lambda^2

  **7**: Algorithmic Source 2. This option chooses between the s\ :sup:`5` Gaussian and Rectangular type sources based on target conditions. It is otherwise similar to option 5.

    The algorithm used is the following:
  
    If n\ :sub:`target` > 10\ :sup:`19` then an s\ :sup:`5` Gaussian Source (option 6 is used).
  
    If Te < 1.3 eV
  
    Start of s\ :sup:`5` Gaussian Source = 0.0
  
    End of s5Gaussian Source = End of 1/2 ring
  
    Width Factor = 14.0 - 10 ( Te
  
    If Te >= 1.3 eV
  
    Start of s\ :sup:`5` Gaussian Source = 0.0
  
    End of s\ :sup:`5` Gaussian Source = End of 1/2 ring
  
    Width Factor = 1.0
  
    If n\ :sub:`target` <= 10\ :sup:`19` then a Rectangular Source (option 4 is used).
  
    If Te < 10 eV
  
    Start of Rectangular Source = 0.0
  
    End of Rectangular Source = 13.0 - Te
  
    If Te >= 10 eV
  
    Start of Rectangular Source = 0.0
  
    End of Rectangular Source = 2.0 meters

  **8**: PIN Source Strength. This option runs PIN to obtain the total ionization on each field-line, affected only by global normalization. It then applies the ionization option specified in the initial ionization option on each subsequent iteration except that instead of the ionization source being normalized to the target flux, it is normalized to the amount of ionization on the 1/2 flux tube as determined from PIN.

  **9**: Offset s\ :sup:`5` Gaussian. This ionization option is almost the same as option 6 except that the width factor for the distribution is also simultaneously used as a zero offset. The form of the ionization source is the following. The source is normalized (by setting the factor A) so that the integral over the source is equal to the target flux for the 1/2 ring.

    .. math::
      S_{iz}(s) = A(s+L)^5 e^{\alpha(s+L)^2
  
    .. math::
      \alpha = 2.5 / \lambda^2
  
    and
  
    .. math::
      L = \lambda / 2

  **10**: Algorithmic Source 3. This option chooses between the Offset s\ :sup:`5` Gaussian source and the Rectangular source based on target conditions.

    The algorithm used is the following:
  
    If ntarget > 10\ :sup:`19` then Offset s\ :sup:`5` Gaussian (option 9 is used).
  
    If Te < 1.3 eV
  
    Start of s\ :sup:`5` Gaussian Source = 0.0
  
    End of s\ :sup:`5` Gaussian Source = End of 1/2 ring
  
    Width Factor = 28.0 - 20 ( Te
  
    If Te >= 1.3 eV
  
    Start of s\ :sup:`5` Gaussian Source = 0.0
  
    End of s\ :sup:`5` Gaussian Source = End of 1/2 ring
  
    Width Factor = 2.0
  
    If ntarget <= 10\ :sup:`19` then a Rectangular Source (option 4 is used).
  
    If Te < 10 eV
  
    Start of Rectangular Source = 0.0
  
    End of Rectangular Source = 13.0 - Te
  
    If Te >= 10 eV
  
    Start of Rectangular Source = 0.0
  
    End of Rectangular Source = 2.0 meters

.. _252:
252 : Initial Ionization Option:
  This ionization option is used to generate the seed plasma solution ("Starter Plasma") for PIN iteration. Since PIN needs a plasma solution before it can determine the hydrogenic ionization and power loss terms, which will be used as the OSM iterates to calculate better estimates of the background plasma, it is necessary to start with an analytically calculated plasma. This option specifies the ionization source to be used when the ionization option is set to 1, 2 or 8.

  Options: Match Previous These are identical to the options above. Options 1, 2 and 8 are invalid as the initial ionization option.

  In addition, the following ionization source options are available on the initial iteration.

  **11**: Edge2D Ionization Source This option is only useful when trying to directly compare OSM model solutions to Edge2D background plasma solutions. This option requires that the Edge2D data be read into DIVIMP so that the ionization source, as it was used by Edge2D in it's calculations, is then available within DIVIMP for use by SOL option 22. This can reduce or eliminate differences in the solutions caused by differences in the source terms.

  **12**: PIN is run once with EDGE2D background in the SOL before SOL 22 is invoked. This generates power terms but does not include puffing.

  **13**: PIN is run twice with an EDGE2D background in the SOL before SOL22 is invoked. This generates power terms and a more correct puffing approximation.

  **14**: PIN is run twice with an EDGE2D background everywhere (SOL+CORE) before SOL 22 is invoked. This option must be used in conjunction with core option -1 or else the EDGE2D solution in the core will be overwritten.

  **15**: Ionization source data is read from the EDGE2D input for the case. The EDGE2D plasma solution is assigned as the "previous" solution so that cross-field gradient dependent perpendicular flux and power terms may be assigned correctly.

  **16**: PIN is run once with previously generated DIVIMP background in the SOL before SOL 22 is invoked. This generates power terms but does not include puffing.

  **17**: PIN is run once with previously generated DIVIMP background in the SOL before SOL 22 is invoked. This generates power terms but does not include puffing. In addition, each subsequent iteration through the solver will ensure that the original solution is used for the core and private plasma regions. Only the main SOL plasma is allowed to evolve. This option pre-dated the piece-wise background plasma options and it is recommended that one use the piece-wise method of combining plasma solutions to obtain the same effect. 

.. _253:
253 : Private Plasma Ionization Option
  This option allows one to specify a different analytic ionization model in the private plasma from what is used in the main SOL. The physics in the private plasma are very different from the main SOL and a very different background plasma may result. This makes using the same initial ionization source option in these two regions somewhat invalid since one might expect quite different behaviour. These options will only be active IF the Tgrad option has been set to "1" so that SOL22 will be applied to both main SOL and PP rings.

  **-6**: A uniform plasma is assigned to each ring in the private flux zone from a listing of temperature and density in the divimp input file. The target values are set from the bulk plasma values.

  **-5**: A uniform plasma is assigned to each ring in the private flux zone from a listing of temperature and density in the divimp input file. The target values are not modified.

  **-4**: Experimental Thomson data is applied to the private plasma. the average value of the Thomson data on a ring is assigned to every cell on the ring. Rings without data are interpolated. The target flux is assigned using the thomson data.

  **-3**: Experimental Thomson data is applied to the private plasma. The average value of the Thomson data on a ring is assigned to every cell on the ring. Rings without data are interpolated. The target flux is specified in the divimp input file.

  **-2**: Specified Plasma This option uses a completely arbitrary, specified background plasma (ref) for the private plasma region. This specification is exactly the same as Trap Temperature Gradient Option 2 and uses all the coefficients that are defined for that option. In practice, this option is identical to using Tgrad Option 2. Historically, this option was developed first and then later generalized into Tgrad Option 2 for use in combination with other SOL options.

  **-1**: Matches Previous Same as Initial Ionization Option if the Main Ionization option has been set to 1.0, 2.0 or 8.0. Otherwise, the private plasma ionization option is set to the value of the Ionization option.

  **0 to 15**: These are identical to the options outlined in the Initial Ionization Option.

.. _254:
254 : 5/2 nv kT Term
  This option turns on/off the first convection term in the fluid equations for both ions and electrons.

  **0**: Off

  **1**: On

.. _255:
255 : 1/2 m v\ :sup:`3` n Term
  This will turn on/off the kinetic convection term for ions.

  **0**: Off

  **1**: On

.. _256:
256 : Prad Option
  This option turns on and off the radiative loss source term. The radiative loss term imposes a certain amount of power loss in each cell that could be due to radiative power losses. There are a couple of options supported.

  **0**: Off - No Radiative Losses

  **1**: On - Exponential Decay radiation source. The radiation falls off exponentially away from the target. The length of the decay source (`209`_), the characteristic decay distance of the exponential (`210`_), and the total integrated power radiated (`211`_) are all specified in parameters described earlier in the document.

    .. math::
      P_{rad}(s) = F_{rr} \times (P_{ae} + P_{ai}) \times e^{-s / \lambda_r}

  **2**: On - The power loss is described by the Garching Model (`212`_ - `215`_). The equations and their parameters are described in the parameter section.

  **3**: On - Radiative losses are proportional to PINQE. The multiplier is specified by this parameter (`249`_).

.. _257:
257 : Phelpi Option
  This option estimates the power loss for electrons due to hydrogenic ionization.

  **0**: Off - No power loss due to hydrogenic ionization.

  **1**: On - Analytic calculation of electron energy losses due to hydrogenic ionization. The analytic formula used to estimate this is the following:

    .. math::
      P_{helpi}(s) = H_{elpi} \times S_{iz}(s)
  
    .. math::
      H_{elpi} = 17.5 + (5.0 + 37.5 / T_e(s)) \times (1.0 + 0.25 / T_e(s)) \times log_{10}(10^{21} / n_e(s))

  **2**: On - PINQE On - PINQE used - Option 1 (ANALYTIC) is used on the seed plasma iteration. For subsequent iterations the electron energy losses as reported by PIN are used.

  **3**: On - PINQE On - PINQE used - Option 0 (OFF) is used on the seed plasma iteration. For subsequent iterations the electron energy losses as reported by PIN are used. 

.. _258:
258 : Pei Option
  Equipartition power term due to electron/ion energy transfer. This term has been found to be quite destabilizing at low temperatures because as one solves the equations moving upstream - the hotter species will get hotter due to this term and the cooler species will cool further. To be stable this term needs adequate compensating flows from other power sources.

  **0**: Off - the term is not calculated and is not included in the power balance equations.

  **1**: On - the term is calculated according to the formula in the DIVIMP Guide - modified by the Pei parameter described previously. (Note: There is a typographical error in the Guide in the Pei equation - eqn. 3.68: the exponent in the numerator should not be 3/2 but instead should be 1):

    .. math::
      P_{ei}(s) = P_{ei-cf} \times \frac{1.14 \times 10^{-32} n^2 (T-e - T_i)}{m_b T_e ^{3/2}}
  
    .. math::
      \lambda = \frac{1.5 \times 10^{13} T_e^{3/2}}{\sqrt{n_e}}

  **3**: Off - (but calculated) Pei is not applied to the solver but the values are calculated using the formula of option 1 and are then stored in the .lim or .SOL output file. 

.. _259:
259 : Pcx Option
  This option instructs the solver to include the designated ion energy loss/gain due to various processes in the ion energy balance equation.

  **0**: Off - No ion power term is used.

  **1**: On - Analytic CX An analytic formulation of ion CX energy losses only is included. (See Guide for more information)

    .. math::
      P_{CX}(S) = CEICF \times (1.5 \times T_i) \times S_{iz}(s)

  Option 2.0: On - PINQI On - For the seed plasma iteration, option 1 is used to calculate the ion power term (keep in mind the CEICF factor specified in the parameter list). On subsequent iterations the ion energy loss/gain from PIN is used.

  Option 3.0: On - PINQI On - For the seed plasma iteration, option 0 (OFF) is used for the ion power term. On subsequent iterations the ion energy loss/gain from PIN is used.

  Option 4.0: On - DIVIMP QI On - The seed plasma iteration is calculated using option 0 (OFF) for the ion power term. On subsequent iterations, DIVIMP uses the quantities returned from PIN to estimate the values of four of the components parts of the ion power source term. These are the Atomic Ionization, Molecular Ionization, Charge Exchange and Recombination. Each of these components to DIVIMP QI can be calculated with different options which are outlined below. Furthermore, each of the contributions has an associated cut-off temperature. If the background temperature in a cell is less than the cutoff then no contribution to the specific source term in that cell is included in the integration. (This term is also referred to as PINQID).

  Option 5.0: On - PINQI On - Internal PCX term (option 1) is on for seed plasma iteration. PINQI term is used for subsequent iterations EXCEPT that the PINQI term is clipped. Plasma heating by the PINQI term is not allowed. 

.. _260:
260 : PINQID - DIVIMP calculated Qi - Atomic Ionization
  This option specifies how the contribution to the ion power source term due to atomic ionization is to be calculated. The values for the densities of neutral hydrogen as well as neutral hydrogen molecules and the ionization source are all generated by PIN. When the temperature in the cell is less than a specified cutoff the contribution is set to zero.

  **0**: Off - no contribution due to atomic ionization

  **1**: On - Ionization from PIN.

    .. math::
      P_{atiz}(s) = \frac{3}{2} kT_{atom} \times \frac{n_H}{n_H + n_{H_2}}

  **2**: On - Ionization rates from ADAS

    .. math::
      P_{atiz}(s) = \frac{3}{2} kT_{atom} \times n_e \times n_H \times <\sigma \nu >_{atiz}

.. _261:
261 : PINQID - DIVIMP calculated Qi - Molecular Ionization
  This option specifies how the contribution to the ion power source term due to molecular ionization is to be calculated. The values for the densities of neutral hydrogen as well as neutral hydrogen molecules and the ionization source are all generated by PIN. When the temperature in the cell is less than a specified cutoff the contribution is set to zero.

  **0**: Off - no contribution due to molecular ionization

  **1**: On - Ionization from PIN - fixed 3 eV contribution/event.

    .. math::
      P_{mliz}(s) = 3.0 \times \frac{n_{H_2}}{n_H + N_{H_2}} \times S_{iz}(s)

  **2**: On - Ionization rates from ADAS - the atomic ionization rate is used as an approximation to the molecular ionization rate.

    .. math::
      P_{mliz}(s) = 3.0 \times n_e \times n_{H_2} \times <\sigma \nu >_{atiz}

.. _262:
262 : PINQID - DIVIMP calculated Qi - Recombination
  This option specifies how the contribution to the ion power source term due to recombination is to be calculated. When the temperature in the cell is less than a specified cutoff the contribution is set to zero.

  **0**: Off - no contribution due to recombination power terms.

  **1**: On - Recombination rates from ADAS.

    .. math::
      P_{rec}(s) = -\frac{3}{2} kT_{i} \times n_i \times n_e \times <\sigma \nu >_{rec}

.. _263:
262 : PINQID - DIVIMP calculated Qi - Charge Exchange
  This option specifies how the contribution to the ion power source term due to charge exchange is to be calculated. The values for the densities of neutral hydrogen and the ionization source are generated by PIN. When the temperature in the cell is less than a specified cutoff the contribution is set to zero.

  **0**: Off - no contribution due to charge exchange.

  **1**: On - Ionization from PIN. Reference Temperature is a specified parameter.

    .. math::
      P_{CX}(s) = \frac{3}{2} kT_{ref} \times R_{CXmult}(s) \times CEICF \times S_{iz}(s)

  **2**: On - Charge Exchange cross-section from analytic formula. The temperature differential is used for calculating contributions.

    .. math::
      P_{CX}(s) = \frac{3}{2}k(T_{atom} - T_i) \times n_e \times n_H \times <\sigma \nu >_{cx}
  
    .. math::
      E_{av} = \frac{3}{2}k \frac{T_{atom} + T_i}{2}
  
    .. math::
      <\sigma \nu >_{cx} = 10^{-13}\ for\ E_{av} > 1000\ eV
      <\sigma \nu >_{cx} = 10^{-14} \times E_{av}^{1/3}\ for\ E_{av} \le 1000\ eV

  **3**: On - Charge Exchange cross-section from ADAS. The temperature differential is used for calculating contributions. ADAS charge exchange cross-sections were considered unreliable at the time of writing.

    .. math::
      P_{CX}(s) = \frac{3}{2} k(T_{atom} - T_i) \times n_e \times n_H \times <\sigma \nu >_{cx}

  **4**: On - This option is identical to Option 2 - except that the contributions are limited to ion cooling only. (i.e. "-" ve - or power loss from ions - if Tatom is greater than Ti - then the contribution (r) 0 for that cell). 

.. _264:
264 : PP ElecLoss

  **0**: Private plasma electron power loss compensation term is OFF.

  **1**: Private Plasma electron power loss compensation term is ON. Electron power lost to each element of the private plasma target is removed from the main SOL rings in the corresponding position relative to the separatrix. This power loss is distrtibuted evenly to the Xpoint on the main SOL rings.

  **2**: Private plasma electron power loss compensation term is ON. Electron power lost to each element of the private plasma target is removed from main SOL rings in the corresponding position relative to the separatrix. This power loss is distrtibuted evenly over Smax * the power distribution parameter (`230`_) from each target on the main SOL rings.

.. _265:
265 : Switch: PP IonLoss

  **0**: Private plasma ion power loss compensation term is OFF.

  **1**: Private Plasma ion power loss compensation term is ON. Ion power lost to each element of the private plasma target is removed from the main SOL rings in the corresponding position relative to the separatrix. This power loss is distrtibuted evenly to the Xpoint on the main SOL rings.

  **2**: Private plasma ion power loss compensation term is ON. Ion power lost to each element of the private plasma target is removed from main SOL rings in the corresponding position relative to the separatrix. This power loss is distrtibuted evenly over Smax * the power distribution parameter (`230`_) from each target on the main SOL rings.

.. _266:
266 : Viscosity Option
  NOT IMPLEMENTED. This was associated with attempts to include parallel viscosity in the transport equations. At present it does nothing but should be specified with an Option of 0 to ensure that it is OFF.

.. _267:
267 : Momentum Loss Option
  This option specifies the form of the momentum loss or pressure term in the equations that are being solved. The pressure loss due to neutral interactions can have a significant effect on the solution. Furthermore, it tends to be a stabilizing factor within the solver since the pressure is an important term in the equation where imaginary values are usually found.

  **0**: Off - pressure loss is turned OFF.

  **1**: On - Rectangular On - rectangular/constant/flat momentum loss source from the target to specified cut off length, with a specified total integral value. (See Guide)

    .. math::
      S_{mom} = S_{mom0}\ for\ S \le F_L \times SMAX
  
    .. math::
      S_{mom} = 0\ for\ S > F_L \times SMAX
  
    .. math::
      S_{mom0} = \frac{P_{target}}{F_L \times SMAX} (\frac{1.0}{F_{fric}} - 1.0)

  **2**: On - Exponential On - Exponential decay momentum loss source - decaying away from the target with specified decay length, maximum length and magnitude specified by the F\ :sub:`fric` factor. (See Guide)

    .. math::
      S_{mom}(s) = S_{mom0} \times e^{-s / (F_\lambda \times SMAX)}\ for\ S \le F_L \times SMAX
  
    .. math::
      S_{mom}(s) = 0.0\ for\ S > F_L \times SMAX
  
    .. math::
      S_{mom0} = \frac{P_{target}}{F_L \times SMAX} (\frac{1.0}{F_{fric}} - 1.0) (\frac{1.0 - e^{-F_L \times SMAX}}{F_\lambda \times SMAX})

  **3**: On - Proportional On - Proportional to ionization source. The magnitude of the momentum loss is defined as in options 2 and 3. However, in this case it is distributed in proportion to the ionization source. The integral is performed over the half field line. This is proportional to the specified ionization source - either analytic or returned from a PIN run.

  .. math::
    S_{mom}(s) = S_{mom0} \times S_{iz}(s)

  .. math::
    S_{mom0} = \frac{P_{target}}{F_L \times SMAX} (\frac{1}{F_fric} - 1) (\frac{1}{\int S_{iz}(s')ds'})

  **4**: On - Proportional On - Proportional to the ionization source. The momentum loss in this case is calculated assuming that the ions are moving at the background velocity found at the given s-position by the solver. In addition it is proportional to the ionization source modified by a number of multiplicative factors defined earlier in the text. R\ :sub:`cxmult` is a function of Te and varies between 1.0 and 1500.0. The other multiplicative factor R\ :sub:`cx/iz` is a constant parameter that is usually set to 1.0 but can be used to increase the effectiveness of the momentum loss term over the entire range.

  .. math::
    S_{mom} = -m_bV_b(s) \times R_{cxmult}(T_e) \times R_{cx/iz} \times S_{iz}(s)

  **5**: Off Off - not used at present. This was an initial test option for PIN related momentum loss. It now sets momentum loss equal to zero if it is selected.

  **6**: On - PIN - Untested On - PIN - This option generates the seed plasma using Momentum Loss Option 0 (OFF) and then reads the momentum source from PIN (PINMP array). The PINMP array may be unreliable for statistical reasons so this option is not recommended unless some verification of the validity of the PINMP values has already been completed.

  **7**: On - PIN - Untested On - PIN - This option generates the seed plasma using Momentum Loss Option 1 (Rectangular) and then reads the momentum source from PIN (PINMP array). The PINMP array may be unreliable for statistical reasons so this option is not recommended unless some verification of the validity of the PINMP values has already been completed.

  **8**: On - PIN - Untested On - PIN - This option generates the seed plasma using Momentum Loss Option 2 (Exponential) and then reads the momentum source from PIN (PINMP array). The PINMP array may be unreliable for statistical reasons so this option is not recommended unless some verification of the validity of the PINMP values has already been completed.

  **9**: On - EDGE2D On - EDGE2D - This option is based on charge exchange momentum loss Cross-sections taken from the Edge2D/NIMBUS implementation. The option needs quantities from a PIN run so for the initial seed plasma calculation the momentum loss option is turned OFF. The rates are multiplied by the R\ :sub:`cx/iz` factor.

  Implemented in DIVIMP by Wojciech Fundamenski.

  **10**: On - EDGE2D On - EDGE2D - This option is similar to option 9. This option factors in the H\ :sub:`0` velocity distributions as returned by PIN.

.. _268:
268 : Iterative Mach Number Option
  This option instructs the solver to attempt to resolve any invalid solutions for density or temperature (e.g. imaginary density values or negative temperatures) by increasing the background flow mach number at the target. Following both the super and sub-sonic solutions to the equations yields a situation in which one follows the super-sonic branch from the target until it can make a smooth transition (at the point where the two solutions just touch) to the sub-sonic branch. The solver steps up the mach number until it finds a solution where there are no imaginary values encountered. It then steps back to the previous mach number that did not work, divides the mach number increment by 10 and starts iterating up from there until it again finds a solution that does not generate invalid results. It continues this process until it reaches the resolution limit for the mach number specified in the input and then stops. At this time, the mach solver finds that the solutions near the critical mach number are very dependent on the exact value of the target mach number. This instability can cause odd effects near the transition point in the plasma solution. If the mach solver is in use, it is important to peruse the results and check them for validity as well as checking the output files for any warnings.

  **0**: Off - the target Mach number is fixed at the input value - usually 1.0. The velocity error setting is used to deal with any imaginary values encountered.

  **1**: On - Target Mach number is changed as indicated to obtain a solution but the density at the target is held fixed. This effectively becomes a modification of the target particle flux.

  **2**: On - The target density is changed as well as the target mach number in such a way that the particle flux onto the target segment is conserved. This is the only option compatible with unnormalized PIN ionization and perpendicular flux correction option 2. Since most target data is actually based on flux measurements - it is believed that this is the better option to select if using mach number iteration.

  **3**: On Fixed Target MACH number option. Target mach number is NOT fixed to 1.0 but is defined to be a constant value calculated from the flow velocity at the target divided by the sound speed. For example, if the flow velocity happens to equal the sound speed the solver will use MACH=1.0. In this option, for a value other than 1.0 to be found - the target flow velocity will have been derived using alternate methods. Perhaps obtained from the ratio of the down flux and the target density. 

.. _269:
269 : Edge 2D Data Compatibility Option
  This is an attempt to improve compatibility between the SOL Option 22 solver and Edge2D for comparison purposes. When this option is active the solver starts at the middle of the first grid cell using data drawn directly from an Edge2D case. There are several different options due to the difficulty of precisely specifying the velocity at the middle of the first cell based on an Edge2D output. In addition, this option was implemented in the first place solely because of the difficulty in precisely specifying the target conditions applicable in an Edge2D case. As such, this option should be turned off when examining experimental data.

  **< 0**: On - EDGE2D compatibility options less than zero instruct the solver to run once using the EDGE2D background plasma solution for the SOL and then to use the absolute value of the compatibility option for subsequent iterations of SOL 22.

  **0**: Off - the solver works from the target with the specified target conditions determining the velocity and the sound speed at the target.

  **1**: On - Pressure On - The velocity at the starting point (middle of the first cell) is calculated from conservation of pressure from the actual target conditions, as reported by Edge2D, to the values reported at the middle of the first cell. This assumes that there is no pressure loss in the first half cell.

  **2**: On - EDGE2D Ghost On - The velocity at the starting point (middle of the first cell) is the average of the velocity at the cell faces (as reported by Edge2D) in the ghost (.g80) plasma background file.

  **3**: On - Parallel Flux On - The velocity at the starting point (middle of the first cell) is calculated by taking the parallel flux into and out of the cell at each cell boundary - obtaining an average flux for the cell and then dividing by the cell density to get an average velocity for the cell centre.

  **4**: Debug - EDGE2D Ghost On - Edge2D Data is read for first cell - the velocity is taken from the EDGE2D value for the first cell centre. SOL 22 is NOT run. The EDGE2D solution is used for the background plasma. This option is used only to produce detailed flux analyses for debugging purposes.

  *5**: Debug - Off Off - fluxes are calculated from the target - optionally taken from the EDGE2D solution. SOL 22 is NOT run. The EDGE2D solution is used for the background plasma. This option is used only to produce detailed flux analyses for debugging purposes.

  **6**: Debug - Pressure On - Edge2D data is read for first cell - EDGE2D pressure is matched at the first cell centre. SOL 22 is NOT run. The EDGE2D solution is used for the background plasma. This option is used only to produce detailed flux analyses for debugging purposes.

  **7**: Debug - Parallel Flux On - Edge2D data is read for first cell - the cell centre velocity is calculated from the cell boundary fluxes and cell density extracted from the EDGE2D solution. SOL 22 is NOT run. The EDGE2D solution is used for the background plasma. This option is used only to produce detailed flux analyses for debugging purposes.

  **8**: On - Down Flux On - Solver will run from the middle of the first cell. Cell centre velocity is obtained by averaging the EDGE2D fluxes into the cell and dividing by the density. The fluxes used in this option are NOT the EDGE2D fluxes from the GHOST file. These are the EDGE2D DOWN fluxes and are read in from an auxiliary input file. These values will be extracted from the EDGE2D down flux listing. If EDGE2D TARGET OPTION 5 is also selected then the down power fluxes as well as particle fluxes will also be used.

  **9**: On - Down Flux (Knot) On - Solver runs from the middle of a SPECIFIED cell. The starting knot index is specified by the Start Knot Index Value described above. Edge2D data is required for the entire ring. EDGE2D DOWN fluxes are extracted from an auxiliary file - the EDGE2D DOWN flux listing. The starting velocity at the cell centre is determined by averaging the cell face down fluxes and dividing by the density. If EDGE2D TARGET OPTION 5 is also specified then the solver will use both the down particle fluxes and down power fluxes. The cells between the target and the cell where the solver begins are filled with values determined by the FILL OPTION described above.

.. _270:
270 : Power Distribution Option
  Many of the older DIVIMP SOL options implicitly assume that all of the power onto the targets plus any other sources must enter at the mid-point between targets and be carried all the way along the SOL by the various transport mechanisms. This option adds the ability to distribute the required input power (target flows and in some options volume sources as well) over various lengths of the whole or half-ring. This models the expected reduction in power transported as the equations are solved towards the mid-point. The effect of this option is to reduce some of the "peakiness" towards the mid-point of the temperature solutions seen in SOL option 12, 13 and 22 (when this option is turned OFF). However, when used in conjunction with PIN ionization and perpendicular flux corrections, it can occur that the convection terms in the heat equation at the mid-point end up carrying all of the heat flux. In some cases, the convection terms carry significantly more than the heat flux required to satisfy the target and volume power sinks on the ring. In order to compensate for this, the conduction term is forced to carry heat in the opposite direction in order to satisfy the conservation equation. This situation can result in temperatures dropping (downward temperature gradient) towards the mid-point as the conduction term tries to counteract the convection term. In particular, when using PINQI (the ion energy term from PIN), the integration of the volume power terms can show a net power gain by the ions, which when combined with a reduced requirement for power flow as one approaches the mid-point, can result in negative ion temperatures being encountered by the solver. This situation is still under investigation. However, errors of this type are flagged and printed in both the output data file and in the .lim file. Also, unless the power distribution option explicitly distributes a specific power source, it will be treated as if the power to supply it was coming in at the mid-point between the targets.

  **Option 0**: Off Off - All power is assumed to come in at the top (mid-point on the ring between targets.)

  **Option 1**: On - Target - Half Ring On - The total TARGET power flux (for both electrons and ions) is evenly distributed over the entire 1/2 ring. This means that the power required to be carried by each species falls linearly from the target power flux at s = 0 to 0.0 at s = SMAX/2. There is no adjustment for volume power terms. These are implicitly assumed to be supplied by a flow from the mid-plane.

  **Option 2**: On - Target - X-point On - The total TARGET power flow is evenly distributed from an S-position approximately equivalent to the X-point for each ring to the mid-point. For the region from the X-point to the target the equations are identical to option 0 where all of the target power flow is being transported by the usual mechanisms.

  **Option 3**: On - Major Radius On - This is identical to option 1 except that the power losses have been modified by a Major Radius correction. This is for use ONLY when the major radius correction option has been selected. In addition, many newer options have not been designed to work with the Major Radius corrected version of the solver and can not be expected to work correctly.

  **Option 4**: On - Targets - Whole Ring On - The power flux to both targets is added together and then evenly distributed over the entire ring. This forces the power being carried by each species to ramp linearly from the value at one target to the value at the other target. This will likely result in power being transported across the mid-point. Keep in mind that the sign of the power flux at the two targets is different because the velocity at the two targets have opposite signs.

  **Option 5**: On - Target + PIN - Half On - The TARGET power flow and power from both the PINQE and PINQI volume terms are added together and distributed over the half ring. This option is equivalent to option 1 for the seed plasma solution and otherwise only works with PIN when either or both of PINQE and PINQI are specified for source power loss terms in the Phelpi and Pcx options mentioned above.

  **Option 6**: On - Target +PIN - Whole On - This is the same as option 4 except that it also includes the PIN based power terms if available and in use.

  **Option 7**: On - Target - Dist On - the target power flow is distributed from a specified position F1 × SMAX to the midpoint of the ring. The factor F1 is described in the parameter section above.

  **Option 8**: On - Target + PIN - Dist On - This is the same as option 7 except that it also includes the PIN based power terms if they are available and in use.

  **Option 9**: On - Target - Dist2 On - the Target power flow is evenly distributed between two given positions on the field-line. (F1 × SMAX to F2 × SMAX).

  **Option 10**: On - Target + PIN - Dist2 On - This is the same as option 9 except that it also includes the PIN based power terms if they are available and in use.

  **Option 11**: On- Target - Whole - Dist On - This is the same as option 4 - the total target power outflux is summed and distributed evenly over the ring starting at a distance F1 ( SMAX from each target.

.. _271:
271 : Private Plasma Power Distribution 
  This option allows a different power distribution option to be specified in the Private Plasma. This feature may be necessary because the physics of the power influx in the private plasma region may be quite different from the processes and behaviour that dominate the heat flux into the main SOL rings.

  **-1**: Special Set the value of this switch equal to the value used for the general power option.

  **0-11**: Off/On - These options are identical to those described above.

.. _272:
272 : Gamma Perp Option
  This is a perpendicular flux correction option. It is used to ensure that the particle balance for each ring is maintained. When using an unnormalized ionization source, with or without recombination particle sources, a condition of over or under-ionization on a flux tube may be encountered. If left uncorrected this will result in a constant background drift velocity beyond the end of the ionization source because there are no additional particle sources or sinks. Furthermore, the only time the velocity will be exactly zero is when the target sink is equal to the ionization source. This is obviously both an unacceptable and unphysical solution since it is clear that cross-field particle sources and sinks will, in steady state, ensure that the sum of sources and sinks on the field line is zero. This option implements a very simple version of a cross-field source - any particle excess or deficit (after considering the target fluxes and the ionization source) is then compensated for by a cross-field source/sink term distributed over the entire field line. This ensures that the velocity will cross through zero at least once and also results in a more realistic evolution of the particle source. This option is only effective with unnormalized sources - typically a PIN result. It has no effect on sources that are normalized to the target flux as there is no excess/deficit that needs to be compensated for in these cases. In general, this option will have no effect on the seed plasma iteration since the seed plasma ionization source will usually be normalized.

  **0**: Off Off - no additional cross-field source or sink is used in the calculations.

  **1**: On - Half-Ring On - The flux correction is calculated for each half ring independently. This causes the flux (and thus also the velocity) to fall to zero at the mid-point of the ring and ensures particle conservation on each half ring.

  **2**: On - Whole Ring On - The flux correction is calculated for the entire flux tube from target to target by summing ionization and target flux for the entire ring and distributing the resulting difference uniformly over the entire ring. The velocity will cross zero somewhere on the ring and in the case of ionization exceeding target fluxes - may cross more than once.

  **3**: On - Whole Ring - N On - The flux correction is calculated for the entire flux tube from target to target by summing ionization and target flux for the entire ring. Additional cross-field sources for under-ionized rings are applied using an evenly distributed perpendicular flux. Cross-field sinks for over-ionized rings are applied using a flux that is proportional to the density in each cell calculated on the previous iteration of the solver.

  **4**: On - Half Ring - N On - This is the same as option 3 except that the particle source excess/deficit is calculated for only one half ring at a time. The distribution of the flux uses the strategy outlined in option 3.

  **5**: On - Half Ring + Rect On - The net flux along the field line goes to zero at the midpoint with a specified fraction of the required perpendicular flux being distributed uniformly over the half-ring with the remainder being distributed uniformly between specified start and end points on the half-ring. The fraction in the rectangular source is specified by the Compound Gperp fraction described previously. The rectangular region is specified from GperpF1 × SMAX to GperpF2 × SMAX from the target.

  **6**: On - Whole Ring + Rect On - This is the same as option 5 except that the net flux goes to zero for the entire flux tube considered as a whole with a specified fraction of the required cross-field flux being included in the flux in two specified regions. One region near each target. The rest of the flux is included uniformly over the entire flux tube. The fraction in the rectangular source is specified by the Compound Gperp fraction described previously. The rectangular region is specified from GperpF1 × SMAX to GperpF2 × SMAX from both targets.

  **7**: On-Whole Ring-Gradient On - The net flux over the entire ring goes to zero with the excess or deficit distributed proportional to the second gradient of the density. This gradient is derived from either an EDGE2D solution or the SOL 22 solution from a previous iteration. This option is changed to a uniform distribution for any rings where the total positive or negative contribution for the ring exceeds five times the integrated value of the whole. The reason for this is that such case are usually unstable and do not produce useful results.

  **8**: On-Whole Ring-Absolute On - The net cross-field flux over the entire ring goes to zero. A cross-field component of the flux is calculated using a second gradient of the density of a previous iteration and a fixed value for the diffusion coefficient. Any remaining excess or deficit after this term is included is imposed as a uniform source or sink over the entire ring.

.. _273:
273 : Private Plasma Gamma Perp Option
  A different Gperp option can be specified for the Private Plasma from that used for the main SOL. However, these options will share the same parameter values if the Gperp option requires them. The options available are the same as those for the regular Gperp option.

.. _274:
274 : Extra Perpendicular Source and Sink Option

  **0**: Extra perpendicular flux term is OFF.

  **1**: Extra perpendicular flux term is ON. An extra source and sink are superimposed on the flux tube. This source and sink exactly cancel but will affect the flow pattern on the flux tube.

  Source and Sink Strength = (total target flux on ring) * Fstr     (`227`_)

  The source is imposed over the region: Smax * [F1,F2]     (`228`_)

  The sink is imposed over the region: Smax * [F3,F4]     (`229`_)

.. _275:
275 : Major Radius Option
  This is an attempt to restructure the entire solver to work using equations that have been adapted to a varying value of Major Radius. Extensive comparisons have been made between this and the standard OSM methods that do not include the major radius effect. Differences are minimal and it is recommended that this option be left turned OFF. Furthermore, a number of the newer features of both the ionization sources and power sources will not work correctly in combination with this option. Some of the options are intended to NOT generate correct major radius solutions but to instead explore the magnitude of the effects of these changes on the solutions.

  **0**: Off OFF - NORMAL operation. This is the recommended setting.

  **1**: On - Target Correction On - All target fluxes are adjusted by Rtarg/R0

  **2**: On - Source Correction On - All ionization sources are adjusted by Rcell/R0

  **3**: On - Source Correction On - All ionization sources are adjusted by R0/Rcell

  **4**: On - General Correction On - Generalized R-correction to both ionization and target fluxes. Both of these quantities must be adjusted in order to correctly include major radius effects.

.. _276:
276 : Core Flux Source
  NOT IMPLEMENTED. The purpose of this option was originally related to compensating the ionization source for ionization occurring inside the innermost ring on the grid. However, this was then superseded by the Perpendicular Flux correction option that ensured particle conservation for each ring. It became unnecessary to try to add to the ion source any particles that would have been ionized within the central core escape region of the grid because the flux of these particles would already be included in the Perpendicular Flux option. In addition, the distribution of such a source would not be better defined using a separate option than would be possible using the available perpendicular flux options. This entry in the input file may be redefined in future releases of DIVIMP

  **0**: Off - This option will do nothing - no matter what value is used. 

.. _277:
277 : Recombination Source Option

  **0**: Off - No recombination particle source is added to the target flux and ionization in the calculation of the spatial particle fluxes and ring particle balance.

  **1**: On - DIVIMP/PIN On - The recombination particle source, as calculated by DIVIMP, is added to the calculation of the net particle flux and is included in the calculation of the particle balance on each ring. The temperature specified for the recombination cut-off limit results in this minimum value being used in the calculation of the recombination rates.

  **2**: On - Special Edge2D On - The recombination particle source calculated from the input Edge2d background is used in the solver. The option specifying the formulae to be used to calculate the recombination is described later in the text.

.. _278:
278 : Smoothing Option:
  This option will apply a smoothing algorithm to the density, Te and Ti profiles across the mid-point of the ring and thus supply a smoother looking background without the peakiness at the mid-point associated with the non-uniform power distribution usually employed by the solver. This is however, a completely ad hoc adjustment to the background and does not reflect any physics at all. The sole purpose is to obtain a smoother looking background plasma. However, since the solutions from each end of the ring would be expected to meet somewhere in the middle, this may not be a bad approximation, just one that is difficult to justify.

  **0**: Off - No smoothing is done

  **1**: On - Smoothing by averaging over a number of adjacent cells is applied to smooth the peaks that usually occur at the mid-point of the rings. 

  **2**: On - Agreement at the midpoints is FORCED by "pivoting" each half-ring solution to meet exactly halfway between where they disagree. Values at the target are left unchanged, while values further away from the target are modified relatively more (think of the solution of each half-ring like two halves of a drawbridge that don't meet at the middle; the tips of the drawbridge halves will move the most so they meet at middle, and the portions that are at the road will not move at all). This algorithm is best used *after one has done everything they can* to obtain agreement with the other SOL 22 options since it will modify the gradients along the entire ring. Implemented by Shawn Zamperini.

.. _279:
279 : Detached Plasma Prescription Option
  This option allows either the inner or outer half-rings to be specified using the detached plasma specification instead of solved using SOL22.

  **0**: Off - Detached Plasma Prescription is OFF

  **1**: On - Detached Plasma Model is used for the first half-ring starting at (IK=1). (OUTER target for JET)

  **2**: On - Detached Plasma Model is used for the second half-ring (IK=NKS(IR)) . (INNER target for JET).

.. _280:
280 : Error Correction Level
  This switch turns on the solver error recovery mechanism. For some rings under certain conditions it can be impossible for the solver to arrive at a consistent and correct solution with the selected options. Usually this will occur when volume power source terms have been turned on and it will occasionally result in negative ion temperatures being found by the solver. Under these circumstances and if the minimum solver temperature option has not been specified, the solver sometimes can not find a valid solution for the ring. If this option is turned on the solver will restart the solution of the ring using a modified set of options that is more conservative than those originally specified. The solver will start at the specified error recovery level and work it's way down from the highest numbered option to the lowest. Each level of error recovery will include all of the error recovery actions taken in all of the higher levels unless the option specifically states otherwise. The eventual solution may be very simple and thus may not reflect the physics that could be included. However, the solution will also not include pathological values like negative temperatures that will interfere with the behaviour of the rest of the particle transport or with the behaviour of the other codes with which DIVIMP interacts. Whenever, this option is turned on, the print-out for SOL 22 should be checked and all rings which demonstrate an error condition should be closely examined to assure validity.

  **0**: Off - No error correction - if the solver dies on a particular ring the contents of the background plasma values will be those from the previous iteration OR the target values if this is the first iteration. The actual contents depend on the plasma decay option selected and the default values that are loaded into the arrays in the "plasma" subroutine.

  **1**: On - Conduction Only On - The ring is solved using an analytic ionization source and only the conduction term - all other options and switches are turned OFF.

  **2**: On - As Level 3.0 + All convective terms are turned off.

  **3**: On - As Level 4.0 + All power terms are turned off.

  **4**: On - As Level 5.0 + turn off convection terms proportional to v\ :sup:`2`

  **5**: On - As Level 6.0 + all power enters at the flux tube at the mid-point.

  **6**: On - As Level 7.0 + replace whole ring uniform particle balance with the equivalent half ring options.

  **7**: On - As Level 8.0 + use half-ring uniform power terms instead of whole ring.

  **8**: On - As Level 9.0 + use only the cooling portion of PINQI if the PINQI term is included.

  **9**: On - As level 10.0 + Use a uniform Gperp source instead of one proportional to d\ :sup:`2` n/dr\ :sup:`2` if this option is active.

  **10**: On - Turn off equipartition if it is on. Pei option is set to zero.

.. _281:
281 : Automatic DEFAULT error correction
  This option reads in a list of ring numbers that are to be solved using error correction level 1.0 default only. It has been found in some cases that only specific rings will experience trouble in the solver and these will always experience problems with the given set of selected options and switches. To save time in the solver for these rings- one can bypass the initial solution and instruct the solver to use the default error recovery initially for this list of specified rings. This option is usually used ONLY for pathological cases. This option is now redundant except for exceptionally pathological cases since the overall error correction level would eventually try to solve the ring using these options in any case if it could not find a solution using more complex options. 
