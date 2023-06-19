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
============ ============================

============ ============================
  `B Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `C Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `D Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `E Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `F Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `G Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `H Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `I Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `N Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `P Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `Q Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `S Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `T Tags`_   
-----------------------------------------
  `B01`_      String
============ ============================

============ ============================
  `Z Tags`_  (SOL 29) 
-----------------------------------------
  `B01`_      String
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

B Tags
------

C Tags
------

D Tags
------

F Tags
------

G Tags
------

H Tags
------

I Tags
------

N Tags
------

P Tags
------

Q Tags
------

S Tags
------

T Tags
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
    \gamma_i = 2.5 * 1/2 M^2 (1.0 + T_{e0} / T_{i0}) + \gamma_{i-corr}

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
  This parameter is used with power distribution (ref here?) options 7, 8 and 9. It specifies the starting position, as a function of SMAX for the ring, of the region where the power being carried by the electrons and ions should start reducing. This can be used for distributing the target power flux and any volume power terms if present. The simplest assumption in SOL 22 occurs when there are no heat sources or sinks along the ring and the power reaching the target is carried by the electrons and ions all the way from the mid-point ("top"). This is an adequate first order assumption - but in fact the power carried to the target or being radiated along the length of the ring would have likely entered the ring with some spatial distribution and not at the top. The power distribution options address the question of distributing this power loading.

.. _223:
223 : Distributed Power End Position
  This specifies the position (times SMAX for the ring) at which all the power has been input. (See previous entry and notes on power distribution options.)

.. _224:
224 : Compound Gperp - Fraction of Gperp in Rectangular Distribution
  Gperp is the amount of additional cross-field source or sink that is necessary to force particle conservation on an individual flux-tube. This additional source or sink can be distributed in a variety of ways - determined by the Gperp option described below. This parameter specifies the fraction of the required amount that will be added in a rectangular distribution in a region defined by the following two parameters. The remainder is added uniformly over the length of the flux tube. This parameter is used with Gperp options 5 and 6 (ref).

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
  If the Start Knot Index for EDGE2D Compatibility Option (ref)...

  This parameter specifies the starting knot number for EDGE2D compatibility option 9 (ref). SOL option 22 will start solving at this knot using the fluxes and background plasma conditions taken from EGE2D at this point. The following option describes what is done with the cells between the target and the starting point of the solver.

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
  PINQID charge exchange option 1 (ref), described in the section on SOL 22 switches, requires a reference temperature in order to calculate the power loss through the ion charge exchange channel.

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
  This option applies to PRAD option 3 (ref). Prad option 3 defines the radiative loss term to be a multiple of the PINQE (electron energy loss term that is calculated by PIN. This parameter specifies the value of the multiplier.

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
    S_{mom0} = \frac{P_target}{F_L \times SMAX} (\frac{1.0}{F_{fric}} - 1.0) (\frac{1.0 - e^{-F_L \times SMAX}}{F_\lambda \times SMAX})

  **3**: On - Proportional On - Proportional to ionization source. The magnitude of the momentum loss is defined as in options 2 and 3. However, in this case it is distributed in proportion to the ionization source. The integral is performed over the half field line. This is proportional to the specified ionization source - either analytic or returned from a PIN run.

  .. math::
    S_{mom}(s) = S_{mom0 \times S_{iz}(s)

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
