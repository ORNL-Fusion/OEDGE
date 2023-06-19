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

============ ========================================================================
  `200 Tags`_   
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
    <\sigma \nu >_{cx} = 10^{-13}\ \mathdefault{for}\ E_{av} > 1000\ \mathdefault{eV}
    <\sigma \nu >_{cx} = 10^{-14} \times E_{av}^{1/3}\ \mathdefault{for}\ E_{av} \le 1000\ \mathdefault{eV}

  **3**: On - Charge Exchange cross-section from ADAS. The temperature differential is used for calculating contributions. ADAS charge exchange cross-sections were considered unreliable at the time of writing.

  .. math::
    P_{CX}(s) = \frac{3}{2} k(T_{atom} - T_i) \times n_e \times n_H \times <\sigma \nu >_{cx}

  **4**: On - This option is identical to Option 2 - except that the contributions are limited to ion cooling only. (i.e. "-" ve - or power loss from ions - if Tatom is greater than Ti - then the contribution (r) 0 for that cell). 




