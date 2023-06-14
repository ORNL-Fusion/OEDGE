module mod_sol22_output



  implicit none

  !public :: prerrdesc


contains




  subroutine echosolorg(spts,npts)
    use mod_solparams
    use mod_solswitch
    use mod_solcommon

    !     This subroutine echoes the input values to standard out - it also
    !     prints the additional calculated values - and is followed by the
    !     tabular output of s,te,ti,n and v at each point S ... which would
    !     be suitable for plotting on a spreadsheet or may be plotted by
    !     calling GHOST routines if the graph option is set.

    !     include 'solparams'
    !     include 'solswitch'
    !     include 'solcommon'

    implicit none
    integer npts,i

    real*8 spts(mxspts)


    CHARACTER  COMENT*77

    ! jdemod assigned in mod_sol22_interface
    !s22_title='DIVIMP S22'
    call prbs
    call prs(s22_title)
    call pris('Ring number: ',ringnum)

    call prbs
    CALL PQS (' Denisty at the target:             ', n0)
    CALL PQS (' Electron temperature at the target:', te0)
    CALL PQS (' Ion temperature at the target:     ', ti0)
    CALL PQS (' Plasma ion mass (amu)              ', mb)
    CALL PQS (' Flow velocity at the target        ', v0)
    CALL PQS (' Electron heat coefficient (K0e)    ', k0e)

    !     Mach number options

    CALL PQS (' Ion heat coefficient      (K0i)    ', k0i)
    call prbs
    if (switch(swmach).eq.0.0) then
       call prs ('Iterative Mach number solution - OFF')
       CALL PQS (' Imposed mach number at the target  ', m0)
    elseif (switch(swmach).eq.1.or.switch(swmach).eq.2) then
       call prs ('Iterative Mach number solution - ON ')
       CALL PQS (' Initial mach number at the target  ',origm0)
       CALL PQS (' Final mach number at the target    ',m0)
       call pqs (' Initial Target Density             ',netarg)
       call pqs (' Final  Target Density              ',nefinal)
       CALL PQS (' Delta m0 for first iteration       ',deltam0)
       CALL PQS (' Ultimate resolution of m0          ',m0res)
       call pris(' Number of iterations to obtain m0  ',miter)
       call pqs (' Supersonic transition point (m)    ',lastiters)
    elseif (switch(swmach).eq.3.0) then
       call prs ('Iterative Mach number solution - OFF')
       CALL PrS (' MACH numbers imposed at target vary.')

    endif
    call prbs
    call prs ('Other calculated quantities:        ')

    call prbs
    call pqs (' Gammae - electron power coefficient', gammae)
    call pqs (' Gammai - ion power coefficient     ', gammai)
    call pqs (' Gamma correction factor in gammai  ', gamcor)
    call pqs (' Electron power to target           ', pae)
    call pqs (' Ion power to target                ', pai)
    call pqs (' Pressure at infinity(w.mach effect)', pinf)

    call pqs (' Value of Smax for ring (m)         ',ringlen)
    call prbs
    call prs ('Set of activated switches           ')

    call prbs
    if (switch(swcond).eq.0.0) then
       call prs(' 5/2 nv * kT     term  is OFF        ')
    else
       call prs(' 5/2 nv * kT     term  is ON         ')

    endif
    if (switch(swconv).eq.0.0) then
       call prs(' 1/2 m v^3 * n   term  is OFF        ')
    else
       call prs(' 1/2 m v^3 * n   term  is ON         ')

    endif
    if (switch(swprad).eq.0.0) then
       call prs(' Prad            term  is OFF        ')
    else
       call prs(' Prad            term  is ON         ')

    endif
    if (switch(swphelp).eq.0.0) then
       call prs(' Phelpi          term  is OFF        ')
    elseif (switch(swphelp).eq.1.0) then
       call prs(' Internal Phelpi term  is ON         ')
    elseif (switch(swphelp).eq.2.0) then
       call prs(' PINQE           term  is ON         ')

    endif
    if (switch(swpei).eq.0.0) then
       call prs(' Pei             term  is OFF        ')
    elseif (switch(swpei).eq.1.0) then
       call prs(' Pei             term  is ON         ')
    elseif (switch(swpei).eq.3.0) then
       call prs(' Pei             term  is Calculated (not used)')

    endif
    if (switch(swpcx).eq.0.0) then
       call prs(' Pcx             term  is OFF        ')
    elseif (switch(swpcx).eq.1.0) then
       call prs(' Internal Pcx    term  is ON         ')
       call pqs(' CX power coefficeint CEICF          ',ceicf)
    elseif (switch(swpcx).eq.2.0) then
       call prs(' PINQI           term  is ON         ')

    endif
    if (switch(swppelec).eq.0.0) then
       call prs(' PPelec          term  is OFF        ')
    elseif (switch(swppelec).eq.1.0) then
       call prs(' PPelec          term  is ON to Xpoint')
    elseif (switch(swppelec).eq.2.0) then
       call pqs(' PPelec          term  is ON to SMAX *',pp_pow_dist)

    endif
    if (switch(swppion).eq.0.0) then
       call prs(' PPion           term  is OFF        ')
    elseif (switch(swppion).eq.1.0) then
       call prs(' PPion           term  is ON to Xpoint')
    elseif (switch(swppion).eq.2.0) then
       call pqs(' PPion           term  is ON to SMAX *',pp_pow_dist)

    endif
    if (switch(swppress).eq.0.0) then
       call prs(' PP Press        term  is OFF        ')
    elseif (switch(swppress).eq.1.0) then
       call prs(' PP Press        term  is ON to Xpoint')
    elseif (switch(swppress).eq.2.0) then
       call pqs(' PP Press        term  is ON to SMAX *',pp_pow_dist)



    endif
    if (switch(swvisc1).eq.0.0) then
       call prs('       Density Viscosity correction term'//'  is NOT IMPLEMENTED')


    endif
    if (switch(swnmom).eq.0.0) then
       call prs('       Neutral BG Momentum loss term      is OFF')
    elseif (switch(swnmom).eq.1.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Smom = Smom0    S < L * Smax        ')
       call prs('              = 0        S > L * Smax        ')
       call prs('         Smom0= Pt/(L * Smax) * (1/Ffric -1)')
       call pqs('         Ffric= ',ffric)
       call pqs('         L    = ',lenmom)
    elseif (switch(swnmom).eq.2.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Smom = Smom0 * exp ( -S / (Lm * Smax))  for'//' S < L * Smax')
       call prs('              = 0        S > L * Smax        ')
       call prs('         Smom0= Pt/(Lm * Smax)*(1/Ffric -1) / expf')
       call prs('         expf = (1 - exp  (-(L*Smax)/(Lm * Smax)))')
       call pqs('         Ffric= ',ffric)
       call pqs('         Lm   = ',lammom)
       call pqs('         L    = ',lenmom)
    elseif (switch(swnmom).eq.3.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Smom = fact * Siz(S) / INT (0 to L) '//'[ Siz(S'') dS'' ] ')
       call prs('         fact = Pt *  (1/Ffric -1) ')
       call pqs('         Ffric= ',ffric)
       call prs('         L    = Smax / 2.0')
    elseif (switch(swnmom).eq.4.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Smom = -m * Vb * Rcxiz *  Siz(S) ')
       call pqs('         Rcxiz= ',rcxmom)
    elseif (switch(swnmom).eq.5.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Smom is read directly from NIMBUS')
    elseif (switch(swnmom).eq.6.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Neutral BG Momentun loss term is taken')
       call prs('         from PIN results except for the first')
       call prs('         iteration, in which case: ')
       call prs('         Neutral BG Momentum loss term is OFF')
    elseif (switch(swnmom).eq.7.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Neutral BG Momentun loss term is taken')
       call prs('         from PIN results except for the first')
       call prs('         iteration, in which case: ')
       call prs('         Initial neutral BG Momentum loss term is:')
       call prs('         Smom = Smom0    S < L * Smax        ')
       call prs('              = 0        S > L * Smax        ')
       call prs('         Smom0= Pt/(L * Smax) * (1/Ffric -1)')
       call pqs('         Ffric= ',ffric)
       call pqs('         L    = ',lenmom)
    elseif (switch(swnmom).eq.8.0) then
       call prs('       Neutral BG Momentum loss term      is ON')
       call prs('         Neutral BG Momentun loss term is taken')
       call prs('         from PIN results except for the first')
       call prs('         iteration, in which case: ')
       call prs('         Initial neutral BG Momentum loss term is:')
       call prs('         Smom = Smom0 * exp ( -S / (Lm * Smax))  for'//' S < L * Smax')
       call prs('              = 0        S > L * Smax        ')
       call prs('         Smom0= Pt/(Lm * Smax)*(1/Ffric -1) / expf')
       call prs('         expf = (1 - exp(-(L * Smax)/(Lm * Smax)))')
       call pqs('         Ffric= ',ffric)
       call pqs('         Lm   = ',lammom)
       call pqs('         L    = ',lenmom)

       !     More switches

    endif
    if (switch(swe2d).eq.0.0) then
       call prs('       Edge2D BG compatibility is OFF')
    elseif (switch(swe2d).eq.1.0) then
       call prs('       Edge2D BG compatibility is ON')
       call prs('       - Edge2D data is read for first cell')
    elseif (switch(swe2d).eq.2.0) then
       call prs('       Edge2D BG compatibility is ON')
       call prs('       - Edge2D data is read for first cell')
       call prs('       - Edge2D Flux is matched at first')
       call prs('         cell centre')

    endif
    if (switch(swpow).eq.0.0) then
       call prs('       Distributed Power Influx is OFF')
    elseif (switch(swpow).eq.1.0) then
       call prs('       Distributed Power Influx is ON')
       call prs('       - Target Power flux is distribited linearly')
       call prs('         along the entire ring to the midpoint')
    elseif (switch(swpow).eq.2.0) then
       call prs('       Distributed Power Influx is ON')
       call prs('       - Target Power flux is distribited linearly')
       call prs('         from the X-point region to the midpoint')

    endif
    if (switch(swgperp).eq.0.0) then
       call prs('       Perpendicular Flux Correction is OFF')
    elseif (switch(swgperp).eq.1.0) then
       call prs('       Perpendicular Flux is ON')
       call prs('       - Net flux along field line ->0 at the'//' midpoint')
       call prs('         using an evenly distributed'//' perpendicular flux.')

    endif
    if (switch(swgperpp).eq.0.0) then
       call prs('       PP Perpendicular Flux Correction is OFF')
    elseif (switch(swgperpp).eq.1.0) then
       call prs('       PP Perpendicular Flux is ON')
       call prs('       - Net flux along field line ->0 at the'//' midpoint')
       call prs('         using an evenly distributed'//' perpendicular flux.')

    endif
    if (recfrac.ne.1.0) then
       call prs('       Recycling Source Fraction is ON')
       call pqs('       - Recycling source fraction = ',recfrac)

       !      if (switch(swcore).eq.0.0) then
       !         call prs('       Cross-field Core Ionization Source is OFF')
       !      elseif (switch(swcore).eq.1.0) then
       !         call prs('       Cross-field Core Ionization Source is ON')
       !         call pqs('       - Core source fraction : ',corefrac)
       !         call prs('       - Core source flux is distribited linearly')
       !         call prs('         along the entire ring to the midpoint')
       !      elseif (switch(swcore).eq.2.0) then
       !         call prs('       Cross-field Core Ionization Source is ON')
       !         call pqs('       - Core source fraction : ',corefrac)
       !         call prs('       - Core source flux is distribited linearly')
       !         call prs('         from the X-point region to the midpoint')
       !      endif



    endif
    if (switch(swmajr).eq.0.0) then
       call prs(' Major Radius Factor is OFF')
    elseif (switch(swmajr).eq.1.0) then
       call prs(' Major Radius Factor is ON')
       call prs('  - All Target Flux adjusted by Rtarg/R0')
    elseif (switch(swmajr).eq.2.0) then
       call prs(' Major Radius Factor is ON')
       call prs('  - All ionization adjusted by Rcell/R0')
    elseif (switch(swmajr).eq.3.0) then
       call prs(' Major Radius Factor is ON')
       call prs('  - All ionization adjusted by R0/Rcell')
    elseif (switch(swmajr).eq.4.0) then
       call prs(' Major Radius Factor is ON')
       call prs('  - Generalized R-correction applied')

       !     Ionization Source characteristics

    endif

    if (switch(swion).eq.0.0) then
       call prbs
       call prs (' Ionization Source Characteristics  S(s)')
       call prs ('       Exponential Decay : ')

       call prbs
       call pqs ('       Length of ionization source        ', lensfi)

       !      call pqs (' Base strength of ionization source ', s0)

       call pqs ('       Decay length of ionization source  ', lams)

    elseif (switch(swion).eq.1) then
       call prbs
       call prs (' Ionization Source Characteristics  S(s)')
       call prbs
       call prs ('       Data for Ionization Source is returned')
       call prs ('       from a PIN/NIMBUS run and is linearly')
       call prs ('       interpolated for values between grid')
       call prs ('       points. The ionization source is')
       call prs ('       integrated over the half-ring and set')
       call prs ('       equal to the flow to the target NoVo')

       call prs ('       as a Normalization factor.')

       call pqs ('       Normalization factor = ',fnorm)

    elseif (switch(swion).eq.2) then
       call prbs
       call prs (' Ionization Source Characteristics  S(s)')
       call prbs
       call prs ('       Data for Ionization Source is returned')
       call prs ('       from a PIN/NIMBUS run and is linearly')
       call prs ('       interpolated for values between grid')
       call prs ('       points. The data is used AS-IS and is')
       call prs ('       NOT Normalized to the target flux for')

       call prs ('       the ring.')

    elseif (switch(swion).eq.3.0) then
       call prs(' Triangular Ionization Source:  S(s)')
       call pqs('          Extending from   :',lensst)
       call pqs('                    to     :',lensfi)

       call prs('          Integral of source normalized'//' to ring target flux.')

    elseif (switch(swion).eq.3.0) then
       call prs(' Rectangular (Flat) Ionization Source: S(s)')
       call pqs('          Extending from   :',lensst)
       call pqs('                    to     :',lensfi)

       call prs('          Integral of source normalized'//' to ring target flux.')

    endif
    if (switch(swprad).eq.1.0) then
       call prbs
       call prs ('Radiation Source Characteristics  Prad')

       call prbs
       call pqs (' Length of radiation source         ', lenr)
       call pqs (' Decay length of radiation  source  ', lamr)
       call pqs (' Source strength fraction (frr)     ', frr)
       call pqs (' Strength of radiation source       ', prad0)
    elseif (switch(swprad).eq.6.0) then
       call prbs
       call prs ('Radiation Source Characteristics  Prad')

       call prbs
       call pqs (' End Length of radiation source     ', lenr)
       call pqs (' Start length of radiation  source  ', lamr)
       call pqs (' Source strength fraction (frr)     ', frr)

       call pqs (' Strength of radiation source       ', prad0)

       !     Output controls

    endif
    call prbs
    if (graph.eq.0) then
       call prs ('Plotting is turned OFF')
    elseif (graph.eq.1) then
       call prs ('Plotting is turned ON')

    endif
    if (graphaux.eq.0) then
       call prs ('Auxilliary plots and tables are turned OFF')
    elseif (graphaux.eq.1) then
       call prs ('Auxilliary plots and tables are turned ON')

    endif
    if (graphvel.eq.0) then
       call prs ('Velocity plots and tables are turned OFF')
    elseif (graphvel.eq.1) then
       call prs ('Velocity plots and tables are turned ON')

    endif
    call prbs
    if (graphran.eq.0.0) then
       call prs ('Close up plots are turned OFF')
    elseif (graphvel.eq.1) then
       call prs ('Close up plots are turned ON')
       call prr ('Range for close up plots is 0.0  to ',graphran)

    endif
    call prbs
    call pris ('Initial number of Runge-Kutta steps between'//' each S:', ndiv)

    call prbs
    call prs ('Set of S-values for Background calculations')
    call prbs
    do i = startn, npts
       call pqs('    S = ', spts(i))
    end do
    call prbs
    call prs ('End of Specification Section')
    call prbs

    call prs ('Results : ')
    return


  end subroutine echosolorg


  SUBROUTINE PRS(STRING)
    use mod_io_units
    !  *********************************************************************
    !  *  PRS:  PRINTS A CHARACTER STRING                                  *
    !  *********************************************************************
    implicit none
    CHARACTER STRING*(*)
    WRITE (s22out,'(1X,A)') STRING
    RETURN

  END SUBROUTINE PRS


  SUBROUTINE PRBS
    use mod_io_units
    !  *********************************************************************
    !  *  PRBS:  PRINTS A BLANK LINE                                       *
    !  *********************************************************************
    implicit none
    WRITE (s22out,'(1X)')
    RETURN


  END SUBROUTINE PRBS


  SUBROUTINE PRIS (NAME, I)
    use mod_io_units
    !  *********************************************************************
    !  *  PRIS:  PRINTS AN INTEGER                                         *
    !  *********************************************************************
    implicit none
    CHARACTER NAME*(*)
    INTEGER   I
    WRITE (s22out,'(1X,A,I7)') NAME,I
    RETURN


  END SUBROUTINE PRIS


  SUBROUTINE PQS (NAME, R)
    use mod_io_units
    !  *********************************************************************
    !  *  PQS:  PRINTS A REAL*8 NUMBER                                     *
    !  *********************************************************************
    implicit none
    CHARACTER NAME*(*)
    real*8      R
    IF (ABS(R).LT.0.1.OR.ABS(R).GE.1000.0) THEN
       WRitE (s22out,'(1X,A,1P,G15.8)') NAME,R
    ELSEIF (ABS(R).LT.1.0) THEN
       WRITE (s22out,'(1X,A,F15.8)') NAME,R
    ELSE
       WRITE (s22out,'(1X,A,F15.8)') NAME,R
    ENDIF
    RETURN

  END SUBROUTINE PQS


  subroutine prerrdesc(errlevel)
    implicit none

    !     PRERRDESC: Prints a description of the ERROR correction
    !                protocol used by the solver.
    real errlevel

    call prb

    CALL PRC('  ERROR CORRECTION IN SOLVER IS TURNED ON')

    CALL PRI('    - INITIAL ERROR CORRECTION LEVEL SET TO : ',NINT(ERRLEVEL))
    CALL PRC('    - ERROR CORRECTION MOVES FROM HIGHEST TO'//' LOWEST LEVEL')

    CALL PRC('      UNTIL SOLUTION IS OBTAINED.')

    call prc('    - ALL CHANGES ARE CUMULATIVE')

    !              errlevel=10- turn off equipartition if it was activated
    !              errlevel=9 - use uniform particles instead of d2n/dr2
    !              errlevel=8 - use only PINQI cooling contributions
    !              errlevel=7 - use 1/2 ring uniform power instead of whole
    !              errlevel=6 - use 1/2 ring uniform power + 1/2 ring particles
    !              errlevel=5 - 1/2 ring uniform particles + power at top
    !              errlevel=4 - 5 + turn off v^2 convection term
    !              errlevel=3 - 4 + no power terms
    !              errlevel=2 - 3 + no convective terms
    !              errlevel=1 - Conduction ONLY

    call prb

    CALL PRC('  LISTING OF ERROR CORRECTION LEVELS:')
    CALL PRC('  10 - TURN OFF EQUIPARTITION IF IT IS ON')
    CALL PRC('   9 - REPLACE DENSITY GRADIENT DEPENDENT'//' CROSS-FIELD TERM WITH UNIFORM')
    call prc('   8 - NO HEATING BY PINQI IS ALLOWED.')
    CALL PRC('   7 - REPLACE WHOLE RING UNIFORM POWER WITH'//' HALF RING UNIFORM.')
    CALL PRC('   6 - HALF RING UNIFORM POWER AND HALF RING UNIFORM'//' PARTICLES')
    CALL PRC('   5 - HALF RING UNIFORM PARTICLES'//' AND POWER IN AT TOP')
    call prc('   4 - 1/2 M V^3 CONVECTIVE TERM TURNED OFF')
    CALL PRC('   3 - ALL ADDITIONAL POWER TERMS TURNED OFF')
    call prc('   2 - ALL CONVECTIVE TERMS TURNED OFF')
    CALL PRC('   1 - CONDUCTION ONLY - ANALYTIC IONIZATION ONLY.')
    call prb
    return
  end subroutine prerrdesc


end module mod_sol22_output
