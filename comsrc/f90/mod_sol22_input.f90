module mod_sol22_input


  implicit none
  !
  ! jdemod - I am creating this data module to start replacing common blocks.
  !
  ! Modules perform the same task but allow for both data sharing and initialization.
  ! In addition, executable code related to variable setup can also be included. 
  !
  ! Ultimately it would be good for all the unstructured input to be organized and added to files 
  ! like this one. 
  !

  integer,public :: debug_sol22 = 0        ! 284 - SOL22 debug flag - default 0 = off
  integer,public :: debug_sol22_ir = 1     ! 285 - SOL22 debug - ring for hi res plasma - default = 1
  integer,public :: debug_sol22_ikopt = 1  ! 286 - SOL22 debug - ikopt (ring end) for detailed plasma - default = 1



contains


  subroutine readsol(ierr)
    !use mod_params
    use mod_io
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

    call rdi(forcet,.TRUE.,0,.true.,1,         'force te=ti'    ,ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER FORCET')

    CALL RDQ(initm0,.TRUE.,0.0d0,.FALSE.,0.0d0,'target mach num',IERR)
    CALL RDQ(deltam0,.TRUE.,0.0d0,.FALSE.,0.0d0,'delta mach num',IERR)

    !     ------------------------------------------------------------------

    !     Ionization Source

    CALL RDQ(m0res,.TRUE.,0.0d0,.FALSE.,0.0d0,'Resolution in m0',IERR)
    CALL RDI(lensind,.TRUE.,0,.TRUE.,1,     'ion source abs/rel',IERR)
    CALL RDQ(lensst,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src start',  IERR)
    CALL RDQ(lensfi,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src finish', IERR)
    CALL RDQ(lams,.TRUE. ,0.0d0,.FALSE.,0.0d0,'ion decay len   ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER IONIZATION SOURCE')
    !     Radiation source

    CALL RDQ(lenri,.TRUE.,0.0d0,.FALSE.,0.0d0,'rad source len  ',IERR)
    CALL RDQ(lamri,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad decay len  ',IERR)
    CALL RDQ(frri ,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad power mult ',IERR)
    CALL RDQ(alfimp,.TRUE.,0.0d0,.FALSE.,0.0d0,'nimp/ne ratio  ',IERR)
    CALL RDQ(talimp,.TRUE. ,0.0d0,.FALSE.,0.0d0,'base Temp     ',IERR)

    call pr_trace('MOD_SOL22_INPUT','AFTER RADIATION SOURCE')
    !     Miscellaneous

    CALL RDQ(ex1imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 1  ',IERR)
    CALL RDQ(ex2imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 2  ',IERR)
    CALL RDQ(gamcor,.false.,0.0d0,.FALSE.,0.0d0,'i power corr. ',IERR)

    CALL RDQ(gamecor,.false.,0.0d0,.FALSE.,0.0d0,'e power corr.',IERR)
    CALL RDQ(ceicf,.TRUE. ,0.0d0,.FALSE.,0.0d0, 'CX power frac ',IERR)

    CALL RDQ(recfrac,.TRUE.,0.0d0,.TRUE.,1.0d0,'Recycle frac ',IERR)
    call rdq(peicf,.true.  ,0.0d0,.false.,0.0d0,'Pei Correction',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MISC')
    !     Power distribution

    call rdi(velsw,.true.  ,0,.true.,3,       'Vel Error Switch',ierr)
    call rdq(spowbeg,.true.,0.0d0,.true.,0.5d0,'Power Dist Beg',ierr)

    !     Gperp Distribution function

    call rdq(spowlen,.true.,0.0d0,.true.,0.5d0,'Power Dist Len',ierr)
    call rdq(gperpfrac,.true.,0.0d0,.true.,1.0d0,'Part Dist Frac',ierr)
    call rdq(gperpbegf,.true.,0.0d0,.true.,0.5d0,'Part Dist Beg',ierr)

    !     Extra Gperp source/sink - start and end positions.

    call rdq(gperpendf,.true.,0.0d0,.true.,0.5d0,'Part Dist Len',ierr)
    call rdq(gextra_mult,.true.,0.0d0,.false.,0.0d0,'Gextra flux mult',ierr)
    call rdq2(gextra_src_start,gextra_src_stop,.true.,0.0d0,.true.,1.0d0, 'Gextra SRC start/stop',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER GPERP')

    !     Field line length fraction for distributing the Private plasma
    !     electron and ion power loads.

    call rdq2(gextra_sink_start,gextra_sink_stop,.true.,0.0d0,.true.,1.0d0,'Gextra SINK start/stop',ierr)

    !     IK index for edge2d compatibility option 9

    call rdq(pp_pow_dist,.true.,0.0d0,.false.,0.0d0,'PP Pow Dist',ierr)
    call rdi(ike2d,.true.,1,.false.,0,'IK start Index for E2D-9',ierr)

    !     Cutoff temperature for PINQID term

    call rdi(fillopt,.true.,0,.true.,3,'Gap fill option - E2D-9',ierr)
    call rdq(tcutqe,.true.,0.0d0,.false.,0.0d0,'Cut T-PINQE',ierr)
    call rdq(tcutatiz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDATIZ',ierr)
    call rdq(tcutmliz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDMLIZ',ierr)
    call rdq(tcutrec,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDREC',ierr)
    call rdq(tcutcx,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDCX',ierr)
    call rdq(trefcx,.true.,0.0d0,.false.,0.0d0,'REF T- QIDCX 1',ierr)
    call rdq(tmin,.false.,0.0d0,.false.,0.0d0,'Min. Allowed T',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER PINQID')

    !     Momentum Source

    call rdq(dropfrac,.true.,0.0d0,.true.,1.0d0,'Allowed T-drop',ierr)
    call rdq(smom_mult,.false.,0.0d0,.false.,0.0d0,'Mom.Loss Multiplier',ierr)
    call rdq(ffric,.false.,0.0d0,.false.,0.0d0,'Mom.loss frac  ',ierr)
    CALL RDQ(lenmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom source len',IERR)
    CALL RDQ(lammom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom decay len ',IERR)
    CALL RDQ(rcxmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'cx/iz ratio   ',IERR)
    call rdq(tcxmom,.TRUE.,1.0001d0,.FALSE.,0.0d0,'T for CXmult',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER MOMENTUM')

    !     Source term multipliers

    call rdq(tcxcut,.true.,0.0d0,.false.,0.0d0,'Cut T -  CXmult',ierr)
    call rdq(qesrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE mult',ierr)

    !     call rdq(qisrc_mult,.false.,0.0d0,.false.,0.0d0,'PINQI mult',ierr)

    !     ------------------------------------------------------------------

    call rdq(radsrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE based PRAD mult',ierr)

    CALL RDI(ndiv,.TRUE., 1,.FALSE., 0,'NUMBER OF STEPS    ',IERR)

    call pr_trace('MOD_SOL22_INPUT','BEFORE SWITCHES')

    !     Read in switches 0.0 is off, 1.0 is on.

    call rdr(switch(swion),.true. ,0.0,.false.,0.0, 'alt ion',ierr)
    call rdr(switch(swioni),.true. ,0.0,.false.,0.0,'init ion',ierr)

    call rdr(switch(swionp),.false.,0.0,.false.,0.0,'pp ion',ierr)
    if (switch(swionp).eq.-1.0) then
       if (switch(swion).eq.1.0.or.switch(swion).eq.2.0.or.switch(swion).eq.8.0) then
          switch(swionp) = switch(swioni)
       else
          switch(swionp) = switch(swion)
       endif

    endif
    call rdr(switch(swcond),.true.,0.0,.false.,0.0, 'cond sw',ierr)
    call rdr(switch(swconv),.true.,0.0,.false.,0.0, 'conv sw',ierr)
    call rdr(switch(swprad),.true.,0.0,.false.,0.0, 'prad sw',ierr)
    call rdr(switch(swphelp),.true.,0.0,.false.,0.0,'phelp sw',ierr)
    call rdr(switch(swpei),.true.,0.0,.false.,0.0,  'pei sw ',ierr)

    call pr_trace('MOD_SOL22_INPUT','AFTER BASIC SWITCHES')

    !     PINQID switches

    call rdr(switch(swpcx),.true.,0.0,.false.,0.0,  'pcx sw ',ierr)
    call rdr(switch(swqidatiz),.true.,0.0,.false.,0.0,'atiz sw',ierr)
    call rdr(switch(swqidmliz),.true.,0.0,.false.,0.0,'mliz sw',ierr)
    call rdr(switch(swqidrec),.true.,0.0,.false.,0.0,'rec sw',ierr)

    call rdr(switch(swqidcx),.true.,0.0,.false.,0.0,'cx sw',ierr)
    call rdr(switch(swppelec),.true.,0.0,.false.,0.0,'pp elec sw',ierr)


    !     jdemod -
    !     switch(swppress) is read in using tag 283 of unstructured input
    !     default value is 0 or OFF

    !      call rdr(switch(swppress),.true.,0.0,.false.,0.0,'pp power sw',ierr)


    call rdr(switch(swppion),.true.,0.0,.false.,0.0,'pp ion sw',ierr)
    call rdr(switch(swvisc1),.true.,0.0,.true.,0.0,'visc1 sw',ierr)
    call rdr(switch(swnmom),.true.,0.0,.false.,0.0, 'N mom sw',ierr)

    call pr_trace('MOD_SOL22_INPUT','HALFWAY THROUGH SWITCHES')


    !     Read in Edge2d compatibility switch and the subsequent values of
    !     ne, Te and Ti at the centre point of the first cell

    call rdr(switch(swmach),.true.,0.0,.false.,0.0, 'mach sw',ierr)
    call rdr(switch(swe2d),.false.,0.0,.false.,0.0, 'e2d sw',ierr)
    call rdr(switch(swpow),.true.,0.0,.false.,0.0, 'power sw',ierr)

    call rdr(switch(swpowp),.false.,0.0,.false.,0.0, 'pp pow',ierr)

    if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)
    call rdr(switch(swgperp),.true.,0.0,.false.,0.0,'GamPerp',ierr)

    call rdr(switch(swgperpp),.false.,0.0,.false.,0.0,'GamPerpP',ierr)

    !     Extra Gperp source/sink term

    if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)

    call rdr(switch(swextra),.true.,0.0,.false.,0.0,'GP Src/Sink',ierr)
    call rdr(switch(swmajr),.true.,0.0,.false.,0.0,'MajorRad',ierr)
    call rdr(switch(swcore),.true.,0.0,.false.,0.0,'Core Src',ierr)
    call rdr(switch(swrecom),.true.,0.0,.false.,0.0,'Recomb.',ierr)
    call rdr(switch(swsmooth),.true.,0.0,.false.,0.0,'Smooth',ierr)
    call rdr(switch(swdetach),.true.,0.0,.false.,0.0,'Detach',ierr)
    call rdr(switch(swerror),.true.,0.0,.false.,0.0,'ERROR',ierr)

    call pr_trace('MOD_SOL22_INPUT','DONE SWITCHES')


    !CALL RDRARN(deflist,ndef,mxspts,0.0,real(maxnrs),.FALSE.,0.0,MACHHI,2,'DEFAULT SOLVER DATA',IERR)
    !
    ! jdemod - The max value should the maximum number of rings but to avoid dependency on params
    !          module I am removing this constraint. 
    CALL RDRARN(deflist,ndef,mxspts,0.0,SOL22_HI,.FALSE.,0.0,SOL22_MACHHI,2,'DEFAULT SOLVER DATA',IERR)

    !     Default plots to off - this can be changed in the calcsol_interface
    !     routine in solascv1.f

    graph = 0
    graphaux = 0
    graphvel = 0

    !      CALL RDI(graph,.TRUE., 0,.TRUE., 1,'GRAPH OPTION       ',IERR)
    !      CALL RDI(graphaux,.TRUE.,0,.TRUE.,1,'AUX GRAPH OPTION  ',IERR)
    !      CALL RDI(graphvel,.TRUE.,0,.TRUE.,1,'VEL GRAPH OPTION  ',IERR)
    !      call rdr(graphran,.true.,0.0,.false.,0.0,'CXMAX VALUE  ',ierr)


    graphran = 0.0

    
    call pr_trace('MOD_SOL22_INPUT','END OF READSOL')


    return

  end subroutine readsol

  subroutine verify_sol22
    use mod_io
    use debug_options
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none


    
    if (switch(swionp).eq.-1.0) then
       if (switch(swion).eq.1.0.or.switch(swion).eq.2.0.or.switch(swion).eq.8.0) then
          switch(swionp) = switch(swioni)
       else
          switch(swionp) = switch(swion)
       endif
    endif

    if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)

    !     Extra Gperp source/sink term
    if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)


    !     Default plots to off - this can be changed in the calcsol_interface
    !     routine in solascv1.f

    graph = 0
    graphaux = 0
    graphvel = 0

    !      CALL RDI(graph,.TRUE., 0,.TRUE., 1,'GRAPH OPTION       ',IERR)
    !      CALL RDI(graphaux,.TRUE.,0,.TRUE.,1,'AUX GRAPH OPTION  ',IERR)
    !      CALL RDI(graphvel,.TRUE.,0,.TRUE.,1,'VEL GRAPH OPTION  ',IERR)
    !      call rdr(graphran,.true.,0.0,.false.,0.0,'CXMAX VALUE  ',ierr)


    graphran = 0.0

    
  end subroutine verify_sol22
    
  
  subroutine sol22_initialize_unstructured_input
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none

    !
    !  Initialization for TAG series 2
    !
    !
    ! TAG:  +201    Initialization
    !
    ! Sample Input
    !
    !  '+201    Force Te=Ti through SOL 22  0=off  1=on         '     0
    !
    ! Input call documentation
    !
    !           'force te=ti'    
    !
    forcet = 0
    !
    ! TAG:  +202    Initialization
    !
    ! Sample Input
    !
    !  '+202    Imposed mach number at the target               '    1.0
    !
    ! Input call documentation
    !
    !  'target mach num'
    !
    initm0 = 1.0
    !
    ! TAG:  +203    Initialization
    !
    ! Comments from input
    !
    !     ------------------------------------------------------------------
    !     Ionization Source
    !
    ! Sample Input
    !
    !  '+203    Delta mach number for initial iterative solution'    0.1
    !
    ! Input call documentation
    !
    !  'delta mach num'
    !
    deltam0 = 0.1
    !
    ! TAG:  +204    Initialization
    !
    ! Sample Input
    !
    !  '+204    Maximum resolution in calculation of m0         ' 0.00001
    !
    ! Input call documentation
    !
    !  'Resolution in m0'
    !
    m0res = 0.00001
    !
    ! TAG:  +205    Initialization
    !
    ! Sample Input
    !
    !  '+205    Ionization Source Lengths  0=Absolute 1=Relative'     1
    !
    ! Input call documentation
    !
    !       'ion source abs/rel'
    !
    lensind = 1
    !
    ! TAG:  +206    Initialization
    !
    ! Sample Input
    !
    !  '+206    Start of Ionization Source (for supported opts) '    0.0
    !
    ! Input call documentation
    !
    !  'ion src start'
    !
    lensst = 0.0
    !
    ! TAG:  +207    Initialization
    !
    ! Sample Input
    !
    !  '+207    End or Length of Ionization Source              '    0.3
    !
    ! Input call documentation
    !
    !  'ion src finish'
    !
    lensfi = 0.3
    !
    ! TAG:  +208    Initialization
    !
    ! Sample Input
    !
    !  '+208    Decay length of ionization source               '    0.03
    !
    ! Input call documentation
    !
    !  'ion decay len   '
    !
    lams = 0.03
    !
    ! TAG:  +209    Initialization
    !
    ! Sample Input
    !
    !  '+209    Length of radiation source                      '    5.0
    !
    ! Input call documentation
    !
    !  'rad source len  '
    !
    lenri = 5.0
    !
    ! TAG:  +210    Initialization
    !
    ! Sample Input
    !
    !  '+210    Decay length of radiation  source               '    0.5
    !
    ! Input call documentation
    !
    !  'rad decay len  '
    !
    lamri = 0.5
    !
    ! TAG:  +211    Initialization
    !
    ! Sample Input
    !
    !  '+211    Source strength fraction (frr)                  '    1.0
    !
    ! Input call documentation
    !
    !  'rad power mult '
    !
    frri  = 1.0
    !
    ! TAG:  +212    Initialization
    !
    ! Sample Input
    !
    !  '+212    Garching Model: Alpha = ni/ne ratio             '    1.0
    !
    ! Input call documentation
    !
    !  'nimp/ne ratio  '
    !
    alfimp = 1.0
    !
    ! TAG:  +213    Initialization
    !
    ! Sample Input
    !
    !  '+213    Garching Model: Tbase = Tratio denominator      '   15.0
    !
    ! Input call documentation
    !
    !  'base Temp     '
    !
    talimp = 15.0
    !
    ! TAG:  +214    Initialization
    !
    ! Sample Input
    !
    !  '+214    Garching Model: Exponent 1                      '    1.5
    !
    ! Input call documentation
    !
    !  'exponent 1  '
    !
    ex1imp  = 1.5
    !
    ! TAG:  +215    Initialization
    !
    ! Sample Input
    !
    !  '+215    Garching Model: Exponent 2                      '   -3.0
    !
    ! Input call documentation
    !
    !  'exponent 2  '
    !
    ex2imp  = -3.0
    !
    ! TAG:  +216    Initialization
    !
    ! Sample Input
    !
    !  '+216    Gamma correction factor in gammai               '    0.0
    !
    ! Input call documentation
    !
    !  'i power corr. '
    !
    gamcor = 0.0
    !
    ! TAG:  +217    Initialization
    !
    ! Sample Input
    !
    !  '+217    Gamma correction factor in gammae               '    0.0
    !
    ! Input call documentation
    !
    !  'e power corr.'
    !
    gamecor = 0.0
    !
    ! TAG:  +218    Initialization
    !
    ! Sample Input
    !
    !  '+218    CX power coefficeint CEICF                      '    1.0
    !
    ! Input call documentation
    !
    !   'CX power frac '
    !
    ceicf = 1.0
    !
    ! TAG:  +219    Initialization
    !
    ! Sample Input
    !
    !  '+219    Recycling  source fraction                      '    1.0
    !
    ! Input call documentation
    !
    !  'Recycle frac '
    !
    recfrac = 1.0
    !
    ! TAG:  +220    Initialization
    !
    ! Sample Input
    !
    !  '+220    Pei Power Transfer Correction Factor            '    1.0
    !
    ! Input call documentation
    !
    !  'Pei Correction'
    !
    peicf = 1.0
    !
    ! TAG:  +221    Initialization
    !
    ! Sample Input
    !
    !  '+221    Velocity Error Switch       0=Cs   1=const      '     1
    !
    ! Input call documentation
    !
    !         'Vel Error Switch'
    !
    velsw = 1
    !
    ! TAG:  +222    Initialization
    !
    ! Comments from input
    !
    !     Gperp Distribution function
    !
    ! Sample Input
    !
    !  '+222    Distributed Power Start position * SMAX         '    0.10
    !
    ! Input call documentation
    !
    !  'Power Dist Beg'
    !
    spowbeg = 0.10
    !
    ! TAG:  +223    Initialization
    !
    ! Sample Input
    !
    !  '+223    Distributed Power End   position * SMAX         '    0.50
    !
    ! Input call documentation
    !
    !  'Power Dist Len'
    !
    spowlen = 0.50
    !
    ! TAG:  +224    Initialization
    !
    ! Sample Input
    !
    !  '+224    Distributed GPERP particle Fraction- non-uniform'    0.8
    !
    ! Input call documentation
    !
    !  'Part Dist Frac'
    !
    gperpfrac = 0.8
    !
    ! TAG:  +225    Initialization
    !
    ! Sample Input
    !
    !  '+225    Distributed GPERP Start position * SMAX         '    0.0
    !
    ! Input call documentation
    !
    !  'Part Dist Beg'
    !
    gperpbegf = 0.0
    !
    ! TAG:  +226    Initialization
    !
    ! Sample Input
    !
    !  '+226    Distributed GPERP End   position * SMAX         '    0.1
    !
    ! Input call documentation
    !
    !  'Part Dist Len'
    !
    gperpendf = 0.1
    !
    ! TAG:  +227    Initialization
    !
    ! Sample Input
    !
    !  '+227    Gextra Source strength - Target flux multiplier '    0.1
    !
    ! Input call documentation
    !
    !  'Gextra flux mult'
    !
    gextra_mult = 0.1
    !
    ! TAG:  +228    Initialization
    !
    ! Sample Input
    !
    !  '+228    Gextra Source Start/Stop * SMAX                 ' 0.2   0.35
    !
    ! Input call documentation
    !
    !   'Gextra SRC start/stop'
    !
    gextra_src_start = 0.2
    gextra_src_stop = 0.35
    !
    ! TAG:  +229    Initialization
    !
    ! Sample Input
    !
    !  '+229    Gextra Sink   Start/Stop * SMAX                 ' 0.65  0.8
    !
    ! Input call documentation
    !
    !  'Gextra SINK start/stop'
    !
    gextra_sink_start = 0.65
    gextra_sink_stop = 0.8
    !
    ! TAG:  +230    Initialization
    !
    ! Sample Input
    !
    !  '+230    PP target power loss redistribution range *SMAX '    0.1
    !
    ! Input call documentation
    !
    !  'PP Pow Dist'
    !
    pp_pow_dist = 0.1
    !
    ! TAG:  +231    Initialization
    !
    ! Sample Input
    !
    !  '+231    Knot Start Index for E2D Option 9               '     8
    !
    ! Input call documentation
    !
    !  'IK start Index for E2D-9'
    !
    ike2d = 8
    !
    ! TAG:  +232    Initialization
    !
    ! Sample Input
    !
    !  '+232    Plasma Fill option for missing knots - E2D Opt 9'     1
    !
    ! Input call documentation
    !
    !  'Gap fill option - E2D-9'
    !
    fillopt = 1
    !
    ! TAG:  +233    Initialization
    !
    ! Sample Input
    !
    !  '+233    Qe Term - Temperature Cutoff (eV)               '    0.0
    !
    ! Input call documentation
    !
    !  'Cut T-PINQE'
    !
    tcutqe = 0.0
    !
    ! TAG:  +234    Initialization
    !
    ! Sample Input
    !
    !  '+234    PINQID - Atomic Ionization    - T cutoff (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'Cut T-QIDATIZ'
    !
    tcutatiz = 0.0
    !
    ! TAG:  +235    Initialization
    !
    ! Sample Input
    !
    !  '+235    PINQID - Molecular Ionization - T cutoff (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'Cut T-QIDMLIZ'
    !
    tcutmliz = 0.0
    !
    ! TAG:  +236    Initialization
    !
    ! Sample Input
    !
    !  '+236    PINQID - Recombination        - T cutoff (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'Cut T- QIDREC'
    !
    tcutrec = 0.0
    !
    ! TAG:  +237    Initialization
    !
    ! Sample Input
    !
    !  '+237    Qi Term/PINQID-Charge Exchange- T cutoff (eV)   '    0.0
    !
    ! Input call documentation
    !
    !  'Cut T- QIDCX'
    !
    tcutcx = 0.0
    !
    ! TAG:  +238    Initialization
    !
    ! Sample Input
    !
    !  '+238    PINQID - CX option 1 - Reference T - (eV)       '    1.0
    !
    ! Input call documentation
    !
    !  'REF T- QIDCX 1'
    !
    trefcx = 1.0
    !
    ! TAG:  +239    Initialization
    !
    ! Sample Input
    !
    !  '+239    Minimum Temperature allowed in Solver (spec<0)  '    0.1
    !
    ! Input call documentation
    !
    !  'Min. Allowed T'
    !
    tmin = 0.1
    !
    ! TAG:  +240    Initialization
    !
    ! Sample Input
    !
    !  '+240    Minimum T allowed as fraction of Tmax reached   '    0.5
    !
    ! Input call documentation
    !
    !  'Allowed T-drop'
    !
    dropfrac = 0.5
    !
    ! TAG:  +241    Initialization
    !
    ! Sample Input
    !
    !  '+241    Momentum loss term multiplier   (Usually 1.0)   '    1.0
    !
    ! Input call documentation
    !
    !  'Mom.Loss Multiplier'
    !
    smom_mult = 1.0
    !
    ! TAG:  +242    Initialization
    !
    ! Sample Input
    !
    !  '+242    Friction factor for Momentum loss formula       '    0.2
    !
    ! Input call documentation
    !
    !  'Mom.loss frac  '
    !
    ffric = 0.2
    !
    ! TAG:  +243    Initialization
    !
    ! Sample Input
    !
    !  '+243    Length of the Momentum loss region * Smax       '    0.1
    !
    ! Input call documentation
    !
    !  'mom source len'
    !
    lenmom = 0.1
    !
    ! TAG:  +244    Initialization
    !
    ! Sample Input
    !
    !  '+244    Decay length of the Momentum loss region * Smax '    0.02
    !
    ! Input call documentation
    !
    !  'mom decay len '
    !
    lammom = 0.02
    !
    ! TAG:  +245    Initialization
    !
    ! Sample Input
    !
    !  '+245    Ratio of CX to IZ events (fixed)                '    1.0
    !
    ! Input call documentation
    !
    !  'cx/iz ratio   '
    !
    rcxmom = 1.0
    !
    ! TAG:  +246    Initialization
    !
    ! Sample Input
    !
    !  '+246    Te cutoff for increased CX multiplier (eV)      '    5.0
    !
    ! Input call documentation
    !
    !  'T for CXmult'
    !
    tcxmom = 5.0
    !
    ! TAG:  +247    Initialization
    !
    ! Sample Input
    !
    !  '+247    Te lower limit cutoff for CX multiplier (eV)    '    1.0
    !
    ! Input call documentation
    !
    !  'Cut T -  CXmult'
    !
    tcxcut = 1.0
    !
    ! TAG:  +248    Initialization
    !
    ! Sample Input
    !
    !  '+248    PINQE multiplier                                '    1.0
    !
    ! Input call documentation
    !
    !  'PINQE mult'
    !
    qesrc_mult = 1.0
    !
    ! TAG:  +249    Initialization
    !
    ! Sample Input
    !
    !  '+249    PRAD option 3 multiplier (x PINQE)              '    0.5
    !
    ! Input call documentation
    !
    !  'PINQE based PRAD mult'
    !
    radsrc_mult = 0.5
    !
    ! TAG:  +250    Initialization
    !
    ! Sample Input
    !
    !  '+250    Initial number of stages for Runge Kutta steps  '    100
    !
    ! Input call documentation
    !
    !  'NUMBER OF STEPS    '
    !
    ndiv = 100
    !
    ! TAG:  +251    Initialization
    !
    ! Sample Input
    !
    !  '+251    Switch: Ionization Opt : 0.0-exp 1.0+ - others  '    2.0
    !
    ! Input call documentation
    !
    !   'alt ion'
    !
    switch(swion) = 2.0
    !
    ! TAG:  +252    Initialization
    !
    ! Sample Input
    !
    !  '+252    Switch: Initial IonOpt : 0.0-exp 3.0+ - others  '    0.0
    !
    ! Input call documentation
    !
    !  'init ion'
    !
    switch(swioni) = 0.0
    !
    ! TAG:  +253    Initialization
    !
    ! Sample Input
    !
    !  '+253    Switch: PPlasma IonOpt : 0.0-exp 3.0+ - others  '   -1.0
    !
    ! Input call documentation
    !
    !  'pp ion'
    !
    switch(swionp) = -1.0
    !
    ! TAG:  +254    Initialization
    !
    ! Sample Input
    !
    !  '+254    Switch: 5/2 nv * kT    : 0.0-off 1.0-on         '    1.0
    !
    ! Input call documentation
    !
    !   'cond sw'
    !
    switch(swcond) = 1.0
    !
    ! TAG:  +255    Initialization
    !
    ! Sample Input
    !
    !  '+255    Switch: 1/2 m v^3 * n  : 0.0-off 1.0-on         '    1.0
    !
    ! Input call documentation
    !
    !   'conv sw'
    !
    switch(swconv) = 1.0
    !
    ! TAG:  +256    Initialization
    !
    ! Sample Input
    !
    !  '+256    Switch: Prad           : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !   'prad sw'
    !
    switch(swprad) = 0.0
    !
    ! TAG:  +257    Initialization
    !
    ! Sample Input
    !
    !  '+257    Switch: Phelpi         : 0.0-off 1.0-on         '    2.0
    !
    ! Input call documentation
    !
    !  'phelp sw'
    !
    switch(swphelp) = 2.0
    !
    ! TAG:  +258    Initialization
    !
    ! Sample Input
    !
    !  '+258    Switch: Pei            : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !    'pei sw '
    !
    switch(swpei) = 0.0
    !
    ! TAG:  +259    Initialization
    !
    ! Sample Input
    !
    !  '+259    Switch: Pcx            : 0.0-off 1.0-on         '    5.0
    !
    ! Input call documentation
    !
    !    'pcx sw '
    !
    switch(swpcx) = 5.0
    !
    ! TAG:  +260    Initialization
    !
    ! Sample Input
    !
    !  '+260    SUB-switch: Pcx Opt 4  : PINQID- Atomic Ioniz.  '    1.0
    !
    ! Input call documentation
    !
    !  'atiz sw'
    !
    switch(swqidatiz) = 1.0
    !
    ! TAG:  +261    Initialization
    !
    ! Sample Input
    !
    !  '+261    SUB-switch: Pcx Opt 4  : PINQID- Molecular Ioniz'    1.0
    !
    ! Input call documentation
    !
    !  'mliz sw'
    !
    switch(swqidmliz) = 1.0
    !
    ! TAG:  +262    Initialization
    !
    ! Sample Input
    !
    !  '+262    SUB-switch: Pcx Opt 4  : PINQID- Recombination  '    1.0
    !
    ! Input call documentation
    !
    !  'rec sw'
    !
    switch(swqidrec) = 1.0
    !
    ! TAG:  +263    Initialization
    !
    ! Sample Input
    !
    !  '+263    SUB-switch: Pcx Opt 4  : PINQID- Charge Exchange'    2.0
    !
    ! Input call documentation
    !
    !  'cx sw'
    !
    switch(swqidcx) = 2.0
    !
    ! TAG:  +264    Initialization
    !
    ! Sample Input
    !
    !  '+264    Switch: PP ElecLoss    : 0.0-off 1.0-XPT 2.0-DIS'    1.0
    !
    ! Input call documentation
    !
    !  'pp elec sw'
    !
    switch(swppelec) = 1.0
    !
    ! TAG:  +265    Initialization
    !
    ! Sample Input
    !
    !  '+265    Switch: PP IonLoss     : 0.0-off 1.0-XPT 2.0-DIS'    1.0
    !
    ! Input call documentation
    !
    !  'pp ion sw'
    !
    switch(swppion) = 1.0
    !
    ! TAG:  +266    Initialization
    !
    ! Sample Input
    !
    !  '+266    Switch: Visc 1 - N calc: 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !  'visc1 sw'
    !
    switch(swvisc1) = 0.0
    !
    ! TAG:  +267    Initialization
    !
    ! Sample Input
    !
    !  '+267    Switch: Momentum loss  : 0.0-off 1.0-on         '    1.0
    !
    ! Input call documentation
    !
    !   'N mom sw'
    !
    switch(swnmom) = 1.0
    !
    ! TAG:  +268    Initialization
    !
    ! Sample Input
    !
    !  '+268    Switch: Iterative Mach : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !   'mach sw'
    !
    switch(swmach) = 0.0
    !
    ! TAG:  +269    Initialization
    !
    ! Sample Input
    !
    !  '+269    Switch: Edge 2D Data   : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !   'e2d sw'
    !
    switch(swe2d) = 0.0
    !
    ! TAG:  +270    Initialization
    !
    ! Sample Input
    !
    !  '+270    Switch: Power Distrib. : 0.0-con 1.0-lin 2.0-xpt'    5.0
    !
    ! Input call documentation
    !
    !   'power sw'
    !
    switch(swpow) = 5.0
    !
    ! TAG:  +271    Initialization
    !
    ! Comments from input
    !
    ! move to after input file read
    !
    ! Sample Input
    !
    !  '+271    Switch: PPlasma PowDist: 0.0-con 1.0-lin 2.0-xpt'    5.0
    !  if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)
    !
    ! Input call documentation
    !
    !   'pp pow'
    !
    switch(swpowp) = 5.0
    !
    ! TAG:  +272    Initialization
    !
    ! Sample Input
    !
    !  '+272    Switch: Gamma Perp     : 0.0-off 1.0-on         '    1.0
    !
    ! Input call documentation
    !
    !  'GamPerp'
    !
    switch(swgperp) = 1.0
    !
    ! TAG:  +273    Initialization
    !
    ! Comments from input
    !
    !     Extra Gperp source/sink term
    ! move to after input file read
    !
    ! Sample Input
    !
    !  '+273    Switch: PP Gamma Perp  : 0.0-off 1.0-on         '    1.0
    !  if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)
    !
    ! Input call documentation
    !
    !  'GamPerpP'
    !
    switch(swgperpp) = 1.0
    !
    ! TAG:  +274    Initialization
    !
    ! Sample Input
    !
    !  '+274    Switch: GPero Src/Sink : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !  'GP Src/Sink'
    !
    switch(swextra) = 0.0
    !
    ! TAG:  +275    Initialization
    !
    ! Sample Input
    !
    !  '+275    Switch: Major Radius   : 0.0-off 1.0-nor 2.0-inv'    0.0
    !
    ! Input call documentation
    !
    !  'MajorRad'
    !
    switch(swmajr) = 0.0
    !
    ! TAG:  +276    Initialization
    !
    ! Sample Input
    !
    !  '+276    Switch: Core Gamma Src : 0.0-off 1.0-all 2.0-xpt'    0.0
    !
    ! Input call documentation
    !
    !  'Core Src'
    !
    switch(swcore) = 0.0
    !
    ! TAG:  +277    Initialization
    !
    ! Sample Input
    !
    !  '+277    Switch: Recomb. Src    : 0.0-off 1.0-PIN        '    1.0
    !
    ! Input call documentation
    !
    !  'Recomb.'
    !
    switch(swrecom) = 1.0
    !
    ! TAG:  +278    Initialization
    !
    ! Sample Input
    !
    !  '+278    Switch: Smooth Centres : 0.0-off 1.0-on         '    0.0
    !
    ! Input call documentation
    !
    !  'Smooth'
    !
    switch(swsmooth) = 0.0
    !
    ! TAG:  +279    Initialization
    !
    ! Sample Input
    !
    !  '+279    Switch: Detached Option: 0.0-off 1.0-out 2.0-in '    0.0
    !
    ! Input call documentation
    !
    !  'Detach'
    !
    switch(swdetach) = 0.0
    !
    ! TAG:  +280    Initialization
    !
    ! Sample Input
    !
    !  '+280    Switch: Error corrected: 0.0-off 1.0-cond       '   10.0
    !
    ! Input call documentation
    !
    !  'ERROR'
    !
    switch(swerror) = 10.0
    !
    ! TAG:  +281    Initialization
    !
    ! Sample Input
    !
    !  '+281 ' 'Automatic DEFAULT Solver condition switches     '
    !  '    DEFAULT applied automatically to these rings'         0
    !
    ! Input call documentation
    !
    !  'DEFAULT SOLVER DATA'
    !
    ndef = 0

    !
    !      TAG 282: SOL22
    !
    !     Initialization of Array input for tag 282 specifying ffric 
    !     values on a ring by ring basis for both targets. 
    !
    n_extffric = 0
    !      extffric = 0.0
    !     
    !     TAG 283: SOL22 - private plasma pressure loss option
    !
    !     Set the default for this value to OFF = 0
    !
    switch(swppress) = 0.0
    !     
    !     TAG 284: SOL22 - debug SOL22 
    !
    !     Set the default for this value to OFF = 0
    !
    debug_sol22 = 0
    !     
    !     TAG 285: SOL22 - debug SOL22 
    !
    !     Set the default for this value to OFF = 0
    !
    debug_sol22_ir = 1
    !     
    !     TAG 286: SOL22 - debug SOL22 
    !
    !     Set the default for this value to OFF = 0
    !
    debug_sol22_ikopt = 1
    !
    !     TAG 287: SOL22 - base ionization source length for algorithmic ionization options
    !
    alg_ion_src_len = 2.0
    !
    !     TAG 288: SOL22 - ring by ring specification of radiation loss parameters
    !
    !     Initialization of Array input for tag 288 specifying radiation parameters
    !     values on a ring by ring basis for both targets. 
    !
    n_extradsrc = 0
    extradsrc = 0.0 
    !
    ! TAG 289: swqperpe - electron perpendicular power flux option
    !          compensates for power terms giving conservation of energy on flux tube
    !          default = OFF
    !
    switch(swqperpe) = 0.0

    !
    ! TAG 290: swqperpi - electron perpendicular power flux option
    !          compensates for power terms giving conservation of energy on flux tube
    !          default = OFF
    !
    switch(swqperpi) = 0.0

    !
    ! TAG 291: DBLSRC_OPT - select option for double source 0=exp+rect  1=rect+rect
    !
    dblsrc_opt = 0

    !
    ! TAG 292: DBLSRC_FRAC - fraction in first source - rest is 1-fraction
    !
    dblsrc_frac = 0.5

    !
    ! TAG 293: DBLSRC 1 - Opt 0 = lams len   Opt 1 = len1 len2
    !
    dblsrc1_p1 = 0.05
    dblsrc1_p2 = 0.1
    !
    ! TAG 294: DBLSRC 2   len1 len2
    !
    dblsrc2_p1 = 0.0
    dblsrc2_p2 = 0.5

    !
    ! TAG 295: switch(swepow) - external electron power term
    !
    switch(swepow) = 0.0

    !
    ! TAG 296: switch(swipow) - external ion power term
    !
    switch(swipow) = 0.0

    !
    ! TAG 297: switch(swepow)=2.0 - external electron power term file name
    !
    ext_epow_fn = 'ext_epow_data.txt'

    !
    ! TAG 298: switch(swipow)=2.0 - external ion power term file name
    !
    ext_ipow_fn = 'ext_ipow_data.txt'

    !
    ! TAG 299: Sol22 print option - used to select more detailed sol22 outputs in the .lim file (fort.6)
    !    sol22_print =0 off
    !
    sol22_print = 0

  end subroutine sol22_initialize_unstructured_input


  subroutine sol22_unstructured_input(tag,line,ierr)
    use mod_io
    use mod_solparams
    use mod_solswitch
    use mod_solcommon
    implicit none
    character*(*) :: tag,line
    integer :: ierr

    !
    !  Inputs for TAG series 2
    !
    if (tag(1:3) .eq. '201') then
       !   Sample input 
       !   '+201    Force Te=Ti through SOL 22  0=off  1=on         '     0
       call divrd(forcet,.TRUE.,0,.true.,1,         'force te=ti'    ,ierr)
    elseif (tag(1:3) .eq. '202') then
       !   Sample input 
       !   '+202    Imposed mach number at the target               '    1.0
       call divrd(initm0,.TRUE.,0.0d0,.FALSE.,0.0d0,'target mach num',IERR)
    elseif (tag(1:3) .eq. '203') then
       !   Sample input 
       !   '+203    Delta mach number for initial iterative solution'    0.1
       !     ------------------------------------------------------------------
       !     Ionization Source
       call divrd(deltam0,.TRUE.,0.0d0,.FALSE.,0.0d0,'delta mach num',IERR)
    elseif (tag(1:3) .eq. '204') then
       !   Sample input 
       !   '+204    Maximum resolution in calculation of m0         ' 0.00001
       call divrd(m0res,.TRUE.,0.0d0,.FALSE.,0.0d0,'Resolution in m0',IERR)
    elseif (tag(1:3) .eq. '205') then
       !   Sample input 
       !   '+205    Ionization Source Lengths  0=Absolute 1=Relative'     1
       call divrd(lensind,.TRUE.,0,.TRUE.,1,     'ion source abs/rel',IERR)
    elseif (tag(1:3) .eq. '206') then
       !   Sample input 
       !   '+206    Start of Ionization Source (for supported opts) '    0.0
       call divrd(lensst,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src start',  IERR)
    elseif (tag(1:3) .eq. '207') then
       !   Sample input 
       !   '+207    End or Length of Ionization Source              '    0.3
       call divrd(lensfi,.TRUE.,0.0d0,.FALSE.,0.0d0,'ion src finish', IERR)
    elseif (tag(1:3) .eq. '208') then
       !   Sample input 
       !   '+208    Decay length of ionization source               '    0.03
       call divrd(lams,.TRUE. ,0.0d0,.FALSE.,0.0d0,'ion decay len   ',IERR)
    elseif (tag(1:3) .eq. '209') then
       !   Sample input 
       !   '+209    Length of radiation source                      '    5.0
       call divrd(lenri,.TRUE.,0.0d0,.FALSE.,0.0d0,'rad source len  ',IERR)
    elseif (tag(1:3) .eq. '210') then
       !   Sample input 
       !   '+210    Decay length of radiation  source               '    0.5
       call divrd(lamri,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad decay len  ',IERR)
    elseif (tag(1:3) .eq. '211') then
       !   Sample input 
       !   '+211    Source strength fraction (frr)                  '    1.0
       call divrd(frri ,.TRUE. ,0.0d0,.FALSE.,0.0d0,'rad power mult ',IERR)
    elseif (tag(1:3) .eq. '212') then
       !   Sample input 
       !   '+212    Garching Model: Alpha = ni/ne ratio             '    1.0
       call divrd(alfimp,.TRUE.,0.0d0,.FALSE.,0.0d0,'nimp/ne ratio  ',IERR)
    elseif (tag(1:3) .eq. '213') then
       !   Sample input 
       !   '+213    Garching Model: Tbase = Tratio denominator      '   15.0
       call divrd(talimp,.TRUE. ,0.0d0,.FALSE.,0.0d0,'base Temp     ',IERR)
    elseif (tag(1:3) .eq. '214') then
       !   Sample input 
       !   '+214    Garching Model: Exponent 1                      '    1.5
       call divrd(ex1imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 1  ',IERR)
    elseif (tag(1:3) .eq. '215') then
       !   Sample input 
       !   '+215    Garching Model: Exponent 2                      '   -3.0
       call divrd(ex2imp ,.FALSE. ,0.0d0,.FALSE.,0.0d0,'exponent 2  ',IERR)
    elseif (tag(1:3) .eq. '216') then
       !   Sample input 
       !   '+216    Gamma correction factor in gammai               '    0.0
       call divrd(gamcor,.false.,0.0d0,.FALSE.,0.0d0,'i power corr. ',IERR)
    elseif (tag(1:3) .eq. '217') then
       !   Sample input 
       !   '+217    Gamma correction factor in gammae               '    0.0
       call divrd(gamecor,.false.,0.0d0,.FALSE.,0.0d0,'e power corr.',IERR)
    elseif (tag(1:3) .eq. '218') then
       !   Sample input 
       !   '+218    CX power coefficeint CEICF                      '    1.0
       call divrd(ceicf,.TRUE. ,0.0d0,.FALSE.,0.0d0, 'CX power frac ',IERR)
    elseif (tag(1:3) .eq. '219') then
       !   Sample input 
       !   '+219    Recycling  source fraction                      '    1.0
       call divrd(recfrac,.TRUE.,0.0d0,.TRUE.,1.0d0,'Recycle frac ',IERR)
    elseif (tag(1:3) .eq. '220') then
       !   Sample input 
       !   '+220    Pei Power Transfer Correction Factor            '    1.0
       call divrd(peicf,.true.  ,0.0d0,.false.,0.0d0,'Pei Correction',ierr)
    elseif (tag(1:3) .eq. '221') then
       !   Sample input 
       !   '+221    Velocity Error Switch       0=Cs   1=const      '     1
       call divrd(velsw,.true.  ,0,.true.,3,       'Vel Error Switch',ierr)
    elseif (tag(1:3) .eq. '222') then
       !   Sample input 
       !   '+222    Distributed Power Start position * SMAX         '    0.10
       !     Gperp Distribution function
       call divrd(spowbeg,.true.,0.0d0,.true.,0.5d0,'Power Dist Beg',ierr)
    elseif (tag(1:3) .eq. '223') then
       !   Sample input 
       !   '+223    Distributed Power End   position * SMAX         '    0.50
       call divrd(spowlen,.true.,0.0d0,.true.,0.5d0,'Power Dist Len',ierr)
    elseif (tag(1:3) .eq. '224') then
       !   Sample input 
       !   '+224    Distributed GPERP particle Fraction- non-uniform'    0.8
       call divrd(gperpfrac,.true.,0.0d0,.true.,1.0d0,'Part Dist Frac',ierr)
    elseif (tag(1:3) .eq. '225') then
       !   Sample input 
       !   '+225    Distributed GPERP Start position * SMAX         '    0.0
       call divrd(gperpbegf,.true.,0.0d0,.true.,0.5d0,'Part Dist Beg',ierr)
    elseif (tag(1:3) .eq. '226') then
       !   Sample input 
       !   '+226    Distributed GPERP End   position * SMAX         '    0.1
       call divrd(gperpendf,.true.,0.0d0,.true.,0.5d0,'Part Dist Len',ierr)
    elseif (tag(1:3) .eq. '227') then
       !   Sample input 
       !   '+227    Gextra Source strength - Target flux multiplier '    0.1
       call divrd(gextra_mult,.true.,0.0d0,.false.,0.0d0,'Gextra flux mult',ierr)
    elseif (tag(1:3) .eq. '228') then
       !   Sample input 
       !   '+228    Gextra Source Start/Stop * SMAX                 ' 0.2   0.35
       call divrd(gextra_src_start,gextra_src_stop,.true.,0.0d0,.true.,1.0d0, 'Gextra SRC start/stop',ierr)
    elseif (tag(1:3) .eq. '229') then
       !   Sample input 
       !   '+229    Gextra Sink   Start/Stop * SMAX                 ' 0.65  0.8
       call divrd(gextra_sink_start,gextra_sink_stop,.true.,0.0d0,.true.,1.0d0,'Gextra SINK start/stop',ierr)
    elseif (tag(1:3) .eq. '230') then
       !   Sample input 
       !   '+230    PP target power loss redistribution range *SMAX '    0.1
       call divrd(pp_pow_dist,.true.,0.0d0,.false.,0.0d0,'PP Pow Dist',ierr)
    elseif (tag(1:3) .eq. '231') then
       !   Sample input 
       !   '+231    Knot Start Index for E2D Option 9               '     8
       call divrd(ike2d,.true.,1,.false.,0,'IK start Index for E2D-9',ierr)
    elseif (tag(1:3) .eq. '232') then
       !   Sample input 
       !   '+232    Plasma Fill option for missing knots - E2D Opt 9'     1
       call divrd(fillopt,.true.,0,.true.,3,'Gap fill option - E2D-9',ierr)
    elseif (tag(1:3) .eq. '233') then
       !   Sample input 
       !   '+233    Qe Term - Temperature Cutoff (eV)               '    0.0
       call divrd(tcutqe,.true.,0.0d0,.false.,0.0d0,'Cut T-PINQE',ierr)
    elseif (tag(1:3) .eq. '234') then
       !   Sample input 
       !   '+234    PINQID - Atomic Ionization    - T cutoff (eV)   '    0.0
       call divrd(tcutatiz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDATIZ',ierr)
    elseif (tag(1:3) .eq. '235') then
       !   Sample input 
       !   '+235    PINQID - Molecular Ionization - T cutoff (eV)   '    0.0
       call divrd(tcutmliz,.true.,0.0d0,.false.,0.0d0,'Cut T-QIDMLIZ',ierr)
    elseif (tag(1:3) .eq. '236') then
       !   Sample input 
       !   '+236    PINQID - Recombination        - T cutoff (eV)   '    0.0
       call divrd(tcutrec,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDREC',ierr)
    elseif (tag(1:3) .eq. '237') then
       !   Sample input 
       !   '+237    Qi Term/PINQID-Charge Exchange- T cutoff (eV)   '    0.0
       call divrd(tcutcx,.true.,0.0d0,.false.,0.0d0,'Cut T- QIDCX',ierr)
    elseif (tag(1:3) .eq. '238') then
       !   Sample input 
       !   '+238    PINQID - CX option 1 - Reference T - (eV)       '    1.0
       call divrd(trefcx,.true.,0.0d0,.false.,0.0d0,'REF T- QIDCX 1',ierr)
    elseif (tag(1:3) .eq. '239') then
       !   Sample input 
       !   '+239    Minimum Temperature allowed in Solver (spec<0)  '    0.1
       call divrd(tmin,.false.,0.0d0,.false.,0.0d0,'Min. Allowed T',ierr)
    elseif (tag(1:3) .eq. '240') then
       !   Sample input 
       !   '+240    Minimum T allowed as fraction of Tmax reached   '    0.5
       call divrd(dropfrac,.true.,0.0d0,.true.,1.0d0,'Allowed T-drop',ierr)
    elseif (tag(1:3) .eq. '241') then
       !   Sample input 
       !   '+241    Momentum loss term multiplier   (Usually 1.0)   '    1.0
       call divrd(smom_mult,.false.,0.0d0,.false.,0.0d0,'Mom.Loss Multiplier',ierr)
    elseif (tag(1:3) .eq. '242') then
       !   Sample input 
       !   '+242    Friction factor for Momentum loss formula       '    0.2
       call divrd(ffric,.false.,0.0d0,.false.,0.0d0,'Mom.loss frac  ',ierr)
    elseif (tag(1:3) .eq. '243') then
       !   Sample input 
       !   '+243    Length of the Momentum loss region * Smax       '    0.1
       call divrd(lenmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom source len',IERR)
    elseif (tag(1:3) .eq. '244') then
       !   Sample input 
       !   '+244    Decay length of the Momentum loss region * Smax '    0.02
       call divrd(lammom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'mom decay len ',IERR)
    elseif (tag(1:3) .eq. '245') then
       !   Sample input 
       !   '+245    Ratio of CX to IZ events (fixed)                '    1.0
       call divrd(rcxmom,.TRUE. ,0.0d0,.FALSE.,0.0d0,'cx/iz ratio   ',IERR)
    elseif (tag(1:3) .eq. '246') then
       !   Sample input 
       !   '+246    Te cutoff for increased CX multiplier (eV)      '    5.0
       call divrd(tcxmom,.TRUE.,1.0001d0,.FALSE.,0.0d0,'T for CXmult',ierr)
    elseif (tag(1:3) .eq. '247') then
       !   Sample input 
       !   '+247    Te lower limit cutoff for CX multiplier (eV)    '    1.0
       call divrd(tcxcut,.true.,0.0d0,.false.,0.0d0,'Cut T -  CXmult',ierr)
    elseif (tag(1:3) .eq. '248') then
       !   Sample input 
       !   '+248    PINQE multiplier                                '    1.0
       call divrd(qesrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE mult',ierr)
    elseif (tag(1:3) .eq. '249') then
       !   Sample input 
       !   '+249    PRAD option 3 multiplier (x PINQE)              '    0.5
       call divrd(radsrc_mult,.false.,0.0d0,.FALSE.,0.0d0,'PINQE based PRAD mult',ierr)
    elseif (tag(1:3) .eq. '250') then
       !   Sample input 
       !   '+250    Initial number of stages for Runge Kutta steps  '    100
       call divrd(ndiv,.TRUE., 1,.FALSE., 0,'NUMBER OF STEPS    ',IERR)
    elseif (tag(1:3) .eq. '251') then
       !   Sample input 
       !   '+251    Switch: Ionization Opt : 0.0-exp 1.0+ - others  '    2.0
       call divrd(switch(swion),.true. ,0.0,.false.,0.0, 'alt ion',ierr)
    elseif (tag(1:3) .eq. '252') then
       !   Sample input 
       !   '+252    Switch: Initial IonOpt : 0.0-exp 3.0+ - others  '    0.0
       call divrd(switch(swioni),.true. ,0.0,.false.,0.0,'init ion',ierr)
    elseif (tag(1:3) .eq. '253') then
       !   Sample input 
       !   '+253    Switch: PPlasma IonOpt : 0.0-exp 3.0+ - others  '   -1.0
       call divrd(switch(swionp),.false.,0.0,.false.,0.0,'pp ion',ierr)
    elseif (tag(1:3) .eq. '254') then
       !   Sample input 
       !   '+254    Switch: 5/2 nv * kT    : 0.0-off 1.0-on         '    1.0
       call divrd(switch(swcond),.true.,0.0,.false.,0.0, 'cond sw',ierr)
    elseif (tag(1:3) .eq. '255') then
       !   Sample input 
       !   '+255    Switch: 1/2 m v^3 * n  : 0.0-off 1.0-on         '    1.0
       call divrd(switch(swconv),.true.,0.0,.false.,0.0, 'conv sw',ierr)
    elseif (tag(1:3) .eq. '256') then
       !   Sample input 
       !   '+256    Switch: Prad           : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swprad),.true.,0.0,.false.,0.0, 'prad sw',ierr)
    elseif (tag(1:3) .eq. '257') then
       !   Sample input 
       !   '+257    Switch: Phelpi         : 0.0-off 1.0-on         '    2.0
       call divrd(switch(swphelp),.true.,0.0,.false.,0.0,'phelp sw',ierr)
    elseif (tag(1:3) .eq. '258') then
       !   Sample input 
       !   '+258    Switch: Pei            : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swpei),.true.,0.0,.false.,0.0,  'pei sw ',ierr)
    elseif (tag(1:3) .eq. '259') then
       !   Sample input 
       !   '+259    Switch: Pcx            : 0.0-off 1.0-on         '    5.0
       call divrd(switch(swpcx),.true.,0.0,.false.,0.0,  'pcx sw ',ierr)
    elseif (tag(1:3) .eq. '260') then
       !   Sample input 
       !   '+260    SUB-switch: Pcx Opt 4  : PINQID- Atomic Ioniz.  '    1.0
       call divrd(switch(swqidatiz),.true.,0.0,.false.,0.0,'atiz sw',ierr)
    elseif (tag(1:3) .eq. '261') then
       !   Sample input 
       !   '+261    SUB-switch: Pcx Opt 4  : PINQID- Molecular Ioniz'    1.0
       call divrd(switch(swqidmliz),.true.,0.0,.false.,0.0,'mliz sw',ierr)
    elseif (tag(1:3) .eq. '262') then
       !   Sample input 
       !   '+262    SUB-switch: Pcx Opt 4  : PINQID- Recombination  '    1.0
       call divrd(switch(swqidrec),.true.,0.0,.false.,0.0,'rec sw',ierr)
    elseif (tag(1:3) .eq. '263') then
       !   Sample input 
       !   '+263    SUB-switch: Pcx Opt 4  : PINQID- Charge Exchange'    2.0
       call divrd(switch(swqidcx),.true.,0.0,.false.,0.0,'cx sw',ierr)
    elseif (tag(1:3) .eq. '264') then
       !   Sample input 
       !   '+264    Switch: PP ElecLoss    : 0.0-off 1.0-XPT 2.0-DIS'    1.0
       call divrd(switch(swppelec),.true.,0.0,.false.,0.0,'pp elec sw',ierr)
    elseif (tag(1:3) .eq. '265') then
       !   Sample input 
       !   '+265    Switch: PP IonLoss     : 0.0-off 1.0-XPT 2.0-DIS'    1.0
       call divrd(switch(swppion),.true.,0.0,.false.,0.0,'pp ion sw',ierr)
    elseif (tag(1:3) .eq. '266') then
       !   Sample input 
       !   '+266    Switch: Visc 1 - N calc: 0.0-off 1.0-on         '    0.0
       call divrd(switch(swvisc1),.true.,0.0,.true.,0.0,'visc1 sw',ierr)
    elseif (tag(1:3) .eq. '267') then
       !   Sample input 
       !   '+267    Switch: Momentum loss  : 0.0-off 1.0-on         '    1.0
       call divrd(switch(swnmom),.true.,0.0,.false.,0.0, 'N mom sw',ierr)
    elseif (tag(1:3) .eq. '268') then
       !   Sample input 
       !   '+268    Switch: Iterative Mach : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swmach),.true.,0.0,.false.,0.0, 'mach sw',ierr)
    elseif (tag(1:3) .eq. '269') then
       !   Sample input 
       !   '+269    Switch: Edge 2D Data   : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swe2d),.false.,0.0,.false.,0.0, 'e2d sw',ierr)
    elseif (tag(1:3) .eq. '270') then
       !   Sample input 
       !   '+270    Switch: Power Distrib. : 0.0-con 1.0-lin 2.0-xpt'    5.0
       call divrd(switch(swpow),.true.,0.0,.false.,0.0, 'power sw',ierr)
    elseif (tag(1:3) .eq. '271') then
       !   Sample input 
       !   '+271    Switch: PPlasma PowDist: 0.0-con 1.0-lin 2.0-xpt'    5.0
       !   if (switch(swpowp).eq.-1.0) switch(swpowp) = switch(swpow)
       ! move to after input file read
       call divrd(switch(swpowp),.false.,0.0,.false.,0.0, 'pp pow',ierr)
    elseif (tag(1:3) .eq. '272') then
       !   Sample input 
       !   '+272    Switch: Gamma Perp     : 0.0-off 1.0-on         '    1.0
       call divrd(switch(swgperp),.true.,0.0,.false.,0.0,'GamPerp',ierr)
    elseif (tag(1:3) .eq. '273') then
       !   Sample input 
       !   '+273    Switch: PP Gamma Perp  : 0.0-off 1.0-on         '    1.0
       !   if (switch(swgperpp).eq.-1.0) switch(swgperpp) = switch(swgperp)
       !     Extra Gperp source/sink term
       ! move to after input file read
       call divrd(switch(swgperpp),.false.,0.0,.false.,0.0,'GamPerpP',ierr)
    elseif (tag(1:3) .eq. '274') then
       !   Sample input 
       !   '+274    Switch: GPero Src/Sink : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swextra),.true.,0.0,.false.,0.0,'GP Src/Sink',ierr)
    elseif (tag(1:3) .eq. '275') then
       !   Sample input 
       !   '+275    Switch: Major Radius   : 0.0-off 1.0-nor 2.0-inv'    0.0
       call divrd(switch(swmajr),.true.,0.0,.false.,0.0,'MajorRad',ierr)
    elseif (tag(1:3) .eq. '276') then
       !   Sample input 
       !   '+276    Switch: Core Gamma Src : 0.0-off 1.0-all 2.0-xpt'    0.0
       call divrd(switch(swcore),.true.,0.0,.false.,0.0,'Core Src',ierr)
    elseif (tag(1:3) .eq. '277') then
       !   Sample input 
       !   '+277    Switch: Recomb. Src    : 0.0-off 1.0-PIN        '    1.0
       call divrd(switch(swrecom),.true.,0.0,.false.,0.0,'Recomb.',ierr)
    elseif (tag(1:3) .eq. '278') then
       !   Sample input 
       !   '+278    Switch: Smooth Centres : 0.0-off 1.0-on         '    0.0
       call divrd(switch(swsmooth),.true.,0.0,.false.,0.0,'Smooth',ierr)
    elseif (tag(1:3) .eq. '279') then
       !   Sample input 
       !   '+279    Switch: Detached Option: 0.0-off 1.0-out 2.0-in '    0.0
       call divrd(switch(swdetach),.true.,0.0,.false.,0.0,'Detach',ierr)
    elseif (tag(1:3) .eq. '280') then
       !   Sample input 
       !   '+280    Switch: Error corrected: 0.0-off 1.0-cond       '   10.0
       call divrd(switch(swerror),.true.,0.0,.false.,0.0,'ERROR',ierr)
    elseif (tag(1:3) .eq. '281') then
       !   Sample input 
       !   '+281 ' 'Automatic DEFAULT Solver condition switches     '
       !   '    DEFAULT applied automatically to these rings'         0
       call divrd(deflist,ndef,mxspts,0.0,SOL22_HI,.FALSE.,0.0,SOL22_MACHHI,2,'DEFAULT SOLVER DATA',IERR)



       !
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     Options added by David for various parts of the code. 
       !
       ! -----------------------------------------------------------------------
       !
       !
       ! -----------------------------------------------------------------------
       !
       !     Trying to keep options in alphabetical order. 
       !
    elseif (tag(1:3).eq.'282') then  
       !     TAG 282 - SOL option 22 - reads in FFRIC values for 
       !               the momentum loss option so that they may vary
       !               from ring to ring and target to target.
       !             - the format is IR FFRIC1 FFRIC2
       !               jdemod - format is IR FFRIC1 LEN1 FFRIC2 LEN2 for option 1
       !                      - format is IR FFRIC1 LAM1 FFRIC2 LAM2 for option 2 - LEN is set to SMAX/2
       !
       !               FFRIC1 etc applies to the first half ring - this would
       !               be the OUTER half ring for X-point up grids and the 
       !               INNER half for X-point down grids. This may also
       !               be designated target 2. 
       !               FFRIC2 etc applies to the second half of the ring  
       !     
       !     Note: the tag line precedes a standard DIVIMP array input of
       !           three lines.  
       !
       CALL RDQARN(extffric,n_extffric,MXSPTS,-sol22_MACHHI,sol22_MACHHI,.FALSE.,&
            -sol22_machhi,sol22_MACHHI,4,'SET OF MOM-LOSS COEF BY RING',IERR)

    elseif (tag(1:3).eq.'283') then  
       !
       !     jdemod
       !     TAG 283 - SOL option 22 - reads the value for switch(swppress) 
       !             - this is the pfz target pressure redistribution option. 
       !     

       CALL ReadR(line,switch(swppress),0.0,2.0,'SOL22 PFZ PRESSURE LOSS OPT')
       !
       !        write(0,*) '*283 - Read switch',swppress,switch(swppress)
       !
    elseif (tag(1:3).eq.'284') then  
       !     jdemod
       !     TAG 284 - SOL option 22 - SOL22 debugging switch
       !     

       CALL ReadI(line,debug_sol22,0,1,'SOL22: DEBUG SWITCH')
       !
    elseif (tag(1:3).eq.'285') then  
       !     jdemod
       !     TAG 285 - SOL option 22 - SOL22 debug IR for detailed profile
       !     

       ! jdemod
       ! maxnrs might not be defined depending on the usage so this is just an index.
       ! Make imax very large for now
       !CALL ReadI(line,debug_sol22_ir,1,maxnrs,'SOL22 DEBUG RING')
       CALL ReadI(line,debug_sol22_ir,1,100000,'SOL22 DEBUG RING')
       !
    elseif (tag(1:3).eq.'286') then  
       !
       !     jdemod
       !     TAG 286 - SOL option 22 - SOL22 debug IKOPT for detailed profile
       !     

       CALL ReadI(line,debug_sol22_ikopt,1,2,'SOL22 DEBUG IKOPT-RING END')
    elseif (tag(1:3).eq.'287') then  
       !
       !     TAG 287: SOL22 - base ionization source length for algorithmic ionization options
       !

       CALL ReadR(line,alg_ion_src_len,0.0,sol22_MACHHI,'DEFAULT IONIZATION SOURCE LENGTH FOR ALGORITHMIC OPTIONS')
    elseif (tag(1:3).eq.'288') then  
       !
       !     TAG 288 - SOL option 22 - reads in radiation values for 
       !               the radiation option so that they may vary
       !               from ring to ring and target to target.
       !             - the format is IR LENR1 LAMR1 FFR1 LENR2 LAMR2 FFR2
       !
       !               "1" applies to the first half ring - this would
       !               be the OUTER half ring for X-point up grids and the 
       !               INNER half for X-point down grids. This may also
       !               be designated target 2. 
       !               "2" applies to the second half of the ring  
       !     
       !     Note: the tag line precedes a standard DIVIMP array input of
       !           three lines.  
       !
       CALL RDQARN(extradsrc,n_extradsrc,MXSPTS,-sol22_MACHHI,sol22_MACHHI,.FALSE.,&
            -sol22_machhi,sol22_MACHHI,6,'SET OF RADIATION COEF BY RING',IERR)
       !
    elseif (tag(1:3).eq.'289') then  
       !
       !     jdemod
       !     TAG 289 - SOL option 22 - reads the value for switch(swqperpe) 
       !             - switch to balance electron power on flux tube
       !             - default is OFF
       !     

       CALL ReadR(line,switch(swqperpe),0.0,1.0,'SOL22: ELECTRON PERPENDICULAR POWER FLUX OPTION')
       !
    elseif (tag(1:3).eq.'290') then  
       !
       !     jdemod
       !     TAG 290 - SOL option 22 - reads the value for switch(swqperpi) 
       !             - switch to balance ion power on flux tube
       !             - default is OFF
       !     

       CALL ReadR(line,switch(swqperpi),0.0,1.0,'SOL22: ION PERPENDICULAR POWER FLUX OPTION')
       !
    elseif (tag(1:3).eq.'291') then  
       !
       !     jdemod  tag 291-294
       !     Double shape ionization source
       !     Option 0 : exp + rect (default)
       !     Option 1 : rect + rect
       !
       !     Other paramters
       !     dblsrc_frac : fraction of particles in first source
       !     dblsrc1_p1, dlbsrc1_p2 : paramters for first source - either lams len or len1 len2
       !     dblsrc2_p1, dlbsrc2_p2 : paramters for second source - len1 len2
       !     
       CALL ReadI(line,dblsrc_opt,0,1,'OPTION for dual profile ionization source')

       !
    elseif (tag(1:3).eq.'292') then  
       !
       !     jdemod
       !     TAG 292
       !     dblsrc_frac : fraction of particles in first source
       !     
       CALL ReadR(line,dblsrc_frac,0.0,1.0,'SOL22: DBLSRC_FRAC : fraction of particles in first source')
       !
    elseif (tag(1:3).eq.'293') then  
       !
       !     jdemod
       !     TAG 293
       !     dblsrc1_p1, dlbsrc1_p2 : paramters for first source - either lams len or len1 len2
       !     
       CALL Read2R(line,dblsrc1_p1,dblsrc1_p2,0.0,1.0,'SOL22: DBLSRC1_P1,2 : parameters for first source')
       !
    elseif (tag(1:3).eq.'294') then  
       !
       !     jdemod
       !     TAG 294
       !     dblsrc2_p1, dlbsrc2_p2 : paramters for second source len1 len2
       !     
       CALL Read2R(line,dblsrc2_p1,dblsrc2_p2,0.0,1.0,'SOL22: DBLSRC2_P1,2 : parameters for second source')
       !
    elseif (tag(1:3).eq.'295') then  
       !
       !     jdemod
       !     TAG 295 - SOL option 22 - reads the value for switch(swepow) for external electron power term 
       !             - default is OFF = 0.0
       !             - option 1 = read data from div aux input file - Looks for tag EXTEPOW:
       !                          Add dataset to the divimp aux input file using the following format
       ! 
       !                          EXTEPOW:
       !                          format(6e18.10)   ((data_array(ik,ir,1),ik=1,nks(ir)),ir=1,nrs)
       ! 
       !             - option 2 = read raw data from R,Z,EPOW file and interpolate onto DIVIMP mesh
       !                          (external data source file must be connected to ext_epow_data.txt)
       !                          External file format is:
       !                          # or $ = Comments ignored on file read
       !                          NROWS: <nr>           
       !                          NCOLS: <nz> 
       !                          DATA:
       !                           R Z RAD       These may be space or comma separated
       !               Any grid coordinates outside the supplied data will be assigned a zero           
       !               Data must be entered from low R to high R and low Z to high Z
       !               i.e. Start from lower left corner of data set
       !               Code could be added to support additional data formats or organization 
       !     

       CALL ReadR(line,switch(swepow),0.0,2.0,'SOL22: EXTERNAL ELECTRON POWER TERM OPTION')
       !
    elseif (tag(1:3).eq.'296') then  
       !
       !     jdemod
       !     TAG 296 - SOL option 22 - reads the value for switch(swipow) for external ion power term 
       !             - default is OFF = 0.0
       !             - option 1 = read data from div aux input file - Looks for tag EXTIPOW:
       ! 
       !                          EXTIPOW:
       !                          format(6e18.10)   ((data_array(ik,ir,1),ik=1,nks(ir)),ir=1,nrs)
       ! 
       !             - option 2 = read raw data from R,Z,IPOW file and interpolate onto DIVIMP mesh
       !                          (external data source file must be connected to the ext_ipow_data.txt)
       !                          External file format is:
       !                          # or $ = Comments ignored on file read
       !                          NROWS: <nr>           
       !                          NCOLS: <nz> 
       !                          DATA:
       !                           R Z RAD       These may be space or comma separated
       !               Any grid coordinates outside the supplied data will be assigned a zero           
       !               Data must be entered from low R to high R and low Z to high Z
       !               i.e. Start from lower left corner of data set
       !               Code could be added to support additional data formats or organization 
       !     

       CALL ReadR(line,switch(swipow),0.0,2.0,'SOL22: EXTERNAL ION POWER TERM OPTION')
       !
    elseif (tag(1:3).eq.'297') then  
       !
       !     jdemod
       !     TAG 297 - SOL option 22 - file name with external electron power source for swepow option = 2.0
       !                               default = ext_epow_data.txt 
       !     

       CALL ReadC(line,ext_epow_fn,'SOL22: EXTERNAL ELECTRON POWER TERM FILE NAME')

       !
    elseif (tag(1:3).eq.'298') then  
       !
       !     jdemod
       !     TAG 298 - SOL option 22 - file name with external ion power source for swipow option = 2.0
       !                               default = ext_ipow_data.txt 
       !     

       CALL ReadC(line,ext_ipow_fn,'SOL22: EXTERNAL ION POWER TERM FILE NAME')

       !

    elseif (tag(1:3).eq.'299') then  

       ! jdemod
       ! TAG 299: Sol22 print option - used to select more detailed sol22 outputs in the .lim file (fort.6)
       !    sol22_print =0 off
       !
       CALL ReadI(line,sol22_print,0,2,'Selects level of detail in SOL22 output to the .lim(fort.6) file ')

    endif

  end subroutine sol22_unstructured_input

end module mod_sol22_input
