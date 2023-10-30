module mod_sol23_input
  use debug_options
  implicit none

  !
  !     common block containing all of the sol23 specifiable input parameters
  !
  !     for explanation contact w. fundamenski
  !
  !     -*-fortran-*-
  ! common /sol23_input/sol23_izlen,sol23_izlam,sol23_izoffset,sol23_momlen,sol23_sm_s,&
  !      sol23_sm_r, sol23_relax, sol23_maxtol,sol23_rmstol, sol23_maxpow,sol23_intopt,&
  !      sol23_bndcond, sol23_seed, sol23_perp,sol23_adaptgrid,sol23_par_ptipte,sol23_par_adaptnr ,&
  !      sol23_par_debugflag, sol23_par_debugnr,sol23_par_refresh , sol23_par_artvisc2 ,&
  !     sol23_par_artvisc4,sol23_par_dtf     , sol23_par_dtg      , sol23_par_grid0  ,&
  !     sol23_par_gridexp , sol23_par_itermax  , sol23_par_ga1    ,sol23_par_ga2     ,&
  !      sol23_par_ga3      , sol23_par_ga4    ,sol23_par_updtqpit, sol23_par_updtqp2  ,&
  !      sol23_par_updtqp3,sol23_par_updtqpte, sol23_par_garelax  , sol23_par_gaiter ,&
  !     sol23_par_limitte , sol23_par_celldte  ,sol23_par_updtbcte,sol23_par_updtdel0,&
  !      sol23_par_updtdel1 ,sol23_par_updtdelm,sol23_par_qbrelax , sol23_par_gridg    ,&
  !     sol23_par_grid_dx0,sol23_par_gae     , sol23_par_gai      , sol23_par_tectrl ,&
  !     sol23_par_drflag  , sol23_par_dsflag   , sol23_par_limrel1,sol23_par_limrel2 ,&
  !      sol23_par_limrel3  , sol23_par_limrel4,sol23_par_tulimit , sol23_par_g0relax  ,&
  !      sol23_par_p0relax,sol23_par_dpdrflag, sol23_par_dpdrtemin,sol23_par_dpdrstep,&
  !     sol23_par_nulimit , sol23_par_nuflag   , sol23_par_vnmult ,sol23_par_pinmom  , sol23_par_emolec   ,&
  !     sol23_par_rec_heat,sol23_par_pinqimult,sol23_par_pinqiflag, sol23_par_prring0,&
  !     sol23_par_prring1 , sol23_par_qperp34  , sol23_par_qeiflag,sol23_par_chie    ,&
  !      sol23_par_joule    , sol23_par_fluxexp,sol23_par_qrec    , sol23_par_fzrad    ,&
  !      sol23_par_dvmrel
  !
  ! save /sol23_input/
  !
  integer,public :: sol23_intopt,sol23_bndcond,sol23_seed,sol23_perp,sol23_adaptgrid,&
       sol23_par_adaptnr,sol23_par_debugflag,sol23_par_debugnr,sol23_par_refresh,sol23_par_itermax,&
       sol23_par_updtqpit,sol23_par_gaiter,sol23_par_limitte,sol23_par_tectrl,&
       sol23_par_drflag,sol23_par_dsflag,sol23_par_dpdrflag,sol23_par_nuflag,sol23_par_pinmom,&
       sol23_par_pinqiflag,sol23_par_prring0,sol23_par_prring1,sol23_par_qperp34,&
       sol23_par_qeiflag,sol23_par_chie,sol23_par_joule,sol23_par_fluxexp,sol23_par_qrec,&
       sol23_par_dvmrel

  integer,public :: readin_sol23_params

  real,public :: sol23_izlen,sol23_izlam,sol23_izoffset,sol23_momlen,sol23_sm_s,sol23_sm_r,&
       sol23_relax,sol23_maxtol,sol23_rmstol,sol23_maxpow,sol23_par_ptipte,sol23_par_artvisc2,&
       sol23_par_artvisc4,sol23_par_dtf,sol23_par_dtg,sol23_par_grid0,sol23_par_gridexp,&
       sol23_par_ga1,sol23_par_ga2,sol23_par_ga3,sol23_par_ga4,sol23_par_updtqp2,&
       sol23_par_updtqp3,sol23_par_updtqpte,sol23_par_garelax,sol23_par_celldte,&
       sol23_par_updtbcte,sol23_par_updtdel0,sol23_par_updtdel1,sol23_par_updtdelm,sol23_par_qbrelax,&
       sol23_par_gridg,sol23_par_grid_dx0,sol23_par_gae,sol23_par_gai,&
       sol23_par_limrel1,sol23_par_limrel2,sol23_par_limrel3,sol23_par_limrel4,sol23_par_tulimit,&
       sol23_par_g0relax,sol23_par_p0relax,sol23_par_dpdrtemin,sol23_par_dpdrstep,&
       sol23_par_nulimit,sol23_par_vnmult,sol23_par_emolec,sol23_par_rec_heat,sol23_par_pinqimult,&
       sol23_par_fzrad

  public :: allocate_mod_sol23_input,deallocate_mod_sol23_input,sol23_initialize_unstructured_input,&
       sol23_read_unstructured_input

contains


  subroutine sol23_initialize_unstructured_input
    use mod_params
    use mod_io
    implicit none
    !
    !  Initialization for TAG series 3
    !
    !
    ! TAG:  +300    Initialization
    !
    ! Sample Input
    !
    !  '+300 SOL23  Parameter read option  0-NO 1-YES            '     0
    !
    ! Input call documentation
    !
    !  'READ SOL23 PARMAS OPT'
    !
    readin_sol23_params = 0
    !
    ! TAG:  +301    Initialization
    !
    ! Sample Input
    !
    !  '+301 SOL23  ptipte ................................      '    1.0
    !
    ! Input call documentation
    !
    !  'par_ptipte       '
    !
    sol23_par_ptipte = 1.0
    !
    ! TAG:  +302    Initialization
    !
    ! Sample Input
    !
    !  '+302 SOL23  adaptnr                                      '     2
    !
    ! Input call documentation
    !
    !  'par_adaptnr      '
    !
    sol23_par_adaptnr = 2
    !
    ! TAG:  +303    Initialization
    !
    ! Sample Input
    !
    !  '+303 SOL23  debugflag                                    '     0
    !
    ! Input call documentation
    !
    !  'par_debugflag    '
    !
    sol23_par_debugflag = 0
    !
    ! TAG:  +304    Initialization
    !
    ! Sample Input
    !
    !  '+304 SOL23  debugnr                                      '     0
    !
    ! Input call documentation
    !
    !  'par_debugnr      '
    !
    sol23_par_debugnr = 0
    !
    ! TAG:  +305    Initialization
    !
    ! Sample Input
    !
    !  '+305 SOL23  refresh                                      '    100
    !
    ! Input call documentation
    !
    !  'par_refresh      '
    !
    sol23_par_refresh = 100
    !
    ! TAG:  +306    Initialization
    !
    ! Sample Input
    !
    !  '+306 SOL23  artvisc2 ..............................      '    0.3
    !
    ! Input call documentation
    !
    !  'par_artvisc2     '
    !
    sol23_par_artvisc2 = 0.3
    !
    ! TAG:  +307    Initialization
    !
    ! Sample Input
    !
    !  '+307 SOL23  artvisc4                                     '    0.01
    !
    ! Input call documentation
    !
    !  'par_artvisc4     '
    !
    sol23_par_artvisc4 = 0.01
    !
    ! TAG:  +308    Initialization
    !
    ! Sample Input
    !
    !  '+308 SOL23  dtf                                          '    0.5
    !
    ! Input call documentation
    !
    !  'par_dtf          '
    !
    sol23_par_dtf = 0.5
    !
    ! TAG:  +309    Initialization
    !
    ! Sample Input
    !
    !  '+309 SOL23  dtg                                          '   1.0E24
    !
    ! Input call documentation
    !
    !  'par_dtg          '
    !
    sol23_par_dtg = 1.0E24
    !
    ! TAG:  +310    Initialization
    !
    ! Sample Input
    !
    !  '+310 SOL23  grid0                                        '    0.1
    !
    ! Input call documentation
    !
    !  'par_grid0        '
    !
    sol23_par_grid0 = 0.1
    !
    ! TAG:  +311    Initialization
    !
    ! Sample Input
    !
    !  '+311 SOL23  gridexp ................................     '    1.25
    !
    ! Input call documentation
    !
    !  'par_gridexp      '
    !
    sol23_par_gridexp = 1.25
    !
    ! TAG:  +312    Initialization
    !
    ! Sample Input
    !
    !  '+312 SOL23  itermax                                      '    1000
    !
    ! Input call documentation
    !
    !  'par_itermax      '
    !
    sol23_par_itermax = 1000
    !
    ! TAG:  +313    Initialization
    !
    ! Sample Input
    !
    !  '+313 SOL23  ga1                                          '    10.0
    !
    ! Input call documentation
    !
    !  'par_ga1          '
    !
    sol23_par_ga1       = 10.0
    !
    ! TAG:  +314    Initialization
    !
    ! Sample Input
    !
    !  '+314 SOL23  ga2                                          '    0.01
    !
    ! Input call documentation
    !
    !  'par_ga2          '
    !
    sol23_par_ga2     = 0.01
    !
    ! TAG:  +315    Initialization
    !
    ! Sample Input
    !
    !  '+315 SOL23  ga3                                          '    0.2
    !
    ! Input call documentation
    !
    !  'par_ga3          '
    !
    sol23_par_ga3     = 0.2
    !
    ! TAG:  +316    Initialization
    !
    ! Sample Input
    !
    !  '+316 SOL23  ga4   .................................      '    10.0
    !
    ! Input call documentation
    !
    !  'par_ga4          '
    !
    sol23_par_ga4      = 10.0
    !
    ! TAG:  +317    Initialization
    !
    ! Sample Input
    !
    !  '+317 SOL23  updtqpit                                     '     0
    !
    ! Input call documentation
    !
    !  'par_updtqpit     '
    !
    sol23_par_updtqpit = 0
    !
    ! TAG:  +318    Initialization
    !
    ! Sample Input
    !
    !  '+318 SOL23  updtqp2                                      '    0.003
    !
    ! Input call documentation
    !
    !  'par_updtqp2      '
    !
    sol23_par_updtqp2  = 0.003
    !
    ! TAG:  +319    Initialization
    !
    ! Sample Input
    !
    !  '+319 SOL23  updtqp3                                      '    0.03
    !
    ! Input call documentation
    !
    !  'par_updtqp3      '
    !
    sol23_par_updtqp3 = 0.03
    !
    ! TAG:  +320    Initialization
    !
    ! Sample Input
    !
    !  '+320 SOL23  updtqpte                                     '    0.01
    !
    ! Input call documentation
    !
    !  'par_updtqpte     '
    !
    sol23_par_updtqpte = 0.01
    !
    ! TAG:  +321    Initialization
    !
    ! Sample Input
    !
    !  '+321 SOL23  garelax ................................     '    0.03
    !
    ! Input call documentation
    !
    !  'par_garelax      '
    !
    sol23_par_garelax = 0.03
    !
    ! TAG:  +322    Initialization
    !
    ! Sample Input
    !
    !  '+322 SOL23  gaiter                                       '     1
    !
    ! Input call documentation
    !
    !  'par_gaiter       '
    !
    sol23_par_gaiter  = 1
    !
    ! TAG:  +323    Initialization
    !
    ! Sample Input
    !
    !  '+323 SOL23  limitte                                      '     1
    !
    ! Input call documentation
    !
    !  'par_limitte      '
    !
    sol23_par_limitte   = 1
    !
    ! TAG:  +324    Initialization
    !
    ! Sample Input
    !
    !  '+324 SOL23  celldte                                      '    0.95
    !
    ! Input call documentation
    !
    !  'par_celldte      '
    !
    sol23_par_celldte = 0.95
    !
    ! TAG:  +325    Initialization
    !
    ! Sample Input
    !
    !  '+325 SOL23  updtbcte                                     '    30.0
    !
    ! Input call documentation
    !
    !  'par_updtbcte     '
    !
    sol23_par_updtbcte = 30.0
    !
    ! TAG:  +326    Initialization
    !
    ! Sample Input
    !
    !  '+326 SOL23  updtdel0 ..............................      '    1.0
    !
    ! Input call documentation
    !
    !  'par_updtdel0     '
    !
    sol23_par_updtdel0 = 1.0
    !
    ! TAG:  +327    Initialization
    !
    ! Sample Input
    !
    !  '+327 SOL23  updtdel1                                     '    1.0
    !
    ! Input call documentation
    !
    !  'par_updtdel1     '
    !
    sol23_par_updtdel1 = 1.0
    !
    ! TAG:  +328    Initialization
    !
    ! Sample Input
    !
    !  '+328 SOL23  updtdelm                                     '    0.1
    !
    ! Input call documentation
    !
    !  'par_updtdelm     '
    !
    sol23_par_updtdelm = 0.1
    !
    ! TAG:  +329    Initialization
    !
    ! Sample Input
    !
    !  '+329 SOL23  qbrelax                                      '    0.02
    !
    ! Input call documentation
    !
    !  'par_qbrelax      '
    !
    sol23_par_qbrelax = 0.02
    !
    ! TAG:  +330    Initialization
    !
    ! Sample Input
    !
    !  '+330 SOL23  gridg                                        '    0.03
    !
    ! Input call documentation
    !
    !  'par_gridg        '
    !
    sol23_par_gridg    = 0.03
    !
    ! TAG:  +331    Initialization
    !
    ! Sample Input
    !
    !  '+331 SOL23  grid_dx0 ...............................     '    0.001
    !
    ! Input call documentation
    !
    !  'par_grid_dx0     '
    !
    sol23_par_grid_dx0 = 0.001
    !
    ! TAG:  +332    Initialization
    !
    ! Sample Input
    !
    !  '+332 SOL23  gae                                          '    5.0
    !
    ! Input call documentation
    !
    !  'par_gae          '
    !
    sol23_par_gae     = 5.0
    !
    ! TAG:  +333    Initialization
    !
    ! Sample Input
    !
    !  '+333 SOL23  gai                                          '    2.5
    !
    ! Input call documentation
    !
    !  'par_gai          '
    !
    sol23_par_gai       = 2.5
    !
    ! TAG:  +334    Initialization
    !
    ! Sample Input
    !
    !  '+334 SOL23  tectrl                                       '     1
    !
    ! Input call documentation
    !
    !  'par_tectrl       '
    !
    sol23_par_tectrl  = 1
    !
    ! TAG:  +335    Initialization
    !
    ! Sample Input
    !
    !  '+335 SOL23  drflag                                       '     0
    !
    ! Input call documentation
    !
    !  'par_drflag       '
    !
    sol23_par_drflag  = 0
    !
    ! TAG:  +336    Initialization
    !
    ! Sample Input
    !
    !  '+336 SOL23  dsflag  ................................     '     0
    !
    ! Input call documentation
    !
    !  'par_dsflag       '
    !
    sol23_par_dsflag   = 0
    !
    ! TAG:  +337    Initialization
    !
    ! Sample Input
    !
    !  '+337 SOL23  limrel1                                      '    0.03
    !
    ! Input call documentation
    !
    !  'par_limrel1      '
    !
    sol23_par_limrel1  = 0.03
    !
    ! TAG:  +338    Initialization
    !
    ! Sample Input
    !
    !  '+338 SOL23  limrel2                                      '    0.1
    !
    ! Input call documentation
    !
    !  'par_limrel2      '
    !
    sol23_par_limrel2  = 0.1
    !
    ! TAG:  +339    Initialization
    !
    ! Sample Input
    !
    !  '+339 SOL23  limrel3                                      '    0.03
    !
    ! Input call documentation
    !
    !  'par_limrel3      '
    !
    sol23_par_limrel3  = 0.03
    !
    ! TAG:  +340    Initialization
    !
    ! Sample Input
    !
    !  '+340 SOL23  limrel4                                      '    0.03
    !
    ! Input call documentation
    !
    !  'par_limrel4      '
    !
    sol23_par_limrel4  = 0.03
    !
    ! TAG:  +341    Initialization
    !
    ! Sample Input
    !
    !  '+341 SOL23  tulimit  ...............................     '    5.0
    !
    ! Input call documentation
    !
    !  'par_tulimit      '
    !
    sol23_par_tulimit = 5.0
    !
    ! TAG:  +342    Initialization
    !
    ! Sample Input
    !
    !  '+342 SOL23  g0relax                                      '    0.1
    !
    ! Input call documentation
    !
    !  'par_g0relax      '
    !
    sol23_par_g0relax = 0.1
    !
    ! TAG:  +343    Initialization
    !
    ! Sample Input
    !
    !  '+343 SOL23  p0relax                                      '    0.02
    !
    ! Input call documentation
    !
    !  'par_p0relax      '
    !
    sol23_par_p0relax   = 0.02
    !
    ! TAG:  +344    Initialization
    !
    ! Sample Input
    !
    !  '+344 SOL23  dpdrflag                                     '     0
    !
    ! Input call documentation
    !
    !  'par_dpdrflag     '
    !
    sol23_par_dpdrflag = 0
    !
    ! TAG:  +345    Initialization
    !
    ! Sample Input
    !
    !  '+345 SOL23  dpdrtemin                                    '    1.0
    !
    ! Input call documentation
    !
    !  'par_dpdrtemin    '
    !
    sol23_par_dpdrtemin = 1.0
    !
    ! TAG:  +346    Initialization
    !
    ! Sample Input
    !
    !  '+346 SOL23  dpdrstep  ..............................     '    0.1
    !
    ! Input call documentation
    !
    !  'par_dpdrstep     '
    !
    sol23_par_dpdrstep = 0.1
    !
    ! TAG:  +347    Initialization
    !
    ! Sample Input
    !
    !  '+347 SOL23  nuflag                                       '     0
    !
    ! Input call documentation
    !
    !  'par_nuflag       '
    !
    sol23_par_nuflag   = 0
    !
    ! TAG:  +348    Initialization
    !
    ! Sample Input
    !
    !  '+348 SOL23  nulimit                                      '   1.0E20
    !
    ! Input call documentation
    !
    !  'par_nulimit      '
    !
    sol23_par_nulimit = 1.0E20
    !
    ! TAG:  +349    Initialization
    !
    ! Sample Input
    !
    !  '+349 SOL23  vnmult                                       '    1.0
    !
    ! Input call documentation
    !
    !  'par_vnmult       '
    !
    sol23_par_vnmult  = 1.0
    !
    ! TAG:  +350    Initialization
    !
    ! Sample Input
    !
    !  '+350 SOL23  pinmom                                       '     0
    !
    ! Input call documentation
    !
    !  'par_pinmom       '
    !
    sol23_par_pinmom   = 0
    !
    ! TAG:  +351    Initialization
    !
    ! Sample Input
    !
    !  '+351 SOL23  emolec                                       '    3.0
    !
    ! Input call documentation
    !
    !  'par_emolec       '
    !
    sol23_par_emolec  = 3.0
    !
    ! TAG:  +352    Initialization
    !
    ! Sample Input
    !
    !  '+352 SOL23  rec_heat  ...............................    '    13.6
    !
    ! Input call documentation
    !
    !  'par_rec_heat     '
    !
    sol23_par_rec_heat = 13.6
    !
    ! TAG:  +353    Initialization
    !
    ! Sample Input
    !
    !  '+353 SOL23  pinqimult                                    '    1.0
    !
    ! Input call documentation
    !
    !  'par_pinqimult    '
    !
    sol23_par_pinqimult = 1.0
    !
    ! TAG:  +354    Initialization
    !
    ! Sample Input
    !
    !  '+354 SOL23  pinqiflag                                    '     0
    !
    ! Input call documentation
    !
    !  'par_pinqiflag    '
    !
    sol23_par_pinqiflag = 0
    !
    ! TAG:  +355    Initialization
    !
    ! Sample Input
    !
    !  '+355 SOL23  prring0                                      '     7
    !
    ! Input call documentation
    !
    !  'par_prring0      '
    !
    sol23_par_prring0 = 7
    !
    ! TAG:  +356    Initialization
    !
    ! Sample Input
    !
    !  '+356 SOL23  prring1                                      '     3
    !
    ! Input call documentation
    !
    !  'par_prring1      '
    !
    sol23_par_prring1  = 3
    !
    ! TAG:  +357    Initialization
    !
    ! Sample Input
    !
    !  '+357 SOL23  qperp34 .................................    '     0
    !
    ! Input call documentation
    !
    !  'par_qperp34      '
    !
    sol23_par_qperp34  = 0
    !
    ! TAG:  +358    Initialization
    !
    ! Sample Input
    !
    !  '+358 SOL23  qeiflag                                      '     1
    !
    ! Input call documentation
    !
    !  'par_qeiflag      '
    !
    sol23_par_qeiflag  = 1
    !
    ! TAG:  +359    Initialization
    !
    ! Sample Input
    !
    !  '+359 SOL23  chie                                         '     1
    !
    ! Input call documentation
    !
    !  'par_chie         '
    !
    sol23_par_chie    = 1
    !
    ! TAG:  +360    Initialization
    !
    ! Sample Input
    !
    !  '+360 SOL23  joule                                        '     1
    !
    ! Input call documentation
    !
    !  'par_joule        '
    !
    sol23_par_joule    = 1
    !
    ! TAG:  +361    Initialization
    !
    ! Sample Input
    !
    !  '+361 SOL23  fluxexp                                      '     1
    !
    ! Input call documentation
    !
    !  'par_fluxexp      '
    !
    sol23_par_fluxexp = 1
    !
    ! TAG:  +362    Initialization
    !
    ! Sample Input
    !
    !  '+362 SOL23  qrec    ..................................   '     1
    !
    ! Input call documentation
    !
    !  'par_qrec         '
    !
    sol23_par_qrec    = 1
    !
    ! TAG:  +363    Initialization
    !
    ! Sample Input
    !
    !  '+363 SOL23  fzrad                                        '    1.0
    !
    ! Input call documentation
    !
    !  'par_fzrad        '
    !
    sol23_par_fzrad = 1.0
    !
    ! TAG:  +364    Initialization
    !
    ! Sample Input
    !
    !  '+364 SOL23  dvmrel                                       '     2
    !
    ! Input call documentation
    !
    !  'par_dvmrel       '
    !
    sol23_par_dvmrel = 2
    !
    ! TAG:  +365    Initialization
    !
    ! Sample Input
    !
    !  '+365 SOL23  intopt                                       '     1
    !
    ! Input call documentation
    !
    !  'intopt'
    !
    sol23_intopt = 1
    !
    ! TAG:  +366    Initialization
    !
    ! Sample Input
    !
    !  '+366 SOL23  adaptgrid                                    '     0
    !
    ! Input call documentation
    !
    !  'adaptgrid'
    !
    sol23_adaptgrid = 0
    !
    ! TAG:  +367    Initialization
    !
    ! Sample Input
    !
    !  '+367 SOL23  seed                                         '     1
    !
    ! Input call documentation
    !
    !  'seed'
    !
    sol23_seed = 1
    !
    ! TAG:  +368    Initialization
    !
    ! Sample Input
    !
    !  '+368 SOL23  izlen   ..................................   '    0.5
    !
    ! Input call documentation
    !
    !  'izlen'
    !
    sol23_izlen = 0.5
    !
    ! TAG:  +369    Initialization
    !
    ! Sample Input
    !
    !  '+369 SOL23  izlam                                        '    0.03
    !
    ! Input call documentation
    !
    !  'izlam'
    !
    sol23_izlam = 0.03
    !
    ! TAG:  +370    Initialization
    !
    ! Sample Input
    !
    !  '+370 SOL23  izoffset                                     '    0.03
    !
    ! Input call documentation
    !
    !  'izoffset'
    !
    sol23_izoffset = 0.03
    !
    ! TAG:  +371    Initialization
    !
    ! Sample Input
    !
    !  '+371 SOL23  momlen                                       '    0.1
    !
    ! Input call documentation
    !
    !  'momlen'
    !
    sol23_momlen = 0.1
    !
    ! TAG:  +372    Initialization
    !
    ! Sample Input
    !
    !  '+372 SOL23  relax                                        '    0.01
    !
    ! Input call documentation
    !
    !  'relax'
    !
    sol23_relax = 0.01
    !
    ! TAG:  +373    Initialization
    !
    ! Sample Input
    !
    !  '+373 SOL23  maxtol                                       '    3.0E-3
    !
    ! Input call documentation
    !
    !  'maxtol'
    !
    sol23_maxtol = 3.0E-3
    !
    ! TAG:  +374    Initialization
    !
    ! Sample Input
    !
    !  '+374 SOL23  rmstol  .................................    '    3.0E-4
    !
    ! Input call documentation
    !
    !  'rmstol'
    !
    sol23_rmstol = 3.0E-4
    !
    ! TAG:  +375    Initialization
    !
    ! Sample Input
    !
    !  '+375 SOL23  Qperp                                        '    1
    !
    ! Input call documentation
    !
    !  'perp'
    !
    sol23_perp = 1


  end subroutine sol23_initialize_unstructured_input



  subroutine sol23_read_unstructured_input(tag,line,ierr)
    use mod_params
    use mod_io
    implicit none
    character*(*) :: tag,line
    integer :: ierr

    !
    !  Inputs for TAG series 3
    !
    if (tag(1:3) .eq. '300') then
       !   Sample input 
       !   '+300 SOL23  Parameter read option  0-NO 1-YES            '     0
       call divrd(readin_sol23_params,.TRUE.,0 ,.true.,1,'READ SOL23 PARMAS OPT',IERR)
    elseif (tag(1:3) .eq. '301') then
       !   Sample input 
       !   '+301 SOL23  ptipte ................................      '    1.0
       call divrd(sol23_par_ptipte,.true.,0.0,.true.,1.0E8,'par_ptipte       ',ierr)
    elseif (tag(1:3) .eq. '302') then
       !   Sample input 
       !   '+302 SOL23  adaptnr                                      '     2
       call divrd(sol23_par_adaptnr,.true.,0  ,.true.,10  ,'par_adaptnr      ',ierr)
    elseif (tag(1:3) .eq. '303') then
       !   Sample input 
       !   '+303 SOL23  debugflag                                    '     0
       call divrd(sol23_par_debugflag,.true.,0 ,.true.,10 ,'par_debugflag    ',ierr)
    elseif (tag(1:3) .eq. '304') then
       !   Sample input 
       !   '+304 SOL23  debugnr                                      '     0
       call divrd(sol23_par_debugnr,.true.,0  ,.true.,10  ,'par_debugnr      ',ierr)
    elseif (tag(1:3) .eq. '305') then
       !   Sample input 
       !   '+305 SOL23  refresh                                      '    100
       call divrd(sol23_par_refresh,.true.,0  ,.true.,10000,'par_refresh      ',ierr)
    elseif (tag(1:3) .eq. '306') then
       !   Sample input 
       !   '+306 SOL23  artvisc2 ..............................      '    0.3
       call divrd(sol23_par_artvisc2,.true.,0.0,.true.,1.0E8,'par_artvisc2     ',ierr)
    elseif (tag(1:3) .eq. '307') then
       !   Sample input 
       !   '+307 SOL23  artvisc4                                     '    0.01
       call divrd(sol23_par_artvisc4,.true.,0.0,.true.,1.0E8,'par_artvisc4     ',ierr)
    elseif (tag(1:3) .eq. '308') then
       !   Sample input 
       !   '+308 SOL23  dtf                                          '    0.5
       call divrd(sol23_par_dtf,.true.,0.0,.true.,1.0E8,'par_dtf          ',ierr)
    elseif (tag(1:3) .eq. '309') then
       !   Sample input 
       !   '+309 SOL23  dtg                                          '   1.0E24
       call divrd(sol23_par_dtg,.true.,0.0,.true.,1.0E36,'par_dtg          ',ierr)
    elseif (tag(1:3) .eq. '310') then
       !   Sample input 
       !   '+310 SOL23  grid0                                        '    0.1
       call divrd(sol23_par_grid0,.true.,0.0,.true.,1.0E8,'par_grid0        ',ierr)
    elseif (tag(1:3) .eq. '311') then
       !   Sample input 
       !   '+311 SOL23  gridexp ................................     '    1.25
       call divrd(sol23_par_gridexp,.true.,0.0,.true.,1.0E8,'par_gridexp      ',ierr)
    elseif (tag(1:3) .eq. '312') then
       !   Sample input 
       !   '+312 SOL23  itermax                                      '    1000
       call divrd(sol23_par_itermax,.true.,0  ,.true.,10000,'par_itermax      ',ierr)
    elseif (tag(1:3) .eq. '313') then
       !   Sample input 
       !   '+313 SOL23  ga1                                          '    10.0
       call divrd(sol23_par_ga1      ,.true.,0.0,.true.,1.0E8,'par_ga1          ',ierr)
    elseif (tag(1:3) .eq. '314') then
       !   Sample input 
       !   '+314 SOL23  ga2                                          '    0.01
       call divrd(sol23_par_ga2    ,.true.,0.0,.true.,1.0E8,'par_ga2          ',ierr)
    elseif (tag(1:3) .eq. '315') then
       !   Sample input 
       !   '+315 SOL23  ga3                                          '    0.2
       call divrd(sol23_par_ga3    ,.true.,0.0,.true.,1.0E8,'par_ga3          ',ierr)
    elseif (tag(1:3) .eq. '316') then
       !   Sample input 
       !   '+316 SOL23  ga4   .................................      '    10.0
       call divrd(sol23_par_ga4     ,.true.,0.0,.true.,1.0E8,'par_ga4          ',ierr)
    elseif (tag(1:3) .eq. '317') then
       !   Sample input 
       !   '+317 SOL23  updtqpit                                     '     0
       call divrd(sol23_par_updtqpit,.true.,0  ,.true.,10   ,'par_updtqpit     ',ierr)
    elseif (tag(1:3) .eq. '318') then
       !   Sample input 
       !   '+318 SOL23  updtqp2                                      '    0.003
       call divrd(sol23_par_updtqp2 ,.true.,0.0,.true.,1.0E8,'par_updtqp2      ',ierr)
    elseif (tag(1:3) .eq. '319') then
       !   Sample input 
       !   '+319 SOL23  updtqp3                                      '    0.03
       call divrd(sol23_par_updtqp3,.true.,0.0,.true.,1.0E8,'par_updtqp3      ',ierr)
    elseif (tag(1:3) .eq. '320') then
       !   Sample input 
       !   '+320 SOL23  updtqpte                                     '    0.01
       call divrd(sol23_par_updtqpte,.true.,0.0,.true.,1.0E8,'par_updtqpte     ',ierr)
    elseif (tag(1:3) .eq. '321') then
       !   Sample input 
       !   '+321 SOL23  garelax ................................     '    0.03
       call divrd(sol23_par_garelax,.true.,0.0,.true.,1.0E8,'par_garelax      ',ierr)
    elseif (tag(1:3) .eq. '322') then
       !   Sample input 
       !   '+322 SOL23  gaiter                                       '     1
       call divrd(sol23_par_gaiter ,.true.,0  ,.true.,10000,'par_gaiter       ',ierr)
    elseif (tag(1:3) .eq. '323') then
       !   Sample input 
       !   '+323 SOL23  limitte                                      '     1
       call divrd(sol23_par_limitte  ,.true.,0  ,.true.,10   ,'par_limitte      ',ierr)
    elseif (tag(1:3) .eq. '324') then
       !   Sample input 
       !   '+324 SOL23  celldte                                      '    0.95
       call divrd(sol23_par_celldte,.true.,0.0,.true.,1.0E8,'par_celldte      ',ierr)
    elseif (tag(1:3) .eq. '325') then
       !   Sample input 
       !   '+325 SOL23  updtbcte                                     '    30.0
       call divrd(sol23_par_updtbcte,.true.,0.0,.true.,1.0E8,'par_updtbcte     ',ierr)
    elseif (tag(1:3) .eq. '326') then
       !   Sample input 
       !   '+326 SOL23  updtdel0 ..............................      '    1.0
       call divrd(sol23_par_updtdel0,.true.,0.0,.true.,1.0E8,'par_updtdel0     ',ierr)
    elseif (tag(1:3) .eq. '327') then
       !   Sample input 
       !   '+327 SOL23  updtdel1                                     '    1.0
       call divrd(sol23_par_updtdel1,.true.,0.0,.true.,1.0E8,'par_updtdel1     ',ierr)
    elseif (tag(1:3) .eq. '328') then
       !   Sample input 
       !   '+328 SOL23  updtdelm                                     '    0.1
       call divrd(sol23_par_updtdelm,.true.,0.0,.true.,1.0E8,'par_updtdelm     ',ierr)
    elseif (tag(1:3) .eq. '329') then
       !   Sample input 
       !   '+329 SOL23  qbrelax                                      '    0.02
       call divrd(sol23_par_qbrelax,.true.,0.0,.true.,1.0E8,'par_qbrelax      ',ierr)
    elseif (tag(1:3) .eq. '330') then
       !   Sample input 
       !   '+330 SOL23  gridg                                        '    0.03
       call divrd(sol23_par_gridg   ,.true.,0.0,.true.,1.0E8,'par_gridg        ',ierr)
    elseif (tag(1:3) .eq. '331') then
       !   Sample input 
       !   '+331 SOL23  grid_dx0 ...............................     '    0.001
       call divrd(sol23_par_grid_dx0,.true.,0.0,.true.,1.0E8,'par_grid_dx0     ',ierr)
    elseif (tag(1:3) .eq. '332') then
       !   Sample input 
       !   '+332 SOL23  gae                                          '    5.0
       call divrd(sol23_par_gae    ,.true.,0.0,.true.,1.0E8,'par_gae          ',ierr)
    elseif (tag(1:3) .eq. '333') then
       !   Sample input 
       !   '+333 SOL23  gai                                          '    2.5
       call divrd(sol23_par_gai      ,.true.,0.0,.true.,1.0E8,'par_gai          ',ierr)
    elseif (tag(1:3) .eq. '334') then
       !   Sample input 
       !   '+334 SOL23  tectrl                                       '     1
       call divrd(sol23_par_tectrl ,.true.,0  ,.true.,10   ,'par_tectrl       ',ierr)
    elseif (tag(1:3) .eq. '335') then
       !   Sample input 
       !   '+335 SOL23  drflag                                       '     0
       call divrd(sol23_par_drflag ,.true.,0  ,.true.,10   ,'par_drflag       ',ierr)
    elseif (tag(1:3) .eq. '336') then
       !   Sample input 
       !   '+336 SOL23  dsflag  ................................     '     0
       call divrd(sol23_par_dsflag  ,.true.,0  ,.true.,10   ,'par_dsflag       ',ierr)
    elseif (tag(1:3) .eq. '337') then
       !   Sample input 
       !   '+337 SOL23  limrel1                                      '    0.03
       call divrd(sol23_par_limrel1 ,.true.,0.0,.true.,1.0E8,'par_limrel1      ',ierr)
    elseif (tag(1:3) .eq. '338') then
       !   Sample input 
       !   '+338 SOL23  limrel2                                      '    0.1
       call divrd(sol23_par_limrel2 ,.true.,0.0,.true.,1.0E8,'par_limrel2      ',ierr)
    elseif (tag(1:3) .eq. '339') then
       !   Sample input 
       !   '+339 SOL23  limrel3                                      '    0.03
       call divrd(sol23_par_limrel3 ,.true.,0.0,.true.,1.0E8,'par_limrel3      ',ierr)
    elseif (tag(1:3) .eq. '340') then
       !   Sample input 
       !   '+340 SOL23  limrel4                                      '    0.03
       call divrd(sol23_par_limrel4 ,.true.,0.0,.true.,1.0E8,'par_limrel4      ',ierr)
    elseif (tag(1:3) .eq. '341') then
       !   Sample input 
       !   '+341 SOL23  tulimit  ...............................     '    5.0
       call divrd(sol23_par_tulimit,.true.,0.0,.true.,1.0E8,'par_tulimit      ',ierr)
    elseif (tag(1:3) .eq. '342') then
       !   Sample input 
       !   '+342 SOL23  g0relax                                      '    0.1
       call divrd(sol23_par_g0relax,.true.,0.0,.true.,1.0E8,'par_g0relax      ',ierr)
    elseif (tag(1:3) .eq. '343') then
       !   Sample input 
       !   '+343 SOL23  p0relax                                      '    0.02
       call divrd(sol23_par_p0relax  ,.true.,0.0,.true.,1.0E8,'par_p0relax      ',ierr)
    elseif (tag(1:3) .eq. '344') then
       !   Sample input 
       !   '+344 SOL23  dpdrflag                                     '     0
       call divrd(sol23_par_dpdrflag,.true.,0  ,.true.,10   ,'par_dpdrflag     ',ierr)
    elseif (tag(1:3) .eq. '345') then
       !   Sample input 
       !   '+345 SOL23  dpdrtemin                                    '    1.0
       call divrd(sol23_par_dpdrtemin,.true.,0.0,.true.,1.0E8,'par_dpdrtemin    ',ierr)
    elseif (tag(1:3) .eq. '346') then
       !   Sample input 
       !   '+346 SOL23  dpdrstep  ..............................     '    0.1
       call divrd(sol23_par_dpdrstep,.true.,0.0,.true.,1.0E8,'par_dpdrstep     ',ierr)
    elseif (tag(1:3) .eq. '347') then
       !   Sample input 
       !   '+347 SOL23  nuflag                                       '     0
       call divrd(sol23_par_nuflag  ,.true.,0  ,.true.,10   ,'par_nuflag       ',ierr)
    elseif (tag(1:3) .eq. '348') then
       !   Sample input 
       !   '+348 SOL23  nulimit                                      '   1.0E20
       call divrd(sol23_par_nulimit,.true.,0.0,.true.,1.0E24,'par_nulimit      ',ierr)
    elseif (tag(1:3) .eq. '349') then
       !   Sample input 
       !   '+349 SOL23  vnmult                                       '    1.0
       call divrd(sol23_par_vnmult ,.true.,0.0,.true.,1.0E8,'par_vnmult       ',ierr)
    elseif (tag(1:3) .eq. '350') then
       !   Sample input 
       !   '+350 SOL23  pinmom                                       '     0
       call divrd(sol23_par_pinmom  ,.true.,0  ,.true.,10   ,'par_pinmom       ',ierr)
    elseif (tag(1:3) .eq. '351') then
       !   Sample input 
       !   '+351 SOL23  emolec                                       '    3.0
       call divrd(sol23_par_emolec ,.true.,0.0,.true.,1.0E8,'par_emolec       ',ierr)
    elseif (tag(1:3) .eq. '352') then
       !   Sample input 
       !   '+352 SOL23  rec_heat  ...............................    '    13.6
       call divrd(sol23_par_rec_heat,.true.,0.0,.true.,1.0E8,'par_rec_heat     ',ierr)
    elseif (tag(1:3) .eq. '353') then
       !   Sample input 
       !   '+353 SOL23  pinqimult                                    '    1.0
       call divrd(sol23_par_pinqimult,.true.,0.0,.true.,1.0E8,'par_pinqimult    ',ierr)
    elseif (tag(1:3) .eq. '354') then
       !   Sample input 
       !   '+354 SOL23  pinqiflag                                    '     0
       call divrd(sol23_par_pinqiflag,.true.,0  ,.true.,10   ,'par_pinqiflag    ',ierr)
    elseif (tag(1:3) .eq. '355') then
       !   Sample input 
       !   '+355 SOL23  prring0                                      '     7
       call divrd(sol23_par_prring0,.true.,0  ,.true.,10   ,'par_prring0      ',ierr)
    elseif (tag(1:3) .eq. '356') then
       !   Sample input 
       !   '+356 SOL23  prring1                                      '     3
       call divrd(sol23_par_prring1 ,.true.,0  ,.true.,10   ,'par_prring1      ',ierr)
    elseif (tag(1:3) .eq. '357') then
       !   Sample input 
       !   '+357 SOL23  qperp34 .................................    '     0
       call divrd(sol23_par_qperp34 ,.true.,0  ,.true.,10   ,'par_qperp34      ',ierr)
    elseif (tag(1:3) .eq. '358') then
       !   Sample input 
       !   '+358 SOL23  qeiflag                                      '     1
       call divrd(sol23_par_qeiflag ,.true.,0  ,.true.,10   ,'par_qeiflag      ',ierr)
    elseif (tag(1:3) .eq. '359') then
       !   Sample input 
       !   '+359 SOL23  chie                                         '     1
       call divrd(sol23_par_chie   ,.true.,0  ,.true.,10   ,'par_chie         ',ierr)
    elseif (tag(1:3) .eq. '360') then
       !   Sample input 
       !   '+360 SOL23  joule                                        '     1
       call divrd(sol23_par_joule   ,.true.,0  ,.true.,10   ,'par_joule        ',ierr)
    elseif (tag(1:3) .eq. '361') then
       !   Sample input 
       !   '+361 SOL23  fluxexp                                      '     1
       call divrd(sol23_par_fluxexp,.true.,0  ,.true.,10   ,'par_fluxexp      ',ierr)
    elseif (tag(1:3) .eq. '362') then
       !   Sample input 
       !   '+362 SOL23  qrec    ..................................   '     1
       call divrd(sol23_par_qrec   ,.true.,0  ,.true.,10   ,'par_qrec         ',ierr)
    elseif (tag(1:3) .eq. '363') then
       !   Sample input 
       !   '+363 SOL23  fzrad                                        '    1.0
       call divrd(sol23_par_fzrad,.true.,0.0,.true.,1.0E8,'par_fzrad        ',ierr)
    elseif (tag(1:3) .eq. '364') then
       !   Sample input 
       !   '+364 SOL23  dvmrel                                       '     2
       call divrd(sol23_par_dvmrel,.true.,1,.true.,10,'par_dvmrel       ',ierr)
    elseif (tag(1:3) .eq. '365') then
       !   Sample input 
       !   '+365 SOL23  intopt                                       '     1
       call divrd(sol23_intopt,.TRUE.,0 ,.true.,1,'intopt',IERR)
    elseif (tag(1:3) .eq. '366') then
       !   Sample input 
       !   '+366 SOL23  adaptgrid                                    '     0
       call divrd(sol23_adaptgrid,.TRUE.,0 ,.true.,5,'adaptgrid',IERR)
    elseif (tag(1:3) .eq. '367') then
       !   Sample input 
       !   '+367 SOL23  seed                                         '     1
       call divrd(sol23_seed,.TRUE.,1 ,.true.,2,'seed',IERR)
    elseif (tag(1:3) .eq. '368') then
       !   Sample input 
       !   '+368 SOL23  izlen   ..................................   '    0.5
       call divrd(sol23_izlen,.true.,0.0,.false.,0.5,'izlen',ierr)
    elseif (tag(1:3) .eq. '369') then
       !   Sample input 
       !   '+369 SOL23  izlam                                        '    0.03
       call divrd(sol23_izlam,.true.,0.0,.false.,0.5,'izlam',ierr)
    elseif (tag(1:3) .eq. '370') then
       !   Sample input 
       !   '+370 SOL23  izoffset                                     '    0.03
       call divrd(sol23_izoffset,.true.,0.0,.false.,0.5,'izoffset',ierr)
    elseif (tag(1:3) .eq. '371') then
       !   Sample input 
       !   '+371 SOL23  momlen                                       '    0.1
       call divrd(sol23_momlen,.true.,0.0,.false.,0.0,'momlen',ierr)
    elseif (tag(1:3) .eq. '372') then
       !   Sample input 
       !   '+372 SOL23  relax                                        '    0.01
       call divrd(sol23_relax,.true.,0.0,.true.,1.0,'relax',ierr)
    elseif (tag(1:3) .eq. '373') then
       !   Sample input 
       !   '+373 SOL23  maxtol                                       '    3.0E-3
       call divrd(sol23_maxtol,.true.,0.0,.true.,1.0,'maxtol',ierr)
    elseif (tag(1:3) .eq. '374') then
       !   Sample input 
       !   '+374 SOL23  rmstol  .................................    '    3.0E-4
       call divrd(sol23_rmstol,.true.,0.0,.true.,1.0,'rmstol',ierr)
    elseif (tag(1:3) .eq. '375') then
       !   Sample input 
       !   '+375 SOL23  Qperp                                        '    1
       call divrd(sol23_perp,.true.,0 ,.true.,10 ,'perp',ierr)
    endif


  end subroutine sol23_read_unstructured_input


  subroutine allocate_mod_sol23_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_sol23_input','ALLOCATE')


  end subroutine allocate_mod_sol23_input


  subroutine deallocate_mod_sol23_input
    implicit none

    call pr_trace('mod_sol23_input','DEALLOCATE')


  end subroutine deallocate_mod_sol23_input

end module mod_sol23_input
