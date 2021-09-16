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

  public :: allocate_mod_sol23_input,deallocate_mod_sol23_input

contains

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