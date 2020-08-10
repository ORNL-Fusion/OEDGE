module mod_sol23_com
  use debug_options
  implicit none

  !
  !     sol23 - interface variables
  !
  !     the following common block contains the information that is
  !     passed to sol option 23 for calculating a cfd solution
  !     for the background plasma on a ring by ring basis.
  !
  !
  !     variable list and descriptions:
  !     - all arrays are declared to be of length:  maxpts - see params block
  !     - l = logical    (.true. or .false.)
  !     - i = integer
  !     - r = real
  !     - a = array
  !
  !     name          type  purpose
  !
  ! general:
  !
  !     seedplasma    i     0 = use passed sources, 1+=generate seed sources
  !     solve_half    i     0 = solve first half, 1=solve second half
  !                         2 = solve entire ring (not supported yet)
  !     intopt        i     integration option (method of calculating integrals)
  !                         0 = stepwise   1 = linearly interpolated
  !     bndcond_opt   i     boundary condition option
  !                         0 = default
  !     s23_adaptgrid i     switch on and off adaptive grid methods
  !                         0 = off
  !                         1 = on
  !
  ! seed plasma:
  !
  !     s23_izlen     r*8   seed ionization length (m)
  !     s23_izlam     r*8   seed ionization lambda - decay length (m)
  !     s23_izoffset  r*8   seed ionization source offset distance (m)
  !     s23_momlen    r*8   seed momentum length (m)
  !
  ! sources:
  !
  !     part_src      r*8 a particle source density for each cell (m-3)
  !     ionpow_src    r*8 a ion energy source for each cell (w/m-3)
  !     elecpow_src   r*8 a electron energy source for each cell (w/m-3)
  !     mom_src       r*8 a momentum source for each cell (not used at present)
  !
  !     h_neut        r*8 a hydrogen neutral density in each cell (m-3)
  !     h2_neut       r*8 a neutral hydrogen molecule density (m-3)
  !     th_neut       r*8 a temperature of atomic hydrogen (ev)
  !
  ! boundary conditions:
  !
  !     mach_0        r*8   mach number at s=0 target
  !     isat_0        r*8   saturation current at s=0 (assuming m=1) (m-2s-1)
  !     ne_0          r*8   density at s=0 target (m-3)
  !     te_0          r*8   electron temperature at s=0 target (ev)
  !     ti_0          r*8   ion temperature at s=0 target (ev)
  !
  !     mach_smax     r*8   mach number at s=smax target
  !     isat_smax     r*8   saturation current at s=smax (assuming m=1) (m-2s-1)
  !     ne_smax       r*8   density at s=smax target  (m-3)
  !     te_smax       r*8   electron temperature at s=smax target (ev)
  !     ti_smax       r*8   ion temperature at s=smax target (ev)
  !
  !     gae_const     r*8   electron sheath power transmission (default=5.0)
  !     gai_const     r*8   ion sheath power transmission (default=2.5)
  !
  !     gtarg         r*8 a target particle fluxes (1->0,2->smax,3->total)
  !     ionptarg      r*8 a ion power fluxes       (1->0,2->smax,3->total)
  !     elecptarg     r*8 a electron power fluxes  (1->0,2->smax,3->total)
  !
  !     curpow        r*8   limit for total power losses on ring
  !
  ! geometric:
  !
  !     smax          r*8   total length of ring (m)
  !     npts          i     number of points on the ring
  !     spts          r*8 a cell centre point array - in terms of s (m)
  !     sbnds         r*8 a cell boundaries - in terms of s (m)
  !     ring_type     i     0 = sol ring,   1 = pp ring
  !
  ! background plasma: (output and input for subsequent iterations)
  !
  !     mb            r*8   mass of background ions  (amu)
  !     cb            r*8   charge of background ions (z-ion)
  !
  !     ne            r*8 a background density (m-3)
  !     te            r*8 a electron tempertaure (ev)
  !     ti            r*8 a background ion temperature (ev)
  !     vb            r*8 a background flow velocity (m/s)
  !
  !     -*-fortran-*-
  ! common /sol23/ divimp_call,seedplasma,solve_half,bndcond_opt,s23_izlen,s23_izlam,&
  !     s23_izoffset,s23_momlen,s23_sm_s, s23_sm_r, s23_relax, s23_maxtol,s23_rmstol, s23_maxpow_0,&
  !      s23_maxpow_n,s23_powlim, s23_momlim, s23_zhi0,s23_masslim, s23_maxmom,&
  !      s23_maxmass,s23_par_ptipte,s23_par_adaptnr , s23_par_debugflag, s23_par_debugnr,&
  !     s23_par_refresh , s23_par_artvisc2 , s23_par_artvisc4,s23_par_dtf     , s23_par_dtg      ,&
  !      s23_par_grid0  ,s23_par_gridexp , s23_par_itermax  , s23_par_updtqpit,&
  !     s23_par_ga1    ,s23_par_ga2     , s23_par_ga3      , s23_par_ga4    ,s23_par_updtqp2  ,&
  !      s23_par_updtqp3,s23_par_updtqpte, s23_par_garelax  , s23_par_gaiter ,&
  !     s23_par_limitte , s23_par_celldte  , s23_par_updtbcte,s23_par_updtdel0, s23_par_updtdel1 ,&
  !      s23_par_updtdelm,s23_par_qbrelax , s23_par_gridg    , s23_par_grid_dx0,&
  !     s23_par_gae     , s23_par_gai      , s23_par_tectrl ,s23_par_drflag  , s23_par_dsflag   ,&
  !      s23_par_dpdrflag,s23_par_limrel1 ,s23_par_limrel2 , s23_par_limrel3  ,&
  !      s23_par_limrel4,s23_par_tulimit , s23_par_g0relax  , s23_par_p0relax,s23_par_dpdrtemin,&
  !      s23_par_dpdrstep,s23_par_nulimit , s23_par_nuflag   , s23_par_pinmom,&
  !     s23_par_vnmult  , s23_par_emolec   , s23_par_rec_heat,s23_par_pinqimult,s23_par_pinqiflag,&
  !      s23_par_prring0,s23_par_prring1 , s23_par_qperp34  , s23_par_qeiflag,&
  !     s23_par_chie    , s23_par_joule    , s23_par_fluxexp ,s23_par_qrec    , s23_par_dvmrel   ,&
  !      s23_par_fzrad,perp_src,part_src,ionpow_src,elecpow_src,mom_src,h_neut,&
  !     h2_neut,th_neut,magn_field,part_in_0, part_in_n,mom_in_0, mom_in_n, mom_tg_0,&
  !      mom_tg_n,pow_in_i0, pow_in_in, pow_in_e0, pow_in_en,pow_in_0,  pow_in_n,pow_tg_i0,&
  !      pow_tg_in, pow_tg_e0, pow_tg_en,pow_tg_0,  pow_tg_n,  power_ei,power_flux, del_power_flux,&
  !     mom_flux, del_mom_flux, part_flux, del_part_flux,isat_0, isat_smax,&
  !      mach_0, mach_smax,ne_0,te_0,ti_0,ne_smax,te_smax,ti_smax,gae_const,gai_const,&
  !     gtarg,momtarg,ionptarg,elecptarg,power_targ, power_lang,poweri_targ, powere_targ,&
  !      mom_targ, gamma_targ,poweri_lang, powere_lang, mom_lang, gamma_lang,te_ctrl, te_lang,&
  !      smax,spts,sbnds,mb, cb, ne,te,ti,vb,xc_cfd,xf_cfd,u_cfd,q_cfd,fract_cfd,&
  !      s23_fract2, s23_fract3, s23_fract1,s23_fract4, s23_fract4_0, s23_fract4_n,s23_fract2_0,&
  !      s23_fract3_0, s23_fract1_0,s23_fract2_n, s23_fract3_n, s23_fract1_n,s23_pinqe,&
  !     intopt,ring_type,npts,ikmax_cfd,ir_cfd,targ_cfd,pin_iter_no, s23_perp, s23_cfdfile,&
  !      s23_symm,s23_2t, finished, npts_mid,s23_adaptgrid
  !
  ! save /sol23/
  !
  real*8,public :: s23_izlen,s23_izlam,s23_izoffset,s23_momlen,s23_sm_s,s23_sm_r,s23_relax,&
       s23_maxtol,s23_rmstol,s23_minpow,s23_maxpow_0,s23_maxpow_n,s23_zhi0,s23_maxmom,&
       s23_maxmass,s23_par_ptipte,s23_par_artvisc2,s23_par_artvisc4,s23_par_dtf,&
       s23_par_dtg,s23_par_grid0,s23_par_gridexp,s23_par_ga1,s23_par_ga2,s23_par_ga3,s23_par_ga4,&
       s23_par_updtqp2,s23_par_updtqp3,s23_par_updtqpte,s23_par_garelax,s23_par_celldte,&
       s23_par_updtbcte,s23_par_updtdel0,s23_par_updtdel1,s23_par_updtdelm,s23_par_qbrelax,&
       s23_par_gridg,s23_par_grid_dx0,s23_par_gae,s23_par_gai,s23_par_limrel1,&
       s23_par_limrel2,s23_par_limrel3,s23_par_limrel4,s23_par_tulimit,s23_par_g0relax,&
       s23_par_p0relax,s23_par_dpdrtemin,s23_par_dpdrstep,s23_par_nulimit,s23_par_vnmult,&
       s23_par_emolec,s23_par_rec_heat,s23_par_pinqimult,s23_par_fzrad,part_in_0,&
       part_in_n,mom_in_0,mom_in_n,mom_tg_0,mom_tg_n,pow_in_i0,pow_in_in,pow_in_e0,pow_in_en,&
       pow_in_0,pow_in_n,pow_tg_i0,pow_tg_in,pow_tg_e0,pow_tg_en,pow_tg_0,pow_tg_n,&
       isat_0,isat_smax,mach_0,mach_smax,ne_0,te_0,ti_0,ne_smax,te_smax,ti_smax,gae_const,&
       gai_const,smax,mb,cb
  real*8,public,allocatable :: s23_powlim(:,:),s23_momlim(:,:),s23_masslim(:,:),perp_src(:),&
       part_src(:),ionpow_src(:),elecpow_src(:),mom_src(:),h_neut(:),h2_neut(:),&
       th_neut(:),magn_field(:),power_ei(:),power_flux(:,:),del_power_flux(:,:),mom_flux(:,:),&
       del_mom_flux(:,:),part_flux(:,:),del_part_flux(:,:),gtarg(:),momtarg(:),&
       ionptarg(:),elecptarg(:),power_targ(:,:),power_lang(:,:),poweri_targ(:,:),powere_targ(:,:),&
       mom_targ(:,:),gamma_targ(:,:),poweri_lang(:,:),powere_lang(:,:),mom_lang(:,:),&
       gamma_lang(:,:),te_ctrl(:,:),te_lang(:,:),spts(:),sbnds(:),ne(:),te(:),&
       ti(:),vb(:),xc_cfd(:,:,:),xf_cfd(:,:,:),u_cfd(:,:,:,:),q_cfd(:,:,:,:),fract_cfd(:,:,:)
  !
  real,public :: s23_fract2,s23_fract3,s23_fract1,s23_fract2_0,s23_fract3_0,s23_fract1_0,&
       s23_fract2_n,s23_fract3_n,s23_fract1_n,s23_pinqe,s23_fract4,s23_fract4_0,s23_fract4_n
  !
  !
  integer,public :: seedplasma,solve_half,ring_type,npts,intopt,divimp_call,bndcond_opt,&
       ir_cfd,targ_cfd,pin_iter_no,s23_perp,s23_cfdfile,s23_symm,s23_2t,npts_mid,s23_adaptgrid,&
       s23_par_adaptnr,s23_par_debugflag,s23_par_debugnr,s23_par_refresh,s23_par_itermax,&
       s23_par_updtqpit,s23_par_gaiter,s23_par_limitte,s23_par_tectrl,s23_par_drflag,&
       s23_par_dsflag,s23_par_dpdrflag,s23_par_nuflag,s23_par_pinmom,s23_par_pinqiflag,&
       s23_par_prring0,s23_par_prring1,s23_par_qperp34,s23_par_qeiflag,s23_par_chie,&
       s23_par_joule,s23_par_fluxexp,s23_par_qrec,s23_par_dvmrel
  integer,public,allocatable :: ikmax_cfd(:,:),finished(:)

  public :: allocate_mod_sol23_com,deallocate_mod_sol23_com

contains

  subroutine allocate_mod_sol23_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_sol23_com','ALLOCATE')

    call allocate_array(s23_powlim,maxpts,3,'s23_powlim',ierr)
    call allocate_array(s23_momlim,maxpts,3,'s23_momlim',ierr)
    call allocate_array(s23_masslim,maxpts,3,'s23_masslim',ierr)
    call allocate_array(perp_src,maxpts,'perp_src',ierr)
    call allocate_array(part_src,maxpts,'part_src',ierr)
    call allocate_array(ionpow_src,maxpts,'ionpow_src',ierr)
    call allocate_array(elecpow_src,maxpts,'elecpow_src',ierr)
    call allocate_array(mom_src,maxpts,'mom_src',ierr)
    call allocate_array(h_neut,maxpts,'h_neut',ierr)
    call allocate_array(h2_neut,maxpts,'h2_neut',ierr)
    call allocate_array(th_neut,maxpts,'th_neut',ierr)
    call allocate_array(magn_field,maxpts,'magn_field',ierr)
    call allocate_array(power_ei,maxpts,'power_ei',ierr)
    call allocate_array(power_flux,maxpts,3,'power_flux',ierr)
    call allocate_array(del_power_flux,maxpts,3,'del_power_flux',ierr)
    call allocate_array(mom_flux,maxpts,3,'mom_flux',ierr)
    call allocate_array(del_mom_flux,maxpts,3,'del_mom_flux',ierr)
    call allocate_array(part_flux,maxpts,3,'part_flux',ierr)
    call allocate_array(del_part_flux,maxpts,3,'del_part_flux',ierr)
    call allocate_array(gtarg,3,'gtarg',ierr)
    call allocate_array(momtarg,3,'momtarg',ierr)
    call allocate_array(ionptarg,3,'ionptarg',ierr)
    call allocate_array(elecptarg,3,'elecptarg',ierr)
    call allocate_array(power_targ,maxpts,3,'power_targ',ierr)
    call allocate_array(power_lang,maxpts,3,'power_lang',ierr)
    call allocate_array(poweri_targ,maxpts,3,'poweri_targ',ierr)
    call allocate_array(powere_targ,maxpts,3,'powere_targ',ierr)
    call allocate_array(mom_targ,maxpts,3,'mom_targ',ierr)
    call allocate_array(gamma_targ,maxpts,3,'gamma_targ',ierr)
    call allocate_array(poweri_lang,maxpts,3,'poweri_lang',ierr)
    call allocate_array(powere_lang,maxpts,3,'powere_lang',ierr)
    call allocate_array(mom_lang,maxpts,3,'mom_lang',ierr)
    call allocate_array(gamma_lang,maxpts,3,'gamma_lang',ierr)
    call allocate_array(te_ctrl,maxpts,3,'te_ctrl',ierr)
    call allocate_array(te_lang,maxpts,3,'te_lang',ierr)
    call allocate_array(spts,maxpts,'spts',ierr)
    call allocate_array(sbnds,0,'sbnds',maxpts,ierr)
    call allocate_array(ne,maxpts,'ne',ierr)
    call allocate_array(te,maxpts,'te',ierr)
    call allocate_array(ti,maxpts,'ti',ierr)
    call allocate_array(vb,maxpts,'vb',ierr)
    call allocate_array(xc_cfd,1,maxpts,1,3,0,maxpts,'xc_cfd',ierr)
    call allocate_array(xf_cfd,1,maxpts,1,3,0,maxpts,'xf_cfd',ierr)
    call allocate_array(u_cfd,1,maxpts,1,3,1,4,0,maxpts,'u_cfd',ierr)
    call allocate_array(q_cfd,1,maxpts,1,3,1,4,0,maxpts,'q_cfd',ierr)
    call allocate_array(fract_cfd,maxpts,2,4,'fract_cfd',ierr)
    call allocate_array(ikmax_cfd,maxpts,3,'ikmax_cfd',ierr)
    call allocate_array(finished,maxpts,'finished',ierr)

  end subroutine allocate_mod_sol23_com


  subroutine deallocate_mod_sol23_com
    implicit none

    call pr_trace('mod_sol23_com','DEALLOCATE')

    if (allocated(s23_powlim)) deallocate(s23_powlim)
    if (allocated(s23_momlim)) deallocate(s23_momlim)
    if (allocated(s23_masslim)) deallocate(s23_masslim)
    if (allocated(perp_src)) deallocate(perp_src)
    if (allocated(part_src)) deallocate(part_src)
    if (allocated(ionpow_src)) deallocate(ionpow_src)
    if (allocated(elecpow_src)) deallocate(elecpow_src)
    if (allocated(mom_src)) deallocate(mom_src)
    if (allocated(h_neut)) deallocate(h_neut)
    if (allocated(h2_neut)) deallocate(h2_neut)
    if (allocated(th_neut)) deallocate(th_neut)
    if (allocated(magn_field)) deallocate(magn_field)
    if (allocated(power_ei)) deallocate(power_ei)
    if (allocated(power_flux)) deallocate(power_flux)
    if (allocated(del_power_flux)) deallocate(del_power_flux)
    if (allocated(mom_flux)) deallocate(mom_flux)
    if (allocated(del_mom_flux)) deallocate(del_mom_flux)
    if (allocated(part_flux)) deallocate(part_flux)
    if (allocated(del_part_flux)) deallocate(del_part_flux)
    if (allocated(gtarg)) deallocate(gtarg)
    if (allocated(momtarg)) deallocate(momtarg)
    if (allocated(ionptarg)) deallocate(ionptarg)
    if (allocated(elecptarg)) deallocate(elecptarg)
    if (allocated(power_targ)) deallocate(power_targ)
    if (allocated(power_lang)) deallocate(power_lang)
    if (allocated(poweri_targ)) deallocate(poweri_targ)
    if (allocated(powere_targ)) deallocate(powere_targ)
    if (allocated(mom_targ)) deallocate(mom_targ)
    if (allocated(gamma_targ)) deallocate(gamma_targ)
    if (allocated(poweri_lang)) deallocate(poweri_lang)
    if (allocated(powere_lang)) deallocate(powere_lang)
    if (allocated(mom_lang)) deallocate(mom_lang)
    if (allocated(gamma_lang)) deallocate(gamma_lang)
    if (allocated(te_ctrl)) deallocate(te_ctrl)
    if (allocated(te_lang)) deallocate(te_lang)
    if (allocated(spts)) deallocate(spts)
    if (allocated(sbnds)) deallocate(sbnds)
    if (allocated(ne)) deallocate(ne)
    if (allocated(te)) deallocate(te)
    if (allocated(ti)) deallocate(ti)
    if (allocated(vb)) deallocate(vb)
    if (allocated(xc_cfd)) deallocate(xc_cfd)
    if (allocated(xf_cfd)) deallocate(xf_cfd)
    if (allocated(u_cfd)) deallocate(u_cfd)
    if (allocated(q_cfd)) deallocate(q_cfd)
    if (allocated(fract_cfd)) deallocate(fract_cfd)
    if (allocated(ikmax_cfd)) deallocate(ikmax_cfd)
    if (allocated(finished)) deallocate(finished)

  end subroutine deallocate_mod_sol23_com

end module mod_sol23_com