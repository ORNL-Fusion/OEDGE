module mod_comtor
  
  use mod_params
  
  implicit none
  
  
  private
  
  !c -*-fortran-*-
  !c
  !c     jdemod - cprint moved to global_options
  !c
  !      common /comtor/ ca,    caw,   cl,    crmb,  ctbin, cgtin1,cltin1,
  !     >                ctboul,cltoul,cnbin, cgnin1,clnin1,cnboul,clnoul,
  !     >                cvin,  crdxo, crdd,  crmi,  cxsc,  cfimp, cqs,
  !     >                cysc,  ctemsc,ctimsc,cebd,  ctgas, cemaxf,cengsc,
  !     >                cizb,  cion,  cizsc, catin, canin, ctresh,cizeff,
  !     >                cnhc,  cnho,  clamhx,clamhy,cxnear,cynear,corect,
  !     >                ctbins,cnbins,ccut,  cpfir, cpsub, cpsc,  csnorm,
  !     >    cdpol, cplsma,ceyin,cvhyin,cstepn,cizset,czenh,debugt,ptracl,
  !     >    ceyout,cvhout,ciseed,cdifop,ctwol, cystag,cbombf,ptracs,
  !     >    cbombz,cirf,  csef,  cgtin2,cgnin2,cltin2,clnin2,cstept,
  !     >    cein2, csolef,cvcx,cymfs, cthetb,csintb,cvpout,cvxmin,
  !     >    cxspls,cnspl, ctboug,cltoug,cnboug,coni  ,cprul,cvxmax,
  !     >    clnoug,clarmr,cki   ,ctbi,  chalfl,debugn,cko, svybar,svyacc,
  !     >    cono,cymflg,cftcut,cnba  ,cgamma,canal ,cap, vpv0,vpalph,
  !     >    cwl,cstepl,debugl,clfact,cyfar ,crdxi ,ctsub, cvpcut,dpbeta,
  !     >    cqpl,cqsl, cioptj, cpco,ctibin,cltiin1,cgtiin1,catiin,cdperp,
  !     >    dpalph,
  !     >    ctiboul,cltioul,ctiboug,cltioug,cgtiin2,cltiin2,cvpopt,
  !     >    ctibins,cprob,ctichg,clpd,lpdion,lpdcum,cvsa,rlc,cdpstp,
  !     >    cbrk,rledge7,calphe,cbetai,qmultp,qmults,ciangn,
  !     >    ycfadd,cspumax,cfbgff,c3halfl,csvymin,cmaxgens,cdcalc,cvpol
  !c
  !      real            cvpout,cvxmin,cvxmax
  !      real            svybar(-maxqxs:maxqxs),svyacc(-maxqxs:maxqxs)
  !      real            ca,caw,cl,crmb,ctbin,cgtin1,cltin1,ctboul,cltoul
  !      real            cnbin,cgnin1,clnin1,cnboul,clnoul,cfimp,crdd,cdpol
  !      real            cvin,crmi,cxsc,cysc,ctemsc,ctimsc,cebd,ctgas
  !      real            cemaxf,cengsc,catin,canin,ctresh,cqs(maxins,2)
  !      real            cnhc,cnho,clamhx,clamhy,cxnear,cynear,ccut,cstepn
  !      real            cpfir,cpsub,cpsc,csnorm,cplsma,ceyin,cvhyin,czenh
  !      real            ceyout,cvhout,ctwol,cystag,cirf,csef,csolef,cvcx
  !      real            cgtin2,ciangn,cgnin2,cltin2,clnin2,cein2,cthetb,
  !     >                csintb
  !      real            ctbins(maxins,2),cnbins(maxins,2),cymfs(maxins,3)
  !      real            cxspls(0:maxins+1),ctboug,cltoug,cnboug,clnoug
  !      real            clarmr,cki,cko,ctbi,cprul,chalfl,coni,cono,cftcut
  !      real            cnba,cgamma,canal,cap,cwl,cstepl,clfact,cyfar
  !      real            crdxo(3),crdxi(3),ctsub,cqpl,cqsl,cpco,cprob
  !      real            ctibin,cltiin1,cgtiin1,catiin,ctiboul,cltioul
  !      real            ctiboug,cltioug,cgtiin2,cltiin2,ctibins(maxins,2)
  !      real            lpdion(maxlpd,2),lpdcum(maxlpd)
  !      real            cvsa(maxins,3),rledge7
  !      real            calphe(maxizs),cbetai(maxizs)
  !      real            ycfadd,cspumax,cfbgff,c3halfl,csvymin
  !      real            qmultp,qmults
  !      real            cvpol,rlc, cdpstp
  !      real            dpalph,dpbeta,vpv0,vpalph,cvpcut
  !      real            ptracs(maxlen,maxt,2)
  !c
  !      integer         cioptj
  !      integer         cmaxgens,cdcalc,cdperp,cvpopt
  !      integer         cizb,cion,cizsc,cizeff,cizset,ciseed,cdifop,corect
  !      integer         cbombf,cbombz,cnspl,cymflg,clpd,cbrk,cstept
  !      integer         ptracl(maxt)
  !      logical         debugn,debugl,debugt,ctichg
  !c
  !c     this file contains declarations for the unstructured
  !c     input values. if it becomes too unwieldy this file will
  !c     be split into separate common blocks for different
  !c     unstructured input values.
  !c
  !      common /sputdat/ csputopt,extra_sputter_angle,
  !     >                 cchemopt, const_yield, impact_energy_opt,
  !     >                 init_y_coord,cselfs,
  !     >                 ss_nymfs,ss_cymfs
  !      integer csputopt,cchemopt,impact_energy_opt,cselfs
  !      integer ss_nymfs
  !      real    extra_sputter_angle,init_y_coord
  !      real    const_yield
  !      real    ss_cymfs(maxins,3)
  !c
  !c     gradient multipliers
  !c
  !      common /gradmult/ ntig,nteg,nnbg,tmig,tmeg,mnbg
  !      integer ntig,nteg,nnbg
  !      real    tmig(maxins,2),tmeg(maxins,2),mnbg(maxins,2)
  !c
  !c     scaling factor
  !c
  !      common /scaling_factor/ absfac
  !c
  !c     in divimp absfac is a real - should lim be the same? since all use of absfac should derive
  !c     from this common block include - the declaration should remain consistent within lim - however
  !c     using real*8 here runs the risk that imported code might have an issue.
  !c
  !      real*8 absfac
  !c      real absfac
  !c
  !c     common block for generic lim related unstructured input
  !c
  !      common /lim_unstruc/ shear_short_circuit_opt,calc_3d_power,
  !     >                     extfluxopt,nextfluxdata,extfluxdata,
  !     >                     vpflow_3d
  !c
  !      integer shear_short_circuit_opt,calc_3d_power,extfluxopt,
  !     >        nextfluxdata
  !      real extfluxdata(maxins,3),vpflow_3d
  !c
  !c     common block for out related unstructured input
  !c
  !      common /out_unstruc/ new_absfac,erosion_scaling_opt
  !c
  !      real*8 new_absfac
  !c
  !c      real new_absfac
  !c
  !      integer erosion_scaling_opt
  !
  
  
  
  !      common /comtor/ ca,    caw,   cl,    crmb,  ctbin, cgtin1,cltin1,
  !     >                ctboul,cltoul,cnbin, cgnin1,clnin1,cnboul,clnoul,
  !     >                cvin,  crdxo, crdd,  crmi,  cxsc,  cfimp, cqs,
  !     >                cysc,  ctemsc,ctimsc,cebd,  ctgas, cemaxf,cengsc,
  !     >                cizb,  cion,  cizsc, catin, canin, ctresh,cizeff,
  !     >                cnhc,  cnho,  clamhx,clamhy,cxnear,cynear,corect,
  !     >                ctbins,cnbins,ccut,  cpfir, cpsub, cpsc,  csnorm,
  !     >    cdpol, cplsma,ceyin,cvhyin,cstepn,cizset,czenh,debugt,ptracl,
  !     >    ceyout,cvhout,ciseed,cdifop,ctwol, cystag,cbombf,ptracs,
  !     >    cbombz,cirf,  csef,  cgtin2,cgnin2,cltin2,clnin2,cstept,
  !     >    cein2, csolef,cvcx,cymfs, cthetb,csintb,cvpout,cvxmin,
  !     >    cxspls,cnspl, ctboug,cltoug,cnboug,coni  ,cprul,cvxmax,
  !     >    clnoug,clarmr,cki   ,ctbi,  chalfl,debugn,cko, svybar,svyacc,
  !     >    cono,cymflg,cftcut,cnba  ,cgamma,canal ,cap, vpv0,vpalph,
  !     >    cwl,cstepl,debugl,clfact,cyfar ,crdxi ,ctsub, cvpcut,dpbeta,
  !     >    cqpl,cqsl, cioptj, cpco,ctibin,cltiin1,cgtiin1,catiin,cdperp,
  !     >    dpalph,
  !     >    ctiboul,cltioul,ctiboug,cltioug,cgtiin2,cltiin2,cvpopt,
  !     >    ctibins,cprob,ctichg,clpd,lpdion,lpdcum,cvsa,rlc,cdpstp,
  !     >    cbrk,rledge7,calphe,cbetai,qmultp,qmults,ciangn,
  !     >    ycfadd,cspumax,cfbgff,c3halfl,csvymin,cmaxgens,cdcalc,cvpol
  !c
  real,public :: cvpout,cvxmin,cvxmax,ca,caw,cl,crmb,ctbin,cgtin1,cltin1,ctboul,cltoul,&
       cnbin,cgnin1,clnin1,cnboul,clnoul,cfimp,crdd,cdpol,cvin,crmi,cxsc,cysc,ctemsc,&
       ctimsc,cebd,ctgas,cemaxf,cengsc,catin,canin,ctresh,cnhc,cnho,clamhx,clamhy,cxnear,&
       cynear,ccut,cstepn,cpfir,cpsub,cpsc,csnorm,cplsma,ceyin,cvhyin,czenh,ceyout,&
       cvhout,ctwol,cystag,cirf,csef,csolef,cvcx,cgtin2,ciangn,cgnin2,cltin2,clnin2,cein2,&
       cthetb,csintb,ctboug,cltoug,cnboug,clnoug,clarmr,cki,cko,ctbi,cprul,chalfl,coni,&
       cono,cftcut,cnba,cgamma,canal,cap,cwl,cstepl,clfact,cyfar,ctsub,cqpl,cqsl,cpco,&
       cprob,ctibin,cltiin1,cgtiin1,catiin,ctiboul,cltioul,ctiboug,cltioug,cgtiin2,cltiin2,&
       rledge7,ycfadd,cspumax,cfbgff,c3halfl,csvymin,qmultp,qmults,cvpol,rlc,cdpstp,&
       dpalph,dpbeta,vpv0,vpalph,cvpcut,ctimsc_win
  real,public ,allocatable:: svybar(:),svyacc(:),cqs(:,:),ctbins(:,:),cnbins(:,:),cymfs(:,:),&
       cxspls(:),crdxo(:),crdxi(:),ctibins(:,:),lpdion(:,:),lpdcum(:),cvsa(:,:),&
       calphe(:),cbetai(:),ptracs(:,:,:)
  
  
  real,public :: te_prof_shift,ti_prof_shift,ne_prof_shift
  real,public :: te_prof_mult,ti_prof_mult,ne_prof_mult
  
  
  integer,public:: cioptj,cmaxgens,cdcalc,cdperp,cvpopt,cizb,cion,cizsc,cizeff,cizset,&
       ciseed,cdifop,corect,cbombf,cbombz,cnspl,cymflg,clpd,cbrk,cstept

  !     jdemod - remove cdwelt_sum option functionality because it isn't
  !              physically meaningful.                  
  !,cdwelt_sum

  integer,public,allocatable:: ptracl(:)
  
  logical,public:: debugn,debugl,debugt,ctichg
  !c
  !c     this file contains declarations for the unstructured
  !c     input values. if it becomes too unwieldy this file will
  !c     be split into separate common blocks for different
  !c     unstructured input values.
  !c
  !      common /sputdat/ csputopt,extra_sputter_angle,
  !     >                 cchemopt, const_yield, impact_energy_opt,
  !     >                 init_y_coord,cselfs,
  !     >                 ss_nymfs,ss_cymfs
  
  integer,public:: csputopt,cchemopt,impact_energy_opt,cselfs,ss_nymfs
  integer, public :: yieldsw
  real,public:: extra_sputter_angle,init_y_coord,const_yield
  real,public,allocatable:: ss_cymfs(:,:)
  !c
  !c     gradient multipliers
  !c
  !      common /gradmult/ ntig,nteg,nnbg,tmig,tmeg,mnbg
  integer,public :: ntig,nteg,nnbg
  real,public,allocatable:: tmig(:,:),tmeg(:,:),mnbg(:,:)
  !c
  !c     scaling factor
  !c
  !      common /scaling_factor/ absfac
  !c
  !c     in divimp absfac is a real - should lim be the same? since all use of absfac should derive
  !c     from this common block include - the declaration should remain consistent within lim - however
  !c     using real*8 here runs the risk that imported code might have an issue.
  !c
  real*8,public :: absfac
  !c      real absfac
  !c
  !c     common block for generic lim related unstructured input
  !c
  !      common /lim_unstruc/ shear_short_circuit_opt,calc_3d_power,
  !     >                     extfluxopt,nextfluxdata,extfluxdata,
  !     >                     vpflow_3d
  !c
  integer,public:: shear_short_circuit_opt,calc_3d_power,extfluxopt,nextfluxdata
  real,public:: vpflow_3d
  real,public,allocatable:: extfluxdata(:,:)
  !c
  !c     common block for out related unstructured input
  !c
  !      common /out_unstruc/ new_absfac,erosion_scaling_opt
  !c
  real*8,public:: new_absfac
  !c
  !c      real new_absfac
  !c
  integer,public:: erosion_scaling_opt
  
  !c    sazmod - switch to vary absorbing boundary that affects plasma solution
  !c             and the accompanying values.
  !integer, public:: vary_absorb, ix_step1, ix_step2
  !real,    public:: yabsorb1a_step, yabsorb2a_step, xabsorb1a_step, xabsorb2a_step
  
  ! add in regions for varying radial dperp values.
  integer, public:: dperp_reg_switch
  real,    public:: dperp_reg1,dperp_reg2,dperp_reg3,dperp_reg4
  
  ! if we don't care about the .raw file save some time and skip it.
  integer, public:: skip_raw
  
  ! modify the velplasma in the step region (right half only right now).
  real,    public:: mod_v_fact
  
  ! option to chose from an exponential distribution in the y direction (3d only).
  integer, public:: choose_exp
  real,    public:: choose_exp_lambda,choose_exp_fact
  
  ! overall scaling factor to apply to the background plasma velocity.
  real,    public:: vel_mod
  
  ! switch to load in fully customizable 2D absorbing boundary.
  integer, public:: vary_2d_bound, bounds_rows, bounds_cols
  
  ! for the DIVIMP injection probability.
  integer, public:: ndivimp_probs
  real, public, allocatable :: divimp_probs(:,:), yinj_cdf(:)
  
  ! For mixed material model.
  real, public :: frac_c
  integer, public :: mm_usage
  
  public :: allocate_mod_comtor, deallocate_mod_comtor
  
  
contains
  
  subroutine allocate_mod_comtor
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(ptes,maxnxs,'ptes',ierr)

    call allocate_array(svybar,-maxqxs,'svybar',maxqxs,ierr)
    call allocate_array(svyacc,-maxqxs,'svyacc',maxqxs,ierr)
    call allocate_array(cqs,maxins,2,'cqs',ierr)
    call allocate_array(ctbins,maxins,2,'ctbins',ierr)
    call allocate_array(cnbins,maxins,2,'cnbins',ierr)
    call allocate_array(cymfs,maxins,3,'cymfs',ierr)
    call allocate_array(cxspls,0,'cxspls',maxins+1,ierr)
    call allocate_array(crdxo,3,'crdxo',ierr)
    call allocate_array(crdxi,3,'crdxi',ierr)
    call allocate_array(ctibins,maxins,2,'ctibins',ierr)
    call allocate_array(lpdion,maxlpd,2,'lpdion',ierr)
    call allocate_array(lpdcum,maxlpd,'lpdcum',ierr)
    call allocate_array(cvsa,maxins,3,'cvsa',ierr)
    call allocate_array(calphe,maxizs,'calphe',ierr)
    call allocate_array(cbetai,maxizs,'cbetai',ierr)
    call allocate_array(ptracs,maxlen,maxt,2,'ptracs',ierr)
    call allocate_array(ptracl,maxt,'ptracl',ierr)
    call allocate_array(ss_cymfs,maxins,3,'ss_cymfs',ierr)
    call allocate_array(tmig,maxins,2,'tmig',ierr)
    call allocate_array(tmeg,maxins,2,'tmeg',ierr)
    call allocate_array(mnbg,maxins,2,'mnbg',ierr)
    call allocate_array(extfluxdata,maxins,3,'extfluxdata',ierr)
    
    ! Just going with maxnys here, since it should be more than large 
    ! enough, but pretty much arbitrary. Ideally it would just 
    ! be ndivimp_probs.
    call allocate_array(divimp_probs, maxnys, 2, 'divimp_probs', ierr)
    call allocate_array(yinj_cdf, maxnys, 'yinj_cdf', ierr)

  end subroutine allocate_mod_comtor
  
  
  subroutine deallocate_mod_comtor
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

    if (allocated(svybar)) deallocate(svybar)
    if (allocated(svyacc)) deallocate(svyacc)
    if (allocated(cqs)) deallocate(cqs)
    if (allocated(ctbins)) deallocate(ctbins)
    if (allocated(cnbins)) deallocate(cnbins)
    if (allocated(cymfs)) deallocate(cymfs)
    if (allocated(cxspls)) deallocate(cxspls)
    if (allocated(crdxo)) deallocate(crdxo)
    if (allocated(crdxi)) deallocate(crdxi)
    if (allocated(ctibins)) deallocate(ctibins)
    if (allocated(lpdion)) deallocate(lpdion)
    if (allocated(lpdcum)) deallocate(lpdcum)
    if (allocated(cvsa)) deallocate(cvsa)
    if (allocated(calphe)) deallocate(calphe)
    if (allocated(cbetai)) deallocate(cbetai)
    if (allocated(ptracs)) deallocate(ptracs)
    if (allocated(ptracl)) deallocate(ptracl)
    if (allocated(ss_cymfs)) deallocate(ss_cymfs)
    if (allocated(tmig)) deallocate(tmig)
    if (allocated(tmeg)) deallocate(tmeg)
    if (allocated(mnbg)) deallocate(mnbg)
    if (allocated(extfluxdata)) deallocate(extfluxdata)

  end subroutine deallocate_mod_comtor
  
  
  
end module mod_comtor
