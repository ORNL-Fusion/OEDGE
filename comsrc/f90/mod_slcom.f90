module mod_slcom
  use debug_options
  implicit none

  
  !
  !   for general use: maxasc  = 5000
  !                    maxnas2 = 2000
  !
  !   for large vacuum grids: maxasc  = 20000 (presently)
  !                           maxnas2 = 10000
  !
  !   the difference in the size of the divimp/out executables for the above
  !   array sets is about 80 mb
  !
  !...maxnas2 should be the same as nlim in parmusr
  !...nmomcha should be the same as nmomcha in eirene (ha! just checked and
  !   nmomcha is not used anywhere in divimp -- really maxdata is 12+nmomcha)
  !
  ! parameters:
  !
  !     -*-fortran-*-
  
  !             number of momentum source channels.  this must be the same
  !             as nmomcha in eirene:
  integer,public :: nonvirtual,virtual,before,after,boundary,iklo,ikhi,permanent,temporary,&
       side14,side23,side34,cell1,total,slout,erout,eirin,eirout,dtout,pinout,pinout2,&
       pinout3,logout,pinout4,tartotar,tartowal,waltotar,waltowal,null,complete,&
       partial,maxstep,maxriter,maxerr,maxnas,maxasd,maxdis,maxdata,maxbgk,maxnas2,maxasd2,&
       maxass,maxtor,maxfnt,maxstrdat,maxnas3,h_balpha,h_bgamma,h_bbeta,nmomcha,maxtally,&
       maxflex,maxstrata,maxpuff,maxvacregion,maxiontime,maxbin,sol1,pfz,core,cmod,&
       diiid,jet,mast
  
  real,public :: nil,center
  
  
  !...    room for 15 momentum loss channels in maxdata, but only reading 1
  !       channel currently:
  parameter (nonvirtual  = 1 , virtual  = 2 , before   = 3 , after     =  4,boundary    = -1,&
       permanent   =  1, temporary = 2 ,iklo        = 1 , ikhi      = 2 ,side14      = 1 ,&
        side23    = 2 , total    = 3 , side34    =  4,cell1       = 3 ,tartotar    = 1 ,&
        tartowal = 2 , waltotar = 3 , waltowal  =  4,slout       = 50 , erout    = 50,&
        eirout   = 81, eirin     = 80,dtout       = 7  , pinout   = 88, pinout2  = 89,&
        logout    = 90,pinout3     = 94 , pinout4  = 95,null        = 0  ,complete    = 1  ,&
        partial   = 2  ,maxstep     = 100, maxriter  = 100, maxerr   = 10,&
       maxasd      = 30 , maxnas    = 100, maxdis   = 16,maxpuff     = 20 ,maxfnt      = 10 ,&
        maxdata   =  26  , nmomcha  =  1,maxbgk      = 30 , maxnas2   = 10000,&
        maxasd2  =  5,maxnas3   = 1000 ,maxass      = 10 , maxtor    = 50   ,h_balpha    =  1 ,&
        h_bgamma  = 2    , h_bbeta = 3      ,maxflex     = 10 , maxstrata = 10   ,&
        maxvacregion = 10,maxiontime  = 40 , maxbin    = 30   , maxtally     = 10,maxstrdat   = 50 ,&
       sol1 = 1, pfz = 2, core = 3,cmod = 1, diiid = 2, jet = 3, mast = 4)
  
  parameter (nil        =  0.0,center     = -1.0)
  integer,public :: h_ion1,h_ion2,h_ion3,h_ion4,h_atm1,h_atm2,h_atm3,h_atm4,h_mol1,&
       h_mol2,h_mol3,h_mol4,h_mp1,h_mp2,h_mp3,h_mp4,h_mp5,h_mp6,h_mp7,h_mp8,h_mp9
  
  !...  block cleared out -sl, 28.07.09
  parameter (h_ion1 =  1,h_ion2 =  2,h_ion3 =  3,h_ion4 =  4,h_atm1 =  5,h_atm2 =  6,&
       h_atm3 =  7,h_atm4 =  8,h_mol1 =  9,h_mol2 = 10,h_mol3 = 11,h_mol4 = 12,h_mp1  = 13,&
       h_mp2  = 14,h_mp3  = 15,h_mp4  = 16,h_mp5  = 17,h_mp6  = 18,h_mp7  = 19,h_mp8  = 20,&
       h_mp9  = 21)
  ! common /slcom1/intpmk,dumpae1,dumpai1,dumpei1,osm_store,osm_preopt,osm_recopt,osm_model,&
  !     osmbulkn,osmbulkte,osmbulkti,osmbulkv,optuser,osmikshift,thesis,tuprati,&
  !     tuprat
  
  ! save /slcom1/
  ! common /slcom2/idring,supflx,ringtype,ikbrks,ikbrke,nbr,nbreak,nbs,irorg2,virtag,&
  !     cgridtype,irbreak,ikbreak,ikti2,ikto2,idtarg,ndsdiv,cgridst,vpolyp,vpolmin,eirgrid,&
  !     eirneut,eirgeom,eirdebug,eirtime,eirdata,eiradd,eirmat1,eirmat2,eirtarg,wallpt2,&
  !     eirtemp1,eirtemp2,eirbgk,eircxd2,eiropacity,eirnstrata,eirstrata,eirph2,eirfuji,&
  !     eiralloc,eiracx,eird2ion,eird2dis,eirspdatmode,eirtrim,eirpgdat,eirnpgdat,eirasdat,&
  !     eirnasdat,eirspdat,eirnspdat,eirres,eirnres,eirpuff,eirnpuff,eirpmode,eirntracks,&
  !     eirniter,eirrinteg,eireinteg,eirermin,eirtrc1,eirtrc2,eird2fu,eirzaa,eiraout,&
  !     eirnaout,eirtorfrac,eirnsdtor,eirsdtor,eirsdvol,eirsrcmul,eirniontime,eiriontime,&
  !     eirntally,eirscale,eirnsection,eirntrans,eirtrans,eirfluxt,eirnttra,eirnstrdat,&
  !     eirstrdat,eirnstrai,eirnatmi,eirnmoli,eirstrlis,eirntorseg,eirnlimi,haddi2,&
  !     eirdtimv,eirtrin,eirtriik,eirtriir,eirpho1,eirpho2,eirtri,eirntri,eirphoton,setmach,&
  !     muldensity,addtemp,osm_ionfnt,fradd,bkadd,cwcopt,cmodopt,stagopt,stopopt,stopopt2,&
  !     stopopt3,rflexopt,iflexopt,virloc,outwall,outplasma,outcell,outpoly,outgeom,&
  !     outtarget,outpin,pincode,outmode,out_source,out_plasma,out_geom,rho
  
  ! save /slcom2/
  ! common /eirenecom/ eirtally
  
  ! save /eirenecom/
  ! common /tagcom/ tagpinatom,tagmulrec
  
  ! save /tagcom/
  !     .  relmode,relreset,relseed,relexpdat,
  ! common /slcom3/prb_num,prb_te,prb_ti,prb_ne,prb_rho,prb_r,prb_z,prb_shift,prb_align,&
  !     prb_tag,dip_te,dip_ti,dip_ne,dip_s,dip_v,dip_ik,prp_te,prp_ti,prp_ne,tarsource,&
  !     tarshift,tarshiftopt,haldata,shot,slice_time,hal_num,hal_ang,hal_val,hal_tag,&
  !     grd_minpl,grd_refine,grd_range,grd_shift,grd_thresh,rel_opt  ,rel_frac,rel_count,&
  !     rel_niter,rel_nstep ,rel_iter  ,rel_step ,rel_pace,rel_data ,rel_ndata,rel_bound1,&
  !     rel_bound2,relmode,relreset,relexpdat,rel_viter,osm_sympt,osm_symopt,osm_powopt,&
  !     osm_dp1,osm_dp2,osm_dp3,osm_dp4,osm_dp5,osm_dp6,osm_code,rel_tol  ,rel_symfr,&
  !     rel_prbfr , rel_mfsp,rel_qemul,osm_peimul,osmikp,rel_deltati,rel_dirtg,rel_hte,&
  !     rel_hti,rel_hne,rel_hproe  ,rel_hproi,osmmock
  
  ! save /slcom3/
  ! common /slcom3_5/ relseed,ringuppr
  
  ! save /slcom3_5/
  !     .  ringuppr,ringmodel,osms28,osmns28,s22pfzion,
  ! common /slcom4/slffric,osm_watch1,osm_watch2,osm_probe,osm_nflag,ikbound,ikfluid,&
  !     lpdati2,lpdato2   ,nlpdati2,nlpdato2,pinion2 ,pinqe2   ,pinrec2 ,pinqi2  ,pinenz2  ,&
  !     pinenm2,pinmol2 ,pinena2  ,pinmp2  ,pinatom2,pinalpha2,pinionz2,pinz02   ,&
  !     pinline ,pinline2 ,pindata ,pindata2,mulrec,mulion,mulqei,mulqer  ,mulrecl,pinior  ,&
  !     pinmpr   ,pinqir  ,pinqer  ,pinior2 ,pinmpr2  ,pinqir2 ,pinqer2 ,pinbgk  ,pinbgk2  ,&
  !     pinasd  ,indasd,pinmoi  ,pinstrata,pinploss,pinalgv,pinioncomp,adp_opt,adp_region,&
  !     adp_upper,adp_lower,osm_mode,osm_range,s21_ndatai,s21_ndatao,s21_datai,&
  !     s21_datao,s22_ir,s22_l1,s22_nr,s22_maxs,s22_tr,s22_midnks,s21_mode,ringmodel,osms28,&
  !     osmns28,s22pfzion,osms28over,osmns28over,s28sym
  
  ! save /slcom4/
  ! common /slcom5/osm_matchp,osm_matchs,osm_sfrac,osm_nfnt,osm_testep,vacnseg,vacseg,&
  !     vacregion,vacnpla,vacpla,firstshape,machine
  
  ! save /slcom5/
  ! common /cioptcom/ orgcioptf,orgcioptg
  ! save /cioptcom/
  
  
  integer,public :: orgcioptf,orgcioptg
  ! common /discom/ disindex,dislev
  ! save /discom/
  integer,public,allocatable :: disindex(:)
  
  
  real,public,allocatable :: dislev(:)
  logical,public :: firstshape
  
  !
  ! constants:
  !
  integer,public :: machine
  real,public :: bzc
  
  parameter (bzc=1.38e-23)
  !
  ! new dperp options:
  !
  integer,public :: outmode
  !
  ! cell width calculation option:
  !
  integer,public :: optuser
  !
  ! broken grid and degnerate grid variables:
  !
  integer,public :: cwcopt
  integer,public :: cgridtype
  integer,public :: nbr
  integer,public,allocatable :: ikbrks(:,:),ikbrke(:,:),nbreak(:)
  integer,public :: nbs
  integer,public,allocatable :: irorg2(:),idring(:),supflx(:,:),ringtype(:)
  integer,public :: irbreak
  integer,public,allocatable :: ikti2(:),ikto2(:),idtarg(:),ikbreak(:)
  integer,public :: ndsdiv
  integer,public,allocatable :: virtag(:,:)
  !
  ! outputgrid related options:
  !
  integer,public :: vpolyp,vpolmin
  !
  ! eirene related options:
  !
  integer,public :: outwall,outplasma,outcell,outpoly,outgeom,outtarget,outpin,out_source,&
       out_plasma,out_geom
  character*128,public :: eirtally(maxtally,10)
  !... replace the 5000:
  integer,public :: cgridst,eirgrid,eirneut,eirgeom,fradd,bkadd,cmodopt,stagopt,stopopt,&
       stopopt2,stopopt3,eirdebug,eirtime,eirdata,eiradd,eirniter,eirtrim,eirmat1,&
       eirmat2,pincode,eirnpgdat,eirbgk,eiropacity,eircxd2,eirnasdat,eirnspdat,eirph2,eirnpuff,&
       eirpmode,eirnstrata,eirnres,eiracx,eird2ion,eird2dis,eirnaout,eirnsdtor,&
       eirfuji,eirniontime,eirntally,eirntrans,eirnstrdat,eirnatmi,eirnmoli,eirnstrai,eirnttra,&
       eirntorseg,eirtrin
  integer,public,allocatable :: eirntracks(:),eiraout(:,:),eirstrlis(:),eirnlimi(:),&
       eirtriik(:),eirtriir(:)
  integer,public :: eirnsection,haddi2,eirntri,eirphoton
  integer,public,allocatable :: iflexopt(:)
  
  
  !      real    rflexopt(maxflex),eirres(20,maxpiniter)
  !
  ! eirpgdat:
  !   1 - gauge code, corresponding to a particular experimental gauge
  !         c-mod:
  !            001
  !         diii-d:
  !            101   pbf1
  !            102   pbf2
  !            103   pbf3
  !            104   pv1
  !            105   pr2
  !            106   vplows
  !            107   pcm105baf
  !            108   pcm240tor
  !   2 - x coordinate for the gauge (m)
  !   3 - y coordinate (m)
  !   4 - toroidal location (degrees)
  !   5 - radius of gauge cylinder (m)
  !   6 - length of gauge cylinder (m)
  !   7 - experimental data
  !   8 - volume of gauge (m3)
  !
  !   9 - d2 energy density, raw data from eirene (ev cm-3)
  !  10 - d2 energy density ( ... mks)
  !  11 - d2 energy denisty, relaxed (... mks)
  !
  !  12 - d energy density, raw data from eirene (ev cm-3)
  !  13 - d energy density ( ... mks)
  !  14 - d energy denisty, relaxed (... mks)
  !
  !
  real,public :: eirrinteg,eireinteg,eirermin,eirsrcmul,eiralloc,setmach,muldensity,&
       addtemp,eirdtimv
  real,public,allocatable :: rflexopt(:),eirres(:,:,:),eirsdtor(:),eirsdvol(:),eirstrata(:,:),&
       eirtrans(:,:),eirstrdat(:,:,:),eiriontime(:,:),eirscale(:),eirfluxt(:),&
       eirpho1(:,:),eirpho2(:,:),eirtri(:,:)
  real,public :: eirtarg,eirtemp1,eirtemp2,eirzaa,eirtorfrac
  real,public,allocatable :: eirpgdat(:,:),eirasdat(:,:),eirspdat(:,:),eirpuff(:,:)
  !
  ! neutral wall stuff:
  !
  integer,public :: eirtrc1,eirtrc2,eird2fu,eirspdatmode
  !
  ! grid:
  !
  real,public,allocatable :: wallpt2(:,:)
  
  integer,public,allocatable :: virloc(:,:)
  integer,public :: in14,out23
  
  
  parameter (in14 = 1, out23 = 2)
  
  real,public,allocatable :: rho(:,:)
  integer,public :: grd_refine
  
  !
  ! cmod experimental data:
  !
  real,public :: grd_minpl,grd_range,grd_shift,grd_thresh
  
  integer,public :: numprb,maxprb,numhal,maxhal,bside_st,btop_st,fside,kbot_st,ofmp,&
       ifmp,fsp1,d3d_rcp
  
  parameter (numprb   = 3, maxprb   = 99,numhal   = 4, maxhal   = 64,fsp1     = 1, ifmp     = 2 ,&
        ofmp  = 3 , d3d_rcp = 4,bside_st = 1, btop_st  = 2 , fside = 3 , kbot_st = 4)
  integer,public,allocatable :: prb_num(:)
  real,public,allocatable :: prb_te(:,:),prb_ti(:,:),prb_ne(:,:),prb_rho(:,:),prb_r(:,:),&
       prb_z(:,:)
  
  character*72,public :: prb_tag(numprb)
  real,public,allocatable :: dip_te(:,:),dip_ti(:,:),dip_ne(:,:),dip_s(:,:),dip_v(:,:)
  
  integer,public,allocatable :: dip_ik(:,:)
  !
  !
  !
  real,public,allocatable :: prp_te(:,:),prp_ti(:,:),prp_ne(:,:)
  integer,public,allocatable :: hal_num(:)
  real,public,allocatable :: hal_ang(:,:),hal_val(:,:)
  
  character*72,public :: hal_tag(numhal)
  integer,public :: tarsource,haldata,shot,slice_time,prb_align,osmnppv,tarshiftopt
  !
  !     relaxation:
  !
  !
  !     rel_viter - number of iterations at each step for pressure comparison
  !                 option
  !
  real,public :: prb_shift
  real,public,allocatable :: tarshift(:)
  
  !     .        osm_dp1(2,maxnrs),osm_dp2(2,maxnrs),osm_dp3(2,maxnrs)
  real,public :: rel_frac,rel_pace,rel_bound1,rel_bound2,rel_tol,osm_testep
  real,public,allocatable :: rel_data(:,:),rel_qemul(:,:),osm_peimul(:,:),osm_dp1(:,:),&
       osm_dp2(:,:),osm_dp3(:,:),osm_dp4(:,:),osm_dp5(:,:),osm_dp6(:,:),osm_ionfnt(:,:),&
       osmppv(:,:)
  
  real,public,allocatable :: rel_dirtg(:),rel_deltati(:),rel_hte(:),rel_hti(:),rel_hne(:),&
       rel_hproe(:,:),rel_hproi(:,:),rel_symfr(:),rel_prbfr(:,:),rel_mfsp(:,:),relexpdat(:,:,:)
  
  !     .        osm_symopt,osm_powopt,osm_preopt,osm_recopt
  integer,public :: rel_opt,rel_count,rel_niter,rel_nstep,rel_iter,rel_step,nlpdati2,&
       nlpdato2,rel_ndata,relmode,relreset,osm_symopt,osm_powopt,osm_preopt,osm_recopt,&
       osm_nfnt,osmmock
  integer,public,allocatable :: rel_viter(:),osm_sympt(:),osm_code(:,:),osm_model(:,:),&
       osmikshift(:,:)
  
  real*8,public :: relseed
  
  real,public,allocatable :: lpdati2(:,:),lpdato2(:,:)
  
  real,public,allocatable :: pinion2(:,:),pinqe2(:,:),pinrec2(:,:),pinqi2(:,:),pinmol2(:,:),&
       pinena2(:,:),pinenz2(:,:),pinenm2(:,:),pinmp2(:,:),pinatom2(:,:),pinalpha2(:,:),&
       pinionz2(:,:),pinz02(:,:)
  
  real,public,allocatable :: pinline(:,:,:,:),pinline2(:,:,:,:),pindata(:,:,:),pindata2(:,:,:),&
       mulrec(:,:),mulion(:,:),mulrecl(:),mulqei(:,:),mulqer(:,:),pinbgk(:,:,:),&
       pinbgk2(:,:,:),pinior(:,:),pinior2(:,:),pinmpr(:,:),pinmpr2(:,:),pinqir(:,:),&
       pinqir2(:,:),pinqer(:,:),pinqer2(:,:),pinmoi(:,:),pinstrata(:,:,:,:),pinploss(:,:,:),&
       pinalgv(:,:,:),pinioncomp(:,:,:)
  
  integer,public :: indasd
  integer,public :: maxasc,maxasc3d,maxascdat
  
  parameter (maxasc=5000,maxasc3d=200,maxascdat=40000)
  ! common /asccom/  asc_cell,asc_link,asc_nvp,asc_rvp,asc_zvp,asc_ncell,asc_grid,asc_region,&
  !     asc_vol,asc_area,asc_nregion,asc_rstart,asc_rend,asccode,ascvertex,ascnvertex,&
  !     ascdata,ascncut,asc_link3d,asc_nvp3d,asc_zmin3d,asc_zmax3d,asc_xvp3d ,asc_yvp3d,&
  !     asc_zvp3d,asc_3dmode
  ! save /asccom/
  integer,public :: asc_ncell,asc_nregion,asccode,ascncut,asc_3dmode
  integer,public,allocatable :: asc_rstart(:),asc_rend(:),asc_cell(:),asc_region(:),&
       asc_link(:,:),asc_grid(:,:),asc_nvp(:),ascnvertex(:),asc_link3d(:,:),asc_nvp3d(:)
  
  real,public,allocatable :: asc_rvp(:,:),asc_zvp(:,:),asc_vol(:),ascvertex(:,:),asc_area(:),&
       ascdata(:,:),asc_zmin3d(:),asc_zmax3d(:),asc_xvp3d(:,:,:),asc_yvp3d(:,:,:),&
       asc_zvp3d(:,:,:)
  
  
  
  !
  !     grid adaptation:
  !
  real,public,allocatable :: pinasd(:,:,:,:)
  integer,public :: adp_opt,adp_region
  !
  ! osm options:
  !
  !
  ! osm_nflag   - flag to see if pressure is within allowed limits
  !
  real,public :: adp_upper,adp_lower
  
  
  real,public :: osm_range,s22_l1,s22_nr,s22_tr,s22_maxs,osmbulkn,osmbulkte,osmbulkti,&
       osmbulkv
  real,public,allocatable :: s21_datai(:,:),s21_datao(:,:),osm_sfrac(:),osms28(:,:),&
       osms28over(:,:)
  
  !.... slfrric is really silly here: it is a copy of ffric in solcommon, but solcommon
  !     cannot be included in caldensitypeak (at the moment) because te1, etc. is declared
  !     in both.
  integer,public :: s21_ndatai,s21_ndatao,osm_store,osm_watch1,osm_watch2,s22_ir,osm_probe,&
       osm_nflag,s22pfzion,osm_mode,s22_midnks,s21_mode,osm_matchs,osm_matchp,ringmodel,&
       osmns28,osmns28over,s28sym
  integer,public,allocatable :: ikbound(:,:),ikfluid(:,:),osmikp(:,:,:)
  
  
  real*8,public :: ringuppr,slffric
  integer,public :: swpow2,swpow3,swion2
  
  !...these 100's are mxspts in solparams... move later...
  parameter (swpow2 = 38,swpow3 = 39,swion2 = 40)
  
  
  !
  !     statistics:
  !
  real*8,public :: dumpei1,dumpae1,dumpai1
  real*8,public,allocatable :: intpmk(:)
  ! common /tempcom2/ err1,err2,tim1,tim2,serr2
  ! save /tempcom2/
  integer,public :: err1,err2,tim1,tim2
  
  
  real,public,allocatable :: serr2(:)
  ! common /imag/ simag1,simag2
  ! save /imag/
  
  real,public :: simag1,simag2
  ! common /errcom/ ierror,ikerror
  ! save /errcom/
  
  
  ! dave's...
  integer,public :: ierror
  integer,public,allocatable :: ikerror(:,:)
  ! common /slcomd/ ne_opt,uedge_bg
  !
  !       new optional input parameters
  !
  ! save /slcomd/
  
  integer,public :: ne_opt,uedge_bg
  ! common /slcom6/osmpar,osmmp  ,osmqe   ,osmqi, osmpot
  
  ! save /slcom6/
  
  
  !
  !     jdemod - moved osmpmk2 to start of common to avoid
  !              alignment issue on some compilers
  !
  real,public,allocatable :: osmpar(:,:),osmmp(:,:),osmqe(:,:),osmqi(:,:),osmpot(:,:)
  !...   sol27:
  ! common /slcom7/ osmpmk2,osmpei ,osmrad,osmmpi,osmmpe,osmnppv,osmppv,osmmp2 ,osmqe2  ,&
  !     osmpei2 ,osmrad2 ,osmmpi2 ,osmmpe2,osmcfp ,osmcfi  ,osmcfe  ,osmcfpflx,osmcve ,&
  !     osmcde  ,osmcvi ,osmcdi,osmion,osmrec   ,osmpmk
  !
  !     jdemod
  !    .  , osmpmk2
  !
  ! save /slcom7/
  
  
  !...   sol27:
  real,public,allocatable :: osmpei(:,:),osmpei2(:,:),osmmp2(:,:),osmqi2(:,:),osmcfp(:,:),&
       osmcfpflx(:,:,:),osmcfi(:,:),osmcfe(:,:),osmcve(:,:),osmcde(:,:),osmcvi(:,:),&
       osmcdi(:,:),osmrad(:,:),osmrad2(:,:),osmmpi(:,:),osmmpi2(:,:),osmmpe(:,:),osmmpe2(:,:),&
       osmqe2(:,:),osmpmk(:,:),osmion(:,:),osmrec(:,:)
  
  real*8,public,allocatable :: osmpmk2(:,:)
  logical,public :: thesis
  
  real,public :: tuprati,tuprat
  integer,public :: vacnseg,vacnpla
  integer,public,allocatable :: vacregion(:)
  
  real,public,allocatable :: vacseg(:,:),vacpla(:,:)
  ! common /intsource/ intion1,intrec1,intqi1,intqe1,intpei1,intcfp1,intcfe1,intcfi1,&
  !     intmp1,intqe2,intrad1,intpmk1
  ! save /intsource/
  
  
  real,public,allocatable :: intion1(:,:),intrec1(:,:),intqi1(:,:),intqe1(:,:),intpei1(:,:),&
       intcfp1(:,:),intcfe1(:,:),intcfi1(:,:),intmp1(:,:),intqe2(:,:),intrad1(:,:),&
       intpmk1(:,:)
  ! common /sol22com1/ ionsum,ionsumt,       recsumt,cfpsum,cfpsumt,cfesum,cfesumt,&
  !     cfisum,cfisumt,       qesumt ,qisumt ,peisum,peisumt,radsum,radsumt,pmksum,pmksumt,&
  !     tifrac
  
  
  !...  tags:
  real*8,public :: ionsum,ionsumt,recsumt,cfpsum,cfpsumt,cfesum,cfesumt,cfisum,cfisumt,&
       qesumt,qisumt,peisum,peisumt,radsum,radsumt,pmksum,pmksumt,tifrac
  
  !...  sol28 options:
  logical,public :: tagpinatom,tagmulrec
  ! common /sol28com/s28fit,s28probe,s28probepfz,s28te,s28ion   ,s28cfp     ,s28rec   ,&
  !     s28mom   ,s28cfm   ,s28ti   ,s28cfpnt   ,s28cfpdrft,s28tepfz,s28ionpfz,s28cfppfz  ,&
  !     s28recpfz,s28mompfz,s28cfmpfz,s28tipfz,s28cfppfznt,s28cfppfzdrft,s28superdet   ,&
  !     s28nemode   ,s28superdetpfz,s28nemodepfz,s28ionset,s28recset,s28momset,s28ionsetpfz,&
  !     s28recsetpfz,s28momsetpfz,s28momfr,s28mode,s28b34,s28b34det,s28tiratio,&
  !     s28ionfrac,s28recfrac,s28momte,s28momfrac,s28ionbound   ,s28ionboundpfz
  
  ! save /sol28com/
  integer,public :: s28fit,s28probe,s28probepfz,s28te,s28ion,s28cfp,s28rec,s28mom,s28cfm,&
       s28ti,s28cfpnt,s28cfpdrft,s28tepfz,s28ionpfz,s28cfppfz,s28recpfz,s28mompfz,&
       s28cfmpfz,s28tipfz,s28cfppfznt,s28cfppfzdrft,s28superdet,s28nemode,s28superdetpfz,&
       s28nemodepfz
  logical,public :: s28ionset,s28recset,s28momset,s28ionsetpfz,s28recsetpfz,s28momsetpfz
  
  
  !...  target data for interpolation:
  real,public :: s28momfr,s28mode,s28b34,s28b34det,s28tiratio,s28momte,s28ionbound,&
       s28ionboundpfz
  real,public,allocatable :: s28ionfrac(:,:),s28recfrac(:,:),s28momfrac(:,:)
  ! common /tarcom/ tarinter,tarninter,tarintermode
  
  ! save /tarcom/
  integer,public,allocatable :: tarninter(:)
  
  !...  grid related data:
  real,public,allocatable :: tarinter(:,:,:),tarintermode(:)
  ! common /grdcom/ grdmod,grdnmod,psindat,psidat,quasidn,ikmidplane,grdntreg,grdntseg,&
  !     grdtseg
  ! save /grdcom/
  integer,public :: grdnmod
  integer,public,allocatable :: psindat(:),ikmidplane(:,:),grdntreg(:),grdntseg(:,:),&
       grdtseg(:,:,:)
  logical,public :: quasidn
  
  real,public,allocatable :: grdmod(:,:),psidat(:,:)
  ! common /debugcom/ sldebug
  ! save /debugcom/
  
  integer,public :: sldebug
  ! common /gridcom/ connected
  ! save /gridcom/
  logical,public :: connected

  public :: allocate_mod_slcom,deallocate_mod_slcom

contains

  subroutine allocate_mod_slcom
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_slcom','ALLOCATE')

    call allocate_array(disindex,maxnrs,'disindex',ierr)
    call allocate_array(dislev,maxnrs,'dislev',ierr)
    call allocate_array(ikbrks,4,maxnrs,'ikbrks',ierr)
    call allocate_array(ikbrke,4,maxnrs,'ikbrke',ierr)
    call allocate_array(nbreak,maxnrs,'nbreak',ierr)
    call allocate_array(irorg2,maxnrs,'irorg2',ierr)
    call allocate_array(idring,maxnrs,'idring',ierr)
    call allocate_array(supflx,2,maxnrs,'supflx',ierr)
    call allocate_array(ringtype,maxnrs,'ringtype',ierr)
    call allocate_array(ikti2,maxnrs,'ikti2',ierr)
    call allocate_array(ikto2,maxnrs,'ikto2',ierr)
    call allocate_array(idtarg,maxnds,'idtarg',ierr)
    call allocate_array(ikbreak,maxnrs,'ikbreak',ierr)
    call allocate_array(virtag,0,maxnks+1,0,maxnrs,'virtag',ierr)
    call allocate_array(eirntracks,1000,'eirntracks',ierr)
    call allocate_array(eiraout,maxasd,10,'eiraout',ierr)
    call allocate_array(eirstrlis,maxstrdat,'eirstrlis',ierr)
    call allocate_array(eirnlimi,5000,'eirnlimi',ierr)
    call allocate_array(eirtriik,maxnks*maxnrs,'eirtriik',ierr)
    call allocate_array(eirtriir,maxnks*maxnrs,'eirtriir',ierr)
    call allocate_array(iflexopt,maxflex,'iflexopt',ierr)
    call allocate_array(rflexopt,maxflex,'rflexopt',ierr)
    call allocate_array(eirres,7,6,maxpiniter,'eirres',ierr)
    call allocate_array(eirsdtor,maxtor,'eirsdtor',ierr)
    call allocate_array(eirsdvol,maxtor,'eirsdvol',ierr)
    call allocate_array(eirstrata,maxstrata,10,'eirstrata',ierr)
    call allocate_array(eirtrans,maxtor,3,'eirtrans',ierr)
    call allocate_array(eirstrdat,maxstrdat,maxstrata,100,'eirstrdat',ierr)
    call allocate_array(eiriontime,maxiontime,20+maxbin*3,'eiriontime',ierr)
    call allocate_array(eirscale,100,'eirscale',ierr)
    call allocate_array(eirfluxt,maxstrata,'eirfluxt',ierr)
    call allocate_array(eirpho1,maxnks,maxnrs,'eirpho1',ierr)
    call allocate_array(eirpho2,maxnks,maxnrs,'eirpho2',ierr)
    call allocate_array(eirtri,maxnrs,20,'eirtri',ierr)
    call allocate_array(eirpgdat,maxnas,maxasd,'eirpgdat',ierr)
    call allocate_array(eirasdat,maxnas2,maxasd,'eirasdat',ierr)
    call allocate_array(eirspdat,maxnas3,maxasd,'eirspdat',ierr)
    call allocate_array(eirpuff,maxnas,maxpuff,'eirpuff',ierr)
    call allocate_array(wallpt2,maxpts,2,'wallpt2',ierr)
    call allocate_array(virloc,maxnrs,2,'virloc',ierr)
    call allocate_array(rho,maxnrs,3,'rho',ierr)
    call allocate_array(prb_num,numprb,'prb_num',ierr)
    call allocate_array(prb_te,maxprb,numprb,'prb_te',ierr)
    call allocate_array(prb_ti,maxprb,numprb,'prb_ti',ierr)
    call allocate_array(prb_ne,maxprb,numprb,'prb_ne',ierr)
    call allocate_array(prb_rho,maxprb,numprb,'prb_rho',ierr)
    call allocate_array(prb_r,maxprb,numprb,'prb_r',ierr)
    call allocate_array(prb_z,maxprb,numprb,'prb_z',ierr)
    call allocate_array(dip_te,maxnrs,numprb,'dip_te',ierr)
    call allocate_array(dip_ti,maxnrs,numprb,'dip_ti',ierr)
    call allocate_array(dip_ne,maxnrs,numprb,'dip_ne',ierr)
    call allocate_array(dip_s,maxnrs,numprb,'dip_s',ierr)
    call allocate_array(dip_v,maxnrs,numprb,'dip_v',ierr)
    call allocate_array(dip_ik,maxnrs,numprb,'dip_ik',ierr)
    call allocate_array(prp_te,maxnrs,numprb,'prp_te',ierr)
    call allocate_array(prp_ti,maxnrs,numprb,'prp_ti',ierr)
    call allocate_array(prp_ne,maxnrs,numprb,'prp_ne',ierr)
    call allocate_array(hal_num,numhal,'hal_num',ierr)
    call allocate_array(hal_ang,maxhal,numhal,'hal_ang',ierr)
    call allocate_array(hal_val,maxhal,numhal,'hal_val',ierr)
    call allocate_array(tarshift,2,'tarshift',ierr)
    call allocate_array(rel_data,maxins,15,'rel_data',ierr)
    call allocate_array(rel_qemul,2,maxnrs,'rel_qemul',ierr)
    call allocate_array(osm_peimul,2,maxnrs,'osm_peimul',ierr)
    call allocate_array(osm_dp1,2,maxnrs,'osm_dp1',ierr)
    call allocate_array(osm_dp2,2,maxnrs,'osm_dp2',ierr)
    call allocate_array(osm_dp3,2,maxnrs,'osm_dp3',ierr)
    call allocate_array(osm_dp4,maxnks,maxnrs,'osm_dp4',ierr)
    call allocate_array(osm_dp5,2,maxnrs,'osm_dp5',ierr)
    call allocate_array(osm_dp6,maxnks,maxnrs,'osm_dp6',ierr)
    call allocate_array(osm_ionfnt,maxfnt,3,'osm_ionfnt',ierr)
    call allocate_array(osmppv,maxnrs,6,'osmppv',ierr)
    call allocate_array(rel_dirtg,maxnrs,'rel_dirtg',ierr)
    call allocate_array(rel_deltati,maxnrs,'rel_deltati',ierr)
    call allocate_array(rel_hte,maxnrs,'rel_hte',ierr)
    call allocate_array(rel_hti,maxnrs,'rel_hti',ierr)
    call allocate_array(rel_hne,maxnrs,'rel_hne',ierr)
    call allocate_array(rel_hproe,3,maxnrs,'rel_hproe',ierr)
    call allocate_array(rel_hproi,3,maxnrs,'rel_hproi',ierr)
    call allocate_array(rel_symfr,maxnrs,'rel_symfr',ierr)
    call allocate_array(rel_prbfr,2,maxnrs,'rel_prbfr',ierr)
    call allocate_array(rel_mfsp,2,maxnrs,'rel_mfsp',ierr)
    call allocate_array(relexpdat,2,3,maxnrs,'relexpdat',ierr)
    call allocate_array(rel_viter,0,'rel_viter',maxstep,ierr)
    call allocate_array(osm_sympt,maxnrs,'osm_sympt',ierr)
    call allocate_array(osm_code,2,maxnrs,'osm_code',ierr)
    call allocate_array(osm_model,2,maxnrs,'osm_model',ierr)
    call allocate_array(osmikshift,2,maxnrs,'osmikshift',ierr)
    call allocate_array(lpdati2,maxins,9,'lpdati2',ierr)
    call allocate_array(lpdato2,maxins,9,'lpdato2',ierr)
    call allocate_array(pinion2,maxnks,maxnrs,'pinion2',ierr)
    call allocate_array(pinqe2,maxnks,maxnrs,'pinqe2',ierr)
    call allocate_array(pinrec2,maxnks,maxnrs,'pinrec2',ierr)
    call allocate_array(pinqi2,maxnks,maxnrs,'pinqi2',ierr)
    call allocate_array(pinmol2,maxnks,maxnrs,'pinmol2',ierr)
    call allocate_array(pinena2,maxnks,maxnrs,'pinena2',ierr)
    call allocate_array(pinenz2,maxnks,maxnrs,'pinenz2',ierr)
    call allocate_array(pinenm2,maxnks,maxnrs,'pinenm2',ierr)
    call allocate_array(pinmp2,maxnks,maxnrs,'pinmp2',ierr)
    call allocate_array(pinatom2,maxnks,maxnrs,'pinatom2',ierr)
    call allocate_array(pinalpha2,maxnks,maxnrs,'pinalpha2',ierr)
    call allocate_array(pinionz2,maxnks,maxnrs,'pinionz2',ierr)
    call allocate_array(pinz02,maxnks,maxnrs,'pinz02',ierr)
    call allocate_array(pinline,1,maxnks,1,maxnrs,1,6,1,3,'pinline',ierr)
    call allocate_array(pinline2,1,maxnks,1,maxnrs,1,6,1,3,'pinline2',ierr)
    call allocate_array(pindata,maxnks,maxnrs,maxdata,'pindata',ierr)
    call allocate_array(pindata2,maxnks,maxnrs,maxdata,'pindata2',ierr)
    call allocate_array(mulrec,maxnks,maxnrs,'mulrec',ierr)
    call allocate_array(mulion,maxnks,maxnrs,'mulion',ierr)
    call allocate_array(mulrecl,maxnrs,'mulrecl',ierr)
    call allocate_array(mulqei,maxnks,maxnrs,'mulqei',ierr)
    call allocate_array(mulqer,maxnks,maxnrs,'mulqer',ierr)
    call allocate_array(pinbgk,maxnks,maxnrs,maxbgk*maxtor,'pinbgk',ierr)
    call allocate_array(pinbgk2,maxnks,maxnrs,maxbgk,'pinbgk2',ierr)
    call allocate_array(pinior,maxnks,maxnrs,'pinior',ierr)
    call allocate_array(pinior2,maxnks,maxnrs,'pinior2',ierr)
    call allocate_array(pinmpr,maxnks,maxnrs,'pinmpr',ierr)
    call allocate_array(pinmpr2,maxnks,maxnrs,'pinmpr2',ierr)
    call allocate_array(pinqir,maxnks,maxnrs,'pinqir',ierr)
    call allocate_array(pinqir2,maxnks,maxnrs,'pinqir2',ierr)
    call allocate_array(pinqer,maxnks,maxnrs,'pinqer',ierr)
    call allocate_array(pinqer2,maxnks,maxnrs,'pinqer2',ierr)
    call allocate_array(pinmoi,maxnks,maxnrs,'pinmoi',ierr)
    call allocate_array(pinstrata,1,maxnks,1,maxnrs,1,3,1,maxstrata,'pinstrata',ierr)
    call allocate_array(pinploss,maxnks,maxnrs,nmomcha,'pinploss',ierr)
    call allocate_array(pinalgv,maxnks,maxnrs,maxtally,'pinalgv',ierr)
    call allocate_array(pinioncomp,maxnks,maxnrs,3,'pinioncomp',ierr)
    call allocate_array(asc_rstart,10,'asc_rstart',ierr)
    call allocate_array(asc_rend,10,'asc_rend',ierr)
    call allocate_array(asc_cell,maxasc,'asc_cell',ierr)
    call allocate_array(asc_region,maxasc,'asc_region',ierr)
    call allocate_array(asc_link,4,maxasc,'asc_link',ierr)
    call allocate_array(asc_grid,2,maxasc,'asc_grid',ierr)
    call allocate_array(asc_nvp,maxasc,'asc_nvp',ierr)
    call allocate_array(ascnvertex,maxasc,'ascnvertex',ierr)
    call allocate_array(asc_link3d,6,maxasc3d,'asc_link3d',ierr)
    call allocate_array(asc_nvp3d,maxasc3d,'asc_nvp3d',ierr)
    call allocate_array(asc_rvp,8,maxasc,'asc_rvp',ierr)
    call allocate_array(asc_zvp,8,maxasc,'asc_zvp',ierr)
    call allocate_array(asc_vol,maxascdat,'asc_vol',ierr)
    call allocate_array(ascvertex,80,maxasc,'ascvertex',ierr)
    call allocate_array(asc_area,maxascdat,'asc_area',ierr)
    call allocate_array(ascdata,maxascdat,5,'ascdata',ierr)
    call allocate_array(asc_zmin3d,maxasc3d,'asc_zmin3d',ierr)
    call allocate_array(asc_zmax3d,maxasc3d,'asc_zmax3d',ierr)
    call allocate_array(asc_xvp3d,6,8,maxasc3d,'asc_xvp3d',ierr)
    call allocate_array(asc_yvp3d,6,8,maxasc3d,'asc_yvp3d',ierr)
    call allocate_array(asc_zvp3d,6,8,maxasc3d,'asc_zvp3d',ierr)
    call allocate_array(pinasd,1,maxascdat,1,maxasd2,1,maxass,1,2,'pinasd',ierr)
    call allocate_array(s21_datai,maxnrs,20,'s21_datai',ierr)
    call allocate_array(s21_datao,maxnrs,20,'s21_datao',ierr)
    call allocate_array(osm_sfrac,maxnrs,'osm_sfrac',ierr)
    call allocate_array(osms28,maxnrs,20,'osms28',ierr)
    call allocate_array(osms28over,maxnrs,20,'osms28over',ierr)
    call allocate_array(ikbound,maxnrs,2,'ikbound',ierr)
    call allocate_array(ikfluid,2,maxnrs,'ikfluid',ierr)
    call allocate_array(osmikp,2,maxnrs,5,'osmikp',ierr)
    call allocate_array(intpmk,100,'intpmk',ierr)
    call allocate_array(serr2,2,'serr2',ierr)
    call allocate_array(ikerror,2,maxnrs,'ikerror',ierr)
    call allocate_array(osmpar,maxnks,maxnrs,'osmpar',ierr)
    call allocate_array(osmmp,maxnks,maxnrs,'osmmp',ierr)
    call allocate_array(osmqe,maxnks,maxnrs,'osmqe',ierr)
    call allocate_array(osmqi,maxnks,maxnrs,'osmqi',ierr)
    call allocate_array(osmpot,maxnks,maxnrs,'osmpot',ierr)
    call allocate_array(osmpei,maxnks,maxnrs,'osmpei',ierr)
    call allocate_array(osmpei2,maxnks,maxnrs,'osmpei2',ierr)
    call allocate_array(osmmp2,maxnks,maxnrs,'osmmp2',ierr)
    call allocate_array(osmqi2,maxnks,maxnrs,'osmqi2',ierr)
    call allocate_array(osmcfp,maxnks,maxnrs,'osmcfp',ierr)
    call allocate_array(osmcfpflx,maxnks,maxnrs,5,'osmcfpflx',ierr)
    call allocate_array(osmcfi,maxnks,maxnrs,'osmcfi',ierr)
    call allocate_array(osmcfe,maxnks,maxnrs,'osmcfe',ierr)
    call allocate_array(osmcve,maxnks,maxnrs,'osmcve',ierr)
    call allocate_array(osmcde,maxnks,maxnrs,'osmcde',ierr)
    call allocate_array(osmcvi,maxnks,maxnrs,'osmcvi',ierr)
    call allocate_array(osmcdi,maxnks,maxnrs,'osmcdi',ierr)
    call allocate_array(osmrad,maxnks,maxnrs,'osmrad',ierr)
    call allocate_array(osmrad2,maxnks,maxnrs,'osmrad2',ierr)
    call allocate_array(osmmpi,maxnks,maxnrs,'osmmpi',ierr)
    call allocate_array(osmmpi2,maxnks,maxnrs,'osmmpi2',ierr)
    call allocate_array(osmmpe,maxnks,maxnrs,'osmmpe',ierr)
    call allocate_array(osmmpe2,maxnks,maxnrs,'osmmpe2',ierr)
    call allocate_array(osmqe2,maxnks,maxnrs,'osmqe2',ierr)
    call allocate_array(osmpmk,0,maxnks,1,maxnrs,'osmpmk',ierr)
    call allocate_array(osmion,maxnks,maxnrs,'osmion',ierr)
    call allocate_array(osmrec,maxnks,maxnrs,'osmrec',ierr)
    call allocate_array(osmpmk2,0,maxnks,1,maxnrs,'osmpmk2',ierr)
    call allocate_array(vacregion,maxvacregion,'vacregion',ierr)
    call allocate_array(vacseg,maxnks,6,'vacseg',ierr)
    call allocate_array(vacpla,maxnks,10,'vacpla',ierr)
    call allocate_array(intion1,maxnks+1,maxnrs,'intion1',ierr)
    call allocate_array(intrec1,maxnks+1,maxnrs,'intrec1',ierr)
    call allocate_array(intqi1,maxnks+1,maxnrs,'intqi1',ierr)
    call allocate_array(intqe1,maxnks+1,maxnrs,'intqe1',ierr)
    call allocate_array(intpei1,maxnks+1,maxnrs,'intpei1',ierr)
    call allocate_array(intcfp1,maxnks+1,maxnrs,'intcfp1',ierr)
    call allocate_array(intcfe1,maxnks+1,maxnrs,'intcfe1',ierr)
    call allocate_array(intcfi1,maxnks+1,maxnrs,'intcfi1',ierr)
    call allocate_array(intmp1,maxnks+1,maxnrs,'intmp1',ierr)
    call allocate_array(intqe2,maxnks+1,maxnrs,'intqe2',ierr)
    call allocate_array(intrad1,maxnks+1,maxnrs,'intrad1',ierr)
    call allocate_array(intpmk1,maxnks+1,maxnrs,'intpmk1',ierr)
    call allocate_array(s28ionfrac,2,maxnrs,'s28ionfrac',ierr)
    call allocate_array(s28recfrac,2,maxnrs,'s28recfrac',ierr)
    call allocate_array(s28momfrac,2,maxnrs,'s28momfrac',ierr)
    call allocate_array(tarninter,2,'tarninter',ierr)
    call allocate_array(tarinter,maxnrs,10,2,'tarinter',ierr)
    call allocate_array(tarintermode,2,'tarintermode',ierr)
    call allocate_array(psindat,2,'psindat',ierr)
    call allocate_array(ikmidplane,maxnrs,2,'ikmidplane',ierr)
    call allocate_array(grdntreg,2,'grdntreg',ierr)
    call allocate_array(grdntseg,maxnrs,2,'grdntseg',ierr)
    call allocate_array(grdtseg,maxnrs,maxnrs,2,'grdtseg',ierr)
    call allocate_array(grdmod,1000,10,'grdmod',ierr)
    call allocate_array(psidat,maxnrs+1,4,'psidat',ierr)

  end subroutine allocate_mod_slcom


  subroutine deallocate_mod_slcom
    implicit none

    call pr_trace('mod_slcom','DEALLOCATE')

    if (allocated(disindex)) deallocate(disindex)
    if (allocated(dislev)) deallocate(dislev)
    if (allocated(ikbrks)) deallocate(ikbrks)
    if (allocated(ikbrke)) deallocate(ikbrke)
    if (allocated(nbreak)) deallocate(nbreak)
    if (allocated(irorg2)) deallocate(irorg2)
    if (allocated(idring)) deallocate(idring)
    if (allocated(supflx)) deallocate(supflx)
    if (allocated(ringtype)) deallocate(ringtype)
    if (allocated(ikti2)) deallocate(ikti2)
    if (allocated(ikto2)) deallocate(ikto2)
    if (allocated(idtarg)) deallocate(idtarg)
    if (allocated(ikbreak)) deallocate(ikbreak)
    if (allocated(virtag)) deallocate(virtag)
    if (allocated(eirntracks)) deallocate(eirntracks)
    if (allocated(eiraout)) deallocate(eiraout)
    if (allocated(eirstrlis)) deallocate(eirstrlis)
    if (allocated(eirnlimi)) deallocate(eirnlimi)
    if (allocated(eirtriik)) deallocate(eirtriik)
    if (allocated(eirtriir)) deallocate(eirtriir)
    if (allocated(iflexopt)) deallocate(iflexopt)
    if (allocated(rflexopt)) deallocate(rflexopt)
    if (allocated(eirres)) deallocate(eirres)
    if (allocated(eirsdtor)) deallocate(eirsdtor)
    if (allocated(eirsdvol)) deallocate(eirsdvol)
    if (allocated(eirstrata)) deallocate(eirstrata)
    if (allocated(eirtrans)) deallocate(eirtrans)
    if (allocated(eirstrdat)) deallocate(eirstrdat)
    if (allocated(eiriontime)) deallocate(eiriontime)
    if (allocated(eirscale)) deallocate(eirscale)
    if (allocated(eirfluxt)) deallocate(eirfluxt)
    if (allocated(eirpho1)) deallocate(eirpho1)
    if (allocated(eirpho2)) deallocate(eirpho2)
    if (allocated(eirtri)) deallocate(eirtri)
    if (allocated(eirpgdat)) deallocate(eirpgdat)
    if (allocated(eirasdat)) deallocate(eirasdat)
    if (allocated(eirspdat)) deallocate(eirspdat)
    if (allocated(eirpuff)) deallocate(eirpuff)
    if (allocated(wallpt2)) deallocate(wallpt2)
    if (allocated(virloc)) deallocate(virloc)
    if (allocated(rho)) deallocate(rho)
    if (allocated(prb_num)) deallocate(prb_num)
    if (allocated(prb_te)) deallocate(prb_te)
    if (allocated(prb_ti)) deallocate(prb_ti)
    if (allocated(prb_ne)) deallocate(prb_ne)
    if (allocated(prb_rho)) deallocate(prb_rho)
    if (allocated(prb_r)) deallocate(prb_r)
    if (allocated(prb_z)) deallocate(prb_z)
    if (allocated(dip_te)) deallocate(dip_te)
    if (allocated(dip_ti)) deallocate(dip_ti)
    if (allocated(dip_ne)) deallocate(dip_ne)
    if (allocated(dip_s)) deallocate(dip_s)
    if (allocated(dip_v)) deallocate(dip_v)
    if (allocated(dip_ik)) deallocate(dip_ik)
    if (allocated(prp_te)) deallocate(prp_te)
    if (allocated(prp_ti)) deallocate(prp_ti)
    if (allocated(prp_ne)) deallocate(prp_ne)
    if (allocated(hal_num)) deallocate(hal_num)
    if (allocated(hal_ang)) deallocate(hal_ang)
    if (allocated(hal_val)) deallocate(hal_val)
    if (allocated(tarshift)) deallocate(tarshift)
    if (allocated(rel_data)) deallocate(rel_data)
    if (allocated(rel_qemul)) deallocate(rel_qemul)
    if (allocated(osm_peimul)) deallocate(osm_peimul)
    if (allocated(osm_dp1)) deallocate(osm_dp1)
    if (allocated(osm_dp2)) deallocate(osm_dp2)
    if (allocated(osm_dp3)) deallocate(osm_dp3)
    if (allocated(osm_dp4)) deallocate(osm_dp4)
    if (allocated(osm_dp5)) deallocate(osm_dp5)
    if (allocated(osm_dp6)) deallocate(osm_dp6)
    if (allocated(osm_ionfnt)) deallocate(osm_ionfnt)
    if (allocated(osmppv)) deallocate(osmppv)
    if (allocated(rel_dirtg)) deallocate(rel_dirtg)
    if (allocated(rel_deltati)) deallocate(rel_deltati)
    if (allocated(rel_hte)) deallocate(rel_hte)
    if (allocated(rel_hti)) deallocate(rel_hti)
    if (allocated(rel_hne)) deallocate(rel_hne)
    if (allocated(rel_hproe)) deallocate(rel_hproe)
    if (allocated(rel_hproi)) deallocate(rel_hproi)
    if (allocated(rel_symfr)) deallocate(rel_symfr)
    if (allocated(rel_prbfr)) deallocate(rel_prbfr)
    if (allocated(rel_mfsp)) deallocate(rel_mfsp)
    if (allocated(relexpdat)) deallocate(relexpdat)
    if (allocated(rel_viter)) deallocate(rel_viter)
    if (allocated(osm_sympt)) deallocate(osm_sympt)
    if (allocated(osm_code)) deallocate(osm_code)
    if (allocated(osm_model)) deallocate(osm_model)
    if (allocated(osmikshift)) deallocate(osmikshift)
    if (allocated(lpdati2)) deallocate(lpdati2)
    if (allocated(lpdato2)) deallocate(lpdato2)
    if (allocated(pinion2)) deallocate(pinion2)
    if (allocated(pinqe2)) deallocate(pinqe2)
    if (allocated(pinrec2)) deallocate(pinrec2)
    if (allocated(pinqi2)) deallocate(pinqi2)
    if (allocated(pinmol2)) deallocate(pinmol2)
    if (allocated(pinena2)) deallocate(pinena2)
    if (allocated(pinenz2)) deallocate(pinenz2)
    if (allocated(pinenm2)) deallocate(pinenm2)
    if (allocated(pinmp2)) deallocate(pinmp2)
    if (allocated(pinatom2)) deallocate(pinatom2)
    if (allocated(pinalpha2)) deallocate(pinalpha2)
    if (allocated(pinionz2)) deallocate(pinionz2)
    if (allocated(pinz02)) deallocate(pinz02)
    if (allocated(pinline)) deallocate(pinline)
    if (allocated(pinline2)) deallocate(pinline2)
    if (allocated(pindata)) deallocate(pindata)
    if (allocated(pindata2)) deallocate(pindata2)
    if (allocated(mulrec)) deallocate(mulrec)
    if (allocated(mulion)) deallocate(mulion)
    if (allocated(mulrecl)) deallocate(mulrecl)
    if (allocated(mulqei)) deallocate(mulqei)
    if (allocated(mulqer)) deallocate(mulqer)
    if (allocated(pinbgk)) deallocate(pinbgk)
    if (allocated(pinbgk2)) deallocate(pinbgk2)
    if (allocated(pinior)) deallocate(pinior)
    if (allocated(pinior2)) deallocate(pinior2)
    if (allocated(pinmpr)) deallocate(pinmpr)
    if (allocated(pinmpr2)) deallocate(pinmpr2)
    if (allocated(pinqir)) deallocate(pinqir)
    if (allocated(pinqir2)) deallocate(pinqir2)
    if (allocated(pinqer)) deallocate(pinqer)
    if (allocated(pinqer2)) deallocate(pinqer2)
    if (allocated(pinmoi)) deallocate(pinmoi)
    if (allocated(pinstrata)) deallocate(pinstrata)
    if (allocated(pinploss)) deallocate(pinploss)
    if (allocated(pinalgv)) deallocate(pinalgv)
    if (allocated(pinioncomp)) deallocate(pinioncomp)
    if (allocated(asc_rstart)) deallocate(asc_rstart)
    if (allocated(asc_rend)) deallocate(asc_rend)
    if (allocated(asc_cell)) deallocate(asc_cell)
    if (allocated(asc_region)) deallocate(asc_region)
    if (allocated(asc_link)) deallocate(asc_link)
    if (allocated(asc_grid)) deallocate(asc_grid)
    if (allocated(asc_nvp)) deallocate(asc_nvp)
    if (allocated(ascnvertex)) deallocate(ascnvertex)
    if (allocated(asc_link3d)) deallocate(asc_link3d)
    if (allocated(asc_nvp3d)) deallocate(asc_nvp3d)
    if (allocated(asc_rvp)) deallocate(asc_rvp)
    if (allocated(asc_zvp)) deallocate(asc_zvp)
    if (allocated(asc_vol)) deallocate(asc_vol)
    if (allocated(ascvertex)) deallocate(ascvertex)
    if (allocated(asc_area)) deallocate(asc_area)
    if (allocated(ascdata)) deallocate(ascdata)
    if (allocated(asc_zmin3d)) deallocate(asc_zmin3d)
    if (allocated(asc_zmax3d)) deallocate(asc_zmax3d)
    if (allocated(asc_xvp3d)) deallocate(asc_xvp3d)
    if (allocated(asc_yvp3d)) deallocate(asc_yvp3d)
    if (allocated(asc_zvp3d)) deallocate(asc_zvp3d)
    if (allocated(pinasd)) deallocate(pinasd)
    if (allocated(s21_datai)) deallocate(s21_datai)
    if (allocated(s21_datao)) deallocate(s21_datao)
    if (allocated(osm_sfrac)) deallocate(osm_sfrac)
    if (allocated(osms28)) deallocate(osms28)
    if (allocated(osms28over)) deallocate(osms28over)
    if (allocated(ikbound)) deallocate(ikbound)
    if (allocated(ikfluid)) deallocate(ikfluid)
    if (allocated(osmikp)) deallocate(osmikp)
    if (allocated(intpmk)) deallocate(intpmk)
    if (allocated(serr2)) deallocate(serr2)
    if (allocated(ikerror)) deallocate(ikerror)
    if (allocated(osmpar)) deallocate(osmpar)
    if (allocated(osmmp)) deallocate(osmmp)
    if (allocated(osmqe)) deallocate(osmqe)
    if (allocated(osmqi)) deallocate(osmqi)
    if (allocated(osmpot)) deallocate(osmpot)
    if (allocated(osmpei)) deallocate(osmpei)
    if (allocated(osmpei2)) deallocate(osmpei2)
    if (allocated(osmmp2)) deallocate(osmmp2)
    if (allocated(osmqi2)) deallocate(osmqi2)
    if (allocated(osmcfp)) deallocate(osmcfp)
    if (allocated(osmcfpflx)) deallocate(osmcfpflx)
    if (allocated(osmcfi)) deallocate(osmcfi)
    if (allocated(osmcfe)) deallocate(osmcfe)
    if (allocated(osmcve)) deallocate(osmcve)
    if (allocated(osmcde)) deallocate(osmcde)
    if (allocated(osmcvi)) deallocate(osmcvi)
    if (allocated(osmcdi)) deallocate(osmcdi)
    if (allocated(osmrad)) deallocate(osmrad)
    if (allocated(osmrad2)) deallocate(osmrad2)
    if (allocated(osmmpi)) deallocate(osmmpi)
    if (allocated(osmmpi2)) deallocate(osmmpi2)
    if (allocated(osmmpe)) deallocate(osmmpe)
    if (allocated(osmmpe2)) deallocate(osmmpe2)
    if (allocated(osmqe2)) deallocate(osmqe2)
    if (allocated(osmpmk)) deallocate(osmpmk)
    if (allocated(osmion)) deallocate(osmion)
    if (allocated(osmrec)) deallocate(osmrec)
    if (allocated(osmpmk2)) deallocate(osmpmk2)
    if (allocated(vacregion)) deallocate(vacregion)
    if (allocated(vacseg)) deallocate(vacseg)
    if (allocated(vacpla)) deallocate(vacpla)
    if (allocated(intion1)) deallocate(intion1)
    if (allocated(intrec1)) deallocate(intrec1)
    if (allocated(intqi1)) deallocate(intqi1)
    if (allocated(intqe1)) deallocate(intqe1)
    if (allocated(intpei1)) deallocate(intpei1)
    if (allocated(intcfp1)) deallocate(intcfp1)
    if (allocated(intcfe1)) deallocate(intcfe1)
    if (allocated(intcfi1)) deallocate(intcfi1)
    if (allocated(intmp1)) deallocate(intmp1)
    if (allocated(intqe2)) deallocate(intqe2)
    if (allocated(intrad1)) deallocate(intrad1)
    if (allocated(intpmk1)) deallocate(intpmk1)
    if (allocated(s28ionfrac)) deallocate(s28ionfrac)
    if (allocated(s28recfrac)) deallocate(s28recfrac)
    if (allocated(s28momfrac)) deallocate(s28momfrac)
    if (allocated(tarninter)) deallocate(tarninter)
    if (allocated(tarinter)) deallocate(tarinter)
    if (allocated(tarintermode)) deallocate(tarintermode)
    if (allocated(psindat)) deallocate(psindat)
    if (allocated(ikmidplane)) deallocate(ikmidplane)
    if (allocated(grdntreg)) deallocate(grdntreg)
    if (allocated(grdntseg)) deallocate(grdntseg)
    if (allocated(grdtseg)) deallocate(grdtseg)
    if (allocated(grdmod)) deallocate(grdmod)
    if (allocated(psidat)) deallocate(psidat)

  end subroutine deallocate_mod_slcom

end module mod_slcom