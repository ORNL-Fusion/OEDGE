module mod_comtor
  use debug_options
  use mod_walls_com
  implicit none

  !
  !     -*-fortran-*-
  !     >  cpdrft,cdrftv,ceimp,fluxpts,fluxropt,fluxinfo,srootopt,csopt2,
  !     >  neut2d_raw,sol22_power_ratio,redefopt,cdrftv_start,cdrftv_end,
  !
  ! common /comtor/ cbphi ,cxtarg,chalfa,crmb  ,cteb0 ,ciopta,cneuta,cnb0  ,cnbout,&
  !     cnbin ,cnbt  ,ctwoa ,ctebin,ctebou,cioptb,cneutb,cstepl,cdperp,ctibou,crmi  ,cxsc  ,&
  !     debugl,ckmax ,cioptc,cneutc,cysc  ,ctem1 ,ctem2 ,cebd  ,clarmr,cemaxf,ctib0 ,&
  !     cioptd,cneutd,cizb  ,cion  ,cizsc ,ctibin,debugn,ctresh,cizeff,ciopte,cneute,cnhc  ,&
  !     cnho  ,clamhx,clamhy,cxnear,cynear,cvcx  ,cioptf,cneutf,cprint,csolef,cein  ,&
  !     cvhin ,cstepn,cizset,czenh ,cioptg,ciopti,ceout ,cvhout,ciseed,cdifop,cbombz,cymfs ,&
  !     cbombf,ciopth,cirf  ,csef  ,qtim  ,fsrate,rizb  ,ckmin ,crect ,csolvf,csnorm,&
  !     cneutd2,cstop ,cksc  ,irspec,ctibt ,ctebt ,cneutg,cfl   ,cfs   ,cioptj,cperip,&
  !     cfrm  ,cfrmin,cfrmax,ckin  ,ckout ,cioptk,cioptl,cfebl1,cfebt ,cfibl1,cfibt ,cioptm,&
  !     cioptn,cfebl2,cfibl2,cfeb2 ,cfib2 ,absfac, absfac_neut,absfac_ion,czo   ,chzo  ,&
  !     cniwa ,tloss,ctebp ,ctibp ,cpa ,ck0,cnebp,lpdati,lpdato,nlpdati,nlpdato,cionr,&
  !     ofield,ceflen,ofield_targ,nrfopt,csolls,csollr,csolpr, ck0i,ciopto,csopt,cpopt,&
  !     czd,injir,injf1,injf2,cneuth,cneuti,ceimp,fluxpts,fluxropt,fluxinfo,srootopt,csopt2,&
  !     cnin,ctemav,lpdatsw,ccoreopt,coredat,ncoredat,rzopt,pdopt,wtsource,wtdep,fgradopt,&
  !     fgradfact,ceffact,sdperpref,sdperppp,const_yield,cfolrec,readaux,cvamult,&
  !     cvrmult,cflatopt,bgplasopt,nbgplas,mtcopt,kelighi,kelighg,ircore,vernum,revnum,&
  !     cneutvel,neut2d_opt,neut2d_vaopt,ngradopt,neut2d_raw,sol22_power_ratio,redefopt,&
  !     cbomb_frac,init_pos_opt,override_bg_velocity_opt,sonnet_grid_sub_type,ext_flx_data_src
  !
  !     >  sol23_izlen,sol23_izlam,sol23_izoffset,sol23_momlen,
  !     >  sol23_intopt,sol23_bndcond,sol23_seed,
  !
  ! common /comtor2/ cpinopt,citersol,csecsol,cpincom,actpin,ctargopt, platco,wallco,&
  !     nplat,nwall,csollt,cfiz,ctrap,csolfr,cgeoopt,wallpol,nitersol,cvolopt,cdiffopt,&
  !     walltemp,nwltemp,trappol,virtgrid,wallswch,injnum,injprob,injrind,injkind,nvaopt,&
  !     cmaxgens,cirhr,csputopt,extra_sputter_angle,xygrid,checkleak,solte,solti,solne,&
  !     solvel,solcor,solpr,solprh,solpcx,solpei,solph,cgridopt,wallco2,nwall2,kprat,kpsiz,&
  !     kprat2,cmiropt,ctimmax,cleq,irstold,ikstold,cleaks,cleakn,terat,nrat,qrat,l1rat,&
  !     l2rat,lvrat,nabsfac,cvbl1,cvbm1,cvbl2,cvbm2,cleakt,cleaksn,northopt,cmachno,&
  !     ctargt,cwallt,cutring,cutpt1,cutpt2,maxrings,maxkpts,cstgrad,ctestsol,cerr,cserr,&
  !     teupstream,tiupstream,cneur,cleakpos,nleakcore,piniter,launchdat,ionizdat,totpintim,&
  !     cellvals,cwallp,cteboup,ctiboup,cnboup,dpmethod,dpsuml,dpouter,dpconv,dpfluxopt,&
  !     dpavopt,dprcopt,dpsmooth,dpnav,dparea,dp_pinqe_mult,fluxes,dpploss,dporth,cosalph,&
  !     sinalph,dppei,dprec,cdeferr,cdefserr,ctes1,ctef1,ctes2,ctef2,ctis1,ctif1,&
  !     ctis2,ctif2,treccut,cnes1,cnef1,cnes2,cnef2,cvbs1,cvbf1,cvbs2,cvbf2,crecopt,kpress,&
  !     kprad,cchemopt,cselfs,s21refsw,te_mult_i,ti_mult_i,n_mult_i,te_mult_o,ti_mult_o,&
  !     n_mult_o,nbupstream,terati,nrati,qrati,l1rati,l2rati,lvrati,vbmult,vbmulti,crploc,&
  !     pincoreiz,pinsoliz,pinppiz,corefv,coreft,corefv2,coreft2,cfdopt,cdperpt,fixtgrad,&
  !     targsrc,targleak,dpxpratio,tirat,tirati,ctegcut,ctigcut,czploc,rlocnum,zlocnum,&
  !     nalph,nalphi,crdivbg,ns21i,s21parmi,ns21o,s21parmo,aux_ns21i,aux_s21parmi,aux_ns21o,&
  !     aux_s21parmo,aux_vel21,tmachine_opt,write_tran,netcdf_opt,refdist,cdeferropt,&
  !     nrat_used,debug_neutv,debug_neutv_einmax,debug_neutv_nbins
  !
  ! save /comtor/
  ! common /shotinfo/ divshotid
  !
  ! save /shotinfo/
  !     >  cdrftv,ceimp,fluxinfo(maxins,4),platco(maxnrs,5),
  !     >  walltemp(maxpts,3),cdrftv_start,cdrftv_end,s21parmi(maxnrs,10),
  !
  !     >  sol23_izlen,sol23_izlam,sol23_izoffset,
  !     >  sol23_momlen,
  !
  !
  real,public :: cbphi,chalfa,crmb,cteb0,cnb0,cnbout,cnbin,ctwoa,csolvf,ctebin,ctebou,&
       cdperp,crmi,cxsc,cysc,ctem1,ctem2,cebd,cemaxf,ctresh,cnhc,cnho,clamhx,clamhy,&
       cxnear,cynear,cein,cvhin,ceout,cvhout,cstepn,czenh,cirf,csef,csolef,cvcx,clarmr,&
       ctib0,cstepl,qtim,fsrate,ctibt,ctibin,cfl,cfs,cnbt,ctebt,cxtarg,ckmax,rizb,ckmin,&
       cksc,cperip,ctibou,cfrm,cfrmin,cfrmax,ckin,ckout,cfebl1,cfebt,cfibl1,cfibt,csnorm,&
       cfebl2,cfibl2,cfeb2,cfib2,absfac,absfac_neut,absfac_ion,chzo,tloss,ctebp,ctibp,&
       cpa,ck0,ck0i,cnebp,csolls,csollr,csolpr,czd,injf1,injf2,ctimmax,ctemav,ceimp,ctes1,&
       ctef1,ctes2,ctef2,ctis1,ctif1,ctis2,ctif2,treccut,cnes1,cnef1,cnes2,cnef2,cvbs1,&
       cvbf1,cvbs2,cvbf2,fgradfact,ceflen,ceffact,cdperpt,sdperpref,sdperppp,aux_vel21,&
       const_yield,cvamult,cvrmult
  real,public,allocatable :: cymfs(:,:),lpdati(:,:),lpdato(:,:),wallco2(:,:),fluxinfo(:,:),&
       platco(:,:),coredat(:,:),kpress(:,:,:),kprad(:,:),walltemp(:,:),s21parmi(:,:),&
       s21parmo(:,:),aux_s21parmi(:,:),aux_s21parmo(:,:),cdeferropt(:,:)
  ! common /pinch_data/ pdf_norm_val,pinchopt,pinch_loc_opt,pinch_npdf,npdf_data,pinch_pdf,&
  !     cvpinch,vpinch,pinch_pdf_data,pinch_correlation_time
  
  ! save /pinch_data/
  real,public :: pdf_norm_val
  integer,public :: pinchopt,pinch_npdf,pinch_loc_opt,npdf_data
  real,public :: cvpinch,vpinch
  real,public,allocatable :: pinch_pdf(:,:)
  real,public,allocatable :: pinch_pdf_data(:,:)
  
  !
  real,public :: pinch_correlation_time
  ! common /debug_pinch_data/ d_pinch_v
  
  ! save /debug_pinch_data/
  integer,public :: max_d_pinch_v
  real,public :: d_pinch_vel
  parameter (max_d_pinch_v=100,d_pinch_vel=10.0)
  !
  real*8,public,allocatable :: d_pinch_v(:)
  ! common  /sol13_options/ sol13_padd,sol13_pdist
  real,public :: sol13_padd,sol13_pdist
  !
  
  !
  real,public :: csollt,cfiz,csolfr,cleakt,ctestsol,cvbl1,cvbm1,cvbl2,cvbm2,cnin,terat,&
       nrat,qrat,l1rat,l2rat,lvrat,nabsfac,cstgrad,ctargt,cwallt,cwallp,totpintim,cteboup,&
       ctiboup,cnboup,dppei,dprec,dp_pinqe_mult,te_mult_i,ti_mult_i,n_mult_i,te_mult_o,&
       ti_mult_o,n_mult_o,crploc,terati,nrati,qrati,l1rati,l2rati,lvrati,vbmult,&
       vbmulti,pincoreiz,pinsoliz,extra_sputter_angle,pinppiz,corefv,coreft,corefv2,coreft2,&
       dpxpratio,tirat,tirati,ctegcut,ctigcut,czploc,kelighi,kelighg,cbomb_frac,nalph,&
       nalphi,debug_neutv_einmax
  real,public,allocatable :: wallco(:,:),injprob(:),cleq(:,:),cleaks(:),cleakn(:,:),&
       cmachno(:,:),solte(:),solti(:),solne(:),solvel(:),solcor(:),solpei(:),solpr(:),&
       solprh(:),solpcx(:),solph(:),tiupstream(:,:),teupstream(:,:),cleakpos(:,:),kprat(:,:),&
       kpsiz(:,:),kprat2(:,:,:),launchdat(:,:),ionizdat(:,:,:,:,:),cellvals(:,:,:),&
       cserr(:,:),fluxes(:,:,:),cosalph(:,:),sinalph(:,:),cdefserr(:,:),nbupstream(:,:),&
       wtsource(:,:,:,:),wtdep(:,:,:),targsrc(:,:),targleak(:,:),refdist(:),bgplasopt(:,:),&
       neut2d_raw(:,:),sol22_power_ratio(:,:,:),nrat_used(:,:)
  !     >  cpdrft,fluxpts,fluxropt,srootopt,cpinopt,citersol,csecsol,
  !
  !     >  sol23_intopt,sol23_bndcond,sol23_seed,
  !
  !
  integer,public :: ciopta,cioptb,cioptc,cioptd,ciopte,cioptf,cioptg,ciopth,ciopti,&
       cneuta,cneutb,cneutc,cneutd,cneute,cneutf,cstop,cneutg,cizb,cion,cizsc,cizeff,cizset,&
       ciseed,cdifop,cbombf,cbombz,cprint,cioptj,cioptk,cioptl,cioptm,cioptn,irspec,&
       crect,czo,cniwa,cionr,ofield,ofield_targ,nrfopt,ciopto,csopt,cpopt,lpdatsw,nlpdati,&
       nlpdato,injir,cneuth,cneuti,fluxpts,fluxropt,srootopt,cpinopt,citersol,csecsol,&
       csopt2,ctargopt,nplat,nwall,ctrap,cgeoopt,trappol,nleakcore,cvolopt,cdiffopt,&
       injnum,nwall2,nvaopt,cmaxgens,cirhr,cgridopt,dpmethod,dpsuml,dpouter,dpconv,csputopt,&
       cmiropt,nitersol,ikstold,irstold,xygrid,cleaksn,northopt,cneutd2,cutring,cutpt1,&
       cutpt2,maxrings,maxkpts,cneur,dpfluxopt,dpavopt,dprcopt,dpsmooth,dpnav,dparea,&
       dpploss,dporth,ccoreopt,ncoredat,cselfs,cchemopt,s21refsw,crecopt,nwltemp,rzopt,&
       pdopt,cfdopt,fgradopt,fixtgrad,rlocnum,zlocnum,ns21i,ns21o,aux_ns21i,aux_ns21o,&
       init_pos_opt,tmachine_opt,write_tran,netcdf_opt,cfolrec,readaux,cflatopt,nbgplas,&
       mtcopt,ircore,vernum,revnum,cneutvel,crdivbg,neut2d_opt,neut2d_vaopt,ngradopt,&
       redefopt,debug_neutv,debug_neutv_nbins,override_bg_velocity_opt,sonnet_grid_sub_type,&
       ext_flx_data_src
  integer,public,allocatable :: wallpol(:),injrind(:),injkind(:),cerr(:,:),cdeferr(:,:)
  !
  logical,public :: debugn,debugl,virtgrid,wallswch,checkleak,piniter
  character*80,public :: cpincom,actpin
  ! slmod begin
  !...  variables for line impurity injection (ciopte=9):
  character*10,public :: divshotid
  ! common /impcom01/ cxsca,cysca,cxscb,cyscb
  ! save /impcom01/
  
  ! slmod end
  !
  !     external plasma overlay option common block
  !
  real,public :: cxsca,cysca,cxscb,cyscb
  ! common /ext_plasma/ external_plasma_overlay,external_plasma_file
  integer,public :: external_plasma_overlay
  !
  !     include the wall common block for backward compatibility for now
  !     need to change later.
  !
  character*1024,public :: external_plasma_file
  
  
  
  !
  ! include 'walls_com'

  public :: allocate_mod_comtor,deallocate_mod_comtor

contains

  subroutine allocate_mod_comtor
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comtor','ALLOCATE')

    call allocate_array(cymfs,maxpts+1,8,'cymfs',ierr)
    call allocate_array(lpdati,maxins,4,'lpdati',ierr)
    call allocate_array(lpdato,maxins,4,'lpdato',ierr)
    call allocate_array(wallco2,maxpts,2,'wallco2',ierr)
    call allocate_array(fluxinfo,maxins,4,'fluxinfo',ierr)
    call allocate_array(platco,maxnrs,5,'platco',ierr)
    call allocate_array(coredat,maxins,5,'coredat',ierr)
    call allocate_array(kpress,maxnks,maxnrs,2,'kpress',ierr)
    call allocate_array(kprad,maxnks,maxnrs,'kprad',ierr)
    call allocate_array(walltemp,maxpts,3,'walltemp',ierr)
    call allocate_array(s21parmi,maxnrs,10,'s21parmi',ierr)
    call allocate_array(s21parmo,maxnrs,10,'s21parmo',ierr)
    call allocate_array(aux_s21parmi,maxnrs,9,'aux_s21parmi',ierr)
    call allocate_array(aux_s21parmo,maxnrs,9,'aux_s21parmo',ierr)
    call allocate_array(cdeferropt,maxnrs,2,'cdeferropt',ierr)
    call allocate_array(pinch_pdf,maxpts,2,'pinch_pdf',ierr)
    call allocate_array(pinch_pdf_data,maxpts,3,'pinch_pdf_data',ierr)
    call allocate_array(d_pinch_v,-max_d_pinch_v,'d_pinch_v',max_d_pinch_v,ierr)
    call allocate_array(wallco,maxpts,2,'wallco',ierr)
    call allocate_array(injprob,maxnks*maxnrs,'injprob',ierr)
    call allocate_array(cleq,maxnrs,2,'cleq',ierr)
    call allocate_array(cleaks,maxpts,'cleaks',ierr)
    call allocate_array(cleakn,maxpts,maxizs+1,'cleakn',ierr)
    call allocate_array(cmachno,maxnrs,2,'cmachno',ierr)
    call allocate_array(solte,0,'solte',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solti,0,'solti',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solne,0,'solne',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solvel,0,'solvel',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solcor,0,'solcor',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solpei,0,'solpei',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solpr,0,'solpr',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solprh,0,'solprh',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solpcx,0,'solpcx',maxnks*msolpt+msolpt,ierr)
    call allocate_array(solph,0,'solph',maxnks*msolpt+msolpt,ierr)
    call allocate_array(tiupstream,maxnrs,2,'tiupstream',ierr)
    call allocate_array(teupstream,maxnrs,2,'teupstream',ierr)
    call allocate_array(cleakpos,maximp,2,'cleakpos',ierr)
    call allocate_array(kprat,maxnrs,2,'kprat',ierr)
    call allocate_array(kpsiz,maxnks,maxnrs,'kpsiz',ierr)
    call allocate_array(kprat2,maxnks,maxnrs,2,'kprat2',ierr)
    call allocate_array(launchdat,maximp,5,'launchdat',ierr)
    call allocate_array(ionizdat,1,2,1,2,1,2,1,2,1,5,'ionizdat',ierr)
    call allocate_array(cellvals,maxnrs,4,2,'cellvals',ierr)
    call allocate_array(cserr,maxnrs,2,'cserr',ierr)
    call allocate_array(fluxes,maxnks,maxnrs,16,'fluxes',ierr)
    call allocate_array(cosalph,maxnks,maxnrs,'cosalph',ierr)
    call allocate_array(sinalph,maxnks,maxnrs,'sinalph',ierr)
    call allocate_array(cdefserr,maxnrs,2,'cdefserr',ierr)
    call allocate_array(nbupstream,maxnrs,2,'nbupstream',ierr)
    call allocate_array(wtsource,1,maxpts,1,maxnrs,1,4,1,6,'wtsource',ierr)
    call allocate_array(wtdep,maxpts,maxpts+1,3,'wtdep',ierr)
    call allocate_array(targsrc,3,4,'targsrc',ierr)
    call allocate_array(targleak,3,4,'targleak',ierr)
    call allocate_array(refdist,maxnrs,'refdist',ierr)
    call allocate_array(bgplasopt,2*maxnrs,12,'bgplasopt',ierr)
    call allocate_array(neut2d_raw,maxnks,maxnrs,'neut2d_raw',ierr)
    call allocate_array(sol22_power_ratio,maxnrs,2,3,'sol22_power_ratio',ierr)
    call allocate_array(nrat_used,maxnrs,2,'nrat_used',ierr)
    call allocate_array(wallpol,maxpts,'wallpol',ierr)
    call allocate_array(injrind,maxnks*maxnrs,'injrind',ierr)
    call allocate_array(injkind,maxnks*maxnrs,'injkind',ierr)
    call allocate_array(cerr,maxnrs,2,'cerr',ierr)
    call allocate_array(cdeferr,maxnrs,2,'cdeferr',ierr)

  end subroutine allocate_mod_comtor


  subroutine deallocate_mod_comtor
    implicit none

    call pr_trace('mod_comtor','DEALLOCATE')

    if (allocated(cymfs)) deallocate(cymfs)
    if (allocated(lpdati)) deallocate(lpdati)
    if (allocated(lpdato)) deallocate(lpdato)
    if (allocated(wallco2)) deallocate(wallco2)
    if (allocated(fluxinfo)) deallocate(fluxinfo)
    if (allocated(platco)) deallocate(platco)
    if (allocated(coredat)) deallocate(coredat)
    if (allocated(kpress)) deallocate(kpress)
    if (allocated(kprad)) deallocate(kprad)
    if (allocated(walltemp)) deallocate(walltemp)
    if (allocated(s21parmi)) deallocate(s21parmi)
    if (allocated(s21parmo)) deallocate(s21parmo)
    if (allocated(aux_s21parmi)) deallocate(aux_s21parmi)
    if (allocated(aux_s21parmo)) deallocate(aux_s21parmo)
    if (allocated(cdeferropt)) deallocate(cdeferropt)
    if (allocated(pinch_pdf)) deallocate(pinch_pdf)
    if (allocated(pinch_pdf_data)) deallocate(pinch_pdf_data)
    if (allocated(d_pinch_v)) deallocate(d_pinch_v)
    if (allocated(wallco)) deallocate(wallco)
    if (allocated(injprob)) deallocate(injprob)
    if (allocated(cleq)) deallocate(cleq)
    if (allocated(cleaks)) deallocate(cleaks)
    if (allocated(cleakn)) deallocate(cleakn)
    if (allocated(cmachno)) deallocate(cmachno)
    if (allocated(solte)) deallocate(solte)
    if (allocated(solti)) deallocate(solti)
    if (allocated(solne)) deallocate(solne)
    if (allocated(solvel)) deallocate(solvel)
    if (allocated(solcor)) deallocate(solcor)
    if (allocated(solpei)) deallocate(solpei)
    if (allocated(solpr)) deallocate(solpr)
    if (allocated(solprh)) deallocate(solprh)
    if (allocated(solpcx)) deallocate(solpcx)
    if (allocated(solph)) deallocate(solph)
    if (allocated(tiupstream)) deallocate(tiupstream)
    if (allocated(teupstream)) deallocate(teupstream)
    if (allocated(cleakpos)) deallocate(cleakpos)
    if (allocated(kprat)) deallocate(kprat)
    if (allocated(kpsiz)) deallocate(kpsiz)
    if (allocated(kprat2)) deallocate(kprat2)
    if (allocated(launchdat)) deallocate(launchdat)
    if (allocated(ionizdat)) deallocate(ionizdat)
    if (allocated(cellvals)) deallocate(cellvals)
    if (allocated(cserr)) deallocate(cserr)
    if (allocated(fluxes)) deallocate(fluxes)
    if (allocated(cosalph)) deallocate(cosalph)
    if (allocated(sinalph)) deallocate(sinalph)
    if (allocated(cdefserr)) deallocate(cdefserr)
    if (allocated(nbupstream)) deallocate(nbupstream)
    if (allocated(wtsource)) deallocate(wtsource)
    if (allocated(wtdep)) deallocate(wtdep)
    if (allocated(targsrc)) deallocate(targsrc)
    if (allocated(targleak)) deallocate(targleak)
    if (allocated(refdist)) deallocate(refdist)
    if (allocated(bgplasopt)) deallocate(bgplasopt)
    if (allocated(neut2d_raw)) deallocate(neut2d_raw)
    if (allocated(sol22_power_ratio)) deallocate(sol22_power_ratio)
    if (allocated(nrat_used)) deallocate(nrat_used)
    if (allocated(wallpol)) deallocate(wallpol)
    if (allocated(injrind)) deallocate(injrind)
    if (allocated(injkind)) deallocate(injkind)
    if (allocated(cerr)) deallocate(cerr)
    if (allocated(cdeferr)) deallocate(cdeferr)

  end subroutine deallocate_mod_comtor

end module mod_comtor
