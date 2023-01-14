module mod_solcommon
  use debug_options
  implicit none

  !
  !     common block containing sol solution parameters
  !     - note - it is in the process of being reorganized so that
  !              related variables will be more closely associated.
  !
  !
  !     -*-fortran-*-
  ! common /solcommon/ n0,te0,ti0,mb,m0,m1,k0e,k0i,v0,pinf,lams,s0,lenr,lamr,prad0,&
  !     gamecor,gammae,gammai,gamcor,pae,pai,lenmom,ceicf,frr,lastm0,machs,vsepmin,deltam0,&
  !     m0res,origm0,ssepmin,ffric,actffric,actlenmom,lastiters,vgradtmp,vgradlast,lammom,&
  !     smom0,ringlen,halfringlen,rcxmom,initm0,imagp,lastimagp,nlast,slast,vgrad,pinf0,&
  !     snegerr,fnorm,sptscopy,intionsrc,momsrc,intmomsrc,hlim,netarg,nefinal,himag,&
  !     n1,te1,ti1,v1e2d,vpe2d,soffset,sxp,lenri,gamma0,gperpcor,targfact,recfrac,gperpcor2,&
  !     slower,gnet,rbnd,sbnd,rconst,r0init,rateperp0,areasum,pstatic0,intarea,pinpute,&
  !     pinputi,gtarg,e2dgtarg,vpg,n1center,n1final,simag,ionptarg,elecptarg,presstarg,&
  !     pae_end,pai_end,pae_start,pai_start,fnorm2,g_pfzsol,pe_pfzsol,pi_pfzsol,pr_pfzsol
  
  ! save /solcommon/
  ! slmod begin - new
  ! slmod end
  ! common /solcommon2/ ionsrc,intioniz,qesrc,intqe,qesum,qisrc,intqi,qisum,e2dm0,lensst,&
  !     lensfi,peicf,ssrcst,ssrcfi,ssrcmid,ssrclen,lastvel,tcxmom,s5gausslen,s5alph,&
  !     s5alph2,s5alph3,pnormfact,s5offset,s5startval,recsrc,intrecsrc,alfimp,talimp,ex1imp,&
  !     ex2imp,gperp,intgperp,tcxcut,nhs,nhs0,ths,ssrcdecay,padd,oldne,oldte,oldti,&
  !     nh2s,qid,intqid,zb,tcutatiz,tcutmliz,tcutrec,tcutcx,trefcx,tcutqe,tfloor,tmin,spow,&
  !     spow2,spowlen,spowbeg,timax,temax,dropfrac,sgperpbeg,sgperpend,gperpfrac,gperpbegf,&
  !     gperpendf,gperprat,smom_mult,qesrc_mult,radsrc_mult,pp_pow_dist,ppelecpow,&
  !     ppionpow,pp_press,radsrc,intrad,pradsum,halflen,croplen,gextra_src_start,gextra_src_stop,&
  !     gextra_sink_start,gextra_sink_stop,start_gextra_src,stop_gextra_src,start_gextra_sink,&
  !     stop_gextra_sink,gextra_mult,gextra_src,gextra_sink
  
  !
  !                        not real*16 (real*8)
  !
  ! save /solcommon2/
  ! common /solcommon3/graph,ndiv,graphaux,nptscopy,forcet,graphvel,founds,newnimag,&
  !     stopimag,lastiter,ringnum,pinavail,lensind,miter,graphran,pinnorm,velsw,irdebug,&
  !     e2dstart,ike2d_start,ike2d,n_extffric,fillopt,startn,debug_s22,title
  !
  ! save /solcommon3/
  ! common /solcommon4/ extffric,alg_ion_src_len
  !
  ! save /solcommon4/
  real*8,public :: n0,te0,ti0,mb,m0,m1,k0e,k0i,v0,pinf,gamecor,lams,s0,machs,imagp,&
       lastimagp,nlast,slast,vgrad,prad0,gammae,gammai,gamcor,pae,pai,ceicf,&
       lastm0,deltam0,m0res,origm0,lastiters,pinf0,vgradtmp,vgradlast,vsepmin,ssepmin,&
       ffric,actffric,actlenmom,actlammom,lammom,lenmom,ringlen,halfringlen,rcxmom,smom0,initm0,&
       snegerr,fnorm,hlim,netarg,nefinal,himag,n1,te1,ti1,v1e2d,vpe2d,soffset,sxp,&
       gamma0,gperpcor,targfact,recfrac,slower,gnet,r0init,rateperp0,pstatic0,pinpute,&
       pinputi,vpg,n1center,n1final,simag,pae_end,pai_end,fnorm2,qesum,qisum,e2dm0,lensst,&
       peicf
  real*8,public :: lenr,lamr,frr,lenri,lamri,frri
  real*8,public :: totprad
  
  real*8,public,allocatable :: ionsrc(:),sptscopy(:),intionsrc(:),momsrc(:),intmomsrc(:),&
       rbnd(:),sbnd(:),rconst(:),intarea(:),gtarg(:,:),areasum(:),ionptarg(:,:),elecptarg(:,:),&
       presstarg(:,:),g_pfzsol(:,:),pe_pfzsol(:,:),pi_pfzsol(:,:),pr_pfzsol(:,:),&
       intioniz(:),intqi(:),intqe(:),qesrc(:),qisrc(:)
  ! slmod begin
  !
  !     >       gperprat(mxspts),pae_start,pai_start,smom_mult,
  ! slmod end
  !
  real*8,public :: ssrcst,ssrcfi,ssrcmid,ssrclen,lensfi,lastvel,tcxmom,s5gausslen,s5alph,&
       s5alph2,s5alph3,pnormfact,s5offset,s5startval,alfimp,talimp,ex1imp,ex2imp,&
       tcxcut,ssrcdecay,padd,zb,tcutatiz,tcutmliz,tcutrec,tcutcx,trefcx,tcutqe,tfloor,tmin,&
       spow,spow2,spowlen,spowbeg,timax,temax,dropfrac,sgperpbeg,sgperpend,gperpfrac,&
       gperpcor2,gperpbegf,gperpendf,pradsum,pae_start,pai_start,smom_mult,halflen,croplen,&
       radsrc_mult,pp_pow_dist,ppelecpow,ppionpow,pp_press,qesrc_mult,gextra_src_start,&
       gextra_src_stop,gextra_sink_start,gextra_sink_stop,start_gextra_src,stop_gextra_src,&
       start_gextra_sink,stop_gextra_sink,gextra_mult,gextra_src,gextra_sink,epowsum,ipowsum
  real*8,public,allocatable :: recsrc(:),intrecsrc(:),gperp(:),intgperp(:),nhs(:),nhs0(:),&
       ths(:),oldne(:),oldte(:),oldti(:),intqid(:),qid(:),nh2s(:),e2dgtarg(:,:),&
       radsrc(:),intrad(:),gperprat(:),epowsrc(:),intepow(:),ipowsrc(:),intipow(:)
  !
  real*8,public,allocatable :: extffric(:,:),extradsrc(:,:)
  !
  real,public :: graphran,alg_ion_src_len
  !
  integer,public :: graph,ndiv,graphaux,miter,graphvel,ringnum,nptscopy,lensind,forcet,&
       velsw,pinnorm,e2dstart,irdebug,ike2d_start,startn,ike2d,fillopt,n_extffric,n_extradsrc
  !
  character,public :: title*80

  ! jdemod - inputs related to double profile particle source option (swion=16)
  !
  integer,public:: dblsrc_opt
  real,public :: dblsrc_frac, dblsrc1_p1, dblsrc1_p2, dblsrc2_p1, dblsrc2_p2
  real*8,public :: ssrcst2, ssrcfi2, s02, ssrclen2, ssrcmid2
  
  ! jdemod - add initial file names for epow and ipow external power data (These can be replaced using *297,*298
  character*256 :: ext_epow_fn = 'ext_epow_data.txt'
  character*256 :: ext_ipow_fn = 'ext_ipow_data.txt'
  
  
  ! Moved from mod_slcom
  real,public :: simag1,simag2
  integer,public :: ierror
  integer,public :: sol22_osm_mode = 0 ! local copy of osm_mode from slcom
  integer,public :: sol22_outmode = 0  ! local copy of outmode from slcom


  integer,public :: sol22_cprint = 0 ! local Sol22 copy of print option selected
  integer,public :: sol22_print  = 0 ! independent sol22 input for selecting data to print as the code runs
  
  real,public,parameter :: sol22_machhi = 1.0e37, sol22_hi = 1.0e37
  
  
  logical,public :: founds,newnimag,stopimag,lastiter,pinavail,debug_s22=.false.
  !logical,public :: founds,newnimag,stopimag,lastiter,pinavail,debug_s22=.true.

  public :: allocate_mod_solcommon,deallocate_mod_solcommon,init_solcommon

contains

  subroutine init_solcommon(osm_mode,outmode)
    integer :: osm_mode,outmode
    sol22_osm_mode = osm_mode
    sol22_outmode = outmode
  end subroutine init_solcommon
  
  subroutine allocate_mod_solcommon
    !use mod_params
    use mod_solparams
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_solcommon','ALLOCATE')

    call allocate_array(ionsrc,mxspts,'ionsrc',ierr)
    call allocate_array(sptscopy,mxspts,'sptscopy',ierr)
    call allocate_array(intionsrc,mxspts,'intionsrc',ierr)
    call allocate_array(momsrc,mxspts,'momsrc',ierr)
    call allocate_array(intmomsrc,mxspts,'intmomsrc',ierr)
    call allocate_array(rbnd,0,'rbnd',mxspts,ierr)
    call allocate_array(sbnd,0,'sbnd',mxspts,ierr)
    call allocate_array(rconst,mxspts,'rconst',ierr)
    call allocate_array(intarea,mxspts,'intarea',ierr)
    call allocate_array(gtarg,mxspts,3,'gtarg',ierr)
    call allocate_array(areasum,mxspts,'areasum',ierr)
    call allocate_array(ionptarg,mxspts,3,'ionptarg',ierr)
    call allocate_array(elecptarg,mxspts,3,'elecptarg',ierr)
    call allocate_array(presstarg,mxspts,3,'presstarg',ierr)
    call allocate_array(g_pfzsol,mxspts,3,'g_pfzsol',ierr)
    call allocate_array(pe_pfzsol,mxspts,3,'pe_pfzsol',ierr)
    call allocate_array(pi_pfzsol,mxspts,3,'pi_pfzsol',ierr)
    call allocate_array(pr_pfzsol,mxspts,3,'pr_pfzsol',ierr)
    call allocate_array(intioniz,mxspts,'intioniz',ierr)
    call allocate_array(intqi,mxspts,'intqi',ierr)
    call allocate_array(intqe,mxspts,'intqe',ierr)
    call allocate_array(qesrc,mxspts,'qesrc',ierr)
    call allocate_array(qisrc,mxspts,'qisrc',ierr)
    call allocate_array(recsrc,mxspts,'recsrc',ierr)
    call allocate_array(intrecsrc,mxspts,'intrecsrc',ierr)
    call allocate_array(gperp,mxspts,'gperp',ierr)
    call allocate_array(intgperp,mxspts,'intgperp',ierr)
    call allocate_array(nhs,mxspts,'nhs',ierr)
    call allocate_array(nhs0,mxspts,'nhs0',ierr)
    call allocate_array(ths,mxspts,'ths',ierr)
    call allocate_array(oldne,mxspts,'oldne',ierr)
    call allocate_array(oldte,mxspts,'oldte',ierr)
    call allocate_array(oldti,mxspts,'oldti',ierr)
    call allocate_array(intqid,mxspts,'intqid',ierr)
    call allocate_array(qid,mxspts,'qid',ierr)
    call allocate_array(nh2s,mxspts,'nh2s',ierr)
    call allocate_array(e2dgtarg,mxspts,3,'e2dgtarg',ierr)
    call allocate_array(radsrc,mxspts,'radsrc',ierr)
    call allocate_array(intrad,mxspts,'intrad',ierr)
    call allocate_array(gperprat,mxspts,'gperprat',ierr)
    call allocate_array(extffric,mxspts,7,'extffric',ierr)
    call allocate_array(extradsrc,mxspts,7,'extradsrc',ierr)
    call allocate_array(epowsrc,mxspts,'epowsrc',ierr)
    call allocate_array(intepow,mxspts,'intepow',ierr)
    call allocate_array(ipowsrc,mxspts,'ipowsrc',ierr)
    call allocate_array(intipow,mxspts,'intipow',ierr)

  end subroutine allocate_mod_solcommon


  subroutine deallocate_mod_solcommon
    implicit none

    call pr_trace('mod_solcommon','DEALLOCATE')

    if (allocated(ionsrc)) deallocate(ionsrc)
    if (allocated(sptscopy)) deallocate(sptscopy)
    if (allocated(intionsrc)) deallocate(intionsrc)
    if (allocated(momsrc)) deallocate(momsrc)
    if (allocated(intmomsrc)) deallocate(intmomsrc)
    if (allocated(rbnd)) deallocate(rbnd)
    if (allocated(sbnd)) deallocate(sbnd)
    if (allocated(rconst)) deallocate(rconst)
    if (allocated(intarea)) deallocate(intarea)
    if (allocated(gtarg)) deallocate(gtarg)
    if (allocated(areasum)) deallocate(areasum)
    if (allocated(ionptarg)) deallocate(ionptarg)
    if (allocated(elecptarg)) deallocate(elecptarg)
    if (allocated(presstarg)) deallocate(presstarg)
    if (allocated(g_pfzsol)) deallocate(g_pfzsol)
    if (allocated(pe_pfzsol)) deallocate(pe_pfzsol)
    if (allocated(pi_pfzsol)) deallocate(pi_pfzsol)
    if (allocated(pr_pfzsol)) deallocate(pr_pfzsol)
    if (allocated(intioniz)) deallocate(intioniz)
    if (allocated(intqi)) deallocate(intqi)
    if (allocated(intqe)) deallocate(intqe)
    if (allocated(qesrc)) deallocate(qesrc)
    if (allocated(qisrc)) deallocate(qisrc)
    if (allocated(recsrc)) deallocate(recsrc)
    if (allocated(intrecsrc)) deallocate(intrecsrc)
    if (allocated(gperp)) deallocate(gperp)
    if (allocated(intgperp)) deallocate(intgperp)
    if (allocated(nhs)) deallocate(nhs)
    if (allocated(nhs0)) deallocate(nhs0)
    if (allocated(ths)) deallocate(ths)
    if (allocated(oldne)) deallocate(oldne)
    if (allocated(oldte)) deallocate(oldte)
    if (allocated(oldti)) deallocate(oldti)
    if (allocated(intqid)) deallocate(intqid)
    if (allocated(qid)) deallocate(qid)
    if (allocated(nh2s)) deallocate(nh2s)
    if (allocated(e2dgtarg)) deallocate(e2dgtarg)
    if (allocated(radsrc)) deallocate(radsrc)
    if (allocated(intrad)) deallocate(intrad)
    if (allocated(gperprat)) deallocate(gperprat)
    if (allocated(extffric)) deallocate(extffric)
    if (allocated(extradsrc)) deallocate(extradsrc)

  end subroutine deallocate_mod_solcommon

end module mod_solcommon
