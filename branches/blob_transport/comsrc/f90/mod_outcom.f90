module mod_outcom
  use debug_options
  use mod_slcom
  use mod_slout
  implicit none

  !
  !---------------------------------------------------------------------
  !
  !     include steve's common blocks - not sure where they are needed
  !
  !     -*-fortran-*-
  ! include 'slcom'
  !
  !---------------------------------------------------------------------
  !
  !     smoothing and graphics print options - to trace routines
  !
  !---------------------------------------------------------------------
  !
  ! include 'slout'
  ! common /nsmooth/ numsmooth,cgrprint
  ! save /nsmooth/
  !
  !---------------------------------------------------------------------
  !
  !     ghost options?
  !
  !---------------------------------------------------------------------
  !
  integer,public :: numsmooth,cgrprint
  ! common /ghostcom/ iopt_ghost
  ! save /ghostcom/
  !
  !---------------------------------------------------------------------
  !
  !     input values to out
  !
  !---------------------------------------------------------------------
  !
  integer,public :: iopt_ghost
  ! common/outcom_input/ alphae,zadj,rxygrid,plrpopt
  ! save /outcom_input/
  real,public :: alphae
  real,public :: zadj
  !
  !---------------------------------------------------------------------
  !
  !     adding extra comments to plots:
  !
  !---------------------------------------------------------------------
  !
  integer,public :: rxygrid,plrpopt
  ! common/outcom_extcom/ ncomments,extra_comments
  ! save /outcom_extcom/
  character*30,public :: extra_comments(10)
  !
  !---------------------------------------------------------------------
  !
  !     colour and other plot characteristics
  !
  !---------------------------------------------------------------------
  !
  integer,public :: ncomments
  ! common/outcom_plotvals/ n_cols,col_opt,ngs,nplots,ismoth,izmin,izmax,iter,niters,&
  !     nizs,clsup,magstep,mgst,mgnd,enldist,nplts,pltmax,pltmin,pltfact,pltmins,pltmaxs,&
  !     mode
  ! save /outcom_plotvals/
  integer,public :: n_cols,col_opt
  integer,public :: ngs,nplots,ismoth
  integer,public :: izmin,izmax
  integer,public :: iter,niters,nizs
  integer,public :: clsup
  real,public :: magstep,mgst,mgnd
  real,public :: enldist
  integer,public :: nplts
  real,public :: pltmax,pltmin,pltfact
  real,public,allocatable :: pltmins(:),pltmaxs(:)
  !
  !     secondary plot related variables - mostly related to los plots
  !
  !
  integer,public :: mode
  ! common/outcom_secondary_plotvals/ ignors,itec,navs,avs,vmin,vmax,numthe,avpts,atype,&
  !     zobs,robs,drad,dthe,themin,themax,themin_start,theres,mfact,pswitch
  ! save /outcom_secondary_plotvals/
  integer,public :: itec,navs
  integer,public,allocatable :: ignors(:)
  real,public :: vmin,vmax
  real,public,allocatable :: avs(:)
  integer,public :: numthe,avpts,atype
  real,public :: zobs,robs,drad,dthe,themin,themax,themin_start
  real,public :: theres,mfact
  !
  !---------------------------------------------------------------------
  !
  !     cngs - option added to allow the number of contours to
  !     be externally specified. slowly convert from fixed value
  !     that is hard-coded to this system.
  !
  !---------------------------------------------------------------------
  !
  logical,public :: pswitch
  ! common/outcom_contour/ cngs,cntropt,global_cngs,localcngs,icntr,minscale,maxscale,&
  !     scalef
  ! save /outcom_contour/
  integer,public :: cngs,cntropt,global_cngs
  integer,public :: localcngs
  real,public :: minscale,maxscale
  integer,public :: icntr
  !
  !     power contour array declarations
  !
  real,public :: scalef
  ! common/outcom_power_contours/ nconts,nclev,conts,clev,pradclev
  ! save /outcom_power_contours/
  integer,public :: nconts,nclev
  real,public,allocatable :: conts(:),clev(:)
  !
  !---------------------------------------------------------------------
  !
  !     global scaling factors
  !
  !---------------------------------------------------------------------
  !
  real,public,allocatable :: pradclev(:)
  ! common/outcom_factors/ facta,factb,ft,fp,fact
  ! save /outcom_factors/
  real,public,allocatable :: facta(:),factb(:)
  real,public :: fp,ft
  !
  !---------------------------------------------------------------------
  !
  !     plrp related data
  !
  !---------------------------------------------------------------------
  !
  real,public :: fact
  ! common/outcom_plrp/ plams,pizs,pind,plrpcnt
  ! save /outcom_plrp/
  real,public,allocatable :: plams(:)
  !
  !---------------------------------------------------------------------
  !
  !     out string declarations - labels and data from the raw file
  !
  !---------------------------------------------------------------------
  !
  !     jdemod - the job string is actually read in from the raw
  !              data file - however, it has become common practice
  !              to use it for other comments since it is directly
  !              used in the call to draw. however, the original
  !              information in the string includes the date, time
  !              and divimp version number for the original run - which
  !              is worth - keeping. the job string is now copied to
  !              job_saved after the call to get.
  !
  integer,public :: plrpcnt
  integer,public,allocatable :: pizs(:),pind(:)
  ! common/outcom_strings/ title,job,job_saved,graph1,equil,xlab,ylab,xpoint,table,&
  !     zlabs,ref,plane,anly,nview,smooth,name,elabs,plabs,klab,datatitle
  ! save /outcom_strings/
  character,public :: title*174,job*72,job_saved*72,graph1*80
  !character desc*1024
  !
  !---- character arrays for draw routines etc
  !
  character,public :: equil*60
  character*36,public :: xlab,ylab,xpoint
  !character*36,public :: table,zlabs(-2:maxizs+1)
  character*36,public :: table
  character*36,public,allocatable :: zlabs(:)
  character*44,public :: ref,plane,anly,nview
  character*72,public :: smooth
  character*36,public :: name,elabs(maxngs)
  character*36,public :: plabs(-2:maxplrp),klab
  !
  !---------------------------------------------------------------------
  !
  !     content summaries
  !
  !---------------------------------------------------------------------
  !
  character*60,public :: datatitle
  ! common/outcom_content/ core_content, core_area,edge_content, edge_area, pp_content,&
  !      pp_area,div_content, div_area, main_content, main_area
  ! save /outcom_content/
  real,public :: core_content,core_area
  real,public :: edge_content,edge_area
  real,public :: pp_content,pp_area
  real,public :: div_content,div_area
  !
  !---------------------------------------------------------------------
  !
  !     adas related data
  !
  !---------------------------------------------------------------------
  !
  real,public :: main_content,main_area
  ! common/outcom_adas/  tadas,dadas,pecae,pecar,pecax,ltrng,ldrng,wlngth,isele,iselr,&
  !     iselx,iseld,iseldef,iadas,npairs,ircode,adasyr,adasyr2,isele2,iselr2,iselx2,iseld2,&
  !     iz_state,z_atom,iz_state2,z_atom2
  ! save  /outcom_adas/
  integer,public :: isele,iselr,iselx,iseld,iseldef
  integer,public :: iadas,npairs,ircode
  real*8,public,allocatable :: tadas(:),dadas(:)
  real*8,public,allocatable :: pecae(:),pecar(:),pecax(:)
  logical*4,public,allocatable :: ltrng(:),ldrng(:)
  real,public :: wlngth
  integer,public :: adasyr
  integer,public :: adasyr2,isele2,iselr2,iselx2,iseld2
  !
  integer,public :: iz_state,z_atom,iz_state2,z_atom2
  ! common/outcom_adas_strings/ graph2,graph3,graph4,hlabs,blabs,adasid,pectitle,plabad,&
  !     adasex,adasgr,adasty,graph5,adasid2,adasex2
  ! save /outcom_adas_strings/
  character,public :: graph2*80,graph3*80,graph4*80
  character*36,public :: hlabs(0:2),blabs
  character,public :: adasid*80,pectitle*120,plabad*36
  character,public :: adasex*3
  character,public :: adasgr*8,adasty*80
  !
  !---------------------------------------------------------------------
  !
  !     mean free path and scale length variables.
  !
  !---------------------------------------------------------------------
  !
  character,public :: graph5*80,adasid2*80,adasex2*3
  ! common/outcom_mfp/ lgradte,lgradti,lmfpii,lmfpee
  ! save /outcom_mfp/
  real,public,allocatable :: lgradte(:,:),lgradti(:,:)
  !
  !---------------------------------------------------------------------
  !
  !     output arrays for passing the calculated results to the
  !     plotting routines.
  !
  !---------------------------------------------------------------------
  !
  real,public,allocatable :: lmfpii(:,:),lmfpee(:,:)
  ! common/outcom_outvals/ touts,twids,tvals,douts,dwids,dvals,kouts,kwids,kvals,louts,&
  !     lwids,lvals,cvalsa,plastmp,xouts,xwids,xvals,youts,ywids,yvals,cvals
  ! save /outcom_outvals/
  real,public,allocatable :: touts(:),twids(:),tvals(:,:)
  real,public,allocatable :: douts(:),dwids(:),dvals(:,:)
  real,public,allocatable :: kouts(:),kwids(:),kvals(:,:)
  real,public,allocatable :: louts(:),lwids(:),lvals(:,:)
  real,public,allocatable :: cvalsa(:,:)
  !
  !     old-style contour plots ...
  !
  real,public,allocatable :: plastmp(:,:)
  real,public,allocatable :: xouts(:),xwids(:),xvals(:,:)
  real,public,allocatable :: youts(:),ywids(:),yvals(:,:)
  !
  !---------------------------------------------------------------------
  !
  !     values related to zoom or portion of 2d plot to be displayed
  !
  !---------------------------------------------------------------------
  !
  real,public,allocatable :: cvals(:,:)
  ! common/outcom_zoom/ xxmin,xxmax,yymin,yymax,xnear2,ynear2,xcen,ycen,zminp,zmaxp,&
  !     zmode,xnear,ynear
  ! save /outcom_zoom/
  real,public :: xxmin,xxmax,yymin,yymax
  real,public :: xnear2,ynear2,xcen,ycen
  real,public :: zminp,zmaxp
  integer,public :: zmode
  !
  !---------------------------------------------------------------------
  !
  !     experimental data related values
  !     - used to hold a list of experimantal dataset indices to be
  !       included on plots if possible
  !
  !---------------------------------------------------------------------
  !
  real,public :: xnear,ynear
  integer,public :: max_expt_datasets
  parameter (max_expt_datasets=10)
  ! common /expt_data_list/ expt_nsets,expt_datasets
  ! save /expt_data_list/
  
  !
  ! ammod begin.
  !---------------------------------------------------------------------
  !
  ! values related to hydrocarbon outout.
  !
  !---------------------------------------------------------------------
  integer,public :: expt_nsets
  integer,public,allocatable :: expt_datasets(:)
  ! common/outcom_hc/ start_hc_species, end_hc_species
  ! save /outcom_hc/
  integer ,public :: start_hc_species
  
  ! ammod end.
  !
  !---------------------------------------------------------------------
  !
  !     local variables
  !
  !---------------------------------------------------------------------
  !
  integer ,public :: end_hc_species
  real,public :: atan2c
  external atan2c
  
  
  logical,public :: griderr
  integer,public :: len,lenstr
  
  external lenstr
  
  !
  !     iplot - move to params
  !
  real,public,allocatable :: ktmp(:,:)
  integer,public :: iplot
  
  parameter (iplot=49)
  character,public :: xfesym*2
  
  external  xfesym

  public :: allocate_mod_outcom,deallocate_mod_outcom

contains

  subroutine allocate_mod_outcom
    use mod_params
    use allocate_arrays
    use error_handling
    implicit none
    integer :: ierr

    call pr_trace('mod_outcom','ALLOCATE')

    call allocate_array(pltmins,maxplts,'pltmins',ierr)
    call allocate_array(pltmaxs,maxplts,'pltmaxs',ierr)
    call allocate_array(ignors,maxngs,'ignors',ierr)
    call allocate_array(avs,0,'avs',100,ierr)
    call allocate_array(conts,maxpts,'conts',ierr)
    call allocate_array(clev,maxpts,'clev',ierr)
    call allocate_array(pradclev,0,'pradclev',maxizs+1,ierr)
    call allocate_array(facta,-1,'facta',maxizs,ierr)
    call allocate_array(factb,-1,'factb',maxizs,ierr)
    call allocate_array(plams,-1,'plams',maxplrp,ierr)
    call allocate_array(pizs,-1,'pizs',maxplrp,ierr)
    call allocate_array(pind,-1,'pind',maxizs+1,ierr)
    call allocate_array(tadas,20,'tadas',ierr)
    call allocate_array(dadas,20,'dadas',ierr)
    call allocate_array(pecae,20,'pecae',ierr)
    call allocate_array(pecar,20,'pecar',ierr)
    call allocate_array(pecax,20,'pecax',ierr)
    call allocate_array(ltrng,20,'ltrng',ierr)
    call allocate_array(ldrng,20,'ldrng',ierr)
    call allocate_array(lgradte,maxnks,maxnrs,'lgradte',ierr)
    call allocate_array(lgradti,maxnks,maxnrs,'lgradti',ierr)
    call allocate_array(lmfpii,maxnks,maxnrs,'lmfpii',ierr)
    call allocate_array(lmfpee,maxnks,maxnrs,'lmfpee',ierr)
    call allocate_array(touts,maxthe,'touts',ierr)
    call allocate_array(twids,maxthe,'twids',ierr)
    call allocate_array(tvals,maxthe,maxngs,'tvals',ierr)
    call allocate_array(douts,maxnds+2,'douts',ierr)
    call allocate_array(dwids,maxnds+2,'dwids',ierr)
    call allocate_array(dvals,maxnds+2,maxngs,'dvals',ierr)
    call allocate_array(kouts,maxnks,'kouts',ierr)
    call allocate_array(kwids,maxnks,'kwids',ierr)
    call allocate_array(kvals,maxnks,maxngs,'kvals',ierr)
    call allocate_array(louts,maxseg,'louts',ierr)
    call allocate_array(lwids,maxseg,'lwids',ierr)
    call allocate_array(lvals,maxseg,maxngs,'lvals',ierr)
    call allocate_array(cvalsa,maxnks,maxnrs,'cvalsa',ierr)
    call allocate_array(plastmp,maxnks,maxnrs,'plastmp',ierr)
    call allocate_array(xouts,maxgxs,'xouts',ierr)
    call allocate_array(xwids,maxgxs,'xwids',ierr)
    call allocate_array(xvals,maxgxs,maxngs,'xvals',ierr)
    call allocate_array(youts,maxgys,'youts',ierr)
    call allocate_array(ywids,maxgys,'ywids',ierr)
    call allocate_array(yvals,maxgys,maxngs,'yvals',ierr)
    call allocate_array(cvals,maxgxs,maxgys,'cvals',ierr)
    call allocate_array(expt_datasets,max_expt_datasets,'expt_datasets',ierr)
    call allocate_array(ktmp,maxnks,maxnrs,'ktmp',ierr)

    ! allocate zlabs explicitly since it is a character array

    if (allocated(zlabs)) deallocate(zlabs)
    allocate(zlabs(-2:maxizs+1),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array ZLABS : IERR =',ierr)
    endif
    
  end subroutine allocate_mod_outcom


  subroutine deallocate_mod_outcom
    implicit none

    call pr_trace('mod_outcom','DEALLOCATE')

    if (allocated(pltmins)) deallocate(pltmins)
    if (allocated(pltmaxs)) deallocate(pltmaxs)
    if (allocated(ignors)) deallocate(ignors)
    if (allocated(avs)) deallocate(avs)
    if (allocated(conts)) deallocate(conts)
    if (allocated(clev)) deallocate(clev)
    if (allocated(pradclev)) deallocate(pradclev)
    if (allocated(facta)) deallocate(facta)
    if (allocated(factb)) deallocate(factb)
    if (allocated(plams)) deallocate(plams)
    if (allocated(pizs)) deallocate(pizs)
    if (allocated(pind)) deallocate(pind)
    if (allocated(tadas)) deallocate(tadas)
    if (allocated(dadas)) deallocate(dadas)
    if (allocated(pecae)) deallocate(pecae)
    if (allocated(pecar)) deallocate(pecar)
    if (allocated(pecax)) deallocate(pecax)
    if (allocated(ltrng)) deallocate(ltrng)
    if (allocated(ldrng)) deallocate(ldrng)
    if (allocated(lgradte)) deallocate(lgradte)
    if (allocated(lgradti)) deallocate(lgradti)
    if (allocated(lmfpii)) deallocate(lmfpii)
    if (allocated(lmfpee)) deallocate(lmfpee)
    if (allocated(touts)) deallocate(touts)
    if (allocated(twids)) deallocate(twids)
    if (allocated(tvals)) deallocate(tvals)
    if (allocated(douts)) deallocate(douts)
    if (allocated(dwids)) deallocate(dwids)
    if (allocated(dvals)) deallocate(dvals)
    if (allocated(kouts)) deallocate(kouts)
    if (allocated(kwids)) deallocate(kwids)
    if (allocated(kvals)) deallocate(kvals)
    if (allocated(louts)) deallocate(louts)
    if (allocated(lwids)) deallocate(lwids)
    if (allocated(lvals)) deallocate(lvals)
    if (allocated(cvalsa)) deallocate(cvalsa)
    if (allocated(plastmp)) deallocate(plastmp)
    if (allocated(xouts)) deallocate(xouts)
    if (allocated(xwids)) deallocate(xwids)
    if (allocated(xvals)) deallocate(xvals)
    if (allocated(youts)) deallocate(youts)
    if (allocated(ywids)) deallocate(ywids)
    if (allocated(yvals)) deallocate(yvals)
    if (allocated(cvals)) deallocate(cvals)
    if (allocated(expt_datasets)) deallocate(expt_datasets)
    if (allocated(ktmp)) deallocate(ktmp)

    ! deallocate zlabs
    if (allocated(zlabs)) deallocate(zlabs)

  end subroutine deallocate_mod_outcom

end module mod_outcom
