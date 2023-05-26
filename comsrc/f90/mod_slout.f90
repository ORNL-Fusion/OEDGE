module mod_slout
  use mod_params
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /drawcom/ plottype,plotnorm,slopt,char,map1x,map2x,map1y,map2y,opt_xscale,&
  !     opt_yscale,restoresolution,slopt2,slopt3,slopt4,slopt5
  
  ! save /drawcom/
  integer,public :: slopt,opt_xscale,opt_yscale,slopt2,slopt3,slopt4,slopt5
  integer,public,allocatable :: plottype(:),plotnorm(:)
  logical,public :: restoresolution
  
  real,public :: map1x,map2x,map1y,map2y
  
  character*100,public :: char(30)
  ! common /grmcom/ grm_opt,ngs2,elabs2,grm_shade,grm_cell,plottype2,ylab2
  
  ! save /grmcom/
  integer,public :: grm_opt
  integer,public,allocatable :: ngs2(:),plottype2(:,:)
  real,public,allocatable :: grm_shade(:,:),grm_cell(:,:)
  !character*36,public,allocatable :: elabs2   (8,maxngs)
  character(len=36),public,allocatable :: elabs2   (:,:)
   
  
  character*128,public :: ylab2    (-30:30)
  integer,public :: maxshow
  
  parameter (maxshow=10000)
  ! common /loscom/ wght0,losopt,nshow,rshow,zshow,ashow
  ! save /loscom/
  integer,public :: losopt,nshow
  
  real,public,allocatable :: wght0(:,:),rshow(:),zshow(:),ashow(:)
  ! common /timecom/ qt
  ! save  /timecom/
  
  real,public :: qt
  ! common /loadcom/ loadstep,stepopt,nsteplist,steplist
  ! save  /loadcom/
  integer,public :: loadstep,stepopt,nsteplist
  integer,public,allocatable :: steplist(:)
  
  ! common /gencom/ sldata
  ! save  /gencom/
  
  
  
  !      common /machplot/ machdat
  !      real              machdat(2,maxnrs)
  real,public :: sldata
  ! common /normcom/ nrmindex,nrmcalculate,nrmdata,nrmvalue,nrmi1,nrmi2,nrmtype,nrmr1,&
  !     nrmr2,nrmnum,nrmstep,nrmcomment
  ! save /normcom/
  integer,public :: nrmindex,nrmtype,nrmi1,nrmi2
  integer,public,allocatable :: nrmnum(:),nrmstep(:)
  logical,public :: nrmcalculate
  real,public :: nrmr1,nrmr2
  real,public,allocatable :: nrmdata(:,:),nrmvalue(:)
  character*1024,public :: nrmcomment(1024)

  public :: allocate_mod_slout,deallocate_mod_slout

contains

  subroutine allocate_mod_slout
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_slout','ALLOCATE')

    call allocate_array(plottype,maxngs,'plottype',ierr)
    call allocate_array(plotnorm,maxngs,'plotnorm',ierr)
    call allocate_array(ngs2,maxngs,'ngs2',ierr)
    call allocate_array(plottype2,8,maxngs,'plottype2',ierr)
    call allocate_array(grm_shade,2,maxngs,'grm_shade',ierr)
    call allocate_array(grm_cell,0,maxnks,1,maxngs,'grm_cell',ierr)
    call allocate_array(wght0,maxnks,maxnrs,'wght0',ierr)
    call allocate_array(rshow,maxshow,'rshow',ierr)
    call allocate_array(zshow,maxshow,'zshow',ierr)
    call allocate_array(ashow,maxshow,'ashow',ierr)
    call allocate_array(steplist,100,'steplist',ierr)
    call allocate_array(nrmnum,1024,'nrmnum',ierr)
    call allocate_array(nrmstep,1024,'nrmstep',ierr)
    call allocate_array(nrmdata,maxthe,2,'nrmdata',ierr)
    call allocate_array(nrmvalue,1024,'nrmvalue',ierr)

    call allocate_array(elabs2,8,maxngs,'E labels',ierr)

    
  end subroutine allocate_mod_slout


  subroutine deallocate_mod_slout
    implicit none

    call pr_trace('mod_slout','DEALLOCATE')

    if (allocated(plottype)) deallocate(plottype)
    if (allocated(plotnorm)) deallocate(plotnorm)
    if (allocated(ngs2)) deallocate(ngs2)
    if (allocated(plottype2)) deallocate(plottype2)
    if (allocated(grm_shade)) deallocate(grm_shade)
    if (allocated(grm_cell)) deallocate(grm_cell)
    if (allocated(wght0)) deallocate(wght0)
    if (allocated(rshow)) deallocate(rshow)
    if (allocated(zshow)) deallocate(zshow)
    if (allocated(ashow)) deallocate(ashow)
    if (allocated(steplist)) deallocate(steplist)
    if (allocated(nrmnum)) deallocate(nrmnum)
    if (allocated(nrmstep)) deallocate(nrmstep)
    if (allocated(nrmdata)) deallocate(nrmdata)
    if (allocated(nrmvalue)) deallocate(nrmvalue)

    if (allocated(elabs2)) deallocate(elabs2)
    
  end subroutine deallocate_mod_slout

end module mod_slout
