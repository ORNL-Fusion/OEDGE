module mod_lim3_local
  use debug_options
  implicit none

  
  
  real,public :: vy0,vy02,partim,p
  real,public :: statim,vfluid,twalln,tdep,rizb,rfail,tfail
  real,public :: spara,timmax,tneut,tatiz,twall,avppos,rdifft,edge2
  real,public :: frqtim,edge1,ssef,yeff,yfact,rres,tres,rres1
  real,public :: avxpos,avypos,ttmax,tcent,tbyond,tcut,ratiz,rneut,quant
  real,public :: tavxpos
  real,public :: fvycol,y,svy,absy,rwalln,rcent,rtmax
  real,public,allocatable :: rdep(:),rwall(:)
  real,public :: absp
  real,public :: porm,oldy
  real,public,allocatable :: temp(:),yldtot(:),yldmax(:)
  real,public :: oldp
  
  real,public :: sputy,rmach,energy,rneut1,ryield,oldalp
  real,public,allocatable :: ythtot(:)
  ! slmod begin
  !
  ! moved to common block slcom:
  !
  !      real      svhins(-maxqxs:maxqxs),seyins(-maxqxs:maxqxs,maxizs)
  !
  ! slmod end
  real,public :: tmp_oldy,tmp_y
  real,public :: gytot1,gtot1,tbelow,spunew,randep
  !
  !     jdemod
  !
  real,public :: rstmax_win,rtime     ! number of time steps to max start time for particles - initial time of particle

  !
  !      real      mat1,mat2
  !
  real,public :: rstruk,tstruk,temold,fact,ran,emax
  !
  integer,public :: mat1,mat2
  real,public,allocatable :: tptrac(:,:)
  integer,public :: kklim,kk,natiz,nprod,ip,ifate,status
  integer,public,allocatable :: icut(:)
  integer,public :: iqx,iqy,ix,iy,iz,maxciz,ic,ii,ioy,iod,io
  integer,public :: imp,implim,matlim,j,jy,jx,it,mput,in
  real,public :: timusd,xm,ym
  real,public,allocatable :: polods(:)
  real,public,allocatable :: svpols(:)
  real,public,allocatable :: rions(:),stots(:)
  real,public :: cistot,cismax,rstmin,tstepl,rconst,svybit,avapos
  real,public :: qfact,yycon,yy,alpha
  real,public,allocatable :: sdtzs(:)
  real,public,allocatable :: tsplit(:),trulet(:),sdyzs(:)
  real,public :: factdeps
  real,public :: svg,svymin,svymod
  real,public :: dpprob
  integer,public :: is,iput
  integer,public,allocatable :: nsplit(:),nrulet(:),iget(:)
  integer,public :: traclen
  !
  logical,public :: diffus,resput,bigtrac
  !
  !     add some local variables related to calculating the scaling of the nerods3 data
  !
  integer,public :: perc
  !
  !     add iqy_tmp to support variable wall location
  !
  real,public :: pbnd1,pbnd2,local_pwid
  
  integer ,public :: iqy_tmp
  !
  !     add logical to record if splitting and rouletting is active to avoid
  !     a bug if alpha > cxspls(is) = 2*ca in one diffusive step
  !
  integer ,public :: ierr
  !
  !
  logical,public :: split
  double precision,public :: dsputy,dtemi,dqfact,deltax
  double precision,public,allocatable :: dtots(:)
  double precision,public :: dact,dwol,dsum4
  double precision,public,allocatable :: demp(:,:)
  !
  !     jdemod - add variables for recording forces
  !
  double precision,public :: dsum1,dsum2,dsum3,diz,dist
  double precision,public,allocatable :: douts(:,:)
  real,public :: ff,fe,feg,fig,fvh,fvel
  
  !
  !      double precision dy1,dy2
  !
  !     jdemod - change the calculation of the yposition
  !              at present, y is recalculated at each time step by combining the
  !              cumulative change in position stored in dy1 and dy2. in order for
  !              reflection to work - a new variable called y_position will hold the
  !              actual y_posiiton and dy1, dy2 -> delta_y1, delta_y2 will be the
  !              change in the current time step
  !
  !
  real,public :: ff2,fe2,fvh2
  !
  ! slmod begin
  double precision ,public :: y_position,old_y_position,delta_y1,delta_y2
  real,public :: ioncnt,ionpnt
  real,public :: ran1,ran2,rgauss,vpara,tpara,vparat
  real,public :: avgtrac
  real,public :: target

  public :: allocate_mod_lim3_local,deallocate_mod_lim3_local

contains

  subroutine allocate_mod_lim3_local
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_lim3_local','ALLOCATE')

    call allocate_array(rdep,2,'rdep',ierr)
    call allocate_array(rwall,2,'rwall',ierr)
    call allocate_array(temp,-maxnys,'temp',maxnys,ierr)
    call allocate_array(yldtot,2,'yldtot',ierr)
    call allocate_array(yldmax,2,'yldmax',ierr)
    call allocate_array(ythtot,2,'ythtot',ierr)
    call allocate_array(tptrac,maxlen,2,'tptrac',ierr)
    call allocate_array(icut,2,'icut',ierr)
    call allocate_array(polods,-maxqxs,'polods',maxqxs,ierr)
    call allocate_array(svpols,-maxqxs,'svpols',maxqxs,ierr)
    call allocate_array(rions,maxizs,'rions',ierr)
    call allocate_array(stots,20,'stots',ierr)
    call allocate_array(sdtzs,maxizs,'sdtzs',ierr)
    call allocate_array(tsplit,maxins,'tsplit',ierr)
    call allocate_array(trulet,maxins,'trulet',ierr)
    call allocate_array(sdyzs,maxizs,'sdyzs',ierr)
    call allocate_array(nsplit,maxins,'nsplit',ierr)
    call allocate_array(nrulet,maxins,'nrulet',ierr)
    call allocate_array(iget,0,'iget',maxput,ierr)
    call allocate_array(dtots,20,'dtots',ierr)
    call allocate_array(demp,-maxnys,maxnys,1,5,'demp',ierr)
    call allocate_array(douts,maxizs,10,'douts',ierr)

  end subroutine allocate_mod_lim3_local


  subroutine deallocate_mod_lim3_local
    implicit none

    call pr_trace('mod_lim3_local','DEALLOCATE')

    if (allocated(rdep)) deallocate(rdep)
    if (allocated(rwall)) deallocate(rwall)
    if (allocated(temp)) deallocate(temp)
    if (allocated(yldtot)) deallocate(yldtot)
    if (allocated(yldmax)) deallocate(yldmax)
    if (allocated(ythtot)) deallocate(ythtot)
    if (allocated(tptrac)) deallocate(tptrac)
    if (allocated(icut)) deallocate(icut)
    if (allocated(polods)) deallocate(polods)
    if (allocated(svpols)) deallocate(svpols)
    if (allocated(rions)) deallocate(rions)
    if (allocated(stots)) deallocate(stots)
    if (allocated(sdtzs)) deallocate(sdtzs)
    if (allocated(tsplit)) deallocate(tsplit)
    if (allocated(trulet)) deallocate(trulet)
    if (allocated(sdyzs)) deallocate(sdyzs)
    if (allocated(nsplit)) deallocate(nsplit)
    if (allocated(nrulet)) deallocate(nrulet)
    if (allocated(iget)) deallocate(iget)
    if (allocated(dtots)) deallocate(dtots)
    if (allocated(demp)) deallocate(demp)
    if (allocated(douts)) deallocate(douts)

  end subroutine deallocate_mod_lim3_local

end module mod_lim3_local
