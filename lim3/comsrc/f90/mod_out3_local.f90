module mod_out3_local
  use debug_options
  implicit none

  integer,public :: maxib,maxqts


  integer ,public,allocatable :: pizs(:),igzs(:),igls(:)
  integer ,public,allocatable :: igss(:),igts(:)
  
  real,public,allocatable :: oyvout(:),oyvwid(:)
  real,public,allocatable :: tptracx(:),tptracy(:)
  real,public,allocatable :: yfuns(:),plams(:),bina(:),binb(:)
  real,public,allocatable :: sum1(:),sum2(:),xfuns(:,:)
  real,public,allocatable :: facta(:),factb(:)
  real,public,allocatable :: qxwids(:),qtwids(:)
  real,public,allocatable :: youts1(:),ywidss(:)
  real,public,allocatable :: tvals(:,:),qts(:)
  real,public,allocatable :: qxfuns(:,:),aux1(:,:)
  real,public,allocatable :: aux2(:,:),yss(:)
  real,public,allocatable :: yjnts(:,:),xyjnts(:,:)
  real,public,allocatable :: xjnts(:,:)
  real,public,allocatable :: xyints(:,:),totals(:,:)
  real,public,allocatable :: pouts(:),yints(:,:),xints(:,:)
  
  !
  !     variables for r,theta graphs
  !
  real,public,allocatable :: rigs(:)
  real,public,allocatable :: sints(:,:)
  real,public,allocatable :: tints(:,:)
  
  
  real,public,allocatable :: rints(:,:)

  character*36,allocatable,public :: tlabs(:),zlabs(:),zlabs1(:),plabs(:),vel_labs(:)
  character*7,allocatable,public :: prinps(:)
  
  public :: allocate_mod_out3_local,deallocate_mod_out3_local

contains

  subroutine allocate_mod_out3_local
    use mod_params
    use mod_rtheta
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_out3_local','ALLOCATE')

    ! calculate derived parameters
    MAXQTS=(MAXIZS+1)*MAXNTS*2
    MAXIB=MAXNXS*2*MAXNYS

    call allocate_array(pizs,maxnls,'pizs',ierr)
    call allocate_array(igzs,-2,'igzs',maxizs+1,ierr)
    call allocate_array(igls,maxnls,'igls',ierr)
    call allocate_array(igss,maxnls,'igss',ierr)
    call allocate_array(igts,maxib,'igts',ierr)
    call allocate_array(oyvout,maxos,'oyvout',ierr)
    call allocate_array(oyvwid,maxos,'oyvwid',ierr)
    call allocate_array(tptracx,maxlen,'tptracx',ierr)
    call allocate_array(tptracy,maxlen,'tptracy',ierr)
    call allocate_array(yfuns,-maxnys,'yfuns',maxnys,ierr)
    call allocate_array(plams,maxnls,'plams',ierr)
    call allocate_array(bina,maxib,'bina',ierr)
    call allocate_array(binb,maxib,'binb',ierr)
    call allocate_array(sum1,-2,'sum1',maxizs,ierr)
    call allocate_array(sum2,-2,'sum2',maxizs,ierr)
    call allocate_array(xfuns,maxnxs,2,'xfuns',ierr)
    call allocate_array(facta,-1,'facta',maxizs,ierr)
    call allocate_array(factb,-1,'factb',maxizs,ierr)
    call allocate_array(qxwids,-maxqxs,'qxwids',maxqxs,ierr)
    call allocate_array(qtwids,0,'qtwids',maxqts,ierr)
    call allocate_array(youts1,-maxnys,'youts1',maxnys,ierr)
    call allocate_array(ywidss,-maxnys,'ywidss',maxnys,ierr)
    call allocate_array(tvals,0,maxqts,-2,maxizs+1,'tvals',ierr)
    call allocate_array(qts,0,'qts',maxqts,ierr)
    call allocate_array(qxfuns,-maxqxs,maxqxs,1,2,'qxfuns',ierr)
    call allocate_array(aux1,1,maxnxs,1-maxnys,maxnys,'aux1',ierr)
    call allocate_array(aux2,1,maxnxs,1-maxnys,maxnys,'aux2',ierr)
    call allocate_array(yss,-maxnys,'yss',maxnys,ierr)
    call allocate_array(yjnts,maxnxs,maxnts+1,'yjnts',ierr)
    call allocate_array(xyjnts,-maxnps,maxnps,1,maxnts+1,'xyjnts',ierr)
    call allocate_array(xjnts,-maxnys,maxnys,1,maxnts+1,'xjnts',ierr)
    call allocate_array(xyints,-maxnps,maxnps,-2,maxizs+1,'xyints',ierr)
    call allocate_array(totals,1,4,-2,maxizs+1,'totals',ierr)
    call allocate_array(pouts,-maxnps,'pouts',maxnps,ierr)
    call allocate_array(yints,1,maxnxs,-2,maxizs+1,'yints',ierr)
    call allocate_array(xints,-maxnys,maxnys,-2,maxizs+1,'xints',ierr)
    call allocate_array(rigs,maxizs+2,'rigs',ierr)
    call allocate_array(sints,1,maxnss,-2,maxizs+1,'sints',ierr)
    call allocate_array(tints,1,maxnrs,-2,maxizs+1,'tints',ierr)
    call allocate_array(rints,1,maxnas,-2,maxizs+1,'rints',ierr)

    ! ALLOCATE character arrays that are not handled by allocate_array
    allocate(tlabs(maxnts+1))
    allocate(zlabs(-2:MAXIZS+1))
    allocate(zlabs1(-2:MAXIZS+1))
    allocate(vel_labs(-2:MAXIZS+1))
    allocate(prinps(-MAXNPS-1:MAXNPS))
    allocate(plabs(maxnls))
          
  end subroutine allocate_mod_out3_local


  subroutine deallocate_mod_out3_local
    implicit none

    call pr_trace('mod_out3_local','DEALLOCATE')

    if (allocated(pizs)) deallocate(pizs)
    if (allocated(igzs)) deallocate(igzs)
    if (allocated(igls)) deallocate(igls)
    if (allocated(igss)) deallocate(igss)
    if (allocated(igts)) deallocate(igts)
    if (allocated(oyvout)) deallocate(oyvout)
    if (allocated(oyvwid)) deallocate(oyvwid)
    if (allocated(tptracx)) deallocate(tptracx)
    if (allocated(tptracy)) deallocate(tptracy)
    if (allocated(yfuns)) deallocate(yfuns)
    if (allocated(plams)) deallocate(plams)
    if (allocated(bina)) deallocate(bina)
    if (allocated(binb)) deallocate(binb)
    if (allocated(sum1)) deallocate(sum1)
    if (allocated(sum2)) deallocate(sum2)
    if (allocated(xfuns)) deallocate(xfuns)
    if (allocated(facta)) deallocate(facta)
    if (allocated(factb)) deallocate(factb)
    if (allocated(qxwids)) deallocate(qxwids)
    if (allocated(qtwids)) deallocate(qtwids)
    if (allocated(youts1)) deallocate(youts1)
    if (allocated(ywidss)) deallocate(ywidss)
    if (allocated(tvals)) deallocate(tvals)
    if (allocated(qts)) deallocate(qts)
    if (allocated(qxfuns)) deallocate(qxfuns)
    if (allocated(aux1)) deallocate(aux1)
    if (allocated(aux2)) deallocate(aux2)
    if (allocated(yss)) deallocate(yss)
    if (allocated(yjnts)) deallocate(yjnts)
    if (allocated(xyjnts)) deallocate(xyjnts)
    if (allocated(xjnts)) deallocate(xjnts)
    if (allocated(xyints)) deallocate(xyints)
    if (allocated(totals)) deallocate(totals)
    if (allocated(pouts)) deallocate(pouts)
    if (allocated(yints)) deallocate(yints)
    if (allocated(xints)) deallocate(xints)
    if (allocated(rigs)) deallocate(rigs)
    if (allocated(sints)) deallocate(sints)
    if (allocated(tints)) deallocate(tints)
    if (allocated(rints)) deallocate(rints)
    
    if (allocated(tlabs)) deallocate(tlabs)
    if (allocated(zlabs)) deallocate(zlabs)
    if (allocated(zlabs1)) deallocate(zlabs1)
    if (allocated(vel_labs)) deallocate(vel_labs)
    if (allocated(prinps)) deallocate(prinps)
    if (allocated(plabs)) deallocate(plabs)
    
  end subroutine deallocate_mod_out3_local

end module mod_out3_local
