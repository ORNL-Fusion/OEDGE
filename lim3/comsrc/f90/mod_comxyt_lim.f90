module mod_comxyt
  
  use mod_params
  
  
  !c
  !      common /comxyt/ nxs,nys,nqxso,nqxsi,nqys,nts,xscalo,xscali,yscale,
  !     >                ny3d,cstmax,xs,ys,qxs,qys,iqxs,iqys,xwids,ywids,
  !     >                ps,xcyls,dwelts,dwelfs,ctimes,pwids,qs,xouts,
  !     >                youts,delps,ctmax,cxdeep
  !      integer         nxs,nys,nqxso,nqxsi,nqys,nts,iqys(maxnys),ny3d
  !      integer         iqxs(maxnxs)
  !      real            ctimes(maxnts+1,0:maxizs),cstmax,xouts(maxnxs)
  !      real            xs(maxnxs),ys(maxnys),qys(maxqys),xscali,yscale
  !      real            qxs(-maxqxs:maxqxs),xscalo
  !      real            xwids(maxnxs),ywids(maxnys),ps(-maxnps:maxnps)
  !      real            xcyls(maxnxs),dwelts(0:maxizs),dwelfs(maxnts)
  !      real            pwids(-maxnps:maxnps),qs(-maxqxs:maxqxs)
  !      real            youts(-maxnys:maxnys),delps(maxnxs,maxnys)
  !      real            ctmax,cxdeep
  !
  
  
  implicit none
  private
  
  
  integer,public:: nxs,nys,nqxso,nqxsi,nqys,nts,ny3d
  integer,public,allocatable:: iqys(:)
  integer,public,allocatable:: iqxs(:)
  real,public:: cstmax
  real,public,allocatable:: ctimes(:,:),xouts(:)
  real,public:: xscali,yscale
  real,public,allocatable:: xs(:),ys(:),qys(:)
  real,public:: xscalo
  real,public,allocatable:: qxs(:)
  real,public,allocatable:: xwids(:),ywids(:)
  real,public,allocatable:: xcyls(:),dwelts(:),dwelfs(:)
  real,public,allocatable:: qs(:)
  real,public,allocatable:: youts(:),delps(:,:)
  real,public:: ctmax,cxdeep
  
  !real,public::            pwids(-maxnps:maxnps),ps(-maxnps:maxnps)
  
  real, allocatable,public :: pwids(:),ps(:)
  
  integer, public:: nsurf,npbins
  real, allocatable,public:: pbin_bnds(:),surf_bnds(:,:)
  integer, allocatable,public:: pzone(:)
  
  public :: allocate_mod_comxyt, deallocate_mod_comxyt
  
  
contains
  
  subroutine allocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(dtev  ,maxnxs,'dtev',ierr)
    ! 1d arrays with both bounds specified have a different argumemnt order to differentiate signatures.
    call allocate_array(ps,-maxnps,'p bin widths',maxnps,ierr)
    call allocate_array(pwids,-maxnps,'p bin widths',maxnps,ierr)
    call allocate_array(pzone,-maxnps,'p bin widths',maxnps,ierr)
    call allocate_array(pbin_bnds,2*maxnps+1,'p bin boundaries',ierr)  ! allow extra one to include 0.0 if desired
    call allocate_array(surf_bnds,max_nsurf,2,'poloidal limiter surface bounds',ierr)  ! surface is defined as p1 to p2 - sets should not overlap
    write(0,*) 'pbin_bnds allocated'
    call allocate_array(iqys,maxnys,'iqys',ierr)
    call allocate_array(iqxs,maxnxs,'iqxs',ierr)
    call allocate_array(ctimes,1,maxnts+1,0,maxizs,'ctimes',ierr)
    call allocate_array(xouts,maxnxs,'xouts',ierr)
    call allocate_array(xs,maxnxs,'xs',ierr)
    call allocate_array(ys,maxnys,'ys',ierr)
    call allocate_array(qys,maxqys,'qys',ierr)
    call allocate_array(qxs,-maxqxs,'qxs',maxqxs,ierr)
    call allocate_array(xwids,maxnxs,'xwids',ierr)
    call allocate_array(ywids,maxnys,'ywids',ierr)
    call allocate_array(xcyls,maxnxs,'xcyls',ierr)
    call allocate_array(dwelts,0,'dwelts',maxizs,ierr)
    call allocate_array(dwelfs,maxnts,'dwelfs',ierr)
    call allocate_array(qs,-maxqxs,'qs',maxqxs,ierr)
    call allocate_array(youts,-maxnys,'youts',maxnys,ierr)
    call allocate_array(delps,maxnxs,maxnys,'delps',ierr)

  end subroutine allocate_mod_comxyt
  
  
  subroutine deallocate_mod_comxyt
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(ps)
    deallocate(pwids)
    deallocate(pzone)
    deallocate(pbin_bnds)
    deallocate(surf_bnds)

    if (allocated(iqys)) deallocate(iqys)
    if (allocated(iqxs)) deallocate(iqxs)
    if (allocated(ctimes)) deallocate(ctimes)
    if (allocated(xouts)) deallocate(xouts)
    if (allocated(xs)) deallocate(xs)
    if (allocated(ys)) deallocate(ys)
    if (allocated(qys)) deallocate(qys)
    if (allocated(qxs)) deallocate(qxs)
    if (allocated(xwids)) deallocate(xwids)
    if (allocated(ywids)) deallocate(ywids)
    if (allocated(xcyls)) deallocate(xcyls)
    if (allocated(dwelts)) deallocate(dwelts)
    if (allocated(dwelfs)) deallocate(dwelfs)
    if (allocated(qs)) deallocate(qs)
    if (allocated(youts)) deallocate(youts)
    if (allocated(delps)) deallocate(delps)

  end subroutine deallocate_mod_comxyt
  
  
  
end module mod_comxyt
