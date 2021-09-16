module mod_div2
  use debug_options
  implicit none

  !
  !
  !     -*-fortran-*-
  ! common /div2c/rstruk,tstruk,fact,ran,emax,ff,fe,feg,fig,fvel,fvh,fvg,stopped_follow,&
  !     iprob,mtcstruk,tmtcstruk,mtcwalln,tmtcwalln,ratio1,ratio2,tau,temold,ran1,ran2,&
  !     rgauss,k11,k12,k13,d11,d12,d13,chipara,lambda2,xkpara,xdparapara,ftotal,coption,&
  !     rions,stots,xxx,sss,sptots,sitots,e2dtots,e2dptots,impurity_content,tstepl,rconst,&
  !     kkk,ssss,sdtzs,sdtzs2
  
  
  ! save /div2c/
  !
  real,public :: rstruk,tstruk,fact,ran,emax,ff,fe,feg,fig,fvel,fvh,fvg
  !
  real,public :: stopped_follow
  !
  !     momentum  transfer collision variables for neut
  !
  real*8,public :: iprob
  !
  real,public :: mtcstruk,tmtcstruk,mtcwalln,tmtcwalln
  !
  ! psmod
  !
  real,public :: ratio1,ratio2,tau,temold,ran1,ran2,rgauss
  real,public :: k11,k12,k13,d11,d12,d13
  !
  !
  !     forceplot
  !
  real,public :: chipara,lambda2,xkpara,xdparapara
  real,public :: ftotal
  
  
  !
  !
  ! psmod
  !
  integer,public :: coption
  real,public,allocatable :: rions(:),stots(:),xxx(:),sss(:)
  real,public,allocatable :: sptots(:,:),sitots(:,:)
  real,public,allocatable :: e2dtots(:),e2dptots(:,:)
  !
  real,public,allocatable :: impurity_content(:,:,:)
  real,public :: tstepl,rconst,kkk
  real,public,allocatable :: ssss(:)
  real,public,allocatable :: sdtzs(:),sdtzs2(:)

  public :: allocate_mod_div2,deallocate_mod_div2

contains

  subroutine allocate_mod_div2
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div2','ALLOCATE')

    call allocate_array(rions,maxizs,'rions',ierr)
    call allocate_array(stots,46,'stots',ierr)
    call allocate_array(xxx,maxizs,'xxx',ierr)
    call allocate_array(sss,maxizs,'sss',ierr)
    call allocate_array(sptots,1,10,0,maxizs+2,'sptots',ierr)
    call allocate_array(sitots,maxizs+2,2,'sitots',ierr)
    call allocate_array(e2dtots,9,'e2dtots',ierr)
    call allocate_array(e2dptots,1,6,0,maxe2dizs,'e2dptots',ierr)
    call allocate_array(impurity_content,0,maxizs,1,4,1,2,'impurity_content',ierr)
    call allocate_array(ssss,maxizs,'ssss',ierr)
    call allocate_array(sdtzs,-1,'sdtzs',maxizs,ierr)
    call allocate_array(sdtzs2,-1,'sdtzs2',maxizs,ierr)

  end subroutine allocate_mod_div2


  subroutine deallocate_mod_div2
    implicit none

    call pr_trace('mod_div2','DEALLOCATE')

    if (allocated(rions)) deallocate(rions)
    if (allocated(stots)) deallocate(stots)
    if (allocated(xxx)) deallocate(xxx)
    if (allocated(sss)) deallocate(sss)
    if (allocated(sptots)) deallocate(sptots)
    if (allocated(sitots)) deallocate(sitots)
    if (allocated(e2dtots)) deallocate(e2dtots)
    if (allocated(e2dptots)) deallocate(e2dptots)
    if (allocated(impurity_content)) deallocate(impurity_content)
    if (allocated(ssss)) deallocate(ssss)
    if (allocated(sdtzs)) deallocate(sdtzs)
    if (allocated(sdtzs2)) deallocate(sdtzs2)

  end subroutine deallocate_mod_div2

end module mod_div2
