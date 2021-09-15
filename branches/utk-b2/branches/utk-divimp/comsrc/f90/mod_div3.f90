module mod_div3
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  !     >  imp_drftvel,bg_drftvel,pol_drftv,
  !     >  sdrft_start,sdrft_end,
  ! common /div3c/xtripp,xtrips,oldz,tntots,stmp,kprob,fptarg,fptart,fpttot,acttarg,&
  !     fpdist,fplosstim,fpent,fpexit,diffr,rfptarg,tautim,wssf,theta1,pintim,rc,reccnt,&
  !     ikorg,irorg,mtcrecstruk,mtcrecwalln,vrec,recloss,dist_travelled
  !
  ! save /div3c/
  real,public :: xtripp,xtrips,oldz
  real,public,allocatable :: tntots(:,:)
  real,public :: stmp,kprob
  real,public :: fptart,fpttot
  real,public,allocatable :: fptarg(:),acttarg(:)
  real,public :: fpdist,fplosstim
  real,public :: fpent,fpexit,diffr,rfptarg
  !
  !     poloidal drift variables
  !
  !      real      imp_drftvel,bg_drftvel,pol_drftv,
  !     >          sdrft_start,sdrft_end
  !
  real,public :: tautim
  real,public :: wssf
  real,public :: theta1
  
  !
  !     variables related to following recombined impurities
  !
  real,public :: pintim
  integer,public :: rc,reccnt,ikorg,irorg
  real,public :: mtcrecstruk,mtcrecwalln,vrec
  real,public :: recloss
  real,public :: dist_travelled

  public :: allocate_mod_div3,deallocate_mod_div3

contains

  subroutine allocate_mod_div3
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div3','ALLOCATE')

    call allocate_array(tntots,0,maxizs+1,1,4,'tntots',ierr)
    call allocate_array(fptarg,maxizs,'fptarg',ierr)
    call allocate_array(acttarg,maxizs,'acttarg',ierr)

  end subroutine allocate_mod_div3


  subroutine deallocate_mod_div3
    implicit none

    call pr_trace('mod_div3','DEALLOCATE')

    if (allocated(tntots)) deallocate(tntots)
    if (allocated(fptarg)) deallocate(fptarg)
    if (allocated(acttarg)) deallocate(acttarg)

  end subroutine deallocate_mod_div3

end module mod_div3