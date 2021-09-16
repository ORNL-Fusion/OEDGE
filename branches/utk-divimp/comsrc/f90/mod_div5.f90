module mod_div5
  use debug_options
  implicit none

  !
  !
  !     -*-fortran-*-
  ! slmod begin - t-dep
  !
  !     >  maxciz,jz,ix,iy,in,ikv,
  ! slmod end
  ! common /div5c/dsputy,dtots,dqtim,ditots,ptots,procterm,lpinz0,best,dsq,rsq,zsq,&
  !     res,kklim,kk,natiz,nprod,status,iw,id,jd,nimps_local,nimps2_local,minciz,maxciz,&
  !     jz,ix,iy,in,ikv,ikstart,irstart,idstart,idtype,iwstart,sstart,imp,mattar,matp,im,&
  !     it,jk,j,m,jr,ik1,ir1,diffus,debug,inmain
  
  ! save /div5c/
  !
  logical,public :: procterm
  !
  logical,public :: lpinz0
  !
  real,public :: best,dsq,rsq,zsq
  !
  !
  integer,public :: res
  integer,public :: kklim,kk,natiz,nprod,status,iw,id,jd
  ! slmod begin - t-dep
  integer,public :: nimps_local,nimps2_local
  !
  !      integer   maxciz,jz,ix,iy,in,ikv
  ! slmod end
  integer,public :: minciz,maxciz,jz,ix,iy,in,ikv
  integer,public :: ikstart,irstart,idstart,idtype,iwstart
  real,public :: sstart
  !
  integer,public :: imp,mattar,matp,im,it,jk,j,m
  !
  !     these are found in the printopt common block
  !
  !      character*5 inner,outer
  !
  integer,public :: jr,ik1,ir1
  logical,public :: diffus,debug,inmain
  double precision,public :: dsputy,dqtim
  double precision,public,allocatable :: dtots(:)
  double precision,public,allocatable :: ditots(:,:)
  double precision,public,allocatable :: ptots(:,:)

  public :: allocate_mod_div5,deallocate_mod_div5

contains

  subroutine allocate_mod_div5
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div5','ALLOCATE')

    call allocate_array(dtots,46,'dtots',ierr)
    call allocate_array(ditots,maxizs+2,2,'ditots',ierr)
    call allocate_array(ptots,1,10,0,maxizs+2,'ptots',ierr)

  end subroutine allocate_mod_div5


  subroutine deallocate_mod_div5
    implicit none

    call pr_trace('mod_div5','DEALLOCATE')

    if (allocated(dtots)) deallocate(dtots)
    if (allocated(ditots)) deallocate(ditots)
    if (allocated(ptots)) deallocate(ptots)

  end subroutine deallocate_mod_div5

end module mod_div5