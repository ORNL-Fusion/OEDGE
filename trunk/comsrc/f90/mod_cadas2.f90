module mod_cadas2
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  !
  ! common /cadas2/ dtev, ddens, dtevd, ddensd, drcofd, zdata,drcofi, titlf,ifail, ievcut,&
  !      itmaxd, idmaxd, izmaxd
  !
  ! save /cadas2/

  integer, public, parameter :: maxads = 100

  real*8,public,allocatable :: dtev(:),ddens(:),dtevd(:),ddensd(:),drcofd(:,:,:),zdata(:),&
       drcofi(:)
  character,public :: titlf*80
  integer,public :: ifail,ievcut,itmaxd,idmaxd,izmaxd

  public :: allocate_mod_cadas2,deallocate_mod_cadas2

contains

  subroutine allocate_mod_cadas2(maxpts)
    !use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr,maxpts

    call pr_trace('mod_cadas2','ALLOCATE')

    call allocate_array(dtev,maxpts,'dtev',ierr)
    call allocate_array(ddens,maxpts,'ddens',ierr)
    call allocate_array(dtevd,maxads,'dtevd',ierr)
    call allocate_array(ddensd,maxads,'ddensd',ierr)
    call allocate_array(drcofd,maxads,maxads,maxads,'drcofd',ierr)
    call allocate_array(zdata,maxads,'zdata',ierr)
    call allocate_array(drcofi,maxpts,'drcofi',ierr)

  end subroutine allocate_mod_cadas2


  subroutine deallocate_mod_cadas2
    implicit none

    call pr_trace('mod_cadas2','DEALLOCATE')

    if (allocated(dtev)) deallocate(dtev)
    if (allocated(ddens)) deallocate(ddens)
    if (allocated(dtevd)) deallocate(dtevd)
    if (allocated(ddensd)) deallocate(ddensd)
    if (allocated(drcofd)) deallocate(drcofd)
    if (allocated(zdata)) deallocate(zdata)
    if (allocated(drcofi)) deallocate(drcofi)

  end subroutine deallocate_mod_cadas2

end module mod_cadas2
