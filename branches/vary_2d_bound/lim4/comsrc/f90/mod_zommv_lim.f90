module mod_zommv
  
  use mod_params
  
  !c
  !      common /zommv/  zicabs,zifabs,zilabs,zisabs,zrtabs,zrvabs,ztbs,
  !     >                zravav,zicay0,zicaaw,zica2l
  !      real            zicabs(maxizs), zifabs(maxizs), zilabs(maxizs)
  !      real            zisabs(maxizs), zicay0, zicaaw, zica2l
  !      real            zrtabs(maxizs), zrvabs(maxizs), zravav(maxizs)
  !      real            ztbs(maxizs)
  
  
  implicit none
  private
  
  real,public,allocatable:: zicabs(:),zifabs(:),zilabs(:)
  real,public:: zicay0,zicaaw,zica2l
  real,public,allocatable:: zisabs(:)
  real,public,allocatable:: zrtabs(:),zrvabs(:),zravav(:)
  real,public,allocatable:: ztbs(:)
  
  
  public :: allocate_mod_zommv, deallocate_mod_zommv
  
  
contains
  
  subroutine allocate_mod_zommv
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(dtev  ,maxnxs,'dtev',ierr)


    call allocate_array(zicabs,maxizs,'zicabs',ierr)
    call allocate_array(zifabs,maxizs,'zifabs',ierr)
    call allocate_array(zilabs,maxizs,'zilabs',ierr)
    call allocate_array(zisabs,maxizs,'zisabs',ierr)
    call allocate_array(zrtabs,maxizs,'zrtabs',ierr)
    call allocate_array(zrvabs,maxizs,'zrvabs',ierr)
    call allocate_array(zravav,maxizs,'zravav',ierr)
    call allocate_array(ztbs,maxizs,'ztbs',ierr)

  end subroutine allocate_mod_zommv
  
  
  subroutine deallocate_mod_zommv
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()

    if (allocated(zicabs)) deallocate(zicabs)
    if (allocated(zifabs)) deallocate(zifabs)
    if (allocated(zilabs)) deallocate(zilabs)
    if (allocated(zisabs)) deallocate(zisabs)
    if (allocated(zrtabs)) deallocate(zrtabs)
    if (allocated(zrvabs)) deallocate(zrvabs)
    if (allocated(zravav)) deallocate(zravav)
    if (allocated(ztbs)) deallocate(ztbs)

  end subroutine deallocate_mod_zommv
  
  
  
end module mod_zommv
