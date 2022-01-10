module mod_div7
  use debug_options
  implicit none

  
  !     -*-fortran-*-
  ! common /div7c/veltot,dsparanorm,dvparanorm,dsparastep,dvparastep,vparastep,vparanorm,&
  !     dsparacnt,dvparacnt,dvmaxv,dvminv,sigma_vel
  
  ! save /div7c/
  !
  !     variables to record parallel diffusive steps
  !
  !     krieger, ipp 12/94
  !
  real,public,allocatable :: veltot(:)
  double precision,public,allocatable :: dsparanorm(:),dvparanorm(:,:)
  double precision,public,allocatable :: dsparastep(:),dvparastep(:,:)
  double precision,public,allocatable :: vparastep(:,:)
  double precision,public,allocatable :: vparanorm(:,:)
  double precision,public,allocatable :: dsparacnt(:),dvparacnt(:,:)
  ! geier ipp/01
  double precision,public,allocatable :: dvmaxv(:,:),dvminv(:,:)
  
  real,public :: sigma_vel

  public :: allocate_mod_div7,deallocate_mod_div7

contains

  subroutine allocate_mod_div7
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div7','ALLOCATE')

    call allocate_array(veltot,maxizs,'veltot',ierr)
    call allocate_array(dsparanorm,6,'dsparanorm',ierr)
    call allocate_array(dvparanorm,6,maxizs+1,'dvparanorm',ierr)
    call allocate_array(dsparastep,6,'dsparastep',ierr)
    call allocate_array(dvparastep,6,maxizs+1,'dvparastep',ierr)
    call allocate_array(vparastep,6,maxizs+1,'vparastep',ierr)
    call allocate_array(vparanorm,6,maxizs+1,'vparanorm',ierr)
    call allocate_array(dsparacnt,6,'dsparacnt',ierr)
    call allocate_array(dvparacnt,6,maxizs+1,'dvparacnt',ierr)
    call allocate_array(dvmaxv,6,maxizs+1,'dvmaxv',ierr)
    call allocate_array(dvminv,6,maxizs+1,'dvminv',ierr)

  end subroutine allocate_mod_div7


  subroutine deallocate_mod_div7
    implicit none

    call pr_trace('mod_div7','DEALLOCATE')

    if (allocated(veltot)) deallocate(veltot)
    if (allocated(dsparanorm)) deallocate(dsparanorm)
    if (allocated(dvparanorm)) deallocate(dvparanorm)
    if (allocated(dsparastep)) deallocate(dsparastep)
    if (allocated(dvparastep)) deallocate(dvparastep)
    if (allocated(vparastep)) deallocate(vparastep)
    if (allocated(vparanorm)) deallocate(vparanorm)
    if (allocated(dsparacnt)) deallocate(dsparacnt)
    if (allocated(dvparacnt)) deallocate(dvparacnt)
    if (allocated(dvmaxv)) deallocate(dvmaxv)
    if (allocated(dvminv)) deallocate(dvminv)

  end subroutine deallocate_mod_div7

end module mod_div7