module mod_parmmod
  use debug_options
  implicit none

  
  !     -*-fortran-*-
  integer,public :: nreac
  
  parameter (nreac=50)
  double precision,public :: eps60
  
  
  parameter       (eps60=1.0d-60)
  ! common /amjuel/iswr,modclf,creac,fparm,rcmx,ifexmx,rcmn,ifexmn
  ! save /amjuel/
  integer,public,allocatable :: iswr(:),modclf(:),ifexmx(:,:),ifexmn(:,:)
  double precision,public,allocatable :: creac(:,:,:),fparm(:,:,:),rcmx(:,:),rcmn(:,:)

  public :: allocate_mod_parmmod,deallocate_mod_parmmod

contains

  subroutine allocate_mod_parmmod
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_parmmod','ALLOCATE')

    call allocate_array(iswr,nreac,'iswr',ierr)
    call allocate_array(modclf,nreac,'modclf',ierr)
    call allocate_array(ifexmx,nreac,2,'ifexmx',ierr)
    call allocate_array(ifexmn,nreac,2,'ifexmn',ierr)
    call allocate_array(creac,1,9,0,9,-10,nreac,'creac',ierr)
    call allocate_array(fparm,nreac,6,2,'fparm',ierr)
    call allocate_array(rcmx,nreac,2,'rcmx',ierr)
    call allocate_array(rcmn,nreac,2,'rcmn',ierr)

  end subroutine allocate_mod_parmmod


  subroutine deallocate_mod_parmmod
    implicit none

    call pr_trace('mod_parmmod','DEALLOCATE')

    if (allocated(iswr)) deallocate(iswr)
    if (allocated(modclf)) deallocate(modclf)
    if (allocated(ifexmx)) deallocate(ifexmx)
    if (allocated(ifexmn)) deallocate(ifexmn)
    if (allocated(creac)) deallocate(creac)
    if (allocated(fparm)) deallocate(fparm)
    if (allocated(rcmx)) deallocate(rcmx)
    if (allocated(rcmn)) deallocate(rcmn)

  end subroutine deallocate_mod_parmmod

end module mod_parmmod