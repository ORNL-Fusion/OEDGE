module mod_psin_data
  use debug_options
  implicit none

  !
  !     this common block is used to hold externally read in psin data
  !     that is used for plotting within the out program
  !
  !     -*-fortran-*-
  ! common /psindata/ surface_data,rpsimin,r_int,z_int,psi_int,dist_int,n_int,n_elements,&
  !     external_psin,psin_filename
  !
  ! save /psindata/
  integer,public :: max_int
  !
  parameter (max_int=10)
  !
  real*8,public :: rpsimin
  real*8,public,allocatable :: surface_data(:,:),r_int(:),z_int(:),psi_int(:),dist_int(:)
  !
  integer,public :: n_int,n_elements
  !
  logical,public :: external_psin
  !
  character*200,public :: psin_filename

  public :: allocate_mod_psin_data,deallocate_mod_psin_data

contains

  subroutine allocate_mod_psin_data
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_psin_data','ALLOCATE')

    call allocate_array(surface_data,maxpts,4,'surface_data',ierr)
    call allocate_array(r_int,max_int,'r_int',ierr)
    call allocate_array(z_int,max_int,'z_int',ierr)
    call allocate_array(psi_int,max_int,'psi_int',ierr)
    call allocate_array(dist_int,max_int,'dist_int',ierr)

  end subroutine allocate_mod_psin_data


  subroutine deallocate_mod_psin_data
    implicit none

    call pr_trace('mod_psin_data','DEALLOCATE')

    if (allocated(surface_data)) deallocate(surface_data)
    if (allocated(r_int)) deallocate(r_int)
    if (allocated(z_int)) deallocate(z_int)
    if (allocated(psi_int)) deallocate(psi_int)
    if (allocated(dist_int)) deallocate(dist_int)

  end subroutine deallocate_mod_psin_data

end module mod_psin_data