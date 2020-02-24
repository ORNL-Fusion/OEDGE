module mod_baffles
  use debug_options
  implicit none

  !
  !     baffles
  !
  !     this common block contains the declarations for the variables that
  !     contain information about any baffles that are defined in the input
  !     grid file. if appropriate and possible - these baffles are incorporated
  !     into the vessel wall for use in divimp.
  !
  !     the variable names were chosen to match those in use in edge2d/nimbus
  !
  !
  !     -*-fortran-*-
  ! common /baffles/ nbufle,rbufle,zbufle,nbufmx,nbufx,rbufx,zbufx,node_origin,nvesorg,&
  !     rvesorg,zvesorg,wallredef
  !
  ! save /baffles/
  integer,public :: nbufle,nbufmx,nvesorg
  integer,public,allocatable :: nbufx(:),node_origin(:,:)
  !
  integer,public :: wallredef
  real,public,allocatable :: rbufle(:),zbufle(:)
  real,public,allocatable :: rbufx(:,:),zbufx(:,:)
  
  
  
  
  real,public,allocatable :: rvesorg(:),zvesorg(:)

  public :: allocate_mod_baffles,deallocate_mod_baffles

contains

  subroutine allocate_mod_baffles
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_baffles','ALLOCATE')

    call allocate_array(nbufx,mbufx,'nbufx',ierr)
    call allocate_array(node_origin,mves,2,'node_origin',ierr)
    call allocate_array(rbufle,mbufle,'rbufle',ierr)
    call allocate_array(zbufle,mbufle,'zbufle',ierr)
    call allocate_array(rbufx,mbufx,mbufle,'rbufx',ierr)
    call allocate_array(zbufx,mbufx,mbufle,'zbufx',ierr)
    call allocate_array(rvesorg,mves,'rvesorg',ierr)
    call allocate_array(zvesorg,mves,'zvesorg',ierr)

  end subroutine allocate_mod_baffles


  subroutine deallocate_mod_baffles
    implicit none

    call pr_trace('mod_baffles','DEALLOCATE')

    if (allocated(nbufx)) deallocate(nbufx)
    if (allocated(node_origin)) deallocate(node_origin)
    if (allocated(rbufle)) deallocate(rbufle)
    if (allocated(zbufle)) deallocate(zbufle)
    if (allocated(rbufx)) deallocate(rbufx)
    if (allocated(zbufx)) deallocate(zbufx)
    if (allocated(rvesorg)) deallocate(rvesorg)
    if (allocated(zvesorg)) deallocate(zvesorg)

  end subroutine deallocate_mod_baffles

end module mod_baffles