module mod_dynam4
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /dynam4/ imode,nts,cstmax,dwelts,dwelfs,ctimes,lims,walks
  ! save /dynam4/
  integer,public :: imode,nts
  real*8,public :: cstmax
  real,public,allocatable :: dwelts(:),dwelfs(:),lims(:,:,:,:),walks(:,:)
  ! jdemod - all timing variables should be double precision since single precision can't handle the possible
  !          number of time steps for long lived ions
  ! ctimes may require conversion to real*8 but should still work as real for most cases since it is just used
  ! for statistical diagnostics and recording density over time. 
  real,public,allocatable :: ctimes(:,:)
  
  public :: allocate_mod_dynam4,deallocate_mod_dynam4,allocate_mod_dynam4_input,allocate_mod_dynam4_input_special

contains

  subroutine allocate_mod_dynam4
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam4','ALLOCATE')

    call allocate_array(ctimes,0,maxnts+1,-1,maxizs,'ctimes',ierr)
    call allocate_array(lims,1,maxnks,1,maxnrs,-1,maxizs,1,maxnts,'lims',ierr)
    call allocate_array(walks,maxnws,2,'walks',ierr)

  end subroutine allocate_mod_dynam4


  subroutine deallocate_mod_dynam4
    implicit none

    call pr_trace('mod_dynam4','DEALLOCATE')

    if (allocated(dwelts)) deallocate(dwelts)
    if (allocated(dwelfs)) deallocate(dwelfs)
    if (allocated(ctimes)) deallocate(ctimes)
    if (allocated(lims)) deallocate(lims)
    if (allocated(walks)) deallocate(walks)

  end subroutine deallocate_mod_dynam4

  subroutine allocate_mod_dynam4_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam4','ALLOCATE INPUT')

    ! jdemod - dwelts is a problem - it is allocated with the default value of maxizs because it is
    !          needed in the input file but maxizs is also read in from the input file
    !        - allocation of this arrays is delayed until after maxizs is read in but before dwelts is
    !          read in
    !call allocate_array(dwelts,-1,'dwelts',maxizs,ierr)
    call allocate_array(dwelfs,maxnts,'dwelfs',ierr)

  end subroutine allocate_mod_dynam4_input

  subroutine allocate_mod_dynam4_input_special
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam4','ALLOCATE INPUT SPECIAL')

    ! jdemod - dwelts is a problem - it is allocated with the default value of maxizs because it is
    !          needed in the input file but 
    call allocate_array(dwelts,-1,'dwelts',maxizs,ierr)

  end subroutine allocate_mod_dynam4_input_special


  
end module mod_dynam4
