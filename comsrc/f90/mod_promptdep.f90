module mod_promptdep
  use debug_options
  implicit none

  !
  !     common block containing variables related to the prompt
  !     deposition option.
  !
  !     -*-fortran-*-
  ! common /prompt_data/ prompt_depopt,promptdeps,mps_thickness,mps_energy
  !
  ! save /prompt_data/
  !
  real,public,allocatable :: promptdeps(:,:),mps_thickness(:,:),mps_energy(:,:)
  !
  integer,public :: prompt_depopt
  
  ! sazmod - adding generalized prompt redeposition coefficients.
  real, public :: prompt_dep_a, prompt_dep_b

  public :: allocate_mod_promptdep,deallocate_mod_promptdep

contains

  subroutine allocate_mod_promptdep
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_promptdep','ALLOCATE')

    call allocate_array(promptdeps,maxnds,9,'promptdeps',ierr)
    call allocate_array(mps_thickness,maxnrs,2,'mps_thickness',ierr)
    call allocate_array(mps_energy,maxnrs,2,'mps_energy',ierr)

  end subroutine allocate_mod_promptdep


  subroutine deallocate_mod_promptdep
    implicit none

    call pr_trace('mod_promptdep','DEALLOCATE')

    if (allocated(promptdeps)) deallocate(promptdeps)
    if (allocated(mps_thickness)) deallocate(mps_thickness)
    if (allocated(mps_energy)) deallocate(mps_energy)

  end subroutine deallocate_mod_promptdep

end module mod_promptdep
