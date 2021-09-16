module mod_temp
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /temp/ rate_calc
  !
  ! save /temp/
  !
  real,public,allocatable :: rate_calc(:,:,:)
  !
  logical,public :: check_rates = .true.
  
  !data check_rates /.true./

  public :: allocate_mod_temp,deallocate_mod_temp

contains

  subroutine allocate_mod_temp
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_temp','ALLOCATE')

    call allocate_array(rate_calc,maxnks,maxnrs,3,'rate_calc',ierr)

  end subroutine allocate_mod_temp


  subroutine deallocate_mod_temp
    implicit none

    call pr_trace('mod_temp','DEALLOCATE')

    if (allocated(rate_calc)) deallocate(rate_calc)

  end subroutine deallocate_mod_temp

end module mod_temp
