module mod_line_profile
  use debug_options
  implicit none

  !
  !     line_profile:
  !
  !     this file contains the common block with the parameters
  !     related to calculating a velocity shifted line profile
  !     function for a specific spectral line - data for which is
  !     loaded from adas.
  !
  !     possible future enhancements:
  !     - calculating for more than one line at a time
  !     - calculating for ion line profiles as well as
  !       neutral ones.
  !
  !     -*-fortran-*-
  integer,public :: max_lp_bins
  
  parameter (max_lp_bins=51)
  ! common /line_profile_data/ line_profile,modified_line_profile,line_profile_opt,&
  !     lp_bin_width,lp_instrument_width,lp_wave,lp_robs,lp_zobs,lp_theta,lp_dtheta,lp_adasyr,&
  !     lp_isele,lp_iselr,lp_iselx,lp_iseld,lp_adasid,lp_adasex
  !
  !     base option
  !
  ! save /line_profile_data/
  !
  !     bin and instrument widths, wavelength
  !
  integer,public :: line_profile_opt
  real,public :: lp_bin_width,lp_instrument_width
  !
  !     instrument view
  !
  real,public :: lp_wave
  !
  !     profile data
  !
  real,public :: lp_robs,lp_zobs,lp_theta,lp_dtheta
  !
  !     adas data - specifiers and selectors
  !
  real*8,public,allocatable :: line_profile(:),modified_line_profile(:)
  integer,public :: lp_adasyr
  integer,public :: lp_isele,lp_iselr,lp_iselx,lp_iseld
  character*80,public :: lp_adasid
  
  !
  character*3,public :: lp_adasex

  public :: allocate_mod_line_profile,deallocate_mod_line_profile

contains

  subroutine allocate_mod_line_profile
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_line_profile','ALLOCATE')

    call allocate_array(line_profile,-max_lp_bins,'line_profile',max_lp_bins,ierr)
    call allocate_array(modified_line_profile,-max_lp_bins,'modified_line_profile',max_lp_bins,&
         ierr)

  end subroutine allocate_mod_line_profile


  subroutine deallocate_mod_line_profile
    implicit none

    call pr_trace('mod_line_profile','DEALLOCATE')

    if (allocated(line_profile)) deallocate(line_profile)
    if (allocated(modified_line_profile)) deallocate(modified_line_profile)

  end subroutine deallocate_mod_line_profile

end module mod_line_profile