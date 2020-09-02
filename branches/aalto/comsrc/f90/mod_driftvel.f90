module mod_driftvel
  use debug_options
  implicit none

  !
  !     this common block holds the data related to the parallel
  !     drift velocity options.
  !
  !
  !     -*-fortran-*-
  !     >       imp_drftvel,bg_drftvel,
  ! common /driftvel_com/cpdrft,cdrftv,drft_region,cdrftv_start,cdrftv_end,pol_drftv,&
  !     sdrft_start,sdrft_end,ndrftvel,drftvel_machopt,drft_distopt,ringdrftvel,ringcs
  !
  !      real imp_drftvel,bg_drftvel,
  ! save /driftvel_com/
  real,public,allocatable :: pol_drftv(:)
  !
  real,public,allocatable :: sdrft_start(:),sdrft_end(:)
  !
  real,public :: cdrftv,cdrftv_start,cdrftv_end
  !
  integer,public :: cpdrft,drft_region,drft_distopt
  !
  integer,public :: ndrftvel,drftvel_machopt
  !
  real,public,allocatable :: ringdrftvel(:,:),ringcs(:)
  !     >      ,osmpotcell
  ! common /drifts/potopt,exb_rad_opt,exb_pol_opt,exb_scale,osmpot2,exb_rad_drft,exb_pol_drft,&
  !     e_pol,e_rad
  integer,public :: potopt,exb_rad_opt,exb_pol_opt
  !     >     ,osmpotcell(5,maxnks*maxnrs)
  real,public :: exb_scale
  real,public,allocatable :: osmpot2(:,:),exb_rad_drft(:,:),exb_pol_drft(:,:),e_pol(:,:),&
       e_rad(:,:)

  public :: allocate_mod_driftvel,deallocate_mod_driftvel,allocate_mod_driftvel_input

contains

  subroutine allocate_mod_driftvel
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_driftvel','ALLOCATE')

    call allocate_array(pol_drftv,maxnrs,'pol_drftv',ierr)
    call allocate_array(sdrft_start,maxnrs,'sdrft_start',ierr)
    call allocate_array(sdrft_end,maxnrs,'sdrft_end',ierr)
    call allocate_array(ringcs,maxnrs,'ringcs',ierr)
    call allocate_array(osmpot2,0,maxnks+1,1,maxnrs,'osmpot2',ierr)
    call allocate_array(exb_rad_drft,maxnks,maxnrs,'exb_rad_drft',ierr)
    call allocate_array(exb_pol_drft,maxnks,maxnrs,'exb_pol_drft',ierr)
    call allocate_array(e_pol,maxnks,maxnrs,'e_pol',ierr)
    call allocate_array(e_rad,maxnks,maxnrs,'e_rad',ierr)

  end subroutine allocate_mod_driftvel


  subroutine deallocate_mod_driftvel
    implicit none

    call pr_trace('mod_driftvel','DEALLOCATE')

    if (allocated(pol_drftv)) deallocate(pol_drftv)
    if (allocated(ringdrftvel)) deallocate(ringdrftvel)
    if (allocated(sdrft_start)) deallocate(sdrft_start)
    if (allocated(sdrft_end)) deallocate(sdrft_end)
    if (allocated(ringcs)) deallocate(ringcs)
    if (allocated(osmpot2)) deallocate(osmpot2)
    if (allocated(exb_rad_drft)) deallocate(exb_rad_drft)
    if (allocated(exb_pol_drft)) deallocate(exb_pol_drft)
    if (allocated(e_pol)) deallocate(e_pol)
    if (allocated(e_rad)) deallocate(e_rad)

  end subroutine deallocate_mod_driftvel

  subroutine allocate_mod_driftvel_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_driftvel','ALLOCATE INPUT')

    call allocate_array(ringdrftvel,maxnrs,2,'ringdrftvel',ierr)
    
  end subroutine allocate_mod_driftvel_input
  
end module mod_driftvel
