module mod_comsol
  use debug_options
  implicit none

  !
  !     common blocks for source functions
  !
  !     -*-fortran-*-
  ! common /setinfop/ plensrc,plamsrc,p0out,p0in,radint,radinti
  ! save /setinfop/
  !
  double precision,public :: plensrc,plamsrc,p0out,p0in,radint,radinti
  ! common /setinfo/ lensrc,lamsrc,facsrc,s0out,s0in,fiz,s0ain,s0aout,s0bin,s0bout,&
  !     ionint,ioninti
  ! save /setinfo/
  !
  double precision,public :: lensrc,lamsrc,facsrc,s0out,s0in,fiz,s0ain,s0aout,s0bin,&
       s0bout,ionint,ioninti
  ! common /srcinfo/ fnorm,fnormi,irn,tss,maxik,tmaxs
  ! save /srcinfo/
  integer,public :: irn,maxik
  real,public :: tmaxs
  real,public,allocatable :: tss(:)
  !
  double precision,public :: fnorm,fnormi
  ! common /lastinfo/ sionl,sradl,iionl,iradl,ipeil,speil,pionl,pradl,ppeil
  ! save /lastinfo/
  double precision,public :: sionl,sradl,iionl,iradl,speil,ipeil
  !
  integer,public :: pionl,pradl,ppeil

  public :: allocate_mod_comsol,deallocate_mod_comsol

contains

  subroutine allocate_mod_comsol
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_comsol','ALLOCATE')

    call allocate_array(tss,maxnks,'tss',ierr)

  end subroutine allocate_mod_comsol


  subroutine deallocate_mod_comsol
    implicit none

    call pr_trace('mod_comsol','DEALLOCATE')

    if (allocated(tss)) deallocate(tss)

  end subroutine deallocate_mod_comsol

end module mod_comsol