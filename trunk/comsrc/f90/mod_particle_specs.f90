module mod_particle_specs
  use debug_options
  implicit none

  !
  !
  !     -*-fortran-*-
  ! common /ion_specs/ ik,ir,iz,s,theta,cross,r,z,ik_last,ir_last,s_last,cross_last,&
  !     theta_last,r_last,z_last,oldtheta,oldcross,ikold,irold,ifate,vel,temi
  !
  ! save /ion_specs/
  integer,public :: ik,ir,iz,ik_last,ir_last
  integer,public :: ikold,irold
  integer,public :: ifate
  real,public :: s,theta,cross,r,z
  real,public :: s_last,theta_last,cross_last,r_last,z_last
  real,public :: oldtheta,oldcross
  !
  real,public :: vel,temi

  public :: allocate_mod_particle_specs,deallocate_mod_particle_specs

contains

  subroutine allocate_mod_particle_specs
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_particle_specs','ALLOCATE')


  end subroutine allocate_mod_particle_specs


  subroutine deallocate_mod_particle_specs
    implicit none

    call pr_trace('mod_particle_specs','DEALLOCATE')


  end subroutine deallocate_mod_particle_specs

end module mod_particle_specs