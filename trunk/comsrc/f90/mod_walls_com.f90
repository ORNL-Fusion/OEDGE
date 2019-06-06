module mod_walls_com
  use debug_options
  implicit none

  !
  !     this common block contains data related to the wall
  !     - including geometry and plasma conditions if any
  !
  !
  !     -*-fortran-*-
  
  integer, parameter ,public :: nwall_data = 32
  ! common /wall_data/ wall_plasma_opt,wall_plasma_fact,wlprob,nwlprob,wlpabs,wlwall1,&
  !     wlwall2,wltrap1,wltrap2,wlwall3,wlwall4,wltrap3,wltrap4,wlind,nwlind,fwlprob,&
  !     totwl,wallpt,wallpts,nimindex,wallindex,wallsrc,wallleak,pcnt,rw,zw
  !
  ! save /wall_data/
  !
  integer,public :: wall_plasma_opt,nwlprob,wallpts,wlwall1,wlwall2,wltrap1,wltrap2,&
       wlwall3,wlwall4,wltrap3,wltrap4,nwlind,wlpabs,pcnt
  integer,public,allocatable :: wlind(:),nimindex(:),wallindex(:)
  !
  real,public :: wall_plasma_fact,totwl
  real,public,allocatable :: wallpt(:,:),fwlprob(:),wlprob(:,:),wallsrc(:,:),wallleak(:,:),&
       rw(:),zw(:)

  public :: allocate_mod_walls_com,deallocate_mod_walls_com,allocate_mod_walls_com_input

contains

  subroutine allocate_mod_walls_com
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_walls_com','ALLOCATE')

    call allocate_array(wlind,maxpts,'wlind',ierr)
    call allocate_array(nimindex,maxnds,'nimindex',ierr)
    call allocate_array(wallindex,maxnds,'wallindex',ierr)
    call allocate_array(wallpt,maxpts,nwall_data,'wallpt',ierr)
    call allocate_array(fwlprob,maxpts,'fwlprob',ierr)
    call allocate_array(wallsrc,5,3,'wallsrc',ierr)
    call allocate_array(wallleak,5,3,'wallleak',ierr)
    call allocate_array(rw,maxpts,'rw',ierr)
    call allocate_array(zw,maxpts,'zw',ierr)

  end subroutine allocate_mod_walls_com


  subroutine deallocate_mod_walls_com
    implicit none

    call pr_trace('mod_walls_com','DEALLOCATE')

    if (allocated(wlind)) deallocate(wlind)
    if (allocated(nimindex)) deallocate(nimindex)
    if (allocated(wallindex)) deallocate(wallindex)
    if (allocated(wallpt)) deallocate(wallpt)
    if (allocated(fwlprob)) deallocate(fwlprob)
    if (allocated(wlprob)) deallocate(wlprob)
    if (allocated(wallsrc)) deallocate(wallsrc)
    if (allocated(wallleak)) deallocate(wallleak)
    if (allocated(rw)) deallocate(rw)
    if (allocated(zw)) deallocate(zw)

  end subroutine deallocate_mod_walls_com

  subroutine allocate_mod_walls_com_input
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_walls_com','ALLOCATE INPUT')

    call allocate_array(wlprob,maxpts,3,'wlprob',ierr)

  end subroutine allocate_mod_walls_com_input


end module mod_walls_com
