module mod_walls_com

  use debug_options

  implicit none

  private


!c     -*-Fortran-*-
!c
!c     This common block contains data related to the wall 
!c     - including geometry and plasma conditions if any
!c
!c
!      integer, parameter :: nwall_data = 32
!
!      common /wall_data/ wall_plasma_opt,wall_plasma_fact,
!     >    WLPROB,NWLPROB,wlpabs,
!    >    WLWALL1,WLWALL2,WLTRAP1,WLTRAP2,
!     >    wlwall3,wlwall4,wltrap3,wltrap4,
!     >    WLIND,NWLIND,FWLPROB,TOTWL,
!     >    WALLPT,WALLPTS,
!     >    nimindex,wallindex,
!     >    wallsrc,wallleak,
!     >    pcnt,rw,zw
!      save /wall_data/
!c
!      integer wall_plasma_opt,
!     >    NWLPROB,WALLPTS,
!     >    WLWALL1,WLWALL2,WLTRAP1,WLTRAP2,
!     >    wlwall3,wlwall4,wltrap3,wltrap4,
!     >    NWLIND,WLIND(MAXPTS),wlpabs,
!     >    nimindex(maxnds),wallindex(maxnds),
!     >    pcnt
!c
!      real wall_plasma_fact,
!     >    WALLPT(MAXPTS,nwall_data),
!     >    FWLPROB(MAXPTS),TOTWL,
!     >    WLPROB(MAXPTS,3),
!     >    wallsrc(5,3),wallleak(5,3),
!     >    rw(maxpts),zw(maxpts)
!c
  !
  !     This common block contains data related to the wall 
  !     - including geometry and plasma conditions if any
  !
  !
  integer, parameter,public :: nwall_data = 32

  integer,public :: wall_plasma_opt,&
       NWLPROB,WALLPTS,WLWALL1,WLWALL2,WLTRAP1,WLTRAP2,&
       wlwall3,wlwall4,wltrap3,wltrap4,NWLIND,wlpabs,pcnt

  integer, allocatable, public :: WLIND(:),nimindex(:),wallindex(:)


  real,public ::  wall_plasma_fact,TOTWL,wallsrc(5,3),wallleak(5,3)


  real, allocatable, public :: WALLPT(:,:),FWLPROB(:),WLPROB(:,:),rw(:),zw(:)


  public :: allocate_walls_com,deallocate_walls_com


contains

  subroutine allocate_walls_com
    use global_parameters
    use allocate_arrays
    integer :: ierr

    call allocate_array(wlind,maxpts,'WLIND',ierr)
    call allocate_array(nimindex,maxnds,'NIMINDEX',ierr)
    call allocate_array(wallindex,maxnds,'WALLINDEX',ierr)

    call allocate_array(wallpt,maxpts,nwall_data,'WALLPT',ierr)
    call allocate_array(fwlprob,maxpts,'FWLPROB',ierr)
    call allocate_array(wlprob,maxpts,3,'WLPROB',ierr)
    call allocate_array(rw,maxpts,'RW',ierr)
    call allocate_array(zw,maxpts,'ZW',ierr)


  end subroutine allocate_walls_com


  subroutine deallocate_walls_com
    implicit none

    deallocate(wlind)
    deallocate(nimindex)
    deallocate(wallindex)

    deallocate(wallpt)
    deallocate(fwlprob)
    deallocate(wlprob)
    deallocate(rw)
    deallocate(zw)

  end subroutine deallocate_walls_com


end module mod_walls_com
