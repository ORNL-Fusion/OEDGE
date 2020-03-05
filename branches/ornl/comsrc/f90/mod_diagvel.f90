module mod_diagvel
  use mod_params
  use debug_options
  implicit none

  !
  !     this file contains the declarations for debugging the
  !     velocity distribution. when this option is not needed
  !     these variables can be minimized so that memory is not
  !     used.
  !
  !     -*-fortran-*-
  ! jdemod - move maxvizs to mod_params
  integer,public :: nvel,nvtime,maxvnks
  !integer,public :: nvel,nvtime,maxvizs,maxvnks
  real,public :: velsep
  parameter (nvel=50,nvtime=10,velsep=0.1)
  !
  ! jdemod - move maxvizs to mod_params
  parameter (maxvnks=maxnks)
  !parameter (maxvizs=maxizs,maxvnks=maxnks)
  ! common /diagvel/ velplate,debugv,cstepv,velspace,velweight,velts,veltsw,ts,velcoord,&
  !     velcell
  !
  ! save /diagvel/
  real,public :: velplate
  real,public,allocatable :: velspace(:,:,:)
  real,public,allocatable :: velweight(:,:,:)
  !
  !     the following 2 are really only used as local variables at the moment
  !
  real,public,allocatable :: velts(:,:),veltsw(:,:),ts(:)
  real,public,allocatable :: velcoord(:)
  !
  real,public,allocatable :: velcell(:)
  !
  integer,public :: cstepv
  !
  !     vtime is a local variable but data values are initialized
  !
  logical,public :: debugv
  real,public,allocatable :: vtime(:)
  integer,public :: ti_calc_opt 

  
  !data vtime /0.1, 0.2, 0.4, 0.6, 0.8, 1.0,1.5, 2.0, 3.0, 5.0/
  
  !
  ! jdemod - moved from mod_div6.f90
  !
  real*8,public,allocatable :: ddvs2(:,:,:)
  real*8,public,allocatable :: ddvs3(:,:,:,:)
  !
  real,public,allocatable :: sdvb(:,:)
  real,public,allocatable :: sdtimp(:,:,:)
  
  

  public :: allocate_mod_diagvel,deallocate_mod_diagvel,allocate_debugv_mod_diagvel

contains

  subroutine allocate_mod_diagvel
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_diagvel','ALLOCATE')

    call allocate_array(velspace,-nvel,nvel+1,1,maxvizs,1,maxvnks,'velspace',ierr)
    call allocate_array(velweight,-nvel,nvel+1,1,maxizs,1,maxvnks,'velweight',ierr)
    call allocate_array(velts,nvel,nvtime,'velts',ierr)
    call allocate_array(veltsw,nvel,nvtime,'veltsw',ierr)
    call allocate_array(ts,nvtime,'ts',ierr)
    call allocate_array(velcoord,-nvel,'velcoord',nvel+1,ierr)
    call allocate_array(velcell,maxnks,'velcell',ierr)
    call allocate_array(vtime,nvtime,'vtime',ierr)

    ! assign initial values to vtime
    vtime = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0,1.5, 2.0, 3.0, 5.0]

    call allocate_debugv_mod_diagvel
    
  end subroutine allocate_mod_diagvel


    
  subroutine allocate_debugv_mod_diagvel
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('debugv_mod_diagvel','ALLOCATE')
    ! jdemod
    ! only allocate if velocity debugging is turned on
    ! all references to these variables should be inside
    ! if debugv statements
    ! the input file has been read before this code is executed
    if (debugv) then
       call allocate_array(ddvs2,1,maxnks,1,maxnrs,-1,maxizs,'ddvs2',ierr)
       call allocate_array(ddvs3,1,maxnks,1,maxnrs,-1,maxizs,1,2,'ddvs3',ierr)
       call allocate_array(sdvb,maxnks,maxnrs,'sdvb',ierr)
       call allocate_array(sdtimp,1,maxnks,1,maxnrs,-1,maxizs,'sdtimp',ierr)
    endif

  end subroutine allocate_debugv_mod_diagvel


  
  subroutine deallocate_mod_diagvel
    implicit none

    call pr_trace('mod_diagvel','DEALLOCATE')

    if (allocated(velspace)) deallocate(velspace)
    if (allocated(velweight)) deallocate(velweight)
    if (allocated(velts)) deallocate(velts)
    if (allocated(veltsw)) deallocate(veltsw)
    if (allocated(ts)) deallocate(ts)
    if (allocated(velcoord)) deallocate(velcoord)
    if (allocated(velcell)) deallocate(velcell)
    if (allocated(vtime)) deallocate(vtime)
    
    if (allocated(ddvs2)) deallocate(ddvs2)
    if (allocated(ddvs3)) deallocate(ddvs3)
    if (allocated(sdvb)) deallocate(sdvb)
    if (allocated(sdtimp)) deallocate(sdtimp)

  end subroutine deallocate_mod_diagvel

end module mod_diagvel
