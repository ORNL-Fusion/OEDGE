module mod_div6
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /div6c/dact,dparas,dsum4,dsum1,dsum2,dsum3,diz,douts,cistiz,cistfp,cist,&
  !     cistot,cismax,cist_recstart,cist_elapsed,coreouts,incore, inedge, inmsol, indiv,&
  !      intrap,tmpncore,tmpiz,num_entered_core,debug0,debug_all,sdvs2,sdvs3,sdvb,sdtimp
  
  ! save /div6c/
  double precision,public :: dact,dsum4
  double precision,public,allocatable :: dparas(:,:)
  !
  !     for some reason - the count of particle time steps was recorded
  !     as double precision - however - the variable cist was the
  !     one used everywhere in the code calculations - this doesn't
  !     make sense - if a dp is required for accuracy then it would
  !     have to be used everywhere for it to make sense - so for now
  !     i am using the single precision cist everywhere.
  !
  !      double precision dista
  !
  !      double precision distacc
  !      real             cistacc,cistlast
  !      logical          cisterr
  !      integer          cisterrcnt
  !
  !     divimp timing variables which record particle timesteps
  !     - switch to double precision for all ion time step accumulation to
  !       enable precise timing for long lived particles
  !     - for a typical case with a time step on the order of 1.0e-7 seconds
  !       a real variable can accurately accumulate just about 1s of
  !       particle lifetime - occasionally this has been found to be
  !       insufficient. since neutrals exist for a shorter length of time
  !       typically - this switch to double precision has not yet been made
  !       for the neutral transport modules.
  !
  !
  !      real      cistiz,cistfp,cist,cistot,cismax,cist_recstart,cist_elapsed
  !
  double precision,public :: dsum1,dsum2,dsum3,diz
  double precision,public,allocatable :: douts(:,:)
  !
  real*8,public :: cistiz,cistfp,cist,cistot,cismax,cist_recstart,cist_elapsed
  !
  !     logical variables for detecting regions
  !
  double precision,public,allocatable :: coreouts(:,:)
  logical,public :: incore,inedge,inmsol,indiv,intrap
  real,public :: tmpncore,tmpiz
  !
  !     debugging flag: hard coded
  !
  real,public :: num_entered_core
  !
  !     array to record ion velocities. (riv)
  !
  !     moved to (dynam3) real   sdvs(maxnks,maxnrs,-1:maxizs)
  !
  logical,public :: debug0,debug_all
  real,public,allocatable :: sdvs2(:,:,:)
  real,public,allocatable :: sdvs3(:,:,:,:)
  real,public,allocatable :: sdvb(:,:)
  real,public,allocatable :: sdtimp(:,:,:)

  public :: allocate_mod_div6,deallocate_mod_div6

contains

  subroutine allocate_mod_div6
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_div6','ALLOCATE')

    call allocate_array(dparas,maxizs,6,'dparas',ierr)
    call allocate_array(douts,maxizs,9,'douts',ierr)
    call allocate_array(coreouts,maxizs,9,'coreouts',ierr)
    call allocate_array(sdvs2,1,maxnks,1,maxnrs,-1,maxizs,'sdvs2',ierr)
    call allocate_array(sdvs3,1,maxnks,1,maxnrs,-1,maxizs,1,2,'sdvs3',ierr)
    call allocate_array(sdvb,maxnks,maxnrs,'sdvb',ierr)
    call allocate_array(sdtimp,1,maxnks,1,maxnrs,-1,maxizs,'sdtimp',ierr)

  end subroutine allocate_mod_div6


  subroutine deallocate_mod_div6
    implicit none

    call pr_trace('mod_div6','DEALLOCATE')

    if (allocated(dparas)) deallocate(dparas)
    if (allocated(douts)) deallocate(douts)
    if (allocated(coreouts)) deallocate(coreouts)
    if (allocated(sdvs2)) deallocate(sdvs2)
    if (allocated(sdvs3)) deallocate(sdvs3)
    if (allocated(sdvb)) deallocate(sdvb)
    if (allocated(sdtimp)) deallocate(sdtimp)

  end subroutine deallocate_mod_div6

end module mod_div6