module mod_div6

  use debug_options

  implicit none
  private


  real*8, public :: DACT,DSUM4
  real*8, public :: DSUM1,DSUM2,DSUM3,DIZ

  real*8, allocatable, public :: DPARAS(:,:),DOUTS(:,:)


  !     For some reason - the count of particle time steps was recorded
  !     as double precision - however - the variable cist was the 
  !     one used everywhere in the code calculations - this doesn't 
  !     make sense - if a DP is required for accuracy then it would
  !     have to be used everywhere for it to make sense - so for now
  !     I am using the single precision cist everywhere. 

  !      double precision dista

  !      double precision distacc
  !      real             cistacc,cistlast
  !      logical          cisterr
  !      integer          cisterrcnt 

  !     DIVIMP timing variables which record particle timesteps 
  !     - switch to double precision for all ion time step accumulation to 
  !       enable precise timing for long lived particles
  !     - for a typical case with a time step on the order of 1.0e-7 seconds
  !       a real variable can accurately accumulate just about 1s of 
  !       particle lifetime - occasionally this has been found to be 
  !       insufficient. Since neutrals exist for a shorter length of time
  !       typically - this switch to double precision has not yet been made
  !       for the neutral transport modules. 

  !      real      cistiz,cistfp,cist,cistot,cismax,cist_recstart,cist_elapsed

  real*8,public ::  cistiz,cistfp,cist,cistot,cismax,&
       cist_recstart,cist_elapsed

  real*8, allocatable,public ::  coreouts(:,:) 

  !     Logical variables for detecting regions

  logical,public :: incore, inedge, inmsol, indiv, intrap  
  real,public ::  tmpncore,tmpiz
  real, public ::  num_entered_core

  !     DEBUGGING FLAG: Hard coded

  logical,public ::  debug0,debug_all

  !     Array to record ion velocities. (RIV)

  !     moved to (dynam3) real   sdvs(maxnks,maxnrs,-1:maxizs)

  real, allocatable, public :: sdvs2(:,:,:),&
       sdvs3(:,:,:,:), sdvb (:,:),&
       sdtimp(:,:,:)


  public :: allocate_div6, deallocate_div6


contains



  subroutine allocate_div6
    use global_parameters
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('MOD_DIV6','ALLOCATE')

    call allocate_array(dparas,maxizs,6,'DPARAS',ierr)
    call allocate_array(douts,maxizs,9,'DOUTS',ierr)

    call allocate_array(coreouts,maxizs,9,'COREOUTS',ierr)

    call allocate_array(sdvs2,1,maxnks,1,maxnrs,-1,maxizs,'SDVS2',ierr)
    call allocate_array(sdvs3,1,maxnks,1,maxnrs,-1,maxizs,1,2,'SDVS3',ierr)
    call allocate_array(sdvb,maxnks,maxnrs,'SDVB',ierr)
    call allocate_array(sdtimp,1,maxnks,1,maxnrs,-1,maxizs,'SDTIMP',ierr)


  end subroutine allocate_div6


  subroutine deallocate_div6
    implicit none

    call pr_trace('MOD_DIV6','DEALLOCATE')

    if (allocated(dparas)) deallocate(dparas)
    if (allocated(douts)) deallocate(douts)

    if (allocated(coreouts)) deallocate(coreouts)

    if (allocated(sdvs2)) deallocate(sdvs2)
    if (allocated(sdvs3)) deallocate(sdvs3)
    if (allocated(sdvb)) deallocate(sdvb)
    if (allocated(sdtimp)) deallocate(sdtimp)

  end subroutine deallocate_div6


end module mod_div6
