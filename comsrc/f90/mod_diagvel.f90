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

  ! for vtig calculations
  real :: integration_const,vtig_crmi,vtig_crmb
  real,public,allocatable :: vtig_array(:,:)   ! holds the calculated values of vtig

  public :: allocate_mod_diagvel,deallocate_mod_diagvel,allocate_debugv_mod_diagvel
  public :: setup_vtig, calc_vtig_array
  
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




  subroutine setup_vtig(crmb,crmi)
    use mod_params
    implicit none
    real :: crmi,crmb

    ! assume ln(lambda) ~= 15.0
    ! note - density is also needed
    ! note - Bi ~= 2.6 Z^2 has been used to simplify the expression which is correct as U -> 1  (U = mz/(mi+mz))
    ! mi = crmb, mz = crmi
    ! note mz in Tau_s cancels the mz in the pre-factor for the vTiG expression

    ! integration constant = e Tau_s Bi/mz Ti^3/2     (  vTiG(s) = e Tau_s Bi/mz  dT/ds)
    ! Tau_s = (1.47e13 mz Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni Z^2 ln_Lam )
    ! mz = cmri, mi = crmb
    ! Bi = 2.6 Z^2
    !
    ! integration_constant = e/(amu mz)  (1.47e13 mz Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni Z^2 ln_Lam ) * 2.6 Z^2
    ! cancel mz term and Z^2 term 
    ! integration_constant = e/amu  (1.47e13 Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni ln_Lam ) * 2.6
    ! pull out Ti^3/2
    ! integration_constant = e/amu  (1.47e13 (1/mi)^(1/2)) / ( (1+mi/mz) ni ln_Lam ) * 2.6    Ti^3/2
    ! pull out n since it isn't known until later
    ! integration_constant = e/amu  (1.47e13 (1/mi)^(1/2)) / ( (1+mi/mz) ln_Lam ) * 2.6    /n     Ti^3/2
    ! integration_constant = e/amu  (1.47e13 (1/crmb)^(1/2)) / ( (1+crmb/crmi) ln_Lam ) * 2.6 
       
    real,parameter :: ln_lambda = 15.0

    integration_const = ech/amu * 1.47e13 * sqrt(1.0/crmb) * 2.6 / ( ( 1.0 + crmb/crmi) * ln_lambda)
    vtig_crmi = crmi
    vtig_crmb = crmb
    
  end subroutine setup_vtig

  subroutine calc_vtig_array(qtim)
    use mod_params
    use mod_cgeom
    use allocate_arrays
    implicit none
    real :: qtim
    integer :: ik,ir,ierr

    ! allocate and calculate the vtig velocity array
    ! use temporary storage
    ierr = 0
    call allocate_array(vtig_array,1,maxnks,1,maxnrs,'vtig_array',ierr)
    
    do ir = 1,nrs
       do ik = 1,nks(ir)
          vtig_array(ik,ir) = calc_vtig(ik,ir,qtim)
       end do
    end do

  end subroutine calc_vtig_array
  
  real function calc_vtig(ik,ir,qtim)
    ! DIVIMP version
    !use mod_io_units
    !use mod_vtig
    !use mod_comtor
    !use mod_comxyt
    !use mod_comt2
    use mod_cgeom

    implicit none
    integer :: ik,ir
    real :: tgscal,qtim


    TGSCAL = (1.6E-19)/(VTIG_CRMI*1.673E-27) * QTIM *QTIM 

    ! print out
    if (knbs(ik,ir).gt.0.0) then 
       calc_vtig = integration_const/knbs(ik,ir) * ktibs(ik,ir)**(1.5) * kfigs(ik,ir) / tgscal
    else
       calc_vtig = 0.0
    endif

  end function calc_vtig

  
    
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
