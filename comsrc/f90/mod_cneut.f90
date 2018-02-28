module mod_cneut
  use debug_options

  implicit none

  private

  !     -*-Fortran-*-                                                                       
  !
  !
  !      COMMON /CNEUT/  SPUTYS,XATIZS,YATIZS,XPRODS,VINS  ,SNEWS ,        
  !     >                YPRODS,RMAXS ,RANVA ,RANVB ,RANVC ,               
  !     >                KATIZS,TEMTIZS,SATIZS,cistizs,
  !     >                IDPRODS,ISPRODS,IDATIZS,
  !     >                mtcinf,mtctotcnt,recinf,rectotcnt,
  !     >                EPRODS
  !      save /cneut/
  !
  !      REAL SPUTYS(MAXIMP),XATIZS(MAXIMP),YATIZS(MAXIMP),KATIZS(MAXIMP), 
  !     >  VINS(MAXIMP),XPRODS(MAXIMP),YPRODS(MAXIMP),TEMTIZS(MAXIMP),
  !     >  RMAXS(MAXIMP),RANVA(MAXIMP),RANVB(MAXIMP),SATIZS(MAXIMP),
  !     >  RANVC(MAXIMP),SNEWS(MAXIMP),cistizs(maximp),
  !     >  mtcinf(7,3),mtctotcnt(0:11,3),recinf(14,7),rectotcnt(0:11),
  !     >  eprods(maximp)
  !                                                                       
  !      INTEGER IDPRODS(MAXIMP),ISPRODS(MAXIMP),IDATIZS(MAXIMP,2)
  !                                                                       
  ! ammod begin.
  ! Added IDATIZS(IMP,3) -> irstart, IDATIZS(IMP,4) -> ikstart
  ! for hydrocarbon launches that reduce to C+.
  !      INTEGER IDPRODS(MAXIMP),ISPRODS(MAXIMP),IDATIZS(MAXIMP,4)
  !
  !     jdemod - bug - travel_locations needs to be in a common 
  !              block so it can be properly accessible in multiple
  !              routines - at least ptravel_locations and div - this 
  !              declaration makes the variable local
  !
  !      common /hc_neut/ travel_locations
  !     ! Travel_Locations contains incore,inedge,inmsol,indiv,intrap
  !      ! for hydrocarbon launches that reduce to C+.
  !      save /hc_neut/
  !
  !      Logical Travel_Locations (MAXIMP,5)
  ! ammod end.


  REAL,allocatable, public ::  SPUTYS(:),XATIZS(:),YATIZS(:),KATIZS(:),& 
       VINS(:),XPRODS(:),YPRODS(:),TEMTIZS(:),&
       RMAXS(:),RANVA(:),RANVB(:),SATIZS(:),&
       RANVC(:),SNEWS(:),cistizs(:),&
       mtcinf(:,:),mtctotcnt(:,:),recinf(:,:),rectotcnt(:),&
       eprods(:)


  INTEGER,allocatable,public ::  IDPRODS(:),ISPRODS(:),IDATIZS(:,:)

  Logical, allocatable,public ::  Travel_Locations (:,:)

  public :: allocate_cneut, deallocate_cneut





contains

  subroutine allocate_cneut
    use global_parameters
    use allocate_arrays
    use error_handling
    implicit none

    integer :: ierr

    call pr_trace('MOD_CNEUT','ALLOCATE')

    call allocate_array(sputys,maximp,'SPUTYS',ierr)
    call allocate_array(snews,maximp,'SNEWS',ierr)

    call allocate_array(xatizs,maximp,'XATIZS',ierr)
    call allocate_array(yatizs,maximp,'YATIZS',ierr)
    call allocate_array(katizs,maximp,'KATIZS',ierr)
    call allocate_array(satizs,maximp,'SATIZS',ierr)
    call allocate_array(cistizs,maximp,'CISTIZS',ierr)

    call allocate_array(vins,maximp,'VINS',ierr)

    call allocate_array(xprods,maximp,'XPRODS',ierr)
    call allocate_array(yprods,maximp,'YPRODS',ierr)
    call allocate_array(eprods,maximp,'EPRODS',ierr)

    call allocate_array(temtizs,maximp,'TEMTIZS',ierr)
    call allocate_array(rmaxs,maximp,'RMAXS',ierr)

    call allocate_array(ranva,maximp,'RANVA',ierr)
    call allocate_array(ranvb,maximp,'RANVB',ierr)
    call allocate_array(ranvc,maximp,'RANVC',ierr)

    call allocate_array(mtcinf,7,3,'MTCINF',ierr)
    call allocate_array(mtctotcnt,0,11,1,3,'MTCTOTCNT',ierr)

    call allocate_array(recinf,14,7,'RECINF',ierr)

    ! The generalized allocate_arrays does not support a single dimension
    ! array with two bounds due to the signature conflict with a 2D array
    ! So it is manually allocated here
    allocate(rectotcnt(0:11),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array '//trim('RECTOTCNT')//' IERR =',ierr)
    endif

    call allocate_array(idprods,maximp,'IDPRODS',ierr)
    call allocate_array(isprods,maximp,'ISPRODS',ierr)
    call allocate_array(idatizs,maximp,4,'IDATIZS',ierr)

    ! no generic procedure for logical arrays ... for now allocate manually

    if (allocated(travel_locations)) deallocate(travel_locations)
    allocate(travel_locations(maximp,5),stat=ierr)
    if (ierr.ne.0) then 
       call errmsg('Error allocating array travel_locations IERR =',ierr)
    endif


  end subroutine allocate_cneut




  subroutine deallocate_cneut
    implicit none

    call pr_trace('MOD_CNEUT','DEALLOCATE')

    if (allocated(sputys))  deallocate(sputys)
    if (allocated(snews))  deallocate(snews)

    if (allocated(xatizs))  deallocate(xatizs)
    if (allocated(yatizs))  deallocate(yatizs)
    if (allocated(katizs))  deallocate(katizs)
    if (allocated(satizs))  deallocate(satizs)
    if (allocated(cistizs))  deallocate(cistizs)

    if (allocated(vins))  deallocate(vins)

    if (allocated(xprods))  deallocate(xprods)
    if (allocated(yprods))  deallocate(yprods)
    if (allocated(eprods))  deallocate(eprods)

    if (allocated(temtizs))  deallocate(temtizs)
    if (allocated(rmaxs))  deallocate(rmaxs)

    if (allocated(ranva))  deallocate(ranva)
    if (allocated(ranvb))  deallocate(ranvb)
    if (allocated(ranvc))  deallocate(ranvc)

    if (allocated(mtcinf))  deallocate(mtcinf)
    if (allocated(mtctotcnt))  deallocate(mtctotcnt)

    if (allocated(recinf))  deallocate(recinf)

    if (allocated(rectotcnt)) deallocate(rectotcnt)

    if (allocated(idprods)) deallocate(idprods)
    if (allocated(isprods)) deallocate(isprods)
    if (allocated(idatizs)) deallocate(idatizs)
    if (allocated(travel_locations)) deallocate(travel_locations)
    

  end subroutine deallocate_cneut



end module mod_cneut
