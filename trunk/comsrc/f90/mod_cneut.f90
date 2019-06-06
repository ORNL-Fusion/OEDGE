module mod_cneut
  use debug_options
  implicit none

  !
  !
  !     -*-fortran-*-
  ! common /cneut/  sputys,xatizs,yatizs,xprods,vins  ,snews ,yprods,rmaxs ,ranva ,&
  !     ranvb ,ranvc ,katizs,temtizs,satizs,cistizs,idprods,isprods,idatizs,mtcinf,mtctotcnt,&
  !     recinf,rectotcnt,eprods
  !
  ! save /cneut/
  !
  !      integer idprods(maximp),isprods(maximp),idatizs(maximp,2)
  !
  ! ammod begin.
  real,public,allocatable :: sputys(:),xatizs(:),yatizs(:),katizs(:),vins(:),xprods(:),&
       yprods(:),temtizs(:),rmaxs(:),ranva(:),ranvb(:),satizs(:),ranvc(:),snews(:),&
       cistizs(:),mtcinf(:,:),mtctotcnt(:,:),recinf(:,:),rectotcnt(:),eprods(:)

  ! moved from comtor
  real,public,allocatable :: cleakpos(:,:),launchdat(:,:)

  ! added idatizs(imp,3) -> irstart, idatizs(imp,4) -> ikstart
  ! for hydrocarbon launches that reduce to c+.
  !
  !     jdemod - bug - travel_locations needs to be in a common
  !              block so it can be properly accessible in multiple
  !              routines - at least ptravel_locations and div - this
  !              declaration makes the variable local
  !
  integer,public,allocatable :: idprods(:),isprods(:),idatizs(:,:)
  ! common /hc_neut/ travel_locations
  ! travel_locations contains incore,inedge,inmsol,indiv,intrap
  ! for hydrocarbon launches that reduce to c+.
  
  ! save /hc_neut/
  ! ammod end.
  logical,public,allocatable :: travel_locations(:,:)

  public :: allocate_mod_cneut,deallocate_mod_cneut

contains

  subroutine allocate_mod_cneut
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cneut','ALLOCATE')

    call allocate_array(sputys,maximp,'sputys',ierr)
    call allocate_array(xatizs,maximp,'xatizs',ierr)
    call allocate_array(yatizs,maximp,'yatizs',ierr)
    call allocate_array(katizs,maximp,'katizs',ierr)
    call allocate_array(vins,maximp,'vins',ierr)
    call allocate_array(xprods,maximp,'xprods',ierr)
    call allocate_array(yprods,maximp,'yprods',ierr)
    call allocate_array(temtizs,maximp,'temtizs',ierr)
    call allocate_array(rmaxs,maximp,'rmaxs',ierr)
    call allocate_array(ranva,maximp,'ranva',ierr)
    call allocate_array(ranvb,maximp,'ranvb',ierr)
    call allocate_array(ranvc,maximp,'ranvc',ierr)
    call allocate_array(satizs,maximp,'satizs',ierr)
    call allocate_array(snews,maximp,'snews',ierr)
    call allocate_array(cistizs,maximp,'cistizs',ierr)
    call allocate_array(mtcinf,7,3,'mtcinf',ierr)
    call allocate_array(mtctotcnt,0,11,1,3,'mtctotcnt',ierr)
    call allocate_array(recinf,14,7,'recinf',ierr)
    call allocate_array(rectotcnt,0,'rectotcnt',11,ierr)
    call allocate_array(eprods,maximp,'eprods',ierr)
    call allocate_array(idprods,maximp,'idprods',ierr)
    call allocate_array(isprods,maximp,'isprods',ierr)
    call allocate_array(idatizs,maximp,4,'idatizs',ierr)
    call allocate_array(travel_locations,maximp,5,'travel_locations',ierr)

    call allocate_array(cleakpos,maximp,2,'cleakpos',ierr)
    call allocate_array(launchdat,maximp,5,'launchdat',ierr)
    
  end subroutine allocate_mod_cneut


  subroutine deallocate_mod_cneut
    implicit none

    call pr_trace('mod_cneut','DEALLOCATE')

    if (allocated(sputys)) deallocate(sputys)
    if (allocated(xatizs)) deallocate(xatizs)
    if (allocated(yatizs)) deallocate(yatizs)
    if (allocated(katizs)) deallocate(katizs)
    if (allocated(vins)) deallocate(vins)
    if (allocated(xprods)) deallocate(xprods)
    if (allocated(yprods)) deallocate(yprods)
    if (allocated(temtizs)) deallocate(temtizs)
    if (allocated(rmaxs)) deallocate(rmaxs)
    if (allocated(ranva)) deallocate(ranva)
    if (allocated(ranvb)) deallocate(ranvb)
    if (allocated(satizs)) deallocate(satizs)
    if (allocated(ranvc)) deallocate(ranvc)
    if (allocated(snews)) deallocate(snews)
    if (allocated(cistizs)) deallocate(cistizs)
    if (allocated(mtcinf)) deallocate(mtcinf)
    if (allocated(mtctotcnt)) deallocate(mtctotcnt)
    if (allocated(recinf)) deallocate(recinf)
    if (allocated(rectotcnt)) deallocate(rectotcnt)
    if (allocated(eprods)) deallocate(eprods)
    if (allocated(idprods)) deallocate(idprods)
    if (allocated(isprods)) deallocate(isprods)
    if (allocated(idatizs)) deallocate(idatizs)
    if (allocated(travel_locations)) deallocate(travel_locations)

    if (allocated(cleakpos)) deallocate(cleakpos)
    if (allocated(launchdat)) deallocate(launchdat)

  end subroutine deallocate_mod_cneut

end module mod_cneut
