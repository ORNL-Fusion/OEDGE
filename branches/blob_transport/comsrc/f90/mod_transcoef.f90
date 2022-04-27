module mod_transcoef
  use debug_options
  implicit none

  !
  !     this file contains the declarations for the transport
  !     coefficient data - including ring coordinates at the
  !     inner and outer midplanes.
  !
  !     -*-fortran-*-
  !     >     rcouter,rcinner
  ! common /transcoef/ inmid,oumid,dperp,odperp,idperp,chiperpi,ochiperpi,ichiperpi,&
  !     chiperpe,ochiperpe,ichiperpe,xperpt,ixperpt,oxperpt,gradte,gradti,gradn,e2dgradte,&
  !     e2dgradti,e2dgradn
  !
  ! save /transcoef/
  !     >     rcouter(maxnrs),rcinner(maxnrs)
  real,public,allocatable :: dperp(:),odperp(:),idperp(:),chiperpi(:),ochiperpi(:),&
       ichiperpi(:),chiperpe(:),ochiperpe(:),ichiperpe(:),gradte(:,:),gradti(:,:),gradn(:,:),&
       xperpt(:),ixperpt(:),oxperpt(:),e2dgradte(:,:),e2dgradti(:,:),e2dgradn(:,:)
  !
  integer,public :: inmid,oumid

  public :: allocate_mod_transcoef,deallocate_mod_transcoef

contains

  subroutine allocate_mod_transcoef
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_transcoef','ALLOCATE')

    call allocate_array(dperp,maxnrs,'dperp',ierr)
    call allocate_array(odperp,maxnrs,'odperp',ierr)
    call allocate_array(idperp,maxnrs,'idperp',ierr)
    call allocate_array(chiperpi,maxnrs,'chiperpi',ierr)
    call allocate_array(ochiperpi,maxnrs,'ochiperpi',ierr)
    call allocate_array(ichiperpi,maxnrs,'ichiperpi',ierr)
    call allocate_array(chiperpe,maxnrs,'chiperpe',ierr)
    call allocate_array(ochiperpe,maxnrs,'ochiperpe',ierr)
    call allocate_array(ichiperpe,maxnrs,'ichiperpe',ierr)
    call allocate_array(gradte,maxnks,maxnrs,'gradte',ierr)
    call allocate_array(gradti,maxnks,maxnrs,'gradti',ierr)
    call allocate_array(gradn,maxnks,maxnrs,'gradn',ierr)
    call allocate_array(xperpt,maxnrs,'xperpt',ierr)
    call allocate_array(ixperpt,maxnrs,'ixperpt',ierr)
    call allocate_array(oxperpt,maxnrs,'oxperpt',ierr)
    call allocate_array(e2dgradte,maxnks,maxnrs,'e2dgradte',ierr)
    call allocate_array(e2dgradti,maxnks,maxnrs,'e2dgradti',ierr)
    call allocate_array(e2dgradn,maxnks,maxnrs,'e2dgradn',ierr)

  end subroutine allocate_mod_transcoef


  subroutine deallocate_mod_transcoef
    implicit none

    call pr_trace('mod_transcoef','DEALLOCATE')

    if (allocated(dperp)) deallocate(dperp)
    if (allocated(odperp)) deallocate(odperp)
    if (allocated(idperp)) deallocate(idperp)
    if (allocated(chiperpi)) deallocate(chiperpi)
    if (allocated(ochiperpi)) deallocate(ochiperpi)
    if (allocated(ichiperpi)) deallocate(ichiperpi)
    if (allocated(chiperpe)) deallocate(chiperpe)
    if (allocated(ochiperpe)) deallocate(ochiperpe)
    if (allocated(ichiperpe)) deallocate(ichiperpe)
    if (allocated(gradte)) deallocate(gradte)
    if (allocated(gradti)) deallocate(gradti)
    if (allocated(gradn)) deallocate(gradn)
    if (allocated(xperpt)) deallocate(xperpt)
    if (allocated(ixperpt)) deallocate(ixperpt)
    if (allocated(oxperpt)) deallocate(oxperpt)
    if (allocated(e2dgradte)) deallocate(e2dgradte)
    if (allocated(e2dgradti)) deallocate(e2dgradti)
    if (allocated(e2dgradn)) deallocate(e2dgradn)

  end subroutine deallocate_mod_transcoef

end module mod_transcoef