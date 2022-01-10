module mod_grbound
  use debug_options
  implicit none

  !
  !       this common block contains the information about
  !       the inner and outer bounds of the defined plasma
  !       grid region. this information is stored for re-use
  !       by the ga15a and ga15b subroutines.
  !
  !     -*-fortran-*-
  ! common /grbound/ ionwpts,ioncpts,neutwpts,riw,ziw,rcw,zcw,rnw,znw,iwindw,iwwork,&
  !     iwtdum,iwxdum,iwydum,icindw,icwork,ictdum,icxdum,icydum,nwindw,nwwork,nwtdum,nwxdum,&
  !     nwydum,lastin,outgrid
  !
  ! save /grbound/
  !
  real,public,allocatable :: riw(:),ziw(:),rcw(:),zcw(:),rnw(:),znw(:),iwwork(:),iwtdum(:),&
       iwxdum(:),iwydum(:),icwork(:),ictdum(:),icxdum(:),icydum(:),nwwork(:),nwtdum(:),&
       nwxdum(:),nwydum(:)
  !
  integer,public :: ionwpts,ioncpts,neutwpts,lastin
  integer,public,allocatable :: iwindw(:,:),icindw(:,:),nwindw(:,:)
  logical,public :: outgrid
  
  
  
  

  public :: allocate_mod_grbound,deallocate_mod_grbound

contains

  subroutine allocate_mod_grbound
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_grbound','ALLOCATE')

    call allocate_array(riw,maxpts,'riw',ierr)
    call allocate_array(ziw,maxpts,'ziw',ierr)
    call allocate_array(rcw,maxpts,'rcw',ierr)
    call allocate_array(zcw,maxpts,'zcw',ierr)
    call allocate_array(rnw,maxpts,'rnw',ierr)
    call allocate_array(znw,maxpts,'znw',ierr)
    call allocate_array(iwwork,4*maxpts,'iwwork',ierr)
    call allocate_array(iwtdum,maxpts,'iwtdum',ierr)
    call allocate_array(iwxdum,maxpts,'iwxdum',ierr)
    call allocate_array(iwydum,maxpts,'iwydum',ierr)
    call allocate_array(icwork,4*maxpts,'icwork',ierr)
    call allocate_array(ictdum,maxpts,'ictdum',ierr)
    call allocate_array(icxdum,maxpts,'icxdum',ierr)
    call allocate_array(icydum,maxpts,'icydum',ierr)
    call allocate_array(nwwork,4*maxpts,'nwwork',ierr)
    call allocate_array(nwtdum,maxpts,'nwtdum',ierr)
    call allocate_array(nwxdum,maxpts,'nwxdum',ierr)
    call allocate_array(nwydum,maxpts,'nwydum',ierr)
    call allocate_array(iwindw,2,maxpts,'iwindw',ierr)
    call allocate_array(icindw,2,maxpts,'icindw',ierr)
    call allocate_array(nwindw,2,maxpts,'nwindw',ierr)

  end subroutine allocate_mod_grbound


  subroutine deallocate_mod_grbound
    implicit none

    call pr_trace('mod_grbound','DEALLOCATE')

    if (allocated(riw)) deallocate(riw)
    if (allocated(ziw)) deallocate(ziw)
    if (allocated(rcw)) deallocate(rcw)
    if (allocated(zcw)) deallocate(zcw)
    if (allocated(rnw)) deallocate(rnw)
    if (allocated(znw)) deallocate(znw)
    if (allocated(iwwork)) deallocate(iwwork)
    if (allocated(iwtdum)) deallocate(iwtdum)
    if (allocated(iwxdum)) deallocate(iwxdum)
    if (allocated(iwydum)) deallocate(iwydum)
    if (allocated(icwork)) deallocate(icwork)
    if (allocated(ictdum)) deallocate(ictdum)
    if (allocated(icxdum)) deallocate(icxdum)
    if (allocated(icydum)) deallocate(icydum)
    if (allocated(nwwork)) deallocate(nwwork)
    if (allocated(nwtdum)) deallocate(nwtdum)
    if (allocated(nwxdum)) deallocate(nwxdum)
    if (allocated(nwydum)) deallocate(nwydum)
    if (allocated(iwindw)) deallocate(iwindw)
    if (allocated(icindw)) deallocate(icindw)
    if (allocated(nwindw)) deallocate(nwindw)

  end subroutine deallocate_mod_grbound

end module mod_grbound