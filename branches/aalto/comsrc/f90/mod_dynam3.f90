module mod_dynam3
  use debug_options
  implicit none

  !
  !     ipp 01/10 - added array wallsiz for charge state resolved
  !     particle exit statistics
  !     jdemod - added wallseiz to get exit energy
  !     add wallsil - to record deposition of divertor leaked particles
  !
  !     -*-fortran-*-
  ! common /dynam3/ tizs,powls,lines,zeffs,walls,deps,neros,elims,wallsn,hpowls,hlines,&
  !     sdvs,chemden,chemizs,ncore, nedge, ntrap,ndivert,nmsol,wallse,wallse_i,wallsi,&
  !     wallsiz,wallseiz,wallsil
  ! save /dynam3/
  real,public,allocatable :: tizs(:,:,:),zeffs(:,:,:),powls(:,:,:),lines(:,:,:),walls(:,:,:),&
       deps(:,:),neros(:,:),elims(:,:,:),wallsn(:),sdvs(:,:,:),chemden(:,:),chemizs(:,:),&
       hpowls(:,:,:),hlines(:,:,:),ncore(:,:),nedge(:,:),ntrap(:,:),ndivert(:,:),&
       nmsol(:,:),wallse(:),wallse_i(:),wallsi(:),wallsiz(:,:),wallseiz(:,:),wallsil(:)

  public :: allocate_mod_dynam3,deallocate_mod_dynam3

contains

  subroutine allocate_mod_dynam3
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam3','ALLOCATE')

    call allocate_array(tizs,1,maxnks,1,maxnrs,-1,maxizs,'tizs',ierr)
    call allocate_array(zeffs,maxnks,maxnrs,3,'zeffs',ierr)
    call allocate_array(powls,1,maxnks,1,maxnrs,-1,maxizs,'powls',ierr)
    call allocate_array(lines,1,maxnks,1,maxnrs,-1,maxizs,'lines',ierr)
    call allocate_array(walls,1,maxnks,1,maxnrs,0,maxizs+1,'walls',ierr)
    call allocate_array(deps,maxnds,maxizs,'deps',ierr)
    call allocate_array(neros,maxnds,5,'neros',ierr)
    call allocate_array(elims,1,maxnks,1,3,-1,maxizs,'elims',ierr)
    call allocate_array(wallsn,maxpts+1,'wallsn',ierr)
    call allocate_array(sdvs,1,maxnks,1,maxnrs,-1,maxizs,'sdvs',ierr)
    call allocate_array(chemden,maxnks,maxnrs,'chemden',ierr)
    call allocate_array(chemizs,maxnks,maxnrs,'chemizs',ierr)
    call allocate_array(hpowls,1,maxnks,1,maxnrs,0,1,'hpowls',ierr)
    call allocate_array(hlines,1,maxnks,1,maxnrs,0,1,'hlines',ierr)
    call allocate_array(ncore,maxnks,maxnrs,'ncore',ierr)
    call allocate_array(nedge,maxnks,maxnrs,'nedge',ierr)
    call allocate_array(ntrap,maxnks,maxnrs,'ntrap',ierr)
    call allocate_array(ndivert,maxnks,maxnrs,'ndivert',ierr)
    call allocate_array(nmsol,maxnks,maxnrs,'nmsol',ierr)
    call allocate_array(wallse,maxpts+1,'wallse',ierr)
    call allocate_array(wallse_i,maxpts+1,'wallse_i',ierr)
    call allocate_array(wallsi,maxpts+1,'wallsi',ierr)
    call allocate_array(wallsiz,maxpts+1,maxizs + 1,'wallsiz',ierr)
    call allocate_array(wallseiz,maxpts+1,maxizs + 1,'wallseiz',ierr)
    call allocate_array(wallsil,maxpts+1,'wallsil',ierr)

  end subroutine allocate_mod_dynam3


  subroutine deallocate_mod_dynam3
    implicit none

    call pr_trace('mod_dynam3','DEALLOCATE')

    if (allocated(tizs)) deallocate(tizs)
    if (allocated(zeffs)) deallocate(zeffs)
    if (allocated(powls)) deallocate(powls)
    if (allocated(lines)) deallocate(lines)
    if (allocated(walls)) deallocate(walls)
    if (allocated(deps)) deallocate(deps)
    if (allocated(neros)) deallocate(neros)
    if (allocated(elims)) deallocate(elims)
    if (allocated(wallsn)) deallocate(wallsn)
    if (allocated(sdvs)) deallocate(sdvs)
    if (allocated(chemden)) deallocate(chemden)
    if (allocated(chemizs)) deallocate(chemizs)
    if (allocated(hpowls)) deallocate(hpowls)
    if (allocated(hlines)) deallocate(hlines)
    if (allocated(ncore)) deallocate(ncore)
    if (allocated(nedge)) deallocate(nedge)
    if (allocated(ntrap)) deallocate(ntrap)
    if (allocated(ndivert)) deallocate(ndivert)
    if (allocated(nmsol)) deallocate(nmsol)
    if (allocated(wallse)) deallocate(wallse)
    if (allocated(wallse_i)) deallocate(wallse_i)
    if (allocated(wallsi)) deallocate(wallsi)
    if (allocated(wallsiz)) deallocate(wallsiz)
    if (allocated(wallseiz)) deallocate(wallseiz)
    if (allocated(wallsil)) deallocate(wallsil)

  end subroutine deallocate_mod_dynam3

end module mod_dynam3