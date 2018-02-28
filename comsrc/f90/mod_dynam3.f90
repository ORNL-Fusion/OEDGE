module mod_dynam3

  use debug_options

  implicit none

  private

  !
  !        The purpose of this module is to allow dynamic allocation of DIVIMP arrays replacing the common blocks
  !        and static allocation originally used in the code. 
  !


  !         include 'dynam3'
  !                                                                       
  !     IPP 01/10 - added array wallsiz for charge state resolved
  !     particle exit statistics
  !     jdemod - added wallseiz to get exit energy
  !     add wallsil - to record deposition of divertor leaked particles
  !
  !      COMMON /DYNAM3/ TIZS,POWLS,LINES,ZEFFS,WALLS,DEPS,NEROS,    
  !     >   ELIMS,wallsn,hpowls,hlines,sdvs,chemden,chemizs,
  !     >   ncore, nedge, ntrap,ndivert,nmsol,wallse,wallse_i,wallsi,
  !     >   wallsiz,wallseiz,wallsil	 
  !      save /DYNAM3/
  !

  !      REAL TIZS(MAXNKS,MAXNRS,-1:MAXIZS),ZEFFS(MAXNKS,MAXNRS,3),        
  !     >   POWLS(MAXNKS,MAXNRS,-1:MAXIZS),LINES(MAXNKS,MAXNRS,-1:MAXIZS), 
  !     >   WALLS(MAXNKS,MAXNRS,0:MAXIZS+1),                               
  !     >   DEPS(MAXNDS,MAXIZS),NEROS(MAXNDS,5),
  !     >   ELIMS(MAXNKS,3,-1:MAXIZS),                                     
  !     >   wallsn(maxpts+1),sdvs(maxnks,maxnrs,-1:maxizs),
  !     >   chemden(maxnks,maxnrs),chemizs(maxnks,maxnrs),
  !     >   HPOWLS(MAXNKS,MAXNRS,0:1), HLINES(MAXNKS,MAXNRS,0:1),           
  !     >   ncore(maxnks,maxnrs),nedge(maxnks,maxnrs),
  !     >   ntrap(maxnks,maxnrs),ndivert(maxnks,maxnrs),
  !     >   nmsol(maxnks,maxnrs),wallse(maxpts+1),wallse_i(maxpts+1),
  !     >   wallsi(maxpts+1),wallsiz(maxpts+1,MAXIZS + 1),
  !     >   wallseiz(maxpts+1,MAXIZS + 1),
  !     >   wallsil(maxpts+1)
  !


  real, allocatable, public :: tizs(:,:,:), zeffs(:,:,:), powls(:,:,:), lines(:,:,:),&
       walls(:,:,:), deps(:,:), neros(:,:), elims(:,:,:), sdvs(:,:,:),&
       chemden(:,:), chemizs(:,:), hpowls(:,:,:), hlines(:,:,:), &
       ncore(:,:), nedge(:,:), ntrap(:,:), ndivert(:,:), nmsol(:,:),& 
       wallsn(:), wallse(:), wallse_i(:), wallsi(:), wallsil(:), &
       wallseiz(:,:), wallsiz(:,:)
  
  public allocate_dynam3, deallocate_dynam3


contains 

  subroutine allocate_dynam3

    use global_parameters
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('MOD_DYNAM3','ALLOCATE')

    call allocate_array(tizs,1,maxnks,1,maxnrs,-1,maxizs,'TIZS',ierr)
    call allocate_array(zeffs,maxnks,maxnrs,3,'ZEFFS',ierr)
    call allocate_array(powls,1,maxnks,1,maxnrs,-1,maxizs,'POWLS',ierr)
    call allocate_array(lines,1,maxnks,1,maxnrs,-1,maxizs,'LINES',ierr)
    call allocate_array(walls,1,maxnks,1,maxnrs,0,maxizs+1,'WALLS',ierr)

    call allocate_array(deps,maxnds,maxizs,'DEPS',ierr)
    call allocate_array(neros,maxnds,5,'NEROS',ierr)
    call allocate_array(elims,1,maxnks,1,3,-1,maxizs,'ELIMS',ierr)

    call allocate_array(sdvs,1,maxnks,1,maxnrs,-1,maxizs,'SDVS',ierr)

    call allocate_array(chemden,maxnks,maxnrs,'CHEMDEN',ierr)
    call allocate_array(chemizs,maxnks,maxnrs,'CHEMIZS',ierr)

    call allocate_array(hpowls,1,maxnks,1,maxnrs,0,1,'HPOWLS',ierr)
    call allocate_array(hlines,1,maxnks,1,maxnrs,0,1,'HLINES',ierr)

    call allocate_array(ncore,maxnks,maxnrs,'NCORE',ierr)
    call allocate_array(nedge,maxnks,maxnrs,'NEDGE',ierr)
    call allocate_array(ntrap,maxnks,maxnrs,'NTRAP',ierr)
    call allocate_array(ndivert,maxnks,maxnrs,'NDIVERT',ierr)
    call allocate_array(nmsol,maxnks,maxnrs,'NMSOL',ierr)

    call allocate_array(wallsn,maxpts+1,'WALLSN',ierr)
    call allocate_array(wallse,maxpts+1,'WALLSE',ierr)
    call allocate_array(wallse_i,maxpts+1,'WALLSE_I',ierr)
    call allocate_array(wallsi,maxpts+1,'WALLSI',ierr)

    call allocate_array(wallseiz,maxpts+1,maxizs+1,'WALLSEIZ',ierr)
    call allocate_array(wallsiz,maxpts+1,maxizs+1,'WALLSIZ',ierr)

    call allocate_array(wallsil,maxpts+1,'WALLSIL',ierr)

  end subroutine allocate_dynam3

  subroutine deallocate_dynam3
    implicit none

    call pr_trace('MOD_DYNAM3','DEALLOCATE')

    deallocate(tizs)
    deallocate(zeffs)
    deallocate(powls)
    deallocate(lines)
    deallocate(walls)

    deallocate(deps)
    deallocate(neros)
    deallocate(elims)

    deallocate(sdvs)

    deallocate(chemden)
    deallocate(chemizs)

    deallocate(hpowls)
    deallocate(hlines)

    deallocate(ncore)
    deallocate(nedge)
    deallocate(ntrap)
    deallocate(ndivert)
    deallocate(nmsol)

    deallocate(wallsn)
    deallocate(wallse)
    deallocate(wallse_i)
    deallocate(wallsi)

    deallocate(wallseiz)
    deallocate(wallsiz)

    deallocate(wallsil)

  end subroutine deallocate_dynam3

end module mod_dynam3
