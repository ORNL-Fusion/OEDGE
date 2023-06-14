module mod_sl_easmesh


implicit none


  ! jdemod - move allocatable arrays from the vacuum.f common block here for now
  !COMMON /EASMESH/ nseg,nwal,rseg,zseg,rwal,zwal
  !    INTEGER          nseg,nwal
  !    DOUBLE PRECISION rseg(MAXNKS),zseg(MAXNKS),
  !   .                 rwal(MAXPTS),zwal(MAXPTS)
      INTEGER          nseg,nwal
      double precision, allocatable :: rseg(:),zseg(:),rwal(:),zwal(:)


contains


  subroutine allocate_mod_sl_easmesh
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call allocate_array(rseg,maxnks,'rseg',ierr)
    call allocate_array(zseg,maxnks,'zseg',ierr)
    call allocate_array(rwal,maxpts,'rwal',ierr)
    call allocate_array(zwal,maxpts,'zwal',ierr)
    
  end subroutine allocate_mod_sl_easmesh

  subroutine deallocate_mod_sl_easmesh
    implicit none

    if (allocated(rseg)) deallocate(rseg)
    if (allocated(zseg)) deallocate(zseg)
    if (allocated(rwal)) deallocate(rwal)
    if (allocated(zwal)) deallocate(zwal)
    

  end subroutine deallocate_mod_sl_easmesh



end module mod_sl_easmesh
