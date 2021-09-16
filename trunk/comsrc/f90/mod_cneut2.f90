module mod_cneut2
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /cneut2/ kmfss ,krmax ,kmfps , kmfpws, kmfcs, kmfcws,kflux ,kener ,kyield,&
  !     kfy,kheat ,kcum,kcum1 ,krmaxw, neut2d_prob, neut2d_index,neut2d_src, neut2d_num
  !
  ! save /cneut2/
  real,public :: neut2d_src
  real,public,allocatable :: kflux(:),kener(:),kyield(:),kmfps(:),kheat(:),kmfcs(:),&
       kfy(:),kcum(:),krmax(:),kmfss(:),kcum1(:),krmaxw(:),kmfpws(:),kmfcws(:),neut2d_prob(:)
  integer,public :: neut2d_num
  integer,public,allocatable :: neut2d_index(:,:)

  public :: allocate_mod_cneut2,deallocate_mod_cneut2

contains

  subroutine allocate_mod_cneut2
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_cneut2','ALLOCATE')

    call allocate_array(kflux,maxnds,'kflux',ierr)
    call allocate_array(kener,maxnds,'kener',ierr)
    call allocate_array(kyield,maxnds,'kyield',ierr)
    call allocate_array(kmfps,maxnds,'kmfps',ierr)
    call allocate_array(kheat,maxnds,'kheat',ierr)
    call allocate_array(kmfcs,maxnds,'kmfcs',ierr)
    call allocate_array(kfy,maxnds,'kfy',ierr)
    call allocate_array(kcum,maxnds,'kcum',ierr)
    call allocate_array(krmax,maxnds,'krmax',ierr)
    call allocate_array(kmfss,maxnds,'kmfss',ierr)
    call allocate_array(kcum1,maxnds,'kcum1',ierr)
    call allocate_array(krmaxw,maxpts,'krmaxw',ierr)
    call allocate_array(kmfpws,maxpts,'kmfpws',ierr)
    call allocate_array(kmfcws,maxpts,'kmfcws',ierr)
    call allocate_array(neut2d_prob,maxnrs*maxnks,'neut2d_prob',ierr)
    call allocate_array(neut2d_index,maxnrs*maxnks,2,'neut2d_index',ierr)

  end subroutine allocate_mod_cneut2


  subroutine deallocate_mod_cneut2
    implicit none

    call pr_trace('mod_cneut2','DEALLOCATE')

    if (allocated(kflux)) deallocate(kflux)
    if (allocated(kener)) deallocate(kener)
    if (allocated(kyield)) deallocate(kyield)
    if (allocated(kmfps)) deallocate(kmfps)
    if (allocated(kheat)) deallocate(kheat)
    if (allocated(kmfcs)) deallocate(kmfcs)
    if (allocated(kfy)) deallocate(kfy)
    if (allocated(kcum)) deallocate(kcum)
    if (allocated(krmax)) deallocate(krmax)
    if (allocated(kmfss)) deallocate(kmfss)
    if (allocated(kcum1)) deallocate(kcum1)
    if (allocated(krmaxw)) deallocate(krmaxw)
    if (allocated(kmfpws)) deallocate(kmfpws)
    if (allocated(kmfcws)) deallocate(kmfcws)
    if (allocated(neut2d_prob)) deallocate(neut2d_prob)
    if (allocated(neut2d_index)) deallocate(neut2d_index)

  end subroutine deallocate_mod_cneut2

end module mod_cneut2