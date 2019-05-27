module mod_dynam1
  use debug_options
  implicit none

  !
  !     -*-fortran-*-
  ! common /dynam1/ ddlims,ddts,ddvoid
  ! save /dynam1/
  double precision,public,allocatable :: ddlims(:,:,:),ddts(:,:,:),ddvoid(:),ddvs(:,:,:)
  !
  
  ! common /dynam1a/ chisq1,chisq2,chisq3,chisq4,chisq5
  ! save /dynam1a/
  !
  double precision,public,allocatable :: chisq1(:),chisq2(:),chisq3(:),chisq4(:),chisq5(:)

  public :: allocate_mod_dynam1,deallocate_mod_dynam1

contains

  subroutine allocate_mod_dynam1
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call pr_trace('mod_dynam1','ALLOCATE')

    call allocate_array(ddlims,1,maxnks,1,maxnrs,-1,maxizs,'ddlims',ierr)
    call allocate_array(ddts,1,maxnks,1,maxnrs,-1,maxizs,'ddts',ierr)
    call allocate_array(ddvs,1,maxnks,1,maxnrs,-1,maxizs,'ddvs',ierr)
    call allocate_array(ddvoid,3,'ddvoid',ierr)
    call allocate_array(chisq1,maxpiniter,'chisq1',ierr)
    call allocate_array(chisq2,maxpiniter,'chisq2',ierr)
    call allocate_array(chisq3,maxpiniter,'chisq3',ierr)
    call allocate_array(chisq4,maxpiniter,'chisq4',ierr)
    call allocate_array(chisq5,maxpiniter,'chisq5',ierr)

  end subroutine allocate_mod_dynam1


  subroutine deallocate_mod_dynam1
    implicit none

    call pr_trace('mod_dynam1','DEALLOCATE')

    if (allocated(ddlims)) deallocate(ddlims)
    if (allocated(ddts)) deallocate(ddts)
    if (allocated(ddvs)) deallocate(ddvs)
    if (allocated(ddvoid)) deallocate(ddvoid)
    if (allocated(chisq1)) deallocate(chisq1)
    if (allocated(chisq2)) deallocate(chisq2)
    if (allocated(chisq3)) deallocate(chisq3)
    if (allocated(chisq4)) deallocate(chisq4)
    if (allocated(chisq5)) deallocate(chisq5)

  end subroutine deallocate_mod_dynam1

end module mod_dynam1
