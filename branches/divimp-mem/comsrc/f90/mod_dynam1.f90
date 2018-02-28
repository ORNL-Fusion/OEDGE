module mod_dynam1

  use debug_options

  implicit none

  private


  !         include 'dynam1'
  !
  !     -*-Fortran-*-
  !                                                                       
  !      COMMON /DYNAM1/ DDLIMS,DDTS,ddvoid
  !      save /dynam1/
  !      DOUBLE PRECISION DDLIMS(MAXNKS,MAXNRS,-1:MAXIZS),                 
  !     >  DDTS(MAXNKS,MAXNRS,-1:MAXIZS),ddvoid(3)
  !
  !
  !      common /dynam1a/ CHISQ1,CHISQ2,CHISQ3,CHISQ4,CHISQ5 
  !      save /dynam1a/
  !      double precision CHISQ1(maxpiniter),CHISQ2(maxpiniter),
  !     >                 CHISQ3(maxpiniter),
  !     >                 CHISQ4(maxpiniter),CHISQ5(maxpiniter)
  !                                                                       


  real*8, allocatable, public :: ddlims(:,:,:), ddts(:,:,:), ddvoid(:)
  real*8, allocatable, public :: chisq1(:),chisq2(:),chisq3(:),chisq4(:),chisq5(:) 

  public allocate_dynam1, deallocate_dynam1
  

contains

  subroutine allocate_dynam1
    use allocate_arrays
    use global_parameters
    implicit none

    integer :: ierr

    call pr_trace('MOD_DYNAM1','ALLOCATE')

    call allocate_array(ddlims,1,maxnks,1,maxnrs,-1,maxizs,'DDLIMS',ierr)
    call allocate_array(ddts,1,maxnks,1,maxnrs,-1,maxizs,'DDTS',ierr)
    call allocate_array(ddvoid,3,'DDVOID',ierr)

    call allocate_array(chisq1,maxpiniter,'CHISQ1',ierr)
    call allocate_array(chisq2,maxpiniter,'CHISQ1',ierr)
    call allocate_array(chisq3,maxpiniter,'CHISQ1',ierr)
    call allocate_array(chisq4,maxpiniter,'CHISQ1',ierr)
    call allocate_array(chisq5,maxpiniter,'CHISQ1',ierr)

  end subroutine allocate_dynam1



  subroutine deallocate_dynam1
    implicit none

    call pr_trace('MOD_DYNAM1','DEALLOCATE')

    deallocate(ddlims)
    deallocate(ddts)
    deallocate(ddvoid)

    deallocate(chisq1)
    deallocate(chisq2)
    deallocate(chisq3)
    deallocate(chisq4)
    deallocate(chisq5)


  end subroutine deallocate_dynam1


end module mod_dynam1
