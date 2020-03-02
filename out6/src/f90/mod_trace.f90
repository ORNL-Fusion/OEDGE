module mod_trace
  use debug_options
  
  private 

  integer,public :: mxxnas, mxxnbs
  
      PARAMETER (MXXNAS=8000)
      
      REAL,allocatable,public :: FACTS(:),NORMS(:),CS(:,:)
      REAL,allocatable,public :: RD(:),FN(:),AREAS(:)
      REAL,allocatable,public :: GN(:),DN(:),THETA(:),XN(:)
      REAL,allocatable,public :: WORKS(:),AENDS(:)

!      COMMON /TRACE/ FACTS,NORMS,AREAS,CS,AENDS,
!     >               RD,XN,FN,GN,DN,THETA,WORKS
!      REAL    FACTS(MXXNBS),NORMS(MXXNAS),CS(MXXNAS,MXXNBS)
!      REAL    RD(MXXNAS),FN(MXXNAS),AREAS(MXXNBS)
!      REAL    GN(MXXNAS),DN(MXXNAS),THETA(MXXNAS),XN(MXXNAS)
!      REAL    WORKS(16*MXXNAS+6),AENDS(MXXNAS)

      
      public:: allocate_mod_trace, deallocate_mod_trace
      
contains


  subroutine allocate_mod_trace
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    mxxnbs = maxizs + 2

    call pr_trace('allocate_mod_trace','ALLOCATE')

    call allocate_array(facts,mxxnbs,'facts',ierr)
    call allocate_array(norms,mxxnas,'norms',ierr)
    call allocate_array(cs,mxxnas,mxxnbs,'cs',ierr)
    
    call allocate_array(rd,mxxnas,'rd',ierr)
    call allocate_array(fn,mxxnas,'fn',ierr)

    call allocate_array(areas,mxxnbs,'areas',ierr)

    call allocate_array(gn,mxxnas,'gn',ierr)
    call allocate_array(dn,mxxnas,'dn',ierr)

    call allocate_array(theta,mxxnas,'theta',ierr)
    call allocate_array(xn,mxxnas,'xn',ierr)

    call allocate_array(aends,mxxnas,'aend',ierr)

    call allocate_array(works,16*mxxnas+6,'works',ierr)
    
  end subroutine allocate_mod_trace




  subroutine deallocate_mod_trace
    implicit none

    call pr_trace('allocate_mod_trace','DEALLOCATE')
    if (allocated(facts)) deallocate(facts)
    if (allocated(norms)) deallocate(norms)
    if (allocated(cs)) deallocate(cs)
    if (allocated(rd)) deallocate(rd)
    if (allocated(fn)) deallocate(fn)
    if (allocated(areas)) deallocate(areas)
    if (allocated(gn)) deallocate(gn)
    if (allocated(dn)) deallocate(dn)
    if (allocated(theta)) deallocate(theta)
    if (allocated(xn)) deallocate(xn)
    if (allocated(works)) deallocate(works)
    if (allocated(aends)) deallocate(aends)
    
  end subroutine deallocate_mod_trace


  
end module mod_trace
