module mod_cadas

  use mod_params

!                                                                       
!     COMMON /CADAS/  YEAR,YEARDF,ICLASS,cdatopt,PTESA,PNESA,PNBS,      
!    >                PNHS,PNZSA,PCOEF,iyearh,iyearz,
!    >                useridh,useridz
!                                                                       
!     INTEGER ICLASS,cdatopt,iyearh,iyearz                       
!     REAL PTESA(MAXNXS),PNESA(MAXNXS),PNBS(MAXNXS),PNHS(MAXNXS)        
!     REAL PNZSA(MAXNXS,0:MAXIZS),PCOEF(MAXNXS,MAXIZS)                  
!     CHARACTER YEAR*2,YEARDF*2,useridh*80,useridz*80             

  implicit none
  private


      INTEGER,public::  ICLASS,cdatopt,iyearh,iyearz                       
      REAL,allocatable,public:: PTESA(:),PNESA(:),PNBS(:),PNHS(:)        
      REAL,allocatable,public::  PNZSA(:,:),PCOEF(:,:)                  
      CHARACTER,public:: YEAR*2,YEARDF*2,useridh*80,useridz*80             
  


  public :: allocate_mod_cadas, deallocate_mod_cadas


contains

  subroutine allocate_mod_cadas
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    call allocate_array(PTESA,maxnxs,'PTESA',ierr)
    call allocate_array(PNESA,maxnxs,'PNESA',ierr)
    call allocate_array(PNBS ,maxnxs,'PNBS',ierr)
    call allocate_array(PNHS ,maxnxs,'PNHS',ierr)
    call allocate_array(PNZSA,1,maxnxs,0,maxizs,'PNZSA',ierr)
    call allocate_array(PCOEF,maxnxs,maxizs,'PCOEF',ierr)

  end subroutine allocate_mod_cadas


  subroutine deallocate_mod_cadas
    use mod_params
    use allocate_arrays
    implicit none

    deallocate(PTESA)
    deallocate(PNESA)
    deallocate(PNBS )
    deallocate(PNHS )
    deallocate(PNZSA)
    deallocate(PCOEF)

  end subroutine deallocate_mod_cadas



end module mod_cadas
