module mod_comvu

  use mod_params

!C                                                                               
!      COMMON /COMVU/ XPPPP,YPPPP,SSQ,XPWID,YPWID                                
!      REAL XPPPP(MAXNXS,-MAXNYS:MAXNYS),YPPPP(MAXNXS,-MAXNYS:MAXNYS)            
!      REAL SSQ  (MAXNXS,-MAXNYS:MAXNYS)                                         
!      REAL XPWID(MAXNXS,-MAXNYS:MAXNYS),YPWID(MAXNXS,-MAXNYS:MAXNYS)            
  
  
  implicit none
  private

      REAL,ALLOCATABLE,PUBLIC:: XPPPP(:,:),YPPPP(:,:)            
      REAL,ALLOCATABLE,PUBLIC:: SSQ  (:,:)                                         
      REAL,ALLOCATABLE,PUBLIC:: XPWID(:,:),YPWID(:,:)                  
  
  public :: allocate_mod_comvu, deallocate_mod_comvu


contains

  subroutine allocate_mod_comvu
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)

      call allocate_array(XPPPP,1,MAXNXS,-MAXNYS,MAXNYS,'XPPPP',ierr)
      call allocate_array(YPPPP,1,MAXNXS,-MAXNYS,MAXNYS,'YPPPP',ierr)
      call allocate_array(SSQ  ,1,MAXNXS,-MAXNYS,MAXNYS,'SSQ  ',ierr)
      call allocate_array(XPWID,1,MAXNXS,-MAXNYS,MAXNYS,'XPWID',ierr)
      call allocate_array(YPWID,1,MAXNXS,-MAXNYS,MAXNYS,'YPWID',ierr)

  end subroutine allocate_mod_comvu


  subroutine deallocate_mod_comvu
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()
    deallocate(XPPPP)
    deallocate(YPPPP)
    deallocate(SSQ  )
    deallocate(XPWID)
    deallocate(YPWID)

  end subroutine deallocate_mod_comvu



end module mod_comvu
