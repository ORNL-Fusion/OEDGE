module mod_limpoly

  use mod_params


!c
!      common /limpoly/ npolyp,korpg,nvertp,xvertp,yvertp,
!     >                 kareas
!      integer npolyp
!      integer korpg(maxnxs,maxnys) 
!      integer nvertp(maxnxs*maxnys)
!      real xvertp(4,maxnxs*maxnys)  
!      real yvertp(4,maxnxs*maxnys)  
!      real kareas(maxnxs,maxnys) 
!c


  implicit none
  private

      integer,public:: npolyp
      integer,allocatable,public:: korpg(:,:) 
      integer,allocatable,public:: nvertp(:)
      real,allocatable,public:: xvertp(:,:)  
      real,allocatable,public:: yvertp(:,:)  
      real,allocatable,public:: kareas(:,:) 
  
  public :: allocate_mod_limpoly, deallocate_mod_limpoly


contains

  subroutine allocate_mod_limpoly
    use mod_params
    use allocate_arrays
    implicit none
    integer :: ierr

    !call allocate_array(DTEV  ,maxnxs,'DTEV',ierr)
     call allocate_array(korpg ,maxnxs,maxnys,  'korpg ',ierr)
     call allocate_array(nvertp,maxnxs*maxnys,  'nvertp',ierr)
     call allocate_array(xvertp,4,maxnxs*maxnys,'xvertp',ierr)
     call allocate_array(yvertp,4,maxnxs*maxnys,'yvertp',ierr)
     call allocate_array(kareas,maxnxs,maxnys,  'kareas',ierr)


  end subroutine allocate_mod_limpoly


  subroutine deallocate_mod_limpoly
    use mod_params
    use allocate_arrays
    implicit none

    !deallocate()
    deallocate(korpg )
    deallocate(nvertp)
    deallocate(xvertp)
    deallocate(yvertp)
    deallocate(kareas)

  end subroutine deallocate_mod_limpoly



end module mod_limpoly
