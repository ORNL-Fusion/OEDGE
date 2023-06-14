module bfield

  use mod_params
  

  implicit none

  ! try using private to block propagation of the visibility of mod_params
  private


  
  real,allocatable,private :: br(:,:),bz(:,:),bt(:,:)
  integer,private,parameter :: veclen=3

  public:: setup_bvectors,get_bvector,write_bvectors,allocate_bfield,deallocate_bfield


  save

contains


  subroutine setup_bvectors
    implicit none


    call set_bcomponents(br,bz,bt)


  end subroutine setup_bvectors

  subroutine get_bvector(ik,ir,b)
    implicit none
    integer :: ik,ir
    real :: b(veclen)

    b(1) = br(ik,ir)
    b(2) = bz(ik,ir)
    b(3) = bt(ik,ir)

  end subroutine get_bvector
  
  subroutine write_bvectors(of,nrs,nks)
    implicit none
    integer :: of,nrs,nks(nrs)
    integer :: ik,ir

      write (of,10) 'BR:'
      write (of,500) ((br(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      write (of,10) 'BZ:'
      write (of,500) ((bz(ik,ir),ik=1,nks(ir)),ir=1,nrs)

      write (of,10) 'BT:'
      write (of,500) ((bt(ik,ir),ik=1,nks(ir)),ir=1,nrs)
    

  10  format(a)
 500  format(6e18.10)


  end subroutine write_bvectors

  subroutine allocate_bfield
    use allocate_arrays
    implicit none
    integer :: ierr
    call allocate_array(br,maxnks,maxnrs,'br',ierr)
    call allocate_array(bz,maxnks,maxnrs,'bz',ierr)
    call allocate_array(bt,maxnks,maxnrs,'bt',ierr)
    
  end subroutine allocate_bfield

  subroutine deallocate_bfield
    implicit none
    if (allocated(br)) deallocate(br)
    if (allocated(bz)) deallocate(bz)
    if (allocated(bt)) deallocate(bt)
    
  end subroutine deallocate_bfield



end module bfield



