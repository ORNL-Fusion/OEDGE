module bfield

  use global_parameters
  

  implicit none

  ! try using private to block propagation of the visibility of global_parameters
  private


  
  real,private :: br(maxnks,maxnrs),bz(maxnks,maxnrs),bt(maxnks,maxnrs)
  integer,private,parameter :: veclen=3

  public:: setup_bvectors,get_bvector,write_bvectors


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




end module bfield



