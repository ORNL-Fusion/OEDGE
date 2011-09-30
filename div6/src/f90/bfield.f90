module bfield

  use global_parameters
  

  implicit none

  ! try using private to block propagation of the visibility of global_parameters
  private


  
  real,private :: br(maxnks,maxnrs),bz(maxnks,maxnrs),bt(maxnks,maxnrs)
  integer,private,parameter :: veclen=3

  public:: setup_bvectors,get_bvector


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



end module bfield



