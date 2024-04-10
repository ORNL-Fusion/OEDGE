module Msort

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A,B)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer :: iq

  if(size(A) > 1) then
     call Partition(A,B, iq)
     call QsortC(A(:iq-1),B(:iq-1))
     call QsortC(A(iq:),B(iq:))
  endif
end subroutine QsortC

subroutine Partition(A,B, marker)
  integer, intent(in out), dimension(:) :: A
  integer, intent(in out), dimension(:) :: B
  integer, intent(out) :: marker
  integer :: i, j
  integer :: temp
  integer tempi
  integer :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
        ! exchange B(i) and B(j)
        tempi = B(i)
        B(i) = B(j)
        B(j) = tempi
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module Msort
