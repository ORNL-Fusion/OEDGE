
c
      function e1 (x)
      USE PRECISION
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: e1, s13aaf
      integer :: ifail
      ifail = 0
      e1 = s13aaf(x, ifail)
      return
      end
