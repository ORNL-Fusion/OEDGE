c
c
      function erf (x)
c  s15aef is a nag routine: error function
c  erf(x)=2/sqrt(pi)*int^x_0 dt exp(-(t**2))
      USE PRECISION
      implicit none
      real(dp), intent(in) :: x
      real(dp) :: erf, s15aef
      integer :: ifail
      erf=s15aef(x,ifail)
      return
      end
