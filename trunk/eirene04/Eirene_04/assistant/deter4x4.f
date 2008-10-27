      function deter4x4 (a)
      use precision
      implicit none

      real(dp), intent(in) :: a(4,4)
      real(dp) :: deter4x4, d11, d12, d13, d14, deter

      d11 = deter(a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d12 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,3),a(3,3),a(4,3),
     .            a(2,4),a(3,4),a(4,4))

      d13 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,4),a(3,4),a(4,4))

      d14 = deter(a(2,1),a(3,1),a(4,1),
     .            a(2,2),a(3,2),a(4,2),
     .            a(2,3),a(3,3),a(4,3))

      deter4x4 = a(1,1)*d11 - a(1,2)*d12 + a(1,3)*d13 - a(1,4)*d14

      return
      end 
      

      
      
      
