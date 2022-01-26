c
c
      function ranset_eirene(ise)
      USE PRECISION
      implicit none
      integer, intent(in) :: ise

      integer :: iseed
      common /cmem/ iseed
      real(dp) :: ranset_eirene
      integer, save :: ifirst=0

      if ((ise <= 0) .and. (ifirst == 0)) then
         iseed = 9876543
      else
         iseed = abs(ise)
      end if

      ifirst = 1
      ranset_eirene=0
      return
      end
