c 
c
      function ranget_eirene(ise)
      USE PRECISION
      implicit none
      integer, intent(inout) :: ise

      integer :: iseed
      common /cmem/ iseed
      integer :: ranget_eirene

      ranget_eirene=iseed
      ise = iseed
      return
      end
