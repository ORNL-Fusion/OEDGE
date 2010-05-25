c
c
      subroutine talusr (ICOUNT,VECTOR,TALTOT,TALAV,
     .              TXTTL,TXTSP,TXTUN,ILAST,*)
      USE PRECISION
      USE PARMMOD
      USE CGRID
      USE CGEOM
      implicit NONE
      integer, intent(in) :: icount
      integer, intent(out) :: ilast
      real(dp), intent(in) :: vector(*), TALTOT, TALAV
      character(len=*) :: txttl,txtsp,txtun
      integer :: i

      ilast=1
      return 1
      end
