      subroutine fxdrchr (iun, t)
      
      integer, intent(in) :: iun
      character(len=*), intent(in) :: t

      return
      end subroutine fxdrchr


      subroutine fxdrcls (iun)
      
      integer, intent(in) :: iun

      return
      end subroutine fxdrcls


      subroutine fxdrdbl (iun, ar, n)
      
      integer, intent(in) :: iun, n
      double precision, intent(in) :: ar(n)

      return
      end subroutine fxdrdbl


      subroutine fxdrint (iun, iar, n)
      
      integer, intent(in) :: iun, n
      integer, intent(in) :: iar(n)

      return
      end subroutine fxdrint


      subroutine fxdrlog (iun, lar, n)
      
      integer, intent(in) :: iun, n
      logical, intent(in) :: lar(n)

      return
      end subroutine fxdrlog


      subroutine fxdropn (name, code, iun)
      
      integer, intent(in) :: iun
      character(len=*), intent(in) :: name, code
      
      return
      end subroutine fxdropn

