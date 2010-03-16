      function intp_adas(ad,p1,p2) result(res)

      use precision
      use comxs, only: adas_data

      implicit none

      type(adas_data), pointer :: ad
      real(dp), intent(in) :: p1, p2
      real(dp) :: res, rx, ry
      integer :: ide, ite
      
      if (p1 <= ad%temp(1)) then
        ite = 1
      else if (p1 >= ad%temp(ad%ntemp)) then
        ite = ad%ntemp-1
      else 
        call binsearch_2 (ad%temp, ad%ntemp, p1, ite)
      end if  

      if (p2 <= ad%dens(1)) then
        ide = 1
      else if (p2 >= ad%dens(ad%ndens)) then
        ide = ad%ndens-1
      else 
        call binsearch_2 (ad%dens, ad%ndens, p2, ide)
      end if  
      
      rx = (ad%temp(ite+1) - p1) * ad%dte(ite)
      ry = (ad%dens(ide+1) - p2) * ad%dde(ide)
      call bilinear_int (ad%fit(ite:ite+1,ide:ide+1), rx, ry, res)

      return
      end function intp_adas
