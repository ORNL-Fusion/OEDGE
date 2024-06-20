!pb  070109  ff(2,2) introduced to avoid array temporaries

      function EIRENE_intp_adas(ad,p1,p2) result(res)
 
      use EIRMOD_precision
      use EIRMOD_comxs, only: adas_data
 
      implicit none
 
      type(adas_data), pointer :: ad
      real(dp), intent(in) :: p1, p2
      real(dp) :: res, rx, ry, ff(2,2)
      integer :: ide, ite
 
      if (p1 <= ad%temp(1)) then
        ite = 1
      else if (p1 >= ad%temp(ad%ntemp)) then
        ite = ad%ntemp-1
      else
        call EIRENE_binsearch_2 (ad%temp, ad%ntemp, p1, ite)
      end if
 
      if (p2 <= ad%dens(1)) then
        ide = 1
      else if (p2 >= ad%dens(ad%ndens)) then
        ide = ad%ndens-1
      else
        call EIRENE_binsearch_2 (ad%dens, ad%ndens, p2, ide)
      end if
 
      rx = (ad%temp(ite+1) - p1) * ad%dte(ite)
      ry = (ad%dens(ide+1) - p2) * ad%dde(ide)
      ff = ad%fit(ite:ite+1,ide:ide+1)
!      call EIRENE_bilinear_int (ad%fit(ite:ite+1,ide:ide+1), rx, ry,
!     .  res)
      call EIRENE_bilinear_int (ff, rx, ry, res)
 
      return
      end function EIRENE_intp_adas
