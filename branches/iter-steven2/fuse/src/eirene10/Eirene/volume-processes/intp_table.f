      function intp_table(tb,p1,p2) result(res)

      use precision
      use comxs, only: hydkin_data

      implicit none

      type(hydkin_data), pointer :: tb
      real(dp), intent(in) :: p1, p2
      real(dp) :: res, rx
      integer :: ite
      
      if (p1 <= tb%temps(1)) then
        ite = 1
      else if (p1 >= tb%temps(tb%ntemps)) then
        ite = tb%ntemps-1
      else 
        call binsearch_2 (tb%temps, tb%ntemps, p1, ite)
      end if  
      
      res = tb%rates(ite) + (p1-tb%temps(ite))*tb%ratio(ite)

      return
      end function intp_table
