module variable_wall
use error_handling
implicit none


  ! LIM wall options allowing for a non-constant wall as a function of Y
  integer :: lim_wall_opt
  real :: ywall_start,caw_min
  real, allocatable :: caw_qys(:)


contains



  real function caw_fnc(y,cl,caw)
    implicit none
    

    real :: y,cl,caw

    real :: yabs,ytmp

    caw_fnc = caw

    if (lim_wall_opt.eq.0) then 
       caw_fnc = caw
    elseif (lim_wall_opt.eq.1) then 

       yabs = abs(y)
       if (yabs.gt.cl) then 
          ytmp= 2.0*cl - yabs
       endif
       
       if (ytmp.lt.ywall_start) then 
          caw_fnc = caw 
       else
          caw_fnc = (ytmp-ywall_start) / (cl-ywall_start) * (caw_min - caw) + caw
       endif
    endif

  end function caw_fnc


  subroutine setup_wall(qys,nqys,cl,caw)
    implicit none
    integer :: nqys
    
    real :: qys(nqys)
    real :: cl,caw

    !real,external :: caw_fnc

    integer :: in
    integer :: ierr

    allocate(caw_qys(nqys),stat=ierr)
    
    if (ierr.ne.0) then 
       call errmsg('SETUP_WALL: ERROR ALLOCATING CAW_QYS:',ierr)
    endif


    do in = 1,nqys
       caw_qys(in) = caw_fnc(qys(in),cl,caw)
    end do


  end subroutine setup_wall




end module variable_wall
