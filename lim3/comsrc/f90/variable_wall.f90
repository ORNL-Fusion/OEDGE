module variable_wall
  use error_handling
  implicit none


  ! LIM wall options allowing for a non-constant wall as a function of Y
  integer :: lim_wall_opt
  real :: ywall_start,caw_min
  real, allocatable :: caw_qys(:),caw_qys_tans(:)

  private :: caw_fnc

contains



  real function caw_fnc(y,cl,caw)
    implicit none


    real :: y,cl,caw

    real :: yabs,ytmp

    caw_fnc = caw

    if (lim_wall_opt.eq.0) then 
       caw_fnc = caw
    elseif (lim_wall_opt.eq.1) then 

       ytmp = abs(y)

       if (ytmp.gt.cl) then 
          ytmp= 2.0*cl - ytmp
       endif


       if (ytmp.lt.ywall_start) then 
          caw_fnc = caw 
       else
          caw_fnc = (ytmp-ywall_start) / (cl-ywall_start) * (caw_min - caw) + caw
       endif
    endif

  end function caw_fnc


  subroutine setup_wall(qys,nqys,cl,caw)
    use constants
    implicit none
    integer :: nqys

    real :: qys(nqys)
    real :: cl,caw

    !real,external :: caw_fnc

    integer :: in
    integer :: ierr
    real :: dx1,dy1,dx2,dy2
    real, external :: atan2c

    allocate(caw_qys(nqys),stat=ierr)

    if (ierr.ne.0) then 
       call errmsg('SETUP_WALL: ERROR ALLOCATING CAW_QYS:',ierr)
    endif

    allocate(caw_qys_tans(nqys),stat=ierr)

    if (ierr.ne.0) then 
       call errmsg('SETUP_WALL: ERROR ALLOCATING CAW_QYS_TANS:',ierr)
    endif

    do in = 1,nqys
       caw_qys(in) = caw_fnc(qys(in),cl,caw)
    end do

    ! Y goes from 0.0 to 2L - caw_qys is assigned over this range
    ! 0.0 to 2L is mirrored onto 0.0 to -2L with the half-separation at CL
    ! So tangents are calculated for 0.0 to 2L and have to be reflected in 
    ! the vertical axis to get the values when Y<0
    ! A normal to a flat horizontal wall surface in the code is PI/2.0 

    
    do in = 1,nqys

       if (in.eq.1) then 
          dx1 = 1.0
          dy1 = 0.0
          dx2 = qys(in+1)-qys(in)
          dy2 = caw_qys(in+1)-caw_qys(in)

       else if (in.eq.nqys) then 
          dx1 = qys(in)-qys(in-1)
          dy1 = caw_qys(in)-caw_qys(in-1)
          dx2 = 1.0
          dy2 = 0.0


       else

          dx1 = qys(in)-qys(in-1)
          dy1 = caw_qys(in)-caw_qys(in-1)
          dx2 = qys(in+1)-qys(in)
          dy2 = caw_qys(in+1)-caw_qys(in)

       endif

       caw_qys_tans(in) = (atan2c(dy1,dx1) + atan2c(dy2,dx2))/2.0 + PI/2.0

       write(debug_message_data,'(i6,1x,g12.5)') in,caw_qys_tans(in)*RADDEG

       call dbgmsg('CAW_QYS_TANS:',debug_message_data)

    enddo



  end subroutine setup_wall



  subroutine cleanup_wall
    implicit none
    
    deallocate(caw_qys)


  end subroutine cleanup_wall


end module variable_wall
