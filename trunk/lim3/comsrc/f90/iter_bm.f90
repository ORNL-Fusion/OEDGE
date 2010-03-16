module iter_bm
  use error_handling
  use common_utilities
  implicit none

  ! ITER BM shaping inputs
  real    :: rtor_setback,rslot_setback,bm_tor_wid,slot_tor_wid,lambda_design

  ! ITER BM inputs for EDGE option 12 - factoring in pitch angle
  real :: bth_bphi_ratio, p_0_value, rho_p_pol, r_ow

  ! ITER BM parameters
  real :: c_lim, y_re, yt_re
  real :: lim_beta, cos_beta, tan_beta
  real :: bm_tor_wid_y,slot_tor_wid_y
  real :: rr_shadow,slot_shadow1,slot_shadow2

  ! ITER BM derived quantities
  real :: xtor_setback, xslot_setback


  ! ITER limiter shape
  integer, parameter :: npts = 1000
  real :: lim_xshift(2)
  integer :: cnts(2)
  real, allocatable :: x_lim(:,:),y_lim(:,:)

  private npts,cnts,x_lim,y_lim

  save

contains


  subroutine calc_iter_limiter_parameters
    use constants
    implicit none
    !
    ! See "ITER Limiter Notes. Ideal Shape of the face of the outer wall BMs for startup/rampdown. 
    !      Set of 18 BMs in a toroidal ring at the outer wall" by P.C. Stangeby for details on the 
    !      formulae and calculations used here. Also see notes and comments in edge.f
    !
    real :: t_rr, p_rr, t_ridge, p_ridge1, p_ridge2


    y_re = ( bm_tor_wid * (1-exp(-rslot_setback/lambda_design)) +slot_tor_wid*(1-exp(-rtor_setback/lambda_design)))/ &
         ( (1-exp(-rslot_setback/lambda_design)) +(1-exp(-rtor_setback/lambda_design)))

    c_lim = (lambda_design * (1-exp(-rtor_setback/lambda_design))) / (bm_tor_wid - y_re)


    ! Calculate the pitch and cosine of the angle based on the input magnetic field ratio

    lim_beta = atan(bth_bphi_ratio)
    cos_beta = cos(lim_beta)
    tan_beta = bth_bphi_ratio

    if (cos_beta.ne.0.0) then 
       yt_re = y_re/cos_beta
       ! Note: the value assigned to bm_tor_wid_y is the Y coordinate of the edge of the limiter with Y=0 located at t_re
       !       the value assigned to slot_tor_wid_y is the Y coordinate of the edge of the slot with Y=0 located at t_re
       !       Normally ... bm_tor_wid_y < 0 ... however since all limiter shapes are recorded with values >0 ... an abs is used
       bm_tor_wid_y = abs((y_re-bm_tor_wid)/cos_beta)
       slot_tor_wid_y = abs((y_re-slot_tor_wid)/cos_beta)
    else
       yt_re = y_re
       bm_tor_wid_y = abs(y_re-bm_tor_wid)
       slot_tor_wid_y = abs(y_re-slot_tor_wid)
    endif

    ! Calculate the shadow distances 

    ! right neighbour ridge
    t_rr = 2.0 * (2.0 * PI * r_ow / 36.0)  - y_re
    p_rr = p_0_value + tan_beta * t_rr

    rr_shadow = p_rr**2 / (2.0 * rho_p_pol)


    ! slot due to self-shadowing in re-entrant region
    t_ridge = y_re
    p_ridge1 = p_0_value + tan_beta * t_ridge
    p_ridge2 = p_0_value - tan_beta * t_ridge

    slot_shadow1 = p_ridge1**2 / (2.0 * rho_p_pol)
    slot_shadow2 = p_ridge2**2 / (2.0 * rho_p_pol)

  end subroutine calc_iter_limiter_parameters



  subroutine calc_iter_limiter_shape
    implicit none

    ! This routine calculates the shape of the limiter using Peter's notes assuming that Y is not the same as t but differs by an angle beta
    ! It also incorporates the effect of p not equal to zero ... p is specified in the input. 
    ! The limiter tip is shifted to 0,0
    ! The limiter shape is calculated numerically as X(Y) and is then interpolated to give the actual required limiter shape in terms of Y(X)

    integer :: ierr
    real :: delta_t, t_start
    real :: t, p, ft, gp
    integer :: in
    real :: maxx
    real :: xtmp,ytmp

    allocate(x_lim(npts,2),stat=ierr) 
    allocate(y_lim(npts,2),stat=ierr)


    !delta_t = (bm_tor_wid - y_re) / real(npts/2)

    delta_t = (bm_tor_wid - slot_tor_wid) / real(npts)

    t_start = slot_tor_wid


    !
    ! Zero point counts for each side of the limiter
    !
    cnts = 0

    !
    ! Note: we have been using the convention in LIM that the reentrant region is Y>0 and the rest Y<0. This is the opposite
    !       of the convention used in Peter's document so care must be taken in mapping the shape. In this case, we will store
    !       the Y<0
    !
    !
    ! As per the convention in edge.f ... side 1 is the Y<0 limiter face and side 2 is the Y>0 face. All Y values stored in the 
    ! limiter arrays are greater than 0.0
    !

    do in = 1,npts+1

       t = t_start + (in-1) * delta_t 

       p = p_0_value + bth_bphi_ratio * t

       if (t.le.y_re) then 
          ft = - lambda_design * log ( 1.0 + c_lim * (t - y_re)/lambda_design)
       else
          ft = - lambda_design * log ( 1.0 - c_lim * (t - y_re)/lambda_design)
       endif

       gp = p**2 / (2.0 * rho_p_pol) 

       ytmp = (y_re - t) / cos_beta

       xtmp = -(ft+gp)

       ! Need to record xtmp for the slot edge and the outer edge of the bm so that the 
       ! edge locations of the limiter shape can be properly detected. These values will need
       ! to be modified by the xshift value found below. 

       if (in.eq.1) then 
          xslot_setback = xtmp
       elseif (in.eq.npts+1) then 
          xtor_setback  = xtmp
       endif


       write(6,'(a,i6,10(1x,g18.6))') 'LIM12:',in,t,p,ft,gp,ytmp,xtmp,y_re,delta_t

       if (ytmp.eq.0.0) then 

          ! Add any point for ytmp=0.0 to both sides of the limiter

          cnts(1) = cnts(1) + 1
          x_lim(cnts(1),1) = xtmp
          y_lim(cnts(1),1) = abs(ytmp)

          cnts(2) = cnts(2) + 1
          x_lim(cnts(2),2) = xtmp
          y_lim(cnts(2),2) = abs(ytmp)

       elseif (ytmp.lt.0.0) then 
          ! side 1 is the Y < 0 limiter face
          cnts(1) = cnts(1) + 1
          x_lim(cnts(1),1) = xtmp
          y_lim(cnts(1),1) = abs(ytmp)


       elseif (ytmp.gt.0.0) then 

          ! side 2 is the Y < 0 limiter face
          cnts(2) = cnts(2) + 1
          x_lim(cnts(2),2) = xtmp
          y_lim(cnts(2),2) = abs(ytmp)

       endif

    end do

    ! Adjust X values so that the maximum value is 0.0

    ! Determine the location of the maximum value in x_lim(:,1) - it should be zero - force the maximum point to be the limiter tip


    in  = maxloc(x_lim(1:cnts(1),1),dim=1) 

    !write(0,*) 'in=',in

    maxx = x_lim(in,1)

    !write(0,*) 'maxx=',maxx

    lim_xshift(1) = maxx

    !write(0,*) 'lim_xshift(1)=',lim_xshift(1)

    x_lim(:,1) = x_lim(:,1) - maxx
    x_lim(in,1) = 0.0
    y_lim(in,1) = 0.0


    ! Determine the location of the maximum value in x_lim(:,1) - it should be zero - force the maximum point to be the limiter tip

    in   = maxloc(x_lim(1:cnts(2),2),dim=1)

    maxx = x_lim(in,2)
    lim_xshift(2) = maxx
    x_lim(:,2) = x_lim(:,2) - maxx
    x_lim(in,2) = 0.0
    y_lim(in,2) = 0.0

    ! Adjust x setback values for maxx
    xslot_setback = xslot_setback - maxx
    xtor_setback  = xtor_setback - maxx


    ! Re-order elements to ascending X values

    if (x_lim(1,1).gt.x_lim(cnts(1),1)) then 

       call sort_arrays(0,cnts(1),x_lim(:,1),y_lim(:,1))       

    endif

    ! Re-order elements to ascending X values

    if (x_lim(1,2).gt.x_lim(cnts(2),2)) then 

       call sort_arrays(0,cnts(2),x_lim(:,2),y_lim(:,2))       

    endif

    write(6,'(a,2(2x,g12.5))') 'X Setback values (slot,tor)',xslot_setback, xtor_setback


    ! Print out the limiter shape 

    write(6,'(a,i6)') 'Y<0 Limiter shape',cnts(1)

    do in = 1,cnts(1)
       write(6,'(i6,2(1x,g18.6))') in, x_lim(in,1),y_lim(in,1)
    end do
    write(6,'(a,i6)') 'Y>0 Limiter shape',cnts(2)

    do in = 1,cnts(2)
       write(6,'(i6,2(1x,g18.6))') in, x_lim(in,2),y_lim(in,2)
    end do

  end subroutine calc_iter_limiter_shape

  real function iter_limiter_shape(x,j)
    implicit none
    real :: tmp
    real :: x
    integer :: j
    ! Returns the Y value from the appropriate side of the limiter shape interpolated to the specified X value

    call fitter(cnts(j),x_lim(:,j),y_lim(:,j),1,x,tmp,'LINEAR')

    write(6,'(a,2i6,3(1x,g18.6))') 'FITTED SHAPE:',j,cnts(j),x,tmp

    iter_limiter_shape = tmp

  end function iter_limiter_shape


  subroutine clean_up_iter_limiter_shape
    implicit none
    ! Desllocates the storage assigned to save the limiter shape
    if (allocated(x_lim))  deallocate(x_lim)
    if (allocated(y_lim))  deallocate(y_lim)

  end subroutine clean_up_iter_limiter_shape



end module iter_bm

