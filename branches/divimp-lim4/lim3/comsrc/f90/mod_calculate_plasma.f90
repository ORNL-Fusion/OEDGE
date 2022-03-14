module mod_calculate_plasma

  use mod_plasma_data

  private


  public:: calculate_plasma,assign_plasma

  real :: qtim_local


contains

  subroutine setup_solvers(qtim)
    implicit none
    ! this is required to bring in some global variables in an easy fashion instead of adding
    ! it to every subroutine in the call stack

    ! It is also used for any other solver setup required
    qtim_local = qtim


  end subroutine setup_solvers
  

  subroutine calculate_plasma
    implicit none
    
    ! code here calls SOLEDGE or SOL22 as appropriate
    ! also allocates storage for the plasma solution

    if (sol22_opt.ne.0.and.nsol22_opt.gt.0) then 
    !
    !     This code calculates the plasma conditions for sections of the simulation
       !     volume using SOL22. There are several scenarios.
       !
!         
! Initialize some output options in SOL22 using values from slcom
         call init_solcommon(0,0)

         call sol22



    endif

       if (soledge_opt.eq.1.and.colprobe3d.eq.1) then 

         ! plasma is calculated from lower absorbing surface to
         ! upper absorbing surface - this allows for
         ! asymmetric placement of the probe
         !call init_soledge(yabsorb1a,yabsorb2a)
         call init_soledge()
         !call init_soledge(-cl,cl)
         
         !if (vary_absorb.eq.1) then
           ! Find the x index where the step happens. I think this is IPOS?
           !ix_step1 = ipos(xabsorb1a_step, xs, nxs-1)
           !ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
           !write(0,*) 'ix_step1 = ',ix_step1,'(x = ',xs(ix_step1),')'
           !write(0,*) 'ix_step2 = ',ix_step2,'(x = ',xs(ix_step2),')'
           
           ! Call soledge for the plasma from the wall to the step.
           !call soledge(1, ix_step1, qtim)
           
           ! Call soledge for the plasma from the step to the top.
           !write(0,*) 'second soledge call'
           !call soledge(ix_step1+1, nxs, qtim)
           
           ! deallocate storage here instead of inside soledge code.
           !call end_soledge
           
         !else
           ! Just do the normal option with one absorbing wall.
           call soledge(1,nxs,qtim)
           !call soledge(1,nxs/2,qtim)
         !endif

      endif

    
    
    
    

  end subroutine calculate_plasma

  subroutine assign_plasma(y_lim,ne_lim,te_lim,ti_lim,vb_lim,ef_lim,teg_lim,tig_lim)
    implicit none
    real :: y_lim,ne_lim,te_lim,ti_lim,vb_lim,ef_lim,teg_lim,tig_lim
    ! This routine interpolates the plasma solution to obtain the
    ! plasma values at the specific coordinate
    real :: y_shift
    
    integer :: iy

    y_shift = y_lim-lower_y_bound ! obtain y coordinate on solved plasma solution
    
    if (y_shift.lt.sd(1)) then ! Y is less than plasma solution space - assign first cell

       ne_lim = ne(1)
       te_lim = te(1)
       ti_lim = ti(1)
       vb_lim = vb(1)
       ef_lim = ef(1)
       teg_lim = teg(1)
       tig_lim = tig(1)
       
    elseif (y_shift.gt.sd(maxn)) then ! Y is greater than plasma solution space - assign last cell

       ne_lim = ne(maxn)
       te_lim = te(maxn)
       ti_lim = ti(maxn)
       vb_lim = vb(maxn)
       ef_lim = ef(maxn)
       teg_lim = teg(maxn)
       tig_lim = tig(maxn)

    else   ! interpolate - note doing the interpolation even if youts=sd without the end points

       iy = ipos(y_shift,s_axis,maxn)
       fact = (y_shift-sd(iy-1))/(sd(iy)-sd(iy-1))

       ne_lim = ne(iy-1) + fact*(ne(iy)-ne(iy-1))
       te_lim = te(iy-1) + fact*(te(iy)-te(iy-1))
       ti_lim = ti(iy-1) + fact*(ti(iy)-ti(iy-1))
       vb_lim = vb(iy-1) + fact*(vb(iy)-vb(iy-1))
       ef_lim = ef(iy-1) + fact*(ef(iy)-ef(iy-1))
       teg_lim = teg(iy-1) + fact*(teg(iy)-teg(iy-1))
       tig_lim = tig(iy-1) + fact*(tig(iy)-tig(iy-1))

    endif       
    
  end subroutine assign_plasma


  
  subroutine calculate_tgrad_e(ix,pz)
    use mod_params
    use mod_comtor
    use mod_comxyt
    ! this routine takes the plasma solution and calculates the temperature gradient
    ! and electric field
    ! the midpoint may have a discontinuity so these are set to zero
    implicit none
    real*8 :: tgscal
    real*8 :: ds1,ds2,dp1,dp2,dt1,dt2,nb1,nb2,e1,e2
    real*8 :: grad1,grad2

    ! proper scaling of teg and e for the simulation (QS is radial time step multiplier if using accelerated time in core or SOL)
    TGSCAL = EMI/CRMI * QTIM_local *QTIM_local * QS(IQXS(IX)) * QS(IQXS(IX)

    !DSTEP = TGSCAL * QS(IQXS(IX)) * QS(IQXS(IX))
    !EFACT  = QTIM_local * QTIM_local * EMI / CRMI -> same as dstep

    ! teg
    do iy = 1,maxn
       if (iy.eq.1) then
          grad1 = (te(iy+1)-te(iy))/(sd(iy+1)-sd(iy))
          teg(iy) = tgscal * grad1

       elseif (iy.eq.maxn) then

          grad2 = (te(iy)-te(iy-1))/(sd(iy)-sd(iy-1))
          teg(iy) = tgscal * grad2 
       else
          grad1 = (te(iy+1)-te(iy))/(sd(iy+1)-sd(iy))
          grad2 = (te(iy)-te(iy-1))/(sd(iy)-sd(iy-1))
          teg(iy) = tgscal * (grad1+grad2)/2.0 

       endif

    end do

    ! tig
    do iy = 1,maxn
       if (iy.eq.1) then
          grad1 = (ti(iy+1)-ti(iy))/(sd(iy+1)-sd(iy))
          tig(iy) = tgscal * grad1
       elseif (iy.eq.maxn) then
          grad2 = (ti(iy)-ti(iy-1))/(sd(iy)-sd(iy-1))
          tig(iy) = tgscal * grad2 
       else
          grad1 = (ti(iy+1)-ti(iy))/(sd(iy+1)-sd(iy))
          grad2 = (ti(iy)-ti(iy-1))/(sd(iy)-sd(iy-1))
          tig(iy) = tgscal * (grad1+grad2)/2.0 
       endif

    end do

    ! ef
    do iy=1,maxn
       if (y.eq.1) then

          ds1 = y(iy+1) - y(iy)
          dp1 = ne(iy+1)*te(iy+1) - ne(iy)*te(iy)
          dt1 = te(iy+1) - te(iy)
          nb1 = 0.5*(ne(iy+1) + ne(iy))

          if (ds1.ne.0.and.nb1.ne.0) then 
             e1 = -(1.0/nb1)*dp1/ds1 - 0.71 * dt1/ds1
          else
             e1 = 0.0
          endif

          ef(iy) = tgscal * e1

       elseif (y.eq.maxn) then

          ds2 = y(iy) - y(iy-1)
          dp2 = ne(iy)*te(iy) - ne(iy-1)*te(iy-1)
          dt2 = te(iy) - te(iy-1)
          nb2 = 0.5*(ne(iy) + ne(iy-1))

          if (ds2.ne.0.and.nb2.ne.0) then 
             e2 = -(1.0/nb2)*dp2/ds2 - 0.71 * dt2/ds2
          else
             e2 = 0.0
          endif

          ef(iy) = tgscal * e2

       else

          ds1 = y(iy+1) - y(iy)
          dp1 = ne(iy+1)*te(iy+1) - ne(iy)*te(iy)
          dt1 = te(iy+1) - te(iy)
          nb1 = 0.5*(ne(iy+1) + ne(iy))

          if (ds1.ne.0.and.nb1.ne.0) then 
             e1 = -(1.0/nb1)*dp1/ds1 - 0.71 * dt1/ds1
          else
             e1 = 0.0
          endif

          ds2 = y(iy) - y(iy-1)
          dp2 = ne(iy)*te(iy) - ne(iy-1)*te(iy-1)
          dt2 = te(iy) - te(iy-1)
          nb2 = 0.5*(ne(iy) + ne(iy-1))

          if (ds2.ne.0.and.nb2.ne.0) then 
             e2 = -(1.0/nb2)*dp2/ds2 - 0.71 * dt2/ds2
          else
             e2 = 0.0
          endif

          ef(iy) = tgscal * 0.5 * (e1+e2)

       endif

    end do




  end subroutine calculate_tgrad_e


  
end module mod_calculate_plasma
