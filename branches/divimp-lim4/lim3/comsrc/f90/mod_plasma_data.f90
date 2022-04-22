module mod_plasma_data

  implicit none
  public

  integer :: maxn = 1000
  integer :: midn        ! calculated in the axis routine
  
  real*8,allocatable:: sd(:),spts(:),ne(:),te(:),ti(:),ef(:),vb(:),teg(:),tig(:),ga(:)
  !real*8,allocatable :: neg(:),ga(:)
  real,allocatable :: s_axis(:)

  real*8 :: n_bnd_lower, te_bnd_lower, ti_bnd_lower, n_bnd_upper, te_bnd_upper, ti_bnd_upper

  real :: lower_y_bound, upper_y_bound , ring_length
  
  !integer :: data_count  ... maxn is used as a local in soledge - may need to use another variable name
  !integer :: maxn = 1000 !initialize to 1000 but adjust depending on simulaton scale size
  ! alternative is to use actual y axes but that has its own challenges. 

  real*8,allocatable :: y_axis(:)
  real*8 :: tg_scale
  real*8 :: ef_scale

  ! LIM variables copied for local use
  real :: qtim_local, crmb_local
  integer :: cizb_local
  integer :: ixout_local
  real :: yscale_local
  integer :: cprint_local = 0
  
contains


  subroutine set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)
    implicit none
    real :: n1,te1,ti1,n2,te2,ti2
    
    n_bnd_lower = n1
    te_bnd_lower = te1
    ti_bnd_lower = ti1

    n_bnd_upper = n2
    te_bnd_upper = te2
    ti_bnd_upper = ti2

  end subroutine set_boundary_conditions

  subroutine set_plasma_data_axis(youts,start_i,end_i,max_i,bnd1,bnd2,axis_opt)
    use allocate_arrays
    implicit none
    integer :: start_i, end_i,max_i,axis_opt   
    real :: youts(-max_i:max_i),bnd1,bnd2

    integer :: pt_cnt,ierr,in
    integer,external :: iposq
    
    ! This routine extracts the Y coordinates for which a plasma solution is required - only ones between the bounds are
    ! included and the bounds are added to the axis if not already present

    ! Determine number of points for axis allocation

    lower_y_bound = bnd1
    upper_y_bound = bnd2

    ring_length = upper_y_bound - lower_y_bound


    
    if (allocated(y_axis)) deallocate(y_axis)

    if (axis_opt.eq.0) then

       ! figure out a value for maxn that gives decent resolution ~ 0.01m?
       ! Minimum 1000 and max 10000 ... 10,000 gives 0.01 resolution for a 100m ring length
       pt_cnt =   max(min(int(ring_length/0.01),1000),10000) 
       
       ! allocate axis array
       call allocate_array(y_axis,pt_cnt,'Local actual Y axis',ierr)

       ! this should give y_axis(1) = bnd1 and y_axis(pt_cnt) = bnd2
       do in = 1,pt_cnt
          y_axis(in) = (in-1) * ring_length/(pt_cnt-1) + bnd1
       end do


    elseif (axis_opt.eq.1) then 

       pt_cnt = 0
       do in = start_i, end_i
          if (youts(in).gt.bnd1.and.youts(in).lt.bnd2) then  ! this explicitly leaves out points=bnd since these are added after as a requirement
             pt_cnt = pt_cnt+1
          endif
       end do
       ! Add 2 points for the end points bnd1,bnd2
       pt_cnt = pt_cnt+2

       ! allocate axis array
       call allocate_array(y_axis,pt_cnt,'Local actual Y axis',ierr)

       y_axis(1) = bnd1
       do in=start_i,end_i
          if (youts(in).gt.bnd1.and.youts(in).lt.bnd2) then
             y_axis(in+1) = youts(in)
          endif
       end do
       y_axis(pt_cnt) = bnd2

    endif

    ! Calculate the internal axis arrays spts and sd (sd is used in soledge and spts in sol22) 
    maxn = pt_cnt
    
    ! set equivalent SOL22 variable for array size - this is done at the start of the SOL22 routine
    !mxspts = maxn
    
    ! y_axis array should start at bnd1 and end at bnd2
    ! calculate shifted y-axis for soledge/sol22 starting at S=0
    ! calculate and allocate both for now - they hold identical values and the only reason there are two is to avoid
    ! having to go through soledge or sol22 and change every occurrence of the axis arrays - that process is likely to generate
    ! bugs

    call allocate_plasma_data(maxn)

    ! Assign axes - explicitly set the end points to avoid the possibility of numerical rounding giving end points that are not 0.0 and ring_length
    sd(1) = 0.0
    do in = 2,maxn-1
       sd(in) = y_axis(in)-bnd1
    end do
    sd(maxn) = ring_length

    !write(6,*) 'Plasma data axis:',lower_y_bound,upper_y_bound,ring_length
    !do in = start_i,end_i
    !   write(6,*) 'youts:',in,youts(in)
    !end do

    !write(6,*) 'Plasma data out axis:',lower_y_bound,upper_y_bound,ring_length
    !do in = 1,maxn
    !   write(6,*) 'sd:',in,sd(in)
    !end do
    

    ! Assign same values to spts (sol22) and s_axis (real version for searches)
    spts = sd
    s_axis = sd

    ! Calculate the midpoint bin - this is needed since axis opt 1 could have unequal numbers of cells on each side
    midn = iposq(dble(ring_length/2.0),sd,maxn)-1

    !write(0,*) 'AXIS setup:', midn,maxn/2,maxn,ring_length/2.0
    
  end subroutine set_plasma_data_axis



  subroutine assign_plasma(y_lim,ne_lim,te_lim,ti_lim,vb_lim,ef_lim,teg_lim,tig_lim)
    implicit none
    real :: y_lim,ne_lim,te_lim,ti_lim,vb_lim,ef_lim,teg_lim,tig_lim
    ! This routine interpolates the plasma solution to obtain the
    ! plasma values at the specific coordinate
    real :: y_shift,fact
    integer,external :: ipos
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


    write(6,'(a,20(1x,g12.5))') 'ASSIGN PLASMA:', y_lim,ne_lim,te_lim,ti_lim,vb_lim,ef_lim,teg_lim,tig_lim,lower_y_bound,y_shift,sd(1),sd(maxn)

    
  end subroutine assign_plasma
       

  subroutine calculate_tgrad_e
    ! this routine takes the plasma solution and calculates the temperature gradient
    ! and electric field
    ! the midpoint may have a discontinuity so these are set to zero
    implicit none

    integer :: iy
    real*8 :: ds1,ds2,dp1,dp2,dt1,dt2,nb1,nb2,e1,e2
    real*8 :: grad1,grad2

    ! proper scaling of teg and e for the simulation (QS is radial time step multiplier if using accelerated time in core or SOL)
    !TGSCAL = EMI/CRMI * QTIM_local *QTIM_local * QS(IQXS(IX)) * QS(IQXS(IX)

    !DSTEP = TGSCAL * QS(IQXS(IX)) * QS(IQXS(IX))
    !EFACT  = QTIM_local * QTIM_local * EMI / CRMI -> same as dstep

    ! scaling for the temperature gradient is included directly in the ctegs and ctigs arrays and so is applied
    ! directly in this routine. However scaling for the electric field is in the cfexzs array already so no scaling
    ! is applied to the efield in this routine (typically ef_scale=1.0)

    
    ! teg
    do iy = 1,maxn
       if (iy.eq.1) then
          grad1 = (te(iy+1)-te(iy))/(sd(iy+1)-sd(iy))
          teg(iy) = tg_scale * grad1
       elseif (iy.eq.maxn) then
          grad2 = (te(iy)-te(iy-1))/(sd(iy)-sd(iy-1))
          teg(iy) = tg_scale * grad2 
       elseif (iy.eq.int(maxn/2).or.iy.eq.int(maxn/2)+1) then
          ! set the midpoint of the solution to zero since it could be disjoint from each target
          teg(iy) = 0.0
       else
          grad1 = (te(iy+1)-te(iy))/(sd(iy+1)-sd(iy))
          grad2 = (te(iy)-te(iy-1))/(sd(iy)-sd(iy-1))
          teg(iy) = tg_scale * (grad1+grad2)/2.0 
       endif

    end do

    ! tig
    do iy = 1,maxn
       if (iy.eq.1) then
          grad1 = (ti(iy+1)-ti(iy))/(sd(iy+1)-sd(iy))
          tig(iy) = tg_scale * grad1
       elseif (iy.eq.maxn) then
          grad2 = (ti(iy)-ti(iy-1))/(sd(iy)-sd(iy-1))
          tig(iy) = tg_scale * grad2 
       elseif (iy.eq.int(maxn/2).or.iy.eq.int(maxn/2)+1) then
          ! set the midpoint of the solution to zero since it could be disjoint from each target
          tig(iy) = 0.0
       else
          grad1 = (ti(iy+1)-ti(iy))/(sd(iy+1)-sd(iy))
          grad2 = (ti(iy)-ti(iy-1))/(sd(iy)-sd(iy-1))
          tig(iy) = tg_scale * (grad1+grad2)/2.0 
       endif

    end do

    ! ef
    do iy=1,maxn
       if (iy.eq.1) then

          ds1 = sd(iy+1) - sd(iy)
          dp1 = ne(iy+1)*te(iy+1) - ne(iy)*te(iy)
          dt1 = te(iy+1) - te(iy)
          nb1 = 0.5*(ne(iy+1) + ne(iy))

          if (ds1.ne.0.and.nb1.ne.0) then 
             e1 = -(1.0/nb1)*dp1/ds1 - 0.71 * dt1/ds1
          else
             e1 = 0.0
          endif

          ef(iy) = ef_scale * e1

       elseif (iy.eq.maxn) then

          ds2 = sd(iy) - sd(iy-1)
          dp2 = ne(iy)*te(iy) - ne(iy-1)*te(iy-1)
          dt2 = te(iy) - te(iy-1)
          nb2 = 0.5*(ne(iy) + ne(iy-1))

          if (ds2.ne.0.and.nb2.ne.0) then 
             e2 = -(1.0/nb2)*dp2/ds2 - 0.71 * dt2/ds2
          else
             e2 = 0.0
          endif

          ef(iy) = ef_scale * e2

       elseif (iy.eq.int(maxn/2).or.iy.eq.int(maxn/2)+1) then
          ! set the midpoint of the solution to zero since it could be disjoint from each target
          ef(iy) = 0.0

       else

          ds1 = sd(iy+1) - sd(iy)
          dp1 = ne(iy+1)*te(iy+1) - ne(iy)*te(iy)
          dt1 = te(iy+1) - te(iy)
          nb1 = 0.5*(ne(iy+1) + ne(iy))

          if (ds1.ne.0.and.nb1.ne.0) then 
             e1 = -(1.0/nb1)*dp1/ds1 - 0.71 * dt1/ds1
          else
             e1 = 0.0
          endif

          ds2 = sd(iy) - sd(iy-1)
          dp2 = ne(iy)*te(iy) - ne(iy-1)*te(iy-1)
          dt2 = te(iy) - te(iy-1)
          nb2 = 0.5*(ne(iy) + ne(iy-1))

          if (ds2.ne.0.and.nb2.ne.0) then 
             e2 = -(1.0/nb2)*dp2/ds2 - 0.71 * dt2/ds2
          else
             e2 = 0.0
          endif

          ef(iy) = ef_scale * 0.5 * (e1+e2)

       endif

    end do

    !do iy = 1,maxn,100
    !   write(6,'(a,i8,20(1x,g12.5))') 'Sample grad,e:',iy,teg(iy),tig(iy),ef(iy)
    !end do
    
  end subroutine calculate_tgrad_e
  
  subroutine prt_plasma(pz,ix,x,debug_step)
    implicit none
    integer :: pz,ix,debug_step
    real :: x
    ! print the plasma solution found to unit 6 for debugging
    ! This will always print the first and last 30 entries and then every debug_step entries
    integer :: in
    integer :: debug_always = 20

    
    write(6,'(2(a,1x,i8),a,1x,g12.5)') 'Plasma solution for zone=',pz,' IX=',ix,' and radius X=',x,' MIDN=',midn,' MAXN=',maxn
    
    do in = 1,maxn
       if ((in.le.debug_always.or.in.ge.maxn-debug_always).or.(int(in/debug_step)*debug_step.eq.in)) then
          write(6,'(3i8,20(1x,g12.5))') in,ix,pz,x,sd(in),sd(in)+lower_y_bound,ne(in),te(in),ti(in),vb(in),ef(in),teg(in),tig(in),ga(in)
       endif
    end do
       
 
  end subroutine prt_plasma
  
  
  subroutine allocate_plasma_data(n)
    use mod_params
    use allocate_arrays
    implicit none
      integer :: n,ierr

      call allocate_array(sd,n,'Local s soledge',ierr)
      call allocate_array(spts,n,'Local s sol22',ierr)
      call allocate_array(s_axis,n,'Local s real for searching',ierr)
      call allocate_array(ne,n,'Local ne',ierr)
      call allocate_array(te,n,'Local te',ierr)
      call allocate_array(ti,n,'Local ti',ierr)
      call allocate_array(ef,n,'Local ef',ierr)
      call allocate_array(vb,n,'Local vb',ierr)
      call allocate_array(teg,n,'Local teg',ierr)
      call allocate_array(tig,n,'Local tig',ierr)
      call allocate_array(ga,n,'Local ga',ierr)

      !call allocate_array(dne,n,'Local dne',ierr)
      !call allocate_array(dte,n,'Local dte',ierr)
      !call allocate_array(dti,n,'Local dti',ierr)

    end subroutine allocate_plasma_data

    subroutine deallocate_plasma_data
      implicit none
      
      if (allocated(sd)) deallocate(sd)
      if (allocated(spts)) deallocate(spts)
      if (allocated(s_axis)) deallocate(s_axis)
      if (allocated(ne)) deallocate(ne)
      if (allocated(te)) deallocate(te)
      if (allocated(ti)) deallocate(ti)
      if (allocated(ef)) deallocate(ef)
      if (allocated(vb)) deallocate(vb)
      if (allocated(teg)) deallocate(teg)
      if (allocated(tig)) deallocate(tig)
      if (allocated(ga)) deallocate(ga)
      !if (allocated(dne)) deallocate(dne)
      !if (allocated(dte)) deallocate(dte)
      !if (allocated(dti)) deallocate(dti)

    end subroutine deallocate_plasma_data
   


end module mod_plasma_data
