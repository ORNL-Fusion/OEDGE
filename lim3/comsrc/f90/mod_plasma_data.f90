module mod_plasma_data

  integer :: maxn = 1000

  real*8,allocatable:: sd(:),spts(:),ne(:),te(:),ti(:),ef(:),vb(:),neg(:),teg(:),tig(:),ga(:)

  real*8 :: n_bnd_lower, te_bnd_lower, ti_bnd_lower, n_bnd_upper, te_bnd_upper, ti_bnd_upper

  real :: lower_y_bound, upper_y_bound , ring_length
  real :: s_axis(:)
  
  !integer :: data_count  ... maxn is used as a local in soledge - may need to use another variable name
  integer :: maxn = 1000 !initialize to 1000 but adjust depending on simulaton scale size
  ! alternative is to use actual y axes but that has its own challenges. 

  real*8,allocatble :: y_axis(:)

  
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
    implicit none
    real :: youts(-max_i:max_i),bnd1,bnd2
    integer :: start_i, end_i,axis_opt

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
       call allocate_array(y_axis,1,pt_cnt,'Local actual Y axis',ierr)
       
       do in = 1,pt_cnt
          y_axis(in) = (in-1) * ring_length/(pt_cnt-1) + bnd1
       end do


    elseif (axis_opt.eq.1) then 

       pt_cnt = 0
       do in = start_i, end_i
          if (youts(in).gt.bnd1.and.youts(in).lt.bnd2) then
             pt_cnt = pt_cnt+1
          endif
       end do
       ! Add 2 points for the end points bnd1,bnd2
       pt_cnt = pt_cnt+2

       ! allocate axis array
       call allocate_array(y_axis,1,pt_cnt,'Local actual Y axis',ierr)

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

    ! Assign same values to spts (sol22) and s_axis (real version for searches)
    spts = sd
    s_axis = sd

  end subroutine set_plasma_data_axis

  subroutine assign_plasma(y,n_out,te_out,ti_out,v_out,e_out,teg_out,tig_out,axis_opt)
    implicit none
    real :: y,n_out,te_out,ti_out,v_out,e_out,teg_out,tig_out
    integer :: axis_opt
    
    ! find Y value mapped into results arrays
    y_map = y-lower_y_bound

    ! interpolate axis value for y_map - this should match a y value almost exactly for axis_opt=0
    ! however, numerical rounding may result in a non-exact match - so use the interpolation code anyway
    ! slightly more computation but it will work for both cases.
    !
    
    ! check for values outside or at boundaries - remember sd=spts so it doesn't matter which is used here
    ! The only reason there are two named is because the axis names in soledge and sol22 are different. 
    if (y_map.le.sd(1)) then 
       
    elseif (y_map.ge.sd(maxn)) then


    else ! interpolate   
    


    endif


       
  end subroutine assign_plasma


  
  
  subroutine allocate_plasma_data(n)
    use mod_params
    use allocate_arrays
    implicit none
      integer :: n,ierr

      call allocate_array(sd,1,n,'Local s soledge',ierr)
      call allocate_array(spts,1,n,'Local s sol22',ierr)
      call allocate_array(ne,1,n,'Local ne',ierr)
      call allocate_array(te,1,n,'Local te',ierr)
      call allocate_array(ti,1,n,'Local ti',ierr)
      call allocate_array(ef,1,n,'Local ef',ierr)
      call allocate_array(vb,1,n,'Local vb',ierr)
      call allocate_array(dne,1,n,'Local dne',ierr)
      call allocate_array(dte,1,n,'Local dte',ierr)
      call allocate_array(dti,1,n,'Local dti',ierr)

    end subroutine allocate_plasma_data

    subroutine deallocate_plasma_data
      implicit none
      
      if (allocated(sd)) deallocate(sd)
      if (allocated(spts)) deallocate(spts)
      if (allocated(ne)) deallocate(ne)
      if (allocated(te)) deallocate(te)
      if (allocated(ti)) deallocate(ti)
      if (allocated(ef)) deallocate(ef)
      if (allocated(vb)) deallocate(vb)
      if (allocated(dne)) deallocate(dne)
      if (allocated(dte)) deallocate(dte)
      if (allocated(dti)) deallocate(dti)

    end subroutine deallocate_plasma_data
   


end module mod_plasma_data
