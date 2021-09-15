module mod_vtig


  implicit none


  private

  integer,public :: n_vtig_blocks, n_vb_blocks
  integer,public :: vtig_opt = 0
  integer,public :: vb_opt = 0
  
  integer,public :: max_data = 10
  
  real,allocatable,public :: vtig_range(:,:)
  real,allocatable,public :: vtig_data(:,:,:)
  integer,allocatable,public :: vtig_ndata(:)
  integer,allocatable,public :: vtig_zones(:,:)

  real,allocatable,public :: vb_range(:,:)
  real,allocatable,public :: vb_data(:,:,:)
  integer,allocatable,public :: vb_ndata(:)
  integer,allocatable,public :: vb_zones(:,:)



  public :: read_v_data,calculate_temperature,setup_vtig,allocate_v_data,deallocate_v_data,calculate_velocity,vtig_tgrad

  character*128 :: buffer
  character*6 :: buff_fmt = '(a128)'

  real,public :: integration_const
  

contains

  subroutine setup_vtig(crmb,crmi,cnbin,ctibin)
    use mod_params
    implicit none
    ! crmi and crmb are part of comtor so they are no longer needed as input parameters
    real :: crmi,crmb,cnbin,ctibin
    
    ! assume ln(lambda) ~= 15.0    - 15.0 is not used in LIM
    ! note - density is also needed
    ! note - Bi ~= 2.6 Z^2 has been used to simplify the expression which is correct as U -> 1  (U = mz/(mi+mz))
    ! mi = crmb, mz = crmi
    ! note mz in Tau_s cancels the mz in the pre-factor for the vTiG expression
    ! integration constant = e Tau_s Bi/mz Ti^3/2     (  vTiG(s) = e Tau_s Bi/mz  dT/ds)
    ! Tau_s = (1.47e13 mz Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni Z^2 ln_Lam )
    ! mz = cmri, mi = crmb
    ! Bi = 2.6 Z^2
    !
    ! integration_constant = e/(amu mz)  (1.47e13 mz Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni Z^2 ln_Lam ) * 2.6 Z^2
    ! cancel mz term and Z^2 term
    ! integration_constant = e/amu  (1.47e13 Ti (Ti/mi)^(1/2)) / ( (1+mi/mz) ni ln_Lam ) * 2.6
    ! pull out Ti^3/2
    ! integration_constant = e/amu  (1.47e13 (1/mi)^(1/2)) / ( (1+mi/mz) ni ln_Lam ) * 2.6    Ti^3/2
    ! pull out n since it isn't known until later
    ! integration_constant = e/amu  (1.47e13 (1/mi)^(1/2)) / ( (1+mi/mz) ln_Lam ) * 2.6    /n     Ti^3/2
    ! integration_constant = e/amu  (1.47e13 (1/crmb)^(1/2)) / ( (1+crmb/crmi) ln_Lam ) * 2.6 
   
    
    !real,parameter :: ln_lambda = 15.0
    real :: ln_lambda
    
    !
    ! LIM uses the following definition of LAMBDA rather than the fixed value of 15. However, CNBIN and CTIBIN are
    ! the inboard base plasma values at the limiter tip so LAMBDA is still fixed for the simulation and doesn't vary
    ! spatially - its value just depends on the values assigned to CNBIN and CTIBIN.
    !
    ! For now, I will change the calculation here but it may be desirable to switch to using a fixed value of 15. 
    !
    ln_LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)

    integration_const = ech/amu * 1.47e13 * sqrt(1.0/crmb) * 2.6 / ( ( 1.0 + crmb/crmi) * ln_lambda)


  end subroutine setup_vtig

  subroutine calculate_temperature(x,y,pz,yz,n,t0,cl,t,nblocks,xrange,ndata,data,zones)
    implicit none

    ! cl = connection length - for scaling the integral
    integer :: nblocks,pz,yz
    real :: x,y,t0,t,n,cl
    integer :: ndata(nblocks)
    integer :: zones(nblocks,2)
    real :: xrange(nblocks,2)
    real :: data(nblocks,max_data,3)
    real :: iv
    integer :: ib,in,ia

    ib = find_block(x,nblocks,xrange,pz,yz,zones)
    
    ! if vTiG is not applied to this X range then return without changing anything - or density is invalid
    if (ib.eq.0.or.n.eq.0.0) return

    ! calculate the corresponding value of ti from integrating the equation.
    ! note - this integral is not correct since it really should be calculated in a way compatible with a spatially
    ! varying value of n - this is fine for constant n. 

    !iv = integral_value(y/cl,ndata(ib),data(ib,:,:))
    
    t = (t0**(5.0/2.0) + 5.0/(2.0*integration_const/n) * integral_value(y/cl,ndata(ib),data(ib,:,:))*cl)**(2.0/5.0)

    !write(6,'(a,3i8,20(1x,g12.5))') 'CT:',ib,pz,yz,x,y,t0,t,n,integration_const,cl,iv,iv*cl,integration_const/n

    !do in = 1,ndata(ib)
    !   write(6,'(a,3i8,20(1x,g12.5))') 'D:',in,ib,nblocks,(data(ib,in,ia),ia=1,3)
    !end do
    
  end subroutine calculate_temperature


  subroutine calculate_velocity(x,y,pz,yz,cl,vel,nblocks,xrange,ndata,data,zones)
    implicit none

    ! cl = connection length - for scaling the integral
    integer :: nblocks,pz,yz
    real :: x,y,cl,vel
    integer :: ndata(nblocks)
    integer :: zones(nblocks,2)
    real :: xrange(nblocks,2)
    real :: data(nblocks,max_data,3)

    integer :: ib

    ib = find_block(x,nblocks,xrange,pz,yz,zones)

    !write(6,'(a,4i8,20(1x,g12.5))') 'CV1:',ib,nblocks,pz,yz,x,y

    
    ! if v is not applied to this X range then return without changing anything
    if (ib.eq.0) return

    vel = real(yz) * interpolate_value(y/cl,ndata(ib),data(ib,:,:),2)

    !write(6,'(a,4i8,20(1x,g12.5))') 'CV2:',ib,nblocks,pz,yz,x,y,vel

  end subroutine calculate_velocity

  real function vtig_tgrad(ix,iy,pz)
    use mod_comtor
    use mod_comxyt
    use mod_comt2
    implicit none
    integer :: ix,iy,pz
    ! This routine assumes that the vtig data entered is NOT velocity data but dTi/ds which
    ! will be directly imposed in the calculation of the temperature gradient forces 

    integer :: iqx,yz
    real :: x,y,y0,yt
    real :: tgrad
    
    iqx = iqxs(ix)
    x = xouts(ix)
    y = youts(iy)

    ! the following gives "+" for -2CL < Y < -CL and 0 < Y < CL
    ! the following gives "-" for -CL < Y < 0 and CL < Y < 2CL
    if (youts(iy).lt.0.0) then
       yz = -int(sign(1.0,youts(iy)+cl))
       if (youts(iy).lt.-cl) then
          yt = 2.0*cl + youts(iy)
       else
          yt = abs(youts(iy))
       endif
    else
       yz = -int(sign(1.0,youts(iy)-cl))
       if (youts(iy).lt.cl) then
          yt = youts(iy)
       else
          yt = 2.0* cl - youts(iy)
       endif
    endif
       
    if (yz.lt.0.0) then 
       y0 = qedges(iqx,1)
    else
       y0 = qedges(iqx,2)
    endif
    
    call calculate_velocity(x,yt-y0,pz,yz,cl-y0,tgrad,n_vtig_blocks,vtig_range,vtig_ndata,vtig_data,vtig_zones)

    vtig_tgrad = tgrad

    write(6,'(a,3i8,20(1x,g12.5))') 'Tgrad:',ix,iy,yz,x,y,tgrad
    
  end function vtig_tgrad



  
  integer function find_block(x,nblocks,xrange,pz,yz,zones)
    implicit none
    real :: x
    integer :: nblocks,pz,yz
    integer :: zones(nblocks,2)
    real :: xrange(nblocks,2)

    ! local 
    integer :: ix

    ! find  - the data block of the profile for this X value - zero if none
    find_block = 0
         
    do ix = 1,nblocks
       ! check for both xrange and poloidal zone
       if (x.ge.xrange(ix,1).and.x.lt.xrange(ix,2).and.(pz.eq.zones(ix,1).or.zones(ix,1).eq.0).and.(yz.eq.zones(ix,2).or.zones(ix,2).eq.0)) then
          find_block = ix
          return
       endif
    end do
  
  end function find_block
  
  subroutine read_v_data(nblocks,xrange,ndata,data,zones)
    use mod_io_units
    implicit none

    integer :: nblocks
    real, allocatable :: xrange(:,:),data(:,:,:)
    integer,allocatable :: ndata(:),zones(:,:)

    ! local
    real temp_profile(max_data,3)
    integer :: in,ix,cnt
    
    !read in data blocks each is of the form 
    !     X1   X2   N  P
    !     X1, X2 are the radial xrange for the profile
    !     N is the number of points in the profile ... only the first 10 points can be read
    !     profiles are lineary interpolated
    
    if (nblocks.gt.0) then
       call allocate_v_data(nblocks,xrange,ndata,data,zones)
       ! read in profile data
       do in = 1,nblocks
10        read(stdin,buff_fmt) buffer
          if (buffer(1:1).eq.'$') goto 10
          !write(0,*) ':',trim(buffer),':'
          read(buffer,*) xrange(in,1),xrange(in,2),ndata(in),zones(in,1),zones(in,2)

          if (ndata(in).gt.max_data-2) then
             ! error - maximum 10 points including end points
             call errmsg('MOD_VTIG:READ_V_DATA:V PROFILES HAVE A MAXIMUM OF MAX_DATA-2 DATA POINTS: LIMIT=',max_data-2)
             stop 'MOD_VTIG: ERROR READING V PROFILE'
          endif

          do ix = 1,ndata(in)
             ! read in profile data
20           read(stdin,buff_fmt) buffer
             if (buffer(1:1).eq.'$') goto 20                     
             !write(0,*) ':',trim(buffer),':'
             read(buffer,*) data(in,ix,1),data(in,ix,2)
          end do
       end do

       ! modify profile to add points at S=0 and S=1 if not present

       do in = 1,nblocks
          
          cnt = 0 
          temp_profile = 0.0
          
          if (data(in,1,1).ne.0.0) then
             cnt = cnt+1
             temp_profile(cnt,1) = 0.0
             temp_profile(cnt,2) = 0.0
          endif

          do ix = 1,ndata(in)
             if ((data(in,ix,1).lt.0.0).or.(data(in,ix,1).gt.1.0)) then
                call errmsg('MOD_VTIG: READ_DATA:','INPUT ERROR: Spatial scale MUST be in the range [0,1]')
                stop 'MOD_VTIG:ERROR READING V DATA'
             endif
             cnt = cnt+1
             temp_profile(cnt,1) = data(in,ix,1)
             temp_profile(cnt,2) = data(in,ix,2)
          end do

          if (data(in,ndata(in),1).ne.1.0) then
             cnt = cnt+1
             temp_profile(cnt,1) = 1.0
             temp_profile(cnt,2) = 0.0
          endif

          ndata(in) = cnt
          data(in,:,:) = temp_profile

          ! calculate integral of profile
          call integrate_profile(ndata(in),data(in,:,:))

       end do

    endif

    !write(0,*) 'READ V DATA:',nblocks
    !do in = 1,nblocks
    !   write(0,*) 'Block:',in,xrange(in,1),xrange(in,2), ndata(in), zones(in,1),zones(in,2)
    !   do ix = 1,ndata(in)
    !      write(0,*) 'Data :', data(in,ix,1), data(in,ix,2),data(in,ix,3)
    !   end do
    !end do


    
  end subroutine read_v_data



  subroutine integrate_profile(ndata,data)
    implicit none
    integer :: ndata
    real :: data(max_data,3)
    integer :: in
    real :: y
    !
    ! Lengths are scaled to L :  0.0 is the limiter surface and 1.0 is L - connection length is 2L
    ! 
    ! Equation needs the integral of vTig ds where s is the parallel to the field line distance
    ! 

    ! calculate the integral of the profile and store it at each point, Integral value at zero is zero
    
    ! first data point is at 0, last at 1 - integral will be scaled by L in later calculations

    data(1,3) = 0.0
    
    do in = 1,ndata-1
       data(in+1,3) = data(in,3) + (data(in,2)+data(in+1,2))/2.0 * (data(in+1,1)-data(in,1))
    end do


    !do in = 1,ndata
    !   write(6,*) 'Data:',in,data(in,1),data(in,2),data(in,3)
    !end do
    
    !do in = 1,100
    !   y = real(in-1)/100.0
    !   write(6,'(a,i8,10(1x,g12.5))') 'Integral:',in,y,integral_value(y,ndata,data)
    !end do
    
  end subroutine integrate_profile
    
    
  real function integral_value(y,ndata,data)
    implicit none
    real :: y
    integer :: ndata
    real :: data(max_data,3)
    integer :: in
    integer,external :: ipos

    if (y.le.0) then
       integral_value = 0.0
       return
    endif
    
    in = ipos(y,data(:,1),ndata)

    ! the lowest bound in the array is at S=0.0
    ! as a result, the above should only find in = 1 for S = 0

    if (in.eq.1) then
       integral_value = data(in,3)
    else    
       ! for in = 2 .. data(in-1,3) should be zero 
       integral_value = (interpolate_value(y,ndata,data,2)+data(in-1,2))/2.0 * (y-data(in-1,1)) + data(in-1,3)
    endif
       
  end function integral_value


  real function interpolate_value(y,ndata,data,index)
    implicit none
    real :: y
    integer :: ndata,index
    real :: data(max_data,3)
    integer :: in,ia,ib
    integer,external :: ipos

    in = ipos(y,data(:,1),ndata)

    ! this should not happen since the first value is S=0 ... but if it does then return the integral of the first element
    if (in.eq.1) then 
       interpolate_value = data(1,index)
    else
       interpolate_value = (y-data(in-1,1)) / (data(in,1)-data(in-1,1)) * (data(in,index)-data(in-1,index)) + data(in-1,index)
    endif

    
    !write(6,'(a,2i8,20(1x,g12.5))') 'IV:',in,index,y,(data(in-1,ia),ia=1,3),(data(in,ib),ib=1,3),interpolate_value
    
    
    
  end function interpolate_value


  
  
  subroutine allocate_v_data(nblocks,xrange,ndata,data,zones)
    use allocate_arrays
    implicit none
    integer :: nblocks
    real, allocatable :: xrange(:,:),data(:,:,:)
    integer,allocatable :: ndata(:),zones(:,:)
    ! local
    integer:: ierr
    
    call allocate_array(ndata,nblocks,'Allocate ndata: Number of points in each profile - max 10 points',ierr)
    call allocate_array(zones,nblocks,2,'Allocate zone: Poloidal zone, and Y zone',ierr)
    call allocate_array(data,nblocks,max_data,3,'Allocate data: v profile data - max 10 ',ierr)
    call allocate_array(xrange,nblocks,2,'Allocate range: Xrange data - X1,X2',ierr)

  end subroutine allocate_v_data

  
  subroutine deallocate_v_data(xrange,ndata,data,zones)
    implicit none
    real, allocatable :: xrange(:,:),data(:,:,:)
    integer,allocatable :: ndata(:),zones(:,:)

    if (allocated(ndata)) deallocate(ndata)
    if (allocated(zones)) deallocate(zones)
    if (allocated(data)) deallocate(data)
    if (allocated(xrange)) deallocate(xrange)

  end subroutine deallocate_v_data


end module mod_vtig
