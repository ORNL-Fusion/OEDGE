module yreflection
  use error_handling
  implicit none

  integer :: yreflection_opt
  real :: cmir_refl_lower, cmir_refl_upper
  real*8 :: cmir_refl_lower_dp, cmir_refl_upper_dp

  real*8 :: relocation_count = 0

  real*8 :: lim_sep,deltay
  real*8 :: minrefl

  real*8 :: yref_upper_cnt = 0.0
  real*8 :: yref_lower_cnt = 0.0

  private :: lim_sep,cmir_refl_lower_dp,cmir_refl_upper_dp,minrefl
  save

  !
  ! Add variables to record reflection statistics
  !
  integer,parameter :: nxbins = 4
  integer,parameter :: none=0, lower = 1, upper = 2

  real*8 :: yreflection_event_count = 0
  real*8 :: x_wall_limit,x_central_limit,dx_inboard,dx_outboard
  real*8 :: refl_cnt(-nxbins:nxbins,2,2),refl_pos(-nxbins:nxbins,2,2)
  real*8 :: refl_cnt_first(-nxbins:nxbins,2,2),refl_pos_first(-nxbins:nxbins,2,2)

  real*8,parameter :: minval = 1.0e-10

contains

  subroutine init_reflection(ctwol,ca,caw)
    implicit none
    real :: ctwol,ca,caw

    ! initialize data
    yreflection_event_count = 0.0
    relocation_count = 0.0

    yref_upper_cnt = 0.0
    yref_lower_cnt = 0.0

    lim_sep = ctwol
    
    cmir_refl_lower_dp = cmir_refl_lower
    cmir_refl_upper_dp = cmir_refl_upper

    minrefl = 1.0d-7
    
    x_central_limit = ca
    x_wall_limit = caw
    dx_inboard = ca/nxbins
    dx_outboard = -caw/nxbins


    refl_cnt = 0.0
    refl_pos = 0.0
    refl_cnt_first = 0.0
    refl_pos_first = 0.0

    !
    ! Check to see if valid input has been specified - if not - turn the option off
    !

    if (yreflection_opt.ne.0.and.(cmir_refl_lower.ge.0.0.or.cmir_refl_upper.le.0.0)) then 
       call errmsg('YREFLECTION:TEST_REFLECTION','Y-REFLECTION OPTION IS ON BUT REFLECTION&
                   & LOCATIONS HAVE NOT BEEN PROPERLY SPECIFIED : Y-REFLECTION DISABLED')
       yreflection_opt=0
       return
    endif



  end subroutine init_reflection

  logical function check_reflected_region(y)
    implicit none
    real  ::  y,absy
    real*8 :: ynew

    ! check to see if the input y value lies in the region with mirrors on either side 
    ! If it does then the particle can't be launched here or does not belong here and 
    ! has leaked somehow.
    
    
    check_reflected_region=.false.

    !
    ! Check to see if the option is active - if it is not then just execute a return
    !

    if (yreflection_opt.eq.0) return

    !
    ! Check to see if the Y position is in a reflected region.
    !

    ynew = y

    if ((ynew.ge.-lim_sep+cmir_refl_upper_dp.and.ynew.le.cmir_refl_lower_dp).or.&
        (ynew.ge.cmir_refl_upper_dp.and.ynew.le.lim_sep+cmir_refl_lower_dp)) then
       check_reflected_region=.true.
    endif


  end function check_reflected_region


  subroutine check_reflection(x,y,oldy,svy,sputy,part_refl_cnt,part_type,debugl,ierr)
    implicit none

    !
    ! part_type indicates the type of particle: 1=neutral  2=ion
    !

    real :: x,y,svy,sputy
    real,intent(in) :: oldy

    real*8 :: part_refl_cnt
    integer,intent(in) :: part_type
    logical :: debugl
    integer :: ierr

    logical :: reflected
    !real :: y_org,oldy_org,svy_org
    real :: y_tmp

    real*8 :: ynew,yprev
    
    integer :: mirror,ix
    real*8 :: last_refl_part_cnt 


    real :: ran
    real,external :: getranf

    !
    ! Just exit if option not active
    !

    if (yreflection_opt.eq.0) return

    mirror = none
    ierr =0
    reflected = .false.

    !if (debugl) then 
    !   y_org = y
    !   oldy_org = oldy
    !   svy_org = svy
    !endif

    !
    ! jdemod - perform intermediate calculations in double precision
    !

    yprev = dble(oldy)
    ynew  = dble(y)

    if (ynew.lt.0) then 
       if ((yprev.gt.cmir_refl_lower_dp).and.(ynew.le.cmir_refl_lower_dp)) then 
          !
          !                   Reflection at lower (<0) mirror
          !

          deltay = (cmir_refl_lower_dp-ynew)
          ynew = cmir_refl_lower_dp + max(deltay,minrefl)

          svy = -svy
          reflected = .true.
          
          mirror = lower

       elseif ((yprev.lt.(-lim_sep+cmir_refl_upper_dp)).and.(ynew.ge.(-lim_sep+cmir_refl_upper_dp))) then 
          !
          !                   Reflection at second upper (< 0) mirror
          !

          deltay = (-lim_sep+cmir_refl_upper_dp-ynew)
          ynew = -lim_sep+cmir_refl_upper_dp +  min(deltay,-minrefl)

          svy = -svy
          reflected = .true.

          mirror = upper

       endif

    elseif (ynew.ge.0) then 
       if ((yprev.lt.cmir_refl_upper_dp).and.(ynew.ge.cmir_refl_upper_dp)) then 
          !
          !                   Reflection at upper (>0) mirror
          !
          deltay = (cmir_refl_upper_dp - ynew)
          ynew = cmir_refl_upper_dp + min(deltay,-minrefl)

          svy = -svy
          reflected = .true.

          mirror = upper

       elseif ((yprev.gt.(lim_sep+cmir_refl_lower_dp)).and.(ynew.le.(lim_sep+cmir_refl_lower_dp))) then 
          !
          !                   Reflection at second lower (>0) mirror
          !
          deltay = (lim_sep +cmir_refl_lower_dp-ynew)
          ynew = lim_sep + cmir_refl_lower_dp + max(deltay,minrefl)

          svy = -svy
          reflected = .true.

          mirror = lower

       endif
    endif

    !
    ! Record statistics if reflection occurred
    !

    if (reflected) then

          yreflection_event_count = yreflection_event_count + 1.0
          part_refl_cnt = part_refl_cnt + 1.0 

          if (x.gt.0.0) then 
             ix = min(int (x/dx_inboard),nxbins)
          elseif (x.le.0.0) then 
             ix = max(int(x/dx_outboard),-nxbins)
          endif

          if (ix.lt.-nxbins.or.ix.gt.nxbins.or.&
             (mirror.ne.1.and.mirror.ne.2).or.&
             (part_type.ne.1.and.part_type.ne.2)) then 

              write(error_message_data,'(a,3i10)') 'DATA=',ix,mirror,part_type
              call errmsg('YREFLECTION:',error_message_data)

          else

             refl_cnt(ix,mirror,part_type) = refl_cnt(ix,mirror,part_type) + sputy
             refl_pos(ix,mirror,part_type) = refl_pos(ix,mirror,part_type) + x * sputy 

             if (part_refl_cnt.eq.1.0) then 
                refl_cnt_first(ix,mirror,part_type) = refl_cnt_first(ix,mirror,part_type) + sputy
                refl_pos_first(ix,mirror,part_type) = refl_pos_first(ix,mirror,part_type) + x *sputy
             endif 

         endif

    endif 


    y = sngl(ynew)

    if (debugl.and.reflected) then
       write(error_message_data,'(a,i10,6(1x,g18.10))') 'REFLECTION OCCURRED:',int(yreflection_event_count),y,oldy,svy
       call dbgmsg('CHECK_REFLECTION',error_message_data)
    endif

    if (check_reflected_region(y)) then 
       
       write(error_message_data,'(a,4(1x,g18.10))') 'REFLECTED PARTICLE HAS ENTERED MIRROR REGION - '//&
                                                  & 'TRY REDUCING SIMULATION TIMESTEPS AND MULTIPLIERS : DATA:', x,y,oldy,svy
       call errmsg('CHECK REFLECTION:WARNING:',error_message_data)
       
       y_tmp = y
       ran = getranf()
       y = cmir_refl_lower + ran * (cmir_refl_upper-cmir_refl_lower)

       write(error_message_data,'(a,i10,l4,10(1x,g20.12))') 'REFLECTION:',int(yreflection_event_count),reflected,x,deltay,ynew,oldy,ynew,yprev,cmir_refl_lower_dp,cmir_refl_upper_dp
       call dbgmsg('CHECK_REFLECTION',error_message_data)

       relocation_count = relocation_count +1.0

       write(error_message_data,'(a,2i12,10(1x,g20.12))') 'REVISED Y COORDINATES:',int(relocation_count),int(yreflection_event_count),oldy,y_tmp,y,svy
       call errmsg('CHECK_REFLECTION PARTICLE RELOCATED:',error_message_data)


       ierr =1
    endif


  end subroutine check_reflection



  subroutine pr_yref_stats(tot_part,tot_refl_part,max_part_refl_cnt,tot_part_refl_cnt)
    implicit none

    real*8 :: tot_part,tot_refl_part,max_part_refl_cnt,tot_part_refl_cnt
    
    character*512 :: output_message
    integer in

    real :: xmin, xmax
!
! Print some statistics on yreflection events
! 
    if (yreflection_opt.ne.0) then 

        call prb
        call prc('  Summary of Y-Reflection events:')
        call prr('  Total number of reflection events  :',sngl(yreflection_event_count))
        call prr('  Total weight of particles reflected:',sngl(sum(refl_cnt(:,:,:))))
        call prr('  Total particles launched           :',sngl(tot_part))
        call prr('  Total particles reflected          :',sngl(tot_refl_part))

        call prr('  Average number of reflection/particle          :',sngl(tot_part_refl_cnt/max(tot_part,minval)))
        call prr('  Average number of reflection/reflected particle:',sngl(tot_part_refl_cnt/max(tot_refl_part,minval)))
        call prr('  Maximum number of reflections for any particle :',sngl(max_part_refl_cnt))

        call prb
        call prc('  Lower Mirror Summary:')
        call prr('  Lower Mirror Location: ',cmir_refl_lower)
        call prr('     - Total Reflections from lower mirror:',sngl(sum(refl_cnt(:,lower,:))))
        call prr('     - Total Ion Reflection Count         :',sngl(sum(refl_cnt(:,lower,2))))
        call prr('     - Total Neutral Reflection Count     :',sngl(sum(refl_cnt(:,lower,1))))
        call prr('     - Total Ion Particles reflecting     :',sngl(sum(refl_cnt_first(:,lower,2))))
        call prr('     - Total Neutral particles reflecting :',sngl(sum(refl_cnt_first(:,lower,1))))
        call prr('     - Average Outboard location of first reflection:', sngl( sum(refl_pos_first(-nxbins:-1,lower,:))/max (sum (refl_cnt_first (-nxbins:-1,lower,:)),minval)))
        call prr('     - Average Inboard  location of first reflection:', sngl( sum(refl_pos_first(0:nxbins,lower,:))/max (sum (refl_cnt_first (0:nxbins,lower,:)),minval)))
        call prc('  Lower mirror Reflection data:')
        call prb

        write(output_message,'(2x,2(a13),6(a21))') 'Bin LBound','Bin Ubound','Particle Count','Average First X','Total Reflections','Average X (all)'
        call prc(output_message)
 

        do in = -nxbins,nxbins-1

          if (in.lt.0) then 
             xmin = in*dx_outboard
             xmax = (in+1) * dx_outboard
          else
             xmin = in * dx_inboard
             xmax = (in+1) * dx_inboard
          endif
          
          write(output_message,'(2x,2(1x,f12.4),6(3x,g18.6))') xmin,xmax,&
                                                    sngl(sum(refl_cnt_first(in,lower,:))),&
                                                    sngl(sum(refl_pos_first(in,lower,:))/max(sum(refl_cnt_first(in,lower,:)),minval)),  &
                                                    sngl(sum(refl_cnt(in,lower,:))),&
                                                    sngl(sum(refl_pos(in,lower,:))/max(sum(refl_cnt(in,lower,:)),minval))
           call prc(output_message)
        end do 


        call prb
        call prc('  Upper Mirror Summary:')
        call prr('  Upper Mirror Location: ',cmir_refl_upper)
        call prr('     - Total Reflections from upper mirror:',sngl(sum(refl_cnt(:,upper,:))))
        call prr('     - Total Ion Reflection Count         :',sngl(sum(refl_cnt(:,upper,2))))
        call prr('     - Total Neutral Reflection Count     :',sngl(sum(refl_cnt(:,upper,1))))
        call prr('     - Total Ion Particles reflecting     :',sngl(sum(refl_cnt_first(:,upper,2))))
        call prr('     - Total Neutral particles reflecting :',sngl(sum(refl_cnt_first(:,upper,1))))
        call prr('     - Average Outboard location of first reflection:', sngl( sum(refl_pos_first(-nxbins:-1,upper,:))/max (sum (refl_cnt_first (-nxbins:-1,upper,:)),minval)))
        call prr('     - Average Inboard  location of first reflection:', sngl( sum(refl_pos_first(0:nxbins,upper,:))/max (sum (refl_cnt_first (0:nxbins,upper,:)),minval)))
        call prc('  Upper mirror Reflection data:')
        call prb

        write(output_message,'(2x,2(a13),6(a21))') 'Bin LBound','Bin Ubound','Particle Count','Average First X','Total Reflections','Average X (all)'
        call prc(output_message)
 
       do in = -nxbins,nxbins-1

          if (in.lt.0) then 
             xmin = in*dx_outboard
             xmax = (in+1) * dx_outboard
          else
             xmin = in * dx_inboard
             xmax = (in+1) * dx_inboard
          endif
          

          write(output_message,'(2x,2(1x,f12.4),6(3x,g18.6))') xmin,xmax,&
                                                    sngl(sum(refl_cnt_first(in,upper,:))),&
                                                    sngl(sum(refl_pos_first(in,upper,:))/max(sum(refl_cnt_first(in,upper,:)),minval)),  &
                                                    sngl(sum(refl_cnt(in,upper,:))),&
                                                    sngl(sum(refl_pos(in,upper,:))/max(sum(refl_cnt(in,upper,:)),minval))
           call prc(output_message)
        end do 

        call prb
        call prr('  Total number of relocated particles: ',sngl(relocation_count))

    endif

  end subroutine pr_yref_stats




end module yreflection
