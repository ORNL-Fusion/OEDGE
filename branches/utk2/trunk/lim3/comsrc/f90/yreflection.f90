module yreflection
  use error_handling
  !  use mod_comtor ! sazmod
  use variable_wall
  implicit none

  !
  ! Absorbtion surfaces options
  !
  integer :: xabsorb_opt, yabsorb_opt

  integer :: yabsorb1_frame, yabsorb2_frame
  real :: yabsorb1a, yabsorb2a, yabsorb1b, yabsorb2b,xabsorb
  !
  !    sazmod - switch to vary absorbing boundary that affects plasma solution
  !             and the accompanying values.
  real :: yabsorb1a_step, yabsorb2a_step, xabsorb1a_step, xabsorb2a_step
  integer, public:: vary_absorb, ix_step1, ix_step2

  !
  ! absorption statistics
  !

  real*8 :: yabsorb1_cnt, yabsorb2_cnt, xabsorb_cnt
  real*8 :: xabsorb_sputy, xabsorb_neut, xabsorb_ion, xabsorb_iz, xabsorb_yavg
  real*8 :: yabsorb1_sputy, yabsorb1_neut, yabsorb1_ion, yabsorb1_iz, yabsorb1_xavg
  real*8 :: yabsorb2_sputy, yabsorb2_neut, yabsorb2_ion, yabsorb2_iz, yabsorb2_xavg


  !
  ! Reflection options
  !
  integer :: yreflection_opt
  real :: cmir_refl_lower, cmir_refl_upper
  real*8 :: cmir_refl_lower_dp, cmir_refl_upper_dp

  real*8 :: relocation_count = 0

  real*8 :: lim_sep  !,deltay
  real*8 :: minrefl

  real*8 :: yref_upper_cnt = 0.0
  real*8 :: yref_lower_cnt = 0.0

  integer :: xreflection_opt
  real :: xreflect_bound

  integer :: preflect_opt
  real :: preflect_bound

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

  real*8,parameter,private :: minval = 1.0e-10


  !
  !     jdemod - add variables for keeping track of individual particle reflections
  !
  real*8 :: part_refl_cnt
  integer :: part_frame
  integer :: part_mirror_state

  integer,parameter :: nframes = 11
  real*8 :: frame_cnt(-nframes:nframes+1)


  real*8 :: max_part_refl_cnt
  real*8 :: tot_part_refl_cnt
  real*8 :: tot_part
  real*8 :: tot_refl_part


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

    ! Can not use a fixed value for minimum reflection distance since the scale size of the system and thus rounding into single precision also change. 
    !minrefl = 1.0d-7
    minrefl = dble(lim_sep)/1.0d7

    x_central_limit = ca
    x_wall_limit = caw
    dx_inboard = ca/nxbins
    dx_outboard = -caw/nxbins


    refl_cnt = 0.0
    refl_pos = 0.0
    refl_cnt_first = 0.0
    refl_pos_first = 0.0



    !
    !     Initialize reflection counters
    !
    max_part_refl_cnt = 0.0
    tot_part_refl_cnt = 0.0
    tot_part = 0.0
    tot_refl_part = 0.0

    !
    !     Initialize the absorbtion counters
    !

    yabsorb1_cnt = 0.0
    yabsorb2_cnt = 0.0
    xabsorb_cnt = 0.0

    !
    !      Initialize absorption statistics      
    !

    xabsorb_sputy = 0.0
    xabsorb_neut  = 0.0
    xabsorb_ion   = 0.0
    xabsorb_iz    = 0.0
    xabsorb_yavg  = 0.0

    yabsorb1_sputy= 0.0
    yabsorb1_neut = 0.0
    yabsorb1_ion  = 0.0 
    yabsorb1_iz   = 0.0
    yabsorb1_xavg = 0.0

    yabsorb2_sputy= 0.0 
    yabsorb2_neut = 0.0
    yabsorb2_ion  = 0.0 
    yabsorb2_iz   = 0.0
    yabsorb2_xavg = 0.0

    !
    !     Set up the secondary absorbing surfaces outside the -L to L range
    !
    if (yabsorb1a.gt.0.0) then 
       yabsorb1b = yabsorb1a-lim_sep
    else
       yabsorb1b = yabsorb1a+lim_sep
    endif

    if (yabsorb2a.gt.0.0) then 
       yabsorb2b = yabsorb2a-lim_sep
    else
       yabsorb2b = yabsorb2a+lim_sep
    endif
	write(6,*)'yabsorb1b = ',yabsorb1b
	write(6,*)'yabsorb2b = ',yabsorb2b


    !
    ! Check to see if valid input has been specified - if not - turn the option off
    !

    if (yreflection_opt.ne.0.and.(cmir_refl_lower.ge.0.0.or.cmir_refl_upper.le.0.0)) then 
       call errmsg('YREFLECTION:TEST_REFLECTION','Y-REFLECTION OPTION IS ON BUT REFLECTION&
            & LOCATIONS HAVE NOT BEEN PROPERLY SPECIFIED : Y-REFLECTION DISABLED')
       yreflection_opt=0
    endif



  end subroutine init_reflection


  subroutine init_part_reflection
    implicit none

    !
    !       jdemod - initialize particle reflection counter
    !     
    !
    !       Note: Frame data and mirror_state is not carried over between ions and neutrals at present. 
    !             This neutrals that start in a different frame will be reverted to the 0 frame when 
    !             started as ions. This effect should not be significant as long as the number of reflected
    !             neutrals is relatively low. If this is not the case then the frame and mirror_state may 
    !             need to be retained for each particle (neutral -> ion) 
    !             At the present time neutrals are not included in the summary data.
    !

    part_refl_cnt = 0.0
    part_frame = 0
    part_mirror_state = 1

  end subroutine init_part_reflection


  subroutine update_part_refl_stats(sputy)
    implicit none
    real :: sputy
    integer :: in
    !
    !       Record some reflection statistics
    !

    if (yreflection_opt.ne.0) then 


       max_part_refl_cnt = max(part_refl_cnt,max_part_refl_cnt)

       tot_part_refl_cnt = tot_part_refl_cnt + part_refl_cnt * sputy

       tot_part = tot_part + sputy

       if (part_refl_cnt.gt.0.0) then 
          tot_refl_part = tot_refl_part + sputy
       endif

       in =   max(min(part_frame,nframes+1),-nframes)
       frame_cnt(in) = frame_cnt(in) + sputy


       !write(0,'(a,3i10,g12.5)') 'FRAME:',in,part_frame,frame_cnt(in),sputy



    endif

  end subroutine update_part_refl_stats



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

    ! Change reflected region test to not include particles exactly on the boundaries (gt,lt and not ge,le)

    if ((ynew.gt.-lim_sep+cmir_refl_upper_dp.and.ynew.lt.cmir_refl_lower_dp).or.&
         (ynew.gt.cmir_refl_upper_dp.and.ynew.lt.lim_sep+cmir_refl_lower_dp)) then
       check_reflected_region=.true.
    endif


  end function check_reflected_region

  subroutine check_reflection(x,y,oldy,svy,sputy,part_type,debugl,ierr)
    implicit none

    !
    ! part_type indicates the type of particle: 1=neutral  2=ion
    !

    real :: x,y,svy,sputy
    real,intent(in) :: oldy

    integer,intent(in) :: part_type

    logical :: debugl
    integer :: ierr

    logical :: reflected
    !real :: y_org,oldy_org,svy_org
    real :: y_tmp

    real*8 :: ynew,yprev,yorg,ytest

    integer :: mirror,ix
    real*8 :: last_refl_part_cnt 
    real*8 :: deltay


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
    yorg  = dble(y)
    ynew  = dble(y)
    deltay = 0.0d0


    !
    ! jdemod - change code to allow for a second reflection option to permit
    !          reflection in either direction from the mirrors. 
    !

    !--------------  Option 1 ---------------------------------

    if (yreflection_opt.eq.1) then 

       if (yorg.lt.0.0) then 
          if ((yprev.gt.cmir_refl_lower_dp).and.(yorg.le.cmir_refl_lower_dp)) then 
             !
             !                   Reflection at lower (<0) mirror
             !

             deltay = (cmir_refl_lower_dp-yorg)
             ynew = cmir_refl_lower_dp + max(deltay,minrefl)

             svy = -svy
             reflected = .true.

             mirror = lower

          elseif ((yprev.lt.(-lim_sep+cmir_refl_upper_dp)).and.(yorg.ge.(-lim_sep+cmir_refl_upper_dp))) then 
             !
             !                   Reflection at second upper (< 0) mirror
             !

             deltay = (-lim_sep+cmir_refl_upper_dp-yorg)
             ynew = -lim_sep+cmir_refl_upper_dp +  min(deltay,-minrefl)

             svy = -svy
             reflected = .true.

             mirror = upper

          endif

       elseif (yorg.ge.0.0) then 
          if ((yprev.lt.cmir_refl_upper_dp).and.(yorg.ge.cmir_refl_upper_dp)) then 
             !
             !                   Reflection at upper (>0) mirror
             !
             deltay = (cmir_refl_upper_dp - yorg)
             ynew = cmir_refl_upper_dp + min(deltay,-minrefl)

             svy = -svy
             reflected = .true.

             mirror = upper

          elseif ((yprev.gt.(lim_sep+cmir_refl_lower_dp)).and.(yorg.le.(lim_sep+cmir_refl_lower_dp))) then 
             !
             !                   Reflection at second lower (>0) mirror
             !
             deltay = (lim_sep +cmir_refl_lower_dp-yorg)
             ynew = lim_sep + cmir_refl_lower_dp + max(deltay,minrefl)

             svy = -svy
             reflected = .true.

             mirror = lower

          endif
       endif


       !--------------  Option 2 ----------------------------------------------------------

    elseif (yreflection_opt.eq.2) then


       ! Need to check for unexpectedly large particle steps - code below should work fine
       ! when the particle doesn't step from one side of the limiter to the other. 



       if (yprev.lt.0.0) then 

          if (abs(yorg-yprev).gt.lim_sep/2.0) then
             ytest = yorg - lim_sep
          else
             ytest = yorg
          endif
          
          
          if ((yprev.gt.cmir_refl_lower_dp).and.(ytest.le.cmir_refl_lower_dp)) then 
             !
             !                   Reflection at lower (<0) mirror
             !

             deltay = max(abs(cmir_refl_lower_dp-ytest),minrefl)
             ynew = cmir_refl_lower_dp + deltay

             svy = -svy
             reflected = .true.

             mirror = lower

          elseif ((yprev.lt.cmir_refl_lower_dp).and.(ytest.ge.cmir_refl_lower_dp)) then

             deltay = max(abs(cmir_refl_lower_dp-ytest),minrefl)
             ynew = cmir_refl_lower_dp - deltay

             svy = -svy
             reflected = .true.

             mirror = lower

          elseif ((yprev.lt.(-lim_sep+cmir_refl_upper_dp)).and.(ytest.ge.(-lim_sep+cmir_refl_upper_dp))) then 
             !
             !                   Reflection at second upper (< 0) mirror
             !

             deltay = max(abs(-lim_sep+cmir_refl_upper_dp-ytest),minrefl)

             ynew = -lim_sep+cmir_refl_upper_dp - deltay

             svy = -svy
             reflected = .true.

             mirror = upper

          elseif ((yprev.gt.(-lim_sep+cmir_refl_upper_dp)).and.(ytest.le.(-lim_sep+cmir_refl_upper_dp))) then 
             !
             !                   Reflection at second upper (< 0) mirror
             !

             deltay = max(abs(-lim_sep+cmir_refl_upper_dp-ytest),minrefl)

             ynew = -lim_sep+cmir_refl_upper_dp +  deltay

             svy = -svy
             reflected = .true.

             mirror = upper

          endif


       elseif (yprev.gt.0.0) then 

          if (abs(yorg-yprev).gt.lim_sep/2.0) then
             ytest = yorg + lim_sep
          else
             ytest = yorg
          endif

          
          if ((yprev.lt.cmir_refl_upper_dp).and.(ytest.ge.cmir_refl_upper_dp)) then 
             !
             !                   Reflection at upper (>0) mirror
             !
             deltay = max(abs(cmir_refl_upper_dp-ytest),minrefl)

             ynew = cmir_refl_upper_dp - deltay

             svy = -svy
             reflected = .true.

             mirror = upper

          elseif ((yprev.gt.cmir_refl_upper_dp).and.(ytest.le.cmir_refl_upper_dp)) then 
             !
             !                   Reflection at upper (>0) mirror
             !
             deltay = max(abs(cmir_refl_upper_dp-ytest),minrefl)

             ynew = cmir_refl_upper_dp + deltay

             svy = -svy
             reflected = .true.

             mirror = upper

          elseif ((yprev.lt.(lim_sep+cmir_refl_lower_dp)).and.(ytest.ge.(lim_sep+cmir_refl_lower_dp))) then 
             !
             !                   Reflection at second lower (>0) mirror
             !

             deltay = max(abs(lim_sep+cmir_refl_lower_dp-ytest),minrefl)

             ynew = lim_sep + cmir_refl_lower_dp - deltay

             svy = -svy
             reflected = .true.

             mirror = lower


          elseif ((yprev.gt.(lim_sep+cmir_refl_lower_dp)).and.(ytest.le.(lim_sep+cmir_refl_lower_dp))) then 
             !
             !                   Reflection at second lower (>0) mirror
             !
             deltay = max(abs(lim_sep+cmir_refl_lower_dp-ytest),minrefl)

             ynew = lim_sep + cmir_refl_lower_dp + deltay

             svy = -svy
             reflected = .true.

             mirror = lower

          endif
       endif

    endif




    !
    ! Record statistics if reflection occurred
    !

    if (reflected) then

       yreflection_event_count = yreflection_event_count + 1.0
       part_refl_cnt = part_refl_cnt + 1.0 

       if (x.gt.0.0) then 
          ix = min(int(x/dx_inboard),nxbins)
       elseif (x.le.0.0) then 
          ix = max(int(x/dx_outboard)-1,-nxbins)
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

       !
       ! Update tracking of particle frame
       !

       if (yreflection_opt.eq.1) then 
          if (mirror.eq.upper) then 
             part_frame = part_frame + part_mirror_state
          elseif (mirror.eq.lower) then 
             part_frame = part_frame - part_mirror_state
          endif
          part_mirror_state = -part_mirror_state
       endif

    endif


    y = sngl(ynew)

    if (debugl.and.reflected) then
       write(error_message_data,'(a,i10,6(1x,g18.10))') 'REFLECTION OCCURRED:',int(yreflection_event_count),y,ynew,yorg,oldy,svy
       call dbgmsg('CHECK_REFLECTION',error_message_data)
    endif

    ! only perform this check for yreflection_option 1
    if (yreflection_opt.eq.1.and.check_reflected_region(y)) then 

       write(error_message_data,'(a,6(1x,g18.10))') 'REFLECTED PARTICLE HAS ENTERED MIRROR REGION '//&
            & ': DATA:', x,y,yorg,oldy,svy
       call errmsg('CHECK REFLECTION:WARNING:',error_message_data)

       y_tmp = y
       ran = getranf()
       y = cmir_refl_lower + ran * (cmir_refl_upper-cmir_refl_lower)

       write(error_message_data,'(a,i10,l4,10(1x,g20.12))') 'REFLECTION:',int(yreflection_event_count),reflected,x,deltay,ynew,oldy,yorg,yprev,cmir_refl_lower_dp,cmir_refl_upper_dp
       call dbgmsg('CHECK_REFLECTION',error_message_data)

       relocation_count = relocation_count +1.0

       write(error_message_data,'(a,2i12,10(1x,g20.12))') 'REVISED Y COORDINATES:',int(relocation_count),int(yreflection_event_count),y,y_tmp,ynew,yorg,oldy,deltay,svy
       call errmsg('CHECK_REFLECTION PARTICLE RELOCATED:',error_message_data)


       ierr =1
    endif


  end subroutine check_reflection

  subroutine check_x_reflection(x,xorg)
    implicit none
    real :: x,xorg
    ! This routine checks to see if an X reflecting surface has been specified and reflects the particle from that surface

    if (xreflection_opt.eq.0) return

    ! reflect the particle X position in the specified surface if it has an X coordinate greater that the specified surface

    if (x.gt.xreflect_bound.and.xorg.le.xreflect_bound) then
       ! new X is less than the reflecting boundary
       x = xreflect_bound - (abs(x-xreflect_bound))
    elseif (x.lt.xreflect_bound.and.xorg.ge.xreflect_bound) then
       ! new X is greater than the reflecting boundary
       x = xreflect_bound + (abs(x-xreflect_bound))
    endif


  end subroutine check_x_reflection




  subroutine check_x_absorption(x,y,sputy,iz,ierr)
    implicit none
    real :: x,y,sputy
    integer :: iz
    integer :: ierr

    !
    !   Given the X coordinate of the particle check for X absorption
    !
    !
    ierr = 0

    if (x.gt.xabsorb) then 
       ierr =1 
       xabsorb_cnt = xabsorb_cnt + 1.0
       xabsorb_sputy = xabsorb_sputy + sputy
       xabsorb_yavg = xabsorb_yavg + y

       if (iz.eq.0) then 
          xabsorb_neut = xabsorb_neut + sputy
       else
          xabsorb_ion = xabsorb_ion + sputy
          xabsorb_iz = xabsorb_iz + sputy*iz
       endif
    endif


  end subroutine check_x_absorption

  subroutine check_y_absorption(x,y,oldy,sputy,iz,ierr)
    implicit none
    real :: x,y,oldy,sputy
    integer :: ierr,iz

    !   Which ierr (either 1 or 2) means what?
    !
    !   Given the change of the Y coordinate of the particle - check for Y 
    !   absorption across surface - frame dependent
    !
    !
    !
    !   There may be 2 absorber surfaces - the number is controlled by the 
    !   value of yabsorb opt ... data for absorb1 is checked first 
    !

    !                 if (y.gt.4.9.or.oldy.gt.4.9) then 
    !                    write(6,'(a,4(1x,g18.6),3i6)') 'Absorb check:',x,y,oldy,sputy,iz,part_frame,yabsorb1_frame
    !                 endif

    !write(0,*) 'yabsorb_opt, part_frame, yabsorb1_frame = ', yabsorb_opt, part_frame, yabsorb1_frame
    ierr = 0
    if (yabsorb_opt.ge.1.and.part_frame.eq.yabsorb1_frame) then 

       !
       ! Check to see if the particle has crossed either the primary 
       ! absorber in -L < Yabsorb < L or the secondary which differs by +/-2L
       !

       if (in_range(y,yabsorb1a,oldy).or.in_range(y,yabsorb1b,oldy)) then 

          ierr = 1
          yabsorb1_cnt = yabsorb1_cnt+1.0
          yabsorb1_sputy = yabsorb1_sputy + sputy
          yabsorb1_xavg = yabsorb1_xavg + x

          if (iz.eq.0) then 
             yabsorb1_neut = yabsorb1_neut + sputy
          else
             yabsorb1_ion = yabsorb1_ion + sputy
             yabsorb1_iz = yabsorb1_iz + sputy*iz
          endif

       endif

       ! sazmod - Will just try and implement the varying boundary here
       !          as a copy/paste kinda thing bc this is confusing.
       ! Note: This is mostly superseced by the completely varying 2D
       ! bound option. Can safely remove as it adds no extra functionality.

       ! Left step.
       if (vary_absorb.eq.1) then

          ! First need to know if the X position is where the step occurs.
          if (x.lt.xabsorb1a_step) then

             ! Then see if it crossed the step boundary.
             if (in_range(y,yabsorb1a_step,oldy)) then

                ! Do the same things as above.
                ierr = 1
                yabsorb1_cnt = yabsorb1_cnt+1.0
                yabsorb1_sputy = yabsorb1_sputy + sputy
                yabsorb1_xavg = yabsorb1_xavg + x 

                if (iz.eq.0) then 
                   yabsorb1_neut = yabsorb1_neut + sputy
                else
                   yabsorb1_ion = yabsorb1_ion + sputy
                   yabsorb1_iz = yabsorb1_iz + sputy*iz
                endif
                !write(0,*) 'Absorbed: Left face' 
             endif

             ! Also need to make sure it isn't "sneaking" around the step.
             ! So like, the top part of the step facing the radial direction
             ! should also be absorbing. This means we have an absorption
             ! when x < x_step (already satisfied if we're in this 
             ! if-statement) and y < y_step (greater than for the right step).
             if (y.lt.yabsorb1a_step) then

                ierr = 1

                ! I would think that we would need to do some statistics here
                ! to count things, but I'm not sure the framework is in place
                ! to count absorptions on the top of the step, and I don't really
                ! care about it right now anyways. But all that stuff would go here.
                !write(0,*) 'Absorbed: Left top'
             endif
          endif
       endif
    endif

    !write(0,*) 'yabsorb_opt, part_frame, yabsorb2_frame = ', yabsorb_opt, part_frame, yabsorb2_frame
    if (yabsorb_opt.ge.2.and.part_frame.eq.yabsorb2_frame) then 

       !
       ! Check to see if the particle has crossed either the primary absorber in -L < Yabsorb < L or the secondary which differs by +/-2L
       !

       if (in_range(y,yabsorb2a,oldy).or.in_range(y,yabsorb2b,oldy)) then 

          ierr = 1
          yabsorb2_cnt = yabsorb2_cnt+1.0
          yabsorb2_sputy = yabsorb2_sputy + sputy
          yabsorb2_xavg = yabsorb2_xavg + x

          if (iz.eq.0) then 
             yabsorb2_neut = yabsorb2_neut + sputy
          else
             yabsorb2_ion = yabsorb2_ion + sputy
             yabsorb2_iz = yabsorb2_iz + iz*sputy
          endif

       endif

       ! Right step.
       if (vary_absorb.eq.1) then

          ! First need to know if the X position is where the step occurs.
          ! I don't think this will work with two different steps on each side?
          if (x.lt.xabsorb2a_step) then

             ! Then see if it crossed the step boundary.
             if (in_range(y,yabsorb2a_step,oldy)) then

                ! Do the same things as above.
                ierr = 1
                yabsorb2_cnt = yabsorb2_cnt+1.0
                yabsorb2_sputy = yabsorb2_sputy + sputy
                yabsorb2_xavg = yabsorb2_xavg + x 

                if (iz.eq.0) then 
                   yabsorb2_neut = yabsorb2_neut + sputy
                else
                   yabsorb2_ion = yabsorb2_ion + sputy
                   yabsorb2_iz = yabsorb2_iz + sputy*iz
                endif

                !write(0,*) 'Absorbed: Right face'

             endif

             ! Check if it passed top of step. See above comments.
             if (y.gt.yabsorb2a_step) then
                ierr = 1
                !write(0,*) 'Absorbed: Right top'
             endif

          endif
       endif
    endif


  end subroutine check_y_absorption

  subroutine check_y_absorption_2d(ip,ix,x,y,oldy,sputy,iz,ierr)
    implicit none
    real :: x, y, oldy, sputy, yabsorb1a_vary, yabsorb2a_vary
    integer :: ierr, iz, ip, ix
    
    ! sazmod
    ! This subroutine is to check for absorption when using the fully
    ! customizable 2D absorbing boundaries.The main difference with the
    ! subroutine check_y_absorption is just that ip and ix are provided 
    ! so that the correct y value from "bounds_1a" and "bounda_2a" can 
    ! be picked out.
    
    ! Get the y value at this ip, ix.
    yabsorb1a_vary = bounds_1a(ix,ip)
    yabsorb2a_vary = bounds_2a(ix,ip)
    
    ierr = 0
    if (yabsorb_opt.ge.1.and.part_frame.eq.yabsorb1_frame) then 

      ! Check to see if the particle has crossed either the primary 
      ! absorber in -L < Yabsorb < L or the secondary which differs 
      ! by +/-2L.
      if (in_range(y,yabsorb1a_vary,oldy).or.in_range(y,yabsorb1b,oldy)) then 
        ierr = 1
        yabsorb1_cnt = yabsorb1_cnt+1.0
        yabsorb1_sputy = yabsorb1_sputy + sputy
        yabsorb1_xavg = yabsorb1_xavg + x

        if (iz.eq.0) then 
           yabsorb1_neut = yabsorb1_neut + sputy
        else
           yabsorb1_ion = yabsorb1_ion + sputy
           yabsorb1_iz = yabsorb1_iz + sputy*iz
        endif

      endif
       
      ! The varying bound requires an extra step to make sure the
      ! does not "slip through the cracks":
      !
      !  |                              |
      !  |                              |
      !  |                O----         |
      !  |                     ---      |
      !  |                        --    |
      !  |                     |    ---
      !  |                     |       ----
      !  |______________|
      !
      ! There are two ways to solve this. Either a) check for absorption
      ! across a technically unspecified X boundary that spans the gap
      ! shown here or b) if the particle is beyond yabsorb1a_vary, 
      ! count it as absorbed. b) is the simplest. Note this probably
      ! doesn't work with whatever the part_frame = 2 is for.
      if (y.ge.yabsorb1a_vary) then
        ierr = 1
        yabsorb1_cnt = yabsorb1_cnt+1.0
        yabsorb1_sputy = yabsorb1_sputy + sputy
        yabsorb1_xavg = yabsorb1_xavg + x
        if (iz.eq.0) then 
          yabsorb1_neut = yabsorb1_neut + sputy
        else
         yabsorb1_ion = yabsorb1_ion + sputy
         yabsorb1_iz = yabsorb1_iz + sputy*iz
        endif
      endif
    endif
   
    if (yabsorb_opt.ge.2.and.part_frame.eq.yabsorb2_frame) then 

      ! Check to see if the particle has crossed either the primary 
      ! absorber in -L < Yabsorb < L or the secondary which differs 
      ! by +/-2L.
      if (in_range(y,yabsorb2a_vary,oldy).or.in_range(y,yabsorb2b,oldy)) then 

        ierr = 1
        yabsorb2_cnt = yabsorb2_cnt+1.0
        yabsorb2_sputy = yabsorb2_sputy + sputy
        yabsorb2_xavg = yabsorb2_xavg + x

        if (iz.eq.0) then 
           yabsorb2_neut = yabsorb2_neut + sputy
        else
           yabsorb2_ion = yabsorb2_ion + sputy
           yabsorb2_iz = yabsorb2_iz + iz*sputy
        endif
      endif
          
      ! Use le here since we are considering the "left", or negative
      ! side of the simulation volume.
      if (y.le.yabsorb2a_vary) then
        ierr = 1
        yabsorb2_cnt = yabsorb2_cnt+1.0
        yabsorb2_sputy = yabsorb2_sputy + sputy
        yabsorb2_xavg = yabsorb2_xavg + x
        if (iz.eq.0) then 
          yabsorb2_neut = yabsorb2_neut + sputy
        else
          yabsorb2_ion = yabsorb2_ion + sputy
          yabsorb2_iz = yabsorb2_iz + sputy*iz
        endif
      endif

    endif
     
  end subroutine check_y_absorption_2d

  subroutine check_p_reflection(p)
    implicit none
    real p
    ! check to see if p has hit a P reflection boundary ... and then reflect it if it has

    if (preflect_opt.eq.1) then 
       if (abs(p).ge.preflect_bound) then 
          ! limit the value of P to the lesser of the reflected particle position or the boundary
          p = sign(min(abs(2.0 * preflect_bound - abs(p)),preflect_bound),p)
       endif
    endif

  end subroutine check_p_reflection

  subroutine check_p_reflection_neut(dp,dpvelf)
    implicit none
    real*8 :: dp, dpvelf
    real :: p
    ! check to see if p has hit a P reflection boundary ... and then reflect it if it has

    p = sngl(dp)
    if (preflect_opt.eq.1) then 
       if (abs(p).ge.preflect_bound) then 
          ! limit the value of P to the lesser of the reflected particle position or the boundary
          dp = sign(min(abs(2.0 * preflect_bound - abs(p)),preflect_bound),p)
          dpvelf = -dpvelf
       endif
    endif

  end subroutine check_p_reflection_neut


  logical function in_range(y1,y0,y2) 
    implicit none
    real :: y0,y1,y2

    if (y2.ge.y1) then 

       if (y0.ge.y1.and.y0.le.y2) then 
          in_range = .true.
       else
          in_range = .false.
       endif

    else

       if (y0.le.y1.and.y0.ge.y2) then 
          in_range = .true.
       else
          in_range = .false.
       endif
    endif

    !    if (y1.gt.4.9.or.y2.gt.4.9) then
    !       write(6,'(a,3(1x,g18.6),l6)') 'In range:',y1,y0,y2,in_range
    !    endif 

  end function in_range

  subroutine pr_yref_stats
    implicit none

    character*512 :: output_message
    integer in,bm

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

       if (yreflection_opt.eq.1) then 

          call prb
          call prr('  Total number of relocated particles: ',sngl(relocation_count))
          call prb

          call prc('Non-local ion transport and deposition summary: Frames [0,1] are the local Limiter')
          call prc('Table of frame index vs. particle weight at end of particle life')


          write(output_message,'(a7,a7,a20)') 'Frame','LIM','Particles in Frame'
          call prc(output_message)

          do in = -nframes,nframes+1

             if (in.lt.0) then 
                bm = int((in-1)/2)
             else
                bm = int(in/2)
             endif

             write(output_message,'(1x,i6,1x,i6,4x,f18.3)') in,bm,frame_cnt(in)
             call prc(output_message)

          end do

       endif
    endif

    !
    !   Print absorption option statistics
    !
    if (xabsorb_opt.ne.0) then 

       call prb
       call prc('  Summary of X-Absortpion events:')
       call prr('  Location of X absorbing surface    :',xabsorb) 
       call prr('  Total number of X-absorptions      :',sngl(xabsorb_cnt))
       call prr('  Total weight of particles absorbed :',sngl(xabsorb_sputy))
       call prr('    - neutrals                       :',sngl(xabsorb_neut))
       call prr('    - ions                           :',sngl(xabsorb_ion))
       call prr('  Average Y value at absorption      :',sngl(xabsorb_yavg/max(xabsorb_sputy,1.0d0)))
       call prr('  Average ion charge at absorption   :',sngl(xabsorb_iz/max(xabsorb_ion,1.0d0)))
       call prb

    endif

    !
    ! Print Y-absorption statistics  - first absorber
    !

    if (yabsorb_opt.ge.1) then 

       call prb
       call prc('  Summary of Y-Absortpion events on first surface:')
       call prr('  Location of Y absorbing surface ( -L < Y <L )      :',yabsorb1a) 
       call pri('  Frame containing first absorbing surface           :',yabsorb1_frame)
       call prr('  Total number of Y-absorptions on this surface      :',sngl(yabsorb1_cnt))
       call prr('  Total weight of particles absorbed on this surface :',sngl(yabsorb1_sputy))
       call prr('    - neutrals                                       :',sngl(yabsorb1_neut))
       call prr('    - ions                                           :',sngl(yabsorb1_ion))
       call prr('  Average X value at absorption                      :',sngl(yabsorb1_xavg/max(yabsorb1_sputy,1.0d0)))
       call prr('  Average ion charge at absorption                   :',sngl(yabsorb1_iz/max(yabsorb1_ion,1.0d0)))
       call prb

       write(6,'(a,5(1x,g18.6))') 'Yabsorb1:',yabsorb1a,yabsorb1_cnt,yabsorb1_sputy,yabsorb1_neut,yabsorb1_ion,yabsorb1_xavg,yabsorb1_iz

    endif

    ! Print Y-absorption statistics  - second absorber

    if (yabsorb_opt.ge.2) then 

       call prb
       call prc('  Summary of Y-Absortpion events on second surface:')
       call prr('  Location of Y absorbing surface ( -L < Y <L )      :',yabsorb2a) 
       call pri('  Frame containing second absorbing surface          :',yabsorb2_frame)
       call prr('  Total number of Y-absorptions on this surface      :',sngl(yabsorb2_cnt))
       call prr('  Total weight of particles absorbed on this surface :',sngl(yabsorb2_sputy))
       call prr('    - neutrals                                       :',sngl(yabsorb2_neut))
       call prr('    - ions                                           :',sngl(yabsorb2_ion))
       call prr('  Average X value at absorption                      :',sngl(yabsorb2_xavg/max(yabsorb2_sputy,1.0d0)))
       call prr('  Average ion charge at absorption                   :',sngl(yabsorb2_iz/max(yabsorb2_ion,1.0d0)))
       call prb

       write(6,'(a,5(1x,g18.6))') 'Yabsorb1:',yabsorb2a,yabsorb2_cnt,yabsorb2_sputy,yabsorb2_neut,yabsorb2_ion,yabsorb2_xavg,yabsorb2_iz

    endif


  end subroutine pr_yref_stats




end module yreflection
