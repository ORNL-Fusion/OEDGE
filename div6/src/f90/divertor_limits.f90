module divertor_limits

  use mod_cgeom

  private


  integer, parameter :: BND=-1,CORE=1, SOL=2, PFZ=3, ISOL=4, OSOL=5, PFZ2 = 6

  ! Number of Xpoints = 0,1,2 for limiter, SN, DN 
  ! Grid  type is LIMITER, SN, DN and CDN (connected Double null)

  integer,parameter :: LIM=0, SN=1, DN=2, CDN=3

  integer :: grid_type, xpoint_cnt

  !  integer :: sep_ring, isep_ring, osep_ring
  !  integer :: sol_bnd, isol_bnd, osol_bnd
  !  integer :: core_bnd, core_sol
  !  integer :: pfz_bnd, pfz2_bnd. pfz_sol_bnd, pfz2_sol_bnd
  !  integer :: ring_cnt

  real,allocatable :: divertor_limit(:,:) 
  real,allocatable :: divertor_leakage(:,:,:)

  integer :: ring_divertor_leaked
  logical :: divertor_leaked


  public :: calculate_divertor_limit,check_divertor_limit,divertor_leaked,ring_divertor_leaked, print_divertor_leakage


contains

  subroutine init_divertor_limits(nizs)
    use allocate_arrays
    implicit none
    integer :: nizs

    integer :: ierr

    call allocate_array(divertor_limit,nrs,2,'Divertor Limit',ierr)
    if (allocated(divertor_limit)) divertor_limit = 0.0

    call allocate_array(divertor_leakage,nrs,nizs+1,3,'Divertor Leakage',ierr)
    if (allocated(divertor_leakage)) divertor_leakage = 0.0

  end subroutine init_divertor_limits


  subroutine calculate_divertor_limit(nizs)
    use error_handling
    implicit none

    integer :: nizs

    ! local variables
    
    integer :: ik,ir,in,ierr
    integer :: irsep1, irsep2, irpfz1, ikpfz1
    logical :: incore,finished

    integer :: irin_last,irout_last,ikin_last,ikout_last,irin,irout,ikin,ikout,iro,iko


    ! analyse the grid for certain characteristics


    !
    !     jdemod - at this point the grid has been read in and a connection map has been built
    !            - perform additional analysis of the grid
    !            - identify number of Xpoints ... SN, DN or Asymmetric DN
    !            - label each ring on the grid for region ... CORE, MAIN SOL, INNER SOL, OUTER SOL, PRIM PFZ, SEC PFZ
    !            - calculate the divertor boundaries from polygon surfaces attachd to the Xpoint
    !            
    !


    call init_divertor_limits(nizs)

    incore = .true.
    finished = .false.

    xpoint_cnt = 0

    do ir = 1,nrs
       divertor_limit(ir,:) = ksmaxs(ir)/2.0
    end do


!    write(0,'(a,i8)') 'NRS:',nrs

    ! set ir to innermost core ring
    ir =1 

    do while (.not.finished)

       irin_last = irins(1,ir)
       irout_last = irouts(1,ir)

       ikin_last = ikins(1,ir)
       ikout_last = ikouts(1,ir)

       do ik = 1,nks(ir)

          irin = irins(ik,ir)
          ikin = ikins(ik,ir)

          irout = irouts(ik,ir)
          ikout = ikouts(ik,ir)

!          write(0,'(a,10i8)') 'INDICES:', ik,ir,irout,irout_last,irin,irin_last
!          write(6,'(a,10i8)') 'INDICES:', ik,ir,irout,irout_last,irin,irin_last

          if (irout.ne.irout_last) then 
             ! found a ring that adjoins multiple rings working outward. This should only happen for a connected double null.
             if (incore) then 
                grid_type = CDN
                irsep1 = irout
                irsep2 = irout_last
             endif
          endif

          if (irin.ne.irin_last) then 
             ! found point where the ring adjoins two rings inward - working outward this should be the first Xpoint
             if (xpoint_cnt.eq.0) then 
                xpoint_cnt = xpoint_cnt + 1

                irpfz1 = irin_last
                ikpfz1 = ikin_last

                finished = .true.
                incore = .false.

!                write(0,*) 'XPoint_cnt:',xpoint_cnt,irpfz1,ikpfz1
!                write(6,*) 'XPoint_cnt:',xpoint_cnt,irpfz1,ikpfz1

                call assign_limit(ik,ir,kss(ik,ir)/ksmaxs(ir))

             elseif (xpoint_cnt.eq.1) then 

!                write(0,*) 'XPoint_cnt1:',xpoint_cnt,irpfz1,ikpfz1
!                write(6,*) 'XPoint_cnt1:',xpoint_cnt,irpfz1,ikpfz1

                call assign_limit(ik,ir,kss(ik,ir)/ksmaxs(ir))

                if (irin.ne.irpfz1) then 
                   ! symmetric double null
                   xpoint_cnt = xpoint_cnt+1
                   ! navigate to Xpoint on other side through PFZ

                   if (grid_type.ne.CDN) then 
                      call errmsg('ANALYSE GRID','ERROR CDN GRID NOT IDENTIFIED')
                   endif

                   write(6,'(a,10i8)') 'CDN RINGS:',irouts(irpfz1,ikpfz1+1),ikouts(irpfz1,ikpfz1+1),&
                        & irouts(irin,ikin-1),ikouts(irin,ikin-1),irsep1,irsep2,irpfz1,irin,ir

                   iro = irouts(irpfz1,ikpfz1+1)
                   iko = ikouts(irpfz1,ikpfz1+1)

                   call assign_limit(iko,iro,kss(iko,iro)/ksmaxs(iro))

                   iro = irouts(irin,ikin-1)
                   iko = ikouts(irin,ikin-1)

                   call assign_limit(iko+1,iro,kss(iko+1,iro)/ksmaxs(iro))
                else
                   grid_type = SN

                endif

             endif

          endif

!          if (.not.finished) then 
!             divertor_limit(ir,1) = 0.0
!             divertor_limit(ir,2) = ksmaxs(ir)
!          endif

          irin_last = irin
          irout_last = irout

          ikin_last = ikin
          ikout_last = ikout

       end do

       ir = ir + 1

       if (ir.eq.nrs) then
          grid_type = LIM
          finished = .true.
       endif


    end do

    finished = .false. 


    ! if there is a single Xpoint ... need to keep scanning in case there are two. 
    if (xpoint_cnt.eq.1.and.ir.lt.nrs) then 

       do while (.not.finished)

          irin_last = irins(1,ir)
          irout_last = irouts(1,ir)

          ikin_last = ikins(1,ir)
          ikout_last = ikouts(1,ir)

          ! scanning for an outward break 

          if (irout_last.ne.ir) then 

          do ik = 1, nks(ir)

             irin = irins(ik,ir)
             ikin = ikins(ik,ir)

             irout = irouts(ik,ir)
             ikout = ikouts(ik,ir)

             !write(0,'(a,10i8)') 'INDICES2:', ik,ir,irout,irout_last,irin,irin_last
             !write(6,'(a,10i8)') 'INDICES2:', ik,ir,irout,irout_last,irin,irin_last

             if (irout.ne.irout_last) then 
                ! found a ring that adjoins multiple rings working outward. This should happen at a wall inflection point on an extended grid OR at a second Xpoint. Either
                ! way this will be the last ring scanned
                finished = .true. 

                ! The test for Xpoint vs wall tangency is whether the adjacent cells are at the ends of their respective rings or not
                if (ikout.ne.1.and.ikout.ne.nks(irout).and.xpoint_cnt.ne.2) then 
                   ! Xpoint
                   xpoint_cnt = xpoint_cnt + 1
                   grid_type = DN

                   call assign_limit(ikout,irout,kss(ikout,irout)/ksmaxs(irout))
                   call assign_limit(ikout_last+1,irout_last,kss(ikout_last,irout_last)/ksmaxs(irout_last))

                   exit

                endif

             endif

             irin_last = irin
             irout_last = irout

             ikin_last = ikin
             ikout_last = ikout

          end do

          else
             ! exit if hit boundary ring
             finished = .true.
          endif

          ! increment loop count
          ir = ir+1
          if (ir.ge.nrs+1) then
             finished = .true.
          endif

       end do

    endif

    write(6,*) 'Divertor Limits:'

    do ir = 1,nrs

       write(6,'(a,i5,3(1x,g18.8))') 'DIVLIM:',ir,divertor_limit(ir,1),divertor_limit(ir,2),ksmaxs(ir)/2.0

    end do

!    write(0,*) 'Divertor Limits:'
!
!    do ir = 1,nrs
!
!       write(0,'(a,i5,3(1x,g18.8))') 'DIVLIM:',ir,divertor_limit(ir,1),divertor_limit(ir,2),ksmaxs(ir)/2.0
!
!    end do


  end subroutine calculate_divertor_limit


  subroutine assign_limit(ik,ir,frac)
    implicit none
    integer :: ik,ir
    real :: frac

    integer :: irout,ikout, iend
    integer :: irlast,iklast
    logical :: finished

    integer :: ikt,irt

    ikt = ik
    irt = ir

    if (frac.lt.0.5) then 
       iend = 1
    else
       iend = 2
    endif

    finished = .false.

!    write(0,*) 'Assign Limit:',ik,ir,iend,frac
!    write(6,*) 'Assign Limit:',ik,ir,iend,frac

    ! use lower KSB for cell since the test will have stepped one cell past the Xpoint

    do while (.not.finished)

       irlast = irt
       divertor_limit(irt,iend) = ksb(ikt-1,irt)
       irt = irouts(ikt,irt)
       ikt = ikouts(ikt,irt)

!       write(0,*) 'Assign:',irt,ikt,irlast,divertor_limit(irt,iend)
!       write(6,*) 'Assign:',irt,ikt,irlast,divertor_limit(irt,iend)

       if (irlast.eq.irt) then ! reached boundary ring - assign and exit
          divertor_limit(irt,iend) = ksb(ikt-1,irt)
          finished = .true.
       endif

    end do


  end subroutine assign_limit


  logical function check_divertor_limit(sputy,ir,iz,time,s)
    implicit none
    integer :: ir,iz
    real :: s,time,sputy

    check_divertor_limit = ((s.gt.divertor_limit(ir,1)).and.(s.lt.divertor_limit(ir,2)))

    if (check_divertor_limit) then 
       divertor_leakage(ir,iz,1) = divertor_leakage(ir,iz,1) + sputy
       divertor_leakage(ir,iz,2) = divertor_leakage(ir,iz,2) + time*sputy
       divertor_leakage(ir,iz,3) = divertor_leakage(ir,iz,3) + s*sputy
       ring_divertor_leaked = ir
    endif

    return
  end function check_divertor_limit


  subroutine print_divertor_leakage(outunit,nizs)
    implicit none
    integer :: outunit,nizs

    real :: tot_divertor_leakage
    real :: ave_divertor_time
    real :: ave_divertor_s
    integer :: ir,iz

    do ir = 1,nrs
       do iz = 1,nizs
          divertor_leakage(ir,nizs+1,1) = divertor_leakage(ir,nizs+1,1) + divertor_leakage(ir,iz,1)
          divertor_leakage(ir,nizs+1,2) = divertor_leakage(ir,nizs+1,2) + divertor_leakage(ir,iz,2)
          divertor_leakage(ir,nizs+1,3) = divertor_leakage(ir,nizs+1,3) + divertor_leakage(ir,iz,3)

          tot_divertor_leakage = tot_divertor_leakage + divertor_leakage(ir,iz,1)
          ave_divertor_time    = ave_divertor_time    + divertor_leakage(ir,iz,2)
          ave_divertor_s       = ave_divertor_s    + divertor_leakage(ir,iz,3)


          if (divertor_leakage(ir,iz,1).gt.0.0) then 
             divertor_leakage(ir,iz,2) = divertor_leakage(ir,iz,2)/divertor_leakage(ir,iz,1)
             divertor_leakage(ir,iz,3) = divertor_leakage(ir,iz,3)/divertor_leakage(ir,iz,1)
          endif
       end do 
       
       if (divertor_leakage(ir,nizs+1,1).gt.0.0) then 
          divertor_leakage(ir,nizs+1,2) = divertor_leakage(ir,nizs+1,2)/divertor_leakage(ir,nizs+1,1)
          divertor_leakage(ir,nizs+1,3) = divertor_leakage(ir,nizs+1,3)/divertor_leakage(ir,nizs+1,1)
       endif

    end do
    
    if (tot_divertor_leakage.gt.0.0) then 
       ave_divertor_time = ave_divertor_time/tot_divertor_leakage
       ave_divertor_s    = ave_divertor_s   /tot_divertor_leakage
    endif


!
! Print out results summary
!

    call prc('SUMMARY OF DIVERTOR LEAKAGE:')

    do ir = 1,nrs
       if (divertor_leakage(ir,nizs+1,1).gt.0.0) then 

          write(outunit,'(a,i8,2(1x,a,g18.8))') 'RING:',ir,'MINS:',divertor_limit(ir,1),'MAXS:',divertor_limit(ir,2)
          
          do iz = 1,nizs

             if (divertor_leakage(ir,iz,1).gt.0.0) then 
                write(outunit,'(4x,a,i8,3(2x,a,g18.8))') 'CHARGE:',iz,'LEAKAGE:', &
                    & divertor_leakage(ir,iz,1),'AVE_S',divertor_leakage(ir,iz,3),'AVE_T:',divertor_leakage(ir,iz,2)
             endif
          end do

          write(outunit,'(4x,a,i8,3(2x,a,g18.8))') 'ALL IZS',ir,'LEAKAGE:', &
                    & divertor_leakage(ir,iz,1),'AVE_S',divertor_leakage(ir,iz,3),'AVE_T:',divertor_leakage(ir,iz,2)


       endif
    end do

    call prc('TOTAL LEAKAGE:')

    write(outunit,'(3(1x,a,g18.8))') 'LEAKAGE:', tot_divertor_leakage, ' AVE_S :',ave_divertor_s, 'AVE_TIME:', ave_divertor_time


  end subroutine print_divertor_leakage


end module divertor_limits
