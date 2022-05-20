module mod_assign_plasma

  use mod_assign_plasma_input
  implicit none

  ! The purpose of this module is to generalize an interface between SOL22 and possibly SOLEDGE
  ! so that these routines can calculate ne, te, ti, ef, vb, dte, dti, dne for a high resolution
  ! grid and then interpolate those values onto the actual simulation mesh being used by the LIM run.
  !
  ! Interpolation is more flexible but likely less accurate.
  ! The alternative is to load the specific set of coordinates desired into an input array and then call this
  ! routine with those coordinates ensuring that these coordinate are among any grid points requested from the
  ! solver routines.
  !

  integer :: maxn = 1000

  !real,allocatable:: s(:),ne(:),te(:),ti(:),ef(:),vb(:),dne(:),dte(:),dti(:)

  integer :: debug_step = 0
  

  public


contains


  subroutine setup_solvers(qtim,crmb,cizb,yscale,ixout,cprint)
    use mod_plasma_data
    implicit none
    real :: qtim, crmb,yscale
    integer :: cizb,ixout,cprint
    ! this is required to bring in some global variables in an easy fashion instead of adding
    ! it to every subroutine in the call stack

    ! It is also used for any other solver setup required
    qtim_local = qtim
    crmb_local = crmb
    cizb_local = cizb
    yscale_local = yscale
    ixout_local = ixout
    cprint_local = cprint   ! set up for future use if there is a need to pass cprint through to the solver code
    
  end subroutine setup_solvers
  
  
  subroutine plasma_solver(ixs,ixe,pzs,pze,solver_opt)
    use mod_params
    use mod_comxyt
    use yreflection
    use mod_plasma_data
    use mod_comtor
    implicit none
    ! solver_opt = 0 - soledge (2PM)  ... = 1 - sol22

    integer :: ixs,ixe,pzs,pze,solver_opt,iqx
    integer :: ring_type
    integer :: ix,pz
    !
    ! This code is called from the plasma_overlay routine in tau.f. The overlay routine handles other sorts of overlays in
    ! addition to the SOLEDGE and SOL22 options.
    !
    !
    ! Enumerating cases:
    !
    ! 1) No absorbing surfaces - within 3D limiter extent - standard limiter style configuration.
    !    Plasma is calculated from QEDGES(-ixqs,2) to 2L-QEDGES(-ixqs,1) over the 
    !    span of Y=0 to 2L. The Y>0 side of the limiter is in QEDGES(-iqxs,2). 
    !
    !    SOLEDGE or SOL22 called - generate plasma solutions for entire array- two
    !    solutions one for each half ring.
    !
    !    Assign plasma to ctembs(ix,iy,pz) for 0 < Y < 2L then copy the values
    !    0 < Y < 2L -> -2L < Y < 0 so that value at -2L+ = values at 0+
    !  
    !    Note: If using 3D probe extents without absorbing surfaces then use base plasma
    !          options to define plasma conditions on flux tubes with no material interactions
    !          present.
    !
    !    Code copies the base plasma option solution onto all flux surfaces before the
    !    assign_plasma pieces are overlaid. Examples would be to assign a constant plasma outside the
    !    probe/limiter area for example. 
    !
    !  2) Absorbing surfaces are present - this code only supports 2 particle absorbing surfaces that
    !     are representative of targets - one in -L < Y < 0 and one in 0 < Y < L. The absorbers do not need
    !     to be at the same Y value. 
    !
    !     A) For 2D or scenarios with a limiter/probe on the field line between the absorbing surfaces
    !     - Two separate calls to the plasma overlays are required. 
    !     - call for -YabsorbA to -QEDGES(iqxs,1)
    !     - call for QEDGES(iqxs,2) to +YabsorbB
    !     - these two solutions full solutions are overlaid onto the ctembs(ix,iy,pz) array first for
    !       iy = -nys/2 to -1 and then for iy = 1 to nys/2 (The solver could alternatively be called for
    !     - 2L-YabsorbA to 2L-QEDGES(iqxs,1) then overlaid on iy = nys/2 to nys
    !     - the entire solution needs to be replicated onto the portions of the remaining plasma background
    !       where it was not directly calculated to get a plasma solution for -2L to +2L
    !
    !     B) scenarios without a probe/limiter surface between the absorbing surfaces.  
    !        - one call to the plasma solver for -YabsorbA -> +YabsorbB
    !        - this solution is then overlaid onto iy = -nys/2 to nys/2
    !        - then the -nys/2,0 segment is copied to nys/2+1 to nys and 
    !               the 0,nys/2 segment is copied to -nys,-nys/2-1
    !
    !     These scenarios cover the most frequent use cases. This code does not support a single
    !     absorbing surface on a field line at the moment. This scenario might be applicable to
    !     a situation with a probe located between limiter surfaces but the same effect can be
    !     obtained by appropriate use of the yabsorb_surf array specifying surfaces at the required
    !     distances upstream and downstream of the probe. 
    !
    !
    !     There are three calculation/plasma assignment cases that are required.
    !
    !     1) Calculation and assignment to an entire field line from limiter/probe -> limiter/probe
    !     2) Calculation and assignment from limiter/probe -> absorbing surface - 2/field line
    !     3) Calculation and assignment from downstream absorbing surface to upstream absorbing surface
    !
    !     The plasma solutions are determined by calls to either SOLEDGE or SOL22
    !
    !     Target conditions on the probe/limiter are determined from standard input values.

    !     Target conditions at the absorbing surfaces are either
    !     - the same as the values at the limiter/probe
    !     - separately specified input
    !
    !     Each poloidal plasma zone for SOLEDGE uses the same parameters but may have different target conditions. 
    !
    !     Each poloidal plasma zone for SOL22 can have different parameters specified through the use
    !     of an external file (name included in the SOL22 input). 
    !     This is part of the SOL22 additional input for LIM.
    !
    ! Note:    
    !
    ! Solve and apply plasma to central region and then copy to Y<-CL and Y>CL regions. 
    !
    ! a) no limiter present - no absorbing surfaces present
    !    - issue error message and exit since use of SOLEDGE and SOL22 assumes the existence of
    !      surfaces on the flux tube
    ! b) no limiter present - two absorbing surfaces - one Y<0 and one Y>0
    !    - solve from Y<0 absorbing surface to Y>0 absorbing surface. 
    !    - possible different boundary conditions at each surface
    !    - apply plasma to the background arrays starting from -nys/2 to 0 then 0 to nys - offset Y to match and interpolate
    ! c) limiter present - no absorbing surfaces
    !    - calculate plasma in 2 sections Y=[0,-CL] Y=[0,CL] surfaces are at Y=0 - or one section with different boundary conditions on each side
    !    - then map onto plasma arrays - might be easiest to do this one for all Y>0 then map back to Y<0
    ! d) limiter present - two absorbing surfaces
    !    - calculate plasma in two sections
    !    - could be either all Y>0 or split Y<0 and Y>0.
    !    - Y<0 limiter surface to Y<0 absorbing surface, Y>0 limiter surface to Y>0 absorbing surface with possibly different boundary conditions on each
    !
    ! Optional input needed for absorbing surface plasma conditions - these are functions that can be interpolated 
    ! nabsorb_plasma
    ! absorb_plasma
    !   PZ, X_abs1, n_abs1, Te_abs1, Ti_abs1, X_abs2, n_abs2, Te_abs2, Ti_abs2   1-> Y<0 surface, 2->Y>0 surface
    !   Support for different plasma zones included but not implemented yet
    !

    ! Determine which type of ring


    do ix = ixs,ixe

       do pz = pzs,pze

          if (plimz(pz).eq.1.and.xouts(ix).le.0.0) then ! limiter present on zone and X is outboard
             ! check for absorbing surfaces
             if (yabsorb_opt.ne.0) then  ! y absorbing surfaces present
                ring_type = 1 ! ring has both limiter surfaces and absorbing surfaces
             else
                ring_type = 2 ! ring has only limiter surfaces
             endif

          else ! no limiter present
             ! check for absorbing surfaces
             if (yabsorb_opt.ne.0) then ! y absorbing surfaces present
                ring_type = 3 ! ring has only absorbing surfaces
             else
                ! ring has neither limiter surfaces nor absorbing surfaces -> NO boundary conditions - neither SOLEDGE or SOL22 is appropriate
                ! no plasma solution calculated for this case - issue warning message and loop
                ring_type = 0 
                call errmsg('MOD_ASSIGN_PLASMA: PLASMA ZONE HAS NEITHER LIMITER NOR ABSORBING SURFACES'//&
                      &' - DEFAULT PLASMA USED FOR ZONE - SOLVERS CAN NOT BE USED FOR X>0  ZONE=',pz)
                cycle ! continue to next plasma zone in loop
             endif

          endif

          ! Calculate plasma for Y>0 then also map to Y<0
          ! Note: YS are bin boundaries and the first bin boundary is implicitly Y=0
          ! YOUTS are the coordinates of the bin centers that should be used for plasma calculation and assignment
          !
          ! Scaling here are for the temperature gradient and electric field arrays
          !
          ! Electric field is later multipled by the cfexzs factor which includes the entire scaling factor
          ! EMI/CRMI * QTIM_local *QTIM_local * QS(IQXS(IX)) * QS(IQXS(IX)) plus the Z dependence and a
          ! temperature relative to the separatrix that needs to be cancelled out. Final scaling of the
          ! electric field efield and the velocity velplasma is handled in this routine so that the default values
          ! for velplasma (cvhys) and efield (ceys) remain correct - these do not need to correct for the temperature scaling
          ! However, any efield or velplasma values calculated using other methods need to remove the temperature dependence. 
          !
          ! CFSS contains QTIM*QS(IQX) so the velocity scaling needs only one factor of this
          ! unless it gets put in cfvhxs in which case cvhys does not need scaling. 
          ! Put all of the scaling in CFEXZS - as for default options then scale
          ! the velocity
          
          iqx = iqxs(ix)

          ! jdemod
          ! note - having the time step multiplier qs depend on iqx rather than ix is a logical inconsistency
          ! since the impurity density arrays are recorded in ix as well as the background plasma conditions so the
          ! other arrays lack the spatial resolution of iqx making changing the timestep based on an iqx division of
          ! the modeling space appear to make little sense
          tg_scale = EMI/CRMI * QTIM_local *QTIM_local * QS(IQX)  * QS(IQX)

          ! ef_scale is set to 1.0 - scaling is applied through the array cfexzs and with inverse temperature factors in
          ! calculate_ring
          ef_scale = 1.0
          
          call calculate_ring(ring_type,youts,nys,ix,xouts(ix),iqx,pz,solver_opt)

          write(0,*) 'Solver Status: Finished plasma for zone = ',pz,' and radial bin = ',ix

       end do
    end do

  end subroutine plasma_solver


  subroutine calculate_ring(ring_type,youts,nys,ix,x,iqx,pz,solver_opt)
    use mod_params
    use mod_plasma_data
    use mod_comt2
    use mod_comtor
    implicit none
    integer :: ix,iqx,iy,pz,solver_opt,ring_type,nys
    real :: x,youts(-maxnys:maxnys)

    real :: n1,te1,ti1,n2,te2,ti2
    real :: bnd1,bnd2
    real :: scale
    
    if (ring_type.eq.1) then
       ! for Y>0
       ! 2 halves
       ! qedge1 -> abs1
       ! use 2cl-abs1 -> 2cl-qedge1 ... map 2cl-abs1 -> 0 (first surface)
       !
       ! Loop through youts to assign results
       !

       ! qedge2 -> abs2 ... map qedge2 -> 0 (first surface)
       ! Y>0 limiter
       call get_boundary_conditions(2,n1,te1,ti1,ix,pz,bnd1)
       ! Y>0 absorbing surface
       call get_boundary_conditions(4,n2,te2,ti2,ix,pz,bnd2)

       call set_plasma_data_axis(youts,1,int(nys/2),maxnys,bnd1,bnd2,solver_axis_opt) 
       call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

       !write(6,'(a,3i8,10(1x,g12.5))') 'Solver1:',ring_type,ix,pz,n1,te1,ti1,n2,te2,ti2,bnd1,bnd2,bnd2-bnd1

       call calculate_plasma(solver_opt,pz,ix,x)

       ! Assign plasma to ring section

       do iy = 1,nys/2
          call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),&
                             velplasma(ix,iy,pz),efield(ix,iy,pz),ctegs(ix,iy,pz),ctigs(ix,iy,pz))
          ! apply temperature scalings to cancel those implicitly in cfvhxs and cfexzs
          ! velocity scale
          if (ctembs(ix,iy,pz).gt.0.0.and.ctembsi(ix,iy,pz).gt.0.0) then 
             scale =sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz)))
             velplasma(ix,iy,pz) = velplasma(ix,iy,pz) * scale
          else
             velplasma(ix,iy,pz) = 0.0
          endif

          ! efield scale
          if (ctembs(ix,iy,pz).gt.0.0) then 
             if (ix.gt.ixout_local) then  ! inboard
                ! efield also has a field length scaling outboard in cfexzs that needs to be removed              
                scale = ctbin/ctembs(ix,iy,pz) 
             else ! outboard
                scale = ctbin/ctembs(ix,iy,pz) * yscale_local/cyscls(iqx)  
             endif
             efield(ix,iy,pz) = efield(ix,iy,pz) * scale
          else
             efield(ix,iy,pz) = 0.0
          endif             
       end do

       ! Y<0 absorbing
       call get_boundary_conditions(3,n1,te1,ti1,ix,pz,bnd1)
       ! Y<0 limiter
       call get_boundary_conditions(1,n2,te2,ti2,ix,pz,bnd2)
       ! Note: by convention qedges is always >0 while the absorbing surface for Y<0 is less than zero
       ! Also - we want to perform this calculation for the CL to 2CL portion of the plasma
       ! so everything is mapped to the second half of the Y>0 plasma
       ! Shift the bounds to CL to 2CL ... note CTWOL is 2 * CL = limiter separation
       bnd1 = CTWOL + bnd1  ! Note value of bnd1 should be <0 for the Y<0 absorbing surface
       bnd2 = CTWOL - bnd2  ! sign of bnd2 is +ve by code convention for qedges

       call set_plasma_data_axis(youts,int(nys/2+1),nys,maxnys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out

       !write(6,'(a,3i8,10(1x,g12.5))') 'Solver2:',ix,pz,ring_type,n1,te1,ti1,n2,te2,ti2,bnd1,bnd2,bnd2-bnd1

       call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

       call calculate_plasma(solver_opt,pz,ix,x)

       ! Assign plasma to ring section - the current plasma conditions calculated and the bounds are stored in the plasma module

       do iy = nys/2+1,nys
          call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),&
                             velplasma(ix,iy,pz),efield(ix,iy,pz),ctegs(ix,iy,pz),ctigs(ix,iy,pz))
          ! apply temperature scalings to cancel those implicitly in cfvhxs and cfexzs
          if (ctembs(ix,iy,pz).gt.0.0.and.ctembsi(ix,iy,pz).gt.0.0) then 
             scale =sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz)))
             velplasma(ix,iy,pz) = velplasma(ix,iy,pz) * scale
          else
             velplasma(ix,iy,pz) = 0.0
          endif
          if (ctembs(ix,iy,pz).gt.0.0) then 
             if (ix.gt.ixout_local) then 
             ! efield also has a field length scaling outboard in cfexzs that needs to be removed              
                scale = ctbin/ctembs(ix,iy,pz) 
             else
                scale = ctbin/ctembs(ix,iy,pz) * yscale_local/cyscls(iqx)  
             endif
             efield(ix,iy,pz) = efield(ix,iy,pz) * scale
          else
             efield(ix,iy,pz) = 0.0
          endif             
       end do

       ! copy plasma to Y<0

       do iy = 1,nys
          crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
          ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
          ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
          velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
          efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
          ctegs(ix,iy-nys-1,pz) = ctegs(ix,iy,pz) 
          ctigs(ix,iy-nys-1,pz) = ctigs(ix,iy,pz) 
       end do

    elseif (ring_type.eq.2) then
       ! qedge2 -> 2cl-qedge1 ... map qedge2 -> 0 (first surface)
       ! This is the easier one to map

       ! Y>0 limiter
       call get_boundary_conditions(2,n1,te1,ti1,ix,pz,bnd1)

       ! Y<0 limiter mapped to 2CL-EDGE
       call get_boundary_conditions(1,n2,te2,ti2,ix,pz,bnd2)

       ! map limiter edge to 2CL-bnd2
       bnd2 = CTWOL - bnd2

       call set_plasma_data_axis(youts,1,nys,maxnys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
       call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

       !write(6,'(a,3i8,10(1x,g12.5))') 'Solver3:',ix,pz,ring_type,n1,te1,ti1,n2,te2,ti2,bnd1,bnd2,bnd2-bnd1

       call calculate_plasma(solver_opt,pz,ix,x)

       ! Assign plasma to ring section - can do the entire ring in one ordered loop

       do iy = 1,nys
          call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),&
                             velplasma(ix,iy,pz),efield(ix,iy,pz),ctegs(ix,iy,pz),ctigs(ix,iy,pz))


          ! apply temperature scalings to cancel those implicitly in cfvhxs and cfexzs
          if (ctembs(ix,iy,pz).gt.0.0.and.ctembsi(ix,iy,pz).gt.0.0) then 
             scale =sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz)))
             velplasma(ix,iy,pz) = velplasma(ix,iy,pz) * scale
          else
             velplasma(ix,iy,pz) = 0.0
          endif
          if (ctembs(ix,iy,pz).gt.0.0) then 
             if (ix.gt.ixout_local) then 
             ! efield also has a field length scaling outboard in cfexzs that needs to be removed              
                scale = ctbin/ctembs(ix,iy,pz)
             else
                scale = ctbin/ctembs(ix,iy,pz) * yscale_local/cyscls(iqx)  
             endif
             efield(ix,iy,pz) = efield(ix,iy,pz) * scale
          else
             efield(ix,iy,pz) = 0.0
          endif             
       end do

       ! copy plasma to Y<0

       do iy = 1,nys
          crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
          ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
          ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
          velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
          efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
          ctegs(ix,iy-nys-1,pz) = ctegs(ix,iy,pz) 
          ctigs(ix,iy-nys-1,pz) = ctigs(ix,iy,pz) 
       end do

    elseif (ring_type.eq.3) then 
       ! -abs1 -> abs2  ... map -abs1-> 0 (first surface)
       ! This is the hardest to map to only Y>0 - easier to map this one from -CL/2 to CL/2 and then copy the missing halves
       ! Y<0 absorber - value of bnd1 has -ve for absorbing surface
       call get_boundary_conditions(3,n1,te1,ti1,ix,pz,bnd1)

       ! Y>0 absorber
       call get_boundary_conditions(4,n2,te2,ti2,ix,pz,bnd2)

       call set_plasma_data_axis(youts,-nys/2,nys/2,maxnys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
       call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

       !write(6,'(a,3i8,10(1x,g12.5))') 'Solver4:',ix,pz,ring_type,n1,te1,ti1,n2,te2,ti2,bnd1,bnd2,bnd2-bnd1
       
       ! This calculates the plasma from -abs2 to +abs1
       call calculate_plasma(solver_opt,pz,ix,x)

       ! Assign plasma to ring section - do the central section

       do iy = -nys/2,nys/2
          call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),&
                             velplasma(ix,iy,pz),efield(ix,iy,pz),ctegs(ix,iy,pz),ctigs(ix,iy,pz))


          ! apply temperature scalings to cancel those implicitly in cfvhxs and cfexzs
          if (ctembs(ix,iy,pz).gt.0.0.and.ctembsi(ix,iy,pz).gt.0.0) then 
             scale =sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz)))
             velplasma(ix,iy,pz) = velplasma(ix,iy,pz) * scale
          else
             velplasma(ix,iy,pz) = 0.0
          endif

          if (ctembs(ix,iy,pz).gt.0.0) then 
             if (ix.gt.ixout_local) then 
             ! efield also has a field length scaling outboard in cfexzs that needs to be removed              
                scale = ctbin/ctembs(ix,iy,pz)
             else
                scale = ctbin/ctembs(ix,iy,pz) * yscale_local/cyscls(iqx)  
             endif
             efield(ix,iy,pz) = efield(ix,iy,pz) * scale
          else
             efield(ix,iy,pz) = 0.0
          endif             
       end do

       ! copy plasma to missing sections

       ! copy to Y < -CL

       do iy = 1,nys/2
          crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
          ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
          ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
          velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
          efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
          ctegs(ix,iy-nys-1,pz) = ctegs(ix,iy,pz) 
          ctigs(ix,iy-nys-1,pz) = ctigs(ix,iy,pz) 
       end do

       ! copy to Y > CL

       do iy = -nys/2,-1
          crnbs(ix,iy+nys+1,pz) = crnbs(ix,iy,pz) 
          ctembs(ix,iy+nys+1,pz) = ctembs(ix,iy,pz) 
          ctembsi(ix,iy+nys+1,pz) = ctembsi(ix,iy,pz) 
          velplasma(ix,iy+nys+1,pz) = velplasma(ix,iy,pz) 
          efield(ix,iy+nys+1,pz) = efield(ix,iy,pz) 
          ctegs(ix,iy+nys+1,pz) = ctegs(ix,iy,pz) 
          ctigs(ix,iy+nys+1,pz) = ctigs(ix,iy,pz) 
       end do

    endif


    !write out for debugging:
    !if (ix.ge.2.and.ix.le.5.and.pz.eq.1) then 

    !do iy=-nys,nys
    !   write(6,'(a,4i8,20(1x,g12.5))') 'Plasma:',ring_type,ix,iy,pz,x,youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz),velplasma(ix,iy,pz)/sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz))) ,efield(ix,iy,pz),ctegs(ix,iy,pz),ctigs(ix,iy,pz)
    !end do
    !endif
    
    ! plasma data allocation is done in the set_plasma_data_axis routine where the number of points
    ! on the ring is algorithmically determined. Each call will deallocate and reallocate the arrays to
    ! the modified size. The last set needs to be explicitly deallocated. 
    call deallocate_plasma_data

    return

  end subroutine calculate_ring



  subroutine get_boundary_conditions(surf,n,te,ti,ix,pz,bnd)
    use mod_params
    use mod_comt2
    use mod_comxyt
    use yreflection
    implicit none

    integer :: surf,ix,pz
    real :: n,te,ti,bnd

    ! this routine returns the plasma conditions at the specified surface to be used as boundary conditions for the solvers
    ! surf - defines which surface the data is for
    !        1 = limiter Y<0 surface  
    !        2 = limiter Y>0 surface
    !        3 = Y<0 absorbing surface
    !        4 = Y>0 absorbing surface
    !
    !        Data are taken first from the already calculated plasma conditions
    !        If data is specified separately for absorbing surfaces then this data is 
    !        is interpolated and returned - otherwise absorbing surface conditions will
    !        be the same as 'limiter surface' conditions for the matching side of the limiter
    !
    real :: x
    integer :: iqx


    iqx = iqxs(ix)
    x = xouts(ix)


    ! this routine should not be called for inboard limiter surfaces but could be called for
    ! inboard absorbing surfaces. 

    if (surf.lt.1.or.surf.gt.4) then
       call errmsg('MOD_ASSIGN_PLASMA:GET_BOUNDARY_CONDITIONS','ERROR Invalid surface specified')
       return
    endif

    if ((surf.eq.1.or.surf.eq.2).or.((surf.eq.3.or.surf.eq.4).and.nabsorb_plasma.eq.0)) then 

       ! Note - a bin center should never be equal to zero. This must be a bin boundary since
       ! it defines the radial region between field lines that hit the probe/limiter and those that do not
       ! If it ever happens - use the second set and post an error message

       if (x.eq.0.0) then
          call errmsg('MOD_ASSIGN_PLASMA:GET_BOUNDARY_CONDITIONS: XOUTS(IX)=0.0 FOR IX=',ix)
       endif
       
       if (x.lt.0.0) then  ! X < 0

          ! no absorption plasma data - absorbtion surface data matches limiter data
          if (surf.eq.1.or.surf.eq.3) then   ! Y < 0
             !tebp = ctembs(ix,0,pz)
             !tibp = ctembsi(ix,0,pz)
             !nbp  = crnbs(ix,0,pz) 

             te = qtembs(iqx,1)  ! Y < 0 
             ti = qtembsi(iqx,1) ! Y < 0 
             n  = qrnbs(iqx,1)   ! Y < 0 

             if (surf.eq.1) then 
                bnd = qedges(iqx,1)  ! Y < 0
             elseif (surf.eq.3) then
                bnd = yabsorb_surf(ix,pz,1)
             endif
          elseif (surf.eq.2.or.surf.eq.4) then  ! Y > 0

             te = qtembs(iqx,2)  ! Y > 0 
             ti = qtembsi(iqx,2) ! Y > 0 
             n  = qrnbs(iqx,2)   ! Y > 0 

             if (surf.eq.2) then 
                bnd = qedges(iqx,2)  ! Y > 0
             elseif (surf.eq.4) then
                bnd = yabsorb_surf(ix,pz,2)
             endif
             !cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward Y>0 surface

          endif

       else  ! X > 0
          
          ! no absorption plasma data - absorbtion surface data matches limiter data
          if (surf.eq.1.or.surf.eq.3) then   ! Y < 0
             !tebp = ctembs(ix,0,pz)
             !tibp = ctembsi(ix,0,pz)
             !nbp  = crnbs(ix,0,pz) 

             te = ctembs(ix,-1,pz)  ! Y < 0 
             ti = ctembsi(ix,-1,pz) ! Y < 0 
             n  = crnbs(ix,-1,pz)   ! Y < 0 

             if (surf.eq.1) then 
                bnd = 0.0   ! Y < 0
             elseif (surf.eq.3) then
                bnd = yabsorb_surf(ix,pz,1)
             endif
          elseif (surf.eq.2.or.surf.eq.4) then  ! Y > 0

             te = ctembs(ix,1,pz)  ! Y > 0 
             ti = ctembsi(ix,1,pz) ! Y > 0 
             n  = crnbs(ix,1,pz)   ! Y > 0 

             if (surf.eq.2) then 
                bnd = 0.0  ! Y > 0
             elseif (surf.eq.4) then
                bnd = yabsorb_surf(ix,pz,2)
             endif
             !cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward Y>0 surface

          endif

       endif

    elseif (nabsorb_plasma.gt.0.and.(surf.eq.3.or.surf.eq.4)) then

       call get_absorb_bounds(surf,x,n,te,ti,ix,pz,bnd)

    endif

    write(6,'(a,5i8,20(1x,g12.5))') 'Boundary conditions:',surf,nabsorb_plasma,pz,iqx,ix,xouts(ix),n,te,ti,bnd,&
         &crnbs(ix,1,pz),crnbs(ix,-1,pz),ctembs(ix,1,pz),ctembs(ix,-1,pz),ctembsi(ix,1,pz),ctembsi(ix,-1,pz)

    
    return
  end subroutine get_boundary_conditions

  subroutine get_absorb_bounds(surf,x,n,te,ti,ix,pz,bnd)
    use yreflection
    use allocatable_input_data
    implicit none
    integer :: surf,ix,pz
    real :: x
    real :: n,te,ti,bnd
    real   :: fact

    ! locals
    integer :: in,iqx
    integer,external :: ipos

    ! interpolate into absorption boundary data which has been loaded into
    ! separate arrays to facilitate interpolation

    ! remember the X axis is positive inboard and negative outboard with most negative value
    ! at the wall
    !
    ! To use IPOS the array of X coordinates must be in ascending order

    if (surf.eq.3) then ! Y<0
       bnd = yabsorb_surf(ix,pz,1)  ! 
       if (x.le.absorb_plasma(1,1)) then  ! assign edge values
          n  = absorb_plasma(1,2)
          te = absorb_plasma(1,3)
          ti = absorb_plasma(1,4)
       elseif (x.ge.absorb_plasma(nabsorb_plasma,1)) then
          n  = absorb_plasma(nabsorb_plasma,2)
          te = absorb_plasma(nabsorb_plasma,3)
          ti = absorb_plasma(nabsorb_plasma,4)
       else ! interpolate

          in = ipos(x,absorb_plasma(:,1),nabsorb_plasma)

          fact = (x-absorb_plasma(in-1,1))/(absorb_plasma(in,1)-absorb_plasma(in-1,1))
          n  = absorb_plasma(in-1,2)+ fact * (absorb_plasma(in,2)-absorb_plasma(in-1,2))
          te = absorb_plasma(in-1,3)+ fact * (absorb_plasma(in,3)-absorb_plasma(in-1,3))
          ti = absorb_plasma(in-1,4)+ fact * (absorb_plasma(in,4)-absorb_plasma(in-1,4))           

       endif
       if (n.eq.0.0.or.te.eq.0.0.or.ti.eq.0.0) then
          write(0,'(a,5i8,20(1x,g12.5))') 'Error in get_absorb_bounds Y<0:',surf,pz,in,nabsorb_plasma,ix,x,n,te,ti,absorb_plasma(1,1),absorb_plasma(nabsorb_plasma,1),fact,&
               &absorb_plasma(in-1,2),absorb_plasma(in,2),(absorb_plasma(in,2)-absorb_plasma(in-1,2)),bnd
       endif

    elseif (surf.eq.4) then ! Y>0
       bnd = yabsorb_surf(ix,pz,2)
       if (x.le.absorb_plasma(1,5)) then ! assign edge values
          n  = absorb_plasma(1,6)
          te = absorb_plasma(1,7)
          ti = absorb_plasma(1,8)
       elseif (x.ge.absorb_plasma(nabsorb_plasma,5)) then
          n  = absorb_plasma(nabsorb_plasma,6)
          te = absorb_plasma(nabsorb_plasma,7)
          ti = absorb_plasma(nabsorb_plasma,8)
       else ! interpolate

          in = ipos(x,absorb_plasma(:,5),nabsorb_plasma)

          fact = (x-absorb_plasma(in-1,5))/(absorb_plasma(in,5)-absorb_plasma(in-1,5))
          n  = absorb_plasma(in-1,6)+ fact * (absorb_plasma(in,6)-absorb_plasma(in-1,6))
          te = absorb_plasma(in-1,7)+ fact * (absorb_plasma(in,7)-absorb_plasma(in-1,7))
          ti = absorb_plasma(in-1,8)+ fact * (absorb_plasma(in,8)-absorb_plasma(in-1,8))

       endif
       if (n.eq.0.0.or.te.eq.0.0.or.ti.eq.0.0) then
          write(0,'(a,5i8,20(1x,g12.5))') 'Error in get_absorb_bounds Y>0:',surf,pz,in,nabsorb_plasma,ix,x,n,te,ti,absorb_plasma(1,5),absorb_plasma(nabsorb_plasma,5),fact,&
               &absorb_plasma(in-1,6),absorb_plasma(in,6),(absorb_plasma(in,6)-absorb_plasma(in-1,6)),bnd
       endif

    endif

    

    return
  end subroutine get_absorb_bounds



  subroutine calculate_plasma(solver_opt,pz,ix,x)
    use mod_soledge
    use mod_solcommon
    use mod_sol22_lim
    use mod_plasma_data
    implicit none
    integer :: solver_opt
    integer :: pz,ix
    real :: x

    if (solver_opt.eq.0) then 

       ! plasma is calculated from lower absorbing surface to
       ! upper absorbing surface - this allows for
       ! asymmetric placement of the probe
       !call init_soledge(yabsorb1a,yabsorb2a)
       !call init_soledge()
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
       call soledge
       !call soledge(1,nxs/2,qtim)
       !endif


    elseif (solver_opt.eq.1) then 
       !
       !     This code calculates the plasma conditions for sections of the simulation
       !     volume using SOL22. There are several scenarios.
       !
       !         
       ! Initialize some output options in SOL22 using values from slcom
       call init_solcommon(0,0)
       call sol22

    endif


    if (debug_step.gt.0) call prt_plasma(pz,ix,x,debug_step)

    
    return


  end subroutine calculate_plasma




!  subroutine allocate_locals(n)
!    use mod_params
!    use allocate_arrays
!    implicit none
!      integer :: n,ierr
!      
!      call allocate_array(s,1,n,'Local s',ierr)
!      call allocate_array(ne,1,n,'Local ne',ierr)
!      call allocate_array(te,1,n,'Local te',ierr)
!      call allocate_array(ti,1,n,'Local ti',ierr)
!      call allocate_array(ef,1,n,'Local ef',ierr)
!      call allocate_array(vb,1,n,'Local vb',ierr)
!      call allocate_array(dne,1,n,'Local dne',ierr)
!      call allocate_array(dte,1,n,'Local dte',ierr)
!      call allocate_array(dti,1,n,'Local dti',ierr)
!
!    end subroutine allocate_locals
!
!    subroutine deallocate_locals
!      implicit none
!      
!      if (allocated(s)) deallocate(s)
!      if (allocated(ne)) deallocate(ne)
!      if (allocated(te)) deallocate(te)
!      if (allocated(ti)) deallocate(ti)
!      if (allocated(ef)) deallocate(ef)
!      if (allocated(vb)) deallocate(vb)
!      if (allocated(dne)) deallocate(dne)
!      if (allocated(dte)) deallocate(dte)
!      if (allocated(dti)) deallocate(dti)
!
!    end subroutine deallocate_locals


end module mod_assign_plasma
