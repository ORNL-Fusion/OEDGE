module mod_assign_plasma

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

  real,allocatable:: s(:),ne(:),te(:),ti(:),ef(:),vb(:),dne(:),dte(:),dti(:)


  public :: assign_plasma
  

contains

  subroutine plasma_solver(ixs,ixe,pzs,pze)
    use mod_soledge
    implicit none

    integer :: ixs,ixe,pzs,pze
    
    integer :: ring_type
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

    
    if (plim(pz).eq.1) then ! limiter present
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
          ! no plasma solution calculated for this case - issue error message and exit
          ring_type = 0 
          call errmsg('MOD_ASSIGN_PLASMA: PLASMA ZONE HAS NEITHER LIMITER NOR ABSORBING SURFACES:',pz)
          cycle ! continue to next plasma zone in loop
       endif

    endif
    
    ! Calculate plasma for Y>0 then also map to Y<0
    ! Note: YS are bin boundaries and the first bin boundary is implicitly Y=0
    !       YOUTS are the coordinates of the bin centers that should be used for plasma calculation and assignment 

    
    call calculate_ring(ring_type,youts,nys,ix,pz)




    
    
    
      

    
    do iy = 1,nys
       crnbs(ix,iy,pz) = ne(iy)  ! Y>0
       ctembs(ix,iy,pz) = te(iy)
       ctembsi(ix,iy,pz) = ti(iy)
       velplasma(ix,iy,pz) = vb(iy)
       efield(ix,iy,pz) = ef(iy)
       
       crnbs(ix,iy-nys-1,pz) = ne(iy) ! Y<0
       ctembs(ix,iy-nys-1,pz) = te(iy)
       ctembsi(ix,iy-nys-1,pz) = ti(iy)
       velplasma(ix,iy-nys-1,pz) = vb(iy)
       efield(ix,iy-nys-1,pz) = ef(iy)
    end do   


    

    ! SOLEDGE code -------------------------------
    

    !
    ! Calling SOLEDGE gives a plasma solution for the entire field line from YMIN to YMAX. Solved in two parts
    ! inside the SOLEDGE routine. Data is stored in arrays shared between assign plasma and the plasma calculation
    ! routines
    !

    ! loop over radial planes
    do ix = ixe,ixs
       ! loop over plasma zones
       do pz = pzs,pze
       
          ! pull out the Y-ccordinates


          


          

      smax = ymax-ymin
	!endif
	
    soffset = ymin



    ! set up s axis
    ! Yaxis is s(in) - ymin
    do in = 1, maxn
    !write(0,*) 'checkpoint6 in = ',in,'/',maxn
       sd(in) = (in-1) * smax/(maxn-1)
       !write(0,*) 'SD:',in,sd(in)
    end do
	!write(0,*) 'checkpoint5'
    yd = sd + soffset





    ikstart = 1
    ikmid = maxn/2
    ikend = maxn

    do ix = ix1,ix2
	   !write(0,*) 'ix,x = ', ix,xs(ix)
       !
       !     jdemod - do not read the starting target conditions from the grid
       !            - they should be loaded from the specified target conditions
       !              for each ring
       !
       ! should use qtembs(iqx)?
       !
       
       tebp = ctembs(ix,0,pz)
       tibp = ctembsi(ix,0,pz)
       nbp  = crnbs(ix,0,pz) 
       V0  = - SQRT(0.5*EMI*(TEBP+TIBP)*(1+RIZB)/CRMB)

       tebpi = ctembs(ix,0,pz)
       tibpi = ctembsi(ix,0,pz)
       nbpi  = crnbs(ix,0,pz) 
       V0i  = - SQRT(0.5*EMI*(TEBPi+TIBPi)*(1+RIZB)/CRMB)


       ! jdemod
       ! the special code applying simple options to field lines attached to probes goes away
       ! when proper plasma calculations are made for every field line

       !
       ! Pull out the boundary conditions for applying Vb and ef on field lines that connect
       ! to the probe. Using the simple model for E and Vb from SOL option 4. (linear decay
       ! to stagnation between surfaces
       !

       ! interpolate at youts(iy) to obtain
       ! fix sd axis

       y_1b = yd(1)
       te_1b = te(1)
       ti_1b = ti(1)
       cs_1b =9.79E+03 * SQRT(((te_1b+ti_1b)/2)* (1.0+REAL(CIZB))/CRMB)  

       !write(0,*) ix, iqxs(ix), qedges(iqxs(ix), 1)
       y_1t = -qedges(iqxs(ix),1)
       in = ipos(y_1t,yd,maxn)
       te_1t = te(in)
       ti_1t = ti(in)
       cs_1t =9.79E+03 * SQRT(((te_1t+ti_1t)/2)* (1.0+REAL(CIZB))/CRMB)  

	   !write(0,*) 'y_1b,y_1t,yabsorb1a_step',y_1b,y_1t,yabsorb1a_step

	   cl_1 = (y_1t+y_1b)/2.0

       ! sazmod - The conn length will change if in step region.
       if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step1 = ipos(xabsorb1a_step, xs, nxs-1)
         
         ! Adjust connection length.
         if (ix.le.ix_step1) then
           cl_1 = (y_1t + yabsorb1a_step) / 2.0  ! Y < 0
           !write(0,*) 'cl_1 adjusted to ',cl_1
         endif
       endif
       !write(0,*) 'cl_1 = ',cl_1
       
       ! max efield
       e_1b = te_1b/cl_1
       e_1t = te_1t/cl_1


       y_2b = qedges(iqxs(ix),2)
       in = ipos(y_2b,yd,maxn)
       te_2b = te(in)
       ti_2b = ti(in)
       cs_2b =9.79E+03 * SQRT(((te_2b+ti_2b)/2)* (1.0+REAL(CIZB))/CRMB)  

       y_2t = yd(maxn)
       te_2t = te(maxn)
       ti_2t = ti(maxn)
       cs_2t =9.79E+03 * SQRT(((te_2t+ti_2t)/2)* (1.0+REAL(CIZB))/CRMB) 
       
       cl_2 = (y_2t+y_2b)/2.0 
       
       !write(0,*)'y_2b,y_2t,yabsorb2a_step',y_2b,y_2t,yabsorb2a_step
       ! sazmod - The conn length will vary if in step region.
       if (vary_absorb.eq.1) then
       
         ! Check if in step region.
         ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
         ! Adjust connection length.
         if (ix.le.ix_step2) then           
           cl_2 = (y_2b + yabsorb2a_step) / 2.0  ! Y > 0   
           !write(0,*) 'cl_2 adjusted to ',cl_2       
         endif
       endif
       !write(0,*) 'ix, x, cl_1, cl_2 = ',ix, xs(ix), cl_1, cl_2
       
       ! max efield
       e_2b = te_2b/cl_2
       e_2t = te_2t/cl_2

       TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 

       do iy = -nys,nys


          in = ipos(youts(iy),yd,maxn)

          ! only overwrite if youts is inside of range
          if (in.gt.1.and.in.lt.maxn) then 


             if (.not.(youts(iy).ge.yd(in-1).and.youts(iy).le.yd(in))) then 
                write(0,*) 'Warning Y not in correct bin:', yd(in),youts(iy),yd(in-1)
                stop 'Y index error'
             endif

                
             dy = youts(iy)-yd(in-1)
             dt = yd(in) - yd(in-1)
             ctembs(ix,iy,pz) = te(in-1) + dy/dt * (te(in)-te(in-1))
             ctembsi(ix,iy,pz)= ti(in-1) + dy/dt * (ti(in)-ti(in-1))

             DSTEP = TGSCAL *  QS(IQXS(IX)) * QS(IQXS(IX))
             ctegs(ix,iy,pz) = (teg(in-1) + dy/dt * (teg(in)-teg(in-1))) * dstep
             ctigs(ix,iy,pz) = (tig(in-1) + dy/dt * (tig(in)-tig(in-1))) * dstep

             crnbs(ix,iy,pz)= ne(in-1) + dy/dt * (ne(in)-ne(in-1))

             ! These weird scaing factors are required due to the way the
             ! force coefficients were originally calculated in LIM.
             ! CFEXZS and CFVHXS both have temperature scaling dependencies since the
             ! original LIM code only calculates the velocity and efield along the field
             ! line at the separatrix and scales everything outboard proportional to the
             ! temperature at the limiter. To use the same coefficients with this option
             ! the inverse temperature scaling is applied here so it cancels out.

             if (ctembs(ix,iy,pz).le.0.0.or.ctembsi(ix,iy,pz).le.0.0.or.crnbs(ix,iy,pz).le.0.0) then
                write(0,'(a,2i8,10(g12.5))') 'Less than zero:',ix,iy, ctembs(ix,iy,pz),ctembsi(ix,iy,pz),crnbs(ix,iy,pz)
                stop
             endif
             
             ! jdemod - the scaling factors for LIM original arrays really make no sense for more complex backgrounds
             ! the vel_efield_opt is intended to transition to using efield and velplasma exclusively and includes revised
             ! coefficients that only include timestep scaling and not scaling to temperature at separatrix (CFVHXS, CFVEXS)
             if (vel_efield_opt.eq.0) then 
                e_scale = ctbin/ctembs(ix,iy,pz)
                if ((ctembs(ix,iy,pz).gt.0.0).and.(ctembsi(ix,iy,pz).gt.0.0).and.ctbin.ge.0.0.and.ctibin.ge.0.0) then 
                   v_scale = sqrt((ctbin+ctibin)/(ctembs(ix,iy,pz)+ctembsi(ix,iy,pz)))
                else
                   v_scale = 0.0
                   write(0,'(a,2i8,10(1x,g12.5))') 'V_scale error:',ix,iy,ctbin,ctibin,ctembs(ix,iy,pz),ctembsi(ix,iy,pz)
                endif
             elseif (vel_efield_opt.eq.1) then 
                e_scale = 1.0
                v_scale = 1.0
             endif
                
             ! sazmod
             ! Apply the overall scaling of the plasma velocity. Obviously
             ! nothing will change if vel_mod = 1.0. Note: mod_v_fact
             ! can still do it's own thing, and is applied after vel_mod.
             v_scale = v_scale * vel_mod
               


             ! This needs to be rewritten for poloidal plasma zones
             
             ! sazmod 
             ! Calculate velplasma differently so the step in the absorbing walls
             ! is factored in and has an appropriate stagnation point. 
             if (vary_absorb.eq.1) then  
                             
                ! This is just the equation for a line with the two points
                ! (cl_1, -cs_1b)
                ! (cl_2,  cs_1b)
                !mult = cs_1b * (2*(youts(iy)-cl_2)/(cl_2-cl_1) + 1)
                !write(0,*) 'youts(iy) =', iy, youts(iy)
                !mult = cs_1b * (2*youts(iy)/(cl_1+cl_2) - 1)
                mult = cs_1b * ((youts(iy)-2*cl_1)/(cl_2-cl_1) - 1)
                
                ! Multiply by the v_scale factor because 3DLIM is weird.
                velplasma(ix,iy,pz) = mult * v_scale
                
                ! Printout some numbers.
                !if (mod(iy, 10).eq.0) then
                !  write(0,'(a,2i6,10(1x,g12.5))')'ix,iy,youts,cs1b,mult,vp1=',ix,iy,youts(iy),cs_1b,mult,velplasma(ix,iy,1)
                !endif
                
             else
	             velplasma(ix,iy,pz)= (vb(in-1) + dy/dt * (vb(in)-vb(in-1))) * v_scale
	             efield(ix,iy,pz)= (ef(in-1) + dy/dt * (ef(in)-ef(in-1))) * e_scale
	     endif
             

             ! calculate second set of velplasma and efield for field lines that end on probe surfaces. 

             ! Y < 0
			 !write(0,*) 'youts, qedges(iqxs(ix),2)', youts(iy), qedges(iqxs(ix),2)
             if (youts(iy).lt.ymin) then
                velplasma(ix,iy,2) = 0.0
                efield(ix,iy,2) = 0.0
             elseif (youts(iy).gt.ymax) then 
                velplasma(ix,iy,2) = 0.0
                efield(ix,iy,2) = 0.0
             elseif (youts(iy).lt.-qedges(iqxs(ix),1)) then 

                ! set sign and scaling ... "-" towards -y surface
                !mult =  youts(iy)/cl_1 + 1.0
                
                ! sazmod - I think this should be -youts(iy)/cl_1, otherwise
                ! there isn't a stagnation point on the left side of the probe.
                mult = -youts(iy)/cl_1 + 1.0

                ! Y < 0

				! y < 0 and y < --5
                if (youts(iy).lt.-cl_1) then
				   !write(0,*) 'Section 1'
                   ! half nearest absorbing surface
                   velplasma(ix,iy,2) = cs_1b * mult * v_scale
                   efield(ix,iy,2) = e_1b * mult * e_scale
                   !write(0,*) 'Section 1 ', youts(iy), velplasma(ix,iy,2)
                   
                ! y < 0 and y > --5 (i.e. never!). This should be fixed but
                ! it doesn't seem to cause errors... yet.
                else
                   !write(0,*) 'Section 2'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_1t * mult * v_scale
                   efield(ix,iy,2) = e_1t * mult * e_scale 
                   !write(0,*) 'Section 2 ', youts(iy), velplasma(ix,iy,2) 
                endif


             elseif (youts(iy).gt.qedges(iqxs(ix),2)) then
                ! Y > 0

                ! set sign and scaling ... "-" towards -y surface
                mult =  youts(iy)/cl_2 - 1.0
                !write(0,*) 'Y>0 mult = ', mult

                ! Y < 0

				! y < 5 and y > 0
                if (youts(iy).lt.cl_2) then
                   !write(0,*) 'Section 3'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_2b * mult * v_scale
                   efield(ix,iy,2) = e_2b * mult * e_scale
                   !write(0,*) 'Section 3 ', youts(iy), velplasma(ix,iy,2)
                   
                   ! sazmod
                   if (vary_absorb.eq.1) then
                     
                     ! Check if in right step region.
					 ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
					 ! Mult. velplasma by factor.
					 if (ix.le.ix_step2) then
					   velplasma(ix,iy,2) = velplasma(ix,iy,2) * mod_v_fact
					   !write(0,*) 'mod_v_fact applied'
					 endif                
                   endif
                   
                ! y > 5 and y > 0
                else
                   !write(0,*) 'Section 4'
                   ! on probe side of midpoint
                   velplasma(ix,iy,2) = cs_2t * mult * v_scale
                   efield(ix,iy,2) = e_2t * mult * e_scale 
                   !write(0,*) 'Section 4 ', youts(iy), velplasma(ix,iy,2) 
                   
                   ! sazmod
                   ! My own little ad-hoc factor just to mess with flows
                   ! in the right step region.
                   if (vary_absorb.eq.1) then
                     
                     ! Check if in right step region.
					 ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
         
					 ! Mult. velplasma by factor.
					 if (ix.le.ix_step2) then
					   velplasma(ix,iy,2) = velplasma(ix,iy,2) * mod_v_fact
					   !write(0,*) 'mod_v_fact applied'
					 endif                
                   endif
                   
                endif
             endif
          endif
       end do


       !write(0,'(a,i8)') 'Plasma data for field line:', ix
       !DO IY = -nys,nys
       !   write(0,'(i8,12(1x,g18.8))') iy,youts(iy),ctembs(ix,iy),ctembsi(ix,iy),crnbs(ix,iy),ctegs(ix,iy),ctigs(ix,iy),velplasma(ix,iy,1),efield(ix,iy,1),velplasma(ix,iy,2),efield(ix,iy,2)
       !end do



       

    end do ! ix




    ! SOL22 code

    
       xbnd1 = sol22_regions(ir,1)
       xbnd2 = sol22_regions(ir,2)

       ! ignore pbnds for now
       ! create two zones to overlay SOL22
       ! pzone = 1 -> |P| =< CPC0 (width of probe)
       ! pzone = 2 -> |P| > CPC0
       ! CPC0 
       ! if maxpzones = 1 OR not 3D (MAXNPS =1) then only
       ! calculate for pzone = 1
       ! 
       ! if there are absorbing surfaces - use them
       !    - calculate plasma over [yabsorb1a,yabsorb2a]
       !    - NOTE: there will need to be duplication of
       !      all plasma arrays for 1:pzone areas
       !      e.g. ctembs(ix,iz,pzone)
       !           ctegs ... and forces by iz and pzone
       !
       ! if not ignore
       !    - calculate plasma over whole range and duplicate
       !    - plasma from limiter surface to limiter surface
       !    - SOL22 is NOT applied to pzone 2 since there are
       !      no surfaces there.
       !    - This means this is only applicable for non-3D cases 
       !          
       !
       !pbnd1 = sol22_regions(ir,3)
       !pbnd2 = sol22_regions(ir,4)
       ! SOL22 input specifies which poloidal zone to use the model 
       pz = sol22_regions(ir,3)
       
       do ix = 1,nxs

          ! calculate y locations where the plasma needs to be calculated and
          ! map them to 0 ... ymax (needs to be mapped back afterward)
          ! Also need starting nb, Te, Ti values extracted from existing plasma
          ! or input. Take limits from SOLEDGE code. 

          ! ctembs(ix,iy) ... iy = -nys,nys 
          ! Two sides 
          ! Y < 0 = nys ... nys/2+1 OR -1 ... -nys

          ! only replace specified regions with SOL22 plasma 

          ! take qedges into account for Y coordinates
          ! qedges(iqx,1) for Y< 0 limiter surface
          ! qedges(iqx,2) for Y> 0 limiter surface
          ! HOWEVER - values are +ve for both so edge for Y<0 is -qedges(iqx,1)
          !

          if (xouts(ix).ge.xbnd1.and.xouts(ix).le.xbnd2) then 

             iqx = iqxs(ix)
             ncnt = 0

             !write(0,*) 'SOL22:',ix,iqx,ixout,xouts(ix),xbnd1,xbnd2


             ! only overlay SOL22 outboard when no absorbing surfaces are present
             if (yabsorb_opt.eq.0.and.ix.le.ixout) then 

                ncnt = 0
                do iy = nys,nys/2+1,-1
                   ! values with spts < 0 will be assigned the value of spts = 0
                   if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      spts(ncnt) = 2.0 * cl - qedges(iqx,1) - youts(iy)
                   endif
                end do

                te0 = qtembs(iqx,1)  ! Y < 0 
                ti0 = qtembsi(iqx,1) ! Y < 0 
                nb0 = qrnbs(iqx,1)   ! Y < 0 
                cs  = getcs_sol22(sngl(te0),sngl(ti0))  ! positive flow toward Y<0 surface
                

                ! assign values back to arrays

                ncnt = 0
                do iy = nys,nys/2+1,-1
                   if ((2.0*cl - qedges(iqx,1) - youts(iy)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      crnbs(ix,iy,pz) = nb(ncnt)
                      ctembs(ix,iy,pz) = te(ncnt)
                      ctembsi(ix,iy,pz) = ti(ncnt)
                      velplasma(ix,iy,pz) = -vb(ncnt)
                   else
                      crnbs(ix,iy,pz) = nb0
                      ctembs(ix,iy,pz) = te0
                      ctembsi(ix,iy,pz) = ti0
                      velplasma(ix,iy,pz) = cs
                   endif
                end do


                ! Y > 0 is 1...nys/2 or -nys...-nys/2-1

                ncnt = 0
                do iy = 1,nys/2
                   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
                      ncnt = ncnt + 1
                      spts(ncnt) = youts(iy) - qedges(iqx,2)
                   endif
                end do


                te0 = qtembs(iqx,2)  ! Y > 0 
                ti0 = qtembsi(iqx,2) ! Y > 0 
                nb0 = qrnbs(iqx,2)   ! Y > 0 
                cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward Y>0 surface

                ! calculate SOL




                

                ncnt = 0
                do iy = 1,nys/2
                   if ((youts(iy)-qedges(iqx,2)) .gt. 0.0) then 
                      ncnt = ncnt +1
                      
                      crnbs(ix,iy,pz) = nb(ncnt)
                      ctembs(ix,iy,pz) = te(ncnt)
                      ctembsi(ix,iy,pz) = ti(ncnt)
                      velplasma(ix,iy,pz) = vb(ncnt)
                   else
                      crnbs(ix,iy,pz) = nb0
                      ctembs(ix,iy,pz) = te0
                      ctembsi(ix,iy,pz) = ti0
                      velplasma(ix,iy,pz) = cs
                   endif

                end do


                ! copy over modified ctembs etc to -nys:-1

                do iy = 1,nys
                   crnbs(ix,-nys+iy-1,pz) = crnbs(ix,iy,pz)
                   ctembs(ix,-nys+iy-1,pz) = ctembs(ix,iy,pz)
                   ctembsi(ix,-nys+iy-1,pz) =  ctembsi(ix,iy,pz)
                   velplasma(ix,-nys+iy-1,pz) =  velplasma(ix,iy,pz)
                end do

                ! copy velplasma from pz = 1 to pz = 2
                !do in = 2,maxpzone
                !   velplasma(:,:,in) = velplasma(:,:,1)
                !end do
                
                


             elseif (yabsorb_opt.ge.2) then

                ! add code for dealing with absorption options - with or without limiter being present
                




             endif


             pz = 1
             ! print out results
             do iy = -nys,nys

                write(6,'(a,2i8,20(1x,g12.5))') 'SOL22:',ix,iy,iqxs(ix),qrnbs(iqx,1),qtembs(iqx,1),qtembs(iqx,1),&
                     qrnbs(iqx,2),qtembs(iqx,2),qtembs(iqx,2),youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz)
             end do


          endif

       end do


    end do


    

    

  end subroutine plasma_solver
  

    subroutine calculate_ring(ring_type,youts,nys,ix,pz)
      implicit none
      integer :: ring_type,ix,pz
      real :: ys
      
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

         call set_plasma_data_axis(youts,1,int(nys/2),nys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
         call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)
         
         call calculate_plasma
         call calculate_tgrad_e(ix,pz)
         
         ! Assign plasma to ring section

         do iy = 1,nys/2
            call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz),efield(ix,iy,pz),solver_axis_opt)
         end do
         
         ! Y<0 absorbing
         call get_boundary_conditions(3,n1,te1,ti1,ix,pz,bnd1)
         ! Y<0 limiter
         call get_boundary_conditions(1,n2,te2,ti2,ix,pz,bnd2)
         ! Note: by convention qedges is always >0 while the absorbing surface for Y<0 is less than zero
         ! Also - we want to perform this calculation for the CL to 2CL portion of the plasma
         ! so everything is mapped to the second half of the Y>0 plasma
         ! Shift the bounds to CL to 2CL
         bnd1 = 2.0 * CL + bnd1  ! Note value of bnd1 should be <0 for the Y<0 absorbing surface
         bnd2 = 2.0 * CL - bnd2  ! sign of bnd2 is +ve by code convention for qedges

         call set_plasma_data_axis(youts,int(nys/2+1),nys,nys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
         call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

         call calculate_plasma
         call calculate_tgrad_e(ix,pz)
         
         ! Assign plasma to ring section - the current plasma conditions calculated and the bounds are stored in the plasma module

         do iy = nys/2+1,nys
            call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz),efield(ix,iy,pz),solver_axis_opt)
         end do
         
         ! copy plasma to Y<0

         do iy = 1,nys
            crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
            ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
            ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
            velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
            efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
         end do
                  
      elseif (ring_type.eq.2) then
         ! qedge2 -> 2cl-qedge1 ... map qedge2 -> 0 (first surface)
         ! This is the easier one to map
         
         ! Y>0 limiter
         call get_boundary_conditions(2,n1,te1,ti1,ix,pz,bnd1)

         ! Y<0 limiter mapped to 2CL-EDGE
         call get_boundary_conditions(1,n2,te2,ti2,ix,pz,bnd2)

         ! map limiter edge to 2CL-bnd2
         bnd2 = 2.0 *CL - bnd2

         call set_plasma_data_axis(youts,1,nys,nys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
         call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)

         call calculate_plasma
         call calculate_tgrad_e(ix,pz)
         
         ! Assign plasma to ring section - can do the entire ring in one ordered loop

         do iy = 1,nys
            call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz),efield(ix,iy,pz),solver_axis_opt)
         end do

         ! copy plasma to Y<0
         
         do iy = 1,nys
            crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
            ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
            ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
            velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
            efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
         end do

      elseif (ring_type.eq.3) then 
         ! -abs1 -> abs2  ... map -abs1-> 0 (first surface)
         ! This is the hardest to map to only Y>0 - easier to map this one from -CL/2 to CL/2 and then copy the missing halves
         ! Y<0 absorber - value of bnd1 has -ve for absorbing surface
         call get_boundary_conditions(3,n1,te1,ti1,ix,pz,bnd1)

         ! Y>0 absorber
         call get_boundary_conditions(4,n2,te2,ti2,ix,pz,bnd2)

         call set_plasma_data_axis(youts,-nys/2,nys/2,nys,bnd1,bnd2,solver_axis_opt) ! alternate approach - may be needed if interpolation doesn't work out
         call set_boundary_conditions(n1,te1,ti1,n2,te2,ti2)
         
         ! This calculates the plasma from -abs2 to +abs1
         call calculate_plasma
         call calculate_tgrad_e(ix,pz)
         
         ! Assign plasma to ring section - do the central section

         do iy = -nys/2,nys/2
            call assign_plasma(youts(iy),crnbs(ix,iy,pz),ctembs(ix,iy,pz),ctembsi(ix,iy,pz),velplasma(ix,iy,pz),efield(ix,iy,pz),solver_axis_opt)
         end do

         ! copy plasma to missing sections
         
         ! copy to Y < -CL

         do iy = 1,nys/2
            crnbs(ix,iy-nys-1,pz) = crnbs(ix,iy,pz) 
            ctembs(ix,iy-nys-1,pz) = ctembs(ix,iy,pz) 
            ctembsi(ix,iy-nys-1,pz) = ctembsi(ix,iy,pz) 
            velplasma(ix,iy-nys-1,pz) = velplasma(ix,iy,pz) 
            efield(ix,iy-nys-1,pz) = efield(ix,iy,pz) 
         end do

         ! copy to Y > CL
         
         do iy = -nys/2,-1
            crnbs(ix,iy+nys+1,pz) = crnbs(ix,iy,pz) 
            ctembs(ix,iy+nys+1,pz) = ctembs(ix,iy,pz) 
            ctembsi(ix,iy+nys+1,pz) = ctembsi(ix,iy,pz) 
            velplasma(ix,iy+nys+1,pz) = velplasma(ix,iy,pz) 
            efield(ix,iy+nys+1,pz) = efield(ix,iy,pz) 
         end do

      endif

    end subroutine calculate_ring



    subroutine get_boundary_conditions(surf,n,te,ti,ix,pz,bnd)
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
      iqx = iqxs(ix)
      
      if (surf.lt.1.or.surf.gt.4) then
         call errmsg('MOD_ASSIGN_PLASMA:GET_BOUNDARY_CONDITIONS','ERROR Invalid surface specified')
         return
      endif
      
      if ((surf.eq.1.or.surf.eq.2).or.((surf.eq.3.or.surf.eq.4).and.nabsorb_plasma.eq.0)) then 
         ! no absorption plasma data - absorbtion surface data matches limiter data
         if (surf.eq.1.or.surf.eq.3) then 
            !tebp = ctembs(ix,0,pz)
            !tibp = ctembsi(ix,0,pz)
            !nbp  = crnbs(ix,0,pz) 

            te0 = qtembs(iqx,1)  ! Y < 0 
            ti0 = qtembsi(iqx,1) ! Y < 0 
            nb0 = qrnbs(iqx,1)   ! Y < 0 

            if (surf.eq.1) then 
               bnd = qedges(iqx,1)  ! Y < 0
            elseif (surf.eq.3) then
               bnd = yabsorb_surf(ix,pz,1)
            endif
         elseif (surf.eq.2.or.surf.eq.4) then

            te0 = qtembs(iqx,2)  ! Y > 0 
            ti0 = qtembsi(iqx,2) ! Y > 0 
            nb0 = qrnbs(iqx,2)   ! Y > 0 

            if (surf.eq.2) then 
               bnd = qedges(iqx,2)  ! Y > 0
            elseif (surf.eq.4) then
               bnd = yabsorb_surf(ix,pz,2)
            endif
            !cs  = -getcs_sol22(sngl(te0),sngl(ti0))  ! negative flow toward Y>0 surface

         endif

      elseif (nabsorb_plasma.gt.0.and.(surf.eq.3.or.surf.eq.4)) then

         call get_absorb_bounds(surf,xs(ix),n,te,ti,ix,pz,bnd)

      endif

      return


    end subroutine get_boundary_conditions

    subroutine get_absorb_bounds(surf,x,n,te,ti,ix,pz,bnd)
      implicit none
      integer :: surf,ix,pz
      real :: x
      real :: n,te,ti
      real   :: fact
      
      ! locals
      integer :: in,iqx

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

            fact = (x-absorb_plasma(in-1,1))/((absorb_plasma(in,1)-absorb_plasma(in-1,1))
            n  = absorb_plasma(in-1,2)+ fact * (absorb_plasma(in,2)-absorb_plasma(in-1,2)
            te = absorb_plasma(in-1,3)+ fact * (absorb_plasma(in,3)-absorb_plasma(in-1,3)
            ti = absorb_plasma(in-1,4)+ fact * (absorb_plasma(in,4)-absorb_plasma(in-1,4)            

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

            fact = (x-absorb_plasma(in-1,5))/((absorb_plasma(in,5)-absorb_plasma(in-1,5))
            n  = absorb_plasma(in-1,6)+ fact * (absorb_plasma(in,6)-absorb_plasma(in-1,6)
            te = absorb_plasma(in-1,7)+ fact * (absorb_plasma(in,7)-absorb_plasma(in-1,7)
            ti = absorb_plasma(in-1,8)+ fact * (absorb_plasma(in,8)-absorb_plasma(in-1,8)

         endif
            
      endif


      return
    end subroutine get_absorb_bounds

      

    



  
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
