c     -*Former Mode Specification*-
c     
      subroutine do_parallel_step(seed,nrand,neutim,
     >                            spara,dspara,vpara,dvpara)
      implicit none
c     
      real*8 seed 
      real neutim
      real spara,dspara,vpara,dvpara
      integer nrand  
c     
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'clocal'
      include    'reiser_com' 
c     
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div5'
      include 'div6'
c     
      include    'particle_specs'
      include    'driftvel'
c     
c     Force functions 
c     
      real force_ff,force_fig,force_feg,force_fvg,force_fe
      external force_ff,force_fig,force_feg,force_fvg,force_fe
      real     delta_s_dperpz,ds_kpinchs
      external delta_s_dperpz,ds_kpinchs

c     
c     Local variables
c     
      real ds_dperpz,ds_pinch
      real bg_drftvel,imp_drftvel 
c     
c     Set ion and/or background drift velocities 
c     
      call set_drift_velocity(s,ir,bg_drftvel,imp_drftvel)
c     
c     Set up dvpara or dspara depending on options - for collisional diffusive transport
c     
      call set_collisional_step(seed,nrand,spara,vpara,
     >                          lfps(ik,ir,iz),lllfps(ik,ir,iz))
c     
c     
c     Set quantites for use in force calculations 
c     
c     Ion velocity
c     
      fvel = vel
c     
c     Plasma flow velocity
c     
      FVH   = KVHS(IK,IR) + bg_drftvel

c     
c     Get forces
c     
c     psmod
c     
c     Select Reiser or normal force calculations
c     
      IF(CIOPTR.EQ.1.or.cioptr.eq.2)THEN



c     
c     Call decision to examine gradients in the cell and determine
c     if it is appropriate to use the Reiser formulation for the 
c     forces.
c     
         CALL DECISION(COPTION,IK,IR)
         IF(COPTION.NE.0)THEN
c     
c     If option 2 - update coeffiecients each time step 
c     
            if (cioptr.eq.2) then  
               call update_reiser_coeff(ik,ir,iz,s,lambda2,fvh)
            else
               LAMBDA2 = LAMBDA1(IZ)*KNBS(IK,IR)
            endif 
c     
c     Calculate Reiser forces  
c     
            CHIpara = ALPHAI(IK,IR)*(FVEL-FVH)
c     
            CALL COULOMB_COLL(XKpara,XDparapara,CHIpara,LAMBDA2,
     >           IK,IR,K11,K12,K13,D11,D12,D13)

            FF  = K11 * QTIM2
            FIG = K12 * QTIM2
            FVG = K13 * QTIM2 
c     
            kk = kk + 1 

c     
c     Reiser code ion transport - gaussian value
c     

 7703       nrand = nrand + 1
            call surand2(seed,1,ran1)
            if (ran1.eq.0.0) goto 7703
            nrand = nrand + 1
            call surand2(seed,1,ran2)
            rgauss = sqrt(-2.0*log(ran1))*cos(2.0*PI*ran2)


c     
c     Reset Vpara for Reiser fomulation
c     
            vpara  = SQRT(XDparapara*QTIM) * RGAUSS *QTIM
            spara  = 0.0
c     
         else
c     
c     Default force behaviour - as standard
c     
            FF    = force_ff(ik,ir,iz,fvh,fvel)
            FIG   = force_fig(ik,ir,iz,s,smax)
            FVG   = force_fvg(ik,ir)

         endif

      else
c     
         FF    = force_ff(ik,ir,iz,fvh,fvel)
         FIG   = force_fig(ik,ir,iz,s,smax)
         FVG   = force_fvg(ik,ir)
      endif
c     
c     psmod
c     


      FEG   = force_feg(ik,ir,iz,s,smax)
      FE    = force_fe(ik,ir,iz,s)
c     
      DS_DPERPZ = delta_s_dperpz(ik,ir,nrand)
      ds_pinch  = ds_kpinchs(ik,ir)
      
      call force_col(dvpara,dspara,vpara,spara,kk)       

c     
c     Calculate change in velocity due to forces
c     
      QUANT = FF + FE + FEG + FIG + FVG
c     
c     
c     changed order of s and vel calculation... Krieger IPP 12/94
c     
c     
c     Record last S value in order to detect parallel reflection situations  
c     
      slast = s
c     
      S     = S + VEL + 0.5 * (QUANT+dvpara) + dspara
     >     + IMP_DRFTVEL + ds_dperpz + ds_pinch
c     
      VEL   = VEL + QUANT + dvpara
c     

      if (debug_all) then 
         write(6,'(a,i4,5(1x,g12.5))') 
     >       'DO  PARA 1:',ifate,s,theta,cross,vel,quant
         write(6,'(a,i4,10(1x,g12.5))') 
     >       'DO  PARA 1:',ifate,ff,fe,fig,feg,dvpara,ds_dperpz,ds_pinch
      endif

c     
c     Accumulate some statistics on forces
c     
      call save_force_data (dvpara)
c     
c     Check to see if the parallel reflection option is active
c     and if parallel ion reflection is possible on this ring.  
c     
      if (s_reflect_opt.eq.1.and.
     >     (s_reflect(ir,1).ne.0.0.or.s_reflect(ir,2).ne.0.0)) then 
c     
c     Check reflection conditions. 
c     
c     Reflection from downstream side 
c     
         if (slast.lt.s_reflect(ir,1).and.
     >        s.ge.s_reflect(ir,1)) then 
c     
            s = s_reflect(ir,1) - (s-s_reflect(ir,1))
            vel = -vel
c     
         elseif (slast.gt.s_reflect(ir,2).and.
     >           s.le.s_reflect(ir,2)) then 
c     
            s = s_reflect(ir,2) + (s_reflect(ir,2)-s)
            vel = -vel
c     
         endif
c     
      endif 
c     
c     For collision option 7 - reverse sign of VEL under
c     certain circumstances.
c     
      if (cioptb.eq.7.and.dspara.ne.0.0) then
         if (vel.gt.0.0.and.dspara.lt.0.0) then
            vel = -vel
         elseif (vel.lt.0.0.and.dspara.gt.0.0) then
            vel = -vel
         endif
      endif

      if (debug_all) write(6,'(a,i4,5(1x,g12.5))') 
     >     'DO  PARA 2:',ifate,s,slast,theta,cross,smax


      if (ir.lt.irsep) then 


C     
C--------LOOPING ROUND MAIN PLASMA CONTOURS
C     
         IF (S.LT.0.0) THEN
 600        S = S + SMAX
            IF (S.LT.0.0) GOTO 600
         ELSEIF (S.GT.SMAX) THEN
 610        S = S - SMAX
            IF (S.GT.SMAX) GOTO 610
         ENDIF


c     
c     Particle not in core 
c     

      else
c     
c     Routine may set IFATE for target impact
c     
         call check_target_impact(seed,nrand,neutim)

      endif



      if (debug_all) write(6,'(a,i4,5(1x,g12.5))') 
     >     'DO  PARA E:',ifate,s,slast,theta,cross




c     
c     Update cell and deperp information based on new S-position
c     

      call update_parallel
c
c     
c     if (ifate.ne.0) return 
c     
c-------------------------------------------------------------------c
c     
      return
      end
c     
c     
c     
      real function force_fe(ik,ir,iz,s)
c     
      implicit none
      integer ik,ir,iz
      real s 
c     
      include 'params'
      include 'cgeom'
      include 'hc_global_opts'
c     
      real local_efield
c     
      if (global_hc_follow_option.ne.0) then 
c     
         call hc_electric_field_mod(ik,ir,iz,s,local_efield)
c     
      else
c     
         local_efield =  KES(IK,IR)
c     
      endif
c     
      force_FE    = real(iz) * local_efield
c     
      return
      end
c     
c     
c     
      real function force_ff(ik,ir,iz,fvh,fvel)
      implicit none
      integer ik,ir,iz
      real fvh,fvel
c     
      include 'params'
      include 'clocal'      
      include 'cioniz'
c     
      force_FF    = KFSSMOD(IK,IR)  * LFSS(IK,IR,IZ) * (FVH-FVEL)
c
      return
      end
c     
c     
c     
      real function force_fig(ik,ir,iz,s,smax)
      implicit none
c     
      integer ik,ir,iz
      real s,smax
c     
      include 'params'
      include 'comtor'
      include 'cgeom'
c     
      if (cioptn.eq.3.and.s.gt.cstgrad*smax
     >     .and.s.lt.smax*(1.0-cstgrad)) then
         force_fig = 0.0

      else
         force_FIG   = KBETAS(IZ) * KFIGS(IK,IR)
      endif 

c     
      return
      end
c     
c     
c     
      real function force_feg(ik,ir,iz,s,smax)
      implicit none
c     
      integer ik,ir,iz
      real s,smax
c     
      include 'params'
      include 'comtor' 
      include 'cgeom'
c     
c     Calculate modifications to forces if any
c     
      if (cioptm.eq.3.and.s.gt.cstgrad*smax
     >     .and.s.lt.smax*(1.0-cstgrad)) then
         force_feg = 0.0
c     
      else 

         force_FEG   = KALPHS(IZ) * KFEGS(IK,IR)

      endif

      return
      end
c     
c     
c     
      real function force_fvg(ik,ir)
      implicit none
      integer ik,ir

      force_FVG   = 0.0

      return
      end
c     
c     
c     
      subroutine force_col(dvpara,dspara,vpara,spara,kk)
      implicit none
c     
      real dvpara,dspara,vpara,spara
      integer kk
c     
      include 'params'
      include 'crand' 

c     
c     Only need one random number since VPARA and SPARA based methods of 
c     diffusion based collisonal transport are mutually exclusive.
c     
      KK   =  KK + 1
      dspara= SIGN(SPARA,RANV(KK)-0.5)
      dvpara= sign(vpara,ranv(kk)-0.5)

      return
      end
c
c
c
      real function ds_kpinchs(ik,ir)
      implicit none
c
      integer ik,ir
c
      include 'params'
      include 'cgeom'
      include 'comtor'

      !
      ! Radial flow options 8 and 9 can result in effective 
      ! parallel displacements by moving particles into adjacent
      ! flux tubes (moving in the P direction) 
      !
      if (pinchopt.eq.8.or.pinchopt.eq.9.or.pinchopt.eq.10) then 

         ds_kpinchs = kpinchs_para(ik,ir)
         
      else
         ds_kpinchs = 0.0

      endif


      return
      end
c     
c     
c     
      real function delta_s_dperpz(ik,ir,nrand)
      implicit none
      integer ik,ir,nrand
      include 'params'
      include 'cgeom'
      include 'dperpz'
c     
c     This routine returns a deltaS displacement that would result
c     from a cross-field step occurring in the Z or P (paramagnetic direction). 
c     
c     
      real fact,ran1
      real*8 seed
c     
      delta_s_dperpz = 0.0
c     
      if (dperpz_opt.ne.0) then
c     
c     Need to calculate Btor/Bpol from KBFS which is Btot/Bpol
c     Btot = sqrt(Bpol**2 + Btor**2)
c     KBFS**2 = (Bpol**2 + Btor**2) / Bpol**2 = 1 + (Btor/Bpol)**2
c     Btor/Bpol = sqrt(kbfs**2 -1)
c     
         fact = sqrt(kbfs(ik,ir)**2-1.0)
c     
         nrand=nrand+1
         call surand2(seed,1,ran1)
c     
         delta_s_dperpz = sign(base_dperpz_step * fact,ran1-0.5)
c     
      endif
c     
      return
      end  
c     
c     
c     
      subroutine init_dperpz
      implicit none
      include 'params'
      include 'comtor'
      include 'dperpz'
c     
c     Initialize the DperpZ Delta S transport option
c     - if this option is active, additional deltaS 
c     steps for the ions will be calculated - these 
c     deltaS steps are actually the result of dperp steps
c     in the orthogonal direction not present in DIVIMP
c     directly. 
c     
c     
c     This routine calculates the basic Dperp step in the additional
c     dimension. This is used in the calls to delta_s_dperpz. Note that
c     this routine does not currently support the spatially varying options
c     for DPERP. 
c     
c     
      if (dperpz_opt.eq.0) then 
         base_dperpz_step = 0.0
      elseif (dperpz_opt.eq.1) then 
         base_dperpz_step = sqrt(2.0 * cdperp * qtim)
      endif
c     
      return
      end
c     
c     
c     
      subroutine update_parallel
      implicit none
      include    'params'
      include    'comtor'
      include    'cgeom'
c     
      include 'div1'
      include 'div6'
c     
      include    'particle_specs'

c     
C     
C--------FIND NEAREST IK CORRESPONDING TO DISTANCE S ALONG CONTOUR IR
C--------FIND THETA VALUE CORRESPONDING TO DISTANCE S ALONG CONTOUR IR
C--------ADJUST CROSS FIELD TERM FOR NEW DISTANCES BETWEEN CONTOURS
C     
      IKOLD = IK
      IROLD = IR

c     
c     Find new cell
c     
      IF (PDOPT.EQ.1) THEN
 620     IF (IK.LT.NKS(IR).AND.S.GT.KSB(IK,IR)) THEN
            IK = IK + 1
            GOTO 620
         ENDIF
 625     IF (IK.GT.1.AND.S.LT.KSB(IK-1,IR)) THEN
            IK = IK - 1
            GOTO 625
         ENDIF
      ELSEIF (PDOPT.eq.0) then 
 630     IF (IK.LT.NKS(IR).AND.S.GT.KSS(IK,IR)) THEN
            IK = IK + 1
            GOTO 630
         ENDIF
 635     IF (IK.GT.1.AND.S.LE.KSS(IK-1,IR)) THEN
            IK = IK - 1
            GOTO 635
         ENDIF
         IF (IK.GT.1.AND.S-KSS(IK-1,IR).LT.KSS(IK,IR)-S) IK = IK - 1   
      ENDIF
c     
c     If using non-orthogonal transport - find theta value
c     
c     
c slmod begin
      if (.not.(ir.lt.irsep.and.cgridopt.eq.LINEAR_GRID))
     .  call calculate_theta(ik,ir,s,theta) 
c
c      call calculate_theta(ik,ir,s,theta) 
c slmod end
c     
c     
c     Adjust cross-field term - if necessary
c     
      call adjust_cross(cross,adjust,ik,ir,ikold,irold,debug_all) 
c     
c     Record some statistics in the core. 
c     
      if (ir.lt.irsep.and.adjust.ne.0.0) then 
         IF (ADJUST.LE.0.0) THEN
            DCROSS(1) = DCROSS(1) + 1.0
            DCROSS(2) = DCROSS(2) + ADJUST
         ELSE
            DCROSS(3) = DCROSS(3) + 1.0
            DCROSS(4) = DCROSS(4) + ADJUST
         ENDIF
      endif
C     
      IF (DEBUG_ALL) WRITE (6,1001) 'D1:',IK,IR,S,K,
     >     THETA,SMAX,CROSS,adjust,
     >     'UPDATED S'


c     
c     Format statements
c     
 1000 format(a,2i4,1p,10(g11.4),1x,a) 
 1001 format(a,2i4,1p,6(g11.4),44x,1x,a) 



      if (debug_all) write(6,'(a,5(1x,g12.5))') 
     >     'UPD PARA  :',s,slast,smax,theta,cross




      return 
      end
c     
c     
c     
      subroutine set_collisional_step(seed,nrand,spara,vpara,
     >                                lfps,lllfps)
      implicit none
c
c     Note: Passing variables as arguments make the routine more accessible from
c           other sections of the code - however, care must be taken with the options
c           selected to make sure that all required values are set properly for the 
c           specified options. For example TEMI is required to be properly set
c           for some options - this is the case from the DIVIMP ion transport routines
c           but may not be the case if this code is called from either the HC or periphery
c           ion transport code.  
c     
      real*8 seed
      integer nrand
      real spara,vpara
      real lfps,lllfps
c     
      include    'params'
      include    'comtor'
c      include    'clocal'
      include    'crand'
      include 'div1'
      include 'div2'
      include 'div5'
      include 'div6'
c     
      include    'particle_specs'

c     
c     calculate SPARA and VPARA
c     
      spara = 0.0
      vpara = 0.0
C     
C------CALCULATE PARALLEL DIFFUSION FACTOR SPARA
C     
      IF (DIFFUS) THEN
         SPARA = TEMI * LLLFPS
      ELSE
         IF (CDIFOP.EQ.2) THEN
            IF (LFPS.GT.0.0) THEN
               RCONST = 2.0 * TEMI / LFPS
            ELSE
               RCONST = HI
            ENDIF
         ENDIF
         IF (CIST.GE.RCONST .OR. RCONST.LT.1.0) THEN
            RDIFFT = RDIFFT + CIST * QTIM * SPUTY
            DIFFUS = .TRUE.
            SPARA  = TEMI * LLLFPS
         ELSE
            SPARA  = 0.0
         ENDIF
      ENDIF
C     
C------COLLISION OPTIONS 3,4,7,9 - RESET SPARA IF NECESSARY
C     
      IF (IR.GE.IRSPEC) THEN
         IF (CIOPTB.EQ.3.OR.CIOPTB.EQ.4.or.cioptb.eq.7
     >        .or.cioptb.eq.9) THEN
            KK = KK + 1
            IF (RANV(KK).GT.LFPS/(2.0*TEMI)) SPARA = 0.0
         ENDIF
      ENDIF
c     
c------Collision option 5 - set parallel diffusive velocity step
c     
      if (cioptb.eq.5.or.cioptb.eq.10) then
         vpara = 0.0
         spara  = 0.0
         KK = KK + 1
         IF (RANV(KK).LE.LFPS/(2.0*TEMI))
     >        vpara = lllfps
c     
      elseif (cioptb.eq.6.or.cioptb.eq.11) then
c     
c     Collision option 6 ... changes on each time-step
c     
         spara  = 0.0
         vpara = lllfps
c     
c     Collision option 14 - set to zero for regions of low 
c     collisionality - specified in input as a fraction of SMAX - 
c     ideally less than 0.5.
c     
      elseif (cioptb.eq.14.and.(s.gt.(cstgrad*smax).and.
     >        s.lt.(1.0-cstgrad)*smax)
     >        ) then
         spara  = 0.0
         vpara  = 0.0
c     
c     Collision options 12, 13 and 14 (when 14 is not set to zero)
c     
c     This sqrt may prove to be computationally too intensive.
c     Previous implementation - see TAU module.
c     
c     >              *sqrt( LFPS/(2.0*TEMI))
c     
c     
      elseif (cioptb.eq.12.or.cioptb.eq.13.or.cioptb.eq.14) then
c     
c     Collision option 12 ... changes on each time-step
c     Collision option 13 ... tau para / (1+Mb/Mi)
c     
         spara  = 0.0
         vpara = lllfps
 7702    nrand = nrand + 1
         call surand2(seed,1,ran1)
         if (ran1.eq.0.0) goto 7702
         nrand = nrand + 1
         call surand2(seed,1,ran2)
         rgauss = sqrt(-2.0* log(ran1))*cos(2.0*PI*ran2)
         vpara = vpara * rgauss
c     
      endif

      return
      end
c     
c     
c     
      subroutine set_drift_velocity(s,ir,bg_drftvel,imp_drftvel)
      implicit none
      real bg_drftvel,imp_drftvel
      include    'params'
c     
      include    'driftvel'
c     
      real s
      integer ir
c     
C-----------------------------------------------------------------------
c     
c     Set the value of the drift velocity displacement 
c     depending on the particles position on the ring.
c     
c     Drifts only apply in the main SOL or PFZ at the moment
c     
c     Assign zero values as default
c     
      imp_drftvel = 0.0
      bg_drftvel = 0.0
c     
c     Replace with assigned values (which may be zero) depending
c     on options - option 1,3 are direct to the impurity particle -
c     option 2 is assigned as a background flow component and coupled
c     to the impurity through friction
c     
c     Region limitations are enforced in the setup_drftv routine
c
c     sdrft_start and sdrft_end are now specified on a ring by 
c     ring basis
c     
      if ((cpdrft.eq.1..or.cpdrft.eq.3).and.
     >     (s.ge.sdrft_start(ir).and.s.le.sdrft_end(ir))) then 
         imp_drftvel = pol_drftv(ir)
         bg_drftvel = 0.0
      elseif (cpdrft.eq.2.and.
     >        (s.ge.sdrft_start(ir).and.s.le.sdrft_end(ir))) then 
         imp_drftvel= 0.0
         bg_drftvel = pol_drftv(ir)
      endif

c     
c     Previous code
c     
c     if (cpdrft.eq.0) then 
c     imp_drftvel = 0.0
c     bg_drftvel = 0.0
c     elseif (ir.ge.irsep.and.cpdrft.eq.1.and.
c     >       (s.ge.sdrft_start.and.s.le.sdrft_end)) then 
c     imp_drftvel = pol_drftv(ir)
c     bg_drftvel = 0.0
c     elseif (ir.ge.irsep.and.cpdrft.eq.2.and.
c     >       (s.ge.sdrft_start.and.s.le.sdrft_end)) then 
c     imp_drftvel= 0.0
c     bg_drftvel = pol_drftv(ir)
c     
c     Limit Drift to PFZ 
c     
c     elseif (cpdrft.eq.3.and.ir.ge.irtrap.and.ir.le.nrs.and.
c     >       (s.ge.sdrft_start.and.s.le.sdrft_end)) then 
c     imp_drftvel = pol_drftv(ir)
c     bg_drftvel = 0.0
c     else
c     imp_drftvel = 0.0
c     bg_drftvel = 0.0
c     endif
c     


      return
      end
c     
c     
c     
      subroutine check_target_impact(seed,nrand,neutim)
      implicit none
c     
      real*8 seed
      integer nrand
      real neutim 
c     
      include    'params'
      include    'comtor'
      include    'cgeom'
c     
      include 'div1'
      include 'div4'
      include 'div5'
c     
      include    'particle_specs'



      integer   verify_id
      external  verify_id



C     
C--------CHECK IF REACHED TARGET
C     
      
      IF (S.LE.0.0 .OR. S.GE.SMAX) THEN
c     
c------------------------------------------------------
c     MIRROR TARGET OFF - TARGET IMPACT FOUND
c------------------------------------------------------
c     
         if ( cmiropt.eq.0.or.
     >        (s.le.0.0.and.cmiropt.eq.3).or.
     >        (s.ge.smax.and.cmiropt.eq.4)) then
c     
c     Set coordinate to zero and determine target element 
c     of impact 
c     
            if (cprint.eq.9) then 
               write(6,'(a,2i6,4(2x,g12.5))') 
     >                   'Update_walldep1: Struck target',
     >                   ik,ir,s,smax,ksmaxs(ir)
            endif
c
            IF (S.LE.0.0) THEN
               S  = 0.0
               IK = 1
               id = verify_id(ik,ir,2)
            ELSEif (s.ge.smax) then 
               S  = SMAX
               IK = NKS(IR)
               id = verify_id(ik,ir,1)
            ENDIF
c     
c     Determine postion on target/initial position options
c     
            if (init_pos_opt.eq.0) then
c     
               R = RP(ID)
               Z = ZP(ID)
c     
            elseif (init_pos_opt.eq.1) then 
c     
               call position_on_target(r,z,cross,id)
c     
            endif

c     
c     Check for direct ion to neutral reflection
c     
            call ion_neutral_reflection(seed,nrand,neutim)
c     
c     
C     IF Ion has not been reflected as a neutral - then call the target strike code
c     
            if (.not.reflect_ion) then


               if (cprint.eq.9) then
                  write(6,'(a,1i6,6x,3(2x,g12.5))') 
     >                   'Update_walldep2: Struck target',
     >                   id,cross,r,z
               endif


               call struck_target

c     Endif for test of target ion/neutral mirror
c     
            endif



c     
c----------------------------------------------------------------------------
c     NIRROR TARGET - ON
c----------------------------------------------------------------------------
c     
c     
c     Mirror target specified - MIRROR TARGET ON
c     
         elseif (cmiropt.eq.1.or.cmiropt.eq.3.or.cmiropt.eq.4) then
c     
c     Reverse sign of the particle velocity
c     To mimic target reflection/impact
c     
            vel = -vel
c     
            if (s.lt.0.and.(cmiropt.eq.1.or.cmiropt.eq.4)) then
               s = -s
            elseif (s.gt.smax.and.(cmiropt.eq.1.or.cmiropt.eq.3)) then
               s = smax - (s-smax)
c     
c     For the mirror target - it is necessary to
c     move the particle just slightly out from the
c     target in order to avoid an infinite loop
c     when S is exactly 0.0 or smax.
c     
            elseif (s.eq.0.0.and.(cmiropt.eq.1.or.cmiropt.eq.4)) then
               if (kss(ik,ir).eq. 0.0) then
                  s = kss(ik+1,ir) / 100.0
               else
                  s = kss(ik,ir) / 100.0
               endif
            elseif (s.eq.smax.and.(cmiropt.eq.1.or.cmiropt.eq.3)) then
               if (kss(ik,ir).eq.smax) then
                  s = smax - (smax-kss(ik-1,ir))/100.0
               else
                  s = smax - (smax-kss(ik,ir))/100.0
               endif
            endif
c     
c     Mirror target specified - confine particle - do not reflect 
c     velocity
c     
         elseif (cmiropt.eq.2) then
c     
c     Reverse sign of the particle velocity
c     To mimic target reflection/impact
c     
c     vel = -vel
c     
            if (s.lt.0) then
               s = -s
            elseif (s.gt.smax) then
               s = smax - (s-smax)
c     
c     For the mirror target - it is necessary to
c     move the particle just slightly out from the
c     target in order to avoid an infinite loop
c     when S is exactly 0.0 or smax.
c     
            elseif (s.eq.0.0) then
               if (kss(ik,ir).eq. 0.0) then
                  s = kss(ik+1,ir) / 100.0
               else
                  s = kss(ik,ir) / 100.0
               endif
            elseif (s.eq.smax) then
               if (kss(ik,ir).eq.smax) then
                  s = smax - (smax-kss(ik-1,ir))/100.0
               else
                  s = smax - (smax-kss(ik,ir))/100.0
               endif
            endif
         endif
      ENDIF


      return
      end
c     
c     
c     
      subroutine ion_neutral_reflection(seed,nrand,neutim)
      implicit none
c     
      real*8  seed
      real    neutim
      integer nrand
c     
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'cneut2'
c slmod begin
      include    'dynam3'
c slmod end
c     
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div4'
      include 'div5'
      include 'div6'
c     
      include    'particle_specs'
c
c     Output velocity along the field line from Launch_one (m/s)
c
      real vout
c slmod begin
      logical code_warning
      data    code_warning / .true. /
c slmod end
c     
c     ---       SPECIAL REFLECTION AS NEUTRAL PARTICLE  ----
c     
c     If the self-sputtering yield for this segment is 
c     set to -99.0 then the ion is converted to a neutral and 
c     then reflected instead of being allowed to strike and
c     sputter at the target. It is treated in the same 
c     manner as recombined neutrals except that the vel/ang
c     flag used for the launcg is different.   
c     
      write (6,*) 'REF:',ik,ir,id,kmfss(id)
c     
c     Set default condition
c     
      reflect_ion = .false.
c     
      if (kmfss(id).le.-99.0) then 
c     
c     Probability of reflection starts at 0.0 for a value
c     of -99.0 and rises to 1.0 for a value of -100.0 or more
c     
         refprob = min(abs(kmfss(id)+99.0),1.0) 
c     
c     Select random number
c     
         NRAND = NRAND + 1
         CALL SURAND2 (SEED, 1, RAN)
c     
         if (ran.le.refprob) reflect_ion = .true.
c     
      endif
c     
c     Test for reflection
c     
      if (reflect_ion) then  
c     
c     Call Launch_one to launch a single reflected neutral -
c     then branch to the appropriate location depending 
c     on the result.    
c     
c slmod begin
        if (code_warning) then 
           write(0,*) 
           write(0,*) 'WARNING ion_neutral_reflection: updating DEPS '//
     .                'for target reflection, which is non-standard'
           write(0,*) 
           code_warning = .false.
        endif
        if (id.lt.1.or.id.gt.nds) then    
           write (6,*) 'DEPS Error (ion_neutral_reflection):',id,iz,
     .                 sputy
        else 
           tneut = tneut + 1
           deps(id,iz) = deps(id,iz) + sputy
        endif
c slmod end
c     
c     Follow reflected impurities
c     
         ikorg = ik
         irorg = ir
c     
         write(6,*) ' debug: launch_one from ion_parallel_transport',id
         call LAUNCH_ONE (IMP,R,Z,RIZPOS,ZIZPOS,id,iwstart,
     >        rc,ctem1,cist,sputy,
     >        refSTRUK,mtcrefstruk,refMAIN,refEXIT,
     >        refATIZ,refNEUT,refWALLN,mtcrefwalln,
     >        refCENT,refTMAX,
     >        SEED,NRAND,
     >        NEUTIM,refFAIL,7,vout,vrec,refloss)
c     
c     
c
c        Set vel to appropriate value by scaling the along the field
c        line vout value returned by launch_one
c
         vel = vout * qtim
c
c slmod begin
         write (6,'(a,2i8,6g12.5)') 'REFLECT IMP:',imp,rc,
     >        r,z,rizpos,zizpos,
     >        temi,cist
c slmod end
c     
c     For all results other than reionization to state 1 - the code
c     will stop following the particle at this point and exit
c     as it would have for the old recombination implementation.
c     
         if (rc.ne.5) then 
            ifate = 11
            return
c     
c     goto 790
c     
         endif

c     
c     Deal with particle that was re-ionized
c     
         r = rizpos
         z = zizpos 
         iz= 1
c     
c     Look near last recorded ik,ir
c     
         call gridpos(ik,ir,r,z,.false.,griderr)
c     
c     If a griderr then exit anyway and issue error message 
c     since this shouldn't happen. 
c     
         if (griderr) then 
            
            write (6,*) 'PTR ION ERROR: NOT ON GRID:',r,z
            ifate = 11
            return
c     
c     goto 790
c     
         else
c     
c     If on grid - re-assign S-value - to cell centre if 
c     it left the cell - otherwise leave as is. 
c     
            if (ik.ne.ikorg.or.ir.ne.irorg) then 
c
c
c              Also need to reset the SMAX value for the new ring
c
               smax = ksmaxs(ir) 
c     
               if (init_pos_opt.eq.0) then
                  s = kss(ik,ir)
                  cross = 0.0
               elseif (init_pos_opt.eq.1) then
c     
                  call getscross_approx(r,z,s,cross,ik,ir)
c     
               endif 
c     
            endif   
c     
         endif 
c slmod begin
        write(6,*) 'debug: done launching reflected neutral'
c slmod end
      endif


      return
      end
c     
c     
c     
      subroutine struck_target
      implicit none
      include    'params'
      include    'dynam3'
      include    'comtor'
      include    'cgeom'
      include    'commv'
      include    'cneut'
      include    'cneut2'
c     
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div5'
      include 'div6'
c     
      include    'particle_specs'



      real      yield
      external  yield

c     
c     ION is NOT reflected - continue as normal
c     
C     
C     WRITE(6,*) 'SPUTTERED:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID)
C     
      CICABS(IZ) = CICABS(IZ) + SPUTY
c     CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c     CILABS(IZ) = MAX (CILABS(IZ), CIST)
      CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
      CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
      CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c     
c     Recording addition to total elapsed time in state IZ
c     
      cieizs(iz) = cieizs(iz) + cistiz * sputy
      citizs(iz) = citizs(iz) + sputy
c     
      CRTABS(IZ) = CRTABS(IZ) + TEMI * SPUTY
      CRVABS(IZ) = CRVABS(IZ) + VEL * SPUTY
      CRAVAV(IZ) = CRAVAV(IZ) + ABS(VEL) * SPUTY
      CTBS  (IZ) = CTBS  (IZ) + KTEBS(IK,IR) * SPUTY
      IM         = MIN (INT(TEMI/(0.2*CTEB0))+1, 10)
      CTEXS(IM)  = CTEXS(IM) + TEMI * SPUTY
c     
      acttarg(iz) = acttarg(iz) + sputy
C     
      RDEP   = RDEP + SPUTY

      if (id.lt.1.or.id.gt.nds) then    
         write (6,*) 'DEPS Error:',id,iz,sputy
      else 
         DEPS(ID,IZ) = DEPS(ID,IZ) + SPUTY
         NEROS(ID,1) = NEROS(ID,1) + SPUTY
      endif

      ENERGY = 3.0 * RIZ * KTEBS(IK,IR) +
     >     5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c     
      if (kmfss(id).ge.0.0) then  
         RYIELD = YIELD (6, MATTAR, ENERGY,
     >        ktebs(ik,ir),ktibs(ik,ir)) * KMFSS(ID)
      elseif (kmfss(id).lt.0.0.and.kmfss(id).ge.-50.0) then 
         RYIELD = abs(KMFSS(ID))
      elseif (kmfss(id).le.-99.0) then 
         RYIELD = YIELD (6, MATTAR, ENERGY,
     >        ktebs(ik,ir),ktibs(ik,ir))
      endif
c     
c     
      SPUNEW = SPUTY * RYIELD
      YLDTOT = YLDTOT + SPUNEW
      YLDMAX = MAX (YLDMAX, SPUNEW)
c     
c     WRITE(6,*) 'SPUTTERED:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID),
c     >                   ryield,kmfss(id),
c     >                    mattar,energy,sputy,spunew,yldtot
C     
      IF (SPUNEW.GT.CTRESH) THEN
         YTHTOT = YTHTOT + SPUNEW
         NPROD  = NPROD + 1
         SNEWS(NPROD) = SPUNEW
         IF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.
     >        CNEUTD.EQ.4) THEN
            EMAX = CEMAXF * ENERGY
            RMAXS(NPROD) = 1.0 / (1.0+CEBD/EMAX)**2
         ELSE
            RMAXS(NPROD) = 1.0
         ENDIF
         XPRODS(NPROD) = R
         YPRODS(NPROD) = Z
c     
c        For segments with a fixed sputtering yield - allow for
c        the energy of the sputtered particle to be set to a
c        specific value. 
c     
         if(kmfss(id).lt.0.0.and.cselfs.eq.2) then
            eprods(nprod) = ctem1
         else
            eprods(nprod) = 0.0
         endif
c     
         IDPRODS(NPROD) = ID
         launchdat(nprod,2) = 0.0
      ENDIF
c     
c     Add ion weight to wall element closest to grid 
c     departure.
c     
      if (cprint.eq.9) then 
         write(6,'(a)') 'Update_walldep3: '//
     >         'struck target'
      endif

      call update_walldep(ik,ir,iz,id,0,iwstart,idtype,sputy,energy)
c     
c     
      IFATE = 2
!     ammod begin.	       
      ! WBC comparison addition for ion prompt deposition.
      call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
!     ammod end.	       
c     
c     GOTO 790
c     

      return
      end
c     
c     
c     
      subroutine save_force_data(dvpara)
      implicit none
      real dvpara
c     
      include    'params'
      include    'comtor'
      include    'reiser_com' 
c     
      include 'div1'
      include 'div2'
c     
      include    'particle_specs'



c     
c     psmod
c     
C-----------------------------------------------------------------------
c     
c     Record data on average forces acting on impurity particles
c     
C-----------------------------------------------------------------------
c     
      Ftotal = ((FE+FEG+FIG+FF+FVG+dvpara)/QTIM2) * CRMI*AMU
c     
      Fcell(IK,IR,IZ) = Fcell(IK,IR,IZ) + Ftotal * SPUTY
c     
c     Examine Frictional K11(FF), Thermal K12(FIG), 
c     and Background Velocity Gradient K13 (FVG) 
c     contributions seperately
c     
      Ffi(IK,IR,IZ)  = Ffi(IK,IR,IZ)  + FF/QTIM2*CRMI*AMU * SPUTY
      Fthi(IK,IR,IZ) = Fthi(IK,IR,IZ) + FIG/QTIM2*CRMI*AMU* SPUTY
      Fvbg(IK,IR,IZ) = Fvbg(IK,IR,IZ) + FVG/QTIM2*CRMI*AMU* SPUTY
      DIFF(IK,IR,IZ) = DIFF(IK,IR,IZ) + dvpara/QTIM2*CRMI*AMU*SPUTY
c     
c     Calculate the total ion velocity per grid cell per charge state
c     
      VELavg(ik,ir,iz) = VELavg(ik,ir,iz) + VEL * SPUTY
c     
c     psmod
c     

      return
      end
c
c
c
      subroutine calculate_theta(ik,ir,s,theta)
      implicit none
      integer ik,ir
      real s,theta

      include 'params'
      include 'cgeom'
c
c     CALCULATE_THETA: Calculate the theta value associated with 
c                      the given S in the specified cell. 
c
c     NOTE: For the core rings the first and last cells are the same.
c     S=0 starts at the cell centre and S=SMAX is at the same
c     cell centre - thus - in the core - the particle should never
c     have an S-value with S < KSS(1,ir) or S > KSS(nks(ir),ir)
c     
      IF (S.GT.KSS(IK,IR)) THEN
         IF (IK.LT.NKS(IR).or.ir.lt.irsep) THEN
            THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
     >           *(THETAG(IK+1,IR)-THETAG(IK,IR))
         ELSE
            THETA = THETAG(IK,IR) + (S-KSS(IK,IR))/KFORDS(IK,IR)
     >           *(THETAT(IDDS(IR,1))-THETAG(IK,IR))
         ENDIF
      ELSE
         IF (IK.GT.1.or.ir.lt.irsep) THEN
c     
c     This was added to cover the case of S=0 occuring for an ion initially
c     injected in the first cell of a core ring.  
c     
            if (s.eq.kss(ik,ir)) then
               theta = thetag(ik,ir)
            else
               THETA = THETAG(IK,IR) - (KSS(IK,IR)-S)/KBACDS(IK,IR)
     >              *(THETAG(IK,IR)-THETAG(IK-1,IR))
            endif  
         ELSE
            THETA = THETAG(IK,IR) - (KSS(IK,IR)-S)/KBACDS(IK,IR)
     >           *(THETAG(IK,IR)-THETAT(IDDS(IR,2)))
         ENDIF
      ENDIF

      return
      end

