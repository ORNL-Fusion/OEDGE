c     -*-Fortran-*-
c
      subroutine execute_transport_step(seed,nrand,neutim,
     >                                  ero_record_data)
      use ero_interface
      implicit none
c
      real*8  seed
      real    neutim
      integer nrand
      logical ero_recorded
c
      include    'params'
      include    'comtor'
      include    'cgeom'
      include    'commv'
c
      include 'div1'
      include 'div3'
      include 'div5'
      include 'div6'
c
      include    'particle_specs'
      include    'driftvel'
c
      logical ero_record_data
      real spara,dspara,vpara,dvpara
c
c     Record current cell and position 
c
      ir_last    = ir
      ik_last    = ik
      s_last     = s
      cross_last = cross
      theta_last = theta
c
c     Execute parallel step
c
      call do_parallel_step(seed,nrand,neutim,spara,dspara,vpara,dvpara)

c
      if (ifate.ne.0) return 
c
c     Update the Cross-field transport term and make 
c     all other related position adjustments.
c
c     Execute cross-field step
c
c
c     Record position
c
      IKOLD = IK
      IROLD = IR
      oldcross = cross
      oldtheta = theta 
c
      call do_crossfield_step(ik,ir,ikold,irold,kk,s,theta,cross,
     >                        oldtheta,oldcross,
     >                        adjust,dcross,ckkmin,smax,k,debug,
     >                        seed,nrand,neutim,cist,imp,debug_all,
     >                        ifate)
c
c
      if (ifate.ne.0) return 
c
c
c     sdrft_start and sdrft_end are now setup on a ring by ring 
c     basis in the setup_drftvel routine
c
c      if (ir.ne.ir_last) then 
c
c           Set range for poloidal drift velocity
c 
c            if (cpdrft.ne.0) then 
c               sdrft_start = cdrftv_start*smax 
c               sdrft_end   = cdrftv_end*smax 
c            endif
c      endif
c
c
c     Check for particle entering the peripheral plasma
c
      call check_reached_grid_edge(seed,nrand)
c
      if (ifate.ne.0) return 
c
c
c     Record statistics and other particle data - execute code specific to individual regions
c

      if (ir.lt.irsep) then  
         call ion_in_main(spara,dspara,vpara,dvpara)
      else 
         call ion_in_sol(spara,dspara,vpara,dvpara)
      endif
c
      if (ifate.ne.0) return
c
c     Find current R,Z location 
c
      call findrz
c
c     Monitor for debug purposes 
c
      if (smax.ne.ksmaxs(ir)) then 
         write(6,'(a,2i6,10(1x,g12.5))') 'SMAX ERROR:',ik,ir,s,smax,
     >         ksmaxs(ir)
         write(0,'(a,2i6,10(1x,g12.5))') 'SMAX ERROR:',ik,ir,s,smax,
     >         ksmaxs(ir)
      endif
c
c     Monitor particle for entering ERO simulation volune - after update to R,Z
c
c     Variable are called: R,Z, SPUTY, IZ, VEL, TEMI - 
c     
      if (ero_record_data) then 
         call ero_check_part_track_ion(ik,ir,iz,r,z,vel,qtim,temi,
     >                             sputy,ero_record_data)
      endif


c        
      return
      end   




c
c
c
      subroutine ion_in_main(spara,dspara,vpara,dvpara)
      implicit none
      real spara,dspara,vpara,dvpara
c
      include    'params'
      include    'dynam3'
      include    'comtor'
      include    'cgeom'
      include    'commv'
      include    'cneut'
      include    'clocal'
c
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div4'
      include 'div5'
      include 'div6'
      include 'div7'
c
      include    'particle_specs'

      real za02as
      external za02as




c
c         Record particle has entered core plasma - count it and save
c         it's starting ionization position. (for possible scatter plot)
c
          if (.not.hasleakedcore) then 
             hasleakedcore = .true.
             nleakcore = nleakcore + 1 
             totleakcore = totleakcore + sputy
             if (nleakcore.le.maximp) then 
                cleakpos(nleakcore,1) = rstart
                cleakpos(nleakcore,2) = zstart
             endif
          endif


c
c         Record average parallel steps for debugging purposes.
c
          if (dvpara.ne.0.0) then
             if (dvpara.lt.0.0) then
                dvparastep(5,iz) = dvparastep(5,iz) + dvpara
                vparastep(5,iz) = vparastep(5,iz) + abs(vel)
                dvparacnt(5,iz)  = dvparacnt(5,iz) + 1.0
                dvmaxv(5,iz)   = max(dvmaxv(5,iz),dble(vel))
                dvminv(5,iz)   = min(dvminv(5,iz),dble(vel))
             elseif (dvpara.gt.0.0) then
                dvparastep(6,iz) = dvparastep(6,iz) + dvpara
                vparastep(6,iz) =  vparastep(6,iz) + abs(vel)
                dvparacnt(6,iz)  = dvparacnt(6,iz) + 1.0
                dvmaxv(6,iz)   = max(dvmaxv(6,iz),dble(vel))
                dvminv(6,iz)   = min(dvminv(6,iz),dble(vel))
             endif
          endif
c
c         Record some quantities in the core
c
          COREOUTS(IZ,1) = COREOUTS(IZ,1) + DSPUTY
          COREOUTS(IZ,2) = COREOUTS(IZ,2) + DSPUTY *
     >                      TEMI / LFPS(IK,IR,IZ)
          COREOUTS(IZ,3) = COREOUTS(IZ,3) + DSPUTY / 
     >                      LFSS(IK,IR,IZ)
          IF (S.LE.0.5*SMAX) THEN
            COREOUTS(IZ,4) = COREOUTS(IZ,4) + DSPUTY * FF
            COREOUTS(IZ,5) = COREOUTS(IZ,5) + DSPUTY * FE
            COREOUTS(IZ,6) = COREOUTS(IZ,6) + DSPUTY * FEG
            COREOUTS(IZ,7) = COREOUTS(IZ,7) + DSPUTY * FIG
            COREOUTS(IZ,8) = COREOUTS(IZ,8) + DSPUTY * FVEL
            COREOUTS(IZ,9) = COREOUTS(IZ,9) + DSPUTY * FVH
          ELSE
            COREOUTS(IZ,4) = COREOUTS(IZ,4) - DSPUTY * FF
            COREOUTS(IZ,5) = COREOUTS(IZ,5) - DSPUTY * FE
            COREOUTS(IZ,6) = COREOUTS(IZ,6) - DSPUTY * FEG
            COREOUTS(IZ,7) = COREOUTS(IZ,7) - DSPUTY * FIG
            COREOUTS(IZ,8) = COREOUTS(IZ,8) - DSPUTY * FVEL
            COREOUTS(IZ,9) = COREOUTS(IZ,9) - DSPUTY * FVH
          ENDIF


c
c          Accumulate statistics if particle has
c          not entered the core plasma previously. 
c
           if (.not.incore) then
              incore = .true. 
c
              num_entered_core = num_entered_core + sputy
c
              ncore(ikstart,irstart) = ncore(ikstart,irstart) + sputy
c
              if (idtype.gt.0) then
                  wtsource(idstart,irstart,3,idtype) =
     >                       wtsource(idstart,irstart,3,idtype) + sputy
                  wtsource(idstart,irstart,4,idtype) =
     >               wtsource(idstart,irstart,4,idtype)+sputy*sstart
              endif
c
c
           endif
c
c
c
           IF (CFLRIN.and.ir_last.ge.irsep) THEN

c nonorth
            IF (DEBUGL) WRITE(6,9003) IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >        THETA,SMAX,VEL,TEMI,SPARA,CROSS,SPUTY,IT,'ENTERED MAIN'
c nonorth
c            IF (DEBUGL) WRITE(6,9003) IMP,CIST,IK,IR,IZ,
c     >        R,Z,S,K,SMAX,VEL,TEMI,SPARA,CROSS,SPUTY,IT,'ENTERED MAIN'

c
c           Current position statistics upon core entry
c
            CNNN(IZ)  = CNNN(IZ)  + SPUTY
            CNNNX(IZ) = CNNNX(IZ) + Z * SPUTY
            CNNNS(IZ) = CNNNS(IZ) + MIN (S,SMAX-S) * SPUTY
c
c           Original ionization position statistics upon core entry
c
            CNNNK(IZ) = CNNNK(IZ) + KATIZS(IMP) * SPUTY
            cnorgs(iz) = cnorgs(iz) + satizs(imp) * sputy
            cnorgr(iz) = cnorgr(iz) + xatizs(imp) * sputy
            cnorgz(iz) = cnorgz(iz) + yatizs(imp) * sputy
c
c           Accumulated statistics on original position 
c
            CVVRM  = CVVRM + XATIZS(IMP) * SPUTY
            CVVZM  = CVVZM + YATIZS(IMP) * SPUTY
            CVVKM  = CVVKM + KATIZS(IMP) * SPUTY
            CVVSM  = CVVSM + SATIZS(IMP) * SPUTY
c
c           Accumulate some statistics on different sources
c
c           Original Neutral from FP launch
c
            if (launchdat(imp,2).eq.1.0) then  
c
c              Original neutral refected.
c
               if (launchdat(imp,3).ne.0.0) then 
                  cvvfpref = cvvfpref + sputy
               else 
                  cvvfpnrf = cvvfpnrf + sputy
               endif 
            else
c
c              Original neutral refected.
c
               if (launchdat(imp,3).ne.0) then 
                  cvvrefm = cvvrefm + sputy
               else
                  cvvnrfm = cvvnrfm + sputy
               endif  
            endif
c
            CNNNT(IZ) = CNNNT(IZ) + QTIM * CIST * SPUTY
            CNNNKT(IZ)= CNNNKT(IZ)+ KKK / CIST * SPUTY
c
            ELIMS(IK,3,IZ) = ELIMS(IK,3,IZ) + SPUTY
c
            CICRIN = CICRIN + SPUTY
            CISRIN = CISRIN + CIST * SPUTY
            CITRIN = CITRIN + TEMI * SPUTY
            CIKRIN = CIKRIN + KATIZS(IMP) * SPUTY
            CKTRIN = CKTRIN + KKK / CIST * SPUTY
c            CIFRIN = MIN (CIFRIN, CIST)
            CIFRIN = MIN (CIFRIN, sngl(CIST))
            CFLRIN = .FALSE.
c
            CVVXE  = CVVXE + Z * SPUTY
            IF (K.GT.KATIZS(IMP)) THEN
              XTRIPP  = XTRIPP + (Z - OLDZ)
            ELSE
              XTRIPS  = XTRIPS + (Z - OLDZ)
            ENDIF
            CVVXP  = CVVXP + XTRIPP * SPUTY
            CVVXS  = CVVXS + XTRIPS * SPUTY
c slmod begin 
            IF (CPRINT.GT.1)
     >        WRITE (6,9022) Z,YATIZS(IMP),XTRIPP,XTRIPS,
     >              100.*XTRIPP/max(1e-8,(Z-YATIZS(IMP))),
     >              100.*XTRIPS/max(1e-8,(Z-YATIZS(IMP))),IMP,SPUTY,
     >           ZA02AS (1) - STATIM
c
c            WRITE (6,9022) Z,YATIZS(IMP),XTRIPP,XTRIPS,
c     >            100.*XTRIPP/max(1e-8,(Z-YATIZS(IMP))),
c     >            100.*XTRIPS/max(1e-8,(Z-YATIZS(IMP))),IMP,SPUTY,
c     >         ZA02AS (1) - STATIM
c slmod end
           ENDIF
c
c          ENDIF
c
c
c        Stop following particle in main - if option is set. 
c
c        If the option is turned on to stop following the particle after core
c        entry then set the ifate value before exiting this routine
c
         IF (CSTOP.EQ.1) THEN

            stopped_follow = stopped_follow+sputy             
            IFATE = 7
c
            return
c
c           GOTO 790
c
         ENDIF
c


c
C
C-------- REFLECT OFF CENTRAL MIRROR
C
          IF (IR.le.ircore) THEN
            CROSS = -ABS(CROSS)
c
            if (ir.lt.ircore) then
c
c               write (6,*) 'WARNING: Adjusting particle'//
c     >           ' found inside central mirror:',ik,ir,cross,r,z,s
c
c              Set to one step outside core mirror ring. 
c
               ir = ircore
               cross = -distout(ik,ircore) - kperps(ik,ircore)
c
            endif 
c
            IF (CFLRXA) THEN
              CICRXA = CICRXA + SPUTY
              CISRXA = CISRXA + CIST * SPUTY
              CITRXA = CITRXA + TEMI * SPUTY
c              CIFRXA = MIN (CIFRXA, CIST)
              CIFRXA = MIN (CIFRXA, sngl(CIST))
              CFLRXA = .FALSE.
c
              IF (DEBUGL) WRITE(6,9003)IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >          THETA,SMAX,VEL,TEMI,SPARA,CROSS,SPUTY,IT,'REFLECTED'
c
            ENDIF
          ENDIF













 9003 FORMAT(1X,I5,F9.1,2I3,I2,2F9.5,F8.3,2F6.2,F8.3,1P,E15.8,
     >  0P,F7.1,1P,E8.1,0P,F8.5,F5.2,I2,:,1X,A,:,F8.5)

 9022 FORMAT(1X,'DIV: ZENTRY',F6.3,', ZCREAT',F6.3,', ZTRIPP',F6.3,
     >  ', ZTRIPS',F6.3,', %P',F7.1,', %S',F7.1,'  (ION',I5,
     >  '  WEIGHT',F5.2,')',' TIME:',f10.2)


      return
      end
c
c
c 
      subroutine ion_in_sol(spara,dspara,vpara,dvpara)
      use divertor_limits
      implicit none
      real  spara,dspara,vpara,dvpara
c      
      include    'params'
      include    'dynam3'
      include    'comtor'
      include    'cgeom'
      include    'commv'
      include    'clocal'
c
      include 'div1'
      include 'div2'
      include 'div4'
      include 'div5'
      include 'div6'
      include 'div7'
c
      include    'particle_specs'
c
      include 'hc_global_opts'
c
      real tmp_time
c
      integer flag
c
c         Accumulate some data if the particle enters various regions
c         of the SOL or trapped plasma. (related to initial ionization
c         positions of particles entering various regions)
c
          if (.not.inedge) then            
             inedge = .true.
             nedge(ikstart,irstart) = nedge(ikstart,irstart) + sputy
          endif      
c
          if (.not.intrap.and.ir.gt.irwall.and.ir.le.nrs) then 
c
c            TRAP region 
c
             intrap = .true.
             ntrap(ikstart,irstart) = ntrap(ikstart,irstart) + sputy
c
          elseif (.not.indiv.and.ir.ge.irsep.and.ir.le.irwall.and. 
     >             ((z.ge.zxp.and.xpoint_up).or.
     >              (z.le.zxp.and.(.not.xpoint_up)))) then  
c
c            Divertor region 
c
             indiv = .true.
             ndivert(ikstart,irstart) = ndivert(ikstart,irstart)+sputy
c
          elseif (.not.inmsol.and.ir.ge.irsep.and.ir.le.irwall.and. 
     >             ((z.lt.zxp.and.xpoint_up).or.
     >              (z.gt.zxp.and.(.not.xpoint_up)))) then  
c
c            Main SOL Region 
c
             inmsol = .true.
             nmsol(ikstart,irstart) = nmsol(ikstart,irstart)+sputy
          endif






C-----------------------------------------------------------------------
c
c
c         Record average parallel steps for debugging purposes.
c
          if (dvpara.ne.0.0) then
             if (s.lt.smax/2.0) then
                if (dvpara.lt.0.0) then
                   dvparastep(1,iz) = dvparastep(1,iz) + dvpara
                   vparastep(1,iz)  = vparastep(1,iz) + abs(vel)
                   dvparacnt(1,iz)  = dvparacnt(1,iz) + 1.0
                   dvmaxv(1,iz)   = max(dvmaxv(1,iz),dble(vel))
                   dvminv(1,iz)   = min(dvminv(1,iz),dble(vel))
                elseif (dvpara.gt.0.0) then
                   dvparastep(2,iz) = dvparastep(2,iz) + dvpara
                   vparastep(2,iz)  = vparastep(2,iz) + abs(vel)
                   dvparacnt(2,iz)  = dvparacnt(2,iz) + 1.0
                   dvmaxv(2,iz)   = max(dvmaxv(2,iz),dble(vel))
                   dvminv(2,iz)   = min(dvminv(2,iz),dble(vel))
                endif
             elseif (s.gt.smax/2.0) then
                if (dvpara.lt.0.0) then
                   dvparastep(3,iz) = dvparastep(3,iz) + dvpara
                   vparastep(3,iz)  = vparastep(3,iz) + abs(vel)
                   dvparacnt(3,iz)  = dvparacnt(3,iz) + 1.0
                   dvmaxv(3,iz)   = max(dvmaxv(3,iz),dble(vel))
                   dvminv(3,iz)   = min(dvminv(3,iz),dble(vel))
                elseif (dvpara.gt.0.0) then
                   dvparastep(4,iz) = dvparastep(4,iz) + dvpara
                   vparastep(4,iz)  = vparastep(4,iz) + abs(vel)
                   dvparacnt(4,iz)  = dvparacnt(4,iz) + 1.0
                   dvmaxv(4,iz)   = max(dvmaxv(4,iz),dble(vel))
                   dvminv(4,iz)   = min(dvminv(4,iz),dble(vel))
                endif
             endif
          endif
c
c         Spara - spatial diffusion.
c
          if (dspara.ne.0.0) then
             if (s.lt.smax/2.0) then
                if (dspara.lt.0.0) then
                   dsparastep(1) = dsparastep(1) + dspara
                   dsparacnt(1)  = dsparacnt(1) + 1.0
                elseif (dspara.gt.0.0) then
                   dsparastep(2) = dsparastep(2) + dspara
                   dsparacnt(2)  = dsparacnt(2) + 1.0
                endif
             elseif (s.gt.smax/2.0) then
                if (dspara.lt.0.0) then
                   dsparastep(3) = dsparastep(3) + dspara
                   dsparacnt(3)  = dsparacnt(3) + 1.0
                elseif (dspara.gt.0.0) then
                   dsparastep(4) = dsparastep(4) + dspara
                   dsparacnt(4)  = dsparacnt(4) + 1.0
                endif
             endif
          endif



C
C-------- EXIT FROM MAIN PLASMA
C
           IF (CFLREX.and.ir_last.lt.irsep) THEN
c
            IF (DEBUGL) WRITE(6,9003) IMP,CIST,IK,IR,IZ,R,Z,S,K,
     >        THETA,SMAX,VEL,TEMI,SPARA,CROSS,SPUTY,IT,'EXITED MAIN'
c
            CFLREX = .FALSE.
            IF (INMAIN) THEN
              CMMM(IZ)  = CMMM(IZ) + SPUTY
              CMMMX(IZ) = CMMMX(IZ) + Z * SPUTY
              CMMMS(IZ) = CMMMS(IZ) + MIN (S, SMAX-S) * SPUTY
              ELIMS(IK,2,IZ) = ELIMS(IK,2,IZ) + SPUTY
            ELSE
              CLLL(IZ)  = CLLL(IZ) + SPUTY
              CLLLX(IZ) = CLLLX(IZ) + Z * SPUTY
              CLLLS(IZ) = CLLLS(IZ) + MIN (S, SMAX-S) * SPUTY
              ELIMS(IK,1,IZ) = ELIMS(IK,1,IZ) + SPUTY
            ENDIF
           ENDIF

 
          DOUTS(IZ,1) = DOUTS(IZ,1) + DSPUTY
          DOUTS(IZ,2) = DOUTS(IZ,2) + DSPUTY * TEMI / LFPS(IK,IR,IZ)
          DOUTS(IZ,3) = DOUTS(IZ,3) + DSPUTY / LFSS(IK,IR,IZ)
          IF (S.LE.0.5*SMAX) THEN
            DOUTS(IZ,4) = DOUTS(IZ,4) + DSPUTY * FF
            DOUTS(IZ,5) = DOUTS(IZ,5) + DSPUTY * FE
            DOUTS(IZ,6) = DOUTS(IZ,6) + DSPUTY * FEG
            DOUTS(IZ,7) = DOUTS(IZ,7) + DSPUTY * FIG
            DOUTS(IZ,8) = DOUTS(IZ,8) + DSPUTY * FVEL
            DOUTS(IZ,9) = DOUTS(IZ,9) + DSPUTY * FVH
          ELSE
            DOUTS(IZ,4) = DOUTS(IZ,4) - DSPUTY * FF
            DOUTS(IZ,5) = DOUTS(IZ,5) - DSPUTY * FE
            DOUTS(IZ,6) = DOUTS(IZ,6) - DSPUTY * FEG
            DOUTS(IZ,7) = DOUTS(IZ,7) - DSPUTY * FIG
            DOUTS(IZ,8) = DOUTS(IZ,8) - DSPUTY * FVEL
            DOUTS(IZ,9) = DOUTS(IZ,9) - DSPUTY * FVH
          ENDIF





c
c------ Check for parallel leakage - record particle once.
c
        if (checkleak.and.ir.ge.irsep
     >      .and.(.not.hasleaked.and.
     >      (s.gt.cleaks(cleakp).and.s.lt.(smax-cleaks(cleakp)))))then
           cleakn(cleakp,iz) = cleakn(cleakp,iz) + sputy
           cleakp = cleakp +1
c           write (6,*) 'Debug leak:',imp,iz,cleakp,s,cleaks(cleakp),
c     >                  smax-cleaks(cleakp),smax,cleakn(cleakp,iz)
           if (cleakp.gt.cleaksn) then
              cleakp = cleaksn
              hasleaked = .true.
              cleakt = cleakt + cist * qtim
           endif
        endif
c
c       Use the same mechanism for divertor leakage checks
c
        if (checkleak.and.(.not.divertor_leaked)) then 
           tmp_time = cist * qtim
           divertor_leaked = check_divertor_limit(sputy,ir,iz,
     >                                       tmp_time,s)
        endif



! ammod begin.
!       Check if ion is inside the WBC geometry boundary and count if not.
        if (global_hc_follow_option.ne.0) then 
c
           flag = 0
c
           call global_hc_check_WBC_Ion_Pos(IK,IR,S,CROSS,IZ,CRMI,
     >                VEL/qtim,TEMI,SPUTY,
     >                NIMPS_local,NIMPS2_local,NATIZ,IMP,flag)
c
!          Add to WBC counting and stop following. 
           if (flag.eq.1) ifate=7

        endif
! ammod end.	       




 9003 FORMAT(1X,I5,F9.1,2I3,I2,2F9.5,F8.3,2F6.2,F8.3,1P,E15.8,
     >  0P,F7.1,1P,E8.1,0P,F8.5,F5.2,I2,:,1X,A,:,F8.5)



      return
      end

c
c
c
      subroutine findrz
      implicit none
      include    'params'
      include    'comtor'
      include    'cgeom'
c
      include    'particle_specs'

c
c     Set new approximate values of R,Z 
c


      IF (rzopt.gt.0) THEN
        CALL GETRZ(IK,IR,S,CROSS,R,Z,rzopt)
      ELSEif (rzopt.eq.0) then 
        R = RS(IK,IR)                                                 
        Z = ZS(IK,IR)                                                 
      ENDIF
c
c      R = RS(IK,IR)
c      Z = ZS(IK,IR)
c sl end
c


      return
      end


      subroutine check_reached_grid_edge(seed,nrand)
      implicit none
c 
      real*8 seed
      integer nrand
c
      include    'params'
      include    'dynam3'
      include    'dynam4'
      include    'comtor'
      include    'cgeom'
      include    'commv'
      include    'cneut'
      include    'cneut2'
      include    'fperiph_com'
c
      include 'div1'
      include 'div2'
      include 'div3'
      include 'div5'
      include 'div6'
c
      include    'particle_specs'

      integer fperiph
      external fperiph

      integer   verify_id
      external  verify_id

      real      yield
      external  yield

c
c     Variables for periphery option 5
c
      integer istate,id_out,is_out
      real rsect,zsect
c      real energy
c
C
C-------- CHECK IF REACHED WALLS
C
C         MODIFY TO SUPPORT ION REFLECTION INSTEAD OF ELIMINATION
C         IF THE OPTION IS ACTIVE.
C
C         OPTION 2 MOVES THE WALL TO THE OUTERMOST RING INSTEAD
C         OF MIDWAY BETWEEN THE LAST AND NEXT TO LAST.
C
          if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
            IKOLD = IK
            IROLD = IR
c nonorth
          endif

c
  770     IF (IR.EQ.IRWALL.OR.IR.EQ.IRTRAP.or.ir.eq.irwall2
     >        .or.ir.eq.irtrap2) THEN
            IF ( ((FPOPT.EQ.0.or.(fpopt.eq.4.and.ir.eq.irwall)).AND.
     >           (CIONR.EQ.0.or.cionr.eq.2)) .OR.
     >           ((FPOPT.EQ.0.or.(fpopt.eq.4.and.ir.eq.irwall)).AND.
     >            CIONR.EQ.1.AND.
     >           ( (CROSS.LE.0.0.AND. (IR.EQ.IRWALL.or.ir.eq.irwall2) )
     >           .OR.
     >           (CROSS.GE.0.0.AND. (IR.EQ.IRTRAP.or.ir.eq.irtrap2) ) )
     >           ))  THEN
              CICABS(IZ) = CICABS(IZ) + SPUTY
c              CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c              CILABS(IZ) = MAX (CILABS(IZ), CIST)
              CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
              CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
              CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c             Recording addition to total elapsed time in state IZ
c
              cieizs(iz) = cieizs(iz) + cistiz * sputy
              citizs(iz) = citizs(iz) + sputy
c
              CRTABS(IZ) = CRTABS(IZ) + TEMI * SPUTY
              CRVABS(IZ) = CRVABS(IZ) + VEL * SPUTY
              CRAVAV(IZ) = CRAVAV(IZ) + ABS(VEL) * SPUTY
              CTBS  (IZ) = CTBS  (IZ) + KTEBS(IK,IR) * SPUTY
              IM         = MIN (INT(TEMI/(0.2*CTEB0))+1, 10)
c slmod begin
c...          This problem has appeared on at least one occasion, but not very
c             often:
c TEMI is the problem
c it is being assigned properly in div.f after the neutral launch, so it's being
c reassigned in a nasty way somewhere along between there and here...
              IF (IM.LT.1) THEN
                CALL WN('check_reached_grid_edge','Array bounds '//
     .                  'violation, setting IM=1')
                 WRITE(0,*) 'TEMI  =',TEMI
                 WRITE(0,*) 'CTEB0 =',CTEB0
                 WRITE(0,*) 'RESULT=',INT(TEMI/(0.2*CTEB0))+1
                 WRITE(6,*) 'TEMI  =',TEMI
                 WRITE(6,*) 'CTEB0 =',CTEB0
                 WRITE(6,*) 'RESULT=',INT(TEMI/(0.2*CTEB0))+1
                 IM = 1
              ENDIF
c slmod end
              CTEXS(IM)  = CTEXS(IM) + TEMI * SPUTY
c
              RWALL      = RWALL + SPUTY
              WALLS(IK,IR,IZ) = WALLS(IK,IR,IZ) + SPUTY
c
              ENERGY = 
     >           5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c
c             Add ion weight to wall element closest to grid 
c             departure.
c
              
c              write(6,*) 'update_walldep: '//
c     >                   'check_reached_grid_edge - hard walls'

              call update_walldep(ik,ir,iz,0,0,iwstart,idtype,sputy,
     >                            energy)
c
              IFATE = 1

! ammod begin.	       
	       ! WBC comparison addition for ion prompt deposition.
               call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
! ammod end.	       

              return

c
c              GOTO 790
c
            ELSEIF ((FPOPT.EQ.1.or.(fpopt.eq.4.and.ir.eq.irtrap))
     >              .AND.CIONR.EQ.1
     >             .AND.
     >         ((CROSS.LE.0.0.AND.(IR.EQ.IRWALL.or.ir.eq.irwall2)
     >           ).OR.
     >          (CROSS.GE.0.0.AND.(IR.EQ.IRTRAP.or.ir.eq.irtrap2)
     >          ))) THEN
C
C             REFLECT ION AND CONTINUE.
C             ALLOW ION TO CONTINUE IT WILL NEVER BE FARTHER THAN
C             THE OUTERMOST RING.
C
              CROSS = -CROSS
            ELSEIF ((FPOPT.EQ.1.or.(fpopt.eq.4.and.ir.eq.irtrap))
     >              .AND.(CIONR.EQ.0.or.cionr.eq.2))THEN
C
C             REFLECT ION AND CONTINUE.
C             ALLOW ION TO CONTINUE IT WILL NEVER BE FARTHER THAN
C             THE INNER EDGE OF THE OUTERMOST RING REGION.
C             FOR THETA ALLOW FOR (REMOTE) POSSIBILITY OF CROSSING
C             SEPARATRIX IN THE INNER DIVERTOR REGION.
C
              CROSS = -CROSS
              IF (IR.EQ.IRWALL) THEN
                JK = IKINS(IK,IR)
                IR = IRINS(IK,IR)
c
                if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
                IF( IR.EQ.NRS
     >              .AND. JK.GE.IKTI )THETA = THETA - DTHETG
c nonorth
                endif
c
              ELSEIF (IR.EQ.IRTRAP) THEN
                JK = IKOUTS(IK,IR)
                IR = IROUTS(IK,IR)
c
                if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
                  IF( IR.EQ.IRSEP
     >                .AND. JK.GT.IKTO )THETA = THETA + DTHETG
c nonorth
                endif
c
              ENDIF
 
c nonorth
              IF (northopt.eq.0.or.northopt.eq.2.or.
     >          ((northopt.eq.1.or.northopt.eq.3).and.
     >          (TAGDV(JK,IR).EQ.0 .AND. TAGDV(IKOLD,IROLD).EQ.0))) THEN
C
C               BOTH POINTS UNSHIFTED FROM ORTHOGONAL GRID
C
c nonorth
                IK = JK
c nonorth
              ELSEif (northopt.eq.1.or.northopt.eq.3) then
C
C               EITHER POINT SHIFTED FROM ORTHOGONAL GRID
C
                IF (IROLD.EQ.IRWALL) THEN
                  JK = IKING(IK,IROLD)
                  IR = IRINS(IK,IROLD)
                ELSEIF (IROLD.EQ.IRTRAP) THEN
                  JK = IKOUTG(IK,IROLD)
                  IR = IROUTS(IK,IROLD)
                ENDIF
c
  771           IF( THETA.LT.THETAG(JK,IR) )THEN
                  IF( JK.GT.1 )THEN
                    THETA1 = (THETAG(JK,IR) - THETA)
     >                       /(THETAG(JK,IR) - THETAG(JK-1,IR))
                    IF( THETA1.GT.0.5 )THEN
                      JK = JK - 1
                      GOTO 771
                    ENDIF
                  ENDIF
                ELSE
                  IF( JK.LT.NKS(IR) )THEN
                    THETA1 = (THETA - THETAG(JK,IR))
     >                       /(THETAG(JK+1,IR) - THETAG(JK,IR))
                    IF( THETA1.GE.0.5 )THEN
                      JK = JK + 1
                      GOTO 771
                    ENDIF
                  ENDIF
                ENDIF
                IK = JK
              endif
c
c             jdemod
c
c             Update SMAX for this periphery option since the value of IR may have changed
c     
              smax = ksmaxs(ir)
c
c nonorth
c
c           Periphery option 3 
c 
            ELSEIF ((FPOPT.EQ.3.or.fpopt.eq.5.or.fpopt.eq.6).AND.
     >              ((CIONR.EQ.0.or.cionr.eq.2).OR.
     >               (CIONR.EQ.1.AND.
     >          ( (CROSS.LE.0.0.AND.(IR.EQ.IRWALL.or.ir.eq.irwall2))
     >          .OR.
     >          (CROSS.GE.0.0.AND.(IR.EQ.IRTRAP.or.ir.eq.irtrap2))
     >          ) ) )) THEN
C
C             FPOPT 3 AND WALL OPT 1 OR 3:
C
C             WALL AT LAST RING
C
C             WALL IS REPLACED BY THE FAR PERIPHERY REGION.
C             IONS HAVE 4 POSSIBLE FATES:
C             1) RE-ENTRY WHERE THEY EXITED - IF THEY DIFFUSE BACK.
C             2) DISCARDED BY REACHING TIME-LIMIT.
C             3) IMPACT WITH WALL AT DISTANCE FPXMAX.
C             4) FP LOSS TO TARGET PLATE IMPACT. USING
C                CHARACTERISTIC TIME FPTIM.
C
C             FPOPT 3 AND WALL OPT 0 OR 2:
C
C             AS ABOVE BUT PLACING WALL/PERIPHERY RING BETWEEN
C             THE LAST TWO RINGS. THIS IS NECESSARY BECAUSE THE
C             JET SHOT DATA DOES NOT CONTAIN INFORMATION
C             ON THE BACKGROUND CHARACTERISTICS FOR THE OUTERMOST
C             RINGS - DESPITE CONTAINING ALL THE GEOMETRY DATA.
C
C
              FPENT = FPENT + SPUTY
c             
c             Outer or Inner FP 
c

              if (s.lt.(ksmaxs(ir)/2.0)) then 
c                            
c                OUTER 
c
                 fpdist = fpxmaxO
                 fplosstim = qtim/fptimO 
c
              else
c
c                INNER          
c
                 fpdist = fpxmaxI
                 fplosstim = qtim/fptimI 
c
              endif
c
c             Initialization for fpopt 5
c
              id_out = 0 
              is_out = 0
              rsect  = 0
              zsect  = 0
c
              if (fpopt.eq.5.or.fpopt.eq.6) then 
c
c                For regular ions istate = iz 
c
                 istate = iz

c
c                call periphery transport routine
c
                 call fp_transport(imp,ik,ir,iz,istate,s,theta,
     >                        cross,vel,temi,
     >                        crmi,nrand,
     >                        cist,cistfp,cstmax,ctemav,rsect,zsect,
     >                        sputy,res)
c
c                If a wall collision has occured then refine the rsect,zsect 
c                values returned to actual wall locations and determine the ID
c                and IS values related to that location. 
c
                 if (res.eq.3) then 
                    call find_nearest_point_on_wall(rsect,zsect,
     >                                              id_out,is_out)
                 endif

              else
c
                 RES = FPERIPH(CIST,cistfp,FPDIST,fplosstim,CSTMAX,
     >                      NRAND,DIFFR,SEED,0.0)
c
              endif
c
c
c             Update ion lifetimes for far periphery excursion
c
              cist  = cist + cistfp
              cistiz= cistiz + cistfp  
c
              if (debug)  
     >           WRITE(6,'(a,i5,1p,6g12.5)') 'FP:',RES,CIST,CSTMAX,
     >                      qtim,fpdist,fplosstim,DIFFR
              if (debug0)  
     >           WRITE(0,'(a,i5,1p,6g12.5)') 'FP:',RES,CIST,CSTMAX,
     >                      qtim,fpdist,fplosstim,DIFFR
c
              IF (RES.EQ.1) THEN
c
c               Particle returns to grid. For FPOPT = 5, the 
c               cross field location of the particle after
c               stepping across the boundary is used to determine
c               the starting position of the particle. 
c
C
C               SET CROSS TO ONE STEP INSIDE THE BOUNDARY.
C               THIS IS THE POSITION IT WILL HAVE REACHED
C               IN ORDER TO EXIT THE FPERIPH ROUTINE
c
c               jdemod - for fpopt 5 ... proper transport
c                        across the FP boundary is needed to 
c                        ensure continuity of the calculated
c                        density profile so the code here and in
c                        fp_transport has been modified to make use
c                        of the actual cross position of the particle
c                        when moving in and out of the FP. 
c
C
                FPEXIT = FPEXIT + SPUTY
c
                IF (CIONR.EQ.1) THEN
                  IF (IR.EQ.IRWALL.or.ir.eq.irwall2) THEN
                    CROSS = KPERPS(IK,IR)
                  ELSEIF (IR.EQ.IRTRAP.or.ir.eq.irtrap2) THEN
                    CROSS = -KPERPS(IK,IR)
                  ENDIF
                ELSEIF (CIONR.EQ.0.or.cionr.eq.2) THEN
                  IF (IR.EQ.IRWALL.or.ir.eq.irwall2) THEN
                    JK = IKINS(IK,IR)
                    IR = IRINS(IK,IR)
c
                    if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
                      IF( IR.EQ.NRS
     >                    .AND. JK.GE.IKTI )THETA = THETA - DTHETG
c nonorth
                    endif
c
                    IK = JK
c
c                   jdemod - handle return from fpopt 5
c
                    if (fpopt.eq.5.or.fpopt.eq.6) then 
                       ! cross is always negative on return from fpopt 5
                       if (abs(cross).gt.abs(distout(ik,ir))) then 
                            write(0,'(a,2i6,10(1x,g12.5)') 
     >                             'FP RE_ENTRY WARNING:',
     >                             ik,ir,cross,distout(ik,ir)
                       endif
                       CROSS = -distout(ik,ir) - cross
                    else
                       CROSS = -distout(ik,ir) + KPERPS(IK,IR)
                    endif


 
                  ELSEIF (IR.EQ.IRTRAP.or.ir.eq.irtrap2) THEN
                    JK = IKOUTS(IK,IR)
                    IR = IROUTS(IK,IR)
c
                    if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
                      IF( IR.EQ.IRSEP
     >                    .AND. JK.GT.IKTO )THETA = THETA + DTHETG
c nonorth
                    endif
c
                    IK = JK
c
c                   jdemod - handle return from fpopt 5
c
                    if (fpopt.eq.5.or.fpopt.eq.6) then 
                       ! cross is always negative on return from fpopt 5
                       if (abs(cross).gt.abs(distin(ik,ir))) then 
                            write(0,'(a,2i6,10(1x,g12.5)') 
     >                             'FP RE_ENTRY WARNING:',
     >                             ik,ir,cross,distin(ik,ir)
                       endif
                       CROSS = distin(ik,ir) + cross
                    else
                       CROSS = distin(ik,ir) - KPERPS(IK,IR)
                    endif
                  else
                    write(6,'(a,2i6)')'CHECK_REACHED_GRID_EDGE:ERROR:'//
     >                     'PARTICLE NOT IN WALL RING AT FP EXIT:',ik,ir
                    write(0,'(a,2i6)')'CHECK_REACHED_GRID_EDGE:ERROR:'//
     >                     'PARTICLE NOT IN WALL RING AT FP EXIT:',ik,ir

                  ENDIF

c
c                 Find the revised theta value
c
                  call calculate_theta(ik,ir,s,theta)
c
c                 Update the value of SMAX
c
                  smax = ksmaxs(ir)
c
c
c                 I don't think this code is either needed or desired upon fp exit
c                 1) Either the particle returns exactly where it left
c                 2) It moves along the ring for fpopt 5 
c                 In either case I don't think this processing is appropriate.  
c     
c
c                  if (northopt.eq.1.or.northopt.eq.3) then
c nonorth
c
c
c                 IF (TAGDV(JK,IR).EQ.0
c     >                .AND. TAGDV(IKOLD,IROLD).EQ.0) THEN
C
C                   BOTH POINTS UNSHIFTED FROM ORTHOGONAL GRID
C
c                    IK = JK
c                  ELSE
C
C                   EITHER POINT SHIFTED FROM ORTHOGONAL GRID
C
c                    IF (IROLD.EQ.IRWALL) THEN
c                      JK = IKING(IK,IROLD)
c                    ELSEIF (IROLD.EQ.IRTRAP) THEN
c                      JK = IKOUTG(IK,IROLD)
c                    ENDIF
c  772               IF( THETA.LT.THETAG(JK,IR) )THEN
c                      IF( JK.GT.1 )THEN
c                        THETA1 = (THETAG(JK,IR) - THETA)
c     >                           /(THETAG(JK,IR) - THETAG(JK-1,IR))
c                        IF( THETA1.GT.0.5 )THEN
c                          JK = JK - 1
c                          GOTO 772
c                        ENDIF
c                      ENDIF
c                    ELSE
c                      IF( JK.LT.NKS(IR) )THEN
c                        THETA1 = (THETA - THETAG(JK,IR))
c     >                           /(THETAG(JK+1,IR) - THETAG(JK,IR))
c                        IF( THETA1.GE.0.5 )THEN
c                          JK = JK + 1
c                          GOTO 772
c                        ENDIF
c                      ENDIF
c                    ENDIF
c                    IK = JK
c                  ENDIF
c
c                  IF (IROLD.EQ.IRWALL) THEN
c                    CROSS = -distout(ik,ir) + KPERPS(IK,IR)
c                  ELSEIF (IROLD.EQ.IRTRAP) THEN
c                    CROSS = distin(ik,ir) - KPERPS(IK,IR)
c                  ENDIF
c nonorth
c                  endif
c

                ENDIF
c
              ELSEIF (RES.EQ.2) THEN
C
C               CIST HAS EXCEEDED CSTMAX - GOTO 780 TO
C               FINISH PROCESSING
C
                FPTTOT = FPTTOT + SPUTY
                ifate = 3
                return
c
c                GOTO 780
c
              ELSEIF (RES.EQ.3) THEN
C
C               TREAT AS NORMAL WALL COLLISION
C
               if (fpropt.eq.0) then
c
c               Normal collision - no recycling
c
                CICABS(IZ) = CICABS(IZ) + SPUTY
c                CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c                CILABS(IZ) = MAX (CILABS(IZ), CIST)
                CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
                CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
                CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c               Recording addition to total elapsed time in state IZ
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
                RWALL      = RWALL + SPUTY
                WALLS(IK,IR,IZ) = WALLS(IK,IR,IZ) + SPUTY
c
                ENERGY = 3.0 * real(IZ) * wallpt(id_out,29) +
     >           5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c
c               Add ion weight to wall element closest to grid 
c               departure.
c     
c               FP option 5 provides the id value of the element as well 
c               as the approximate  R,Z value where the particle strikes
c
c
c                write(6,*) 'update_walldep: '//
c     >                   'check_reached_grid_edge - FP res=3'

                call update_walldep(ik,ir,iz,0,id_out,
     >                              iwstart,idtype,sputy,energy)
c
                IFATE = 1

! ammod begin.	       
	       ! WBC comparison addition for ion prompt deposition.
               call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
! ammod end.	       


                return
c 
c                GOTO 790
c
               elseif (fpropt.eq.1) then
c
c               Recycle particle from edge of nearest plate
c
c               For the outer wall the particle will
c               be launched from either ID=1 or ID=NDS
c               For the trap wall it will be launched
c               from ID = NDSIN or ID = NDSIN+1
c
c               Unless DDS(ID) = 0 in which case it shifts
c               accross to the first target segment with non-zero
c               size.
c
c               Note: To have reached this point IR must be equal
c               to IRWALL or IRTRAP
c
c
c               Record statistics as if for normal wall collision
c
c         WRITE(6,*) 'FP RELAUNCH BW:',IK,IR,R,Z,iz
c


c ---      REFLECT ---- IONS may need to be reflected here


                CICABS(IZ) = CICABS(IZ) + SPUTY
c                CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c                CILABS(IZ) = MAX (CILABS(IZ), CIST)
                CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
                CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
                CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c               Recording addition to total elapsed time in state IZ
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
                RWALL      = RWALL + SPUTY
                WALLS(IK,IR,IZ) = WALLS(IK,IR,IZ) + SPUTY
c
c               Add ion weight to wall element closest to grid 
c               departure.
c

c
c NOTE:!!! This code needs adjustments to support relaunch of fp particles from the walls
c          if fp recycle is turned on.This code is also set up to calculate a self-sputtering
c          event which may or may not be correct depending on the context.  
c
c


c
c               Find target segment for re-launch
c
                if (ik.gt.nks(ir)/2) then 
                   ik = nks(ir)
                   id = verify_id(ik,ir,1)
                else
                   ik = 1
                   id = verify_id(ik,ir,2)
                endif
c
c
c               Postion on target/initial position options
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
c               Do not record the statistics of this
c               as a standard relaunched particle.
c
c
C
C
c            CICABS(IZ) = CICABS(IZ) + SPUTY
c            CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c            CILABS(IZ) = MAX (CILABS(IZ), CIST)
c            CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
c            CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
c            CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c            CRTABS(IZ) = CRTABS(IZ) + TEMI * SPUTY
c            CRVABS(IZ) = CRVABS(IZ) + VEL * SPUTY
c            CRAVAV(IZ) = CRAVAV(IZ) + ABS(VEL) * SPUTY
c            CTBS  (IZ) = CTBS  (IZ) + KTEBS(IK,IR) * SPUTY
c            IM         = MIN (INT(TEMI/(0.2*CTEB0))+1, 10)
c            CTEXS(IM)  = CTEXS(IM) + TEMI * SPUTY
C
c            RDEP   = RDEP + SPUTY
c            DEPS(ID,IZ) = DEPS(ID,IZ) + SPUTY
c            NEROS(ID,1) = NEROS(ID,1) + SPUTY
c
c
                ENERGY = 3.0 * RIZ * KTEBS(IK,IR) +
     >            5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c
c
c                write(6,*) 'update_walldep: '//
c     >                   'check_reached_grid_edge - FP Recycle'

                call update_walldep(ik,ir,iz,0,id_out,
     >                              iwstart,idtype,sputy,energy)
c
c
                if (kmfss(id).ge.0.0) then  
                   RYIELD = YIELD (6, MATTAR, ENERGY,
     >                      ktebs(ik,ir),ktibs(ik,ir)) * KMFSS(ID)
                elseif (kmfss(id).lt.0.0.and.kmfss(id).ge.-50.0) then 
                   RYIELD = abs(KMFSS(ID))
                elseif (kmfss(id).le.-99.0) then 
                   RYIELD = YIELD (6, MATTAR, ENERGY,
     >                       ktebs(ik,ir),ktibs(ik,ir))
                endif
c
                SPUNEW = SPUTY * RYIELD
                YLDTOT = YLDTOT + SPUNEW
                YLDMAX = MAX (YLDMAX, SPUNEW)
C
                IF (SPUNEW.GT.CTRESH) THEN
                  YTHTOT = YTHTOT + SPUNEW
                  NPROD  = NPROD + 1
                  SNEWS(NPROD) = SPUNEW
                  IF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.
     >                CNEUTD.EQ.4) THEN
                    EMAX = CEMAXF * ENERGY
                    RMAXS(NPROD) = 1.0 / (1.0+CEBD/EMAX)**2
                  ELSE
                    RMAXS(NPROD) = 1.0
                  ENDIF
                  XPRODS(NPROD) = R
                  YPRODS(NPROD) = Z
c
c                 For segments with a fixed sputtering yield - allow for
c                 the energy of the sputtered particle to be set to a
c                 specific value. 
c
                  if(kmfss(id).lt.0.0.and.cselfs.eq.2) then
                     eprods(nprod) = ctem1
                  else
                     eprods(nprod) = 0.0
                  endif
c
                  IDPRODS(NPROD) = ID
                  launchdat(nprod,2) = 1.0
                ENDIF
c
c           WRITE(6,*) 'FP RELAUNCH:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID),
c     >                 spunew,nprod,energy
c 
                IFATE = 2


! ammod begin.	       
	       ! WBC comparison addition for ion prompt deposition.
               call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
! ammod end.	       

c
                return 
c
c                GOTO 790
c
               endif
c
              ELSEIF (RES.EQ.4) THEN
c
c               FP opt 5 - may need to add code to generate a better estimate of where the particle 
c               struck the target. 
c
C
C               FP TARGET IMPACT - RECORD AND GOTO NEXT ION
C               TREAT AS NORMAL TARGET LOSS FOR STATISTICS BUT
C               NOT SELF-SPUTTERING
C



               if (fpropt.eq.0) then
                CICABS(IZ) = CICABS(IZ) + SPUTY
c                CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c                CILABS(IZ) = MAX (CILABS(IZ), CIST)
                CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
                CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
                CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c               Recording addition to total elapsed time in state IZ
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
C
                FPTARG(IZ) = FPTARG(IZ) + SPUTY
                FPTART     = FPTART + SPUTY
                rfptarg    = rfptarg + sputy 
c
c               Find closest target segment for loss
c
                if (ik.gt.nks(ir)/2) then 
                   ik = nks(ir)
                   id = verify_id(ik,ir,1)
                else
                   ik = 1
                   id = verify_id(ik,ir,2)
                endif
c
c               Update SMAX for actual exit ring
c
                smax = ksmaxs(ir)
c
c               Add ion weight to wall element 
c               closest to grid departure.
c
c                write (6,*) 'FPTARG:',ik,ir,id,wallindex(id)
c
c                write(6,*) 'update_walldep: '//
c     >             'check_reached_grid_edge - FP res=4 '
c
                ENERGY = 3.0 * real(IZ) * KTEDS(ID) +
     >            5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI

c
                call update_walldep(ik,ir,iz,id, 0,
     >                              iwstart,idtype,sputy,energy)
c
                IFATE = 9
                return
c
c                GOTO 790
c
               elseif (fpropt.eq.1) then
c
c               Recycle particle from edge of nearest plate
c
c               Particles trigger a self-sputtering event upon recycling
c
c               For the outer fp target the particle will
c               be launched from either ID=1 or ID=NDS
c               For the trap fp target it will be launched
c               from ID = NDSIN or ID = NDSIN+1
c
c               Unless DDS(ID) = 0 in which case it shifts
c               accross to the first target segment with non-zero
c               size.
c
c               Note: To have reached this point IR must be equal
c               to IRWALL or IRTRAP
c
c               Record statistics as if for normal FP-TARGET collision
c
c             WRITE(6,*) 'FP RELAUNCH BT:',IK,IR,R,Z,iz
c 

c ---      REFLECT ---- IONS may need to be reflected here

                CICABS(IZ) = CICABS(IZ) + SPUTY
c                CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c                CILABS(IZ) = MAX (CILABS(IZ), CIST)
                CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
                CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
                CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c
c               Recording addition to total elapsed time in state IZ
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
                FPTARG(IZ) = FPTARG(IZ) + SPUTY
                FPTART     = FPTART + SPUTY
                rfptarg    = rfptarg + sputy
c
c
c               Find target segment for re-launch
c
                if (ik.gt.nks(ir)/2) then 
                   ik = nks(ir)
                   id = verify_id(ik,ir,1)
                else
                   ik = 1
                   id = verify_id(ik,ir,2)
                endif
c
c               Update SMAX for actual exit ring
c
                smax = ksmaxs(ir)
c
c
c               Postion on target/initial position options
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
c
c               Add ion weight to wall element closest to grid 
c               departure.
c

                ENERGY = 3.0 * real(IZ) * KTEDS(ID) +
     >            5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI


c                write(6,*) 'update_walldep: '//
c     >             'check_reached_grid_edge - FP res=4 Recycle'

                call update_walldep(ik,ir,iz,id, 0,
     >                              iwstart,idtype,sputy,energy)
c
C
C            WRITE(6,*) 'SPUTTERED:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID)
C
c            CICABS(IZ) = CICABS(IZ) + SPUTY
c            CIFABS(IZ) = MIN (CIFABS(IZ), CIST)
c            CILABS(IZ) = MAX (CILABS(IZ), CIST)
c            CIFABS(IZ) = MIN (CIFABS(IZ), sngl(CIST))
c            CILABS(IZ) = MAX (CILABS(IZ), sngl(CIST))
c            CISABS(IZ) = CISABS(IZ) + CIST * SPUTY
c            CRTABS(IZ) = CRTABS(IZ) + TEMI * SPUTY
c            CRVABS(IZ) = CRVABS(IZ) + VEL * SPUTY
c            CRAVAV(IZ) = CRAVAV(IZ) + ABS(VEL) * SPUTY
c            CTBS  (IZ) = CTBS  (IZ) + KTEBS(IK,IR) * SPUTY
c            IM         = MIN (INT(TEMI/(0.2*CTEB0))+1, 10)
c            CTEXS(IM)  = CTEXS(IM) + TEMI * SPUTY
C
c            RDEP   = RDEP + SPUTY
c            DEPS(ID,IZ) = DEPS(ID,IZ) + SPUTY
c            NEROS(ID,1) = NEROS(ID,1) + SPUTY
c
c
 
                ENERGY = 3.0 * RIZ * KTEBS(IK,IR) +
     >            5.22E-9 * CRMI * VEL/QTIM * VEL/QTIM + 2.0 * TEMI
c
                if (kmfss(id).ge.0.0) then  
                   RYIELD = YIELD (6, MATTAR, ENERGY,
     >                      ktebs(ik,ir),ktibs(ik,ir)) * KMFSS(ID)
                elseif (kmfss(id).lt.0.0.and.kmfss(id).ge.-50.0) then 
                   RYIELD = abs(KMFSS(ID))
                elseif (kmfss(id).le.-99.0) then 
                   RYIELD = YIELD (6, MATTAR, ENERGY,
     >                      ktebs(ik,ir),ktibs(ik,ir))
                endif
c
                SPUNEW = SPUTY * RYIELD
                YLDTOT = YLDTOT + SPUNEW
                YLDMAX = MAX (YLDMAX, SPUNEW)
C
                IF (SPUNEW.GT.CTRESH) THEN
                  YTHTOT = YTHTOT + SPUNEW
                  NPROD  = NPROD + 1
                  SNEWS(NPROD) = SPUNEW
                  IF (CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5.OR.
     >                CNEUTD.EQ.4) THEN
                    EMAX = CEMAXF * ENERGY
                    RMAXS(NPROD) = 1.0 / (1.0+CEBD/EMAX)**2
                  ELSE
                    RMAXS(NPROD) = 1.0
                  ENDIF
                  XPRODS(NPROD) = R
                  YPRODS(NPROD) = Z
c
c                 For segments with a fixed sputtering yield - allow for
c                 the energy of the sputtered particle to be set to a
c                 specific value. 
c
                  if(kmfss(id).lt.0.0.and.cselfs.eq.2) then
                     eprods(nprod) = ctem1
                  else
                     eprods(nprod) = 0.0
                  endif
c
                  IDPRODS(NPROD) = ID
                  launchdat(nprod,2) = 1.0
                ENDIF
c
c           WRITE(6,*) 'FP RELAUNCH:',IK,IR,ID,R,Z,IKDS(ID),IRDS(ID),
c     >                 spunew,nprod,energy
c 
                IFATE = 2

! ammod begin.	       
	       ! WBC comparison addition for ion prompt deposition.
               call global_hc_wbc_comp(iz,crmi,vel,temi,sputy)
! ammod end.	       


                return
c
c                GOTO 790
c
               endif
              ENDIF
C
C           ELSEIF (FPOPT.EQ.2) THEN
C
C             THIS OPTION ALLOWS THE ION TO DRIFT OUT AS FAR AS
C             IT CAN. IT WILL ALWAYS BE ASSOCIATED WITH THE
C             CHARACTERISTICS OF THE OUTERMOST RING.
C
            ENDIF
          ENDIF


          return




      end
