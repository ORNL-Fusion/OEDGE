c     -*Fortran*-
c
c
      subroutine do_crossfield_step(ik,ir,ikold,irold,kk,s,theta,cross,
     >                        oldtheta,oldcross,
     >                        adjust,dcross,ckkmin,smax,k,debug,
     >                        seed,nrand,neutim,cist,imp,debug_all,
     >                        ifate)
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_crand
        implicit none
c
      integer ik,ir,ikold,irold,kk,nrand,imp,ifate
      real s,theta,adjust,dcross(4),cross,ckkmin,smax,k
      real oldtheta,oldcross,neutim
      real*8 seed,cist
      logical debug,debug_all
c
c*************************************************************************
c
c     DO_CROSSFIELD_STEP: formerly UPDATE_CROSS
c                   This routine replaces the code that was originally
c                   in the DIV module to implement the cross-field
c                   transport. The base option uses the code that was 
c                   in DIV - the other options implement different 
c                   algorithms or methods for calculating the cross-field
c                   transport. The routine takes all of the particle 
c                   position related quantites from DIV as it's input. In
c                   addition some variables that are used to accumulate 
c                   statistics are also passed in order to keep the 
c                   existing code as intact as possible.          
c
c
c     David Elder, Jan 22, 1997
c
c*************************************************************************
c
c     include    'params'
c     include    'comtor'
c     include    'cgeom'
c     include    'crand' 
c
c     Local variables
c
      real kprob,theta1,pinchvel
      real kdrefin,kdrefout
      real kratin,kratout
      real sdperptmp,tmpran 
      real kperpstepin,kperpstepout
      integer jk,ikreftmp,in
      integer flag
c
c     Available to trap error conditions if desired 
c
      integer ierr
      logical vr_assigned
c
C      real find_vr
C      external find_vr
c
      real tmpout,tmpin,crossfrac,tmptheta
c      real oldtheta,oldcross
      real tmpcross,tmpcross2
      integer iktmp,irtmp,jktmp
c
      real exb_rad_drftvel
c
c
c

c
c     Initialization - set vr_assigned to off for routine
c
      vr_assigned = .false.

C
C-------- UPDATE CROSS FIELD DIFFUSION. 
c
c         Record position
c
          IKOLD = IK
          IROLD = IR
          oldcross = cross
          oldtheta = theta 
c
c         Perform Cross-field diffusive step. - various methods
c
c         Applying different methods depending on options. 
c
c
c        if (debug) write(6,*) '3d:',ik,ir,s,cross,theta,adjust
          

c
c         Calculate the probability of inward step and perform 
c         the appropriate cross-field step.
c
c---------------------------------------------------------------
c
c         Regular cross-field diffusion - cioptj = 0,1 ...
c
       if (cioptj.eq.0.or.cioptj.eq.1) then 
c
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.tdistin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.tdistout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.2) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.tdistin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.tdistout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI.or.distin(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI.or.distout(ik,ir).eq.0.0) then 
                 kprob = kins(ik,ir)
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
c              write(6,*) 'Kpout:',ik,ir,kprat2(ik,ir,2),distout(ik,ir)
c

c
            endif 
c
          endif
c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)

c
c
c
c         Update Cross
c
c
c          write (6,*) 'Kprob:',ik,ir,kprob,cdiffopt,cross
c
c
          KK = KK +1 

          if (pinchopt.eq.4.and.vr_assigned) then 
             CROSS = CROSS + PINCHVEL
          else
             CROSS = CROSS + SIGN (KPERPS(IK,IR), kprob-ranv(kk))
     >                      + PINCHVEL 
          endif
c
          tmpran = ranv(kk)
c
c---------------------------------------------------------------
c
c         Cioptj = 2 - reference line proportional Dperp
c
c         This option provides an alternative way of calculating 
c         "CROSS" - the cross-field displacement. 
c
c
       elseif (cioptj.eq.2) then
c
c         Calculate step direction probability
c
c
c         Calculate the ACTUAL size of the Dperp step that will
c         be taken.
c
            if (ir.lt.irsep) then 
c
               kdrefin = tdistin(ikrefcore,ir)
               kdrefout= tdistout(ikrefcore,ir)
               sdperptmp = sdperpref
c
            elseif (ir.ge.irsep.and.ir.le.irwall) then 
c
               kdrefin =  tdistin(ikrefsol,ir)
               kdrefout= tdistout(ikrefsol,ir)
               sdperptmp  = sdperpref
c
            elseif (ir.ge.irtrap.and.ir.le.nrs) then 
c
               kdrefin = tdistin(ikrefpp,ir)
               kdrefout= tdistout(ikrefpp,ir)
c
c              Keep the Dperp in the Private Plasma constant on
c              each set of knots matching the corresponding cell 
c              on the separatrix so that a pinch does not develop - 
c              This matches the procedure that appears to have been
c              used in EDGE2D.
c
               ikreftmp= ikouts(ik,nrs) 
c
               sdperptmp=sdperpref
     >               *(distin(ikreftmp,irsep)
     >               /(distin(ikrefsol,irsep)))
c
c               sdperptmp  = sdperppp
c
            endif
c
            if (kdrefout.eq.0.0) then 
               kratout = 1.0
            else   
               kratout = tdistout(ik,ir)/kdrefout
            endif
c
            if (kdrefin.eq.0.0) then 
               kratin = 1.0
            else   
               kratin = tdistin(ik,ir)/kdrefin
            endif
c
c           These are estimates necessary for approximating the step probability
c
c           The estimate process in the Private plasma keeps the Dperp constant
c           across a set of knots starting on the separatrix. This matches what
c           is done in EDGE2D and will avoid pinch effects but does not seem very 
c           physical.  
c   
            if (ir.ge.irtrap.and.ir.le.nrs) then 

               kperpstepin = sdperptmp 
               kperpstepout = sdperptmp 
               kratin = 1.0
               kratout = 1.0

            else  

               kperpstepin = sdperptmp * kratin
               kperpstepout = sdperptmp* kratout

            endif
c
c 
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistin(ik,ir)))
              endif
c
            endif 
c
         elseif (cdiffopt.eq.2) then 
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperpstepin/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperpstepout/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
            endif 
c
          endif
c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)

c
c  
c         Map to cross on reference line
c
          if (cross.gt.0.0) then 
c
             tmpcross = cross / kratin
c
          elseif (cross.le.0.0) then 
c
             tmpcross = cross / kratout
c
          endif
c
c
c         Update TmpCross
c
          kk = kk +1
          if (pinchopt.eq.4.and.vr_assigned) then 
             tmpCROSS = tmpCROSS + PINCHVEL
          else
             tmpCROSS = tmpCROSS + SIGN (sdperptmp,kprob-ranv(kk))
     >                      + PINCHVEL
          endif
c
c         Map cross back using the reference line ratios again
c
          if (cross.gt.0.0) then 
c
             cross = tmpcross * kratin
c
          elseif (cross.le.0.0) then 
c
             cross = tmpcross * kratout
c
          endif
c
c
c
c     Spatially varying Dperp - more complicated 
c
      elseif (cioptj.eq.3.or.cioptj.eq.4) then 
c
c         Uses Kperps(ik,ir) with spatially
c         varying values.  
c
c
          if (cdiffopt.eq.0.or.
     >        (cdiffopt.eq.1.and.ir.ge.irsep)) then 
c
            kprob = kins(ik,ir) 
c
c
          elseif (cdiffopt.eq.1.and.ir.lt.irsep) then 
c
c
            if ((ir.eq.irsep-1.and.cross.le.0.0).or.ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir) ))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.2) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/tdistin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/tdistout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/tdistout(ik,ir)))
              endif
c
            endif 
c
          elseif (cdiffopt.eq.3) then 
c
c
            if (ir.eq.1) then 
c
              kprob = kins(ik,ir)
c
            elseif (cross.gt.0.0) then 
c
              if (kprat2(ik,ir,1).eq.HI) then 
                 kprob = 0.5
              else
                 kprob = (kprat2(ik,ir,1)+
     >                ((cross+kperps(ik,ir)/2.0)/distin(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,1)+cross/distin(ik,ir)))
              endif 
c
            elseif (cross.le.0.0) then 
c
              if (kprat2(ik,ir,2).eq.HI) then 
                 kprob = 0.5
              else 
                 kprob = (kprat2(ik,ir,2)+
     >                ((cross+kperps(ik,ir)/2.0)/distout(ik,ir)))
     >                / (2.0 * (kprat2(ik,ir,2)+cross/distout(ik,ir)))
              endif
c
            endif 
c
          endif
c
          call set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)

c
c
c
c         This is where the process gets fancy - check to
c         see if the Dperp step will carry cross into the 
c         next cell - if it does - use a portion of the
c         step-size in the next cell to calculate the 
c         final position of the particle. I.E. Use the
c         differeing cross-field step sizes in the 
c         adjoining cells in proportion to amounts already
c         used.
c
c         e.g. for example if a dperp step carries a particle
c              across the boundary of cell A at a distance
c              that is 40% of a cell A step - then beyond the
c              boundary - the particle will be moved 60% of
c              a cell B step.  
c
c         The Dperp step is done first and THEN a pinch velocity
c         (if any) is added on. 
c          
          oldcross = cross
c
          kk = kk +1
          tmpCROSS = CROSS + SIGN (KPERPS(IK,IR), kprob-RANV(KK))
c
c
c         Check to see if has crossed a cell boundary. 
c
c         Check for only ONE boundary for now - if the particle
c         is crossing multiple rings with one cross-field step
c         then the time-step is too large. However, there is 
c         code that checks this condition and prints an error
c         message. 
c
c
c         First - find out exactly what cell the particle will
c         be in - so we can use the appropriate Dperp for the 
c         second portion of the step.
c          
c
          tmpcross2 = tmpcross
          tmptheta = theta
          iktmp = ik
          irtmp = ir
          jktmp = jk
c
c         Use the DO_CFSTEP routine to find the new cell into which the 
c         particle will be transported with it's current value of cross. 
c
          call do_cfstep(jktmp,iktmp,irtmp,irold,tmpcross2,adjust,
     >                   tmptheta,flag,debug_all)
c
          if (flag.gt.0) then 
c
c            Particle moved inward - use the Dperp of the cell 
c            it would have wound up in to calculate the final 
c            value of cross. Can double check this if one wants
c            to by calling do-cfstep again with the corrected
c            value of cross. 
c      
c            Note: cosali and cosalo have been set to 1.0 for all
c            cells on grids that have been assumed orthogonal. 
c
             crossfrac = (distin(ik,ir)- oldcross) / kperps(ik,ir)
c
             tmpcross =  distin(ik,ir)
     >                  + (1.0-crossfrac) * kperps(iktmp,irtmp)
c 
             if (flag.gt.1) then 
                write (6,*) 'ERROR: Cross-field step error -' 
                write (6,*) '       Crossed multiple rings = ',flag
                write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
             end if 
c
c            Info for debugging purposes ...
c             
             tmpout = distout(iktmp,irtmp)
             tmpin = distin(ik,ir)
c
c
c
c             write (6,'(a,5g13.6)') 'DperpIN:',tmpcross,tmpin,
c     >                 tmpcross2,
c     >                 tmpcross-tdistin(ik,ir),
c     >                 tmpout
c             write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
c
          elseif (flag.lt.0) then  
c      
c            Note: cosali and cosalo have been set to 1.0 for all
c            cells on grids that have been assumed orthogonal. 
c
             crossfrac = abs(-distout(ik,ir)-
     >                     oldcross) / kperps(ik,ir)
c
             tmpcross = -distout(ik,ir)
     >                  - (1.0-crossfrac) * kperps(iktmp,irtmp)
c 
             if (flag.lt.-1) then 
                write (6,*) 'ERROR: Cross-field step error -' 
                write (6,*) '       Crossed multiple rings = ',flag
                write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
             end if 
c
c            Info for debugging purposes ...
c             
             tmpout = distout(iktmp,irtmp)
             tmpin =  distin(ik,ir)
c
c             write (6,'(a,5g13.6)') 'DperpIN:',tmpcross,tmpin,
c     >                 tmpcross2,
c     >                 tmpcross-tdistin(ik,ir),
c     >                 tmpout
c             write (6,*) 'IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp
c
          endif
c
          if (flag.eq.0.and.(ir.ne.irtmp.or.ik.ne.iktmp)) then 
             
             write (6,*) 'FLAG=0:IK,IR,IKTMP,IRTMP:',ik,ir,iktmp,irtmp

          endif    
c
c         Finalize the updating of cross by adding the pinch velocity
c
          if (pinchopt.eq.4.and.vr_assigned) then 
             CROSS = cross + PINCHVEL
          else
             CROSS = tmpcross + PINCHVEL
          endif
c
c
c         End of cioptj IF statement
c
       endif
C
c      jdemod - all cross updates continue here - add in the radial exb drift term at this point
c
c
       call set_exb_rad_drift(ik,ir,exb_rad_drftvel)
c
c     update cross with the exb term
c
       cross = cross + exb_rad_drftvel
c


          IF (DEBUG_ALL) WRITE (6,1000) 'D2:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,kprob,
     >      kprob-tmpran,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'




      if (debug_all) write(6,'(a,i4,3(1x,g12.5))') 
     >             'DO  CROSS :',ifate,s,theta,cross


c
c     Update S and other values based on new cross-field location
c
c

      call update_crossfield(ik,ir,ikold,irold,kk,s,theta,cross,
     >                        oldtheta,oldcross,
     >                        adjust,dcross,ckkmin,smax,k,debug,
     >                        seed,nrand,neutim,cist,imp,debug_all,
     >                        ifate)
c
c      if (ifate.ne.0) return 
c


c
c     Format statements
c
 1000 format(a,2i4,1p,10(g11.4),1x,a) 
c 1001 format(a,2i4,1p,6(g11.4),44x,1x,a) 


        return
        end



      subroutine update_crossfield(ik,ir,ikold,irold,kk,s,theta,cross,
     >                        oldtheta,oldcross,
     >                        adjust,dcross,ckkmin,smax,k,debug,
     >                        seed,nrand,neutim,cist,imp,debug_all,
     >                        ifate)
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_crand
      implicit none
c
      integer ik,ir,ikold,irold,kk,nrand,imp,ifate
      real s,theta,adjust,dcross(4),cross,ckkmin,smax,k
      real oldtheta,oldcross
      real*8 seed,cist
      logical debug,debug_all
      real neutim
c   
c     include    'params'
c     include    'comtor'
c     include    'cgeom'
c     include    'crand' 
c
c     Local variables 
c
      real theta1
      integer flag,jk
c

c
c         Move the particle across rings if the CROSS value
c         has exceeded the distance to the next ring. Recalculate
c         THETA if using NON-orthogonal transport
C
          call do_cfstep(jk,ik,ir,irold,cross,adjust,theta,
     >                   flag,debug_all)
c
c         Record some statistics in the core. 
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
          IF (DEBUG_ALL) WRITE (6,1000) 'D3:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,0.0,
     >      0.0,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'





c
c         Re-calculate S if the particle has changed rings.  
c
c         There is a bug caused by adjusting S when entering 
c         a virtual ring. The FP code does not properly deal with
c         it - in order to fix this - try not adjusting S for 
c         values in the virtual rings.
c
c
          IF (IR.NE.IROLD.and.
     >       (ir.ne.irwall.and.ir.ne.irtrap)) THEN
            K      = KKS(IR)
            CKKMIN = MIN (CKKMIN, K)
            SMAX   = KSMAXS(IR)
c 
c           ITER grid 
c
            if (cgridopt.eq.2) then
               if ( (((ir.ge.irsep.and.ir.le.irwall2).or.
     >               (ir.ge.irsep2.and.ir.le.irwall)).and.
     >               (irold.ge.irtrap)).or.
     >              (((irold.ge.irsep.and.irold.le.irwall2).or.
     >               (irold.ge.irsep2.and.irold.le.irwall)).and.
     >               (ir.ge.irtrap))) then
                  IF (S.GT.KSS(IKOLD,IROLD)) THEN
                    S = KSS(IK,IR) - (S-KSS(IKOLD,IROLD)) *
     >                         (KBACDS(IK,IR)/KFORDS(IKOLD,IROLD))
                  ELSE
                    S = KSS(IK,IR) + (KSS(IKOLD,IROLD)-S) *
     >                         (KFORDS(IK,IR)/KBACDS(IKOLD,IROLD))
                  ENDIF
               elseIF (S.GT.KSS(IKOLD,IROLD)) THEN
                    S = KSS(IK,IR) + (S-KSS(IKOLD,IROLD)) *
     >                         (KFORDS(IK,IR)/KFORDS(IKOLD,IROLD))
               ELSE
                  S = KSS(IK,IR) - (KSS(IKOLD,IROLD)-S) *
     >                         (KBACDS(IK,IR)/KBACDS(IKOLD,IROLD))
               ENDIF
c
c           Orthogonal Transport 
c
            ELSEIF (NORTHOPT.EQ.0.or.northopt.eq.2) THEN
c
              IF (S.GT.KSS(IKOLD,IROLD)) THEN
                  S = KSS(IK,IR) + (S-KSS(IKOLD,IROLD)) *
     >                         (KFORDS(IK,IR)/KFORDS(IKOLD,IROLD))
              ELSE
                 S = KSS(IK,IR) - (KSS(IKOLD,IROLD)-S) *
     >                         (KBACDS(IK,IR)/KBACDS(IKOLD,IROLD))
              ENDIF
c
c
c           Non-orthogonal Transport
c
c           This block generates a new value for S after non-orthogonal 
c           cross-field diffusion (that results in changing rings).
c
            ELSEIF (NORTHOPT.EQ.1.or.northopt.eq.3) THEN
c
c             Adjustment to Theta moved to DO_CFSTEP where
c             code decides new IK value of cell.
c
c             Special for the last knot on a core ring:
c
c              IF (IKOLD.EQ.NKS(IROLD).AND.IR.GE.IRSEP
c     >            .and.irold.lt.irsep) THEN
c                THETA = THETAG(IK,IR) - 
c     >                  (THETAG(IK,IR) - THETAG(IK-1,IR)            ) *
c     >                  (THETAG(IKOLD,IROLD) - THETA                ) /
c     >                  (THETAG(IKOLD,IROLD) - THETAG(IKOLD-1,IROLD))
c              ENDIF
c
              IF (IR.LT.IRSEP.AND.IK.EQ.1.AND.
     +            THETA.LT.THETAG(IK,IR)) THEN

                IK = NKS(IR)
c
c               Handle an error condition when IKOLD was also 1.
c               The spacing from nks(ir) to nks(ir) -1 on core rings
c               is the same as cell 1 to a mythical cell 0. since cell
c               1 and nks(ir) coincide. IKOLD should not be 1 for a 
c               non-core 
c
c
                if (ikold.eq.1) then 
                   THETA = THETAG(IK,IR) - 
     +                  (THETAG(IK,IR)       - THETAG(IK-1,IR)      ) *
     +                  (THETAG(IKOLD,IROLD) - THETA                ) /
     +          (THETAG(nks(irold),IROLD) - THETAG(nks(irold)-1,IROLD))
c
c               Regular case
c
                else 
                   THETA = THETAG(IK,IR) - 
     +                  (THETAG(IK,IR)       - THETAG(IK-1,IR)      ) *
     +                  (THETAG(IKOLD,IROLD) - THETA                ) /
     +                  (THETAG(IKOLD,IROLD) - THETAG(IKOLD-1,IROLD))
                endif               


          IF (DEBUG_ALL) WRITE (6,1000) 'D4:',IKold,IRold,S,K,
     >      THETA,thetag(ik,ir),thetag(ik-1,ir),
     >      thetag(ikold,irold),thetag(ikold-1,irold),
     >      0.0,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'



              ENDIF


          IF (DEBUG_ALL) WRITE (6,1000) 'D5:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,0.0,
     >      0.0,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'

c
c             Re-calculate THETA and S.
c
c
c             First Half of cell
c
              IF (THETA.LT.THETAG(IK,IR)) THEN
c
                IF (IK.EQ.1) THEN
                  IF (IR.LT.IRSEP) THEN
                    IK = NKS(IR)             
 
                    THETA = THETA + (THETAG(IK,IR) - THETAG(1,IR))

                    THETA1 = (THETAG(IK,IR) - THETA          ) /
     >                       (THETAG(IK,IR) - THETAG(IK-1,IR))
                    S = KSS(IK,IR) - KBACDS(IK,IR) * THETA1
                  ELSE
                    IF( THETA.LE.THETAT(IDDS(IR,2)) )THEN
c
c                     Particle has struck target cross-field 
c
c                     If Mirror target option is ON - place particle
c                     back in old ring with cross set to it's
c                     previous value.
c
c                     Otherwise do the usual. 
c
                      if (cmiropt.eq.0.or.cmiropt.eq.3) then 

                         S = 0.0

                      elseif (cmiropt.eq.1.or.
     >                        cmiropt.eq.2.or.
     >                        cmiropt.eq.4) then                      
c
c                        Do not change S - reset CROSS and IR
c
                         cross = oldcross
                         theta = oldtheta  
                         ir    = irold
                         ik    = ikold
c slmod begin
c...                     BUG: Need to restore this as well -SL 24.10.06
                         smax  = ksmaxs(ir)
c slmod end
c
                      endif
c
                    ELSE
                      S = KSS(IK,IR) * 
     >                    (THETA         - THETAT(IDDS(IR,2))) /
     >                    (THETAG(IK,IR) - THETAT(IDDS(IR,2)))
                    ENDIF
                  ENDIF
                ELSE
                  THETA1 = (THETAG(IK,IR) - THETA          ) /
     >                     (THETAG(IK,IR) - THETAG(IK-1,IR))
                  S = KSS(IK,IR) - KBACDS(IK,IR) * THETA1
                ENDIF
c
c             Particle in second half of cell. 
c
              ELSE
c
                IF( IK.EQ.NKS(IR) )THEN
                  IF (IR.LT.IRSEP) THEN
                    IK = 1
                    THETA = THETA - (THETAG(NKS(IR),IR) - THETAG(1,IR))
  
                    THETA1 = (THETA           - THETAG(IK,IR)) /
     >                       (THETAG(IK+1,IR) - THETAG(IK,IR))
                    S = KSS(IK,IR) + KFORDS(IK,IR) * THETA1
                  ELSE
                    IF (THETA.GE.THETAT(IDDS(IR,1))) THEN

c
c                     Particle has struck target cross-field 
c
c                     If Mirror target option is ON - place particle
c                     back in old ring with cross set to it's
c                     previous value.
c
c                     Otherwise do the usual. 
c
                      if (cmiropt.eq.0.or.cmiropt.eq.4) then 

                         S = SMAX

                      elseif (cmiropt.eq.1.or.
     >                        cmiropt.eq.2.or.
     >                        cmiropt.eq.3) then                      
c
c                        Do not change S - reset CROSS and IR
c
                         cross = oldcross 
                         theta = oldtheta
                         ir    = irold
                         ik    = ikold 
c slmod begin
c...                     BUG: Need to restore this as well -SL 24.10.06
                         smax  = ksmaxs(ir)
c slmod end
c
                      endif
c
                    ELSE
                      S = KSS(IK,IR) + KFORDS(IK,IR) *
     >                    (THETAT(IDDS(IR,1)) - THETA        ) /
     >                    (THETAT(IDDS(IR,1)) - THETAG(IK,IR))
                    ENDIF
                  ENDIF
                ELSE
                  THETA1 = (THETA           - THETAG(IK,IR)) /
     >                     (THETAG(IK+1,IR) - THETAG(IK,IR))
                  S = KSS(IK,IR) + KFORDS(IK,IR) * THETA1
                ENDIF
c
              ENDIF
c    
            ENDIF


          IF (DEBUG_ALL) WRITE (6,1000) 'D6:',IK,IR,S,K,
     >      THETA,SMAX,CROSS,adjust,0.0,
     >      0.0,
     >      distin(ik,ir),
     >      -distout(ik,ir),
     >      'UPDATED CROSS'



          ENDIF

c        if (debug) write(6,*) '3d:',ik,ir,cross,s

c
c     Adjust particle S value if required for the new ring 
c
      if (ir.lt.irsep) then 


C
C-------- LOOPING ROUND MAIN PLASMA CONTOURS
C
          IF (S.LT.0.0) THEN
  600       S = S + SMAX
            IF (S.LT.0.0) GOTO 600
          ELSEIF (S.GT.SMAX) THEN
  610       S = S - SMAX
            IF (S.GT.SMAX) GOTO 610
          ENDIF


c
c     Particle not in core 
c

      else
c
c        Routine may set IFATE for target impact
c
c        NOTE!: This code works properly in DIVIMP because the 
c               arguments passed in are actually stored in the 
c               common block "particle_specs" - as a result when
c               this subroutine is called - all the relevant 
c               data values have been updated. This is not 
c               immediately obvious. The intention is to 
c               eventually remove the dependence on the common
c               block or rewrite all the transport code as a module.
c               At the moment things are in a mixed state which 
c               makes it difficult to modify this code for the HC
c               routines. 
c

         call check_target_impact(seed,nrand,neutim)

      endif

      if (debug_all)  write(6,'(a,i4,3(1x,g12.5))') 
     >           'UPD CROSS :',ifate,s,theta,cross



c
c     End of routine 
c
      return
c
c     Format statements
c
 1000 format(a,2i4,1p,10(g11.4),1x,a) 
c 1001 format(a,2i4,1p,6(g11.4),44x,1x,a) 








      
      return
      end
c
c
c

      subroutine set_exb_rad_drift(ik,ir,exb_rad_drftvel)
      use mod_params
      use mod_driftvel
      implicit none
      integer :: ik,ir
      real :: exb_rad_drftvel
c     include 'params'
c     include 'driftvel'

c     assign exb radial drift if the option is ON
      if (exb_rad_opt .eq. 1) then 
         ! note the exb_rad_drft is defined so that the drift is "+" for radially outward
         ! similar to E-rad which is "+" for radially outward. 
         ! However, the code uses "+" for radially inward transport which requires the 
         ! sign of the radial drift velocity be changed for particle transport
         exb_rad_drftvel = -exb_rad_drft(ik,ir)
      else
         exb_rad_drftvel = 0.0
      endif

      return
      end
c
c
c
      subroutine do_cfstep(jk,ik,ir,irold,cross,adjust,theta,flag,debug)
      use mod_fp_data
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_fperiph_com
      implicit none
      integer jk,ik,ir,flag,irold
      real cross,theta,adjust
      logical debug
c
c     DO_CFSTEP: This routine executes the cross-field step and determines
c                which ik,ir cell the particle will end up in - this 
c                depends on the transport options chosen among other things. 
c
c                The flag value indicates how many rings were crossed with the
c                current Dperp step - if none then flag = 0 .
c
c                David Elder              1997, Feb 17 
c
c     include    'params'
c     include    'comtor'
c     include    'cgeom'
c     include    'fperiph_com'
c
c     Local variable,ik
c
      integer ikold,iroldtmp
      real theta1
c
c     Set flag to zero 
c
c     Flag is incremented by 1 for inward steps and decremented by 1 for 
c     outward steps.
c
      flag = 0 
      ikold = ik
      iroldtmp = ir

c        if (debug) write(6,*) '3ba:',ik,ir,cross,flag

c
c         Move the particle across rings if the CROSS value
c         has exceeded the distance to the next ring. Recalculate
c         THETA if using NON-orthogonal transport
C
c
c         Regular Transport
c
          if (northopt.eq.0.or.northopt.eq.2) then
c
c           Move IN 
c   
 640        IF (CROSS.GT.distin(ik,ir)) THEN
              CROSS = CROSS - tdistin(ik,ir)
              JK    = IKINS(IK,IR)
              IR    = IRINS(IK,IR)
              IK = JK
              flag = flag + 1
c
              if (distin(ik,ir).eq.0.0.and.distout(ik,ir).eq.0.0) then 
c
c                 write (6,'(a,2i4,1p,4g12.5)')
c     >                   'ERROR: Particle INWARD transport into'//
c     >                  ' zero volume ring.',ik,ir,cross,distin(ik,ir),
c     >                    distout(ik,ir)
c
c                jdemod - this condition is an indicator that the particle has
c                         stepped off the grid edge onto the boundary ring and potentially
c                         into the periphery. Depending on the options in play it is necessary
c                         to keep the actual cross location of the particle for 
c                         initialization of the periphery particle. 
c                         fp_cross_tmp should always be positive for a position in the FP
c                         given the current fp transport sign conventions
c
                 fp_cross_tmp = cross
c
                 cross = 0.0     
c
              else 
                 GOTO 640
              endif 
            ENDIF
c
c           Move OUT
c
 650        IF (CROSS.LT.-distout(ik,ir)) THEN
              CROSS = CROSS + tdistout(ik,ir)
              JK    = IKOUTS(IK,IR)
              IR    = IROUTS(IK,IR)
              IK    = JK
              flag  = flag -1
c
C
              if (distin(ik,ir).eq.0.0.and.distout(ik,ir).eq.0.0) then 
c
c                 write (6,'(a,2i4,1p,4g12.5)')
c     >                   'ERROR: Particle OUTWARD transport into'//
c     >                  ' zero volume ring.',ik,ir,cross,distin(ik,ir),
c     >                    distout(ik,ir)
c
c                jdemod - record cross value for use in FP
c
                 fp_cross_tmp = cross
c
                 cross = 0.0     
c
              else 
                 GOTO 650
              endif
            ENDIF
c
c         Non-orthogonal Transport
c
          elseif (northopt.eq.1.or.northopt.eq.3) then
c
c           Move IN
c
c
  641       IF (CROSS.GT.distin(ik,ir)) THEN
c
c
c              Set cross to the corrected value AFTER the new cell is 
c              known
c
c              CROSS = CROSS - tdistin(ik,ir)
c              JK    = IKING(IK,IR)
c
              JK    = IKINs(IK,IR)
              ir    = irins(ik,ir)
              flag = flag + 1
c
C             (IF SEPARATRIX HAS BEEN CROSSED IN INNER DIVERTOR REGION)
c
c             IKTI is the last cell, counting up from the INNER target,
c             on the separatrix, which is adjacent to a private plasma cell.
c
              IF( IR.gt.irtrap .AND. IK.GE.IKTI
     >            .and.iroldtmp.lt.irwall) THETA = THETA - DTHETG
c
c             Recalculate THETA
c
  645         IF( THETA.LT.THETAG(JK,IR) )THEN
                IF( JK.GT.1 )THEN
                  THETA1 = (THETAG(JK,IR) - THETA)
     >                   /(THETAG(JK,IR) - THETAG(JK-1,IR))
c
                  IF ((PDOPT.EQ.0.AND.THETA1.GT.0.5).OR.
     +                (PDOPT.EQ.1.AND.
     +                 KSS(JK,IR) - THETA1 * KBACDS(JK,IR).LT.
     +                 KSB(JK-1,IR))) THEN
c
                    JK = JK - 1
                    GOTO 645
                  ENDIF
                ENDIF
              ELSEIF( THETA.ge.THETAG(JK,IR) )THEN

                IF( JK.LT.NKS(IR) )THEN
                  THETA1 = (THETA - THETAG(JK,IR))
     >                   /(THETAG(JK+1,IR) - THETAG(JK,IR))
c
                  IF ((PDOPT.EQ.0.AND.THETA1.GE.0.5).OR.
     +                (PDOPT.EQ.1.AND.
     +               KSS(JK,IR) + THETA1 * KFORDS(JK,IR).GT.
     +               KSB(JK,IR))) THEN
                    JK = JK + 1
                    GOTO 645
                  ENDIF
                ENDIF
              ENDIF
              IK = JK

c
         if (debug) write(6,'(a,5i5,6g12.5)') 'DO_CF:1',ik,ir,
     >                  ikold,iroldtmp,
     >                  flag,cross,
     >                 cross-distin(ikold,iroldtmp)-distout(ik,ir),
     >                  distout(ik,ir),
     >                  distin(ikold,iroldtmp),theta
c
c             Set revised value of cross - based on cell moved into
c
              cross = cross - distin(ikold,iroldtmp)-distout(ik,ir)

c
c             Adjust cross if necessary for changes in IK
c
c              call adjust_cross(cross,adjust,ik,ir,ikold,
c     >                          iroldtmp,debug) 
              ikold = ik
              iroldtmp = ir
C
              if (distin(ik,ir).eq.0.0.and.distout(ik,ir).eq.0.0) then 
c
c                 write (6,'(a,2i4,1p,4g12.5)')
c     >                   'ERROR: Particle INWARD transport into'//
c     >                  ' zero volume ring.',ik,ir,cross,distin(ik,ir),
c     >                    distout(ik,ir)
c
c
c                jdemod - record cross value for use in FP
c
                 fp_cross_tmp = cross
c
                 cross = 0.0     
c
              else   
c
                 GOTO 641
c
              endif
            ENDIF
c
c           Move OUT
c
c
 651        IF (CROSS.LT.-distout(ik,ir))
     >                           THEN
c
c
c              Set cross to the corrected value AFTER the new cell is 
c              known
c
c
c              CROSS = CROSS + tdistout(ik,ir)
c              JK    = IKOUTG(IK,IR)
c
              JK    = IKOUTs(IK,IR)
              ir    = irouts(ik,ir)
              flag = flag - 1
c
C             (IF SEPARATRIX HAS BEEN CROSSED IN INNER DIVERTOR REGION)
c
c             IKTO is the first cell, counting up from the INNER target,
c             on the separatrix, after looping around the core plasma, 
c             which is adjacent to a private plasma cell.
c
              IF( IR.lt.IRwall .AND. IK.GT.IKTO
     >            .and.iroldtmp.gt.irtrap)   THETA = THETA + DTHETG
c
c             Re-calculate THETA
c
c             Modify THETA value if stepping out of the last cell on a core
c             ring into the main SOL so that the THETA value is appropriate
c             for that side of the SOL. Without this - the particle can
c             step out of the core into a SOL cell on the other side of the 
c             plasma adjacent to the PP. Note - TAU also revised so
c             that thetag(nks(ir),ir) for ir < irsep is not equal to 
c             thetag(1,ir)
c
              if (theta.le.thetag(ikold,iroldtmp).and.
     >            iroldtmp.eq.irsep-1.and.ikold.eq.nks(iroldtmp)) then   
c
                  theta = theta  
     >                  - (thetag(ikold,iroldtmp)-thetag(1,iroldtmp))
c
              endif
c

  655         IF( THETA.LT.THETAG(JK,IR) )THEN

                IF( JK.GT.1 )THEN

                  THETA1 = (THETAG(JK,IR) - THETA)
     >                     /(THETAG(JK,IR) - THETAG(JK-1,IR))
c
                   IF ((PDOPT.EQ.0.AND.THETA1.GT.0.5).OR.
     +                 (PDOPT.EQ.1.AND.
     +                  KSS(JK,IR)-THETA1*KBACDS(JK,IR).LT.KSB(JK-1,IR)
     +                  )) THEN
                    JK = JK - 1
                    goto 655
                  ENDIF
                ENDIF
              ELSEif ( THETA.ge.THETAG(JK,IR) ) THEN
                IF( JK.LT.NKS(IR) )THEN
                  THETA1 = (THETA - THETAG(JK,IR))
     >                   /(THETAG(JK+1,IR) - THETAG(JK,IR))
c
                   IF ((PDOPT.EQ.0.AND.THETA1.GE.0.5).OR.
     +                (PDOPT.EQ.1.AND.KSS(JK,IR)+THETA1*KFORDS(JK,IR)
     +                 .GT.KSb(JK,IR))) THEN
                    JK = JK + 1
                    GOTO 655
                  ENDIF
                ENDIF
              ENDIF
              IK    = JK
c
         if (debug) write(6,'(a,5i5,6g12.5)') 'DO_CF:2',ik,ir,
     >                  ikold,iroldtmp,
     >                  flag,cross,
     >                 cross + distout(ikold,iroldtmp)+distin(ik,ir),
     >                  distin(ik,ir),
     >                  distout(ikold,iroldtmp),theta


c
c             Set revised value of cross - based on cell moved into
c
              cross = cross + distout(ikold,iroldtmp)+distin(ik,ir)
c
c             Adjust cross if necessary for changes in IK
c
c              call adjust_cross(cross,adjust,ik,ir,ikold,
c     >                          iroldtmp,debug) 
              ikold = ik
              iroldtmp = ir
c
              if (distin(ik,ir).eq.0.0.and.distout(ik,ir).eq.0.0) then 
c
c                 write (6,'(a,2i4,1p,4g12.5)')
c     >                   'ERROR: Particle OUTWARD transport into'//
c     >                  ' zero volume ring.',ik,ir,cross,distin(ik,ir),
c     >                    distout(ik,ir)
c
c
c                jdemod - record cross value for use in FP
c
                 fp_cross_tmp = cross
c
                 cross = 0.0     
c
              else   
c
                 GOTO 651
c
              endif 
c
            ENDIF
  


          endif
c
c          if (debug) 
c     >        write (6,*) 'Return from CFSTEP:' 
c
      return
      end 
c
c
c
      subroutine adjust_cross(cross,adjust,ik,ir,ikold,irold,debug)
      use mod_params
      use mod_comtor
      use mod_cgeom
      use mod_crand
      implicit none
      real cross,adjust
      integer ik,ir,ikold,irold
      logical debug
c
c     ADJUST_CROSS: The purpose of this routine is to adjust the value
c                   of cross whenever the IK index of the particle changes.
c                   This can be due to either parallel or perpendicular 
c                   transport - as such - this routine is called after 
c                   both stages of the motion - one at the beginning of 
c                   cfield and the other at the end of do_cfstep.
c
c                   1997, July 4 - David Elder
c
c     include    'params'
c     include    'comtor'
c     include    'cgeom'
c     include    'crand' 
c
c     Local variables
c
c
c     Zero adjustment  
c
      adjust = 0.0

c
c     If in new cell - adjust cross-field term
c     
      IF (IK.NE.IKOLD.AND.CROSS.NE.0.0) THEN
c
c
c     
        IF (CROSS.GT.0.0) THEN
c     
           if (distin(ikold,ir).ne.0.0) then  
              ADJUST = CROSS * (distin(ik,ir)
     >                 / distin(ikold,ir)-1.0)
           else
              adjust = 0.0
           endif
c
        ELSEif (cross.lt.0.0) then 
c     
           if (distout(ikold,ir).ne.0.0) then  
      
              ADJUST = CROSS * (distout(ik,ir)
     >                   /distout(ikold,ir) -1.0)
           else
              adjust = 0.0
           endif
c     
        ELSE
c     
           adjust = 0.0
c     
        ENDIF
c
      
        CROSS = CROSS + ADJUST

         if (debug) write(6,'(a,2i5,5g12.5)') 'DO_ADJ:',ik,ir,
     >                  cross,adjust,cross-adjust
c
      
      ENDIF

c
c     End of cross field adjustment routine
c
      return
      end   
c
c
c
      real function find_vr(ik,ir,nrand,vr_assigned,
     >                      cist,imp,ierr)
      use mod_params
      use mod_cgeom
      use mod_comtor
      implicit none
      integer ik,ir,imp 
      logical vr_assigned
      real*8 cist
      integer nrand,ierr
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c
c     FIND_VR: 
c
c     NOTE: DIVIMP SIGN CONVENTION!!!
c
c     In DIVIMP a positive velocity is toward the core plasma
c     except in the PFZ - general usage has a positive velocity
c     outward - this code switches the sign on the applied
c     velocity appropriately depending on region. 
c
c
c
c     This routine randomly selects a VR/PINCH value from 
c     a specified probability distribution function for this
c     quantity. The Function VR_PDF will return a "probability"
c     weight value for any given input value of V_IN. By integrating
c     over the PDF and normalizing the result - it is then possible
c     to select a random number on [0,1] and then find a 
c     corresponding velocity. The same procedure is used in the 
c     find_thomson_velocity routines. 
c
c     CIST contains the total number of ion time steps spent 
c     following the current particle. 
c
c     IMP is the index of the current particle being followed - when 
c     IMP changes the time window is reset. 
c
c     Local declarations
c
      real       vr_pdf_random,result_val,vr_direction,ran,getranf
      external   vr_pdf_random,getranf
      integer    in
c
c     Set up local saved data
c
      !data       current_particle /0/    
      integer,save :: current_particle = 0
      real*8,save ::  last_time_chosen = 0.0,
     >                vr_last_assigned=0.0, 
     >                current_time=0.0
c
c     Initialization
c
c     Set VR_assigned equal to false to indicate that a Radial Velocity has NOT
c     been assigned by this routine. This will be set to true if a radial velocity 
c     is later assigned. This is necessary to allow DPERP transport and VR transport
c     to be exclusive when that option is selected. 
c
      vr_assigned = .false.
      find_vr = 0.0
c
c     Check to see if a value needs to be calculated based on the
c     section of the grid.
c
c     Option 2: Entire grid - no test necessary but sign change reguired
c                             for ir > irwall since PFZ direction 
c                             definition is reversed.
c
c     Option 0: Entire grid except PFZ
c
      if (pinch_loc_opt.eq.0.and.ir.gt.irwall) then 
         find_vr = 0.0
         return
c
c     Option 1: Only main SOL
c
      elseif (pinch_loc_opt.eq.1.and.
     >        (ir.gt.irwall.or.ir.lt.irsep)) then 
         find_vr = 0.0
         return
c
c     Option 3: Only main SOL above approximate Xpoint region.  
c     
      elseif (pinch_loc_opt.eq.3.and.(
     >        (ir.gt.irwall.or.ir.lt.irsep).or.
     >        ((ir.ge.irsep.and.ir.le.irwall).and.
     >        (ik.lt.(ikto+1).or.ik.gt.(ikti-1)))
     >       )) then 
         find_vr = 0.0
         return
c
c     Option 4: Main SOL above approximate Xpoint region + core
c     
      elseif (pinch_loc_opt.eq.4.and.(
     >        (ir.gt.irwall).or.
     >        ((ir.ge.irsep.and.ir.le.irwall).and.
     >        (ik.lt.(ikto+1).or.ik.gt.(ikti-1)))
     >       )) then 
         find_vr = 0.0
         return

      endif
c
c     VR will be assigned - set variable
c
      vr_assigned = .true.
c
c     In the code - a negative radial velocity is outward while 
c     a positive one is inward. (This sign convention reverses in 
c     the PFZ where a positive is outwards towards the wall.)
c
c     The input vr_pdf is entered using the convention that positive
c     velocities are radially outward. As a result, the vr_direction
c     variable is used to change the direction of the vr_pdf 
c     result to match the transport conventions in use on the grid. 
c
      if (ir.ge.irtrap.and.ir.le.nrs) then 
         vr_direction =  1.0
      else 
         vr_direction =  -1.0 
      endif

      ! Set the current particle time 
      current_time = cist * qtim 

      ! Check to see if a new value needs to be chosen - or if we are 
      ! still within the correlation time. This is only really intended
      ! and physically makes sense in the SOL (where blobs are), so 
      ! also make sure we are less than irsep.
      if (imp.ne.current_particle.or.current_time.ge.
     >  (last_time_chosen+dble(pinch_correlation_time)).and.(ir.ge.
     >  irsep)) then 
     
         ! At this point, we are not in the blob yet.
         in_blob = .false.

         ! sazmod
         ! Only pick a new velocity if a blob passes by. This is modeled
         ! via the probability fblob * qtim, i.e., the number of blobs
         ! seen between each timestep. This should be less than 1, 
         ! otherwise a new velocity is always chosen. If fblob = 0 then
         ! result_val will always = 0 (no blobs, no transport).
         nrand = nrand + 1
         ran = getranf()
         !write(0,*) 'ran, fblob*qtim = ',ran,fblob*qtim
         if (ran.le.(fblob*qtim)) then

           ! Get random number for fraction of integration
           nrand = nrand + 1 
           ran = getranf()

           ! Returns velocity from input distribution. Multiply by 
           ! direction factor.
           result_val = vr_pdf_random(ran) * vr_direction
           !write(0,*) ' choose vr = ',result_val

           ! Reset the last assigned and time chosen values
           current_particle = imp
           last_time_chosen = current_time

           ! Save the resultant velocity (including direction factor),
           ! set in_blob flag to true.
           vr_last_assigned = result_val
           in_blob = .true.
           
         ! If a velocity was not chosen from the distribution, then 
         ! assign zero radial transport.
         else
           !result_val = vr_last_assigned
           result_val = 0.0
         endif
      
      ! Set to zero if in core. The original constant pinch option
      ! happens below, so it's unaffected by this. 
      else if (ir.lt.irsep) then
        result_val = 0.0     
    
      else
         result_val = real(vr_last_assigned)
      endif
c
c      write(0,'(a,i6,6(1x,g15.8))')
c     >  'DEBUG VR:',imp,result_val,current_time,last_time_chosen,
c     >              pinch_correlation_time, 
c     >              last_time_chosen+pinch_correlation_time,
c     >              vr_last_assigned
c
c
c     Add debugging code for the value of pinchvel. 
c
      in = int(sign(
     >      min(( (abs(result_val)+d_pinch_vel/2.0)
     >           / d_pinch_vel),
     >           real(max_d_pinch_v)),
     >           result_val/vr_direction))
      d_pinch_v(in) = d_pinch_v(in) + 1.0 
      
!      write(0,*) 'current_particle, time, result_val = ', 
!     > current_particle, current_time, result_val

      ! sazmod - A constant pinch drift (cvpinch) can be specified
      ! on top of everything. cvpinch assumed to be in DIVIMP convention
      ! already (negative = outwards). 
      ! Multiply by unit ion time step.
      !find_vr = result_val * dble(qtim)
      find_vr = real((result_val + cvpinch) * dble(qtim))
      
      ! The limited number of blob measurements available suggest lower 
      ! radial velocities in the divertor, which here is approximated
      ! as above/below the X-point (USN/LSN). Pass in a factor to force lower 
      ! radial velocities here.
      if (.not.xpoint_up.and.zs(ik,ir).lt.zxp) then
        find_vr = find_vr * div_vr_fact
!        write(0,*) 'xpoint_up, zs, div_vr_fact = ',xpoint_up,
!     >    zs(ik,ir),div_vr_fact
      elseif (xpoint_up.and.zs(ik,ir).gt.zxp) then
        find_vr = find_vr * div_vr_fact
      endif

      ! Pathetic attempt at ballooning nature by multiplying by the
      ! factor (1/B^2 @ OMP) / (1/B^2) = B^2 / B^2 @ OMP. Toroidal
      ! field is fine since it's easily available. midplane_b can
      ! be zero if the ring does not cross the midplane (e.g. extended
      ! grids this could happen). The PFZ is given the corresponding
      ! core values. 
      if (balloon_opt.eq.1) then
        if (midplane_b(ir).ne.0.0) then
          find_vr = find_vr * (bts(ik,ir)*bts(ik,ir)) / 
     >      (midplane_b(ir) * midplane_b(ir))
!          write(0,*) 'ir,ik,bts,midplane_b,ratio',ir,ik,bts(ik,ir),
!     >      midplane_b(ir),bts(ik,ir)/midplane_b(ir) 
        endif
      endif

      return
      end 
c
c
c
      real function vr_pdf_random(ran)
      use mod_params
      use mod_comtor
      implicit none
      real ran
c     include 'params'
c     include 'comtor'
c
c     VR_PDF_RANDOM: 
c
c     Given a random number on [0,1] this code
c     calculates the value of velocity from the 
c     PDF which would correspond to this random 
c     number.
c
c     pinch_pdf_data(*,2) contains the normalized
c     integral of the PDF which is used for this
c     calculation - the random number is linearly 
c     interpolated on this array and a corresponding
c     velocity is returned. 
c
c     Local variables 
c 
c      real a,b,c
c      real v1,v2,x1,x2,i1,i2
c
      real*8 vr_pdf_int
      external vr_pdf_int
c
      real*8 eps
      parameter(eps=1.0d-8)
c
      integer iter_cnt,max_iter_cnt
      parameter(max_iter_cnt=10000)
c
      real*8 vmin,vmax,vtest,res_frac
c
      integer in,ipos
      external ipos  
c
c     Check that RAN is in range
c
      if (ran.le.0) then 
         vr_pdf_random = pinch_pdf_data(1,1)
      elseif (ran.ge.1.0) then     
         vr_pdf_random = pinch_pdf_data(npdf_data,1)
      endif    
c
c     Perform search loop - include maximum iteration test - set at 10,000 for now
c
      iter_cnt = 0
c
c     Find appropriate position in pinch_pdf_data array
c
      in = ipos(ran,pinch_pdf_data(1,3),npdf_data)
c
      res_frac = -1.0 
      vmin = pinch_pdf_data(in-1,1) 
      vmax = pinch_pdf_data(in,1) 
      vtest = (vmin+vmax)/2.0
c
      do while ((abs(res_frac-ran).gt.eps).and.
     >          (iter_cnt.lt.max_iter_cnt))
c
         iter_cnt = iter_cnt + 1
c
         res_frac = vr_pdf_int(vtest,in)/pdf_norm_val
c
c        Find new test value and iterate
c        - increase test value
c
         if (res_frac.lt.ran) then 
c
            vmin = vtest
            vtest = 0.5*(vmin+vmax)
c
c        Decrease test value
c
         else
c
            vmax = vtest
            vtest = 0.5*(vmin+vmax)
c
         endif
c
      end do

c
c     The value of the range for which the fractional condition is 
c     satisfied is stored in last_r
c
      vr_pdf_random = real(vtest)
c
c     Linearly interpolating the integrated probability results in an 
c     equal probability of all velocities in a given bin as opposed to 
c     one which varies linearly with the probability itself. 
c
c     vr_pdf_random = pinch_pdf_data(in-1,1) + 
c     >         ((ran-pinch_pdf_data(in-1,2))/
c     >         (pinch_pdf_data(in,2)-pinch_pdf_data(in-1,2)))*
c     >         (pinch_pdf_data(in,1)-pinch_pdf_data(in-1,1))
c
c      
c
c      write(6,'(a,2i6,6(1x,f14.7))') 'VR_PDF:',in,npdf_data,
c     >    pinch_pdf_data(in-1,3),ran, pinch_pdf_data(in,3),
c     >    pinch_pdf_data(in-1,1),vr_pdf_random,pinch_pdf_data(in,1)
c
      return
c
      end
      
c
c
c
      real*8 function vr_pdf(v_in)
      use mod_params
      use mod_comtor
      implicit none
      real*8 v_in
c     include 'params'
c     include 'comtor'
c
c     VR_PDF: This function returns the PDF 
c             value for the given value of 
c             V_IN. The PDF is specified as 
c             an ascending set of values 
c             
c             VEL_1   PDF_VALUE_1
c             VEL_2   PDF_VALUE_2
c                  ...
c             VEL_N   PDF_VALUE_N
c
c
c     The code linearly interpolates on the VEL
c     axis to obtain the PDF value for the specified 
c     velocity. For values less than VEL_1 or greater
c     than VEL_N the PDF value is returned as 0.0       
c
c     Local variables 
c 
      integer in,ipos
      external ipos  
c
c     Check V_IN in range
c
      if (v_in.lt.pinch_pdf_data(1,1).or.
     >    v_in.gt.pinch_pdf_data(npdf_data,1)) then 
         vr_pdf = 0.0d0
         return
      endif   
c
c     Find appropriate position in pinch_pdf_data array
c
      in = ipos(real(v_in),pinch_pdf_data(1,1),npdf_data)
c
c     Linearly interpolate
c
      vr_pdf = pinch_pdf_data(in-1,2) + 
     >        (v_in - pinch_pdf_data(in-1,1))/
     >        (pinch_pdf_data(in,1)-pinch_pdf_data(in-1,1)) *
     >        (pinch_pdf_data(in,2) - pinch_pdf_data(in-1,2))
c
      return
c
      end
c
c
c
      real*8 function vr_pdf_int(v_in,test_in)
      use mod_params
      use mod_comtor
      implicit none
      integer test_in 
      real*8     v_in
c     include 'params'
c     include 'comtor'
c
c     VR_PDF_INT: 
c  
c     This function calculates the integral of
c     the PDF at the velocity value V_in. 
c
c     Local variables 
c 
      real*8   pdfn
      integer  ipos,in
      external ipos  
c
c     Check V_IN in range
c
      if (v_in.le.pinch_pdf_data(1,1)) then 
         vr_pdf_int = 0.0
         return
      elseif (v_in.ge.pinch_pdf_data(npdf_data,1)) then 
         vr_pdf_int = pdf_norm_val
         return  
      endif 
c
c     Find appropriate position in pinch_pdf array
c
      if (test_in.le.0) then 

         in = ipos(real(v_in),pinch_pdf_data(1,1),npdf_data)

      else 
         
         in = test_in

      endif
c
c      write(6,'(a,2i4,3(1x,f14.8))') 'VR_PDF_INT 1:',
c     >       in,pinch_npdf,
c     >       pinch_pdf(in,1),v_in,pinch_pdf(in-1,1)
c
c
c     Calculate the integral at Vin
c     - find value of PDF at V_in then 
c       calculate the integral to that point
c
      pdfn =  pinch_pdf_data(in-1,2) + 
     >      ((v_in - pinch_pdf_data(in-1,1))/
     >      (pinch_pdf_data(in,1)-pinch_pdf_data(in-1,1))) *
     >      (pinch_pdf_data(in,2)-pinch_pdf_data(in-1,2))
c
      vr_pdf_int = pinch_pdf_data(in-1,3)*pdf_norm_val 
     >         + 0.5 *
     >        (v_in - pinch_pdf_data(in-1,1)) *
     >        (pdfn+pinch_pdf_data(in-1,2))
c
c
      return
c
      end
c
c
c
      subroutine set_pinch_velocity(ik,ir,nrand,imp,cist,pinchvel,
     >                              vr_assigned,ierr)
      use mod_params
      use mod_comtor
      use mod_cgeom
      implicit none
c
      real*8 cist
      integer ik,ir,nrand,imp,ierr
      real pinchvel
      logical vr_assigned
c
c     SET_PINCH_VELOCITY: This routine sets a cross field velocity based on
c                         options specified in the input file.
c
c     
c
c
      real find_vr
      external find_vr
c
c     include 'params'
c     include 'comtor'
c     include 'cgeom'
c
      vr_assigned = .false.
c
c         Pinch Velocity
c
          pinchvel = 0.0 
c
          if (pinchopt.eq.0) then 
c
             pinchvel = 0.0
c
          elseif (pinchopt.eq.1.or.
     >           (pinchopt.eq.2.and.
     >           (ir.ge.irsep.and.ir.le.irwall))) then 
c
             pinchvel = vpinch
c 
c         Pinch options 3 or 6 - load from grid pinch array
c
          elseif (pinchopt.eq.3.or.pinchopt.eq.6.or.pinchopt.eq.7.or.
     >            pinchopt.eq.8.or.pinchopt.eq.9.or.pinchopt.eq.10.or.
     >            pinchopt.eq.11.or.pinchopt.eq.12.or.pinchopt.eq.13.or.
     >            pinchopt.eq.14.or.pinchopt.eq.15)then 
c
             pinchvel = kpinchs(ik,ir)
c
          elseif (pinchopt.eq.4.or.pinchopt.eq.5) then 
c
             pinchvel=find_vr(ik,ir,nrand,vr_assigned,
     >                        cist,imp,ierr)
c
          endif


      return
      end
