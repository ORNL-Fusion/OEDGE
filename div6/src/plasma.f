      SUBROUTINE INITPLASMA(irstart,irend,ikopt)
      use debug_options
      IMPLICIT  none
c
      integer irstart,irend,ikopt
c
      include 'params'
c
      include 'comtor'
C
C  *********************************************************************
C  *                                                                   *
C  *  INITPLASMA:  THIS ROUTINE SETS THE INITIAL TEMPERATURES AND      * 
c  *               DENSITIES.                                          *
C  *                                                                   *
C  *  CHRIS FARRELL  (HUNTERSKIL)  JUNE     1989                       *
C  *  DAVID ELDER                  AUG/SEPT 1991                       *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IR
C
C     CALCULATE AND/OR ASSIGN THE PLASMA BACKGROUND TEMPERATURES
C     AND DENSITIES. FOR OPTIONS 2 AND 3 THE SOL DATA IS READ FROM
C     THE INPUT FILE AND THEN ASSIGNED. THE TEMPERATURES MAY BE MODIFIED
C     AND OR CALCULATED AGAIN IN THE TEMPERATURE GRADIENT SECTION.
C     THUS THIS ROUTINE IS EXECUTED TWICE FOR CERTAIN OPTIONS. THIS
C     IS NECESSARY IN ORDER TO DEFINE THE BASE VALUES FOR INBOARD
C     CHARACTERISTICS.   CTEB0,CTIB0 AND CNB0
c
c
c     Modified - use the base values for cteb0,ctibo and cnb0 from 
c                the input file for initial processing. These
c                will be overwritten depending on the later core
c                option selected. 
c
C
c      IF (CIOPTK.GE.2.AND.CIOPTK.LE.7) CTEB0 = 1.0E+02
c      IF (CIOPTL.GE.2.AND.CIOPTL.LE.7) CTIB0 = 1.0E+02
c      IF (CIOPTG.EQ.2.OR.CIOPTG.EQ.3.OR.CIOPTG.EQ.4.or.cioptg.eq.5
c     >     .or.cioptg.eq.7 
c     >     .OR.CIOPTF.EQ.12.OR.CIOPTF.EQ.13.OR.CIOPTF.EQ.14
c     >     .OR.CIOPTF.EQ.15.or.cioptf.eq.16.or.cioptf.eq.17
c     >     .or.cioptf.eq.18.or.cioptf.eq.19.or.cioptf.eq.20
c     >     .or.cioptf.eq.21.or.cioptf.eq.22.or.cioptf.eq.23) THEN
c                CTEb0 = 1.0E02
c                CTIb0 = 1.0E02
c                CNb0  = 1.0E19
c      ENDIF
C
      call pr_trace('PLASMA','START OF INITPLASMA')

      DO IR = irstart, irend
c
        call set_initplasma(ir,ikopt)
c
      end do  
c
      call pr_trace('PLASMA','END OF INITPLASMA')
c
      return 
      end 
c
c
c
      subroutine sol_plasma(irstart,irend,ikopt)
      use debug_options
      IMPLICIT  none
c
      integer irstart,irend,ikopt
c
      include 'params'
c
      include 'cgeom'
c
      include 'comtor'
C
C  *********************************************************************
C  *                                                                   *
C  *  SOL_PLASMA:  THIS ROUTINE SETS THE EVOLUTION OF TEMPERATURES AND *
c  *               DENSITIES ALONG THE FIELD LINES IN THE SOL AND PP.  *
C  *                                                                   *
C  *  DAVID ELDER                  JULY  1997                          *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IK,IR,IKMID,irlim1,irlim2,in,ierr
c
c
C-----------------------------------------------------------------------
C     CALCULATE TEMPERATURE GRADIENTS IN SOL ACCORDING TO OPTIONS
C
C     FOR OPTIONS 2,3,4 AND 5 THIS MEANS CALCULATING THE TEMPERATURES
C     IN THE SOL BASED ON GIVEN PLATE VALUES. THIS MEANS THE THE
C     TEMPERATURE USED FOR CALCULATING INBOARD VALUES CTEB0,CTIB0 ARE
C     INITIALLY NOT DEFINED AND THE SUBROUTINE NEEDS TO BE INVOKED
C     ITERATIVELY.
C
C-----------------------------------------------------------------------
C
C     NOTE : THE OPTION CIOPTO CONTROLS WHETHER OR NOT THE
C            GRADIENT OPTIONS ARE APPLIED TO THE TRAPPED PLASMA
C
C           NOTE: THERE WAS SOME CONFUSION ABOUT THE DESIGNATION OF
C           INNER AND OUTER PLATES. THE ARRAYS LPDATO - OUTER AND LPDATI
C           INNER CONTAIN THE SPECIFIC PLATE INFORMATION. FOR UNIFORM CA
C           LPDATI CONTAINS THE INFORMATION FOR THE BOTH TARGETS.
C
C
      call pr_trace('PLASMA','START OF SOL_PLASMA')
c

      irlim1 = irstart
c     
      if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
        IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >      ciopto.eq.3.or.ciopto.eq.4) THEN
          IRLIM2 = min(IRWALL,irend)
        ELSEIF (CIOPTO.EQ.1) THEN
          IRLIM2 = min(NRS,irend)
        ENDIF
      elseif (cgridopt.eq.2) then
        IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >      ciopto.eq.3.or.ciopto.eq.4) THEN
          IRLIM2 = IRWALL2
        ELSEIF (CIOPTO.EQ.1) THEN
          IRLIM2 = NRS2
        ENDIF
c slmod begin
      elseif (cgridopt.eq.LINEAR_GRID.or.cgridopt.eq.RIBBON_GRID) then
        IRLIM2 = min(nrs,irend)
c slmod end
      endif
c
c     For testing SOL options on one ring only.
c
      if (ctestsol.gt.0.0) then 
         irlim1 = int(ctestsol)   
         irlim2 = int(ctestsol) 
      endif

c
      IF (CIOPTF.EQ.12.OR.CIOPTF.EQ.13.OR.CIOPTF.EQ.14.OR.CIOPTF.EQ.15
     >    .or.cioptf.eq.16.or.cioptf.eq.17.or.cioptf.eq.18.or.
     >    cioptf.eq.19.or.cioptf.eq.20.or.cioptf.eq.21)
     >        THEN
C
C        CALL THE ROUTINE TO SET UP N,TE,TI,V,E FOR SOL 12,13,14 OR 15
C
C
         CALL soledge(IRlim1,IRLIM2,ikopt)
c

         if (cgridopt.eq.2) then
           irlim1 = irsep2
           IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >         ciopto.eq.3.or.ciopto.eq.4) THEN
             IRLIM2 = IRWALL
           ELSEIF (CIOPTO.EQ.1) THEN
             IRLIM2 = NRS
           ENDIF
           CALL soledge(IRlim1,IRLIM2,ikopt)
         endif
C
      ELSEif (cioptf.eq.22) then
c
         call calcsol_interface(irlim1,irlim2,ikopt)
c
         if (cgridopt.eq.2) then
           irlim1 = irsep2
           IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >         ciopto.eq.3.or.ciopto.eq.4) THEN
             IRLIM2 = IRWALL
           ELSEIF (CIOPTO.EQ.1) THEN
             IRLIM2 = NRS
           ENDIF
           CALL calcsol_interface(IRlim1,IRLIM2,ikopt)
         endif
c
      ELSEif (cioptf.eq.23) then
c
         call sol23_interface(irlim1,irlim2)
c
      ELSEif (cioptf.eq.24) then

c
         call specplas(irlim1,irlim2,ikopt)
c
c slmod begin
      ELSEIF (cioptf.EQ.28) THEN
         CALL SOL28(irlim1,irlim2,ikopt)
c slmod end
      else
C
         CALL TGRAD(IRlim1,IRLIM2,ikopt)
         call Ngrad(irlim1,irlim2,ikopt)
c
         if (cgridopt.eq.2) then
           irlim1 = irsep2
           IF (CIOPTO.EQ.0.or.ciopto.eq.2.or.
     >         ciopto.eq.3.or.ciopto.eq.4) THEN
             IRLIM2 = IRWALL
           ELSEIF (CIOPTO.EQ.1) THEN
             IRLIM2 = NRS
           ENDIF
           CALL TGRAD(IRlim1,IRLIM2,ikopt)
           call Ngrad(irlim1,irlim2,ikopt)
         endif
      endif

c
c     If specified private plasma is ON - then it calls the 
c     routine to calculate the private plasma values.
c
      if (ciopto.eq.2) then 
         call specplas(irtrap,nrs,ikopt) 
      elseif (ciopto.eq.3.or.ciopto.eq.4) then 
         call thompp(irtrap,nrs,ikopt,ciopto) 
      endif
c
      call pr_trace('PLASMA','END OF SOL_PLASMA')


      return 
      end
c
c
c
      subroutine core_plasma(irstart,irend,ikopt)
      use debug_options
      IMPLICIT  none
c
      integer irstart,irend,ikopt
c
      include 'params'
c
      include 'cgeom'
c
      include 'comtor'
C
C  *********************************************************************
C  *                                                                   *
C  *  CORE_PLASMA:  THIS ROUTINE SETS THE EVOLUTION OF TEMPERATURES AND*
c  *                DENSITIES ALONG THE FIELD LINES IN THE CORE.        *
C  *                                                                   *
C  *  DAVID ELDER                  JULY  1997                          *
C  *                                                                   *
C  *********************************************************************
C
      INTEGER IK,IR,IKMID,in,ierr
c
c
C
C     This is the processing that needs to be repeated. This
C     is generally only the core plasma conditions after
C     the SOL has been calculated. So for the appropriate
C     options repeat the calculations.
C
c
c     NOTE: The following may need to be changed for ITER grids.
c           The implementation below will use the midpoint
c           of the inner separatrix ring as the base temperature
c           for the core plasma ... it might be desirable to
c           have the average of inner and outer separatrix
c           rings or some other formula in the future.
c
c     Local variables
c
      real nstart,testart,tistart,vstart,nmid,temid,timid,stot,s
      integer ringno,nr,ikstart,ikend
      external ringno
C
c     Find mid-point of separatrix.      
c

      call pr_trace('PLASMA','START OF CORE_PLASMA')

      ikmid = ikmids(irsep)+1
c
c     The following turns off all core option processing.
c
      if (ccoreopt.eq.-1) return 
c
c     Continue with CORE option processing 
c
c     Calculate core conditions - start with mid-separatrix conditions
c
      IF (CIOPTK.GE.2.AND.CIOPTK.LE.7)
     >          CTEB0 = (KTEBS(IKMID,IRSEP)+KTEBS(IKMID-1,IRSEP))/2.0
      IF (CIOPTL.GE.2.AND.CIOPTL.LE.7)
     >          CTIB0 = (KTIBS(IKMID,IRSEP)+KTIBS(IKMID-1,IRSEP))/2.0
      IF (CIOPTG.EQ.2.OR.CIOPTG.EQ.3.OR.CIOPTG.EQ.4
     >           .or.cioptg.eq.7
     >           .OR.CIOPTF.EQ.12.OR.CIOPTF.EQ.13.OR.CIOPTF.EQ.14
     >           .OR.CIOPTF.EQ.15.or.cioptf.eq.16.or.cioptf.eq.17
     >           .or.cioptf.eq.18.or.cioptf.eq.19.or.cioptf.eq.20
     >           .or.cioptf.eq.21.or.cioptf.eq.22.or.cioptf.eq.23) THEN
         CTEB0 = (KTEBS(IKMID,IRSEP)+KTEBS(IKMID-1,IRSEP))/2.0
         CTIB0 = (KTIBS(IKMID,IRSEP)+KTIBS(IKMID-1,IRSEP))/2.0
         CNB0  = (KNBS(IKMID,IRSEP)+KNBS(IKMID-1,IRSEP))/2.0
      ENDIF
c
c     
c
c slmod begin
      if (ccoreopt.eq.28) then
c...     Allow SOL28 to be called for core rings:
         call sol28(irstart,irend,ikopt)
      elseif (ccoreopt.eq.0.and.cioptg.ne.7) then 
c
c      IF (ccoreopt.eq.0.and.cioptg.ne.7) then 
c slmod end

c            .and.(
c     >       CIOPTK.EQ.2.OR.CIOPTK.EQ.3.OR.CIOPTK.EQ.4.OR.CIOPTK.EQ.5
c     >      .OR.CIOPTK.EQ.6.OR.CIOPTK.EQ.7
c     >      .OR.CIOPTL.EQ.2.OR.CIOPTL.EQ.3.OR.CIOPTL.EQ.4.OR.CIOPTL.EQ.5
c     >      .OR.CIOPTL.EQ.6.OR.CIOPTL.EQ.7
c     >      .OR.CIOPTG.EQ.2.OR.CIOPTG.EQ.3.OR.CIOPTG.EQ.4
c     >      .OR.CIOPTF.EQ.12.OR.CIOPTF.EQ.13.OR.CIOPTF.EQ.14
c     >      .OR.CIOPTF.EQ.15.or.cioptf.eq.16.or.cioptf.eq.17
c     >      .or.cioptf.eq.18.or.cioptf.eq.19.or.cioptf.eq.20
c     >      .or.cioptf.eq.21.or.cioptf.eq.22)) THEN
C
         DO 500 IR = irstart, irend
c
            call set_ikvals(ir,ikstart,ikend,ikopt)
c
            DO 500 IK = ikstart, ikend
C
C              Recalculate core conditions
C
               KTEBS(IK,IR) = CTEB0 + REAL(IRSEP-IR) * CTEBIN
               KTIBS(IK,IR) = CTIB0 + REAL(IRSEP-IR) * CTIBIN
               KNBS (IK,IR) = CNB0  + REAL(IRSEP-IR) * CNBIN
c
 500     continue
c
c     Assign constant values to the core plasma from input file.     
c
      elseif (ccoreopt.eq.1) then
c
c        Loop through rings               
c
         do ir = irstart,irend
c
            call set_ikvals(ir,ikstart,ikend,ikopt)
c
            NR = RINGNO(IR,coreDAT,NcoreDAT,MAXINS,5,ierr)
c
            do ik = ikstart,ikend
c
               ktebs(ik,ir) = coredat(nr,2)
               ktibs(ik,ir) = coredat(nr,3)
               knbs(ik,ir)  = coredat(nr,4)
               kvhs(ik,ir)  = coredat(nr,5)

            end do
c
         end do
c
c     Marfe Core plasma description - options 2 and 3
c
      elseif (ccoreopt.eq.2.or.ccoreopt.eq.3) then       
c
c        Loop through rings 
c
         do ir = irstart,irend   
c
            ikmid = ikmids(ir)
c
            call set_ikvals(ir,ikstart,ikend,ikopt)
c
c            if (ikend.eq.nks(ir)) ikend = ikend-1
c          
            NR = RINGNO(IR,coreDAT,NcoreDAT,MAXINS,5,ierr)
c          
            stot = kss(nks(ir)-1,ir)
c
            testart = coreDAT(NR,2)
            tistart = coreDAT(NR,3)
            nstart  = coreDAT(NR,4)
            vstart  = coredat(nr,5) 
c
            nmid = knbs(ikmid,ir)
            temid = ktebs(ikmid,ir)
            timid = ktibs(ikmid,ir)
c            
            do ik = ikstart,ikend
c           
               s = kss(ik,ir) 
c           
c              Constant from top to midplane - ramps down to 
c              entered value at the X-point. 
c           
c              NOTE: Data entered for INNER and OUTER should be 
c                    the same.
c            
c              At the present time there are two cells adjacent
c              at the X-point - these are index 1 and index
c              nks(ir)-1 ... the index nks(ir) points to the same
c              cell as 1 - in order to give a closed ring.
c            
               if (ierr.eq.0.and.kss(ik,ir).lt.0.25*stot) then 
c            
                  ktebs(ik,ir) = s/(0.25*stot)  
     >              * (temid - testart)
     >              + testart
                  ktibs(ik,ir) = s/(0.25*stot)  
     >              * (timid - tistart)
     >              + tistart
                  kvhs(ik,ir) = (1.0- s/(0.25*stot) ) * vstart
                  kes(ik,ir) = 0.0
c
                  if (ccoreopt.eq.2) then 

                     knbs(ik,ir) = s/(0.25*stot)  
     >                 * (nmid - nstart)
     >                 + nstart
c
                  elseif (ccoreopt.eq.3) then 
c
                     knbs(ik,ir) = nmid * (temid+timid) * ech /
     >                      ( (ktebs(ik,ir)+ktibs(ik,ir)) * ech + 
     >                            crmb * amu * kvhs(ik,ir)**2)
c
                  endif 
c            
               elseif (ierr.eq.0.and.kss(ik,ir).gt.0.75*stot) then 
c           
                  ktebs(ik,ir) = (stot-s)/(0.25*stot)  
     >              * (temid - testart)
     >              + testart
                  ktibs(ik,ir) = (stot-s)/(0.25*stot)  
     >              * (timid - tistart)
     >              + tistart
                  kvhs(ik,ir) = -(1.0-(stot-s)/(0.25*stot))
     >              * vstart
                  kes(ik,ir) = 0.0
c
                  if (ccoreopt.eq.2) then 
c
                     knbs(ik,ir) = (stot-s)/(0.25*stot)  
     >                 * (nmid - nstart)
     >                 + nstart
c
                  elseif (ccoreopt.eq.3) then 
c
                     knbs(ik,ir) = nmid * (temid+timid) * ech /
     >                      ( (ktebs(ik,ir)+ktibs(ik,ir)) * ech + 
     >                            crmb * amu * kvhs(ik,ir)**2)
c
                  endif 
c	    
               else 
c
                  KTEBS(IK,IR) = temid
                  KTIBS(IK,IR) = timid
                  KNBS(IK,IR)  = nmid
                  kvhs(ik,ir)  = 0.0
                  kes(ik,ir)   = 0.0
c
               endif 
c
            end do
c
         end do
c
c
      elseif (ccoreopt.eq.4.or.ccoreopt.eq.5.or.ccoreopt.eq.6) then 
c
c        Core Option 4 - Marfe - variable distances
c
         do ir = irstart,irend
c
            ikmid = ikmids(ir)
c
            call set_ikvals(ir,ikstart,ikend,ikopt)
c          
            NR = RINGNO(IR,coreDAT,NcoreDAT,MAXINS,5,ierr)
c          
            stot = kss(nks(ir)-1,ir)
c            
            testart = coreDAT(NR,2)
            tistart = coreDAT(NR,3)
            nstart  = coreDAT(NR,4)
            vstart  = coredat(nr,5) 
c
            nmid = knbs(ikmid,ir)
            temid = ktebs(ikmid,ir)
            timid = ktibs(ikmid,ir)
c
            do ik = ikstart,ikend
c           
               s = kss(ik,ir) 
c           
c              Constant from top to specified points - ramps down to 
c              entered value at the X-point. 
c           
c              NOTE: Data entered for INNER and OUTER should be 
c                    the same.
c            
c              At the present time there are two cells adjacent
c              at the X-point - these are index 1 and index
c              nks(ir)-1 ... the index nks(ir) points to the same
c              cell as 1 - in order to give a closed ring.
c            
               if (ierr.eq.0) then 
c
c                 Electric Field 
c
                  kes(ik,ir) = 0.0 
c
c                 Velocity
c
                  if (kss(ik,ir).lt.corefv*stot.or.
     >                kss(ik,ir).gt.(1.0-corefv)*stot) then   
c
                     kvhs(ik,ir) = 0.0
c
                  elseif (kss(ik,ir).ge.corefv*stot.and.
     >                    kss(ik,ir).lt.corefv2*stot) then  
c
                      kvhs(ik,ir)=(1.0-(s-corefv*stot)
     >                              / ((corefv2-corefv)*stot))
     >                              * vstart
c
                  elseif (kss(ik,ir).gt.(1.0-corefv2)*stot.and.
     >                    kss(ik,ir).le.(1.0-corefv)*stot) then 
c
                     kvhs(ik,ir) = -(1.0-(stot-(s+corefv*stot))
     >                              /((corefv2-corefv)*stot))
     >                              * vstart
c
                  else
c
                     kvhs(ik,ir) = 0.0
              
                  endif
c
c                 Temperature
c
                  if (kss(ik,ir).lt.coreft*stot.or.
     >                kss(ik,ir).gt.(1.0-coreft)*stot) then   

                     ktebs(ik,ir) = testart
                     ktibs(ik,ir) = tistart 

                  elseif (kss(ik,ir).ge.coreft*stot.and.
     >                    kss(ik,ir).lt.coreft2*stot) then  

                     ktebs(ik,ir) = (s-coreft*stot)
     >                 / ((coreft2-coreft)*stot)  
     >                 * (temid - testart)
     >                 + testart
                     ktibs(ik,ir) = (s-coreft*stot)
     >                 / ((coreft2-coreft)*stot)  
     >                 * (timid - tistart)
     >                 + tistart
c

                  elseif (kss(ik,ir).gt.(1.0-coreft2)*stot.and.
     >                    kss(ik,ir).le.(1.0-coreft)*stot) then 
c
                     ktebs(ik,ir) = (stot-(s+coreft*stot))
     >                 / ((coreft2-coreft)*stot)  
     >                 * (temid  - testart)
     >                 + testart
                     ktibs(ik,ir) = (stot-(s+coreft*stot))
     >                 / ((coreft2-coreft)*stot)  
     >                 * (timid  - tistart)
     >                 + tistart
c
                  else
c
                     KTEBS(IK,IR) = temid
                     KTIBS(IK,IR) = timid
c
                  endif
c            
c                 Density
c
c                 Add an option for ramp of plasma density as well as pressure conservation 
c                 The density ramp uses the same parameters as the temperature for 
c                 along ring locations
c
                  if (ccoreopt.eq.6) then 
c
                     if (kss(ik,ir).lt.coreft*stot.or.
     >                   kss(ik,ir).gt.(1.0-coreft)*stot) then   
                                          
                        knbs(ik,ir) = nstart
                                          
                     elseif (kss(ik,ir).ge.coreft*stot.and.
     >                       kss(ik,ir).lt.coreft2*stot) then  
                                          
                        knbs(ik,ir) = (s-coreft*stot)
     >                    / ((coreft2-coreft)*stot)  
     >                    * (nmid - nstart)
     >                    + nstart
c                                       
                                          
                     elseif (kss(ik,ir).gt.(1.0-coreft2)*stot.and.
     >                       kss(ik,ir).le.(1.0-coreft)*stot) then 
c                                       
                        knbs(ik,ir) = (stot-(s+coreft*stot))
     >                    / ((coreft2-coreft)*stot)  
     >                    * (nmid  - nstart)
     >                    + nstart
c                                       
                     else
c                                       
                        KNBS(IK,IR) = nmid
c                                       
                     endif
c

                  else
                     knbs(ik,ir) = nmid * (temid+timid) * ech /
     >                      ( (ktebs(ik,ir)+ktibs(ik,ir)) * ech + 
     >                            crmb * amu * kvhs(ik,ir)**2)
                  endif


c
               else 
c
                  KTEBS(IK,IR) = temid
                  KTIBS(IK,IR) = timid
                  KNBS(IK,IR)  = nmid
                  kvhs(ik,ir)  = 0.0
                  kes(ik,ir)   = 0.0
c
               endif 
c
            end do
c
         end do
c
      endif
c     
c     Set last point values to first point since ring is closed.
c
      do ir = irstart,irend
         KTEBS(nks(ir),IR) = ktebs(1,ir)
         KTIBS(nks(ir),IR) = ktibs(1,ir)
         KNBS(nks(ir),IR)  = knbs(1,ir)
         kvhs(nks(ir),ir)  = kvhs(1,ir)
         kes(nks(ir),ir)   = kes(1,ir)
      end do 
c
      call pr_trace('PLASMA','END OF CORE_PLASMA')
C
C
      RETURN
      END
C
C
C
      INTEGER FUNCTION RINGNO(IR,LPDAT,NLPDAT,MAXN1,MAXN2,ierr)
      implicit none
      INTEGER IR,NLPDAT
      INTEGER MAXN1,MAXN2,ierr
      REAL LPDAT(MAXN1,MAXN2)
C
      INTEGER I
C
C     THIS ROUTINE CHECKS THROUGH THE DATA ARRAY TO SEE IF THE
C     FIRST VALUE MATCHES THE RING NUMBER BEING SOUGHT. IF IT MATCHES
C     IT RETURNS THE INDEX INTO THE ARRAY, IF IT DOESN'T MATCH IT RETURN
C     THE INDEX OF THE LAST SET OF DATA.
C
      ierr = 0
c
      DO 100 I = 1,NLPDAT
         IF (INT(LPDAT(I,1)).EQ.IR) THEN
            GOTO 200
         ENDIF
100   CONTINUE
      I=NLPDAT
      WRITE(6,'(a,i5,a,i5,a)') 'DATA NOT FOUND FOR RING - ',IR,
     >           ' RETURNED DATA FOR RING ',
     >           INT(LPDAT(I,1)),' INSTEAD'
      ierr = 1
200   RINGNO = I
      RETURN
      END
c
c
c
      subroutine NGRAD(IRLIM1,IRLIM2,ikopt)
      IMPLICIT NONE
      INTEGER IRLIM1,IRLIM2,ikopt
C     INCLUDE   "PARAMS"
      include 'params'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C     INCLUDE   "COMTOR"
      include 'comtor'
C
c
C  *********************************************************************
c  *                                                                   *
c  *   NGRAD:                                                          *
c  *                                                                   *
c  *   This subroutine applies the selected density gradient           *
c  *   option to the range of rings from irlim1 to irlim2              *
c  *   inclusive.                                                      *
c  *                                                                   *
C  *********************************************************************
c
      INTEGER IK,IR,ikstart,ikend
      REAL    S,SMAX,TEB,TIB,nb
      real    tebi,tibi,tebo,tibo,nbo,nbi
c
c     Start looping through rings
c
       DO IR = irlim1, IRLIM2
c
c         Set IK values for inner loops
c
          call set_ikvals(ir,ikstart,ikend,ikopt)
c
          SMAX = KSMAXS(IR)
c 
c         Pull out target conditions locally
c
          TEBo  = KTeds(idds(ir,2))
          TIBo  = KTIds(idds(ir,2))
          tebi  = kteds(idds(ir,1))
          tibi  = ktids(idds(ir,1))
          NBO   = KndS(idds(ir,2))
          NBI   = KNdS(idds(ir,1))
c
C         Loop over specified section of ring
c 
          DO IK = ikstart, ikend
c
             S  = KSS(IK,IR)
c
             IF (S.GE.0.5*SMAX) THEN
                S = SMAX - S
                NB  = NBI
                TEB = tebi
                TIB = tibi 
             else
                NB  = NBO
                TEB = tebo
                TIB = tibo 
             endif 
c
c            Calculate density in given cell 
c
             if (ngradopt.eq.1) then

                knbs(ik,ir) = 4.0 * nb * teb / 
     >                     (ktebs(ik,ir)+ktibs(ik,ir)) 

             endif  
c
          end do 
c     
c         Feb 2008 - jde - not sure why this code is needed. At the 
c                    beginning of this routine - nbo and nbi are 
c                    assigned from the contents of knds. The conditions
c                    at the target surface are expected to be different
c                    from those at the center of the first cell - so it 
c                    seems counter intuitive to assign the first cell
c                    center contents back to the target condition array
c                    at the end of this routine. So I am commenting out 
c                    this code for the time being.  
c
c         IPP/08 Krieger - introduced code to ensure that densities
c         calculated by NGRAD for SOL also show up at target plates
c         Actually not clear if still needed. Similar add-on for
c         routine TGRAD by Klaus Schmid has been commented out.
c          if (ikopt.eq.1.or.ikopt.eq.3) then            
c             KNDS(IDDS(IR,2)) =  KNBS(1,IR)
c          endif
c          if (ikopt.eq.2.or.ikopt.eq.3) then  
c             KNDS(IDDS(IR,1)) =  KNBS(NKS(IR),IR)
c          endif
c
      end do 
c


      return
      end
c
c
c
c
      SUBROUTINE TGRAD(IRLIM1,IRLIM2,ikopt)
      IMPLICIT NONE
      INTEGER IRLIM1,IRLIM2,ikopt
C     INCLUDE   "PARAMS"
      include 'params'
C     INCLUDE   "CGEOM"
      include 'cgeom'
C     INCLUDE   "COMTOR"
      include 'comtor'
C
c
C  *********************************************************************
c  *                                                                   *
c  *   TGRAD:                                                          *
c  *                                                                   *
c  *   This subroutine applies the selected temperature gradient       *
c  *   options to the range of rings from irlim1 to irlim2             *
c  *   inclusive. Extracting this routine into a separate              *
c  *   module became necessary when separate ranges of                 *
c  *   rings needed to be treated as SOL and others as                 *
c  *   TRAP ... this situation arises when dealing with ITER           *
c  *   grids.                                                          *
c  *                                                                   *
C  *********************************************************************
c
      INTEGER IK,IR,NR,NRI, IKMID,ikstart,ikend
      REAL    S,SMAX,TEB,TIB,SE1,SE2,SI1,SI2,TEBP,TIBP,NBP
      REAL    LPPA,TEBPI,NBPI,TIBPI,LPPAI
      real    tebi,tibi,tebo,tibo
      REAL    LPPAEB,LPPAEI,LPPAIB,LPPAII
      LOGICAL INPLATE
      INTEGER RINGNO
      EXTERNAL RINGNO
c
c     Start looping through rings
c
       DO 400 IR = irlim1, IRLIM2
c
c        Set IK values for inner loops
c
         call set_ikvals(ir,ikstart,ikend,ikopt)
c
        SMAX = KSMAXS(IR)
        SE1  = CFEBL1 * SMAX
        SE2  = CFEBL2 * SMAX
        SI1  = CFIBL1 * SMAX
        SI2  = CFIBL2 * SMAX
c
        TEBo  = KTEdS(idds(ir,2))
        TIBo  = KTIds(idds(ir,2))
        tebi = kteds(idds(ir,1))
        tibi = ktids(idds(ir,1))
        NBP  = KNDS(idds(ir,2))
c
c        TEBo  = KTEBS(1,IR)
c        TIBo  = KTIBS(1,IR)
c        tebi = ktebs(nks(ir),ir)
c        tibi = ktibs(nks(ir),ir)
c        NBP  = KNBS(1,IR)
c
        IF (CIOPTG.EQ.3.OR.CIOPTG.EQ.4) THEN
c           NBPI = KNBS(NKS(IR),IR)
           NBPI = KNdS(idds(ir,1))
        ELSE
           NBPI = NBP
        ENDIF
C
C       SET UP ELECTRON PLATE TEMPERATURES
C
        IF (CIOPTK.EQ.0.OR.CIOPTK.EQ.1) THEN
           TEBP = Tebo * CFEBT
           TEBPI = tebi * CFEBT
        ELSEIF (CIOPTK.EQ.2) THEN
           if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
             IF (IR.LE.IRWALL) THEN
                TEBP = CTEBP - FLOAT(IR-IRSEP)*CTEBOU
             ELSE
                TEBP = CTEBP - FLOAT(NRS-IR+1)*CTEBOU
             ENDIF
           elseif (cgridopt.eq.2) then
             IF (IR.LE.IRWALL2) THEN
                TEBP = CTEBP - FLOAT(IR-IRSEP)*CTEBOU
             ELSEIF(IR.LE.NRS2) then
                TEBP = CTEBP - FLOAT(NRS2-IR+1)*CTEBOU
             ElseIF(IR.LE.IRWALL2) then
                TEBP = CTEBP - FLOAT(IR-IRSEP2)*CTEBOU
             ELSEIF(IR.LE.NRS) then
                TEBP = CTEBP - FLOAT(NRS-IR+1)*CTEBOU
             ENDIF
           endif
           TEBPI = TEBP
           IF (TEBP.LT.1.0E+00) TEBP = 1.0E+00
           IF (TEBPI.LT.1.0E+00) TEBPI = 1.0E+00
        ELSEIF (CIOPTK.EQ.3) THEN
           TEBP = tebo
           TEBPI = TEBP
        ELSEIF (CIOPTK.EQ.4.OR.CIOPTK.EQ.5
     >          .OR.CIOPTK.EQ.6.OR.CIOPTK.EQ.7) THEN
           TEBP = tebo
           TEBPI = tebi
        ENDIF
C
C       SET UP ION PLATE TEMPERATURES
C
        IF (CIOPTL.EQ.0.OR.CIOPTL.EQ.1) THEN
           TIBP = tibo * CFIBT
           TIBPI = tibi * CFIBT
        ELSEIF (CIOPTL.EQ.2) THEN
           if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3) then
             IF (IR.LE.IRWALL) THEN
                TIBP = CTIBP - FLOAT(IR-IRSEP)*CTIBOU
             ELSE
                TIBP = CTIBP - FLOAT(NRS-IR+1)*CTIBOU
             ENDIF
           elseif (cgridopt.eq.2) then
             IF (IR.LE.IRWALL2) THEN
                TIBP = CTIBP - FLOAT(IR-IRSEP)*CTIBOU
             ELSEIF (IR.LE.NRS2) then
                TIBP = CTIBP - FLOAT(NRS2-IR+1)*CTIBOU
             ELSEIF (IR.LE.IRWALL) then
                TIBP = CTIBP - FLOAT(IR-IRSEP2)*CTIBOU
             ELSEIF (IR.LE.NRS) then
                TIBP = CTIBP - FLOAT(NRS-IR+1)*CTIBOU
             ENDIF
           endif
           TIBPI = TIBP
           IF (TIBP.LT.1.0E+00) TIBP = 1.0E+00
           IF (TIBPI.LT.1.0E+00) TIBPI = 1.0E+00
        ELSEIF (CIOPTL.EQ.3) THEN
           TIBP = tibo
           TIBPI = TIBP
        ELSEIF (CIOPTL.EQ.4.OR.CIOPTL.EQ.5
     >          .OR.CIOPTL.EQ.6.OR.CIOPTL.EQ.7) THEN
           TIBP = tibo
           TIBPI = tibi
        ENDIF
C
c
c       NOTE: PLATE TEMPERATURES ARE NOW STORED IN SET_INITPLASMA
c
C       SAVE PLATE TEMPERATURES AND DENSITY FOR EACH LAUNCH POSITION
C
C       THE ARRAY IDDS GIVES THE NDS INDICES OF THE PLATE POINTS FOR
C       EACH RING. THE INNER PLATE IS 1 AND THE OUTER IS 2.
C
C
c       IPP/08 Krieger - introduced code to ensure that temperatures
c       calculated by TGRAD for SOL also show up at target plates
c       This mod by Klaus Schmid has been commented out. Not needed
c       anymore ?
C
c        if (ikopt.eq.1.or.ikopt.eq.3) then  
c
c           KTEDS(IDDS(IR,2)) = TEBP
c           KTIDS(IDDS(IR,2)) = TIBP
c           KNDS(IDDS(IR,2)) = NBP
c
c        endif
c
c        if (ikopt.eq.2.or.ikopt.eq.3) then  
c
c           KTEDS(IDDS(IR,1)) = TEBPI
c           KTIDS(IDDS(IR,1)) = TIBPI
c           KNDS(IDDS(IR,1)) = NBPI
c
c        endif
c
c       Calculate power coefficients
C
        IF (CIOPTK.EQ.3.OR.CIOPTL.EQ.3) THEN
           LPPA = (2.5*TIBP+5.0*TEBP)*1.602192E-19*NBP*
     >            SQRT((TIBP+TEBP)*EMI/CRMB)
        ENDIF
        IF (CIOPTK.EQ.4.OR.CIOPTL.EQ.4
     >      .OR.CIOPTK.EQ.6.OR.CIOPTL.EQ.6) THEN
           LPPA = (2.5*TIBP+5.0*TEBP)*1.602192E-19*NBP*
     >            SQRT((TIBP+TEBP)*EMI/CRMB)
           LPPAI = (2.5*TIBPI+5.0*TEBPI)*1.602192E-19*NBPI*
     >            SQRT((TIBPI+TEBPI)*EMI/CRMB)
        ENDIF
        IF (CIOPTK.EQ.5.OR.CIOPTL.EQ.5
     >      .OR.CIOPTK.EQ.7.OR.CIOPTL.EQ.7) THEN
           LPPAEB = 5.0*TEBP*1.602192E-19*NBP*
     >            SQRT((TIBP+TEBP)*EMI/CRMB)
           LPPAEI = 5.0*TEBPI*1.602192E-19*NBPI*
     >            SQRT((TIBPI+TEBPI)*EMI/CRMB)
           LPPAIB = 2.5*TIBP*1.602192E-19*NBP*
     >            SQRT((TIBP+TEBP)*EMI/CRMB)
           LPPAII = 2.5*TIBPI*1.602192E-19*NBPI*
     >            SQRT((TIBPI+TEBPI)*EMI/CRMB)
        ENDIF
C
C       Loop over specified section of ring
c 
        DO 300 IK = ikstart, ikend
          INPLATE = .FALSE.
          S  = KSS(IK,IR)
          TEB = TEBo
          TIB = TIBo
          IF (S.GE.0.5*SMAX) THEN
             S = SMAX - S
             INPLATE = .TRUE.
             TEB = TEBi
             TIB = TIBi
          ENDIF
C
          IF     (CIOPTK.EQ.0) THEN
            IF     (S.LE.0.0) THEN
              KTEBS(IK,IR) = TEB * CFEBT
            ELSEIF (S.LE.SE1) THEN
              KTEBS(IK,IR) = TEB * (CFEBT+(1.0-CFEBT)*S/SE1)
            ELSE
              KTEBS(IK,IR) = TEB
            ENDIF
          ELSEIF (CIOPTK.EQ.1) THEN
            IF     (S.LE.0.0) THEN
              KTEBS(IK,IR) = TEB * CFEBT
            ELSEIF (S.LE.SE1) THEN
              KTEBS(IK,IR) = TEB * (CFEBT+(CFEB2-CFEBT)*S/SE1)
            ELSEIF (S.LE.SE2) THEN
              KTEBS(IK,IR) = TEB * (CFEB2+(1.0-CFEB2)*(S-SE1)/(SE2-SE1))
            ELSE
              KTEBS(IK,IR) = TEB
            ENDIF
          ELSEIF (CIOPTK.EQ.2) THEN
            KTEBS(IK,IR) = ( TEBP**3.5 + 1.75*CPA*S/CK0 ) ** (2.0/7.0)
          ELSEIF (CIOPTK.EQ.3) THEN
            KTEBS(IK,IR) = ( TEBP**3.5 + 1.75*LPPA*S/CK0 ) ** (2.0/7.0)
          ELSEIF (CIOPTK.EQ.4) THEN
            IF (INPLATE) THEN
              KTEBS(IK,IR) = (TEBPI**3.5+1.75*LPPAI*S/CK0)**(2.0/7.0)
            ELSE
              KTEBS(IK,IR) = (TEBP**3.5 + 1.75*LPPA*S/CK0 )**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTK.EQ.5) THEN
            IF (INPLATE) THEN
              KTEBS(IK,IR) = (TEBPI**3.5+1.75*LPPAEI*S/CK0)**(2.0/7.0)
            ELSE
              KTEBS(IK,IR) = (TEBP**3.5+1.75*LPPAEB*S/CK0 )**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTK.EQ.6) THEN
            IF (INPLATE) THEN
              KTEBS(IK,IR) = (TEBPI**3.5+3.5*LPPAI*S/CK0)**(2.0/7.0)
            ELSE
              KTEBS(IK,IR) = (TEBP**3.5 +3.5*LPPA*S/CK0 )**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTK.EQ.7) THEN
            IF (INPLATE) THEN
              KTEBS(IK,IR) = (TEBPI**3.5+3.5*LPPAEI*S/CK0)**(2.0/7.0)
            ELSE
              KTEBS(IK,IR) = (TEBP**3.5+3.5*LPPAEB*S/CK0 )**(2.0/7.0)
            ENDIF
          ENDIF
C
C
          IF     (CIOPTL.EQ.0) THEN
            IF     (S.LE.0.0) THEN
              KTIBS(IK,IR) = TIB * CFIBT
            ELSEIF (S.LE.SI1) THEN
              KTIBS(IK,IR) = TIB * (CFIBT+(1.0-CFIBT)*S/SI1)
            ELSE
              KTIBS(IK,IR) = TIB
            ENDIF
          ELSEIF (CIOPTL.EQ.1) THEN
            IF     (S.LE.0.0) THEN
              KTIBS(IK,IR) = TIB * CFIBT
            ELSEIF (S.LE.SI1) THEN
              KTIBS(IK,IR) = TIB * (CFIBT+(CFIB2-CFIBT)*S/SI1)
            ELSEIF (S.LE.SI2) THEN
              KTIBS(IK,IR) = TIB * (CFIB2+(1.0-CFIB2)*(S-SI1)/(SI2-SI1))
            ELSE
              KTIBS(IK,IR) = TIB
            ENDIF
          ELSEIF (CIOPTL.EQ.2) THEN
            KTIBS(IK,IR) = ( TIBP**3.5 + 1.75*CPA*S/CK0 ) ** (2.0/7.0)
          ELSEIF (CIOPTL.EQ.3) THEN
            KTIBS(IK,IR) = ( TIBP**3.5 + 1.75*LPPA*S/CK0 ) ** (2.0/7.0)
          ELSEIF (CIOPTL.EQ.4) THEN
            IF (INPLATE) THEN
              KTIBS(IK,IR) = (TIBPI**3.5+1.75*LPPAI*S/CK0)**(2.0/7.0)
            ELSE
              KTIBS(IK,IR) = (TIBP**3.5 + 1.75*LPPA*S/CK0 )**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTL.EQ.5) THEN
            IF (INPLATE) THEN
              KTIBS(IK,IR) =(TIBPI**3.5+1.75*LPPAII*S/CK0I)**(2.0/7.0)
            ELSE
              KTIBS(IK,IR) = (TIBP**3.5+1.75*LPPAIB*S/CK0I)**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTL.EQ.6) THEN
            IF (INPLATE) THEN
              KTIBS(IK,IR) = (TIBPI**3.5+3.5*LPPAI*S/CK0)**(2.0/7.0)
            ELSE
              KTIBS(IK,IR) = (TIBP**3.5 +3.5*LPPA*S/CK0 )**(2.0/7.0)
            ENDIF
          ELSEIF (CIOPTL.EQ.7) THEN
            IF (INPLATE) THEN
              KTIBS(IK,IR) =(TIBPI**3.5+3.5*LPPAII*S/CK0I)**(2.0/7.0)
            ELSE
              KTIBS(IK,IR) = (TIBP**3.5+3.5*LPPAIB*S/CK0I)**(2.0/7.0)
            ENDIF
          ENDIF
C
  300   CONTINUE
  400  CONTINUE
C
 
      return
      end
c
c
c
      subroutine set_initplasma(ir,ikopt)
      implicit none
      integer ir,ikopt
c
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'cedge2d'
c
c     SET_INITPLASMA: This routine sets the initial plasma
c                     values on the specified ring depending 
c                     on the plasma decay option specified.
c
c
c
c     Local variables 
c
      integer ik,in,ierr,nr,ikstart,ikend
      INTEGER RINGNO
      EXTERNAL RINGNO
c
      call set_ikvals(ir,ikstart,ikend,ikopt)
c
c     Loop through ring assigning initial values to each cell 
c
        DO 100 IK = ikstart,ikend
C
          IF     (IR.LT.IRSEP) THEN
C
            if (cioptg.eq.7.or.ccoreopt.eq.5) then 

              IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                NR = RINGNO(IR,LPDATO,NLPDATO,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATO(NR,2)
                KTIBS(IK,IR) = LPDATO(NR,3)
                KNBS (IK,IR) = LPDATO(NR,4)
                kvhs (ik,ir) = 0.0
              ELSE
                NR = RINGNO(IR,LPDATI,NLPDATI,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATI(NR,2)
                KTIBS(IK,IR) = LPDATI(NR,3)
                KNBS (IK,IR) = LPDATI(NR,4)
                kvhs (ik,ir) = 0.0
              ENDIF
c
            else
               KTEBS(IK,IR) = CTEb0 + REAL(IRSEP-IR) * CTEBIN
               KTIBS(IK,IR) = CTIb0 + REAL(IRSEP-IR) * CTIBIN
               KNBS (IK,IR) = CNb0  + REAL(IRSEP-IR) * CNBIN
               kvhs(ik,ir) = 0.0
            endif
C
          ELSEIF (((cgridopt.eq.0.or.cgridopt.eq.1
     >       .or.cgridopt.eq.3).and.IR.LE.IRWALL)
     >      .or.(cgridopt.eq.2.and.((ir.ge.irsep.and.ir.le.irwall2)
     >          .or.(ir.ge.irsep2.and.ir.le.irwall)))) THEN
C
            if (e2dtargopt.eq.1.or.e2dtargopt.eq.3.or.
     >          e2dtargopt.eq.4.or.e2dtargopt.eq.6) then 
               IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                  KNBS (IK,IR) = e2dtarg(ir,1,2)
                  KTEBS(IK,IR) = e2dtarg(ir,2,2)
                  KTIBS(IK,IR) = e2dtarg(ir,3,2)
               else
                  KNBS (IK,IR) = e2dtarg(ir,1,1)
                  KTEBS(IK,IR) = e2dtarg(ir,2,1)
                  KTIBS(IK,IR) = e2dtarg(ir,3,1)
               endif
            elseif (e2dtargopt.eq.2.or.e2dtargopt.eq.5) then 
               IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                  KNBS (IK,IR) = e2dtarg(ir,7,2)
                  KTEBS(IK,IR) = e2dtarg(ir,2,2)
                  KTIBS(IK,IR) = e2dtarg(ir,3,2)
               else
                  KNBS (IK,IR) = e2dtarg(ir,7,1)
                  KTEBS(IK,IR) = e2dtarg(ir,2,1)
                  KTIBS(IK,IR) = e2dtarg(ir,3,1)
               endif
            elseif (CIOPTG.EQ.0) THEN
              KTEBS(IK,IR) = MAX (LO, CTEB0 - REAL(IR-IRSEP) * CTEBOU)
              KTIBS(IK,IR) = MAX (LO, CTIB0 - REAL(IR-IRSEP) * CTIBOU)
              KNBS (IK,IR) = MAX (1.E10, CNB0-REAL(IR-IRSEP) * CNBOUT)
            ELSEIF (CIOPTG.EQ.1) THEN
c              IF (IR.NE.IRSEP.AND.IK.EQ.1) DIST = DIST + KINDS(IKREF,IR)
              KTEBS(IK,IR) = CTEB0 * EXP (-refDIST(ir) / CTEBOU)
              KTIBS(IK,IR) = CTIB0 * EXP (-refDIST(ir) / CTIBOU)
              KNBS (IK,IR) = CNB0  * EXP (-refDIST(ir) / CNBOUT)
            ELSEIF (CIOPTG.EQ.2) THEN
              NR = RINGNO(IR,LPDATI,NLPDATI,MAXINS,4,ierr)
              KTEBS(IK,IR) = LPDATI(NR,2)
              KTIBS(IK,IR) = LPDATI(NR,3)
              if (lpdatsw.eq.0) then 
                 KNBS (IK,IR) = LPDATI(NR,4)
              elseif (lpdatsw.eq.1) then  
c
c                This initially assumes that the Mach number 
c                at the target = 1.0
c
                 KNBS (IK,IR) = abs(lpdati(nr,4) / (ech *
     >                 (9.79E+03 * 
     >                 SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
     >                 *(1.0+RIZB)/CRMB))))
              endif
            ELSEIF (CIOPTG.EQ.3.OR.CIOPTG.EQ.4.or.cioptg.eq.7) THEN
              IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                NR = RINGNO(IR,LPDATO,NLPDATO,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATO(NR,2)
                KTIBS(IK,IR) = LPDATO(NR,3)
                if (lpdatsw.eq.0) then 
                   KNBS (IK,IR) = LPDATO(NR,4)
                elseif (lpdatsw.eq.1) then 
c
c                  This initially assumes that the Mach number 
c                  at the target = 1.0
c
                   KNBS (IK,IR) = abs(lpdato(nr,4) / (ech *
     >                 (9.79E+03 * 
     >                 SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
     >                 *(1.0+RIZB)/CRMB))))
                endif
              ELSE
                NR = RINGNO(IR,LPDATI,NLPDATI,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATI(NR,2)
                KTIBS(IK,IR) = LPDATI(NR,3)
                if (lpdatsw.eq.0) then 
                   KNBS (IK,IR) = LPDATI(NR,4)
                elseif (lpdatsw.eq.1) then  
c
c                  This initially assumes that the Mach number 
c                  at the target = 1.0
c
                   KNBS (IK,IR) = abs(lpdati(nr,4) / (ech *
     >                   (9.79E+03 * 
     >                   SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
     >                   *(1.0+RIZB)/CRMB))))
                endif
              ENDIF
            ELSEIF (CIOPTG.EQ.5.or.cioptg.eq.6) THEN
c
              if (kss(ik,ir).lt.0.5*ksmaxs(ir)) then
c
c                Outer target
c
                 in = 2
              elseif (kss(ik,ir).ge.0.5*ksmaxs(ir)) then
c
c                Inner target
c
                 in = 1
              endif
c
              KTEBS(IK,IR) = CTEBP*EXP(-sepdist(idds(ir,in))/CTEBOU)
              KTIBS(IK,IR) = CTIBP*EXP(-sepdist(idds(ir,in))/CTIBOU)
              KNBS (IK,IR) = CNEBP*EXP(-sepdist(idds(ir,in))/CNBOUT)
c
            ENDIF
C
c         At this point in time any remaining regions are trap
c         regions and the test on the ELSEIF is really not
c         required.
c
          ELSEIF (IR.LE.NRS) THEN
C
            if (e2dtargopt.eq.1.or.e2dtargopt.eq.3.or.
     >          e2dtargopt.eq.4.or.e2dtargopt.eq.6) then 
               IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                  KNBS (IK,IR) = e2dtarg(ir,1,2)
                  KTEBS(IK,IR) = e2dtarg(ir,2,2)
                  KTIBS(IK,IR) = e2dtarg(ir,3,2)
               else
                  KNBS (IK,IR) = e2dtarg(ir,1,1)
                  KTEBS(IK,IR) = e2dtarg(ir,2,1)
                  KTIBS(IK,IR) = e2dtarg(ir,3,1)
               endif
            elseif (e2dtargopt.eq.2.or.e2dtargopt.eq.5) then 
               IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                  KNBS (IK,IR) = e2dtarg(ir,7,2)
                  KTEBS(IK,IR) = e2dtarg(ir,2,2)
                  KTIBS(IK,IR) = e2dtarg(ir,3,2)
               else
                  KNBS (IK,IR) = e2dtarg(ir,7,1)
                  KTEBS(IK,IR) = e2dtarg(ir,2,1)
                  KTIBS(IK,IR) = e2dtarg(ir,3,1)
               endif
            elseIF (CIOPTG.EQ.4.or.cioptg.eq.7) THEN
              IF (KSS(IK,IR).LE.0.5*KSMAXS(IR)) THEN
                NR = RINGNO(IR,LPDATO,NLPDATO,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATO(NR,2)
                KTIBS(IK,IR) = LPDATO(NR,3)
                if (lpdatsw.eq.0) then 
                   KNBS (IK,IR) = LPDATO(NR,4)
                elseif (lpdatsw.eq.1) then 
c
c                  This initially assumes that the Mach number 
c                  at the target = 1.0
c
                   KNBS (IK,IR) = abs(lpdato(nr,4) / (ech *
     >                 (9.79E+03 * 
     >                 SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
     >                 *(1.0+RIZB)/CRMB))))
                endif
              ELSE
                NR = RINGNO(IR,LPDATI,NLPDATI,MAXINS,4,ierr)
                KTEBS(IK,IR) = LPDATI(NR,2)
                KTIBS(IK,IR) = LPDATI(NR,3)
                if (lpdatsw.eq.0) then 
                   KNBS (IK,IR) = LPDATI(NR,4)
                elseif (lpdatsw.eq.1) then  
c
c                  This initially assumes that the Mach number 
c                  at the target = 1.0
c
                   KNBS (IK,IR) = abs(lpdati(nr,4) / (ech *
     >                   (9.79E+03 * 
     >                   SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
     >                   *(1.0+RIZB)/CRMB))))

c                   write(6,'(a,2i4,8g8.2)') 'N:',ik,ir,
c     >                   (9.79E+03 * 
c     >                   SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
c     >                   *(1.0+RIZB)/CRMB)),ktebs(ik,ir),
c     >                   ktibs(ik,ir),rizb,crmb,knbs(ik,ir)
c

                endif
              ENDIF
            elseif (cioptg.eq.6) then
c
              if (kss(ik,ir).lt.0.5*ksmaxs(ir)) then
c
c               Outer target
c
                in = 2
              elseif (kss(ik,ir).ge.0.5*ksmaxs(ir)) then
c
c               Inner target
c
                in = 1
              endif
c
              KTEBS(IK,IR) = CTEBP*EXP(-sepdist(idds(ir,in))/CTEBOUP)
              KTIBS(IK,IR) = CTIBP*EXP(-sepdist(idds(ir,in))/CTIBOUP)
              KNBS (IK,IR) = CNEBP*EXP(-sepdist(idds(ir,in))/CNBOUP)
c
            ELSE
              KTEBS(IK,IR) = CTEBT
              KTIBS(IK,IR) = CTIBT
              KNBS (IK,IR) = CNBT
            ENDIF
C
          ENDIF
C
  100   CONTINUE
C
C       ASSIGN INITIAL PLATE VALUES OF TEMPERATURE AND
C       DENSITY - THESE MAY BE OVERWRITTEN AT A LATER
C       POINT WITH CORRECTED VALUES.
c
c       Assume sonic initially ... 
C
        IF (IR.GE.IRSEP) THEN
C
C         OUTER PLATE
C 
          if (ikopt.eq.1.or.ikopt.eq.3) then  
c
             KTEDS(IDDS(IR,2)) = KTEBS(1,IR)
             KTIDS(IDDS(IR,2)) = KTIBS(1,IR)
             KNDS(IDDS(IR,2)) =  KNBS(1,IR)
c
             if (e2dtargopt.eq.1.or.e2dtargopt.eq.2.or.
     >           e2dtargopt.eq.3.or.e2dtargopt.eq.5) then 
c  
                kvds(idds(ir,2)) = e2dtarg(ir,4,2) 
c
             elseif (e2dtargopt.eq.4) then  
c
                kvds(idds(ir,2)) = e2dtarg(ir,8,2) 
c
             else
c     >                 (9.79E+03 * 
c     >                 SQRT(0.5*(KTEBs(ik,ir)+KTIBS(IK,ir))
c     >                 *(1.0+RIZB)/CRMB)))
c
c     >         SQRT(0.5*EMI*(KTEdS(idds(i,1))+KTIdS(idds(i,1)))
c     >         *(1+RIZB)/CRMB)

                kvds(idds(ir,2)) = 
     >          -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))/crmb 
     >          * emi)
             endif
c
          endif 

C
C         INNER PLATE
C
          if (ikopt.eq.2.or.ikopt.eq.3) then  
c
             KTEDS(IDDS(IR,1)) = KTEBS(NKS(IR),IR)
             KTIDS(IDDS(IR,1)) = KTIBS(NKS(IR),IR)
             KNDS(IDDS(IR,1)) =  KNBS(NKS(IR),IR)
c
             if (e2dtargopt.eq.1.or.e2dtargopt.eq.2.or.
     >           e2dtargopt.eq.3.or.e2dtargopt.eq.5) then 
c  
                kvds(idds(ir,1)) = e2dtarg(ir,4,1) 
c
             elseif (e2dtargopt.eq.4) then  
c
                kvds(idds(ir,1)) = e2dtarg(ir,8,1) 
c
             else

                kvds(idds(ir,1)) = 
     >          sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))/crmb 
     >          * emi)
c slmod begin
                IF (cvhout.EQ.99.0.AND.ir.EQ.irsep) 
     .            cvhout = -kvds(idds(ir,1))
c slmdo end                
c                write (0,'(a,i4,10g9.2)') 'T:',ir,
c     >            kvds(idds(ir,1)),kteds(idds(ir,1)),
c     >            ktids(idds(ir,1)),
c     >            ktebs(nks(ir),ir),ktibs(nks(ir),ir), 
c     >                 crmb,emi,rizb,
c     >                 (9.79E+03 * 
c     >                 SQRT(0.5*(KTEBs(nks(ir),ir)+KTIBS(nks(ir),ir))
c     >                 *(1.0+RIZB)/CRMB)),
c     >         SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
c     >         *(1+RIZB)/CRMB)
c
c                write (6,'(a,i4,10g9.2)') 'T:',ir,
c     >            kvds(idds(ir,1)),kteds(idds(ir,1)),
c     >            ktids(idds(ir,1)),
c     >            ktebs(nks(ir),ir),ktibs(nks(ir),ir), 
c     >                 crmb,emi,rizb,
c     >                 (9.79E+03 * 
c     >                 SQRT(0.5*(KTEBs(nks(ir),ir)+KTIBS(nks(ir),ir))
c     >                 *(1.0+RIZB)/CRMB)),
c     >         SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
c     >         *(1+RIZB)/CRMB)
c

c
             endif
c
          endif
c
          write (6,'(a,1x,2i4,4g14.6)') 
     >                    'TargO:',ir,idds(ir,2),kteds(idds(ir,2)),
     >                    ktids(idds(ir,2)),knds(idds(ir,2)),
     >                    kvds(idds(ir,2)) 
c
          write (6,'(a,1x,2i4,4g14.6)') 
     >                    'TargI:',ir,idds(ir,1),kteds(idds(ir,1)),
     >                    ktids(idds(ir,1)),knds(idds(ir,1)),
     >                    kvds(idds(ir,1)) 
c
        elseif (ir.lt.irsep) then 
c
c             Initialize last cell on core rings
c            
              KTEBS(nks(ir),IR) = ktebs(1,ir)
              KTIBS(nks(ir),IR) = ktibs(1,ir)
              KNBS (nks(ir),IR) = knbs(1,ir)
              KVHS (nks(ir),ir) = kvhs(1,ir)
c
        ENDIF
c


      return
      end 



