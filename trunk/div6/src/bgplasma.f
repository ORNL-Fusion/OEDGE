c     -*Fortran*-
c
      subroutine bgplasma(title,equil)
      use error_handling
      use debug_options
c slmod begin
      USE mod_sol28_global
c slmod end
      implicit none
C
      include 'params'
c
      include 'cgeom'
c
      include 'cedge2d'
c
      include 'comtor'
c
      include 'cioniz'
c
c      include 'reader'
c
      include 'dynam1'
      include 'dynam5'
c
      include 'pindata'
c slmod begin - new
      include 'slcom'
     
      INTEGER status
      LOGICAL callsol28,message_reverse

      DATA message_reverse /.TRUE./
c slmod end
c
      character*(*) title,equil
c
c     BGPLASMA: This subroutine utilizes the various
c     plasma decay, SOL and temperature gradient options
c     to specify or calculate a background plasma for
c     the simulation. This code was written to facilitate
c     the use of different options for different half-rings
c     and to allow EDGE2D solutions and specified solutions
c     to be combined in the same plasma background.
c
c     David Elder, July 26, 1997.
c
c
c     The following code divides into two sections
c
c     1) The original logic BG plasma code for each option
c     2) The newer BG plasma options allowing for combinations
c        of options to be specified for each ring.
c
c     The code that performs the basic calculations for the
c     plasma decay, SOL and temperature gradient options has
c     been parameterized to take the range of ring numbers and
c     which portion of each ring as arguments. The argument
c     that specifies which portion of each ring currently can
c     take three values -
c
c     1 - first half of ring - corresponds to OUTER for JET grids
c     2 - inner half of ring - corresponds to INNER for JET grids
c     3 - full ring
c
c     NOTE: The target designations themselves differ from the
c           1/2 ring designations. The first target - indexed by 1 -
c           is the INNER target, while the second target - indexed by 2-
c           is the OUTER target for JET grids.
c
c
c
c     Local variables
c
      integer ik,ir,retcode,iitersol,ir1,ir2,id
      logical lpinopt,lpinavail,litersol,liter
      real    pintim
c
c     Temporary variables - to hold original variable values
c
      integer tmpcioptg,tmpcioptf,tmpcsopt
c
      integer  ikscut1,ikscut2,sfind,ikopt,in,num
      real     scut1,scut2
      external sfind
c
      call pr_trace('BGPLASMA','START OF BGPLASMA')
c
c     Initialization
c
      lpinopt = .false.
      lpinavail = .false.
      litersol = .false.
      liter = .false.
c
c     jdemod - init chisq records
c
      chisq1 = 0.0
      chisq2 = 0.0
      chisq3 = 0.0
      chisq4 = 0.0
      chisq5 = 0.0
c
c     Zero PIN time accounting variables
c
      pintim = 0.0
      totpintim = 0.0
c slmod begin 

c...  Assign IKBOUNDS array:
      CALL SetBounds

      opt%cflukin = 0

c...  Determine if SOL28 is going to be called:
      callsol28 = .FALSE.
      if (cioptf.eq.28) callsol28 = .TRUE.
      if (cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) then
        do id = 1, nbgplas
          if (bgplasopt(id,5).eq.28.0.or.bgplasopt(id,8).eq.28.0)
     .      callsol28 = .TRUE.
        enddo
      endif
      if (callsol28.AND.s28mode.GE.4.0) THEN
        CALL SetupSOL28
c...    If these are set to zero, on account of SOL28 and the fact that it
c       doesn't use the standard target data arrays, then initialize:
        IF (ncoredat.EQ.0) THEN
          coredat = 0.0
          DO ir = 2, irsep
            ncoredat = ncoredat + 1
            coredat(ncoredat,1) = REAL(ir)
          ENDDO
        ENDIF
        IF (nlpdato.EQ.0) THEN
          lpdato = 0.0
          DO ir = irsep, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            nlpdato = nlpdato + 1
            lpdato(nlpdato,1) = REAL(ir)
          ENDDO
        ENDIF
        IF (nlpdati.EQ.0) THEN
          lpdati = 0.0
          DO ir = irsep, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            nlpdati = nlpdati + 1
            lpdati(nlpdati,1) = REAL(ir)
          ENDDO
        ENDIF
      ENDIF
c slmod end
c
c     Piece-meal SOL options
c

      if (cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) then


         call pr_trace('BGPLASMA','PIECEWISE SOL OPT START')

         write (6,*) 'PD90 and 91 Under development ...'
c         write (0,*) 'PD90 and 91 Under development ...'
c
C
C-----------------------------------------------------------------------
C     INITIALISE COUNTER FOR USE IF SOL IS TO BE CALCULATED ITERATIVELY
C-----------------------------------------------------------------------
C
      IF ((CPINOPT.EQ.1.OR.CPINOPT.EQ.4).or.
     .    (cneuta.eq.1.and.ciopte.eq.4)) THEN
        LPINOPT = .TRUE.
        IITERSOL = 1
        IITERPIN = 0
c slmod begin
c...    Avoid the 0th call to PIN if boundary condition relaxation is
c       being used:
        IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
          CALL SetupIteration(iitersol)
          CALL SaveSolution
          IITERSOL = 2
          IITERPIN = 1
        elseif (citersol.eq.2) then
          IF (sloutput) THEN
            WRITE(0,*) 
            WRITE(0,*) '*********************************************'
            WRITE(0,*) '*  WARNING: CITERSOL.EQ.2 logic has changed *'
            WRITE(0,*) '*********************************************'
            WRITE(0,*) 
          ENDIF
          IITERSOL = 1
          IITERPIN = 1
        ENDIF
c slmod end
      ENDIF
c
      if (citersol.eq.1.or.citersol.eq.2) then
         litersol = .true.
      endif

C
C-----------------------------------------------------------------------
C     OVERRIDE PLASMA DECAY AND SOL OPTIONS IF REQUIRED.
C-----------------------------------------------------------------------
C
  360 CONTINUE
c slmod begin - new
      opt%cflukin = opt%cflukin + 1
      CALL SetupIteration(iitersol)
c slmod end
c
c      solstart = 1
c
c
C-----------------------------------------------------------------------
c
c        SET up targets for base FLUID CODE plasma for option 90
c
c        If option 90 is in use the basis for the BG is
c        an EDGE2D solution.
c
c        If option 91 is in use then the basis for the BG
c        assumes plasma decay option 4 and uses all of the
c        other options normally.
c
c        After the background has been calculated as stated
c        above - it will then be overlaid in the specific
c        regions as requested in the BGPLASOPT input array.
c
c        The core is recalculated using the base option after
c        SOL modifications have been applied and then any core
c        ring modifications are applied.
c

      IF (cioptg.eq.90) THEN
         if (cgridopt.eq.0) then

c
c           For JET grids - calculate target quantities appropriately
c           based on Edge2d values.
c
            DO IR = IRSEP, NRS
c
              KTEDS(IDDS(IR,1)) = e2dtarg(ir,2,1)
              KTEDS(IDDS(IR,2)) = e2dtarg(ir,2,2)
              KTIDS(IDDS(IR,1)) = e2dtarg(ir,3,1)
              KTIDS(IDDS(IR,2)) = e2dtarg(ir,3,2)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
              if (e2dtargopt.eq.2.or.e2dtargopt.eq.5) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,7,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,7,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)
c
              elseif (e2dtargopt.eq.4) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,8,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,8,2)
c
              elseif (e2dtargopt.eq.6) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)

                 KVDS(IDDS(IR,1))  =
     >               sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))
     >                /crmb *  emi)

                 KVDS(IDDS(IR,2))  =
     >               sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))
     >                /crmb *  emi)
c
              else
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)
c
              endif
c
            end do
c
         elseif (cgridopt.eq.3) then
c
            DO IR = IRSEP, NRS
c
c             Use loaded fluid code target conditions from the 
c             B2REPL routine in TAU
c           
              KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
              KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
              KTEDS(IDDS(IR,1)) = e2dtarg(ir,2,1)
              KTEDS(IDDS(IR,2)) = e2dtarg(ir,2,2)
              KTIDS(IDDS(IR,1)) = e2dtarg(ir,3,1)
              KTIDS(IDDS(IR,2)) = e2dtarg(ir,3,2)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
              if (e2dtargopt.eq.6) then 
c
                 KVDS(IDDS(IR,1))  =
     >               sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))
     >                /crmb *  emi)

                 KVDS(IDDS(IR,2))  =
     >               sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))
     >                /crmb *  emi)
c
              else

                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)

              endif
c
c
c              KTEDS(IDDS(IR,1)) = KTEBS(NKS(IR),IR)
c              KTEDS(IDDS(IR,2)) = KTEBS(1      ,IR)
c              KTIDS(IDDS(IR,1)) = KTIBS(NKS(IR),IR)
c              KTIDS(IDDS(IR,2)) = KTIBS(1      ,IR)
c              KNDS(IDDS(IR,1))  = KNBS(NKS(IR),IR)
c              KNDS(IDDS(IR,2))  = KNBS(1      ,IR)
c
c             Set to sound speed for now -
c
c              KVDS(IDDS(IR,1)) = 9.79E3 *
c     >                 SQRT (0.5*(kteds(idds(ir,1))
c     >                           +ktids(idds(ir,1)))
c     >                       *(1.0+RIZB)/CRMB)
c
c              KVDS(IDDS(IR,2))  = 9.79E3 *
c     >                 SQRT (0.5*(kteds(idds(ir,2))
c     >                           +ktids(idds(ir,2)))
c     >                       *(1.0+RIZB)/CRMB)
c
c              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
c              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
            ENDDO
         else
            DO IR = IRSEP, NRS
              KTEDS(IDDS(IR,1)) = KTEBS(NKS(IR),IR)
              KTEDS(IDDS(IR,2)) = KTEBS(1      ,IR)
              KTIDS(IDDS(IR,1)) = KTIBS(NKS(IR),IR)
              KTIDS(IDDS(IR,2)) = KTIBS(1      ,IR)
              KNDS(IDDS(IR,1))  = KNBS(NKS(IR),IR)
              KNDS(IDDS(IR,2))  = KNBS(1      ,IR)
              KVDS(IDDS(IR,1))  = KVHS(NKS(IR),IR)
              KVDS(IDDS(IR,2))  = KVHS(1      ,IR)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
            ENDDO
         endif
c
c        Overwrite the private plasma if that option has
c        been selected.
c
c        If specified private plasma is ON - then it calls the
c        routine to calculate the private plasma values.
c
         if (ciopto.eq.2.or.ciopto.eq.3.or.ciopto.eq.4) then
c
            tmpcioptg = cioptg
            cioptg = 4
c
            write (6,*) 'Spec:',irtrap,nrs,ciopto
c
            do ir = irtrap,nrs
c
               call set_initplasma(ir,3)
c
            end do
c
            cioptg = tmpcioptg
c
            if (ciopto.eq.2) then

               call specplas(irtrap,nrs,3)
c
            elseif (ciopto.eq.3.or.ciopto.eq.4) then

               call thompp(irtrap,nrs,3,ciopto)

            endif
c
         endif
c
c        Print out background for reference
c
c         write (6,*)
c         write (6,*) 'Background Plasma:'
c         write (6,*)
c
c         do ir = 1,nrs
c            do ik = 1,nks(ir)
c               write(6,'(2i4,5g15.5)') ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                    knbs(ik,ir),kvhs(ik,ir)
c            end do
c         end do
c
c     Read in a background plasma format written by a previous
c     DIVIMP run - in a DIVIMP format.
c
      elseif (cioptg.eq.92) then
c
          call readdivbg
c
      elseif (cioptg.eq.91) then
c
c         Calculate the base plasma based on the set input quantities and
c         using CIOPTG=4
c
          tmpcioptg=cioptg
          cioptg = 4
c
          CALL INITPLASMA(1,nrs,3)
          CALL SOL_PLASMA(irsep,nrs,3)
          CALL CORE_PLASMA(1,irsep-1,3)
          CALL SOL(irsep,nrs,3)
c
          cioptg = tmpcioptg
c
      ENDIF

      call pr_trace('BGPLASMA','BASIC PLASMA COMPLETE')
c
C-----------------------------------------------------------------------
c     BGPLASOPT:
c
c     At this point the basic plasma background has been either
c     loaded or calculated. It is now time to start overwriting with
c     the various options outlined in the BGPLASOPT array.
c
c     The BGPLASOPT array contains the data in the following order
c
c     BGPLASOPT(N,1) = Ring number Start
c                ,2  = Ring number End     2 >= 1
c                ,3  = Section specifer 1=first 2=second 3=all
c                ,4  = Plasma Decay option to use
c                ,5  = Sol option for section
c                ,6  = Te gradient option for section
c                ,7  = Ti gradient option for section
c                ,8  = Core option for ring
c                ,9  = E-field option for ring
c
c
c     Loop through over-writing SOL and PP rings
c
c     Save original BG opts.
c
      call load_bgopts(-2)
c
c     Set ciopto = 1 so that PP rings will be properly affected
c     by all options.
c
c     Process SOL and PP rings first
c
      do in = 1,nbgplas
         num  = in
         ir1   = bgplasopt(num,1)
         ir2   = bgplasopt(num,2)
c
         if (ir1.ge.irsep.and.ir1.le.nrs.and.
     >       ir2.ge.irsep.and.ir2.le.nrs) then
c
            ikopt = bgplasopt(num,3)
c
            call load_bgopts(num)
c
c           Calculate the BG for this plasma entry
c
c slmod begin
            if (ir1.le.ir2) then 
              CALL INITPLASMA(ir1,ir2,ikopt)
              CALL SOL_PLASMA(ir1,ir2,ikopt)
              CALL SOL(ir1,ir2,ikopt)
            else
c              Process the rings individually, in reverse order.  This is required
c              for SOL28 for PFZ rings that reference the neighbouring ring that
c              is closer to the separatrix.  In order for the neighbouring ring to
c              be defined (have a plasma solution already), eventhough it has a 
c              higher ring index, it needs to be sent to the solver first, hence 
c              the need to call the rings in reverse order. -SL, 04/04/2012
               if (message_reverse) then
                 CALL WN('bgplasma','Solving some rings in '//
     .                   'reverse order')
                 message_reverse = .FALSE.
               endif
               do ir = ir1, ir2, -1
                  CALL INITPLASMA(ir,ir,ikopt)
                  CALL SOL_PLASMA(ir,ir,ikopt)
                  CALL SOL(ir,ir,ikopt)
               enddo
            endif
c
c            CALL INITPLASMA(ir1,ir2,ikopt)
c            CALL SOL_PLASMA(ir1,ir2,ikopt)
c            CALL SOL(ir1,ir2,ikopt)
c slmod end
c
         endif
c
      end do
c
c     Reload original options
c
      call load_bgopts(-1)
c
c     Recalculate the core plasma based on the original options
c     factoring in any changes to the SOL - only do this if cioptg=91
c     was selected - option 90 will already have a valid core plasma.
c
      if (cioptg.eq.91) CALL CORE_PLASMA(1,irsep-1,3)
c
c     Run through the BGPLSOPT array and process any core rings listed.
c
      do in = 1,nbgplas
         num   = in
         ir1   = bgplasopt(num,1)
         ir2   = bgplasopt(num,2)
c
         if (ir1.ge.1.and.ir1.lt.irsep.and.
     >       ir2.ge.1.and.ir2.lt.irsep) then
c
            ikopt = bgplasopt(num,3)
c
            call load_bgopts(num)
c
c           Calculate the BG for this plasma entry
c
            CALL CORE_PLASMA(ir1,ir2,ikopt)
c
         endif
c
      end do
c
c     This should have patched all of the changes into the background
c     plasma - now add other change options and iterate through PIN
c     if necessary.
c
      call pr_trace('BGPLASMA','PLASMA PATCHING COMPLETE')
c
c     Reload original options
c
      call load_bgopts(-1)
c
C-----------------------------------------------------------------------
c     If Te and Ti over-rides are in effect - then impose flat Te and Ti
c     for S values greater than the specified Cutoff*SMAX.
c
c     Applied only to main SOL rings.
c
c     Te is flat for S > Fact * SMAX for each half ring
c
      call pr_trace('BGPLASMA','BEFORE FLATTEN T PROFILES')
c
      call flat_t

c
c     Overlay a separate interpolated plasma file if one is provided. This is ONLY ne, Te and will not be consistent
c     with any other plasma options unless some further development is done. This option over-writes the ne and Te by
c     reading in a file and interpolating for grid cells that lie within the region covered by the file. 
c
      if (external_plasma_overlay.eq.1) then 
         call overlay_plasma
      endif
c
C-----------------------------------------------------------------------
c     If the OFIELD option has been set to 1 ... over-ride the
c     calculated EFIELD and replace it by zeroes.
c
      call pr_trace('BGPLASMA','CALCULATE EFIELD')
c
      if (ofield.eq.0) then
c
c        Calculate E-field for B2 Background case - in AsdexU geometry
c
         if ((cgridopt.eq.3.or.cgridopt.eq.4.or.cgridopt.eq.5)
     >       .and.(cioptf.eq.99.or.cioptg.eq.99.or.cioptg.eq.90)) then
            call calcefb2
         endif
c
      elseif (ofield.eq.1) then
         call rzero(kes,maxnks*maxnrs)
c slmod end
c...  Added OFIELD=5 which is the same as 3 but also calculates the e-field 
c     in the PFZ:
      elseif (ofield.eq.3.or.ofield.eq.5) then
c
c      elseif (ofield.eq.3) then
c slmod begin
c
c        Use calculated BG electric field to overwrite anything read in.
c
         call calcefb2
c
      elseif (ofield.eq.4) then
c
c        Use calculated BG electric field to overwrite anything read in.
c
         call calcef4(lpinavail)
c
      endif
c
C-----------------------------------------------------------------------
c
C     MULTIPLY BY ENHANCEMENT FACTORS (TIME-STEP DONE LATER)
C
c     Properly normalize the velocity and electric field.      
c
      DO  IR = 1,NRS
        DO IK = 1,NKS(IR)
          KVHS(IK,IR) = KVHS(IK,IR) * CSOLVF
          KES(IK,IR)  = KES(IK,IR)  * CSOLEF
        end do
      end do
c
C-----------------------------------------------------------------------
c
c     Assign value of Ne
c
c     Assign the effective electron density for the case based
c     on the defined input option. (NE_OPT).
c
      if (ne_opt.eq.0) then
c
c        ne = ni
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               knes(ik,ir) = knbs(ik,ir)
            end do
         end do
      elseif (ne_opt.eq.1) then
c
c        ne = ni + Sigma (iz * nz)  (from fluid code result)
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               knes(ik,ir) = e2dnes(ik,ir)
            end do
         end do
      endif
c
c     Mirror target options 3 or 4 - zero the recorded target density
c     on these target elements in order to create a zero flux value 
c     calculated by the code for these elements.
c
      if (cmiropt.eq.3) then 
         do id = 1,ndsin
            knds(id) = 0.0
         enddo
      elseif (cmiropt.eq.4) then 
         do id = ndsin+1,nds
            knds(id) = 0.0
         enddo
      endif
c
c
c slmod begin 
      if (override_bg_velocity_opt.eq.1.and.osmns28.gt.0) then 
        CALL PrescribeFlow 
      endif

c...  Generate some standard analysis on the state of the solution:
      CALL OutputAnalysis
c slmod end
C
C-----------------------------------------------------------------------
c
C     INSERT CALL THAT WILL EXECUTE PIN :- IF IT IS NECESSARY TO ITERATE
C     SOL OPTION - FLOW WILL THEN BRANCH TO RE-EXECUTE THE PLASMA AND
C     SOL SUBROUTINES.
C

      call pr_trace('BGPLASMA','BEFORE EXECUTE PIN')
c
c     Execute the HYDROGENIC neutral code  PIN/NIMBUS or EIRENE
c
      call PINEXE(title,equil,lpinopt,litersol,iitersol,liter,
     >            tmpcsopt,tmpcioptf,lpinavail,iiterpin)
c
c     Iterate after PIN call if required
c
c slmod begin
      if (liter) then
        goto 360
      elseif (cpinopt.ne.0.and.citersol.eq.2.and.lpinopt) then
c      elseif (cpinopt.ne.0.and.citersol.eq.2.and.lpinopt.and.
c     .        rel_opt.ne.0) then
c...    Calculate plasma solution one more time, so that it is consistent
c       with the current PIN source terms:
        lpinopt = .false.
        goto 360
      endif
c
c      if (liter) goto 360
c slmod end
c
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c
c     REGULAR SOL OPTIONS ...
c
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
c

      else
c
         call pr_trace('BGPLASMA','STANDARD SOL OPTIONS')
c
         write (6,*) 'STANDARD SOL OPTIONS:'
C
C-----------------------------------------------------------------------
C     INITIALISE COUNTER FOR USE IF SOL IS TO BE CALCULATED ITERATIVELY
C-----------------------------------------------------------------------
C
      IF ((CPINOPT.EQ.1.OR.CPINOPT.EQ.4).or.
     .    (cneuta.eq.1.and.ciopte.eq.4)) THEN
        LPINOPT = .TRUE.
        IITERSOL = 1
        IITERPIN = 0
c slmod begin
c...    Avoid the 0th call to PIN if boundary condition relaxation is 
c       being used:
        IF (rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
          CALL SetupIteration(iitersol)
c          CALL SaveSolution
          IITERSOL = 2
          IITERPIN = 1
        elseif (citersol.eq.2) then
          if (sloutput) then
            WRITE(0,*) 
            WRITE(0,*) '*********************************************'
            WRITE(0,*) '*  WARNING: CITERSOL.EQ.2 logic has changed *'
            WRITE(0,*) '*********************************************'
            WRITE(0,*)
          endif
          IITERSOL = 1
          IITERPIN = 1
        ENDIF 
c slmod end
      ENDIF
c
      if (citersol.eq.1.or.citersol.eq.2) then
         litersol = .true.
      endif


  361 CONTINUE
c slmod begin - new
      opt%cflukin = opt%cflukin + 1
      CALL SetupIteration(iitersol)
c slmod end
c
c
C
C-----------------------------------------------------------------------
C     IF EDGE PLASMA IS TO BE USED, MAP NEAREST PLASMA POINT ONTO TARGET
C-----------------------------------------------------------------------
C
      IF (CIOPTG.EQ.99.or.cioptg.eq.90) THEN

         call pr_trace('BGPLASMA','LOAD EXTERNAL PLASMA')
         
         if (cgridopt.eq.0) then
c
c           For JET grids - calculate target quantities appropriately
c           based on Edge2d values.
c
            DO IR = IRSEP, NRS
c
              KTEDS(IDDS(IR,1)) = e2dtarg(ir,2,1)
              KTEDS(IDDS(IR,2)) = e2dtarg(ir,2,2)
              KTIDS(IDDS(IR,1)) = e2dtarg(ir,3,1)
              KTIDS(IDDS(IR,2)) = e2dtarg(ir,3,2)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
              if (e2dtargopt.eq.2.or.e2dtargopt.eq.5) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,7,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,7,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)
c
              elseif (e2dtargopt.eq.4) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,8,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,8,2)
c
              elseif (e2dtargopt.eq.6) then
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)

                 KVDS(IDDS(IR,1))  =
     >               sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))
     >                /crmb *  emi)

                 KVDS(IDDS(IR,2))  =
     >               sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))
     >                /crmb *  emi)
c
c
              else
c
                 KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
                 KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)
c
              endif
c
           end do
c
         elseif (cgridopt.eq.3) then
            DO IR = IRSEP, NRS
c
c             Use loaded fluid code target conditions from the 
c             B2REPL routine in TAU
c           
              KNDS(IDDS(IR,1))  = e2dtarg(ir,1,1)
              KNDS(IDDS(IR,2))  = e2dtarg(ir,1,2)
              KTEDS(IDDS(IR,1)) = e2dtarg(ir,2,1)
              KTEDS(IDDS(IR,2)) = e2dtarg(ir,2,2)
              KTIDS(IDDS(IR,1)) = e2dtarg(ir,3,1)
              KTIDS(IDDS(IR,2)) = e2dtarg(ir,3,2)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
              if (e2dtargopt.eq.6) then 
c
                 KVDS(IDDS(IR,1))  =
     >               sqrt((kteds(idds(ir,1))+ktids(idds(ir,1)))
     >                /crmb *  emi)

                 KVDS(IDDS(IR,2))  =
     >              -sqrt((kteds(idds(ir,2))+ktids(idds(ir,2)))
     >                /crmb *  emi)
c
              else

                 KVDS(IDDS(IR,1))  = e2dtarg(ir,4,1)
c slmod begin
c...             I had a problem with a plasma file that didn't have
c                the correct sign on the target velocities. This is not a 
c                general problem, but I put this check in to be sure:
                 IF (e2dtarg(ir,4,2).GT.0.0) THEN
                   WRITE(0,*) 'WARNING: LOW INDEX TARGET FLUID '//
     .                        'VELOCITY +VE, CHANGING SIGN',ir
                   KVDS(IDDS(IR,2)) = -e2dtarg(ir,4,2)
                 ELSE
                   KVDS(IDDS(IR,2)) = e2dtarg(ir,4,2)
                 ENDIF
c
c                 KVDS(IDDS(IR,2))  = e2dtarg(ir,4,2)
c slmod end

              endif
c
c
c              KTEDS(IDDS(IR,1)) = KTEBS(NKS(IR),IR)
c              KTEDS(IDDS(IR,2)) = KTEBS(1      ,IR)
c              KTIDS(IDDS(IR,1)) = KTIBS(NKS(IR),IR)
c              KTIDS(IDDS(IR,2)) = KTIBS(1      ,IR)
c              KNDS(IDDS(IR,1))  = KNBS(NKS(IR),IR)
c              KNDS(IDDS(IR,2))  = KNBS(1      ,IR)
c
c             Set to sound speed for now -
c
c              KVDS(IDDS(IR,1)) = 9.79E3 *
c     >                 SQRT (0.5*(kteds(idds(ir,1))
c     >                           +ktids(idds(ir,1)))
c     >                       *(1.0+RIZB)/CRMB)
c
c              KVDS(IDDS(IR,2))  = 9.79E3 *
c     >                 SQRT (0.5*(kteds(idds(ir,2))
c     >                           +ktids(idds(ir,2)))
c     >                       *(1.0+RIZB)/CRMB)
c
c              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
c              KEDS(IDDS(IR,2))  = KES(1      ,IR)
c
            ENDDO
         else
            DO IR = IRSEP, NRS
              KTEDS(IDDS(IR,1)) = KTEBS(NKS(IR),IR)
              KTEDS(IDDS(IR,2)) = KTEBS(1      ,IR)
              KTIDS(IDDS(IR,1)) = KTIBS(NKS(IR),IR)
              KTIDS(IDDS(IR,2)) = KTIBS(1      ,IR)
              KNDS(IDDS(IR,1))  = KNBS(NKS(IR),IR)
              KNDS(IDDS(IR,2))  = KNBS(1      ,IR)
              KVDS(IDDS(IR,1))  = KVHS(NKS(IR),IR)
              KVDS(IDDS(IR,2))  = KVHS(1      ,IR)
              KEDS(IDDS(IR,1))  = KES(NKS(IR),IR)
              KEDS(IDDS(IR,2))  = KES(1      ,IR)
            ENDDO
         endif
c
c
c        If specified private plasma is ON - then it calls the
c        routine to calculate the private plasma values.
c
         if (ciopto.eq.2.or.ciopto.eq.3.or.ciopto.eq.4) then
c
            tmpcioptg = cioptg
            cioptg = 4
c
            write (6,*) 'Spec:',irtrap,nrs,ciopto
c
            do ir = irtrap,nrs
c
               call set_initplasma(ir,3)
c
            end do
c
            cioptg = tmpcioptg
c
            if (ciopto.eq.2) then
c
               call specplas(irtrap,nrs,3)

            elseif (ciopto.eq.3.or.ciopto.eq.4) then

               call thompp(irtrap,nrs,3,ciopto)

            endif
c
         endif
c
c        If ccoreopt in greater than zero - overwrite the core plasma supplied 
c        by the fluid code
c
         if (ccoreopt.gt.0) then  

            CALL CORE_PLASMA(1,irsep-1,3)
         endif           
c
c        Print out background for reference
c
c         write (6,*)
c         write (6,*) 'Background Plasma:'
c         write (6,*)
c
c         do ir = 1,nrs
c            do ik = 1,nks(ir)
c               write(6,'(2i4,5g15.5)') ik,ir,ktebs(ik,ir),ktibs(ik,ir),
c     >                    knbs(ik,ir),kvhs(ik,ir)
c            end do
c         end do
c
c     Read in a background plasma format written by a previous
c     DIVIMP run - in a DIVIMP format.
c
      elseif (cioptg.eq.98) then
c
         call readdivbg
c
      ENDIF


C
C-----------------------------------------------------------------------
C     OVERRIDE PLASMA DECAY AND SOL OPTIONS IF REQUIRED.
C-----------------------------------------------------------------------
C
c
c     Generate background plasma
c
      IF (CIOPTG.NE.99.and.cioptg.ne.98) THEN
         CALL INITPLASMA(1,nrs,3)
         CALL SOL_PLASMA(irsep,nrs,3)
         CALL CORE_PLASMA(1,irsep-1,3)
      endif

      call pr_trace('BGPLASMA','BASIC PLASMA LOAED')
      
c
      IF (CIOPTF.NE.99.and.cioptf.ne.98) CALL SOL(irsep,nrs,3)
c
C-----------------------------------------------------------------------
c
c     If Te and Ti over-rides are in effect - then impose flat Te and Ti
c     for S values greater than the specified Cutoff*SMAX.
c
c     Applied only to main SOL rings.
c
c     Te is flat for S > Fact * SMAX for each half ring
c
      call pr_trace('BGPLASMA','FLATTEN T PROFILES 2')
      call flat_t

c
c     Overlay a separate interpolated plasma file if one is provided. This is ONLY ne, Te and will not be consistent
c     with any other plasma options unless some further development is done. This option over-writes the ne and Te by
c     reading in a file and interpolating for grid cells that lie within the region covered by the file. 
c
      if (external_plasma_overlay.eq.1) then 
         call overlay_plasma
      endif
c
C-----------------------------------------------------------------------
c
c     If the OFIELD option has been set to 1 ... over-ride the
c     calculated EFIELD and replace it by zeroes.
c
      call pr_trace('BGPLASMA','CALCULATE EFIELD 2')

      if (ofield.eq.0) then
c
c        Calculate E-field for B2 Background case - in AsdexU geometry
c
         if ((cgridopt.eq.3.or.cgridopt.eq.4.or.cgridopt.eq.5)
     >       .and.(cioptf.eq.99.or.cioptg.eq.99.or.cioptg.eq.90)) then
            call calcefb2
         endif
c
      elseif (ofield.eq.1) then
         call rzero(kes,maxnks*maxnrs)
      elseif (ofield.eq.3) then
c
c        Use calculated BG electric field to overwrite anything read in.
c
         call calcefb2
c
      elseif (ofield.eq.4) then
c
c        Use calculated BG electric field to overwrite anything read in.
c
         call calcef4(lpinavail)
c
      endif
c
      call pr_trace('BGPLASMA','BEFORE SCALING VB,EF')
C-----------------------------------------------------------------------
c
C     MULTIPLY BY ENHANCEMENT FACTORS (TIME-STEP DONE LATER)
C
c     Properly normalize the velocity and electric field.
c
      DO  IR = 1,NRS
        DO IK = 1,NKS(IR)
          KVHS(IK,IR) = KVHS(IK,IR) * CSOLVF
          KES(IK,IR)  = KES(IK,IR)  * CSOLEF
c
c         Correct for velocity orientation problem from
c         read-in background plasma for ITER B2-EIRENE runs
c
          if (cgridopt.eq.2.and.cioptf.eq.99) then
             if ((ir.ge.irsep.and.ir.le.irwall2).and.
     >           ik.gt.(nks(ir)/2)) then
                kvhs(ik,ir) = -kvhs(ik,ir)
             elseif ((ir.ge.irsep2.and.ir.le.irwall).and.
     >           ik.le.(nks(ir)/2)) then
                kvhs(ik,ir) = -kvhs(ik,ir)
             endif
          endif
c
        end do
      end do
c
C-----------------------------------------------------------------------
c
c     Assign value of Ne
c
c     Assign the effective electron density for the case based
c     on the defined input option. (NE_OPT).
c
      if (ne_opt.eq.0) then
c
c        ne = ni
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               knes(ik,ir) = knbs(ik,ir)
            end do
         end do
      elseif (ne_opt.eq.1) then
c
c        ne = ni + Sigma (iz * nz)  (from fluid code result)
c
         do ir = 1,nrs
            do ik = 1,nks(ir)
               knes(ik,ir) = e2dnes(ik,ir)
            end do
         end do
      endif
c
c     Mirror target options 3 or 4 - zero the recorded target density
c     on these target elements in order to create a zero flux value 
c     calculated by the code for these elements.
c
      if (cmiropt.eq.3) then 
         do id = 1,ndsin
            knds(id) = 0.0 
         enddo
      elseif (cmiropt.eq.4) then 
         do id = ndsin+1,nds
            knds(id) = 0.0 
         enddo
      endif
c
      call pr_trace('BGPLASMA','BEFORE PrescribeFlow')
c
c
c slmod begin 
      if (override_bg_velocity_opt.eq.1.and.osmns28.gt.0) then 
        CALL PrescribeFlow 
      endif


c...  Generate some standard analysis on the state of the solution:
      CALL OutputAnalysis

      IF (rel_opt.EQ.1) THEN
        CALL LoadObjects('osm_geometry.raw',status)
        IF (status.NE.0) 
     .    CALL ER('bgplasma','Unable to load geometry data',*99)
        CALL GenerateOutputFiles(iitersol-1)
        CALL geoClean
      ENDIF
c slmod end
C
C-----------------------------------------------------------------------
c
C     INSERT CALL THAT WILL EXECUTE PIN :- IF IT IS NECESSARY TO ITERATE
C     SOL OPTION - FLOW WILL THEN BRANCH TO RE-EXECUTE THE PLASMA AND
C     SOL SUBROUTINES.
C
      call pr_trace('BGPLASMA','EXECUTE PIN 2')
c
c     Execute the HYDROGENIC neutral code  PIN/NIMBUS or EIRENE
c
      call PINEXE(title,equil,lpinopt,litersol,iitersol,liter,
     >            tmpcsopt,tmpcioptf,lpinavail,iiterpin)
c
c     Iterate after PIN call if required
c
c slmod begin
      if (liter) then
        goto 361
      elseif (cpinopt.ne.0.and.citersol.eq.2.and.lpinopt) then
c      elseif (cpinopt.ne.0.and.citersol.eq.2.and.lpinopt.and.
c     .        rel_opt.ne.0) then
c...    Calculate plasma solution one more time, so that it is consistent
c       with the current PIN source terms:
        lpinopt = .false.
        goto 361
      endif
c
c      if (liter) goto 361
c slmod end
c
c     Endif for CIOPTG = 90 or 91
c
      endif
c
C-----------------------------------------------------------------------
c
c     If the Velocity Override option has been specified - impose or
c     calculate a new background velocity array. 
c
C-----------------------------------------------------------------------
c
c     Specified flow field 
c
      if (override_bg_velocity_opt.eq.1.and.osmns28.gt.0) then 

c slmod begin 
c        Adding some capacity to over-ride the SOL21 calcualted hydrogenic flow
c        velocity.  This should be moved in the future, preferentially to some
c        generic "over-ride" routine that is called after the solution has been
c        calculated, irrespective of the solver that is in use: 
c
           WRITE(0,*) 'CALLING vb OVER-RIDE CODE'

           CALL PrescribeFlow 
c slmod end 

c
c     Recalculate velocity field using EIRENE input - only
c     valid if EIRENE has been run - half ring
c
      elseif (override_bg_velocity_opt.eq.2.and.lpinopt) then 
c
         call recalculate_bg_velocity(1)
c
c
c     Recalculate velocity field using EIRENE input - only
c     valid if EIRENE has been run - full ring
c
      elseif (override_bg_velocity_opt.eq.3.and.lpinopt) then 
c
         call recalculate_bg_velocity(2)
c
c     Force the background flow velocity to zero everywhere
c
      elseif (override_bg_velocity_opt.eq.4) then 
c
         call set_bg_velocity(0)
c
      elseif (override_bg_velocity_opt.eq.5) then 
c
         call set_bg_velocity(1)
c
      elseif (override_bg_velocity_opt.eq.6) then 
c
         call set_bg_velocity(2)
c
      elseif (override_bg_velocity_opt.eq.7) then 
c
         call set_bg_velocity(3)
c
      elseif (override_bg_velocity_opt.eq.8) then 
c
         call set_bg_velocity(4)
c
      elseif (override_bg_velocity_opt.eq.9) then 
c
         call set_bg_velocity(5)
c
      elseif (override_bg_velocity_opt.eq.10) then 
c
         call set_bg_velocity(6)
c
      elseif (override_bg_velocity_opt.eq.11) then 
c
         call set_bg_velocity(7)
c
      endif
c
C-----------------------------------------------------------------------
c
c     Wrap up processing after the BG plasma has been calculated
c
C-----------------------------------------------------------------------
c
      call pr_trace('BGPLASMA','FINISH PROCESSING BGP')
      
c
c     IF Ofield = 4 and PIN has been run - call the E-field
c     calculation routine again to revise the values based on the
c     PIN generated Leq values.
c
      if ((cpinopt.eq.1.or.cpinopt.eq.4).and.ofield.eq.4) then
c
         call calcef4(lpinavail)
c
         DO  IR = 1,NRS
           DO IK = 1,NKS(IR) 
              KES(IK,IR)  = KES(IK,IR)  * CSOLEF
           end do
         end do
c
      endif
c
c     Reset the SOL and iteration ionization options
c
      if (piniter) then
         cioptf  = tmpcioptf
         csopt   = tmpcsopt
      endif
c
C     Calculate the upstream conditions - both inner and outer
C
      do ir = irsep,nrs
         teupstream(ir,2) = ktebs(ikmids(ir),ir)
         teupstream(ir,1) = ktebs(ikmids(ir)+1,ir)
c
         tiupstream(ir,2) = ktibs(ikmids(ir),ir)
         tiupstream(ir,1) = ktibs(ikmids(ir)+1,ir)
c
         nbupstream(ir,2) = knbs(ikmids(ir),ir)
         nbupstream(ir,1) = knbs(ikmids(ir)+1,ir)
c
      end do

c
c     End of BG plasma calculation - finish and clean up is in TAU
c
c slmod begin - new
c...Add print option:
      CALL AnalyseSolution(PINOUT)

c...TEMP:
      WRITE(67,'(A,I4,A)') 'ITERATION ',rel_count,' (after update)'
      WRITE(67,*)
      CALL OutputEIRENE(67,'DONE CALCULATING PLASMA SOLUTION')
      WRITE(67,*)

      IF (cpinopt.EQ.0.OR.citersol.EQ.0.OR.rel_opt.EQ.0.OR.
     .    citersol.EQ.2)
     .  CALL SaveSolution

      IF (callsol28.AND.s28mode.GE.4.0) CALL CloseSOL28

c...  This is needed, i.e. KNDS is NANQ for some cases:
      IF (.NOT.nopriv.AND.cgridopt.NE.LINEAR_GRID.AND.
     .                    cgridopt.NE.RIBBON_GRID) THEN
        knds(idds(irwall,1:2)) = 0.0
        knds(idds(irtrap,1:2)) = 0.0
      ENDIF

      call pr_trace('BGPLASMA','END OF BGPLASMA')


c      CALL SaveSolution
c slmod end
      return
c slmod begin
99    stop
c slmod end
      end
c     
c     
c     
      subroutine calcef4(lpinavail)
      implicit none
      logical lpinavail
      include 'params'
      include 'cgeom'
      include 'comtor'
c     
c     CALCEF4:
c     
c     This subroutine uses a simple prescription to estimate
c     the electric field from the background density and
c     temperature and the distances between grid points.
c     
c     This routine is modelled on the calcefb2 routine and
c     implements (OFIELD option = 4) where the collisionality of
c     the ring is checked before applying different E-field
c     models.
c     
c     
c     
      real ds1,dp1,dt1,nb1,ds2,dp2,dt2,nb2,irlimit
      real ef1,efnks,dist1,distnks,tav,smax
      integer ik,ir
C     
C     CALCULATE ELECTRIC FIELD
C     
C     IN THE FOLLOWING EQUATIONS THE FACTOR E CANCELS WITH THE
C     SAME FACTOR USED IN CONVERTING T IN EV TO KT.
C     
      call rzero(kes,maxnks*maxnrs)
c     
c     Calculate for Private Plasma only if the entire SOL is being
c     calculatd and not specified for the PP.
c     
      if (ciopto.eq.0.or.ciopto.eq.2.or.
     >     ciopto.eq.3.or.ciopto.eq.4) then
         irlimit = irwall
      elseif (ciopto.eq.1) then
         irlimit = nrs
      endif
c     
c     Loop through rings
c     
      do 600 ir = irsep,irlimit
c     
         smax = ksmaxs(ir)

c     
c     The DIST1 and DISTNKS contain the length of the source region
c     if the half ring is found to be non-collisional. This region
c     is the distance over which the electric field will be non-zero.
c     dist1 is for the first half of the ring with IK =1 (usually the
c     OUTER half on JET grids) and distnks is for the IK =NKS(IR) end
c     of the ring. DIST = -1.0 is used to identify a collisional ring
c     for which the other formula should be used.
c     
c     The E-field is -ve for toward the target at the IK=1 end
c     and +ve for toward the target at the IK=NKS(IR) end.
c     
c     
c     IK = 1
c     
c     write (6,*) 'Debug E-field option 4:',ir,lpinavail
c     
         if (teupstream(ir,2).lt.ceffact*kteds(idds(ir,2))) then
            tav = 0.5 * (teupstream(ir,2) +kteds(idds(ir,2)))
c     
            if (lpinavail.and.cleq(ir,2).ne.0.0) then
               dist1 = cleq(ir,2)
               ef1 = - tav / (2.0*dist1)
            elseif (ceflen.ne.0.0) then
               dist1 =  ceflen * smax
               ef1 = -tav / (2.0 * dist1)
            else
               dist1 = 0.0
               ef1 = 0.0
            endif
c     
         else
            dist1 = -1.0
            ef1 =   0.0
         endif
c     
c     write (6,'(a,6(g13.5,2x))') 'IK=1  ',
c     >           teupstream(ir,2),kteds(idds(ir,2)),cleq(ir,2),
c     >           tav,ef1,dist1
c     
c     IK = NKS(IR)
c     
         if (teupstream(ir,1).lt.ceffact*kteds(idds(ir,1))) then
            tav = 0.5 * (teupstream(ir,1) +kteds(idds(ir,1)))
c     
            if (lpinavail.and.cleq(ir,1).ne.0.0) then
               distnks = cleq(ir,1)
               efnks =   tav / (2.0*distnks)
            elseif (ceflen.ne.0.0) then
               distnks =  ceflen * smax
               efnks =  tav / (2.0 * distnks)
            else
               distnks = 0.0
               efnks = 0.0
            endif
c     
         else
            distnks = -1.0
            efnks = 0.0
         endif
c     
c     write (6,'(a,6(g13.5,2x))') 'IK=NKS',
c     >     teupstream(ir,1),kteds(idds(ir,1)),cleq(ir,1),
c     >     tav,efnks,distnks
c     
c     
c     Use the values calculated above in the evaluation for each point
c     
c     Start with end-points and then loop through the ring.
c     
         if (dist1.ge.0.0) then
            if (kss(1,ir).le.dist1) then
               kes(1,ir) = ef1
            else
               kes(1,ir) = 0.0
            endif
         else
            DS1 = KSS(2,IR) - KSS(1,IR)
            DP1 = (KNBS(2,IR)*KTEBS(2,IR)-KNBS(1,IR)*KTEBS(1,IR))
            DT1 = (KTEBS(2,IR)-KTEBS(1,IR))
            NB1 = 0.5*(KNBS(2,IR)+KNBS(1,IR))
c     
c     jdemod - modified to include and set target data
c     
c     KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1

            DS2 = KSS(1,IR) - KSb(0,IR)
            DP2 = (KNBS(1,IR)*KTEBS(1,IR)
     >            -  KNdS(idds(IR,2))*KTEdS(idds(IR,2)))
            DT2 = (KTEBS(1,IR)-KTEdS(idds(ir,2)))
            NB2 = 0.5*(KNBS(1,IR)+KNdS(idds(ir,2)))
c     

            if ((ds1 .ne. 0.0) .and. (nb1.ne.0.0).and.
     >           (ds2 .ne. 0.0) .and. (nb2.ne.0.0)) then 

               if (ofield_targ.eq.1) then 
                  KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                  keds(idds(ir,2)) =  0.0
               elseif (ofield_targ.eq.2) then 
                  KES(1,IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                  keds(idds(ir,2)) =  kes(1,ir)
               elseif (ofield_targ.eq.3) then 
                  KES(1,IR) = 0.5*((-(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1)
     >                        + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
                  keds(idds(ir,2)) = (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2)
               endif


            else
C     
               kes(1,ir) = 0.0
               keds(idds(ir,2)) = 0.0
               write(6,'(a,2i8,4(1x,g12.5))') 
     >              'KES CALCULATION: IK,IR,DS1,NB1:',
     >              1,ir,ds1,nb1
            endif


         endif
C     
         if (distnks.ge.0.0) then
            if ((smax-kss(nks(ir),ir)).le.distnks) then
               kes(nks(ir),ir) = efnks
            else
               kes(nks(ir),ir) = 0.0
            endif
         else
            DS1 = KSS(NKS(IR),IR) - KSS(NKS(IR)-1,IR)
            DP1 = (KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR)
     >           -KNBS(NKS(IR)-1,IR)*KTEBS(NKS(IR)-1,IR))
            DT1 = (KTEBS(NKS(IR),IR)-KTEBS(NKS(IR)-1,IR))
            NB1 = 0.5*(KNBS(NKS(IR),IR)+KNBS(NKS(IR)-1,IR))

c     
c     jdemod - modified to include and set target data
c     
c     KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1


            DS2 = KSb(NKS(IR),IR) - KSS(NKS(IR),IR)
            DP2 = (KNdS(idds(ir,1))*KTEdS(idds(ir,1))
     >           -KNBS(NKS(IR),IR)*KTEBS(NKS(IR),IR))
            DT2 = (KTEdS(idds(ir,1))-KTEBS(NKS(IR),IR))
            NB2 = 0.5*(KNdS(idds(ir,1))+KNBS(NKS(IR),IR))

            if ((ds1 .ne. 0.0) .and. (nb1.ne.0.0).and.
     >           (ds2 .ne. 0.0) .and. (nb2.ne.0.0)) then 

               if (ofield_targ.eq.1) then 
                  KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                  keds(idds(ir,1))= 0.0
               elseif (ofield_targ.eq.1) then 
                  KES(NKS(IR),IR) = -(1/NB1)*DP1/DS1 - 0.71 * DT1/DS1
                  keds(idds(ir,1)) =  kes(nks(ir),ir)
               elseif (ofield_targ.eq.1) then 
                  KES(nks(ir),IR)=0.5*((-(1/NB1)*DP1/DS1-0.71*DT1/DS1)
     >                    + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
                  keds(idds(ir,1)) = (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2)
               endif

c               KES(nks(ir),IR) =0.5*((-(1/NB1)*DP1/DS1 -0.71 * DT1/DS1)
c     >                             + (-(1/NB2)*DP2/DS2 -0.71 * DT2/DS2))
c               keds(idds(ir,1)) =  (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2)

            else
C     
               kes(nks(ir),ir) = 0.0
               keds(idds(ir,1)) = 0.0
               write(6,'(a,2i8,4(1x,g12.5))') 
     >              'KES CALCULATION: IK,IR,DS1,NB1:',
     >              nks(ir),ir,ds1,nb1
            endif



         endif
c     
c     Loop through all the other knots
c     
         DO 500 IK = 2,NKS(IR)-1
c     
            if (dist1.ge.0.0.and.
     >           kss(ik,ir).le.smax/2.0) then
c     
               if (kss(ik,ir).le.dist1) then
                  kes(ik,ir) = ef1
               else
                  kes(ik,ir) = 0.0
               endif
c     
            elseif (distnks.ge.0.0.and.
     >              kss(ik,ir).gt.smax/2.0) then
c     
               if ((smax-kss(ik,ir)).le.distnks) then
                  kes(ik,ir) = efnks
               else
                  kes(ik,ir) = 0.0
               endif
c     
            else

               DS1 = KSS(IK,IR) - KSS(IK-1,IR)
               DP1 = KNBS(IK,IR)*KTEBS(IK,IR)
     >              -KNBS(IK-1,IR)*KTEBS(IK-1,IR)
               DT1 = (KTEBS(IK,IR)-KTEBS(IK-1,IR))
               NB1 = 0.5*(KNBS(IK,IR)+KNBS(IK-1,IR))
               DS2 = KSS(IK+1,IR) - KSS(IK,IR)
               DP2 = KNBS(IK+1,IR)*KTEBS(IK+1,IR)
     >              -KNBS(IK,IR)*KTEBS(IK,IR)
               DT2 = (KTEBS(IK+1,IR)-KTEBS(IK,IR))
               NB2 = 0.5*(KNBS(IK+1,IR)+KNBS(IK,IR))
               KES(IK,IR) = 0.5*((-(1/NB1)*DP1/DS1-0.71*DT1/DS1)
     >              + (-(1/NB2)*DP2/DS2 - 0.71 * DT2/DS2))
            endif

 500     CONTINUE
c     
 600  CONTINUE
c     
      return
      end
c
c
c
      subroutine flat_t
      implicit none
c
c     FLAT_T: If the upstream temperature flattening options are
c             on - then this routine will modify the upstream
c             temperatures to make them flat using one of  a
c             number of options.
c
      include 'params'
c
      include 'cgeom'
c
      include 'comtor'
c
c     Local variables
c
      integer ir,ik,ikscut1,ikscut2,sfind
      real scut1,scut2
      external sfind
c
      if ((cflatopt.eq.1.or.cflatopt.eq.2.or.cflatopt.eq.3)
     >     .and.ctegcut.gt.0.0) then
c
         do ir = irsep,irwall
c
            if (cflatopt.eq.1.or.cflatopt.eq.2) then 
               scut1 = ctegcut * ksmaxs(ir)
               scut2 = (1.0-ctegcut)*ksmaxs(ir)
            elseif (cflatopt.eq.3) then 
               scut1 = ksb(ikmids(ir),ir) - ctegcut * ksmaxs(ir)
               scut2 = ksb(ikmids(ir),ir) + ctegcut * ksmaxs(ir)
            endif
c
            ikscut1 = sfind(scut1,ir)
            ikscut2 = sfind(scut2,ir)
c
            do ik = 1,nks(ir)
c
               if (cflatopt.eq.1.and.
     >             kss(ik,ir).gt.scut1.and.
     >             kss(ik,ir).lt.scut2) then
c
                  if (kss(ik,ir).gt.(ksmaxs(ir)/2.0)) then
                     ktebs(ik,ir) = ktebs(ikscut2,ir)
                  else
                     ktebs(ik,ir) = ktebs(ikscut1,ir)
                  endif
c
               elseif (cflatopt.eq.2.and.
     >                 kss(ik,ir).gt.scut1.and.
     >                 ktebs(ik,ir).gt.ktebs(ikscut1,ir)) then
c
                   ktebs(ik,ir) = ktebs(ikscut1,ir)
c
               elseif (cflatopt.eq.3.and.
     >                 kss(ik,ir).gt.scut1.and.
     >                 kss(ik,ir).lt.scut2) then
                  if (kss(ik,ir).gt.((scut2+scut1)/2.0)) then
                     ktebs(ik,ir) = ktebs(ikscut2,ir)
                  else
                     ktebs(ik,ir) = ktebs(ikscut1,ir)
                  endif

                  

               endif
c
            end do
c
         end do
c
      endif
c
c     Ti is flat for S > Fact * SMAX for each half ring
c
      if ((cflatopt.eq.1.or.cflatopt.eq.2.or.cflatopt.eq.3)
     >          .and.ctigcut.gt.0.0) then
c
         do ir = irsep,irwall
c
            if (cflatopt.eq.1.or.cflatopt.eq.2) then 
               scut1 = ctigcut * ksmaxs(ir)
               scut2 = (1.0-ctigcut)*ksmaxs(ir)
            elseif (cflatopt.eq.3) then 
               scut1 = ksb(ikmids(ir),ir) - ctigcut * ksmaxs(ir)
               scut2 = ksb(ikmids(ir),ir) + ctigcut * ksmaxs(ir)
            endif
c
            ikscut1 = sfind(scut1,ir)
            ikscut2 = sfind(scut2,ir)
c
            do ik = 1,nks(ir)
c
               if (cflatopt.eq.1.and.
     >             kss(ik,ir).gt.scut1.and.
     >             kss(ik,ir).lt.scut2) then
c
                  if (kss(ik,ir).gt.(ksmaxs(ir)/2.0)) then
                     ktibs(ik,ir) = ktibs(ikscut2,ir)
                  else
                     ktibs(ik,ir) = ktibs(ikscut1,ir)
                  endif
c
               elseif (cflatopt.eq.2.and.
     >                 kss(ik,ir).gt.scut1.and.
     >                 ktibs(ik,ir).gt.ktibs(ikscut1,ir)) then
c
                   ktibs(ik,ir) = ktibs(ikscut1,ir)
c
               elseif (cflatopt.eq.3.and.
     >                 kss(ik,ir).gt.scut1.and.
     >                 kss(ik,ir).lt.scut2) then
                  if (kss(ik,ir).gt.((scut2+scut1)/2.0)) then
                     ktibs(ik,ir) = ktibs(ikscut2,ir)
                  else
                     ktibs(ik,ir) = ktibs(ikscut1,ir)
                  endif
c
               endif
c
            end do
c
         end do
c
      endif
c
c     RETURN
c
      return
      end
c
c
c

      subroutine load_bgopts(in)
      implicit none
      integer in
c
      include 'params'
      include 'comtor'
c
c     LOAD_BGOPTS: This routine handles all the manipulation of the
c                  global variables that must be correctly set for
c                  calling all of the various BG Plasma routines.
c
c     The argument IN is an index into the BGPLASOPT array AND an
c     indicator of how to load the tmp values.
c
c     It is also necessary to save the value of ciopto - the private
c     plasma option switch because for all over-written rings it will
c     need to be set to 1.
c
c     IN = -2 load the current values of cioptg ... into the tmp copies
c     IN = -1 load the values in the tmp copies back into cioptg ...
c     IN =  1,.. load the values from the index given into the cioptg ...
c     IN =  0 Scan the BGPLASOPT array for a 0 ring listing - store the
c             values of cioptg into the tmp variables and load the 0
c             ring set of values.
c
c
c     The BGPLASOPT array contains the data in the following order
c
c     BGPLASOPT(N,1) = Ring number Start
c                ,2  = Ring number End     2 >= 1
c                ,3  = Section specifer 1=first 2=second 3=all
c                ,4  = Plasma Decay option to use
c                ,5  = Sol option for section
c                ,6  = Te gradient option for section
c                ,7  = Ti gradient option for section
c                ,8  = Core option for ring
c                ,9  = E-field option for ring
c
c     CIOPTG   = plasma decay option
c     CIOPTF   = SOL option
c     CIOPTK   = Te gradient option
c     CIOPTL   = Ti Gradient option
c     CCOREOPT = Core option
c     OFIELD   = Efield over-ride option
c
      integer first
      data first /0/
c
      integer tmpcioptg,tmpcioptf,tmpcioptk,tmpcioptl,
     >        tmpccoreopt,tmpofield,tmpciopto
      integer ic,iopt
c slmod begin
      save
c slmod end
c
c     Save original values on first invocation of this routine
c
      if (first.eq.0) then
         tmpcioptg = cioptg
         tmpcioptf = cioptf
         tmpcioptk = cioptk
         tmpcioptl = cioptl
         tmpccoreopt = ccoreopt
         tmpofield = ofield
         tmpciopto = ciopto
c
         first = first +1
c
      endif
c
      if (in.eq.-2) then
         tmpcioptg = cioptg
         tmpcioptf = cioptf
         tmpcioptk = cioptk
         tmpcioptl = cioptl
         tmpccoreopt = ccoreopt
         tmpofield = ofield
         tmpciopto = ciopto
      elseif (in.eq.-1) then
         cioptg  =   tmpcioptg
         cioptf  =   tmpcioptf
         cioptk  =   tmpcioptk
         cioptl  =   tmpcioptl
         ccoreopt=   tmpccoreopt
         ofield  =   tmpofield
         ciopto  =   tmpciopto
      elseif (in.eq.0) then
c
         iopt = 0
c
         do ic = 1,nbgplas
             if (bgplasopt(ic,1).eq.0) then
               iopt = ic
               goto 100
            endif
         end do
 100     continue
c
c        Copy defaults into tmp values
c
         tmpcioptg = cioptg
         tmpcioptf = cioptf
         tmpcioptk = cioptk
         tmpcioptl = cioptl
         tmpccoreopt = ccoreopt
         tmpofield = ofield
         tmpciopto = ciopto
c
c        If a set of zero ring values was found then
c        copy these into the in use values.
c
         if (iopt.gt.0) then
            cioptg  =  bgplasopt(iopt,4)
            cioptf  =  bgplasopt(iopt,5)
            cioptk  =  bgplasopt(iopt,6)
            cioptl  =  bgplasopt(iopt,7)
            ccoreopt=  bgplasopt(iopt,8)
            ofield  =  bgplasopt(iopt,9)
         endif
c
c     For any other LEGAL value of IN
c
      elseif (in.le.nbgplas) then
c
         cioptg  =  bgplasopt(in,4)
         cioptf  =  bgplasopt(in,5)
         cioptk  =  bgplasopt(in,6)
         cioptl  =  bgplasopt(in,7)
         ccoreopt=  bgplasopt(in,8)
         ofield  =  bgplasopt(in,9)
c
      else
c
         write (6,*) 'ERROR: Illegal value of Index: ',in
         write (6,*) '       Found in LOAD_BGOPTS from BGPLASMA'
c
         write (0,*) 'ERROR: Illegal value of Index: ',in
         write (0,*) '       Found in LOAD_BGOPTS from BGPLASMA'
c
      endif
c
      return
      end
c
c
c
      subroutine PINEXE(title,equil,lpinopt,litersol,iitersol,
     >                  liter,tmpcsopt,tmpcioptf,lpinavail,
     >                  iiterpin)
      use error_handling
c slmod begin
      USE mod_sol28_global
c slmod end
      implicit none
      logical lpinopt,litersol,liter,lpinavail
      integer iitersol,tmpcsopt,tmpcioptf,iiterpin
      character*(*) title,equil
c
      include 'params'
      include 'dynam1'
      include 'cgeom'
      include 'comtor'
      INCLUDE 'slcom'
c
c     PINEXE: This routine contains the code that will set up and
c             call the appropriate hydrogenic neutral code.
c
c
c     Local variables
c
      integer retcode,ik,ir
      real pintim
c slmod begin
      integer iparam
      character pin_command*1024
c slmod end
c
c     Initialization
c
      pintim = 0.0
      liter = .false.
c

      IF (LPINOPT) THEN
c
c       Set up for PIN call - only for JET grids
c
        retcode = 0
c
c slmod begin - new
        IF (pincode.NE.4.AND.pincode.NE.5.AND.
     .      (rel_opt.EQ.1.OR.rel_opt.EQ.3)) CALL StoreSources(-1)

        IF (rel_opt.EQ.2.AND.rel_count.EQ.0) CALL SaveSolution
c *TEMP*
c         CALL SaveSolution
        IF (rel_opt.EQ.1.AND.
     .      (rel_iter.EQ.rel_niter.OR.rel_step.EQ.0)) CALL SaveSolution

        IF (s28ion.GT.0) s28ionset = .TRUE.
        IF (s28rec.GT.0) s28recset = .TRUE.
        IF (s28mom.GT.0) s28momset = .TRUE.
        IF (s28ionpfz.GT.0) s28ionsetpfz = .TRUE.
        IF (s28recpfz.GT.0) s28recsetpfz = .TRUE.
        IF (s28mompfz.GT.0) s28momsetpfz = .TRUE.

        if (pincode.EQ.0) then
c
c        if (pin_code.eq.0) then
c slmod end
c
           CALL WRTPIN(title,equil,17)
c
           write(0,*) 'Calling PIN for iteration ',iitersol
           CALL INVOKEPIN(ACTPIN,pintim,retcode)
           write(0,*) 'Return from PIN after ',pintim,' (seconds)'
           write(6,*) 'AFTER PIN:'
c
c          Load PIN results
c
           CALL READPIN
c
c       Set up for Eirene call - only for SONNET grids
c
c slmod begin - new
        elseif (pincode.EQ.1.OR.pincode.EQ.2.OR.pincode.EQ.3) then
c
c        elseif (pin_code.EQ.2) then
c slmod end
c
c          changed by Krieger IPP 12/94
c          file descriptor 17 hidden in b2wrpl
c
c slmod begin - new
           IF     (pincode.EQ.1) THEN
             call wrteirene_97
             CALL WriteGeometryFile_97
             CALL WriteInputFile_97
             write(0     ,*) 'Calling EIRENE97 for iteration ',iitersol
             write(PINOUT,*) 'Calling EIRENE97 for iteration ',iitersol
           ELSEIF (pincode.EQ.2.OR.pincode.EQ.3) THEN
             CALL WrtEIRENE
             IF (rel_opt.EQ.1.OR.rel_opt.EQ.2.OR.rel_opt.EQ.3) THEN
               WRITE(0     ,'(A,I3,A,I2,A,I2,A)')
     .           ' Calling EIRENE for iteration ',iitersol - 1,
     .           ' (step ',rel_step,'  subiteration ',rel_iter,
     .           ')'
               WRITE(PINOUT,'(A,I3,A,I2,A,I2,A)')
     .           ' Calling EIRENE for iteration ',iitersol - 1,
     .           ' (step ',rel_step,'  subiteration ',rel_iter,
     .           ')'
             ELSE
               WRITE(0,*) 'Calling EIRENE99 for iteration ',iitersol
               WRITE(PINOUT,*) 'Calling EIRENE99 for iter ',iitersol
             ENDIF
           ELSE
             CALL ER('PINEXE','Invalid PINCODE value',*99)
99           STOP
           ENDIF

           CALL InvokePIN(actpin,pintim,retcode)
           WRITE(0     ,'(A,I6,A)') ' Return from EIRENE after ',
     .                              NINT(pintim),' s'
           WRITE(PINOUT,'(A,I6,A)') ' Return from EIRENE after ',
     .                              NINT(pintim),' s'
           WRITE(6,*) 'AFTER PIN:'
c
c          Load EIRENE results
c
           IF     (pincode.EQ.1) THEN
             call readeire_97
           ELSEIF (pincode.EQ.2.OR.pincode.EQ.3) THEN
             call readeire
           ENDIF

        ELSEIF (pincode.EQ.4.OR.pincode.EQ.5) THEN
c...       Calling EIRENE04/06/07:
           SELECTCASE (pincode)
             CASE(4)
               CALL WriteEireneFiles_04
             CASE(5)
               CALL WriteEireneFiles_06(iitersol)
           ENDSELECT

           IF (rel_opt.NE.0) THEN
             WRITE(0     ,'(A,I3,A,I2,A,I2,A)')
     .         ' Calling EIRENE0x, iteration ',iitersol - 1,
     .         ' (step ',rel_step,'  subiteration ',rel_iter,
     .         ')'
             WRITE(PINOUT,'(A,I3,A,I2,A,I2,A)')
     .         ' Calling EIRENE0x, iteration ',iitersol - 1,
     .         ' (step ',rel_step,'  subiteration ',rel_iter,
     .         ')'
           ELSE
             WRITE(0     ,*) 'Calling EIRENE0x, iteration',iitersol,
     .                                                     nitersol
             WRITE(PINOUT,*) 'Calling EIRENE0x, iteration',iitersol
           ENDIF
c...       Pass iteration number to the EIRENE script if tetrahedrons
c          are being called:  *** HACK *** special for filaments at the moment...
           IF (citersol.GT.0.AND.
     .         (eirgeom.EQ.3.OR.opt_eir%ntime.NE.0)) THEN
             iparam = iitersol    
           ELSE
             iparam = -1
           ENDIF
           WRITE(pin_command,'(A,I5.3)') TRIM(actpin),iparam
           WRITE(0,*) 'PIN_COMMAND:',TRIM(pin_command)
           CALL InvokePIN(pin_command,pintim,retcode)
c           CALL InvokePIN(actpin,pintim,retcode)

           WRITE(0     ,'(A,I6,A)') ' Return from EIRENE after ',
     .                              NINT(pintim),' s'
           WRITE(PINOUT,'(A,I6,A)') ' Return from EIRENE after ',
     .                              NINT(pintim),' s'
           WRITE(6,*) 'AFTER PIN:'

c...       Read PIN results:
           SELECTCASE (pincode)
             CASE(4)
               CALL ReadEireneResults_04
             CASE(5)
               CALL ReadEireneResults_06(iitersol)
           ENDSELECT

c
c             call wrteirene
c             write(0,*) 'Calling EIRENE for iteration ',iitersol
c
c             CALL WriteGeometryFile
c             CALL WriteInputFile
c
c           STOP
c
c
c           call invokepin(actpin,pintim,retcode)
c           write(0,*) 'Return from EIRENE after ',pintim,' (seconds)'
c           write(6,*) 'AFTER PIN:'
c
c          Load EIRENE results
c
c           call readeire
c slmod end
c
c          end of changes by Krieger IPP 12/94
c
        endif
c
c       Check for error in PIN - use time for run as first indicator
c       until return codes are available
c
c       Choose 30 seconds as time limit for now.
c
c       jdemod - reduce the time limit to 2 seconds for exit and 10 seconds for warning
c
c slmod begin
        if ((retcode.ne.0.or.pintim.lt.10.0).and.
     .      pincode.ne.4.and.pincode.ne.5) then
c
           if (pintim.lt.10.0) then 
c
c             jdemod - change to warning for 10s
c
              write(error_message_data,'(a,f8.3,a)')
     >         'PIN WARNING: LESS THAN 10 SECONDS USED IN PIN'//
     >         ': TIME USED = ',pintim,' (S) : CHECK FOR'//
     >         ' INPUT ERROR'
              call errmsg(error_message_data)


           elseif (pintim.lt.2.0) then
c
c             jdemod - exit if less than 2 seconds
c
              write(error_message_data,'(a,f8.3,a)')
     >         'PIN ERROR: LESS THAN 2 SECONDS USED IN PIN'//
     >         ': TIME USED = ',pintim,' (S) : CHECK FOR'//
     >         ' INPUT ERROR: PROGRAM HALT'
              call errmsg(error_message_data)
              stop 'PIN < 2 SECONDS'

           endif
c
c          check PIN return code
c
           if (retcode.ne.0) then
c
              retcode = retcode/256
c
c             Issue error message for PIN return code 
c             jdemod - not sure about the divide by 256 unless the 
c                      retcode starts off large
c
              call errmsg('ERROR: PIN DID NOT EXECUTE CORRECTLY:'//
     >                    ' PROGRAM HALT : PIN RETURN CODE =',retcode)
c
              stop 'PIN EXECUTION ERROR'
c
           endif
c
        endif
c
        iiterpin =  iiterpin + 1
c
c       Calculate equivalent lengths
c
        call calcleq
c
        totpintim = totpintim + pintim
c
        lpinavail = .true.
c
        IF (litersol) THEN
c
C
C          SET IONIZATION OPTION TO BE PIN DATA
C
           piniter = .true.
c
           tmpcioptf = cioptf
           tmpcsopt   = csopt
           CIOPTF = CSECSOL
           CSOPT = CSOPT2
C
C          ALLOW FLUX RECIRCULATION FOR IONIZATION OPTION 2
C
           IF (CSOPT.EQ.2) FLUXROPT = 1
C
C          AFTER PRESCRIBED NUMBER OF ITERATIONS, SWITCH OFF CALL TO PIN
C
           IF (IITERSOL.GE.NITERSOL) litersol = .false.
c slmod begin - tmp
           WRITE(PINOUT,*) 'IITER,NITER,LITER=',
     .       iitersol,nitersol,litersol,lpinopt
c slmod end
C
C          STORE VALUES OF PLASMA PARAMETERS FOR CONVERGENCE STUDY
C          AND CALCULATE CHI-SQUARED VARIATION FROM LAST ITERATION
C
           DO 362 IR=IRSEP,IRWALL-1,1
C
             DO 363 IK=1,NKS(IR),1
C
c
c            jdemod - chisq doesn't work since it isn't initialized from 
c                     what I can see 
c
             IF (IITERSOL.GT.1) THEN
C
             CHISQ1(IITERSOL-1)=CHISQ1(IITERSOL-1)+(KTEBS(IK,IR)-
     >                          OKTEBS(IK,IR))**2
             CHISQ2(IITERSOL-1)=CHISQ2(IITERSOL-1)+(KTIBS(IK,IR)-
     >                          OKTIBS(IK,IR))**2
c
             if (oknbs(ik,ir).gt.0.0) then
                CHISQ3(IITERSOL-1)=CHISQ3(IITERSOL-1)+
     >                  DBLE((KNBS(IK,IR)- OKNBS(IK,IR)))**2
             endif
c
             CHISQ4(IITERSOL-1)=CHISQ4(IITERSOL-1)+(KVHS(IK,IR)-
     >                          OKVHS(IK,IR))**2
             CHISQ5(IITERSOL-1)=CHISQ5(IITERSOL-1)+(KES(IK,IR)-
     >                          OKES(IK,IR))**2
C
             ENDIF
C
             OKTEBS(IK,IR) = KTEBS(IK,IR)
             OKTIBS(IK,IR) = KTIBS(IK,IR)
             OKNBS(IK,IR) = KNBS(IK,IR)
             OKVHS(IK,IR) = KVHS(IK,IR)
             OKES(IK,IR) = KES(IK,IR)
C
 363         CONTINUE
C
 362       CONTINUE
C
C          ITERATE!!!
C
           IITERSOL = IITERSOL + 1
           liter = .true.
c slmod begin
           IF (citersol.EQ.2.AND..NOT.litersol) liter = .FALSE.
c slmod end
        ENDIF
c
c slmod begin - new
c...    Relax PIN volume sources (if EIRENE04 not being used, relaxation done as sources
c       are loaded): 
        IF (pincode.NE.4.AND.pincode.NE.5.AND.
     .      (rel_opt.EQ.1.OR.rel_opt.EQ.3)) CALL UpdateSources(-1)

c...    I do not particularly like this here since the plasma solution
c       and the sources are not quaranteed to be consistent:
        IF (rel_opt.EQ.2.AND.rel_count.GE.1.AND.rel_iter.EQ.rel_niter)
     .    CALL SaveSolution

c...    Change output option:
        IF (outmode.GE.2) THEN
          WRITE(67,'(A,I4,A)') 'ITERATION ',rel_count,' (after update)'
          WRITE(67,*)
          CALL OutputEIRENE(67,'DONE EIRENE ITERATION')
          WRITE(67,*)
          CALL OutputAdditionalCellData
        ENDIF

        IF (asc_ncell.GT.0) CALL OutputAdditionalCellData

      ELSE
c...    Save solution in supplimental .RAW file(s):
c        WRITE(0,*) 'NEW SAVE SOLUITON'
c        IF (rel_opt.GT.0) CALL SaveSolution
      ENDIF
c
c      ENDIF
c slmod end
c
c
      return
      end
c
c
c
c slmod begin
c
c ======================================================================
c
c subroutine: MapParameters
c
c
      SUBROUTINE MapParameters(ir,type,s,v,n,MAXVAL)
      IMPLICIT none
       
      INTEGER ir,type,n,MAXVAL
      REAL    s(MAXVAL),v(MAXVAL)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      INTEGER i0,i1,ik,id
      REAL    frac,v0,v1
      REAL*8  a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd

      n = 0
      s = 0.0
      v = 0.0

c...  Step-based interpolation fraction:
      frac = 0.0
c      frac = MAX(0.0,REAL(rel_step - 1) / REAL(rel_nstep - 1))
c      frac = rel_bound1 + frac * (rel_bound2 - rel_bound1) 
      
      DO i1 = 2, osmns28
        i0 = i1 - 1
        IF (osms28(i0,1).NE.osms28(i1,1)) CYCLE

        IF (s28mode.GE.2.0.AND.
     .      ((osms28(i1,9) .NE.0.0.AND.REAL(ir).LT.osms28(i1,9) ).OR.
     .       (osms28(i1,10).NE.0.0.AND.REAL(ir).GT.osms28(i1,10))))
     .    CYCLE

        a1 = DBLE(osms28(i1-1,3))
        a2 = DBLE(osms28(i1-1,4))
        b1 = DBLE(osms28(i1  ,3))
        b2 = DBLE(osms28(i1  ,4))
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          c1 = 0.5 * DBLE(rvertp(1,id) + rvertp(2,id))
          c2 = 0.5 * DBLE(zvertp(1,id) + zvertp(2,id))
          d1 = 0.5 * DBLE(rvertp(3,id) + rvertp(4,id))
          d2 = 0.5 * DBLE(zvertp(3,id) + zvertp(4,id))
      
          CALL CalcInter(a1,a2,b1,b2,c1,c2,d1,d2,tab,tcd)

          IF (tab.GE.0.0.AND.tab.LT.1.0.AND.
     .        tcd.GE.0.0.AND.tcd.LT.1.0) THEN

            IF     (osms28(i1,1).EQ.REAL(type)) THEN

              n = n + 1

              IF (n.GT.MAXVAL) CALL ER('MapParameters','Array bounds '//
     .                                 'exceeded',*99)

c...          Assign the field-line coordinate:
              s(n) = ksb(ik-1,ir)+SNGL(tcd)*(ksb(ik,ir)-ksb(ik-1,ir))

c...          Assign the quantity of interest, with 2D interpolation and interpolation
c             with solution interation:
              v0 = osms28(i0,5) + frac*(osms28(i0,7)-osms28(i0,5))
              v1 = osms28(i1,5) + frac*(osms28(i1,7)-osms28(i1,5))
              v(n) = v0 + SNGL(tab)**osms28(i1,2) * (v1 - v0)

            ELSEIF (osms28(i1,1).GT.12.0) THEN
              CALL ER('MapParameters','Invalid parameter '//
     .                'type',*99)

            ENDIF
          ENDIF
        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: PrescribeFlow
c
c
      SUBROUTINE PrescribeFlow
      IMPLICIT none
       
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

      REAL GetCs

      INTEGER MAXVAL
      PARAMETER (MAXVAL=100)

      INTEGER ik,ir,n,i1,type,ik1,ik2,oldte(MAXNKS)
      LOGICAL status
      REAL    s(MAXVAL),v(MAXVAL),shold,vhold,machno


c...  Limit calculation to SOL for now:
      DO ir = irsep, nrs
c      DO ir = irsep, irwall-1

c...    Skip virtual rings:
        IF (idring(ir).EQ.-1) CYCLE

        IF (.TRUE.) THEN
c...      Parameter identifier in the data list in the input file:
          type = 10

c...      Check for 2D parameter specification:
          CALL MapParameters(ir,type,s,v,n,MAXVAL)

c...      Data found:
          IF (n.GT.0) THEN
c...        Add end points, assuming sonic flow to both targets (fix later -- base 
c           on plasma solution):            
            n = n + 2
            s(n-1) =  0.0
            v(n-1) = -1.0
            s(n)   =  ksmaxs(ir)
            v(n)   =  1.0

c...        Sort data points in ascending order:
            status = .TRUE.
            DO WHILE (status)
              status = .FALSE.
              DO i1 = 1, n-1
                IF (s(i1).GT.s(i1+1)) THEN
                  shold = s(i1)
                  vhold = v(i1)
                  s(i1) = s(i1+1)
                  v(i1) = v(i1+1)
                  s(i1+1) = shold
                  v(i1+1) = vhold
                  status = .TRUE. 
                ENDIF
              ENDDO
            ENDDO

c...        Quick check on if all entries were valid:
            DO i1 = 1, n
              IF (s(i1).LT. 0.0.OR.s(i1).GT.ksmaxs(ir)) 
     .          CALL ER('PrescribeFlow','s-value out of range',*99)
              IF (v(i1).LT.-1.0.OR.v(i1).GT.1.0) 
     .          CALL ER('PrescribeFliow','v-value out of range',*99)
c              WRITE(0,*) '--->',ir,i1,s(i1),v(i1)
            ENDDO

c...        Interpolate for each cell on the ring and calculate
c           flow from specified Mach number and the Te,Ti values 
c           (already calculated elsewhere):
            DO ik = 1, nks(ir)
              CALL Fitter(n,s,v,1,kss(ik,ir),machno,'LINEAR')

              kvhs(ik,ir) = machno * GetCs(ktebs(ik,ir),ktibs(ik,ir))
            ENDDO

          ENDIF

         ELSE
           CALL ER('PrescribeFlow','Invalid option',*99)
         ENDIF

      ENDDO




c...  Also process the temperature profile to get rid of
c     discontinuity:
      DO ir = irsep, irwall-1

        IF (.FALSE.) THEN
c        IF (.TRUE.) THEN

c...      Find midplane bounds:
          ik1 = 0
          ik2 = 0
          DO ik = 1, nks(ir)-1
            IF (zs(ik,ir).LT.0.0.AND.zs(ik+1,ir).GT.0.0.OR.
     .          zs(ik,ir).GT.0.0.AND.zs(ik+1,ir).LT.0.0) THEN
              IF (ik1.EQ.0) ik1 = ik 
              IF (ik1.NE.0) ik2 = ik
            ENDIF
          ENDDO
          IF (ik1.EQ.0.OR.ik2.EQ.0) 
     .      CALL ER('...','Midplane cell(s) not found',*99)

c...      Copy Te temperature profile:
          DO ik = 1, nks(ir)
            oldte(ik) = ktebs(ik,ir)
          ENDDO
          DO ik = ik1+2, ik2-2
            ktebs(ik,ir) = (0.15 * oldte(ik-2) + 0.20 * oldte(ik-1) +
     .                      0.30 * oldte(ik  ) + 0.20 * oldte(ik+1) + 
     .                      0.15 * oldte(ik+2)) * 1.0
            ktibs(ik,ir) = ktebs(ik,ir)
          ENDDO 

        ENDIF

      ENDDO



      RETURN
99    STOP
      END
c
c ======================================================================
c
c slmod end
c
c
      subroutine recalculate_bg_velocity(halfopt)
      implicit none
      integer halfopt
c
      include 'params'
      include 'comtor'
      include 'cgeom'
      include 'pindata'
c
c     RECALCULATE_BG_VELOCITY:
c
c     This routine recalculates the background velocity - overwriting 
c     whatever was calculated by the specified SOL option. 
c
c     It does this by gathering the particle source data for the
c     ring - ionization, recombination, target out fluxes - and 
c     assumes an evenly distributed cross-field outflux to make up
c     any differences.
c
c     Using the input target flux and the density from the SOL option
c     solution it then calculates the parallel flux in every cell along
c     the ring and uses G=nv to find the velocity. 
c
c     Some data used in this routine:
c     KNBS    = cell ion density
c     KNDS    = target density
c     KVDS    = target velocity 
c     PINION  = ionization data
c     PINREC  = recombination data
c      
c     This routine will calculate the velocity field on either a full
c     ring or a half-ring basis. In the half-ring method the velocity is 
c     forced to zero at the mid-point of the ring and the plasma solution
c     for the inner target will not affect the outer target flow
c     field. For the full ring basis, the out flux at both targets as 
c     well as the sources is taken into account in the calculations. 
c
c     The code is written so as to calculate the data in half-ring segments
c     - the difference in the treatment is in the precalculation of the 
c     magnitude of the cross-field source/sink term.
c
c     Local variables:
c
      integer ik,ir
      real gamma(maxnrs,2)
      real cfterm(maxnrs,2)
      real src1,src2,src,g,sl,su
c
      real srci1, srcr1,srci2,srcr2
c
c     Calculate the target out fluxes
c     Target 1 = INNER X-pt up / OUTER Xpt down = target at IK=NKS(IR) index
c     Target 2 = OUTER X-pt up / INNER Xpt down = target at IK=1 index
c
      do ir = irsep,nrs
c
         gamma(ir,1) = knds(idds(ir,1)) * abs(kvds(idds(ir,1)))
         gamma(ir,2) = knds(idds(ir,2)) * abs(kvds(idds(ir,2)))
c
         src1 = 0.0
         src2 = 0.0
         srci1 = 0.0
         srci2 = 0.0
         srcr1 = 0.0
         srcr2 = 0.0         
c
         write(6,'(a,1x,i6,1p,2(1x,g12.5))') 'GAMMA:',ir,
     >                 gamma(ir,1),gamma(ir,2)
c
         do ik = 1,nks(ir)
c
            if (ksb(ik-1,ir).lt.ksmaxs(ir)/2.0.and.
     >          ksb(ik,ir).lt.ksmaxs(ir)/2.0) then 
c
                srci2 = srci2 + pinion(ik,ir)*(ksb(ik,ir)-ksb(ik-1,ir))
                srcr2 = srcr2 + pinrec(ik,ir)*(ksb(ik,ir)-ksb(ik-1,ir))
c
                src2 = src2 + pinion(ik,ir) * (ksb(ik,ir)-ksb(ik-1,ir))
     >                      - pinrec(ik,ir) * (ksb(ik,ir)-ksb(ik-1,ir))
                

            elseif (ksb(ik-1,ir).gt.ksmaxs(ir)/2.0.and.
     >              ksb(ik,ir).gt.ksmaxs(ir)/2.0) then 

                srci1 = srci1 + pinion(ik,ir)*(ksb(ik,ir)-ksb(ik-1,ir))
                srcr1 = srcr1 + pinrec(ik,ir)*(ksb(ik,ir)-ksb(ik-1,ir))
c
                src1 = src1 + pinion(ik,ir) * (ksb(ik,ir)-ksb(ik-1,ir))
     >                      - pinrec(ik,ir) * (ksb(ik,ir)-ksb(ik-1,ir))
            else
c
                srci2 = srci2 + 
     >                  pinion(ik,ir)*(ksmaxs(ir)/2.0-ksb(ik-1,ir))
c
                srcr2 = srcr2 + 
     >                  pinrec(ik,ir)*(ksmaxs(ir)/2.0-ksb(ik-1,ir))
c
                src2 = src2+pinion(ik,ir)*(ksmaxs(ir)/2.0-ksb(ik-1,ir))
     >                     -pinrec(ik,ir)*(ksmaxs(ir)/2.0-ksb(ik-1,ir))

c
                srci1 = srci1 + 
     >                  pinion(ik,ir)*(ksb(ik,ir)-ksmaxs(ir)/2.0)
                srcr1 = srcr1 + 
     >                  pinrec(ik,ir)*(ksb(ik,ir)-ksmaxs(ir)/2.0)
c
                src1= src1 + pinion(ik,ir)*(ksb(ik,ir)-ksmaxs(ir)/2.0)
     >                      - pinrec(ik,ir)*(ksb(ik,ir)-ksmaxs(ir)/2.0)
c
            endif
         end do 
c
         write(6,'(a,1x,i6,1p,3(1x,g12.5))') 'SRC 1:',ir,
     >                src1,srci1,srcr1
         write(6,'(a,1x,i6,1p,3(1x,g12.5))') 'SRC 2:',ir,
     >                src2,srci2,srcr2
c
c        Calculate cross-field term
c
c        Half-ring option
c
         if (halfopt.eq.1) then

            cfterm(ir,1) = (gamma(ir,1) - src1) / (ksmaxs(ir)/2.0)  
            cfterm(ir,2) = (gamma(ir,2) - src2) / (ksmaxs(ir)/2.0)  

c
c        Full-ring option
c
         elseif (halfopt.eq.2) then  

           cfterm(ir,1) = (gamma(ir,1)+gamma(ir,2)
     >                     -(src1+src2))/ksmaxs(ir)
           cfterm(ir,2) = cfterm(ir,1)

         endif 

         write(6,'(a,i6,1p,2(1x,g12.5))') 'CFTERM:',ir,
     >                cfterm(ir,1)*ksmaxs(ir)/2.0,
     >                cfterm(ir,2)*ksmaxs(ir)/2.0


      end do

c
c     Recaculate the velocity field
c
      do ir = irsep,nrs

         src = 0.0
c
c        target 2 - velocity is negative so -gamma is used to start
c
         do ik = 1,ikmids(ir)
c
            su = ksb(ik,ir) - kss(ik,ir)
            sl = kss(ik,ir) - ksb(ik-1,ir)
c
            src = src + 
     >           (pinion(ik,ir) - pinrec(ik,ir) + cfterm(ir,2)) * sl
c
            g = -gamma(ir,2) + src
c
c           Calculate velocity 
c
            kvhs(ik,ir) = g / knbs(ik,ir)
c
            write(6,'(a,1x,i6,1x,i6,4(1x,g12.5))') 'VBO:',
     >             ik,ir,kvhs(ik,ir),knbs(ik,ir),g,src
c
            src = src + 
     >           (pinion(ik,ir) - pinrec(ik,ir) + cfterm(ir,2)) * su
c
         end do

c
c        Target 1 - velocity is positive so signs change in expressions
c
         src = 0.0
c
         do ik = nks(ir),ikmids(ir)+1,-1
c
c           Cell sizes
c
            sl = ksb(ik,ir) - kss(ik,ir)
            su = kss(ik,ir) - ksb(ik-1,ir)
c
            src = src + 
     >           (pinion(ik,ir) - pinrec(ik,ir) + cfterm(ir,1)) * sl
c
            g = gamma(ir,1) - src
c
c           Calculate velocity 
c
            kvhs(ik,ir) = g / knbs(ik,ir)
c
            write(6,'(a,1x,i6,1x,i6,4(1x,g12.5))') 'VBO:',
     >             ik,ir,kvhs(ik,ir),knbs(ik,ir),g,src
c
            src = src + 
     >           (pinion(ik,ir) - pinrec(ik,ir) + cfterm(ir,1)) * su
c
         end do 

 
      end do
c
c
c
      return
      end
c
c
c
      subroutine set_bg_velocity(velopt)
      implicit none
      integer velopt
      include 'params'
      include 'cgeom'
      include 'comtor'
c
c     SET_BG_VELOCITY:
c
c     This routine overwrites the velocity values stored in the 
c     kvhs(ik,ir) array. It is intended as a method to quickly 
c     apply a simple velocity distribution for testing purposes. 
c
      integer ik,ir
      real ssfact
c
c     VELOPT 0: Overwrite the velocity array with 0.0
c
      if (velopt.eq.0) then 
         write(0,*) 'Velocity override 4 - velopt 0'
         call rzero(kvhs,maxnks*maxnrs)
      elseif (velopt.eq.1) then 
         write(0,*) 'Velocity override 5 - velopt 1'
         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(4.0e3,ssfact)
               else
                  kvhs(ik,ir) = sign(2.0e4,ssfact)
               endif
            end do
         end do
      elseif (velopt.eq.2) then 
         write(0,*) 'Velocity override 6 - velopt 2'
         do ir=irsep,nrs 
            do ik = 1,nks(ir)
               
               if (ik.gt.nks(ir)/2) then 
                  kvhs(ik,ir) =  cvhout
               else 
                  kvhs(ik,ir) = -cvhout
               endif

            end do 
         end do
      elseif (velopt.eq.3) then 
         write(0,*) 'Velocity override 7 - velopt 3'
         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(4.0e3,ssfact)
               elseif (ssfact.gt.0.04.and.ir.le.irsep+2) then 
                  kvhs(ik,ir) = -sign(1.0e4,ssfact)
               else
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               endif
            end do
         end do
      
      elseif (velopt.eq.4) then 
         write(0,*) 'Velocity override 8 - velopt 4'
         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(1.0e3,ssfact)
               elseif (ssfact.gt.0.04.and.ir.le.irsep+2) then 
                  kvhs(ik,ir) = -sign(5.0e3,ssfact)
               else
                  kvhs(ik,ir) = sign(1.0e3,ssfact)
               endif
            end do
         end do
      elseif (velopt.eq.5) then 
         write(0,*) 'Velocity override 9 - velopt 5'
         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(0.0,ssfact)
               elseif (ssfact.gt.0.03.and.ir.le.irsep+2) then 
                  kvhs(ik,ir) = -sign(5.0e3,ssfact)
               else
                  kvhs(ik,ir) = sign(5.0e4,ssfact)
               endif
            end do
         end do
      elseif (velopt.eq.6) then 
         write(0,*) 'Velocity override 10 - velopt 6'

         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(0.0,ssfact)
               elseif (ssfact.gt.0.03.and.ir.le.irsep+2) then 
                  kvhs(ik,ir) = -sign(5.0e3,ssfact)
               else
                  kvhs(ik,ir) = sign(2.0e4,ssfact)
               endif
            end do
         end do
      
      elseif (velopt.eq.7) then 
         write(0,*) 'Velocity override 11 - velopt 7'

         do ir = irsep,nrs
            do ik = 1,nks(ir)

               if (ik.gt.nks(ir)/2) then 
                  ssfact = -(kss(ik,ir)-ksmaxs(ir))/ksmaxs(ir)
               else 
                  ssfact = kss(ik,ir)/ksmaxs(ir)
               endif
c
               if (ssfact.lt.0.005) then 
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               elseif (ssfact.gt.0.005.and.ssfact.lt.0.015) then
                  kvhs(ik,ir) = sign(0.0,ssfact)
               elseif (ssfact.gt.0.03.and.ir.le.irsep+2) then 
                  kvhs(ik,ir) = -sign(5.0e3,ssfact)
               else
                  kvhs(ik,ir) = sign(1.0e4,ssfact)
               endif
            end do
         end do
      
      endif

      return
      end
c
c
c
      subroutine overlay_plasma
      use error_handling
      use plasma_overlay
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      
! read in the overlay plasma file including array sizes and bounds
! Loop through grid and for any cell with NON-ZERO values within
! the region covered - replace the values with the ones interpolated from 
! the file. This includes target data. This is only done for ne,Te at the
! present time but the code could be enhanced to interpolate any number
! of plasma parameters for any kind of plasma solution. 

      integer ierr
      integer ik,ir,id
      real r,z

! exit if option is not set
      if (external_plasma_overlay.ne.1) return
      
      call load_plasma_overlay(external_plasma_file,ierr)

      if (ierr.ne.0) return

      ! overlay grid first
      do ir = 1,nrs
         do ik = 1,nks(ir)
            r = rs(ik,ir)
            z = zs(ik,ir)
            call interpolate_overlay(r,z,ktebs(ik,ir),knbs(ik,ir))
         end do
      end do

      ! Now do the same for the target values

      do id = 1,nds
         r = rp(id)
         z = zp(id)
         call interpolate_overlay(r,z,kteds(id),knds(id))
      end do

      ! That should have replaced ALL grid cells in the region
      ! with data from the DTS fit

      return
      end








