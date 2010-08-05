c     -*-Fortran-*-
c     @PROCESS NOOPT
      SUBROUTINE PRDATA (NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
C     THIS ROUTINE PRINTS ALL THE OPTION FLAGS, TORUS ATTRIBUTES ETC
C     STORED IN COMMONS COMTOR AND COMTAU
C
C     C.M.FARRELL    NOVEMBER 1987
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
      include    'printopt' 
c
      INTEGER  IT,IZ,IN
c
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
c      external lenstr
c
c
c
c     Initialization
c
c      sp = '                        '
c      sa = '  '
c      s1 = '      '
c      s2 = '            '
c
      prsol21 = .false.
      prsol22 = .false.


c
c     Set switches for turning on SOL 21 and SOL 22 print outs
c
      if (cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) then

         do in = 1,nbgplas
c
c           SET Flags for printing SOL 21 and SOL 22 if they are
c           specified as part of the background plasma solution
c
            if (bgplasopt(in,5).eq.21) prsol21 = .true.
            if (bgplasopt(in,5).eq.22) prsol22 = .true.

         end do
c
      endif
c
C-----------------------------------------------------------------------
c
c     Print Table of Content Information
c
      call pr_toc


c
C-----------------------------------------------------------------------
c
c     Print Geometry and grid information 
c
      call pr_geom 
c
C-----------------------------------------------------------------------
c
c     Print Background Plasma Information
c
      call pr_bg(NIZS,NIMPS,NIMPS2,nymfs) 
c
C-----------------------------------------------------------------------
c
c     Print Impurity Simulation Information
c
      if (ctestsol.eq.0.0) then 

         call pr_sim(NIZS,NIMPS,NIMPS2,nymfs)  

      endif   
c
C-----------------------------------------------------------------------
c
c     Print Specific Option Selections
c
      call pr_options(NIZS,NIMPS,NIMPS2,nymfs) 
c
c
C-----------------------------------------------------------------------
c     Print Dwell times as debugging output
c
      IF (IMODE.EQ.1) THEN
        WRITE (6,9020) (DWELFS(IT),IT=1,NTS)
        DO 30 IZ = 0, NIZS
          IF (IZ.EQ.0 .AND. CNEUTA.NE.0) GOTO 30
          WRITE (6,9021) IZ,(DWELTS(IZ)*DWELFS(IT),IT=1,NTS)
   30   CONTINUE
      ENDIF


C-----------------------------------------------------------------------
C
C---- CHECK FOR DUBIOUS / NON-IMPLEMENTED  COMBINATIONS ...
C
      CONTINUE
      IF (CIOPTI.NE.0.AND.CIZB.NE.1) THEN
       CALL PRB
       CALL PRC ('*** WARNING *** C-X RECOM SPECIFIED FOR NON-HYDROGENIC
     > PLASMA...')
      ENDIF
c
C
c

      RETURN
 9010 FORMAT(1X,A,F9.2)
 9020 FORMAT(/1X,'DWELL TIME FACTORS AND OUTPUT TIMES (S)',
     >  /1X,'IZ ',7F8.1,/,(4X,7F8.1))
 9021 FORMAT(I3,1P,1X,7G8.1,/,(4X,7G8.1))
 9030 FORMAT(/1X,'FUNCTIONS    ELECTRON TEMP     ION TEMP      PLASMA ',
     >  '  SMAX   DPERP        ',
     >       /1X,'  OF K       (TARG)  (REF)   (TARG)  (REF)   DENSITY',
     >  '  FACTOR              ')
 9031 FORMAT(1X,(' K =',F6.3,4F8.1,1P,1X,E9.2,0P,F7.2,F9.5,:))
      END
c
c
c
      subroutine pr_toc
      implicit none
      include    'params'
      include    'comtor'
c      include    'cadas'
c      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
      include    'promptdep'
c      include    'slcom'  
      include    'printopt' 
c
c     PR_TOC: This routine adds a section listing to the
c             beginning of the DIVIMP .dat file. In the
c             hypertext version it also includes links to 
c             each section.
c
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new

      new = .TRUE.
c slmod end
      call prb
      call prc('--- TABLE OF CONTENTS ---')
      call prb 

C-----------------------------------------------------------------------


      call prchtml('GEOMETRY GRID AND WALL SECTION',
     >             '0','pr_geom','B')

      call prchtml(' - GEOMETRY GRID AND WALL OPTIONS',
     >             '0','pr_geom_options','0')

      call prchtml(' - GEOMETRY GRID AND WALL INFORMATION',
     >             '0','pr_geom_info','0')

C-----------------------------------------------------------------------
      call prchtml('BACKGROUND PLASMA SECTION','0','pr_bg','B')
      call prchtml(' - BACKGROUND PLASMA OPTIONS',
     >             '0','pr_bg_options','0')
      call prchtml(' - PIN OPTIONS',
     >             '0','pr_pin_options','0')
      call prchtml(' - INPUT PLASMA CONDITIONS','0','pr_plasmain','0')
      call prchtml(' - SUMMARY OF TARGET CONDITIONS',
     >             '0','pr_target','0')

      if (cpinopt.ne.0.or.(cneuta.eq.1.and.ciopte.eq.4)) then  
         call prchtml(' - PIN ITERATION PRINT OUT','0','pr_pinprn','0')
      endif 

C-----------------------------------------------------------------------

      if (ctestsol.eq.0.0) then 

         call prchtml('IMPURITY SIMULATION INPUT SECTION',
     >             '0','pr_sim','B')

         CALL PRChtml (' - ATOMIC PROCESS AND DATA OPTIONS',
     >              '0','pr_proc_options','B')

         call prchtml(' - IMPURITY SIMULATION OPTIONS',
     >                 '0','pr_sim_options','0')

         call prchtml(' - IMPURITY TRANSPORT OPTIONS',
     >                 '0','pr_transport_options','0')

         call prchtml (' - IMPURITY SIMULATION INFORMATION',
     >              '0','pr_siminfo','B')
      endif  

C-----------------------------------------------------------------------

      call prchtml('OTHER OPTIONS LISTING SECTION',
     >             '0','pr_options','B')

      call prchtml(' - MISCELLANEOUS OPTIONS',
     >             '0','pr_misc_options','0')

      call prchtml(' - OPTIONAL INPUT VALUES',
     >             '0','pr_optional_values','0')

!ammod begin
      call prchtml(' - HYDROCARBON OPTIONS',
     >             '0','pr_hydrocarbon_options','0')
!ammod end

C-----------------------------------------------------------------------
     
      if (cprint.eq.1.or.cprint.eq.9) then 

         call prchtml('TRANSPORT COEFFICIENTS EXTRACTOR'//
     >                ' OPTIONS AND RESULTS',
     >             '0','pr_transport_coeff','B')

      endif
c
      if (cprint.eq.1.or.cprint.eq.9) then 

         call prchtml('RECIPROCATING PROBE MODEL RESULTS',
     >             '0','pr_rcp','B')

      endif
c 
      if (ctestsol.eq.0) then 

         call prchtml('DIVIMP IMPURITY SIMULATION INFORMATION',
     >             '0','pr_runtime','B')
         call prc('   - INCLUDING FLUX AND YIELD DATA (IF ANY)') 
c

         call prchtml('DIVIMP SUMMARY OF RESULTS',
     >             '0','pr_summary','B')
         call prchtml(' - NEUTRAL SUMMARY',
     >             '0','pr_neut','0')
         call prchtml(' - ION SUMMARY',
     >             '0','pr_ion','0')
         call prchtml(' - INITIAL IMPURITY NEUTRAL IONIZATION',
     >             '0','pr_ioniz','0')
         call prchtml(' - PARTICLE DEPOSITION/EROSION SUMMARY',
     >             '0','pr_depero','0')
         call prchtml(' - IMPURITY CHARGE STATE SUMMARY',
     >             '0','pr_chargestate','0')  
         CALL PRChtml(' - ALL IONIZATION STATE SUMMARY',
     >              '0','pr_allstates','0')
         if (prompt_depopt.eq.1.or.prompt_depopt.eq.2) then 

            call prchtml(' - PROMPT ION DEPOSITION SUMMARY',
     >                '0','pr_promptdep','0')

         endif 
	
         call prchtml(' - CORE LEAKAGE SUMMARY',
     >                '0','pr_leakage','0')  

         call prchtml(' - ADDITIONAL DIVIMP CALCULATED VALUES',
     >                '0','pr_other','0')
c
c         call prchtml('   - MORE MAIN LEAKAGE INFORMATION',
c     >             '0','pr_addleak','0')
c
         call prchtml(' - PLASMA IMPURITY CONTENT SUMMARY',
     >                '0','pr_content','0')

         call prchtml(' - ESTIMATED DIVERTOR RETENTION',
     >                '0','pr_divrent','0')
         call prchtml(' - RADIATION SUMMARY',
     >                '0','pr_rad','0')

         call prchtml(' - FORCE AND DIFFUSION SUMMARY',
     >                '0','pr_force','0')

      endif
c
      call prchtml('GLOBAL POWER SUMMARY','0','pr_power','0')  
c
      call prchtml('CASE EPILOGUE','0','pr_end','B')

      call prb

      return
      end 
c
c
c
      SUBROUTINE PR_GEOM
      use subgrid_options
      IMPLICIT none
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS GEOMETRY AND TORUS ATTRIBUTES ETC
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
      include    'printopt'
c slmod begin - new
      INTEGER i,j

      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      INTEGER  I,IT,IZ,irlim,ir,in,len2,lenstr,id,ind
      INTEGER   ir,len2,lenstr,ik,id,in
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
      real     tmp_bratio 
c
c     Temporary unit number
c
      integer tmp_unit

      external lenstr
c
      
c      call prc('GEOMETRY AND GRID INFORMATION SECTION')  
      call prb
      call prchtml('--- GEOMETRY AND GRID INFORMATION SECTION ---',
     >             'pr_geom','0','B')
      call prb
c


c
C-----------------------------------------------------------------------
c
c     Print GEOMETRY Options
c
      call pr_geom_options

c
      call prb
      call prchtml('--- GEOMETRY AND GRID INFORMATION ---',
     >             'pr_geom_info','0','B')  
      call prb
C-----------------------------------------------------------------------
      call prb
      call prc ('GRID CHARACTERISTICS:')
      CALL PRI ('  NO OF RINGS                       ', NRS)
      CALL PRI ('  NO OF R PTS FOR RECTANGULAR GRID  ', NXS)
      CALL PRI ('  NO OF Z PTS FOR RECTANGULAR GRID  ', NYS)
      CALL PRI ('  NO OF T PTS FOR RESULTS           ', NTS)
      call prb

      CALL PRB
      CALL PRC  ('TORUS GEOMETRY')
      len2 = lenstr(crun)
      CALL PRC  ('  RUN                  '//CRUN(1:len2))
      CALL PRI  ('  SHOT NUMBER          ', ISHOT)
      CALL PRR  ('  TIME SLICE           ', TSLICE)
      call prb
      CALL PRR2 ('  LOWER LEFT CORNER           (RMIN,ZMIN)', RMIN,ZMIN)
      CALL PRR2 ('  UPPER RIGHT CORNER          (RMAX,ZMAX)', RMAX,ZMAX)
      CALL PRR2 ('  CENTRE                      (RO  ,ZO  )', R0  ,Z0  )
      CALL PRR2 ('  X POINT                     (RXP ,ZXP )', RXP ,ZXP )
      CALL PRR  ('  TOROIDAL FIELD               BPHI      ', CBPHI)
      call prb
c
      if (zxp.gt.z0) then
        call prc('WARNING: THIS GRID HAS THE X-POINT AT THE TOP')
      else
        call prc('WARNING: THIS GRID HAS THE X-POINT AT THE BOTTOM')
      endif
c
      CALL PRC  ('MAGNETIC GEOMETRY - RING NUMBERS')
      CALL PRI  ('  TOTAL NUMBER OF RINGS                  ',NRS)
      CALL PRI2 ('  MAIN PLASMA RINGS                      ',1,IRSEP-1)
c
      if (cgridopt.eq.0.or.cgridopt.eq.1.or.cgridopt.eq.3.or.
     >    cgridopt.eq.LINEAR_GRID.or.cgridopt.eq.GEN_GRID) then
c...    cgridopt=6 for LINEAR_GRID
c..     cgridopt=7 for GEN_GRID
c
      CALL PRI2 ('  SOL RINGS                              ',IRSEP,
     >           IRWALL)
      CALL PRI2 ('  TRAP RINGS                             ',IRTRAP,NRS)
      CALL PRI  ('  SEPARATRIX                             ',IRSEP)
      CALL PRI  ('  WALL RING                              ',IRWALL)
      CALL PRI  ('  TRAP WALL RING                         ',IRTRAP)
      call prc  ('  RING CHARACTERISTICS: ')
c slmod begin - new
      IF (new) THEN
         WRITE(coment,'(2A12)') 'Ring No.','Length (m)'
         CALL PRC(coment)
         DO i = 1, (nrs + 1) / 2
           j = i + (nrs + 1) / 2
           WRITE(coment,1000) i,ksmaxs(i)
           IF (j.LE.nrs)
     .       WRITE(coment(LEN_TRIM(coment)+2:LEN(coment)),1000)
     .         j,ksmaxs(j)
           CALL PRC(coment)
         ENDDO
1000     FORMAT(I12,F12.4)
      ELSE
        do ir = 1,nrs
           write(coment,
     >           '(''   Ring No : '',i6,''  Length : '',f9.5,'' (m)'')')
     >            ir,ksmaxs(ir)
           len2 = lenstr(coment)
           call prc (coment(1:len2))
        end do
      ENDIF
c
c      do ir = 1,nrs
c         write(coment,
c     >         '(''     Ring No : '',i6,''  Length : '',f9.5,'' (m)'')')
c     >          ir,ksmaxs(ir)
c         len2 = lenstr(coment)
c         call prc (coment(1:len2))
c      end do
c slmod end      do ir = 1,nrs
c


c---------------- SUBGRID ------------------------
c
c     Print out the subgrid option if it is active
c
      if (subgrid_opt.eq.0) then 
         call prc('SUBGRID OVERLAY OPTION 0 (*G38) : SUBGRID OFF')
      elseif (subgrid_opt.eq.1) then 
         call prc ('SUBGRID OVERLAY OPTION 1 (*G38) : SUBGRID ON')
         call pri2(' - SUBGRID DIMENSIONS    (*G39) : RDIM,ZDIM:',
     >             sg_rdim,sg_zdim)
         call prr2(' - SUBGRID R-RANGE       (*G40) : RMIN,RMAX:',
     >             sg_rmin,sg_rmax)
         call prr2(' - SUBGRID Z-RANGE       (*G41) : ZMIN,ZMAX:',
     >             sg_zmin,sg_zmax)
      endif

c
c     Calculate and print the average distance for SOL rings from the separatrix
c     - assumes orthogonal grid 
c    
      call calc_average_crossfield_dist
c

      call prb
      CALL PRC  ('WALL CHRACTERISTICS - INDEX NUMBERS:     ')
      CALL PRI2 ('  FIRST AND LAST FOR OUTER WALL          ',
     >           WLWALL1,WLWALL2)
      CALL PRI2 ('  FIRST PLATE                            ',
     >           WLWALL2+1,WLTRAP1-1)
      CALL PRI2 ('  FIRST AND LAST FOR PRIVATE PLASMA WALL ',
     >           WLTRAP1,WLTRAP2)
      CALL PRI2 ('  SECOND PLATE                           ',
     >           WLTRAP2+1,WALLPTS)
c
      elseif (cgridopt.eq.2) then
c
      CALL PRI2 ('  SOL RINGS (INNER)                      ',IRSEP,
     >           IRWALL2)
      CALL PRI2 ('  SOL RINGS (OUTER)                      ',IRSEP2,
     >           IRWALL)
      CALL PRI2 ('  TRAP RINGS (LOWER)                     ',IRTRAP2,
     >          NRS2)
      CALL PRI2 ('  TRAP RINGS (UPPER)                     ',IRTRAP,
     >          NRS)
      CALL PRI2 ('  SEPARATRIX RINGS                       ',IRSEP,
     >          IRSEP2)
      CALL PRI2 ('  WALL RINGS                             ',IRWALL2,
     >           IRWALL)
      CALL PRI2 ('  TRAP WALL RINGS                        ',IRTRAP2,
     >          IRTRAP)
      call prb
      CALL PRC  ('WALL CHRACTERISTICS - INDEX NUMBERS      ')
      CALL PRI2 ('  FIRST AND LAST FOR OUTER WALL          ',
     >           WLWALL1,WLWALL2)
      CALL PRI2 ('  FIRST AND LAST FOR INNER WALL          ',
     >           WLWALL3,WLWALL4)
      CALL PRI2 ('  FIRST PLATE  (Lower Outer)             ',
     >           WLWALL2+1,WLTRAP1-1)
      CALL PRI2 ('  FIRST AND LAST FOR FIRST TRAP WALL     ',
     >           WLTRAP1,WLTRAP2)
      CALL PRI2 ('  SECOND PLATE (Lower Inner)             ',
     >           WLTRAP2+1,WLWALL3-1)
      CALL PRI2 ('  THIRD PLATE  (Upper Inner)             ',
     >           WLWALL4+1,WLTRAP3-1)
      CALL PRI2 ('  FIRST AND LAST FOR SECOND TRAP WALL    ',
     >           WLTRAP3,WLTRAP4)
      CALL PRI2 ('  FOURTH PLATE (Upper Outer)             ',
     >           WLTRAP4+1,WALLPTS)
c
      endif
c
c     Print out separatrix areas 
c
      call prb
      call prc ('SUMMARY OF SEPARATRIX GEOMETRY:')
      call prr ('  POLOIDAL LENGTH OF SEPARATRIX           = ',
     >                                                      asep)   
      call prr ('  TOROIDAL AREA OF SEPARATRIX             = ',
     >                                                  asep_tor)   
      call prr ('  EFFECTIVE POLOIDAL LENGTH OF SEPARATRIX = ',
     >                                                  asep_eff)   
      call prr ('  EFFECTIVE TOROIDAL AREA OF SEPARATRIX   = ',
     >                                              asep_tor_eff)   
 
c
c     Print out target geometry information  
c
      call prb
      CALL PRC  ('TARGET ELEMENT GEOMETRY:')
      call prb
c
      write(coment,'(1x,'' ID '','' IK '','' IR '',6x,''R'',3x,'//
     >            '6x,''Z'',3x,5x,''PSI'',2x,3x,''LENGTH'','//
     >            '6x,''Bth/B'',3x,''SEP DIST'',3x,''MID DIST'')')
c
      call prc(coment)
c
      do id = 1,nds
c
c        ID, IK, IR, R, Z, PSI, LENGTH 
c
         if (ikds(id).gt.nks(irds(id))/2) then 
            in = 1
         else
            in = 2
         endif 
c
         if (kbfst(irds(id),in).ne.0.0) then
            tmp_bratio = 1.0 / kbfst(irds(id),in)
         else
            tmp_bratio = 0.0
         endif
c
         write(coment,
     >       '(3i4,2(1x,f9.5),1x,f9.6,1x,f9.5,2x,g13.5,2(1x,g12.5))')
     >       id,ikds(id),irds(id),rp(id),zp(id),psitarg(irds(id),in),
     >       dds(id),tmp_bratio,sepdist2(id),middist(irds(id),in)
c
         call prc(coment)
c
      end do



c
c     Print out target geometry information to debug file in ring order
c
      
      call find_free_unit_number(tmp_unit)

      open(tmp_unit,file='target_geometry.dat',form='formatted')
     
      write(tmp_unit,*)
      write(tmp_unit,'(a)') 'TARGET ELEMENT GEOMETRY:'
      write(tmp_unit,*)
c
c
      do in = 1,2

         if (in.eq.1) then 
            write(tmp_unit,'(a)') INNER//' TARGET:'
            
         elseif (in.eq.2) then 
            write(tmp_unit,'(a)') OUTER//' TARGET:'

         endif


         write(tmp_unit,'(1x,'' ID '','' IK '','' IR '',6x,''R'',3x,'//
     >            '6x,''Z'',3x,5x,''PSI'',2x,3x,''LENGTH'','//
     >            '6x,''Bth/B'',3x,''SEP DIST'',3x,''MID DIST'')')


         do ir = irsep,nrs
            id = idds(ir,in)

            if (kbfst(irds(id),in).ne.0.0) then
               tmp_bratio = 1.0 / kbfst(irds(id),in)
            else
               tmp_bratio = 0.0
            endif

            write(tmp_unit,
     >       '(3i4,2(1x,f9.5),1x,f9.6,1x,f9.5,2x,g13.5,2(1x,g12.5))')
     >       id,ikds(id),irds(id),rp(id),zp(id),psitarg(irds(id),in),
     >       dds(id),tmp_bratio,sepdist2(id),middist(irds(id),in)

         end do

      end do
      write(tmp_unit,*)

      close(tmp_unit)
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_BG  (NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS INFORMATION ABOUT THE CALCULATED 
c     BACKGROUND PLASMA
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
      include    'cadas'
      include    'cgeom'
      include    'cioniz'
      include    'cedge2d'
      include    'dynam4'
      include    'dynam5'
      include    'pindata'
      include    'adpak_com'
      include    'promptdep'
      include    'slcom'  
      include    'printopt'
c
      CHARACTER prtype*4
      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
      real     totmhpin
      external lenstr
c slmod begin - new
      INTEGER j

      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
      CALL PRB
c
c      CALL PRC  ('BACKGROUND PLASMA CHARACTERISTICS:')
c
      call prchtml('--- BACKGROUND PLASMA SECTION ---',
     >                'pr_bg','0','B')
      call prb 

      CALL PRR  ('  PLASMA ION MASS              MB           ', CRMB)
      CALL PRI  ('  PLASMA ION CHARGE            ZB           ', CIZB)
      call prb
c
C-----------------------------------------------------------------------
c
c     Print BACKGROUND PLASMA Options
c
      call pr_bg_options(NIZS,NIMPS,NIMPS2,nymfs)
c
C-----------------------------------------------------------------------
c
c     Print PIN Options
c
      call pr_pin_options 
c
      call prb
c
c     Electron Temperatures
c
      call prchtml('--- INPUT PLASMA CONDITIONS ---',
     >                'pr_plasmain','0','B')
c
C-----------------------------------------------------------------------
      IF (CIOPTG.NE.99.and.cioptg.ne.98.and.
     >    (.not.((cre2d.eq.1.or.cre2d.eq.2).and.
     >            e2dtargopt.ne.0.0))) THEN
       if (Cioptg.eq.5.or.cioptg.eq.6) then
        CALL PRR ('  ELECTRON TEMPERATURE AT SEPARATRIX (EV)   ', CTEBP)
       else
        CALL PRR ('  TEMPERATURE OF ELECTRONS     TEB0  (EV)   ', CTEB0)
       endif
       IF( CIOPTK.EQ.2 .and.cioptf.lt.12)
     > CALL PRR ('  TEMP OF ELECTRONS AT PLATES  TEBP  (EV)   ', CTEBP)
       IF ((CIOPTK.EQ.0.OR.CIOPTK.EQ.1).and.cioptf.lt.12) THEN
       CALL PRR ('    GRADIENT DISTANCE FACTOR   FEBL1        ', CFEBL1)
       CALL PRR ('    GRADIENT DISTANCE FACTOR   FEBL2        ', CFEBL2)
       CALL PRR ('    GRADIENT TARGET FRACTION   FEBT         ', CFEBT)
       CALL PRR ('    GRADIENT TARGET FRACTION   FEB2         ', CFEB2)
       endif
c
       IF (CIOPTG.EQ.0)
     > CALL PRR ('    OUTER TEB STEP             TEBOUT(EV)   ', CTEBOU)
c
       IF (CIOPTG.EQ.1.or.cioptg.eq.5.or.cioptg.eq.6)
     > CALL PRR ('    SOL TARGET TEB EXP DECAY FACTOR     (M)', CTEBOU)
c
       IF (cioptg.eq.6)
     > CALL PRR ('    PRIVATE PLASMA TEB EXP DECAY FACTOR (M) ',CTEBOUP)
c
       if (cioptg.ne.4.and.cioptg.ne.6)
     >  CALL PRR ('    PRIVATE PLASMA CONSTANT    TEBT  (EV)   ', CTEBT)
c
       if (ccoreopt.eq.0.or.ccoreopt.eq.2.or.ccoreopt.eq.3.or.
     >     ccoreopt.eq.4) then
       CALL PRR ('    INNER TEB STEP             TEBIN (EV)   ', CTEBIN)
       endif
c

       IF (CIOPTG.EQ.2) THEN
         CALL PRC ('    TEB SPECIFIED BY INPUT DATA/RING')
         call prr ('    TEB INPUT MULTIPLIER:  ',te_mult_i)
c slmod begin - new
         IF (new) THEN

         ELSE
           DO 100 I=1,NLPDATI
             CALL PRR2('       RING, TEB   : ',LPDATI(I,1),LPDATI(I,2))
 100       CONTINUE
         ENDIF
c
c           DO 100 I=1,NLPDATI
c             CALL PRR2('       RING, TEB   : ',LPDATI(I,1),LPDATI(I,2))
c 100       CONTINUE
c slmod end
       ENDIF


       IF (CIOPTG.EQ.3.OR.CIOPTG.EQ.4
     >     .or.cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) THEN
c
         CALL PRC ('    TEB SPECIFIED BY INPUT DATA/RING '
     >                  //INner//'/'//OUter)
         call prr ('    TEB '//INner//' INPUT MULTIPLIER:  ',te_mult_i)
c slmod begin - new
         IF (new) THEN

         ELSE
           CALL PRC ('    '//INNER//' PLATE DATA: ')
           DO 200 I=1,NLPDATI
             CALL PRR2('       RING, TEB: ',LPDATI(I,1),LPDATI(I,2))
 200       CONTINUE
         ENDIF
c
c         CALL PRC ('    '//INner//' PLATE DATA: ')
c
c         DO 200 I=1,NLPDATI
c           CALL PRR2('       RING, TEB: '
c     >         ,LPDATI(I,1),LPDATI(I,2))
c 200     CONTINUE
c slmod end
c
         call prr ('    TEB '//OUTER//' INPUT MULTIPLIER:  ',te_mult_o)
c slmod begin - new
         IF (new) THEN

         ELSE
           CALL PRC ('    '//OUTER//' PLATE DATA: ')
           DO 300 I=1,NLPDATO
             CALL PRR2('       RING, TEB: '
     >           ,LPDATO(I,1),LPDATO(I,2))
 300       CONTINUE
         ENDIF
c
c         CALL PRC ('    '//OUTER//' PLATE DATA: ')
cc
c         DO 300 I=1,NLPDATO
c           CALL PRR2('       RING, TEB: '
c     >         ,LPDATO(I,1),LPDATO(I,2))
c 300     CONTINUE
c slmod end
       ENDIF


c
      if (ctegcut.gt.0.0) then
c
         if (cflatopt.eq.1) then

            call prr('* ELECTRON TEMPERATURE FORCED FLAT FOR'//
     >                 ' S OR (SMAX-S) > SMAX * ',
     >             ctegcut)
c
         elseif (cflatopt.eq.2) then
c
            call prr('* ELECTRON TEMPERATURE FORCED FLAT FOR'//
     >                 ' S  > SMAX * ',
     >             ctegcut)
            call prc('  THE VALUE OF Te AT THIS POSITION IS THE'//
     >               ' MAXIMUM Te')
            call prc('  ALLOWED ON THE RING.')
c
         endif
c
      endif
C
C-----------------------------------------------------------------------
c
c      Ion temperatures
c
c
       if (Cioptg.eq.5.or.cioptg.eq.6) then
        CALL PRR ('  ION TEMPERATURE AT SEPARATRIX      (EV)   ', CTIBP)
       else
        CALL PRR ('  TEMPERATURE OF IONS          TIB0  (EV)   ', CTIB0)
       endif
       IF( CIOPTL.EQ.2 .and.cioptf.lt.12)
     > CALL PRR ('  TEMP OF IONS AT PLATES       TIBP  (EV)   ', CTIBP)
       IF ((CIOPTL.EQ.0.OR.CIOPTL.EQ.1).and.cioptf.lt.12) THEN
       CALL PRR ('    GRADIENT DISTANCE FACTOR   FIBL1        ', CFIBL1)
       CALL PRR ('    GRADIENT DISTANCE FACTOR   FIBL2        ', CFIBL2)
       CALL PRR ('    GRADIENT TARGET FRACTION   FIBT         ', CFIBT)
       CALL PRR ('    GRADIENT TARGET FRACTION   FIB2         ', CFIB2)
       endif
c
       IF (CIOPTG.EQ.0)
     > CALL PRR ('    OUTER TIB STEP             TIBOUT(EV)   ', CTIBOU)
c
       IF (CIOPTG.EQ.1.or.cioptg.eq.5.or.cioptg.eq.6)
     > CALL PRR ('    SOL TARGET TIB EXP DECAY FACTOR  (M)    ', CTIBOU)
c
       IF (cioptg.eq.6)
     > CALL PRR ('    PRIVATE PLASMA TIB EXP DECAY FACTOR (M) ',CTIBOUP)
c
       if (cioptg.ne.4.and.cioptg.ne.6)
     >  CALL PRR ('    PRIVATE PLASMA CONSTANT    TIBT  (EV)   ', CTIBT)
c
       if (ccoreopt.eq.0.or.ccoreopt.eq.2.or.ccoreopt.eq.3.or.
     >     ccoreopt.eq.4) then
       CALL PRR ('    INNER TIB STEP             TIBIN (EV)   ', CTIBIN)
       endif
c
       IF (CIOPTG.EQ.2) THEN
         CALL PRC ('    TIB SPECIFIED BY INPUT DATA/RING')
         call prr ('    TIB INPUT MULTIPLIER:  ',ti_mult_i)
c slmod begin - new
         IF (new) THEN

         ELSE
           DO 110 I=1,NLPDATI
             CALL PRR2('       RING, TIB   : ',LPDATI(I,1),LPDATI(I,3))
 110       CONTINUE
         ENDIF
c
c         DO 110 I=1,NLPDATI
c           CALL PRR2('       RING, TIB   : ',LPDATI(I,1),LPDATI(I,3))
c 110     CONTINUE
c slmod end
       ENDIF
       IF (CIOPTG.EQ.3.OR.CIOPTG.EQ.4
     >     .or.cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) THEN
c
         CALL PRC ('    TIB SPECIFIED BY INPUT DATA/RING '
     >                //INNER//'/'//OUTER)
         call prr ('    TIB '//INNER//' INPUT MULTIPLIER:  ',ti_mult_i)
c slmod begin - new
         IF (new) THEN
           call prr ('    TIB '//OUTER//' INPUT MULTIPLIER: ',ti_mult_o)
         ELSE
           CALL PRC ('    '//INNER//' PLATE DATA: ')
           DO 210 I=1,NLPDATI
             CALL PRR2('       RING, TIB: '
     >           ,LPDATI(I,1),LPDATI(I,3))
 210       CONTINUE
           call prr ('    TIB '//OUTER//' INPUT MULTIPLIER: ',ti_mult_o)
           CALL PRC ('    '//OUTER//' PLATE DATA: ')
           DO 310 I=1,NLPDATO
             CALL PRR2('       RING, TIB: '
     >           ,LPDATO(I,1),LPDATO(I,3))
 310       CONTINUE
         ENDIF
c
c         CALL PRC ('    '//INNER//' PLATE DATA: ')
cc
c         DO 210 I=1,NLPDATI
c           CALL PRR2('       RING, TIB: '
c     >         ,LPDATI(I,1),LPDATI(I,3))
c 210     CONTINUE
cc
c         call prr ('    TIB '//OUTER//' INPUT MULTIPLIER:  ',ti_mult_o)
c         CALL PRC ('    '//OUTER//' PLATE DATA: ')
c         DO 310 I=1,NLPDATO
c           CALL PRR2('       RING, TIB: '
c     >         ,LPDATO(I,1),LPDATO(I,3))
c 310     CONTINUE
c slmod end
c
       ENDIF
c
      if (ctigcut.gt.0.0) then
c
         if (cflatopt.eq.1) then

            call prr('* ION TEMPERATURE      FORCED FLAT FOR'//
     >                 ' S OR (SMAX-S) > SMAX * ',
     >             ctigcut)
         elseif (cflatopt.eq.2) then
c
            call prr('* ION TEMPERATURE      FORCED FLAT FOR'//
     >                 ' S  > SMAX * ',
     >             ctegcut)
            call prc('  THE VALUE OF Ti AT THIS POSITION IS THE'//
     >               ' MAXIMUM Ti')
            call prc('  ALLOWED ON THE RING.')
c
         endif
c
      endif
c
      if (ctegcut.gt.0.0.or.ctigcut.gt.0.0.and.cflatopt.ne.0) call prb

c
C-----------------------------------------------------------------------
C
c      Densities
c
       if (Cioptg.eq.5.or.cioptg.eq.6) then
        CALL PRR ('  DENSITY AT SEPARATRIX TARGET NEBP (M**-3)', CNEBP)
       elseif (cioptg.eq.0.or.cioptg.eq.1) then
        CALL PRR ('  DENSITY AT 0                 NB0   (M**-3)', CNB0)
       endif
c
       IF (CIOPTG.EQ.0)
     > CALL PRR ('    OUTER NB STEP              NBOUT (M**-3)', CNBOUT)
c
       IF (CIOPTG.EQ.1.or.cioptg.eq.5.or.cioptg.eq.6)
     > CALL PRR ('    SOL TARGET NB EXP DECAY FACTOR      (M) ', CNBOUT)
c
       IF (cioptg.eq.6)
     > CALL PRR ('    PRIVATE PLASMA NB EXP DECAY FACTOR  (M) ',CNBOUP)
c
       if (cioptg.ne.4.and.cioptg.ne.6)
     >  CALL PRR ('    PRIVATE PLASMA CONSTANT     NBT   (M**-3)', CNBT)
c
       if (ccoreopt.eq.0.or.ccoreopt.eq.2.or.ccoreopt.eq.3.or.
     >     ccoreopt.eq.4) then
       CALL PRR ('    INNER NB STEP              NBIN  (M**-3)', CNBIN)
       endif
c
c      Set type of quantity being printed from the LPDAT arrays
c
       if (lpdatsw.eq.0) then
          prtype = 'NB  '
       elseif (lpdatsw.eq.1) then
          prtype = 'Isat'
       endif

c
       IF (CIOPTG.EQ.2) THEN
c
         if (lpdatsw.eq.0) then
            CALL PRC ('    '//prtype//'SPECIFIED BY INPUT DATA/RING')
         elseif (lpdatsw.eq.1) then
            CALL PRC ('    '//prtype//'SPECIFIED BY INPUT DATA/RING'//
     >                ' - CONVERTED TO NB')
         endif
c
         call prr ('    NB INPUT MULTIPLIER:  ',n_mult_i)
c
c slmod begin - new
         IF (new) THEN

         ELSE
           DO 120 I=1,NLPDATI
             if (lpdatsw.eq.0) then
                CALL PRR2('       RING, '//prtype//
     >                 ' : ',LPDATI(I,1),LPDATI(I,4))
             elseif (lpdatsw.eq.1) then
                CALL PRR3('       RING, '//prtype//
     >                 ' : ',LPDATI(I,1),LPDATI(I,4),
     >                 knds(idds(int(lpdati(i,1)),1)))
             endif
 120       CONTINUE
         ENDIF
c
c         DO 120 I=1,NLPDATI
c           if (lpdatsw.eq.0) then
c              CALL PRR2('       RING, '//prtype//
c     >               ' : ',LPDATI(I,1),LPDATI(I,4))
c           elseif (lpdatsw.eq.1) then
c              CALL PRR3('       RING, '//prtype//
c     >               ' : ',LPDATI(I,1),LPDATI(I,4),
c     >               knds(idds(int(lpdati(i,1)),1)))
c           endif
c 120     CONTINUE
c slmod end
       ENDIF
c
       IF (CIOPTG.EQ.3.OR.CIOPTG.EQ.4
     >     .or.cioptg.eq.90.or.cioptg.eq.91.or.cioptg.eq.92) THEN
c
c
c
         if (lpdatsw.eq.0) then
            CALL PRC ('    '//prtype//
     >             ' SPECIFIED BY INPUT DATA/RING '
     >            //INNER//'/'//OUTER)
         elseif (lpdatsw.eq.1) then
            CALL PRC ('    '//prtype//
     >             ' SPECIFIED BY INPUT DATA/RING '
     >             // INNER//'/'//OUTER//
     >             ' - CONVERTED TO NB')
         endif
c
         call prr ('    NB '//INNER//' INPUT MULTIPLIER:  ',n_mult_i)
c slmod begin - new
         IF (new) THEN

         ELSE
           CALL PRC ('    '//INNER//' PLATE DATA: ')
           DO 220 I=1,NLPDATI
             if (lpdatsw.eq.0) then
                CALL PRR2('       RING, '//prtype//': '
     >              ,LPDATI(I,1),LPDATI(I,4))
             elseif (lpdatsw.eq.1) then
                CALL PRR3('       RING, '//prtype//': '
     >              ,LPDATI(I,1),LPDATI(I,4),
     >              knds(idds(int(lpdati(i,1)),1)))
             endif
 220       CONTINUE
         ENDIF
c
c         CALL PRC ('    '//INNER//' PLATE DATA: ')
c         DO 220 I=1,NLPDATI
c
c           if (lpdatsw.eq.0) then
c              CALL PRR2('       RING, '//prtype//': '
c     >            ,LPDATI(I,1),LPDATI(I,4))
c           elseif (lpdatsw.eq.1) then
c              CALL PRR3('       RING, '//prtype//': '
c     >            ,LPDATI(I,1),LPDATI(I,4),
c     >            knds(idds(int(lpdati(i,1)),1)))
c           endif
cc
c 220     CONTINUE
c slmod end         CALL PRC ('    '//INNER//' PLATE DATA: ')
c
         call prr ('    NB '//OUTER//' INPUT MULTIPLIER:  ',n_mult_o)
c slmod begin - new
         IF (new) THEN

         ELSE
           CALL PRC ('    '//OUTER//' PLATE DATA: ')
           DO 320 I=1,NLPDATO
             if (lpdatsw.eq.0) then
                CALL PRR2('       RING, '//prtype//': '
     >              ,LPDATO(I,1),LPDATO(I,4))
             elseif (lpdatsw.eq.1) then
                CALL PRR3('       RING, '//prtype//': '
     >              ,LPDATO(I,1),LPDATO(I,4),
     >              knds(idds(int(lpdato(i,1)),2)))
             endif
 320       CONTINUE
         ENDIF
c
c         CALL PRC ('    '//OUTER//' PLATE DATA: ')
cc
c         DO 320 I=1,NLPDATO
c
c           if (lpdatsw.eq.0) then
c              CALL PRR2('       RING, '//prtype//': '
c     >            ,LPDATO(I,1),LPDATO(I,4))
c           elseif (lpdatsw.eq.1) then
c              CALL PRR3('       RING, '//prtype//': '
c     >            ,LPDATO(I,1),LPDATO(I,4),
c     >            knds(idds(int(lpdato(i,1)),2)))
c           endif
cc
c 320     CONTINUE
c slmod end
      ENDIF
c slmod begin - new
      IF (new.AND.(cioptg.EQ. 2.OR.cioptg.EQ. 3.OR.cioptg.EQ. 4.OR.
     .             cioptg.EQ.90.OR.cioptg.EQ.91.OR.cioptg.EQ.92)) THEN

        CALL PRC ('  '//inner//' PLATE DATA:                   '//
     .            '  '//outer//' PLATE DATA: ')

        IF     (lpdatsw.EQ.0) THEN
          WRITE(coment,321)
        ELSEIF (lpdatsw.eq.1) THEN
          WRITE(coment,322)
        ENDIF
        CALL PRC(coment)

        DO i = 1, MAX(nlpdati,nlpdato)
          IF (i.LE.nlpdati)
     .      WRITE(coment                                ,323)
     .        INT(lpdati(i,1)),(lpdati(i,j),j=2,4)
          IF (i.LE.nlpdato)
     .      WRITE(coment(LEN_TRIM(coment)+3:LEN(coment)),323)
     .        INT(lpdato(i,1)),(lpdato(i,j),j=2,4)
          CALL PRC(coment)
        ENDDO

321     FORMAT(2(2X,'Ring  Te(eV)  Ti(eV)       Ne(m-2)  '))
322     FORMAT(2(2X,'Ring  Te(eV)  Ti(eV) Jsat(s-1 m-2)  '))
323     FORMAT(2X,I4,2F8.2,1P,E14.2,0P)
      ENDIF

      CALL SLOPT01(7)
c slmod endc
C-----------------------------------------------------------------------
c
c     Print out plasma characteristics across target for plasma decay
c     option 5
c
      if (cioptg.eq.5.or.cioptg.eq.6) then
c
       if (cioptg.eq.5) then
          irlim = irwall
       elseif (cioptg.eq.6) then
          irlim = nrs
       endif
c
       call prb
       call prchtml('  PLASMA VARIATION ACROSS TARGET:',
     >              'pr_target','0','B')
c
c       CALL PRC ('  PLASMA VARIATION ACROSS TARGET:')
c
       call prr2('  '//INNER//' TARGET  (R,Z)  = ',rp(idds(irsep,1)),
     >                               zp(idds(irsep,1)))
       call prc ('  RING    DIST '//
     >           '       Ne (m**3)     Te (eV)    Ti (eV)')
c
       do i = irsep,irlim
          write (coment,'(3x,i4,1x,f8.5,1x,g15.6,1x,2f10.4)') i,
     >           sepdist(idds(i,1)),knds(idds(i,1)),kteds(idds(i,1)),
     >           ktids(idds(i,1))
          call prc(coment)
       end do
c
       call prr2('  '//OUTER//' TARGET  (R,Z)  =  ',rp(idds(irsep,2)),
     >                              zp(idds(irsep,2)))
       call prc ('  RING    DIST    Ne (m**3)     Te (eV)    Ti (eV)')
c
       do i = irsep,irlim
          write (coment,'(3x,i4,1x,f8.5,1x,g15.6,1x,2f10.4)') i,
     >           sepdist(idds(i,2)),knds(idds(i,2)),kteds(idds(i,2)),
     >           ktids(idds(i,2))
          call prc(coment)
       end do
c
c     For cases where cioptg not equal to 5 or 6
c
      else
c
c      Print out the target conditions that have
c      been extracted from the plasma file or calculated from inputs.

       call prb 
       call prchtml(' SUMMARY OF ACTUAL TARGET CONDITIONS:',
     >              'pr_target','0','B')
c
c       call prc(' SUMMARY OF ACTUAL TARGET CONDITIONS:')
c
       call prr2('   '//INNER//' TARGET (R,Z) = ',rp(idds(irsep,1)),
     >                              zp(idds(irsep,1)))
       call prc('   RING   Ne (m**3)   Te (eV)  TeCRIT(eV) Ti (eV)'
     >           //'     Vb         Cs        Isat')
       do i = irsep,nrs
          write (coment,
     >           '(3x,i4,1x,1p,g12.4,1x,0p,
     >             3(1x,f7.3),1p,3(1x,g10.3))')
     >           i,
     >           knds(idds(i,1)),kteds(idds(i,1)),
     >           (1.0e-18*knds(idds(i,1))*ksmaxs(i)/2.0
     >            *sqrt(ktids(idds(i,1))))**0.4,
     >           ktids(idds(i,1)),kvds(idds(i,1)),
     >         SQRT(0.5*EMI*(KTEdS(idds(i,1))+KTIdS(idds(i,1)))
     >         *(1+RIZB)/CRMB),
     >            knds(idds(i,1))*kvds(idds(i,1))*ech
          call prc(coment)
       end do
c
       if (nsheath_vali.gt.0) then 
          call prb 
          call prc('SPECIFIED '//INNER//' TARGET SHEATH TEMPERATURES:')
          do i = 1,nsheath_vali
             write(coment,'(a,1x,i5,6x,a,1x,f10.2)') 'RING :',
     >         int(sheath_vali(i,1)), 'Te:',sheath_vali(i,2)
             call prc(coment)
          end do
       endif
 

       call prr2('   '//OUTER//' TARGET (R,Z) = ',rp(idds(irsep,2)),
     >                              zp(idds(irsep,2)))
       call prc('   RING   Ne (m**3)   Te (eV)  TeCRIT(eV) Ti (eV)'
     >           //'     Vb         Cs      Isat')
c
c       call prc('   RING       Ne (m**3)     Te (eV)    Ti (eV)'
c     >           //'    Vb              Cs')
c
       do i = irsep,nrs
          write (coment,
     >           '(3x,i4,1x,1p,g12.4,1x,0p,
     >             3(1x,f7.3),1p,3(1x,g10.3))')
c
c     >           '(3x,i4,1x,1p,g15.6,0p,1x,1p,
c     >             2g12.4,0p,1x,g13.5,1x,g13.5)')
c
     >           i,
     >           knds(idds(i,2)),kteds(idds(i,2)),
     >           (1.0e-18*knds(idds(i,2))*ksmaxs(i)/2.0
     >            *sqrt(ktids(idds(i,2))))**0.4,
     >           ktids(idds(i,2)),kvds(idds(i,2)),
     >         SQRT(0.5*EMI*(KTEdS(idds(i,2))+KTIdS(idds(i,2)))
     >         *(1+RIZB)/CRMB),
     >              knds(idds(i,2))*kvds(idds(i,2))*ech
          call prc(coment)
       end do
c
       if (nsheath_valo.gt.0) then 
          call prb
          call prc('SPECIFIED '//OUTER//' TARGET SHEATH TEMPERATURES:')
          do i = 1,nsheath_valo
             write(coment,'(a,1x,i5,6x,a,1x,f10.2)') 'RING :',
     >         int(sheath_valo(i,1)), 'Te:',sheath_valo(i,2)
             call prc(coment)
          end do
       endif
c
       call prb
c
      endif
c
c     For EDGE2D background plasmas on JET grids - print out the
c     target conditions - as read in or calculated. Or for cases that
c     use the EDGE2D target conditions.
c
      else
c
c
c      Print out the target conditions that have
c      been extracted from the plasma file or calculated from inputs.
c
       call prb  
       call prchtml(' SUMMARY OF TARGET CONDITIONS:',
     >              'pr_target','0','B')
c
c       call prc(' SUMMARY OF TARGET CONDITIONS:')
c
       call prr2('   '//INNER//' TARGET (R,Z) = ',rp(idds(irsep,1)),
     >                              zp(idds(irsep,1)))
c
c  IPP/01 - Krieger - changed formats for table and added R,Z 
c
       call prc('  RING  R       Z      Ne(m**3)    Te(eV)  Ti(eV)'
     >           //' Vb        Cs        Isat')
c
       do i = irsep,nrs
          write (coment,
     >           '(2x,i4,1x,f6.3,2x,f6.3,1x,1p,g10.3,0p,1x,
     >             f7.3,1x,f7.3,1p,3(1x,g9.2))')
     >           i,rp(idds(i,1)),zp(idds(i,1)),
     >           knds(idds(i,1)),kteds(idds(i,1)),
     >           ktids(idds(i,1)),kvds(idds(i,1)),
     >         SQRT(0.5*EMI*(KTEdS(idds(i,1))+KTIdS(idds(i,1)))
     >         *(1+RIZB)/CRMB),
     >            knds(idds(i,1))*kvds(idds(i,1))*ech
          call prc(coment)
       end do
c
c       call prc('   RING       Ne (m**3)     Te (eV)    Ti (eV)'
c     >           //'    Vb        Cs       Isat ')
c
c       do i = irsep,nrs
c          write (coment,
c     >           '(3x,i4,1x,1p,g12.4,0p,1x,1p,
c     >             2g12.4,3(1x,g10.3))')
c     >           i,
c     >           knds(idds(i,1)),kteds(idds(i,1)),
c     >           ktids(idds(i,1)),kvds(idds(i,1)),
c     >         SQRT(0.5*EMI*(KTEdS(idds(i,1))+KTIdS(idds(i,1)))
c     >         *(1+RIZB)/CRMB),
c     >            knds(idds(i,1))*kvds(idds(i,1))*ech
c          call prc(coment)
c       end do
c
       call prr2('   '//OUTER//' TARGET (R,Z) = ',rp(idds(irsep,2)),
     >                              zp(idds(irsep,2)))

c
c  IPP/01 - Krieger - changed formats for table and added R,Z 
c

       call prc('  RING  R       Z      Ne(m**3)    Te(eV)  Ti(eV)'
     >           //' Vb        Cs        Isat')
c
       do i = irsep,nrs
          write (coment,
     >           '(2x,i4,1x,f6.3,2x,f6.3,1x,1p,g10.3,0p,1x,
     >             f7.3,1x,f7.3,1p,3(1x,g9.2))')
     >           i,rp(idds(i,2)),zp(idds(i,2)),
     >           knds(idds(i,2)),kteds(idds(i,2)),
     >           ktids(idds(i,2)),kvds(idds(i,2)),
     >         SQRT(0.5*EMI*(KTEdS(idds(i,2))+KTIdS(idds(i,2)))
     >         *(1+RIZB)/CRMB),
     >            knds(idds(i,2))*kvds(idds(i,2))*ech
          call prc(coment)
       end do
c
c       call prc('   RING       Ne (m**3)     Te (eV)    Ti (eV)'
c     >           //'    Vb        Cs       Isat')
c
c       do i = irsep,nrs
c          write (coment,
c     >           '(3x,i4,1x,1p,g12.4,0p,1x,1p,
c     >             2g12.4,3(1x,g10.3))')
c     >           i,
c     >           knds(idds(i,2)),kteds(idds(i,2)),
c     >           ktids(idds(i,2)),kvds(idds(i,2)),
c     >         SQRT(0.5*EMI*(KTEdS(idds(i,2))+KTIdS(idds(i,2)))
c     >         *(1+RIZB)/CRMB),
c     >            knds(idds(i,2))*kvds(idds(i,2))*ech
c          call prc(coment)
c       end do
c
       call prb
c
c     Endif from cioptg.ne.99 -
c
      endif
c
c     Print out the target heat flux summary (MOVED)
c
      call pr_targfluxdata
c
C-----------------------------------------------------------------------
c
      if (cpinopt.ne.0.or.(cneuta.eq.1.and.ciopte.eq.4)) then 
         call prb
         call prchtml('PIN ITERATION SUMMARIES:',
     >              'pr_pinprn','0','B')
c
c         call prc('PIN ITERATION SUMMARIES:')
c
c        Insert the PIN iteration summary information at this point in
c        the .dat file
c
         call appendfile (datunit, pinunit) 
         call prb
c
c        Print out summary of Johnson-Hinov factor information  
c       
         call pr_jhfact
c
      endif
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_SIM(NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS OUT SIMULATION RELATED QUANTITIES NOT
c     COVERED AS OPTIONS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
c      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
      include    'dynam4'
c      include    'dynam5'
      include    'pindata'
c      include    'adpak_com'
      include    'promptdep'
      include    'slcom'  
      include    'printopt'
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
      integer in 
c
c      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
c      external lenstr
c
c      call prc ('IMPURITY SIMULATION CHARACTERISTICS:')
      call prb 
      call prchtml ('--- IMPURITY SIMULATION SECTION ---',
     >              'pr_sim','0','B')
      call prb 
c
C-----------------------------------------------------------------------
c
c     Print ATOMIC PROCESS AND DATA OPTIONS Options
c
      call pr_procdata_options(NIZS)
c
C-----------------------------------------------------------------------
c
c     Print SIMULATION Options
c
      call pr_sim_options(NIZS,NIMPS,NIMPS2,nymfs)
c
C-----------------------------------------------------------------------
c
c     Print SIMULATION TRANSPORT Options
c
      call pr_transport_options
c

      call prb 
      call prchtml ('--- IMPURITY SIMULATION INFORMATION ---',
     >              'pr_siminfo','0','B')
      call prb 
c
C-----------------------------------------------------------------------
      IF     (IMODE.EQ.1) THEN
        CALL PRC ('  OPERATION MODE  TIME DEPENDENT')
      ELSEIF (IMODE.EQ.2) THEN
        CALL PRC ('  OPERATION MODE  STEADY STATE')
      ENDIF

      CALL PRR  ('  IMPURITY ION MASS            MI           ', CRMI)
      CALL PRI  ('  IMPURITY ATOMIC NUMBER       ZI           ', CION)


      CALL PRI ('  MAXIMUM IONIZATION STATE          ', NIZS)
      CALL PRI ('  NO OF IMPURITY PARTICLES TO FOLLOW     ', NIMPS)
      IF (NIMPS2.GT.0.and.cneuta.eq.0)
     >  CALL PRI ('  NO OF SUPPLEMENTARY PARTICLES TO FOLLOW', NIMPS2)
      if (cneutd.eq.6.and.nimps2.gt.0.0) then
        call prc ('  NOTE : FOR THE PHYSICAL+CHEMICAL SPUTTERING OPTION'
     >)
        call pri ('  CHOSEN THE TOTAL NUMBER OF PARTICLES IS  : ',
     >  nimps+nimps2)
        call prc ('  DISTRIBUTED APPROPRIATELY BETWEEN THE TWO SPUTTER')
        call prc ('  MECHANISMS. ')
      endif

      CALL PRR ('  TIMESTEP FOR ATOMS (S)            ', FSRATE)
      CALL PRR ('  TIMESTEP FOR IONS  (S)            ', QTIM)
c
      if (nabsfac.gt.0.0.and.cioptf.ne.21) then
      CALL PRR ('* ABSOLUTE FACTOR EXTERNALLY IMPOSED',nabsfac)
      elseif (nabsfac.gt.0.0.and.cioptf.eq.21) then
      CALL PRR ('  SPECIFIED PLATE POWER FLUX        ',nabsfac)
      endif
c
      CALL PRR ('  ION REMOVAL LOSS TIME (S)         ', TLOSS)
      call prr ('  MAXIMUM ION DWELL TIME (S)        ', CSTMAX*qtim)
      WRITE (7,'(1X,''  DIVIMP RANDOM NUMBER SEED'',10X,I15)') CISEED
      WRITE (7,'(1X,''  PIN/NIMBUS RANDOM NUMBER SEED'',10X,I15)')
     >              piniseed
c
      if (ircore.ne.1) then
         call prb
         call pri('* WARNING: THE CORE MIRROR FOR IMPURITY TRANSPORT'//
     >             ' HAS BEEN MOVED TO: ',ircore)
         call prb
      endif
c

      IF (CIZSET.GE.1.AND.CIZSET.LE.NIZS)
     >CALL PRI  ('  SET TI=TB FOR IONS REACHING STATE         ', CIZSET)
C-----------------------------------------------------------------------

      CALL PRB
      CALL PRC  ('SOURCE CHARACTERISTICS')
      IF     (CNEUTA.EQ.0) THEN
       CIZSC = 1
       CALL PRC ('  INITIAL ION TEMPERATURE  (EV)     PASSED BACK FROM N
     >EUT')
       IF (CNEUTB.EQ.0.or.cneutb.eq.3) THEN
        CALL PRC ('  SOURCE LOCATION FOR NEUT          DISTRIBUTED ALONG
     > TARGET ')
       ELSEIF (CNEUTB.EQ.1) THEN
        CALL PRR ('  SOURCE R POSN FOR NEUT   (M)     ', CXSC)
        CALL PRR ('  SOURCE Z POSN FOR NEUT   (M)     ', CYSC)
       ELSEIF (CNEUTB.EQ.2.or.cneutb.eq.4) THEN
        CALL PRC ('  SOURCE LOCATION FOR NEUT          DISTRIBUTED ALONG
     > WALLS ')
       elseif (cneutb.eq.6) then 
c
        CALL PRR2('  SOURCE R AND Z POSITIONS      P1=', CXSCA,CYSCA)
        CALL PRR2('  LIE ON LINE JOINING P1 TO P2  P2=', CXSCB,CYSCB)
c
       elseif (cneutb.eq.7) then 
c
        CALL PRR2('  SOURCE R POSN FOR NEUT IN RANGE  ', CXSCA,CXSCB)
        CALL PRR2('  SOURCE Z POSN FOR NEUT IN RANGE  ', CYSCA,CYSCB)
c
       ENDIF
      ELSEIF (CNEUTA.EQ.1) THEN
       CALL PRR ('  INITIAL ION TEMPERATURE  (EV)    ', CTEM1)
       if (ciopte.eq.1) then
       CALL PRR ('  SOURCE R POSITION        (M)     ', CXSC)
       CALL PRR ('  SOURCE Z POSITION        (M)     ', CYSC)
c slmod begin
       elseif (ciopte.eq.9) then
c       CALL PRR ('  SOURCE R POSITION        (M)     ', 0.0)
c       CALL PRR ('  SOURCE Z POSITION        (M)     ', 1.0)
c slmod end
       endif
      ENDIF
      CALL PRI ('  INITIAL IONIZATION STATE         ', CIZSC)


      IF (CSTOP.EQ.1) CALL PRC ('  IONS NOT FOLLOWED AFTER REACHING MAIN
     > PLASMA')

c
C-----------------------------------------------------------------------
      CALL PRR  ('  PERPENDICULAR DIFFUSION      DPERP (M*M/S)', CDPERP)
      if (cdperpt.ne.cdperp)
     >CALL PRR  ('  PRIVATE PLASMA PERPENDICULAR DIFFUSION    ',CDPERPT)
c
      if (pinchopt.eq.0) then
         call prc  ('  PINCH OPTION 0:'// 
     >              '  NO PERPENDICULAR PINCH VELOCITY APPLIED')
      elseif (pinchopt.eq.1) then 
         call prc  ('  PINCH OPTION 1:')
         call prr  ('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
c
      elseif (pinchopt.eq.2) then
         call prc  ('  PINCH OPTION 2:')
         call prr('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc('  - PINCH IS APPLIED ONLY IN THE MAIN SOL')
c
c     added pinch velocity Krieger IPP/97
c
      elseif (pinchopt.eq.3) then
         call prc  ('  PINCH OPTION 3:')
         CALL PRR  ('  PERPENDICULAR PINCH AT SEP.'//
     >             '  PINCH (  M/S)', CVPINCH)
         call prc('  - PINCH IS APPLIED ONLY IN THE CORE')
         call prc('  - PINCH IS ADJUSTED BY POLOIDAL LENGTH FACTOR')
c
      elseif (pinchopt.eq.4) then
c
         call prc  ('  PINCH OPTION 4:')
         call prc('     - PDF BASED PERPENDICULAR VELOCITY CALCULATED')
         if (pinch_loc_opt.eq.0) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED EVERYWHERE'//
     >               'EXCEPT PFZ')
         elseif (pinch_loc_opt.eq.1) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED ONLY'//
     >               ' IN MAIN SOL')
         elseif (pinch_loc_opt.eq.2) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED EVERYWHERE')
         elseif (pinch_loc_opt.eq.3) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED ONLY IN'//
     >               ' MAIN SOL ABOVE X-POINT')
         elseif (pinch_loc_opt.eq.4) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED IN MAIN SOL'// 
     >               ' ABOVE X-POINT AND IN CORE REGION')
         endif
         call prc ('     - DPERP TRANSPORT TURNED OFF IN'//
     >                     ' PINCH REGION')
c
      elseif (pinchopt.eq.5) then
c
         call prc  ('  PINCH OPTION 5:')
         call prc('     - PDF BASED PERPENDICULAR VELOCITY CALCULATED')
         if (pinch_loc_opt.eq.0) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED EVERYWHERE'//
     >               'EXCEPT PFZ')
         elseif (pinch_loc_opt.eq.1) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED ONLY'//
     >               ' IN MAIN SOL')
         elseif (pinch_loc_opt.eq.2) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED EVERYWHERE')
         elseif (pinch_loc_opt.eq.3) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED ONLY IN'//
     >               ' MAIN SOL ABOVE X-POINT')
         elseif (pinch_loc_opt.eq.4) then 
           CALL PRC ('     - PDF BASED PINCH APPLIED IN MAIN SOL'// 
     >               ' ABOVE X-POINT AND IN CORE REGION')
         endif
         call prc ('     - DPERP DIFFUSIVE TRANSPORT ALSO'//
     >               ' ACTIVE IN PINCH REGION')
c
      elseif (pinchopt.eq.6) then 
         call prc  ('  PINCH OPTION 6:')
         call prr  ('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc  ('  - PINCH HAS OPPOSITE SIGN FOR R < Rxpt'//
     >              ' THAN FOR R > Rxpt')
         call prc  ('  - THIS APPROXIMATES AN OVERALL INWARD OR'//
     >              ' OUTWARD PINCH')
      elseif (pinchopt.eq.7) then 
         call prc  ('  PINCH OPTION 7:')
         call prr  ('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc  ('  - PINCH HAS OPPOSITE SIGN FOR R < Rxpt'//
     >              ' THAN FOR R > Rxpt')
         call prc  ('  - THIS APPROXIMATES AN OVERALL INWARD OR'//
     >              ' OUTWARD PINCH')
         call prc  ('  - PINCH IS NOT PRESENT INSIDE CONFINED PLASMA')
      elseif (pinchopt.eq.8) then 
         call prc  ('  PINCH OPTION 8:')
         call prr  ('  RADIAL PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc  ('  - PINCH IS PROJECTED ONTO S AND CROSS AXES')
         call prc  ('    BASED ON THE CELL GEOMETRY')
      elseif (pinchopt.eq.9) then 
         call prc  ('  PINCH OPTION 9:')
         call prr  ('  RADIAL PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc  ('  - PINCH IS PROJECTED ONTO S AND CROSS AXES')
         call prc  ('    BASED ON THE CELL GEOMETRY')
         call prc  ('  - CORE PLASMA EXCLUDED')
      elseif (pinchopt.eq.10) then 
         call prc('  PINCH OPTION 10:')
         call prr('  RADIAL PINCH VELOCITY APPLIED (M/S)',
     >               cvpinch)
         call prc('  - PINCH IS PROJECTED ONTO S AND CROSS AXES')
         call prc('    BASED ON THE CELL GEOMETRY')
         call prc('  - PINCH IS APPLIED ONLY IN THE MAIN SOL')
         call prc('  - PINCH IS APPLIED ONLY IN CELLS ABOVE THE XPOINT')
      elseif (pinchopt.eq.11) then
         call prc  ('  PINCH OPTION 11:')
         call prr('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc('  - PINCH IS APPLIED ONLY IN THE MAIN SOL')
         call prc('  - PINCH IS APPLIED ONLY IN CELLS ABOVE THE XPOINT')
      elseif (pinchopt.eq.12) then 
         call prc('  PINCH OPTION 12:')
         call prr('  RADIAL PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc('  - PINCH IS PROJECTED ONTO S AND CROSS AXES')
         call prc('    BASED ON THE CELL GEOMETRY')
         call prc('  - PINCH IS APPLIED ONLY IN THE INNER SOL : R < R0')
         call prc('  - PINCH IS APPLIED ONLY IN CELLS ABOVE THE XPOINT')
      elseif (pinchopt.eq.13) then
         call prc  ('  PINCH OPTION 13:')
         call prr('  PERPENDICULAR PINCH VELOCITY APPLIED (M/S)',
     >                 cvpinch)
         call prc('  - PINCH IS APPLIED ONLY IN THE INNER SOL : R < R0')
         call prc('  - PINCH IS APPLIED ONLY IN CELLS ABOVE THE XPOINT')
      endif
c

      RETURN
      END
c
c
c
      SUBROUTINE PR_options (NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
c
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
c
C     THIS ROUTINE DIRECTS PRINTING OF ALL THE OPTION FLAGS
c
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
c      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
      include    'printopt'
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
c      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
c      external lenstr
c
c
C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('----------------  OPTIONS ------------------- ')
      call prb 
      CALL PRChtml ('---  OPTIONS  ---',
     >               'pr_options','0','B')
      CALL PRB
c
C-----------------------------------------------------------------------
c
c     Print MISCELLANEOUS Options
c


       call pr_misc_options

!ammod begin
c
c     Print HYDROCARBON Options
c
       

       call pr_hydrocarbon_options
c
!ammod end

c
C-----------------------------------------------------------------------
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_misc_options
c
      IMPLICIT none
c
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS ALL THE MISCELLANEOUS OPTION FLAGS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
c      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
       include    'printopt' 
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
c      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
c      external lenstr
c
C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('------------ MISC  OPTIONS ------------------- ')
      CALL PRChtml ('--- MISC OPTIONS ---',
     >              'pr_misc_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------

      if (cprint.eq.0) then
        CALL PRC ('  PRINT OPTION 0     : STANDARD PRINT OPTION')
      elseif (cprint.eq.1) then
        CALL PRC ('  PRINT OPTION 1     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - Dperp EXTRACTOR')
        CALL PRC ('                       - SCANNING PROBE')
        CALL PRC ('                       - TRAP CONTENT DATA')
        CALL PRC ('                       IF THESE ARE SWITCHED ON')
      elseif (cprint.eq.2) then
        CALL PRC ('  PRINT OPTION 2     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - CORE LEAKAGE AND SOURCE ')
        CALL PRC ('                         DESCRIPTION INFORMATION')
      elseif (cprint.eq.3) then
        CALL PRC ('  PRINT OPTION 3     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - DEBUG - GEOMETRICAL DATA ABO
     >UT')
        CALL PRC ('                         GRID, TARGETS AND WALLS')
      elseif (cprint.eq.4) then
        CALL PRC ('  PRINT OPTION 4     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - INFORMATION ON EDGE2D TARGET
     >')
        CALL PRC ('                         CONDITIONS')
      elseif (cprint.eq.5) then
        CALL PRC ('  PRINT OPTION 5     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - ADDITIONAL INFORMATION ABOUT
     >')
        CALL PRC ('                         BACKGROUND CONDITIONS ')
        CALL PRC ('                       - SOME CHARACTERISTIC TIMES DA
     >TA')
        CALL PRC ('                       - RETENTION PREDICTOR VALUES')
      elseif (cprint.eq.6) then
        CALL PRC ('  PRINT OPTION 6     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - WRITES THE GRID IN A SONNET'
     >)
        CALL PRC ('                         STYLE FORMAT')
      elseif (cprint.eq.7) then
        CALL PRC ('  PRINT OPTION 7     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - EXTRA PIN RELATED DATA FOR')
        CALL PRC ('                         DEBUGGING')
      elseif (cprint.eq.8) then
        CALL PRC ('  PRINT OPTION 8     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - NOTHING EXTRA AT PRESENT')
      elseif (cprint.eq.9) then
        CALL PRC ('  PRINT OPTION 9     : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - EVERYTHING. COMBINES ALL PRI
     >NT')
        CALL PRC ('                         OPTIONS FROM 1 TO 8')
      elseif (cprint.eq.10) then
        CALL PRC ('  PRINT OPTION 10    : STANDARD PRINT OPTION + ')
        CALL PRC ('                       - DIVIMP BACKGROUND PLASMA FIL
     >E')
      endif
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_geom_options
      IMPLICIT none
c
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS ALL THE GEOMETRY OPTION FLAGS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
      include    'printopt'    
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new

      INTEGER j,k
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
      INTEGER  I,ir
c      logical  prsol21,prsol22
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
      integer lenval,lenstr
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('--------- GEOMETRY GRID AND WALL OPTIONS -------')
      CALL PRChtml ('--- GEOMETRY GRID AND WALL OPTIONS ---',
     >              'pr_geom_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      IF (CGRIDOPT.EQ.0) then
        CALL PRC ('  GRID OPTION 0:  STANDARD JET GRID FILES')
      ELSEIF (CGRIDOPT.EQ.1) then
        CALL PRC ('  GRID OPTION 1:  STANDARD ASDEX GRID FILES')
      ELSEIF (CGRIDOPT.EQ.2) then
        CALL PRC ('  GRID OPTION 2:  STANDARD ITER GRID FILES')
      ELSEIF (CGRIDOPT.EQ.3) then
        CALL PRC ('  GRID OPTION 3:  STANDARD SONNET (UPGRADE, CMOD, DII
     >ID) GRID FILES')

        if (sonnet_grid_sub_type.eq.1) then 
          call prc(sp//'SONNET GRID SUB-TYPE: 1')
          call prc(sp//'- USE FRC VERSION 1 SONNET GRID '
     >               //' CUSTOMIZATIONS')
        endif 
c
        if (sonnet_grid_sub_type.eq.2) then 
          call prc(sp//'SONNET GRID SUB-TYPE: 2')
          call prc(sp//'- SONNET GRID WITHOUT BOUNDARY CELLS - '
     >               //'CELLS ARE ADDED FOR INTERNAL PROCESSING')
        endif 


      ELSEIF (CGRIDOPT.EQ.LINEAR_GRID) then
        CALL PRC ('  GRID OPTION 6:  LINEAR DEVICE GRID BUILT')
      ELSEIF (CGRIDOPT.EQ.GEN_GRID) then
        CALL PRC ('  GRID OPTION 6:  GENERALIZED GRID FILE :'//
     >            ' SONNET BASE')
      ENDIF
c
c     TMACHINE_OPT is an optional input value used to set the 
c     machine type in the JET post-processor TRAN file. The 
c     default value is set to JET. 
c

      if (tmachine_opt.eq.0) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  0:  JET')
      elseif (tmachine_opt.eq.1) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  1:  DIIID')
      elseif (tmachine_opt.eq.2) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  2:  ALCATOR (CMO
     >D)')
      elseif (tmachine_opt.eq.3) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  3:  AUG (ASDEX U
     >PGRADE)')
      elseif (tmachine_opt.eq.4) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  4:  ITER')
      elseif (tmachine_opt.eq.5) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  5:  IGNITOR')
      elseif (tmachine_opt.eq.6) then  
        CALL PRC ('  MACHINE TYPE (*G34) FOR TRAN FILE  6:  FIRE')
      endif
c
c     DIVSHOTID is an optional input that will specify the shot ID 
c     explicitly for the TRAN file. If not specified then the SHOT ID
c     (a shot number) is extracted from the case name by the cataloguing
c     routines. 
c
      call prb 
      if (divshotid.eq.' ') then 
        CALL PRC ('  SHOT ID (*G35) NOT SPECIFIED: SHOT'//
     >                 ' ID WILL BE EXTRACTED FROM CASE NAME')
      else
        lenval = lenstr(divshotid)  
        CALL PRC ('  SHOT ID (*G35) SPECIFIED: SHOT'//
     >                 ' ID WILL BE '//divshotid(1:lenval))
      endif 
c
c     PARALLEL ION REFLECTION
c 
c     This option is used for broken grids in order to allow for
c     ion reflection at the breaks. This allows for simulation of 
c     1/2 of a double null configuration.
c
      call prb
c
      if (s_reflect_opt.eq.0) then      
        call prc('  PARALLEL ION REFLECTION OPT 0: (*G36) THIS'//
     >              ' OPTION IS OFF')
        call prc('                                 IONS ARE NOT'//
     >           ' REFLECTED AT GRID PARALLEL GEOMETRY ERRORS')
      elseif (s_reflect_opt.eq.1) then      
        call prc('  PARALLEL ION REFLECTION OPT 1: (*G36) THIS'//
     >              ' OPTION IS ON')
        call prc('                                 IONS ARE'//
     >           ' REFLECTED AT GRID PARALLEL GEOMETRY ERRORS')
        call prc('       TABLE OF CALCULATED S-REFLECTION VALUES:')
        call prc('       (VALUES OF 0.0 MEAN NO REFLECTION)')
        call prc('    RING     REFLECT S<     REFLECT S>  ')
        do ir = 1,nrs
           write(coment,'(4x,i4,2(3x,f10.3))') ir,
     >            s_reflect(ir,1),s_reflect(ir,2)
           call prc(coment)
        end do
c
      endif 
c
c 
      call prb
c
      if (northopt.eq.0) then
        CALL PRC ('  NON-OTHOGONAL OPT 0: GRIDS ARE TREATED AS ORTHOGONA
     >L EVERYWHERE')
      elseif (northopt.eq.1) then
        CALL PRC ('  NON-OTHOGONAL OPT 1: JET GRID NON-ORTHOGONAL TREATM
     >ENT')
      elseif (northopt.eq.2) then
        CALL PRC ('  NON-OTHOGONAL OPT 2: TARGET NON-ORTHOGONALITY IS US
     >ED')
        call prc ('                       TO CALCULATE FLUXES ONLY. PART
     >ICLE')
        call prc ('                       TRANSPORT IS ORTHOGONAL')
c
      elseif (northopt.eq.3) then
        CALL PRC ('  NON-OTHOGONAL OPT 3: NON-ORTHOGONAL TREATMENT')
        call prc ('                       BASED ON CALCULATING AN INDEPE
     >NDENT')
        call prc ('                       ORTHOGONAL CO-ORDINATE (THETAG
     >)')
        call prc ('                       FROM THE GRID POLYGONS')
        call prc ('                       USABLE FOR JET OR SONNET GRIDS
     >')
        call prc ('                       AND BASED ON THE METHODS USED
     >IN')
        call prc ('                       JET GRIDS')
        call prc ('                       SHOULD BE USED WITH:')
        call prc ('                       PARALLEL DISTANCE OPTION 1')
      endif
      call prb
C-----------------------------------------------------------------------
      if (rzopt.eq.0) then
         call prc('  R,Z OPTION        0: CELL CENTER R,Z VALUES ARE USE
     >D')
         CALL PRC('                       TO ESTIMATE PARTICLE POSITION'
     >)
      ELSEif (rzopt.eq.1) then
         call prc('  R,Z OPTION        1: SIMPLE GETRZ ROUTINE IS'//
     >            ' USED TO CALCULATE')
         CALL PRC('                       ACTUAL PARTICLE POSITION.')
         call prc('                       - CROSS VALUE IS NOT USED.')
         call prc('                       - END CELLS NOT ACCURATE')
      ELSEif (rzopt.eq.2) then
         call prc('  R,Z OPTION        2: MODIFIED GETRZ ROUTINE IS'//
     >            ' USED TO CALCULATE')
         CALL PRC('                       ACTUAL PARTICLE POSITION')
         call prc('                       - CROSS IS USED.')
         call prc('                       - CALCULATED R,Z WILL'//
     >            ' BE IN CELL')
      ELSEif (rzopt.eq.3) then
         call prc('  R,Z OPTION        2: MODIFIED GETRZ ROUTINE IS'//
     >            ' USED TO CALCULATE')
         CALL PRC('                       ACTUAL PARTICLE POSITION')
         call prc('                       - S AND CROSS ARE'//
     >            ' ORTHOGONAL AT THE CELL CENTER')
         call prc('                       - RESULT MAY NOT BE'//
     >            ' IN CELL OR ON GRID')
      ENDIF
C-----------------------------------------------------------------------
      IF     (cvolopt.EQ.0) THEN
       CALL PRC ('  CELL AREAS OPTION 0: APPROXIMATED FROM CELL CENTRES'
     >)
       call prc ('                       ****** IMPORTANT NOTE ******')
       call prc ('                       USE PARALLEL DISTANCE OPTION 0'
     >)
       call prc ('                       AND CROSS-FIELD DISTANCE OPTION
     > 0')
       CALL PRC ('                       IN ORDER TO GET TRANSPORT AND P
     >ARTICLE')
       CALL PRC ('                       DENSITIES COMPATIBLE WITH THIS
     >METHOD')
       CALL PRC ('                       OF CALCULATING CELL AREAS.')
      ELSEIF (cvolopt.EQ.1) THEN
       CALL PRC ('  CELL AREAS OPTION 1: CALCULATED FROM POLYGON VERTICE
     >S')
       call prc ('                       ****** IMPORTANT NOTE ******')
       call prc ('                       USE PARALLEL DISTANCE OPTION 1'
     >)
       call prc ('                       AND CROSS-FIELD DISTANCE OPTION
     > 1+')
       CALL PRC ('                       IN ORDER TO GET TRANSPORT AND P
     >ARTICLE')
       CALL PRC ('                       DENSITIES COMPATIBLE WITH THIS
     >METHOD')
       CALL PRC ('                       OF CALCULATING CELL AREAS.')
      endif
C-----------------------------------------------------------------------
      if (pdopt.eq.0) then
        call prc('  PARALLEL DIST OPT 0: STANDARD (ORIGINAL) OPTION')
        call prc('                       S DISTANCES ALONG THE FIELD LIN
     >ES')
        CALL PRC('                       ARE CALCULATED BY TAKING THE DI
     >STANCES')
        CALL PRC('                       BETWEEN CELL CENTRES.')
      ELSEIF (PDOPT.EQ.1) THEN
        call prc('  PARALLEL DIST OPT 1: POLYGON GRID OPTION')
        call prc('                       S DISTANCES ALONG THE FIELD LIN
     >ES')
        CALL PRC('                       ARE CALCULATED BY TAKING THE DI
     >STANCES')
        CALL PRC('                       BETWEEN THE MID-POINTS OF OPPOS
     >ING.')
        CALL PRC('                       CELL POLYGON SIDES AND PASSING
     >THROUGH')
        CALL PRC('                       THE CENTRE OF THE CELL.')
      endif
C-----------------------------------------------------------------------
      if (cfdopt.eq.0) then
        call prc('  CROSS-FIELD DIST  0: STANDARD (ORIGINAL) OPTION')
        call prc('                       CROSS-FIELD DIFFUSIVE DISTANCES
     > ARE')
        CALL PRC('                       MEASURED TO THE MID-POINT BETWE
     >EN')
        CALL PRC('                       ADJACENT CELL CENTRES')
      ELSEIF (cfdOPT.EQ.1) THEN
        call prc('  CROSS-FIELD DIST  1: POLYGON GRID OPTION')
        call prc('                       THE CROSS-FIELD DISTANCES REQUI
     >RED')
        CALL PRC('                       FOR DIFFUSION BETWEEN RINGS ARE
     > CALCULATED')
        CALL PRC('                       BY USING THE ACTUAL POLYGON CEL
     >L SIZES')
      endif
C-----------------------------------------------------------------------
      IF     (xygrid.EQ.0) THEN
       CALL PRC ('  XY GRID OPTION   0 : OFF - XY GRID IS NOT GENERATED'
     >)
      ELSEIF (xygrid.EQ.1) THEN
       CALL PRC ('  XY GRID OPTION   1 : ON - XY GRID MAPPING CREATED AN
     > STORED')
      endif
C-----------------------------------------------------------------------
      IF (CTARGOPT.EQ.0) THEN
       CALL PRC ('  TARGET OPTION    0 : TARGET IS LOCATED AT SECOND GRI
     >D POINTS ON')
       CALL PRC ('                       THE SOL AND TRAP RINGS.')
       call prc ('                       VIRTUAL POINTS ARE DISCARDED')
      ELSEIF (CTARGOPT.EQ.1) THEN
       CALL PRC ('  TARGET OPTION    1 : TARGET IS LOCATED MID-WAY BETWE
     >EN THE')
       CALL PRC ('                       VIRTUAL POINT AND FIRST REAL PO
     >INT ON')
       call prc ('                       THE SOL AND TRAP RINGS.')
       call prc ('                       VIRTUAL POINTS ARE THEN DISCARD
     >ED')
      ELSEIF (CTARGOPT.EQ.2) THEN
       CALL PRC ('  TARGET OPTION    2 : TARGET IS SPECIFIED BY A SET OF
     > POINTS')
       CALL PRC ('                       ENTERED IN THE DATA FILE. ONE P
     >OINT')
       CALL PRC ('                       FOR EACH END OF EACH RING.')
       call prc ('                       VIRTUAL POINTS ARE DISCARDED')
      ELSEIF (CTARGOPT.EQ.3) THEN
       CALL PRC ('  TARGET OPTION    3 : TARGET IS SPECIFIED BY A SET OF
     > POINTS')
       CALL PRC ('                       THAT ARE HARD-CODED. ONE POINT'
     >)
       CALL PRC ('                       FOR EACH END OF EACH RING.')
       call prc ('                       THE SET OF POINTS IS SELECTED B
     >Y THE')
       call prc ('                       GEOMETRY OPTION')
       call prc ('                       VIRTUAL POINTS ARE DISCARDED')
      ELSEIF (CTARGOPT.EQ.4) THEN
       CALL PRC ('  TARGET OPTION    4 : TARGET IS LOCATED AT THE FIRST
     >GRID POINTS ON')
       CALL PRC ('                       THE SOL AND TRAP RINGS.')
       call prc ('                       VIRTUAL POINTS ARE NOT DISCARDE
     >D')
      ELSEIF (CTARGOPT.EQ.5) THEN
       CALL PRC ('  TARGET OPTION    5 : TARGET IS SPECIFIED BY A SET OF
     > POINTS')
       CALL PRC ('                       ENTERED IN THE DATA FILE. ONE P
     >OINT')
       CALL PRC ('                       FOR EACH END OF EACH RING.')
       call prc ('                       VIRTUAL POINTS ARE NOT DISCARDE
     >D')
      ELSEIF (CTARGOPT.EQ.6) THEN
       CALL PRC ('  TARGET OPTION    6 : TARGET IS SPECIFIED BY THE OUTE
     >R')
       CALL PRC ('                       POLYGON BOUNDARIES OF THE LAST
     >SET')
       CALL PRC ('                       OF REAL PLASMA POLYGONS ON EACH
     > RING')
       call prc ('                       VIRTUAL POINTS ARE DISCARDED')
      ENDIF
C-----------------------------------------------------------------------
 999  IF     (CIONR.EQ.0) THEN
       CALL PRC ('  ION WALL OPTION  0 : ION WALLS MIDWAY BETWEEN LAST 2
     > RINGS')
      ELSEIF (CIONR.EQ.1) THEN
       CALL PRC ('  ION WALL OPTION  1 : ION WALLS AT LAST RING')
      ELSEIF (CIONR.EQ.2) THEN
       CALL PRC ('  ION WALL OPTION  2 : ION WALLS AT THE POLYGON EDGE O
     >F THE')
       call prc ('                       OUTERMOST RINGS OF POLYGONS. IO
     >N')
       call prc ('                       TRANSPORT TO WALLS IS THE SAME'
     >)
       call prc ('                       OPTION 0.')
      ENDIF
C-----------------------------------------------------------------------
1001  IF     (CNEUR.EQ.0) THEN
       CALL PRC ('  NEUT WALL OPTION 0 : NEUTRAL WALLS MIDWAY BETWEEN LA
     >ST 2 RINGS')
      ELSEIF (CNEUR.EQ.1) THEN
       CALL PRC ('  NEUT WALL OPTION 1 : NEUTRAL WALLS AT LAST RING')
      ELSEIF (CNEUR.EQ.2) THEN
       CALL PRC ('  NEUT WALL OPTION 2 : MAIN NEUTRAL WALLS SPECIFIED BY
     >')
       call prc ('                       THE FOLLOWING VALUES IN THE INP
     >UT FILE')
c slmod begin - new
       IF (new) THEN
         DO i = 1, (nwall + 1) / 3
           j = i + (nwall + 1) / 3
           k = j + (nwall + 1) / 3
           WRITE(coment,1009) i,wallco(i,1),wallco(i,2)
           IF (j.LE.nwall)
     .       WRITE(coment(LEN_TRIM(coment)+1:LEN(coment)),1009)
     .         j,wallco(j,1),wallco(j,2)
           IF (k.LE.nwall)
     .       WRITE(coment(LEN_TRIM(coment)+1:LEN(coment)),1009)
     .         k,wallco(k,1),wallco(k,2)
           CALL PRC(coment)
         ENDDO
         CALL PRB
1009     FORMAT(4X,I4,2F8.4)
       ELSE
         DO 1010 I = 1,NWALL
           CALL PRR2('                      ',WALLCO(I,1),WALLCO(I,2))
1010     CONTINUE
       ENDIF
c
c       DO 1010 I = 1,NWALL
c       CALL PRR2('                      ',WALLCO(I,1),WALLCO(I,2))
c1010   CONTINUE
c slmod end
      ELSEIF (CNEUR.EQ.3) THEN
       CALL PRC ('  NEUT WALL OPTION 3 : NEUTRAL WALLS SPECIFIED BY')
       call prc ('                       HARD-CODED DATA')
c
c       DO 1020 I = 1,NWALL
c       CALL PRR3('                      ',WALLCO(I,1),WALLCO(I,2),
c     >           real(wallpol(i)))
c1020   CONTINUE
c
      elseif (cneur.eq.4) then
       CALL PRC ('  NEUT WALL OPTION 4 : NEUTRAL WALLS SPECIFIED BY')
       call prc ('                       THE VESSEL COORIDNATES IN THE')
       call prc ('                       GRID FILE')
c
c       DO  I = 1,NWALL
c       CALL PRR2('                      ',WALLCO(I,1),WALLCO(I,2))
c       end do
c
      elseif (cneur.eq.5) then
       CALL PRC ('  NEUT WALL OPTION 5 : NEUTRAL WALLS ARE READ FROM')
       call prc ('                       THE PIN/NIMBUS TRANSFER FILE')
c
c       DO  I = 1,NWALL
c       CALL PRR2('                      ',WALLCO(I,1),WALLCO(I,2))
c       end do
c
      ELSEIF (CNEUR.EQ.7) THEN
       CALL PRC ('  NEUT WALL OPTION 7 : NEUTRAL WALLS AT OUTER POLYGON
     >BOUNDARY')
       call prc ('                       OF LAST REAL RING.')
      ENDIF
C-----------------------------------------------------------------------
      IF     (CTRAP.EQ.0) THEN
       CALL PRC ('  TRAP WALL OPTION 0 : NEUTRAL TRAP WALLS MIDWAY BETWE
     >EN LAST')
       CALL PRC ('                       2 RINGS')
       CALL PRC ('                       IONS ARE TREATED AS SPECIFIED I
     >N THE')
       CALL PRC ('                       SELECTED WALL OPTION')
      ELSEIF (CTRAP.EQ.1) THEN
       CALL PRC ('  TRAP WALL OPTION 1 : NEUTRAL TRAP WALLS AT LAST RING
     >')
       CALL PRC ('                       IONS ARE TREATED AS SPECIFIED I
     >N THE')
       CALL PRC ('                       SELECTED WALL OPTION')
      ELSEIF (CTRAP.EQ.2) THEN
       CALL PRC ('  TRAP WALL OPTION 2 : NEUTRAL TRAP WALL CREATED BY JO
     >INING THE')
       CALL PRC ('                       ENDS OF THE INNER AND OUTER TAR
     >GETS')
       CALL PRC ('                       IONS ARE TREATED AS SPECIFIED I
     >N THE')
       CALL PRC ('                       SELECTED WALL OPTION')
      ELSEIF (ctrap.eq.3) THEN
       CALL PRC ('  TRAP WALL OPTION 3 : NEUTRAL TRAP WALL SPECIFIED BY
     >SET OF   ')
       CALL PRC ('                       ADDITIONAL COORDINATES. IONS AR
     >E   ')
       CALL PRC ('                       TREATED AS SPECIFIED IN THE SEL
     >ECTED')
       CALL PRC ('                       WALL OPTION. TRAP WALL SPECIFIE
     >D BY')
c slmod begin - new
       IF (new) THEN
         DO i = 1, (nwall2 + 1) / 3
           j = i + (nwall2 + 1) / 3
           k = j + (nwall2 + 1) / 3
           WRITE(coment,1009) i,wallco2(i,1),wallco2(i,2)
           IF (j.LE.nwall)
     .       WRITE(coment(LEN_TRIM(coment)+1:LEN(coment)),1009)
     .         j,wallco2(j,1),wallco2(j,2)
           IF (k.LE.nwall)
     .       WRITE(coment(LEN_TRIM(coment)+1:LEN(coment)),1009)
     .         k,wallco2(k,1),wallco2(k,2)
           CALL PRC(coment)
         ENDDO
         CALL PRB
       ELSE
         DO I = 1,NWALL2
           CALL PRR2('                      ',WALLCO2(I,1),WALLCO2(I,2))
         end do
       ENDIF
c
c       DO I = 1,NWALL2
c       CALL PRR2('                      ',WALLCO2(I,1),WALLCO2(I,2))
c       end do
c slmod end
      ELSEIF (ctrap.eq.4) THEN
       CALL PRC ('  TRAP WALL OPTION 4 : NEUTRAL TRAP WALL SPECIFIED BY
     >SET OF   ')
       CALL PRC ('                       ADDITIONAL COORDINATES TAKEN F
     >ROM THE')
       call prc ('                       GRID GEOMETRY FILE. IONS ARE')
       CALL PRC ('                       TREATED AS SPECIFIED IN THE SEL
     >ECTED')
       CALL PRC ('                       WALL OPTION.')
c
c       DO I = 1,NWALL2
c       CALL PRR2('                      ',WALLCO2(I,1),WALLCO2(I,2))
c       end do
c
      elseif (ctrap.eq.5) then
       CALL PRC ('  TRAP WALL OPTION 5 : NEUTRAL WALLS ARE READ FROM')
       call prc ('                       THE PIN/NIMBUS TRANSFER FILE')
c
c       DO  I = 1,NWALL2
c       CALL PRR2('                      ',WALLCO2(I,1),WALLCO2(I,2))
c       end do
c
      ELSEIF (CTRAP.EQ.7) THEN
       CALL PRC ('  TRAP WALL OPTION 7 : NEUTRAL WALLS AT OUTER POLYGON
     >BOUNDARY')
       call prc ('                       OF LAST REAL PRIVATE PLASMA RIN
     >G.')
      ELSEIF (CTRAP.EQ.8) THEN
       CALL PRC ('  TRAP WALL OPTION 8 : THERE IS NO WALL DEFINITION '
     >//'USED FROM THE PFZ')
       call prc ('                       THIS OPTION IS ONLY USED WITH'
     >//' LIMITER GRIDS')
       call prc ('                       WHICH DO NOT HAVE A PFZ.')
      ENDIF

c
C-----------------------------------------------------------------------
c      REDFEOPT  - Vessel Wall redefinition option
C-----------------------------------------------------------------------
c
      IF     (redefopt.EQ.0) THEN
       CALL PRC ('  VESSEL REDEF OPT 0 : VESSEL REDEFINITION IS OFF')
      ELSEIF (redefopt.EQ.1) THEN
       CALL PRC ('  VESSEL REDEF OPT 1 : VESSEL REDEFINITION IS ON')
       CALL PRC ('                       VESSEL WALL IS MODIFIED TO INCL
     >UDE BAFFLES')
      endif
c
C-----------------------------------------------------------------------
c      List the elements in the entire neutral wall - with index numbers
C-----------------------------------------------------------------------
c
c
c     Print out wall temperatures for chemical sputtering cases.
c
      if (cneutd.eq.2.or.cneutd.eq.5.or.cneutd.eq.6.or.
     >    cneutd2.eq.2.or.cneutd2.eq.5.or.cneutd2.eq.6)
     >   call prtemp
c
c      Print out wall elements 
c
       call prb
c slmod begin
       IF (new) THEN
         call prc ('  NEUTRAL WALL ELEMENT LISTING:')
         WRITE(coment,'(2X,2A8,A7,A9)') 'R (m)','Z (m)','Wall #',
     .                                  'Target #'
         CALL PRC(coment)
         DO i = 1, (wallpts+1)/2
          j = i + (wallpts + 1) / 2
          WRITE(coment,1011) wallpt(i,1),wallpt(i,2),i,INT(wallpt(i,18))
          IF (j.LE.wallpts)
     .      WRITE(coment(LEN_TRIM(coment)+5:LEN(coment)),1011)
     .        wallpt(j,1),wallpt(j,2),j,INT(wallpt(j,18))
          CALL PRC(coment)
         ENDDO
1011     FORMAT(2X,2F8.4,I7,I9)
       ELSE
         call prc ('                     NEUTRAL WALL ELEMENT LISTING:')
         write (coment,'(22x,5x,''R (m)'',13x,''Z (m)'',6x,''Wall#'',2x,
     >          ''Targ#'')')
         call prc (coment)
         do i = 1,wallpts
            write (coment,'(16x,2f18.10,2(1x,i5))') wallpt(i,1),
     >         wallpt(i,2),i,int(wallpt(i,18))
          call prc(coment)
         end do
       ENDIF
c
c       call prc ('                       NEUTRAL WALL ELEMENT LISTING:')
c
c       write (coment,'(24x,5x,''R (m)'',13x,''Z (m)'',6x,''Wall#'',2x,
c     >        ''Targ#'')')
c       call prc (coment)
c
c       do i = 1,wallpts
c          write (coment,'(18x,2f18.10,2(1x,i5))') wallpt(i,1),
c     >             wallpt(i,2),i,int(wallpt(i,18))
c
c          call prc(coment)
c       end do
c slmod end
       call prb

C-----------------------------------------------------------------------
      IF (CGEOOPT.EQ.-1) THEN
       CALL PRC ('  GEOMETRY OPTION -1 : SHOT WALL AND TARGET DATA')
       CALL PRC ('                       IS NOT PRE-LOADED')
      ELSEIF (CGEOOPT.EQ.0) THEN
       CALL PRC ('  GEOMETRY OPTION  0 : SHOT 24719 WALL AND TARGET')
       CALL PRC ('                       DATA AVAILABLE')
      ELSEIF (CGEOOPT.EQ.1) THEN
       CALL PRC ('  GEOMETRY OPTION  1 : SHOT 26308 WALL AND TARGET') 
       CALL PRC ('                       DATA AVAILABLE')
      ENDIF
 
      RETURN
      END
c
c
c
      SUBROUTINE PR_procdata_options(NIZS)
      IMPLICIT none
c
      INTEGER NIZS
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS ALL ATOMIC PROCESS AND ATOMIC DATA OPTIONS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
      include    'cadas'
      include    'cgeom'
      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
      include    'adpak_com'
      include    'reiser_com'
c      include    'promptdep'
c      include    'slcom'  

      include    'dperpz'
      include    'printopt'
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod endc
c      CHARACTER  COMENT*80,prtype*4
c      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
      INTEGER  len,lenstr
c      logical  prsol21,prsol22
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('----- ATOMIC PROCESS AND DATA OPTIONS ------- ')
      CALL PRChtml ('--- ATOMIC PROCESS AND DATA OPTIONS ---',
     >              'pr_proc_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------
c
c  SIMULATION - ATOMIC DATA and PROCESSES
c
      IF     (cdatopt.EQ.0) THEN
       CALL PRC ('  DATA SOURCE OPT  0 : NOCORONA DATABASE USED')
      ELSEIF (cdatopt.EQ.1) THEN
       CALL PRC ('  DATA SOURCE OPT  1 : ADAS DATABASE USED')
       CALL PRC ('                       USERID FOR HYDROGEN:')
       LEN = lenstr(useridh)
       CALL PRC ('                       '//USERIDH(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR HYDROGEN : '',I2)') IYEARH
       CALL PRC (COMENT)
       CALL PRC ('                       USERID FOR IMPURITY:')
       LEN = lenstr(useridz)
       CALL PRC ('                       '//USERIDZ(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR IMPURITY : '',I2)') IYEARZ
       CALL PRC (COMENT)
       call prr ('               ADAS IONIZATION    RATE MULTIPLIER = ',
     >                              adas_iz_rate_mult)
       call prr ('               ADAS RECOMBINATION RATE MULTIPLIER = ',
     >                              adas_rec_rate_mult)
      ELSEIF (cdatopt.EQ.2) THEN
       CALL PRC ('  DATA SOURCE OPT  2 : B2-FRATES DATABASE USED')
       call prc ('                       APPLIES TO IONIZATION AND RECOM
     >BINATION')
       call prc ('                       OF CARBON ONLY (AT PRESENT)')
       call prc ('                       SPECIFIC DATASET USED:')
       len = lenstr(mcfile)
       call prc ('                       '//mcfile(1:len))
c
       call prb
       CALL PRC ('                       ADAS DATABASE IS USED')
       call prc ('                       FOR ALL OTHER CALCULATIONS')
       CALL PRC ('                       USERID FOR HYDROGEN:')
       LEN = lenstr(useridh)
       CALL PRC ('                       '//USERIDH(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR HYDROGEN : '',I2)') IYEARH
       CALL PRC (COMENT)
       CALL PRC ('                       USERID FOR IMPURITY:')
       LEN = lenstr(useridz)
       CALL PRC ('                       '//USERIDZ(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR IMPURITY : '',I2)') IYEARZ
       CALL PRC (COMENT)
c
      ELSEIF (cdatopt.EQ.3) THEN
       CALL PRC ('  DATA SOURCE OPT  3 : INEL DATABASE USED')
       call prc ('                       APPLIES TO IONIZATION AND RECOM
     >BINATION')
       call prc ('                       OF CARBON ONLY (AT PRESENT)')
       call prc ('                       SPECIFIC DATASET USED:')
       len = lenstr(mcfile)
       call prc ('                       '//mcfile(1:len))
c
       call prb
       CALL PRC ('                       ADAS DATABASE IS USED')
       call prc ('                       FOR ALL OTHER CALCULATIONS')
       CALL PRC ('                       USERID FOR HYDROGEN:')
       LEN = lenstr(useridh)
       CALL PRC ('                       '//USERIDH(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR HYDROGEN : '',I2)') IYEARH
       CALL PRC (COMENT)
       CALL PRC ('                       USERID FOR IMPURITY:')
       LEN = lenstr(useridz)
       CALL PRC ('                       '//USERIDZ(1:LEN))
       WRITE (COMENT,'(23X,''YEAR FOR IMPURITY : '',I2)') IYEARZ
       CALL PRC (COMENT)
c
      endif
C-----------------------------------------------------------------------
      IF     (CIOPTA.EQ.0) THEN
       CALL PRC ('  IONISATION OPT   0 : RATES FROM S(Z,TE) DATA')
      ELSEIF (CIOPTA.EQ.1) THEN
       CALL PRC ('  IONISATION OPT   1 : RATES FROM S(Z,TE) DATA')
      ELSEIF (CIOPTA.EQ.2) THEN
       CALL PRC ('  IONISATION OPT   2 : RATES TAKEN AS MAX (S(Z,TE))')
      ELSEIF (CIOPTA.EQ.3) THEN
       CALL PRC ('  IONISATION OPT   3 : IONISATION PLUS E-I RECOMBINATI
     >ON')
      ELSEIF (CIOPTA.EQ.4) THEN
       CALL PRC ('  IONISATION OPT   4 : IMPURITY RECOMBINATION OFF')
      ELSEIF (CIOPTA.EQ.5) THEN
       CALL PRC ('  IONISATION OPT   5 : IONISATION PLUS E-I RECOMBINATI
     >ON')
      ELSEIF (CIOPTA.EQ.6) THEN
       CALL PRC ('  IONISATION OPT   6 : IMPURITY RECOMBINATION OFF')
      ENDIF
C
      IF (NIZS.LT.CION) THEN
       IF (CIOPTA.EQ.0.OR.CIOPTA.EQ.3.OR.CIOPTA.EQ.4) THEN
        WRITE (COMENT,'(23X,''IONS NOT FOLLOWED AFTER STATE'',I3)') NIZS
        CALL PRC (COMENT)
       ELSE
        WRITE (COMENT,'(23X,''IONISATION DISABLED AFTER STATE'',I3)')
     >   NIZS
        CALL PRC (COMENT)
       ENDIF
      ENDIF
      IF (CNEUTA.EQ.0)
     > CALL PRR ('                       NEUTRAL IONISATION RATE FACTOR'
     >, CIRF)
C-----------------------------------------------------------------------
      WRITE (COMENT,'(23X,A,I3,2X,A)') 'FOR RINGS >=', IRSPEC,
     >  '.  ELSEWHERE'
C-----------------------------------------------------------------------

      if (cioptr.eq.0) then
c
c       Reiser is OFF
c
        CALL PRC('  REISER COLL. OPT 0: REISER COULOMB COLLISION OPTION'
     >)
        CALL PRC('                      HAS BEEN TURNED OFF')
c
c       Print NORMAL collision and friction options
c
c
      IF     (CIOPTB.EQ.0) THEN
       CALL PRC ('  COLLISION OPTION 0 : TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ELSEIF (CIOPTB.EQ.1) THEN
       CALL PRC ('  COLLISION OPTION 1 : TAU PARA = INFINITY, NO DIFFUSI
     >ON')
      ELSEIF (CIOPTB.EQ.2) THEN
       CALL PRC ('  COLLISION OPTION 2 : TAU PARA = SQRT(MI.TI).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZEFF.ZI.ZI
     >.LAM)')
       CALL PRI ('                       WHERE ZEFF(SELF) =',CIZEFF)
      ELSEIF (CIOPTB.EQ.3) THEN
       CALL PRC ('  COLLISION OPTION 3 : TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = TAU
     > PARA')
       IF (IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = DEL
     >TA T')
       ENDIF
      ELSEIF (CIOPTB.EQ.4) THEN
       CALL PRC ('  COLLISION OPTION 4 : TAU PARA = MB.SQRT(TI/MI).TI.TI
     >/')
       CALL PRC ('                            (9.0E4.TB.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       CALL PRC ('                       WHEN TI > TB.MI/MB AND')
       CALL PRC (COMENT)
       CALL PRC ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = TAU
     > PARA')
       IF (IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = DEL
     >TA T')
       ENDIF
      ELSEIF (CIOPTB.EQ.5) THEN
       CALL PRC ('  COLLISION OPTION 5 : PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = SQRT(8kTI/PI.MI)')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                       Diffusion occurs if : ')
       call prc ('                       0 < (random) < dt/(TAU PARA)')
      ELSEIF (CIOPTB.EQ.6) THEN
       CALL PRC ('  COLLISION OPTION 6 : PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = SQRT(8kTI/PI.MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ELSEIF (CIOPTB.EQ.7) THEN
       CALL PRC ('  COLLISION OPTION 7 : TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = TAU
     > PARA')
       IF (IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = DEL
     >TA T')
       ENDIF
       call prc ('                       Diffusive steps in the directio
     >n')
       call prc ('                       opposite of the particles veloc
     >ity')
       call prc ('                       reverse the sign of that V. ')
      ELSEIF (CIOPTB.EQ.8) THEN
       CALL PRC ('  COLLISION OPTION 8 : TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                       S Diffusive steps are based on:
     >')
       call prc ('                         sqrt(2.0 * kTi/Mi)*TauPara')
      ELSEIF (CIOPTB.EQ.9) THEN
       CALL PRC ('  COLLISION OPTION 9 : TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                       S Diffusive steps are based on:
     >')
       call prc ('                         sqrt(2.0 * kTi/Mi)*TauPara')
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = TAU
     > PARA')
       IF (IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TIME BETWEEN S DIFF STEPS = DEL
     >TA T')
       ENDIF
      ELSEIF (CIOPTB.EQ.10) THEN
       CALL PRC ('  COLLISION OPTION 10: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = SQRT(2kTI/MI)')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                       Diffusion occurs if : ')
       call prc ('                       0 < (random) < dt/(TAU PARA)')
      ELSEIF (CIOPTB.EQ.11) THEN
       CALL PRC ('  COLLISION OPTION 11: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = SQRT(2kTI/MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ELSEIF (CIOPTB.EQ.12) THEN
       CALL PRC ('  COLLISION OPTION 12: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = Rg * SQRT(2kTI/MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ELSEIF (CIOPTB.EQ.13) THEN
       CALL PRC ('  COLLISION OPTION 13: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = Rg * SQRT(2kTI/MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                               / (1.0+MB/MI) ')
      ELSEIF (CIOPTB.EQ.14) THEN
       CALL PRC ('  COLLISION OPTION 14: PARALLEL VELOCITY DIFFUSION')
       call prc ('                       DELTAV = Rg * SQRT(2kTI/MI)')
       call prc ('                           * sqrt( dt / (tau para))')
       call prc ('                       At every time step       ')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       call prc ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
       call prc ('                               / (1.0+MB/MI) ')
       call prr ('                       DELTAV = 0.0 for S > SMAX*',
     >cstgrad)
       call prc ('                                         from target')
      ENDIF
c
      CALL PRR  ('  COLLISION ENHANCEMENT FACTOR ZENH         ', CZENH)
      CALL PRB
C
      IF ((CIOPTB.EQ.1.OR.CIOPTB.EQ.2) .AND. IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TAU PARA = MI.SQRT(TB/MB).TI/')
       CALL PRC ('                               (6.8E4.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ENDIF
C-----------------------------------------------------------------------
      IF     (CIOPTC.EQ.0) THEN
       CALL PRC ('  FRICTION OPTION  0 : TAU STOP = MI.TB.SQRT(TB/MB)/')
       CALL PRC ('                       (6.8E4.(1+MB/MI).NB.ZB.ZB.ZI.ZI
     >.ZENH.LAM)')
      ELSEIF (CIOPTC.EQ.1) THEN
       CALL PRC ('  FRICTION OPTION  1 : TAU STOP = INFINITY')
      ELSEIF (CIOPTC.EQ.2) THEN
       CALL PRC ('  FRICTION OPTION  2 : TAU STOP = TAU PARALLEL')
      ELSEIF (CIOPTC.EQ.3) THEN
       CALL PRC ('  FRICTION OPTION  3 : TAU STOP = TI.SQRT(TI.MI)/')
       CALL PRC ('                       (9.0E4.(1+MI/MB).NB.ZB.ZB.ZI.ZI
     >.ZENH.LAM)')
       CALL PRC ('                       WHEN TI > TB.MI/MB AND')
      ELSEIF     (CIOPTC.EQ.4) THEN
       CALL PRC ('  FRICTION OPTION  4 : TAU STOP = MI.TB.SQRT(TB/MB)/')
       CALL PRC ('                       (6.8E4.(1+MB/MI).NB.ZB.ZB.ZI.ZI
     >.ZENH.LAM)')
       call prc ('                       Friction -> 0.0 for a cell if')
       call prc ('                       the target is within one mean f
     >ree path')
      ENDIF
C
      IF ((CIOPTC.NE.0.and.cioptc.ne.4).AND.IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TAU STOP = MI.TB.SQRT(TB/MB)/')
       CALL PRC ('                       (6.8E4.(1+MB/MI).NB.ZB.ZB.ZI.ZI
     >.ZENH.LAM)')
       ENDIF
c
c
c     Reiser's Coulomb Collision Options
c
c
      elseif (cioptr.eq.1) then 

        CALL PRC('  REISER COLL. OPT 1: REISER COULOMB COLLISION OPTION'
     >)
        CALL PRC('                      HAS BEEN TURNED ON')
c
      ELSEIF (CIOPTR.EQ.2) THEN
c
        CALL PRC('  REISER COLL. OPT 2: REISER COULOMB COLLISION OPTION'
     >)
        CALL PRC('                      HAS BEEN TURNED ON')
        call prc('                      REISER COEFFICIENTS ARE RECALCU'
     >//'LATED')
        call prc('                      AT EVERY TIME STEP.')
        call prc('                      WARNING!: MAY BE COMPUTAT'//
     >                                  'IONALLY INTENSIVE')
      endif

c
C-----------------------------------------------------------------------
c
c     Only print out if the Reiser force option is in use 
c
      if (cioptr.eq.1.or.cioptr.eq.2) then 

      IF (Coulomb_log.NE.15.0) THEN

        CALL PRC('  COULOMB LOGARITHM IS NOT ITS USUAL VALUE OF 15')

      ENDIF
c
      IF (ASWITCH.EQ.0) THEN

        CALL PRC('  MASTER SWITCH TO ANALYSE THE INDIVIDUAL FORCE AND')
        CALL PRC('  VELOCITY DIFFUSION COMPONENTS OF THE DRIFT-KINETIC')
        CALL PRC('  MODEL HAS BEEN TURNED "OFF"')

      ELSEIF (ASWITCH.EQ.1) THEN

        CALL PRC('  MASTER SWITCH TO ANALYSE THE INDIVIDUAL FORCE AND')
        CALL PRC('  VELOCITY DIFFUSION COMPONENTS OF THE DRIFT-KINETIC')
        CALL PRC('  MODEL HAS BEEN TURNED "ON". THE FOLLOWING HAVE BEEN'
     >)
        CALL PRC('  SELECTED:')

       IF (SK11.EQ.1) CALL PRC('            FF(DK)')
       IF (SK12.EQ.1) CALL PRC('            FIG(DK)')
       IF (SK13.EQ.1) CALL PRC('            FIV(DK)')
       IF (SD11.EQ.1) CALL PRC('            DIFFUSION: PRESSURE GRADIENT
     >')
       IF (SD12.EQ.1) CALL PRC('            DIFFUSION: TEMPERATURE GRADI
     >ENT')
       IF (SD13.EQ.1) CALL PRC('            DIFFUSION: VELOCITY GRADIENT
     >')

      ENDIF
c
      IF (COMBINE.EQ.1) THEN

        CALL PRC('  -THE FLUID APPROXIMATION AND DRIFT-KINETIC MODELS')
        CALL PRC('   HAVE BEEN COMBINED')

      ENDIF
c
c      IF (VELOVER.EQ.1) THEN
c
c        CALL PRC('  -IMPURITY IONS ARE NO LONGER TRACKED WHEN THEY HEAD'
c     >)
c        CALL PRC('   TOWARD THE OUTER TARGET (DECREASING "S")')
c
c      ELSEIF (VELOVER.EQ.2) THEN
c
c        CALL PRC('  -IMPURITY IONS ARE NO LONGER TRACKED WHEN THEY HEAD'
c     >)
c        CALL PRC('   TOWARD THE INNER TARGET (INCREASING "S")')
c
c      ENDIF
c
c      IF (SWITCHB.EQ.1) CALL PRC('  -PLASMA VALUES HAVE BEEN WRITTEN TO 
c     >THE .LIM FILE')
c      IF (SWITCHZ.EQ.1) CALL PRC('  -FORCE VALUES HAVE BEEN WRITTEN TO T 
c     >HE .LIM FILE')
c      IF (SWITCHV.EQ.1) CALL PRC('  -BINNED VELOCITY VALUES HAVE BEEN WR
c     >ITTEN TO THE .LIM FILE')
c      IF (LINEARPEAK.EQ.1) CALL PRC('  -T-GRAD IS LINEAR ON BOTH SIDES O
c     >F THE MIDPOINT')
c
c
      endif       
c psmod end    

c
C-----------------------------------------------------------------------
c
c      TAG T34 - Delta S steps resulting from 3D Dperp
c
      if (dperpz_opt.eq.0) then 
        call prc('  DPERPZ DELTA S OPTION 0: OFF  (TAG T34)')
        call prc('                           NO S TRANSPORT AS A'//
     >         ' RESULT OF MOTION ALONG THE ORTHOGONAL MAGNETIC AXIS')

      elseif (dperpz_opt.eq.1) then 
        call prc('  DPERPZ DELTA S OPTION 1: ON  (TAG T34)')
        call prc('                           DELTA S STEPS ARE USED'//
     >         ' TO APPROXIMATE THE EFFECT')
        call prc('                           OF DPERP STEPS ALONG THE'//
     >           '"Z" ORTHOGONAL MAGNETIC AXIS')
      endif

C-----------------------------------------------------------------------
      IF     (CIOPTD.EQ.0) THEN
       CALL PRC ('  HEATING OPTION   0 : TAU HEAT = MI.TB.SQRT(TB/MB)/')
       CALL PRC ('                               (1.4E5.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ELSEIF (CIOPTD.EQ.1) THEN
       CALL PRC ('  HEATING OPTION   1 : TAU HEAT = INFINITY')
      ELSEIF (CIOPTD.EQ.2) THEN
       CALL PRC ('  HEATING OPTION   2 : TAU HEAT = ZERO')
      ELSEIF (CIOPTD.EQ.3) THEN
       CALL PRC ('  HEATING OPTION   3 : TAU HEAT =(MI.TB+MB.TI)**3/2 /'
     >)
       CALL PRC ('                       (1.4E5SQRT(MI.MB)NB.ZB.ZB.ZI.ZI
     >.ZENH.LAM)')
      ENDIF
C
      IF (CIOPTD.NE.0.AND.CIOPTD.NE.3.AND.IRSPEC.GT.0) THEN
       CALL PRC (COMENT)
       CALL PRC ('                       TAU HEAT = MI.TB.SQRT(TB/MB)/')
       CALL PRC ('                               (1.4E5.NB.ZB.ZB.ZI.ZI.Z
     >ENH.LAM)')
      ENDIF
c
C-----------------------------------------------------------------------
c
      IF     (CIOPTI.EQ.0) THEN
       CALL PRC ('  CX RECOMB OPTION 0 : NO CHARGE EXCHANGE RECOMBINATIO
     >N')
      ELSEIF (CIOPTI.EQ.1) THEN
       CALL PRC ('  CX RECOMB OPTION 1 : NH=NHO (CONSTANT EVERYWHERE)')
C
C     +NHO.EXP(-X/LAMHX)EXP(-Y/
C     >LAMHY) X>0')
C
C       CALL PRC ('                          NHO.EXP(-Y/LAMHY)
C     >       X<0')
C
       CALL PRC ('                       AND VCX=SQRT(2TB/MB)')
       WRITE (COMENT,'(23X,''WHERE  NHO='',1P,G11.4)') CNHO
C
C       ,'',   NHO='',
C     >   G11.4)') CNHC,CNHO
C
       CALL PRC (COMENT)
C
C       WRITE (COMENT,'(23X,''AND  LAMHX='',1P,G11.4,'', LAMHY='',
C     >   G11.4)') CLAMHX,CLAMHY
C       CALL PRC (COMENT)
C
      ELSEIF (CIOPTI.EQ.2) THEN
       CALL PRC ('  CX RECOMB OPTION 2 : NH=NHO (CONSTANT EVERYWHERE)')
C
C          +NHO.EXP(-X/LAMHX)EXP(-Y/
C          >LAMHY) X>0')
C       CALL PRC ('                          NHO.EXP(-Y/LAMHY)
C     >       X<0')
C
       CALL PRR ('                       WITH CONSTANT VCX =', CVCX)
       WRITE (COMENT,'(23X,''WHERE  NHO='',1P,G11.4)') CNHO
C
C        ,'',   NHO='',
C     >   G11.4)') CNHC,CNHO
C
       CALL PRC (COMENT)
C
C       WRITE (COMENT,'(23X,''AND  LAMHX='',1P,G11.4,'', LAMHY='',
C     >   G11.4)') CLAMHX,CLAMHY
C       CALL PRC (COMENT)
C
      ELSEIF (CIOPTI.EQ.3) THEN
       CALL PRC ('  CX RECOMB OPTION 3 : NH=NHC (CONSTANT IN CORE)')
       CALL PRC ('                       NH=NHO.EXP(-S/LAMHS) (FOR SOL R
     >INGS)')
C
       CALL PRC ('                       AND VCX=SQRT(2TB/MB)')
       WRITE (COMENT,'(23X,''WHERE  NHC='',1P,G11.4,'',   NHO='',
     >   G11.4)') CNHC,CNHO
C
       CALL PRC (COMENT)
C
       WRITE (COMENT,'(23X,''AND  LAMHS='',1P,G11.4)') CLAMHX
       CALL PRC (COMENT)
C
C         ,'', LAMHY='',
C     >   G11.4)') CLAMHX,CLAMHY
C
      ELSEIF (CIOPTI.EQ.4) THEN
       CALL PRC ('  CX RECOMB OPTION 4 : NH FROM NIMBUS')
       CALL PRC ('                       AND VCX=SQRT(2TB/MB)')
C
      ELSEIF (CIOPTI.EQ.5) THEN
       CALL PRC ('  CX RECOMB OPTION 5 : NH FROM NIMBUS')
       CALL PRC ('                       WITH CONSTANT VCX')
      ELSEIF (CIOPTI.EQ.6) THEN
       CALL PRC ('  CX RECOMB OPTION 6 : NH FROM NIMBUS')
       CALL PRC ('                       WITH CCD FROM ADAS')
      ELSEIF (CIOPTI.EQ.7) THEN
       CALL PRC ('  CX RECOMB OPTION 7 : NH FROM OTHER CODE')
       CALL PRC ('                       <SIG V>cx FROM ADPAK/INEL FILES
     >')
       call prc ('                       RATE FILE USED:')
       len = lenstr(mcfile)
       call prc ('                       '//mcfile(1:len))
      ELSEIF (CIOPTI.EQ.8) THEN
       CALL PRC ('  CX RECOMB OPTION 8 : NH FROM NIMBUS OR FLUID CODE')
       CALL PRC ('                       <SIG V>cx FROM CFM PHD THESIS')
       call prc ('                       RATES HAVE BEEN FITTED TO A')
       call prc ('                       THREE PARAMETER EXPONENTIAL FOR
     >MULA')
       call prc ('                       BY TOM ROGNLIEN (LLNL).')
      ELSEIF (CIOPTI.EQ.9) THEN
       CALL PRC ('  CX RECOMB OPTION 9 : NH FROM NIMBUS OR FLUID CODE')
       CALL PRC ('                       <SIG V>cx FROM CFM PHD THESIS')
       call prc ('                       RATES HAVE BEEN FITTED TO A')
       call prc ('                       THREE PARAMETER EXPONENTIAL FOR
     >MULA')
       call prc ('                       BY TOM ROGNLIEN (LLNL).')
       call prc ('                       FITS FOR CII RECOMBINATION AT')
       call prc ('                       LOW TEMPERATURE HAVE BEEN CHANG
     >ED')
       call prc ('                       TO REDUCE THE RATE AND MATCH NE
     >WER')
       call prc ('                       EXPERIMENTAL DATA.')
      ENDIF

      RETURN
      END
c
c
c
      SUBROUTINE PR_sim_options(NIZS,NIMPS,NIMPS2,nymfs)
      use eckstein_2002_yield_data
      use eckstein_2007_yield_data
      IMPLICIT none
c
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS ALL THE SIMULATION OPTION FLAGS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
      include    'pindata'
c      include    'adpak_com'
      include    'promptdep'
      include    'slcom'  
      include    'printopt' 
c
      include    'fperiph_com' 
c
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
      real     totmhpin
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('----- IMPURITY SIMULATION OPTIONS --------- ')
      CALL PRChtml ('--- IMPURITY SIMULATION OPTIONS ---',
     >                'pr_sim_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
      IF     (CNEUTA.EQ.0) THEN
       CALL PRC ('  CONTROL SWITCH   0 : NEUT ON: FOLLOW ATOMS TO IONISA
     >TION POSNS')
      ELSEIF (CNEUTA.EQ.1) THEN
       CALL PRC ('  CONTROL SWITCH   1 : NEUT OFF: INJECT IONS AS "INITI
     >AL STATE"')
      ENDIF
c
      IF     (CNEUTA.EQ.0) THEN
       CALL PRC ('  INJECTION OPT    * : FROM NEUT')
      ELSEIF (CIOPTE.EQ.1) THEN
       CALL PRR ('  INJECTION OPT    1 : INJECT AT GIVEN (R,Z), V0=',
     >   CTEM1)
      ELSEIF (CIOPTE.EQ.2) THEN
       CALL PRI ('  INJECTION OPT    2 : INJECT ON RING      : ',INJIR)
       CALL PRR2('                       BETWEEN FACTORS     : ',INJF1,
     >       INJF2)
       CALL PRR ('                       TIMES SMAX FOR RING : ',
     >       KSMAXS(INJIR))
       CALL PRC ('                       FROM BOTH PLATES. ')
       CALL PRR ('                       VALUE FOR V0 (TEMP) = ',
     >     CTEM1)
c       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
c     >     1.56E4*SQRT(CTEM1/CRMI))
       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
     >     9.79e3*SQRT(CTEM1/CRMI))
      ELSEIF (CIOPTE.EQ.3) THEN
       CALL PRI ('  INJECTION OPT    3 : INJECT ON RING      : ',INJIR)
       CALL PRR2('                       BETWEEN FACTORS     : ',INJF1,
     >       INJF2)
       CALL PRR ('                       TIMES SMAX FOR RING : ',
     >       KSMAXS(INJIR))
       CALL PRR ('                       VALUE FOR V0 (TEMP) = ',
     >     CTEM1)
c       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
c     >     1.56E4*SQRT(CTEM1/CRMI))
       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
     >     9.79e3*SQRT(CTEM1/CRMI))
      ELSEIF (CIOPTE.EQ.4) THEN
       CALL PRC ('  INJECTION OPT    4 : Neutral impurity ionization pro
     >files')
       call prc ('                       taken from a HYDROGEN CODE run
     >are used')
       call prc ('                       to generate a probability map f
     >or ion')
       call prc ('                       injection. The initial ion ener
     >gy is')
       call prc ('                       taken from the HYDROGEN CODE re
     >sults')


       if ((pincode.eq.0.and.cneur.eq.5.and.ctrap.eq.5).or.
     >      pincode.eq.2.or.pincode.eq.3.or.pincode.eq.4.or.
     >      pincode.eq.5) then
c
c       Print out information about the PIN wall-target source.
c
        call prb
        call prc('   PIN Source characteristics:')

        CALL PRC('   SAMPLE PRIMARY FLUX AND YIELD DATA')
c
        call prb
c
        totfpin = 0.0
        totfypin= 0.0
        totzfpin= 0.0
        tothpin = 0.0
        totmhpin= 0.0
        totfapin= 0.0
c
        if (wallpts.gt.0) write(7,9001)
c
        DO ID = 1, wallpts
c
           IN = WALLPT(id,17)
c
           if (flxhw2(in).le.0.0) then
              tmpy = 0.0
           else
              tmpy = flxhw3(in)/flxhw2(in)
           endif
c
           WRITE (7,9003) id,
     >         wallpt(ID,1),wallpt(id,2),fluxhw(in),flxhw2(in),
     >         flxhw5(in),flxhw3(in),tmpy,
     >         flxhw4(in),wallpt(id,16)
c
           totfpin = totfpin + flxhw2(in) * wallpt(id,7)
           totfypin = totfypin + flxhw3(in) * wallpt(id,7)
           totzfpin = totzfpin + flxhw4(in) * wallpt(id,7)
           totfapin = totfapin + flxhw6(in) * wallpt(id,7)
           tothpin = tothpin + flxhw6(in) * flxhw5(in) * wallpt(id,7)
           totmhpin =totmhpin + (fluxhw(in)-flxhw6(in))* kboltz
     >                          * wallpt(id,19) * wallpt(id,7)
c
        end do
c
         call prb
c
         call prc(' CALCULATED FROM SEGMENT DATA:')
         CALL PRR('   TOTAL PRIMARY INTEGRATED ATOM+ION FLUX',
     >                totfpin)
         CALL PRR('   TOTAL PRIMARY INTEGRATED FLUX*YIELD   ',
     >                totfypin)
         CALL PRR('   TOTAL INTEGRATED Z-REDEPOSITION       ',
     >                totzfpin)
         call prr('   TOTAL ATOMIC H HEAT FLUX TO SURFACES  ',
     >                tothpin*ech)
         call prr('   ATOM HEAT+REC POWER FLUX TO SURFACES  ',
     >                tothpin*ech+2.2*ech*totfapin)
         call prr('   MOL H HEAT FLUX TO SURFACES (E=SURF T)',
     >                totmhpin)
c
         if (totfpin.gt.0.0) then
         call prr('   TOTAL EFFECTIVE YIELD               ',
     >                totfypin/totfpin)
         endif
c
         if (pincode.eq.0) then

            call prb
            call prc(' PIN/NIMBUS REPORTED TOTALS: ')
            call prr('   TOTAL NUMBER OF SPUTTERED IMPURITIES',
     >                zsput)
            call prr('   TOTAL NEUTRAL SPUTTERED IMPURITIES  ',
     >                zsputn)
            call prr('   IMPURITIES LOST TO CENTRAL ESCAPE   ',
     >                zescpd)
            call prr('   IMPURITIES LOST TO ALBEDO REGIONS   ',
     >                zescal)
            call prr('   IMPURITIES LOST TO LEAK REGIONS     ',
     >                zesclk)
         endif
c
         call prb
c
c        Format statements - for PIN output
c

 9001 FORMAT(' IND',3x,'R',5x,'Z',5x,'AT+MOL.F',2x,
     >      'AT+ION.F',1x,'AT.ENERGY',1x,'Z-FLUX',2x,
     >      'EFF-YIELD',3x,'Z-REDEP',3x,'TYPE')
 9003 FORMAT(I4,2(1x,F5.2),2(1x,g9.2),1x,f8.3,1x,G9.2,
     >       1x,g9.2,1x,g9.2,1x,f5.2)
c
c
c
       else
        call prb
        call prc('                       PIN  wall options NOT specifi
     >ed')
        call prc('                       Detailed information for WALL
     >source')
        call prc('                       is NOT printed. Specify - ')
        call prc('                       NEUTRAL WALL OPTION = 5 and')
        call prc('                       TRAP WALL OPTION    = 5')
        call prc('                       To get more detail.')
        call prb
       endif


      ELSEIF (CIOPTE.EQ.5) THEN
       CALL PRI ('  INJECTION OPT    5 : INJECT ON RING      : ',INJIR)
       CALL PRR2('                       BETWEEN FACTORS     : ',INJF1,
     >       INJF2)
       CALL PRR ('                       TIMES SMAX FOR RING : ',
     >       KSMAXS(INJIR))
       CALL PRC ('                       FROM BOTH PLATES. ')
       call prc ('                       Initial Velocity is calculated
     >from:')
       call prc ('                       Vinit =  Rg * V0')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       CALL PRR ('                       VALUE FOR V0 (TEMP) = ',
     >     CTEM1)
c       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
c     >     1.56E4*SQRT(CTEM1/CRMI))
       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
     >     9.79e3*SQRT(CTEM1/CRMI))
      ELSEIF (CIOPTE.EQ.6) THEN
       CALL PRI ('  INJECTION OPT    6 : INJECT ON RING      : ',INJIR)
       CALL PRR2('                       BETWEEN FACTORS     : ',INJF1,
     >       INJF2)
       CALL PRR ('                       TIMES SMAX FOR RING : ',
     >       KSMAXS(INJIR))
       call prc ('                       Initial Velocity is calculated
     >from:')
       call prc ('                       Vinit =  Rg * V0')
       call prc ('                       Where Rg = SQRT(-2*ln(X1))* cos
     >(2PI*X2) ')
       call prc ('                       X1,X2 are uniform on [0,1] ')
       CALL PRR ('                       VALUE FOR V0 (TEMP) = ',
     >     CTEM1)
c       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
c     >     1.56E4*SQRT(CTEM1/CRMI))
       CALL PRR ('                       VALUE FOR V0 (M/S)  = ',
     >     9.79e3*SQRT(CTEM1/CRMI))
c
      ELSEIF (CIOPTE.EQ.7) THEN
c
       CALL PRC ('  INJECTION OPT    7 : Neutral impurity ionization sou
     >rce profiles')
       call prc ('                       taken from a FLUID CODE run are
     > used')
       call prc ('                       to generate a probability map f
     >or ion')
       call prc ('                       injection. The initial ion ener
     >gy is')
       call prc ('                       taken from the FLUID CODE data'
     >)
c
      ELSEIF (CIOPTE.EQ.8) THEN
c
       CALL PRC ('  INJECTION OPT    8 : Neutral impurity ionization sou
     >rce profiles')
       call prc ('                       taken from a FLUID CODE run are
     > used')
       call prc ('                       to generate a probability map f
     >or ion')
       call prc ('                       injection. The initial ion ener
     >gy is')
       call prc ('                       taken from the FLUID CODE data.
     >')
       call prc ('                       Initial particle velocity is Ma
     >xwellian')
c slmod begin
      ELSEIF (CIOPTE.EQ.9) THEN
       CALL PRR ('  INJECTION OPT    9 : LINE INJECTION (UNIFORM)  V0=',
     >   CTEM1)

       WRITE(DATUNIT,'(7X,A,2F10.3)') 'END POINT A (R,Z)=',CXSCA,CYSCA
       WRITE(DATUNIT,'(7X,A,2F10.3)') 'END POINT B (R,Z)=',CXSCB,CYSCB
c slmod end
c jdemod begin
      ELSEIF (CIOPTE.EQ.10) THEN
       CALL PRR ('  INJECTION OPT   10 : BOX INJECTION (UNIFORM)  V0=',
     >   CTEM1)
       call prc ('                       PARTICLES RANDOMLY DISTRIBUTED
     > OVER THE BOX') 
       call prc ('                       DEFINED BY THE CORNER POINTS')

       WRITE(DATUNIT,'(7X,A,2F10.3)') 'END POINT A (R,Z)=',CXSCA,CYSCA
       WRITE(DATUNIT,'(7X,A,2F10.3)') 'END POINT B (R,Z)=',CXSCB,CYSCB
c jdemod end
      ENDIF

C-----------------------------------------------------------------------
c
      if (cneuta.eq.1.and.ciopte.eq.4) then
c
        call prc ('  PIN RUN          * : PIN IS RUN FOR INJECTION'//
     >                                  ' OPTION 4.')
        call prb
c
        if (cpinopt.eq.0) then
c
c
c         Print analysis of the PIN input passed through DIVIMP
c
          if (pincode.eq.0) then
c
               call prc ('                       '//
     >             'NIMBUS called to generate impurity neutral'//
     >                  ' distribution')
             call prnimbin(cnimbin,nlines)
          elseif (pincode.eq.1.or.pincode.eq.2.or.pincode.eq.3.or.
     >            pincode.eq.4.or.pincode.eq.5) then
                 call prc ('                       '//
     >            'EIRENE called to generate impurity neutral'//
     >                  ' distribution')
          endif
c
        endif
c
      endif

c
c     Branch around launch and sputter code for injection cases  
c
      if (cneuta.eq.1) goto 997  

c
C-----------------------------------------------------------------------
      IF     (CNEUTB.EQ.0) THEN
       CALL PRC ('  LAUNCH OPTION    0 : DISTRIBUTED LAUNCH ALONG TARGET
     >')
      ELSEIF (CNEUTB.EQ.1) THEN
       CALL PRC ('  LAUNCH OPTION    1 : AT GIVEN (R,Z)')
       CALL PRR2('                       COORDINATES: ',CXSC,CYSC)
      ELSEIF (CNEUTB.EQ.2) THEN
       CALL PRC ('  LAUNCH OPTION    2 : HOMOGENEOUSLY ALONG WALLS')
       CALL PRR ('                       IMPACT NEUTRAL ENERGY (EV) = ',
     >             CEIMP)
       CALL PRC ('                       WALL SEGMENT LAUNCH PROB'//
     >           'ABILITY')
       IF (WLPABS.EQ.0) THEN
         CALL PRC ('                       MODIFIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.1) THEN
         CALL PRC ('                       SET TO THESE VALUES')
       ELSEIF (WLPABS.EQ.2) THEN
         CALL PRC ('                       READ FROM PIN DATA AND')
         CALL PRC ('                       MULTIPLIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.3) THEN
         CALL PRC ('                       READ FROM PIN DATA')
       ENDIF
       IF (WLPABS.NE.3) THEN
         CALL PRC ('                       INDEX 1    INDEX 2   PROB MUL
     >T.')
         DO 400 I = 1,NWLPROB
         CALL PRR3 ('                      ',WLPROB(I,1),WLPROB(I,2),
     >             WLPROB(I,3))
400      CONTINUE
       ENDIF
      elseIF    (CNEUTB.EQ.3) THEN
       CALL PRC ('  LAUNCH OPTION    3 : DISTRIBUTED LAUNCH ALONG TARGET
     >')
       call prc ('                       DUE TO ION IMPACT.')
       call prc ('                       USING PIN/NIMBUS TARGET DATA')
      ELSEIF (CNEUTB.EQ.4) THEN
       CALL PRC ('  LAUNCH OPTION    4 : DISTRIBUTED ALONG ''WALLS'' ')
       call prc ('                       DUE TO ATOM IMPACT')
       call prc ('                       WALLS ARE ALL VESSEL SURFACES')
       call prc ('                       INCLUDING TARGET SEGMENTS')
       call prc ('                       USING PIN/NIMBUS WALL DATA')
       call prb
c
       if (nwlprob.gt.0) then
c
       CALL PRC ('                       WALL SEGMENT LAUNCH PROB'//
     >           'ABILITY')
       IF (WLPABS.EQ.0) THEN
         CALL PRC ('                       MODIFIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.1) THEN
         CALL PRC ('                       SET TO THESE VALUES')
       ELSEIF (WLPABS.EQ.2) THEN
         CALL PRC ('                       READ FROM PIN DATA AND')
         CALL PRC ('                       MULTIPLIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.3) THEN
         CALL PRC ('                       READ FROM PIN DATA')
       ENDIF
       IF (WLPABS.NE.3) THEN
         CALL PRC ('                       INDEX 1    INDEX 2   PROB MUL
     >T.')
         DO  I = 1,NWLPROB
         CALL PRR3 ('                      ',WLPROB(I,1),WLPROB(I,2),
     >             WLPROB(I,3))
         end do
       ENDIF
       endif
      ELSEIF (CNEUTB.EQ.5) THEN
       CALL PRC ('  LAUNCH OPTION    5 : 2D NEUTRAL LAUNCH - NEUTRALS')
       CALL PRC ('                       ARE LAUNCHED FROM THE CENTRE OF
     >')
       CALL PRC ('                       CELLS ON THE (IK,IR) GRID. THE
     >')
       CALL PRC ('                       PROBABILITY FOR LAUNCH IS GIVEN
     >')
       CALL PRC ('                       BY AN EXTERNALLY SUPPLIED SOURC
     >E.')
       CALL PRC ('                       PROBABILITY IS PROPORTIONAL TO
     >THE')
       CALL PRC ('                       SOURCE RATE IN EACH CELL.')

      ELSEIF (CNEUTB.EQ.6) THEN
       CALL PRC ('  LAUNCH OPTION    6 : FREESPACE LAUNCH')
       call prc ('                       ON LINE JOINING P1 P2')
       CALL PRR2('                       P1 = ',CXSCA,CYSCA)
       CALL PRR2('                       P2 = ',CXSCA,CYSCA)
      ELSEIF (CNEUTB.EQ.7) THEN
       CALL PRC ('  LAUNCH OPTION    7 : FREESPACE LAUNCH')
       call prc ('                       WITHIN THE BOX DEFINED'//
     >                                 ' BY P1 AND P2')
       CALL PRR2('                       P1 = ',CXSCA,CYSCA)
       CALL PRR2('                       P2 = ',CXSCA,CYSCA)
      ENDIF
C-----------------------------------------------------------------------
c    
c     Print notification if Far Periphery neutral ionization option is
c     ON. 
c
      if (fp_neut_opt.gt.0) then 
        call prb   
        call prc('  FAR PERIPHERY NEUTRAL IONIZATION OPTION : ON')
        call prc('  - SEE FAR PERIPHERY DATA PRINT OUT FOR DETAILS')
        call prb
      endif

C-----------------------------------------------------------------------
C
C------- INITIAL NEUTRAL V/A FLAG --------------------------------------
C
C  THIS OPTION IS SET EQUAL TO CNEUTC, IF REQ'D, IN READIN
C  IMMEDIATELY AFTER THE OTHER LAUNCH OPTIONS HAVE BEEN ASSIGNED
C
c      IF (NVAOPT.NE.CNEUTC)  THEN
      if (cneuta.eq.0) then
c
c     Initial Neutral V/A
c
       call prb
       CALL PRC ('  INITAL V/A FLAG    : THIS OPTION IS USED TO SPECIFY'
     >)
       CALL PRC ('                       THE V/A FLAG THAT WILL BE USED
     >FOR THE')
       CALL PRC ('                       VERY FIRST BATCH OF INITIALLY R
     >ELEASED')
       CALL PRC ('                       NEUTRAL PARTICLES. THE REGULAR
     >V/A FLAG')
       CALL PRC ('                       IS THEN USED FOR EACH SUBSEQUEN
     >T LAUNCH')
       CALL PRC ('                       OF SELF-SPUTTERED NEUTRALS.')
       call prb
c
      IF     (NVAOPT.EQ.0) THEN
       CALL PRC ('  INITIAL V/A FLAG 0 : THETA =+/-ASIN($), $ IN (0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (NVAOPT.EQ.1) THEN
       CALL PRC ('  INITIAL V/A FLAG 1 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NVAOPT.EQ.2) THEN
       CALL PRC ('  INITIAL V/A FLAG 2 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(
     >1-$))),')
       CALL PRR ('                             $ IN (0,1),  TG=',CTEM1)
      ELSEIF (NVAOPT.EQ.3) THEN
       CALL PRC ('  INITIAL V/A FLAG 3 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.4) THEN
       CALL PRC ('  INITIAL V/A FLAG 4 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (NVAOPT.EQ.5) THEN
       CALL PRC ('  INITIAL V/A FLAG 5 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (NVAOPT.EQ.6) THEN
       CALL PRC ('  INITIAL V/A FLAG 6 : THETA = 0')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.7) THEN
       CALL PRC ('  INITIAL V/A FLAG 7 : THETA =+/-ACOS((1-$)**1/3), $ I
     >N (0,1)   "FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.8) THEN
       CALL PRC ('  INITIAL V/A FLAG 8 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.9) THEN
       CALL PRC ('  INITIAL V/A FLAG 9 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO
     >VALUES')
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-
     >')
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')
     >   CTEM1,CTEM2
       CALL PRC (COMENT)
      ELSEIF (NVAOPT.EQ.10) THEN
       CALL PRC ('  INITIAL V/A FLAG 10: BETA = ACOS((1-$)**1/3)  "3D FR
     >EE JET"')
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.11) THEN
       CALL PRC ('  INITIAL V/A FLAG 11: THETA =+/-ACOS((1-$)**1/3), 0<$
     ><1   "2.5D FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NVAOPT.EQ.12) THEN
       CALL PRC ('  INITIAL V/A FLAG 12: EMISSION AT CONSTANT ENERGY, EI
     >N,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       EIN=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (NVAOPT.EQ.13) THEN
       CALL PRC ('  INITIAL V/A FLAG 13: EMISSION AT TEMPERATURE, TG,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       TG=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (NVAOPT.EQ.14) THEN
       CALL PRC ('  INITIAL V/A FLAG 14: THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = 1.38E4 * SQRT(Etarg/MI)')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (NVAOPT.EQ.15) THEN
       CALL PRC ('  INITIAL V/A FLAG 15: THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRc ('                       VIN = SQRT(2TiB/MI)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (NVAOPT.EQ.16) THEN
       CALL PRC ('  INITIAL V/A FLAG 16: THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NVAOPT.EQ.17) THEN
       CALL PRC ('  INITIAL V/A FLAG 17: THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (NVAOPT.EQ.18) THEN
       CALL PRC ('  INITIAL V/A FLAG 18: THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NVAOPT.EQ.19) THEN
       CALL PRC ('  INITIAL V/A FLAG 19: THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       BETA = 2PI$, $ in [0,1]')
       call prc ('                       PHI  = 2PI$ - PI, $ in [0,1]')
      ENDIF
      ENDIF


997   continue


C-----------------------------------------------------------------------
C
C BRANCH AROUND THE SPUTTER OPTION AND V/A FOR LAUNCH OPTION 1 AND
C INJECTION CASES IF THE
C OPTION SPECIFIED DOES NOT INCLUDE SELF-SPUTTERING
C
      IF (CNEUTB.EQ.1.OR.CNEUTH.EQ.1.OR.CNEUTA.EQ.1.or.
     >    cneutb.eq.6.or.cneutb.eq.7) THEN
         IF (cselfs.eq.0) THEN
            GOTO 998
         ELSE
           CALL PRB
           CALL PRC ('  SPUTTER OPTION AND V/A FLAG APPLY ONLY TO ')
           CALL PRC ('  PLATE LAUNCHED NEUTRALS (IF ANY) AND ')
           CALL PRC ('  SELF-SPUTTERING.')
           CALL PRB
         ENDIF
      ENDIF
C-----------------------------------------------------------------------
      IF     (CNEUTC.EQ.0) THEN
       CALL PRC ('  VEL/ANGLE FLAG   0 : THETA =+/-ASIN($), $ IN (0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (CNEUTC.EQ.1) THEN
       CALL PRC ('  VEL/ANGLE FLAG   1 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTC.EQ.2) THEN
       CALL PRC ('  VEL/ANGLE FLAG   2 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(
     >1-$))),')
       CALL PRR ('                             $ IN (0,1),  TG=',CTEM1)
      ELSEIF (CNEUTC.EQ.3) THEN
       CALL PRC ('  VEL/ANGLE FLAG   3 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.4) THEN
       CALL PRC ('  VEL/ANGLE FLAG   4 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (CNEUTC.EQ.5) THEN
       CALL PRC ('  VEL/ANGLE FLAG   5 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (CNEUTC.EQ.6) THEN
       CALL PRC ('  VEL/ANGLE FLAG   6 : THETA = 0')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.7) THEN
       CALL PRC ('  VEL/ANGLE FLAG   7 : THETA =+/-ACOS((1-$)**1/3), $ I
     >N (0,1)   "FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.8) THEN
       CALL PRC ('  VEL/ANGLE FLAG   8 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.9) THEN
       CALL PRC ('  VEL/ANGLE FLAG   9 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO
     >VALUES')
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-
     >')
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')
     >   CTEM1,CTEM2
       CALL PRC (COMENT)
      ELSEIF (CNEUTC.EQ.10) THEN
       CALL PRC ('  VEL/ANGLE FLAG  10 : BETA = ACOS((1-$)**1/3)  "3D FR
     >EE JET"')
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.11) THEN
       CALL PRC ('  VEL/ANGLE FLAG  11 : THETA =+/-ACOS((1-$)**1/3), 0<$
     ><1   "2.5D FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTC.EQ.12) THEN
       CALL PRC ('  VEL/ANGLE FLAG  12 : EMISSION AT CONSTANT ENERGY, EI
     >N,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       EIN=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (CNEUTC.EQ.13) THEN
       CALL PRC ('  VEL/ANGLE FLAG  13 : EMISSION AT TEMPERATURE, TG,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       TG=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (CNEUTC.EQ.14) THEN
       CALL PRC ('  VEL/ANGLE FLAG  14 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = 1.38E4 * SQRT(Etarg/MI)')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (CNEUTC.EQ.15) THEN
       CALL PRC ('  VEL/ANGLE FLAG  15 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRc ('                       VIN = SQRT(2TiB/MI)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (CNEUTC.EQ.16) THEN
       CALL PRC ('  VEL/ANGLE FLAG  16 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTC.EQ.17) THEN
       CALL PRC ('  VEL/ANGLE FLAG  17 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (CNEUTC.EQ.18) THEN
       CALL PRC ('  VEL/ANGLE FLAG  18 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTC.EQ.19) THEN
       CALL PRC ('  VEL/ANGLE FLAG  19 : 3D ISOTROPIC')
       call prc ('                       THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRR ('                       VIN = SQRT(2EIN/MI)'//
     >                '*SQRT(|COS(B)**2+SIN(B)**2.COS(PHI)**2|)')
       call prr ('                       EIN = ',ctem1)

       call prc ('                       BETA = 2PI$, $ in [0,1]')
       call prc ('                       PHI  = 2PI$ - PI, $ in [0,1]')
      ENDIF
C-----------------------------------------------------------------------
C
C     SUPPLEMENTAL LAUNCH OPTIONS ARE ONLY PRINTED IF SOME NUMBER
C     OF SUPPLEMENTAL LAUNCHES (NIMPS2 > 0) HAVE BEEN SPECIFIED
C
C-----------------------------------------------------------------------
      IF (NIMPS2.GT.0) THEN
C
      IF     (CNEUTH.EQ.0) THEN
       CALL PRC ('  SUP LAUNCH OPTION 0: DISTRIBUTED LAUNCH ALONG TARGET
     >')
      ELSEIF (CNEUTH.EQ.1) THEN
       CALL PRC ('  SUP LAUNCH OPTION 1: AT GIVEN (R,Z)')
       CALL PRR2('                       COORDINATES: ',CXSC,CYSC)
      ELSEIF (CNEUTH.EQ.2) THEN
       CALL PRC ('  SUP LAUNCH OPTION 2: HOMOGENEOUSLY ALONG WALLS')
       CALL PRR ('                       IMPACT NEUTRAL ENERGY (EV) = ',
     >             CEIMP)
       CALL PRC ('                       WALL SEGMENT LAUNCH PROB'//
     >           'ABILITY')
       IF (WLPABS.EQ.0) THEN
         CALL PRC ('                       MODIFIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.1) THEN
         CALL PRC ('                       SET TO THESE VALUES')
       ELSEIF (WLPABS.EQ.2) THEN
         CALL PRC ('                       READ FROM PIN DATA AND')
         CALL PRC ('                       MULTIPLIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.3) THEN
         CALL PRC ('                       READ FROM PIN DATA')
       ENDIF
       IF (WLPABS.NE.3) THEN
         CALL PRC ('                       INDEX 1    INDEX 2    PROB MU
     >LT.')
         DO 500 I = 1,NWLPROB
         CALL PRR3 ('                      ',WLPROB(I,1),WLPROB(I,2),
     >             WLPROB(I,3))
500      CONTINUE
       ENDIF
      elseIF    (CNEUTH.EQ.3) THEN
       CALL PRC ('  SUP LAUNCH OPTION 3: DISTRIBUTED LAUNCH ALONG TARGET
     >')
       call prc ('                       DUE TO ION IMPACT.')
       call prc ('                       USING PIN/NIMBUS TARGET DATA')
      ELSEIF (CNEUTH.EQ.4) THEN
       CALL PRC ('  SUP LAUNCH OPTION 4: DISTRIBUTED ALONG ''WALLS'' ')
       call prc ('                       DUE TO ATOM IMPACT')
       call prc ('                       WALLS ARE ALL VESSEL SURFACES')
       call prc ('                       INCLUDING TARGET SEGMENTS')
       call prc ('                       USING PIN/NIMBUS WALL DATA')
c
       if (nwlprob.gt.0) then
c
       call prb
       CALL PRC ('                       WALL SEGMENT LAUNCH PROB'//
     >           'ABILITY')
       IF (WLPABS.EQ.0) THEN
         CALL PRC ('                       MODIFIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.1) THEN
         CALL PRC ('                       SET TO THESE VALUES')
       ELSEIF (WLPABS.EQ.2) THEN
         CALL PRC ('                       READ FROM PIN DATA AND')
         CALL PRC ('                       MULTIPLIED BY THE FOLLOWING')
       ELSEIF (WLPABS.EQ.3) THEN
         CALL PRC ('                       READ FROM PIN DATA')
       ENDIF
       IF (WLPABS.NE.3) THEN
         CALL PRC ('                       INDEX 1    INDEX 2   PROB MUL
     >T.')
         DO  I = 1,NWLPROB
         CALL PRR3 ('                      ',WLPROB(I,1),WLPROB(I,2),
     >             WLPROB(I,3))
         end do
       ENDIF
       endif
C
      ELSEIF (CNEUTH.EQ.5) THEN
       CALL PRC ('  SUP LAUNCH OPTION 5: 2D NEUTRAL LAUNCH - NEUTRALS')
       CALL PRC ('                       ARE LAUNCHED FROM THE CENTRE OF
     >')
       CALL PRC ('                       CELLS ON THE (IK,IR) GRID. THE
     >')
       CALL PRC ('                       PROBABILITY FOR LAUNCH IS GIVEN
     >')
       CALL PRC ('                       BY AN EXTERNALLY SUPPLIED SOURC
     >E.')
       CALL PRC ('                       PROBABILITY IS PROPORTIONAL TO
     >THE')
       CALL PRC ('                       SOURCE RATE IN EACH CELL.')
c
      ENDIF
C-----------------------------------------------------------------------
       call prb
       CALL PRC ('  SUP V/A FLAG       : THIS OPTION IS USED TO SPECIFY'
     >)
       CALL PRC ('                       THE V/A FLAG THAT WILL BE USED
     >FOR THE')
       CALL PRC ('                       FIRST BATCH OF SUPPLEMENTAL NEU
     >TRAL')
       CALL PRC ('                       PARTICLES LAUNCHED. THE REGULAR
     > V/A FLAG')
       CALL PRC ('                       IS THEN USED FOR EACH SUBSEQUEN
     >T LAUNCH')
       CALL PRC ('                       OF SELF-SPUTTERED NEUTRALS.')
       call prb
c
      IF     (CNEUTI.EQ.0) THEN
       CALL PRC ('  SUP V/A FLAG     0 : THETA =+/-ASIN($), $ IN (0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (CNEUTI.EQ.1) THEN
       CALL PRC ('  SUP V/A FLAG     1 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTI.EQ.2) THEN
       CALL PRC ('  SUP V/A FLAG     2 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(
     >1-$))),')
       CALL PRR ('                             $ IN (0,1),  TG=',CTEM1)
      ELSEIF (CNEUTI.EQ.3) THEN
       CALL PRC ('  SUP V/A FLAG     3 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.4) THEN
       CALL PRC ('  SUP V/A FLAG     4 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (CNEUTI.EQ.5) THEN
       CALL PRC ('  SUP V/A FLAG     5 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (CNEUTI.EQ.6) THEN
       CALL PRC ('  SUP V/A FLAG     6 : THETA = 0')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.7) THEN
       CALL PRC ('  SUP V/A FLAG     7 : THETA =+/-ACOS((1-$)**1/3), $ I
     >N (0,1)   "FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.8) THEN
       CALL PRC ('  SUP V/A FLAG     8 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.9) THEN
       CALL PRC ('  SUP V/A FLAG     9 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO
     >VALUES')
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-
     >')
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')
     >   CTEM1,CTEM2
       CALL PRC (COMENT)
      ELSEIF (CNEUTI.EQ.10) THEN
       CALL PRC ('  SUP V/A FLAG    10 : BETA = ACOS((1-$)**1/3)  "3D FR
     >EE JET"')
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.11) THEN
       CALL PRC ('  SUP V/A FLAG    11 : THETA =+/-ACOS((1-$)**1/3), 0<$
     ><1   "2.5D FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (CNEUTI.EQ.12) THEN
       CALL PRC ('  SUP V/A FLAG    12 : EMISSION AT CONSTANT ENERGY, EI
     >N,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       EIN=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (CNEUTI.EQ.13) THEN
       CALL PRC ('  SUP V/A FLAG    13 : EMISSION AT TEMPERATURE, TG,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       TG=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (CNEUTI.EQ.14) THEN
       CALL PRC ('  SUP V/A FLAG    14 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = 1.38E4 * SQRT(Etarg/MI)')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (CNEUTI.EQ.15) THEN
       CALL PRC ('  SUP V/A FLAG    15 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRc ('                       VIN = SQRT(2TiB/MI)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (CNEUTI.EQ.16) THEN
       CALL PRC ('  SUP V/A FLAG    16 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTI.EQ.17) THEN
       CALL PRC ('  SUP V/A FLAG    17 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (CNEUTI.EQ.18) THEN
       CALL PRC ('  SUP V/A FLAG    18 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (CNEUTI.EQ.19) THEN
       CALL PRC ('  SUP V/A FLAG    19 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       BETA = 2PI$, $ in [0,1]')
       call prc ('                       PHI  = 2PI$ - PI, $ in [0,1]')
      ENDIF
C-----------------------------------------------------------------------
      if (neut2d_opt.eq.0) then
       CALL PRC ('  2D NEUTRAL SOURCE 0: THE 2D NEUTRAL SOURCE OPTION IS
     > OFF')
       CALL PRC ('                       NO NEUTRALS ARE LAUNCHED PROPOR
     >TIONAL TO')
       CALL PRC ('                       A GRID BASED 2D DISTRIBUTION.')
      ELSEIF(NEUT2D_OPT.EQ.1) THEN
       CALL PRC ('  2D NEUTRAL SOURCE 1: THE 2D NEUTRAL SOURCE OPTION IS
     > ON')
       CALL PRC ('                       NEUTRALS ARE LAUNCHED PROPORTIO
     >NAL TO')
       CALL PRC ('                       TO A GIVEN 2D DISTRIBUTION. THE
     >')
       CALL PRC ('                       DISTRIBUTION AND STRENGTH IN TH
     >IS')
       CALL PRC ('                       CASE ARE CALCULATED FROM THE UE
     >DGE RECOMBINATION')
       CALL PRC ('                       SOURCE IN EACH CELL.')
      ENDIF
C
C     PRINT OUT CORRESPONDING V/A FLAG
C
      IF (NEUT2D_OPT.NE.0.0) THEN
C
      IF     (NEUT2D_VAOPT.EQ.0) THEN
       CALL PRC ('  2D NEUTRAL V/A   0 : THETA =+/-ASIN($), $ IN (0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (NEUT2D_VAOPT.EQ.1) THEN
       CALL PRC ('  2D NEUTRAL V/A   1 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NEUT2D_VAOPT.EQ.2) THEN
       CALL PRC ('  2D NEUTRAL V/A   2 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2TG/MI).SQRT(ABS(LN(
     >1-$))),')
       CALL PRR ('                             $ IN (0,1),  TG=',CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.3) THEN
       CALL PRC ('  2D NEUTRAL V/A   3 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.4) THEN
       CALL PRC ('  2D NEUTRAL V/A   4 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (NEUT2D_VAOPT.EQ.5) THEN
       CALL PRC ('  2D NEUTRAL V/A   5 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EBD/MI(1/SQRT($)-1)
     >) $<$MAX')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (NEUT2D_VAOPT.EQ.6) THEN
       CALL PRC ('  2D NEUTRAL V/A   6 : THETA = 0')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.7) THEN
       CALL PRC ('  2D NEUTRAL V/A   7 : THETA =+/-ACOS((1-$)**1/3), $ I
     >N (0,1)   "FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.8) THEN
       CALL PRC ('  2D NEUTRAL V/A   8 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.9) THEN
       CALL PRC ('  2D NEUTRAL V/A   9 : THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRC ('                       VIN = SQRT(2EIN/MI), WHERE TWO
     >VALUES')
       CALL PRC ('                       ARE USED ALTERNATELY FOR EIN :-
     >')
       WRITE (COMENT,'(23X,''EIN1 ='',1P,G11.4,''  AND EIN2 ='',G11.4)')
     >   CTEM1,CTEM2
       CALL PRC (COMENT)
      ELSEIF (NEUT2D_VAOPT.EQ.10) THEN
       CALL PRC ('  2D NEUTRAL V/A  10 : BETA = ACOS((1-$)**1/3)  "3D FR
     >EE JET"')
       CALL PRC ('                       PSI = 2PI$,  $ IN (0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.11) THEN
       CALL PRC ('  2D NEUTRAL V/A  11 : THETA =+/-ACOS((1-$)**1/3), 0<$
     ><1   "2.5D FREE JET"')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
      ELSEIF (NEUT2D_VAOPT.EQ.12) THEN
       CALL PRC ('  2D NEUTRAL V/A  12 : EMISSION AT CONSTANT ENERGY, EI
     >N,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       EIN=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (NEUT2D_VAOPT.EQ.13) THEN
       CALL PRC ('  2D NEUTRAL V/A  13 : EMISSION AT TEMPERATURE, TG,')
       CALL PRC ('                       INTO A COS**N DISTRIBUTION')
       CALL PRR ('                       TG=',CTEM1)
       CALL PRR ('                       N=',CNIN)
      ELSEIF (NEUT2D_VAOPT.EQ.14) THEN
       CALL PRC ('  2D NEUTRAL V/A  14 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       VIN = 1.38E4 * SQRT(Etarg/MI)')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (NEUT2D_VAOPT.EQ.15) THEN
       CALL PRC ('  2D NEUTRAL V/A  15 : THETA = 2PI$,  $ IN (0,1)  "ISO
     >TROPIC"')
       CALL PRc ('                       VIN = SQRT(2TiB/MI)')
       call prc ('                       IN ADDITION VIN IS MULTIPLIED B
     >Y:')
       call prr ('                       VIN = VIN * ',cvamult)
      ELSEIF (NEUT2D_VAOPT.EQ.16) THEN
       CALL PRC ('  2D NEUTRAL V/A  16 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NEUT2D_VAOPT.EQ.17) THEN
       CALL PRC ('  2D NEUTRAL V/A  17 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       VIN = SQRT(2EBD/(MI.(1/SQRT($)-
     >1))),')
       CALL PRC ('                             $ IN (0,1)')
      ELSEIF (NEUT2D_VAOPT.EQ.18) THEN
       CALL PRC ('  2D NEUTRAL V/A  18 : THETA = 0.0 - NORMAL TO SURFACE
     >')
       CALL PRC ('                       Y(E) = E/(E+Ebd)**3 * [1-SQRT((
     >E+Ebd)/(G(1-G)Eimp))]')
       call prc ('                       G = 4 mi*mb/(mi+mb)**2')
       call prc ('                       E selected randomly from Y(E)')
       call prc ('                       VIN = SQRT(2 * E /mi) ')
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
      ELSEIF (NEUT2D_VAOPT.EQ.19) THEN
       CALL PRC ('  2D NEUTRAL V/A  19 : THETA = ATAN(TAN(BETA)COS(PHI))
     >')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
       CALL PRC ('                       *SQRT(|COS(B)**2+SIN(B)**2.COS(
     >PHI)**2|)')
       call prc ('                       BETA = 2PI$, $ in [0,1]')
       call prc ('                       PHI  = 2PI$ - PI, $ in [0,1]')
      ENDIF
C
      ENDIF
C
c     End of options for neutrals - not applicable for ion injections
c
      ENDIF

C-----------------------------------------------------------------------
      IF     (CSPUTOPT.EQ.1) THEN
       CALL PRC ('  SPUTTER SOURCE   1 : Formulation due to Bohdansky')
       CALL PRC ('                       (modified coefficients)')
      ELSEIF (CSPUTOPT.EQ.2) THEN
       CALL PRC ('  SPUTTER SOURCE   2 : Eckstein IPP9/82 (1993)')
      ELSEIF (CSPUTOPT.EQ.3) THEN
       CALL PRC ('  SPUTTER SOURCE   3 : Based on Eckstein IPP9/82 (1993
     >)')
       call prc ('                       Slight changes to H,D,T on C co
     >efficients')
       call prc ('                       From Garcia/Rosales-Roth 1996')
       call prc ('                       Slight changes to D on W co'//
     >'efficients')
      elseif (csputopt.eq.4) then
       CALL PRC ('  SPUTTER SOURCE   4 : Specified CONSTANT Yield Value'
     >)
       call prr ('                       Yield = ',const_yield)
      ELSEIF (CSPUTOPT.EQ.5) THEN
       CALL PRC ('  SPUTTER SOURCE   5 : Defaults to SPUTTER SOURCE'//
     >       ' 3 except for cases with custom data specified instead')
       call prc ('                       Default:'//
     >                            'Based on Eckstein IPP9/82 (1993)')
       call prc ('                       Slight changes to H,D,T on C co
     >efficients')
       call prc ('          List of custom data sets:')
       call prc ('          W Custom Yield Data : K. Krieger')
       call prc ('          D-> Be and Be-> Be  : Ecsktein 2002') 
       call prc ('          D->C and C->C       : Eckstein 2002')
       if ((cion.eq.4.or.cion.eq.6).and.crmb.eq.2.0) then 
       ! Be or C selected
          if (cion.eq.4) then 
             call prc('     BERYLLIUM SPUTTERING YIELD DATA SELECTED:')
          elseif (cion.eq.6) then 
             call prc('     CARBON SPUTTERING YIELD DATA SELECTED:')
          endif
          if (extra_sputter_angle.lt.0.0) then
             call prc('     - ANGLE AVERAGED DATA SELECTED')
          else 
             call prr('     - DATA SELECTED FOR INCIDENT ANGLE =',
     >                  extra_sputter_angle)
          endif
          call print_eck2002_yields(datunit)
       elseif (cion.eq.74) then 
          ! W selected 
             CALL PRC ('    TUNGSTEN SPUTTERING DATA SELECTED:')
       endif
      ELSEIF (CSPUTOPT.EQ.6) THEN
       CALL PRC ('  SPUTTER SOURCE   6 : Based on Eckstein'//
     >                               ' "Sputtering Yields" 2007')
       call prc ('                       Defaults to '//
     >      'Sputter data option 3 (modified  Eckstein IPP9/82 (1993))') 
       call prc ('                       if 2007 data is unavailable')

       call print_eck2007_yields(datunit)

      ENDIF
C-----------------------------------------------------------------------
      IF     (CCHEMOPT.EQ.1) THEN
       CALL PRC ('  CHEM SPUTTER SRC 1 : DIVIMP - Garcia-Rosales/Roth 19
     >94')
      ELSEIF (CCHEMOPT.EQ.2) THEN
       CALL PRC ('  CHEM SPUTTER SRC 2 : DIVIMP - Garcia-Rosales/Roth 19
     >96')
      ELSEIF (CCHEMOPT.EQ.3) THEN
       CALL PRC ('  CHEM SPUTTER SRC 3 : JET 1 - Garcia-Rosales formula(
     >EPS94)')
      ELSEIF (CCHEMOPT.EQ.4) THEN
       CALL PRC ('  CHEM SPUTTER SRC 4 : JET 2 - Pospieszczyk (EPS95)')
      ELSEIF (CCHEMOPT.EQ.5) THEN
       CALL PRC ('  CHEM SPUTTER SRC 5 : JET 3 - Vietzke (in Physical pr
     >ocesses')
       call prc ('                       of the inter-action of Fusion P
     >lasmas')
       call prc ('                       with Solids)')
      ELSEIF (CCHEMOPT.EQ.6) THEN
       CALL PRC ('  CHEM SPUTTER SRC 6 : JET 4 - Haasz (Submitted to J.N
     >ucl.Mater.,')
       call prc ('                               Dec. 1995)')
      ELSEIF (CCHEMOPT.EQ.7) THEN
       CALL PRC ('  CHEM SPUTTER SRC 7 : JET 5 - Roth & Garcia-Rosales')
       call prc ('                       (Submitted to Nucl.Fusion, Marc
     >h 1996)')
      ELSEIF (CCHEMOPT.EQ.8) THEN
       CALL PRC ('  CHEM SPUTTER SRC 8 : JET 6 - Haasz 1997 (Brian Mech
     >PhD Thesis')
       call prc ('                       Fits are based on D.')
       call prc ('                       Approximate mass dependence'//
     >                                ' included for H.') 
      ELSEIF (CCHEMOPT.EQ.9) THEN
       CALL PRC ('  CHEM SPUTTER SRC 9 : DIVIMP - FIXED YIELD')
       call prr ('                       Ychem = ',const_yield)
      ELSEIF (CCHEMOPT.EQ.10) THEN
       CALL PRC ('  CHEM SPUTTER SRC 10: DIVIMP- Haasz 1997 (Brian Mech
     >PhD Thesis')
       call prc ('                       Fits are based on D.')
       call prc ('                       Approximate mass dependence'//
     >                                ' included for H.') 
       call prc ('                       Modified to reduce the yield')
       call prc ('                       to 1/5 of its value as the')
       call prc ('                       temperature drops from 10eV to
     >5eV')
       call prc ('                       Constant at 1/5 below 5eV.')
      ELSEIF (CCHEMOPT.EQ.11) THEN
       CALL PRC ('  CHEM SPUTTER SRC 11: DIVIMP- Haasz 1997 (Brian Mech
     >PhD Thesis')
       call prc ('                       Fits are based on D.')
       call prc ('                       Approximate mass dependence'//
     >                                ' included for H.') 
       call prc ('                       Modified to reduce the yield')
       call prc ('                       to 1/5 of its value as the')
       call prc ('                       temperature drops from 10eV to
     >5eV')
       call prc ('                       Constant at 1/5 below 5eV.')
       call prc ('                       Also modified to reduce the'//
     >                                 ' effective')
       call prc ('                       surface temperature as flux'//
     >                                 ' increases.')
       call prc ('                       Flux < 1.0e19 - Teff = Tsurf')
       call prc ('                       Flux > 1.0e20 - Teff = '//
     >                                  'Tsurf-100.0')
       call prc ('                       Linear adjustment going from'//
     >                                 ' 0 to 100 as flux')
       call prc ('                       increases in that range.')
      ENDIF

C-----------------------------------------------------------------------
      IF     (CNEUTD.EQ.0) THEN
       CALL PRC ('  SPUTTER OPTION   0 : SPUTTERING BY BACKGROUND IONS O
     >NLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
      ELSEIF (CNEUTD.EQ.1) THEN
       CALL PRC ('  SPUTTER OPTION   1 : SPUTTERING BY SPECIFIED ION TYP
     >E ONLY')
       CALL PRI ('                       EIMP=TB(2+3ZIMP),  ZIMP=',
     >   CBOMBZ)
       CALL PRC ('                       EMAX=EIMP')
      ELSEIF (CNEUTD.EQ.2) THEN
       CALL PRC ('  SPUTTER OPTION   2 : INITIAL CHEMICAL SPUTTERING SOU
     >RCE ONLY')
      ELSEIF (CNEUTD.EQ.3) THEN
       CALL PRC ('  SPUTTER OPTION   3 : INITIAL SPUTTERING BY BACKGND I
     >ONS ONLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
      ELSEIF (CNEUTD.EQ.4) THEN
       CALL PRC ('  SPUTTER OPTION   4 : INITIAL SPUTTERING BY BACKGND I
     >ONS ONLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
       call prc ('                       SELF-SPUTTERED PARTICLE ENERGY'
     >)
       call prc ('                       CUT-OFFS ARE MULTIPLIED BY THE'
     >)
       call prc ('                       EMAX FACTOR BELOW')
       CALL PRR ('                       EMAX-FACTOR =',CEMAXF)
      ELSEIF (CNEUTD.EQ.5) THEN
       CALL PRC ('  SPUTTER OPTION   5 : INITIAL CHEMICAL SPUTTERING SOU
     >RCE ONLY.')
      ELSEIF (CNEUTD.EQ.6) THEN
       CALL PRC ('  SPUTTER OPTION   6 : Combined Physical and Chemical
     >Sputtering')
       call prc ('                       Two groups of atoms are launche
     >d')
       call prc ('                       The first uses PHYSICAL sputter
     >ing')
       call prc ('                       The second uses CHEMICAL sputte
     >ring')
       call prc ('                       The physically sputtered group'
     >)
       call prc ('                       uses the specified Initial Neut
     >ral')
       call prc ('                       V/A flag. The chemically sputte
     >red')
       call prc ('                       group uses a V/A flag of 3.')
       call prc ('                       The ratio of particles launched
     >')
       call prc ('                       through each method is proporti
     >onal to')
       call prc ('                       the total BG FLUX*YIELD for eac
     >h type')
       call prc ('                       of sputter source')
       call prc ('                       Total number of particles launc
     >hed')
       call pri ('                       is NIMPS + NIMPS2 = ',nimps+nim
     >ps2)
       call prb
       call prc ('                       PHYSICAL SPUTTERING: ')
       call prc ('                       INITIAL SPUTTERING BY BACKGND I
     >ONS')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
       call prb
       call prc ('                       CHEMICAL SPUTTERING: ')
       CALL PRC ('                       INITIAL CHEMICAL SPUTTERING SOU
     >RCE')
       CALL PRC ('                       THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
       call prb
       call prc ('                       The number of supplementary par
     >ticles')
       call prc ('                       specified, if any, is overwritt
     >en.')
c
c IPP/02 Geier option 7 for sputtering from external flux file      
      ELSEIF (CNEUTD.EQ.7) THEN
       CALL PRC ('  SPUTTER OPTION   7 : SPUTTERING ALONG TARGET BY ION 
     > IMPACT. FLUX TAKEN FROM EXTERNAL FILE.')
c
      ELSEIF (CNEUTD.EQ.8) THEN
       CALL PRC ('  SPUTTER OPTION   8 : PHYSICAL SPUTTERING BY WALL PLA
     >SMA CONTACT')
       call prc ('                       WALL PLASMA CONDITIONS SPECIFIE
     >D BY WALL PLASMA OPTION')
       call prc ('                       WALL ELEMENTS ACROSS THE TARGET
     >S HAVE TARGET')
       call prc ('                       CONDITIONS ASSIGNED')
      ENDIF
c
c
c
      CALL PRR ('                       SPUTTERING ENHANCEMENT FACTOR',
     >   CSEF)
      CALL PRI ('                       MAXIMUM GENERATIONS: ',CMAXGENS)
C-----------------------------------------------------------------------
c
c     Print yield modifier functions
c
      if (nymfs.gt.0) then

         call prb
         call prc ('        DESCRIPTION OF YMFS :')
         call prb
         call prc ('                   TARGET          WALLS')
         call prc ('        ELEMENTS   PHYS    CHEM    PHYS    CHEM')
      endif
c
      do in = 1,nymfs
c
         if (cymfs(in,1).eq.0.0) then
c      
            write(coment,'(1x,''DEFAULTS: '',
     >             i5,'' TO '',i5,4(2x,f7.3))')
     >          INT(cymfs(in,1)), INT(cymfs(in,2)),
     >          cymfs(in,3),cymfs(in,5),cymfs(in,6),cymfs(in,7)
       
            call prc (coment)
c      
         else 
            write(coment,'(12x,
     >             i5,'' TO '',i5,4(2x,f7.3))')
     >          INT(cymfs(in,1)), INT(cymfs(in,2)),
     >          cymfs(in,3),cymfs(in,5),cymfs(in,6),cymfs(in,7)
       
            call prc (coment)
c      
         endif
c
      end do 
c
      call prb   
C-----------------------------------------------------------------------

      IF (CEBD.EQ.0.0) THEN
       CALL PRC  ('  CHARACTERISTIC ENERGY (EBD) READ FROM LOADED DATA')
      ELSEIF (cebd.eq.-1.0) then
        CALL PRC  ('  CHARACTERISTIC ENERGY (EBD) READ FROM DATA FILE')
      else
        CALL PRR  ('  CHARACTERISTIC ENERGY        EBD   (EV)   ', CEBD)
      ENDIF

C-----------------------------------------------------------------------
c
c     Self-sputtering
c
      if (cselfs.eq.0) then
        call prc ('  SELF-SPUTTER OPT 0 : SELF-SPUTTERING IS OFF')
        call prc ('                       SELF-SPUTTERING IS CONTROLLED'
     >)
        call prc ('                       THROUGH THIS OPTION ONLY!')
        call prc ('                       COMMENTS IN OTHER SPUTTER OPTI
     >ONS')
        call prc ('                       MAY BE MISLEADING AND WILL BE'
     >)
        call prc ('                       CORRECTED AT A FUTURE DATE.')
      elseif (cselfs.eq.1) then
        call prc ('  SELF-SPUTTER OPT 1 : SELF-SPUTTERING IS ON')
        call prc ('                       SELF-SPUTTERING IS CONTROLLED'
     >)
        call prc ('                       THROUGH THIS OPTION ONLY!')
        call prc ('                       COMMENTS IN OTHER SPUTTER OPTI
     >ONS')
        call prc ('                       MAY BE MISLEADING AND WILL BE'
     >)
        call prc ('                       CORRECTED AT A FUTURE DATE.')
       CALL PRC ('                       PROPER SELF-SPUTTERING USING AC
     >TUAL ZIMP')
       CALL PRC ('                       VALUES ON EXIT AND SAME VEL/ANG
     > FLAG.')
       CALL PRC ('                       EIMP=3TB.ZIMP+5.22E-9.MI.VEXIT.
     >VEXIT+2TI')
       CALL PRR ('                       EMAX=EIMP.  THRESHOLD YIELD=',
     >   CTRESH)


       if (nymfs.gt.1) then

          call prb
          call prc ('         DESCRIPTION OF SELF-SPUTTERING YMFS :')
          call prb

       endif
c
       do in = 1,nymfs
          if (cymfs(in,4).le.-99.0) then

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' IMPURITY REFLECTION PROB = '',f5.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >           max(cymfs(in,4)+100.0,0.0)
             call prc (coment)

             write(coment,'(10x,
     >          ''- V/A FLAG 3 APPLIES TO REFLECTED IONS'',
     >          ''  ENERGY = '',f10.3)') ctem1

          elseif (cymfs(in,4).lt.0.0) then

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' THE SELF-SPUTTERING YIELD IS '',f10.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >            abs(cymfs(in,4))

          else

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' THE Y.M.F. IS '',f10.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >            cymfs(in,4)

          endif

          call prc(coment)

       end do
       call prb

      elseif (cselfs.eq.2) then
        call prc ('  SELF-SPUTTER OPT 2 : SELF-SPUTTERING IS ON')
        call prc ('                       SELF-SPUTTERING IS CONTROLLED'
     >)
        call prc ('                       THROUGH THIS OPTION ONLY!')
        call prc ('                       COMMENTS IN OTHER SPUTTER OPTI
     >ONS')
        call prc ('                       MAY BE MISLEADING AND WILL BE'
     >)
        call prc ('                       CORRECTED AT A FUTURE DATE.')
       CALL PRC ('                       SELF-SPUTTERING ENERGY IS FIXED
     >')
       CALL PRC ('                       FOR EACH SEGMENT WITH A FIXED Y
     >IELD.')
       call prr ('                       ENERGY FOR FIXED YIELD SPUTTERI
     >NG (eV) =',ctem1)
       CALL PRR ('                       THRESHOLD YIELD=',
     >   CTRESH)


       if (nymfs.gt.1) then

          call prb
          call prc ('         DESCRIPTION OF SELF-SPUTTERING YMFS :')
          call prb

       endif

       do in = 1,nymfs

          if (cymfs(in,4).le.-99.0) then

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' IMPURITY REFLECTION PROB = '',f5.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >           max(cymfs(in,4)+100.0,0.0)

             call prc (coment)

             write(coment,'(10x,
     >          ''- V/A FLAG 3 APPLIES TO REFLECTED IONS'',
     >          ''  ENERGY = '',f10.3)') ctem1

          elseif (cymfs(in,4).lt.0.0) then

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' THE SELF-SPUTTERING YIELD IS '',f10.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >            abs(cymfs(in,4))

          else

             write(coment,'(10x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >       '' THE Y.M.F. IS '',f10.3)')
     >           INT(cymfs(in,1)), INT(cymfs(in,2)),
     >            cymfs(in,4)

          endif

          call prc(coment)

       end do
       call prb


      endif
c
C-----------------------------------------------------------------------
c     Supplementary sputter source option
C-----------------------------------------------------------------------
c
      if (nimps2.gt.0) then
c
      IF     (CNEUTD2.EQ.0) THEN
       CALL PRC ('  SUP. SPUTTER OPT 0 : SPUTTERING BY BACKGROUND IONS O
     >NLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
      ELSEIF (CNEUTD2.EQ.1) THEN
       CALL PRC ('  SUP. SPUTTER OPT 1 : SPUTTERING BY SPECIFIED ION TYP
     >E ONLY')
       CALL PRI ('                       EIMP=TB(2+3ZIMP),  ZIMP=',
     >   CBOMBZ)
       CALL PRC ('                       EMAX=EIMP')
      ELSEIF (CNEUTD2.EQ.2) THEN
       CALL PRC ('  SUP. SPUTTER OPT 2 : INITIAL CHEMICAL SPUTTERING SOU
     >RCE ONLY')
      ELSEIF (CNEUTD2.EQ.3) THEN
       CALL PRC ('  SUP. SPUTTER OPT 3 : INITIAL SPUTTERING BY BACKGND I
     >ONS ONLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
      ELSEIF (CNEUTD2.EQ.4) THEN
       CALL PRC ('  SUP. SPUTTER OPT 4 : INITIAL SPUTTERING BY BACKGND I
     >ONS ONLY')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
      ELSEIF (CNEUTD2.EQ.5) THEN
       CALL PRC ('  SUP. SPUTTER OPT 5 : INITIAL CHEMICAL SPUTTERING SOU
     >RCE ONLY')
      ELSEIF (CNEUTD2.EQ.6) THEN
       CALL PRC ('  SUP. SPUTTER OPT 6 : Combined Physical and Chemical
     >Sputtering')
       call prc ('                       Two groups of atoms are launche
     >d')
       call prc ('                       The first uses PHYSICAL sputter
     >ing')
       call prc ('                       The second uses CHEMICAL sputte
     >ring')
       call prc ('                       The physically sputtered group'
     >)
       call prc ('                       uses the specified Supplementar
     >y')
       call prc ('                       V/A flag. The chemically sputte
     >red')
       call prc ('                       group uses a V/A flag of 3.')
       call prc ('                       The ratio of particles launched
     >')
       call prc ('                       The ratio of particles launched
     >')
       call prc ('                       through each method is proporti
     >onal to')
       call prc ('                       the total BG FLUX*YIELD for eac
     >h type')
       call prc ('                       of sputter source')
       call prc ('                       Total number of particles launc
     >hed')
       call pri ('                       is NIMPS + NIMPS2 = ',nimps+nim
     >ps2)
       call prb
       call prc ('                       PHYSICAL SPUTTERING: ')
       call prc ('                       INITIAL SPUTTERING BY BACKGND I
     >ONS')
       CALL PRC ('                       EIMP=TB(2+3ZB)')
       CALL PRC ('                       EMAX=EIMP.GAMMA(1-GAMMA)-EBD')
       call prb
       call prc ('                       CHEMICAL SPUTTERING: ')
       CALL PRC ('                       INITIAL CHEMICAL SPUTTERING SOU
     >RCE')
       CALL PRC ('                       THETA =+/-ASIN(SQRT($)), $ IN (
     >0,1)')
       CALL PRR ('                       VIN = SQRT(2EIN/MI), EIN=',
     >    CTEM1)
       call prb
       call prc ('                       The number of supplementary par
     >ticles')
       call prc ('                       specified, if any, is overwritt
     >en.')
c
c IPP/02 Geier option 7 for sputtering from external flux file      
      ELSEIF (CNEUTD2.EQ.7) THEN
       CALL PRC ('  SPUTTER OPTION   7 : SPUTTERING ALONG TARGET BY ION 
     > IMPACT. FLUX TAKEN FROM EXTERNAL FILE.')
c
      ELSEIF (CNEUTD2.EQ.8) THEN
       CALL PRC ('  SUP. SPUTTER OPT 8 : PHYSICAL SPUTTERING BY WALL PLA
     >SMA CONTACT')
       call prc ('                       WALL PLASMA CONDITIONS SPECIFIE
     >D BY WALL PLASMA OPTION')
       call prc ('                       WALL ELEMENTS ACROSS THE TARGET
     >S HAVE TARGET')
       call prc ('                       CONDITIONS ASSIGNED')
      ENDIF
c
      endif
c
C-----------------------------------------------------------------------
C
c     Print out wall plasma options for a wall launch using the wall
c     plasma sputtering options.
c
      if ((cneutb.eq.2.or.cneutb.eq.4).and. 
     >    (cneutd.eq.8.or.cneutd2.eq.8)) then 
c
c      Print WALL_PLASMA_OPT options  
c
       if (wall_plasma_opt.eq.0) then   
c
        call prc('  WALL PLASMA OPT 0  : PLASMA CONDITIONS AT WALLS ARE'
     >   //' TAKEN FROM CLOSEST PLASMA CELL')
        call prc('                       WALL ELEMENTS AT THE TARGETS'
     >   //' ARE ASSIGNED ACTUAL TARGET CONDITIONS') 
c
       elseif (wall_plasma_opt.eq.1) then   
c
        call prc('  WALL PLASMA OPT 1  : PLASMA CONDITIONS AT WALLS ARE'
     >   //' CALCULATED AS FOLLOWS:')
        call prc('                       -FIND NEAREST PLASMA CELL FOR'
     >   //' BASE CONDITIONS')
        call prc('                       -FIND S-CORRDINATE ALONG OUTER'
     >   //' RING THAT CORRESPONDS TO WALL ELEMENT LOCATION')
        call prc('                       -INTERPOLATE ALONG RING TO FIN'
     >   //'D PLASMA CONDITIONS FOR WALL ELEMENT')
        call prc('                       -TAKE CLOSEST DISTANCE TO'
     >   //' OUTER RING FROM WALL ELEMENT')    
        call prc('                       -APPLY A LINEAR DECAY OVER'
     >   //' THIS DISTANCE WITH')
        call prr('                        SCALE LENGTH =',
     >             wall_plasma_fact )
        call prc('                       -MINIMUM VALUES ARE'//
     >           ' HARD CODED')
c
       elseif (wall_plasma_opt.eq.2) then   
c
        call prc('  WALL PLASMA OPT 2  : PLASMA CONDITIONS AT WALLS ARE'
     >   //' CALCULATED AS FOLLOWS:')
        call prc('                       -FIND NEAREST PLASMA CELL FOR'
     >   //' BASE CONDITIONS')
        call prc('                       -FIND S-CORRDINATE ALONG OUTER'
     >   //' RING THAT CORRESPONDS TO WALL ELEMENT LOCATION')
        call prc('                       -INTERPOLATE ALONG RING TO FIN'
     >   //'D PLASMA CONDITIONS FOR WALL ELEMENT')
        call prc('                       -TAKE CLOSEST DISTANCE TO'
     >   //' OUTER RING FROM WALL ELEMENT')    
        call prc('                       -APPLY AN EXPONENTIAL DECAY '
     >   //' WITH')
        call prr('                        SCALE LENGTH =',
     >             wall_plasma_fact )
        call prc('                       -MINIMUM VALUES ARE'//
     >           ' HARD CODED')
c
       endif
c
      endif
c
C
C-----------------------------------------------------------------------
C
C    Modify normal option printing for a free-space launch
C
 998  IF     (CNEUTE.EQ.0) THEN
       CALL PRC ('  NORMAL OPTION    0 : MEASURE THETA FROM SURFACE NORM
     >AL')
       if (cneutb.eq.1.or.cneutb.eq.6.or.cneutb.eq.7
     >                .or.cneuth.eq.1) then
       call prc ('                       "SURFACE" NORMAL IS THE +ve R A
     >XIS')
         call prc ('                       FOR A FREE SPACE LAUNCH')
       endif
      ELSEIF (CNEUTE.EQ.1.OR.CNEUTE.EQ.2) THEN
       WRITE (COMENT,'(A,I1,A,F7.2,A)') '  NORMAL OPTION    ',CNEUTE,
     >   ' : MEASURE THETA FROM ',CSNORM*RADDEG,' DEGS TO X=0'
       CALL PRC (COMENT)
       IF (CNEUTE.EQ.1) WRITE (7,'(23X,A)') 'FOR PRIMARIES ONLY'
       IF (CNEUTE.EQ.2) WRITE (7,'(23X,A)') 'FOR PRIMARIES & SELF-SPUT.'
      ENDIF
C-----------------------------------------------------------------------
      IF     (CNEUTF.EQ.0) THEN
       CALL PRC ('  NEUT SPREADING   0 : OFF (LAUNCH AT MESH PTS ONLY)')
      ELSEIF (CNEUTF.EQ.1) THEN
C      CALL PRC ('  NEUT SPREADING   1 : ON  (LAUNCH BETWEEN MESH PTS)')
       CALL PRC ('  NEUT SPREADING   1 : NOT IMPLEMENTED')
      ENDIF
C-----------------------------------------------------------------------
      IF     (CNEUTG.EQ.0) THEN
       CALL PRC ('  INITIAL ION VEL  0 : VI = 0.0')
      elseIF     (CNEUTG.EQ.1) THEN
       CALL PRC ('  INITIAL ION VEL  1 : +/-0.5VN ALONG S')
      ELSEIF (CNEUTG.EQ.2) THEN
       CALL PRC ('  INITIAL ION VEL  2 : VN ALONG S AWAY FROM TARGET')
      ELSEIF (CNEUTG.EQ.3) THEN
       CALL PRC ('  INITIAL ION VEL  3 : +/-SQRT($).VN ALONG S, $ IN (0,
     >1)')
      ENDIF
c
C-----------------------------------------------------------------------
C  
c     Initial position option
C
      if (init_pos_opt.eq.0) then 
       CALL PRC ('  I24-INIT POSITION 0: IONS '//
     >                   'RECYCLE/SPUTTER AS NEUTRALS FROM THE CENTER')
       CALL PRC ('                       OF  THE TARGET'//
     >                   ' ELEMENT OF IMPACT')
       CALL PRC ('                       THE INITIAL S,CROSS'//
     >                   ' COORDINATES FOR AN ION')
       CALL PRC ('                       ARE AT THE CENTER OF THE'//
     >                             ' CELL WHERE THEY ARE CREATED')
      elseif (init_pos_opt.eq.1) then 
       CALL PRC ('  I24-INIT POSITION 1: IONS'//
     >               'RECYCLE/SPUTTER AS NEUTRALS FROM THE POSITION')
       CALL PRC ('                       OF IMPACT ON THE TARGET'//
     >                   ' ELEMENT (BASED ON CROSS COORDINATE) ')
       CALL PRC ('                       THE INITIAL S,CROSS'//
     >                   ' COORDINATES FOR AN ION')
       CALL PRC ('                       ARE APPROXIMATED BASED'//
     >                             ' ON THE ACTUAL R,Z LOCATION')
       CALL PRC ('                       OF CREATION')
      elseif (init_pos_opt.eq.2) then 
c
c      NOTE - NOT YET IMPLEMENTED
c
       CALL PRC ('  I24-INIT POSITION 2: IONS'//
     >               'RECYCLE/SPUTTER AS NEUTRALS FROM THE POSITION')
       CALL PRC ('                       OF IMPACT ON THE TARGET'//
     >                   ' ELEMENT (BASED ON CROSS COORDINATE) ')
       CALL PRC ('                       THE INITIAL S,CROSS'//
     >                   ' COORDINATES FOR AN ION')
       CALL PRC ('                       ARE APPROXIMATED BASED'//
     >                             ' ON THE ACTUAL R,Z LOCATION')
       CALL PRC ('                       OF CREATION')
       call prc ('                       PRIMARY NEUTRALS ARE ALSO'//
     >                   ' DISTRIBUTED ACROSS THE ELEMENT OF ORIGIN')
      endif

c
C-----------------------------------------------------------------------
c
      IF (cneutvel.EQ.0) THEN
       CALL PRC ('  IMP.NEUT VEL OPT 0 : IMPURITY NEUTRAL VELOCITY OPTIO
     >N 0')
       call prc ('                       IMPURITY NEUTRAL WILL HAVE A CO
     >NSTANT VELOCITY')
       call prc ('                       FROM CREATION TO REMOVAL.')
      ELSEIF (cneutvel.eq.1) THEN
       CALL PRC ('  IMP.NEUT VEL OPT 1 : IMPURITY NEUTRAL VELOCITY OPTIO
     >N 1')
       call prc ('                       IMPURITY NEUTRAL WILL HAVE A VE
     >LOCITY')
       call prc ('                       BASED ON THE LOCAL PLASMA ION')
       call prc ('                       TEMPERATURE. THIS ASSIGNED VELO
     >CITY WILL CHANGE')
       call prc ('                       AS THE NEUTRAL MOVES INTO REGIO
     >NS ')
       call prc ('                       WITH DIFFERENT PLASMA ION TEMPE
     >RATURES.')
      ELSEIF (cneutvel.eq.2) THEN
       CALL PRC ('  IMP.NEUT VEL OPT 2 : IMPURITY NEUTRAL VELOCITY OPTIO
     >N 2')
       call prc ('                       IMPURITY NEUTRAL WILL HAVE A VE
     >LOCITY')
       call prc ('                       THAT REMAINS CONSTANT BETWEEN M
     >TC EVENTS.')
       call prc ('                       WHEN AN MTC EVENT OCCURS THE VE
     >LOCITY')
       call prc ('                       OF THE IMPURITY NEUTRAL WILL BE
     > SET ')
       call prc ('                       TO A VALUE BASED ON THE LOCAL P
     >ASMA ION')
       call prc ('                       TEMPERATURE AND THE MASS OF THE
     > IMPURITY.')
      ENDIF
C-----------------------------------------------------------------------
      IF (NRFOPT.EQ.0) THEN
       CALL PRC ('  REFLECTION OPT   0 : NEUTRAL REFLECTION - OFF')
      ELSEIF (NRFOPT.EQ.1) THEN
       CALL PRC ('  REFLECTION OPT   1 : NEUTRAL REFLECTION - ON')
       CALL PRC ('                       ATOMS ARE REFLECTED SPECULARLY'
     >)
       CALL PRC ('                       AT WALL IMPACT')
      ELSEIF (NRFOPT.EQ.2) THEN
       CALL PRC ('  REFLECTION OPT   2 : NEUTRAL REFLECTION - ON')
       CALL PRC ('                       ATOMS ARE REFLECTED WITH A COSI
     >NE')
       CALL PRC ('                       DISTRIBUTION AT WALL IMPACT')
      ENDIF
c
      if (nrfopt.eq.1.or.nrfopt.eq.2) then
c
c     Print out individual segment reflection coefficients if active
c
c
       ind = 0
c
       do in = 1,nymfs
c
          if (ind.eq.0) then

             call prb
             call prc ('         SEGMENT WALL REFLECTION IS ACTIVE:')
             call prb
             call prc ('     THE FOLLOWING REFLECTION COEFFICIENTS'//
     >                    ' WILL BE USED')
             call prc ('     - LISTING FOR SEGMENT 0 IS THE'//
     >                    ' DEFAULT VALUE')
             ind = 1
          endif
c
          if (cymfs(in,8).gt.0.0) then
c
             write(coment,'(5x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >                      '' NEUTRALS REFELCTED WITH WEIGHT:'',
     >                         f10.3)')
     >                         INT(cymfs(in,1)), INT(cymfs(in,2)),
     >                         cymfs(in,8)
c
             call prc(coment)
c
          elseif (cymfs(in,8).lt.0.0) then
c
             write(coment,'(5x,''FOR ELEMENTS: '',i5,'' TO '',i5,
     >                      '' NEUTRALS PTR WITH WEIGHT:'',f10.3)')
     >                         INT(cymfs(in,1)), INT(cymfs(in,2)),
     >                         abs(cymfs(in,8))
             call prc(coment)
c
             write (coment,'(5x,'' - ENERGY FOR PTR IS : '',f10.3)')
     >                        ctem1
c
             call prc(coment)
c

          endif
c
       end do
c
      endif
c
C-----------------------------------------------------------------------
      IF (cfolrec.eq.0) then
       CALL PRC ('  FOLLOW REC. OPT  0 : OFF - IMPURITY IONS ARE NOT FOL
     >LOWED')
       call prc ('                       AFTER RECOMBINING TO FORM NEUTR
     >ALS')
      ELSEIF (cfolrec.EQ.1) THEN
       CALL PRC ('  FOLLOW REC. OPT  1 : ON - IMPURITY IONS THAT RECOMBI
     >NE')
       call prc ('                       TO FORM NEUTRALS ARE FOLLOWED A
     >S A')
       call prc ('                       NEUTRAL - THESE NEUTRALS START
     >FROM')
       call prc ('                       THE POSITION OF RECOMBINATION W
     >ITH')
       call prc ('                       AN ISOTROPIC 3-D VELOCITY DISTR
     >IBUTION')
       call prc ('                       THE MAGNITUDE OF THIS VELOCITY
     >IS')
       call prc ('                       CALCULATED FROM THE ACTUAL ION
     >TEMPERATURE')
       call prc ('                       RECOMBINED NEUTRAL VELOCITY IS
     >MULTIPLIED')
       call prr ('                       BY: ',cvrmult)
      ENDIF
C-----------------------------------------------------------------------
      IF (mtcopt.eq.0) then
       CALL PRC ('  NEUT.MOM.COL OPT 0 : OFF - IMPURITY NEUTRAL MOMENTUM
     >')
       call prc ('                       TRANSFER COLLISON OPTION IS NOT
     >')
       CALL PRC ('                       ACTIVE.')
      ELSEIF (MTCOPT.EQ.1) THEN
       CALL PRC ('  NEUT.MOM.COL OPT 1 : ON - IMPURITY NEUTRAL MOMENTUM'
     >)
       call prc ('                       TRANSFER COLLISON OPTION IS ACT
     >IVATED.')
       call prc ('                       IMPURITY NEUTRALS WILL UNDERGO
     >A')
       call prc ('                       +/-90 DEGREE PATH CHANGE WITH A
     > FREQUENCY')
       CALL PRC ('                       CALCULATED FROM THE FOLLOWING')
       CALL PRC ('                       FORMULA:')
       CALL PRC ('                       Vfr = 16mb/3(mz+mb) * ')
       call prc ('                             (kelighi * nb+  + kelighg
     >* nb0)')
       call prr ('                       Where KELIGHI = ',kelighi)
       call prr ('                         and KELIGHG = ',kelighg)
       call prc ('                       Nb0 = local background neutral
     >density')
       call prc ('                       Nb+ = local background ion dens
     >ity')
c
      ELSEIF (MTCOPT.EQ.2) THEN
       CALL PRC ('  NEUT.MOM.COL OPT 2 : ON - IMPURITY NEUTRAL MOMENTUM'
     >)
       call prc ('                       TRANSFER COLLISON OPTION IS ACT
     >IVATED.')
       call prc ('                       IMPURITY NEUTRALS WILL UNDERGO
     >A')
       call prc ('                       3D 90 DEGREE PATH CHANGE IF A
     >3D VELOCITY IS SPECIFIED')  
       call prc ('                       WITH A FREQUENCY')
       CALL PRC ('                       CALCULATED FROM THE FOLLOWING')
       CALL PRC ('                       FORMULA:')
       CALL PRC ('                       Vfr = 16mb/3(mz+mb) * ')
       call prc ('                             (kelighi * nb+  + kelighg
     >* nb0)')
       call prr ('                       Where KELIGHI = ',kelighi)
       call prr ('                         and KELIGHG = ',kelighg)
       call prc ('                       Nb0 = local background neutral
     >density')
       call prc ('                       Nb+ = local background ion dens
     >ity')
c
      endif
C-----------------------------------------------------------------------
      IF (prompt_depopt.EQ.0) THEN
       CALL PRC ('  PROMPT DEP OPT   0 : PROMPT ION REDEPOSITION OPTION
     >- OFF')
      ELSEIF (prompt_depopt.EQ.1) THEN
       CALL PRC ('  PROMPT DEP OPT   1 : PROMPT ION REDEPOSITION OPTION
     >- ON')
       CALL PRC ('                       INITIAL IONS WITHIN ONE LARMOR
     >RADIUS')
       CALL PRC ('                       OF THE TARGET ARE REMOVED. IF')
       CALL PRC ('                       SELF-SPUTTERING IS ON THESE ION
     >S')
       CALL PRC ('                       MAY CAUSE SPUTTERING EVENTS')
      ELSEIF (prompt_depopt.EQ.2) THEN
       CALL PRC ('  PROMPT DEP OPT   2 : PROMPT ION REDEPOSITION OPTION
     >- ON')
       CALL PRC ('                       ALL IONS WITHIN ONE LARMOR'//
     >' RADIUS')
       CALL PRC ('                       OF THE TARGET ARE REMOVED. IF')
       CALL PRC ('                       SELF-SPUTTERING IS ON THESE ION
     >S')
       CALL PRC ('                       MAY CAUSE SPUTTERING EVENTS')
      ENDIF

c
C-----------------------------------------------------------------------
c
      IF (FPOPT.EQ.0) THEN
       CALL PRC ('  PERIPHERY OPTION 0 : HARD WALL - IONS LOST AT WALL')
      ELSEIF (FPOPT.EQ.1) THEN
       CALL PRC ('  PERIPHERY OPTION 1 : REFLECTING WALL - IONS RETURN T
     >O PLASMA')
      ELSEIF (FPOPT.EQ.2) THEN
       CALL PRC ('  PERIPHERY OPTION 2 : NO WALL - IONS DIFFUSE CROSS-FI
     >ELD')
       CALL PRC ('                       INDEFINITELY - ASSOCIATED WITH
     >OUTERMOST')
       CALL PRC ('                       RING')
      ELSEIF (FPOPT.EQ.3) THEN
       CALL PRC ('  PERIPHERY OPTION 3 : FAR PERIPHERY MODEL AT WALL')
       CALL PRC ('                       IONS TRACKED BEYOND LAST RING')
       CALL PRC ('                       USING FAR PERIPHERY MODEL')
       CALL PRR ('                       DISTANCE TO WALLS   ('
     >                                   //OUTER//'): ',FPXMAXO)
       CALL PRR ('                       DISTANCE TO WALLS   ('
     >                                   //INNER//'): ',FPXMAXI)
       CALL PRR ('                       FP TARGET LOSS TIME ('
     >                                   //OUTER//'): ',FPTIMO)
       CALL PRR ('                       FP TARGET LOSS TIME ('
     >                                   //INNER//'): ',FPTIMI)
       CALL PRR ('                       FP DIFFUSION COEFF  : ',CDPERPF
     >P)
       if (fpropt.eq.0.or.(fpropt.eq.1.and.cselfs.eq.0)) then
       call prc ('  FP RECYCLE OPTION 0: IONS ARE LOST AT FP TARGET OR W
     >ALL IMPACT')
       elseif (fpropt.eq.1) then
       call prc ('  FP RECYCLE OPTION 1: FP TARGET AND WALL LOSSES ARE R
     >EINJECTED')
       call prc ('                       FROM THE EDGE OF THE TARGET')
       endif
      ELSEIF (FPOPT.EQ.4) THEN
       CALL PRC ('  PERIPHERY OPTION 4 : MAIN VESSEL WALL - TREATED AS H
     >ARD')
       call prc ('                                        - IONS LOST AT
     > WALL')
       call prc ('                                        - AS OPTION 0'
     >)
       CALL PRC ('                       TRAP PLASMA WALL - TREATED AS R
     >EFLECTING')
       call prc ('                                        - IONS REFLECT
     >ED AT WALL')
       call prc ('                                        - AS OPTION 1'
     >)
      ELSEIF (FPOPT.EQ.5) THEN
       CALL PRC (sa//'PERIPHERY OPTION 5 :'// 
     >               ' FAR PERIPHERY MODEL AT WALL')
       CALL PRC (sp//'IONS TRACKED BEYOND LAST RING')
       CALL PRC (sp//'USING FAR PERIPHERY MODEL')
       CALL PRc (sp//'DISTANCE TO WALLS IS CALCULATED FOR EACH CELL'
     >             //' ON RING')
       CALL PRc (sp//'PARRALLEL TRANSPORT IS ON')
       call prc (sp//'PERIPHERAL PLASMA IS SPECIFIED BY THE'//
     >               ' FP_PLASMA OPTION')
       call prc (sp//'PERIPHERAL GEOMETRY TAKEN FROM LAST REAL RING')
       CALL PRR (sp//'FP DIFFUSION COEFF  : ',CDPERPFP)
c
       if (fpropt.eq.0.or.(fpropt.eq.1.and.cselfs.eq.0)) then
       call prc ('  FP RECYCLE OPTION 0: IONS ARE LOST AT FP TARGET OR W
     >ALL IMPACT')
       elseif (fpropt.eq.1) then
       call prc ('  FP RECYCLE OPTION 1: FP TARGET AND WALL LOSSES ARE R
     >EINJECTED')
       call prc ('                       FROM THE EDGE OF THE TARGET')
       endif

      ENDIF
C-----------------------------------------------------------------------
      call prb
      if (fp_neut_opt.eq.0) then 
       call prc('  FP NEUTRAL IONIZATION OPT  0: OFF ') 
      elseif (fp_neut_opt.eq.1) then 
       call prc('  FP NEUTRAL IONIZATION OPT  1: ON ') 
      endif
c
c     FP Plasma option ---------------------------------------------
c
      if (fpopt.eq.5.or.fp_neut_opt.ne.0) then 
c
       call prb
       if (fp_plasma_opt.eq.0) then  
c
        call prc('  FP PLASMA OPTION 0 : PLASMA DATA FROM NEAREST'//
     >          ' GRID CELL') 
c
       elseif (fp_plasma_opt.eq.1) then 
c
        call prc('  FP PLASMA OPTION 1 : FP DENSITY FROM NEAREST'//
     >          ' GRID CELL') 
        call prr('                       FP TEMPERATURE = ',fp_te)
c
       elseif (fp_plasma_opt.eq.2) then 
c
        call prc('  FP PLASMA OPTION 1 : FP PLASMA SPECIFIED'//
     >          ' GRID CELL') 
        call prr('                       FP TEMPERATURE = ',fp_te)
        call prr('                       FP DENSITY     = ',fp_ne)
c
       elseif (fp_plasma_opt.eq.3) then 
c
        call prc('  FP PLASMA OPTION 1 : FP PLASMA SPECIFIED'//
     >          ' GRID CELL') 
        call prr('                       FP TEMPERATURE = ',fp_te)
        call prr('                       FP DENSITY     = ',fp_ne)
        call prc('                       Vb and E       = 0.0')
c
       endif

       call prb

       if (fp_flow_opt.eq.0) then
          call prc('   FP FLOW OPTION 0: POLOIDAL DRIFT'//
     >             ' FLOW IS OFF IN PERIPHERY')
       elseif (fp_flow_opt.eq.1) then
          call prc('   FP FLOW OPTION 1: POLOIDAL DRIFT'//
     >             ' FLOW MATCHES ASSOCIATED GRID RING')
          do in = 1,num_fp_regions
             write(coment,'(a,i6,a,g12.5)') '     FLOW VEL.'//
     >                ' IN FP_REG: ',in,' = ',fp_flow_velocity(in)
             call prc(coment)
          end do

       elseif (fp_flow_opt.eq.2) then
          call prc('   FP FLOW OPTION 2: POLOIDAL DRIFT'//
     >             ' FLOW IS SPECIFIED = ',fp_flow_velocity_input)
       elseif (fp_flow_opt.eq.3) then
          call prc('   FP FLOW OPTION 3: POLOIDAL DRIFT'//
     >             ' FLOW MATCHES ASSOCIATED BOUNDARY RING')
          do in = 1,num_fp_regions
             write(coment,'(a,i6,a,g12.5)') '     FLOW VEL.'//
     >                ' IN FP_REG: ',in,' = ',fp_flow_velocity(in)
             call prc(coment)
          end do

       endif

      endif

C-----------------------------------------------------------------------
      IF (CmirOPT.EQ.0) THEN
       CALL PRC ('  TARGET MIRROR OPT 0 : TARGET IS STANDARD. SPUTTERING
     > OCCURS')
       CALL PRC ('                       NORMALLY.')
       call prc ('                       VIRTUAL POINTS ARE DISCARDED')
      ELSEIF (CmirOPT.EQ.1) THEN
       CALL PRC ('  TARGET MIRROR OPT 1 : TARGET IS A MIRROR. IONS ARE R
     >EFLECTED.')
       call prc ('                       SIGN OF VELOCITY IS CHANGED.')
       CALL PRC ('                       PARTICLES ARE ELIMINATED AFTER
     >REACHING')
       call prc ('                       THE TIME LIMIT.')
      ELSEIF (CmirOPT.EQ.2) THEN
       CALL PRC ('  TARGET MIRROR OPT 2 : TARGET IS A MIRROR. IONS ARE R
     >EFLECTED.')
       call prc ('                       SIGN OF VELOCITY IS UNCHANGED.'
     >)
       CALL PRC ('                       PARTICLES ARE ELIMINATED AFTER
     >REACHING')
       call prc ('                       THE TIME LIMIT.')
      ELSEIF (CmirOPT.EQ.3) THEN
       CALL PRC ('  TARGET MIRROR OPT 3 : '//INNER//
     >              ' TARGET IS A MIRROR. IONS ARE REFLECTED.')
       call prc ('                       SIGN OF VELOCITY IS CHANGED.')
       call prc (sp//INNER//' TARGET DENSITY SET TO ZERO TO'//
     >                      ' FORCE PARTICLE FLUXES TO ZERO')
      ELSEIF (CmirOPT.EQ.4) THEN
       CALL PRC ('  TARGET MIRROR OPT 4 : '//OUTER//
     >              ' TARGET IS A MIRROR. IONS ARE REFLECTED.')
       call prc ('                       SIGN OF VELOCITY IS CHANGED.')
       call prc (sp//OUTER//' TARGET DENSITY SET TO ZERO TO'//
     >                      ' FORCE PARTICLE FLUXES TO ZERO')
      endif


      RETURN
      END
c
c
c
      SUBROUTINE PR_transport_options
c                (NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
c
c      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
C     THIS ROUTINE PRINTS ALL THE OPTION FLAGS, TORUS ATTRIBUTES ETC
C     STORED IN COMMONS COMTOR AND COMTAU
C
C     C.M.FARRELL    NOVEMBER 1987
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
c      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
c      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
c      include    'slcom'  
      include   'printopt' 
      include   'driftvel'
c slmod begin - new
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind,irstart,irend
c      logical  prsol21,prsol22
c      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
c      real     totmhpin
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('--- IMPURITY SIMULATION TRANSPORT OPTIONS ---- ')
      CALL PRChtml ('--- IMPURITY SIMULATION TRANSPORT OPTIONS ---',
     >               'pr_transport_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------
      IF (CIOPTB.EQ.1) THEN
       CALL PRC ('  FIRST DIFFUSION  * : IGNORED')
      ELSEIF (CDIFOP.EQ.0) THEN
       CALL PRC ('  FIRST DIFFUSION  0 : DIFFUSION STARTS IMMEDIATELY')
      ELSEIF (CDIFOP.EQ.1) THEN
       CALL PRC ('  FIRST DIFFUSION  1 : AFTER T=-TAUPARA.LOG$/2, BASED
     >ON TI(0)')
      ELSEIF (CDIFOP.EQ.2) THEN
       CALL PRC ('  FIRST DIFFUSION  2 : AFTER T=TAUPARA, ITSELF CHANGIN
     >G WITH T')
      ENDIF
C-----------------------------------------------------------------------
      IF     (CIOPTJ.EQ.0) THEN
       CALL PRC ('  DPERP OPTION     0 : CONSTANT')
      ELSEIF (CIOPTJ.EQ.1) THEN
       CALL PRC ('  DPERP OPTION     1 : DPERP = DPERP0.NB0/NB IN SOL, T
     >RAP')
       CALL PRC ('                               CONSTANT DPERP0 IN MAIN
     >')
      elseIF (CIOPTJ.EQ.2) THEN
       CALL PRC ('  DPERP OPTION     2 : DPERP IS CONSTANT ALONG THE REF
     >ERENCE')
       CALL PRC ('                       LINE LOCATED AT NKS(IRSEP)/2 +1
     >')
       CALL PRC ('                       TRANSPORT ELSEWHERE IS BASED ON
     >')
       CALL PRC ('                       MOVING PARTICLES IN PROPORTION
     >TO')
       CALL PRC ('                       THEIR EQUIVALENT POSITION ON TH
     >E')
       CALL PRC ('                       REFERENCE. TRANSPORT IN THE PRI
     >VATE PLASMA')
       CALL PRC ('                       IS MAPPED RELATIVE TO THE ADJAC
     >ENT')
       CALL PRC ('                       CELL ON THE SEPARATRIX AT THE I
     >K=1 INDEX')
       CALL PRC ('                       WHICH IS IN TURN MAPPED TO THE
     >REFERENCE')
      ENDIF
C-----------------------------------------------------------------------
      IF     (cdiffopt.EQ.0) THEN
       CALL PRC ('  CROSS-FIELD STEP 0 : CONSTANT - The Probability of m
     >aking')
       call prc ('                       inward and outward diffusive st
     >eps')
       call prc ('                       equals 0.5')
      ELSEIF (Cdiffopt.EQ.1) THEN
       CALL PRC ('  CROSS-FIELD STEP 1 : The probability of making inwar
     >d and')
       call prc ('                       outward diffusive steps is equa
     >l')
       call prc ('                       to the ratio of the length of t
     >he sides')
       call prc ('                       parallel to the field lines of
     >a small')
       call prc ('                       cell located at the current cro
     >ss-field')
       call prc ('                       position of the particle.')
       call prc ('                       This applies ONLY in the core')
      ELSEIF (Cdiffopt.EQ.2) THEN
       CALL PRC ('  CROSS-FIELD STEP 2 : The probability of making inwar
     >d and')
       call prc ('                       outward diffusive steps is equa
     >l')
       call prc ('                       to the ratio of the length of t
     >he sides')
       call prc ('                       parallel to the field lines of
     >a small')
       call prc ('                       cell located at the current cro
     >ss-field')
       call prc ('                       position of the particle.')
       call prc ('                       This applies EVERYWHERE')
      ELSEIF (Cdiffopt.EQ.3) THEN
       CALL PRC ('  CROSS-FIELD STEP 3 : The probability of making inwar
     >d and')
       call prc ('                       outward diffusive steps is equa
     >l')
       call prc ('                       to the ratio of the length of t
     >he sides')
       call prc ('                       parallel to the field lines of
     >a small')
       call prc ('                       cell located at the current cro
     >ss-field')
       call prc ('                       position of the particle.')
       call prc ('                       The ratios are calculated for')
       call prc ('                       each 1/2 cell independently')
       call prc ('                       This option can ONLY be used wi
     >th')
       call prc ('                       CROSS FIELD DISTANCE OPTION 1')
       call prc ('                       This option applies EVERYWHERE'
     >)
      ENDIF
C-----------------------------------------------------------------------
      IF     (pinchopt.EQ.0) THEN
       CALL PRC ('  PINCH VEL. OPT   0 : OFF - NO PINCH VELOCITY APPLIED
     >')
      ELSEIF (pinchopt.EQ.1) THEN
       CALL PRC ('  PINCH VEL. OPT   1 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
      ELSEIF (pinchopt.EQ.2) THEN
       CALL PRC ('  PINCH VEL. OPT   2 : ON - PINCH VELOCITY IS APPLIED
     >ONLY IN MAIN SOL')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
      elseif (pinchopt.eq.3) then
       CALL PRC ('  PINCH VEL. OPT   3 : ON - PINCH VELOCITY IS APPLIED
     >ONLY IN CORE')
       call prr ('                       PERPENDICULAR PINCH AT SEP',
     >                                   cvpinch)
       call prc ('                       PINCH IS ADJUSTED BY POLOIDAL 
     >LENGTH FACTOR')

      elseif (pinchopt.eq.4) then
       CALL PRC ('  PINCH VEL. OPT   4 : ON - RADIAL VELOCITY IS DETERMI
     >NED')
       call prc (sp//'FROM A SPECIFIED PROBABILITY DISTRIBUTION FUNCTION
     >')
       if (pinch_loc_opt.eq.0) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED EVERYWHERE EXCEPT'//
     >                                ' PFZ')
       elseif (pinch_loc_opt.eq.1) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED ONLY IN MAIN SOL')
       elseif (pinch_loc_opt.eq.2) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED EVERYWHERE')
       elseif (pinch_loc_opt.eq.3) then 
         CALL PRC ('     - PDF BASED PINCH APPLIED ONLY IN'//
     >               ' MAIN SOL ABOVE X-POINT')
       elseif (pinch_loc_opt.eq.4) then 
         CALL PRC ('     - PDF BASED PINCH APPLIED IN MAIN SOL'// 
     >               ' ABOVE X-POINT AND IN CORE REGION')
       endif
       call prc (sp//'- DPERP TRANSPORT TURNED OFF IN'//
     >                      ' PINCH REGION')
c
         call prb
         call prc ('     ** PDF FUNCTION USED **')
         call prc ('      VELOCITY (M/S)     PROBABILITY')
c
         do in = 1,pinch_npdf
            write(coment,'(4x,f13.2,2x,1p,e12.5)') 
     >          pinch_pdf(in,1),pinch_pdf(in,2)
            call prc(coment)
         end do
c
      elseif (pinchopt.eq.5) then
       CALL PRC ('  PINCH VEL. OPT   5 : ON - RADIAL VELOCITY IS DETERMI
     >NED')
       call prc (sp//'FROM A SPECIFIED PROBABILITY DISTRIBUTION FUNCTION
     >')
       if (pinch_loc_opt.eq.0) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED EVERYWHERE EXCEPT'//
     >                                ' PFZ')
       elseif (pinch_loc_opt.eq.1) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED ONLY IN MAIN SOL')
       elseif (pinch_loc_opt.eq.2) then 
         CALL PRC (sp//'- PDF BASED PINCH APPLIED EVERYWHERE')
       elseif (pinch_loc_opt.eq.3) then 
         CALL PRC ('     - PDF BASED PINCH APPLIED ONLY IN'//
     >               ' MAIN SOL ABOVE X-POINT')
       elseif (pinch_loc_opt.eq.4) then 
         CALL PRC ('     - PDF BASED PINCH APPLIED IN MAIN SOL'// 
     >               ' ABOVE X-POINT AND IN CORE REGION')
       endif
       call prc (sp//'- DPERP DIFFUSIVE TRANSPORT ALSO ACTIVE IN'
     >           //' PINCH REGION')
c
         call prb
         call prc ('     ** PDF FUNCTION USED **')
         call prc ('      VELOCITY (M/S)     PROBABILITY')
c
         do in = 1,pinch_npdf
            write(coment,'(4x,f13.2,2x,1p,e12.5)') 
     >          pinch_pdf(in,1),pinch_pdf(in,2)
            call prc(coment)
         end do
c
      ELSEIF (pinchopt.EQ.6) THEN
       CALL PRC ('  PINCH VEL. OPT   6 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH DIRECTION IS DEPENDENT'//
     >           ' ON CELL LOCATION')
       call prc ('                       CELLS WITH R<Rxpt HAVE '//
     >           ' OPPOSITE SIGN OF PINCH IN R>Rxpt CELLS')
      ELSEIF (pinchopt.EQ.7) THEN
       CALL PRC ('  PINCH VEL. OPT   7 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH DIRECTION IS DEPENDENT'//
     >           ' ON CELL LOCATION')
       call prc ('                       CELLS WITH R<Rxpt HAVE '//
     >           ' OPPOSITE SIGN OF PINCH IN R>Rxpt CELLS')
       call prc ('                       PINCH IS ZERO IN CORE')
      ELSEIF (pinchopt.EQ.8) THEN
       CALL PRC ('  PINCH VEL. OPT   8 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH IS RADIAL')
       call prc ('                       PINCH IS PROJECTED ONTO'//
     >                                   ' S AND CROSS AXES')
       call prc ('                       BASED ON LOCAL CELL GEOMETRY')
      ELSEIF (pinchopt.EQ.9) THEN
       CALL PRC ('  PINCH VEL. OPT   9 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH IS RADIAL')
       call prc ('                       PINCH IS PROJECTED ONTO'//
     >                                   ' S AND CROSS AXES')
       call prc ('                       BASED ON LOCAL CELL GEOMETRY')
       call prc ('                       CORE PLASMA EXCLUDED')
      ELSEIF (pinchopt.EQ.10) THEN
       CALL PRC ('  PINCH VEL. OPT  10 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH IS RADIAL')
       call prc ('                       PINCH IS PROJECTED ONTO'//
     >                                   ' S AND CROSS AXES')
       call prc ('                       BASED ON LOCAL CELL GEOMETRY')
       call prc ('                       CORE PLASMA EXCLUDED')
       call prc ('                       REGION BELOW X-POINT EXCLUDED')
      ELSEIF (pinchopt.EQ.11) THEN
       CALL PRC ('  PINCH VEL. OPT  11 : ON - PINCH VELOCITY IS APPLIED
     >ONLY IN MAIN SOL ABOVE XPOINT')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
      ELSEIF (pinchopt.EQ.12) THEN
       CALL PRC ('  PINCH VEL. OPT  12 : ON - PINCH VELOCITY IS APPLIED
     >EVERYWHERE')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
       call prc ('                       PINCH IS RADIAL')
       call prc ('                       PINCH IS PROJECTED ONTO'//
     >                                   ' S AND CROSS AXES')
       call prc ('                       BASED ON LOCAL CELL GEOMETRY')
       call prc ('                       CORE PLASMA EXCLUDED')
       call prc ('                       REGION BELOW X-POINT EXCLUDED')
       call prc ('                       OUTER SOL EXCLUDED')
      ELSEIF (pinchopt.EQ.13) THEN
       CALL PRC ('  PINCH VEL. OPT  13 : ON - PINCH VELOCITY IS APPLIED
     >ONLY IN INNER SOL ABOVE XPOINT')
       call prr ('                       PINCH VELOCITY = ',cvpinch)
      ENDIF
C-----------------------------------------------------------------------
      IF (CPDRFT.EQ.0) THEN
       CALL PRC ('  POL DRIFT OPT    0 : OFF - NO ADDITIONAL POLOIDAL DR
     >IFT')
      ELSEIF (CPDRFT.EQ.1) THEN
c
       CALL PRC ('  POL DRIFT OPT    1 : ON  - ADDITIONAL DRIFT VELOCITY
     >')
       call prc ('                       APPLIED DIRECTLY TO ION TRANSPO
     >RT')
c
      ELSEIF (CPDRFT.EQ.2) THEN
c
       CALL PRC ('  POL DRIFT OPT    2 : ON  - ADDITIONAL DRIFT VELOCITY
     >')
       call prc ('                       APPLIED TO BACKGROUND PLASMA FL
     >OW')

c
      ELSEIF (CPDRFT.EQ.3) THEN
c
       CALL PRC ('  POL DRIFT OPT    3 : ON  - ADDITIONAL DRIFT VELOCITY
     >')
       call prc ('                       APPLIED DIRECTLY TO ION TRANSPO
     >RT')
       call prc('                        PFZ ONLY!')
       call prc('                        THIS OPTION HAS BEEN'//
     >                       ' DEPRECATED BY THE DRIFT REGION OPTION')
c
      ENDIF
c
c     Print out drift region interpretation and S values for each ring
c
      if (cpdrft.ne.0) then 

         if (drft_distopt.eq.0) then 
            write (coment,'(25x,''EXTENDING FROM: '',f6.2,'' * SMAX'')')
     >                   cdrftv_start
           call prc(coment)
           write (coment,'(25x,''            TO: '',f6.2,'' * SMAX'')')
     >                   cdrftv_end
           call prc(coment)
         elseif (drft_distopt.eq.1) then 

            write (coment,'(25x,''EXTENDING FROM: '',f6.2,'' * PMAX'')')
     >                   cdrftv_start
           call prc(coment)
           write (coment,'(25x,''            TO: '',f6.2,'' * PMAX'')')
     >                   cdrftv_end
           call prc(coment)


         elseif (drft_distopt.eq.2) then 

            write (coment,'(25x,''STARTING FROM:'',f6.2,'' (M) '',a)')
     >                   cdrftv_start,outer
           call prc(coment)
           write (coment,'(25x,''ENDING AT    :'',f6.2,'' (M) '',a)')
     >                   cdrftv_end,inner
           call prc(coment)

         endif

c
c        Write out drift region S values
c         
         call prc(s2//'TABLE OF DRIFT REGION BY RING - RINGS'//
     >                ' WITHOUT FLOW ARE NOT LISTED')
         call get_drftv_rings(irstart,irend)

         call prc(s2//'  IR      Vdrift (m/s)       '//
     >             ' S_START (m)       S_END (m)')

         do ir=irstart,irend

            write(coment,'(a,i6,1(1x,g12.5),2(1x,f16.5))')  s2,
     >           ir,pol_drftv(ir)/qtim,
     >           sdrft_start(ir),sdrft_end(ir)
            call prc(coment)

         end do 

      endif

c
c     Print out flow region information
c
      if (cpdrft.ne.0) then 
         if (drft_region.eq.1.and.cpdrft.ne.3) then 
            call prc(s2//'DRIFT REGION OPT 1: DRIFTS'//
     >                    ' APPLIED TO SOL+PFZ')
         elseif (drft_region.eq.2) then 
            call prc(s2//'DRIFT REGION OPT 2: DRIFTS'//
     >                    ' APPLIED TO SOL ONLY')
         elseif (drft_region.eq.3.or.cpdrft.eq.3) then 
            call prc(s2//'DRIFT REGION OPT 3: DRIFTS'//
     >                    ' APPLIED TO PFZ ONLY')
         elseif (drft_region.eq.4) then 
            call prc(s2//'DRIFT REGION OPT 4: DRIFTS'//
     >                    ' APPLIED TO CORE ONLY')
         endif
      endif
c
c     Print out flow data
c
      if (cpdrft.ne.0.and.ndrftvel.ne.0) then 
         call prc(sp//'DETAILED DRIFT VELOCITY SPECIFIED')
         call prc(sp//'DIFFERENT FLOW MAY BE SPECIFIED ON EVERY RING')
         call prc(sp//'VELOCITY IS ZERO ON UNLISTED RINGS:')
         call get_drftv_rings(irstart,irend)

         if (drftvel_machopt.gt.0) then
            call prr(sp//'DEFAULT FLOW MACH NUMBER = ',cdrftv)
            call pri('DRFIT VELOCITY SPECIFIED BY MACH VALUE. OPTION = '
     >               ,drftvel_machopt)
            call prc(sp//' IR      VEL (M/S)       CS(MID)     MACH') 
            do ir = irstart,irend
               write(coment,'(a,i5,3(2x,g12.5))') sp,ir,
     >                       pol_drftv(ir)/qtim,
     >                       ringcs(ir),pol_drftv(ir)/ringcs(ir)/qtim
               call prc(coment) 
            end do
         else
            call prr(sp//'DEFAULT FLOW RATE (M/S) = ',cdrftv)
            call prc(sp//' IR      VEL (M/S)') 
            do ir = irstart,irend
               write(coment,'(a,i5,3x,g12.5)') sp,ir,pol_drftv(ir)/qtim
               call prc(coment) 
            end do
         endif

         call prb

      elseif (cpdrft.ne.0.and.drftvel_machopt.gt.0) then 
c
         call prr(sp//'DEFAULT FLOW MACH NUMBER = ',cdrftv)
         call pri('DRFIT VELOCITY SPECIFIED BY MACH VALUE. OPTION = '
     >            ,drftvel_machopt)
         call prc(sp//' IR      VEL (M/S)       CS(MID)     MACH') 
c
         call get_drftv_rings(irstart,irend)
c
         do ir = irstart,irend
            write(coment,'(a,i5,3(2x,g12.5))') sp,ir,
     >                    pol_drftv(ir)/qtim,
     >                    ringcs(ir),pol_drftv(ir)/ringcs(ir)/qtim
            call prc(coment) 
         end do
         call prb
      elseif (cpdrft.ne.0.and.drftvel_machopt.eq.0) then 
         CALL PRR (sp//'DEFAULT FLOW VELOCITY (m/s) = ',CDRFTV)
      endif
c
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_bg_options(NIZS,NIMPS,NIMPS2,nymfs)
      IMPLICIT none
c
      INTEGER NIZS,NIMPS,NIMPS2,nymfs
C
C***********************************************************************
C     THIS ROUTINE PRINTS ALL THE OPTION FLAGS, TORUS ATTRIBUTES ETC
C     STORED IN COMMONS COMTOR AND COMTAU
C
C     C.M.FARRELL    NOVEMBER 1987
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
      include    'cedge2d'
c      include    'dynam4'
      include    'dynam5'
      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
      include    'slcom'  
      include    'printopt'
c slmod begin - new
      INTEGER i1,i2
      LOGICAL status
      COMMON /NEWCOM/ new
      LOGICAL         new
c slmod end
c
c      CHARACTER  COMENT*80,prtype*4
      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
      real     totmhpin
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('--- BACKGROUND PLASMA OPTIONS --- ')
      CALL PRChtml ('--- BACKGROUND PLASMA OPTIONS ---',
     >              'pr_bg_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------
c
c     PLASMA DECAY OPTIONS
c
      IF     (CIOPTG.EQ.0) THEN
       CALL PRC ('  PLASMA DECAY OPT 0 : STANDARD')
      ELSEIF (CIOPTG.EQ.1) THEN
       CALL PRC ('  PLASMA DECAY OPT 1 : EXPONENTIAL DECAY OUTBOARD')
      ELSEIF (CIOPTG.EQ.2) THEN
       CALL PRC ('  PLASMA DECAY OPT 2 : FROM INPUT DATA FOR RINGS IN SO
     >L')
      ELSEIF (CIOPTG.EQ.3) THEN
       CALL PRC ('  PLASMA DECAY OPT 3 : FROM INPUT DATA FOR RINGS IN SO
     >L')
       CALL PRC ('                       INNER AND OUTER PLATE DIFFER')
      ELSEIF (CIOPTG.EQ.4) THEN
       CALL PRC ('  PLASMA DECAY OPT 4 : FROM INPUT DATA FOR RINGS IN SO
     >L AND TRAP')
       CALL PRC ('                       INNER AND OUTER PLATE DIFFER')
      ELSEIF (CIOPTG.EQ.5) THEN
       CALL PRC ('  PLASMA DECAY OPT 5 : EXPONENTIAL DECAY OUTBOARD')
      ELSEIF (CIOPTG.EQ.6) THEN
       CALL PRC ('  PLASMA DECAY OPT 6 : EXPONENTIAL DECAY OUTBOARD')
       call prc ('                       DIFFERENT EXPONENTIAL FACTORS')
       call prc ('                       IN SOL AND PRIVATE PLASMA')
      ELSEIF (CIOPTG.EQ.7) THEN
       CALL PRC ('  PLASMA DECAY OPT 7 : FROM INPUT DATA FOR RINGS IN SO
     >L AND TRAP AND CORE')
       CALL PRC ('                       INNER AND OUTER PLATE DIFFER')
       call prc ('                       VALUES SPECIFIED FOR CORE MAY')
       CALL PRC ('                       OVERRIDE THE PRESCRIBED CONDITI
     >ONS')
       CALL PRC ('                       LISTED ABOVE FOR CTEBIN,CTIBIN
     >AND')
       CALL PRC ('                       CNBIN. THESE VALUES ARE STILL')
       CALL PRC ('                       USED IF CORE DECAY OPTION 1 IS'
     >)
       CALL PRC ('                       SPECIFIED')
      ELSEIF (CIOPTG.EQ.90) THEN
       CALL PRC ('  PLASMA DECAY OPT 90: PIECE-WISE BACKGROUND PLASMA')
       call prc ('                       BASE PLASMA FROM EQUILIBRIUM DA
     >TA')
       call prc ('                       SPECIFIC RING SECTIONS WERE OVE
     >RLAID WITH')
       call prc ('                       THE FOLLOWING OPTIONS.')
       if (nbgplas.gt.0) then
c
         call prb
         call prc ('  Section: 1 = '//OUTER//' 2 = '//INNER//
     >              ' 3 = WHOLE RING')
         call prc ('  Ring1 to Ring2  Sect   PlasDec   SOL'//
     >                  '   TeGrad  TiGrad   Core   E-field')
       endif
       do in = 1,nbgplas
          write(coment,'(3x,i5,3x,i5,3x,i2,2x,6(3x,i5))')
     >                 (int(bgplasopt(in,i)),i=1,9)
          call prc (coment)
       end do
       call prb
      ELSEIF (CIOPTG.EQ.91) THEN
       CALL PRC ('  PLASMA DECAY OPT 91: PIECE-WISE BACKGROUND PLASMA')
       call prc ('                       BASE PLASMA FROM REGULAR OPTION
     >S.')
       call pri ('                       PLASMA DECAY = ',4)
       call pri ('                       SOL OPTION   = ',cioptf)
       call pri ('                       Te Grad Opt  = ',cioptk)
       call pri ('                       Ti Grad Opt  = ',cioptl)
       call pri ('                       CORE OPTION  = ',ccoreopt)
       call pri ('                       E-FIELD OPT  = ',ofield)
       call prb
       call prc ('                       SPECIFIC RING SECTIONS WERE OVE
     >RLAID WITH')
       call prc ('                       THE FOLLOWING OPTIONS.')
       if (nbgplas.gt.0) then
         call prb
         call prc ('  Section: 1 = '//OUTER//' 2 = '//INNER//
     >              ' 3 = WHOLE RING')
         call prc ('  Ring1 to Ring2  Sect   PlasDec   SOL'//
     >                  '   TeGrad  TiGrad   Core   E-field')
       endif
       do in = 1,nbgplas
          write(coment,'(3x,i5,3x,i5,3x,i2,2x,6(3x,i5))')
     >                 (int(bgplasopt(in,i)),i=1,9)
          call prc (coment)
       end do
       call prb
      ELSEIF (CIOPTG.EQ.92) THEN
       CALL PRC ('  PLASMA DECAY OPT 92: PIECE-WISE BACKGROUND PLASMA')
       call prc ('                       BASE PLASMA FROM PREVIOUS DIVIM
     >P RUN')
       call prc ('                       SPECIFIC RING SECTIONS WERE OVE
     >RLAID WITH')
       call prc ('                       THE FOLLOWING OPTIONS.')
       if (nbgplas.gt.0) then
         call prb
         call prc ('  Section: 1 = '//OUTER//' 2 = '//INNER//
     >              ' 3 = WHOLE RING')
         call prc ('  Ring1 to Ring2  Sect   PlasDec   SOL'//
     >                  '   TeGrad  TiGrad   Core   E-field')
       endif
       do in = 1,nbgplas
          write(coment,'(3x,i5,3x,i5,3x,i2,2x,6(3x,i5))')
     >                 (int(bgplasopt(in,i)),i=1,9)
          call prc (coment)
       end do
       call prb
      ELSEIF (CIOPTG.EQ.98) THEN
       CALL PRC ('  PLASMA DECAY OPT 98: FROM DIVIMP PLASMA INPUT FILE')
      ELSEIF (CIOPTG.EQ.99) THEN
       CALL PRC ('  PLASMA DECAY OPT 99: FROM EQUILIBRIUM DATA')
       if (ciopto.eq.2.or.ciopto.eq.3.or.ciopto.eq.4) then
        call prc('                       WARNING!!!! - THE PRIVATE PLASM
     >A')
        call prc('                       CONDITIONS FROM EQUILIBRIUM DAT
     >A')
        call prc('                       HAVE BEEN REPLACED BY A PLASMA'
     >)
        call pri('                       SPECIFIED BY PRIVATE PLASMA OPT
     >ION : ',ciopto)
       endif
      ENDIF
c
C-----------------------------------------------------------------------
c
c     FLUID CODE Data Options - Calculation of target conditions
c
      if (cioptg.eq.90.or.cioptg.eq.99.or.cre2d.ne.0) then  
c

      if (fc_target_calc_option.eq.0) then 
       CALL PRC ('  FC BASE TARG OPT 0 : FLUID CODE - BASE Target'//
     >' Options: EDGE2D')
       call prc ('                       NE CALC OPT = 2') 
       call prc ('                       TE CALC OPT = 1') 
       call prc ('                       TI CALC OPT = 2') 
       call prc ('                       VB CALC OPT = 0') 
      elseif (fc_target_calc_option.eq.1) then 
       CALL PRC ('  FC BASE TARG OPT 1 : FLUID CODE - BASE Target'//
     >' Options: UEDGE')
       call prc ('                       NE CALC OPT = 2') 
       call prc ('                       TE CALC OPT = 0') 
       call prc ('                       TI CALC OPT = 0') 
       call prc ('                       VB CALC OPT = 0') 
      elseif (fc_target_calc_option.eq.2) then 
       CALL PRC ('  FC BASE TARG OPT 2 : FLUID CODE - BASE Target'//
     >' Options: AUG/DIV-B2')
       call prc ('                       NE CALC OPT = 2') 
       call prc ('                       TE CALC OPT = 2') 
       call prc ('                       TI CALC OPT = 2') 
       call prc ('                       VB CALC OPT = 1') 
      elseif (fc_target_calc_option.eq.3) then 
       CALL PRC ('  FC BASE TARG OPT 3 : FLUID CODE - BASE Target'//
     >' Options: B2/B2.5')
       call prc ('                       NE CALC OPT = 2') 
       call prc ('                       TE CALC OPT = 1') 
       call prc ('                       TI CALC OPT = 1') 
       call prc ('                       VB CALC OPT = 1') 
      endif
c
c     Detailed printouts of the target data interpretation options 
c
c     Ne
c
      if (fc_ne_calc_opt.eq.0) then 
       CALL PRC ('  FC NE CALC OPT 0   : FLUID CODE Target Density taken
     > from value in guard cell')
      elseif (fc_ne_calc_opt.eq.1) then 
       CALL PRC ('  FC NE CALC OPT 1   : FLUID CODE Target Density is ar
     >ithmetic mean of guard and first real cell')
      elseif (fc_ne_calc_opt.eq.2) then 
       CALL PRC ('  FC NE CALC OPT 2   : FLUID CODE Target Density taken
     > from first real cell')
      endif
c 
c     Te
c
      if (fc_te_calc_opt.eq.0) then 
       CALL PRC ('  FC TE CALC OPT 0   : FLUID CODE Target Te taken'//
     >' from value in guard cell')
      elseif (fc_te_calc_opt.eq.1) then 
       CALL PRC ('  FC TE CALC OPT 1   : FLUID CODE Target Te is ar'//
     >'ithmetic mean of guard and first real cell')
      elseif (fc_te_calc_opt.eq.2) then 
       CALL PRC ('  FC TE CALC OPT 2   : FLUID CODE Target Te taken'//
     >' from first real cell')
      endif
c
c     Ti
c
      if (fc_ti_calc_opt.eq.0) then 
       CALL PRC ('  FC TI CALC OPT 0   : FLUID CODE Target Ti taken'//
     >' from value in guard cell')
      elseif (fc_ti_calc_opt.eq.1) then 
       CALL PRC ('  FC TI CALC OPT 1   : FLUID CODE Target Ti is ar'//
     >'ithmetic mean of guard and first real cell')
      elseif (fc_ti_calc_opt.eq.2) then 
       CALL PRC ('  FC TI CALC OPT 2   : FLUID CODE Target Ti taken'//
     >' from first real cell')
      endif
c
c     Vb
c
      if (fc_v_calc_opt.eq.0) then 
       CALL PRC ('  FC VB CALC OPT 0   : FLUID CODE Target Velocity'//
     >' is set to FC boundary value')
      elseif (fc_v_calc_opt.eq.1) then 
       CALL PRC ('  FC VB CALC OPT 1   : FLUID CODE Target Velocity'//
     >' is set to sound speed')
      endif
c
c
c     Interpretation of velocity in the fort.31 file for B2/B2.5
c     formatted plasmas.
c
      if (cgridopt.eq.3.or.cgridopt.eq.GEN_GRID) then 
c
       if (fc_v_interp_opt.eq.0) then
       CALL PRC ('  FC VEL OPTION 0    : FLUID CODE Velocity'//
     >' Interpretation Option 0')
       call prc ('                       Velocity in plasma file is a'//
     >' cell boundary quantity')
       elseif (fc_v_interp_opt.eq.1) then
       CALL PRC ('  FC VEL OPTION 1    : FLUID CODE Velocity'//
     >' Interpretation Option 1')
       call prc ('                       Velocity in plasma file is a'//
     >' cell center quantity')
       endif 
c
      endif 

c
c     ENDIF for (cioptg.eq.90.or.cioptg.eq.99.or.cre2d.ne.0)  
c
      endif
c
C-----------------------------------------------------------------------
c
c     FLUID CODE options related to OSM target conditions 
c
      if (cre2d.eq.0) then
       CALL PRC ('  FC READ OPT 0      : FLUID CODE Background is NOT re
     >ad')
       call prc ('                       in separately for reference.')
      elseif (cre2d.eq.1) then
       CALL PRC ('  FC READ OPT 1      : EDGE2D Background is read')
       call prc ('                       in separately for reference.')
      elseif (cre2d.eq.2) then
       CALL PRC ('  FC READ OPT 2      : UEDGE DATA is read')
       call prc ('                       in separately for reference.')
c slmod begin
      elseif (cre2d.eq.3) then
       CALL PRC ('  FC READ OPT 2      : B2 DATA (Rhozansky) is read')
       call prc ('                       in separately for reference.')
c slmod end
c
      endif
c
C-----------------------------------------------------------------------
c     Target condition options
c
      if (cre2d.eq.1.or.cre2d.eq.2) then 
c
c      USE FLUID CODE DATA for OSM TARGET conditions
c
       if (e2dtargopt.eq.0) then
c
        CALL PRC('  FC TARG OPT 0      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS NOT USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS.')
c
       elseif (e2dtargopt.eq.1) then
c
        CALL PRC('  FC TARG OPT 1      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS.')
c
       elseif (e2dtargopt.eq.2) then
c
        CALL PRC('  FC TARG OPT 2      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS')
        CALL PRC('                       AND ADDITIONAL TARGET FLUX INFO
     >RMATION.')
        CALL PRC('                       THE TARGET DENSITY IS CALCULATE
     >D FROM THE')
        CALL PRC('                       GIVEN FLUXES: Nt= G / CS ')
        call prb
        call prc('          Imposed Target Fluxes:')
        call prc('          Ring         '//INNER//'         '//OUTER)
        do in = 1,fluxpts
           write (coment,'(11x,i4,6x,g16.8,6x,g16.8)')
     >                         int(fluxinfo(in,1)),
     >                         fluxinfo(in,2),fluxinfo(in,3)
           call prc(coment)
        end do
c
       elseif (e2dtargopt.eq.3) then
c
        CALL PRC('  FC TARG OPT 3      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS.')
        CALL PRC('                       THE ADDITIONAL TARGET FLUX INFO
     >RMATION')
        CALL PRC('                       IS PASSED DIRECTLY TO PIN TO DE
     >FINE THE')
        CALL PRC('                       ACTUAL TARGET FLUXES. THE TARGE
     >T CONDITIONS')
        CALL PRC('                       WILL NOT MATCH THESE FLUXES.')
        CALL PRB
        call prc('          Imposed Target Fluxes:')
        call prc('          Ring         '//INNER//'         '//OUTER)
        do in = 1,fluxpts
           write (coment,'(11x,i4,6x,g16.8,6x,g16.8)')
     >                         int(fluxinfo(in,1)),
     >                         fluxinfo(in,2),fluxinfo(in,3)
           call prc(coment)
        end do
c
       elseif (e2dtargopt.eq.4) then
c
        CALL PRC('  FC TARG OPT 4      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS')
        CALL PRC('                       AND ADDITIONAL TARGET FLUX INFO
     >RMATION.')
        CALL PRC('                       THE TARGET VELOCITY IS CALCULAT
     >ED FROM THE')
        CALL PRC('                       GIVEN FLUXES: Vt= G / NB ')
        call prb
        call prc('          Imposed Target Fluxes:')
        call prc('          Ring         '//INNER//'         '//OUTER)
        do in = 1,fluxpts
           write (coment,'(11x,i4,6x,g16.8,6x,g16.8)')
     >                         int(fluxinfo(in,1)),
     >                         fluxinfo(in,2),fluxinfo(in,3)
           call prc(coment)
        end do
c
c
       elseif (e2dtargopt.eq.5) then
c
        CALL PRC('  FC TARG OPT 5      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS')
        CALL PRC('                       AND ADDITIONAL TARGET FLUX INFO
     >RMATION.')
        CALL PRC('                       THE TARGET DENSITY IS CALCULATE
     >D FROM THE')
        CALL PRC('                       GIVEN FLUXES: Nt= G / CS ')
        CALL PRC('                       EDGE2D TARGET POWER FLUXES FOR'
     >)
        CALL PRC('                       ELECTRONS AND IONS ARE ALSO IMP
     >OSED')
        CALL PRC('                       WHERE RELEVANT (E.G.SOL22)')
        call prb
        call prc('          Imposed Target Particle Fluxes:')
        call prc('          Ring         '//INNER//'         '//OUTER)
        do in = 1,fluxpts
           write (coment,'(11x,i4,6x,g16.8,6x,g16.8)')
     >                         int(fluxinfo(in,1)),
     >                         fluxinfo(in,2),fluxinfo(in,3)
           call prc(coment)
        end do
c
       elseif (e2dtargopt.eq.6) then
c
        CALL PRC('  FC TARG OPT 6      : EDGE2D DATA THAT HAS BEEN READ'
     >)
        CALL PRC('                       IS USED TO ASSIGN INITIAL')
        CALL PRC('                       TARGET CONDITIONS BASED ON THE'
     >)
        CALL PRC('                       BACKGROUND PLASMA CONDITIONS.')
        call prc('                       Vb=Cs AT TARGET IS FORCED')
c
       endif

c
c      Endif for cre2d = 1 or 2
c
       endif
c
C-----------------------------------------------------------------------
c
      if (cioptg.eq.99.and.cgridopt.eq.0) then 
c
      call prb  
      call prc('EDGE2D EQUILIBRIUM SOLUTION OPTIONS:')

C-----------------------------------------------------------------------
      IF     (CNIWA.EQ.0) THEN
       CALL PRC('  LOST SOL RINGS   0 : PLASMA SET TO MINIMUM VALUES')
      ELSEIF (CNIWA.EQ.1) THEN
       CALL PRC('  LOST SOL RINGS   1 : PLASMA SET TO OUTER RING VALUE')
      ENDIF
C-----------------------------------------------------------------------
c
      IF(cprint.eq.9 ) CALL PRVMF
c
       IF( CNVMF.LE.0 ) THEN
           CALL PRC('  VMF                : 1.000 ')
       ELSE
          DO 900 I = 1 , CNVMF
             CALL PRI2('  VMF RING RANGE     :' ,CIRNG0(I) ,CIRNG1(I) )
             CALL PRI2('  J0 & J1            :' ,CJ0(I)    ,CJ1(I)    )
             CALL PRR3('  VMF0,VMF1,VMF2     :'
     >                ,CVMF0(I) ,CVMF1(I) ,CVMF2(I) )
  900     CONTINUE
       END IF

      endif 
C-----------------------------------------------------------------------
      IF     (CIOPTF.EQ.-1) THEN
       WRITE (COMENT,'(''  SOL OPTION      -1 : SOL1A,  (FL,FS) = ('',
     >   F6.3,'','',F6.3,'')'')') CFL,CFS
       CALL PRC (COMENT)
      endif

      IF (CIOPTF.EQ.0) THEN
       CALL PRC ('  SOL OPTION       0 : SOL0,  VHOUT=0,  EOUT=0')
      endif

      IF (CIOPTF.EQ.1) THEN
       CALL PRC ('  SOL OPTION       1 : SOL1')
      endif

      IF (CIOPTF.EQ.2) THEN
       CALL PRC ('  SOL OPTION       2 : SOL2')
      endif

      IF (CIOPTF.EQ.3) THEN
       CALL PRC ('  SOL OPTION       3 : SOL3')
      endif

      IF (CIOPTF.EQ.4) THEN
       CALL PRC ('  SOL OPTION       4 : SOL4')
      endif

      IF (CIOPTF.EQ.5) THEN
       WRITE (COMENT,'(''  SOL OPTION       5 : SOL5,  VHOUT='',
     >   1P,G11.4,'',  EOUT='',G11.4)') CVHOUT,CEOUT
       CALL PRC (COMENT)
      endif

      if (cioptf.eq.6) then
       CALL PRc ('  SOL OPTION       6 : SOL6, FIXED E-field, Variable V
     >')
       call prr ('                       Eout = ',ceout)
       call prr ('                       Base Velocity Vhout (V0) = ',
     >           cvhout)
       write (coment,'(24x,''Fixed at V0 to '',g8.3,
     >'' * SMAX'')')   cvbl1
       call prc(coment)
       write (coment,'(24x,''Then fixed at '',g8.3,'' *V0 to '',g8.3,
     >'' * SMAX'')') cvbm1,cvbl2
       call prc(coment)
       write (coment,'(24x,''Then fixed at '',g8.3,'' *V0 to SMAX/2.0'')
     > ')       cvbm2
       call prc(coment)
      endif

      if (cioptf.eq.7) then
       CALL PRc ('  SOL OPTION       7 : SOL7, FIXED E-field, Linear Vb
     >')
       call prr ('                       Eout = ',ceout)
       call prr ('                       Base Velocity Vhout (V0) = ',
     >           cvhout)
       write (coment,'(24x,''Falling to '',g8.3,'' *V0 at '',g8.3,
     >'' * SMAX'')') cvbm1,cvbl1
       call prc(coment)
       write (coment,'(24x,''Then to '',g8.3,'' *V0 at '',g8.3,
     >'' * SMAX'')') cvbm2,cvbl2
       call prc(coment)
       call prc ('                       Then to Vb = 0 at SMAX/2 ')
      endif

      IF (CIOPTF.EQ.9) THEN
       CALL PRC ('  SOL OPTION       9 : SOL9')
      endif

      IF (CIOPTF.EQ.10) THEN
       CALL PRR ('  SOL OPTION      10 : SOL10, REVERSAL MACH NUMBER',
     >    CFRM)
       WRITE (COMENT,'(''                       KIN'',F5.2,'',KOUT'',
     >   F5.2,'',FRMIN'',F6.3,'',FRMAX'',F6.3)')CKIN,CKOUT,CFRMIN,CFRMAX
       CALL PRC (COMENT)
       WRITE (7,'(24X,''(FL,FS) = ('',F6.3,'','',F6.3,'')'')') CFL,CFS
      endif

      IF (CIOPTF.EQ.12) THEN
       CALL PRC ('  SOL OPTION      12 : SOL12 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/2K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(2KT+MV**2) = 4NOKTO')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      endif

      IF (CIOPTF.EQ.13) THEN
       CALL PRC ('  SOL OPTION      13 : SOL13 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/2K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(KTI+KTE+MV**2) = 2N0K(T0i+T0
     >e)')
       CALL PRC ('                       P/A = (5KTEBP)*NBP*CSBP FOR EL
     >ECTRONS')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRC ('                       P/A = (2KTIBP)*NBP*CSBP FOR IO
     >NS')
       CALL PRC ('                       PR/A= 0 FOR IONS')
       CALL PRR ('                       K0I = ', CK0I)
      endif

      IF (CIOPTF.EQ.14) THEN
       CALL PRC ('  SOL OPTION      14 : SOL14 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       HEAT CONDUCTION AND CONVECTION'
     >)
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       D/DS5NVKT-K0*T**(5/2)*(DT/DS)
     >=-PR(S)/A')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(2KT+MV**2) = 4NOKTO')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      endif

      IF (CIOPTF.EQ.15) THEN
       CALL PRC ('  SOL OPTION      15 : SOL15 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/4K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(KTI+KTE+MV**2) = 2N0K(T0i+T0
     >e)')
       CALL PRC ('                       P/A = (5KTEBP)*NBP*CSBP FOR EL
     >ECTRONS')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRC ('                       P/A = (2KTIBP)*NBP*CSBP FOR IO
     >NS')
       CALL PRC ('                       PR/A= 0 FOR IONS')
       CALL PRR ('                       K0I = ', CK0I)
      endif

      if (cioptf.eq.16) then
       CALL PRC ('  SOL OPTION      16 : SOL16 - AS SOL 14 - Linear Solu
     >tions')
      endif

      if (cioptf.eq.17) then
       CALL PRC ('  SOL OPTION      17 : SOL17 - AS SOL 16 - Te and Ti')
       call prc ('                       are calculated independently.')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRR ('                       K0I = ', CK0I)
      endif

      if (cioptf.eq.18) then
       CALL PRC ('  SOL OPTION      18 : SOL18 - AS SOL 16 - except')
       call prc ('                       the 1/2mv**2*Gamma term has')
       call prc ('                       been added.')
      endif

      if (cioptf.eq.19) then
       CALL PRC ('  SOL OPTION      19 : SOL19 - AS SOL 18 - except')
       call prc ('                       the hydrogenic rad term has')
       call prc ('                       been added.')
      endif

      if (cioptf.eq.20) then
       CALL PRC ('  SOL OPTION      20 : SOL20 - AS SOL 17 - except')
       call prc ('                       the hydrogenic rad term and')
       call prc ('                       Pei terms have been added.')
      endif

      if (cioptf.eq.21.or.prsol21) then
       CALL PRC ('  SOL OPTION      21 : Detached plasma model')
       if (s21refsw.eq.0) then
        call prc('                       NOTE: REF OPTION 0 - ')
        call prc('                       All distances are specified')
        call prc('                       in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc('                       NOTE: REF OPTION 1 - ')
        call prc('                       All distances are specified')
        call prc('                       in units of PMAX. (MF=PMAX)')
        call prc('                       Converted to S for each'//
     >                                 ' ring.')
       elseif (s21refsw.eq.2) then
        call prc('                       NOTE: REF OPTION 2 - ')
        call prc('                       All distances are specified')
        call prc('                       in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc('                       NOTE: REF OPTION 3 - ')
        call prc('                       All distances are specified')
        call prc('                       in meters along P.(MF=1.0)')
        call prc('                       Converted to S for each'//
     >                                 ' ring.')
       endif
c
c      OUTER
c
       call prc ('                       '//OUTER//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rat)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rat,l2rat
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rat)
       call prr ('                       Te at end of A: Te0 * ',
     >            terat)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirat)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrat)
c
       if (nrat.lt.0.0) then

        call prc('                       Nratio < 0 - variable N may'//
     >' be used to match pressure.')
        call prc('             IR          Nratio  for '//OUTER)
        do ir = irsep,nrs
           write(coment,'(12x,i5,10x,f10.2)') ir,nrat_used(ir,2)
           call prc(coment)
        end do
        call prb
       endif
c
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalph)
       call prr ('                       Radiated power in B (Qrad): Q0
     >* ',qrat)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases lineraly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrat)
c
c      INNER
c
       call prb
       call prc ('                       '//INNER//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rati)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rati,l2rati
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rati)
       call prr ('                       Te at end of A: Te0 * ',
     >            terati)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirati)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrati)
c
       if (nrati.lt.0.0) then

        call prc('                       Nratio < 0 - variable N may'//
     >' be used to match pressure.')
        call prc('             IR          Nratio  for '//INNER)
        do ir = irsep,nrs
           write(coment,'(12x,i5,10x,f10.2)') ir,nrat_used(ir,1)
           call prc(coment)
        end do
        call prb
       endif
c
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalphi)
       call prr ('                       Radiated power in B (Qrad): Q0
     >* ',qrati)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases linearly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrati)

c
c     Check for arrays of parameters specified
c
       if (ns21i.gt.0) then
c
         call prc('    SOL21 - Per ring parameters - '//INNER)
         call prc('    These parameters over ride defaults above')
         call prc('  IR   Te/Te0  Ti/Ti0  Ne/Ne0  Qr/Q0  L1/SM  '
     >            //'L2/SM  LV/SM  Vmult')
         do in = 1,ns21i
            write(coment,'(i4,8(1x,f8.4))') int(s21parmi(in,1)),
     >              (s21parmi(in,i),i=2,9)
            call prc(coment)
         end do
       end if
c
       if (ns21o.gt.0) then
c
         call prc('    SOL21 - Per ring parameters - '//OUTER)
         call prc('    These parameters over ride defaults above')
         call prc('  IR   Te/Te0  Ti/Ti0  Ne/Ne0  Qr/Q0  L1/SM  '
     >            //'L2/SM  LV/SM  Vmult')
         do in = 1,ns21o
            write(coment,'(i4,8(1x,f8.4))') int(s21parmo(in,1)),
     >              (s21parmo(in,i),i=2,9)
            call prc(coment)
         end do
       end if
c
c
c     Check for arrays of extra parameters specified
c
       if (aux_ns21i.gt.0) then
c
         call prc('    SOL21 - EXTRA Per ring parameters - '//INNER)
         call prc('    These parameters over ride defaults above')
         call prc('  IR   LR1A LR1B  NeR(A) NeR(B) TeR(A) TeR(B)'//
     >            ' TiR(A) TiR(B)')
c
         do in = 1,aux_ns21i
            write(coment,'(i4,8(1x,f8.4))') int(aux_s21parmi(in,1)),
     >              (aux_s21parmi(in,i),i=2,9)
            call prc(coment)
         end do
       end if
c
       if (aux_ns21o.gt.0) then
c
         call prc('    SOL21 - EXTRA Per ring parameters - '//OUTER)
         call prc('    These parameters over ride defaults above')
         call prc('  IR   LR1A LR1B  NeR(A) NeR(B) TeR(A) TeR(B)'//
     >            ' TiR(A) TiR(B)')
c
         do in = 1,aux_ns21o
            write(coment,'(i4,8(1x,f8.4))') int(aux_s21parmo(in,1)),
     >              (aux_s21parmo(in,i),i=2,9)
            call prc(coment)
         end do
c
       end if
c
c
      endif

      if (cioptf.eq.22.or.prsol22) then
       CALL PRC ('  SOL OPTION      22 : Runge Kutta SOL equation solver
     >')
       call echosol
      endif

      if (cioptf.eq.23) then
       CALL PRC ('  SOL OPTION      23 : CFD ring by ring plasma solver'
     >)
       call echosol23
      endif

      IF (CIOPTF.EQ.98) THEN
       CALL PRC ('  SOL OPTION      98 : FROM DIVIMP PLASMA INPUT FILE')
       call prc ('                       THIS ONLY READS IN THE ENTIRE')
       call prc ('                       BACKGROUND AND CAN ONLY BE SPECI
     >FIED')
       call prc ('                       IN COMBINATION WITH PLASMA DECAY
     >')
       call prc ('                       OPTION 98.')
      endif

      IF (CIOPTF.EQ.99) THEN
       CALL PRC ('  SOL OPTION       99: FROM EQUILIBRIUM DATA')
      ENDIF

      IF (CSOLEF.NE.1.0.OR.CSOLVF.NE.1.0) THEN
       WRITE (COMENT,'(23X,2(A,F10.4))') 'E FACTOR',CSOLEF,
     >                                ',  V FACTOR',CSOLVF
       CALL PRC (COMENT)
      ENDIF



C-----------------------------------------------------------------------
      IF (CIOPTF.EQ.12.OR.CIOPTF.EQ.13.OR.CIOPTF.EQ.14.OR.
     >    CIOPTF.EQ.15.or.cioptf.eq.16.or.cioptf.eq.17.or.
     >    cioptf.eq.18.or.cioptf.eq.19.or.cioptf.eq.20) THEN


       IF (CSOPT.EQ.0) THEN
       CALL PRR ('  SOL  :IONIZATION 0 : CONSTANT OVER 0 < S < SMAX *'
     >           ,CSOLLS)
       CALL PRC ('                       WITH SO = -NOVO')
       ELSEIF (CSOPT.EQ.1) THEN
       CALL PRC ('  SOL  :IONIZATION 1 : EXPONENTIAL DECAY DESCRIBED BY:
     > ')
       CALL PRC ('                       SI(S) = SO * EXP(-S/LS)')
       CALL PRR ('                       WITH DECAY LENGTH LS = SMAX * '
     >         ,CSOLLS)
       CALL PRC ('                       SO= -NOVO/(LS*(1-EXP(-L/LS)))')
       CALL PRC ('                       L = SMAX/2.0 FOR FIELD LINE')
       ELSEIF (CSOPT.EQ.4) THEN
       CALL PRC ('  SOL  :IONIZATION 4 : COMBINATION OF TWO CONSTANT SOU
     >RCES')
       CALL PRC ('                       S(S) = -(1-F)NOVO/LS - FNOVO/L'
     >)
       CALL PRC ('                       FOR 0 < S < LS')
       CALL PRC ('                       S(S) = - FNOVO/L')
       CALL PRC ('                       FOR LS < S < L')
       CALL PRR ('                       WHERE LS IS SMAX * ',CSOLLS)
       CALL PRR ('                       WHERE L  IS SMAX * ',CSOLLT)
       CALL PRR ('                       WHERE F  IS        ',CFIZ)
       ELSEIF (CSOPT.EQ.5) THEN
       CALL PRC ('  SOL  :IONIZATION 5 : COMBINATION OF CONSTANT AND EXP
     >ONENTIAL')
       CALL PRC ('                       S(S) = SO *EXP(-S/LS) - FNOVO/L
     >')
       CALL PRC ('                       SO= -(1-F)NOVO/(LS*(1-EXP(-L/LS
     >)))')
       CALL PRR ('                       WHERE LS IS SMAX * ',CSOLLS)
       CALL PRR ('                       WHERE L  IS SMAX * ',CSOLLT)
       CALL PRR ('                       WHERE F  IS        ',CFIZ)
       ENDIF
c
c
c
       IF (CPOPT.EQ.0) THEN
       CALL PRR ('  SOL  :RADIATION  0 : CONSTANT OVER 0 < S < SMAX * '
     >          ,CSOLLR)
       CALL PRR ('                       WITH PR/A = ',CSOLPR)
       ELSEIF (CPOPT.EQ.1) THEN
       CALL PRC ('  SOL  :RADIATION  1 : EXPONENTIAL DECAY DESCRIBED BY:
     > ')
       CALL PRC ('                       PR(S)/A = PRO/A * EXP(S/LR)')
       CALL PRR ('                       WITH DECAY LENGTH LR = SMAX * '
     >          ,CSOLLR)
       CALL PRR ('                       PRO/A = ',CSOLPR)

       ELSEIF (CPOPT.EQ.2) THEN
       CALL PRR ('  SOL  :RADIATION  2 : CONSTANT OVER 0 < S < SMAX * '
     >          ,CSOLLR)
       CALL PRC ('                       WITH PR/A = FRR * ( P/A ) /LR'
     >          )
       CALL PRR ('                       WHERE:   FRR = ',CSOLFR)
       CALL PRR ('                       AND      LR  = SMAX * ',CSOLLR)
       ELSEIF (CPOPT.EQ.3) THEN
       CALL PRC ('  SOL  :RADIATION  3 : EXPONENTIAL DECAY DESCRIBED BY:
     > ')
       CALL PRC ('                       PR(S)/A = PRO/A * EXP(-S/LR)')
       CALL PRR ('                       WITH DECAY LENGTH LR = SMAX * '
     >          ,CSOLLR)
       CALL PRC ('                       PRO/A=FRR * (P/A)/(LR*(1-EXP(-S
     >MAX/2LR)))')
       CALL PRR ('                       WHERE:   FRR = ',CSOLFR)
       ENDIF
c
c
c
       IF (FLUXROPT.EQ.0) THEN
       CALL PRC ('  SOL  : FLUX OPT  0 : NO FLUX RECIRCULATION')
       ELSEIF (FLUXROPT.EQ.1) THEN
       CALL PRC ('  SOL  : FLUX OPT  1 : IONIZATION SOURCE MODIFIED')
       CALL PRC ('                       FOR FLUX RECIRCULATION.')
       CALL PRC
     >('                       RING  SRCMULT  SRCLEN*SMAX  SRCLAM*SMAX')
       DO 370 I = 1,FLUXPTS
         WRITE(COMENT,'(23X,I4,3(2X,G9.3))') INT(FLUXINFO(I,1)),
     >           FLUXINFO(I,2),FLUXINFO(I,3),FLUXINFO(I,4)
         CALL PRC(COMENT)
 370   CONTINUE
       ENDIF
c
c
c
       IF (SROOTOPT.EQ.0) THEN
       CALL PRC('  SOL  : ROOT OPT  0 : IMAGINARY VALUES ARE SET TO ZERO
     >')
       ELSEIF (SROOTOPT.EQ.1) THEN
       CALL PRC('  SOL  : ROOT OPT  1 : IMAGINARY VALUES FORCE THE FLOW'
     >)
       CALL PRC('                       VELOCITY TO BE SET TO THE LOCAL'
     >)
       CALL PRC('                       SOUND SPEED.')
       ENDIF
c
c
c
       if (ciopto.eq.0.or.ciopto.eq.2.or.
     >     ciopto.eq.3.or.ciopto.eq.4) then
          irlim = irwall
       elseif (ciopto.eq.1) then
          irlim = nrs
       endif
c
c      Only print for debugging purposes - code will be deleted later
c      if it will never be used.
c
       if (cprint.eq.5.or.cprint.eq.9) then
c
       call prc('  SOL  : NEAR-PLATE RETENTION PREDICTOR VALUES (BARRIER
     > < -7 )')
       call prc('         '//Outer//' Plate  (S=0) ')
       call prc('         Ring         FF-factor       Fig-factor
     >Total  ')
       do 380 i = irsep,irlim
          write(coment,'(10X,i4,3(5x,g12.6))') i,
     >        kpredbar(idds(i,2),1,1),kpredbar(idds(i,2),2,1),
     >             kpredbar(idds(i,2),3,1)
          call prc(coment)
 380   continue
       call prc('                      (S=Sinj) ')
       do 381 i = irsep,irlim
          write(coment,'(10X,i4,3(5x,g12.6))') i,
     >        kpredbar(idds(i,2),1,2),kpredbar(idds(i,2),2,2),
     >             kpredbar(idds(i,2),3,2)
          call prc(coment)
 381   continue
       call prc('         '//Inner//' Plate  (S=0) ')
       call prc('         Ring         FF-factor       Fig-factor
     >Total  ')
       do 382 i = irsep,irlim
          write(coment,'(10X,i4,3(5x,g12.6))') i,
     >        kpredbar(idds(i,1),1,1),kpredbar(idds(i,1),2,1),
     >             kpredbar(idds(i,1),3,1)
          call prc(coment)
 382   continue
       call prc('                      (S=Sinj) ')
       do 383 i = irsep,irlim
          write(coment,'(10X,i4,3(5x,g12.6))') i,
     >        kpredbar(idds(i,1),1,2),kpredbar(idds(i,1),2,2),
     >             kpredbar(idds(i,1),3,2)
          call prc(coment)
 383   continue
c
c      endif for cprint = 5
c
       endif
c

      ENDIF

c
c     OVERRIDE Background Velocity Option
c

      call prb
      if (override_bg_velocity_opt.eq.0) then 
         call prc(sa//'FLOW VELOCITY OVERRIDE OPT 0: OFF')

      elseif (override_bg_velocity_opt.eq.1) then 
         call prc(sa//'FLOW VELOCITY OVERRIDE OPT 1:'//
     >                 ' PRESCRIBED FLOW')
c slmod begin
         IF (osmns28.GT.0) THEN
           CALL HD(7,
     >       '       HYDROGENIC FLOW PRESCRIPTION','HFLOWPRE-HD',
     >             8,67)
           status = .FALSE.
           DO i1 = 1, osmns28
             IF (osms28(i1,1).EQ.10.0) status = .TRUE.
           ENDDO
           IF (status) THEN
             WRITE(7,*)
             WRITE(7,'(7X,4A8)') 'Mode','x (m)','y (m)','Mach1'
             DO i1 = 1, osmns28
               IF (osms28(i1,1).EQ.10.0) THEN
                 WRITE(7,'(7X,7F8.2)') (osms28(i1,i2),i2=2,8)
               ENDIF
             ENDDO
           ELSE
             WRITE(7,*)
             WRITE(7,*) '       Not in use.'
           ENDIF
           WRITE(7,*)
         ENDIF
c slmod end

      elseif (override_bg_velocity_opt.eq.2) then 
         call prc(sa//'FLOW VELOCITY OVERRIDE OPT 2:'//
     >              ' RECALCULATE')
         call prc(sp//'FLOW IS RECALCULATED USING TARGET FLUXES')
         call prc(sp//'AND SOURCE DATA CALCULATED BY EIRENE')
         call prc(sp//'PLASMA DENSITY IS UNCHANGED')
      endif
      call prb 

c
C-----------------------------------------------------------------------
c
c     Iterative SOL Option
c
c slmod begin
      if (cpinopt.eq.1.or.cpinopt.eq.4) then  
c
c      if (cpinopt.eq.1) then  
c slmod end
c
       IF (CITERSOL.EQ.0) THEN
       CALL PRC ('  PIN ITERATION OPT 0: THE SOL IS NOT CALCULATED ITERA
     >TIVELY.')
c slmod begin
       ELSEIF (CITERSOL.EQ.1.OR.CITERSOL.EQ.2) THEN
c
c       ELSEIF (CITERSOL.EQ.1) THEN
c slmod end
       CALL PRC ('  PIN ITERATION OPT 1: THE SOL IS CALCULATED ITERATIVE
     >LY.')
       CALL PRI ('                       USING SOL OPT ',CSECSOL)
c
       CALL PRC ('                       WITH PIN DATA AVAILABLE')
c
       CALL PRC ('                       OTHER PRINTED BACKGROUND OPTION
     >S')
       CALL PRC ('                       APPLY TO THE FIRST ITERATION')
c slmod begin
       IF (CITERSOL.EQ.2) THEN
         CALL PRB
         CALL PRC ('    * SOLVER CALLED AFTER LAST CALL TO PIN *')
         CALL PRB
       ENDIF
c slmod end
c
c      Print out SOL option used for ITERATIONS if different from main
c      SOL option.
c
       if (csecsol.ne.cioptf) then

       call prc ('  SECONDARY SOL OPTION USED FOR PLASMA ITERATIONS: ')
c
c
      IF     (CSECSOL.EQ.-1) THEN
       WRITE (COMENT,'(''  SEC SOL OPTION  -1 : SOL1A,  (FL,FS) = ('',
     >   F6.3,'','',F6.3,'')'')') CFL,CFS
       CALL PRC (COMENT)
      ELSEIF (CSECSOL.EQ.0) THEN
       CALL PRC ('  SEC SOL OPTION   0 : SOL0,  VHOUT=0,  EOUT=0')
      ELSEIF (CSECSOL.EQ.1) THEN
       CALL PRC ('  SEC SOL OPTION   1 : SOL1')
      ELSEIF (CSECSOL.EQ.2) THEN
       CALL PRC ('  SEC SOL OPTION   2 : SOL2')
      ELSEIF (CSECSOL.EQ.3) THEN
       CALL PRC ('  SEC SOL OPTION   3 : SOL3')
      ELSEIF (CSECSOL.EQ.4) THEN
       CALL PRC ('  SEC SOL OPTION   4 : SOL4')
      ELSEIF (CSECSOL.EQ.5) THEN
       WRITE (COMENT,'(''  SEC SOL OPTION   5 : SOL5,  VHOUT='',
     >   1P,G11.4,'',  EOUT='',G11.4)') CVHOUT,CEOUT
       CALL PRC (COMENT)
      elseif (csecsol.eq.6) then
       CALL PRc ('  SEC SOL OPTION   6 : SOL6, FIXED E-field, Variable V
     >')
       call prr ('                       Eout = ',ceout)
       call prr ('                       Base Velocity Vhout (V0) = ',
     >           cvhout)
       write (coment,'(24x,''Fixed at V0 to '',g8.3,
     >'' * SMAX'')')   cvbl1
       call prc(coment)
       write (coment,'(24x,''Then fixed at '',g8.3,'' *V0 to '',g8.3,
     >'' * SMAX'')') cvbm1,cvbl2
       call prc(coment)
       write (coment,'(24x,''Then fixed at '',g8.3,'' *V0 to SMAX/2.0'')
     > ')       cvbm2
       call prc(coment)
      elseif (csecsol.eq.7) then
       CALL PRc ('  SEC SOL OPTION   7 : SOL7, FIXED E-field, Linear Vb
     >')
       call prr ('                       Eout = ',ceout)
       call prr ('                       Base Velocity Vhout (V0) = ',
     >           cvhout)
       write (coment,'(24x,''Falling to '',g8.3,'' *V0 at '',g8.3,
     >'' * SMAX'')') cvbm1,cvbl1
       call prc(coment)
       write (coment,'(24x,''Then to '',g8.3,'' *V0 at '',g8.3,
     >'' * SMAX'')') cvbm2,cvbl2
       call prc(coment)
       call prc ('                       Then to Vb = 0 at SMAX/2 ')
      ELSEIF (CSECSOL.EQ.9) THEN
       CALL PRC ('  SEC SOL OPTION   9 : SOL9')
      ELSEIF (CSECSOL.EQ.10) THEN
       CALL PRR ('  SEC SOL OPTION  10 : SOL10, REVERSAL MACH NUMBER',
     >    CFRM)
       WRITE (COMENT,'(''                       KIN'',F5.2,'',KOUT'',
     >   F5.2,'',FRMIN'',F6.3,'',FRMAX'',F6.3)')CKIN,CKOUT,CFRMIN,CFRMAX
       CALL PRC (COMENT)
       WRITE (7,'(24X,''(FL,FS) = ('',F6.3,'','',F6.3,'')'')') CFL,CFS
      ELSEIF (CSECSOL.EQ.12) THEN
       CALL PRC ('  SEC SOL OPTION  12 : SOL12 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/2K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(2KT+MV**2) = 4NOKTO')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CSECSOL.EQ.13) THEN
       CALL PRC ('  SEC SOL OPTION  13 : SOL13 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/2K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(KTI+KTE+MV**2) = 2N0K(T0i+T0
     >e)')
       CALL PRC ('                       P/A = (5KTEBP)*NBP*CSBP FOR EL
     >ECTRONS')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRC ('                       P/A = (2KTIBP)*NBP*CSBP FOR IO
     >NS')
       CALL PRC ('                       PR/A= 0 FOR IONS')
       CALL PRR ('                       K0I = ', CK0I)
      ELSEIF (CSECSOL.EQ.14) THEN
       CALL PRC ('  SEC SOL OPTION  14 : SOL14 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       HEAT CONDUCTION AND CONVECTION'
     >)
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       D/DS5NVKT-K0*T**(5/2)*(DT/DS)
     >=-PR(S)/A')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(2KT+MV**2) = 4NOKTO')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CSECSOL.EQ.15) THEN
       CALL PRC ('  SEC SOL OPTION  15 : SOL15 - PSEUDO SELF-CONSISTENT'
     >)
       CALL PRC ('                       OPTION OVERRIDES TGRAD OPTIONS
     >AND PLASMA')
       CALL PRC ('                       OPTION IN THE SOL.')
       CALL PRC ('                       TI,TE,NB,VH,E ARE CALCULATED')
       CALL PRC ('                       FROM THE FOLLOWING EQUATIONS')
       CALL PRC ('                       T=(TO**3.5+7/4K0*(P/A*S+INT2(P
     >R/A)))**(2/7)')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRC ('                       NV = NOVO + INT1(S(S)DS)')
       CALL PRC ('                       N(KTI+KTE+MV**2) = 2N0K(T0i+T0
     >e)')
       CALL PRC ('                       P/A = (5KTEBP)*NBP*CSBP FOR EL
     >ECTRONS')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRC ('                       P/A = (2KTIBP)*NBP*CSBP FOR IO
     >NS')
       CALL PRC ('                       PR/A= 0 FOR IONS')
       CALL PRR ('                       K0I = ', CK0I)
      elseif (csecsol.eq.16) then
       CALL PRC ('  SEC SOL OPTION  16 : SOL16 - AS SOL 14 - Linear Solu
     >tions')
      elseif (csecsol.eq.17) then
       CALL PRC ('  SEC SOL OPTION  17 : SOL17 - AS SOL 16 - Te and Ti')
       call prc ('                       are calculated independently.')
       CALL PRR ('                       K0  = ', CK0)
       CALL PRR ('                       K0I = ', CK0I)
      elseif (csecsol.eq.18) then
       CALL PRC ('  SEC SOL OPTION  18 : SOL18 - AS SOL 16 - except')
       call prc ('                       the 1/2mv**2*Gamma term has')
       call prc ('                       been added.')
      elseif (csecsol.eq.19) then
       CALL PRC ('  SEC SOL OPTION  19 : SOL19 - AS SOL 18 - except')
       call prc ('                       the hydrogenic rad term has')
       call prc ('                       been added.')
      elseif (csecsol.eq.20) then
       CALL PRC ('  SEC SOL OPTION  20 : SOL20 - AS SOL 17 - except')
       call prc ('                       the hydrogenic rad term and')
       call prc ('                       Pei terms have been added.')
      elseif (csecsol.eq.21) then
       CALL PRC ('  SEC SOL OPTION  21 : Detached plasma model')
       if (s21refsw.eq.0) then
        call prc('                       NOTE: REF OPTION 0 - ')
        call prc('                       All distances are specified')
        call prc('                       in units of SMAX. (MF=SMAX)')
       elseif (s21refsw.eq.1) then
        call prc('                       NOTE: REF OPTION 1 - ')
        call prc('                       All distances are specified')
        call prc('                       in units of PMAX. (MF=PMAX)')
        call prc('                       Converted to S for each'//
     >                                 ' ring.')
       elseif (s21refsw.eq.2) then
        call prc('                       NOTE: REF OPTION 2 - ')
        call prc('                       All distances are specified')
        call prc('                       in meters along S.(MF=1.0)')
       elseif (s21refsw.eq.3) then
        call prc('                       NOTE: REF OPTION 3 - ')
        call prc('                       All distances are specified')
        call prc('                       in meters along P.(MF=1.0)')
        call prc('                       Converted to S for each'//
     >                                 ' ring.')
       endif
c
c      OUTER
c
       call prc ('                       '//OUTER//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rat)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rat,l2rat
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rat)
       call prr ('                       Te at end of A: Te0 * ',
     >            terat)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirat)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrat)
c
       if (nrat.lt.0) then

        call prc('                       Nratio < 0 - variable N may'//
     >' be used to match pressure.')
        call prc('             IR          Nratio  for '//OUTER)
        do ir = irsep,nrs
           write(coment,'(12x,i5,10x,f10.2)') in,nrat_used(ir,2)
           call prc(coment)
        end do
        call prb
       endif
c
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalph)
       call prr ('                       Radiated power in B (Qrad): Q0
     >* ',qrat)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases lineraly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrat)
c
c      INNER
c
       call prb
       call prc ('                       '//INNER//' HALF OF SOL:')
c
       call prc ('                       Three regions :   ')
       call prr ('                       A :  0 < MF * ',l1rati)
       write (coment,'(''                       B :'',g8.3,
     >          '' * MF < '',g8.3, '' * MF '')') l1rati,l2rati
       call prc(coment)
       call prr ('                       C : S  >  MF * ',l2rati)
       call prr ('                       Te at end of A: Te0 * ',
     >            terati)
       call prr ('                       Ti at end of A: Ti0 * ',
     >            tirati)
       call prr ('                       Ne at end of A: Ne0 * '
     >           ,nrati)
c
       if (nrati.lt.0) then

        call prc('                       Nratio < 0 - variable N may'//
     >' be used to match pressure.')
        call prc('             IR          Nratio  for '//OUTER)
        do ir = irsep,nrs
           write(coment,'(12x,i5,10x,f10.2)') in,nrat_used(ir,1)
           call prc(coment)
        end do
        call prb
       endif
c
       call prc ('                       Ne(s) = Ne0 + (Ne1-Ne0)'//
     >                                    ' * (s/L)**ALPHA')
       call prr ('                       ALPHA = ',nalphi)
       call prr ('                       Radiated power in B (Qrad): Q0
     >* ',qrati)
       call prc ('                       Te increases linearly in A')
       call prc ('                       Ti increases linearly in A')
       call prc ('                       Ne increases linearly in A')
       call prc ('                       Velocity in A: v(s)=N0V0/n(s)')
       CALL PRC ('                       B: T=(T1**3.5+7/2K0*(Q0(s-L1)+'
     >)
       call prc ('                         1/2(s-L1)**2*(Qrad/Lrad)))**(
     >2/7)')
       CALL PRC ('                       C: T=(T2**3.5+7/2K0*(Qtot(s-L2)
     >))**(2/7)')
       call prc ('                          Qtot = Q0 + Qrad')
       call prc ('                       B,C: N(s) = N1 * (T(s)/T1)**(-1
     >)')
       call prr ('                       Veocity in B,C linearly->0 at M
     >F * ',lvrati)
c
      elseif (csecsol.eq.22) then
       CALL PRC ('  SEC SOL OPTION  22 : Runge Kutta SOL equation solver
     >')
       call echosol
      elseif (csecsol.eq.23) then
       CALL PRC ('  SEC SOL OPTION  23 : CFD ring by ring plasma solver'
     >)
       call echosol23
      ELSEIF (CSECSOL.EQ.98) THEN
       CALL PRC ('  SEC SOL OPTION  98 : FROM DIVIMP PLASMA INPUT FILE')
       call prc ('                       THIS ONLY READS IN THE ENTIRE')
       call prc ('                       BACKGROUND AND CAN ONLY BE SPECI
     >FIED')
       call prc ('                       IN COMBINATION WITH PLASMA DECAY
     >')
       call prc ('                       OPTION 98.')
      ELSEIF (CSECSOL.EQ.99) THEN
       CALL PRC ('  SEC SOL OPTION   99: FROM EQUILIBRIUM DATA')
      ENDIF

      IF (CSOLEF.NE.1.0.OR.CSOLVF.NE.1.0) THEN
       WRITE (COMENT,'(23X,2(A,F10.4))') 'E FACTOR',CSOLEF,
     >                                ',  V FACTOR',CSOLVF
       CALL PRC (COMENT)
      ENDIF

c
c     END OF CSECSOL.ne.CIOPTF
c
      call prb
c
      endif
c
c     Endif for Iteration Option (citersol.gt.0)
c
      ENDIF
c
c     Endif for PIN being run  (cpinopt.eq.1)
c       
      endif
c


c
C-----------------------------------------------------------------------
c
c     ELECTRON and ION Temperature gradient Options
c
      IF (CIOPTF.NE.12.AND.CIOPTF.NE.13.AND.CIOPTF.NE.14.and.
     >    cioptf.ne.15.and.cioptf.ne.16.and.cioptf.ne.17.and.
     >    cioptf.ne.18.and.cioptf.ne.19.and.cioptf.ne.20.and.
     >    cioptf.ne.21.and.cioptf.ne.22.and.cioptf.ne.23)
     >    THEN
      IF     (CIOPTK.EQ.0) THEN
       CALL PRC ('  TEB GRADIENT OPT 0 : LINEAR, FROM TEB0.FEBT AT TARGE
     >T TO TEB0')
       CALL PRC ('                       AT FEBL1.SMAX, THEN CONSTANT')
      ELSEIF (CIOPTK.EQ.1) THEN
       CALL PRC ('  TEB GRADIENT OPT 1 : LINEAR, FROM TEB0.FEBT AT TARGE
     >T, TO')
       CALL PRC ('                       TEB0.FEB2 AT FEBL1.SMAX, TO TEB
     >0')
       CALL PRC ('                       AT FEBL2.SMAX, THEN CONSTANT')
      ELSEIF (CIOPTK.EQ.2) THEN
       CALL PRR ('  TEB GRADIENT OPT 2 : P/A DRIVEN GRADIENTS. P/A =',
     >CPA)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.3) THEN
       CALL PRC ('  TEB GRADIENT OPT 3 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.4) THEN
       CALL PRC ('  TEB GRADIENT OPT 4 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.5) THEN
       CALL PRC ('  TEB GRADIENT OPT 5 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = 5KTEBP*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.6) THEN
       CALL PRC ('  TEB GRADIENT OPT 6 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.7) THEN
       CALL PRC ('  TEB GRADIENT OPT 7 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = 5KTEBP*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTK.EQ.98) THEN
       CALL PRC ('  TEB GRADIENT OPT 98: FROM DIVIMP PLASMA INPUT FILE')
       call prc ('                       THIS ONLY READS IN THE ENTIRE')
       call prc ('                       BACKGROUND AND CAN ONLY BE SPECI
     >FIED')
       call prc ('                       IN COMBINATION WITH PLASMA DECAY
     >')
       call prc ('                       OPTION 98.')
      ELSEIF (CIOPTK.EQ.99) THEN
       CALL PRC ('  TEB GRADIENT OPT 99: FROM EQUILIBRIUM DATA')
      ENDIF
C-----------------------------------------------------------------------
c
c     ION Temperature gradient Option
c
      IF     (CIOPTL.EQ.0) THEN
       CALL PRC ('  TIB GRADIENT OPT 0 : LINEAR, FROM TIB0.FIBT AT TARGE
     >T TO TIB0')
       CALL PRC ('                       AT FIBL1.SMAX, THEN CONSTANT')
      ELSEIF (CIOPTL.EQ.1) THEN
       CALL PRC ('  TIB GRADIENT OPT 1 : LINEAR, FROM TIB0.FIBT AT TARGE
     >T, TO')
       CALL PRC ('                       TIB0.FIB2 AT FIBL1.SMAX, TO TIB
     >0')
       CALL PRC ('                       AT FIBL2.SMAX, THEN CONSTANT')
      ELSEIF (CIOPTL.EQ.2) THEN
       CALL PRR ('  TIB GRADIENT OPT 2 : P/A DRIVEN GRADIENTS. P/A =',
     >CPA)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTL.EQ.3) THEN
       CALL PRC ('  TIB GRADIENT OPT 3 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTL.EQ.4) THEN
       CALL PRC ('  TIB GRADIENT OPT 4 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTL.EQ.5) THEN
       CALL PRC ('  TIB GRADIENT OPT 5 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = 2KTIBP*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/4 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0I)
      ELSEIF (CIOPTL.EQ.6) THEN
       CALL PRC ('  TIB GRADIENT OPT 6 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = (2KTIBP+5KTEBP)*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0)
      ELSEIF (CIOPTL.EQ.7) THEN
       CALL PRC ('  TIB GRADIENT OPT 7 : P/A DRIVEN GRADIENTS. P/A '//
     > 'CALCULATED')
       CALL PRC ('                       BASED ON INPUT DATA/RING')
       CALL PRC ('                       P/A = 2KTIBP*NBP*CSBP'
     >)
       CALL PRC ('                       DATA FOR INNER AND OUTER PLATES
     >')
       CALL PRC ('                       7/2 FACTOR ASSUMED')
       CALL PRR ('                       K0 = ', CK0I)
      ELSEIF (CIOPTL.EQ.98) THEN
       CALL PRC ('  TEB GRADIENT OPT 98: FROM DIVIMP PLASMA INPUT FILE')
       call prc ('                       THIS ONLY READS IN THE ENTIRE')
       call prc ('                       BACKGROUND AND CAN ONLY BE SPECI
     >FIED')
       call prc ('                       IN COMBINATION WITH PLASMA DECAY
     >')
       call prc ('                       OPTION 98.')
      ELSEIF (CIOPTL.EQ.99) THEN
       CALL PRC ('  TIB GRADIENT OPT 99: FROM EQUILIBRIUM DATA')
      ENDIF

C-----------------------------------------------------------------------
c
c     Density gradient Option
c
      IF     (ngradopt.EQ.0) THEN
       CALL PRC ('  DENSITY GRAD OPT 0 : CONSTANT NE - NO GRADIENT APPLI
     >ED')
      ELSEIF (ngradopt.EQ.1) THEN
       CALL PRC ('  DENSITY GRAD OPT 1 : DENSITY IS CALCULATED FROM THE
     >FOLLOWING:')
       call prc ('                       N(S) = 4 Nt * TEt / (Te(S)+Ti(S
     >))')
       call prc ('                       Tet = Tit IS ASSUMED') 
      endif
c
c     Endif for SOL (cioptf) option if-statement
c
      ENDIF
C-----------------------------------------------------------------------
      IF     (CIOPTM.EQ.0) THEN
       CALL PRC ('  TEB GRAD COEFF   0 : ALPHAE = 0')
      ELSEIF (CIOPTM.EQ.1) THEN
       CALL PRC ('  TEB GRAD COEFF   1 : ALPHAE = 0.71 ZI.ZI')
c       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,min(6,NIZS))
      ELSEIF (CIOPTM.EQ.2) THEN
       CALL PRC ('  TEB GRAD COEFF   2 : ALPHAE = 1.5(1-0.6934(1.3167**-
     >ZI))ZI.ZI')
c       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,min(6,nizs))
      ELSEIF (CIOPTM.EQ.3) THEN
       CALL PRC ('  TEB GRAD COEFF   3 : ALPHAE = 0.71 ZI.ZI')
c       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KALPHS(IZ),IZ=1,min(6,nizs))
       call prr ('                       FEG -> 0 for S(from target)'//
     >           ' > SMAX*',cstgrad)
      ENDIF
C-----------------------------------------------------------------------
      IF     (CIOPTN.EQ.0) THEN
       CALL PRC ('  TIB GRAD COEFF   0 : BETAI = 0')
      ELSEIF (CIOPTN.EQ.1) THEN
       CALL PRC ('  TIB GRAD COEFF   1 : BETAI=-3(1-MU-5.ZI.ZI.SQRT(2MU)
     >MU(1.1MU-')
       CALL PRC ('                       0.35))/(2.6-2MU+5.4MU.MU), MU=M
     >I/(MI+MB)')
c       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,min(6,NIZS))
      ELSEIF (CIOPTN.EQ.2) THEN
       CALL PRC ('  TIB GRAD COEFF   2 : BETAI=H(ZO)ZI.ZI/(ZO+SQRT(0.5(1
     >+MB/MI)))')
       CALL PRC ('                       WHERE H(ZO)=1.56(1+1.41ZO)(1+0.
     >52ZO)/')
       CALL PRC ('                                   ((1+2.65ZO)(1+0.285
     >ZO))')
       WRITE (7,'(23X,A,I3,A,F6.2)') 'ZO=', CZO, ',  H(ZO)=', CHZO
c       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,min(6,NIZS))
      ELSEIF (CIOPTN.EQ.3) THEN
       CALL PRC ('  TIB GRAD COEFF   3 : BETAI=-3(1-MU-5.ZI.ZI.SQRT(2MU)
     >MU(1.1MU-')
       CALL PRC ('                       0.35))/(2.6-2MU+5.4MU.MU), MU=M
     >I/(MI+MB)')
c       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,NIZS)
       WRITE (7,'(23X,6F6.2)') (KBETAS(IZ),IZ=1,min(6,NIZS))
       call prr ('                       FIG -> 0 for S(from target)'//
     >           ' > SMAX *',cstgrad)
      ENDIF
C-----------------------------------------------------------------------
      if (fgradopt.eq.0) then
       call prc ('  F-GRAD MOD OPTION 0: OFF - The temperature gradient
     >forces')
       call prc ('                       are NOT corrected for kinetic e
     >ffects.')
       call prc ('                       The force modifier is set equal
     > to')
       call prc ('                       1.0 everywhere.')
      elseif (fgradopt.eq.1) then
       call prc ('  F-GRAD MOD OPTION 1: ON - The temperature gradient f
     >orces')
       call prc ('                       are corrected for kinetic effec
     >ts')
       call prc ('                       using an ad-hoc formula adapted
     > from')
       call prc ('                       the UEDGE code. The force in ea
     >ch cell')
       call prc ('                       is modified by the following fa
     >ctor')
       call prc ('                       which is calculated for each ce
     >ll.')
       call prc ('                       Cuedge = 1 / (1+Fact*(Lmfp/Lgra
     >d)**2)')
       call prr ('                       Where: Fact = ',fgradfact)
       call prc ('                       Lmfp = sum of mean free paths')
       call prc ('                       Lgrad = minimum of the scale le
     >ngths')
       call prc ('                               for Te,Ti, n and Pressu
     >re')
      elseif (fgradopt.eq.2) then
       call prc ('  F-GRAD MOD OPTION 2: ON - The temperature gradient f
     >orces')
       call prc ('                       are corrected for kinetic effec
     >ts')
       call prc ('                       using an ad-hoc formula adapted
     > from')
       call prc ('                       the UEDGE code. The force in ea
     >ch cell')
       call prc ('                       is modified by the following fa
     >ctor')
       call prc ('                       which is calculated for each ce
     >ll.')
       call prc ('                       Cuedge = 1 / (1+Fact*(Lmfp/Lgra
     >d)**2)')
       call prr ('                       Where: Fact = ',fgradfact)
       call prc ('                       Lmfp = sum of mean free paths')
       call prc ('                       Lgrad = minimum of the scale le
     >ngths')
       call prc ('                               for Te,Ti, n and Pressu
     >re')
       call prc ('                       FiG -> 0.0 within one ion mean'
     >)
       call prc ('                       free path of the target')
      elseif (fgradopt.eq.3) then
       call prc ('  F-GRAD MOD OPTION 3: ON - The temperature gradient f
     >orces')
       call prc ('                       are corrected for kinetic effec
     >ts')
       call prc ('                       using an ad-hoc formula adapted
     > from')
       call prc ('                       the GARCHING-B2 code. The force
     > in each')
       call prc ('                       cell is modified by the followi
     >ng')
       call prc ('                       which is calculated for each ce
     >ll.')
       call prc ('                       T-GRAD = min(T-GRAD(norm),0.3 T
     >/Lam')
       call prc ('                       Lam =  Fact * 1.5e16 * T**2/n')
       call prr ('                       Where: Fact = ',fgradfact)
      elseif (fgradopt.eq.4) then
       call prc ('  F-GRAD MOD OPTION 4: ON - The temperature gradient f
     >orces')
       call prc ('                       are corrected for kinetic effec
     >ts')
       call prc ('                       using an ad-hoc formula adapted
     > from')
       call prc ('                       the UEDGE code. The force in ea
     >ch cell')
       call prc ('                       is modified by the following fa
     >ctor')
       call prc ('                       which is calculated for each ce
     >ll.')
       call prc ('                       Cuedge = 1 / (1+Fact*(Lmfp/Lgra
     >d)**2)')
       call prr ('                       Where: Fact = ',fgradfact)
       call prc ('                       Lmfp = sum of mean free paths')
       call prc ('                       Lgrad = minimum of the scale le
     >ngths')
       call prc ('                               for Te,Ti, n and Pressu
     >re')
       call prc ('                       FiG -> 0.0 within one ion mean'
     >)
       call prc ('                       free path of the target')
       call prc ('                       FeG -> 0.0 within one electron
     >mean')
       call prc ('                       free path of the target')
      endif


C-----------------------------------------------------------------------
      if (ofield.eq.0) then
       CALL PRC ('  E-FIELD OPTION   0 : ELECTRIC FIELD IS NOT OVER-WRIT
     >TEN')
      elseif (ofield.eq.1) then
       CALL PRC ('  E-FIELD OPTION   1 : ELECTRIC FIELD SET TO ZERO')
      elseif (ofield.eq.2) then
       CALL PRC ('  E-FIELD OPTION   2 : ELECTRIC FIELD LIMITED TO CONST
     >ANT')
       CALL PRC ('                       INSIDE FIRST INTERVAL FOR SOL O
     >PTS 6 & 7,')
       call prc ('                       ZERO ELSEWHERE')
      elseif (ofield.eq.3) then
       CALL PRC ('  E-FIELD OPTION   3 : ELECTRIC FIELD CALCULATED USING
     >')
       call prc ('                       THE DENSITY AND TEMPERATURES AN
     >D')
       call prc ('                       A SIMPLE FORMULA BASED ON PRESS
     >URE')
       call prc ('                       AND TEMPERATURE GRADIENTS.')
      elseif (ofield.eq.4) then
       CALL PRC ('  E-FIELD OPTION   4 : ELECTRIC FIELD CALCULATED USING
     >')
       call prc ('                       THE DENSITY AND TEMPERATURES AN
     >D')
       call prc ('                       A SIMPLE FORMULA FOR RINGS')
       call prc ('                       WHICH ARE DEEMED COLLISONAL.')
       call prc ('                       A HALF-RING IS DEEMED NON-COLLI
     >SIONAL')
       call prc ('                       IF Te-mid < FACT * Te-target')
       call prr ('                       WHERE FACT = ',ceffact)
       call prc ('                       THE EFIELD IS THEN CALULATED TO
     > BE')
       call prc ('                       Ef = -(Te-average)/2*Lsource')
       call prc ('                       Te-average= 0.5 * (Te-mid+Te-ta
     >rget')
c
       if (cpinopt.eq.0) then
          call prr ('                       Lsource = SMAX * ', ceflen)
c slmod begin
       elseif (cpinopt.eq.1.or.cpinopt.eq.4) then
c
c       elseif (cpinopt.eq.1) then
c slmod end
          call prc ('                       Lsource = Lequiv from PIN')
       endif
c
      endif
c
      if (ofield.eq.1.or.ofield.eq.2.or.ofield.eq.3.or.ofield.eq.4) then
         if (ciopto.eq.0.or.ciopto.eq.2.or.
     >       ciopto.eq.3.or.ciopto.eq.4) then

       CALL PRC ('                       ELECTRIC FIELD IS NOT OVERWRITT
     >EN')
       call prc ('                       IN THE PRIVATE PLASMA')

         elseif (ciopto.eq.1) then

       CALL PRC ('                       ELECTRIC FIELD OPTION ALSO APPL
     >IES')
       call prc ('                       IN THE PRIVATE PLASMA')

         endif

      endif

C-----------------------------------------------------------------------
      if (ccoreopt.eq.-1) then
       CALL PRC ('  CORE DECAY OPT -1  : CORE SPECIFICIATION IS OFF')
       call prc ('                       OPTION IS USED IN COMBINATION')
       call prc ('                       WITH OTHER OPTIONS THAT SPECIFY
     >')
       call prc ('                       THE ENTIRE PLASMA.')
      elseif (ccoreopt.eq.0) then
       CALL PRC ('  CORE DECAY OPT 0   : STANDARD - PROFILES FLAT IN COR
     >E')
       call prc ('                       Specified by quantities above')

      elseif (ccoreopt.eq.1) then
       CALL PRC ('  CORE DECAY OPT 1   : CORE QUANTITIES ARE SPECIFIED I
     >N THE INPUT')
       call prc ('                       ARRAY. CORE PROFILES ARE FLAT.'
     >)
       call prc ('                       DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') (coredat(in,i),i=1,5)
          call prc(coment)
       end do
c
      elseif (ccoreopt.eq.2) then
c
       CALL PRC ('  CORE DECAY OPT 2   : MARFE PRESCRIPTION - CORE PROFI
     >LES')
       call prc ('                       ARE FLAT FOR THE HALF OF RINGS'
     >)
       CALL PRC ('                       NEAREST THE MID-PLANE - FALLING
     >')
       CALL PRC ('                       TO VALUES ENTERED FOR THE CORE'
     >)
       CALL PRC ('                       PLASMA DATA AT THE X-POINT')
       CALL PRC ('                       VALUES IN THE TOP HALF-RING ARE
     >')
       CALL PRC ('                       CALCULATED USING STANDARD METHO
     >DS')
       call prc ('                       X-POINT DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') (coredat(in,i),i=1,5)
          call prc(coment)
       end do
c
c      Core option 3
c
      elseif (ccoreopt.eq.3) then
c
       CALL PRC ('  CORE DECAY OPT 3   : MARFE PRESCRIPTION - CORE PROFI
     >LES')
       call prc ('                       ARE FLAT FOR THE HALF OF RINGS'
     >)
       CALL PRC ('                       NEAREST THE MID-PLANE - FALLING
     >')
       CALL PRC ('                       TO VALUES ENTERED FOR THE CORE'
     >)
       CALL PRC ('                       PLASMA DATA AT THE X-POINT')
       CALL PRC ('                       VALUES IN THE TOP HALF-RING ARE
     >')
       CALL PRC ('                       CALCULATED USING STANDARD METHO
     >DS')
       call prc ('                       DENSITY IS NOT USED ON SPECIFIC
     >ATION ')
       call prc ('                       IT IS CALCULATED BASED ON PRESS
     >URE')
       call prc ('                       CONSERVATION.')
       call prc ('                       X-POINT DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') coredat(in,1),coredat(in,2),
     >           coredat(in,3),knbs(1,int(coredat(in,1))),coredat(in,5)
          call prc(coment)
       end do
c
c
      elseif (ccoreopt.eq.4) then
c
c
       CALL PRC ('  CORE DECAY OPT 4   : MARFE PRESCRIPTION - CORE PROFI
     >LES')
       call prc ('                       ARE FLAT FOR SPECIFIED LENGTHS'
     >)
       CALL PRC ('                       FALLING TO VALUES ENTERED'//
     >                                         ' FOR THE CORE')
       CALL PRC ('                       PLASMA DATA AT THE X-POINT')
       CALL PRC ('                       VALUES IN THE FLAT PROFILE ARE'
     >                                   //'A ARE')
       CALL PRC ('                       CALCULATED USING STANDARD METHO
     >DS')
       call prc ('                       DENSITY IS NOT USED ON SPECIFIC
     >ATION ')
       call prc ('                       IT IS CALCULATED BASED ON PRESS
     >URE')
       call prc ('                       CONSERVATION.')
c
       call prr ('                       VELOCITY -> 0     FOR Sv < SMA'
     >                                   //'X * ',corefv)
       call prr ('                       VELOCITY -> 0     FOR Sv > SMA'
     >                                   //'X * ',corefv2)
       call prc ('                       VELOCITY -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
       call prr ('                       Te,Ti -> XPT T FOR St < SMAX *'
     >                                   ,coreft)
       call prr ('                       Te,Ti -> MID T FOR St > SMAX *'
     >                                   ,coreft2)
       call prc ('                       Te,Ti -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
c
       call prc ('                       X-POINT DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') coredat(in,1),coredat(in,2),
     >           coredat(in,3),knbs(1,int(coredat(in,1))),coredat(in,5)
          call prc(coment)
       end do
c
c
      elseif (ccoreopt.eq.5) then
c
c
       CALL PRC ('  CORE DECAY OPT 5   : MARFE PRESCRIPTION - CORE PROFI
     >LES')
       call prc ('                       ARE FLAT FOR SPECIFIED LENGTHS'
     >)
       CALL PRC ('                       FALLING TO VALUES ENTERED'//
     >                                        ' FOR THE CORE')
       CALL PRC ('                       PLASMA DATA AT THE X-POINT')
       CALL PRC ('                       VALUES IN THE FLAT PROFILE ARE'
     >                                   //'A ARE')
       call prc ('                       READ FROM INPUT DATA FOR EACH R
     >ING')
       call prc ('                       DENSITY IS NOT USED ON SPECIFIC
     >ATION ')
       call prc ('                       IT IS CALCULATED BASED ON PRESS
     >URE')
       call prc ('                       CONSERVATION.')
c
       call prr ('                       VELOCITY -> 0     FOR Sv < SMA'
     >                                   //'X * ',corefv)
       call prr ('                       VELOCITY -> 0     FOR Sv > SMA'
     >                                   //'X * ',corefv2)
       call prc ('                       VELOCITY -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
       call prr ('                       Te,Ti -> XPT T FOR St < SMAX *'
     >                                   ,coreft)
       call prr ('                       Te,Ti -> MID T FOR St > SMAX *'
     >                                   ,coreft2)
       call prc ('                       Te,Ti -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
c
       call prc ('                       X-POINT DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') coredat(in,1),coredat(in,2),
     >           coredat(in,3),knbs(1,int(coredat(in,1))),coredat(in,5)
          call prc(coment)
       end do
c
       call prc ('                       MID-POINT DATA FOR CORE RINGS:'
     >)
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do ir = 1,irsep-1
          write (coment,'(10x,i10,4g10.3)') ir,ktebs(nks(ir)/2,ir),
     >            ktibs(nks(ir)/2,ir),knbs(nks(ir)/2,ir),
     >            kvhs(nks(ir)/2,ir)
          call prc(coment)
       end do
c
c
      elseif (ccoreopt.eq.6) then
c
c
       CALL PRC ('  CORE DECAY OPT 6   : MARFE PRESCRIPTION - CORE PROFI
     >LES')
       call prc ('                       ARE FLAT FOR SPECIFIED LENGTHS'
     >)
       CALL PRC ('                       FALLING TO VALUES ENTERED'//
     >                                         ' FOR THE CORE')
       CALL PRC ('                       PLASMA DATA AT THE X-POINT')
       CALL PRC ('                       VALUES IN THE FLAT PROFILE ARE'
     >                                   //'A ARE')
       CALL PRC ('                       CALCULATED USING STANDARD METHO
     >DS')
       call prc ('                       DENSITY IS RAMPED TO THE'//
     >                               ' INPUT X-POINT VALUES USING')
       call prc ('                       THE SAME DISTANCE'//
     >                               ' SPECIFICATIONS AS TEMPERTATURE')

c
       call prr ('                       VELOCITY -> 0     FOR Sv < SMA'
     >                                   //'X * ',corefv)
       call prr ('                       VELOCITY -> 0     FOR Sv > SMA'
     >                                   //'X * ',corefv2)
       call prc ('                       VELOCITY -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
       call prr ('                       Te,Ti -> XPT T FOR St < SMAX *'
     >                                   ,coreft)
       call prr ('                       Te,Ti -> MID T FOR St > SMAX *'
     >                                   ,coreft2)
       call prc ('                       Te,Ti -> LINEAR RAMP'//
     >                                                ' IN BETWEEN')
c
       call prc ('                       X-POINT DATA FOR CORE RINGS:')
       write(coment,
     >     '(6x,6x,''RING'',8x,''Te'',8x,''Ti'',8x,''Nb'',8x,''Vb'')')
       call prc(coment)
       do in = 1,ncoredat
          write (coment,'(10x,5g10.3)') coredat(in,1),coredat(in,2),
     >           coredat(in,3),knbs(1,int(coredat(in,1))),coredat(in,5)
          call prc(coment)
       end do
c
c
      ENDIF
c
C-----------------------------------------------------------------------
      IF     (CIOPTO.EQ.0.and.cioptg.ne.99) THEN
       CALL PRC ('  PRIVATE PLASMA BG 0: SOL BACKGROUND PLASMA OPTIONS O
     >FF IN TRAP')
       call prc ('                       PRIVATE PLASMA CONDITIONS ARE C
     >ONSTANT')
       call prc ('                       AT THE TARGET VALUES')
      ELSEIF (CIOPTO.EQ.1.and.cioptg.ne.99) THEN
       CALL PRC ('  PRIVATE PLASMA BG 1: SOL BACKGROUND PLASMA OPTIONS O
     >N IN TRAP')
       call prc ('                       THE PRIVATE PLASMA IS CALCULATE
     >D')
       CALL PRC ('                       USING THE SPECIFIED SOL OPTION.
     >')
      ELSEIF (CIOPTO.EQ.2) THEN
       CALL PRC ('  PRIVATE PLASMA BG 2: SPECIFIED PLASMA CONDITIONS')
       call prc ('                       IN THE PRIVATE PLASMA REGION')
       call prb
       call prc('      All distances are proportions of SMAX'//
     >          ' for the ring')
c
       call prc('  Quantities are linearly interpolated between points')
       call prc('  - constant outside range at end-point values')
       write(coment,1161)  ctes1,ctef1,ctes2,ctef2
       call prc(coment)
       write(coment,1171)  ctis1,ctif1,ctis2,ctif2
       call prc(coment)
       write(coment,1181)  cnes1,cnef1,cnes2,cnef2
       call prc(coment)
       write(coment,1191)  cvbs1,cvbf1,cvbs2,cvbf2
       call prc(coment)
c
 1161  format('  (S,Te) : (0.0,Tet) -> (',f6.3,',',f6.3,'*Tet) -> (',
     >            f6.3,',',f6.3,'*Tet)')
 1171  format('  (S,Ti) : (0.0,Tit) -> (',f6.3,',',f6.3,'*Tit) -> (',
     >            f6.3,',',f6.3,'*Tit)')
 1181  format('  (S,Ne) : (0.0,Net) -> (',f6.3,',',f6.3,'*Net) -> (',
     >            f6.3,',',f6.3,'*Net)')
 1191  format('  (S,Vb) : (0.0,Vbt) -> (',f6.3,',',f6.3,'*Vbt) -> (',
     >            f6.3,',',f6.3,'*Vbt)')
      ELSEIF (CIOPTO.EQ.3) THEN
       CALL PRC (sa//'PRIVATE PLASMA BG 2: PLASMA CONDITIONS')
       call prc (sp//'IN THE PRIVATE PLASMA REGION')
       call prc (sp//'ARE CALCULATED FROM THOMSON MEASUREMENTS.') 
       call prc (sp//'ALL DATA FOR EACH FLUX TUBE ARE AVERAGED.') 
       call prc (sp//'TARGET CONDITIONS ARE TAKEN FROM INPUT') 
       call prc (sp//'LANGMUIR PROBE VALUES.') 
       call prb
      ELSEIF (CIOPTO.EQ.4) THEN
       CALL PRC (sa//'PRIVATE PLASMA BG 2: PLASMA CONDITIONS')
       call prc (sp//'IN THE PRIVATE PLASMA REGION')
       call prc (sp//'ARE CALCULATED FROM THOMSON MEASUREMENTS.') 
       call prc (sp//'ALL DATA FOR EACH FLUX TUBE ARE AVERAGED.') 
       call prc (sp//'TARGET CONDITIONS ARE ASSIGNED TO EQUAL THE') 
       call prc (sp//'FLUX TUBE VALUES.') 
       call prb
      ENDIF


c
C-----------------------------------------------------------------------
c
c     If full print add this summary
c
      if (cprint.eq.1.or.cprint.eq.9) then
c
         call prb
         call prc('  SUMMARY OF TARGET AND UPSTREAM CONDITIONS FOR'//
     >             ' THE SELECTED SOL OPTION:')
c
         write(coment,1030)
         call prc(coment)
c
         do ir = irsep,irwall
c
c           Adjust print-out for INNER/OUTER
c
            write(coment,1040) ir,inner,kteds(idds(ir,1)),
     >            ktids(idds(ir,1)),knds(idds(ir,1)),teupstream(ir,1),
     >            tiupstream(ir,1),nbupstream(ir,1)
            call prc(coment)
            write(coment,1040) ir,outer,kteds(idds(ir,2)),
     >            ktids(idds(ir,2)),knds(idds(ir,2)),teupstream(ir,2),
     >            tiupstream(ir,2),nbupstream(ir,2)
            call prc(coment)
c
         end do
c
         call prb
c
      endif

1030  format(2x,'Ring',12x,'Te targ',2x,'Ti targ',3x,'Nb Targ',3x,
     >      'Te mid',2x,'Ti mid',4x,'Nb Mid')
1040  format(1x,i4,3x,a5,3x,2(1x,f8.3),1x,g9.2,2(1x,f8.3),1x,
     >       1P,g9.2)
c
c1050  format(1x,i4,3x,'OUTER',3x,2(1x,f8.3),1x,g9.2,2(1x,f8.3),1x,
c     >       1P,g9.2)
c
c
c
      RETURN
      END
c
c
c
      SUBROUTINE PR_pin_options
      IMPLICIT none
c
C
C***********************************************************************
c
C     THIS ROUTINE PRINTS ALL THE PIN OPTION FLAGS
C
C***********************************************************************
C
      include    'params'
      include    'comtor'
c      include    'cadas'
      include    'cgeom'
c      include    'cioniz'
c      include    'cedge2d'
c      include    'dynam4'
c      include    'dynam5'
      include    'pindata'
c      include    'adpak_com'
c      include    'promptdep'
      include    'slcom'  
      include    'printopt'
c
c      CHARACTER  COMENT*80,prtype*4
      INTEGER  I,IT,IZ,irlim,ir,in,len,lenstr,id,ind
c      logical  prsol21,prsol22
      real     totfpin,totfypin,tmpy,totzfpin,tothpin,totfapin
      real     totmhpin
      external lenstr


C-----------------------------------------------------------------------
      CALL PRB
c      CALL PRC ('--- PIN OPTIONS --- ')
      CALL PRChtml ('--- PIN OPTIONS ---','pr_pin_options','0','B')
      CALL PRB
C-----------------------------------------------------------------------
c
c     PIN CODE OPTIONS
c
c
      IF (CPINOPT.EQ.0) THEN
       CALL PRC ('  PIN RUN OPT      0 : PIN NOT EXECUTED FOR BG')
       call prc ('                       PLASMA SOLUTION')

      ELSEIF (CPINOPT.EQ.1.OR.CPINOPT.EQ.4) THEN
       IF (CPINOPT.EQ.1) THEN
         CALL PRC ('  PIN RUN OPT      1 : PIN IS EXECUTED. SOME '// 
     >             'RESULTS RETAINED.')
       ELSE
         CALL PRC ('  PIN RUN OPT      4 : PIN IS EXECUTED. SOME '// 
     >             'RESULTS RETAINED.  PIN IMPURITY SOURCE USED.')
       ENDIF
       CALL PRC ('                       THE PIN INVOCATION '//
     >           'COMMAND IS:')
       LEN = LENSTR(actpin)
       CALL PRC ('                      '//ACTPIN(1:LEN))
c
c
c      Which PIN code is running?
c
       if (pincode.eq.0) then
          call prb
          call prc ('NIMBUS called for hydrogenic data')
          call prnimbin(cnimbin,nlines)
       elseif (pincode.eq.1) then
          call prb
          call prc ('EIRENE called for hydrogenic data')
c slmod begin
c...      Need designated stream for .dat file:
          call slopt02(7)
c slmod end
       elseif (pincode.eq.2) then
          call prb
          call prc ('EIRENE-99 called for hydrogenic data')
c slmod begin
c...      Need designated stream for .dat file:
          call slopt02(7)
       elseif (pincode.eq.3) then
          call prb
          call prc ('EIRENE-02 called for hydrogenic data')
          call slopt02(7)
       elseif (pincode.eq.4) then
          call prb
          call prc ('EIRENE-04 called for hydrogenic data')
          call slopt02(7)
       elseif (pincode.eq.5) then
          call prb
          call prc ('EIRENE-06 called for hydrogenic data')
          call slopt02(7)
c slmod end
       endif
c
       call prb
c
       if (pinprint.eq.0) then
        CALL PRC ('  PIN PRINT OPTION 0 : STANDARD PIN DATA'//
     >            ' PRINT OPTION')
       elseif (pinprint.eq.1) then
        CALL PRC ('  PIN PRINT OPTION 1 : EXTENDED PIN DATA'//
     >                                    ' PRINT OPTION')
        CALL PRC ('                       - Includes STANDARD')
        CALL PRC ('                       - ITEMIZES IONIZATION AND'//
     >            ' NEUTRAL CONTENT BY RING')
        CALL PRC ('                       - ITEMIZES AVERAGE'//
     >            ' IONIZATION AND NEUTRAL DENSITIES BY RING')
       endif
c
       if (pincode.eq.0.or.pincode.eq.2.or.pincode.eq.3.or.
     >     pincode.eq.4.or.pincode.eq.5) then
c
c       Print out information about the PIN wall-target source.
c
        call prb
        call prc('   PIN Source Characteristics:')
c
c        CALL PRC('   SAMPLE PRIMARY FLUX AND YIELD DATA')
c
        call prb
c
        totfpin = 0.0
        totfypin= 0.0
        totzfpin= 0.0
        tothpin = 0.0
        totmhpin= 0.0
        totfapin= 0.0
c
c        if (wallpts.gt.0) write(7,9001)
c
        DO ID = 1, wallpts
c
           IN = WALLPT(id,17)
c
c           if (flxhw2(in).le.0.0) then
c              tmpy = 0.0
c           else
c              tmpy = flxhw3(in)/flxhw2(in)
c           endif
c
c           WRITE (7,9003) id,
c     >         wallpt(ID,1),wallpt(id,2),fluxhw(in),flxhw2(in),
c     >         flxhw5(in),flxhw3(in),tmpy,
c     >         flxhw4(in),wallpt(id,16)
c
           totfpin = totfpin + flxhw2(in) * wallpt(id,7)
           totfypin = totfypin + flxhw3(in) * wallpt(id,7)
           totzfpin = totzfpin + flxhw4(in) * wallpt(id,7)
           totfapin = totfapin + flxhw6(in) * wallpt(id,7)
           tothpin = tothpin + flxhw6(in) * flxhw5(in) * wallpt(id,7)
           totmhpin =totmhpin + (fluxhw(in)-flxhw6(in))* kboltz
     >                          * wallpt(id,19) * wallpt(id,7)
c
        end do
c
         call prb
c
         call prc('   CALCULATED FROM INDIVIDUAL SEGMENT DATA:')
         CALL PRR('     TOTAL PRIMARY INTEGRATED ATOM+ION FLUX',
     >                totfpin)
         CALL PRR('     TOTAL PRIMARY INTEGRATED FLUX*YIELD   ',
     >                totfypin)
         CALL PRR('     TOTAL INTEGRATED Z-REDEPOSITION       ',
     >                totzfpin)
         call prr('     TOTAL ATOMIC H HEAT FLUX TO SURFACES  ',
     >                tothpin*ech)
         call prr('     ATOM HEAT+REC POWER FLUX TO SURFACES  ',
     >                tothpin*ech+2.2*ech*totfapin)
         call prr('     MOL H HEAT FLUX TO SURFACES (E=SURF T)',
     >                totmhpin)
c
         call prb
c
c
       endif
c
c      HYBRID WALL OPTION - only for JET
c
       if (ihybrid.eq.0) then
        call prc('  PIN HYBRID WALL  0 : STANDARD WALL FROM THE EQUILIBR
     >IUM')
        call prc('                       FILE IS USED IN PIN')
       elseif (ihybrid.eq.1) then
        call prc('  PIN HYBRID WALL  1 : MARK I HYBRID WALL FILE')
        call prc('                       IS USED IN PIN')
       elseif (ihybrid.eq.2) then
        call prc('  PIN HYBRID WALL  2 : MARK IIA HYBRID WALL FILE')
        call prc('                       IS USED IN PIN')
       elseif (ihybrid.eq.3) then
        call prc('  PIN HYBRID WALL  3 : CORRECTED MARK IIA HYBRID'//
     >                      ' WALL FILE')
        call prc('                       IS USED IN PIN')
       elseif (ihybrid.eq.4) then
        call prc('  PIN HYBRID WALL  4 : MODIFIED CORRECTED MARK IIA'//
     >                      ' HYBRID'//
     >                      ' WALL FILE')
        call prc('                       IS USED IN PIN')
       endif
c
C-----------------------------------------------------------------------
c
c      Check for options applicable only to JET/NIMBUS
c
       if (pincode.eq.0) then
c
C-----------------------------------------------------------------------
c
c     Recombination Calculation Option
c
      if (crecopt.eq.0) then
       CALL PRC ('  RECOMBINATION OPT 0: HYDROGENIC RECOMBINATION IS OFF
     >')
c
      elseif (crecopt.eq.1) then
c
       CALL PRC ('  RECOMBINATION OPT 1: HYDROGENIC RECOMBINATION IS CAL
     >CULATED')
       call prc ('                       USING GORDEEV COEFFICIENTS')
       call prc ('                       WITH AN IMPOSED LOWER TEMPERATU
     >RE')
       call prr ('                       LIMIT OF (eV) :',treccut)
c
      elseif (crecopt.eq.2) then
c
       CALL PRC ('  RECOMBINATION OPT 2: HYDROGENIC RECOMBINATION IS CAL
     >CULATED')
       call prc ('                       USING JANEV COEFFICIENTS')
       call prc ('                       WITH AN IMPOSED LOWER TEMPERATU
     >RE')
       call prr ('                       LIMIT OF (eV) :',treccut)
c
      elseif (crecopt.eq.3) then
c
       CALL PRC ('  RECOMBINATION OPT 3: HYDROGENIC RECOMBINATION IS CAL
     >CULATED')
       call prc ('                       USING NRL COEFFICIENTS')
       call prc ('                       WITH AN IMPOSED LOWER TEMPERATU
     >RE')
       call prr ('                       LIMIT OF (eV) :',treccut)
c
      elseif (crecopt.eq.4) then
       CALL PRC ('  RECOMBINATION OPT 4: HYDROGENIC RECOMBINATION IS CAL
     >CULATED')
       call prc ('                       USING ADAS COEFFICIENTS')
       call prc ('                       WITH AN IMPOSED LOWER TEMPERATU
     >RE')
       call prr ('                       LIMIT OF (eV) :',treccut)
c
      endif
c
c      PIN CELL AREA OPTION
c
       if (ihcorr.eq.0) then
        call prc('  PIN IHCORR OPT   0 : EDGE2D COMPATIBILITY - USES EDG
     >2D CELL')
        call prc('                       AREAS IN NIMBUS. OTHER EFFECTS
     >MAY ALSO')
        call prc('                       RESULT.')
       elseif (ihcorr.eq.1) then
        call prc('  PIN IHCORR OPT   1 : STANDARD - USES POLYGON AREAS F
     >OR CELLS')
        call prc('                       IN NIMBUS. OTHER EFFECTS MAY AL
     >SO')
        call prc('                       RESULT.')
       endif
c
c      PIN puffing options - option 0 - OFF
c
       if (pinpuff.eq.0) then
c
       call prc ('  PIN PUFF OPT     0 : PIN PUFFING IS OFF')
c
c         Option 1 - ON - recycle pumped gas
c
       elseif (pinpuff.eq.1) then
c
       call prc ('  PIN PUFF OPT     1 : PIN PUFFING IS ON')
c
       call prc ('                       PARTICLES LOST TO ALBEDO (PUMP)
     > ESCAPE')
       call prc ('                       REGIONS ARE REINJECTED WITH THE
     >')
       call prc ('                       FOLLOWING CHARACTERISTICS')
c
       call prr ('                       RE-PUFF FRACTION OF ',hpcpuf)
       call prr ('                       RE-PUFF TEMPERATURE (eV) OF ',
     >                                   tpufh)
c
       if (swpvhpf.eq.0) then
c
         call prc('      LOCATION OPT 0 : FROM MAIN SOL WALLS')
c
       elseif (swpvhpf.eq.1) then
c
         call prc('      LOCATION OPT 1 : FROM PRIVATE VOID WALLS')
c
       endif
       call prc ('                       WITH SEGMENTS SPECIFIED BY:')
       call pri2('                       JHPUF1(1) , JHPUF1(2) = ',
     >    jhpuf1(1),jhpuf1(2))
       call pri2('                       JHPUF2(1) , JHPUF2(2) = ',
     >    jhpuf2(1),jhpuf2(2))

c
c         Option 2 - ON - recycle extra gas
c
       elseif (pinpuff.eq.2) then
c
       call prc ('  PIN PUFF OPT     2 : PIN PUFFING IS ON')
c
       call prc ('                       FRACTION OF RECYLING IONS')
       CALL PRC ('                       ARE NOT RECYLCED BUT ARE ')
       CALL PRC ('                       RATHER RE-INSERTED AS A ')
       CALL PRR ('                       PUFF: PPCPUF = ',ppcpuf)
       call prr ('                       PUFF TEMPERATURE (eV) OF ',
     >                                   tpufh)
c
       if (swpvhpf.eq.0) then
c
        call prc('      LOCATION OPT 0 : FROM MAIN SOL WALLS')
c
       elseif (swpvhpf.eq.1) then
c
        call prc('      LOCATION OPT 1 : FROM PRIVATE VOID WALLS')
c
       endif
       call prc ('                       WITH SEGMENTS SPECIFIED BY:')
       call pri2('                       JHPUF1(1) , JHPUF1(2) = ',
     >    jhpuf1(1),jhpuf1(2))
       call pri2('                       JHPUF2(1) , JHPUF2(2) = ',
     >    jhpuf2(1),jhpuf2(2))

c
c      Endif for pinpuff
c
       endif

c
c     Endif for pincode option above
c
      endif
c
c      Print the Lequiv values 
c
       call prb
c
       call prc('  PIN  : L-Equiv VALUES BASED ON PIN IONIZATION DATA')
       call prc('         Ring        '//Outer//' Plate      '
     >                                 //Inner//' Plate')
       do 390 i = irsep,nrs
          write(coment,'(10X,i4,2(5x,g12.6))') i,cleq(i,2),cleq(i,1)
          call prc(coment)
 390   continue
c
c     END IF for CPINOPT
c
      ENDIF
c
c     Need extra endif?
c
c     endif
c
      RETURN
      END
C
C
C
      CHARACTER*9 FUNCTION FACTOR (A,OPTION)
      implicit none
      REAL A,B
      INTEGER OPTION
      CHARACTER*9 VAL
C     INCLUDE "PARAMS"
      include    'params'
C     INCLUDE "CGEOM"
      include    'cgeom'
C     INCLUDE "COMTOR"
      include    'comtor'
C
      B = A
C
      IF (OPTION.EQ.2) THEN
        IF (B.GT.0.0) THEN
          WRITE (VAL,'(1P,E9.2)') 2.0 * CTEMAV * QTIM / B
        ELSE
          VAL = ' INFIN   '
        ENDIF
C
      ELSEIF (OPTION.EQ.3.OR.OPTION.EQ.4) THEN
        IF (B.GE.1.0) THEN
          VAL = ' INSTANT '
        ELSEIF (B.GT.1.E-3) THEN
          WRITE (VAL,'(1P,E9.2)') -QTIM / LOG (1.0-B)
        ELSEIF (B.GT.0.0) THEN
          WRITE (VAL,'(1P,E9.2)') QTIM / B
        ELSE
          VAL = ' INFIN   '
        ENDIF
C
      ELSE
c        IF (OPTION.EQ.6) B = B / CA
        IF (OPTION.EQ.7) B = B / (QTIM*QTIM*EMI/CRMI)
        IF (OPTION.EQ.8) B = B / QTIM
c        IF (OPTION.EQ.9) B = B / (CA*CA)
        IF (OPTION.EQ.11)B = B / (QTIM*QTIM*EMI/CRMI)
        IF (B.GE.HI) THEN
          VAL = ' INFIN   '
        ELSEIF (B.LE.-HI) THEN
          VAL = '-INFIN   '
        ELSEIF (ABS(B).LE.LO) THEN
          VAL = '   0.0   '
        ELSE
          WRITE (VAL,'(1P,E9.2)') B
        ENDIF
      ENDIF
C
      FACTOR = VAL
      RETURN
      END
C
C
C
      subroutine prnimbin(nimbin,nlines)
      implicit none
      integer nlines
      character*(*) nimbin(nlines)
c
c     PRNIMBIN:
c
c     This routine examines and prints out a meaning
c     when one is avaialble for each of the quantities
c     in the NIMBUS input block that is passed through
c     the DIVIMP input file to NIMBUS.
c
      integer in,rc,len,taglen,valstart
      character*80 tag
      character*16 sp
      integer lenstr,extract_tag
      external lenstr,extract_tag
c
      sp = '                '
c
      call prb
      call prc (sp//'THE NIMBUS INPUT PARAMETERS WERE:')
      call prb
      call prc (sp//'ANALYSIS OF INPUT BLOCK:')
      call prc (sp//'DEFAULT VALUES ARE LISTED IN [ ] ')
      call prb
c
c     Analyse input block - one line at a time
c
      do in = 1,nlines
c
c        Pass the input line and extarct the namelist tag
c
         len = lenstr(nimbin(in))
c
         rc = extract_tag(nimbin(in),tag,taglen,valstart)
c
c         write(6,*) 'NB:',in,':',nimbin(in),':',tag(1:taglen),':',
c     >                   taglen,valstart
c
         if (rc.eq.0) then
c
c           If a valid namelist tag was found.
c
c           Search list to see if it matches a known tag and
c           print the associated meaning and the current value
c           of the tag.
c

C-----------------------------------------------------------------------
            if (tag(1:taglen).eq.'NHIST') then

               call prc(sp//'NHIST: Number of Neutral Histories [2000]')

            elseif (tag(1:taglen).eq.'IFCHAN') then

               call prc(sp//'IFCHAN: 0=No Channels  1=Yes  [1]')

            elseif (tag(1:taglen).eq.'IFWALD') then

               call prc(sp//'IFWALD: 0=off                     [0]')
               call prc(sp//'        1=request distributions of'//
     >                               ' sputtering and power')
               call prc(sp//'          along walls in print file.')

            elseif (tag(1:taglen).eq.'IFPRIM') then

               call prc(sp//'IFPRIM: 0=do not follow impurity'//
     >                               ' neutrals  1=do  [1]')

            elseif (tag(1:taglen).eq.'IZWALL') then

               call prc(sp//'IZWALL: Atomic number of wall')

            elseif (tag(1:taglen).eq.'IAEMIS') then

               call prc(sp//'IAEMIS: 0=Mol. reemission          [0]')
               call prc(sp//'        1=Atomic Reemission')
               call prc(sp//'        -1/-2=INUTPG @ EATMR + AT./MOL')

            elseif (tag(1:taglen).eq.'KINDPR') then

               call prc(sp//'KINDPR: Print Switch 0=minimum  [0]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'TWALL') then

               call prc(sp//'TWALL: Vessel Wall Temperature (C)  [300]')

            elseif (tag(1:taglen).eq.'ZESCUT') then

               call prc(sp//'ZESCUT: Gap polygons above or equal'//
     >                              ' are wall  [INF]')

            elseif (tag(1:taglen).eq.'JXLM') then

               call prc(sp//'JXLM: Knot for projection beyond X-pt.'//
     >                               ' (Use default) [0]')

            elseif (tag(1:taglen).eq.'JXRM') then

               call prc(sp//'JXRM: Knot for projection beyond X-pt.'//
     >                                 ' (Use default) [0]')

            elseif (tag(1:taglen).eq.'XC1') then

               call prc(sp//'XC1: Point for projection beyond X-pt.'//
     >                                  ' (Use default) [RPX]')

            elseif (tag(1:taglen).eq.'YC1') then

               call prc(sp//'YC1: Point for projection beyond X-pt.'//
     >                                 ' (Use default) [ZPX]')

            elseif (tag(1:taglen).eq.'IALB') then

               call prc(sp//'IALB: Albedo condition 0=wall 1=albedo'//
     >                              ' 2=void for P.P. void  [0]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'LWALL') then

               call prc(sp//'LWALL: Use Actual Vessel As Wall'//
     >                                  ' (True/False)  [T]')

            elseif (tag(1:taglen).eq.'LBUFLE') then

               call prc(sp//'LBUFLE: Use Baffle (Set False)'//
     >                                  ' (True/False)  [T]')

            elseif (tag(1:taglen).eq.'LPWALL') then

               call prc(sp//'LPWALL: Use Vessel Wall for private region'
     >)
               call prc(sp//'       (True/False) (forces LBUFLE=F) [F]')

            elseif (tag(1:taglen).eq.'LPSEG') then

               call prc(sp//'LPSEG: Use Explicit Source segments'//
     >                              ' around private void.')
               call prc(sp//'       (Set True!) (True/False)  [F]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'IALBPG') then

               call prc(sp//'IALBPG: Switch for turning specific wall'//
     >                         ' segments into')
               call prc(sp//'albedo regions - rely on pump files'//
     >                        ' - [off]')

            elseif (tag(1:taglen).eq.'ALBPG') then

               call prc(sp//'ALBPG: Switch for turning specific wall'//
     >                         ' segments into')
               call prc(sp//'albedo regions - rely on pump files'//
     >                        ' - [off]')

            elseif (tag(1:taglen).eq.'ALBEPG') then

               call prc(sp//'ALBEPG: Switch for turning specific wall'//
     >                         ' segments into')
               call prc(sp//'albedo regions - rely on pump files'//
     >                        ' - [off]')

            elseif (tag(1:taglen).eq.'ALBATO') then

               call prc(sp//'ALATO: Switch for turning specific wall'//
     >                         ' segments into')
               call prc(sp//'albedo regions - rely on pump files'//
     >                        ' - [off]')

            elseif (tag(1:taglen).eq.'LNWESC') then

               call prc(sp//'LNWESC: Use New Escape Figure Method'//
     >                          ' (Use T) (True/False) [T]')

            elseif (tag(1:taglen).eq.'EATMD') then

               call prc(sp//'EATMD: Energy (eV) of neutrals'//
     >                             ' re-emitted as atoms.')
               call prc(sp//'       0.0=Franck-Condon  [0.025]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'MCX') then

               call prc(sp//'MCX: Model for energy after CX -'//
     >                        ' use default  [0]')

            elseif (tag(1:taglen).eq.'NTSPUT') then

               call prc(sp//'NTSPUT: Turn on neutral sputtering'//
     >                      ' of impurities? 0=off 1=on [1]')

            elseif (tag(1:taglen).eq.'IHOR') then

               call prc(sp//'IHOR: Switch for multi-group vel.'//
     >                        ' distributions. 0=off 1=on [0]')

            elseif (tag(1:taglen).eq.'DECIMA') then

               call prc(sp//'DECIMA: Decimation probability'//
     >                          ' - leave as is -  [0.7]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'MODATM') then

               call prc(sp//'MODATM: Model for atomic CX Losses [1]')

            elseif (tag(1:taglen).eq.'NCOLP') then

               call prc(sp//'NCOLP: Max # of collisons before R.R.'//
     >                             ' analog game (0 to 100,00)')
               call prc(sp//'       (use default) [0]')

            elseif (tag(1:taglen).eq.'ISEHHE') then

               call prc(sp//'ISEHHE: Model for elastic scattering.')
               call prc(sp//'        0=none 1to7=diff. comb.'//
     >                                  ' of H,HZ,HE (use 0) [0]')

            elseif (tag(1:taglen).eq.'RNLITE') then

               call prc(sp//'RNLITE: Reflected fraction'//
     >                             ' of light impurity (?) [1]')

            elseif (tag(1:taglen).eq.'EWLITE') then

               call prc(sp//'EWLITE: Energy (eV) of reflected'//
     >                              ' light impurity (?) [0]')

            elseif (tag(1:taglen).eq.'INUTPG') then

               call prc(sp//'INUTPG: Regions to be set as'//
     >                    ' recyclers (use default) [All Wall]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'EATMR') then

               call prc(sp//'EATMR: Enrgy (eV) of forced'//
     >                               ' reflected neutrals [EATMD]')

            elseif (tag(1:taglen).eq.'TDIV') then

               call prc(sp//'TDIV: Divertor wall'//
     >                                 ' temperature (C) [TWALL]')

            elseif (tag(1:taglen).eq.'LPUMP') then

               call prc(sp//'LPUMP: Switch on pump -'//
     >                          ' in pump structure - (T/F) [T]')

            elseif (tag(1:taglen).eq.'INPUMP') then

               call prc(sp//'INPUMP: Channel for reading pumpfile'//
     >                           ' - SET=18!! - [LDUMIO]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'FPUMP') then

               call prc(sp//'FPUMP: Pump structure file name'//
     >                        ' '' ''=none   ['' '']')

            elseif (tag(1:taglen).eq.'ALBPMP') then

               call prc(sp//'ALBPMP: Albedo for pump'//
     >                            ' <0=use pump value  [-1e30]')

            elseif (tag(1:taglen).eq.'PSEMPO') then

               call prc(sp//'PSEMPO: Transparency of Outer SOL DIV'//
     >                              ' <0= use pump file [-1e30]')

            elseif (tag(1:taglen).eq.'PSEMPT') then

               call prc(sp//'PSEMPT: Transparency of TARGET DIV'//
     >                              ' <0=use pump file [-1e30]')

            elseif (tag(1:taglen).eq.'PSEMPI') then

               call prc(sp//'PSEMPI: Transparency of Inner SOL DIV'//
     >                                 ' <0= use pump file [-1e30]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'PSEMPX') then

               call prc(sp//'PSEMPO: Transparency of CHEVRON'//
     >                                ' <0= use pump file [-1e30]')

            elseif (tag(1:taglen).eq.'IPSEMP') then

               call prc(sp//'IPSEMP: Define wall regions to be'//
     >                               ' semi-transparent')
               call prc(sp//'        Rely on pump files - [none]')

            elseif (tag(1:taglen).eq.'PSEMP') then

               call prc(sp//'PSEMP: Define wall regions to be'//
     >                               ' semi-transparent')
               call prc(sp//'       Rely on pump files - [0]')

            elseif (tag(1:taglen).eq.'PSEMPB') then

               call prc(sp//'PSEMPB: Transparency of baffle'//
     >                              '  1e30 = do not use  [1e30]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'IPVOID') then

               call prc(sp//'IPVOID: Flag for treatment of'//
     >                          ' pump void walls - use -1! - [1]')

            elseif (tag(1:taglen).eq.'IVIEW') then

               call prc(sp//'IVIEW: 0=std. Nimbus geom. map'//
     >                             '  -1=user defined window  [0]')

            elseif (tag(1:taglen).eq.'VIEW') then

               call prc(sp//'VIEW: (Rmin,Zmin,Rlen, Zlen) of'//
     >                         ' user defined window')
               call prc(sp//'      GEOM from GRID2D  -  [GEOM]')

            elseif (tag(1:taglen).eq.'ITRIM') then

               call prc(sp//'ITRIM: 0=no TRIM files '//
     >                              ' 1=use TRIM files -use 1!- [0]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'FTRIM') then

               call prc(sp//'FTRIM: TRIM file prefix -'//
     >                                 ' set in /NIMBIN/ - [CTRIMF]')

            elseif (tag(1:taglen).eq.'LFULL') then

               call prc(sp//'LFULL: Full setup for NIMBUS'//
     >                          ' at every call (T/F) - use T -  [T]')

            elseif (tag(1:taglen).eq.'ITARHZ') then

               call prc(sp//'ITARHZ: Switches (2) to'//
     >                                 ' determine when to use')
               call prc(sp//'        horizontal escape figure'//
     >                                     ' - use default - [2*MX]')

            elseif (tag(1:taglen).eq.'ICHKP') then

               call prc(sp//'ICHKP: 0=Stop on polygon'//
     >                                 ' problems in NIMBUS')
               call prc(sp//'       1=warn - use default - [1]')

            elseif (tag(1:taglen).eq.'LTIME') then

               call prc(sp//'LTIME: Time dependent Monte Carlo'//
     >                            ' (T/F) - use F - [F]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'TWIDTH') then

               call prc(sp//'TWIDTH: Time slice width -'//
     >                          ' only for time-dependent mode - [0]')

            elseif (tag(1:taglen).eq.'TWDMIN') then

               call prc(sp//'TWDMIN: Minimum time slice width'//
     >                      ' - only for time dependent mode - [0]')

            elseif (tag(1:taglen).eq.'AYIZ') then

               call prc(sp//'AYIZ: Enhanced Yield from H: Y''=AY+B [1]')

            elseif (tag(1:taglen).eq.'BYIZ') then

               call prc(sp//'BYIZ: Enhanced Yield from H: Y''=AY+B [0]')

            elseif (tag(1:taglen).eq.'ICORRN') then

               call prc(sp//'ICORRN: Random number correlation flag'//
     >                                   ' - use default - [UNDEF]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'ALBLK') then

               call prc(sp//'ALBLK: Albedo of pump leaks'//
     >                               '  <0=use pump file  - [-1e30]')

            elseif (tag(1:taglen).eq.'DCUTCX') then

               call prc(sp//'DCUTCX: Maximum density for CX  [1e30]')

            elseif (tag(1:taglen).eq.'TCUTCX') then

               call prc(sp//'TCUTCX: Minimum temperature for CX'//
     >                                                  '  [-1e30]')

            elseif (tag(1:taglen).eq.'ITAU') then

               call prc(sp//'ITAU: Model for flux estimator       [1]')
               call prc(sp//'      <=1=estimated dist.  >1=dist.-'//
     >                                          ' use default - [1]')

            elseif (tag(1:taglen).eq.'IYCHEM') then

               call prc(sp//'IYCHEM: Model for chemical sputtering')
               call prc(sp//'        0=off 1+=model selected - [0]')

C-----------------------------------------------------------------------
            elseif (tag(1:taglen).eq.'EYCHEM') then

               call prc(sp//'EYCHEM: Energy (eV) of chemically'//
     >                                ' sputtered C ''atom'' - [0.0]')

            elseif (tag(1:taglen).eq.'IDBHST') then

               call prc(sp//'IDBHST: Number of histories to store'//
     >                            ' trajectories - use default - [0]')

            elseif (tag(1:taglen).eq.'XGAUGE') then

               call prc(sp//'XGAUGE: R position -'//
     >              ' to override gauge location in pump file - [1e30]')

            elseif (tag(1:taglen).eq.'YGAUGE') then

               call prc(sp//'YGAUGE: Z position -'//
     >              ' to override gauge location in pump file - [1e30]')

            elseif (tag(1:taglen).eq.'RGAUGE') then

               call prc(sp//'RGAUGE: Radius -'//
     >              ' to override gauge location in pump file - [1e30]')

            elseif (tag(1:taglen).eq.'LGAUGE') then

               call prc(sp//'LGAUGE: Label - ''G'' for pump gauge'//
     >                                        ' (leave as is)')
               call prc(sp//'        ''K'' for vessel - [''G'',''K'']')

            elseif (tag(1:taglen).eq.'MODEZR') then

               call prc(sp//'MODEZR: Model for impurity ion'//
     >                               ' recycling - use default - [1]')

            elseif (tag(1:taglen).eq.'ISPOFF') then

               call prc(sp//'ISPOFF: Switch off recycling in'//
     >                                    ' specifed macrozones - [0]')

            elseif (tag(1:taglen).eq.'MIMP') then

               call prc(sp//'MIMP: # imp switch -'//
     >              ' still being developed - use 0 for now - [0]')

            elseif (tag(1:taglen).eq.'LRS') then

               call prc(sp//'LRS: Override leak recycling segments'//
     >                          ' in pump file :')
               call prc(sp//'     ''<ID> X1 Y1 X2 Y2'' - [60*'' '']')

            elseif (tag(1:taglen).eq.'GAP') then

               call prc(sp//'GAP: Override Gap segments in pump file :')
               call prc(sp//'     ''<ID> X1 Y1 X2 Y2 T'' '//
     >                                 '(transparency) - [60*'' '']')

c
c           NEW Input values added for NIMBUS/PIN release 5.0
c

            elseif (tag(1:taglen).eq.'YNOSP') then

               call prc(sp//'YNOSP: Switch off sputtering for'//
     >                      ' Y < YNOSP (cm) - [-1.0e30 = OFF]')

            elseif (tag(1:taglen).eq.'LMC') then

               call prc(sp//'LMC: Test switch - if TRUE run link'//
     >                      ' but do not run NIMBUS - [FALSE]')

            elseif (tag(1:taglen).eq.'INIZ') then

               call prc(sp//'INIZ: Random number seed for NIMBUS'//
     >                      ' - [INIZ=1]')

            elseif (tag(1:taglen).eq.'NOZPMP') then

               call prc(sp//'NOZPMP: Allow pumping of impurities'//
     >                      ' - 0=yes 1=no  [0]')

            elseif (tag(1:taglen).eq.'AYCHEM') then

               call prc(sp//'AYCHEM: Multiplicative factor for'//
     >                      ' Chemical Sputtering'//
     >                      ' - [1.0]')

            elseif (tag(1:taglen).eq.'PWMAT') then

               call prc(sp//'PWMAT: Pump wall material '//
     >                      ' - set to "FE" as default  ["  "]')
               call prc(sp//'       - NIMBUS default is the same'//
     >                      'effect as IPVOID=-1')
               call prc(sp//'       - use this instead of IPVOID=-1')

            elseif (tag(1:taglen).eq.'TDIVW') then

               call prc(sp//'TDIVW: Divertor Wall Temperature (C)'//
     >                      ' - [TWALL]')

            elseif (tag(1:taglen).eq.'TTARG') then

               call prc(sp//'TTARG: Target Temperature (C)'//
     >                      ' - [TWALL]')

            elseif (tag(1:taglen).eq.'TPRIV') then

               call prc(sp//'TPRIV: Private Region Wall'//
     >                      ' Temperature (C)'//
     >                      ' - [TWALL]')

            elseif (tag(1:taglen).eq.'TPRIV') then

               call prc(sp//'TPRIV: Sub-divertor Region Wall'//
     >                      ' Temperature (C)'//
     >                      ' - [TWALL]')

            elseif (tag(1:taglen).eq.'TEVGAP') then

               call prc(sp//'TEVGAP: Temperature of Void Gaps (eV)'//
     >                      ' - [0.05]')
               call prc(sp//'        Not used at present - could be '
     >                    //'used for sputtering from gaps')

            elseif (tag(1:taglen).eq.'RECMAT') then

               call prc(sp//'RECMAT: Reflection matrix for fuel')
               call prc(sp//'        - Reemission coefficients'//
     >                      ' for implanted atoms.')
               call prc(sp//'        - [Identity Matrix is default]')
               call prc(sp//'        - used for multi-species hydrogen'
     >                    //' runs.')

            elseif (tag(1:taglen).eq.'ANGNOR') then

               call prc(sp//'ANGNOR: Angle with wall normal for '//
     >                      ' H-flux scores - [90.0]')
               call prc(sp//'        Angle not equal to default may '
     >                    //'not be compatible with')
               call prc(sp//'        flux dependent chemical'//
     >                      ' sputtering.')
               call prc(sp//'        Default value is recommended.')

            elseif (tag(1:taglen).eq.'NPWALL') then

               call prc(sp//'NPWALL: Number of extra points inserted '//
     >                      ' between private region')
               call prc(sp//'        and end-point extrapolations. [0]')

            elseif (tag(1:taglen).eq.'FPWALL') then

               call prc(sp//'FPWALL: Array of positions of  extra'//
     >                      ' points inserted between private region')
               call prc(sp//'        end-point extrapolations,'//
     >                       ' expressed as a fraction of the')
               call prc(sp//'        distance along the wall.'
     >                    //' (monotonic and in the'//
     >                      ' range (0,1) exclusive). ')

            elseif (tag(1:taglen).eq.'ALFCH7') then

               call prc(sp//'ALFCH7: Exponent for flux dependency of '//
     >                      ' Haasz97 chemical sputtering yields [0.1]')

c
C-----------------------------------------------------------------------
c
            else
c
               call prc(sp//tag(1:taglen)//': Unidenitified Tag')

            endif

            call prc(sp//nimbin(in)(1:len))
            call prb

C-----------------------------------------------------------------------
         else
c
c           If there was no tag identified - just echo the line to
c           the output
c
            call prc(sp//nimbin(in)(1:len))
c
         endif
c
      end do
c
      call prc (sp//'NOTE: Any parameters not listed were')
      call prc (sp//'      used with their built-in NIMBUS')
      call prc (sp//'      default values.')
      call prb
c
      return
      end
c
c
c
      integer function extract_tag(string,tag,taglen,valstart)
      implicit none
      character*(*) string
      character*(*) tag
      integer taglen,valstart
c
c     EXTRACT_TAG: This routine processes a line of input - it assumes
c                  that only one Identifier tag will be entered on a
c                  line - it assumes that this tag will be immediately
c                  before the first "=" on the line - after removing
c                  any whitespace - it then strips any whitespace
c                  around this identifier and copies this substring
c                  into the tag character string. It returns this
c                  tag - the number of characters in it in taglen and
c                  the location of the "=" +1 in the original
c                  string in the valstart - so that other code
c                  might be able to extract the actual value entered
c                  on the line at a later date. It returns a value of
c                  0 if a tag is successfully extracted - a value of
c                  1 if the line contains no "=" and a value of 2 if
c                  the tag length is found to be zero.
c
      integer in,rc,pos,endn
c
      extract_tag = 0
c
c     Find '='
c
      pos = index(string,'=')
c
      if (pos.eq.0) then
         extract_tag = 1
         return
      endif
c
      valstart = pos+1
c
c     Find first non-space before '='
c
      endn = pos
c
 10   endn = endn -1
c
      if (endn.le.0) then
         extract_tag=2
         return
      endif
c
      if (string(endn:endn).eq.' ') goto 10
c
      extract_tag = 0
c
      do in = endn,1,-1

         if (in.eq.1) then
            tag=string(in:endn)
            taglen = endn-in+1
            return
         elseif (string(in-1:in-1).eq.' ') then
            tag=string(in:endn)
            taglen = endn-in+1
            return
         endif
c
      end do
c
c     Code should never get here.
c
      return
      end
c
c
c
      subroutine prtemp
      implicit none
      include 'params'
      include 'comtor'
      include 'cgeom'
c
c     PRTEMP: This routine prints out the wall temperatures that
c             were specifed for this case.
c
      character*77 comment
      integer in
c
      call prb
      call prc('  WALL AND TARGET TEMPERATURES IN USE:')
      call prr('    Default Target Temperature  (K) =',ctargt)
      call prr('    Default Vessel Wall Temp.   (K) =',cwallt)
      call prr('    Default PP Wall Temperature (K) =',cwallp)
      if (nwltemp.gt.0) then
         call prc('  THESE VALUES ARE OVERRIDDEN FOR THE'//
     >                            ' FOLLOWING SEGMENTS:')
         call prc('    Index #1   to   Index #2     '//
     >                            'Temperature (K)')
         do in = 1,nwltemp
            write(comment,'(4x,i6,10x,i6,8x,f11.2)')
     >                        int(walltemp(in,1)),
     >                        int(walltemp(in,2)),walltemp(in,3)
            call prc(comment)
         end do
c
      endif
c
      return
      end
c
c
c
      subroutine appendfile(outunit, inunit)
      implicit none
c
      integer outunit, inunit,len,lenstr
      character line*512 
      external lenstr
c
      write(6,*) 'APPENDFILE:',outunit,inunit
c
c     This routine appends the contents of one file to the other. 
c
      rewind(inunit)
c
 10   read(inunit,'(A512)',end=20,err=20) line
c
      len = lenstr(line)
c
      write(outunit,'(a)') line(1:len) 
c
      goto 10
c
 20   return
c
      end 
c
c
c
      subroutine pr_jhfact
      implicit none
      include 'printopt'
c
c     Print out a summary of the Johnson-Hinov ratios over a variety
c     of regions.      
c   
      real jhtots(10,3)
      character comment*72
c
      call calc_jhfactors(jhtots)
c
      call prb
      call prc('  SUMMARY OF JOHNSON-HINOV FACTORS BY REGION:')
      call prc('         (EXTRACTED FROM PIN RESULTS)')
      call prb
c
c     Title
c
      write(comment,'(5x,a6,5x,4x,a7,4x,2x,a11,2x,3x,a9,3x)')
     >          'REGION','PHOTONS','IONIZATIONS','JH-FACTOR'
      call prc(comment)
c
c     Regions
c

      write(comment,10) 'CORE',
     >                  jhtots(1,1),jhtots(1,2),jhtots(1,3)
      call prc(comment)
c
      write(comment,10) OUTER//' DIVERTOR',
     >                  jhtots(2,1),jhtots(2,2),jhtots(2,3)
      call prc(comment)
c
      write(comment,10) OUTER//' MAIN VES',
     >                  jhtots(3,1),jhtots(3,2),jhtots(3,3)
      call prc(comment)
c
      write(comment,10) OUTER//' SOL',
     >                  jhtots(4,1),jhtots(4,2),jhtots(4,3)
      call prc(comment)
c
      write(comment,10) INNER//' DIVERTOR',
     >                  jhtots(5,1),jhtots(5,2),jhtots(5,3)
      call prc(comment)
c
      write(comment,10) INNER//' MAIN VES',
     >                  jhtots(6,1),jhtots(6,2),jhtots(6,3)
      call prc(comment)
c
      write(comment,10) INNER//' SOL',
     >                  jhtots(7,1),jhtots(7,2),jhtots(7,3)
      call prc(comment)
c
      write(comment,10) 'ENTIRE SOL',
     >                  jhtots(8,1),jhtots(8,2),jhtots(8,3)
      call prc(comment)
c
      write(comment,10) 'PFZ',
     >                  jhtots(9,1),jhtots(9,2),jhtots(9,3)
      call prc(comment)
c
      write(comment,10) 'ENTIRE PLASMA',
     >                  jhtots(10,1),jhtots(10,2),jhtots(10,3)
      call prc(comment)
c
      call prb
c
c     Format statements 
c
 10   format(1x,a14,1x,1x,1p,g13.5,2x,g13.5,2x,0p,f13.5)
c
      return
      end
c
c
c
      subroutine pr_power_summary
      implicit none
c
      include 'params'
      include 'slcom'           
c
      call prb  
      call prchtml('GLOBAL POWER SUMMARY','pr_power','0','B')
      call prb 



      call pr_heatfluxdata




      return
      end
c
c
c
      subroutine pr_heatfluxdata
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'printopt'
c
c     PR_HEATFLUXDATA: 
c
c     This routine prints a table that contains a summary of the 
c     heat flux information for each element of the wall.
c
      integer id,in,is,it
c
      call prb
      call prc('SUMMARY OF WALL HEAT FLUXES (W/M^2):')       
c
c     Inner target
c
      call prb
      call prc('ELEM TYPE   R      Z     ION     ATOM'//
     >      '     MOL   PTOTAL IMP RAD  H RAD  TOT RAD  GTOTAL')
c
      do in= 1,wallpts
c         
         id = wallpt(in,18)
c
         write(coment,'(i4,i4,2f7.3,1p,8(1x,e7.1))') 
     .          in,id,wallpt(in,1),wallpt(in,2),
     >          (wallfluxdata(in,is,4),is=1,4),
     >          (wallprad(in,is),is=1,3),
     >          wallfluxdata(in,4,4)+wallprad(in,3)
c
         call prc(coment)
c
      end do
c
c     Totals:
c
      call prb
      call prc('WALL REGION TOTALS [W]:')
      call prb 
      call prc('REGION        ION     ATOM   '//
     >      '   MOL    PTOTAL   IMP RAD'//
     >      '   H RAD   TOT RAD  GTOTAL')
c
c     First Target 
c
      write(coment,'(a9,1x,1p,8(e9.2))') 'TOT '//inner,
     >          (wallfluxdata(maxpts+1,is,4),is=1,4),
     >          (wallprad(maxpts+1,is),is=1,3),
     >          wallfluxdata(maxpts+1,4,4)+wallprad(maxpts+1,3)
      call prc(coment)
c
c     Second target 
c
      write(coment,'(a9,1x,1p,8(e9.2))') 'TOT '//outer,
     >          (wallfluxdata(maxpts+2,is,4),is=1,4),
     >          (wallprad(maxpts+2,is),is=1,3),
     >          wallfluxdata(maxpts+2,4,4)+wallprad(maxpts+2,3)
      call prc(coment)
c
c     MAIN Wall
c
      write(coment,'(a9,1x,1p,8(e9.2))') 'TOT  MAIN',
     >          (wallfluxdata(maxpts+3,is,4),is=1,4),
     >          (wallprad(maxpts+3,is),is=1,3),
     >          wallfluxdata(maxpts+3,4,4)+wallprad(maxpts+3,3)
      call prc(coment)
c
c     PFZ Wall
c
      write(coment,'(a9,1x,1p,8(e9.2))') 'TOT   PFZ',
     >          (wallfluxdata(maxpts+4,is,4),is=1,4),
     >          (wallprad(maxpts+4,is),is=1,3),
     >          wallfluxdata(maxpts+4,4,4)+wallprad(maxpts+4,3)
      call prc(coment)
c
c     Grand Totals for wall
c
      write(coment,'(a9,1x,1p,8(e9.2))') 'TOT   ALL',
     >          (wallfluxdata(maxpts+5,is,4),is=1,4),
     >          (wallprad(maxpts+5,is),is=1,3),
     >          wallfluxdata(maxpts+5,4,4)+wallprad(maxpts+5,3)
      call prc(coment)
c
      call prb
c
c     Radiation source strength
c
      write(coment,'(a9,1x,36x,1p,8(e9.2))') 'RAD   SRC',
     >          (wallprad(maxpts+6,is),is=1,3)
      call prc(coment)
c
c
      return 
      end




      subroutine pr_targfluxdata
      implicit none
      include 'params'
      include 'cgeom'
      include 'comtor'
      include 'printopt'
c
c     PR_TARGFLUXDATA: 
c
c     This routine prints a table that contains a summary of the 
c     heat flux information for each element of the target.
c
c
      integer id,in
c
      call prb
      call prc('SUMMARY OF TARGET PARTCLE'//
     >         ' HEAT FLUXES (W/M^2):')       
      call prc('  - FOR THE RADIATIVE SUMMARY AND FULL WALL'//
     >         ' HEAT FLUXES SEE THE')
      CALL PRC('    GLOBAL HEAT FLUX SUMMARY TABLE'//
     >         ' TOWARD THE END OF THIS FILE.')
c
c     Inner target
c
      call prb
      call prc(inner//' Target:')
      call prc('ELEMENT     R           ION        ATOM'//
     >         '         MOLECULE        TOTAL')
c
      do id= 1,ndsin
c         
         write(coment,'(2x,i4,1x,f8.4,4x,1p,4(1x,e12.4))') 
     .     id,rp(id),
     >          (targfluxdata(id,in,4),in=1,4)
c
         call prc(coment)
c
      end do
c
      write(coment,'(a9,10x,1p,4(1x,e12.4),2x,a)') 'TOT '//inner,
     >          (targfluxdata(maxnds+1,in,4),in=1,4),'(W)'
      call prc(coment)
c
c     Outer target 
c
      call prb
      call prc(outer//' Target:')
      call prc('ELEMENT     R           ION        ATOM'//
     >         '         MOLECULE        TOTAL')
c
      do id= ndsin+1,nds
c         
         write(coment,'(2x,i4,1x,f8.4,4x,1p,4(1x,e12.4))') id,rp(id),
     >          (targfluxdata(id,in,4),in=1,4)
c
         call prc(coment)
c
      end do
c
      write(coment,'(a9,10x,1p,4(1x,e12.4),2x,a)') 'TOT '//outer,
     >          (targfluxdata(maxnds+2,in,4),in=1,4),'(W)'
      call prc(coment)
c
      call prb      
c
      write(coment,'(a9,10x,1p,4(1x,e12.4),2x,a)') 'GRAND TOT',
     >          (targfluxdata(maxnds+3,in,4),in=1,4),'(W)'
      call prc(coment)
c
      call prb
c
      return 
      end
c
c
c
! ammod begin.
      Subroutine pr_hydrocarbon_options
c
      Use ComHC ! Contains optional hydrocarbon options
                ! selected in input file (H15-H51).

      ! Every good Fortran routine has...
      Implicit none
      
      Include 'comtor'
      
      ! Declare local variables.
      Integer I
      Logical higher_hcs_available
      
      ! Initialize local variables.
      higher_hcs_available = .false.
      
c
c
c***********************************************************************
c
c     THIS ROUTINE PRINTS ALL THE HYDROCARBON OPTION FLAGS
c
c***********************************************************************
c
c
c-----------------------------------------------------------------------
      Call PRB
      Call PRChtml ('--- HYDROCARBON OPTIONS ---',
     >              'pr_hydrocarbon_options','0','B')
      Call PRB
c-----------------------------------------------------------------------
      ! (H15)
      If (hc_follow_option.eq.0) Then
         ! Hydrocarbon following not activated.
         Call PRC ('  Hydrocarbon following off.')
      ElseIf (hc_follow_option .eq. 1) Then
         ! Hydrocarbon following activated.
         Call PRC ('  Hydrocarbon following on.')

         call prb
         ! (H64) jdemod
         call prr('H64: Default H Isotope Mass in HC molecules has'//
     >           ' been set to:',input_HC_H_mass) 
         call prb

         ! Print all other options.
         ! (H16)
	 If (hc_higher_hcs_option .eq. 0) Then
            Call PRC ('  Hydrocarbons beyond CH4 ignored.')
         ElseIf (hc_higher_hcs_option .eq. 1) Then
            Call PRC ('  Hydrocarbons beyond CH4 (up to C3H8)'//
     >                ' followed.')
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        'Error:  Option for hc_higher_hcs_option'//
     >        ' not supported:', hc_higher_hcs_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H17)
	 If (hc_wbc_comp_option .eq. 0) Then
            Call PRC ('  WBC comparison output not printed.')
         ElseIf (hc_wbc_comp_option .eq. 1) Then
            Call PRC ('  WBC following and tracking enabled'//
     > ' followed.')
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_wbc_comp_option'//
     >        ' not supported:", hc_wbc_comp_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H20)
         If (hc_sputtering_model .eq. 0) Then
            Call PRC ('  Preset single HC species will be emitted.')
            ! (H21) dependent on H20.
            Call PRI ('  HC species launched: ',
     > hc_sputtered_hc_species)
         ElseIf (hc_sputtering_model .eq. 1) Then
            Call PRC ('  Mech, Davis, Haasz (Nucl.Mat 1997) HC launch'//
     > ' distribution used (with possibilty of launch of CH4 - C3H8)')
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        'Error:  Option for hc_sputtering_model'//
     >        ' not supported:', hc_sputtering_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H22)
         If (hc_evolution_model_primary .eq. 0) Then
	    ! No data.
	    Write (Output_Unit_HC_Alert,*)
     >        "Error: No hydrocarbon reaction data specified"//
     >        " in input file (H22)."
	    Write (Output_Unit_HC_Alert,*) "Program stop."
            Stop
	 ElseIf (hc_evolution_model_primary .eq. 1) Then
	    ! E&L data selected as primary.
	    Call PRC ('  Ehrnhardt and Langer (PPPL, 1987) data is'//
     > ' used as primary.')
            higher_hcs_available = .false.
	 ElseIf (hc_evolution_model_primary .eq. 2) Then
	    ! Alman/Ruzic/Brooks data selected as primary.
	    Call PRC ('  Alman, Ruzic, Brooks (Phys. Plasma, 2000)'//
     > ' data is used as primary.')
            higher_hcs_available = .true.
	 ElseIf (hc_evolution_model_primary .eq. 3) Then
	    ! Janev data selected as primary.
            Call PRC ('  Janev, et al. (NIFS, 2001) data is used as'//
     > ' primary.')
            higher_hcs_available = .true.
	 End If

         ! (H23)
         ! jdemod - this option doesn't work at present and is being 
         !          commented out for now
c         If (hc_evolution_model_secondary .eq. 0) Then
c	    ! E&L data will be used t.  Check that only lower hydrocarbons are being followed.
c	    Call PRC ('  No data will be used as a secondary set.')
c         ElseIf (hc_evolution_model_secondary .eq. 1) Then
c	    ! E&L data will be used as secondary.
c	    Call PRC ('  Ehrnhardt and Langer (PPPL, 1987) data is'//
c     > ' used as secondary (to compliment primary reaction data).')
c	 ElseIf (hc_evolution_model_secondary .eq. 2) Then
c            ! A,R,B data will be used second.
c            Call PRC ('  Alman, Ruzic, Brooks (Phys. Plasma, 2000)'//
c     > ' data is used as secondary (to compliment primary reaction'//
c     > ' data).')
c            higher_hcs_available = .true.
c         ElseIf (hc_evolution_model_secondary .eq. 3) Then
c            ! Janev data will be used second.
c            Call PRC ('  Janev, et al. (NIFS, 2001) data is used as'//
c     > ' secondary (if primary reaction data missing).')
c            higher_hcs_available = .true.
c         Else
c            ! Unsupported reaction or order specified.
c            Write (Output_Unit_HC_Alert,*) 
c     >        "Error: Unsupported reaction or order"//
c     >        " specified (H23)."
c            Write (Output_Unit_HC_Alert,*) "Program stop."
c            Stop
c     	 End If
       
         ! Check to be sure if higher hcs are being followed that the data (Janev or Alman/Ruzic/Brooks) is available.
c slmod begin
         If ((hc_higher_hcs_option .eq. 1) .And. 
     .       (.not.higher_hcs_available)) Then
c
c         If ((hc_higher_hcs_option .eq. 1) .And. (higher_hcs_available
c     >     .eq. .false.)) Then
c slmod end
            ! Error - must use additional reaction rate date.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  More data required for higher"//
     >        " hydrocarbons."
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop         
         End If

         ! (H24)
         If (hc_launch_location .eq. -1) Then
	    ! Defaults to DIVIMP Launch option CNEUTB.
	    Call PRI ('  HC launch option set to DIVIMP launch option: '
     > ,hc_launch_location)
         ElseIf (hc_launch_location .ge. 0 .and. 
     > hc_launch_location .le. 6) Then
            ! DIVIMP launch location option used.
            Call PRI ('  HC launch option set independent of DIVIMP'//
     > ' launch option CNEUTB as: ',hc_launch_location)
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_launch_location"//
     >        " not supported:", hc_launch_location
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H25)
         If (hc_launch_angle_velocity .eq. -1) Then
	    ! Defaults to DIVIMP velocity/angle option CNEUTC.
	    Call PRI ('  HC vel/ang option set to DIVIMP V/A option: '
     > ,hc_launch_angle_velocity)
         ElseIf (hc_launch_angle_velocity .ge. 0 .and.
     > hc_launch_angle_velocity .le. 20) Then
            ! DIVIMP velocity/angle option used.
            Call PRI ('  HC vel/ang option set independent of DIVIMP'//
     > ' launch option CNEUTC as: ',hc_launch_angle_velocity)
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_launch_angle_velocity"//
     >        " not supported:", hc_launch_angle_velocity
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H26)
	 
	 ! (H27)
	 
	 ! (H28)
	 


         ! (H30)
         If (hc_neut_ion_velocity .eq. -1) Then
	    ! Defaults to DIVIMP initial ion velocity option CNEUTG.
	    Call PRI ('  HC initial ion velocity option set to DIVIMP'//
     > ' initial ion velocity option CNEUTG: ',hc_neut_ion_velocity)
         ElseIf (hc_neut_ion_velocity .ge. 0 .and.
     > hc_neut_ion_velocity .le. 3) Then
            ! DIVIMP initial ion velocity option used.
            Call PRI ('  HC initial ion velocity option set'//
     > ' independent of DIVIMP ion velocity option CNEUTG as: ',
     > hc_neut_ion_velocity)
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_neut_ion_velocity"//
     >        " not supported:", hc_neut_ion_velocity
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H31)
         If (hc_ion_neut_angle .eq. 0) Then
	    ! Isotropic ejection distribution after neutralization.
	    Call PRC ('  Isotropic ejection distribution after'//
     > ' neutralization is used.')
         ElseIf (hc_ion_neut_angle .eq. 1) Then
            ! Sine bias forward distribution after neutralization.
            Call PRC ('  Sine-biased forward distribution after'//
     > ' neutralization is used.')
         ElseIf (hc_ion_neut_angle .eq. 2) Then
            ! S-direction at point of neutralization.
            Call PRC ('  S-direction at point of neutralization'//
     > ' is used.')
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >       "Error:  Option for hc_ion_neut_angle"//
     >       " not supported:", hc_ion_neut_angle
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H32)
         If (hc_ion_neut_velocity .eq. 0) Then
	    ! Original ion energy is used for neutral energy upon neutralization.
	    Call PRC ('  Ion energy is used for energy of particle'//
     > ' after neutralization.')
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_ion_neut_velocity"//
     >        " not supported:", hc_ion_neut_velocity
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H33)
	 If (hc_lambda_calc .eq. 0) Then
	    ! Using standard L=15 throughout grid for ion transport.
	    Call PRC ('  Using L=15 for all ion transport.')
	 ElseIf (hc_lambda_calc .eq. 1) Then
	    ! Using enhanced calculation of Sivukhin for L.
	    Call PRC ('  Using temperature dependant calculation'//
     >        ' of Sivukhin for Lambda.')
	 Else
	    ! Error: Option not supported.
	    Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_lambda_calc"//
     >        " not supported:", hc_lambda_calc
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If
	 
	 ! (H34)
	 If (hc_disable_transitions .eq. 0) Then
	    ! Hydrocarbon transitions are not disabled.
	    Call PRC ('  Hydrocarbons are free to evolve between'//
     >        ' species')
	 ElseIf (hc_disable_transitions .eq. 1) Then
	    ! Hydrocarbon transitions are disabled.
	    Call PRC ('  Hydrocarbon transitions are disabled. '//
     >       'Note: This should only be used for debugging '//
     >       'purposes.')
	 Else
	    ! Error: Option not supported.
	    Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_disable_transitions"//
     >        " not supported:", hc_disable_transitions
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If	 
	
	 ! (H35) and (H36) and (H37)
         If (hc_presheath_efield .eq. 0) Then
            ! Standard sheath treatment is kept.
	    Call PRC ('  Using standard sheath e-field treatment')    
         ElseIf (hc_presheath_efield .eq. 1) Then
            ! Using enhanced calculation for sheath e-field.
            Call PRC ('  Using enhanced sheath e-field of Brooks.')
            ! Print Debye sheath drop fraction.
            Call PRR ('  Voltage drop fraction in the Debye region:',
     >        hc_efield_drop_fraction)
            ! Print range of applicable cells from target.
	    Call PRI ('  Applied to number of cells from target:',
     >        hc_efield_cells)
	 Else
	    ! Error: Option not supported.
	    Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_presheath_efield"//
     >        " not supported:", hc_presheath_efield
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If	 
	 
         ! (H40) and (H43)
	 If (hc_neutral_reflection_option .eq. 0) Then
	    ! Reflection of neutrals is off.
	    Call PRC ('  Neutral hydrocarbons are not reflected from'//
     >        ' the vessel.')
	 ElseIf (hc_neutral_reflection_option .eq. 1) Then
	    ! Reflection of neutrals is on.
	    Call PRC ('  Neutral hydrocarbons are reflected from the'//
     >        'vessel.')
            ! (H40) dependent on H43.
            Call PRR ('  Preset neutral reflection coefficient',
     > hc_reflection_coef_preset)
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_neutral_reflection"//
     >        "_option not supported:", hc_neutral_reflection_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (41) and (H43)
	 If (hc_ion_reflection_option .eq. 0) Then
	    ! Reflection of ions is off.
	    Call PRC ('  Ionized hydrocarbons are not reflected from'//
     >        ' the vessel.')
	 ElseIf (hc_ion_reflection_option .eq. 1) Then
	    ! Reflection of ions is on.
	    Call PRC ('  Ionized hydrocarbons are reflected from the'//
     >        'vessel.')
            ! (H41) dependent on H43.
            Call PRR ('  Preset ion reflection coefficient',
     > hc_reflection_coef_preset)
         Else
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_ion_reflection"//
     >        "_option not supported:", hc_ion_reflection_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
         End If

         ! (H42)
         If (hc_reflection_coef_model .eq. 0) Then
            ! Use preset reflection coefficient in all cases (typically 0.5).
            Call PRC ('  Preset global sticking coefficient is used.')
         ElseIf (hc_reflection_coef_model .eq. 1) Then
            Call PRC ('  Janev data used.'//
     > ' Not currently supported.')
         ElseIf (hc_reflection_coef_model .eq. 2) Then
            Call PRC ('  Alman and Ruzic data used.'//
     > ' Results based on average wall impact energy and'//
     > ' should be considered preliminary.')
         ElseIf (hc_reflection_coef_model .eq. 3) Then
            Call PRC ('  Refl coef for CH4 = 1.0, and C0=0.0.'//
     > ' Other HC fragment coeffs will use that in option H43.')
         ElseIf (hc_reflection_coef_model .eq. 4) Then
            Call PRC ('  Refl coef from Jacob JNM 337-339.'//
     > ' (2005) 839-846 (table 1 p844)')
         ElseIf (hc_reflection_coef_model .eq. 5) Then
            Call PRC ('  All molecular HC species are reflected'//
     > ' at the preset value. C and C+ are NOT reflected.')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_reflection_coef_model"//
     >        " not supported:", hc_reflection_coef_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
            Write (0,*) 
     >        "Error:  Option for hc_reflection_coef_model"//
     >        " not supported:", hc_reflection_coef_model
            Write (0,*) "Program stop."
	    Stop
	 End If

         ! (H44)
         If (hc_reflection_species_model .eq. 0) Then
            ! Preset reflection table used for reflected hydrocarbon species.
            Call PRC ('  Preset sputtering table used for reflected'//
     >        ' hydocarbon species.')
         ElseIf (hc_reflection_species_model .eq. 1) Then
	    ! Alman and Ruzic data used for reflected hydrocarbon species.
            Call PRC ('  Alman and Ruzic data used for reflected'//
     >        ' hydrocarbon species.')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) "Error:  Option for "//
     >        " hc_reflection_species_model not supported:", 
     >        hc_reflection_species_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H45) and (H46) and (H47)
         If (hc_reflection_energy_model .eq. 0) Then
            ! Reflected hydrocarbon will always have preset energy set in option H31 (typically a few ev total).
	    Call PRR ('  Reflected hydrocarbon from neutral impact '//
     > 'will have energy of (eV): ', hc_refl_energy_neutral_preset) ! H46.
	    Call PRR ('  Reflected hydrocarbon from ion impact '//
     > 'will have energy of (eV): ', hc_refl_energy_ion_preset) ! H47.
         ElseIf (hc_reflection_energy_model .eq. 1) Then
	    ! Reflected hydrocarbon energy is the impact energy.
            Call PRC ('  Reflected hydrocarbon energy equals'//
     > ' impact energy.')
         ElseIf (hc_reflection_energy_model .eq. 2) Then
	    ! Reflected hydrcarbon energy is calculated from the thermal substrate temperature.
            Call PRC ('  Reflected hydrocarbon energy calculated from'//
     > ' thermal substrate temperature.')
         ElseIf (hc_reflection_energy_model .eq. 3) Then
	    ! Reflected hydrcarbon energy calculated from Alman and Ruzic reflection data.
            Call PRC ('  Reflected hydrocarbon energy calculated from'//
     > ' Alman and Ruzic data.')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*)
     > "Error:  Option for hc_reflection_energy_"//
     > "model not supported:", hc_reflection_energy_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H48)
	 If (hc_reflection_angle_model .eq. -1) Then
	    ! Defaults to DIVIMP neutral reflectinonangle option NRFOPT.
	    Call PRI ('  HC reflection angle model option set to'//
     >        ' DIVIMP neutral reflection angle option NRFOPT: ',
     >        nrfopt)
         ElseIf (hc_reflection_angle_model .eq. 0) Then
	    ! Perfrom a reflection with a normal emittance angle.
	    Call PRC ('  Reflection angle normal to vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 1) Then
	    ! Perfrom a reflection with an specular emittance angle.
	    Call PRC ('  Reflection angle specular to vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 2) Then
	    ! Perfrom a reflection with a COS off normal emittance angle.
	    Call PRC ('  Reflection angle distributed with a COS off'//
     >        'the normal of the vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 3) Then
	    ! Perfrom a reflection with a SIN off normal emittance angle.
	    Call PRC ('  Reflection angle distributed with a SIN off'//
     >        'the normal of the vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 4) Then
	    ! Perfrom a reflection with a  SIN(SQRT) off normal emittance angle.
	    Call PRC ('  Reflection angle distributed with a SIN(SQRT'//
     >        ') off the normal of the vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 5) Then
	    ! Perfrom a reflection with a proj SIN(SQRT) off normal emittance angle.
	    Call PRC ('  Reflection angle distributed with a proj SIN('//
     >        'SQRT) off the normal of the vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 10) Then
	    ! Perfrom a reflection with a normal emittance angle.
	    Call PRC ('  Reflection angle distributed'//
     >        ' normal to the vessel wall.')
	 ElseIf (hc_reflection_angle_model .eq. 11) Then
	    ! Perfrom a reflection Janev emittance angle.
	    Call PRC ('  Reflection angle distributed'//
     >        ' normal to the vessel wall.')
         Else   
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_reflection_angle_model"//
     >        " not supported:", hc_reflection_angle_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H50)
	 If (hc_sputtering_option .eq. 0) Then
	    ! Sputtering option is off.
	    Call PRC ('  Sputtering option is off.')
	 ElseIf (hc_sputtering_option .eq. 1) Then
	    ! Sputtering option is on.
	    Call PRC ('  Sputtering option is on.')
	 Else
	    ! Error: option not supported.
	    Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_sputtering_option"//
     >        " not supported:", hc_sputtering_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H51) and (H52)
         If (hc_sticking_coef_model .eq. 0) Then
            ! Use preset sticking coefficient in all cases (typically 0.5).
            Call PRC ('  Preset global sticking coefficient is used.')
            ! (H52) dependent on H51.
            Call PRR ('  Preset neutral sticking coefficient',
     > hc_sticking_coef_preset)
         ElseIf (hc_sticking_coef_model .eq. 1) Then
            Call PRC ('  Mech, Davis, Haasz (Nucl.Mat 1997) HC launch'//
     > ' distribution used (with possibilty of launch of CH4 - C3H8).')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_sticking_coef_model"//
     > " not supported:", hc_sticking_coef_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H53)
         If (hc_sputtering_species_model .eq. 0) Then
            ! Preset sputtering table used for sputtered hydrocarbon species.
            Call PRC ('  Preset sputtering table used for sputtered'//
     >        ' hydocarbon species.')
         ElseIf (hc_sputtering_species_model .eq. 1) Then
	    ! Alman and Ruzic data used for sputtered hydrocarbon species.
            Call PRC ('  Alman and Ruzic data used for sputtered'//
     >        ' hydrocarbon species.')
         ElseIf (hc_sputtering_species_model .eq. 2) Then
	    ! Sputtered species will be equal to incoming species.
            Call PRC ('  Sputtered species will be equal to incomng'//
     >        ' hydrocarbon species.')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) "Error:  Option for "//
     >        " hc_sputtering_species_model not supported:", 
     >        hc_sputtering_species_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H54) and (H55) and (H56)
         If (hc_sputtering_energy_model .eq. 0) Then
            ! Sputtered hydrocarbon will always have preset energy set in option H55/56 (typically a few ev total).
	    Call PRR ('  Sputtered hydrocarbon from neutral impact '//
     > 'will have energy of (eV): ', hc_sput_energy_neutral_preset) ! H55.
	    Call PRR ('  Sputtered hydrocarbon from ion impact '//
     > 'will have energy of (eV): ', hc_sput_energy_ion_preset) ! H56.
         ElseIf (hc_sputtering_energy_model .eq. 1) Then
	    ! Sputtered hydrocarbon energy is the impact energy.
            Call PRC ('  Sputtered hydrocarbon energy equals'//
     > ' impact energy.')
         ElseIf (hc_sputtering_energy_model .eq. 2) Then
	    ! Sputtered hydrcarbon energy is calculated from the thermal substrate temperature.
            Call PRC ('  Sputtered hydrocarbon energy calculated from'//
     > ' thermal substrate temperature.')
         ElseIf (hc_sputtering_energy_model .eq. 3) Then
	    ! Sputtered hydrcarbon energy calculated from Alman and Ruzic reflection data.
            Call PRC ('  Sputtered hydrocarbon energy calculated from'//
     > ' Alman and Ruzic data.')
         Else
            ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_sputtering_energy_"//
     >        "model not supported:", hc_sputtering_energy_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H57)
	 If (hc_sputtering_angle_model .eq. -1) Then
	    ! Defaults to DIVIMP neutral sputtering angle option NRFOPT.
	    Call PRI ('  HC sputtering angle model option set to'//
     >        ' DIVIMP neutral sputtering angle option NRFOPT: ',
     >        nrfopt)
         ElseIf (hc_sputtering_angle_model .eq. 0) Then
	    ! Perfrom a sputter with a normal emittance angle.
	    Call PRC ('  Sputtering angle normal to vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 1) Then
	    ! Perfrom a sputter with an specular emittance angle.
	    Call PRC ('  Sputtering angle specular to vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 2) Then
	    ! Perfrom a sputter with a COS off normal emittance angle.
	    Call PRC ('  Sputtering angle distributed with a COS off'//
     >        'the normal of the vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 3) Then
	    ! Perfrom a sputter with a SIN off normal emittance angle.
	    Call PRC ('  Sputtering angle distributed with a SIN off'//
     >        'the normal of the vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 4) Then
	    ! Perfrom a sputter with a SIN(SQRT) off normal emittance angle.
	    Call PRC ('  Sputtering angle distributed with a SIN(SQRT'//
     >        ') off the normal of the vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 5) Then
	    ! Perfrom a sputter with a proj SIN(SQRT) off normal emittance angle.
	    Call PRC ('  Sputtering angle distributed with a proj'//
     >        'SIN(SQRT) off the normal of the vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 10) Then
	    ! Perfrom a sputter with a normal emittance angle.
	    Call PRC ('  Sputtering angle distributed'//
     >        ' normal to the vessel wall.')
	 ElseIf (hc_sputtering_angle_model .eq. 11) Then
	    ! Perfrom a sputter Janev emittance angle.
	    Call PRC ('  Sputtering angle distributed'//
     >        ' normal to the vessel wall.')
         Else   
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_sputtering_angle_model"//
     >        " not supported:", hc_reflection_angle_model
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! Hydrocarbon optional output.
         ! (H90)
         If (hc_coord_print_option .eq. 0) Then
            ! Decide not to print the r,z position data and hydocarbon species at each timestep.
	    Call PRC ('  hc_coord_print ascii file for hydrocarbon'//
     > ' r,z,species will not be saved.')
         ElseIf (hc_coord_print_option .eq. 1) Then
            ! Decide to print the r,z position data and hydocarbon species at each timestep.
            Call PRC ('  hc_coord_print ascii file for hydrocarbon'//
     > ' r,z,species will be saved.')
         Else   
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_coord_print_option"//
     >        " not supported:", hc_coord_print_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

         ! (H91)
         If (hc_evolve_print_option .eq. 0) Then
            ! Decide not to print the r,z position data and hydocarbon species at each hydrocarbon evolution.
	    Call PRC ('  hc_evolve_print ascii file for hydrocarbon'//
     > ' r,z,species evolution events will not be saved.')
         ElseIf (hc_evolve_print_option .eq. 1) Then
            ! Decide to print the r,z position data and hydocarbon species at each evolution event.
            Call PRC ('  hc_evolve_print ascii file for hydrocarbon'//
     > ' r,z,species evolution event will be saved.')
         Else   
	    ! Error: option not supported.
            Write (Output_Unit_HC_Alert,*) 
     >        "Error:  Option for hc_evolve_print_option"//
     >        " not supported:", hc_evolve_print_option
            Write (Output_Unit_HC_Alert,*) "Program stop."
	    Stop
	 End If

      
      End If
c
      Return
      End
! ammod end.
c
c
c
      subroutine calc_average_crossfield_dist
      implicit none
      include 'params'
      include 'cgeom'
c
c     This routine calculates the following:
c     1) For each ring in the main SOL - the average
c        poloidal-length weighted distance from the
c        separatrix for cells above and adjacent to the 
c        X-point.
c     2) The poloidal length/area of the separatrix
c        around the confined plasma and the total
c        toroidal surface area of the confined plasma. 
c
      character*80 comment
c
      real tot_spara, spara, tot_cfdist, cfdist
      real rsc,zsc
      real tot_dpara,dpara,tot_area,area
c
      integer ik,ir,nr,nc
c    
      
      WRITE(comment,'(A12,2(5x,A20))') 
     >          'Ring No.',
     >          ' CF-Dist to Sep ',
     >          ' Ring Xpt to Xpt'
      CALL PRC(comment)
c
      do ir = irsep,irwall-1
         tot_spara=0.0
         tot_cfdist=0.0 

         do ik = ikouts(1,irsep-1),ikouts(nks(irsep-1)-1,irsep-1)

            nr = korpg(ik,irsep)
            nc = korpg(ik,ir)
       
            if (nr.ne.0.and.nc.ne.0) then 
c                 
               rsc = (rvertp(4,nr) + rvertp(1,nr))/2.0
               zsc = (zvertp(4,nr) + zvertp(1,nr))/2.0

               cfdist=sqrt((rsc-rs(ik,ir))**2+(zsc-zs(ik,ir))**2)
               spara = ksb(ik,ir)-ksb(ik-1,ir)

               tot_spara = tot_spara+spara
               tot_cfdist = tot_cfdist+cfdist*spara
c
            endif
c                 
         end do
c
         if (tot_spara.gt.0.0) then
            WRITE(comment,1000) ir,tot_cfdist/tot_spara,tot_spara
            call prc(comment)
         endif
c
      end do
c
      call prb
c
c     Calculate the length and total area of the separatrix 
c
      ir = irsep-1
c
      tot_dpara = 0.0

c     IPP/08 Krieger - tot_area was not initialized
      tot_area = 0.0

      do ik = 1,nks(ir)-1
c
         nr = korpg(ik,ir)

         if (nr.ne.0) then 
c                 
            dpara = sqrt((rvertp(2,nr)-rvertp(3,nr))**2+
     >                   (zvertp(2,nr)-zvertp(3,nr))**2)
c
            tot_dpara = tot_dpara + dpara
c
            area = dpara * 2.0 * PI * (rvertp(2,nr)+rvertp(3,nr))/2.0
c
            tot_area=tot_area + area
c
         endif
c
      end do
c
      write(comment,'(a,f15.4)') 'Total poloidal length'//
     >                           ' of core boundary (m):',
     >                                 tot_dpara
      call prc(comment)
      write(comment,'(a,f15.4)') 'Total toroidal surfac'//
     >                           'e area of Core (m2)  :',
     >                                 tot_area
      call prc(comment)
      call prb





1000  FORMAT(I12,5x,F15.6,10x,f15.6)
      return
      end 
      
