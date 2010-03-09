c     -*-Fortran-*-
c
      SUBROUTINE MONPRI (FACT,VFLUID,NIZS,SDTZS,sdtzs2,
     >            STOTS,DOUTS,RIONS,TDEP,TWALL,DPARAS,DCROSS,
     >            TNTOTS,FPTARG,acttarg,SPTOTS,sitots,coreouts,
     >            e2dtots,e2dptots,cre2d,cre2dizs,
     >            impurity_content)
      implicit none 
      include    'params'
      include    'cadas'
      include    'commv'
      include    'cgeom'
      include    'comtor'
      include    'dynam1'
      include    'dynam3' 
      include    'promptdep' 
c
      include    'fperiph_com' 
c
      INTEGER   NIZS,cre2d,cre2dizs
      REAL      FACT,VFLUID,STOTS(46),TDEP,TWALL
      real      sptots(10,0:maxizs+2),sitots(maxizs+2,2)
      REAL      TNTOTS(0:MAXIZS+1,2),FPTARG(MAXIZS),acttarg(maxizs)
      REAL      RIONS(MAXIZS),DCROSS(4)
      real      SDTZS(-1:MAXIZS),sdtzs2(-1:maxizs)
      real      e2dtots(9),e2dptots(6,0:maxe2dizs)
      real      impurity_content(0:maxizs,4,2)
c
      DOUBLE PRECISION DOUTS(MAXIZS,9),DPARAS(MAXIZS,6)
      DOUBLE PRECISION coreOUTS(MAXIZS,9)
C
C***********************************************************************
C
C       PRINTS MONITORING VARIABLES STORED IN COMMON BLOCK COMMV
C
C         ARGUMENTS -
C     QTIM   : TIMESTEP
C     FACT   : RECIPROCAL OF THE TOTAL NUMBER OF IONS INJECTED
C     VFLUID : FLOW VELOCITY AT Y=0
C     NIZS   : NUMBER OF IONIZATION STATES BEING FOLLOWED
C     SDTZS  : CLOUD TEMPERATURES
C     STOTS  : ZEFF AVERAGES, TOTAL POWER, ETC, ETC ...
C
C***********************************************************************
C
C
      INTEGER   IZ,IM,II,IN,i,j,id
      REAL      RCR,RCAB,RSAB,TEXT
      REAL      RTR,RTAB,RTAV,RAVA,RMACH,RAVMCH,RENEGY
      REAL      RAVEGY,RZ0,RZ1,RZ2,DRZ,RTBS,FORCE
      REAL      VEXIT,MTSO,CICOUT,mtsc
      real      entexdat(3,3),mfact
c
      real      tmpval
      real      carea,main_content,div_content,div_reten
      real      totsum
c     IPP/08 Krieger - extended comment string length
      character*256 comment 
c
c     Core and divertor content calculations - local variables  
c
      integer ik,ir
      real neutc_all,neutc_prim,neutc_area,ionc_all,ionc_area
      real neutc_temp,neutc_ave_temp
      real neutc_div_density,neutc_divedge_density
      real ionc_core_density,ionc_coreedge_density
      real tot_incore,tot_inpfz      

C
C-----------------------------------------------------------------------
C          INITIALISATION
C-----------------------------------------------------------------------
C
c     IPP/08 Krieger - RTBS was not initialized
      RTBS = 0.0
      RCR  = 0.0
      RTR  = 0.0
      RCAB = 0.0
      RSAB = 0.0
      RTAB = 0.0
      RTAV = 0.0
      RAVA = 0.0
      RMACH = 0.0
      RENEGY = 0.0
      RAVMCH = 0.0
      RAVEGY = 0.0
      RZ0  = 0.0
      RZ1  = 0.0
      RZ2  = 0.0
      DRZ  = 0.0
      CICOUT = 0.0
c
c
C-----------------------------------------------------------------------
C
      call prb
      call prchtml('--- SUMMARY BY CHARGE STATE ---',
     >             'pr_chargestate','0','B')  
      call prb  
C
       DO  100  IZ = 1, NIZS
C
C------- CALCULATE MEAN TIME SPENT OUTBOARD IN EACH CHARGE STATE.
C
         IF (RIONS(IZ).GT.0.0) THEN
           MTSO = QTIM * SNGL(DOUTS(IZ,1)) / RIONS(IZ)
           mtsc = QTIM * SNGL(coreOUTS(IZ,1)) / RIONS(IZ)
         ELSE
           MTSO = 0.0
           mtsc = 0.0
         ENDIF
         CICOUT = CICOUT + QTIM * SNGL(DOUTS(IZ,1)) * FACT
C
C-----------------------------------------------------------------------
C        PRINT FIRST LOT OF DETAILS
C-----------------------------------------------------------------------
C
         CALL PRB
         WRITE (7,'('' ***     IONIZATION   STATE '',I3,6X,''***'')') IZ
         CALL PRB
         CALL PRR ('CLOUD TEMPERATURE  (EV)                      ',
     >                                 SDTZS(IZ))
         CALL PRR ('MEAN VELOCITY BASED CLOUD TEMPERATURE  (EV)  ',
     >                                 SDTZS2(IZ))
         CALL PRR ('MEAN TIME SPENT OUTBOARD (S)                 ',
     >                                 MTSO)
         CALL PRR ('MEAN TIME SPENT IN CORE  (S)                 ',
     >                                 MTSC)
c
         if (citizs(iz).eq.0.0) then
            tmpval = 0.0
         else
            tmpval = qtim * cieizs(iz)/citizs(iz)
         endif   
         call prr ('MEAN ELAPSED TIME SPENT IN STATE (S)         ',
     >                 tmpval )
c
         CALL PRR0('NUMBER OF IONS REACHING THIS CHARGE STATE  ',
     >                                  RIONS(IZ) )
c
         call prc ('Radiation statistics (normalized to ion launch)')
         call prr ('  Power radiated from MAIN    :  ',sptots(1,iz))
         call prr ('  Line radiation from MAIN    :  ',sptots(2,iz))
         call prr ('  Power radiated from SOL     :  ',sptots(3,iz))
         call prr ('  -Divertor Region "below" Xpt:  ',sptots(7,iz))  
         call prr ('  -Main Region     "above" Xpt:  ',sptots(9,iz))  
         call prr ('  Line radiation from SOL     :  ',sptots(4,iz))
         call prr ('  -Divertor Region "below" Xpt:  ',sptots(8,iz))  
         call prr ('  -Main Region     "above" Xpt:  ',sptots(10,iz))  
         call prr ('  Power radiated from TRAP    :  ',sptots(5,iz))
         call prr ('  Line radiation from TRAP    :  ',sptots(6,iz))
         call prr ('  Total power radiated        :  ',sptots(1,iz)
     >               + sptots(3,iz)+ sptots(5,iz))
         call prr ('  Total line radiation        :  ',sptots(2,iz)
     >               + sptots(4,iz)+sptots(6,iz))
c       
         if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) then 
c
         call prb
         call prc ('FOLLOWING BASED ON FLUID CODE IMPURITY DATA:')
         call prc ('Radiation data (absolute) E2D/DIVIMP')
         call prr2 ('  Power radiated from MAIN    :  ',e2dptots(1,iz),
     >                            sptots(1,iz) * absfac)
         call prr2 ('  Line radiation from MAIN    :  ',e2dptots(2,iz),
     >                            sptots(2,iz) * absfac)
         call prr2 ('  Power radiated from SOL     :  ',e2dptots(3,iz),
     >                            sptots(3,iz) * absfac)
         call prr2 ('  Line radiation from SOL     :  ',e2dptots(4,iz),
     >                            sptots(4,iz) * absfac)
         call prr2 ('  Power radiated from TRAP    :  ',e2dptots(5,iz),
     >                            sptots(5,iz) * absfac)
         call prr2 ('  Line radiation from TRAP    :  ',e2dptots(6,iz),
     >                            sptots(6,iz) * absfac)
         call prr2 ('  Total EDGE2D power radiated :  ',e2dptots(1,iz)
     >        + e2dptots(3,iz)+e2dptots(5,iz),
     >        absfac*(sptots(1,iz)+sptots(3,iz)+sptots(5,iz)))
         call prr2 ('  Total EDGE2D line radiation :  ',e2dptots(2,iz)
     >        + e2dptots(4,iz)+e2dptots(6,iz),
     >        absfac*(sptots(2,iz)+sptots(4,iz)+sptots(6,iz)))
c
         call prb
c
         endif
c
         IF (CICUTS(IZ).GT.0.0) THEN
           CALL PRB
           CALL PRI ('NUMBER OF IONS STILL IN PLASMA AT CUTOFF     ',
     >                             NINT(CICUTS(IZ)))
           CALL PRR ('  MEAN TEMPERATURE  (EV)                     ',
     >                             (CRTRCS(IZ)/CICUTS(IZ)))
         ENDIF
         IF (RIONS(IZ).GT.0.0) THEN
           CALL PRR ('AVERAGE MINIMUM Z VALUE REACHED              ',
     >                              CXXX(IZ)/RIONS(IZ))
           CALL PRR ('AVERAGE MAXIMUM S OR SMAX-S VALUE REACHED    ',
     >                              CSSS(IZ)/RIONS(IZ))
      CALL PRI ('NO OF EXITING IONS ORIG. IONIZ IN MAIN PLASMA      ',
     >                                    NINT(CMMM(IZ)))
           IF (CMMM(IZ).GT.0.0) THEN
      CALL PRR ('AVERAGE Z VALUE FOR THESE EXITS                    ',
     >                            CMMMX(IZ)/MAX(LO,CMMM(IZ)))
      CALL PRR ('AVERAGE S OR SMAX-S FOR THESE EXITS                ',
     >                            CMMMS(IZ)/MAX(LO,CMMM(IZ)))
           ENDIF
           CALL PRI2 ('NO OF OTHER IONS ENTERING / EXITING MAIN',
     >                    NINT(CNNN(IZ)), NINT(CLLL(IZ)))
           IF (CNNN(IZ).GT.0.0.OR.CLLL(IZ).GT.0.0) THEN
           CALL PRR2 ('AVERAGE Z VALUE FOR THESE ENTRIES/EXITS ',
     >        CNNNX(IZ)/MAX(LO,CNNN(IZ)), CLLLX(IZ)/MAX(LO,CLLL(IZ)))
           CALL PRR2 ('AVERAGE S OR SMAX-S FOR ENTRIES/EXITS   ',
     >        CNNNS(IZ)/MAX(LO,CNNN(IZ)), CLLLS(IZ)/MAX(LO,CLLL(IZ)))
           ENDIF
           IF (CNNN(IZ).GT.0.0) THEN
      WRITE(7,'(1X,A,F11.6)')'AV. ORIG. R AT NEUT IONIZ. FOR ENTRIES  ',
     >                            Cnorgr(IZ)/MAX(LO,CNNN(IZ))
      WRITE(7,'(1X,A,F11.6)')'AV. ORIG. Z AT NEUT IONIZ. FOR ENTRIES  ',
     >                            Cnorgz(IZ)/MAX(LO,CNNN(IZ))
      WRITE(7,'(1X,A,F11.6)')'AV. ORIG. |S| AT NEUT IONIZ. FOR ENTRIES',
     >                            CNorgs(IZ)/MAX(LO,CNNN(IZ))
      WRITE(7,'(1X,A,F11.6)')'AV. ORIG. K AT NEUT IONIZ. FOR ENTRIES  ',
     >                            CNNNK(IZ)/MAX(LO,CNNN(IZ))
           CALL PRR ('AV. TIME FROM NEUT IONIZ. TO ENTRY      ',
     >                            CNNNT(IZ)/MAX(LO,CNNN(IZ)))
      WRITE(7,'(1X,A,F11.6)')'AV. K OVER TIME FROM ORIGINAL IONIZ     ',
     >                            CNNNKT(IZ)/MAX(LO,CNNN(IZ))
           ENDIF
         ENDIF
C
C-----------------------------------------------------------------------
C        PRINT ABSORPTION DETAILS FOR TARGET
C-----------------------------------------------------------------------
C
         CALL PRB
         CALL PRR0('TOTAL NUMBER ABSORBED ON TARGET & WALLS    ',
     >                                          CICABS(IZ) )
         IF (FPOPT.eq.3) THEN
c           CALL PRR ('NUMBERS ABSORBED ON ACTUAL TARGET          ',
c     >         CICABS(IZ)-FPTARG(IZ)-TNTOTS(IZ,1)-TNTOTS(IZ,2) )
           CALL PRR ('NUMBERS ABSORBED ON ACTUAL TARGET          ',
     >               acttarg(iz) )
           CALL PRR ('NUMBERS ABSORBED ON FP TARGET              ',
     >              FPTARG(IZ))
         ELSE
           CALL PRR ('NUMBERS ABSORBED ON TARGET PLATES          ',
     >                acttarg(iz))
         ENDIF
         CALL PRR ('NUMBERS ABSORBED ON OUTER WALL             ',
     >               TNTOTS(IZ,1))
         CALL PRR ('NUMBERS ABSORBED ON TRAP WALL              ',
     >               TNTOTS(IZ,2))
C
         IF (CICABS(IZ).GT.0.0) THEN
           CICABS(IZ) = MAX (LO,CICABS(IZ))
           CALL PRR  ('  TIME FIRST ION ABSORBED  (S)               ',
     >                                    QTIM*CIFABS(IZ))
           CALL PRR  ('  TIME LAST ION ABSORBED  (S)                ',
     >                                    QTIM*CILABS(IZ))
           CALL PRR  ('  MEAN ABSORPTION TIME  (S)                  ',
     >                                    QTIM*CISABS(IZ)/CICABS(IZ))
           CALL PRR  ('  MEAN TEMPERATURE AT ABSORPTION (EV)        ',
     >                                    CRTABS(IZ)/CICABS(IZ))
           CALL PRR  ('  MEAN VELOCITY AT ABSORPTION (M/S)          ',
     >                                    CRVABS(IZ)/(CICABS(IZ)*QTIM))
           CALL PRR  ('  MEAN ABS(VEL) AT ABSORPTION (M/S)          ',
     >                                    CRAVAV(IZ)/(CICABS(IZ)*QTIM))
C
           IF (VFLUID.LE.0.0 .OR. CRAVAV(IZ).LE.0.0) THEN
             RMACH = 0.0
             DRZ   = 0.0
             VEXIT = 0.0
           ELSE
             RMACH = CRAVAV(IZ) / (CICABS(IZ)*QTIM) / VFLUID
             DRZ   = CICABS(IZ) / RMACH
             VEXIT = CRAVAV(IZ) / (CICABS(IZ)*QTIM)
           ENDIF
           RENEGY = 3.0 * REAL(IZ) * CTBS(IZ) / CICABS(IZ) +
     >              5.22E-9 * CRMI * VEXIT * VEXIT +
     >              2.0 * CRTABS(IZ) / CICABS(IZ)
C
           CALL PRR  ('  IMPURITY MACH NUMBER AT ABSORPTION         ',
     >                                    RMACH)
           CALL PRR  ('  IMPACT ENERGY AT ABSORPTION (EV)           ',
     >                                    RENEGY)
c slmod begin
c DUMP WALL IMPURITY FLUX DATA HERE...     
c slmod end
         ENDIF
C
C-----------------------------------------------------------------------
C        PRINT IONISATION & RECOMBINATION DETAILS
C-----------------------------------------------------------------------
C
         CALL PRB
         CALL PRC ('NUMBER OF         (FROM T=0)  FIRST TIME  LAST TIME
     > MEAN TIME')
         IF (CICIZS(IZ).GT.0.0) THEN
           WRITE (7,9001) '  IONIZATIONS   ',     CICIZS(IZ) ,
     >       QTIM*CIFIZS(IZ),QTIM*CILIZS(IZ),QTIM*CISIZS(IZ)/CICIZS(IZ)
         ELSE
           WRITE (7,9001) '  IONIZATIONS         0'
         ENDIF
         IF (CICRCS(IZ).GT.0.0) THEN
           WRITE (7,9001) '  RECOMBINATIONS',     CICRCS(IZ) ,
     >       QTIM*CIFRCS(IZ),QTIM*CILRCS(IZ),QTIM*CISRCS(IZ)/CICRCS(IZ)
         ELSE
           WRITE (7,9001) '  RECOMBINATIONS      0'
         ENDIF
 9001 FORMAT(1X,A,1X,F9.2:,1P,3X,3(2X,G9.2))
C
C-----------------------------------------------------------------------
C        PRINT ION REMOVAL LOSS DETAILS
C-----------------------------------------------------------------------
C
         CALL PRB
         CALL PRR0('NUMBERS OF IONS LOST BY REMOVAL            ',
     >                   CICLOS(IZ)  )
C
         IF( CICLOS(IZ).GT.0.0 ) THEN
             CALL PRR ('  TIME FIRST ION REMOVED (S)                 ',
     >                 QTIM*CIFLOS(IZ) )
             CALL PRR ('  TIME LAST  ION REMOVED (S)                 ',
     >                 QTIM*CILLOS(IZ) )
             CALL PRR ('  MEAN LOSS TIME         (S)                 ',
     >                 QTIM*CISLOS(IZ)/CICLOS(IZ) )
             CALL PRR ('  AVERAGE S OR SMAX-S VALUE REACHED (M)      ',
     >                 CSSSS(IZ)/CICLOS(IZ) )
         END IF
C
C-----------------------------------------------------------------------
C        UPDATE TOTALS
C-----------------------------------------------------------------------
C
         RCR  = RCR  + CICUTS(IZ)
         RTR  = RTR  + CRTRCS(IZ)
C
         RCAB = RCAB + CICABS(IZ)
         RSAB = RSAB + CISABS(IZ)
         RTAB = RTAB + CRTABS(IZ)
         RTAV = RTAV + CRVABS(IZ)
         RAVA = RAVA + CRAVAV(IZ)
         RTBS = RTBS + CTBS(IZ)
         RAVMCH = RAVMCH +  RMACH*CICABS(IZ)
         RAVEGY = RAVEGY + RENEGY*CICABS(IZ)
         RZ0  = RZ0 + DRZ
         RZ1  = RZ1 + DRZ*FLOAT(IZ)
         RZ2  = RZ2 + DRZ*FLOAT(IZ)*FLOAT(IZ)
  100  CONTINUE
C
C-----------------------------------------------------------------------
C     PRINT SUMMARY DETAILS OVER ALL IONISATION STATES
C-----------------------------------------------------------------------
C
      CALL PRB
c
c      CALL PRC ('***   ALL   IONIZATION   STATES     ***')
c
      CALL PRChtml ('***   ALL   IONIZATION   STATES     ***',
     >              'pr_allstates','0','B')
      CALL PRB
      CALL PRR ('FRACTION OF IONS STILL IN PLASMA AT CUTOFF   ',
     >                             (FACT*RCR))
      IF (RCR .GT. 0.0)
     >   CALL PRR ('  MEAN TEMPERATURE  (EV)                     ',
     >                             (RTR/RCR))
      CALL PRB
      CALL PRR ('FRACTION OF IONS IONIZED BEYOND LIMIT        ',
     >                             (FACT*CICIZS(NIZS)))
      CALL PRB
c
      if (checkleak) then
         call prr ('Total number of ions exceeding maximum specified S',
     >              cleakn(cleaksn,nizs+1))
         call prr ('Maximum value of S specified for leakage (m)      ',
     >             cleaks(cleaksn))
         call prb
         call prc ('Detailed Summary of leakage: ')
c
         do 180 in = 1,cleaksn
c
           if (cleakn(in,nizs+1).gt.0.0) then
             call prr('  Leakage by charge state:  S  >  ',cleaks(in))
             do 190 iz = 1,nizs
               write(7,'(''     Ionization state: '',
     >            I4,''  Leakage:    '',G12.5)' ) iz,cleakn(in,iz)
 190         continue
             if (in.eq.cleaksn)
     >         call prr('Mean time needed to exceed S       ',
     >              cleakt/cleakn(in,nizs+1))
           endif
 180     continue
         call prb
      endif
c
      CALL PRR ('FRACTION OF IONS ABSORBED                    ',
     >                             (FACT*RCAB))
      IF (RCAB .GT. 0.0)  THEN
         CALL PRR ('  MEAN ABSORPTION TIME  (S)                  ',
     >                             ((QTIM*RSAB)/RCAB))
         CALL PRR ('  MEAN ION TEMPERATURE AT ABSORPTION (EV)    ',
     >                             (RTAB/RCAB))
         CALL PRR ('  MEAN PLASMA TEMPERATURE AT ABSORPTION (EV) ',
     >                             (RTBS/RCAB))
         CALL PRR ('  MEAN VELOCITY AT ABSORPTION (M/S)          ',
     >                             (RTAV/(RCAB*QTIM)))
         CALL PRR ('  MEAN ABSOLUTE VELOCITY AT ABSORPTION (M/S) ',
     >                             (RAVA/(RCAB*QTIM)))
         CALL PRR ('  MEAN WEIGHTED MACH NUMBER AT ABSORPTION    ',
     >                            (RAVMCH/RCAB))
         CALL PRR ('  MEAN WEIGHTED IMPACT ENERGY (EV)           ',
     >                            (RAVEGY/RCAB))
         CALL PRR ('  EXIT IMPURITY DENSITY                      ',
     >                            (RZ0/RCAB))
         CALL PRR ('  EXIT IMPURITY  Z                           ',
     >                            (RZ1/RCAB))
         CALL PRR ('  EXIT IMPURITY  Z**2                        ',
     >                            (RZ2/RCAB))
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED IN EACH STATE')
         WRITE (7,'((3X,6(I3,F7.4)))')
     >     (IZ,  CICABS(IZ)/RCAB, IZ=1,NIZS)
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED ON')
         WRITE (7,'(3X,2(2X,A7,1P,G9.2))')
     >            'TARGET ',TDEP/RCAB,'WALLS  ',TWALL/RCAB
C
         TEXT = 0.0
         DO 200 IM = 1, 10
           TEXT = TEXT + CTEXS(IM)
  200    CONTINUE
         CALL PRC ('  RATIO OF ABSORPTIONS BY EXIT TEMPS AS MULTIPLES OF
     > TB(0)')
         WRITE (7,'(3X,9F6.2,'' HIGHER'')') (REAL(IM)*0.2,IM=1,9)
         WRITE (7,'(4X,10F6.3)') (CTEXS(IM)/TEXT,IM=1,10)
      ENDIF
c
C-----------------------------------------------------------------------
c 
c     PROMPT REDEPOSITION SUMMARY 
c
c     If the prompt deposition option was ON -print a summary.
c
C-----------------------------------------------------------------------
c 
      if (prompt_depopt.eq.1.or.prompt_depopt.eq.2) then   
c
         totsum = 0.0 
         do id = 1,nds
            totsum = totsum + promptdeps(id,1)
c
            if (promptdeps(id,1).gt.0.0) then  
               promptdeps(id,2) = promptdeps(id,2) / promptdeps(id,1) 
               promptdeps(id,3) = promptdeps(id,3) / promptdeps(id,1) 
               promptdeps(id,4) = promptdeps(id,4) / promptdeps(id,1) 
               promptdeps(id,5) = promptdeps(id,5) / promptdeps(id,1) 
            endif 
c
         end do
c 
c
C-----------------------------------------------------------------------
C

         call prb
c
c         call prc('   PROMPT ION DEPOSITION WAS TURNED ON -')
c
         call prchtml('   PROMPT ION DEPOSITION SUMMARY',
     >                'pr_promptdep','0','B')
         CALL PRR('   - NUMBER OF PARTICLES DEPOSITED BY PROMPT'//
     >                     ' DEPOSITION: ',TOTSUM)
c
         IF (TOTSUM.GT.0.0) THEN 
            CALL PRC ('   - ANALYSIS BY TARGET SEGMENT :')
            call prc ('   TARGET   WEIGHT       SELF      AV-SH '//
     >                '       IMPACT'//
     >                '   AVERAGE DIST     AVERAGE      SHEATH')
            call prc ('     IND    REDEP        SPUT      V-DROP'//
     >                '       ENERGY'//
     >                '     AT IONIZ    LARMOR RADIUS   THICKNESS')
            DO ID = 1,NDS
c
               i = irds(id)
c
               if (id.le.ndsin) then 
                  j = 1
               else
                  j = 2
               endif   
c
               IF (promptdeps(id,1).gt.0.0) then   
c                 IPP/08 Krieger - changed formatting
                  WRITE (COMMENT,'(3X,I4,1X,g12.5,1x,g12.3,
     >                  1x,g12.5,1x,g12.5,
     >                  3(1x,g12.5))') id,promptdeps(id,1),
c                  WRITE (COMMENT,'(3X,I4,1X,2(1x,F7.2),1x,f6.2,1x,
c     >                  f9.2,
c     >                  3(1x,g12.5))') id,promptdeps(id,1),
     >                  promptdeps(id,6),
     >                  promptdeps(id,4),promptdeps(id,5),
     >                  promptdeps(id,2), 
     >                  promptdeps(id,3),mps_thickness(i,j)
                  call prc(comment)
c
               endif
            end do
         endif 
c
      endif


C
C-----------------------------------------------------------------------
C     PRINT CORE LEAKAGE SECTION
C-----------------------------------------------------------------------
C

      call prb 
      call prchtml('--- ANALYSIS OF CORE LEAKAGE ---',
     >           'pr_leakage','0','B')  
      call prb

      CALL PRR ('FRACTION OF IONS ENTERING MAIN PLASMA    ',
     >                             (FACT*CICRIN))
c
c     Calculate totals of all ions entering/exiting main and 
c     print the statistics.
c
      do i = 1,3
         do j = 1,3
            entexdat(i,j) = 0.0 
         end do 
      end do
   
      do iz = 1,nizs 
         entexdat(1,1) = entexdat(1,1) + cmmm(iz)
         entexdat(1,2) = entexdat(1,2) + cmmmx(iz)
         entexdat(1,3) = entexdat(1,3) + cmmms(iz)
         entexdat(2,1) = entexdat(2,1) + cnnn(iz)
         entexdat(2,2) = entexdat(2,2) + cnnnx(iz)
         entexdat(2,3) = entexdat(2,3) + cnnns(iz)
         entexdat(3,1) = entexdat(3,1) + clll(iz)
         entexdat(3,2) = entexdat(3,2) + clllx(iz)
         entexdat(3,3) = entexdat(3,3) + cllls(iz)
      end do
c
      CALL PRI ('  TOTAL IONS ORIG. IONIZ IN MAIN PLASMA      ',
     >                                    NINT(entexdat(1,1)))
      IF (entexdat(1,1).GT.0.0) THEN
        CALL PRR ('  AVERAGE Z VALUE FOR THESE EXITS            ',
     >                entexdat(1,2)/entexdat(1,1))
        CALL PRR ('  AVERAGE S OR SMAX-S FOR THESE EXITS        ',
     >                entexdat(1,3)/entexdat(1,1))
      ENDIF
      CALL PRI2('  TOTAL OF OTHER IONS ENTERING / EXITING MAIN',
     >                   nint(entexdat(2,1)), nint(entexdat(3,1)))
      IF (entexdat(2,1).GT.0.0.or.entexdat(3,1).GT.0.0) THEN
        CALL PRR2 ('  AVERAGE Z VALUE FOR THESE ENTRIES/EXITS   ',
     >         entexdat(2,2)/max(lo,entexdat(2,1)), 
     >         entexdat(3,2)/max(lo,entexdat(3,1)))
        CALL PRR2 ('  AVERAGE S OR SMAX-S FOR ENTRIES/EXITS     ',
     >         entexdat(2,3)/max(lo,entexdat(2,1)), 
     >         entexdat(3,3)/max(lo,entexdat(3,1)))
      ENDIF
c
c     Sources of particles entering main
c       
      call prb 
      call prc('ORIGINAL SOURCE OF PARTICLES ENTERING MAIN:')
c
c      call prchtml('ORIGINAL SOURCE OF PARTICLES ENTERING MAIN',
c     >             'pr_addleak','0','B')
c
      call prr ('  NUMBER FROM REGULAR LAUNCH (NO REFLECTIONS):',
     >          cvvnrfm)     
      call prr ('  NUMBER FROM FP      LAUNCH (NO REFLECTIONS):',
     >          cvvfpnrf)     
      if (nrfopt.eq.1.or.nrfopt.eq.2) then 
      call prr ('  NUMBER FROM REGULAR LAUNCH (REFLECTION >=1):',
     >          cvvrefm)     
      call prr ('  NUMBER FROM FP      LAUNCH (REFLECTION >=1):',
     >          cvvfpref)     
      endif
      call prb
c
      IF (CICRIN .GT. 0.0) THEN
         CALL PRR ('  TIME FIRST ION ENTERED MAIN PLASMA (S) ',
     >                             (QTIM*CIFRIN))
         CALL PRR ('  MEAN TIME TO ENTER MAIN PLASMA (S)     ',
     >                             (QTIM*CISRIN/CICRIN))
         CALL PRR ('  MEAN ION TEMP AT MAIN PLASMA ENTRY (EV)',
     >                             (CITRIN/CICRIN))
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. ORIG. K AT NEUT IONIZ FOR ENTRIES   ',CIKRIN/CICRIN
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. K OVER TIME FROM ORIGINAL IONIZ     ',CKTRIN/CICRIN
         CALL PRR ('  AV. R VALUE AT CREATION OF ION ENTRIES  ',
     >                             CVVRM/CICRIN)
         CALL PRR ('  AV. Z VALUE AT CREATION OF ION ENTRIES  ',
     >                             CVVZM/CICRIN)
         CALL PRR ('  AV. K VALUE AT CREATION OF ION ENTRIES  ',
     >                             CVVKM/CICRIN)
         CALL PRR ('  AV. |S| VALUE AT CREATION OF ION ENTRIES',
     >                             CVVSM/CICRIN)
         CALL PRR ('  AVERAGE Z VALUE AT ENTRY               ',
     >                             CVVXE/CICRIN)
         CALL PRR ('  AVERAGE Z TRIP ON PERIPHERAL K''S      ',
     >                             CVVXP/CICRIN)
         CALL PRR ('  AVERAGE Z TRIP ON K''S NEAR SEPTRIX    ',
     >                             CVVXS/CICRIN)
      ENDIF
      CALL PRR ('FRACTION OF IONS NOT ENTERING MAIN PLASMA',
     >                             (FACT*CICRNO))
      IF (CICRNO .GT. 0.0) THEN
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. ORIG. R AT NEUT IONIZ FOR NON-ENTS ',CIRRNO/CICRNO
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. ORIG. Z AT NEUT IONIZ FOR NON-ENTS ',CIZRNO/CICRNO
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. ORIG. K AT NEUT IONIZ FOR NON-ENTS ',CIKRNO/CICRNO
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. ORIG.|S| AT NEUT IONIZ FOR NON-ENTS',CISRNO/CICRNO
         WRITE (7,'(1X,A,F11.6)')
     >     '  AV. K OVER TIME FROM ORIGINAL IONIZ    ',CKTRNO/CICRNO
      ENDIF
      CALL PRR ('FRACTION OF IONS STARTING IN MAIN PLASMA ',
     >                             (FACT*CICRNJ))
      CALL PRR ('FRACTION OF IONS REACHING CENTRAL MIRROR ',
     >                             (FACT*CICRXA))
      IF (CICRXA .GT. 0.0) THEN
         CALL PRR ('  TIME FIRST ION REACHED MIRROR (S)      ',
     >                             (QTIM*CIFRXA))
         CALL PRR ('  MEAN TIME TO REACH MIRROR (S)          ',
     >                             (QTIM*CISRXA/CICRXA))
         CALL PRR ('  MEAN ION TEMPERATURE AT MIRROR (EV)    ',
     >                             (CITRXA/CICRXA))
      ENDIF
      CALL PRR ('MEAN TIME SPENT OUTBOARD BY EACH ION (S) ',
     >                               CICOUT)
      IF (CICOUT .GT. 0.0)
     >   CALL PRR ('DEEPEST OUTBOARD PENETRATION  KMAX       ',
     >                                 CKKMAX)
      IF (CICRIN .GT. 0.0)
     >   CALL PRR ('FARTHEST INBOARD PENETRATION  KMIN       ',
     >                                 CKKMIN)

c
c     Print out the summary of average S position and field line
c     of ionization for all neutrals from all target elements.   
c
c
 
      call prleakage


C
C-----------------------------------------------------------------------
C     PRINT "OTHER INFORMATION" SECTION
C-----------------------------------------------------------------------
C
      CALL PRB
      call prchtml('--- ADDITIONAL CALCULATED VALUES ---',
     >             'pr_other','0','B')
      call prb

      CALL PRI ('MEAN NUMBER OF COLLISIONS PER ION        ',
     >                             NINT(FACT*CICCOL))
      CALL PRR ('AVERAGE Z VALUE AT CREATION              ',
     >                             FACT*CVVXC)
c
      if (cgridopt.eq.0) then
       CALL PRC ('"NEAR TARGET" MEANS Z > ZXP, |R-RXP| < 0.4')
      elseif (cgridopt.eq.1.or.cgridopt.eq.2.or.cgridopt.eq.3) then
       CALL PRC ('"NEAR TARGET" MEANS Z < ZXP, |R-RXP| < 0.4')
      endif
c
      CALL PRR ('AREA OF "NEAR TARGET" REGION     (M**2)  ',STOTS(9))
      CALL PRR ('MEAN NIE NEAR TARGET             (M**-3) ',STOTS(10))
      CALL PRR ('MEAN ZB.NBT NEAR TARGET          (M**-3) ',STOTS(11))
      CALL PRR ('MEAN ZEFF NEAR TARGET                    ',STOTS(12))
      CALL PRR ('AREA OF SOL+TRAP REGION          (M**2)  ',STOTS(15))
      CALL PRR ('TOTAL PLASMA   IN SOL+TRAP               ',STOTS(16))
      CALL PRR ('TOTAL IMPURITY IN SOL                    ',STOTS(6))
      CALL PRR ('TOTAL IMPURITY IN TRAP                   ',STOTS(37))
      CALL PRR ('TOTAL IMPURITY IN SOL+TRAP               ',STOTS(37)+
     >                                                      stots(6))
      if (cgridopt.eq.0) then
       CALL PRR ('TOTAL IMPURITY  NEAR DIVERTOR Z > ZXP   ',STOTS(21))
       CALL PRR ('TOTAL GRID AREA NEAR DIVERTOR Z > ZXP   ',STOTS(44))
       CALL PRR ('AVERAGE IMPURITY DENSITY      Z > ZXP   ',STOTS(21)/
     >                          stots(44))
       IF (CZD.NE.99.0) THEN
        CALL PRR ('TOTAL IMPURITY IN DIVERTOR   Z > ZD     ',STOTS(22))
        CALL PRR ('                                 ZD   = ',CZD)
       ENDIF
      elseif (cgridopt.eq.1.or.cgridopt.eq.2.or.cgridopt.eq.3) then
       CALL PRR ('TOTAL IMPURITY  NEAR DIVERTOR Z < ZXP   ',STOTS(21))
       CALL PRR ('TOTAL GRID AREA NEAR DIVERTOR Z < ZXP   ',STOTS(44))
       CALL PRR ('AVERAGE IMPURITY DENSITY      Z < ZXP   ',STOTS(21)/
     >                          stots(44))
       IF (CZD.NE.99.0) THEN
        CALL PRR ('TOTAL IMPURITY IN DIVERTOR   Z < ZD     ',STOTS(22))
        CALL PRR ('                                 ZD   = ',CZD)
       ENDIF
      endif
c
      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) 
     >CALL PRR2 ('TOTAL *FLUID CODE* IMPURITY IN SOL      (E/D)',
     >                                           e2dTOTS(4),
     >                                     stots(6) * absfac)
      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) 
     >CALL PRR2 ('TOTAL *FLUID CODE* IMPURITY IN TRAP     (E/D)',
     >                                           e2dTOTS(7),
     >                                    stots(37) * absfac)
c
      CALL PRR ('AREA OF MAIN PLASMA, EXC MIRROR  (M**2)  ',STOTS(1))
      CALL PRR ('TOTAL PLASMA   IN MAIN     NBTOT         ',STOTS(2))
      CALL PRR ('TOTAL IMPURITY IN MAIN     NITOT         ',STOTS(3))
c
      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) 
     >CALL PRR2 ('TOTAL *FLUID CODE* IMPURITY IN MAIN     (E/D)',
     >                                               e2dTOTS(1),
     >                                       stots(3) * absfac)
c
      CALL PRR ('PLASMA IN MAIN / MAIN AREA NBAVG         ',STOTS(13))
      CALL PRR ('IMPURITY IN MAIN/MAIN AREA NIAVG         ',STOTS(17))
      call prc ('IMPURITY DENSITY JUST INSIDE SEPARATRIX: ')
      do 300 iz = 1,nizs
        write(7,'(''     IONIZATION STATE: '',
     >            I4,''  NI        '',G12.5)' ) iz,sitots(iz,1)
 300  continue
      call prr ('    ALL STATES                           ',
     >          sitots(maxizs+1,1))
c
c     For top half of core - if core marfe option is on
c
      if (ccoreopt.eq.2.or.ccoreopt.eq.3.or.ccoreopt.eq.4
     >    .or.ccoreopt.eq.5.or.ccoreopt.eq.6) then 

         CALL PRR ('IMPURITY IN MAIN(TOP)/MAIN(TOP)AREA NIAVG',
     >                   STOTS(35)/stots(36))
         call prc ('IMPURITY DENSITY JUST INSIDE SEPARATRIX(TOP HALF):')
         do 305 iz = 1,nizs
           write(7,'(''     IONIZATION STATE: '',
     >            I4,''  NI        '',G12.5)' ) iz,sitots(iz,2)
 305     continue
         call prr ('    ALL STATES                           ',
     >          sitots(maxizs+1,2))
c
      endif 
c
C-----------------------------------------------------------------------
C
c
c     Content Summary 
c
      call prb
c      call prc ('PLASMA IMPURITY CONTENT SUMMARY:')
      call prchtml('--- PLASMA IMPURITY CONTENT SUMMARY ---',
     >             'pr_content','0','B')
      call prb
c
      if (absfac.eq.0.0) then 
         mfact = 1.0
      else
         mfact = absfac
      endif 
c
      do iz = 0,nizs
        call pri('CHARGE STATE  :  ',iz)    

        call prr2('   - IN PLASMA GRID - DIVERTOR REGION (I/O):',
     >      impurity_content(iz,1,1)*mfact,
     >      impurity_content(iz,2,1)*mfact)
        call prr ('   - IN PLASMA GRID - MAIN SOL             :',
     >      impurity_content(iz,3,1)*mfact)
        call prr ('   - IN PLASMA GRID - CORE                 :',
     >      impurity_content(iz,4,1)*mfact)
        call prr ('   - IN PLASMA GRID - TOTAL                :',
     >      (impurity_content(iz,1,1) +    
     >       impurity_content(iz,2,1) +
     >       impurity_content(iz,3,1) +
     >       impurity_content(iz,4,1)) * mfact )
c
        if (iz.eq.0) then       

          call prr ('   - IN VOID        - MAIN REGION          :',
     >                                                sngl(ddvoid(1)))
          call prr ('   - IN VOID        - PRIVATE PLASMA       :',
     >                                                sngl(ddvoid(2)))
          call prr ('   - IN VOID        - DIVERTOR (NOT PP)    :',
     >                                                sngl(ddvoid(3)))
        endif
      end do         
c
      call prb
c
c
c     If auxiliary FLUID CODE data has been read in then print
c     the integrated impurity content 
c
      if (cre2d.eq.1.or.cre2d.eq.2) then       
c
         call prc ('FLUID CODE RESULT - IMPURITY CONTENT SUMMARY:')
c
         do iz = 0,cre2dizs
           call pri('CHARGE STATE  :  ',iz)    

           call prr2('   - IN PLASMA GRID - DIVERTOR REGION (I/O):',
     >         impurity_content(iz,1,2),
     >         impurity_content(iz,2,2))
           call prr ('   - IN PLASMA GRID - MAIN  SOL            :',
     >         impurity_content(iz,3,2))
           call prr ('   - IN PLASMA GRID - CORE                 :',
     >         impurity_content(iz,4,2))
           call prr ('   - IN PLASMA GRID - TOTAL                :',
     >      (impurity_content(iz,1,2) +    
     >       impurity_content(iz,2,2) +
     >       impurity_content(iz,3,2) +
     >       impurity_content(iz,4,2)))
c
         end do         
c
         call prb
c
      endif 
c
C-----------------------------------------------------------------------
C
c
c      call prc ('CALCUATION OF ESTIMATED DIVERTOR RETENTION:')
      call prb
      call prchtml('CALCUATION OF ESTIMATED DIVERTOR RETENTION:',
     >             'pr_divrent','0','B')
      call prb
      call prc ('        Impurity content of Divertor Region')
      call prc ('   DR = -----------------------------------')
      call prc ('        Impurity content of Main Region    ')
      call prb
      call prc (' - The Line Z=Zxp is used to divide these two'//
     >                ' regions.')
      call prb
      call prc (' Impurity Content of Divertor Region: ')
      call prr ('      Neutral Content (Void+Grid)     : ',
     >                          stots(43)+sngl(ddvoid(2)+ddvoid(3)))
      call prr ('      Ion Content     (Grid)          : ',
     >                                      stots(21))
c 
      div_content =   stots(21)+stots(43)+sngl(ddvoid(2)+ddvoid(3))
c
      call prr ('    Total Divertor Impurity Content   : ',div_content)
      call prc (' Impurity Content of Main Region:')
      call prr ('      Neutral Content (Void+Grid)     : ',
     >                        sngl(ddvoid(1))+stots(42))
      call prr ('      Ion Content in Main SOL         : ',
     >                   stots(6)-(stots(21)-stots(37)))
c
      carea = stots(18)**2 / (4.0*PI)
c
      call prr ('      - Estimated Area of Entire Core : ',carea)
      call prr ('      - Total Density Just Inside Core: ',
     >                          sitots(maxizs+1,1))
      call prr ('      Estimated Ion Content of Core   : ',
     >                          carea * sitots(maxizs+1,1))
c
      main_content = carea * sitots(maxizs+1,1)
     >    +sngl(ddvoid(1))+stots(42)+stots(6)-(stots(21)-stots(37))
c
      call prr ('   Total Main Impurity Content        : ',main_content)
c 
      if (main_content.gt.0.0) then 
         div_reten = div_content/main_content
      else
         div_reten = 0.0
      endif
c
      call prb
      call prr ('ESTIMATED VALUED FOR DIVERTOR RETENTION (DR) = ',
     >               div_reten)
      call prb
c
c     Additional Retention Data - Neutral content - various areas
c
      call prc('ADDITIONAL DATA ON PLASMA CONTENT:')
c    
      neutc_all = 0.0
      neutc_prim = 0.0
      neutc_temp = 0.0 
      neutc_area = 0.0
      neutc_div_density = 0.0 
c
      ir = irtrap+1
      do ik = 1,nks(ir)
c
         neutc_all = neutc_all 
     >                     + ddlims(ik,ir,0) * kareas(ik,ir)
         neutc_prim = neutc_prim
     >                     + ddlims(ik,ir,-1) * kareas(ik,ir)
         neutc_temp = neutc_temp
     >         + ddts(ik,ir,0) * ddlims(ik,ir,0) * kareas(ik,ir) 
c
         neutc_area = neutc_area + kareas(ik,ir) 
c
      end do
c      
      if (neutc_area.gt.0.0) then
       call prb
       call prr(' - TOTAL NEUTRAL CONTENT AT PFZ EDGE   ',neutc_all)  
       call prr(' - PRIMARY NEUTRAL CONTENT AT PFZ EDGE ',neutc_prim)  
       call prr(' - AREA OF PFZ EDGE                    ',neutc_area)
       neutc_divedge_density = neutc_all/neutc_area

       if (neutc_all.gt.0.0) then 
          neutc_ave_temp = neutc_temp/neutc_all
       else
          neutc_ave_temp = 0.0
       endif

       call prr(' - NEUTRAL DENSITY AT PFZ EDGE         ',
     >                                      neutc_divedge_density)
       call prr(' - AVER. NEUTRAL TEMP AT PFZ EDGE  (eV)',
     >                                      neutc_ave_temp)
c
      endif
c
      neutc_all = 0.0
      neutc_prim = 0.0
      neutc_area = 0.0
      neutc_div_density = 0.0 
c
      do ir = irtrap+1, nrs
         do ik = 1,nks(ir)
c
            neutc_all = neutc_all 
     >                     + ddlims(ik,ir,0) * kareas(ik,ir)
            neutc_prim = neutc_prim
     >                     + ddlims(ik,ir,-1) * kareas(ik,ir)
            neutc_area = neutc_area + kareas(ik,ir) 
c
         end do 
c
      end do
c   
      if (neutc_area.gt.0.0) then 
       call prb      
       call prr(' - TOTAL NEUTRAL CONTENT OF PFZ        ',neutc_all)  
       call prr(' - PRIMARY NEUTRAL CONTENT OF PFZ      ',neutc_prim)  
       call prr(' - AREA OF PFZ                         ',neutc_area)
       neutc_div_density = neutc_all/neutc_area
       call prr(' - NEUTRAL DENSITY FOR ALL OF PFZ      ',
     >                                      neutc_div_density)
      endif
c
      ionc_all   = 0.0
      ionc_area = 0.0
      ionc_coreedge_density = 0.0 
c
      ir = irsep-1
      do ik = 1,nks(ir)-1
         do iz = 1,nizs
c
            ionc_all = ionc_all 
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
         end do
c
         ionc_area = ionc_area + kareas(ik,ir) 
c
      end do
c      
      if (ionc_area.gt.0.0) then 
       call prb
       call prr(' - TOTAL ION CONTENT AT CORE EDGE      ',ionc_all)  
       call prr(' - AREA OF CORE EDGE                   ',ionc_area)
c
       ionc_coreedge_density = ionc_all/ionc_area
       call prr(' - ION DENSITY AT CORE EDGE            ',
     >                                      ionc_coreedge_density)
      endif
   
c
      ionc_all   = 0.0
      ionc_area = 0.0
      ionc_core_density = 0.0 
c
      do ir = ircore+1, irsep-1
         do ik = 1,nks(ir)-1
            do iz = 1,nizs
c
               ionc_all = ionc_all 
     >                     + ddlims(ik,ir,iz) * kareas(ik,ir)
c
            end do
c
            ionc_area = ionc_area + kareas(ik,ir) 
c
         end do
c
      end do
c      
      if (ionc_area.gt.0.0) then
       call prb
       call prr(' - TOTAL ION CONTENT OF CORE RINGS     ',ionc_all)  
       call prr(' - TOTAL AREA OF CORE RINGS            ',ionc_area)
c
       ionc_core_density = ionc_all/ionc_area
c
       call prr(' - ION DENSITY IN CORE                 ',
     >                                      ionc_core_density)
c
       endif
c
       call prb
c
      if (neutc_divedge_density.ne.0.0) then 
         call prr(' - EDGE DENSITY RATIO  (CORE/DIVERTOR) ',
     >              ionc_coreedge_density/neutc_divedge_density)
      else 
         call prr(' - EDGE DENSITY RATIO  (CORE/DIVERTOR) ',
     >              0.0)
      endif
c
      if  (ionc_coreedge_density.ne.0.0) then
         call prr(' - EDGE DENSITY RATIO  (DIVERTOR/CORE) ',
     >              neutc_divedge_density/ionc_coreedge_density)
      else
         call prr(' - EDGE DENSITY RATIO  (DIVERTOR/CORE) ',
     >              0.0)
      endif
c
      if (neutc_div_density.ne.0.0) then
         call prr(' - TOTAL DENSITY RATIO  (CORE/DIVERTOR)',
     >              ionc_core_density/neutc_div_density)
      else 
         call prr(' - TOTAL DENSITY RATIO  (CORE/DIVERTOR)',
     >              0.0)
      endif
c
      if (ionc_core_density.ne.0) then 
         call prr(' - TOTAL DENSITY RATIO  (DIVERTOR/CORE)',
     >              neutc_div_density/ionc_core_density)
      else 
         call prr(' - TOTAL DENSITY RATIO  (DIVERTOR/CORE)',
     >              0.0)
      endif 

      call prb 
      call prr('OVERALL AVERAGE NEUTRAL IMP TEMP (eV)',
     >                  sdtzs(0))
      call prb
c
      tot_incore = 0.0
      tot_inpfz  = 0.0
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
            tot_incore = tot_incore + ncore(ik,ir) 
            tot_inpfz  = tot_inpfz + ntrap(ik,ir)
         end do 
      end do
c
      call prr('TOTAL NUMBER OF IONS ENTERING PFZ  :',tot_inpfz)
      call prr('TOTAL NUMBER OF IONS ENTERING CORE :',tot_incore)
c
      call prb
c
      CALL PRR ('RADIATION CONSTANT IN MAIN LZ    (W.M**3)',STOTS(14))
      call prb 
      CALL PRR ('ION "ABSOLUTE" FACTOR      (IONS/m-tor/s)',ABSFAC_ion)
      CALL PRR ('NEUTRAL "ABSOLUTE" FACTOR  (NEUT/m-tor/s)',ABSFAC_neut)
      call prb
      CALL PRR ('POL. PERIPHERAL CIRCUM OF MAIN.  (M)     ',STOTS(18))
      call prr ('ESTIMATED TOTAL AREA OF MAIN     (M**2)  ',carea)
      WRITE (7,'('' TOTAL NO. OF IN/OUT C.F ADJS IN MAIN'',2I11)')
     >                    NINT(DCROSS(3)),NINT(DCROSS(1))
      CALL PRR2 ('AV.  SIZE OF IN/OUT C.F ADJS IN MAIN',
     >        DCROSS(4)/MAX(DCROSS(3),LO),DCROSS(2)/MAX(DCROSS(1),LO))
c
C-----------------------------------------------------------------------
C
c
c     Summary of radiation information
c
      call prb
c
c      call prc ('Summary Radiation statistics: ')
c
      call prchtml ('Summary Radiation statistics: ',
     >              'pr_rad','0','B')
      call prb
c
      call prc ('  - Ion values are normalized to net ion creation'//
     >          ' rate')
      call prc ('  - Neutral values are normalized to the net neutral'//
     >          ' creation rate')
      call prc ('  - Totals use the neutral creation rate for'//
     >          ' normalization')
c
      if (cdatopt.eq.1.or.cdatopt.eq.2.or.cdatopt.eq.3) then
         call prc ('  - Hydrogen values are absolute')
      endif
c
c     SOL + TRAP region
c
      call prb 
      call prc ('SOL + TRAP REGIONS:')
      CALL PRR ('  Power from ions        IN SOL       (W)  ',sTOTS(7))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                              sptots(7,maxizs+1))
      call prr ('      - Main SOL Region "above" Xpt   (W)  ',
     >                              sptots(9,maxizs+1))
      CALL PRR ('  Power from ions        IN TRAP      (W)  ',sTOTS(38))
      CALL PRR ('  Power from neutrals    IN SOL       (W)  ',STOTS(25))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                              sptots(7,0))
      call prr ('      - Main SOL Region "above" Xpt   (W)  ',
     >                              sptots(9,0))
      CALL PRR ('  Power from neutrals    IN TRAP      (W)  ',STOTS(40))
      CALL PRR ('  Line rad from ions     IN SOL       (W)  ',STOTS(8))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                              sptots(8,maxizs+1))
      call prr ('      - Main SOL Region "above" Xpt   (W)  ',
     >                              sptots(10,maxizs+1))
      CALL PRR ('  Line rad from ions     IN TRAP      (W)  ',STOTS(39))
      CALL PRR ('  Line rad from neutrals IN SOL       (W)  ',STOTS(26))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                              sptots(8,0))
      call prr ('      - Main SOL Region "above" Xpt   (W)  ',
     >                              sptots(10,0))
      CALL PRR ('  Line rad from neutrals IN TRAP      (W)  ',STOTS(41))
      call prr ('  Total Power            IN SOL+TRAP  (W)  ',stots(29))
      call prr ('  Total Line Rad         IN SOL+TRAP  (W)  ',stots(30))
c
c     Hydrogen data only if ADAS used.
c
      if (cdatopt.eq.1.or.cdatopt.eq.2.or.cdatopt.eq.3) then
      call prb
      call prr ('  Hydrogen Power         IN SOL+TRAP  (W)  ',stots(31))
      call prr ('  Hydrogen Line Rad      IN SOL+TRAP  (W)  ',stots(32))
      endif
c
c     MAIN plasma
c
      call prb
      call prc ('MAIN PLASMA REGION (CORE):')
      CALL PRR ('  Power from ions        IN MAIN      (W)  ',STOTS(4))
      CALL PRR ('  Power from neutrals    IN MAIN      (W)  ',STOTS(23))
      CALL PRR ('  Line rad from ions     IN MAIN      (W)  ',STOTS(5))
      CALL PRR ('  Line rad from neutrals IN MAIN      (W)  ',STOTS(24))
      call prr ('  Total Power            IN MAIN      (W)  ',stots(27))
      call prr ('  Total Line Rad         IN MAIN      (W)  ',stots(28))
c
c     Hydrogen data only if ADAS used.
c
      if (cdatopt.eq.1.or.cdatopt.eq.2.or.cdatopt.eq.3) then
      call prb
      call prr ('  Hydrogen Power         IN MAIN      (W)  ',stots(33))
      call prr ('  Hydrogen Line Rad      IN MAIN      (W)  ',stots(34))
      endif
c
c     Total for both regions
c
      call prb
      call prc ('ALL REGIONS:')
      CALL PRR ('  Total Power from ions               (W)  ',STOTS(4)
     >     +stots(7)+stots(38))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(7,maxizs+1)+sptots(5,maxizs+1))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(9,maxizs+1)+sptots(1,maxizs+1))
      CALL PRR ('  Total Power from neutrals           (W)  ',STOTS(23)
     >     +stots(25)+stots(40))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(7,0)+sptots(5,0))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(9,0)+sptots(1,0))
      CALL PRR ('  Total Line rad from ions            (W)  ',STOTS(5)
     >     +stots(8)+stots(39))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(8,maxizs+1)+sptots(6,maxizs+1))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(10,maxizs+1)+sptots(2,maxizs+1))
      CALL PRR ('  Total Line rad from neutrals        (W)  ',STOTS(24)
     >     +stots(26)+stots(41))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(8,0)+sptots(6,0))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(10,0)+sptots(2,0))
      call prr ('  Total Power                   PRAD  (W)  ',stots(27)
     >     +stots(29))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(7,maxizs+2)+sptots(5,maxizs+2))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(9,maxizs+2)+sptots(1,maxizs+2))
      call prr ('  Total Line Rad                      (W)  ',stots(28)
     >     +stots(30))
      call prr ('      - Divertor Region "below" Xpt   (W)  ',
     >                    sptots(8,maxizs+2)+sptots(6,maxizs+2))
      call prr ('      - Main Region     "above" Xpt   (W)  ',
     >                    sptots(10,maxizs+2)+sptots(2,maxizs+2))
c
c     Hydrogen data - only if ADAS used.
c
      if (cdatopt.eq.1.or.cdatopt.eq.2.or.cdatopt.eq.3) then
      call prr ('  Hydrogen Power                PRAD  (W)  ',stots(31)
     >     +stots(33))
      call prr ('  Hydrogen Line Rad                   (W)  ',stots(32)
     >     +stots(34))
      endif
c
      call prb
c
c     EDGE2D impurity radiation summary  
c
      if ((cre2d.eq.1.or.cre2d.eq.2).and.cre2dizs.gt.-1) then 
c
      call prb
      call prc ('*FLUID CODE* IMPURITY - Radiation Summary: ')
      call prc ('  - Units are absolute  (*FLUID CODE*/DIVIMP).')
c
c     SOL + TRAP region
c
      CALL PRR2('  Power from ions        IN SOL       (W)  ',
     >                     e2dtots(5),stots(7)*absfac)   
      CALL PRR2('  Power from ions        IN TRAP      (W)  ',
     >                     e2dtots(8),stots(38)*absfac)   
c
      CALL PRR2('  Power from neutrals    IN SOL       (W)  ',
     >                     e2dptots(3,0),stots(25)*absfac)
      CALL PRR2('  Power from neutrals    IN TRAP      (W)  ',
     >                     e2dptots(5,0),stots(40)*absfac)
c
      CALL PRR2('  Line rad from ions     IN SOL       (W)  ',
     >                     e2dtots(6),stots(8)*absfac)   
      CALL PRR2('  Line rad from ions     IN TRAP      (W)  ',
     >                     e2dtots(9),stots(39)*absfac)   
c
      CALL PRR2('  Line rad from neutrals IN SOL       (W)  ',
     >                     e2dptots(4,0),stots(26)*absfac)
      CALL PRR2('  Line rad from neutrals IN TRAP      (W)  ',
     >                     e2dptots(6,0),stots(41)*absfac)
c
      call prr2('  Total *FC*   Power     IN SOL+TRAP  (W)  ',e2dtots(5)     
     >   +  e2dptots(3,0)+e2dtots(8)+e2dptots(5,0),
     >                stots(29)*absfac)
      call prr2('  Total *FC*   Line Rad  IN SOL+TRAP  (W)  ',e2dtots(6)
     >   +  e2dptots(4,0)+e2dtots(9)+e2dptots(6,0),
     >                stots(30)*absfac)
c
c     MAIN plasma
c
      CALL PRR2('  Power from ions        IN MAIN      (W)  ',
     >                       e2dtots(2),stots(4)*absfac)
      CALL PRR2('  Power from neutrals    IN MAIN      (W)  ',
     >                       e2dptots(1,0),stots(23)*absfac)
      CALL PRR2('  Line rad from ions     IN MAIN      (W)  ',
     >                       e2dtots(3),stots(5)*absfac)
      CALL PRR2('  Line rad from neutrals IN MAIN      (W)  ',
     >                       e2dptots(2,0),stots(24)*absfac)
      call prr2('  Total *FC*   Power     IN MAIN      (W)  ',e2dtots(2)
     >               +  e2dptots(1,0),stots(27)*absfac)
      call prr2('  Total *FC*   Line Rad  IN MAIN      (W)  ',e2dtots(3)
     >               +  e2dptots(2,0),stots(28)*absfac)
c
c     Total for both regions
c
      CALL PRR2('  Total Power from ions               (W)  ',e2dTOTS(2)
     >         +e2dtots(5)+e2dtots(8),
     >          absfac*(stots(4)+stots(7)+stots(38)))
      CALL PRR2('  Total Power from neutrals           (W)  ',
     >     e2dpTOTS(1,0)+e2dptots(3,0)+e2dptots(5,0),
     >     absfac*(stots(23)+stots(25)+stots(40)))
      CALL PRR2('  Total Line rad from ions            (W)  ',e2dTOTS(3)
     >          +e2dtots(6)+e2dtots(9),
     >            absfac*(stots(5)+stots(8)+stots(39)))
      CALL PRR2('  Total Line rad from neutrals        (W)  ',
     >      e2dpTOTS(2,0)+e2dptots(4,0)+e2dptots(6,0),
     >          absfac*(stots(24)+stots(26)+stots(41)))
      call prr2('  Total *FC*   Power            PRAD  (W)  ',
     >          e2dtots(2)+e2dtots(5)+e2dptots(1,0)+e2dptots(3,0)+
     >          e2dtots(8)+e2dptots(5,0),
     >          absfac*(stots(27)+stots(29))) 
      call prr2('  Total *FC*   Line Rad               (W)  ',
     >          e2dtots(3)+e2dtots(6)+e2dptots(2,0)+e2dptots(4,0)+
     >          e2dtots(9)+e2dptots(6,0), 
     >          absfac*(stots(28)+stots(30))) 
      call prb
c
      endif
c
c
c---------------------------------------------------------------
C
c     Print out power summary
c
      call pr_power_summary
c
c---------------------------------------------------------------
C

      call prb  
      call prchtml('SUMMARY OF DIFFUSION AND FORCES',
     >             'pr_force','0','B')
      call prb 

      FORCE = CRMI * 1.67E-27 / (QTIM*QTIM)
      DO 889 IZ = 1, MAXIZS
        IF (DOUTS(IZ,1).GT.0.0) THEN
          DOUTS(IZ,2) = 2.0 * QTIM * DOUTS(IZ,2) / DOUTS(IZ,1)
          DOUTS(IZ,3) = QTIM * DOUTS(IZ,3) / DOUTS(IZ,1)
          DOUTS(IZ,4) = FORCE * DOUTS(IZ,4) / DOUTS(IZ,1)
          DOUTS(IZ,5) = FORCE * DOUTS(IZ,5) / DOUTS(IZ,1)
          DOUTS(IZ,6) = FORCE * DOUTS(IZ,6) / DOUTS(IZ,1)
          DOUTS(IZ,7) = FORCE * DOUTS(IZ,7) / DOUTS(IZ,1)
          DOUTS(IZ,8) = DOUTS(IZ,8) / DOUTS(IZ,1) / QTIM
          DOUTS(IZ,9) = DOUTS(IZ,9) / DOUTS(IZ,1) / QTIM
        ENDIF
c
        IF (COREOUTS(IZ,1).GT.0.0) THEN
          COREOUTS(IZ,2) = 2.0 * QTIM * COREOUTS(IZ,2) / COREOUTS(IZ,1)
          COREOUTS(IZ,3) = QTIM * COREOUTS(IZ,3) / COREOUTS(IZ,1)
          COREOUTS(IZ,4) = FORCE * COREOUTS(IZ,4) / COREOUTS(IZ,1)
          COREOUTS(IZ,5) = FORCE * COREOUTS(IZ,5) / COREOUTS(IZ,1)
          COREOUTS(IZ,6) = FORCE * COREOUTS(IZ,6) / COREOUTS(IZ,1)
          COREOUTS(IZ,7) = FORCE * COREOUTS(IZ,7) / COREOUTS(IZ,1)
          COREOUTS(IZ,8) = COREOUTS(IZ,8) / COREOUTS(IZ,1) / QTIM
          COREOUTS(IZ,9) = COREOUTS(IZ,9) / COREOUTS(IZ,1) / QTIM
        ENDIF
c
        DO 888 II = 1, 5
          IF (RIONS(IZ).GT.0.0) DPARAS(IZ,II) = DPARAS(IZ,II)/RIONS(IZ)
  888   CONTINUE
        IF (DPARAS(IZ,3).GT.0.0) DPARAS(IZ,6)=DPARAS(IZ,4)/DPARAS(IZ,3)
  889 CONTINUE


      CALL PRC ('AV. NO. OF DELTAS = 0 DIFF STEPS TAKEN BY EACH STATE')
      WRITE(7,9010)(IZ,NINT(SNGL(DPARAS(IZ,1))),IZ=1,NIZS)
      CALL PRC ('AV. NO. OF DELTAS = 0 DIFF STEPS TAKEN IN SOL+TRAP  ')
      WRITE(7,9010)(IZ,NINT(SNGL(DPARAS(IZ,5))),IZ=1,NIZS)
      CALL PRC ('AV. NO. OF DELTAS > 0 DIFF STEPS TAKEN BY EACH STATE')
      WRITE(7,9010)(IZ,NINT(SNGL(DPARAS(IZ,2))),IZ=1,NIZS)
      CALL PRC ('AV. NO. OF DELTAS > 0 DIFF STEPS TAKEN IN SOL+TRAP  ')
      WRITE(7,9010)(IZ,NINT(SNGL(DPARAS(IZ,3))),IZ=1,NIZS)
c
      CALL PRC ('AV. SIZE OF DELTAS > 0 STEPS IN SOL+TRAP REGION (M) ')
      WRITE(7,9011) (IZ,DPARAS(IZ,6),IZ=1,NIZS)
c
c     SOL + TRAP
c
      CALL PRC ('AVERAGE TAU PARALLEL IN SOL+TRAP                    ')
      WRITE(7,9011) (IZ,DOUTS(IZ,2),IZ=1,NIZS)
      CALL PRC ('AVERAGE TAU STOPPING IN SOL+TRAP                    ')
      WRITE(7,9011) (IZ,DOUTS(IZ,3),IZ=1,NIZS)
      CALL PRC ('AVERAGE FRICTION FORCE IN SOL+TRAP   FF             ')
      WRITE(7,9011) (IZ,DOUTS(IZ,4),IZ=1,NIZS)
      CALL PRC ('AVERAGE ELECTRIC FORCE IN SOL+TRAP   FE             ')
      WRITE(7,9011) (IZ,DOUTS(IZ,5),IZ=1,NIZS)
      CALL PRC ('AVERAGE GRADIENT FORCE IN SOL+TRAP   FEG            ')
      WRITE(7,9011) (IZ,DOUTS(IZ,6),IZ=1,NIZS)
      CALL PRC ('AVERAGE GRADIENT FORCE IN SOL+TRAP   FIG            ')
      WRITE(7,9011) (IZ,DOUTS(IZ,7),IZ=1,NIZS)
      CALL PRC ('AVERAGE ION VELOCITY IN SOL+TRAP     VIMP           ')
      WRITE(7,9011) (IZ,DOUTS(IZ,8),IZ=1,NIZS)
      CALL PRC ('AVERAGE BACKGROUND VEL IN SOL+TRAP   VB             ')
      WRITE(7,9011) (IZ,DOUTS(IZ,9),IZ=1,NIZS)
c
c
c     CORE 
c
      CALL PRC ('AVERAGE TAU PARALLEL IN CORE                    ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,2),IZ=1,NIZS)
      CALL PRC ('AVERAGE TAU STOPPING IN CORE                    ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,3),IZ=1,NIZS)
      CALL PRC ('AVERAGE FRICTION FORCE IN CORE   FF             ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,4),IZ=1,NIZS)
      CALL PRC ('AVERAGE ELECTRIC FORCE IN CORE   FE             ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,5),IZ=1,NIZS)
      CALL PRC ('AVERAGE GRADIENT FORCE IN CORE   FEG            ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,6),IZ=1,NIZS)
      CALL PRC ('AVERAGE GRADIENT FORCE IN CORE   FIG            ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,7),IZ=1,NIZS)
      CALL PRC ('AVERAGE ION VELOCITY IN CORE     VIMP           ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,8),IZ=1,NIZS)
      CALL PRC ('AVERAGE BACKGROUND VEL IN CORE   VB             ')
      WRITE(7,9011) (IZ,COREOUTS(IZ,9),IZ=1,NIZS)
c
c     Write out temporary case specific data 
c
      call write_tmp_data
c
c     Formats for the above write statements  
c
 9010 format (2x,6(i3,i8))     
 9011 format (2x,6(i3,1p,d8.1))     
c
      RETURN
      END
c
c
c
      subroutine write_tmp_data
      implicit none
      include 'params' 
      include 'comtor'
      include 'cgeom' 
      include 'dynam1'
c
      integer iktmp,iktmpc
      integer ik,ir
c
c     Write out the line density along the S=SMAX/2 line 
c
c     
c     jdemod - change ir=20 to ir=irsep so it will work for more grids
c            - could also change to (irsep+irwall-1)/2 if one wants the middle of the SOL
c
      iktmp = 0
      iktmpc = 0
c
      ir = irsep
c
      do ik = 1,nks(ir)
c
         if (ksb(ik-1,ir).le.ksmaxs(ir)/2.0.and.
     >       ksb(ik,ir).ge.ksmaxs(ir)/2.0) then 
           iktmp = ik 
           exit
         endif
      end do
c
c
c     jdemod - why this is behind a ippchange flag I do not know - since 
c              it is just identifying an out of bounds condition. The 
c              problem is probably the fixed use of ir=20 - which doesn't
c              make any sense - it can only be guaranteed to work for one grid
c            - Also -    
c
c slmod begin
c      IF (ippchange) THEN
c        IF (ik.EQ.nks(ir)+1) THEN
c
c     jdemod - test value of iktmp
c     
        IF (iktmp.EQ.0) THEN
          WRITE(6,*) 'WARNING: IKTMP NOT DEFINED',iktmp,ik
          RETURN
        ENDIF
c      ENDIF
c sldod end
c
c     Assign core: iktmpc value 
c
      iktmpc = ikins(iktmp,irsep)
c
      do ir = 1,irwall 
c
         if (ir.lt.irsep) then 
            ik = iktmpc
         else
            ik = iktmp
         endif
c
         if ((distin(ik,ir)+distout(ik,ir)).gt.0.0) then 
            write(6,'(a,2i6,6(1x,g20.8))')  'LINE DENSITY:',
     >       ik,ir,ddlims(ik,ir,1),kareas(ik,ir),
     >       distin(ik,ir),distout(ik,ir),
     >       ddlims(ik,ir,1)
     >       *kareas(ik,ir)/(distin(ik,ir)+distout(ik,ir))  
         else
            write(6,'(a,2i6,6(1x,g20.8))')  'LINE DENSITY:',
     >       ik,ir,ddlims(ik,ir,1),kareas(ik,ir),
     >       distin(ik,ir),distout(ik,ir)
         endif
c
      end do
      return
      end 
c
c
c
      subroutine pr_line_profile 
      implicit none
c
      include    'params'
      include    'line_profile'
c
c     Print out data related to the line_profile_opt calculation
c     if this was turned on
c
c
      real lam,lam0
      real*8 maxraw,maxinst
      integer in,ii
c
      real*8 instf,isler_inst_funct,dlam
      external isler_inst_funct
c
c     Initialization
c
      maxraw  = -HI
      maxinst = -HI
c
      call prb
c
      if (line_profile_opt.ne.0) then 
c
         call prc('LINE PROFILE CALCULATION OPTION WAS ACTIVATED:')
         call prr(' - WAVELENGTH OF LINE TRACKED         = ',
     >                                        lp_wave)
         call prr(' - SIZE OF WAVELENGTH BINS USED       = ',
     >                                        lp_bin_width)   
         if (lp_instrument_width.eq.-1.0) then 
            call prc(' - ISLER INSTRUMENT WIDTH FORMULA APPLIED')
            call prc('   IF = 0.95*EXP(-(DLAM/(1.201*0.1))**2 +')
            call prc('        0.05*EXP(-(DLAM/(1.201*0.35))**2')
         else
            call prr(' - EFFECTIVE INSTRUMENT WIDTH APPLIED = ',
     >                                        lp_instrument_width)
         endif
         call pri(' - NUMBER OF BINS IN USE          +/-   ',
     >                                        max_lp_bins)       
c
c        Calculate the modified line profile for the 
c        instrument width. 
c
         call dzero(mod_line_profile,2*max_lp_bins+1) 
c
         do in = -max_lp_bins+1, max_lp_bins-1
c
            do ii =  -max_lp_bins+1, max_lp_bins-1
c
               lam0 = in * lp_bin_width
               lam  = ii * lp_bin_width
c 
               dlam = lam-lam0
c
               if (lp_instrument_width.eq.-1.0) then 
                  instf = isler_inst_funct(dlam) 
               else
                  instf =  exp(-((lam-lam0)/lp_instrument_width)**2)
               endif

               mod_line_profile(in)=mod_line_profile(in) +
     >                          line_profile(ii) * instf

            end do 
c
            maxraw = max(maxraw,line_profile(in))
            maxinst = max(maxinst,mod_line_profile(in))
c
        end do
c
c       Write these profiles to the UNIT 6 file:
c
        write(6,'(a)')  
        write(6,'(a)')  'LINE PROFILE CALCULATION:'
        write(6,'(a)')  
        write(6,'(a,f12.5)')  'WAVELENGTH:',lp_wave
        write(6,'(a)')  
        do in =  -max_lp_bins, max_lp_bins
           write(6,'(i5,1x,f10.5,4(1x,g20.12))')
     >             in,in*lp_bin_width,line_profile(in),
     >             mod_line_profile(in),
     >             line_profile(in)/maxraw,
     >             mod_line_profile(in)/maxinst
        end do
c 
      endif

      return 
      end

      real*8 function const_inst_funct(dlam,lp_instrument_width)
      implicit none
      real*8 dlam
      real*8 lp_instrument_width 

      const_inst_funct = exp(-(dlam/lp_instrument_width)**2)

      return
      end

      real*8 function isler_inst_funct(dlam)
      implicit none
      real*8 dlam
      
      isler_inst_funct = 0.95 * exp(-(dlam/(1.201*0.1))**2) +
     >                   0.05 * exp(-(dlam/(1.201*0.35))**2)

      return
      end
