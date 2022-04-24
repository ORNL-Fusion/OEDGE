      SUBROUTINE TAUIN1 (QTIM,NIZS,ICUT,                                        
     >     FSRATE,IGEOM,NTBS,NTIBS,NNBS)                         
      use mod_params
      use mod_comt2
      use variable_wall
      use mod_comtor
      use mod_cadas
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_global_options
      use mod_slcom
      use yreflection
      use mod_vtig
      use mod_lambda
      use debug_options
      IMPLICIT  none
      REAL      QTIM,FSRATE                                                     
      INTEGER   NIZS,ICUT(2),IGEOM,NTBS,NTIBS,NNBS,IQXBRK
C     
C***********************************************************************
C     
C     SETS UP VALUES IN COMMON BLOCKS COMTAU/COMT2                            
C     
C***********************************************************************
C     
      INTEGER   IPOS,IXOUT,IZ,IQX,LIMIZ,IX,IY,J                                 
      CHARACTER MESAGE*80
      REAL      TMPION      
      real      width
c     
      integer :: solver_opt
      integer :: pz
      pz = 1 

      call pr_trace('TAUIN1','Start')

!     jdemod - use common code for calculating lambda. LIM default is
! option lambda_opt = 2
      lambda = coulomb_lambda(cnbin,ctibin)
      
C     
      LIMIZ = MIN (CION, NIZS)                                                  
c     write(0,*) 'LIMIZ:', cion,nizs,limiz
      write(0,"(A,I2,A,I2,A,I2)") ' CION = ',cion, ' NIZS = ',nizs, 
     >     ' --> LIMZ = ',limiz
C     
C-----------------------------------------------------------------------
C     SET UP QEDGES, QTANS AND QDISTS                           
c     
c     jdemod - setup geometry BEFORE calculating plasma since some of
c     the plasma calculations may depend on the location of the
c     limiter edges. 
C-----------------------------------------------------------------------
C     
      call pr_trace('TAUIN1','Before EDGE')

      WRITE (6,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      WRITE (0,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      CALL EDGE (QXS,QEDGES,QTANS,QDISTS,NQXSO,CAW,CL,ICUT,CCUT,XSCALO,         
     >     WEDGAN,XL1,YL1,XL2,YL2,TC,SC,TO,SO,GC,RP,CIOPTH,CORECT,
     >     XST1,YST1,XST2,YST2,XST3,YST3,RLEDGE7,CA,RLC)        

C     
C-----------------------------------------------------------------------
C     SET UP CTEMBS AND CRNBS, QTEMBS AND QRNBS                 
C-----------------------------------------------------------------------
C     
      call pr_trace('TAUIN1','Before PLASMA')

      WRITE (6,'('' TAU: CALLING PLASMA OPTIONS'',I3,I3)') 
     >     CIOPTG,CIOPTK                   
      WRITE (0,'('' TAU: CALLING PLASMA OPTIONS'',I3,I3)') 
     >     CIOPTG,CIOPTK                   

      CALL PLASMA (NTBS,NTIBS,NNBS,CIOPTG,CIOPTK,QTIM)                     
C     
C-----------------------------------------------------------------------
C     SET UP VARIABLE CAW
C-----------------------------------------------------------------------
c     
      call pr_trace('TAUIN1','Before Setup_wall')
      call setup_wall (qys,nqys,cl,caw)
C     
C-----------------------------------------------------------------------
C     SET UP CEYS AND CVHYS                                     
C-----------------------------------------------------------------------
C     
      call pr_trace('TAUIN1','Before SOL')
      
      WRITE (6,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      WRITE (0,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      CALL SOL (QYS,CEYS,CVHYS,NQYS,CTBIN,CTIBIN,CRMB,CL,CIZB,                 
     >     CEYOUT,CVHOUT,CYSTAG,CRMI,CSOLEF,CIOPTF)                        
c
c     jdemod - after base plasma complete - calculate temperature gradients
c     
c     jdemod
c     moved this before the plasma overlay code - the plasma overlay code will recalculate the temperature
c     gradients as required to accomodate the fact that the midpoints may not be at +/-nys/2 among other
c     complicating factors      
c     
      call pr_trace('TAUIN1','Before calculate_tgrad')
      call calculate_tgrad(qtim)
c     
c     
C-----------------------------------------------------------------------
C     SET UP CYSCLS, CFEXZS, CFVHXS                                    
C     ALL VALUES DEFINED FOR OUTBOARD X POSITIONS ONLY                 
c     
c     These need to be calculated inboard for colprobe3d on
c     
C-----------------------------------------------------------------------
C     
      DO 200  IQX = 1-NQXSO, 0                                                  
         WIDTH      = CTWOL - QEDGES(IQX,1) - QEDGES(IQX,2)                     
         CYSCLS(IQX)= REAL(NQYS) / WIDTH                                        
 200  CONTINUE                                                                  

! jdemod - efield calculations use cyscls
! cyscls appears to be intended to compress an efield calculated for all cells
! from 0 to 2CL onto an actual plasma length of 2CL-EDGE1-EDGE2
! these factors are only used with the base plasma options - the overlays
! calculate the efield and other quantities from EDGE1 to EDGE2 and set any cells
! outside this range to the boundary values at the respective limiter surfaces

      call pr_trace('TAUIN1','Before calculate_efield')

      ! the following routine sets the base values of the efield and velplasma arrays to the
      ! contents of cvhys and ceys - so that base plasma options will be able to work on sections
      ! without more realistic plasma specifications. The velplasma and efield overlays will
      ! include inverse temperature scaling factors to cancel those implicitly included in
      ! cfexzs (efield) and cfvhxs (velocity). 
      call calculate_efield(qtim,limiz)

c     
c     jdemod - at this point the LIM plasma has been fully calculated
c     using the base options which are simple and relatively efficient
c     - both SOLEDGE and SOL22 are then implemented as overlays on top of
c     this existing plasma description. velplasma and efield can be initialized
c     using this so that it contains valid values for inboard and other regions 
c     where SOL22 is not appropriate
c     - these routines include re-calculating and overwriting the temperature gradients
c     as required 
c     
      call pr_trace('TAUIN1','Before plasma_overlay')
      call plasma_overlay(qtim)
      
c     write(0,*) 'crmi,qtim:',crmi,qtim
c     
c     Depending on the plasma overlay option specified - rewrite the ion temperature to
c     be constant at the target value.       
c     
      if (vtig_opt.eq.2) then 
         do pz = 1,maxpzone
            do ix = 1,nxs
               iqx = iqxs(ix)
               do iy = -nys/2,-1
                  ctembsi(ix,iy,pz) = qtembsi(iqx,1)
                  ctembsi(ix,iy+nys+1,pz) =  qtembsi(iqx,1)
               end do   
               do iy = -nys,-nys/2-1
                  ctembsi(ix,iy,pz) = qtembsi(iqx,2)
                  ctembsi(ix,iy+nys+1,pz) =  qtembsi(iqx,2)
               end do   
            end do
         end do
      endif

      
! jdemod - write out background plasma
      if (cprint.eq.9) then 
         do pz = 1, maxpzone
            write(6,*) 'PLASMA BACKGROUND AFTER OVERLAY:',pz
            do ix = 1,nxs
               do iy = -nys,nys
                  write(6,'(3i8,10(1x,g12.5))')ix,iy,pz,
     >                 xouts(ix),youts(iy),crnbs(ix,iy,pz),
     >                 ctembs(ix,iy,pz),
     >                 ctembsi(ix,iy,pz), velplasma(ix,iy,pz),
     >                 efield(ix,iy,pz),             
     >                 ctegs(ix,iy,pz),ctigs(ix,iy,pz)
                     
               end do
            end do
         end do
      endif

      
C     
c     jdemod - pull out code to calculate varying dperp - must be run
c     after plasma density is finalized       
c     
      call pr_trace('TAUIN1','Before setup_dperp')
      call setup_dperp(qtim,igeom)

c     
c     jdemod: Scale the outboard flow velocity outside the limiters
c     (only used in 3D) 
c     (QX(IQX) scaling is done in line when the new location
c     is calculated in lim3.f)
c     
      vpflow_3d = vpflow_3d * qtim

C     
C-----------------------------------------------------------------------
C     SET UP CMIZS                                              
C-----------------------------------------------------------------------
C---- NOTE:                                                                     
C---- THE NORMAL ROUTE IS TO DISABLE IONISATION BEYOND LEVEL NIZS.              
C---- A SPECIAL CASE IS ALLOWED WHEN OPTA=0 SUCH THAT IONISATION FROM           
C---- LEVEL NIZS TO NIZS+1 IS ALLOWED, BUT WHEN IT OCCURS THE ION IS            
C---- NO LONGER FOLLOWED   (LEAP TO LABEL 790 IN LIM ROUTINE)                   
C     
      IF (CIOPTA.EQ.0 .OR. CIOPTA.EQ.3 .OR. CIOPTA.EQ.4) THEN                   
         CMIZS = MIN (CION, NIZS+1)                                             
      ELSE                                                                      
         CMIZS = MIN (CION, NIZS)                                               
      ENDIF                                                                     
C     
C-----------------------------------------------------------------------
C     SET IONISATION / E-I RECOMBINATION TIME INTERVALS    CFIZS,CFRCS           
C-----------------------------------------------------------------------
C     
      call pr_trace('TAUIN1','Before IZTAU')
      WRITE (0,'('' TAU: CALLING IZTAU  OPTION'',I3)') CIOPTA                   
      CALL IZTAU (CRMI,NXS,NYS,CION,CIZB,CIOPTA)                                
C     
C-----------------------------------------------------------------------
C     SET COMBINED C-X AND E-I RECOMBINATION TIMES         CFCXS                 
C-----------------------------------------------------------------------
C     
      call pr_trace('TAUIN1','Before CXREC')
      WRITE (0,'('' TAU: CALLING CXREC  OPTION'',I3)') CIOPTI                   
      CALL CXREC (NIZS,CION,CIOPTI,CIZB,CL,CRMB,CVCX,                           
     >     CNHC,CNHO,CLAMHX,CLAMHY)                                  
C     
C-----------------------------------------------------------------------
C     SET PROBABILITY OF EITHER AN IONISATION OR A RECOMBINATION                
C     SET PROPORTION OF THESE WHICH WILL BE RECOMBINATIONS                      
C     PREVENT ANY IONISATION BEYOND MAXIMUM LIMIT SPECIFIED IF REQUIRED         
C-----------------------------------------------------------------------
C     
c     slmod tmp
c     
c     IF (CIOPTE.EQ.10.AND.CDATOPT.EQ.0) THEN
c     DO IZ = 1, CMIZS-1
c     DO IX = 1, NXS
c     DO IY = -NYS,NYS
c     CFIZS(IX,IY,IZ) = 1.0 / (TMPION*CNBIN)
c     ENDDO
c     ENDDO
c     ENDDO
c     ENDIF
c     slmod end

      call pr_trace('TAUIN1','Before calculating state change')


      IF (CMIZS.GT.1) THEN                                                      
         do pz = 1,maxpzone
            DO 370 IZ = 1, CMIZS-1                                                  
               DO 360 IX = 1, NXS                                                     
                  IQX = IQXS(IX)                                                        
                  DO 350 IY = -NYS, NYS                                                 
                     IF (CFCXS(IX,IY,IZ,pz).LE.0.0) THEN                                    
                        CPCHS(IX,IY,IZ,pz) = QTIM * QS(IQX)
     >                                      / CFIZS(IX,IY,IZ,pz)                
                        CPRCS(IX,IY,IZ,pz) = 0.0                                             
                     ELSE                                                                
                        CPCHS(IX,IY,IZ,pz) =
     >                       (CFIZS(IX,IY,IZ,pz) + CFCXS(IX,IY,IZ,pz)) *           
     >                       QTIM * QS(IQX)
     >                        /(CFIZS(IX,IY,IZ,pz) * CFCXS(IX,IY,IZ,pz))             

! jdemod - cfizs and cfcxs contain characteristic TIMES (see print outs below)
! fast ionization is a SMALLER time meaning less chance for recombination. 
                        
                        CPRCS(IX,IY,IZ,pz) = CFIZS(IX,IY,IZ,pz) /                               
     >                       (CFCXS(IX,IY,IZ,pz) + CFIZS(IX,IY,IZ,pz))             
                     ENDIF                                                               
                     CPCHS(IX,IY,IZ,pz) = MIN (1.0, CPCHS(IX,IY,IZ,pz))                        
 350              CONTINUE                                                              
 360           CONTINUE                                                               
 370        CONTINUE                                                                
         end do
      ENDIF                                                                     
C     
      IF (CMIZS .LE. LIMIZ) THEN                                                
         do pz = 1,maxpzone
            DO 390 IX = 1, NXS                                                      
               IQX = IQXS(IX)                                                        
               DO 380 IY = -NYS, NYS                                                 
                  IF (CFCXS(IX,IY,IZ,pz).LE.0.0) THEN                                    
                     CPCHS(IX,IY,CMIZS,pz) = 0.0                                          
                  ELSE                                                                
                     CPCHS(IX,IY,CMIZS,pz) =QTIM * QS(IQX)
     >                     / CFCXS(IX,IY,IZ,pz)             
                  ENDIF                                                               
                  CPRCS(IX,IY,CMIZS,pz) = 1.0                                            
                  CPCHS(IX,IY,IZ,pz) = MIN (1.0, CPCHS(IX,IY,IZ,pz))                        
 380           CONTINUE                                                              
 390        CONTINUE                                                                
         end do
      ENDIF                                                                     
C     
C---- SET IONISATION PROBABILITIES FOR NEUT ...                                 
C---- (SAVES REPEATED CALCULATION EVERY ITERATION)                              
C---- THEY ARE ALL MULTIPLIED BY THE "IONISATION RATE FACTOR" IRF               
C---- WHICH TYPICALLY MIGHT BE 0.2 TO GIVE DEEPER IONISATION.                   
C     
      do pz = 1,maxpzone
         DO IY = -NYS, NYS                                                     
            DO IX = 1, NXS                                                      
               CPCHS(IX,IY,0,pz) = MIN (1.0, CIRF*FSRATE
     >               / CFIZS(IX,IY,0,pz))            
            end do
         end do
      end do
C     
c     slmod begin - N2 break
! N2 currently not 3D poloidal zone compatible
      pz = 1
      IF (N2OPT.EQ.1) THEN
         WRITE(63,*) ' '
         WRITE(63,*) 'Ionisation probabilities:'
         WRITE(63,*) ' '
         WRITE(63,'(A3,3A12)')
     +        'IX','NNCPCHS','N2CPCHS','CFIZS0'

         DO IX = 1, NXS
            N2CPCHS(IX) = MIN (1.0, CIRF * FSRATE / N2RATE(IX))            
            NNCPCHS(IX) = MIN (1.0, CIRF * FSRATE / NNRATE(IX))            

            WRITE(63,'(I3,3E12.5)') 
     +           IX,NNCPCHS(IX),N2CPCHS(IX),CPCHS(IX,1,0,pz)
         ENDDO
         WRITE(63,*) ' '
      ENDIF

c     slmod end


c
c     end of tau - write out entire plasma solution if degugging is on

      call wrt_plasma
      

      call pr_trace('TAUIN1','End')
      WRITE(0,*) 'TAU finished'   
      
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE TAUIN2 (QTIM,NIZS)                                             
      use mod_params
      use mod_comt2
      use mod_comtor
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_lambda
      IMPLICIT  none
      REAL      QTIM                                                            
      INTEGER   NIZS                                                            
C                                                                               
C***********************************************************************        
C                                                                               
C       SETS UP VALUES IN COMMON BLOCKS COMTAU/COMT2                            
C       THESE VALUES MAY DEPEND ON CTEMSC, AVERAGE IONISATION TEMP              
C       RETURNED FROM A CALL TO NEUT.                                           
C                                                                               
C***********************************************************************        
C                                                                               
C                                                                               
      INTEGER   IPOS,IZ,IQX,LIMIZ,JX,IX,IY                                      
!      real      tmp1
      
      ! test for issues with precision - especially when TAU values are multiplied
      REAL*8    TEMP,FTAU,FTAUP,FTAUS,FTAUT,RIZSQR,STAU,TAU                     
      real*8    lambda
      REAL*8    ROOTMI,ROOTTT                                            

      integer :: pz

c      
      if (lambda_vary_opt.eq.0) then 
         lambda = coulomb_lambda(cnbin,ctibin)
      else
         lambda = 1.0
      endif
         
c
C         LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)
c       ELSE
c         LAMBDA = 15.0
c       ENDIF
c slmod end

      LIMIZ  = MIN (CION, NIZS)                                                 
      ROOTMI = SQRT (CRMI)                                                      
      ROOTTT = SQRT (CTEMSC)                                                    
C                                                                               
C---- DETERMINE CROSSOVER POINT WHERE NON-ZERO CHARACTERISTICS WILL             
C---- CHANGE TO NORMAL PLASMA CHARACTERISTICS.  WE WANT THE SPECIALS            
C---- OUTSIDE OF X = CPLSMA ONLY.  FOR UNIFORM SPECIAL CHARACTERISTICS,         
C---- SET CPLSMA > CA  (SAY TO 99.0).                                           
C                                                                               
      JX = IPOS (CPLSMA*0.999, XS, NXS-1)                                       
C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP CFPS, CFSS, CFTS                                   
C  NOTE 212: BACKGROUND ION TEMPERATURE CAN BE OVERRIDDEN WITH CONSTANT         
C  VALUE "CTBI" IF SPECIFIED AS > 0.0                                           
C  NOTE 215: EXTRA HEATING OPTION, SET UP CONSTANTS C215A AND C215B.            
C  SET CTOLDS ARRAY TO ION TEMPERATURES USED IN THESE CALCULATIONS.             
C-----------------------------------------------------------------------        
C                                                                               
      FTAU  = CZENH * SQRT(CRMB) * CIZB * CIZB * LAMBDA * QTIM * sf_tau                 
      FTAUP = FTAU * 6.8E-14                                                    
      FTAUS = FTAU * 6.8E-14 * (1.0 + CRMB/CRMI)                                
      FTAUT = FTAU * 1.4E-13                                                    
      C215A = FTAU * 1.4E-13 * SQRT(CRMI)                                       
      C215B = (CRMI * CTBI + CRMB * CTEMSC) ** 1.5                              
c
      write(6,'(a,3i8,10(1x,g12.5))') 'TAU CONSTANTS:',
     >       cioptb,cioptc,cioptd,ftau,
     >       ftaup,ftaus, ftaut,c215a,c215b
      
C
      do pz = 1,maxpzone
         
      DO 540  IZ = 1, LIMIZ                                                     
         RIZSQR = REAL (IZ) * REAL (IZ)                                         
         DO 520 IY = -NYS, NYS                                                  
          DO 520  IX = 1, NXS                                                   
            IF (CTBI.GT.0.0) THEN                                               
             TEMP  = CRNBS(IX,IY,pz) / (CRMI * CTBI**1.5)                          
            ELSE                                                                
             TEMP  = CRNBS(IX,IY,pz) / (CRMI * CTEMBSI(IX,IY,pz)**1.5)                
            ENDIF                                                               
            IQX = IQXS(IX)                                                      
!
            ! jdemod - move lambda into STAU if lambda varies with plasma conditions
            if (lambda_vary_opt.eq.1) then
               lambda =coulomb_lambda(crnbs(ix,iy,pz),ctembsi(ix,iy,pz))
            else
               lambda = 1.0
            endif

            STAU = TEMP * RIZSQR * QS(IQX) * LAMBDA                                     
            CTOLDS(IX,IZ) = CTEMSC                                              
C                                                                               
C-----------------------------------------------------------------------        
C           TAU PARALLEL           NOTES 3,50,103                               
C-----------------------------------------------------------------------        
C                                                                               
C---------- STANDARD CASE:-                                                     
C---------- CFPS = 2.DELTAT.TI/TAUPARALLEL.                                     
C---------- FOR NON-ZERO OPTIONS, SPECIAL CASE APPLIES FOR X OUTSIDE OF         
C---------- CPLSMA ONLY.  EACH TIME AN ION ENTERS THIS REGION ITS               
C---------- TEMPERATURE IS COMPARED AGAINST A PREVIOUS VALUE TO SEE             
C---------- WHETHER ANY VALUES NEED RECALCULATING.                              
C                                                                               
            IF (CTBI.GT.0.0) THEN                                               
              CFPS(IX,IY,IZ,pz) = STAU * CTBI * FTAUP * 2.0                        
            ELSE                                                                
              CFPS(IX,IY,IZ,pz) = STAU * CTEMBSI(IX,IY,pz) * FTAUP * 2.0             
            ENDIF                                                               
C                                                                               
            IF     (CIOPTB.EQ.1.AND.IX.LE.JX) THEN                              
               CFPS(IX,IY,IZ,pz) = 0.0                                             
            ELSEIF (CIOPTB.EQ.2.AND.IX.LE.JX) THEN                              
               CFPS(IX,IY,IZ,pz) = 2.0 * CRNBS(IX,IY,pz) * 6.8E-14 *                  
     >                         REAL (CIZB * CIZEFF) * RIZSQR * LAMBDA /         
     >                         (ROOTMI * ROOTTT) * QTIM * QS(IQX)               
            ENDIF                                                               
C                                                                               
C---------- SET ADDITIONAL ARRAY FOR USE IN INNERMOST LOOP OF LIM2              
C---------- EQUIVALENT TO SOME CONSTANT TIMES CFPS VALUES                       
C---------- THIS SAVES 20% OF CPU TIME BY ELIMINATING A SQUARE ROOT             
C---------- CCCFPS = SQRT (4.88E8/(CFPS.MI)) . DELTAT. SIN(THETAB)              
C---------- FOR NOTE 284, SET CCCFPS = AS ABOVE . SQRT(2TI/CFPS)                
C                                                                               
            IF (CFPS(IX,IY,IZ,pz).EQ.0.0) THEN                                     
               CCCFPS(IX,IY,IZ,pz) = 0.0                                           
            ELSEIF (CIOPTB.EQ.3 .AND. IX.LE.JX) THEN 
               CCCFPS(IX,IY,IZ,pz) = SQRT (9.76E8 * CTEMSC / CRMI) *               
     >           QTIM * QS(IQX) / CFPS(IX,IY,IZ,pz)                                
            ELSEIF (CIOPTB.EQ.4 .AND. IX.LE.JX) THEN
               CFPS(IX,IY,IZ,pz) = 2.0E0 * CFPS(IX,IY,IZ,pz)
               CCCFPS(IX,IY,IZ,pz) = SQRT(9.76E8 * CTEMSC / CRMI ) *
     >            QTIM * QS(IQX) / CFPS(IX,IY,IZ,pz)
c slmod
            ELSEIF (CIOPTB.EQ.13) THEN

c CRMB    - plasma ion mass
c CRMI    - impurity ion mass
c CIZB    - plasma ion charge
c CTEMBSI - local background temperature
c CRNBS   - local background density

c               TPARA = DTEMI*12.0*SQRT(CTEMBSI(IX,IY)/2.0)/6.8E-14
c     +                 /LAMBDA/CNBIN/(1+2.0/12.0)
c               VPARAT = SQRT(2.0*1.6E-19*DTEMI/1.67E-27/12.0)*
c     +                  SQRT(QTIM/TPARA)

c               TPARA = CRMI*SQRT(CTEMBSI(IX,IY)/CRMB)/6.8E-14
c     +                 /LAMBDA/CNBIN/(1+CRMB/CRMI)/CIZB**2/
c     +                 REAL(IZ)*REAL(IZ)
c               VPARAT = SQRT(2.0*1.6E-19/1.67E-27/CRMB)*
c     +                  SQRT(QTIM/TPARA)

c               CCCFPS(IX,IY,IZ) = 
c               tmp1  = 
c     +                  SQRT(2.0*1.6E-19/1.67E-27/CRMI)*
c     +                  SQRT(QTIM/
c     +                   (CRMI*SQRT(CTEMBSI(IX,IY)/CRMB)/6.8E-14
c     +                    /LAMBDA/CRNBS(IX,IY)/(1+CRMB/CRMI)/
c     +                    (REAL(CIZB)*REAL(CIZB))/
c     +                    (REAL(IZ)*REAL(IZ))) ) * qtim

               CCCFPS(ix,iy,iz,pz) = 1.56e4 * SQRT(PI/4.0 * 1.0/CRMI
     >             * (cfps(ix,iy,iz,pz) *(1.0+CRMB/CRMI))  /2.0) * qtim

               
c               WRITE (78,'(a,3i8,20(1x,g12.5))')
c     +                  'CCCFPS:',ix,iy,iz,tmp1,CCCFPS(IX,IY,IZ),
c     +                   CRNBS(IX,IY),CTEMBSI(IX,IY),
c     +                   QTIM,CRMB,CRMI,CIZB,REAL(IZ),LAMBDA 
c slmod end
            ELSE                                                                
              CCCFPS(IX,IY,IZ,pz)=SQRT(4.88E8/(CFPS(IX,IY,IZ,pz)*CRMI))*         
     >           QTIM * QS(IQX)                                                 
            ENDIF                                                               
c     Apply the scaling factor to the parallel diffusive transport 
c     sf_vdiff default value is 1.0
c     This quantity is used to scale either spatial or velocity diffusive
c     step sizes depending on which is in use      
            cccfps(ix,iy,iz,pz) = cccfps(ix,iy,iz,pz) * sf_vdiff
            
C     
C-----------------------------------------------------------------------        
C           TAU STOPPING           NOTES 3,50,103                               
C-----------------------------------------------------------------------        
C                                                                               
C---------- STANDARD CASE:-                                                     
C---------- CFSS = 1 - EXP (-DELTAT/TAUSTOPPING)                                
C---------- NON ZERO OPTIONS APPLY OUTSIDE OF X = CPLSMA ONLY ...               
C---------- FOR OPTION 2, TAUSTOPPING = TAUPARALLEL                             
C----------                           = 2.DELTAT.TI/CFPS                        
C                                                                               
            TAU = STAU * FTAUS                                                  
            IF (TAU.GT.1.E-3) THEN                                              
              CFSS(IX,IY,IZ,pz) = 1.0 - EXP(-TAU)                                  
            ELSE                                                                
              CFSS(IX,IY,IZ,pz) = TAU                                              
            ENDIF                                                               

!            write(6,'(a,3i8,10(1x,g12.5))') 'CFSS:',ix,iy,iz,
!     >            cfss(ix,iy,iz),tau,stau,ftaus            
C
            
            IF     (CIOPTC.EQ.1.AND.IX.LE.JX) THEN                              
               CFSS(IX,IY,IZ,pz) = 0.0                                             
            ELSEIF (CIOPTC.EQ.2.AND.IX.LE.JX) THEN                              
               TAU = CFPS(IX,IY,IZ,pz) / (2.0*CTEMSC)                              
               IF (TAU.GT.1.E-3) THEN                                           
                 CFSS(IX,IY,IZ,pz) = 1.0 - EXP(-TAU)                               
               ELSE                                                             
                 CFSS(IX,IY,IZ,pz) = TAU                                           
               ENDIF                                                            
            ELSEIF (CIOPTC.EQ.3 .AND. IX.LE.JX) THEN
                 CFSS(IX,IY,IZ,pz) = 1.0E20 
            ENDIF                                                               
C                                                                               
C-----------------------------------------------------------------------        
C           TAU HEATING            NOTES 3,50,103,215                           
C-----------------------------------------------------------------------        
C                                                                               
C---------- STANDARD CASE:-                                                     
C---------- CFTS = 1 - EXP (-DELTAT/TAUHEATING)                                 
C---------- NON ZERO OPTIONS APPLY OUTSIDE OF X = CPLSMA ONLY ...               
C                                                                               
            TAU = STAU * FTAUT                                                  
            IF (TAU.GT.1.E-3) THEN                                              
              CFTS(IX,IY,IZ,pz) = 1.0 - EXP(-TAU)                                  
            ELSE                                                                
              CFTS(IX,IY,IZ,pz) = TAU                                              
            ENDIF                                                               
C                                                                               
            IF     (CIOPTD.EQ.1.AND.IX.LE.JX) THEN                              
               CFTS(IX,IY,IZ,pz) = 0.0                                             
            ELSEIF (CIOPTD.EQ.2.AND.IX.LE.JX) THEN                              
               CFTS(IX,IY,IZ,pz) = 1.0                                             
            ELSEIF (CIOPTD.EQ.3.AND.IX.LE.JX) THEN                              
               IF (CTBI.LE.0.0)                                                 
     >           C215B = (CRMI * CTEMBsI(IX,IY,pz)+ CRMB*CTEMSC) ** 1.5      
               ! jdemod - include lambda in C215A (from FTAU) for cases
               ! where spatially varying lambda is in use
               TAU = C215A * RIZSQR * QS(IQX) * CRNBS(IX,IY,pz) / C215B            
     >                  * LAMBDA
               IF (TAU.GT.1.E-3) THEN                                           
                 CFTS(IX,IY,IZ,pz) = 1.0 - EXP(-TAU)                               
               ELSE                                                             
                 CFTS(IX,IY,IZ,pz) = TAU                                           
               ENDIF                                                            
            ENDIF                                                               
520      CONTINUE                                                               
540   CONTINUE                                                                  

      end do ! end of pz loop
C

      
      !write(6,*) 'CFSS,CFTS:'
      !do ix = 1,nxs
      !   do iy = -nys,nys
      !      write(6,'(2i8,10(1x,g12.5))') 
     >!        ix,iy,cfss(ix,iy,iz),cfts(ix,iy,iz) 
      !   end do
      !end do

      call check_tau(cion,nizs)


      
      RETURN                                                                    
      END                                                                       

C
      subroutine check_tau(cion,nizs)

      use mod_params
      use mod_comt2
      !use mod_comtor
      !use mod_comtau
      use mod_comxyt
      !use mod_coords
      
      implicit none
      integer :: cion,nizs
      
      real :: tau_warn(3,4,maxizs+1,maxpzone)
      real :: tau_ave(3,4,maxizs+1,maxpzone)
      real :: tau_cnt
      integer :: ix,iy,iz

      integer :: pz
      !pz = 1
      tau_warn= 0.0
      tau_ave = 0.0
      tau_cnt = 0.0

      do pz = 1,maxpzone
      DO  IZ = 1,  MIN (CION, NIZS)
        DO IX = 1, NXS
          DO IY = -NYS,NYS

            tau_cnt = tau_cnt + 1.0
             
            ! jdemod
            ! add diagnostic checks on the values of cfss, cfts and cfps
            ! These are QTIM/TAU where TAU is TAU_Stopping, TAU_Heating and
            ! TAU_parallel ... these should all be << 1
            if (cfts(ix,iy,iz,pz).ge.1.0) then 
               write(6,'(a,4i8,20(1x,g12.5))')
     >              'CFTS > 1:',ix,iy,iz,pz,
     >              crnbs(ix,iy,pz),ctembsi(ix,iy,pz),cfts(ix,iy,iz,pz)
               tau_warn(1,1,iz,pz) = tau_warn(1,1,iz,pz) +1.0
               tau_ave(1,1,iz,pz) = tau_ave(1,1,iz,pz)+cfts(ix,iy,iz,pz)
            elseif (cfts(ix,iy,iz,pz).ge.0.1) then
               tau_warn(1,2,iz,pz) = tau_warn(1,2,iz,pz) +1.0
               tau_ave(1,2,iz,pz) = tau_ave(1,2,iz,pz)+cfts(ix,iy,iz,pz)
            elseif (cfts(ix,iy,iz,pz).ge.0.01) then
               tau_warn(1,3,iz,pz) = tau_warn(1,3,iz,pz) +1.0
               tau_ave(1,3,iz,pz) = tau_ave(1,3,iz,pz)+cfts(ix,iy,iz,pz)
            else
               tau_warn(1,4,iz,pz) = tau_warn(1,4,iz,pz) +1.0
               tau_ave(1,4,iz,pz) = tau_ave(1,4,iz,pz)+cfts(ix,iy,iz,pz)
            endif   
c
            if (cfss(ix,iy,iz,pz).ge.1.0) then 
               write(6,'(a,4i8,20(1x,g12.5))')
     >              'CFSS > 1:',ix,iy,iz,pz,
     >              crnbs(ix,iy,pz),ctembsi(ix,iy,pz),cfss(ix,iy,iz,pz)
               tau_warn(2,1,iz,pz) = tau_warn(2,1,iz,pz) +1.0
               tau_ave(2,1,iz,pz) = tau_ave(2,1,iz,pz)+cfss(ix,iy,iz,pz)
            elseif (cfss(ix,iy,iz,pz).ge.0.1) then
               tau_warn(2,2,iz,pz) = tau_warn(2,2,iz,pz) +1.0
               tau_ave(2,2,iz,pz) = tau_ave(2,2,iz,pz)+cfss(ix,iy,iz,pz)
            elseif (cfss(ix,iy,iz,pz).ge.0.01) then
               tau_warn(2,3,iz,pz) = tau_warn(2,3,iz,pz) +1.0
               tau_ave(2,3,iz,pz) = tau_ave(2,3,iz,pz)+cfss(ix,iy,iz,pz)
            else
               tau_warn(2,4,iz,pz) = tau_warn(2,4,iz,pz) +1.0
               tau_ave(2,4,iz,pz) = tau_ave(2,4,iz,pz)+cfss(ix,iy,iz,pz)
            endif   
c
            if (cfps(ix,iy,iz,pz).ge.1.0) then 
               write(6,'(a,4i8,20(1x,g12.5))')
     >              'CFPS > 1:',ix,iy,iz,pz,
     >              crnbs(ix,iy,pz),ctembsi(ix,iy,pz),cfps(ix,iy,iz,pz)
               tau_warn(3,1,iz,pz) = tau_warn(3,1,iz,pz) +1.0
               tau_ave(3,1,iz,pz) = tau_ave(3,1,iz,pz)+cfps(ix,iy,iz,pz)
            elseif (cfps(ix,iy,iz,pz).ge.0.1) then
               tau_warn(3,2,iz,pz) = tau_warn(3,2,iz,pz) +1.0
               tau_ave(3,2,iz,pz) = tau_ave(3,2,iz,pz)+cfps(ix,iy,iz,pz)
            elseif (cfps(ix,iy,iz,pz).ge.0.01) then
               tau_warn(3,3,iz,pz) = tau_warn(3,3,iz,pz) +1.0
               tau_ave(3,3,iz,pz) = tau_ave(3,3,iz,pz)+cfps(ix,iy,iz,pz)
            else
               tau_warn(3,4,iz,pz) = tau_warn(3,4,iz,pz) +1.0
               tau_ave(3,4,iz,pz) = tau_ave(3,4,iz,pz)+cfps(ix,iy,iz,pz)
            endif   
c            
              
         enddo
        enddo
      enddo
      end do 
      
      do pz = 1,maxpzone
      do iz = 1,  MIN (CION, NIZS)
         do iy = 1,3
            do ix = 1,4
              tau_warn(iy,ix,maxizs+1,pz) = tau_warn(iy,ix,maxizs+1,pz)
     >                                 + tau_warn(iy,ix,iz,pz)
              tau_ave(iy,ix,maxizs+1,pz) = tau_ave(iy,ix,maxizs+1,pz)
     >                                + tau_ave(iy,ix,iz,pz)
           end do
        end do
      end do
      end do

      do pz = 1,maxpzone
      do iz = 1, MAXIZS
         do iy = 1,3
            do ix = 1,4
              tau_warn(iy,ix,maxizs+1,pz) = tau_warn(iy,ix,maxizs+1,pz)
     >              + tau_warn(iy,ix,iz,pz)
              if (tau_warn(iy,ix,iz,pz).ne.0.0) then 
                 tau_ave(iy,ix,iz,pz)=tau_ave(iy,ix,iz,pz)
     >                              /tau_warn(iy,ix,iz,pz)
              else
                 tau_ave(iy,ix,iz,pz)=0.0
              endif
           end do
        end do
      end do
      end do 

      do pz = 1,maxpzone
!     issue tau warnings
      if (tau_warn(1,1,maxizs+1,pz).ne.0.0.or.
     >    tau_warn(2,1,maxizs+1,pz).ne.0.0.or.
     >    tau_warn(3,1,maxizs+1,pz).ne.0.0.or.
     >    tau_warn(1,2,maxizs+1,pz).ne.0.0.or.
     >    tau_warn(2,2,maxizs+1,pz).ne.0.0.or.
     >    tau_warn(3,2,maxizs+1,pz).ne.0.0) then
         write(0,*) 'WARNING: Time step may be too large in some'//
     >               ' cells for some charge states' ,pz
c         write(0,*) 'Tau Para warnings are only for'//
c     >              ' spatial diffusion options'
         write(0,*) 'Total ix,iy,iz checked = ', tau_cnt
         write(0,'(14x,6x,a,5x,5x,a,4x,4x,a,4x,5x,a)')
     >           'dt/Tau','>1','>0.1','>0.01','rest'
         write(0,'(a,8(1x,g12.5))')
     >        'Tau_t warn   :',tau_warn(1,1,maxizs+1,pz),
     >                         tau_warn(1,2,maxizs+1,pz),
     >                         tau_warn(1,3,maxizs+1,pz),
     >                         tau_warn(1,4,maxizs+1,pz)
         write(0,'(a,8(1x,g12.5))')
     >        'dt/Tau_t ave :',tau_ave(1,1,maxizs+1,pz),
     >                         tau_ave(1,2,maxizs+1,pz),
     >                         tau_ave(1,3,maxizs+1,pz),
     >                         tau_ave(1,4,maxizs+1,pz)

         write(0,'(a,8(1x,g12.5))')
     >        'Tau_s warn   :',tau_warn(2,1,maxizs+1,pz),
     >                         tau_warn(2,2,maxizs+1,pz),
     >                         tau_warn(2,3,maxizs+1,pz),
     >                         tau_warn(2,4,maxizs+1,pz)
         write(0,'(a,8(1x,g12.5))')
     >        'dt/Tau_s ave :',tau_ave(2,1,maxizs+1,pz),
     >                         tau_ave(2,2,maxizs+1,pz),
     >                         tau_ave(2,3,maxizs+1,pz),
     >                         tau_ave(2,4,maxizs+1,pz)

c         write(0,'(a,8(1x,g12.5))')
c     >        'Tau_p warn   :',tau_warn(3,1,maxizs+1),
c     >                         tau_warn(3,2,maxizs+1),
c     >                         tau_warn(3,3,maxizs+1),
c     >                         tau_warn(3,4,maxizs+1)
c         write(0,'(a,8(1x,g12.5))')
c     >        'dt/Tau_p ave :',tau_ave(3,1,maxizs+1),
c     >                         tau_ave(3,2,maxizs+1),
c     >                         tau_ave(3,3,maxizs+1),
c     >                         tau_ave(3,4,maxizs+1)


      endif

      write(6,*) 'TAU testing results >1, >0.1, >0.01 (not inclusive):' 
      write(6,*) 'Tau Para warnings are for spatial'//
     >              ' diffusion options only'
      write(6,*) 'Results for zone = ',pz
      write(6,*) 'Total ix,iy,iz checked = ', tau_cnt
         write(6,'(a,8(1x,g12.5))')
     >        'Tau_t warn   :',tau_warn(1,1,maxizs+1,pz),
     >                         tau_warn(1,2,maxizs+1,pz),
     >                         tau_warn(1,3,maxizs+1,pz),
     >                         tau_warn(1,4,maxizs+1,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_t ave :',tau_ave(1,1,maxizs+1,pz),
     >                         tau_ave(1,2,maxizs+1,pz),
     >                         tau_ave(1,3,maxizs+1,pz),
     >                         tau_ave(1,4,maxizs+1,pz)

         write(6,'(a,8(1x,g12.5))')
     >        'Tau_s warn   :',tau_warn(2,1,maxizs+1,pz),
     >                         tau_warn(2,2,maxizs+1,pz),
     >                         tau_warn(2,3,maxizs+1,pz),
     >                         tau_warn(2,4,maxizs+1,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_s ave :',tau_ave(2,1,maxizs+1,pz),
     >                         tau_ave(2,2,maxizs+1,pz),
     >                         tau_ave(2,3,maxizs+1,pz),
     >                         tau_ave(2,4,maxizs+1,pz)

         write(6,'(a,8(1x,g12.5))')
     >        'Tau_p warn   :',tau_warn(3,1,maxizs+1,pz),
     >                         tau_warn(3,2,maxizs+1,pz),
     >                         tau_warn(3,3,maxizs+1,pz),
     >                         tau_warn(3,4,maxizs+1,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_p ave :',tau_ave(3,1,maxizs+1,pz),
     >                         tau_ave(3,2,maxizs+1,pz),
     >                         tau_ave(3,3,maxizs+1,pz),
     >                         tau_ave(3,4,maxizs+1,pz)


      write(6,*) 'TAU testing results >1,>0.1,>0.01 (by charge state):' 


      do iz = 1, MIN (CION,NIZS)

         write(6,*) 'TAU testing results >1,>0.1,>0.01'//
     >             ' (by charge state) IZ=:',iz 

         write(6,'(a,8(1x,g12.5))')
     >        'Tau_t warn   :',tau_warn(1,1,iz,pz),
     >                         tau_warn(1,2,iz,pz),
     >                         tau_warn(1,3,iz,pz),
     >                         tau_warn(1,4,iz,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_t ave :',tau_ave(1,1,iz,pz),
     >                         tau_ave(1,2,iz,pz),
     >                         tau_ave(1,3,iz,pz),
     >                         tau_ave(1,4,iz,pz)

         write(6,'(a,8(1x,g12.5))')
     >        'Tau_s warn   :',tau_warn(2,1,iz,pz),
     >                         tau_warn(2,2,iz,pz),
     >                         tau_warn(2,3,iz,pz),
     >                         tau_warn(2,4,iz,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_s ave :',tau_ave(2,1,iz,pz),
     >                         tau_ave(2,2,iz,pz),
     >                         tau_ave(2,3,iz,pz),
     >                         tau_ave(2,4,iz,pz)

         write(6,'(a,8(1x,g12.5))')
     >        'Tau_p warn   :',tau_warn(3,1,iz,pz),
     >                         tau_warn(3,2,iz,pz),
     >                         tau_warn(3,3,iz,pz),
     >                         tau_warn(3,4,iz,pz)
         write(6,'(a,8(1x,g12.5))')
     >        'dt/Tau_p ave :',tau_ave(3,1,iz,pz),
     >                         tau_ave(3,2,iz,pz),
     >                         tau_ave(3,3,iz,pz),
     >                         tau_ave(3,4,iz,pz)

       end do

       end do  ! end of pz loop


      end
C     
      SUBROUTINE TAUPR1 (QTIM,NIZS)                                             
      use mod_params
      use mod_comt2
      use mod_comtor
      use mod_comtau 
      use mod_comxyt
      use mod_printr
      use mod_global_options
      IMPLICIT  none
      INTEGER   NIZS                                                            
      REAL      QTIM                                                            
C     
C***********************************************************************
C     
C     PRINTS REPRESENTATIVE VALUES FROM THE ARRAYS OF                         
C     FACTORS SET UP IN COMMONS COMTAU & COMT2                                  
C     CHRIS FARRELL    DECEMBER 1987                             
C     
C***********************************************************************
C     
c     INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c     INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c     INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c     INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c     INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c     INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c     
c     include   'global_options'
C     
      INTEGER   IQX,IZ,PMIZS,MAXIZ,J                                            

      integer :: pz
! use poloidal zone near P=0
      pz = pzones(1)
C     
C-----------------------------------------------------------------------
      MAXIZ  = MIN (NIZS, CION)                                                 
      PMIZS  = CMIZS                                                            
      IF ((CIOPTA.EQ.0 .OR. CIOPTA.EQ.3 .OR. CIOPTA.EQ.4) .AND.                 
     >     NIZS.LT.CION) PMIZS = CMIZS - 1                                      
C-----------------------------------------------------------------------
      CALL PRB                                                                  
      IF (CTBI.GT.0.0) THEN                                                     
         CALL PRC ('FUNCTIONS    LIMITER    PLASMA TEMPS(TBI) DENSITIES            
     >        TIMESTEP  ')                                                              
      ELSE                                                                      
         CALL PRC ('FUNCTIONS    LIMITER    PLASMA TEMPS      DENSITIES            
     >        TIMESTEP  ')                                                              
      ENDIF                                                                     

      CALL PRI ('  OF X      Y<0   Y>0   Y<0     Y>0      Y<0      Y>0          
     >     FACTOR    PZ=',pz)                                                              

      WRITE (7,9107) 'X = A   ',                                                
     >     CTEMBS(IXA  ,-1,pz),CTEMBS(IXA  ,1,pz),                      
     >     CRNBS (IXA  ,-1,pz),CRNBS (IXA  ,1,pz),QS(IQXA  )                  
      WRITE (7,9107) 'X = A/2 ',                                                
     >     CTEMBS(IXA2 ,-1,pz),CTEMBS(IXA2 ,1,pz),                      
     >     CRNBS (IXA2 ,-1,pz),CRNBS (IXA2 ,1,pz),QS(IQXA2 )                  
      WRITE (7,9107) 'X = A/4 ',                                                
     >     CTEMBS(IXA4 ,-1,pz),CTEMBS(IXA4 ,1,pz),                      
     >     CRNBS (IXA4 ,-1,pz),CRNBS (IXA4 ,1,pz),QS(IQXA4 )                  
      WRITE (7,9107) 'INBOARD ',                                                
     >     CTEMBS(IXIN ,-1,pz),CTEMBS(IXIN ,1,pz),                      
     >     CRNBS (IXIN ,-1,pz),CRNBS (IXIN ,1,pz),QS(IQXIN )                  
      WRITE (7,9108) 'OUTBOARD',-QEDGES(IQXOUT,1),                              
     >     QEDGES(IQXOUT,2),CTEMBS(IXOUT,-1,pz),CTEMBS(IXOUT,1,pz),                      
     >     CRNBS (IXOUT,-1,pz),CRNBS (IXOUT,1,pz),QS(IQXOUT)                  
      WRITE (7,9108) 'X =-AW/4',-QEDGES(IQXAW4,1),                              
     >     QEDGES(IQXAW4,2),CTEMBS(IXAW4,-1,pz),CTEMBS(IXAW4,1,pz),                      
     >     CRNBS (IXAW4,-1,pz),CRNBS (IXAW4,1,pz),QS(IQXAW4)                  
      WRITE (7,9108) 'X =-AW/2',-QEDGES(IQXAW2,1),                              
     >     QEDGES(IQXAW2,2),CTEMBS(IXAW2,-1,pz),CTEMBS(IXAW2,1,pz),                      
     >     CRNBS (IXAW2,-1,pz),CRNBS (IXAW2,1,pz),QS(IQXAW2)                  
      WRITE (7,9108) 'X =-AW  ',-QEDGES(IQXAW ,1),                              
     >     QEDGES(IQXAW ,2),CTEMBS(IXAW ,-1,pz),CTEMBS(IXAW ,1,pz),                      
     >     CRNBS (IXAW ,-1,pz),CRNBS (IXAW ,1,pz),QS(IQXAW )                  
      IF (IQXFAC.LT.0)                                                          
     >     WRITE (7,9108) 'X =-LAM ',-QEDGES(IQXFAC,1),                            
     >     QEDGES(IQXFAC,2),CTEMBS(IXFAC,-1,pz),CTEMBS(IXFAC,1,pz),                      
     >     CRNBS (IXFAC,-1,pz),CRNBS (IXFAC,1,pz),pz,QS(IQXFAC)                  
C     
      IF (CTBI.GT.0.0) THEN                                                     
         CALL PRB                                                                
         CALL PRR ('*** BACKGROUND ION TEMP OVERRIDDEN TO VALUE TBI =',          
     >        CTBI)                                                                 
      ELSE
C     
C     IF THE BACKGROUND ION TEMPERATURE IS NOT OVERRIDDEN BY   
C     A CONSTATNT, PRINT OUT REPRESENTATIVE VALUES. 
C     
         CALL PRC ('FUNCTIONS    LIMITER    PLASMA ION TEMPS                       
     >        TIMESTEP  ')                                                              
         CALL PRI ('  OF X      Y<0   Y>0   Y<0     Y>0                            
     >        FACTOR    PZ=',pz)                                                              
         WRITE (7,9109) 'X = A   ',                                                
     >        CTEMBSI(IXA  ,-1,pz),CTEMBSI(IXA  ,1,pz),QS(IQXA  )                  
         WRITE (7,9109) 'X = A/2 ',                                                
     >        CTEMBSI(IXA2 ,-1,pz),CTEMBSI(IXA2 ,1,pz), QS(IQXA2 )                  
         WRITE (7,9109) 'X = A/4 ',                                                
     >        CTEMBSI(IXA4 ,-1,pz),CTEMBSI(IXA4 ,1,pz),QS(IQXA4 )                  
         WRITE (7,9109) 'INBOARD ',                                                
     >        CTEMBSI(IXIN ,-1,pz),CTEMBSI(IXIN ,1,pz),QS(IQXIN )                  
         WRITE (7,9110) 'OUTBOARD',-QEDGES(IQXOUT,1),                              
     >        QEDGES(IQXOUT,2),CTEMBSI(IXOUT,-1,pz),CTEMBSI(IXOUT,1,pz),                   
     >        QS(IQXOUT)                  
         WRITE (7,9110) 'X =-AW/4',-QEDGES(IQXAW4,1),                              
     >        QEDGES(IQXAW4,2),CTEMBSI(IXAW4,-1,pz),CTEMBSI(IXAW4,1,pz),                   
     >        QS(IQXAW4)                  
         WRITE (7,9110) 'X =-AW/2',-QEDGES(IQXAW2,1),                              
     >        QEDGES(IQXAW2,2),CTEMBSI(IXAW2,-1,pz),CTEMBSI(IXAW2,1,pz),                    
     >        QS(IQXAW2)                  
         WRITE (7,9110) 'X =-AW  ',-QEDGES(IQXAW ,1),                              
     >        QEDGES(IQXAW ,2),CTEMBSI(IXAW ,-1,pz),CTEMBSI(IXAW ,1,pz),                   
     >        QS(IQXAW )                  
         IF (IQXFAC.LT.0)                                                          
     >        WRITE (7,9110) 'X =-LAM ',-QEDGES(IQXFAC,1),                            
     >        QEDGES(IQXFAC,2),CTEMBSI(IXFAC,-1,pz),CTEMBSI(IXFAC,1,pz),                    
     >        QS(IQXFAC)                  
C     
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF (CPRINT.EQ.1.or.cprint.eq.9) THEN
         DO 222 J = 1, 3
            CALL PRB                                                                  
            IF (J.EQ.1) THEN                                                          
               CALL PRC ('X DRIFT & INWARD PINCH (M),'//
     >              ' & INWARD STEP PROBS -L/2 < Y < 0')          
            ELSEIF (J.EQ.2) THEN                                                
               CALL PRC ('X DRIFT & INWARD PINCH (M),'//
     >               ' & INWARD STEP PROBS  0 < Y < L/2')          
            ELSE                                                                      
               CALL PRC ('X DRIFT & INWARD PINCH (M),'//
     >                   ' & INWARD STEP PROBS  L/2< AY <3L/2')          
            ENDIF                                                                     
            WRITE (7,9003) 'NEAR X = A     ', CXBFS(IQXA,J) /                       
     >           SQRT(QS(IQXA  )),CXAFS(IQXA,J)/QS(IQXA)
     >           ,CXCFS(IQXA  ,J)             
            WRITE (7,9003) 'NEAR X = A/2   ', CXBFS(IQXA2,J) /                       
     >           SQRT(QS(IQXA2 )),CXAFS(IQXA2,J)/QS(IQXA2)
     >           ,CXCFS(IQXA2 ,J)             
            WRITE (7,9003) 'NEAR X = A/4   ', CXBFS(IQXA4,J) /                       
     >           SQRT(QS(IQXA4 )),CXAFS(IQXA4,J)/QS(IQXA4)
     >            ,CXCFS(IQXA4 ,J)             
            WRITE (7,9003) 'JUST INBOARD   ', CXBFS(IQXIN,J) /                       
     >           SQRT(QS(IQXIN )),CXAFS(IQXIN,J)/QS(IQXIN)
     >            ,CXCFS(IQXIN ,J)             
            WRITE (7,9003) 'JUST OUTBOARD  ', CXBFS(IQXOUT,J) /                       
     >           SQRT(QS(IQXOUT)),CXAFS(IQXOUT,J)/QS(IQXOUT)
     >            ,CXCFS(IQXOUT,J)             
            WRITE (7,9003) 'NEAR X =-AW/4  ', CXBFS(IQXAW4,J) /                       
     >           SQRT(QS(IQXAW4)),CXAFS(IQXAW4,J)/QS(IQXAW4)
     >            ,CXCFS(IQXAW4,J)             
            WRITE (7,9003) 'NEAR X =-AW/2  ', CXBFS(IQXAW2,J) /                       
     >           SQRT(QS(IQXAW2)),CXAFS(IQXAW2,J)/QS(IQXAW2)
     >            ,CXCFS(IQXAW2,J)             
            WRITE (7,9003) 'NEAR X =-AW    ', CXBFS(IQXAW,J) /                       
     >           SQRT(QS(IQXAW )),CXAFS(IQXAW,J)/QS(IQXAW)
     >            ,CXCFS(IQXAW ,J)             
            IF (IQXFAC.LT.0)                                                          
     >           WRITE (7,9003) 'NEAR X =-LAMBDA', CXBFS(IQXFAC,J) /                       
     >           SQRT(QS(IQXFAC)),CXAFS(IQXFAC,J)/QS(IQXFAC)
     >            ,CXCFS(IQXFAC,J)             
 222     CONTINUE                                                                  
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF (MAXIZ.LT.1) GOTO 999                                                  
C-----------------------------------------------------------------------
      IF (CIOPTI.NE.0)                                                          
     >     CALL TAUPRA (7,CNHS,'NEUTRAL HYDROGEN DENSITY (M**-3)',-1)              
C-----------------------------------------------------------------------
      IF (CPRINT.EQ.1.or.cprint.eq.9) THEN
         do pz = 1,maxpzone
            CALL PRB                                                                  
            CALL PRI ('IONISATION AND RECOMBINATION TIMES   PZ=',pz)                           
            WRITE (7,9101)                                                            
            CALL PRB                                                                  
            CALL PRC ('  VALUES NEAR X = A/2  (Y>0)')                                        
            WRITE (7,9102) 0,CFIZS(IXA2,IY0,0,pz)                                        
            DO IZ = 1, PMIZS                                                       
               WRITE (7,9102) IZ,CFIZS(IXA2 ,IY0,IZ,pz)
     >               ,CFRCS(IXA2 ,IY0,IZ,pz),               
     >              CFCXS(IXA2 ,IY0,IZ,pz),CFCXS(IXA2 ,IYL8,IZ,pz),
     >              CFCXS(IXA2 ,IYL4,IZ,pz)          
            end do
c     call prb
c     CALL PRC ('  VALUES NEAR X = A/2  (Y<0)')                                        
c     WRITE (7,9102) 0,CFIZS(IXA2,IY0LT,0)                                        
c     DO IZ = 1, PMIZS                                                       
c     WRITE (7,9102) IZ,CFIZS(IXA2 ,IY0LT,IZ),CFRCS(IXA2 ,IY0LT,IZ),               
c     >   CFCXS(IXA2,IY0LT,IZ),CFCXS(IXA2,-IYL8,IZ),CFCXS(IXA2 ,-IYL4,IZ)          
c     end do
C-----------------------------------------------------------------------
            CALL PRB                                                                  
            CALL PRC ('  VALUES JUST OUTBOARD (Y>0)')                                       
            WRITE (7,9102) 0,CFIZS(IXOUT,IY0,0,pz)                                       
            DO IZ = 1, PMIZS                                                       
               WRITE (7,9102) IZ,CFIZS(IXOUT,IY0,IZ,pz)
     >               ,CFRCS(IXOUT,IY0,IZ,pz),               
     >              CFCXS(IXOUT,IY0,IZ,pz),CFCXS(IXOUT,IYL8,IZ,pz),
     >              CFCXS(IXOUT,IYL4,IZ,pz)          
            end do
            CALL PRB                                                                  
            CALL PRC ('  VALUES JUST OUTBOARD (Y<0)')                                       
            WRITE (7,9102) 0,CFIZS(IXOUT,IY0LT,0,pz)                                       
            DO IZ = 1, PMIZS                                                       
               WRITE (7,9102) IZ,CFIZS(IXOUT,IY0LT,IZ,pz),
     >              CFRCS(IXOUT,IY0LT,IZ,pz),               
     >              CFCXS(IXOUT,IY0LT,IZ,pz),CFCXS(IXOUT,-IYL8,IZ,pz),
     >              CFCXS(IXOUT,-IYL4,IZ,pz)          
            end do
C-----------------------------------------------------------------------
            IF (IQXFAC.LT.0) THEN                                                     
               CALL PRB                                                                
               CALL PRR ('  VALUES NEAR X = -LAMBDA (Y>0) = '
     >                   ,QXS(IQXFAC))                   
               WRITE (7,9102) 0,CFIZS(IXFAC,IY0,0,pz)                                     
               DO IZ = 1, PMIZS                                                     
                  WRITE (7,9102) IZ,CFIZS(IXFAC,IY0,IZ,pz),
     >                 CFRCS(IXFAC,IY0,IZ,pz),            
     >                 CFCXS(IXFAC,IY0,IZ,pz),CFCXS(IXFAC,IYL8,IZ,pz),
     >                 CFCXS(IXFAC,IYL4,IZ,pz)        
               end do
               CALL PRB                                                                
               CALL PRR ('  VALUES NEAR X = -LAMBDA (Y<0) = '
     >                  ,QXS(IQXFAC))                   
               WRITE (7,9102) 0,CFIZS(IXFAC,IY0LT,0,pz)                                     
               DO IZ = 1, PMIZS                                                     
                  WRITE (7,9102) IZ,CFIZS(IXFAC,IY0LT,IZ,pz),
     >                 CFRCS(IXFAC,IY0LT,IZ,pz),            
     >                 CFCXS(IXFAC,IY0LT,IZ,pz)
     >                 ,CFCXS(IXFAC,-IYL8,IZ,pz),
     >                 CFCXS(IXFAC,-IYL4,IZ,pz)        
               end do
            ENDIF                                                                     
C-----------------------------------------------------------------------
C     
C     NOTE ON THIS SECTION WHEN VARIABLE DELTAT IS IN USE:-  IF A PROB. IS          
C     SAY 0.22 FOR THE STANDARD DT, THEN WHEN A TIMESTEP MULTIPLE OF SAY            
C     1000 IS USED THIS SHOOTS UP TO "220", WHICH IS ADJUSTED TO 1 IN               
C     TAUIN1.  THEN WHEN WE ATTEMPT TO RECREATE THE GENUINE PROB. FOR PRINT         
C     PURPOSES BELOW, IT BECOMES 0.001 WHICH IS OBVIOUSLY WRONG.  HOWEVER,          
C     WHAT ELSE CAN ONE DO?  AT LEAST THE CALCULATION USES CORRECT VALUES.          
C     
            CALL PRB                                                                  
            CALL PRC ('CHANGE OF STATE PROBABILITIES (USING LOCAL TIMESTEPS FO        
     >           R IONS)')                                                                 
            WRITE (7,9103)                                                            
            CALL PRB                                                                  
            CALL PRC ('  VALUES NEAR X = A/2 (Y>0)')                                        
            DO IZ = 0, PMIZS                                                       
               WRITE (7,9104) IZ,                                                      
     >              CPCHS(IXA2 ,IY0 ,IZ,pz)
     >               ,100.0*CPRCS(IXA2 ,IY0 ,IZ,pz),            
     >              CPCHS(IXA2 ,IYL8,IZ,pz)
     >               ,100.0*CPRCS(IXA2 ,IYL8,IZ,pz),            
     >              CPCHS(IXA2 ,IYL4,IZ,pz)
     >               ,100.0*CPRCS(IXA2 ,IYL4,IZ,pz)             
            end do
c     CALL PRB                                                                  
c     CALL PRC ('  VALUES NEAR X = A/2 (Y<0)')                                        
c     DO IZ = 0, PMIZS                                                       
c     WRITE (7,9104) IZ,                                                      
c     >    CPCHS(IXA2 ,IY0LT ,IZ)         ,100.0*CPRCS(IXA2 ,IY0LT,IZ),            
c     >    CPCHS(IXA2 ,-IYL8,IZ)          ,100.0*CPRCS(IXA2 ,-IYL8,IZ),            
c     >    CPCHS(IXA2 ,-IYL4,IZ)          ,100.0*CPRCS(IXA2 ,-IYL4,IZ)             
c     end do
C-----------------------------------------------------------------------
            CALL PRB                                                                  
            CALL PRC ('  VALUES JUST OUTBOARD (Y>0)')                                       
            DO IZ = 0, PMIZS                                                       
               WRITE (7,9104) IZ,                                                      
     >            CPCHS(IXOUT,IY0 ,IZ,pz),100.0*CPRCS(IXOUT,IY0 ,IZ,pz),           
     >            CPCHS(IXOUT,IYL8,IZ,pz),100.0*CPRCS(IXOUT,IYL8,IZ,pz),           
     >            CPCHS(IXOUT,IYL4,IZ,pz),100.0*CPRCS(IXOUT,IYL4,IZ,pz)            
            end do
            CALL PRB                                                                  
            CALL PRC ('  VALUES JUST OUTBOARD (Y<0)')                                       
            DO IZ = 0, PMIZS                                                       
               WRITE (7,9104) IZ,                                                      
     >          CPCHS(IXOUT,IY0LT,IZ,pz),100.0*CPRCS(IXOUT,IY0LT,IZ,pz),           
     >          CPCHS(IXOUT,-IYL8,IZ,pz),100.0*CPRCS(IXOUT,-IYL8,IZ,pz),           
     >          CPCHS(IXOUT,-IYL4,IZ,pz),100.0*CPRCS(IXOUT,-IYL4,IZ,pz)            
            end do
C-----------------------------------------------------------------------
            IF (IQXFAC.LT.0) THEN                                                     
               CALL PRB                                                                
               CALL PRR ('  VALUES NEAR X = -LAMBDA (Y>0) = '
     >                    ,QXS(IQXFAC))                   
               DO IZ = 0, PMIZS                                                     
                  WRITE (7,9104) IZ,                                                    
     >            CPCHS(IXFAC,IY0 ,IZ,pz),100.0*CPRCS(IXFAC,IY0 ,IZ,pz),         
     >            CPCHS(IXFAC,IYL8,IZ,pz),100.0*CPRCS(IXFAC,IYL8,IZ,pz),         
     >            CPCHS(IXFAC,IYL4,IZ,pz),100.0*CPRCS(IXFAC,IYL4,IZ,pz)          
               end do
               CALL PRB                                                                
               CALL PRR ('  VALUES NEAR X = -LAMBDA (Y<0) = '
     >                     ,QXS(IQXFAC))                   
               DO IZ = 0, PMIZS                                                     
                  WRITE (7,9104) IZ,                                                    
     >                 CPCHS(IXFAC,IY0LT,IZ,pz)
     >                  ,100.0*CPRCS(IXFAC,IY0LT,IZ,pz),         
     >                 CPCHS(IXFAC,-IYL8,IZ,pz)
     >                  ,100.0*CPRCS(IXFAC,-IYL8,IZ,pz),         
     >                 CPCHS(IXFAC,-IYL4,IZ,pz)
     >                  ,100.0*CPRCS(IXFAC,-IYL4,IZ,pz)          
               end do

            ENDIF                                                                     
         end do              ! end of pz loop

      ENDIF                                                                     
C-----------------------------------------------------------------------
      CALL PRB                                                                  
      CALL PRC ('ELECTRIC FIELD AND DRIFT VELOCITY')                            
      CALL PRC ('  Y POSITION FACTORS  ELECTRIC    DRIFT')                      
      CALL PRC ('                        FIELD    VELOCITY')                    
      WRITE (7,9106) 'NEAR Y = 0      ', CEYS(IQY0  ),CVHYS(IQY0  )             
      IF ((CIOPTF.LE.5).OR.(CIOPTF.EQ.9)) THEN                                 
         WRITE (7,9106) 'NEAR Y = L/8    ', CEYS(IQYL8 ),CVHYS(IQYL8 )             
         WRITE (7,9106) 'NEAR Y = L/4    ', CEYS(IQYL4 ),CVHYS(IQYL4 )             
         WRITE (7,9106) 'NEAR Y = L/2    ', CEYS(IQYL2 ),CVHYS(IQYL2 )             
         WRITE (7,9106) 'NEAR Y = L      ', CEYS(IQYL  ),CVHYS(IQYL  )             
         WRITE (7,9106) 'NEAR Y = 3L/2   ', CEYS(IQY3L2),CVHYS(IQY3L2)             
         WRITE (7,9106) 'NEAR Y = 7L/4   ', CEYS(IQY7L4),CVHYS(IQY7L4)             
         WRITE (7,9106) 'NEAR Y =15L/8   ', CEYS(IQY158),CVHYS(IQY158)             
         WRITE (7,9106) 'NEAR Y = 2L     ', CEYS(IQY2L ),CVHYS(IQY2L )             
      ELSEIF (CIOPTF.EQ.6) THEN                                                 
         WRITE (7,9106) 'NEAR Y = YSTAG/2', CEYS(IQYS2 ),CVHYS(IQYS2 )             
         WRITE (7,9106) 'NEAR Y =3YSTAG/2', CEYS(IQY3S2),CVHYS(IQY3S2)             
         WRITE (7,9106) 'NEAR Y =5YSTAG/2', CEYS(IQY5S2),CVHYS(IQY5S2)             
      ELSEIF (CIOPTF.EQ.7) THEN                                                 
         WRITE (7,9106) 'NEAR Y = YSTAG/2', CEYS(IQYS2 ),CVHYS(IQYS2 )             
         WRITE (7,9106) 'NEAR Y = 5YSTAG ', CEYS(IQY5S ),CVHYS(IQY5S )             
         WRITE (7,9106) 'NEAR Y =15YSTAG ', CEYS(IQY15S),CVHYS(IQY15S)             
      ENDIF                                                                     
C-----------------------------------------------------------------------
      IF (CPRINT.EQ.1.or.cprint.eq.9) THEN

         do pz = 1,maxpzone
            CALL PRB                                                                  
            CALL PRC ('X POSITION FACTORS ')                                          
            CALL PRI ('ELECTRIC FIELD, IONIZATION STATE ', 1)                         
            CALL PRR (' JUST OUTBOARD   ',CFEXZS(IXOUT,1,1,pz)
     >            /QS(IQXOUT)**2)            
            CALL PRR (' NEAR X =-AW/4   ',CFEXZS(IXAW4,1,1,pz)
     >            /QS(IQXAW4)**2)            
            CALL PRR (' NEAR X =-AW/2   ',CFEXZS(IXAW2,1,1,pz)
     >            /QS(IQXAW2)**2)            
            CALL PRR (' NEAR X =-AW     ',CFEXZS(IXAW ,1,1,pz)
     >            /QS(IQXAW )**2)            
            IF (IQXFAC.LT.0)                                                          
     >           CALL PRR (' NEAR X =-LAMBDA ',CFEXZS(IXFAC,1,1,pz)
     >            /QS(IQXFAC)**2)            
            IF (MAXIZ.GT.1) THEN                                                      
               CALL PRI ('ELECTRIC FIELD, IONIZATION STATE', MAXIZ)                      
               CALL PRR (' JUST OUTBOARD   ',CFEXZS(IXOUT,1,MAXIZ,pz)
     >              /QS(IQXOUT)**2)        
               CALL PRR (' NEAR X =-AW/4   ',CFEXZS(IXAW4,1,MAXIZ,pz)
     >              /QS(IQXAW4)**2)        
               CALL PRR (' NEAR X =-AW/2   ',CFEXZS(IXAW2,1,MAXIZ,pz)
     >              /QS(IQXAW2)**2)        
               CALL PRR (' NEAR X =-AW     ',CFEXZS(IXAW ,1,MAXIZ,pz)
     >              /QS(IQXAW )**2)        
               IF (IQXFAC.LT.0)                                                          
     >              CALL PRR (' NEAR X =-LAMBDA ' 
     >              ,CFEXZS(IXFAC,1,MAXIZ,pz)
     >              /QS(IQXFAC)**2)        
            ENDIF                                                                     
C-----------------------------------------------------------------------
            CALL PRB                                                                  
            CALL PRC ('DRIFT VELOCITY, ALL IONIZATION STATES')                        
            CALL PRR (' JUST OUTBOARD   ',CFVHXS(IXOUT,1,pz)/QS(IQXOUT))                
            CALL PRR (' NEAR X =-AW/4   ',CFVHXS(IXAW4,1,pz)/QS(IQXAW4))                
            CALL PRR (' NEAR X =-AW/2   ',CFVHXS(IXAW2,1,pz)/QS(IQXAW2))                
            CALL PRR (' NEAR X =-AW     ',CFVHXS(IXAW ,1,pz)/QS(IQXAW ))                
            IF (IQXFAC.LT.0)                                                          
     >           CALL PRR (' NEAR X =-LAMBDA ', CFVHXS(IXFAC,1,pz)
     >                       /QS(IQXFAC))                

         end do                 ! end of pz loop
         
      ENDIF                                                                     
C-----------------------------------------------------------------------
C     
C---- MISCELLANEOUS EXTRA TABLES TO PRINT.  COMMENT OUT AS REQUIRED.            
C---- SET OUTPUT CHANNELS TO 6 OR 7 ALSO AS REQUIRED.                           
C     
C--   WRITE (6,9200) 'QXS          ',(QXS  (IQX),IQX=1-NQXSO,NQXSI,100)         
      WRITE (6,9200) 'CXBFS Y<0  ',(CXBFS(IQX,1),IQX=1-NQXSO,NQXSI,100)         
      WRITE (6,9200) 'CXBFS Y>0  ',(CXBFS(IQX,2),IQX=1-NQXSO,NQXSI,100)         
C--   WRITE (6,9200) 'CXAFS Y<0  ',(CXAFS(IQX,1),IQX=1-NQXSO,NQXSI,100)         
C--   WRITE (6,9200) 'CXAFS Y>0  ',(CXAFS(IQX,2),IQX=1-NQXSO,NQXSI,100)         
      do pz = 1,maxpzone
         CALL TAUPRA (6,CTEMBS(1,-maxnys,pz),
     >         'PLASMA TEMPERATURE (EV)',-1)              
         CALL TAUPRA (6,CTEMBSI(1,-maxnys,pz),
     >         'PLASMA ION TEMPS (EV) ',-1) 
         IF (NTIG.NE.0.OR.NTEG.NE.0) THEN
            CALL TAUPRA (7,CTEMBS(1,-maxnys,pz),
     >            'PLASMA TEMPERATURE (E)',-1)              
            CALL TAUPRA (7,CTEMBSI(1,-maxnys,pz),
     >            'PLASMA ION TEMPS (EV)',-1) 
         ENDIF
         CALL TAUPRA (6,CRNBS(1,-maxnys,pz),
     >         'PLASMA DENSITY    (M**-3)',-1)                    
      end do
      
C--   WRITE (6,9200) 'TEMPERATURE <',(QTEMBS(IQX,1),IQX=1-NQXSO,1,100)          
C--   WRITE (6,9200) 'TEMPERATURE >',(QTEMBS(IQX,2),IQX=1-NQXSO,1,100)          
C--   WRITE (6,9200) 'DENSITY     <',(QRNBS (IQX,1),IQX=1-NQXSO,1,100)          
C--   WRITE (6,9200) 'DENSITY     >',(QRNBS (IQX,2),IQX=1-NQXSO,1,100)          
C--   CALL TAUPRA (6,CNHS  ,'NEUTRAL HYDROGEN   (M**-3)',-1)                    
C--   DO 505 IZ = 0, NIZS, 2                                                    
C--   CALL TAUPRA (6,CFIZS(1,-MAXNYS,IZ),'IONISATION TIMES (S)  ',IZ)           
C--   CALL TAUPRA (6,CFRCS(1,-MAXNYS,IZ),'E-I RECOMB TIMES (S)  ',IZ)           
C--   CALL TAUPRA (6,CFCXS(1,-MAXNYS,IZ),'C-X AND E-I RECOMB (S)',IZ)           
C--   CALL TAUPRA (6,CPCHS(1,-MAXNYS,IZ),'CHANGE OF STATE PROBS ',IZ)           
C     -505 CONTINUE                                                                  
C--   DO 510 IZ = 1, NIZS                                                       
C--   CALL TAUPRA (6,CFEXZS(1,-MAXNYS,IZ),'ELECTRIC FIELD FACTOR',IZ)           
C     -510 CONTINUE                                                                  
C--   CALL TAUPRA (6,CFVHXS              ,'DRIFT VELOCITY FACTOR',-1)           
C-----------------------------------------------------------------------
 9003 FORMAT(1X,'  ',A15,1P,' +/-',G11.2,' +',G11.2,',  PROB',0P,F8.5)          
 9101 FORMAT( 1X,'       IONISATION  E-I RECOMB  MODIFIED RECOM',               
     >     'BINATION TIMES DUE           ',           
     >     /1X,'        TIME FOR    TIME FOR   TO CHARGE EXCH',               
     >     'ANGE RECOMBINATION           ',           
     >     /1X,'  IZ   IZ -> IZ+1  IZ -> IZ-1  NEAR Y = 0  NE',               
     >     'AR Y=L/8 NEAR Y=L/4          ')           
 9102 FORMAT(1X,I4,1P,5(1X,G11.2))                                              
 9103 FORMAT( 1X,'      PROBABILITY  RECOMB PROBABILITY  RECOMB ',              
     >     'PROBABILITY RECOMB       ',            
     >     /1X,'  IZ  NEAR Y = 0    FRAC  NEAR Y=L/8    FRAC  ',              
     >     'NEAR Y=L/4   FRAC        ')            
 9104 FORMAT(1X,I4,3(1X,1P,G11.2,0P,F7.2,'%'))                                  
 9106 FORMAT(2X,A,1P,5(1X,G11.2))                                               
 9107 FORMAT(2X,A,12X,2F8.1,1X,1P,2G9.2,0P,F8.1)                                
 9108 FORMAT(2X,A,1X,2F6.3,F7.1,F8.1,1X,1P,2G9.2,0P,F8.1)                       
 9109 FORMAT(2X,A,12X,2F8.1,1X,1P,18X,0P,F8.1)             
 9110 FORMAT(2X,A,1X,2F6.3,F7.1,F8.1,1X,1P,18X,0P,F8.1)                       
 9200 FORMAT(/1X,A,1P,/,(1X,14G9.2))                                            
C     
 999  RETURN                                                                    
      END                                                                       
C     
C     
C     
      SUBROUTINE TAUPR2 (QTIM, NIZS)                                            
      use mod_params
      use mod_comt2
      use mod_comtor
      use mod_comtau
      use mod_comxyt
      use mod_printr
      use mod_global_options
      IMPLICIT  none
      INTEGER   NIZS                                                            
      REAL      QTIM                                                            
C     
C***********************************************************************
C     
C     PRINTS REPRESENTATIVE VALUES OF CHARACTERISTIC TIMES, DERIVED           
C     FROM FACTORS CFPS,CFSS,CFTS SET UP IN COMMON COMT2, USING THE             
C     INITIAL ION TEMPERATURE RETURNED FROM NEUT/LAUNCH.                        
C     
C     NOTE "TAUPRF" IS INCLUDED IN MEMBER "RUNLIM3" SINCE IT USES F77           
C     CHARACTER FEATURES INCOMPATIBLE WITH THE CRAY CFT2 COMPILER.              
C     HENCE THE REST OF TAU CAN BE COMPILED WITH CFT2.                          
C     
C     CHRIS FARRELL    DECEMBER 1987                        
C     
C***********************************************************************
C     
c     INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c     INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c     INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c     INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c     INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c     INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c     
c     include   'global_options'
C     
      INTEGER   MAXIZ,IZ                                                        
C-----------------------------------------------------------------------
      MAXIZ  = MIN (NIZS, CION)                                                 
      IF (CPRINT.EQ.0) GOTO 999                                                 
      CALL PRB                                                                  
      CALL PRC ('SAMPLE CHARACTERISTIC TIMES BASED ON INITIAL ION TEMPER        
     >     ATURE')                                                                   
      call PRI ('DATA PRINTED FOR PLASMA ZONE NEAR P=0: PZ=',pzones(1))
      CALL PRC ('                        TAU PARALLEL  TAU STOPPING   TA        
     >     U HEATING')                                                               
      CALL PRI ('IONISATION STATE', 1)                                          
      CALL TAUPRF ('NEAR X = A     ',IXA,  1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('NEAR X = A/2   ',IXA2, 1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('NEAR X = A/4   ',IXA4, 1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('JUST INBOARD   ',IXIN, 1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('JUST OUTBOARD  ',IXOUT,1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('NEAR X =-AW/4  ',IXAW4,1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('NEAR X =-AW/2  ',IXAW2,1,1,CTEMSC,QTIM,CIOPTC)              
      CALL TAUPRF ('NEAR X =-AW    ',IXAW, 1,1,CTEMSC,QTIM,CIOPTC)              
      IF (IQXFAC.LT.0)                                                          
     >     CALL TAUPRF ('NEAR X =-LAMBDA',IXFAC,1,1,CTEMSC,QTIM,CIOPTC)              
C     
      IF (MAXIZ.GT.1) THEN                                                      
        CALL PRB                                                                  
        CALL PRI ('IONISATION STATE', MAXIZ)                                      
        CALL TAUPRF ('NEAR X = A     ',IXA,  1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('NEAR X = A/2   ',IXA2, 1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('NEAR X = A/4   ',IXA4, 1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('JUST INBOARD   ',IXIN, 1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('JUST OUTBOARD  ',IXOUT,1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('NEAR X =-AW/4  ',IXAW4,1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('NEAR X =-AW/2  ',IXAW2,1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        CALL TAUPRF ('NEAR X =-AW    ',IXAW, 1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
        IF (IQXFAC.LT.0)                                                          
     >       CALL TAUPRF ('NEAR X =-LAMBDA',IXFAC,1,MAXIZ,CTEMSC,QTIM
     >                   ,CIOPTC)          
      ENDIF                                                                     
C     
C     DO 500 IZ = 1, MAXIZ                                                      
C     CALL TAUPRA (6,CFPS  (1,-MAXNYS,IZ),'CFPS FACTORS',IZ)                  
C     CALL TAUPRA (6,CFSS  (1,-MAXNYS,IZ),'CFSS FACTORS',IZ)                  
C     CALL TAUPRA (6,CFTS  (1,-MAXNYS,IZ),'CFTS FACTORS',IZ)                  
C     CALL TAUPRA (6,CCCFPS(1,-MAXNYS,IZ),'CCCFPS FACTORS',IZ)                
C     500 CONTINUE                                                                  
C     
 999  RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE TAUFIX (IX,TEMOLD,TEMNEW)                                      
      use mod_params
      use mod_comt2
      use mod_comtor
      use mod_comtau
      use mod_comxyt
      use mod_lambda
      IMPLICIT none
      INTEGER IX                                                                
      REAL    TEMOLD,TEMNEW                                                     
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *   TAUFIX:  USED TO QUICKLY ADJUST RELEVANT CHARACTERISTIC TIMES   *        
C  *   FOR NON-STANDARD PLASMA OPTIONS.  BECAUSE OF THE SQUARE ROOTS   *        
C  *   THIS ROUTINE WANTS TO BE CALLED SPARINGLY TO PREVENT DRAMATIC   *        
C  *   SLOWING OF PROGRAM EXECUTION - FOR EXAMPLE, WHENEVER THE        *        
C  *   TEMPERATURE CHANGES BY 10% OR MORE.                             *        
C  *   CALLED FROM LIM3.                                               *        
C  *                                                                   *        
C  *                        CHRIS FARRELL (HUNTERSKIL)  SEPT 1988      *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
C                                                                               
c      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
c      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
c      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
c      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      REAL :: RIZSQR,RATIO1,RATIO2,TAU                                             
      INTEGER :: IQX,IY                                                            
      real :: lambda
      integer :: pz
C     
C     WRITE (6,9001) TEMOLD,                                                    
C    >  CFPS(IX,1,CIZ),CCCFPS(IX,1,CIZ),CFSS(IX,1,CIZ),CFTS(IX,1,CIZ)           
C                                                                               
      do pz = 1,maxpzone

       IF (CIOPTB.EQ.2) THEN                                                     
        RATIO1 = SQRT (TEMOLD / TEMNEW)                                         
        RATIO2 = 1.0 / SQRT (RATIO1)                                            
        DO 100 IY = -NYS, NYS                                                   
          CFPS(IX,IY,CIZ,pz)   = RATIO1 * CFPS(IX,IY,CIZ,pz)                          
          CCCFPS(IX,IY,CIZ,pz) = RATIO2 * CCCFPS(IX,IY,CIZ,pz)                        
  100   CONTINUE                                                                
C                                                                               
      ELSEIF (CIOPTB.EQ.3) THEN                                                 
        RATIO2 = SQRT (TEMNEW / TEMOLD)                                         
        DO 110 IY = -NYS, NYS                                                   
          CCCFPS(IX,IY,CIZ,pz) = RATIO2 * CCCFPS(IX,IY,CIZ,pz)                        
  110   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      IF (CIOPTC.EQ.2) THEN                                                     
        DO 200 IY = -NYS, NYS                                                   
          TAU = CFPS(IX,IY,CIZ,pz) / (2.0*TEMNEW)                                  
          IF (TAU.GT.1.E-3) THEN                                                
            CFSS(IX,IY,CIZ,pz) = 1.0 - EXP(-TAU)                                   
          ELSE                                                                  
            CFSS(IX,IY,CIZ,pz) = TAU                                               
          ENDIF                                                                 
  200   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      IF (CIOPTD.EQ.3) THEN                                                     
        RIZSQR = REAL(CIZ) * REAL(CIZ)                                          
        IQX = IQXS(IX)                                                          
        C215B = (CRMI * CTBI + CRMB * TEMNEW) ** 1.5                            
        DO 300 IY = -NYS, NYS                                                   
          IF (CTBI.LE.0.0)                                                      
     >      C215B = (CRMI * CTEMBSI(IX,IY,pz) + CRMB * TEMNEW)** 1.5             

            ! jdemod - if lambda is spatially varying change it in C215A
            !   calculation of tau
            if (lambda_vary_opt.eq.1) then
              lambda = coulomb_lambda(crnbs(ix,iy,pz),ctembsi(ix,iy,pz))
            else
               lambda = 1.0
            endif
          
          TAU = C215A * RIZSQR*QS(IQX)*CRNBS(IX,IY,pz) / C215B * LAMBDA
          IF (TAU.GT.1.E-3) THEN                                                
            CFTS(IX,IY,CIZ,pz) = 1.0 - EXP (-TAU)                                  
          ELSE                                                                  
            CFTS(IX,IY,CIZ,pz) = TAU                                               
          ENDIF                                                                 
  300   CONTINUE                                                                
      ENDIF                                                                     

      end do ! end of pz loop
      
C


      
C     WRITE (6,9002) TEMNEW,                                                    
C    >  CFPS(IX,1,CIZ),CCCFPS(IX,1,CIZ),CFSS(IX,1,CIZ),CFTS(IX,1,CIZ)           
C                                                                               
 9001 FORMAT(10X,'TAUFIX: TEMOLD=',F7.2,' OLDVALS=',1P,4G11.4)                  
 9002 FORMAT(10X,'        TEMNEW=',F7.2,' NEWVALS=',1P,4G11.4)                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE TAUPRA (JW,AS,NAME,ISTATE)                                     
      use mod_params
      use mod_printr
      IMPLICIT none
c      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      REAL AS(MAXNXS,-MAXNYS:MAXNYS)                                            
      INTEGER JW,ISTATE                                                         
      CHARACTER*(*) NAME                                                        
c      INCLUDE 'printr'                                                          
C     INCLUDE (PRINTR)                                                          
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *   TAUPRA:  ROUTINE TO PRINT A 2-D ARRAY GIVEN IN ARGUMENT LIST    *        
C  *   AT A SET OF X AND Y POSITIONS.                                  *        
C  *                                                                   *        
C  *                        CHRIS FARRELL (HUNTERSKIL)  AUGUST 1988    *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      IF (ISTATE.GE.0) THEN                                                     
        WRITE (JW,9104) NAME,ISTATE                                             
      ELSE                                                                      
        WRITE (JW,9104) NAME                                                    
      ENDIF                                                                     
      WRITE (JW,9105)                                                           
C                                                                               
      WRITE (JW,9106) 'X = A   ',AS(IXA  ,-IYL8),AS(IXA  ,IY0),                 
     >  AS(IXA  ,IYL4),AS(IXA  ,IYL2),AS(IXA  ,IYL),AS(IXA  ,IY2L)              
      WRITE (JW,9106) 'X = A/2 ',AS(IXA2 ,-IYL8),AS(IXA2 ,IY0),                 
     >  AS(IXA2 ,IYL4),AS(IXA2 ,IYL2),AS(IXA2 ,IYL),AS(IXA2 ,IY2L)              
      WRITE (JW,9106) 'X = A/4 ',AS(IXA4 ,-IYL8),AS(IXA4 ,IY0),                 
     >  AS(IXA4 ,IYL4),AS(IXA4 ,IYL2),AS(IXA4 ,IYL),AS(IXA4 ,IY2L)              
      WRITE (JW,9106) 'INBOARD ',AS(IXIN ,-IYL8),AS(IXIN ,IY0),                 
     >  AS(IXIN ,IYL4),AS(IXIN ,IYL2),AS(IXIN ,IYL),AS(IXIN, IY2L)              
      WRITE (JW,9106) 'OUTBOARD',AS(IXOUT,-IYL8),AS(IXOUT,IY0),                 
     >  AS(IXOUT,IYL4),AS(IXOUT,IYL2),AS(IXOUT,IYL),AS(IXOUT,IY2L)              
      WRITE (JW,9106) 'X =-AW/4',AS(IXAW4,-IYL8),AS(IXAW4,IY0),                 
     >  AS(IXAW4,IYL4),AS(IXAW4,IYL2),AS(IXAW4,IYL),AS(IXAW4,IY2L)              
      WRITE (JW,9106) 'X =-AW/2',AS(IXAW2,-IYL8),AS(IXAW2,IY0),                 
     >  AS(IXAW2,IYL4),AS(IXAW2,IYL2),AS(IXAW2,IYL),AS(IXAW2,IY2L)              
      WRITE (JW,9106) 'X =-AW  ',AS(IXAW ,-IYL8),AS(IXAW ,IY0),                 
     >  AS(IXAW ,IYL4),AS(IXAW ,IYL2),AS(IXAW ,IYL),AS(IXAW ,IY2L)              
      IF (IQXFAC.LT.0)                                                          
     >WRITE (JW,9106) 'X =-LAM ',AS(IXFAC,-IYL8),AS(IXFAC,IY0),                 
     >  AS(IXFAC,IYL4),AS(IXFAC,IYL2),AS(IXFAC,IYL),AS(IXFAC,IY2L)              
C                                                                               
 9104 FORMAT(/1X,A,:,'  FOR STATE',I2)                                          
 9105 FORMAT(10X,'  Y=-L/8   Y =+0    Y=L/4    Y=L/2    Y = L    Y =2L')        
 9106 FORMAT(2X,A8,1P,6G9.2)                                                    
      RETURN                                                                    
      END                                                                       


      
C     
C     
      SUBROUTINE TAUIN3 (QTIM,NIZS,DEFACT)                                      
      use mod_params
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
      use mod_comtor
      use mod_comtau
      use mod_comxyt
      use mod_coords
      use mod_lambda
      IMPLICIT  none
      REAL      QTIM                                                            
      INTEGER   NIZS                                                            
      DOUBLE PRECISION DEFACT                                                   
C     
C     *********************************************************************        
C     *                                                                   *        
C     *  TAUIN3:  THIS ROUTINE RECALCULATES THE CHARACTERISTIC TIMES      *        
C     *  ACCORDING TO THE STEADY STATE ION DISTRIBUTIONS AT THE END OF    *        
C     *  A LIM RUN.  THESE TIMES CAN THEN BE USED FOR ANOTHER ITERATION   *        
C     *  OF LIM, UNTIL EVENTUALLY WE MAY ACHIEVE A SELF CONSISTENT PLASMA *        
C     *  AS DESCRIBED IN NOTE 207.                                        *        
C     *  NOTE THAT TO PREVENT DIVIDE BY ZEROES IN INNER LOOPS THIS        *        
C     *  ROUTINE HAS TO BE COMPILED WITH VECTORISATION OFF.               *        
C     *  NOTE 297: ZEFFS RELATED QUANTITIES ARE NOW MULTIPLIED BY         *        
C     *  SIN(THETAB) AND THE RATIO LP+/LP-                                *        
C     *                                                                   *        
C     *                          CHRIS FARRELL  (HUNTERSKIL)  SEPT 1988   *        
C     *                                                                   *        
C     *********************************************************************        
C     
c     INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c     INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c     INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c     INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c     INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
c     INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c     INCLUDE   'dynam1'                                                        
C     INCLUDE   (DYNAM1)                                                        
c     INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
C     
      INTEGER   IX,IY,IZ,IQX,JZ                                                 
      REAL      LAMBDA,TEMP,SQRTMB,SQRTMI,TOTALP,ARG,RIZB,TOTALT                
      REAL      TAUP(0:MAXIZS),TAUS(0:MAXIZS),TAUT(0:MAXIZS),TOTALS,RIZ         

      integer :: pz
      pz = 1
      
c     slmod
c     PARAMETER (LAMBDA=15.0)                                                   
C     
c     IF (CIOPTE.EQ.10) THEN
c     lambda = 17.3 - .5*ln(n/1e20) + 1.5*ln(T/1000.)
c     
c     jdemod - replace with global lambda options
c     LIM baseline is lambda_opt = 2      
      if (lambda_vary_opt.eq.0) then
         lambda = coulomb_lambda(cnbin,ctibin)
      else
         lambda = 1.0
      endif


c     LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)
c     ELSE
c     LAMBDA = 15.0
c     ENDIF
c     slmod end
C     
      SQRTMB = SQRT (CRMB)                                                      
      SQRTMI = SQRT (CRMI)                                                      
      RIZB   = REAL (CIZB)                                                      
C     
c     
c     code uses ddlims and is not 3D poloidal zone compatible     
      pz = 1
!     do pz = 1,maxpzone

      DO 530 IZ = 1, NIZS                                                       
         WRITE (6,9001) IZ,(JZ,JZ=1,9)                                            
         RIZ = REAL (IZ)                                                          
         DO 520 IY = -NYS, NYS                                                    
            IF (IY.EQ.0) GOTO 520                                                   
            DO 510 IX = 1, NXS                                                      
               IF (CTBI.GT.0.0) THEN                                                  
                  TEMP = CTBI                                                          
               ELSE                                                                   
                  TEMP = CTEMBSI(IX,IY,pz)                                             
               ENDIF                                                                  
               IQX = IQXS(IX)                                                         
C_______________________________________________________________________
C     
C     TAU PARALLEL     CFPS = 2.DELTAT.TI/TAUPARA                            
C_______________________________________________________________________
C     

! jdemod - if lambda is spatially varying change it in C215A
!   calculation of tau
               if (lambda_vary_opt.eq.1) then
                  lambda = coulomb_lambda(crnbs(ix,iy,pz)
     >                                   ,ctembsi(ix,iy,pz))
               else
                  lambda = 1.0
               endif

               TOTALP = 0.0                                                           
               IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
                  TAUP(0) = CRMI * SQRT(TEMP) /                                        
     >                 (6.8E-14 * SQRTMB * ZEFFS(IX,IY,5) *                                
     >                 RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
               ELSE                                                                   
                  TAUP(0) = 0.0                                                        
               ENDIF                                                                  
               IF (TAUP(0).GT.0.0) TOTALP = 1.0 / TAUP(0)                             
C     
               DO 100 JZ = 1, NIZS                                                    
                  IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
                     TAUP(JZ) = CRMI * SQRT(SNGL(DDTS(IX,IY,JZ))) /                     
     >                    (6.8E-14*SQRTMI*SNGL(DDLIMS(IX,IY,JZ)*DEFACT)*               
     >                    CSINTB*REAL(JZ)**2 *RIZ*RIZ*LAMBDA*CZENH)               
                  ELSE                                                                 
                     TAUP(JZ) = 0.0                                                     
                  ENDIF                                                                
                  IF (TAUP(JZ).GT.0.0) TOTALP = TOTALP + 1.0 / TAUP(JZ)                
 100           CONTINUE                                                               
C     
               CFPS(IX,IY,IZ,pz) = 2.0 * QTIM * QS(IQX) * TOTALP                         
C     
               IF (CFPS(IX,IY,IZ,pz).LE.0.0) THEN                                        
                  CCCFPS(IX,IY,IZ,pz) = 0.0                                               
               ELSE                                                                   
                  CCCFPS(IX,IY,IZ,pz) = SQRT(4.88E8
     >                  /(CFPS(IX,IY,IZ,pz)*CRMI))*QTIM * QS(IQX)                                    
               ENDIF                                                                  
C_______________________________________________________________________
C     
C     TAU STOPPING     CFSS = 1 - EXP (-DELTAT/TAUSTOP)                      
C_______________________________________________________________________
C     
               TOTALS = 0.0                                                           
               IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
                  TAUS(0) = CRMI * (TEMP**1.5) /                                       
     >                 (6.8E-14*SQRTMB*(1.0+CRMB/CRMI)*ZEFFS(IX,IY,5) *              
     >                 RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
               ELSE                                                                   
                  TAUS(0) = 0.0                                                        
               ENDIF                                                                  
               IF (TAUS(0).GT.0.0) TOTALS = 1.0 / TAUS(0)                             
C     
               DO 200 JZ = 1, NIZS                                                    
                  IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
                     TAUS(JZ) = CRMI * (SNGL(DDTS(IX,IY,JZ))**1.5) /                    
     >                    (6.8E-14 * SQRTMI * 2.0 *
     >                     SNGL(DDLIMS(IX,IY,JZ)*DEFACT) *         
     >                    CSINTB * REAL(JZ)**2 * RIZ*RIZ*LAMBDA*CZENH)               
                  ELSE                                                                 
                     TAUS(JZ) = 0.0                                                     
                  ENDIF                                                                
                  IF (TAUS(JZ).GT.0.0) TOTALS = TOTALS + 1.0 / TAUS(JZ)                
 200           CONTINUE                                                               
C     
               ARG = QTIM * QS(IQX) * TOTALS                                          
               IF (ARG.GT.1.E-3) THEN                                                 
                  CFSS(IX,IY,IZ,pz) = 1.0 - EXP(-ARG)                                     
               ELSE                                                                   
                  CFSS(IX,IY,IZ,pz) = ARG                                                 
               ENDIF                                                                  
C_______________________________________________________________________
C     
C     TAU HEATING      CFTS = 1 - EXP (-DELTAT/TAUHEAT)                      
C_______________________________________________________________________
C     
               TOTALT = 0.0                                                           
               IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
                  TAUT(0) = CRMI * (TEMP**1.5) /                                       
     >                 (1.4E-13 * SQRTMB * ZEFFS(IX,IY,5) *                                
     >                 RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
               ELSE                                                                   
                  TAUT(0) = 0.0                                                        
               ENDIF                                                                  
               IF (TAUT(0).GT.0.0) TOTALT = 1.0 / TAUT(0)                             
C     
               DO 300 JZ = 1, NIZS                                                    
                  IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
                     TAUT(JZ) = CRMI * (SNGL(DDTS(IX,IY,JZ))**1.5) /                    
     >                    (1.4E-13 * SQRTMI *
     >                     SNGL(DDLIMS(IX,IY,JZ)*DEFACT) *               
     >                    CSINTB * REAL(JZ)**2 * RIZ*RIZ*LAMBDA*CZENH)               
                  ELSE                                                                 
                     TAUT(JZ) = 0.0                                                     
                  ENDIF                                                                
                  IF (TAUT(JZ).GT.0.0) TOTALT = TOTALT + 1.0 / TAUT(JZ)                
 300           CONTINUE                                                               
C     
               ARG = QTIM * QS(IQX) * TOTALT                                          
               IF (ARG.GT.1.E-3) THEN                                                 
                  CFTS(IX,IY,IZ,pz) = 1.0 - EXP(-ARG)                                     
               ELSE                                                                   
                  CFTS(IX,IY,IZ,pz) = ARG                                                 
               ENDIF                                                                  
C     
               IF (10*(IY/10).EQ.IABS(IY) .AND. IY.LE.50 .AND.                        
     >              10*(IX/10).EQ.IX .AND. 2*(IZ/2).NE.IZ) THEN                        
                  WRITE (6,9002) IY,IX,1.0/TOTALP,(TAUP(JZ),JZ=0,NIZS)                  
                  WRITE (6,9003)       1.0/TOTALS,(TAUS(JZ),JZ=0,NIZS)                  
                  WRITE (6,9004)       1.0/TOTALT,(TAUT(JZ),JZ=0,NIZS)                  
               ENDIF                                                                  
 510        CONTINUE                                                                
 520     CONTINUE                                                                 
 530  CONTINUE                                                                  

!     end do                    ! end of pz loop
      
C     
      CALL TAUPRA (7,ZEFFS(1,-MAXNYS,5),'ZB.NBT (NT BASED) (M**-3)',-1)         
      RETURN                                                                    
 9001 FORMAT(//1X,'TAUIN3: IONIZATION STATE',I3,///1X,                          
     >     '  IY  IX        TOTAL    TAUB',9(5X,'TAU',I1),/1X,130('-'))            
 9002 FORMAT(1X,2I4,' PARL',1P,10E9.2)                                          
 9003 FORMAT(9X,    ' STOP',1P,10E9.2)                                          
 9004 FORMAT(9X,    ' HEAT',1P,10E9.2)                                          
      END                                                                       
C     
C     
C     
      SUBROUTINE TAUPRF (STRING, IX, IY, IZ, CTEMSC, QTIM, CIOPTC)              
      use mod_params
      use mod_comt2
      use mod_comxyt
      IMPLICIT none
      CHARACTER*15 STRING                                                       
      INTEGER      IX,IY,IZ,CIOPTC,IQX                                          
      REAL         CTEMSC,QTIM                                                  

      integer :: pz
C     
C***********************************************************************
C     SHORT ROUTINE TO PRINT A LINE OF THE "CHARACTERISTIC TIMES" TABLE           
C     OF TAU PARA, TAU STOP AND TAU HEAT.                                         
C     CHRIS FARRELL     DECEMBER 1987                      
C***********************************************************************
C     
c     INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c     INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
c     INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      CHARACTER*11 WTAUP,WTAUS,WTAUT                                            
      REAL         TAUP,TAUS,TAUT                                               
C     
C     WRITE (6,'('' TP,TS,TT='',3G11.4)') CFPS(IX,1,IZ),                        
C     >  CFSS(IX,1,IZ),CFTS(IX,1,IZ)                                             
C     
c     Use poloidal zone near P=0
      pz = pzones(1)
      IQX = IQXS(IX)                                                            
C     
      IF (CFPS(IX,IY,IZ,pz).GT.0.0) THEN                                           
         TAUP = 2.0 * CTEMSC * QTIM * QS(IQX) / CFPS(IX,IY,IZ,pz)                   
         WRITE (WTAUP,'(1P,G11.2)') TAUP                                         
      ELSE                                                                      
         WTAUP = ' INFINITE  '                                                   
      ENDIF                                                                     
C     
      IF     (CIOPTC.EQ.2) THEN                                                 
         WTAUS = WTAUP                                                           
      ELSEIF (CFSS(IX,IY,IZ,pz).GE.1.0) THEN                                       
         WTAUS = ' INSTANT   '                                                   
      ELSEIF (CFSS(IX,IY,IZ,pz).GT.1.E-3) THEN                                     
         TAUS = -QTIM * QS(IQX) / LOG (1.0 - CFSS(IX,IY,IZ,pz))                     
         WRITE (WTAUS,'(1P,G11.2)') TAUS                                         
      ELSEIF (CFSS(IX,IY,IZ,pz).GT.0.0) THEN                                       
         TAUS = QTIM * QS(IQX) / CFSS(IX,IY,IZ,pz)                                  
         WRITE (WTAUS,'(1P,G11.2)') TAUS                                         
      ELSE                                                                      
         WTAUS = ' INFINITE  '                                                   
      ENDIF                                                                     
C     
      IF     (CFTS(IX,IY,IZ,pz).GE.1.0) THEN                                       
         WTAUT = ' INSTANT   '                                                   
      ELSEIF (CFTS(IX,IY,IZ,pz).GT.1.E-3) THEN                                     
         TAUT = -QTIM * QS(IQX) / LOG (1.0 - CFTS(IX,IY,IZ,pz))                     
         WRITE (WTAUT,'(1P,G11.2)') TAUT                                         
      ELSEIF (CFTS(IX,IY,IZ,pz).GT.0.0) THEN                                       
         TAUT = QTIM * QS(IQX) / CFTS(IX,IY,IZ,pz)                                  
         WRITE (WTAUT,'(1P,G11.2)') TAUT                                         
      ELSE                                                                      
         WTAUT = ' INFINITE  '                                                   
      ENDIF                                                                     
C     
      WRITE (7,9003) STRING,WTAUP,WTAUS,WTAUT                                   
 9003 FORMAT(5X,A15,5X,A11,3X,A11,3X,A11)                                       
      RETURN                                                                    
      END                                                                       
c     
c     
c     
      subroutine plasma_overlay(qtim)
      use mod_params
      use error_handling
      use mod_global_options
      use mod_soledge_input
      use mod_sol22_input_lim
      use mod_sol22_lim
      use mod_solcommon
      use mod_comxyt
      use mod_comt2
      use yreflection
      use mod_vtig
      use mod_comtor
      use mod_assign_plasma
      implicit none
      real :: qtim

c     
      integer :: ix,iy,iqx,iqy
      real :: t0,vel,n,x,y,y0,ti
      integer :: pz,yz
      integer :: ixout
      integer,external :: ipos

      integer :: solver_opt
      integer :: ir, new_unit,ierr,pz1,pz2,ixs,ixe,pzs,pze
      real :: xbnd1,xbnd2
      
      call pr_trace('TAU:plasma_overlay','Start')

      ! calculate X index of the separatrix
      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
      
!     set quantities needed in solvers that won't be passed as arguments down the call stack
      call setup_solvers(qtim,crmb,cizb,yscale,ixout,cprint)
      
!     the setup_vtig routine assigns masses and calculates the integration constant
!     and should be called for all vtig options - this is needed to calculate an estimate
!     of vtig from the temperature gradients and should be called in all cases 
      call setup_vtig(crmb,crmi,cnbin,ctibin)
c     
      
      call pr_trace('TAU:plasma_overlay:soledge_opt=',soledge_opt)
      
c     
c     If the collector probe 3D plamsa options are in effect then call the
c     code to set up the modified plamsa, efield and plasma velocity arrays
c     
c     sazmod - Maybe use a separate switch for this statement to allow
c     only setting up forces in lim3.f without prescribing a 
c     complex SOL (like SOL12, 13, etc.). 
c     
c      if (soledge_opt.eq.1.and.colprobe3d.eq.1) then 
      if (soledge_opt.eq.1) then 

!     plasma is calculated from lower absorbing surface to
!     upper absorbing surface - this allows for
!     asymmetric placement of the probe
!     call init_soledge(yabsorb1a,yabsorb2a)
!     call init_soledge(-cl,cl)
         
!     if (vary_absorb.eq.1) then
!     Find the x index where the step happens. I think this is IPOS?
!     ix_step1 = ipos(xabsorb1a_step, xs, nxs-1)
!     ix_step2 = ipos(xabsorb2a_step, xs, nxs-1)
!     write(0,*) 'ix_step1 = ',ix_step1,'(x = ',xs(ix_step1),')'
!     write(0,*) 'ix_step2 = ',ix_step2,'(x = ',xs(ix_step2),')'
         
!     Call soledge for the plasma from the wall to the step.
!     call soledge(1, ix_step1, qtim)
         
!     Call soledge for the plasma from the step to the top.
!     write(0,*) 'second soledge call'
!     call soledge(ix_step1+1, nxs, qtim)
         
!     deallocate storage here instead of inside soledge code.
!     call end_soledge
         
!     else
!     Just do the normal option with one absorbing wall.
!     call soledge(1,nxs,qtim)
!     call soledge(1,nxs/2,qtim)
!     endif

!     The default soledge option applies it to the entire plasma region
!     select the two point model solver using last argument 0
         call plasma_solver(1,nxs,1,maxpzone,0)

      endif



      call pr_trace('TAU:plasma_overlay:sol22_opt=',sol22_opt)

      
      if (sol22_opt.eq.1) then 
c     
c     This code calculates the plasma conditions for sections of the simulation
c     volume using SOL22. There are several scenarios.
c     
c     1) No absorbing surfaces - standard limiter or probe simulations
c     SOL22 is used to calculate the background plasma on field lines
c     that connect to the probe/limiter. (i.e. PZONE = 1)
c     
c     
!     Initialize some output options in SOL22 using values from slcom
!     call init_solcommon(0,0)

!     call sol22

!     Loop through the SOL22 input blocks at this level so that parameter files
!     only need to be loaded once



         if (nsol22_opt.gt.0) then 

            
            do ir = 1,nsol22_opt

!     jdemod - customized parameters can't be loaded here since storage has not
!              yet been allocated - delay this until after storage allocation.              
!     This requires storing the filename to a location that can be
!     accessed later               
!
!     call load_sol22_parameter_file(sol22_filenames(ir),ierr)

               call set_sol22_parameter_file(sol22_filenames(ir))

               
               xbnd1 = sol22_regions(ir,1)
               xbnd2 = sol22_regions(ir,2)

!     ignore pbnds for now
!     create two zones to overlay SOL22
!     pzone = 1 -> |P| =< CPC0 (width of probe)
!     pzone = 2 -> |P| > CPC0
!     CPC0 
!     if maxpzones = 1 OR not 3D (MAXNPS =1) then only
!     calculate for pzone = 1
!     
!     if there are absorbing surfaces - use them
!     - calculate plasma over [yabsorb1a,yabsorb2a]
!     - NOTE: there will need to be duplication of
!     all plasma arrays for 1:pzone areas
!     e.g. ctembs(ix,iz,pzone)
!     ctegs ... and forces by iz and pzone
!     
!     if not ignore
!     - calculate plasma over whole range and duplicate
!     - plasma from limiter surface to limiter surface
!     - SOL22 is NOT applied to pzone 2 since there are
!     no surfaces there.
!     - This means this is only applicable for non-3D cases 
!     
!     
               pz1 = sol22_regions(ir,3)
               pz2 = sol22_regions(ir,4)

!     SOL22 input specifies which poloidal zone to use the model 
!     pzs = sol22_regions(ir,3)
!     pze = sol22_regions(ir,4)

               write(0,'(a,1x,g12.5,a,1x,g12.5,a,1x,i8,a,1x,i8,a,a)')
     >           'Solver applied from X=',xbnd1,' to X=',xbnd2,
     >           ' for Zone 1=',pz1,' to Zone 2=',pz2,
     >           ' using SOL22 option file:',trim(sol22_filenames(ir))

               ixs=ipos(xbnd1,xouts,nxs)

               ! ipos returns the last element even if the value is greater than the maximum value of the
               ! array - however, to avoid overlapping solver ranges we need to exclude the bin that is greater
               ! than the xbnd2 not include it - except if it is the last bin
               if (xbnd2.gt.xouts(nxs)) then
                  ixe = nxs
               else
                  ixe=ipos(xbnd2,xouts,nxs)-1
               endif
                  
               pzs=pz1
               pze=pz2
               ! jdemod
               ! This option allows selection of either SOL22 or SOLEDGE 
               ! to be applied to a particular range (or any other solver
               ! that is implemented in the future). It isn't limited to
               ! SOL22. Solver_opt=1 specifies SOL22.
               solver_opt = sol22_regions(ir,5)   
               ! check validity of solver_opt
               if (solver_opt.lt.0.or.solver_opt.gt.1) then
                  ! set the solver_opt to sol22
                  solver_opt = 1
               endif
                  
               call plasma_solver(ixs,ixe,pzs,pze,solver_opt)

            end do

         elseif (nsol22_opt.eq.0) then 

!           jdemod            
!           This option applies SOL22 to the entire plasma 
!           space using the default SOL22 input parameters

             call set_sol22_parameter_file(sol22_default_filename)
!            call load_sol22_parameter_file('sol22_default.txt',ierr)


!     Apply SOL22 to the entire plasma region if applicable
!     SOL22 will only be applied on flux tubes with surfaces, either limiter or
!     absorbing surfaces (ie target analogue)
            ixs=1
            ixe=nxs
            pzs=1
            pze=maxpzone
            solver_opt = 1
            
            call plasma_solver(ixs,ixe,pzs,pze,solver_opt)

         endif
         
      endif

c     

      
c     Calculate Ti profiles from input vTiG profiles
c     

      
      if ((vtig_opt.eq.1.or.vtig_opt.eq.2).and.n_vtig_blocks.gt.0) then

         do pz = 1,maxpzone
            do ix = 1,ixout
               do iy =  1,nys/2
                  iqx = iqxs(ix)
                  
                  x = xouts(ix)
                  y = youts(iy)
                  yz = int(sign(1.0,youts(iy)))

                  n = crnbs(ix,iy,pz)
                  y0 = qedges(iqx,2)
                  t0 = qtembsi(iqx,2)

                  if (y.lt.y0) then
                     ctembsi(ix,iy,pz) = t0
                  else
                     call calculate_temperature(x,y-y0,pz,yz,n,t0,
     >                    cl-y0,ti,
     >                    n_vtig_blocks,vtig_range,vtig_ndata,
     >                    vtig_data,vtig_zones)
                     ctembsi(ix,iy,pz) = ti
                  endif 
               end do
               do iy =  -nys/2,-1
                  iqx = iqxs(ix)
                  
                  x = xouts(ix)
                  yz = int(sign(1.0,youts(iy)))

                  y = abs(youts(iy))
                  n = crnbs(ix,iy,pz)
                  y0 = qedges(iqx,1)
                  t0 = qtembsi(iqx,1)
                  
                  call calculate_temperature(x,y-y0,pz,yz,n,t0,
     >                 cl-y0,ti,
     >                 n_vtig_blocks,vtig_range,vtig_ndata,
     >                 vtig_data,vtig_zones)
                  ctembsi(ix,iy,pz) = ti
               end do

!     Copy the central portion to the rest of the range
               do iy = -nys,-nys/2-1
                  ctembsi(ix,iy,pz) = ctembsi(ix,iy+nys+1,pz)
               end do
               do iy = nys/2+1,nys
                  ctembsi(ix,iy,pz) = ctembsi(ix,iy-nys-1,pz)
               end do
!     average zero value between +/- 1
               ctembsi(ix,0,pz)=(ctembsi(ix,1,pz)+ctembsi(ix,-1,pz))/2.0
            end do
         end do
      endif
c     
c     Impose background vb profile 
c     
      if (vb_opt.eq.1.and.n_vb_blocks.gt.0) then 

         do pz = 1,maxpzone
            do ix = 1,ixout
               do iy =  1,nys/2
                  iqx = iqxs(ix)
                  
                  x = xouts(ix)
                  y = youts(iy)
                  yz = int(sign(1.0,youts(iy)))
                  y0 = qedges(iqx,2)
                  
                  call calculate_velocity(x,y-y0,pz,
     >                 yz,cl-y0,vel,
     >                 n_vb_blocks,vb_range,vb_ndata,
     >                 vb_data,vb_zones)
                  velplasma(ix,iy,pz) = vel
               end do
               do iy =  -nys/2,-1
                  iqx = iqxs(ix)
                  
                  x = xouts(ix)
                  yz = int(sign(1.0,youts(iy)))
                  y = abs(youts(iy))
                  y0 = qedges(iqx,1)
                  
                  call calculate_velocity(x,y-y0,pz,
     >                 yz,cl-y0,vel,
     >                 n_vb_blocks,vb_range,vb_ndata,
     >                 vb_data,vb_zones)
!     velocity profiles are entered as if for the first Y>0, Y<CL section
!     so the sign is swapped for the Y<0 Y>-CL section                    
!     the code in calculate_velocity uses the value of yz to determine the sign
!     of the returned velocity. vel = yz * vel (flow towards the target in Y>0, Y< L
!     is negative
!     velplasma(ix,iy,pz) = -vel
                  velplasma(ix,iy,pz) = vel
               end do

!     Copy the central portion to the rest of the range
               do iy = -nys,-nys/2-1
                  velplasma(ix,iy,pz) = velplasma(ix,iy+nys+1,pz) 
               end do
               do iy = nys/2+1,nys
                  velplasma(ix,iy,pz) = velplasma(ix,iy-nys-1,pz) 
               end do
!     average zero value between +/- 1
               velplasma(ix,0,pz) = (velplasma(ix,1,pz)
     >              +velplasma(ix,-1,pz))/2.0
            end do
         end do
         
      endif 
c     

      end 
      

      
      
      subroutine setup_dperp(qtim,igeom)
      use mod_params
      use mod_comt2
      use mod_comxyt
      use mod_comtor
      use yreflection
      implicit none
      real :: qtim
      integer :: igeom
      
!     locals
      integer :: j,ix,iy,iqx
      INTEGER   IQXCV1,IQXCV2,iqxbrk
      integer,external :: ipos
      REAL      RDX,NX,dnx
      integer :: pz
      pz = 1 
      
C-----------------------------------------------------------------------
C     SET DIFFUSION DECAY ETC, USING DPERP FACTORS                                 
C     ENSURE INWARD STEPPING PROBABILITIES ARE WITHIN RANGE (0,1)                  
C-----------------------------------------------------------------------
C     
      IQXBRK = CBRK
      IQXCV1 = IPOS(CVXMIN,QXS(1-NQXSO),(NQXSI+NQXSO)) -NQXSO 
      IQXCV2 = IPOS(CVXMAX,QXS(1-NQXSO),(NQXSI+NQXSO)) -NQXSO
      write (6,*) 'range for arb v:',cvxmin,cvxmax,iqxcv1,iqxcv2,cvpout
C     
      DO 130 J = 1, 3                                                           
         DO 100 IQX = 1-NQXSO, IQXBRK                                           

            IF (CVPOPT.EQ.0) THEN
               CXAFS(IQX,J) = 2.0 * CRDXO(J) * CVIN * (CA-QXS(IQX)) *
     >              QTIM * QS(IQX) / (CA*CA)                               
            ELSEIF (CVPOPT.EQ.1) THEN 
               CXAFS(IQX,J) = VPV0*(((CA-QXS(IQX))/CA)**VPALPH) *
     >              QTIM * QS(IQX)
            ELSEIF (CVPOPT.EQ.2) THEN
               IF (QXS(IQX).GE.CVPCUT) THEN 
                  CXAFS(IQX,J) = 0.0
               ELSEIF (IQX.EQ.NQXSI) THEN                
                  CXAFS(IQX,J) = 0.0
               ELSE
                  IX = IPOS(QXS(IQX),XS,NXS)
                  IF (IX.LE.NXS) THEN 
                     IY = 1
                     DNX = (CRNBS(IX,IY,pz)-CRNBS(IX-1,IY,pz))
     >                    /(XS(IX)-XS(IX-1))
                     NX = CRNBS(IX-1,IY,pz) + DNX * (QXS(IQX)-XS(IX-1))               
                     CXAFS(IQX,J) = CVIN * CRDXO(J) * DNX / NX * 
     >                    QTIM * QS(IQX)
                  ELSE
                     CXAFS(IQX,J) = 0.0
                  ENDIF 
               ENDIF
            ENDIF 
C     
C     CXAFS(IQX,J) = 2.0 * CRDXO(J) * CVIN * (CA-QXS(IQX)) *                
C     >                   QTIM * QS(IQX) / (CA*CA)                               
C     

c     sazmod
c     Will try and fit in the code for just specifying the radial 
c     diffusion into regions here, overwriting whatever happens 
c     above. Essentially all we want is to swap CRDXO(J) out with the
c     correct dperp_reg. Look at the picture in unstructured_input.f
c     if someone besides me is looking at this (hello from the past!).

c     Fortunately, we already know that if J=1 or 3 it's in the left
c     side of things, so either region 1 or 3. We just need to know
c     the X (i.e. radial) value to see if it's in the step region 
c     or not.
            
!     write(0,*) 'IQX, QXS = ', IQX, QXS(IQX)
            
            if (dperp_reg_switch.eq.1) then
               
!     Left region.
               if ((J.eq.1).or.(J.eq.3)) then
                  
!     See if in step part or not.
                  if (QXS(IQX).lt.xabsorb1a_step) then
                     
!     If in step part, region 3.
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg3*QTIM*QS(IQX))                 
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg3 
     >                    / (CDPSTP * CDPSTP) 
                  else
                     
!     If not in step part, region 1.
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg1*QTIM*QS(IQX))                 
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg1 
     >                    / (CDPSTP * CDPSTP) 
                  endif
                  
!     Right region.
               else
                  
!     If in step part, region 4.
                  if (QXS(IQX).lt.xabsorb2a_step) then
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg4*QTIM*QS(IQX))                 
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg4 
     >                    / (CDPSTP * CDPSTP) 
                     
!     If not in step part, region 2.
                  else
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg2*QTIM*QS(IQX))                 
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg2 
     >                    / (CDPSTP * CDPSTP) 
                  endif
               endif
               
            else

               CXBFS(IQX,J) = SQRT (2.0 * CRDXO(J) * QTIM * QS(IQX))                 
               CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * CRDXO(J) 
     >              / (CDPSTP * CDPSTP) 
            endif

            
 100     CONTINUE                                                                
C     
         DO 110 IQX = IQXBRK+1, NQXSI                               
            IF (CDPERP.EQ.0) THEN 
               RDX = CRDXI(J) + CRDD * QXS(IQX) / CA                                 
            ELSEIF (CDPERP.EQ.1) THEN 
               RDX = CRDXI(J)* (1.0+DPALPH*((CA-QXS(IQX))/CA)**DPBETA)
            ENDIF
            IF (CVPOPT.EQ.0) THEN
               CXAFS(IQX,J) = 2.0 * RDX * CVIN * (CA-QXS(IQX)) *                     
     >              QTIM * QS(IQX) / (CA*CA)                               
            ELSEIF (CVPOPT.EQ.1) THEN 
               CXAFS(IQX,J) = VPV0*(((CA-QXS(IQX))/CA)**VPALPH) *
     >              QTIM * QS(IQX)
            ELSEIF (CVPOPT.EQ.2) THEN
               IF (QXS(IQX).GE.CVPCUT) THEN 
                  CXAFS(IQX,J) = 0.0
               ELSEIF (IQX.EQ.NQXSI) THEN                
                  CXAFS(IQX,J) = 0.0
               ELSE
                  IX = IPOS(QXS(IQX),XS,NXS)
                  IF (IX.LE.NXS) THEN 
                     IY = 1
                     DNX = (CRNBS(IX,IY,pz)-CRNBS(IX-1,IY,pz))
     >                    /(XS(IX)-XS(IX-1))
                     NX = CRNBS(IX-1,IY,pz) + DNX * (QXS(IQX)-XS(IX-1))                     
                     CXAFS(IQX,J) = CVIN * RDX * DNX / NX*QTIM*QS(IQX)
                  ELSE
                     CXAFS(IQX,J) = 0.0
                  ENDIF 
               ENDIF
            ENDIF 
            
!     sazmod
!     A little different for the inboard side compared to the
!     outboard, but still simple enough. This time we are just
!     replacing RDX. 
!     write(0,*) 'IQX, QXS = ', IQX, QXS(IQX)
            if (dperp_reg_switch.eq.1) then
               
!     Left region.
               if ((J.eq.1).or.(J.eq.3)) then
                  
!     If in step part, region 3.
                  if (QXS(IQX).lt.xabsorb2a_step) then
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg3*QTIM*QS(IQX))                      
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg3 
     >                    / (CDPSTP * CDPSTP) 
                     
!     If not in step part, region 1.
                  else
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg1*QTIM*QS(IQX))                      
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg1 
     >                    / (CDPSTP * CDPSTP) 
                  endif
                  
!     Right region.
               else  
                  
!     If in step part, region 4.
                  if (QXS(IQX).lt.xabsorb2a_step) then
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg4*QTIM*QS(IQX))                      
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg4 
     >                    / (CDPSTP * CDPSTP) 
                     
!     If not in step part, region 2.
                  else
                     CXBFS(IQX,J) = SQRT (2.0*dperp_reg2*QTIM*QS(IQX))                      
                     CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg2 
     >                    / (CDPSTP * CDPSTP) 
                  endif
                  
               endif
               
            else
               CXBFS(IQX,J) = SQRT (2.0 * RDX * QTIM * QS(IQX))                      
               CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * RDX 
     >              / (CDPSTP * CDPSTP) 
            endif
            
 110     CONTINUE                                                                
C     
C     An outward pinch velocity or an arbitrary pinch over a specified 
C     region may also be entered. This is simply added to the contents
C     of CXAFS - the array containing the step due to the pinch velocity
C     component of the motion.
C     
C     Note: The definition of signs in the code differs from that used 
C     in the literature. Normally, an inward pinnch is descibed by
C     Vpinch = -2.0 Dperp * r / a**2 - in the code the "+" ve sign
C     is used for motion towards the core - the "-" ve sign is used 
C     for drift towards the SOL. Thus an arbitrary pinch of +0.5 m/s 
C     in the specification list -  is a drift towards the walls and
C     is applied to the CXAFS array as   "-" Vin * QTIM * QS(IQX) 
C     
         DO 115 IQX = IQXCV1,IQXCV2
            CXAFS(IQX,J) = CXAFS(IQX,J) + CVPOUT * QTIM * QS(IQX)
 115     CONTINUE
C     
         DO 120 IQX = 1-NQXSO, NQXSI                                             
            IF     (IGEOM.EQ.0) THEN                                              
               CXCFS(IQX,J) = 0.5                                                  
            ELSEIF (IGEOM.EQ.1) THEN                                              
               CXCFS(IQX,J) = (CA-QXS(IQX)-0.5*CXBFS(IQX,J)) /                     
     >              (2.0*(CA-QXS(IQX)))                                  
               CXCFS(IQX,J) = MIN (1.0, MAX (0.0, CXCFS(IQX,J)))                   
            ENDIF                                                                 
 120     CONTINUE                                                                
 130  CONTINUE                                                                  


      return
      end



      subroutine wrt_plasma
      use mod_params
      use mod_global_options
      use mod_comt2
      use mod_comxyt
      use mod_io_units
      implicit none
      integer :: pz
      ! write out the plasma solution zone by zone to stdout

      if (cprint.ne.9) return
      
      do pz = 1,maxpzone
         call prt_bg_array(pz,crnbs,'Density')
         call prt_bg_array(pz,ctembs,'Te')
         call prt_bg_array(pz,ctembsi,'Ti')
         call prt_bg_array(pz,velplasma,'Vb')
         call prt_bg_array(pz,efield,'Efield')
      end do

      return
      end
      
      subroutine prt_bg_array(pz,bgarray,name)
      use mod_params
      use mod_comt2
      use mod_comxyt
      use mod_io_units
      use yreflection
      use debug_options
      implicit none

      integer :: pz,ix,iy
      real:: bgarray(maxnxs,-maxnys:maxnys,maxpzone)
      character*(*) :: name
      
      write(stdout,*) 'Plasma '//trim(name)//' for zone:',pz

      write(stdout,'(5a13,1000(1x,g12.5))') '   X   ','   YABS1   ',
     >     '   YABS2   ','   YABS1_EXT  ','   YABS2_EXT  ',    
     >     (youts(iy),iy=-nys,nys)
      
      do ix = 1,nxs

         ! try to write it on one line
         write(stdout,'(1000(1x,g12.5))') xouts(ix),
     >        yabsorb_surf(ix,pz,1),yabsorb_surf(ix,pz,2),
     >        yabsorb_surf_ext(ix,pz,1),yabsorb_surf_ext(ix,pz,2),
     >        (bgarray(ix,iy,pz),iy=-nys,nys)

      end do
      
      return
      end
