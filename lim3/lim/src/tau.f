      SUBROUTINE TAUIN1 (QTIM,NIZS,ICUT,                                        
     >                   FSRATE,IGEOM,NTBS,NTIBS,NNBS)                         
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
      IMPLICIT  none
      REAL      QTIM,FSRATE                                                     
      INTEGER   NIZS,ICUT(2),IGEOM,NTBS,NTIBS,NNBS,IQXBRK
C                                                                               
C***********************************************************************        
C                                                                               
C       SETS UP VALUES IN COMMON BLOCKS COMTAU/COMT2                            
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c
c      include   'global_options'
c
c slmod
c     INCLUDE   'cadas'
c      INCLUDE   'slcom'
c slmod end
C                                              ,ICNT
c slmdo begin
      INTEGER   IPOS,IXOUT,IZ,IQX,LIMIZ,IX,IY,J,ip                                
      !INTEGER   IQXCV1,IQXCV2 
      REAL      FEX,WIDTH,FEXZ                              
c slmod
      CHARACTER MESAGE*80
      REAL      TMPION      
      
      
c      PARAMETER (LAMBDA=10.06)
c      PARAMETER (LAMBDA=15.0)

c       IF (CIOPTE.EQ.10) THEN

c         OPEN (UNIT=51,FILE='../data/limin/tmp.dli',STATUS='OLD')
 
c         WRITE (0,*) ' '
c         WRITE (0,*) '--------------------------------------------'

c         READ (51,'(A33,F5.2)')    MESAGE,LAMBDA

c Just calculate LAMBDA here ya idiot.
c lambda = 17.3 - .5*ln(n/1e20) + 1.5*ln(T/1000.)
         LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)
c         WRITE(0,*) 'Calculating LAMBDA ...', LAMBDA

c         WRITE(0 ,'(1X,A33,1X,F4.1)') MESAGE,LAMBDA

c         READ (51,'(A33,G8.2)')    MESAGE,TMPION
c         WRITE(0 ,'(1X,A33,2X,G9.3)') MESAGE,TMPION

c         WRITE(0,*) '--------------------------------------------'
c         WRITE(0,*) ' '

c         CLOSE (UNIT=51)

c       ELSE
c         LAMBDA = 15.0
c       ENDIF


c        WRITE(0,*) 'Starting TAU'

c slmod end
C                                                                               
         LIMIZ = MIN (CION, NIZS)                                                  
c         write(0,*) 'LIMIZ:', cion,nizs,limiz
         write(0,"(A,I2,A,I2,A,I2)") ' CION = ',cion, ' NIZS = ',nizs, 
     >      ' --> LIMZ = ',limiz
C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP QEDGES, QTANS AND QDISTS                           
c
c     jdemod - setup geometry BEFORE calculating plasma since some of
c     the plasma calculations may depend on the location of the
c     limiter edges. 
C-----------------------------------------------------------------------
C                                                                               
      WRITE (6,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      WRITE (0,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      CALL EDGE (QXS,QEDGES,QTANS,QDISTS,NQXSO,CAW,CL,ICUT,CCUT,XSCALO,         
     >           WEDGAN,XL1,YL1,XL2,YL2,TC,SC,TO,SO,GC,RP,CIOPTH,CORECT,
     >           XST1,YST1,XST2,YST2,XST3,YST3,RLEDGE7,CA,RLC)        

C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP CTEMBS AND CRNBS, QTEMBS AND QRNBS                 
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (6,'('' TAU: CALLING PLASMA OPTIONS'',I3,I3)') 
     >           CIOPTG,CIOPTK                   
      WRITE (0,'('' TAU: CALLING PLASMA OPTIONS'',I3,I3)') 
     >           CIOPTG,CIOPTK                   

      CALL PLASMA (NTBS,NTIBS,NNBS,CIOPTG,CIOPTK,QTIM)                     
C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP VARIABLE CAW
C-----------------------------------------------------------------------        
c
      call setup_wall (qys,nqys,cl,caw)

c     sazmod - Load in 2D data for fully customizable absorbing boundary.
c     The data gets stored in the array "bounds". Only for yabsorb1a 
c     side of simulation currently.
      if (vary_2d_bound.eq.1) then
        write(0,*) 'Loading in varying boundary...'
        write(6,*) 'Loading in varying boundary...'
        call load_varying_boundary_1a
      endif

C     
C-----------------------------------------------------------------------        
C                     SET UP CEYS AND CVHYS                                     
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (6,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      WRITE (0,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      CALL SOL (QYS,CEYS,CVHYS,NQYS,CTBIN,CTIBIN,CRMB,CL,CIZB,                 
     >          CEYOUT,CVHOUT,CYSTAG,CRMI,CSOLEF,CIOPTF)                        

c     Whatever has been done with these old LIM options to the plasma,
c     go through and assign it to the 3D/4D arrays just to stick to
c     the convention where the real solution is overlaid on top of the
c     initial LIM plasma. Fun little Fortran tidbit: order of loops
c     matters for speed since Fortran stores arrays as column-major 
c     (look it up).
      do iy = -nys, nys
        do ix = 1, nxs
          do ip = 1, npbins
            ctembs_3d(ip, ix, iy) = ctembs(ix, iy)
            ctembsi_3d(ip, ix, iy) = ctembsi(ix, iy)
            crnbs_3d(ip, ix, iy) = crnbs(ix, iy)
          end do
        end do
      end do

c
c     jdemod - at this point the LIM plasma has been fully calculated
c              using the base options which are simple and relatively efficient
c     - both SOLEDGE and SOL22 are then implemented as overlays on top of
c       this existing plasma description. velplasma and efield can be initialized
c       using this so that it contains valid values for inboard and other regions 
c       where SOL22 is not appropriate
c
      call plasma_overlay(qtim)
c
c     jdemod - after plasma has been finalized - calculate temperature gradients
c      
      
      write(0,*) 'crmi,qtim:',crmi,qtim
      call calculate_tgrad(qtim)
c
c     Depending on the plasma overlay option specified - rewrite the ion temperature to
c     be constant at the target value.       
c
c     These seems to be new and currently maybe still being implemented,
c     so will not worry about vtig_mod for vary_2d_bound (unless I see
c     vtig_opt = 0, the default value).
      if (vtig_opt.eq.2) then 
         do ix = 1,nxs
            iqx = iqxs(ix)
            do iy = -nys/2,-1
               ctembsi(ix,iy) = qtembsi(iqx,1)
               ctembsi(ix,iy+nys+1) =  qtembsi(iqx,1)
            end do   
            do iy = -nys,-nys/2-1
               ctembsi(ix,iy) = qtembsi(iqx,2)
               ctembsi(ix,iy+nys+1) =  qtembsi(iqx,2)
            end do   
         end do
      endif

         
      ! jdemod - write out background plasma
      if (cprint.eq.9) then 
         write(6,*) 'PLASMA BACKGROUND AFTER OVERLAY:'      
         do ix = 1,nxs
            do iy = -nys,nys
               write(6,'(2i8,10(1x,g12.5))') ix,iy,xouts(ix),youts(iy),
     >           ctembs(ix,iy),
     >           ctembsi(ix,iy),crnbs(ix,iy),ctigs(ix,iy),ctegs(ix,iy),
     >           velplasma(ix,iy,1)
            end do
         end do
      endif

      
C
c     jdemod - pull out code to calculate varying dperp - must be run
c     after plasma density is finalized       
c 
      call setup_dperp(qtim,igeom)

c
c     jdemod: Scale the outboard flow velocity outside the limiters
c             (only used in 3D) 
c             (QX(IQX) scaling is done in line when the new location
c              is calculated in lim3.f)
c
      vpflow_3d = vpflow_3d * qtim

C                                                                               
C-----------------------------------------------------------------------        
C              SET UP CYSCLS, CFEXZS, CFVHXS                                    
C              ALL VALUES DEFINED FOR OUTBOARD X POSITIONS ONLY                 
c
c              These need to be calculated inboard for colprobe3d on
c     
C-----------------------------------------------------------------------
C                                                                               
      DO 200  IQX = 1-NQXSO, 0                                                  
         WIDTH      = CTWOL - QEDGES(IQX,1) - QEDGES(IQX,2)                     
         CYSCLS(IQX)= REAL(NQYS) / WIDTH                                        
  200 CONTINUE                                                                  
C                                                                               
      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
C                                                                               
      FEX = QTIM * QTIM * (1.602192E-19 / (1.672614E-27 * CRMI))                
      IF (LIMIZ.GT.0) THEN                                                      
        DO 300  IZ = 1, LIMIZ                                                   
          FEXZ = FEX * REAL (IZ)                                                
          DO 250 IY = -NYS, NYS     
                                                      
           ! Changed to NXS to support transport forces inboard of the 
           ! probe tip                                               
           DO 250 IX = 1, NXS
            IQX = IQXS(IX)                                                      
            
            ! Varying 2D boundary uses different arrays.
            if (vary_2d_bound.eq.1) then
              do ip = 1, npbins
                if (vel_efield_opt.eq.0) then
                  if (ix.gt.ixout) then 
                    cfexzs_4d(ip,ix,iy,iz) = fexz * ctembs_3d(ip,ix,iy)
     >                 / ctbin * qs(iqx) * qs(iqx)           
                   else
                     cfexzs_4d(ip,ix,iy,iz) = fexz * ctembs_3d(ip,ix,iy)
     >                  / ctbin * cyscls(iqx)/yscale * qs(iqx) * qs(iqx)           
                   endif
                elseif (vel_efield_opt.eq.1) then 
                
                   ! If using velplasma/efield values then the CFVHXS 
                   ! contains only timestep scaling and not temperature 
                   ! relative to the separatrix
                   if (ix.gt.ixout) then 
                      cfexzs_4d(ip,ix,iy,iz) = fexz * qs(iqx) * qs(iqx)           
                   else
                   
                      ! Not sure about the cyscls/yscale factor for 
                      ! efield - leave for now
                      cfexzs_4d(ip,ix,iy,iz) = fexz *                      
     >                         cyscls(iqx)/yscale * qs(iqx) * qs(iqx)           
                   endif
                endif
              end do
              
            ! Normal stuff.
            else
            
                if (vel_efield_opt.eq.0) then
                   if (ix.gt.ixout) then 
                      CFEXZS(IX,IY,IZ) = FEXZ * CTEMBS(IX,IY)/CTBIN 
     >                             * QS(IQX) * QS(IQX)           
                   else
                      CFEXZS(IX,IY,IZ) = FEXZ * CTEMBS(IX,IY)/CTBIN *                     
     >                         CYSCLS(IQX)/YSCALE * QS(IQX) * QS(IQX)           
                   endif
                elseif (vel_efield_opt.eq.1) then 
                   !  if using velplasma/efield values then the CFVHXS contains only timestep
                   ! scaling and not temperature relative to the separatrix
                   if (ix.gt.ixout) then 
                      CFEXZS(IX,IY,IZ) = FEXZ * QS(IQX) * QS(IQX)           
                   else
                      ! not sure about the cyscls/yscale factor for efield - leave for now
                      CFEXZS(IX,IY,IZ) = FEXZ *                      
     >                         CYSCLS(IQX)/YSCALE * QS(IQX) * QS(IQX)           
                   endif
                endif
            
            endif
               
 250        CONTINUE                                                              
  300   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      !do 310 iy = -nys, nys  
      do iy = -nys, nys                     
                                     
        ! Changed to NXS to support transport forces inboard of the 
        ! probe tip.
        !do 310 ix = 1, nxs 
        do ix = 1, nxs                                                     
          iqx = iqxs(ix)          
          if (vary_2d_bound.eq.1) then
            do ip = 1, npbins
              if (vel_efield_opt.eq.0) then 
                cfvhxs_3d(ip, ix,iy) = sqrt((ctembs_3d(ip, ix,iy) + 
     >            ctembsi_3d(ip,ix,iy))/(ctbin+ctibin)) * qtim * qs(iqx) 
                 
              elseif (vel_efield_opt.eq.1) then
            
                ! If using velplasma/efield values then the cfvhxs 
                ! contains only timestep scaling and not temperature 
                ! relative to the separatrix
                cfvhxs_3d(ip,ix,iy) = qtim * qs(iqx)
              endif
            end do
        
          ! Normal stuff.
          else                                            
            if (vel_efield_opt.eq.0) then 
              cfvhxs(ix,iy) = sqrt((ctembs(ix,iy) + ctembsi(ix,iy)) / 
     >          (ctbin+ctibin)) * qtim * qs(iqx)             
            elseif (vel_efield_opt.eq.1) then
            
              ! If using velplasma/efield values then the cfvhxs 
              ! contains only timestep scaling and not temperature 
              ! relative to the separatrix
              cfvhxs(ix,iy) = qtim * qs(iqx)             
            endif 
          endif
        end do
      end do
c 310  continue                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP CMIZS                                              
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
C    SET IONISATION / E-I RECOMBINATION TIME INTERVALS    CFIZS,CFRCS           
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (0,'('' TAU: CALLING IZTAU  OPTION'',I3)') CIOPTA 
      WRITE (6,'('' TAU: CALLING IZTAU  OPTION'',I3)') CIOPTA                     
c      CALL IZTAU (CRMI,NXS,NYS,CION,CIZB,CIOPTA)   
      CALL IZTAU (CIOPTA)                             
C                                                                               
C-----------------------------------------------------------------------        
C    SET COMBINED C-X AND E-I RECOMBINATION TIMES         CFCXS                 
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (0,'('' TAU: CALLING CXREC  OPTION'',I3)') CIOPTI 
      WRITE (6,'('' TAU: CALLING CXREC  OPTION'',I3)') CIOPTI                  
      CALL CXREC (NIZS,CION,CIOPTI,CIZB,CL,CRMB,CVCX,                           
     >            CNHC,CNHO,CLAMHX,CLAMHY)          
     
c     I've never used charge-exchange in my simulations before, and I'm
c     not sure if it's even relevant really, so unless otherwise needed
c     I'm just gonna create an array of zeros.
      if (vary_2d_bound.eq.1) then
        do iy = -nys, nys
          do ix = 1, nxs
            do ip = 1, npbins
              do iz = 0, cion
                cfcxs_4d(ip,ix,iy,iz) = 0.0
              end do
            end do
          end do
        end do
      endif     
                             
C                                                                               
C-----------------------------------------------------------------------        
C     SET PROBABILITY OF EITHER AN IONISATION OR A RECOMBINATION                
C     SET PROPORTION OF THESE WHICH WILL BE RECOMBINATIONS                      
C     PREVENT ANY IONISATION BEYOND MAXIMUM LIMIT SPECIFIED IF REQUIRED         
C-----------------------------------------------------------------------        
C                                                                              
      if (cmizs.gt.1) then                                                      
        do 370 iz = 1, cmizs-1                                                  
         do 360 ix = 1, nxs                                                     
          iqx = iqxs(ix)                                                        
          do 350 iy = -nys, nys 
          
            ! 2D boundary routine.    
            if (vary_2d_bound.eq.1) then
              do ip = 1, npbins
                if (cfcxs_4d(ip,ix,iy,iz).le.0.0) then                                    
                  cpchs_4d(ip,ix,iy,iz) = qtim * qs(iqx) / 
     >              cfizs_4d(ip,ix,iy,iz)                
                  cprcs_4d(ip,ix,iy,iz) = 0.0                                             
                else                                                                
                  cpchs_4d(ip,ix,iy,iz) = (cfizs_4d(ip,ix,iy,iz) + 
     >               cfcxs_4d(ip,ix,iy,iz)) * qtim * qs(iqx) / 
     >              (cfizs_4d(ip,ix,iy,iz) * cfcxs_4d(ip,ix,iy,iz))             
     
                  ! jdemod - cfizs and cfcxs contain characteristic 
                  ! TIMES (see print outs below) fast ionization is a 
                  ! SMALLER time meaning less chance for recombination. 
                  cprcs_4d(ip,ix,iy,iz) = cfizs_4d(ip,ix,iy,iz) /                               
     >              (cfcxs_4d(ip,ix,iy,iz) + cfizs_4d(ip,ix,iy,iz))             
                endif                                                               
                cpchs_4d(ip,ix,iy,iz) = min(1.0, cpchs_4d(ip,ix,iy,iz)) 
              end do
              
            ! Normal routine.
            else                                            
              if (cfcxs(ix,iy,iz).le.0.0) then                                    
                cpchs(ix,iy,iz) = qtim * qs(iqx) / cfizs(ix,iy,iz)                
                cprcs(ix,iy,iz) = 0.0                                             
              else                                                                
                cpchs(ix,iy,iz) = (cfizs(ix,iy,iz) + cfcxs(ix,iy,iz)) *           
     >            qtim * qs(iqx) /(cfizs(ix,iy,iz) * cfcxs(ix,iy,iz))             

                ! jdemod - cfizs and cfcxs contain characteristic TIMES 
                ! (see print outs below) fast ionization is a SMALLER 
                ! time meaning less chance for recombination. 
                cprcs(ix,iy,iz) = cfizs(ix,iy,iz) /                               
     >            (cfcxs(ix,iy,iz) + cfizs(ix,iy,iz))             
              endif                                                               
              cpchs(ix,iy,iz) = min (1.0, cpchs(ix,iy,iz)) 
            endif                       
  350     continue                                                              
  360    continue                                                               
  370   continue                                                                
      endif                                                                     
                                                                               
      if (cmizs .le. limiz) then                                                
        do 390 ix = 1, nxs                                                      
          iqx = iqxs(ix)                                                        
          do 380 iy = -nys, nys 
          
            ! 2D boundary routine.
            if (vary_2d_bound.eq.1) then
              do ip = 1, npbins
                if (cfcxs_4d(ip,ix,iy,iz).le.0.0) then                                    
                  cpchs_4d(ip,ix,iy,cmizs) = 0.0                                          
                else                                                                
                  cpchs_4d(ip,ix,iy,cmizs) = qtim * qs(iqx) / 
     >              cfcxs_4d(ip,ix,iy,iz)             
                endif                                                               
                cprcs_4d(ip,ix,iy,cmizs) = 1.0                                            
                cpchs_4d(ip,ix,iy,iz) = min(1.0, cpchs_4d(ip,ix,iy,iz))  
              end do
              
            ! Normal routine.
            else                                                
              if (cfcxs(ix,iy,iz).le.0.0) then                                    
                cpchs(ix,iy,cmizs) = 0.0                                          
              else                                                                
                cpchs(ix,iy,cmizs) = qtim * qs(iqx) / cfcxs(ix,iy,iz)             
              endif                                                               
              cprcs(ix,iy,cmizs) = 1.0                                            
              cpchs(ix,iy,iz) = min (1.0, cpchs(ix,iy,iz))  
            endif                      
  380     continue                                                              
  390   continue                                                                
      endif                                                                     
C                                                                               
C---- SET IONISATION PROBABILITIES FOR NEUT ...                                 
C---- (SAVES REPEATED CALCULATION EVERY ITERATION)                              
C---- THEY ARE ALL MULTIPLIED BY THE "IONISATION RATE FACTOR" IRF               
C---- WHICH TYPICALLY MIGHT BE 0.2 TO GIVE DEEPER IONISATION.                   
C                                                                               
      do 500 iy = -nys, nys                                                     
        do 500 ix = 1, nxs      
          if (vary_2d_bound.eq.1) then
            do ip = 1, npbins
              cpchs_4d(ip,ix,iy,0) = min(1.0, cirf * fsrate / 
     >          cfizs_4d(ip,ix,iy,0))
            end do
          else
            cpchs(ix,iy,0) = min (1.0, cirf * fsrate / cfizs(ix,iy,0)) 
          endif
  500 continue                                                                  
C                                                                               
c slmod begin - N2 break
      IF (N2OPT.EQ.1) THEN
        WRITE(63,*) ' '
        WRITE(63,*) 'Ionisation probabilities:'
        WRITE(63,*) ' '
        WRITE(63,'(A3,3A12)')
     +    'IX','NNCPCHS','N2CPCHS','CFIZS0'

        DO IX = 1, NXS
          N2CPCHS(IX) = MIN (1.0, CIRF * FSRATE / N2RATE(IX))            
          NNCPCHS(IX) = MIN (1.0, CIRF * FSRATE / NNRATE(IX))            

          WRITE(63,'(I3,3E12.5)') 
     +      IX,NNCPCHS(IX),N2CPCHS(IX),CPCHS(IX,1,0)
        ENDDO
        WRITE(63,*) ' '
      ENDIF
c slmod end

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
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
C                                                                               
      INTEGER   IPOS,IZ,IQX,LIMIZ,JX,IX,IY,ip                                      
!      real      tmp1
      
      ! test for issues with precision - especially when TAU values are multiplied
      REAL*8    TEMP,FTAU,FTAUP,FTAUS,FTAUT,RIZSQR,STAU,TAU                     
      real*8    lambda
      REAL*8    ROOTMI,ROOTTT                                            

c     slmod
c      PARAMETER (LAMBDA=15.0)                                                   
C                                                                               
c       IF (CIOPTE.EQ.10) THEN
c lambda = 17.3 - .5*ln(n/1e20) + 1.5*ln(T/1000.)
         LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)
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
      ftau  = czenh * sqrt(crmb) * cizb * cizb * lambda * qtim * sf_tau                 
      ftaup = ftau * 6.8e-14                                                    
      ftaus = ftau * 6.8e-14 * (1.0 + crmb/crmi)                                
      ftaut = ftau * 1.4e-13                                                    
      c215a = ftau * 1.4e-13 * sqrt(crmi)                                       
      c215b = (crmi * ctbi + crmb * ctemsc) ** 1.5                              
                                                                               
      do 540  iz = 1, limiz                                                     
         rizsqr = real (iz) * real (iz)                                         
         do 520 iy = -nys, nys                                                  
          do 520  ix = 1, nxs     
          
            ! 2D varying boundary routine. Comments removed since this 
            ! code is copied from the normal routine (scroll down), just
            ! the arrays are swapped with the 3D/4D ones that include ip.
            if (vary_2d_bound.eq.1) then
              do ip = 1, npbins
                if (ctbi.gt.0.0) then                                               
                  temp  = crnbs_3d(ip,ix,iy) / (crmi * ctbi**1.5)                          
                else                                                                
                  temp  = crnbs_3d(ip,ix,iy) / (crmi * 
     >              ctembsi_3d(ip,ix,iy)**1.5)                
                endif                                                               
                iqx = iqxs(ix)                                                      
                stau = temp * rizsqr * qs(iqx)                                      
                ctolds_3d(ip,ix,iz) = ctemsc                                              
                       
                ! Tau parallel                                                                
                if (ctbi.gt.0.0) then                                               
                  cfps_4d(ip,ix,iy,iz) = stau * ctbi * ftaup * 2.0                        
                else                                                                
                  cfps_4d(ip,ix,iy,iz) = stau * ctembsi_3d(ip,ix,iy) * 
     >              ftaup * 2.0             
                endif                                                               
                                                                               
                if (cioptb.eq.1.and.ix.le.jx) then                              
                  cfps_4d(ip,ix,iy,iz) = 0.0                                             
                elseif (cioptb.eq.2.and.ix.le.jx) then                              
                  cfps_4d(ip,ix,iy,iz) = 2.0 * crnbs_3d(ip,ix,iy) 
     >              * 6.8e-14 * real (cizb * cizeff) * rizsqr * lambda /         
     >              (rootmi * roottt) * qtim * qs(iqx)               
                endif                                                               
            
                ! Innermost loop something.                                                                              
                if (cfps_4d(ip,ix,iy,iz).eq.0.0) then                                     
                   cccfps_4d(ip,ix,iy,iz) = 0.0                                           
                elseif (cioptb.eq.3 .and. ix.le.jx) then 
                   cccfps_4d(ip,ix,iy,iz) = sqrt (9.76e8 * ctemsc / 
     >                crmi) * qtim * qs(iqx) / cfps_4d(ip,ix,iy,iz)                                
                elseif (cioptb.eq.4 .and. ix.le.jx) then
                   cfps_4d(ip,ix,iy,iz) = 2.0e0 * cfps_4d(ip,ix,iy,iz)
                   cccfps_4d(ip,ix,iy,iz) = sqrt(9.76e8 * ctemsc / crmi) 
     >                * qtim * qs(iqx) / cfps_4d(ip,ix,iy,iz)

                ! Something Lisgo did with minimal comments, as 
                ! is tradition.
                elseif (cioptb.eq.13) then
                   cccfps_4d(ip,ix,iy,iz) = 1.56e4 * sqrt(pi/4.0 * 
     >               1.0/crmi * (cfps_4d(ip,ix,iy,iz) * (1.0+crmb/crmi))
     >               / 2.0) * qtim
                else                                                                
                   cccfps_4d(ip,ix,iy,iz) = sqrt (4.88e8 / 
     >               (cfps_4d(ip,ix,iy,iz) * crmi)) * qtim * qs(iqx)                                                 
                endif  
                                                                             
                ! Scaling factor to the parallel diffusive transport.
                cccfps_4d(ip,ix,iy,iz) = cccfps_4d(ip,ix,iy,iz)*sf_vdiff
            
                ! Tau stopping.                                                                                                                                  
                tau = stau * ftaus                                                  
                if (tau.gt.1.e-3) then                                              
                  cfss_4d(ip,ix,iy,iz) = 1.0 - exp(-tau)                                  
                else                                                                
                  cfss_4d(ip,ix,iy,iz) = tau                                              
                endif                                                               
            
                if (cioptc.eq.1.and.ix.le.jx) then                              
                  cfss_4d(ip,ix,iy,iz) = 0.0                                             
                elseif (cioptc.eq.2.and.ix.le.jx) then                              
                  tau = cfps_4d(ip,ix,iy,iz) / (2.0 * ctemsc)                              
                  if (tau.gt.1.e-3) then                                           
                    cfss_4d(ip,ix,iy,iz) = 1.0 - exp(-tau)                               
                  else                                                             
                    cfss_4d(ip,ix,iy,iz) = tau                                           
                  endif                                                            
                elseif (cioptc.eq.3 .and. ix.le.jx) then
                  cfss_4d(ip,ix,iy,iz) = 1.0e20 
                endif                                                               

                ! Tau heating.                                                                     
                tau = stau * ftaut                                                  
                if (tau.gt.1.e-3) then                                              
                  cfts_4d(ip,ix,iy,iz) = 1.0 - exp(-tau)                                  
                else                                                                
                  cfts_4d(ip,ix,iy,iz) = tau                                              
                endif                                                               
                                                                               
                if (cioptd.eq.1.and.ix.le.jx) then                              
                  cfts_4d(ip,ix,iy,iz) = 0.0                                             
                elseif (cioptd.eq.2.and.ix.le.jx) then                              
                  cfts_4d(ip,ix,iy,iz) = 1.0                                             
                elseif (cioptd.eq.3.and.ix.le.jx) then                              
                  if (ctbi.le.0.0) then                                                 
                    c215b = (crmi * ctembsi_3d(ip,ix,iy) + crmb * 
     >                ctemsc) ** 1.5   
                  endif   
                  tau = c215a * rizsqr * qs(iqx) * crnbs_3d(ip,ix,iy) 
     >              / c215b            
                  if (tau.gt.1.e-3) then                                           
                    cfts_4d(ip,ix,iy,iz) = 1.0 - exp(-tau)                               
                  else                                                             
                    cfts_4d(ip,ix,iy,iz) = tau                                           
                  endif                                                            
                endif
              end do
              
            ! Normal routine, cleaned up some.
            else
              if (ctbi.gt.0.0) then                                               
                temp  = crnbs(ix,iy) / (crmi * ctbi**1.5)                          
              else                                                                
                temp  = crnbs(ix,iy) / (crmi * ctembsi(ix,iy)**1.5)                
              endif                                                               
              iqx = iqxs(ix)                                                      
              stau = temp * rizsqr * qs(iqx)                                      
              ctolds(ix,iz) = ctemsc                                              
       
              ! Tau parallel - Notes 3, 50, 103                                                                                      
              ! Standard case: cfps = 2 * deltat * ti / tauparallel.                                     
              ! For non-zero options, special case applies for X outside 
              ! of cplsma only. Each time an ion enters this region its               
              ! temperature is compared against a previous value to see             
              ! whether any values need recalculating.                              
              if (ctbi.gt.0.0) then                                               
                cfps(ix,iy,iz) = stau * ctbi * ftaup * 2.0                        
              else                                                                
                cfps(ix,iy,iz) = stau * ctembsi(ix,iy) * ftaup * 2.0             
              endif                                                               
                                                                               
              if (cioptb.eq.1.and.ix.le.jx) then                              
                cfps(ix,iy,iz) = 0.0                                             
              elseif (cioptb.eq.2.and.ix.le.jx) then                              
                cfps(ix,iy,iz) = 2.0 * crnbs(ix,iy) * 6.8e-14 *                  
     >            real (cizb * cizeff) * rizsqr * lambda /         
     >            (rootmi * roottt) * qtim * qs(iqx)               
              endif                                                               
                                                                              
              ! Set additional array for use in innermost loop of LIM2              
              ! Equivalent to some constant times cfps values.                      
              ! This saves 20% of cpu time by eliminating a square root             
              ! cccfps = sqrt (4.88e8/(cfps.mi)) . deltat. sin(thetab)              
              ! for note 284, set cccfps = as above . sqrt(2ti/cfps)                
                                                                             
              if (cfps(ix,iy,iz).eq.0.0) then                                     
                cccfps(ix,iy,iz) = 0.0                                           
              elseif (cioptb.eq.3 .and. ix.le.jx) then 
                cccfps(ix,iy,iz) = sqrt (9.76e8 * ctemsc / crmi) *               
     >          qtim * qs(iqx) / cfps(ix,iy,iz)                                
              elseif (cioptb.eq.4 .and. ix.le.jx) then
                cfps(ix,iy,iz) = 2.0e0 * cfps(ix,iy,iz)
                cccfps(ix,iy,iz) = sqrt(9.76e8 * ctemsc / crmi ) *
     >            qtim * qs(iqx) / cfps(ix,iy,iz)
     
              ! slmod
              elseif (cioptb.eq.13) then

                ! crmb    - plasma ion mass
                ! crmi    - impurity ion mass
                ! cizb    - plasma ion charge
                ! ctembsi - local background temperature
                ! crnbs   - local background density
                cccfps(ix,iy,iz) = 1.56e4 * sqrt(pi/4.0 * 1.0/crmi
     >            * (cfps(ix,iy,iz) *(1.0+crmb/crmi))  /2.0) * qtim

              else                                                                
                cccfps(ix,iy,iz) = sqrt (4.88e8 /(cfps(ix,iy,iz)*crmi))*         
     >            qtim * qs(iqx)                                                 
              endif
                                                                             
              ! Apply the scaling factor to the parallel diffusive 
              ! transport. sf_vdiff default value is 1.0. This quantity 
              ! is used to scale either spatial or velocity diffusive
              ! step sizes depending on which is in use.
              cccfps(ix,iy,iz) = cccfps(ix,iy,iz) * sf_vdiff
                   
              ! Tau stopping - Notes 3, 50, 103                               
              ! Standard case: cfss = 1 - exp (-deltat/taustopping)                                
              ! Non zero options apply outside of X = cplsma only ...               
              ! For option 2, taustopping = tauparallel                             
              ! = 2.deltat.ti/cfps                                                                           
              tau = stau * ftaus                                                  
              if (tau.gt.1.e-3) then                                              
                cfss(ix,iy,iz) = 1.0 - exp(-tau)                                  
              else                                                                
                cfss(ix,iy,iz) = tau                                              
              endif                                                               
            
              if (cioptc.eq.1.and.ix.le.jx) then                              
                cfss(ix,iy,iz) = 0.0                                             
              elseif (cioptc.eq.2.and.ix.le.jx) then                              
                tau = cfps(ix,iy,iz) / (2.0*ctemsc)                              
                if (tau.gt.1.e-3) then                                           
                  cfss(ix,iy,iz) = 1.0 - exp(-tau)                               
                else                                                             
                  cfss(ix,iy,iz) = tau                                           
                endif                                                            
              elseif (cioptc.eq.3 .and. ix.le.jx) then
                cfss(ix,iy,iz) = 1.0e20 
              endif                                                               
       
              ! Tau heating: Notes 3, 50, 103, 215                           
              ! Standard case: cfts = 1 - exp (-deltat/tauheating)                                 
              ! Non zero options apply outside of x = cplsma only ...                                                                        
              tau = stau * ftaut                                                  
              if (tau.gt.1.e-3) then                                              
                cfts(ix,iy,iz) = 1.0 - exp(-tau)                                  
              else                                                                
                cfts(ix,iy,iz) = tau                                              
              endif                                                               
                                                                               
              if (cioptd.eq.1.and.ix.le.jx) then                              
                cfts(ix,iy,iz) = 0.0                                             
              elseif (cioptd.eq.2.and.ix.le.jx) then                              
                cfts(ix,iy,iz) = 1.0                                             
              elseif (cioptd.eq.3.and.ix.le.jx) then                              
                if (ctbi.le.0.0)                                                 
     >            c215b = (crmi * ctembsi(ix,iy)+ crmb*ctemsc) ** 1.5      
                tau = c215a * rizsqr * qs(iqx) * crnbs(ix,iy) / c215b            
                if (tau.gt.1.e-3) then                                           
                  cfts(ix,iy,iz) = 1.0 - exp(-tau)                               
                else                                                             
                  cfts(ix,iy,iz) = tau                                           
                endif                                                            
              endif
            
            endif
520      continue                                                               
540   continue                                                                  
                
      !write(6,*) 'CFSS,CFTS:'
      !do ix = 1,nxs
      !   do iy = -nys,nys
      !      write(6,'(2i8,10(1x,g12.5))') 
     >!        ix,iy,cfss(ix,iy,iz),cfts(ix,iy,iz) 
      !   end do
      !end do

      RETURN                                                                    
      END                                                                       
C                                                                               
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
C       PRINTS REPRESENTATIVE VALUES FROM THE ARRAYS OF                         
C     FACTORS SET UP IN COMMONS COMTAU & COMT2                                  
C                    CHRIS FARRELL    DECEMBER 1987                             
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c      INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c
c      include   'global_options'
C                                                                               
      INTEGER   IQX,IZ,PMIZS,MAXIZ,J                                            
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
     >TIMESTEP  ')                                                              
      ELSE                                                                      
      CALL PRC ('FUNCTIONS    LIMITER    PLASMA TEMPS      DENSITIES            
     >TIMESTEP  ')                                                              
      ENDIF                                                                     
      CALL PRC ('  OF X      Y<0   Y>0   Y<0     Y>0      Y<0      Y>0          
     >  FACTOR  ')                                                              

      WRITE (7,9107) 'X = A   ',                                                
     >                   CTEMBS(IXA  ,-1),CTEMBS(IXA  ,1),                      
     >             CRNBS (IXA  ,-1),CRNBS (IXA  ,1),QS(IQXA  )                  
      WRITE (7,9107) 'X = A/2 ',                                                
     >                   CTEMBS(IXA2 ,-1),CTEMBS(IXA2 ,1),                      
     >             CRNBS (IXA2 ,-1),CRNBS (IXA2 ,1),QS(IQXA2 )                  
      WRITE (7,9107) 'X = A/4 ',                                                
     >                   CTEMBS(IXA4 ,-1),CTEMBS(IXA4 ,1),                      
     >             CRNBS (IXA4 ,-1),CRNBS (IXA4 ,1),QS(IQXA4 )                  
      WRITE (7,9107) 'INBOARD ',                                                
     >                   CTEMBS(IXIN ,-1),CTEMBS(IXIN ,1),                      
     >             CRNBS (IXIN ,-1),CRNBS (IXIN ,1),QS(IQXIN )                  
      WRITE (7,9108) 'OUTBOARD',-QEDGES(IQXOUT,1),                              
     >  QEDGES(IQXOUT,2),CTEMBS(IXOUT,-1),CTEMBS(IXOUT,1),                      
     >             CRNBS (IXOUT,-1),CRNBS (IXOUT,1),QS(IQXOUT)                  
      WRITE (7,9108) 'X =-AW/4',-QEDGES(IQXAW4,1),                              
     >  QEDGES(IQXAW4,2),CTEMBS(IXAW4,-1),CTEMBS(IXAW4,1),                      
     >             CRNBS (IXAW4,-1),CRNBS (IXAW4,1),QS(IQXAW4)                  
      WRITE (7,9108) 'X =-AW/2',-QEDGES(IQXAW2,1),                              
     >  QEDGES(IQXAW2,2),CTEMBS(IXAW2,-1),CTEMBS(IXAW2,1),                      
     >             CRNBS (IXAW2,-1),CRNBS (IXAW2,1),QS(IQXAW2)                  
      WRITE (7,9108) 'X =-AW  ',-QEDGES(IQXAW ,1),                              
     >  QEDGES(IQXAW ,2),CTEMBS(IXAW ,-1),CTEMBS(IXAW ,1),                      
     >             CRNBS (IXAW ,-1),CRNBS (IXAW ,1),QS(IQXAW )                  
      IF (IQXFAC.LT.0)                                                          
     >  WRITE (7,9108) 'X =-LAM ',-QEDGES(IQXFAC,1),                            
     >  QEDGES(IQXFAC,2),CTEMBS(IXFAC,-1),CTEMBS(IXFAC,1),                      
     >             CRNBS (IXFAC,-1),CRNBS (IXFAC,1),QS(IQXFAC)                  
C                                                                               
      IF (CTBI.GT.0.0) THEN                                                     
        CALL PRB                                                                
        CALL PRR ('*** BACKGROUND ION TEMP OVERRIDDEN TO VALUE TBI =',          
     >    CTBI)                                                                 
      ELSE
C
C     IF THE BACKGROUND ION TEMPERATURE IS NOT OVERRIDDEN BY   
C     A CONSTATNT, PRINT OUT REPRESENTATIVE VALUES. 
C
      CALL PRC ('FUNCTIONS    LIMITER    PLASMA ION TEMPS                       
     >TIMESTEP  ')                                                              
      CALL PRC ('  OF X      Y<0   Y>0   Y<0     Y>0                            
     >  FACTOR  ')                                                              
      WRITE (7,9109) 'X = A   ',                                                
     >     CTEMBSI(IXA  ,-1),CTEMBSI(IXA  ,1),QS(IQXA  )                  
      WRITE (7,9109) 'X = A/2 ',                                                
     >     CTEMBSI(IXA2 ,-1),CTEMBSI(IXA2 ,1), QS(IQXA2 )                  
      WRITE (7,9109) 'X = A/4 ',                                                
     >     CTEMBSI(IXA4 ,-1),CTEMBSI(IXA4 ,1),QS(IQXA4 )                  
      WRITE (7,9109) 'INBOARD ',                                                
     >     CTEMBSI(IXIN ,-1),CTEMBSI(IXIN ,1),QS(IQXIN )                  
      WRITE (7,9110) 'OUTBOARD',-QEDGES(IQXOUT,1),                              
     >  QEDGES(IQXOUT,2),CTEMBSI(IXOUT,-1),CTEMBSI(IXOUT,1),                   
     >  QS(IQXOUT)                  
      WRITE (7,9110) 'X =-AW/4',-QEDGES(IQXAW4,1),                              
     >  QEDGES(IQXAW4,2),CTEMBSI(IXAW4,-1),CTEMBSI(IXAW4,1),                   
     >  QS(IQXAW4)                  
      WRITE (7,9110) 'X =-AW/2',-QEDGES(IQXAW2,1),                              
     >  QEDGES(IQXAW2,2),CTEMBSI(IXAW2,-1),CTEMBSI(IXAW2,1),                    
     >  QS(IQXAW2)                  
      WRITE (7,9110) 'X =-AW  ',-QEDGES(IQXAW ,1),                              
     >  QEDGES(IQXAW ,2),CTEMBSI(IXAW ,-1),CTEMBSI(IXAW ,1),                   
     >  QS(IQXAW )                  
      IF (IQXFAC.LT.0)                                                          
     >  WRITE (7,9110) 'X =-LAM ',-QEDGES(IQXFAC,1),                            
     >  QEDGES(IQXFAC,2),CTEMBSI(IXFAC,-1),CTEMBSI(IXFAC,1),                    
     >  QS(IQXFAC)                  
C                                                                               
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (CPRINT.EQ.1.or.cprint.eq.9) THEN
      DO 222 J = 1, 3
      CALL PRB                                                                  
      IF (J.EQ.1) THEN                                                          
      CALL PRC ('X DRIFT & INWARD PINCH (M), & INWARD STEP PROBS'//
     >          ' -L/2 < Y < 0')          
      ELSEIF (J.EQ.2) THEN                                                
      CALL PRC ('X DRIFT & INWARD PINCH (M), & INWARD STEP PROBS'// 
     >          '  0 < Y < L/2')          
      ELSE                                                                      
      CALL PRC ('X DRIFT & INWARD PINCH (M), & INWARD STEP PROBS'//
     >          ' L/2< AY <3L/2')          
      ENDIF                                                                     
      WRITE (7,9003) 'NEAR X = A     ', CXBFS(IQXA  ,J) /                       
     >  SQRT(QS(IQXA  )),CXAFS(IQXA  ,J)/QS(IQXA  ),CXCFS(IQXA  ,J)             
      WRITE (7,9003) 'NEAR X = A/2   ', CXBFS(IQXA2 ,J) /                       
     >  SQRT(QS(IQXA2 )),CXAFS(IQXA2 ,J)/QS(IQXA2 ),CXCFS(IQXA2 ,J)             
      WRITE (7,9003) 'NEAR X = A/4   ', CXBFS(IQXA4 ,J) /                       
     >  SQRT(QS(IQXA4 )),CXAFS(IQXA4 ,J)/QS(IQXA4 ),CXCFS(IQXA4 ,J)             
      WRITE (7,9003) 'JUST INBOARD   ', CXBFS(IQXIN ,J) /                       
     >  SQRT(QS(IQXIN )),CXAFS(IQXIN ,J)/QS(IQXIN ),CXCFS(IQXIN ,J)             
      WRITE (7,9003) 'JUST OUTBOARD  ', CXBFS(IQXOUT,J) /                       
     >  SQRT(QS(IQXOUT)),CXAFS(IQXOUT,J)/QS(IQXOUT),CXCFS(IQXOUT,J)             
      WRITE (7,9003) 'NEAR X =-AW/4  ', CXBFS(IQXAW4,J) /                       
     >  SQRT(QS(IQXAW4)),CXAFS(IQXAW4,J)/QS(IQXAW4),CXCFS(IQXAW4,J)             
      WRITE (7,9003) 'NEAR X =-AW/2  ', CXBFS(IQXAW2,J) /                       
     >  SQRT(QS(IQXAW2)),CXAFS(IQXAW2,J)/QS(IQXAW2),CXCFS(IQXAW2,J)             
      WRITE (7,9003) 'NEAR X =-AW    ', CXBFS(IQXAW ,J) /                       
     >  SQRT(QS(IQXAW )),CXAFS(IQXAW ,J)/QS(IQXAW ),CXCFS(IQXAW ,J)             
      IF (IQXFAC.LT.0)                                                          
     >WRITE (7,9003) 'NEAR X =-LAMBDA', CXBFS(IQXFAC,J) /                       
     >  SQRT(QS(IQXFAC)),CXAFS(IQXFAC,J)/QS(IQXFAC),CXCFS(IQXFAC,J)             
  222 CONTINUE                                                                  
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      IF (MAXIZ.LT.1) GOTO 999                                                  
C-----------------------------------------------------------------------        
      IF (CIOPTI.NE.0)                                                          
     >  CALL TAUPRA (7,CNHS,'NEUTRAL HYDROGEN DENSITY (M**-3)',-1)              
C-----------------------------------------------------------------------        
      IF (CPRINT.EQ.1.or.cprint.eq.9) THEN
      CALL PRB                                                                  
      CALL PRC ('IONISATION AND RECOMBINATION TIMES')                           
      WRITE (7,9101)                                                            
      CALL PRB                                                                  
      CALL PRC ('  VALUES NEAR X = A/2  (Y>0)')                                        
      WRITE (7,9102) 0,CFIZS(IXA2,IY0,0)                                        
      DO IZ = 1, PMIZS                                                       
       WRITE (7,9102) IZ,CFIZS(IXA2 ,IY0,IZ),CFRCS(IXA2 ,IY0,IZ),               
     >   CFCXS(IXA2 ,IY0,IZ),CFCXS(IXA2 ,IYL8,IZ),CFCXS(IXA2 ,IYL4,IZ)          
      end do
c      call prb
c      CALL PRC ('  VALUES NEAR X = A/2  (Y<0)')                                        
c      WRITE (7,9102) 0,CFIZS(IXA2,IY0LT,0)                                        
c      DO IZ = 1, PMIZS                                                       
c       WRITE (7,9102) IZ,CFIZS(IXA2 ,IY0LT,IZ),CFRCS(IXA2 ,IY0LT,IZ),               
c     >   CFCXS(IXA2,IY0LT,IZ),CFCXS(IXA2,-IYL8,IZ),CFCXS(IXA2 ,-IYL4,IZ)          
c      end do
C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC ('  VALUES JUST OUTBOARD (Y>0)')                                       
      WRITE (7,9102) 0,CFIZS(IXOUT,IY0,0)                                       
      DO IZ = 1, PMIZS                                                       
       WRITE (7,9102) IZ,CFIZS(IXOUT,IY0,IZ),CFRCS(IXOUT,IY0,IZ),               
     >   CFCXS(IXOUT,IY0,IZ),CFCXS(IXOUT,IYL8,IZ),CFCXS(IXOUT,IYL4,IZ)          
      end do
      CALL PRB                                                                  
      CALL PRC ('  VALUES JUST OUTBOARD (Y<0)')                                       
      WRITE (7,9102) 0,CFIZS(IXOUT,IY0LT,0)                                       
      DO IZ = 1, PMIZS                                                       
       WRITE (7,9102) IZ,CFIZS(IXOUT,IY0LT,IZ),CFRCS(IXOUT,IY0LT,IZ),               
     > CFCXS(IXOUT,IY0LT,IZ),CFCXS(IXOUT,-IYL8,IZ),CFCXS(IXOUT,-IYL4,IZ)          
      end do
C-----------------------------------------------------------------------        
      IF (IQXFAC.LT.0) THEN                                                     
        CALL PRB                                                                
        CALL PRR ('  VALUES NEAR X = -LAMBDA (Y>0) = ',QXS(IQXFAC))                   
        WRITE (7,9102) 0,CFIZS(IXFAC,IY0,0)                                     
        DO IZ = 1, PMIZS                                                     
          WRITE (7,9102) IZ,CFIZS(IXFAC,IY0,IZ),CFRCS(IXFAC,IY0,IZ),            
     >     CFCXS(IXFAC,IY0,IZ),CFCXS(IXFAC,IYL8,IZ),CFCXS(IXFAC,IYL4,IZ)        
        end do
        CALL PRB                                                                
        CALL PRR ('  VALUES NEAR X = -LAMBDA (Y<0) = ',QXS(IQXFAC))                   
        WRITE (7,9102) 0,CFIZS(IXFAC,IY0LT,0)                                     
        DO IZ = 1, PMIZS                                                     
          WRITE (7,9102) IZ,CFIZS(IXFAC,IY0LT,IZ),CFRCS(IXFAC,IY0LT,IZ),            
     > CFCXS(IXFAC,IY0LT,IZ),CFCXS(IXFAC,-IYL8,IZ),CFCXS(IXFAC,-IYL4,IZ)        
        end do
      ENDIF                                                                     
C-----------------------------------------------------------------------        
C                                                                               
C NOTE ON THIS SECTION WHEN VARIABLE DELTAT IS IN USE:-  IF A PROB. IS          
C SAY 0.22 FOR THE STANDARD DT, THEN WHEN A TIMESTEP MULTIPLE OF SAY            
C 1000 IS USED THIS SHOOTS UP TO "220", WHICH IS ADJUSTED TO 1 IN               
C TAUIN1.  THEN WHEN WE ATTEMPT TO RECREATE THE GENUINE PROB. FOR PRINT         
C PURPOSES BELOW, IT BECOMES 0.001 WHICH IS OBVIOUSLY WRONG.  HOWEVER,          
C WHAT ELSE CAN ONE DO?  AT LEAST THE CALCULATION USES CORRECT VALUES.          
C                                                                               
      CALL PRB                                                                  
      CALL PRC ('CHANGE OF STATE PROBABILITIES (USING LOCAL TIMESTEPS FO        
     >R IONS)')                                                                 
      WRITE (7,9103)                                                            
      CALL PRB                                                                  
      CALL PRC ('  VALUES NEAR X = A/2 (Y>0)')                                        
      DO IZ = 0, PMIZS                                                       
        WRITE (7,9104) IZ,                                                      
     >    CPCHS(IXA2 ,IY0 ,IZ)          ,100.0*CPRCS(IXA2 ,IY0 ,IZ),            
     >    CPCHS(IXA2 ,IYL8,IZ)          ,100.0*CPRCS(IXA2 ,IYL8,IZ),            
     >    CPCHS(IXA2 ,IYL4,IZ)          ,100.0*CPRCS(IXA2 ,IYL4,IZ)             
      end do
c      CALL PRB                                                                  
c      CALL PRC ('  VALUES NEAR X = A/2 (Y<0)')                                        
c      DO IZ = 0, PMIZS                                                       
c        WRITE (7,9104) IZ,                                                      
c     >    CPCHS(IXA2 ,IY0LT ,IZ)         ,100.0*CPRCS(IXA2 ,IY0LT,IZ),            
c     >    CPCHS(IXA2 ,-IYL8,IZ)          ,100.0*CPRCS(IXA2 ,-IYL8,IZ),            
c     >    CPCHS(IXA2 ,-IYL4,IZ)          ,100.0*CPRCS(IXA2 ,-IYL4,IZ)             
c      end do
C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC ('  VALUES JUST OUTBOARD (Y>0)')                                       
      DO IZ = 0, PMIZS                                                       
        WRITE (7,9104) IZ,                                                      
     >    CPCHS(IXOUT,IY0 ,IZ)           ,100.0*CPRCS(IXOUT,IY0 ,IZ),           
     >    CPCHS(IXOUT,IYL8,IZ)           ,100.0*CPRCS(IXOUT,IYL8,IZ),           
     >    CPCHS(IXOUT,IYL4,IZ)           ,100.0*CPRCS(IXOUT,IYL4,IZ)            
      end do
      CALL PRB                                                                  
      CALL PRC ('  VALUES JUST OUTBOARD (Y<0)')                                       
      DO IZ = 0, PMIZS                                                       
        WRITE (7,9104) IZ,                                                      
     >    CPCHS(IXOUT,IY0LT ,IZ)         ,100.0*CPRCS(IXOUT,IY0LT,IZ),           
     >    CPCHS(IXOUT,-IYL8,IZ)          ,100.0*CPRCS(IXOUT,-IYL8,IZ),           
     >    CPCHS(IXOUT,-IYL4,IZ)          ,100.0*CPRCS(IXOUT,-IYL4,IZ)            
      end do
C-----------------------------------------------------------------------        
      IF (IQXFAC.LT.0) THEN                                                     
        CALL PRB                                                                
        CALL PRR ('  VALUES NEAR X = -LAMBDA (Y>0) = ',QXS(IQXFAC))                   
        DO IZ = 0, PMIZS                                                     
          WRITE (7,9104) IZ,                                                    
     >      CPCHS(IXFAC,IY0 ,IZ)           ,100.0*CPRCS(IXFAC,IY0 ,IZ),         
     >      CPCHS(IXFAC,IYL8,IZ)           ,100.0*CPRCS(IXFAC,IYL8,IZ),         
     >      CPCHS(IXFAC,IYL4,IZ)           ,100.0*CPRCS(IXFAC,IYL4,IZ)          
        end do
        CALL PRB                                                                
        CALL PRR ('  VALUES NEAR X = -LAMBDA (Y<0) = ',QXS(IQXFAC))                   
        DO IZ = 0, PMIZS                                                     
          WRITE (7,9104) IZ,                                                    
     >      CPCHS(IXFAC,IY0LT,IZ)         ,100.0*CPRCS(IXFAC,IY0LT,IZ),         
     >      CPCHS(IXFAC,-IYL8,IZ)         ,100.0*CPRCS(IXFAC,-IYL8,IZ),         
     >      CPCHS(IXFAC,-IYL4,IZ)         ,100.0*CPRCS(IXFAC,-IYL4,IZ)          
        end do
      ENDIF                                                                     


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
      CALL PRB                                                                  
      CALL PRC ('X POSITION FACTORS ')                                          
      CALL PRI ('ELECTRIC FIELD, IONIZATION STATE ', 1)                         
      CALL PRR (' JUST OUTBOARD   ',CFEXZS(IXOUT,1,1)/QS(IQXOUT)**2)            
      CALL PRR (' NEAR X =-AW/4   ',CFEXZS(IXAW4,1,1)/QS(IQXAW4)**2)            
      CALL PRR (' NEAR X =-AW/2   ',CFEXZS(IXAW2,1,1)/QS(IQXAW2)**2)            
      CALL PRR (' NEAR X =-AW     ',CFEXZS(IXAW ,1,1)/QS(IQXAW )**2)            
      IF (IQXFAC.LT.0)                                                          
     >CALL PRR (' NEAR X =-LAMBDA ',CFEXZS(IXFAC,1,1)/QS(IQXFAC)**2)            
      IF (MAXIZ.GT.1) THEN                                                      
      CALL PRI ('ELECTRIC FIELD, IONIZATION STATE', MAXIZ)                      
      CALL PRR (' JUST OUTBOARD   ',CFEXZS(IXOUT,1,MAXIZ)/QS(IQXOUT)**2)        
      CALL PRR (' NEAR X =-AW/4   ',CFEXZS(IXAW4,1,MAXIZ)/QS(IQXAW4)**2)        
      CALL PRR (' NEAR X =-AW/2   ',CFEXZS(IXAW2,1,MAXIZ)/QS(IQXAW2)**2)        
      CALL PRR (' NEAR X =-AW     ',CFEXZS(IXAW ,1,MAXIZ)/QS(IQXAW )**2)        
      IF (IQXFAC.LT.0)                                                          
     >CALL PRR (' NEAR X =-LAMBDA ',CFEXZS(IXFAC,1,MAXIZ)/QS(IQXFAC)**2)        
      ENDIF                                                                     
C-----------------------------------------------------------------------        
      CALL PRB                                                                  
      CALL PRC ('DRIFT VELOCITY, ALL IONIZATION STATES')                        
      CALL PRR (' JUST OUTBOARD   ', CFVHXS(IXOUT,1)/QS(IQXOUT))                
      CALL PRR (' NEAR X =-AW/4   ', CFVHXS(IXAW4,1)/QS(IQXAW4))                
      CALL PRR (' NEAR X =-AW/2   ', CFVHXS(IXAW2,1)/QS(IQXAW2))                
      CALL PRR (' NEAR X =-AW     ', CFVHXS(IXAW ,1)/QS(IQXAW ))                
      IF (IQXFAC.LT.0)                                                          
     >CALL PRR (' NEAR X =-LAMBDA ', CFVHXS(IXFAC,1)/QS(IQXFAC))                
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
      CALL TAUPRA (6,CTEMBS, 'PLASMA TEMPERATURE (EV)   ',-1)              
      CALL TAUPRA (6,CTEMBSI,'PLASMA ION TEMPS (EV)     ',-1) 
      IF (NTIG.NE.0.OR.NTEG.NE.0) THEN
        CALL TAUPRA (7,CTEMBS, 'PLASMA TEMPERATURE (EV)   ',-1)              
        CALL TAUPRA (7,CTEMBSI,'PLASMA ION TEMPS (EV)     ',-1) 
      ENDIF
      CALL TAUPRA (6,CRNBS ,'PLASMA DENSITY     (M**-3)',-1)                    
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
C-505 CONTINUE                                                                  
C--   DO 510 IZ = 1, NIZS                                                       
C--   CALL TAUPRA (6,CFEXZS(1,-MAXNYS,IZ),'ELECTRIC FIELD FACTOR',IZ)           
C-510 CONTINUE                                                                  
C--   CALL TAUPRA (6,CFVHXS              ,'DRIFT VELOCITY FACTOR',-1)           
C-----------------------------------------------------------------------        
 9003 FORMAT(1X,'  ',A15,1P,' +/-',G11.2,' +',G11.2,',  PROB',0P,F8.5)          
 9101 FORMAT( 1X,'       IONISATION  E-I RECOMB  MODIFIED RECOM',               
     >                               'BINATION TIMES DUE           ',           
     >       /1X,'        TIME FOR    TIME FOR   TO CHARGE EXCH',               
     >                               'ANGE RECOMBINATION           ',           
     >       /1X,'  IZ   IZ -> IZ+1  IZ -> IZ-1  NEAR Y = 0  NE',               
     >                               'AR Y=L/8 NEAR Y=L/4          ')           
 9102 FORMAT(1X,I4,1P,5(1X,G11.2))                                              
 9103 FORMAT( 1X,'      PROBABILITY  RECOMB PROBABILITY  RECOMB ',              
     >                                  'PROBABILITY RECOMB       ',            
     >       /1X,'  IZ  NEAR Y = 0    FRAC  NEAR Y=L/8    FRAC  ',              
     >                                  'NEAR Y=L/4   FRAC        ')            
 9104 FORMAT(1X,I4,3(1X,1P,G11.2,0P,F7.2,'%'))                                  
 9106 FORMAT(2X,A,1P,5(1X,G11.2))                                               
 9107 FORMAT(2X,A,12X,2F8.1,1X,1P,2G9.2,0P,F8.1)                                
 9108 FORMAT(2X,A,1X,2F6.3,F7.1,F8.1,1X,1P,2G9.2,0P,F8.1)                       
 9109 FORMAT(2X,A,12X,2F8.1,1X,1P,18X,0P,F8.1)             
 9110 FORMAT(2X,A,1X,2F6.3,F7.1,F8.1,1X,1P,18X,0P,F8.1)                       
 9200 FORMAT(/1X,A,1P,/,(1X,14G9.2))                                            
C                                                                               
  999 RETURN                                                                    
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
C       PRINTS REPRESENTATIVE VALUES OF CHARACTERISTIC TIMES, DERIVED           
C     FROM FACTORS CFPS,CFSS,CFTS SET UP IN COMMON COMT2, USING THE             
C     INITIAL ION TEMPERATURE RETURNED FROM NEUT/LAUNCH.                        
C                                                                               
C     NOTE "TAUPRF" IS INCLUDED IN MEMBER "RUNLIM3" SINCE IT USES F77           
C     CHARACTER FEATURES INCOMPATIBLE WITH THE CRAY CFT2 COMPILER.              
C     HENCE THE REST OF TAU CAN BE COMPILED WITH CFT2.                          
C                                                                               
C                         CHRIS FARRELL    DECEMBER 1987                        
C                                                                               
C***********************************************************************        
C                                                                               
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c      INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c
c      include   'global_options'
C                                                                               
      INTEGER   MAXIZ,IZ                                                        
C-----------------------------------------------------------------------        
      MAXIZ  = MIN (NIZS, CION)                                                 
      IF (CPRINT.EQ.0) GOTO 999                                                 
      CALL PRB                                                                  
      CALL PRC ('SAMPLE CHARACTERISTIC TIMES BASED ON INITIAL ION TEMPER        
     >ATURE')                                                                   
      CALL PRC ('                        TAU PARALLEL  TAU STOPPING   TA        
     >U HEATING')                                                               
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
     >CALL TAUPRF ('NEAR X =-LAMBDA',IXFAC,1,1,CTEMSC,QTIM,CIOPTC)              
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
     >CALL TAUPRF ('NEAR X =-LAMBDA',IXFAC,1,MAXIZ,CTEMSC,QTIM,CIOPTC)          
      ENDIF                                                                     
C                                                                               
C     DO 500 IZ = 1, MAXIZ                                                      
C       CALL TAUPRA (6,CFPS  (1,-MAXNYS,IZ),'CFPS FACTORS',IZ)                  
C       CALL TAUPRA (6,CFSS  (1,-MAXNYS,IZ),'CFSS FACTORS',IZ)                  
C       CALL TAUPRA (6,CFTS  (1,-MAXNYS,IZ),'CFTS FACTORS',IZ)                  
C       CALL TAUPRA (6,CCCFPS(1,-MAXNYS,IZ),'CCCFPS FACTORS',IZ)                
C 500 CONTINUE                                                                  
C                                                                               
  999 RETURN                                                                    
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
      REAL RIZSQR,RATIO1,RATIO2,TAU                                             
      INTEGER IQX,IY                                                            
C                                                                               
C     WRITE (6,9001) TEMOLD,                                                    
C    >  CFPS(IX,1,CIZ),CCCFPS(IX,1,CIZ),CFSS(IX,1,CIZ),CFTS(IX,1,CIZ)           
C                                                                               
      IF (CIOPTB.EQ.2) THEN                                                     
        RATIO1 = SQRT (TEMOLD / TEMNEW)                                         
        RATIO2 = 1.0 / SQRT (RATIO1)                                            
        DO 100 IY = -NYS, NYS                                                   
          CFPS(IX,IY,CIZ)   = RATIO1 * CFPS(IX,IY,CIZ)                          
          CCCFPS(IX,IY,CIZ) = RATIO2 * CCCFPS(IX,IY,CIZ)                        
  100   CONTINUE                                                                
C                                                                               
      ELSEIF (CIOPTB.EQ.3) THEN                                                 
        RATIO2 = SQRT (TEMNEW / TEMOLD)                                         
        DO 110 IY = -NYS, NYS                                                   
          CCCFPS(IX,IY,CIZ) = RATIO2 * CCCFPS(IX,IY,CIZ)                        
  110   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      IF (CIOPTC.EQ.2) THEN                                                     
        DO 200 IY = -NYS, NYS                                                   
          TAU = CFPS(IX,IY,CIZ) / (2.0*TEMNEW)                                  
          IF (TAU.GT.1.E-3) THEN                                                
            CFSS(IX,IY,CIZ) = 1.0 - EXP(-TAU)                                   
          ELSE                                                                  
            CFSS(IX,IY,CIZ) = TAU                                               
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
     >      C215B = (CRMI * CTEMBSI(IX,IY) + CRMB * TEMNEW)** 1.5             
          TAU = C215A * RIZSQR * QS(IQX) * CRNBS(IX,IY) / C215B                 
          IF (TAU.GT.1.E-3) THEN                                                
            CFTS(IX,IY,CIZ) = 1.0 - EXP (-TAU)                                  
          ELSE                                                                  
            CFTS(IX,IY,CIZ) = TAU                                               
          ENDIF                                                                 
  300   CONTINUE                                                                
      ENDIF                                                                     
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
C@PROCESS OPT(0),VECTOR(LEV(0))                                                 
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
      IMPLICIT  none
      REAL      QTIM                                                            
      INTEGER   NIZS                                                            
      DOUBLE PRECISION DEFACT                                                   
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  TAUIN3:  THIS ROUTINE RECALCULATES THE CHARACTERISTIC TIMES      *        
C  *  ACCORDING TO THE STEADY STATE ION DISTRIBUTIONS AT THE END OF    *        
C  *  A LIM RUN.  THESE TIMES CAN THEN BE USED FOR ANOTHER ITERATION   *        
C  *  OF LIM, UNTIL EVENTUALLY WE MAY ACHIEVE A SELF CONSISTENT PLASMA *        
C  *  AS DESCRIBED IN NOTE 207.                                        *        
C  *  NOTE THAT TO PREVENT DIVIDE BY ZEROES IN INNER LOOPS THIS        *        
C  *  ROUTINE HAS TO BE COMPILED WITH VECTORISATION OFF.               *        
C  *  NOTE 297: ZEFFS RELATED QUANTITIES ARE NOW MULTIPLIED BY         *        
C  *  SIN(THETAB) AND THE RATIO LP+/LP-                                *        
C  *                                                                   *        
C  *                          CHRIS FARRELL  (HUNTERSKIL)  SEPT 1988   *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
c     INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c     INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
c      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c      INCLUDE   'dynam1'                                                        
C     INCLUDE   (DYNAM1)                                                        
c      INCLUDE   'dynam3'                                                        
C     INCLUDE   (DYNAM3)                                                        
C                                                                               
      INTEGER   IX,IY,IZ,IQX,JZ                                                 
      REAL      LAMBDA,TEMP,SQRTMB,SQRTMI,TOTALP,ARG,RIZB,TOTALT                
      REAL      TAUP(0:MAXIZS),TAUS(0:MAXIZS),TAUT(0:MAXIZS),TOTALS,RIZ         
c slmod
c      PARAMETER (LAMBDA=15.0)                                                   
C                                                                               
c       IF (CIOPTE.EQ.10) THEN
c lambda = 17.3 - .5*ln(n/1e20) + 1.5*ln(T/1000.)
         LAMBDA = 17.3 - 0.5*LOG(CNBIN/1.0E20) + 1.5*LOG(CTIBIN/1000.0)
c       ELSE
c         LAMBDA = 15.0
c       ENDIF
c slmod end
C                                                                               
      SQRTMB = SQRT (CRMB)                                                      
      SQRTMI = SQRT (CRMI)                                                      
      RIZB   = REAL (CIZB)                                                      
C                                                                               
      DO 530 IZ = 1, NIZS                                                       
       WRITE (6,9001) IZ,(JZ,JZ=1,9)                                            
       RIZ = REAL (IZ)                                                          
       DO 520 IY = -NYS, NYS                                                    
        IF (IY.EQ.0) GOTO 520                                                   
        DO 510 IX = 1, NXS                                                      
         IF (CTBI.GT.0.0) THEN                                                  
           TEMP = CTBI                                                          
         ELSE                                                                   
           TEMP = CTEMBSI(IX,IY)                                             
         ENDIF                                                                  
         IQX = IQXS(IX)                                                         
C_______________________________________________________________________        
C                                                                               
C        TAU PARALLEL     CFPS = 2.DELTAT.TI/TAUPARA                            
C_______________________________________________________________________        
C                                                                               
         TOTALP = 0.0                                                           
         IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
           TAUP(0) = CRMI * SQRT(TEMP) /                                        
     >      (6.8E-14 * SQRTMB * ZEFFS(IX,IY,5) *                                
     >       RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
         ELSE                                                                   
           TAUP(0) = 0.0                                                        
         ENDIF                                                                  
         IF (TAUP(0).GT.0.0) TOTALP = 1.0 / TAUP(0)                             
C                                                                               
         DO 100 JZ = 1, NIZS                                                    
           IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
             TAUP(JZ) = CRMI * SQRT(SNGL(DDTS(IX,IY,JZ))) /                     
     >        (6.8E-14 * SQRTMI * SNGL(DDLIMS(IX,IY,JZ)*DEFACT) *               
     >         CSINTB * REAL(JZ)**2 * RIZ * RIZ * LAMBDA * CZENH)               
           ELSE                                                                 
             TAUP(JZ) = 0.0                                                     
           ENDIF                                                                
           IF (TAUP(JZ).GT.0.0) TOTALP = TOTALP + 1.0 / TAUP(JZ)                
  100    CONTINUE                                                               
C                                                                               
         CFPS(IX,IY,IZ) = 2.0 * QTIM * QS(IQX) * TOTALP                         
C                                                                               
         IF (CFPS(IX,IY,IZ).LE.0.0) THEN                                        
           CCCFPS(IX,IY,IZ) = 0.0                                               
         ELSE                                                                   
           CCCFPS(IX,IY,IZ) = SQRT (4.88E8 /(CFPS(IX,IY,IZ)*CRMI))*             
     >                        QTIM * QS(IQX)                                    
         ENDIF                                                                  
C_______________________________________________________________________        
C                                                                               
C        TAU STOPPING     CFSS = 1 - EXP (-DELTAT/TAUSTOP)                      
C_______________________________________________________________________        
C                                                                               
         TOTALS = 0.0                                                           
         IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
           TAUS(0) = CRMI * (TEMP**1.5) /                                       
     >      (6.8E-14 * SQRTMB * (1.0+CRMB/CRMI) * ZEFFS(IX,IY,5) *              
     >       RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
         ELSE                                                                   
           TAUS(0) = 0.0                                                        
         ENDIF                                                                  
         IF (TAUS(0).GT.0.0) TOTALS = 1.0 / TAUS(0)                             
C                                                                               
         DO 200 JZ = 1, NIZS                                                    
           IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
             TAUS(JZ) = CRMI * (SNGL(DDTS(IX,IY,JZ))**1.5) /                    
     >        (6.8E-14 * SQRTMI * 2.0 * SNGL(DDLIMS(IX,IY,JZ)*DEFACT) *         
     >         CSINTB * REAL(JZ)**2 * RIZ * RIZ * LAMBDA * CZENH)               
           ELSE                                                                 
             TAUS(JZ) = 0.0                                                     
           ENDIF                                                                
           IF (TAUS(JZ).GT.0.0) TOTALS = TOTALS + 1.0 / TAUS(JZ)                
  200    CONTINUE                                                               
C                                                                               
         ARG = QTIM * QS(IQX) * TOTALS                                          
         IF (ARG.GT.1.E-3) THEN                                                 
           CFSS(IX,IY,IZ) = 1.0 - EXP(-ARG)                                     
         ELSE                                                                   
           CFSS(IX,IY,IZ) = ARG                                                 
         ENDIF                                                                  
C_______________________________________________________________________        
C                                                                               
C        TAU HEATING      CFTS = 1 - EXP (-DELTAT/TAUHEAT)                      
C_______________________________________________________________________        
C                                                                               
         TOTALT = 0.0                                                           
         IF (ZEFFS(IX,IY,5).GT.0.0) THEN                                        
           TAUT(0) = CRMI * (TEMP**1.5) /                                       
     >      (1.4E-13 * SQRTMB * ZEFFS(IX,IY,5) *                                
     >       RIZB * RIZ * RIZ * LAMBDA * CZENH)                                 
         ELSE                                                                   
           TAUT(0) = 0.0                                                        
         ENDIF                                                                  
         IF (TAUT(0).GT.0.0) TOTALT = 1.0 / TAUT(0)                             
C                                                                               
         DO 300 JZ = 1, NIZS                                                    
           IF (DDLIMS(IX,IY,JZ).GT.0.0D0) THEN                                  
             TAUT(JZ) = CRMI * (SNGL(DDTS(IX,IY,JZ))**1.5) /                    
     >        (1.4E-13 * SQRTMI * SNGL(DDLIMS(IX,IY,JZ)*DEFACT) *               
     >         CSINTB * REAL(JZ)**2 * RIZ * RIZ * LAMBDA * CZENH)               
           ELSE                                                                 
             TAUT(JZ) = 0.0                                                     
           ENDIF                                                                
           IF (TAUT(JZ).GT.0.0) TOTALT = TOTALT + 1.0 / TAUT(JZ)                
  300    CONTINUE                                                               
C                                                                               
         ARG = QTIM * QS(IQX) * TOTALT                                          
         IF (ARG.GT.1.E-3) THEN                                                 
           CFTS(IX,IY,IZ) = 1.0 - EXP(-ARG)                                     
         ELSE                                                                   
           CFTS(IX,IY,IZ) = ARG                                                 
         ENDIF                                                                  
C                                                                               
         IF (10*(IY/10).EQ.IABS(IY) .AND. IY.LE.50 .AND.                        
     >       10*(IX/10).EQ.IX .AND. 2*(IZ/2).NE.IZ) THEN                        
          WRITE (6,9002) IY,IX,1.0/TOTALP,(TAUP(JZ),JZ=0,NIZS)                  
          WRITE (6,9003)       1.0/TOTALS,(TAUS(JZ),JZ=0,NIZS)                  
          WRITE (6,9004)       1.0/TOTALT,(TAUT(JZ),JZ=0,NIZS)                  
         ENDIF                                                                  
  510   CONTINUE                                                                
  520  CONTINUE                                                                 
  530 CONTINUE                                                                  
C                                                                               
      CALL TAUPRA (7,ZEFFS(1,-MAXNYS,5),'ZB.NBT (NT BASED) (M**-3)',-1)         
      RETURN                                                                    
 9001 FORMAT(//1X,'TAUIN3: IONIZATION STATE',I3,///1X,                          
     >  '  IY  IX        TOTAL    TAUB',9(5X,'TAU',I1),/1X,130('-'))            
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
C                                                                               
C***********************************************************************        
C   SHORT ROUTINE TO PRINT A LINE OF THE "CHARACTERISTIC TIMES" TABLE           
C   OF TAU PARA, TAU STOP AND TAU HEAT.                                         
C                          CHRIS FARRELL     DECEMBER 1987                      
C***********************************************************************        
C                                                                               
c      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
c      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      CHARACTER*11 WTAUP,WTAUS,WTAUT                                            
      REAL         TAUP,TAUS,TAUT                                               
C                                                                               
C     WRITE (6,'('' TP,TS,TT='',3G11.4)') CFPS(IX,1,IZ),                        
C    >  CFSS(IX,1,IZ),CFTS(IX,1,IZ)                                             
C                                                                               
      IQX = IQXS(IX)                                                            
C                                                                               
      IF (CFPS(IX,IY,IZ).GT.0.0) THEN                                           
        TAUP = 2.0 * CTEMSC * QTIM * QS(IQX) / CFPS(IX,IY,IZ)                   
        WRITE (WTAUP,'(1P,G11.2)') TAUP                                         
      ELSE                                                                      
        WTAUP = ' INFINITE  '                                                   
      ENDIF                                                                     
C                                                                               
      IF     (CIOPTC.EQ.2) THEN                                                 
        WTAUS = WTAUP                                                           
      ELSEIF (CFSS(IX,IY,IZ).GE.1.0) THEN                                       
        WTAUS = ' INSTANT   '                                                   
      ELSEIF (CFSS(IX,IY,IZ).GT.1.E-3) THEN                                     
        TAUS = -QTIM * QS(IQX) / LOG (1.0 - CFSS(IX,IY,IZ))                     
        WRITE (WTAUS,'(1P,G11.2)') TAUS                                         
      ELSEIF (CFSS(IX,IY,IZ).GT.0.0) THEN                                       
        TAUS = QTIM * QS(IQX) / CFSS(IX,IY,IZ)                                  
        WRITE (WTAUS,'(1P,G11.2)') TAUS                                         
      ELSE                                                                      
        WTAUS = ' INFINITE  '                                                   
      ENDIF                                                                     
C                                                                               
      IF     (CFTS(IX,IY,IZ).GE.1.0) THEN                                       
        WTAUT = ' INSTANT   '                                                   
      ELSEIF (CFTS(IX,IY,IZ).GT.1.E-3) THEN                                     
        TAUT = -QTIM * QS(IQX) / LOG (1.0 - CFTS(IX,IY,IZ))                     
        WRITE (WTAUT,'(1P,G11.2)') TAUT                                         
      ELSEIF (CFTS(IX,IY,IZ).GT.0.0) THEN                                       
        TAUT = QTIM * QS(IQX) / CFTS(IX,IY,IZ)                                  
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
      use mod_soledge
      use mod_sol22_input_lim
      use mod_sol22_lim
      use mod_solcommon
      use mod_comxyt
      use mod_comt2
      use yreflection
      use mod_vtig
      use mod_comtor
      use variable_wall
      implicit none
      real :: qtim

c
      integer :: ix,iy,iqx,iqy,ip
      real :: ti,t0,vel,n,x,y,y0
      integer :: pz,yz
      integer :: ixout
      integer,external :: ipos


      
      !     the setup_vtig routine assigns masses and calculates the integration constant
      ! and should be called for all vtig options - this is needed to calculate an estimate
      ! of vtig from the temperature gradients and should be called in all cases 
      call setup_vtig(crmb,crmi,cnbin,ctibin)

c
      
      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
c     
c       
c     IQYS and IQXS should be setup to map IX and IY to IQX ad IQY
c     Need to include multiplying by the radial scale factors  ?
c     Consider removing CVHXS and CVEXZS      
c      
      do pz = 1,maxpzone
         do ix = 1,ixout
            do iy = -nys,nys
               if (iy.lt.0) then
                  iqy = iqys(iy+nys+1)
               elseif (iy.eq.0) then
                  iqy = 1
               else
                  iqy = iqys(iy)
               endif
               efield(ix,iy,pz) = ceys(iqy)
               velplasma(ix,iy,pz) = cvhys(iqy)
            end do
         end do
         do ix = ixout+1,nxs
            do iy = -nys,nys
               iqx = iqxs(ix)
               if (iy.lt.0) then
                  iqy = iqys(iy+nys+1)
               elseif (iy.eq.0) then
                  iqy = 1
               else
                  iqy = iqys(iy)
               endif
               efield(ix,iy,pz) = ceyin              !ceys(iqy)  
               velplasma(ix,iy,pz) = cvhyin          !cvhys(iqy)
            end do
         end do
      end do

c
c     If the collector probe 3D plamsa options are in effect then call the
c     code to set up the modified plasma, efield and plasma velocity arrays
c     
c     sazmod - Maybe use a separate switch for this statement to allow
c              only setting up forces in lim3.f without prescribing a 
c              complex SOL (like SOL12, 13, etc.). 
       if (soledge_opt.eq.1.and.colprobe3d.eq.1) then 

         if (vary_2d_bound.eq.1) then
         
           write(0,*) 'SOLEDGE1D: Plasma with 2D boundary...'
         
           ! sazmod - the 2D array of the fully customizable boundary
           ! which as of this writing is only applied to yabsorb1a side,
           ! is implemented "on top" of yabsorb1a. It in effect adds
           ! absorbing surfaces on top of it to extend the boundary and
           ! decrease the connection length to the specified value.
           
           ! The way PS is filled in is that the specified bin 
           ! boundaries in the input file are essentially placed in the
           ! middle of the PS array, and the unassigned cells on the 
           ! edges of the array are filled in via smaller steps. I.e.,
           ! if the lowest input P bin is -5, then we may see the array
           ! looking like [-5.04, -5.03, -5.02, -5.01, 5.00, ...]
           ! depending on how many extra indices need to be filled (PS is
           ! 2*MAXNPS+1 big, so 2*MAXNPS+1-NPBINS indices are filled in, 
           ! split evenly among each edge of PS). All this is to say,
           ! the 2D customizable bound option will ONLY work when NPBINS
           ! = 2*MAXNPS+1, otherwise the indexing between the two becomes
           ! a mess. So if you get this error, go into mod_params_lim
           ! and adjust MAXNPS and recompile (make clean, make).
           if (npbins.ne.(2*maxnps+1)) then
             write(0,*) 'WARNING: NPBINS != 2*MAXNPS+1. The varying 2D'
             write(0,*) 'boundary option will not work correctly!'
             write(0,*) 'Change MAXNPS to ',(npbins-1) / 2,'.'
           endif
                      
           ! This is probably pretty expensive memory-wise, but we will 
           ! need to go through one flux tube at a time to account for 
           ! each individual connection length and use soledge for each 
           ! tube (hence the 1d suffix here). Each ip, ix has a 
           ! corresponding absorbing boundary distance.
           !write(6,*) 'bounds begin'
           write(6,*) 'maxnps = ',maxnps
           write(6,*) 'npbins = ',npbins
           write(6,*) 'yabsorb2a = ',yabsorb2a
           do ip=1, npbins 
             do ix=1, nxs   
               
               ! From bounds pull out the location of the absorbing 
               ! boundary for this flux tube. Initialize soledge arrays.
!               write(6,*) 'ix, ip, X, P, bound = ',ix,ip,xs(ix),
!     >           ps(-maxnps+ip-1),bounds(ix,ip)
               call init_soledge(bounds(ix, ip), yabsorb2a)
               call soledge_1d(ix, ip, qtim)
               
             enddo
           enddo
           !write(6,*) 'bounds end'
           
!           write(6,*) 'ctembs_3d begin at iy = 0'
!           do ix=1, nxs
!             do ip=1, npbins
!               write(6,*) 'ix, ip, X, P, ctembs_3d = ',ix,ip,xs(ix),
!     >           ps(-maxnps+ip),ctembs_3d(ip,ix,0)
!             end do
!           end do
!           write(6,*) 'ctembs_3d end'
    
         else
             ! plasma is calculated from lower absorbing surface to
             ! upper absorbing surface - this allows for
             ! asymmetric placement of the probe
             call init_soledge(yabsorb1a,yabsorb2a)
             call soledge(1,nxs,qtim)
         endif
         
         
         
      endif

      
      if (sol22_opt.gt.0.and.nsol22_opt.gt.0) then 
c
c     This code calculates the plasma conditions for sections of the simulation
c     volume using SOL22. There are several scenarios.
c
c     1) No absorbing surfaces - standard limiter or probe simulations
c     SOL22 is used to calculate the background plasma on field lines
c     that connect to the probe/limiter. (i.e. PZONE = 1)
c              
c         
! Initialize some output options in SOL22 using values from slcom
         call init_solcommon(0,0)

         call sol22

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

                  n = crnbs(ix,iy)
                  y0 = qedges(iqx,2)
                  t0 = qtembsi(iqx,2)

                  if (y.lt.y0) then
                     ctembsi(ix,iy) = t0
                  else
                     call calculate_temperature(x,y-y0,pz,yz,n,t0,
     >                 cl-y0,ti,
     >                 n_vtig_blocks,vtig_range,vtig_ndata,
     >                 vtig_data,vtig_zones)
                     ctembsi(ix,iy) = ti
                  endif 
               end do
               do iy =  -nys/2,-1
                  iqx = iqxs(ix)
                  
                  x = xouts(ix)
                  yz = int(sign(1.0,youts(iy)))

                  y = abs(youts(iy))
                  n = crnbs(ix,iy)
                  y0 = qedges(iqx,1)
                  t0 = qtembsi(iqx,1)
                  
                  call calculate_temperature(x,y-y0,pz,yz,n,t0,
     >                 cl-y0,ti,
     >                 n_vtig_blocks,vtig_range,vtig_ndata,
     >                 vtig_data,vtig_zones)
                  ctembsi(ix,iy) = ti
               end do

               ! Copy the central portion to the rest of the range
               do iy = -nys,-nys/2-1
                  ctembsi(ix,iy) = ctembsi(ix,iy+nys+1)
               end do
               do iy = nys/2+1,nys
                  ctembsi(ix,iy) = ctembsi(ix,iy-nys-1)
               end do
               ! average zero value between +/- 1
               ctembsi(ix,0) = (ctembsi(ix,1)+ctembsi(ix,-1))/2.0
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
!                 velocity profiles are entered as if for the first Y>0, Y<CL section
!                 so the sign is swapped for the Y<0 Y>-CL section                    
!                 the code in calculate_velocity uses the value of yz to determine the sign
!                 of the returned velocity. vel = yz * vel (flow towards the target in Y>0, Y< L
!                 is negative
!                  velplasma(ix,iy,pz) = -vel
                  velplasma(ix,iy,pz) = vel
               end do

               ! Copy the central portion to the rest of the range
               do iy = -nys,-nys/2-1
                  velplasma(ix,iy,pz) = velplasma(ix,iy+nys+1,pz) 
               end do
               do iy = nys/2+1,nys
                  velplasma(ix,iy,pz) = velplasma(ix,iy-nys-1,pz) 
               end do
               ! average zero value between +/- 1
               velplasma(ix,0,pz) = (velplasma(ix,1,pz)
     >                             +velplasma(ix,-1,pz))/2.0
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
      
      ! locals
      integer :: j,ix,iy,iqx
      INTEGER   IQXCV1,IQXCV2,iqxbrk
      integer,external :: ipos
      REAL      RDX,NX,dnx
      
      
C-----------------------------------------------------------------------
C  SET DIFFUSION DECAY ETC, USING DPERP FACTORS                                 
C  ENSURE INWARD STEPPING PROBABILITIES ARE WITHIN RANGE (0,1)                  
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
     >                     QTIM * QS(IQX) / (CA*CA)                               
          ELSEIF (CVPOPT.EQ.1) THEN 
            CXAFS(IQX,J) = VPV0*(((CA-QXS(IQX))/CA)**VPALPH) *
     >                     QTIM * QS(IQX)
          ELSEIF (CVPOPT.EQ.2) THEN
            IF (QXS(IQX).GE.CVPCUT) THEN 
               CXAFS(IQX,J) = 0.0
            ELSEIF (IQX.EQ.NQXSI) THEN                
               CXAFS(IQX,J) = 0.0
            ELSE
              IX = IPOS(QXS(IQX),XS,NXS)
              IF (IX.LE.NXS) THEN 
                IY = 1
                DNX = (CRNBS(IX,IY)-CRNBS(IX-1,IY))/(XS(IX)-XS(IX-1))
                NX = CRNBS(IX-1,IY) + DNX * (QXS(IQX)-XS(IX-1))               
                CXAFS(IQX,J) = CVIN * CRDXO(J) * DNX / NX * 
     >                         QTIM * QS(IQX)
              ELSE
                CXAFS(IQX,J) = 0.0
              ENDIF 
            ENDIF
          ENDIF 
C
C         CXAFS(IQX,J) = 2.0 * CRDXO(J) * CVIN * (CA-QXS(IQX)) *                
C    >                   QTIM * QS(IQX) / (CA*CA)                               
C

c         sazmod
c         Will try and fit in the code for just specifying the radial 
c         diffusion into regions here, overwriting whatever happens 
c         above. Essentially all we want is to swap CRDXO(J) out with the
c         correct dperp_reg. Look at the picture in unstructured_input.f
c         if someone besides me is looking at this (hello from the past!).

c         Fortunately, we already know that if J=1 or 3 it's in the left
c         side of things, so either region 1 or 3. We just need to know
c         the X (i.e. radial) value to see if it's in the step region 
c         or not.
         
          !write(0,*) 'IQX, QXS = ', IQX, QXS(IQX)
          
          if (dperp_reg_switch.eq.1) then
            
            ! Left region.
            if ((J.eq.1).or.(J.eq.3)) then
            
              ! See if in step part or not.
              if (QXS(IQX).lt.xabsorb1a_step) then
              
                ! If in step part, region 3.
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg3 * QTIM * QS(IQX))                 
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg3 
     >                         / (CDPSTP * CDPSTP) 
              else
              
                ! If not in step part, region 1.
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg1 * QTIM * QS(IQX))                 
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg1 
     >                         / (CDPSTP * CDPSTP) 
              endif
             
            ! Right region.
            else
              
              ! If in step part, region 4.
              if (QXS(IQX).lt.xabsorb2a_step) then
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg4 * QTIM * QS(IQX))                 
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg4 
     >                         / (CDPSTP * CDPSTP) 
     
              ! If not in step part, region 2.
              else
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg2 * QTIM * QS(IQX))                 
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg2 
     >                         / (CDPSTP * CDPSTP) 
              endif
            endif
              
          else

            CXBFS(IQX,J) = SQRT (2.0 * CRDXO(J) * QTIM * QS(IQX))                 
            CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * CRDXO(J) 
     >                      / (CDPSTP * CDPSTP) 
         endif

     
  100   CONTINUE                                                                
C                                                                               
        DO 110 IQX = IQXBRK+1, NQXSI                               
          IF (CDPERP.EQ.0) THEN 
            RDX = CRDXI(J) + CRDD * QXS(IQX) / CA                                 
          ELSEIF (CDPERP.EQ.1) THEN 
            RDX = CRDXI(J)* (1.0+DPALPH*((CA-QXS(IQX))/CA)**DPBETA)
          ENDIF
          IF (CVPOPT.EQ.0) THEN
            CXAFS(IQX,J) = 2.0 * RDX * CVIN * (CA-QXS(IQX)) *                     
     >                     QTIM * QS(IQX) / (CA*CA)                               
          ELSEIF (CVPOPT.EQ.1) THEN 
            CXAFS(IQX,J) = VPV0*(((CA-QXS(IQX))/CA)**VPALPH) *
     >                     QTIM * QS(IQX)
          ELSEIF (CVPOPT.EQ.2) THEN
            IF (QXS(IQX).GE.CVPCUT) THEN 
               CXAFS(IQX,J) = 0.0
            ELSEIF (IQX.EQ.NQXSI) THEN                
               CXAFS(IQX,J) = 0.0
            ELSE
              IX = IPOS(QXS(IQX),XS,NXS)
              IF (IX.LE.NXS) THEN 
                IY = 1
                DNX = (CRNBS(IX,IY)-CRNBS(IX-1,IY))/(XS(IX)-XS(IX-1))
                NX = CRNBS(IX-1,IY) + DNX * (QXS(IQX)-XS(IX-1))                     
                CXAFS(IQX,J) = CVIN * RDX * DNX / NX * QTIM * QS(IQX)
              ELSE
                CXAFS(IQX,J) = 0.0
              ENDIF 
            ENDIF
          ENDIF 
          
          !sazmod
          ! A little different for the inboard side compared to the
          ! outboard, but still simple enough. This time we are just
          ! replacing RDX. 
          !write(0,*) 'IQX, QXS = ', IQX, QXS(IQX)
          if (dperp_reg_switch.eq.1) then
          
            ! Left region.
            if ((J.eq.1).or.(J.eq.3)) then
            
              ! If in step part, region 3.
              if (QXS(IQX).lt.xabsorb2a_step) then
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg3 * QTIM * QS(IQX))                      
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg3 
     >                         / (CDPSTP * CDPSTP) 
     
              ! If not in step part, region 1.
              else
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg1 * QTIM * QS(IQX))                      
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg1 
     >                         / (CDPSTP * CDPSTP) 
              endif
            
            ! Right region.
            else  
            
              ! If in step part, region 4.
              if (QXS(IQX).lt.xabsorb2a_step) then
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg4 * QTIM * QS(IQX))                      
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg4 
     >                         / (CDPSTP * CDPSTP) 
     
              ! If not in step part, region 2.
              else
                CXBFS(IQX,J) = SQRT (2.0 * dperp_reg2 * QTIM * QS(IQX))                      
                CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * dperp_reg2 
     >                         / (CDPSTP * CDPSTP) 
              endif
              
            endif
            
          else
            CXBFS(IQX,J) = SQRT (2.0 * RDX * QTIM * QS(IQX))                      
            CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * RDX 
     >                      / (CDPSTP * CDPSTP) 
          endif
       
  110   CONTINUE                                                                
C                                                                               
C       An outward pinch velocity or an arbitrary pinch over a specified 
C       region may also be entered. This is simply added to the contents
C       of CXAFS - the array containing the step due to the pinch velocity
C       component of the motion.
C
C       Note: The definition of signs in the code differs from that used 
C             in the literature. Normally, an inward pinnch is descibed by
C             Vpinch = -2.0 Dperp * r / a**2 - in the code the "+" ve sign
C             is used for motion towards the core - the "-" ve sign is used 
C             for drift towards the SOL. Thus an arbitrary pinch of +0.5 m/s 
C             in the specification list -  is a drift towards the walls and
C             is applied to the CXAFS array as   "-" Vin * QTIM * QS(IQX) 
C  
        DO 115 IQX = IQXCV1,IQXCV2
           CXAFS(IQX,J) = CXAFS(IQX,J) + CVPOUT * QTIM * QS(IQX)
 115    CONTINUE
C
        DO 120 IQX = 1-NQXSO, NQXSI                                             
          IF     (IGEOM.EQ.0) THEN                                              
            CXCFS(IQX,J) = 0.5                                                  
          ELSEIF (IGEOM.EQ.1) THEN                                              
            CXCFS(IQX,J) = (CA-QXS(IQX)-0.5*CXBFS(IQX,J)) /                     
     >                     (2.0*(CA-QXS(IQX)))                                  
            CXCFS(IQX,J) = MIN (1.0, MAX (0.0, CXCFS(IQX,J)))                   
          ENDIF                                                                 
  120   CONTINUE                                                                
  130 CONTINUE                                                                  


      return
      end
