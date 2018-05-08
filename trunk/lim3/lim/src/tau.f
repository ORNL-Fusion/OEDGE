      SUBROUTINE TAUIN1 (QTIM,NIZS,ICUT,                                        
     >                   FSRATE,IGEOM,NTBS,NTIBS,NNBS)                         
      use mod_comt2
      use variable_wall
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
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c
      include   'global_options'
c
c slmod
      INCLUDE   'cadas'
      INCLUDE   'slcom'
c slmod end
C                                              ,ICNT
c slmdo begin
      INTEGER   IPOS,IXOUT,IZ,IQX,LIMIZ,IX,IY,J                                 
      INTEGER   IQXCV1,IQXCV2 
      REAL      RDX,FEX,WIDTH,FEXZ,DNX,NX                              
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
C                     SET UP QEDGES, QTANS AND QDISTS                           
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (6,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      WRITE (0,'('' TAU: CALLING EDGE   OPTION'',I3)') CIOPTH                   
      CALL EDGE (QXS,QEDGES,QTANS,QDISTS,NQXSO,CAW,CL,ICUT,CCUT,XSCALO,         
     >           WEDGAN,XL1,YL1,XL2,YL2,TC,SC,TO,SO,GC,RP,CIOPTH,CORECT,
     >           XST1,YST1,XST2,YST2,XST3,YST3,RLEDGE7,CA,RLC)        

C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP VARIABLE CAW
C-----------------------------------------------------------------------        
c
      call setup_wall (qys,nqys,cl,caw)
C                                                                               
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
          CXBFS(IQX,J) = SQRT (2.0 * CRDXO(J) * QTIM * QS(IQX))                 
          CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * CRDXO(J) 
     >                   / (CDPSTP * CDPSTP) 
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
          CXBFS(IQX,J) = SQRT (2.0 * RDX * QTIM * QS(IQX))                      
          CXDPS(IQX,J) = 2.0 * QTIM * QS(IQX) * RDX
     >                   / (CDPSTP * CDPSTP) 
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
C                                                                               
C-----------------------------------------------------------------------        
C                     SET UP CEYS AND CVHYS                                     
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (6,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      WRITE (0,'('' TAU: CALLING SOL    OPTION'',I3)') CIOPTF                   
      CALL SOL (QYS,CEYS,CVHYS,NQYS,CTBIN,CTIBIN,CRMB,CL,CIZB,                 
     >          CEYOUT,CVHOUT,CYSTAG,CRMI,CSOLEF,CIOPTF)                        
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
           DO 250 IX = 1, IXOUT                                                 
            IQX = IQXS(IX)                                                      
            CFEXZS(IX,IY,IZ) = FEXZ * CTEMBS(IX,IY)/CTBIN *                     
     >                         CYSCLS(IQX)/YSCALE * QS(IQX) * QS(IQX)           
  250     CONTINUE                                                              
  300   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 310 IY = -NYS, NYS                                                     
       DO 310 IX = 1, IXOUT                                                     
        IQX = IQXS(IX)                                                          
        CFVHXS(IX,IY) = 
     >     SQRT((CTEMBS(IX,IY)+CTEMBSI(IX,IY))/(CTBIN+CTIBIN))
     >     * QTIM * QS(IQX)             
  310 CONTINUE                                                                  
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
      CALL IZTAU (CRMI,NXS,NYS,CION,CIZB,CIOPTA)                                
C                                                                               
C-----------------------------------------------------------------------        
C    SET COMBINED C-X AND E-I RECOMBINATION TIMES         CFCXS                 
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (0,'('' TAU: CALLING CXREC  OPTION'',I3)') CIOPTI                   
      CALL CXREC (NIZS,CION,CIOPTI,CIZB,CL,CRMB,CVCX,                           
     >            CNHC,CNHO,CLAMHX,CLAMHY)                                      
C                                                                               
C-----------------------------------------------------------------------        
C     SET PROBABILITY OF EITHER AN IONISATION OR A RECOMBINATION                
C     SET PROPORTION OF THESE WHICH WILL BE RECOMBINATIONS                      
C     PREVENT ANY IONISATION BEYOND MAXIMUM LIMIT SPECIFIED IF REQUIRED         
C-----------------------------------------------------------------------        
C                                                                               
c slmod tmp
c
c      IF (CIOPTE.EQ.10.AND.CDATOPT.EQ.0) THEN
c        DO IZ = 1, CMIZS-1
c          DO IX = 1, NXS
c            DO IY = -NYS,NYS
c              CFIZS(IX,IY,IZ) = 1.0 / (TMPION*CNBIN)
c            ENDDO
c          ENDDO
c        ENDDO
c      ENDIF
c slmod end
      IF (CMIZS.GT.1) THEN                                                      
        DO 370 IZ = 1, CMIZS-1                                                  
         DO 360 IX = 1, NXS                                                     
          IQX = IQXS(IX)                                                        
          DO 350 IY = -NYS, NYS                                                 
            IF (CFCXS(IX,IY,IZ).LE.0.0) THEN                                    
              CPCHS(IX,IY,IZ) = QTIM * QS(IQX) / CFIZS(IX,IY,IZ)                
              CPRCS(IX,IY,IZ) = 0.0                                             
            ELSE                                                                
              CPCHS(IX,IY,IZ) = (CFIZS(IX,IY,IZ) + CFCXS(IX,IY,IZ)) *           
     >          QTIM * QS(IQX) /(CFIZS(IX,IY,IZ) * CFCXS(IX,IY,IZ))             
              CPRCS(IX,IY,IZ) = CFIZS(IX,IY,IZ) /                               
     >                          (CFCXS(IX,IY,IZ) + CFIZS(IX,IY,IZ))             
            ENDIF                                                               
            CPCHS(IX,IY,IZ) = MIN (1.0, CPCHS(IX,IY,IZ))                        
  350     CONTINUE                                                              
  360    CONTINUE                                                               
  370   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      IF (CMIZS .LE. LIMIZ) THEN                                                
        DO 390 IX = 1, NXS                                                      
          IQX = IQXS(IX)                                                        
          DO 380 IY = -NYS, NYS                                                 
            IF (CFCXS(IX,IY,IZ).LE.0.0) THEN                                    
              CPCHS(IX,IY,CMIZS) = 0.0                                          
            ELSE                                                                
              CPCHS(IX,IY,CMIZS) = QTIM * QS(IQX) / CFCXS(IX,IY,IZ)             
            ENDIF                                                               
            CPRCS(IX,IY,CMIZS) = 1.0                                            
            CPCHS(IX,IY,IZ) = MIN (1.0, CPCHS(IX,IY,IZ))                        
  380     CONTINUE                                                              
  390   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- SET IONISATION PROBABILITIES FOR NEUT ...                                 
C---- (SAVES REPEATED CALCULATION EVERY ITERATION)                              
C---- THEY ARE ALL MULTIPLIED BY THE "IONISATION RATE FACTOR" IRF               
C---- WHICH TYPICALLY MIGHT BE 0.2 TO GIVE DEEPER IONISATION.                   
C                                                                               
      DO 500 IY = -NYS, NYS                                                     
        DO 500 IX = 1, NXS                                                      
          CPCHS(IX,IY,0) = MIN (1.0, CIRF * FSRATE / CFIZS(IX,IY,0))            
  500 CONTINUE                                                                  
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
      use mod_comt2
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
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
C                                                                               
      INTEGER   IPOS,IZ,IQX,LIMIZ,JX,IX,IY                                      
      REAL      LAMBDA,ROOTMI,ROOTTT                                            

      REAL      TEMP,FTAU,FTAUP,FTAUS,FTAUT,RIZSQR,STAU,TAU                     
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
      FTAU  = CZENH * SQRT(CRMB) * CIZB * CIZB * LAMBDA * QTIM                  
      FTAUP = FTAU * 6.8E-14                                                    
      FTAUS = FTAU * 6.8E-14 * (1.0 + CRMB/CRMI)                                
      FTAUT = FTAU * 1.4E-13                                                    
      C215A = FTAU * 1.4E-13 * SQRT(CRMI)                                       
      C215B = (CRMI * CTBI + CRMB * CTEMSC) ** 1.5                              
C                                                                               
      DO 540  IZ = 1, LIMIZ                                                     
         RIZSQR = REAL (IZ) * REAL (IZ)                                         
         DO 520 IY = -NYS, NYS                                                  
          DO 520  IX = 1, NXS                                                   
            IF (CTBI.GT.0.0) THEN                                               
             TEMP  = CRNBS(IX,IY) / (CRMI * CTBI**1.5)                          
            ELSE                                                                
             TEMP  = CRNBS(IX,IY) / (CRMI * CTEMBSI(IX,IY)**1.5)                
            ENDIF                                                               
            IQX = IQXS(IX)                                                      
            STAU = TEMP * RIZSQR * QS(IQX)                                      
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
              CFPS(IX,IY,IZ) = STAU * CTBI * FTAUP * 2.0                        
            ELSE                                                                
              CFPS(IX,IY,IZ) = STAU * CTEMBSI(IX,IY) * FTAUP * 2.0             
            ENDIF                                                               
C                                                                               
            IF     (CIOPTB.EQ.1.AND.IX.LE.JX) THEN                              
               CFPS(IX,IY,IZ) = 0.0                                             
            ELSEIF (CIOPTB.EQ.2.AND.IX.LE.JX) THEN                              
               CFPS(IX,IY,IZ) = 2.0 * CRNBS(IX,IY) * 6.8E-14 *                  
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
            IF (CFPS(IX,IY,IZ).EQ.0.0) THEN                                     
               CCCFPS(IX,IY,IZ) = 0.0                                           
            ELSEIF (CIOPTB.EQ.3 .AND. IX.LE.JX) THEN 
               CCCFPS(IX,IY,IZ) = SQRT (9.76E8 * CTEMSC / CRMI) *               
     >           QTIM * QS(IQX) / CFPS(IX,IY,IZ)                                
            ELSEIF (CIOPTB.EQ.4 .AND. IX.LE.JX) THEN
               CFPS(IX,IY,IZ) = 2.0E0 * CFPS(IX,IY,IZ)
               CCCFPS(IX,IY,IZ) = SQRT(9.76E8 * CTEMSC / CRMI ) *
     >            QTIM * QS(IQX) / CFPS(IX,IY,IZ)
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

               CCCFPS(IX,IY,IZ) = 
     +                  SQRT(2.0*1.6E-19/1.67E-27/CRMI)*
     +                  SQRT(QTIM/
     +                   (CRMI*SQRT(CTEMBSI(IX,IY)/CRMB)/6.8E-14
     +                    /LAMBDA/CRNBS(IX,IY)/(1+CRMB/CRMI)/
     +                    (REAL(CIZB)*REAL(CIZB))/
     +                    (REAL(IZ)*REAL(IZ))) )

c            WRITE (78,*) CCCFPS(IX,IY,IZ),CRNBS(IX,IY),CTEMBSI(IX,IY),
c     +                   QTIM,CRMB,CRMI,CIZB,REAL(IZ),LAMBDA 
c slmod end
            ELSE                                                                
               CCCFPS(IX,IY,IZ) = SQRT (4.88E8 /(CFPS(IX,IY,IZ)*CRMI))*         
     >           QTIM * QS(IQX)                                                 
            ENDIF                                                               
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
              CFSS(IX,IY,IZ) = 1.0 - EXP(-TAU)                                  
            ELSE                                                                
              CFSS(IX,IY,IZ) = TAU                                              
            ENDIF                                                               
C                                                                               
            IF     (CIOPTC.EQ.1.AND.IX.LE.JX) THEN                              
               CFSS(IX,IY,IZ) = 0.0                                             
            ELSEIF (CIOPTC.EQ.2.AND.IX.LE.JX) THEN                              
               TAU = CFPS(IX,IY,IZ) / (2.0*CTEMSC)                              
               IF (TAU.GT.1.E-3) THEN                                           
                 CFSS(IX,IY,IZ) = 1.0 - EXP(-TAU)                               
               ELSE                                                             
                 CFSS(IX,IY,IZ) = TAU                                           
               ENDIF                                                            
            ELSEIF (CIOPTC.EQ.3 .AND. IX.LE.JX) THEN
                 CFSS(IX,IY,IZ) = 1.0E20 
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
              CFTS(IX,IY,IZ) = 1.0 - EXP(-TAU)                                  
            ELSE                                                                
              CFTS(IX,IY,IZ) = TAU                                              
            ENDIF                                                               
C                                                                               
            IF     (CIOPTD.EQ.1.AND.IX.LE.JX) THEN                              
               CFTS(IX,IY,IZ) = 0.0                                             
            ELSEIF (CIOPTD.EQ.2.AND.IX.LE.JX) THEN                              
               CFTS(IX,IY,IZ) = 1.0                                             
            ELSEIF (CIOPTD.EQ.3.AND.IX.LE.JX) THEN                              
               IF (CTBI.LE.0.0)                                                 
     >           C215B = (CRMI * CTEMBSI(IX,IY)+ CRMB*CTEMSC) ** 1.5      
               TAU = C215A * RIZSQR * QS(IQX) * CRNBS(IX,IY) / C215B            
               IF (TAU.GT.1.E-3) THEN                                           
                 CFTS(IX,IY,IZ) = 1.0 - EXP(-TAU)                               
               ELSE                                                             
                 CFTS(IX,IY,IZ) = TAU                                           
               ENDIF                                                            
            ENDIF                                                               
520      CONTINUE                                                               
540   CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE TAUPR1 (QTIM,NIZS)                                             
      use mod_comt2
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
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
      INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c
      include   'global_options'
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
      IF (CPRINT.EQ.1) THEN                                                     
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
      IF (CPRINT.EQ.1) THEN                                                     
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
      IF (CPRINT.EQ.1) THEN                                                     
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
      use mod_comt2
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
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
      INCLUDE   'printr'                                                        
C     INCLUDE   (PRINTR)                                                        
c
      include   'global_options'
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
      use mod_comt2
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
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
c      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
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
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      REAL AS(MAXNXS,-MAXNYS:MAXNYS)                                            
      INTEGER JW,ISTATE                                                         
      CHARACTER*(*) NAME                                                        
      INCLUDE 'printr'                                                          
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
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
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
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
c      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'coords'                                                        
C     INCLUDE   (COORDS)                                                        
      INCLUDE   'comxyt'                                                        
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
      use mod_comt2
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
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
c      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
      INCLUDE 'comxyt'                                                          
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
