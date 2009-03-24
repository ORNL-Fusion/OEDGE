c     -*-Fortran-*-
c
      SUBROUTINE NEUT (NATIZ,FSRATE,RRES,                                       
     >                 ICUT,MATLIM,MAT1,MAT2,NPROD,IMPADD,IMPCF,
     >                 QTIM,GTOT1,GYTOT1,
     >                 RSTRUK,              
     >                 RATIZ,RNEUT,RWALLN,RCENT,RTMAX,SEED,NRAND,               
     >                 NEUTIM,RFAIL,NYMFS,NCVS,STATUS)                          
      use variable_wall
      implicit none                                                    
      DOUBLE PRECISION SEED                                                     
      INTEGER   NRAND,NATIZ,ICUT(2),MATLIM,NPROD,NYMFS,STATUS                   
      INTEGER   IMPADD,IMPCF
      INTEGER   MAT1,MAT2
      integer   ncvs 
      REAL      RFAIL,RSTRUK,NEUTIM,GTOT1,GYTOT1                                
      REAL      FSRATE,QTIM,RATIZ,RNEUT,RWALLN,RCENT,RTMAX,RRES                 
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  NEUT: CONTROL ROUTINE FOR SETTING UP PRIMARY NEUTRALS            *        
C  *  -----------------------------------------------------            *        
C  *                                                                   *        
C  *    THIS ROUTINE CREATES A SET OF PRIMARY NEUTRALS DETAILS TO BE   *        
C  *  PASSED TO LAUNCH ROUTINE.  DETAILS CREATED INCLUDE X POSITIONS,  *        
C  *  Y POSITIONS, P POSITIONS, MAXIMUM RANDOM NUMBERS TO BE USED IN   *        
C  *  VELOCITY CALCULATIONS AT LAUNCH POINTS, ETC.  LAUNCH MAY BE      *        
C  *  CALLED A SECOND TIME WHERE THE "SIMPLE SELF-SPUTTERING" OPTION   *        
C  *  HAS BEEN SPECIFIED, WITH FIMP > 0.  IN THIS CASE A SECOND SET OF *        
C  *  NEUTRALS DETAILS IS CREATED AND PASSED TO LAUNCH.                *        
C  *                                                                   *        
C  * INPUT ARGUMENTS :-                                                *        
C  *   ICUT : CUTOFF POINT FOR NEUTRAL INJECTION (POINTER TO QXS ARRAY)*        
C  *   QTIM : TIMESTEP IN LIM2                                         *        
C  *  NPROD : NUMBER OF NEUTRALS TO LAUNCH                             *        
C  *   SEED : RANDOM GENERATOR SEED  (ALSO PASSED BACK TO LIM3)        *        
C  *  NRAND : COUNTS TOTAL RANDOMS USED   (PASSED BACK TO LIM3 TOO)    *        
C  *  NYMFS : NUMBER OF POINTS FOR INTERPOLATING YIELD MODIFIER FUNCTN *        
C  *        : NYMFS IS NO LONGER USED IN NEUT - INITIALIZATION IN LIM  *
C  *   NCVS : NUMBER OF POINTS FOR INTERPOLATING VS(X)                 *
C  * IMPADD : NUMBER OF ADDITIONAL NEUTRALS TO LAUNCH USING ADDITIONAL *
C  *          NEUTRAL LAUNCH OPTIONS 0-FLAT 1-NORMAL DIST.             * 
C  *                                                                   *        
C  * OUTPUTS :-                                                       *        
C  * XATIZS : ARRAY WITH X COORDINATES OF INJECTED IONS - IN /CNEUT/   *        
C  * YATIZS : ARRAY WITH Y COORDINATES OF INJECTED IONS - IN /CNEUT/   *        
C  * PATIZS : ARRAY WITH P COORDINATES OF INJECTED IONS - IN /CNEUT/   *        
C  *   VINS : ARRAY WITH VELOCITIES OF INJECTED IONS    - IN /CNEUT/   *        
C  * SPUTYS : ARRAY WITH FRAGMENT SIZES OF NEUTRALS  (ALL 1'S) /CNEUT/ *        
C  * FSRATE : TIMESTEP USED IN NEUT  (MAY BE SMALLER THAN LIM TIMESTEP)*        
C  *  NATIZ : NUMBER OF NEUTRALS WHICH IONISED                         *        
C  * MATLIM : MATERIAL OF LIMITER REFERENCE FOR "YIELD" FUNCTION       *        
C  *  GTOT1 : TOT INTEGRATED PRIMARY FLUX = 0.5 * (FTOT1(1)+FTOT1(2))  *        
C  * GYTOT1 : TOT INTEGRATED PRIMARY FLUX*YIELD = MEAN OF (1) AND (2)  *        
C  *  RATIZ : TOTAL OF IONISED NEUTRAL FRAGMENTS                       *        
C  *  RNEUT : TOTAL OF LAUNCHED NEUTRAL FRAGMENTS                      *        
C  * RWALLN : TOTAL OF FRAGMENTS PLATING OUT ON WALLS                  *        
C  *  RCENT : TOTAL OF FRAGMENTS REACHING CENTRE                       *        
C  *  RTMAX : TOTAL OF FRAGMENTS EXISTING AT TMAX                      *        
C  * RSTRUK : TOTAL OF FRAGMENTS STRIKING LIMITER                      *        
C  * NEUTIM : TIME SPENT TRACKING NEUTRALS SO FAR (CALC. IN LAUNCH)    *        
C  *  RFAIL : NUMBER OF FAILED NEUTRAL LAUNCHES (V>VMAX TOO MANY TIMES)*        
C  *  IMPCF : NUMBER OF BACKGROUND SPUTTERED ADDITIONAL CROSS-FIELD    *
C  *          NEUTRALS LAUNCHED                                        *
C  *                                                                   *        
C  *                    CHRIS FARRELL (HUNTERSKIL) MARCH 1988          *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'dynam1'                                                          
C     INCLUDE (DYNAM1)                                                          
      INCLUDE 'dynam3'                                                          
C     INCLUDE (DYNAM3)                                                          
      INCLUDE 'cyield'                                                          
C     INCLUDE (CYIELD)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comtau'                                                          
C     INCLUDE (COMTAU)                                                          
      INCLUDE 'comt2'                                                           
C     INCLUDE (COMT2)                                                           
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'printr'                                                          
C     INCLUDE (PRINTR)                                                          
      INCLUDE 'cneut'                                                           
C     INCLUDE (CNEUT)                                                           
      INCLUDE 'comnet'                                                          
C     INCLUDE (COMNET)                                                          
      INCLUDE 'coords'
      EXTERNAL VLAN                                                             
C                                                                               
c      REAL      RADDEG,PI,GAMMA,GAMBL,DELTAX,CS,YIELD,RYIELD                    
      REAL      GAMMA,GAMBL,DELTAX,CS,YIELD,RYIELD                    
      REAL      FTOT2(2),FYTOT2(2),FYTOT(2),FRAC1(2),FTOT1(2),FYTOT1(2)         
      REAL      EMAX1,EMAX2,RAN,X0,Y0,P0,RATIZ1,RNEUT1,RWALL1,RCENT1            
      REAL      RTMAX1,RATIZ2,RNEUT2,RWALL2,RCENT2,RTMAX2,VLAN,RRES1            
      REAL      RSTRK1,RSTRK2,RFAIL1,RFAIL2,THETA,OLDX0,OLDY0,D0,RRES2          
      REAL      RAN1,RAN2,CFBGC,IMPNJ,NINIT
      LOGICAL   PRIME,SPREAD,RESPUT,RES                                         
      INTEGER   IQX,NPROD1,NPROD2,IPROD,IX,IY,NATIZ1,NATIZ2,IPOS,J              
      INTEGER   IP,IT,KQX,IO,IOY,IMPCNT,IERR,ITER
      REAL      TEMP1,TEMP2,TEMPOPT
c slmod
      REAL      PHI
c slmod end
c
c     jdemod  - for wall options
c
      integer :: iqy_tmp 

C                                                                               
c      DATA      RADDEG /57.29577952/, PI /3.141592654/                          
C
      IERR = 0

C                                                                               
C-----------------------------------------------------------------------        
C       SET UP SECTION                                                          
C-----------------------------------------------------------------------        
C                                                                               
C---- CHECK VALUE FOR FIMP IS OK FOR ALTERNATE SPUTTERING OPTIONS.              
C---- SET UP DATA IN COMMON BLOCK CYIELD: SPUTTERING DATA.                      
C                                                                               
      CALL PRB                                                                  
      IF (CNEUTD.NE.2) CFIMP = 0.0                                              
C
C     SET CNEUTC TO THE INITIAL NEUTRAL V/A FLAG FOR USE IN LAUNCH 
C     THEN CHANGE IT BACK TO ITS ORIGINAL VALUE AT THE END OF NEUT
C     SO THAT IT HAS THE CORRECT VALUE FOR LATER SELF-SPUTTERING
C
      TEMPOPT = CNEUTC
      CNEUTC = NVAOPT
C
C
C     MOVE INITIALIZATION OF YIELD DATA TO LIM3 SUBROUTINE SO THAT 
C     IT IS DONE FOR BOTH NEUTRAL AND ION LAUNCHES
C
C     CALL SYIELD (MATLIM,MAT1,MAT2,CNEUTD,                                     
C    >             CBOMBF,CBOMBZ,CION,CIZB,CRMB,CTSUB)                          
C                                                                               
C---- CALCULATE GAMMA, FACTOR DETERMINING MAXIMUM ENERGY EXCHANGE.              
C---- CONSTANT X SPACING IN OUTBOARD MESH FOR IQX = 1-NQXSO,0                   
C---- (DIFFERENT) CONSTANT X SPACING IN INBOARD MESH FOR IQX=1,NQXSI            
C                                                                               
      GAMMA  = 4.0 * CRMB * CRMI / ((CRMB+CRMI) * (CRMB+CRMI))                  
      GAMBL  = GAMMA * (1.0 - GAMMA)                                            
      DELTAX = 1.0 / XSCALO                                                     
C
C     THE ROUTINE TO CALCULATE THE YMF FUNCTIONS HAS BEEN 
C     MOVED TO LIM SO THAT SELF-SPUTTERING BY ION LAUNCHES MAY 
C     TAKE PLACE. NYMFS IS NO LONGER NEEDED AS A PARAMETER, HOWEVER,
C     IT IS LEFT IN CASE THIS SECTION IS EVER RE-INSTATED.
C
C     D.ELDER , NOV 29 , 1990
C
C                                                                               
C---- FIT INTERPOLATING CURVE TO SET OF YIELD MODIFIER VALUES                   
C---- CALCULATE INTERPOLATED/EXTRAPOLATED YMF AT EACH OUTBOARD X POSN.          
C---- DIFFERENT YIELD MODIFIERS FOR EACH SIDE OF Y = 0                          
C---- FLAG DETERMINES WHETHER TO APPLY TO PRIMARIES, SECONDARIES, BOTH          
C                                                                               
      DO 5 J = 1, 2                                                             
        DO 5 IQX = 1-NQXSO, 0                                                   
C
C         CYMFPS(IQX,J) = 1.0                                                   
C         CYMFSS(IQX,J) = 1.0                                                   
C
          CVS(IQX,J) = 1.0
    5 CONTINUE                                                                  
C                                                                              
C     IF (CYMFLG.NE.-2) THEN                                                    
C       CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,2),                               
C    >               NQXSO,QXS(1-NQXSO),CYMFPS(1-NQXSO,1),'LINEAR')             
C       CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,3),                               
C    >               NQXSO,QXS(1-NQXSO),CYMFPS(1-NQXSO,2),'LINEAR')             
C     ENDIF                                                                     
C                                                                              
C     IF (CYMFLG.NE.-1) THEN                                                    
C       CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,2),                               
C    >               NQXSO,QXS(1-NQXSO),CYMFSS(1-NQXSO,1),'LINEAR')             
C       CALL FITTER (NYMFS,CYMFS(1,1),CYMFS(1,3),                               
C    >               NQXSO,QXS(1-NQXSO),CYMFSS(1-NQXSO,2),'LINEAR')             
C     ENDIF                                                                     
C
      IF (CNEUTD.EQ.8) THEN 
         CALL FITTER (NCVS,CVSA(1,1),CVSA(1,2),
     >                NQXSO,QXS(1-NQXSO),CVS(1-NQXSO,1),'LINEAR')
         CALL FITTER (NCVS,CVSA(1,1),CVSA(1,3),
     >                NQXSO,QXS(1-NQXSO),CVS(1-NQXSO,2),'LINEAR')
      ENDIF
C                                                                               
C---- LARMOR RADIUS EFFECT :-  ALLOW BOMBARDING FOR IONS FROM X = +RL           
C---- OUTWARD.  THE LIMITER SURFACE IMPACT POINT FOR THESE IONS IS              
C---- ASSIGNED IN A RANDOM FASHION BETWEEN -RL < X < 0.  A SINGLE BIN           
C---- (IQX=1) IS USED FOR THE +RL AREA - WE HAVE TO ALLOW FOR THE               
C---- DIFFERENCE IN SIZE WITH DELTAX, WHICH IS SIMPLY DONE BY SETTING           
C---- THE YMF VALUE APPROPRIATELY.  THE VALUE KQX IDENTIFIES THE                
C---- BIN IN WHICH -RL LIES, AND IS USED AS THE BASIS FOR SELECTING             
C---- A NEW LAUNCH POINT, WITH SPREADING TEMPORARILY TURNED ON TO GIVE          
C---- A UNIFORM DISTRIBUTION BETWEEN  -Y(RL) < Y < +Y(RL)                       
C                                                                               
      
c      write(6,*) 'Before clarmr ipos'
      KQX = IPOS (-CLARMR, QXS(-NQXSO), NQXSO) - NQXSO - 1                      
c      write(6,*) 'After clarmr ipos'
      IF (CNEUTB.NE.0) CLARMR = 0.0                                             
      IF (CLARMR.GT.0.0) THEN                                                   
        CYMFPS(1,1) = CLARMR / DELTAX                                           
        CYMFPS(1,2) = CLARMR / DELTAX                                           
      ELSE                                                                      
        CYMFPS(1,1) = 0.0                                                       
        CYMFPS(1,2) = 0.0                                                       
      ENDIF                                                                     
      WRITE (6,'('' NEUT: LARMOR KQX,QXS(KQX)'',I5,G11.4)') KQX,QXS(KQX)        
C                                                                               
C                                                                               
C********** LOOP FOR EACH LIMITER SURFACE  1= Y < 0     2= Y > 0   *****        
C                                                                               
C                                                                               
      DO 666 J = 1, 2                                                           
C                                                                               
C---- CALCULATE FLUXES AND YIELDS FOR PRIMARY AND SECONDARY NEUTRALS            
C---- 26/7/88:  ADDED CSINTB FACTOR FOR TOROIDAL CASES.                         
C                                                                               
      DO 10 IQX = 1-NQXSO, 1                                                    
C
C       INCLUDE DEPENDENCE ON BOTH ELECTRON AND ION TEMPERATURES
C       1990, FEB 8,  DAVID ELDER 
C
        CS = 9.79E3 * SQRT (((QTEMBS(IQX,J)+QTEMBSI(IQX,J))/2)*
     >       (1.0+REAL(CIZB))/CRMB)                
        FLUX1(IQX,J)  = QRNBS(IQX,J) * CS * CSINTB                              
        FLUX2(IQX,J)  = FLUX1(IQX,J) * CFIMP                                    
        ENEGY1(IQX,J) = (2.0*QTEMBSI(IQX,J)) + 
     >                  (3.0*REAL(CIZB) * QTEMBS(IQX,J))             
        ENEGY2(IQX,J) = (2.0*QTEMBSI(IQX,J)) + 
     >                  (3.0*REAL(CBOMBZ) * QTEMBS(IQX,J))                
        IF (CNEUTD.EQ.8) 
     >     ENEGY1(IQX,J) = 2.0*QTEMBS(IQX,J) + REAL(CIZB) * CVS(IQX,J)
        IF (CNEUTD.EQ.1) ENEGY1(IQX,J) = ENEGY2(IQX,J)                          
        YIELD1(IQX,J) = YIELD (MAT1,MATLIM,ENEGY1(IQX,J),0.0,0.0)
     >                      *QMULTP*CYMFPS(IQX,J)           
        YIELD2(IQX,J) = YIELD (MAT2,MATLIM,ENEGY2(IQX,J),0.0,0.0) 
     >                     *QMULTS*CYMFSS(IQX,J)        
        IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN       
          YIELD1(IQX,J) = YIELD1(IQX,J)*(1.0+CQPL/
     >                    (CQ(MAT1,MATLIM)*QMULTP))            
          YIELD2(IQX,J) = YIELD2(IQX,J)*(1.0+CQSL/
     >                    (CQ(MAT2,MATLIM)*QMULTS))            
        ENDIF                                                                   
        FY1(IQX,J)    = FLUX1(IQX,J) * YIELD1(IQX,J)                            
        FY2(IQX,J)    = FLUX2(IQX,J) * YIELD2(IQX,J)                            

c        write(6,'(a,2i9,12(1x,g12.5))') 'Y:',iqx,j,
c     >      flux1(iqx,j),flux2(iqx,j),
c     >      enegy1(iqx,j),enegy2(iqx,j),
c     >      yield1(iqx,j),yield2(iqx,j),
c     >      fy1(iqx,j),fy2(iqx,j)
c

   10 CONTINUE                                                                  
C                                                                               
C---- GENERATE CUMULATIVE DISTRIBUTION FUNCTIONS FOR NEUTRAL PRODUCTION         
C                                                                               
      FTOT1(J)  = 0.0                                                           
      FTOT2(J)  = 0.0                                                           
      FYTOT1(J) = 0.0                                                           
      FYTOT2(J) = 0.0                                                           
      DO 20 IQX = ICUT(J), 1                                                    
        FTOT1(J)  = FTOT1(J)  + FLUX1(IQX,J)                                    
        FTOT2(J)  = FTOT2(J)  + FLUX2(IQX,J)                                    
        FYTOT1(J) = FYTOT1(J) + FY1(IQX,J)                                      
        FYTOT2(J) = FYTOT2(J) + FY2(IQX,J)                                      
   20 CONTINUE                                                                  
      FTOT1(J)  = FTOT1(J)  * DELTAX                                            
      FTOT2(J)  = FTOT2(J)  * DELTAX                                            
      FYTOT1(J) = FYTOT1(J) * DELTAX                                            
      FYTOT2(J) = FYTOT2(J) * DELTAX                                            
      FYTOT(J)  = FYTOT1(J) + FYTOT2(J)                                         
      FRAC1(J)  = FYTOT1(J) / FYTOT(J)                                          
      WRITE (6,'('' NEUT: J,FYTOT1,FYTOT2,FYTOT'',I2,3G11.4)')                  
     >  J,FYTOT1(J),FYTOT2(J),FYTOT(J)                                          
C                                                                               
      DO 30 IQX = -NQXSO, 1                                                     
        IF (IQX.LT.ICUT(J)) THEN                                                
          FYCUM(IQX,J)  = 0.0                                                   
          FSPLIT(IQX,J) = 0.0                                                   
        ELSE                                                                    
          FYCUM(IQX,J)  = FYCUM(IQX-1,J) + FY1(IQX,J) + FY2(IQX,J)              
          FSPLIT(IQX,J) = FYCUM(IQX-1,J) + FY1(IQX,J)                           
        ENDIF                                                                   
   30 CONTINUE                                                                  
C                                                                               
      DO 35 IQX = ICUT(J), 1                                                    
        FYCUM(IQX,J)  = FYCUM(IQX,J)  * DELTAX / FYTOT(J)                       
        FSPLIT(IQX,J) = FSPLIT(IQX,J) * DELTAX / FYTOT(J)                       

       WRITE (6,'('' NEUT: IQX,FYCUM,FSPLIT'',I5,1P,2G15.6)') IQX,
     >       FYCUM(IQX,J),FSPLIT(IQX,J)
                                            
   35 CONTINUE                                                                  
      FYCUM(1,J) = 1.0                                                          
C                                                                               
C------ EXTRA CHECK IN CASE OF ROUNDING ERRORS: NO SECONDARIES IF               
C------ FIMP SPECIFIED AS 0.0.                                                  
C                                                                               
      IF (CFIMP.LE.0.0) THEN                                                    
        DO 36 IQX = ICUT(J), 1                                                  
          FSPLIT(IQX,J) = FYCUM(IQX,J)                                          
   36   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- SET UP LIMITING RANDOM NUMBERS FOR CALCULATING LAUNCH VELOCITY.           
C---- THIS HAS TO BE DONE OUTWARDS FROM X=0 SO THAT ONCE WE FIND                
C---- EMAX <= 0,  THE VALUES FROM THE PREVIOUS X POSITION CAN BE USED           
C---- TO PROVIDE A CONSTANT, LOW VALUE FOR RMAX1 OUT TO THE WALL.               
C                                                                               
      IF (CNEUTC.EQ.1) CEMAXF = 1.0                                             
C                                                                               
      DO 15 IQX = 1, 1-NQXSO, -1                                                
        RMAX1(IQX,J) = 1.0                                                      
        RMAX2(IQX,J) = 1.0                                                      
        IF (CNEUTC.EQ.1 .OR. CNEUTC.EQ.4 .OR. CNEUTC.EQ.5
     >      .OR.CNEUTC.EQ.12.OR.CNEUTC.EQ.14) THEN                 
          EMAX1 = CEMAXF * (ENEGY1(IQX,J) * GAMBL-CEBD)                         
          EMAX2 = CEMAXF * ENEGY2(IQX,J)                                        
          IF (CNEUTD.EQ.1) EMAX1 = EMAX2                                        
          IF (EMAX1.GT.0.0) THEN                                                
            RMAX1(IQX,J) = 1.0 /((1.0+CEBD/EMAX1) * (1.0+CEBD/EMAX1))           
          ELSE                                                                  
            RMAX1(IQX,J) = 0.0                                                  
          ENDIF                                                                 
          IF (EMAX2.GT.0.0) THEN                                                
            RMAX2(IQX,J) = 1.0 /((1.0+CEBD/EMAX2) * (1.0+CEBD/EMAX2))           
          ELSE                                                                  
            RMAX2(IQX,J) = 0.0                                                  
          ENDIF                                                                 
        ENDIF                                                                   
   15 CONTINUE                                                                  
C                                                                               
C---- PRINT TABLE OF LAUNCH DATA.  DON'T PRINT BOTH SIDES IF IDENTICAL!         
C                                                                               
C     BRANCH AROUND THIS FOR NEUTRAL LAUNCH OPTIONS WHICH DO NOT INVOLVE
C     SPUTTERING
C 
      CALL PRB                                                                  
      IF ((CNEUTB.EQ.6).OR.(CNEUTB.EQ.7).OR.(CNEUTB.EQ.8)) GOTO 666
C
      IF (J.EQ.1) THEN                                                          
       CALL PRC ('SAMPLE PRIMARY FLUX AND YIELD DATA FOR Y < 0 SURFACE')        
       CALL PRC ('                                       *****        ')        
      ELSEIF (J.EQ.2 .AND. CTBOUL.EQ.CTBOUG .AND. CLTOUL.EQ.CLTOUG .AND.        
     >   CPRINT.EQ.0 .AND. CNBOUL.EQ.CNBOUG .AND. CLNOUL.EQ.CLNOUG .AND.        
     >   CYMFPS(0,1).EQ.CYMFPS(0,2)) THEN                                       
       CALL PRC ('SAMPLE PRIMARY FLUX AND YIELD DATA FOR Y > 0 SURFACE A        
     >S ABOVE')                                                                 
       CALL PRC ('                                       *****        ')        
       GOTO 666                                                                 
      ELSE                                                                      
       CALL PRC ('SAMPLE PRIMARY FLUX AND YIELD DATA FOR Y > 0 SURFACE')        
       CALL PRC ('                                       *****        ')        
      ENDIF                                                                     
      WRITE (7,9000)                                                            
      IF (CLARMR.GT.0.0) WRITE (7,9002)                                         
     >  QXS(1),CYMFPS(1,J),FLUX1(1,J),ENEGY1(1,J),                              
     >  YIELD1(1,J),FY1(1,J),FY1(1,J)*DELTAX/FYTOT1(J)                          
      DO 100 IQX = 0, -NQXSO/2, MIN (-1,-NQXSO/25)                              
        IF (YIELD1(IQX,J).GT.1.E-3) WRITE (7,9002)                              
     >    QXS(IQX),CYMFPS(IQX,J),FLUX1(IQX,J),ENEGY1(IQX,J),                    
     >    YIELD1(IQX,J),FY1(IQX,J),FY1(IQX,J)*DELTAX/FYTOT1(J)                  
  100 CONTINUE                                                                  
C                                                                               
      CALL PRB                                                                  
      CALL PRR ('TOTAL PRIMARY INTEGRATED FLUX*SINTB      ',FTOT1(J))           
      IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN      
        CALL PRR ('PRIMARY (HIGH) INTEGRATED FLUX*YPH*SINTB ',                  
     >   FYTOT1(J) / (1.0+CQPL/(CQ(MAT1,MATLIM)*QMULTP)))          
        CALL PRR ('PRIMARY (LOW)  INTEGRATED FLUX*YPL*SINTB ',                  
     >   FYTOT1(J) - FYTOT1(J) / (1.0+CQPL/(CQ(MAT1,MATLIM)*QMULTP)))         
      ENDIF                                                                     
      CALL PRR ('TOTAL PRIMARY INTEGRATED FLUX*YIELD*SINTB',FYTOT1(J))          
      IF (CNEUTD.EQ.2)                                                          
     >CALL PRR ('EXPECTED PROPORTION OF LAUNCHES          ',FRAC1(J))           
C                                                                               
C---- PRINT MAXIMUM LAUNCH VELOCITIES ...                                       
C---- VLAN IS A FUNCTION TO CALCULATE THIS, DEPENDING ON NEUT OPTIONS           
C                                                                               
      CALL PRB                                                                  
      CALL PRC ('LIMITING RANDOM NUMBERS AND MAXIMUM LAUNCH VELOCITIES          
     >(EV)')                                                                    
      WRITE (7,9001) 'JUST OUTBOARD   ', RMAX1(IQXOUT,J),                       
     >           VLAN (CNEUTC,RMAX1(IQXOUT,J))                                  
      WRITE (7,9001) 'NEAR X =-AW/4   ', RMAX1(IQXAW4,J),                       
     >           VLAN (CNEUTC,RMAX1(IQXAW4,J))                                  
      WRITE (7,9001) 'NEAR X =-AW/2   ', RMAX1(IQXAW2,J),                       
     >           VLAN (CNEUTC,RMAX1(IQXAW2,J))                                  
      WRITE (7,9001) 'NEAR X =-AW     ', RMAX1(IQXAW,J),                        
     >           VLAN (CNEUTC,RMAX1(IQXAW ,J))                                  
      IF (IQXFAC.LT.0)                                                          
     >WRITE (7,9001) 'NEAR X =-LAMBDA ', RMAX1(IQXFAC,J),                       
     >           VLAN (CNEUTC,RMAX1(IQXFAC,J))                                  
C                                                                               
  666 CONTINUE                                                                  
      GTOT1  = 0.5 * (FTOT1 (1) + FTOT1 (2))                                    
      GYTOT1 = 0.5 * (FYTOT1(1) + FYTOT1(2))                                    
C                                                                               
C-----------------------------------------------------------------------        
C   GENERATE LAUNCH POSITIONS FOR PRIMARY AND SECONDARY LAUNCHES                
C-----------------------------------------------------------------------        
C                                                                               
C---- "NPROD1" IS TOTAL NUMBER OF PRIMARY NEUTRALS PRODUCED                     
C---- "NPROD2" IS TOTAL NUMBER OF SECONDARY NEUTRALS PRODUCED                   
C---- WHICH WE COUNT BACKWARDS FROM NPROD, UNTIL THE TWO MEET                   
C---- SOMEWHERE IN THE MIDDLE.                                                  
C---- GENERATE VECTOR OF RANDOM NUMBERS IN RANGE 0 TO 1                         
C                                                                               
      NPROD1 = 0                                                                
      NPROD2 = NPROD + 1                                                        
C                                                                               
      CALL SURAND (SEED, NPROD, RANVA)                                          
      NRAND = NRAND + NPROD                                                     
      CALL SURAND (SEED, NPROD, RANVB)                                          
      NRAND = NRAND + NPROD                                                     
      CALL SURAND (SEED, NPROD, RANVC)                                          
      NRAND = NRAND + NPROD                                                     
      CALL SURAND (SEED, NPROD, RANVD)                                          
      NRAND = NRAND + NPROD                                                     
C                                                                               
C---- FOR EACH NEUTRAL, FIRST ASSUME IT IS PRIMARY AND GET RANDOM NO.           
C                                                                               
      DO 300 IPROD = 1, NPROD                                                   
        PRIME  = .TRUE.                                                         
        SPREAD = .FALSE.                                                        
        RAN    = RANVA (IPROD)                                                  
        IF (RAN.LE.FYTOT1(1)/(2.0*GYTOT1)) THEN                                 
          J = 1                                                                 
        ELSE                                                                    
          J = 2                                                                 
        ENDIF                                                                   
        RAN = RANVB (IPROD)                                                     
C                                                                               
C------ DEPENDING ON OPTION CHOSEN,  SELECT LAUNCH POSITION ON LIMITER          
C------ EDGE AND DETERMINE IF LAUNCH IS A SECONDARY LAUNCH ...                  
C------ NOTE: CNEUTB=3 CASE.  SET IQX=0 EVEN IF LAUNCH IS AT X0>0, SINCE        
C------ WE LATER REQUIRE RMAX1(IQX) WHICH IS UNDEFINED FOR IQX0>...             
C------ IF RMAX=0 TO PREVENT A FAILED LAUNCH SELECT A NEW X POSITION ..         
C------ FOR LARMOR RADIUS EFFECT, IQX=1 BECOMES A VALUE WITHIN {KQX,0}.         
C------ THE 1.E-10 CORRECTION IS NEEDED FOR SLAB LIMITER, SO THAT "J"           
C------ WILL BE ALTERNATELY 1 AND 2 IN LAUNCH.                                  
C                                                                               
        IF (CNEUTB.EQ.0.OR.CNEUTB.EQ.5) THEN                                 
c          write(6,*) 'before fycum ipos' 
  285     IQX = IPOS (RAN, FYCUM(ICUT(J),J), 1-ICUT(J)) + ICUT(J) - 1           
          IF (RAN.GT.FSPLIT(IQX,J)) PRIME = .FALSE.                             
          IF (((   PRIME  ).AND.(RMAX1(IQX,J).LE.0.0)) .OR.                     
     >        ((.NOT.PRIME).AND.(RMAX2(IQX,J).LE.0.0))) THEN                    
            CALL SURAND (SEED, 1, RAN)                                          
            GOTO 285                                                            
          ENDIF                                                                 
          IF (IQX.EQ.1) THEN                                                    
            NRAND  = NRAND + 1                                                  
            CALL SURAND (SEED, 1, RAN)                                          
            IQX    = NINT (REAL(KQX) * RAN)                                     
            SPREAD = .TRUE.                                                     
          ENDIF                                                                 
          X0  = QXS(IQX)                                                        
          IF (J.EQ.1) THEN                                                      
            Y0    =-QEDGES(IQX,1) - 1.E-10                                      
            THETA = PI - QTANS(IQX,1)                                           
          ELSE                                                                  
            Y0    = QEDGES(IQX,2) + 1.E-10                                      
            THETA = QTANS(IQX,2)                                                
          ENDIF                                                                 
C                                                                               
        ELSEIF (CNEUTB.EQ.1 .OR. CNEUTB.EQ.2) THEN                              
          IQX = INT (CXSC * XSCALO)                                             
          X0  = CXSC                                                            
          IF ((CNEUTB.EQ.1) .AND. J.EQ.1) THEN                                  
            Y0    =-QEDGES(IQX,1) - 1.E-10                                      
            THETA = PI - QTANS(IQX,1)                                           
          ELSE                                                                  
            Y0    = QEDGES(IQX,2) + 1.E-10                                      
            THETA = QTANS(IQX,2)                                                
          ENDIF                                                                 
          IF (RAN.GT.FRAC1(J)) PRIME = .FALSE.                                  
C                                                                               
        ELSEIF (CNEUTB.EQ.3) THEN                                               
          X0    = CXSC                                                          
          IQX   = 0                                                             
          Y0    = CYSC                                                  
          THETA = PI/2.0                                                        
          IF (RAN.GT.FRAC1(J)) PRIME = .FALSE.                                  
C
c slmod
        ELSEIF (CNEUTB.EQ.9) THEN
c
c 3D Gaussian launch about the point (X0,Y0,P0).
c
          X0    = CXSC                                                          
          IQX   = 0                                                             

250       CALL SURAND(SEED,1,RAN)
          NRAND = NRAND +1
          Y0 = 6.0E-2 * (2.0 * RAN - 1.0) 
          
          CALL SURAND(SEED,1,RAN)
          NRAND = NRAND + 1
          IF (RAN.GT.EXP(-((Y0/0.02)**2))) GOTO 250

          Y0 = Y0 + CYSC

252       CALL SURAND(SEED,1,RAN)
          NRAND = NRAND +1
          P0 = 6.0E-2 * (2.0 * RAN - 1.0)
          
          CALL SURAND(SEED,1,RAN)
          NRAND = NRAND + 1
          IF (RAN.GT.EXP(-((P0/0.02)**2))) GOTO 252

          P0 = P0 + CPSC

          THETA = PI/2.0                                                        
          IF (RAN.GT.FRAC1(J)) PRIME = .FALSE.                           

        ELSEIF (CNEUTB.EQ.10) THEN
c
c True 3D Gaussian launch about (X0,Y0,P0).
c
          X0    = CXSC                                                          
          IQX   = 0                                                             

255       CALL SURAND(SEED,1,RAN)
          NRAND = NRAND +1
          Y0 = 6.0E-2 * RAN 
          
          CALL SURAND(SEED,1,RAN)
          NRAND = NRAND + 1
          IF (RAN.GT.EXP(-((Y0/0.02)**2))) GOTO 255

          CALL SURAND(SEED,1,RAN)
          NRAND = NRAND + 1
          PHI = 2.0 * PI * RAN
          
          P0 = Y0 * SIN (PHI) + CPSC
          Y0 = Y0 * COS (PHI) + CYSC

          THETA = PI/2.0                                                        
          IF (RAN.GT.FRAC1(J)) PRIME = .FALSE.                           
c slmod end
C                                                                               
        ELSEIF (CNEUTB.EQ.4) THEN                                               
c
c         Wall launch - need to adjust for new wall options
c
          if (lim_wall_opt.eq.0) then  

             IQX = 1-NQXSO                                                         
             X0  = QXS(IQX)                                                        
  290        Y0  = RAN * CTWOL                                                     
             IF (J.EQ.1) Y0 = -Y0                                                  
             IF ((Y0.LE.QEDGES(IQX,2)-CTWOL) .OR.                                  
     >           (Y0.GE.-QEDGES(IQX,1).AND.Y0.LE.QEDGES(IQX,2)) .OR.               
     >           (Y0.GE.CTWOL-QEDGES(IQX,1))) THEN                                 
               NRAND = NRAND + 1                                                   
               CALL SURAND (SEED, 1, RAN)                                          
               GOTO 290                                                            
             ENDIF                                                                 
             THETA = PI/2.0                                                        

          elseif (lim_wall_opt.eq.1) then 


 291         Y0  = RAN * CTWOL                                                     

             IQY_TMP = max(min(int(y0*yscale)+1,nqys),1)

             X0 = caw_qys(iqy_tmp)

             IQX = INT (X0 * XSCALO)                                             

             IF (J.EQ.1) Y0 = -Y0                                                  
             IF ((Y0.LE.QEDGES(IQX,2)-CTWOL) .OR.                                  
     >           (Y0.GE.-QEDGES(IQX,1).AND.Y0.LE.QEDGES(IQX,2)) .OR.               
     >           (Y0.GE.CTWOL-QEDGES(IQX,1))) THEN                                 
               NRAND = NRAND + 1                                                   
               CALL SURAND (SEED, 1, RAN)                                          
               GOTO 291                                                            
             ENDIF                                                                 
c
c            Normal incidence needs to be corrected at some point for these wall options
c
             THETA = PI/2.0                                                        

          endif


        ELSEIF (CNEUTB.EQ.6) THEN
          X0 = CXSC
          IQX = 0
          NRAND = NRAND+1
          CALL SURAND (SEED,1,RAN)
          Y0 = CYSC*2.0*(RAN-0.5)
          THETA = 0.0
        ELSEIF (CNEUTB.EQ.7) THEN
          NRAND = NRAND +2
          CALL SURAND(SEED,1,RAN)
          IQX = IPOS(RAN,LPDCUM,CLPD-1)
          X0 = LPDION(IQX,1)
          CALL SURAND (SEED,1,RAN)
          Y0 = Y0S + RAN * (Y0L-Y0S)
          IQX = 0
          THETA = 0.0
        ELSEIF (CNEUTB.EQ.8) THEN
          NRAND = NRAND +2
          CALL SURAND(SEED,1,RAN)
          IQX = IPOS(RAN,LPDCUM,CLPD-1)
          Y0 = LPDION(IQX,1)
          CALL SURAND (SEED,1,RAN)
          X0 = X0S + RAN * (X0L-X0S)
          IQX = 0
          THETA = 0.0
        ENDIF                                                                   
C                                                                               
        IF (CNEUTC.EQ.11.OR.CNEUTB.EQ.5) THEN                             
          NRAND = NRAND + 1                                                     
          CALL SURAND (SEED, 1, RAN)                                            
          P0 = CPSC * (2.0 * RAN - 1.0)                                         
c slmod
        ELSEIF (CNEUTB.EQ.9.OR.CNEUTB.EQ.10) THEN
          P0 = P0
          WRITE(78,*) Y0,P0,X0,PHI*180.0/PI
c slmod end
        ELSE
          P0 = CPSC                                                             
        ENDIF                                                                   
C                                                                               
C------ NEUTRAL SPREADING OPTION ... CHANGE LAUNCH POINT TO ANYWHERE            
C------ WITHIN X0-DX/2 TO X0+DX/2, AND UPDATE CORRESPONDING EDGE POINT.         
C------ ALSO APPLIES TO LARMOR RADIUS EFFECT!                                   
C                                                                               
        IF ((CNEUTF.EQ.1) .OR. SPREAD) THEN                                     
          RAN   = RANVC(IPROD)                                                  
          OLDX0 = X0                                                            
          OLDY0 = Y0                                                            
          X0    = OLDX0 + (RAN-0.5) * DELTAX                                    
          CALL EDGINT (X0,IQX,J,Y0,D0)                                          
          IF (J.EQ.1) Y0 = -Y0                                                  
C         WRITE (6,'('' NEUT: X0,Y0,X0`,Y0`,T'',4F11.6,F9.2)')                  
C    >      OLDX0,OLDY0,X0,Y0,THETA*RADDEG                                      
        ENDIF                                                                   
C                                                                               
C------ STORE LAUNCH DETAILS IN PRODUCTION ARRAYS, COUNTING UP FOR              
C------ PRIMARIES AND DOWN FOR SECONDARIES.                                     
C                                                                               
        IF (PRIME) THEN                                                         
          NPROD1 = NPROD1 + 1                                                   
          XPRODS(NPROD1) = X0                                                   
          YPRODS(NPROD1) = Y0                                                   
          PPRODS(NPROD1) = P0                                                   
          RMAXS(NPROD1)  = RMAX1(IQX,J)                                         
          IF (MATLIM.EQ.13.OR.MATLIM.EQ.14
     >         .OR.MATLIM.EQ.15.OR.MATLIM.EQ.16
     >         .OR. MATLIM.EQ.17.OR.MATLIM.EQ.18) THEN
            RYIELD = 1.0
            RESPUT = .FALSE.
          ELSE
            RYIELD = YIELD1(IQX,J)/(1.0+CQPL/(CQ(MAT1,MATLIM)*QMULTP))       
            RESPUT = RES(MAT1,MATLIM,RYIELD,PRIME,CNEUTD,RANVD(IPROD),
     >                    QMULTP)        
          ENDIF 
          IF (RESPUT) RMAXS(NPROD1) = -1.0                                      
          SPUTYS(NPROD1) = 1.0                                                  
        ELSE                                                                    
          NPROD2 = NPROD2 - 1                                                   
          XPRODS(NPROD2) = X0                                                   
          YPRODS(NPROD2) = Y0                                                   
          PPRODS(NPROD2) = P0                                                   
          RMAXS(NPROD2)  = RMAX2(IQX,J)                                         
          IF (MATLIM.EQ.13.OR.MATLIM.EQ.14
     >        .OR.MATLIM.EQ.15.OR.MATLIM.EQ.16
     >        .OR.MATLIM.EQ.17.OR.MATLIM.EQ.18) THEN
            RYIELD = 1.0
            RESPUT = .FALSE.
          ELSE
            RYIELD = YIELD2(IQX,J)/(1.0+CQSL/(CQ(MAT2,MATLIM)*QMULTS))        
            RESPUT = RES(MAT2,MATLIM,RYIELD,PRIME,CNEUTD,RANVD(IPROD),
     >                    QMULTS)        
          ENDIF
          IF (RESPUT) RMAXS(NPROD2) = -1.0                                      
          SPUTYS(NPROD2) = 1.0                                                  
        ENDIF                                                                   
  300 CONTINUE                                                                  
      NPROD2 = NPROD - NPROD1                                                   
C
C     ALL OF THE REGULAR NEUTRALS HAVE BEEN LAUNCHED USING
C     THE SPECIFIED LAUNCH OPTION. THE FOLLOWING ADDITIONAL 
C     NEUTRAL PRIMARY LAUNCHES ARE TO SIMULATE CROSS-FIELD SPUTTERING.
C
C     TWO POSSIBLE METHODS ARE IN USE. THE FIRST CALCULATES AN
C     EXPECTED FLUX CONTRIBUTION FROM EACH Y-EROSION BIN BASED ON 
C     CONTINUITY (PARTICLE CONSERVATION) AN AN ARBITRARY STRENGTH FACTOR.
C     CFBGFF. THE SECOND RELEASES A FIXED NUMBER OF PARTICLES WITH A SPECIFIED 
C     DISTRIBUTION IN ADDITION TO THE ORIGINAL NUMBER SPECIFIED.
C
C     DAVID ELDER, NOV 22 1990 
C
      IF (CFBGFF.GT.0.0) THEN
C
C        CALCULATE THE NUMBER AND APPROXIMATE LAUNCH POSITION OF 
C        THE ADDITIONAL PRIMARIES. THEN AFTER THE TOTAL NUMBER IS KNOWN
C        MAKE ROOM FOR THEM AND THEN GENERATE THEIR PRECISE POSITIONS.
C
C        NJ = CFBGFF*NIMPS/(2*LP-) * NB(X)/NB(0) * DELTA Y
C
         IMPCF = 0 
         CFBGC = CFBGFF*NPROD/(CAP/(CWL*CLFACT))
         DO 400 IOY = 1,MAXOS
           IF (OYOUTS(IOY).GE.-OYMAX2(1).AND.OYOUTS(IOY).LE.OYMAX2(2))
     >       THEN   
             CALL YEDGINT(OYOUTS(IOY),X0,IQX,J,IERR)
C
C            WRITE(6,*) 'OYO,X0,IQX,J,IERR:',OYOUTS(IOY),X0,IQX,J,IERR 
C
             IF (IERR.EQ.0) THEN                   
               IMPNJ = CFBGC*QRNBS(IQX,J)/QRNBS(0,J)*OYWIDS(IOY)  
               NEROYS(IOY,6) = AINT(IMPNJ)
               NRAND = NRAND+1
               CALL SURAND(SEED,1,RAN)
               IF (RAN.LE.(IMPNJ-NEROYS(IOY,6))) 
     >           NEROYS(IOY,6) = NEROYS(IOY,6) + 1.0  
               IMPCF = IMPCF + INT(NEROYS(IOY,6))
             ENDIF
           ELSE
             NEROYS(IOY,6) = 0.0
           ENDIF 
400      CONTINUE  
C
C        MAKE SURE THAT THE TOTAL NUMBER IS LESS THAN THE MAXIMUM
C        IF IT IS NOT ISSUE A WARNING MESSAGE AND DON'T LAUNCH ANY 
C        ADDITIONAL
C
         IF ((NPROD+IMPCF).GT.MAXIMP) THEN 

            CALL PRB
            CALL PRI('NEUT ERROR: TOO MANY EXTRA CF PRIMARIES :'
     >        //' NONE LAUNCHED ',IMPCF)
            CALL PRB      
            WRITE(6,9100) IMPCF  
9100        FORMAT('NEUT ERROR: CREATING TOO MANY CF PRIM. ',G12.2)
         ELSE
C
         WRITE(6,*) 'IMPCF:' ,IMPCF
C
          DO 405 IPROD = 1,IMPCF
C
C         COPY THE FIRST BLOCK OF SECONDARIES TO THE END TO MAKE ROOM
C         FOR THE ADDITIONAL PRIMARIES
C
            XPRODS(NPROD+IPROD) = XPRODS(NPROD1+IPROD)
            YPRODS(NPROD+IPROD) = YPRODS(NPROD1+IPROD)
            PPRODS(NPROD+IPROD) = PPRODS(NPROD1+IPROD)
            RMAXS(NPROD+IPROD)  = RMAXS(NPROD1+IPROD)
            SPUTYS(NPROD+IPROD) = SPUTYS(NPROD1+IPROD)  

405      CONTINUE    

C
C         CALCULATE LAUNCH COORDINATES
C
       IPROD = 0
       DO 410 IOY = 1,MAXOS
          DO 420 IMPCNT = 1,INT(NEROYS(IOY,6))      
            IPROD = IPROD +1   

C            WRITE(6,*) 'IOY,NER:',IOY,NEROYS(IOY,6)
            
            ITER = 0
425         CONTINUE
C
C           DETERMINE A Y LAUNCH POSITION ON THE LIMITER
C           AT LEAST HALF OF THE BIN MUST BE VALID FOR IT TO HAVE 
C           ANY PARTICLES LAUNCHED FROM IT.
C
            NRAND = NRAND + 1
            CALL SURAND(SEED,1,RAN)
            Y0 = OYOUTS(IOY) + (RAN-0.5)*OYWIDS(IOY)               
            ITER = ITER + 1
            IF( ITER.GT.100) THEN 
              WRITE(6,*) 'ERROR POSITION NOT FOUND'
              STOP
            ENDIF
            IF (Y0.GT.OYMAX2(2).OR.Y0.LT.-OYMAX2(1)) GOTO 425  
C
C                  
            IF (CNEUTC.EQ.11.OR.CNEUTB.EQ.5) THEN                             
              NRAND = NRAND + 1                                               
              CALL SURAND (SEED, 1, RAN)                                      
              P0 = CPSC * (2.0 * RAN - 1.0)                                   
            ELSE                                                              
c slmod tmp
               WRITE(0,*) 'P0 set elsewhere: A'
c slmod end
               P0 = CPSC                                                       
            ENDIF                                                            
C
C           CALCULATE THE ACTUAL EDGE POSITION AND IQX BIN FROM THE 
C           GIVEN Y POSITION FOR LAUNCH.
C 
            CALL YEDGINT(Y0,X0,IQX,J,IERR)
            XPRODS(NPROD1+IPROD) = X0                                         
            YPRODS(NPROD1+IPROD) = Y0                                          
            PPRODS(NPROD1+IPROD) = P0                                       
            RMAXS(NPROD1+IPROD)  = RMAX1(IQX,J)                             
            IF (MATLIM.EQ.13.OR.MATLIM.EQ.14
     >          .OR.MATLIM.EQ.15.OR.MATLIM.EQ.16
     >          .OR.MATLIM.EQ.17.OR.MATLIM.EQ.18) THEN
              RYIELD = 1.0
              RESPUT = .FALSE.
            ELSE
              NRAND = NRAND + 1                                               
              CALL SURAND (SEED, 1, RAN)                                      
              RYIELD = YIELD1(IQX,J)/(1.0+CQPL/(CQ(MAT1,MATLIM)*QMULTP))    
              RESPUT = RES(MAT1,MATLIM,RYIELD,PRIME,CNEUTD,RAN,QMULTP)        
            ENDIF 
            IF (RESPUT) RMAXS(NPROD1+IPROD) = -1.0                           
            SPUTYS(NPROD1+IPROD) = 1.0                                       
420       CONTINUE
410     CONTINUE
        NPROD = NPROD + IMPCF
        NPROD1 = NPROD1 + IMPCF
C
        CALL PRB
        CALL PRC('ADDITIONAL CROSS-FIELD PRIMARY LAUNCH OPTION'
     >          //' : ON')
        CALL PRR('   CROSS-FIELD BACKGROUND FLUX FACTOR  = ',CFBGFF)
        CALL PRI('   TOTAL NUMBER OF ADDITIONAL PRIMARIES= ',IMPCF) 
        CALL PRB
       ENDIF 

      ENDIF
C
C     METHOD 2:
C
C     THE VARIABLE IMPADD CONTAINS THE NUMBER OF ADDITIONAL 
C     IMPURITIES TO BE LAUNCHED. THEY WILL BE LAUNCHED FROM THE 
C     LIMITER SURFACE UNIFORMLY ON A RANGE Y = +/- YCFADD.
C
C     DAVID ELDER, OCT 17 1990
C
      IF (IMPADD.GT.0) THEN
        DO 305 IPROD = 1,IMPADD
C
C       COPY THE FIRST BLOCK OF SECONDARIES TO THE END TO MAKE ROOM
C       FOR THE ADDITIONAL PRIMARIES
C
          XPRODS(NPROD+IPROD) = XPRODS(NPROD1+IPROD)
          YPRODS(NPROD+IPROD) = YPRODS(NPROD1+IPROD)
          PPRODS(NPROD+IPROD) = PPRODS(NPROD1+IPROD)
          RMAXS(NPROD+IPROD)  = RMAXS(NPROD1+IPROD)
          SPUTYS(NPROD+IPROD) = SPUTYS(NPROD1+IPROD)  
C
C         CALCULATE LAUNCH COORDINATES
C
          IF (CEXNEUT.EQ.0) THEN
C
C         FLAT - RECTANGULAR DISTRIBUTION OVER Y = +/- YCFADD
C       
307          CALL SURAND(SEED,1,RAN)
             NRAND = NRAND+1
             Y0 = 2.0 * YCFADD * (RAN - 0.5)     
             IF (Y0.GT.OYMAX2(2).OR.Y0.LT.-OYMAX2(1)) GOTO 307    
          ELSEIF (CEXNEUT.EQ.1) THEN 
C
C         NORMAL OPTION - YCFADD USED AS HALF-WIDTH   
C
309          CALL SURAND(SEED,1,RAN1)
             NRAND = NRAND+1
             IF (RAN1.EQ.0.0) GOTO 309       
             CALL SURAND(SEED,1,RAN2)
             NRAND = NRAND +1
             TEMP1 = SQRT(-2.0*LOG(RAN1))
             TEMP2 = COS(2.0*PI*RAN2)
             Y0 =YCFADD*TEMP1*TEMP2        
             IF (Y0.GT.OYMAX2(2).OR.Y0.LT.-OYMAX2(1)) GOTO 309    
C         WRITE(6,*) 'Y0,T1,T2,R1,R2 : ',Y0,TEMP1,TEMP2,RAN1,RAN2
          ENDIF  
C
          IF (CNEUTC.EQ.11.OR.CNEUTB.EQ.5) THEN                             
            NRAND = NRAND + 1                                               
            CALL SURAND (SEED, 1, RAN)                                      
            P0 = CPSC * (2.0 * RAN - 1.0)                                   
          ELSE                                                              
c slmod tmp
               WRITE(0,*) 'P0 set elsewhere: A'
c slmod end
            P0 = CPSC                                                       
          ENDIF                                                            
C
C         CALCULATE THE ACTUAL EDGE POSITION AND IQX BIN FROM THE 
C         GIVEN Y POSITION FOR LAUNCH.
C 
          CALL YEDGINT(Y0,X0,IQX,J,IERR) 
c          write(6,*) 'after yedgint ipos'
          IOY = IPOS(Y0,OYS,MAXOS-1)
          NEROYS(IOY,6) = NEROYS(IOY,6) + 1.0   
          XPRODS(NPROD1+IPROD) = X0                                         
          YPRODS(NPROD1+IPROD) = Y0                                          
          PPRODS(NPROD1+IPROD) = P0                                             
          RMAXS(NPROD1+IPROD)  = RMAX1(IQX,J)                                   
          IF (MATLIM.EQ.13.OR.MATLIM.EQ.14
     >        .OR.MATLIM.EQ.15.OR.MATLIM.EQ.16
     >        .OR.MATLIM.EQ.17.OR.MATLIM.EQ.18) THEN
            RYIELD = 1.0
            RESPUT = .FALSE.
          ELSE
            NRAND = NRAND + 1                                               
            CALL SURAND (SEED, 1, RAN)                                      
            RYIELD = YIELD1(IQX,J)/(1.0+CQPL/(CQ(MAT1,MATLIM)*QMULTP))         
            RESPUT = RES (MAT1,MATLIM,RYIELD,PRIME,CNEUTD,RAN,QMULTP)        
          ENDIF 
          IF (RESPUT) RMAXS(NPROD1+IPROD) = -1.0                              
          SPUTYS(NPROD1+IPROD) = 1.0                                         
305     CONTINUE
        NPROD = NPROD + IMPADD
        NPROD1 = NPROD1 + IMPADD
      ENDIF        
C
C     CALCULATE THE DEUTERIUM FLUX - PARALLEL AND CROSSFIELD FOR EACH 
C     EROSION BIN. STORE THESE AND THE TOTAL IN THE ARRAY CDFLUX. IF NO
C     CROSS-FIELD FLUXES HAVE BEEN SPECIFIED THESE ARRAY WILL BE ZERO.
C
C     THIS IS CALCULATED ONLY IF ASKED FOR I.E. CDCALC = 1
C
      CALL RZERO(CDFLUX,MAXOS*3)
      WRITE (6,*) 'NEUT:CDCALC : ',CDCALC
      IF (CDCALC.EQ.1) THEN 
        DO 310 IOY = 1,MAXOS
          CALL YEDGINT(OYOUTS(IOY),X0,IQX,J,IERR)
          IF (IERR.NE.0) THEN 
            OYCOORD(IOY,1) = CAW   
            OYIQX(IOY)  = 1-NQXSO    
          ELSE 
            OYCOORD(IOY,1) = X0         
            OYIQX(IOY) = IQX
          ENDIF
          CALL YEDGINT(OYS(IOY),X0,IQX,J,IERR)
          IF (IERR.EQ.1) THEN    
            OYCOORD(IOY,3) = CAW
          ELSE
            OYCOORD(IOY,3) = X0
          ENDIF
310     CONTINUE
C
C     CALCULATE DELTA X VALUES
C
        OYCOORD(1,2) = OYCOORD(1,3) - CAW
        OYCOORD(MAXOS,2) = OYCOORD(MAXOS,3) -CAW
        DO 312 IOY = 2,MAXOS-1
          IF (OYOUTS(IOY+1).GT.0.0) THEN 
            IF (OYOUTS(IOY).LE.0.0) THEN 
               OYCOORD(IOY,2) = -OYCOORD(IOY-1,3)
               X0 = 0.0 
            ELSE
               OYCOORD(IOY,2) = X0 -OYCOORD(IOY,3)
               X0 = OYCOORD(IOY,3)
            ENDIF
          ELSE
            OYCOORD(IOY,2) = OYCOORD(IOY,3) -OYCOORD(IOY-1,3)
          ENDIF
312     CONTINUE 
C
C     CALCULATE FLUXES
C
        NINIT = FLOAT(NPROD) - FLOAT(IMPADD) -FLOAT(IMPCF)
        WRITE(6,*) 'NINIT:',NINIT 
        J = 1
        DO 314 IOY = 1,MAXOS
          IF (OYOUTS(IOY).LT.-OYMAX2(1).OR.OYOUTS(IOY).GT.OYMAX2(2)) 
     >      THEN   
            CDFLUX(IOY,1) = 0.0
            CDFLUX(IOY,2) = 0.0
            CDFLUX(IOY,3) = 0.0
          ELSE
            IF (OYOUTS(IOY).GT.0.0) J=2
            CDFLUX(IOY,1) = FLUX1(OYIQX(IOY),J)*OYCOORD(IOY,2)
            IF (YIELD1(OYIQX(IOY),J).LE.0.0) THEN
              CDFLUX(IOY,2) = 1.0E25
            ELSE      
              CDFLUX(IOY,2) = NEROYS(IOY,6)/NINIT * 
     >                 FYTOT1(J)/YIELD1(OYIQX(IOY),J)
            ENDIF 
            CDFLUX(IOY,3) = CDFLUX(IOY,1)+CDFLUX(IOY,2)
          ENDIF
          WRITE(6,'(a,10(1x,g12.5))') 
     >      'IOY:FLUX1,OYCOORD2,CDF1:',IOY,
     >      FLUX1(OYIQX(IOY),J),OYCOORD(IOY,2),CDFLUX(IOY,1)
          WRITE(6,'(a,2i8,10(1x,g12.5))') 
     >      'IOY:J,YIELD1,NER,FYT,N,CDF2:',IOY,J,
     >      YIELD1(OYIQX(IOY),J),NEROYS(IOY,6),FYTOT(J),NINIT,
     >      CDFLUX(IOY,2)
314     CONTINUE 
        WRITE(6,*) 'NINIT:' ,NINIT
        DO 316 IOY = 1,MAXOS
          WRITE(6,'(i8,a,5(1x,g12.5))') 
     >             IOY,':OYC:',OYCOORD(IOY,1),OYCOORD(IOY,2),
     >              OYCOORD(IOY,3),OYIQX(IOY)
          WRITE(6,'(i8,a,5(1x,g12.5))') 
     >             IOY,':CDF:',CDFLUX(IOY,1),CDFLUX(IOY,2),
     >              CDFLUX(IOY,3)
316     CONTINUE
      ENDIF
C             
C                                                                               
C-----------------------------------------------------------------------        
C       LAUNCH AND FOLLOW PRIMARY NEUTRALS                                      
C-----------------------------------------------------------------------        
C                                                                               
      RATIZ1 = 0.0                                                              
      RNEUT1 = 0.0                                                              
      RWALL1 = 0.0                                                              
      RCENT1 = 0.0                                                              
      RTMAX1 = 0.0                                                              
      RSTRK1 = 0.0                                                              
      RFAIL1 = 0.0                                                              
      RRES1  = 0.0                                                              
      IF (NPROD1.GT.0)                                                          
     >CALL LAUNCH (FSRATE,1,NPROD1,1,NATIZ1,RSTRK1,RRES1,                       
     >             RATIZ1,RNEUT1,RWALL1,RCENT1,RTMAX1,SEED,NRAND,               
     >             NEUTIM,RFAIL1,STATUS,MAT1,MATLIM,QMULTP)                     

C                                                                               
C-----------------------------------------------------------------------        
C       STORE PRIMARY DETAILS IN "IONISATION -1" POSITION IN ARRAYS             
C-----------------------------------------------------------------------        
C                                                                               
C---- WHEN PLOTTING GRAPHS, THE "-1" POSITION WILL BE CALLED THE                
C---- "PRIMARY NEUTRAL" LINES, AND THE "0" POSITION WILL BE CALLED THE          
C---- "TOTAL NEUTRAL" LINES.  HENCE NOW THAT ALL PRIMARIES HAVE BEEN            
C---- LAUNCHED, BOTH POSITIONS SHOULD BE IDENTICAL.  FOR SUBSEQUENT             
C---- LAUNCHES OF SECONDARY, TERTIARY, ..ETC NEUTRALS ONLY THE "0"              
C---- POSITION WILL BE AFFECTED IN LAUNCH.                                      
C                                                                               
      DO 320 J = 1, 2                                                           
       DO 320 IX = 1, NXS                                                       
        NEROXS(IX,2,J) = NEROXS(IX,3,J)                                         
  320 CONTINUE                                                                  
      DO 325 IO = 1, MAXOS                                                      
        NEROYS(IO,2) = NEROYS(IO,3)                                             
        NERODS(IO,2) = NERODS(IO,3)                                             
        do ip=-maxnps,maxnps
           NERODS3(IO,IP,2) = NERODS3(IO,IP,3)
        end do
  325 CONTINUE                                                                  
C                                                                               
      DO 360 IX = 1, NXS                                                        
        DO 330 IY = -NYS, NYS                                                   
          DDLIMS(IX,IY,-1) = DDLIMS(IX,IY,0)                                    
          TIZS  (IX,IY,-1) = TIZS  (IX,IY,0)                                    
  330   CONTINUE                                                                
        DO 350 IY = -NY3D, NY3D                                                 
          IF (IY.EQ.0) GOTO 350                                                 
          DO 340 IP = -MAXNPS, MAXNPS                                           
            DDLIM3(IX,IY,-1,IP) = DDLIM3(IX,IY,0,IP)                            
            TIZ3  (IX,IY,-1,IP) = TIZ3  (IX,IY,0,IP)                            
            DO 335 IT = 1, NTS                                                  
              LIM5(IX,IY,-1,IP,IT) = LIM5(IX,IY,0,IP,IT)                        
  335       CONTINUE                                                            
  340     CONTINUE                                                              
  350   CONTINUE                                                                
  360 CONTINUE                                                                  
C                                                                               
      DO 370 IY = -NYS, NYS                                                     
        WALLS(IY,-1) = WALLS(IY,0)                                              
  370 CONTINUE                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C       SET UP FOR LAUNCH OF SECONDARY NEUTRALS                                 
C-----------------------------------------------------------------------        
C                                                                               
      RATIZ2 = 0.0                                                              
      RNEUT2 = 0.0                                                              
      RWALL2 = 0.0                                                              
      RCENT2 = 0.0                                                              
      RTMAX2 = 0.0                                                              
      RSTRK2 = 0.0                                                              
      RFAIL2 = 0.0                                                              
      RRES2  = 0.0                                                              
      IF (NPROD2.EQ.0) GOTO 999                                                 
      STATUS = 2                                                                
      CALL PRB                                                                  
      CALL PRC ('***  LAUNCHING  SECONDARY NEUTRALS  ***')                      
      CALL PRC ('      (USING SPUTTERING OPTION 2)')                            
C                                                                               
C                                                                               
C********** LOOP FOR EACH LIMITER SURFACE  1= Y < 0     2= Y > 0   *****        
C                                                                               
C                                                                               
      DO 777 J = 1, 2                                                           
      CALL PRB                                                                  
      IF (J.EQ.1) THEN                                                          
      CALL PRC('SAMPLE SECONDARY FLUX AND YIELD DATA FOR Y < 0 SURFACE')        
      CALL PRC('                                         *****        ')        
      ELSEIF (J.EQ.2 .AND. CTBOUL.EQ.CTBOUG .AND. CLTOUL.EQ.CLTOUG .AND.        
     >   CPRINT.EQ.0 .AND. CNBOUL.EQ.CNBOUG .AND. CLNOUL.EQ.CLNOUG .AND.        
     >   CYMFSS(0,1).EQ.CYMFSS(0,2)) THEN                                       
      CALL PRC('SAMPLE SECONDARY FLUX AND YIELD DATA FOR Y > 0 SURFACE A        
     >S ABOVE')                                                                 
      CALL PRC('                                         *****        ')        
      GOTO 777                                                                  
      ELSE                                                                      
      CALL PRC('SAMPLE SECONDARY FLUX AND YIELD DATA FOR Y > 0 SURFACE')        
      CALL PRC('                                         *****        ')        
      ENDIF                                                                     
      WRITE (7,9000)                                                            
      IF (CLARMR.GT.0.0) WRITE (7,9002)                                         
     >  QXS(1),CYMFSS(1,J),FLUX2(1,J),ENEGY2(1,J),                              
     >  YIELD2(1,J),FY2(1,J),FY2(1,J)*DELTAX/FYTOT2(J)                          
      DO 390 IQX = 0, -NQXSO/2, MIN (-1,-NQXSO/25)                              
        IF (YIELD2(IQX,J).GT.1.E-3) WRITE (7,9002)                              
     >    QXS(IQX),CYMFSS(IQX,J),FLUX2(IQX,J),ENEGY2(IQX,J),                    
     >    YIELD2(IQX,J),FY2(IQX,J),FY2(IQX,J)*DELTAX/FYTOT2(J)                  
  390 CONTINUE                                                                  
C                                                                               
      CALL PRB                                                                  
      CALL PRR ('TOTAL SECONDARY INTEGRATED FLUX*SINTB      ',FTOT2(J))         
      IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN                 
        CALL PRR ('SECNDRY (HIGH) INTEGRATED FLUX*YSH*SINTB ',                  
     >    FYTOT2(J) / (1.0+CQSL/(CQ(MAT2,MATLIM)*QMULTS)))         
        CALL PRR ('SECNDRY (LOW)  INTEGRATED FLUX*YSL*SINTB ',                  
     >    FYTOT2(J) - FYTOT2(J) / (1.0+CQSL/(CQ(MAT2,MATLIM)*QMULTS)))         
      ENDIF                                                                     
      CALL PRR ('TOTAL SECONDARY INTEGRATED FLUX*YIELD*SINTB',FYTOT2(J))        
      CALL PRR ('EXPECTED PROPORTION OF LAUNCHES            ',                  
     >                                                    1.0-FRAC1(J))         
      CALL PRB                                                                  
      CALL PRC ('LIMITING RANDOM NUMBERS AND MAXIMUM LAUNCH VELOCITIES          
     >(EV)')                                                                    
      WRITE (7,9001) 'JUST OUTBOARD   ', RMAX2(IQXOUT,J),                       
     >           VLAN (CNEUTC,RMAX2(IQXOUT,J))                                  
      WRITE (7,9001) 'NEAR X =-AW/4   ', RMAX2(IQXAW4,J),                       
     >           VLAN (CNEUTC,RMAX2(IQXAW4,J))                                  
      WRITE (7,9001) 'NEAR X =-AW/2   ', RMAX2(IQXAW2,J),                       
     >           VLAN (CNEUTC,RMAX2(IQXAW2,J))                                  
      WRITE (7,9001) 'NEAR X =-AW     ', RMAX2(IQXAW,J),                        
     >           VLAN (CNEUTC,RMAX2(IQXAW ,J))                                  
      IF (IQXFAC.LT.0)                                                          
     >WRITE (7,9001) 'NEAR X =-LAMBDA ', RMAX2(IQXFAC,J),                       
     >           VLAN (CNEUTC,RMAX2(IQXFAC,J))                                  
C                                                                               
  777 CONTINUE                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C       LAUNCH AND FOLLOW SECONDARY NEUTRALS                                    
C-----------------------------------------------------------------------        
C                                                                               
      CALL LAUNCH (FSRATE,NPROD1+1,NPROD2,NATIZ1+1,NATIZ2,RSTRK2,RRES2,         
     >             RATIZ2,RNEUT2,RWALL2,RCENT2,RTMAX2,SEED,NRAND,               
     >             NEUTIM,RFAIL2,STATUS,MAT2,MATLIM,QMULTS)               
      STATUS = 51                                               
C                                                                               
  999 CONTINUE                                                                  
      NATIZ = NATIZ1 + NATIZ2                                                   
      RATIZ = RATIZ1 + RATIZ2                                                   
      RNEUT = RNEUT1 + RNEUT2                                                   
      RWALLN= RWALL1 + RWALL2                                                   
      RCENT = RCENT1 + RCENT2                                                   
      RTMAX = RTMAX1 + RTMAX2                                                   
      RSTRUK= RSTRK1 + RSTRK2                                                   
      RFAIL = RFAIL1 + RFAIL2                                                   
      RRES  = RRES1  + RRES2                                                    
      WRITE (6,'('' NEUT: NPROD1 NPROD2 NPROD  NATIZ1 NATIZ2 NATIZ '',          
     >         /5X,6I7)') NPROD1,NPROD2,NPROD ,NATIZ1,NATIZ2,NATIZ              
C
C     CHANGE THE VEL/ANG FLAG BACK TO ITS ORIGINAL VALUE
C
      CNEUTC = TEMPOPT

C
      RETURN                                                                    
C                                                                               
 9000 FORMAT(3X,' X POSN   YMF    FLUX        ENERGY      YIELD    ',           
     >  ' FLUX*YIELD    FRACTION')                                              
 9001 FORMAT(1X,'    ',A16,1P,2(5X,G11.2))                                      
 9002 FORMAT(2X,F9.5,F7.3,1P,5(G11.2,1X))                                       
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LAUNCH (FSRATE,LPROD,NPROD,LATIZ,NATIZ,SSTRUK,SRES,            
     >                   SATIZ,SNEUT,SWALLN,SCENT,STMAX,SEED,NRAND,             
     >                   NEUTIM,SFAIL,STATUS,MAT,MATLIM,QMULT)              
      use variable_wall
      implicit none                                                    
      INTEGER    NPROD,NATIZ,NRAND,LPROD,LATIZ,STATUS,MAT,MATLIM                
      REAL       SATIZ,SNEUT,SWALLN,SCENT,STMAX,SSTRUK,NEUTIM,SFAIL             
      REAL       FSRATE,SRES,QMULT                                       
      DOUBLE PRECISION SEED                                                     
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  LAUNCH:  THIS ROUTINE CREATED TO IMPLEMENT NOTE 87 ON SELF -     *        
C  *  SPUTTERING.  FOLLOWS NEUTRALS FROM LAUNCH POSITIONS INPUT        *        
C  *  VIA THE ARGUMENT LIST UNTIL EVENTUAL IONISATION.  THE IONISATION *        
C  *  POSITIONS ARE RETURNED, TOGETHER WITH LAUNCH VELOCITIES, ETC.    *        
C  *    ALL NEUTRALS ARE FOLLOWED USING THE SAME TIMESTEP, FSRATE,     *        
C  *  NO MATTER WHAT THE X POSITION OF THE NEUTRAL HAPPENS TO BE.      *        
C  *  THUS ARRAY CPROBS HOLDS THE IONISING PROBABILITIES BASED ON      *        
C  *  FSRATE ALONE AND IS NOT AFFECTED BY TIMESTEP VARIATIONS.         *        
C  *                                                                   *        
C  *                                                                   *        
C  *  INPUT:  XPRODS ) LAUNCH COORDINATES FROM /CNEUT/                 *        
C  *          YPRODS )                                                 *        
C  *          PPRODS )                                                 *        
C  *          RMAXS    MAXIMUM RANDOM NUMBERS ALLOWED IN VIN CALCS.    *        
C  *                   USED IN SELF-SPUTTER CASES TO HOLD EIMPS /CNEUT/*        
C  *          SPUTYS   FRACTIONS OF ATOMS TO BE LAUNCHED - IN /CNEUT/  *        
C  *          FSRATE   TIMESTEP TO BE USED IN FOLLOWING NEUTRALS       *        
C  *          NPROD    NUMBER OF NEUTRAL FRAGMENTS TO START WITH       *        
C  *          SEED     RANDOM NUMBER GENERATOR SEED  (PASSED BACK TOO) *        
C  *          NRAND    COUNTS TOTAL RANDOMS USED     (PASSED BACK TOO) *        
C  *          NEUTIM   TIME USED TRACKING NEUTRALS   (PASSED BACK TOO) *        
C  *          MAT      BOMBARDING ION TYPE, USED FOR YIELDS            *        
C  *          MATLIM   LIMITER TYPE, USED FOR YIELDS                   *        
C  *                                                                   *        
C  *  OUTPUT: XATIZS ) FINAL POSITION COORDS AT IONISATION IN /CNEUT/  *        
C  *          YATIZS )                                                 *        
C  *          PATIZS )                                                 *        
C  *          VINS     LAUNCH VELOCITY W.R.T Y DIRECTION - IN /CNEUT/  *        
C  *          NATIZ    NUMBER OF IONISED NEUTRAL FRAGMENTS ON RETURN   *        
C  *          SSTRUK   TOTAL OF FRAGMENTS STRIKING LIMITER SURFACE     *        
C  *          SATIZ    TOTAL OF IONISED NEUTRAL FRAGMENTS              *        
C  *          SNEUT    TOTAL OF LAUNCHED NEUTRAL FRAGMENTS             *        
C  *          SWALLN   TOTAL OF FRAGMENTS PLATING OUT ON WALLS         *        
C  *          SCENT    TOTAL OF FRAGMENTS REACHING CENTRE              *        
C  *          STMAX    TOTAL OF FRAGMENTS EXISTING AT TMAX             *        
C  *          SFAIL    TOTAL OF FRAGMENTS WITH LAUNCH FAILURES         *        
C  *                     (TOTALLED OVER L AND R SIDES)                 *        
C  *                                      C.M.FARRELL   NOVEMBER 1987  *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE    'params'                                                       
C     INCLUDE    (PARAMS)                                                       
      INCLUDE    'dynam1'                                                       
C     INCLUDE    (DYNAM1)                                                       
      INCLUDE    'dynam3'                                                       
C     INCLUDE    (DYNAM3)                                                       
      INCLUDE    'comxyt'                                                       
C     INCLUDE    (COMXYT)                                                       
      INCLUDE    'comtor'                                                       
C     INCLUDE    (COMTOR)                                                       
      INCLUDE    'comtau'                                                       
C     INCLUDE    (COMTAU)                                                       
      INCLUDE    'comt2'                                                        
C     INCLUDE    (COMT2)                                                        
      INCLUDE    'cneut'                                                        
C     INCLUDE    (CNEUT)                                                        
      INCLUDE    'comnet'                                                       
C     INCLUDE    (COMNET)                                                       
      INCLUDE    'cyield'                                                       
C     INCLUDE    (CYIELD)                                                       
      INCLUDE    'crand'                                                        
C     INCLUDE    (CRAND)                                                        
c slmod begin
      INCLUDE    'slcom'
c slmod end      
C                                                                               
      REAL    XTOT(2),XMIN(2),XMAX(2),YTOT(2),AATIZ(2),PTOT(2),RRES(2)          
      REAL    ATOT(2),VTOT(2),VTOTM(2),VTOTA(2),VTOTAM(2),PATIZ(2)              
      REAL    VMULTT(2),VMULTM(2),ETOT(2),FRACIN(2),XATIZ(2),YATIZ(2)           
      REAL    VATIZ(2),VATIZM(2),TATIZ(2),EATIZ(2),REJECT(2),RATIZ(2)           
      REAL    RNEUT(2),RWALLN(2),RCENT(2),RTMAX(2),RSTRUK(2),RFAIL(2)           
C                                                                               
c      REAL    X,Y,ABSY,RMAX,SPUTY,PI,RADDEG                                     
      REAL    X,Y,ABSY,RMAX,SPUTY
      REAL    ANGLE,BETA,PSI,VMULT,TANLAN,TANGNT,ANGLAN,VIN                     
      REAL    TEMN,RAN,STATIM,ZA02AS,P,TANTRU,RPROD                             
      REAL    CISTOT,CISMAX,RSTMAX,TSTEPN,E,D                                   
      REAL    VY0,VY02,OLDY
      REAL    MAXXP(2)
      INTEGER IPROD,IQX,IQY,IX,IY,IFATE,IP,IPOS,IT,NREJEC,KK,KKLIM              
      integer :: iqy_tmp
      INTEGER JY,J,IOY,IOD                                                      
      CHARACTER FATE(6)*14                                                      
      DOUBLE PRECISION DSPUTY,DX,DY,DP,DXVELF,DYVELF,DPVELF,DWOL                
      LOGICAL RESPUT,FREEL                                                    

c
c     jdemod minimum rvalue to prevent division by zero errors
c
      real, parameter :: minrval = 1.0e-10
      
c slmod begin - N2 break
      LOGICAL N2STAT
      REAL    ENEW, VNEW,NFAST,N2Time
c slmod end
C                                                                               
c      DATA PI   /3.141592654/,     RADDEG / 57.29577952 /                       
      DATA FATE /'REACHED WALL',           'REACHED CENTRE',                    
     >           'TIME = TMAX',            'STRUCK LIMITER',                    
     >           'IONISED TO 1',           'FAILED LAUNCH'/                     
C                                                                               
C
C     IF THE NEUTRAL LAUNCH OCCURS IN FREESPACE THEN SET THE
C     VARIABLE FREEL TO TRUE SO THAT THE LIMITER COLLISION 
C     CHECKING WILL BE TURNED ON AND THE INITIAL LAUNCH CHECK
C     TURNED OFF.
C     D. ELDER , MAY 2 1990
C
c slmod begin - tmp
c      WRITE(0,*) 'Starting LAUNCH'

      IF (N2OPT.EQ.1) THEN
        CRMI = 2.0 * CRMI
        WRITE(0,*) 'Impurity mass set to ',CRMI,'for LAUNCH routine.'
      ENDIF
c slmod end
      IF ((CNEUTB.EQ.3).OR.(CNEUTB.EQ.6).OR.(CNEUTB.EQ.7).OR.
c slmod
     >     (CNEUTB.EQ.8).OR.(CNEUTB.EQ.9).OR.(CNEUTB.EQ.10)) THEN
c     >     (CNEUTB.EQ.8)) THEN
c slmod end
        FREEL = .TRUE.
      ELSE
        FREEL = .FALSE.
      ENDIF 
C
      DWOL = DBLE (CTWOL)                                                       
      STATIM = ZA02AS (1)                                                       
      DEBUGN = .FALSE.                                                          
      IF (CSTEPN.GT.0.0) THEN                                                   
        WRITE (6,9004) NINT(CSTEPN),FSRATE                                      
        DEBUGN = .TRUE.                                                         
      ENDIF                                                                   
      CISTOT = 0.0                                                              
      CISMAX = 0.0                                                              
      IQY    = 0                                                                
      iqy_tmp = 0
C                                                                               
C-----------------------------------------------------------------------        
C        SET INITIAL VALUES FOR DIAGNOSTICS                                     
C-----------------------------------------------------------------------        
C                                                                               
C---- XTOT  : TOTAL OF X POSITIONS AT PRODUCTION                                
C---- XMIN  : MINIMUM X FROM WHICH A NEUTRAL WAS LAUNCHED                       
C---- XMAX  : MAXIMUM    "     "        "      "     "                          
C---- MAXXP : MAXIMUM X AT WHICH A NEUTRAL WAS IONIZED  
C---- YTOT  : TOTAL OF ABS(Y) AT PRODUCTION                                     
C---- ATOT  : TOTAL OF ANGLES AT PRODUCTION, CALCULATED AS IF ALL               
C----         NEUTRALS WERE LAUNCHED FROM +Y REGION.                            
C---- VTOTA : TOTAL OF VELOCITIES AT PRODUCTION WITHOUT VMULT FACTORS           
C---- VTOTAM: MAXIMUM VELOCITY AT PRODUCTION WITHOUT VMULT FACTOR               
C---- VTOT  : TOTAL OF VELOCITIES AT PRODUCTION                                 
C---- VTOTM : MAXIMUM VELOCITY AT PRODUCTION                                    
C---- VMULTT: TOTAL OF VELOCITY ANGLE MULTIPLIERS "VMULT"                       
C---- VMULTM: MAXIMUM OF VELOCITY ANGLE MULTIPLIERS "VMULT"                     
C---- ETOT  : TOTAL OF TEMPERATURES AT PRODUCTION                               
C---- PTOT  : TOTAL OF ABS(P) POSITIONS AT PRODUCTION                           
C---- RWALLN: NO OF NEUTRALS REACHING X=-AW                                     
C---- RCENT : NO OF NEUTRALS REACHING X=A                                       
C---- RTMAX : NO OF NEUTRALS EXISTING AT T=TMAX                                 
C---- FRACIN: NO OF IONISATIONS THAT OCCUR INBOARD                              
C---- XATIZ : TOTAL OF X POSITIONS AT IONISATION                                
C---- YATIZ : TOTAL OF ABS(Y) AT IONISATION                                     
C---- AATIZ : TOTAL OF ANGLES AT IONISATION                                     
C---- VATIZ : TOTAL OF VELOCITIES AT IONISATION                                 
C---- VATIZM: MAXIMUM      "      "       "                                     
C---- EATIZ : TOTAL OF TEMPERATURES AT IONISATION                               
C---- TATIZ : TOTAL OF TIMES TO IONISATION                                      
C---- PATIZ : TOTAL OF ABS(P) AT IONISATION                                     
C---- REJECT: NUMBER OF VELOCITIES GREATER THAN MAX RANDOMS                     
C---- RNEUT : SUM OF FRAGMENTS LAUNCHED ) SET TO SMALL NO. TO PREVENT           
C---- RATIZ : SUM OF FRAGMENTS IONISED  ) POSSIBLE DIVIDE BY ZERO LATER         
C---- RSTRUK: NUMBER OF FRAGMENTS STRIKING LIMITER SURFACE                      
C---- RFAIL : NUMBER OF FAILED LAUNCHES  V > VMAX AT LEAST 1000                 
C----         TIMES, WHICH SUGGESTS THAT RMAX IS VERY NEAR ZERO!                
C---- RRES  : NO. OF NEUTRALS DUE TO RADIATION ENHANCED SUBLIMATION             
C                                                                               
      DO 50 J = 1, 2                                                            
        XTOT  (J) = 0.0                                                         
        XMIN  (J) = CA                                                          
        XMAX  (J) = CAW                                                         
        MAXXP (J) = CAW 
        YTOT  (J) = 0.0                                                         
        ATOT  (J) = 0.0                                                         
        VTOT  (J) = 0.0                                                         
        VTOTM (J) = 0.0                                                         
        VTOTA (J) = 0.0                                                         
        VTOTAM(J) = 0.0                                                         
        VMULTT(J) = 0.0                                                         
        VMULTM(J) = 0.0                                                         
        ETOT  (J) = 0.0                                                         
        PTOT  (J) = 0.0                                                         
        RWALLN(J) = 0.0                                                         
        RCENT (J) = 0.0                                                         
        RTMAX (J) = 0.0                                                         
        FRACIN(J) = 0.0                                                         
        XATIZ (J) = 0.0                                                         
        YATIZ (J) = 0.0                                                         
        AATIZ (J) = 0.0                                                         
        VATIZ (J) = 0.0                                                         
        VATIZM(J) = 0.0                                                         
        TATIZ (J) = 0.0                                                         
        EATIZ (J) = 0.0                                                         
        PATIZ (J) = 0.0                                                         
        REJECT(J) = 0.0                                                         
        ! jdemod - these fixed constants are too large for R*4 - should use the parameter LO=1.0e-37 for R4
        !RNEUT (J) = 1.0E-50                                                     
        !RATIZ (J) = 1.0E-50                                                     
        RNEUT (J) = LO                                                    
        RATIZ (J) = LO                                                     
        RSTRUK(J) = 0.0                                                         
        RFAIL (J) = 0.0                                                         
        RRES  (J) = 0.0                                                         
   50 CONTINUE                                                                  
C                                                                               
C---- SET  RSTMAX: MAX NUMBER OF ITERATIONS UP TO TIME = TMAX (0.1S)            
C---- (DIFFERS FROM CSTMAX WHICH APPLIES TO IONS,  BY FSRATE/QTIM AND           
C---- BECAUSE TMAX IS 0.1 HERE BUT 10 FOR IONS).                                
C---- SET NATIZ : NO OF NEUTRAL "FRAGMENTS" THAT GET IONISED                    
C---- SET RPROD : SUM OF FRAGMENTS PRODUCED                                     
c
c     NOTE: apply cutoff time for ions to neutrals as well - if it is less than 
c           0.1s - may need to change this if it causes problems later.  
C                                                                               
      if (ctmax.lt.0.1) then 
         RSTMAX = ctmax / FSRATE                                                     
      else
         RSTMAX = 0.1 / FSRATE                                                     
      endif  
      NATIZ  = 0                                                                
      RPROD  = 0.0                                                              
c slmod begin - N2 break
c
c Initialize the N2 break-up statistics:
c
      n2_t    = 0.0
      n_t     = 0.0
      n_n     = 0.0
      ni_n    = 0.0
      n2_x    = 0.0
      n2_y    = 0.0
      n2_p    = 0.0
      n_x     = 0.0
      n_y     = 0.0
      n_p     = 0.0
      n_e     = 0.0
      ni_e    = 0.0
      n_x_max = 0.0
      n_y_max = 0.0
      n_p_max = 0.0
      n_lim   = 0.0
      n2_lim  = 0.0
      n2_n    = 0.0
      n2_e    = 0.0
      n_beta  = 0.0
      n_psi   = 0.0
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C      LOOP FOR EACH NEUTRAL FRAGMENT TO BE LAUNCHED .....                      
C-----------------------------------------------------------------------        
C                                                                               
C---- SET START COORDINATES, POINTERS, ETC...  UPDATE EROSION DISTRIB'N         
C                                                                               
      CALL SURAND (SEED, NPROD, RANVA)                                          
      CALL SURAND (SEED, NPROD, RANVB)                                          
      CALL SURAND (SEED, NPROD, RANVC)                                          
      NRAND = NRAND + 3 * NPROD                                                 
      KK    = 1000 * ISECT                                                      
      KKLIM = KK - 10                                                           

C
C     SET UP VELOCITIES FOR USE IN THE NEUTRAL INJECTION CASE 6
C     - PARTICLES ARE INJECTED AT X0 IN +/-Y0 WITH ONE OF TWO
C       VELOCITIES
C
      IF (CNEUTC.EQ.13.OR.CNEUTC.EQ.16.OR.CNEUTC.EQ.17) THEN 
         VY0 = 1.38E4*SQRT(CENGSC/CRMI)
         VY02 = 1.38E4*SQRT(CEIN2/CRMI)
      ENDIF 
C
C                                                                               
      DO 900 IPROD = LPROD, LPROD+NPROD-1                                       
c slmod begin - N2 break
        IF (N2OPT.EQ.1) THEN
          N2STAT = .FALSE.
        ELSE
          N2STAT = .TRUE.
        ENDIF
c slmod end
        TSTEPN= CSTEPN                                                          
        IF (DEBUGN) WRITE (6,9005)                                              
        X     = XPRODS(IPROD)                                                   
        Y     = YPRODS(IPROD)                                                   
        P     = PPRODS(IPROD)                                                   
        DX    = DBLE (X)                                                        
        DY    = DBLE (Y)                                                        
        DP    = DBLE (P)                                                        
        ABSY  = ABS (Y)                                                         
        SPUTY = SPUTYS(IPROD)                                                   
        DSPUTY= DBLE (SPUTY)                                                    
        RPROD = RPROD + SPUTY                                                   
        J     = 1                                                               
        IF (Y.GE.0.0) J = 2                                                     
        IF (X.GE.0.0) THEN                                                      
          IQX = INT (X * XSCALI) + 1                                            
        ELSE                                                                    
          IQX = INT (X * XSCALO)                                                
        ENDIF                                                                   
        
c        write(0,'(a,i8,5g18.10)') 'IQX:',iqx,
c     >       INT (X * XSCALI) + 1,INT (X * XSCALO),x,xscalo,xscali


c       Calculate an IQY_TMP value to access CAW_QYS which gives the wall distance
c       at a specific value of QYS - this IQY_TMP has a different meaning than the 
c       IQY calculated below (which is relative to the limiter faces). 
c
        if (y.lt.0.0) then 
           IQY_TMP = max(min(int((y+ctwol)*yscale)+1,nqys),1)
        else
           IQY_TMP = max(min(int(y*yscale)+1,nqys),1)
        endif
c
        IX    = IPOS (X, XS, NXS-1)                                             
        IY    = IPOS (ABSY, YS, NYS-1)                                          
        IF (Y.LT.0.0) IY = -IY                                                  
        JY    = IABS (IY)                                                       
        IP    = IPOS (P, PS, 2*MAXNPS) - MAXNPS - 1                             
        CIST  = CTIMSC / FSRATE                                                 
        IT    = IPOS (CIST, CTIMES(1,0), NTS)                                   
c
c        write(0,'(a,4i7,3(1x,g12.5))') 'LAUNCH:',iprod,ix,iy,ip,x,y,p
C                                                                               

        RMAX  = RMAXS(IPROD)                                                    
        RESPUT= .FALSE.                                                         
        IF (RMAX.LT.0.0) THEN                                                   
          RMAX = 1.0                                                            
          RESPUT = .TRUE.                                                       
          RRES(J) = RRES(J) + SPUTY                                             
        ENDIF                                                                   
C                                                                               
        CALL EDGINT (X,IQX,J,E,D)                                               
        NEROXS(IX,3,J) = NEROXS(IX,3,J) + SPUTY                                 
        IF (J.EQ.1) THEN                                                        
          IOY = IPOS (-E, OYS, MAXOS-1)                                         
          IOD = IPOS (-D, ODS, MAXOS-1)                                         
        ELSE                                                                    
          IOY = IPOS (E, OYS, MAXOS-1)                                          
          IOD = IPOS (D, ODS, MAXOS-1)                                          
        ENDIF                                                                   
        NEROYS(IOY,3) = NEROYS(IOY,3) + SPUTY                                   
        NERODS(IOD,3) = NERODS(IOD,3) + SPUTY                                   
        NERODS3(IOD,IP,3) = NERODS3(IOD,IP,3) + SPUTY                                   
C                                                                               
C-----------------------------------------------------------------------        
C      SELECT ANGLE FOR LAUNCH ...                                              
C      PETER'S NOTES 38,65,83,93,109,218                                        
C-----------------------------------------------------------------------        
C                                                                               
C      ADD TRUE 3D FOR LIMITED EXTENT (TN3)
C
        IF     (RESPUT) THEN                                                    
          ANGLAN = SIGN (ASIN (SQRT(RANVA(IPROD))), RANVB(IPROD)-0.5)           
        ELSEIF (CNEUTC.EQ.13.OR.CNEUTC.EQ.15) THEN
          ANGLAN = SIGN(PI/2.0,RANVB(IPROD)-0.5)
        ELSEIF (CNEUTC.EQ.17) THEN
          ANGLAN = CIANGN
        ELSEIF (CNEUTC.EQ.0) THEN                                               
          ANGLAN = SIGN (ASIN (RANVA(IPROD)), RANVB(IPROD)-0.5)                 
        ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.5.OR.CNEUTC.EQ.9
     >          .OR.CNEUTC.EQ.16) THEN                 
          ANGLAN = SIGN (ASIN (SQRT(RANVA(IPROD))), RANVB(IPROD)-0.5)           
        ELSEIF (CNEUTC.EQ.1.OR.CNEUTC.EQ.2.OR.CNEUTC.EQ.4) THEN                 
          BETA   = ASIN (SQRT (RANVA(IPROD)))                                   
          PSI    = 2.0 * PI * RANVB(IPROD)                                      
          ANGLAN = ATAN (TAN (BETA) * COS (PSI))                                
        ELSEIF (CNEUTC.EQ.6) THEN                                               
          ANGLAN = 0.0                                                          
        ELSEIF (CNEUTC.EQ.7.OR.CNEUTC.EQ.11) THEN                               
          ANGLAN = SIGN (ACOS ((1.0-RANVA(IPROD)) ** (1.0/3.0)),                
     >                   RANVB(IPROD)-0.5)                                      
        ELSEIF (CNEUTC.EQ.8.OR.CNEUTC.EQ.14) THEN                        
          ANGLAN = 2.0 * PI * RANVA(IPROD) - PI                                 
        ELSEIF (CNEUTC.EQ.10) THEN                                              
          BETA   = ACOS ((1.0-RANVA(IPROD)) ** (1.0/3.0))                       
          PSI    = 2.0 * PI * RANVB(IPROD)                                      
          ANGLAN = BETA                                                         
          IF (PSI.LT.PI/2.0 .OR. PSI.GT.3.0*PI/2.0) ANGLAN = -BETA              
        ELSEIF (CNEUTC.EQ.12) THEN
          BETA = ASIN(SQRT(RANVA(IPROD)))
          PSI = 2.0*PI*RANVB(IPROD)
C         TRUE 3D - ANGLAN REPRESENTS PROJECTION ONTO 2D (NOT USED) 
          ANGLAN = ATAN(TAN(BETA)*COS(PSI))
        ENDIF                                                                   
C                                                                               
        IF ((CNEUTC.EQ.1.OR.CNEUTC.EQ.4).AND.(.NOT.RESPUT)) THEN                
          VMULT = SQRT (ABS(COS(BETA)**2+(SIN(BETA)**2*COS(PSI)**2)))           
        ELSE                                                                    
          VMULT = 1.0                                                           
        ENDIF                                                                   

C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE NORMAL TO LIMITER SURFACE AND ADD TO ANGLE                      
C     FOR SIDE PUFF GAS INJECTION, PUFF ALONG Y=0, IE RESET TANLAN TO 0         
C     THIS MAY RESULT IN SOME NEUTRALS BEING LAUNCHED STRAIGHT AT THE           
C     LIMITER SURFACE: COUNT THESE BUT DO NOT FOLLOW THEM.                      
C     VARIABLES ANGLAN & TANLAN REFER TO +Y REGION ONLY.  TRUE VALUES           
C     USED ARE ANGLE & TANGNT WHICH CORRECT FOR -Y OR +Y REGIONS.               
C-----------------------------------------------------------------------        
C                                                                               
c slmod

        if (iqx.ge.0) then 
           write(6,'(a,4i8,5(1x,g12.5))') 'IQX>=0:',
     >                iqx,ix,iy,iprod,x,y,sputy
           write(0,'(a,4i8,5(1x,g12.5))') 'IQX>=0:',
     >                iqx,ix,iy,iprod,x,y,sputy
        endif

        IF (CNEUTB.EQ.3 .OR. CNEUTB.EQ.4 .OR. CNEUTB.EQ.9 .OR.
     +      CNEUTB.EQ.10) THEN   
c        IF (CNEUTB.EQ.3 .OR. CNEUTB.EQ.4) THEN                                  
c slmod end
          TANTRU = PI / 2.0                                                     
        ELSEIF (FREEL) THEN
          TANTRU = 0.0 
        ELSEIF (Y.GE.0.0) THEN                                                  
          TANTRU = QTANS(IQX,2)                                                 
        ELSE                                                                    
          TANTRU = QTANS(IQX,1)                                                 
        ENDIF                                                                   
C                                                                               
C       TANLAN IS USED DIRECTLY IN VEL/ANG 12
C
        IF (CNEUTE.EQ.0) THEN                                                   
          TANLAN = TANTRU                                                       
        ELSE                                                                    
          TANLAN = CSNORM                                                       
        ENDIF                                                                   
C                                                                               
        IF (Y.GE.0.0) THEN                                                      
          ANGLE  = ANGLAN                                                       
          TANGNT = TANLAN                                                      
        ELSE                                   
          ANGLE  = -ANGLAN                                                      
          TANGNT = PI - TANLAN                                                  
        ENDIF                                                                   
C                                                                               
        XTOT(J)  = XTOT(J) + X * SPUTY                                          
        XMIN(J)  = MIN (XMIN(J), X)                                             
        XMAX(J)  = MAX (XMAX(J), X)                                             
        YTOT(J)  = YTOT(J) + ABSY * SPUTY                                       
        PTOT(J)  = PTOT(J) + ABS(P) * SPUTY                                     
        ATOT(J)  = ATOT(J) + (ANGLAN+TANLAN) * SPUTY                            
        RNEUT(J) = RNEUT(J) + SPUTY                                             

C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE LAUNCH VELOCITY.  FIRST SELECT RANDOM NUMBER IN THE             
C     RANGE 0 < RAN < RANMAX.  IF NOT FOUND, REJECT NUMBER AND FIND             
C     ANOTHER ONE ... TO PREVENT A POTENTIALLY INFINITE LOOP HERE, COUNT        
C     NUMBER OF REJECTED VELOCITIES IN NREJEC AND LIMIT TO (SAY) 1000;          
C     IF THIS LIMIT EXCEEDED THEN NEUTRAL IS CALLED A "FAILED LAUNCH"           
C-----------------------------------------------------------------------        
C                                                                               
        RAN = RANVC(IPROD)                                                      
        NREJEC = 0                                                              
  100   CONTINUE                                                                
        IF (RAN.GT.RMAX) THEN                                                   
          REJECT(J) = REJECT(J) + SPUTY                                         
          NREJEC = NREJEC + 1                                                   
          CALL SURAND (SEED, 1, RAN)                                            
          NRAND = NRAND + 1                                                     
          IF (NREJEC.LT.1000) GOTO 100                                          
          VIN   = 0.0                                                           
          TEMN  = 0.0                                                           
          RFAIL(J) = RFAIL(J) + SPUTY                                           
          IFATE = 6                                                             
          GOTO 899                                                              
        ENDIF                                                                   
C                                                                               
        IF     (RESPUT) THEN                                                    
          VIN = 1.38E4 * SQRT (0.15/CRMI)                                       
        ELSEIF (CNEUTC.EQ.13.OR.CNEUTC.EQ.16.OR.CNEUTC.EQ.17) THEN
          IF (RAN.GT.CPROB) THEN 
            VIN = VY02   
          ELSE
            VIN = VY0
          ENDIF 
        ELSEIF (CNEUTC.EQ.0.OR.CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.                   
     >          CNEUTC.EQ.5.OR.CNEUTC.EQ.12.OR.CNEUTC.EQ.14.OR.
     >          CNEUTC.EQ.15) THEN           
          VIN = 1.38E4 * SQRT (CEBD/(1.0/SQRT(RAN)-1.0)/CRMI) * VMULT           
        ELSEIF (CNEUTC.EQ.2) THEN                                               
          IF (RAN.GE.1.0) RAN = 0.999999                                        
          VIN = 1.38E4 * SQRT (CTGAS*ABS(LOG(1.0-RAN))/CRMI)                    
        ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.6.OR.CNEUTC.EQ.7.OR.                   
     >          CNEUTC.EQ.8.OR.CNEUTC.EQ.10.OR.CNEUTC.EQ.11) THEN               
          VIN = 1.38E4 * SQRT (CENGSC/CRMI)                                     
C                                                                               
C------ NOTE 156 VEL/ANG9 FLAG.  NOTICE THAT LAUNCHES ALTERNATELY AT            
C------ EIN1, EIN2,  BUT WE ALSO HAVE THE +/-Y THING.  HENCE HERE               
C------ WE NEED TO LAUNCH 2 PARTICLES AT EIN1, FOLLOWED BY 2 AT EIN2,           
C------ ETC TO ENSURE THAT SOME OF EACH ARE LAUNCHED ON EACH SIDE               
C------ OF Y = 0.                                                               
C                                                                               
        ELSEIF (CNEUTC.EQ.9) THEN                                               
          IF (2*(IPROD/4).EQ.IPROD/2) THEN                                      
            VIN = 1.38E4 * SQRT (CENGSC/CRMI)                                   
          ELSE                                                                  
            VIN = 1.38E4 * SQRT (CEIN2/CRMI)                                    
          ENDIF                                                                 
        ENDIF                                                                   
C                                                                               
        VTOTA(J)  = VTOTA(J) + VIN/VMULT * SPUTY                                
        VTOTAM(J) = MAX (VTOTAM(J), VIN/VMULT)                                  
        VMULTT(J) = VMULTT(J) + VMULT * SPUTY                                   
        VMULTM(J) = MAX (VMULTM(J), VMULT)                                      
        VTOT(J)   = VTOT(J) + VIN * SPUTY                                       
        VTOTM(J)  = MAX (VTOTM(J), VIN)                                         
C                                                                               
C------ CALCULATE X,Y COMPONENTS OF VELOCITY, PROB OF IONISATION, ETC           
C                                                                               
        IF (CNEUTC.EQ.13.OR.CNEUTC.EQ.15) THEN 
          DYVELF = VIN*FSRATE*SIN(ANGLE)
          DXVELF = 0.0
          DPVELF = 0.0
        ELSEIF (CNEUTC.EQ.17) THEN
          DYVELF = VIN*FSRATE*sin(angle)
          DXVELF = vin*fsrate*cos(angle)
          DPVELF = 0.0
        ELSEIF ((CNEUTC.EQ.10).AND.(.NOT.RESPUT)) THEN                      
          DXVELF = DBLE (VIN * COS(BETA) * FSRATE)                              
          DYVELF = DBLE (VIN * SIN(BETA) * COS(PSI) * FSRATE)                   
          DPVELF = DBLE (VIN * SIN(BETA) * SIN(PSI) * FSRATE)                   
        ELSEIF (CNEUTC.EQ.12) THEN 
          DXVELF = DBLE (VIN * FSRATE * (COS(BETA)*SIN(TANLAN)-SIN(BETA)
     >     * COS(PSI) * COS(TANLAN) ))
          DYVELF = DBLE (VIN * FSRATE * (COS(BETA) * COS(TANLAN) +
     >      SIN(BETA) * COS(PSI) * SIN(TANLAN)))
          DPVELF = DBLE( VIN * FSRATE * SIN(BETA) * SIN(PSI) )
C
C         IF LAUNCHING ON NEGATIVE Y SIDE REVERSE SIGN OF Y VELOCITY
C         ALL OTHERS ARE RANDOMLY SYMMETRIC AND DON'T NEED 
C         ADJUSTMENT
C         NOTE: WORK OUT A PROPER FORMULA FOR TANGNT IN THIS CASE AND
C               THIS TEST IS UNNECESSARY 
C
          IF (Y.LT.0) THEN
             DYVELF = -DYVELF
          ENDIF
          IF (DEBUGN) THEN
             WRITE(6,'(4(A10,G12.5))') 'V/AF12-VX:',DXVELF,
     >       '   VY:',DYVELF,'   VP:',DPVELF,' TANLAN:',TANLAN
          ENDIF 
c slmod begin
        ELSEIF (cneutc.eq.6) then
          DXVELF = DBLE (VIN * FSRATE)                      
          DYVELF = 0.0
          DPVELF = 0.0                                                         
c slmod end
        ELSE                      
          DXVELF = DBLE (VIN * SIN(ANGLE+TANGNT) * FSRATE)                      
          DYVELF = DBLE (VIN * COS(ANGLE+TANGNT) * FSRATE)                      
          DPVELF = 0.0                                                         
c slmod tmp
          WRITE(78,*) 'Inj vel:',DXVELF,DYVELF,DPVELF
c slmod end
          IF (DEBUGN) THEN
             WRITE(6,'(4(A10,G12.5))') 'V/AF  VX:',DXVELF,
     >       '   VY:',DYVELF,'   VP:',DPVELF,' TANLAN:',TANLAN
          ENDIF 
        ENDIF                                                                   

        TEMN   = CRMI * (VIN/1.38E4) * (VIN/1.38E4)                             
c slmod begin - N2 Break-up
        IF (n2opt.EQ.0) THEN
          ETOT(J)= ETOT(J) + TEMN * SPUTY                                         
        ELSE
          n2_e = n2_e + temn
          n2_v = n2_v + vin
        ENDIF
c   
c        ETOT(J)= ETOT(J) + TEMN * SPUTY                                         
c slmod end
        IF (DEBUGN) THEN                                                        
          WRITE (6,9003) IPROD,CIST,IQX,IQY,IX,IY,X,Y,VIN,TEMN,                 
     >      SPUTY,(ANGLE+TANGNT)*RADDEG,IP,P,IT,'NEUTRAL LAUNCH'                
        ENDIF                                                                   
C                                                                               
C------ CHECK IF NEUTRAL IS GOING TO STRIKE LIMITER SURFACE                     
C                                                                               
C------ CANNOT STRIKE IN NEUTRAL OPTION 6, 7 OR 8 I.E. FREE LAUNCH 
C
        IF ((.NOT.(FREEL))
     >    .AND.(((ANGLAN+TANLAN).LE.(TANTRU-PI/2.0))
     >    .OR. ((ANGLAN+TANLAN).GE.(TANTRU+PI/2.0)))) THEN                
          RSTRUK(J) = RSTRUK(J) + SPUTY                                         
          IFATE  = 4                                                            
          GOTO 899                                                              
        ENDIF                                                                   
C                                                                               
C-----------------------------------------------------------------------        
C  ITERATE AROUND MOTION FOR EACH TIMESTEP UNTIL IONISATION OCCURS.             
C  GENERATE NEW SET OF RANDOMS WHEN OLD LOT ARE ALL USED UP                     
C-----------------------------------------------------------------------        
C                                                                               

  200   CONTINUE                                                                
        IF (KK.GT.KKLIM) THEN                                                   
          CALL SURAND (SEED, KK, RANV)                                          
          NRAND = NRAND + KK                                                    
          KK = 0                                                                
        ENDIF                                                                   
C                                                                               
C------ UPDATE (X,Y,P) COORDINATES, IQX POINTER                                 
C------ CHECK Y VALUE ACCEPTABLE, IE  -2L < Y < 2L                              
C                                                                               
        OLDY = SNGL(DY)
        DX = DX + DXVELF                                                        
        DY = DY + DYVELF                                                        
        DP = DP + DPVELF                                                        
        IF (DY.GE.DWOL) THEN                                                    
          DY = DY - 2.0D0 * DWOL                                                
        ELSEIF (DY.LE.-DWOL) THEN                                               
          DY = DY + 2.0D0 * DWOL                                                
        ENDIF                                                                   
        X = SNGL (DX)                                                           
        Y = SNGL (DY)                                                           
        P = SNGL (DP)                                                           
        ABSY = ABS (Y)                                                          
        IF (X.GE.0.0) THEN                                                      
          IQX = INT (X * XSCALI) + 1                                            
        ELSE                                                                    
          IQX = INT (X * XSCALO)                                                
        ENDIF                                                                   
c
c       Update IQY_TMP
c
        if (y.lt.0.0) then 
           IQY_TMP = max(min(int((y+ctwol)*yscale)+1,nqys),1)
c           write(6,'(a,2i8,6g18.10)'),'IQY_TMP-:',iqy_tmp,
c     >               int((y+ctwol)*yscale)+1,
c     >               y,ctwol,yscale,y+ctwol,(y+ctwol)*yscale
        else
           IQY_TMP = max(min(int(y*yscale)+1,nqys),1)
c           write(6,'(a,2i8,6g18.10)'),'IQY_TMP+:',iqy_tmp,
c     >               int((y+ctwol)*yscale)+1,
c     >               y,ctwol,yscale,y*yscale
        endif
C                                                                               
C------ CHECK IF NEUTRAL HAS REACHED WALL, REACHED CENTRE,                      
C------ OR SURVIVED UNTIL TIME CUTOFF POINT ...                                 
C------ WALLS ARRAY IS UPDATED                                                  
C------ THE "IFATE" FLAG CAN BE USED WHEN PRINTING DEBUG MESSAGES WHICH         
C------ ALLOW THE NEUTRALS TO BE TRACKED.  (WRITE 9003 STATEMENTS)              
C                                                                               
        IF (X.LE.CAW_qys(iqy_tmp)) THEN                                                      

c
C         REFLECT NEUTRALS AT WALL IMPACT - IF OPTION ACTIVE
C
          IF (NRFOPT.EQ.1) THEN 
c
c           jdemod - this is only correct for a flat wall - if an alternate wall option
c                    is in use this should ideally be modified.  
c
            DXVELF = -DXVELF 
          ELSE
            RWALLN(J) = RWALLN(J) + SPUTY                                         
            WALLS(IY,0) = WALLS(IY,0) + SPUTY                                     
            IFATE = 1                                                             
            GOTO 899                                                              
          ENDIF  
        ELSEIF (X.GE.CA) THEN                                                   
          RCENT(J) = RCENT(J) + SPUTY                                           
          IFATE = 2                                                             
          GOTO 899                                                              
        ELSEIF (CIST.GT.RSTMAX) THEN                                            
          RTMAX(J) = RTMAX(J) + SPUTY                                           
          IFATE = 3                                                             
          GOTO 899                                                              
        ELSEIF ((X.LT.0.0).AND.(FREEL)) THEN
          IF  ( ( (OLDY.GT.0.0) .AND.
     >        ( (Y.LT.QEDGES(IQX,2)) .OR.
     >        (Y.GT. (CTWOL-QEDGES(IQX,1)) ) ) ) .OR.
     >        ( (OLDY.LE.0.0) .AND.
     >        ( (Y.LT. (-CTWOL+QEDGES(IQX,2)) ) .OR.
     >        (Y.GT. (-QEDGES(IQX,1)) ) ) ) ) THEN
C            WRITE(6,*) 'NEUT LIM:',X,Y,OLDY,QEDGES(IQX,2),
C     >                 -QEDGES(IQX,1),CTWOL,DWOL,IQX
            RSTRUK(J) = RSTRUK(J)  + SPUTY
c slmod begin - N2 break
            IF (N2OPT.EQ.1) THEN
              IF (N2STAT) THEN
                n_lim = n_lim + 1.0 
              ELSE
                n2_lim = n2_lim + 1.0
              ENDIF
            ENDIF
c slmod end
            IFATE = 4
            GOTO 899     
          ENDIF
        ENDIF                                                                   

C                                                                               
C------ SCORE PARTICLE IN DDLIMS ARRAY IN "IONISATION 0" POSITION.              
C------ CHECKING THE IP BIN IS NOT NEEDED FOR THE STANDARD LIM VERSION.         
C                                                                               
        IF (.NOT.BIG) GOTO 450                                                  
  430   CONTINUE                                                                
          IF ((IP.LE.-MAXNPS).OR.(PS(IP-1).LT.P))GOTO 440                       
          IP = IP - 1                                                           
          GOTO 430                                                              
  440   CONTINUE                                                                
          IF ((IP.GE.MAXNPS) .OR. (PS(IP).GE.P)) GOTO 450                       
          IP = IP + 1                                                           
          GOTO 440                                                              
  450   CONTINUE                                                                
          IF ((IX.LE.1) .OR. (XS(IX-1).LT.X))    GOTO 460                       
          IX = IX - 1                                                           
          GOTO 450                                                              
  460   CONTINUE                                                                
          IF ((IX.GE.NXS) .OR. (XS(IX).GE.X))    GOTO 470                       
          IX = IX + 1                                                           
          GOTO 460                                                              
  470   CONTINUE                                                                
          IF ((JY.LE.1) .OR. (YS(JY-1).LT.ABSY)) GOTO 480                       
          JY = JY - 1                                                           
          GOTO 470                                                              
  480   CONTINUE                                                                
          IF ((JY.GE.NYS) .OR. (YS(JY).GE.ABSY)) GOTO 490                       
          JY = JY + 1                                                           
          GOTO 480                                                              
  490   CONTINUE                                                                
        IY = JY                                                                 
        IF (Y.LT.0.0) IY = -IY                                                  
C                                                                               
C------ SCORE POSITION (IN D.P. ARRAYS DDLIMS/DDLIM3 FOR ACCURACY)              
C------ IF TIME POINT REACHED SCORE TIME POSITION ALSO, AND INCREMENT.          
C                                                                               

        DDLIMS(IX,IY,0) = DDLIMS(IX,IY,0) + DSPUTY                              
        IF (JY.LE.NY3D)                                                         
     >    DDLIM3(IX,IY,0,IP) = DDLIM3(IX,IY,0,IP) + DSPUTY                      
        IF (CIST.GE.CTIMES(IT,0)) THEN                                          
          IF (JY.LE.NY3D)                                                       
     >      LIM5(IX,IY,0,IP,IT) = LIM5(IX,IY,0,IP,IT) + SPUTY                   
          IT = IT + 1                                                           
        ENDIF                                                                   

C                                                                               
C------ SET NEW IONISATION PROBABILITY                                          
C------ CHECK FOR IONISATION: IF IT HAS NOT OCCURED JUMP BACK FOR               
C------ ANOTHER ITERATION.                                                      
C                                                                               
        CIST = CIST + 1.0                                                       
        IF (DEBUGN) THEN                                                        
          IF (CIST.GE.TSTEPN) THEN                                              
  495       TSTEPN = TSTEPN + CSTEPN                                            
            IF (TSTEPN.LE.CIST) GOTO 495                                        
            WRITE (6,9003) IPROD,CIST,IQX,IQY,IX,IY,X,Y,VIN,TEMN,               
     >        SPUTY,(ANGLE+TANGNT)*RADDEG,IP,P,IT,' '                           
          ENDIF                                                                 
        ENDIF                                                                   

C                                                                               
        KK = KK + 1                                                             
c slmod begin - N2 break
        IF (((     N2STAT).AND.(RANV(KK).GE.CPCHS  (IX,IY,0))).OR.
     +      ((.NOT.N2STAT).AND.(RANV(KK).GE.N2CPCHS(IX)     ))) GOTO 200
C                                                                               
C  IONISATION HAS OCCURED : STORE PARTICLE DETAILS IN ARRAYS / TOTALS           
C  ------------------------------------------------------------------           
C                                                                               
c 
c Check to see which path to molecular break-up this neutral has taken:
c
c 1) dissociation in to neutrals (N2 + e -> N + N + e)
c 2) ionisation dissociation (N2 + e -> N+ + N + 2e)
c 3) molecular ionisation (N2 + e -> N2+ + 2e ... + e -> N+ + N+ + 3e)
c
c Check to see if a high energy neutral is produced, and if it is then 
c reset neutral velocity and direction:
c
c *** OUTPUT N2 IONISATION DATA ***
c
c

         CALL SURAND (SEED, 1, RAN)                                          
         NRAND = NRAND + 1
c
c Check to see if the N2 event generated a fast neutral:
c
         IF (.NOT.N2STAT.AND.RAN.LT.(NNCPCHS(IX)/N2CPCHS(IX))) THEN
c
c Get a better ENEW: *** IMPORTANT ***
c
           IF (neopt.EQ.1) then
             enew = neng
             temn = enew
           ENDIF

           VNEW = 1.38E4 * SQRT (ENEW / (0.5 * CRMI))                                    

           CALL SURAND (SEED, 1, RAN)                                          
           NRAND = NRAND + 1
           beta = PI * RAN

           CALL SURAND (SEED, 1, RAN)                                          
           NRAND = NRAND + 1
           psi = 2.0 * PI * RAN
c
c *** ARGH! *** Was I making this horrible mistake all along?  Should it really be VNEW?
c
c *** ARGH! *** Add some good ****** statistics to catch this kind of thing!

           DXVELF = DBLE (VNEW * SIN(beta) * COS(psi) * FSRATE)                      
           DYVELF = DBLE (VNEW * SIN(beta) * SIN(psi) * FSRATE)                      
           DPVELF = DBLE (VNEW * COS(beta) * FSRATE)    
c
c Old:
c
c           DXVELF = DBLE (VIN * SIN(ANGLE) * FSRATE)                      
c           DYVELF = DBLE (VIN * COS(ANGLE) * FSRATE)                      
c           DPVELF = 0.0                                                         
c
c Set flag to indicate that a fast neutral is now being followed:
c
           N2STAT = .TRUE.

c           WRITE(63,'(A,6G10.3)') 'Fast neutral: ',
c     +       X,Y,VNEW,ANGLE * 180.0 / PI, DXVELF, DYVELF

c           WRITE( 0,'(A,6G10.3)') 'Fast neutral: ',
c     +       X,Y,VNEW,ANGLE * 180.0 / PI, DXVELF, DYVELF

c
c Statistics for fast neutral production:
c
c
c I think that some of the original neutral statistics (a little ways below) 
c are screwed up by this method - or maybe not since they only track the
c slow neutrals, not the fast ones.  Think on this...
c
           N2Time = CIST * FSRATE * SPUTY

           n2_t  = n2_t + n2time 
           n2_x  = n2_x + x
           n2_y  = n2_y + y
           n2_p  = n2_p + p
           n_e   = n_e  + enew
           n_v   = n_v  + vnew

           n_beta = n_beta + beta
           n_psi  = n_psi  + psi

           n_n   = n_n  + 1.0
           n2_n  = n2_n + 1.0

           GOTO 200
         ENDIF

c
c Record more statistics if a fast neutral was just ionized:
c
         IF (N2OPT.EQ.1) THEN
           IF (N2STAT) THEN
             n_t     = n_t + CIST * FSRATE * SPUTY - N2TIME
             n_x     = n_x + x
             n_y     = n_y + y
             n_p     = n_p + p
             n_x_max = MAX(ABS(x),n_x_max)
             n_y_max = MAX(ABS(y),n_y_max)
             n_p_max = MAX(ABS(p),n_p_max)
           ELSE
c
c A high energy primary (not ionisation of a fast neutral) ion was produced:
c
             IF (nieopt.EQ.1) THEN
               temn = nieng
             ENDIF

             ni_e = ni_e + temn
             n2_x = n2_x + x
             n2_y = n2_y + y
             n2_p = n2_p + p

             ni_n = ni_n + 1.0
             n2_n = n2_n + 1.0
           ENDIF

           N2STAT = .FALSE.
         ENDIF
c
c Original code:
c
c        IF (RANV(KK).GE.CPCHS(IX,IY,0)) GOTO 200                                
cC                                                                               
cC  IONISATION HAS OCCURED : STORE PARTICLE DETAILS IN ARRAYS / TOTALS           
cC  ------------------------------------------------------------------           
cC                                                                               

        IF (n2opt.EQ.1) THEN
          ETOT(J)= ETOT(J) + TEMN * SPUTY                                         
        ENDIF
c slmod end
         TIZS(IX,IY,0) = TIZS(IX,IY,0) + SPUTY                                  
         IF (JY.LE.NY3D)                                                        
     >     TIZ3(IX,IY,0,IP) = TIZ3(IX,IY,0,IP) + SPUTY                          
         NATIZ = NATIZ + 1                                                      
         XATIZS(LATIZ+NATIZ-1) = X                                              
         IF (X.GT.0.0) FRACIN(J) = FRACIN(J) + SPUTY                            
         YATIZS(LATIZ+NATIZ-1) = Y                                              
         PATIZS(LATIZ+NATIZ-1) = P                                              
C
C        THE Y-COMPONENT OF THE VELOCITY IS NEEDED AT THIS POINT TO 
C        PASS TO THE ION FOLLOWING ROUTINE. INSTEAD OF RECALCULATING 
C        THE Y-COMPONENT. SIMPLY EXTRACT IT FROM DYVELF BY DIVIDING 
C        BY THE TIME CONSTANT. THIS ALLOWS SPECIALIZED INITIAL 
C        CALCULATIONS FOR THE VELOCITY WITHOUT THE REQUIREMENT OF 
C        REPEATING THEM. SECOND, IT ALLOWS FOR THE POSSIBILITY THAT 
C        FUTURE CODE CHANGES COULD MODIFY THE Y-COMPONENT OF THE 
C        NEUTRAL VELOCITY. IN THIS CASE ONLY THE CURRENT VELOCITY 
C        OF THE NEUTRAL WOULD BE RETURNED TO LIM.
C
C        IF ((CNEUTC.EQ.10).OR.(CNEUTC.EQ.12)) THEN                          
C          VINS(LATIZ+NATIZ-1) = VIN * SIN (BETA) * COS (PSI)                   
C        ELSE                                                                   
C          VINS(LATIZ+NATIZ-1) = VIN * COS (ANGLE+TANGNT)                       
C        ENDIF                                                                  
C
         IF (FSRATE.GE.0.0) THEN 
            VINS(LATIZ+NATIZ-1) = SNGL(DYVELF)/FSRATE
         ELSE
            VINS(LATIZ+NATIZ-1) = 0.0
         ENDIF
C
         TATIZ(J) = TATIZ(J) + CIST * FSRATE * SPUTY                            
         XATIZ(J) = XATIZ(J) + X * SPUTY                                        
         YATIZ(J) = YATIZ(J) + ABSY * SPUTY                                     
         AATIZ(J) = AATIZ(J) + (ANGLAN+TANLAN) * SPUTY                          
         VATIZ(J) = VATIZ(J) + VIN * SPUTY                                      
         EATIZ(J) = EATIZ(J) + TEMN * SPUTY                                     
         PATIZ(J) = PATIZ(J) + ABS(P) * SPUTY                                   
         VATIZM(J)= MAX (VATIZM(J), VIN)                                        
         RATIZ(J) = RATIZ(J) + SPUTY                                            
         IFATE = 5                                                              
C                                                                               
  899  CONTINUE                                                                 

       MAXXP(J) = MAX(MAXXP(J),X) 
       IF (DEBUGN.OR.(IFATE.EQ.6))                                              
     >   WRITE (6,9003) IPROD,CIST,IQX,IQY,IX,IY,X,Y,VIN,TEMN,                  
     >   SPUTY,(ANGLE+TANGNT)*RADDEG,IP,P,IT,FATE(IFATE)                        
       CISTOT = CISTOT + CIST * SPUTY                                           
       CISMAX = MAX (CISMAX, CIST)                                              
  900  CONTINUE                                                                 


C                                                                               
C-----------------------------------------------------------------------        
C      PRINT DIAGNOSTICS                                                        
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (7,9001) NINT(CETH(MAT,MATLIM)),NINT(CETF(MAT,MATLIM)),             
     >                                    CQ(MAT,MATLIM)*QMULT         
      IF (CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7) THEN              
      CALL PRI2 ('NUMBER OF (HIGH) NEUTRALS LAUNCHED ',                         
     >           NINT(RNEUT(1)-RRES(1)),NINT(RNEUT(2)-RRES(2)))                 
      CALL PRI2 ('NUMBER OF (LOW)  NEUTRALS LAUNCHED ',                         
     >           NINT(RRES(1)),      NINT(RRES(2)))                             
      ELSE                                                                      
      CALL PRI2 ('NUMBER OF NEUTRALS LAUNCHED        ',                         
     >           NINT(RNEUT (1)),    NINT(RNEUT (2)))                           
      ENDIF                                                                     


      IF (STATUS.LE.2) THEN                                                     
      CALL PRI2 ('NUMBER OF NEUTRALS REACHING WALL   ',                         
     >           NINT(RWALLN(1)),    NINT(RWALLN(2)))                           
      CALL PRI2 ('NUMBER OF NEUTRALS REACHING CENTRE ',                         
     >           NINT(RCENT (1)),    NINT(RCENT (2)))                           
      CALL PRI2 ('NUMBER OF NEUTRALS EXISTING AT 0.1S',                         
     >           NINT(RTMAX (1)),    NINT(RTMAX (2)))                           
      CALL PRI2 ('NUMBER OF NEUTRALS STRIKING LIMITER',                         
     >           NINT(RSTRUK(1)),    NINT(RSTRUK(2)))                           
      CALL PRI2 ('NO FAILED LAUNCHES (1000 DISCARDS) ',                         
     >           NINT(RFAIL (1)),    NINT(RFAIL (2)))                           
      CALL PRI2 ('NO OF VELOCITIES > VMAX (DISCARDED)',                         
     >           NINT(REJECT(1)),    NINT(REJECT(2)))                           


      CALL PRR2 ('AVERAGE X PRODUCTION POSITION (M)  ',                         
     >           XTOT(1)/max(RNEUT(1),minrval),
     >           XTOT(2)/max(RNEUT(2),minrval))                          

      CALL PRR2 ('MINIMUM X PRODUCTION POSITION (M)  ',                         
     >           XMIN(1),            XMIN(2))                                   
      CALL PRR2 ('MAXIMUM X PRODUCTION POSITION (M)  ',                         
     >           XMAX(1),            XMAX(2))                                   
      CALL PRR2 ('AVERAGE ABS(Y) PRODN POSITION (M)  ',                         
     >           YTOT(1)/max(RNEUT(1),minrval),   
     >           YTOT(2)/max(RNEUT(2),minrval))                          

      CALL PRR2 ('AVERAGE ABS(P) PRODN POSITION (M)  ',                         
     >           PTOT(1)/max(RNEUT(1),minrval),
     >           PTOT(2)/max(RNEUT(2),minrval))                          
      CALL PRR2 ('AVERAGE ANGLE AT PRODUCTION (DEGS) ',                         
     >    RADDEG*ATOT(1)/max(RNEUT(1),minrval),
     >    RADDEG*ATOT(2)/max(RNEUT(2),minrval))                   

      CALL PRR2 ('AVERAGE VELOCITY WITHOUT ANGLE FACT',                         
     >           VTOTA(1)/max(RNEUT(1),minrval),
     >           VTOTA(2)/max(RNEUT(2),minrval))                         
      CALL PRR2 ('MAXIMUM VELOCITY WITHOUT ANGLE FACT',                         
     >           VTOTAM(1),          VTOTAM(2))                                 
      CALL PRR2 ('AVERAGE ANGULAR VEL CORRECTION FACT',                         
     >           VMULTT(1)/max(RNEUT(1),minrval),
     >           VMULTT(2)/max(RNEUT(2),minrval))                        

      CALL PRR2 ('MAXIMUM ANGULAR VEL CORRECTION FACT',                         
     >           VMULTM(1),          VMULTM(2))                                 
      CALL PRR2 ('AVERAGE VELOCITY AT PRODUCTION(M/S)',                         
     >           VTOT(1)/max(RNEUT(1),minrval),
     >           VTOT(2)/max(RNEUT(2),minrval))                          

      CALL PRR2 ('MAXIMUM VELOCITY AT PRODUCTION(M/S)',                         
     >           VTOTM(1),           VTOTM(2))                                  

      CALL PRR2 ('AVERAGE TEMP AT PRODUCTION (EV)    ',                         
     >           ETOT(1)/max(RNEUT(1),minrval),
     >           ETOT(2)/max(RNEUT(2),minrval))                          
      ENDIF                                                                     
C                                                                               
      SNEUT  = RNEUT (1) + RNEUT (2)                                            
      SWALLN = RWALLN(1) + RWALLN(2)                                            
      SCENT  = RCENT (1) + RCENT (2)                                            
      STMAX  = RTMAX (1) + RTMAX (2)                                            
      SSTRUK = RSTRUK(1) + RSTRUK(2)                                            
      SFAIL  = RFAIL (1) + RFAIL (2)                                            
      SATIZ  = RATIZ (1) + RATIZ (2)                                            
      SRES   = RRES  (1) + RRES  (2)                                            

      CTEMSC =(EATIZ (1) + EATIZ (2)) / max(SATIZ,minrval)                                   

      CALL PRB                                                                  
      CALL PRI2 ('NUMBER OF NEUTRALS IONISED         ',                         
     >           NINT(RATIZ(1)),     NINT(RATIZ(2)))                            
C                                                                               
      IF (SATIZ.GT.0.0 .AND. STATUS.LE.2) THEN                                  
         
        CALL PRR2 ('AVERAGE X IONISATION POSITION (M)  ',                       
     >             XATIZ(1)/max(RATIZ(1),minrval),
     >             XATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('MAXIMUM X IONIZATION POSITION (M)  ',
     >              MAXXP(1) , MAXXP(2) )                 
        CALL PRR2 ('AVERAGE ABS(Y) IONIS POSITION (M)  ',                       
     >             YATIZ(1)/max(RATIZ(1),minrval),
     >             YATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('AVERAGE ABS(P) IONIS POSITION (M)  ',                       
     >             PATIZ(1)/max(RATIZ(1),minrval),
     >             PATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('AVERAGE ANGLE AT IONISATION (DEGS) ',                       
     >      RADDEG*AATIZ(1)/max(RATIZ(1),minrval),
     >      RADDEG*AATIZ(2)/max(RATIZ(2),minrval))               
        CALL PRR2 ('AVERAGE VELOCITY AT IONISATION(M/S)',                       
     >             VATIZ(1)/max(RATIZ(1),minrval),
     >             VATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('MAXIMUM VELOCITY AT IONISATION(M/S)',                       
     >             VATIZM(1),           VATIZM(2))                              
        CALL PRR2 ('AVERAGE TEMP AT IONISATION (EV)    ',                       
     >             EATIZ(1)/max(RATIZ(1),minrval),
     >             EATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('AVERAGE TIME TO IONISATION (S)     ',                       
     >             TATIZ(1)/max(RATIZ(1),minrval),
     >             TATIZ(2)/max(RATIZ(2),minrval))                      
        CALL PRR2 ('FRACTION OF IONISATIONS INBOARD    ',                       
     >             FRACIN(1)/max(RATIZ(1),minrval),
     >             FRACIN(2)/max(RATIZ(2),minrval))                     
      ENDIF                                                                     
c slmod begin - N2 break
      IF (N2OPT.EQ.1) THEN
        WRITE(7,*) ' '
        WRITE(7,*) 'MOLECULAR NITROGEN BREAK-UP DIAGNOSTICS'        
        WRITE(7,*) ' '
 
        nfast = n_n - n_lim

        CALL PRI('Number of N2 events                 ', INT(n2_n))
        CALL PRR('  Mean X                            ', n2_x / n2_n)
        CALL PRR('  Mean Y                            ', n2_y / n2_n)
        CALL PRR('  Mean P                            ', n2_p / n2_n)
        CALL PRR('  Average N2 launch energy          ', n2_e / n2_n)
        CALL PRR('  Average N2 velocity               ', n2_v / n2_n)
        CALL PRR('  Fraction of events giving fast N  ', n_n / n2_n)
        CALL PRR('  Number of N2 striking limiter     ', n2_lim)

        IF (n_n.GT.0.0) THEN
          CALL PRI('Number of fast N produced           ', INT(n_n))
          CALL PRR('  Mean time to fast N production    ', n2_t / n_n)
          CALL PRR('  Mean energy                       ', n_e / n_n)
          CALL PRR('  Mean velocity                     ', n_v / n_n)
          CALL PRR('  Mean beta (1.571)                 ', n_beta / n_n)
          CALL PRR('  Mean psi  (3.142)                 ', n_psi / n_n)
          CALL PRR('  Mean lifetime of fast N           ', n_t / nfast)
          CALL PRR('  Mean X at ionisation              ', n_x / nfast)   
          CALL PRR('  Mean Y                            ', n_y / nfast)   
          CALL PRR('  Mean P                            ', n_p / nfast) 
          CALL PRR('  Maximum absolute X                ', n_x_max) 
          CALL PRR('  Maximum Y                         ', n_y_max) 
          CALL PRR('  Maximum P                         ', n_p_max) 
        ENDIF

        CALL PRR('  Number of fast N striking limiter ', n_lim)
        CALL PRI('Number of primary N+ produced       ', INT(ni_n))
        CALL PRR('  Mean energy                       ', ni_e / ni_n)
        CALL PRR('Average injection energy for N+     ',ctemsc)

        WRITE(7,*) ' '
      ENDIF
c slmod end
C                                                                               


      WRITE (6,'(1X,A,I15)') 'AV. ITERS PER NEUT ',NINT(CISTOT/RPROD)           
      WRITE (6,'(1X,A,I15)') 'MAX ITERS ANY NEUT ',NINT(CISMAX)                 
      WRITE (6,'(1X,A,I15)') 'TOTAL IONS CREATED ',NINT(SATIZ)                  
C                                                                               
      NEUTIM = NEUTIM + ZA02AS (1) - STATIM                                     
c slmod begin - tmp
      IF (N2OPT.EQ.1) THEN
        CRMI = 0.5 * CRMI
        WRITE(0,*) 'Impurity mass reset to ',CRMI
      ENDIF

c      WRITE(0,*) 'LAUNCH done'
c slmod end
      RETURN                                                                    
C                                                                               
 9001 FORMAT(/1X,'ETH,ETF,Q',I6,',',I7,',',F8.4,7X,'Y < 0     Y > 0')           
 9003 FORMAT(1X,I5,F8.1,I6,I6,I4,I4,' X',F9.6,' Y',F9.6,                        
     >  ' VIN',1P,G10.3,0P,' TEMP',F8.2,' FRAC',F5.2,                           
     >  ' ANG',F7.3,' P',I2,F7.4,I2,1X,A)                                       
 9004 FORMAT(//1X,'NEUT DEBUG: DIAGNOSTICS TO BE PRINTED EVERY',I6,             
     >  ' TIMESTEPS  (DELTA T =',G10.3,' SECONDS).',//)                         
 9005 FORMAT(1X,'-NEUT----TIME---IQX---IQY--IX--IY',72('-'),                    
     >  'IP-------IT',16('-'))                                                  
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      REAL FUNCTION VLAN (CNEUTC,RAN)                                           
      implicit none                                                    
      INTEGER CNEUTC                                                            
      REAL    RAN                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  VLAN:  THIS FUNCTION RETURNS THE MAXIMUM POSSIBLE NEUTRAL LAUNCH *        
C  *  VELOCITY FOR A GIVEN MAXIMUM RANDOM NUMBER "RAN" <= 1.0.         *        
C  *                                                                   *        
C  *                                      C.M.FARRELL   NOVEMBER 1987  *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      IF (CNEUTC.EQ.0.OR.CNEUTC.EQ.1.OR.CNEUTC.EQ.4.OR.CNEUTC.EQ.5
     >  .OR.CNEUTC.EQ.12.OR.CNEUTC.EQ.14.OR.CNEUTC.EQ.15) THEN       
        IF (RAN.LE.0.0) THEN                                                    
          VLAN = 0.0                                                            
        ELSE                                                                    
          IF (RAN.GE.1.0) RAN = 0.999999
          VLAN = 1.38E4 * SQRT (CEBD / (1.0/SQRT(RAN)-1.0) / CRMI)              
        ENDIF                                                                   
      ELSEIF (CNEUTC.EQ.2) THEN                                                 
        IF (RAN.GE.1.0) RAN = 0.999999                                          
        VLAN = 1.38E4 * SQRT (CTGAS * ABS(LOG(1.0-RAN)) / CRMI)                 
      ELSEIF (CNEUTC.EQ.3.OR.CNEUTC.EQ.6.OR.CNEUTC.EQ.7.OR.                     
     >        CNEUTC.EQ.8.OR.CNEUTC.EQ.10.OR.CNEUTC.EQ.11) THEN                 
        VLAN = 1.38E4 * SQRT (CENGSC / CRMI)                                    
      ELSEIF (CNEUTC.EQ.9.OR.CNEUTC.EQ.13.OR.CNEUTC.EQ.16.OR.
     >        CNEUTC.EQ.17) THEN
        VLAN = 1.38E4 * SQRT (MAX (CENGSC,CEIN2) / CRMI)                        
      ENDIF                                                                     
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
       LOGICAL FUNCTION RES (ION,LIMITR,RYIELD,PRIME,CNEUTD,RAN,QMULT)        
       implicit none                                                   
       REAL RYIELD,YL,RAN,YH,QMULT                                              
       INTEGER ION,LIMITR,CNEUTD                                                
       LOGICAL PRIME                                                            
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RES:   ADJUSTS YIELD IF NECESSARY FOR RADIATION ENHANCED         *        
C  *  SUBLIMATION CASE  (SPUTTER OPTION 5).                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'cyield'                                                          
C     INCLUDE (CYIELD)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      YH  = RYIELD                                                              
      YL  = 1.0                                                                 
      RES = .FALSE.                                                             
C
      IF (LIMITR.EQ.13.OR.LIMITR.EQ.14
     >     .OR.LIMITR.EQ.15.OR.LIMITR.EQ.16
     >     .OR.LIMITR.EQ.17.OR.LIMITR.EQ.18) RETURN
C
      IF ((CNEUTD.EQ.5.OR.CNEUTD.EQ.6.OR.CNEUTD.EQ.7)
     >    .AND.YH.GT.0.0) THEN                                       
        IF (PRIME) THEN                                                         
          YL = YH / (CQ(ION,LIMITR)*QMULT) * CQPL                           
        ELSE                                                                    
          YL = YH / (CQ(ION,LIMITR)*QMULT) * CQSL                           
          RYIELD = YH + YL                                                      
        ENDIF                                                                   
        IF ((CNEUTD.EQ.5.AND.RAN.LE.YL/(YH+YL)).OR.CNEUTD.EQ.7) THEN            
          IF (PRIME) RYIELD = YL                                                
          RES = .TRUE.                                                          
        ENDIF                                                                   
      ENDIF                                                                     
C     WRITE (6,9001) YH,YL,RAN,YL/(YH+YL),RYIELD                                
 9001 FORMAT(1X,'RES: YH',1P,E9.2,', YL',E9.2,', RAN',0P,F7.4,                  
     >  ', YL/(YH+YL)',F7.4,' RYIELD',1P,E9.2)                                  
      RETURN                                                                    
      END                                                                       
c
c
c
      SUBROUTINE SYIELD_set_mat2 (MAT2,CNEUTD,CBOMBF,CBOMBZ)
      IMPLICIT NONE
      INTEGER MAT2,CNEUTD,CBOMBF,CBOMBZ                   
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  SYIELD_SET_MAT2:  SETS THE VALUE OF MAT2 FOR SPUT OPT 2          *
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
      CHARACTER*6  BACMAT(7)                                                    
c                                                                                
      DATA BACMAT/                                                              
     &  ' H    ',' D    ',' T    ',' HE4  ',' C    ',' SELF ',' O    '/         
c
c     Set MAT2 
c
      MAT2 = CBOMBF                                                             
      IF (MAT2.EQ.0) MAT2 = 6                                                   
c
      IF (CNEUTD.EQ.2) THEN                                                     
        call prc ('SECOND BOMBARDING IONS '// BACMAT(MAT2))
        CALL PRI ('         WITH ZIMP', CBOMBZ)                                 
      ENDIF                                                                     
c
      RETURN                                                                    
      END                                                                       

