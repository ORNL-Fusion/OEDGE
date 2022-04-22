      SUBROUTINE PLASMA (NTBS,NTIBS,NNBS,CIOPTG,CIOPTK,QTIM)
      use mod_params
      use mod_comt2
      use mod_comtor
      use mod_comxyt
      use mod_slcom
      use yreflection
      use mod_soledge
      use allocatable_input_data
c     use mod_sol22_input_lim
c      use mod_sol22_lim
      IMPLICIT  none
      INTEGER   NTBS,NTIBS,NNBS,CIOPTG,CIOPTK 
      REAL      QTIM                            
c     INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c     INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
c     INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
c     INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c     slmod begin
c     INCLUDE   'slcom'
c     slmod end      
C     
C     *********************************************************************        
C     *                                                                   *        
C     *  PLASMA                                                           *        
C     *  ------                                                           *        
C     *  ROUTINE TO CALCULATE THE TEMPERATURE AND ION DENSITY OF THE      *        
C     *  BACKGROUND PLASMA AT A SET OF POINTS GIVEN BY QXS.  VARIOUS      *        
C     *  OPTIONS ARE ALLOWED TO GIVE COMBINATIONS OF EXPONENTIAL AND      *        
C     *  LINEAR DECAY.                                                    *        
C     *                                                                   *        
C     *  AMENDMENTS 16/8/88 ... ARRAYS CTEMBS, CRNBS NOW DEPENDENT ON     *        
C     *  X AND Y, EVERYTHING HAS TO BE CHANGED.  RETAIN QTEMBS AND QRNBS  *        
C     *  DIMENSIONED FOR QX AND RIGHT & LEFT SIDES OF LIMITER JUST FOR    *        
C     *  USE IN NEUT WHEN CALCULATING LAUNCH DISTRIBUTION FUNCTIONS.      *        
C     *                                                                   *
C     *  AMENDED 06/02/90 ... ARRAY CTEMBSI - CONTAINS THE ION TEMPERATURE*
C     *  DISTRIBUTION AS A FUNCTION OF X,Y. THIS IS THE FIRST STEP IN     *
C     *  CREATING COMPLETELY INDEPENDENT DIST'NS FOR IONS AND ELECTRONS.  *
C     *  USED FOR CALCULATING A NUMBER OF THE TAU VALUES (EXCEPT TIZ)     *
C     *  DAVID ELDER , FEB 6 , 1990      (08/02/90 - ALSO QTEMBSI)        *
C     *                                                                   *
C     *  AMENDED 25/09/90 ... ARRAYS CTIGS,CTEGS CONTAINING THE Y AXIS    *
C     *  TEMPERATURE GRADIENTS. THE TEMPERATURES ARE MODIFIED BY GIVEN    *
C     *  INPUT FUNCTIONS AND THE GRADIENT FOR EACH BIN IS CALCULATED      *
C     *  FROM THE RESULTING TEMPERATURES.                                 *
C     *  DAVID ELDER , SEPT 25 , 1990                                     *
C     *                                                                   *        
C     *  ARGUMENTS (ARRAYS FROM COMMONS) :-                               *        
C     *  QXS    : SET OF X POSITIONS TEMP AND DENSITY TO BE CALCULATED AT *        
C     *  CTEMBS : ARRAY TEMPERATURES TO BE RETURNED IN X,Y                *        
C     *  CTEMBSI: ARRAY OF ION TEMPERATURES INDEXED BY IX,IY              *
C     *  CRNBS  : ARRAY DENSITIES TO BE RETURNED IN X,Y                   *        
C     *  QTEMBS : ARRAY TEMPERATURES TO BE RETURNED IN QX,RIGHT & LEFT    *        
C     *  QTEMBSI: ARRAY OF ION TEMPERATURES INDEXED BY QX,RIGHT & LEFT    *        
C     *  QRNBS  : ARRAY DENSITIES TO BE RETURNED IN QX, RIGHT & LEFT      *        
C     *  NQXSO  : NUMBER OF POINTS (1-NQXSO:0) OUTBOARD                   *        
C     *  NQXSI  : NUMBER OF POINTS (1:NQXSI)   INBOARD                    *        
C     *  CA     : POSITION OF CENTRE "A"                                  *        
C     *  ELECTRON TEMPERATURE VARIABLES                                   *
C     *  CTBINS : SET OF (X,TB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C     *  CTBIN  : INBOARD  TEMPERATURE OF PLASMA BASE VALUE               *        
C     *  CLTIN1 : INBOARD  TEMP EXPONENTIAL DECAY   FOR   0   < X < CATIN *        
C     *  CLTIN2 : INBOARD  TEMP EXPONENTIAL DECAY   FOR CATIN < X <  CA   *        
C     *  CGTIN1 : INBOARD  TEMP LINEAR DECAY FACTOR FOR   0   < X < CATIN *        
C     *  CGTIN2 : INBOARD  TEMP LINEAR DECAY FACTOR FOR CATIN < X <  CA   *        
C     *  CATIN  : INBOARD  TEMP CROSSOVER POINT                           *        
C     *  CTBOUL : OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y < 0         *        
C     *  CTBOUG : OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y > 0         *        
C     *  NTEG   : NUMBER OF DIFFERENT Te MULTIPLIERS/REGIONS              *  
C     *  TMEG   : ARRAY SPECIFYING THE DIFFERENT REGIONS AND MULTIPLIERS  *
C     *  CTEGS  : CALCULATED Te TEMPERATURE GRADIENT FOR EACH BIN         *
C     *  ION TEMPERATURE PARAMETERS                                       *
C     *  CTIBINS: SET OF (X,TB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C     *  CTIBIN : INBOARD  TEMPERATURE OF PLASMA BASE VALUE               *        
C     *  CLTIIN1: INBOARD  TEMP EXPONENTIAL DECAY   FOR   0   < X < CATIIN*        
C     *  CLTIIN2: INBOARD  TEMP EXPONENTIAL DECAY   FOR CATIN < X <  CA   *        
C     *  CGTIIN1: INBOARD  TEMP LINEAR DECAY FACTOR FOR   0   < X < CATIIN*        
C     *  CGTIIN2: INBOARD  TEMP LINEAR DECAY FACTOR FOR CATIN < X <  CA   *        
C     *  CATIIN : INBOARD  TEMP CROSSOVER POINT                           *        
C     *  CTIBOUL: OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y < 0         *        
C     *  CTIBOUG: OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y > 0         *        
C     *  NTIG   : NUMBER OF DIFFERENT Ti MULTIPLIERS/REGIONS              *  
C     *  TMIG   : ARRAY SPECIFYING THE DIFFERENT REGIONS AND MULTIPLIERS  *
C     *  CTIGS  : CALCULATED Ti TEMPERATURE GRADIENT FOR EACH BIN         *
C     *  DENSITY PARAMETERS                                               *
C     *  CNBINS : SET OF (X,NB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C     *  CNBIN  : INBOARD  ION DENSITY OF PLASMA BASE VALUE               *        
C     *  CLNIN1 : INBOARD  DENS EXPONENTIAL DECAY   FOR   0   < X < CANIN *        
C     *  CLNIN2 : INBOARD  DENS EXPONENTIAL DECAY   FOR CANIN < X <  CA   *        
C     *  CGNIN1 : INBOARD  DENS LINEAR DECAY FACTOR FOR   0   < X < CANIN *        
C     *  CGNIN2 : INBOARD  DENS LINEAR DECAY FACTOR FOR CANIN < X <  CA   *        
C     *  CANIN  : INBOARD  DENSITY CROSSOVER POINT                        *        
C     *  CNBOUL : OUTBOARD ION DENSITY OF PLASMA BASE VALUE Y < 0         *        
C     *  CNBOUG : OUTBOARD ION DENSITY OF PLASMA BASE VALUE Y > 0         *        
C     *  NTBS   : NUMBER OF PAIRED VALUES IN CTBINS ARRAY                 *        
C     *  NTIBS  : NUMBER OF PAIRED VALUES IN CTIBINS ARRAY                *        
C     *  NNBS   : NUMBER OF PAIRED VALUES IN CNBINS ARRAY                 *        
C     *                                                                   *        
C     *  OPTION    OUTBOARD      INBOARD (1)          INBOARD (2)         *        
C     *  ------------------------------------------------------------     *        
C     *   0     CONSTANT          CONSTANT             CONSTANT           *        
C     *   1     EXPONENTIAL       LINEAR               LINEAR             *        
C     *   2     EXPONENTIAL       EXPONENTIAL          EXPONENTIAL        *        
C     *   3     EXPONENTIAL       FITTED TO SET OF GIVEN VALUES           *        
C     *   4     EXPONENTIAL       LINEAR               EXPONENTIAL        *        
C     *   5     EXPONENTIAL       EXPONENTIAL          LINEAR             *        
C     *   6     EXPONENTIAL       STANDARD JET  (NOTE 274)                *        
C     *                                                                   *        
C     *  CHRIS FARRELL  (HUNTERSKIL)  APRIL 1988                          *        
C     *  DAVID ELDER                  OCTOBER 1990
c
c     jdemod - the plasma calculated here is applied as the default
c     background for all poloidal plasma slices. It is too much
c     work at this time to replicate the functionality for multiple poloidal
c     slices. Instead the SOL22 code can be applied to each poloidal plasma
c     slice as required using the SOL22 inputs.      
c      
C     *                                                                   *        
C     *********************************************************************        
C     
      INTEGER IPOS,IXAT,IXAN,IXOUT,IX,IY,IQX,J                                  
      INTEGER COUNTI,COUNTJ
      INTEGER NTEMP,IN
      REAL TEGOUL,TEGOUG,TIGOUL,TIGOUG
      REAL TGPOS(2*MAXINS+2),TGMULT(2*MAXINS+1)
      REAL LEN1,LEN2,MULT,GRAD1,GRAD2,DSTEP
      REAL TGSCAL
c
      integer :: pz
C     
c     Apply plasma solution to default plasma poloidal zone and copy to
c     all other zones at the end of the routine.       
c      
      pz = 1

      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
C     
C     IF CIOPTK = -1 THEN USE Ti = Te
C     DO THIS BY SETTING ALL OF THE ION CONSTANTS EQUAL TO THOSE 
C     ENTERED FOR ELECTRONS.
C     THERE ARE CURRENTLY ONLY 10 OF THEM PLUS THE INTERPOLATION ARRAY.
C     DAVID ELDER 12/02/90
C     
      CTICHG = .FALSE.
      IF (CIOPTK.EQ.-1) THEN 
         CTICHG = .TRUE.
         CIOPTK = CIOPTG
         NTIBS = NTBS
         CTIBIN = CTBIN
         CLTIIN1 = CLTIN1
         CLTIIN2 = CLTIN2
         CGTIIN1 = CGTIN1
         CGTIIN2 = CGTIN2   
         CATIIN = CATIN
         CTIBOUL = CTBOUL
         CTIBOUG = CTBOUG
         CLTIOUL = CLTOUL
         CLTIOUG = CLTOUG
         DO COUNTI = 1,2
            DO COUNTJ = 1,MAXINS
               CTIBINS(COUNTJ,COUNTI) = CTBINS(COUNTJ,COUNTI)
            end do
         end do 
C     
C     ALSO COPY ANY GRADIENT INFORMATION
C     
         NTIG = NTEG
         DO COUNTI = 1,2
            DO COUNTJ = 1,MAXINS
               TMIG(COUNTJ,COUNTI) = TMEG(COUNTJ,COUNTI)     
            end do
         end do
      ENDIF
C     
C     
C     
      DO  IY = -NYS, NYS                                                    
C     
C-----------------------------------------------------------------------
C     INBOARD CONSTANT TEMPERATURE AND DENSITY (OPTION 0)                       
C-----------------------------------------------------------------------
C     
C     
         IF (CIOPTG.EQ.0) THEN                                                     
            DO  IX = IXOUT+1, NXS                                                
               CTEMBS(IX,IY,PZ) = CTBIN                                                 
               CRNBS(IX,IY,PZ)  = CNBIN                                                 
            end do
         ENDIF                                                                     
C     ION TEMPERATURE OPTION 0 
         IF (CIOPTK.EQ.0) THEN                                                     
            DO  IX = IXOUT+1, NXS                                                
               CTEMBSI(IX,IY,PZ) = CTIBIN                                               
            end do                                                                
         ENDIF                                                                     
C     
C     
C-----------------------------------------------------------------------
C     INBOARD TEMPERATURE (OPTIONS 1,2,4,5,6)                                   
C-----------------------------------------------------------------------
C     
C---- ATIN SPECIFIED AS 0 MEANS THAT THE SECOND DECAY FACTORS                   
C---- WILL APPLY TO THE ENTIRE INBOARD REGION.                                  
C     
         IF (CATIN.LE.0.0) THEN                                                    
            DO IX = IXOUT+1, NXS                                                
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.5.OR.CIOPTG.EQ.6) then 
                  CTEMBS(IX,IY,PZ) = CTBIN + CGTIN2 * XOUTS(IX)                          
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
                  CTEMBS(IX,IY,PZ) = CTBIN * EXP (XOUTS(IX) / CLTIN2)                    
               ENDIF                                                                 
            end do                                                                
C     
C---- ATIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL               
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C     
         ELSEIF (CATIN.GE.CA) THEN                                                 
            DO  IX = IXOUT+1, NXS                                                
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.4.OR.CIOPTG.EQ.6) THEN               
                  CTEMBS(IX,IY,PZ) = CTBIN + CGTIN1 * XOUTS(IX)                          
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
                  CTEMBS(IX,IY,PZ) = CTBIN * EXP (XOUTS(IX) / CLTIN1)                    
               ENDIF                                                                 
            end do
C     
C---- ATIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.  IN            
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C     
         ELSE                                                                      
            IXAT = IPOS (CATIN*0.999, XS, NXS-1)                                    
            DO  IX = IXOUT+1, IXAT                                               
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.4.OR.CIOPTG.EQ.6) THEN               
                  CTEMBS(IX,IY,PZ) = CTBIN + CGTIN1 * XOUTS(IX)                          
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
                  CTEMBS(IX,IY,PZ) = CTBIN * EXP (XOUTS(IX) / CLTIN1)                    
               ENDIF                                                                 
            end do                                                                
C     
            DO IX = IXAT+1, NXS                                                 
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.5.OR.CIOPTG.EQ.6) THEN               
                  CTEMBS(IX,IY,PZ) = CTEMBS(IXAT,IY,pz) +                                   
     >                 CGTIN2 * (XOUTS(IX)-XOUTS(IXAT))                    
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
                  CTEMBS(IX,IY,PZ) = CTEMBS(IXAT,IY,pz) *                                   
     >                 EXP ((XOUTS(IX)-XOUTS(IXAT))/CLTIN2)                
               ENDIF                                                                 
            end do                                                                
         ENDIF                                                                     
C     
C     
C-----------------------------------------------------------------------
C     ION INBOARD TEMPERATURE (OPTIONS 1,2,4,5,6)                      
C-----------------------------------------------------------------------
C     
C---- ATIIN SPECIFIED AS 0 MEANS THAT THE SECOND DECAY FACTORS                  
C---- WILL APPLY TO THE ENTIRE INBOARD REGION.                                  
C     
         IF (CATIIN.LE.0.0) THEN                                                  
            DO IX = IXOUT+1, NXS                                                
               IF (CIOPTK.EQ.1.OR.CIOPTK.EQ.5.OR.CIOPTK.EQ.6) THEN              
                  CTEMBSI(IX,IY,pz) = CTIBIN + CGTIIN2 * XOUTS(IX)                      
               ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.4) THEN                              
                  CTEMBSI(IX,IY,pz) = CTIBIN * EXP (XOUTS(IX) / CLTIIN2)               
               ENDIF                                                               
            end do
C     
C---- ATIIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL            
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C     
         ELSEIF (CATIIN.GE.CA) THEN                                               
            DO IX = IXOUT+1, NXS                                                
               IF (CIOPTK.EQ.1.OR.CIOPTK.EQ.4.OR.CIOPTK.EQ.6) THEN               
                  CTEMBSI(IX,IY,pz) = CTIBIN + CGTIIN1 * XOUTS(IX)                      
               ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.5) THEN                              
                  CTEMBSI(IX,IY,pz) = CTIBIN * EXP (XOUTS(IX) / CLTIIN1)               
               ENDIF                                                                 
            end do
C     
C---- ATIIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.IN          
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C     
         ELSE                                                                      
            IXAT = IPOS (CATIIN*0.999, XS, NXS-1)                                  
            DO  IX = IXOUT+1, IXAT                                               
               IF (CIOPTK.EQ.1.OR.CIOPTK.EQ.4.OR.CIOPTK.EQ.6) THEN               
                  CTEMBSI(IX,IY,pz) = CTIBIN + CGTIIN1 * XOUTS(IX)                      
               ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.5) THEN                              
                  CTEMBSI(IX,IY,pz) = CTIBIN * EXP (XOUTS(IX) / CLTIIN1)               
               ENDIF                                                                 
            end do
C     
            DO IX = IXAT+1, NXS                                                 
               IF (CIOPTK.EQ.1.OR.CIOPTK.EQ.5.OR.CIOPTK.EQ.6) THEN               
                  CTEMBSI(IX,IY,pz) = CTEMBSI(IXAT,IY,pz) +                              
     >                 CGTIIN2 * (XOUTS(IX)-XOUTS(IXAT))               
               ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.4) THEN                              
                  CTEMBSI(IX,IY,pz) = CTEMBSI(IXAT,IY,pz) *                               
     >                 EXP ((XOUTS(IX)-XOUTS(IXAT))/CLTIIN2)            
               ENDIF                                                                 
            end do
         ENDIF                                                                     
C     
C-----------------------------------------------------------------------
C     INBOARD DENSITY (OPTIONS 1,2,4,5)                                         
C-----------------------------------------------------------------------
C     
C---- ANIN SPECIFIED AS 0 MEANS THAT THE SECOND DECAY FACTORS                   
C---- WILL APPLY TO THE ENTIRE INBOARD REGION.                                  
C     
         IF (CANIN.LE.0.0) THEN                                                    
            DO  IX = IXOUT+1, NXS                                                
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.5) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN + CGNIN2 * XOUTS(IX)                           
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN * EXP (XOUTS(IX) / CLNIN2)                     
               ENDIF                                                                 
            end do
C     
C---- ANIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL               
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C     
         ELSEIF (CANIN.GE.CA) THEN                                                 
            DO  IX = IXOUT+1, NXS                                                
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.4) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN + CGNIN1 * XOUTS(IX)                           
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN * EXP (XOUTS(IX) / CLNIN1)                     
               ENDIF                                                                 
            end do
C     
C---- ANIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.  IN            
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C     
         ELSE                                                                      
            IXAN = IPOS (CANIN*0.999, XS, NXS-1)                                    
            DO  IX = IXOUT+1, IXAN                                               
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.4) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN + CGNIN1 * XOUTS(IX)                           
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
                  CRNBS(IX,IY,pz) = CNBIN * EXP (XOUTS(IX) / CLNIN1)                     
               ENDIF                                                                 
            end do
C     
            DO IX = IXAN+1, NXS                                                 
               IF (CIOPTG.EQ.1.OR.CIOPTG.EQ.5) THEN                              
                  CRNBS(IX,IY,pz) = CRNBS(IXAN,IY,pz) +                                     
     >                 CGNIN2 * (XOUTS(IX)-XOUTS(IXAN))                     
               ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
                  CRNBS(IX,IY,pz) = CRNBS(IXAN,IY,pz) *                                     
     >                 EXP ((XOUTS(IX)-XOUTS(IXAN))/CLNIN2)                 
               ENDIF                                                                 
            end do
         ENDIF                                                                     
C     
C-----------------------------------------------------------------------
C     INBOARD DENSITY (OPTION 6) NOTE 274    "STANDARD JET"                     
C-----------------------------------------------------------------------
C     
         IF (CIOPTG.EQ.6) THEN                                                     
            DO  IX = IXOUT+1, NXS                                                
               CRNBS(IX,IY,pz) = CNBIN + (CNBA-CNBIN) *                                 
     >          (1.0 - (1.0-XOUTS(IX)/CA)*(1.0-XOUTS(IX)/CA)) ** CGAMMA             
            end do
         ENDIF                                                                     
C     
C-----------------------------------------------------------------------
C     INBOARD CURVE FITTING OPTION FOR TEMPERATURE AND DENSITY                  
C-----------------------------------------------------------------------
C     
C---- OPTION 3 OR 7: CURVE FITTING INBOARD                                          
C     
         IF (CIOPTG.EQ.3.or.cioptg.eq.7) THEN                                                     
C     
C------INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
            CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >         NXS-IXOUT,XOUTS(IXOUT+1),CTEMBS(IXOUT+1,IY,pz),'LINEAR')          
C     
C------INTERPOLATE DENSITIES INBOARD FROM CUBIC SPLINE FIT.                    
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
            CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >         NXS-IXOUT,XOUTS(IXOUT+1),CRNBS (IXOUT+1,IY,pz),'LINEAR')          
         ENDIF                                                                     
C     
C---- OPTION 3 or 7: CURVE FITTING INBOARD                                          
C     
         IF (CIOPTK.EQ.3.or.cioptk.eq.7) THEN                                                     
C     
C------INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
            CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >         NXS-IXOUT,XOUTS(IXOUT+1),CTEMBSI(IXOUT+1,IY,pz),'LINEAR')        
C     
         ENDIF
C     
c     slmod tmp
c     
c     jdemod - steve has an option to override the plasma options with constant values inboard      
c     
         IF (SLOPT.EQ.1) THEN
            DO IX = IXOUT + 1, NXS      
               IF (XS(IX).GT.0.02) THEN
                  CTEMBSI(IX,IY,pz) = 100.0
                  CRNBS  (IX,IY,pz) = 1.0E20
               ENDIF  
            ENDDO
         ENDIF
c     slmod end
      end do
C     
C-----------------------------------------------------------------------
C     OUTBOARD TEMPERATURE AND DENSITY                                          
C     SET UP TWO SEPARATE SETS OF DETAILS :                                     
C     (A) OUTBOARD DETAILS IN QTEMBS,QRNBS FOR USE WITH NEUT LAUNCHING          
C     (B) FULL DETAILS OVER X,Y BINS IN CTEMBS,CRNBS FOR USE ELSEWHERE.         
C     
C-----------------------------------------------------------------------
C     
C     ALSO MODIFY THE BACKGROUND TEMPERATURES TO ACCOUNT FOR ANY 
C     APPLIED GRADIENTS AND THEN CALCULATE THE SPECIFIC GRADIENTS.
C     THE SPECIFIC GRADIENTS ARE CALCULATED BASED ON THE TEMPERATURE
C     DIFFERENCES BETWEEN ADJOINING BINS AND NOT ON THE INPUT 
C     TEMPERATURE MULTIPLIER DATA. THIS WAS DONE SO THAT AN 
C     ARBITRARY VARIATION IN TEMPERATURE COULD BE DEALT WITH, HOWEVER,
C     IT SHOULD BE RELATIVELY SIMPLE TO CONVERT TO USING THE INPUT
C     MULTIPLIER DATA TO CALCULATE GRADIENTS, IF NECESSARY.
C     
C     D.ELDER OCT 3, 1990
C     
C     
C---- OPTIONS 0           : CONSTANT                                            
C---- OPTIONS 1,2,3,4,5,6 : EXPONENTIAL DECAY                                   
c---- Options 7           : FITTED
C     
C     
C     SET UP GRADIENT MULTIPLIER FOR LIMITER EDGE POSITIONS
C     MULTIPLY INITIAL TEMPERATURES BY GRADIENT VALUE
C     
c     jdemod - July 2009 
c     - allow the base outboard density profile to be specified
c     in the input and fitted instead of just using an 
c     exponential 
c     
c     
c     
C     
      IF (NTEG.NE.0) THEN 
         TEGOUL = CTBOUL *  TMEG(1,2)
         TEGOUG = CTBOUG *  TMEG(1,2)
      ELSE
         TEGOUL = CTBOUL 
         TEGOUG = CTBOUG 
      ENDIF
      IF (NTIG.NE.0) THEN 
         TIGOUL = CTIBOUL * TMIG(1,2)
         TIGOUG = CTIBOUG * TMIG(1,2)
      ELSE
         TIGOUL = CTIBOUL
         TIGOUG = CTIBOUG
      ENDIF

c     
c     Base Density and temperature profiles along the limiter surface
c     
c     
      if (cioptg.eq.7) then 

C     
C-----------------------------------------------------------------------
C     OUTBOARD CURVE FITTING OPTION FOR TEMPERATURE AND DENSITY                  
C-----------------------------------------------------------------------
C     
C---- OPTION 7 : CURVE FITTING OUTBOARD
C     
C     
C------INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
         CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >        NQXSO-1,qxs(1-NQXSO),QTEMBS(1-nqxso,1),'LINEAR')          
C     
C------INTERPOLATE DENSITIES INBOARD FROM CUBIC SPLINE FIT.                    
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
         CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >        nqxso-1,qxs(1-nqxso),qrnbs(1-nqxso,1),'LINEAR')          

         qtembs(:,2) = qtembs(:,1)
         qrnbs(:,2) = qrnbs(:,1)

      else

         DO IQX = 1-NQXSO, 0                                                  
            IF (CIOPTG.EQ.0) THEN                                                   
               QTEMBS(IQX,1) = TEGOUL                             
               QTEMBS(IQX,2) = TEGOUG  
               QRNBS(IQX,1)  = CNBOUL                                                
               QRNBS(IQX,2)  = CNBOUG                                                
            ELSE                                                                    
               QTEMBS(IQX,1) = TEGOUL * EXP (QXS(IQX) / CLTOUL)        
               QTEMBS(IQX,2) = TEGOUG * EXP (QXS(IQX) / CLTOUG)            
               QRNBS(IQX,1)  = CNBOUL * EXP (QXS(IQX) / CLNOUL)                      
               QRNBS(IQX,2)  = CNBOUG * EXP (QXS(IQX) / CLNOUG)                      
            ENDIF                                                                   

         end do
C     
      endif

C     
C---- OPTION 7 : CURVE FITTING OUTBOARD
C     
      IF (CIOPTK.EQ.7) THEN                                                     
C     
C------INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
         CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >        nqxso-1,qxs(1-nqxso),qtembsi(1-nqxso,1),'LINEAR')        
C     
         qtembsi(:,2) = qtembsi(:,1)
c     
      else

         DO IQX = 1-NQXSO, 0                                                  
C     
C     BACKGROUND ION TEMPERATURES ALONG THE LIMITER FOR 
C     USE IN THE CALCULATION OF NEUTRAL PARTICLE INJECTION 
C     ENERGIES  .08/02/90. DAVID ELDER
C     
            IF (CIOPTK.EQ.0) THEN                                                   
               QTEMBSI(IQX,1) = TIGOUL                                     
               QTEMBSI(IQX,2) = TIGOUG                                    
            ELSE                                                                    
               QTEMBSI(IQX,1) = TIGOUL * EXP (QXS(IQX) / CLTIOUL)      
               QTEMBSI(IQX,2) = TIGOUG * EXP (QXS(IQX) / CLTIOUG)       
            ENDIF                                                                   
         end do

      endif
c     
c     Calculate base plasma conditions on the general grid
c     


      IF (CIOPTG.EQ.3.or.cioptg.eq.7) THEN                                                     
         CTBIN = CTBINS(1,2)                                      
         CNBIN = CNBINS(1,2)                                                     
         WRITE (6,'('' PLASMA: CTBIN,CNBIN='',1P,2G12.4)') CTBIN,CNBIN           
      ENDIF                                                                     
C     
      IF (CIOPTK.EQ.3.or.cioptk.eq.7) THEN                                                     
         CTIBIN = CTIBINS(1,2)                                                   
         WRITE (6,'('' PLASMA: CTIBIN='',1P,G12.4)') CTIBIN           
      ENDIF                                                                     
C     
      DO J = 1, 2                                                          
         QTEMBSI(1,J) = CTIBIN 
         QTEMBS(1,J) = CTBIN                                           
         QRNBS (1,J) = CNBIN                                                     
      end do
C     




      DO IY = -NYS, NYS                                                    
c     
c     Density and electron temperature
c     
         if (cioptg.eq.7) then 
c     
C     INTERPOLATE TEMPERATURES OUTBOARD FROM FITTER ROUTINE.                   
C     
            CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >           IXOUT,XOUTS(1),CTEMBS(1,IY,pz),'LINEAR')          
C     
C     INTERPOLATE DENSITIES OUTBOARD FROM FITTER ROUTINE.                   
c     
            CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >           IXOUT,XOUTS(1),CRNBS(1,IY,pz),'LINEAR')          
         else

            DO IX = 1, IXOUT                                                   
               IF (CIOPTG.EQ.0) THEN                                                 
                  IF (IY.GT.0) THEN                                                   
                     CTEMBS(IX,IY,PZ) = CTBOUG                                            
                     CRNBS(IX,IY,pz)  = CNBOUG                                            
                  ELSE                                                                
                     CTEMBS(IX,IY,PZ) = CTBOUL                                            
                     CRNBS(IX,IY,pz)  = CNBOUL                                            
                  ENDIF                                                               
               ELSE                                                                  
                  IF (IY.GT.0) THEN                                                   
                    CTEMBS(IX,IY,PZ) = CTBOUG * EXP (XOUTS(IX) / CLTOUG)                 
                    CRNBS(IX,IY,pz)  = CNBOUG * EXP (XOUTS(IX) / CLNOUG)                 
                  ELSE                                                                
                    CTEMBS(IX,IY,PZ) = CTBOUL * EXP (XOUTS(IX) / CLTOUL)                 
                    CRNBS(IX,IY,pz)  = CNBOUL * EXP (XOUTS(IX) / CLNOUL)                 
                  ENDIF                                                               
               ENDIF                                                                 
            end do
         endif
c     
c     Ion temperature
c     
         if (cioptk.eq.7) then

C------INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C     
            CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >           IXOUT,XOUTS(1),CTEMBSI(1,IY,pz),'LINEAR')        
C     

         else
            

            DO IX = 1, IXOUT                                                   

               IF (CIOPTK.EQ.0) THEN                                                 
                  IF (IY.GT.0) THEN                                                   
                     CTEMBSI(IX,IY,pz) = CTIBOUG                                      
                  ELSE                                                                
                     CTEMBSI(IX,IY,pz) = CTIBOUL                                       
                  ENDIF                                                               
               ELSE                                                                  
                  IF (IY.GT.0) THEN                                                   
                     CTEMBSI(IX,IY,pz) = CTIBOUG*EXP(XOUTS(IX)/CLTIOUG)              
                  ELSE                                                                
                     CTEMBSI(IX,IY,pz) = CTIBOUL*EXP(XOUTS(IX)/CLTIOUL)              
                  ENDIF                                                               
               ENDIF                                                                 
            end do

         endif
         
      end do




      
c     
c------------------------------------------------------------
c     Apply density and temperature gradients outboard
c------------------------------------------------------------
c     
C     
C     UNFORTUNATELY, THE METHOD THAT ALLOWED SIMPLE CALCULATION 
C     OF THE TEMPERATURES AT THE LIMITER EDGE WILL NOT WORK 
C     FARTHER OUT BECAUSE THE LIMITER SHAPE CHANGES AND THUS THE
C     CONNECTION LENGTH (LIMITER SEPARATION) AND THE POSITIONS 
C     OF THE CHANGING GRADIENT. CALCULATE AN ARRAY FROM 0 TO 2L
C     CONTAINING THE VARIOUS POINTS AT WHICH THE MULTIPLIERS ARE
C     SPECIFIED. A PARALLEL ARRAY CONTAINS THE ACTUAL MULTIPLIERS.
C     APPLY THESE TO THE VALUES IN THE CTEMBS,CTEMBSI ARRAYS. USE 
C     THE CALCULATED TEMPERATURES TO ESTIMATE BIN TO BIN GRADIENTS.
C     
C     D.ELDER, OCT 4, 1990
C     
C     ELECTRON TEMPERATURE GRADIENTS
C     
      IF (NTEG.NE.0) THEN
         DO 1130 IX = 1,IXOUT  
            TGPOS(1) = QEDGES(IQXS(IX),2)
            TGMULT(1) = TMEG(1,2)
            LEN1 = CL - TGPOS(1)
            DO 1140 IN = 2,NTEG
               TGPOS(IN) = TGPOS(1) + TMEG(IN,1) * LEN1
               TGMULT(IN) = TMEG(IN,2)
 1140       CONTINUE
            NTEMP = NTEG + 1
            TGMULT (NTEMP) = 1.0
            TGPOS (NTEMP) = CL
            LEN2 = CL - QEDGES(IQXS(IX),1)
            DO 1150 IN = NTEG,1,-1
               TGPOS(NTEMP+NTEG-IN+1) = CL+(1.0-TMEG(IN,1)) * LEN2
               TGMULT(NTEMP+NTEG-IN+1) = TMEG(IN,2)
 1150       CONTINUE
            NTEMP = NTEMP +NTEG
            TGMULT(NTEMP+1) = TMEG(1,2)               
C     
C     WRITE(6,*) 'NTEMP:',NTEMP 
C     WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C     WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C     
            CTEMBS(IX,0,pz) = CTEMBS(IX,0,pz) * TGMULT(1)   
            DO 1160 IY = 1,NYS
               IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
               IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
                  MULT = TGMULT(IN)
               ELSE
                  MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >                 * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
               ENDIF
               CTEMBS(IX,IY,PZ) = MULT * CTEMBS(IX,IY,PZ)
               CTEMBS(IX,IY-NYS-1,pz) = MULT * CTEMBS(IX,IY-NYS-1,pz)
 1160       CONTINUE 
C     WRITE(6,'(A,I5,2X,100F7.2)') 'CTEMBS: ',IX,
C     >                  (CTEMBS(IX,IY,PZ),IY=1,NYS)    
C     
C     END OF IX LOOP
C     
 1130    CONTINUE 
      ENDIF

C     
C     
C     ION TEMPERATURE GRADIENTS
C     
      IF (NTIG.NE.0) THEN
         DO 1330 IX = 1,IXOUT
            TGPOS(1) = QEDGES(IQXS(IX),2)
            TGMULT(1) = TMIG(1,2)
            LEN1 = CL - TGPOS(1)
            DO 1340 IN = 2,NTIG
               TGPOS(IN) = TGPOS(1) + TMIG(IN,1) * LEN1
               TGMULT(IN) = TMIG(IN,2)
 1340       CONTINUE
            NTEMP = NTIG + 1
            TGMULT (NTEMP) = 1.0
            TGPOS (NTEMP) = CL
            LEN2 = CL - QEDGES(IQXS(IX),1)
            DO 1350 IN = NTIG,1,-1
               TGPOS(NTEMP+NTIG-IN+1) = CL+(1.0-TMIG(IN,1)) * LEN2
               TGMULT(NTEMP+NTIG-IN+1) = TMIG(IN,2)
 1350       CONTINUE
            NTEMP = NTEMP +NTIG
            TGMULT(NTEMP+1) = TMIG(1,2)               
C     
C     WRITE(6,*) 'NTEMP:',NTEMP 
C     WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C     WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C     
            CTEMBS(IX,0,pz) = CTEMBS(IX,0,pz) * TGMULT(1)   
            DO 1360 IY = 1,NYS
               IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
               IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
                  MULT = TGMULT(IN)
               ELSE
                  MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >                 * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
               ENDIF
               CTEMBSI(IX,IY,pz) = MULT * CTEMBSI(IX,IY,pz)
               CTEMBSI(IX,IY-NYS-1,pz) = MULT * CTEMBSI(IX,IY-NYS-1,pz)
 1360       CONTINUE 
C     
C     END OF IX LOOP
C     
 1330    CONTINUE 
      ENDIF
c     
c     Apply a gradient (multiplier) to the outboard density as well. 
c     
      
      IF (NNBG.NE.0) THEN
         DO IX = 1,IXOUT  
            TGPOS(1) = QEDGES(IQXS(IX),2)
            TGMULT(1) = MNBG(1,2)
            LEN1 = CL - TGPOS(1)
            do IN = 2,NNBG
               TGPOS(IN) = TGPOS(1) + MNBG(IN,1) * LEN1
               TGMULT(IN) = MNBG(IN,2)
            end do
            NTEMP = NNBG + 1
            TGMULT (NTEMP) = 1.0
            TGPOS (NTEMP) = CL
            LEN2 = CL - QEDGES(IQXS(IX),1)
c     
            DO IN = NNBG,1,-1
               TGPOS(NTEMP+NTEG-IN+1) = CL+(1.0-MNBG(IN,1)) * LEN2
               TGMULT(NTEMP+NTEG-IN+1) = MNBG(IN,2)
            end do
c     
            NTEMP = NTEMP +NTEG
            TGMULT(NTEMP+1) = MNBG(1,2)               
C     
C     WRITE(6,*) 'NTEMP:',NTEMP 
C     WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C     WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C     
            CRNBS(IX,0,pz) = CRNBS(IX,0,pz) * TGMULT(1)   
            DO IY = 1,NYS
               IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
               IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
                  MULT = TGMULT(IN)
               ELSE
                  MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >                 * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
               ENDIF
               CRNBS(IX,IY,pz) = MULT * CRNBS(IX,IY,pz)
               CRNBS(IX,IY-NYS-1,pz) = MULT * CRNBS(IX,IY-NYS-1,pz)
            end do
C     WRITE(6,'(A,I5,2X,100F7.2)') 'CRNBS: ',IX,
C     >                  (CRNBS(IX,IY),IY=1,NYS)    
C     
C     END OF IX LOOP
C     
         end do
c     
      endif 

!     jdemod - the code to calculate the temperature gradients has
!     separated to make it possible to call after other background
!     modifications have been made.
!     This code will need to be enhanced when multiple plasma slices
!     are in use. 
      
c      call calculate_tgrad(qtim)



C     
C-----------------------------------------------------------------------
C     NORMAL AND ERROR RETURN POINTS ...                                        
C-----------------------------------------------------------------------
C     
C     DO 2000 IQX = 1-NQXSO, NQXSI, 20                                          
C     WRITE (6,'('' IQX='',I5,''  X='',1P,G11.4,''  TB,NB='',4G11.4)')        
C     >    IQX,QXS(IQX),(QTEMBS(IQX,J),QRNBS(IQX,J),J=1,2)                       
C     2000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       


      subroutine calculate_tgrad(qtim)
      use mod_params
      use mod_comt2
      use mod_comxyt
      use mod_comtor
      use mod_vtig
!     jdemod - this routine puts the code to calculate the temperature gradients
! in one routine- these are CTEGS and CTIGS
      implicit none
      real :: qtim

      ! locals
      integer :: ix,iy,IXOUT
      integer, external :: ipos
      real :: dstep,grad1,grad2,tgscal
      ! set pz = 1 - these routines assign a baseline plasma condition that
      ! will be copied to all plasma zones at the end of the routine
      integer :: pz
      pz = 1
c     
      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
C     
      ctigs = 0.0
      ctegs = 0.0

c      cALL RZERO(CTIGS,MAXNXS*(2*MAXNYS+1)) 
c      CALL RZERO(CTEGS,MAXNXS*(2*MAXNYS+1))
C     
C     SET SCALING FACTOR FOR TEMPERATURE GRADIENT FORCES
C     

      TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 

C     
C     CALCULATE THE TEMPERATURE GRADIENT FOR A BIN.
C     IT WILL BE THE AVERAGE OF THE GRADIENTS FROM THE CURRENT 
C     BIN CENTRE TO THE ADJACENT BIN CENTRES. FOR BINS AT THE  
C     ENDS THE GRADIENT IS ONLY BASED ON THE ONE NEIGHBOURING
C     BIN.
C     
      DO IX = 1,IXOUT
C     
C     DSTEP CONVERTS THE FORCE GRADIENT QUANTITIES TO 
C     A DISTANCE STEP FOR 1 TIMESTEP IN LIM (AT A GIVEN X POSITION)
C     THE QUANTITY IS deltaT^2 SINCE ONE dT COMES FROM THE VELOCITY 
C     EQUATION AND THE SECOND IN CONVERTING THE dV TO A DISTANCE.
C     
         DSTEP = TGSCAL * QS(IQXS(IX)) * QS(IQXS(IX))
         DO  IY = -NYS+1,-2
            GRAD1 = (CTEMBS(IX,IY+1,PZ)-CTEMBS(IX,IY,PZ)) 
     >           / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBS(IX,IY,PZ)-CTEMBS(IX,IY-1,PZ)) 
     >           / ( YOUTS(IY) -YOUTS(IY-1))
            CTEGS(IX,IY,pz) = DSTEP * (GRAD1+GRAD2)/2.0
         end do
         DO IY = 2,NYS-1
            GRAD1 = (CTEMBS(IX,IY+1,PZ)-CTEMBS(IX,IY,PZ)) 
     >           / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBS(IX,IY,PZ)-CTEMBS(IX,IY-1,PZ)) 
     >           / ( YOUTS(IY) -YOUTS(IY-1))
            CTEGS(IX,IY,pz) = DSTEP * (GRAD1+GRAD2)/2.0
         end do
C     
C     SPECIAL CASES - THE END POINTS 
C     THE NULL POINT 
C     
         CTEGS(IX,-NYS,pz) = (CTEMBS(IX,-NYS+1,PZ)-CTEMBS(IX,-NYS,pz))   
     >        / (YOUTS(-NYS+1)+YOUTS(-NYS)) * DSTEP
         CTEGS(IX,NYS,pz) = (CTEMBS(IX,NYS,pz)-CTEMBS(IX,NYS-1,pz))   
     >        / (YOUTS(NYS)-YOUTS(NYS-1)) * 1.6E-19 *DSTEP
         CTEGS(IX,1,pz) = (CTEMBS(IX,2,pz)-CTEMBS(IX,1,pz))   
     >        / (YOUTS(2)-YOUTS(1)) * DSTEP
         CTEGS(IX,-1,pz) = (CTEMBS(IX,-1,pz)-CTEMBS(IX,-2,pz))   
     >        / (YOUTS(-1)-YOUTS(-2)) * DSTEP
         CTEGS(IX,0,pz) = 0.0
      end do



C     
C     CALCULATE THE TEMPERATURE GRADIENT FOR A BIN.
C     IT WILL BE THE AVERAGE OF THE GRADIENTS FROM THE CURRENT 
C     BIN CENTRE TO THE ADJACENT BIN CENTRES. FOR BINS AT THE  
C     ENDS THE GRADIENT IS ONLY BASED ON THE ONE NEIGHBOURING
C     BIN.
C     
      DO  IX = 1,IXOUT
         DSTEP = TGSCAL *  QS(IQXS(IX)) * QS(IQXS(IX))
         DO  IY = -NYS+1,-2
            if (vtig_opt.eq.3) then 
               grad1 = vtig_tgrad(ix,iy,1)
               grad2 = grad1
            else
               GRAD1 = (CTEMBSI(IX,IY+1,pz)-CTEMBSI(IX,IY,pz)) 
     >           / ( YOUTS(IY+1) -YOUTS(IY))
               GRAD2 = (CTEMBSI(IX,IY,pz)-CTEMBSI(IX,IY-1,pz)) 
     >           / ( YOUTS(IY) -YOUTS(IY-1))
            endif
            CTIGS(IX,IY,pz) = DSTEP * (GRAD1+GRAD2)/2.0
         end do
         DO  IY = 2,NYS-1
            if (vtig_opt.eq.3) then 
               grad1 = vtig_tgrad(ix,iy,1)
               grad2 = grad1
            else
               GRAD1 = (CTEMBSI(IX,IY+1,pz)-CTEMBSI(IX,IY,pz)) 
     >           / ( YOUTS(IY+1) -YOUTS(IY))
               GRAD2 = (CTEMBSI(IX,IY,pz)-CTEMBSI(IX,IY-1,pz)) 
     >           / ( YOUTS(IY) -YOUTS(IY-1))
            endif
            CTIGS(IX,IY,pz) = DSTEP * (GRAD1+GRAD2)/2.0
         end do
C     
C     SPECIAL CASES - THE END POINTS 
C     THE NULL POINT 
C     
         
         if (vtig_opt.eq.3) then 
            CTIGS(IX,-NYS,pz) = vtig_tgrad(ix,-nys,1) * dstep
            CTIGS(IX,NYS,pz) = vtig_tgrad(ix,nys,1) * dstep
            CTIGS(IX,1,pz) = vtig_tgrad(ix,1,1) * dstep
            CTIGS(IX,-1,pz) = vtig_tgrad(ix,-1,1) * dstep
            CTIGS(IX,0,pz) = 0.0
         else
           CTIGS(IX,-NYS,pz)=(CTEMBSI(IX,-NYS+1,pz)-CTEMBSI(IX,-NYS,pz))   
     >        / (YOUTS(-NYS+1)+YOUTS(-NYS)) *DSTEP
           CTIGS(IX,NYS,pz) =(CTEMBSI(IX,NYS,pz)-CTEMBSI(IX,NYS-1,pz))   
     >        / (YOUTS(NYS)-YOUTS(NYS-1))  * DSTEP
            CTIGS(IX,1,pz) = (CTEMBSI(IX,2,pz)-CTEMBSI(IX,1,pz))   
     >        / (YOUTS(2)-YOUTS(1)) * 1.6E-19 * DSTEP
            CTIGS(IX,-1,pz) = (CTEMBSI(IX,-1,pz)-CTEMBSI(IX,-2,pz))   
     >        / (YOUTS(-1)-YOUTS(-2)) * DSTEP
            CTIGS(IX,0,pz) = 0.0
         endif
      end do

      ! jdemod - after the temperature gradients have been calculated the default LIM
      ! background is complete - copy to all poloidal zones

      ! Loop over pz and assign the background plasma properties calculated here
      ! to every poloidal zone. 
      do pz = 2,maxpzone
         crnbs(:,:,pz) = crnbs(:,:,1)
         ctembs(:,:,pz) = ctembs(:,:,1)
         ctembsi(:,:,pz) = ctembsi(:,:,1)
         ctegs(:,:,pz) = ctigs(:,:,1)
         ctigs(:,:,pz) = ctigs(:,:,1)
      end do
      
      return
      end


      subroutine calculate_efield(qtim,limiz)
      use mod_params
      use mod_comt2
      use mod_comxyt
      use mod_comtor
      use mod_vtig
      implicit none

      ! jdemod - this routine puts the code to calculate the electric fields for the base plasma options
      !          in one routine
      ! sets default values of efield and velplasma arrays
      ! this routine is called BEFORE any plasma overlays have been applied.
      ! As a result the pz loops could be replaced with a calculation for pz=1 and assignment
      ! from 1 to any other plasma zones. However, in case there are any changes to the code elswhere
      ! that might change the background values upstream from this code - leave the pz loops for now. 
      real :: qtim
      integer :: limiz
      
      ! locals
      integer :: ix,iy,IXOUT,iqx,iqy
      integer, external :: ipos
      real :: dstep,grad1,grad2,tgscal

      real :: fex,fexz
      integer :: pz,iz

      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
C                                                                               
c     Assign default values to efield and velplasma
c     
c     IQYS and IQXS should be setup to map IX and IY to IQX ad IQY
c     Need to include multiplying by the radial scale factors  ?
c     Consider removing CVHXS and CVEXZS      
c      
!
      do pz = 1,maxpzone
         ! outboard
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
         ! inboard
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
c
      FEX = QTIM * QTIM * (1.602192E-19 / (1.672614E-27 * CRMI))                

      IF (LIMIZ.GT.0) THEN                                                      

        do pz = 1,maxpzone

        DO 300  IZ = 1, LIMIZ                                                   
          FEXZ = FEX * REAL (IZ)                                                
          DO 250 IY = -NYS, NYS                                                 
           ! changed to NXS to support transport forces inboard of the probe tip
           !DO 250 IX = 1, IXOUT                                                 
           DO 250 IX = 1, NXS
            IQX = IQXS(IX)                                                      

        !    if (vel_efield_opt.eq.0) then
        ! The velocity and efield from vel_efield_opt =1 are scaled to cancel the temperature factor in
        ! CFVHXS and CFEXZS - this lets plasma backgrounds calculated with the old options and those
        ! calculated with the new ones (soledge/sol22) to coexist in the same simulation.        
               if (ix.gt.ixout) then 
                  CFEXZS(IX,IY,IZ,pz) = FEXZ * CTEMBS(IX,IY,pz)/CTBIN 
     >                             * QS(IQX) * QS(IQX)           
               else
                  CFEXZS(IX,IY,IZ,pz) = FEXZ * CTEMBS(IX,IY,pz)/CTBIN *                     
     >                         CYSCLS(IQX)/YSCALE * QS(IQX) * QS(IQX)           
               endif
         !   elseif (vel_efield_opt.eq.1) then 
               !  if using velplasma/efield values then the CFEHXS contains only timestep and charge state
               ! scaling and not temperature relative to the separatrix
         !      if (ix.gt.ixout) then 
         !         CFEXZS(IX,IY,IZ,pz) = FEXZ * QS(IQX) * QS(IQX)           
         !      else
                  ! not sure about the cyscls/yscale factor for efield - leave for now
         !         CFEXZS(IX,IY,IZ,pz) = FEXZ *                      
!     >                         CYSCLS(IQX)/YSCALE * QS(IQX) * QS(IQX)           
         !      endif
         !   endif
               
 250        CONTINUE                                                              
  300   CONTINUE                                                                

        end do
      ENDIF                                                                     
C                                                                               
      do pz = 1,maxpzone
      DO  IY = -NYS, NYS                                                     
       ! changed to NXS to support transport forces inboard of the probe tip
       DO  IX = 1, NXS                                                     
       !DO 310 IX = 1, IXOUT
        IQX = IQXS(IX)                                                          
        ! The velocity and efield from vel_efield_opt =1 are scaled to cancel the temperature factor in
        ! CFVHXS and CFEXZS - this lets plasma backgrounds calculated with the old options and those
        ! calculated with the new ones (soledge/sol22) to coexist in the same simulation.        

        !if (vel_efield_opt.eq.0) then 
           CFVHXS(IX,IY,pz) = 
     >        SQRT((CTEMBS(IX,IY,pz)+CTEMBSI(IX,IY,pz))/(CTBIN+CTIBIN))
     >        * QTIM * QS(IQX)             
        !elseif (vel_efield_opt.eq.1) then
           !  if using velplasma/efield values then the CFVHXS contains only timestep
           ! scaling and not temperature relative to the separatrix
        !   CFVHXS(IX,IY,pz) = QTIM * QS(IQX)             
        !endif   
           
          end do
        end do
      end do

      
      return
      end
