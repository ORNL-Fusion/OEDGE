      SUBROUTINE PLASMA (NTBS,NTIBS,NNBS,CIOPTG,CIOPTK,QTIM)
      IMPLICIT  none
      INTEGER   NTBS,NTIBS,NNBS,CIOPTG,CIOPTK 
      REAL      QTIM                            
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comxyt'                                                        
C     INCLUDE   (COMXYT)                                                        
      INCLUDE   'comtor'                                                        
C     INCLUDE   (COMTOR)                                                        
      INCLUDE   'comt2'                                                         
C     INCLUDE   (COMT2)                                                         
c slmod begin
      INCLUDE   'slcom'
c slmod end      
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  PLASMA                                                           *        
C  *  ------                                                           *        
C  *  ROUTINE TO CALCULATE THE TEMPERATURE AND ION DENSITY OF THE      *        
C  *  BACKGROUND PLASMA AT A SET OF POINTS GIVEN BY QXS.  VARIOUS      *        
C  *  OPTIONS ARE ALLOWED TO GIVE COMBINATIONS OF EXPONENTIAL AND      *        
C  *  LINEAR DECAY.                                                    *        
C  *                                                                   *        
C  *  AMENDMENTS 16/8/88 ... ARRAYS CTEMBS, CRNBS NOW DEPENDENT ON     *        
C  *  X AND Y, EVERYTHING HAS TO BE CHANGED.  RETAIN QTEMBS AND QRNBS  *        
C  *  DIMENSIONED FOR QX AND RIGHT & LEFT SIDES OF LIMITER JUST FOR    *        
C  *  USE IN NEUT WHEN CALCULATING LAUNCH DISTRIBUTION FUNCTIONS.      *        
C  *                                                                   *
C  *  AMENDED 06/02/90 ... ARRAY CTEMBSI - CONTAINS THE ION TEMPERATURE*
C  *  DISTRIBUTION AS A FUNCTION OF X,Y. THIS IS THE FIRST STEP IN     *
C  *  CREATING COMPLETELY INDEPENDENT DIST'NS FOR IONS AND ELECTRONS.  *
C  *  USED FOR CALCULATING A NUMBER OF THE TAU VALUES (EXCEPT TIZ)     *
C  *  DAVID ELDER , FEB 6 , 1990      (08/02/90 - ALSO QTEMBSI)        *
C  *                                                                   *
C  *  AMENDED 25/09/90 ... ARRAYS CTIGS,CTEGS CONTAINING THE Y AXIS    *
C  *  TEMPERATURE GRADIENTS. THE TEMPERATURES ARE MODIFIED BY GIVEN    *
C  *  INPUT FUNCTIONS AND THE GRADIENT FOR EACH BIN IS CALCULATED      *
C  *  FROM THE RESULTING TEMPERATURES.                                 *
C  *  DAVID ELDER , SEPT 25 , 1990                                     *
C  *                                                                   *        
C  *  ARGUMENTS (ARRAYS FROM COMMONS) :-                               *        
C  *  QXS    : SET OF X POSITIONS TEMP AND DENSITY TO BE CALCULATED AT *        
C  *  CTEMBS : ARRAY TEMPERATURES TO BE RETURNED IN X,Y                *        
C  *  CTEMBSI: ARRAY OF ION TEMPERATURES INDEXED BY IX,IY              *
C  *  CRNBS  : ARRAY DENSITIES TO BE RETURNED IN X,Y                   *        
C  *  QTEMBS : ARRAY TEMPERATURES TO BE RETURNED IN QX,RIGHT & LEFT    *        
C  *  QTEMBSI: ARRAY OF ION TEMPERATURES INDEXED BY QX,RIGHT & LEFT    *        
C  *  QRNBS  : ARRAY DENSITIES TO BE RETURNED IN QX, RIGHT & LEFT      *        
C  *  NQXSO  : NUMBER OF POINTS (1-NQXSO:0) OUTBOARD                   *        
C  *  NQXSI  : NUMBER OF POINTS (1:NQXSI)   INBOARD                    *        
C  *  CA     : POSITION OF CENTRE "A"                                  *        
C  *  ELECTRON TEMPERATURE VARIABLES                                   *
C  *  CTBINS : SET OF (X,TB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C  *  CTBIN  : INBOARD  TEMPERATURE OF PLASMA BASE VALUE               *        
C  *  CLTIN1 : INBOARD  TEMP EXPONENTIAL DECAY   FOR   0   < X < CATIN *        
C  *  CLTIN2 : INBOARD  TEMP EXPONENTIAL DECAY   FOR CATIN < X <  CA   *        
C  *  CGTIN1 : INBOARD  TEMP LINEAR DECAY FACTOR FOR   0   < X < CATIN *        
C  *  CGTIN2 : INBOARD  TEMP LINEAR DECAY FACTOR FOR CATIN < X <  CA   *        
C  *  CATIN  : INBOARD  TEMP CROSSOVER POINT                           *        
C  *  CTBOUL : OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y < 0         *        
C  *  CTBOUG : OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y > 0         *        
C  *  NTEG   : NUMBER OF DIFFERENT Te MULTIPLIERS/REGIONS              *  
C  *  TMEG   : ARRAY SPECIFYING THE DIFFERENT REGIONS AND MULTIPLIERS  *
C  *  CTEGS  : CALCULATED Te TEMPERATURE GRADIENT FOR EACH BIN         *
C  *  ION TEMPERATURE PARAMETERS                                       *
C  *  CTIBINS: SET OF (X,TB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C  *  CTIBIN : INBOARD  TEMPERATURE OF PLASMA BASE VALUE               *        
C  *  CLTIIN1: INBOARD  TEMP EXPONENTIAL DECAY   FOR   0   < X < CATIIN*        
C  *  CLTIIN2: INBOARD  TEMP EXPONENTIAL DECAY   FOR CATIN < X <  CA   *        
C  *  CGTIIN1: INBOARD  TEMP LINEAR DECAY FACTOR FOR   0   < X < CATIIN*        
C  *  CGTIIN2: INBOARD  TEMP LINEAR DECAY FACTOR FOR CATIN < X <  CA   *        
C  *  CATIIN : INBOARD  TEMP CROSSOVER POINT                           *        
C  *  CTIBOUL: OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y < 0         *        
C  *  CTIBOUG: OUTBOARD TEMPERATURE OF PLASMA BASE VALUE Y > 0         *        
C  *  NTIG   : NUMBER OF DIFFERENT Ti MULTIPLIERS/REGIONS              *  
C  *  TMIG   : ARRAY SPECIFYING THE DIFFERENT REGIONS AND MULTIPLIERS  *
C  *  CTIGS  : CALCULATED Ti TEMPERATURE GRADIENT FOR EACH BIN         *
C  *  DENSITY PARAMETERS                                               *
C  *  CNBINS : SET OF (X,NB) POINTS FOR CURVE FITTING  (OPTION 3)      *        
C  *  CNBIN  : INBOARD  ION DENSITY OF PLASMA BASE VALUE               *        
C  *  CLNIN1 : INBOARD  DENS EXPONENTIAL DECAY   FOR   0   < X < CANIN *        
C  *  CLNIN2 : INBOARD  DENS EXPONENTIAL DECAY   FOR CANIN < X <  CA   *        
C  *  CGNIN1 : INBOARD  DENS LINEAR DECAY FACTOR FOR   0   < X < CANIN *        
C  *  CGNIN2 : INBOARD  DENS LINEAR DECAY FACTOR FOR CANIN < X <  CA   *        
C  *  CANIN  : INBOARD  DENSITY CROSSOVER POINT                        *        
C  *  CNBOUL : OUTBOARD ION DENSITY OF PLASMA BASE VALUE Y < 0         *        
C  *  CNBOUG : OUTBOARD ION DENSITY OF PLASMA BASE VALUE Y > 0         *        
C  *  NTBS   : NUMBER OF PAIRED VALUES IN CTBINS ARRAY                 *        
C  *  NTIBS  : NUMBER OF PAIRED VALUES IN CTIBINS ARRAY                *        
C  *  NNBS   : NUMBER OF PAIRED VALUES IN CNBINS ARRAY                 *        
C  *                                                                   *        
C  *  OPTION    OUTBOARD      INBOARD (1)          INBOARD (2)         *        
C  *  ------------------------------------------------------------     *        
C  *   0     CONSTANT          CONSTANT             CONSTANT           *        
C  *   1     EXPONENTIAL       LINEAR               LINEAR             *        
C  *   2     EXPONENTIAL       EXPONENTIAL          EXPONENTIAL        *        
C  *   3     EXPONENTIAL       FITTED TO SET OF GIVEN VALUES           *        
C  *   4     EXPONENTIAL       LINEAR               EXPONENTIAL        *        
C  *   5     EXPONENTIAL       EXPONENTIAL          LINEAR             *        
C  *   6     EXPONENTIAL       STANDARD JET  (NOTE 274)                *        
C  *                                                                   *        
C  *  CHRIS FARRELL  (HUNTERSKIL)  APRIL 1988                          *        
C  *  DAVID ELDER                  OCTOBER 1990
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IPOS,IXAT,IXAN,IXOUT,IX,IY,IQX,J                                  
      INTEGER COUNTI,COUNTJ
      INTEGER NTEMP,IN
      REAL TEGOUL,TEGOUG,TIGOUL,TIGOUG
      REAL TGPOS(2*MAXINS+2),TGMULT(2*MAXINS+1)
      REAL LEN1,LEN2,MULT,GRAD1,GRAD2,DSTEP
      REAL TGSCAL
C                                                                               
      IXOUT = IPOS (-1.E-10, XS, NXS-1)                                         
C  
C    IF CIOPTK = -1 THEN USE Ti = Te
C    DO THIS BY SETTING ALL OF THE ION CONSTANTS EQUAL TO THOSE 
C    ENTERED FOR ELECTRONS.
C    THERE ARE CURRENTLY ONLY 10 OF THEM PLUS THE INTERPOLATION ARRAY.
C    DAVID ELDER 12/02/90
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
         DO 50 COUNTI = 1,2
            DO 50 COUNTJ = 1,MAXINS
               CTIBINS(COUNTJ,COUNTI) = CTBINS(COUNTJ,COUNTI)
  50     CONTINUE
C
C        ALSO COPY ANY GRADIENT INFORMATION
C
         NTIG = NTEG
         DO 60 COUNTI = 1,2
            DO 60 COUNTJ = 1,MAXINS
               TMIG(COUNTJ,COUNTI) = TMEG(COUNTJ,COUNTI)     
  60     CONTINUE
      ENDIF
C
C
C                                                                               
      DO 1000 IY = -NYS, NYS                                                    
C                                                                               
C-----------------------------------------------------------------------        
C     INBOARD CONSTANT TEMPERATURE AND DENSITY (OPTION 0)                       
C-----------------------------------------------------------------------        
C                                                                               
C                                                                              
      IF (CIOPTG.EQ.0) THEN                                                     
        DO 100 IX = IXOUT+1, NXS                                                
          CTEMBS(IX,IY) = CTBIN                                                 
          CRNBS(IX,IY)  = CNBIN                                                 
  100   CONTINUE                                                                
      ENDIF                                                                     
C     ION TEMPERATURE OPTION 0 
      IF (CIOPTK.EQ.0) THEN                                                     
        DO 112 IX = IXOUT+1, NXS                                                
          CTEMBSI(IX,IY) = CTIBIN                                               
  112   CONTINUE                                                                
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
        DO 200 IX = IXOUT+1, NXS                                                
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.5.OR.CIOPTG.EQ.6) THEN               
            CTEMBS(IX,IY) = CTBIN + CGTIN2 * XOUTS(IX)                          
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
            CTEMBS(IX,IY) = CTBIN * EXP (XOUTS(IX) / CLTIN2)                    
          ENDIF                                                                 
  200   CONTINUE                                                                
C                                                                               
C---- ATIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL               
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C                                                                               
      ELSEIF (CATIN.GE.CA) THEN                                                 
        DO 300 IX = IXOUT+1, NXS                                                
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.4.OR.CIOPTG.EQ.6) THEN               
            CTEMBS(IX,IY) = CTBIN + CGTIN1 * XOUTS(IX)                          
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
            CTEMBS(IX,IY) = CTBIN * EXP (XOUTS(IX) / CLTIN1)                    
          ENDIF                                                                 
  300   CONTINUE                                                                
C                                                                               
C---- ATIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.  IN            
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C                                                                               
      ELSE                                                                      
        IXAT = IPOS (CATIN*0.999, XS, NXS-1)                                    
        DO 400 IX = IXOUT+1, IXAT                                               
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.4.OR.CIOPTG.EQ.6) THEN               
            CTEMBS(IX,IY) = CTBIN + CGTIN1 * XOUTS(IX)                          
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
            CTEMBS(IX,IY) = CTBIN * EXP (XOUTS(IX) / CLTIN1)                    
          ENDIF                                                                 
  400   CONTINUE                                                                
C                                                                               
        DO 500 IX = IXAT+1, NXS                                                 
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.5.OR.CIOPTG.EQ.6) THEN               
            CTEMBS(IX,IY) = CTEMBS(IXAT,IY) +                                   
     >                      CGTIN2 * (XOUTS(IX)-XOUTS(IXAT))                    
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
            CTEMBS(IX,IY) = CTEMBS(IXAT,IY) *                                   
     >                      EXP ((XOUTS(IX)-XOUTS(IXAT))/CLTIN2)                
          ENDIF                                                                 
  500   CONTINUE                                                                
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
        DO 210 IX = IXOUT+1, NXS                                                
          IF     (CIOPTK.EQ.1.OR.CIOPTK.EQ.5.OR.CIOPTK.EQ.6) THEN              
            CTEMBSI(IX,IY) = CTIBIN + CGTIIN2 * XOUTS(IX)                      
          ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.4) THEN                              
            CTEMBSI(IX,IY) = CTIBIN * EXP (XOUTS(IX) / CLTIIN2)               
          ENDIF                                                               
  210   CONTINUE                                                                
C                                                                               
C---- ATIIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL            
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C                                                                               
      ELSEIF (CATIIN.GE.CA) THEN                                               
        DO 310 IX = IXOUT+1, NXS                                                
          IF     (CIOPTK.EQ.1.OR.CIOPTK.EQ.4.OR.CIOPTK.EQ.6) THEN               
            CTEMBSI(IX,IY) = CTIBIN + CGTIIN1 * XOUTS(IX)                      
          ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.5) THEN                              
            CTEMBSI(IX,IY) = CTIBIN * EXP (XOUTS(IX) / CLTIIN1)               
          ENDIF                                                                 
  310   CONTINUE                                                                
C                                                                               
C---- ATIIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.IN          
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C                                                                               
      ELSE                                                                      
        IXAT = IPOS (CATIIN*0.999, XS, NXS-1)                                  
        DO 410 IX = IXOUT+1, IXAT                                               
          IF     (CIOPTK.EQ.1.OR.CIOPTK.EQ.4.OR.CIOPTK.EQ.6) THEN               
            CTEMBSI(IX,IY) = CTIBIN + CGTIIN1 * XOUTS(IX)                      
          ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.5) THEN                              
            CTEMBSI(IX,IY) = CTIBIN * EXP (XOUTS(IX) / CLTIIN1)               
          ENDIF                                                                 
  410   CONTINUE                                                                
C                                                                               
        DO 510 IX = IXAT+1, NXS                                                 
          IF     (CIOPTK.EQ.1.OR.CIOPTK.EQ.5.OR.CIOPTK.EQ.6) THEN               
            CTEMBSI(IX,IY) = CTEMBSI(IXAT,IY) +                              
     >                      CGTIIN2 * (XOUTS(IX)-XOUTS(IXAT))               
          ELSEIF (CIOPTK.EQ.2.OR.CIOPTK.EQ.4) THEN                              
            CTEMBSI(IX,IY) = CTEMBSI(IXAT,IY) *                               
     >                      EXP ((XOUTS(IX)-XOUTS(IXAT))/CLTIIN2)            
          ENDIF                                                                 
  510   CONTINUE                                                                
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
        DO 600 IX = IXOUT+1, NXS                                                
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.5) THEN                              
            CRNBS(IX,IY) = CNBIN + CGNIN2 * XOUTS(IX)                           
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
            CRNBS(IX,IY) = CNBIN * EXP (XOUTS(IX) / CLNIN2)                     
          ENDIF                                                                 
  600   CONTINUE                                                                
C                                                                               
C---- ANIN SPECIFIED AS A MEANS THAT THE FIRST DECAY FACTORS WILL               
C---- APPLY TO THE ENTIRE INBOARD REGION.                                       
C                                                                               
      ELSEIF (CANIN.GE.CA) THEN                                                 
        DO 700 IX = IXOUT+1, NXS                                                
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.4) THEN                              
            CRNBS(IX,IY) = CNBIN + CGNIN1 * XOUTS(IX)                           
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
            CRNBS(IX,IY) = CNBIN * EXP (XOUTS(IX) / CLNIN1)                     
          ENDIF                                                                 
  700   CONTINUE                                                                
C                                                                               
C---- ANIN BETWEEN 0 AND A MEANS WE WILL HAVE A CROSSOVER POINT.  IN            
C---- THE FIRST SECTION, THE FIRST DECAY FACTORS WILL APPLY, AND IN             
C---- THE SECOND SECTION, THE SECOND SET OF FACTORS WILL APPLY.                 
C                                                                               
      ELSE                                                                      
        IXAN = IPOS (CANIN*0.999, XS, NXS-1)                                    
        DO 800 IX = IXOUT+1, IXAN                                               
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.4) THEN                              
            CRNBS(IX,IY) = CNBIN + CGNIN1 * XOUTS(IX)                           
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.5) THEN                              
            CRNBS(IX,IY) = CNBIN * EXP (XOUTS(IX) / CLNIN1)                     
          ENDIF                                                                 
  800   CONTINUE                                                                
C                                                                               
        DO 900 IX = IXAN+1, NXS                                                 
          IF     (CIOPTG.EQ.1.OR.CIOPTG.EQ.5) THEN                              
            CRNBS(IX,IY) = CRNBS(IXAN,IY) +                                     
     >                     CGNIN2 * (XOUTS(IX)-XOUTS(IXAN))                     
          ELSEIF (CIOPTG.EQ.2.OR.CIOPTG.EQ.4) THEN                              
            CRNBS(IX,IY) = CRNBS(IXAN,IY) *                                     
     >                     EXP ((XOUTS(IX)-XOUTS(IXAN))/CLNIN2)                 
          ENDIF                                                                 
  900   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     INBOARD DENSITY (OPTION 6) NOTE 274    "STANDARD JET"                     
C-----------------------------------------------------------------------        
C                                                                               
      IF (CIOPTG.EQ.6) THEN                                                     
        DO 950 IX = IXOUT+1, NXS                                                
          CRNBS(IX,IY) = CNBIN + (CNBA-CNBIN) *                                 
     >      (1.0 - (1.0-XOUTS(IX)/CA)*(1.0-XOUTS(IX)/CA)) ** CGAMMA             
  950   CONTINUE                                                                
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
C------ INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
        CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >           NXS-IXOUT,XOUTS(IXOUT+1),CTEMBS(IXOUT+1,IY),'LINEAR')          
C                                                                               
C------ INTERPOLATE DENSITIES INBOARD FROM CUBIC SPLINE FIT.                    
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
        CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >           NXS-IXOUT,XOUTS(IXOUT+1),CRNBS (IXOUT+1,IY),'LINEAR')          
      ENDIF                                                                     
C                                                                               
C---- OPTION 3 or 7: CURVE FITTING INBOARD                                          
C                                                                               
      IF (CIOPTK.EQ.3.or.cioptk.eq.7) THEN                                                     
C                                                                               
C------ INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
        CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >           NXS-IXOUT,XOUTS(IXOUT+1),CTEMBSI(IXOUT+1,IY),'LINEAR')        
C                                                                               
      ENDIF
C
c slmod tmp
      IF (SLOPT.EQ.1) THEN
        DO IX = IXOUT + 1, NXS      
          IF (XS(IX).GT.0.02) THEN
            CTEMBSI(IX,IY) = 100.0
            CRNBS  (IX,IY) = 1.0E20
          ENDIF  
        ENDDO
      ENDIF
c slmod end
 1000 CONTINUE                                                                  
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
c            - allow the base outboard density profile to be specified
c              in the input and fitted instead of just using an 
c              exponential 
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
C------ INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
          CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >           NQXSO-1,qxs(1-NQXSO),QTEMBS(1-nqxso,1),'LINEAR')          
C                                                                               
C------ INTERPOLATE DENSITIES INBOARD FROM CUBIC SPLINE FIT.                    
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
          CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >           nqxso-1,qxs(1-nqxso),qrnbs(1-nqxso,1),'LINEAR')          

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
C------ INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
          CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >           nqxso-1,qxs(1-nqxso),qtembsi(1-nqxso,1),'LINEAR')        
C                                                                               
          qtembsi(:,2) = qtembsi(:,1)
c
      else

         DO IQX = 1-NQXSO, 0                                                  
C
C          BACKGROUND ION TEMPERATURES ALONG THE LIMITER FOR 
C          USE IN THE CALCULATION OF NEUTRAL PARTICLE INJECTION 
C          ENERGIES  .08/02/90. DAVID ELDER
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
c       Density and electron temperature
c
         if (cioptg.eq.7) then 
c
C           INTERPOLATE TEMPERATURES OUTBOARD FROM FITTER ROUTINE.                   
C                                                                               
           CALL FITTER (NTBS,CTBINS(1,1),CTBINS(1,2),                              
     >             IXOUT,XOUTS(1),CTEMBS(1,IY),'LINEAR')          
C                                                                               
C           INTERPOLATE DENSITIES OUTBOARD FROM FITTER ROUTINE.                   
c
           CALL FITTER (NNBS,CNBINS(1,1),CNBINS(1,2),                              
     >           IXOUT,XOUTS(1),CRNBS(1,IY),'LINEAR')          
        else

           DO IX = 1, IXOUT                                                   
             IF (CIOPTG.EQ.0) THEN                                                 
               IF (IY.GT.0) THEN                                                   
                 CTEMBS(IX,IY) = CTBOUG                                            
                 CRNBS(IX,IY)  = CNBOUG                                            
               ELSE                                                                
                 CTEMBS(IX,IY) = CTBOUL                                            
                 CRNBS(IX,IY)  = CNBOUL                                            
               ENDIF                                                               
             ELSE                                                                  
               IF (IY.GT.0) THEN                                                   
                 CTEMBS(IX,IY) = CTBOUG * EXP (XOUTS(IX) / CLTOUG)                 
                 CRNBS(IX,IY)  = CNBOUG * EXP (XOUTS(IX) / CLNOUG)                 
               ELSE                                                                
                 CTEMBS(IX,IY) = CTBOUL * EXP (XOUTS(IX) / CLTOUL)                 
                 CRNBS(IX,IY)  = CNBOUL * EXP (XOUTS(IX) / CLNOUL)                 
               ENDIF                                                               
             ENDIF                                                                 
          end do
       endif
c
c     Ion temperature
c
       if (cioptk.eq.7) then

C------ INTERPOLATE TEMPERATURES INBOARD FROM FITTER ROUTINE.                   
C------ EXTRAPOLATE OUTER VALUES AS CONSTANTS IF REQUIRED                       
C                                                                               
        CALL FITTER (NTIBS,CTIBINS(1,1),CTIBINS(1,2),                           
     >            IXOUT,XOUTS(1),CTEMBSI(1,IY),'LINEAR')        
C                                                                               

      else
       

         DO IX = 1, IXOUT                                                   

          IF (CIOPTK.EQ.0) THEN                                                 
            IF (IY.GT.0) THEN                                                   
              CTEMBSI(IX,IY) = CTIBOUG                                      
            ELSE                                                                
              CTEMBSI(IX,IY) = CTIBOUL                                       
            ENDIF                                                               
          ELSE                                                                  
            IF (IY.GT.0) THEN                                                   
              CTEMBSI(IX,IY) = CTIBOUG * EXP (XOUTS(IX) / CLTIOUG)              
            ELSE                                                                
              CTEMBSI(IX,IY) = CTIBOUL * EXP (XOUTS(IX) / CLTIOUL)              
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
      CALL RZERO(CTIGS,MAXNXS*(2*MAXNYS+1)) 
      CALL RZERO(CTEGS,MAXNXS*(2*MAXNYS+1))
C
C     SET SCALING FACTOR FOR TEMPERATURE GRADIENT FORCES
C
      TGSCAL = (1.6E-19)/(CRMI*1.673E-27) * QTIM *QTIM 
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
1140      CONTINUE
          NTEMP = NTEG + 1
          TGMULT (NTEMP) = 1.0
          TGPOS (NTEMP) = CL
          LEN2 = CL - QEDGES(IQXS(IX),1)
          DO 1150 IN = NTEG,1,-1
            TGPOS(NTEMP+NTEG-IN+1) = CL+(1.0-TMEG(IN,1)) * LEN2
            TGMULT(NTEMP+NTEG-IN+1) = TMEG(IN,2)
 1150     CONTINUE
          NTEMP = NTEMP +NTEG
          TGMULT(NTEMP+1) = TMEG(1,2)               
C
C         WRITE(6,*) 'NTEMP:',NTEMP 
C         WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C         WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C
          CTEMBS(IX,0) = CTEMBS(IX,0) * TGMULT(1)   
          DO 1160 IY = 1,NYS
            IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
            IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
               MULT = TGMULT(IN)
            ELSE
               MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >           * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
            ENDIF
            CTEMBS(IX,IY) = MULT * CTEMBS(IX,IY)
            CTEMBS(IX,IY-NYS-1) = MULT * CTEMBS(IX,IY-NYS-1)
1160      CONTINUE 
C         WRITE(6,'(A,I5,2X,100F7.2)') 'CTEMBS: ',IX,
C    >                  (CTEMBS(IX,IY),IY=1,NYS)    
C
C     END OF IX LOOP
C
1130    CONTINUE 
C
C       CALCULATE THE TEMPERATURE GRADIENT FOR A BIN.
C       IT WILL BE THE AVERAGE OF THE GRADIENTS FROM THE CURRENT 
C       BIN CENTRE TO THE ADJACENT BIN CENTRES. FOR BINS AT THE  
C       ENDS THE GRADIENT IS ONLY BASED ON THE ONE NEIGHBOURING
C       BIN.
C      
        DO 1180 IX = 1,IXOUT
C
C       DSTEP CONVERTS THE FORCE GRADIENT QUANTITIES TO 
C       A DISTANCE STEP FOR 1 TIMESTEP IN LIM (AT A GIVEN X POSITION)
C       THE QUANTITY IS deltaT^2 SINCE ONE dT COMES FROM THE VELOCITY 
C       EQUATION AND THE SECOND IN CONVERTING THE dV TO A DISTANCE.
C
          DSTEP = TGSCAL * QS(IQXS(IX)) * QS(IQXS(IX))
          DO 1190 IY = -NYS+1,-2
            GRAD1 = (CTEMBS(IX,IY+1)-CTEMBS(IX,IY)) 
     >              / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBS(IX,IY)-CTEMBS(IX,IY-1)) 
     >              / ( YOUTS(IY) -YOUTS(IY-1))
            CTEGS(IX,IY) = DSTEP * (GRAD1+GRAD2)/2.0
1190      CONTINUE          
          DO 1200 IY = 2,NYS-1
            GRAD1 = (CTEMBS(IX,IY+1)-CTEMBS(IX,IY)) 
     >              / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBS(IX,IY)-CTEMBS(IX,IY-1)) 
     >              / ( YOUTS(IY) -YOUTS(IY-1))
            CTEGS(IX,IY) = DSTEP * (GRAD1+GRAD2)/2.0
1200      CONTINUE          
C
C         SPECIAL CASES - THE END POINTS 
C                         THE NULL POINT 
C
          CTEGS(IX,-NYS) = (CTEMBS(IX,-NYS+1)-CTEMBS(IX,-NYS))   
     >            / (YOUTS(-NYS+1)+YOUTS(-NYS)) * DSTEP
          CTEGS(IX,NYS) = (CTEMBS(IX,NYS)-CTEMBS(IX,NYS-1))   
     >            / (YOUTS(NYS)-YOUTS(NYS-1)) * 1.6E-19 *DSTEP
          CTEGS(IX,1) = (CTEMBS(IX,2)-CTEMBS(IX,1))   
     >            / (YOUTS(2)-YOUTS(1)) * DSTEP
          CTEGS(IX,-1) = (CTEMBS(IX,-1)-CTEMBS(IX,-2))   
     >            / (YOUTS(-1)-YOUTS(-2)) * DSTEP
          CTEGS(IX,0) = 0.0
1180    CONTINUE
      ENDIF
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
1340      CONTINUE
          NTEMP = NTIG + 1
          TGMULT (NTEMP) = 1.0
          TGPOS (NTEMP) = CL
          LEN2 = CL - QEDGES(IQXS(IX),1)
          DO 1350 IN = NTIG,1,-1
            TGPOS(NTEMP+NTIG-IN+1) = CL+(1.0-TMIG(IN,1)) * LEN2
            TGMULT(NTEMP+NTIG-IN+1) = TMIG(IN,2)
 1350     CONTINUE
          NTEMP = NTEMP +NTIG
          TGMULT(NTEMP+1) = TMIG(1,2)               
C
C         WRITE(6,*) 'NTEMP:',NTEMP 
C         WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C         WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C
          CTEMBS(IX,0) = CTEMBS(IX,0) * TGMULT(1)   
          DO 1360 IY = 1,NYS
            IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
            IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
               MULT = TGMULT(IN)
            ELSE
               MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >           * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
            ENDIF
            CTEMBSI(IX,IY) = MULT * CTEMBSI(IX,IY)
            CTEMBSI(IX,IY-NYS-1) = MULT * CTEMBSI(IX,IY-NYS-1)
1360      CONTINUE 
C
C     END OF IX LOOP
C
1330    CONTINUE 
C
C       CALCULATE THE TEMPERATURE GRADIENT FOR A BIN.
C       IT WILL BE THE AVERAGE OF THE GRADIENTS FROM THE CURRENT 
C       BIN CENTRE TO THE ADJACENT BIN CENTRES. FOR BINS AT THE  
C       ENDS THE GRADIENT IS ONLY BASED ON THE ONE NEIGHBOURING
C       BIN.
C      
        DO 1380 IX = 1,IXOUT
          DSTEP = TGSCAL *  QS(IQXS(IX)) * QS(IQXS(IX))
          DO 1390 IY = -NYS+1,-2
            GRAD1 = (CTEMBSI(IX,IY+1)-CTEMBSI(IX,IY)) 
     >              / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBSI(IX,IY)-CTEMBSI(IX,IY-1)) 
     >              / ( YOUTS(IY) -YOUTS(IY-1))
            CTIGS(IX,IY) = DSTEP * (GRAD1+GRAD2)/2.0
1390      CONTINUE          
          DO 1400 IY = 2,NYS-1
            GRAD1 = (CTEMBSI(IX,IY+1)-CTEMBSI(IX,IY)) 
     >              / ( YOUTS(IY+1) -YOUTS(IY))
            GRAD2 = (CTEMBSI(IX,IY)-CTEMBSI(IX,IY-1)) 
     >              / ( YOUTS(IY) -YOUTS(IY-1))
            CTIGS(IX,IY) = DSTEP * (GRAD1+GRAD2)/2.0
1400      CONTINUE          
C
C         SPECIAL CASES - THE END POINTS 
C                         THE NULL POINT 
C
          CTIGS(IX,-NYS) =(CTEMBSI(IX,-NYS+1)-CTEMBSI(IX,-NYS))   
     >           / (YOUTS(-NYS+1)+YOUTS(-NYS)) *DSTEP
          CTIGS(IX,NYS) = (CTEMBSI(IX,NYS)-CTEMBSI(IX,NYS-1))   
     >           / (YOUTS(NYS)-YOUTS(NYS-1))  * DSTEP
          CTIGS(IX,1) = (CTEMBSI(IX,2)-CTEMBSI(IX,1))   
     >           / (YOUTS(2)-YOUTS(1)) * 1.6E-19 * DSTEP
          CTIGS(IX,-1) = (CTEMBSI(IX,-1)-CTEMBSI(IX,-2))   
     >           / (YOUTS(-1)-YOUTS(-2)) * DSTEP
          CTIGS(IX,0) = 0.0
1380    CONTINUE
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
C         WRITE(6,*) 'NTEMP:',NTEMP 
C         WRITE(6,*) 'POS:',(TGPOS(IN),IN=1,NTEMP+1)
C         WRITE(6,*) 'MULT:',(TGMULT(IN),IN=1,NTEMP+1)  
C
          CRNBS(IX,0) = CRNBS(IX,0) * TGMULT(1)   
          DO IY = 1,NYS
            IN = IPOS(YOUTS(IY),TGPOS,NTEMP)
            IF (IN.EQ.1.OR.IN.EQ.NTEMP+1) THEN
               MULT = TGMULT(IN)
            ELSE
               MULT=(TGMULT(IN)-TGMULT(IN-1))/(TGPOS(IN)-TGPOS(IN-1))  
     >           * (YOUTS(IY)-TGPOS(IN-1)) + TGMULT(IN-1)
            ENDIF
            CRNBS(IX,IY) = MULT * CRNBS(IX,IY)
            CRNBS(IX,IY-NYS-1) = MULT * CRNBS(IX,IY-NYS-1)
          end do
C         WRITE(6,'(A,I5,2X,100F7.2)') 'CRNBS: ',IX,
C    >                  (CRNBS(IX,IY),IY=1,NYS)    
C
C         END OF IX LOOP
C
        end do
c
      endif 
C
C                                                                               
C-----------------------------------------------------------------------        
C     NORMAL AND ERROR RETURN POINTS ...                                        
C-----------------------------------------------------------------------        
C                                                                               
C     DO 2000 IQX = 1-NQXSO, NQXSI, 20                                          
C       WRITE (6,'('' IQX='',I5,''  X='',1P,G11.4,''  TB,NB='',4G11.4)')        
C    >    IQX,QXS(IQX),(QTEMBS(IQX,J),QRNBS(IQX,J),J=1,2)                       
C2000 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
