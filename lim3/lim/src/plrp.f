      SUBROUTINE PLRP (NIZS,PLAMS,PIZS,NLS,CION)                  
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
      IMPLICIT    none
      INCLUDE     'params'                                                      
C     INCLUDE     (PARAMS)                                                      
c      INCLUDE     'dynam1'                                                      
C     INCLUDE     (DYNAM1)                                                      
c      INCLUDE     'dynam3'                                                      
C     INCLUDE     (DYNAM3)                                                      
      INTEGER     NIZS,PIZS(MAXNLS),NLS,CION                                    
      REAL        PLAMS(MAXNLS)                           
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  PLRP:  THIS ROUTINE DERIVED FROM NOTES 81,84. EXTRACTS PARTICULAR*        
C  *  LINE RADIATION DATA FROM TABLES BELOW AND APPLIES TO THE         *        
C  *  NUMBER DENSITY ION DISTRIBUTION PASSED IN AS DDLIMS ARRAY.       *        
C  *                                                                   *        
C  *  INPUT:  NIZS     NUMBER OF IONISATION STATES USED                *        
C  *          DDLIMS   NUMBER DENSITY ION DISTRIBUTION FROM /DYNAM1/   *        
C  *       OR TIZS     IONISATION DENSITY DISTRIBUTION FROM /DYNAM3/   *        
C  *          DDLIM3   NUMBER DENSITY ION 3D DISTB'N   FROM /DYNAM1/   *        
C  *       OR TIZ3     IONISATION DENSITY 3D DISTB'N   FROM /DYNAM3/   *        
C  *                                                                   *        
C  *  OUTPUT: PLAMS    WAVELENGTH OF EACH LINE                         *        
C  *          PIZS     IONISATION STATE EACH LINE REFERS TO            *        
C  *          NLS      NUMBER OF LINES RETURNED                        *        
C  *                                                                   *        
C  *                                      C.M.FARRELL   NOVEMBER 1987  *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER        L,IX,IY,IZ,IT,IP,JY                                        
      REAL           ETAS(MAXNXS),VAL                                           
c      INCLUDE        'comt2'                                                    
C     INCLUDE        (COMT2)                                                    
      INCLUDE        'comxyt'                                                   
C     INCLUDE        (COMXYT)                                                   
C                                                                               
C---- PARTICULAR LINE RADIATION DATA                                            
C---- PARAMETERS NL= NUMBER OF LINES OF DATA IN FOLLOWING TABLE                 
C----            NT= NUMBER OF TEMPERATURE POINTS DATA GIVEN FOR                
C---- NOTE IBM ONLY ALLOWS 19 CONTINUATION LINES IN ANY ONE STATEMENT.          
C---- EACH LINE HAS AT NO., IZ STATE, LAMBDA, 5 VALUES, FLAG                    
C---- THE FLAG INDICATES 0: IGNORE THIS LINE  1: USE THIS LINE                  
C                                                                               
      INTEGER     NL,NT                                                         
c slmod begin - N PEC data
      PARAMETER   (NL=37,NT=5)                                                  
c      PARAMETER   (NL=36,NT=5)                                                  
c slmod end
      REAL        RD(9*NL),LLAMS(NL),LETAS(NT,NL),LTBS(NT)                      
      INTEGER     LIONS(NL),LIZS(NL),LFLAGS(NL)                                 
      DATA        LTBS   / 10.,20.,50.,100.,200. /                              
C                                                                               
C         *************** HELIUM **********************                         
C         NOTE 150: THERE IS NO HE0 DATA, USE HE+ DATA                          
C         NOTE 186: NOW THERE IS, LAMBDA FROM NOTE 187.                         
C         NOTE 303: NEW HE0 DATA                                                
C                                                                               
      DATA        (RD(L), L=1,19*9) /                                           
     >     2, 0, 587.5, 0.01, .054, 0.22, 0.46, 0.78, 0                         
     >,    2, 0, 667.8, 24. , 52. , 127., 211., 330., 1                         
     >,    2, 1, 468.6, 7.9 , 23. , 50. , 74. , 90.,  1                         
C         
C         *************** BERYLLIUM ******************* 
C
     >,    4, 0, 440.7,359.0,393.0,423.0,437.0,455.0, 0
     >,    4, 0, 825.4, 12.8, 15.1, 17.4, 19.3, 21.5, 1
     >,    4, 1, 527.1, 10.2, 15.9, 23.3, 28.8, 34.0, 1
     >,    4, 1, 436.0,  6.2, 11.4, 20.3, 28.5, 37.0, 0
     >,    4, 2, 372.0,.0056,0.066, 0.58, 2.04, 6.07, 1 
     >,    4, 3, 465.9, 6.92, 17.4, 45.7, 83.6, 140.0,1 
C                                                                               
C         *************** CARBON **********************                         
C                                                                               
     >,    6, 0, 909.5, 7.4 , 16. , 35. , 56. , 85.,  1                         
     >,    6, 1, 657.8, 4.16, 11.3, 26.6, 38.1, 65.0, 1
     >,    6, 1, 514.5, 3.57, 8.42, 17.1, 23.9, 34.0, 0
     >,    6, 1,90.409, 0.14, 0.31, 0.55, 0.69, 0.77, 0                          
     >,    6, 1, 687.2, 0.85, 1.61, 2.54, 3.01, 3.22, 0 
     >,    6, 2, 569.6, 44.1,120.0,245.0,342.0,460.0, 0                         
     >,    6, 2, 464.7, 0.31, 1.21, 3.50, 5.68, 9.2,  1                         
     >,    6, 2, 977.0,.0014, .013, .059, 0.11, 0.16, 0                         
     >,    6, 2, 459.6, .028, 0.14, 0.44, 0.72, 0.97, 0                         
     >,    6, 3, 31.24, .046, 0.26, 0.83 , 1.3 ,1.64, 1  /                      
c slmod begin - N PEC data - May 5, 97
C                                                                               
C         *************** NITROGEN ********************                         
C                                                                               
      DATA        (RD(L), L=19*9+1, 20*9) /                                     
     >     7, 1, 463.7, 30.0, 83.5,190.0,269.0,326.0, 1  /
C                                                                               
C         *************** OXYGEN **********************                         
C         NOTE 287: THERE IS NO O0 DATA, USE C0 DATA                            
C                                                                               
      DATA        (RD(L), L=20*9+1, 9*NL) /                                     
c
c      DATA        (RD(L), L=19*9+1, 9*NL) /                                     
c slmod end
     >     8, 0, 909.5, 7.4 , 16. , 35. , 56. , 85.,  1                         
     >,    8, 1, 374.9, 5.88, 15.4, 34.8, 53.1, 81.0, 0                         
     >,    8, 1, 435.1, 1.3 , 4.6 , 15. , 27. , 42.,  0
     >,    8, 1, 441.5, 3.29, 8.49, 18.9, 28.2, 43.0, 1
     >,    8, 1, 391.2, 2.1 , 8.2 , 27. , 50. , 78.,  0                         
     >,    8, 1, 397.3, 4.34, 10.8, 23.7, 34.8, 54.0, 0                         
     >,    8, 2, 376.0, 2.44, 10.4, 33.0, 57.8, 96.0, 1                         
     >,    8, 2, 559.2, 3.91, 15.4, 40.6, 57.6, 93.0, 0                         
     >,    8, 2, 370.3, 3.0 , 15. , 57. , 172., 188., 0                         
     >,    8, 2, 703.3,.0074, .093, 0.56, 1.2 , 1.9,  0                         
     >,    8, 3, 789.4,1.9E-4,.0084,0.12, 0.32, 0.58, 1                         
     >,    8, 3, 233.6, .012, .087, 0.40, 0.79, 1.1,  0                         
     >,    8, 4, 629.7,1.5E-6,3.E-4,.0089,.033, 0.07, 1                         
     >,    8, 4, 760.4,2.1E-6,5.E-4,.017, .066, 0.14, 0                         
     >,    8, 5,1032.0,3.4E-8,3.3E-5,.0028,.015,0.04, 1                         
C                                                                               
C         ************** CHLORINE *********************                         
C                                                                               
     >,   17, 1, 489.7, 39. , 92. , 205., 312., 430., 1                         
C                                                                               
C         ************** CHROMIUM *********************                         
C                                                                               
     >,   24, 0, 425.4, 0.34, 0.62, 1.  , 1.3 , 1.6,  1  /                      
C                                                                               
C-----------------------------------------------------------------------        
C     INITIALISE LOCAL ARRAYS.  DIAGNOSTICS TO CHANNEL 6                        
C     SET PHOTON EFFICIENCY TO CONSTANT VALUES WHEN EXTRAPOLATING BEYOND        
C     LAST SPECIFIED POINTS.                                                    
C-----------------------------------------------------------------------        
C                                                                               
      WRITE (6,9000)                                                            
      DO 20 L = 1, NL                                                           
        LIONS(L) = NINT (RD(9*(L-1)+1))                                         
        LIZS(L)  = NINT (RD(9*(L-1)+2))                                         
        LLAMS(L) = RD(9*(L-1)+3)                                                
        DO 10 IT = 1, 5                                                         
          LETAS(IT,L) = RD(9*(L-1)+3+IT)                                        
   10   CONTINUE                                                                
        LFLAGS(L) = NINT(RD(9*(L-1)+9))                                         
        IF (CION.EQ.LIONS(L))                                                   
     >    WRITE (6,9001) LIONS(L),LIZS(L),LLAMS(L),(LETAS(IT,L),IT=1,5)         
   20 CONTINUE                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C     FOR EACH IONISATION STATE CHECK EACH ROW OF DATA TABLE UNTIL ...          
C-----------------------------------------------------------------------        
C                                                                               
      NLS = 0                                                                   
      DO 500 IZ = 0, NIZS                                                       
        DO 400 L = 1, NL                                                        
          IF (CION.EQ.LIONS(L) .AND. IZ.EQ.LIZS(L)) THEN                        
            IF (LFLAGS(L).EQ.0) THEN                                            
              WRITE (6,9004) LLAMS(L)                                           
              GOTO 400                                                          
            ENDIF                                                               
            IF (NLS.EQ.MAXNLS) THEN                                             
              WRITE (6,9003) LLAMS(L),MAXNLS                                    
              GOTO 400                                                          
            ENDIF                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     ... WE FIND A MATCH.  INCREMENT COUNT, STORE LAMBDA/STATE DETAILS         
C     AND SET UP FITTING DATA ON TEMPERATURE VS ETA VALUES USING THE            
C     FITTER ROUTINE, AND THEN ...                                              
C-----------------------------------------------------------------------        
C                                                                               
            NLS = NLS + 1                                                       
            PLAMS(NLS) = LLAMS(L)                                               
            PIZS(NLS)  = LIZS(L)                                                
C                                                                               
C-----------------------------------------------------------------------        
C     ... FOR EACH X POINT, INTERPOLATE ON TEMPERATURE TO OBTAIN AN             
C     EXACT ETA FROM THE CUBIC SPLINE, USING HARWELL ROUTINE                    
C     TG01B (FLAG,N,X,F,F',X0) WITHIN RANGE 10 < TEMP < 200, OR                 
C     USING BORDER VALUE BEYOND THIS RANGE, AND ...                             
C-----------------------------------------------------------------------        
C                                                                               
            DO 310 IY = -NYS, NYS                                               
             JY = IABS(IY)                                                      
             CALL FITTER (5,LTBS,LETAS(1,L),                                    
     >                    NXS,CTEMBS(1,IY),ETAS,'SPLINE')                       
             DO 300 IX = 1, NXS                                                 
C             IF (10*(IY/10).EQ.IABS(IY).AND.10*(IX/10).EQ.IX)                  
C    >          WRITE (6,9002) IZ,IY,L,NLS,IX,CTEMBS(IX,IY),                    
C    >                         ETAS(IX),CFIZS(IX,IY,IZ)                         
C                                                                               
C-----------------------------------------------------------------------        
C     ... FOR EACH Y POINT CALCULATE PLRP VALUE FROM                            
C     (A)  ION DENSITY, IONISATION TIME AND PHOTON EFFICIENCY,  OR              
C     (B)  IONISATION DENSITY AND PHOTON EFFICIENCY                             
C          (COMMENT OUT ONE OPTION)                                             
C-----------------------------------------------------------------------        
C                                                                               
              VAL = SNGL(DDLIMS(IX,IY,IZ)) / CFIZS(IX,IY,IZ) / ETAS(IX)      
C#              VAL = TIZS  (IX,IY,IZ) / ETAS(IX)                              
              PLRPS(IX,IY,NLS) = VAL                                            
              IF (JY.LE.NY3D) THEN                                              
                DO 210 IP = -MAXNPS, MAXNPS                                     
                  VAL=SNGL(DDLIM3(IX,IY,IZ,IP))/CFIZS(IX,IY,IZ)/ETAS(IX)   
C#                  VAL = TIZ3  (IX,IY,IZ,IP) / ETAS(IX)                        
                  PLRP3(IX,IY,NLS,IP) = VAL                                     
  210           CONTINUE                                                        
              ENDIF                                                             
  300        CONTINUE                                                           
  310       CONTINUE                                                            

c slmod begin 
            WRITE(63,*) ' '
            WRITE(63,*) 'Not sure this output is right...'
            WRITE(63,*) ' '

            IY = 1

            WRITE(63,'(2A4,A10,A18)')
     +        'IX','IY','Te','Photon Efficiency'

            DO IX = 1, NXS
             WRITE(63,'(2I4,F10.4,3E12.4)') 
     +        IX,IY,CTEMBSI(IX,1),
c     +        PLRP3(IX,IY,0,1) * CFIZS(IX,IY,0) / SNGL(DDLIMS(IX,IY,0)),
     +        PLRP3(IX,IY,1,1) * CFIZS(IX,IY,1) / SNGL(DDLIMS(IX,IY,1)),
     +        PLRP3(IX,IY,2,1) * CFIZS(IX,IY,2) / SNGL(DDLIMS(IX,IY,2)),
     +        PLRP3(IX,IY,3,1) * CFIZS(IX,IY,3) / SNGL(DDLIMS(IX,IY,3))
            ENDDO

c            WRITE(0,*) 'DEBUG: Photon efficiency print out complete'
c slmod end
C                                                                               
C-----------------------------------------------------------------------        
C     EXTRA SECTION FOR NEUTRALS.  CONVERT PREVIOUS ENTRY INTO ENTRY            
C     FOR PRIMARY NEUTRALS, AND CALCULATE ADDITIONAL ENTRY FOR TOTAL            
C     NEUTRALS.                                                                 
C-----------------------------------------------------------------------        
C                                                                               
            IF (IZ.EQ.0) THEN                                                   
              NLS = NLS + 1                                                     
              PLAMS(NLS) = LLAMS(L)                                             
              PIZS(NLS)  = 0                                                    
              PIZS(NLS-1)= -1                                                   
              DO 330 IY = -NYS, NYS                                             
               JY = IABS (IY)                                                   
               DO 320 IX = 1, NXS                                               
                PLRPS(IX,IY,NLS) = PLRPS(IX,IY,NLS-1)                           
                VAL = 0.0                                                       
                IF (DDLIMS(IX,IY,0).GT.0.0D0)                              
     >            VAL = SNGL (DDLIMS(IX,IY,-1) / DDLIMS(IX,IY,0))           
C#                IF (TIZS  (IX,IY,0).GT.0.0)                                 
C#     >            VAL = TIZS(IX,IY,-1) / TIZS(IX,IY,0)                     
                PLRPS(IX,IY,NLS-1) = PLRPS(IX,IY,NLS-1) * VAL                   
                IF (JY.LE.NY3D) THEN                                            
                 DO 311 IP = -MAXNPS, MAXNPS                                    
                  PLRP3(IX,IY,NLS,IP) = PLRP3(IX,IY,NLS-1,IP)                   
                  VAL = 0.0                                                     
                  IF (DDLIM3(IX,IY,0,IP).GT.0.0D0)                           
     >              VAL = SNGL(DDLIM3(IX,IY,-1,IP)/DDLIM3(IX,IY,0,IP))     
C#                  IF (TIZ3  (IX,IY,0,IP).GT.0.0)                             
C#     >              VAL = TIZ3(IX,IY,-1,IP) / TIZ3(IX,IY,0,IP)               
                  PLRP3(IX,IY,NLS-1,IP) = PLRP3(IX,IY,NLS-1,IP) * VAL           
  311            CONTINUE                                                       
                ENDIF                                                           
  320          CONTINUE                                                         
  330         CONTINUE                                                          
            ENDIF                                                               
C-----------------------------------------------------------------------        
          ENDIF                                                                 
  400   CONTINUE                                                                
  500 CONTINUE                                                                  
C
C     INDICATE IN THE OUTPUT FILE WHICH METHOD WAS USED TO CALCULATE
C     THE PLRP'S.
C
      CALL PRB
      CALL PRC('  NOTE: PLRPS WERE CALCULATED BASED ON ION DENSITY')  
C#      CALL PRC('  NOTE: PLRPS WERE CALCULATED BASED ON IONIZATION')
      CALL PRB
C
      RETURN                                                                    
C                                                                               
 9000 FORMAT(/1X,'PLRP:  PHOTON EFFICIENCY DATA',//1X,                          
     >  '      SPECIES          LINE                  ELECTRON TEMPER',         
     >  'ATURE  (EV)',/1X,                                                      
     >  '    ATOM   CHARGE    WAVELENGTH      10        20        50 ',         
     >  '      100       200',/1X,81('-'))                                      
 9001 FORMAT(1X,2I8,3X,F9.1,6X,1P,5G10.3)                                       
 9002 FORMAT(1X,'IZ=',I2,' IY=',I4,' L=',I2,' NLS=',I2,' IX=',I3,               
     >  ' TEMP=',F9.2,' ETA=',F9.2,' CFIZS=',1P,G9.2)                           
 9003 FORMAT(1X,'PLRP: LINE LAMBDA =',F7.1,' IGNORED - MAX OF',I3,              
     >  ' LINES REACHED.')                                                      
 9004 FORMAT(1X,'PLRP: LINE LAMBDA =',F7.1,' IS NOT USED AT PRESENT.')          
      END                                                                       
c
c
c
      SUBROUTINE SPECTEMP (PLAMS,PIZS,NLS)      
      use mod_dynam1
      use mod_dynam3
      use mod_comt2
C
C     THIS SUBROUTINE CALCULATES THE SPECTROSCOPIC TEMPERATURE
C     OF THE PLASMA BY WEIGHTING THE BIN TEMPERATURES BY THE VALUE
C     OF THE PARTICULAR LINE RADIATION IN THAT BIN. THIS VALUE IS 
C     CALCULATED INBOARD AND OUTBOARD OF A SPECIFIED X VALUE.
C     IT IS ALSO CALCULATED FOR THE OVERALL PLASMA AND AS AN
C     ALTERNATIVE AS A FUNCTION OF THE BACKGROUND ELECTRON 
C     TEMPERATURE INSTEAD OF THE ION TEMPERATURE.
C     THESE RESULTS ARE THEN PRINTED OUT.
C
C     THIS ROUTINE ALSO CALCULATES THE SPECTROSCOPICALLY WEIGHTED
C     TEMPERATURES AS A FUNCTION OF X AND Y.  APRIL 6,1990
C
C     D. ELDER , 1990, MARCH 27
C
      IMPLICIT none
C
      INCLUDE 'params'
C
c      INCLUDE 'dynam1'
C
c      INCLUDE 'dynam3'
C
      INCLUDE 'comxyt'
C
c      INCLUDE 'comt2'
C
      REAL PLAMS(MAXNLS)
      INTEGER NLS,PIZS(MAXNLS)
C
C     LOCAL VARIABLES
C
C     4 ELEMENTS SHALLOW, DEEP, TOTAL, AND ELECTRON
      REAL SCT(MAXNLS,4)
C     2 ELEMENTS SHALLOW AND DEEP 
      REAL INTPLRP(MAXNLS,3),FACTA,FACTB,SUM1,SUM2,SUM3,SUM4
      INTEGER I,J,K,IX,IPOS,ABSJ
C
C
      CALL RZERO (SCT,MAXNLS*4)
      CALL RZERO (INTPLRP,MAXNLS*3)
      CALL RZERO (SCTXS,MAXNXS*MAXNLS)
      CALL RZERO (SCTYS,(2*MAXNYS+1)*MAXNLS)
      SUM1 = 0.0
      SUM2 = 0.0
      SUM3 = 0.0
      SUM4 = 0.0
C 
      IX = IPOS(CXDEEP,XS,NXS-1)
      WRITE(6,*) 'SPECTEMP-IXDEEP:',IX 
      WRITE(6,1000) 
C
      DO 700 K = 1,NLS
        IF (PIZS(K).NE.0.AND.PIZS(K).NE.-1) THEN
          DO 150 J = -NYS,NYS
            ABSJ = ABS(J)

            ! jdemod ywids has no J=0 element - need to avoid executing this
            !        code then
            if (j.ne.0) then 

            DO 100 I = 1, IX
              FACTA = XWIDS(I)*XCYLS(I)*YWIDS(ABSJ)*DELPS(I,ABSJ)
              SCT(K,1)=SCT(K,1)+SNGL(DDTS(I,J,PIZS(K)))
     >                   *PLRPS(I,J,K)*FACTA 
              SCT(K,4) = SCT(K,4)+CTEMBS(I,J)*PLRPS(I,J,K)*FACTA 
              INTPLRP(K,1) = INTPLRP(K,1) + PLRPS(I,J,K) *FACTA
100         CONTINUE 
            DO 50 I = IX+1,NXS
              FACTB = XWIDS(I)*XCYLS(I)*YWIDS(ABSJ)*DELPS(I,ABSJ)
              SCT(K,2)=SCT(K,2)+SNGL(DDTS(I,J,PIZS(K)))
     >                    *PLRPS(I,J,K)*FACTB
              SCT(K,4) = SCT(K,4)+CTEMBS(I,J)*PLRPS(I,J,K)*FACTB
              INTPLRP(K,2) = INTPLRP(K,2) + PLRPS(I,J,K)*FACTB
 50         CONTINUE

            endif

 150      CONTINUE
          INTPLRP(K,3) = INTPLRP(K,1)+INTPLRP(K,2) 
          SCT(K,3) = SCT(K,1)+SCT(K,2)
C
C
C         WRITE(6,*) 'SCTS :',SCT(K,1),SCT(K,2),SCT(K,3),SCT(K,4)
C         WRITE(6,*) 'IPLRP:',INTPLRP(K,1),INTPLRP(K,2),
C    >                INTPLRP(K,3)
          IF (INTPLRP(K,1).GT.0.0) THEN
            SCT(K,1) = SCT(K,1) / INTPLRP(K,1)
          ENDIF
          IF (INTPLRP(K,2).GT.0.0) THEN
            SCT(K,2) = SCT(K,2) / INTPLRP(K,2)
          ENDIF 
          IF (INTPLRP(K,3).GT.0.0) THEN
            SCT(K,3)=SCT(K,3)/INTPLRP(K,3)
            SCT(K,4)=SCT(K,4)/INTPLRP(K,3)
          ENDIF 
C         WRITE(6,*) 'SCTF :',SCT(K,1),SCT(K,2),SCT(K,3),SCT(K,4)
C
C
C     CALCULATE THE AVERAGE TEMPERATURES AS A FUNCTION OF X AND Y
C     THE DELPS ARRAY HAS BEEN DROPPED FROM THE NORMALISATION
C     CORRECTION BECAUSE THE FEATURE IS NOT USED AND SHOULD BE 
C     COMPLETELY EXCISED IN TIME. (THE VALUE IS ALWAYS 1.0)
C
C     THE ZERO INDEX Y-SLOT IN THE ARRAYS FOR AVERAGE TEMPERTATURE 
C     CREATES A PROBLEM SINCE IN THE CODE IT IS NOT USED BUT HERE  
C     IT NEEDS TO BE SUMMED ACCROSS THE EASIEST SOLUTION IS TO 
C     STORE IN THIS VARIABLE SCTYS(0,?) THE AVERAGE OF
C     SCTYS(1,?) AND SCTYS(-1,?) 
C
C
          DO 450 J = 1,NYS
            SUM1 = 0.0
            SUM2 = 0.0
            SUM3 = 0.0
            SUM4 = 0.0
            DO 500 I = 1, NXS
              FACTA = XWIDS(I)*XCYLS(I)
              SUM1 = SUM1 + SNGL(DDTS(I,J,PIZS(K)))
     >               *PLRPS(I,J,K)*FACTA 
              SUM2 = SUM2 + PLRPS(I,J,K) *FACTA
              SUM3 = SUM3 + SNGL(DDTS(I,-J,PIZS(K)))
     >               *PLRPS(I,-J,K)*FACTA 
              SUM4 = SUM4 + PLRPS(I,-J,K) *FACTA
 500        CONTINUE 
            IF (SUM2.GT.0.0) SCTYS(J,K) = SUM1 / SUM2
            IF (SUM4.GT.0.0) SCTYS(-J,K) = SUM3 / SUM4           
 450      CONTINUE
          SCTYS(0,K) = (SCTYS(1,K)+SCTYS(-1,K))/2.0
          DO 550 I = 1,NXS
            SUM1 = 0.0
            SUM2 = 0.0
            DO 600 J = 1, NYS
              FACTA = YWIDS(J)
              SUM1 = SUM1 + FACTA*( SNGL(DDTS(I,J,PIZS(K)))
     >         *PLRPS(I,J,K) + SNGL(DDTS(I,-J,PIZS(K)))
     >         *PLRPS(I,-J,K)) 
              SUM2 =SUM2+FACTA*(PLRPS(I,J,K)+PLRPS(I,-J,K))
 600        CONTINUE 
            IF (SUM2.GT.0.0) SCTXS(I,K) = SUM1 / SUM2
 550      CONTINUE
        ENDIF    
 700  CONTINUE
C
C
C
C     PRINT OUT SPECTROSCOPIC TEMPERATURES CALCULATED EARLIER
C
C
      CALL PRB
      WRITE(7,1000)  
      WRITE(7,1002) CXDEEP
      WRITE(7,1004)
      DO 400 K = 1,NLS
        IF (PIZS(K).NE.0.AND.PIZS(K).NE.-1) 
     >     WRITE(7,1006) PIZS(K),PLAMS(K),SCT(K,1),SCT(K,2),SCT(K,3),
     >                  SCT(K,4)
         
 400  CONTINUE
      CALL PRB
      RETURN
 1000 FORMAT('*********    SPECTROSCOPIC TEMPERATURES   *********')
 1002 FORMAT('             SHALLOW TEMPERATURES FOR X < ',G9.2)
 1004 FORMAT('STATE    LAMBDA    Tshallow      Tdeep      Ttotal',
     >       '     Telectron')
 1006 FORMAT(1X,I3,5X,F7.2,5X,G9.3,3X,G9.3,3X,G9.3,3X,G9.3)
      END 
