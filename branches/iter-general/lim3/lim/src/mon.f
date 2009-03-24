      SUBROUTINE MONINI (RTMIN, RTMAX, NIZS, CTBIN)                             
      IMPLICIT  none
      INTEGER   NIZS                                                            
      REAL      RTMIN,RTMAX,CTBIN                                               
C                                                                               
C***********************************************************************        
C                                                                               
C       INITIALIZES MONITORING VARIABLES IN COMMON BLOCKS COMMV/ZOMMV           
C                                                                               
C         ARGUMENTS -                                                           
C     RTMIN  : SCALED VALUE OF START TIME POINT                                 
C     RTMAX  : SCALED VALUE OF CUTOFF TIME                                      
C     NIZS   : NUMBER OF IONIZATION STATES BEING FOLLOWED                       
C     CTBIN  : TEMPERATURE BASE VALUE INBOARD.                                  
C                                                                               
C***********************************************************************        
C                                                                               
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'commv'                                                         
C     INCLUDE   (COMMV)                                                         
      INCLUDE   'zommv'                                                         
C     INCLUDE   (ZOMMV)                                                         
C                                                                               
      INTEGER   IZ                                                              
C                                                                               
C-----------------------------------------------------------------------        
C     INITIALISE ARRAYS                                                         
C-----------------------------------------------------------------------        
C                                                                               
      DO 100  IZ = 1, NIZS                                                      
         CICIZS(IZ) = 0.0                                                       
         CIFIZS(IZ) = RTMAX                                                     
         CILIZS(IZ) = RTMIN                                                     
         CISIZS(IZ) = 0.0                                                       
C                                                                               
         CICABS(IZ) = 0.0                                                       
         CIFABS(IZ) = RTMAX                                                     
         CILABS(IZ) = RTMIN                                                     
         CISABS(IZ) = 0.0                                                       
         CRTABS(IZ) = 0.0                                                       
         CICUTS(IZ) = 0.0                                                       
         CRTRCS(IZ) = 0.0                                                       
         CRVABS(IZ) = 0.0                                                       
         CRAVAV(IZ) = 0.0                                                       
         CTBS  (IZ) = 0.0                                                       
C                                                                               
         ZICABS(IZ) = 0.0                                                       
         ZIFABS(IZ) = RTMAX                                                     
         ZILABS(IZ) = RTMIN                                                     
         ZISABS(IZ) = 0.0                                                       
         ZRTABS(IZ) = 0.0                                                       
         ZRVABS(IZ) = 0.0                                                       
         ZRAVAV(IZ) = 0.0                                                       
         ZTBS  (IZ) = 0.0                                                       
C                                                                               
         CICRCS(IZ) = 0.0                                                       
         CIFRCS(IZ) = RTMAX                                                     
         CILRCS(IZ) = RTMIN                                                     
         CISRCS(IZ) = 0.0                                                       
         COUTS(IZ)  = 0.0                                                       
C
         CIY0IZ(IZ) = 0.0
         CI2LIZ(IZ) = 0.0
         ZIY0IZ(IZ) = 0.0
         ZI2LIZ(IZ) = 0.0  
C
         CTY0IZ(IZ) = 0.0
         CT2LIZ(IZ) = 0.0
         ZTY0IZ(IZ) = 0.0
         ZT2LIZ(IZ) = 0.0  
  100 CONTINUE                                                                  
C                                                                               
      CALL RZERO (CTEXS, 10)                                                    
      CTTT = 0.2 * CTBIN                                                        
C                                                                               
C-----------------------------------------------------------------------        
C     INITIALISE SCALARS                                                        
C-----------------------------------------------------------------------        
C                                                                               
      CICCOL = 0.0                                                              
      CICSTX = 0.0                                                              
      CICSTZ = 0.0                                                              
      CICRXA = 0.0                                                              
      CIFRXA = RTMAX                                                            
      CICY2L = 0.0                                                              
      CIFY2L = RTMAX                                                            
      CRXMIN = 0.0                                                              
      CICAY0 = 0.0                                                              
      CICAAW = 0.0                                                              
      CICA2L = 0.0                                                              
      ZICAY0 = 0.0                                                              
      ZICAAW = 0.0                                                              
      ZICA2L = 0.0                                                              
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE MONUP (ICODE,SPUTY)                                            
      IMPLICIT  none
      INTEGER   ICODE                                                           
      REAL      SPUTY                                                           
C                                                                               
C***********************************************************************        
C                                                                               
C    NOTE:  SOME CALLS TO MONUP REPLACED WITH INLINE CODE IN LIM3               
C           FOR SPEED IMPROVEMENTS                                              
C                                                                               
C                                                                               
C       RECORDS IN COMMON BLOCK COMMV AN EVENT HAPPENING TO AN ION.             
C                                                                               
C         ARGUMENTS -                                                           
C     ICODE  : CODE FOR EVENT HAPPENING                                         
C                                                                               
C         READ FROM COMMON COMTAU -                                             
C     CIZ    : CURRENT IONIZATION STATE                                         
C     CIST   : SCALED CURRENT TIME OF EVENT                                     
C     CX     : CURRENT X POSITION OF ION                                        
C     CIAB   : CODE FOR LIMITER SURFACE; 0 WALL,  +1 Y=0 FROM Y>0,              
C              +2 Y=2L,  -1 Y=0 FROM Y<0,  -2 Y=-2L                             
C     CTEMI  : CURRENT ION TEMPERATURE                                          
C     CTBIQX : BACKGROUND TEMPERATURE AT ABSORPTION POINT                       
C                                                                               
C         VALID ICODE VALUES              VARIABLES READ BY MONUP               
C                                    CIZ   CIST   CX   CIAB CTEMI CTBIQX        
C  1  ION MADE OUTBOARD ITERATION                 X                             
C  2  ION HAD A COLLISION                                                       
C  3  ION FURTHER IONIZED             X     X                                   
C  4  ION REACHED X = A                     X                 X                 
C  5  ION REACHED Y = 2L, X > 0             X                                   
C  6  ION ABSORBED FROM Y>0 REGION    X     X           X     X      X          
C  7  ION STILL IN PLASMA AT CUTOFF   X                       X                 
C  8  NEW ION BEING INJECTED                                                    
C  9  ION RECOMBINED                  X     X                                   
C 10  ION ABSORBED FROM Y<0 REGION    X     X           X     X      X          
C                                    CIZ   CIST   CX   CIAB CTEMI CTBIQX        
C                                                                               
C***********************************************************************        
C                                                                               
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INCLUDE   'comtau'                                                        
C     INCLUDE   (COMTAU)                                                        
      INCLUDE   'commv'                                                         
C     INCLUDE   (COMMV)                                                         
      INCLUDE   'zommv'                                                         
C     INCLUDE   (ZOMMV)                                                         
      INTEGER   IM                                                              
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS COMPLETED OUTBOARD STEP                   
C-----------------------------------------------------------------------        
C                                                                               
      IF      (ICODE .EQ. 1)  THEN                                              
         WRITE (6,'('' MONUP(1):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** COUTS(CIZ) = COUTS(CIZ) + SPUTY                                        
C ****** SUPERCEDED BY D.P. VERSION INLINE IN LIM3 ROUTINE !!!                  
C ****** IF (CX .LT. CRXMIN) CRXMIN = CX                                        
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS HAD A COLLISION                           
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 2)  THEN                                              
         WRITE (6,'('' MONUP(2):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** CICCOL = CICCOL + SPUTY                                                
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS IONIZED FURTHER                           
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 3)  THEN                                              
         WRITE (6,'('' MONUP(3):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** CICIZS(CIZ) = CICIZS(CIZ) + SPUTY                                      
C ****** IF (CIST .LT. CIFIZS(CIZ)) CIFIZS(CIZ) = CIST                          
C ****** IF (CIST .GT. CILIZS(CIZ)) CILIZS(CIZ) = CIST                          
C ****** CISIZS(CIZ) = CISIZS(CIZ) + CIST * SPUTY                               
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS REACHED X = A                             
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 4)  THEN                                              
         WRITE (6,'('' MONUP(4):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** IF (CFLRXA) THEN                                                       
C ******    CICRXA = CICRXA + SPUTY                                             
C ******    CISRXA = CISRXA + CIST                                              
C ******    CITRXA = CITRXA + CTEMI * SPUTY                                     
C ******    IF (CIST .LT. CIFRXA) CIFRXA = CIST                                 
C ******    CFLRXA = .FALSE.                                                    
C ****** ENDIF                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS REACHED Y = 2L, X > 0                     
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 5)  THEN                                              
         WRITE (6,'('' MONUP(5):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** IF (CFLY2L) THEN                                                       
C ******    CICY2L = CICY2L + SPUTY                                             
C ******    IF (CIST .LT. CIFY2L) CIFY2L = CIST                                 
C ******    CFLY2L = .FALSE.                                                    
C ****** ENDIF                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS BEEN ABSORBED FROM Y>0                    
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 6)  THEN                                              
         CICABS(CIZ) = CICABS(CIZ) + SPUTY                                      
         IF (CIST .LT. CIFABS(CIZ)) CIFABS(CIZ) = CIST                          
         IF (CIST .GT. CILABS(CIZ)) CILABS(CIZ) = CIST                          
         CISABS(CIZ) = CISABS(CIZ) + CIST * SPUTY                               
         CRTABS(CIZ) = CRTABS(CIZ) + CTEMI * SPUTY                              
         CRVABS(CIZ) = CRVABS(CIZ) + CVABS * SPUTY                              
         CRAVAV(CIZ) = CRAVAV(CIZ) + ABS(CVABS) * SPUTY                         
         CTBS  (CIZ) = CTBS  (CIZ) + CTBIQX * SPUTY                             
         IF     (CIAB.EQ.0) THEN                                                
            CICAAW = CICAAW + SPUTY                                             
         ELSEIF (CIAB.EQ.1) THEN                                                
            CICAY0 = CICAY0 + SPUTY                                             
            CIY0IZ(CIZ) = CIY0IZ(CIZ) + SPUTY
            CTY0IZ(CIZ) = CTY0IZ(CIZ) + CTEMI * SPUTY                              
         ELSEIF (CIAB.EQ.-2) THEN                                               
            ZICA2L = ZICA2L + SPUTY          
            ZI2LIZ(CIZ) = ZI2LIZ(CIZ) + SPUTY                                   
            ZT2LIZ(CIZ) = ZT2LIZ(CIZ) + CTEMI * SPUTY                              
         ENDIF                                                                  
         IM = MIN (INT(CTEMI/CTTT)+1, 10)                                       
         CTEXS(IM) = CTEXS(IM) + CTEMI * SPUTY                                  
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION REACHED CUTOFF TIME STILL IN PLASMA           
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 7)  THEN                                              
         CICUTS(CIZ) = CICUTS(CIZ) + SPUTY                                      
         CRTRCS(CIZ) = CRTRCS(CIZ) + CTEMI * SPUTY                              
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION JUST BEING INJECTED IN                        
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 8)  THEN                                              
         CFLRXA = .TRUE.                                                        
         CFLY2L = .TRUE.                                                        
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS RECOMBINED                                
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ. 9)  THEN                                              
         WRITE (6,'('' MONUP(9):  REPLACED WITH INLINE CODE !'')')              
         STOP                                                                   
C ****** CICRCS(CIZ) = CICRCS(CIZ) + SPUTY                                      
C ****** IF (CIST .LT. CIFRCS(CIZ)) CIFRCS(CIZ) = CIST                          
C ****** IF (CIST .GT. CILRCS(CIZ)) CILRCS(CIZ) = CIST                          
C ****** CISRCS(CIZ) = CISRCS(CIZ) + CIST * SPUTY                               
C                                                                               
C-----------------------------------------------------------------------        
C                     TEST IF ION HAS BEEN ABSORBED FROM Y<0                    
C-----------------------------------------------------------------------        
C                                                                               
      ELSE IF (ICODE .EQ.10)  THEN                                              
         ZICABS(CIZ) = ZICABS(CIZ) + SPUTY                                      
         IF (CIST .LT. ZIFABS(CIZ)) ZIFABS(CIZ) = CIST                          
         IF (CIST .GT. ZILABS(CIZ)) ZILABS(CIZ) = CIST                          
         ZISABS(CIZ) = ZISABS(CIZ) + CIST * SPUTY                               
         ZRTABS(CIZ) = ZRTABS(CIZ) + CTEMI * SPUTY                              
         ZRVABS(CIZ) = ZRVABS(CIZ) + CVABS * SPUTY                              
         ZRAVAV(CIZ) = ZRAVAV(CIZ) + ABS(CVABS) * SPUTY                         
         ZTBS  (CIZ) = ZTBS  (CIZ) + CTBIQX * SPUTY                             
         IF     (CIAB.EQ.0) THEN                                                
            ZICAAW = ZICAAW + SPUTY                                             
         ELSEIF (CIAB.EQ.-1) THEN                                               
            ZICAY0 = ZICAY0 + SPUTY 
            ZIY0IZ(CIZ) = ZIY0IZ(CIZ) + SPUTY                                            
            ZTY0IZ(CIZ) = ZTY0IZ(CIZ) + CTEMI * SPUTY                              
         ELSEIF (CIAB.EQ.2) THEN                                                
            CICA2L = CICA2L + SPUTY          
            CI2LIZ(CIZ) = CI2LIZ(CIZ) + SPUTY                                   
            CT2LIZ(CIZ) = CT2LIZ(CIZ) + CTEMI * SPUTY                              
         ENDIF                                                                  
         IM = MIN (INT(CTEMI/CTTT)+1, 10)                                       
         CTEXS(IM) = CTEXS(IM) + CTEMI * SPUTY                                  
C                                                                               
C-----------------------------------------------------------------------        
C                     IF ALL TESTS FAILED THEN IS AN ERROR                      
C-----------------------------------------------------------------------        
C                                                                               
      ELSE                                                                      
         WRITE(6, *)  '  ERROR - SUBROUTINE MONUP '                             
         WRITE(6, *)  '  INVALID EVENT CODE ', ICODE                            
      ENDIF                                                                     
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
      SUBROUTINE MONPRI (QTIM,FACT,VFLUID,NIZS,SDTZS,SDYZS,                     
     >                   STOTS,DOUTS,RIONS,CTBIN,CRMI)                          
      IMPLICIT  none
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INTEGER   NIZS                                                            
      REAL      QTIM,FACT,VFLUID,STOTS(20)                                      
      REAL      RIONS(MAXIZS),CTBIN,CRMI,SDTZS(MAXIZS),SDYZS(MAXIZS)            
      DOUBLE PRECISION DOUTS(MAXIZS)                                            
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
C     SDYZS  : AVERAGE Y STEPS DUE TO PARALLEL DIFFUSION                        
C     STOTS  : ZEFF AVERAGES, TOTAL POWER, ETC, ETC ...                         
C                                                                               
C***********************************************************************        
C                                                                               
      INCLUDE   'commv'                                                         
C     INCLUDE   (COMMV)                                                         
      INCLUDE   'zommv'                                                         
C     INCLUDE   (ZOMMV)                                                         
C                                                                               
      INTEGER   IZ,IM                                                           
      REAL      RCR,RCAB,RSAB,TEXT                                              
      REAL      RTR, RTAB, RTAV, RAVA, RMACH, RAVMCH, RENEGY                    
      REAL      RAVEGY, RZ0, RZ1, RZ2, DRZ, ZMACH, ZRZ, ZENEGY, RTBS            
      REAL      VEXIT,ZEXIT,MTSO,CICOUT                                         
C                                                                               
C-----------------------------------------------------------------------        
C          INITIALISATION                                                       
C-----------------------------------------------------------------------        
C                                                                               
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
      ZMACH = 0.0                                                               
      ZRZ = 0.0                                                                 
      ZENEGY = 0.0                                                              
      CICOUT = 0.0                                                              
      WRITE (6,'(A,(8G11.4))') ' MONPRI: SDTZS ',(SDTZS (IZ),IZ=1,NIZS)         
      WRITE (6,'(A,(8G11.4))') ' MONPRI: SDYZS ',(SDYZS (IZ),IZ=1,NIZS)         
      WRITE (6,'(A,(8G11.4))') ' MONPRI: DOUTS ',(DOUTS (IZ),IZ=1,NIZS)         
      WRITE (6,'(A,(8G11.4))') ' MONPRI: RIONS ',(RIONS (IZ),IZ=1,NIZS)         
      WRITE (6,'(A,(8G11.4))') ' MONPRI: STOTS ',(STOTS (IM),IM=1,20)           
C                                                                               
       DO  100  IZ = 1, NIZS                                                    
C                                                                               
C------- CALCULATE MEAN TIME SPENT OUTBOARD IN EACH CHARGE STATE.               
C                                                                               
         IF (RIONS(IZ).GT.0.0) THEN                                             
           MTSO = QTIM * SNGL(DOUTS(IZ)) / RIONS(IZ)                            
         ELSE                                                                   
           MTSO = 0.0                                                           
         ENDIF                                                                  
         CICOUT = CICOUT + QTIM * SNGL(DOUTS(IZ)) * FACT                        
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
         CALL PRR ('MEAN TIME SPENT OUTBOARD (S)                 ',             
     >                                 MTSO)                                    
         CALL PRI ('NUMBER OF IONS REACHING THIS CHARGE STATE    ',             
     >                             NINT(RIONS(IZ)))                             
         IF (CICUTS(IZ).GT.0.0) THEN                                            
           CALL PRB                                                             
           CALL PRI ('NUMBER OF IONS STILL IN PLASMA AT CUTOFF     ',           
     >                             NINT(CICUTS(IZ)))                            
           CALL PRR ('  MEAN TEMPERATURE  (EV)                     ',           
     >                             (CRTRCS(IZ)/CICUTS(IZ)))                     
         ENDIF                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C        PRINT ABSORPTION DETAILS FOR Y>0 AND Y<0 REGIONS                       
C-----------------------------------------------------------------------        
C                                                                               
         CALL PRB                                                               
         CALL PRI2 ('NUMBERS ABSORBED ON Y < 0 AND Y > 0 SIDES ',               
     >                   NINT(ZICABS(IZ)), NINT(CICABS(IZ)))                    
C
         CALL PRI2 ('NUMBERS ABSORBED ON Y < 0 FOR Y=0, Y= -2L ',
     >                   NINT(ZIY0IZ(IZ)),NINT(ZI2LIZ(IZ)))
         CALL PRI2 ('NUMBERS ABSORBED ON Y > 0 FOR Y=0, Y=  2L ',
     >                   NINT(CIY0IZ(IZ)),NINT(CI2LIZ(IZ)))
C                                                                               
         IF (CICABS(IZ).GT.0.0.OR.ZICABS(IZ).GT.0.0) THEN                       
           ! jdemod - these fixed constants are too large for R*4 - should use the parameter LO=1.0e-37 for R4
           !CICABS(IZ) = MAX (CICABS(IZ), 1.E-50)                                
           !ZICABS(IZ) = MAX (ZICABS(IZ), 1.E-50)                                 
           CICABS(IZ) = MAX (CICABS(IZ), LO)                                
           ZICABS(IZ) = MAX (ZICABS(IZ), LO)                                 
           CALL PRR2 ('  TIME FIRST ION ABSORBED  (S)          ',               
     >                   QTIM*ZIFABS(IZ), QTIM*CIFABS(IZ))                      
           CALL PRR2 ('  TIME LAST ION ABSORBED  (S)           ',               
     >                   QTIM*ZILABS(IZ), QTIM*CILABS(IZ))                      
           CALL PRR2 ('  MEAN ABSORPTION TIME  (S)             ',               
     >        QTIM*ZISABS(IZ)/ZICABS(IZ), QTIM*CISABS(IZ)/CICABS(IZ))           
           CALL PRR2 ('  MEAN TEMPERATURE AT ABSORPTION (EV)   ',               
     >             ZRTABS(IZ)/ZICABS(IZ), CRTABS(IZ)/CICABS(IZ))                
C
           IF (ZIY0IZ(IZ).GT.0.0) 
     >       CALL PRR ('  MEAN TEMP. AT ABSORPTION (EV) Y<0 FOR Y=  0',               
     >             ZTY0IZ(IZ)/ZIY0IZ(IZ))                
           IF (ZI2LIZ(IZ).GT.0.0) 
     >       CALL PRR ('  MEAN TEMP. AT ABSORPTION (EV) Y<0 FOR Y=-2L',               
     >             ZT2LIZ(IZ)/ZI2LIZ(IZ))                
           IF (CIY0IZ(IZ).GT.0.0) 
     >       CALL PRR ('  MEAN TEMP. AT ABSORPTION (EV) Y>0 FOR Y=  0',               
     >             CTY0IZ(IZ)/CIY0IZ(IZ))                
           IF (CI2LIZ(IZ).GT.0.0) 
     >       CALL PRR ('  MEAN TEMP. AT ABSORPTION (EV) Y>0 FOR Y= 2L',               
     >             CT2LIZ(IZ)/CI2LIZ(IZ))                
C                                                                               
           CALL PRR2 ('  MEAN VELOCITY AT ABSORPTION (M/S)     ',               
     >      ZRVABS(IZ)/(ZICABS(IZ)*QTIM), CRVABS(IZ)/(CICABS(IZ)*QTIM))         
           CALL PRR2 ('  MEAN ABS(VEL) AT ABSORPTION (M/S)     ',               
     >      ZRAVAV(IZ)/(ZICABS(IZ)*QTIM), CRAVAV(IZ)/(CICABS(IZ)*QTIM))         
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
           IF (VFLUID.LE.0.0 .OR. ZRAVAV(IZ).LE.0.0) THEN                       
             ZMACH = 0.0                                                        
             ZRZ   = 0.0                                                        
             ZEXIT = 0.0                                                        
           ELSE                                                                 
             ZMACH = ZRAVAV(IZ) / (ZICABS(IZ)*QTIM) / VFLUID                    
             ZRZ   = ZICABS(IZ) / ZMACH                                         
             ZEXIT = ZRAVAV(IZ) / (ZICABS(IZ)*QTIM)                             
           ENDIF                                                                
           ZENEGY = 3.0 * REAL(IZ) * ZTBS(IZ) / ZICABS(IZ) +                    
     >            5.22E-9 * CRMI * ZEXIT * ZEXIT +                              
     >            2.0 * ZRTABS(IZ) / ZICABS(IZ)                                 
C                                                                               
           CALL PRR2 ('  IMPURITY MACH NUMBER AT ABSORPTION    ',               
     >                             ZMACH, RMACH)                                
           CALL PRR2 ('  IMPACT ENERGY AT ABSORPTION (EV)      ',               
     >                            ZENEGY, RENEGY)                               
         ENDIF                                                                  
C                                                                               
C-----------------------------------------------------------------------        
C        PRINT IONISATION & RECOMBINATION DETAILS                               
C-----------------------------------------------------------------------        
C                                                                               
         CALL PRB                                                               
         CALL PRC ('NUMBER OF                     FIRST TIME  LAST TIME         
     > MEAN TIME')                                                              
         IF (CICIZS(IZ).GT.0.0) THEN                                            
           WRITE (7,9001) '  IONIZATIONS   ',NINT(CICIZS(IZ)),                  
     >       QTIM*CIFIZS(IZ),QTIM*CILIZS(IZ),QTIM*CISIZS(IZ)/CICIZS(IZ)         
         ELSE                                                                   
           WRITE (7,9001) '  IONIZATIONS   ',NINT(CICIZS(IZ))                   
         ENDIF                                                                  
         IF (CICRCS(IZ).GT.0.0) THEN                                            
           WRITE (7,9001) '  RECOMBINATIONS',NINT(CICRCS(IZ)),                  
     >       QTIM*CIFRCS(IZ),QTIM*CILRCS(IZ),QTIM*CISRCS(IZ)/CICRCS(IZ)         
         ELSE                                                                   
           WRITE (7,9001) '  RECOMBINATIONS',NINT(CICRCS(IZ))                   
         ENDIF                                                                  
 9001 FORMAT(1X,A16,3X,I7,:,1P,3X,3(2X,G9.2))                                   
C                                                                               
C-----------------------------------------------------------------------        
C        UPDATE TOTALS                                                          
C-----------------------------------------------------------------------        
C                                                                               
         RCR  = RCR  + CICUTS(IZ)                                               
         RTR  = RTR  + CRTRCS(IZ)                                               
C                                                                               
         RCAB = RCAB + CICABS(IZ) + ZICABS(IZ)                                  
         RSAB = RSAB + CISABS(IZ) + ZISABS(IZ)                                  
         RTAB = RTAB + CRTABS(IZ) + ZRTABS(IZ)                                  
         RTAV = RTAV + CRVABS(IZ) + ZRVABS(IZ)                                  
         RAVA = RAVA + CRAVAV(IZ) + ZRAVAV(IZ)                                  
         RTBS = RTBS + CTBS(IZ)   + ZTBS(IZ)                                    
         RAVMCH = RAVMCH +  RMACH*CICABS(IZ) + ZMACH*ZICABS(IZ)                 
         RAVEGY = RAVEGY + RENEGY*CICABS(IZ) + ZENEGY*ZICABS(IZ)                
         RZ0  = RZ0 + DRZ + ZRZ                                                 
         RZ1  = RZ1 + (DRZ+ZRZ)*FLOAT(IZ)                                       
         RZ2  = RZ2 + (DRZ+ZRZ)*FLOAT(IZ)*FLOAT(IZ)                             
  100  CONTINUE                                                                 
C                                                                               
C-----------------------------------------------------------------------        
C     PRINT SUMMARY DETAILS OVER ALL IONISATION STATES                          
C-----------------------------------------------------------------------        
C                                                                               
      CALL PRB                                                                  
      CALL PRC ('***   ALL   IONIZATION   STATES     ***')                      
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
     >     (IZ, (CICABS(IZ)+ZICABS(IZ))/RCAB, IZ=1,NIZS)                        
C
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED ON Y < 0 SIDE AT')            
         WRITE (7,'(3X,3(2X,A7,1P,G9.2))') 'Y = 0  ',ZICAY0/RCAB,               
     >            'Y = -2L',ZICA2L/RCAB,'X = AW ',ZICAAW/RCAB                   
C
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED ON Y < 0 SIDE'//
     >             ' IN EACH STATE')            
         DO 150 IZ = 1,NIZS           
           WRITE (7,'(5X,A7,I3,2(2X,A7,G12.5))') '  IZ = ',IZ,
     >     'Y = 0  ',ZIY0IZ(IZ)/RCAB, 'Y = -2L',ZI2LIZ(IZ)/RCAB
 150     CONTINUE 
C
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED ON Y > 0 SIDE AT')            
         WRITE (7,'(3X,3(2X,A7,1P,G9.2))') 'Y = 0  ',CICAY0/RCAB,               
     >            'Y = 2L ',CICA2L/RCAB,'X = AW ',CICAAW/RCAB                   
C
         CALL PRC ('  RATIO OF FRACTIONS ABSORBED ON Y > 0 SIDE'//
     >             ' IN EACH STATE')            
         DO 160 IZ = 1,NIZS           
           WRITE (7,'(5X,A7,I3,2(2X,A7,G12.5))') '  IZ = ',IZ,
     >     'Y = 0  ',CIY0IZ(IZ)/RCAB, 'Y =  2L',CI2LIZ(IZ)/RCAB
 160     CONTINUE 
C                                                                               
         TEXT = 0.0                                                             
         DO 200 IM = 1, 10                                                      
           TEXT = TEXT + CTEXS(IM)                                              
  200    CONTINUE                                                               
         CALL PRC ('  RATIO OF ABSORPTIONS BY EXIT TEMPS AS MULTIPLES OF        
     > TB(0)')                                                                  
         WRITE (7,'(3X,9F6.2,'' HIGHER'')') (REAL(IM)*CTTT/CTBIN,IM=1,9)        
         WRITE (7,'(4X,10F6.3)') (CTEXS(IM)/TEXT,IM=1,10)                       
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     PRINT "OTHER INFORMATION" SECTION                                         
C-----------------------------------------------------------------------        
C                                                                               
      CALL PRB                                                                  
      CALL PRI ('MEAN NUMBER OF COLLISIONS PER ION        ',                    
     >                         NINT(FACT*CICCOL))                               
      CALL PRR ('FRACTION OF IONS PASSING Y = 2L OR -2L   ',                    
     >                             (FACT*CICY2L))                               
      IF (CICY2L .GT. 0.0)                                                      
     >   CALL PRR ('  TIME FIRST ION PASSED Y = 2L OR -2L (S)',                 
     >                             (QTIM*CIFY2L))                               
      CALL PRR ('FRACTION OF IONS REACHING X = A          ',                    
     >                             (FACT*CICRXA))                               
      IF (CICRXA .GT. 0.0) THEN                                                 
         CALL PRR ('  TIME FIRST ION REACHED X = A  (S)      ',                 
     >                             (QTIM*CIFRXA))                               
         CALL PRR ('  MEAN TIME TO REACH X = A  (S)          ',                 
     >                             (QTIM*CISRXA/CICRXA))                        
         CALL PRR ('  MEAN ION TEMPERATURE AT X = A  (EV)    ',                 
     >                             (CITRXA/CICRXA))                             
      ENDIF                                                                     
      CALL PRR ('MEAN TIME SPENT OUTBOARD BY EACH ION (S) ',                    
     >                               CICOUT)                                    
      IF (CICOUT .GT. 0.0)                                                      
     >   CALL PRR ('DEEPEST OUTBOARD PENETRATION  (M)        ',                 
     >                             ABS(CRXMIN))                                 
      CALL PRR ('MEAN NIE NEAR LIMITER (NT BASED) (M**-3) ',STOTS(10))          
      CALL PRR ('MEAN ZB.NBT NEAR LIM  (NT BASED) (M**-3) ',STOTS(11))          
      CALL PRR ('MEAN ZEFF NEAR LIMITER(NT BASED)         ',STOTS(12))          
      CALL PRR ('TOTAL POWER RADIATED       PRAD  (W)     ',STOTS(4))           
      CALL PRR ('TOTAL POWER RADIATED IN SOL      (W)     ',STOTS(7))           
      CALL PRR ('TOTAL LINE RADIATION             (W)     ',STOTS(5))           
      CALL PRR ('TOTAL LINE RADIATION IN SOL      (W)     ',STOTS(8))           
      CALL PRR ('TOTAL IMPURITY CONTENT     NITOT         ',STOTS(3))           
      CALL PRR ('TOTAL PLASMA CONTENT       NBTOT         ',STOTS(2))           
      CALL PRR ('TOTAL PLASMA / VOLUME      NBAVG         ',STOTS(13))          
      CALL PRR ('RADIATION CONSTANT         LZ    (W.M**3)',STOTS(14))          
C                                                                               
      CALL PRC ('AVERAGE DELTA Y STEPS OUTBOARD (AS USED) (M)')                 
      WRITE (7,'((2X,6(I2,1P,E8.1)))') (IZ,SDYZS(IZ),IZ=1,NIZS)                 
      RETURN                                                                    
      END                                                                       
