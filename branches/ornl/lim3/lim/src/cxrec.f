C@PROCESS OPT(1),VECTOR(LEV(0))                                               
C                                                                               
C                                                                               
      SUBROUTINE CXREC (NIZS,CION,CIOPTI,CIZB,CL,CRMB,CVCX,                     
     >                  CNHC,CNHO,CLAMHX,CLAMHY)                                
      use mod_params
      use mod_comt2
      use mod_comxyt
      IMPLICIT    none                                                
      INTEGER     NIZS,CION,CIOPTI,CIZB                                         
      REAL        CL,CRMB,CVCX,CNHC,CNHO,CLAMHX,CLAMHY                          
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  CXREC: THIS ROUTINE DERIVED FROM NOTE 89                         *        
C  *         VCX MODIFICATION GIVEN IN NOTE 173 APRIL 1988             *        
C  *  DEALS WITH CHARGE EXCHANGE RECOMBINATION, IMPORTANT FOR          *        
C  *  HYDROGENIC PLASMAS.                                              *        
C  *                                                                   *        
C  *  VECTORISATION SWITCHED OFF DUE TO POSSIBILITY OF DIVIDE BY ZERO  *        
C  *  ERRORS WITHIN FINAL DO-LOOPS.  SHOULDN'T HAPPEN IN SCALAR SINCE  *        
C  *  IF STATEMENTS TAKE CARE OF IT.  HOWEVER, THE IBM OPTIMISER       *        
C  *  SHUFFLES THE LOOP AND THE ERRORS DO OCCUR.  HENCE THE SPECIAL    *        
C  *  COMPILER OPTIONS ABOVE FOR THIS ROUTINE ONLY.                    *        
C  *                                                                   *        
C  *                                      C.M.FARRELL   JANUARY 1988   *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
c      INCLUDE     'params'                                                      
C     INCLUDE     (PARAMS)                                                      
c      INCLUDE     'comt2'                                                       
C     INCLUDE     (COMT2)                                                       
c      INCLUDE     'comxyt'                                                      
C     INCLUDE     (COMXYT)                                                      
C                                                                               
      INTEGER IX,IY,IZ                                                          
      REAL    X,Y,VCX,V,Q(28),RIZB,SIGCX                                        
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE NEUTRAL HYDROGEN ATOM DENSITY  (FUNCTION OF X AND Y)            
C     VARIOUS OPTIONS ALLOWED, KEYED WITH CIOPTI                                
C-----------------------------------------------------------------------        
C                                                                               
C  OPTION 0                                                                     
C  --------                                                                     
      IF (CIOPTI.EQ.0) THEN                                                     
        CALL RZERO (CNHS, MAXNXS*(2*MAXNYS+1))                                  
C                                                                               
C  OPTIONS 1,2                                                                  
C  -----------                                                                  
      ELSEIF (CIOPTI.EQ.1.OR.CIOPTI.EQ.2) THEN                                  
        DO 130 IX = 1, NXS                                                      
          X = XOUTS(IX)                                                         
          DO 120 IY = 1, NYS                                                    
            Y = YOUTS(IY)                                                       
            IF     (Y.LE.CL .AND. X.GE.0.0) THEN                                
              CNHS(IX,IY) = CNHC + CNHO * EXP(-X/CLAMHX) *EXP(-Y/CLAMHY)        
            ELSEIF (Y.LE.CL .AND. X.LT.0.0) THEN                                
              CNHS(IX,IY) = CNHO * EXP (-Y/CLAMHY)                              
            ELSEIF (Y.GT.CL) THEN                                               
              CNHS(IX,IY) = CNHS(IX,NYS+1-IY)                                   
            ENDIF                                                               
            CNHS(IX,-IY) = CNHS(IX,IY)                                          
  120     CONTINUE                                                              
          CNHS(IX,0) = CNHS(IX,1)                                               
  130   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     SET UP OUTER LOOPS BEFORE CALCULATING Q VALUES                            
C-----------------------------------------------------------------------        
C                                                                               
      DO 240 IY = -NYS, NYS                                                     
C                                                                               
      DO 230 IX = 1, NXS                                                        
        IF (CIOPTI.EQ.0.OR.CIOPTI.EQ.1) THEN                                    
          VCX = 1.56E4 * SQRT (CTEMBSI(IX,IY)/CRMB)                            
        ELSEIF (CIOPTI.EQ.2) THEN                                               
          VCX = CVCX                                                            
        ENDIF                                                                   
        V   = 1.E-4 * VCX                                                       
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE Q VALUES (FORMULAE DIFFERENT FOR EACH SPECIES)                  
C-----------------------------------------------------------------------        
C                                                                               
        DO 200 IZ = 1, NIZS                                                     
          Q(IZ) = 0.0                                                           
  200   CONTINUE                                                                
C                                                                               
C       HELIUM                                                                  
C       ------                                                                  
        IF     (CION.EQ.2) THEN                                                 
          Q(1) =-0.884 + 0.021*V                                                
          Q(2) = 0.114 - 0.038*V + 0.001929*V*V                                 
C                                                                               
C       BERYLLIUM                                                               
C       ---------                                                               
        ELSEIF (CION.EQ.4) THEN                                                 
          Q(1) = 0.061 * V * (1.0 + V / 61.5)                                   
          Q(2) = 9.16 * EXP (-18.68 / V)                                        
          Q(3) = 1.2
          Q(4) = 2.109 + 1.048*V - 0.00718*V*V                                  
C                                                                               
C       CARBON                                                                  
C       ------                                                                  
        ELSEIF (CION.EQ.6) THEN                                                 
          Q(1) = 0.061 * V * (1.0 + V / 61.5)                                   
          Q(2) = 9.16 * EXP (-18.68 / V)                                        
          Q(3) = 21.93 - 1.93*V + 0.07*V*V - 0.000753*V*V*V                     
          Q(4) = 2.66 -1.76*LOG(ABS(V))/V + 17.14/V + 37.87*EXP(-8.51/V)        
          Q(5) = 3.24-377.3*LOG(ABS(V))/V/V+276.6/V +0.32*V*EXP(-6.52/V)        
          Q(6) = V * EXP (0.676 - 15.05/V)                                      
C                                                                               
C       OXYGEN                                                                  
C       ------                                                                  
        ELSEIF (CION.EQ.8) THEN                                                 
          Q(1) = 15.83  - 0.374 *V + 0.003913*V*V - 0.0000138*V*V*V             
          Q(2) =  8.84  - 0.29  *V + 0.00283 *V*V                               
          Q(3) = 27.22  + 1.013 *V - 0.0203  *V*V                               
          Q(4) = 41.53  - 0.994 *V + 0.0147  *V*V                               
          Q(5) = 58.75  - 0.685 *V + 0.00927 *V*V                               
          Q(6) = 25.772 + 0.755 *V - 0.0143  *V*V + 0.0000794*V*V*V             
          Q(7) = 31.25  - 0.0553*V                                              
          Q(8) =  1.493 + 1.809 *V - 0.0116  *V*V                               
C                                                                               
C       IRON   (NOTE 114:  SIGCX = 1.E-19 * Z, HENCE SET Q AS BELOW)            
C       ----                                                                    
        ELSEIF (CION.EQ.26) THEN                                                
          DO 205 IZ = 1, NIZS                                                   
            Q(IZ) = 10.0 * REAL (IZ)                                            
  205     CONTINUE                                                              
C                                                                               
C       OTHER IMPURITIES                                                        
C       ----------------                                                        
        ELSEIF (IX.EQ.1.AND.CIOPTI.NE.0) THEN                                   
          CALL PRC ('WARNING: NO CHARGE EXCHANGE RECOMBINATION DATA AVAI        
     >LABLE FOR THIS IMPURITY')                                                 
        ENDIF                                                                   
C                                                                               
C       WRITE (6,'('' CXREC: IY,IX,V,NB,Q'',2I5,(1X,6G11.4))')                  
C    >    IY,IX,V,CRNBS(IX,IY),(Q(IZ),IZ=1,NIZS)                                
C                                                                               
C-----------------------------------------------------------------------        
C     CALCULATE NEW RECOMBINATION TIMES, STORE IN CFCXS ARRAY                   
C-----------------------------------------------------------------------        
C                                                                               
        RIZB = REAL (CIZB)                                                      
        DO 220 IZ = 1, NIZS                                                     
            SIGCX = CNHS(IX,IY) * Q(IZ) * 1.E-20 * VCX / CRNBS(IX,IY)           
            IF (CFRCS(IX,IY,IZ).GT.0.0) SIGCX = SIGCX +                         
     >        1.0 / (RIZB * CRNBS(IX,IY) * CFRCS(IX,IY,IZ))                     
            IF (SIGCX.GT.0.0) THEN                                              
              CFCXS(IX,IY,IZ) = 1.0 / (RIZB * CRNBS(IX,IY) * SIGCX)             
            ELSE                                                                
              CFCXS(IX,IY,IZ) = 0.0                                             
            ENDIF                                                               
  220   CONTINUE                                                                
  230 CONTINUE                                                                  
C                                                                               
  240 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
