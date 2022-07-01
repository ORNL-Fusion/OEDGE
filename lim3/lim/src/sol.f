      SUBROUTINE SOL (QYS,CEYS,CVHYS,NQYS,CTBIN,CTIBIN,CRMB,CL,
     >           CIZB,CEYOUT,CVHOUT,CYSTAG,CRMI,CSOLEF,CIOPTF)                  
      use mod_params
      IMPLICIT  none
c      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
      INTEGER   NQYS,CIZB,CIOPTF                                                
      REAL      QYS(MAXQYS),CEYS(MAXQYS),CVHYS(MAXQYS)                          
      REAL      CTBIN,CRMB,CL,CEYOUT,CVHOUT,CYSTAG,CRMI,CSOLEF                  
      REAL      CTIBIN
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  SOL:   THIS ROUTINE CALCULATES E(Y) AND VH(Y) IN THE SCRAPE OFF  *        
C  *         LAYER.                                                    *        
C  *                                                                   *        
C  *  ARGUMENT LIST :-                                                 *        
C  *  QYS    : ARRAY OF Y VALUES E AND VH TO BE CALCULATED FOR         *        
C  *  CEYS   : ARRAY FOR ELECTRIC FIELDS (V/M)                         *        
C  *  CVHYS  : ARRAY FOR BACKGROUND DRIFT VELOCITY (M/S)               *        
C  *  NQYS   : NUMBER OF Y POSITIONS EY AND VHY WANTED FOR             *        
C  *  CTBIN  : TEMPERATURE OF BACKGROUND ELECTRONS AT X=0 (EV)         *
C  *  CTIBIN : TEMPERATURE OF BACKGROUND IONS AT X=0 POSITION (EV)     *        
C  *  CRMB   : MASS OF BACKGROUND IONS (PROTON MASSES)                 *        
C  *  CL     : HALF THE SEPARATION BETWEEN LIMITERS (M)                *        
C  *  CIZB   : PLASMA IONS CHARGE                                      *        
C  *                                                                   *        
C  *  CHRIS FARRELL  (HUNTERSKIL)  MARCH 1988                          *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IQY,IYSTAG,IPOS                                                   
      REAL    RCS,EFAC,VFAC,REPL,TEMP,M,MSQR,QY,VHY                             
C                                                                               
C  *********************************************************************        
C  *  SOL0:    E(Y) AND VH(Y) ARE SET TO 0 FOR ALL Y                   *        
C  *********************************************************************        
C                                                                               
      IF (CIOPTF.EQ.0) THEN                                                     
        DO 10 IQY = 1, NQYS                                                     
          CVHYS(IQY) = 0.0                                                      
          CEYS(IQY)  = 0.0                                                      
   10   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL1: E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING :-        *        
C  *       FOR Y < L                                                   *        
C  *   E(Y)  = -(CTBIN/L).M.(1+M.M)/(1-M.M)                            *        
C  *   VH(Y) = -CS.M                                                   *        
C  *   CS    = 9.79E+03.SQRT(((CTBIN+CTIBIN)/2)*(1+CIZB)/CRMB)         *        
C  *   M     = TEMP-SQRT(TEMP.TEMP-1)                                  *        
C  *   TEMP  = 1/(1-Y/L)                                               *        
C  *        AND FOR Y > L                                              *        
C  *   E(Y)  = -E(2L-Y)                                                *        
C  *   VH(Y) = -VH(2L-Y)                                               *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.1) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)
     >          *(1.0+REAL(CIZB))/CRMB)                      
        REPL = 1.0 / CL                                                         
        EFAC = CTBIN * REPL                                                     
        DO 100  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                     
          ELSE                                                                  
            TEMP = 1.0 / ABS(QY)                                                
            M    = SQRT(TEMP * TEMP - 1.0) - TEMP                               
            MSQR = M * M                                                        
            CVHYS(IQY) = SIGN((RCS * M), QY)                                    
            CEYS(IQY) = SIGN(((EFAC*M*(1.0+MSQR)) / (1.0-MSQR)), QY)            
          ENDIF                                                                 
  100   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL2: E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING FORMULAE *         
C  *       FOR Y < L                                                   *        
C  *  E(Y)  = -CTBIN/2L                                                *        
C  *  VH(Y) = CS.(1-Y/L)                                               *        
C  *  CS    = 9.79E+03.SQRT(((CTBIN+CTIBIN)/2)*(1+CIZB)/CRMB)          *        
C  *       AND FOR Y > L                                               *        
C  *  E(Y)  = -E(2L-Y)                                                 *        
C  *  VH(Y) = -VH(2L-Y)                                                *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.2) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)*
     >         ( 1.0+REAL(CIZB))/CRMB)                      

        WRITE(6,*) 'sol:2:',RCS,CTBIN,CTIBIN,CIZB,CRMB
        EFAC = -CTBIN / (2.0 * CL)                                              
        REPL = 1.0 / CL                                                         
        DO 200  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                     
          ELSE                                                                  
            CVHYS(IQY) = RCS * QY                                               
            CEYS(IQY) = SIGN(EFAC, QY)                                          
          ENDIF                                                                 
  200   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL3:  E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING FORMULAE *        
C  *        FOR Y < L                                                  *        
C  *   E(Y)  = -(CTBIN/(PI.L)).(1-Y/L)                                 *        
C  *   VH(Y) = (4/3).CS.(1-Y/L)                                        *        
C  *   CS    = 9.79E+03.SQRT(((CTBIN+CTIBIN)/2)*(1+CIZB)/CRMB)         *        
C  *        AND FOR Y > L                                              *        
C  *   E(Y)  = -E(2L-Y)                                                *        
C  *   VH(Y) = -VH(2L-Y)                                               *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.3) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)*
     >         (1.0+REAL(CIZB))/CRMB)                      
        EFAC = CTBIN / (3.14159265 * CL)                                        
        VFAC = 1.33 * RCS                                                       
        REPL = 1.0 / CL                                                         
        DO 300  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                     
          ELSE                                                                  
            CVHYS(IQY) = VFAC * QY                                              
            CEYS(IQY) = EFAC * QY                                               
          ENDIF                                                                 
  300   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *   SOL4: E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING FORMULAE *        
C  *        FOR Y < L                                                  *        
C  *   E(Y)  = -(CTBIN/(2.L)).(1-Y/L)                                  *        
C  *   VH(Y) = CS.(1-Y/L)                                              *        
C  *   CS    = 9.79E+03.SQRT(((CTBIN+CTIBIN)/2)*(1+CIZB)/CRMB)         *        
C  *        AND FOR Y > L                                              *        
C  *   E(Y)  = -E(2L-Y)                                                *        
C  *   VH(Y) = -VH(2L-Y)                                               *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.4) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)*
     >         (1.0+REAL(CIZB))/CRMB)                      
        EFAC = CTBIN / (2.0 * CL)                                               
        VFAC = RCS                                                              
        REPL = 1.0 / CL                                                         
        DO 400  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                     
          ELSE                                                                  
            CVHYS(IQY) = VFAC * QY                                              
            CEYS(IQY) = EFAC * QY                                               
          ENDIF                                                                 
  400   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL5:E(Y) AND VH(Y) ARE CONSTANT, BEARING IN MIND ANY CHANGES OF *        
C  *      SIGN IN THE FOUR REGIONS ...                                 *        
C  *                                                                   *        
C  *   |  +VH(Y)   |   -VH(Y)    |   +VH(Y)    |   -VH(Y)   |          *        
C  *   |   +E(Y)   |    -E(Y)    |    +E(Y)    |    -E(Y)   |          *        
C  *   |-----------|-------------|-------------|------------|----> Y   *        
C  * -2L          -L             0             L            2L         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.5) THEN                                                 
        DO 500 IQY = 1, NQYS/2                                                  
          CEYS(IQY)  = CEYOUT                                                   
          CVHYS(IQY) = CVHOUT                                                   
  500   CONTINUE                                                                
        DO 510 IQY = (NQYS/2)+1, NQYS                                           
         CEYS(IQY)  =-CEYOUT                                                   
          CVHYS(IQY) =-CVHOUT                                                   
  510   CONTINUE                                                                
c slmod end
C                                                                               
C  *********************************************************************        
C  *  SOL6:  STAGNATION PT SPECIFIED NEAR THE MIDDLE OF THE IONISATION *        
C  *  DISTRIBUTIONS  (SEE NOTE 132).  SPECIFIED YSTAG GIVES 3 ZONES IN *        
C  *  +Y REGION, WHICH ARE MIMICKED FOR -Y REGION  (ALTHOUGH SOL6 WILL *        
C  *  BE USED FOR GAS PUFF SIMULATIONS WHEN -Y REGION DOESN'T REALLY   *        
C  *  MATTER).                                                         *        
C  *                                                                   *        
C  *  |                                REGION1   REGION2   REGION3     *        
C  *  |                                                                *        
C  *  |--------|----------|----------|---------|---------|---------|->Y*        
C  * -2L    -2YSTAG     -YSTAG       0        YSTAG    2YSTAG     2L   *        
C  *                                                                   *        
C  *  IN REGION1:  E(Y) = -TB/2YSTAG   VH(Y) = -9788.SQRT(TB/MI)       *        
C  *     REGION2:  E(Y) = +TB/2YSTAG   VH(Y) = +9788.SQRT(TB/MI)       *        
C  *     REGION3:  E(Y) = 0            VH(Y) = 0                       *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.6) THEN                                                 
        IYSTAG = IPOS (CYSTAG, QYS, NQYS-1)                                     
        VHY    = 9788.0 * SQRT (((CTBIN+CTIBIN)/2) / CRMI)                     
        DO 600 IQY = 1, IYSTAG                                                  
          CEYS(IQY)  = -CTBIN / (2.0 * CYSTAG)                                  
          CVHYS(IQY) = -VHY                                                     
  600   CONTINUE                                                                
        DO 610 IQY = IYSTAG+1, 2*IYSTAG                                         
          CEYS(IQY)  = CTBIN / (2.0 * CYSTAG)                                   
          CVHYS(IQY) = VHY                                                      
  610   CONTINUE                                                                
        DO 620 IQY = 2*IYSTAG+1, NQYS                                           
          CEYS(IQY)  = 0.0                                                      
          CVHYS(IQY) = 0.0                                                      
  620   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL7:  A VARIATION ON SOL6 (NOTE 132), SOL7 DEFINED BY NOTE 142. *        
C  *                                                                   *        
C  *  |                                REGION1   REGION2   REGION3     *        
C  *  |                                                                *        
C  *  |--------|----------|----------|---------|---------|---------|->Y*        
C  * -2L    -2YSTAG     -YSTAG       0        YSTAG   10YSTAG     2L   *        
C  *                                                                   *        
C  *  IN REGION1:  E(Y) = -TB/2YSTAG   VH(Y) = -9788.SQRT(TB/MI)       *        
C  *     REGION2:  E(Y) = +TB/20YSTAG  VH(Y) = +9788.SQRT(TB/MI).      *        
C  *     REGION3:  E(Y) = 0            VH(Y) = 0     (Y-YSTAG)/10YSTAG *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.7) THEN                                                 
        IYSTAG = IPOS (CYSTAG, QYS, NQYS-1)                                     
        VHY    = 9788.0 * SQRT (((CTBIN+CTIBIN)/2) / CRMI)                   
        DO 700 IQY = 1, IYSTAG                                                  
          CEYS(IQY)  = -CTBIN / (2.0 * CYSTAG)                                  
          CVHYS(IQY) = -VHY                                                     
  700   CONTINUE                                                                
        DO 710 IQY = IYSTAG+1, 10*IYSTAG                                        
          CEYS(IQY)  = CTBIN / (20.0 * CYSTAG)                                  
          CVHYS(IQY) = VHY * (QYS(IQY) - CYSTAG) / (10.0 * CYSTAG)              
  710   CONTINUE                                                                
        DO 720 IQY = 10*IYSTAG+1, NQYS                                          
          CEYS(IQY)  = 0.0                                                      
          CVHYS(IQY) = 0.0                                                      
  720   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL8:       (BASED ON SOL4 AND NOTE 168)                         *        
C  *      SIMILAR TO SOL4 WITH A SOL ENHANCEMENT FACTOR THROWN IN.     *        
C  *      E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING FORMULAE    *        
C  *        FOR Y < L                                                  *        
C  *   E(Y)  = -(CTBIN/(2.L)).(1-Y/L).SOLEF                            *        
C  *   VH(Y) = CS.(1-Y/L).SOLEF                                        *        
C  *   CS    = 9.79E+03.SQRT(CTBIN*(1+CIZB)/CRMB)                      *        
C  *        AND FOR Y > L                                              *        
C  *   E(Y)  = -E(2L-Y)                                                *        
C  *   VH(Y) = -VH(2L-Y)                                               *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.8) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)*
     >         (1.0+REAL(CIZB))/CRMB)                      
        EFAC = CTBIN / (2.0 * CL) * CSOLEF                                      
        VFAC = RCS * CSOLEF                                                     
        REPL = 1.0 / CL                                                         
        DO 800  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                     
          ELSE                                                                  
            CVHYS(IQY) = VFAC * QY                                              
            CEYS(IQY) = EFAC * QY                                               
          ENDIF                                                                 
  800   CONTINUE                                                                
C                                                                               
C  *********************************************************************        
C  *  SOL9 : E(Y) AND VH(Y) ARE CALCULATED FROM THE FOLLOWING FORMULAE *       
C  *       FOR Y < L                                                   *        
C  *  E(Y)  = -CTBIN/2L                                                *        
C  *  VH(Y) = CS                                                       *        
C  *  CS    = 9.79E+03.SQRT(((CTBIN+CTIBIN)/2)*(1+CIZB)/CRMB)          *        
C  *       AND FOR Y > L                                               *        
C  *  E(Y)  = -E(2L-Y)                                                 *        
C  *  VH(Y) = -VH(2L-Y)                                                *        
C  *********************************************************************        
C                                                                               
      ELSEIF (CIOPTF.EQ.9) THEN                                                 
        RCS = 9.79E+03 * SQRT(((CTBIN+CTIBIN)/2)*
     >         ( 1.0+REAL(CIZB))/CRMB)                      
        EFAC = -CTBIN / (2.0 * CL)                                              
        REPL = 1.0 / CL                                                         
        DO 900  IQY = 1, NQYS                                                   
          QY = QYS(IQY) * REPL - 1.0                                            
          IF (QY .EQ. 0.0)  THEN                                                
            CVHYS(IQY) = 0.0                                                    
            CEYS(IQY) = 0.0                                                   
          ELSE                                                                  
            CVHYS(IQY) = SIGN(RCS, QY)                                          
            CEYS(IQY) = SIGN(EFAC, QY)                                          
          ENDIF                                                                 
  900   CONTINUE                                                                
      ENDIF                                                                     
      RETURN                                                                    
      END                                                                       
