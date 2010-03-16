      SUBROUTINE RINTX (VS,V3,IPLANE,MAXIZ,XINTS,IFOLD,RV,XV,YV,IVU,            
     >                  SSS,NYSLIM,FP,FT,AVER,IXMIN,IXMAX)                
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    V3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)                
      REAL    XINTS(-MAXNYS:MAXNYS,-2:MAXIZS+1),RV,XV,YV,FP,FT,AVER             
      INTEGER IPLANE,MAXIZ,IFOLD,IVU,NYSLIM,IXMIN,IXMAX                       
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTX:  TAKE REAL ARRAY VS  (OR V3(,,,IPLANE) IF IPLANE<>99) AND *        
C  *  INTEGRATE OVER X DIMENSION TO OBTAIN RESULTS ARRAY XINTS.  THE   *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO XINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTXS.          *        
C  *                                                                   *        
C  *  VIEWPOINT CALCULATIONS ADDED JULY 88 FOR BELT LIMITER CASES.     *        
C  *                                                                   *        
C  *  ARGUMENT "AVER" ADDED 14/11/88 FOR AVERAGING TEMPS, DELTAY STEPS *        
C  *                                                                   *        
C  *  ARGUMENTS "IXMIN,IXMAX" ADDED 29/10/90 FOR INTEGRATION RANGE SPEC*
C  *                                                                   *
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH   1988                          *        
C  *  DAVID ELDER                OCT 29  1990                          *
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTXS(-MAXNYS:MAXNYS,-2:MAXIZS+1),DXSUM1,DXSUM2          
      INTEGER IZ,IY,IX,JY,KY,LY                                                 
      REAL    FRAC,TMP                                                          
C                                                                               
      CALL DZERO (DNTXS, (2*MAXNYS+1)*(MAXIZS+4))                               
C                                                                               
      IF (IVU.GT.0) THEN                                                        
        CALL VIEW2 (RV,XV,YV,SSS,NYSLIM,IVU)                                    
        DO 230 IY = -NYSLIM, NYSLIM                                             
          IF (IY.EQ.0) GOTO 230                                                 
          JY = IABS (IY)                                                        
          DO 220 IX = IXMIN,IXMAX                                          
            IF (XPPPP(IX,IY).LT.0.0.AND.IY.LT.0) GOTO 220                       
            DO 210 KY = -NYS, NYS                                               
              IF (KY.EQ.0) GOTO 210                                             
              LY = IABS (KY)                                                    
              IF (KY.LT.0) THEN                                                 
                TMP = FRAC (YWIDS(LY)-YS(LY),YWIDS(LY),                         
     >                                       YPPPP(IX,IY),YPWID(IX,IY))         
              ELSE                                                              
                TMP = FRAC (YS(LY),YWIDS(LY),YPPPP(IX,IY),YPWID(IX,IY))         
              ENDIF                                                             
              IF (TMP.LE.0.0) GOTO 210                                          
              TMP = TMP * XWIDS(IX) * YWIDS(JY) / YWIDS(LY) / SSQ(IX,IY)        
              DO 200 IZ = -1, MAXIZ                                             
               IF (IPLANE.EQ.99) THEN                                           
                DNTXS(KY,IZ) = DNTXS(KY,IZ) + DBLE (TMP * VS(IX,IY,IZ))         
               ELSEIF (JY.LE.NY3D) THEN                                         
                DNTXS(KY,IZ) = DNTXS(KY,IZ) +                                   
     >                                 DBLE (TMP * V3(IX,IY,IZ,IPLANE))         
               ENDIF                                                            
  200         CONTINUE                                                          
  210       CONTINUE                                                            
  220     CONTINUE                                                              
  230   CONTINUE                                                                
C                                                                               
      ELSEIF (IVU.EQ.0) THEN                                                    
        DO 130 IZ = -1, MAXIZ                                                   
         DO 120 IY = 1, NYS                                                     
          DXSUM1 = 0.0D0                                                        
          DXSUM2 = 0.0D0                                                        
          DO 110 IX = IXMIN,IXMAX                                          
            IF (IPLANE.EQ.99) THEN                                              
              DXSUM1 = DXSUM1 + DBLE (VS(IX,IY,IZ) * XWIDS(IX))                 
              DXSUM2 = DXSUM2 + DBLE (VS(IX,-IY,IZ) * XWIDS(IX))                
            ELSEIF (IY.LE.NY3D) THEN                                            
              DXSUM1 = DXSUM1 + DBLE (V3(IX,IY,IZ,IPLANE) * XWIDS(IX))          
              DXSUM2 = DXSUM2 + DBLE (V3(IX,-IY,IZ,IPLANE) * XWIDS(IX))         
            ENDIF                                                               
  110     CONTINUE                                                              
          IF (IFOLD.EQ.0) THEN                                                  
            DNTXS(IY,IZ)  = DXSUM1                                              
            DNTXS(-IY,IZ) = DXSUM2                                              
          ELSEIF (IFOLD.EQ.1) THEN                                              
            DNTXS(IY,IZ)  = 0.5 * (DXSUM1 + DXSUM2)                             
            DNTXS(-IY,IZ) = 0.5 * (DXSUM1 + DXSUM2)                             
          ENDIF                                                                 
  120    CONTINUE                                                               
  130   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 300 IZ = -1, MAXIZ                                                     
        DNTXS(0,IZ) = 0.5 * (DNTXS(1,IZ) + DNTXS(-1,IZ))                        
  300 CONTINUE                                                                  
      DO 330 IY = -NYS, NYS                                                     
        DNTXS(IY,-2) = FT * DNTXS(IY,0) - FP * DNTXS(IY,-1)                     
        DO 310 IZ = 1, MAXIZ                                                    
          DNTXS(IY,MAXIZ+1) = DNTXS(IY,MAXIZ+1) + DNTXS(IY,IZ)                  
  310   CONTINUE                                                                
        DO 320 IZ = -2, MAXIZ+1                                                 
          XINTS(IY,IZ) = SNGL (DNTXS(IY,IZ)) / AVER                             
  320   CONTINUE                                                                
  330 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE RINTY (VS,V3,IPLANE,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,            
     >                  SSS,NYSLIM,FP,FT,AVER,IYMIN,IYMAX)                  
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    V3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)                
      REAL    YINTS(MAXNXS,-2:MAXIZS+1),RV,XV,YV,FP,FT,AVER                     
      INTEGER IPLANE,MAXIZ,IFOLD,IVU,NYSLIM,IYMIN,IYMAX
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTY:  TAKE REAL ARRAY VS  (OR V3(,,,IPLANE) IF IPLANE<>99) AND *        
C  *  INTEGRATE OVER Y DIMENSION TO OBTAIN RESULTS ARRAY YINTS.  THE   *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO YINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTYS.          *        
C  *                                                                   *        
C  *  INTEGRATION IS EITHER STRAIGHTFORWARD, OR BY METHOD OF NOTE 91   *        
C  *  WHICH TAKES INTO ACCOUNT AN OBSERVER'S VIEWPOINT RELATIVE TO THE *        
C  *  Y AXIS.  WE ONLY USE VALUES FROM -RV/3 TO RV/3  (OR FROM -L TO L *        
C  *  IF SMALLER).  POINTS WHICH WOULD BE PROJECTED ON TO THE X AXIS   *        
C  *  BEYOND THE WALL OR INSIDE THE CENTRE ARE IGNORED.  CONTRIBUTIONS *        
C  *  FROM -Y PROJECTED ONTO X' AXIS <0 ARE IGNORED  (SINCE THE LINE-  *        
C  *  OF-SIGHT IS BLOCKED BY THE LIMITER).  TAKE INTO ACCOUNT ANY      *        
C  *  DIFFERENCES IN BIN-WIDTHS FROM X --> X' SYSTEMS.                 *        
C  *  MODIFICATION FROM NOTE 172: SCALE PLRPS BY DISTANCE SQUARED SSQ  *        
C  *                                                                   *        
C  *  QUICK FIX: NOTE 240.  CHOOSE IVU=-1 FOR REGION1, -2 FOR REGION2  *        
C  *  RESULTS WHEN RUNNING ELONGATION CASE.  SEE ALSO NOTE 211.        *        
C  *                                                                   *        
C  *  ARGUMENT "AVER" ADDED 14/11/88 FOR AVERAGING TEMPS, DELTAY STEPS *        
C  *                                                                   *        
C  *  ARGUMENTS "IYMIN,IYMAX" ADDED 29/10/90 FOR INTEGRATION RANGE SPEC*
C  *                                                                   *
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH   1988                          *        
C  *  DAVID ELDER                OCT 29  1990                          *
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTYS(MAXNXS,-2:MAXIZS+1),DYSUM                          
      INTEGER IZ,IY,IX,JY,JX                                       
      REAL    FRAC,TMP                                                          
C                                                                               
      CALL DZERO (DNTYS, MAXNXS*(MAXIZS+4))                             
C   
C
      IF (IVU.GT.0) THEN                                                        
        CALL VIEW2 (RV,XV,YV,SSS,NYSLIM,IVU)                                    
C
C       MODIFY LIMITS OF INTEGRATION TO MATCH ANY SPECIFIED 
C       VALUES SMALLER THAN YLIMIT (VIEW RANGE- SPECIFIED IN 
C       INPUT.  
C
        IYMIN = MAX(-NYSLIM,IYMIN)
        IYMAX = MIN( NYSLIM,IYMAX)         
        DO 230 IY = IYMIN, IYMAX                                            
          IF (IY.EQ.0) GOTO 230                                                 
          JY = IABS (IY)                                                        
          DO 220 IX = 1, NXS                                                    
            IF (XPPPP(IX,IY).LT.0.0.AND.IY.LT.0) GOTO 220                       
            DO 210 JX = 1, NXS                                                  
              TMP = FRAC (XS(JX),XWIDS(JX),XPPPP(IX,IY),XPWID(IX,IY))           
              IF (TMP.LE.0.0) GOTO 210                                          
              TMP = TMP * YWIDS(JY) * XWIDS(IX) / XWIDS(JX) / SSQ(IX,IY)        
              DO 200 IZ = -1, MAXIZ                                             
               IF (IPLANE.EQ.99) THEN                                           
                DNTYS(JX,IZ) = DNTYS(JX,IZ) + DBLE (TMP * VS(IX,IY,IZ))         
               ELSEIF (JY.LE.NY3D) THEN                                         
                DNTYS(JX,IZ) = DNTYS(JX,IZ) +                                   
     >                                 DBLE (TMP * V3(IX,IY,IZ,IPLANE))         
               ENDIF                                                            
  200         CONTINUE                                                          
  210       CONTINUE                                                            
  220     CONTINUE                                                              
  230   CONTINUE                                                                
C                                                                               
      ELSE                                                                      
        DO 270 IZ = -1, MAXIZ                                                   
          DO 260 IX = 1, NXS                                                    
            DYSUM = 0.0                                                         
            DO 250 IY = IYMIN, IYMAX
              IF (IY.EQ.0) GOTO 250 
              JY = IABS(IY)                                       
              IF (YS(JY).GT.0.2501*CL.AND.YS(JY).LE.0.7501*CL) THEN             
                IF (IVU.EQ.-1) GOTO 250                                         
              ELSE                                                              
                IF (IVU.EQ.-2) GOTO 250                                         
              ENDIF                                                             
              IF (IPLANE.EQ.99) THEN                                            
                DYSUM = DYSUM + DBLE (YWIDS(JY) *                          
     >                  VS(IX,IY,IZ))                           
              ELSEIF (JY.LE.NY3D) THEN                                      
                DYSUM = DYSUM + DBLE (YWIDS(JY) *                           
     >                  V3(IX,IY,IZ,IPLANE))            
              ENDIF                                                             
  250       CONTINUE                                                            
            DNTYS(IX,IZ) = DYSUM                                                
  260     CONTINUE                                                              
  270   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 320 IX = 1, NXS                                                        
        DNTYS(IX,-2) = FT * DNTYS(IX,0) - FP * DNTYS(IX,-1)                     
        DO 300 IZ = 1, MAXIZ                                                    
          DNTYS(IX,MAXIZ+1) = DNTYS(IX,MAXIZ+1) + DNTYS(IX,IZ)                  
  300   CONTINUE                                                                
        DO 310 IZ = -2, MAXIZ+1                                                 
          YINTS(IX,IZ) = SNGL (DNTYS(IX,IZ)) / AVER                             
  310   CONTINUE                                                                
  320 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       

      SUBROUTINE RINTM (V3,ISTATE,IPLANE,IFOLD,NPTS,MPTS,SURFAS,            
     >               XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >               COORD1,COORD2,POUTS,ARRLIM)
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      INTEGER ARRLIM
      REAL    V3(MAXNXS,-MAXY3D:MAXY3D,-1:ARRLIM,-MAXNPS:MAXNPS)                
      REAL    COORD1(90),COORD2(90),POUTS(-MAXNPS:MAXNPS)
      REAL    SURFAS(192,192),XMIN,XMAX,YMIN,YMAX,PMIN,PMAX                  
      INTEGER ISTATE,IPLANE,IFOLD,NPTS,NIZS,MPTS,IPLOT                        
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTM: THIS ROUTINE PREPARES AN ARRAY OF SURFACE HEIGHTS IN      *        
C  *  ARRAY SURFAS, FOR A GIVEN IONISATION STATE/PLANE COMBINATION.    *        
C  *  WRITTEN FOR USE WITH MESH/CONTOUR GRAPH OPTION. FOR MESH PLOTS   *
C  *  IT REQUIRES THAT THE IRREGULAR GRID OF X BINS AND Y BINS BE      *
C  *  CONVERTED INTO A REGULAR GRID STRUCTURE, NPTS BY MPTS IN THE TWO *
C  *  NON-INTEGRATED DIRECTIONS. THE LIMIT                             *        
C  *  ON NPTS IS 192 SINCE THIS IS A GHOST80 LIMITATION.               *        
C  *  MATLAB LIMITS NPTS,MPTS TO 90                                    *
C  *  THE SUMMING IS DONE IN D.P. FOR ACCURACY.                        *        
C  *                                                                   *
C  *  THE LIMITER EDGE SUPERPOSITION IS EXCLUDED BECAUSE OF THE WAY IT *
C  *  WOULD VARY DEPENDING ON WHICH PROJECTION WAS TO BE DISPLAYED     *
C  *                                                                   *
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  FEBRUARY 1988                         *        
C  *  D. ELDER , MARCH 15, 1990                                        *
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER  IPT,IPOS,IRS1(90),IRS2(90),IRX,IRY,IZ1,IZ2,IZ,IQX,IMOD 
      INTEGER  IXMIN,IXMAX,IYMIN,IYMAX,IPMIN,IPMAX,IRP,IN1,IN2,J      
      REAL     RATIO,RS1(90),RS2(90),VAL1,VAL2,X,Y(4)            
      DOUBLE PRECISION DURFAS(90,90)                                          

      WRITE(6,*) 'ARGS TO RINTM:',ISTATE,IPLANE,IFOLD,NPTS,MPTS,
     >    XMIN,XMAX,YMIN,YMAX,PMIN,PMAX,NIZS,IPLOT,
     >    (POUTS(J),J=-MAXNPS,MAXNPS),SURFAS(1,1)      
C
C---- FIRST INTEGRATE AND THEN MOVE TO A REGULAR GRID IF NECESSARY
C
C
c jdemod - comment out print out for now - takes up a lot of space
c slmod
c      WRITE(6,*) 'Output V3:',NPTS,MPTS
c slmod end      
c      WRITE(6,'(10G13.5)') ((V3(1,IN1,1,IN2),IN1=-MAXY3D,MAXY3D),
c     >                     IN2 = -MAXNPS,MAXNPS )
      IXMIN=IPOS(XMIN,XS,NXS-1)
      IXMAX=IPOS(XMAX,XS,NXS-1)
      IF (YMIN.GE.0) THEN 
        IYMIN=IPOS(YMIN,YS,NY3D-1)
      ELSE
        IYMIN=-IPOS(-YMIN,YS,NY3D-1)
      ENDIF 
      IF (YMAX.GE.0) THEN 
        IYMAX=IPOS(YMAX,YS,NY3D-1)
      ELSE
        IYMAX=-IPOS(-YMAX,YS,NY3D-1)
      ENDIF 
      IF (PMIN.GE.0) THEN
        IPMIN=IPOS(PMIN,PS,MAXNPS-1)
      ELSE
        IPMIN=-IPOS(-PMIN,PS,MAXNPS-1)
      ENDIF
      IF (PMAX.GE.0) THEN
        IPMAX=IPOS(PMAX,PS,MAXNPS-1)
      ELSE
        IPMAX=-IPOS(-PMAX,PS,MAXNPS-1)
      ENDIF 

      WRITE(6,*) 'IX IY IP:',IXMIN,IXMAX,IYMIN,IYMAX,IPMIN,IPMAX

      CALL DZERO (DURFAS, 90*90)                                              
      IF (ISTATE.EQ.NIZS+1) THEN                                                
        IZ1 = 1                                                                 
        IZ2 = NIZS                                                              
      ELSE                                                                      
        IZ1 = ISTATE                                                            
        IZ2 = ISTATE                                                            
      ENDIF                                                                     

      WRITE(6,*) 'IZS:',IZ1,IZ2
      WRITE(6,*) 'MAX:', MAXNXS,MAXY3D,MAXIZS,MAXNPS 
C                                                                               
C     THE ARRAYS OF DATA HAVE A ZERO INDEXED ROW IN THE Y-DIRECTION
C     THAT IS NOT USED. THIS MUST BE EDITED OUT OF THE DATA TO BE 
C     PLOTTED. HOWEVER, THE INPUT ALLOWS SPECIFYING ASSYMMETRIC RANGES
C     AS WELL AS ONES WHICH DO NOT CROSS ZERO. THIS MAKES SEPARATE LOOP 
C     IMPLEMENTATIONS CLUMSY. HERE A VARIABLE IMOD IS USED, IF THE 
C     INDEX CROSSES THROUGH ZERO THEN FURTHER ENTRIES ARE SHIFTED BY 
C     ONE POSITION AND THE TOTAL NUMBER REDUCED BY ONE.
C     UNFORTUNATELY, THIS MAY INHIBIT VECTORIZATION.
C
C     D.ELDER, APRIL 15 1990
C
      IMOD=0
      DO 220 IZ = IZ1, IZ2                                                      
        DO 210 IRY = IYMIN,IYMAX                                           
          IF (IRY.EQ.0) THEN 
             IMOD=-1
             GOTO 210
          ENDIF  
c slmod
c          WRITE(6,*) ' '
c slmod end
          DO 200 IRX = IXMIN,IXMAX                                             
            DO 190 IRP = IPMIN,IPMAX
              VAL1 = V3(IRX,IRY,IZ,IRP)                              
              VAL2 = V3(IRX,-IRY,IZ,IRP)  
              IF (IPLANE.EQ.0) THEN                                          
                IN1 = IRY-IYMIN+1+IMOD
                IN2 = IRP-IPMIN+1                
              ELSEIF (IPLANE.EQ.1) THEN                                      
                IN1 = IRX-IXMIN+1
                IN2 = IRP-IPMIN+1
              ELSEIF (IPLANE.EQ.2) THEN
                IN1 = IRX-IXMIN+1
                IN2 = IRY-IYMIN+1+IMOD
              ENDIF
              IF (IFOLD.EQ.0) THEN                                            
                DURFAS(IN1,IN2) = DURFAS(IN1,IN2) + DBLE (VAL1)                 
              ELSE                                                              
                DURFAS(IN1,IN2)=DURFAS(IN1,IN2)+DBLE(0.5*(VAL1+VAL2))      
              ENDIF                                                             
c slmod 
c             WRITE(6,*) 'DURFAS:',IN1,IN2,DURFAS(IN1,IN2),IRX,IRY,IRP
c slmod
  190       CONTINUE 
  200     CONTINUE                                                              
  210   CONTINUE                                                                
  220 CONTINUE                                                                       

      IF (IPLOT.EQ.2) THEN
C                                                                               
C---- SET INDEXING ARRAYS FOR GOING FROM IRREGULAR TO REGULAR GRIDS ...         
C                                                                               
C
C     FIRST COORDINATE  0 = Y , 1 = X, 2 = X
C
      DO 100 IPT = 1, NPTS                                                      
        RATIO = REAL (IPT) / REAL (NPTS+1)                                      
        IF (IPLANE.EQ.0) THEN
          RS1(IPT) = YMIN + RATIO * (YMAX - YMIN)                              
          IF (RS1(IPT).GE.0.0) THEN                                            
            IRS1(IPT) = IPOS (RS1(IPT), YS, NY3D-1) + IMOD                     
          ELSE                                                                  
            IRS1(IPT) =-IPOS(-RS1(IPT), YS, NY3D-1)                             
          ENDIF                                                                 
        ELSE
          RS1(IPT) = XMIN + RATIO * (XMAX - XMIN)                              
          IRS1(IPT) = IPOS (RS1(IPT), XS, NXS-1) 
        ENDIF  
        COORD1(IPT) = RS1(IPT)
  100 CONTINUE                                                                  
C
C     SECOND COORDINATE MPTS - 0 =P , 1 = P , 2 =Y 
C

      DO 110 IPT = 1, MPTS                                                      
        RATIO = REAL (IPT) / REAL (MPTS+1)                                      
        IF (IPLANE.EQ.2) THEN
          RS2(IPT) = YMIN + RATIO * (YMAX - YMIN)                              
          IF (RS2(IPT).GE.0.0) THEN                                            
            IRS2(IPT) = IPOS (RS2(IPT), YS, NY3D-1) + IMOD
          ELSE                                                                  
            IRS2(IPT) =-IPOS(-RS2(IPT), YS, NY3D-1)                             
          ENDIF                                                                
        ELSE
          RS2(IPT) = PMIN + RATIO * (PMAX - PMIN)                              
          IF (RS2(IPT).GE.0) THEN        
            IRS2(IPT) = IPOS (RS2(IPT), PS, MAXNPS-1) 
          ELSE
            IRS2(IPT) = -IPOS(-RS2(IPT),PS,MAXNPS-1)
          ENDIF
        ENDIF  
        COORD2(IPT) = RS2(IPT)
  110 CONTINUE                                                                  
c
c jdemod - comment out print outs for now
c     WRITE(6,*) 'PS:',(PS(IPT),IPT=-MAXNPS,MAXNPS)
c     WRITE(6,*) 'RS1:',(RS1(IPT),IPT=1,NPTS)  
c     WRITE(6,*) 'IRS1:',(IRS1(IPT),IPT=1,NPTS) 
c     WRITE(6,*) 'RS2:',(RS2(IPT),IPT=1,NPTS)     
c     WRITE(6,*) 'IRS2:',(IRS1(IPT),IPT=1,MPTS) 
C
C     MAP THE DURFAS ARRAY ONTO AN EQUALLY SPACED GRID SPECIFIED BY
C     THE PRECEDING CALCULATIONS AND STORE IN THE SURFAS ARRAY
C        
C     THE NEGATIVE SIGN REPRESENTS ONE ADDED ON BECAUSE THE 
C     QUANTITIES ARE SUBTRACTED
C
        IF (IPLANE.EQ.0) THEN                                          
          IN1 = IYMIN-1
          IN2 = IPMIN-1                
        ELSEIF (IPLANE.EQ.1) THEN                                      
          IN1 = IXMIN-1
          IN2 = IPMIN-1
        ELSEIF (IPLANE.EQ.2) THEN
          IN1 = IXMIN-1
          IN2 = IYMIN-1
        ENDIF
c jdemod - comment out print out for now
c       WRITE(6,*) 'IN1:',IN1,' IN2:',IN2
C                                                                               
        DO 120 IRX = 1, NPTS                                                    
          DO 125 IRY = 1, MPTS  
            SURFAS(IRX,IRY)=SNGL(
     >           DURFAS(IRS1(IRX)-IN1,IRS2(IRY)-IN2))                
 125      CONTINUE
 120    CONTINUE
C
C
      ELSE     
C
C     IPLOT = 1 -- CONTOUR PLOT (SET NPTS,MPTS = NUMBER OF COORDINATES)
C                  COPY RESULTS TO SURFAS
C     
        IF (IPLANE.EQ.0) THEN 
           NPTS = IYMAX-IYMIN+1 + IMOD
           MPTS = IPMAX-IPMIN+1  
        ELSEIF (IPLANE.EQ.1) THEN 
           NPTS = IXMAX-IXMIN+1
           MPTS = IPMAX-IPMIN+1
        ELSEIF (IPLANE.EQ.2) THEN 
           NPTS = IXMAX-IXMIN+1
           MPTS = IYMAX-IYMIN+1 + IMOD
        ENDIF
C
c jdemod - comment out print outs for now 
c
c       WRITE(6,*) 'N,M:',NPTS,MPTS
C
C
c       WRITE(6,*) 'XOUTS:',(J,':',XOUTS(J),J=IXMIN,IXMAX)
c       WRITE(6,*) 'YOUTS:',(J,':',YOUTS(J),J=IYMIN,IYMAX)
c       WRITE(6,*) 'IPMIN,IPMAX:',IPMIN,IPMAX
c       WRITE(6,*) 'POUTS:',(J,':',POUTS(J),J=IPMIN,IPMAX)
c
        DO 130 IRX = 1,NPTS
            IF (IPLANE.EQ.0) THEN 
               IF (((IRX+IYMIN-1).LT.0).OR.(IMOD.EQ.0)) THEN 
                  COORD1(IRX) = YOUTS(IRX+IYMIN-1)
               ELSE
                  COORD1(IRX) = YOUTS(IRX+IYMIN)
               ENDIF 
            ELSEIF (IPLANE.EQ.1) THEN 
               COORD1(IRX) = XOUTS(IRX+IXMIN-1)
            ELSEIF (IPLANE.EQ.2) THEN 
               COORD1(IRX) = XOUTS(IRX+IXMIN-1)
            ENDIF  
c
c jdemod - comment out print out for now 
c            WRITE(6,*) 'I,Y:',IRX+IYMIN-1,YOUTS(IRX+IYMIN-1)
c            WRITE(6,*) 'C1:',IRX,COORD1(IRX)
 130    CONTINUE     


        DO 140 IRY = 1,MPTS
            IF (IPLANE.EQ.0) THEN 
               COORD2(IRY) = POUTS(IRY+IPMIN-1)
            ELSEIF (IPLANE.EQ.1) THEN 
               COORD2(IRY) = POUTS(IRY+IPMIN-1)
            ELSEIF (IPLANE.EQ.2) THEN 
               IF (((IRY+IYMIN-1).LT.0) .OR. (IMOD.EQ.0)) THEN
                  COORD2(IRY) = YOUTS(IRY+IYMIN-1)
               ELSE
                  COORD2(IRY) = YOUTS(IRY+IYMIN) 
               ENDIF
            ENDIF  
c
c jdemod - comment out print out for now
c            WRITE(6,*) 'I,P:',IRY+IPMIN-1,POUTS(IRY+IPMIN-1)
c            WRITE(6,*) 'C2: ',IRY,COORD2(IRY) 
c
 140    CONTINUE     

        DO 155 IRX = 1,NPTS
          DO 150 IRY = 1,MPTS
            SURFAS(IRX,IRY) = SNGL (DURFAS(IRX,IRY))
 150      CONTINUE
 155    CONTINUE
c
c jdemod - comment out print out for now
c slmod
c      WRITE(6,*) 'Output SURFAS:',NPTS,MPTS
c slmod end
c      WRITE(6,'(10G13.5)') ((SURFAS(IRX,IRY),
c     >             IRX=1,NPTS),IRY=1,MPTS)
     
      ENDIF

      RETURN
                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE RINTXY (VS,V3,IPLANE,MAXIZ,XYINTS,IFOLD,RV,XV,YV,IVU,          
     >                   SSS,NYSLIM,FP,FT)                                      
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    V3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)                
      REAL    XYINTS(-MAXNPS:MAXNPS,-2:MAXIZS+1),RV,XV,YV,FP,FT                 
      INTEGER IPLANE,MAXIZ,IFOLD,IVU,NYSLIM                                     
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTXY: TAKE REAL ARRAY V3 (NOT SUITABLE FOR VS ARRAY) AND       *        
C  *  INTEGRATE OVER X AND Y TO OBTAIN RESULTS ARRAY XYINTS.  THE      *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO XYINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTXYS.         *        
C  *                                                                   *        
C  *  INTEGRATION OVER Y IS PERFORMED BY ROUTINE RINTY  (AND TAKES     *        
C  *  OBSERVER'S VIEWPOINT INTO ACCOUNT IF REQUIRED).                  *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTXYS(-MAXNPS:MAXNPS,-2:MAXIZS+1),DXYSUM                
      REAL    YINTS(MAXNXS,-2:MAXIZS+1)                                         
      INTEGER IP,IZ,IX                                                          
C                                                                               
C      WRITE(6,*) 'RINTXY:START IFOLD,SSS,IVU,NYSLIM,FP',
C     >         IFOLD,SSS,IVU,NYSLIM,FP
C
      CALL DZERO (DNTXYS, (2*MAXNPS+1)*(MAXIZS+4))                              
C                                                                               
      IF (IPLANE.EQ.99) THEN                                                    
C     ======================                                                    
      DO 200 IP = -MAXNPS, MAXNPS                                               
        CALL RINTY (VS,V3,IP,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,SSS,NYSLIM,         
     >              FP,FT,1.0,-NYS/2,NYS/2)                               
        DO 190 IZ = -1, MAXIZ                                                   
          DXYSUM = 0.0D0                                                        
          DO 180 IX = 1, NXS                                                    
            DXYSUM = DXYSUM + DBLE (YINTS(IX,IZ) * XWIDS(IX))                   
  180     CONTINUE                                                              
          DXYSUM = DXYSUM / DBLE (PWIDS(IP))                                    
          DNTXYS(IP,IZ) = DXYSUM                                                
          IF (IZ.GT.0) DNTXYS(IP,MAXIZ+1) = DNTXYS(IP,MAXIZ+1) + DXYSUM         
  190   CONTINUE                                                                
  200 CONTINUE                                                                  
C                                                                               
      ELSEIF (IPLANE.EQ.0) THEN                                                 
C     =========================                                                 
        CALL RINTY (VS,V3,99,MAXIZ,YINTS,IFOLD,RV,XV,YV,IVU,SSS,NYSLIM,         
     >              FP,FT,1.0,-NYS/2,NYS/2)                                   
        DO 290 IZ = -1, MAXIZ                                                   
          DXYSUM = 0.0D0                                                        
          DO 280 IX = 1, NXS                                                    
            DXYSUM = DXYSUM + DBLE (YINTS(IX,IZ) * XWIDS(IX))                   
  280     CONTINUE                                                              
          DNTXYS(0,IZ) = DXYSUM                                                 
          IF (IZ.GT.0) DNTXYS(0,MAXIZ+1) = DNTXYS(0,MAXIZ+1) + DXYSUM           
  290   CONTINUE                                                                
C                                                                               
      ENDIF                                                                     
C     =====                                                                     
      DO 310 IP = -MAXNPS, MAXNPS                                               
        XYINTS(IP,-2) =                                                         
     >            SNGL (DBLE(FT)*DNTXYS(IP,0) - DBLE(FP)*DNTXYS(IP,-1))         
        DO 300 IZ = -1, MAXIZ+1                                                 
          XYINTS(IP,IZ) = SNGL (DNTXYS(IP,IZ))                                  
  300   CONTINUE                                                                
  310 CONTINUE                                                                  
C                                                                               
      WRITE(6,*) 'RINTXY: END'
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE RINT3D (VS,V3,ISTATE,IPLANE,IFOLD,NPTS,                        
     >                  SURFAS,XMIN,XMAX,YMIN,YMAX,NIZS,SUREDG,QEDGES,
     >                  TFTRG)          
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    V3(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS)                
      REAL    SURFAS(192,192),XMIN,XMAX,YMIN,YMAX                               
      REAL    SUREDG(192,192),QEDGES(-MAXQXS:0,2)                               
      INTEGER ISTATE,IPLANE,IFOLD,NPTS,NIZS,TFTRG
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RSET3D: THIS ROUTINE PREPARES AN ARRAY OF SURFACE HEIGHTS IN     *        
C  *  ARRAY SURFAS, FOR A GIVEN IONISATION STATE/P PLANE COMBINATION.  *        
C  *  WRITTEN FOR USE WITH 3D GRAPH OPTION, IT REQUIRES THAT THE       *        
C  *  IRREGULAR GRID OF X BINS AND Y BINS BE CONVERTED INTO A REGULAR  *        
C  *  GRID STRUCTURE, NPTS BY NPTS IN X AND Y DIRECTIONS.  THE LIMIT   *        
C  *  ON NPTS IS 192 SINCE THIS IS A GHOST80 LIMITATION.               *        
C  *  THE SUMMING IS DONE IN D.P. FOR ACCURACY.                        *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  FEBRUARY 1988                         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER  IPT,IPOS,IRXS(192),IRYS(192),IRX,IRY,IZ1,IZ2,IZ,IQX,J            
      integer  irxspq,iryspq 
      REAL     RATIO,RXS(192),RYS(192),VAL1,VAL2,X,Y(4)                         
      REAL     P,Q,RATIO1,RATIO2,XTMP,YTMP
      DOUBLE PRECISION DURFAS(192,192)                                          
C
C     THE FOLLOWING IS DONE IN TWO SEGMENTS - 1 FOR TFTR GRID MAPPED
C     PLOTS AND ONE FOR THE STANDARD CASE. 
C
      write (6,*) 'tftrg:', tftrg,xmin,xmax,ymin,ymax,ifold

      if (tftrg.eq.0) then 
 
C                                                                               
C---- SET INDEXING ARRAYS FOR GOING FROM IRREGULAR TO REGULAR GRIDS ...         
C                                                                                      
      DO 100 IPT = 1, NPTS                                                      
        RATIO = REAL (IPT) / REAL (NPTS+1)                                      
        RXS(IPT) = XMIN + RATIO * (XMAX - XMIN)                                 
        RYS(IPT) = YMIN + RATIO * (YMAX - YMIN)                                 
        IRXS(IPT)= IPOS (RXS(IPT), XS, NXS-1)                                   
        IF (RYS(IPT).GE.0.0) THEN                                               
          IRYS(IPT) = IPOS (RYS(IPT), YS, NYS-1)                                
        ELSE                                                                    
          IRYS(IPT) =-IPOS(-RYS(IPT), YS, NYS-1)                                
        ENDIF                                                                   
  100 CONTINUE                                                                  
C                                                                               
C---- CREATE LIMITER EDGE TO BE SUPERIMPOSED ...                                
C                                                                               
      CALL RZERO (SUREDG, 192*192)                                              
      DO 150 IQX = 1-NQXSO, 0                                                   
        X = QXS(IQX)                                                            
        IF (X.GE.XMIN.AND.X.LE.XMAX) THEN                                       
          IRX  = IPOS (X, RXS, NPTS-1)                                          
          Y(1) = QEDGES(IQX,2)                                                  
          Y(2) =-QEDGES(IQX,1)                                                  
          Y(3) = CL + CL - QEDGES(IQX,1)                                        
          Y(4) =-CL - CL + QEDGES(IQX,2)                                        
          DO 140 J = 1, 4                                                       
            IF (Y(J).GE.YMIN.AND.Y(J).LE.YMAX) THEN                             
              IRY = IPOS (Y(J), RYS, NPTS-1)                                    
              SUREDG(IRX,IRY) = 1.0                                             
            ENDIF                                                               
  140     CONTINUE                                                              
        ENDIF                                                                   
  150 CONTINUE                                                                  
C                                                                               
C---- CALCULATE SURFACE VALUES TO BE PLOTTED.  FOUR POSSIBILITIES:-             
C---- (1) INTEGRATED OVER P  (IPLANE=99),  Y FOLDING OFF                        
C---- (2) INTEGRATED OVER P  (IPLANE=99),  Y FOLDING ON  (IE ADD -Y             
C----     REGION RESULTS INTO +Y REGION AND VICE-VERSA)                         
C---- (3) PARTICULAR P PLANE  (-4 TO 4), Y FOLDING OFF                          
C---- (4) PARTICULAR P PLANE  (-4 TO 4), Y FOLDING ON                           
C---- ON TOP OF THIS, TO GET "ALL IONISATIONS" ONE MUST ADD FOR                 
C---- IZ = 1:NIZS, OTHERWISE WE JUST WANT STATE "ISTATE".                       
C---- THE SUMMATION IS DONE IN DOUBLE PRECISION FOR ACCURACY ...                
C                                                                               
      CALL DZERO (DURFAS, 192*192)                                              
      IF (ISTATE.EQ.NIZS+1) THEN                                                
        IZ1 = 1                                                                 
        IZ2 = NIZS                                                              
      ELSE                                                                      
        IZ1 = ISTATE                                                            
        IZ2 = ISTATE                                                            
      ENDIF                                                                     
C                                                                               
      DO 220 IZ = IZ1, IZ2                                                      
        DO 210 IRX = 1, NPTS                                                    
          DO 200 IRY = 1, NPTS                                                  
            IF (IPLANE.EQ.99) THEN                                              
              VAL1 = VS(IRXS(IRX),IRYS(IRY),IZ)                                 
              VAL2 = VS(IRXS(IRX),-IRYS(IRY),IZ)                                
            ELSE                                                                
              VAL1 = V3(IRXS(IRX),IRYS(IRY),IZ,IPLANE)                          
              VAL2 = V3(IRXS(IRX),-IRYS(IRY),IZ,IPLANE)                         
            ENDIF                                                               
            IF (IFOLD.EQ.0) THEN                                                
              DURFAS(IRX,IRY) = DURFAS(IRX,IRY) + DBLE (VAL1)                   
            ELSE                                                                
              DURFAS(IRX,IRY) = DURFAS(IRX,IRY) + DBLE (0.5*(VAL1+VAL2))        
            ENDIF                                                               
  200     CONTINUE                                                              
  210   CONTINUE                                                                
  220 CONTINUE                                                                  
C
      ELSEIF (TFTRG.EQ.1) THEN       
C
C                                                                               
C---- CALCULATE SURFACE VALUES TO BE PLOTTED.  FOUR POSSIBILITIES:-             
C---- (1) INTEGRATED OVER P  (IPLANE=99),  Y FOLDING OFF                        
C---- (2) INTEGRATED OVER P  (IPLANE=99),  Y FOLDING ON  (IE ADD -Y             
C----     REGION RESULTS INTO +Y REGION AND VICE-VERSA)                         
C---- (3) PARTICULAR P PLANE  (-4 TO 4), Y FOLDING OFF                          
C---- (4) PARTICULAR P PLANE  (-4 TO 4), Y FOLDING ON                           
C---- ON TOP OF THIS, TO GET "ALL IONISATIONS" ONE MUST ADD FOR                 
C---- IZ = 1:NIZS, OTHERWISE WE JUST WANT STATE "ISTATE".                       
C---- THE SUMMATION IS DONE IN DOUBLE PRECISION FOR ACCURACY ...                
C                                                                               
      CALL RZERO (SUREDG, 192*192)                                              
      CALL DZERO (DURFAS, 192*192)                                              
      IF (ISTATE.EQ.NIZS+1) THEN                                                
        IZ1 = 1                                                                 
        IZ2 = NIZS                                                              
      ELSE                                                                      
        IZ1 = ISTATE                                                            
        IZ2 = ISTATE                                                            
      ENDIF                                                                     
C                                                                               
      DO 520 IRX = 1, NPTS                                                    
        RATIO1 = REAL (IRX) / REAL (NPTS+1)                                      
        P = XMIN + RATIO1 * (XMAX - XMIN)                                 
        DO 510 IRY = 1, NPTS                                                  
          RATIO2 = REAL (IRY) / REAL (NPTS+1)
          Q = YMIN + RATIO2 * (YMAX - YMIN)                                                     
          RXS(IRX) = P
          RYS(IRY) = Q 
          XTMP = CA - SQRT((CA-P)**2+Q**2)
          YTMP = CA * ATAN2(Q,CA-P)                  
          IRXSPQ= IPOS (XTMP, XS, NXS-1)                                   
          IF (YTMP.GE.0.0) THEN                                               
            IRYSPQ = IPOS (YTMP, YS, NYS-1)                                
          ELSE                                                                    
            IRYSPQ =-IPOS(-YTMP, YS, NYS-1)                                
          ENDIF                                                                   

          write (6,*) irx,iry,p,q,xtmp,ytmp,
     >                irxspq,iryspq

          DO 500 IZ = IZ1, IZ2                                                      
            IF (IPLANE.EQ.99) THEN                                              
              VAL1 = VS(IRXSPQ,IRYSPQ,IZ)                                 
              VAL2 = VS(IRXSPQ,-IRYSPQ,IZ)                                
            ELSE                                                                
              VAL1 = V3(IRXSPQ,IRYSPQ,IZ,IPLANE)                          
              VAL2 = V3(IRXSPQ,-IRYSPQ,IZ,IPLANE)                         
            ENDIF                                                               
            IF (IFOLD.EQ.0) THEN                                                
              DURFAS(IRX,IRY) = DURFAS(IRX,IRY) + DBLE (VAL1)                   
            ELSE                                                                
              DURFAS(IRX,IRY) = DURFAS(IRX,IRY) + DBLE (0.5*(VAL1+VAL2))        
            ENDIF                                                               
  500     CONTINUE                                                              
  510   CONTINUE                                                                
  520 CONTINUE                                                                  
C
      ENDIF
C                                                                               
      DO 310 IRX = 1, NPTS                                                      
        DO 300 IRY = 1, NPTS                                                    
          SURFAS(IRX,IRY) = SNGL (DURFAS(IRX,IRY))                              
  300   CONTINUE                                                                
  310 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C
C=======================================================================
C
      SUBROUTINE RSMO3D(SURFAS,NPTS)
      IMPLICIT none
      INTEGER NPTS
      REAL SURFAS(192,192)  
C
C     RSMO3D: This routine applies a 2-D 3 X 3 matrix convolution to 
C             the input array surfas. The values of this function 
C             are set in the data statement and are currently all
C             set to 1.0. This will generate a simple weighted 
C             average smoothing of the surface.
C
C             The temporary array is used to store the modified 
C             graph data. A more memory efficient implementation
C             would be to use only one or two rows to hold the 
C             modified data as it workd through the array.
C
C             David Elder   Mar 2, 1992
C
c      automatic surtmp
c
      real surtmp(192,192) 
      real mw(3,3),totmw
      data mw/1.0, 1.0, 1.0,
     >        1.0, 1.0, 1.0,
     >        1.0, 1.0, 1.0/
      integer i,j    
C
C     Loop through centre segments of surfas - deal with edges 
C     separately.  
C
      totmw = 0
      do 10 i = 1,3
        do 10 j = 1,3
          totmw = totmw + mw(i,j)
10    continue

C
      do 100 i = 2,npts-1
C
        do 100 j = 2,npts-1
C
C       Leave the edges for later - see how it looks - if this 
C       is adequate then it may not be worthwhile smoothing the 
C       outermost entries.
C
        surtmp(i,j) =(surfas(i-1,j-1)*mw(1,1)+
     >                surfas(i-1,j)*mw(1,2)+
     >                surfas(i-1,j+1)*mw(1,3)+
     >                surfas(i,j-1)*mw(2,1)+
     >                surfas(i,j)*mw(2,2)+
     >                surfas(i,j+1)*mw(2,3)+
     >                surfas(i+1,j-1)*mw(3,1)+
     >                surfas(i+1,j)*mw(3,2)+
     >                surfas(i+1,j+1)*mw(3,3))/totmw
100   continue     

      do 200 i = 2,npts-1
        do 200 j = 2,npts-1
          surfas(i,j) = surtmp(i,j)
200   continue
      return
      end     
C                                                                               
C=======================================================================        
C                                                                               
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE TINTX (V5,ISTATE,IPLANE,MAXIT,XJNTS,IFOLD,RV,XV,YV,IVU,        
     >                  SSS,NYSLIM)                                             
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                                                                               
      REAL    V5(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS,MAXNTS)         
      REAL    XJNTS(-MAXNYS:MAXNYS,MAXNTS+1),RV,XV,YV                           
      INTEGER ISTATE,IPLANE,MAXIT,IFOLD,IVU,NYSLIM                              
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  TINTX:  TAKE THE REAL ARRAY V5.  IF IPLANE=99, THEN INTEGRATE    *        
C  *  OVER P, OTHERWISE FOCUS ON V5(,,,IPLANE,) VALUES.  THEN          *        
C  *  INTEGRATE OVER X DIMENSION TO OBTAIN RESULTS ARRAY XINTS.  THE   *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  ALL SUMMING FOR A GIVEN CHARGE STATE  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTXS.          *        
C  *                                                                   *        
C  *  VIEWPOINTS ADDED - JULY 88.                                      *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTXS(-MAXNYS:MAXNYS,MAXNTS+1),DXSUM1,DXSUM2             
      INTEGER IP,IT,IY,IX,JY,KY,LY,MYSLIM                                       
      REAL    FRAC,TMP                                                          
C                                                                               
      CALL DZERO (DNTXS, (2*MAXNYS+1)*(MAXNTS+1))                               
C                                                                               
      MYSLIM = MIN (NYSLIM, NY3D)                                               
      IF (IVU.GT.0) THEN                                                        
        CALL VIEW2 (RV,XV,YV,SSS,MYSLIM,IVU)                                    
        DO 230 IY = -MYSLIM, MYSLIM                                             
          IF (IY.EQ.0) GOTO 230                                                 
          JY = IABS (IY)                                                        
          DO 220 IX = 1, NXS                                                    
            IF (YPPPP(IX,IY).GT.CTWOL .OR. YPPPP(IX,IY).LT.-CTWOL .OR.          
     >         (XPPPP(IX,IY).LT.0.0.AND.IY.LT.0)) GOTO 220                      
            DO 210 KY = -NY3D, NY3D                                             
              IF (KY.EQ.0) GOTO 210                                             
              LY = IABS (KY)                                                    
              IF (KY.LT.0) THEN                                                 
                TMP = FRAC (YWIDS(LY)-YS(LY),YWIDS(LY),                         
     >                                       YPPPP(IX,IY),YPWID(IX,IY))         
              ELSE                                                              
                TMP = FRAC (YS(LY),YWIDS(LY),YPPPP(IX,IY),YPWID(IX,IY))         
              ENDIF                                                             
              IF (TMP.LE.0.0) GOTO 210                                          
              TMP = TMP * XWIDS(IX) * YWIDS(JY) / YWIDS(LY) / SSQ(IX,IY)        
              DO 200 IT = 1, MAXIT                                              
               IF (IPLANE.EQ.99) THEN                                           
                DO 190 IP = -MAXNPS, MAXNPS                                     
                 DNTXS(KY,IT) = DNTXS(KY,IT) +                                  
     >                              DBLE (TMP * V5(IX,IY,ISTATE,IP,IT))         
  190           CONTINUE                                                        
               ELSE                                                             
                DNTXS(KY,IT) = DNTXS(KY,IT) +                                   
     >                           DBLE (TMP * V5(IX,IY,ISTATE,IPLANE,IT))        
               ENDIF                                                            
  200         CONTINUE                                                          
  210       CONTINUE                                                            
  220     CONTINUE                                                              
  230   CONTINUE                                                                
C                                                                               
      ELSEIF (IVU.EQ.0) THEN                                                    
      DO 130 IT = 1, MAXIT                                                      
       DO 120 IY = 1, NY3D                                                      
        DXSUM1 = 0.0D0                                                          
        DXSUM2 = 0.0D0                                                          
        DO 110 IX = 1, NXS                                                      
         IF (IPLANE.EQ.99) THEN                                                 
          DO 100 IP = -MAXNPS, MAXNPS                                           
           DXSUM1 = DXSUM1 + DBLE (V5(IX,IY,ISTATE,IP,IT) * XWIDS(IX))          
           DXSUM2 = DXSUM2 + DBLE (V5(IX,-IY,ISTATE,IP,IT) * XWIDS(IX))         
  100     CONTINUE                                                              
         ELSE                                                                   
          DXSUM1 = DXSUM1 + DBLE (V5(IX,IY,ISTATE,IPLANE,IT)* XWIDS(IX))        
          DXSUM2 = DXSUM2 + DBLE (V5(IX,-IY,ISTATE,IPLANE,IT)*XWIDS(IX))        
         ENDIF                                                                  
  110   CONTINUE                                                                
        IF (IFOLD.EQ.0) THEN                                                    
         DNTXS(IY,IT)  = DXSUM1                                                 
         DNTXS(-IY,IT) = DXSUM2                                                 
        ELSEIF (IFOLD.EQ.1) THEN                                                
         DNTXS(IY,IT)  = 0.5 * (DXSUM1 + DXSUM2)                                
         DNTXS(-IY,IT) = 0.5 * (DXSUM1 + DXSUM2)                                
        ENDIF                                                                   
  120  CONTINUE                                                                 
  130 CONTINUE                                                                  
      ENDIF                                                                     
C                                                                               
      DO 300 IT = 1, MAXIT                                                      
        DNTXS(0,IT) = 0.5 * (DNTXS(1,IT) + DNTXS(-1,IT))                        
  300 CONTINUE                                                                  
      DO 330 IY = -NYS, NYS                                                     
        DO 310 IT = 1, MAXIT                                                    
          DNTXS(IY,MAXIT+1) = DNTXS(IY,MAXIT+1) + DNTXS(IY,IT)                  
  310   CONTINUE                                                                
        DO 320 IT = 1, MAXIT+1                                                  
          XJNTS(IY,IT) = SNGL (DNTXS(IY,IT))                                    
  320   CONTINUE                                                                
  330 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE TINTY (V5,ISTATE,IPLANE,MAXIT,YJNTS,IFOLD,RV,XV,YV,IVU,        
     >                  SSS,NYSLIM)                                             
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                                                                               
      REAL    V5(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS,MAXNTS)         
      REAL    YJNTS(MAXNXS,MAXNTS+1),RV,XV,YV                                   
      INTEGER ISTATE,IPLANE,MAXIT,IFOLD,IVU,NYSLIM                              
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  TINTY:  TAKE THE REAL ARRAY V5.  IF IPLANE=99, THEN INTEGRATE    *        
C  *  OVER P, OTHERWISE FOCUS ON V5(,,,IPLANE,) VALUES.  THEN          *        
C  *  INTEGRATE OVER Y DIMENSION TO OBTAIN RESULTS ARRAY YJNTS.  THE   *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL TIME     *        
C  *  POINTS AND WRITTEN TO YJNTS(,MAXIT+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTYS.          *        
C  *                                                                   *        
C  *  INTEGRATION IS EITHER STRAIGHTFORWARD, OR BY METHOD OF NOTE 91   *        
C  *  WHICH TAKES INTO ACCOUNT AN OBSERVER'S VIEWPOINT RELATIVE TO THE *        
C  *  Y AXIS.  WE ONLY USE VALUES FROM -RV/3 TO RV/3  (OR FROM -L TO L *        
C  *  IF SMALLER).  POINTS WHICH WOULD BE PROJECTED ON TO THE X AXIS   *        
C  *  BEYOND THE WALL OR INSIDE THE CENTRE ARE IGNORED.  CONTRIBUTIONS *        
C  *  FROM -Y PROJECTED ONTO X' AXIS <0 ARE IGNORED  (SINCE THE LINE-  *        
C  *  OF-SIGHT IS BLOCKED BY THE LIMITER).  TAKE INTO ACCOUNT ANY      *        
C  *  DIFFERENCES IN BIN-WIDTHS FROM X --> X' SYSTEMS.                 *        
C  *  MODIFICATION FROM NOTE 172: SCALE PLRPS BY DISTANCE SQUARED SSQ  *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTYS(MAXNXS,MAXNTS+1),DYSUM                             
      INTEGER IT,IP,IY,IX,JY,JX,MYSLIM                                          
      REAL    FRAC,TMP                                                          
C                                                                               
      CALL DZERO (DNTYS, MAXNXS*(MAXNTS+1))                                     
C                                                                               
      IF (IVU.GT.0) THEN                                                        
        MYSLIM = MIN (NYSLIM, NY3D)                                             
        CALL VIEW2 (RV,XV,YV,SSS,MYSLIM,IVU)                                    
        DO 230 IY = -MYSLIM, MYSLIM                                             
          IF (IY.EQ.0) GOTO 230                                                 
          JY = IABS (IY)                                                        
          DO 220 IX = 1, NXS                                                    
            IF (XPPPP(IX,IY).GT.CA .OR. XPPPP(IX,IY).LT.CAW .OR.                
     >         (XPPPP(IX,IY).LT.0.0.AND.IY.LT.0)) GOTO 220                      
            DO 210 JX = 1, NXS                                                  
C             TMP = FRAC (XS(JX),XWIDS(JX),XPPPP(IX,IY),XWIDS(IX))              
              TMP = FRAC (XS(JX),XWIDS(JX),XPPPP(IX,IY),XPWID(IX,IY))           
              IF (TMP.LE.0.0) GOTO 210                                          
              TMP = TMP * YWIDS(JY) * XWIDS(IX) / XWIDS(JX) / SSQ(IX,IY)        
              DO 205 IT = 1, MAXIT                                              
               IF (IPLANE.EQ.99) THEN                                           
                DO 200 IP = -MAXNPS, MAXNPS                                     
                  DNTYS(JX,IT) = DNTYS(JX,IT) +                                 
     >                           DBLE (TMP * V5(IX,IY,ISTATE,IP,IT))            
  200           CONTINUE                                                        
               ELSE                                                             
                DNTYS(JX,IT) = DNTYS(JX,IT) +                                   
     >                           DBLE (TMP * V5(IX,IY,ISTATE,IPLANE,IT))        
               ENDIF                                                            
  205         CONTINUE                                                          
  210       CONTINUE                                                            
  220     CONTINUE                                                              
  230   CONTINUE                                                                
C                                                                               
      ELSEIF (IVU.EQ.0) THEN                                                    
        DO 270 IT = 1, MAXIT                                                    
          DO 260 IX = 1, NXS                                                    
            DYSUM = 0.0                                                         
            DO 250 IY = 1, NY3D                                                 
              IF (IPLANE.EQ.99) THEN                                            
                DO 240 IP = -MAXNPS, MAXNPS                                     
                  DYSUM = DYSUM + DBLE (YWIDS(IY) *                             
     >              (V5(IX,IY,ISTATE,IP,IT) + V5(IX,-IY,ISTATE,IP,IT)))         
  240           CONTINUE                                                        
              ELSE                                                              
                DYSUM = DYSUM + DBLE (YWIDS(IY) *                               
     >       (V5(IX,IY,ISTATE,IPLANE,IT) + V5(IX,-IY,ISTATE,IPLANE,IT)))        
              ENDIF                                                             
  250       CONTINUE                                                            
            DNTYS(IX,IT) = DYSUM                                                
  260     CONTINUE                                                              
  270   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
      DO 320 IX = 1, NXS                                                        
        DO 300 IT = 1, MAXIT                                                    
          DNTYS(IX,MAXIT+1) = DNTYS(IX,MAXIT+1) + DNTYS(IX,IT)                  
  300   CONTINUE                                                                
        DO 310 IT = 1, MAXIT+1                                                  
          YJNTS(IX,IT) = SNGL (DNTYS(IX,IT))                                    
  310   CONTINUE                                                                
  320 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE TINTXY (V5,ISTATE,IPLANE,MAXIT,XYJNTS,IFOLD,RV,XV,YV,          
     >                   IVU,SSS,NYSLIM)                                        
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      REAL    V5(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS,MAXNTS)         
      REAL    XYJNTS(-MAXNPS:MAXNPS,MAXNTS+1),RV,XV,YV                          
      INTEGER ISTATE,IPLANE,MAXIT,IFOLD,IVU,NYSLIM                              
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  TINTXY: TAKE THE REAL ARRAY V5.  IF IPLANE=99, THEN INTEGRATE    *        
C  *  OVER P, OTHERWISE FOCUS ON V5(,,,IPLANE,) VALUES.  THEN          *        
C  *  INTEGRATE OVER X AND Y TO OBTAIN RESULTS ARRAY XYJNTS.  THE      *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL TIME     *        
C  *  POINTS AND WRITTEN TO XYJNTS(,MAXIT+1) LOCATIONS. ALL SUMMING IS *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTXYS.         *        
C  *                                                                   *        
C  *  INTEGRATION OVER Y IS PERFORMED BY ROUTINE TINTY  (AND TAKES     *        
C  *  OBSERVER'S VIEWPOINT INTO ACCOUNT IF REQUIRED).                  *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTXYS(-MAXNPS:MAXNPS,MAXNTS+1),DXYSUM                   
      REAL    YJNTS(MAXNXS,MAXNTS+1)                                            
      INTEGER IT,IP,IX                                                          
C                                                                               
      CALL DZERO (DNTXYS, (2*MAXNPS+1)*(MAXNTS+1))                              
C                                                                               
      DO 200 IP = -MAXNPS, MAXNPS                                               
        CALL TINTY (V5,ISTATE,IP,MAXIT,YJNTS,IFOLD,RV,XV,YV,IVU,                
     >              SSS,NYSLIM)                                                 
        DO 190 IT = 1, MAXIT                                                    
          DXYSUM = 0.0D0                                                        
          DO 180 IX = 1, NXS                                                    
            DXYSUM = DXYSUM + DBLE (YJNTS(IX,IT) * XWIDS(IX))                   
  180     CONTINUE                                                              
          DXYSUM = DXYSUM / DBLE (PWIDS(IP))                                    
          DNTXYS(IP,IT) = DXYSUM                                                
          DNTXYS(IP,MAXIT+1) = DNTXYS(IP,MAXIT+1) + DXYSUM                      
  190   CONTINUE                                                                
  200 CONTINUE                                                                  
C                                                                               
      DO 310 IT = 1, MAXIT+1                                                    
        DO 300 IP = -MAXNPS, MAXNPS                                             
          XYJNTS(IP,IT) = SNGL (DNTXYS(IP,IT))                                  
  300   CONTINUE                                                                
  310 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE TINTT (V5,IPLANE,MAXIZ,QTS,TVALS,IFOLD,RV,XV,YV,IVU,           
     >                  SSS,NYSLIM,FP,FT,MAXQTS,NQTS)                           
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
C                                                                               
      INTEGER IPLANE,MAXIZ,IFOLD,IVU,NYSLIM,MAXQTS,NQTS                         
      REAL    V5(MAXNXS,-MAXY3D:MAXY3D,-1:MAXIZS,-MAXNPS:MAXNPS,MAXNTS)         
      REAL    QTS(0:MAXQTS),TVALS(0:MAXQTS,-2:MAXIZS+1),RV,XV,YV,FP,FT          
      LOGICAL SSS                                                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  TINTT:  TAKE THE REAL ARRAY V5.  IF IPLANE=99, THEN INTEGRATE    *        
C  *  OVER P, OTHERWISE FOCUS ON V5(,,,IPLANE,) VALUES.  THEN          *        
C  *  INTEGRATE OVER X AND Y TO OBTAIN RESULTS ARRAY XYJNTS USING      *        
C  *  ROUTINE TINTXY.  THE TWO LOCAL ARRAYS TS AND TLIMS ARE SET TO    *        
C  *  HOLD THE TIMEPOINTS AND INTEGRATED DENSITIES AT THESE TIMEPOINTS *        
C  *  FOR THE CURRENT IONISATION STATE.  THESE FIGURES ARE THEN        *        
C  *  LINEARLY INTERPOLATED ON TO THE REGULAR GRID QTS TO GIVE THE     *        
C  *  RESULTS ARRAY TVALS, WHERE THE TVALS(,MAXIZ+1) ENTRIES RELATE TO *        
C  *  THE TOTAL OVER ALL IONISATION STATES.  IF NEUT WAS USED, THEN    *        
C  *  AT TIME 0 WE WOULD HAVE HAD 1.0 NEUTRALS IN THE PLASMA; OTHER-   *        
C  *  WISE IF IONS WERE INJECTED IN STATE IZ THEN WE WOULD HAVE HAD    *        
C  *  1.0 OF THESE AT TIME 0.                                          *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)  MARCH 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      REAL    XYJNTS(-MAXNPS:MAXNPS,MAXNTS+1)                                   
      REAL    TS(0:MAXNTS+1),TLIMS(0:MAXNTS+1),VAL                              
      INTEGER IZ,IT,IP,IQT,IPOS                                                 
      LOGICAL FLAG                                                              
C                                                                               
      FLAG = .FALSE.                                                            
      CALL RZERO (TVALS, (MAXQTS+1)*(MAXIZS+4))                                 
C     WRITE (6,*) ' TINTT CALLED ...'                                           
C                                                                               
      DO 500 IZ = -1, MAXIZ                                                     
        CALL TINTXY (V5,IZ,IPLANE,NTS,XYJNTS,IFOLD,RV,XV,YV,IVU,                
     >               SSS,NYSLIM)                                                
C                                                                               
        DO 200 IT = 1, NTS                                                      
          IF (IPLANE.EQ.99) THEN                                                
            TLIMS(IT) = 0.0                                                     
            DO 100 IP = -MAXNPS, MAXNPS                                         
              TLIMS(IT) = TLIMS(IT) + XYJNTS(IP,IT) * PWIDS(IP)                 
  100       CONTINUE                                                            
          ELSE                                                                  
            TLIMS(IT) = XYJNTS(IPLANE,IT) * PWIDS(IPLANE)                       
          ENDIF                                                                 
          IF (IZ.EQ.-1) THEN                                                    
            TS(IT) = DWELTS(0) * DWELFS(IT)                                     
          ELSE                                                                  
            TS(IT) = DWELTS(IZ) * DWELFS(IT)                                    
          ENDIF                                                                 
  200   CONTINUE                                                                
C                                                                               
        TS(0) = 0.0                                                             
        IF (IZ.LE.0.AND.TLIMS(1).GT.0.0) THEN                                   
          TLIMS(0) = 1.0                                                        
          FLAG = .TRUE.                                                         
        ELSEIF ((IZ.EQ.CIZSC).AND.(.NOT.FLAG)) THEN                             
          TLIMS(0) = 1.0                                                        
        ELSE                                                                    
          TLIMS(0) = 0.0                                                        
        ENDIF                                                                   
        TS(NTS+1)    = 100.0                                                    
        TLIMS(NTS+1) = 0.0                                                      
C       WRITE (6,*) ' IZ,TS   ',IZ,(TS(IT),IT=0,NTS+1)                          
C       WRITE (6,*) ' IZ,TLIMS',IZ,(TLIMS(IT),IT=0,NTS+1)                       
C                                                                               
        TVALS(0,IZ) = TLIMS(0)                                                  
        IF (IZ.GT.0) TVALS(0,MAXIZ+1) = TVALS(0,MAXIZ+1) + TLIMS(0)             
        DO 300 IQT = 1, NQTS                                                    
          IT  = IPOS (QTS(IQT), TS, NTS+1) - 1                                  
          VAL = TLIMS(IT-1) + (QTS(IQT)-TS(IT-1)) *                             
     >          (TLIMS(IT)-TLIMS(IT-1)) / (TS(IT)-TS(IT-1))                     
C         WRITE (6,*) ' IZ,IQT,IT,VAL',IZ,IQT,IT,VAL                            
          TVALS(IQT,IZ) = VAL                                                   
          IF (IZ.GT.0) TVALS(IQT,MAXIZ+1) = TVALS(IQT,MAXIZ+1) + VAL            
  300   CONTINUE                                                                
C                                                                               
  500 CONTINUE                                                                  
C                                                                               
      DO 600 IQT = 0, NQTS                                                      
        TVALS(IQT,-2) = FT * TVALS(IQT,0) - FP * TVALS(IQT,-1)                  
  600 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      REAL FUNCTION FRAC (UP,BINWID,PPPP,WID)                                   
      REAL UP,BINWID,PPPP,WID                                                   
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  FRAC:   WORKS OUT WHAT PROPORTION OF THE PARCEL (HIGH,LOW)       *        
C  *  IS CONTAINED WITHIN THE BIN  (UP,DOWN).  THE PARCEL IS DEFINED   *        
C  *  BY ITS MIDPOINT PPPP AND WIDTH WID;  THE BIN IS DEFINED BY ITS   *        
C  *  UPPER BOUND UP AND WIDTH BINWID.                                 *        
C  *  THE REGIONS ARE DEFINED IN THIS WAY TO ENHANCE FLEXIBILITY.      *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)   JULY 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      REAL DOWN,HIGH,LOW                                                        
C                                                                               
      DOWN = UP - BINWID                                                        
      HIGH = PPPP + 0.5 * WID                                                   
      LOW  = PPPP - 0.5 * WID                                                   
C                                                                               
      IF (UP.GT.HIGH) THEN                                                      
C     ~~~~~~~~~~~~~~~~~~~~                                                      
        IF (DOWN.GT.LOW) THEN                                                   
          FRAC = (HIGH-DOWN) / (HIGH-LOW)                                       
        ELSE                                                                    
          FRAC = 1.0                                                            
        ENDIF                                                                   
C                                                                               
      ELSEIF (UP.GT.LOW) THEN                                                   
C     ~~~~~~~~~~~~~~~~~~~~~~~                                                   
        IF (DOWN.GT.LOW) THEN                                                   
          FRAC = (UP-DOWN) / (HIGH-LOW)                                         
        ELSE                                                                    
          FRAC = (UP-LOW) / (HIGH-LOW)                                          
        ENDIF                                                                   
C                                                                               
      ELSE                                                                      
C     ~~~~                                                                      
        FRAC = 0.0                                                              
      ENDIF                                                                     
C                                                                               
C     IF (FRAC.GT.0.0)                                                          
C    >  WRITE (6,'('' FRAC: UP,DOWN,HIGH,LOW,FRAC'',5F9.4)')                    
C    >    UP,DOWN,HIGH,LOW,FRAC                                                 
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE VIEW2 (RV,XV,YV,SSS,MYSLIM,IVU)                                
      REAL    RV,XV,YV                                                          
      LOGICAL SSS                                                               
      INTEGER MYSLIM,IVU                                                        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  VIEW2:  CALCULATES PROJECTION OF EACH OF THE FOUR CORNERS OF A   *        
C  *  PARCEL ONTO THE X AXIS AND THE Y AXIS.  THE VALUE XPPPP IS       *        
C  *  ASSIGNED TO THE MIDPOINT OF THE PROJECTION ONTO THE X AXIS, AND  *        
C  *  XPWID IS SET TO THE WIDTH OF THE PARCEL WHEN PROJECTED ON TO THE *        
C  *  X AXIS.  SIMILAR CONSIDERATIONS APPLY TO THE CALCULATION OF      *        
C  *  YPPPP AND YPWID.   SSQ IS A DISTANCE-SQUARED FACTOR TO BE        *        
C  *  APPLIED WHEN DEALING WITH PLRPS.  IT IS CALCULATED AT EACH OF    *        
C  *  THE FOUR CORNERS AND THEN AVERAGED.                              *        
C  *  PERCNT IS A MODIFYING FACTOR FOR THE CALCULATED MAXIMUM EXTENT   *        
C  *  OF THE PROJECTION ONTO THE AXES - NORMALLY 100%.                 *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)   JULY 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IY,IX,J                                                           
      REAL    YP(4),XP(4),XPP,YPP,XPPP,YPPP,XMIN,XMAX,YMIN,YMAX,SSUM            
      REAL    PERCNT                                                            
C                                                                               
      PERCNT = 1.0                                                              
      IF (IVU.EQ.3) PERCNT = 0.9                                                
      IF (IVU.EQ.4) PERCNT = 0.8                                                
      IF (IVU.EQ.5) PERCNT = 0.7                                                
C                                                                               
      DO 230 IY = -MYSLIM, MYSLIM                                               
        IF (IY.LT.0) THEN                                                       
          YP(1) = -YS(-IY)                                                      
          YP(2) = -YS(-IY) + YWIDS(-IY)                                         
        ELSEIF (IY.EQ.0) THEN                                                   
          GOTO 230                                                              
        ELSE                                                                    
          YP(1) = YS(IY) - YWIDS(IY)                                            
          YP(2) = YS(IY)                                                        
        ENDIF                                                                   
        YP(3) = YP(2)                                                           
        YP(4) = YP(1)                                                           
C                                                                               
        DO 210 IX = 1, NXS                                                      
          XP(1) = XS(IX)                                                        
          XP(2) = XP(1)                                                         
          XP(3) = XS(IX) - XWIDS(IX)                                            
          XP(4) = XP(3)                                                         
          XMIN = 1.E10                                                          
          XMAX =-1.E10                                                          
          YMIN = 1.E10                                                          
          YMAX =-1.E10                                                          
          SSUM = 0.0                                                            
C                                                                               
          DO 200 J = 1, 4                                                       
            XPP = (RV + XP(J)) * COS (YP(J) / RV)                               
            YPP = (RV + XP(J)) * SIN (YP(J) / RV)                               
            XPPP = (XPP-RV) - YPP * (XV-XPP) / (YV-YPP)                         
            YPPP = YPP - (XPP-RV) * (YV-YPP) / (XV-XPP)                         
            XMIN = MIN (XMIN, XPPP)                                             
            XMAX = MAX (XMAX, XPPP)                                             
            YMIN = MIN (YMIN, YPPP)                                             
            YMAX = MAX (YMAX, YPPP)                                             
            SSUM = SSUM + (YV-YP(J))**2 + (XV-RV-XP(J))**2                      
  200     CONTINUE                                                              
C                                                                               
          XPPPP(IX,IY) = 0.5 * (XMIN + XMAX)                                    
          YPPPP(IX,IY) = 0.5 * (YMIN + YMAX)                                    
          XPWID(IX,IY) = PERCNT * (XMAX - XMIN)                                 
          YPWID(IX,IY) = PERCNT * (YMAX - YMIN)                                 
          IF (SSS) THEN                                                         
            SSQ(IX,IY) = 0.25 * SSUM                                            
          ELSE                                                                  
            SSQ(IX,IY) = 1.0                                                    
          ENDIF                                                                 
  210   CONTINUE                                                                
  230 CONTINUE                                                                  
      RETURN                                                                    
      END                                                                       
C                                                                               
C=======================================================================        
C                                                                               
      SUBROUTINE PROJEC (MAXIZ,XINTS,YOUTS,YWIDSS,TV,SV,TC,SC,GC,IVU,           
     >                   CSINTB)                                                
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      REAL    XINTS(-MAXNYS:MAXNYS,-2:MAXIZS+1),YWIDSS(-MAXNYS:MAXNYS)          
      REAL    YOUTS(-MAXNYS:MAXNYS),TV,SV,TC,SC,GC,CSINTB                       
      INTEGER MAXIZ,IVU                                                         
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  PROJEC:  CALCULATES THE BIN-WIDTHS CORRECTION TERMS TO BE        *        
C  *  APPLIED IN THE TRANSFORMATION BETWEEN THE Y AXIS AND THE T AXIS  *        
C  *  AND ROTATES AND TRANSLATES THE YOUTS ARRAY TO CORRESPOND TO THE  *        
C  *  T AXIS.  ENSURE THAT DESCENDING ORDER OF T POINTS IS PRESERVED   *        
C  *  DURING THIS PROCESS.                                             *        
C  *  AS A FINAL STEP, INVERTS ALL ARRAYS SO THE T POINTS REVERT TO    *        
C  *  ASCENDING ORDER - THIS IS REQUIRED IN CASE SMOOTHING HAS BEEN    *        
C  *  SPECIFIED, OTHERWISE THE SMOOTHING ROUTINE CALLED FROM DRAW      *        
C  *  FAILS.                                                           *        
C  *  NOTE 299 - MULTIPLIES Y VALUES BY THE SIN(THETAB) FACTOR.        *        
C  *                                                                   *        
C  *                                                                   *        
C  *  C.M.FARRELL  (HUNTERSKIL)   JULY 1988                            *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IY,JY,IZ                                                          
      REAL    COSGC,SINGC,FACTOR,YL,YU,TL,TU,TEMP                               
C                                                                               
      DO 60 IY = -MAXNYS, MAXNYS                                                
        DO 50 IZ = -2, MAXIZ+1                                                  
          XINTS(IY,IZ) = XINTS(IY,IZ) * CSINTB                                  
   50   CONTINUE                                                                
   60 CONTINUE                                                                  
C                                                                               
      IF (IVU.LT.2) RETURN                                                      
C                                                                               
      COSGC = COS (GC)                                                          
      SINGC = SIN (GC)                                                          
      YOUTS(0) = SV * (TV - TC) / (SC - SV) + TV                                
C                                                                               
      DO 110 IY = -MAXNYS, MAXNYS                                               
        IF (IY.EQ.0) GOTO 110                                                   
        YL = YOUTS(IY) - 0.5 * YWIDSS(IY)                                       
        YU = YOUTS(IY) + 0.5 * YWIDSS(IY)                                       
        TL = SV * (TV - TC + YL*COSGC) / (SC - YL*SINGC - SV) + TV              
        TU = SV * (TV - TC + YU*COSGC) / (SC - YU*SINGC - SV) + TV              
        YOUTS(IY) = 0.5 * (TU + TL)                                             
        IF (IY.GT.-MAXNYS .AND. YOUTS(IY).GE.YOUTS(IY-1)) THEN                  
          YOUTS(IY) = YOUTS(IY-1) - 0.001                                       
          FACTOR = 0.0                                                          
          YWIDSS(IY) = 0.001                                                    
        ELSE                                                                    
          FACTOR = (YU - YL) / (TL - TU)                                        
          YWIDSS(IY) = TL - TU                                                  
        ENDIF                                                                   
        DO 100 IZ = -2, MAXIZ+1                                                 
          XINTS(IY,IZ) = FACTOR * XINTS(IY,IZ)                                  
  100   CONTINUE                                                                
C       WRITE (6,'('' PROJEC: IY,YL,YU,TL,TU,FACT,YOUT'',I5,6F10.5)')           
C    >    IY,YL,YU,TL,TU,FACTOR,YOUTS(IY)                                       
  110 CONTINUE                                                                  
C                                                                               
      DO 120 IZ = -2, MAXIZ+1                                                   
        XINTS(0,IZ) = 0.5 * (XINTS(-1,IZ) + XINTS(1,IZ))                        
  120 CONTINUE                                                                  
C                                                                               
      DO 210 IY = 1, MAXNYS                                                     
        TEMP       = YOUTS(IY)                                                  
        YOUTS(IY)  = YOUTS(-IY)                                                 
        YOUTS(-IY) = TEMP                                                       
        TEMP       = YWIDSS(IY)                                                 
        YWIDSS(IY) = YWIDSS(-IY)                                                
        YWIDSS(-IY)= TEMP                                                       
        DO 200 IZ = -2, MAXIZ+1                                                 
          TEMP          = XINTS(IY,IZ)                                          
          XINTS(IY,IZ)  = XINTS(-IY,IZ)                                         
          XINTS(-IY,IZ) = TEMP                                                  
  200   CONTINUE                                                                
  210 CONTINUE                                                                  
C                                                                               
C     WRITE (6,'('' PROJEC:  YOUTS ='',/,(1X,8F9.5))')                          
C    >      (YOUTS(IY),IY=-MAXNYS,MAXNYS)                                       
C     WRITE (6,'('' PROJEC: YWIDSS ='',/,(1X,8F9.5))')                          
C    >      (YWIDSS(IY),IY=-MAXNYS,MAXNYS)                                      
      RETURN                                                                    
      END                                                                       
C
      SUBROUTINE RINTR (VS,IPLANE,MAXIZ,RINTS,IFOLD,            
     >                  SSS,FP,FT,AVER,SMIN,SMAX)                             
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C
      INCLUDE 'rtheta'
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    RINTS(MAXNAS,-2:MAXIZS+1),FP,FT,AVER
      REAL    SMIN,SMAX
      INTEGER MAXIZ,IFOLD,IPLANE   
      LOGICAL SSS                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTR:  TAKE REAL ARRAY VS, MAP A SPECIFIED R,THETA GRID ONTO IT *
C  *  AND INTEGRATE OVER THE R DIMENSION. NO 3D IS POSSIBLE RIGHT NOW. *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO RINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTRS.          *        
C  *                                                                   *        
C  *                                                                   *        
C  *  D. ELDER , JUNE 19, 1990                                         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTRS(MAXNAS,-2:MAXIZS+1),DRSUM1          
      INTEGER IZ,IR,IT,IPOS,SMIN1,SMAX1                                     
      EXTERNAL IPOS
C                                                                               
      CALL DZERO (DNTRS, (MAXNAS)*(MAXIZS+4))                               
C                                                                               
C     SET UP LIMITS OF INTEGRATION
C
      SMIN1 = IPOS(SMIN,RS,NRS-1)
      SMAX1 = IPOS(SMAX,RS,NRS-1)
C
        DO 130 IZ = -1, MAXIZ                                                   
         CALL RTXY(VS,NRS,NAS,IFOLD,IZ,MAXIZ)
         DO 120 IT = 1, NAS                                                     
          DRSUM1 = 0.0D0                                                        
          DO 110 IR = SMIN1, SMAX1                                              
              DRSUM1 = DRSUM1 + 
     >          DBLE (RTMESH(IR,IT) * RWIDS(IR))                 
  110     CONTINUE                                                              
          DNTRS(IT,IZ) = DRSUM1                             
  120    CONTINUE                                                               
  130   CONTINUE                                                                
C                        
      DO 330 IT = 1,MAXNAS                                                     
        DNTRS(IT,-2) = FT * DNTRS(IT,0) - FP * DNTRS(IT,-1)                     
        DO 310 IZ = 1, MAXIZ                                                    
          DNTRS(IT,MAXIZ+1) = DNTRS(IT,MAXIZ+1) + DNTRS(IT,IZ)                  
  310   CONTINUE                                                                
        DO 320 IZ = -2, MAXIZ+1                                                 
          RINTS(IT,IZ) = SNGL (DNTRS(IT,IZ)) / AVER                             
  320   CONTINUE                                                                
  330 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C                                                                               
C
      SUBROUTINE RINTT (VS,IPLANE,MAXIZ,TINTS,IFOLD,            
     >                  SSS,FP,FT,AVER,SMIN,SMAX)                             
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)                                                           
C                           
      INCLUDE 'rtheta'
C                                                    
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    TINTS(MAXNRS,-2:MAXIZS+1),FP,FT,AVER             
      REAL    SMIN,SMAX
      INTEGER MAXIZ,IFOLD,IPLANE
      LOGICAL SSS
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTT:  TAKE REAL ARRAY VS, MAP A SPECIFIED R,THETA GRID ONTO IT *
C  *  AND INTEGRATE OVER THE T DIMENSION. NO 3D IS POSSIBLE RIGHT NOW. *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO RINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTRS.          *        
C  *                                                                   *        
C  *                                                                   *        
C  *  D. ELDER , JUNE 19, 1990                                         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTTS(MAXNRS,-2:MAXIZS+1),DTSUM1          
      INTEGER IZ,IR,IT,SMIN1,SMAX1,SMIN2,SMAX2,IPOS
      EXTERNAL IPOS
C                                                                               
      CALL DZERO (DNTTS, (MAXNRS)*(MAXIZS+4))                               
C                                                                               
C     SET UP BOUNDARIES OF INTEGRATION. SINCE THE INTEGRATION RANGE CAN BE 
C     SPECIFIED IN TWO WAYS IT IS NECESSARY TO BE ABLE TO SPLIT THE 
C     INTEGRATION. IN SOME CASES, ONLY THE FIRST AND LAST PORTIONS OF 
C     THE ARRAY MUST BE INTEGRATED.
C
      SMIN1 = 0
      SMAX1 = 0
      SMIN2 = 0
      SMAX2 = 0
      WRITE(6,*) 'SMIN,2*PI+SMIN:',SMIN,2*PI+SMIN
      WRITE(6,*) 'SMAX,2*PI+SMAX:',SMAX,2*PI+SMAX

      IF (SMIN.LT.0.0) THEN 
         IF (SMAX.LT.0.0) THEN    
            SMIN1 = IPOS(2*PI+SMIN,TS,NAS-1)
            SMAX1 = IPOS(2*PI+SMAX,TS,NAS-1)
         ELSE
            SMIN1 = 1
            SMAX1 = IPOS(SMAX,TS,NAS-1)
            SMIN2 = IPOS(2*PI+SMIN,TS,NAS-1)
            SMAX2 = NAS
         ENDIF
      ELSE
         SMIN1 = IPOS(SMIN,TS,NAS-1)
         SMAX1 = IPOS(SMAX,TS,NAS-1)
         WRITE(6,*) 'SMIN,SMIN1,SMAX,SMAX1:',SMIN,SMIN1,
     >             SMAX,SMAX1 
      ENDIF 
C
C
        WRITE(6,*) 'SMIN1,SMAX1,SMIN2,SMAX2:',SMIN1,SMAX1,
     >              SMIN2,SMAX2
        DO 130 IZ = -1, MAXIZ                                                   
         CALL RTXY(VS,NRS,NAS,IFOLD,IZ,MAXIZ)
         DO 120 IR = 1, NRS                                                     
          DTSUM1 = 0.0D0                                                        
          DO 110 IT = SMIN1, SMAX1                                              
              DTSUM1 = DTSUM1 + 
     >          DBLE (RTMESH(IR,IT) * TWIDS(IT)* ROUTS(IR))                 
  110     CONTINUE                                                              
          IF (SMIN2.NE.0) THEN 
            DO 105 IT = SMIN2, SMAX2                                        
                DTSUM1 = DTSUM1 + 
     >            DBLE (RTMESH(IR,IT) * TWIDS(IT)* ROUTS(IR))                 
  105       CONTINUE                                                         
          ENDIF  
          DNTTS(IR,IZ) = DTSUM1                             
  120    CONTINUE                                                               
  130   CONTINUE                                                                
C                                                                               
      DO 330 IR = 1,MAXNRS                                                     
        DNTTS(IR,-2) = FT * DNTTS(IR,0) - FP * DNTTS(IR,-1)                     
        DO 310 IZ = 1, MAXIZ                                                    
          DNTTS(IR,MAXIZ+1) = DNTTS(IR,MAXIZ+1) + DNTTS(IR,IZ)                  
  310   CONTINUE                                                                
        DO 320 IZ = -2, MAXIZ+1                                                 
          TINTS(IR,IZ) = SNGL (DNTTS(IR,IZ)) / AVER                             
  320   CONTINUE                                                                
  330 CONTINUE                                                                  
C                                                                               
      RETURN                                                                    
      END                                                                       
C
C
      SUBROUTINE RINTPSI (VS,IPLANE,MAXIZ,SINTS,IFOLD,            
     >                  SSS,FP,FT,AVER,VMIN,VMAX,SMIN,SMAX)             
      IMPLICIT none
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
      INCLUDE 'comtor'                                                          
C     INCLUDE (COMTOR)                                                          
      INCLUDE 'comvu'                                                           
C     INCLUDE (COMVU)     
C                                                      
      INCLUDE 'rtheta'
C                                                                               
      REAL    VS(MAXNXS,-MAXNYS:MAXNYS,-1:MAXIZS)                               
      REAL    SINTS(MAXNSS,-2:MAXIZS+1),FP,FT,AVER,VMIN,VMAX
      REAL    SMIN,SMAX   
      INTEGER MAXIZ,IFOLD,IPLANE   
      LOGICAL SSS                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  RINTPSI:TAKE REAL ARRAY VS, MAP A SPECIFIED R,THETA GRID ONTO IT *
C  *  AND THEN INTEGRATE FROM A GIVEN OBSERVATION POSITION.            *
C  *  THE SPACE IS DIVIDED INTO EQUAL ARCS FROM THE OBSERVATION POINT  *
C  *  AND ALL BIN CENTRES LYING IN AN ARC CONTRIBUTE TO THE INTEGRATED *
C  *  FLUX IN THAT ARC.                                                *        
C  *  IFOLD PARAMETER WHEN SET TO 1 INDICATES THAT -Y AND +Y RESULTS   *        
C  *  ARE TO BE ADDED TOGETHER.  RESULTS ARE TOTALLED FOR ALL CHARGE   *        
C  *  STATES AND WRITTEN TO RINTS(,MAXIZ+1) LOCATIONS. ALL SUMMING IS  *        
C  *  DONE IN D.P. FOR ACCURACY, USING TEMPORARY ARRAY DNTSS.          *        
C  *                                                                   *        
C  *                                                                   *        
C  *  D. ELDER , JUNE 19, 1990                                         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      DOUBLE PRECISION DNTSS(MAXNSS,-2:MAXIZS+1)
      INTEGER IZ,IR,IN,IPOS
      EXTERNAL IPOS
      REAL RACT,PRWID2,ATMP,ROBS,AACT
      REAL RDISTS(MAXNRDIV,2)
      REAL VAL,XP,YP,SCOS,SSIN    
C
      REAL DX,DY,AT,A11,A12,A21,A22
      REAL VAL11,VAL12,VAL21,VAL22
      INTEGER IX1,IY1,IX2,IY2
C
C                                                                               
      CALL DZERO (DNTSS, (MAXNSS)*(MAXIZS+4))                               
C                                                          
C     CALCULATE THE X,Y COORDINATES OF THE R,THETA BIN CENTRES IN THE 
C     PSI COORDINATE SYSTEM CENTRED ON THE OBSERVATION POSITION WITH 
C     THE POSITIVE X-AXIS CONNECTING THE OBSERVATION POSITION TO THE 
C     PLASMA CENTRE
C
C     SET UP RANGE OF PSI VALUES
C
      WRITE(6,*) 'RPSI:1'
      IF (.NOT.RTDONE) THEN
        WRITE(6,*) '1A:'
        IF (RVPOS.GE.(CA-CAW)) THEN
           WRITE(6,*) 'CA-CAW,RVPOS,RATIO:',CA-CAW,RVPOS,
     >          (CA-CAW)/RVPOS
           WRITE(6,*) 'NSS:',NSS 
           PSIMAX = ASIN((CA-CAW)/RVPOS)
           PSIMIN = -PSIMAX
           PSISTEP = (PSIMAX-PSIMIN)/NSS
           WRITE(6,*) '1B>:',RVPOS 
        ELSE
           PSIMAX = PI
           PSIMIN = -PI
           PSISTEP = 2*PI / NSS
           WRITE(6,*) '1B<:',RVPOS
         ENDIF
        PSIHALF = 0.5 * PSISTEP  
C
C       SET UP PSI COORDINATES OF ARC CENTRES
C
        WRITE(6,*) '1C,PSIMIN,PSIMAX:',PSIMIN,PSIMAX
        DO 20 IN = 1,NSS
C
C          CALCULATE SOUTS IN RADIANS AND SOUTS1 IN DEGREES FOR 
C          PLOTTING PURPOSES 
C
C          D. ELDER , AUG 3 1990
C
           SOUTS(IN) = PSIMIN + (IN-1) * PSISTEP + PSIHALF
           SOUTS1(IN) = SOUTS(IN) * RADDEG
           SWIDS(IN) = PSISTEP
20      CONTINUE
C
C       SET UP COORDINATES OF THE RPSI,THETAPSI BINS IN X,Y SYSTEM 
C       BASED ON THE VIEWING POSITION AND TORUS CHARACTERISTICS
C
        WRITE(6,*) 'RPSI:2.5 : ',NRDIV
        PSIRWID = (RVPOS+CA-CAW)/NRDIV
        PRWID2 = PSIRWID / 2.0
        ROBS = RVPOS**2
C        WRITE(6,*) 'PSIRWID,PRWID2,ROBS: ',PSIRWID,PRWID2,ROBS    
C        WRITE(6,*) 'RVPOS,TVPOS : ', RVPOS, TVPOS
        DO 45 IR = 1 , NRDIV
          RDISTS(IR,1) = REAL(IR) * PSIRWID+PRWID2
          RDISTS(IR,2) = RDISTS(IR,1)**2 
45      CONTINUE
        DO 50 IN = 1,NSS
           SCOS = COS(SOUTS(IN))
           SSIN = SIN(SOUTS(IN))
C           WRITE(6,*) 'IN,SIN,COS,S:',IN,SSIN,SCOS,SOUTS(IN)  
           DO 60 IR = 1,NRDIV
C
C          CALCULATE THE ACTUAL R,THETA POSITIONS CORRESPONDING
C          TO THE NEXT PSI BIN TO BE INDEXED. USE THESE TO FIND
C          THE INDICES OF THE CORRESPONDING X,Y BIN IN WHICH THE
C          POINT FALLS.
C
C          D.ELDER AUG 28, 1990 
C
             RACT = SQRT(RDISTS(IR,2)+ROBS-2*RDISTS(IR,1)*RVPOS*SCOS)
             ATMP = ACOS((RACT**2+ROBS-RDISTS(IR,2))/(2.0*RACT*RVPOS))
             AACT = TVPOS + SIGN(ATMP,SOUTS(IN))  
C             WRITE(6,*) 'AACT1: ', AACT
C
C            ADJUST THE THETA RANGE TO LIE IN [0,2PI]
C
             IF (AACT.GT.2.0*PI) THEN
                AACT = AACT - 2.0*PI
             ELSEIF (AACT.LT.0.0) THEN
                AACT = AACT + 2.0*PI
             ENDIF      
C             WRITE(6,*) 'RD,SSIN,RACT,RD*SS/RA,ASIN(),AACT:',
C     >        RDISTS(IR,1),SSIN,RACT,RDISTS(IR,1)*SSIN/RACT,
C     >        ATMP ,AACT
             XP = CA - RACT
             YP = CL * AACT / PI
C             WRITE(6,*) 'XP,YP:',XP,YP
             IF (XP.LT.CAW) THEN 
                RTSXY(IR,IN,1) = -1
                RTSXY(IR,IN,2) =  0
             ELSE     
                RTSXY(IR,IN,1) = IPOS(XP,XS,MAXNXS-1)
                RTSXY(IR,IN,2) = IPOS(YP,YS,2*MAXNYS-1)
                RTAXY(IR,IN,1) = XP
                RTAXY(IR,IN,2) = YP
             ENDIF
C             WRITE(6,*) 'IR,IN,RTSXY 1, 2: ', IN, IR,
C     >            RTSXY(IR,IN,1),RTSXY(IR,IN,2)     
60         CONTINUE
50      CONTINUE  
      RTDONE = .TRUE.
      ENDIF
      WRITE(6,*) 'RPSI:2: PSIMIN,PSIMAX,PSISTEP',PSIMIN,PSIMAX,
     >               PSISTEP    
C
C      ADJUST THE PLOTTING RANGES VMIN,VMAX AS REQUIRED TO 
C      CONFORM WITH THE CALCULATED LIMITS PSIMIN,PSIMAX 
C      
C      THEN CONVERT THE VALUES TO DEGREES FOR USE IN THE PLOTTING
C      ROUTINE
C
       VMIN = MAX(VMIN,PSIMIN) * RADDEG
       VMAX = MIN(VMAX,PSIMAX) * RADDEG
C
C      LOOP THROUGH THE RPSI, THETAPSI GRID DETERMINING THE 
C      DENSITY IN EACH REGION AND THEN ADD THE FLUX 
C      CONTRIBUTION TO SUM FOR THAT VIEW ANGLE.
C       
      WRITE(6,*) 'VMIN,VMAX: ',VMIN,VMAX
      WRITE(6,*) 'SMIN,SMAX: ',SMIN,SMAX
C
C     ADJUST THE X-INTEGRATION RANGE FOR THE PSI PLOTS. THIS ALLOWS
C     SPECIFIC CENTRAL OR EDGE CONTRIBUTIONS TO BE REMOVED FORM THE 
C     PSI PLOT INTEGRATIONS. OVER ESTIMATES ARE ALLOWED AND A RANGE
C     OF 0.0 TO 0.0 IS INTERPRETED AS -99.0 TO 99.0 TO SAVE SPACE IN 
C     THE INPUT FILE.
C
      IF (SMIN.EQ.0.0.AND.SMAX.EQ.0.0) THEN 
         SMIN = -99.0
         SMAX =  99.0
      ENDIF
C
C
      DO 70 IZ = -1,MAXIZ
        DO 80 IN = 1,NSS
          DO 90 IR = 1,NRDIV
C
C           IF IT IS NOT OUTSIDE THE PLASMA AND IS WITHIN THE 
C           X-INTEGRATION RANGE THEN ADD TO THE BIN CONTENTS
C
            IF (RTSXY(IR,IN,1).NE.-1.AND.
     >          RTAXY(IR,IN,1).GE.SMIN.AND.
     >          RTAXY(IR,IN,1).LE.SMAX) THEN
C
C            CALCULATE THE OTHER QUANTITIES NECESSARY FOR
C            4 BIN AVERAGING. RAW PSI PLOTS APPEAR TO BE 
C            EXCESSIVELY NOISY AND STATISTICS DEPENDENT.
C            IT IS HOPED THAT BY FORMING THE POINT DENSITY VALUES 
C            FROM AN AREA WEIGTHED MEAN OF THE FOUR NEAREST BINS
C            A SMOOTHER AND MORE REPRESENTATIVE PLOT MAY BE 
C            FORMED.
C
              IX1= RTSXY(IR,IN,1)
              IY1 =RTSXY(IR,IN,2)
              XP = RTAXY(IR,IN,1)
              YP = RTAXY(IR,IN,2)          
              DX = XP - XOUTS(IX1)
              DY = YP - YOUTS(IY1)
              IF (DX.LE.0.0D0) THEN
                IX2 = IX1-1
              ELSE
                IX2= IX1+1
              ENDIF
              IF (DY.LE.0.0D0) THEN
                IY2= IY1-1
              ELSE
                IY2 = IY1+1
              ENDIF
              IF ((DX.EQ.0.0) .OR.
     >            (DY.EQ.0.0) .OR.
     >            (IX2.EQ.0)  .OR.
     >            (IX2.EQ.NXS+1) .OR.
     >            (IY2.EQ.0)  .OR.
     >            (IY2.EQ.NYS+1)) THEN
                 A11 = 0.0
                 A12 = 0.0
                 A21 = 0.0
                 A22 = 1.0
                 AT  = 1.0
               ELSE
                 AT = ABS((XOUTS(IX1)-XOUTS(IX2))*
     >                      (YOUTS(IY1)-YOUTS(IY2)))        
                 A11 = ABS(DX*DY)
                 A12 = ABS(DX*(YP-YOUTS(IY2)))        
                 A21 = ABS((XP-XOUTS(IX2))*DY)
                 A22 = ABS((XP-XOUTS(IX2))*(YP-YOUTS(IY2)))        
               ENDIF 
C
               VAL11 = VS(IX1,IY1,IZ)
               VAL12 = VS(IX1,IY2,IZ)
               VAL21 = VS(IX2,IY1,IZ)
               VAL22 = VS(IX2,IY2,IZ)
               IF (IFOLD.EQ.1) THEN 
                 VAL11= 0.5*(VAL11+
     >             VS(IX1,-IY1,IZ))
                 VAL12 = 0.5*(VAL12+
     >             VS(IX1,-IY2,IZ))
                 VAL21 = 0.5*(VAL21+
     >             VS(IX2,-IY1,IZ))
                 VAL22 = 0.5*(VAL22+
     >             VS(IX2,-IY2,IZ))
               ENDIF
               VAL = (A11 * VAL22  +
     >                A12 * VAL21  +
     >                A21 * VAL12  +
     >                A22 * VAL11) / AT
               DNTSS(IN,IZ) = DNTSS(IN,IZ) + 
     >                DBLE(  VAL*PSISTEP*PSIRWID/(2.0*PI))
            ENDIF 
90        CONTINUE          
80      CONTINUE
70    CONTINUE

      WRITE(6,*) 'RPSI:4'  
 
C                        
      DO 330 IN = 1,NSS                                                     
        DNTSS(IN,-2) = FT * DNTSS(IN,0) - FP * DNTSS(IN,-1)                     
        DO 310 IZ = 1, MAXIZ                                                    
          DNTSS(IN,MAXIZ+1) = DNTSS(IN,MAXIZ+1) + DNTSS(IN,IZ)              
  310   CONTINUE                                                                
        DO 320 IZ = -2, MAXIZ+1                                                 
          SINTS(IN,IZ) = SNGL (DNTSS(IN,IZ)) / AVER                             
  320   CONTINUE                                                                
  330 CONTINUE                                                         
C                                                                               
      RETURN                                                                    
      END                                                                    
