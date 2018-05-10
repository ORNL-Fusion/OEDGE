      SUBROUTINE EDGE (QXS,QEDGES,QTANS,QDISTS,NQXSO,CAW,CL,ICUT,CCUT,          
     >     XSCALO,WEDGAN,XL1,YL1,XL2,YL2,TC,SC,TO,SO,GC,RP,CIOPTH,
     >     CORECT,
     >     XST1,YST1,XST2,YST2,XST3,YST3,RLEDGE7,CA,RLC)
      use iter_bm
      use mod_comnet
      IMPLICIT  none
      INCLUDE   'params'                                                        
C     INCLUDE   (PARAMS)                                                        
c      INCLUDE   'comnet'                                                        
C     INCLUDE   (COMNET)                                                        
c
      include   'global_options'

      INTEGER   NQXSO,ICUT(2),CIOPTH,CORECT                                     
      REAL      QXS(-MAXQXS:MAXQXS),QEDGES(-MAXQXS:0,2),CAW,CL,CCUT             
      REAL      QTANS(-MAXQXS:0,2),XL1(2),YL1(2),XL2(2),YL2(2),XSCALO           
      REAL      WEDGAN(2),TC,SC,TO,SO,GC,RP,QDISTS(-MAXQXS:0,2)                 
      REAL      XST1(2),YST1(2),XST2(2),YST2(2),XST3(2),YST3(2)
      REAL      RLEDGE7,CA,RLC
C     
C     *********************************************************************        
C     *                                                                   *        
C     *  EDGE:  THIS ROUTINE CALCULATES THE LIMITER SHAPE AND STORES IT   *        
C     *  AS A SET OF Y VALUES CORRESPONDING TO THE X POSITIONS IN "QXS".  *        
C     *  THE SHAPE CAN BE DIFFERENT ON EACH SIDE OF Y = 0.  THE EDGES     *        
C     *  FOR THE Y < 0 REGION ARE STORED IN QEDGES(,1) AND THE SURFACE    *        
C     *  NORMALS (WITH RESPECT TO THE X AXIS) ARE STORED IN QTANS(,1).    *        
C     *  THE CORRESPONDING VALUES FOR THE Y > 0 REGION ARE STORED IN      *        
C     *  THE (,2) POSITIONS.  NOTE THAT FOR THE Y < 0 REGION, THE STORED  *        
C     *  VALUES ARE CALCULATED AS IF THE LIMITER EDGE WERE FLIPPED OVER   *        
C     *  INTO THE Y > 0 REGION, SO THAT WHEN ACTUALLY LAUNCHING A NEUTRAL *        
C     *  WE WILL USE -QEDGES AND PI-QTANS FOR THE Y < 0 REGION.           *        
C     *  PARAMETER ICUT IS SET TO THE MINIMUM X POINT FROM WHICH WE       *        
C     *  WILL LAUNCH.  IT CAN BE OVER-RIDDEN BY THE INPUT VALUE "CCUT".   *        
C     *  THE LIMITER EDGE IS "CORRECTED" BY A SIMPLE ALGORITHM IF         *        
C     *  REQUIRED, TO ACCOUNT FOR THE PLASMA CURVATURE.                   *        
C     *  DISTANCES ALONG LIMITER SURFACE ARE STORED IN QDISTS.            *        
C     *                                                                   *        
C     *  CHRIS FARRELL   (HUNTERSKIL)  SEPT 1988                          *        
C     *  D. ELDER    MAY 1990 - ADDED BELT LIMITER FROM J. SPENCE (JET)   *
C     *                                                                   *
C     *                                                                   *        
C     *********************************************************************        
C     
      INTEGER IQX,J,NJ,NK,INTVL,K,JC,IO,IN,JPOS
      EXTERNAL JPOS
      PARAMETER (NJ=246,NK=15)                                                  
      REAL    C(0:NK+1),TK(0:NK+1),T(NJ),S(NJ),X(NJ),Y(NJ),DS(NJ),HH,H          
c     REAL    DEGRAD,PI,X0,Y0,RADIUS,XMAX,RADS2,YMAX,COSGC,SINGC,AA,B           
      REAL    X0,Y0,RADIUS,XMAX,RADS2,YMAX,COSGC,SINGC,AA,B           
      REAL    TEMP(-MAXQXS:1),DR(-MAXQXS:1),A,DA,GAM,THETA,FACT                 
      REAL    X08,Y08
C     
C     VARIABLES FOR EDGE OPTION 7 
C     
      REAL NUM1,NUM2,NUM3,NUM4,NUM5,NUM6
      REAL ALPHA,BETA,VAL,INTVAL0,INTVAL1,INTVAL2,INTVAL3
      REAL SLOPE
C     
C     Edge option 9
C     
      INTEGER   MAXNP,NP
      PARAMETER (MAXNP=18,NP = 18)   
      REAL      XP(MAXNP),YP(MAXNP) 
C     
      DATA XP /0.0, -0.001, -0.002, -0.003, -0.004, -0.005, -0.006,
     >     -0.007, -0.008, -0.009, -0.010, -0.011, -0.012, -0.013,
     >     -0.014, -0.015, -0.016, -0.017/ 
      DATA YP/ 0.0893, 0.0956, 0.1012, 0.1063, 0.1109, 0.1152, 0.1190,
     >     0.1224, 0.1255, 0.1283, 0.1309, 0.1332, 0.1353, 0.1371,
     >     0.1386, 0.1401, 0.1413, 0.1422/

C     
c     DATA    DEGRAD /.017453292/,       PI /3.141592654/                       
      DATA    X0 /-0.065/,  Y0 /0.0/,    RADIUS /0.065/,  XMAX /-0.1/           
      DATA    C /119.4256, -20.5781, -36.6792, -37.4321, -32.6931,              
     >     -22.9295, -11.1687, 1.625282, 15.10991, 28.54260,              
     >     43.68845, 58.98647, 74.96905, 89.81274, 110.4915,              
     >     126.4284, 240.4677/                                            
c     
c     Edge Option 10
c     
      real c1,c2,xtmp,ytmp
      real theta1,theta2
      real atan2c
      external atan2c


C     
C     *********************************************************************        
C     *  EDGE0:  SLAB LIMITER HAS 0 THICKNESS                             *        
C     *********************************************************************        
C     
      IF (CIOPTH.EQ.0) THEN                                                     
         DO J = 1, 2                                                          
            ICUT(J) = 1-NQXSO                                                     
            DO IQX = 1-NQXSO, 0                                                
               QEDGES(IQX,J) = 0.0                                                 
               QTANS(IQX,J)  = 0.0                                                 
            end do
         end do
C     
C     *********************************************************************        
C     *  EDGE1: THE LIMITER IS WEDGE SHAPED WITH AN ANGLE THETA IN RADIANS*        
C     *  THE TWO ANGLES GIVEN DEFINE THE WEDGE SHAPE ON EACH SIDE OF Y=0. *        
C     *                                                                   *        
C     *  :.                                                               *        
C     *  : :.                                                             *        
C     *  :   :.                                                           *        
C     *  :     :.                                                         *        
C     *  :.....(.:. <-- WEDGAN MULTIPLIED BY DEGREES-->RADIANS FACTOR     *        
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.1) THEN                                                 
         DO J = 1, 2                                                         
            ICUT(J) = 1-NQXSO                                                     
            DO IQX = 1-NQXSO, 0                                               
               QEDGES(IQX,J) = -QXS(IQX) / TAN (DEGRAD*WEDGAN(J))                  
               QTANS(IQX,J)  = PI/2.0 - DEGRAD*WEDGAN(J)                           
            end do
         end do 
C     
C     *********************************************************************        
C     *  EDGE2: THE LIMITER IS SHAPED LIKE THE ARC OF A CIRCLE WITH CENTRE*        
C     *  (X0,Y0) AND RADIUS = RADIUS.                                     *        
C     *  THE ARC IS STOPPED AT X=XMAX AND Y IS HELD CONST FOR |X|<|XMAX|  *        
C     *       0....                                                       *        
C     *        :   :..                                                    *        
C     *        :      :.                                                  *        
C     *    XMAX:       :.                                                 *        
C     *        :        :                                                 *        
C     *        :        :                                                 *        
C     *        :        :                                                 *        
C     *      AW:........:                                                 *        
C     *        0       YMAX                                               *        
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.2) THEN                                                 
         RADS2 = RADIUS * RADIUS                                                 
         YMAX  = SQRT (RADS2-(XMAX-X0)**2) + Y0                                  
         DO  J = 1, 2                                                         
            ICUT(J) = 1-NQXSO                                                     
            DO  IQX = 1-NQXSO, 0                                               
               IF (ABS(QXS(IQX)).LE.ABS(XMAX)) THEN                                
                  QEDGES(IQX,J) = SQRT (RADS2-(QXS(IQX)-X0)**2) + Y0                
                  QTANS(IQX,J) = ATAN ((QXS(IQX)-X0)/(QEDGES(IQX,J)-Y0))          
               ELSE                                                                
                  QEDGES(IQX,J) = YMAX                                              
                  QTANS(IQX,J)  = 0.0                                               
                  ICUT(J) = ICUT(J) + 1                                             
               ENDIF                                                               
            end do 
         end do
C     
C     *********************************************************************        
C     *  EDGE3: THE LIMITER IS A "BLUNT NOSE" SPECIFIED BY TWO POINTS     *        
C     *  THE TWO POINTS ARE DEFINED SEPARATELY FOR EACH SIDE OF Y=0.      *        
C     *                                                                   *        
C     *  0........(XL1,YL1)                                               *        
C     *   :       :.                                                      *        
C     *   :         :.                                                    *        
C     *   :           :.(XL2,YL2)                                         *        
C     *   :             :                                                 *        
C     *   :             :                                                 *        
C     *   :             :                                                 *        
C     *   :.............:                                                 *        
C     *   0                                                               *        
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.3) THEN                                                 
         DO J = 1, 2                                                         
            ICUT(J) = 1-NQXSO                                                     
            DO IQX = 1-NQXSO, 0                                               
               IF (QXS(IQX).LT.XL2(J)) THEN                                        
                  QEDGES(IQX,J) = YL2(J)                                            
                  QTANS(IQX,J)  = 0.0                                               
               ELSEIF (QXS(IQX).LT.XL1(J)) THEN                                    
                  QEDGES(IQX,J) = YL1(J) + (XL1(J)-QXS(IQX)) *                      
     >                 (YL2(J)-YL1(J)) / (XL1(J)-XL2(J))                 
                  QTANS(IQX,J)  = PI/2.0 -                                          
     >                 ATAN ((XL1(J)-XL2(J)) / (YL2(J)-YL1(J)))          
               ELSEIF (QXS(IQX).GE.XL1(J)) THEN                                    
                  QEDGES(IQX,J) = QXS(IQX) * YL1(J) / XL1(J)                        
                  QTANS(IQX,J)  = PI/2.0 - ATAN (-XL1(J)/YL1(J))                    
               ENDIF                                                               
            end do 
         end do

C     
C     *********************************************************************        
C     *  EDGE4:  THE LIMITER IS A "BELT LIMITER", DESCRIBED IN NOTE 196.  *        
C     *  BEWARE CONVERSIONS FROM MM TO M UNITS.                           *        
C     *  IT IS NOT SYMMETRICAL ABOUT Y = 0                                *        
C     *                                                                   *        
C     *             .........                                             *        
C     *          ..:     :   :.....                                       *        
C     *       ..:        :         :...                                   *        
C     *     .:           :             :..                                *        
C     *     :            :                :                               *        
C     *     :            :               :                                *        
C     *     :            :              :                                 *        
C     *     :............:..............:                                 *        
C     *                  0                                                *        
C     *  EDGE5:  SYMMETRICAL "BELT LIMITER", USING Y>0 SHAPE REFLECTED    *        
C     *  ABOUT X AXIS.                                                    *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.4.OR.CIOPTH.EQ.5) THEN                                  
C     
C------CALCULATE POINT OF INFLECTION OF CURVED LIMITER SURFACE                 
C------INITIALISE T POINTS, ETC                                                
C     
         T(1) = .000186                                                          
         T(NJ)= 261.0627                                                         
         HH   = (T(NJ)-T(1)) / REAL(NJ-1)                                        
         DO J = 2, NJ-1                                                      
            T(J) = T(1) + HH * REAL(J-1)                                          
            IF (TC.GE.T(J)*1.E-3) JC = J                                          
         end do
         COSGC = COS (GC)                                                        
         SINGC = SIN (GC)                                                        
         ICUT(1) = 1-NQXSO                                                       
         ICUT(2) = 1-NQXSO                                                       
C     
C------SET KNOT POINTS FOR PIECEWISE CUBIC SPLINE APPROXIMATION                
C     
         H = (T(NJ)-T(1)) * 1.00002 / REAL(NK-1)                                 
         AA= T(1)  - (T(NJ)-T(1)) * .00001                                       
         B = T(NJ) + (T(NJ)-T(1)) * .00001                                       
         DO K = 1, NK                                                        
            TK(K) = H * REAL(K-1) + AA                                            
         end do 
C     
C------CALCULATE BETA FUNCTIONS DEFINING SPLINE, AND THEREFORE SPLINE          
C------VALUES S AT EACH POINT T.  CONVERT TO X,Y COORDINATE SYSTEM             
C     
         DO 430 J = 1, NJ                                                        
            INTVL = INT((T(J)-AA)/H) + 1                                          
            S(J)  = 0.0                                                           
            DS(J) = 0.0                                                           
            DO 420 K = 0, NK+1                                                    
               A = 0.0                                                             
               DA= 0.0                                                             
               IF     (K.EQ.INTVL+2) THEN                                          
                  GAM = (T(J)-TK(K-2)) / H                                          
                  A = GAM*GAM*GAM / 4.0                                             
                  DA= 3.0*GAM*GAM / (4.0*H)                                         
               ELSEIF (K.EQ.INTVL+1) THEN                                          
                  GAM = (T(J)-TK(K-1)) / H                                          
                  A = (1.0+3.0*GAM*(1.0+GAM-GAM*GAM)) /4.0                          
                  DA= (3.0+6.0*GAM-9.0*GAM*GAM) / (4.0*H)                           
               ELSEIF (K.EQ.INTVL  ) THEN                                          
                  GAM = (TK(K+1)-T(J)) / H                                          
                  A = (1.0+3.0*GAM*(1.0+GAM-GAM*GAM)) /4.0                          
                  DA=-(3.0+6.0*GAM-9.0*GAM*GAM) / (4.0*H)                           
               ELSEIF (K.EQ.INTVL-1) THEN                                          
                  GAM = (TK(K+2)-T(J)) / H                                          
                  A = GAM*GAM*GAM/4.0                                               
                  DA=-3.0*GAM*GAM / (4.0*H)                                         
               ENDIF                                                               
               S(J) = S(J) + A * C(K)                                              
               DS(J)= DS(J)+ DA* C(K)                                              
 420        CONTINUE                                                              
            Y(J) = (TC-T(J)*1.E-3)*COSGC + (SC-S(J)*1.E-3)*SINGC                  
            X(J) = (SC-S(J)*1.E-3)*COSGC - (TC-T(J)*1.E-3)*SINGC                  
C     WRITE (6,9002) J,INTVL,X(J),Y(J),                                     
C     >      S(J),T(J),DS(J)                                                     
 430     CONTINUE                                                                
C     
C------CALCULATE EDGES AND NORMALS FOR Y < 0 REGION                            
C------QEDGES(,1) AND QTANS(,1) ARE SET UP AS IF THE LIMITER WERE              
C------REALLY IN THE Y > 0 REGION.                                             
C     
         DO 450 IQX = 1-NQXSO, 0                                                 
            IF (QXS(IQX).LE.X(NJ)) THEN                                           
               QEDGES(IQX,1) = -Y(NJ)                                              
               QTANS(IQX,1)  = 0.0                                                 
            ELSE                                                                  
               J = NJ                                                              
 440           CONTINUE                                                            
               IF (QXS(IQX).GT.X(J) .AND. J.GT.1) THEN                             
                  J = J - 1                                                         
                  GOTO 440                                                          
               ENDIF                                                               
               QEDGES(IQX,1) =-Y(J) + (QXS(IQX)-X(J)) * (Y(J)-Y(J+1))/            
     >              ( X(J+1) -X(J))                              
               QEDGES(IQX,1) = MAX (QEDGES(IQX,1),0.0)                             
               THETA = ATAN (DS(J)- (QXS(IQX)-X(J)) * (DS(J)-DS(J+1))/            
     >              ( X(J+1) -X(J)))                               
               QTANS(IQX,1)  = PI/2.0 + GC - THETA                                 
C     WRITE (6,9001) IQX,QXS(IQX),-QEDGES(IQX,1),THETA/DEGRAD,            
C     >        (PI-QTANS(IQX,1))/DEGRAD,J,T(J),S(J),X(J),Y(J)                    
            ENDIF                                                                 
 450     CONTINUE                                                                
C     
C------CALCULATE EDGES AND NORMALS FOR Y > 0 REGION                            
C     
         DO 470 IQX = 1-NQXSO, 0                                                 
            IF (QXS(IQX).LE.X(1)) THEN                                            
               QEDGES(IQX,2) = Y(1)                                                
               QTANS(IQX,2)  = 0.0                                                 
            ELSE                                                                  
               J = 1                                                               
 460           CONTINUE                                                            
               IF (QXS(IQX).GT.X(J) .AND. J.LT.NJ) THEN                            
                  J = J + 1                                                         
                  GOTO 460                                                          
               ENDIF                                                               
               QEDGES(IQX,2) = Y(J) - (QXS(IQX)-X(J)) * (Y(J)-Y(J-1))/            
     >              ( X(J-1) -X(J))                              
               QEDGES(IQX,2) = MAX (QEDGES(IQX,2),0.0)                             
               THETA = ATAN (DS(J)- (QXS(IQX)-X(J)) * (DS(J)-DS(J-1))/            
     >              ( X(J-1) -X(J)))                               
               QTANS(IQX,2)  = PI/2.0 - GC + THETA                                 
C     WRITE (6,9001) IQX,QXS(IQX),QEDGES(IQX,2),THETA/DEGRAD,             
C     >        QTANS(IQX,2)/DEGRAD,J,T(J),S(J),X(J),Y(J)                         
            ENDIF                                                                 
 470     CONTINUE                                                                
C     WRITE (6,9000) JC,SC,TC,GC/DEGRAD,RP,ICUT(1),ICUT(2)                    
C     
C------OVERRIDE Y<0 SHAPE WITH REFLECTED Y>0 SHAPE FOR OPTION 5.               
C     
         IF (CIOPTH.EQ.5) THEN                                                   
            DO 500 IQX = 1-NQXSO, 0                                               
               QEDGES(IQX,1) = QEDGES(IQX,2)                                       
               QTANS (IQX,1) = QTANS (IQX,2)                                       
 500        CONTINUE                                                              
            ICUT(1) = ICUT(2)                                                     
         ENDIF                                                                   

C     
C     *********************************************************************
C     *  EDGE6: THE LIMITER IS A "STEP LIMITER" SPECIFIED BY THREE        *
C     *  POINTS. THE THREE POINTS ARE DEFINED SEPARATELY FOR EACH SIDE    *
C     *  OF Y=0.                                                          *
C     *                                                                   *
C     *  0.........(XST1,YST1)                                            *
C     *   :       :                                                       *
C     *   :       :                                                       *
C     *   :       :                                                       *
C     *   :       :..........(XST3,YST3)                                  *
C     *   :  (XST2,YST2)    :                                             *
C     *   :                 :                                             *
C     *   :.................:                                             *
C     *                                                                   *
C     *   NOTE: THIS IS OVER SPECIFIED. TOO MUCH INFORMATION IS PROVIDED  *
C     *   THE CODE SHOULD BE MODIFIED TO ALLOW A GENERALIZED 3-POINT      *
C     *   DEFINITION OF THE LIMITER SHAPE WITH A "STEP LIMITER" BEING     *
C     *   A DEGENERATE FORM OF IT.                                        *
C     *   D. ELDER   MAY 30 1990                                          *
C     *   UPGRADED TO GENERAL SPEC. JUNE 27,1990                          * 
C     *********************************************************************
C     
      ELSEIF (CIOPTH.EQ.6) THEN
         DO J = 1,2
            ICUT(J) = 1-NQXSO
            DO IQX = 1-NQXSO,0
               IF (QXS(IQX).LT.XST3(J)) THEN 
                  QEDGES(IQX,J) = YST3(J)
                  QTANS(IQX,J) = 0.0
               ELSEIF (QXS(IQX).LT.XST2(J)) THEN
                  QEDGES(IQX,J) = YST2(J) + (XST2(J)-QXS(IQX)) 
     >                 * (YST3(J)-YST2(J)) / (XST2(J)-XST3(J))
                  QTANS(IQX,J) = PI/2.0 - 
     >                 ATAN( (XST2(J)-XST3(J))/(YST3(J)-YST2(J)))
               ELSEIF (QXS(IQX).LT.XST1(J)) THEN
                  QEDGES(IQX,J) = YST1(J) + (XST1(J)-QXS(IQX)) 
     >                 * (YST2(J)-YST1(J)) / (XST1(J)-XST2(J))
                  QTANS(IQX,J) = PI/2.0 - 
     >                 ATAN( (XST1(J)-XST2(J))/(YST2(J)-YST1(J)))
               ELSEIF (QXS(IQX).GE.XST1(J)) THEN 
                  QEDGES(IQX,J) =  QXS(IQX) 
     >                 * YST1(J) / XST1(J)
                  QTANS(IQX,J) = PI/2.0 - 
     >                 ATAN(-XST1(J)/YST1(J))
               ENDIF
            end do 
         end do

C     
C     
C     *********************************************************************
C     * EDGE7: THE LIMITER IS A "BELT LIMITER" WITH THE A RADIUS OF       *
C     * CURVATURE IN THE SAME DIRECTION AS THAT FOR THE LCFS. IT IS       *
C     * TANGENTIAL TO THE LCFS AT THE POINT OF CONTACT AND IS SYMMETRIC   *
C     * FOR THE Y>0 AND Y<0 SIDES.                                        *
C     *                           Y=0                                     *
C     *    -----------------------|--------------------------- LCFS       *
C     *                ....... """"""".......                             *
C     *       ...."""""                      """""....                    * 
C     *    =================================================== WALL       *
C     *                                                                   *
C     * THE SHAPE IS DEFINED BY THE DISTANCE BETWEEN THE LCFS GIVEN BY    *
C     * RP(CA) AND THE LIMITER SURFACE GIVEN BY RL.                       *
C     *                                                                   *
C     * DAVID ELDER , JUNE 18 , 1990                                      *
C     *                                                                   *
C     *********************************************************************
C     
      ELSEIF (CIOPTH.EQ.7) THEN
C     
C     CALCULATE CONSTANT EXPRESSIONS
C     
         NUM1 = 2* (CA**2-RLEDGE7*CA)
         NUM2 = RLEDGE7-CA
         NUM3 = 1-(RLEDGE7/CA)
         NUM4 = NUM1/2
         NUM5 = 2*CA
         NUM6 = NUM2 * 2
         WRITE(6,*) 'N1,N2,N3,N4,N5,N6,XSO:',NUM1,NUM2,NUM3,NUM4,
     >        NUM5,NUM6,XSCALO
C     
         J = 1
         ICUT(J) = 1-NQXSO                                                     
         DO 850 IQX = 1-NQXSO, 0                                               
C     
C     NOTE : FOR CALCULATIONS VAL MUST BE > 0 
C     XSCALO < 0 AND IQX < 0 SO OK
C     
C     CALCULATE EDGE POSITIONS
C     

            VAL = -IQX/XSCALO
            INTVAL0 = (VAL**2 + NUM5*VAL + NUM1)/(NUM1-NUM6*VAL)
            QEDGES(IQX,J) = CA * ACOS(INTVAL0)
C     
C     CALCULATE THE NORMALS
C     
            INTVAL1 = COS ( QEDGES(IQX,J)/CA)
            INTVAL2 = SIN ( QEDGES(IQX,J)/CA)
            ALPHA = CA + NUM2 * INTVAL1
            BETA  = NUM4 * ( 1 - INTVAL1)
C     
            INTVAL3 = SQRT(ALPHA**2 - 2* BETA)   
C     
C     CALCULATE SLOPE          
C     
            SLOPE = (ALPHA/INTVAL3-1)*NUM3*INTVAL2 +
     >           INTVAL2*NUM2/INTVAL3
            QTANS(IQX,J) = PI/2.0 - ATAN(SLOPE)
C     
C     REFLECT Y<0 ONTO Y>0 
C     
            QEDGES(IQX,2) = QEDGES(IQX,1)
            QTANS(IQX,2)  = QTANS(IQX,1)
 850     CONTINUE               
         ICUT(2) = ICUT(1)

C     
C     *********************************************************************        
C     *  EDGE8: THE LIMITER IS SHAPED LIKE THE ARC OF A CIRCLE WITH CENTRE*        
C     *  AT (-RLC,0.0) AND RADIUS = RLC.                                  *        
C     *  THE ARC IS STOPPED AT XMAX= -RLC AND Y IS HELD CONSTANT          *   
C     *  FOR |X|>|XMAX|                                                   *
C     *  THE COORDINATES OF THE CENTRE OF THE ARC ARE GENERALIZED SO THAT *
C     *  FUTURE ADDITIONS CAN EXPAND UPON THIS LIMITER SHAPE.             *
C     *      
C     *       0....                                                       *        
C     *        :   :..                                                    *        
C     *        :      :.                                                  *        
C     *    XMAX:       :.                                                 *        
C     *        :        :                                                 *        
C     *        :        :                                                 *        
C     *        :        :                                                 *        
C     *      AW:........:                                                 *        
C     *        0       YMAX                                               *        
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.8) THEN                                                 
         X08 = -RLC
         Y08 = 0.0
         RADS2 = RLC * RLC                                                 
         DO J = 1, 2                                                         
            ICUT(J) = 1-NQXSO                                                     
            DO IQX = 1-NQXSO, 0                                               
               IF (QXS(IQX).GE.-RLC) THEN                                
                  QEDGES(IQX,J)=SQRT (RADS2-(QXS(IQX)-X08)**2) + Y08            
                  QTANS(IQX,J) =ATAN((QXS(IQX)-X08)/(QEDGES(IQX,J)-Y08))       
               ELSE                                                                
                  QEDGES(IQX,J) = RLC                                              
                  QTANS(IQX,J)  = 0.0                                               
               ENDIF                                                               
            end do 
         end do

C     
C     *********************************************************************        
C     *  EDGE9: The Limiter is composed of a set of 18 points which are   *
C     *  currently hrad coded. This can be changed to allow a more general*
C     *  specification of an N-point limiter shape. In addition, the set  *
C     *  of points could be interpolated with a spline or linear fit. At  *
C     *  this time a linear fit is used for the purposes of simplicity.   *
C     *  However, this could also be easily changed if necessary. The     *
C     *  limiter extends straight to the wall after the last point.       *
C     *                                                                   * 
C     *       0...X...                                                    *        
C     *        :     :. X...                5-point example               *        
C     *        :           :... X.                                        *        
C     *    XMAX:                 :..                                      *        
C     *        :                   :.X...                                 *        
C     *        :                        :..X                              *        
C     *        :                           :                              *        
C     *      AW:...........................:                              *        
C     *        0                         YMAX                             *        
C     *                                                                   *
C     *                                                                   *
C     *   David Elder Dec 18/1991                                         *
C     *   Improvements to be made : Generalized N -points and spline      *
C     *                             fitting ans asymmetry                 *
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.9) THEN 
C     
C     USE THE J VARIABLE INSTEAD OF HARD-CODING 1 OR 2 SO THAT 
C     MOVING TO AN ASYMMETRIC SITUATION WILL BE EASIER
C     
         J=1 
         ICUT(J) = 1-NQXSO
         ICUT(2) = 1-NQXSO
         DO 1000 IQX = 1-NQXSO,0
            IN = JPOS (QXS(IQX),XP,NP)  
            IF (IN.EQ.NP) THEN 
               QEDGES(IQX,J) = YP(IN)
               QTANS(IQX,J) = 0.0
            ELSE
               QEDGES(IQX,J) = YP(IN) + (XP(IN)-QXS(IQX)) 
     >              * (YP(IN+1)-YP(IN)) / (XP(IN)-XP(IN+1))
               QTANS(IQX,J) = PI/2.0 - 
     >              ATAN( (XP(IN)-XP(IN+1))/(YP(IN+1)-YP(IN)))
            ENDIF    
C     
C     MAKE SYMMETRIC FOR NOW
C     
            QEDGES(IQX,2) = QEDGES(IQX,1)
            QTANS(IQX,2) = QTANS(IQX,1)  
C     
 1000    CONTINUE 

C     
C     *********************************************************************        
C     *  EDGE10: ITER Startup limiter shape                               *
c     *          The separation between the LCFS and the limiter surface  *
c     *          is defined by the equation:                              * 
c     *          X = -(0.04 * Y**2 + 0.09 * Y**4)                         *
c     *          It is symmetric about Y=0                                *
C     *                                                                   * 
C     *       0...X...                                                    *        
C     *        :     :. X...                                              *        
C     *        :           :... X.                                        *        
C     *    XMAX:                 :..                                      *        
C     *        :                   :.X...                                 *        
C     *        :                        :..X                              *        
C     *        :                           :                              *        
C     *      AW:...........................:                              *        
C     *        0                         YMAX                             *        
C     *                                                                   *
C     *                                                                   *
C     *   David Elder Apr 6/2007                                          *
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.10) THEN 
C     
c     
c     In this option the distance between the LCFS and the limiter surface is defined by
c     the relation p = 0.04 * q**2 + 0.09 * q**4
c     where p is the radial distance to the limiter surface and q is the distance along
c     the field line from the limiter tip defined as q=0.
c     
c     These map to the LIM coordinates as q == Y and p == -X since the limiter surface 
c     is located at negative values of X.
c     
c     Keep in mind that the Y<0 and Y>0 limiter surfaces are specified with the same sign
c     thus a symmetric limiter has the same data stored in QEDGES and QTANS for each side 
c     of the limiter. 
c     
c     To calculate the limiter surface we need to invert the given expression to get q=f(p).
c     since QEDGES and QTANS are indexed by the "X" coordinate variable - thus Y=f(X) is required.
c     
c     The equation here has been parameterized with fixed constants though these could be 
c     changed to input values for greater flexibility if such is ever needed.
c     
         c1 = 0.04
         c2 = 0.09
C     
C     USE THE J VARIABLE INSTEAD OF HARD-CODING 1 OR 2 SO THAT 
C     MOVING TO AN ASYMMETRIC SITUATION WILL BE EASIER
C     
         J=1 
         ICUT(J) = 1-NQXSO
         ICUT(2) = 1-NQXSO
         DO IQX = 1-NQXSO,0
            xtmp = qxs(iqx)
            ytmp = sqrt((-c1 + sqrt(c1**2 + 4.0*c2*abs(xtmp)))/(2.0*c2))
c     write(6,'(a,i8,5(1x,g12.5))') 'calc:',iqx,xtmp,ytmp,
c     >         sqrt(c1**2+4.0*c2*abs(xtmp)) 
            qedges(iqx,j) = ytmp
         end do
c     
c     Calculate QTANS from QEDGES data
c     
c     Note: sign convention and calculation are copied from 
c     other entries to remain consistent.
c     
         DO IQX = 1-NQXSO,0
            if (iqx.eq.(1-nqxso)) then 
               qtans(iqx,j) = PI/2.0-atan((qxs(iqx)-qxs(iqx+1))/
     >              (qedges(iqx+1,j)-qedges(iqx,j)))
            elseif (iqx.eq.0) then 
c     qtans(iqx,j) = PI/2.0-atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                    qedges(iqx,j)-qedges(iqx-1,j))
               qtans(iqx,j) = PI/2.0-atan((qxs(iqx-1)-qxs(iqx))/
     >              (qedges(iqx,j)-qedges(iqx-1,j)))
            else 
c     theta1 = PI/2.0 - atan2c(qxs(iqx)-qxs(iqx+1),
c     >                                qedges(iqx+1,j)-qedges(iqx,j))
c     theta2 = PI/2.0 - atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                qedges(iqx,j)-qedges(iqx-1,j))
               theta1 = PI/2.0 - atan((qxs(iqx)-qxs(iqx+1))/
     >              (qedges(iqx+1,j)-qedges(iqx,j)))
               theta2 = PI/2.0 - atan((qxs(iqx-1)-qxs(iqx))/
     >              (qedges(iqx,j)-qedges(iqx-1,j)))
               qtans(iqx,j) = (theta1+theta2)/2.0
            endif
         end do
c     
c     Copy data over to other elements - make symmetric
c     
         qedges(:,2) = qedges(:,1)
         qtans(:,2) = qtans(:,1)
c     
c     
C     
C     *********************************************************************        
C     *  EDGE11: ITER Blanket module  shape                               *
c     *          The separation between the LCFS and the BM surface       *
c     *          is defined by the equation:                              * 
c     *                        
c     * X = - lambda_design ln(1 -/+ C*(t - t_re-entrant)/lambda_design ) *
c     *
c     * This limiter specification has two limiter tips at t=+/-t_re-entrant
c     * which configuration is not supported in LIM at the moment. In order
c     * to work around this the limiter tip is mapped to (0, t_re-entrant)
c     * and one half the limiter will become a slab at t=slot_tor_wid
c     *
c     * Code to handle a fully symmetric double tangency limiter shape is
c     * planned. 
c     * 
C     *                                                                   * 
C     *       0...X...                                                    *        
C     *        :     :. X...                                              *        
C     *        :           :... X.                                        *        
C     *    XMAX:                 :..                                      *        
C     *        :                   :.X...                                 *        
C     *        :                        :..X                              *        
C     *        :                           :                              *        
C     *      AW:...........................:                              *        
C     *        0                         YMAX                             *        
C     *                                                                   *
C     *                                                                   *
C     *   David Elder Apr 6/2007                                          *
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.11) THEN 
C     
c     
c     Note: In Peter's notes and the ITER nomenclature - the variables are called
c     t for the parallel to field line distance and x for the radial distance. 
c     
c     In this option the distance between the LCFS and the limiter surface is defined by
c
c     x(t) = - lambda_design ln(1 -/+ C*(t - t_re-entrant)/lambda_design ) 
c     
c     where x is the radial distance to the limiter surface and t is the distance along
c     the field line from the limiter tip defined as t=0.
c     
c     These map to the LIM coordinates as t == Y and x == X since the limiter surface 
c     is located at negative values of X.
c     
c     Keep in mind that the Y<0 and Y>0 limiter surfaces are specified with the same sign
c     thus a symmetric limiter has the same data stored in QEDGES and QTANS for each side 
c     of the limiter. 
c     
c     To calculate the limiter surface we need to invert the given expression to get t=f(x).
c     since QEDGES and QTANS are indexed by the "X" coordinate variable - thus Y=f(X) is required.
c     
c     Inverting the equation for x and t gives:
c     
c     t =  -lam/C  * ( exp(-x/lam) -1 ) + t_re-entrant  for t > t_re-entrant
c     t =   lam/C  * ( exp(-x/lam) -1 ) + t_re-entrant  for tslot < t < t_re_entrant
c     
c     So:
c     
c     Y =  -lam/C  * ( exp(-X/lam) -1 ) + Y_re-entrant  for Y > Y_re-entrant
c     Y =   lam/C  * ( exp(-X/lam) -1 ) + Y_re-entrant  for Yslot < Y < Y_re_entrant
c     
c     NOTE: limiter shape is symmetric about Y_re-entrant except for different cutoffs.           
c
c     NOTE: on the slot setback side the limiter shape becomes flat for Yslot > Y > 0
c           Ideally the slot at Y=0 should be closer to the separatrix than at the
c           start of the slot setback - however, for now the code won't support this
c           type of geometry so we will make it go almost flat relative to the LCFS.
c     
c     Calculate C and t_re-entrant parameters for the limiter shape from the following formulae:
c     
c     
c     y_re = t_re-entrant = bm_tor_wid * (1-exp(-rtor_setback/lambda_design)) + slot_tor_wid * (1-exp(-rslot_setback/lambda_design))/
c     (( 1-exp(-rtor_setback/lambda_design)) + (1-exp(-rslot_setback/lambda_design)))
c     
c     C = C_inner = C_outer = lambda_design * (1-exp(-rtor_setback/lambda_design)) / (bm_tor_wid - t_re-entrant)
c     
c     Since LIM does not support two limiter tips at present we can only model 1/2 of the limiter. We do this by mapping 
c     one limiter tip at y=y_re to y=0 and then we cut off the limiter at y=bm_tor_wid-y_re on one side and y=y_re-slot_tor_wid on the other.
c     A fully symmetric limiter option might be possible if we can address the issue of the multiple surface degeneracy as a function of X. 
c     
c     
         call calc_iter_limiter_parameters
c        
c         y_re = ( bm_tor_wid * (1-exp(-rtor_setback/lambda_design)) 
c     >        +slot_tor_wid*(1-exp(-rslot_setback/lambda_design)))/
c     >        ( (1-exp(-rtor_setback/lambda_design)) 
c     >         +(1-exp(-rslot_setback/lambda_design)))
c         c_lim = (lambda_design * (1-exp(-rtor_setback/lambda_design)))
c     >        / (bm_tor_wid - y_re)
c     
         write(6,'(a,10(1x,g12.5))') 'EDGE OPTION 11: BM Parameters:',
     >        rtor_setback,rslot_setback,bm_tor_wid,slot_tor_wid,
     >        lambda_design,y_re,c_lim
c     
c     Calculate QEDGES and QTANS: Extend limiter to Y=bm_tor_wid-y_re for "+"
c     Extend limiter to Y=slot_tor_wid-y_re for "-" 
c     

C     
C     USE THE J VARIABLE INSTEAD OF HARD-CODING 1 OR 2 SO THAT 
C     MOVING TO AN ASYMMETRIC SITUATION WILL BE EASIER
C     
c     
         ICUT(1) = 1-NQXSO
         ICUT(2) = 1-NQXSO
c     
c     First half
c     
         J=1    
         DO IQX = 1-NQXSO,0
            xtmp = qxs(iqx)
            ytmp = lambda_design/c_lim  * 
     >             (1.0-exp(-abs(xtmp)/lambda_design))
c
            if (xtmp.lt.-rtor_setback) then 
c            if (ytmp.gt.bm_tor_wid-y_re) then 
               qedges(iqx,j) = bm_tor_wid-y_re
            else
               qedges(iqx,j) = ytmp
            endif
c     
            write(6,'(a,2i8,5(1x,g12.5))') 'edge calc:',j,iqx,xtmp,ytmp,
     >                                  qedges(iqx,j),-rtor_setback
c     
         end do
c     
c     Second half
c     
         J=2    
         DO IQX = 1-NQXSO,0
            xtmp = qxs(iqx)
            ytmp = lambda_design/c_lim  * 
     >      (1.0-exp(-abs(xtmp)/lambda_design))
c
            if (xtmp.lt.-rslot_setback) then 
c
c            if (ytmp.gt.(y_re-slot_toar_wid)) then 
c
c               Adjust edge to vertical at y_re - this should give an almost
c               flat top to the slot - however, need to check how the code
c               works since I think qtans is calculated as the average of 
c               adjacent segments and the launches take place at points which 
c               could be problematic. 
c
c               qedges(iqx,j) = y_re-slot_tor_wid
c
               qedges(iqx,j) = y_re  
            else
               qedges(iqx,j) = ytmp
            endif
c     
            write(6,'(a,2i8,5(1x,g12.5))') 'edge calc:',j,iqx,xtmp,ytmp,
     >                                  qedges(iqx,j),-rslot_setback
c     
         end do

c     
c     Calculate QTANS from QEDGES data
c     
c     Note: sign convention and calculation are copied from 
c     other entries to remain consistent.
c     
         do j = 1,2
            DO IQX = 1-NQXSO,0
               if (iqx.eq.(1-nqxso)) then 
                  if ((qedges(iqx+1,j)-qedges(iqx,j)).eq.0.0) then 
                     qtans(iqx,j) = 0.0
                  else
                     qtans(iqx,j) = PI/2.0-atan((qxs(iqx)-qxs(iqx+1))/
     >                 (qedges(iqx+1,j)-qedges(iqx,j)))
                  endif

               elseif (iqx.eq.0) then 
c     qtans(iqx,j) = PI/2.0-atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                    qedges(iqx,j)-qedges(iqx-1,j))
                  if ((qedges(iqx,j)-qedges(iqx-1,j)).eq.0.0) then 
                     qtans(iqx,j) = 0.0
                  else
                     qtans(iqx,j) = PI/2.0-atan((qxs(iqx-1)-qxs(iqx))/
     >                 (qedges(iqx,j)-qedges(iqx-1,j)))
                  endif
c
               else 
c     theta1 = PI/2.0 - atan2c(qxs(iqx)-qxs(iqx+1),
c     >                                qedges(iqx+1,j)-qedges(iqx,j))
c     theta2 = PI/2.0 - atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                qedges(iqx,j)-qedges(iqx-1,j))
                  if ((qedges(iqx+1,j)-qedges(iqx,j)).eq.0.0) then 
                     theta1 = 0.0
                  else
                     theta1 = PI/2.0 - atan((qxs(iqx)-qxs(iqx+1))/
     >                 (qedges(iqx+1,j)-qedges(iqx,j)))
                  endif
                  if ((qedges(iqx,j)-qedges(iqx-1,j)).eq.0.0) then 
                     theta2 = 0.0
                  else
                     theta2 = PI/2.0 - atan((qxs(iqx-1)-qxs(iqx))/
     >                 (qedges(iqx,j)-qedges(iqx-1,j)))
                  endif
                  qtans(iqx,j) = (theta1+theta2)/2.0
               endif
            end do
         end do


c     
C     
C     *********************************************************************        
C     *  EDGE12: ITER Blanket module  shape                               *
c     *          The separation between the LCFS and the BM surface       *
c     *          is defined by the equation:                              * 
c     *                        
c     *  X = f(t) + g(p)
c     *  Y = (t - t_reentrant) / cos (beta)
c     * f(t) = - lambda_design ln(1 -/+ C*(t - t_re-entrant)/lambda_design )  
c     * g(p) = p**2 / (2* rho_pp)
c     *
c     * This limiter specification has two limiter tips at t=+/-t_re-entrant
c     * which configuration is not supported in LIM at the moment. In order
c     * to work around this the limiter tip is mapped to (0, t_re-entrant)
c     * and one half the limiter will become a slab at t=slot_tor_wid
c     *
c     * Code to handle a fully symmetric double tangency limiter shape is
c     * planned. 
c     * 
C     *                                                                   *
C     *   David Elder Oct 28, 2009                                        *
C     *                                                                   *        
C     *********************************************************************        
C     
      ELSEIF (CIOPTH.EQ.12) THEN 
C     
c     
c     Note: In Peter's notes and the ITER nomenclature - the variables are called
c     Y for the parallel to field line distance and X for the radial distance. 
c
c     The limiter surface is defined in terms of t and p
c     
c     In this option the distance between the LCFS and the limiter surface is defined by
c     the set of X,Y points determined from the following equations.
c
c      X = f(t) + g(p)
c      Y = (t - t_reentrant) / cos (beta)
c      f(t) = - lambda_design ln(1 -/+ C*(t - t_re-entrant)/lambda_design )  
c      g(p) = p**2 / (2* rho_pp)
c     
c     where X is the radial distance to the limiter surface and Y is the distance along
c     the field line from the limiter tip defined as Y=0.
c     
c     The Y axis stretches across the limiter surface which is defined in the t,p plane
c     
c     Keep in mind that the Y<0 and Y>0 limiter surfaces are specified with the same sign
c     thus a symmetric limiter has the same 'sign' data stored in QEDGES and QTANS for each side 
c     of the limiter. 
c     
c     It is not possible to analytically invert the X(Y) expression to obtain Y(X) needed by the 
c     code. As a result, the code calculates a function X,Y numerically defining the surface
c     and then uses interpolation to find appropriate values at the needed X coordinates on each side
c     of the limiter. 
c
c
c     NOTE: on the slot setback side the limiter shape becomes flat for +Yslot_setback < Y < +Y_reentrant
c     
c     Calculate C and t_re-entrant parameters for the limiter shape from the following formulae:
c     
c     
c     y_re = t_re-entrant = bm_tor_wid * (1-exp(-rtor_setback/lambda_design)) + slot_tor_wid * (1-exp(-rslot_setback/lambda_design))/
c     (( 1-exp(-rtor_setback/lambda_design)) + (1-exp(-rslot_setback/lambda_design)))
c     
c     C = C_inner = C_outer = lambda_design * (1-exp(-rtor_setback/lambda_design)) / (bm_tor_wid - t_re-entrant)
c     
c     Since LIM does not support two limiter tips at present we can only model 1/2 of the limiter. We do this by mapping 
c     one limiter tip at y=y_re to y=0 and then we cut off the limiter at y=bm_tor_wid-y_re on one side and y=y_re-slot_tor_wid on the other.
c     A fully symmetric limiter option might be possible if we can address the issue of the multiple surface degeneracy as a function of X. 
c     
c     

         call calc_iter_limiter_parameters
c        
c         y_re = ( bm_tor_wid * (1-exp(-rtor_setback/lambda_design)) 
c     >        +slot_tor_wid*(1-exp(-rslot_setback/lambda_design)))/
c     >        ( (1-exp(-rtor_setback/lambda_design)) 
c     >         +(1-exp(-rslot_setback/lambda_design)))
c         c_lim = (lambda_design * (1-exp(-rtor_setback/lambda_design)))
c     >        / (bm_tor_wid - y_re)
c     
         write(6,'(a,10(1x,g12.5))') 'EDGE OPTION 12: BM Parameters:',
     >        rtor_setback,rslot_setback,bm_tor_wid,slot_tor_wid,
     >        lambda_design,y_re,c_lim

c
c        Calculate the limiter shape function
c
c     
c     Calculate QEDGES and QTANS: Extend limiter to Y=bm_tor_wid-y_re for "+"
c     Extend limiter to Y=slot_tor_wid-y_re for "-" 
c     
         call calc_iter_limiter_shape

C     
C     USE THE J VARIABLE INSTEAD OF HARD-CODING 1 OR 2 SO THAT 
C     MOVING TO AN ASYMMETRIC SITUATION WILL BE EASIER
C     
c     
         ICUT(1) = 1-NQXSO
         ICUT(2) = 1-NQXSO
c     
c     Apply the limiter shape to the First half Y<0 face though all Y values used are >0
c     
         J=1    
         DO IQX = 1-NQXSO,0
c
            xtmp = qxs(iqx)           
c
c            ytmp = lambda_design/c_lim  * (1.0-exp(-xtmp/lambda_design))
c
            if (xtmp.lt.xtor_setback) then 
c
c            if (ytmp.gt.bm_tor_wid-y_re) then 
c               
               qedges(iqx,j) = bm_tor_wid_y
            else
               qedges(iqx,j) = iter_limiter_shape(xtmp,j)
            endif
c     
c     write(6,'(a,2i8,5(1x,g12.5))') 'edge calc:',j,iqx,xtmp,ytmp
c     
         end do

c     
c     Apply the limiter shape to the Second half
c     
         J=2    
         DO IQX = 1-NQXSO,0
c
            xtmp = qxs(iqx)
c
c            ytmp = lambda_design/c_lim  * (1.0-exp(-xtmp/lambda_design))
c
            if (xtmp.lt.xslot_setback) then 
c
c            if (ytmp.gt.(y_re-slot_tor_wid)) then 
c
c               Adjust edge to vertical at y_re-slot_tor_wid - this should give an almost
c               flat top to the slot - however, need to check how the code
c               works since I think qtans is calculated as the average of 
c               adjacent segments and the launches take place at points which 
c               could be problematic. 
c
c               qedges(iqx,j) = y_re-slot_tor_wid
c
               qedges(iqx,j) = slot_tor_wid_y
            else
               qedges(iqx,j) = iter_limiter_shape(xtmp,j)
            endif
c     
c     write(6,'(a,2i8,5(1x,g12.5))') 'edge calc:',j,iqx,xtmp,ytmp
c     
         end do

c
c       After the limiter shape has been stored in QEDGES - clean up the limiter module
c
         call clean_up_iter_limiter_shape

c     
c     Calculate QTANS from QEDGES data
c     
c     Note: sign convention and calculation are copied from 
c     other entries to remain consistent.
c     
         do j = 1,2
            DO IQX = 1-NQXSO,0
               if (iqx.eq.(1-nqxso)) then 
                  if ((qedges(iqx+1,j)-qedges(iqx,j)).eq.0.0) then 
                     qtans(iqx,j) = 0.0
                  else
                     qtans(iqx,j) = PI/2.0-atan((qxs(iqx)-qxs(iqx+1))/
     >                 (qedges(iqx+1,j)-qedges(iqx,j)))
                  endif

               elseif (iqx.eq.0) then 
c     qtans(iqx,j) = PI/2.0-atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                    qedges(iqx,j)-qedges(iqx-1,j))
                  if ((qedges(iqx,j)-qedges(iqx-1,j)).eq.0.0) then 
                     qtans(iqx,j) = 0.0
                  else
                     qtans(iqx,j) = PI/2.0-atan((qxs(iqx-1)-qxs(iqx))/
     >                 (qedges(iqx,j)-qedges(iqx-1,j)))
                  endif
c
               else 
c     theta1 = PI/2.0 - atan2c(qxs(iqx)-qxs(iqx+1),
c     >                                qedges(iqx+1,j)-qedges(iqx,j))
c     theta2 = PI/2.0 - atan2c(qxs(iqx-1)-qxs(iqx),
c     >                                qedges(iqx,j)-qedges(iqx-1,j))
                  if ((qedges(iqx+1,j)-qedges(iqx,j)).eq.0.0) then 
                     theta1 = 0.0
                  else
                     theta1 = PI/2.0 - atan((qxs(iqx)-qxs(iqx+1))/
     >                 (qedges(iqx+1,j)-qedges(iqx,j)))
                  endif
                  if ((qedges(iqx,j)-qedges(iqx-1,j)).eq.0.0) then 
                     theta2 = 0.0
                  else
                     theta2 = PI/2.0 - atan((qxs(iqx-1)-qxs(iqx))/
     >                 (qedges(iqx,j)-qedges(iqx-1,j)))
                  endif
                  qtans(iqx,j) = (theta1+theta2)/2.0
               endif
            end do
         end do



      ENDIF
c     

c     
C     
C     *********************************************************************        
C     *  CHANGE ICUT IF CCUT PARAMETER SPECIFIED ON INPUT STREAM REQUIRES *        
C     *  A DIFFERENT CUTOFF POSITION ... RESET CCUT TO HOLD EXACT VALUE.  *        
C     *********************************************************************        
C     
      DO  J = 1, 2                                                           
         DO IQX = ICUT(J), 0                                                 
            IF (QXS(IQX).LT.CCUT) ICUT(J) = ICUT(J) + 1                           
         end do
      end do
C     
C     *********************************************************************        
C     *  CORRECT LIMITER SURFACE FOR PLASMA CURVATURE - NOTE 197.         *        
C     *  MAKE TEMPORARY COPY OF OLD LIMITER EDGE POSITIONS IN +Y REGION.  *        
C     *  CALCULATE DISTANCES DR FROM LCFS                                 *        
C     *  REPEAT ENTIRE PROCESS FOR OTHER HALF OF LIMITER IN -Y REGION     *        
C     *  18/11/88 ADD TEMP(1),DR(1) DUMMY POINT SO THAT CURVATURE EFFECT  *        
C     *  IS CORRECTLY APPLIED AT THE LIMITER TIP.                         *        
C     *********************************************************************        
C     
      IF (CORECT.EQ.1) THEN                                                     
         DO 930 J = 1, 2                                                         
            DO 910 IQX = 1-NQXSO, 0                                               
               TEMP(IQX) = QEDGES(IQX,J)                                           
               DR  (IQX) =-SQRT (QEDGES(IQX,J)**2 
     >                     + (QXS(IQX)-RP)**2) + RP         
 910        CONTINUE                                                              
            TEMP(1) = 0.0                                                         
            DR  (1) = 0.0                                                         
            CALL FITTER (NQXSO+1,DR(1-NQXSO),TEMP(1-NQXSO),                       
     >           NQXSO,QXS(1-NQXSO),QEDGES(1-NQXSO,J),'LINEAR')           
 930     CONTINUE                                                                
      ENDIF                                                                     
C     
C     *********************************************************************        
C     *  CALCULATE DISTANCES ALONG LIMITER SURFACE                        *        
C     *********************************************************************        
C     
      FACT = (1.0/XSCALO)**2                                                    
      QDISTS(0,1) = SQRT (QEDGES(0,1)**2 + FACT/4.0)                            
      QDISTS(0,2) = SQRT (QEDGES(0,2)**2 + FACT/4.0)                            
      DO 950 IQX = -1, 1-NQXSO, -1                                              
         QDISTS(IQX,1) = QDISTS(IQX+1,1) +                                       
     >        SQRT ((QEDGES(IQX,1)-QEDGES(IQX+1,1))**2 + FACT)                      
         QDISTS(IQX,2) = QDISTS(IQX+1,2) +                                       
     >        SQRT ((QEDGES(IQX,2)-QEDGES(IQX+1,2))**2 + FACT)                      
 950  CONTINUE                                                                  
C     
C     *********************************************************************        
C     *  CALCULATE OYS,ODS ARRAYS FOR NET EROSION OUTPUTS                 *        
C     *********************************************************************        
C     
      OYMAX = 0.0                                                               
      DO J = 1, 2                                                           
         DO IQX = 1-NQXSO, 0                                                 
            OYMAX2(J) = MAX (OYMAX2(J), QEDGES(IQX,J))                        
         end do
      end do

      OYMAX = MAX(OYMAX2(1),OYMAX2(2))
      WRITE(6,'(a,3(1x,g12.5))') 
     >     'OYMAX,OYM2(1),OYM2(2):',OYMAX,OYMAX2(1),OYMAX2(2) 

      ODMAX = MAX (QDISTS(1-NQXSO,1), QDISTS(1-NQXSO,2))                        
      WRITE(6,'(a,3(1x,g12.5))') 
     >     'ODMAX,QD(1-N,1),QD(1-N,2):',ODMAX,
     >       QDISTS(1-NQXSO,1), QDISTS(1-NQXSO,2)

      DO IO = 1, MAXOS                                                      

c
c       jdemod - issues here
c       1) Why the 1.05 factor? Makes no sense and simply
c          expands both Y and Distance along surface by a constant factor
c       2) ODS/OYS are bin boundaries while OYOUTS and ODOUTS should be bin 
c          centers. However, for some reason, the code offsets the bounds
c          by half a cell width and the cell center by a second 1/2 cell
c       3) OYWIDS and ODWIDS are loop constants but are recalculated on every
c          loop interation? Maybe just to more easily assign it? 
c
c        Cell bounds should extend from -OYMAX to OYMAX and -ODMAX to ODMAX
c        - need to check IPOS to see what the value in the first element needs
c          to be but probably ODS(-MAXOS) = -ODMAX+ODWIDS
c
c        OYS, OYOUTS may also not be correct but are used in NEUT primarily and
c        there may be code in NEUT to compensate for the defintion (not re-writing
c        NEUT at the moment). 
c
c         OYWIDS(IO) = 1.05 * OYMAX / REAL(MAXOS/2)                               
c
         OYWIDS(IO) = OYMAX / REAL(MAXOS/2)                               
         OYS(IO)    = (REAL(IO-MAXOS/2)-0.5) * OYWIDS(IO)                        
         OYOUTS(IO) = OYS(IO) - 0.5 * OYWIDS(IO)                                 

c         ODWIDS(IO) = 1.05 * ODMAX / REAL(MAXOS/2)                               
c         ODS(IO)    = (REAL(IO-MAXOS/2)-0.5) * ODWIDS(IO)                        
c         ODOUTS(IO) = ODS(IO) - 0.5 * ODWIDS(IO)                                 

         ODWIDS(IO) = ODMAX / REAL(MAXOS/2)                               
c
c        Set ODS to upper bin bounds and ODOUTS to bin center
c
         ODS(IO)    = REAL(IO-MAXOS/2) * ODWIDS(IO)                        
         ODOUTS(IO) = ODS(IO) - 0.5 * ODWIDS(IO)                                 
         
      end do

c
C     
      if (cprint.eq.1.or.cprint.ge.9) then 
         WRITE (6,9004) OYMAX,ODMAX,(IO,OYWIDS(IO),OYS(IO),OYOUTS(IO),             
     >      ODWIDS(IO),ODS(IO),ODOUTS(IO),IO=1,MAXOS)                               
         WRITE (6,9003) (QXS(IQX),(PI-QTANS(IQX,1))/DEGRAD,                        
     >     -QDISTS(IQX,1),-QEDGES(IQX,1),QEDGES(IQX,2),QDISTS(IQX,2),              
     >      QTANS(IQX,2)/DEGRAD,IQX=0,1-nqxso,-1)                                       
      endif
c
      RETURN                                                                    
C     
c     I/O Formatting
c
 9000 FORMAT(1X,'EDGE: JC,SC,TC,GC,RP,ICUT()=',I5,1P,4G11.4,2I5)                
 9001 FORMAT(1X,'EDGE: IQX',I6,' QXS',F7.4,' EDGE',F7.4,                        
     >     ' THETA',F7.2,' TAN',F7.2,' J',I4,' T',F7.2,' S',F7.2,                  
     >     ' X',F7.4,' Y',F7.4)                                                    
 9002 FORMAT(1X,'EDGE: J',I4,' INT',I3,                                         
     >     ' X',F7.4,' Y',F7.4,' S',F7.2,' T',F7.2,' DS',F7.2)                     
 9003 FORMAT(//1X,                                                              
     >'       X      THETA1      DIST1       EDGE1       EDGE2  ',
     >     '     DIST2    THETA2',/1X,                                             
     >     79('-'),/,(1X,F11.7,F9.3,4F12.7,F9.3))                                  
 9004 FORMAT(//1X,'EDGE: OYMAX=',F9.4,'  ODMAX=',F9.4,//1X,                     
     >  '  IO     OYWIDS   OYS    OYOUTS    ODWIDS   ODS    ODOUTS',/1X,        
     >           65('-'),/,(1X,I5,6F9.4))                                                
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE EDGINT (X,IQX,J,E,D)                                           
      use error_handling
      use mod_comt2
      IMPLICIT none
      REAL    X,E,D                                                             
      INTEGER IQX,J                                                             
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  EDGINT:  ROUTINE TO RETURN AN "EXACT" EDGE Y COORDINATE FROM A   *        
C  *  GIVEN X COORDINATE WHICH LIES SOMEWHERE WITHIN THE REGION        *        
C  *  DEFINED BY THE INDEX IQX, ON SIDE J  (J=1 Y<0,  J=2 Y>0).        *        
C  *  THE ROUTINE INTERPOLATES THROUGH THE SET OF CALCULATED EDGE      *        
C  *  POSITIONS STORED IN QEDGES ARRAY.                                *        
C  *  THE CORRESPONDING DISTANCE ALONG THE LIMITER SURFACE IS ALSO     *        
C  *  RETURNED.                                                        *        
C  *  NOTE: ALWAYS RETURNS POSITIVE VALUES,  EVEN FOR SIDE Y<0.        *        
C  *                                                                   *        
C  *            CHRIS FARRELL  (HUNTERSKIL)  FEBRUARY 1989             *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'params'                                                          
C     INCLUDE (PARAMS)                                                          
      INCLUDE 'comxyt'                                                          
C     INCLUDE (COMXYT)                                                          
c      INCLUDE 'comt2'                                                            
C     INCLUDE (COMT2)                                                           
      REAL    QT,QB,ET,EB,DT,DB                                                 
C                                                                               
c      write(0,'(a,2i8,10(1x,g12.5))') 'EDGINT1:',iqx,j,x,e,d
c
      IF (X.GE.0.0) THEN                                                        
        E = 0.0                                                                 
        D = 0.0                                                                 
        GOTO 999                                                                
      ELSEIF (X.GT.QXS(0)) THEN                                                 
        QT = 0.0                                                                
        QB = QXS(0)                                                             
        ET = 0.0                                                                
        EB = QEDGES(0,J)                                                        
        DT = 0.0                                                                
        DB = QDISTS(0,J)                                                        
      ELSEIF (X.LT.QXS(1-NQXSO)) THEN                                           
        QT = QXS(1-NQXSO)                                                       
        QB = QXS(1-NQXSO) - 1.0/XSCALO                                          
        ET = QEDGES(1-NQXSO,J)                                                  
        EB = QEDGES(1-NQXSO,J)                                                  
        DT = QDISTS(1-NQXSO,J)                                                  
        DB = QDISTS(1-NQXSO,J)                                                  
      ELSEIF (X.GT.QXS(IQX)) THEN                                               
        QT = QXS(IQX+1)                                                         
        QB = QXS(IQX)                                                           
        ET = QEDGES(IQX+1,J)                                                    
        EB = QEDGES(IQX,J)                                                      
        DT = QDISTS(IQX+1,J)                                                    
        DB = QDISTS(IQX,J)                                                      
      ELSE                                                                      
        QT = QXS(IQX)                                                           
        QB = QXS(IQX-1)                                                         
        ET = QEDGES(IQX,J)                                                      
        EB = QEDGES(IQX-1,J)                                                    
        DT = QDISTS(IQX,J)                                                      
        DB = QDISTS(IQX-1,J)                                                    
      ENDIF                                                                     
      
      E = ET + (QT-X) / (QT-QB) * (EB-ET)                                       
      D = DT + (QT-X) / (QT-QB) * (DB-DT)                                       

c      write(0,'(a,20(1x,g12.5))') 'EDGINT2:',x,qxs(0),qxs(iqx),
c     >     qt,qb,et,eb,dt,db

C                                                                               
  999 RETURN                                                                    
      END                                                                       
C
C
C
      SUBROUTINE YEDGINT(Y,X,IQX0,J,IERR)
      use mod_comt2
      use mod_comnet
      implicit none
      REAL     Y,X
      INTEGER  IQX0,J,IERR
C
      INCLUDE 'params'
C
      INCLUDE 'comxyt'
C
c      INCLUDE 'comt2'
C
c      INCLUDE 'comnet'
C
C**************************************************************
C
C     YEDGINT: THIS ROUTINE CALCULATES THE X-POSITION ON THE 
C     LIMITER GIVEN A STARTING Y-POSITION. IT MAKES A BINARY 
C     SEARCH ON THE QEDGES ARRAY UNTIL IT INTERSECTS 
C     THE LIMITER. IF IT DOES NOT INTERSECT, A LAUNCH POSITION 
C     OF Y,0.0 AND AN IQX VALUE OF 0 ARE RETURNED. IT IS 
C     ASSUMED THAT, IN GENERAL, THE LIMITER SHAPE WILL BE TAKEN 
C     INTO CONSIDERATION WHEN THE +/- Ycf FOR LAUNCH IS SPECIFIED
C     AND THUS A NON-INTERSECTING CONDITION WILL NOT OCCUR. 
C     THE ROUTINE USES THE LIMITER COORDINATES STORED IN THE 
C     QEDGES ARRAY AND THEN LINEARLY INTERPOLATES BETWEEN THE 
C     POINTS TO OBTAIN A MORE ACCURATE ESTIMATE OF THE X-LAUNCH 
C     POSITION.
C     THIS ROUTINE ASSUMES THAT THE QEDGES VALUES ARE MONOTONICALLY
C     INCREASING AS ONE MOVES FROM IQX=0 TO IQX = 1-NQXSO.
C
C     DAVID ELDER, OCT 17 1990, JAN 30 1991
C
C***************************************************************
C
      INTEGER  ISTRT,IEND,IMID,ITER,LIMIT 
      REAL     ABSY
      real     frac
C
C      WRITE(6,*) 'Y: ',Y
C      WRITE(6,*) 'OYMAX,OYM2(1),OYM2(2):',OYMAX,OYMAX2(1),
C     >            OYMAX2(2)
C 
      IQX0 = 0
      IERR = 0
      ISTRT = 0 
      LIMIT = 1-NQXSO 
      IEND = LIMIT -1 

      IF (Y.GE.0.0) THEN 
         J =2
      ELSE 
         J = 1
      ENDIF
C
      ABSY = ABS(Y)
      IF (ABSY.GT.OYMAX2(J)) THEN 
         IERR = 1
         IQX0 = 0
         X= 0.0
         WRITE(6,'(''YEDGINT ERROR: X LAUNCH POSITION'//
     >             ' FOR Y='',G14.4,
     >    '' NOT FOUND: (X,Y) = (0,0, '',G14.4,
     >    '') ASSUMED'')') Y,Y
      ELSE
         ITER = 0
100      IMID = INT((ISTRT+IEND)/2)            
         IF (QEDGES(IMID,J).GT.ABSY) THEN 
            IEND =  IMID
         ELSE
            ISTRT = IMID
         ENDIF
         ITER = ITER + 1    
         IF ((ISTRT-1).NE.IEND.AND.ITER.LT.200) GOTO 100 
C
         IQX0 = ISTRT
c
c         WRITE(6,*) 'YEDGINT: IQX0:',IQX0  
C 
         IF (ITER.GE.200) THEN            
            X = QEDGES(IQX0,J)             
            WRITE(6,'(''YEDGINT ERROR: TOO MANY ITERATIONS''//
     >            '' CLOSEST X ASSUMED : COORDS SUPPLIED ('',G14.4,
     >            '','',G14.4,'')'')')   X,Y          
         ELSE
            IF (ISTRT.LE.LIMIT) THEN
               IQX0 = LIMIT 
               X = QEDGES(IQX0,J)
               WRITE(6,'(a,1x,2i8,1x,g12.5)') 
     >               'DEBUG YINT: LIMIT : ISTRT,X :',LIMIT,
     >                    ISTRT,X 
           ELSE
               FRAC = 
     >          (ABSY-QEDGES(ISTRT,J))/(QEDGES(IEND,J)-QEDGES(ISTRT,J))
               X = (IQX0 - FRAC) / XSCALO
               WRITE (6,'(a,2x,g12.5,2i8,5(1x,g12.5))') 
     >          'DEBUG YINT: Y,IQX0,J,QE(IQX0),QE(IQX0-1),FRAC,X:' 
     >            ,Y,IQX0,J,QEDGES(IQX0,J),QEDGES(IQX0-1,J),FRAC,X
            ENDIF
         ENDIF    
      ENDIF
C
      RETURN
      END
