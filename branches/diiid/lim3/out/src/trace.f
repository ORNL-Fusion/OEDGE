C@PROCESS VECTOR(LEV(0)),OPT(1)                                                
C                                                                               
C                                                                               
      SUBROUTINE LIM_DRAW (AS,WS,BS,MAXNAS,NAS,ANLY,                                
     >  NBS,ISMOTH,ASTART,AEND,BSTART,BEND,IGS,ITEC,AVS,NAVS,                   
     >  JOB,TITLE,AAXLAB,BAXLAB,BLABS,REF,VIEW,PLANE,TABLE,IDRAW,IFLAG)         
      IMPLICIT none
      include 'params'
      INTEGER   MAXNAS,NAS,IBS,NBS,IDRAW,ISMOTH,IFLAG,IGS(*),ITEC,NAVS          
      REAL      AS(*),WS(*),BS(MAXNAS,*),ASTART,AEND,BSTART,BEND                
      REAL      AVS(0:NAVS)                                                     
      CHARACTER JOB*72,TITLE*80,TABLE*36,ANLY*36                                
      CHARACTER AAXLAB*24,BAXLAB*24,BLABS(*)*36,REF*36,VIEW*36,PLANE*36         
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  * DRAW:    DRAW A GRAPH OF VALUES IN BS AGAINST VALUES IN AS, FOR   *        
C  *          A VALUES WITHIN THE RANGE ASTART TO AEND                 *        
C  *                                                                   *        
C  * ARGUMENTS :-                                                      *        
C  *                                                                   *        
C  * AS     - A VALUES (TO GO ON X AXIS OF GRAPH)                      *        
C  * WS     - A BIN WIDTHS  (FOR SMOOTHING PURPOSES)                   *        
C  * BS     - B VALUES (TO GO ON Y AXIS OF GRAPH)                      *        
C  * MAXNAS - DIMENSION FOR AS ARRAY, FIRST DIMENSION OF BS ARRAY      *        
C  * NAS    - NUMBER OF AS VALUES USED                                 *        
C  * NBS    - LAST BS VALUE TO BE PLOTTED                              *        
C  * ASTART - MINIMUM AS VALUE TO BE PLOTTED                           *        
C  * AEND   - MAXIMUM AS VALUE TO BE PLOTTED                           *        
C  * BSTART - MINIMUM BS VALUE FOR AXIS: ACTUALLY USE MAX (BSTART, BS) *        
C  * BEND   - MAXIMUM BS VALUE FOR AXIS: ACTUALLY USE MIN (BEND, BS)   *        
C  *          BSTART AND BEND ONLY APPLY WHEN IDRAW=1, IE UNNORMALISED.*        
C  * IGS    - LIST OF CURVES TO BE IGNORED (NOT PLOTTED)               *        
C  * ITEC   - SMOOTHING TECHNIQUE  0:SPLINE  1:WEIGHTED AVERAGE        *        
C  * AVS    - SET OF WEIGHTS FOR WEIGHTED AVERAGE SMOOTHING            *        
C  * NAVS   - DIMENSION OF AVS.  "AVS(-1)" ASSUMED EQUAL TO AVS(1) ETC *        
C  * JOB    - JOB REFERENCE STRING INCORPORATING TIME/DATE ETC         *        
C  * TITLE  - TITLE ABOVE GRAPH                                        *        
C  * AAXLAB - X-AXIS LABEL                                             *        
C  * BAXLAB - Y-AXIS LABEL                                             *        
C  * BLABS  - LABELS FOR EACH SET OF BS VALUES                         *        
C  * REF    - LABEL FOR GRAPH  (BOTTOM RIGHT HAND CORNER)              *        
C  * VIEW   - ANY EXTRA NOTE TO BE PRINTED IN SYMBOL TABLE (TYPICALLY  *        
C  *          NOTE ON OBSERVATION POINT ETC)                           *        
C  * PLANE  - DITTO  (TYPICALLY NOTE ON 3D PLANE ETC)                  *        
C  * TABLE  - USUALLY "SYMBOL TABLE", BUT CAN BE CHANGED IF REQUIRED   *        
C  * IDRAW  - 1 NORMAL CASE, 2 NORMALISE GRAPHS TO UNITY AND PRINT     *        
C  *          A SCALING FACTOR FOR EACH IONISATION LEVEL.  NOT ALLOWED *        
C  *          IDRAW=2 FOR CASES WHERE FUNCTION VALUE < 0 OCCURS.       *        
C  *          IDRAW=3 AS 2, EXCEPT PARTIAL AREA CALCULATED, NOTE 238.  *        
C  * ISMOTH - INDICATES WHICH GRAPHS TO SMOOTH, IE GRAPHS "ISMOTH" AND *        
C  *          ABOVE OUT OF "NBS" GRAPHS IN ALL.                        *        
C  *          SET ISMOTH > NBS FOR NO SMOOTHING ANYWHERE.              *        
C  * IFLAG  - 1  CALL GRTSET (AND DRAW AXES), PLOT GRAPHS, IGNORE FRAME*        
C  *          2  CALL GRTSET (AND DRAW AXES), PLOT GRAPHS, CALL FRAME  *        
C  *          3  CALL GRTSET (IGNORING AXES), PLOT GRAPHS, CALL FRAME  *        
C  *          4  SUPERIMPOSE GRAPH ON PREVIOUS FRAME, DONT CALL FRAME  *        
C  *          5  SUPERIMPOSE GRAPH ON PREVIOUS FRAME, THEN CALL FRAME  *        
C  *          6  AS 2, BUT DRAW EXTRA LINE ALONG F = 0  (USED EG       *        
C  *             FOR NET EROSION GRAPHS)                               *        
C  *          7  AS 2, BUT DIVIDE TOTAL AREA BY 2 (FOR Y PLOTS)        *        
C  *                                                                   *        
C  * THE PARAMETERS MXXNAS,MXXNBS ARE USED FOR DIMENSIONING LOCAL ARRAY*        
C  * USED WITH THE CURVE SMOOTHING OPTION.  THEY ARE SET LARGE TO ALLOW*        
C  * ROUTINE DRAW TO BE CALLED FROM PROGRAMS OTHER THAN LIM FOR EXAMPLE*        
C  * FOR LIM USE, MXXNAS >= MAXNAS,   MXXNBS >= MAXNBS                 *        
C  *                                                                   *        
C  *                     CHRIS FARRELL    SEPTEMBER 1988               *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INTEGER IA,IB,NKNOTS,MXXNAS,MXXNBS,IA1,IA2,IA3,IA4,KOUNT,J,JA,IPOS        
      INTEGER NBBS                                                              
      PARAMETER (MXXNAS=10200,MXXNBS=20)                                        
      REAL    BMIN,BMAX,FACTS(MXXNBS),NORMS(MXXNAS),CS(MXXNAS,MXXNBS)           
      REAL    RD(MXXNAS),FN(MXXNAS),AREAS(MXXNBS)                               
      REAL    GN(MXXNAS),DN(MXXNAS),THETA(MXXNAS),XN(MXXNAS),TG01B              
      REAL    WORKS(16*MXXNAS+6),WMIN,TOTAV,AENDS(MXXNAS)                       
      INTEGER RANGE,ICOUNT
      CHARACTER SMOOTH*72                                                       
      COMMON /TRACE/ FACTS,NORMS,AREAS,CS,AENDS,                                
     >               RD,XN,FN,GN,DN,THETA,WORKS                                 
C                                                                               
      COMMON /NSMOOTH/ NUMSMOOTH
      INTEGER NUMSMOOTH
C
      WRITE (6,'(/,'' DRAW: '',A36,/1X,42(''-''))') REF                         
c
      write(6,*) 'REF   :',trim(ref),':'
      write(6,*) 'TABLE :',trim(table),':'
      write(6,*) 'VIEW  :',trim(view),':'
      write(6,*) 'PLANE :',trim(plane),':'
      write(6,*) 'ANLY  :',trim(anly),':'
      write(6,*) 'JOB   :',trim(job),':'
      write(6,*) 'AAXLAB:',trim(aaxlab),':'
      write(6,*) 'BAXLAB:',trim(baxlab),':'
      write(6,*) 'BLABS1:',trim(blabs(1)),':'
c
c     Print the plot data to the .plt file ... outunit = 49
c
      WRITE (outunit,'(/,'' DRAW: '',A36,/1X,42(''-''))') REF                         
c
      write(outunit,*) 'REF   :',trim(ref),':'
      write(outunit,*) 'TABLE :',trim(table),':'
      write(outunit,*) 'VIEW  :',trim(view),':'
      write(outunit,*) 'PLANE :',trim(plane),':'
      write(outunit,*) 'ANLY  :',trim(anly),':'
      write(outunit,*) 'JOB   :',trim(job),':'
      write(outunit,*) 'AAXLAB:',trim(aaxlab),':'
      write(outunit,*) 'BAXLAB:',trim(baxlab),':'
      write(outunit,*) 'BLABS1:',trim(blabs(1)),':'


c
c      
c jdemod - removing printout
c
c      WRITE (6,'('' CHECK PARMS: MAXNAS    NAS    NBS MXXNAS MXXNBS'',          
c     >  /13X,5I7)') MAXNAS,NAS,NBS,MXXNAS,MXXNBS                                
C                                                                               
C-----------------------------------------------------------------------        
C     GENERATE CS ARRAY TO BE USED FOR THE ACTUAL PLOTTING                      
C-----------------------------------------------------------------------        
C                                                                               
C---- SET RANGE OF AS VALUES FOR PLOTTING WINDOW (IA1,IA2)                      
C---- SET RANGE OF AS VALUES FOR CUBIC SPLINE INTERPOLATION (IA3,IA4)           
C---- SLIGHTLY LARGER SO CURVES ARE SMOOTH AT EDGES OF VISIBLE REGION.          
C                                                                               
      IA1 = 1                                                                   
      IA2 = 0                                                                   
      DO 10 IA = 1, NAS                                                         
        IF (AS(IA).LT.ASTART) IA1 = IA1 + 1                                     
        IF (AS(IA).LT.AEND  ) IA2 = IA2 + 1                                     
        AENDS(IA) = AS(IA) + 0.5 * WS(IA)                                       
        DO 10 IB = 1, NBS                                                       
          CS(IA,IB) = BS(IA,IB)                                                 
   10 CONTINUE                                                                  
      IA1 = MAX (1,   IA1-1)                                                    
      IA2 = MIN (NAS, IA2+1)                                                    
      IA3 = MAX (1,   IA1-5)                                                    
      IA4 = MIN (NAS, IA2+5)                                                    
      WMIN = 1.E10                                                              
      DO 15 IA = IA3, IA4                                                       
        IF (WS(IA).GT.0.0) WMIN = MIN (WMIN, WS(IA))                            
   15 CONTINUE                                                                  
      TOTAV = 0.0                                                               
      DO 20 J = -NAVS, NAVS                                                     
        TOTAV = TOTAV + AVS(IABS(J))                                            
   20 CONTINUE                                                                  
c
c jdemod - remove useless print out of data
c
c      WRITE (6,'('' ASTART,AEND,IA1,IA2,WMIN,TOTAV,IFLAG='',                    
c     >  2G10.3,2I6,2G10.3,I3)') ASTART,AEND,IA1,IA2,WMIN,TOTAV,IFLAG            
c      WRITE (6,'('' WS'',/,(1X,8F9.5))') (WS(IA),IA=IA1,IA2)                    
c      WRITE (6,'('' AENDS'',/,(1X,8F9.5))') (AENDS(IA),IA=IA1,IA2)              
c      WRITE (6,'('' AS'',/,(1X,8F9.5))') (AS(IA),IA=IA1,IA2)                    
c
      NBBS = 0                                                                  
      DO 22 IB = 1, NBS                                                         
        IF (IGS(IB).GT.0) NBBS = NBBS + 1                                       
   22 CONTINUE                                                                  
C                                                                               
C---- CALCULATE VALUES TO BE PLOTTED AND STORE IN CS ARRAY                      
C---- FOR SMOOTHING, CALL HARWELL ROUTINE TO CREATE CUBIC SPLINE                
C---- OR USE WEIGHTED AVERAGING TECHNIQUE  (ADDED DEC 88)                       
C---- FOR WEIGHTED AVERAGES, NOTE THE MESSY FIX IN CASE WS=0  (THIS             
C---- REPRESENTS THE CASE WHERE A DUMMY POINT IS ADDED AT Y=0 BY TAKING         
C---- THE AVERAGE OF THE TWO VALUES ON EITHER SIDE OF Y=0)                      
C---- FIRST CHECK THERE ARE SOME NON-ZERO VALUES - DON'T SMOOTH IF ALL          
C---- VALUES ARE ZERO.                                                          
C---- WEIGHTS ARRAY FOR SMOOTHING IS SIMPLY WS, THE BIN WIDTHS                  
C---- OTHERWISE, JUST COPY APPROPRIATE VALUES OVER TO CS ARRAY                  
C                                                                               
      SMOOTH = ' '                                                              
      DO 50 IB = 1, NBS                                                         
        IF (IGS(IB).LE.0) GOTO 50                                               
        KOUNT = 0                                                               
        DO 25 IA = IA3, IA4                                                     
          IF (BS(IA,IB).EQ.0.0) KOUNT = KOUNT + 1                               
   25   CONTINUE                                                                
        BLABS(IB)(31:31) = ' '                                                  

c
c jdemod - remove additional useless print out       
c
c        WRITE (6,'(1X,A32)') BLABS(IB)(5:36)                                    
c        WRITE (6,'('' BS'',/,1P,(1X,8E9.2))') (BS(IA,IB),IA=IA1,IA2)            
C                                                                               
        IF (IB.GE.ISMOTH .AND. KOUNT.LT.IA4-IA3+1) THEN                         
C                                                                               
          IF (ITEC.EQ.0) THEN                                                   
C         ===================                                                   
            NKNOTS = IA4-IA3+1                                                  
            CALL VC03A  (IA4-IA3+1,NKNOTS,AS(IA3),BS(IA3,IB),WS(IA3),           
     >                                    RD,XN,FN,GN,DN,THETA,0,WORKS)         
            CS(IA3,IB) = TG01B (-1,NKNOTS,XN,FN,GN,AS(IA3))                     
            IF (IA4.GT.IA3) THEN                                                
              DO 30 IA = IA3+1, IA4                                             
                CS(IA,IB) = TG01B (1,NKNOTS,XN,FN,GN,AS(IA))                    
   30         CONTINUE                                                          
            ENDIF                                                               
            SMOOTH = 'SMOOTHED BY CUBIC SPLINE        *'                        
C                                                                               
          ELSEIF (ITEC.EQ.1) THEN                                               
C         =======================                                               
            DO 35 IA = IA3, IA4                                                 
              CS(IA,IB) = 0.0                                                   
              DO 35 J = -NAVS, NAVS                                             
                JA = IPOS (AS(IA)+REAL(J)*WMIN, AENDS, NAS-1)                   
                IF (WS(JA).LE.0.0) JA = JA - 1                                  
                CS(IA,IB) = CS(IA,IB) + BS(JA,IB) * AVS(IABS(J)) / TOTAV        
   35       CONTINUE                                                            
            DO 36 IA = IA3, IA4                                                 
              IF (WS(IA).LE.0.0) CS(IA,IB) = 0.5 *                              
     >                                       (CS(IA-1,IB) + CS(IA+1,IB))        
   36       CONTINUE                                                            
            WRITE (SMOOTH,'(''SMOOTHED,'',I3,'' WEIGHTS, D'',A1,                
     >        F7.4,''  *'')') 2*NAVS+1,AAXLAB(4:4),WMIN                         
C
          ELSEIF (ITEC.EQ.2) THEN
C         ======================= 
C
C         SUM AVERAGE SMOOTHING - ADD CONSECUTIVE POINTS AND THEN
C         AVERAGE TO OBTAIN THE VALUE FOR THE MIDPOINT. SIMILAR TO 
C         WEIGHTED AVERAGE SMOOTHING WITH A WEIGHT OF 1.0 AND WIDTH OF 
C         AN ODD NUMBER. 
C 
       
            range = int(numsmooth/2)
            write(6,*) 'n,r',numsmooth,range
c
c           deal with the bulk of the range of values
c           
            do 510 ia = ia1+range,ia2-range
              cs(ia,ib) = 0.0
              do 520 j = ia-range,ia+range
                cs(ia,ib) = cs(ia,ib) + bs(j,ib)
520           continue
              cs(ia,ib) = cs(ia,ib)/numsmooth
510         continue
c
c           deal with the end points
c           for the missing end points simply weight the point of 
c           interest heavier - by the missing number of bins
c
c           bottom end
c
            icount = range
            do 530 ia = ia1,ia1+range-1
              cs(ia,ib) = 0.0
              do 540 j = ia-range+icount,ia+range
                cs(ia,ib) = cs(ia,ib) + bs(j,ib)
540           continue
              cs(ia,ib) = (cs(ia,ib) + icount * bs(ia,ib))/numsmooth
              icount = icount -1
530         continue
c
c           top end
c
            icount = range
            do 550 ia = ia2,ia2-range,-1
              cs(ia,ib) = 0.0
              do 560 j = ia-range,ia+range-icount
                cs(ia,ib) = cs(ia,ib) + bs(j,ib)
560           continue               
              cs(ia,ib) = (cs(ia,ib) + icount*bs(ia,ib))/numsmooth
              icount = icount -1
550         continue
c
             write(smooth,'(i4,'' Num. Average Smoothing'')') numsmooth
c
          ENDIF                                                                 
C                                                                               
          BLABS(IB)(31:31) = '*'                                                
          WRITE (6,'('' CS'',/,1P,(1X,8E9.2))') (CS(IA,IB),IA=IA1,IA2)          
        ENDIF                                                                   
C                                                                               
        AREAS(IB) = 0.0                                                         
        DO 45 IA = 1, NAS                                                       
          IF (IDRAW.EQ.3 .AND. (AS(IA).LT.0.0.OR.AS(IA).GT.0.4)) GOTO 45        
          AREAS(IB) = AREAS(IB) + CS(IA,IB) * WS(IA)                            
   45   CONTINUE                                                                
        IF (IFLAG.EQ.7.AND.IDRAW.EQ.2) AREAS(IB) = 0.5 * AREAS(IB)              
   50 CONTINUE                                                                  
c
C-----------------------------------------------------------------------        
c     Print out plot data in spreadsheet friendly format 
C-----------------------------------------------------------------------        
c
      write(outunit,*)
      write(outunit,'(a)') 'PLOT DATA:',REF
      write(outunit,'(a,2(1x,g15.6))') 'PLOT RANGE:',astart,aend
c      write(outunit,'(a,5i6)') 'IA VALUES:',ia1,ia2,ia3,ia4
      
      write(outunit,'(31(1x,a12))') 'AXIS:',(blabs(ib)(1:12),ib=1,nbs)
c
      do ia = ia1,ia2
         write(outunit,'(30(1x,g12.5))')
     >       as(ia),(cs(ia,ib),ib=1,nbs)
      end do 
      write(outunit,*)

C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE NO PLOT POSSIBLE                                               
C-----------------------------------------------------------------------        
C                                                                               
      IF (IDRAW.EQ.0 .OR. ASTART.EQ.AEND .OR. IA2.LT.IA1) THEN                  
        GOTO 9999                                                               
C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE AN ORDINARY, UNNORMALISED PLOT IS REQUIRED                     
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (IDRAW.EQ.1) THEN                                                  
C                                                                               
C------ CALCULATE MIN AND MAX FUNCTION VALUES WITHIN THE A RANGE.               
C------ DRAW FRAMEWORK FOR GRAPH.                                               
C------ PLOT RESULTS, TRUNCATED TO REGION BSTART TO BEND ...                    
C                                                                               
        BMIN = BEND                                                             
        BMAX = BSTART                                                           
        DO 110 IB = 1, NBS                                                      
          DO 100 IA = IA1, IA2                                                  
            BMAX = MAX (BMAX, CS(IA,IB))                                        
            BMIN = MIN (BMIN, CS(IA,IB))                                        
  100     CONTINUE                                                              
  110   CONTINUE                                                                
        BMIN = MAX (BMIN, BSTART)                                               
        BMAX = MIN (BMAX, BEND)                                                 
        CALL LIM_GRTSET (TITLE,REF,VIEW,PLANE,JOB,ASTART,AEND,
     >               BMIN,BMAX,            
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)          
        DO 120 IB = 1, NBS                                                      
          IF (IGS(IB).GT.0)                                                     
     >      CALL LIM_GRTRAC (AS(IA1),CS(IA1,IB),IA2-IA1+1,BLABS(IB),
     >                       'LINE')         
  120   CONTINUE                                                                
C                                                                               
C-----------------------------------------------------------------------        
C     CASE WHERE A NORMALISED PLOT IS REQUIRED                                  
C-----------------------------------------------------------------------        
C                                                                               
      ELSEIF (IDRAW.GE.2) THEN                                                  
C                                                                               
C------ CALCULATE MIN AND MAX FUNCTION VALUES WITHIN THE A RANGE.               
C------ VALUES LESS THAN 0 WILL NOT BE PLOTTED  (GRAPH FROM 0 TO 1)             
C------ DRAW FRAMEWORK FOR GRAPH.                                               
C------ SCALE VALUES AND PLOT RESULTS                                           
C                                                                               
        DO 210 IB = 1, NBS                                                      
          BMAX = 0.0                                                            
          DO 200 IA = IA1, IA2                                                  
            BMAX = MAX (BMAX, CS(IA,IB))                                        
  200     CONTINUE                                                              
          FACTS(IB) = BMAX                                                      
  210   CONTINUE                                                                
        CALL LIM_GRTSET (TITLE,REF,VIEW,PLANE,JOB,ASTART,AEND,0.0,1.0,              
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)          
        DO 230 IB = 1, NBS                                                      
          DO 220 IA = IA1, IA2                                                  
            NORMS(IA) = 0.0                                                     
            IF (FACTS(IB).NE.0.0) NORMS(IA) = CS(IA,IB) / FACTS(IB)             
  220     CONTINUE                                                              
          IF (AREAS(IB).NE.0.0)                                                 
     >      WRITE (BLABS(IB)(12:20),'(1P,E9.2)') AREAS(IB)                      
          IF (FACTS(IB).NE.0.0)                                                 
     >      WRITE (BLABS(IB)(22:30),'(1P,E9.2)') FACTS(IB)                      
          IF (IGS(IB).GT.0)                                                     
     >      CALL LIM_GRTRAC (AS(IA1),NORMS(IA1),IA2-IA1+1,BLABS(IB),
     >                       'LINE')         
  230   CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C-----------------------------------------------------------------------        
C     ALL CASES: FINISH GHOST PICTURE IF REQUIRED                               
C-----------------------------------------------------------------------        
C                                                                               
      IF (IFLAG.EQ.2.OR.IFLAG.EQ.3.OR.IFLAG.EQ.5.OR.IFLAG.EQ.6.OR.              
     >    IFLAG.EQ.7) CALL FRAME                                                
 9999 RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LIM_GRTSET (TITLE,REF,VIEW,PLANE,JOB,XMIN,XMAX,                    
     >    YMIN,YMAX,TABLE,XLABEL,YLABEL,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)           
      REAL      XMIN,XMAX,YMIN,YMAX                                             
      INTEGER   IFLAG,IDRAW,NBBS                                                
      CHARACTER YLABEL*24,XLABEL*24,REF*36,VIEW*72,PLANE*36,TABLE*36            
      CHARACTER TITLE*80,JOB*72,SMOOTH*72,ANLY*36                               
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  * GRTSET:  DRAWS FRAMES, TITLE, X AND Y AXES AND ASSOCIATED LABELS  *        
C  *          FOR A NEW PICTURE AND INITIALISES COMMON BLOCK LIM_COMGRA*        
C  *          FOR ROUTINE GRTRAC.                                      *        
C  *                                                                   *        
C  * ARGUMENTS :-                                                      *        
C  *                                                                   *        
C  * TITLE  - TITLE OF PICTURE                                         *        
C  * REF    - STRING OF REFERENCE INFO TO BE DISPLAYED                 *        
C  * VIEW   - EXTRA STRING OF DATA TO BE DISPLAYED                     *        
C  * PLANE  - DITTO                                                    *        
C  * TABLE  - NORMALLY "SYMBOL TABLE"                                  *        
C  * XMIN   - MINIMUM X VALUE                                          *        
C  * XMAX   - MAXIMUM X VALUE                                          *        
C  * YMIN   - MINIMUM Y VALUE                                          *        
C  * YMAX   - MAXIMUM Y VALUE                                          *        
C  * XLABEL - X-AXIS LABEL                                             *        
C  * YLABEL - Y-AXIS LABEL                                             *        
C  * IFLAG  - PASSED FROM DRAW, 3 MEANS DONT DRAW AXES ON GRAPH        *        
C  * SMOOTH - SOME RESULTS ARE SMOOTHED, PRINT MESSAGE IN SYMBOL TABLE.*        
C  * NBBS   - Number of plots                                          *
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'gcom1'

      COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
     >                    NPLOTS,ISPOT          
      REAL            CXMIN,CXMAX,CYMIN,CYMAX                                   
      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT                                  
C                                                                               
      IF (IFLAG.EQ.4.OR.IFLAG.EQ.5) RETURN                                      
C     ====================================                                      
C                                                                               
      CXMIN  = XMIN                                                             
      CXMAX  = XMAX                                                             
      CYMIN  = YMIN                                                             
      CYMAX  = YMAX                                                             
      IPLOTS = 0                                                                
      ICOL   = 1                                                                
      NPLOTS = NBBS                                                             
      ISPOT  = 12                                                               
      IF (NPLOTS.GT.10) ISPOT = 10                                              
      IF (NPLOTS.GT.15) ISPOT = 8                                              
C                                                                               
C---- DRAW TITLES                                                               
C                                                                               
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)                                         
      CALL MAP    (0.0, 1.35, 0.0, 1.0)                                         
      CALL CTRMAG (20)                                                         
      CALL LINCOL (1)                                                          
      CALL THICK  (2) 
      L = LENSTR(TITLE)                                                         
      CALL PCSCEN (0.8, 0.95, TITLE(:L))
      L = LENSTR(TABLE)                                           
      CALL PCSCEN (1.18, 0.87, TABLE(:L))                                          
      L = LENSTR (XLABEL)                                                       
      IF (IFLAG.NE.3) CALL PCSCEN (0.5,0.025,XLABEL(:L))                       
      CALL THICK  (1)                                                          
      CALL CTRMAG (10)                                                         
      L = LENSTR (REF)                                                          
      CALL PCSCEN (1.16, 0.12, REF(:L))                                        
      CALL CTRMAG (10)                                                         
      CALL PCSCEN (1.20, 0.17, JOB(37:72))                                     
      CALL PCSCEN (1.20, 0.22, JOB( 1:36))                                     
      CALL CTRMAG (10)                                                         
      L = LENSTR(VIEW)
      CALL PCSCEN (1.18, 0.27, VIEW(:L))
      L = LENSTR(PLANE)                                           
      CALL PCSCEN (1.18, 0.31, PLANE(:L))                                          
      L = LENSTR(ANLY)
      CALL PCSCEN (1.18, 0.35, ANLY(:L))                                           
      CALL PCSCEN (1.18, 0.39, SMOOTH(37:72))                                  
      CALL PCSCEN (1.18, 0.43, SMOOTH( 1:36))                                  
      IF     (IDRAW.EQ.2) THEN                                                 
        CALL PLOTST (1.05 ,0.818, '        TOT AREA   SCALE        ')          
      ELSEIF (IDRAW.EQ.3) THEN                                                 
        CALL PLOTST (1.05 ,0.818, '      0:0.4 AREA   SCALE        ')          
      ENDIF                                                                    
C                                                                               
C---- DRAW FRAMES                                                               
C                                                                               
C      THE FOLLOWING USES GHOST TO DRAW FRAMES AROUND THE GRAPHS
C      COMMENT THIS OUT FOR NOW. SINCE CURRENTLY WE ARE NOT USING GHOST 
C      AND THESE WOULD ONLY GET IN THE WAY.
C 
C      DAVID ELDER, 1989 NOV 21 
C
C
      CALL LINCOL (3)                                                         
      CALL POSITN (0.1, 0.1)                                                  
      CALL   JOIN (0.1, 0.9)                                                  
      CALL   JOIN (0.9, 0.9)                                                  
      CALL   JOIN (0.9, 0.1)                                                  
      CALL   JOIN (0.1, 0.1)                                                 
      CALL POSITN (0.93, 0.1)                                                 
      CALL   JOIN (0.93, 0.9)                                                 
      CALL   JOIN (1.35, 0.9)                                                 
      CALL   JOIN (1.35, 0.1)                                                 
      CALL   JOIN (0.93, 0.1)                                                 
      CALL POSITN (0.93, 0.85)                                                 
      CALL   JOIN (1.35, 0.85)                                                
      CALL POSITN (0.93, 0.15)                                                
      CALL   JOIN (1.35, 0.15)                                                
      CALL POSITN (0.93, 0.20)                                                
      CALL   JOIN (1.35, 0.20)                                                
      CALL POSITN (0.93, 0.25)                                                
      CALL   JOIN (1.35, 0.25)                                                
C                                                                               
C---- FOR 3D CASE IGNORE AXES CODE                                              
C---- DRAW X AXIS                                                               
C                                                                               
      IF (IFLAG.EQ.3) RETURN                                                    
      ITEN = IEXP(XMIN, XMAX)                                                   
      IF (ITEN .NE. 0) THEN                                                     
        CALL POSITN (0.8,0.04)                                                 
        CALL CTRMAG (14)                                                       
        CALL TYPECS ('X10')                                                    
        CALL CTRMAG (12)                                                       
        CALL TYPENI ((ITEN))                                                   
      ENDIF                                                                     
      POWER = 10.0**(-ITEN)                                                     
      TMIN = XMIN * POWER                                                       
      TMAX = XMAX * POWER                                                       
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)                                         
      CALL MAP    (TMIN, TMAX, 0.1, 0.9)                                       
      CALL CTRMAG (10)                                                         
      CALL XSCALE                                                              
C                                                                               
C---- DRAW Y AXIS AND LABELS                                                    
C                                                                               
      IF (YMAX .EQ. YMIN) THEN                                                  
         WRITE(6, *)  ' ERROR - ROUTINE GRTSET '                                
         WRITE(6, *)  ' SCALE ',YLABEL,' FMIN = FMAX ',YMIN                     
      ENDIF                                                                     
      ITEN = IEXP(YMIN, YMAX)                                                   
      POWER = 10.0**(-ITEN)                                                     
      TMIN = YMIN * POWER                                                       
      TMAX = YMAX * POWER                                                       
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)                                       
      CALL MAP    (0.0, 1.0, TMIN, TMAX)                                       
      CALL LINCOL (3)                                                          
      CALL CTRMAG (10)                                                         
      CALL YSCALE                                                              
      CALL PSPACE (0.0, 1.35, 0.11, 0.89)                                       
      CALL MAP    (0.0, 1.35, 0.0, 1.0)                                         
      CALL CTRORI (90.0)                                                       
      IF (ITEN .NE. 0) THEN                                                     
         CALL POSITN (0.02, 0.01)                                              
         CALL CTRMAG (14)                                                      
         CALL TYPECS ('X10')                                                   
         CALL CTRMAG (12)                                                      
         CALL TYPENI ((ITEN))                                                  
      ENDIF                                                                     
      CALL LINCOL (1)                                                          
      CALL CTRMAG (14)                                                         
      CALL THICK  (2)                                                          
      L = LENSTR (YLABEL)                                                       
      CALL PCSEND (0.02,0.99,YLABEL(:L))                                       
      CALL CTRORI (0.0)                                                        
      CALL THICK  (1)                                                          
C                                                                               
C---- DRAW LINE FOR FUNCTION=0 IN CASES WHERE IFLAG=6  (EG NET EROSION)         
C                                                                               
      IF (IFLAG.EQ.6) THEN                                                     
        CALL LINCOL (3)                                                        
        CALL PSPACE (0.1, 0.9, 0.11, 0.89)                                     
        CALL MAP    (CXMIN,CXMAX,CYMIN,CYMAX)                                  
        CALL POSITN (CXMIN, 0.0)                                               
        CALL JOIN   (CXMAX, 0.0)                                               
      ENDIF           

 1235 FORMAT (A6,A74)
 1236 FORMAT (A6,A36)
 1237 FORMAT (A6,A24)
 1238 FORMAT (A6,I7)
                                                         
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LIM_GRTRAC (X ,Y ,NPTS ,NAME, CURVE)                               
      CHARACTER NAME*36,CURVE*(*)                                               
      INTEGER   NPTS                                                            
      REAL      X(NPTS),Y(NPTS)                                                 
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  * GRTRAC: DRAWS GRAPHS OF THE DATA POINTS GIVEN BY X AND Y.         *        
C  *         N.B. AN INITIALISING CALL TO GRTSET MUST BE MADE          *        
C  *         BEFORE USING THIS SUBROUTINE.                             *        
C  *                                                                   *        
C  * ARGUMENTS :-                                                      *        
C  *                                                                   *        
C  *  X     - X DATA VALUES                                            *        
C  *  Y     - Y DATA VALUES                                            *        
C  *  NPTS  - TOTAL NUMBER OF DATA POINTS                              *        
C  *  NAME  - IDENTIFYING LABEL TO GO IN SYMBOL TABLE                  *        
C  *                                                                   *        
C  *********************************************************************        
C                    
      INCLUDE 'gcom1'
                                                           
      COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
     >                    NPLOTS,ISPOT          
      REAL            CXMIN,CXMAX,CYMIN,CYMAX                                   
      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT                                  
C                                                                               
      INTEGER COLOUR(8)                                                         
      DATA COLOUR /2,4,6,5,7,3,6,8/                                             
C                                                                               
C     WRITE (6,'('' GRTRAC: X ='',/,(1X,8F9.5))')    (X(I),I=1,NPTS)            
C     WRITE (6,'(''     AND Y ='',/,1P,(1X,8E9.2))') (Y(I),I=1,NPTS)            
C                                                                               
C---- INCREMENT COUNTS, SET COLOUR AND LINE PATTERN                             
C---- ARGUMENTS TO "BROKEN" GIVE SIZES OF  (DASH1, GAP1, DASH2, GAP2)           
C                                                                               
      IBROK  = IPLOTS                                                           
      IPLOTS = IPLOTS + 1                                                       
      CALL LINCOL (COLOUR(ICOL))                                               
      ICOL   = ICOL + 1                                                         
      IF (ICOL.GT.8) ICOL = 1                                                   
      IF (IPLOTS.LE.1) THEN                                                    
         CALL FULL                                                             
      ELSEIF (NPLOTS.LE.5) THEN                                                
         CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)                         
      ELSE                                                                     
         CALL BROKEN (2*IBROK,1*IBROK,2*IBROK,1*IBROK)                         
      ENDIF                                                                    
C                                                                               
C---- DRAW PLOT                                                                 
C                                                                               
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)                                       
      CALL MAP    (CXMIN,CXMAX,CYMIN,CYMAX)                                    
C                                                                               
      IF (CURVE.EQ.'POINT') THEN                                                
        DO 10 I = 1, NPTS                                                       
          IF (X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.Y(I).GE.CYMIN.AND.            
     >        Y(I).LE.CYMAX) call point ( X(I),Y(I) )                     
   10   CONTINUE                                                                
      ELSE                                                                      
        CALL POSITN (X(1), Y(1))                                               
        DO 20 I = 2, NPTS                                                       
          CALL JOIN (X(I), Y(I))                                               
          IF ((I.EQ.3.OR.I.EQ.NPTS-2).and.(
     >        X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.
     >        Y(I).GE.CYMIN.AND.Y(I).LE.CYMAX            
     >        )) CALL PLOTST (X(I),Y(I),NAME(1:4))         
 20     CONTINUE                                                                
      ENDIF                                                                     
C                                                                               
C---- WRITE ENTRY IN SYMBOL TABLE                                               
C                                                                               
      CALL CTRMAG (ISPOT)                                                      
      CALL PSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL CSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL MAP    (0.0, 1.35, 0.0,1.0)                                          
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0                               
      IF (CURVE.EQ.'POINT') THEN                                               
        CALL PLOTST (0.96,SPOT,'+')                                            
        CALL PLOTST (0.99,SPOT,'+')                                            
        CALL PLOTST (1.02,SPOT,'+')                                            
      ELSE                                                                     
        CALL POSITN (0.95,SPOT)                                                
        CALL JOIN   (1.03,SPOT)                                                
      ENDIF                                                                    
      CALL FULL                                                                
      CALL PLOTST (1.05 ,SPOT, NAME(5:36))                                     
 1236 FORMAT(A6,A36)
 1250 FORMAT(2G20.8)
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LIM_GR3D (SURFAS,NPTS,NAME,IVEW3D,PROJ3D,IBAS3D,                   
     >                 SUREDG,LIMEDG)                                           
      INTEGER  IBOX,NPTS,IVEW3D,IBAS3D,LIMEDG                                   
      REAL     SURFAS(192,192),PROJ3D,SUREDG(192,192)                           
      CHARACTER*36 NAME                                                         
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  GR3D  THIS ROUTINE DRAWS GRAPHS FROM A 3D VANTAGE POINT.  THE    *        
C  *  VALUES IN SURFAS ARE ASSUMED TO LIE ON A REGULAR GRID, NPTS*NPTS,*        
C  *  AND A FEW PARAMETERS LIKE "IVEW3D" ALLOW THE VIEWPOINT TO BE     *        
C  *  CHANGED.  AN INITIALISING CALL TO GRTSET WITH LAST ARG = 3       *        
C  *  SHOULD BE MADE BEFORE CALLING THIS ROUTINE.                      *        
C  *                                                                   *        
C  *  C.M.FARRELL   FEBRUARY 1988                                      *        
C  *                                                                   *        
C  *********************************************************************        
C                    
      INCLUDE 'gcom1' 
                                                           
      COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
     >                    NPLOTS,ISPOT          
      REAL            CXMIN,CXMAX,CYMIN,CYMAX                                   
      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT                                  
C                                                                               
C---- SURCOL: TOP COLOUR, UNDERSIDE COLOUR, BASE COLOUR                         
C---- (1:8 REPRESENT BLACK,RED,GREEN,BLUE,WHITE,CYAN,MAGENTA,YELLOW)            
C---- SURDIR: 0 TO 3 GIVES 0 DEGREES, 90,180,270 DEGREES VIEW                   
C---- SURANG: 35.2644 ISOMETRIC, OTHERWISE IN RANGE -90 TO 90 DEGREES           
C---- SURBAS: 0/1 UNDERSIDE OFF/ON, 0/1/2 BASE OFF/ON/ON AT LAST ARG            
C---- SURPLT: DIM ARRAY(192,192), BUT ONLY PLOT ARRAY(1:NPTS,1:NPTS)            
C                                                                               
      IPLOTS = IPLOTS + 1                                                       
      CALL PSPACE (0.105, 0.895, 0.11, 0.89)                                   
      CALL SURCOL (2,7,2)                                                      
      CALL SURDIR (IVEW3D)                                                     
      CALL SURANG (PROJ3D)                                                     
      CALL SURBAS (0,IBAS3D,0.0)                                               
      CALL SURPLT (SURFAS,1,NPTS,192,1,NPTS,192)                                
      IF (LIMEDG.EQ.1) THEN                                                     
        CALL SURCOL (4,6,2)                                                    
        CALL SURBAS (0,0,0.0)                                                  
        CALL SURPLT (SUREDG,1,NPTS,192,1,NPTS,192)                              
      ENDIF                                                                     
C                                                                               
C---- WRITE ENTRY IN SYMBOL TABLE                                               
C                                                                               
      CALL LINCOL (2)                                                          
      CALL CTRMAG (ISPOT)                                                      
      CALL PSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL CSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL MAP    (0.0, 1.35, 0.0,1.0)                                          
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0                               
      SPOT = 0.818 - 0.03 * IPLOTS                                             
      CALL POSITN (0.95,SPOT)                                                  
      CALL JOIN   (1.03,SPOT)                                                  
      CALL FULL                                                                
      CALL PLOTST (1.05 ,SPOT, NAME)                                           
 1236 FORMAT (A6,A36)
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LIM_GRM (SURFAS,NPTS,MPTS,IPLOT,IPLANE,NAME,
     >               IVEW3D,PROJ3D,IBAS3D,COORD1,COORD2)
      INTEGER  NPTS,MPTS,IPLOT,IPLANE,IVEW3D,IBAS3D                      
      REAL     SURFAS(192,192),PROJ3D,COORD1(90),COORD2(90)
      CHARACTER*(*) NAME                                                      
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  GRM   THIS ROUTINE DRAWS GRAPHS FROM A 3D VANTAGE POINT OR A     *
C  *  CONTOUR PLOT OF THE SAME DATA, DEPENDING ON WHICH PLOT OPTION    *
C  *  CHOSEN. FOR A MESH PLOT THE                                      *        
C  *  VALUES IN SURFAS ARE ASSUMED TO LIE ON A REGULAR GRID, NPTS*MPTS,*        
C  *  AND A FEW PARAMETERS LIKE "IVEW3D" ALLOW THE VIEWPOINT TO BE     *        
C  *  CHANGED.  AN INITIALISING CALL TO GRTSET WITH LAST ARG = 3       *        
C  *  SHOULD BE MADE BEFORE CALLING THIS ROUTINE.                      *        
C  *  FOR A CONTOUR PLOT THE POINTS ARE IRREGULARLY SPACED AND THE     *
C  *  SPACING MUST ALSO BE PASSED TO THE PLOTTING ROUTINE.             *
C  *                                                                   *        
C  *  C.M.FARRELL   FEBRUARY 1988                                      *        
C  *  D. ELDER      MARCH 15 1990                                      *
C  *                                                                   *        
C  *********************************************************************        
C                    
      INCLUDE 'gcom1' 
                                                           
      INTEGER   I,J 

C     COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
C     >                    NPLOTS,ISPOT          
C     REAL            CXMIN,CXMAX,CYMIN,CYMAX                                   
C     INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT                                  
C                                                                               
C---- SURCOL: TOP COLOUR, UNDERSIDE COLOUR, BASE COLOUR                         
C---- (1:8 REPRESENT BLACK,RED,GREEN,BLUE,WHITE,CYAN,MAGENTA,YELLOW)            
C---- SURDIR: 0 TO 3 GIVES 0 DEGREES, 90,180,270 DEGREES VIEW                   
C---- SURANG: 35.2644 ISOMETRIC, OTHERWISE IN RANGE -90 TO 90 DEGREES           
C---- SURBAS: 0/1 UNDERSIDE OFF/ON, 0/1/2 BASE OFF/ON/ON AT LAST ARG            
C---- SURPLT: DIM ARRAY(192,192), BUT ONLY PLOT ARRAY(1:NPTS,1:NPTS)            
C                                                                               
      IPLOTS = IPLOTS + 1                                                    
      CALL PSPACE (0.105, 0.895, 0.11, 0.89)                                   
      CALL SURCOL (2,7,2)                                                      
      CALL SURDIR (IVEW3D)                                                     
      CALL SURANG (PROJ3D)                                                     
      CALL SURBAS (0,IBAS3D,0.0)                                               

C      IF (IPLOT.EQ.2) THEN  
C        WRITE(GNOUT,'(A6,I7,I7)') 'SURFA:',NPTS,MPTS
C        WRITE(GNOUT,'(10G13.5)') ((SURFAS(I,J),I=1,NPTS),J=1,MPTS)
C        WRITE(GNOUT,'(A6)') 'XGRID:'
C        WRITE(GNOUT,'(10G13.5)') (COORD1(I),I=1,NPTS)
C        WRITE(GNOUT,'(A6)') 'YGRID:'
C        WRITE(GNOUT,'(10G13.5)') (COORD2(I),I=1,MPTS)               
C      ELSE
C        WRITE(GNOUT,'(A6,I7,I7,I7)') 'CONTR:',NPTS,MPTS,10
C        WRITE(GNOUT,'(10G13.5)') ((SURFAS(I,J),I=1,NPTS),J=1,MPTS)
C        WRITE(GNOUT,'(A6)') 'XGRID:'
C        WRITE(GNOUT,'(10G13.5)') (COORD1(I),I=1,NPTS)
C        WRITE(GNOUT,'(A6)') 'YGRID:'
C        WRITE(GNOUT,'(10G13.5)') (COORD2(I),I=1,MPTS)               
C      ENDIF
      IF (LIMEDG.EQ.1) THEN                                                 
        CALL SURCOL (4,6,2)                                                    
        CALL SURBAS (0,0,0.0)                                                  
        CALL SURPLT (SUREDG,1,NPTS,192,1,NPTS,192)                           
C        WRITE(GNOUT,1236) 'NAME0:','LIMITER EDGE'     
      ENDIF                                                                   
C                                                                               
C---- WRITE ENTRY IN SYMBOL TABLE                                               
C                                                                               
      CALL LINCOL (2)                                                          
      CALL CTRMAG (ISPOT)                                                      
      CALL PSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL CSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL MAP    (0.0, 1.35, 0.0,1.0)                                          
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0                               
      SPOT = 0.818 - 0.03 * IPLOTS                                             
      CALL POSITN (0.95,SPOT)                                                  
      CALL JOIN   (1.03,SPOT)                                                  
      CALL FULL                                                                
      CALL PLOTST (1.05 ,SPOT, NAME)                                           
 1236 FORMAT (A6,A36)
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
      SUBROUTINE LIM_GRCONT (VALS,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,                   
     >                   MAXNYS,CLEVEL,XOUTS,YOUTS,NAME)                        
      INTEGER  IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS                            
      REAL     CLEVEL,XOUTS(MAXNXS),YOUTS(2*MAXNYS)                             
      REAL     VALS(MAXNXS,2*MAXNYS)                                            
      CHARACTER*36 NAME                                                         
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  GRCONT  THIS ROUTINE DRAW A CONTOUR PLOT.                        *        
C  *          THE ARRAYS VALS AND YOUTS ARE PECULIARLY DIMENSIONED FOR *        
C  *          PASSING TO ROUTINE CONTIL; THE Y RANGE IS SPLIT INTO TWO *        
C  *          IN CASE IT MIGHT EXCEED THE GHOST MAXIMUM OF 192.        *        
C  *                                                                   *        
C  *  C.M.FARRELL   MARCH 1988                                         *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      INCLUDE 'gcom1'

      COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
     >                    NPLOTS,ISPOT          
      REAL            CXMIN,CXMAX,CYMIN,CYMAX                                   
      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT                                  
C                                                                               
      INTEGER COLOUR(8)                                                         
      DATA COLOUR /2,4,6,5,7,3,6,8/                                             
      WRITE (6,'('' GRCONT: IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS'',            
     >  /7X,6I7)') IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS                        
C                                                                               
C---- INCREMENT COUNTS, SET COLOUR AND LINE PATTERN                             
C                                                                               
      IBROK  = IPLOTS                                                           
      IPLOTS = IPLOTS + 1                                                       
      CALL LINCOL (COLOUR(ICOL))                                               
      ICOL   = ICOL + 1                                                         
      IF (ICOL.GT.8) ICOL = 1                                                   
      IF (IPLOTS.LE.1) THEN                                                    
         CALL FULL                                                             
      ELSEIF (NPLOTS.LE.5) THEN                                                
         CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)                         
      ELSE                                                                     
         CALL BROKEN (2*IBROK,1*IBROK,2*IBROK,1*IBROK)                         
      ENDIF                                                                    
C                                                                               
C---- DRAW CONTOUR PLOT  (DISABLE ANNOTATION)                                   
C---- SPLIT UP Y RANGE SINCE LIMIT OF 192 APPLIES, AND Y RANGE IS 200.          
C---- (USE 1/3 : 2/3 SPLIT TO ENSURE REGION NEAR Y=0 ISN'T AFFECTED).           
C                                                                               
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)                                      
      CALL MAP    (CXMIN, CXMAX, CYMIN, CYMAX)                                 
      CALL CONOTA (0)                                                           
      LYMIN = IYMIN + MAXNYS                                                    
      LYMAX = IYMAX + MAXNYS                                                    
      LYMID = LYMIN + (LYMAX - LYMIN) / 3                                       
c slmod begin
c      WRITE(0,*) '   Here in GRCONT'
c
c      CALL CONTIL (VALS,IXMIN,IXMAX,MAXNXS,LYMIN,LYMID,2*MAXNYS,                
c     >             CLEVEL,1,1,XOUTS,YOUTS)                                      
c      CALL CONTIL (VALS,IXMIN,IXMAX,MAXNXS,LYMID,LYMAX,2*MAXNYS,                
c     >             CLEVEL,1,1,XOUTS,YOUTS)                                      
c
      CALL CONTIL (VALS,IXMIN,IXMAX,MAXNXS,IYMIN,IYMID,MAXNYS,                
     >             CLEVEL,1,1,XOUTS,YOUTS)                                      
c slmod end
C                                                                               
C---- WRITE ENTRY IN SYMBOL TABLE                                               
C                                                                               
      CALL CTRMAG (ISPOT)                                                      
      CALL PSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL CSPACE (0.0, 1.35, 0.0,1.0)                                          
      CALL MAP    (0.0, 1.35, 0.0,1.0)                                          
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0                              
      CALL POSITN (0.95,SPOT)                                                  
      CALL JOIN   (1.03,SPOT)                                                  
      CALL FULL                                                                
      CALL PLOTST (1.05 ,SPOT, NAME(5:36))                                     
 1236 FORMAT(A6,A36)
      RETURN                                                                    
      END    
c
c slmod begin
c
      SUBROUTINE LIM_GRCONT95 (VALS,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,
     >                     MAXNYS,CLEVEL,XOUTS,YOUTS,NAME)
      INTEGER  IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS
      REAL     CLEVEL,XOUTS(MAXNXS),YOUTS(MAXNYS)
      REAL     VALS(MAXNXS,MAXNYS)
      CHARACTER*36 NAME
C
C  *********************************************************************
C  *                                                                   *
C  *  GRCONT  THIS ROUTINE DRAW A CONTOUR PLOT.                        *
C  *          THE ARRAYS VALS AND YOUTS ARE PECULIARLY DIMENSIONED FOR *
C  *          PASSING TO ROUTINE CONTIL; THE X,Y RANGES ARE SPLIT IF   *
C  *          THEY EXCEED THE GHOST LIMIT OF 192.                      *
C  *                                                                   *
C  *  C.M.FARRELL   FEBRUARY 1989                                      *
C  *                                                                   *
C  *********************************************************************
C
      COMMON /LIM_COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,
     >                    NPLOTS,ISPOT          

      REAL            CXMIN,CXMAX,CYMIN,CYMAX
      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT
C
      INTEGER COLOUR(8),IXB,IXE,IYB,IYE
      DATA COLOUR /2,4,6,5,7,3,6,8/
      WRITE (6,'(A6,1X,A32)') 'NAME: ',NAME(5:36)
      WRITE (6,'('' GRCONT: IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS'',
     >  /7X,6I7)') IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS
C
C---- INCREMENT COUNTS, SET COLOUR AND LINE PATTERN
C
      IBROK  = IPLOTS
      IPLOTS = IPLOTS + 1
      CALL LINCOL (COLOUR(ICOL))
      ICOL   = ICOL + 1
      IF (ICOL.GT.8) ICOL = 1
      IF (IPLOTS.LE.1) THEN
         CALL FULL
      ELSEIF (NPLOTS.LE.5) THEN
         CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
      ELSE
         CALL BROKEN (2*IBROK,1*IBROK,2*IBROK,1*IBROK)
      ENDIF
c slmod tmp
      CALL FULL
c slmod end
C
C---- DRAW CONTOUR PLOT  (DISABLE ANNOTATION)
C---- SPLIT UP X,Y RANGES IF REQUIRED ...
C
      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
      CALL MAP    (CXMIN, CXMAX, CYMIN, CYMAX)
      CALL CONOTA (0)
      DO 110 IXB = IXMIN, IXMAX, 192
        IXE = MIN (IXMAX, IXB+192)
        DO 100 IYB = IYMIN, IYMAX, 192
          IYE = MIN (IYMAX, IYB+192)
          CALL CONTIL (VALS,IXB,IXE,MAXNXS,IYB,IYE,MAXNYS,
     >                 CLEVEL,1,1,XOUTS,YOUTS)
  100   CONTINUE
  110 CONTINUE
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
      CALL CTRMAG (ISPOT)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0
      CALL POSITN (0.95,SPOT)
      CALL JOIN   (1.03,SPOT)
      CALL FULL
      CALL PLOTST (1.05 ,SPOT, NAME(5:36))
      RETURN
      END
c slmod end
C
      INTEGER FUNCTION IEXP (SMIN, SMAX)                                        
      implicit none 
      REAL SMIN,SMAX,power
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  IEXP: RETURNS SUITABLE EXPONENT FOR A SCALE SMIN TO SMAX.        *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
      POWER = MAX(ABS(SMIN), ABS(SMAX))                                         
      IF (POWER.EQ.0.0) then
        iexp = 0
      ELSEIF (POWER .GE. 100.0) THEN
        IEXP = INT(LOG10(POWER) - 0.3)
      ELSE IF (POWER .LT. 0.5) THEN
        IEXP = INT(LOG10(POWER) - 0.3)
      ELSE
        IEXP = 0
      ENDIF
      RETURN                                                                    
      END                                                                       
C                                                                               
C                                                                               
C                                                                               
c      INTEGER FUNCTION LENSTR (ASTR)                                            
c      CHARACTER*(*) ASTR                                                        
C                                                                               
C  *********************************************************************        
C  *                                                                   *        
C  *  LENSTR: RETURNS EFFECTIVE LENGTH OF STRING ASTR IGNORING         *        
C  *          ANY TRAILING BLANKS.                                     *        
C  *                                                                   *        
C  *********************************************************************        
C                                                                               
c      DO 10 I = LEN(ASTR),1,-1                                                  
c         IF (ASTR(I:I) .NE. ' ') THEN                                           
c            LENSTR = I                                                          
c            RETURN                                                              
c         ENDIF                                                                  
c   10 CONTINUE                                                                  
c      LENSTR = 1                                                                
c      RETURN                                                                    
c      END                                                                       
c
