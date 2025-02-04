C@PROCESS VECTOR(LEV(0)),OPT(1)
C
C
      SUBROUTINE DRAW (AS,WS,BS,MAXNAS,NAS,ANLY,
     >  NBS,ISMOTH,ASTART,AEND,BSTART,BEND,IGS,ITEC,AVS,NAVS,JOB,
     >  TITLE,AAXLAB,BAXLAB,BLABS,REF,VIEW,PLANE,TABLE,IDRAW,IFLAG,AAA,
     >  IEXPT)
      use mod_params
      use mod_slout
      use allocate_arrays
      implicit none
c
c      IMPLICIT LOGICAL (A-Z)
c
      INTEGER   MAXNAS,NAS,IBS,NBS,IDRAW,ISMOTH,IFLAG,IGS(*),ITEC,NAVS
      integer   iexpt  
      REAL      AS(*),WS(*),BS(MAXNAS,*),ASTART,AEND,BSTART,BEND
      REAL      AVS(0:NAVS),AAA
      CHARACTER JOB*72,TITLE*80,TABLE*36
c slmod begin - new
c...  Allow for longer y-axis labels:
      CHARACTER AAXLAB*36,BAXLAB*(*),BLABS(NBS)*(*)
c...  Array size should be set by some parameter, but 1000 ought to
c     be good enough:
      REAL      RV(1000),ZV(1000)
c
c      CHARACTER AAXLAB*36,BAXLAB*36,BLABS(*)*36
c slmod end
      character REF*(*),VIEW*(*),PLANE*(*),anly*(*)
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
c  *        - for ASTART,AEND - factor in experimental data as well    *
c  *          to get AMIN,AMAX                                         *
C  * BSTART - MINIMUM BS VALUE FOR AXIS: ACTUALLY USE MAX (BSTART, BS) *
C  * BEND   - MAXIMUM BS VALUE FOR AXIS: ACTUALLY USE MIN (BEND, BS)   *
c  *        - for BSTART, BEND also factor in experimental data        *
c  *          to get BMIN,BMAX                                         *
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
C  *          IDRAW=5 - DO NOT PRINT AREAS                             *
C  *          IDRAW=7 - unnormalized - MIN SET TO 0.0                  *
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
C  * IEXPT  - INDEX OF EXPERIMENTAL DATA SET TO BE INCLUDED IN THE     *
C  *          PLOT.                                                    *
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
      INTEGER IA,IB,in,NKNOTS,IA1,IA2,IA3,IA4,KOUNT,J,JA
      INTEGER NBBS,IPOS,ierr
C     INCLUDE "PARAMS"
c      include 'params'
c     mxxnbs was too small - must be larger than maximum of ring number
c     and charge state ; Krieger IPP/98
      integer :: mxxnas,mxxnbs
      PARAMETER (MXXNAS=8000)
      !PARAMETER (MXXNAS=8000,MXXNBS=maxizs+2)
      real    amin,amax,bmin,bmax,wmin,totav,tg01b

!      REAL    BMIN,BMAX,FACTS(MXXNBS),NORMS(MXXNAS),CS(MXXNAS,MXXNBS)
!      real    amin,amax
!      REAL    RD(MXXNAS),FN(MXXNAS),AREAS(MXXNBS)
!      REAL    GN(MXXNAS),DN(MXXNAS),THETA(MXXNAS),XN(MXXNAS),TG01B
!      REAL    WORKS(16*MXXNAS+6),WMIN,TOTAV,AENDS(MXXNAS)

      REAL    NORMS(MXXNAS),RD(MXXNAS),FN(MXXNAS)
      REAL    GN(MXXNAS),DN(MXXNAS),THETA(MXXNAS),XN(MXXNAS)
      REAL    WORKS(16*MXXNAS+6),AENDS(MXXNAS)

      !REAL    FACTS(MXXNBS),CS(MXXNAS,MXXNBS),AREAS(MXXNBS)
      REAL,allocatable::  FACTS(:),CS(:,:),AREAS(:)
      

      INTEGER RANGE,ICOUNT
      CHARACTER SMOOTH*72
c
      integer len,lenstr
      external lenstr
c
      LOGICAL NEGTIV
c      COMMON /TRACE/ FACTS,NORMS,AREAS,CS,AENDS,
c     >               RD,XN,FN,GN,DN,THETA,WORKS
c slmod begin
c      include 'slout'

      INTEGER i,i1,i2
c slmod end
C
c      COMMON /NSMOOTH/ NUMSMOOTH,cgrprint
c
      INTEGER NUMSMOOTH,cgrprint
c
c     Local variables for experimental data
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=1)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=1000,dataunit=13,maxcols=1)
c
c      integer axis_type,num_expt,ncols
c      real expt_axis(maxdatx),expt_data(maxdatx)
      real expt_min,expt_max
      integer iexpt_plot 
c      character*100 datatitle
c
c     Allocate local storage
c
      mxxnbs = maxizs + 2
      call allocate_array(facts,mxxnbs,'FACTS: Draw local',ierr)
      call allocate_array(cs,mxxnas,mxxnbs,'CS: Draw local',ierr)
      call allocate_array(areas,mxxnbs,'AREAS: Draw local',ierr)      
C
c     Plot start ...
c
      WRITE (6,'(/,'' DRAW: '',A36,i4,/1X,42(''-''))') REF,nbs
c      write(6,'(a,2(1x,g12.5))')  'ASTART,AEND:',astart,aend
c      write(6,'(a,2(1x,g12.5))')  'BSTART,BEND:',bstart,bend

      IF (NAS.GT.MXXNAS.OR.NBS.GT.MXXNBS) THEN
        WRITE (6,'('' DRAW: ERROR! CHECK PARMS: NAS,MAXNAS,MXXNAS'',3I5,
     >    '',  NBS,MXXNBS'',2I3)') NAS,MAXNAS,MXXNAS,NBS,MXXNBS
        STOP
      ENDIF
C-----------------------------------------------------------------------
c
c     Load and initialize any experimental data required for the plot
c
C-----------------------------------------------------------------------

c     Initialize experimental plot variable. 
c
      iexpt_plot = iexpt
      amin = astart
      amax = aend
c
c     This is called for all values of IEXPT since EXPT_NSETS could be
c     non-zero 
c
      call find_scales(amin,amax,expt_min,expt_max,iexpt_plot)
c
c      write(6,'(a,2(1x,g12.5))')  'AMIN,AMAX:',amin,amax
c      write(6,'(a,2(1x,g12.5))')  'EMIN,EMAX:',expt_min,expt_max
c
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
        IF (AS(IA).LT.Amin) IA1 = IA1 + 1
        IF (AS(IA).LE.Amax  ) IA2 = IA2 + 1
        AENDS(IA) = AS(IA) + 0.5 * WS(IA)
        DO 10 IB = 1, NBS
          CS(IA,IB) = BS(IA,IB)
   10 CONTINUE
c
      IA1 = MAX (1,   IA1-1)
      IA2 = MIN (NAS, IA2+1)
      IA3 = MAX (1,   IA1-5)
      IA4 = MIN (NAS, IA2+5)
c
      WMIN = 1.E10
      DO 15 IA = IA3, IA4
        IF (WS(IA).GT.0.0) WMIN = MIN (WMIN, WS(IA))
   15 CONTINUE
      TOTAV = 0.0
      DO 20 J = -NAVS, NAVS
        TOTAV = TOTAV + AVS(IABS(J))
   20 CONTINUE
C     WRITE (6,'('' AS,AE,IA1/2,WMIN,TOTAV,IFL '',
C    >  2G9.2,2I4,2G9.2,I2)') amin,Amax,IA1,IA2,WMIN,TOTAV,IFLAG
C     WRITE (6,'('' WS'',/,(1X,8F9.5))') (WS(IA),IA=IA1,IA2)
C     WRITE (6,'('' AENDS'',/,(1X,8F9.5))') (AENDS(IA),IA=IA1,IA2)
C     WRITE (6,'('' AS'',/,(1X,8F9.5))') (AS(IA),IA=IA1,IA2)
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
c slmod begin
c        BLABS(IB)(31:31) = ' '
c slmod end
        WRITE (6,'(1X,A32)') BLABS(IB)(5:36)
C       WRITE (6,'('' BS'',/,1P,(1X,8E9.2))') (BS(IA,IB),IA=IA1,IA2)
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
C
        ENDIF
C
        AREAS(IB) = 0.0
        DO 45 IA = 1, NAS
          IF (IDRAW.EQ.3 .AND. (AS(IA).LT.0.0.OR.AS(IA).GT.0.4)) GOTO 45
          AREAS(IB) = AREAS(IB) + CS(IA,IB) * WS(IA) * AAA
   45   CONTINUE
        IF (IFLAG.EQ.7.AND.IDRAW.EQ.2) AREAS(IB) = 0.5 * AREAS(IB)
   50 CONTINUE
C
C-----------------------------------------------------------------------
C     CASE WHERE NO PLOT POSSIBLE
C-----------------------------------------------------------------------
C
      IF (IDRAW.EQ.0 .OR. amin.EQ.amax .OR. IA2.LT.IA1) THEN
        write(6,*) 'DRAW: ERROR: NO PLOT POSSIBLE'
        write(6,*) 'DRAW: ERROR: AMIN=',amin,' AMAX=',amax,
     >       ' IA1,IA2=',ia1,ia2,' IDRAW=',idraw
        GOTO 9999
C
C-----------------------------------------------------------------------
C     CASE WHERE AN ORDINARY, UNNORMALISED PLOT IS REQUIRED
C-----------------------------------------------------------------------
c
c     IDRAW=9 - plot as points instead of connected lines 
C
      ELSEIF (IDRAW.EQ.1.or.idraw.eq.7.or.idraw.eq.9) THEN
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
c
c        write(6,'(a,4(1x,g12.5))')  'BMIN,BMAX0:',bmin,bmax,bstart,bend
c

        BMIN = MAX (BMIN, BSTART)
        BMAX = MIN (BMAX, BEND)
c slmod begin
        IF (slopt2.EQ.1.OR.slopt2.EQ.2) THEN
          IF (bmin.LT.0.0) THEN
            bmin = bmin - MAX(-0.1*bmin,0.1*bmax)
          ELSE
            bmin = bmin * 0.95
          ENDIF
          bmax = bmax * 1.10


c...      Super-power override:
          IF (bstart.NE.-HI) bmin = bstart
          IF (bend  .NE. HI) bmax = bend

        ENDIF
c slmod end
c
c        write(6,'(a,4(1x,g12.5))')  'BMIN,BMAX1:',bmin,bmax,bstart,bend
c
c       Modify BMIN and BMAX for experimental data 
c
        if (iexpt_plot.ne.0) then 
           bmin = min(bmin,expt_min)
           bmax = max(bmax,expt_max)
        endif
c
c        write(6,'(a,4(1x,g12.5))')  'BMIN,BMAX2:',bmin,bmax,expt_min,
c     >                                  expt_max
c
c
c       Set minimum of scale to zero for idraw option = 7
c
        if (idraw.eq.7) bmin = 0.0
c
        CALL GRTSET (TITLE,REF,VIEW,PLANE,JOB,AMIN,AMAX,BMIN,BMAX,
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
 
        DO 120 IB = 1, NBS
          IF (IGS(IB).GT.0) then 
c
c           Plot as LINE 
c
            if (idraw.eq.1.or.idraw.eq.7) then  
              CALL GRTRAC (AS(IA1),CS(IA1,IB),IA2-IA1+1,BLABS(IB),
     >                     'LINE',1)
c
c           Plot as POINT
c
            elseif (idraw.eq.9) then   

              CALL GRTRAC (AS(IA1),CS(IA1,IB),IA2-IA1+1,BLABS(IB),
     >                     'POINT',1)

            endif
c
          endif
c
  120   CONTINUE
C
C-----------------------------------------------------------------------
C     CASE WHERE A NORMALISED PLOT IS REQUIRED
C-----------------------------------------------------------------------
C
c     IDRAW = 10 - plot as points instead of connected lines
c
      ELSEIF (IDRAW.eq.2.or.idraw.eq.3.or.idraw.eq.5.or.
     >        idraw.eq.6.or.idraw.eq.10) THEN
C
C------ CALCULATE MIN AND MAX FUNCTION VALUES WITHIN THE A RANGE.
C------ VALUES LESS THAN 0 WILL NOT BE PLOTTED  (GRAPH FROM 0 TO 1)
C------ DRAW FRAMEWORK FOR GRAPH.
C------ SCALE VALUES AND PLOT RESULTS
C
        NEGTIV = .FALSE.
        DO 210 IB = 1, NBS
          BMAX = 0.0
          BMIN = 0.0
          DO 200 IA = IA1, IA2
            BMAX = MAX (BMAX, CS(IA,IB))
            BMIN = MIN (BMIN, CS(IA,IB))
  200     CONTINUE
          FACTS(IB) = MAX (BMAX, ABS(BMIN))
          IF (BMIN.LT.0.0.AND.BSTART.LT.0.0) NEGTIV = .TRUE.
  210   CONTINUE

c
c       One off option to get Te and Ti on the same scale
c
        if (idraw.eq.6) then
c
c          IB = 3 = Te   IB = 4 = Ti
c
           facts(3) = max(facts(3),facts(4))
           facts(4) = facts(3)
c
        endif

        IF (NEGTIV) THEN
          CALL GRTSET (TITLE,REF,VIEW,PLANE,JOB,AMIN,AMAX,-1.0,1.0,
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
        ELSE
          CALL GRTSET (TITLE,REF,VIEW,PLANE,JOB,AMIN,AMAX,0.0,1.0,
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
        ENDIF
        DO 230 IB = 1, NBS
          DO 220 IA = IA1, IA2
            NORMS(IA) = 0.0
            IF (FACTS(IB).NE.0.0) NORMS(IA) = CS(IA,IB) / FACTS(IB)
  220     CONTINUE
          IF (AREAS(IB).NE.0.0.and.(idraw.ne.5.and.idraw.ne.6))
     >      WRITE (BLABS(IB)(12:20),'(1P,E9.2)') AREAS(IB)
          IF (FACTS(IB).NE.0.0)
     >      WRITE (BLABS(IB)(22:30),'(1P,E9.2)') FACTS(IB)


          
          IF (IGS(IB).GT.0) then 

             if (idraw.eq.10) then  
                CALL GRTRAC (AS(IA1),NORMS(IA1),IA2-IA1+1,BLABS(IB),
     >                       'POINT',1)
             else 
                CALL GRTRAC (AS(IA1),NORMS(IA1),IA2-IA1+1,BLABS(IB),
     >                       'LINE',1)
             endif
          endif

  230   CONTINUE
c
c     Plot points as well as lines!
c
      ELSEIF (IDRAW.EQ.4) THEN
C
C------ CALCULATE MIN AND MAX FUNCTION VALUES WITHIN THE A RANGE.
C------ DRAW FRAMEWORK FOR GRAPH.
C------ PLOT RESULTS, TRUNCATED TO REGION BSTART TO BEND ...
C
        BMIN = BEND
        BMAX = BSTART
        DO  IB = 1, NBS
          DO  IA = IA1, IA2
            BMAX = MAX (BMAX, CS(IA,IB))
            BMIN = MIN (BMIN, CS(IA,IB))
          end do
        end do
        BMIN = MAX (BMIN, BSTART)
        BMAX = MIN (BMAX, BEND)
c
c       Modify BMIN and BMAX for experimental data 
c
        if (iexpt_plot.ne.0) then 
           bmin = min(bmin,expt_min)
           bmax = max(bmax,expt_max)
        endif
c
        CALL GRTSET (TITLE,REF,VIEW,PLANE,JOB,AMIN,AMAX,BMIN,BMAX,
     >               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
        DO  IB = 1, NBS
          IF (IGS(IB).GT.0) then
           CALL GRTRAC (AS(IA1),CS(IA1,IB),IA2-IA1+1,BLABS(IB),'LINE',1)
           CALL GRTRAC(AS(IA1),CS(IA1,IB),IA2-IA1+1,BLABS(IB),'POINT',1)
          endif
        end do
c slmod begin - new
      ELSEIF (IDRAW.EQ.8) THEN
        BMIN = BEND
        BMAX = BSTART
        DO IB = 1, NBS
          DO IA = IA1, IA2
            BMAX = MAX (BMAX, CS(IA,IB))
            BMIN = MIN (BMIN, CS(IA,IB))
          ENDDO
        ENDDO
        BMIN = MAX (BMIN, BSTART)
        BMAX = MIN (BMAX, BEND)

        IF (slopt2.EQ.2) bmax = bmax * 1.10
c
        CALL GRTSET (TITLE,REF,VIEW,PLANE,JOB,AMIN,AMAX,BMIN,BMAX,
     .               TABLE,AAXLAB,BAXLAB,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
        DO IB = 1, NBS
c...      Build large polygon:          
           I1 = 0
           IF (IB.EQ.1) THEN
             I1 = I1 + 2
             RV(1) = as(IA2)
             ZV(1) = bmin
             RV(2) = as(IA1)
             ZV(2) = bmin
           ELSE
             DO IA = IA2, IA1, -1
             I1 = I1 + 1
               RV(I1) = AS(IA)
               ZV(I1) = CS(IA,IB-1)
             ENDDO
           ENDIF
           DO IA = IA1, IA2
             IF (I1.EQ.999) THEN
               WRITE(6,*) 'DRAW: RV,ZV ARRAY BOUNDS EXCEEDED'
               WRITE(0,*) 'DRAW: RV,ZV ARRAY BOUNDS EXCEEDED'
               RETURN
             ENDIF
             I1 = I1 + 1
             RV(I1) = AS(IA)
             ZV(I1) = CS(IA,IB)
           ENDDO
c           DO I2 = 1, i1
c             WRITE(6,*) 'MARK: SOLID: ',i2,rv(i1),zv(i1)
c           ENDDO
           CALL GRTRAC (RV,ZV,I1,BLABS(IB),'SOLID',1)
         ENDDO
c slmod end


      ENDIF

c
c     If the print option is turned on - then print the graph values
c     opposite their indices in the print file ... set for now
c     to fort.24
c

      if (cgrprint.eq.1) then
         call grprint(as,cs,ia1,ia2,mxxnas,nbs,igs,title,blabs,ref)
      endif
c slmod begin - new
c...  Print comments:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL CTRMAG (12)
      DO i = 1, 10
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.590+(i-1)*0.02,char(i))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST (1.04,0.250-(j-1)*0.02,char(i))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
      ENDDO
c slmod end

c
C-----------------------------------------------------------------------
c     ADD EXPERIMENTAL DATA TO THE PLOT IF IEXPT IS NON-ZERO 
C-----------------------------------------------------------------------
c
      if (iexpt_plot.ne.0) then 
c
         call plot_allexpt(iexpt_plot) 
c
      endif
C
C-----------------------------------------------------------------------
C     ALL CASES: FINISH GHOST PICTURE IF REQUIRED
C-----------------------------------------------------------------------
C
      IF (IFLAG.EQ.2.OR.IFLAG.EQ.3.OR.IFLAG.EQ.5.OR.IFLAG.EQ.6.OR.
     >    IFLAG.EQ.7) CALL FRAME
 9999 RETURN
      END
c
c
c
      subroutine grprint (as,cs,ia1,ia2,mxxnas,nbs,igs,title,blabs,ref)
      implicit  none
      INTEGER   MXXNAS,NBS,IGS(*),ia1,ia2
      REAL      AS(*),CS(MXXNAS,*)
      CHARACTER TITLE*80,ref*36
      CHARACTER BLABS(*)*36
c
c     This routine prints the values that are plotted in the graph in a
c     two column format ... the first column is the "X" axis and the
c     second is the "Y" axis. These columns of data can then be
c     extracted and plotted using a spreadsheet program. The purpose
c     of this is to allow some "easy" collation of results from
c     different cases without incurring a lot of overhead since
c     one would expect to need this type of collation primarily
c     for presentations and reports.
c
c     David Elder,   1994 April 19.
c
      integer in,ia,ib,iout
      parameter (iout=26)

      do 10 ib = 1,nbs
         write(iout,100)
         write(iout,*) title
         write(iout,*) ref,blabs(ib)
         write(iout,100)
         do 10 ia = ia1,ia2
            write (iout,110) as(ia),cs(ia,ib)
 10   continue
      write (iout,100)
 100  format (10x)
 110  format (e14.8,',',e14.8)
      return
      end
C
C
C
      SUBROUTINE GRTSET (TITLE,REF,VIEW,PLANE,JOB,XMIN,XMAX,
     >    YMIN,YMAX,TABLE,XLABEL,YLABEL,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
      use mod_params
      use mod_colours
      use mod_comgra
      use mod_slout
      implicit none
      REAL      XMIN,XMAX,YMIN,YMAX
      INTEGER   IFLAG,IDRAW,NBBS
c slmod begin - new
c...  Allow for longer y-axis labels:
      CHARACTER YLABEL*(*),XLABEL*36,TABLE*36
      CHARACTER TITLE*(*),JOB*72,SMOOTH*72
c
c      CHARACTER TITLE*80,JOB*72,SMOOTH*72
c      CHARACTER YLABEL*36,XLABEL*36,TABLE*36
c slmod end
      character REF*(*),VIEW*(*),PLANE*(*),anly*(*)
C
C  *********************************************************************
C  *                                                                   *
C  * GRTSET:  DRAWS FRAMES, TITLE, X AND Y AXES AND ASSOCIATED LABELS  *
C  *          FOR A NEW PICTURE AND INITIALISES COMMON BLOCK COMGRA FOR*
C  *          ROUTINE GRTRAC.                                          *
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
C  *                                                                   *
C  *********************************************************************
C
c      include 'params'
c      include 'slout'
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c      include 'comgra'
c slmod begin
      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost
c slmod end
c
c     Declare variables
c
      real power,tmin,tmax
      integer iexp,l,lenstr,iten
c slmod begin
      INTEGER L1

c...
c      IF (slopt4.EQ.1) THEN
c        WRITE(0,*) 'CALLING GRTSET_TRIM'
c        CALL GRTSET_TRIM (TITLE,REF,VIEW,PLANE,JOB,XMIN,XMAX,
c     .    YMIN,YMAX,TABLE,XLABEL,YLABEL,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
c        RETURN
c      ENDIF 
c
c      CALL THICK2(4)
c slmod end
C
      IF (IFLAG.EQ.4.OR.IFLAG.EQ.5) RETURN
C     ====================================
C
      CXMIN  = XMIN
      CXMAX  = XMAX
      CYMIN  = YMIN
      CYMAX  = YMAX
      IPLOTS = 0
      COL   = init_col()
      NPLOTS = NBBS
      ISPOT  = 12
      IF (NPLOTS.GT.10) ISPOT = 10
      IF (NPLOTS.GT.15) ISPOT = 8
c 
      write(6,*) 'NPLOTS:',nplots,ispot
C
C---- DRAW TITLES
C
      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRMAG (20)
c
c     added for AUG
c
c      CALL BACCOL(5)
c
      CALL LINCOL (defcol)
      CALL THICK  (2)
c slmod begin - new
      L = LENSTR(TITLE)
      IF (L.GT.76) THEN
        DO L1 = 76, 1, -1
          IF (TITLE(L1:L1).EQ.' ') EXIT
        ENDDO
c...    Maybe there are no spaces, so set L1 arbitrarily:
        IF (L1.EQ.1) L1 = 76
        CALL PCSCEN (0.68, 0.97, TITLE(1:L1-1))   
        CALL PCSCEN (0.68, 0.94, TITLE(L1+1:L))           
      ELSE
        CALL PCSCEN (0.68, 0.95, TITLE(:L))
      ENDIF
c
c      L = LENSTR(TITLE)
c      CALL PCSCEN (0.8, 0.95, TITLE(:L))
c slmod end
      L = LENSTR(TABLE)
      CALL PCSCEN (1.18, 0.87, TABLE(:L))
      L = LENSTR (XLABEL)
c slmod begin
      IF (iopt_ghost.EQ.0) THEN
        IF (IFLAG.NE.3) CALL PCSCEN (0.5,0.025,XLABEL(:L))
      ELSE
        CALL CTRMAG(12)
        IF (IFLAG.NE.3) CALL PCSCEN (0.5*(MAP1X+MAP2X),MAP1Y-0.05,
     .                               XLABEL(:L))
      ENDIF
c
c      IF (IFLAG.NE.3) CALL PCSCEN (0.5,0.025,XLABEL(:L))
c slmod end
      CALL THICK  (1)
      CALL CTRMAG (12)
      L = LENSTR (REF)
      CALL PCSCEN (1.16, 0.12, REF(:L))
      CALL CTRMAG (12)
      CALL PCSCEN (1.16, 0.17, JOB(37:72))
      CALL PCSCEN (1.16, 0.22, JOB( 1:36))
      CALL CTRMAG (12)
      L=LENSTR(VIEW)
      CALL PCSCEN (1.16, 0.27, VIEW(:L))
      L=LENSTR(PLANE)
      CALL PCSCEN (1.16, 0.31, PLANE(:L))
      L=LENSTR(ANLY)
      CALL PCSCEN (1.16, 0.35, ANLY(:L))
      CALL PCSCEN (1.16, 0.39, SMOOTH(37:72))
      CALL PCSCEN (1.16, 0.43, SMOOTH( 1:36))
      CALL CTRMAG (ISPOT)
c
      IF     (IDRAW.EQ.2.or.idraw.eq.10) THEN
        CALL PLOTST (1.05 ,0.838, '        TOT AREA   SCALE        ')
      ELSEIF (IDRAW.EQ.3) THEN
        CALL PLOTST (1.05 ,0.838, '      0:0.4 AREA   SCALE        ')
      ELSEIF (IDRAW.EQ.5) THEN
        CALL PLOTST (1.05 ,0.838, '                   SCALE        ')
      ENDIF
C
C---- DRAW FRAMES
C
      CALL LINCOL (defcol)
c
c      CALL LINCOL (3)
c
c slmod begin
c...  Allow non-standard sized plots:
      CALL POSITN (MAP1X,MAP1Y)
      CALL   JOIN (MAP1X,MAP2Y)
      CALL   JOIN (MAP2X,MAP2Y)
      CALL   JOIN (MAP2X,MAP1Y)
      CALL   JOIN (MAP1X,MAP1Y)
c
c      CALL POSITN (0.1, 0.1)
c      CALL   JOIN (0.1, 0.9)
c      CALL   JOIN (0.9, 0.9)
c      CALL   JOIN (0.9, 0.1)
c      CALL   JOIN (0.1, 0.1)
c slmod end
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
c slmod begin
c...  Allow non-standard sized plots:
c      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
c
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c slmod end
      CALL MAP    (TMIN, TMAX, 0.1, 0.9)
      CALL CTRMAG (10)
c slmod begin
      IF (iopt_ghost.NE.2) THEN
        IF (OPT_XSCALE.EQ.2) THEN
          CALL XLOGSCALE
        ELSE
          CALL XSCALE
        ENDIF
      ENDIF
c
c      CALL XSCALE
c slmod end
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
c slmod begin
c...True shape
c      CALL PSPACE (MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c
c      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod begin
      CALL MAP    (0.0, 1.0, TMIN, TMAX)
      CALL LINCOL (defcol)
      CALL CTRMAG (10)
c slmod begin
      IF (OPT_YSCALE.EQ.2) THEN
        CALL YLOGSCALE
      ELSE
        CALL YSCALE
      ENDIF

c...  True shape
      CALL PSPACE (0.0, 1.35, 0.1, 0.9)
c
c      CALL YSCALE
c      CALL PSPACE (0.0, 1.35, 0.11, 0.89)
c slmod end
      CALL MAP    (0.0, 1.35, 0.0, 1.0)
      CALL CTRORI (90.0)
c slmod begin
      IF (ITEN .NE. 0.AND.OPT_YSCALE.NE.2) THEN
c      IF (ITEN .NE. 0) THEN
c slmod end
         CALL POSITN (0.02, 0.01)
         CALL CTRMAG (14)
         CALL TYPECS ('X10')
         CALL CTRMAG (12)
         CALL TYPENI ((ITEN))
      ENDIF
      CALL LINCOL (defcol)
      CALL CTRMAG (14)
      CALL THICK  (2)
      L = LENSTR (YLABEL)
c slmod begin
      IF (iopt_ghost.EQ.0) THEN
        CALL PCSEND (0.02,0.99,YLABEL(:L))
      ELSE
        CALL CTRMAG(11)
        CALL PSPACE (0.0,1.35,map1y,map2y)
        CALL MAP    (0.0,1.35,0.0  ,1.0  )
        CALL PLOTST (0.04,0.1,YLABEL(:L))
        CALL PSPACE (0.0, 1.35, 0.1, 0.9)
        CALL MAP    (0.0, 1.35, 0.0, 1.0)         

c        CALL PCSCEN (MAP1X-0.06,0.5*(MAP1Y+MAP2Y),YLABEL(:L))   
      ENDIF
c slmod end
      CALL CTRORI (0.0)
      CALL THICK  (1)
C
C---- DRAW LINE FOR FUNCTION=0 IN CASES WHERE IFLAG=6  (EG NET EROSION)
C
      IF (IFLAG.EQ.6) THEN
        CALL LINCOL (defcol)
c slmod begin
        CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c
c        CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod end
        CALL MAP    (CXMIN,CXMAX,CYMIN,CYMAX)
        CALL POSITN (CXMIN, 0.0)
        CALL JOIN   (CXMAX, 0.0)
      ENDIF
      RETURN
      END
C
C
C
      SUBROUTINE GRTRAC (X ,Y ,NPTS ,NAME, CURVE, INC)
      use mod_params
      use mod_colours
      use mod_comgra
      use mod_slout
      implicit  none
c slmod begin - new
      CHARACTER NAME*(*),CURVE*(*)
c
c      CHARACTER NAME*36,CURVE*(*)
c slmod end
      INTEGER   NPTS,INC
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
c  *  INC   - INCREMENT PLOT COUNT AND INCLUDE IN SYMBOL TABLE -> 1=ON *
C  *                                                                   *
C  *********************************************************************
C
c      include 'params'
c      include 'slout'
c
c      include 'comgra'
c      include 'colours'
c      integer init_col,get_col,next_col,save_col
c      external init_col,get_col,next_col
c
c slmod begin
      LOGICAL firstpoint
 
      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER i1,mark
      REAL    spot,dspot,x1,y1,xwid,ywid
c
      integer nsymbols,nsym
      character*4 symbols
c
c      real spot
c slmod end
      integer i,j,ibrok
C
c      INTEGER COLOUR(8)
c
c      DATA COLOUR /2,4,6,5,7,3,6,8/
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
c      CALL RGB
c      CALL COLSET(0.7,1.0,0.,8)
c      CALL COLSET(1.0,0.7,0.,9)
C
C     WRITE (6,'('' GRTRAC: X ='',/,(1X,8F9.5))')    (X(I),I=1,NPTS)
C     WRITE (6,'(''     AND Y ='',/,1P,(1X,8E9.2))') (Y(I),I=1,NPTS)
C
C---- INCREMENT COUNTS, SET COLOUR AND LINE PATTERN
C---- ARGUMENTS TO "BROKEN" GIVE SIZES OF  (DASH1, GAP1, DASH2, GAP2)
C
      nsymbols = 3
      symbols = '+XO'
c
      nsym = mod(iplots,nsymbols) +1
c
      IBROK  = 5 * IPLOTS
c
      IF (INC.EQ.1) then 
         IPLOTS = IPLOTS + 1
      endif  
c
      if (inc.eq.-1) then
         call lincol(defcol)
      else
         CALL LINCOL (col)
         IF (INC.EQ.1) COL   = next_col()
      endif
c
c     Done automatically in next_col()
c
c      IF (ICOL.LT.1) ICOL = NCOLS
c
      IF (IPLOTS.LE.1) THEN
         CALL FULL
c         write(6,*) 'FULL:',iplots,ibrok,nplots
      ELSEIF (NPLOTS.LE.5) THEN
         CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
c         write(6,*) 'BROKEN1:',ibrok,3*IBROK,2*IBROK,3*IBROK,2*IBROK
      ELSE
         CALL BROKEN (2*IBROK,1*IBROK,2*IBROK,1*IBROK)
c         write(6,*) 'BROKEN2:',ibrok,2*IBROK,1*IBROK,2*IBROK,1*IBROK
      ENDIF
c slmod begin
      IF ((slopt2.EQ.1.OR.slopt2.EQ.2).AND.plottype(iplots).GT.1) THEN
        CALL FULL
        col = plottype(iplots)+ncols
        CALL LINCOL(col)
      ENDIF
c slmod end

C
C---- DRAW PLOT
C
c slmod begin
c...  Allows for plots to be drawn in a subregion of the 
c     page:
c      IF (slopt.GT.0) THEN
c        CALL PSPACE (MAP1X,MAP2X,MAP1Y,MAP2Y)
c      ELSE
        CALL PSPACE (0.102, 0.898, 0.102, 0.898)
c      ENDIF
c
c      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod end
      CALL MAP    (CXMIN,CXMAX,CYMIN,CYMAX)
C
c slmod begin - new
      IF (iplots.GT.0.AND.plottype(iplots).LE.-1) THEN

        xwid = map2x - map1x
        ywid = map2y - map1y

        CALL HATOPT(1)
        CALL FILCOL(ABS(plottype(iplots))+ncols)
        CALL LINCOL(ABS(plottype(iplots))+ncols)
c        CALL LINCOL(1)
        CALL FULL
        firstpoint = .TRUE.
        DO I = 1, NPTS
          IF (X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.Y(I).GE.CYMIN.AND.
     .        Y(I).LE.CYMAX.AND.Y(I).NE.LO)
     .        CALL BOX (x(i)-0.005*(CXMAX-CXMIN),
     .                  x(i)+0.005*(CXMAX-CXMIN),
     .                  y(i)-0.005*(CYMAX-CYMIN)*xwid/ywid,
     .                  y(i)+0.005*(CYMAX-CYMIN)*xwid/ywid)

c          IF (y(i).NE.LO) THEN
c            IF (firstpoint) THEN
c              CALL POSITN(x(i),y(i))
c              firstpoint = .FALSE.
c            ELSE
c              CALL JOIN(x(i),y(i))
c            ENDIF
c          ENDIF

c     .     CALL BOX (x(i)+0.015-0.040,x(i)+0.015+0.040,
c     .               y(i)+0.002-0.140,y(i)+0.002+0.140)

        ENDDO
      ELSEIF (CURVE.EQ.'POINT') THEN
        DO 10 I = 1, NPTS
c          IF (X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.Y(I).GE.CYMIN.AND.
c     >        Y(I).LE.CYMAX.AND.Y(I).NE.LO) CALL PLOTST (X(I),Y(I),'+')
          IF (X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.Y(I).GE.CYMIN.AND.
     >        Y(I).LE.CYMAX.AND.Y(I).NE.LO) CALL PLOTST (X(I),Y(I),
     >        symbols(nsym:nsym))
   10   CONTINUE
      ELSEIF (CURVE.EQ.'SOLID') THEN
c...    Plot a solid region for the data.  This is presently only
c       useable with plot 980, sub-plots 4 and 10:
        CALL FULL
        CALL FILCOL(col)
        CALL LINCOL(col)
        CALL PTPLOT(x,y,1,npts,1)         
      ELSE
        IF (OPT_XSCALE.EQ.2.OR.OPT_YSCALE.EQ.2) THEN
          MARK = 0
          CALL CTRMAG (ISPOT)
          DO I1 = 1, NPTS
            x1 = x(i1)
            y1 = y(i1)
            IF (opt_xscale.EQ.2) x1 = LOG10(x1)
            IF (opt_yscale.EQ.2) y1 = LOG10(y1)
            IF (Y1.NE.LO.AND.Y1.NE.-37) THEN
              CALL POSITN (X1, Y1)
              GOTO 15
            ENDIF
          ENDDO
 15       DO I = I1+1, NPTS
            x1 = x(i)
            y1 = y(i)
            IF (opt_xscale.EQ.2) x1 = LOG10(x1)
            IF (opt_yscale.EQ.2) y1 = LOG10(y1)
            IF (y1.NE.LO.AND.y1.NE.-37.0) THEN
              CALL JOIN (X1,Y1)
              IF (mark.EQ.0) THEN
                CALL PLOTST (X1,Y1,NAME(1:4))
                mark = 1
              ENDIF
            ENDIF
          ENDDO
        ELSE
c...      Original line plotting code:
          I = 1
          DO WHILE(Y(I).EQ.LO.AND.I.LT.NPTS) 
            I = I + 1
          ENDDO
          CALL POSITN (X(I), Y(I))
          DO 20 J = I+1, NPTS
            IF (Y(J).NE.LO) THEN
              CALL JOIN (X(J),Y(J))
              IF (J.EQ.I+3.OR.J.EQ.NPTS-5) CALL PLOTST(X(J),Y(J),
     .                                                 NAME(1:4))
            ENDIF
 20       CONTINUE
        ENDIF
      ENDIF
c
c      IF (CURVE.EQ.'POINT') THEN
c        DO 10 I = 1, NPTS
c          IF (X(I).GE.CXMIN.AND.X(I).LE.CXMAX.AND.Y(I).GE.CYMIN.AND.
c     >        Y(I).LE.CYMAX) CALL PLOTST (X(I), Y(I), '+')
c   10   CONTINUE
c      ELSE
c        CALL POSITN (X(1), Y(1))
c        DO 20 I = 2, NPTS
c          CALL JOIN (X(I), Y(I))
c          IF (I.EQ.3.OR.I.EQ.NPTS-5) CALL PLOTST (X(I),Y(I),NAME(1:4))
c 20     CONTINUE
c      ENDIF
c slmod end
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
c slmod begin
      IF (iopt_ghost.EQ.2) RETURN

      IF (INC.NE.1) RETURN

      IF (slopt2.EQ.2) THEN
        ispot = 20
        CALL THICK(2)
      ENDIF

      IF (plottype(iplots).LE.-1) THEN

        CALL THICK(1)

        DSPOT = REAL(ISPOT) / 500.0
        SPOT = 0.838 - IPLOTS * DSPOT

        CALL CTRMAG (ISPOT)
        CALL PSPACE (0.0, 1.35, 0.0,1.0)
        CALL CSPACE (0.0, 1.35, 0.0,1.0)
        CALL MAP    (0.0, 1.35, 0.0,1.0)
        SPOT = 0.838 - REAL (ISPOT*IPLOTS) / 500.0

        CALL HATOPT(1)
        CALL FILCOL(ABS(plottype(iplots))+ncols)
        CALL LINCOL(ABS(plottype(iplots))+ncols)
        CALL FULL
        CALL BOX (0.99-0.006,
     .            0.99+0.006,
     .            spot-0.006*1.1+0.002,
     .            spot+0.006*1.1+0.002)

      ELSEIF (CURVE.EQ.'SOLID') THEN
        CALL CTRMAG (ISPOT)
        CALL PSPACE (0.0, 1.35, 0.0,1.0)
        CALL CSPACE (0.0, 1.35, 0.0,1.0)
        CALL MAP    (0.0, 1.35, 0.0,1.0)
        DSPOT = REAL(ISPOT) / 500.0
        SPOT = 0.838 - IPLOTS * DSPOT
        CALL BOX (1.00-DSPOT/2.0,1.00+DSPOT/2.0,
     .            SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
        CALL FULL
      ELSE
        CALL CTRMAG (ISPOT)
        CALL PSPACE (0.0, 1.35, 0.0,1.0)
        CALL CSPACE (0.0, 1.35, 0.0,1.0)
        CALL MAP    (0.0, 1.35, 0.0,1.0)
        SPOT = 0.838 - REAL (ISPOT*IPLOTS) / 500.0
        IF (CURVE.EQ.'POINT') THEN
c          CALL PLOTST (0.96,SPOT,'+')
c          CALL PLOTST (0.99,SPOT,'+')
c          CALL PLOTST (1.02,SPOT,'+')
          CALL PLOTST (0.96,SPOT,symbols(nsym:nsym))
          CALL PLOTST (0.99,SPOT,symbols(nsym:nsym))
          CALL PLOTST (1.02,SPOT,symbols(nsym:nsym))
        ELSE
          CALL POSITN (0.95,SPOT)
          CALL JOIN   (1.03,SPOT)
        ENDIF
      ENDIF

      CALL THICK(1)
c
c      IF (INC.NE.1) RETURN
c      CALL CTRMAG (ISPOT)
c      CALL PSPACE (0.0, 1.35, 0.0,1.0)
c      CALL CSPACE (0.0, 1.35, 0.0,1.0)
c      CALL MAP    (0.0, 1.35, 0.0,1.0)
c      SPOT = 0.838 - REAL (ISPOT*IPLOTS) / 500.0
c      IF (CURVE.EQ.'POINT') THEN
c        CALL PLOTST (0.96,SPOT,'+')
c        CALL PLOTST (0.99,SPOT,'+')
c        CALL PLOTST (1.02,SPOT,'+')
c      ELSE
c        CALL POSITN (0.95,SPOT)
c        CALL JOIN   (1.03,SPOT)
c      ENDIF
c slmod end
      CALL FULL
c
      call lincol(defcol)
c
c slmod begin - new
      CALL PLOTST (1.05 ,SPOT, NAME(5:LEN_TRIM(NAME)))
c
c       CALL PLOTST (1.05 ,SPOT, NAME(5:36))
c slmod end
      RETURN
      END
C
C
C
      SUBROUTINE GR3D (SURFAS,NPTS,NAME,IVEW3D,PROJ3D,IBAS3D,
     >                 SUREDG,LIMEDG)
      use mod_colours
      use mod_comgra
      implicit none
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
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c      include 'comgra'
c
      real spot
c
C
C---- SURCOL: TOP COLOUR, UNDERSIDE COLOUR, BASE COLOUR
C---- (1:8 REPRESENT BLACK,RED,GREEN,BLUE,WHITE,CYAN,MAGENTA,YELLOW)
C---- SURDIR: 0 TO 3 GIVES 0 DEGREES, 90,180,270 DEGREES VIEW
C---- SURANG: 35.2644 ISOMETRIC, OTHERWISE IN RANGE -90 TO 90 DEGREES
C---- SURBAS: 0/1 UNDERSIDE OFF/ON, 0/1/2 BASE OFF/ON/ON AT LAST ARG
C---- SURPLT: DIM ARRAY(192,192), BUT ONLY PLOT ARRAY(1:NPTS,1:NPTS)
C
      CALL BACCOL(5)
c
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
      CALL LINCOL (defcol)
      CALL CTRMAG (ISPOT)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0
C     SPOT = 0.818 - 0.03 * IPLOTS
      CALL POSITN (0.95,SPOT)
      CALL JOIN   (1.03,SPOT)
      CALL FULL
      CALL PLOTST (1.05 ,SPOT, NAME)
      RETURN
      END
C
C
C
      SUBROUTINE GRCONT (VALS,IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,
     >                   MAXNYS,CLEVEL,XOUTS,YOUTS,NAME)
      use mod_colours
      use mod_comgra
      implicit none
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
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c      include 'comgra'
c
      real clev(1)
c
      real spot
c
      integer         iye,iyb,ixe,ixb,ibrok
c slmod begin
      COMMON /DUMCOM/ map1x_d,map2x_d,map1y_d,map2y_d,tag_d
      INTEGER         tag_d
      REAL            map1x_d,map2x_d,map1y_d,map2y_d

      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

c
c      slmod end
C
c
c      INTEGER COLOUR(8),IXB,IXE,IYB,IYE
c
c      DATA COLOUR /2,4,6,5,7,3,6,8/
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
      clev(1) = clevel

      WRITE (6,'(1X,A32)') NAME(5:36)
C     WRITE (6,'('' GRCONT: IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS'',
C    >  /7X,6I7)') IXMIN,IXMAX,MAXNXS,IYMIN,IYMAX,MAXNYS
C
C---- INCREMENT COUNTS, SET COLOUR AND LINE PATTERN
C
      IBROK  = IPLOTS
c
      IPLOTS = IPLOTS + 1
c
c
      COL = get_col()
c
      CALL LINCOL (COL)
c
      col = next_col()
c
c      ICOL   = ICOL - 1
c      IF (ICOL.LT.1) ICOL = NCOLS
c
      IF (IPLOTS.LE.1) THEN
         CALL FULL
      ELSEIF (NPLOTS.LE.5) THEN
         CALL BROKEN (3*IBROK,2*IBROK,3*IBROK,2*IBROK)
      ELSE
         CALL BROKEN (2*IBROK,1*IBROK,2*IBROK,1*IBROK)
      ENDIF
C
C---- DRAW CONTOUR PLOT  (DISABLE ANNOTATION)
C---- SPLIT UP X,Y RANGES IF REQUIRED ...
C
c slmod begin
c...True space
      IF (tag_d.NE.0) THEN
        CALL PSPACE (MAP1X_D,MAP2X_D,MAP1Y_D,MAP2Y_D)
      ELSE
        CALL PSPACE (0.1, 0.9, 0.1, 0.9)
      ENDIF

      iconcol = iconcol + 1
c
c      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod end
      CALL MAP    (CXMIN, CXMAX, CYMIN, CYMAX)
      CALL CONOTA (0)
      DO 110 IXB = IXMIN, IXMAX, 1024
        IXE = MIN (IXMAX, IXB+1024)
        DO 100 IYB = IYMIN, IYMAX, 1024
          IYE = MIN (IYMAX, IYB+1024)
c...dev (dummy needed)
c          WRITE(0,*) 'IXB,IXY=',ixmin,ixmax
c          CALL CONTIL2(VALS,IXB,IXE,MAXNXS,IYB,IYE,MAXNYS,
c     >                 CLEVEL,1,1,XOUTS,YOUTS)
          CALL CONTIL2(VALS,IXB,IXE,MAXNXS,IYB,IYE,MAXNYS,
     >                 CLEV,1,1,XOUTS,YOUTS)
c          CALL CONTIL (VALS,IXB,IXE,MAXNXS,IYB,IYE,MAXNYS,
c     >                 CLEVEL,1,1,XOUTS,YOUTS)
  100   CONTINUE
  110 CONTINUE


c       IF (icon.GT.0) THEN
c          CALL FILCOL(2)
c          CALL LINCOL(2)
c          CALL PTPLOT(xcon,ycon,1,icon,1)
c       ENDIF


c      DO 110 IXB = IXMIN, IXMAX, 192
c        IXE = MIN (IXMAX, IXB+192)
c        DO 100 IYB = IYMIN, IYMAX, 192
c          IYE = MIN (IYMAX, IYB+192)
c          CALL CONTIL (VALS,IXB,IXE,MAXNXS,IYB,IYE,MAXNYS,
c     >                 CLEVEL,1,1,XOUTS,YOUTS)
c  100   CONTINUE
c  110 CONTINUE


C
C---- WRITE ENTRY IN SYMBOL TABLE
C
c slmod begin
      IF (iopt_ghost.NE.0) THEN


      ELSE
        CALL CTRMAG (ISPOT)
        CALL PSPACE (0.0, 1.35, 0.0,1.0)
        CALL CSPACE (0.0, 1.35, 0.0,1.0)
        CALL MAP    (0.0, 1.35, 0.0,1.0)
        SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0
        CALL POSITN (0.95,SPOT)
        CALL JOIN   (1.03,SPOT)
        CALL FULL
c
        call lincol(defcol)
c
        CALL PLOTST (1.05 ,SPOT, NAME(5:36))
      ENDIF
c
c      CALL CTRMAG (ISPOT)
c      CALL PSPACE (0.0, 1.35, 0.0,1.0)
c      CALL CSPACE (0.0, 1.35, 0.0,1.0)
c      CALL MAP    (0.0, 1.35, 0.0,1.0)
c      SPOT = 0.818 - REAL (ISPOT*IPLOTS) / 500.0
c      CALL POSITN (0.95,SPOT)
c      CALL JOIN   (1.03,SPOT)
c      CALL FULL
c
c      call lincol(defcol)
c
c      CALL PLOTST (1.05 ,SPOT, NAME(5:36))
c slmod end
      RETURN
      END
C
C
C
      SUBROUTINE GRCOLR (VS,VLO,VHI,NAME)
      use mod_params
      use mod_colours
      use mod_comgra
      use mod_comxyt
      use mod_limpoly
      use mod_slout
      IMPLICIT NONE
c
c      include 'params'
c
      REAL VS(MAXNXS,MAXNYS),VLO,VHI
      CHARACTER*36 NAME
c
c      include 'comxyt'
c      include 'limpoly' 
C
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c
c      include 'comgra'
c
c      include 'slout'
C
C  *********************************************************************
C  *                                                                   *
C  *  GRCOLR  THIS ROUTINE DRAWS A FALSE COLOUR PLOT BY CHECKING       *
C  *          TO SEE IF THE VALUE OF A GIVEN PLASMA CELL IS IN THE     *
C  *          PRESENT RANGE.                                           *
C  *                                                                   *
C  *  L.D.HORTON    JULY 1993                                          *
C  *                                                                   *
C  *********************************************************************
C
C
      integer IX,IY,K,IRI,IKI,IRO,IKO
c
      REAL SPOT,DSPOT
      real xoffset
      integer yoffset
C
      WRITE (6,'(1X,A32)') NAME(5:36)
C
C---- INCREMENT COUNTS, SET COLOUR
C
      IPLOTS = IPLOTS + 1
      CALL HATOPT (1)
c
      col = get_col()
c
      CALL LINCOL (COL)
      CALL FILCOL (COL)
c
c      CALL LINCOL (COLOUR(ICOL+1))
c      CALL FILCOL (COLOUR(ICOL+1))
c
      col = next_col()
c
c      ICOL   = ICOL - 1
c      IF (ICOL.LT.1) ICOL = NCOLS
C
C---- DRAW PLOT ONE CELL AT A TIME
C
c slmod begin
c...True space
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c
c      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod eng
      CALL MAP    (CXMIN, CXMAX, CYMIN, CYMAX)
      DO 210 IX = 1, NXS
        DO 200 IY = 1, NYS
          K = KORPG(IX,IY)
          IF (K.EQ.0) GOTO 200
C
C---- IF CELL VALUE IS IN RANGE, DRAW FILLED BOX
C
          IF (VS(IX,IY).GT.VLO .AND. VS(IX,IY).LE.VHI) THEN
            CALL PTPLOT(XVERTP(1,K),YVERTP(1,K),1,NVERTP(K),1)
C           CALL PTJOIN(XVERTP(1,K),YVERTP(1,K),1,NVERTP(K),-1)
c slmod begin - temp (avoid drawing cell boundaries because):
c...        For plot 972:
            IF (slopt.EQ.3) CYCLE
c slmod end
C
C---- DRAW EDGES OF CONTOURS FOR PATTERN MAPS TO B&W
C
C---- FORWARD NEIGHBOUR IS BETWEEN VERTICES 3 AND 4
C
            IF (IY.EQ.NYS) GOTO 100
            IF (KORPG(IX,IY+1).EQ.0) GOTO 100
            IF (VS(IX,IY+1).LE.VLO) THEN
              CALL POSITN (XVERTP(3,K),YVERTP(3,K))
              CALL JOIN   (XVERTP(4,K),YVERTP(4,K))
            ENDIF
C
C---- INNER NEIGHBOUR IS BETWEEN VERTICES 4 AND 1
C
  100       IF (IX.EQ.1) GOTO 110
            IF (KORPG(IX-1,IY).EQ.0) GOTO 110
            IF (VS(IX-1,IY).LE.VLO) THEN
              CALL POSITN (XVERTP(4,K),YVERTP(4,K))
              CALL JOIN   (XVERTP(1,K),YVERTP(1,K))
            ENDIF
C
C---- BACKWARD NEIGHBOUR IS BETWEEN VERTICES 1 AND 2
C
  110       IF (IY.EQ.1) GOTO 120
            IF (KORPG(IX,IY-1).EQ.0) GOTO 120
            IF (VS(IX,IY-1).LE.VLO) THEN
              CALL POSITN (XVERTP(1,K),YVERTP(1,K))
              CALL JOIN   (XVERTP(2,K),YVERTP(2,K))
            ENDIF
C
C---- OUTER NEIGHBOUR IS BETWEEN VERTICES 2 AND 3
C
  120       IF (IX.eq.NXS) GOTO 200
            IF (KORPG(IX+1,IY).EQ.0) GOTO 200
            IF (VS(IX+1,IY).LE.VLO) THEN
              CALL POSITN (XVERTP(2,K),YVERTP(2,K))
              CALL JOIN   (XVERTP(3,K),YVERTP(3,K))
            ENDIF
          ENDIF
C
  200   CONTINUE
  210 CONTINUE
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
c
c     If there are more than 20 contours - allow for two columns 
c     in the symbol table. 
c
      if (iplots.gt.20) then
        xoffset=0.24
        yoffset=20
      else
        xoffset=0.0
        yoffset=0
      endif
c
      CALL CTRMAG (ISPOT)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      DSPOT = REAL(ISPOT) / 500.0

      SPOT = 0.818 - (IPLOTS-yoffset) * DSPOT
      CALL BOX (0.98+xoffset-DSPOT/2.0,0.98+xoffset+DSPOT/2.0,
     >          SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
      CALL FULL
      call lincol(defcol)
      CALL PLOTST (1.01+xoffset ,SPOT, NAME(5:36))
c
c      SPOT = 0.818 - IPLOTS * DSPOT
c      CALL BOX (1.00-DSPOT/2.0,1.00+DSPOT/2.0,
c     >          SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
c      CALL FULL
c
c
c      CALL PLOTST (1.05 ,SPOT, NAME(5:36))
c
      RETURN
      END
C
C
C
c slmod begin
      SUBROUTINE GRCOLRXY (VS,maxix,maxiy,nix,niy,raxis,zaxis,
     >                     VLO,VHI,vmin,vmax,NAME)
      use mod_params
      use mod_colours
      use mod_comgra
      use mod_slout
c      SUBROUTINE GRCOLRXY (VS,maxix,maxiy,nix,niy,raxis,zaxis,
c     >                     VLO,VHI,NAME)
c slmod end
      IMPLICIT NONE
c
c      include 'params'
c
      integer maxix,maxiy,nix,niy,drawcelledges
      real raxis(maxix),zaxis(maxiy)
      REAL VS(MAXix,maxiy),VLO,VHI
      CHARACTER*36 NAME
c
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c      include 'comgra'
c
c slmod begin
      REAL vmin,vmax,val,frac,frac5,mod5
c slmod end
c      include 'slout'
C
C  *********************************************************************
C  *                                                                   *
C  *  GRCOLRXY:  THIS ROUTINE DRAWS A FALSE COLOUR PLOT BY CHECKING    *
C  *          TO SEE IF THE VALUE OF A GIVEN CELL IS IN THE            *
C  *          PRESENT RANGE.                                           *
C  *                                                                   *
C  *********************************************************************
C
C
      integer Ix,Iy
c
      REAL SPOT,DSPOT
      real xoffset
      integer yoffset
      real rvert(4),zvert(4),r1,r2,z1,z2
      integer nvert
c slmod begin
c...  From plot 982:

      REAL bright,pastel

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale

c slmod end
C
      WRITE (6,'(1X,A32)') NAME(5:36)
C
C---- INCREMENT COUNTS, SET COLOUR
C
      IPLOTS = IPLOTS + 1
      CALL HATOPT (1)
c
c     Get next colour
c
      col = get_col()
c
c     Set line colour 
c
      CALL LINCOL (COL)
      CALL FILCOL (COL)
c
      col = next_col()
c slmod begin

      IF (colscale) THEN
        CALL HSV      

        val = 0.5 * (vlo + vhi)
        frac = MIN(0.98,(val - vmin) / (vmax - vmin))

        CALL ColSet(frac,1.0,1.0-(1.0-frac)**20,255)

        IF (hardscale) THEN
          IF (frac.LT.0.07) THEN
            bright = 1.0
            frac   = 1.0
            pastel = 0.0
          ELSE
            IF (frac.LE.0.27) THEN
              bright = 1.0-((0.27-frac)/(0.27-0.07))**2
              bright = MAX(0.0,bright)
            ELSE
              bright = 1.0
            ENDIF
            frac = (1.0 - frac) * 0.90
            frac = frac + 0.34
            IF (frac.GT.1.0) frac = frac - 1.0
            pastel = 1.0
          ENDIF
          CALL ColSet(frac,pastel,bright,255)
        ENDIF

c        IF (hardscale) THEN
c          IF (frac.LT.vmin) THEN
c            bright = 1.0
c            frac   = 1.0
c            pastel = 0.0
c          ELSE
c            IF (frac.LE.0.25) THEN
c              bright = 1.0-((0.25-frac)/(0.25-vmin))**2
c            ELSE
c              bright = 1.0
c            ENDIF
c            frac = (1.0 - frac) * 0.90
c            frac = frac + 0.34
c            IF (frac.GT.1.0) frac = frac - 1.0
c            pastel = 1.0
c          ENDIF
c          CALL ColSet(frac,pastel,bright,255)
c        ENDIF

        CALL LINCOL (255)
        CALL FILCOL (255)
      ENDIF

c slmod end
c
C
C---- DRAW PLOT ONE CELL AT A TIME
C
c slmod begin
c...True space
      CALL PSPACE (0.1, 0.9, 0.1, 0.9)
c
c      CALL PSPACE (0.1, 0.9, 0.11, 0.89)
c slmod end
      CALL MAP    (CXMIN, CXMAX, CYMIN, CYMAX)
      DO 210 Ix = 1, nix
        DO 200 Iy = 1, niy
C
C---- IF CELL VALUE IS IN RANGE, DRAW FILLED BOX
C
          IF (VS(Ix,Iy).GT.VLO .AND. VS(Ix,Iy).LE.VHI) THEN
c
c            Calculate and record corner points of cell.  
c	   
c            Set up cell for INTERS routine 
c	   
             nvert = 4
c	   
c            Calculate R,Z vertices around the centre of the ix,iy cell.
c	   
c            R
c	   
             if (ix.eq.1) then 
c               
                r2 = (raxis(ix+1)+raxis(ix))/2.0
                r1 = raxis(ix) + (raxis(ix)-r2)
c	   
             elseif (ix.eq.nix) then 
c	   
                r1 = (raxis(ix)+raxis(ix-1))/2.0
                r2 = raxis(ix) + (raxis(ix)-r1)
c	   
             else
c	   
                r1 = (raxis(ix)+raxis(ix-1))/2.0
                r2 = (raxis(ix+1)+raxis(ix))/2.0
             endif 
c	   
c            Z
c	   
             if (iy.eq.1) then 
c               
                z2 = (zaxis(iy+1)+zaxis(iy))/2.0
                z1 = zaxis(iy) + (zaxis(iy)-z2) 
c	   
             elseif (iy.eq.niy) then 
c	   
                z1 = (zaxis(iy)+zaxis(iy-1))/2.0
                z2 = zaxis(iy) + (zaxis(iy)-z1)
c	   
             else
c	   
                z1 = (zaxis(iy)+zaxis(iy-1))/2.0
                z2 = (zaxis(iy+1)+zaxis(iy))/2.0
c	   
             endif 
	   
c	   
c            Combine these values in the correct order to define the
c            vertices of the polygon.
c	   
             rvert(1) = r1
             zvert(1) = z1 
c	   
             rvert(2) = r1 
             zvert(2) = z2
c	   
             rvert(3) = r2 
             zvert(3) = z2
c	   
             rvert(4) = r2 
             zvert(4) = z1
c
             CALL PTPLOT(RVERT,ZVERT,1,NVERT,1)
c
C            CALL PTJOIN(RVERTP(1,K),ZVERTP(1,K),1,NVERTP(K),-1)
c slmod begin - temp (avoid drawing cell boundaries because):
c...        For plot 972,982:
            IF (slopt.EQ.3.or.slopt.eq.42.OR.slopt.EQ.4) CYCLE
c slmod end
C
C---- DRAW EDGES OF CONTOURS FOR PATTERN MAPS TO B&W
C
c           Next Higher indexed Neighbour IN R    
C
            if (ix.eq.nix) goto 100  
            IF (VS(Ix+1,Iy).LE.VLO) THEN
              CALL POSITN (R2,Z1)
              CALL JOIN   (R2,Z2)
            ENDIF
C
c           Next Lower indexed Neighbour IN R    
C
  100       IF (ix.eq.1) goto 110
            IF (VS(ix-1,iy).LE.VLO) THEN
              CALL POSITN (R1,Z1)
              CALL JOIN   (R1,Z2)
            ENDIF
C
c           Next Higher indexed in Z
C
  110       IF (iy.eq.niy) GOTO 120
            IF (VS(ix,Iy+1).LE.VLO) THEN
              CALL POSITN (R1,Z2)
              CALL JOIN   (R2,Z2)
            ENDIF
C
c           Next lower indexed in Z
C
  120       IF (iy.eq.1) GOTO 200
            IF (VS(ix,iy-1).LE.VLO) THEN
              CALL POSITN (R1,Z1)
              CALL JOIN   (R2,Z1)
            ENDIF
          ENDIF
C
  200   CONTINUE
  210 CONTINUE
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
C
C---- WRITE ENTRY IN SYMBOL TABLE
C
c
c     If there are more than 20 contours - allow for two columns 
c     in the symbol table. 
c
      if (iplots.gt.20) then
        xoffset=0.24
        yoffset=20
      else
        xoffset=0.0
        yoffset=0
      endif
c
c slmod begin
      IF (slopt.EQ.4) RETURN
c slmod end
      CALL CTRMAG (ISPOT)
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      DSPOT = REAL(ISPOT) / 500.0

      SPOT = 0.818 - (IPLOTS-yoffset) * DSPOT
      CALL BOX (0.98+xoffset-DSPOT/2.0,0.98+xoffset+DSPOT/2.0,
     >          SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
      CALL FULL
      call lincol(defcol)
      CALL PLOTST (1.01+xoffset ,SPOT, NAME(5:36))
c
c      SPOT = 0.818 - IPLOTS * DSPOT
c      CALL BOX (1.00-DSPOT/2.0,1.00+DSPOT/2.0,
c     >          SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
c      CALL FULL
c
c
c      CALL PLOTST (1.05 ,SPOT, NAME(5:36))
c
      RETURN
      END
C
C
C
c      INTEGER FUNCTION IEXP (SMIN, SMAX)
c      implicit none
c      REAL SMIN,SMAX,power
C
C  *********************************************************************
C  *                                                                   *
C  *  IEXP: RETURNS SUITABLE EXPONENT FOR A SCALE SMIN TO SMAX.        *
C  *                                                                   *
C  *********************************************************************
C
c      POWER = MAX(ABS(SMIN), ABS(SMAX))
c      IF (POWER.EQ.0.0) then
c        iexp = 0
c      elseIF (POWER .GE. 100.0) THEN
c        IEXP = INT(LOG10(POWER) - 0.3)
c      ELSE IF (POWER .LT. 0.5) THEN
c        IEXP = INT(LOG10(POWER) - 0.3)
c      ELSE
c        IEXP = 0
c      ENDIF
c      RETURN
c      END
c
c
c
      subroutine DRAWM (MOUTS,MWIDS,MVALS,MAXNMS,maxplts,maxngs,pnks,
     >                  nplts,pngs,pltlabs,mlabs,xlab,ylab,ref,title,
     >                  sctype,numplots,pltmins,pltmaxs,pltfact,naxes,
     >                  mdrawtype,drawtypesw)
      use mod_grminfo
      implicit none
      integer nplts,pngs(maxplts),sctype,drawtypesw
      integer maxnms,maxplts,maxngs,numplots,naxes
      integer mdrawtype(maxplts,maxngs)
      integer pnks(maxplts,maxngs)
      real mouts(maxnms,maxplts,maxngs)
      real mvals(maxnms,maxplts,maxngs)
      real mwids(maxnms,maxplts,maxngs)
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact
      character*(*) pltlabs(maxplts)
      character*(*) mlabs(maxplts,maxngs)
      character*(*) xlab,ylab,ref,title
c
c     Use elabs for local labeling
c
      character*36 elabs(maxngs)
c
c     Use drawtype for local drawing of plots
c
      integer drawtype(maxngs) 
c
c      include 'grminfo'
c
      real hi
      parameter(hi=1.0e37)
C
C  *********************************************************************
C  *                                                                   *
C  * DRAWM:   DRAW A GRAPH OF VALUES IN MVALS AGAINST VALUES IN MOUTS  *
C  *          OVER THE RANGE OF VALUES IN MOUTS. PUTTING UP TO 8 PLOTS *
C  *          ON ONE PAGE                                              *
C  *                                                                   *
C  * ARGUMENTS :-                                                      *
C  *                                                                   *
C  * MOUTS  - VALUES TO GO ON X AXIS OF GRAPHS (FOR EACH DATA SET)     *
C  * MWIDS  - BIN WIDTHS  (FOR SMOOTHING PURPOSES)                     *
C  * MVALS  - VALUES TO GO ON Y AXIS OF GRAPHS (FOR EACH DATA SET)     *
C  * MAXNMS - FIRST DIMENSION OF M ARRAYS (IK INDEX)                   *
C  * MAXPLTS- SECOND DIMENSION OF M ARRAYS - MAX NUM OF PLOTS          *
C  * MAXNGS - MAX NUMBER OF DATA SETS FOR EACH PLOT                    *
C  * PNKS   - NUMBER OF DATA VALUES IN EACH SET OF DATA                *
C  * NPLTS  - NUMBER OF PLOTS TO BE MADE                               *
C  * PNGS   - NUMBER OF DATA SETS ON EACH PLOT                         *
C  * PLTLABS- LABELS SPECIFIC TO EACH PLOT                             *
C  * MLABS  - LABELS SPECIFIC TO EACH DATA SET ON EACH PLOT            *
C  * XLAB   - X-AXIS LABEL                                             *
C  * YLAB   - Y-AXIS LABEL                                             *
C  * REF    - REFERENCE LABEL FOR ALL PLOTS                            *
C  * TITLE  - TITLE OF CASE                                            *
c  * PLTMINS- EXTERNAL PLOT MINIMUMS                                   *
c  * PLTMAXS- EXTERNAL PLOT MAXIMUMS                                   *
C  * SCTYPE - TYPE OF PLOT SCALING TO BE EMPLOYED              *
C  *          1 - ALL PLOTS ON PAGE SHARE THE SAME SCALING - MIN AND   *
C  *              MAX CALCULATED FOR ALL DATA                          *
C  *          2 - MINIMUM AND MAXIMUM CALCULATED FOR EACH PLOT         *
C  *              SEPARATELY ON THE PAGE.                              *
C  *          3 - MAX CALCULATED FOR ALL PLOTS ON PAGE - MIN IS ZERO   *
C  *              UNLESS PLOT CROSSES X-AXIS.                  *
C  *          4 - MAXIMUM CALCULATED FOR EACH PLOT SEPARATELY - MIN IS *
C  *              SET TO 0.0 UNLESS PLOT CROSSES X-AXIS            *
C  *          5 - Logarithmic Scaling                                  *
C  *              MAX CALCULATED FOR ALL PLOTS ON PAGE - MIN IS ZERO   *
C  *              UNLESS PLOT CROSSES X-AXIS.                  *
C  *          6 - Logarithmic Scaling                                  *
C  *              MAXIMUM CALCULATED FOR EACH PLOT SEPARATELY - MIN IS *
C  *              SET TO 0.0 UNLESS PLOT CROSSES X-AXIS            *
c  *          7 - Externally provided plot minimums and maximums.      *
c  * NUMPLOTS-NUMBER OF PLOTS ON EACH PAGE                             *
c  * NAXES  - IF SET TO 1 THIS INDICATES THAT ONLY ONE SET OF AXIS     *
C  *          COORDINATES HAVE BEEN SET FOR EACH PLOT. IN THIS CASE    *
C  *          THE FIRST SET ARE COPIED INTO THE AXIS SPACE FOR EACH OF *
C  *          THE OTHER DATA SETS. THIS OPTION WAS ADDED TO MINIMIZE   *
C  *          THE CHANGES REQUIRED TO EXISTING CODE WHEN MULTIPLE      *
C  *          AXIS CAPABILITY WAS ADDED TO THE PLOTTING CODE.          *
C  *                                           *
C  *                                           *
C  *          David Elder, NOV 21, 1995                    *
C  *                                                                   *
C  *********************************************************************
C
c
c     Local Variables
c
      integer in,ip,ik,ir,ipmax,id,posin
      real pltmax,pltmin,axmax,axmin,axfact
c
c     Set up axis factor based on value of pltfact
c
      if (pltfact.gt.1.0) then
         axfact = pltfact-1.0
         pltfact= 0.0
      else
         axfact = 0.0
      endif 
C
C     Copy the specified axis information for all plots if information
C     for only one horizontal axis has been specified.
C
      if (naxes.eq.1) then
         do ip = 1,nplts
            do in = 2,pngs(ip)
               pnks(ip,in) = pnks(ip,1)
               do ik = 1,pnks(ip,1)
                  mouts(ik,ip,in) = mouts(ik,ip,1)
               end do
            end do
         end do
      endif
c
c     If the drawtype array is not switched on then load a
c     regular unmarked line as the default
c
      if (drawtypesw.eq.0) then
         do in = 1,maxngs
            drawtype(in) = 1
         end do
      endif
c
c     Set up number of plots / page
c
      pageplots = numplots
c
c     Adjust size of text depending on the number of plots on the page
c
      if (pageplots.le.8) then 
         textsize = 10
         axistextsize = 8
      elseif (pageplots.le.12) then 
         textsize = 9
         axistextsize = 7
      elseif (pageplots.le.20) then 
         textsize = 8
         axistextsize = 6
      else 
         textsize = 4
         axistextsize = 2
      endif 
c
c     Add overall Page labels and references
c
      call grmtitle(title,ref)
c
c     Loop through all the plots on the page
c
c
c     Find the Max and Min for the plots - the axis scale
c     Max and Min are always calculated separately for each plot.
c     The axis max and min assume that the array of axis coordinates
c     is in ascending order.
c
c     Global scaling
c
      write (6,*) 'sctype:',sctype,':',xlab,':',ylab,':'
c
c
c     Convert all values to logarithmic for sctype 5 and 6.
c
      if (sctype.eq.5.or.sctype.eq.6) then
         do ip = 1,nplts
            do in = 1,pngs(ip)
               do ik = 1,pnks(ip,in)
                  if (mvals(ik,ip,in).ne.0.0) then
                     mvals(ik,ip,in) = log10(mvals(ik,ip,in))
                  else
                     mvals(ik,ip,in) = 0.0
                  endif
               end do
            end do
         end do
      endif
c
      pltmax = -hi
      pltmin =  hi
      if (sctype.eq.1.or.sctype.eq.3.or.sctype.eq.5) then
         do ip = 1,nplts
            do in = 1,pngs(ip)
               write (6,*) 'ip,in:',ip,in,pnks(ip,in),pngs(ip)
               do ik = 1,pnks(ip,in)
c
                  pltmax = max(pltmax,mvals(ik,ip,in))
                  pltmin = min(pltmin,mvals(ik,ip,in))
c
c                 write (6,*) 'plt:',ik,ip,in,pltmax,pltmin,
c     >                     mvals(ik,ip,in)
c
               end do
            end do
         end do
         if (sctype.eq.3) then
            if (pltmax.ge.0.0.and.pltmin.ge.0.0) then
               pltmin = 0.0
            elseif (pltmax.le.0.0.and.pltmin.le.0.0) then
               pltmax = 0.0
            endif
         endif
c
c        Rescale plot max and min -
c        just so the plots don't hit the edges
c
         if (pltmax.ge.0.0) then
             pltmax = pltmax * 1.05
         else
             pltmax = pltmax * 0.95
         endif
c
         if (pltmin.ge.0.0) then
            pltmin = pltmin * 0.95
         else
            pltmin = pltmin * 1.05
         endif
c
      endif
c
      do id = 1,nplts,pageplots
c
        ipmax = min(nplts,id+pageplots-1)
c
        write (6,*)'id:',id,nplts,pageplots,ipmax
c
        do ip = id,ipmax
c
c          Copy graph labels for individual plot into local data
c
           do in = 1,pngs(ip)
c
              elabs(in) = mlabs(ip,in) 
c
c             Load drawing instructions if needed
c
              if (drawtypesw.ne.0) then  
                 drawtype(in) = mdrawtype(ip,in)
              endif 
c
           end do
c
c         Set externally imposed scales for scale option 7
c
          if (sctype.eq.7) then
             pltmin = pltmins(ip)
             pltmax = pltmaxs(ip)
          endif
c
c         SET UP BOX FOR EACH PLOT - BASED ON POSITION/PLOT NUMBER
c
          call grmbox(ip)
c
c         Add labels for plot (individual plot titles and annotation.
c
          call grmlabels(ip,pltlabs(ip),elabs,pngs(ip),drawtype)
c
c         Label and plot the axes for the given plot position
c
          axmax = -hi
          axmin = hi
c
c         Set axis to cover the entire range from all data sets
c
c         Scan entire data sets because they may no longer be ordered.
c
          do in = 1,pngs(ip)
c
             do ik = 1,pnks(ip,in)
c
                axmin = min(mouts(ik,ip,in),axmin)
                axmax = max(mouts(ik,ip,in),axmax)
c
             end do
c
          end do
c
c         Adjust the ends of the axis to expand viewing region if 
c         requested. 
c
          if (axfact.gt.0.0) then
             axmax = axmax + abs(axmax-axmin)*axfact 
             axmin = axmin - abs(axmax-axmin)*axfact 
          endif
c
c         Allow for closeup plots at either end of the ring by passing
c         another parameter - if 0.0 it plots the whole scale.
c         if < 0.0 it plots from axmin to abs(pltfact)*(axmax-axmin)
c         and if greater than 0.0 it plots from
c         (1.0-pltfact)*(axmax-axmin) to axmax.
c
          if (pltfact.lt.0.0) then
             if (abs(pltfact).gt.1.0) pltfact = 0.1
             axmax = abs(pltfact) * (axmax-axmin) + axmin
          elseif (pltfact.gt.0.0) then
             if (pltfact.gt.1.0) pltfact = 0.1
             axmin = (1.0-pltfact) * (axmax-axmin) + axmin
          endif
c
          if (sctype.eq.2.or.sctype.eq.4.or.sctype.eq.6) then
c
            pltmax = -hi
            pltmin =  hi

            do in = 1,pngs(ip)
              do ik = 1, pnks(ip,in)
                pltmax = max(pltmax,mvals(ik,ip,in))
                pltmin = min(pltmin,mvals(ik,ip,in))
c
c                write(6,'(a,3i4,5(1x,g12.5))') 'PLT:',ip,in,ik,
c     >            mouts(ik,ip,in),mvals(ik,ip,in),pltmax,pltmin
c
              end do
            end do
c
            if (sctype.eq.4) then
              if (pltmax.ge.0.0.and.pltmin.ge.0.0) then
                pltmin = 0.0
              elseif (pltmax.le.0.0.and.pltmin.le.0.0) then
                pltmax = 0.0
              endif
            endif
c
c           Rescale plot max and min - just so the plots
c                                      don't hit the edges
c
            if (pltmax.ge.0.0) then
               pltmax = pltmax * 1.05
            else
               pltmax = pltmax * 0.95
            endif
c
            if (pltmin.ge.0.0) then
               pltmin = pltmin * 0.95
            else
               pltmin = pltmin * 1.05
            endif
c
          endif
c
          write (6,*) 'range:',axmin,axmax,pltmin,pltmax
c
          call grmaxes(ip,pltmin,pltmax,axmin,axmax,
     >                 xlab,ylab,sctype)
c
c         Plot data for this plot
c
          call grmdata(mvals,mouts,pnks,maxnms,maxplts,
     >               maxngs,ip,pngs(ip),
     >               pltmin,pltmax,axmin,axmax,elabs,sctype,
     >               drawtype)
c
        end do
c
        call frame
c
      end do
c
      return
      end
c
c
c
      subroutine grmtitle(title,ref)
      use mod_params
      use mod_slout
      use mod_grminfo
      implicit none
      character*(*) title,ref
c slmod begin - new
c      INCLUDE 'params'
c      INCLUDE 'slout'
c slmod end
c
c     GRMTITLE: This function places the titles on the page.
c
      integer i,j,l,l2,lenstr
      external lenstr
c      include 'grminfo'
c
      CALL PSPACE (0.0, 1.0, 0.0, 1.0)
      CALL MAP    (0.0, 1.0, 0.0, 1.0)
      CALL CTRMAG (textsize)

      if (pageplots.le.8) then 
c
         L = LENSTR(TITLE)
         CALL PCSCEN (0.5,0.99, title(:l))
c
         L = LENSTR(REF)
         CALL PCSCEN (0.5,0.975, ref(:l))
c
      else
c
c slmod begin - new
         IF (slopt3.EQ.2) THEN
           CALL CTRMAG(14)
           L = LENSTR(TITLE)
           CALL PCSCEN (0.5,0.998, title(:l))
         ELSE
c
c          For many plots/page - put title and ref on one line in larger text
c
           L = LENSTR(TITLE)
           L2 = LENSTR(REF)
           call ctrmag(textsize+2)
           CALL PCSCEN (0.5,0.995, ref(:l2)//' '//title(:l))
           call ctrmag(textsize)  

c
c            L = LENSTR(TITLE)
c           CALL PCSCEN (0.5,0.998, title(:l))
c
c           L = LENSTR(REF)
c           CALL PCSCEN (0.5,0.992, ref(:l))
         ENDIF
c
c         CALL CTRMAG(14)
c         L = LENSTR(TITLE)
c         CALL PCSCEN (0.5,0.998, title(:l))
c
c         L = LENSTR(REF)
c         CALL PCSCEN (0.5,0.992, ref(:l))
c slmod end
      endif  

c
c     Writing Additional comments 
c
c slmod begin - new
      CALL CTRMAG(10)
      DO l = 1, 10
        IF (char(l).NE.' ') CALL PLOTST (1.04,0.737+(l-1)*0.02,char(l))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST (1.04,0.250-(j-1)*0.02,char(i))
      ENDDO
c slmod end
      return
      end
c
c
c
      subroutine grmbox(boxindex)
      use mod_colours
      implicit none
      integer boxindex
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c
c     GRMBOX: DRAWS THE OUTSIDE BOX FOR THE PLOT INDICATED BY
c             THE BOXINDEX.
c
c
      real xpt,ypt,xwid,ywid,xsep,ysep
c
      call grmfindxy(boxindex,xpt,ypt,xwid,ywid,xsep,ysep)
c
c     Set up vector/ND space mapping
c
      call pspace(0.0,1.0,0.0,1.0)
      call map (0.0,1.0,0.0,1.0)
      call full
      CALL LINCOL (defcol)
      CALL THICK  (1)
c
c     Start drawing box
c
      call positn(xpt,ypt)
c
c     Draw box
c
      call join(xpt+xwid,ypt)
      call join(xpt+xwid,ypt+ywid)
      call join(xpt,ypt+ywid)
      call join(xpt,ypt)
c
      return
      end
c
c
c
      subroutine grmlabels(ip,pltlabs,elabs,ngs,drawtype)
      use mod_params
      use mod_colours
      use mod_grminfo
      use mod_slout
      implicit none
      integer ip,ngs,drawtype(ngs)
      character*(*) pltlabs
      character*(*) elabs(ngs)
c
c     GRMLABELS : This routine places the labels on
c                 each individual plot.
c
c
      integer in,l,lenstr
      external lenstr
      real xpt,ypt,xwid,ywid,xsep,ysep,xp,yp
c
c      include 'colours'
c      include 'grminfo'
c slmod begin - new
c      INCLUDE 'params'
c      INCLUDE 'slout'
c slmod end
c     integer init_col,get_col,next_col,
      integer :: drawcount
c      external init_col,get_col,next_col
c
c      INTEGER COLOUR(8)
c
c      integer icol
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
      call grmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
      CALL PSPACE (0.0, 1.0, 0.0, 1.0)
      CALL MAP    (0.0, 1.0, 0.0, 1.0)
      CALL CTRMAG (textsize)
      call full
      call thick(1)
c
      L = LENSTR(PLTLABS)
c slmod begin - new
      IF (pageplots.LE.8.OR.(slopt3.NE.2.AND.slopt3.NE.3)) THEN
        CALL PCSCEN (xpt+xwid/2.0,ypt+ywid+0.75*ysep, pltlabs(:L))
      ELSEIF (slopt3.EQ.3) THEN
        CALL CTRMAG(12)
        CALL PCSCEN(xpt+xwid/2.0,ypt+ywid-1.20*ysep, pltlabs(:L))
        CALL CTRMAG(textsize)
        RETURN
      ELSEIF (slopt3.EQ.2) THEN
        CALL CTRMAG(12)
        CALL PCSCEN(xpt+xwid/2.0,ypt+ywid-1.50*ysep, pltlabs(:L))
        RETURN
      ENDIF
c
c      CALL PCSCEN (xpt+xwid/2.0,ypt+ywid+0.75*ysep, pltlabs(:L))
c slmod end

c
      col = init_col()
      drawcount = 0
c
c      icol = ncols
c
      DO IN = 1,NGS
c
c       Set line colour
c
c
c        icol = icol - 1
c        IF (ICOL.LT.1) ICOL = NCOLS
c
        CALL LINCOL (COL)

        col = next_col()
c
c       Set line type
c
        IF (in.LE.1) THEN
           CALL FULL
        ELSEIF (in.LE.5) THEN
           CALL BROKEN (3*in,2*in,3*in,2*in)
        ELSE
           CALL BROKEN (2*in,1*in,2*in,1*in)
        ENDIF
c
c       Draw line segment
c
        if (drawtype(in).eq.1.or.drawtype(in).eq.3) then

           xp = xpt + in * (xwid/(ngs+1)) - 0.05
           yp = ypt + ywid + 0.25 * ysep
           call positn(xp,yp)
           call join(xp+0.03,yp)

        endif
c
c       Draw marker
c
        if (drawtype(in).eq.2.or.drawtype(in).eq.3) then

           drawcount = drawcount + 1

           xp = xpt + in * (xwid/(ngs+1)) - 0.05
           yp = ypt + ywid + 0.25 * ysep

           call plotnc(xp+0.015,yp,marker(drawcount))

        endif
c
        call full
        L = LENSTR(ELABS(IN))
        CALL PLOTCS(xp+0.04,yp,elabs(in)(5:L))
c
      end do
c
      return
      end
c
c
c
      subroutine grmaxes(ip,pltmin,pltmax,axmin,axmax,xlab,ylab,
     >                   sctype)
      use mod_params
      use mod_colours
      use mod_slout
      use mod_grminfo
      implicit none
      integer ip,iten,sctype
      real pltmin,pltmax,axmin,axmax
      character*(*) xlab,ylab
c      include 'colours'
c      include 'grminfo'
c slmod begin - new
c      INCLUDE 'params'
c      INCLUDE 'slout'
c slmod end
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c
c     GRMAXES: This routine plots the axes on the
c              particular box and sets the mapping between
c              the ND-space and vector-space in preparation
c              for the plotting.
c
c
      real xpt,ypt,xwid,ywid,xsep,ysep
      real tmin,tmax,tdiv,tdiv_tmp,mult,power
      integer in,ik,l,lenstr,iexp
      external lenstr,iexp
c
      call grmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
c     Draw Axes scales as appropriate.
c
      CALL LINCOL (defcol)
      CALL THICK  (2)
      call full
c
c     Xscale
c
      ITEN = IEXP(AXMIN,AXMAX)
c
c     Xaxis scale
c
      POWER = 10.0**(-ITEN)
      TMIN = AXMIN * POWER
      TMAX = AXMAX * POWER
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      write(6,*) 'xscale:',axmin,axmax,iten,power,tmin,tmax
      CALL MAP    (tmin, tmax, ypt, ypt+ywid)
c
      CALL CTRMAG (axistextsize)
      if (pageplots.le.8) then 
         CALL XSCALE
      else
c slmod begin - new        
          IF (slopt3.NE.2) THEN
           tdiv = abs(tmax-tmin)/5.0  
           mult = 1.0
           tdiv_tmp = int(tdiv)

           do while (tdiv_tmp.le.0.0)
              tdiv = tdiv * 10.0
              mult = mult * 10.0 
              tdiv_tmp = int(tdiv)
           end do            
c
           tdiv = tdiv_tmp/mult
c
           CALL XSCALI(tdiv)
         ENDIF
c
c           CALL XSCALI(0.0)
c slmod end
      endif
      CALL CTRMAG (textsize)
c
c     Xaxis labels
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      call map    (xpt, xpt+xwid,ypt, ypt+ywid)
c
      IF (ITEN .NE. 0) THEN
        CALL CTRMAG (axistextsize)

        if (pageplots.le.8) then  
           CALL PLOTCS (xpt+0.85*xwid,ypt-0.6*ysep,'X10')
           CALL TYPENI (ITEN)
        else
c slmod begin - new
           IF (slopt3.NE.2) THEN
             CALL PLOTCS (xpt+0.85*xwid,ypt-1.4*ysep,'X10')
c             CALL PLOTCS (xpt+0.85*xwid,ypt-0.8*ysep,'X10')
             CALL TYPENI (ITEN)
           ENDIF  
c
c           CALL PLOTCS (xpt+0.85*xwid,ypt-0.8*ysep,'X10')
c           CALL TYPENI (ITEN)
c slmod end
        endif

        CALL CTRMAG (textsize)
      ENDIF
c
      CALL CTRMAG (axistextsize)
      L=LENSTR(xlab)
      if (pageplots.le.8) then 
         CALL PCSCEN (xpt+0.5*xwid,ypt-0.6*ysep,xlab(1:L))
      else
c slmod begin - new
         IF (slopt3.NE.2) THEN
c
c           call ctrmag(axistextsize+2) 
           CALL PCSCEN (xpt+0.5*xwid,ypt-1.4*ysep,xlab(1:L))
c           call ctrmag(axistextsize) 
c
c           CALL PCSCEN (xpt+0.5*xwid,ypt-0.8*ysep,xlab(1:L))
c
         ENDIF
c
c         CALL PCSCEN (xpt+0.5*xwid,ypt-0.8*ysep,xlab(1:L))
c slmod end
      endif

      CALL CTRMAG (textsize)
C
c     yscales
c
      if (sctype.eq.1.or.sctype.eq.2
     >    .or.sctype.eq.3.or.sctype.eq.4.or.sctype.eq.7) then
        ITEN = IEXP(PLTMIN, PLTMAX)
c
        POWER = 10.0**(-ITEN)
        TMIN = pltMIN * POWER
        TMAX = pltMAX * POWER
        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
        write(6,*) 'yscale:',pltmin,pltmax,iten,power,tmin,tmax
c
        CALL CTRMAG (axistextsize)
        if (pageplots.le.8) then 
           CALL YSCALE
        else
c slmod begin - new
           IF (slopt3.NE.2) THEN
c
              tdiv = abs(tmax-tmin)
c
              if (tdiv.gt.50) then 
                 tdiv=10   
              elseif (tdiv.gt.20) then 
                 tdiv=5   
              elseif (tdiv.gt.10) then 
                 tdiv=2   
              else 
                 tdiv=1   
              endif
c
c             mult = 1.0
c             tdiv_tmp = int(tdiv)
c
c             do while (tdiv_tmp.le.0.0)
c                tdiv = tdiv * 10.0
c                mult = mult * 10.0 
c                tdiv_tmp = int(tdiv)
c             end do            
c
c             tdiv = tdiv_tmp/mult
c
             CALL YSCALI(tdiv) 
           ENDIF
c
c           CALL YSCALI(0.0) 
c slmod end
        endif
        CALL CTRMAG (textsize)
c
c     Y = 0 line if the scale crosses zero
c
        if (pltmin*pltmax.lt.0.0) then
           call full
           call thick(1)
           CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
           CALL MAP    (axmin,axmax,pltmin,pltmax)
           call positn(axmin,0.0)
           call join(axmax,0.0)
        endif
      elseif (sctype.eq.5.or.sctype.eq.6) then
c
c     Logarithmic Y-axis
c
c        ITEN = IEXP(PLTMIN, PLTMAX)
c
c        POWER = 10.0**(-ITEN)
c        TMIN = pltMIN * POWER
c        TMAX = pltMAX * POWER
c        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
c        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
c        write(6,*) 'yscall:',pltmin,pltmax,iten,power,tmin,tmax
c        CALL CTRMAG (textsize)
c        CALL YSCALL
c
        ITEN = IEXP(PLTMIN, PLTMAX)
c
        POWER = 10.0**(-ITEN)
        TMIN = pltMAX - 4.0
        TMAX = pltMAX
        CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
        CALL MAP    (xpt, xpt+xwid,tmin,tmax)
        write(6,*) 'yscall:',pltmin,pltmax,iten,power,tmin,tmax
        CALL CTRMAG (textsize)
        CALL YSCALE
      endif
c
c     Y axis labels
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
      CALL MAP    (xpt, xpt+xwid,ypt, ypt+ywid)
c
c     Rotate text 90 degrees
c
      CALL CTRORI (90.0)
      IF (ITEN .NE. 0) THEN
         CALL CTRMAG (axistextsize)
         if (pageplots.le.8) then 
            CALL PLOTCS (xpt-0.6*xsep,ypt+0.90*ywid,'X10')
            CALL TYPENI (ITEN)
         else
c slmod begin - new
            IF (slopt3.NE.2) THEN
              CALL PLOTCS (xpt-1.2*xsep,ypt+0.90*ywid,'X10')
c              CALL PLOTCS (xpt-1.0*xsep,ypt+0.90*ywid,'X10')
              CALL TYPENI (ITEN)
            ENDIF
c
c            CALL PLOTCS (xpt-1.0*xsep,ypt+0.90*ywid,'X10')
c            CALL TYPENI (ITEN)
c slmod end
         endif
         CALL CTRMAG (textsize)
      ENDIF
      CALL CTRMAG (axistextsize)
      L = LENSTR (YLAB)

      if (pageplots.le.8) then
         CALL PCSCEN (xpt-0.6*xsep,ypt+0.5*ywid,YLAB(:L))
      else
c slmod begin - new
         IF (slopt3.NE.2) THEN
           CALL PCSCEN (xpt-1.2*xsep,ypt+0.5*ywid,YLAB(:L))
c           CALL PCSCEN (xpt-1.0*xsep,ypt+0.5*ywid,YLAB(:L))
         ENDIF
c
c         CALL PCSCEN (xpt-1.0*xsep,ypt+0.5*ywid,YLAB(:L))
c slmod end
      endif        

      CALL CTRMAG (textsize)
c
      if (sctype.eq.5.or.sctype.eq.6) then
         CALL PCSCEN (xpt-0.6*xsep,ypt+0.1*ywid,'(LOG)')
      endif
c
      CALL CTRORI (0.0)
c
      return
      end
c
c
c
      subroutine grmdata (mvals,mouts,pnks,maxnms,
     >                    maxplts,maxngs,ip,ngs,
     >                    pltmin,pltmax,axmin,axmax,elabs,
     >                    sctype,drawtype)
      use mod_colours
      use mod_grminfo
      implicit none
      integer maxnms,maxplts,maxngs,ngs,ip,sctype
      integer pnks(maxplts,maxngs),drawtype(ngs)
      real pltmin,pltmax,axmin,axmax
      real mvals(maxnms,maxplts,maxngs)
      real mouts(maxnms,maxplts,maxngs)
      character*(*) elabs(maxngs)
c
c     GRMDATA: This routine plots the sepcified data
c              on one of the small plots on the page
c              specified by the ip value.
c
      real xpt,ypt,xwid,ywid,xsep,ysep
      real tmin,tmax
      integer in,ik
c
c      include 'colours'
c      include 'grminfo'
c slmod begin
      REAL x1,y1,LO
      PARAMETER (LO=1.E-37 )
c slmod end

c      integer init_col,get_col,next_col,drawcount
      integer :: drawcount
c      external init_col,get_col,next_col
c
c      integer icol
c
c      INTEGER COLOUR(8)
c
c      DATA COLOUR /1,2,3,4,9,7,6,8/
c
      call grmfindxy(ip,xpt,ypt,xwid,ywid,xsep,ysep)
c
      CALL PSPACE (xpt, xpt+xwid,ypt, ypt+ywid)
c
      if (sctype.eq.1.or.sctype.eq.2
     >    .or.sctype.eq.3.or.sctype.eq.4.or.sctype.eq.7) then
         CALL MAP  (axmin,axmax,pltmin,pltmax)
      elseif (sctype.eq.5.or.sctype.eq.6) then
         CALL MAPYL(axmin,axmax,pltmax-6.0,pltmax)
      endif
c
c      CALL RGB
c      CALL COLSET(0.7,1.0,0.,8)
c      CALL COLSET(1.0,0.7,0.,9)
C
C
      call ctrmag(axistextsize)
      call thick(1)
c
      col = init_col()
c
      drawcount = 0
c
c      icol = ncols
c
      do in = 1,ngs
c
c       Set line colour
c
c        icol = icol - 1
c        IF (ICOL.lT.1) ICOL = ncols
c
        CALL LINCOL (COL)
c
        col = next_col()
c
c       Set line type
c
        IF (in.LE.1) THEN
           CALL FULL
        ELSEIF (in.LE.5) THEN
           CALL BROKEN (3*in,2*in,3*in,2*in)
        ELSE
           CALL BROKEN (2*in,1*in,2*in,1*in)
        ENDIF
c
c       Draw each plot using lines if requested
c
        if (drawtype(in).eq.1.or.drawtype(in).eq.3) then

           call positn(mouts(1,ip,in),mvals(1,ip,in))

c
           do ik = 1,pnks(ip,in)
c
              call join(mouts(ik,ip,in),mvals(ik,ip,in))
c
              IF (IK.EQ.3.OR.IK.EQ.pnks(ip,in)-3)
     >           CALL PLOTST (mouts(ik,ip,in),mvals(ik,ip,in),
     >                     elabs(in)(1:4))
c
c              write(6,*) 'plt:',in,ip,ik,
c     >              mouts(ik,ip,in),mvals(ik,ip,in)
c
           end do
c
        endif
c
c       Draw each line using markers if requested
c
        if (drawtype(in).eq.2.or.drawtype(in).eq.3) then
c
           drawcount = drawcount + 1
c
           do ik = 1,pnks(ip,in)
c

            x1 = mouts(ik,ip,in)
            y1 = mvals(ik,ip,in)
            IF (X1.GE.AXMIN.AND.X1.LE.AXMAX.AND.Y1.GE.pltMIN.AND.
     >          Y1.LE.PLTMAX.AND.Y1.NE.LO)
     .        CALL BOX (x1-0.008*(AXMAX-AXMIN),
     .                  x1+0.008*(AXMAX-AXMIN),
     .                  y1-0.008*(PLTMAX-PLTMIN)*xwid/ywid,
     .                  y1+0.008*(PLTMAX-PLTMIN)*xwid/ywid)
c              call plotnc(mouts(ik,ip,in),mvals(ik,ip,in),
c     >                    marker(drawcount))


c
           end do
c
        endif


      end do
c
      call full
c
      return
      end
c
c
c
      subroutine grmfindxy (boxindex,xpt,ypt,xwid,ywid,
     >                      xsep,ysep)
      use mod_params
      use mod_slout
      use mod_grminfo
      implicit none
      integer boxindex
      real xpt,ypt,xwid,ywid,xsep,ysep
c
c      include 'grminfo'
c slmod begin - new
c      INCLUDE 'params'
c      INCLUDE 'slout'
c slmod end
c
c     GRMFINDXY: This routine finds the X,Y co-ordinates and
c                size of the box to be used for the particular
c                plot.
c
      integer ix,iy,n,xplts,yplts,plotindex
      real    xpts(5),ypts(5),xwid0,ywid0,xsep0,ysep0
c
c     The following describes the box location in ND space
c     X,Y are the location of the lower left corner of
c     the box.
c     Xwid and Ywid are the X and Y extent of the box
c     Xsep and Ysep are the anount of X and Y border
c     around each box.
c
      integer in
c
      data    xsep0 /0.05/
      data    ysep0 /0.05/
c
c     Draw Box - Boxindex is 1 to pageplots
c     Modify so that this is so.
c
      n = (boxindex -1)/ pageplots
      plotindex = boxindex - n * pageplots
c
c      xplts = 2
c      yplts = (pageplots-1) /2 + 1
c
      if (pageplots.le.8) then 
         xplts = 2
         yplts = (pageplots-1) /xplts + 1
         xsep0 = 0.05
         ysep0 = 0.05
      elseif (pageplots.le.12) then
         xplts = 3
         yplts = (pageplots-1) /xplts + 1
         xsep0 = 0.03
         ysep0 = 0.03
      elseif (pageplots.le.20) then
         xplts = 4
         yplts = (pageplots-1) /xplts + 1
c slmod begin - new
         IF (slopt3.EQ.2) THEN
           xsep0 = 0.01
           ysep0 = 0.01
         ELSE
           xsep0 = 0.02
           ysep0 = 0.02
         ENDIF
c
c         xsep0 = 0.02
c         ysep0 = 0.02
c slmod end
      else
         xplts = 5
         yplts = (pageplots-1) /xplts + 1
         xsep0 = 0.01
         ysep0 = 0.01
      endif      
c
      xwid0 = (1.0 - xplts * 2.0 * xsep0 - 2.0*xsep0) / xplts
      ywid0 = (1.0 - yplts * 2.0 * ysep0 - 2.0*ysep0) / yplts
c
      do in = 1,xplts
         xpts(in) = 2.0*xsep0 + (in-1) * (xwid0+2.0*xsep0)
      end do
c
      do in = yplts,1,-1
         ypts(in) = 2.0*ysep0+(yplts-in)*(ywid0+2.0*ysep0)
      end do
c
c     iy is 1 to pageplots/2
c     ix is 1 to 2
c
      iy = (plotindex-1)/xplts  +1 
c
      ix = plotindex - (iy-1) * xplts   
c
      xpt = xpts(ix)
      ypt = ypts(iy)
c
      xwid = xwid0
      ywid = ywid0
c
      xsep = xsep0
      ysep = ysep0
c
c      write (6,*) 'GRMFind:',plotindex,boxindex,n,iy,ix,
c     >                   xpt,ypt,xwid,
c     >                   ywid,xsep,ysep,pageplots
c
      return
      end
c
c
c
      subroutine grbar(valsts,istart,istop,novals,nosets,
     >                 ymin,ymax,iflag,pnames1,pnames2,grtitle,
     >                 cnames,ylab)
      use mod_params
      use mod_colours
      use mod_slout
      implicit none
      integer istart,istop,novals,nosets,iflag
      character*(*) pnames1(istop),pnames2(istop),
     >              grtitle,cnames(novals)
      character*(*) ylab
      real valsts (novals,nosets),ymin,ymax
c
c      include 'params'
c      include 'colours'
c      integer init_col,get_col,next_col
c      external init_col,get_col,next_col
c
c     GRBAR: This routine draws a multi-part incremental
c              bar chart in one of three position on the page.
c              Top, bottom or full page. It also annotates
c              using the labels provided if they are not equal to
c              a null string.
c
c              IFLAG =1 - use upper half of drawing space - bar chart
c              IFLAG =2 - use lower half of drawing space - bar chart
c              IFLAG =3 - use complete drawing space - bar chart
c              IFLAG =4 - use complete drawing space -histogram
c
c
c slmod begin
c              IFLAG =5
c
c      INCLUDE 'slout'

      INTEGER i1
c slmod end
      integer in,npts
      real xmin,xmax,xpt,ypt,ysep,spot,dspot
c
      integer lstlen,l,ispot,lenstr
      integer ncolrs(8)
      external lenstr
c
      real baspos,barwid,barpos(maxpts)
c
c     Check for ymin = ymax (no plot) - if this is the case then return
c
      if (ymin.eq.ymax) return 
c
      npts = istop-istart+1
c
c     Set up a selection of fill patterns for the bar charts
c
      lstlen = 8
c
      ncolrs(1) = -3
      ncolrs(2) = -22
      ncolrs(3) = -33
      ncolrs(4) = -4
      ncolrs(5) = -14
      ncolrs(6) = -34
      ncolrs(7) = -40
      ncolrs(8) = -38
c
      call lstcol(ncolrs,lstlen)
c
c     Set up MIN and MAX on Horizontal axis
c
      xmin = 0.0
      xmax = npts
c
c     Set up plotting space in universal coordinates
c
      if (IFLAG.eq.1) then
         CALL PSPACE (0.15, 0.85, 0.55, 0.85)
c
c        Set up value in coordinate space representing 0.05 in real space
c
         ysep = (ymax-ymin)/6.0
c
      elseif (iflag.eq.2) then
         CALL PSPACE (0.15, 0.85, 0.15, 0.45)
c
c        Set up value in coordinate space representing 0.05 in real space
c
         ysep = (ymax-ymin)/6.0
c
      elseif (iflag.eq.3.or.iflag.eq.4) then
         CALL PSPACE (0.15, 0.85, 0.15, 0.85)
c
c        Set up value in coordinate space representing 0.05 in real space
c
         ysep = (ymax-ymin)/14.0
c
c slmod begin
      ELSEIF (iflag.EQ.5) THEN
         CALL PSPACE (map1x,map2x,map1y,map2y)
         ysep = (ymax-ymin)/3.0
c slmod end
      endif
c
c     Map plotting space onto local coordinates
c
      CALL MAP    (xmin, xmax , YMIN, YMAX)
c
c     Draw Y-AXIS
c
      CALL LINCOL (defcol)
c slmod begin
      IF (iflag.EQ.5) THEN
      ELSE
        CALL THICK  (2)
        CALL CTRMAG (10)
      ENDIF
c
c        CALL CTRMAG (10)
c slmod end
      CALL YAXIS
c
      call thick(1)
      l = lenstr(ylab)
c
      if (l.gt.0) then
         call ctrori(90.0)
c slmod begin
         IF (iflag.EQ.5) THEN
           call PLOTST(xmin-1.7*(xmax-xmin)/20.0,
     >                ymin+0.1*(ymax+ymin),ylab(:l))
         ELSE
           call pcscen(xmin-1.5*(xmax-xmin)/28.0,
     >                (ymax+ymin)/2.0-(ymax+ymin)/4.0,ylab(:l))
         ENDIF
c
c         call pcscen(xmin-1.5*(xmax-xmin)/28.0,
c   >                (ymax+ymin)/2.0,ylab(:l))
c slmod end
         call ctrori(0.0)
      endif
c
c     Set up the plotting positions and the widths of each
c     barchart.
c
      if (npts.gt.100) then
         barwid = (xmax-xmin) / npts
      else
         barwid = (xmax-xmin) / (npts*1.5)
      endif
c
      baspos = 0.0
c
      do in = istart,istop
         barpos(in) = in - istart+ 0.5
      end do
c
      call full
c
c     Draw incremantal bar chart for optiosn 1 to 3
c
      if (iflag.eq.1.or.iflag.eq.2.or.iflag.eq.3) then

         call incbar(baspos,barwid,barpos,valsts,istart,istop,
     >               novals,nosets)

      elseif (iflag.eq.4) then

         call mulhis(0.0,0.0,barwid,valsts,istart,istop,
     >               novals,nosets)
c slmod begin
      ELSEIF (iflag.EQ.5) THEN

         barwid = (xmax-xmin) / (npts*2.0)

         ncolrs(1) = 1
         DO i1 = 2, 8
           ncolrs(i1) = ncols + i1
         ENDDO

c         IF (nosets.EQ.1) ncolrs(1) = ncols

         call lstcol(ncolrs,lstlen)

         call mulbar(baspos,barwid,barpos,valsts,istart,istop,
     >               novals,nosets)

c slmod end
      endif
c
c     Add annotation from pnames1, pnames2 and grtitle
c
c
c     Title First
c
c slmod begin
      IF (iflag.EQ.5) THEN
      ELSE
        CALL CTRMAG (12)
      ENDIF
c
c      CALL CTRMAG (12)
c slmod end
      L = lenstr(grtitle)
      call plotcs(xmin,ymax+ysep/2.0,grtitle(:L))
c
c     Add comments for each bar chart - in two lines
c     pnames1 - first line
c     pnames2 - second line
c
c     but only if string is non-null
c
c slmod begin
      IF (iflag.EQ.5) THEN
      ELSE
        call ctrmag(8)
      ENDIF
c
c      call ctrmag(8)
c slmod end
c
      do in = istart,istop
c
         L = lenstr(pnames1(in))
c
         if (l.gt.0) then
            xpt = in - istart + 0.5
            ypt = ymin - 0.333 * ysep
            call pcscen(xpt,ypt,pnames1(in)(:l))
         endif
c
         L = lenstr(pnames2(in))
c
         if (l.gt.0) then
            xpt = in - istart + 0.5
            ypt = ymin - 0.667 * ysep
            call pcscen(xpt,ypt,pnames2(in)(:l))
         endif
c
      end do
c
c     Put entry in symbol table for each catagory - with fill pattern
c
      do in = 1,nosets
c
         L = lenstr(cnames(in))
c
         if (l.gt.0) then
C
C            WRITE ENTRY IN SYMBOL TABLE
C
             ispot = 12
c
c            Calculate position
c
             CALL FULL
c slmod begin
             IF (iflag.EQ.5) THEN
               CALL PSPACE (0.0, 1.35, 0.0,1.0)
               CALL CSPACE (0.0, 1.35, 0.0,1.0)
               CALL MAP    (0.0, 1.35, 0.0,1.0)
               dspot = ispot * 1.5 / 500.0
               SPOT = 0.818 - dspot * in
               call filcol(ncolrs(in))
               call box(1.03-dspot*0.667,1.03,spot-0.333*dspot,
     .                                       spot+0.333*dspot)
               call filcol(0)
               CALL PLOTST (1.05 ,SPOT, cnames(in)(:l))
             ELSE
               CALL CTRMAG (ISPOT)
               CALL PSPACE (0.0, 1.35, 0.0,1.0)
               CALL CSPACE (0.0, 1.35, 0.0,1.0)
               CALL MAP    (0.0, 1.35, 0.0,1.0)
c
               dspot = ispot * 4.0 /500.0
c
               SPOT = 0.818 - dspot * in
c
               call filcol(ncolrs(in))
c
               call box(0.96,1.03,spot-0.333*dspot,spot+0.333*dspot)
c
               call filcol(0)
c
               CALL PLOTST (1.05 ,SPOT, cnames(in)(:l))
             ENDIF
c
c             CALL CTRMAG (ISPOT)
c             CALL PSPACE (0.0, 1.35, 0.0,1.0)
c             CALL CSPACE (0.0, 1.35, 0.0,1.0)
c             CALL MAP    (0.0, 1.35, 0.0,1.0)
c
c             dspot = ispot * 4.0 /500.0
c
c             SPOT = 0.818 - dspot * in
c
c             call filcol(ncolrs(in))
c
c             call box(0.96,1.03,spot-0.333*dspot,spot+0.333*dspot)
c
c             call filcol(0)
c
c             CALL PLOTST (1.05 ,SPOT, cnames(in)(:l))
c slmod end    
c
          end if
      end do
c
c
c     Turn off fill patterns after plotting
c
      call lstcol(ncolrs,0)
c
      return
      end
c
c
c
      subroutine find_scales(axis_min,axis_max,
     >                            expt_min,expt_max,
     >                            iexpt_plot)
      use mod_params
      use mod_expt_data
      implicit none
      real axis_min,axis_max,expt_min,expt_max
      integer iexpt_plot
c 
c      include 'params'
c      include 'expt_data'  
c      include 'outcom'
c
c     FIND_SCALES
c
c     This routine scans all of the requested datasets for the plot
c     and then calculates the axis range based on the entire set of 
c     experimental data.
c
c     If IEXPT_PLOT is specified and the extended experimental
c     dataset input is not present (EXPT_NSETS=0) then the code 
c     uses only the IEXPT_PLOT data. However, if EXPT_NSETS is
c     non-zero then the code will load all of the data sets and 
c     calculate the max and min over all datasets. 
c
c
      integer in
      real expt_data_min,expt_data_max      
      real expt_axis_min,expt_axis_max      
c
      expt_min =  HI
      expt_max = -HI  
c
c     Extended experimental datasets specified. 
c
      if (expt_nsets.gt.0) then
c
         do in = 1,expt_nsets
c
            call load_expt_scales(expt_axis_min,expt_axis_max,
     >               expt_data_min,expt_data_max,expt_datasets(in))  
c
            axis_min = min(axis_min,expt_axis_min)      
            axis_max = max(axis_max,expt_axis_max)      
            expt_min = min(expt_min,expt_data_min)      
            expt_max = max(expt_max,expt_data_max)      
c
         end do
c
c        Set iexpt plot to a non-zero but invalid dataset so that
c        the experimental data plot routine is properly called. 
c
         iexpt_plot = -1
c
c     Single experimental data set specified in plot input line. 
c
      elseif (iexpt_plot.ne.0) then 
c
            call load_expt_scales(expt_axis_min,expt_axis_max,
     >               expt_data_min,expt_data_max,abs(iexpt_plot))  
c
            axis_min = min(axis_min,expt_axis_min)      
            axis_max = max(axis_max,expt_axis_max)      
            expt_min = min(expt_min,expt_data_min)      
            expt_max = max(expt_max,expt_data_max)      
c
      endif
c
      return
      end      
c
c
c
      subroutine load_expt_scales(axis_min,axis_max,
     >                     expt_min,expt_max,iexpt)
      use mod_params
      use mod_expt_data
      implicit none
c
      real axis_min,axis_max,expt_min,expt_max
      integer iexpt
c
c      include 'params'
c      include 'expt_data'
c
c     load_expt_scales: 
c
c     This routine loads a set of experimental data and finds the min and 
c     max values over the data range. 
c
c
c     Local variables to hold experimental data
c
      integer local_iexpt 
      integer dataunit,in
c
c      integer maxcols
c      parameter(maxcols=1)
c
c      integer axis_type,num_expt,ncols,dataunit,in
c      real expt_axis(maxdatx),expt_data(maxdatx,maxcols)
c      character*100 datatitle
c
c     Initialize
c
      dataunit = exptunit 
      axis_min = HI 
      axis_max =-HI 
      expt_min = HI 
      expt_max =-HI 
c
c     Check for valid dataset 
c
c
c     if (iexpt.le.0) return
c
c
c     Load experimental data - if it is not already loaded to the
c     common block. 
c     
      if (expt_data_available.eq.0) then 
c
         local_iexpt = abs(iexpt)
c
         call load_expt_data(dataunit,local_iexpt,
     >                 expt_data_axis,expt_data_axis_type,
     >                 expt_data_values,
     >                 expt_data_maxcols,expt_data_maxdatx,
     >                 expt_data_num,expt_data_ncols,
     >                 expt_data_title)

      endif

c     
c     
c     Check to see that data was loaded correctly         
c     
      if (expt_data_num.le.0) then
c     
         write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 local_iexpt, ' - NO ELEMENTS FOUND'
         write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 local_iexpt, ' - NO ELEMENTS FOUND'
c      
c        On error - set actual iexpt to zero
c
         iexpt = 0
c     
c     Calculate minimum and maximum values over the experimental data
c     
      else  
c     
c        Find min and max
c     
         do in = 1,expt_data_num
c     
            expt_min = min(expt_min,expt_data_values(in,1))              
            expt_max = max(expt_max,expt_data_values(in,1))              
            axis_min = min(axis_min,expt_data_axis(in))              
            axis_max = max(axis_max,expt_data_axis(in))              
c     
         end do 
c     
      endif
c
      return 
      end
c
c
c
      subroutine plot_allexpt(iexpt_plot)
      use mod_params
      use mod_expt_data
      implicit none
      integer iexpt_plot
c
c      include 'params'
c      include 'expt_data'
c
c     PLOT_ALLEXPT:
c 
c     This add all the required experimental data sets to the
c     plot using a POINT method.
c
      integer in
c
c     Extended experimental datasets specified. 
c
      if (expt_nsets.gt.0) then
c
         do in = 1,expt_nsets
c
            call plot_expt(expt_datasets(in))
c
         end do
c
c     Single experimental data set specified in plot input line. 
c
      elseif (iexpt_plot.ne.0) then 
c
         call plot_expt(iexpt_plot)
c
      endif
c
      return 
c
      end
c
c
c 
      subroutine plot_expt(iexpt)
      use mod_params
      use mod_expt_data
      implicit none
      integer iexpt
c
c      include 'params'
c      include 'expt_data'  
c
c     PLOT_EXPT:
c
c     This routine loads and plots the one set of 
c     experimental data referenced by IEXPT
c
c
c     Local variables to hold experimental data
c
      integer local_iexpt
      integer dataunit
c
c      integer maxcols
c      parameter(maxcols=1)
c
c      integer axis_type,num_expt,ncols,dataunit
c      real expt_axis(maxdatx),expt_data(maxdatx,maxcols)
c      real expt_min,expt_max
c      character*100 datatitle

c
c     Other local variables 
c
      integer in,len,lenstr
      external lenstr 
c
c     Check for valid dataset 
c
c      if (iexpt.le.0) return
c
      local_iexpt = abs(iexpt)
c
c     Initialize 
c
      dataunit = exptunit 
c
c     Load experimental data  
c     
      if (expt_data_available.eq.0) then  

         call load_expt_data(dataunit,local_iexpt,
     >                 expt_data_axis,expt_data_axis_type,
     >                 expt_data_values,
     >                 expt_data_maxcols,expt_data_maxdatx,
     >                 expt_data_num,expt_data_ncols,
     >                 expt_data_title)
      endif
c     
c     Check to see that data was loaded correctly         
c     
      if (expt_data_num.le.0) then
c     
         write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 local_iexpt, ' - NO ELEMENTS FOUND'
         write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                 local_iexpt, ' - NO ELEMENTS FOUND'
c      
c       In the case of an error - set the actual IEXPT to zero
c
         iexpt = 0
c     
c     Add the data to the plot
c     
      else  
c     
         len = min(lenstr(expt_data_title),25)
c        
c        Plot experimental data point-wise on the graph
c 
         if (iexpt.gt.0) then 

            CALL GRTRAC(expt_data_axis,expt_data_values,
     >               expt_data_num,
     >              'EXPT'//expt_data_title(1:len),'POINT',1)

         elseif (iexpt.lt.0) then 

            CALL GRTRAC(expt_data_axis,expt_data_values,
     >               expt_data_num,
     >              'EXPT'//expt_data_title(1:len),'LINE',1)

         endif

c       
      endif
c
      return
      end

c
c
c
      SUBROUTINE CONTIL2(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,NPTSY,
     &                   CLEVLS,ISTRTL,ISTOPL,XGRIDS,YGRIDS)
      use mod_colours
      ! jdemod - this routine has a lot of implicitly declared variables and
      ! switching to explicit declaration will be too time consuming at the moment
      ! Code has been ok for ddecades and should remain so :) 
      !implicit none
C     
C          ------------------------------------------------
C          ROUTINE NO. ( 124)   VERSION (A8.8)    15:MAY:87
C          ------------------------------------------------
C
C          THIS DRAWS STRAIGHT-ELEMENT CONTOURS ON AN IRREGULAR GRID.
C
C
C          THE ARGUMENTS ARE AS FOLLOWS:
C
C          [SURFAS]  IS THE ARRAY OF SURFACE HEIGHT VALUES,
C          <ISTRTX>  IS THE LOWER X-EXTENT,
C          <ISTOPX>  IS THE UPPER X-EXTENT, WHILE
C          <ISTRTY>  AND
C          <ISTOPY>  ARE THE CORRESPONDIMG Y-BOUNDS.
C          <NPTSX>   IS THE ACTUAL ARRAY X-EXTENT, AND
C          <NPTSY>   IS THE ACTUAL ARRAY Y-EXTENT.
C          [CLEVLS]  CONTAINS THE VALUES OF CONTOUR HEIGHTS,
C          <ISTRTL>  IS THE STARTING POINT, AND
C          <ISTOPL>  IS THE END POINT OF THIS ARRAY.
C          [XGRIDS]  ARE THE GRID X-POSITIONS.
C          [YGRIDS]  ARE THE GRID Y-POSITIONS.
C
C
      REAL    SURFAS(NPTSX,NPTSY),CLEVLS(ISTOPL),
     &        XGRIDS(NPTSX),YGRIDS(NPTSY)
      LOGICAL OPENCO,ERRON,CURVED
C
      COMMON /T0ACON/ ANGCON
      COMMON /T0CANG/ STANG0,CRANG0
      COMMON /T0CLBL/ LLBCON
      COMMON /T0CLEV/ LEVEL,HEIGHT
      COMMON /T0CMAP/ MAPBIT(6,2048),ISTPTX,ISTPTY
      COMMON /T0CTYP/ OPENCO,CURVED
      COMMON /T0CWID/ CWIDX,CWIDY
      COMMON /T0PPOS/ XPLOT0,YPLOT0
      COMMON /T0TRAC/ IPRINT
      COMMON /T0TRAI/ ITRAC1,ITRAC2,ITRAC3,ITRAC4
      COMMON /T0WNDO/ X1WND0,X2WND0,Y1WND0,Y2WND0
      COMMON /T3ERRS/ ERRON,NUMERR
c slmod begin
      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

c      include 'colours'
c slmod end
C
      DATA NBITSW /32/
C
C
      ITRAC1= ISTRTX
      ITRAC2= ISTOPX
      ITRAC3= ISTRTY
      ITRAC4= ISTOPY
      IF (IPRINT.EQ.1) CALL G0MESG(48,8)
C
      ILENX= ISTOPX-ISTRTX
      ILENY= ISTOPY-ISTRTY
      IF (ISTRTX.LT.1.OR.ISTRTY.LT.1)         GO TO 901
      IF (ILENX.LT.1.OR.ILENX.GT.2048)         GO TO 901
      IF (ILENY.LT.1.OR.ILENY.GT.2048)         GO TO 901
      IF (ISTOPX.GT.NPTSX.OR.ISTOPY.GT.NPTSY) GO TO 901
      IF (ISTRTL.LT.1.OR.ISTRTL.GT.ISTOPL)    RETURN
C
C          THE ROUTE-TRACING FLAG, THE CURVE METHOD AND THE
C          CHARACTER MAGNIFICATION ARE SAVED, THEN TRACE IS
C          DISABLED, CURVE METHOD (1) SET AND MAGN. (6) SET.
C
      IPRSAV= IPRINT
      IPRINT= 0
      CURVED= .FALSE.
      ANGSAV= STANG0
      UNSAVE= ANGCON
      ANGCON= 1.0
      IF (LLBCON.GT.0) CALL CTRORI(0.5)
C
      CALL G0AUTO(XGRIDS,YGRIDS,ISTRTX,ISTOPX,ISTRTY,ISTOPY,1)
      XHERE= XPLOT0
      YHERE= YPLOT0
      CWIDX= (X2WND0-X1WND0)/(ISTOPX-ISTRTX)
      CWIDY= (Y2WND0-Y1WND0)/(ISTOPY-ISTRTY)
C
C          TO BEGIN WITH, SOME CONSTANTS ARE SET UP.
C
      ISTRX1= ISTRTX+1
      ISTRY1= ISTRTY+1
      ISTPX1= ISTOPX-1
      ISTPY1= ISTOPY-1
C
C          LOOP-100 DOES ALL THE CONTOURS AT EACH HEIGHT IN TURN.
C          THE BIT ARRAY <MAPBIT> IS INITIALISED FOR EACH HEIGHT,
C          ELEMENTS BEING SET TO '1' WHEN A CONTOUR LINE OF THE
C          CURRENT HEIGHT CROSSES AN EDGE IN AN UPWARDS DIRECTION
C          PROVIDED 'HIGH GROUND' LIES TO THE RIGHT OF THE LINE.
C

c...dev
      DO 100 LEVELH= ISTRTL,ISTOPL
        HEIGHT= CLEVLS(LEVELH)

c...dev
        WRITE(6,*) 'HEIGHT=',height,ISTOPL,ISTRTL,iconcol

        LEVEL= LEVELH
C
        DO 200 ISTPTY= ISTRY1,ISTPY1
          DO 200 ISTPTX= ISTRX1,ISTOPX
            IBIT= 0
            IF (SURFAS(ISTPTX-1,ISTPTY).LT.HEIGHT.AND.
     &          SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)     IBIT= 1
C
            IXWRD= (ISTPTX-ISTRX1)/NBITSW+1
            IXBIT= MOD((ISTPTX-ISTRX1),NBITSW)+1
            CALL G4PUTB(MAPBIT(IXWRD,(ISTPTY-ISTRTY)),IXBIT,IBIT)
  200   CONTINUE
C
C          <OPENCO> IS INITIALISED AND A SEARCH FOR OPEN
C          CONTOURS BEGINNING ON EACH EDGE IS MADE IN TURN:
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTY>= <ISTRTY>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        OPENCO= .TRUE.
        ISTPTY=  ISTRTY
C
        DO 300 ISTPTX= ISTRX1,ISTOPX

       icon = 0

          IF (SURFAS(ISTPTX-1,ISTPTY).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,-1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'A'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  300   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTX>= <ISTOPX>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTX= ISTOPX
C
        DO 400 ISTPTY= ISTRY1,ISTOPY

       icon = 0

          IF (SURFAS(ISTPTX,ISTPTY-1).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,0,-1)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'B'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  400   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTY>= <ISTOPY>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTY= ISTOPY
C
        DO 500 NEGIX= ISTRTX,ISTPX1
          ISTPTX= ISTPX1+ISTRTX-NEGIX

       icon = 0

          IF (SURFAS(ISTPTX+1,ISTPTY).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour            
          ENDIF
          WRITE(6,*) 'C'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  500   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTX>= <ISTRTX>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTX= ISTRTX
C
        DO 600 NEGIY= ISTRTY,ISTPY1
          ISTPTY= ISTPY1+ISTRTY-NEGIY

       icon = 0

          IF (SURFAS(ISTPTX,ISTPTY+1).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,0,1)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'D'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  600   CONTINUE
C
C          A SEARCH IS MADE FOR CLOSED CONTOURS. <OPENCO> IS
C          SET TO .FALSE. AND <MAPBIT> IS SCANNED: IF THE
C          ELEMENT (ISTPTX,ISTPTY) IS '1', A CLOSED CONTOUR
C          STARTS FROM THE LINE JOINING (ISTPTX-1,ISTPTY) TO
C          (ISTPTX,ISTPTY) AND IT IS FOLLOWED USING <G0CTIA>.
C
        OPENCO= .FALSE.
C
c slmod begin
        IF (.FALSE..AND.icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
          ENDIF
          WRITE(6,*) 'D.5'
          CALL PlotContour
          icon = 0
        ENDIF
c slmod begin
        DO 700 NEGIY= ISTRY1,ISTPY1
          ISTPTY= ISTOPY+ISTRTY-NEGIY
C
          DO 700 NEGIX= ISTRTX,ISTPX1
            ISTPTX= ISTOPX+ISTRTX-NEGIX
            IXWRD= (ISTPTX-ISTRX1)/NBITSW+1
            IXBIT= MOD((ISTPTX-ISTRX1),NBITSW)+1
            CALL G4GETB(MAPBIT(IXWRD,(ISTPTY-ISTRTY)),IXBIT,IBIT)

       icon = 0

            IF (IBIT.EQ.1)
     &              CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,
     &                          ISTOPY,NPTSY,XGRIDS,YGRIDS,-1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048)THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'E'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  700   CONTINUE
  100 CONTINUE
C
C          THE ENTRY STATE IS ALWAYS RESTORED BEFORE ENDING.
C
      CALL POSITN(XHERE,YHERE)
      CALL CTRORI(ANGSAV)
      ANGCON= UNSAVE
      IPRINT= IPRSAV
      RETURN
C
 901  NUMERR= 11
      IF (ERRON) CALL G0ERMS
C
      RETURN
      END
c
c ======================================================================
c
c subroutine: PlotContour
c
c
      SUBROUTINE PlotContour
      use mod_params
      IMPLICIT none
c
c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'comgra'
c      INCLUDE 'colours'
c      INCLUDE 'slcom'
c      INCLUDE 'slout'     
c
      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

      REAL       TOL
      PARAMETER (TOL=1.0E-5)

      INTEGER i1,i2,nver(2),imin,icut,ngrp,imid
      LOGICAL farenough
      REAL    xtmp(2048),ytmp(2048),xver(2048,2),yver(2048,2),xmin,xmid,
     .        xmax

c...  Check if end points of the contour make sense:
      IF (SQRT((xcon(1)-xcon(icon))**2+
     .         (ycon(1)-ycon(icon))**2).GT.0.01) THEN

        IF     (ABS(xcon(1)-xcon(icon)).LT.TOL) THEN
        ELSEIF (ABS(ycon(1)-ycon(icon)).LT.TOL) THEN
        ELSEIF (xcon(1).GT.xcon(icon).AND.
     .          ycon(1).LT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(1)
          ycon(icon) = ycon(icon-1)
        ELSEIF (xcon(1).LT.xcon(icon).AND.
     .          ycon(1).LT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(icon-1)
          ycon(icon) = ycon(1)
c        ELSEIF (xcon(1).LT.xcon(icon).AND.
c     .          ycon(1).LT.ycon(icon)) THEN
c          icon = icon + 1
c          xcon(icon) = xcon(icon-1)
c          ycon(icon) = ycon(1)
        ELSEIF (xcon(1).LT.xcon(icon).AND.
     .          ycon(1).GT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(icon-1)
          ycon(icon) = ycon(1)
        ELSE
c          WRITE(0,*) 'xCON=',
c     .      xcon(1),ycon(1),xcon(icon),ycon(icon)
          CALL ER('PlotContour','Unknown situation',*99)
        ENDIF

      ENDIF
 

c      IF (.FALSE..AND.icon.GT.999) THEN
c
      IF (icon.GT.2000) THEN

        WRITE(0,*) 'ICON= ',icon
        STOP 'CONTOUR ARRAY TOO LARGE, INCREASE ARRAY SIZES IN GHOST'
c...    Split array of verticies

        ngrp = 2

c...    Find the mid-point of the contour:
        xmin =  HI
        xmax = -HI
        DO i1 = 1, icon
          IF (xcon(i1).LT.xmin) THEN
            imin = i1
            xmin = xcon(i1)
          ENDIF
          xmax = MAX(xmax,xcon(i1))
        ENDDO
c...    Copy contour array starting at imin:
        i1 = imin
        DO i2 = 1, icon
          xtmp(i2) = xcon(i1)
          ytmp(i2) = ycon(i1)
          i1 = i1 + 1
          IF (i1.GT.icon) i1 = 1
        ENDDO

c...    Find the x-axis midpoint and cut the array there:
        xmid = 0.5 * (xmin + xmax)
        DO i1 = 1,icon
          IF (xtmp(i1).GT.xmid) THEN
            imid = i1
            EXIT
          ENDIF
        ENDDO
        farenough = .FALSE.
        icut = 0
        DO i1 = imid+1,icon
          IF (xtmp(i1).GT.1.1*xmid) farenough = .TRUE.
          IF (farenough.AND.xtmp(i1).LT.xmid) THEN
            icut = i1
            EXIT
          ENDIF
        ENDDO
        IF (imid.EQ.icon) STOP 'SHIT'

        ngrp = 2

        nver(1) = 0
        DO i1 = 1, imid
          nver(1) = nver(1) + 1
          xver(nver(1),1) = xtmp(i1)
          yver(nver(1),1) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      1,nver(1),xver(nver(1),1),yver(nver(1),1)
        ENDDO
        DO i1 = icut, icon
          nver(1) = nver(1) + 1
          xver(nver(1),1) = xtmp(i1)
          yver(nver(1),1) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      1,nver(1),xver(nver(1),1),yver(nver(1),1)
        ENDDO

        nver(2) = 0
        DO i1 = imid, icut
          nver(2) = nver(2) + 1
          xver(nver(2),2) = xtmp(i1)
          yver(nver(2),2) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      2,nver(2),xver(nver(2),2),yver(nver(2),2)
        ENDDO

c        WRITE(0,*) '--->',nver(1),nver(2),icon

      ELSE
c...    Copy array of verticies as-is:
        ngrp = 1
        nver(1) = icon
        DO i1 = 1, icon
          WRITE(6,*) 'CNTR:',i1,xcon(i1),ycon(i1)
          xver(i1,1) = xcon(i1)
          yver(i1,1) = ycon(i1)
        ENDDO
      ENDIF

      DO i1 = 1, ngrp
        IF (nver(i1).GT.2048) WRITE(0,*) 'WARNING: NVER TOO BIG'
        CALL PTPLOT(xver(1,i1),yver(1,i1),1,nver(i1),1)
      ENDDO

      RETURN
99    STOP
      END
