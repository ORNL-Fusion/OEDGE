c     -*-Fortran-*-
      INTEGER FUNCTION CH1(STRING)
      IMPLICIT none
      CHARACTER STRING*(*)
      INTEGER   I1
      CH1=1
      DO I1 = 1, LEN_TRIM(STRING)
        IF (STRING(I1:I1).NE.' ') EXIT
        CH1=CH1+1
      ENDDO
      RETURN
      END
c
c ======================================================================
c
      SUBROUTINE GRTSET_TRIM (TITLE,REF,VIEW,PLANE,JOB,XMIN,XMAX,
     >    YMIN,YMAX,TABLE,XLABEL,YLABEL,IFLAG,SMOOTH,IDRAW,ANLY,NBBS)
      implicit none
      REAL      XMIN,XMAX,YMIN,YMAX
      INTEGER   IFLAG,IDRAW,NBBS
c slmod begin - new
c...  Allow for longer y-axis labels:
      CHARACTER YLABEL*(*),XLABEL*36,TABLE*36
c
c      CHARACTER YLABEL*36,XLABEL*36,TABLE*36

      INTEGER CH1

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

c slmod end
      CHARACTER TITLE*(*),JOB*(*),SMOOTH*(*)
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
      include 'params'
      include 'slout'
      include 'colours'
      integer init_col,get_col,next_col
      external init_col,get_col,next_col
      include 'comgra'
c
c     Declare variables
c
      real power,tmin,tmax,rdata(10)
      integer iexp,l,lenstr,iten,l1,idata(10)
C
c      CALL THICK(2)

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
c      CALL CTRMAG (20)
c
c     added for AUG
c
c      CALL BACCOL(5)
c
      CALL LINCOL (defcol)
c      CALL THICK  (2)
c slmod begin - new
c      IF (iopt_ghost.EQ.2) iflag = 3

      CALL ENQCHR(idata,rdata)
      CALL CTRMAG(20)  

      L = LENSTR(TITLE)
      IF (L.GT.76) THEN
        DO L1 = 76, 1, -1
          IF (TITLE(L1:L1).EQ.' ') EXIT
        ENDDO
c...    Maybe there are no spaces, so set L1 arbitrarily:
        IF (L1.EQ.1) L1 = 76
        CALL PCSCEN (0.68, 0.98, TITLE(1:L1-1))   
        CALL PCSCEN (0.68, 0.95, TITLE(L1+1:L))           
      ELSE
        CALL PCSCEN (0.68, 0.98, TITLE(:L))
c        CALL PCSCEN (0.68, 0.97, TITLE(1:L1-1))   
c        CALL PCSCEN (0.68, 0.94, TITLE(L1+1:L))           
c      ELSE
c        CALL PCSCEN (0.68, 0.95, TITLE(:L))
      ENDIF

      CALL CTRMAG(idata(2))

c      L = MIN(LENSTR(TITLE),75)
c      CALL PCSCEN (0.68, 0.95, TITLE(:L))

c
c      L = LENSTR(TITLE)
c      CALL PCSCEN (0.8, 0.95, TITLE(:L))
c slmod end
      CALL THICK  (1)
c      CALL CTRMAG (12)

c      CALL CTRMAG (12)
      CALL PLOTCS (0.01, 0.01, JOB(CH1(job):LEN_TRIM(job)))

      IF (.FALSE.) THEN
        L = LENSTR (REF)
        CALL PCSCEN (1.16, 0.12, REF(:L))
c        CALL CTRMAG (12)
        L=LENSTR(VIEW)
        CALL PCSCEN (1.16, 0.27, VIEW(:L))
        L=LENSTR(PLANE)
        CALL PCSCEN (1.16, 0.31, PLANE(:L))
        L=LENSTR(ANLY)
        CALL PCSCEN (1.16, 0.35, ANLY(:L))
        CALL PCSCEN (1.16, 0.39, SMOOTH(37:72))
        CALL PCSCEN (1.16, 0.43, SMOOTH( 1:36))
c        CALL CTRMAG (ISPOT)
      ENDIF
C
C---- DRAW FRAMES
C
      CALL LINCOL (defcol)

c slmod begin
      IF (ABS(xmin).LT.0.01.OR.ABS(xmax).GT.1.0E+3) THEN
        ITEN = IEXP(XMIN, XMAX)
        POWER = 10.0**(-ITEN)
        TMIN = XMIN * POWER
        TMAX = XMAX * POWER
      ELSE
        ITEN = 0
        TMIN = XMIN
        TMAX = XMAX
      ENDIF
c
c        ITEN = IEXP(XMIN, XMAX)
c        POWER = 10.0**(-ITEN)
c        TMIN = XMIN * POWER
c        TMAX = XMAX * POWER
c slmod end

      CALL PSPACE (map1x,map2x,map1y,map2y)
      CALL MAP    (TMIN , TMAX,0.1  ,0.9  )

      IF (iopt_ghost.EQ.0) CALL CTRMAG (10)

c...  Use IFLAG here, not IOP_GHOST2:
      
      IF (xlabel(1:LEN_TRIM(xlabel)).EQ.'none') GOTO 50

      IF (OPT_XSCALE.EQ.2) THEN
        CALL XLOGSCALE
      ELSE
c        CALL CTRMAG (12)
        CALL XSCALE
      ENDIF
      
      CALL PSPACE (map1x,map2x,map1y,map2y)
c      CALL PSPACE (map1x,map2x,map1y-0.07,map2y)
      CALL MAP    (0.0 ,1.0 ,0.0 , 1.0)
      
      IF (ITEN .NE. 0) THEN
        CALL POSITN (0.8,0.04)
c        CALL CTRMAG (14)
        CALL TYPECS ('X10')
c        CALL CTRMAG (12)
        CALL TYPENI ((ITEN))
      ENDIF
      
      L = LENSTR (XLABEL)
      IF (iopt_ghost.EQ.0) THEN
        IF (IFLAG.NE.3) CALL PCSCEN (0.5,0.025,XLABEL(:L))
      ELSE
        CALL PSPACE (0.0, 1.35, 0.0, 1.0)
        CALL MAP    (0.0, 1.35, 0.0, 1.0)
        IF (IFLAG.NE.3) CALL PCSCEN (0.5*(MAP1X+MAP2X),MAP1Y-0.042,  ! -0.05
     .                               XLABEL(:L))
      ENDIF
c      IF (IFLAG.NE.3) CALL PCSCEN (0.5,0.035,XLABEL(:L))
C
C---- DRAW Y AXIS AND LABELS
C
50    CONTINUE

      IF (YMAX .EQ. YMIN) THEN
         WRITE(6, *)  ' ERROR - ROUTINE GRTSET '
         WRITE(6, *)  ' SCALE ',YLABEL,' FMIN = FMAX ',YMIN
      ENDIF
      ITEN = IEXP(YMIN, YMAX)
      POWER = 10.0**(-ITEN)
      TMIN = YMIN * POWER
      TMAX = YMAX * POWER

      CALL PSPACE (map1x,map2x,map1y,map2y)
      CALL MAP    (0.0  ,1.0  ,TMIN ,TMAX )

      CALL LINCOL (defcol)
      IF (iopt_ghost.EQ.0) CALL CTRMAG (10)

      IF (ylabel(1:LEN_TRIM(ylabel)).NE.'none') THEN
        CALL YSCALE
      ELSE
        CALL POSITN(0.0,TMIN)
        CALL JOIN(0.0,TMAX)
      ENDIF

c...  Clear out the ticks inside the box:
c      CALL PSPACE(map1x,map1y,map1y,map2y)
c      CALL MAP(0.0,1.0,0.0,1.0)
c      CALL BOX(0.0,1.0,0.0,1.0)

      

      CALL PSPACE (0.0,1.35,map1y,map2y)
      CALL MAP    (0.0,1.35,0.0  ,1.0  )


      IF (ylabel(1:LEN_TRIM(ylabel)).NE.'none') THEN
        CALL CTRORI (90.0)
        IF (ITEN .NE. 0) THEN
           CALL POSITN (0.02, 0.03)
c           CALL CTRMAG (14)
           CALL TYPECS ('X10')
c           CALL CTRMAG (12)
           CALL TYPENI ((ITEN))
        ENDIF
        CALL LINCOL (defcol)
        IF (iopt_ghost.EQ.0) CALL CTRMAG (14)
        CALL THICK  (2)
        L = LENSTR (YLABEL)
c        CALL CTRMAG(11)

        CALL PSPACE (0.0,1.35,0.0,1.0)
        CALL MAP    (0.0,1.35,0.0,1.0)
        CALL PLOTST (map1x-0.07,map1y+0.1*(map2y-map1y),YLABEL(:L))

c        CALL PCSEND (0.02,0.99,YLABEL(:L))
        CALL CTRORI (0.0)
        CALL THICK  (1)
      ENDIF

      RETURN
      END
c
c ======================================================================
c


      SUBROUTINE Plot981(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
c      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure

      integer cngs,cntropt

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256,graph6*128

      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,XPOINT
      CHARACTER*72 YLAB
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,PLABS(-2:MAXPLRP),KLAB

      CHARACTER*128 elabs(MAXNGS)

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact

      REAL tauii,pii
      INTEGER nenum,tenum,opt_const,plot_mode(30),iter1,iter2,
     .        xaxis,ring,mode,inorm(MAXNGS)

      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
     .        max1,max2,ynorm(MAXNGS)

      REAL cdata(MAXNKS,MAXNRS)
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,SUM(10),VAL,VALTOT
      CHARACTER*72 SMOOTH
      integer icntr
      integer nconts,nclev
      real conts(maxpts),clev(maxpts)

      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)


      CHARACTER*128 dataline,cdum1,cdum2

c...630:
      INTEGER i,k
      REAL    PLTMAX,PLTMIN
      REAL LOUTS(MAXSEG),LWIDS(MAXSEG),LVALS(MAXSEG,MAXNGS)
      REAL ydata(MAXSEG)
      integer llabs(maxseg)

c...980:
      INTEGER NUMTHE,AVPTS,ATYPE
      INTEGER numth2,numthe1(MAXNGS)
      INTEGER IGNORS(MAXNGS),ITEC,NAVS
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS),
     .     touts1(MAXTHE,MAXNGS)
      REAL TOUTS2(MAXTHE),TWIDS2(MAXTHE),TVALS2(MAXTHE,MAXNGS),
     .     dum1(MAXTHE),dum2(MAXTHE),dum3(MAXTHE),dum4(MAXTHE),
     .     dum5(MAXTHE),dum6(MAXTHE),dum7(MAXTHE),dum8(MAXTHE),
     .     den1,teb1
      REAL ZOBS,ROBS,DRAD,DTHE,THEMIN,THEMAX,themin_start
      real theres,mfact
      REAL WLNGTH
      real    zsuma(maxizs),cvalsa(maxnks,maxnrs)
      REAL LEVEL,AVS(0:100),VMIN,VMAX
      CHARACTER ADASID*80,PECTITLE*120,PLABAD*36
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      CHARACTER ADASGR*8,ADASTY*80
C
C     Second sets of ADAS data for RATIO plots
C
      character graph5*80,adasid2*80,adasex2*3
      integer   adasyr2,isele2,iselr2,iselx2,iseld2
      integer   iz_state,z_atom,iz_state2,z_atom2
      INTEGER plot
      character graph2*80,graph3*80,graph4*80
      INTEGER ISELE,ISELR,ISELX,iseld,iseldef
      INTEGER IADAS,NPAIRS,IRCODE,IKK,LEN,LENSTR
      INTEGER IK,II,IT,LT,UT,IREF,IYMIN,IYMAX,IR,JD,cnt
      INTEGER IZMIN,IZMAX,IW,LW,UW
      REAL PLRPAD(MAXNKS,MAXNRS)
      LOGICAL plotr
      real zadj
      REAL peak1,peak2,array(MAXNKS,MAXNRS)
      INTEGER plotcode,line,iopt1,iopt2,iopt3
      CHARACTER cname*256,resdir*256,cmnd*256
      LOGICAL status,oldraw


      COMMON /NSMOOTH/ NUMSMOOTH,cgrprint
      INTEGER NUMSMOOTH,cgrprint

c...  For reading from the .experiment file (UNIT=13):
      INTEGER    MAXTDAT     ,MAXCOLS   
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)
      INTEGER   etype,ndata,ncol
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS)
      CHARACTER datatitle*128

c...  981:
      INTEGER count,nlines,cline(100),idum1,multiplot,size
      REAL    yycen,yydis,lines(100),deltax,deltay,xpos,ypos
      CHARACTER dummy*5000,caption*5000

      DATA multiplot /0/

      nview = ' '
      plane = ' '
      anly  = ' '
      TABLE = 'SYMBOL TABLE'

c...  Use GRTSET_TRIM:
      slopt4 = 1

c...  Trigger multi-plot code:
      IF (iopt.EQ.99) THEN
        iopt = 1
        multiplot = 1
      ELSEIF (multiplot.GT.0) THEN
        multiplot = multiplot + 1
c...    Spagetti (sorry, just want to get this going):
        GOTO 40
      ELSE
        multiplot = 0
      ENDIF

      CALL THICK2(4)

      yycen = 0.5 * (yymin + yymax)
      yydis = yymax - yymin

      deltax = 0.80
      deltay = 0.18
c      deltay = 0.20

      yymin = yycen - 0.5 * (deltay / deltax) * (xxmax - xxmin)
      yymax = yycen + 0.5 * (deltay / deltax) * (xxmax - xxmin)
c      yymin = yycen - 0.25 * yydis
c      yymax = yycen + 0.25 * yydis

      CALL SLSET (0.10,0.10+deltax,0.80,0.80+deltay,
     .            xxmin,xxmax,yymin,yymax)
c      CALL SLSET (0.10,0.10+deltax,0.70+0.008,0.70+deltay,
c     .            xxmin,xxmax,yymin,yymax)
c      CALL SLSET (0.10,0.10+deltax,0.05,0.05+deltay,
c     .            xxmin,xxmax,yymin,yymax)
c      CALL SUPIMP('FULL')
      CALL SUPIMP('PARTIAL')

      CALL FULL
      CALL THICK(1)

      CALL MAP    (0.0,1.0,0.0,1.0)

      CALL POSITN (0.0,0.0)
      CALL JOIN   (0.0,1.0)              
      CALL POSITN (0.0,1.0)
      CALL JOIN   (1.0,1.0)              
      CALL POSITN (1.0,1.0)
      CALL JOIN   (1.0,0.0)              
      CALL POSITN (1.0,0.0)
      CALL JOIN   (0.0,0.0)              


      count  = 0
      nlines = 0

40    READ(5,'(A256)',END=50) graph1
      BACKSPACE 5
      IF (graph1(8:11).EQ.'View'.OR.graph1(8:11).EQ.'VIEW'.OR.
     .    graph1(8:11).EQ.'view') THEN

c...    Store views:
        READ(graph1,*) cdum1,idum1,rdum1,rdum1,
     .                 idum1,(lines(i1+nlines),i1=1,idum1)        
        DO i1 = 1, idum1
          cline(i1+nlines) = count + 1
        ENDDO
        nlines = nlines + idum1

        count = count + 1
c        CALL LinCol(ncols+count)
c        CALL FilCol(ncols+count)

        IF (multiplot.LE.1) CALL DrawHalpha(1,1)
        GOTO 40
      ELSE
        CALL LinCol(1)
        CALL FilCol(1)
      ENDIF
50    CONTINUE


c       DO i1 = 1, 1


c...True space
c         CALL SLSET (0.10,0.90,0.50,0.90,0.0,0.0,0.0,0.0)

c...dev
         CALL CTRMAG(12)

         IF (multiplot.GT.0) THEN
           i1 = multiplot
           CALL SLSET (0.10,0.90,0.80-(i1  )*deltay,
     .                           0.80-(i1-1)*deltay,0,0.0,0.0,0.0)
c           CALL SLSET (0.10,0.90,0.70-(i1  )*0.20,
c     .                           0.70-(i1-1)*0.20,0.0,0.0,0.0,0.0)
c           CALL SLSET (0.10,0.90,0.70-(i1  )*0.20+0.008,
c     .                           0.70-(i1-1)*0.20,0.0,0.0,0.0,0.0)
           IF (multiplot.EQ.3) THEN
             iopt_ghost = 3
           ELSE
             iopt_ghost = 2   
           ENDIF
         ELSE
           CALL SLSET (0.10,0.90,0.25,0.65,0.0,0.0,0.0,0.0)
c           CALL SLSET (0.10,0.90,0.15,0.65,0.0,0.0,0.0,0.0)
         ENDIF

c...     Blank ref for now:
         REF = '                                   '

         CALL Plot980(job,graph,ref,title,iopt,
     .               xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .               ismoth,ignors,itec,avs,navs)

c...     Add picture to legend:
c        slopt = 2
c        CALL SLSET(0.95,1.31,0.38,0.73,xxmin,xxmax,yymin,yymax)
c        CALL DrawHalpha(ret,0)
c        CALL SUPIMP ('PARTIAL')


c...     Draw the views on the LOS plot:
         CALL PSPACE (map1x,map2x,map1y,map2y)
         CALL MAP    (cxmin,cxmax,0.0,1.0)
         CALL BROKEN(6,6,6,6)
         CALL LinCol(1)
         DO i2 = 1, nlines
c           CALL LinCol(ncols+cline(i2))
c           CALL FilCol(ncols+cline(i2))
           CALL POSITN (lines(i2),0.0)
           CALL JOIN   (lines(i2),1.0)        
         ENDDO


c...    Add a comment to the plot:
49      READ(5,'(A5000)') dummy
        IF   (dummy(8:11).EQ.'Cmnt'.OR.dummy(8:11).EQ.'cmnt'.OR.
     .        dummy(8:11).EQ.'CMNT') THEN
      
          READ (dummy,*) cdum1,xpos,ypos,size,caption
      
c...      Annotate graph:
          CALL PSPACE (map1x,map2x,map1y,map2y)
          CALL MAP    (0.0,1.0,0.0,1.0)
          CALL LinCol(1)
          CALL CTRMAG(size)
          CALL PLOTST(xpos,ypos,caption(1:LEN_TRIM(caption)))
c...      Another comment:        
          GOTO 49
        ELSE
          BACKSPACE 5
        ENDIF




c...    Complete the plot frame that isn't quite done right
c       in GRTET_TRIM:
        CALL PSPACE (map1x,map2x,map1y,map2y)
        CALL MAP    (0.0,1.0,0.0,1.0)

        CALL FULL
        CALL THICK(1)

        CALL POSITN (0.0,1.0)
        CALL JOIN   (1.0,1.0)              
        CALL POSITN (1.0,1.0)
        CALL JOIN   (1.0,0.0)              

c      ENDDO





c...  Add a caption to the plot:
      READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Note'.OR.dummy(8:11).EQ.'note'.OR.
     .      dummy(8:11).EQ.'NOTE') THEN
        READ(dummy,*) cdum1,xpos,ypos,size,caption
        CALL AddCaption(caption,xpos,ypos,size)
      ELSE
        BACKSPACE 5
      ENDIF


c      IF (multiplot.EQ.0.OR.multiplot.EQ.3) THEN
c        CALL FRAME
c        multiplot = 0
c      ENDIF

c...  Add a caption to the plot:
      READ(5,'(A5000)') dummy
      IF   (dummy(8:14).EQ.'Noframe'.OR.dummy(8:12).EQ.'noframe'.OR.
     .      dummy(8:14).EQ.'NOFRAME') THEN
      ELSE
        BACKSPACE 5
        CALL FRAME
        multiplot = 0
      ENDIF


      RETURN
99    STOP
      END
