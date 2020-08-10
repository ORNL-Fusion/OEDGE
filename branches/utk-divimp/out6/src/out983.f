c
c ======================================================================
c
c
   

c
c ======================================================================
c
c subroutine: Plot982
c
c 2D neutral pressure plot
c
c
      SUBROUTINE Plot983(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_comgra
      use mod_colours
      use mod_printopt
      use mod_slcom
      use mod_slout
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'comgra'
c     INCLUDE 'colours'
c     include 'printopt'

c     INCLUDE 'slcom'
c     INCLUDE 'slout'
c      INCLUDE 'outcom'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure

      integer cngs,cntropt

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*(*),TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256,
     .          graph6*128

      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,XPOINT
      CHARACTER*72 YLAB
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,PLABS(-2:MAXPLRP),KLAB

      CHARACTER*128 elabs(MAXNGS)

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,i3,ret,plt_opt1,osm_cfs,npts,ntime,
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
      REAL    nemin,nestep,temin,temax,neTe,xrange1,xrange2,
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


c...  983:
      INTEGER calc_random_seed

      INTEGER MAXSOLID
      PARAMETER (MAXSOLID=10000)

      INTEGER axis,fp1,index,axis1(3),idum1,idum2,idum3,idum4,
     .        trackno,seed,colind,lasttrack,nseg,nsur,ndiv1,ndiv2
      LOGICAL skipnext,skiptrack,axisdata,solid
      REAL    x1,y1,z1,x2,y2,z2,angle1(3),frac,ran,dangle,ang,phi,
     .        xcen(MAXSOLID),ycen(MAXSOLID),zcen(MAXSOLID),theta,
     .        xsur(4,MAXSOLID),ysur(4,MAXSOLID),zsur(4,MAXSOLID),
     .        dsur(MAXSOLID),
     .        xtemp,ytemp,ztemp,cxtmp,cytmp,cztmp,dtemp,
     .        xsur2(4),ysur2(4),zsur2(4),
     .        v1x,v2x,v3x,v4x,v1y,v2y,v3y,v4y,v1z,v2z,v3z,v4z,
     .        frac1,frac2,lightv(MAXSOLID),lightvtmp
      REAL*8 p1(3,6),p2(3,6),xshift,yshift,zshift
      REAL*8 mat(3,3),angle,posx,posy,posz,
     .       light(3),normv(3),a1,a2,a3,b1,b2,b3


    


      XLAB   = '   R  (M)'
      YLAB   = '   Z  (M)'

c      title  = title(1:77)//'   '
      NVIEW  = '                                    '
      PLANE  = '                                    '
      SMOOTH = '                                    '//
     .         '                                    '
      ANLY   = '                                    '
      REF    = graph(14:LEN_TRIM(graph))

      CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     .             YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)


c...  Draw polygons:
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)


      axis1 (1) = 1
      axis1 (2) = 2
      axis1 (3) = 3
      angle1(1) = 0.0
      angle1(2) = 0.0
      angle1(3) = 0.0

      axisdata = .FALSE.

      READ(5,'(A256)') graph1
      IF (graph1(8:11).EQ.'Axis'.OR.graph1(8:11).EQ.'AXIS'.OR.
     .    graph1(8:11).EQ.'axis') THEN
        READ(graph1,*) cdum1,(axis1(i1),angle1(i1),i1=1,3)
        axisdata = .TRUE.
      ELSEIF (iopt.GT.0) THEN
        CALL ER('983','Rotation data not found',*99)
      ELSE
        BACKSPACE 5
      ENDIF

c...  Setup transformation matrix:
      CALL Calc_Transform(mat,0.0D0,1,0)
      DO i1 = 1, 3
        angle = DBLE( angle1(i1)*PI/180.0)
        CALL Calc_Transform(mat,angle,axis1(i1),1)
      ENDDO


c      axis  = 1
c      angle = DBLE(-45.0*PI/180.0)
c      CALL Calc_Transform(mat,angle,axis,1)
c      axis  = 2
c      angle = DBLE( 45.0*PI/180.0)
c      CALL Calc_Transform(mat,angle,axis,0)



c...  Move origin:
      IF (.FALSE..AND.(iopt.GT.0.OR.axisdata)) THEN
c        xshift = DBLE(-rxp)
c        yshift = DBLE(-zp(idds(irsep,1)))
c        zshift =  0.0D0
        xshift = -0.650D0
        yshift = +0.450D0
        zshift = -0.235D0
      ELSE
        xshift = 0.0D0
        yshift = 0.0D0
        zshift = 0.0D0
      ENDIF

      CALL LINCOL(ncols+3)

      DO ir = 2, nrs
        IF (ir.EQ.irtrap) CYCLE
c        IF (ir.EQ.irwall.OR.ir.EQ.irtrap) CYCLE

        DO ik = 1, nks(ir)

          id = korpg(ik,ir)
         
          DO in1 = 1, nvertp(id)
            in2 = in1 + 1
            IF (in2.EQ.nvertp(id)+1) in2 = 1

            IF (.NOT.(ir.EQ.irwall  .AND.in1.EQ.2.OR.
c            IF (.NOT.(ir.EQ.irwall-1.AND.in1.EQ.2.OR.
     .                ir.EQ.irtrap+1.AND.in1.EQ.4.OR.
     .                ir.EQ.irsep   .AND.in1.EQ.4.OR.
     .                ir.EQ.2       .AND.in1.EQ.2)) CYCLE

            IF (ir.EQ.irwall) id = korpg(ikins(ik,ir),irins(ik,ir))


            p1(1,1) = DBLE(rvertp(in1,id)) + xshift
            p1(2,1) = DBLE(zvertp(in1,id)) + yshift
            p1(3,1) = 0.0D0                + zshift
            p2(1,1) = DBLE(rvertp(in2,id)) + xshift
            p2(2,1) = DBLE(zvertp(in2,id)) + yshift
            p2(3,1) = 0.0D0                + zshift

            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))

            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 

          ENDDO
        ENDDO
      ENDDO

      CALL LINCOL(55)

      DO i1 = 1, nvesm+nvesp

        p1(1,1) = DBLE(rvesm(i1,1)) + xshift
        p1(2,1) = DBLE(zvesm(i1,1)) + yshift
        p1(3,1) = 0.0D0             + zshift
        p2(1,1) = DBLE(rvesm(i1,2)) + xshift
        p2(2,1) = DBLE(zvesm(i1,2)) + yshift
        p2(3,1) = 0.0D0             + zshift
	
        call transform_vect(mat,p1(1,1))
        call transform_vect(mat,p2(1,1))
	
        CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
        CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 

      ENDDO


c...  Draw toroidally extended surfaces:
      solid = .TRUE.
      nsur = 0

      IF (eirntorseg.NE.0) THEN

        dangle = 360.0 / REAL(eirntorseg) / RADDEG

        DO i1 = 1, nvesm+nvesp

c...      Don't draw anything above the x-point:
          IF (zvesm(i1,1).GT.zxp.OR.zvesm(i1,2).GT.zxp) CYCLE

c...      Don't draw other things:
          IF (jvesm(i1).NE.7.0.AND.jvesm(i1).NE.8.0.AND.
     .        jvesm(i1).NE.0.0) CYCLE

          DO ang = 0.0, 2.0*PI-dangle, dangle

            p1(1,1) = DBLE(rvesm(i1,1))
            p1(2,1) = DBLE(zvesm(i1,1))
            p1(3,1) = DBLE(rvesm(i1,1) * TAN(-0.5 * dangle))
            p2(1,1) = DBLE(rvesm(i1,2))
            p2(2,1) = DBLE(zvesm(i1,2))
            p2(3,1) = DBLE(rvesm(i1,2) * TAN(-0.5 * dangle))
	
            p1(1,2) = DBLE(rvesm(i1,2))
            p1(2,2) = DBLE(zvesm(i1,2))
            p1(3,2) = DBLE(rvesm(i1,2) * TAN(-0.5 * dangle))
            p2(1,2) = DBLE(rvesm(i1,2))
            p2(2,2) = DBLE(zvesm(i1,2))
            p2(3,2) = DBLE(rvesm(i1,2) * TAN(+0.5 * dangle))

            p1(1,3) = DBLE(rvesm(i1,2))
            p1(2,3) = DBLE(zvesm(i1,2))
            p1(3,3) = DBLE(rvesm(i1,2) * TAN(+0.5 * dangle))
            p2(1,3) = DBLE(rvesm(i1,1))
            p2(2,3) = DBLE(zvesm(i1,1))
            p2(3,3) = DBLE(rvesm(i1,1) * TAN(+0.5 * dangle))
	  
	    p1(1,4) = DBLE(rvesm(i1,1))
            p1(2,4) = DBLE(zvesm(i1,1))
            p1(3,4) = DBLE(rvesm(i1,1) * TAN(+0.5 * dangle))
            p2(1,4) = DBLE(rvesm(i1,1))
            p2(2,4) = DBLE(zvesm(i1,1))
            p2(3,4) = DBLE(rvesm(i1,1) * TAN(-0.5 * dangle))

c            IF (solid) THEN
            IF (solid.AND.ang.LT.36.0*PI/180.0) THEN
              nsur = nsur + 1
            ENDIF
	  
            DO i2 = 1, 4
              x1 = p1(1,i2)
              z1 = p1(3,i2)
              x2 = p2(1,i2)
              z2 = p2(3,i2)
	    
              p1(1,i2) = DBLE( COS(ang) * x1 - SIN(ang) * z1)
              p1(3,i2) = DBLE(+SIN(ang) * x1 + COS(ang) * z1)
              p2(1,i2) = DBLE( COS(ang) * x2 - SIN(ang) * z2)
              p2(3,i2) = DBLE(+SIN(ang) * x2 + COS(ang) * z2)

c...          Store polygon for solid surface plot:              
c              IF (solid) THEN     
              IF (solid.AND.ang.LT.36.0*PI/180.0) THEN     
                xsur(i2,nsur) = SNGL(p1(1,i2))
                ysur(i2,nsur) = SNGL(p1(2,i2))
                zsur(i2,nsur) = SNGL(p1(3,i2))
              ENDIF

              call transform_vect(mat,p1(1,i2))
              call transform_vect(mat,p2(1,i2))
              CALL POSITN(SNGL(p1(1,i2)),SNGL(p1(2,i2)))
              CALL JOIN  (SNGL(p2(1,i2)),SNGL(p2(2,i2))) 

            ENDDO

          ENDDO

        ENDDO
      ENDIF



      fp1=98
      OPEN(UNIT=fp1,FILE='dump.dat',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=15)

      status = .TRUE.
      DO WHILE (status)

        colind = 0

        READ(fp1,*,END=20) index,x1,y1,z1,x2,y2,z2,colind

        IF (colind.EQ.0) THEN
          IF (index.LT.ascncut) THEN
            CALL LINCOL(ncols+3)
          ELSE
            CALL LINCOL(ncols+1)
          ENDIF
        ELSE
          CALL LINCOL(ncols+colind)
        ENDIF

        p1(1,1) = DBLE(x1) + xshift
        p1(2,1) = DBLE(y1) + yshift
        p1(3,1) = DBLE(z1) + zshift
        p2(1,1) = DBLE(x2) + xshift
        p2(2,1) = DBLE(y2) + yshift
        p2(3,1) = DBLE(z2) + zshift
	
        call transform_vect(mat,p1(1,1))
        call transform_vect(mat,p2(1,1))
	
        CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
        CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 

      ENDDO
      GOTO 20

15    CALL WN('983','Can''t find dump file')

20    CONTINUE

      CLOSE(fp1)


      frac = 1.0
      READ(5,'(A256)') graph1
      IF (graph1(8:11).EQ.'Samp'.OR.graph1(8:11).EQ.'SAMP'.OR.
     .    graph1(8:11).EQ.'samp') THEN
        READ(graph1,*) cdum1,frac
        frac = MIN(frac,1.0,MAX(frac,0.0))
        seed = calc_random_seed(-1)
      ELSE
        BACKSPACE 5
      ENDIF


c...  Plot particle track:
      IF (iopt.NE.0) THEN

        fp1=98
        OPEN(UNIT=fp1,FILE='eirtrac.dat',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=25)

        trackno = 0
        lasttrack = 1

        CALL LINCOL(ncols+6)

        skipnext  = .FALSE.
        skiptrack = .FALSE.

        DO WHILE (.TRUE.)

          READ(fp1,*,ERR=25,END=25) idum1,x1,y1,z1,idum2,idum3,idum4,
     .                              nseg

          IF (idum2.EQ.0) trackno = trackno + 1

          p1(1,1) = x2
          p1(2,1) = y2
          p1(3,1) = z2

          p2(1,1) = x1/100.0D0 + xshift
          p2(2,1) = y1/100.0D0 + yshift
          p2(3,1) = z1/100.0D0 + zshift

c...dev
c          p2(1,1) = x1/100.0D0 * COS(phi / RADDEG)
c          p2(3,1) = x1/100.0D0 * SIN(phi / RADDEG)

c          WRITE(0,*) 'phi=',phi

          phi = dangle * REAL(nseg - 1)

          p2(1,1) = ( COS(phi) * x1 - SIN(phi) * z1) / 100.0D0
          p2(3,1) = (+SIN(phi) * x1 + COS(phi) * z1) / 100.0D0

          x2 = p2(1,1)
          y2 = p2(2,1)
          z2 = p2(3,1)

          IF ((trackno.EQ.iopt.OR.
     .         iopt.EQ.-99.OR.iopt.EQ.-98.OR.
     .         (iopt.LE.-1.AND.idum4.EQ.-iopt)).AND.
     .        idum2.NE.0.AND..NOT.skipnext.AND.
     .                       .NOT.skiptrack) THEN

            IF (idum2.EQ.98) skipnext = .TRUE.

            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))

            IF (iopt.EQ.-98.AND.lasttrack.EQ.trackno) THEN      
              lasttrack = lasttrack + 1
c
c             Reflect particle tracks if the grid has been reflected.  
c
              if (refct.eq.1) then 
                 CALL POSITN (SNGL(p1(1,1))      ,-SNGL(p1(2,1)))
                 CALL JOIN   (SNGL(p1(1,1))+0.001,-SNGL(p1(2,1)))  
              else
                 CALL POSITN (SNGL(p1(1,1))      ,SNGL(p1(2,1)))
                 CALL JOIN   (SNGL(p1(1,1))+0.001,SNGL(p1(2,1)))  
              endif

            ELSEIF (iopt.NE.-98) THEN

c
c             Reflect particle tracks if the grid has been reflected.  
c
              if (refct.eq.1) then 
                 CALL POSITN (SNGL(p1(1,1)),-SNGL(p1(2,1)))
                 CALL JOIN   (SNGL(p2(1,1)),-SNGL(p2(2,1)))  
              else
                 CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                 CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1)))  
              endif

            ENDIF

          ELSEIF (iopt.LT.0.AND.idum2.EQ.0.AND.frac.LT.1.0) THEN

            skiptrack = .FALSE.
            CALL SURAND2 (SEED, 1, RAN)
            IF (ran.GT.frac) skiptrack = .TRUE.
c            WRITE(0,*) 'SKIP:',seed,ran,frac,skiptrack

          ELSEIF (iopt.GE.1.AND.trackno.GT.iopt) THEN

            GOTO 25

          ELSE

            skipnext = .FALSE.

          ENDIF

        ENDDO

      ENDIF

25    CONTINUE

      CLOSE(fp1)


c...  Some comments:     
      IF (iopt.LT.0) THEN
        CALL LINCOL(1)
        WRITE(char(20),'(A,I7  )') 'STRATUM INDEX   =',ABS(iopt)
        if (abs(iopt).eq.1) then 
           WRITE(char(21),'(A)') OUTER//' TARGET SOURCE'
        elseif (abs(iopt).eq.2) then 
           WRITE(char(21),'(A)') INNER//' TARGET SOURCE'
        elseif (abs(iopt).eq.3) then 
           WRITE(char(21),'(A)') 'RECOMBINATION SOURCE'
        elseif (abs(iopt).eq.99) then 
           WRITE(char(21),'(A)') 'ALL PARTICLES'
        endif
        WRITE(char(22),'(A,F7.3)') 'SAMPLE FRACTION =',frac
      ENDIF


c          CALL FILCOL(shade)
c          CALL LINCOL(shade)
c          CALL PTPLOT(rv(1,i1),zv(1,i1),1,nv(i1),1)






c...  Draw scale:
c      CALL PSPACE (0.0, 1.35, 0.0,1.0)
c      CALL CSPACE (0.0, 1.35, 0.0,1.0)
c      CALL MAP    (0.0, 1.35, 0.0,1.0)
c      CALL FULL

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
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
      ENDDO



      CALL FRAME
      CALL LINCOL(1)


c.... Solid surface plot:
      IF (solid) THEN
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     .               YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)

c...    Draw polygons:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)

c...    Position vector:
c        p1(1,1) = 0.0D0
c        p1(2,1) = 0.0D0
c        p1(3,1) = 1.0D0
c        call transform_vect(mat,p1(1,1))
c        posx = p1(1,1)
c        posy = p1(2,1)
c        posz = p1(3,1)

        posx = 0.5 * (cxmin + cxmax)
        posy = 0.5 * (cymin + cymax)
c...    Assumes 90 degree field of view:
        posz = 10.0*0.5 * (cxmax - cxmin)


c...    Determine light shading:
        light(1) = 0.0D0
        light(2) = 1.0D0
        light(3) = 0.0D0

        DO i1 = 1, nsur 
c...      Surface normal:
          a1 = -DBLE(xsur(1,i1) - xsur(2,i1))
          a2 = -DBLE(ysur(1,i1) - ysur(2,i1))
          a3 = -DBLE(zsur(1,i1) - zsur(2,i1))
          b1 =  DBLE(xsur(3,i1) - xsur(2,i1))
          b2 =  DBLE(ysur(3,i1) - ysur(2,i1))
          b3 =  DBLE(zsur(3,i1) - zsur(2,i1))
          normv(1) = a2 * b3 - a3 * b2                
          normv(2) = a3 * b1 - a1 * b3
          normv(3) = a1 * b2 - a2 * b1                

c...      Find angle between surface and light vector:
          theta = SNGL(DACOS((light(1) * normv(1) +
     .                        light(2) * normv(2) +
     .                        light(3) * normv(3)) /
     .            DSQRT(light(1)**2+light(2)**2+light(3)**2) /
     .            DSQRT(normv(1)**2+normv(2)**2+normv(3)**2)))
 
          theta = ABS(theta)

          IF (theta.GE.PI-1.0E-04) THEN
            theta = ABS(theta - PI)
          ENDIF

          lightv(i1) = 0.90 * (1.0 - theta / PI / 2.0) + 0.05

        ENDDO


c...    Rotate all of the surfaces:
        DO i1 = 1, nsur
          DO i2 = 1, 4
            p1(1,1) = DBLE(xsur(i2,i1))
            p1(2,1) = DBLE(ysur(i2,i1))
            p1(3,1) = DBLE(zsur(i2,i1))
            call transform_vect(mat,p1(1,1))
            xsur(i2,i1) = SNGL(p1(1,1))
            ysur(i2,i1) = SNGL(p1(2,1))
            zsur(i2,i1) = SNGL(p1(3,1))
          ENDDO
        ENDDO


c...    Sub-divide surfaces to specified resolution:
        ndiv1 = 2
        ndiv2 = 2
        i1 = 1

        WRITE(0,*) 'NSUR:',nsur

        DO WHILE (.FALSE..AND.i1.LE.nsur)
          
c...      Make room for new surfaces:
          DO i2 = nsur, i1+1, -1
            DO i3 = 1, 4
              xsur(i3,i2+(ndiv1*ndiv2-1)) = xsur(i3,i2)
              ysur(i3,i2+(ndiv1*ndiv2-1)) = ysur(i3,i2)
              zsur(i3,i2+(ndiv1*ndiv2-1)) = zsur(i3,i2)
            ENDDO
            lightv(i2+(ndiv1*ndiv2-1)) = lightv(i2)
          ENDDO

c...      Replicate light vector:
          DO i2 = i1+1, i1+(ndiv1*ndiv2-1)
            lightv(i2) = lightv(i2-1)
          ENDDO

          nsur = nsur + (ndiv1*ndiv2-1)

          IF (nsur.GT.MAXSOLID) 
     .      CALL ER('Solid','Need more space in array',*99)

          WRITE(0,*)

c...      Dupe surface:
          DO i3 = 1, 4
            xsur2(i3) = xsur(i3,i1)
            ysur2(i3) = ysur(i3,i1)
            zsur2(i3) = zsur(i3,i1)

            WRITE(0,*) '--> SURF:',xsur2(i3),ysur2(i3),zsur2(i3)
          ENDDO

c...      Process surfaces:
          DO i2 = 1, ndiv1

            frac1 = DBLE(i2-1) / DBLE(ndiv1)
            frac2 = DBLE(i2  ) / DBLE(ndiv1)

c...        Find intermediate vectors:
            v1x = xsur2(1) + frac1 * (xsur2(2) - xsur2(1))
            v1y = ysur2(1) + frac1 * (ysur2(2) - ysur2(1))
            v1z = zsur2(1) + frac1 * (zsur2(2) - zsur2(1))
            v2x = xsur2(4) + frac1 * (xsur2(3) - xsur2(4))
            v2y = ysur2(4) + frac1 * (ysur2(3) - ysur2(4))
            v2z = zsur2(4) + frac1 * (zsur2(3) - zsur2(4))
					     		
            v3x = xsur2(1) + frac2 * (xsur2(2) - xsur2(1))
            v3y = ysur2(1) + frac2 * (ysur2(2) - ysur2(1))
            v3z = zsur2(1) + frac2 * (zsur2(2) - zsur2(1))
            v4x = xsur2(4) + frac2 * (xsur2(3) - xsur2(4))
            v4y = ysur2(4) + frac2 * (ysur2(3) - ysur2(4))
            v4z = zsur2(4) + frac2 * (zsur2(3) - zsur2(4))

            WRITE(0,*)
            WRITE(0,*) '--> FR12:',frac1,frac2
            WRITE(0,*) '--> VECS:',v1x,v1y,v1z
            WRITE(0,*) '-->      ',v2x,v2y,v2z
            WRITE(0,*) '--> VECS:',v3x,v3y,v3z
            WRITE(0,*) '-->      ',v4x,v4y,v4z

            DO i3 = 1, ndiv2 
              frac1 = DBLE(i3-1) / DBLE(ndiv2)
              frac2 = DBLE(i3  ) / DBLE(ndiv2)

c...          Assign sub-surfaces:
              xsur(1,i1) = v1x + frac1 * (v2x - v1x)
              ysur(1,i1) = v1y + frac1 * (v2y - v1y)
              zsur(1,i1) = v1z + frac1 * (v2z - v1z)

              xsur(2,i1) = v3x + frac1 * (v4x - v3x)
              ysur(2,i1) = v3y + frac1 * (v4y - v3y)
              zsur(2,i1) = v3z + frac1 * (v4z - v3z)

              xsur(3,i1) = v3x + frac2 * (v4x - v3x)
              ysur(3,i1) = v3y + frac2 * (v4y - v3y)
              zsur(3,i1) = v3z + frac2 * (v4z - v3z)

              xsur(4,i1) = v1x + frac2 * (v2x - v1x)
              ysur(4,i1) = v1y + frac2 * (v2y - v1y)
              zsur(4,i1) = v1z + frac2 * (v2z - v1z)

              WRITE(0,*)
              WRITE(0,*) '--> FR12:',frac1,frac2,lightv(i1)
              WRITE(0,*) '--> SUR2:',xsur(1,i1),ysur(1,i1),zsur(1,i1)
              WRITE(0,*) '-->     :',xsur(2,i1),ysur(2,i1),zsur(2,i1)
              WRITE(0,*) '-->     :',xsur(3,i1),ysur(3,i1),zsur(3,i1)
              WRITE(0,*) '-->     :',xsur(4,i1),ysur(4,i1),zsur(4,i1)
              
              i1 = i1 + 1
            ENDDO



          ENDDO

c          STOP 'sdfsdf'

        ENDDO

        WRITE(0,*) 'NSUR:',nsur


c...    Estimate surface center-point:
        DO i1 = 1, nsur 
          xcen(i1) = 0.0
          ycen(i1) = 0.0
          zcen(i1) = 0.0
          DO i2 = 1, 4
            xcen(i1) = xcen(i1) + 0.25 * xsur(i2,i1)
            ycen(i1) = ycen(i1) + 0.25 * ysur(i2,i1)
            zcen(i1) = zcen(i1) + 0.25 * zsur(i2,i1)
          ENDDO
        ENDDO

c...    Sort surfaces:
        DO i1 = 1, nsur
          dsur(i1) = SQRT((posx-xcen(i1))**2 + (posy-ycen(i1))**2 + 
     .                    (posz-zcen(i1))**2)
        ENDDO

        status = .TRUE.
        DO WHILE (status)
          status = .FALSE.
          DO i1 = 1, nsur-1
            IF (dsur(i1).GT.dsur(i1+1)) THEN
              dtemp = dsur(i1)
              cxtmp = xcen(i1)
              cytmp = ycen(i1)
              cztmp = zcen(i1)
              lightvtmp = lightv(i1)
              dsur(i1) = dsur(i1+1)
              xcen(i1) = xcen(i1+1)
              ycen(i1) = ycen(i1+1)
              zcen(i1) = zcen(i1+1)
              lightv(i1) = lightv(i1+1)
              dsur(i1+1) = dtemp
              xcen(i1+1) = cxtmp
              ycen(i1+1) = cytmp
              zcen(i1+1) = cztmp
              lightv(i1+1) = lightvtmp
              DO i2 = 1, 4
                xtemp = xsur(i2,i1)
                ytemp = ysur(i2,i1)
                ztemp = zsur(i2,i1)
                xsur(i2,i1) = xsur(i2,i1+1)
                ysur(i2,i1) = ysur(i2,i1+1)
                zsur(i2,i1) = zsur(i2,i1+1)
                xsur(i2,i1+1) = xtemp
                ysur(i2,i1+1) = ytemp
                zsur(i2,i1+1) = ztemp
             ENDDO
             status = .TRUE.
            ENDIF
          ENDDO
        ENDDO



c...    Plot surfaces:
        DO i1 = 1, nsur 
          CALL RGB
          CALL ColSet(lightv(i1),lightv(i1),lightv(i1),255)          
          CALL LINCOL(255) 
          CALL FILCOL(255)
          CALL PTPLOT(xsur(1,i1),ysur(1,i1),1,4,1)
        ENDDO


c...    Find vector of light source:
        light(1) = 0.0D0
        light(2) = 1.0D0
        light(3) = 0.0D0

c        DO i1 = 1, nsur 
        DO i1 = 1, 0 

c...      Surface normal:
          a1 = -DBLE(xsur(1,i1) - xsur(2,i1))
          a2 = -DBLE(ysur(1,i1) - ysur(2,i1))
          a3 = -DBLE(zsur(1,i1) - zsur(2,i1))
          b1 = DBLE(xsur(3,i1) - xsur(2,i1))
          b2 = DBLE(ysur(3,i1) - ysur(2,i1))
          b3 = DBLE(zsur(3,i1) - zsur(2,i1))
          normv(1) = a2 * b3 - a3 * b2                
          normv(2) = a3 * b1 - a1 * b3
          normv(3) = a1 * b2 - a2 * b1                

c...      Find angle between surface and light vector:
          theta = SNGL(DACOS((light(1) * normv(1) +
     .                        light(2) * normv(2) +
     .                        light(3) * normv(3)) /
     .            DSQRT(light(1)**2+light(2)**2+light(3)**2) /
     .            DSQRT(normv(1)**2+normv(2)**2+normv(3)**2)))
 
          theta = ABS(theta)

          IF (theta.GE.PI-1.0E-04) THEN
            theta = ABS(theta - PI)
          ENDIF

          frac = 0.90 * (1.0 - theta / PI / 2.0) + 0.05
        
c          WRITE(0,'(A,2F14.8,2(2X,3F10.4),2X,6F10.4)') 
c     .       'FRAC:',frac,theta*180.0/PI,
c     .                    normv(1),normv(2),normv(3),
c     .                    light(1),light(2),light(3),
c     .              a1,a2,a3,b1,b2,b3
 
c           WRITE(0,'(3(A,3F10.3))')
c     .       'VECT:',xsur(1,i1),xsur(2,i1),xsur(3,i1),
c     .       'VECT:',ysur(1,i1),ysur(2,i1),ysur(3,i1),
c     .       'VECT:',zsur(1,i1),zsur(2,i1),zsur(3,i1)
         
          CALL RGB
          WRITE(0,*) '???',i1,frac,lightv(i1)

          CALL ColSet(frac,frac,frac,255)          
          CALL LINCOL(255) 
          CALL FILCOL(255)
c          CALL LINCOL(55) 
c          CALL FILCOL(55)


c          DO i2 = 1, 4
c            p1(1,1) = DBLE(xsur(i2,i1))
c            p1(2,1) = DBLE(ysur(i2,i1))
c            p1(3,1) = DBLE(zsur(i2,i1))
c            call transform_vect(mat,p1(1,1))
c            xsur(i2,i1) = SNGL(p1(1,1))
c            ysur(i2,i1) = SNGL(p1(2,1))
c            zsur(i2,i1) = SNGL(p1(3,1))
c          ENDDO

          CALL PTPLOT(xsur(1,i1),ysur(1,i1),1,4,1)

        ENDDO



        CALL FRAME
        CALL LINCOL(1)
      ENDIF

      RETURN
 9012 FORMAT(1X,'PLOT',I3,4X,A)
99    STOP
      END
