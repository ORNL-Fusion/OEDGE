

c
c ======================================================================
c
c subroutine: BasicPlot01
c
      SUBROUTINE BasicPlot01(ndata,xdata,ydata,pdata,xwid,ywid,
     .                       title,xlab,ylab,col,ncon)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slcom'
      INCLUDE 'slout'
      include 'comgra'
c
c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT
c
      INTEGER   ndata,col,ncon
      REAL      xdata(ndata),ydata(ndata),pdata(ndata),xwid,ywid
      CHARACTER title*(*),xlab*(*),ylab*(*)

      INTEGER i1
      REAL    r1,x1,x2,y1,y2,x(4),y(4),cpval1,cpval2,cpmin,cpmax

      CALL PSPACE (map1x,map2x,map1y,map2y)
      CALL MAP    (cxmin,cxmax,cymin,cymax)


      CALL XScale
      CALL YScale


      cpmin =  HI
      cpmax = -HI
      DO i1 = 1, ndata
        IF (pdata(i1).NE.0.0.AND.pdata(i1).LT.cpmin) cpmin = pdata(i1)
        IF (pdata(i1).NE.0.0.AND.pdata(i1).GT.cpmax) cpmax = pdata(i1)
      ENDDO

      DO i1 = 1, ndata
        x1 = xdata(i1) - xwid
        x2 = xdata(i1) + xwid
        y1 = ydata(i1) - ywid
        y2 = ydata(i1) + ywid

        x(1) = x1
        y(1) = y1
        x(2) = x1
        y(2) = y2
        x(3) = x2
        y(3) = y2
        x(4) = x2
        y(4) = y1

        DO r1 = 0.0, REAL(ncon)-1.0
          cpval1 = cpmax - (REAL(ncon) - r1        ) / REAL(ncon) *
     .                     (cpmax - cpmin)
          cpval2 = cpmax - (REAL(ncon) - (r1 + 1.0)) / REAL(ncon) *
     .                     (cpmax - cpmin)

          IF (pdata(i1).GE.cpval1.AND.pdata(i1).LE.cpval2) THEN
            CALL FILCOL(INT(r1)+1)
            CALL LINCOL(INT(r1)+1)
          ENDIF
        ENDDO

        CALL PTPLOT(x(1),y(1),1,4,1)
      ENDDO


      CALL PCSCEN(0.5*(cxmin+cxmax),cymax*1.10,title(1:LEN_TRIM(title)))
      CALL PCSCEN(0.5*(cxmin+cxmax),cymin*1.10,xlab(1:LEN_TRIM(xlab)))
      CALL PCSCEN(1.10*cxmin,0.5*(cymin+cymax),ylab(1:LEN_TRIM(ylab)))

      CALL PSPACE (0.0, 1.35, 0.0, 1.0)
      CALL MAP    (0.0, 1.35, 0.0, 1.0)

      DO r1 = 0.0, REAL(ncon)-1.0
        x1 = map2x + 0.05
        x2 = map2x + 0.08
        y1 = map2y - (REAL(ncon) - r1) / REAL(ncon) * (map2y - map1y)
        y2 = y1 + (map2y - map1y) / REAL(ncon)

        x(1) = x1
        y(1) = y1
        x(2) = x1
        y(2) = y2
        x(3) = x2
        y(3) = y2
        x(4) = x2
        y(4) = y1

        WRITE(0,*) 'basic01',INT(r1)+1,x1,x2,y1,y2

        CALL FILCOL(INT(r1)+1)
        CALL LINCOL(INT(r1)+1)

        CALL PTPLOT(x(1),y(1),1,4,1)
      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c

c
c ======================================================================
c
      SUBROUTINE Plot990(nplts,ringnos,graph,nplots,ref,title,iopt,
     .                   ngrm,pltmins,pltmaxs,pltfact,iplot)
      use mod_comtor
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure

      INTEGER ik,ir

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*80,JOB*72,GRAPH*80,GRAPH1*80

      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,YLAB,XPOINT
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,ELABS(MAXNGS),PLABS(-2:MAXPLRP),KLAB

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,NPLOTS,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact

      REAL tauii,pii
      LOGICAL status
      INTEGER nenum,tenum,opt_const,plot_mode(30),array,iter1,iter2,
     .        xaxis,ring,mode,inorm(MAXNGS)

      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
     .        max1,max2,ynorm(MAXNGS)


      REAL CalcDist
      INTEGER iks,ike,ik1,ndata
      REAL    deltas,dist(MAXNKS,MAXNRS),sum,magcfp,sumn,sump,factp,
     .        factn,pm,pl,xdata(MAXPTS),ydata(MAXPTS),
     .        pdata(MAXPTS,MAXNGS)

      INTEGER    MAXTDAT     ,MAXCOLS   ,MAXNPS     ,MAXTYP
      PARAMETER (MAXTDAT=1000,MAXCOLS=10,MAXNPS=1000,MAXTYP=10)

      INTEGER nraw,ind
      REAL   thodat(MAXNKS,MAXNRS,MAXTYP),raw(MAXTDAT,MAXCOLS)

      INTEGER gndata(MAXNRS,MAXTYP)
      REAL    gydata(MAXNPS,MAXNRS,MAXTYP),gxdata(MAXNPS,MAXNRS,MAXTYP)



      real mvals(maxnps,maxplts,maxngs)
      real mouts (maxnps,maxplts),mwids (maxnps,maxplts)

      REAL xxmin,xxmax,yymin,yymax

      INTEGER xshiftn,yshiftn,qnt,nsum
      REAL    xshifts,xshifte,xshiftd,yshifts,yshifte,yshiftd,x1,y1,
     .        vsum

      iopt_ghost = 1
      slopt      = 1

      CALL RZero(mvals,MAXNPS*MAXPLTS*MAXNGS)
      CALL RZero(mouts,MAXNPS*MAXPLTS)
      CALL RZero(mwids,MAXNPS*MAXPLTS)

      CALL RZero(xdata,MAXPTS)
      CALL RZero(ydata,MAXPTS)
      CALL RZero(pdata,MAXPTS*MAXNGS)


      xshifts = -15.0
      xshifte =  15.0
      xshiftd =   5.0
      xshiftn =  INT((xshifte - xshifts) / xshiftd) + 1

      yshifts = -15.0
      yshifte =  15.0
      yshiftd =   5.0
      yshiftn =  INT((yshifte - yshifts) / yshiftd) + 1

      WRITE(0,*) '990: x,ySHIFT= ',xshiftn,yshiftn

      mode = 1

      qnt  = 1

      DO y1 = yshifts, yshifte, yshiftd
        DO x1 = yshifts, yshifte, yshiftd

         ndata = ndata + 1
         xdata(ndata) = x1
         ydata(ndata) = y1

         CALL RSet(thodat,MAXNKS*MAXNRS*MAXTYP,0.0)

c         CALL LoadThomsonData(thodat(1,1,1),nraw,raw,
c     .                        MAXTDAT,MAXCOLS,1,x1*0.001,y1*0.001,1)
         CALL LoadThomsonData(nraw,raw,MAXTDAT,MAXCOLS,
     .                        x1*0.001,y1*0.001,1)
         STOP '990: NEED TO USE ASSIGN THOMSON DATA'
c         CALL LoadThomsonData(thodat(1,1,2),nraw,raw,
c     .                        MAXTDAT,MAXCOLS,2,x1*0.001,y1*0.001,1)


c...     Convert units of postion data to millimeters:
         DO i1 = 1, nraw
           DO i2 = 1, MAXCOLS
             raw(i1,i2) = raw(i1,i2) * 1000.0
           ENDDO
         ENDDO

         IF     (mode.EQ.1) THEN
c...       Generate statistics for cell averaged Thomson values:
c...       Make sure only compare sol22 rings -- or sol12 or whatever, but
c          not bogus rings like sol6
c...only outer sol for now
           nsum = 0
           vsum = 0.0
           DO ir = irsep, irwall-1
             DO ik = 1, nks(ir)
c...           Temperature:
               IF (thodat(ik,ir,2).NE.0.0) THEN
                 nsum = nsum + 1
c...do some proper statistics:
                 vsum = vsum + ABS(ktebs(ik,ir) - thodat(ik,ir,2))
               ENDIF
             ENDDO
           ENDDO
           IF (nsum.GT.0) THEN
             pdata(ndata,1) = vsum / REAL(nsum)
             pdata(ndata,2) = REAL(nsum)
           ENDIF

c         ELSEIF (mode.EQ.1) THEN
c...        Generate statistics using raw Thomson data:
c           CALL AsgnThomsonData(nraw,raw,qnt,gxdata,gydata,gndata,
c     .                          MAXNPS,MAXTYP,MAXTDAT,MAXCOLS)
         ELSE
           CALL ER('990','Invalid MODE',*99)
         ENDIF
        ENDDO
      ENDDO


      DO i1 = 1, ndata
        WRITE(6,90) i1,xdata(i1),ydata(i1),(pdata(i1,i2),i2=1,4)
90      FORMAT(I4,2F10.4,1P,10E12.4,0P)
      ENDDO


      xxmin = xshifts - xshiftd
      xxmax = xshifte + xshiftd
      yymin = yshifts - yshiftd
      yymax = yshifte + yshiftd

      CALL SLSet(0.15,0.55,0.60,1.00,xxmin,xxmax,yymin,yymax)
      CALL BasicPlot01(ndata,xdata,ydata,pdata(1,1),xshiftd,yshiftd,
     .                 'Plot 1','x','y',1,10)

      CALL SLSet(0.75,1.15,0.60,1.00,xxmin,xxmax,yymin,yymax)
      CALL BasicPlot01(ndata,xdata,ydata,pdata(1,2),xshiftd,yshiftd,
     .                 'Plot 2','x','y',1,10)

      CALL SLSet(0.15,0.55,0.05,0.45,xxmin,xxmax,yymin,yymax)
      CALL BasicPlot01(ndata,xdata,ydata,pdata(1,3),xshiftd,yshiftd,
     .                 'Plot 3','x','y',1,10)

      CALL SLSet(0.75,1.15,0.05,0.45,xxmin,xxmax,yymin,yymax)
      CALL BasicPlot01(ndata,xdata,ydata,pdata(1,4),xshiftd,yshiftd,
     .                 'Plot 4','x','y',1,10)

      CALL FRAME

      RETURN
99    STOP
      END


















      SUBROUTINE PlotXXX(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .                   xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .                   xouts,youts,nconts,conts,cntropt)
      use mod_comtor
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'



      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure

      INTEGER ik,ir

      integer cngs,cntropt

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*80,JOB*72,GRAPH*80,GRAPH1*80

      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,YLAB,XPOINT
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,ELABS(MAXNGS),PLABS(-2:MAXPLRP),KLAB

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,NPLOTS,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact

      REAL tauii,pii
      LOGICAL status
      INTEGER nenum,tenum,opt_const,plot_mode(30),array,iter1,iter2,
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
      real conts(maxpts),clev(maxpts),offset
      INTEGER ipeak

      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)


      CHARACTER*128 dataline,cdum1,cdum2
      CHARACTER*7   pinlabel(-1:2)

      INTEGER plotcode,line,iopt1,iopt2,iopt3
      CHARACTER cname*256,resdir*256,cmnd*256
      LOGICAL oldraw
      CHARACTER graph6*128


c...630:
      INTEGER i,k,ii
      REAL    PLTMAX,PLTMIN
      REAL LOUTS(MAXSEG),LWIDS(MAXSEG),LVALS(MAXSEG,MAXNGS)
      REAL ydata(MAXSEG,MAXNGS)
      integer llabs(maxseg)

      pinlabel(-1) = 'NIMBUS'
      pinlabel( 0) = 'NIMBUS'
      pinlabel( 1) = 'EIRENE'
      pinlabel( 2) = 'EIRENE'


      TABLE  = 'SYMBOL TABLE'

      ngs = 1

c...  IOPT = 1 - 630   PIN - Major radius of vessel wall
c            2 - 631   PIN - Major radius of divertor wall
c            3 - 632   PIN - Height of vessel wall
c            4 - 633   PIN - Height of divertor wall
c            5 - 634   PIN - Neutral flux to vessel wall
c            6 - 635   PIN - Neutral flux to divertor wall
c            7 -
c            8 -

      IF     (iopt.EQ.1) THEN
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'MAJOR RADIUS (M)'
        REF  = 'MAJOR RADIUS OF VESSEL WALL'
        NAME = '    RVESM'
        DO i1 = 1, nvesm
          WRITE(6,'(A,I4,I4,2F10.4,2X,2F10.4)')
     .      'xVESM',i1,jvesm(i1),rvesm(i1,1),zvesm(i1,1),rxp,zxp
        ENDDO
      ELSEIF (iopt.EQ.2) THEN
        NAME = '    RVESM'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'MAJOR RADIUS (M)'
        REF  = 'MAJOR RADIUS OF DIVERTOR WALL'
      ELSEIF (iopt.EQ.3) THEN
      ELSEIF (iopt.EQ.4) THEN
      ELSEIF (iopt.EQ.5) THEN
c...    634:
        NAME = '    FLUXHW'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'HYDROGEN FLUX (M-2 S-1)'
        REF  = 'NEUTRAL FLUX TO VESSEL WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = fluxhw(i1)
        ENDDO
      ELSEIF (iopt.EQ.6) THEN
c...    635:
        NAME     = '    FLUXHW'
        elabs(1) = '    FLUXHW '//pinlabel(pincode)
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'HYDROGEN FLUX (M-2 S-1)'
        REF  = 'NEUTRAL FLUX TO DIVERTOR WALL'

        status = .TRUE.
        oldraw = .FALSE.
        DO WHILE (status)
          DO i1 = 1, MAXSEG
            ydata(i1,ngs) = fluxhw(i1)
          ENDDO
          status = .FALSE.
          READ(5,'(A128)',END=10) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,cname,iopt1,cmnd
            CALL LoadData(cname,cmnd,NULL,iopt1)
            status = .TRUE.
            oldraw = .TRUE.
            ngs = ngs + 1
            elabs(ngs) = '    FLUXHW '//pinlabel(pincode)
          ELSE
            BACKSPACE 5
          ENDIF
10        CONTINUE
        ENDDO
        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

      ELSEIF (iopt.EQ.7) THEN
c...    636:
        NAME = '    FLXHW2'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ION+ATOM FLUX (M-2 S-1)'
        REF  = 'ION+ATOM FLUX TO VESSEL WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw2(i1)
        ENDDO
      ELSEIF (iopt.EQ.8) THEN
c...    637:
        NAME = '    FLXHW2'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ION+ATOM FLUX (M-2 S-1)'
        REF  = 'ION+ATOM FLUX TO DIVERTOR WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw2(i1)
        ENDDO
      ELSEIF (iopt.EQ.9) THEN
c...    642:
        NAME = '    FLXHW5'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'H0 IMPACT ENERGY (EV)'
        REF  = 'ATOM ENERGY ON VESSEL WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw5(i1)
        ENDDO
      ELSEIF (iopt.EQ.10) THEN
c...    643:
        NAME = '    FLXHW5'
        elabs(1) = '    FLXHW5 '//pinlabel(pincode)
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'H0 IMPACT ENERGY (EV)'
        REF  = 'ATOM ENERGY ON DIVERTOR WALL'

        status = .TRUE.
        oldraw = .FALSE.
        DO WHILE (status)
          DO i1 = 1, MAXSEG
            ydata(i1,ngs) = flxhw5(i1)
          ENDDO
          status = .FALSE.
          READ(5,'(A128)',END=30) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,cname,iopt1,cmnd
            CALL LoadData(cname,cmnd,NULL,iopt1)
            status = .TRUE.
            oldraw = .TRUE.
            ngs = ngs + 1
            elabs(ngs) = '    FLUXHW '//pinlabel(pincode)
          ELSE
            BACKSPACE 5
          ENDIF
30        CONTINUE
        ENDDO
        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

      ELSEIF (iopt.EQ.11) THEN
        NAME = '    FLXHW6'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ATOM FLUX TO WALL'
        REF  = 'ATOM FLUX ON VESSEL WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw6(i1)
        ENDDO
      ELSEIF (iopt.EQ.12) THEN
        NAME = '    FLXHW6'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ATOM FLUX TO WALL'
        REF  = 'ATOM FLUX ON DIVERTOR WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw6(i1)
        ENDDO
      ELSEIF (iopt.EQ.13) THEN
        NAME = '    FLXHW2-FLXHW6'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ION FLUX TO WALL'
        REF  = 'ION FLUX ON VESSEL WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw2(i1) - flxhw6(i1)
        ENDDO
      ELSEIF (iopt.EQ.14) THEN
        NAME = '    FLXHW2-FLXHW6'
        XLAB = 'DISTANCE ALONG WALL (M)'
        YLAB = 'ION FLUX TO WALL'
        REF  = 'ION FLUX ON DIVERTOR WALL'
        DO i1 = 1, MAXSEG
          ydata(i1,1) = flxhw2(i1) - flxhw6(i1)
        ENDDO
      ELSEIF (iopt.EQ.15) THEN
      ELSEIF (iopt.EQ.16) THEN
      ELSEIF (iopt.EQ.17) THEN
      ENDIF

      IF     (iopt.EQ.1) THEN
c...    Based on plot 630:
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NPLOTS = NPLOTS + 1

          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = 0.5*(rvesm(1,2)+rvesm(1,1))
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = 0.5*(rvesm(i,2)+rvesm(i,1))
c
            pltmax = max(pltmax,lvals(i,1))
            pltmin = min(pltmin,lvals(i,1))
          enddo
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,
     >      louts(1)-0.5*lwids(1),louts(nvesm)+0.5*lwids(nvesm),
     >      pltmin,pltmax,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL GRTRAC (louts,lvals,nvesm,NAME,'LINE',1)
          call region (nvesm,louts,lwids,jvesm,pltmin,pltmax)
          CALL FRAME
        endif

      ELSEIF (iopt.EQ.2) THEN
c...    Based on plot 631:
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
c  count the number of wall segments in the divertor
c  separating outer and inner
c
          j = 0
          k = 0
          do i = 1, nvesm
c: inner (including private void region)
            if (jvesm(i).ge.4.and.jvesm(i).le.6.or.jvesm(i).eq.8)
     >        j = j + 1
c: outer
            if (jvesm(i).ge.1.and.jvesm(i).le.3)
     >        k = k + 1
          enddo
c
c  work forward through outer divertor segments
c
          lwids(j+1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                     +(zvesm(1,2)-zvesm(1,1))**2)
          louts(j+1) = lwids(j+1) * 0.5
          lvals(j+1,1) = 0.5*(rvesm(1,2)+rvesm(1,1))
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = 0.5*(rvesm(i,2)+rvesm(i,1))
            llabs(ii) = jvesm(i)
c
            pltmax = max(pltmax,lvals(ii,1))
            pltmin = min(pltmin,lvals(ii,1))
          enddo
c
c  work backward through inner divertor segments
c
          lwids(j) = sqrt((rvesm(nvesm,2)-rvesm(nvesm,1))**2
     >                   +(zvesm(nvesm,2)-zvesm(nvesm,1))**2)
          louts(j) = -lwids(j) * 0.5
          lvals(j,1) = 0.5*(rvesm(nvesm,2)+rvesm(nvesm,1))
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = 0.5*(rvesm(i,2)+rvesm(i,1))
            llabs(ii) = jvesm(i)
c
            pltmax = max(pltmax,lvals(ii,1))
            pltmin = min(pltmin,lvals(ii,1))
          enddo
c
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,
     >      louts(1)-0.5*lwids(1),louts(j+k)+0.5*lwids(j+k),
     >      pltmin,pltmax,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL GRTRAC (louts,lvals,j+k,NAME,'LINE',1)
          call region (j+k,louts,lwids,llabs,pltmin,pltmax)
          CALL FRAME
        endif

      ELSEIF (iopt.EQ. 5.OR.iopt.EQ.7.OR.iopt.EQ.9.OR.iopt.EQ.11.OR.
     .        iopt.EQ.13      ) THEN
c...    Based on 634:
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI


          ipeak = 1
          DO i = 1, nvesm
            IF (zvesm(i,1).GT.zvesm(ipeak,1)) ipeak = i
          ENDDO

          offset = 0.0

          i  = ipeak
          ii = 0
          DO j = 1, nvesm
              ii = ii + 1

              lwids(ii) = SQRT((rvesm(i,2) - rvesm(i,1))**2 +
     >                         (zvesm(i,2) - zvesm(i,1))**2)
              IF (ii.GT.1) THEN
                louts(ii) = louts(ii-1)+0.5*(lwids(ii)+lwids(ii-1))
              ELSE
                louts(ii) = 0.5 * lwids(ii) + offset
              ENDIF
              lvals(ii,1) = ydata(i,1)
              llabs(ii)   = jvesm(i)

              pltmax = MAX(pltmax,lvals(ii,1))
              pltmin = MIN(pltmin,lvals(ii,1))

            IF (i.EQ.nvesm) THEN
              i = 1
            ELSE
              i = i + 1
            ENDIF
          ENDDO

          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,
     >      louts(1)-0.5*lwids(1),louts(nvesm)+0.5*lwids(nvesm),
     >      pltmin,pltmax,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL GRTRAC (louts,lvals,nvesm,NAME,'LINE',1)
          call region (nvesm,louts,lwids,llabs,pltmin,pltmax)
          CALL FRAME
        endif

      ELSEIF (iopt.EQ. 6.OR.iopt.EQ.8.OR.iopt.EQ.10.OR.iopt.EQ.12.OR.
     .        iopt.EQ.14      ) THEN
c...    Based on 635:
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI

          ipeak = 1
          DO i = 1, nvesm
            IF (zvesm(i,1).GT.zvesm(ipeak,1)) ipeak = i
          ENDDO


          i  = ipeak
          ii = 0
          DO j = 1, nvesm
c...Condition -- half-way between the x-point and the center point:
            IF (0.5*(zvesm(i,1)+zvesm(i,2)).LT. 0.5*(zxp+z0)) THEN
              ii = ii + 1

              lwids(ii) = SQRT((rvesm(i,2) - rvesm(i,1))**2 +
     >                         (zvesm(i,2) - zvesm(i,1))**2)
              IF (ii.GT.1) THEN
                louts(ii) = louts(ii-1)+0.5*(lwids(ii)+lwids(ii-1))
              ELSE
c...      Find offset, so that distance around the wall coordinate
c         matches full wall plots:
                offset = 0.0
                DO i1 = ipeak, i
                offset = offset + SQRT((rvesm(i1,2) - rvesm(i1,1))**2 +
     >                                 (zvesm(i1,2) - zvesm(i1,1))**2)
                ENDDO

                louts(ii) = 0.5 * lwids(ii) + offset
c                louts(ii) =
              ENDIF
              llabs(ii)   = jvesm(i)

              DO i1 = 1, ngs
                lvals(ii,i1) = ydata(i,i1)
                pltmax = MAX(pltmax,lvals(ii,i1))
                pltmin = MIN(pltmin,lvals(ii,i1))
              ENDDO
            ENDIF

            IF (i.EQ.nvesm) THEN
              i = 1
            ELSE
              i = i + 1
            ENDIF
          ENDDO


          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,
     >      louts(1)-0.5*lwids(1),louts(ii)+0.5*lwids(ii),
     >      pltmin,pltmax,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)

          IF (ngs.GT.1) THEN
            DO i1 = 1, ngs
              CALL GRTRAC (louts,lvals(1,i1),ii,elabs(i1),'LINE',1)
            ENDDO
          ELSE
            CALL GRTRAC (louts,lvals,ii,NAME,'LINE',1)
          ENDIF

          call region (ii,louts,lwids,llabs,pltmin,pltmax)
          CALL FRAME
        endif

      ELSE
        CALL ER('PlotXXX','Unknown plot option',*99)
      ENDIF

      RETURN
99    STOP
9012  FORMAT(1X,'PLOT',I3,4X,A)
      END


