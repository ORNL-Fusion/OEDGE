c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE ShiftInversion(inv,xpt,ypt,ir,direction)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_comgra
      use mod_slcom
      use mod_slout
      IMPLICIT none   
c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'comgra'
c     INCLUDE 'slcom'
c     INCLUDE 'slout'      

      INTEGER, PARAMETER :: MAXINVDATA = 1000000
      TYPE :: type_inversion
        INTEGER :: n
        REAL*8  :: x   (MAXINVDATA)
        REAL*8  :: y   (MAXINVDATA)
        REAL*8  :: data(MAXINVDATA)
        INTEGER :: npts
        REAL    :: xpts(4,MAXINVDATA)
        REAL    :: ypts(4,MAXINVDATA)
        REAL    :: xmin
        REAL    :: xmax
        REAL    :: ymin
        REAL    :: ymax
      ENDTYPE type_inversion
      TYPE(type_inversion) :: inv 
      INTEGER ir,direction
      REAL    xpt,ypt

      INTEGER n,i1,i2,nx,ny,npts,imax,id,ik
      LOGICAL tag(MAXINVDATA)
      REAL    x1,x2,y1,y2,xsep,ysep,frac,deltax,deltay,yval,xval,
     .        fmax,xdist,ydist
      REAL*8  x(MAXNKS),y(MAXNKS),f(MAXNKS)

      tag = .FALSE.

      n    = inv%n
      npts = inv%npts


      SELECTCASE (direction)
        CASE (1)  ! Horizontal (and so outboard of x-point)

          nx = MIN(50,MAXNKS)
          ny = 100

          y1 = MAX(inv%ymin,ypt)
          y2 = inv%ymax
          deltay = (y2 - y1) / REAL(ny)

          DO i1 = 1, ny
            x1 = MAX(inv%xmin,xpt)
            x2 = inv%xmax

c            WRITE(0,*) 'i1:',i1
            yval = (REAL(i1-1) + 0.5) * deltay + y1
            y(1:nx) = DBLE(yval)                              
c...        Find distance between x-location of peak and x-value
c           of separatrix:
            DO ik = 1, nks(ir)
              IF (rs(ik,ir).LT.x1.OR.rs(ik,ir).GT.x2.OR.  ! Build list above to save some time...?
     .            zs(ik,ir).LT.y1.OR.zs(ik,ir).GT.y2) CYCLE
              id = korpg(ik,ir)
              IF ((yval.LT.zvertp(1,id).AND.yval.GE.zvertp(4,id)).OR.
     .            (yval.GE.zvertp(1,id).AND.yval.LT.zvertp(4,id))) THEN
                frac = (yval         - zvertp(1,id)) / 
     .                 (zvertp(4,id) - zvertp(1,id))
                xsep = (1.0 - frac) * rvertp(1,id) + frac * rvertp(4,id)
c                WRITE(0,*) 'ik,ir,frac,xsep:',ik,ir,frac,xsep
              ENDIF
            ENDDO
c...        Adjustment is not performed if no magnetic information found:
            IF (xsep.EQ.0.0) CYCLE
c...        Narrow x-range when looking for peak in radial slice of inversion profile:
            x1 = MAX(x1,xsep-0.10)
            x2 = MIN(x2,xsep+0.10)
            deltax = (x2 - x1) / REAL(nx)
            DO i2 = 1, nx
              x(i2) = DBLE((REAL(i2-1) + 0.5) * deltax) + x1  
            ENDDO                                             
c...        Interpolate:
c            WRITE(0,*) 'x:',x(1:nx)
c            WRITE(0,*) 'y:',y(1:nx)
            CALL Interpolate(inv%n,inv%x(1:n),inv%y(1:n),inv%data(1:n),
     .                       nx,x,y,f,0)
c            WRITE(0,*) 'f:',f(1:nx)
c...        Find peak in interpolation:     
            imax = 0
            fmax = -1.0E+20
            DO i2 = 1, nx
              IF (f(i2).GT.fmax) THEN  ! Shorthand..?
                fmax = SNGL(f(i2))
                imax = i2
              ENDIF
            ENDDO
c            WRITE(0,*) 'PEAK:',imax,fmax
c...        Scan through inversion data and shift so that the peak in the
c           radial profile corresponds to the separatrix:
            xdist = x(imax) - xsep
c            WRITE(0,*) 'XDIST:',x(imax),xsep,xdist
            DO i2 = 1, n
              IF (SNGL(inv%y(i2)).GE.yval-0.5*deltay.AND.
     .            SNGL(inv%y(i2)).LT.yval+0.5*deltay.AND.
     .            SNGL(inv%x(i2)).GT.xpt+xdist) THEN

                IF (tag(i2)) CALL ER('LoadInversionData','Bad TAG',*99)
                inv%x   (       i2) = inv%x   (       i2) - DBLE(xdist)
                inv%xpts(1:npts,i2) = inv%xpts(1:npts,i2) -      xdist
c                WRITE(0,*) 'moving:',i2
                tag(i2) = .TRUE.
              ENDIF
            ENDDO

          ENDDO
c
c       ----------------------------------------------------------------
        CASE (2)  ! Vertical (upper, inner)

          nx = 100 ! 100
          ny = MIN(50,MAXNKS)

          x1 = inv%xmin
c          x1 = MAX(inv%xmin,rp(idds(ir,1))-0.01)
          x2 = MIN(inv%xmax,xpt)
          deltax = (x2 - x1) / REAL(nx)

          DO i1 = 1, nx
            y1 = (inv%ymin + ypt) / 2.0  ! Try for now...
            y2 =  inv%ymax

c            WRITE(0,*) 'i1:',i1
            xval = (REAL(i1-1) + 0.5) * deltax + x1
            x(1:ny) = DBLE(xval)                              
c            WRITE(0,*) 'xval:',xval
c...        Find distance between x-location of peak and x-value
c           of separatrix:
            ysep = 0.0
            frac = 0.0
            DO ik = 1, nks(ir)
              IF (rs(ik,ir).LT.x1.OR.rs(ik,ir).GT.x2.OR.  ! Build list above to save some time...?
     .            zs(ik,ir).LT.y1.OR.zs(ik,ir).GT.y2) CYCLE
              id = korpg(ik,ir)
              IF ((xval.LT.rvertp(1,id).AND.xval.GE.rvertp(4,id)).OR.
     .            (xval.GE.rvertp(1,id).AND.xval.LT.rvertp(4,id))) THEN
                frac = (xval         - rvertp(1,id)) / 
     .                 (rvertp(4,id) - rvertp(1,id))
                ysep = (1.0 - frac) * zvertp(1,id) + frac * zvertp(4,id)
c                WRITE(0,*) 'ik,ir,frac,ysep:',ik,ir,frac,ysep
              ENDIF
            ENDDO
c            WRITE(0,*) 'final ik,ir,frac,ysep:',ik,ir,frac,ysep
c...        Adjustment is not performed if no magnetic information found:
            IF (ysep.EQ.0.0) CYCLE
c...        Narrow y-range when looking for peak in radial slice of inversion profile:
            y1 = MAX(y1,ysep-0.10)
            y2 = MIN(y2,ysep+0.10)
            deltay = (y2 - y1) / REAL(ny)
            DO i2 = 1, ny
              y(i2) = DBLE((REAL(i2-1) + 0.5) * deltay) + y1  
            ENDDO                                             
c...        Interpolate:
c            IF (i1.EQ.25) THEN
c              WRITE(0,*) 'x:',x(1:ny)
c              WRITE(0,*) 'y:',y(1:ny)
c            ENDIF
            CALL Interpolate(inv%n,inv%x(1:n),inv%y(1:n),inv%data(1:n),
     .                       ny,x,y,f,0)
c            IF (i1.EQ.25) WRITE(0,*) 'f:',f(1:ny)
c...        Find peak in interpolation:     
            imax = 0
            fmax = -1.0E+20
            DO i2 = 1, ny
              IF (f(i2).GT.fmax) THEN  ! Shorthand..?
                fmax = SNGL(f(i2))
                imax = i2
              ENDIF
            ENDDO
c            WRITE(0,*) 'PEAK:',imax,fmax
c...        Scan through inversion data and shift so that the peak in the
c           radial profile corresponds to the separatrix:
            ydist = y(imax) - ysep
c            WRITE(0,*) 'YDIST:',y(imax),ysep,ydist
            DO i2 = 1, n
              IF (SNGL(inv%x(i2)).GE.xval-0.5*deltax.AND.
     .            SNGL(inv%x(i2)).LT.xval+0.5*deltax.AND.
     .            SNGL(inv%y(i2)).GT.(inv%ymin+ypt)/2.0) THEN  

                IF (tag(i2)) CALL ER('LoadInversionData','Bad TAG',*99)
                inv%y   (       i2) = inv%y   (       i2) - DBLE(ydist)
                inv%ypts(1:npts,i2) = inv%ypts(1:npts,i2) -      ydist
c                WRITE(0,*) 'moving:',i2
                tag(i2) = .TRUE.
              ENDIF
            ENDDO

          ENDDO


        CASE DEFAULT
          CALL ER('ShiftInversion','Unknown orientation',*99)
      ENDSELECT


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadInversionData(val,fname,scale,location)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_comgra
      use mod_slcom
      use mod_slout
      IMPLICIT none   
c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'comgra'
c     INCLUDE 'slcom'
c     INCLUDE 'slout'

      INTEGER       location
      REAL          val(MAXNKS,MAXNRS),scale
      CHARACTER*(*) fname

      INTEGER   fp,i1,i2,i3,n,npts,ik,ir,iklist(MAXNKS),nx
      REAL      rdum1
      REAL*8    x(MAXNKS),y(MAXNKS),f(MAXNKS)
      CHARACTER buffer*512

      REAL   , ALLOCATABLE :: xpts(:,:),ypts(:,:),
     .                        save_xpts(:,:),save_ypts(:,:)
      REAL*8 , ALLOCATABLE :: dinv(:),xcen(:),ycen(:)
      
      INTEGER, PARAMETER :: MAXINVDATA = 1000000
      TYPE :: type_inversion
        INTEGER :: n
        REAL*8  :: x   (MAXINVDATA)
        REAL*8  :: y   (MAXINVDATA)
        REAL*8  :: data(MAXINVDATA)
        INTEGER :: npts
        REAL    :: xpts(4,MAXINVDATA)
        REAL    :: ypts(4,MAXINVDATA)
        REAL    :: xmin
        REAL    :: xmax
        REAL    :: ymin
        REAL    :: ymax
      ENDTYPE type_inversion
      TYPE(type_inversion) :: inv 

c...  Plotting:
      REAL qmin,qmax


c...  Load inversion data:
      fp = 99
      OPEN(fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)
      READ(fp,*,ERR=97) n
      ALLOCATE(dinv(n))
      ALLOCATE(xpts(4,n))
      ALLOCATE(ypts(4,n))
      ALLOCATE(save_xpts(4,n))
      ALLOCATE(save_ypts(4,n))
      ALLOCATE(xcen(n))
      ALLOCATE(ycen(n))
      DO i1 = 1, n
        READ(fp,'(A512)') buffer
        READ(buffer,*) i2,rdum1,npts
        READ(buffer,*) i2,rdum1,npts,(xpts(i3,i2),ypts(i3,i2),i3=1,npts)
        dinv(i2) = DBLE(rdum1)
        xcen(i2) = DBLE(SUM(xpts(1:npts,i2))/REAL(npts))
        ycen(i2) = DBLE(SUM(ypts(1:npts,i2))/REAL(npts))
        save_xpts(1:npts,i2) = xpts(1:npts,i2)
        save_ypts(1:npts,i2) = ypts(1:npts,i2)
      ENDDO
      CLOSE(fp)

      inv%n         = n 
      inv%x   (1:n) = xcen(1:n)
      inv%y   (1:n) = ycen(1:n)
      inv%data(1:n) = dinv(1:n)
      inv%npts      = npts
      DO i1 = 1, npts
        inv%xpts(i1,1:n) = xpts(i1,1:n)
        inv%ypts(i1,1:n) = ypts(i1,1:n)
      ENDDO

c...  Determine spatial domain of the inversion:
      inv%xmin =  1.0E+20
      inv%xmax = -1.0E+20
      inv%ymin =  1.0E+20
      inv%ymax = -1.0E+20
      DO i1 = 1,inv%n
        DO i2 = 1, inv%npts
          inv%xmin = MIN(inv%xmin,inv%xpts(i2,i1))  ! Some shorthand please...
          inv%xmax = MAX(inv%xmax,inv%xpts(i2,i1))
          inv%ymin = MIN(inv%ymin,inv%ypts(i2,i1))
          inv%ymax = MAX(inv%ymax,inv%ypts(i2,i1))
        ENDDO
      ENDDO
c...  Adjust inversion to match equilibrium (that is, the peak emission
c     is assumed to correspond to the separatrix):
      SELECTCASE (location)
        CASE (3)  ! HU12 -- assuming UDND for now...
          ir = irsep

          CALL ShiftInversion(inv,rxp,zxp,ir,1)
          CALL ShiftInversion(inv,rxp,zxp,ir,2)

        CASE DEFAULT
          CALL ER('LoadInversionData','Unknown location',*99)
      ENDSELECT

c...  Interpolate:     
      DO ir = irsep, nrs  ! SOL only for now...
        IF (idring(ir).EQ.BOUNDARY) CYCLE
        nx = 0
        DO ik = 1, nks(ir)
          IF (rs(ik,ir).GT.inv%xmin.AND.rs(ik,ir).LT.inv%xmax.AND.
     .        zs(ik,ir).GT.inv%ymin.AND.zs(ik,ir).LT.inv%ymax) THEN
            nx = nx + 1
            x(nx) = DBLE(rs(ik,ir))
            y(nx) = DBLE(zs(ik,ir))
            iklist(nx) = ik
          ENDIF
        ENDDO

        IF (ir.EQ.irsep+1) THEN
c          WRITE(0,*) 'ir:',ir
c          WRITE(0,*) 'min/max:',inv%xmin,inv%xmax,inv%ymin,inv%ymax
c          WRITE(0,*) 'iklist:',iklist(1:nx)
c          WRITE(0,*) 'x:',x(1:nx)
c          WRITE(0,*) 'y:',y(1:nx)
        ENDIF

        IF (nx.GT.0) THEN
          CALL Interpolate(inv%n,inv%x(1:n),inv%y(1:n),inv%data(1:n),
     .                   nx,x,y,f,0)

          IF (ir.EQ.irsep+1) THEN
c            WRITE(0,*) 'f:',f(1:nx)
          ENDIF
        ENDIF

        DO i1 = 1, nx
           val(iklist(i1),ir) = MAX(0.0,SNGL(f(i1))) * scale
c           IF (ir.EQ.irsep) WRITE(0,*) 
c     .      'val:',iklist(i1),val(iklist(i1),ir)
        ENDDO        

      ENDDO
   








c...  Plot original and shifted inversions:
      IF (.TRUE.) THEN
        qmin =  0.0
        qmax = -HI
        cxmin = inv%xmin
        cxmax = inv%xmax
        cymin = inv%ymin
        cymax = inv%ymax
        DO i1 = 1, n
          qmax = MAX(0.0,MAX(qmax,SNGL(inv%data(i1))))
        ENDDO
        cxmin = cxmin - 0.05 * (cxmax - cxmin)
        cxmax = cxmax + 0.05 * (cxmax - cxmin)
        cymin = cymin - 0.05 * (cymax - cymin)
        cymax = cymax + 0.05 * (cymax - cymin)
        map1x = 0.10             
        map1y = 0.52
        IF (cxmax-cxmin.GT.cymax-cymin) THEN
          map2x = map1x + 0.45
          map2y = map1y + 0.45 * (cymax - cymin) / (cxmax - cxmin)
        ELSE
          map2x = map1x + 0.45 * (cxmax - cxmin) / (cymax - cymin)
          map2y = map1y + 0.45 
        ENDIF
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
c...    Original inversion:
        DO i1 = 1, n
          CALL SetCol255_04(2,MAX(0.0,SNGL(inv%data(i1))),qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255) 
          CALL PTPLOT(save_xpts(1,i1),save_ypts(1,i1),1,inv%npts,1)
        ENDDO
        CALL Supimp('PARTIAL')
        CALL DrawFrame
c...    Shifted:
        map1y = map1y - 0.50
        map2y = map2y - 0.50
        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (cxmin,cxmax,cymin,cymax)
        DO i1 = 1, n
          CALL SetCol255_04(2,MAX(0.0,SNGL(inv%data(i1))),qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255) 
          CALL PTPLOT(inv%xpts(1,i1),inv%ypts(1,i1),1,inv%npts,1)
        ENDDO
c...    Frame:
        CALL Supimp('PARTIAL')
        CALL DrawFrame
        CALL Frame
      ENDIF

c...  Clear arrays:
      IF (ALLOCATED(dinv) ) DEALLOCATE(dinv )
      IF (ALLOCATED(xpts)) DEALLOCATE(xpts)
      IF (ALLOCATED(ypts)) DEALLOCATE(ypts)
      IF (ALLOCATED(save_xpts)) DEALLOCATE(save_xpts)
      IF (ALLOCATED(save_ypts)) DEALLOCATE(save_ypts)
      IF (ALLOCATED(xcen)) DEALLOCATE(xcen)
      IF (ALLOCATED(ycen)) DEALLOCATE(ycen)

      RETURN
 97   CALL ER('LoadInversionData','Problem reading file',*99)
 98   CALL ER('LoadInversionData','Unable to open file',*99)
 99   WRITE(0,*) '    FILE:',fname(1:LEN_TRIM(fname))
      STOP
      END
c
c
c
      SUBROUTINE Plot978(graph,nplots,ref,title,iopt,pltmins,pltmaxs,
     .                   iplot,job)
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_slcom
      use mod_slout
      IMPLICIT none

c     INCLUDE 'params'
c     INCLUDE 'cgeom'
c     INCLUDE 'comtor'
c     INCLUDE 'pindata'
c     INCLUDE 'slcom'
c     INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure

      INTEGER ik,ir

      integer iplot
c      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*80,JOB*72,GRAPH*80,GRAPH1*80
c
      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36  pltlabs(maxplts)
c
      CHARACTER*36 XLAB,YLAB

      CHARACTER*36 REF
      CHARACTER*36 NAME,ELABS(MAXNGS)
c
      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,NPLOTS,M,ID,JR
c
      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact
c
c      REAL tauii,pii
      LOGICAL status
c      INTEGER nenum,tenum,opt_const,plot_mode(30),array,iter1,iter2,
c     .        xaxis,ring,mode,inorm(MAXNGS)
c
c      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
c     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
c     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
c      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
c     .        max1,max2,ynorm(MAXNGS)
c
c
      character graph2*80,graph3*80,graph4*80,cdum1*128


      REAL GetCs,GetJsat

      INTEGER probe,idum1
      REAL    pu
      INTEGER ikp,iks,ike,nfram,ikm

      INTEGER    MAXTDAT     ,MAXCOLS
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)

      INTEGER thnnraw
      REAL thndat(MAXNKS,MAXNRS),thtdat(MAXNKS,MAXNRS),
     .     thnraw(MAXTDAT,MAXCOLS)



       REAL psin(2,MAXNRS)



c...  978: 
      REAL IonViscosity

      INTEGER   offset,ik1,ik2,optflow,ndata,i,line,optscale,ir1
      LOGICAL   shiftthetav
      CHARACTER spec*5,fname*1025,dummy*256
      REAL      intg1,intg2,dmax,gmax,thetav,xdata(MAXNRS),scale,
     .          dmax2,gmax2,fact,
     .          ydata(MAXNRS,5),t,rdum1,osmtmp (MAXNKS,MAXNRS),
     .                                  osmtmp2(MAXNKS,MAXNRS),
     .                                  osmtmp3(MAXNKS,MAXNRS)

      REAL, ALLOCATABLE :: tmppinion(:,:)
      REAL, ALLOCATABLE :: tmposmcfp(:,:)

      COMMON /RSEPCOM/ rsep
      REAL             rsep

c
c     Array for extra comments 
c
      character*20 extra_comments(1)
c
     


c      WRITE(0,*) '978: ',iopt
      WRITE(6,*) '978: ',iopt

      iopt_ghost = 1
      grm_opt    = 1
      slopt      = 1
      slopt4     = 0
      optscale   = 0

c...  Plot the symbol as well as the line:
      READ(5,'(A80)',END=10) graph1
      IF (graph1(8:11).EQ.'Mark'.OR.graph1(8:11).EQ.'MARK'.OR.
     .    graph1(8:11).EQ.'makr') THEN
        READ(graph1,*) cdum1,slopt4
      ELSE
        BACKSPACE 5
      ENDIF
10    CONTINUE

      CALL RSet(thndat,MAXNKS*MAXNRS,LO)
      CALL RSet(thtdat,MAXNKS*MAXNRS,LO)

      CALL RZero(mvals,MAXNKS*MAXPLTS*MAXNGS)
      CALL RZero(mouts,MAXNKS*MAXPLTS)
      CALL RZero(mwids,MAXNKS*MAXPLTS)
c      CALL RZero(machdat,2*MAXNRS)

      WRITE(0,*) 'NPLTS:',nplts,maxplts
      nplts = 0
      DO ip = 1, maxplts
        WRITE (pltlabs(ip),*) ' '
      ENDDO

c      WRITE(0,*) '"',graph,'"'

      READ(graph(14:15),*) xtype
      READ(graph(17:18),*) ytype

      probe = FSP1

      IF (ytype.NE.3) THEN
        IF (prb_num(probe).EQ.0) THEN
          CALL ReadProbeData(status)
          IF (status) THEN
            CALL InterpolateProbeData(probe)
          ELSE
            DO ir = irsep, irwall-1
              ikm = nks(ir) / 2 + 1
              prp_ne(ir,probe) = knbs (ikm,ir)
              prp_te(ir,probe) = ktebs(ikm,ir)
              prp_ti(ir,probe) = ktibs(ikm,ir)
              dip_v (ir,probe) = kvhs (ikm,ir) / qt
            ENDDO
          ENDIF
        ELSE
          CALL InterpolateProbeData(probe)
        ENDIF
      ENDIF



      WRITE(0,*) '978: XTYPE,YTYPE= ',xtype,ytype 
 
c ...change to an array...
      IF     (xtype.EQ.1) THEN
        XLAB = '   s (m) '
      ELSEIF (xtype.EQ.2) THEN
        XLAB = '   theta  '
      ELSEIF (xtype.EQ.3) THEN
        XLAB = '   poloidal dist (m)'
      ELSEIF (xtype.EQ.4) THEN
        XLAB = '   s / smax'
      ELSEIF (xtype.EQ.5) THEN
        XLAB = '   iteration'
      ELSEIF (xtype.EQ.6) THEN
        XLAB = '   psi-n'

c...    Load psi-n data from OUT input file (this will over-ride data from 
c       the .RAW file):
        READ(5,'(A80)') graph1
        IF (graph1(8:11).EQ.'Psin'.OR.graph1(8:11).EQ.'PSIN'.OR.
     .      graph1(8:11).EQ.'psin') THEN
          READ(graph1,*) cdum1,idum1
          IF (idum1.EQ.irwall-2) THEN
            DO i1 = 1, irwall-2
              READ(5,*) cdum1,idum1,psin(IKLO,idum1),psin(IKHI,idum1)
c              WRITE(0,*) '978: PSIN=',idum1,psin(IKLO,idum1),
c     .                                      psin(IKHI,idum1)
            ENDDO
          ELSE
            CALL ER('978','Incorrect number of PSI-N entries',*99)
          ENDIF
        ELSE
          BACKSPACE 5
c...      Set PSIn from PSITARG:
          DO ir = irsep, nrs
            IF (idring(ir).EQ.BOUNDARY) CYCLE
            psin(IKLO,ir) = psitarg(ir,2)
            psin(IKHI,ir) = psitarg(ir,1)
          ENDDO            
        ENDIF

      ELSEIF (xtype.EQ.7) THEN
        XLAB = '   rho (m)'
        CALL CalculateRho

c...    Load rho data from OUT input file (this will over-ride data from 
c       the .RAW file):
        READ(5,'(A80)') graph1
        IF (graph1(8:11).EQ.'Rho '.OR.graph1(8:11).EQ.'RHO '.OR.
     .      graph1(8:11).EQ.'rho ') THEN
          READ(graph1,*) cdum1,idum1
          IF (idum1.EQ.irwall-2) THEN
            DO i1 = 1, irwall-2
              READ(5,*) cdum1,idum1,rho(idum1,CELL1)
            ENDDO
          ELSE
            CALL ER('978','Incorrect number of RHO entries',*99)
          ENDIF
        ELSE
          BACKSPACE 5
        ENDIF

      ELSEIF (xtype.EQ.8) THEN
        WRITE(0,*) '978: 3 XTYPE=8 WORKS FOR HIGH IK TARGET ONLY'
        XLAB = '   r (m)'
c
c     jdemod 
c
      ELSEIF (xtype.EQ.9) THEN
c
        XLAB = 'Separatrix Distance (R-Rsep if flat)'

c...    Load sepdist2 data into the appropriate locations in the PSIN
c       array - this option reuses the code for the psin plots - shouldn't
c       make things any more confusing than they already are ...
c  
c
c       WARNING: ALL OF THESE PLOTS WILL ONLY WORK PROPERLY FOR X-POINT
c                DOWN SONNET STYLE GRIDS - The OUTER/INNER definitions
c                are messed up for other grids. 
c

        do ir = irsep,nrs
c
           if (ir.le.irwall) then 
              psin(IKLO,ir) = sepdist2(idds(ir,2))
              psin(IKHI,ir) = sepdist2(idds(ir,1))
           else
              psin(IKLO,ir) = -sepdist2(idds(ir,2))
              psin(IKHI,ir) = -sepdist2(idds(ir,1))
           endif
c
        end do
c
c     jdemod
c
      ELSEIF (xtype.EQ.10) THEN
c...    10 - MapArray called, s inverted to outer target is s=0
        XLAB = '   s (m) '
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)
c            kss(ik,ir) = ksmaxs(ir) - kss(ik,ir)
c          ENDDO
c        ENDDO

      ELSE
        CALL ER('978','Invalid XTYPE',*99)
      ENDIF

      IF     (ytype.EQ.1) THEN
        call rdg4(graph4,ngrm,nplts,ringnos,maxplts,pltfact,ierr)
        ylab = ' '

        elabs2(1,1) = '    Te'
        elabs2(1,2) = '    Ti'
        ngs2  (1)   = 2

        elabs2(2,1) = '    n'
        ngs2  (2)   = 1

        elabs2(3,1) = '    g'
        ngs2  (3)   = 1
c        elabs2(4,1) = '    M'
c        ngs2  (4)   = 1

        elabs2(4,1) = '    ion'
        elabs2(4,2) = '    rec'
        elabs2(4,3) = '    cfp'
        elabs2(4,4) = '    net'
        ngs2  (4)   = 4

        elabs2(5,1) = '    cx '
        elabs2(5,2) = '    rec'
        elabs2(5,3) = '    osm'
        elabs2(5,4) = '    p  '
        ngs2  (5)   = 4

        elabs2(6,1) = '    Qe '
        elabs2(6,2) = '    -Pei'
        elabs2(6,3) = '    cfe'
        elabs2(6,4) = '    net'
        ngs2  (6)   = 4

        elabs2(7,1) = '    Qi '
        elabs2(7,2) = '    cfi'
        elabs2(7,3) = '    net'
        ngs2  (7)   = 3

        elabs2(8,1) = '    Qecv'
        elabs2(8,2) = '    Qecd'
        elabs2(8,3) = '    Qicv'
        elabs2(8,4) = '    Qicd'
        elabs2(8,5) = '    net'
        ngs2  (8)   = 5

        btype = 1
        nfram = 8

        CALL SetPlotComments(978,job,extra_comments,0,0.0)

c...    Turn on coloured lines on plots:
        slopt2 = 2
        DO i1 = 1, 8
          plottype2(i1,1) = 2
          plottype2(i1,2) = 3
          plottype2(i1,3) = 4
          plottype2(i1,4) = 6
          plottype2(i1,5) = 7
        ENDDO

      ELSEIF (ytype.EQ.2.OR.(ytype.GE.4.AND.ytype.LE.7)) THEN
c...    EIRENE BGK relaxation:

        IF (ytype.EQ.2) THEN
          WRITE(0,*) '978:2 SUB-OPTION NOT SUPPORTED (REASSIGNED TO '//
     .               'SUB-OPTION 4)'
          ytype = 4
        ENDIF

        IF (ytype.EQ.4) SPEC = 'D'
        IF (ytype.EQ.5) SPEC = 'D2'   
        IF (ytype.EQ.6) SPEC = 'DD2'
        IF (ytype.EQ.7) SPEC = 'D2D'

        offset = ytype - 4

        ylab = ' '
        elabs2(1,1) = '    RATN ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (1)   = 1
        elabs2(2,1) = '    RATM1 ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (2)   = 1
        elabs2(3,1) = '    RATM2 ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (3)   = 1
        elabs2(4,1) = '    RATM3 ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (4)   = 1
        elabs2(5,1) = '    RATE ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (5)   = 1
        elabs2(6,1) = '    RESN ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (6)   = 1
        elabs2(7,1) = '    RESE ('//SPEC(1:LEN_TRIM(SPEC))//')'
        ngs2  (7)   = 1
        elabs2(8,1) = '    blank'
        ngs2  (8)   = 1
        btype = 2
        nfram = 8
        nplts = 1
        ngrm  = 8
        pltfact = 0.0

      ELSEIF (ytype.EQ.3) THEN
c...    DIVIMP target data:
        ylab = ' '

        ylab2 (1)   = 'Inner target T (eV)'
        elabs2(1,1) = '    Te'
        elabs2(1,2) = '    Ti'
        ngs2  (1)   = 1

        ylab2 (2)   = 'Outer target T (eV)'
        elabs2(2,1) = '    Te'
        elabs2(2,2) = '    Ti'
        ngs2  (2)   = 1

        ylab2 (3)   = 'Inner target Jsat (A m-2)'
        elabs2(3,1) = '    Jsat'
        ngs2  (3)   = 1

        ylab2 (4)   = 'Outer target Jsat (A m-2)'
        elabs2(4,1) = '    Jsat'
        ngs2  (4)   = 1

        btype   = 2
        nfram   = 4
        nplts   = 1
        ngrm    = 4
        grm_opt = 1

c...    Plot the location of the separatrix:
        IF     (xtype.EQ.6) THEN
          slopt3 = 4
        ELSEIF (xtype.EQ.7.or.xtype.eq.9) THEN
          slopt3 = 5
        ELSEIF (xtype.EQ.8) THEN
          slopt3 = 6
          rsep = rvertp(4,korpg(nks(irsep),irsep))
        ENDIF
        pltfact = 0.0

        CALL SetPlotComments(978,job,extra_comments,0,0.0)

      ELSEIF (ytype.EQ.8) THEN
c...    Average neutral density on the half-ring:
        ylab = ' '

        ylab2 (1)   = 'Inner target T (eV)'
        elabs2(1,1) = '    Inner SOL'
        ngs2  (1)   = 1

        ylab2 (2)   = 'Outer target T (eV)'
        elabs2(2,1) = '    Outer SOL'
        ngs2  (2)   = 1
        btype   = 2
        nfram   = 2
        nplts   = 1
        ngrm    = 2
        grm_opt = 1
        IF     (xtype.EQ.6) THEN
          slopt3 = 4
        ELSEIF (xtype.EQ.7.or.xtype.eq.9) THEN
          slopt3 = 5
        ELSEIF (xtype.EQ.8) THEN
          slopt3 = 6
          rsep = rvertp(4,korpg(nks(irsep),irsep))
        ENDIF
        pltfact = 0.0

        CALL SetPlotComments(978,job,extra_comments,0,0.0)

      ELSEIF (ytype.EQ.9.OR.ytype.EQ.11) THEN
        call rdg4(graph4,ngrm,nplts,ringnos,maxplts,pltfact,ierr)
        ylab = ' '

        ylab2 (1)   = '(eV)'
        elabs2(1,1) = '    Te'
        elabs2(1,2) = '    Ti'
        ngs2  (1)   = 2

        ylab2 (2)   = '(m-3)'
        elabs2(2,1) = '    n'
        elabs2(2,2) = '    Bgamma(norm)'
        ngs2  (2)   = 2

        ylab2 (3)   = '(s-1 m-2)'
        elabs2(3,1) = '    parallel flux'
        ngs2  (3)   = 1

        elabs2(4,1) = '    Mach no.'
        ngs2  (4)   = 1

        ylab2 (5)   = '(m-3 s-1)'
        elabs2(5,1) = '    ion'
        elabs2(5,2) = '    osm'
        elabs2(5,3) = '    rec'
        elabs2(5,4) = '    cfp'
        elabs2(5,5) = '    net'
        ngs2  (5)   = 5


c        ylab2 (6)   = '(Pa)'
c        elabs2(6,1) = '    osm'
c        elabs2(6,2) = '    p  '
c        ngs2  (6)   = 2
        ylab2 (6)   = '(eV m-3)'
c        ylab2 (6)   = '(Pa)'
        elabs2(6,1) = '    pin'
        elabs2(6,2) = '    osm'
        elabs2(6,3) = '    rec'
        elabs2(6,4) = '    p  '
        elabs2(6,5) = '    int'
        ngs2  (6)   = 5


c        ngs2  (6)   = 0

        IF (ytype.EQ.11) THEN
          elabs2(2,2) = '    nD'
          elabs2(2,3) = '    nD2'
          ngs2  (2)   = 3

          elabs2(5,1) = '    ion'
          elabs2(5,2) = '    osm'
          elabs2(5,3) = '    rec'
          elabs2(5,4) = '    osm'
          elabs2(5,5) = '    cfp'
          elabs2(5,6) = '    net'
          elabs2(5,7) = '    dft'
          ngs2  (5)   = 6
          CALL CalcRadialFlux(osmtmp3)
 
          ylab2 (7)   = ' Dalpha'
          elabs2(7,1) = '    Dalpha'
c          ylab2 (7)   = '(m s-1)'
c          elabs2(7,1) = '    v_b'
          ngs2  (7)   = 1

c          ylab2 (8)   = ' Dgamma'
c          elabs2(8,1) = '    EIRENE'
c          ngs2  (8)   = 1

          ylab2 (8)   = ' power'
          elabs2(8,1) = '    AN'
          elabs2(8,2) = '    QE'
          ngs2  (8)   = 2

          scale = 0.0
          osmtmp  = 0.0
          osmtmp2 = 0.0

 11       READ(5,'(A256)') dummy
          IF (dummy(8:10).EQ.'Inv'.OR.dummy(8:10).EQ.'inv'.OR.
     .        dummy(8:10).EQ.'INV') THEN
            WRITE(fname,'(1024X)')
            READ(dummy,*) cdum1,line,optscale,scale,fname
            SELECTCASE (line)
              CASE (1) ! Dalpha
                CALL LoadInversionData(osmtmp2,fname,scale,3)
                ngs2(7) = ngs2(7) + 1
                elabs2(7,2) = '    camera (scaled)'
              CASE (2) ! Dgamma
                CALL LoadInversionData(osmtmp ,fname,scale,3)
                ngs2(8) = ngs2(8) + 1
c                DO ir = 1, nrs
c                 IF (idring(ir).EQ.BOUNDARY) CYCLE
c                  DO ik = 1, nks(ir)
c                    IF (osmtmp (ik,ir).NE.0.0.AND.
c     .                  osmtmp2(ik,ir).NE.0.0) THEN
c                      osmtmp(ik,ir) = osmtmp(ik,ir) / osmtmp2(ik,ir)
c                    ELSE
c                      osmtmp(ik,ir) = 0.0
c                    ENDIF
c                  ENDDO
c                ENDDO
                elabs2(8,2) = '    camera (scaled)'
              CASE DEFAULT
                CALL ER('987','Unknown line',*99)
            ENDSELECT
            GOTO 11
          ELSE
            BACKSPACE 5
          ENDIF
c...      Load and map toroidal inversion data, if specified:
          READ(5,'(A256)') dummy
          IF     (dummy(8:11).EQ.'File'.OR.dummy(8:11).EQ.'file'.OR.
     .            dummy(8:11).EQ.'FILE') THEN
            WRITE(fname,'(1024X)')
            READ(dummy,*) cdum1,line,optscale,scale,fname
            CALL LoadCameraData(osmtmp,fname,scale)
            ngs2(8) = ngs2(8) + 1
            elabs2(8,2) = '    camera (scaled)'
          ELSE
            BACKSPACE 5
          ENDIF
c...      Load and map toroidal inversion data, if specified:
          READ(5,'(A256)') dummy
          IF     (dummy(8:11).EQ.'File'.OR.dummy(8:11).EQ.'file'.OR.
     .            dummy(8:11).EQ.'FILE') THEN
            WRITE(fname,'(1024X)')
            READ(dummy,*) cdum1,line,optscale,scale,fname
            CALL LoadCameraData(osmtmp2,fname,scale)
            elabs2(7,1) = '    EIRENE'
            ngs2(7) = ngs2(7) + 1
            elabs2(8,2) = '    camera (scaled)'
          ELSE
            BACKSPACE 5
          ENDIF

c...      Check if parallel plasma flow should be over-written:
          optflow = 0
          READ(5,'(A80)') graph1
          IF (graph1(8:11).EQ.'Flow'.OR.graph1(8:11).EQ.'FLOW'.OR.
     .        graph1(8:11).EQ.'flow') THEN
            READ(graph1,*) cdum1,optflow

c...        Store ionisation source, PINION:
            ALLOCATE(tmppinion(MAXNKS,MAXNRS))
            ALLOCATE(tmposmcfp(MAXNKS,MAXNRS))
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                tmppinion(ik,ir) = pinion(ik,ir)
                tmposmcfp(ik,ir) = osmcfp(ik,ir)
                pinion(ik,ir) = osmion(ik,ir)
              ENDDO
            ENDDO
            osmcfp = 0.0
            CALL CalcFlow(optflow)
          ELSE
            BACKSPACE 5
          ENDIF


        ELSE
          elabs2(7,1) = '    Qe '
          elabs2(7,2) = '    Qi '
          IF (osmpmk(1,1).NE.LO) THEN
            elabs2(7,3) = '    -mock'
          ELSE
            elabs2(7,3) = '    -Pei'
          ENDIF
          elabs2(7,4) = '    cf'
          elabs2(7,5) = '    Oe'
          elabs2(7,6) = '    net'
          ngs2  (7)   = 6

          elabs2(8,1) = '    Qcd'
          elabs2(8,2) = '    Qcve'
          elabs2(8,3) = '    Qcvi'
          elabs2(8,4) = '    net'
          ngs2  (8)   = 4
        ENDIF

        btype = 1
        nfram = 8

        CALL SetPlotComments(978,job,extra_comments,0,0.0)

c...    Check if parallel plasma flow should be over-written:
        optflow = 0
        READ(5,'(A80)') graph1
        IF (graph1(8:11).EQ.'Flow'.OR.graph1(8:11).EQ.'FLOW'.OR.
     .      graph1(8:11).EQ.'flow') THEN
          READ(graph1,*) cdum1,optflow

c...      Store ionisation source, PINION:
          ALLOCATE(tmppinion(MAXNKS,MAXNRS))
          ALLOCATE(tmposmcfp(MAXNKS,MAXNRS))
          DO ir = 1, nrs
            DO ik = 1, nks(ir)
              tmppinion(ik,ir) = pinion(ik,ir)
              tmposmcfp(ik,ir) = osmcfp(ik,ir)
c              pinion(ik,ir) = osmion(ik,ir)
            ENDDO
          ENDDO
          CALL CalcFlow(optflow)
          osmcfp(1,1) = 0.0
        ELSE
          BACKSPACE 5
        ENDIF

c...    Turn on coloured lines on plots:
        slopt2 = 2
        DO i1 = 1, 8
          plottype2(i1,1) = 2
          plottype2(i1,2) = 3
          plottype2(i1,3) = 4
          plottype2(i1,4) = 6
          plottype2(i1,5) = 7
        ENDDO

      ELSEIF (ytype.EQ.10) THEN
c...    Radial plasma profiles:
        ylab = ' '

        ylab2 (1)   = ' ... '
        elabs2(1,1) = '    Te'
        elabs2(1,2) = '    Ti'
        ngs2  (1)   = 2

        ylab2 (2)   = ' ... '
        elabs2(2,1) = '    n'
        ngs2  (2)   = 1

        ylab2 (3)   = ' velocity'
        elabs2(3,1) = '    v'
        ngs2  (3)   = 1

        ylab2 (4)   = ' pressure'
        elabs2(4,1) = '    p'
        ngs2  (4)   = 1

        btype   = 2
        nfram   = 4
        nplts   = 1
        ngrm    = 4
        grm_opt = 1

c...    Assemble data to plot:
        ndata = 0
        ir = irtrap 
        DO WHILE (ir.NE.irwall-1) 
c...      Assign THETAV:
          thetav = 201.5
c          thetav = 197.5

          ir = irouts(nks(ir)/2,ir)
          ik2 = 0
          shiftthetav = .FALSE.
          DO ik1 = 1, nks(ir)
            IF (.NOT.shiftthetav.AND.
     .          ir.GT.irtrap.AND.ik1.GE.ikti2(ir)) THEN
c...          Shift THETAV because the intersection is in the high 
c             index PFZ:
              thetav = thetav - dthetg
              shiftthetav = .TRUE.
            ENDIF
            IF (thetav.GT.thetag(ik1,ir)) ik2 = ik1
          ENDDO

c...      Assume that THETAV is not valid for this ring:
          IF (ik2.EQ.0.OR.ik2.EQ.nks(ir)) CYCLE

          t = (thetav           - thetag(ik2,ir)) / 
     .        (thetag(ik2+1,ir) - thetag(ik2,ir))

c...      Assign data:
          ndata = ndata + 1
          xdata(ndata) = rho(ir,CELL1)
          ydata(ndata,1) = (1.0-t) * ktebs(ik2,ir) + t * ktebs(ik2+1,ir) 
          ydata(ndata,2) = (1.0-t) * ktibs(ik2,ir) + t * ktibs(ik2+1,ir) 
          ydata(ndata,3) = (1.0-t) * knbs (ik2,ir) + t * knbs (ik2+1,ir) 
          ydata(ndata,4) = (1.0-t) * kvhs (ik2,ir) + t * kvhs (ik2+1,ir) 
          ydata(ndata,4) = ydata(ndata,4) / qt
          ydata(ndata,5) = CalcPressure
     .                       (ydata(ndata,3),ydata(ndata,1),
     .                        ydata(ndata,2),ydata(ndata,4)) * ECH
        ENDDO

        WRITE(0,*) 'NDATA:',ndata


c...    Plot the location of the separatrix:
c        IF     (xtype.EQ.6) THEN
c          slopt3 = 4
c        ELSEIF (xtype.EQ.7.or.xtype.eq.9) THEN
c          slopt3 = 5
c        ELSEIF (xtype.EQ.8) THEN
c          slopt3 = 6
c          rsep = rvertp(4,korpg(nks(irsep),irsep))
c        ENDIF
c        pltfact = 0.0

        CALL SetPlotComments(978,job,extra_comments,0,0.0)

      ELSE
        CALL ER('978','Invalid YTYPE',*9997)
      ENDIF

c
      sctype = iopt
      if (sctype.lt.1.or.sctype.gt.4) sctype = 1

      NPLOTS = NPLOTS + 1
      WRITE (IPLOT,9012) NPLOTS,REF

      CALL rzero (mvals,MAXNKS*MAXNGS*MAXPLTS)

      DO i2 = 1, nplts

        IF     (ytype.EQ.1.OR.ytype.EQ.9.OR.ytype.EQ.11) THEN
          
c          WRITE(0,*) '978: IR= ',ir

          ir = ringnos(i2)
c...FSP1 hardcoded:
          IF (ir.EQ.0) 
     .      CALL ER('PLOT978','IR error, likely plot count too '//
     .              'high',*99)
          pu = CalcPressure(prp_ne(ir,probe),prp_te(ir,probe),
     .                      prp_ti(ir,probe),dip_v (ir,probe))
          ikp = dip_ik(ir,probe)
          iks = 1
          ike = nks(ir)
          IF (ir.LT.irsep) ike = ike - 1

           IF (ir.LT.irsep) THEN
             WRITE(0,*) 
     .         'IKE:',ike,nks(ir),pinion(ike,ir),pinion(nks(ir),ir)
          ENDIF

          IF (ytype.EQ.9) THEN
c...        Get scaling for Dgamma:
            dmax = LO
            gmax = LO
            DO ik = 1, nks(ir)
              dmax = MAX(dmax,knbs(ik,ir))
              gmax = MAX(gmax,pinline(ik,ir,6,H_BGAMMA))
            ENDDO
          ELSEIF (ytype.EQ.11) THEN
c...        Get scaling for mapped camera data (OUTER TARGET BASED):
            dmax = 1.0
            gmax = 1.0
            dmax2 = 1.0
            gmax2 = 1.0
            IF     (optscale.EQ.0) THEN
            ELSEIF (optscale.EQ.1) THEN
c...          Scale accoring to local peak on outside of ring:
              dmax = LO
              gmax = LO
              DO ik = nks(ir)/2, nks(ir)
                IF (machine.EQ.CMOD) THEN
                  dmax = MAX(dmax,pinline(ik,ir,6,H_BGAMMA)*
     .                            (6.63E-34*3.0E+08)/(4340.0*1.0E-10))
                ELSE
                  dmax = MAX(dmax,pinline(ik,ir,6,H_BGAMMA))
                ENDIF
                gmax = MAX(gmax,osmtmp(ik,ir))
              ENDDO
            ELSEIF (optscale.EQ.2) THEN
c...          Scale accoring to local peak on outside separatrix:
              dmax = LO
              gmax = LO
              ir1 = irsep + 2
              DO ik = nks(ir1)/2, nks(ir1)
                IF (machine.EQ.CMOD) THEN
                  dmax = MAX(dmax,pinline(ik,ir1,6,H_BGAMMA)*
     .                            (6.63E-34*3.0E+08)/(4340.0*1.0E-10))
                ELSE
                  dmax = MAX(dmax,pinline(ik,ir1,6,H_BGAMMA))
                ENDIF
                gmax = MAX(gmax,osmtmp(ik,ir1))
              ENDDO
c *TEMP*
              dmax2 = LO
              gmax2 = LO
              DO ik = nks(ir1)/2, nks(ir1)
c              DO ik = nks(ir1)/2, nks(ir1)-10
                dmax2 = MAX(dmax2,pinline(ik,ir1,6,H_BALPHA))
                gmax2 = MAX(gmax2,osmtmp2(ik,ir1))
              ENDDO
            ELSEIF (optscale.EQ.3) THEN
              DO ik = 1, nks(ir)
                osmtmp(ik,ir) = 2.0 * osmtmp(ik,ir)
              ENDDO
            ELSE
              CALL ER('978','UNKNOWN CAMERA SCALING OPTION',*99)
            ENDIF
          ENDIF
        ELSEIF (ytype.GE.4.AND.ytype.LE.7) THEN
          ir  = irsep
          iks = 1
          ike = eirnres
        ELSEIF (ytype.EQ.3) THEN
          IF (xtype.EQ.8) THEN
c...        Only plot out to irbreak:
            STOP 'DOES NOT WORK'
            WRITE(0,*) '978: 3 IRBREAK= ',irbreak
            iks = 1
            ike = irbreak - 2
          ELSE
            iks = 1
            ike = irwall - 7
c            ike = irwall - 2
          ENDIF
          ir = irsep
        ELSEIF (ytype.EQ.8) THEN
          iks = 1
          ike = irwall - irsep 
          ir = irsep
        ELSEIF (ytype.EQ.10) THEN
          iks = 1
          ike = ndata
        ENDIF

        IF (idring(ir).EQ.-1) CYCLE

        DO ip = 1, nfram

          IF (ytype.EQ.1.OR.ytype.EQ.9) THEN
            IF (osm_model(IKLO,ir).EQ.24) THEN
              grm_shade(IKLO,ip) = kss(ikfluid(IKLO,ir),ir)
c              grm_shade(IKLO,ip) = kss(ikbound(ir,IKLO),ir)
              DO ik = 0, ike
                grm_cell(ik,ip) = ksb(ik,ir)
              ENDDO
            ELSE
              grm_shade(IKLO,ip) = 0.0
            ENDIF
            IF (osm_model(IKHI,ir).EQ.24) THEN
              grm_shade(IKHI,ip) = kss(ikfluid(IKHI,ir),ir)
c              grm_shade(IKHI,ip) = kss(ikbound(ir,IKHI),ir)
              DO ik = 0, ike
                grm_cell(ik,ip) = ksb(ik,ir)
              ENDDO
            ELSE
              grm_shade(IKHI,ip) = 0.0
            ENDIF

          ELSEIF (ytype.EQ.3) THEN
            ir = irtrap
          ELSEIF (ytype.EQ.8) THEN
            ir = nrs
          ENDIF

          IF     (ytype.EQ.1.OR.ytype.EQ.9.OR.ytype.EQ.11) THEN
            WRITE(pltlabs(ip),'(I3)') ir
          ELSEIF (ytype.NE.2) THEN
            pltlabs(ip) = ' '
          ENDIF

          IF     (btype.EQ.2.OR.ir.LT.irsep) THEN
c...        Boundary data is at cell centers, not the target locations:
            in  = 0
            inc = 0

          ELSEIF (btype.EQ.1) THEN
c...        Boundary data includes the target locations:
            in = 1
            inc = 2
            mouts(1,ip) = 0.0

            IF     (xtype.EQ.1.OR.xtype.EQ.10) THEN
              mwids(1      ,ip) = kss(1,ir)
              mouts(ike+inc,ip) = ksmaxs(ir)
              mwids(ike+inc,ip) = ksmaxs(ir) - kss(ike,ir)
            ELSEIF (xtype.EQ.2) THEN
              mouts(1      ,ip) = thetat(idds(ir,2))
              mWIDS(1      ,ip) = thetag(1,ir) - thetat(idds(ir,2))
              mouts(ike+inc,ip) = thetat(idds(ir,1))
              mwids(ike+inc,ip) = thetat(idds(ir,1)) - thetag(ike,ir)
            ELSEIF (xtype.EQ.3) THEN
              mwids(1      ,ip) = kps(1,ir)
              mouts(ike+inc,ip) = kpmaxs(ir)
              mwids(ike+inc,ip) = kpmaxs(ir) - kps(ike,ir)
            ELSEIF (xtype.EQ.4) THEN
              mwids(      1,ip) = kss(1,ir) / ksmaxs(ir)
              mouts(ike+inc,ip) = 1.0
              mwids(ike+inc,ip) = 1.0 - kss(ike,ir) / ksmaxs(ir)
c            ELSEIF (xtype.EQ.10) THEN
c              mouts(1      ,ip) = ksmaxs(ir)
c              mouts(ike+inc,ip) = 0.0
            ELSE
              CALL ER('964','Invalid x-axis data',*9997)
            ENDIF
c
c           IF (xtype.EQ.10) THEN
c             id1 = idds(ir,1)
c             id2 = idds(ir,2)
c           ELSE
             id1 = idds(ir,2)
             id2 = idds(ir,1)
c           ENDIF
c            id1 = idds(ir,2)
c            id2 = idds(ir,1)
c
c...        Assign boundary data (where available):
            IF     (ytype.EQ.1) THEN
              IF     (ip.EQ.1) THEN
                mvals(1      ,ip,1) = kteds(id1)
                mvals(ike+inc,ip,1) = kteds(id2)
                mvals(1      ,ip,2) = ktids(id1)
                mvals(ike+inc,ip,2) = ktids(id2)
                mvals(1      ,ip,3) = LO
                mvals(ike+inc,ip,3) = LO
              ELSEIF (ip.EQ.2) THEN
                mvals(1      ,ip,1) = knds(id1)
                mvals(ike+inc,ip,1) = knds(id2)
                mvals(1      ,ip,2) = LO
                mvals(ike+inc,ip,2) = LO
              ELSEIF (ip.EQ.3) THEN
                mvals(1      ,ip,1) = kvds(id1) * knds(id1)
                mvals(ike+inc,ip,1) = kvds(id2) * knds(id2)
c              ELSEIF (ip.EQ.4) THEN
c                mvals(1      ,ip,1) = kvds(id1) /
c     .            GetCs(kteds(id1),ktids(id1))
c                mvals(ike+inc,ip,1) = kvds(id2) /
c     .            GetCs(kteds(id2),ktids(id2))
              ELSEIF (ip.EQ.4) THEN
                mvals(1      ,ip,1) = pinion(1      ,ir)
                mvals(ike+inc,ip,1) = pinion(ike,ir)
                mvals(1      ,ip,2) = -pinrec(1      ,ir)
                mvals(ike+inc,ip,2) = -pinrec(ike,ir)
                mvals(1      ,ip,3) = osmcfp(1      ,ir)
                mvals(ike+inc,ip,3) = osmcfp(ike,ir)
              ELSEIF (ip.EQ.5) THEN
                mvals(1      ,ip,1) = LO
                mvals(ike+inc,ip,1) = LO
                mvals(1      ,ip,2) = LO
                mvals(ike+inc,ip,2) = LO
                mvals(1      ,ip,3) = LO
                mvals(ike+inc,ip,3) = LO
                mvals(1      ,ip,4) = LO
                mvals(ike+inc,ip,4) = LO
              ELSEIF (ip.EQ.6) THEN
                mvals(1      ,ip,1) = pinqe (1      ,ir)
                mvals(ike+inc,ip,1) = pinqe (ike,ir)
                mvals(1      ,ip,2) =-osmpei(1      ,ir)
                mvals(ike+inc,ip,2) =-osmpei(ike,ir)
                mvals(1      ,ip,3) = osmcfe(1      ,ir)
                mvals(ike+inc,ip,3) = osmcfe(ike,ir)
              ELSEIF (ip.EQ.7) THEN
                mvals(1      ,ip,1) = pinqi (1      ,ir)
                mvals(ike+inc,ip,1) = pinqi (ike,ir)
                mvals(1      ,ip,2) = osmcfi(1      ,ir)
                mvals(ike+inc,ip,2) = osmcfi(ike,ir)
              ELSEIF (ip.EQ.8) THEN
                mvals(1      ,ip,1) = osmcve(1  ,ir)
                mvals(ike+inc,ip,1) = osmcve(ike,ir)
                mvals(1      ,ip,2) = osmcde(1  ,ir)
                mvals(ike+inc,ip,2) = osmcde(ike,ir)
                mvals(1      ,ip,3) = osmcvi(1  ,ir)
                mvals(ike+inc,ip,3) = osmcvi(ike,ir)
                mvals(1      ,ip,4) = osmcdi(1  ,ir)
                mvals(ike+inc,ip,4) = osmcdi(ike,ir)
              ENDIF

            ELSEIF (ytype.EQ.9.OR.ytype.EQ.11) THEN
              IF     (ip.EQ.1) THEN
                mvals(1      ,ip,1) = kteds(id1)
                mvals(ike+inc,ip,1) = kteds(id2)
                mvals(1      ,ip,2) = ktids(id1)
                mvals(ike+inc,ip,2) = ktids(id2)
              ELSEIF (ip.EQ.2) THEN
                mvals(1      ,ip,1) = knds(id1)
                mvals(ike+inc,ip,1) = knds(id2)
                mvals(1      ,ip,2) = LO
                mvals(ike+inc,ip,2) = LO
                mvals(1      ,ip,3) = LO
                mvals(ike+inc,ip,3) = LO
              ELSEIF (ip.EQ.3) THEN
                mvals(1      ,ip,1) = kvds(id1) * knds(id1)
                mvals(ike+inc,ip,1) = kvds(id2) * knds(id2)
              ELSEIF (ip.EQ.4) THEN
                mvals(1      ,ip,1) = kvds(id1) /
     .                          GetCs(kteds(id1),ktids(id1))
                mvals(ike+inc,ip,1) = kvds(id2) /
     .                          GetCs(kteds(id2),ktids(id2))
              ELSEIF (ip.EQ.5) THEN
                mvals(1      ,ip,1) = pinion(1      ,ir)
                mvals(ike+inc,ip,1) = pinion(nks(ir),ir)
                mvals(1      ,ip,2) = osmion(1      ,ir)
                mvals(ike+inc,ip,2) = osmion(nks(ir),ir)
                mvals(1      ,ip,3) = -pinrec(1      ,ir)
                mvals(ike+inc,ip,3) = -pinrec(nks(ir),ir)
                mvals(1      ,ip,4) = -osmrec(1      ,ir)
                mvals(ike+inc,ip,4) = -osmrec(nks(ir),ir)
                mvals(1      ,ip,5) = osmcfp(1      ,ir)
                mvals(ike+inc,ip,5) = osmcfp(nks(ir),ir)
c                mvals(1      ,ip,4) = mvals(1      ,ip,1) +
c     .                                mvals(1      ,ip,2) +
c     .                                mvals(1      ,ip,3)
c                mvals(ike+inc,ip,4) = mvals(nks(ir),ip,1) +
c     .                                mvals(nks(ir),ip,2) +
c     .                                mvals(nks(ir),ip,3)
                mvals(1      ,ip,6) = LO
                mvals(ike+inc,ip,6) = LO
                mvals(1      ,ip,7) = osmtmp3(1      ,ir)-
     .                                osmtmp3(ikouts(1,ir),
     .                                        irouts(1,ir))
                mvals(ike+inc,ip,7) = osmtmp3(nks(ir),ir)-
     .                                osmtmp3(ikouts(nks(ir),ir),
     .                                        irouts(nks(ir),ir))
              ELSEIF (ip.EQ.6) THEN
c                mvals(1      ,ip,1) = LO
c                mvals(ike+inc,ip,1) = LO
c                id = id1
c                mvals(1      ,ip,2) = 
c     .            CalcPressure(knds (id),kteds(id),
c     .                         ktids(id),kvds (id))*ECH
c                id = id2
c                mvals(ike+inc,ip,2) = 
c     .            CalcPressure(knds (id),kteds(id),
c     .                         ktids(id),kvds (id))*ECH

                mvals(1      ,ip,1) = LO
                mvals(ike+inc,ip,1) = LO
                mvals(1      ,ip,2) = LO
                mvals(ike+inc,ip,2) = LO
                mvals(1      ,ip,3) = LO
                mvals(ike+inc,ip,3) = LO

c                fact = ECH
                fact = 1.0
                mvals(1      ,ip,4) = 
     .            CalcPressure(knds (id1),kteds(id1),
     .                         ktids(id1),kvds (id1))*fact
                mvals(ike+inc,ip,4) = 
     .            CalcPressure(knds (id2),kteds(id2),
     .                         ktids(id2),kvds (id2))*fact

                mvals(1      ,ip,5) = 
     .            CalcPressure(knds (id1),kteds(id1),
     .                         ktids(id1),kvds (id1))
                CALL CalcIntegral3(osmmp,1,nks(ir),ir,
     .                             mvals(ike+inc,ip,5),4)
                mvals(ike+inc,ip,5) = mvals(ike+inc,ip,5) / ECH
                mvals(ike+inc,ip,5) = (mvals(ike+inc,ip,5) + 
     .                                 mvals(1      ,ip,5))*fact
              ELSEIF (ip.EQ.7) THEN
                IF (ytype.EQ.11) THEN
c                  mvals(1      ,ip,1) = kvds(id1)
c                  mvals(ike+inc,ip,1) = kvds(id2)
                ENDIF
              ELSEIF (ip.EQ.8) THEN
                IF (ytype.EQ.11) THEN
                  IF (machine.EQ.CMOD) THEN
                   mvals(1      ,ip,1) = pinline(1      ,ir,6,H_BGAMMA)*
     .                               (6.63E-34*3.0E+08)/(4340.0*1.0E-10)
                   mvals(ike+inc,ip,1) = pinline(nks(ir),ir,6,H_BGAMMA)*
     .                               (6.63E-34*3.0E+08)/(4340.0*1.0E-10)
                  ELSE
                    mvals(1      ,ip,1) = osmcfe(1      ,ir)
                    mvals(ike+inc,ip,1) = osmcfe(nks(ir),ir)
                    mvals(1      ,ip,2) = osmqe (1      ,ir)
                    mvals(ike+inc,ip,2) = osmqe (nks(ir),ir)
c                    mvals(1      ,ip,1) = pinline(1      ,ir,6,H_BGAMMA)
c                    mvals(ike+inc,ip,1) = pinline(nks(ir),ir,6,H_BGAMMA)
                  ENDIF
c                  mvals(1      ,ip,2) = osmtmp(1      ,ir) * dmax / gmax
c                  mvals(ike+inc,ip,2) = osmtmp(nks(ir),ir) * dmax / gmax
                ENDIF
              ENDIF
            ELSE
              CALL ER('978','Invalid YTYPE',*9997)
            ENDIF
          ELSE
            CALL ER('964','Invalid boundary data type',*9997)
          ENDIF

          pnks(ip) = (ike - iks + 1) + inc
c
c
c
c
c
          i1 = 0

          DO ik = iks, ike

            i1 = i1 + 1

            IF (ytype.EQ.3.OR.ytype.EQ.8) THEN
              ir = irouts(1,ir)
              id1 = idds(ir,2)
              id2 = idds(ir,1)
            ENDIF

c...        Assign x-axis data:

            IF     (xtype.EQ.1.OR.xtype.EQ.10) THEN
              mOUTS(i1+in,ip) = KSS(IK,IR)
c              WRITE(0,*) 'MOUTS:',in,mouts(1:10,ip)
c              WRITE(0,*) '     :',ir,kss  (1:10,ir)
              mWIDS(i1+in,ip) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
            ELSEIF (xtype.EQ.2) THEN
              mOUTS(i1+in,ip) = thetag(IK,IR)
              mWIDS(i1+in,ip) = 0.0
              IF (IK.GT.1) mWIDS(i1+in,ip)=0.5*(thetag(IK,IR)-
     .                                          thetag(IK-1,IR))
              IF (IK.LT.ike) mWIDS(i1+in,ip) = mWIDS(i1+in,ip) +
     >                           0.5 * (thetag(IK+1,IR)-thetag(IK,IR))
            ELSEIF (xtype.EQ.3) THEN
              mOUTS(i1+in,ip) = KPS(IK,IR)
              mWIDS(i1+in,ip) = 0.0
              IF (IK.GT.1) mWIDS(i1+in,ip)=0.5*(KPS(IK,IR)-KPS(IK-1,IR))
              IF (IK.LT.ike) mWIDS(i1+in,ip) = mWIDS(i1+in,ip) +
     .                             0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
            ELSEIF (xtype.EQ.4) THEN
              mOUTS(i1+in,ip) = KSS(IK,IR) / ksmaxs(ir)
              mWIDS(i1+in,ip) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR)) /
     .                          ksmaxs(ir)
            ELSEIF (xtype.EQ.5) THEN
              mOUTS(i1+in,ip) = REAL(ik)
              mWIDS(i1+in,ip) = 1.0
            ELSEIF (xtype.EQ.6.or.xtype.eq.9) THEN
              IF (ytype.EQ.3) THEN
                IF     (ip.EQ.1.OR.ip.EQ.3) THEN
                  mouts(i1+in,ip) = psin(IKLO,ir)
                  mwids(i1+in,ip) = 1.0
                ELSEIF (ip.EQ.2.OR.ip.EQ.4) THEN
                  mouts(i1+in,ip) = psin(IKHI,ir)
                  mwids(i1+in,ip) = 1.0
                ENDIF
              ELSE
                CALL ER('978','Invalid YTYPE for XTYPE=6',*99)
              ENDIF
            ELSEIF (xtype.EQ.7) THEN
c...          RHO C-Mod coordinate:
              mouts(i1+in,ip) = rho(ir,CELL1)
              mwids(i1+in,ip) = 1.0
            ELSEIF (xtype.EQ.8) THEN
c...          R *OUTER* target plate coordinate:
              mouts(i1+in,ip) = rp(id2)
              mwids(i1+in,ip) = 1.0
            ENDIF
c
c           Assign y-series data:
            IF     (ytype.EQ.1) THEN
c...          Plasma solution combination plot:
              IF     (ip.EQ.1) THEN
                mvals(i1+in,ip,1) = ktebs (ik,ir)
                mvals(i1+in,ip,2) = ktibs (ik,ir)
                mvals(i1+in,ip,3) = thtdat(ik,ir)
              ELSEIF (ip.EQ.2) THEN
                mvals(i1+in,ip,1) = knbs  (ik,ir)
                mvals(i1+in,ip,2) = thndat(ik,ir)
              ELSEIF (ip.EQ.3) THEN
                mvals(i1+in,ip,1) = knbs(ik,ir) * kvhs(ik,ir) / qt
c              ELSEIF (ip.EQ.4) THEN
c                mvals(i1+in,ip,1) = kvhs(ik,ir) / qt /
c     .                        GetCs(ktebs(ik,ir),ktibs(ik,ir))
              ELSEIF (ip.EQ.4) THEN
                IF (pinion(1,1).NE.LO) mvals(i1+in,ip,1) = pinion(ik,ir)
                IF (pinrec(1,1).NE.LO) mvals(i1+in,ip,2) =-pinrec(ik,ir)
                IF (osmcfp(1,1).NE.LO) mvals(i1+in,ip,3) = osmcfp(ik,ir)
                mvals(i1+in,ip,4) = mvals(i1+in,ip,1) +
     .                              mvals(i1+in,ip,2) +
     .                              mvals(i1+in,ip,3)
              ELSEIF (ip.EQ.5) THEN
                mvals(i1+in,ip,1) = pinmp(ik,ir)
                mvals(i1+in,ip,2) = -2.0 * AMU * kvhs(ik,ir) / qt *
     .                               pinrec(ik,ir)
                mvals(i1+in,ip,3) = osmmp(ik,ir)
                mvals(i1+in,ip,4) =
     .            CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                         ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH
              ELSEIF (ip.EQ.6) THEN
                IF (pinqe (1,1).NE.LO) mvals(i1+in,ip,1) = pinqe (ik,ir)
                IF (osmpei(1,1).NE.LO) mvals(i1+in,ip,2) =-osmpei(ik,ir)
                IF (osmcfe(1,1).NE.LO) mvals(i1+in,ip,3) = osmcfe(ik,ir)
                mvals(i1+in,ip,4) = mvals(i1+in,ip,1) +
     .                              mvals(i1+in,ip,2) +
     .                              mvals(i1+in,ip,3)
              ELSEIF (ip.EQ.7) THEN
                IF (pinqi (1,1).NE.LO) mvals(i1+in,ip,1) = pinqi (ik,ir)
                IF (osmcfi(1,1).NE.LO) mvals(i1+in,ip,2) = osmcfi(ik,ir)
                mvals(i1+in,ip,3) = mvals(i1+in,ip,1) +
     .                              mvals(i1+in,ip,2) 
              ELSEIF (ip.EQ.8) THEN
                IF (osmcve(1,1).NE.LO) mvals(i1+in,ip,1) = osmcve(ik,ir)
                IF (osmcde(1,1).NE.LO) mvals(i1+in,ip,2) = osmcde(ik,ir)
                IF (osmcvi(1,1).NE.LO) mvals(i1+in,ip,3) = osmcvi(ik,ir)
                IF (osmcdi(1,1).NE.LO) mvals(i1+in,ip,4) = osmcdi(ik,ir)
c     .                                                   + osmcde(ik,ir)
                mvals(i1+in,ip,5) = mvals(i1+in,ip,1) +
     .                              mvals(i1+in,ip,2) +
     .                              mvals(i1+in,ip,3) + 
     .                              mvals(i1+in,ip,4)
              ENDIF

            ELSEIF (ytype.GE.4.AND.ytype.LE.7) THEN
c...          EIRENE BGK relaxation:
              IF     (ip.EQ.1) THEN
                mvals(i1+in,ip,1) = eirres(1,3+offset,ik)
              ELSEIF (ip.EQ.2) THEN
                mvals(i1+in,ip,1) = eirres(2,3+offset,ik)
              ELSEIF (ip.EQ.3) THEN
                mvals(i1+in,ip,1) = eirres(3,3+offset,ik)
              ELSEIF (ip.EQ.4) THEN
                mvals(i1+in,ip,1) = eirres(4,3+offset,ik)
              ELSEIF (ip.EQ.5) THEN
                mvals(i1+in,ip,1) = eirres(5,3+offset,ik)
              ELSEIF (ip.EQ.6) THEN
                mvals(i1+in,ip,1) = eirres(6,3+offset,ik)
              ELSEIF (ip.EQ.7) THEN
                mvals(i1+in,ip,1) = eirres(7,3+offset,ik)
              ELSEIF (ip.EQ.8) THEN
                mvals(i1+in,ip,1) = 0.0
              ENDIF

            ELSEIF (ytype.EQ.3) THEN
c...          DIVIMP target data:
              IF     (ip.EQ.1) THEN
                mvals(i1+in,ip,1) = kteds(id1)
                mvals(i1+in,ip,2) = ktids(id1)
              ELSEIF (ip.EQ.2) THEN
                mvals(i1+in,ip,1) = kteds(id2)
                mvals(i1+in,ip,2) = ktids(id2)
              ELSEIF (ip.EQ.3) THEN
                mvals(i1+in,ip,1) = -GetJsat(0.0,0.0,knds(id1),
     .                                               kvds(id1))
              ELSEIF (ip.EQ.4) THEN
                mvals(i1+in,ip,1) =  GetJsat(0.0,0.0,knds(id2),
     .                                               kvds(id2))
              ENDIF

            ELSEIF (ytype.EQ.8) THEN
c...          Half-ring average neutral densities:
              IF     (ip.EQ.1) THEN
                ik1 = 1
                ik2 = nks(ir) / 2
              ELSEIF (ip.EQ.2) THEN
                ik1 = nks(ir) / 2 + 1
                ik2 = nks(ir)
              ENDIF

              CALL CalcIntegral3(pinatom,ik1,ik2,ir,intg1,6)
              CALL CalcIntegral3(pinmol ,ik1,ik2,ir,intg2,6)

              mvals(i1+in,ip,1) = (intg1 + 2.0 * intg2) / 
     .                            (kpb(ik2,ir) - kpb(ik1-1,ir))

              WRITE(6,*) '978 8 : ',ip,ir,ik1,ik2,intg1,intg2

            ELSEIF (ytype.EQ.9.OR.ytype.EQ.11) THEN
c...          Plasma solution combination plot for Ti=Te:
              IF     (ip.EQ.1) THEN
                mvals(i1+in,ip,1) = ktebs(ik,ir)
                mvals(i1+in,ip,2) = ktibs(ik,ir)
              ELSEIF (ip.EQ.2) THEN
                mvals(i1+in,ip,1) = knbs(ik,ir)
                IF (ytype.EQ.11) THEN
c                  mvals(i1+in,ip,2) = LO
c                  mvals(i1+in,ip,3) = LO
                  mvals(i1+in,ip,2) = pinatom(ik,ir)
                  mvals(i1+in,ip,3) = pinmol(ik,ir)
                ELSE
                  mvals(i1+in,ip,2) = pinline(ik,ir,6,H_BGAMMA) * 
     .                                dmax / gmax
                ENDIF
              ELSEIF (ip.EQ.3) THEN
                mvals(i1+in,ip,1) = knbs(ik,ir) * kvhs(ik,ir) / qt
              ELSEIF (ip.EQ.4) THEN
                mvals(i1+in,ip,1) = kvhs(ik,ir) / qt /
     .                        GetCs(ktebs(ik,ir),ktibs(ik,ir))
              ELSEIF (ip.EQ.5) THEN
c                IF (pinion(1,1).NE.LO) mvals(i1+in,ip,1) = pinion(ik,ir)
c                IF (pinrec(1,1).NE.LO) mvals(i1+in,ip,2) =-pinrec(ik,ir)
c                IF (osmcfp(1,1).NE.LO) mvals(i1+in,ip,3) = osmcfp(ik,ir)
                mvals(i1+in,ip,1) = pinion(ik,ir)
                mvals(i1+in,ip,2) = osmion(ik,ir)
                mvals(i1+in,ip,3) =-pinrec(ik,ir)
                mvals(i1+in,ip,4) =-osmrec(ik,ir)
                mvals(i1+in,ip,5) = osmcfp(ik,ir)


            IF (ir.EQ.irsep.AND.ik.EQ.nks(ir))  
     .        WRITE(0,*) 'FLUX:',osmcfp(nks(irsep),irsep)
c                mvals(i1+in,ip,5) = mvals(i1+in,ip,2) +
c     .                              mvals(i1+in,ip,3) +
c     .                              mvals(i1+in,ip,4)
                mvals(i1+in,ip,6) = LO

                mvals(i1+in,ip,7) = osmcfpflx(ik,ir,1)
c                mvals(i1+in,ip,7) = osmtmp3(ik,ir)
c                mvals(i1+in,ip,7) = osmtmp3(ik,ir)-
c     .                              osmtmp3(ikouts(ik,ir),
c     .                                      irouts(ik,ir))
              ELSEIF (ip.EQ.6) THEN
c                mvals(i1+in,ip,1) = osmmp(ik,ir)
c                mvals(i1+in,ip,2) =
c     .            CalcPressure(knbs (ik,ir),ktebs(ik,ir),
c     .                         ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH

c                fact = ECH
                fact = 1.0
                mvals(i1+in,ip,1) = pinmp(ik,ir)/ECH*fact
                mvals(i1+in,ip,2) = osmmp(ik,ir)/ECH*fact
                mvals(i1+in,ip,3) = -2.0 * AMU * kvhs(ik,ir) / qt *
     .                               pinrec(ik,ir)/ECH*fact
                mvals(i1+in,ip,4) =
     .            CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                         ktibs(ik,ir),kvhs (ik,ir)/qt) * fact

                IF (.TRUE.) THEN
                  CALL CalcIntegral3(osmmp,ik-1,ik,ir,rdum1,9)
                  mvals(i1+in,ip,5) = mvals(MAX(1,i1+in-1),ip,5) + 
     .                                rdum1/ECH*fact
c                  IF (i1.EQ.1) THEN
c                    WRITE(0,*) 'I!=1:',i1+in,mvals(i1+in,ip,5),
c     .                   mvals(i1+in-1,ip,5),rdum1
c                   STOP 'sfsd'
c                  ENDIF

                ELSEIF (.FALSE..AND.
     .              ik.GE.ikbound(ir,IKLO).AND.
     .              ik.LE.ikbound(ir,IKHI)) THEN
                  CALL CalcIntegral3(osmmp,1,ik,ir,mvals(i1+in,ip,5),1)
	      
                  ik1 = ikbound(ir,IKLO)
	      
                  mvals(i1+in,ip,5) = mvals(i1+in,ip,5) +
     .            CalcPressure(knbs (ik1,ir),ktebs(ik1,ir),
     .                         ktibs(ik1,ir),kvhs (ik1,ir)/qt) * ECH
                ELSE
                  mvals(i1+in,ip,5) = 0.0
                ENDIF

                rdum1 = IonViscosity(ik,ir,knbs(ik,ir),ktibs(ik,ir),
     .                  (crmb*AMU),1.0)           


              ELSEIF (ip.EQ.7) THEN
                IF (ytype.EQ.11) THEN
                  mvals(i1+in,ip,1) = pinline(ik,ir,6,H_BALPHA)
                  mvals(i1+in,ip,2) = osmtmp2(ik,ir) * dmax2 / gmax2
c                  mvals(i1+in,ip,1) = kvhs(ik,ir)
                ELSE
c                  IF (pinqe (1,1).NE.LO) mvals(i1+in,ip,1) = pinqe (ik,ir)
c                  IF (pinqi (1,1).NE.LO) mvals(i1+in,ip,2) = pinqi (ik,ir)
c                  IF (osmpei(1,1).NE.LO) mvals(i1+in,ip,3) =-osmpei(ik,ir)
c                  IF (osmpmk(1,1).NE.LO) mvals(i1+in,ip,3) =-osmpmk(ik,ir)
c                  IF (osmcfe(1,1).NE.LO) mvals(i1+in,ip,4) = osmcfe(ik,ir)
c                  IF (osmcfi(1,1).NE.LO) mvals(i1+in,ip,4) = 
c     .                                 mvals(i1+in,ip,4) + osmcfi(ik,ir)
c                  IF (pinqe (1,1).NE.LO) mvals(i1+in,ip,5) = osmqe (ik,ir)

                  mvals(i1+in,ip,1) = pinqe (ik,ir)
                  mvals(i1+in,ip,2) = pinqi (ik,ir)
                  mvals(i1+in,ip,3) =-osmpei(ik,ir)
                  mvals(i1+in,ip,3) =-osmpmk(ik,ir)
                  mvals(i1+in,ip,4) = osmcfe(ik,ir)
                  mvals(i1+in,ip,4) = mvals(i1+in,ip,4) + osmcfi(ik,ir)
                  mvals(i1+in,ip,5) = osmqe (ik,ir)
                  mvals(i1+in,ip,6) = mvals(i1+in,ip,1) +
     .                                mvals(i1+in,ip,2) +
     .                                mvals(i1+in,ip,3) +
     .                                mvals(i1+in,ip,4) +
     .                                mvals(i1+in,ip,5) 
                  ENDIF
              ELSEIF (ip.EQ.8) THEN
                IF (ytype.EQ.11) THEN
                  IF (machine.EQ.CMOD) THEN
                    mvals(i1+in,ip,1) = pinline(ik,ir,6,H_BGAMMA)*
     .                              (6.63E-34*3.0E+08)/(4340.0*1.0E-10)
                  ELSE
                    mvals(i1+in,ip,1) = osmcfe(ik,ir)
                    mvals(i1+in,ip,2) = osmqe (ik,ir)
c                    mvals(i1+in,ip,1) = pinline(ik,ir,6,H_BGAMMA) 
c                    mvals(i1+in,ip,1) = pinline(ik,ir,6,H_BGAMMA) /
c     .                                  pinline(ik,ir,6,H_BALPHA) 
                  ENDIF
c                  mvals(i1+in,ip,2) = osmtmp(ik,ir) * dmax / gmax
                ELSE
c                  IF (osmcde(1,1).NE.LO) mvals(i1+in,ip,1) = osmcde(ik,ir)
c                  IF (osmcdi(1,1).NE.LO) mvals(i1+in,ip,1) = 
c     .                                   mvals(i1+in,ip,1) + osmcdi(ik,ir)
c                  IF (osmcve(1,1).NE.LO) mvals(i1+in,ip,2) = osmcve(ik,ir)
c                  IF (osmcvi(1,1).NE.LO) mvals(i1+in,ip,3) = osmcvi(ik,ir)
                  mvals(i1+in,ip,1) = osmcde(ik,ir)
                  mvals(i1+in,ip,1) = mvals(i1+in,ip,1) + osmcdi(ik,ir)
                  mvals(i1+in,ip,2) = osmcve(ik,ir)
                  mvals(i1+in,ip,3) = osmcvi(ik,ir)
                  mvals(i1+in,ip,4) = mvals(i1+in,ip,1) +
     .                                mvals(i1+in,ip,2) +
     .                                mvals(i1+in,ip,3)
                ENDIF
              ENDIF

            ELSEIF (ytype.EQ.10) THEN

c...          Assign independent variable here instead of above (what a mess):
              mouts(i1+in,ip) = xdata(i1)
              mwids(i1+in,ip) = 1.0

              IF     (ip.EQ.1) THEN
                mvals(i1+in,ip,1) = ydata(i1,1)
                mvals(i1+in,ip,2) = ydata(i1,2)
              ELSEIF (ip.EQ.2) THEN
                mvals(i1+in,ip,1) = ydata(i1,3)
              ELSEIF (ip.EQ.3) THEN
                mvals(i1+in,ip,1) = ydata(i1,4)
              ELSEIF (ip.EQ.4) THEN
                mvals(i1+in,ip,1) = ydata(i1,5)
              ENDIF

            ELSE
              CALL ER('978','Invalid YTYPE',*9997)
            ENDIF

          ENDDO
        ENDDO
c
c
c
        WRITE(6,*) 'PLOT 978, YTYPE=',ytype
        DO ip = 1, nfram
          DO i1 = 1, ngs2(ip)
            status = .FALSE.
            DO ik = iks, ike+inc
              IF (mvals(ik,ip,i1).NE.0.0) status = .TRUE.
            ENDDO
            IF (.NOT.status) THEN
              WRITE(6,*) '978: NULL DATA ARRAY DETECTED: IP,I1 = ',ip,i1
              CALL RSet(mvals(1,ip,i1),MAXNKS,LO)
            ENDIF
          ENDDO

c...      Invert order:
          IF (xtype.EQ.10) THEN
c...        Invert x-coordinate:
            DO ik = iks, ike+inc
              mouts(ik,ip) = ksmaxs(ir) - mouts(ik,ip)
            ENDDO
c...        Reorder data points: 
            status = .TRUE.
            DO WHILE (status)
              status = .FALSE.
              DO ik = iks, ike+inc-1
                IF (mouts(ik,ip).GT.mouts(ik+1,ip)) THEN
                  mouts(MAXNKS,ip) = mouts(ik,ip)
                  mwids(MAXNKS,ip) = mwids(ik,ip)
                  DO i1 = 1, ngs2(ip)
                    mvals(MAXNKS,ip,i1) = mvals(ik,ip,i1)
                  ENDDO
                  mouts(ik,ip) = mouts(ik+1,ip)
                  mwids(ik,ip) = mwids(ik+1,ip)
                  DO i1 = 1, ngs2(ip)
                    mvals(ik,ip,i1) = mvals(ik+1,ip,i1)
                  ENDDO
                  mouts(ik+1,ip) = mouts(MAXNKS,ip)
                  mwids(ik+1,ip) = mwids(MAXNKS,ip)
                  DO i1 = 1, ngs2(ip)
                    mvals(ik+1,ip,i1) = mvals(MAXNKS,ip,i1)
                  ENDDO                  
                  status = .TRUE.
                ENDIF
              ENDDO

            ENDDO
          ENDIF


          DO ik = iks, ike+inc
            WRITE(6,'(2I4,2X,2I4,2F8.4,2X,2F8.4,2X,1P,10(E11.3:))')
     .        ik,ir,ip,inc,
     .        mouts(ike+inc,ip),mwids(ike+inc,ip),
     .        mouts(ik,ip),mwids(ik,ip),
     .        (mvals(ik,ip,i1),i1=1,ngs2(ip))
          ENDDO
        ENDDO

c        WRITE(0,*) '978: XTYPE,YTYPE= ',xtype,TRIM(xlab)


        CALL SLDRAWM (mouts,mwids,mvals,MAXNKS,pnks,
     >                nfram,ngs,pltlabs,elabs,xlab,ylab,ref,title,
     >                sctype,nfram,pltmins,pltmaxs,pltfact)


      ENDDO


      IF (optflow.NE.0) THEN
c...    Restore PINION, OSMCFP:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            pinion(ik,ir) = tmppinion(ik,ir)
            osmcfp(ik,ir) = tmposmcfp(ik,ir)
          ENDDO
        ENDDO
        DEALLOCATE(tmppinion,STAT=i)
        DEALLOCATE(tmposmcfp,STAT=i)
      ENDIF


      IF (xtype.EQ.10) THEN
c...    Restore KSS:
c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)
c            kss(ik,ir) = ksmaxs(ir) - kss(ik,ir)
c          ENDDO
c        ENDDO
      ENDIF


 9997 CONTINUE
      RETURN
99    STOP
 9012 FORMAT(1X,'PLOT',I3,4X,A)
      END
