

      SUBROUTINE Plot972(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .                  xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .                  xouts,youts,nconts,conts,cntropt,n_cols,col_opt)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'
      INCLUDE 'colours'


      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure,GetCs

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
      real conts(maxpts),clev(maxpts)
      INTEGER ii
      integer n_cols,col_opt
      real saturation,value,hue

      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)

      INTEGER count

      CHARACTER*128 dataline,cdum1,cdum2





        WRITE(0,*) '972'

        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
c slmod begin
        ref(1:31) = graph(11:41)
c
c        REF = 'PIN - H-ALPHA '// XPOINT
c slmod end
        WRITE (IPLOT,9012) NPLOTS,REF

        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c slmod begin
c        IF (graph(41:44).NE.'    ') THEN
c          READ(graph(41:42),'(I2)') i1
c          READ(graph(43:44),'(I2)') i2
c
c          DO ii = i1, i2
c            CALL DrawHalpha(ii,1)
c          ENDDO
c        ENDIF
c slmod end

        CALL StoreGrid(1)

        CALL LoadCamera(cdata)
       
        slopt = 3

c        CALL slCONTOUR (ICNTR,NGS,cdata,1,1,1,FT,FP,1.0,
c     >                XOUTS,1,NXS,YOUTS,1,NYS,
c     >                XXMIN,XXMAX,YYMIN,YYMAX,
c     >                nconts,conts,cntropt)

      CALL CONTOUR(ICNTR,NGS,cdata,1,1,1,FT,FP,1.0,
     .             XOUTS,1,NXS,YOUTS,1,NYS,
     .             XXMIN,XXMAX,YYMIN,YYMAX,
     .             nconts,conts,cntropt,0.0,0.0)

        CALL StoreGrid(2)

        slopt = 0

        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)

        CALL SUPIMP2('SELECT')

        CALL FRAME





      RETURN
99    STOP
9012  FORMAT(1X,'PLOT',I3,4X,A)
      END

c
c ======================================================================
c
c subroutine: StoreGrid
c
      SUBROUTINE StoreGrid(mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      INTEGER mode

      INTEGER ik,ir,i1,i2

      INTEGER nrs1,irsep1,irwall1,irtrap1,nbr1,nks1(MAXNRS),
     .        nv(  MAXNKS*MAXNRS),ko(  MAXNKS,MAXNRS)
      REAL    rv(5,MAXNKS*MAXNRS),zv(5,MAXNKS*MAXNRS),
     .        r (  MAXNKS,MAXNRS),z (  MAXNKS,MAXNRS)

      SAVE

      IF (mode.EQ.1) THEN

        irsep1  = irsep
        irwall1 = irwall
        irtrap1 = irtrap
        nbr1    = nbr
        nrs1    = nrs

        DO ir = 1, MAXNRS
          nks1(ir) = nks(ir)

          DO ik = 1, MAXNKS
            ko(ik,ir) = korpg(ik,ir)
            r (ik,ir) = rs(ik,ir)
            z (ik,ir) = zs(ik,ir)
          ENDDO
        ENDDO

        DO i1 = 1, MAXNKS*MAXNRS
          nv(i1) = nvertp(i1)

          DO i2 = 1, 5
            rv(i2,i1) = rvertp(i2,i1)
            zv(i2,i1) = zvertp(i2,i1)
          ENDDO
        ENDDO

      ELSE

        irsep  = irsep1
        irwall = irwall1
        irtrap = irtrap1
        nbr    = nbr1
        nrs    = nrs1

        DO ir = 1, MAXNRS
          nks(ir) = nks1(ir)

          DO ik = 1, MAXNKS
            korpg(ik,ir) = ko(ik,ir)
            rs(ik,ir) = r (ik,ir)
            zs(ik,ir) = z (ik,ir)
          ENDDO
        ENDDO

        DO i1 = 1, MAXNKS*MAXNRS
          nvertp(i1) = nv(i1)

          DO i2 = 1, 5
            rvertp(i2,i1) = rv(i2,i1)
            zvertp(i2,i1) = zv(i2,i1)
          ENDDO
        ENDDO

      ENDIF

      RETURN
99    STOP
      END
c
c subroutine: LoadCamera
c
      SUBROUTINE LoadCamera(cdata)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'

      REAL cdata(MAXNKS,MAXNRS)

      INTEGER   fp,in,ik,ir,i1,nrs_camera
      REAL      deltar,deltaz,rdum1,zlimit
      CHARACTER buffer*200

      fp = 69

      OPEN(UNIT=fp,FILE='camera.dat',ACCESS='SEQUENTIAL')

      CALL TransferLine(fp,6,buffer,5)

      READ(buffer,*,END=97,ERR=98) nrs_camera

      nrs = MIN(nrs_camera,MAXNRS)

      IF (nrs.GT.MAXNRS) CALL ER('LoadCamera','R bounds error',*99)

      DO ir = 1, nrs_camera
        READ(fp,*,END=97,ERR=98) rdum1
c...new
c        IF (ir.LE.nrs) rs(1,ir)= rdum1
        IF (ir.LE.nrs) rs(1,ir)= rdum1+0.002
c        IF (ir.LE.nrs) rs(1,ir)= rdum1+0.005
      ENDDO

      CALL TransferLine(fp,6,buffer,2)

      READ(buffer,*) nks(1)

      IF (nks(1).GT.MAXNKS) CALL ER('LoadCamera','K bounds error',*99)

      DO ik = 1, nks(1)
        READ(fp,*,END=97,ERR=98) zs(ik,1)
c...new
c        zs(ik,1) = zs(ik,1) 
        zs(ik,1) = zs(ik,1) + 0.01
      ENDDO

      CALL TransferLine(fp,6,buffer,1)
      DO ir = 1, nrs_camera
        DO ik = 1, nks(1)
          READ(fp,*,END=97,ERR=98) rdum1
c...new
          IF (ir.LE.nrs) cdata(ik,ir) = MAX(0.0,rdum1)
        ENDDO
      ENDDO

      IF (nrs.LT.nrs_camera) THEN
        WRITE(0,*) 'WARNING LoadCamera: Dropping the last',
     .             nrs_camera-nrs,' columns of camera data'
      ENDIF

      CLOSE(fp)
c
c
c
      DO ir = 1, nrs
        nks(ir) = nks(1)
        DO ik = 1, nks(ir)
          rs(ik,ir) = rs(1,ir) 
          zs(ik,ir) = zs(ik,1) 
        ENDDO
      ENDDO
c
c     ... Assumes a regular grid:
c
      deltar = 0.5 * ABS(rs(1,1) - rs(1,2))
      deltaz = 0.5 * ABS(zs(1,1) - zs(2,1))

      in = 0

      DO ir = 1, nrs
        DO ik = 1, nks(ir)

          in = in + 1

          korpg(ik,ir) = in

          nvertp(in) = 4

          rvertp(1,in) = rs(ik,ir) + deltar
          rvertp(2,in) = rs(ik,ir) - deltar
          rvertp(3,in) = rs(ik,ir) - deltar
          rvertp(4,in) = rs(ik,ir) + deltar

          zvertp(1,in) = zs(ik,ir) - deltaz
          zvertp(2,in) = zs(ik,ir) - deltaz
          zvertp(3,in) = zs(ik,ir) + deltaz
          zvertp(4,in) = zs(ik,ir) + deltaz

        ENDDO
      ENDDO

c...  Blank data outside the grid (an estimate for grid a39):
c      WRITE(0,*) 'Blanking camera data outside the grid'
c      DO ir = 1, nrs
c        DO ik = 1, nks(ir)
c          zlimit = -11.0 * (rs(ik,ir) - 0.525)**2.0 - 0.440
c          IF (zs(ik,ir).LT.zlimit) cdata(ik,ir) = 0.0
c        ENDDO
c      ENDDO           


      DO ir = 1, nrs
        DO ik = 1, nks(ir)

          in = korpg(ik,ir)

          WRITE(6,'(2I3,2F10.4,1X,2I5,1X,1P,E15.7,0P)')
     .       ik,ir,rs(ik,ir),zs(ik,ir),in,nvertp(in),cdata(ik,ir)

        ENDDO
      ENDDO


      DO ir = 1, nrs
        DO ik = 1, nks(ir)

          in = korpg(ik,ir)

          WRITE(6,'(2I3,4(2F8.3,1X))')
     .       ik,ir,(rvertp(i1,in),zvertp(i1,in),i1=1,4)

        ENDDO
      ENDDO




      RETURN
97    STOP 'END OF FILE ERROR'
98    STOP 'FILE ACCESS ERROR'
99    STOP
      END


