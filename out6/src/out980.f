c     -*-Fortran-*-
      subroutine fluidprobe(lvals,louts,osmvals,osmplots,
     >                 r1p,z1p,r2p,z2p,int_type,crmb,qtim)
      implicit none
c
      include 'params'
      include 'cgeom'
      include 'cedge2d'
      include 'slcom'
c
      integer osmvals,osmplots,int_type,i1
      real lvals(maxseg,maxngs),louts(maxseg)
      real r1p,z1p,r2p,z2p,crmb,qtim
c
c     OSMPROBE: This routine extract the OSM results from the background
c               plasma along a specific line intended to model the location
c               of a probe - this line is defined from P(R1,Z1) to P(R2,Z2).
c
c               The INT_TYPE parameter will allow for future enhancement for
c               calculation of cuts through any part of the grid - for
c               now - intersections are only calculated for the main SOL
c               and only for the first intersection between the line
c               segment and the rings composing the grid.
c
c               **NOTE**: INT_TYPE now defines the type of axis to be used
c                         for the rcp/osm plot
c                         1 = mid-plane distance
c                         2 = R-coordinate of intersections
c                         3 = Z-coordinate of intersections
c
c
c               The probe values returned are found by interpolating
c               between the nearest grid points. This code assumes that
c               a properly defined polygonal grid is in use - it can
c               be modified to work for a cell center only grid in necessary
c               but since these are not in common use at this time it
c               is not a useful enhancement.
c
c               The extracted values are recorded relative to the outer
c               outer mid-plane coordinates.
C
C               osmvals - output index in lvals
C                     1 - Ne
c                     2 - Te
C                     3 - Ti
c                     4 - pressure
c                     5 - Vb
c
c
c     Local Variables
c
      integer in,ip,ik,ir,icnt
      real sint,psin
      real tmp,rsect,zsect
c
c     Loop through the main SOL rings starting at the separatrix and the
c     outer target - find the first intersection with the specified LOS
c     and extract/interpolate the values of Ne, Te, Ti and pressure.
c
c     Counter for number of intersections found
c
      icnt = 0

c...  Move to out.o6a:
      outmode = 1     
      CALL WN('FluidProbe','Fluid velocity is taken to be at the '//
     .                     'cell center')

c
      do ir = irsep,irwall-1
c
         do ik = 1,nks(ir)
c
            call find_intsect(ik,ir,r1p,z1p,r2p,z2p,rsect,zsect,
     >                        sint,psin)
c
            if (sint.gt.0.0) then
c
c              Increment counter
c
               icnt = icnt + 1
c
c              Assign axis value depending on option
c
               if (int_type.eq.1) then

                  louts(icnt) = middist(ir,2)

               elseif (int_type.eq.2) then

                  louts(icnt) = rsect

               elseif (int_type.eq.3) then

                  louts(icnt) = zsect

               endif
c
c              Intersection greater than cell center
c
               if (sint.ge.kss(ik,ir)) then
c
c                 Grid point at target
c
                  if (ik.eq.nks(ir)) then

                    CALL ER('FluidProbe','Intersection between the '//
     .                      'target and the last cell center - '//
     .                      'development needed.',*99)
c
c                 Grid points away from target
c
                  else
c
c                    Ne
c
                     lvals(icnt,1) = e2dnbs(ik,ir)     +
     >                 (e2dnbs(ik+1,ir)-e2dnbs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Te
c

                     lvals(icnt,2) = e2dtebs(ik,ir)     +
     >                 (e2dtebs(ik+1,ir)-e2dtebs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Ti
c

                     lvals(icnt,3) = e2dtibs(ik,ir)     +
     >                 (e2dtibs(ik+1,ir)-e2dtibs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (e2dvhs(ik,ir)     +
     >                 (e2dvhs(ik+1,ir)-e2dvhs(ik,ir)) *
     >                 (sint-kss(ik,ir))/(kss(ik+1,ir)-kss(ik,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
                  endif
c
c              Intersection less than cell center
c

               else
c
c                 Grid points at the target
c
                  if (ik.eq.1) then

                    CALL ER('FluidProbe','Intersection between the '//
     .                      'target and the first cell center - '//
     .                      'development needed.',*99)
c
c                 Grid points away from the target
c
                  else

c
c                    Ne
c
                     lvals(icnt,1) = e2dnbs(ik,ir)     +
     >                 (e2dnbs(ik-1,ir)-e2dnbs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Te
c

                     lvals(icnt,2) = e2dtebs(ik,ir)     +
     >                 (e2dtebs(ik-1,ir)-e2dtebs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Ti
c

                     lvals(icnt,3) = e2dtibs(ik,ir)     +
     >                 (e2dtibs(ik-1,ir)-e2dtibs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir))
c
c                    Vb
c
                     lvals(icnt,5) = (e2dvhs(ik,ir)     +
     >                 (e2dvhs(ik-1,ir)-e2dvhs(ik,ir)) *
     >                 (kss(ik,ir)-sint)/(kss(ik,ir)-kss(ik-1,ir)))
     >                 /qtim
c
c                    Pressure
c
                     lvals(icnt,4) = lvals(icnt,1) * (
     >                         (lvals(icnt,2)+lvals(icnt,3)) *ech +
     >                         crmb * amu * lvals(icnt,5)**2)
c
                  endif

               endif
c slmod begin - temp
c         WRITE(6,'(5X,A,3I4,2F8.3,1P,E10.2,0P,2F8.2,1P,2E10.2,0P)') 
c     .     '--> ',
c     .     ik,ir,icnt,louts(icnt),zsect,(lvals(icnt,i1),i1=1,5)
c slmod end
               goto 100
c
            endif

         end do
c
c
c        Proceed with next ring
c
100      continue

c
      end do
c
c     Wrap up processing
c
c     Set number of points on each plot
c
      osmvals = icnt
c
c     Set number of sets of data available
c
      osmplots= 5
c
c     Verify order of scaling on axis and reorder if necessary
c
      if (louts(1).gt.louts(osmvals) ) then
c
c        Reorder from lowest to highest
c
         do in = 1,osmvals
c
            tmp = louts(osmvals-in+1)
            louts(osmvals-in+1) = louts(in)
            louts(in) = tmp
c
            do ip = 1,osmplots
c
               tmp = lvals(osmvals-in+1,ip)
               lvals(osmvals-in+1,ip) = lvals(in,ip)
               lvals(in,ip) = tmp
c
            end do
c
         end do
c
      end if
c
c
c     Print results
c
      write (6,*) 'FLUID CODE probe results:',osmvals
c
      do in = 1,osmvals
c
         write(6,'(a,i4,6(1x,g15.6))') 'FC :',in,louts(in),
     >        lvals(in,1),lvals(in,2),lvals(in,3),lvals(in,4),
     >        lvals(in,5)
c
      end do
c
      return
99    STOP
      end
c
c
c
c
c
c
c
c
      subroutine load_expt(iseld,touts,tvals,numthe,MAXTHE,MAXNGS,
     .                     MAXDATX,datatitle,colselect)
      implicit none
c
      integer iseld,numthe,MAXNGS,MAXTHE,MAXDATX,colselect
      real    touts(maxthe),tvals(maxthe)
      character datatitle*(*)
c
c     Local variables
c
c     jdemod - MAXDATX now in common block so as to standardize the 
c              interface to the experimental data - at present it was 
c              set to different values in a number of separate routines. 
c
      integer dataunit,maxcols
      parameter (dataunit=13,maxcols=5)
c
c      integer maxdatx,dataunit,maxcols
c      parameter (maxdatx=1000,dataunit=13,maxcols=1)
c
c     jdemod 
c

      integer axis_type,in,cur_index,num_expt,ipos,ncols,i1
      external ipos

c
c
c     Experimental data
c
      real expt_axis(maxdatx),expt_data(maxdatx,maxcols),theta
c
c     Write out input -
c
      write(6,*) 'Input:',iseld,numthe,maxngs
c
      call load_expt_data(dataunit,iseld,expt_axis,axis_type,expt_data,
     >                    maxcols,maxdatx,num_expt,ncols,
     >                    datatitle)
c
      if (num_expt.le.0) then
c
         write(6,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
         write(0,*) 'ERROR IN EXPERIMENTAL DATA: CAN NOT'//
     >                    ' LOAD DATASET # ',
     >                    iseld, ' - NO ELEMENTS FOUND'
c
         return
c
      endif
c
c
c
      DO i1 = 1, num_expt
        touts(i1) = expt_axis(i1)
      ENDDO

      IF     (colselect.EQ.-1) THEN
c...    Calculate upstream pressure (zero velocity assumed):
        DO i1 = 1, num_expt
          tvals(i1) = 2.0 * expt_data(i1,2) * expt_data(i1,1)*1.602E-19
        ENDDO
      ELSEIF (colselect.EQ.-2) THEN
c...    Calculate target pressure (M = 1 assumed):
        DO i1 = 1, num_expt
          tvals(i1) = 2.0 * 2.0 * expt_data(i1,2) * expt_data(i1,1) *
     .                1.602E-19
        ENDDO
      ELSEIF (colselect.EQ.-3) THEN
c...    
c       tiratio = 2.0
        DO i1 = 1, num_expt
          tvals(i1) = 2.0 * expt_data(i1,2) * expt_data(i1,1)*1.602E-19
        ENDDO

      ELSE
        DO i1 = 1, num_expt
          tvals(i1) = expt_data(i1,colselect)
        ENDDO
      ENDIF

      numthe = num_expt

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: LoadData
c
      SUBROUTINE LoadData(fname,cmnd,step,mode)
      IMPLICIT   none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER   NIZS,ITER,NITERS
      REAL      FACTA(-1:MAXIZS),FACTB(-1:MAXIZS)
      CHARACTER TITLE*80,JOB*72,EQUIL*60,graph*128

      INTEGER       step,mode
      CHARACTER*(*) fname,cmnd
      character*256 extn,resdir
      CHARACTER*128 blank

      INTEGER retcode
      INTEGER idum1,idum2
      LOGICAL mainraw

      DATA mainraw /.FALSE./

      CHARACTER command*256

      integer len,lenstr
      external lenstr
c
      integer savestep
      data savestep /-1/

c
c      len=lenstr(fname)
c      write(0,'(a,a)') 'FNAME:',fname(1:len)
c      len=lenstr(cmnd)
c      write(0,'(a,a)') 'CMND :',cmnd(1:len)
c      len=lenstr(fname)
c      write(0,'(a,4i5)') 'S/M  :',step,mode,loadstep,savestep
c

      WRITE(blank,'(128X)')

      extn = '.zip'

      CALL GetEnv('DIVRESDIR',resdir)
c...ERROR CONDITION

      IF (mode.EQ.-2) THEN
c...    Restore original fort.8:
        IF (mainraw) THEN
          CLOSE(8)

          WRITE(0,*) 'MOVING FILES'

          CALL CIssue('mvc temp.8  fort.8 ',retcode)
          CALL CIssue('mvc temp.94 plasma.dat',retcode)
          CALL CIssue('mvc temp.89 source.dat',retcode)
          CALL CIssue('mvc temp.95 geomty.dat',retcode)
          CALL CIssue('mvc temp.79 fort.79',retcode)

          WRITE(0,*) 'LOADING'

          CALL GET(TITLE,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)

          WRITE(0,*) 'LOADING SUP.'

c...      Load step if required:
          IF (savestep.NE.-1) THEN
            WRITE(graph,'(A,I2)') '      999 ',savestep    
            WRITE(6,*) 'GRAPH-2: "'//graph(1:LEN_TRIM(graph))//'"'
            CALL SetupSourcePlot(idum1,graph,idum2)
            WRITE(6,*) 'IDUM2-2: ',idum2
          ENDIF

c...      For JET grids:
          IF (REFCT.EQ.1) CALL REFLECT
        ELSE
          CALL ER('LoadData','Main .RAW file unavailable',*99)
        ENDIF
        mainraw = .FALSE.

      ELSEIF (mode.EQ.-1) THEN

        savestep = loadstep

        CLOSE(8)
        CLOSE(94)
        CLOSE(89)
        CLOSE(95)
        CLOSE(79)

        WRITE(6,*) 'ATTEMPTING TO OPEN '//fname(1:LEN_TRIM(fname))//
     .             '.raw'

c...    If this is the first call, move the existing fort.8 to temp.8:
        IF (.NOT.mainraw) THEN
          CALL CIssue('mvc fort.8  temp.8  ',retcode)
          CALL CIssue('mvc plasma.dat temp.94 ',retcode)
          CALL CIssue('mvc source.dat temp.89 ',retcode)
          CALL CIssue('mvc geomty.dat temp.95 ',retcode)
          CALL CIssue('mvc fort.79 temp.79 ',retcode)
          mainraw = .TRUE.
        ENDIF

c...    Copy file from results directory:
        command = cmnd  (1:LEN_TRIM(cmnd  ))//' '//
     .            resdir(1:LEN_TRIM(resdir))//' '//
     .            fname (1:LEN_TRIM(fname ))//blank
        CALL CISSUE(command,retcode)

c...    Load solution:
        CALL GET(TITLE,NIZS,JOB,EQUIL,FACTA,FACTB,ITER,NITERS)

c...    Load step if required:
        IF (step.NE.-1) THEN
          WRITE(graph,'(A,I2)') '      999 ',step    
          WRITE(6,*) 'GRAPH: "'//graph(1:LEN_TRIM(graph))//'"'
          CALL SetupSourcePlot(idum1,graph,idum2)
          WRITE(6,*) 'IDUM2: ',idum2
        ENDIF

c        STOP 'sdfsd'

c...    For JET grids:
        IF (REFCT.EQ.1) CALL REFLECT
      ELSE
        WRITE(0,*) 'MODE=',mode
        CALL ER('LoadData','Invalid mode',*99)
      ENDIF

      RETURN
98    CALL ER('LoadData','Unable to open '//fname(1:LEN_TRIM(fname)),
     .        *99)
99    CONTINUE
      STOP
      RETURN
      END
c
c ======================================================================
c


      SUBROUTINE Plot980(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'comgra'
      INCLUDE 'pindata'
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
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*512,graph6*128
c
c      real mvals(maxnks,maxplts,maxngs)
c      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
c      character*36 pltlabs(maxplts)
c
      CHARACTER*36 XLAB
c    >            ,XPOINT
      CHARACTER*72 YLAB
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE
c    >            ,ZLABS(-2:MAXIZS+1)
c      CHARACTER*36 NAME,PLABS(-2:MAXPLRP),KLAB
c
      CHARACTER*128 elabs(MAXNGS)

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
c
c      integer  sctype,ngrm
c      real pltmins(maxplts),pltmaxs(maxplts)
c      real pltfact
c
c      REAL tauii,pii
c
c      INTEGER nenum,tenum,opt_const,plot_mode(30),iter1,iter2,
c     .        xaxis,ring,mode,inorm(MAXNGS)
c
c
c      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
c     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
c     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
c
c      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
c     .        max1,max2,ynorm(MAXNGS)
c
c
      REAL cdata(MAXNKS,MAXNRS)
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,SUM(10),VAL,VALTOT
      CHARACTER*72 SMOOTH
      integer icntr
      integer nconts,nclev
      real conts(maxpts),clev(maxpts)
c
c      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
c      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)
c
c
c      CHARACTER*128 dataline,cdum1,cdum2
c
      CHARACTER*128 cdum1,cdum2

c...630:
c      INTEGER i,k
c      REAL    PLTMAX,PLTMIN
c      REAL LOUTS(MAXSEG),LWIDS(MAXSEG),LVALS(MAXSEG,MAXNGS)
c      REAL ydata(MAXSEG)
c      integer llabs(maxseg)
c
c...980:
      INTEGER NUMTHE,AVPTS,ATYPE,step,iflag,ln1,ln2
      INTEGER numth2,numthe1(MAXNGS),idum1,idum2,i3,i4,i5
      INTEGER IGNORS(MAXNGS),ITEC,NAVS,oldline,optscale
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS),
     .     touts1(MAXTHE,MAXNGS),maxwght0
      REAL TOUTS2(MAXTHE),TWIDS2(MAXTHE),TVALS2(MAXTHE,MAXNGS),
     .     dum1(MAXTHE),dum2(MAXTHE),dum3(MAXTHE),dum4(MAXTHE),
     .     dum5(MAXTHE),dum6(MAXTHE),dum7(MAXTHE),dum8(MAXTHE),
     .     den1,teb1,frac,rdum1,rdum2,osmtmp(MAXNKS,MAXNRS),
     .     osmhld(MAXNKS,MAXNRS)
      REAL ZOBS,ROBS,DRAD,DTHE,THEMIN,THEMAX,themin_start
      real theres,mfact,scale,qmin,qmax
      REAL WLNGTH
      real    zsuma(maxizs),cvalsa(maxnks,maxnrs)
      REAL LEVEL,AVS(0:100),VMIN,VMAX
      CHARACTER ADASID*80,PECTITLE*120,PLABAD*36
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      CHARACTER ADASGR*8,ADASTY*80,units*20,fname*1024
C
C     Second sets of ADAS data for RATIO plots
C
      character graph5*80,adasid2*80,adasex2*3
      integer   adasyr2,isele2,iselr2,iselx2,iseld2
      integer   iz_state,z_atom,iz_state2,z_atom2
      INTEGER plot,cion2
      character graph2*80,graph3*80,graph4*80
      INTEGER ISELE,ISELR,ISELX,iseld,iseldef
      INTEGER IADAS,NPAIRS,IRCODE,IKK,LEN,LENSTR
      INTEGER IK,II,IT,LT,UT,IREF,IYMIN,IYMAX,IR,JD,cnt
      INTEGER IZMIN,IZMAX,IW,LW,UW
      REAL PLRPAD(MAXNKS,MAXNRS)
      LOGICAL plotr
      real zadj
      REAL peak1,peak2,array(MAXNKS,MAXNRS),delta1,delta2,
     .     barerror1,barerror2
      INTEGER plotcode,line,iopt1,iopt2,iopt3
      CHARACTER cname*256,resdir*256,cmnd*256
      LOGICAL status,oldraw
      INTEGER barmode,barseries,barindex


      COMMON /NSMOOTH/ NUMSMOOTH,cgrprint
      INTEGER NUMSMOOTH,cgrprint

c...  For reading from the .experiment file (UNIT=13):
      INTEGER    MAXTDAT     ,MAXCOLS   
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)
      INTEGER   etype,ndata,ncol
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS)
      CHARACTER datatitle*128

c...  Mapping data from PSIn to azimuthal coordinates:
      REAL ATAN3C
      INTEGER ndat,subopt
      LOGICAL success
      REAL ang(MAXNDS),psi(MAXNDS),r,z
c
c     Array for extra comments 
c
      character*20 extra_comments(1)
c
c jdemod - note this is a compiled in option for now - mostly hacked
c
c     Variables for mapping to PSIN for cumulative plots 
c
      real psin
      integer targ 
      logical plot_psin  
c
c jdemod   
c
      nview = ' '
      plane = ' '
      anly  = ' '
      TABLE = 'SYMBOL TABLE'

      ngs = 0


      iflag = 1
      IF (ismoth.EQ.-1) THEN
        ismoth = 99
        iflag = 1
      ENDIF

      READ(graph(14:15),*) plot
      WRITE(6,*) 'OUT980: PLOT=',plot
      WRITE(0,*) 'OUT980: PLOT=',plot


c...temp!
      CALL RSet (tvals  ,MAXTHE*MAXNGS,LO)
      CALL RSet (tvals2 ,MAXTHE*MAXNGS,LO)
      CALL IZero(numthe1,MAXNGS)
      CALL RZero(touts2 ,MAXTHE)
      CALL RZero(touts1 ,MAXTHE*MAXNGS)
      CALL RZero(array  ,MAXNKS*MAXNRS)

      xlab = 'degrees'

c...  Read ADAS data specifications from input file:
      IF (.NOT.(plot.GE.4.AND.plot.LE.15).AND.plot.NE.24) THEN
        WRITE(6,*) 'READING LINE 1'

c...    Making ADAS line optional, since I have been sloppy about
c       what graphs actually need this data, and I don't want to
c       have to clean up all my old OUT files:
        READ(5,'(A128)',END=6) graph6
        IF (graph6(8:11).EQ.'Adas'.OR.graph6(8:11).EQ.'ADAS'.OR.
     .      graph6(8:11).EQ.'adas') THEN
          BACKSPACE 5

          CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .               ISELD,IERR)
          IF (IERR.NE.0) THEN
            WRITE(6,*) 'ERROR READING ADAS DETAILS, IERR = ',IERR
            IERR = 0
            GOTO 99
          ENDIF

        ELSE
          BACKSPACE 5
        ENDIF
 6      CONTINUE

      ENDIF

c...  Read 2nd line of ADAS data for line ratio plot:
      IF (plot.EQ.3) THEN
        WRITE(6,*) 'READING LINE 2'
        CALL RDG1(graph3,adasid2,adasyr2,adasex2,
     >            isele2,iselr2,iselx2,iseld2,ierr)
        IF (IERR.NE.0) THEN
          WRITE(6,*) 'ERROR READING ADAS DETAILS, IERR = ',IERR
          IERR = 0
          GOTO 99
        ENDIF
      ENDIF

c...  Read plots specifications:
      CALL RDG2 (GRAPH2,ROBS,ZOBS,THEMIN,DTHE,THERES,NUMTHE,
     >           IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
      IF (IERR.EQ.1) THEN
         WRITE(6,*) 'RDG2 ERROR READING 980 SERIES- GRAPH '//
     .              'DETAILS'
         WRITE(6,*) graph2
         IERR = 0
         GOTO 99
      ENDIF

c...  Shift viewing vertex if necessary (this blows at the moment
c     because it is based on an simple z-shift, which is only
c     'good' for DIII-D since the target plate is horizontal -- and
c     even then, it is a poor approximation farther out on the targets
c     where they are not.  Also, the inner and outer shift could be
c     different in theory, which isn't good):
      IF (tarshift(IKHI).NE.0.0) THEN
        CALL WN('980','R-shift applied to viewing vertex')
        robs = robs + tarshift(IKHI)        
        WRITE(6,*) 'ROBS= ',robs
      ENDIF

c...  Quick and dirty method of making R or theta or PSIN plots optional:
c
c     jdemod - use also for plots vs. PSIN
c            - if numsmooth <= -10 then set EXTERNAL PSIN axis instead
c
      plotr = .false.
      plot_psin = .false.
      if (numsmooth.lt.0.and.numsmooth.gt.-10) then
         plotr = .true.
         numsmooth = abs(numsmooth)
         numsmooth = 0
         WRITE(6,*) '980:   SETTING NUMSMOOTH TO 0 -'//
     >              ' R AXIS SELECTED'
      elseif (numsmooth.le.-10) then
         plot_psin = .true.
         numsmooth = 0
         WRITE(6,*) '980:   SETTING NUMSMOOTH TO 0 -'//
     >              ' PSIN AXIS SELECTED'
      endif
c
c     jdemod
c

c...For C-Mod:
c     if (cgridopt.eq.3) themin = themin + 180.0

      slopt = 1


c...  Create THETA vector:
      idum2  = numthe
      numthe = 0
5     READ(5,'(A512)') graph1
      IF (graph1(8:11).EQ.'View'.OR.graph1(8:11).EQ.'VIEW'.OR.
     .    graph1(8:11).EQ.'view') THEN
        READ(graph1,*) cdum1,idum1,(touts(ii),ii=1+numthe,idum1+numthe)
        numthe = numthe + idum1
        THEMIN = TOUTS(1)
        THEMAX = TOUTS(NUMTHE)
        GOTO 5
      ELSE
        IF (numthe.EQ.0) THEN
          numthe = idum2
          DO II = 1, NUMTHE
            TOUTS(II) = THEMIN + (II-1)*DTHE
          ENDDO
          THEMAX = TOUTS(NUMTHE)
        ENDIF
        BACKSPACE 5
      ENDIF

c...  Specify width of viewing cone (this rather convoluted method
c     of reading in the data allows more than one 'CONE' data line
c     to be specified in the OUT input file):
      idum1 = 0
      idum2 = 0
7     READ(5,'(A512)') graph1
      IF (graph1(8:11).EQ.'Cone'.OR.graph1(8:11).EQ.'CONE'.OR.
     .    graph1(8:11).EQ.'cone') THEN
        READ(graph1,*) cdum1,idum1,(twids(ii),ii=1+idum2,idum1+idum2)
        IF (idum1.EQ.1.AND.idum2.EQ.0) THEN
c...      Assume that a single value is provided that is to be applied
c         to all viewing cones:
          DO ii = 2, numthe
            twids(ii) = twids(1)
          ENDDO
        ENDIF
        THERES = -99.0
        idum2 = idum2 + idum1
        GOTO 7
c...     Old code.  I can't figure why it was done this way. - SL, Feb, 2004
c        READ(graph1,*) cdum1,idum1,(twids(ii),ii=1+numthe,idum1+numthe)
c        numthe = numthe + idum1
c        THERES = -99.0
c        GOTO 7
      ELSE
        IF (idum1.EQ.0) THEN
c...      Individual cone widths were not found, so set TWIDS 
c         from THERES:
          DO II = 1, NUMTHE
            TWIDS(II) = THERES
          ENDDO
        ENDIF
        BACKSPACE 5
      ENDIF


c...  Set scale factor:
      units = '                    '
      IF (ATYPE.EQ.0) THEN
        MFACT = 1.0
        WRITE(char(1),'(''NO SCALE FACTOR APPLIED'')')
        units = '(PH M-2 S-1)'
      ELSEIF (ATYPE.EQ.1) THEN
        MFACT = DTHE / (2.0 * PI)
        WRITE(char(1),'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >             DTHE
      ELSEIF (ATYPE.EQ.2) THEN
        MFACT = 1.0 / (2.0 * PI)
        WRITE(char(1),'(''SCALE FACTOR = 1 / (2*PI)'') ')
      ELSEIF (ATYPE.EQ.3) THEN
        MFACT = 1.0 / ( 4.0 * PI)
        WRITE(char(1),'(''SCALE FACTOR = 1 / (4*PI)'')')
        units = '(PH M-2 S-1 SR-1)'
c        YLAB = 'PLRP (PH M-2 S-1 SR-1)'
      ENDIF
      WRITE(0,*) '980:',mfact,units

c...  Decide what to plot - PINCODE is passed from DIVIMP 
c     and is the parameter specifiying which neutral code
c     was run:
c
c     0 - NIMBUS
c     1 - EIRENE97
c     2 - EIRENE99
c
      IF (pincode.EQ.-1.OR.pincode.EQ.5) THEN
        IF (graph(5:5).EQ.'0'.OR.graph(5:5).EQ.'1'.OR.
     .      graph(5:5).EQ.'2'.OR.graph(5:5).EQ.'9') THEN
          READ(graph(5:5),*) plotcode
        ELSE
          plotcode = -1
        ENDIF
      ELSE
        plotcode = pincode
      ENDIF 

      WRITE(0,*) 'PLOTCODE:',plotcode
c
c
c     PLOTS
c
c
      oldraw = .FALSE.
      oldline = 0

      IF (plot.EQ.1.OR.plot.EQ.2.OR.plot.EQ.23) THEN

c...    Turn on colour plot lines:
        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 
        ENDDO
        plottype(1) = 55

        IF (graph(5:5).EQ.'3') slopt5 = 1
c...    Glorious hack!  Needed for Kbottom array on C-Mod, which is located in the plenum 
c       chamber.  This confuses the LOS plot code that counts wall intersections, since it assumes
c       that the array connot see outside the chamber, which is not the case of course:
        IF (graph(5:5).EQ.'4') slopt5 = 2

        ref      = graph(7:46)
        IF     (plot.EQ.1) THEN
c...      Dalpha:
          YLAB = 'Dalpha '//units
          line = H_BALPHA
        ELSEIF (plot.EQ.2) THEN
c...      Dgamma:
          YLAB = 'Dgamma '//units
          line = H_BGAMMA
        ELSEIF (plot.EQ.23) THEN
c...      Dbeta:
          YLAB = 'Dbeta '//units
          line = H_BBETA
        ENDIF

c...    Store the x-axis data:
        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

        elabs(1) = '                                                  '
        elabs(2) = '                                                  '

c...    Load experimental data:
        if (machine.EQ.DIIID.AND.iseld.ne.0) then
          ngs = ngs + 1
          numthe1(ngs) = 0
          CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                   MAXTHE,MAXNGS,MAXDATX,datatitle,1)
          elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))

          IF (.NOT..TRUE.) THEN
c...        FOR OUTER TARGET, MAP PSIn to ANGLE:
            WRITE(0,*) '980: PSIn to ANGLE mapping',plot,iseld

            DO i1 = 1, numthe1(ngs)
              WRITE(6,*) 'DATA:',touts1(i1,ngs),tvals(i1,ngs)
            ENDDO
            i1 = 0
            DO ir = irtrap+1, nrs
              r = rp(idds(ir,1))
              z = zp(idds(ir,1))
              IF (z.NE.zp(idds(irsep,1))) CYCLE
              i1 = i1 + 1
              psi(i1) = psitarg(ir,1)
              ang(i1) = ATAN3C(z-1.533,r-2.027)
c              WRITE(0,*) 'STUFF:',ir,z,psi(i1),ang(i1)
            ENDDO
            DO ir = irsep,irwall-1
              r = rp(idds(ir,1))
              z = zp(idds(ir,1))
              IF (z.NE.zp(idds(irsep,1))) CYCLE
              i1 = i1 + 1
              psi(i1) = psitarg(ir,1)
              ang(i1) = ATAN3C(z-1.533,r-2.027)
c              WRITE(0,*) 'STUFF:',ir,psi(i1),ang(i1)
            ENDDO      
            ndat = i1
c...        Interpolate x-axis data:
            i1 = 1
            DO WHILE (i1.LE.numthe1(ngs))
c              WRITE(0,*) 'PROGRESS:',i1,numthe1(ngs)
              success = .FALSE.
              DO i2 = 1, ndat-1
                IF (touts1(i1,ngs).GT.psi(i2).AND.
     .              touts1(i1,ngs).LT.psi(i2+1)) THEN
                  frac = (touts1(i1,ngs) - psi(i2)) / 
     .                   (psi(i2+1) - psi(i2))
                  touts1(i1,ngs) = ang(i2) + frac * (ang(i2+1)-ang(i2))
                  success = .TRUE.                  
                  EXIT
                ENDIF
              ENDDO
              IF (.NOT.success) THEN
c...            Drop data point:
                DO i2 = i1, numthe1(ngs)
                  touts1(i2,ngs) = touts1(i2+1,ngs)
                  tvals (i2,ngs) = tvals (i2+1,ngs)
                ENDDO
                numthe1(ngs) = numthe1(ngs) - 1
c                WRITE(0,*) 'CUSTOM HACK: INTERPOLATION FAILED, '//
c     .                     'HALTING OUT',i1,touts1(i1,ngs)
c                STOP
              ELSE
                i1 = i1 + 1
              ENDIF
            ENDDO
            DO i1 = 1, numthe1(ngs)
              WRITE(6,*) 'DATA 2:',touts1(i1,ngs),tvals(i1,ngs)
            ENDDO
          ENDIF

          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
c...      Plot black squares:
          plottype(ngs) = -55




c...SPECIAL
          IF (.NOT..TRUE..AND.plot.EQ.1) THEN
            ngs = ngs + 1
            numthe1(ngs) = 0
            CALL load_expt(82,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                     MAXTHE,MAXNGS,MAXDATX,datatitle,1)
            elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))

c...        FOR OUTER TARGET, MAP PSIn to ANGLE:
            DO i1 = 1, numthe1(ngs)
              WRITE(6,*) 'DATA:',touts1(i1,ngs),tvals(i1,ngs)
            ENDDO
            i1 = 0
            DO ir = irtrap+1, nrs
              r = rp(idds(ir,1))
              z = zp(idds(ir,1))
              IF (z.NE.zp(idds(irsep,1))) CYCLE
              i1 = i1 + 1
              psi(i1) = psitarg(ir,1)
              ang(i1) = ATAN3C(z-1.533,r-2.027)
c              WRITE(0,*) 'STUFF:',ir,z,psi(i1),ang(i1)
            ENDDO
            DO ir = irsep,irwall-1
              r = rp(idds(ir,1))
              z = zp(idds(ir,1))
              IF (z.NE.zp(idds(irsep,1))) CYCLE
              i1 = i1 + 1
              psi(i1) = psitarg(ir,1)
              ang(i1) = ATAN3C(z-1.533,r-2.027)
c              WRITE(0,*) 'STUFF:',ir,psi(i1),ang(i1)
            ENDDO      
            ndat = i1
c...        Interpolate x-axis data:
            i1 = 1
            DO WHILE (i1.LE.numthe1(ngs))
c              WRITE(0,*) 'PROGRESS:',i1,numthe1(ngs)
              success = .FALSE.
              DO i2 = 1, ndat-1
                IF (touts1(i1,ngs).GT.psi(i2).AND.
     .              touts1(i1,ngs).LT.psi(i2+1)) THEN
                  frac = (touts1(i1,ngs) - psi(i2)) / 
     .                   (psi(i2+1) - psi(i2))
                  touts1(i1,ngs) = ang(i2) + frac * (ang(i2+1)-ang(i2))
                  success = .TRUE.                  
                  EXIT
                ENDIF
              ENDDO
              IF (.NOT.success) THEN
c...            Drop data point:
                DO i2 = i1, numthe1(ngs)
                  touts1(i2,ngs) = touts1(i2+1,ngs)
                  tvals (i2,ngs) = tvals (i2+1,ngs)
                ENDDO
                numthe1(ngs) = numthe1(ngs) - 1
c                WRITE(0,*) 'CUSTOM HACK: INTERPOLATION FAILED, '//
c     .                     'HALTING OUT',i1,touts1(i1,ngs)
c                STOP
              ELSE
                i1 = i1 + 1
              ENDIF
            ENDDO
            DO i1 = 1, numthe1(ngs)
              WRITE(6,*) 'DATA 2:',touts1(i1,ngs),tvals(i1,ngs)
            ENDDO

            CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
c...        Plot black squares:
            plottype(ngs) = -2
          ENDIF
c...      END OF DIII-D BUSINESS
        ENDIF

c...    A loop is required in case results are plotted for more than
c       one case:
        status = .TRUE.
        DO WHILE (status)

c...      Conduct the LOS integrations for the ADAS data:
          IF ((plotcode.EQ.-1.OR.(plotcode.EQ.0.AND.plot.EQ.1).OR.
     .                           (plotcode.EQ.0.AND.plot.EQ.2).OR.
     .                           (plotcode.EQ.9.AND.plot.EQ.1)).OR.
     .        .FALSE.) THEN
            ngs = ngs + 1
            WRITE(0,*) 'PLOT 980: STARTING ADAS',ngs
            IF     (plotcode.EQ.-1) THEN
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    ADAS'
            ELSEIF (plotcode.EQ.0) THEN
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    NIMBUS-ADAS'
            ELSE
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    EIRENE-ADAS'
            ENDIF

            CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
  
            CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)

            WRITE(0,*) 'CVALSA:',1,irsep,cvalsa(1,irsep)

c            WRITE(0,*) '980: USING slLOSINT'
            CALL slLOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,
     .                    AVPTS,CVALSA,0.0)

        

c...        Store y-axis data for this data set:
            numthe1(ngs) = numthe
            DO i1 = 1, numthe
              touts1(i1,ngs) = touts(i1)
            ENDDO

c            WRITE(0,*) 'PLOT 980: DONE'
          ENDIF

c...      Conduct the LOS integrations for the EIRENE data:
          IF ((plotcode.EQ.-1.OR.(plotcode.EQ.0.AND.plot.EQ.1).OR.
     .         plotcode.EQ. 1.OR. plotcode.EQ.2.OR.
     .         plotcode.EQ. 4).AND.
     .        .TRUE.) THEN
            ngs = ngs + 1
            WRITE(0,*) 'PLOT 980: STARTING EIRENE',ngs,plotcode,line
            IF     (plotcode.EQ.-1) THEN
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    OSM+PIN'
            ELSEIF (plotcode.EQ. 0) THEN
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    NIMBUS'
            ELSE
              IF (elabs(ngs)(5:5).EQ.' ') elabs(ngs) = '    EIRENE'
c              IF (elabs(1)(5:5).EQ.' ') elabs(ngs) = '    OSM+EIRENE'
            ENDIF
            IF (plotcode.EQ.0) THEN
              CALL RVALKR(CVALSA,pinalpha           ,1,1,1,FT,FP,
     >                    MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            ELSE
              CALL RVALKR(CVALSA,pinline(1,1,6,line),1,1,1,FT,FP,
     >                    MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            ENDIF

            WRITE(0,*) 'CVALSA:',1,irsep,cvalsa(1,irsep)

c            WRITE(0,*) '980: USING slLOSINT'
            CALL slLOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,
     .                    AVPTS,CVALSA,0.0)




c...        Store y-axis data for this data set:
            numthe1(ngs) = numthe
            DO i1 = 1, numthe
              touts1(i1,ngs) = touts(i1)
            ENDDO

            IF (.FALSE.) THEN
              IF (machine.EQ.CMOD.AND.iseld.EQ.3) THEN
                WRITE(0,*) '980: **************************' 
                WRITE(0,*) '980: ***SCALING MODEL DALPHA***' 
                WRITE(0,*) '980: **************************' 
                DO i1 = 1, numthe
                  tvals(i1,ngs) = tvals(i1,ngs) * 0.4
                ENDDO
              ENDIF
              IF (machine.EQ.CMOD.AND.iseld.EQ.6) THEN
                WRITE(0,*) '980: **************************' 
                WRITE(0,*) '980: ***SCALING MODEL DALPHA***' 
                WRITE(0,*) '980: **************************' 
                DO i1 = 1, numthe
                  tvals(i1,ngs) = tvals(i1,ngs) * 0.4
                ENDDO
              ENDIF  
            ENDIF


          ENDIF

c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          status = .FALSE.
          READ(5,'(A128)',END=10) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,elabs(ngs+1),cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
c...        TEMPORARY!
            IF (pincode.EQ.-1) THEN
              plotcode = iopt2
            ELSE
              plotcode = pincode
            ENDIF
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
c...      Do LOS intergration for camera data:
          READ(5,'(A128)',END=10) graph6
          IF (graph6(8:11).EQ.'File'.OR.graph6(8:11).EQ.'FILE'.OR.
     .        graph6(8:11).EQ.'file') THEN
            scale = 0.0
            WRITE(fname,'(1024X)')
            READ(graph6,*) cdum1,line,optscale,scale,fname
            osmtmp = 0.0
            CALL LoadCameraData(osmtmp,fname,scale)
            IF     (machine.EQ.DIIID) THEN
              WRITE(0,*) 'WARNING 980: HARDCODED DIII-D SCALING OF'//
     .                   ' CAMERA DATA'
            ELSEIF (machine.EQ.CMOD) THEN                
              WRITE(0,*) 'WARNING 980: HARDCODED C-MOD SCALING OF'//
     .                   ' DGAMMA DATA'
            ENDIF
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                IF (oldline.EQ.0) osmhld(ik,ir) = pinline(ik,ir,6,line)
                pinline(ik,ir,6,line) = osmtmp(ik,ir)
c *HARDCODED SCALING*
                IF     (machine.EQ.DIIID) THEN
                  IF (line.EQ.1) pinline(ik,ir,6,line) = 
     .                           pinline(ik,ir,6,line) * 1.5E+07
c     .                           pinline(ik,ir,6,line) * 1.3E+13
                  IF (line.EQ.2) pinline(ik,ir,6,line) = 
     .                           pinline(ik,ir,6,line) * 2.0E+19
c     .                           pinline(ik,ir,6,line) * 0.4*3.0E+19  ! Eric's MAR
c     .                           pinline(ik,ir,6,line) * 3.0E+12
c     .                           pinline(ik,ir,6,line) * 6.5E+12
                ELSEIF (machine.EQ.CMOD) THEN                
                  IF (line.EQ.2) pinline(ik,ir,6,line) = 
     .                           pinline(ik,ir,6,line) * 0.19E+20
                ENDIF
              ENDDO
            ENDDO
            plottype(ngs+1) = -55
            oldline = line
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
c...      Cycle again, but for reflections:
          READ(5,'(A128)',END=10) graph6
          IF (graph6(8:11).EQ.'Cycl'.OR.graph6(8:11).EQ.'CYCL'.OR.
     .        graph6(8:11).EQ.'cycl') THEN
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
10        CONTINUE

        ENDDO

c...    Restore main .raw file if necessary:
        IF (oldline.NE.0) THEN
          DO ir = 1, nrs
            DO ik = 1, nks(ir)
              pinline(ik,ir,6,oldline) = osmhld(ik,ir)
            ENDDO
          ENDDO
        ENDIF



c...    Restore main .raw file if necessary:
        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

c...    Load experimental data:
        if (machine.EQ.CMOD.AND.iseld.ne.0) then
          ngs = ngs + 1
          numthe1(ngs) = 0
          CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                   MAXTHE,MAXNGS,MAXDATX,datatitle,1)
          elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))

          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
c...      Plot black squares:
          plottype(ngs) = -55
       
          IF (iseld.EQ.3) plottype(ngs) = 2
          IF (iseld.EQ.6) plottype(ngs) = 2

          IF (.NOT..TRUE.) THEN
c...       Load Kbot straight: 
           ngs = ngs + 1
           numthe1(ngs) = 0
           CALL load_expt(5    ,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                    MAXTHE,MAXNGS,MAXDATX,datatitle,1)
           elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))
           CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
c...       Plot black squares:
           plottype(ngs) = -55
          ENDIF

        ENDIF
c...    Map individual y-axis data sets to the master y-axis data array:
        DO i1 = 1, ngs
          CALL MapArray(touts2      ,tvals2(1,i1),1,numth2,
     .                  touts1(1,i1),tvals (1,i1),1,numthe1(i1))
        ENDDO


c...    Adjust the x-axis data if required:
        if (plotr) then
          WRITE(6,*) 'PLOTR: ',zadj
          call adjustout(touts2,numth2,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numth2)
          write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
        endif

      ELSEIF (plot.EQ.3) THEN
c       
c...    Line ratios:
c
c  subopt 0: thesis
c  subopt 1: Dbeta/Dalpha for DIII-D
c

        subopt = 1
         
        ngs  = 0
        ref  = graph(7:46)
        YLAB = 'Dbeta / Dalpha'
c        YLAB = 'Dalpha / Dgamma'
        XLAB = 'Viewing angle (degrees)'

        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 
        ENDDO
        plottype(1) = 55

c        WRITE(0,*) '980:03 USING slLOSINT'

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

c...    Experiment : HARDCODED!
        IF (subopt.EQ.0) THEN
          CALL Load_Expt_Data(13,7,xdata,etype,edata,
     .                        MAXCOLS,MAXTDAT,ndata,ncol,datatitle)
          CALL LoadArray(touts2,numth2,xdata,1,ndata)

        ELSE

c...      Load experimental data:
          if (iseld.ne.0) then
            ngs = ngs + 1
            numthe1(ngs) = 0
            CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs)
     .                     ,MAXTHE,MAXNGS,MAXDATX,datatitle,1)
            elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))

            IF (.TRUE.) THEN
c... FOR OUTER TARGET, MAP PSIn to ANGLE:
              DO i1 = 1, numthe1(ngs)
                WRITE(6,*) 'DATA:',touts1(i1,ngs),tvals(i1,ngs)
              ENDDO
              i1 = 0
              DO ir = irtrap+1, nrs
                r = rp(idds(ir,1))
                z = zp(idds(ir,1))
                IF (z.NE.zp(idds(irsep,1))) CYCLE
                i1 = i1 + 1
                psi(i1) = psitarg(ir,1)
                ang(i1) = ATAN3C(z-1.533,r-2.027)
              ENDDO
              DO ir = irsep,irwall-1
                r = rp(idds(ir,1))
                z = zp(idds(ir,1))
                IF (z.NE.zp(idds(irsep,1))) CYCLE
                i1 = i1 + 1
                psi(i1) = psitarg(ir,1)
                ang(i1) = ATAN3C(z-1.533,r-2.027)
              ENDDO      
              ndat = i1
c...          Interpolate x-axis data:
              i1 = 1
              DO WHILE (i1.LE.numthe1(ngs))
                success = .FALSE.
                DO i2 = 1, ndat-1
                  IF (touts1(i1,ngs).GT.psi(i2).AND.
     .                touts1(i1,ngs).LT.psi(i2+1)) THEN
                    frac = (touts1(i1,ngs) - psi(i2)) / 
     .                     (psi(i2+1) - psi(i2))
                    touts1(i1,ngs) = ang(i2) + frac*(ang(i2+1)-ang(i2))
                    success = .TRUE.                  
                    EXIT
                  ENDIF
                ENDDO
                IF (.NOT.success) THEN
c...              Drop data point:
                  DO i2 = i1, numthe1(ngs)
                    touts1(i2,ngs) = touts1(i2+1,ngs)
                    tvals (i2,ngs) = tvals (i2+1,ngs)
                  ENDDO
                  numthe1(ngs) = numthe1(ngs) - 1
                ELSE
                  i1 = i1 + 1
                ENDIF
              ENDDO
              DO i1 = 1, numthe1(ngs)
                WRITE(6,*) 'DATA 2:',touts1(i1,ngs),tvals(i1,ngs)
              ENDDO
            ENDIF


            CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
            CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                    touts1(1,ngs),tvals(1,ngs),1,numthe1(ngs))
c...        Plot black squares:
            plottype(ngs) = -55
          ENDIF



        ENDIF

        status = .TRUE.
        oldraw = .FALSE.
        DO WHILE (status)

c...TEMP!          
c          plotcode = -1

          IF (plotcode.EQ.-1.OR.plotcode.EQ.0.OR.
     .        plotcode.EQ.999) THEN
c...        ADAS:
            ngs = ngs + 1
            IF (plotcode.EQ.-1.OR.plotcode.EQ.2) THEN
              elabs(ngs) = '    ADAS'
            ELSE
              elabs(ngs) = '    NIMBUS-ADAS'
            ENDIF

            CALL RSet(tvals,MAXTHE*MAXNGS,0.0)

            CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
            CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                    CVALSA,0.0)

            CALL LDADAS(1,IZMIN,ADASID2,ADASYR2,ADASEX2,ISELE2,ISELR2,
     >                  ISELX2,plrpad,Wlngth,IRCODE)
            CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,2),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                    CVALSA,0.0)

            DO i1 = 1, numthe
              tvals(i1,1) = tvals(i1,1) / tvals(i1,2)
            ENDDO

            CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                    touts ,tvals (1,1  ),1,numthe)
          ENDIF

          IF (plotcode.EQ.-1.OR.plotcode.EQ.1.OR.plotcode.EQ.2.OR.
     .        plotcode.EQ. 4) THEN
c...        PIN:
            ngs = ngs + 1
            IF (plotcode.EQ.-1) THEN
              elabs(ngs) = '    OSM+PIN'
            ELSE
              elabs(ngs) = '    OSM+EIRENE'
            ENDIF

            CALL RSet(tvals,MAXTHE*MAXNGS,0.0)

c            CALL RVALKR(CVALSA,pinline(1,1,2,H_BBETA),1,1,1,FT,FP,
            CALL RVALKR(CVALSA,pinline(1,1,6,H_BBETA),1,1,1,FT,FP,
c            CALL RVALKR(CVALSA,pinline(1,1,6,H_BALPHA),1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                    CVALSA,0.0)

c            CALL RVALKR(CVALSA,pinline(1,1,2,H_BALPHA),1,1,1,FT,FP,
            CALL RVALKR(CVALSA,pinline(1,1,6,H_BALPHA),1,1,1,FT,FP,
c            CALL RVALKR(CVALSA,pinline(1,1,6,H_BGAMMA),1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,2),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                    CVALSA,0.0)

            DO i1 = 1, numthe
              tvals(i1,1) = tvals(i1,1) / tvals(i1,2)
            ENDDO

            CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                    touts ,tvals (1,1  ),1,numthe)
          ENDIF

c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          status = .FALSE.
          READ(5,'(A128)',END=20) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,cname,iopt1,iopt2,cmnd
            CALL LoadData(cname,cmnd,-1,iopt1)
c...        TEMPORARY!
            IF (pincode.EQ.-1) THEN
              plotcode = iopt2
            ELSE
              plotcode = pincode
            ENDIF
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
c...      Cycle again, but for reflections:
c
c          jdemod - I think this should be end=20 
c                 - otherwise the branch is outside the loop
c
c          READ(5,'(A128)',END=10) graph6
c
          READ(5,'(A128)',END=20) graph6
          IF (graph6(8:11).EQ.'Cycl'.OR.graph6(8:11).EQ.'CYCL'.OR.
     .        graph6(8:11).EQ.'cycl') THEN
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
20        CONTINUE
        ENDDO

        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

c...    Modify independent variable:
        if (plotr) then
           call adjustout(touts2,numthe,zadj,robs,zobs)
           themin = touts2(1)
           themax = touts2(numthe)
           write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
        endif

c...    Experiment : HARDCODED!
        IF (subopt.EQ.0) THEN
          WRITE(0,*) '980:03 EXPERIMENTAL DATA HARDCODED'
          ngs = ngs + 1
          plottype(ngs) = -55
          elabs(ngs) = '    experiment'
          CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                  xdata ,edata (1,1  ),1,ndata)
        ENDIF


c        call calc_expt(7,touts2,tvals2,maxthe,numthe,
c     >                 themin,themax,maxngs,ngs,datatitle)
c        ngs = ngs + 1
c        call calc_expt(2,touts2,tvals2,maxthe,numthe,
c     >                 themin,themax,maxngs,ngs,datatitle)
c        ngs = ngs - 1
c...    Experimental data less than a significant fraction of the peak
c       is very noisy, so edit it out:
c        peak1 = -HI
c        peak2 = -HI
c        DO i1 = 1, numth2
c          peak1 = MAX(peak1,tvals2(i1,ngs  ))
c          peak2 = MAX(peak2,tvals2(i1,ngs+1))
c        ENDDO
c        DO i1 = 1, numth2
c          IF (tvals2(i1,ngs  ).LT.0.05*peak1.OR.
c     .        tvals2(i1,ngs+1).LT.0.05*peak2) THEN
c            tvals2(i1,ngs) = LO
c          ELSE
c            tvals2(i1,ngs) = tvals2(i1,ngs) / tvals2(i1,ngs+1)
c          ENDIF
c        ENDDO

      ELSEIF (plot.EQ.4.OR.plot.EQ.10.OR.plot.EQ.24) THEN
c...    Dalpha (4) or Dgamma (10) components (cumulative):
        ngs      = 5
        ref      = graph(7:46)
        elabs(1) = '    H (excitation)'
        elabs(2) = '    H+ (recombination)'
        elabs(3) = '    H2'
        elabs(4) = '    H2+'
        elabs(5) = '    H-'

        IF (graph(5:5).EQ.'3') slopt5 = 1
c...    Glorious hack!  Needed for Kbottom array on C-Mod, which is located in the plenum 
c       chamber.  This confuses the LOS plot code that counts wall intersections, since it assumes
c       that the array connot see outside the chamber, which is not the case of course:
        IF (graph(5:5).EQ.'4') slopt5 = 2

c...    Turn on coloured lines on plots:
        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 + 1
        ENDDO


        WRITE(cdum1,'(A,I2,A)') '980:',plot,' USING slLOSINT'
c        WRITE(0,*) cdum1

        IF     (plot.EQ. 4) THEN
          YLAB = 'Dalpha (cumulative) '//units
          line = H_BALPHA
        ELSEIF (plot.EQ.10) THEN
          YLAB = 'Dgamma (cumulative) '//units
          line = H_BGAMMA
        ELSEIF (plot.EQ.24) THEN
          YLAB = 'Dbeta (cumulative) '//units
          line = H_BBETA
        ENDIF

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts         ,1,numthe)

        DO i1 = 1, 5
          CALL RVALKR(CVALSA,pinline(1,1,i1,line),1,1,1,FT,FP,
     >                MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
          CALL slLOSINT(TVALS(1,i1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                CVALSA,0.0)
c          CALL LOSINT(TVALS(1,i1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
c     .                CVALSA,0.0,0)
          CALL MapArray(touts2,tvals2(1,i1),1,numth2,
     .                  touts ,tvals (1,i1),1,numthe)
        ENDDO

        DO i1 = 2, 5
          DO i2 = 1, numth2
            tvals2(i2,i1) = tvals2(i2,i1) + tvals2(i2,i1-1)
          ENDDO
        ENDDO

        if (plotr) then
          call adjustout(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
        endif
c
c jdemod - added as hard-coded (for now) option to map as psin 
c        - changed to use negative value of numsmooth which is reset to 0
c
c       Load PSIN data from file - needs input line - and then 
c       map the angle to a PSIN intersection - the PSIN load data
c       could be made optional
c
c       Get filename and load data for an externally specified PSIN mapping
c
        if (plot_psin) then 
c
            call get_psin_filename(ierr)
            if (ierr.ne.0) return 
c
c           Load PSIN reference data
c
            call read_psin_data
c
            do in = 1,numthe
c
               call calcpsin(robs,zobs,touts2(in),psin,targ)
c
               touts2(in) = psin
c
            end do 
c
            themin = touts2(1)
            themax = touts2(numthe)
            write(XLAB,'(''PSIN'')') 
            if (atype.eq.3) then 
               IF     (plot.EQ. 4) THEN
                 YLAB = 'Dalpha components (cumulative) (PH/M2/S/SR)'
               ELSEIF (plot.EQ.10) THEN
                 YLAB = 'Dgamma components (cumulative) (PH/M2/S/SR)'
               ENDIF
            endif
c
         endif
c
c jdemod
c
c-----------
c
c jdemod - added plot 25 
c
      ELSEIF (plot.EQ.25.OR.plot.EQ.26) THEN

        IF (plot.EQ.25) THEN
          cion2 = cion
        ELSE
c...      Hydrogenic data:  
          cion2 = 1
        ENDIF
c
c...    ADAS Line (25) components (cumulative):
c
        len = lenstr(graph)
        ref = graph(17:len)
c
c       Read in ADAS data input line - already read for plots outside 4 to 15. 
c
c        CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
c                   ISELD,IERR)
c
        IF (graph(5:5).EQ.'3') slopt5 = 1
c...    Glorious hack!  Needed for Kbottom array on C-Mod, which is located in the plenum 
c       chamber.  This confuses the LOS plot code that counts wall intersections, since it assumes
c       that the array connot see outside the chamber, which is not the case of course:
        IF (graph(5:5).EQ.'4') slopt5 = 2

c...    Turn on coloured lines on plots:
        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 + 1
        ENDDO


        WRITE(cdum1,'(A,I2,A)') '980:',plot,' USING slLOSINT'
c        WRITE(0,*) cdum1
c
         if (atype.eq.3) then 

            ylab = 'ADAS Spectral Line components '//
     >             '(cumulative) (PH/M2/S/SR)'
         else

            ylab = 'ADAS Spectral Line components (cumulative)'//
     >             ' (PH/M2/S)'
         endif

c
c        Scale by ABSFAC
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC

c
        numth2 = 0
c
        CALL LoadArray(touts2,numth2,touts         ,1,numthe)
c
        ngs      = 0
c
c       Excitation
c
         
        if (isele.ne.0) then 

           ngs = ngs + 1  
           elabs(ngs) = 'EXC Excitation'
c
           CALL LDADAS(CION2,IZMIN,ADASID,ADASYR,ADASEX,ISELE,0,0,
     >                 plrpad,Wlngth,IRCODE)
           CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
           CALL slLOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,
     .                AVPTS,CVALSA,0.0)
c
c           CALL LOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
c    .                 CVALSA,0.0,0)
c
           CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                   touts ,tvals (1,ngs),1,numthe)
        endif 
c
c       Recombination
c
        if (iselr.ne.0) then 

           ngs=ngs+1
           elabs(ngs) = 'REC Recombination'

c
           CALL LDADAS(CION2,IZMIN,ADASID,ADASYR,ADASEX,0,ISELR,0,
     >                  plrpad,Wlngth,IRCODE)
           CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
           CALL slLOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,
     .                AVPTS,CVALSA,0.0)
c
c           CALL LOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
c     .                CVALSA,0.0,0)
c
           CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                  touts ,tvals (1,ngs),1,numthe)
c
        endif
c
c       CX 
c
        if (iselx.ne.0) then 

           ngs=ngs+1
           elabs(ngs) = 'CX  CX'
c
           CALL LDADAS(CION2,IZMIN,ADASID,ADASYR,ADASEX,0,0,ISELX,
     >                  plrpad,Wlngth,IRCODE)
           CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >                MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
           CALL slLOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,
     .                AVPTS,CVALSA,0.0)
c
c           CALL LOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
c     .                CVALSA,0.0,0)
           CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                  touts ,tvals (1,ngs),1,numthe)
c
        endif
c
c       Create cumulative data
c
        if (ngs.gt.1) then 
           DO i1 = 2,ngs
             DO i2 = 1, numth2
               tvals2(i2,i1) = tvals2(i2,i1) + tvals2(i2,i1-1)
             ENDDO
           ENDDO
        endif 
c
        if (plotr) then
          call adjustout(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
        endif
c
c jdemod - added as hard-coded (for now) option to map as psin 
c        - changed to use numsmooth < -9 as indicator - set numsmooth to 0
c
c       Load PSIN data from file - needs input line - and then 
c       map the angle to a PSIN intersection - the PSIN load data
c       could be made optional
c
c       Get filename and load data for an externally specified PSIN mapping
c
        if (plot_psin) then 
c
            call get_psin_filename(ierr)
            if (ierr.ne.0) return 
c
c           Load PSIN reference data
c
            call read_psin_data
c
            do in = 1,numthe
c
               call calcpsin(robs,zobs,touts2(in),psin,targ)
c
               touts2(in) = psin
c
            end do 
c
            themin = touts2(1)
            themax = touts2(numthe)
            write(XLAB,'(''PSIN'')') 
c
c
         endif

c
c jdemod - end of plot 25
c
c

      ELSEIF ((plot.GE. 5.AND.plot.LE. 9).OR.
     .        (plot.GE.11.AND.plot.LE.15)) THEN
c...    Dalpha or Dgamma components (individual):
        ngs      = 0
        ref      = graph(7:46)
        elabs(1) = '    H (excitation)'
        elabs(2) = '    H+ (recombination)'
        elabs(3) = '    H2'
        elabs(4) = '    H2+'
        elabs(5) = '    H-'

        IF     (plot.GE. 5.AND.plot.LE. 9) THEN
          elabs(1) = elabs(plot-4)
          ylab = elabs(plot-4)
          ylab = ylab(5:LEN_TRIM(ylab))//'  Dalpha   PLRP'//units
          line = H_BALPHA
          i1   = plot - 4
        ELSEIF (plot.GE.11.AND.plot.LE.15) THEN
          elabs(1) = elabs(plot-10)
          ylab = elabs(plot-10)
          ylab = ylab(5:LEN_TRIM(ylab))//'  Dgamma   PLRP'//units
          line = H_BGAMMA
          i1   = plot - 10
        ENDIF

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

        ngs = ngs + 1
        CALL RVALKR(CVALSA,pinline(1,1,i1,line),1,1,1,FT,FP,
     >              MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
        CALL LOSINT(TVALS(1,ngs),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .              CVALSA,0.0,0)
        numthe1(ngs) = numthe
        DO i1 = 1, numthe
          touts1(i1,ngs) = touts(i1)
        ENDDO

c...    Plot experimental data:

        WRITE(0,*) '980:   TEMP! SETTING ISELD TO 3'
        iseld = 3

        if (iseld.ne.0) then
          ngs = ngs + 1
          numthe1(ngs) = 0
          CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                   MAXTHE,MAXNGS,MAXDATX,datatitle,1)
          elabs(ngs) = '    '//datatitle(1:LEN_TRIM(datatitle))
          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
        ENDIF

c...    Map data to master array:
        DO i1 = 1, ngs
          CALL MapArray(touts2      ,tvals2(1,i1),1,numth2,
     .                  touts1(1,i1),tvals (1,i1),1,numthe1(i1))
        ENDDO



        if (plotr) then

          STOP '980:10 NO LONGER GUARANTEED TO WORK'
 
          call adjustout(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
           write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
        endif

      ELSEIF (plot.EQ.16) THEN
c
c...    Dalpha without core contributions:
c
        ngs      = 1
        ref      = graph(7:46)
        elabs(1) = '    Dalpha (no core emission)'
        elabs(2) = '    experiment'

        ylab = 'Dalpha PLRP without core emission (ph m-2 s-1)'

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

        DO ir = irsep, nrs
          DO ik = 1, nks(ir)
            array(ik,ir) = pinalpha(ik,ir)
c            array(ik,ir) = pinline(ik,ir,6,H_BALPHA)			
          ENDDO
		ENDDO
         
        CALL RValKR(cvalsa,array,1,1,1,ft,fp,
     .              mfact,xxmin,xxmax,yymin,yymax,vmin,vmax)
        CALL LosInt(tvals(1,1),touts,twids,numthe,robs,zobs,avpts,
     .              CVALSA,0.0,0)
        CALL MapArray(touts2,tvals2(1,1),1,numth2,
     .                touts ,tvals (1,1),1,numthe)

        IF (plotr) THEN
          CALL AdjustOut(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          WRITE(xlab,'(''RADIUS (M) AT Z='',f6.3)') zadj
        ENDIF

c...    Plot experimental data:
        if (iseld.ne.0) then
          ngs = ngs + 1
          elabs(ngs) = '    experiment'
          call calc_expt(iseld,touts2,tvals2,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
        ENDIF

      ELSEIF (plot.EQ.17) THEN
c...    Generic LOS plot that uses the information in the 'ADAS'
c       data line in the OUT input file:
        ngs      = 1
        ref      = graph(7:46)
        elabs(1) = '    line'

        ylab = 'PLRP (ph m-2 s-1)'

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

        CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >              plrpad,Wlngth,IRCODE)
        CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     >              MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
        CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .              CVALSA,0.0,0)
        CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                touts ,tvals (1,1  ),1,numthe)
         
        IF (plotr) THEN
          CALL AdjustOut(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          WRITE(xlab,'(''RADIUS (M) AT Z='',f6.3)') zadj
        ENDIF

c...    Plot experimental data:
        if (iseld.ne.0) then
          ngs = ngs + 1
          elabs(ngs) = '    experiment'
          call calc_expt(iseld,touts2,tvals2,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
        ENDIF

      ELSEIF (plot.EQ.18) THEN
c...    Toroidal camera LOS plot:
        ngs      = 1
        ref      = graph(7:46)
        elabs(1) = '    line'

c        WRITE(0,*) '980:18 USING slLOSINT'

        ylab = 'PLRP (xxx)'

        numth2 = 0
        CALL LoadArray(touts2,numth2,touts,1,numthe)

        CALL StoreGrid(1)

        CALL LoadCamera(cdata)
       
        CALL RVALKR(CVALSA,cdata,1,1,1,FT,FP,
     .              MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
        CALL slLOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,ROBS,ZOBS,AVPTS,
     .                CVALSA,0.0)

        CALL MapArray(touts2,tvals2(1,ngs),1,numth2,
     .                touts ,tvals (1,1  ),1,numthe)

        CALL StoreGrid(2)
         
        IF (plotr) THEN
          CALL AdjustOut(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          WRITE(xlab,'(''RADIUS (M) AT Z='',f6.3)') zadj
        ENDIF

c...    Plot experimental data:
        if (iseld.ne.0) then
          ngs = ngs + 1
          elabs(ngs) = '    experiment'
          call calc_expt(iseld,touts2,tvals2,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
        ENDIF

      ELSEIF (plot.EQ.19.OR.plot.EQ.20.OR.plot.EQ.21.OR.plot.EQ.22.OR.
     .        plot.EQ.27) THEN
c...    Recomination weighted density integral:
        ngs      = 1
        ref      = graph(7:46)
        IF (elabs(1)(5:5).EQ.' ') elabs(1) = '    line'

        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 
        ENDDO
        plottype(1) = 55

        IF (graph(5:5).EQ.'3') slopt5 = 1
c...    Glorious hack!  Needed for Kbottom array on C-Mod, which is located in the plenum 
c       chamber.  This confuses the LOS plot code that counts wall intersections, since it assumes
c       that the array connot see outside the chamber, which is not the case of course:
        IF (graph(5:5).EQ.'4') slopt5 = 2

c        WRITE(0,*) '980:19 USING slLOSINT'

        xlab = 'degrees'

        status = .TRUE.
        DO WHILE (status)

c...      Process x-axis data:
          numthe1(ngs) = numthe
          DO i1 = 1, numthe
            touts1(i1,ngs) = touts(i1)
          ENDDO
          numth2 = 0
          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))

          IF (plot.EQ.27) THEN
c...        Load ADAS data:
            CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .                 ISELD,IERR)
            CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
            osmtmp = 0.0
            CALL RVALKR(osmtmp,plrpad,1,1,1,FT,FP,
     >                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
          ENDIF

c         DO ir = 1, nrs
c           DO ik = 1, nks(ir)
c             knbs(ik,ir) = 0.0
c             ktebs(ik,ir) = 0.0
c             pinrec(ik,ir) = 0.0
c           ENDDO
c         ENDDO
c         ir = 36
c         DO ik = 1, nks(ir)
c           knbs(ik,ir) = 1.0E+20
c           ktebs(ik,ir) = 1.0
c           pinrec(ik,ir) = 1.0 
c         ENDDO
c         ir = 30
c         DO ik = 1, nks(ir)
c           knbs(ik,ir) = 1.0E+21
c           ktebs(ik,ir) = 1.0
c           pinrec(ik,ir) = 0.1
c         ENDDO

c          DO ir = 14, 24
c            DO ik = 1, nks(ir)
c              pinline(ik,ir,6,H_BGAMMA) = 0.0
c            ENDDO
c          ENDDO



c          DO ir = 1, 26
c            DO ik = 1, nks(ir)
c              pinline(ik,ir,6,H_BGAMMA) = 0.0
c            ENDDO
c          ENDDO

c          DO ir = 27, 34
c            DO ik = 1, nks(ir)
c              pinline(ik,ir,6,H_BGAMMA) = 0.0
c            ENDDO
c          ENDDO
c          DO ir = 36, 38
c            DO ik = 1, nks(ir)
c              pinline(ik,ir,6,H_BGAMMA) = 0.0
c            ENDDO
c          ENDDO



c         ir = 29
c         DO ik = 1, nks(ir)
c           pinline(ik,ir,6,H_BGAMMA) = 0.0
c         ENDDO

c         ir = 35
c         DO ik = 1, nks(ir)
c           knbs(ik,ir) = 1.0E+21
c         ENDDO
       



c...      Apply weight function:
          losopt = 1
          maxwght0 = 0.0
c          DO ir = irtrap+1, nrs
          DO ir = irsep, nrs
c          DO ir = 1, nrs
            DO ik = 1, nks(ir)

c               ktebs(ik,ir) = 1.0
c               knbs(ik,ir) = 1.0

              IF     (plot.EQ.21.OR.plot.EQ.22) THEN
                maxwght0 = MAX(maxwght0,pinline(ik,ir,6,H_BGAMMA))
              ELSEIF (plot.EQ.27) THEN
                maxwght0 = MAX(maxwght0,osmtmp(ik,ir))
              ELSE
c                maxwght0 = MAX(maxwght0,(1.0E-20*pinrec(ik,ir))**2)
                maxwght0 = MAX(maxwght0,pinrec(ik,ir))
              ENDIF
            ENDDO
          ENDDO

          CALL RZero(wght0,MAXNKS*MAXNRS)

c          DO ir = irtrap+1, nrs
          DO ir = irsep, nrs
c          DO ir = 1, nrs
            DO ik = 1, nks(ir)
              IF     (plot.EQ.21.OR.plot.EQ.22) THEN
                wght0(ik,ir) = pinline(ik,ir,6,H_BGAMMA) / maxwght0
              ELSEIF (plot.EQ.27) THEN
                wght0(ik,ir) = osmtmp(ik,ir) / maxwght0
              ELSE
                wght0(ik,ir) = pinrec(ik,ir) / maxwght0
              ENDIF
            ENDDO
          ENDDO


          IF (plot.EQ.19.OR.plot.EQ.21.OR.plot.EQ.27) THEN
            ylab = 'n (m-3)'
 
            CALL RVALKR(CVALSA,knbs,1,1,1,FT,FP,
     .                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,ngs),TOUTS1(1,ngs),TWIDS,NUMTHe1(ngs),
     .                    ROBS,ZOBS,AVPTS,CVALSA,0.0)
          ELSE
            ylab = 'Te (eV)'

            CALL RVALKR(CVALSA,ktebs,1,1,1,FT,FP,
     .                  MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
            CALL slLOSINT(TVALS(1,ngs),TOUTS1(1,ngs),TWIDS,NUMTHe1(ngs),
     .                    ROBS,ZOBS,AVPTS,CVALSA,0.0)
          ENDIF

c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          status = .FALSE.
          READ(5,'(A128)',END=100) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,elabs(ngs+1),cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
            ngs = ngs + 1
          ELSE
            BACKSPACE 5
          ENDIF
c...      Cycle again, but for reflections:
c
c          jdemod - I think this should be end=100
c                 - otherwise the branch is outside the loop
c
c          READ(5,'(A128)',END=10) graph6
c
          READ(5,'(A128)',END=100) graph6
          IF (graph6(8:11).EQ.'Cycl'.OR.graph6(8:11).EQ.'CYCL'.OR.
     .        graph6(8:11).EQ.'cycl') THEN
            ngs = ngs + 1
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
100       CONTINUE

        ENDDO

c...    Restore main .raw file if necessary:
        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

c...    Convert to radial coordinate on x-axis (this assumes that the experimental:
        IF (plotr) THEN
          CALL AdjustOut(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          WRITE(xlab,'(''RADIUS (M) AT Z='',f6.3)') zadj
        ENDIF

c...    Plot experimental data:
        if (iseld.ne.0) then
          ngs = ngs + 1
          elabs(ngs) = '    experiment'
          CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                   MAXTHE,MAXNGS,MAXDATX,datatitle,1)
          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
          plottype(ngs) = -55
        ENDIF

c...    Copy the data to the plot arrays:
        DO i1 = 1, ngs
          CALL MapArray(touts2      ,tvals2(1,i1),1,numth2,
     .                  touts1(1,i1),tvals (1,i1),1,numthe1(i1))
        ENDDO

      ELSEIF (plot.EQ.98) THEN
c...    Recomination weighted density integral:
        ngs      = 1
        ref      = graph(7:46)
        IF (elabs(1)(5:5).EQ.' ') elabs(1) = '    line'

        slopt2 = 1
c        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 
        ENDDO
        plottype(1) = 55

        IF (graph(5:5).EQ.'3') slopt5 = 1
c...    Glorious hack!  Needed for Kbottom array on C-Mod, which is located in the plenum 
c       chamber.  This confuses the LOS plot code that counts wall intersections, since it assumes
c       that the array connot see outside the chamber, which is not the case of course:
        IF (graph(5:5).EQ.'4') slopt5 = 2

c        WRITE(0,*) '980:19 USING slLOSINT'

        xlab = 'degrees'
        ylab = '..........'


      WRITE(0,*) 'XLAB:'//xlab

c       DO ir = 1, nrs
c         DO ik = 1, nks(ir)
c           ktebs(ik,ir) = 0.5
c           ktibs(ik,ir) = 0.5
c           knbs(ik,ir) = 1.0E+21
c           pinatom(ik,ir) = 0.0
c         ENDDO
c       ENDDO


        status = .TRUE.
        DO WHILE (status)

c...      Process x-axis data:
          numthe1(ngs) = numthe
          DO i1 = 1, numthe
            touts1(i1,ngs) = touts(i1)
          ENDDO
          numth2 = 0
          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))


c...      Load ADAS data:
          plrpad = 0.0
          cvalsa = 0.0
          CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .               ISELD,IERR)
          CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                plrpad,Wlngth,IRCODE)
          CALL RVALKR(CVALSA,plrpad,1,1,1,FT,FP,
     .                MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)

          WRITE(0,*) 'ADAS DATA LOADED:',isele,iselr

c*CHECK*
c TRY DALPHA AND MAKE SURE IT IS THE SAME AS THE OTHER PLOTS

c...      Integrate:
          CALL slLOSINT(TVALS(1,ngs),TOUTS1(1,ngs),TWIDS,NUMTHe1(ngs),
     .                  ROBS,ZOBS,AVPTS,CVALSA,0.0)

c...      Check if additional data to be plotted from another DIVIMP case:
          status = .FALSE.
          READ(5,'(A128)',END=102) graph6
          IF (graph6(8:11).EQ.'Case'.OR.graph6(8:11).EQ.'CASE'.OR.
     .        graph6(8:11).EQ.'case') THEN
            READ(graph6,*) cdum1,elabs(ngs+1),cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
            ngs = ngs + 1
          ELSE
            BACKSPACE 5
          ENDIF

c...      More ADAS data specified, so keep going:
          READ(5,'(A128)',END=102) graph6
          IF (graph6(8:11).EQ.'Adas'.OR.graph6(8:11).EQ.'ADAS'.OR.
     .        graph6(8:11).EQ.'adas') THEN
            BACKSPACE 5
            ngs = ngs + 1
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF

c...      Cycle again, but for reflections:
          READ(5,'(A128)',END=102) graph6
          IF (graph6(8:11).EQ.'Cycl'.OR.graph6(8:11).EQ.'CYCL'.OR.
     .        graph6(8:11).EQ.'cycl') THEN
            ngs = ngs + 1
            status = .TRUE.
          ELSE
            BACKSPACE 5
          ENDIF
102       CONTINUE

        ENDDO

        READ(5,'(A128)',END=103) graph6
        IF (graph6(8:11).EQ.'Saha'.OR.graph6(8:11).EQ.'SAHA'.OR.
     .      graph6(8:11).EQ.'saha') THEN
          READ(graph6,*) cdum1,ln1,ln2
          WRITE(0,*) tvals(1,1)
          CALL FitTeSaha(MAXTHE,MAXNGS,tvals,numthe1(1),ngs,ln1,ln2)
          ngs = 1
        ELSE
          BACKSPACE 5
        ENDIF
 103    CONTINUE

c...    Restore main .raw file if necessary:
        IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

c...    Convert to radial coordinate on x-axis (this assumes that the experimental...):
        IF (plotr) THEN
          CALL AdjustOut(touts2,numthe,zadj,robs,zobs)
          themin = touts2(1)
          themax = touts2(numthe)
          WRITE(xlab,'(''RADIUS (M) AT Z='',f6.3)') zadj
        ENDIF

c...    Plot experimental data:
        if (iseld.ne.0) then
          ngs = ngs + 1
          elabs(ngs) = '    experiment'
          CALL load_expt(iseld,touts1(1,ngs),tvals(1,ngs),numthe1(ngs),
     .                   MAXTHE,MAXNGS,MAXDATX,datatitle,1)
          CALL LoadArray(touts2,numth2,touts1(1,ngs),1,numthe1(ngs))
          plottype(ngs) = -55
        ENDIF

c...    Copy the data to the plot arrays:
        DO i1 = 1, ngs
          CALL MapArray(touts2      ,tvals2(1,i1),1,numth2,
     .                  touts1(1,i1),tvals (1,i1),1,numthe1(i1))
        ENDDO

      ELSEIF (plot.EQ.99) THEN
        ref      = graph(7:46)
        YLAB     = 'emission ratio             '
        elabs(1) = '    Hg / Ha '
        elabs(2) = '    Hgamma  '
        elabs(3) = '    Halpha  '
        elabs(4) = '    Camera (scaled)'

      ELSE
        CALL ER('980','Unknown plot selected',*99)
      ENDIF

c...
      IF (slopt2.NE.1.AND.slopt2.NE.2) THEN
        WRITE (char(2),'(A,F6.2,A)') 'Z   = ',zobs,' m'
        WRITE (char(3),'(A,F6.2,A)') 'R   = ',robs,' m'
        WRITE (char(4),'(A,I6    )') 'NPT = ',avpts
        IF (theres.EQ.-99.0) THEN
          WRITE (char(5),'(A)'       ) 'INDIVIDUAL CHORD WIDTHS'
        ELSE 
          WRITE (char(5),'(A,F6.2,A)') 'RES = ',theres,' degrees'
        ENDIF
        WRITE (char(6),'(A,F6.2,A)') 'WL  = ',wlngth/10.0,' nm'
      ELSE
        char(1) = '                                   '
      ENDIF



c...  Normalize data, as required:
 104  READ(5,'(A80)',END=50) graph6
      IF   (graph6(8:11).EQ.'Norm'.OR.graph6(8:11).EQ.'NORM'.OR.
     .      graph6(8:11).EQ.'norm') THEN
        READ(graph6,*) cdum1,idum1,idum2

        CALL NormalizeData(numth2,touts2,tvals2(1,idum1),
     .                                   tvals2(1,idum2))

        GOTO 104
      ELSE
        BACKSPACE 5
      ENDIF


c      WRITE(0,*) INT(wlngth/10.0)

c      WRITE(6,*)'NGS= ',ngs
c      WRITE(0,*)'NGS= ',ngs

      DO ii = 1, numth2
        WRITE(6,'(A,I4,F10.4,1P,E10.2,10(E12.4:))')
     .    '980:   ',ii,touts2(ii),twids(ii),
     .    (tvals2(ii,i1),i1=1,ngs)
c        WRITE(0,'(A,I4,F10.4,1P,E10.2,10(E12.4:))')
c     .    '980???:   ',ii,touts2(ii),twids(ii),
c     .    (tvals2(ii,i1),i1=1,ngs)
      ENDDO

      WRITE(6,*) 'EXCEL OUTPUT:'
      DO i1 = 1, ngs
        DO ii = 1, numth2
          IF (tvals2(ii,i1).NE.LO)
     .      WRITE(6,'(A,2I4,F10.4,1P,E10.2,E12.4)')
     .       '980: ',i1,ii,touts2(ii),tvals2(ii,i1)
        ENDDO
      ENDDO


c...  Set the y-axis scale, and scale the data to be plotted:
      qmin = -HI
      qmax =  HI
      scale = 1.0
      READ(5,'(A80)',END=50) graph6
      IF   (graph6(8:11).EQ.'Scal'.OR.graph6(8:11).EQ.'SCAL'.OR.
     .      graph6(8:11).EQ.'scal') THEN
        READ(graph6,*) cdum1,qmin,qmax,scale

        IF (qmin.EQ.99.0) qmin = -HI
        IF (qmax.EQ.99.0) qmax =  HI

        WRITE(6,*) 
        WRITE(6,*) 'SCALED DATA:'

        DO i1 = 1, ngs
          DO ii = 1, numth2
            IF (tvals2(ii,i1).NE.LO) 
     .        tvals2(ii,i1) = tvals2(ii,i1) * scale

            IF (tvals2(ii,i1).NE.LO)
     .        WRITE(6,'(A,2I4,F10.4,1P,E10.2,E12.4)')
     .         '980: ',i1,ii,touts2(ii),tvals2(ii,i1)
          ENDDO
        ENDDO

      ELSE
        BACKSPACE 5
      ENDIF
50    CONTINUE

      title2 = title

      WRITE(0,*) 'XLAB:'//xlab

      CALL CustomisePlot(title2,xlab,ylab,elabs)

      WRITE(0,*) 'XLAB:'//xlab
     
      IF (slopt2.NE.1.AND.slopt2.NE.2)
     .  CALL SetPlotComments(980,job,extra_comments,0,tarshift(IKHI))

c...  Clear REF for now:
      IF (iopt_ghost.NE.0) REF = '                                   '


c. *************************** TEMP ********************
c        DO i1 = 1, numth2
c          tvals2(i1,1) = LO
c        ENDDO


      WRITE(0,*) 'XLAB:'//xlab


      IF (plot.EQ.4.OR.plot.EQ.10.OR.plot.EQ.24.or.plot.eq.25.OR.
     .    plot.EQ.26) THEN
c...    Cumulative plots:
        CALL DRAW(TOUTS2,TWIDS,TVALS2,MAXTHE,NUMTH2,ANLY,ngs,
     >            ISMOTH,THEMIN,THEMAX,qmin,qmax,IGNORS,ITEC,AVS,NAVS,
     >            JOB,TITLE2,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >            TABLE,8,iflag,1.0,0)

      ELSE

        CALL DRAW(TOUTS2,TWIDS,TVALS2,MAXTHE,NUMTH2,ANLY,ngs,
     >            ISMOTH,THEMIN,THEMAX,qmin,qmax,IGNORS,ITEC,AVS,NAVS,
     >            JOB,TITLE2,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >            TABLE,IOPT,iflag,1.0,0)
  
      ENDIF



c...  Plot error bars:
75    READ(5,'(A256)') graph1
      IF (graph1(8:11).EQ.'Bars'.OR.graph1(8:11).EQ.'BARS'.OR.
     .    graph1(8:11).EQ.'bars') THEN
        READ(graph1,*) cdum1,barmode,barseries,barindex,barerror1
        i2 = 0
        DO i1 = 1, numth2
          IF (tvals2(i1,barseries).NE.LO.AND.
     .        touts2(i1).GE.cxmin.AND.touts2(i1).LE.cxmax) THEN
            i2 = i2 + 1
            IF (barindex.EQ.-1.OR.barindex.EQ.i2) THEN
              IF     (barmode.EQ.1) THEN
c...            Relative error specified:
                delta1 = barerror1 * tvals2(i1,barseries)
                delta2 = delta1
              ELSEIF (barmode.EQ.2) THEN
c...            Absolute error specified:     
                delta1 = barerror1 * scale
                delta2 = delta1
              ELSEIF (barmode.EQ.3) THEN
                READ(graph1,*) cdum1,barmode,barseries,barindex,
     .                         barerror1,barerror2
c...            Relative error specified:
                delta1 = barerror1 * tvals2(i1,barseries)
                delta2 = barerror2 * tvals2(i1,barseries)
              ELSE
                CALL ER('982','Invalid BARMODE',*99)
              ENDIF 
c...          Draw error bar:
              CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
              CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)
              CALL POSITN(touts2(i1),tvals2(i1,barseries)+delta2)
              CALL JOIN  (touts2(i1),tvals2(i1,barseries)-delta1)        
              CALL POSITN(touts2(i1)-0.005*(CXMAX-CXMIN),
     .                    tvals2(i1,barseries)+delta2)
              CALL JOIN  (touts2(i1)+0.005*(CXMAX-CXMIN),
     .                    tvals2(i1,barseries)+delta2)        
              CALL POSITN(touts2(i1)-0.005*(CXMAX-CXMIN),
     .                    tvals2(i1,barseries)-delta1)
              CALL JOIN  (touts2(i1)+0.005*(CXMAX-CXMIN),
     .                    tvals2(i1,barseries)-delta1)        
            ENDIF
          ENDIF
        ENDDO
        GOTO 75
      ELSE
        BACKSPACE 5
      ENDIF






c...   Add picture to legend:
c      slopt = 2
c      CALL SLSET(0.95,1.31,0.38,0.73,xxmin,xxmax,yymin,yymax)
c      CALL DrawHalpha(ret,0)
c      CALL SUPIMP ('PARTIAL')



c...  Check if calculation of a norm is requested:

      nrmcalculate = .FALSE.
      READ(5,'(A256)',END=55) graph1
      IF (graph1(8:11).EQ.'Norm'.OR.graph1(8:11).EQ.'NORM'.OR.
     .    graph1(8:11).EQ.'norm') THEN
        nrmindex = nrmindex + 1
        READ(graph1,*) cdum1,nrmtype,nrmi1,nrmi2,nrmr1,nrmr2,
     .                 nrmcomment(nrmindex)
        nrmcalculate = .TRUE.
      ELSE
        BACKSPACE 5
      ENDIF
55    CONTINUE
      
c '000 Norm'  TYPE IND1 IND2 XRANGE1 XRANGE2 COMMENT
      
      IF (nrmcalculate) THEN
c...    Scan the selected range of plot data and store relevant data:
        CALL RZero(nrmdata,MAXTHE*2)

        i3 = 0

        DO i1 = 1, numth2 
          IF ((touts2(i1).GE.nrmr1.AND.touts2(i1).LE.nrmr2).OR.
     .        (nrmr1.EQ.0.0.AND.nrmr2.EQ.0.0)) THEN 
c...        Search for data in the specified range (x-axis data):            

            IF (tvals2(i1,nrmi1).NE.LO) THEN
c...          Primary data for norm found:
              i3 = i3 + 1
              nrmdata(i3,1) = tvals2(i1,nrmi1) 

c              WRITE(0,*) 'Searching:',i1,numth2,touts2(i1)

              IF (tvals2(i1,nrmi2).NE.LO) THEN 
c...            Secondary data assigned directly from plot data:
                nrmdata(i3,2) = tvals2(i1,nrmi2) 
              ELSE
c...            Secondary data linearly interpolated:
                DO i4 = i1-1, 1, -1
                  IF (tvals2(i4,nrmi2).NE.LO) THEN 
                    rdum1 = tvals2(i4,nrmi2)
                    EXIT
                  ENDIF
                ENDDO
                DO i5 = i1+1, numth2
                  IF (tvals2(i5,nrmi2).NE.LO) THEN 
                    rdum2 = tvals2(i5,nrmi2)
                    EXIT
                  ENDIF
                ENDDO
                IF (i4.EQ.0.OR.i5.EQ.numth2+1)
     .            CALL ER('982','Unable to interpolate secondary '//
     .                          'norm data',*99)

                frac = (touts2(i1) - touts2(i4)) /
     .                 (touts2(i5) - touts2(i4))
                nrmdata(i3,2) = rdum1 + frac * (rdum2 - rdum1)
              ENDIF
            ENDIF

          ENDIF
        ENDDO
      
c...    Store the number of data sets recorded:
        nrmnum(nrmindex) = i3

c...    Generate the norm:
        CALL CalculateNorm
      ENDIF


      RETURN
99    RETURN
c99    STOP
      END



c
c
c
c
c ======================================================================
c
c subroutine: CalculateNorm
c 
c
c
      SUBROUTINE CalculateNorm
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'

      INTEGER i1
      REAL    sum,count

      
      IF (nrmtype.EQ.1) THEN 
c...    ...
        nrmstep(nrmindex) = loadstep
 
        sum = 0.0
        count = 0.0
        DO i1 = 1, nrmnum(nrmindex)
          IF (nrmdata(i1,1).NE.LO.AND.nrmdata(i1,2).NE.LO) THEN
            count = count + 1.0
            sum = sum + ABS(nrmdata(i1,2) - nrmdata(i1,1)) / 
     .                  ABS(nrmdata(i1,1))
          ENDIF
        ENDDO
        IF (count.GT.0.0) THEN
          nrmvalue(nrmindex) = sum / count
        ELSE
          nrmvalue(nrmindex) = 0.0
        ENDIF

        nrmnum(nrmindex) = NINT(count)

        WRITE(6,'(A,3I6,1P,E12.4,0P,A)') 'NORM:',
     .    nrmindex,NINT(count),nrmstep(nrmindex),
     .    nrmvalue(nrmindex),
     .    nrmcomment(nrmindex)(1:LEN_TRIM(nrmcomment(nrmindex)))

      ELSE
        CALL ER('CalculateNorm','Invalid norm type',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE FitTeSaha(PARAM1,PARAM2,linedat,ndat,nline,line1,line2)
      IMPLICIT none
 
      INTEGER PARAM1,PARAM2,ndat,nline,line1,line2
      REAL    linedat(PARAM1,PARAM2)

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'cedge2d'
      INCLUDE 'slcom'

      INTEGER i1,i2
      REAL    A(5:9),tedat(5:9),tesaha(5:9),energy(5:9),
     .        minte,te,minlsd,lsd

c...  Einstein coefficients for Balmer lines:
      A(5) = 2.53E+06 
      A(6) = 9.73E+05
      A(7) = 4.39E+05
      A(8) = 2.22E+05
      A(9) = 1.22E+05      

c          WRITE(0,*) linedat(1,1)

      IF (line2-line1+1.NE.nline)
     .  CALL ER('FitTeSaha','Inconsistent line indices',*99)

      DO i1 = 1, ndat

        tedat = 0.0
        DO i2 = line1, line2
          tedat(i2) = linedat(i1,i2-line1+1) / A(i2)
        ENDDO
        DO i2 = line2, line1, -1
          tedat(i2) = tedat(i2) / tedat(line1)
        ENDDO



c       WRITE(0,'(A,5E10.2)') 'DAT:',(linedat(i1,i2-line1+1),i2=line1,line2)
c       WRITE(0,'(A,5F10.2)') 'DAT:',(tedat(i2),i2=line1,line2)

        minte  = 0.0
        minlsd = HI
        DO te = 0.10, 2.0, 0.01

c...      Find Saha distribution for TE:
          energy = 0.0
          tesaha = 0.0
          DO i2 = line1, line2
            energy(i2) = 13.6 / REAL(i2**2)
            tesaha(i2) = (REAL(i2) / REAL(line1))**2 * 
     .                   EXP((energy(i2) - energy(line1)) / te)
          ENDDO

c...      Find least squares difference:
          lsd = 0.0
          DO i2 = line1, line2
            lsd = lsd + (tedat(i2) - tesaha(i2))**2
          ENDDO
          lsd = SQRT(lsd)

c...      See if this TE is the best fit or not:
          IF (lsd.LT.minlsd) THEN
            minlsd = lsd
            minte  = te
          ENDIF



        ENDDO

        WRITE(0,'(A,I6,6(F10.3:))') 
     .    '-->',i1,minte,(tedat(i2),i2=line1,line2)
c      STOP 'gfshgr'

        linedat(i1,1) = minte

      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c
c Scale the data array VDAT1 so that its peak value equals the 
c corresponding value in VDAT2.
c
c
c
      SUBROUTINE NormalizeData(ndat,xdat,vdat1,vdat2)
      IMPLICIT none

      INCLUDE 'params'

      INTEGER ndat
      REAL    xdat(ndat),vdat1(ndat),vdat2(ndat)

      INTEGER i1,i1max
      REAL    x2a,x2b,v2a,v2b,factor,xval,val2

      i1max = 1
      DO i1 = 2, ndat
        IF (vdat1(i1).GT.vdat1(i1max)) i1max = i1
      ENDDO

c...  Find corresponding VDAT2 value:
      IF (vdat2(i1max).NE.LO) THEN
        factor = vdat2(i1max) / vdat1(i1max)
          WRITE(0,*) 'XVAL A:',factor
      ELSE
        xval = xdat(i1max)

        x2a = -HI
        x2b = -HI
        DO i1 = 1, ndat
          IF (xdat(i1).LT.xval.AND.vdat2(i1).NE.LO) THEN
            x2a = xdat(i1)
            v2a = vdat2(i1)
          ENDIF
        ENDDO
        DO i1 = ndat, 1, -1
          IF (xdat(i1).GT.xval.AND.vdat2(i1).NE.LO) THEN
            x2b = xdat(i1)
            v2b = vdat2(i1)
          ENDIF
        ENDDO

        IF (x2a.NE.-HI.AND.x2b.NE.-HI) THEN
          val2 = v2a + (xval - x2a) / (x2b - x2a) * (v2b - v2a)
          factor = val2 / vdat1(i1max)
          WRITE(0,*) 'XVAL B:',factor
        ELSE
          WRITE(0,*) 'XVAL:',xval,x2a,x2b
          CALL ER('NormalizeData','Scale factor could not be '//
     .            'determined',*99)
        ENDIF
      ENDIF


      DO i1 = 1, ndat
        IF (vdat1(i1).NE.LO) vdat1(i1) = vdat1(i1) * factor
      ENDDO

c      STOP 'sdfdssadfsd'


      RETURN
 99   STOP
      END
c
c ======================================================================
c
