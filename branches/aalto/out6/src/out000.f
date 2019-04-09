c     -*-Fortran-*-
c
      subroutine out000(iref,graph,iopt,ierr)
      use mod_params
      use mod_outcom
      use mod_cgeom
      use mod_comtor
      use mod_cneut2
      use mod_dynam2
      use mod_dynam3
      use mod_pindata
      use mod_outxy
      use mod_cedge2d
      use mod_printopt
      implicit none  
c
      integer iref,iopt,ierr
      character*(*) graph
c slmod begin
      INTEGER       idum1
      REAL          rdum1
      CHARACTER*256 dummy,cdum1

c slmod end
c
c     include 'params'
c     include 'outcom'
c
c     Other common blocks
c
c     include 'cgeom'
c     include 'comtor'
c     include 'cneut2'
c     include 'dynam2'
c     include 'dynam3'
c      include 'dynam4'
c     include 'pindata'
c      include 'cadas'
c      include 'grbound'
c     include 'outxy'
c     include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c     include 'printopt' 
c
c     Local Variables
c
      integer itarg
      integer startid,endid,stepid,switchid
      integer startin,endin,stepin
      integer inc
      integer iexpt
c
c     Counters
c

      integer i,j,k
      integer ix,iy
      integer ik,ir,iz,irref,ikref2
      integer jk,jr,jz
      integer id,jd,ii,in
      integer m
      real r,z




      real point 
      real SUM(10)

      real    ntotal


      real tdist(maxnds+2)
      real cwids(maxpts)

      real mach,lamii,lamti

      REAL KVALSPEC(MAXNKS*MSOLPT+msolpt+1,MAXNGS)
      REAL KWIDSPEC(MAXNKS*MSOLPT+msolpt+1)
      REAL KOUTSPEC(MAXNKS*MSOLPT+msolpt+1)

      REAL ROUTS(MAXNRS),RWIDS(MAXNRS),RVALS(MAXNRS,MAXNGS)

      real drn
      real xnvals(maxpts,maxngs)
      real xnouts(maxpts)

      integer i1

      real totsrc,tmpsum(maxnrs),tottmpsum,tmpsumiz



c...  Temp : looking at location of Thomson data:
      INTEGER    MAXCOLS   ,MAXDATA2,    DATAUNIT
      PARAMETER (MAXCOLS=10,MAXDATA2=500,DATAUNIT=13)
      INTEGER   eindex,etype,ndata,ncol
      LOGICAL   outofgrid
      REAL      xdata(MAXDATA2),edata(MAXDATA2,MAXCOLS)
c



        IF (IOPT.EQ.0) return
c
c       This is just a sanity check and not really necessary - especially since 
c       we keep adding plots that interpret various values of IOPT to perform
c       different functions.    
c
c        IF (IREF.NE.82.AND.IREF.NE.83.and.iref.ne.24.AND.
c     .      iref.NE.11.AND.iref.NE.12
c     >      .AND.IOPT.NE.1) return
c



      call init_plot(iref,graph,iopt)


C
C-----------------------------------------------------------------------
C     K CONTOUR PLOTS (FULL VIEW, CLOSE UPS, PLUS GRIDS, ETC)
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.11.OR.IREF.EQ.12) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        REF    = 'K CONTOURS' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
c slmod begin - new
        READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .        dummy(8:11).EQ.'SIZE') THEN
          WRITE(title,'(174X)')
          slopt4 = 1
          iopt_ghost = 1
c          CALL THICK2(4)
          CALL CTRMAG(10)
          READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
        ELSE
          CALL THICK2(1)
          BACKSPACE 5
        ENDIF

        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >     YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)

        IF     (iopt.EQ.20) THEN
          CALL SUPIMP ('TRIANGLES')
        ELSEIF (iopt.EQ.3.OR.iopt.EQ.4.OR.iopt.EQ.6.OR.iopt.EQ.8) THEN
          CALL SUPIMP ('PARTIAL')
        ELSE
          CALL SUPIMP ('FULL')
        ENDIF

        IF (iopt.EQ.2.OR.iopt.EQ.4.OR.iopt.EQ.5.OR.iopt.EQ.6.OR.
     .      iopt.EQ.7.OR.iopt.EQ.8) 
     .    CALL DrawAdditionalSurfaces(iopt)

        IF (graph(5:5).EQ.'1') THEN
c          DO i = 5, 20

c            i = 1
c
c           Read the index of the dataset from the GRAPH string 
c 
            read(graph(7:9),*) i

            eindex = i
 
            CALL Load_Expt_Data(DATAUNIT,eindex,xdata,etype,edata,
     .                          MAXCOLS,MAXDATA2,ndata,ncol,datatitle)

            WRITE(6,*) '== PLOT: 11 ================================='
            WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
            WRITE(6,*) 'INDEX  = ',eindex
            WRITE(6,*) 'TYPE   = ',etype
            WRITE(6,*) 'NDATA  = ',ndata
            WRITE(6,*) 'NCOL   = ',ncol

            call ctrmag(10)
            DO i1 = 1, ndata
              r = xdata(i1)
              z = edata(i1,1) 
              CALL PCSCEN(r,z,'+')
            ENDDO

c          ENDDO
        ENDIF

        IF (graph(5:5).EQ.'2') THEN

15        READ(5,'(A80)',END=110) graph1
          IF (graph1(8:11).EQ.'Thom'.OR.graph1(8:11).EQ.'THOM'.OR.
     .        graph1(8:11).EQ.'thom') THEN
 
c '000   Thom'  'DTS 500  '  99.00   0.00  01  94
            READ(graph1,*) cdum1,cdum1,rdum1,rdum1,idum1,eindex
 
            ndata = 0
            CALL Load_Expt_Data(DATAUNIT,eindex,xdata,etype,edata,
     .                          MAXCOLS,MAXDATA2,ndata,ncol,datatitle)

            WRITE(6,*) '== PLOT: 11 ================================='
            WRITE(6,*) 'TITLE  = ',datatitle(1:LEN_TRIM(datatitle))
            WRITE(6,*) 'INDEX  = ',eindex
            WRITE(6,*) 'TYPE   = ',etype
            WRITE(6,*) 'NDATA  = ',ndata
            WRITE(6,*) 'NCOL   = ',ncol

            call ctrmag(10)
            DO i1 = 1, ndata
              r = xdata(i1)
              z = edata(i1,1) 
              CALL PCSCEN(r,z,'+')
            ENDDO

c...        Look for more data:
            GOTO 15

          ELSE
            BACKSPACE 5
          ENDIF
        ENDIF


110     READ(5,'(A80)',END=115) graph1
        BACKSPACE 5
        IF (graph1(8:11).EQ.'View'.OR.graph1(8:11).EQ.'VIEW'.OR.
     .      graph1(8:11).EQ.'view') THEN
          call setup_col(1,5)
          CALL DrawHalpha(1,1)
          call setup_col(n_cols,col_opt) 
          GOTO 110
        ENDIF
115     CONTINUE

c       IPP/09 - Krieger - I think that has to go there ...
        call setup_col(n_cols,col_opt) 

        READ(5,'(A256)') dummy
        IF   (dummy(8:14).EQ.'Noframe'.OR.dummy(8:14).EQ.'noframe'.OR.
     .        dummy(8:14).EQ.'NOFRAME') THEN
        ELSE
          CALL FRAME
          BACKSPACE 5
        ENDIF

        IF (slopt4.EQ.1) THEN
c...      Finish off the plot:
          CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
          CALL MAP(0.0,1.0,0.0,1.0)
          CALL FULL
          CALL POSITN (0.0,1.0)
          CALL JOIN   (1.0,1.0)
          CALL JOIN   (1.0,0.0)
          CALL JOIN   (0.0,0.0)
          CALL JOIN   (0.0,1.0)
        ENDIF
c
c        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
c     >     YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)
c        CALL SUPIMP ('FULL')
c        CALL FRAME
c slmod end
      ENDIF
C
      IF (IREF.EQ.13.OR.IREF.EQ.14) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        ELABS(1) = ' '
        REF    = 'K CONTOURS & GRID' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >     YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)
        CALL SUPIMP ('FULL')
        DO 130 IX = 1, NXS
          DO 130 IY = 1, NYS
            IF (IFXYS(IX,IY).EQ.1) THEN
              R = (RMAX-RMIN) * REAL(IX)/REAL(NXS) + RMIN - DR
              Z = (ZMAX-ZMIN) * REAL(IY)/REAL(NYS) + ZMIN - DZ
              XVALS(1,1) = R
              XVALS(2,1) = R
              XVALS(3,1) = R + DR
              XVALS(4,1) = R + DR
              XVALS(5,1) = R
              XVALS(1,2) = Z
              XVALS(2,2) = Z + DZ
              XVALS(3,2) = Z + DZ
              XVALS(4,2) = Z
              XVALS(5,2) = Z
              CALL GRTRAC (XVALS(1,1),XVALS(1,2),5,ELABS(1),'LINE',0)
            ENDIF
  130   CONTINUE
        CALL FRAME
      ENDIF
c
      IF (IREF.EQ.15.OR.IREF.EQ.16) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        REF    = 'K CONTOURS' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >     YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)
        CALL SUPIMPOLD ('FULL')
        CALL FRAME
c       IPP/10 - Krieger - reset colors required
        call setup_col(n_cols,col_opt) 
      ENDIF
c
c     Contour Plot with wall segments numbered.
c
      IF (IREF.EQ.17.OR.IREF.EQ.18) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        REF    = 'K CONTOURS' // XPOINT
        PLANE  = 'WALL SEGMENTS LABELLED'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >     YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)
        CALL SUPIMP ('FULL')
        call label_wall
        CALL FRAME
c       IPP/10 - Krieger - reset colors required
        call setup_col(n_cols,col_opt) 
      ENDIF
C
C-----------------------------------------------------------------------
C     TEMP, DENSITY, K, SMAX ALONG REFERENCE LINE AND ALONG TARGET
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.21.OR.IREF.EQ.22.OR.IREF.EQ.23) THEN
        ELABS(1) = 'TE  TEB    '
        ELABS(2) = '  TITIB    '
        ELABS(3) = 'NB  NB     '
        ELABS(4) = '  K K      '
        if (iref.eq.23) then 
           ELABS(5) = 'FLUXFLUX   '
           ELABS(6) = ' Y  YIELD Y'
        else
           ELABS(5) = 'SMAXSMAX   '
           ELABS(6) = ' V  TOTAL V'
        endif
      ENDIF
C
      IF (IREF.EQ.21) THEN
        XLAB   = '   Z  (M) '
        YLAB   = '          '
        REF    = 'ALONG REFERENCE LINE'
        WRITE (IPLOT,9012) NPLOTS,REF
        JR = 0
        IK = IKREF
        CALL RZERO (RVALS, MAXNRS*MAXNGS)
c
        if (zs(ikref,irsep).lt.zs(ikref,irwall)) then
           startin = 1
           endin   = irwall
           stepin  = 1
        else
           startin = irwall
           endin   = 1
           stepin  = -1
        endif
c
        DO 210 IR = startin, endin, stepin
          JR = JR + 1
c
          jk = nks(ir)/2 +1
c
c          IF (IR.EQ.IRSEP-1) IK = IK - IKTO
c
          ROUTS(JR)   = ZS(IK,IR)
          RWIDS(JR)   = 0.5 * (KINDS(IK,IR) + KOUTDS(IK,IR))
          RVALS(JR,1) = KTEBS(IK,IR)
          RVALS(JR,2) = KTIBS(IK,IR)
          RVALS(JR,3) = KNBS(IK,IR)
          RVALS(JR,4) = KKS(IR)
          RVALS(JR,5) = KSMAXS(IR)
          RVALS(JR,6) = KTOTAS(IR)
  210   CONTINUE
c
        CALL DRAW (ROUTS,RWIDS,RVALS,MAXNRS,IRWALL,ANLY,6,99,
     >    ROUTS(1),ROUTS(IRWALL),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
      IF (IREF.EQ.22) THEN
        XLAB   = '   R  (M)'
        YLAB   = '           '
        REF    = 'ALONG TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin +1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif
c
        DO 220 ID = startid, endid, stepid
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          DVALS(JD,1) = KTEDS(ID)
          DVALS(JD,2) = KTIDS(ID)
          DVALS(JD,3) = KNDS(ID)
          DVALS(JD,4) = KKS(IR)
          DVALS(JD,5) = KSMAXS(IR)
          DVALS(JD,6) = KTOTAS(IR)
          IF (ID.EQ.switchid) JD = JD + 2
  220   CONTINUE
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    6,99,DOUTS(1),DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
c
c     Asdex targets
c
      IF (IREF.EQ.23) THEN

c       inner target

        XLAB   = ' SLUNT(M)'
        YLAB   = '           '
        REF    = 'ALONG INNER TARGET'
        WRITE (6,9012) NPLOTS,REF
c
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        write(6,*) 'Dist. along target, ne, Te, Ti'
        DO ID = nds,ndsin+1,-1
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          tdist(jd)=sqrt((rs(ik,ir)-1.132)**2+(zs(ik,ir)+0.822)**2)+.01
          DVALS(JD,1) = KTEDS(ID)
          DVALS(JD,2) = KTIDS(ID)
          DVALS(JD,3) = KNDS(ID)
          DVALS(JD,4) = KKS(IR)
          write(6,*) tdist(jd),knds(id),kteds(id),ktids(id)
        ENDDO
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,nds-ndsin,ANLY,
     >    4,99,tdist(1),tdist(nds-ndsin),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)

c       outer target

        XLAB   = ' SRUNT(M)'
        YLAB   = '           '
        REF    = 'ALONG OUTER TARGET'
        WRITE (6,9012) NPLOTS,REF
c
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        ELABS(5) = 'FLUXFLUX   '
        ELABS(6) = 'Y   YIELD  '
        write(6,*) 'Dist. along target, ne, Te, Ti'
        DO ID = ndsin,1,-1
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          tdist(jd)   = sqrt((rs(ik,ir)-1.487)**2+(zs(ik,ir)+0.954)**2)
          DVALS(JD,1) = KTEDS(ID)
          DVALS(JD,2) = KTIDS(ID)
          DVALS(JD,3) = KNDS(ID)
          DVALS(JD,4) = KKS(IR)
          DVALS(JD,5) = KFLUX(ID)
          DVALS(JD,6) = KYIELD(ID)
          write(6,*) tdist(jd),knds(id),kteds(id),ktids(id)
        ENDDO
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,ndsin,ANLY,
     >    6,99,tdist(1),tdist(ndsin),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
c
c     Heat fluxes across the target 
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.24) THEN
c
        if (iopt.lt.0) then
           iexpt = abs(iopt)
           ngs = 5
        else
           iexpt = 0
           ngs = 4
        endif 
c
        ELABS(1) = 'I   ION     '
        ELABS(2) = ' A  ATOM    '
        ELABS(3) = '  M MOLECULE'
        ELABS(4) = 'TOT TOTAL   '
        ELABS(5) = 'EXPTEXPT    '
        XLAB = '   R  (M) '
        YLAB = '          '
        REF  = 'HEAT FLUXES ACROSS THE TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
c
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin + 1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif
c
        DO ID = startid, endid, stepid
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          DVALS(JD,1) = targfluxdata(id,1,4)
          DVALS(JD,2) = targfluxdata(id,2,4)
          DVALS(JD,3) = targfluxdata(id,3,4)
          DVALS(JD,4) = targfluxdata(id,4,4)
          IF (ID.EQ.switchid) JD = JD + 2
        end do
c
c        if (iexpt.gt.0) then 
c           call calc_expt(iexpt,douts,dvals,maxnds+2,nds+2,
c     >               douts(1),douts(nds+2),maxngs,ngs,datatitle)
c        endif
c
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    ngs,99,DOUTS(1),DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,iexpt)
      ENDIF
C
C-----------------------------------------------------------------------
c
c     Ion Particle flux across the target 
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.25) THEN
c
c       Initialize
c
        call rzero(tdist,maxnds+2)  
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
c       Set scaling and target depending on iopt values 
c        
        if (iopt.lt.0) then
           ref   = 'PARTICLE FLUX ON '//inner//' TARGET'
           itarg = 1  
        else
           ref   = 'PARTICLE FLUX ON '//outer//' TARGET'
           itarg = 2 
        endif
c
        iopt = abs(iopt)
c
c       Axis - IOPT = 1 - PSIN, IOPT=2 - R-TARG
c      
        if (iopt.eq.1) then 
           XLAB = 'PSIN'
        elseif (iopt.eq.2) then
           XLAB = 'R ALONG TARGET (M)'
        endif 
c
        ELABS(1) = 'FLUXFLUX DENS'
        YLAB = 'FLUX DENSITY'

        WRITE (IPLOT,9012) NPLOTS,REF
c
c
c       Set up looping range for array based on selected axis
c        
c       PSIN
c        
        if (iopt.eq.1) then 
c
           if (itarg.eq.1) then
              startid = ndsin -1
              endid   = 2     
              stepid  = -1
           elseif (itarg.eq.2) then 
              startid = ndsin +2
              endid   = nds -1
              stepid  = 1
           endif
c
c       R axis
c
        elseif (iopt.eq.2) then 
c
           if (itarg.eq.1) then
              startid = 2
              endid   = ndsin-1     
              stepid  = 1
           elseif (itarg.eq.2) then 
              startid = ndsin +2
              endid   = nds-1
              stepid  = 1
           endif
c               
c          Switch ordering if necessary to get increasing R values   
c
           if (rp(startid).gt.rp(endid)) then
             switchid = startid
             startid  = endid
             endid    = switchid
             stepid   = -stepid
           endif
c
        endif 
c
c       Loop across target
c
        JD = 0
c
        write(6,*) 'PLOT 25:'
c
        DO ID = startid, endid, stepid
          JD = JD + 1
c
          IK = IKDS(ID)
          IR = IRDS(ID)
c
          if (iopt.eq.1) then 
             tdist(jd) = psitarg(ir,itarg)
          elseif(iopt.eq.2) then   
             tdist(jd) = rp(id)
          endif
          dvals(jd,1) = targfluxdata(id,1,1)
c
          write(6,'(i4,2(1x,g12.5))') jd,tdist(jd),dvals(jd,1)
c
        end do
c
c        if (iexpt.gt.0) then 
c           call calc_expt(iexpt,douts,dvals,maxnds+2,nds+2,
c     >               douts(1),douts(nds+2),maxngs,ngs,datatitle)
c        endif
c
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,abs(endid-startid)+1,ANLY,
     >    1,99,tdist(1),tdist(abs(endid-startid)+1),0.0,HI,IGNORS,ITEC,
     >    AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,iexpt)
      ENDIF
C
C-----------------------------------------------------------------------
c
C
C-----------------------------------------------------------------------
C     TEMPERATURES, DENSITIES & K AGAINST R, FOR PARTICULAR Z VALUES
C-----------------------------------------------------------------------
C
C     increased number of points, Krieger IPP/98
C
      IF (IREF.EQ.31) THEN
        ELABS(1) = '    TEB    '
        ELABS(2) = '    TIB    '
        ELABS(3) = '    NB     '
        ELABS(4) = '    K      '
        XLAB   = '   R  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(F7.3)') POINT
        WRITE (REF,'(''SECTION ALONG Z ='',F7.4)') POINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (TVALS, MAXTHE*MAXNGS)
        ik = 0
        ir = 0
        drn = max((rmax-rmin)/(maxthe-1),0.001)
c
*       added output to file for comparison with ivp; Krieger, IPP/95
        write(IPLOT,'(1x,a     )')
     >    'R       Te        Ti        Ne          K      IR'
c
        DO 310 IX = 1,maxthe
          R = rmin + drn * (ix-1)
          touts(ix) = r
          twids(ix) = 0.0
          call gridpos (ik,ir,r,point,.false.,griderr)
          if (griderr) goto 310
c
          TVALS(IX,1) = KTEBS(IK,IR)
          TVALS(IX,2) = KTIBS(IK,IR)
          TVALS(IX,3) = KNBS (IK,IR)
          TVALS(IX,4) = KKS  (IR)
c
*         added output; Krieger, IPP/95
          write(iplot,'(1x,f5.3,2(3x,f7.2),3x,e9.3,3x,f4.2,3x,i2)')
     >      touts(ix),tvals(ix,1),tvals(ix,2),tvals(ix,3),tvals(ix,4),
     >      ir
c
  310   CONTINUE
        CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,MAXTHE,ANLY,
     >    4,99,TOUTS(1),TOUTS(maxthe),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
c
      ENDIF
C
C-----------------------------------------------------------------------
C     TEMPERATURES, DENSITIES & K AGAINST R, FOR IVP z-value, IPP 6/95
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.32) THEN
        ELABS(1) = '    TEB    '
        ELABS(2) = '    TIB    '
        ELABS(3) = '    NB     '
        ELABS(4) = '    K      '
        XLAB   = '   R  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(F7.3)') POINT
        IY   = (POINT-ZMIN) / DZ + 1
        WRITE (REF,'(''SECTION ALONG Z ='',F7.4)') POINT
        WRITE (6,9012) NPLOTS,REF
        CALL RZERO (XNVALS, MAXpts*MAXNGS)
        ik = 0
        ir = 0
        drn = (rmax-rmin)/(maxpts-1)
*       added output to file for comparison with ivp; Krieger, IPP/95
        write(6,'(1x,a     )') 'r        Te        Ti     ne'
        DO 311 IX = 1,maxpts
           R = rmin + drn * (ix-1)
           xnouts(ix) = r
c
c          IK = IKXYS(IX,IY)
c          IR = IRXYS(IX,IY)
c          IF (IFXYS(IX,IY).EQ.0) GOTO 310
c
          call gridpos (ik,ir,r,point,.false.,griderr)
          if (griderr) goto 311
c
          XnVALS(IX,1) = KTEBS(IK,IR)
          XnVALS(IX,2) = KTIBS(IK,IR)
          XnVALS(IX,3) = KNBS (IK,IR)
          XnVALS(IX,4) = KKS  (IR)
*         added output to file for comparison with ivp; Krieger, IPP/95
          write(6,'(1x,f5.3,2(3x,f6.2),3x,e9.3)')
     >      r,xnvals(ix,1),xnvals(ix,2),xnvals(ix,3)
  311   CONTINUE
        CALL DRAW (XnOUTS,XWIDS,XnVALS,MAXpts,maxpts,ANLY,
     >    4,99,1.2,1.8,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
c
c kkmod
c
C                                                                       
C-----------------------------------------------------------------------
C     Impurity density against R at given z-Value and along AUG target
C     Krieger IPP 2/98
C     (use method of plot option 36)
C-----------------------------------------------------------------------
C                                                                       
      IF (IREF.EQ.33) THEN                                              
        ELABS(1) = '    NZ     '                                        
        ELABS(2) = '    FZ     '
        ELABS(3) = '    K      '                                        
        XLAB   = '   R  (M)'                                            
        YLAB   = '           '                                          
        READ (GRAPH(38:44),'(F7.3)') POINT                              
C
C       LOAD 2-D ARRAY WITH RESPECTIVE VALUES AT EACH LOCATION
C
c       We want here the "real" density, not per injected
c       particle -> multiply by absfac
c
        call rzero(plastmp,maxnks*maxnrs)
        do iz = 0,nizs
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 if (knbs(ik,ir).ne.0.0) then 
                    plastmp(ik,ir) = plastmp(ik,ir) 
     >                    + sdlims(ik,ir,iz)*absfac
                 endif
              end do 
           end do
        end do 
        CALL CUT(TOUTS,TVALS(1,1),NUMTHE,MAXTHE,PLASTMP,
     >           RMIN,POINT,RMAX,POINT)
c
c       We want here the "real" concentration, not per injected
c       particle -> multiply by absfac
c
        call rzero(plastmp,maxnks*maxnrs)
        do iz = 0,nizs
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 if (knbs(ik,ir).ne.0.0) then 
                    plastmp(ik,ir) = plastmp(ik,ir) 
     >                    + sdlims(ik,ir,iz)/knbs(ik,ir)*absfac
                 endif
              end do 
           end do
        end do 
        CALL CUT(TOUTS,TVALS(1,2),NUMTHE,MAXTHE,PLASTMP,
     >           RMIN,POINT,RMAX,POINT)

        call rzero(plastmp,maxnks*maxnrs)
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = KKS(IR)
          ENDDO
        ENDDO
        CALL CUT(TOUTS,TVALS(1,3),NUMTHE,MAXTHE,PLASTMP,
     >           RMIN,POINT,RMAX,POINT)

        WRITE (REF,'(''SECTION ALONG Z ='',F7.4)') POINT                
        WRITE (IPLOT,9012) NPLOTS,REF                                       
*       added output to file; Krieger, IPP/95
        write(iplot,'(1x,a     )') 'R       Nz          Fz          K'
        DO II = 1, NUMTHE
          TOUTS(II) = RMIN + TOUTS(II)
          TWIDS(II) = 0.0
*         added output; Krieger, IPP/98
          call gridpos (ik,ir,touts(ii),point,.false.,griderr)
          write(iplot,'(1x,f5.3,1p,2(3x,e9.3),3x,0p,f4.2)')
     >      touts(ii),tvals(ii,1),tvals(ii,2),tvals(ii,3)
        ENDDO
        CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,
     >    3,99,RMIN,RMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF                                                             
c
c kkmod
c
c
C-----------------------------------------------------------------------
c
c     JET version of plot 31 using the CUT subroutine - from LDH
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.36) THEN
        ELABS(1) = '    TEB    '
        ELABS(2) = '    TIB    '
        ELABS(3) = '    NB     '
        ELABS(4) = '    K      '
        XLAB   = '   R  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(F7.3)') POINT
        CALL CUT(TOUTS,TVALS(1,1),NUMTHE,MAXTHE,KTEBS,
     >           RMIN,POINT,RMAX,POINT)
        CALL CUT(TOUTS,TVALS(1,2),NUMTHE,MAXTHE,KTIBS,
     >           RMIN,POINT,RMAX,POINT)
        CALL CUT(TOUTS,TVALS(1,3),NUMTHE,MAXTHE,KNBS,
     >           RMIN,POINT,RMAX,POINT)
C
C        LOAD 2-D ARRAY WITH K-VALUES AT EACH LOCATION
C
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = KKS(IR)
          ENDDO
        ENDDO
        CALL CUT(TOUTS,TVALS(1,4),NUMTHE,MAXTHE,PLASTMP,
     >           RMIN,POINT,RMAX,POINT)
C
        WRITE (REF,'(''SECTION ALONG Z ='',F7.4)') POINT
        WRITE (IPLOT,9012) NPLOTS,REF
        DO 316 II = 1, NUMTHE
          TOUTS(II) = RMIN + TOUTS(II)
          TWIDS(II) = 0.0
  316   CONTINUE
        CALL DRAW (TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,
     >    4,99,RMIN,RMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     DISTANCES, Y/A, SLOPE FACTORS ETC AGAINST S ALONG CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.41.OR.IREF.EQ.42) THEN
        ELABS(1) = 'R   R      '
        ELABS(2) = ' Z  Z      '
        ELABS(3) = '  B BRATIO '
        ELABS(4) = '   VVOLUME '
        IF (IREF.EQ.41) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        DO 410 IK = 1, NKS(IR)
          IF (IREF.EQ.41) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          KVALS(IK,1) = RS(IK,IR)
          KVALS(IK,2) = ZS(IK,IR)
          KVALS(IK,3) = KBFS(IK,IR)
          KVALS(IK,4) = KAREAS(IK,IR)
  410   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     ELECTRON and ION MEAN FREE PATH AGAINST S ALONG CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.43.OR.IREF.EQ.44) THEN
c
c       Produces two plots - one for ions and one for electrons.
c
c
c       PLOT The ION Mean free path
c
        ELABS(1) = 'LTI LAMTI  '
        ELABS(2) = '5LTI5 LAMTI'
        ELABS(3) = 'SCTITI SLEN'
c
c        ELABS(4) = '500X500 LTI'
c
        IF (IREF.EQ.43) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'ION MEAN FREE PATH/SCALE LENGTH (M)'
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ION MFP: '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (kvals, maxnks*maxngs)
c
        DO IK = 1, NKS(IR)
          IF (IREF.EQ.43) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = lmfpii(ik,ir)
          KVALS(IK,2) = 5.0 * lmfpii(ik,ir)
          KVALS(IK,3) = lgradti(ik,ir)
c
c          KVALS(IK,4) = 500.0 * lmfpii(ik,ir)
c
        end do
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    3,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >           3,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,
     >           IGNORS,ITEC,AVS,NAVS,
     >           JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,
     >           PLANE,TABLE,1,6,1.0,0)
c
           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >           3,99,mgnd*kouts(nks(ir)),kouts(nks(ir)),-HI,HI,
     >           IGNORS,ITEC,AVS,NAVS,
     >           JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,
     >           PLANE,TABLE,1,6,1.0,0)
c
        endif

c
c       Plot the Electron MEAN FREE PATH
c
        ELABS(1) = 'LTE LAMTE  '
        ELABS(2) = '5LTE5 LAMTE'
        ELABS(3) = 'SCTETE SLEN'
c
c        ELABS(4) = '500X500 LTE'
c
        IF (IREF.EQ.43) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = 'ELEC.MEAN FREE PATH/SCALE LENGTH (M)'
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ELECTRON MFP:'',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (kvals, maxnks*maxngs)
c
        DO IK = 1, NKS(IR)
          IF (IREF.EQ.43) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = lmfpee(ik,ir)
          KVALS(IK,2) = 5.0 * lmfpee(ik,ir)
          KVALS(IK,3) = lgradte(ik,ir)
c
c          KVALS(IK,4) = 500.0 * lmfpee(ik,ir)
c
        end do
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    3,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >           3,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,
     >           IGNORS,ITEC,AVS,NAVS,
     >           JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,
     >           PLANE,TABLE,1,6,1.0,0)
c
           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >           3,99,mgnd*kouts(nks(ir)),kouts(nks(ir)),-HI,HI,
     >           IGNORS,ITEC,AVS,NAVS,
     >           JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,
     >           PLANE,TABLE,1,6,1.0,0)
c
        endif

      ENDIF
c
C-----------------------------------------------------------------------
c     Comparative plots of DIVIMP and EDGE2D values along the
c     field lines vs. S and P.
C-----------------------------------------------------------------------
c
c     Plot of fluxes - IREF = 45
c
c
      IF (IREF.EQ.45) THEN
        ELABS(1) = 'GaE GaE'
        ELABS(2) = 'GaD GaD'
c        ELABS(3) = 'TEE TEE    '
c        ELABS(4) = 'TED TED    '
c        ELABS(5) = 'NBE NBE    '
c        ELABS(6) = 'NBD NBD    '
c
        IF (IREF.EQ.45) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c        write (6,*) 'Plot 45:',iref,iopt,ir
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.45) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
           KVALS(1,1) = KvdS(idds(ir,2)) * knds(idds(ir,2))
           KVALS(1,2) = KvdS(idds(ir,2)) * knds(idds(ir,2))
c
c           KVALS(1,3) = KTEdS(idds(ir,2))
c           KVALS(1,4) = KTEdS(idds(ir,2))
c           KVALS(1,5) = KNdS (idds(ir,2))
c           KVALS(1,6) = KNdS (idds(ir,2))
c
           KVALS(nks(ir)+2,1) = KvdS(idds(ir,1)) * knds(idds(ir,1))
           KVALS(nks(ir)+2,2) = KvdS(idds(ir,1)) * knds(idds(ir,1))

c
        endif
c
c
        DO IK = 1, NKS(IR)
          IF (IREF.eq.45) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK+in,1) = e2dvhs(ik,ir) * e2dnbs(ik,ir)
          KVALS(IK+in,2) = (kvhs(ik,ir)/qtim) * knbs(ik,ir)
c
        enddo
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
      ENDIF

c
c     Plot 46 - Only Ionization Density - non-normalised - vs. S
c

C
C-----------------------------------------------------------------------
C     PIN IONIZATION FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.46) THEN
        IF (IREF.EQ.46) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'H IONIZ'
c
        ELABS(1) = 'HizEHizE   '
        ELABS(2) = 'HizDHizD   '
c
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        PLTMAX = -HI
        PLTMIN = HI
        DO 461 IK = 1, NKS(IR)
          IF (IREF.EQ.46) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = E2DION(IK,IR)
          KVALS(IK,2) = PINION(IK,IR)
c
c          KVALS(IK,1) = PINION(IK,IR)*KVOLS(IK,IR)
c
C          IF (KVALS(IK,1).GT.0.0) THEN
C             KVALS(IK,1) = LOG10(KVALS(IK,1))
C             PLTMAX = MAX(PLTMAX,KVALS(IK,1))
C             PLTMIN = MIN(PLTMIN,KVALS(IK,1))
C          ENDIF
c
           PLTMAX = MAX(PLTMAX,KVALS(IK,1))
           PLTMIN = MIN(PLTMIN,KVALS(IK,1))
  461   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)),PLTMIN,PLTMAX,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,KOUTS(1),mgst*kouts(nks(ir)),pltmin,pltmax,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),pltmin,pltmax,
     >    IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
       ENDIF

c
c     Plot 47 - Only Temperatures - non-normalised - vs. S
c
      IF (IREF.EQ.47) THEN
        ELABS(1) = 'TIE TIE'
        ELABS(2) = 'TID TID'
        ELABS(3) = 'TEE TEE    '
        ELABS(4) = 'TED TED    '
c        ELABS(5) = 'NBE NBE    '
c        ELABS(6) = 'NBD NBD    '
c
        IF (IREF.EQ.47) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c        write (6,*) 'Plot 45:',iref,iopt,ir
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.47) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
           KVALS(1,1) = KTIdS(idds(ir,2))
           KVALS(1,2) = KTIdS(idds(ir,2))
           KVALS(1,3) = KTEdS(idds(ir,2))
           KVALS(1,4) = KTEdS(idds(ir,2))
c           KVALS(1,5) = KNdS (idds(ir,2))
c           KVALS(1,6) = KNdS (idds(ir,2))
c
           KVALS(nks(ir)+2,1) = KTIdS(idds(ir,1))
           KVALS(nks(ir)+2,2) = KTIdS(idds(ir,1))
           KVALS(nks(ir)+2,3) = KTEdS(idds(ir,1))
           KVALS(nks(ir)+2,4) = KTEdS(idds(ir,1))
c           KVALS(nks(ir)+2,5) = KNdS (idds(ir,1))
c           KVALS(nks(ir)+2,6) = KNdS (idds(ir,1))
c
        endif
c
c
        DO IK = 1, NKS(IR)
          IF (IREF.eq.47) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK+in,1) = e2dtibs(ik,ir)
          KVALS(IK+in,2) = ktibs(ik,ir)
          KVALS(IK+in,3) = e2dtebs(ik,ir)
          KVALS(IK+in,4) = ktebs(ik,ir)
c          KVALS(IK+in,5) = e2dnbs(ik,ir)
c          KVALS(IK+in,6) = KNBS (IK,IR)
c
        enddo
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    4,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
      ENDIF
c
c     Plot 48 - Only Densities - non-normalised - vs. S
c
      IF (IREF.EQ.48) THEN
c        ELABS(1) = 'TIE TIE'
c        ELABS(2) = 'TID TID'
c        ELABS(3) = 'TEE TEE    '
c        ELABS(4) = 'TED TED    '
        ELABS(1) = 'NBE NBE    '
        ELABS(2) = 'NBD NBD    '
c
        IF (IREF.EQ.48) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c        write (6,*) 'Plot 45:',iref,iopt,ir
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.48) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
c           KVALS(1,1) = KTIdS(idds(ir,2))
c           KVALS(1,2) = KTIdS(idds(ir,2))
c           KVALS(1,3) = KTEdS(idds(ir,2))
c           KVALS(1,4) = KTEdS(idds(ir,2))
           KVALS(1,1) = KNdS (idds(ir,2))
           KVALS(1,2) = KNdS (idds(ir,2))
c
c           KVALS(nks(ir)+2,1) = KTIdS(idds(ir,1))
c           KVALS(nks(ir)+2,2) = KTIdS(idds(ir,1))
c           KVALS(nks(ir)+2,3) = KTEdS(idds(ir,1))
c           KVALS(nks(ir)+2,4) = KTEdS(idds(ir,1))
           KVALS(nks(ir)+2,1) = KNdS (idds(ir,1))
           KVALS(nks(ir)+2,2) = KNdS (idds(ir,1))
c
        endif
c
c
        DO IK = 1, NKS(IR)
          IF (IREF.eq.48) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
c          KVALS(IK+in,1) = e2dtibs(ik,ir)
c          KVALS(IK+in,2) = ktibs(ik,ir)
c          KVALS(IK+in,3) = e2dtebs(ik,ir)
c          KVALS(IK+in,4) = ktebs(ik,ir)
          KVALS(IK+in,1) = e2dnbs(ik,ir)
          KVALS(IK+in,2) = KNBS (IK,IR)
c
        enddo
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
      ENDIF

c
c     Plot 49 - Only Velocities - non-normalised - vs. S
c
      IF (IREF.EQ.49) THEN
c        ELABS(1) = 'TIE TIE'
c        ELABS(2) = 'TID TID'
c        ELABS(3) = 'TEE TEE    '
c        ELABS(4) = 'TED TED    '
        ELABS(1) = 'VBE VBE    '
        ELABS(2) = 'VBD VBD    '
c
        IF (IREF.EQ.49) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.49) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
c           KVALS(1,1) = KTIdS(idds(ir,2))
c           KVALS(1,2) = KTIdS(idds(ir,2))
c           KVALS(1,3) = KTEdS(idds(ir,2))
c           KVALS(1,4) = KTEdS(idds(ir,2))
           KVALS(1,1) = KvdS (idds(ir,2))
           KVALS(1,2) = KvdS (idds(ir,2))
c
c           KVALS(nks(ir)+2,1) = KTIdS(idds(ir,1))
c           KVALS(nks(ir)+2,2) = KTIdS(idds(ir,1))
c           KVALS(nks(ir)+2,3) = KTEdS(idds(ir,1))
c           KVALS(nks(ir)+2,4) = KTEdS(idds(ir,1))
           KVALS(nks(ir)+2,1) = KvdS (idds(ir,1))
           KVALS(nks(ir)+2,2) = KvdS (idds(ir,1))
c
        endif
c
c
        DO IK = 1, NKS(IR)
          IF (IREF.eq.49) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
c          KVALS(IK+in,1) = e2dtibs(ik,ir)
c          KVALS(IK+in,2) = ktibs(ik,ir)
c          KVALS(IK+in,3) = e2dtebs(ik,ir)
c          KVALS(IK+in,4) = ktebs(ik,ir)
          KVALS(IK+in,1) = e2dvhs(ik,ir)
          KVALS(IK+in,2) = KvhS (IK,IR) / qtim
c
        enddo
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
      ENDIF


C
C-----------------------------------------------------------------------
C     CONDUCTION AND CONVECTION TERMS IN BACKGROUND PLASMA CALCULATIONS
C-----------------------------------------------------------------------
c
c     Power Balance
C
      IF (IREF.EQ.51) THEN
        ELABS(1) = '1   5NVKT  '
        ELABS(2) = ' 2  CONDUCT'
        ELABS(3) = '  3 1/2MV3N'
        ELABS(4) = '   4(1)+(2)'
        ELABS(5) = '5   (1)/(4)'
        ELABS(6) = ' 6  (3)/(4)'
        XLAB = '   S  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,
     >     '(''Power Balance - Ring '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        DO 450 IK = 1, NKS(IR)
          KOUTS(IK) = KSS(IK,IR)
          KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          KVALS(IK,1) = 5.0*KNBS(IK,IR)*KVHS(IK,IR)/QTIM
     >                  *ECH*KTEBS(IK,IR)
          IF (IK.EQ.1) THEN
            KVALS(IK,2)=-CK0*KTEBS(IK,IR)**2.5
     >                    *(KTEBS(1,IR)-KTEBS(2,IR))
     >                    /(KSS(1,IR)-KSS(2,IR))
          ELSEIF (IK.EQ.NKS(IR)) THEN
            KVALS(IK,2)=-CK0*KTEBS(IK,IR)**2.5*(KTEBS(NKS(IR)-1,IR)
     >                    -KTEBS(NKS(IR),IR))
     >                    /(KSS (NKS(IR)-1,IR)-KSS(NKS(IR),IR) )
          ELSE
            KVALS(IK,2) = -CK0*KTEBS(IK,IR)**2.5*( (KTEBS(IK,IR)
     >                  -KTEBS(IK+1,IR)) / (KSS(IK,IR)-KSS(IK+1,IR))
     >                  + (KTEBS(IK-1,IR) - KTEBS(IK,IR)) /
     >                  (KSS(IK-1,IR)-KSS(IK,IR)) )/2.0
          ENDIF
          KVALS(IK,3) = 0.5*CRMB*AMU*(KVHS(IK,IR)/QTIM)**3*knbs(ik,ir)
          KVALS(IK,4) = KVALS(IK,1) + KVALS(IK,2)
          KVALS(IK,5) = KVALS(IK,1)/KVALS(ik,4)
          kvals(ik,6) = kvals(ik,3)/kvals(ik,4)
  450   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    6,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then

           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >       6,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,
     >       IGNORS,ITEC,AVS,NAVS,
     >       JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
           CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >       6,99,mgnd*kouts(nks(ir)),KOUTS(NKS(IR)),-HI,HI,
     >       IGNORS,ITEC,AVS,NAVS,
     >       JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
         endif

c
      ENDIF
C
c     Electron Power Balance
c

      IF (IREF.EQ.52) THEN
        ELABS(1) = '1   5/2NVKT'
        ELABS(2) = ' 2  CONDUCT'
C        ELABS(3) = '  3 1/2MV3N'
C        ELABS(4) = '   4(1)+(2)'
C        ELABS(5) = '5   (1)/(4)'
C        ELABS(6) = ' 6  (3)/(4)'
        XLAB = '   S  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,
     >    '(''E-Power Balance-Ring '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        DO 455 IK = 1, NKS(IR)
          KOUTS(IK) = KSS(IK,IR)
          KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          KVALS(IK,1) = 2.5*KNBS(IK,IR)*KVHS(IK,IR)/QTIM
     >                  *ECH*KTEBS(IK,IR)
          IF (IK.EQ.1) THEN
            KVALS(IK,2)=-CK0*KTEBS(IK,IR)**2.5*(KTEBS(1,IR)-KTEBS(2,IR))
     >                    /(KSS(1,IR)-KSS(2,IR))
          ELSEIF (IK.EQ.NKS(IR)) THEN
            KVALS(IK,2)=-CK0*KTEBS(IK,IR)**2.5*(KTEBS(NKS(IR)-1,IR)
     >                    -KTEBS(NKS(IR),IR))
     >                    /(KSS (NKS(IR)-1,IR)-KSS(NKS(IR),IR) )
          ELSE
            KVALS(IK,2) = -CK0*KTEBS(IK,IR)**2.5*( (KTEBS(IK,IR)
     >                  -KTEBS(IK+1,IR)) / (KSS(IK,IR)-KSS(IK+1,IR))
     >                  + (KTEBS(IK-1,IR) - KTEBS(IK,IR)) /
     >                  (KSS(IK-1,IR)-KSS(IK,IR)) )/2.0
          ENDIF
C          KVALS(IK,3) = 0.5*CRMB*AMU*(KVHS(IK,IR)/QTIM)**3*knbs(ik,ir)
C          KVALS(IK,4) = KVALS(IK,1) + KVALS(IK,2)
C          KVALS(IK,5) = KVALS(IK,1)/KVALS(ik,4)
C          kvals(ik,6) = kvals(ik,3)/kvals(ik,4)
  455   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,mgnd*kouts(nks(ir)),KOUTS(NKS(IR)),-HI,HI,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
      ENDIF
C
c     Ion Power Balance
C
      IF (IREF.EQ.53) THEN
        ELABS(1) = '1   5/2NVKT'
        ELABS(2) = ' 2  CONDUCT'
        ELABS(3) = '  3 1/2MV3N'
        ELABS(4) = '   4SUM    '
c        ELABS(5) = '5   (1)/(4)'
c        ELABS(6) = ' 6  (2)/(4)'
c        ELABS(7) = '  7 (3)/(4)'
        XLAB = '   S  (M)'
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,
     >   '(''I-Power Balance-Ring '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        DO 460 IK = 1, NKS(IR)
          KOUTS(IK) = KSS(IK,IR)
          KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          KVALS(IK,1) = 2.5*KNBS(IK,IR)*KVHS(IK,IR)/QTIM
     >                  *ECH*KTIBS(IK,IR)
          IF (IK.EQ.1) THEN
            KVALS(IK,2)=-CK0i*KTiBS(IK,IR)**2.5
     >                    *(KTiBS(1,IR)-KTiBS(2,IR))
     >                    /(KSS(1,IR)-KSS(2,IR))
          ELSEIF (IK.EQ.NKS(IR)) THEN
            KVALS(IK,2)=-CK0i*KTiBS(IK,IR)**2.5*(KTiBS(NKS(IR)-1,IR)
     >                    -KTiBS(NKS(IR),IR))
     >                    /(KSS (NKS(IR)-1,IR)-KSS(NKS(IR),IR) )
          ELSE
            KVALS(IK,2) = -CK0i*KTiBS(IK,IR)**2.5*( (KTiBS(IK,IR)
     >                  -KTiBS(IK+1,IR)) / (KSS(IK,IR)-KSS(IK+1,IR))
     >                  + (KTiBS(IK-1,IR) - KTiBS(IK,IR)) /
     >                  (KSS(IK-1,IR)-KSS(IK,IR)) )/2.0
          ENDIF
          KVALS(IK,3) = 0.5*CRMB*AMU*(KVHS(IK,IR)/QTIM)**3*knbs(ik,ir)
          KVALS(IK,4) = KVALS(IK,1) + KVALS(IK,2) + kvals(ik,3)
C          KVALS(IK,5) = KVALS(IK,1)/KVALS(ik,4)
C          kvals(ik,6) = kvals(ik,3)/kvals(ik,4)
  460   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    4,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    4,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    4,99,mgnd*kouts(nks(ir)),KOUTS(NKS(IR)),-HI,HI,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
      ENDIF

C
C-----------------------------------------------------------------------
C     SPOL AND Z AS A FUNCTION OF S FOR A GIVEN RING.
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.61) THEN
        ELABS(1) = 'SPOLSPOL  M'
        ELABS(2) = 'Z   Z     M'
        IF (IREF.EQ.61) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = '   POSITION (M)     '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        VMAX = 0.0
        CALL rzero (kvals, maxnks*maxngs)
        DO 610 IK = 1, NKS(IR)
          IF (IREF.EQ.61) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          KVALS(IK,1) = KPS(IK,IR)
          KVALS(IK,2) = ZS(IK,IR)
  610   CONTINUE
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
      ENDIF
c
C-----------------------------------------------------------------------
c     Plots 65 and 66, 67 and 68 are plots of retention criterion
c     functions for the specified ring. These retention functions are
c     defined in TN959.  66 and 68 for plots vs. P are roughed-in but
c     not completely supported ... the dT/dS term would need to be
c     recalculated.
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.65) THEN
c        ELABS(1) = 'LII  LAM II'
c        ELABS(2) = 'LTI  LAM TI'
        ELABS(1) = 'I/T LII/LTI'
        ELABS(2) = 'MACHMACH   '
c
c       Include Plot name information
c
        plane = 'Neuhauser Retention Criterion Plot'
c
        IF (IREF.EQ.65) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.65) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
           lamii  = 1.0e16 * KTIDS(idds(ir,2))**2/KNDS(idds(ir,2))
           lamti  = 1.0 / abs(1.0/KTIDS(idds(ir,2))  *
     >                  (ktibs(1,ir)-ktids(idds(ir,2)))/
     >                        (kss(1,ir)))
           KVALS(1,1) = lamii/lamti
           KVALS(1,2) = abs(kvds(idds(ir,2)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,2))+KTIdS(idds(ir,2)))
     >           *(1+RIZB)/CRMB))
c
           lamii = 1.0e16 * KTIDS(idds(ir,1))**2
     >                          /  KNDS(idds(ir,1))
           lamti = 1.0/abs(1.0/KTIDS(idds(ir,1))  *
     >                        (ktids(idds(ir,1))-ktibs(nks(ir),ir))/
     >                        (ksmaxs(ir)-kss(nks(ir),ir)))
           KVALS(nks(ir)+2,1) = lamii/lamti
           KVALS(nks(ir)+2,2) = abs(kvds(idds(ir,1)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
     >              *(1+RIZB)/CRMB))
c
        endif
c
c
        DO 650 IK = 1, NKS(IR)
          IF (IREF.EQ.65) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          lamii = 1.0e16 * KTIBS(ik,ir)**2/KNBS(ik,ir)
c
c         Calculate the LamTi term which involves a required estimate
c         for dT/dS ... this results in the complicated code which
c         follows.
c
          if (ik.eq.1) then
             if (kss(1,ir).eq.0.0) then
                lamti = abs( KTIBS(ik,ir)  *
     >          ((kss(2,ir)-kss(1,ir))/(ktibs(2,ir)-ktibs(1,ir))))
             else
                 lamti  = abs( KTIBS(ik,ir)  *
     >            0.5*((kss(2,ir)-kss(1,ir))/(ktibs(2,ir)-ktibs(1,ir))+
     >               (kss(1,ir))/(ktibs(1,ir)-ktids(idds(ir,2))) ))
             endif
          elseif (ik.eq.nks(ir)) then
             if (kss(1,ir).eq.0.0) then
                lamti = abs(KTIBS(ik,ir)  *
     >              ( (kss(nks(ir),ir)-kss(nks(ir)-1,ir))
     >                /(ktibs(nks(ir),ir)-ktibs(nks(ir)-1,ir))))
             else
                lamti = abs(KTIBS(ik,ir)  *
     >           0.5 * (  (kss(nks(ir),ir)-kss(nks(ir)-1,ir))
     >                  /(ktibs(nks(ir),ir)-ktibs(nks(ir)-1,ir) )
     >                 +
     >                  (ksmaxs(ir)-kss(nks(ir),ir))
     >                / (ktids(idds(ir,1))-ktibs(nks(ir),ir)) ))
             endif
          else
             lamti = abs( KTIBS(ik,ir)  *
     >          0.5 *( (kss(ik,ir)-kss(ik-1,ir))
     >                 / (ktibs(ik,ir)-ktibs(ik-1,ir) )
     >              +  (kss(ik+1,ir)-kss(ik,ir))
     >                 /  (ktibs(ik+1,ir)-ktibs(ik,ir)) ))
          endif
c
          KVALS(ik+in,1) = lamii/lamti
          KVALS(ik+in,2) = abs((kvhs(ik,ir)/qtim) /
     >      SQRT(0.5*EMI*(KTEBS(ik,ir)+KTIBS(ik,ir))*(1+RIZB)/CRMB))
c
  650   CONTINUE
C
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),mgst*kouts(nks(ir)+inc),-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,mgnd*kouts(nks(ir)+inc),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
        plane = ' '
c
      ENDIF
c
C-----------------------------------------------------------------------
c
c     Plot of LamFig and LamFF.
c
C-----------------------------------------------------------------------
c

      IF (IREF.EQ.67) THEN
        ELABS(1) = 'LFF 2LAM FF'
        ELABS(2) = 'LIG LAM FIG '
c
c       Include Plot name information
c
        plane =  'Neuhauser Retention Criterion Plot'
c
        IF (IREF.EQ.67) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
c
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.67) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           else
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
           mach = abs(kvds(idds(ir,2)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,2))+KTIdS(idds(ir,2)))
     >          * (1+RIZB)/CRMB))
           if (mach.le.5.0e-2) mach = 5.0e-2
c
           KVALS(1,1) = 3.7e14 * KTIDS(idds(ir,2))**2/KNDS(idds(ir,2))
     >                              / mach  * 2.0
           KVALS(1,2) = 0.031 * abs(KTIDS(idds(ir,2))
     >                        * (kss(1,ir))
     >                     / (ktibs(1,ir)-ktids(idds(ir,2))) )
c
c
           mach = kvds(idds(ir,1)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
     >          *(1+RIZB)/CRMB)
           if (mach.le.5.0e-2) mach = 5.0e-2
c
           KVALS(nks(ir)+2,1) = 3.7e14 * KTIDS(idds(ir,1))**2
     >                          / KNDS(idds(ir,1)) / mach * 2.0
           KVALS(nks(ir)+2,2) = 0.031 * abs(KTIDS(idds(ir,1))  *
     >                         (ksmaxs(ir)-kss(nks(ir),ir))
     >                / (ktids(idds(ir,1))-ktibs(nks(ir),ir)) )
c
        endif
c
c      Note: the Mach number is limited to a value of .01 to
c      prevent LAM TI from becoming INF when M -> 0
c
        DO 670 IK = 1, NKS(IR)
          IF (IREF.EQ.67) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          mach = abs((kvhs(ik,ir)/qtim) /
     >      SQRT(0.5*EMI*(KTEBS(ik,ir)+KTIBS(ik,ir))*(1+RIZB)/CRMB))
          if (mach.le.5.0e-2) mach = 5.0e-2
c
          KVALS(ik+in,1) = 3.7e14 * KTIBS(ik,ir)**2/KNBS(ik,ir) / mach
     >                       * 2.0
c
c         Calculate the LamFig term which involves a required estimate
c         for dT/dS ... this results in the complicated code which
c         follows.
c
          if (ik.eq.1) then
             if (kss(1,ir).eq.0.0) then
                KVALS(ik+in,2) = abs( KTIBS(ik,ir)  *
     >          ((kss(2,ir)-kss(1,ir))/(ktibs(2,ir)-ktibs(1,ir))))
             else
                KVALS(ik+in,2) = abs( KTIBS(ik,ir)  *
     >            0.5*((kss(2,ir)-kss(1,ir))/(ktibs(2,ir)-ktibs(1,ir))+
     >               (kss(1,ir))/(ktibs(1,ir)-ktids(idds(ir,2))) ))
             endif
          elseif (ik.eq.nks(ir)) then
             if (kss(1,ir).eq.0.0) then
                KVALS(ik+in,2) = abs(KTIBS(ik,ir)  *
     >              ( (kss(nks(ir),ir)-kss(nks(ir)-1,ir))
     >                /(ktibs(nks(ir),ir)-ktibs(nks(ir)-1,ir))))
             else
                KVALS(ik+in,2) = abs(KTIBS(ik,ir)  *
     >           0.5 * (  (kss(nks(ir),ir)-kss(nks(ir)-1,ir))
     >                  /(ktibs(nks(ir),ir)-ktibs(nks(ir)-1,ir) )
     >                 +
     >                  (ksmaxs(ir)-kss(nks(ir),ir))
     >                / (ktids(idds(ir,1))-ktibs(nks(ir),ir)) ))
             endif
          else
             KVALS(ik+in,2) = abs( KTIBS(ik,ir)  *
     >          0.5 *( (kss(ik,ir)-kss(ik-1,ir))
     >                 / (ktibs(ik,ir)-ktibs(ik-1,ir) )
     >              +  (kss(ik+1,ir)-kss(ik,ir))
     >                 /  (ktibs(ik+1,ir)-ktibs(ik,ir)) ))
          endif
c
c         Since LAM Fig is simply 0.031 X LAM Ti ... LAM Ti is calculate
c         and the constant multiplication is done below.
c
          kvals(ik+in,2) = 0.031 * kvals(ik+in,2)

c
  670   CONTINUE
C
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,KOUTS(1),mgst*kouts(nks(ir)+inc),-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,maxnks,NKS(IR)+inc,ANLY,
     >    2,99,mgnd*kouts(nks(ir)+inc),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        endif
c
        plane = ' '
c
      ENDIF
C
C-----------------------------------------------------------------------
C     ELECTRIC FIELD AND DRIFT VELOCITY FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.71.OR.IREF.EQ.72) THEN
        ELABS(1) = 'VB  VB  M/S'
        ELABS(2) = '   EE   V/M'
C
C        ELABS(3) = '100EE*100  '
C
        ELABS(3) = 'TE  TEB    '
        ELABS(4) = '  TITIB    '
        ELABS(5) = 'NB  NB     '
        ELABS(6) = 'FE  FEG1   '
        ELABS(7) = '  FIFIG1   '
        ELABS(8) = 'MACHMACH   '
        IF (IREF.EQ.71) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0
           if (iref.eq.71) then
              KWIDS(1) = kss(1,ir)
              kouts(nks(ir)+2) = ksmaxs(ir)
              kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
           elseif (iref.eq.72) then
              KWIDS(1) = kps(1,ir)
              kouts(nks(ir)+2) = kpmaxs(ir)
              kwids(nks(ir)+2) = kpmaxs(ir) - kps(nks(ir),ir)
           endif
c
           KVALS(1,1) = KVDS(idds(ir,2))
           KVALS(1,2) = KEDS(idds(ir,2))
           KVALS(1,3) = KTEdS(idds(ir,2))
           KVALS(1,4) = KTIdS(idds(ir,2))
           KVALS(1,5) = KNdS (idds(ir,2))
           KVALS(1,6) = KFEdS(idds(ir,2))*KALPHS(1)/(QTIM*QTIM*EMI/CRMI)
           KVALS(1,7) = KFIdS(idds(ir,2))*KBETAS(1)/(QTIM*QTIM*EMI/CRMI)
           KVALS(1,8) = KVALS(1,1) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,2))+KTIdS(idds(ir,2)))
     >         *(1+RIZB)/CRMB)
c
           KVALS(nks(ir)+2,1) = KVDS(idds(ir,1))
           KVALS(nks(ir)+2,2) = KEDS(idds(ir,1))
           KVALS(nks(ir)+2,3) = KTEdS(idds(ir,1))
           KVALS(nks(ir)+2,4) = KTIdS(idds(ir,1))
           KVALS(nks(ir)+2,5) = KNdS (idds(ir,1))
           KVALS(nks(ir)+2,6) = -KFEdS(idds(ir,1))*KALPHS(1)
     >                          / (QTIM*QTIM*EMI/CRMI)
           KVALS(nks(ir)+2,7) = -KFIdS(idds(ir,1))*KBETAS(1)
     >                          / (QTIM*QTIM*EMI/CRMI)
           KVALS(nks(ir)+2,8) = KVALS(nks(ir)+2,1) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
     >           *(1+RIZB)/CRMB)
c
        endif
c
c
        DO 710 IK = 1, NKS(IR)
          IF (IREF.EQ.71) THEN
            KOUTS(IK+in) = KSS(IK,IR)
            KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK+in) = KPS(IK,IR)
            KWIDS(IK+in) = 0.0
            IF (IK.GT.1) KWIDS(IK+in) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK+in) = KWIDS(IK+in) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          KVALS(IK+in,1) = KVHS(IK,IR) / QTIM
          KVALS(IK+in,2) = KES(IK,IR) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK+in,3) = KTEBS(IK,IR)
          KVALS(IK+in,4) = KTIBS(IK,IR)
          KVALS(IK+in,5) = KNBS (IK,IR)
          KVALS(IK+in,6) = KFEGS(IK,IR)*KALPHS(1) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK+in,7) = KFIGS(IK,IR)*KBETAS(1) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK+in,8) = KVALS(IK+in,1) /
     >    SQRT(0.5*EMI*(KTEBS(IK,IR)+KTIBS(IK,IR))*(1+RIZB)/CRMB)
  710   CONTINUE
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    8,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    8,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    8,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     ELECTRIC FIELD AND DRIFT VELOCITY - WIDE SCALE ONLY
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.73.OR.IREF.EQ.74) THEN
        ELABS(1) = 'VB  VB  M/S'
        ELABS(2) = '   EE   V/M'
C
C        ELABS(3) = '100EE*100  '
C
        ELABS(3) = 'TE  TEB    '
        ELABS(4) = '  TITIB    '
        ELABS(5) = 'NB  NB     '
        ELABS(6) = 'FE  FEG1   '
        ELABS(7) = '  FIFIG1   '
        ELABS(8) = 'MACHMACH   '
        IF (IREF.EQ.73) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        VMAX = 0.0
        CALL rzero (kvals, maxnks*maxngs)
        DO 711 IK = 1, NKS(IR)
          IF (IREF.EQ.73) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
          KVALS(IK,1) = KVHS(IK,IR) / QTIM
          KVALS(IK,2) = KES(IK,IR) / (QTIM*QTIM*EMI/CRMI)
          VMAX = MAX (VMAX, ABS(KVALS(IK,2)))
          KVALS(IK,3) = KTEBS(IK,IR)
          KVALS(IK,4) = KTIBS(IK,IR)
          KVALS(IK,5) = KNBS (IK,IR)
          KVALS(IK,6) = KFEGS(IK,IR)*KALPHS(1) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK,7) = KFIGS(IK,IR)*KBETAS(1) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK,8) = KVALS(IK,1) /
     >    SQRT(0.5*EMI*(KTEBS(IK,IR)+KTIBS(IK,IR))*(1+RIZB)/CRMB)
  711   CONTINUE
C
C        DO 716 IK = 1, NKS(IR)
C          KVALS(IK,3) = MAX (-VMAX, MIN (VMAX, 100.0*KVALS(IK,2)))
C  716   CONTINUE
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    8,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN IONIZATION FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.75) THEN
        ELABS(1) = 'H+  H+ IONIZATION'
        IF (IREF.EQ.75) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'PIN H IONIZ'
C        YLAB   = 'LOG10(IONIZATION)'
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        PLTMAX = -HI
        PLTMIN = HI
        DO 750 IK = 1, NKS(IR)
          IF (IREF.EQ.75) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = PINION(IK,IR)
c
c          KVALS(IK,1) = PINION(IK,IR)*KVOLS(IK,IR)
c
C          IF (KVALS(IK,1).GT.0.0) THEN
C             KVALS(IK,1) = LOG10(KVALS(IK,1))
C             PLTMAX = MAX(PLTMAX,KVALS(IK,1))
C             PLTMIN = MIN(PLTMIN,KVALS(IK,1))
C          ENDIF
c
           PLTMAX = MAX(PLTMAX,KVALS(IK,1))
           PLTMIN = MIN(PLTMIN,KVALS(IK,1))
  750   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),PLTMIN,PLTMAX,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN Hydrogen Atom Density FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.76) THEN
        ELABS(1) = 'H0  H0 DENSITY'
        IF (IREF.EQ.76) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'PIN H NEUTRAL'
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        PLTMAX = -HI
        PLTMIN = HI
        DO 760 IK = 1, NKS(IR)
          IF (IREF.EQ.76) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = PINATOM(IK,IR)
c
           PLTMAX = MAX(PLTMAX,KVALS(IK,1))
           PLTMIN = MIN(PLTMIN,KVALS(IK,1))
  760   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),PLTMIN,PLTMAX,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN Hydrogen Molecular Density FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.77) THEN
        ELABS(1) = 'H2  H2 DENSITY'
        IF (IREF.EQ.77) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'PIN H2 DENSITY'
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        PLTMAX = -HI
        PLTMIN = HI
        DO 770 IK = 1, NKS(IR)
          IF (IREF.EQ.77) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = PINMOL(IK,IR)
c
          PLTMAX = MAX(PLTMAX,KVALS(IK,1))
          PLTMIN = MIN(PLTMIN,KVALS(IK,1))
  770   CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),PLTMIN,PLTMAX,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
       ENDIF
C
C-----------------------------------------------------------------------
C     ELECTRIC FIELD AND DRIFT VELOCITY FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
C     High-resolution data
C
C
      IF (IREF.EQ.78) THEN
        ELABS(1) = 'VB  VB  M/S'
C
C        ELABS(2) = '   EE   V/M'
C
C        ELABS(3) = '100EE*100  '
C
        ELABS(2) = 'TE  TEB    '
        ELABS(3) = '  TITIB    '
        ELABS(4) = 'NB  NB     '
C        ELABS(5) = 'FE  FEG1   '
C        ELABS(7) = '  FIFIG1   '
        ELABS(5) = 'MACHMACH   '
        IF (IREF.EQ.78) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = '           '
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''HI-RES BG : RING '',I3,'', K'',F7.4)')
     >         CIRHR,KKS(CIRHR)
        WRITE (IPLOT,9012) NPLOTS,REF
        VMAX = 0.0
        CALL RZERO (KVALSPEC,(MAXNKS*MSOLPT+MSOLPT+1)*maxngs)
        write(6,*) 'CIR',cirhr
        DO 720 IK = 1, NKS(CIRHR)*msolpt+msolpt+1
          IF (IREF.EQ.78) THEN
            KOUTSPEC(IK) = solcor(ik-1)
            if (ik.eq.1) then
              KWIDSPEC(IK) = (solcor(ik)-solcor(IK-1))
            elseif (ik.eq.nks(cirhr)*msolpt+msolpt+1) then
              KWIDSPEC(IK) = (solcor(ik-1)-solcor(IK-2))
            else
              KWIDSPEC(IK) = 0.5 * (solcor(ik)-solcor(IK-2))
            endif
          ENDIF
          KVALSPEC(IK,1) = solvel(IK-1)
          KVALSPEC(IK,2) = solte(IK-1)
          KVALSPEC(IK,3) = solti(IK-1)
          KVALSPEC(IK,4) = solne (IK-1)
          KVALSPEC(IK,5) = KVALSPEC(IK,1) /
     >    SQRT(0.5*EMI*(solte(IK-1)+solti(IK-1))*(1+RIZB)/CRMB)
c          write(6,*) 'r:',ik,kouts(ik),kwids(ik),kvalspec(ik,1),
c     > kvalspec(ik,2),kvalspec(ik,3),kvalspec(ik,4),kvalspec(ik,5)
  720   CONTINUE
C
c       introduced for AUG
c
c        enldist=0.5
c
        CALL DRAW (KOUTSPEC,KWIDSPEC,KVALSPEC,MAXNKS*MSOLPT+msolpt+1,
     >    NKS(cIRhr)*MSOLPT+msolpt+1,ANLY,
     >    5,99,KOUTSPEC(1),KOUTSPEC(NKS(cIRhr)*MSOLPT+msolpt+1),-HI,HI,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then
c


        CALL DRAW (KOUTSPEC,KWIDSPEC,KVALSPEC,MAXNKS*MSOLPT+msolpt+1,
     >    NKS(cIRhr)*MSOLPT+msolpt+1,ANLY,
     >    5,99,KOUTSPEC(1),mgst*koutspec(nks(cirhr)*msolpt+msolpt+1),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTSPEC,KWIDSPEC,KVALSPEC,MAXNKS*MSOLPT+msolpt+1,
     >    NKS(cirhr)*MSOLPT+msolpt+1,ANLY,
     >    5,99,mgnd*KOUTSPEC(NKS(cIRhr)*MSOLPT+msolpt+1),
     >    KOUTSPEC(NKS(cIRhr)*MSOLPT+msolpt+1),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     Te,Ti,N,Vb,E for selected contours
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.79) THEN
c        ELABS(1) = '    VB  M/S'
        ELABS(1) = '    MACH   '

        ELABS(2) = '    E   V/M'
        ELABS(3) = '    TEB    '
        ELABS(4) = '    TIB    '
        ELABS(5) = '    NB     '

        XLAB = '   S  (M)'
        YLAB   = 'Scaled Plasma Quantities'

        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in = 0
           inc = 0
        else
           in = 1
           inc = 2
           kouts(1) = 0.0

           KWIDS(1) = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
c
c           KVALS(1,1) = KVDS(idds(ir,2))
c
           KVALS(1,1) = KVDS(idds(ir,2)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,2))+KTIdS(idds(ir,2)))
     >         *(1+RIZB)/CRMB)
c

           KVALS(1,2) = KEDS(idds(ir,2))
           KVALS(1,3) = KTEdS(idds(ir,2))
           KVALS(1,4) = KTIdS(idds(ir,2))
           KVALS(1,5) = KNdS (idds(ir,2))

c           KVALS(nks(ir)+2,1) = KVDS(idds(ir,1))
           KVALS(nks(ir)+2,1) = KVDS(idds(ir,1)) /
     >      SQRT(0.5*EMI*(KTEdS(idds(ir,1))+KTIdS(idds(ir,1)))
     >           *(1+RIZB)/CRMB)
c
           KVALS(nks(ir)+2,2) = KEDS(idds(ir,1))
           KVALS(nks(ir)+2,3) = KTEdS(idds(ir,1))
           KVALS(nks(ir)+2,4) = KTIdS(idds(ir,1))
           KVALS(nks(ir)+2,5) = KNdS (idds(ir,1))
        endif
c
c
        DO 714 IK = 1, NKS(IR)
c
          KOUTS(IK+in) = KSS(IK,IR)
          KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
c
c          KVALS(IK+in,1) = KVHS(IK,IR) / QTIM
          KVALS(IK+in,1) = (KVHS(IK,IR) / QTIM) /
     >    SQRT(0.5*EMI*(KTEBS(IK,IR)+KTIBS(IK,IR))*(1+RIZB)/CRMB)
c
          KVALS(IK+in,2) = KES(IK,IR) / (QTIM*QTIM*EMI/CRMI)
          KVALS(IK+in,3) = KTEBS(IK,IR)
          KVALS(IK+in,4) = KTIBS(IK,IR)
          KVALS(IK+in,5) = KNBS (IK,IR)
  714   CONTINUE
C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    5,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,6,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    5,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,6,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    5,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >    -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,6,6,1.0,0)
c
        endif
c
      ENDIF
C
C-----------------------------------------------------------------------
C     DEPOSITION AND NET EROSION QUANTITIES, ETC
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.81) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
c
        NGS    = cngs
c
        REF    = 'WALL DEPOSITION  (TOTAL)'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
c       IPP/09 - Krieger - better zero array before use
        call rzero(plastmp,maxnks*maxnrs)
        DO 801 IR = 1,NRS
          DO 801 IK = 1,NKS(IR)
            PLASTMP(IK,IR) = WALLS(IK,IR,MAXIZS+1)
801     CONTINUE
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.82) THEN
        XLAB = '   R  (M)  '
        YLAB = '           '
        REF  = 'DEPOSITION'
        JD   = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin + 1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif

C       Schmid IPP/06 - Dump DEPS to ascii
        open(99, file='TargetDeposition.dat', status='replace')
        write(6,*) 'Storing TargetDeposition.dat'
        write(99,*) 'ID, IZ, DEPS(ID,IZ)'
        DO 825 ID = startid, endid, stepid
          JD = JD + 1
c
c         bug if nizs > maxngs ; Krieger IPP/97
c         DO 820 IZ = 1, NIZS
c
          DO 820 IZ = 1, min(NIZS,maxngs)
            DVALS(JD,IZ) = DEPS(ID,IZ)
            write(99,1000) ID,'',IZ,'',DEPS(ID,IZ)
  820     CONTINUE
          IF (ID.EQ.switchid) JD = JD + 2
  825   CONTINUE
        close(99)
 1000   Format(I4,A1,I4,A1,E12.4)        
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    NIZS,99,DOUTS(1),DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
CX      CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
CX   >    NIZS,99,2.95   ,DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
CX   >    JOB,TITLE,XLAB,YLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.83) THEN
        XLAB = '   Y/A  (M)'
        YLAB = '           '
        JD   = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        CALL RZERO (SUM, 10)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin +1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif
c
        DO 835 ID = startid, endid, stepid
          JD = JD + 1
          DO 830 II = 1, 5
            DVALS(JD,II) = NEROS(ID,II)
            IF (NEROS(ID,II).GT.0.0) THEN
              SUM(II)   = SUM(II)   + NEROS(ID,II) * DWIDS(JD)
            ELSE
              SUM(II+5) = SUM(II+5) + NEROS(ID,II) * DWIDS(JD)
            ENDIF
  830     CONTINUE
          IF (ID.EQ.switchid) JD = JD + 2
  835   CONTINUE
        WRITE(ELABS(1),'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
        WRITE(ELABS(2),'(A,F8.4)')'    PRIMARY REMOVAL  =',SUM(2)+SUM(7)
        WRITE(ELABS(3),'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
        WRITE(ELABS(4),'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
        WRITE(ELABS(5),'(A,2F7.4)')'    NENNL      =',    SUM(5),SUM(10)
        REF = 'NET EROSION'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    5,99,DOUTS(1),DOUTS(NDSIN+1),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    5,99,DOUTS(NDSIN+2),DOUTS(NDS+2),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.84) THEN
        ELABS(1) = 'F   FLUX   '
        ELABS(2) = ' E  ENERGY '
        ELABS(3) = '  Y YIELD  '
        ELABS(4) = 'FY  F * Y  '
        ELABS(5) = '  MRMAX RAN'
        ELABS(6) = 'FRACFY FRAC'
        ELABS(7) = '  T THETA  '
        ELABS(8) = '   BBRATIO '
        XLAB = '   R  (M) '
        YLAB = '          '
        REF  = 'NEUTRAL LAUNCH FUNCTIONS'
        WRITE (IPLOT,9012) NPLOTS,REF
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin + 1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif
c
        DO 840 ID = startid, endid, stepid
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          DVALS(JD,1) = KFLUX (ID)
          DVALS(JD,2) = KENER (ID)
          DVALS(JD,3) = KYIELD(ID)
          DVALS(JD,4) = KFY   (ID)
          DVALS(JD,5) = KRMAX (ID)
          DVALS(JD,6) = KCUM  (ID)
          DVALS(JD,7) = THETAS(ID) * RADDEG
          DVALS(JD,8) = KBFS(IK,IR)
          IF (ID.EQ.switchid) JD = JD + 2
  840   CONTINUE
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    8,99,DOUTS(1),DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.85.OR.IREF.EQ.86.OR.IREF.EQ.87.OR.IREF.EQ.88) THEN
        IZMIN = NIZS + 1
        IZMAX = NIZS + 1
        IF     (IREF.EQ.85) THEN
          M   = 1
          IR  = IRSEP
          REF = 'EXITS ORIG. IONIZ IN SOL,TRAP'
        ELSEIF (IREF.EQ.86) THEN
          M   = 2
          IR  = IRSEP
          REF = 'EXITS ORIG. IONIZ IN MAIN'
        ELSEIF (IREF.EQ.87) THEN
          M   = 3
          IR  = IRSEP - 1
          REF = 'NEUTRALS ENTERING MAIN'
          IZMIN = 0
          IZMAX = 0
        ELSEIF (IREF.EQ.88) THEN
          M   = 3
          IR  = IRSEP - 1
          REF = 'IONS ENTERING MAIN'
        ENDIF
        CALL rzero (kvals, maxnks*maxngs)
        XLAB   = '   S  (M)  '
        YLAB   = '           '
        WRITE (KLAB,'(''CONTOUR K ='',F10.6)') KKS(IR)
        DO 852 IK = 1, NKS(IR)
          KOUTS(IK) = KSS(IK,IR)
          KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          KVALS(IK,1) = (FT*ELIMS(IK,M,0)-FP*ELIMS(IK,M,-1))/KWIDS(IK)
          FACT = 0.0
          DO 850 IZ = -1, NIZS
            KVALS(IK,IZ+3) = ELIMS(IK,M,IZ) / KWIDS(IK)
            IF (IZ.GT.0) FACT = FACT + KVALS(IK,IZ+3)
  850     CONTINUE
          KVALS(IK,NIZS+4) = FACT
  852   CONTINUE
        NPLOTS = NPLOTS + 1
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL DRAW (KOUTS,KWIDS,KVALS(1,IZMIN+3),MAXNKS,NKS(IR),ANLY,
     >   IZMAX-IZMIN+1,99,KOUTS(1),KOUTS(NKS(IR)),
     >   0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >   JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,KLAB,TABLE,2,
     >   2,1.0,0)
      ENDIF
c
C-----------------------------------------------------------------------
c
c     Plot of net erosion/deposition/impurity quantities 
c     for Asdex Upgrade targets
c
C-----------------------------------------------------------------------
c
c
c     Plot of target impurity quantities for Asdex Upgrade target
c     Krieger IPP/98
c
      if (iref.eq.89) then
c
c       erosion/deposition inner AUG target
c
        XLAB   = ' SLUNT(M)'                                            
        YLAB   = '         '                                          
        REF = 'EROS./DEPOS.-INNER AUG TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        CALL RZERO (SUM, 10)
        write(iplot,*)
     >       'Dist. along target, Depos., Gross eros., Net eros., K' 

        JD = 0                                                          
        DO ID = nds,ndsin+1,-1
          JD = JD + 1                                                   
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          tdist(jd) = sqrt((rs(ik,ir)-1.1300)**2+(zs(ik,ir)+0.8200)**2)
          DO II = 1, 5
            DVALS(JD,II) = NEROS(ID,II)
            IF (NEROS(ID,II).GT.0.0) THEN
              SUM(II)   = SUM(II)   + NEROS(ID,II) * DWIDS(JD)
            ELSE
              SUM(II+5) = SUM(II+5) + NEROS(ID,II) * DWIDS(JD)
            ENDIF
          enddo    
          write(iplot,'(1X,1P,E11.4,3(3X,E11.4),3x,0p,f4.2)')
     >      tdist(jd),dvals(jd,1),dvals(jd,3),dvals(jd,4),kks(ir)
        ENDDO
c
        WRITE(ELABS(1),'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
        WRITE(ELABS(2),'(A,F8.4)')'    PRIMARY REMOVAL  =',SUM(2)+SUM(7)
        WRITE(ELABS(3),'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
        WRITE(ELABS(4),'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
        WRITE(ELABS(5),'(A,2F7.4)')'    NENNL      =',    SUM(5),SUM(10)
        WRITE(iplot,'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
        WRITE(iplot,'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
        WRITE(iplot,'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
c
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,nds-ndsin,ANLY,                 
     >    5,99,tdist(1),tdist(nds-ndsin),-HI,HI,IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)      
c
c       erosion/deposition outer AUG target
c
        XLAB   = ' SRUNT(M)'                                            
        YLAB   = '           '                                          
        REF = 'EROS./DEPOS.-OUTER AUG TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        CALL RZERO (SUM, 10)
        write(iplot,*)
     >       'Dist. along target, Depos., Gross eros., Net eros., K' 

        JD = 0                                                          
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)                               
        DO ID = ndsin,1,-1
          JD = JD + 1                                                   
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          tdist(jd) = sqrt((rs(ik,ir)-1.4820)**2+(zs(ik,ir)+0.9526)**2)
          DO II = 1, 5
            DVALS(JD,II) = NEROS(ID,II)
            IF (NEROS(ID,II).GT.0.0) THEN
              SUM(II)   = SUM(II)   + NEROS(ID,II) * DWIDS(JD)
            ELSE
              SUM(II+5) = SUM(II+5) + NEROS(ID,II) * DWIDS(JD)
            ENDIF
          enddo    
          write(iplot,'(1X,1P,E11.4,3(3X,E11.4),3x,0p,f4.2)')
     >      tdist(jd),dvals(jd,1),dvals(jd,3),dvals(jd,4),kks(ir)
        ENDDO
c
        WRITE(ELABS(1),'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
        WRITE(ELABS(2),'(A,F8.4)')'    PRIMARY REMOVAL  =',SUM(2)+SUM(7)
        WRITE(ELABS(3),'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
        WRITE(ELABS(4),'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
        WRITE(ELABS(5),'(A,2F7.4)')'    NENNL      =',    SUM(5),SUM(10)
        WRITE(iplot,'(A,F8.4)')'    TOTAL DEPOSITION =',SUM(1)+SUM(6)
        WRITE(iplot,'(A,F8.4)')'    TOTAL REMOVAL    =',SUM(3)+SUM(8)
        WRITE(iplot,'(A,2F7.4)')'    NET EROSION=',     SUM(4),SUM(9)
c
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,ndsin,ANLY,                 
     >    5,99,tdist(1),tdist(ndsin),-HI,HI,IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,2,1.0,0)      
c
c       plot deposition along inner AUG target, Krieger IPP/98
c
        XLAB = ' SLUNT(M)'                                            
        YLAB = '         '
        REF  = 'AUG DEPOSITION - INNER TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        write(iplot,*)
     >       'Dist. along target, Deposition iz=1,8' 
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        JD   = 0
        DO ID = nds,ndsin+1,-1
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          JD = JD + 1
          tdist(jd) = sqrt((rs(ik,ir)-1.1300)**2+(zs(ik,ir)+0.8200)**2)
          DO IZ = 1,min(8,nizs)
            DVALS(JD,IZ) = DEPS(ID,IZ)
          END DO
          write(iplot,'(1X,1P,E10.3,8(2X,E10.3))')
     >         tdist(jd),(dvals(jd,iz),iz=1,min(nizs,maxngs,8))
        END DO
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,NDS-ndsin,ANLY,
     >    min(8,nizs),99,tdist(1),tdist(nds-ndsin),0.0,HI,
     >    IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,1,2,1.0,0)
c
c       plot deposition along outer AUG target, Krieger IPP/98
c
        XLAB = ' SRUNT(M)'                                            
        YLAB = '         '
        REF  = 'AUG DEPOSITION - OUTER TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        write(iplot,*)
     >       'Dist. along target, Deposition iz=1,8' 
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        JD   = 0
        DO ID = ndsin,1,-1
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          JD = JD + 1
          tdist(jd) = sqrt((rs(ik,ir)-1.4820)**2+(zs(ik,ir)+0.9526)**2)
          DO IZ = 1,min(8,nizs)
            DVALS(JD,IZ) = DEPS(ID,IZ)
          END DO
          write(iplot,'(1X,1P,E10.3,8(2X,E10.3))')
     >         tdist(jd),(dvals(jd,iz),iz=1,min(nizs,maxngs,8))
        END DO
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,ndsin,ANLY,
     >    min(8,nizs),99,tdist(1),tdist(ndsin),0.0,HI,
     >    IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,1,2,1.0,0)

c       first calculate total densities/concentrations

c       We want here the "real" density, not per injected
c       particle -> multiply by absfac

        call rzero(plastmp,maxnks*maxnrs)
        do iz = 0,nizs
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 if (knbs(ik,ir).ne.0.0) then 
                    plastmp(ik,ir) = plastmp(ik,ir) 
     >                    + sdlims(ik,ir,iz)*absfac
                 endif
              end do 
           end do
        end do 

c       We want here the "real" concentration, not per injected
c       particle -> multiply by absfac

        call rzero(ktmp,maxnks*maxnrs)
        do iz = 0,nizs
           do ir = 1,nrs
              do ik = 1,nks(ir)
                 if (knbs(ik,ir).ne.0.0) then 
                    ktmp(ik,ir) = ktmp(ik,ir) 
     >                    + sdlims(ik,ir,iz)/knbs(ik,ir)*absfac
                 endif
              end do 
           end do
        end do 
c
c       impurity density/concentration along inner AUG target
c
        XLAB   = ' SLUNT(M)'                                            
        YLAB   = '         '                                          
        REF = 'Nz/Fz - INNER AUG TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        write(iplot,*)
     >       'Dist. along target, Nz, Fz, K' 
        JD = 0                                                          
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)                               
        ELABS(1) = '      Nz '
        ELABS(2) = '      Fz '
        DO ID = nds,ndsin+1,-1
          JD = JD + 1                                                   
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          tdist(jd) = sqrt((rs(ik,ir)-1.1300)**2+(zs(ik,ir)+0.8200)**2)
          DVALS(JD,1) = plastmp(ik,ir)
          DVALS(JD,2) = ktmp(ik,ir)
          write(iplot,'(1X,1P,E11.4,2(3X,E11.4),3x,0p,f4.2)')
     >      tdist(jd),dvals(jd,1),dvals(jd,2),kks(ir)
        ENDDO
c
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,nds-ndsin,ANLY,                 
     >    2,99,tdist(1),tdist(nds-ndsin),0.0,HI,IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ELABS(1),REF,NVIEW,PLANE,TABLE,2,2,1.0,0)      
c
c       impurity density along outer AUG target, Krieger, IPP/98
c
        XLAB   = ' SRUNT(M)'                                            
        YLAB   = '           '                                          
        REF = 'Nz/Fz - OUTER AUG TARGET'
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
        write(iplot,*)
     >       'Dist. along target, Nz, Fz, K' 
        JD = 0                                                          
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)                               
        ELABS(1) = '      Nz '
        ELABS(2) = '      Fz '
        DO ID = ndsin,1,-1
          JD = JD + 1                                                   
          IK = IKDS(ID)                                                 
          IR = IRDS(ID)                                                 
          tdist(jd) = sqrt((rs(ik,ir)-1.4820)**2+(zs(ik,ir)+0.9526)**2)
          DVALS(JD,1) = plastmp(ik,ir)
          DVALS(JD,2) = ktmp(ik,ir)
          write(iplot,'(1X,1P,E11.4,2(3X,E11.4),3x,0p,f4.2)')
     >      tdist(jd),dvals(jd,1),dvals(jd,2),kks(ir)
        ENDDO
c
        CALL DRAW (tdist,DWIDS,DVALS,MAXNDS+2,ndsin,ANLY,                 
     >    2,99,tdist(1),tdist(ndsin),0.0,HI,IGNORS,ITEC,AVS,NAVS,       
     >    JOB,TITLE,XLAB,YLAB,ELABS(1),REF,NVIEW,PLANE,TABLE,2,2,1.0,0)      
c
      endif
c
C-----------------------------------------------------------------------
c
c     Plot of particles leaving grid towards wall  Krieger IPP/97
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.90) THEN
c
c       SOL wall ring
c

         write(6,*) 'WALLS Integration:'
         tmpsum = 0.0
        do ir = 1,nrs
           do ik = 1,nks(ir)
              tmpsumiz = 0.0
              do iz = 1,nizs
                 tmpsumiz = tmpsumiz + walls(ik,ir,iz)
              end do

              if (tmpsumiz.ne.walls(ik,ir,maxizs+1)) then 
                 write(6,'(a,2i8,3(1x,g18.8))') 
     >        'WALLS SUM MISMATCH:',ik,ir,tmpsumiz,walls(ik,ir,maxizs+1)
              endif 
              tmpsum(ir) = tmpsum(ir) + tmpsumiz
              tottmpsum = tottmpsum + tmpsum(ir)
          end do 
       end do

       do ir = 1,nrs
          write(6,'(a,i8,1x,g18.8)') 'WALLS Int:',ir,tmpsum(ir)
       end do
       write(6,'(a,i8,1x,g18.8)') 'WALLS Int:',ir,tottmpsum



        ELABS(1) = 'N   Nexits '
        XLAB = '   POLOIDAL DIST (M)'
        YLAB = '           '
        IR = IRWALL
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        

c
c     
c       Get total number of particles launched
c
        totsrc = (targsrc(3,4) + wallsrc(5,3))
        if (totsrc.le.1.0) totsrc = 1.0
c
        irref = irwall

        write(iplot,'(a,1x,i8,2(1x,g18.8))') 
     >             ' Wall exits - SOL wall ring:', irref,totsrc,absfac
        write(iplot,'(a)') ' IK  L-Pol  N   WID   RS  ZS'//
     >                     '   S   P  WALLS   WALLS-FLUX'
c
c       jdemod: Note that the axis data stored in IR=IRWALL may not be 
c               correct since the otuermost ring is used for accounting and
c               may not have polygon data or the polygons may be zero volume.
c
c               Use the axis data from one ring inward.
c
c
        DO IKref2 = 1, NKS(IrRef)
c
c         Overwrite IR in loop 
c
          ir = irins(ikref2,irref)
          ik = ikins(ikref2,irref)
c
          KOUTS(IK) = KPS(IK,IR)
          KWIDS(IK) = 0.0
          IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
          IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
c
c         jdemod: walls data must still reference actual wall ring
c
          KVALS(IK,1) = WALLS(IKref2,IRref,MAXIZS+1) / kwids(ik)
c
c         jdemod - print the actual particles out
c
          write(iplot,'(1x,i6,1x,f12.5,2x,10(1p,e18.8))')
     >        ik, kouts(ik), kvals(ik,1),kwids(ik),
     >        rs(ik,ir),zs(ik,ir),kss(ik,ir),kps(ik,ir), 
     >        walls(ikref2,irref,maxizs+1),
     >        kvals(ik,1)/totsrc*absfac
        END DO
        

*       KVALS(1,1) = 0.0
*       KVALS(NKS(IR),1) = 0.0
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IRref),ANLY,
     >    1,1,KOUTS(1),KOUTS(NKS(IRref)),0.,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
c
c       Trap wall ring
c
c       Krieger IPP/09 - do this plot only if trap exists
        if (irtrap.le.nrs) then
          ELABS(1) = 'N   Nexits '
          XLAB = '   POLOIDAL DIST (M)'
          YLAB = '           '
          IRref = IRTRAP
          WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (kvals, maxnks*maxngs)
c
c
          write(iplot,'(a,1x,i8,1x,2(1x,g18.8))')
     >        ' Wall exits - PFZ wall ring:', irref,totsrc,absfac
          write(iplot,'(a)') ' IK  L-Pol  N   WID   RS  ZS'//
     >                     '   S   P  WALLS   WALLS-FLUX'
          DO IKref2 =  1, NKS(IRref)

          ir = irouts(ikref2,irref)
          ik = ikouts(ikref2,irref)
c

            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
            KVALS(IK,1) = WALLS(IKref2,IRref,MAXIZS+1) / kwids(ik)

            write(iplot,'(1x,i6,1x,f12.5,2x,10(1p,e18.8))')
     >        ik, kouts(ik), kvals(ik,1),kwids(ik),
     >        rs(ik,ir),zs(ik,ir),kss(ik,ir),kps(ik,ir), 
     >                    walls(ikref2,irref,maxizs+1),
     >        kvals(ik,1)/totsrc*absfac

          END DO
*         KVALS(1,1) = 0.0
*         KVALS(NKS(IR),1) = 0.0
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IRref),ANLY,
     >      1,1,KOUTS(1),KOUTS(NKS(IRref)),0.,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
        endif
c
c       Neutral exits on wall elements
c
        ELABS(1) = 'N   Nexits '
        XLAB = '   # of Wall Element'
        YLAB = '           '
        IR = IRTRAP
        WRITE (REF,'(''ALONG NEUTRAL WALL'')')
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
		write(iplot,*) 'Neutral wall exits'
		write(iplot,'(a)') ' IN  N'
		IK=0
        DO IN = 1, MAXPTS
          if (wallsn(in).gt.0.0) then
            ik = ik + 1
            KOUTS(IK) = IN
            KWIDS(IK) = 1.0
            KVALS(IK,1) = WALLSN(IN)
            write(iplot,'(i3,1x,e12.5)') in, wallsn(in)
          endif
        END DO
*       CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,IK,ANLY,
*    >    1,1,1.,KOUTS(IK),0.,HI,IGNORS,ITEC,AVS,NAVS,
*    >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
C
C-----------------------------------------------------------------------
C     Z EFFECTIVES ETC CONTOUR PLOTS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.91.OR.IREF.EQ.92.OR.IREF.EQ.93) THEN
c
c test for recycling impurities (fudge factor)
c
        fact = 2.9e-3
        do ir = 1, nrs
          do ik = 1, nks(ir)
            cvalsa(ik,ir) = 1.0
            do iz = 0, nizs
              cvalsa(ik,ir) = cvalsa(ik,ir)
     >          + iz*(iz-1)*fact*absfac*sdlims(ik,ir,iz)
     >               /(rizb*knbs(ik,ir))
            enddo
          enddo
        enddo
        call prdmat(cvalsa,maxnks,nks(irsep),nrs,irsep,irwall,
     >              ikto+1,ikti-1,6,'zeff:')
c
c end of test
c
        ELABS(1) = '    NIE    '
        ELABS(2) = '    ZB.NBT '
        ELABS(3) = '    ZEFF   '
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
        NGS    = CNGS
        II     = IREF - 90
        REF    = ELABS(II)(5:11)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,ZEFFS,II,3,3,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
c
C-----------------------------------------------------------------------
c
c     ZEFFS ALONG A CONTOUR LINE
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.94.or.IREF.eq.95) THEN
        write (elabs(1),'(''ZeffZeffs'')')
        IF (IREF.EQ.94) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'IMPURITY ZEFFS'
        READ (GRAPH(38:44),'(I4)') IR
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        DO 940 IK = 1, NKS(IR)
          IF (IREF.EQ.94) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = zeffs(ik,ir,3)
c
 940    CONTINUE
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*KOUTS(NKS(IR)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
      ENDIF
c
C-----------------------------------------------------------------------
c
c     LEAKAGE PLOTS
c
c     Plot type 96: Plot of ion leakage as a function of charge
c     state and distance S along the field line. The S-values are
c     specified as input to DIVIMP.
c
C-----------------------------------------------------------------------
c
      IF (IREF.EQ.96) THEN
        XLAB = '   S  (M)'
        YLAB   = 'IMPURITY LEAKAGE'
        call rzero (cwids,maxpts)
        CALL DRAW (cleaks,cWIDS,cleakn,MAXPTS,cleaksn,ANLY,
     >    nizs+1,99,cleaks(1),cleaks(cleaksn),-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,
     >    JOB,TITLE,XLAB,YLAB,ZLABS(1),REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
      ENDIF

C
C-----------------------------------------------------------------------
C     PIN MOMENTUM SOURCE FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.97) THEN
        ELABS(1) = '    MOMENTUM SRC'
        IF (IREF.EQ.97) THEN
          XLAB = '   S  (M)'
        ELSE
          XLAB = '   POLOIDAL DIST (M)'
        ENDIF
        YLAB   = 'PIN MOM SRC (PA/M)'
        READ (GRAPH(38:44),'(I4)') IR
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''ALONG RING '',I3,'', K'',F7.4)') IR,KKS(IR)
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL rzero (kvals, maxnks*maxngs)
        PLTMAX = -HI
        PLTMIN = HI
        DO IK = 1, NKS(IR)
          IF (IREF.EQ.97) THEN
            KOUTS(IK) = KSS(IK,IR)
            KWIDS(IK) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))
          ELSE
            KOUTS(IK) = KPS(IK,IR)
            KWIDS(IK) = 0.0
            IF (IK.GT.1) KWIDS(IK) = 0.5 * (KPS(IK,IR)-KPS(IK-1,IR))
            IF (IK.LT.NKS(IR)) KWIDS(IK) = KWIDS(IK) +
     >                               0.5 * (KPS(IK+1,IR)-KPS(IK,IR))
          ENDIF
c
          KVALS(IK,1) = PINMP(IK,IR)
c
          PLTMAX = MAX(PLTMAX,KVALS(IK,1))
          PLTMIN = MIN(PLTMIN,KVALS(IK,1))
        ENDDO
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),KOUTS(NKS(IR)),PLTMIN,PLTMAX,
     >    IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,1,6,1.0,0)
c
        if (clsup.eq.1) then
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,KOUTS(1),mgst*kouts(nks(ir)),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR),ANLY,
     >    1,99,mgnd*KOUTS(NKS(IR)),KOUTS(NKS(IR)),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
        endif
c
       ENDIF

C
C-----------------------------------------------------------------------
C     TARGET ELEMENT LEAKAGE PLOTS - TARGSOURCE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.98) THEN
        ELABS(1) = 'TAU Tau Leak'
        ELABS(2) = 'DN  DelatN'
        ELABS(3) = '%-LK%-Leakage'
        REF  = 'TARGET TO CORE LEAKAGE'
        WRITE (IPLOT,9012) NPLOTS,REF
c
        JD = 0
        CALL RZERO (DVALS, (MAXNDS+2)*MAXNGS)
c
        if (rp(1).gt.rp(nds)) then
          startid = nds
          endid   = 1
          stepid  = -1
          switchid = ndsin + 1
        else
          startid = 1
          endid   = nds
          stepid = 1
          switchid = ndsin
        endif
c
        ntotal = 0.0
c
        DO ID = startid, endid, stepid
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
          do ir = 1,maxnrs
             do in = 1,3
                dvals(jd,1) = dvals(jd,1) + wtsource(id,ir,1,in)
                dvals(jd,2) = dvals(jd,2) + wtsource(id,ir,3,in)
             end do
          end do
c
          if (dvals(jd,1) .gt.0.0) then
             dvals(jd,3) = dvals(jd,2)/dvals(jd,1)
          else
             dvals(jd,3) = 0.0
          endif
c
          ntotal = ntotal + dvals(jd,2)
c
          IF (ID.EQ.switchid) JD = JD + 2
c
        end do
c
c
        DO ID = startid, endid, stepid
          JD = JD + 1
          IK = IKDS(ID)
          IR = IRDS(ID)
c
          if (ntotal.gt.0.0) then
             dvals(jd,2) = dvals(jd,2)/ntotal * core_content
          else
             dvals(jd,2) = 0.0
          endif
c
          if (dvals(jd,1).gt.0.0) then
             dvals(jd,1) =  dvals(jd,2)/(dvals(jd,1) * facta(1))
          else
             dvals(jd,1) = 0.0
          endif
c
          IF (ID.EQ.switchid) JD = JD + 2
c
        end do
c
        CALL DRAW (DOUTS,DWIDS,DVALS,MAXNDS+2,NDS+2,ANLY,
     >    3,99,DOUTS(1),DOUTS(NDS+2),0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,2,1.0,0)
      ENDIF
c


      return 
c
c     Fortmat statements
c

 9012 FORMAT(1X,'PLOT',I3,4X,A)


      end

