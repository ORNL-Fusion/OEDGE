c
c ======================================================================
c


      SUBROUTINE Plot970(job,graph,ref,title,iopt)

      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'slout'



      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      REAL CalcPressure,GetCs

      integer cngs,cntropt

      integer iplot
c
c      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,POINT
c
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*(*),TITLE2*80,JOB*72,GRAPH*80,GRAPH1*80,graph6*128

      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB
c    >            ,XPOINT
      CHARACTER*72 YLAB
      CHARACTER*36 REF
c    >            ,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,PLABS(-2:MAXPLRP),KLAB

      CHARACTER*128 elabs(MAXNGS)

      REAL xrange1,xrange2

      integer ik,ir,i1,i2

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,NPLOTS,M,ID,JR

c
c      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
c     .        in2,xtype,ytype,btype,id1,id2
c
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact
c
c      REAL tauii,pii
c
c      INTEGER nenum,tenum,opt_const,plot_mode(30),iter1,iter2,
c     .        xaxis,ring,mode,inorm(MAXNGS)
c
c
c      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
c     .        max1,max2,ynorm(MAXNGS)
c
c      REAL cdata(MAXNKS,MAXNRS)
c
c      REAL    XXMIN,XXMAX,YYMIN,YYMAX,SUM(10),VAL,VALTOT
c
c      CHARACTER*72 SMOOTH
c
c      integer icntr
c      integer nconts,nclev
c      real conts(maxpts),clev(maxpts)
c
c      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
c      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)
c
c
c      CHARACTER*128 dataline,cdum1,cdum2
c
      character*128 cdum1,cdum2

      INTEGER   size
      REAL      xpos,ypos
      CHARACTER dummy*5000,caption*5000

c...980:
c      INTEGER NUMTHE,AVPTS,ATYPE
c      INTEGER numth2,numthe1(MAXNGS)
c      INTEGER IGNORS(MAXNGS),ITEC,NAVS
c      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS),
c     .     touts1(MAXTHE,MAXNGS)
c      REAL TOUTS2(MAXTHE),TWIDS2(MAXTHE),TVALS2(MAXTHE,MAXNGS),
c     .     dum1(MAXTHE),dum2(MAXTHE),dum3(MAXTHE),dum4(MAXTHE),
c     .     dum5(MAXTHE),dum6(MAXTHE),dum7(MAXTHE),dum8(MAXTHE),
c     .     den1,teb1
c      REAL ZOBS,ROBS,DRAD,DTHE,THEMIN,THEMAX,themin_start
c      real theres,mfact
c      REAL WLNGTH
c      real    zsuma(maxizs),cvalsa(maxnks,maxnrs)
c      REAL LEVEL,AVS(0:100),VMIN,VMAX
c      CHARACTER ADASID*80,PECTITLE*120,PLABAD*36
c      CHARACTER XFESYM*2
c      character adasex*3
c      integer   adasyr
c      CHARACTER ADASGR*8,ADASTY*80
c      character datatitle*128
C
C     Second sets of ADAS data for RATIO plots
C
c      character graph5*80,adasid2*80,adasex2*3
c      integer   adasyr2,isele2,iselr2,iselx2,iseld2
c      integer   iz_state,z_atom,iz_state2,z_atom2
c      INTEGER plot
c      character graph2*80,graph3*80,
c
      character graph4*80
c
c      INTEGER ISELE,ISELR,ISELX,iseld,iseldef
c      INTEGER IADAS,NPAIRS,IRCODE,IKK,LEN,LENSTR
c      INTEGER IK,II,IT,LT,UT,IREF,IYMIN,IYMAX,IR,JD,cnt
c      INTEGER IZMIN,IZMAX,IW,LW,UW
c      REAL PLRPAD(MAXNKS,MAXNRS)
c      LOGICAL plotr
c      real zadj
c      REAL peak1,peak2,array(MAXNKS,MAXNRS)
c      INTEGER plotcode,line,iopt1,iopt2,iopt3
c      CHARACTER cname*256,resdir*256,cmnd*256
       LOGICAL status
c     >       ,oldraw
c
c...970:
      CHARACTER*50   tag
      CHARACTER*2048 ldata
      INTEGER       ifp,xmode,xbnd1,xbnd2,ymode,step,count,match,ir1,
     .              lstep,liter,lcount,in1,eircount,leircount,
     .              idum1,idum2,idum3,idum4,idum5,probe,iliin
      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,fact,
     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS),
     .              atm(6),mol(6)
c
c     Array for extra comments 
c
      character*20 extra_comments(1)
c


      CHARACTER*30 gaugename(110),surmat(20)

c      gaugename(001) = 'Midplane gauge'
c      gaugename(002) = 'Divertor gas box gauge'

      WRITE(0,*) '970:'

      gaugename(101) = 'PBF1 pump chamber'
      gaugename(102) = 'PBF2'
      gaugename(103) = 'PBF3'
      gaugename(104) = 'PV1  private zone'
      gaugename(105) = 'PR2'
      gaugename(106) = 'VPLOWS'
      gaugename(107) = 'PCM105BAF'
      gaugename(108) = 'PCM240TOR'

        ref = graph( 7:80)

        iopt_ghost = 1

        CALL RZero(mvals,MAXNKS*MAXPLTS*MAXNGS)
        CALL RZero(mouts,MAXNKS*MAXPLTS)
        CALL RZero(mwids,MAXNKS*MAXPLTS)
c
c       Relaxation history plots:
c
        READ(5,*) cdum1,xmode,xbnd1,xbnd2,ymode

c        WRITE(0,*) '970: YMODE = ',ymode
c
c       Load plots per page info:
c
        CALL RDG4(graph4,ngrm,nplts,ringnos,maxplts,pltfact,ierr)

        DO ip = 1, nplts
          IF ((ymode.EQ.14.OR.ymode.EQ.15.OR.ymode.EQ.18).AND.
     .        ringnos(ip).LE.0.AND.ringnos(ip).NE.-99) THEN
c...        Use the gauge name instead of a number index.  This 
c           will only work if the EIRPGDAT array is being assigned,
c           which shouldn't be a problem:
            pltlabs(ip) = 
     .        gaugename(NINT(eirpgdat(eirnpgdat+ringnos(ip),1)))
          ELSE
            IF     (ringnos(ip).LT.1000) THEN
              WRITE (pltlabs(ip),'(I3)') ringnos(ip)
            ELSEIF (ringnos(ip).LT.10000) THEN 
              WRITE (pltlabs(ip),'(I4)') ringnos(ip)
            ELSEIF (ringnos(ip).LT.100000) THEN
              WRITE (pltlabs(ip),'(I5)') ringnos(ip)
            ENDIF
          ENDIF
        ENDDO

        sctype = iopt
        IF (sctype.LT.1.OR.sctype.GT.4) sctype = 1

        IF     (xmode.EQ.1) THEN
          xlab     = 'iteration'
        ELSEIF (xmode.EQ.2) THEN
          xlab     = 'step'
        ELSEIF (xmode.EQ.3) THEN
          xlab     = 'EIRENE self-iteration'
        ELSE
          CALL ER('970','Invalid XMODE value',*99)
        ENDIF

        IF     (ymode.EQ.1) THEN
          ylab     = ' '
          elabs(1) = '    Te1'
          elabs(2) = '    Ti1'
          elabs(3) = '    Te2'
          elabs(4) = '    Ti2'
          ngs      = 4
        ELSEIF (ymode.EQ.2) THEN
          ylab     = ' '
          elabs(1) = '    RMSTe'
          elabs(2) = '    RMSTi'
          elabs(3) = '    RMSne'
c          elabs(4) = '    PeakTe'
c          elabs(5) = '    PeakTi'
c          elabs(6) = '    Peakne'
          ngs      = 3
        ELSEIF (ymode.EQ.3) THEN
          ylab     = ' '
          elabs(1) = '    QeMlo'
          elabs(2) = '    QeMhi'
          ngs      = 2
        ELSEIF (ymode.EQ.4) THEN
          ylab     = ' '
          elabs(1) = '    ne1'
          elabs(2) = '    ne2'
          ngs      = 2
        ELSEIF (ymode.EQ.5) THEN
          ylab     = ' '
          elabs(1) = '    PeiMlo'
          elabs(2) = '    PeiMhi'
          ngs      = 2
        ELSEIF (ymode.EQ.6) THEN
          ylab     = ' '
          elabs(1) = '    Te1'
          elabs(2) = '    Te2'
          elabs(3) = '    Ti1'
          elabs(4) = '    Ti2'
          ngs      = 4
        ELSEIF (ymode.EQ.7) THEN
c...      Momentum peak:
          ylab     = ' '
          elabs(1) = '    Ne1'
          elabs(2) = '    Ne2'
          ngs      = 2
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
        ELSEIF (ymode.EQ.8) THEN
          ylab     = ' '
          elabs(1) = '    Mach1'
          elabs(2) = '    Mach2'
          ngs      = 2
        ELSEIF (ymode.GE.9.AND.ymode.LE.12) THEN
          ylab     = ' '

c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 5
          plottype(5) = 6

          elabs(1) = '    1'
          elabs(2) = '    2'
          elabs(3) = '    3'
          elabs(4) = '    4'
          elabs(5) = '    5'
          ngs      = 5
        ELSEIF (ymode.EQ.13) THEN
          ylab     = ' '
          elabs(1) = '    1'
          elabs(2) = '    2'
          ngs      = 2
        ELSEIF (ymode.EQ.14) THEN
          ylab     = ' T (eV)'
          elabs(1) = '    D'
          elabs(2) = '    Davg'
          elabs(3) = '    D2'
          elabs(4) = '    D2avg'
          ngs      = 4
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 6
        ELSEIF (ymode.EQ.15) THEN
          ylab     = ' n (m-3)'
          elabs(1) = '    D'
          elabs(2) = '    Davg'
          elabs(3) = '    D2'
          elabs(4) = '    D2avg'
          ngs      = 4
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 6
        ELSEIF (ymode.EQ.16) THEN
c...      Momentum loss:
          ylab     = ' '
          elabs(1) = '    1'
          elabs(2) = '    2'
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          ngs      = 2
        ELSEIF (ymode.EQ.17) THEN
c...      Pressure gauge data:
          ylab     = 'p (mTorr)'
          elabs(1) = '    D'
          elabs(2) = '    D2'
          elabs(3) = '    Total'
          elabs(4) = '    Exp'
          ngs      = 4
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 6
        ELSEIF (ymode.EQ.18) THEN
c...      Pressure gauge data (from additional cell data):
          ylab     = 'p (mTorr)'
          elabs(1) = '    D'
          elabs(2) = '    Davg'
          elabs(3) = '    D2'
          elabs(4) = '    D2avg'
          ngs      = 4
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 6
        ELSEIF (ymode.EQ.19) THEN
c...      Additional cell T_D2 from ... and ...:
          ylab     = ' T (eV)'
          elabs(1) = '    1'
          elabs(2) = '    2'
          ngs      = 2
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
        ELSEIF (ymode.EQ.20) THEN
c...      Particle flux through EIRENE additional surface:
c          ylab     = ' atomic flux / source'
          ylab     = ' flux (atoms s-1)'
          elabs(1) = '    +atm'
          elabs(2) = '    -atm'
          elabs(3) = '    +mol'
          elabs(4) = '    -mol'
          elabs(5) = '    net'
          ngs      = 5
c...      Turn on coloured lines on plots:
          slopt2 = 2
          plottype(1) = 2
          plottype(2) = 3
          plottype(3) = 4
          plottype(4) = 6
          plottype(5) = 7
c          WRITE(0,*) '970 20: DIVIDING BY 5.0E+21'
        ELSE
          CALL ER('970','Invalid YMODE value',*99)
        ENDIF


        ifp = 79
        OPEN(UNIT=ifp,ACCESS='SEQUENTIAL',STATUS='OLD',ERR=99)

        DO ip = 1, nplts

          REWIND(ifp)

          IF (ymode.EQ.14.OR.ymode.EQ.15.OR.ymode.EQ.18.OR.
     .        ymode.EQ.19) THEN
c...        This is a little nasty, but since the additional cells are stored
c           with the large vacuum cell and the pressure gauge cells, and
c           are better referred to as additional cell 1,2,etc. for the user (rather
c           than 1+1+1...1+2+eirnpgdat,etc., especially since changing eirnpgdat
c           from case to case would mess up the indexing in the OUT input file) I
c           have added the required index shift here.  Pressure gauge data can
c           still be accessed as 1-x, where x is the pressure gauge data index in
c           reverse, but it is still preferred that pressure gauge data be
c           accessed via plot 17:
            ir = ringnos(ip) + 1 + eirnpgdat
          ELSE
            ir = ringnos(ip)
          ENDIF
          in     =  0
          lstep  = -1
          liter  = -1
          lcount = -1

          DO WHILE (1.EQ.1)
            READ(ifp,'(A2048)',END=2500) ldata

            READ(ldata,*) tag

            IF     (tag(1:10).EQ.'SYMMETRY  '.AND.
     .              (ymode.EQ.1.OR.ymode.EQ.4)) THEN

              IF     (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*) tag,step,iter,match,ir1,
     .                        idum1,idum2,idum3,
     .                        te1,ti1,ne1,norm,
     .                        rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7
                count = 0
                te2   = te1
                ti2   = ti1
                ne2   = ne1
              ELSEIF (tag(11:14).EQ.'1.01') THEN
                READ(ldata,*) tag,step,iter,count,match,ir1,
     .                        idum1,idum2,idum3,
     .                        te1,ti1,ne1,norm,
     .                        rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7
                te2 = te1
                ti2 = ti1
                ne2 = ne1
              ELSEIF (tag(11:14).EQ.'1.02') THEN
                READ(ldata,*) tag,step,iter,count,match,ir1,
     .                        idum1,idum2,idum3,
     .                        te1,ti1,ne1,te2,ti2,ne2,norm,
     .                        rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7
              ELSEIF (tag(11:14).EQ.'1.03') THEN
                READ(ldata,*) tag,step,iter,count,match,ir1,
     .                        idum1,idum2,idum3,
     .                        rdum1,rdum2,rdum3,rdum4,rdum5
                READ(ifp  ,*) tag,te1,ti1,ne1,te2,ti2,ne2
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF

              IF (ir.EQ.ir1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  IF     (ymode.EQ.1) THEN
                    mvals(in,ip,1) = te1
                    mvals(in,ip,2) = ti1
                    mvals(in,ip,3) = te2
                    mvals(in,ip,4) = ti2
                  ELSEIF (ymode.EQ.4) THEN
                    mvals(in,ip,1) = ne1
                    mvals(in,ip,2) = ne2
                  ENDIF
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF

            ELSEIF (tag(1:10).EQ.'PLASMAVAR '.AND.ymode.EQ.2) THEN

              IF     (tag(11:14).EQ.'1.01') THEN
                READ(ldata,*) tag,step,iter,count,ir1,rdum1,
     .                        (radum1(i1),radum2(i1),i1=1,3)
              ELSE
                CALL ER('970','RMS y-data version unsupported',*99)
              ENDIF

              IF (ir.EQ.ir1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  mvals(in,ip,1) = radum1(1)
                  mvals(in,ip,2) = radum1(2)
                  mvals(in,ip,3) = radum1(3)

c                  mvals(in,ip,4) = radum2(1)
c                  mvals(in,ip,5) = radum2(2)
c                  mvals(in,ip,6) = radum3(3)
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF

            ELSEIF (tag(1:10).EQ.'QeMUL     '.AND.ymode.EQ.3) THEN

              IF     (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*)
     .            tag,step,iter,count,ir1,radum1(1),radum2(1)
              ELSE
                CALL ER('970','QeM y-data version unsupported',*99)
              ENDIF

              IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

              IF (ir.EQ.ir1) THEN
                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  mvals(in,ip,1) = radum1(1)
                  mvals(in,ip,2) = radum2(1)
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF


            ELSEIF (tag(1:10).EQ.'PeiMUL    '.AND.ymode.EQ.5) THEN

              IF (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*)
     .            tag,step,iter,count,ir1,radum1(1),radum2(1)
              ELSE
                CALL ER('970','PeiM y-data version unsupported',*99)
              ENDIF

              IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

              IF (ir.EQ.ir1) THEN
                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  mvals(in,ip,1) = radum1(1)
                  mvals(in,ip,2) = radum2(1)
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF

            ELSEIF (tag(1:10).EQ.'TARGET    '.AND.
     .              (ymode.EQ.6.OR.ymode.EQ.7.OR.ymode.EQ.8)) THEN
c
c 6 - temperatures
c 7 - densities
c 8 - Mach number at target
c
              IF     (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*) tag,step,iter,count,ir1,
     .                        te1,te2,ti1,ti2,ne1,ne2,
     .                        rdum1,rdum2
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF

              IF (ir.EQ.ir1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  IF     (ymode.EQ.6) THEN
                    mvals(in,ip,1) = te1
                    mvals(in,ip,2) = ti1
                    mvals(in,ip,3) = te2
                    mvals(in,ip,4) = ti2
                  ELSEIF (ymode.EQ.7) THEN
                    mvals(in,ip,1) = ne1
                    mvals(in,ip,2) = ne2
                  ELSEIF (ymode.EQ.8) THEN
                    mvals(in,ip,1) = rdum1
                    mvals(in,ip,2) = rdum2
                  ENDIF
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF

            ELSEIF ((tag(1:10).EQ.'FLUXES    '.OR.
     .               tag(1:10).EQ.'ION       '.OR.
     .               tag(1:10).EQ.'REC       '.OR.
     .               tag(1:10).EQ.'CFP       '.OR.
     .               tag(1:10).EQ.'MOM LOSS  ').AND.
     .              (ymode.EQ. 9.OR.ymode.EQ.10.OR.ymode.EQ.11.OR.
     .               ymode.EQ.12.OR.ymode.EQ.16)) THEN
c
c  9 - fluxes
c 10 - integrated ionisation
c 11 - integrated recombination
c 12 - integrated cross-field particle source
c 16 - integrated momentum loss - CX only?
c
              IF     (tag(11:14).EQ.'1.00') THEN
                IF (tag(1:10).EQ.'MOM LOSS  ') THEN
                  READ(ldata,*) tag,step,iter,count,ir1,
     .                          rdum1,rdum2,rdum3,rdum4
                ELSE
                  READ(ldata,*) tag,step,iter,count,ir1,
     .                          rdum1,rdum2,rdum3,rdum4,rdum5
                ENDIF
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF

              IF (ir.EQ.ir1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                 IF   ((tag(1:10).EQ.'FLUXES    '.AND.ymode.EQ. 9).OR.
     .                 (tag(1:10).EQ.'ION       '.AND.ymode.EQ.10).OR.
     .                 (tag(1:10).EQ.'REC       '.AND.ymode.EQ.11).OR.
     .                 (tag(1:10).EQ.'CFP       '.AND.ymode.EQ.12)) THEN
                   mvals(in,ip,1) = rdum1
                   mvals(in,ip,2) = rdum2
                   mvals(in,ip,3) = rdum3
                   mvals(in,ip,4) = rdum4
                   mvals(in,ip,5) = rdum5
                 ELSEIF (tag(1:10).EQ.'MOM LOSS  '.AND.ymode.EQ.16) THEN
                   mvals(in,ip,1) = rdum1
                   mvals(in,ip,2) = rdum2
                 ENDIF
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF


            ELSEIF (tag(1:10).EQ.'RELAX     '.AND.
     .              (ymode.EQ.13)) THEN
c
c 13 - peimul
c
              IF     (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*) tag,step,iter,count,ir1,
     .                        rdum1,rdum2,rdum3,rdum4
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF

              IF (ir.EQ.ir1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  IF     (ymode.EQ.13) THEN
                    mvals(in,ip,1) = rdum3
                    mvals(in,ip,2) = rdum4
                  ENDIF
                ENDIF
                liter  = iter
                lstep  = step
                lcount = count
              ENDIF


            ELSEIF ((tag(1:10).EQ.'ACD 01    '.OR.
     .               tag(1:10).EQ.'ACD 03    ').AND.
     .              (ymode.EQ.14.OR.ymode.EQ.15.OR.ymode.EQ.18.OR.
     .               ymode.EQ.19)) THEN
c
c 14 - D and D2 temperature in additional cells
c 15 - D and D2 density
c 18 - additional cell pressure
c 19 - D2 temperature from ... and ...
c
              IF     (tag(11:14).EQ.'1.01') THEN
                READ(ldata,*) tag,step,iter,count,eircount,idum1,idum2,
     .                        idum3
              ELSEIF (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*) tag,step,iter,count,idum1,idum2,idum3
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF

              IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

              IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                in           = in + 1
                mouts(in,ip) = REAL(count)
              ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                in           = in + 1
                mouts(in,ip) = REAL(step)
              ELSEIF (xmode.EQ.3.AND.eircount.NE.leircount.AND.
     .                eircount.GE.xbnd1.AND.eircount.LE.xbnd2) THEN
                in           = in + 1
                mouts(in,ip) = REAL(eircount)
              ENDIF

              IF (in.GT.0) THEN
c...search for data
c                DO i1 = 1, idum3
c                  DO i2 = 1, idum2
                DO WHILE (.TRUE.) 
 
                  READ (ifp,'(A2048)',END=1000) ldata
		
c...              Loop exit condition - another tag is found:
                  IF (ldata(1:1).EQ.'''') THEN
                    BACKSPACE ifp
                    EXIT
                  ENDIF
                  READ (ldata,*) idum4,idum5,rdum1,rdum2,rdum3,rdum4,
     .                           rdum5
		
c                  READ (ifp,*) idum4,idum5,rdum1,rdum2,rdum3,rdum4,
c     .                         rdum5
		
                  IF     (ymode.EQ.14.AND.idum5.EQ.ir) THEN
c...                Temperature (...):
                    rdum1 = rdum1 + 1.0E-10
                    IF (tag(1:6).EQ.'ACD 01') THEN
                      IF (idum4.EQ.1) mvals(in,ip,1) = rdum2/rdum1*0.667
                      IF (idum4.EQ.3) mvals(in,ip,3) = rdum2/rdum1*0.667
                    ENDIF
                    IF (tag(1:6).EQ.'ACD 03') THEN
                      IF (idum4.EQ.1) mvals(in,ip,2) = rdum2/rdum1*0.667
                      IF (idum4.EQ.3) mvals(in,ip,4) = rdum2/rdum1*0.667
                    ENDIF
                  ELSEIF (ymode.EQ.15.AND.idum5.EQ.ir) THEN
c...                Density (convert from cm-3 to m-3):
                    IF (tag(1:6).EQ.'ACD 01') THEN
                      IF (idum4.EQ.1) mvals(in,ip,1) = rdum1 * 1.0E+06
                      IF (idum4.EQ.3) mvals(in,ip,3) = rdum1 * 1.0E+06
                    ENDIF
                    IF (tag(1:6).EQ.'ACD 03') THEN
                      IF (idum4.EQ.1) mvals(in,ip,2) = rdum1 * 1.0E+06
                      IF (idum4.EQ.3) mvals(in,ip,4) = rdum1 * 1.0E+06
                    ENDIF
                  ELSEIF (ymode.EQ.18.AND.idum5.EQ.ir) THEN
c                  ELSEIF (ymode.EQ.18.AND.i2.EQ.ir) THEN
c...                Pressure (convert eV cm-3 to mTorr):
                    fact = ECH * 1.0E+06 * 0.67 * 7.502
                    IF (tag(1:6).EQ.'ACD 01') THEN
                      IF (idum4.EQ.1) mvals(in,ip,1) = rdum2 * fact
                      IF (idum4.EQ.3) mvals(in,ip,3) = rdum2 * fact
c                      IF (i1.EQ.1) mvals(in,ip,1) = rdum2 * fact
c                      IF (i1.EQ.3) mvals(in,ip,3) = rdum2 * fact
                    ENDIF
                    IF (tag(1:6).EQ.'ACD 03') THEN
                      IF (idum4.EQ.1) mvals(in,ip,2) = rdum2 * fact
                      IF (idum4.EQ.3) mvals(in,ip,4) = rdum2 * fact
c                      IF (i1.EQ.1) mvals(in,ip,2) = rdum2 * fact
c                      IF (i1.EQ.3) mvals(in,ip,4) = rdum2 * fact
                    ENDIF
                  ELSEIF (ymode.EQ.19.AND.i2.EQ.ir) THEN
c...                D2 temperature from ... and ...:
                    rdum1 = rdum1 + 1.0E-10
                    IF (tag(1:6).EQ.'ACD 01') THEN
                      IF (i1.EQ.3) mvals(in,ip,1) = rdum2/rdum1*0.667
                      IF (i1.EQ.8) mvals(in,ip,2) = rdum2
                    ENDIF
		
                  ENDIF

                ENDDO
1000            CONTINUE

c                  ENDDO
c                ENDDO
              ENDIF
              liter     = iter
              lstep     = step
              lcount    = count
              leircount = eircount
c
c
c
            ELSEIF ((tag(1:10).EQ.'ASDATA    '.OR.
     .               tag(1:10).EQ.'PG-DATA   ').AND.ymode.EQ.17) THEN
c...          17 - pressure gauge data

              IF     (tag(11:14).EQ.'1.01') THEN
                READ(ldata,*) tag,step,iter,count,in1,
     .                        (radum1(i1),i1=1, 7)
                READ(ifp  ,*) cdum1,(radum1(i1),i1=8,14)
              ELSEIF (tag(11:14).EQ.'1.00') THEN
                READ(ldata,*) tag,step,iter,count,in1,
     .                        (radum1(i1),i1=1,14)
              ELSE
                CALL ER('970','Pressure gauge data version not '//
     .                        'supported',*99)
              ENDIF

              IF (ir.EQ.in1) THEN
                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)
c...            Assign independent variable (either iteration number 
c               or step number):
                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
c...              Multiply by 7.502 to convert from Pa to mTorr:
                  mvals(in,ip,1) = radum1(13) * 7.502
                  mvals(in,ip,2) = radum1(10) * 7.502
c                  mvals(in,ip,1) = radum1(14) * 7.502
c                  mvals(in,ip,2) = radum1(11) * 7.502
                  mvals(in,ip,3) = mvals(in,ip,1) + mvals(in,ip,2)
c....             Experimental data is expected to be in mTorr:
                  mvals(in,ip,4) = radum1( 7) 
                ENDIF

c...            Assign plot labels:
                IF (in.EQ.1) THEN
                  pltlabs(ip) = gaugename(INT(radum1(1)))
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF


            ELSEIF (tag(1:10).EQ.'PINFLUX   '.AND.
     .              (ymode.EQ.20)) THEN
c
c 20 - particle flux
c
c              IF     (tag(11:14).EQ.'1.00') THEN
c                READ(ldata,*) tag,step,iter,count,i1,i2,iliin,
c     .                        atm1,atm2,atm3,mol1,mol2,mol3
              IF (tag(11:14).EQ.'1.01') THEN
                READ(ldata,*) tag,step,iter,count,i1,i2,iliin,
     .                        atm(1),atm(2),atm(3),mol(1),mol(2),mol(3),
     .                        atm(4),atm(5),atm(6),mol(4),mol(5),mol(6)
              ELSE
                CALL ER('970','Rel. y-data version unsupported',*99)
              ENDIF


              IF (ir.EQ.i1) THEN

                IF (in.GE.MAXNKS-1) CALL ER('970','Overflow',*99)

                IF     (xmode.EQ.1.AND.count.NE.lcount.AND.
     .                  count.GE.xbnd1.AND.count.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(count)
                ELSEIF (xmode.EQ.2.AND.step.NE.lstep.AND.
     .                  step.GE.xbnd1.AND.step.LE.xbnd2) THEN
                  in           = in + 1
                  mouts(in,ip) = REAL(step)
                ENDIF

                IF (in.GT.0) THEN
                  IF     (ymode.EQ.20) THEN
                    mvals(in,ip,1) = atm(1)
                    mvals(in,ip,2) = atm(4)
                    mvals(in,ip,3) = mol(1)
                    mvals(in,ip,4) = mol(4)
                    mvals(in,ip,5) = atm(1)+atm(4)+mol(1)+mol(4)
c                    mvals(in,ip,1) = atm1 / 5.0E+21
c                    mvals(in,ip,2) = mol1 / 5.0E+21
                  ENDIF
                ENDIF

                liter  = iter
                lstep  = step
                lcount = count
              ENDIF



            ENDIF

          ENDDO
2500      CONTINUE
c
c
c
c
          pnks(ip) = in

          DO ik = 1, pnks(ip)
            WRITE(6,'(2I4,2X,2I4,2X,2F8.4,2X,1P,10(E11.3:))')
     .        ik,ir,ip,in,
     .        mouts(ik,ip),mwids(ik,ip),(mvals(ik,ip,i1),i1=1,ngs)
          ENDDO


        ENDDO

        CLOSE(ifp)

        DO ip = 1, nplts
          i2 = 0
          DO i1 = 1, ngs
            status = .FALSE.
            DO ik = 1, pnks(ip)
              IF (mvals(ik,ip,i1).NE.0.0) status = .TRUE.
            ENDDO
            IF (.NOT.status) THEN
c              i2 = i2 + 1
              WRITE(6,*) '970: NULL DATA ARRAY DETECTED: IP,I1 = ',ip,i1
              CALL RSet(mvals(1,ip,i1),MAXNKS,LO)
            ENDIF
          ENDDO
        ENDDO


c...    Change plot labels:    
        CALL CustomisePlot(title,xlab,ylab,elabs)

      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Xran'.OR.dummy(8:11).EQ.'xran'.OR.
     .      dummy(8:11).EQ.'XRAN') THEN
        READ(dummy,*) cdum1,xrange1,xrange2
        DO ip = 1, nplts
c...      Shift data points to make room for x-axis range points:
          DO ik = pnks(ip)+1, 2
            mouts(ik,ip) = mouts(ik-1,ip)
            DO i1 = 1, ngs
              mvals(ik,ip,i1) = mvals(ik-1,ip,i1)
            ENDDO
          ENDDO
          pnks(ip) = pnks(ip) + 2
c...      Add new data points that mark x-axis range:
          mouts(1       ,ip) = xrange1
          mouts(pnks(ip),ip) = xrange2
          DO i1 = 1, ngs
            mvals(1       ,ip,i1) = LO
            mvals(pnks(ip),ip,i1) = LO
          ENDDO          
        ENDDO
      ELSE
        BACKSPACE 5
      ENDIF


c...  Add a caption to the plot - THIS READ SHOULD BE INSIDE ADDCAPTION:
      READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Note'.OR.dummy(8:11).EQ.'note'.OR.
     .      dummy(8:11).EQ.'NOTE') THEN
        READ(dummy,*) cdum1,xpos,ypos,size,caption
        CALL AddCaption(caption,xpos,ypos,size)
      ELSE
        BACKSPACE 5
      ENDIF


        IF (i2.LT.ngs) THEN

          CALL SetPlotComments(970,job,extra_comments,0,tarshift(IKHI))

          CALL SLDRAWM (mouts,mwids,mvals,MAXNKS,pnks,
     .                nplts,ngs,pltlabs,elabs,xlab,ylab,ref,title,
     >                sctype,ngrm,pltmins,pltmaxs,pltfact)
c     .                grm_opt,ngs2,grm_shade,elabs2)
        ELSE
          WRITE(0,*) '970: No data, skipping plot'
        ENDIF
c        WRITE(0,*) '970: DONE'


      RETURN
99    STOP
      END
