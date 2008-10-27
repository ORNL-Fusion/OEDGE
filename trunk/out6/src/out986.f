c     -*-Fortran-*-
c
c
c
      SUBROUTINE Plot986(graph,nplots2,ref,title,iopt,pltmins,pltmaxs,
     .                   iplot,job)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'
      INCLUDE 'comgra'
c      INCLUDE 'outcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER    MAXEXPNDATA
      PARAMETER (MAXEXPNDATA=5000)

      REAL CalcPressure,GetJsat,GetCs,GetIonSrc

      INTEGER ik,ir

      integer iplot,NPLOTS2,iopt,ndata(MAXNGS),i1,i2,i3,ngs,loc1,loc2,
     .        type,plotformat,plotmode,ndata1,i1start,plotframe,
     .        idum1,idum2,idum3,id,ik1,size,anotemode,fudge,step,
     .        index,colselect,ecount,navg,stratum,region,expndata
      LOGICAL oldraw,status,cumulative
      real pltmins(maxplts),pltmaxs(maxplts),ydata(MAXGXS,MAXNGS),
     .     xdata(MAXGXS),TWIDS(MAXGXS),ydata1(MAXGXS,MAXNGS),rdum1,
     .     expxdata(MAXEXPNDATA),expydata(MAXEXPNDATA),
     .     xdata1(MAXGXS,MAXNGS),probfact(MAXNGS),
     .     probescaledat(MAXGXS,MAXNGS),probescaledat2(MAXGXS,MAXNGS),
     .     ydata2(MAXGXS,MAXNGS,10),flux,area,vavg,fact,mass,
     .     xavg(MAXGXS),yavg(MAXGXS),psin,deltapsin
      REAL AVS(0:100),maxy,scale,xrange1,xrange2,xpos,ypos,sum1,
     .     yrange1,yrange2
      INTEGER IGNORS(MAXNGS),ITEC,NAVS,nydata2,ndataydata2(MAXNGS,10),
     .        ngsydata2(10)
      INTEGER probescale
      CHARACTER TITLE*80,JOB*72,GRAPH*80,GRAPH1*80,pattern*128,file*128,
     .          dummy*1024,cdum1*32,mark*256,caption*5000
      CHARACTER*36 XLAB
      CHARACTER*128 cnames(MAXGXS)
      CHARACTER*32 tag(100)
      CHARACTER*64 ylab
      CHARACTER*36 REF,VIEW,ANLY,PLANE,table
      CHARACTER*36 NAME,pnames1(MAXGXS),pnames2(MAXGXS)
      CHARACTER*26 alphabet
      CHARACTER*128 elabs(MAXNGS)
      CHARACTER cname*256,cmnd*256

      INTEGER i,optflow
     
      REAL, ALLOCATABLE :: tmppinion(:,:)
      REAL, ALLOCATABLE :: tmposmcfp(:,:)
      REAL, ALLOCATABLE :: tmpkvhs(:,:)

      DATA alphabet /'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      DATA probescale,nydata2 /0,0/


      REAL totalleakage1,totalleakage2,totalleakage3,totalleakage4


c FSP
      integer osmvals,osmplots
      integer int_type,icnt
      real    r1p,z1p,r2p,z2p,rsect,zsect
      REAL LOUTS(MAXSEG),LVALS(MAXSEG,MAXNGS)
      SAVE

c EXPT
      INTEGER ntmp
      REAL    xtmp(5000),ytmp(5000)

      CALL THICK2(4)

      DO i1 = 1, MAXGXS
        WRITE(cnames(i1),'(128X)') 
        WRITE(pnames1(i1),'(36X)') 
        WRITE(pnames2(i1),'(36X)') 
      ENDDO

      iplots = 0
      ngs = 0
      fudge = 0
      xrange1 = HI
      xrange2 = HI
      yrange1 = -HI
      yrange2 =  HI
      oldraw = .FALSE.
      cumulative = .FALSE.

      CALL RZero(twids,MAXGXS)
      CALL ISet(ignors,MAXNGS,1)

      CALL RSet(ydata ,MAXGXS*MAXNGS,LO)
      CALL RSet(ydata1,MAXGXS*MAXNGS,LO)

      plotframe = 1


c...  Load plot parameters:
c       PLOTFORMAT	1 - Each data line is a separate group of plot data
c			2 - Each data line is a data point on the same line

c       PLOTMODE	1 - line chat
c   			2 - bar chart


      CALL CustomisePlot(title,xlab,ylab,elabs)

10    READ(5,'(A256)',END=15) dummy
      IF   (dummy(8:11).EQ.'Plot'.OR.dummy(8:11).EQ.'plot'.OR.
     .      dummy(8:11).EQ.'PLOT') THEN
        READ(dummy,*) cdum1,plotformat,plotmode,plotframe
        GOTO 10
      ELSE
        BACKSPACE 5
      ENDIF
15    CONTINUE

      READ(5,'(A256)',END=25) dummy
      IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .      dummy(8:11).EQ.'SIZE') THEN
        READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
      ELSE
        map1x = -1.0
        BACKSPACE 5
      ENDIF
25    CONTINUE

c...  X-axis labels:
      IF (plotframe.EQ.1.OR.plotframe.EQ.2) THEN
        DO i1 = 1, 26
          pnames1(i1) = alphabet(i1:i1)
        ENDDO
      ENDIF

30    READ(5,'(A256)',END=35) dummy
      IF   (dummy(8:11).EQ.'Labs'.OR.dummy(8:11).EQ.'labs'.OR.
     .      dummy(8:11).EQ.'LABS') THEN
        READ(dummy,*) cdum1,idum1,pnames1(idum1),pnames2(idum1)
        GOTO 30
      ELSE
        BACKSPACE 5
      ENDIF
35    CONTINUE

40    READ(5,'(A256)',END=35) dummy
      IF   (dummy(8:11).EQ.'Clab'.OR.dummy(8:11).EQ.'clab'.OR.
     .      dummy(8:11).EQ.'CLAB') THEN
        READ(dummy,*) cdum1,idum1,cnames(idum1)
        GOTO 30
      ELSE
        BACKSPACE 5
      ENDIF
45    CONTINUE

      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Xran'.OR.dummy(8:11).EQ.'xran'.OR.
     .      dummy(8:11).EQ.'XRAN') THEN
        READ(dummy,*) cdum1,xrange1,xrange2
      ELSE
        BACKSPACE 5
      ENDIF

      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Yran'.OR.dummy(8:11).EQ.'yran'.OR.
     .      dummy(8:11).EQ.'YRAN') THEN
        READ(dummy,*) cdum1,yrange1,yrange2
        IF (yrange1.EQ.99.0) yrange1 = -HI
        IF (yrange2.EQ.99.0) yrange2 =  HI
      ELSE
        BACKSPACE 5
      ENDIF

c...  Load plot data from data files in the results directory:
c
      i1start = 1
      scale = 1.0
      WRITE(mark,'(256X)') 
      CALL IZero(ndata,MAXNGS)
50    READ(5,'(A256)',END=55) dummy

      IF   (dummy(8:11).EQ.'Scal'.OR.dummy(8:11).EQ.'scal'.OR.
     .      dummy(8:11).EQ.'SCAL') THEN
        READ(dummy,*) cdum1,idum1
        IF     (idum1.EQ.-1) THEN
c..HACK-ATTACK!
            probescale = 1
        ELSEIF (idum1.EQ.0) THEN
c...      Reset scale:
          READ(dummy,*) cdum1,idum1,rdum1
          scale = rdum1
          fudge = 0
        ELSEIF (idum1.EQ.1) THEN
c...      Multiplication:
          READ(dummy,*) cdum1,idum1,rdum1
          scale = scale * rdum1
        ELSEIF (idum1.EQ.2) THEN
c...      Multiplication, but load scaling factor from data file:
          WRITE(mark,'(256X)') 
          READ(5,'(A256)') dummy
          IF   (dummy(8:11).EQ.'Mark'.OR.dummy(8:11).EQ.'mark'.OR.
     .          dummy(8:11).EQ.'MARK') THEN
            READ(dummy,*) cdum1,mark
          ELSE
            BACKSPACE 5
          ENDIF
          READ(5,'(A256)') dummy
          READ(dummy,*) cdum1,file,pattern,loc1,loc2,type
          CALL SnatchData(file,mark,pattern,loc1,loc2,type,
     .                    ndata1,ydata1(1,1))
          scale = scale * ydata1(1,1)
c...      Terrible:
          IF (ylab(1:3).EQ.'p (') THEN
            probescale = 1
          ENDIF
        ELSEIF (idum1.EQ.3) THEN
c...      Oh golly:
        ELSEIF (idum1.EQ.4) THEN
c...      Multiplication, but load scaling factor from data file:
          WRITE(mark,'(256X)') 
          READ(5,'(A256)') dummy
          IF   (dummy(8:11).EQ.'Mark'.OR.dummy(8:11).EQ.'mark'.OR.
     .          dummy(8:11).EQ.'MARK') THEN
            READ(dummy,*) cdum1,mark
          ELSE
            BACKSPACE 5
          ENDIF
          READ(5,'(A256)') dummy
          READ(dummy,*) cdum1,file,pattern,loc1,loc2,type
          CALL SnatchData(file,mark,pattern,loc1,loc2,type,
     .                    ndata1,ydata1(1,1))
          scale = scale * 1.0 / (1.0 - ydata1(1,1))
        ELSEIF (idum1.EQ.5) THEN
c...      Oh golly, fudge:
          fudge = 1
        ELSEIF (idum1.EQ.6) THEN
c...      Oh golly, fudge:
          fudge = 2
        ELSEIF (idum1.EQ.7) THEN
c...      Oh golly, fudge:
          fudge = 3
        ELSEIF (idum1.EQ.8) THEN
c...      Calculating momentum flux into plenum:
          WRITE(0,*) 'SUPERSIZE!'
          fudge = 4
        ELSEIF (idum1.EQ.9) THEN
c...      Calculating momentum flux into plenum:
          WRITE(0,*) 'SUPERSIZE! TO GO!'
          fudge = 5
        ELSEIF (idum1.EQ.10) THEN
c...      Calculating leakage per section:
          WRITE(0,*) 'NET LEAKAGE!'
          fudge = 6
        ELSEIF (idum1.EQ.11) THEN
c...      Oh golly, fudge:
          fudge = 7
        ELSE
          CALL ER('986','Invalid scaling option',*99)
        ENDIF
        GOTO 50

      ELSEIF   (dummy(8:11).EQ.'Mark'.OR.dummy(8:11).EQ.'mark'.OR.
     .          dummy(8:11).EQ.'MARK') THEN
        WRITE(mark,'(256X)') 
        READ(dummy,*) cdum1,mark
        GOTO 50

      ELSEIF   (dummy(8:11).EQ.'Data'.OR.dummy(8:11).EQ.'data'.OR.
     .          dummy(8:11).EQ.'DATA') THEN
        READ(dummy,*) cdum1,file,pattern,loc1,loc2,type

        CALL SnatchData(file,mark,pattern,loc1,loc2,type,
     .                  ndata1,ydata1(1,1))

        IF     (plotformat.EQ.1) THEN
          IF (ngs.EQ.0) THEN
            DO i1 = 1, ndata1
              xdata(i1) = REAL(i1)
              pnames1(i1) = alphabet(i1:i1)
            ENDDO
          ENDIF
          ngs = ngs + 1
          ndata(ngs) = ndata1
          DO i1 = 1, ndata1
            ydata(i1,ngs) = ydata(i1,ngs) + ydata1(i1,1) * scale
          ENDDO
        ELSEIF (plotformat.EQ.2) THEN
          IF (ngs.EQ.0) ngs = ndata1
          DO i1 = i1start, ngs
            IF (i1start.EQ.1) THEN
              ndata(i1) = ndata(i1) + 1
              ydata(ndata(i1),i1) = ydata(ndata(i1),i1) + 
     .                              ydata1(i1,1) * scale
            ELSE
              DO i2 = 1, ndata1
                ndata(i1) = ndata(i1) + 1
                ydata(ndata(i1),i1) = ydata(ndata(i1),i1) + 
     .                                ydata1(i2,1) * scale
              ENDDO
            ENDIF
            xdata(ndata(i1)) = REAL(ndata(i1))
          ENDDO

          READ(5,'(A256)',END=55) dummy
          IF   (dummy(8:11).EQ.'Next'.OR.dummy(8:11).EQ.'next'.OR.
     .          dummy(8:11).EQ.'NEXT') THEN
            ngs = ngs + 1
            i1start = ngs
            scale = 1.0
            fudge = 0
c            WRITE(0,*) 'INCREASING NGS:',ngs
          ELSE
            BACKSPACE 5
          ENDIF

          READ(5,'(A256)',END=55) dummy
          IF   (dummy(8:10).EQ.'Sum'.OR.dummy(8:10).EQ.'sum'.OR.
     .          dummy(8:10).EQ.'SUM') THEN
            ndata(ngs) = 0
            i1start = ngs

c            WRITE(0,*) 'SUM ACTION'
          ELSE
            BACKSPACE 5
          ENDIF

        ELSE
          CALL ER('986','Invalid plot format',*99)
        ENDIF       

        GOTO 50
      ELSE
        BACKSPACE 5
      ENDIF
55    CONTINUE




    
      IF     (iopt.EQ.1) THEN
c...    Plot the fraction of each stratum ionised in 
c       each plasma regon:

c        nregions = 3
c        nstrata = 3




      ELSEIF (iopt.EQ.2) THEN
c...

        IF (map1x.EQ.-1.0) THEN
          map1x=0.1
          map2x=0.9
          map1y=0.1
          map2y=0.4
        ENDIF

        iopt_ghost = 2

        CALL CTRMAG(12)

      ELSEIF (iopt.EQ.3) THEN
c...

        map1x=0.1
        map2x=0.9
        map1y=0.1
        map2y=0.4

        iopt_ghost = 2

        CALL CTRMAG(12)

        ngs = 4
                
        CALL SnatchData('tor0101?.dat',mark,'2443   0.604 -0.661',81,87,
     .                  1,ndata1,ydata1(1,1))				
									
        CALL SnatchData('tor0102?.dat',mark,'2443   0.604 -0.661',81,87,
     .                  1,ndata1,ydata1(1,1))				
									
        CALL SnatchData('tor0103?.dat',mark,'2443   0.604 -0.661',81,87,
     .                  1,ndata1,ydata1(1,1))				
									
        CALL SnatchData('tor0104?.dat',mark,'2443   0.604 -0.661',81,87,
     .                  1,ndata1,ydata1(1,1))


      ELSEIF (iopt.EQ.4) THEN
c...    TEMP: Assigning rho:

        IF (map1x.EQ.-1.0) THEN
          map1x=0.1
          map2x=0.9
          map1y=0.1
          map2y=0.4
        ENDIF

        iopt_ghost = 1

        CALL CTRMAG(12)

c...    Assign RHO to y coordinate:
        CALL CalculateRho  ! remove eventually, data now passed from DIVIMP...
        ir = irtrap+1
        i1 = 1
        DO WHILE(ir.NE.irwall)        
          xdata(i1) = rho(ir,CELL1) * 1000.0
          ir = irouts(1,ir)
          i1 = i1 + 1
        ENDDO


      ELSEIF (iopt.EQ.10) THEN
c...    Inner target Te:
        CALL CalculateRho
        ngs = ngs + 1
        ir = irtrap+1
        ndata(ngs) = 0
        DO WHILE(ir.NE.irwall)        
          ndata(ngs) = ndata(ngs) + 1
          xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
          ydata1(ndata(ngs),ngs) = kteds(idds(ir,2)) * scale
          ir = irouts(1,ir)
        ENDDO

      ELSEIF (iopt.EQ.11) THEN
c...    Inner target Isat:
        CALL CalculateRho
        ngs = ngs + 1
        ir = irtrap+1
        ndata(ngs) = 0
        DO WHILE(ir.NE.irwall)        
          ndata(ngs) = ndata(ngs) + 1
          xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
          id = idds(ir,2)
          ydata1(ndata(ngs),ngs) = scale *
     .      ABS(GetJsat(kteds(id),ktids(id),knds(id),kvds(id))) 
          ir = irouts(1,ir)
        ENDDO

      ELSEIF (iopt.EQ.12) THEN
c...    Outer target Te:
        status = .TRUE.
        DO WHILE(status)
          status = .FALSE.
          CALL CalculateRho
          ngs = ngs + 1
          ir = irtrap+1
          ndata(ngs) = 0
          DO WHILE(ir.NE.irwall)        
            ndata(ngs) = ndata(ngs) + 1
            xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
            ydata1(ndata(ngs),ngs) = kteds(idds(ir,1)) * scale
            ir = irouts(1,ir)
          ENDDO
c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          READ(5,'(A128)') dummy
          IF (dummy(8:11).EQ.'Case'.OR.dummy(8:11).EQ.'CASE'.OR.
     .        dummy(8:11).EQ.'case') THEN
            READ(dummy,*) cdum1,cdum1,cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
c...        Restore main .raw file if necessary:
            IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

            BACKSPACE 5
          ENDIF
        ENDDO

      ELSEIF (iopt.EQ.13) THEN
c...    Outer target Isat:
        CALL CalculateRho
        ngs = ngs + 1
        ir = irtrap+1
        ndata(ngs) = 0
        DO WHILE(ir.NE.irwall)        
          ndata(ngs) = ndata(ngs) + 1
          xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
          id = idds(ir,1)
          ydata1(ndata(ngs),ngs) = scale *
     .      ABS(GetJsat(kteds(id),ktids(id),knds(id),kvds(id)))
          ir = irouts(1,ir)
        ENDDO

      ELSEIF (iopt.EQ.14) THEN
c...    Inner target nmax:

c... REALLY NEED OSM_MODEL FOR THIS PLOT, BUT IT IS NOT PASSED IN .RAW

        CALL CalculateRho
        ngs = ngs + 1
        ir = irtrap+1
        ndata(ngs) = 0
        DO WHILE(ir.NE.irwall)        
          ndata(ngs) = ndata(ngs) + 1
          xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
          IF (ir.LT.irwall) THEN
c...        SOL:
            ydata1(ndata(ngs),ngs) = knbs(ikbound(ir,IKLO),ir) * scale
          ELSE
c...        PFZ:
            DO ik = 2, nks(ir)
              IF (ktebs(ik,ir).GT.ktebs(ik-1,ir)) EXIT
              ydata1(ndata(ngs),ngs) = knbs(ik,ir) * scale
            ENDDO
          ENDIF
          ir = irouts(1,ir)
        ENDDO

      ELSEIF (iopt.EQ.15) THEN
c...    Tmax:

c... REALLY NEED OSM_MODEL FOR THIS PLOT, BUT IT IS NOT PASSED IN .RAW
        status = .TRUE.
        DO WHILE(status)
          status = .FALSE.
          CALL CalculateRho
          ngs = ngs + 1
          ir = irtrap+1
          ndata(ngs) = 0
          DO WHILE(ir.NE.irwall)        
            ndata(ngs) = ndata(ngs) + 1
            xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0


c            ydata1(ndata(ngs),ngs) = kteds(idds(ir,1)) * scale

            IF (ir.LT.irwall) THEN
c...          SOL:
              IF (ir.EQ.irsep) THEN
                ydata1(ndata(ngs),ngs) = ktebs(ik1,ir) * scale
               ELSE
                ydata1(ndata(ngs),ngs) = LO
              ENDIF
            ELSE
c...          PFZ:
              ydata1(ndata(ngs),ngs) = -HI
              DO ik = 1, nks(ir)
                IF (ktebs(ik,ir)*scale.GT.ydata1(ndata(ngs),ngs)) THEN
                  ydata1(ndata(ngs),ngs) = ktebs(ik,ir) * scale
                  ik1 = ikouts(ik,ir)
                ENDIF
              ENDDO
            ENDIF


            ir = irouts(1,ir)
          ENDDO
c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          READ(5,'(A128)') dummy
          IF (dummy(8:11).EQ.'Case'.OR.dummy(8:11).EQ.'CASE'.OR.
     .        dummy(8:11).EQ.'case') THEN
            READ(dummy,*) cdum1,cdum1,cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
c...        Restore main .raw file if necessary:
            IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

            BACKSPACE 5
          ENDIF
        ENDDO


c        CALL CalculateRho
c        ngs = ngs + 1
c        ir = irtrap+1
c        ndata(ngs) = 0
c        DO WHILE(ir.NE.irwall)        
c          ndata(ngs) = ndata(ngs) + 1
c          xdata(ndata(ngs)    ) = rho(ir,CELL1) * 1000.0
c          IF (ir.LT.irwall) THEN
cc...        SOL:
c            IF (ir.EQ.irsep) THEN
c              ydata1(ndata(ngs),ngs) = ktebs(ik1,ir) * scale
c            ELSE
c              ydata1(ndata(ngs),ngs) = LO
c            ENDIF
c          ELSE
cc...        PFZ:
c            ydata1(ndata(ngs),ngs) = -HI
c            DO ik = 1, nks(ir)
c              IF (ktebs(ik,ir)*scale.GT.ydata1(ndata(ngs),ngs)) THEN
c                ydata1(ndata(ngs),ngs) = ktebs(ik,ir) * scale
c                ik1 = ikouts(ik,ir)
c              ENDIF
c            ENDDO
c          ENDIF
c          ir = irouts(1,ir)
c        ENDDO

      ELSEIF (iopt.EQ.16) THEN
c...    PFZ and outer SOL ionisation:
        cumulative = .TRUE.
        CALL CalculateRho
        DO i1 = 1, ngs
          ndata(i1) = 0
          ir = irtrap+1
          DO WHILE(ir.NE.irwall)        
            ndata(i1) = ndata(i1) + 1
            xdata1(ndata(i1),i1) = rho(ir,CELL1) * 1000.0
            ydata1(ndata(i1),i1) = MAX(0.0,ydata(ndata(i1),i1))
            ir = irouts(1,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.17) THEN
c...    Outer target pressure:
        status = .TRUE.
        DO WHILE(status)
          status = .FALSE.
          CALL CalculateRho
          ngs = ngs + 1
          ir = irtrap+1
          ndata(ngs) = 0
          DO WHILE(ir.NE.irwall)        
            ndata(ngs) = ndata(ngs) + 1
            xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
            ydata1(ndata(ngs),ngs) = scale * ECH *
     .        CalcPressure(knds (idds(ir,1)),kteds(idds(ir,1)),
     .                     ktids(idds(ir,1)),kvds (idds(ir,1)))
            ir = irouts(1,ir)
          ENDDO
c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          READ(5,'(A128)') dummy
          IF (dummy(8:11).EQ.'Case'.OR.dummy(8:11).EQ.'CASE'.OR.
     .        dummy(8:11).EQ.'case') THEN
            READ(dummy,*) cdum1,cdum1,cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
c...        Restore main .raw file if necessary:
            IF (oldraw) CALL LoadData(' ',' ',NULL,-2)

            BACKSPACE 5
          ENDIF
        ENDDO

      ELSEIF (iopt.GE.18.AND.iopt.LE.22.OR.
     .        iopt.GE.40.AND.iopt.LE.45) THEN
c...    Reciprocating probe data:
        status = .TRUE.
        DO WHILE(status)
          status = .FALSE.
          CALL CalculateRho

c...      Check if parallel plasma flow should be over-written:
          optflow = 0
          READ(5,'(A80)') graph1
          IF (graph1(8:11).EQ.'Flow'.OR.graph1(8:11).EQ.'FLOW'.OR.
     .        graph1(8:11).EQ.'flow') THEN
            READ(graph1,*) cdum1,optflow

c...        Store ionisation source, PINION:
            ALLOCATE(tmppinion(MAXNKS,MAXNRS))
            ALLOCATE(tmposmcfp(MAXNKS,MAXNRS))
            ALLOCATE(tmpkvhs(MAXNKS,MAXNRS))
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                tmppinion(ik,ir) = pinion(ik,ir)
                tmposmcfp(ik,ir) = osmcfp(ik,ir)
                tmpkvhs(ik,ir) = kvhs(ik,ir)

                pinion(ik,ir) = osmion(ik,ir)
              ENDDO
            ENDDO
            CALL CalcFlow(optflow)

c            DO ir = 1, nrs
c              DO ik = 1, nks(ir)
c                IF (ir.GE.irwall) kvhs(ik,ir) = 2.0
c                WRITE(0,*) 'KVHS:',ik,ir,kvhs(ik,ir)
c              ENDDO
c            ENDDO


          ELSE
            BACKSPACE 5
          ENDIF

          IF     (machine.EQ.CMOD) THEN
c...        C-Mod vertical scanning probe:
            r1p =  0.732
            z1p = -0.325
            r2p =  0.732
            z2p = -0.280 
c...        C-Mod horizontal scanning probe:
c            r1p =  0.8500
c            z1p =  0.1088 
c            r2p =  0.8900
c            z2p =  0.1088

          ELSEIF (machine.EQ.DIIID) THEN
c...        DIII-D horizontal scanning probe:
            r1p =  1.500
            z1p = -0.188
            r2p =  3.000
            z2p = -0.188
          ELSE
            CALL ER('986:18-22','Invalid MACHINE specified',*99)
          ENDIF

          int_type = 3
          call osmprobe(lvals,louts,osmvals,osmplots,
     >                  r1p,z1p,r2p,z2p,int_type,crmb,qtim)

c          WRITE(0,*) '***',kvhs(63,22)/
c     .       GETCS(ktebs(63,22),ktebs(63,22))/qt

          ngs = ngs + 1
          IF (iopt.GE.40.AND.iopt.LE.45) THEN
            ir = irtrap + 1
          ELSE
            ir = irsep
          ENDIF
          ndata(ngs) = 0
          elabs(ngs) = '    OSM'
          WRITE(0,*) '*** 982: IGNORING IRWALL-1'
          DO WHILE(ir.NE.irwall-1)        
c          DO WHILE(ir.NE.irwall)        
            ndata(ngs) = ndata(ngs) + 1
            IF (machine.EQ.CMOD) THEN
c...          C-Mod vertical scanning probe, versus RHO:
              xdata1(ndata(ngs),ngs) = rho(ir,CELL1)
c              xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
            ELSE
c...          DIII-D horizontal scanning probe, versus PSIn:
              xdata1(ndata(ngs),ngs) = psitarg(ir,1)
            ENDIF

            ik = ndata(ngs)

            WRITE(0,*) '---',lvals(ik,1),lvals(ik,2),lvals(ik,3),
     .                       lvals(ik,5)

            IF    (iopt.EQ.18) THEN
c...          Pressure
              IF (lvals(ik,5).LT.0.0)
     .           lvals(ik,5) = MAX(lvals(ik,5),
     .                    -0.9*GetCs(lvals(ik,2),lvals(ik,3)))

              ydata1(ndata(ngs),ngs) = 
     .          CalcPressure(lvals(ik,1),lvals(ik,2),
     .                       lvals(ik,3),lvals(ik,5)) * ECH * scale
            ELSEIF (iopt.EQ.19) THEN
c...          n:
              ydata1(ndata(ngs),ngs) = lvals(ik,1) * scale
            ELSEIF (iopt.EQ.20) THEN
c...          Te:
              ydata1(ndata(ngs),ngs) = lvals(ik,2) * scale
            ELSEIF (iopt.EQ.21) THEN
c...          Ti:
              ydata1(ndata(ngs),ngs) = lvals(ik,3) * scale
            ELSEIF (iopt.EQ.22) THEN
c...          vb:
              ydata1(ndata(ngs),ngs) = lvals(ik,5) * scale /
     .                          GetCs(lvals(ik,2),lvals(ik,3))

              ydata1(ndata(ngs),ngs)=MAX(-0.9,ydata1(ndata(ngs),ngs))

c              ydata1(ndata(ngs),ngs)=1.0

              IF (optflow.EQ.5.AND.ir.EQ.irsep.OR.
     .                             ir.EQ.irsep+1) 
     .          ydata1(ndata(ngs),ngs)=MAX(-0.1,ydata1(ndata(ngs),ngs))


c              WRITE(0,*) '___:',ydata1(ndata(ngs),ngs),lvals(ik,5),
c     .                        GetCs(lvals(ik,2),lvals(ik,3))
            ELSEIF (iopt.EQ.40.OR.iopt.EQ.41) THEN
c...          Target Isat:
              IF (iopt.EQ.40) id = idds(ir,2)
              IF (iopt.EQ.41) id = idds(ir,1)
              ydata1(ndata(ngs),ngs) = scale *
     .          ABS(GetJsat(kteds(id),ktids(id),knds(id),kvds(id))) 

            ELSEIF (iopt.EQ.42.OR.iopt.EQ.43) THEN
c...          Target Isat:
              IF (iopt.EQ.42) id = idds(ir,2)
              IF (iopt.EQ.43) id = idds(ir,1)
              ydata1(ndata(ngs),ngs) = scale * kteds(id)

            ELSEIF (iopt.EQ.44.OR.iopt.EQ.45) THEN
c...          Target Isat:
              IF (iopt.EQ.44) id = idds(ir,2)
              IF (iopt.EQ.45) id = idds(ir,1)
              ydata1(ndata(ngs),ngs) = scale * ECH * 
     .          ABS(CalcPressure(knds(id),kteds(id),ktids(id),kvds(id))) 

            ENDIF
            ir = irouts(1,ir)
          ENDDO

          IF (optflow.NE.0) THEN
c...        Restore PINION, OSMCFP:
            DO ir = 1, nrs
              DO ik = 1, nks(ir)
                pinion(ik,ir) = tmppinion(ik,ir)
                osmcfp(ik,ir) = tmposmcfp(ik,ir)
                kvhs(ik,ir) = tmpkvhs(ik,ir)
              ENDDO
            ENDDO
            DEALLOCATE(tmppinion,STAT=i)
            DEALLOCATE(tmposmcfp,STAT=i)
          ENDIF

c...      Check if additional data is to be plotted from another DIVIMP
c         case:
          READ(5,'(A128)') dummy
          IF (dummy(8:11).EQ.'Case'.OR.dummy(8:11).EQ.'CASE'.OR.
     .        dummy(8:11).EQ.'case') THEN
            READ(dummy,*) cdum1,cdum1,cname,step,cmnd
            CALL LoadData(cname,cmnd,step,-1)
            status = .TRUE.
            oldraw = .TRUE.
          ELSE
c...        Restore main .raw file if necessary:
            IF (oldraw) THEN
              CALL LoadData(' ',' ',NULL,-2)
              WRITE(0,*) 'RESTORING BASE SOLUTION'
            ENDIF
            BACKSPACE 5
          ENDIF
        ENDDO

      ELSEIF (iopt.EQ.23.OR.iopt.EQ.24) THEN
c...    Integrated ionisation source as a fraction of the total ion sink
c       for the half-ring:
        cumulative = .TRUE.
        IF (iopt.EQ.23) region = IKLO
        IF (iopt.EQ.24) region = IKHI
        stratum = 6
        ngs = stratum
        CALL CalculateRho

        sum1 = 0.0

        DO i1 = 1, stratum
          ndata(i1) = 0
          ir = irtrap+1
          DO WHILE(ir.NE.irwall)        
            ndata(i1) = ndata(i1) + 1
c            xdata1(ndata(i1),i1) = REAL(ndata(i1))
            IF (region.EQ.IKLO) xdata1(ndata(i1),i1) = psitarg(ir,2)
            IF (region.EQ.IKHI) xdata1(ndata(i1),i1) = psitarg(ir,1)
c            xdata1(ndata(i1),i1) = rho(ir,CELL1) * 1000.0
            ydata1(ndata(i1),i1) = GetIonSrc(region,ir,i1,1)

c            IF (stratum.EQ.6) sum1 = sum1 + GetIonSrc(region,ir,0,0)

            ir = irouts(1,ir)
          ENDDO
        ENDDO
c        WRITE(0,*) 'SUM1=',sum1


      ELSEIF (iopt.EQ.25) THEN
c...    Flux of atoms to wall surfaces:

        CALL CalculateRho

        IF     (machine.EQ.CMOD) THEN
        ELSEIF (machine.EQ.DIIID) THEN
        ENDIF

        ngs = ngs + 1
        ir = irtrap+1
        ndata(ngs) = nvesm+nvesp
        DO i1 = 1, nvesm+nvesp
          xdata1(i1,ngs) = REAL(i1)
          ydata1(i1,ngs) = flxhw5(i1)
        ENDDO

c        ngs = ngs + 1
c        ir = irtrap+1
c        ndata(ngs) = 0
c        DO WHILE(ir.NE.irwall)        
c          ndata(ngs) = ndata(ngs) + 1
c          IF (machine.EQ.CMOD) THEN
c            xdata1(ndata(ngs),ngs) = rho(ir,CELL1) * 1000.0
c          ELSE
c            xdata1(ndata(ngs),ngs) = psitarg(ir,1)
c          ENDIF
c          ydata1(ndata(ngs),ngs) = ...
c        ENDDO

      ELSE
        CALL ER('986','Invalid plot option',*99)
      ENDIF



      IF (iopt.GE.10) THEN
c...    Process x-y scatter plots:

        slopt2 = 1
        DO i1 = 1, 9
          plottype(i1) = i1 + 1
        ENDDO
        plottype(1) = 55

        IF (plotframe.EQ.1.OR.plotframe.EQ.2) THEN
          iopt_ghost = 1
        ELSE
          iopt_ghost = 2
        ENDIF
        CALL CTRMAG(11)

        ndata1 = 0
        DO i1 = 1, ngs
          CALL LoadArray(xdata,ndata1,xdata1(1,i1),1,ndata(i1))
        ENDDO


        ecount = 0
75      READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'Expt'.OR.dummy(8:11).EQ.'expt'.OR.
     .        dummy(8:11).EQ.'EXPT') THEN

          ecount = ecount + 1
          ngs = ngs + 1

          dummy(1020:1022) = '999'
          READ(dummy,*) cdum1,index,colselect,plottype(ngs)

          IF (plottype(ngs).EQ.999) THEN
            IF (ecount.EQ.1) THEN
              plottype(ngs) = -55
            ELSE
              plottype(ngs) = -(ecount)
            ENDIF
          ENDIF

c          ndata(ngs) = 0
          expndata = 0
          CALL load_expt(index,expxdata,expydata,expndata,
     .                   MAXEXPNDATA,0,MAXDATX,dummy,colselect)
c          CALL load_expt(index,xdata1(1,ngs),ydata1(1,ngs),ndata(ngs),
c     .                   MAXGXS,0,MAXDATX,dummy,colselect)

c...      Only keep data in plot range:
          ndata(ngs) = 0
          DO i1 = 1, expndata
            IF (expxdata(i1).GE.xrange1.AND.
     .          expxdata(i1).LE.xrange2) THEN
              ndata(ngs) = ndata(ngs) + 1
              xdata1(ndata(ngs),ngs) = expxdata(i1)
              ydata1(ndata(ngs),ngs) = expydata(i1)
c              WRITE(0,*) 'KEEP:'
            ELSE 
c              WRITE(0,*) 'CHUCK:',expxdata(i1),xrange1,xrange2
            ENDIF
          ENDDO


c          elabs(ngs) = '    experiment'
          elabs(ngs) = dummy(1:20)
          CALL LoadArray(xdata,ndata1,xdata1(1,ngs),1,ndata(ngs))
c...      Scale experimental data:
          DO i1 = 1, ndata(ngs)
            ydata1(i1,ngs) = ydata1(i1,ngs) * scale
c...        Make sure -ve density values don't mess things up:
            IF (iopt.EQ.19.AND.ydata1(i1,ngs).LE.0.0) ydata1(i1,ngs)=LO
          ENDDO

c...HARDCODED NASTY!
          IF (index.EQ.24.OR.index.EQ.37) THEN
            WRITE(0,*) '*********************'
            WRITE(0,*) ' NASTY MACH HARDCODE! '
            WRITE(0,*) '*********************'
            DO i1 = 1, ndata(ngs)
              ydata1(i1,ngs) = 0.4*LOG(ydata1(i1,ngs))
            ENDDO
          ELSEIF (iopt.EQ.18.OR.iopt.EQ.19.OR.iopt.EQ.20) THEN
            IF (index.EQ.55) plottype(ngs) = 2 
            IF (index.EQ.56) plottype(ngs) = 2 
          ENDIF
          IF (index.EQ.26) plottype(ngs) = -1
          IF (index.EQ.28) plottype(ngs) = -3

c...      Check for more experimental data to plot:
          GOTO 75
        ELSE
c...      Calculate averages of experimental data:

          IF (.NOT..TRUE.) THEN
            WRITE(6,*) 'PLOT AVERAGE FOR 986',iopt

            deltapsin = 0.018   !  Core Thomson 
c            deltapsin = 0.010   ! RCP

            i1 = 0

c            DO psin = -2.0*deltapsin+1.0, 1.34, deltapsin  ! RCP
            DO psin = -0.196, 1.50, deltapsin              ! Core Thomson
              i1 = i1 + 1
              IF (i1.GT.MAXGXS)  
     .          CALL ER('PLOT986','xAVG array bounds error',*99)
              xavg(i1) = psin + 0.5 * deltapsin
              yavg(i1) = 0.0
              navg = 0
              DO i2 = ngs-ecount+1, ngs
                DO i3 = 1, ndata(i2)
                  IF (ydata1(i3,i2).NE.LO.AND.
     .                ABS(xdata1(i3,i2)-xavg(i1)).LT.0.5*deltapsin) THEN
                    navg = navg + 1
                    yavg(i1) = yavg(i1) + ydata1(i3,i2) 
                  ENDIF
                ENDDO
              ENDDO
              IF (navg.GT.0) THEN
                yavg(i1) = yavg(i1) / REAL(navg)
                WRITE(0,*) 'DATA:',xavg(i1),yavg(i1),navg
                WRITE(6,*) i1,xavg(i1),yavg(i1),navg
              ELSE
                i1 = i1 - 1
              ENDIF
            ENDDO 

c...        Plot average:
            ngs = ngs + 1
            ndata(ngs) = i1
            elabs(ngs) = '    AVERAGE'
            DO i1 = 1, ndata(ngs)
              xdata1(i1,ngs) = xavg(i1)
              ydata1(i1,ngs) = yavg(i1)
            ENDDO
            plottype(ngs) = 3 
            CALL LoadArray(xdata,ndata1,xdata1(1,ngs),1,ndata(ngs))
          ENDIF

          IF (.NOT..TRUE.) THEN

            WRITE(0,*) 
            WRITE(0,*) '***********************'
            WRITE(0,*) ' SWAPING ARRAYS IN 986'
            WRITE(0,*) '***********************'
            WRITE(0,*) 

            elabs(MAXNGS) = elabs(1)
            ndata(MAXNGS) = ndata(1)
            plottype(MAXNGS) = plottype(1)
            DO i1 = 1, ndata(MAXNGS)   
              xdata1(i1,MAXNGS) = xdata1(i1,1)
              ydata1(i1,MAXNGS) = ydata1(i1,1)
            ENDDO

            DO i2 = 1, ngs-1
              elabs(i2) = elabs(i2+1)
              ndata(i2) = ndata(i2+1)
              plottype(i2) = plottype(i2+1)
              DO i1 = 1, ndata(i2)   
                xdata1(i1,i2) = xdata1(i1,i2+1)
                ydata1(i1,i2) = ydata1(i1,i2+1)
              ENDDO
            ENDDO

            elabs(ngs) = '    OSM'
            ndata(ngs) = ndata(MAXNGS)
            plottype(ngs) = plottype(MAXNGS)
            DO i1 = 1, ndata(ngs)   
              xdata1(i1,ngs) = xdata1(i1,MAXNGS)
              ydata1(i1,ngs) = ydata1(i1,MAXNGS)
            ENDDO

c *HACK FOR PLOTTING THOMSON AND RCP DATA*
            IF (.FALSE..AND.ecount.EQ.1) THEN
              ngs = ngs - 1           
              plottype(1) = -55
            ELSEIF (.NOT..TRUE.) THEN
              ngs = ngs - 1

c              plottype(2) = 3

c              DO i1 = 1, ngs
c                plottype(i1) = i1 
c              ENDDO
c              plottype(1) = 55

              IF (.FALSE.) THEN
                ngs = ngs + 1
                elabs(ngs) = '    Thomson PSIn +0.015'
                ndata(ngs) = ndata(ngs-1)
                plottype(ngs) = 4
                DO i1 = 1, ndata(ngs)   
                  xdata1(i1,ngs) = xdata1(i1,ngs-1) + 0.015
                  ydata1(i1,ngs) = ydata1(i1,ngs-1)
                ENDDO
              ENDIF

c...          DIII-D upstream data shifting:
              IF (.NOT..TRUE..AND.
     .            machine.EQ.DIIID.AND.iopt.GE.18.AND.iopt.LE.20) THEN
                WRITE(0,*)
                WRITE(0,*) '*** SHIFTING UPSTREAM DIII-D DATA ***'
                WRITE(0,*)
                i2 = 1
                DO i1 = 1, ndata(i2)
                  xdata1(i1,i2) = xdata1(i1,i2) + 0.00  !0.01
c                  xdata1(i1,i2) = xdata1(i1,i2) + 0.014 ! ?
c                  xdata1(i1,i2) = xdata1(i1,i2) + 0.018 ! Te match
                ENDDO
                i2 = 2
                DO i1 = 1, ndata(i2)
c                  xdata1(i1,i2) = xdata1(i1,i2) + 0.002
                ENDDO
              ENDIF

c...          C-Mod upstream data shifting:
              IF (.FALSE..AND.
     .            machine.EQ.CMOD.AND.iopt.GE.18.AND.iopt.LE.20) THEN
                WRITE(0,*)
                WRITE(0,*) '*** SHIFTING UPSTREAM C-MOD DATA ***'
                WRITE(0,*)
                i2 = 1
                DO i1 = 1, ndata(i2)
                  xdata1(i1,i2) = xdata1(i1,i2) - 0.002
                ENDDO
                i2 = 2
                DO i1 = 1, ndata(i2)
c                  xdata1(i1,i2) = xdata1(i1,i2) + 0.002
                ENDDO
              ENDIF


c...          Estimate Ti:
              IF (.NOT..TRUE..AND.index.EQ.70.AND.ngs.EQ.2) THEN

                IF (iopt.EQ.18) THEN
c...              Load Thomson Te data:
                  tirat = 2.0
                  expndata = 0
                  CALL load_expt(77,expxdata,expydata,expndata,
     .                           MAXEXPNDATA,0,MAXDATX,dummy,2)
                  ntmp = 0
                  DO i1 = 1, expndata
                    IF (expxdata(i1).GE.xrange1.AND.
     .                  expxdata(i1).LE.xrange2) THEN
                      ntmp = ntmp + 1
                      xtmp(ntmp) = expxdata(i1) + 0.00 !0.01
                      ytmp(ntmp) = expydata(i1)
                    ENDIF
                  ENDDO
                  WRITE(0,*) ntmp,ndata(1)
                  WRITE(0,*) xrange1,xrange2
                  DO i1 = 1, ntmp
                    WRITE(0,*) i1,xtmp(i1),ytmp(i1),xdata1(i1,1)
                  ENDDO
                  IF (ntmp.NE.ndata(1))
     .              CALL ER('986:19','Te data not matching up-core',*99)

                  DO i1 = 1, ndata(1)   
                    IF (xtmp(i1).NE.xdata1(i1,1))
     .                CALL ER('986:19','X data no matching up-core',*99)

                    IF (xtmp(i1).GT.1.20) CYCLE

                    ydata1(i1,1) = ydata1(i1,1) * 
     .                (ytmp(i1) * (1.0 + tirat)) / 
     .                (ytmp(i1) * (1.0 + 1.0  ))
                  ENDDO
c...              Load RCP Te data:
                  tirat = 2.0
                  expndata = 0
                  CALL load_expt(index,expxdata,expydata,expndata,
     .                           MAXEXPNDATA,0,MAXDATX,dummy,2)
                  ntmp = 0
                  DO i1 = 1, expndata
                    IF (expxdata(i1).GE.xrange1.AND.
     .                  expxdata(i1).LE.xrange2) THEN
                      ntmp = ntmp + 1
                      xtmp(ntmp) = expxdata(i1)
                      ytmp(ntmp) = expydata(i1)
                    ENDIF
                  ENDDO
                  IF (ntmp.NE.ndata(ngs))
     .              CALL ER('986:19','Te data not matching up',*99)

                  DO i1 = 1, ndata(ngs)   
                    IF (xtmp(i1).NE.xdata1(i1,ngs))
     .                CALL ER('986:19','X data not matching up',*99)

                    IF (xtmp(i1).GT.1.20) CYCLE

                    ydata1(i1,ngs) = ydata1(i1,ngs) * 
     .                SQRT(ytmp(i1) * (1.0 + tirat)) / 
     .                SQRT(ytmp(i1) * (1.0 + 1.0  ))
                  ENDDO

                ENDIF

                IF (iopt.EQ.19) THEN
c...              Load RCP Te data:
                  tirat = 2.0
                  expndata = 0
                  CALL load_expt(index,expxdata,expydata,expndata,
     .                           MAXEXPNDATA,0,MAXDATX,dummy,2)
                  ntmp = 0
                  DO i1 = 1, expndata
                    IF (expxdata(i1).GE.xrange1.AND.
     .                  expxdata(i1).LE.xrange2) THEN
                      ntmp = ntmp + 1
                      xtmp(ntmp) = expxdata(i1)
                      ytmp(ntmp) = expydata(i1)
                    ENDIF
                  ENDDO
                  IF (ntmp.NE.ndata(ngs))
     .              CALL ER('986:19','Te data not matching up',*99)

                  DO i1 = 1, ndata(ngs)   

                    IF (xtmp(i1).NE.xdata1(i1,ngs))
     .                CALL ER('986:19','X data not matching up',*99)

                    IF (xtmp(i1).GT.1.20) CYCLE

                    ydata1(i1,ngs) = ydata1(i1,ngs) * 
     .                SQRT(ytmp(i1) * (1.0 + 1.0  )) / 
     .                SQRT(ytmp(i1) * (1.0 + tirat))
                  ENDDO

                ENDIF
                IF (iopt.EQ.20) THEN

                ENDIF
              ENDIF


c...          Reset x-axis data:
              ndata1 = 0
              DO i1 = 1, ngs
                CALL LoadArray(xdata,ndata1,xdata1(1,i1),1,ndata(i1))
              ENDDO

            ENDIF

          ENDIF



          BACKSPACE 5
        ENDIF


        DO i1 = 1, ngs
          CALL MapArray(xdata       ,ydata (1,i1),1,ndata1   ,
     .                  xdata1(1,i1),ydata1(1,i1),1,ndata(i1))          
        ENDDO


      ENDIF

      REF = '                                   '
      JOB = '                                   '

c...
c
c
c
c
      IF     (plotmode.EQ.1) THEN
c...    Line plot:

        slopt2 = 2
        DO i1 = 1, 9
          plottype(i1) = i1 
        ENDDO

c          CALL DRAW(xdata,TWIDS,ydata,MAXgXs,ndata1,ANLY,ngs,
c     >              100,xrange1,xrange2,-HI,HI,IGNORS,0,AVS,0,
c     >              JOB,TITLE,XLAB,YLAB,ELABS,REF,VIEW,
c     >              PLANE,TABLE,1,1,1.0,0)

        IF (iopt.EQ.4) THEN

          WRITE(0,*) 'XADATA:',(xdata(i1),i1=1,ndata(1))

          IF (xrange1.EQ.HI) xrange1 = xdata(1)
          IF (xrange2.EQ.HI) xrange2 = xdata(ndata(1))

          CALL DRAW(xdata,TWIDS,ydata,MAXgXs,ndata(1),ANLY,ngs,
     >              100,xrange1,xrange2,-HI,HI,IGNORS,0,AVS,0,
     >              JOB,TITLE,XLAB,YLAB,ELABS,REF,VIEW,
     >              PLANE,TABLE,1,1,1.0,0)

c...      Annotate graph:
          CALL PSPACE (map1x,map2x,map1y,map2y)
          CALL MAP    (cxmin,cxmax,cymin,cymax)
          CALL BROKEN(6,6,6,6)
          CALL LinCol(1)

            tag(1) = 'separatrix'
            tag(2) = 'inner nose'
            tag(3) = 'outer nose'
            tag(4) = 'floor'

            idum2 = 20

c...        Separatrix:
            CALL POSITN (0.0,cymin)
            CALL JOIN   (0.0,cymax)        
            CALL CTRORI(90.0)
            CALL CTRMAG(10)
            CALL PLOTST(0.0+0.015*(cxmax-cxmin),
     .                  cymin+0.05*(cymax-cymin),tag(1))
c...        Outer nose:
            CALL POSITN(rho(idum2,OUT23)*1.0E+03,cymin)
            CALL JOIN  (rho(idum2,OUT23)*1.0E+03,cymax)        
            CALL PLOTST(rho(idum2,OUT23)*1.0E+03+0.015*(cxmax-cxmin),
     .                cymin+0.05*(cymax-cymin),tag(3))
c...        Bottom of vertical target (HARDCODED!):
            CALL POSITN(rho(29,OUT23)*1.0E+03,cymin)
            CALL JOIN  (rho(29,OUT23)*1.0E+03,cymax)        
            CALL PLOTST(rho(29,OUT23)*1.0E+03-0.010*(cxmax-cxmin),
     .                cymin+0.05*(cymax-cymin),tag(4))

        ELSE

          DO i1 = 1, ndata1
            WRITE(6,'(A,I6,1P,E10.2,3X,50(E10.2))') 
     .        'LINE PLOT: ',i1,xdata(i1),(ydata(i1,i2),i2=1,ngs)
          ENDDO

          CALL DRAW(xdata,TWIDS,ydata,MAXgXs,ndata(1),ANLY,ngs,
     >              100,0.5,REAL(ndata(1))+0.5,-HI,HI,IGNORS,0,AVS,0,
     >              JOB,TITLE,XLAB,YLAB,ELABS,REF,VIEW,
     >              PLANE,TABLE,1,1,1.0,0)

c...      Draw lablels on x-axis:
          CALL PSPACE(MAP1X,MAP2X,MAP1Y-0.05,MAP1Y)
          CALL MAP(0.5,REAL(ndata)+0.5,0.0,1.0)
          DO i1 = 1, ndata(1)
            CALL PCSCEN(REAL(i1),map2y-0.04,
     .                  pnames1(i1)(1:LEN_TRIM(pnames1(i1))))
          ENDDO
        ENDIF

      ELSEIF (plotmode.EQ.2.OR.plotmode.EQ.4) THEN
c...    Bar chart:

        CALL RGB
        CALL ColSet(1.0,1.0,1.0,255)
        CALL FilCol(255)

c...HACK ALERT!
        IF (probescale.EQ.1) THEN
c...      Record scaling data:
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              probescaledat(i2,i1) = ydata(i2,i1)
            ENDDO
          ENDDO
          probescale = 2
        ELSEIF (probescale.EQ.2) THEN
c...      Record scaling data:
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              probescaledat2(i2,i1) = ydata(i2,i1)
            ENDDO
          ENDDO
          probescale = 3
        ELSEIF (probescale.EQ.3) THEN
c...      Scale data:
          ngs = ngsydata2(1)
          DO i1 = 1, ngs
          
            ndata(i1) = ndataydata2(i1,1)          

            DO i2 = 2, ndata(i1)
c              WRITE(0,*) 'DATA:',i1,i2,probescaledat (i2,i1),
c     .                                 probescaledat2(i2,i1)
              IF (probescaledat2(i2,i1).LE.LO) CYCLE
              ydata(i2,i1) = probescaledat(i2,i1) / 
     .                       probescaledat2(i2,i1)
c              ydata(i2,i1) = ydata(i2,i1) * probescaledat(i2,i1) /
c     .                                      probescaledat(i2,1)
            ENDDO
          ENDDO

          maxy = -HI
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              maxy = MAX(maxy,ydata(i2,i1))
              WRITE(0,*) 'MAX:',ydata(i2,i1),maxy
            ENDDO
          ENDDO
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              ydata(i2,i1) = ydata(i2,i1) / maxy 
              WRITE(0,*) 'MAX:',ydata(i2,i1),maxy
            ENDDO
          ENDDO

          probescale = 0
        ENDIF

c...    Fudge the data to include n-n collisions:
        IF (fudge.EQ.1) THEN
          WRITE(0,*) 'FUDGECICLE!'
          DO i1 = 1, ngs 
            rdum1 = ydata(1,i1)
            ydata(1,i1) = 0.718*ydata(1,i1)**1.3832
            DO i2 = 2, ndata(i1)
              ydata(i2,i1) = ydata(i2,i1) / rdum1 * ydata(1,i1)
            ENDDO
          ENDDO

        ELSEIF (fudge.EQ.7) THEN
          WRITE(0,*) 'FUDGECICLE! 7'
          DO i1 = 1, ngs 
            DO i2 = 1, ndata(i1)
              ydata(i2,i1) = 0.718*ydata(i2,i1)**1.3832
            ENDDO
          ENDDO

        ELSEIF (fudge.EQ.4) THEN
c...      Calculate the momentum flux into the plenum (lots of assumptions
c         about plots that came before):

          ngs = ngsydata2(3)
          DO i1 = 1, ngs 
            ndata(i1) = ndataydata2(i1,3)
            DO i2 = 1, ndata(i1)
              mass = crmb * REAL(i1) * amu


c...          THIS IS WRONG!
              vavg = 0.9*SQRT(2.0 * ydata2(i2,i1,3) * ECH / mass)
c              area = PI * 2.0 * 0.63 * 0.0233
              area = 1.1 * PI * 2.0 * 0.639 * 0.026
              fact = 760.0E+03 / 101.3E+03 
              flux = ydata2(i2,i1,2) * 1.0E+21 / REAL(i1)

c              WRITE(0,'(A,1P,5E10.2)') 
c     .          'DATA:',mass,vavg,area,fact,flux

              ydata(i2,i1) = flux * mass * vavg * fact / area
c              ydata(i2,i1) = 2.0/3.0*flux * mass * vavg * fact / area



c... NEW
              flux = ydata2(i2,i1,2) * 1.0E+20 / REAL(i1)
              area = 0.1 * PI * 2.0 * 0.639 * 0.026

              IF (i1.EQ.2) THEN
                ydata(i2,i1) = (flux / area) * 4.0 * SQRT(PI / 8.0) * 
     .                         SQRT(mass) * 
     .                         SQRT(0.67 * 0.5 * ydata2(i2,i1,3) * ECH)
     .                         * fact
              ELSE
                ydata(i2,i1) = (flux / area) * 4.0 * SQRT(PI / 8.0) * 
     .                         SQRT(mass) * 
     .                         SQRT(0.67 * ydata2(i2,i1,3) * ECH)
     .                         * fact
              ENDIF

c              ydata(i2,i1) = (flux / area) * 4.0 * SQRT(PI / 8.0) * 
c     .                       SQRT(mass) * 
c     .                       SQRT(0.67 * 0.038 * ECH) * fact


              WRITE(0,*) 'MOM FLUX:',ydata(i2,i1),ydata2(i2,i1,3),mass

            ENDDO
          ENDDO

        ELSEIF (fudge.EQ.5) THEN
c...      Calculate the momentum flux into the plenum (lots of assumptions
c         about plots that came before):

          ngs = ngsydata2(1)
          DO i1 = 1, ngs 
            ndata(i1) = ndataydata2(i1,1)
            DO i2 = 1, ndata(i1)
              ydata(i2,i1) = ydata2(i2,i1,1) / 
     .                      (ydata2(i2,1,4) + ydata2(i2,2,4))
          WRITE(0,*) 'MOM FRACTION:',ydata(i2,i1)
            ENDDO
          ENDDO



        ELSEIF (fudge.EQ.6) THEN
c...      Calculate the leakage flux:
          ngs = 2
          ndata(1) = ndataydata2(1,1)
          ndata(2) = ndataydata2(1,1)

          DO i2 = 1, ndata(1)
            ydata(i2,1) = 0.0
            ydata(i2,2) = 0.0          
          ENDDO
          totalleakage1 = 0.0
          totalleakage2 = 0.0
          totalleakage3 = 0.0
          totalleakage4 = 0.0
          DO i1 = 1, ngsydata2(1) 
c          DO i1 = 2, ngsydata2(1) 
            DO i2 = 1, ndata(1)
              ydata(i2,1) = ydata(i2,1) + ydata2(i2,i1,1) * scale
              ydata(i2,2) = ydata(i2,2) + ydata2(i2,i1,2) * scale
              totalleakage1 = totalleakage1 + ydata2(i2,i1,1) * scale
              totalleakage2 = totalleakage2 + ydata2(i2,i1,2) * scale
              IF (i1.GT.2) THEN
c...            Flux through the outer divertor substructure:              
                totalleakage3 = totalleakage3 + ydata2(i2,i1,1) * scale
                totalleakage4 = totalleakage4 + ydata2(i2,i1,2) * scale
              ENDIF
            ENDDO
          ENDDO
       
          WRITE(0,*) 'TOTAL LEAK 1=',totalleakage1/10.0,
     .                               totalleakage1/scale
          WRITE(0,*) 'TOTAL LEAK 2=',totalleakage2/10.0,
     .                               totalleakage2/scale
          WRITE(0,*) 'TOTAL LEAK 3=',totalleakage3/10.0,
     .                               totalleakage3/scale
          WRITE(0,*) 'TOTAL LEAK 4=',totalleakage4/10.0,
     .                               totalleakage4/scale
        ENDIF

c...    Store plotted data for use in future processing:
        nydata2 = nydata2 + 1
        ngsydata2(nydata2) = ngs
        DO i1 = 1, ngsydata2(nydata2) 
          ndataydata2(i1,nydata2) = ndata(i1)
          DO i2 = 1, ndataydata2(i1,nydata2)
c            WRITE(0,*) 'DA:',i2,i1,ydata(i2,i1)
            ydata2(i2,i1,nydata2) = ydata(i2,i1)
          ENDDO
        ENDDO

        IF (plotmode.EQ.4) THEN
c...      Stacked bar chart:
          iopt_ghost = 2
          maxy = -HI
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              maxy = MAX(maxy,ydata(i2,i1))
            ENDDO
          ENDDO
          maxy = maxy * 2.0
        ELSE
          iopt_ghost = 3
          maxy = -HI
          DO i1 = 1, ngs
            DO i2 = 1, ndata(i1)
              maxy = MAX(maxy,ydata(i2,i1))
            ENDDO
          ENDDO
        maxy = maxy * 1.10
        ENDIF

        READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'Maxy'.OR.dummy(8:11).EQ.'maxy'.OR.
     .        dummy(8:11).EQ.'MAXY') THEN
          READ(dummy,*) cdum1,maxy
        ELSE
          BACKSPACE 5
        ENDIF

        CALL GrBar(ydata,1,ndata(1),MAXGXS,ngs,0.0,maxy,5,
     .             pnames1,pnames2,
     .             title,cnames,ylab)

        READ(5,'(A5000)') dummy
        IF   (dummy(8:11).EQ.'Line'.OR.dummy(8:11).EQ.'line'.OR.
     .        dummy(8:11).EQ.'LINE') THEN
c...      Annotate graph:
          CALL PSPACE (map1x,map2x,map1y,map2y)
          CALL MAP    (0.0,1.0,0.0,maxy)
          CALL BROKEN(6,6,6,6)
          CALL LinCol(1)
c...      Flow reversal line:
          CALL POSITN (0.0,1.0)
          CALL JOIN   (1.0,1.0)        
        ELSE
          BACKSPACE 5
        ENDIF

c...    Finish off the plot:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP(0.0,1.0,0.0,1.0)
        CALL FULL
        CALL POSITN (0.0,1.0)
        CALL JOIN   (1.0,1.0)
        CALL JOIN   (1.0,0.0)
        CALL JOIN   (0.0,0.0)
        CALL JOIN   (0.0,1.0)

      ELSEIF (plotmode.EQ.3) THEN
c...    X-Y scatter:
        slopt4 = 1 
c        plottype(1) = 55

        IF (xrange1.EQ.HI) xrange1 = xdata(1)
        IF (xrange2.EQ.HI) xrange2 = xdata(ndata1)

        IF (cumulative) THEN
c...      Cumulative plot:
          slopt2 = 2
          DO i1 = 1, 9
            plottype(i1+1) = i1 + 1
          ENDDO
          DO i1 = 2, ngs
            DO i2 = 1, ndata1
              ydata(i2,i1) = ydata(i2,i1) + ydata(i2,i1-1)
            ENDDO
          ENDDO
 
          DO i1 = 1, ndata1
            WRITE(6,'(A,I6,1P,E10.2,3X,50(E10.2))') 
     .        'X-Y PLOT: ',i1,xdata(i1),(ydata(i1,i2),i2=1,ngs)
          ENDDO

          CALL DRAW(xdata,TWIDS,ydata,MAXgXs,ndata1,ANLY,ngs,
     >              100,xrange1,xrange2,-HI,HI,IGNORS,0,AVS,0,
     >              JOB,TITLE,XLAB,YLAB,ELABS,REF,VIEW,
     >              PLANE,TABLE,8,1,1.0,0)
        ELSE
c...      Standard line plot:

          DO i1 = 1, ndata1
            WRITE(6,'(A,I6,1P,E16.8,3X,50(E10.2))') 
     .        'X-Y PLOT: ',i1,xdata(i1),(ydata(i1,i2),i2=1,ngs)
          ENDDO

          CALL DRAW(xdata,TWIDS,ydata,MAXgXs,ndata1,ANLY,ngs,
     >              100,xrange1,xrange2,yrange1,yrange2,IGNORS,0,AVS,0,
c     >              100,xrange1,xrange2,-HI,HI,IGNORS,0,AVS,0,
     >              JOB,TITLE,XLAB,YLAB,ELABS,REF,VIEW,
     >              PLANE,TABLE,1,1,1.0,0)

        ENDIF

47      READ(5,'(A5000)') dummy
        IF   (dummy(8:11).EQ.'Line'.OR.dummy(8:11).EQ.'line'.OR.
     .        dummy(8:11).EQ.'LINE') THEN

          READ (dummy,*) cdum1,anotemode,idum1,idum2

c...      Annotate graph:
          CALL PSPACE (map1x,map2x,map1y,map2y)
          CALL MAP    (cxmin,cxmax,cymin,cymax)
          CALL BROKEN(6,6,6,6)
          CALL LinCol(1)

          IF (idum1.EQ.1) THEN
            tag(1) = 'separatrix'
            tag(2) = 'inner nose'
            tag(3) = 'outer nose'
            tag(4) = 'floor'
          ELSE
            tag(1) = '                                '
            tag(2) = '                                '
            tag(3) = '                                '
            tag(4) = '                                '
          ENDIF

          IF (anotemode.EQ.1) THEN
c...        Separatrix:
            CALL POSITN (0.0,cymin)
            CALL JOIN   (0.0,cymax)        
            CALL CTRORI(90.0)
            CALL CTRMAG(10)
            CALL PLOTST(0.0+0.015*(cxmax-cxmin),
     .                  cymin+0.05*(cymax-cymin),tag(1))
c...        Inner nose:
            CALL POSITN(rho(idum2,OUT23)*1.0E+03,cymin)
            CALL JOIN  (rho(idum2,OUT23)*1.0E+03,cymax)        
            CALL PLOTST(rho(idum2,OUT23)*1.0E+03+0.015*(cxmax-cxmin),
     .                cymin+0.05*(cymax-cymin),tag(2))
          ELSEIF (anotemode.EQ.2) THEN
c...        Separatrix:
            CALL POSITN (0.0,cymin)
            CALL JOIN   (0.0,cymax)        
            CALL CTRORI(90.0)
            CALL CTRMAG(10)
            CALL PLOTST(0.0+0.015*(cxmax-cxmin),
     .                  cymin+0.05*(cymax-cymin),tag(1))
c...        Outer nose:
            CALL POSITN(rho(idum2,OUT23)*1.0E+03,cymin)
            CALL JOIN  (rho(idum2,OUT23)*1.0E+03,cymax)        
            CALL PLOTST(rho(idum2,OUT23)*1.0E+03+0.015*(cxmax-cxmin),
     .                cymin+0.05*(cymax-cymin),tag(3))
c...        Bottom of vertical target (HARDCODED!):
            CALL POSITN(rho(29,OUT23)*1.0E+03,cymin)
            CALL JOIN  (rho(29,OUT23)*1.0E+03,cymax)        
            CALL PLOTST(rho(29,OUT23)*1.0E+03-0.010*(cxmax-cxmin),
     .                cymin+0.05*(cymax-cymin),tag(4))
          ELSEIF (anotemode.EQ.3) THEN

c...        Flow reversal line:
            CALL POSITN (cxmin,1.0)
            CALL JOIN   (cxmax,1.0)        

          ELSEIF (anotemode.EQ.4) THEN

c...        Zero line:
            CALL POSITN (cxmin,0.0)
            CALL JOIN   (cxmax,0.0)        

          ELSE
            CALL ER('986','Invalid ANOTEMODE value',*99)
          ENDIF
          CALL CTRORI (0.0)

c...      See if there are any more lines to be plotted:
          GOTO 47
        ELSE
          BACKSPACE 5
        ENDIF

c...    Finish off the plot:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP(0.0,1.0,0.0,1.0)
        CALL FULL
        CALL POSITN (0.0,1.0)
        CALL JOIN   (1.0,1.0)
        CALL JOIN   (1.0,0.0)
        CALL JOIN   (0.0,0.0)
        CALL JOIN   (0.0,1.0)

      ELSE
        CALL ER('986','Invalid mode',*99)
      ENDIF

c...  Add a comment to the plot:
49    READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Cmnt'.OR.dummy(8:11).EQ.'cmnt'.OR.
     .      dummy(8:11).EQ.'CMNT') THEN
      
        READ (dummy,*) cdum1,xpos,ypos,size,caption
      
c...    Annotate graph:
        CALL PSPACE (map1x,map2x,map1y,map2y)
        CALL MAP    (0.0,1.0,0.0,1.0)
        CALL LinCol(1)
        CALL CTRMAG(size)
        CALL PLOTST(xpos,ypos,caption(1:LEN_TRIM(caption)))
c...    Another comment:        
        GOTO 49
      ELSE
        BACKSPACE 5
      ENDIF



c...  Add a caption to the plot:
      READ(5,'(A5000)') dummy
      IF   (dummy(8:11).EQ.'Note'.OR.dummy(8:11).EQ.'note'.OR.
     .      dummy(8:11).EQ.'NOTE') THEN
        READ(dummy,*) cdum1,xpos,ypos,size,caption
        CALL AddCaption(caption,xpos,ypos,size)
      ELSE
        BACKSPACE 5
      ENDIF


      IF (plotframe.EQ.1.OR.plotframe.EQ.3) THEN
        nydata2 = 0
        CALL RSet(ydata2,MAXGXS*MAXNGS*10,0.0)
        CALL FRAME
      ENDIF


      RETURN
99    CONTINUE
      WRITE(6,*) 'PLOT OPTION=',iopt
      STOP
      END


      SUBROUTINE AddCaption(buffer,xpos,ypos,size)
      IMPLICIT none

      INTEGER size,lenbuf,i2,i3,linepos,indent,maxlen
      REAL xpos,ypos,xpos1,ypos1
      CHARACTER*(*) buffer
      CHARACTER*1024 line,sp,word

      CALL PSPACE(0.0, 1.35, 0.0, 1.0)
      CALL MAP(0.0, 1.35, 0.0, 1.0)

      WRITE(sp,'(1024X)')

      lenbuf = LEN_TRIM(buffer)
      
      buffer(lenbuf+1:lenbuf+1) = ' '
      
      lenbuf  = lenbuf + 1
      i3      = 0
      linepos = 0
      indent  = 0      
      maxlen  = NINT((1.35 - xpos) * 1100.0) / size

      CALL CTRMAG(size)

      WRITE(line,'(1024X)')
      
      xpos1 = xpos
      ypos1 = ypos

      DO i2 = 1, lenbuf
        IF (buffer(i2:i2).EQ.' ') THEN

          IF (i2.EQ.lenbuf) THEN
            CALL PLOTST(xpos1,ypos1,line(1:LEN_TRIM(line))//' '//
     .                              word(1:LEN_TRIM(word)))

          ELSEIF (linepos+i3+1.GT.maxlen) THEN
            IF (i3.LT.maxlen-indent) THEN
              CALL PLOTST(xpos1,ypos1,line(1:LEN_TRIM(line)))
              ypos1 = ypos1 - 0.015 * 0.1 * REAL(size)
              WRITE(line,'(1024X)')
              WRITE(line,'(A)') sp(1:indent)
              linepos = indent
            ELSE
              CALL ER('GetInfo','Word too long',*99)
            ENDIF
          ENDIF

          WRITE(line(linepos+1:1024),'(A)') word(1:i3)//' '
          WRITE(word,'(1024X)')
          linepos = linepos + i3 + 1
          i3 = 0
        ELSE
          i3 = i3 + 1
          word(i3:i3) = buffer(i2:i2)
        ENDIF
      ENDDO      

      RETURN
99    STOP
      END
c
c
c
c
c
c
c
c
      SUBROUTINE SnatchData(file,mark,pattern,loc1,loc2,type,ndata,
     .                      vdata)
      IMPLICIT none

      INTEGER type,ndata,loc1,loc2
      REAL    vdata(*)
      CHARACTER*(*) pattern,mark
      CHARACTER*128 file

      INCLUDE 'params'

      INTEGER MAXFILES
      PARAMETER (MAXFILES = 100)


      INTEGER retcode,fp,i1,i2,len1,len2,idum1
      LOGICAL foundmark
      CHARACTER command*256,dummy*256
      CHARACTER*128 fname(MAXFILES)


c...  Return nil:
      IF (pattern(1:3).EQ.'Nil'.OR.pattern(1:3).EQ.'NIL'.OR.
     .    pattern(1:3).EQ.'nil') THEN
        ndata = 1
        vdata(1) = 0.0
        RETURN
      ENDIF        

c...  Return data in the input file:
      IF (file(1:3).EQ.'Dir'.OR.file(1:3).EQ.'DIR'.OR.
     .    file(1:3).EQ.'dir') THEN
        ndata = 1
        IF (type.EQ.1) THEN
c...      REAL:
          READ(pattern,*) vdata(1)
        ELSE
          CALL ER('SnatchData','Unsupported data type',*99)
        ENDIF
        RETURN
      ENDIF        



c...  Check if data files for the current case are to be searched:
      WRITE(6,*) '>'//file(1:9)//'<'  
      IF (file(1:9).EQ.'<current>') THEN
        WRITE(dummy,'(256X)')
        WRITE(dummy,'(A)') file(10:LEN_TRIM(file))
        WRITE(file,'(128X)')
        CALL GetEnv('CASENAME',file)

        ndata = 1
        fname(ndata) = '/home/steven/divimp/results/'//
     .                 file(1:LEN_TRIM(file))//
     .                 dummy(1:LEN_TRIM(dummy))
        WRITE(6,*) 'FNAME:',fname(ndata)(1:LEN_TRIM(fname(ndata)))
      ELSE
c...    Get list of files:
        WRITE(command,'(256X)')
        command = 'rm -rf temp.dat'
        CALL CISSUE(command,retcode)

        WRITE(command,'(256X)')
        command = 'ls ~/divimp/results/'//file(1:LEN_TRIM(file))//
     .            ' > temp.dat'
        CALL CISSUE(command,retcode)

        fp = 99
        OPEN(UNIT=fp,FILE='temp.dat',STATUS='OLD',ERR=99)
        ndata = 0
        DO WHILE (.TRUE.)
          READ(fp,'(A128)',END=10,ERR=99) fname(ndata+1)
          ndata = ndata + 1
          WRITE(6,*) 'FNAME:',fname(ndata)(1:LEN_TRIM(fname(ndata)))
        ENDDO
10      CONTINUE
      ENDIF

c...  Scan the file and extract the requested data:

      DO i1 = 1, ndata

        fp = 99
        OPEN(UNIT=fp,FILE=fname(i1),STATUS='OLD',ERR=98)

        vdata(i1) = LO

        len1 = LEN_TRIM(mark)  
	
        IF (len1.EQ.0) THEN
          len1 = LEN_TRIM(pattern)  
          foundmark = .TRUE.
        ELSE
          foundmark = .FALSE.
        ENDIF

        DO WHILE(.TRUE.)

          READ(fp,'(A256)',END=20,ERR=99) dummy

          len2 = LEN_TRIM(dummy)  

          IF (foundmark) THEN
c...        Mark found, search for data:

            DO i2 = 1, len2-len1+1
              IF (pattern(1:len1).EQ.dummy(i2:i2+len1-1)) THEN
                IF     (type.EQ.1) THEN
                  READ(dummy(loc1:loc2),*,ERR=99) vdata(i1)
                  WRITE(6,*) 'VDATA:',vdata(i1)
                ELSEIF (type.EQ.2) THEN
                  READ(dummy(loc1:loc2),*,ERR=99) idum1
                  vdata(i1) = REAL(idum1)
                ENDIF
                GOTO 20
              ENDIF
            ENDDO
          ELSE
c...        Search for mark:
            DO i2 = 1, len2-len1+1
              IF (mark(1:len1).EQ.dummy(i2:i2+len1-1)) THEN
                len1 = LEN_TRIM(pattern)  
                foundmark = .TRUE.
              ENDIF
            ENDDO

          ENDIF

        ENDDO
20      CONTINUE

        IF (vdata(i1).EQ.LO) THEN
          WRITE(6,*) 'PATTERN: >'//pattern//'<'
          CALL ER('SnatchData','Pattern not found',*99)
        ENDIF

        CLOSE(fp)

      ENDDO




      RETURN
98    CONTINUE
      WRITE(0,*) 'SnatchData: Cannot open data file >'//fname(i1)//'<'
      WRITE(6,*) 'SnatchData: Cannot open data file >'//fname(i1)//'<'
      STOP
99    CONTINUE
      WRITE(0,*) 'PATTERN:',fname(i1),pattern
      WRITE(0,*) 'SnatchData: Problem'
      STOP
      END
c
c ======================================================================
c
c
c

