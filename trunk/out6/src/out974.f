
c
c ======================================================================
c
c subroutine: DrawHalpha
c
      SUBROUTINE DrawHalpha(array,opt_lines)

      IMPLICIT none

      INTEGER array,opt_lines

      INCLUDE 'params'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'
      include 'comgra'
c
c      COMMON /COMGRA/ CXMIN,CXMAX,CYMIN,CYMAX,IPLOTS,ICOL,NPLOTS,ISPOT
c      REAL            CXMIN,CXMAX,CYMIN,CYMAX
c      INTEGER         IPLOTS,ICOL,NPLOTS,ISPOT
c
      REAL FARAWAY,CLOSEBY
c      PARAMETER (FARAWAY=1.0,CLOSEBY=0.03)
      PARAMETER (FARAWAY=100.0,CLOSEBY=0.03)

      INTEGER       id,in,in1,in2,i1,i2,count,num_lines,idum1
      REAL          deltar,deltaz,angle,ang_lines(100),dfact(4),
     .              xpos,ypos,xpos2,ypos2,linelen,siz
      integer len,lenstr
      external lenstr
      character*100 desc
      CHARACTER*100 cdum1

      REAL hal_r(NUMHAL),hal_z(NUMHAL)
      DATA (hal_r(in),in=1,NUMHAL) /1.173, 0.754, 1.953, 0.833/
      DATA (hal_z(in),in=1,NUMHAL) /0.000, 0.465, 0.050,-0.488/

      CHARACTER*20 hal_name(NUMHAL)
      DATA (hal_name(in),in=1,NUMHAL) /'B side straight     ',
     .                                 'B top straight      ',
     .                                 'F side              ',
     .                                 'K bottom straight   '/



      dfact(1) = 0.60
      dfact(2) = 1.07
      dfact(3) = 1.0
      dfact(4) = 1.0
c
c
c
      IF (array.EQ.NUMHAL+1) THEN
        in1 = 1
        in2 = NUMHAL
      ELSE
        in1 = array
        in2 = array
      ENDIF

      CALL THICK  (1)
      CALL CTRMAG (7)

c      CALL LINCOL(1)
      CALL PSPACE(map1x,map2x,map1y,map2y)
      CALL MAP   (cxmin,cxmax,cymin,cymax)


      DO in = in1, in2
c
c       Cross-hairs:
c

c        CALL FULL

c        DO angle = 0.0, 270.0, 90.0
c          deltar = COS(angle*PI/180.0) * CLOSEBY
c          deltaz = SIN(angle*PI/180.0) * CLOSEBY
c
c          CALL POSITN(hal_r(in),hal_z(in))
c          CALL LINE  (deltar,deltaz)
c        ENDDO

c
c       View bounds:
c
c        CALL BROKEN(3,3,3,3)
c
c        deltar = COS(hal_ang(1,in)*PI/180.0) * FARAWAY
c        deltaz = SIN(hal_ang(1,in)*PI/180.0) * FARAWAY
c
c        CALL POSITN(hal_r(in),hal_z(in))
c        CALL LINE  (deltar,deltaz)
c
c
c        id = hal_num(in)
c        deltar = COS(hal_ang(id,in)*PI/180.0) * FARAWAY
c        deltaz = SIN(hal_ang(id,in)*PI/180.0) * FARAWAY
c
c        CALL POSITN(hal_r(in),hal_z(in))
c        CALL LINE  (deltar,deltaz)
c
c       Draw peaks on view if requested:
c

        IF (opt_lines.EQ.1) THEN
          CALL BROKEN(6,6,6,6)
          CALL LINCOL(ncols+1)
          CALL LINCOL(1)

          READ(5,*) cdum1,num_lines
          BACKSPACE 5

          IF (num_lines.EQ.-2) THEN
            READ(5,*) cdum1,idum1,xpos,ypos,
     .                num_lines,(ang_lines(i1),i1=1,num_lines)

            DO i1 = 1, num_lines
              deltar = COS(ang_lines(i1)*PI/180.0)
              deltaz = SIN(ang_lines(i1)*PI/180.0)

              CALL POSITN(xpos,ypos)
              CALL LINE  (deltar*FARAWAY,deltaz*FARAWAY)

              CALL POSITN(xpos+deltar*dfact(in),
     .                    ypos+deltaz*dfact(in))
              CALL TYPENF(ang_lines(i1),1)
            ENDDO

          ELSEIF (num_lines.EQ.-5) THEN
            READ(5,*) cdum1,idum1,xpos,ypos,
     .                xpos2,ypos2,desc
c            DO i1 = 1, num_lines
c              deltar = COS(ang_lines(i1)*PI/180.0)
c              deltaz = SIN(ang_lines(i1)*PI/180.0)

              call full   
              call lincol(4)  
              CALL POSITN(xpos,ypos)
              CALL JOIN  (xpos2,ypos2)
              call lincol(1)
              CALL POSITN(xpos+CLOSEBY,
     .                    ypos)
              len = lenstr(desc)
              call ctrmag(16)
              CALL TYPECS(desc(1:len))
              call ctrmag(10) 
c            ENDDO

          elseIF (num_lines.EQ.-6) THEN
            READ(5,*) cdum1,idum1,xpos,ypos,linelen,
     .                num_lines,(ang_lines(i1),i1=1,num_lines),desc
            DO i1 = 1, num_lines
              deltar = COS(ang_lines(i1)*PI/180.0)
              deltaz = SIN(ang_lines(i1)*PI/180.0)

              call full   
              call lincol(3)
              CALL POSITN(xpos,ypos)
              CALL LINE  (deltar*linelen,deltaz*linelen)
              call lincol(1)
 
              CALL POSITN(xpos+deltar*dfact(in),
     .                    ypos+deltaz*dfact(in))
              call ctrmag(12)
c              CALL TYPENF(ang_lines(i1),1)
              CALL TYPENI(INT(ang_lines(i1)))
              call ctrmag(10)
            ENDDO
            call lincol(1)
            CALL POSITN(xpos+CLOSEBY,
     .                    ypos)
            len = lenstr(desc)
            call ctrmag(16)
            CALL TYPECS(desc(1:len))
            call ctrmag(10)
          ELSEIF (num_lines.EQ.-7) THEN
            desc = ' '  
            READ(5,*) cdum1,idum1,xpos,ypos,siz,desc
c            DO i1 = 1, num_lines
c              deltar = COS(ang_lines(i1)*PI/180.0)
c              deltaz = SIN(ang_lines(i1)*PI/180.0)
              call full  
              call lincol(6)  
              call filcol(6)   
              CALL POSITN(xpos,ypos)
              CALL circle(siz)
              call lincol(1)
              call filcol(0)   
              call ctrmag(16)
              len = lenstr(desc)
              if (len.gt.1) then 
                 call positn(xpos+CLOSEBY,ypos)
                 CALL TYPECS(desc(1:len))
              endif
              call ctrmag(10) 

          ELSEIF (num_lines.EQ.-3) THEN

            READ(5,*) cdum1,idum1,xpos,ypos,
     .                num_lines,(ang_lines(i1),i1=1,num_lines+1)

            DO i1 = 2, num_lines+1
              CALL POSITN(xpos,ypos)
              CALL JOIN  (ang_lines(i1),ang_lines(1))
            ENDDO

          ELSEIF (num_lines.EQ.-4) THEN

            READ(5,*) cdum1,idum1,xpos,ypos,
     .                num_lines,(ang_lines(i1),i1=1,num_lines)

            DO i1 = 1, num_lines
              deltar = FARAWAY * COS(ang_lines(i1)*PI/180.0)
              deltaz = FARAWAY * SIN(ang_lines(i1)*PI/180.0)

              CALL POSITN(xpos       ,ypos       )
              CALL JOIN  (xpos+deltar,ypos+deltaz)
            ENDDO

          ELSEIF (num_lines.EQ.-1) THEN
            READ(5,*) cdum1,num_lines,(ang_lines(i1),i1=1,num_lines)

            DO i1 = 11, id-11
              count = 0
              DO i2 = i1-10, i1+10
                IF (hal_val(i2,in).GT.hal_val(i1,in)) count = 1
              ENDDO

              IF (count.EQ.0) THEN
                deltar = COS(hal_ang(i1,in)*PI/180.0)
                deltaz = SIN(hal_ang(i1,in)*PI/180.0)

                CALL POSITN(hal_r(in),hal_z(in))
                CALL LINE  (deltar*FARAWAY,deltaz*FARAWAY)

                CALL POSITN(hal_r(in)+deltar,hal_z(in)+deltaz)
                CALL TYPENF(hal_ang(i1,in),1)
              ENDIF
            ENDDO
          ELSE
            READ(5,*) cdum1,num_lines,(ang_lines(i1),i1=1,num_lines)
c
c         Display chords to plot from DATA line:
c
            DO i1 = 1, num_lines

              deltar = COS(ang_lines(i1)*PI/180.0)
              deltaz = SIN(ang_lines(i1)*PI/180.0)

              CALL POSITN(hal_r(in),hal_z(in))
              CALL LINE  (deltar*FARAWAY,deltaz*FARAWAY)

              CALL POSITN(hal_r(in)+deltar*dfact(in),
     .                    hal_z(in)+deltaz*dfact(in))
              CALL TYPENF(ang_lines(i1),1)
            ENDDO

          ENDIF
        ENDIF

      ENDDO

      CALL FULL
c
      call lincol(1)
      call ctrmag(10)
c
      RETURN
99    STOP
      END
c
c ======================================================================
c








      SUBROUTINE CustomisePlot(title,xlab,ylab,elabs)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'colours'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      CHARACTER TITLE*80,XLAB*(*),YLAB*(*)
      CHARACTER*(*) elabs(MAXNGS)


      INTEGER       i1,in,in1,idum1
      REAL r,g,b,red,green,blue,frac
      CHARACTER*128 dataline,cdum1,cdum2

10    READ(5,'(A128)') dataline
      IF (dataline(8:11).EQ.'info'.OR.dataline(8:11).EQ.'Info'.OR.
     .    dataline(8:11).EQ.'INFO') THEN
        READ(dataline,*) cdum1,cdum2
        IF     (cdum2.EQ.'title') THEN
          READ(dataline,*) cdum1,cdum2,title
        ELSEIF (cdum2.EQ.'xaxis') THEN
          READ(dataline,*) cdum1,cdum2,xlab
        ELSEIF (cdum2.EQ.'yaxis') THEN
          READ(dataline,*) cdum1,cdum2,ylab
        ELSEIF (cdum2.EQ.'elabs') THEN
          READ(dataline,*) cdum1,cdum2,idum1,elabs(idum1)
        ELSEIF (cdum2.EQ.'type') THEN
          READ(dataline,*) cdum1,cdum2,slopt,slopt2,iopt_ghost
        ELSEIF (cdum2.EQ.'note') THEN
          READ(dataline,*) cdum1,cdum2,i1,char(i1)
        ELSEIF (cdum2.EQ.'colour') THEN
          READ(dataline,*) cdum1,cdum2,r,g,b
          CALL RGB
          ncols = ncols + 1
          DO in = 2, ncols
            in1 = ncols - in + 1
            frac = (REAL(in1) / REAL(ncols-1)) ** 1.5
            red   = 1.0 + frac * (r - 1.0)
            green = 1.0 + frac * (g - 1.0)
            blue  = 1.0 + frac * (b - 1.0)
c            red   = r * REAL(in-1) / REAL(ncols-1)
c            green = g * REAL(in-1) / REAL(ncols-1)
c            blue  = b * REAL(in-1) / REAL(ncols-1)
            CALL ColSet(red,green,blue,in)
            colour(in) = in
          ENDDO
          ncols = ncols - 1
        ELSE
          CALL ER('INFO','Unrecognized label : '//
     .                   cdum2(1:LEN_TRIM(cdum2)),*99)
        ENDIF
        GOTO 10
      ELSE
        BACKSPACE 5
      ENDIF

      RETURN
99    STOP
      END



      SUBROUTINE Plot974(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .                  xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .                  xouts,youts,nconts,conts,cntropt,n_cols,col_opt)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      include 'dynam2'
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

      INTEGER ndata
      REAL vdata(1000)
c
c     Array for extra comments 
c
      character*20 extra_comments(1)
c
c
c...Colours/shades:
      IF (n_cols.NE.-1.AND..FALSE.) THEN
        ncols = n_cols
        CALL HSI
        CALL ColSet(0.0,0.0,0.0,1)
        colour(1) = 1
        ncols = ncols + 1
        DO in = 2, ncols
          hue = 0.0
          saturation = 0.0
          value = REAL(in-1) / REAL(ncols-1)
          CALL ColSet(hue,saturation,value,in)
          colour(in) = in
        ENDDO

        DO i1 = in, in + 6
          colour(i1) = i1
        ENDDO
        CALL RGB
        CALL ColSet(1.0,0.0,0.0,ncols+1)
        CALL ColSet(0.0,1.0,0.0,ncols+2)
        CALL ColSet(0.0,0.0,1.0,ncols+3)
        CALL ColSet(0.5,0.7,0.2,ncols+4)
        CALL ColSet(0.5,0.0,0.5,ncols+5)
        CALL ColSet(0.0,0.5,0.5,ncols+6)
        ncols = ncols - 1
      ENDIF


      nview  = ' '
      plane  = ' '
      smooth = ' '
      anly   = ' '
      TABLE  = 'SYMBOL TABLE'

      XLAB   = '   X/A  (M)'
      YLAB   = '   Y/A  (M)'
      NGS    = CNGS
      IZ     = IOPT

c      WRITE(0,*) '974: ',iopt

      CALL RZero(cdata,MAXNKS*MAXNRS)

      ref(1:31) = graph(14:41)

      WRITE (IPLOT,9012) NPLOTS,REF


c      CALL CustomisePlot(title,xlab,ylab,elabs)


      IF     (iopt.EQ.1) THEN
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinion(ik,ir)
          ENDDO
         ENDDO
      ELSEIF (iopt.EQ.2) THEN
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinrec(ik,ir)
          ENDDO
         ENDDO
      ELSEIF (iopt.EQ.3) THEN
        DO ir = irsep,irwall-1
          DO ik = 1, nks(ir)
            cdata(ik,ir) = knbs(ik,ir)
          ENDDO
         ENDDO
      ELSEIF (iopt.EQ.4) THEN
        DO ir = 1, nrs
          IF (idring(ir).EQ.-1) CYCLE

          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinatom(ik,ir)
          ENDDO
        ENDDO
      ELSEIF (iopt.EQ.5) THEN
c
c       Ratio of PIN Hgamma to Halpha (excluding values of Halpha
c       less than 10% of peak):
c
        rdum1 = LO

        DO ir = 2, nrs
          DO ik = 1, nks(ir)
            rdum1 = MAX(rdum1,pinline(ik,ir,6,H_BALPHA))
          ENDDO
        ENDDO

        DO ir = 2, nrs
          DO ik = 1, nks(ir)
            IF (pinline(ik,ir,6,H_BALPHA).GT.0.10*rdum1) THEN
              cdata(ik,ir) = pinline(ik,ir,6,H_BGAMMA) /
     .                       pinline(ik,ir,6,H_BALPHA)
            ELSE
              cdata(ik,ir) = 0.0
            ENDIF
          ENDDO
         ENDDO
      ELSEIF (iopt.EQ.6) THEN
        DO ir = irsep, irwall-1
          DO ik = 1, nks(ir)
            cdata(ik,ir) = ktebs(ik,ir)
          ENDDO
         ENDDO

      ELSEIF (iopt.GE.7.AND.iopt.LE.11) THEN
c...    Contributions to Balmer alpha emission:
c          7 - H ionisation
c          8 - H+ recombination
c          9 - H2 dissociation
c         10 - H2+ dissociation
c         11 - CX of H and H+
c
        DO ir = 2, nrs
          IF (idring(ir).EQ.-1) CYCLE
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinline(ik,ir,iopt-6,H_BALPHA)
          ENDDO
        ENDDO

      ELSEIF (iopt.GE.12.AND.iopt.LE.20) THEN
c...    Stratum data:
c         12 - Ionistaion : stratum 1
c         13 -                      2
c         14 -                      3
c         15 - H density  : stratum 1
c         16 -                      2
c         17 -                      3
c         18 - H2 denisty : stratum 1
c         19 -                      2
c         20 -                      3
c
        STOP 'STOP 974: PINDATA ARRAY NOT AVAILABLE'

        DO ir = 1, nrs
          IF (idring(ir).EQ.-1) CYCLE

          DO ik = 1, nks(ir)
            cdata(ik,ir) = pindata(ik,ir,iopt-11)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.21) THEN
        DO ir = 1, nrs
          IF (idring(ir).EQ.-1) CYCLE

          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinmol(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.22) THEN

        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinrec(ik,ir) * kvols(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.23) THEN
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = kvhs(ik,ir) / qt /
     .                     GetCs(ktebs(ik,ir),ktibs(ik,ir))
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.24) THEN
c...    Dalpha:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinalpha(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.25) THEN
c...    D2+ density:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinmoi(ik,ir)
          ENDDO
        ENDDO
      ELSEIF (iopt.EQ.26) THEN
c...    Dalpha - core emission only:
        DO ir = 2, irsep-1
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinalpha(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.27) THEN
c...    Dgamma:

        WRITE(0,*) '974:27 Scaling Dgamma to kW m-3'
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinline(ik,ir,6,H_BGAMMA)*1.0E-3*
     .                    (6.63E-34*3.0E+08)/(4340.0*1.0E-10)
c            cdata(ik,ir) = pinline(ik,ir,6,H_BGAMMA) 
          ENDDO
        ENDDO
 
        char(30) = 'UNITS= kW m-3'

      ELSEIF (iopt.EQ.28) THEN
c...    Dgamma:
        DO ir = 1, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinline(ik,ir,6,H_BGAMMA) 
          ENDDO
        ENDDO

      ELSEIF (iopt.EQ.29) THEN
c...    NIU

      ELSEIF (iopt.GE.30.AND.iopt.LE.36) THEN
c...    Impurity radiation:
        DO ir = 2, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = sdlims(ik,ir,iopt-30)
          ENDDO
        ENDDO

      ELSEIF (iopt.GE.37) THEN
c...    Impurity neutral density:
        DO ir = 2, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinz0(ik,ir)
          ENDDO
        ENDDO

      ELSEIF (iopt.GE.38) THEN
c...    Impurity ion density:
        DO ir = 2, nrs
          DO ik = 1, nks(ir)
            cdata(ik,ir) = pinionz(ik,ir)
          ENDDO
        ENDDO

      ENDIF

      CALL GRTSET(TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     .  YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)


      count = 0
40    READ(5,'(A80)',END=50) graph1
      BACKSPACE 5
      IF (graph1(8:11).EQ.'View'.OR.graph1(8:11).EQ.'VIEW'.OR.
     .    graph1(8:11).EQ.'view') THEN

        count = count + 1
        CALL LinCol(ncols+count)
        CALL FilCol(ncols+count)

        CALL DrawHalpha(1,1)
        GOTO 40
      ELSE
        CALL LinCol(1)
        CALL FilCol(1)
      ENDIF
50    CONTINUE

c      IF (graph(41:44).NE.'    '.AND.graph(41:44).NE.' 0 0') THEN
c
c        READ(graph(41:42),'(I2)') i1
c        READ(graph(43:44),'(I2)') i2
c
c        DO ii = i1, i2
c          CALL DrawHalpha(ii,1)
c        ENDDO
c      ENDIF

      CALL SetPlotComments(974,job,extra_comments,0,0.0)
c
c jdemod - added 0.0, 0.0 for min/max scaling values for contours
c
      CALL CONTOUR(0,NGS,cdata,1,1,1,FT,FP,1.0,
     .             XOUTS,1,NXS,YOUTS,1,NYS,
     .             XXMIN,XXMAX,YYMIN,YYMAX,
     .             nconts,conts,cntropt,-0.01,-0.01)

c      CALL CONTOUR(ICNTR,NGS,cdata,1,1,1,FT,FP,1.0,
c     .             XOUTS,1,NXS,YOUTS,1,NYS,
c     .             XXMIN,XXMAX,YYMIN,YYMAX,
c     .             nconts,conts,cntropt,-0.01,-0.01)

c     .             nconts,conts,cntropt,0.0,0.0)
c
c jdemod
c

c...  Restore colours:
      IF (n_cols.NE.-1.AND..FALSE.) call setup_col(n_cols,col_opt)

      RETURN
99    STOP
9012  FORMAT(1X,'PLOT',I3,4X,A)
      END
