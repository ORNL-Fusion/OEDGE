      subroutine out200(iref,graph,iopt,ierr)
      use mod_params
      use mod_outcom
      use mod_cgeom
      use mod_comtor
      use mod_dynam2
      use mod_dynam3
      use mod_dynam4
      use mod_pindata
      implicit none
      integer iref,iopt,ierr
      character*(*) graph

c
c     include 'params'
c     include 'outcom'
c
c     Other common blocks
c
c     include 'cgeom'
c     include 'comtor'
c      include 'cneut2'
c     include 'dynam2'
c     include 'dynam3'
c     include 'dynam4'
c     include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c      include 'cedge2d'
c      include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Local Variables
c

      integer ik,ir,iz,it
      integer jk,jr,jz
      integer id,in,ii,jd



      integer startin,endin,stepin


c
c     Variables for time dependent LOS plots
c
      logical time_plot
      integer time_cnt


      integer iz1
      REAL WVALS(MAXTHE,MAXNGS)

      integer ikk

      REAL PLRPAD(MAXNKS,MAXNRS)
      logical plotr

      REAL GFF


      character*80 augform



      REAL TVALS2(MAXTHE,MAXNGS)



c
c     Exit if no plot requested 
c
      IF (IOPT.EQ.0) return


c
c     All plots in this module should be less than 310 - in the 
c     range 101 to 200 to be specific.  
c
c
c      IF (IREF.LT.310) THEN
c
        ref = graph(5:41)
        CALL RDG2 (GRAPH2,ROBS,ZOBS,THEMIN,DTHE,THERES,NUMTHE,
     >              IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
c
c       Make adjustment for ASDEX UPGRADE plots ... their angular
c       coordinate system is 180 degrees out of phase with JET.
c       i.e. Angles measured from -ve R-axis ... so add 180 to THEMIN
c       for ASDEX plots.
c
c       Not all SONNET grids are Asdex Upgrade grids - remove this for
c       now - replace with input parameter that will allow for
c       an arbitrary offset. (implement when needed)
c
c        if (cgridopt.eq.3) themin = themin + 180.0
c
        IF (IERR.EQ.1) THEN
           WRITE(6,*) 'RDG2 ERROR READING 200+SERIES- GRAPH DETAILS'
           IERR = 0
           return
        ENDIF
c
c      endif  
c

c


      call init_plot(iref,graph,iopt) 




c
C-----------------------------------------------------------------------
C
C     LINE OF SIGHT PLOTS FOR DENSITY AND PLRPS... FROM A GIVEN POSITION
C
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.201) THEN
C
C     LINE-INTEGRATED IMPURITY DENSITIES
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SCALED DENSITY'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
c
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 200 PLOT'
         ENDIF
c
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN+2,IZMAX+2,SDLIMS,MAXIZS+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2200 IZ = IZMIN,IZMAX
           DO 2200 IT = 1,NUMTHE
             IF (ABSFAC.GT.0.0)
     >         TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) *ABSFAC
2200     CONTINUE
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
c
c
      ELSEIF (IREF.EQ.206) THEN
C
C     LINE-INTEGRATED IMPURITY DENSITIES
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'DENSITY (M-2)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
         ENDIF
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
C
C  LOOP OVER IONISATION STATES
C
         DO IZ = IZMIN, IZMAX
           IZ1 = IZ - IZMIN + 1
C
C  LOAD INTEGRAND
C
           CALL RVALKR (CVALSA,SDLIMS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         ENDDO
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           DO IZ = IZMIN,IZMAX
             WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,IZ-IZMIN+1)
           ENDDO
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
c
c
      ELSEIF (IREF.EQ.208) THEN
C
C     LINE-INTEGRATED BACKGROUND PLASMA DENSITIES
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
c
         ELABS(1) = 'BG-NBG-NE LOS'
c
c
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'DENSITY (M-2)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
c
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
         ENDIF
c
c         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,knbs,1,1,1,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
c
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
c
c
      ELSEIF (IREF.EQ.211) THEN
C
C     CALCULATION BASED ON PLRPS
C
         YLAB = 'SCALED PLRPS'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'PLRP LOS PLOT'
c
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 210 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         IF (IZMIN.NE.-2) THEN
            IZMIN = PIND(IZMIN)
         ENDIF
         IZMAX = PIND(IZMAX+1)-1

         PSWITCH = .TRUE.
         PSHIFT = PNCNT
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,
     >               MAXPLRP+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
         PSWITCH = .FALSE.
         PSHIFT = 1
C
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2210 IZ = IZMIN,IZMAX
           DO 2210 IT = 1,NUMTHE
             IF (ABSFAC.GT.0.0)
     >         TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) *ABSFAC
2210     CONTINUE
c
c        Print for Asdex Upgrade only
c
         if (cgridopt.eq.3) then
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(50,'(1x,a)')
     >   'Data from Plot Option #211 (BLS - Wide Angle Overview)'
         write(50,'(1x,a,i4.4)') 'Number of points: ',numthe
         write(50,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a11,00(3x,a6,4x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
*        write(50,'(1x,a)') augform
         write(50,augform) 'Angle (Deg)',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f7.3,4x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
*        write(50,'(1x,a)') augform
         do 2211 it=1,numthe
           write(50,augform) touts(it),(tvals(it,iz),iz=1,izmax-izmin+1)
 2211    continue
c
         endif
c
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0,0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
c
      ELSEIF (IREF.EQ.212) THEN
C
C     LINE-INTEGRATED PHOTON EMISSION (PRLPS) FROM ADAS
c
c     FOR TIME DEPENDENT DATA.
C
C     VERIFY PARAMETERS
C
c
c     Time_plot is used to determine if one time dependent plot is
c     in use or a series of 2D LOS plots that evolve over time.
c
         time_plot = .false.
         time_cnt= 0
c
      write (6,*) 'PLOT 212:',nts
c
c
c     Quick and dirty method of making R or theta plots optional
c
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
c
C        ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
c
         IF (IZMIN.LT.0) IZMIN = 0
c
C        IZMAX IGNORED FOR ADAS PLRP INTEGRALS
c
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        Record starting THETA value
c
         themin_start = themin
C
C        LABELS AND SCALE FACTORS
C
         YLAB = 'PLRP (PH M-2 S-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'PLRP (PH M-2 S-1 SR-1)'
         ENDIF
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
c        LOOP OVER ALL OF THE STAGES OF THE TIME DEPENDENT RESULTS
c
         do it = 1,nts
C
C           CREATE THETA VECTOR
C
            themin = themin_start
c
            DO II = 1, NUMTHE
               TOUTS(II) = THEMIN + (II-1)*DTHE
               TWIDS(II) = THERES
            ENDDO
            THEMAX = TOUTS(NUMTHE)
c
c           Zero radiation array
c
            call rzero(plrpad,maxnks*maxnrs)
c
            WRITE (NAME,'(1P,E8.1,A)') DWELTS(IZMIN)*DWELFS(IT),'S'
c
c           Load radiation array
c
            call LDADAS_TIMEDEP(CION,IZMIN,IT,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
            IF (IRCODE.NE.0) THEN
               WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
               RETURN
            ENDIF
C
C  LOAD INTEGRAND
C
            CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >                   XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
            CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >                  ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
            PLABAD = '    XX XXXXX ('
            WRITE(PLABAD(5:6),'(I2)') IZMIN
            WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
            LEN = LENSTR(PLABAD)
            IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
            LEN = LENSTR(PLABAD)
            IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
            LEN = LENSTR(PLABAD)
            IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
            LEN = LENSTR(PLABAD)
            PLABAD = PLABAD(1:LEN) // ') '
C
            WRITE (iplot,9012) NPLOTS, REF
            WRITE (iplot,9040) PLABAD
            WRITE (iplot,*) NVIEW
            WRITE (iplot,*) PLANE
            WRITE (iplot,*) ANLY
c
c           Adjust TOUTS to against R instead of Theta if required.
c
            if (plotr) then
               call adjustout(touts,numthe,zadj,robs,zobs)
               themin = touts(1)
               themax = touts(numthe)
               write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
            endif
c
            write (6,*) '212:NUMTHE:',numthe,name,themin,themax
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
            IF (NUMTHE.GT.1) THEN
c
               ELABS(1) = PLABAD
c
               if (iseld.ne.0) then

                  ngs = 2

                  call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                      themin,themax,maxngs,ngs,datatitle)

                  ELABS(2) = 'EXPT'//DATATITLE
c
                  do in = 1,numthe
                     write(6,*) 'TVALS:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
                  end do
c
                  CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,NGS,
     >                 ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >                 JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >                 TABLE,IOPT,2,1.0,0)
c
               else
c
                  ngs = 1
c
                  CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >                 ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >                 JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >                 TABLE,IOPT,2,1.0,0)
c
               endif
c
            ELSE
c
c              Single LOS
c
               time_plot=.true.
               time_cnt = time_cnt+1
c
c              Record the time dependent data
c
               louts(time_cnt) = (ctimes(time_cnt,izmin)
     >                            +ctimes(time_cnt-1,izmin))/2.0
c
               lvals(time_cnt,1) = tvals(1,1)
c
               WRITE(iplot,'(a,i4,a,f8.3,3(a,g12.5))')
     >                   'RESULT 212: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   ' TIME = ',louts(time_cnt),
     >                   ' VAL = ',TVALS(1,1)
c
            endif
c
         end do
c
c        Create a time evolving plot of single point LOS's
c
         if (time_plot) then
c
            themin = 0.0
            themax = louts(time_cnt)
            XLAB = 'TIME (SECONDS)'
            ELABS(1) = PLABAD
c
            if (iseld.gt.0) then
c
               ngs = 2
c
               call calc_expt(iseld,louts,lvals,maxseg,time_cnt,
     >                      themin,themax,maxngs,ngs,datatitle)

               ELABS(2) = 'EXPT'//DATATITLE
c
               do in = 1,time_cnt
                  write(6,*) 'LVALS:',in,louts(in),lvals(in,1),
     >                                  lvals(in,2)
               end do
c
               CALL DRAW(LOUTS,LWIDS,LVALS,MAXSEG,time_cnt,ANLY,NGS,
     >                 ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >                 JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >                 TABLE,IOPT,2,1.0,0)
c
            else
c
               ngs = 1
c
               do in = 1,time_cnt
                  write(6,*) 'LVALS:',in,louts(in),lvals(in,1)
               end do
c
               CALL DRAW(LOUTS,LWIDS,LVALS,MAXSEG,time_cnt,ANLY,NGS,
     >                 ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >                 JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >                 TABLE,IOPT,2,1.0,0)
c
            endif
c
         endif
c
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
c
c
      ELSEIF (IREF.EQ.216) THEN
C
C     LINE-INTEGRATED PHOTON EMISSION (PRLPS) FROM ADAS
C
C     VERIFY PARAMETERS
C
c
c        Quick and dirty method of making R or theta plots optional
c

         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'PLRP (PH M-2 S-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'PLRP (PH M-2 S-1 SR-1)'
         ENDIF
c
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
c
          call LDADAS(CION,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         PLABAD = '    XX XXXXX ('
         WRITE(PLABAD(5:6),'(I2)') IZMIN
         WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
         LEN = LENSTR(PLABAD)
         IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
         LEN = LENSTR(PLABAD)
         IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
         LEN = LENSTR(PLABAD)
         IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
         LEN = LENSTR(PLABAD)
         PLABAD = PLABAD(1:LEN) // ') '
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            ELABS(1) = PLABAD
c
            if (iseld.ne.0) then

               ngs = 2

               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(2) = 'EXPT'//DATATITLE
c
               do in = 1,numthe
                  write(6,*) 'TVALS:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
               end do
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,NGS,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif

CPRINT
c           DO II = 1, NUMTHE
c             WRITE(6,*) TOUTS(II),TVALS(II,1)
c           ENDDO
CPRINT
         ELSE
c
c           Get experimental data
c
            if (iseld.gt.0) then
c
               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 216: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1),' EXPT = ',tvals(1,2)

c
            else
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 216: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1)

c
            endif
c
c           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,1)
c
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.217) THEN
C
C     LINE-INTEGRATED PHOTON EMISSION (PLRPS) RATIOS FROM ADAS
C
C     VERIFY PARAMETERS
C
c
c        Quick and dirty method of making R or theta plots optional
c

         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
c
c        Both IZMIN and IZMAX are ignored for ratio plots
c

C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS


         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS

         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'RATIO'
         XLAB = 'THETA (DEGREES)'
c
c         REF  = GRAPH(5:41)
c
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
c
c
c        ANLY: Scale factor not relevant for ratios
c
c
c        Fill in ANLY
c
c         IF (ATYPE.EQ.0) THEN
c           MFACT = 1.0
c           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
c         ELSEIF (ATYPE.EQ.1) THEN
c           MFACT = DTHE / (2.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
c     >                DTHE
c         ELSEIF (ATYPE.EQ.2) THEN
c           MFACT = 1.0 / (2.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
c         ELSEIF (ATYPE.EQ.3) THEN
c           MFACT = 1.0 / ( 4.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
c           YLAB = 'PLRP (PH M-2 S-1 SR-1)'
c         ENDIF
c
c
c
c       Load data for the first line
c
c       Verify that it is an acceptable line
c
        if ((z_atom.eq.cion.and.iz_state.ge.0.and.iz_state.le.nizs).or.
     >      (z_atom.eq.cizb.and.iz_state.ge.0.and.iz_state.le.cizb))
     >       then
c
          call LDADAS(z_atom,IZ_state,ADASID,ADASYR,ADASEX,
     >                  ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
          REF = 'ADAS Z=XX XX XXXXX ('
          WRITE(ref(8:9),'(I2)') Z_atom
          WRITE(REF(11:12),'(I2)') IZ_state
          WRITE(REF(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(REF)
          IF (ISELE.GT.0) REF = REF(1:LEN) // 'E'
          LEN = LENSTR(REF)
          IF (ISELR.GT.0) REF = REF(1:LEN) // 'R'
          LEN = LENSTR(REF)
          IF (ISELX.GT.0) REF = REF(1:LEN) // 'C'
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // ') '
          LEN = LENSTR(REF)
          REF = REF(1:LEN) // XPOINT
          ref = 'NUMER:'//ref
          WRITE (IPLOT,9012) NPLOTS,REF
c
c         calculate LOS plot
c
C
C  LOAD INTEGRAND
C
          MFACT = 1.0
c
c         Set absolute scaling factor only for impurity species
c
          IF (ABSFAC.GT.0.0.and.z_atom.eq.cion)
     >               MFACT = MFACT * ABSFAC
c
          CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >                 XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
c
C  INTEGRATE
C
          CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >                ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
        else
          WRITE(6,*) 'PLOT 217 - INVALID ION AND CHARGE STATE'//
     >               ' SPECIFIED WITH FIRST ADAS INPUT LINE'
          RETURN
        endif
c
c
c       Load data for the second line
c
c       Verify that it is an acceptable line
c
        if ((z_atom2.eq.cion.and.iz_state2.ge.0.and.iz_state2.le.nizs)
     >               .or.
     >      (z_atom2.eq.cizb.and.iz_state2.ge.0.and.iz_state2.le.cizb))
     >             then
c
          call LDADAS(z_atom2,IZ_state2,ADASID2,ADASYR2,ADASEX2,
     >                  ISELE2,ISELR2,ISELX2,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
          anly = 'ADAS Z=XX XX XXXXX ('
          WRITE(anly(8:9),'(I2)') Z_atom2
          WRITE(anly(11:12),'(I2)') IZ_state2
          WRITE(anly(14:18),'(I5)') NINT(WLNGTH)
          LEN = LENSTR(anly)
          IF (ISELE2.GT.0) anly = anly(1:LEN) // 'E'
          LEN = LENSTR(anly)
          IF (ISELR2.GT.0) anly = anly(1:LEN) // 'R'
          LEN = LENSTR(anly)
          IF (ISELX2.GT.0) anly = anly(1:LEN) // 'C'
          LEN = LENSTR(anly)
          anly = anly(1:len) // ') '
          LEN = LENSTR(anly)
          anly = anly(1:LEN) // XPOINT
          anly = 'DENOM:'//anly
          WRITE (IPLOT,9012) NPLOTS,REF
c
C
C  LOAD INTEGRAND
C
          MFACT = 1.0
c
c         Set absolute scaling factor only for impurity species
c
          IF (ABSFAC.GT.0.0.and.z_atom.eq.cion)
     >               MFACT = MFACT * ABSFAC
c
          CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >                 XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
c
C  INTEGRATE
C
          CALL LOSINT(TVALS2,TOUTS,TWIDS,NUMTHE,
     >                ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
c
        else
          WRITE(6,*) 'PLOT 217 - INVALID ION AND CHARGE STATE'//
     >               ' SPECIFIED WITH SECOND ADAS INPUT LINE'
          RETURN
        endif
c
c
        PLABAD = 'DIV DIVIMP RATIO'
c
c       Calculate RATIO and call the drawing routine to plot it -
c       if the denominator contribution is 0.0 then the ratio is not calculated
c       for that cell.
c
        do in = 1,numthe

           if (tvals2(in,1).gt.0.0) then
              tvals(in,1) = tvals(in,1) / tvals2(in,1)
           else
              tvals(in,1) = 0.0
           endif

        end do
C
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            ELABS(1) = PLABAD
c
            if (iseld.ne.0.and.iseld2.ne.0) then

               ngs = 2

               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               call calc_expt(iseld2,touts,tvals2,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
c
               do in = 1,numthe
c
                  if (tvals2(in,2).gt.0.0) then
                    tvals(in,2) = tvals(in,2) / tvals2(in,2)
                  else
                    tvals(in,2) = 0.0
                  endif
c
               end do
c
               ELABS(2) = 'EXPTEXPT RATIO'
c
               do in = 1,numthe
                  write(6,*) 'TVALS:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
               end do
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,NGS,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif

CPRINT
c           DO II = 1, NUMTHE
c             WRITE(6,*) TOUTS(II),TVALS(II,1)
c           ENDDO
CPRINT
         ELSE
c
c           Get experimental data
c
            if (iseld.gt.0.and.iseld2.gt.0) then
c
               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
c
               call calc_expt(iseld2,touts,tvals2,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
c
c              Find ratio
c
               if (tvals2(1,2).gt.0.0) then
                  tvals(1,2) = tvals(1,2)/tvals2(1,2)
               else
                  tvals(1,2) = 0.0
               endif
c
c              Write out results
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 217: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1),' EXPT = ',tvals(1,2)

c
            else
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 217: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1)

c
            endif
c
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.218) THEN
C
C     LINE-INTEGRATED HYDROGEN PHOTON EMISSION (PRLPS) FROM ADAS
C
C     VERIFY PARAMETERS
C
c
c        Quick and dirty method of making R or theta plots optional
c
c         write (6,*) 'Doing plot 218:'
c
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'PLRP (PH M-2 S-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
c
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'PLRP (PH M-2 S-1 SR-1)'
         ENDIF
c
          call LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)

C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
c
         PLABAD = '    XX XXXXX ('
         WRITE(PLABAD(5:6),'(I2)') IZMIN
         WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
         LEN = LENSTR(PLABAD)
         IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
         LEN = LENSTR(PLABAD)
         IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
         LEN = LENSTR(PLABAD)
         IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
         LEN = LENSTR(PLABAD)
         PLABAD = PLABAD(1:LEN) // ') '
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
c
c        Adjust Theta to R plotting coordinates.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            ELABS(1) = PLABAD
c
            if (iseld.ne.0) then

               ngs = 2

               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(2) = 'EXPT'//DATATITLE
c
               do in = 1,numthe
                  write(6,*) 'TVALS:',in,touts(in),tvals(in,1),
     >                                  tvals(in,2)
               end do
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)

            else
c
                CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif
c
         ELSE
c
c           Get experimental data
c
            if (iseld.gt.0) then
c
               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 218: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1),' EXPT = ',tvals(1,2)

               WRITE(6,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 218: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1),' EXPT = ',tvals(1,2)

c
            else
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 218: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1)

               WRITE(6,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT 218: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1)

c
            endif
c
c           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,1)
            WRITE(6,*) 'Result 218: ', IZMIN,TOUTS(1),TVALS(1,1)
c
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
c
      ELSEIF (IREF.EQ.231) THEN
C
C     CALCULATION BASED ON ION POWER LOSS
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SCALED POWER LOSS'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'POWER LOSS LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 230 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN+2,IZMAX+2,POWLS,MAXIZS+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2220 IZ = IZMIN,IZMAX
           DO 2220 IT = 1,NUMTHE
             IF (ABSFAC.GT.0.0)
     >         TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) *ABSFAC
2220     CONTINUE
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
      ELSEIF (IREF.EQ.236) THEN
C
C     LINE-INTEGRATED RADIATED POWER
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'PRAD (W M-2)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
         ENDIF
         IF (ABSFAC.GT.0.0) MFACT = MFACT * ABSFAC
C
C  LOOP OVER IONISATION STATES
C
         DO IZ = IZMIN, IZMAX
           IZ1 = IZ - IZMIN + 1
C
C  LOAD INTEGRAND
C
           CALL RVALKR (CVALSA,POWLS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         ENDDO
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
CPRINT
         DO IZ = IZMIN,IZMAX
           DO II = 1,NUMTHE
             WRITE(iplot,*) IZ,TOUTS(II),TVALS(II,IZ-IZMIN+1)
           ENDDO
         ENDDO
CPRINT
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           DO IZ = IZMIN,IZMAX
             WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,IZ-IZMIN+1)
           ENDDO
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.238) THEN
C
C     LINE-INTEGRATED HYDROGEN RADIATED POWER
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.0) IZMIN = 0
         IF (IZMAX.GT.2) IZMAX = 2
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
         write (6,*) '238:', numthe,avpts

C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'HYDROGEN PRAD (W M-2)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
         ENDIF
C
C  LOOP OVER IONISATION STATES
C
         DO IZ = IZMIN, IZMAX
           IZ1 = IZ - IZMIN + 1
C
C  LOAD INTEGRAND
C
           CALL HVALKR (CVALSA,HPOWLS,IZ+1,2,2,FT,FP,1.0,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         ENDDO
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY
         write (iplot,*) 'IZMIN=',izmin,' IZMAX=',izmax,
     >                   ' NUMTHE=',numthe,' AVPTS=',avpts
CPRINT
         DO IZ = IZMIN,IZMAX
           DO II = 1,NUMTHE
             WRITE(iplot,*) IZ,TOUTS(II),TVALS(II,IZ-IZMIN+1)
           ENDDO
         ENDDO
CPRINT
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,HLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           DO IZ = IZMIN,IZMAX
             WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,IZ-IZMIN+1)
           ENDDO
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.241) THEN
C
C     CALCULATION BASED ON SPECTROSCOPICALLY WEIGHTED TEMPERATURES
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SPECT. TEMP.'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'SPECTOSCOPIC T LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 240 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         IF (IZMIN.NE.-2) THEN
            IZMIN = PIND(IZMIN)
         ENDIF
         IZMAX = PIND(IZMAX+1)-1

         PSWITCH = .TRUE.
         PSHIFT = PNCNT
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,MAXPLRP+2,
     >               SDTS,KTEBS,1,ANLY,PSWITCH,PIZS,FT,FP)
         PSWITCH = .FALSE.
         PSHIFT = 1
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.246) THEN
C
C     EMISSION-WEIGHTED, LINE-INTEGRATED IMPURITY TEMPERATURE
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'EMIS-WGHT IMP TEMP (EV)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         ANLY = ' '
         MFACT = 1.0
c
          call LDADAS(CION,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
c
c --- delete ---
C
C  LOAD ADAS-BASED PLRP
C
c          CALL XXUID(ADASID)
c          DO IR = 1,NRS
c            DO IK = 1,NKS(IR),20
c              NPAIRS = MIN0(20,NKS(IR)-(IK-1))
c              DO IADAS = 1,NPAIRS
c                TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
c                DADAS(IADAS) = DBLE(1.E-6*KNBS(IK+(IADAS-1),IR)*RIZB)
c              ENDDO
c              CALL DZERO(PECA1,NPAIRS)
c              IF (ISELE.GT.0) THEN
c                CALL SPEC(ISELE,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA1,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELE.EQ.-1) THEN
c                CALL DINIT(PECA1,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              CALL DZERO(PECA2,NPAIRS)
c              IF (ISELR.GT.0) THEN
c                CALL SPEC(ISELR,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA2,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELR.EQ.-1) THEN
c                CALL DINIT(PECA2,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
cldh - use ion temperature for CX rate
c              DO IADAS = 1,NPAIRS
c                TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
c              ENDDO
cldh
c              CALL DZERO(PECA3,NPAIRS)
c              IF (ISELX.GT.0) THEN
c                CALL SPEC(ISELX,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA3,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELX.EQ.-1) THEN
c                CALL DINIT(PECA3,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              DO IADAS = 1,NPAIRS
c                IKK = IK + (IADAS-1)
c                PLRPAD(IKK,IR) = 1.E-6*
c     >    ( KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN  ) *PECA1(IADAS)
c     >    + KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN+1) *PECA2(IADAS)
c     >    + PINATOM(IKK,IR)   *SDLIMS(IKK,IR,IZMIN+1) *PECA3(IADAS))
c              ENDDO
c            ENDDO
c          ENDDO


C
C  LOAD INTEGRAND FOR DENOMINATOR
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
C  LOAD INTEGRAND FOR NUMERATOR
C
         DO IR = 1,NRS
           DO IK = 1,NKS(IR)
             PLASTMP(IK,IR) = SDTS(IK,IR,IZMIN) * PLRPAD(IK,IR)
           ENDDO
         ENDDO
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(WVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         DO II = 1,NUMTHE
           IF (TVALS(II,1).GT.0.0) THEN
             TVALS(II,1) = WVALS(II,1) / TVALS(II,1)
           ELSE
             TVALS(II,1) = 0.0
           ENDIF
         ENDDO
C
         PLABAD = '    XX XXXXX ('
         WRITE(PLABAD(5:6),'(I2)') IZMIN
         WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
         LEN = LENSTR(PLABAD)
         IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
         LEN = LENSTR(PLABAD)
         IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
         LEN = LENSTR(PLABAD)
         IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
         LEN = LENSTR(PLABAD)
         PLABAD = PLABAD(1:LEN) // ') '
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,1)
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.251) THEN
C
C     CALCULATION BASED ON DENSITY WEIGHTED TEMPERATURES
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'DENSITY TEMP.'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'DENSITY T LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 250 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN+2,IZMAX+2,SDLIMS,MAXIZS+2,
     >               SDTS,KTEBS,1,ANLY,PSWITCH,PIZS,FT,FP)
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.256) THEN
C
C     DENSITY-WEIGHTED, LINE-INTEGRATED IMPURITY TEMPERATURE
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  NO SECONDARY NEUTRAL WEIGHTING
         IF (IZMIN.LT.-1) IZMIN = -1
C  NO SUMS OVER IONISATION STATES
         IF (IZMAX.GT.NIZS) IZMAX = NIZS
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'DENS-WGHT IMP TEMP (EV)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         MFACT = 1.0
         ANLY = ' '
C
C  LOOP OVER IONISATION STATES
C
         DO IZ = IZMIN, IZMAX
           IZ1 = IZ - IZMIN + 1
C
C  LOAD INTEGRAND FOR DENOMINATOR
C
           CALL RVALKR (CVALSA,SDLIMS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
C  LOAD INTEGRAND FOR NUMERATOR
C
           DO IR = 1,NRS
             DO IK = 1,NKS(IR)
               PLASTMP(IK,IR) = SDTS(IK,IR,IZ) * SDLIMS(IK,IR,IZ)
             ENDDO
           ENDDO
           CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(WVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
           DO II = 1,NUMTHE
             IF (TVALS(II,IZ1).GT.0.0) THEN
               TVALS(II,IZ1) = WVALS(II,IZ1) / TVALS(II,IZ1)
             ELSE
               TVALS(II,IZ1) = 0.0
             ENDIF
           ENDDO
C
         ENDDO
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           DO IZ = IZMIN,IZMAX
             WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,IZ-IZMIN+1)
           ENDDO
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.261) THEN
C
C     CALCULATION BASED ON SPECTROSCOPICALLY WEIGHTED
C     BACKGROUND (H or D) ION TEMPERATURE
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SPECT. ION TEMP.'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'SPECTOSCOPIC T LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 260 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         IF (IZMIN.NE.-2) THEN
            IZMIN = PIND(IZMIN)
         ENDIF
         IZMAX = PIND(IZMAX+1)-1

         PSWITCH = .TRUE.
         PSHIFT = PNCNT
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,MAXPLRP+2,
     >               SDTS,KTIBS,2,ANLY,PSWITCH,PIZS,FT,FP)
         PSWITCH = .FALSE.
         PSHIFT = 1
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.266) THEN
C
C     EMISSION-WEIGHTED, LINE-INTEGRATED ION TEMPERATURE
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'EMIS-WGHT ION TEMP (EV)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         ANLY = ' '
         MFACT = 1.0
c
          call LDADAS(CION,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
c
c --- delete ---
C
C  LOAD ADAS-BASED PLRP
C
c          CALL XXUID(ADASID)
c          DO IR = 1,NRS
c            DO IK = 1,NKS(IR),20
c              NPAIRS = MIN0(20,NKS(IR)-(IK-1))
c              DO IADAS = 1,NPAIRS
c                TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
c                DADAS(IADAS) = DBLE(1.E-6*KNBS(IK+(IADAS-1),IR)*RIZB)
c              ENDDO
c              CALL DZERO(PECA1,NPAIRS)
c              IF (ISELE.GT.0) THEN
c                CALL SPEC(ISELE,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA1,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELE.EQ.-1) THEN
c                CALL DINIT(PECA1,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              CALL DZERO(PECA2,NPAIRS)
c              IF (ISELR.GT.0) THEN
c                CALL SPEC(ISELR,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA2,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELR.EQ.-1) THEN
c                CALL DINIT(PECA2,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
cldh - use ion temperature for CX rate
c              DO IADAS = 1,NPAIRS
c               TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
c              ENDDO
cldh
c              CALL DZERO(PECA3,NPAIRS)
c              IF (ISELX.GT.0) THEN
c                CALL SPEC(ISELX,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA3,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELX.EQ.-1) THEN
c                CALL DINIT(PECA3,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              DO IADAS = 1,NPAIRS
c                IKK = IK + (IADAS-1)
c                PLRPAD(IKK,IR) = 1.E-6*
c     >    ( KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN  ) *PECA1(IADAS)
c     >    + KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN+1) *PECA2(IADAS)
c     >    + PINATOM(IKK,IR)   *SDLIMS(IKK,IR,IZMIN+1) *PECA3(IADAS))
c              ENDDO
c            ENDDO
c          ENDDO



C
C  LOAD INTEGRAND FOR DENOMINATOR
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
C  LOAD INTEGRAND FOR NUMERATOR
C
         DO IR = 1,NRS
           DO IK = 1,NKS(IR)
             PLASTMP(IK,IR) = KTIBS(IK,IR) * PLRPAD(IK,IR)
           ENDDO
         ENDDO
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(WVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         DO II = 1,NUMTHE
           IF (TVALS(II,1).GT.0.0) THEN
             TVALS(II,1) = WVALS(II,1) / TVALS(II,1)
           ELSE
             TVALS(II,1) = 0.0
           ENDIF
         ENDDO
C
         PLABAD = '    XX XXXXX ('
         WRITE(PLABAD(5:6),'(I2)') IZMIN
         WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
         LEN = LENSTR(PLABAD)
         IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
         LEN = LENSTR(PLABAD)
         IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
         LEN = LENSTR(PLABAD)
         IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
         LEN = LENSTR(PLABAD)
         PLABAD = PLABAD(1:LEN) // ') '
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,1)
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.271) THEN
C
C     CALCULATION BASED ON SPECTROSCOPICALLY WEIGHTED
C     BACKGROUND ELECTRON TEMPERATURE
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SPECT. ELECTRON TEMP.'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'SPECTOSCOPIC T LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 240 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         IF (IZMIN.NE.-2) THEN
            IZMIN = PIND(IZMIN)
         ENDIF
         IZMAX = PIND(IZMAX+1)-1

         PSWITCH = .TRUE.
         PSHIFT = PNCNT
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,MAXPLRP+2,
     >               SDTS,KTEBS,2,ANLY,PSWITCH,PIZS,FT,FP)
         PSWITCH = .FALSE.
         PSHIFT = 1
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.276) THEN
C
C     EMISSION-WEIGHTED, LINE-INTEGRATED ELECTRON TEMPERATURE
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
         IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'EMIS-WGHT ELEC TEMP (EV)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         ANLY = ' '
         MFACT = 1.0
c
          call LDADAS(CION,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >                  plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             RETURN
          ENDIF
c
c
c --- delete ---
C
C  LOAD ADAS-BASED PLRP
C
c          CALL XXUID(ADASID)
c          DO IR = 1,NRS
c            DO IK = 1,NKS(IR),20
c              NPAIRS = MIN0(20,NKS(IR)-(IK-1))
c              DO IADAS = 1,NPAIRS
c                TADAS(IADAS) = DBLE(KTEBS(IK+(IADAS-1),IR))
c                DADAS(IADAS) = DBLE(1.E-6*KNBS(IK+(IADAS-1),IR)*RIZB)
c              ENDDO
c              CALL DZERO(PECA1,NPAIRS)
c              IF (ISELE.GT.0) THEN
c                CALL SPEC(ISELE,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA1,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELE.EQ.-1) THEN
c                CALL DINIT(PECA1,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              CALL DZERO(PECA2,NPAIRS)
c              IF (ISELR.GT.0) THEN
c                CALL SPEC(ISELR,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA2,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELR.EQ.-1) THEN
c                CALL DINIT(PECA2,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
ccldh - use ion temperature for CX rate
c              DO IADAS = 1,NPAIRS
c                TADAS(IADAS) = DBLE(KTIBS(IK+(IADAS-1),IR))
c              ENDDO
ccldh
c              CALL DZERO(PECA3,NPAIRS)
c              IF (ISELX.GT.0) THEN
c                CALL SPEC(ISELX,IZMIN,CION,NPAIRS,TADAS,DADAS,
c     >                    WLNGTH,PECA3,LTRNG,LDRNG,PECTITLE,IRCODE)
c                IF (IRCODE.NE.0) THEN
c                  WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
c                  RETURN
c                ENDIF
c              ELSE IF (ISELX.EQ.-1) THEN
c                CALL DINIT(PECA3,NPAIRS,1.D6)
c                WLNGTH = 0.0
c              ENDIF
c              DO IADAS = 1,NPAIRS
c                IKK = IK + (IADAS-1)
c                PLRPAD(IKK,IR) = 1.E-6*
c     >    ( KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN  ) *PECA1(IADAS)
c     >    + KNBS(IKK,IR)*RIZB *SDLIMS(IKK,IR,IZMIN+1) *PECA2(IADAS)
c     >    + PINATOM(IKK,IR)   *SDLIMS(IKK,IR,IZMIN+1) *PECA3(IADAS))
c              ENDDO
c            ENDDO
c          ENDDO


C
C  LOAD INTEGRAND FOR DENOMINATOR
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
C  LOAD INTEGRAND FOR NUMERATOR
C
         DO IR = 1,NRS
           DO IK = 1,NKS(IR)
             PLASTMP(IK,IR) = KTEBS(IK,IR) * PLRPAD(IK,IR)
           ENDDO
         ENDDO
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(WVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         DO II = 1,NUMTHE
           IF (TVALS(II,1).GT.0.0) THEN
             TVALS(II,1) = WVALS(II,1) / TVALS(II,1)
           ELSE
             TVALS(II,1) = 0.0
           ENDIF
         ENDDO
C
         PLABAD = '    XX XXXXX ('
         WRITE(PLABAD(5:6),'(I2)') IZMIN
         WRITE(PLABAD(8:12),'(I5)') NINT(WLNGTH)
         LEN = LENSTR(PLABAD)
         IF (ISELE.GT.0) PLABAD = PLABAD(1:LEN) // 'E'
         LEN = LENSTR(PLABAD)
         IF (ISELR.GT.0) PLABAD = PLABAD(1:LEN) // 'R'
         LEN = LENSTR(PLABAD)
         IF (ISELX.GT.0) PLABAD = PLABAD(1:LEN) // 'C'
         LEN = LENSTR(PLABAD)
         PLABAD = PLABAD(1:LEN) // ') '
C
         WRITE (iplot,9012) NPLOTS, REF
         WRITE (iplot,9040) PLABAD
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,1)
         ENDIF
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.281) THEN
C
C     CALCULATION BASED ON DENSITY WEIGHTED K-VALUES
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'DENSITY WEIGHTED K-VALUE'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
C        REF  = 'DENSITY K LOS PLOT'
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
         IF (IZMIN.LT.-2) IZMIN = -2
         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 250 PLOT'
         ENDIF
C
C        LOAD 2-D ARRAY WITH K-VALUES AT EACH LOCATION
C
         DO 2600 IR = 1,NRS
           DO 2600 IK = 1,NKS(IR)
             KTMP(IK,IR) = KKS(IR)
 2600    CONTINUE
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN+2,IZMAX+2,SDLIMS,MAXIZS+2,
     >               SDTS,KTMP,2,ANLY,PSWITCH,PIZS,FT,FP)
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.286) THEN
C
C     DENSITY-WEIGHTED, LINE-INTEGRATED K-VALUES
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  NO SECONDARY NEUTRAL WEIGHTING
         IF (IZMIN.LT.-1) IZMIN = -1
C  NO SUMS OVER IONISATION STATES
         IF (IZMAX.GT.NIZS) IZMAX = NIZS
         IF (IZMIN.GT.IZMAX) IZMIN=IZMAX
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'DENS-WGHT K-VALUE'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         MFACT = 1.0
         ANLY = ' '
C
C  LOOP OVER IONISATION STATES
C
         DO IZ = IZMIN, IZMAX
           IZ1 = IZ - IZMIN + 1
C
C  LOAD INTEGRAND FOR DENOMINATOR
C
           CALL RVALKR (CVALSA,SDLIMS,IZ+2,NIZS+2,MAXIZS+2,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
C  LOAD INTEGRAND FOR NUMERATOR
C
           DO IR = 1,NRS
             DO IK = 1,NKS(IR)
               PLASTMP(IK,IR) = KKS(IR) * SDLIMS(IK,IR,IZ)
             ENDDO
           ENDDO
           CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >       XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
           CALL LOSINT(WVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >                 ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
           DO II = 1,NUMTHE
             IF (TVALS(II,IZ1).GT.0.0) THEN
               TVALS(II,IZ1) = WVALS(II,IZ1) / TVALS(II,IZ1)
             ELSE
               TVALS(II,IZ1) = 0.0
             ENDIF
           ENDDO
C
         ENDDO
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           DO IZ = IZMIN,IZMAX
             WRITE(iplot,*) IZ,TOUTS(1),TVALS(1,IZ-IZMIN+1)
           ENDDO
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.291) THEN
C
C     LOS PLOT OF PIN H-ALPHA DATA
C
C     CALCULATION BASED ON PIN H-ALPHA DATA
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SCALED H-ALPHA'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
            WRITE(iplot,*)'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C
C        OVERWRITE INPUT SINCE ONLY ONE SET OF DATA FOR H-ALPHA IS
C        AVAILABLE
C
         IZMIN = 1
         IZMAX = 1
C
C         IF (IZMIN.LT.-2) IZMIN = -2
C         IF (IZMAX.GT.NIZS+1) IZMAX = NIZS+1
C
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 290 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN,IZMAX,PINALPHA,1,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
C
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)

C
c        we want all data written in a file as ASCII to compare with
c        experimental values; added by Krieger, IPP 6/95

         write(6,'(1x,a,i4.4)') 'Number of points: ',numthe
         augform='(1x,a11,3x,a6)'
         write(6,augform)       'Angle (Deg)','Halpha'
         augform='(1x,f7.3,4x,1p,2x,e11.4)'
         do 2221 it=1,numthe
           write(6,augform) touts(it),tvals(it,1)
 2221    continue
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
      ELSEIF (IREF.EQ.296) THEN
C
C     LINE-INTEGRATED PIN H-ALPHA EMISSION
C
C     VERIFY PARAMETERS
C
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
c
c        Use numsmooth as a quick way to pass a different parameter.
c        Experimental data set selector.
c
         iseld = 0
         plotr = .false.
c
         if (numsmooth.lt.0) then
c
            iseld = abs(numsmooth)
            numsmooth = 1
c
         endif
c
C  NEUTRAL!
         IF (IZMIN.NE.0) IZMIN = 0
C  IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
         IZ1 = 1
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'PIN H-ALPHA (PH M-2 S-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'PIN H-ALPHA (PH M-2 S-1 SR-1)'
         ENDIF
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,PINALPHA,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,IZ1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (iplot,9012) NPLOTS,REF
         WRITE (iplot,*) NVIEW
         WRITE (iplot,*) PLANE
         WRITE (iplot,*) ANLY

         write (6,*) '296:',iseld,numsmooth,numthe ,plotr

c slmod begin - new
         CALL SetPlotComments(iref,job,extra_comments,0,0.0)
c slmod end
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN

            ELABS(1) = 'PIN PIN H-ALPHA'

            if (iseld.ne.0) then

               ngs = 2

               call calc_expt(iseld,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(2) = 'EXPT'//DATATITLE
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)

            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)

c
            endif

         ELSE
           WRITE(iplot,*) IZMIN,TOUTS(1),TVALS(1,IZ1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
c
      ELSEIF (IREF.EQ.297) THEN
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION
C
C     VERIFY PARAMETERS
C
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  IZMIN AND IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
         ENDIF
c
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
         WLNGTH = 5235.0
C
C  LOAD INTEGRAND
C
c         CALL LDBREM(plastmp,CVALS,IRCODE,NIZS)
c
c delete
c
c
         DO IR = 1,NRS
           DO IK = 1,NKS(IR)
             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
     >                       *(RIZB*KNBS(IK,IR))**2
     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
           ENDDO
         ENDDO
c
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
         ELSE
           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
C
      ELSEIF (IREF.EQ.298) THEN
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION
C
C     VERIFY PARAMETERS
C
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  IZMIN AND IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)// '- NEW METHOD'
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
         ENDIF
c
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
         WLNGTH = 5235.0
C
C  LOAD INTEGRAND
C
         CALL LDBREM(wlngth,plastmp,IRCODE,NIZS)
c
c delete
c
c
c         DO IR = 1,NRS
c           DO IK = 1,NKS(IR)
c             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
c             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c           ENDDO
c         ENDDO
c
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            if (iseldef.ne.0) then

               ngs = 2

               call calc_expt(iseldef,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(1) = BLABS
               ELABS(2) = 'EXPT'//DATATITLE
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif
c
         ELSE
c
               WRITE(iplot,'(a,i4,a,f8.3,a,g12.5,a,g12.5)')
     >                   'RESULT: IZ =',IZMIN,' POS = ',TOUTS(1),
     >                   'VAL = ',TVALS(1,1)
c
c           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
c
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
C
      ELSEIF (IREF.EQ.299) THEN
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION - ASSUMING NO IMPURITIES
C
C     VERIFY PARAMETERS
C
         plotr = .false.
         if (numsmooth.lt.0) then
            plotr = .true.
            numsmooth = abs(numsmooth)
         endif
c
         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
            ITEC = 2
            ISMOTH = 0
         ENDIF
C  IZMIN AND IZMAX IGNORED
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
         XLAB = 'THETA (DEGREES)'
         REF  = GRAPH(5:41)// '- NEW METHOD'
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
         IF (ATYPE.EQ.0) THEN
           MFACT = 1.0
           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
         ELSEIF (ATYPE.EQ.1) THEN
           MFACT = DTHE / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
     >                DTHE
         ELSEIF (ATYPE.EQ.2) THEN
           MFACT = 1.0 / (2.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
         ELSEIF (ATYPE.EQ.3) THEN
           MFACT = 1.0 / ( 4.0 * PI)
           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
         ENDIF
c
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 A (523.5 nm) FOR NOW
c
         WLNGTH = 5235.0
C
C  LOAD INTEGRAND - NO IMPURITIES - USE NIZS=-1
C
         CALL LDBREM(wlngth,plastmp,IRCODE,-1)
c
c delete
c
c
c         DO IR = 1,NRS
c           DO IK = 1,NKS(IR)
c             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
c             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c           ENDDO
c         ENDDO
c
         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
         WRITE (IPLOT,9012) NPLOTS,REF
         WRITE (IPLOT,*) NVIEW
         WRITE (IPLOT,*) PLANE
         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
         if (plotr) then
            call adjustout(touts,numthe,zadj,robs,zobs)
            themin = touts(1)
            themax = touts(numthe)
            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
         IF (NUMTHE.GT.1) THEN
c
            if (iseldef.ne.0) then

               ngs = 2

               call calc_expt(iseldef,touts,tvals,maxthe,numthe,
     >                   themin,themax,maxngs,ngs,datatitle)

               ELABS(1) = BLABS
               ELABS(2) = 'EXPT'//DATATITLE
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,ngs,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            else
c
               CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
     >               TABLE,IOPT,2,1.0,0)
c
            endif
c
         ELSE
           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
         ENDIF
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
C
C
c     OLD BREM
C
C     LINE-INTEGRATED BREMSSTRAHLUNG EMISSION - ZEFF = 1.0
C
C     VERIFY PARAMETERS
C
c         plotr = .false.
c         if (numsmooth.lt.0) then
c            plotr = .true.
c            numsmooth = abs(numsmooth)
c         endif
c
c         IF (NUMSMOOTH.GE.3) THEN
C
C           ALLOW ONLY ODD NUMBERS FOR SMOOTHING
C
c            IF (FLOAT(NUMSMOOTH/2).EQ.NUMSMOOTH/2.0)
c     >         NUMSMOOTH = NUMSMOOTH+1
C           WRITE(6,*) 'NUMBER FOR SUM-AVERAGE SMOOTHING',NUMSMOOTH
c            ITEC = 2
c            ISMOTH = 0
c         ENDIF
C  IZMIN AND IZMAX IGNORED
c         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
C
C  CREATE THETA VECTOR
C
c         DO II = 1, NUMTHE
c           TOUTS(II) = THEMIN + (II-1)*DTHE
c           TWIDS(II) = THERES
c         ENDDO
c         THEMAX = TOUTS(NUMTHE)
C
C  LABELS AND SCALE FACTORS
C
c         YLAB = 'BREMS. (PH M-2 S-1 NM-1)'
c         XLAB = 'THETA (DEGREES)'
c         REF  = GRAPH(5:41)
c         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
c         WRITE (PLANE,'(''RESOLUTION='',F12.6,'' DEGREES'')') THERES
c         IF (ATYPE.EQ.0) THEN
c           MFACT = 1.0
c           WRITE(ANLY,'(''NO SCALE FACTOR APPLIED'')')
c         ELSEIF (ATYPE.EQ.1) THEN
c           MFACT = DTHE / (2.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = '',G12.6,'' / (2*PI)'')')
c     >                DTHE
c         ELSEIF (ATYPE.EQ.2) THEN
c           MFACT = 1.0 / (2.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = 1 / (2*PI)'') ')
c         ELSEIF (ATYPE.EQ.3) THEN
c           MFACT = 1.0 / ( 4.0 * PI)
c           WRITE(ANLY,'(''SCALE FACTOR = 1 / (4*PI)'')')
c           YLAB = 'BREMS. (PH M-2 S-1 NM-1 SR-1)'
c         ENDIF
C  ASSUME BREMSSTRAHLUNG IS MEASURED AT 5235 NM FOR NOW
c         WLNGTH = 523.5
C
C  LOAD INTEGRAND
C
c         DO IR = 1,NRS
c           DO IK = 1,NKS(IR)
c             GFF = 1.456*(2.368/KTEBS(IK,IR))**(-0.149)
c
c            Assumed ZEFF = 1.0 everywhere
c
c             PLASTMP(IK,IR) = 9.55E-21*   1.0   *GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c
c             PLASTMP(IK,IR) = 9.55E-21*ZEFFS(IK,IR,3)*GFF
c     >                       *(RIZB*KNBS(IK,IR))**2
c     >                       /(SQRT(KTEBS(IK,IR))*WLNGTH)
c
c           ENDDO
c         ENDDO
c         CALL RVALKR (CVALSA,PLASTMP,1,1,1,FT,FP,MFACT,
c     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
C
C  INTEGRATE
C
c         CALL LOSINT(TVALS(1,1),TOUTS,TWIDS,NUMTHE,
c     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
C
c         WRITE (IPLOT,9012) NPLOTS,REF
c         WRITE (IPLOT,*) NVIEW
c         WRITE (IPLOT,*) PLANE
c         WRITE (IPLOT,*) ANLY
c
c         Adjust TOUTS to against R instead of Theta if required.
c
c         if (plotr) then
c            call adjustout(touts,numthe,zadj,robs,zobs)
c            themin = touts(1)
c            themax = touts(numthe)
c            write(XLAB,'(''RADIUS (M) AT Z='',f6.3)') zadj
c         endif
C
C  DON'T PLOT SINGLE LINES OF SIGHT
C
c         IF (NUMTHE.GT.1) THEN
c           CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,1,
c     >               ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
c     >               JOB,TITLE,XLAB,YLAB,BLABS,REF,NVIEW,PLANE,
c     >               TABLE,IOPT,2,1.0,0)
c         ELSE
c           WRITE(iplot,*) TOUTS(1),TVALS(1,1)
c         ENDIF
C
c         ITEC = 1
c         ISMOTH = 99
c         ANLY = ' '
C
C
      endif

      return 
c
c     Format statements    
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)
 9040 FORMAT(1X,'IZ, WAVELENGTH: ',A)


      end

