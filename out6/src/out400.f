      subroutine out400(iref,graph,iopt,ierr)
      use mod_params
      use mod_outcom
      use mod_cgeom
      use mod_comtor
      use mod_dynam2
      use mod_dynam3
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
c      include 'dynam4'
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
c
      integer ik,ir,iz,it
      integer in,ig,ii
c
c     Local Variables
c


c
      real plrpad(maxnks,maxnrs)
c
c     AUG variables and data
c
c
      REAL augOUTS(MAXTHE),augWIDS(MAXTHE),augVALS(MAXTHE,MAXNGS)
C     IPP/00 - Krieger: added for new plot option
      REAL augOUTS2a(MAXTHE)
      real mslope,bint,mslope2,bint2,augplato(4),deltar,deltaz
      real rintpt,zintpt,tdisp,augwid,rpt,zpt
C     IPP/00 - Krieger: added for new plot option
      real augsite(16,4), augsite2a(39,4)

      character*80 augform

      integer init_augdata,i
      data init_augdata /0/

C IPP/00 - Krieger: added for new plot option
C the variables augsite and augouts refer to a special spectrometer
C with 16 lines of sight, which was installed in the ASDEX Upgrade
C Divertor I configuration

C for the Divertor IIa configuration a similar spectrometer provides
C 19 LOS for the inner half and 20 LOS for the outer half of the
C Lyra shaped divertor.
C We introduce a new plot option for these and store the coordinates in
C the variables augsite2a and augouts2a respectively (dimension 39).
C New coordinates for the Divertor IIb configuration will be kept
C in variables augsite2b and augouts2b

      data ((augsite(ig,in),in=1,4),ig=1,16)
     >   /1.84420, -1.05390, 1.02100, -0.81790,
     >    1.84620, -1.05470, 1.01660, -0.80990,
     >    1.84660, -1.05480, 1.01700, -0.80400,
     >    1.89450, -1.05910, 0.87850, -0.77350,
     >    1.89860, -1.06160, 0.87380, -0.73920,
     >    1.89900, -1.06220, 0.88250, -0.70300,
     >    1.95270, -1.06740, 0.73350, -0.64500,
     >    1.95440, -1.07190, 0.74280, -0.58150,
     >    1.95490, -1.07330, 0.75770, -0.51970,
     >    1.83330, -1.05520, 1.14130, -0.50840,
     >    1.83310, -1.05590, 1.16310, -0.43350,
     >    1.85650, -1.06300, 1.11570, -0.33140,
     >    1.85580, -1.06500, 1.14340, -0.23700,
     >    1.88190, -1.07350, 1.09270, -0.11510,
     >    1.88000, -1.07660, 1.12800, -0.00220,
     >    1.91360, -1.08840, 1.05920,  0.14480/
c

      data ((augsite2a(ig,in),in=1,4),ig=1,39)
     >   /1420.3, -1176.9, 1242.0, -1147.0,   ! inner divertor LOS
     >    1420.3, -1176.9, 1244.5, -1138.5,
     >    1420.3, -1176.9, 1247.0, -1130.0,
     >    1420.3, -1176.9, 1252.0, -1113.0,
     >    1420.3, -1176.9, 1254.5, -1104.6,
     >    1420.3, -1176.9, 1257.0, -1096.1,
     >    1420.3, -1176.9, 1259.5, -1087.6,
     >    1420.3, -1176.9, 1262.0, -1079.1,
     >    1420.3, -1176.9, 1264.5, -1070.6,
     >    1445.0, -1127.5, 1267.0, -1062.1,
     >    1445.0, -1127.5, 1269.6, -1053.6,
     >    1445.0, -1127.5, 1272.1, -1045.2,
     >    1445.0, -1127.5, 1274.6, -1036.7,
     >    1445.0, -1127.5, 1277.1, -1028.2,
     >    1445.0, -1127.5, 1279.6, -1019.7,
     >    1445.0, -1127.5, 1282.1, -1011.2,
     >    1445.0, -1127.5, 1284.6, -1002.7,
     >    1445.0, -1127.5, 1287.1, -994.2,
     >    1445.0, -1127.5, 1289.6, -985.7,
     >    1523.9, -1181.6, 1723.6, -1193.0,   ! outer divertor LOS
     >    1523.9, -1181.6, 1721.5, -1184.0,
     >    1523.9, -1181.6, 1719.6, -1176.0,
     >    1523.9, -1181.6, 1717.4, -1166.0,
     >    1523.9, -1181.6, 1715.6, -1159.0,
     >    1523.9, -1181.6, 1713.5, -1150.0,
     >    1523.9, -1181.6, 1711.6, -1141.5,
     >    1523.9, -1181.6, 1709.8, -1134.0,
     >    1523.9, -1181.6, 1707.9, -1126.0,
     >    1523.9, -1181.6, 1705.8, -1117.0,
     >    1504.7, -1140.5, 1703.7, -1108.0,
     >    1504.7, -1140.5, 1701.6, -1099.0,
     >    1504.7, -1140.5, 1699.7, -1091.0,
     >    1504.7, -1140.5, 1697.6, -1082.0,
     >    1504.7, -1140.5, 1695.7, -1074.0,
     >    1504.7, -1140.5, 1693.8, -1066.0,
     >    1504.7, -1140.5, 1691.6, -1056.5,
     >    1504.7, -1140.5, 1689.7, -1048.5,
     >    1504.7, -1140.5, 1687.5, -1039.0,
     >    1504.7, -1140.5, 1685.5, -1030.5/
c
      data (augouts(in),in=1,16) / 1.5, 3.0, 4.5, 9.0,
c
c     The rest are (n-1)**2 mm where n is the number of the LOS
c
     >  16.0, 25.0, 36.0, 49.0, 64.0, 81.0, 100.0, 121.0, 144.0,
     > 169.0, 196.0, 225.0/
c
c
c     the values below are old and possibly not accurate, I replace
c     them by values obtained from Stefan Bosch
c     augplato(3/4)=r,z of beginning of outer target plate (left edge)
c     augplato(1/2)=r,z of another point on the target plate, I think
c                   Dave used a point projected to the radius of the
c                   BLS mirror, I used the last point in Steve Boschs
c                   file...
c     data (augplato(ig),ig=1,4) /2.7490, -1.3125, 1.4815, -0.9525/
c     in daves notation:          2.7490, -1.3164, 1.4870, -0.9540
      data (augplato(ig),ig=1,4) /1.8234, -1.0506, 1.4870, -0.9540/
C

C   IPP/00 - Krieger: added for new plot option
C - scale augsite2a from mm to m -> divide by 1000
C - preset augouts2a with incremented numbers

      if (init_augdata.eq.0) then 
         do i=1,39
           augouts2a(i) = i
           do in=1,4
             augsite2a(i, in) = augsite2a(i, in)/1000.0
           enddo
         enddo
         init_augdata = 1
      endif
c
      IF (IREF.LT.500) THEN
        ref = graph(5:41)
c
c       Series 400 are for ASDEX specific plots ... the meanings of
c       some of the variables may differ.
c
        CALL RDG2 (GRAPH2,ROBS,ZOBS,THEMIN,DTHE,THERES,NUMTHE,
     >              IZMIN,IZMAX,AVPTS,NUMSMOOTH,ATYPE,IERR)
c
c       Make adjustment for ASDEX UPGRADE plots ... their angular
c       coordinate system is 180 degrees out of phase with JET.
c       i.e. Angles measured from -ve R-axis ... so add 180 to THEMIN
c       for ASDEX plots.
c
        if (cgridopt.eq.3) themin = themin + 180.0
c
        IF (IERR.EQ.1) THEN
           WRITE(6,*) 'RDG2 ERROR READING 400 SERIES- GRAPH DETAILS'
           IERR = 0
           RETURN
        ENDIF
        IF (IOPT.EQ.0) RETURN
c
c       Load ADAS data for 400+ series plots that require it
c
c      IPP/01 Krieger - added plot option 446 for AUG DIVII spectr.
c
        IF (IREF.EQ.446) then 
           CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     >             ISELE,ISELR,ISELX,ISELD,IERR)
           IF (IERR.NE.0) THEN
             WRITE(6,*) 'ERROR READING ADAS DETAILS, IREF = ',IREF
             IERR = 0
             return
           ENDIF
        ENDIF
      endif

      call init_plot(iref,graph,iopt) 

c
c     Reserve the 400 series plots for ASDEX UPGRADE related plots - i.e
c     specific diagnostics as opposed to differing lines of sight.
c
      IF (IREF.EQ.401) THEN
C
C     CALCULATION BASED ON ION DENSITIES
C
C     VERIFY PARAMETERS THEN CALL INTEGRATION ROUTINE
C
         YLAB = 'SCALED DENSITY'
         XLAB = 'DISTANCE ALONG TARGET(M)'
         REF  = GRAPH(5:41)
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 400 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R='',F12.6,'' Z='',F12.6)') ROBS,ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN+2,IZMAX+2,SDLIMS,MAXIZS+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2400 IZ = IZMIN,IZMAX
           DO 2400 IT = 1,NUMTHE
             IF (ABSFAC.GT.0.0)
     >         TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) *ABSFAC
2400     CONTINUE
C
c
c        Adjust the coordinate system of the result to correspond to dis
c        along the ASDEX UPGRADE outer target - starting from the innerm
c        point as S=0.
c
         deltar = (augplato(3) - augplato(1))
         deltaz = (augplato(4) - augplato(2))
         rpt = augplato(3)
         zpt = augplato(4)
         mslope= deltaz/deltar
     >
         bint = augplato(4) - mslope * augplato(3)
c
         do 2405 it = 1,numthe
            if (touts(it).eq.90.0.or.touts(it).eq.-90.0) then
c
c              Deal with vertical LOS
c
               mslope2 = tan(degrad* touts(it))
               bint2 = zobs - mslope2 * robs
               rintpt = robs
               zintpt = mslope * rintpt + bint
               tdisp = sqrt((rpt-rintpt)**2+(zpt-zintpt)**2)
            else
c
c              Deal with all other cases
c
               mslope2 = tan(degrad* touts(it))
               bint2 = zobs - mslope2 * robs
               rintpt = (bint2-bint)/(mslope-mslope2)
               zintpt = mslope * rintpt + bint
               tdisp = sqrt((rpt-rintpt)**2+(zpt-zintpt)**2)
            endif
            if ((rintpt-rpt).lt.0.0) tdisp = -tdisp
            touts(it) = tdisp
c            write(6,*) 'out400:',it,tvals(it,1),touts(it),twids(it),
c     >         rpt,zpt,rintpt,zintpt,mslope,mslope2,bint,bint2
 2405    continue
c
c
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,IZMAX-IZMIN+1,
     >             ISMOTH,touts(1),touts(numthe),0.0,HI,IGNORS,ITEC,
     >             AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,ZLABS(IZMIN),REF,NVIEW,PLANE,
     >             TABLE,IOPT,2,1.0,0)
C
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
      ELSEIF (IREF.EQ.411) THEN
C
C     CALCULATION BASED ON PLRPS
C
         YLAB = 'SCALED PLRPS'
         XLAB = 'DISTANCE ALONG TARGET(M)'
         REF  = GRAPH(5:41)
C        REF  = 'PLRP LOS PLOT'
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 410 PLOT'
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
     >               ATYPE,themin,themax,
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
         DO 2410 IZ = IZMIN,IZMAX
           DO 2410 IT = 1,NUMTHE
             IF (ABSFAC.GT.0.0)
     >         TVALS(IT,IZ-IZMIN+1) = TVALS(IT,IZ-IZMIN+1) *ABSFAC
2410     CONTINUE
C
c
c        Adjust the coordinate system of the result to correspond to dis
c        along the ASDEX UPGRADE outer target - starting from the innerm
c        point as S=0.
c
         deltar = (augplato(3) - augplato(1))
         deltaz = (augplato(4) - augplato(2))
         rpt = augplato(3)
         zpt = augplato(4)
         mslope= deltaz/deltar
     >
         bint = augplato(4) - mslope * augplato(3)
c
         do 2415 it = 1,numthe
            if (touts(it).eq.90.0.or.touts(it).eq.-90.0) then
c
c              Deal with vertical LOS
c
               mslope2 = tan(degrad* touts(it))
               bint2 = zobs - mslope2 * robs
               rintpt = robs
               zintpt = mslope * rintpt + bint
               tdisp = sqrt((rpt-rintpt)**2+(zpt-zintpt)**2)
            else
c
c              Deal with all other cases
c
               mslope2 = tan(degrad* touts(it))
               bint2 = zobs - mslope2 * robs
               rintpt = (bint2-bint)/(mslope-mslope2)
               zintpt = mslope * rintpt + bint
               tdisp = sqrt((rpt-rintpt)**2+(zpt-zintpt)**2)
            endif
            if ((rintpt-rpt).lt.0.0) tdisp = -tdisp
            touts(it) = tdisp
c            write(6,*) 'out410:',it,tvals(it,1),touts(it),twids(it),
c     >         rpt,zpt,rintpt,zintpt,mslope,mslope2,bint,bint2
 2415    continue
c
c        For Asdex Upgrade only - though unlikely anyone else would
c        use this plot
c
         if (cgridopt.eq.3) then
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(51,'(1x,a)')
     >   'Data from Plot Option #411 (BLS - Fine Resolution)'
         write(51,'(1x,a,i4.4)') 'Number of points: ',numthe
         write(51,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a11,00(2x,a6,5x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
*        write(51,'(1x,a)') augform
         write(51,augform) 'Target  (m)',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f7.4,3x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
*        write(51,'(1x,a)') augform
         do 2416 it=1,numthe
           write(51,augform) touts(it),(tvals(it,iz),iz=1,izmax-izmin+1)
 2416    continue
c
         endif
c
c
         CALL DRAW(TOUTS,TWIDS,TVALS,MAXTHE,NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,touts(1),touts(numthe),0.0,HI,IGNORS,
     >             ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0,0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
      ELSEIF (IREF.EQ.421) THEN
C
C     CALCULATION BASED ON PLRPS
C
c     This plot is an odd ASDEX diagnostic involving 16 separate
c     lines of sight, which cross the volume in front of the outer
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
c
         YLAB = 'SCALED PLRPS'
         XLAB = 'VERT. DISP. Y (MM)'
         REF  = GRAPH(5:41)
C        REF  = 'PLRP LOS PLOT'
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 420 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
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
c
         PSWITCH = .TRUE.
         PSHIFT = PNCNT
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts" coordinates for each chord.
c
c
         do 2420 in = 1,16
c
c         Calculate line
c
          deltar = (augsite(in,3) - augsite(in,1))
          deltaz = (augsite(in,4) - augsite(in,2))
          mslope= deltaz/deltar
          bint = augsite(in,2) - mslope * augsite(in,1)
c          write(6,*) 'in:',in,(augsite(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 420: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif
c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
c          write(6,*) 'out420a:', robs,zobs,themin*raddeg,bint,mslope
c          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,
     >               MAXPLRP+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
c
c        Copy necessary values into aug arrays
c
         do 2425 iz = izmin,izmax
           do 2425 ig = 1,numthe
             augvals((in-1)*numthe+ig,iz-izmin+1) = tvals(ig,iz-izmin+1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
 2425    continue
c
c        Need to add a section that calculates the actual perpendicular
c        positions based on the position of the plate.
c
 2420   continue
c
        PSWITCH = .FALSE.
        PSHIFT = 1
C
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2430 IZ = IZMIN,IZMAX
           DO 2430 IT = 1,16*NUMTHE
             IF (ABSFAC.GT.0.0)
     >         augVALS(IT,IZ-IZMIN+1) = augVALS(IT,IZ-IZMIN+1) *ABSFAC
c             write (6,*) 'aug:',it,iz,augvals(it,iz),absfac
2430     CONTINUE
c
         do 2435 it = 1,16
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts(it) - augouts(it-1)
           endif

           do 2440 in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
2440       continue
c           write(6,*) 'out421:',it,numthe,augouts(it),augwids(it)
2435     continue
c
         themin = augouts(1)
         themax = augouts(16)
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(52,'(1x,a)')
     >   'Data from Plot Option #421 (DIV - Spectrometer)'
         write(52,'(1x,a,i4.4)') 'Number of points: ',16*numthe
         write(52,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a12,00(3x,a6,4x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         write(52,augform) 'Pl-Dist (mm)',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f6.2,6x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         do 2436 it=1,16*numthe
           write(52,augform) augOUTS(it),
     >                      (augVALS(it,iz),iz=1,izmax-izmin+1)
 2436    continue

C
         CALL DRAW(augOUTS,augWIDS,augvals,MAXTHE,16*NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0,0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
      ELSEIF (IREF.EQ.422) THEN
C IPP/00 - Krieger:
C added this new plot option to plot results of DivIIa divertor
C spectrometer
C
C     CALCULATION BASED ON PLRPS
C
c     This plot is an even odder ASDEX diagnostic involving 39 separate
c     lines of sight, which cross the volume perpendicularly in front of
C     the inner and outer (now vertically oriented)
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
c
         YLAB = 'SCALED PLRPS'
         XLAB = 'Spectrometer Chord #'
         REF  = GRAPH(5:41)
C        REF  = 'PLRP LOS PLOT'
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 420 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
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
c
         PSWITCH = .TRUE.
         PSHIFT = PNCNT
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts2a" coordinates for each chord.
c
c
         do 2421 in = 1,39
c
c         Calculate line
c
          deltar = (augsite2a(in,3) - augsite2a(in,1))
          deltaz = (augsite2a(in,4) - augsite2a(in,2))
          mslope= deltaz/deltar
          bint = augsite2a(in,2) - mslope * augsite2a(in,1)
c          write(6,*) 'in:',in,(augsite2a(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 422: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif
c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
c          write(6,*) 'out422a:', robs,zobs,themin*raddeg,bint,mslope
c          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,
     >               MAXPLRP+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
c
c        Copy necessary values into aug arrays
c
         do iz = izmin,izmax
           do ig = 1,numthe
             augvals((in-1)*numthe+ig,iz-izmin+1) = tvals(ig,iz-izmin+1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
           enddo
         enddo 
c
c        Need to add a section that calculates the actual perpendicular
c        positions based on the position of the plate.
c
 2421   continue
c
        PSWITCH = .FALSE.
        PSHIFT = 1
C
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO IZ = IZMIN,IZMAX
           DO IT = 1,39*NUMTHE
             IF (ABSFAC.GT.0.0)
     >         augVALS(IT,IZ-IZMIN+1) = augVALS(IT,IZ-IZMIN+1) *ABSFAC
c             write (6,*) 'aug:',it,iz,augvals(it,iz),absfac
           ENDDO
         ENDDO
c
         do it = 1,39
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts2a(it) - augouts2a(it-1)
           endif

           do in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
           enddo
c           write(6,*) 'out422:',it,numthe,augouts2a(it),augwids(it)
         enddo
c
         themin = augouts2a(1)
         themax = augouts2a(39)
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(52,'(1x,a)')
     >   'Data from Plot Option #422 (DIV - Spectrometer)'
         write(52,'(1x,a,i4.4)') 'Number of points: ',39*numthe
         write(52,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a12,00(3x,a6,4x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         write(52,augform) 'Channel Num.',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f6.2,6x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         do it=1,39*numthe
           write(52,augform) augouts2a(it),
     >                      (augVALS(it,iz),iz=1,izmax-izmin+1)
         enddo

C
         CALL DRAW(augouts2a,augWIDS,augvals,MAXTHE,39*NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
      ELSEIF (IREF.EQ.431) THEN
C
C     CALCULATION BASED ON DENSITY
C
c     This plot is an odd ASDEX diagnostic involving 16 separate
c     lines of sight, which cross the volume in front of the outer
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
c
         YLAB = 'SCALED DENSITY'
         XLAB = 'VERT. DISP. Y (MM)'
         REF  = GRAPH(5:41)
C        REF  = 'DENSITY LOS PLOT'
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 430 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts" coordinates for each chord.
c
c
         do 2450 in = 1,16
c
c         Calculate line
c
          deltar = (augsite(in,3) - augsite(in,1))
          deltaz = (augsite(in,4) - augsite(in,2))
          mslope= deltaz/deltar
          bint = augsite(in,2) - mslope * augsite(in,1)
c          write(6,*) 'in:',in,(augsite(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 431: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif
c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
c          write(6,*) 'out420a:', robs,zobs,themin*raddeg,bint,mslope
c          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,nizs+3,
     >               IZMIN+2,IZMAX+2,SDLIMS,
     >               MAXizs+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
c
c        Copy necessary values into aug arrays
c
         do 2455 iz = izmin,izmax
           do 2455 ig = 1,numthe
             augvals((in-1)*numthe+ig,iz-izmin+1) = tvals(ig,iz-izmin+1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
 2455    continue
c
c        Need to add a section that calculates the actual perpendicular
c        positions based on the position of the plate.
c
 2450   continue
c
C
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO 2460 IZ = IZMIN,IZMAX
           DO 2460 IT = 1,16*NUMTHE
             IF (ABSFAC.GT.0.0)
     >         augVALS(IT,IZ-IZMIN+1) = augVALS(IT,IZ-IZMIN+1) *ABSFAC
c             write (6,*) 'aug:',it,iz,augvals(it,iz),absfac
2460     CONTINUE
c
         do 2465 it = 1,16
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts(it) - augouts(it-1)
           endif

           do 2470 in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
2470       continue
c           write(6,*) 'out420:',it,numthe,augouts(it),augwids(it)
2465     continue
c
         themin = augouts(1)
         themax = augouts(16)
C
         CALL DRAW(augOUTS,augWIDS,augvals,MAXTHE,16*NUMTHE,ANLY,
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



      ELSEIF (IREF.EQ.441) THEN
C IPP/01 - Krieger:
C added this new plot option to plot results of DivIIa divertor
C spectrometer
C
C     CALCULATION BASED ON PLRPS
C
c     This plot is an even odder ASDEX diagnostic involving 39 separate
c     lines of sight, which cross the volume perpendicularly in front of
C     the inner and outer (now vertically oriented)
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
c
		 YLAB = 'SCALED PLRPS'
         XLAB = 'Spectrometer Chord #'
         REF  = GRAPH(5:41)
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
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 420 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
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
c
         PSWITCH = .TRUE.
         PSHIFT = PNCNT
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts2a" coordinates for each chord.
c
c
         do in = 1,39
c
c         Calculate line
c
          deltar = (augsite2a(in,3) - augsite2a(in,1))
          deltaz = (augsite2a(in,4) - augsite2a(in,2))
          mslope= deltaz/deltar
          bint = augsite2a(in,2) - mslope * augsite2a(in,1)
c          write(6,*) 'in:',in,(augsite2a(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 441: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif
c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
c          write(6,*) 'out441a:', robs,zobs,themin*raddeg,bint,mslope
c          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,PLRPCNT+3,
     >               IZMIN+2,IZMAX+2,PLRPS,
     >               MAXPLRP+2,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
c
c        Copy necessary values into aug arrays
c
         do iz = izmin,izmax
           do ig = 1,numthe
             augvals((in-1)*numthe+ig,iz-izmin+1) = tvals(ig,iz-izmin+1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
           enddo
         enddo 
c
        enddo
c
        PSWITCH = .FALSE.
        PSHIFT = 1
C
C
C        MULTIPLY BY MAGNITUDE FACTOR
C
         DO IZ = IZMIN,IZMAX
           DO IT = 1,39*NUMTHE
             IF (ABSFAC.GT.0.0)
     >         augVALS(IT,IZ-IZMIN+1) = augVALS(IT,IZ-IZMIN+1) *ABSFAC
c             write (6,*) 'aug:',it,iz,augvals(it,iz),absfac
           ENDDO
         ENDDO
c
         do it = 1,39
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts2a(it) - augouts2a(it-1)
           endif

           do in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
           enddo
c           write(6,*) 'out421:',it,numthe,augouts2a(it),augwids(it)
         enddo
c
         themin = augouts2a(1)
         themax = augouts2a(39)
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(59,'(1x,a)')
     >   'Data from Plot Option #441 (DIV II - Spectrometer)'
         write(59,'(1x,a,i4.4)') 'Number of points: ',39*numthe
         write(59,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a12,00(3x,a6,4x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
         write(59,augform) 'Channel Num.',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f6.2,6x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
         do it=1,39*numthe
           write(59,augform) augouts2a(it),
     >                      (augVALS(it,iz),iz=1,izmax-izmin+1)
         enddo

C
         CALL DRAW(augouts2a,augWIDS,augvals,MAXTHE,39*NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABS(IZMIN),REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '

      ELSEIF (IREF.EQ.446) THEN
c
C IPP/01 - Krieger:
C added this new plot option to plot results of DivIIa divertor
C spectrometer
C
C     CALCULATION BASED ON PLRPS from ADAS
C
c     This plot is an even odder ASDEX diagnostic involving 39 separate
c     lines of sight, which cross the volume perpendicularly in front of
C     the inner and outer (now vertically oriented)
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
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
C  ONLY TOTAL NEUTRAL CONTRIBUTIONS ALLOWED FOR ADAS PRLP INTEGRALS
                 IF (IZMIN.LT.0) IZMIN = 0
C  IZMAX IGNORED FOR ADAS PLRP INTEGRALS
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0

         XLAB = 'Spectrometer Chord #'
         REF  = GRAPH(5:41)
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
          call LDADAS(CION,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,
     >                ISELX,plrpad,Wlngth,IRCODE)
c
          IF (IRCODE.NE.0) THEN
             WRITE(6,*) 'SPEC ERROR, IRCODE = ',IRCODE
             return
          ENDIF
C
C  LOAD INTEGRAND
C
         CALL RVALKR (CVALSA,PLRPAD,1,1,1,FT,FP,MFACT,
     >     XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 440 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         PSWITCH = .TRUE.

         PSHIFT = PNCNT
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts2a" coordinates for each chord.
c
c
         do in = 1,39
c
c         Calculate line
c
          deltar = (augsite2a(in,3) - augsite2a(in,1))
          deltaz = (augsite2a(in,4) - augsite2a(in,2))
          mslope= deltaz/deltar
          bint = augsite2a(in,2) - mslope * augsite2a(in,1)
c          write(6,*) 'in:',in,(augsite2a(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 446: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif

c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
          write(6,*) 'out446a:', robs,zobs,themin*raddeg,bint,mslope
          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
C  INTEGRATE
C
C  CREATE THETA VECTOR
C
         DO II = 1, NUMTHE
           TOUTS(II) = THEMIN + (II-1)*DTHE
           TWIDS(II) = THERES
         ENDDO
         THEMAX = TOUTS(NUMTHE)
C
         CALL LOSINT(TVALS,TOUTS,TWIDS,NUMTHE,
     >               ROBS,ZOBS,AVPTS,CVALSA,0.0,0)
c
c        Copy necessary values into aug arrays
c
           do ig = 1,numthe
             augvals((in-1)*numthe+ig,1) = tvals(ig,1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
           enddo
c
        enddo
c
        PSWITCH = .FALSE.
        PSHIFT = 1

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
C
         do it = 1,39
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts2a(it) - augouts2a(it-1)
           endif

           do in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
           enddo
c           write(6,*) 'out446:',it,numthe,augouts2a(it),augwids(it)
         enddo
c
         themin = augouts2a(1)
         themax = augouts2a(39)
c
c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(59,'(1x,a)')
     >   'Data from Plot Option #446 (DIV II - Spectrometer)'
         write(59,'(1x,a,i4.4)') 'Number of points: ',39*numthe
         write(59,'(1x,a,i4.4)') 'Number of curves: ',1
         augform='(1x,a12,00(3x,f6.1,4x))'
         write(augform( 9:10),'(i2)') 1
         write(59,augform) 'Wave length ',wlngth/10.
         augform='(1x,f6.2,6x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') 1
         do it=1,39*numthe
           write(59,augform) augouts2a(it), augVALS(it,1)


         enddo

C
         CALL DRAW(augouts2a,augWIDS,augvals,MAXTHE,39*NUMTHE,ANLY,
     >             1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,PLABAD,REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0,0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
      ELSEIF (IREF.EQ.491) THEN
C
C     CALCULATION BASED ON PINALPHA
C
c     This plot is an odd ASDEX diagnostic involving 16 separate
c     lines of sight, which cross the volume in front of the outer
c     target. The data for this lines ... i.e. start and end points
c     is loaded in the beginning ... these are then used to calculate
c     individual line integrals that are combined onto one plot.
c
         YLAB = 'SCALED H-alpha'
         XLAB = 'VERT. DISP. Y (MM)'
         REF  = GRAPH(5:41)
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
         IZMIN = 1
         IZMAX = 1
         IF (ATYPE.LT.0 .OR. ATYPE.GT.3) ATYPE = 0
c
c        For these plots THERES is interpreted as DRAD
c
         DRAD = THERES
c
         IF (DRAD.EQ.0.0) THEN
            DRAD = 0.001
            WRITE(6,*) 'ERROR: DRAD SPECIFIED AS ZERO FOR 420 PLOT'
         ENDIF
         WRITE (NVIEW,'(''R= VARIES     Z='',F12.6)') ZOBS
         WRITE (PLANE,'(I4,'' POINTS WITH '',I4,'' PT. AVERAGING'')')
     >          NUMTHE,AVPTS
C
C        NOTE: IF SECONDARY NEUTRAL PLRP PLOTS ARE REQUESTED
C              THE PLOT WILL BE GENERATED ONLY FOR THE FIRST
C              NEUTRAL LINE.
C
         PSWITCH = .TRUE.
         PSHIFT = PNCNT
c
c        These plots are implemented based on the pre-loaded
c        data from the chord end-points. Currently, an arbitrary
c        effective viewing position based on the intersection of
c        viewing chord with the line Z = -ZOBS is calculated
c        and then the angle of observation (THEMIN) is calculated
c        from the equation of the line - taking into account DTHE
c        which is assumed to be the effective instrument angular
c        viewing cone. Then the value of the LOS is calculated for
c        a number of points - usually equal to 1, with averaging over
c        smaller points in between. These points are then linked togethe
c        and plotted versus the "augouts" coordinates for each chord.
c
c
         do in = 1,16
c
c         Calculate line
c
          deltar = (augsite(in,3) - augsite(in,1))
          deltaz = (augsite(in,4) - augsite(in,2))
          mslope= deltaz/deltar
          bint = augsite(in,2) - mslope * augsite(in,1)
c          write(6,*) 'in:',in,(augsite(in,ig),ig=1,4),deltar,deltaz,
c     >            atan2c(deltaz,deltar)*raddeg,mslope,bint
c
c         If MSLOPE = 0 this is obviously an error since the observation
c         line will never intersect the viewing position ... or will
c         do so at an infinite number of points. Set ROBS = 3.0 ... i.e.
c         large but reasonably near the torus under these conditions.
c         Issue an error message.
c
          if (mslope.eq.0.0) then
             robs = 3.0
             write (6,*) 'ERROR in plot 491: observation'//
     >                   ' line slope is ZERO'
          else
             robs = (zobs - bint) / mslope
          endif
c
          themin = atan2c(deltaz,deltar) * raddeg
          themin = themin - 0.5 * (numthe-1) * dthe
c
c          write(6,*) 'out420a:', robs,zobs,themin*raddeg,bint,mslope
c          write(6,*) 'num:',numthe,avpts,drad,atype,izmin,izmax
c
         CALL INTLOS(TVALS,TOUTS,TWIDS,ROBS,ZOBS,DRAD,NUMTHE,AVPTS,
     >               ATYPE,THEMIN,THEMAX,
     >               DTHE,NIZS+3,
     >               IZMIN,IZMAX,PINALPHA,1,
     >               SDTS,KTEBS,0,ANLY,PSWITCH,PIZS,FT,FP)
c
c        Copy necessary values into aug arrays
c
         do 2525 iz = izmin,izmax
           do 2525 ig = 1,numthe
             augvals((in-1)*numthe+ig,iz-izmin+1) = tvals(ig,iz-izmin+1)
c             write (6,*) in,ig,iz,(in-1)*numthe+ig,tvals(ig,iz-izmin+1)
c     >                   augvals((in-1)*numthe+ig,iz-izmin+1)
 2525    continue
c
c        Need to add a section that calculates the actual perpendicular
c        positions based on the position of the plate.
c
        enddo
c
        PSWITCH = .FALSE.
        PSHIFT = 1
C
C
         do it = 1,16
           if (it.eq.1) then
             augwid = 1.0
           else
             augwid = augouts(it) - augouts(it-1)
           endif

           do in = 1,numthe
             augwids((in-1)*numthe+it) = augwid/numthe
           enddo
c           write(6,*) 'out421:',it,numthe,augouts(it),augwids(it)
         enddo
c
         themin = augouts(1)
         themax = augouts(16)

c        we want all data written in a file as ASCII to compare with
c        experimental values

         write(52,'(1x,a)')
     >   'Data from Plot Option #491 (DIV - Spectrometer)'
         write(52,'(1x,a,i4.4)') 'Number of points: ',16*numthe
         write(52,'(1x,a,i4.4)') 'Number of curves: ',izmax-izmin+1
         augform='(1x,a12,00(3x,a6,4x))'
         write(augform( 9:10),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         write(52,augform) 'Pl-Dist (mm)',
     >                     (plabs(izmin+iz)(6:11),iz=0,izmax-izmin)
         augform='(1x,f6.2,6x,1p,00(2x,e11.4))'
         write(augform(16:17),'(i2)') izmax-izmin+1
*        write(52,'(1x,a)') augform
         do 2437 it=1,16*numthe
           write(52,augform) augOUTS(it),
     >                      (augVALS(it,iz),iz=1,izmax-izmin+1)
 2437    continue
c
         klab='    HAL 656'
         CALL DRAW(augOUTS,augWIDS,augvals,MAXTHE,16*NUMTHE,ANLY,
     >             IZMAX-IZMIN+1,
     >             ISMOTH,THEMIN,THEMAX,0.0,HI,IGNORS,ITEC,AVS,NAVS,
     >             JOB,TITLE,XLAB,YLAB,klab,REF,NVIEW,
     >             PLANE,
     >             TABLE,IOPT,2,1.0,0)
         ITEC = 1
         ISMOTH = 99
         ANLY = ' '
c
      endif 


 9012 FORMAT(1X,'PLOT',I3,4X,A)
 9040 FORMAT(1X,'IZ, WAVELENGTH: ',A)


      return
      end

