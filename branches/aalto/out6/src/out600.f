      subroutine out600(iref,graph,iopt,ierr)
      use mod_params
      use mod_outcom
      use mod_cgeom
      use mod_comtor
      use mod_pindata
      use mod_cedge2d
      use mod_transcoef
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
c      include 'dynam2'
c      include 'dynam3'
c      include 'dynam4'
c     include 'pindata'
c      include 'cadas'
c      include 'grbound'
c      include 'outxy'
c     include 'cedge2d'
c     include 'transcoef'
c      include 'cioniz'
c      include 'reiser' 
c      include 'printopt' 
c
c     Local Variables
c
c
      integer i,j,k
      integer ik,ir,iz,it
      integer in,ii,ig
      integer lw,uw,iw
c
c     Local Variables
c

      CHARACTER*9  VLABEL(14)
      DATA VLABEL/'DENSITY','IONIS','CX',
     >            'VX','VY','VZ','VX VX','VY VY','VZ VZ',
     >            'VX VY','VX VZ','VY VZ',
     >            'H2 IONIS','H2 DISSOC'/

      integer inc


      real tmpsum

      integer llabs(maxseg)


      REAL    TAUS,SLVALS(MAXNRS,2)
      integer midnks

c slmod begin
      integer count
c slmod end

c
        if ((iref.lt.669.or.iref.gt.680).and.iopt.eq.0) then 
           return
c slmod begin - new
c...    For the neutral particle track plot (650), the IF statement
c       following this one voided the plot if more than NIZS+1
c       particle tracks were requested:
        elseif (iref.EQ.650.AND.iopt.GT.0) then
c slmod end
        elseif (IOPT.LT.-2.OR.IOPT.GT.NIZS+1) then 
           return
        endif 
c

      call init_plot(iref,graph,iopt) 


C
C-----------------------------------------------------------------------
C
C     CONTOUR PLOTS OF PIN DATA - IONIZATION, DENSITY, HALPHA, ETC
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C     PIN - BACKGROUND IONISATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.601.OR.IREF.EQ.602) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - BG IONISATION '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINION,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - BG NEUTRAL ATOM DENSITY CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.603.OR.IREF.EQ.604) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - BG ATOM DENSITY '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINATOM,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - H-ALPHA CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.605.OR.IREF.EQ.606) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - H-ALPHA '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINALPHA,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - MOLECULAR HYDROGEN DENSITY CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.607.OR.IREF.EQ.608) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - MOL. DENSITY '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINMOL,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - NEUTRAL IMPURITY DENSITY CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.609.OR.IREF.EQ.610) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - NEUTRAL IMP. DENSITY '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINZ0,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - NEUTRAL IMPURITY IONISATION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.611 .OR.IREF.EQ.612) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'PIN - IMP. IONISATION '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINIONZ,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF

C
C-----------------------------------------------------------------------
C     PIN - MOMENTUM SOURCE CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.613 .OR.IREF.EQ.614) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = 10
        IZ     = IOPT
        REF = 'PIN - MOMENTUM SOURCE '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - ION ENERGY SOURCE CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.615 .OR.IREF.EQ.616) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = 10
        IZ     = IOPT
        REF = 'PIN - ION ENERGY SOURCE '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINQI,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - ELECTRON ENERGY SOURCE CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.617 .OR.IREF.EQ.618) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = 10
        IZ     = IOPT
        REF = 'PIN - ELEC ENERGY SOURCE '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PINQE,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF


C
C-----------------------------------------------------------------------
C     PIN - COLD NEUTRAL VELOCITY DISTRIBUTION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.621 .OR.IREF.EQ.622) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = PINVDIST(1,IOPT,IK,IR)
          ENDDO
        ENDDO
        REF = 'PIN - COLD NEUTRAL '// VLABEL(IOPT) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - HOT NEUTRAL VELOCITY DISTRIBUTION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.623 .OR.IREF.EQ.624) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = PINVDIST(2,IOPT,IK,IR)
          ENDDO
        ENDDO
        REF = 'PIN - HOT NEUTRAL '// VLABEL(IOPT) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - REFLECTED NEUTRAL VELOCITY DISTRIBUTION CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.625 .OR.IREF.EQ.626) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = PINVDIST(3,IOPT,IK,IR)
          ENDDO
        ENDDO
        REF = 'PIN - REFL. NEUTRAL '// VLABEL(IOPT) // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - HOT NEUTRAL TEMPERATURE RELATIVE TO ION TEMPERATURE
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.627 .OR.IREF.EQ.628) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        DO IR = 1,NRS
          DO IK = 1,NKS(IR)
            PLASTMP(IK,IR) = (CRMB/EMI/3.0*(PINVDIST(2,7,IK,IR)
     >                        +PINVDIST(2,8,IK,IR)+PINVDIST(2,9,IK,IR))
     >                       - KTIBS(IK,IR)) / KTIBS(IK,IR)
          ENDDO
        ENDDO
        REF = 'PIN - (THOT-TION)/TION ' // XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,PLASTMP,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - CORE PROFILES OF IONIZATION AND NEUTRALS
C-----------------------------------------------------------------------
C
      if (iref.eq.629) then
c
        ELABS(1) =  'IZ  IONIZATION'
        ELABS(2) =  'IZD IONIZ.DENS'
        ELABS(3) =  'NEUTNEUTRALS'
        ELABS(3) =  'ND  NEUT.DENS'
c
        XLAB = 'R (M) OUTER MIDPLANE'
c
        YLAB   = 'PIN DATA'
        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''PLOT OF PIN IONIZ AND NEUT FOR CORE'')')
        WRITE (IPLOT,9012) NPLOTS,REF
c
        CALL rzero (kvals, maxnks*maxngs)
c
        DO Ir = 1,irsep-1
c
          in = ir
c
          KOUTS(in) = rcouter(ir)
c
c         Approximate width
c
          kwids(in) = (kinds(oumid,ir)+koutds(oumid,ir) )/2.0
c
          KVALS(In,1) = piniz_info(ir,1)
          KVALS(In,2) = piniz_info(ir,2)
          KVALS(In,3) = piniz_info(ir,3)
          KVALS(In,4) = piniz_info(ir,4)
c
        enddo
C
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,irsep-1,ANLY,
     >    4,99,KOUTS(1),KOUTS(irsep-1),-HI,HI,IGNORS,
     >    ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
c
      ENDIF
C
C-----------------------------------------------------------------------
C
C     PLOTS OF PIN WALL PROFILES
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C     PIN - MAJOR RADIUS OF THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.630) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    RVESM'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'MAJOR RADIUS (M)'
          NPLOTS = NPLOTS + 1
          REF = 'MAJOR RADIUS OF VESSEL WALL'
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - MAJOR RADIUS OF THE WALL NEAR THE DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.631) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    RVESM'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'MAJOR RADIUS (M)'
          NPLOTS = NPLOTS + 1
          REF = 'MAJOR RADIUS OF DIVERTOR WALL'
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - HEIGHT OF THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.632) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    ZVESM'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'HEIGHT (M)'
          NPLOTS = NPLOTS + 1
          REF = 'HEIGHT OF VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = 0.5*(zvesm(1,2)+zvesm(1,1))
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = 0.5*(zvesm(i,2)+zvesm(i,1))
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - HEIGHT OF THE WALL NEAR THE DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.633) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    ZVESM'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'HEIGHT (M)'
          NPLOTS = NPLOTS + 1
          REF = 'HEIGHT OF DIVERTOR WALL'
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
          lvals(j+1,1) = 0.5*(zvesm(1,2)+zvesm(1,1))
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = 0.5*(zvesm(i,2)+zvesm(i,1))
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
          lvals(j,1) = 0.5*(zvesm(nvesm,2)+zvesm(nvesm,1))
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = 0.5*(zvesm(i,2)+zvesm(i,1))
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF HYDROGEN TO THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.634) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLUXHW'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'HYDROGEN FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'NEUTRAL FLUX TO VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = FLUXHW(1)
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = FLUXHW(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF HYDROGEN TO THE WALL NEAR THE DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.635) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLUXHW'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'HYDROGEN FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'NEUTRAL FLUX TO DIVERTOR WALL'
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
          lvals(j+1,1) = FLUXHW(1)
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = FLUXHW(I)
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
          lvals(j,1) = FLUXHW(nvesm)
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = FLUXHW(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF ATOMS+IONS TO THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.636) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW2'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'ION+ATOM FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'ION+ATOM FLUX TO VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = FLXHW2(1)
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = FLXHW2(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF ATOMS+IONS TO THE WALL NEAR THE DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.637) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW2'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'ION+ATOM FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'ION+ATOM FLUX TO DIVERTOR WALL'
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
          lvals(j+1,1) = FLUXHW(1)
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = FLXHW2(I)
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
          lvals(j,1) = FLXHW2(nvesm)
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = FLXHW2(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF IMPURITIES SPUTTERED FROM THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.638) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW3'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'Z SPUT. FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'Z SPUTTERING FROM VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = FLXHW3(1)
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = FLXHW3(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF IMPURITIES SPUTTERED FROM THE WALL NEAR DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.639) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW3'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'Z SPUT. FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'Z SPUTTERING FROM VESSEL WALL'
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
          lvals(j+1,1) = FLXHW3(1)
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = FLXHW3(I)
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
          lvals(j,1) = FLXHW3(nvesm)
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = FLXHW3(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF IMPURITIES REDEPOSITED ON THE WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.640) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW4'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'Z REDEP FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'Z REDEPOSITION ON VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = FLXHW4(1)
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = FLXHW4(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - FLUX OF IMPURITIES REDEPOSITED ON THE WALL NEAR DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.641) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW4'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'Z REDEP FLUX (M-2 S-1)'
          NPLOTS = NPLOTS + 1
          REF = 'Z REDEPOSITION ON VESSEL WALL'
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
          lvals(j+1,1) = FLXHW4(1)
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = FLXHW4(I)
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
          lvals(j,1) = FLXHW4(nvesm)
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = FLXHW4(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - AVERAGE ENERGY OF ATOMS HITTING WALL
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.642) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW5'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'H0 IMPACT ENERGY (EV)'
          NPLOTS = NPLOTS + 1
          REF = 'ATOM ENERGY ON VESSEL WALL'
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL rzero (lvals, maxseg*maxngs)
          PLTMAX = -HI
          PLTMIN = HI
c
          lwids(1) = sqrt((rvesm(1,2)-rvesm(1,1))**2
     >                   +(zvesm(1,2)-zvesm(1,1))**2)
          louts(1) = lwids(1) * 0.5
          lvals(1,1) = FLXHW5(1)
          pltmax = max(pltmax,lvals(1,1))
          pltmin = min(pltmin,lvals(1,1))
          do i = 2, nvesm
            lwids(i) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                     +(zvesm(i,2)-zvesm(i,1))**2)
            louts(i) = louts(i-1) + 0.5 * (lwids(i)+lwids(i-1))
            lvals(i,1) = FLXHW5(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - AVERAGE ENERGY OF ATOMS HITTING WALL NEAR DIVERTOR
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.643) THEN
        if (nvesm.le.0) then
          write(6,*) ' no wall data to plot!'
        else
          NAME = '    FLXHW5'
          XLAB = 'DISTANCE ALONG WALL (M)'
          YLAB   = 'H0 IMPACT ENERGY (EV)'
          NPLOTS = NPLOTS + 1
          REF = 'ATOM ENERGY ON VESSEL WALL'
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
          lvals(j+1,1) = FLXHW5(1)
          llabs(j+1) = jvesm(1)
          pltmax = max(pltmax,lvals(j+1,1))
          pltmin = min(pltmin,lvals(j+1,1))
          do i = 2, k
            ii = j + i
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                       +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii-1) + 0.5 * (lwids(ii)+lwids(ii-1))
            lvals(ii,1) = FLXHW5(I)
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
          lvals(j,1) = FLXHW5(nvesm)
          llabs(j) = jvesm(nvesm)
          pltmax = max(pltmax,lvals(j,1))
          pltmin = min(pltmin,lvals(j,1))
          do i = nvesm-1, nvesm-j+1, -1
            ii = j - (nvesm-i)
            lwids(ii) = sqrt((rvesm(i,2)-rvesm(i,1))**2
     >                      +(zvesm(i,2)-zvesm(i,1))**2)
            louts(ii) = louts(ii+1) - 0.5 * (lwids(ii)+lwids(ii+1))
            lvals(ii,1) = FLXHW5(I)
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
       ENDIF
C
C-----------------------------------------------------------------------
C     PIN - TRACKS OF SELECTED NEUTRALS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.650) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NPLOTS = NPLOTS - 1
        LW     = 1
        UW     = 0
        IG     = 1
        NAME = '     TRACK'
C
        rmax = 100.0

        count = 0
 6500   CONTINUE
        count = count + 1
        DO IW = LW, MAXNWS
c          write(0,*) 'HWALKS(IW,1)',iw,HWALKS(IW,1),rmax
          IF (HWALKS(IW,1).GE.10.0*RMAX) GOTO 6505
          UW = IW
        ENDDO
 6505   CONTINUE
        write(0,*) '650:', iw,maxnws,uw
        WRITE (REF,'(''NEUTRAL TRACK'',I3,''  ('',I5,'' PTS)'')')
     >      IG,UW-LW+1
        IF (LW.LT.MAXNWS.AND.UW.GT.LW.AND.IG.LE.IOPT) THEN
          NPLOTS = NPLOTS + 1
          WRITE (IPLOT,9012) NPLOTS,REF
          CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
          CALL SUPIMP ('SELECT')
c          IF (count.eq.1) CALL SUPIMP ('SELECT')
          CALL GRTRAC (HWALKS(LW,1),HWALKS(LW,2),UW-LW+1,NAME,'LINE',0)
          CALL FRAME
          IG = IG + 1
        ENDIF
        LW = UW + 2
c slmod begin
c Added the UW check to avoid an infinite loop when EIRENE not run. -SL, 20/01/12
        IF (UW.GT.0.AND.LW.LT.MAXNWS.AND.count.LT.1e+6) GOTO 6500
c        IF (UW.GT.0.AND.LW.LT.MAXNWS) GOTO 6500
        CALL FRAME
c
c        IF (LW.LT.MAXNWS) GOTO 6500
c slmod end
      ENDIF
C
C-----------------------------------------------------------------------
C     PIN - TRACKS OF ALL RECORDED NEUTRALS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.651) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NPLOTS = NPLOTS - 1
        LW     = 1
        UW     = 0
        IG     = 1
        NAME = '    TRACKS'
c
        WRITE (REF,'(''NEUTRAL TRACKS: '',I3)') iopt 
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >      YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL SUPIMP ('SELECT')
C
 6501   CONTINUE
        DO IW = LW, MAXNWS
          IF (HWALKS(IW,1).GE.10.0*RMAX) GOTO 6506
          UW = IW
        ENDDO
 6506   CONTINUE
        IF (LW.LT.MAXNWS.AND.UW.GT.LW.AND.IG.LE.IOPT) THEN
          CALL GRTRAC (HWALKS(LW,1),HWALKS(LW,2),UW-LW+1,NAME,'LINE',1)
          IG = IG + 1
        ENDIF
        LW = UW + 2
        IF (LW.LT.MAXNWS) GOTO 6501
        CALL FRAME
      ENDIF
c
c  660 series plots - plotting force balances
c
c
c slmod begin
C
C-----------------------------------------------------------------------
C     FORCE BALANCE FOR SELECTED CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.660.OR.IREF.EQ.661) THEN
        ELABS(1) = 'FF  FF     '
        ELABS(2) = ' TiGFTiG   '
        ELABS(3) = 'TeG FTeG   '
        ELABS(4) = '  FEFE     '
        ELABS(5) = 'NET NET    '

        XLAB = '   S  (M)'
        YLAB = '           '

        FACT = QTIM**2 * EMI / CRMI
c
c Add ionisation state info here
c
        READ (GRAPH(30:33),'(I4)') IZ
        READ (GRAPH(38:44),'(I4)') IR

        NPLOTS = NPLOTS + 1
        WRITE (REF,'(''Ring: '',I3,'', K: '',F7.4,'' ,State: '',I3)')
     +        IR,KKS(IR),IZ
c
        WRITE (IPLOT,9012) NPLOTS,REF
        write(iplot,'(2x,''IK'',2x,''IR'',5x,''FF'',5x,5x'//
     >                   ',''FIG'',4x,5x,''FEG'',4x,5x'//
     >                    ',''FE'',5x,5x,''NET'',4x,6x,''S'',5x'//
     >                    ',6x,''R'',5x,6x,''Z'')')

c
        CALL rzero (kvals, maxnks*maxngs)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2
           kouts(1) = 0.0
c
c Fix constant for density units
c
           TAUS = CRMI * KTIDS(IDDS(IR,2))**1.5 * SQRT(1.0/CRMB) /
     +            (6.8E-14 * (1 + CRMB / CRMI) * KNDS(IDDS(IR,2)) *
     +             REAL(IZ)**2.0 * RIZB**2 * 15.0)

           KWIDS(1)         = kss(1,ir)
           kouts(nks(ir)+2) = ksmaxs(ir)
           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
c
c------------
*psmod
*
c  If one wants a non-zero impurity velocity to modify the force equation
c  then you add an input for Vimp and use the modified force equation below. 
c
*  Add Impurity Velocity to the equation
*
c A)        KVALS(1,1) = AMU * CRMI * (KVDS(idds(ir,2))-Vimp) / TAUS
* B)        KVALS(nks(ir)+2,1) = AMU*CRMI*(KVDS(idds(ir,1))-Vimp) / TAUS
c C)        KVALS(IK+in,1) = AMU * CRMI * (KVHS(IK,IR)-(Vimp*QTIM))
c     +                                       / QTIM / TAUS
*psmod
c-------------
c

           KVALS(1,1) = AMU * CRMI * KVDS(idds(ir,2)) / TAUS

           KVALS(1,2) = KFIdS(idds(ir,2)) * KBETAS(IZ) * ECH / FACT

           KVALS(1,3) = KFEdS(idds(ir,2)) * KALPHS(IZ) * ECH / FACT
           KVALS(1,4) = IZ * KEDS(idds(ir,2)) * ECH

           KVALS(1,5) = KVALS(1,1) + KVALS(1,2)
     +                + KVALS(1,3) + KVALS(1,4)

c
           TAUS = CRMI * KTIDS(IDDS(IR,1))**1.5 * SQRT(1.0/CRMB) /
     +            (6.8E-14 * (1 + CRMB / CRMI) * KNDS(IDDS(IR,1)) *
     +             REAL(IZ)**2.0 * RIZB**2 * 15.0)

c           WRITE(0,*) TAUS

           KVALS(nks(ir)+2,1) = AMU * CRMI * KVDS (idds(ir,1)) / TAUS
           KVALS(nks(ir)+2,2) = -KFIdS(idds(ir,1)) * KBETAS(IZ) * ECH /
     +                           FACT

           KVALS(nks(ir)+2,3) = -KFEdS(idds(ir,1)) * KALPHS(IZ) * ECH /
     +                           FACT
           KVALS(nks(ir)+2,4) = IZ * ECH * KEDS (idds(ir,1))

           KVALS(nks(ir)+2,5) = KVALS(NKS(IR)+2,1) + KVALS(NKS(IR)+2,2)
     +                        + KVALS(NKS(IR)+2,3) + KVALS(NKS(IR)+2,4)
c
        endif
c
c
        DO 715 IK = 1, NKS(IR)

          KOUTS(IK+in) = KSS(IK,IR)
          KWIDS(IK+in) = 0.5 * (KBACDS(IK,IR) + KFORDS(IK,IR))

          TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +           (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +            REAL(IZ)**2.0 * RIZB**2 * 15.0)
c
c            WRITE(49,'(A,5G11.3,2I11,G11.3)') 'Data: ',
c     +       TAUS,CRMI,KTIBS(IK,IR),CRMB,KNBS(IK,IR),IZ,CIZB,
c     +       KVHS(IK,IR)/QTIM
c
c            WRITE(49,'(A,7G11.3)') '     ',
c     +       KBETAS(IZ),KFIGS(IK,IR),KTIBS(IK,IR),
c     +       KBACDS(IK,IR),KTIBS(IK-1,IR),KFORDS(IK,IR),KTIBS(IK+1,IR)
c
c            WRITE(49,'(A,2G11.3)') '     ',
c     +       KTEBS(IK,IR),(KVHS(IK,IR)/QTIM)
c
          KVALS(IK+in,1) = AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
c     >                     * KFSSMOD(ik,ir)
          KVALS(IK+in,2) = KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT

          KVALS(IK+in,3) = KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
          KVALS(IK+in,4) = IZ * KES(IK,IR) * ECH / FACT


          KVALS(IK+in,5) = KVALS(IK+IN,1) + KVALS(IK+IN,2)
     +                   + KVALS(IK+IN,3) + KVALS(IK+IN,4)
c
          WRITE(iplot,'(2I4,8G12.4)') IK,IR,KVALS(IK+IN,1),
     +      KVALS(IK+IN,2),
     +      KVALS(IK+IN,3),KVALS(IK+IN,4),KVALS(IK+IN,5),KSS(IK,IR),
     +      RS(IK,IR),ZS(IK,IR)
  715   CONTINUE

C
        enldist = mgst*kouts(nks(ir)+inc)
c
        CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >    5,99,KOUTS(1),KOUTS(NKS(IR)+inc),-HI,HI,IGNORS,ITEC,AVS,NAVS,
     >    JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)

        IF (IREF.EQ.660) THEN
          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      5,99,KOUTS(1),enldist   ,-HI,HI,IGNORS,ITEC,AVS,
     >      NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)

          CALL DRAW (KOUTS,KWIDS,KVALS,MAXNKS,NKS(IR)+inc,ANLY,
     >      5,99,KOUTS(NKS(IR)+inc)-enldist,KOUTS(NKS(IR)+inc),
     >      -HI,HI,IGNORS,ITEC,AVS,NAVS,
     >      JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,2,6,1.0,0)
        ENDIF
c
      ENDIF
c
C
C-----------------------------------------------------------------------
C     LEGRANGE POINT FOR INNER (IK = 1, 2, ...) SOL FORCES
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.662) THEN
        ELABS(1) = 'TRAPTRAP PT'
        ELABS(2) = 'X PTX PT   '

        XLAB = '  Ring     '
        YLAB = '   S       '

        FACT = QTIM**2 * EMI / CRMI

c
c Add ionisation state info here
c
        READ (GRAPH(38:44),'(I4)') IZ

        NPLOTS = NPLOTS + 1
c
        if (cgridopt.eq.0) then
           WRITE (REF,'(''Zi:'',I3,
     +             '', Outer SOL - IK=1'')') IZ
        else
           WRITE (REF,'(''Zi:'',I3,
     +             '', Inner SOL - IK=1'')') IZ
        endif
c
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL RZERO (KVALS, MAXNKS*MAXNGS)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
c
c Start of IR loop:
c
        DO IR = IRSEP, IRWALL - 1


        DO IK = 1, NKS(IR)
          IF (KSS(IK,IR).GT.0.5*KSMAXS(IR)) THEN
               MIDNKS = IK - 1
               GOTO 725
          ENDIF
        ENDDO
725     CONTINUE

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2

           TAUS = CRMI * KTIDS(IDDS(IR,2))**1.5 * SQRT(1.0/CRMB) /
     +            (6.8E-14 * (1 + CRMB / CRMI) * KNDS(IDDS(IR,2)) *
     +             REAL(IZ)**2.0 * RIZB**2 * 15.0)

c           WRITE(0,*) TAUS

c           KWIDS(1)         = kss(1,ir)
c           kouts(nks(ir)+2) = ksmaxs(ir)
c           kwids(nks(ir)+2) = ksmaxs(ir) - kss(nks(ir),ir)
c
           KVALS(1,1) = AMU * CRMI * KVDS(idds(ir,2)) / TAUS
           KVALS(1,2) = KFIdS(idds(ir,2)) * KBETAS(IZ) * ECH / FACT

           KVALS(1,3) = KFEdS(idds(ir,2)) * KALPHS(IZ) * ECH / FACT
           KVALS(1,4) = IZ * KEDS(idds(ir,2)) * ECH

           KVALS(1,5) = KVALS(1,1) + KVALS(1,2)
     +                + KVALS(1,3) + KVALS(1,4)

        ENDIF

        DO IK = 1, MIDNKS

          TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +           (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +            REAL(IZ)**2.0 * RIZB**2 * 15.0)

          KVALS(IK+in,1) = AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
c     >                     * KFSSMOD(ik,ir)
          KVALS(IK+in,2) = KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT

          KVALS(IK+in,3) = KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
          KVALS(IK+in,4) = IZ * KES(IK,IR) * ECH / FACT

          KVALS(IK+in,5) = KVALS(IK+IN,1) + KVALS(IK+IN,2)
     +                   + KVALS(IK+IN,3) + KVALS(IK+IN,4)

        ENDDO
C
        KOUTS(IR-IRSEP+1)  = REAL(IR)
        KWIDS(IR-IRSEP+1)  = 0.5

        SLVALS(IR-IRSEP+1,1) = ksmaxs(ir)/2.0

        DO IK = 1, MIDNKS
          IF (KVALS(IK+IN,5).GT.0.0) THEN

            SLVALS(IR-IRSEP+1,1) = - KVALS(IK-1+IN,5) *
     +        (KSS(IK,IR) - KSS(IK-1,IR)) /
     +         (KVALS(IK+IN,5) - KVALS(IK-1+IN,5)) + KSS(IK-1,IR)

            WRITE(49,'(I3,A,I4,5G10.3)') IR,' Zero pt: ',
     +            IK,SLVALS(IR-IRSEP+1,1),KSS(IK,IR),KSS(IK-1,IR),
     +            KVALS(IK+IN,5),KVALS(IK-1+IN,5)

            GOTO 730
          ENDIF
        ENDDO

730     CONTINUE

        SLVALS(IR-IRSEP+1,2) = 0.0

        DO IK = 2, MIDNKS
c
c         I am not sure how this will work for inverted grids.
c
          IF ((ZS(IK,IR).GT.ZXP.AND.ZS(IK-1,IR).LT.ZXP).OR.
     +        (ZS(IK,IR).LT.ZXP.AND.ZS(IK-1,IR).GT.ZXP)) THEN

            SLVALS(IR-IRSEP+1,2) = (ZXP - ZS(IK-1,IR)) *
     +        (KSS(IK,IR) - KSS(IK-1,IR)) /
     +         (ZS(IK,IR) - ZS(IK-1,IR)) + KSS(IK-1,IR)

               WRITE(49,'(A,I4,6G10.3)') '       X-pt: ',
     +           IK,SLVALS(IR-IRSEP+1,2),ZXP,KSS(IK,IR),KSS(IK-1,IR),
     +           ZS(IK,IR),ZS(IK-1,IR)

            GOTO 735
          ENDIF
        ENDDO

735     CONTINUE

        WRITE(49,'(A,I4,2G10.3,I4)') '    Results: ',
     +    IR,SLVALS(IR-IRSEP+1,1),SLVALS(IR-IRSEP+1,2),IK

c
c End of IR loop:
c
        ENDDO
c
        CALL DRAW (KOUTS,KWIDS,SLVALS,MAXNRS,IRWALL-IRSEP,ANLY,
     >    2,99,KOUTS(1),KOUTS(IRWALL-IRSEP),-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,
     >    1,6,1.0,0)

      ENDIF
c
C
C-----------------------------------------------------------------------
C     LEGRANGE POINT FOR OUTER (IK = NKS(IR),NKS(IR)-1,...) SOL FORCES
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.664) THEN
        ELABS(1) = 'TRAPTRAP PT'
        ELABS(2) = 'MID MID PL '

        XLAB = '  Ring     '
        YLAB = '  S - Smax '

        FACT = QTIM**2 * EMI / CRMI
c
c Add ionisation state info here
c
        READ (GRAPH(38:44),'(I4)') IZ

        write(0,*) 'graph: >'//graph//'<'
        write(0,*) 'iz:     ',iz

        NPLOTS = NPLOTS + 1

        if (cgridopt.eq.0) then
           WRITE (REF,'(''Zi:'',I3,
     +             '', Inner SOL - IK=NKS(IR)'')') IZ
        else
           WRITE (REF,'(''Zi:'',I3,
     +             '', Outer SOL - IK=NKS(IR)'')') IZ
        endif
c
        WRITE (IPLOT,9012) NPLOTS,REF

        CALL RZERO (KVALS, MAXNKS*MAXNGS)
c
c       Modify the following to include values for S=0 ... even if this
c       is not part of the grid. Assume that if one end is not a grid
c       point then neither is the other.
c
c
c Start of IR loop:
c
        DO IR = IRSEP, IRWALL - 1

        DO IK = 1, NKS(IR)
          IF (KSS(IK,IR).GT.0.5*KSMAXS(IR)) THEN
               MIDNKS = IK
               GOTO 737
          ENDIF
        ENDDO
737     CONTINUE

        if (kss(1,ir).eq.0.0.or.ir.lt.irsep) then
           in  = 0
           inc = 0
        else
           in       = 1
           inc      = 2

           TAUS = CRMI * KTIDS(IDDS(IR,1))**1.5 * SQRT(1.0/CRMB) /
     +            (6.8E-14 * (1 + CRMB / CRMI) * KNDS(IDDS(IR,1)) *
     +             REAL(IZ)**2.0 * RIZB**2 * 15.0)


           KVALS(nks(ir)+1,1) = AMU * CRMI * KVDS (idds(ir,1)) / TAUS
           KVALS(nks(ir)+1,2) = -KFIdS(idds(ir,1)) * KBETAS(IZ)* ECH /
     +                           FACT


           KVALS(nks(ir)+1,3) = -KFEdS(idds(ir,1)) * KALPHS(IZ) * ECH /
     +                           FACT
           KVALS(nks(ir)+1,4) = IZ * KEDS (idds(ir,1)) * ECH

           KVALS(nks(ir)+1,5) = KVALS(NKS(IR)+2,1) + KVALS(NKS(IR)+2,2)
     +                        + KVALS(NKS(IR)+2,3) + KVALS(NKS(IR)+2,4)
c
        endif

        DO IK = NKS(IR), MIDNKS, -1

          TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +           (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +            REAL(IZ)**2.0 * RIZB**2 * 15.0)

          KVALS(IK,1) = AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
c     >                  * KFSSMOD(ik,ir)
          KVALS(IK,2) = KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT

          KVALS(IK,3) = KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
          KVALS(IK,4) = IZ * KES(IK,IR) * ECH / FACT


          KVALS(IK,5) = KVALS(IK+IN,1) + KVALS(IK+IN,2)
     +                + KVALS(IK+IN,3) + KVALS(IK+IN,4)

        ENDDO
C
        KOUTS(IR-IRSEP+1)  = REAL(IR)
        KWIDS(IR-IRSEP+1)  = 0.5

        SLVALS(IR-IRSEP+1,1) = ksmaxs(ir)/2.0

        DO IK = NKS(IR), MIDNKS, -1
          IF (KVALS(IK,5).LT.0.0) THEN
            SLVALS(IR-IRSEP+1,1) = -KVALS(IK,5) *
     +        (KSS(IK+1,IR) - KSS(IK,IR)) /
     +         (KVALS(IK+1,5) - KVALS(IK,5)) + KSS(IK,IR)

            SLVALS(IR-IRSEP+1,1) = KSMAXS(IR) - SLVALS(IR-IRSEP+1,1)

               WRITE(49,'(I3,A,I4,5G10.3)') IR,' Zero pt: ',
     +           IK,SLVALS(IR-IRSEP+1,1),KSS(IK,IR),KSS(IK+1,IR),
     +           KVALS(IK+IN,5),KVALS(IK+1+IN,5)

            GOTO 740
          ENDIF
        ENDDO

740     CONTINUE

        SLVALS(IR-IRSEP+1,2) = 0.0

        DO IK = NKS(IR) - 1, MIDNKS, -1
c
c         I am not sure how this will work for inverted grids.
c
          IF ((ZS(IK,IR).GT.0.0.AND.ZS(IK+1,IR).LT.0.0).OR.
     +        (ZS(IK,IR).LT.0.0.AND.ZS(IK+1,IR).GT.0.0)) THEN

            SLVALS(IR-IRSEP+1,2) = (0.0 - ZS(IK,IR))*
     +        (KSS(IK+1,IR) - KSS(IK,IR)) /
     +         (ZS(IK+1,IR) - ZS(IK,IR)) + KSS(IK,IR)

            SLVALS(IR-IRSEP+1,2) = KSMAXS(IR) - SLVALS(IR-IRSEP+1,2)

               WRITE(49,'(A,I4,6G10.3)') '  Mid-plane: ',
     +           IK,SLVALS(IR-IRSEP+1,2),0.0,ZS(IK,IR),
     +           ZS(IK+1,IR),KSS(IK,IR),KSS(IK+1,IR)

            GOTO 745
          ENDIF
        ENDDO

745     CONTINUE

        WRITE(49,'(A,2G10.3)') '    Results: ',
     +    SLVALS(IR-IRSEP+1,1),SLVALS(IR-IRSEP+1,2)
c
c End of IR loop:
c
        ENDDO
c
        CALL DRAW (KOUTS,KWIDS,SLVALS,MAXNRS,IRWALL-IRSEP,ANLY,
     >    2,99,KOUTS(1),KOUTS(IRWALL-IRSEP),-HI,HI,IGNORS,ITEC,AVS,
     >    NAVS,JOB,TITLE,XLAB,YLAB,ELABS,REF,NVIEW,PLANE,TABLE,
     >    1,6,1.0,0)

      ENDIF
c
c slmod end


c
c kkmod
c
C
C-----------------------------------------------------------------------
C     CONTOURS OF NET FORCE ON IMPURITY CHARGE STATE
C-----------------------------------------------------------------------
C
c      IF (IREF.EQ.109.OR.IREF.EQ.110) THEN
c
      IF (IREF.EQ.669.OR.IREF.EQ.670) THEN
        XLAB   = '   R  (M)'
        YLAB   = '   Z  (M)'
c
        IZ     = IOPT
c
        REF    = 'IMPURITY NET FORCE ' // XPOINT
        CALL RZERO (KTMP, MAXNRS*MAXNKS)
        WRITE (IPLOT,9012) NPLOTS,REF
        FACT = QTIM**2 * EMI / CRMI
        DO 2601 IR = 1,NRS
          DO 2601 IK = 1,NKS(IR)
            TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +             (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +             REAL(IZ)**2.0 * RIZB**2 * 15.0)
            TMPSUM =          AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
            TMPSUM = TMPSUM + KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT

            TMPSUM = TMPSUM + KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
            TMPSUM = TMPSUM + IZ * KES(IK,IR) * ECH / FACT
            KTMP(IK,IR) = TMPSUM
 2601   CONTINUE
c       normalization added by Krieger IPP/97
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)

c       if filled plot selected, use special color table for
c       directed quantities  Krieger, IPP/97
        if (icntr.eq.0) then
          CALL CONTOUR (ICNTR,NGS,KTMP,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,cntropt,minscale,maxscale)
        else
c
c         Setup colours
c
          call setup_col(41,4) 
c
          CALL CONTOUR (ICNTR,NGS,KTMP,1,1,1,FT,FP,1.0,
     >                  XOUTS,1,NXS,YOUTS,1,NYS,
     >                  XXMIN,XXMAX,YYMIN,YYMAX,
     >                  nconts,conts,5,minscale,maxscale)
c
c         restore colours
c
          call setup_col(n_cols,col_opt) 
c
        endif
      ENDIF
c
c kkmod
c
c
C-----------------------------------------------------------------------
c
c     Plots from 671 to 699 are for EDGE2D/Fluid code quantities.
c
C-----------------------------------------------------------------------
C     FLUID CODE - BG NEUTRAL ATOM DENSITY CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.671.OR.IREF.EQ.672) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'FLUID CODE-BG ATOM DEN '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,E2DATOM,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
C
C-----------------------------------------------------------------------
C     FLUID CODE - NEUTRAL IMPURITY DENSITY CONTOURS
C-----------------------------------------------------------------------
C
      IF (IREF.EQ.673.OR.IREF.EQ.674) THEN
        XLAB   = '   X/A  (M)'
        YLAB   = '   Y/A  (M)'
        NGS    = CNGS
        IZ     = IOPT
        REF = 'FLUID CODE-NEUTRAL IMP.DEN '// XPOINT
        WRITE (IPLOT,9012) NPLOTS,REF
        CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     >    YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NGS)
        CALL CONTOUR (ICNTR,NGS,E2DZ0,1,1,1,FT,FP,1.0,
     >                XOUTS,1,NXS,YOUTS,1,NYS,
     >                XXMIN,XXMAX,YYMIN,YYMAX,
     >                nconts,conts,cntropt,minscale,maxscale)
      ENDIF
c
      IF (iref.EQ.699) THEN
c...    Thomson shift analysis plot:
        call setup_col(n_cols,5)
        CALL PlotXXX(cngs,job,graph,nplots,ref,title,iopt,iplot,
     .               xxmin,xxmax,yymin,yymax,icntr,ft,fp,
     .               xouts,youts,nconts,conts,cntropt)
        call setup_col(n_cols,col_opt)
      endif
 
      return
c
c     Format statements    
c
 9012 FORMAT(1X,'PLOT',I3,4X,A)

      end

