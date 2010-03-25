c     -*-Fortran-*-


      SUBROUTINE LoadCameraData(osmtmp,fname,scale)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'slcom'

c...OSMTMP not looking dynamic here...
      REAL      osmtmp(MAXNKS,MAXNRS),scale
      CHARACTER fname*(*)

      INTEGER i1,i2,ik,ir,nc,npt,fp1,ic
      REAL    xcen,ycen,dpt(100),vpt(100),dist,totdpt,tmpdpt,tmpvpt
      CHARACTER dummy*1024

      REAL,    ALLOCATABLE :: cq(:),rc(:),zc(:)

      INTEGER    MAXCELL
      PARAMETER (MAXCELL=MAXNKS*MAXNRS+MAXASC)

      ALLOCATE(cq(MAXCELL))
      ALLOCATE(rc(MAXCELL))
      ALLOCATE(zc(MAXCELL))
      nc = 0
      cq = 0.0

      IF (scale.EQ.0.0) 
     .  CALL ER('LoadCameraData','SCALE not defined',*99)

      fp1 = 99
      OPEN(fp1,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     STATUS='OLD')

c...  Read header:
      DO WHILE (.TRUE.)
        ! jdemod - * by itself in the end position isn't meaningful 
        !          assume it is a * format specifier and move it to 2nd argument
        !READ(fp1,ERR=98,END=98,*) dummy
        READ(fp1,*,ERR=98,END=98) dummy
        IF (dummy(1:1).NE.'#'.AND.dummy(1:1).NE.'*') THEN
          BACKSPACE fp1
          EXIT
        ENDIF
      ENDDO

c...  Load data:
      nc = 0        
      DO WHILE(.TRUE.) 
        ! jdemod - * by itself in the end position isn't meaningful 
        !          assume it is a * format specifier and move it to 2nd argument
        !READ(fp1,ERR=98,END=9,*) rc(nc+1),zc(nc+1),cq(nc+1)
        READ(fp1,*,ERR=98,END=9) rc(nc+1),zc(nc+1),cq(nc+1)
        nc = nc + 1
        IF (nc.GT.MAXCELL) 
     .    CALL ER('LoadCameraData','MAXCELL exceeded, increase '//
     .            'value',*99)
      ENDDO
 9    CLOSE(fp1)

c...  Process data:
      DO i1 = 1, nc
        cq(i1) = MAX(0.0,cq(i1)) 
      ENDDO

c...  Map data to grid:
      DO ir = 2, nrs
        IF (idring(ir).EQ.BOUNDARY) CYCLE

        DO ik = 1, nks(ir)
          npt = 0
          vpt = 0.0
          dpt = 0.0
          xcen = rs(ik,ir)
          ycen = zs(ik,ir)
          DO i1 = 1, nc
c *BUG* (SORT OF)
c            IF (cq(i1).LE.0.0) CYCLE
            dist = SQRT((xcen - rc(i1))**2 + (ycen - zc(i1))**2)
            IF (dist.LT.scale) THEN
c...          Store up to 100 data points:
              IF (npt.LT.100) npt = npt + 1
              IF (dpt(npt).EQ.0.0.OR.dist.LT.dpt(npt)) THEN
                dpt(npt) = MAX(dist,1.0E-12) 
                vpt(npt) = MAX(0.0,cq(i1))
c...            Sort data points from closest to farthest:                 
                i2 = 1
                DO WHILE(i2.LT.npt)
                  i2 = i2 + 1
                  IF (dpt(i2-1).GT.dpt(i2)) THEN
                    tmpdpt = dpt(i2-1)
                    tmpvpt = vpt(i2-1)
                    dpt(i2-1) = dpt(i2)
                    vpt(i2-1) = vpt(i2)
                    dpt(i2) = tmpdpt
                    vpt(i2) = tmpvpt
                    i2 = 1   
                  ENDIF
                ENDDO
              ENDIF
            ENDIF
          ENDDO

c...      Calculate weight function for averaging:
          totdpt = 0.0
          DO i1 = 1, npt
c *TEMP*
c            dpt(i1) = 1.0 / dpt(i1)
c *BUG* (SORT OF)
            dpt(i1) = 1.0 / (dpt(i1)**2)
            totdpt = totdpt + dpt(i1)
          ENDDO

c...      Average the data points:
          IF (totdpt.GT.0.0) THEN
            DO i1 = 1, npt
              osmtmp(ik,ir) = osmtmp(ik,ir) + dpt(i1) / totdpt * vpt(i1)
c DEBUG:
c *TEMP*
              IF (ik.EQ.134.AND.ir.EQ.20) THEN
C              IF (ik.GT.128..AND.ik.LT.139.AND.ir.EQ.20) THEN
                WRITE(6,'(A,2I6,1P,5E10.2,0P)') 
     .            'VALUE:',ik,ir,dpt(i1)/totdpt,dpt(i1),vpt(i1),
     .            osmtmp(ik,ir)
              ENDIF

            ENDDO
          ENDIF

        ENDDO
      ENDDO

      DEALLOCATE(cq)
      DEALLOCATE(rc)
      DEALLOCATE(zc)

      RETURN
 98   CALL ER('LoadCameraData','Problems accessing data file',*99)
      WRITE(0,*) '>'//fname(1:LEN_TRIM(fname))//'<'
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SetCol255
c
c
      SUBROUTINE SetCol255(frac1,qmin,qmax)
      IMPLICIT none
      REAL frac1,frac,qmin,qmax,hue,bright,pastel,midfrac,scale

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale

      LOGICAL grayscale

      grayscale = .FALSE.

      frac = frac1
      IF (qmin.LT.-1.0E-6.AND.qmax.GT.1.0E-06) THEN
c...
        midfrac = -qmin / (qmax - qmin) 
        scale = MAX(midfrac,1.0-midfrac)

        IF (frac.LT.midfrac) THEN
          hue    = 0.34
          pastel = (midfrac - frac) / scale
          bright = 1.0
        ELSE
          hue    = 0.0
          pastel = (frac - midfrac) / scale
          bright = 1.0
        ENDIF
        CALL ColSet(hue,pastel,bright,255)

      ELSEIF (grayscale) THEN
        IF (hardscale.AND.frac.LT.0.07) THEN
          CALL ColSet(0.0,0.0,1.0,255)
        ELSE
          IF (frac.GT.1.0) STOP 'sdfsdsgsd'
          IF (frac.LT.0.0) STOP 'sdfsdsgsd adsfasd'
          frac = frac * 0.9 + 0.1
          CALL ColSet(0.0,0.0,1.0-frac,255)
        ENDIF

      ELSEIF (hardscale) THEN

        IF (.TRUE.) THEN
c        IF (.NOT..TRUE.) THEN
          IF (.FALSE..AND.frac.LT.0.01) THEN
            bright = 1.0
            frac   = 1.0
            pastel = 0.0
          ELSE
            frac = frac * 0.93 + 0.07
            IF (frac.LE.0.27) THEN
              bright = 1.0-((0.27-frac)/(0.27-0.07))**2
              bright = MAX(0.0,bright)
            ELSE
              bright = 1.0
            ENDIF
            frac = (1.0 - frac) * 0.90
            frac = frac + 0.34
            IF (frac.GT.1.0) frac = frac - 1.0
            pastel = 1.0
          ENDIF
        ELSE
          IF (frac.LT.0.07) THEN
            bright = 1.0
            frac   = 1.0
            pastel = 0.0
          ELSE
            IF (frac.LE.0.27) THEN
              bright = 1.0-((0.27-frac)/(0.27-0.07))**2
              bright = MAX(0.0,bright)
            ELSE
              bright = 1.0
            ENDIF
            frac = (1.0 - frac) * 0.90
            frac = frac + 0.34
            IF (frac.GT.1.0) frac = frac - 1.0
            pastel = 1.0
          ENDIF
        ENDIF
        CALL ColSet(frac,pastel,bright,255)

      ELSE
        bright = 1.0-(0.98-frac)**20
        CALL ColSet(frac,1.0,bright,255)

      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: PlotContour
c
c
      SUBROUTINE PlotContour
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'     

      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

      REAL       TOL
      PARAMETER (TOL=1.0E-5)

      INTEGER i1,i2,nver(2),imin,icut,ngrp,imid
      LOGICAL farenough
      REAL    xtmp(2048),ytmp(2048),xver(2048,2),yver(2048,2),xmin,xmid,
     .        xmax

c...  Check if end points of the contour make sense:
      IF (SQRT((xcon(1)-xcon(icon))**2+
     .         (ycon(1)-ycon(icon))**2).GT.0.01) THEN

        IF     (ABS(xcon(1)-xcon(icon)).LT.TOL) THEN
        ELSEIF (ABS(ycon(1)-ycon(icon)).LT.TOL) THEN
        ELSEIF (xcon(1).GT.xcon(icon).AND.
     .          ycon(1).LT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(1)
          ycon(icon) = ycon(icon-1)
        ELSEIF (xcon(1).LT.xcon(icon).AND.
     .          ycon(1).LT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(icon-1)
          ycon(icon) = ycon(1)
c        ELSEIF (xcon(1).LT.xcon(icon).AND.
c     .          ycon(1).LT.ycon(icon)) THEN
c          icon = icon + 1
c          xcon(icon) = xcon(icon-1)
c          ycon(icon) = ycon(1)
        ELSEIF (xcon(1).LT.xcon(icon).AND.
     .          ycon(1).GT.ycon(icon)) THEN
          icon = icon + 1
          xcon(icon) = xcon(icon-1)
          ycon(icon) = ycon(1)
        ELSE
c          WRITE(0,*) 'xCON=',
c     .      xcon(1),ycon(1),xcon(icon),ycon(icon)
          CALL ER('PlotContour','Unknown situation',*99)
        ENDIF

      ENDIF
 

c      IF (.FALSE..AND.icon.GT.999) THEN
c
      IF (icon.GT.2000) THEN

        WRITE(0,*) 'ICON= ',icon
        STOP 'CONTOUR ARRAY TOO LARGE, INCREASE ARRAY SIZES IN GHOST'
c...    Split array of verticies

        ngrp = 2

c...    Find the mid-point of the contour:
        xmin =  HI
        xmax = -HI
        DO i1 = 1, icon
          IF (xcon(i1).LT.xmin) THEN
            imin = i1
            xmin = xcon(i1)
          ENDIF
          xmax = MAX(xmax,xcon(i1))
        ENDDO
c...    Copy contour array starting at imin:
        i1 = imin
        DO i2 = 1, icon
          xtmp(i2) = xcon(i1)
          ytmp(i2) = ycon(i1)
          i1 = i1 + 1
          IF (i1.GT.icon) i1 = 1
        ENDDO

c...    Find the x-axis midpoint and cut the array there:
        xmid = 0.5 * (xmin + xmax)
        DO i1 = 1,icon
          IF (xtmp(i1).GT.xmid) THEN
            imid = i1
            EXIT
          ENDIF
        ENDDO
        farenough = .FALSE.
        icut = 0
        DO i1 = imid+1,icon
          IF (xtmp(i1).GT.1.1*xmid) farenough = .TRUE.
          IF (farenough.AND.xtmp(i1).LT.xmid) THEN
            icut = i1
            EXIT
          ENDIF
        ENDDO
        IF (imid.EQ.icon) STOP 'SHIT'

        ngrp = 2

        nver(1) = 0
        DO i1 = 1, imid
          nver(1) = nver(1) + 1
          xver(nver(1),1) = xtmp(i1)
          yver(nver(1),1) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      1,nver(1),xver(nver(1),1),yver(nver(1),1)
        ENDDO
        DO i1 = icut, icon
          nver(1) = nver(1) + 1
          xver(nver(1),1) = xtmp(i1)
          yver(nver(1),1) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      1,nver(1),xver(nver(1),1),yver(nver(1),1)
        ENDDO

        nver(2) = 0
        DO i1 = imid, icut
          nver(2) = nver(2) + 1
          xver(nver(2),2) = xtmp(i1)
          yver(nver(2),2) = ytmp(i1)
          WRITE(6,*) 'CNTR:',
     .      2,nver(2),xver(nver(2),2),yver(nver(2),2)
        ENDDO

c        WRITE(0,*) '--->',nver(1),nver(2),icon

      ELSE
c...    Copy array of verticies as-is:
        ngrp = 1
        nver(1) = icon
        DO i1 = 1, icon
          WRITE(6,*) 'CNTR:',i1,xcon(i1),ycon(i1)
          xver(i1,1) = xcon(i1)
          yver(i1,1) = ycon(i1)
        ENDDO
      ENDIF

      DO i1 = 1, ngrp
        IF (nver(i1).GT.2048) WRITE(0,*) 'WARNING: NVER TOO BIG'
        CALL PTPLOT(xver(1,i1),yver(1,i1),1,nver(i1),1)
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: FilterContour
c
c
      SUBROUTINE FilterContour
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'     

      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

      INTEGER i1,i2
      LOGICAL status
      REAL    t1,t2

      RETURN

c...  Scan through the contour verticices and see if any are along the same trajectory, 
c     and delete those:      

      status = .TRUE.
      DO WHILE(status)
        status = .FALSE.
        DO i1 = 1, icon-2
          t1 = (xcon(i1+1) - xcon(i1)) /
     .         (xcon(i1+2) - xcon(i1))
          t2 = (ycon(i1+1) - ycon(i1)) /
     .         (ycon(i1+2) - ycon(i1))
          IF (xcon(i1).EQ.xcon(i2)) THEN
            WRITE(0,*) xcon(i1),xcon(i1+1),xcon(i1+2)
            WRITE(0,*) ycon(i1),ycon(i1+1),ycon(i1+2)
            WRITE(0,*) t1,t2,i1,icon
          ENDIF
          IF (ABS(t1-t2).LT.1.0E-05) THEN
c...        Delete vertex:
            DO i2 = i1+1, icon-1
              xcon(i2) = xcon(i2+1)
              ycon(i2) = ycon(i2+1)
            ENDDO
            icon = icon - 1
            status = .TRUE.
          ENDIF
        ENDDO
      ENDDO

      WRITE(0,*) 'REDUCED ICON=',icon

      RETURN
99    STOP
      END
c
c
c
      SUBROUTINE CONTIL2(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,NPTSY,
     &                   CLEVLS,ISTRTL,ISTOPL,XGRIDS,YGRIDS)
      implicit none
C
C          ------------------------------------------------
C          ROUTINE NO. ( 124)   VERSION (A8.8)    15:MAY:87
C          ------------------------------------------------
C
C          THIS DRAWS STRAIGHT-ELEMENT CONTOURS ON AN IRREGULAR GRID.
C
C
C          THE ARGUMENTS ARE AS FOLLOWS:
C
C          [SURFAS]  IS THE ARRAY OF SURFACE HEIGHT VALUES,
C          <ISTRTX>  IS THE LOWER X-EXTENT,
C          <ISTOPX>  IS THE UPPER X-EXTENT, WHILE
C          <ISTRTY>  AND
C          <ISTOPY>  ARE THE CORRESPONDIMG Y-BOUNDS.
C          <NPTSX>   IS THE ACTUAL ARRAY X-EXTENT, AND
C          <NPTSY>   IS THE ACTUAL ARRAY Y-EXTENT.
C          [CLEVLS]  CONTAINS THE VALUES OF CONTOUR HEIGHTS,
C          <ISTRTL>  IS THE STARTING POINT, AND
C          <ISTOPL>  IS THE END POINT OF THIS ARRAY.
C          [XGRIDS]  ARE THE GRID X-POSITIONS.
C          [YGRIDS]  ARE THE GRID Y-POSITIONS.
C
C
      integer nptsx,nptsy,istopl,llbcon,level,istptx,istpty,
     >        iprint,itrac1,itrac2,itrac3,itrac4,numerr,nbitsw,
     >        istrtx,istopx,istrty,istopy,istrtl,ilenx,ileny,
     >        iprsav,istrx1,istry1,istpx1,istpy1,levelh,ibit,
     >        ixwrd,ixbit,negix,negiy
      real angcon,stang0,crang0,height,mapbit,cwidx,cwidy,xplot0,
     >     yplot0,x1wnd0,x2wnd0,y1wnd0,y2wnd0,angsav,unsave,xhere,
     >     yhere

      REAL    SURFAS(NPTSX,NPTSY),CLEVLS(ISTOPL),
     &        XGRIDS(NPTSX),YGRIDS(NPTSY)
      LOGICAL OPENCO,ERRON,CURVED
C
      COMMON /T0ACON/ ANGCON
      COMMON /T0CANG/ STANG0,CRANG0
      COMMON /T0CLBL/ LLBCON
      COMMON /T0CLEV/ LEVEL,HEIGHT
      COMMON /T0CMAP/ MAPBIT(6,2048),ISTPTX,ISTPTY
      COMMON /T0CTYP/ OPENCO,CURVED
      COMMON /T0CWID/ CWIDX,CWIDY
      COMMON /T0PPOS/ XPLOT0,YPLOT0
      COMMON /T0TRAC/ IPRINT
      COMMON /T0TRAI/ ITRAC1,ITRAC2,ITRAC3,ITRAC4
      COMMON /T0WNDO/ X1WND0,X2WND0,Y1WND0,Y2WND0
      COMMON /T3ERRS/ ERRON,NUMERR
c slmod begin
      COMMON /CONCON/ icon,iconcol,xcon,ycon
      INTEGER icon,iconcol
      REAL    xcon(2048),ycon(2048)

      include 'colours'
c slmod end
C
      DATA NBITSW /32/
C
C
      ITRAC1= ISTRTX
      ITRAC2= ISTOPX
      ITRAC3= ISTRTY
      ITRAC4= ISTOPY
      IF (IPRINT.EQ.1) CALL G0MESG(48,8)
C
      ILENX= ISTOPX-ISTRTX
      ILENY= ISTOPY-ISTRTY
      IF (ISTRTX.LT.1.OR.ISTRTY.LT.1)         GO TO 901
      IF (ILENX.LT.1.OR.ILENX.GT.2048)         GO TO 901
      IF (ILENY.LT.1.OR.ILENY.GT.2048)         GO TO 901
      IF (ISTOPX.GT.NPTSX.OR.ISTOPY.GT.NPTSY) GO TO 901
      IF (ISTRTL.LT.1.OR.ISTRTL.GT.ISTOPL)    RETURN
C
C          THE ROUTE-TRACING FLAG, THE CURVE METHOD AND THE
C          CHARACTER MAGNIFICATION ARE SAVED, THEN TRACE IS
C          DISABLED, CURVE METHOD (1) SET AND MAGN. (6) SET.
C
      IPRSAV= IPRINT
      IPRINT= 0
      CURVED= .FALSE.
      ANGSAV= STANG0
      UNSAVE= ANGCON
      ANGCON= 1.0
      IF (LLBCON.GT.0) CALL CTRORI(0.5)
C
      CALL G0AUTO(XGRIDS,YGRIDS,ISTRTX,ISTOPX,ISTRTY,ISTOPY,1)
      XHERE= XPLOT0
      YHERE= YPLOT0
      CWIDX= (X2WND0-X1WND0)/(ISTOPX-ISTRTX)
      CWIDY= (Y2WND0-Y1WND0)/(ISTOPY-ISTRTY)
C
C          TO BEGIN WITH, SOME CONSTANTS ARE SET UP.
C
      ISTRX1= ISTRTX+1
      ISTRY1= ISTRTY+1
      ISTPX1= ISTOPX-1
      ISTPY1= ISTOPY-1
C
C          LOOP-100 DOES ALL THE CONTOURS AT EACH HEIGHT IN TURN.
C          THE BIT ARRAY <MAPBIT> IS INITIALISED FOR EACH HEIGHT,
C          ELEMENTS BEING SET TO '1' WHEN A CONTOUR LINE OF THE
C          CURRENT HEIGHT CROSSES AN EDGE IN AN UPWARDS DIRECTION
C          PROVIDED 'HIGH GROUND' LIES TO THE RIGHT OF THE LINE.
C

c...dev
      DO 100 LEVELH= ISTRTL,ISTOPL
        HEIGHT= CLEVLS(LEVELH)

c...dev
        WRITE(6,*) 'HEIGHT=',height,ISTOPL,ISTRTL,iconcol

        LEVEL= LEVELH
C
        DO 200 ISTPTY= ISTRY1,ISTPY1
          DO 200 ISTPTX= ISTRX1,ISTOPX
            IBIT= 0
            IF (SURFAS(ISTPTX-1,ISTPTY).LT.HEIGHT.AND.
     &          SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)     IBIT= 1
C
            IXWRD= (ISTPTX-ISTRX1)/NBITSW+1
            IXBIT= MOD((ISTPTX-ISTRX1),NBITSW)+1
            CALL G4PUTB(MAPBIT(IXWRD,(ISTPTY-ISTRTY)),IXBIT,IBIT)
  200   CONTINUE
C
C          <OPENCO> IS INITIALISED AND A SEARCH FOR OPEN
C          CONTOURS BEGINNING ON EACH EDGE IS MADE IN TURN:
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTY>= <ISTRTY>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        OPENCO= .TRUE.
        ISTPTY=  ISTRTY
C
        DO 300 ISTPTX= ISTRX1,ISTOPX

       icon = 0

          IF (SURFAS(ISTPTX-1,ISTPTY).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,-1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'A'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  300   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTX>= <ISTOPX>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTX= ISTOPX
C
        DO 400 ISTPTY= ISTRY1,ISTOPY

       icon = 0

          IF (SURFAS(ISTPTX,ISTPTY-1).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,0,-1)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'B'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  400   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTY>= <ISTOPY>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTY= ISTOPY
C
        DO 500 NEGIX= ISTRTX,ISTPX1
          ISTPTX= ISTPX1+ISTRTX-NEGIX

       icon = 0

          IF (SURFAS(ISTPTX+1,ISTPTY).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour            
          ENDIF
          WRITE(6,*) 'C'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  500   CONTINUE
C
C          CONTOURS BEGINNING ON THE EDGE <ISTPTX>= <ISTRTX>
C          ARE FOUND, THEN FOLLOWED AND DRAWN USING <G0CTIA>.
C
        ISTPTX= ISTRTX
C
        DO 600 NEGIY= ISTRTY,ISTPY1
          ISTPTY= ISTPY1+ISTRTY-NEGIY

       icon = 0

          IF (SURFAS(ISTPTX,ISTPTY+1).LT.HEIGHT.AND.
     &        SURFAS(ISTPTX,  ISTPTY).GE.HEIGHT)
     &        CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,ISTOPY,
     &                    NPTSY,XGRIDS,YGRIDS,0,1)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'D'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  600   CONTINUE
C
C          A SEARCH IS MADE FOR CLOSED CONTOURS. <OPENCO> IS
C          SET TO .FALSE. AND <MAPBIT> IS SCANNED: IF THE
C          ELEMENT (ISTPTX,ISTPTY) IS '1', A CLOSED CONTOUR
C          STARTS FROM THE LINE JOINING (ISTPTX-1,ISTPTY) TO
C          (ISTPTX,ISTPTY) AND IT IS FOLLOWED USING <G0CTIA>.
C
        OPENCO= .FALSE.
C
c slmod begin
        IF (.FALSE..AND.icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048) THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
          ENDIF
          WRITE(6,*) 'D.5'
          CALL PlotContour
          icon = 0
        ENDIF
c slmod begin
        DO 700 NEGIY= ISTRY1,ISTPY1
          ISTPTY= ISTOPY+ISTRTY-NEGIY
C
          DO 700 NEGIX= ISTRTX,ISTPX1
            ISTPTX= ISTOPX+ISTRTX-NEGIX
            IXWRD= (ISTPTX-ISTRX1)/NBITSW+1
            IXBIT= MOD((ISTPTX-ISTRX1),NBITSW)+1
            CALL G4GETB(MAPBIT(IXWRD,(ISTPTY-ISTRTY)),IXBIT,IBIT)

       icon = 0

            IF (IBIT.EQ.1)
     &              CALL G0CTIA(SURFAS,ISTRTX,ISTOPX,NPTSX,ISTRTY,
     &                          ISTOPY,NPTSY,XGRIDS,YGRIDS,-1,0)

       IF (icon.GT.0) THEN
          CALL FILCOL(255)
          CALL LINCOL(0)
          IF (icon.GT.2048)THEN
            WRITE(0,*) 'WARNING: MAX ICON EXCEEDED',ICON
c            CALL FilterContour
          ENDIF
          WRITE(6,*) 'E'
          CALL PlotContour
c          CALL PTPLOT(xcon,ycon,1,MIN(icon,2048),1)
       ENDIF

  700   CONTINUE
  100 CONTINUE
C
C          THE ENTRY STATE IS ALWAYS RESTORED BEFORE ENDING.
C
      CALL POSITN(XHERE,YHERE)
      CALL CTRORI(ANGSAV)
      ANGCON= UNSAVE
      IPRINT= IPRSAV
      RETURN
C
 901  NUMERR= 11
      IF (ERRON) CALL G0ERMS
C
      RETURN
      END


c
c ======================================================================
c
c subroutine: Plot982
c
c 2D neutral pressure plot
c
c
      SUBROUTINE Plot982(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs,zval)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'dynam2'
      INCLUDE 'cedge2d'
      INCLUDE 'pindata'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      COMMON /DUMCOM/ map1x_d,map2x_d,map1y_d,map2y_d,tag_d
      INTEGER         tag_d
      REAL            map1x_d,map2x_d,map1y_d,map2y_d

      REAL CalcPressure,GetCs,GetEAD

      integer cngs,cntropt

      integer iplot
      REAL    NP,NS,NT,FP,FT,TIME,TIME1,ZA02AS,ASPECT,FACT,POINT
      INTEGER IOPT,J,NIZS,ITER,NITERS,IZ,JZ,ISMOTH,JK,IN,inc
      integer  nplts,ringnos(maxplts),ip,pnks(maxplts)
      CHARACTER TITLE*(*),TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256,
     .          graph6*128

      real mvals(maxnks,maxplts,maxngs)
      real mouts (maxnks,maxplts),mwids (maxnks,maxplts)
      character*36 pltlabs(maxplts)

      CHARACTER*36 XLAB,XPOINT
      CHARACTER*72 YLAB
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 NAME,PLABS(-2:MAXPLRP),KLAB

      CHARACTER*128 elabs(MAXNGS)

      INTEGER IX,IY,IG,NGS,IERR,IXMIN,IXMAX,M,ID,JR

      INTEGER ii1,ii2,midnks,i1,i2,ret,plt_opt1,osm_cfs,npts,ntime,
     .        in1,in2,xtype,ytype,btype,id1,id2
      integer  sctype,ngrm
      real pltmins(maxplts),pltmaxs(maxplts)
      real pltfact,i1r

      REAL tauii,pii
      INTEGER nenum,tenum,opt_const,plot_mode(30),iter1,iter2,
     .        xaxis,ring,mode,inorm(MAXNGS)

      REAL          te1,ti1,ne1,te2,ti2,ne2,norm,
     .              rdum1,rdum2,rdum3,rdum4,rdum5,rdum6,rdum7,
     .              radum1(MAXNRS),radum2(MAXNRS),radum3(MAXNRS)
      REAL    nemin,nestep,temin,temax,neTe,frac1,xrange1,xrange2,
     .        max1,max2,ynorm(MAXNGS),frac2,bright,pastel,miny,maxy

      REAL cdata(MAXNKS,MAXNRS)
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,SUM(10),VAL,VALTOT
      CHARACTER*72 SMOOTH
      integer nconts,nclev
c      real conts(maxpts),clev(maxpts)

      REAL XOUTS(MAXGXS),XWIDS(MAXGXS),XVALS(MAXGXS,MAXNGS)
      REAL YOUTS(MAXGYS),YWIDS(MAXGYS),YVALS(MAXGYS,MAXNGS)


      CHARACTER*128 dataline,cdum1,cdum2

c...630:
      INTEGER i,k
      REAL    PLTMAX,PLTMIN
      REAL LOUTS(MAXSEG),LWIDS(MAXSEG),LVALS(MAXSEG,MAXNGS)
      REAL ydata(MAXSEG)
      integer llabs(maxseg)

c...980:
      INTEGER NUMTHE,AVPTS,ATYPE
      INTEGER numth2,numthe1(MAXNGS)
      INTEGER IGNORS(MAXNGS),ITEC,NAVS
      REAL TOUTS(MAXTHE),TWIDS(MAXTHE),TVALS(MAXTHE,MAXNGS),
     .     touts1(MAXTHE,MAXNGS)
      REAL TOUTS2(MAXTHE),TWIDS2(MAXTHE),TVALS2(MAXTHE,MAXNGS),
     .     dum1(MAXTHE),dum2(MAXTHE),dum3(MAXTHE),dum4(MAXTHE),
     .     dum5(MAXTHE),dum6(MAXTHE),dum7(MAXTHE),dum8(MAXTHE),
     .     den1,teb1
      REAL ZOBS,ROBS,DRAD,DTHE,THEMIN,THEMAX,themin_start
      real theres,mfact
      REAL WLNGTH
      real    zsuma(maxizs)
      REAL LEVEL,AVS(0:100),VMIN,VMAX
      CHARACTER ADASID*80,PECTITLE*120,PLABAD*36
      CHARACTER XFESYM*2
      character adasex*3
      integer   adasyr
      CHARACTER ADASGR*8,ADASTY*80
C
C     Second sets of ADAS data for RATIO plots
C
      character graph5*80,adasid2*80,adasex2*3
      integer   adasyr2,isele2,iselr2,iselx2,iseld2
      integer   iz_state,z_atom,iz_state2,z_atom2
      INTEGER plot
      character graph2*80,graph3*80,graph4*80
      INTEGER ISELE,ISELR,ISELX,iseld,iseldef
      INTEGER IADAS,NPAIRS,IRCODE,IKK,LEN,LENSTR
      INTEGER IK,II,IT,LT,UT,IREF,IYMIN,IYMAX,IR,JD,cnt
      INTEGER IZMIN,IZMAX,IW,LW,UW
      REAL PLRPAD(MAXNKS,MAXNRS)
      LOGICAL plotr
      real zadj
      REAL peak1,peak2,array(MAXNKS,MAXNRS)
      INTEGER plotcode,line,iopt1,iopt2,iopt3
      CHARACTER cname*256,resdir*256,cmnd*256
      LOGICAL status,oldraw


      COMMON /NSMOOTH/ NUMSMOOTH,cgrprint
      INTEGER NUMSMOOTH,cgrprint

c...  For reading from the .experiment file (UNIT=13):
      INTEGER    MAXTDAT     ,MAXCOLS   
      PARAMETER (MAXTDAT=1000,MAXCOLS=10)
      INTEGER   etype,ndata,ncol
      REAL      xdata(MAXTDAT),edata(MAXTDAT,MAXCOLS)
      CHARACTER datatitle*128


c...  982:
      INTEGER    MAXCELL                     ,MAXPOINT
      PARAMETER (MAXCELL=MAXNKS*MAXNRS+MAXASC,MAXPOINT=100)

      REAL, ALLOCATABLE :: osmtmp(:,:)

      INTEGER, ALLOCATABLE :: nv(:),rq(:)
      REAL, ALLOCATABLE :: rv(:,:),zv(:,:),cq(:),rc(:),zc(:)

      INTEGER nc,shade,shift,set,avoid(MAXASCDAT),idum1,fp1,
     .        istart,iend,i3,cut1,optflow,size,irs,ire,izstate
      LOGICAL vacuum,standard,inside,setqmin,setqmax
      REAL    x1,y1,qmin,qmax,frac,spot,dspot,range1,range2,zval,xval,
     .        xmin,ymin,xmax,ymax,red,green,blue,frac5,fmod5,zlim,
     .        rlim,stepsize,scale,deltar,deltaz,
     .        pmax(2,MAXNRS),pval,xpos,ypos,p1,p2,pavg,vtot,Ti,nD,
     .        ne,sigmav,sigmav1,sigmav2,marsum1,marsum2,cellvol,nd2p
      CHARACTER label*128,dummy*5000,caption*5000

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale


      REAL, ALLOCATABLE :: tmppinion(:,:)
      REAL, ALLOCATABLE :: tmposmcfp(:,:)


c...  Contour plot option:
      LOGICAL centerinside
      REAL*8 t1,t2
c
c     jdemod - added the allocatable attribute
c
      real,allocatable:: image(:,:),image1(:,:),raxis(:),zaxis(:)  
      integer,allocatable :: tpt(:,:),tpr(:,:)
c
      real uconts(maxpts)
      integer icntr,ncntr,xres,yres,maxix,maxiy,nix,niy,ic,
     .        icount,npt
      real minscale,maxscale,dpt(100),vpt(100),tmpdpt,tmpvpt,dist,totdpt
      real xcen,ycen,xnear,ynear
      character*36 blabs(2)
      character*1024 fname

c      CALL THICK2(4)

      iopt_ghost = 0

      call setup_col(10,5)

c      WRITE(0,*) '982:'

      CALL IZero(avoid,MAXASC)

      vacuum   = .FALSE.
      standard = .FALSE.

      marsum1 = 0.0
      marsum2 = 0.0

      ALLOCATE(osmtmp(MAXNKS,MAXNRS))

      ALLOCATE(nv(MAXCELL))
      ALLOCATE(rv(MAXPOINT,MAXCELL))
      ALLOCATE(zv(MAXPOINT,MAXCELL))
      ALLOCATE(rq(MAXCELL))
      ALLOCATE(cq(MAXCELL))
      ALLOCATE(rc(MAXCELL))
      ALLOCATE(zc(MAXCELL))
      nc = 0
      cq = 0.0

      IF     (iopt.EQ. 1) THEN
c...    D2 density plot for standard and vacuum grid:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 2) THEN
c...    D2 density plot for vacuum grid:
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 3) THEN
c...    D2 density plot for standard grid:
        standard = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 4) THEN
c...    D density plot for standard and vacuum grid:
        WRITE(0,*) 'WARNING 982: 4 - USING PINATOM!'
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 5) THEN
c...    D density plot for vacuum grid:
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 6) THEN
c...    D density plot for standard grid:
        standard = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ. 7) THEN
c...    D2 temperature (eV):
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ. 8) THEN
c...    D temperature (eV): 
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ. 9) THEN
c...    D2 pressure (mTorr):
        vacuum   = .TRUE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.10) THEN
c...    D pressure (mTorr):
        vacuum   = .TRUE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.11) THEN
c...    D2/D denisity on the standard and vacuum grid:
        standard = .TRUE.
        vacuum   = .TRUE.
      ELSEIF (iopt.EQ.12) THEN
c...    D temperature everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ.13) THEN
c...    D2 temperature everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ.14) THEN
c...    D pressure everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.15) THEN
c...    D2 pressure everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.16) THEN
c...    Ionisation on the vacuum grid:
        standard = .FALSE.
        vacuum   = .TRUE.
      ELSEIF (iopt.EQ.17) THEN
c...    Ionisation everywhere:
        WRITE(0,*) 'WARNING 982:17 - STRANGE NORMALIZATION ON '//
     .             'STANDARD GRID'
        standard = .TRUE.
        vacuum   = .TRUE.
      ELSEIF (iopt.EQ.18) THEN
c...    D pressure on standard grid:
        standard = .TRUE.
        vacuum   = .FALSE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.19) THEN
c...    D2 pressure on standard grid:
        standard = .TRUE.
        vacuum   = .FALSE.
        char(30) = 'UNITS = mTorr'
      ELSEIF (iopt.EQ.20) THEN
c...    Ion temperature in vacuum grid:
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ.21) THEN
c...    Ion density in vacuum grid:
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ.22) THEN
c...    Ion temperature everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.EQ.23) THEN
c...    Ion density everywhere:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = m-3'
      ELSEIF (iopt.EQ.24) THEN
c...    Dgamma:
        standard = .TRUE.
        vacuum   = .FALSE.
        char(30) = 'UNITS = kW m-3'
      ELSEIF (iopt.EQ.25) THEN
c...    D2 temperature on the standard grid:
        standard = .TRUE.
        char(30) = 'UNITS = eV'
      ELSEIF (iopt.GE.26.AND.iopt.LE.31) THEN
c...    Dalpha contributions:
        standard = .TRUE.
        char(30) = 'UNITS = !!!'
c        char(30) = 'UNITS = kW m-3'
      ELSEIF (iopt.EQ.32) THEN
c...    Volume recombination:
        standard = .TRUE.
        char(30) = 'UNITS = m-3 s-1'
      ELSEIF (iopt.GE.33.AND.iopt.LE.44) THEN
c...    Ionisation per stratum:
        standard = .TRUE.
        WRITE(char(28),'(A,I3)') 'STRATUM ',iopt-32
        char(30) = 'UNITS = m-3 s-1'
      ELSEIF (iopt.EQ.45) THEN
c...    ABS(PINQE):
        standard = .TRUE.
        char(30) = 'UNITS = ?'
      ELSEIF (iopt.EQ.46) THEN
c...    Local Mach no.:
        standard = .TRUE.
        char(30) = 'UNITS = ?'

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
            ENDDO
          ENDDO
          CALL CalcFlow(optflow)
        ELSE
          BACKSPACE 5
        ENDIF

      ELSEIF (iopt.EQ.47) THEN
c...    Pressure:
        standard = .TRUE.
        char(30) = 'UNITS = ?'

c...    Find the maximum pressure on each ring:
        DO ir = irsep, nrs
          IF (ir.GT.irtrap) THEN
            pmax(IKLO,ir) = -HI
            DO ik = 1, nks(ir)
              pval = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                            ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH   
              pmax(IKLO,ir) = MAX(pmax(IKLO,ir),pval)
            ENDDO
          ELSE
            ik = ikbound(ir,IKLO)
            pmax(IKLO,ir) = 
     .        CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                     ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH   
            ik = nks(ir) / 2 + 1
            pmax(IKHI,ir) = 
     .        CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                     ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH   
          ENDIF
        ENDDO

      ELSEIF (iopt.EQ.48) THEN
c...    Momentum loss from PIN:
        standard = .TRUE.
        char(30) = 'UNITS = ?'

      ELSEIF (iopt.EQ.49) THEN
c...    Momentum loss from OSM:
        standard = .TRUE.
        char(30) = 'UNITS = ?'

      ELSEIF (iopt.EQ.50) THEN
c...    D / D+:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = ' '

      ELSEIF (iopt.EQ.51) THEN
c...    D2 / D+:
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = ' '

      ELSEIF (iopt.EQ.52) THEN
c...    D2+ / D+:
        standard = .TRUE.
        char(30) = ' '

      ELSEIF (iopt.EQ.53) THEN
c...    Dalpha MFP (D density based):
        standard = .TRUE.
        vacuum   = .TRUE.
        char(30) = 'UNITS = m'

      ELSEIF (iopt.EQ.54) THEN
c...    Impurity density:
        standard = .TRUE.
        vacuum   = .FALSE.
c...    Read ionisation state to plot:       
        READ(graph(14:15),*) izstate
        WRITE(char(29),'(A,I6)') 'ISTATE = ',izstate
        char(30) = 'UNITS  = m-3 s-1'

      ELSEIF (iopt.EQ.55) THEN
c...    Volume recombination + MAR:
        standard = .TRUE.
        vacuum = .TRUE.
        char(30) = 'UNITS = m-3 s-1'

      ELSEIF (iopt.EQ.56) THEN
c...    Facelift photons - transparent:
        standard = .TRUE.
        vacuum = .FALSE.
        char(30) = 'UNITS = m-3 s-1'

      ELSEIF (iopt.EQ.57) THEN
c...    Facelift photons - opaque:
        standard = .TRUE.
        vacuum = .FALSE.
        char(30) = 'UNITS = m-3 s-1'

      ELSEIF (iopt.EQ.60) THEN
c...    SOLPS data:
        standard = .TRUE.
        vacuum = .FALSE.
        char(30) = 'UNITS = m-3 s-1'

      ELSEIF (iopt.GE.80.AND.iopt.LE.89) THEN
c...    Impurity temperature to background ion temperature ratio:
        standard = .TRUE.
        vacuum   = .FALSE.

      ELSEIF (iopt.EQ.90) THEN
c...    ADAS hydrogenic line:
        standard = .TRUE.
        vacuum   = .FALSE.

        mfact = 1.0
        char(30) = 'NO SCALE FACTOR APPLIED'

        CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .             ISELD,IERR)
        CALL LDADAS(1,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >              plrpad,Wlngth,IRCODE)
        osmtmp = 0.0
        CALL RVALKR(osmtmp,plrpad,1,1,1,FT,FP,
     >              MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)

c        DO ir = 1, nrs
c          DO ik = 1, nks(ir)
c            osmtmp(ik,ir) = pinline(ik,ir,6,H_BGAMMA)
c          ENDDO
c        ENDDO

      ELSEIF (iopt.EQ.91) THEN
c...    ADAS impurity line:
        standard = .TRUE.
        vacuum   = .FALSE.

        mfact = 1.0
        char(30) = 'NO SCALE FACTOR APPLIED'

        READ (graph(14:15),*) izmin

        WRITE(0,*) 'CHARGE STATE=',izmin

        IF (izmin.LT.0.OR.izmin.GT.cion) 
     .    CALL ER('982:91','Invalid IZMIN',*99)

c '000 CIII 4650' '*' 96  'pju' 30   80   130  0
c '000 CII  5140' '*' 96  'pju' 34   84   134  0

c '000 CII  6580' '*' 96  'pju' 16   66   116  0     432.3?-16,66,116   
c
c CII 426.6-18,68,118
c CI  911.7- 5,58
c

        CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,ISELE,ISELR,ISELX,
     .             ISELD,IERR)
        CALL LDADAS(cion,IZMIN,ADASID,ADASYR,ADASEX,ISELE,ISELR,ISELX,
     >              plrpad,Wlngth,IRCODE)
        osmtmp = 0.0
        CALL RVALKR(osmtmp,plrpad,1,1,1,FT,FP,
     >              MFACT,XXMIN,XXMAX,YYMIN,YYMAX,VMIN,VMAX)

        WRITE(char(28),'(A,I10  )') 'CHARGE STATE=',izmin
        WRITE(char(29),'(A,F10.2)') 'WAVELENGTH=  ',Wlngth


      ELSEIF (iopt.EQ.200.OR.iopt.EQ.201) THEN
c...    Load external data file and plot:

        scale = 0.0
        WRITE(fname,'(1024X)')
c...    Get the size of the plot, if non-standard:
        READ(5,'(A256)') dummy
        IF   (dummy(8:11).EQ.'File'.OR.dummy(8:11).EQ.'file'.OR.
     .        dummy(8:11).EQ.'FILE') THEN
          READ(dummy,*) cdum1,scale,fname
        ELSE
          CALL ER('982:200','Data file not specified',*99)
        ENDIF
c...
c        nv
c        rv
c        zv
c        nc
c        cq

        IF     (machine.EQ.CMOD) THEN
c...      C-Mod camera:
          WRITE(0,*) 'GRID RESOLUTION SET FOR C-MOD CAMERAS (HARDCODED)'
          deltar = 0.005
          deltaz = 0.005
c          deltar = 0.0025
c          deltaz = 0.0025
        ELSEIF (machine.EQ.DIIID) THEN
c...      DIII-D camera:
          WRITE(0,*) 'GRID RESOLUTION SET FOR DIII-D CAMERAS '//
     .               '(HARDCODED)'
c          deltar = 0.005 / 2.0
c          deltaz = 0.005 / 2.0
c          deltar = 0.01270 / 2.0
c          deltaz = 0.01016 / 2.0
          deltar = 0.02041 / 2.0
          deltaz = 0.02069 / 2.0
        ENDIF

        fp1 = 99
        OPEN(fp1,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .       STATUS='OLD')
c        OPEN(fp1,FILE='/home/steven/work/people/chris/cmod_990429019_'//
c    .               '0950.txt',ACCESS='SEQUENTIAL',STATUS='OLD')
c...    Read header
        DO i1 = 1, 5
          ! jdemod - * by itself in the end position isn't meaningful 
          !          assume it is a * format specifier and move it to 2nd argument
          !READ(fp1,ERR=98,END=98,*)
          READ(fp1,*,ERR=98,END=98)
        ENDDO
c...    Load data:
        nc = 0        
        DO WHILE(.TRUE.) 
          ! jdemod - * by itself in the end position isn't meaningful 
          !          assume it is a * format specifier and move it to 2nd argument
          !READ(fp1,ERR=98,END=9,*) rc(nc+1),zc(nc+1),cq(nc+1)
          READ(fp1,*,ERR=98,END=9) rc(nc+1),zc(nc+1),cq(nc+1)
          nc = nc + 1
c...      Determine cell verticies:
          nv(nc) = 4
          rv(1,nc) = rc(nc) + deltar
          zv(1,nc) = zc(nc) - deltaz
          rv(2,nc) = rc(nc) - deltar
          zv(2,nc) = zc(nc) - deltaz
          rv(3,nc) = rc(nc) - deltar
          zv(3,nc) = zc(nc) + deltaz
          rv(4,nc) = rc(nc) + deltar
          zv(4,nc) = zc(nc) + deltaz
          WRITE(6,*) 'DATA:',nc,rc(nc),zc(nc),cq(nc)
        ENDDO
 9      CLOSE(fp1)
c...    Process data:
        DO i1 = 1, nc
          cq(i1) = MAX(0.0,cq(i1)) 
        ENDDO

        IF     (iopt.EQ.200) THEN
          standard = .FALSE.
          vacuum   = .FALSE.
        ELSEIF (iopt.EQ.201) THEN
          standard = .TRUE.
          vacuum   = .FALSE.
          osmtmp = 0.0
          CALL LoadCameraData(osmtmp,fname,scale)
          nc = 0

          IF (.NOT..TRUE.) THEN
c...        Only keep data in the SOL:
            DO ir = 2, nrs
              DO ik = 1, nks(ir) 
                IF (ir.LT.irsep.OR.ir.GT.irwall-1) osmtmp(ik,ir) = 0.0
              ENDDO
            ENDDO
          ENDIF

        ENDIF

      ELSE
        CALL ER('982','Unknown plot option',*99)
      ENDIF



      XLAB   = '   R (m)'
      YLAB   = '   Z (m)'

c      title  = title(1:77)//'   '
      NVIEW  = '                                    '
      PLANE  = '                                    '
      JOB    = '                                    '
      SMOOTH = '                                    '//
     .         '                                    '
      ANLY   = '                                    '
      REF    = graph(14:LEN_TRIM(graph))


c...  Check if data is to be plotted below a specified Z value only:
c      zlim = -HI
      zlim = HI
      READ(5,'(A80)',END=10) graph6
      IF   (graph6(8:11).EQ.'Zlim'.OR.graph6(8:11).EQ.'ZLIM'.OR.
     .      graph6(8:11).EQ.'zlim') THEN
        READ(graph6,*) cdum1,zlim
      ELSE
        BACKSPACE 5
      ENDIF
10    CONTINUE

c...  Check if data is to be plotted below a specified Z value only:
      rlim = 0.0
      READ(5,'(A80)',END=12) graph6
      IF   (graph6(8:11).EQ.'Rlim'.OR.graph6(8:11).EQ.'RLIM'.OR.
     .      graph6(8:11).EQ.'rlim') THEN
        READ(graph6,*) cdum1,rlim
      ELSE
        BACKSPACE 5
      ENDIF
12    CONTINUE

      IF ((asc_3dmode.EQ.1.OR.asc_3dmode.EQ.2).AND.zval.EQ.-99.0) THEN
c...    Look for zval data in OUT input file:
        READ(5,'(A80)',END=15) graph6
        IF   (graph6(8:11).EQ.'Zval'.OR.graph6(8:11).EQ.'ZVAL'.OR.
     .        graph6(8:11).EQ.'zval') THEN
          READ(graph6,*) cdum1,zval
          IF (eirzaa.LT.0.0.AND.eirzaa.NE.-1.0) THEN
            zval = -zval / eirzaa * 360.0
            WRITE(char(29),'(A,F6.1,A)') 'Tval  = ',zval,' degrees'
          ELSE
            WRITE(char(29),'(A,F6.3,A)') 'Tval  = ',zval,' m'
          ENDIF
        ELSE
          BACKSPACE 5
        ENDIF
15      CONTINUE
      ENDIF

      xval = -99.0
      IF (asc_3dmode.EQ.2) THEN
c...    Look for xval data in OUT input file:
        READ(5,'(A80)',END=20) graph6
        IF   (graph6(8:11).EQ.'Xval'.OR.graph6(8:11).EQ.'XVAL'.OR.
     .        graph6(8:11).EQ.'xval') THEN
          READ(graph6,*) cdum1,xval
          WRITE(char(29),'(A,F6.3,A)') 'Rval  = ',xval,' m'
        ELSE
          BACKSPACE 5
        ENDIF
20      CONTINUE
      ENDIF

c...  Check for vacuum cells that are to be avoided:
25    READ(5,'(A80)',END=30) graph6
      IF   (graph6(8:11).EQ.'Kill'.OR.graph6(8:11).EQ.'KILL'.OR.
     .      graph6(8:11).EQ.'kill') THEN
        READ(graph6,*) cdum1,idum1
        avoid(idum1) = 1
        GOTO 25
      ELSE
        BACKSPACE 5
      ENDIF
30    CONTINUE


      qmin =  HI
      qmax = -HI
      scale = 1.0
      hardscale = .FALSE.
c...  Look for scale data in OUT input file:
      READ(5,'(A80)',END=33) graph6
      IF   (graph6(8:11).EQ.'Scal'.OR.graph6(8:11).EQ.'SCAL'.OR.
     .      graph6(8:11).EQ.'scal') THEN
        READ(graph6,*) cdum1,qmin,qmax,scale
c        WRITE(0,*)
c        WRITE(0,*) '982: NOT USING HARDSCALE WITH SCALING'
c        WRITE(0,*)
        hardscale = .TRUE.
        IF (qmin.EQ.99.0) qmin =  HI
        IF (qmax.EQ.99.0) qmax = -HI
      ELSE
        BACKSPACE 5
      ENDIF
33    CONTINUE


c...  Assign polygons from standard grid:

      IF (standard) THEN



c        DO ir = 1, nrs
c          IF (idring(ir).EQ.BOUNDARY) CYCLE
c          DO ik = 1, nks(ir)
c            IF (virtag(ik,ir).EQ.1) WRITE(0,*) 'DEBUG: VIRTUAL',ik,ir
c            IF (zs(ik,ir).GT.1.5) WRITE(0,*) 'DEBUG: HIGH',ik,ir
c            id = korpg(ik,ir)
c            DO in = 1, 4
c              IF (zvertp(in,id).GT.1.5) WRITE(0,*) 'DEBUG: HIGH P',ik,ir
c            ENDDO
c          ENDDO
c        ENDDO

        IF (eirnsdtor.GT.1.AND.xval.NE.-99.0) THEN
          istart = 1
          iend   = eirnsdtor
        ELSE
          istart = 1
          iend   = 1
        ENDIF

        shift = 0

        IF (eirnsdtor.GT.1.AND.zval.NE.-99.0) THEN
          DO i1 = 1, eirnsdtor
            IF (zval.GE.eirsdtor(i1)) THEN
              shift = MAXBGK * (i1 - 1) 
              WRITE(6,*) '982: STD SHIFT BOOST ',zval,i1,shift
              WRITE(0,*) '982: STD SHIFT BOOST ',zval,i1,shift,
     .          eirsdtor(i1),eirsdtor(i1+1)
            ENDIF
          ENDDO
        ENDIF

        DO cut1 = istart, iend

          irs = 2
          ire = nrs
c...      Just plot data for the core:
          READ(5,'(A80)',END=33) graph6
          IF   (graph6(8:11).EQ.'Core'.OR.graph6(8:11).EQ.'CORE'.OR.
     .          graph6(8:11).EQ.'core') THEN
            ire = irsep -1
          ELSE
            BACKSPACE 5
          ENDIF

          DO ir = irs, ire
            IF (idring(ir).EQ.-1) CYCLE
            DO ik = 1, nks(ir)

              id = korpg(ik,ir)
              status = .FALSE.

              IF (id.EQ.0.OR.nvertp(id).EQ.0) CYCLE

c              IF (zs(ik,ir).LT.zlim) CYCLE
c              IF (zs(ik,ir).GT.zlim) CYCLE
c              IF (rs(ik,ir).LT.rlim) CYCLE

              IF (xval.NE.-99.0) THEN
                xmin =  HI
                ymin =  HI
                xmax = -HI
                ymax = -HI
                DO i1 = 1, nvertp(id)
                  xmin = MIN(xmin,rvertp(i1,id))
                  ymin = MIN(ymin,zvertp(i1,id))
                  xmax = MAX(xmax,rvertp(i1,id))
                  ymax = MAX(ymax,zvertp(i1,id))
                ENDDO

c...            Decide if cell is in the plotting region:
                IF (xval.GE.xmin.AND.xval.LT.xmax) status = .TRUE.

                zmin = eirsdtor(cut1)
                zmax = eirsdtor(MIN(cut1+1,iend))
                IF (cut1.EQ.iend) zmax = eirzaa

                shift = MAXBGK * (cut1 - 1) 
              ELSE
                DO i1 = 1, nvertp(id)
                  x1 = rvertp(i1,id)
                  y1 = zvertp(i1,id)
c...              Decide if cell is in the plotting region:
                  IF (x1.GE.xxmin.AND.x1.LE.xxmax.AND.
     .                y1.GE.yymin.AND.y1.LE.yymax) status = .TRUE.
                ENDDO
              ENDIF

c...          If the cell is in the plotting area, then add
c             the appropriate data to the plot arrays:
              IF (status) THEN

c...            Add cell geometry data:
                IF (xval.NE.-99.0) THEN
                  nc = nc + 1
                  nv(nc) = 4
                  rv(1,nc) = zmax
                  zv(1,nc) = ymin
                  rv(2,nc) = zmin
                  zv(2,nc) = ymin
                  rv(3,nc) = zmin
                  zv(3,nc) = ymax
                  rv(4,nc) = zmax
                  zv(4,nc) = ymax
                ELSE
                  nc = nc + 1
                  nv(nc) = nvertp(id)
                  DO i1 = 1, nvertp(id)
                    rv(i1,nc) = rvertp(i1,id)
                    zv(i1,nc) = zvertp(i1,id)
                  ENDDO
                ENDIF

c...            Assign quantity to be plotted:
                IF     (iopt.EQ.1) THEN
c                IF     (iopt.EQ.1.OR.iopt.EQ.3) THEN
c                  fact = 0.47 / (2.0 * PI * rs(ik,ir) * 0.1)
c                  cq(nc) = pinmol(ik,ir) * fact
                  cq(nc) = pinbgk(ik,ir,15+1+shift)*1.0E+6
c                  IF (ir.EQ.30)
c     .            WRITE(6,*) '-->',ik,ir,pinmol(ik,ir),
c     .             pinbgk(ik,ir,15+1+shift)*1.0E+6
                ELSEIF (iopt.EQ.3) THEN
                  cq(nc) = pinmol(ik,ir) 

                ELSEIF (iopt.EQ.4 .OR.iopt.EQ.6) THEN
                  cq(nc) = pinatom(ik,ir) 

                ELSEIF (iopt.EQ.11) THEN
c...              Need to filter out the noise.
                  cq(nc) = pinmol(ik,ir) / (pinatom(ik,ir)+1.0E-10)
c                  cq(nc) = pinatom(ik,ir)/(pinmol(ik,ir)+1.0E-10)

                ELSEIF (iopt.EQ.12) THEN
c...              D temperature:
c                  IF (ir.GT.irtrap) cq(nc) = pinbgk(ik,ir,10+2+shift)
                  cq(nc) = pinbgk(ik,ir,10+2+shift)

                ELSEIF (iopt.EQ.13.OR.iopt.EQ.25) THEN
c...              D2 temperature:
c                  IF (ir.GT.irtrap) cq(nc) = pinbgk(ik,ir,15+2+shift)
                  cq(nc) = pinbgk(ik,ir,15+2+shift)

                ELSEIF (iopt.EQ.14.OR.iopt.EQ.18) THEN
c...             D pressure:
                 fact = ECH * 1.0E+06 * 7.502
                 cq(nc)=pinbgk(ik,ir,10+2+shift) * 
     .                  pinbgk(ik,ir,10+1+shift) * fact

                ELSEIF (iopt.EQ.15.OR.iopt.EQ.19) THEN
c...             D2 pressure:
                 fact = ECH * 1.0E+06 * 7.502
                 cq(nc)=pinbgk(ik,ir,15+2+shift) *
     .                  pinbgk(ik,ir,15+1+shift) * fact

                ELSEIF (iopt.EQ.17) THEN
c...              Ionisation:
c                  IF (ir.LT.22.OR.ir.GT.24.OR.ik.LT.nks(ir)/2) CYCLE                  
c                  IF (ir.GE.irsep) CYCLE                  
c                  fact = 0.47 / (2.0 * PI * rs(ik,ir) * 0.1)
                  fact = 1.0
                  cq(nc) = ABS(pinion(ik,ir)) * fact

                ELSEIF (iopt.EQ.22) THEN
c...              Ion temperature (below x-point):
c                  IF (zs(ik,ir).LT.zxp) cq(nc) = ktibs(ik,ir)
c                  IF (ir.GT.irtrap.AND.ir.LE.nrs) cq(nc) = ktibs(ik,ir)
                   IF (ir.GE.irsep) cq(nc) = ktibs(ik,ir)

                ELSEIF (iopt.EQ.23) THEN
c...              Ion density:
                  cq(nc) = knbs(ik,ir)

                ELSEIF (iopt.EQ.24) THEN
c...              Dgamma emission:

c...HARDCODED BELOW!

                  cq(nc) = pinline(ik,ir,6,H_BGAMMA)
     .                    
c                  cq(nc) = pinline(ik,ir,6,H_BGAMMA)*1.0E-3*
c     .                    (6.63E-34*3.0E+08)/(4340.0*1.0E-10)



c... THIS SHOULD NOT HAVE BEEN HERE! -MAR 17,2004
c                  IF (ir.GT.irtrap.AND.ik.GT.nks(ir)/1.5) 
c     .              cq(nc)=0.5*cq(nc)

                ELSEIF (iopt.GE.26.AND.iopt.LE.31) THEN
c...              Dalpha emission:
                  cq(nc) = pinline(ik,ir,iopt-25,H_BALPHA)
c                  cq(nc) = pinline(ik,ir,iopt-25,H_BALPHA)*1.0E-3*
c     .                    (6.63E-34*3.0E+08)/(6564.0*1.0E-10)

                ELSEIF (iopt.EQ.32) THEN
c...              Volume recombination:
                  cq(nc) = pinrec(ik,ir)

                ELSEIF (iopt.GE.33.AND.iopt.LE.44) THEN
c...              Ionisation per stratum:
c                  IF (ir.LT.22.OR.ir.GT.24.OR.ik.LT.nks(ir)/2) CYCLE                  
                  cq(nc) =MAX(0.0,pindata(ik,ir,iopt-32))
c                  IF (iopt.EQ.33) cq(nc) =MAX(0.0,pindata(ik,ir,H_ION1))
c                  IF (iopt.EQ.34) cq(nc) =MAX(0.0,pindata(ik,ir,H_ION2))
c                  IF (iopt.EQ.35) cq(nc) =MAX(0.0,pindata(ik,ir,H_ION3))

                ELSEIF (iopt.EQ.45) THEN
c...              Absolute value of power loss from electron channel:
                  cq(nc) = ABS(pinqe(ik,ir))

                ELSEIF (iopt.EQ.46) THEN
c...              Mach number:
                  cq(nc) = kvhs(ik,ir) /
     .                GetCs(ktebs(ik,ir),ktibs(ik,ir)) / qt

                ELSEIF (iopt.EQ.47) THEN
c...              Mach number:
                  IF (ir.LT.irsep) CYCLE
c                  IF (ir.LT.irsep.OR.
c     .                ir.LT.irtrap.AND.ik.GT.nks(ir)/2) CYCLE
                  cq(nc) = 
     .              CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                           ktibs(ik,ir),kvhs (ik,ir)/qt) * ECH 
c...              Scale pressure:
                  IF ((ir.LT.irwall.AND.ik.LE.nks(ir)/2).OR.
     .                 ir.GT.irtrap) THEN
                    cq(nc) = cq(nc) / pmax(IKLO,ir)
                  ELSE
                    cq(nc) = cq(nc) / pmax(IKHI,ir)                         
                  ENDIF              
                  cq(nc) = MIN(cq(nc),1.0)

                ELSEIF (iopt.EQ.48) THEN
c...              Momentum loss from EIRENE:
                  cq(nc) = pinmp(ik,ir)

                ELSEIF (iopt.EQ.49) THEN
c...              Momentum loss from OSM:
                  IF     (ir.LT.irwall.AND.ik.LE.ikbound(ir,IKLO)) THEN
                    p1 = CalcPressure(knbs (ik,ir),ktebs(ik,ir),
     .                                ktibs(ik,ir),kvhs (ik,ir)/qt)
                    p2 = CalcPressure(knbs (ik+1,ir),ktebs(ik+1,ir),
     .                                ktibs(ik+1,ir),kvhs (ik+1,ir)/qt)
                    cq(nc) = (p2 - p1) * ECH / (kss(ik+1,ir)-kss(ik,ir))
c                    IF (ir.EQ.irsep) WRITE(0,*) 'MOM:',cq(nc)
c     .                
                  ELSEIF (ir.LT.irwall.AND.ik.GT.ikbound(ir,IKLO).AND.
     .                    ik.LT.nks(ir)/2) THEN
                  ELSE
                    cq(nc) = osmmp(ik,ir) 
                  ENDIF

                ELSEIF (iopt.EQ.50) THEN
c...               D / D+:
                   cq(nc)=pinbgk(ik,ir,10+1+shift)*1.0E+6 / knbs(ik,ir)

                ELSEIF (iopt.EQ.51) THEN
c...               D2 / D+:
                   cq(nc)=pinbgk(ik,ir,15+1+shift)*1.0E+6 / knbs(ik,ir)

                ELSEIF (iopt.EQ.52) THEN
c...               D2 / D+:
                   cq(nc)=pinmoi(ik,ir) / knbs(ik,ir)

                ELSEIF (iopt.EQ.53) THEN
c...               Dalpha MFP (D density based):
                   IF (ir.LT.irsep) CYCLE
                   Ti = ktibs(ik,ir)
                   nD = pinbgk(ik,ir,10+1+shift)*1.0E+6/1.0E+20
                   IF (Ti.EQ.0.0.OR.nD.LT.0.001) CYCLE
                   cq(nc) = 2.2E-4 * SQRT(Ti / 1.0) / nD

                ELSEIF (iopt.EQ.54) THEN
c...              Impurity density:
                  cq(nc) = MAX(0.0,sdlims(ik,ir,izstate))

                ELSEIF (iopt.EQ.55) THEN
c...              Volume recombination + MAR:

                  ti = ktibs(ik,ir)
                  ne = knbs(ik,ir)
c...              D2+ equilibrium density:
                  sigmav = GetEAD(ti,ne,22,'H.11')
                  nD2p = sigmav * pinmol(ik,ir)
c...              MAR: e + H2+ -> e + H + H
                  sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...              MAD: e + H2+ -> e + H + H+
                  sigmav2 = GetEAD(ti,ne,15,'H.4 ')
                  sigmav = sigmav1**2 / (sigmav1 + sigmav2)
                  cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0

                  IF (ir.LT.irwall) THEN
                    MARsum1 = MARsum1 + cq(nc) * kvols(ik,ir)
                  ELSE
                    MARsum2 = MARsum2 + cq(nc) * kvols(ik,ir)
                  ENDIF

c                  WRITE(0,*) 'MARSUM:',marsum

c...              Direct recombination:
                  cq(nc) = cq(nc) + pinrec(ik,ir)

                ELSEIF (iopt.EQ.56) THEN
c...              Facelift - transparent:
                  cq(nc) = eirpho1(ik,ir)

                ELSEIF (iopt.EQ.57) THEN
c...              Facelift - opaque:
c                  IF (eirpho2(ik,ir).LE.0.0) THEN
c                    cq(nc) = -eirpho2(ik,ir)
c                  ELSE 
c                    cq(nc) = 0.0
c                    WRITE(0,*) 'PROBLEM:',ik,ir,eirpho2(ik,ir)
c                  ENDIF
                   cq(nc) = eirpho2(ik,ir)

                ELSEIF (iopt.EQ.60) THEN
c...              SOLPS data:
                  cq(nc) = e2dion(ik,ir)

                ELSEIF (iopt.GE.80.AND.iopt.LE.89) THEN
c...              Impurity temperature ratio for states 0 through 9:
                  cq(nc) = sdts(ik,ir,iopt-80) / (ktebs(ik,ir) + LO)

                ELSEIF (iopt.EQ.90.OR.iopt.EQ.91) THEN
c...              ADAS hydrogenic and impurity lines:
                  cq(nc) = osmtmp(ik,ir)

                ELSEIF (iopt.EQ.201) THEN
c...              Interpolated Dgamma data from Chris:
                  cq(nc) = osmtmp(ik,ir)

                ELSE
                  CALL ER('982','Unknown quantity to be plotted',*99)
                ENDIF
              ENDIF
            ENDDO    

c...      END of IR loop:
          ENDDO

c...    END of CUT1 loop:
        ENDDO

      ENDIF

c...  Assign polygons from vacuum grid:

      set = 2

      IF (vacuum.AND.asc_ncell.GT.0) THEN
c      IF (vacuum) THEN

        shift = 1 + eirnpgdat

        IF (asc_3dmode.EQ.2.AND.zval.NE.-99.0) THEN
          DO i1 = 1, ascncut
c            WRITE(0,*) 'CHECKING:',asc_zmin3D(i1),
c     .                             asc_zmax3D(i1),eirzaa
            IF (zval.GE.asc_zmin3D(i1).AND.zval.LT.asc_zmax3D(i1)) THEN
              shift = shift + asc_ncell * (i1-1)
              WRITE(6,*) '982: VAC SHIFT BOOST ',i1-1
              EXIT
            ENDIF
          ENDDO
        ENDIF

        IF (xval.NE.-99.0) THEN
          istart = 1
          iend   = asc_ncell * ascncut
        ELSE
          istart = 1
          iend   = asc_ncell
        ENDIF

        pavg = 0.0
        vtot = 0.0

        DO i1 = istart, iend

c...      If the toroidal extent of the cell is outside the specified z-value,
c         then don't plot the cell:
c          IF (asc_3dmode.EQ.1.AND.(asc_zmin3D(i1).GT.zval.OR.
c     .                             asc_zmax3D(i1).LE.zval)) CYCLE

c          fact = ECH * 1.0E+06 * 0.67 * 7.502
c          WRITE(6,*) '--?',i1,i1+shift,pinasd(i1+shift,2,3,set) * fact
            
c...      See if cell is blacklisted:
          IF (avoid(i1).EQ.1) CYCLE

c...      Find x,ymin and x,ymax for the cell:
          xmin =  HI
          ymin =  HI
          xmax = -HI
          ymax = -HI
          i3 = MOD(i1,asc_ncell)            
          IF (i3.EQ.0) i3 = asc_ncell
          DO i2 = 1, ascnvertex(i3)
            xmin = MIN(xmin,ascvertex(i2*2-1,i3))
            ymin = MIN(ymin,ascvertex(i2*2  ,i3))
            xmax = MAX(xmax,ascvertex(i2*2-1,i3))
            ymax = MAX(ymax,ascvertex(i2*2  ,i3))
          ENDDO

c...      If the cell is above the specified Z limit (DIVIMP
c         coordinate system now), then skip the cell:
          IF (0.5*(ymin+ymax).GT.zlim) CYCLE

          IF (zlim.NE.-HI.AND.0.5*(ymin+ymax).GT.zlim) CYCLE
          IF (0.5*(xmin+xmax).LT.rlim) CYCLE

          IF (xval.NE.-99.0) THEN
c...        Toroidal plot specified:
            IF (.NOT.(xmin.LE.xval.AND.xmax.GT.xval)) CYCLE

c...        Find zmin and zmax for the cell:
            i3 = INT(REAL(i1) / (REAL(asc_ncell) + 0.5)) + 1
            zmin = asc_zmin3D(i3)
            zmax = asc_zmax3D(i3)

            nc = nc + 1
            nv(nc) = 4
            rv(1,nc) = zmax
            zv(1,nc) = ymin
            rv(2,nc) = zmin
            zv(2,nc) = ymin
            rv(3,nc) = zmin
            zv(3,nc) = ymax
            rv(4,nc) = zmax
            zv(4,nc) = ymax
          ELSE
            nc = nc + 1
            nv(nc) = ascnvertex(i1)
            xcen = 0.0
            ycen = 0.0
            DO i2 = 1, ascnvertex(i1)
              rv(i2,nc) = ascvertex(i2*2-1,i1)
              zv(i2,nc) = ascvertex(i2*2  ,i1)
c...          Estimate the center of the cell:
              xcen = xcen + ascvertex(i2*2-1,i1)
              ycen = ycen + ascvertex(i2*2  ,i1)
            ENDDO
            xcen = xcen / REAL(ascnvertex(i1))
            ycen = ycen / REAL(ascnvertex(i1))
          ENDIF

c...      Store the additional cell region of the cell for averaging limitations
c         in the contour plot:
          rq(nc) = asc_region(i3)

c...      PINASD(index,quantity,species,relaxation):        
          IF     (iopt.EQ.1 .OR.iopt.EQ.2 ) THEN
c...        D2 density:
            cq(nc) = pinasd(i1+shift,1,3,set)*1.0E+06

          ELSEIF (iopt.EQ.4 .OR.iopt.EQ.5 ) THEN
c...        D density:
            cq(nc) = pinasd(i1+shift,1,1,set)*1.0E+06

          ELSEIF (iopt.EQ.7 .OR.iopt.EQ.13) THEN
c...        D2 temperature:
c            IF (rq(nc).NE.3.0)  
c     .        cq(nc) = pinasd(i1+shift,2,8,set)
            cq(nc) = pinasd(i1+shift,2,8,set)
c            cq(nc) = pinasd(i1+shift,2,3,set)/
c     .              (pinasd(i1+shift,1,3,set)+1.0E-10)*0.67

          ELSEIF (iopt.EQ.8 .OR.iopt.EQ.12) THEN
c...        D temperature:
c            IF (rq(nc).NE.3.0)  
c     .        cq(nc) = pinasd(i1+shift,2,7,set)
            cq(nc) = pinasd(i1+shift,2,7,set)
c            cq(nc) = pinasd(i1+shift,2,1,set)/
c     .              (pinasd(i1+shift,1,1,set)+1.0E-10)*0.67

          ELSEIF (iopt.EQ. 9.OR.iopt.EQ.15) THEN
c...        D2 pressure (mTorr):
            fact = ECH * 1.0E+06 * 7.502
            cq(nc) = pinasd(i1+shift,2,8,set) *
     .               pinasd(i1+shift,1,8,set) * fact
            IF (xcen.GT. 0.57.AND.xcen.LT. 0.64.AND.
     .          ycen.GT.-0.68.AND.ycen.LT.-0.62.AND.
     .          cq(nc).GT.1.0) THEN
              pavg = pavg + cq(nc) * asc_vol(i1)
              vtot = vtot + asc_vol(i1)
c              WRITE(0,*) '>>>',xcen,ycen,cq(nc),asc_vol(i1)
            ENDIF
c            fact = ECH * 1.0E+06 * 0.67 * 7.502
c            cq(nc) = pinasd(i1+shift,2,3,set) * fact

          ELSEIF (iopt.EQ.10.OR.iopt.EQ.14) THEN
c...        D pressure:
            fact = ECH * 1.0E+06 * 7.502
            cq(nc) = pinasd(i1+shift,2,7,set) *
     .               pinasd(i1+shift,1,7,set) * fact
c            fact = ECH * 1.0E+06 * 0.67 * 7.502
c            cq(nc) = pinasd(i1+shift,2,1,set) * fact

          ELSEIF (iopt.EQ.11) THEN
c...        D2:D ratio:
            cq(nc) =  pinasd(i1+shift,1,3,set)*1.0E+06 /
     .               (pinasd(i1+shift,1,1,set)*1.0E+06 + 1.0E-10)
c            cq(nc) =  pinasd(i1+shift,1,1,set)*1.0E+06 / 
c     .               (pinasd(i1+shift,1,3,set)*1.0E+06 + 1.0E-10)

          ELSEIF (iopt.EQ.16.OR.iopt.EQ.17) THEN
c...        Ionisation:
c            CYCLE
            cq(nc) = ABS(ascdata(i1+shift-1-eirnpgdat,1))
            IF (i1.EQ.788.OR.i1.EQ.787) cq(nc) = cq(nc) * 0.2

          ELSEIF (iopt.EQ.20.OR.iopt.EQ.22) THEN
c...        Ion temperature in additional cells:
            cq(nc) = pinasd(i1+shift,2,5,1)

          ELSEIF (iopt.EQ.21.OR.iopt.EQ.23) THEN
c...        Ion density in additional cells:
            cq(nc) = pinasd(i1+shift,1,5,1)*1.0E+06

          ELSEIF (iopt.EQ.50) THEN
c...        D / D+:
            IF (pinasd(i1+shift,1,5,1).GT.1.0E+10) THEN
              cq(nc)=pinasd(i1+shift,1,1,set) / pinasd(i1+shift,1,5,1)
            ENDIF
          ELSEIF (iopt.EQ.51) THEN
c...        D2 / D+:
            IF (pinasd(i1+shift,1,5,1).GT.1.0E+10) THEN
              cq(nc)=pinasd(i1+shift,1,3,set) / pinasd(i1+shift,1,5,1)
            ENDIF

          ELSEIF (iopt.EQ.53) THEN
c...        Dalpha MFP (D density based):
            Ti = pinasd(i1+shift,2,5,1)
            nD = pinasd(i1+shift,1,1,set)*1.0E+06/1.0E+20 
            IF (Ti.EQ.0.0.OR.nD.LT.0.001) CYCLE
            cq(nc) = 2.2E-4 * SQRT(Ti / 1.0) / nD

          ELSEIF (iopt.EQ.55) THEN
c...        Volume recombination + MAR:

            ti = pinasd(i1+shift,2,5,1)
            ne = pinasd(i1+shift,1,5,1)*1.0E+06

            IF (ne.LT.1.0E+15) CYCLE

c...        D2+ equilibrium density:
            sigmav = GetEAD(ti,ne,22,'H.11')
            nD2p = sigmav * pinasd(i1+shift,1,3,set) * 1.0E+06
c...        MAR: e + H2+ -> e + H + H
            sigmav1 = GetEAD(ti,ne,17,'H.4 ')
c...        MAD: e + H2+ -> e + H + H+
            sigmav2 = GetEAD(ti,ne,15,'H.4 ')
            sigmav = sigmav1**2 / (sigmav1 + sigmav2)
            cq(nc) = sigmav * 1.0E-06 * ne * nD2p * 2.0
c...        * Estimated! *
            cellvol = (xmax - xmin) * (ymax - ymin) *
     .                2.0 * PI * 0.5 * (xmin + xmax) 
c            WRITE(0,*) 'estimated cellvol:',cellvol

            MARsum2 = MARsum2 + cq(nc) * cellvol

c...        Direct recombination:
            sigmav = GetEAD(ti,ne,3,'H.4 ')            
            cq(nc) = cq(nc) + (sigmav * 1.0E-06 * ne) * ne

          ENDIF

c...      Output data to .OUT file:
          WRITE(6,'(A,3I6,1P,E12.4,0P)') 'VAC:',i1,i1+shift,nc,cq(nc)

        ENDDO

        IF (pavg.NE.0.0) THEN
          WRITE(0,*) 'AVERAGE PRESSURE',pavg/vtot
        ENDIF 

      ENDIF

c...  Scale the data:
      DO i1 = 1, nc
        cq(i1) = cq(i1) * scale
      ENDDO

c...  Find plotting quantity minimum and maximum:
      nconts = 10
      setqmin = .TRUE.
      setqmax = .TRUE.
      IF (qmin.NE. HI) setqmin = .FALSE.
      IF (qmax.NE.-HI) setqmax = .FALSE.
      IF (setqmin.OR.setqmax) THEN
        DO i1 = 1, nc
c...      Decide if the cell is within the viewing range:
          inside = .FALSE.
          DO i2 = 1, nv(i1)
            IF (rv(i2,i1).GE.xxmin.AND.rv(i2,i1).LE.xxmax.AND.
     .          zv(i2,i1).GE.yymin.AND.zv(i2,i1).LE.yymax) inside=.TRUE.
          ENDDO
          IF (inside) THEN
           IF (setqmin.AND.cq(i1).NE.0.0.AND.qmin.GT.cq(i1)) qmin=cq(i1)
           IF (setqmax.AND.cq(i1).NE.0.0.AND.qmax.LT.cq(i1)) qmax=cq(i1)
          ENDIF
        ENDDO
        IF (setqmin.AND.qmin.GT.-1.0E-10) qmin = MAX(qmin,0.01*qmax)
      ENDIF

      miny =  HI
      maxy = -HI

c...  Get the size of the plot, if non-standard:
      READ(5,'(A256)') dummy
      IF   (dummy(8:11).EQ.'Size'.OR.dummy(8:11).EQ.'size'.OR.
     .      dummy(8:11).EQ.'SIZE') THEN
        READ(dummy,*) cdum1,map1x,map2x,map1y,map2y
      ELSE
        BACKSPACE 5
      ENDIF


c...  Check if a contour plot is to be plotted:
      rlim = 0.0
      READ(5,'(A80)') graph6
      IF   (graph6(8:11).EQ.'Cont'.OR.graph6(8:11).EQ.'CONT'.OR.
     .      graph6(8:11).EQ.'cont') THEN
        BACKSPACE 5
        GOTO 100
      ELSE
        BACKSPACE 5
        CALL CustomisePlot(title,xlab,ylab,elabs)
      ENDIF


      CALL GRTSET (TITLE,REF,NVIEW,PLANE,JOB,XXMIN,XXMAX,
     .             YYMIN,YYMAX,TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,NRS)


c...  Draw polygons:
      WRITE(6,*) 'QMIN:',qmin,qmax
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)
      CALL HSV
      DO i1 = 1, nc
        IF (cq(i1).GE.qmin) THEN
          frac = (cq(i1) - qmin) / (qmax - qmin)
          frac5 = 100.0*frac
          fmod5 = AMOD(frac5,2.0)
          frac = MIN(0.98,(frac5-fmod5)/100.0)
          CALL SetCol255(frac,qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255) 
          WRITE(6,*) 'PLOT:',rv(1,i1),zv(1,i1),cq(i1),nv(i1)
          WRITE(6,*) '    :',rv(2,i1),zv(2,i1)
          WRITE(6,*) '    :',rv(3,i1),zv(3,i1)
          WRITE(6,*) '    :',rv(4,i1),zv(4,i1)
          CALL PTPLOT(rv(1,i1),zv(1,i1),1,nv(i1),1)
        ENDIF
      ENDDO


      GOTO 200


100   CONTINUE


c
c        Read in contour options if contour plot was specified  
c
         call rdg_contopts(graph,icntr,ncntr,uconts,maxpts,
     >                     xcen,ycen,xnear,ynear,ierr)

      

      colscale = .TRUE.

c        xres = 19
c        yres = 19

        xres = 99
        yres = 99

c        xres = 79
c        yres = 79

c        xres = 129
c        yres = 129

        
        READ(5,'(A80)') graph6
        IF   (graph6(8:11).EQ.'Cres'.OR.graph6(8:11).EQ.'CRES'.OR.
     .        graph6(8:11).EQ.'cres') THEN
          READ(graph6,*) cdum1,xres,yres
        ELSE
          BACKSPACE 5
        ENDIF


c
c        Plot the image as a contour plot
c
         maxix = xres
         maxiy = yres
         nix   = xres
         niy   = yres       
c
c         write(0,'(a)') 'Plotting contours'
         if (ierr.ne.0) 
     .     CALL ER('982','Problem loading contour plot options',*99)

         IF (ncntr.EQ.1.AND.uconts(1).EQ.-1.0) THEN
c           ncntr = 12
           ncntr = 29
c           ncntr = 19
c           ncntr = 79
           DO i1 = 1, ncntr

             frac = REAL(i1-1) / REAL(ncntr)

c             frac = REAL(i1) / REAL(ncntr + 1)
c             uconts(i1) = REAL(i1) / REAL(ncntr + 1)

             uconts(i1) = qmin + frac * (qmax - qmin)
             WRITE(6,*) 'uconts:',uconts(i1)
           ENDDO
         ENDIF

c
c        Set up colours based on contour options  
c
c         if (icntr.eq.0.or.icntr.eq.2) then
c             call setup_col(ncntr+1,2)
c         elseif (icntr.eq.1) then 
c             ncntr = 10
c             call setup_col(ncntr+1,2)
c         elseif (icntr.eq.3.or.icntr.eq.4) then
c             ncntr = ncntr
c             call setup_col(ncntr+1,2)
c         endif
c
c         write(0,'(a,2i4,8(1x,g12.5))') 'PLOT_CONT:',icntr,ncntr,
c     >         xcen,ycen,xnear,ynear,maxval,minval
c
c        Set NVIEW, PLANE ...
c
c         NVIEW  = 'IMAGE CONTOUR PLOT'
c         write(plane,'(a,f5.2,a)') 'IMAGE PLANE AT ',
c     >                  vect_mag(direction),' (M)'
         PLANE = '                                    '

         smooth = ' '
         anly   = ' '
c
c        Set min and max levels for contour plot
c
         maxscale = qmax
         minscale = qmin

c         WRITE(0,*) 'MIN,MAXSCALE:',minscale,maxscale

c
c
c        Set up axes for the contour plot call - use the 
c        base camera definition space at the given distance.
c
c
c        Allocate space for axes and check to see if space available
c

         allocate (image(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'IMAGE array could not be allocated: ',xres
            write(6,*) 'IMAGE array could not be allocated: ',xres
            goto 99 
         endif        

         allocate (image1(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'IMAGE1 array could not be allocated: ',xres
            write(6,*) 'IMAGE1 array could not be allocated: ',xres
            goto 99 
         endif        

         allocate (tpt(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'TPT array could not be allocated: ',xres
            write(6,*) 'TPT array could not be allocated: ',xres
            goto 99 
         endif        

         allocate (tpr(xres,yres),STAT=ierr)
         if (ierr.ne.0) then 
            write(0,*) 'TPR array could not be allocated: ',xres
            write(6,*) 'TPR array could not be allocated: ',xres
            goto 99 
         endif        

         allocate (raxis(xres),STAT=ierr)
c
         if (ierr.ne.0) then 
            write(0,*) 'RAXIS array could not be allocated: ',xres
            write(6,*) 'RAXIS array could not be allocated: ',xres
            goto 99 
         endif        
c
         allocate (zaxis(yres),STAT=ierr)
c
         if (ierr.ne.0) then 
            write(0,*) 'ZAXIS array could not be allocated: ',yres
            write(6,*) 'ZAXIS array could not be allocated: ',yres
            goto 99 
         endif        


         CALL IZero(tpt,xres*yres)
         CALL IZero(tpr,xres*yres)
 
c
c        Assign axis values
c

         WRITE(6,*) 'XXMIN,MAX:',xxmin,xxmax
         WRITE(6,*) 'YYMIN,MAX:',yymin,yymax

         dr = (xxmax - xxmin) / xres

c         view_width = xxmax - xxmin
c         dr = view_width/xres
c
c        Column center coordinates
c
         do ic = 1, xres
           raxis(ic) = xxmin + REAL(ic - 1) * dr + 0.5 * dr

c            raxis(ic) = -(view_width/2.0) + dr/2.0 + (ic-1)*dr 
            write(6,'(a,i5,2(1x,g12.5))') 'RAXIS:',ic,raxis(ic),dr 
         end do 
c
c        Row center coordinates
c

         dz = (yymax - yymin) / REAL(yres)

c         view_height = yymax - yymin
c         dz = view_height/yres

         do ir = 1,yres  
            zaxis(ir) = yymin + REAL(ir - 1) * dz + 0.5 * dz 

c            zaxis(ir) = (view_height/2.0) 
c     >                         - dz/2.0 - (ir-1)*dz 

            write(6,'(a,i5,2(1x,g12.5))') 'ZAXIS:',ir,zaxis(ir),dz 
         end do 


c...     Assign data to image array from 2D cell plot:



c...     Method2: Some weighted average of nearby data points.  If the XY cell center is not 
c        inside the grid, then assign a null value:

         DO ir = 1, yres
           DO ic = 1, xres
 
             tpt(ic,ir) = 0
             
             xcen = raxis(ic)
             ycen = zaxis(ir)
             centerinside = .FALSE.
             i1 = 0
             DO WHILE(.NOT.centerinside.AND.i1.LT.nc)
               i1 = i1 + 1

               IF (cq(i1).EQ.0.0) CYCLE
c               IF (cq(i1).LT.qmin) CYCLE
c...           Decide if cell center is inside or outside the parent cell polygon:
               icount = 0
               DO i2 = 1, nv(i1)
                 i3 = i2 + 1
                 IF (i2.EQ.nv(i1)) i3 = 1
                 CALL CalcInter
     .             (DBLE(rv(i2,i1)),DBLE(zv(i2,i1)),
     .              DBLE(rv(i3,i1)),DBLE(zv(i3,i1)),
     .              DBLE(xcen),DBLE(ycen),DBLE(xcen+1000.0),DBLE(ycen),
     .              t1,t2)
                 IF (t1.GT.0.0D0.AND.t1.LT.1.0D0.AND.
     .               t2.GT.0.0D0.AND.t2.LT.1.0D0) icount = icount + 1
               ENDDO
               IF (icount.EQ.0.OR.MOD(icount,2).EQ.0) THEN
               ELSE
                 centerinside = .TRUE.
                 tpr(ic,ir) = rq(i1)
               ENDIF
             ENDDO

             IF (centerinside) tpt(ic,ir) = i1

           ENDDO
         ENDDO


c...     Find the center point of each cell:
         DO i1 = 1, nc
           rc(i1) = 0.0
           zc(i1) = 0.0
           DO i2 = 1, nv(i1)
             rc(i1) = rc(i1) + 1.0 / REAL(nv(i1)) * rv(i2,i1)
             zc(i1) = zc(i1) + 1.0 / REAL(nv(i1)) * zv(i2,i1)
           ENDDO
         ENDDO           

c...     Find the 5 closest points with CQ greater than QMIN and 
c        average:

c         WRITE(0,*) '***********************'
c         WRITE(0,*) ' CONTOUR PLOT MODIFIED'
c         WRITE(0,*) '***********************'

         DO ir = 1, yres
           DO ic = 1, xres

             image(ic,ir) = 0.0

             WRITE(6,'(A,4I6)') 'IMAGE:',ic,ir,tpt(ic,ir),tpr(ic,ir)
c             WRITE(0,'(A,4I6)') 'IMAGE:',ic,ir,tpt(ic,ir),tpr(ic,ir)
   
             IF (tpt(ic,ir).EQ.0) CYCLE

             xcen = raxis(ic)
             ycen = zaxis(ir)

c...TEMP!
c             IF (ycen.GT.-0.35) CYCLE
c             IF (ycen.LT.-0.65) CYCLE
c             IF (xcen.GT. 0.80) CYCLE
c             IF (xcen.LT. 0.70) CYCLE

             npt = 1
             vpt(npt) = 0.0
             dpt(npt) = HI

             DO i1 = 1, nc

c *TEMP*
c               IF (cq(i1).LT.0.01*qmin) CYCLE
c               IF (cq(i1).GT.qmax) CYCLE

c...           Don't average cells across different additional cell
c              regions:
               IF (rq(i1).NE.0.AND.tpr(ic,ir).NE.0.AND.
     .             rq(i1).NE.tpr(ic,ir)) CYCLE

               dist = SQRT((xcen - rc(i1))**2 + (ycen - zc(i1))**2)

c IS THIS OKAY?
              IF (dist.LT.dr+dz.OR.
c *TEMP*
     .            (machine.EQ.DIIID.AND.dist.LT.0.03).OR.
     .             (rq(i1).NE.0.AND.dist.LT.0.01)) THEN

c               IF (dist.LT.dr+dz) THEN
                 IF (npt.LT.100) npt = npt + 1
c               IF (dist.LT.dpt(npt)) THEN
c                 IF (npt.LT.5) npt = npt + 1
                 dpt(npt) = dist 
                 vpt(npt) = MAX(0.0,cq(i1))


c...             Sort data points from closest to farthest:                 
                 i2 = 1
                 DO WHILE(i2.LT.npt)
                   i2 = i2 + 1
                   IF (dpt(i2-1).GT.dpt(i2)) THEN
                     tmpdpt = dpt(i2-1)
                     tmpvpt = vpt(i2-1)
                     dpt(i2-1) = dpt(i2)
                     vpt(i2-1) = vpt(i2)
                     dpt(i2) = tmpdpt
                     vpt(i2) = tmpvpt
                     i2 = 1   
                   ENDIF
                 ENDDO
               ENDIF

             ENDDO

c...         No value identified for cell:
             IF (xval.NE.-99.OR.(npt.EQ.1.AND.vpt(1).EQ.0.0)) THEN
               npt = 1
               vpt(1) = cq(tpt(ic,ir))
               dpt(1) = 1.0
             ENDIF

c...         Average the data points:

             totdpt = 0.0
             DO i1 = 1, npt
c               dpt(i1) = 1.0 / (dpt(i1)**5)

               IF (dpt(i1).LT.1.0E-10) THEN
                 dpt(i1) = 1.0E+20
               ELSE
c*TEMP*
c               dpt(i1) = 1.0 / dpt(i1)
                 dpt(i1) = 1.0 / (dpt(i1)**2)
               ENDIF                 

               totdpt = totdpt + dpt(i1)
               WRITE(6,*) '     :',dpt(i1),vpt(i1)
             ENDDO

c            STOP 'dfsdsd'

             IF (totdpt.GT.0.0) THEN
               DO i1 = 1, npt
                 WRITE(6,'(A,2I6,1P,3E10.2,0P)') 
     .             'VALUE:',ic,ir,dpt(i1)/totdpt,vpt(i1),image(ic,ir)

                 image(ic,ir) = image(ic,ir) + dpt(i1)/totdpt * vpt(i1)
               ENDDO
             ELSE
               image(ic,ir) = 0.0
                 WRITE(6,'(A,2I6,A)')
     .             'VALUE:',ic,ir,' NO DATA?'
             ENDIF


           ENDDO                          
         ENDDO

175      CONTINUE


c
c        Draw contour plot 
c
         CALL CustomisePlot(title,xlab,ylab,elabs)

         CALL CTRMAG(10)
 
c         IF (.TRUE.) THEN
c           map1x = 0.1
c           map2x = 0.1 + 0.3
c           map1y = 0.4
c           map2y = 0.4 + 0.3
c         ENDIF

c...     Have GRTSET_TRIM called:

         slopt4 = 1
         char(30) = ' '
         REF    = '                                    '

         CALL GRTSET (TITLE,REF,NVIEW,PLANE,blabs(1),
     .         xXMIN,xXMAX,
     >         yYMIN,yYMAX,
     .         TABLE,XLAB,YLAB,2,SMOOTH,1,ANLY,Ncntr)
 
c
c        Avoid drawing cell boundaries 
c
c         slopt = 4
c


c...dev
          CALL HSV

c...      Clean up IMAGE array:
          DO i1 = 1, xres
            DO i2 = 1, yres
              IF (image(i1,i2).GT.-1.0E-10.AND.
     .            image(i1,i2).LT. 1.0E-10) image(i1,i2) = 0.0
            ENDDO
          ENDDO

c...      Erase boarder:
          DO i1 = 1, xres
            image(i1,1) = 0.0
c            image(i1,2) = 0.0
c            image(i1,3) = 0.0
            image(i1,yres) = 0.0
c            image(i1,yres-1) = 0.0
c            image(i1,yres-2) = 0.0
          ENDDO
          DO i2 = 1, yres
            image(1,i2) = 0.0
c            image(2,i2) = 0.0
c            image(3,i2) = 0.0
            image(xres,i2) = 0.0
c            image(xres-1,i2) = 0.0
c            image(xres-2,i2) = 0.0
          ENDDO

c...      Erase boarder:
          DO i2 = 1, yres
            DO i1 = 1, xres
              WRITE(6,*) 'CDATA:',i1,i2,image(i1,i2) 
            ENDDO
          ENDDO


c...      Sort contours in increasing absolute value:
          status = .TRUE.
          DO WHILE (status)
            status = .FALSE.
            DO i1 = 1, ncntr-1
              IF (ABS(uconts(i1)).GT.ABS(uconts(i1+1))) THEN
                rdum1 = uconts(i1+1)
                uconts(i1+1) = uconts(i1)
                uconts(i1) = rdum1
                status = .TRUE.
              ENDIF
            ENDDO
          ENDDO

          DO i1 = 1, ncntr

            frac = (uconts(i1) - qmin) / (qmax - qmin)

            frac5 = 100.0*frac
            fmod5 = AMOD(frac5,2.0)

            frac = MIN(0.98,(frac5-fmod5)/100.0)

            CALL SetCol255(frac,qmin,qmax)

            map1x_d = map1x
            map2x_d = map2x
            map1y_d = map1y
            map2y_d = map2y

            tag_d = 1

c...        Copy IMAGE to IMAGE1 and account for 
c           negative contours:
            IF (uconts(i1).LT.0.0) THEN
              DO i2 = 1, xres
                DO i3 = 1, yres
                  IF (image(i2,i3).LT.-1.0E-10) THEN
                    image1(i2,i3) = -image(i2,i3)
                  ELSE
                    image1(i2,i3) = 0.0
                  ENDIF
                ENDDO
              ENDDO
            ELSE       
              DO i2 = 1, xres
                DO i3 = 1, yres
                  image1(i2,i3) = image(i2,i3)
                ENDDO
              ENDDO
            ENDIF

c...        Plot contour:
            CALL GRCONT (image1,1,xres,xres,1,yres,yres,
     >                   ABS(uconts(i1)),raxis,zaxis,'hello')

            tag_d = 0
          ENDDO

c...      Finish off the plot boarder:
          CALL FULL
          CALL LINCOL(defcol)
          CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
          CALL MAP(0.0,1.0,0.0,1.0)
          CALL POSITN (0.0,1.0)
          CALL JOIN   (1.0,1.0)
          CALL JOIN   (1.0,0.0)

         deallocate(raxis)
         deallocate(zaxis)
         deallocate(image)
         deallocate(image1)
         deallocate(tpt)



200   CONTINUE





      IF (elabs(1)(1:4).EQ.'none') GOTO 210




c...  Draw scale:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL HSV
      DO i1 = 98, 0, -2
c      DO i1 = 95, 0, -4

        frac = REAL(i1) / 100.0

        CALL SetCol255(frac,qmin,qmax)
       
        CALL FILCOL(255)
        CALL LINCOL(255)

        IF (iopt_ghost.EQ.0) THEN
          DSPOT = 0.016
          SPOT  = 0.817 - (98.0 - REAL(i1)) / 4.0 * DSPOT
          CALL BOX (0.98-DSPOT*1.0,0.98+DSPOT*1.0,
     .              SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
          miny = MIN(miny,SPOT-DSPOT/2.0)
          maxy = MAX(maxy,SPOT+DSPOT/2.0)
        ELSE
c          DSPOT = 0.0115
          DSPOT = (MAP2Y - MAP1Y) * 0.0373
c          DSPOT = (MAP2Y - MAP1Y) * 0.0383
          SPOT  = MAP2Y - DSPOT - (100.0 - REAL(i1)) / 4.0 * DSPOT
          CALL BOX (MAP2X+0.02*(MAP2X-MAP1X),
     .              MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT,
     .              SPOT+0.5*DSPOT,SPOT+DSPOT)
          miny = MIN(miny,SPOT+0.5*DSPOT)
          maxy = MAX(maxy,SPOT+DSPOT)
c          CALL BOX (MAP2X+0.02-DSPOT*1.0,MAP2X+0.02+DSPOT*1.0,
c     .              SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
c          miny = MIN(miny,SPOT-DSPOT/2.0)
c          maxy = MAX(maxy,SPOT+DSPOT/2.0)
        ENDIF


      ENDDO

c...  Box:
      CALL LINCOL(1)
      IF (iopt_ghost.EQ.0) THEN
        CALL POSITN (0.98-DSPOT*1.0,miny)
        CALL JOIN   (0.98-DSPOT*1.0,maxy)
        CALL POSITN (0.98-DSPOT*1.0,maxy)
        CALL JOIN   (0.98+DSPOT*1.0,maxy)
        CALL POSITN (0.98+DSPOT*1.0,maxy)
        CALL JOIN   (0.98+DSPOT*1.0,miny)
        CALL POSITN (0.98+DSPOT*1.0,miny)
        CALL JOIN   (0.98-DSPOT*1.0,miny)
      ELSE

        CALL POSITN (MAP2X+0.02*(MAP2X-MAP1X),miny)
        CALL JOIN   (MAP2X+0.02*(MAP2X-MAP1X),maxy)
        CALL POSITN (MAP2X+0.02*(MAP2X-MAP1X),maxy)
        CALL JOIN   (MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT,maxy)
        CALL POSITN (MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT,maxy)
        CALL JOIN   (MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT,miny)
        CALL POSITN (MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT,miny)
        CALL JOIN   (MAP2X+0.02*(MAP2X-MAP1X),miny)
c        CALL POSITN (MAP2X+0.02-DSPOT*1.0,miny)
c        CALL JOIN   (MAP2X+0.02-DSPOT*1.0,maxy)
c        CALL POSITN (MAP2X+0.02-DSPOT*1.0,maxy)
c        CALL JOIN   (MAP2X+0.02+DSPOT*1.0,maxy)
c        CALL POSITN (MAP2X+0.02+DSPOT*1.0,maxy)
c        CALL JOIN   (MAP2X+0.02+DSPOT*1.0,miny)
c        CALL POSITN (MAP2X+0.02+DSPOT*1.0,miny)
c        CALL JOIN   (MAP2X+0.02-DSPOT*1.0,miny)
      ENDIF


      CALL LINCOL(1)

      IF (hardscale) THEN
        stepsize = -25.0
      ELSE
        stepsize = -20.0
      ENDIF

      DO i1r = 100.0, 0.0, stepsize

        IF (iopt_ghost.EQ.0) THEN
          DSPOT = 0.016
          SPOT  = 0.820 - (100.0 - REAL(i1r)) / 4.0 * DSPOT
          IF (hardscale) THEN 
            WRITE(label,'(1P,E10.2,0P,A,I3,A)') 
     .        qmin + (qmax - qmin) * i1r/100.0 ,' (',
     .        NINT(1.0+(100.0-1.0)*i1r/100.0),'%)'
          ELSE
            WRITE(label,'(1P,E10.2,0P,A,I3,A)') 
     .        qmin + (qmax - qmin) * i1r/100.0 ,' (',
     .        NINT(1.0+(100.0-1.0)*i1r/100.0),'%)'
          ENDIF
          CALL PLOTST(0.99,SPOT-0.004,label(1:LEN_TRIM(label)))       
        ELSE
          CALL CTRMAG(12)
c          DSPOT = 0.0115
          DSPOT = (MAP2Y - MAP1Y) * 0.0373
c          DSPOT = (MAP2Y - MAP1Y) * 0.0383
          SPOT  = MAP2Y - DSPOT - (100.0 - REAL(i1r)) / 4.0 * DSPOT
          IF (qmax.LT.100.0) THEN
            IF (hardscale) THEN 
              WRITE(label,'(F4.1)') 
     .          qmin + (qmax - qmin) * i1r/100.0
            ELSE
              WRITE(label,'(F4.1)') 
     .          qmin + (qmax - qmin) * i1r/100.0
            ENDIF
            DIST = 0.045*(MAP2X-MAP1X)+2.0*DSPOT
c            DIST = 0.04*(MAP2X-MAP1X)+2.0*DSPOT
            CALL PLOTST(MAP2X+DIST,SPOT,label(1:LEN_TRIM(label)))
          ELSE							 
            IF (hardscale) THEN 				 
              WRITE(label,'(1P,E10.2,0P,A,I3,A)') 		 
     .          qmin + (qmax - qmin) * i1r/100.0 ,' (',		 
     .          NINT(1.0+(100.0-1.0)*i1r/100.0),'%)'		 
            ELSE						 
              WRITE(label,'(1P,E10.2,0P,A,I3,A)') 		 
     .          qmin + (qmax - qmin) * i1r/100.0 ,' (',		 
     .          NINT(1.0+(100.0-1.0)*i1r/100.0),'%)'		 
            ENDIF						 
            CALL PLOTST(MAP2X+0.03,SPOT,label(1:LEN_TRIM(label)))
          ENDIF
        ENDIF
      ENDDO

c...  Draw contour plot label:
      CALL CTRORI(90.0)
      CALL PLOTST(MAP2X+0.02*(MAP2X-MAP1X)+2.0*DSPOT+0.07,miny,
     .            elabs(1)(1:LEN_TRIM(elabs(1))))
      CALL CTRORI(0.0)





210   CONTINUE




c...  Print comments:
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL CTRMAG (12)
      DO i = 1, 10
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.590+(i-1)*0.02,char(i))
      ENDDO
      DO i = 20, 30
        j = i - 19
        IF (char(i).NE.' ') CALL PLOTST(1.00,0.550-(j-1)*0.02,char(i))
      ENDDO


c...  Add a comment to the plot:
79    READ(5,'(A5000)') dummy
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
        GOTO 79
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

c...  Plot grid outline:
      IF (xval.NE.-99.0) THEN

      ELSE
        iplots = 0
        CALL FULL
        CALL LINCOL(1)
        CALL SUPIMP('PARTIAL')      
      ENDIF

c *TEMP*
      READ(5,'(A256)') dummy
      IF   (dummy(8:13).EQ.'Addsur'.OR.dummy(8:13).EQ.'addsur'.OR.
     .      dummy(8:13).EQ.'ADDSUR') THEN
        CALL DrawAdditionalSurfaces(6)
      ELSE
        CALL FRAME
        BACKSPACE 5
      ENDIF

      READ(5,'(A256)') dummy
      IF   (dummy(8:14).EQ.'Noframe'.OR.dummy(8:14).EQ.'noframe'.OR.
     .      dummy(8:14).EQ.'NOFRAME') THEN
      ELSE
        CALL FRAME
        BACKSPACE 5
      ENDIF


c...  Write plot data to the OUT data file:
c      WRITE(6,*) '982: DATA'
c      DO i1 = 1, nc
c        WRITE(6,'(I6,1P,E12.4,0P)') i1,cq(i1)
c      ENDDO

      zval = -99.0

      hardscale = .FALSE.
      colscale = .FALSE.

      IF (optflow.NE.0.AND.
     .    ALLOCATED(tmppinion).AND.
     .    ALLOCATED(tmposmcfp)) THEN
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


c...  Free arrays:
      DEALLOCATE(osmtmp)
      DEALLOCATE(nv)
      DEALLOCATE(rq)
      DEALLOCATE(rv)
      DEALLOCATE(cq)
      DEALLOCATE(zv)
      DEALLOCATE(rc)
      DEALLOCATE(zc)

      IF (MARsum1.NE.0.0) WRITE(0,*) 'MARSUM SOL 982:',MARsum1
      IF (MARsum2.NE.0.0) WRITE(0,*) 'MARSUM PFZ 982:',MARsum2

c     IPP/09 - Krieger - should reset colors here
      call setup_col(ncntr,icntr)

      RETURN
 9012 FORMAT(1X,'PLOT',I3,4X,A)
 98   CALL ER('PLOT982','Problem reading external datafile',*99)
 99   WRITE(0,*) '>'//fname(1:LEN_TRIM(fname))//'<'
      STOP
      END



c
c ======================================================================
c
c
c

