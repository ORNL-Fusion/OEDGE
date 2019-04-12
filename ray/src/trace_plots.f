c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE Output985(iopt,MAXPIXEL,npixel,pixel,image)
      USE mod_out985
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1000,1000)

c...  Initialize GHOST:
      CALL GPSTOP (100)
      CALL PAPER  (1)
      CALL PrinterInit

      CALL WireFramePlot(MAXPIXEL,npixel,pixel,image)       

      CALL DetectorPlot(MAXPIXEL,npixel,pixel,image)       


c      IF (.TRUE.) THEN
c....   Solid:
c        WRITE(0,*) 'DRAWING SOLID 3D PLOT'
c        CALL FRAME
c        CALL DrawSolidPlot(opt,nobj,obj)
c        WRITE(0,*) 'DRAWING FRAME'
c        CALL DrawFrame
c        WRITE(0,*) 'DONE'
c      ENDIF


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: DetectorPlots
c
      SUBROUTINE DetectorPlot(MAXPIXEL,npixel,pixel,image)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_plots
      IMPLICIT none

      INTEGER MAXPIXEL,npixel         ! Put this into _variables... 
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1000,1000)

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER, PARAMETER :: MAXIZS = 10

      CHARACTER glabel*512,caption*1024
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
      CHARACTER*36 XLAB
      CHARACTER*72 YLAB
      CHARACTER*72 SMOOTH
      CHARACTER*128 cdum1

      INTEGER axis1(3)
      LOGICAL axisdata
      REAL    angle1(3)

      INTEGER CH1

      INTEGER   nxbin,nybin,ipixel,nbin,isrf,isrf1,isrf2,isid,
     .          idet,iplot,iobj,i1,i2,i3,i4,ix,iy
      REAL      deltar,deltaz,qmin,qmax,qval,ang1,ang2,ang,
     .          frac1,frac2,xcen,ycen,xnear,ynear,count
      REAL      dangle
      REAL*8    angle,p1(3,6),p2(3,6),x1,x2,z1,mat(3,3)

      INTEGER dat1,dat2
      REAL, ALLOCATABLE :: xdat(:),ydat(:)
      CHARACTER xlabel*256,ylabel*256

      INTEGER, PARAMETER :: ncols = 1
      REAL, PARAMETER :: PI = 3.141596
      REAL cxmin,cxmax,cymin,cymax

c *TEMP*
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)





c...  PLOT B: 3D LOS integration through the 3D objects
c     -build list of primary chords
c     -build list of secondary chords (reflections)
c     -fast routine for determining the intersection between a line and a surface
c     -ray trace (using connection map and wedge index to speed things up)



      IF (npixel.GT.1) THEN
c...    Image:

        DO idet = 1, opt%ndet

          CALL FRAME

          nxbin = 0
          nybin = 0
          qmin =  0.0
c          qmin =  HI
          qmax = -1.0E+30
          DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
c          DO ipixel = 1, npixel
            ix = pixel(ipixel)%xindex
            iy = pixel(ipixel)%yindex
            nxbin = MAX(nxbin,ix)
            nybin = MAX(nybin,iy)
c            qmin = MIN(qmin,pixel(ipixel)%integral)
            qmax = MAX(qmax,SNGL(pixel(ipixel)%integral(1)))
          ENDDO
c          qmin = MIN(qmin,0.0)

          WRITE(0,*) 'QMAD:',qmin,qmax

          IF (qmin.EQ.qmax) THEN
            WRITE(0,*) '...PROBLEM, CYCLING'
            CYCLE
          ENDIF

          IF (nxbin.GE.nybin) THEN
            map1x = 0.70          
            map2x = map1x + 0.60
            map2y = 0.70
            map1y = map2y - 0.60 * REAL(nybin) / REAL(nxbin)
          ELSE
            map1x = 0.70          
            map2x = map1x + 0.60 * REAL(nxbin) / REAL(nybin)
            map1y = 0.15 
            map2y = map1y + 0.60 
          ENDIF 

          CALL PSPACE(map1x,map2x,map1y,map2y)
          CALL MAP   (0.0,1.0,1.0,0.0)

          ALLOCATE(nv(1))
          ALLOCATE(rv(4,1))  
          ALLOCATE(zv(4,1))
          ALLOCATE(cq(1))

          IF (.NOT..TRUE.) THEN
            DO frac1 = 0.0, 1.0*0.99999, 0.01                  ! Needs more work, grouping regons of common colour. 
              CALL SetCol255_04(2,frac1+0.005,0.0,1.0)         ! Otherwise, not much savings... 
              CALL FILCOL(255)                         
              CALL LINCOL(255)                         
              DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
c              DO ipixel = 1, npixel
                cq(1) = SNGL(pixel(ipixel)%integral(1))
                frac2 = (cq(1) - qmin) / (qmax - qmin)
                IF (frac2.LT.frac1.OR.frac2.GE.frac1+0.01) CYCLE

                ix = pixel(ipixel)%xindex
                iy = pixel(ipixel)%yindex
                deltar = 1.0 / REAL(nxbin)
                deltaz = 1.0 / REAL(nybin)
                i1 = 1
                rv(1,i1) = (ix - 1) * deltar 
                zv(1,i1) = (iy - 1) * deltaz
                rv(2,i1) = (ix - 1) * deltar 
                zv(2,i1) = (iy    ) * deltaz
                rv(3,i1) = (ix    ) * deltar 
                zv(3,i1) = (iy    ) * deltaz
                rv(4,i1) = (ix    ) * deltar 
                zv(4,i1) = (iy - 1) * deltaz
                CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
              ENDDO
            ENDDO
          ELSE
c            DO ipixel = 1, npixel
            DO ipixel = opt%det_istart(idet), opt%det_iend(idet)
              ix = pixel(ipixel)%xindex
              iy = pixel(ipixel)%yindex
              deltar = 1.0 / REAL(nxbin)
              deltaz = 1.0 / REAL(nybin)
              i1 = 1
              rv(1,i1) = (ix - 1) * deltar 
              zv(1,i1) = (iy - 1) * deltaz
              rv(2,i1) = (ix - 1) * deltar 
              zv(2,i1) = (iy    ) * deltaz
              rv(3,i1) = (ix    ) * deltar 
              zv(3,i1) = (iy    ) * deltaz
              rv(4,i1) = (ix    ) * deltar 
              zv(4,i1) = (iy - 1) * deltaz
              cq(1) = SNGL(pixel(ipixel)%integral(1))
              CALL SetCol255_04(2,cq(i1),qmin,qmax)    ! Should really reorganize this loop, to draw all pixels of a 
              CALL FILCOL(255)                         ! particular shade, from a limited set of contours, to 
              CALL LINCOL(255)                         ! try and shrink the size of the .ps file... 
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
            ENDDO
          ENDIF
c...      Clear arrays:
          DEALLOCATE(nv)
          DEALLOCATE(rv)
          DEALLOCATE(zv)
          DEALLOCATE(cq)
          CALL DrawColourScale(1,2,qmin,qmax,'none')
c...      Frame the plot:
          CALL DrawFrame

c...      Line plot:
c          slopt2 = 1
c          iopt_ghost = 1  ! 2
c          plottype(1) = 2
c          plottype(2) = 3
c          map1x = 0.65             
c          map2x = map1x + 0.55
c          map1y = 0.77 
c          map2y = map1y + 0.15
cc...      Assign data:
c          dat1 = 35600+1
c          dat2 = 35600+200
c          ALLOCATE(xdat(dat2-dat1+1))
c          ALLOCATE(ydat(dat2-dat1+1))
c          xlabel = 'pixel    '
c          ylabel = 'signal   '
c          DO i1 = dat1, dat2
c            xdat(i1-dat1+1) = REAL(i1)
c          ENDDO
c          ydat(1:dat2-dat1+1) = SNGL(pixel(dat1:dat2)%integral(1))
c          IF (idet.EQ.1) THEN
c            WRITE(6,*) '*PIXEL DATA'
c            WRITE(6,*) '*'
c            WRITE(6,*) dat2-dat1+1
c            DO i1 = 1, dat2-dat1+1
c              WRITE(6,*) REAL(i1),xdat(i1),ydat(i1)
c            ENDDO
c          ENDIF
c          cxmin = 1.0
c          cxmax = REAL(dat2-dat1+1)
c          cymin =  HI
c          cymax = -HI
c          DO i1 = 1, dat2-dat1+1
c            cymin = MIN(cymin,ydat(i1))
c            cymax = MAX(cymax,ydat(i1))
c          ENDDO
c          CALL GRTSET_TRIM(' ',' ',' ',' ',' ',      ! TITLE,REF,nVIEW,PLANE,glabel,
c     .                     cxmin,cxmax,cymin,cymax,
c     .                     ' ',xlabel,ylabel,        ! TABLE,XLAB,YLAB,
c     .                     0,' ',0,' ',1)            ! 0,smooth,0,ANLY,1)
c          CALL GRTRAC(xdat,ydat,dat2-dat1+1,'ref ','LINE',1)        
c          DEALLOCATE(xdat)
c          DEALLOCATE(ydat)
c          CALL DrawFrame

        ENDDO
      ENDIF







      IF (.TRUE..AND.opt%img_opt.NE.0) THEN
c...    Plot loaded camera image (if there is one):

c        IF (opt%img_nxbin.GE.100) THEN
c          nbin = opt%img_nxbin / 50          ! Should also check opt%nybin
c          nbin = 1
c        ELSE
          nbin = 1
c        ENDIF

c...    Plot image:
c        nbin = 1 !6  ! Binning (NBIN=1 is no binning)
        IF (nbin.GT.1.AND.
     .      (nbin.GE.opt%img_nxbin.OR.
     .       nbin.GE.opt%img_nybin))
     .    CALL ER('985','Bin request greater than image resolution',
     .            *99)

        nxbin = opt%img_nxbin / nbin
        nybin = opt%img_nybin / nbin

        WRITE(0,*) ' *** NXBIN',nxbin,nybin

        qmin =  0.0
        qmax = -1.0E+30
        DO ix = 1, nxbin
          DO iy = 1, nybin
            qval = opt%img_image(ix,iy)
c            qval = 0.0
c            DO i2 = 0, nbin-1
c              DO i3 = 0, nbin-1
c                qval = qval + SNGL(opt%img_image(nbin*(ix-1)+1+i2,
c     .                                           nbin*(iy-1)+1+i3))
c              ENDDO
c            ENDDO           
            qmax = MAX(qmax,qval)  ! Need to average here, and below?  
          ENDDO
        ENDDO

c        WRITE(0,*) 'QVAL:',qmax,nxbin,nybin
c        STOP 'sdfsd'

        IF (nxbin.GE.nybin) THEN
          map1x = 0.05            
          map2x = map1x + 0.60
          map2y = 0.70
          map1y = map2y - 0.60 * REAL(nybin) / REAL(nxbin)
        ELSE
          map1x = 0.05            
          map2x = map1x + 0.60 * REAL(nxbin) / REAL(nybin)
          map1y = 0.15 
          map2y = map1y + 0.60 
        ENDIF 

        CALL PSPACE(map1x,map2x,map1y,map2y)
        CALL MAP   (0.0,1.0,1.0,0.0)
        ALLOCATE(nv(1))
        ALLOCATE(rv(4,1))  
        ALLOCATE(zv(4,1))
        ALLOCATE(cq(1))
        deltar = 1.0 / REAL(nxbin)
        deltaz = 1.0 / REAL(nybin)
        DO ix = 1, nxbin
          DO iy = 1, nybin
            i1 = 1
            rv(1,i1) = REAL(ix - 1) * deltar 
            zv(1,i1) = REAL(iy - 1) * deltaz
            rv(2,i1) = REAL(ix - 1) * deltar 
            zv(2,i1) = REAL(iy    ) * deltaz
            rv(3,i1) = REAL(ix    ) * deltar 
            zv(3,i1) = REAL(iy    ) * deltaz
            rv(4,i1) = REAL(ix    ) * deltar 
            zv(4,i1) = REAL(iy - 1) * deltaz
            cq(i1) = opt%img_image(ix,iy)
c            cq(i1) = 0.0
c            DO i2 = 0, nbin-1
c              DO i3 = 0, nbin-1
c                cq(i1) = cq(i1) + SNGL(opt%img_image(nbin*(ix-1)+1+i2,
c     .                                               nbin*(iy-1)+1+i3))
c              ENDDO
c            ENDDO
            CALL SetCol255_04(2,cq(i1),qmin,qmax)   ! See above note on reducing size of .ps file...
            CALL FILCOL(255)
            CALL LINCOL(255) 
            CALL PTPLOT(rv(1,i1),zv(1,i1),1,4,1)
          ENDDO
        ENDDO
c...    Clear arrays:
        IF (ALLOCATED(nv)) DEALLOCATE(nv)
        IF (ALLOCATED(rv)) DEALLOCATE(rv)
        IF (ALLOCATED(zv)) DEALLOCATE(zv)
        IF (ALLOCATED(cq)) DEALLOCATE(cq)
c...    Annotate:
        CALL LinCol(ncols+1)
        CALL CTRMAG(12)
        WRITE(caption,'(A,2I5)') 'XBIN,YBIN:',nxbin,nybin
        CALL PLOTST(0.02,0.03,caption(1:LEN_TRIM(caption)))
        CALL DrawFrame
c...    Clear memory:
c        IF (ALLOCATED(opt%img_image)) DEALLOCATE(opt%img_image)     
      ENDIF


      IF (.TRUE.) THEN
c...    Spectra:
c        slopt2 = 1
        iopt_ghost = 1
        DO idet = 1, opt%ndet
          WRITE(0,*) 'PLOTTING LINE SHAPES',
     .        opt%det_istart(idet),
     .        opt%det_iend  (idet)
          CALL FRAME
          CALL PlotLineShapes(opt%det_istart(idet),
     .                        opt%det_iend  (idet),
     .                        MAXPIXEL,npixel,pixel)
        ENDDO
        WRITE(0,*) 'DONE'
      ENDIF


      CALL Frame



      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WireFramePlot
c
      SUBROUTINE WireFramePlot(MAXPIXEL,npixel,pixel,image)
      USE mod_out985
      USE mod_out985_variables
      USE mod_out985_plots
      IMPLICIT none

      INTEGER MAXPIXEL,npixel         ! Put this into _variables... 
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1000,1000)

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      INTEGER, PARAMETER :: MAXIZS = 10

      CHARACTER GRAPH1*256
c      CHARACTER glabel*512,caption*1024
c      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256
c      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:MAXIZS+1)
c      CHARACTER*36 XLAB
c      CHARACTER*72 YLAB
c      CHARACTER*72 SMOOTH
      CHARACTER*128 cdum1

      INTEGER axis1(3)
      LOGICAL axisdata
      REAL    angle1(3)

      INTEGER CH1

      INTEGER   nxbin,nybin,ipixel,nbin,isrf,isrf1,isrf2,isid,
     .          idet,iplot,iobj,i1,i2,i3,i4,ix,iy
      REAL      deltar,deltaz,qmin,qmax,qval,ang1,ang2,ang,
     .          frac1,frac2,xcen,ycen,xnear,ynear,count
      REAL      dangle
      REAL*8    angle,p1(3,6),p2(3,6),x1,x2,z1,mat(3,3)

      INTEGER dat1,dat2
      REAL, ALLOCATABLE :: xdat(:),ydat(:)
c      CHARACTER xlabel*256,ylabel*256

      INTEGER, PARAMETER :: ncols = 1
      REAL, PARAMETER :: PI = 3.141596
      REAL cxmin,cxmax,cymin,cymax

c *TEMP*
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)



c...  PLOT A: 3D surface geometry plot for a given X,Y,Z along a particular view
c     -Translate and rotate as necessary:
c     -Project onto the R,Z plane:

c      glabel = job  (CH1(job  ):LEN_TRIM(job  ))//'     '//
c     .         graph(CH1(graph):LEN_TRIM(graph))

      ! jdemod - changed the : to ) to balance the parentheses in the format specifiers
      !          NOTE: these write statements will produce lines of 512 characters each
      !                right-padded with blanks 
c      WRITE(nview ,'(512(A))') (' ',i1=1,LEN(nview )) 
c      WRITE(plane ,'(512(A))') (' ',i1=1,LEN(plane )) 
c      WRITE(anly  ,'(512(A))') (' ',i1=1,LEN(anly  )) 
c      WRITE(table ,'(512(A))') (' ',i1=1,LEN(table )) 
c      WRITE(xlab  ,'(512(A))') (' ',i1=1,LEN(xlab  )) 
c      WRITE(ylab  ,'(512(A))') (' ',i1=1,LEN(ylab  )) 
c      WRITE(smooth,'(512(A))') (' ',i1=1,LEN(smooth)) 
c      WRITE(graph ,'(512(A))') (' ',i1=1,LEN(graph )) 
c      WRITE(job   ,'(512(A))') (' ',i1=1,LEN(job   )) 

c      XLAB = '   R  (M)'
c      YLAB = '   Z  (M)'



c      slopt = 1



c...  Look for zoom data:
      READ(5,'(A80)',END=150) graph1
      IF (graph1(8:11).EQ.'Zoom'.OR.graph1(8:11).EQ.'ZOOM'.OR.
     .    graph1(8:11).EQ.'zoom') THEN
        READ(graph1,*) cdum1,xcen,ycen,xnear,ynear
        cxmin = xcen - xnear
        cxmax = xcen + xnear
        cymin = ycen - ynear
        cymax = ycen + ynear
      ELSE
        cxmin =  0.0
        cxmax =  6.0
        cymin = -3.0
        cymax =  3.0
        BACKSPACE 5
      ENDIF
150   CONTINUE

      axis1 (1) = 1
      axis1 (2) = 2
      axis1 (3) = 3
      angle1(1) = 0.0
      angle1(2) = 0.0
      angle1(3) = 0.0

      axisdata = .FALSE.

      iplot = 0
      DO WHILE (.TRUE.) 
c...    Wireframe:

        READ(5,'(A256)') graph1
        IF (graph1(8:11).EQ.'Axis'.OR.graph1(8:11).EQ.'AXIS'.OR.
     .      graph1(8:11).EQ.'axis') THEN
          READ(graph1,*) cdum1,(axis1(i1),angle1(i1),i1=1,3)

          iplot = iplot + 1
          IF (iplot.EQ.5) THEN
           iplot = 1
           CALL FRAME
          ENDIF
c          IF (axisdata) CALL FRAME
        ELSE
          BACKSPACE 5
          IF (axisdata) EXIT
        ENDIF

        axisdata = .TRUE.  ! Improper use, but handy at the moment...

c...    Setup transformation matrix:
        CALL Calc_Transform2(mat,0.0D0,1,0)
        DO i1 = 1, 3
          angle = DBLE(angle1(i1)*PI/180.0)
          CALL Calc_Transform2(mat,angle,axis1(i1),1)
        ENDDO

c...    Use GRTSET_TRIM:
c        slopt4 = 1
c...    Stopping resizing of scale font in ghost1.o6a:
        iopt_ghost = 1

        SELECTCASE (iplot)
          CASE(1)
            map1x = 0.07
            map1y = 0.55
          CASE(2)
            map1x = 0.57
            map1y = 0.55
          CASE(3)
            map1x = 0.07
            map1y = 0.08
          CASE(4)
            map1x = 0.57
            map1y = 0.08
          CASEDEFAULT
        ENDSELECT

        map2x = map1x + 0.40
        map2y = map1y + 0.40

c...    Draw polygons:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)

        count = 0.0

        DO iobj = 1, nobj

c          IF (iobj.NE.1035) CYCLE  ! *TEMP*
c          IF (iobj.NE.471) CYCLE  ! *TEMP*
c          IF (iobj.NE.471.AND.iobj.NE.472.AND.iobj.NE.458) CYCLE  ! *TEMP*


          CALL LINCOL(ncols+obj(iobj)%colour) 

          DO isid = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)

c             IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE

c            IF (isid.NE.2) CYCLE

c            IF (obj(iobj)%flag(isid).NE.-1.AND.
c     .          obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE  ! *TEMP*

            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL.AND.
     .          obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY.AND. 
     .          obj(iobj)%tsur(isid).NE.SP_GRID_SURFACE) CYCLE

c              WRITE(6,*) 'BOUNDARY?',obj(iobj)%ik,obj(iobj)%ir

c            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE   ! Recent filter...

c            IF (obj(iobj)%imap(1,isid).NE.0) CYCLE
c            IF (obj(iobj)%ik.NE.16) CYCLE


            IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN
              count = count + 
     .                1.0 / REAL(MAX(obj(iobj)%nsur,obj(iobj)%nside))
c             IF (count.GT.100.0) CYCLE
              IF (count.GT.15000.0) CYCLE
            ENDIF




c *** NEED TO STORE LINE SEGMENTS AND DELETE DUPLICATES... 


            IF (obj(iobj)%nside.NE.0) THEN
              isrf1 = obj(iobj)%iside(isid,1)
              isrf2 = obj(iobj)%iside(isid,2)
            ELSE
              isrf1 = 1
              isrf2 = 1
            ENDIF


            IF     (obj(iobj)%gsur(isid).EQ.GT_TC) THEN

c              count = count + 20

              dangle =  15.0 * PI / 180.0

c              IF (.FALSE.) THEN
              IF (.TRUE..OR.iplot.EQ.1) THEN
                ang1 = 0.0
                ang2 = 1.0 * PI / 180.0
              ELSE
                ang1 = 0.0
                ang2 = 359.0 * PI / 180.0
              ENDIF

c              DO ang = 0.0, 1.0 * PI / 180.0, dangle
              DO ang = ang1, ang2, dangle

c                IF (ang.GT.0.0.AND.
c     .              obj(iobj)%type.NE.OP_INTEGRATION_VOLUME) CYCLE

c                DO isur = 1, 1  ! THIS IS HERE FOR WHEN SIDES ARE USED!

                DO isrf = isrf1, isrf2

                  IF (obj(iobj)%nside.NE.0) THEN
                    p1(1,1) = vtx(1,srf(isrf)%ivtx(1))
                    p1(2,1) = vtx(2,srf(isrf)%ivtx(1))
                    p1(3,1) = p1(1,1) * DTAN(-0.5D0*dangle)
                    p2(1,1) = p1(1,1)
                    p2(2,1) = p1(2,1)
                    p2(3,1) = p2(1,1) * DTAN(+0.5D0*dangle)

                    p1(1,2) = vtx(1,srf(isrf)%ivtx(2))
                    p1(2,2) = vtx(2,srf(isrf)%ivtx(2))
                    p1(3,2) = p1(1,2) * DTAN(-0.5D0*dangle)
                    p2(1,2) = p1(1,2)
                    p2(2,2) = p1(2,2)
                    p2(3,2) = p2(1,2) * DTAN(+0.5D0*dangle)
                  ELSE
                    p1(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(1,isid)) 
                    p1(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(1,isid)) 
                    p1(3,1) = DBLE(p1(1,1))*DTAN(DBLE(-0.5*dangle))
                    p2(1,1) = p1(1,1)
                    p2(2,1) = p1(2,1)
                    p2(3,1) = DBLE(p2(1,1))*DTAN(DBLE(+0.5*dangle))

                    p1(1,2) = obj(iobj)%v(1,obj(iobj)%ipts(2,isid)) 
                    p1(2,2) = obj(iobj)%v(2,obj(iobj)%ipts(2,isid)) 
                    p1(3,2) = DBLE(p1(1,2))*DTAN(DBLE(-0.5*dangle))
                    p2(1,2) = p1(1,2)
                    p2(2,2) = p1(2,2)
                    p2(3,2) = DBLE(p2(1,2))*DTAN(DBLE(+0.5*dangle))
                  ENDIF
c...              Rotate vertices:
                  DO i3 = 1, 2
                    x1 = p1(1,i3)
                    z1 = p1(3,i3)
                    p1(1,i3) = DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang))*z1
                    p1(3,i3) = DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang))*z1
                    x1 = p2(1,i3)
                    z1 = p2(3,i3)
                    p2(1,i3) = DCOS(DBLE(ang)) * x1 - DSIN(DBLE(ang))*z1
                    p2(3,i3) = DSIN(DBLE(ang)) * x1 + DCOS(DBLE(ang))*z1
                    call transform_vect(mat,p1(1,i3))
                    call transform_vect(mat,p2(1,i3))
                  ENDDO

                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p1(1,2)),SNGL(p1(2,2))) 
                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                  CALL POSITN (SNGL(p1(1,2)),SNGL(p1(2,2)))
                  CALL JOIN   (SNGL(p2(1,2)),SNGL(p2(2,2))) 
                  CALL POSITN (SNGL(p2(1,1)),SNGL(p2(2,1)))
                  CALL JOIN   (SNGL(p2(1,2)),SNGL(p2(2,2))) 
                ENDDO
              ENDDO

            ELSEIF (obj(iobj)%gsur(isid).EQ.GT_TD) THEN

              IF (obj(iobj)%nside.NE.0) THEN
                DO isrf = isrf1, isrf2
                  DO i3 = 1, srf(isrf)%nvtx
                    i4 = i3 + 1
                    IF (i4.GT.srf(isrf)%nvtx) i4 = 1
                    p1(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i3))
                    p2(1:3,1) = vtx(1:3,srf(isrf)%ivtx(i4))
                    call transform_vect(mat,p1(1,1))
                    call transform_vect(mat,p2(1,1))
                    CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                    CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                  ENDDO
                ENDDO
              ELSE
                DO i3 = 1, obj(iobj)%npts(isid)
                  i4 = i3 + 1
                  IF (i4.EQ.obj(iobj)%npts(isid)+1) i4 = 1
                  p1(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(i3,isid))
                  p1(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(i3,isid))
                  p1(3,1) = obj(iobj)%v(3,obj(iobj)%ipts(i3,isid))
                  p2(1,1) = obj(iobj)%v(1,obj(iobj)%ipts(i4,isid))
                  p2(2,1) = obj(iobj)%v(2,obj(iobj)%ipts(i4,isid))
                  p2(3,1) = obj(iobj)%v(3,obj(iobj)%ipts(i4,isid))
                  call transform_vect(mat,p1(1,1))
                  call transform_vect(mat,p2(1,1))
                  CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
                  CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
                ENDDO
              ENDIF
            ELSE
              CALL ER('Plot985','Unknown surface geometry type',*98)
            ENDIF

          ENDDO
        ENDDO 
c...    Draw pixel views:
        IF (.FALSE.) THEN
          CALL LINCOL(ncols+55) 
          DO i1 = 1, npixel
            IF (pixel(i1)%yindex.NE.1) CYCLE
            p1(1,1) = pixel(i1)%v1(1)
            p1(2,1) = pixel(i1)%v1(2)
            p1(3,1) = pixel(i1)%v1(3)
            p2(1,1) = pixel(i1)%v2(1)
            p2(2,1) = pixel(i1)%v2(2)
            p2(3,1) = pixel(i1)%v2(3)

            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))

            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ELSEIF (.TRUE.) THEN
c        ELSEIF (.TRUE..AND.iplot.GT.1) THEN
          CALL LINCOL(ncols+55) 
c          DO i1 = 1, MIN(1,nchord)
          DO i1 = 1, MIN(500,nchord)
            p1(1,1) = s_chord(i1)%v1(1)
            p1(2,1) = s_chord(i1)%v1(2)
            p1(3,1) = s_chord(i1)%v1(3)
            p2(1,1) = s_chord(i1)%v2(1)
            p2(2,1) = s_chord(i1)%v2(2)
            p2(3,1) = s_chord(i1)%v2(3)
            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))
            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ENDIF
c...    Frame:
        CALL LINCOL(1)
        CALL DrawFrame
      ENDDO
 
      CALL Frame

      RETURN
98    WRITE(0,*) 'OBJECT,SIDE,TYPE=',i1,i2,obj(i1)%gsur(i2)
99    STOP
      END
c
c ======================================================================
c
      SUBROUTINE DrawFrame
      USE mod_out985_plots
      IMPLICIT none
c...  Finish plot:
      CALL LINCOL(1)
      CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
      CALL MAP   (0.0,1.0,0.0,1.0)
      CALL POSITN(0.0,0.0)
      CALL JOIN  (0.0,1.0)              
      CALL POSITN(0.0,1.0)
      CALL JOIN  (1.0,1.0)              
      CALL POSITN(1.0,1.0)
      CALL JOIN  (1.0,0.0)              
      CALL POSITN(1.0,0.0)
      CALL JOIN  (0.0,0.0)                      
      RETURN
      STOP
      END
c
c ======================================================================
c
c subroutine: SetCol255
c
c
      SUBROUTINE SetCol255_04(mode,qval,qmin,qmax)
      IMPLICIT none

      INTEGER mode 
      REAL    qval,qmin,qmax
    

c * OLD *
      REAL frac1,hue,pastel,midfrac,scale
      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale
      LOGICAL grayscale

      INTEGER last_mode
      REAL bright
      REAL frac,frac5,fmod5,last_frac,last_bright   ! mode.eq.1

      DATA last_mode,last_frac,last_bright /-1, -1.0, -1.0/    

      SAVE


      IF     (mode.EQ.1) THEN

        frac = (qval - qmin) / (qmax - qmin)
        frac5 = 100.0*frac
        fmod5 = AMOD(frac5,2.0)
        frac = MIN(0.98,(frac5-fmod5)/100.0)

        bright = 1.0-(0.98-frac)**20

        IF (mode.NE.last_mode.OR.frac.NE.last_frac.OR.
     .      bright.NE.last_bright) 
     .    CALL ColSet(0.75*frac+0.25,1.0,bright,255)

c      ELSEIF (mode.EQ.2) THEN

      ELSEIF (mode.EQ.2) THEN
c...    Grayscale:

        frac = (qval - qmin) / (qmax - qmin)

c        frac = frac * 0.9 + 0.1

        IF (mode.NE.last_mode.OR.frac.NE.last_frac) 
c     .    CALL ColSet(0.0,0.0,frac,255)  ! Default
c     .    CALL ColSet(0.50,1.0-frac,frac, ! Original DIVCAM
c     .                255)
     .    CALL ColSet(0.30+0.230*frac,
     .                MIN(1.0,4.0*(1.0-frac)    ),
     .                MIN(1.0,4.0*     frac     ),
     .                255)
c     .    CALL ColSet(0.25+1.000*frac,
c     .                MIN(1.0,4.0*(1.0-frac)    ),
c     .                MIN(1.0,4.0*     frac     ),
c     .                255)
c     .    CALL ColSet(0.50,1.0-2.0*(MAX(0.0,frac**0.5-0.5)),frac**0.5,
c     .                255)

c        CALL ColSet(0.0,0.0,frac,255)
c        CALL ColSet(0.0,0.0,1.0-frac,255)


      ELSE
        CALL ER('SetCol255_04','Invalid mode',*99)
      ENDIF


      last_mode = mode
      last_frac = frac
      last_bright = bright

      RETURN
c
c * OLD *
c

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

        IF (.NOT..TRUE.) THEN
          IF (frac.LT.0.01) THEN
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
c subroutine: DrawColourScale
c
c
      SUBROUTINE DrawColourScale(mode,colmode,qmin,qmax,label)
      USE mod_out985_plots
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'slout'

      INTEGER   mode,colmode 
      REAL      qmin,qmax
      CHARACTER label*(*)

      INTEGER CH1

      REAL      qval,dspot,spot,dist,minx,maxx,miny,maxy,stepsize,i1r,
     .          rscale,dscale,xmove,ymove,xpos,ypos,nchar
      LOGICAL   yshift
      CHARACTER nums*256      
 
      CALL PSPACE (0.0, 1.35, 0.0,1.0)
      CALL CSPACE (0.0, 1.35, 0.0,1.0)
      CALL MAP    (0.0, 1.35, 0.0,1.0)
      CALL FULL
      CALL HSV  ! This should be moved...

      IF     (mode.EQ.1) THEN  ! Vertical scale drawn to the left of the plot
        dspot = 0.016         
        minx = map2x + 0.02
        maxx = map2x + 0.04
        miny = 1.0E+30
        maxy = map2y - dspot/2.0
        dscale = 2.0
        DO rscale = 100.0-dscale, 0.0, -dscale
          qval = (rscale + 0.5 * dscale) / 100.0 * (qmax - qmin) + qmin
          CALL SetCol255_04(colmode,qval,qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255)
          SPOT  = maxy - (100.0 - dscale - rscale) / 4.0 * DSPOT
          CALL BOX (minx,maxx,SPOT-DSPOT/2.0,SPOT+DSPOT/2.0)
          miny = MIN(miny,SPOT-DSPOT/2.0)
          maxy = MAX(maxy,SPOT+DSPOT/2.0)
        ENDDO
c...    Box:
        CALL LINCOL(1)
        CALL POSITN (minx,miny)
        CALL JOIN   (minx,maxy)
        CALL POSITN (minx,maxy)
        CALL JOIN   (maxx,maxy)
        CALL POSITN (maxx,maxy)
        CALL JOIN   (maxx,miny)
        CALL POSITN (maxx,miny)
        CALL JOIN   (minx,miny)
c...    Text:
        CALL CTRMAG(12)
        dscale = 20.0
        DO rscale = 100.0, 0.0, -dscale
          qval = qmin + (qmax - qmin) * rscale / 100.0
          IF (qmax.GT.1.0.AND.qmax.LT.100.0) THEN
            WRITE(nums,'(F4.1)') 
     .        qval
c            WRITE(nums,'(F4.1,A,I3,A)') 
c     .        qval,' (',NINT(qval/qmax*100.0),'%)'
          ELSE
            WRITE(nums,'(1P,E10.2,0P,A,I3,A)') 
     .        qval,' (',NINT(qval/qmax*100.0),'%)'
          ENDIF
          spot = rscale / 100.0 * (maxy - miny) + miny
          IF (label.NE.'none') THEN
c...        Tick:
            CALL POSITN(maxx      ,spot)
            CALL JOIN  (maxx+0.005,spot)
c...        Label:
            CALL PLOTST(maxx+0.015,spot,nums(1:LEN_TRIM(nums)))
c            CALL PLOTST(maxx+0.005,spot,nums(1:LEN_TRIM(nums)))
          ENDIF
        ENDDO
        IF (label.NE.'none') THEN
          CALL CTRORI(90.0)
          CALL PLOTST(maxx+0.070,0.5*(miny+maxy),
     .                label(1:LEN_TRIM(label)))
          CALL CTRORI(0.0)
        ENDIF
c        CALL PLOTST(maxx+0.010,0.5*(miny+maxy),'T')


      ELSEIF (mode.EQ.2) THEN  ! Horizontal scale below plot
        dspot = 0.010 ! for dist = 0.70 987 plots                
c        dspot = 0.012 ! for dist = 0.80 987 plots                
c        dspot = 0.0134 
c        dspot = 0.016
        IF (label.EQ.'none') THEN
          minx = 1.0E+30           ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.030     ! map1y - 0.04      ! HI
          maxy = map1y - 0.010     ! map1y - 0.02      ! map2y - dspot/2.0
        ELSE
          minx = 1.0E+30           ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.075     ! map1y - 0.04      ! HI
          maxy = map1y - 0.055     ! map1y - 0.02      ! map2y - dspot/2.0
        ENDIF
        dscale = 2.0
        DO rscale = 100.0-dscale, 0.0, -dscale
          qval = (rscale + 0.5 * dscale) / 100.0 * (qmax - qmin) + qmin
          CALL SetCol255_04(colmode,qval,qmin,qmax)
          CALL FILCOL(255)
          CALL LINCOL(255)
          SPOT  = maxx - (100.0 - dscale - rscale) / 4.0 * DSPOT
          CALL BOX (spot-dspot/2.0,spot+dspot/2.0,miny,maxy)
          minx = MIN(minx,SPOT-DSPOT/2.0)
          maxx = MAX(maxx,SPOT+DSPOT/2.0)
        ENDDO
c...    Box:
        CALL LINCOL(1)
        CALL POSITN (minx,miny)
        CALL JOIN   (minx,maxy)
        CALL POSITN (minx,maxy)
        CALL JOIN   (maxx,maxy)
        CALL POSITN (maxx,maxy)
        CALL JOIN   (maxx,miny)
        CALL POSITN (maxx,miny)
        CALL JOIN   (minx,miny)
c...    Text:
        CALL CTRMAG(12)
        dscale = 20.0
        yshift = .FALSE.
        ymove = 0.0
        DO rscale = 100.0, 0.0, -dscale
          qval = qmin + (qmax - qmin) * rscale / 100.0
c          qval = qval / qmax
          IF     (qmax.GT.0.1.AND.qmax.LE.1.0) THEN
            WRITE(nums,'(F3.1)') qval
          ELSEIF ((qmax.GT. 1.0.AND.qmax.LE. 999.9).OR.
     .            (qmax.LE.-1.0.AND.qmax.GT.-999.9)) THEN
            WRITE(nums,'(F5.1)') qval
          ELSE
            WRITE(nums,'(1P,E9.1)') qval
            yshift = .TRUE.
          ENDIF
          spot = rscale / 100.0 * (maxx - minx) + minx

          IF (label.NE.'none') THEN
c...        Tick:
            CALL POSITN(spot,miny      )
            CALL JOIN  (spot,miny-0.005)
c...        Label:
            ypos = miny - 0.020 + ymove
            CALL PCSCEN(spot,ypos,nums(CH1(nums):LEN_TRIM(nums)))
            IF (yshift.AND.ymove.EQ.0.0) THEN
              ymove = -0.02
            ELSE
              ymove = 0.0
            ENDIF       
          ENDIF
        ENDDO

        IF (label.NE.'none') THEN
          IF (yshift) ymove = -0.02
          IF (.TRUE.) THEN
            nchar = REAL(LEN_TRIM(label) - CH1(label) + 1) 
            xpos = 0.5 * (minx + maxx) - 0.5 * nchar * 0.0090
          ELSE
            xpos = 0.5 * (minx + maxx)  !need to not count formatting characters...
          ENDIF
          ypos = miny - 0.045 + ymove 
          CALL PLOTST(xpos,ypos,label(CH1(label):LEN_TRIM(label)))
        ENDIF

      ELSE
        CALL ER('DrawColourScale','Invalid mode',*99)
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
