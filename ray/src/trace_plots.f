c
c ======================================================================
c
c subroutine: PlotLineShapes
c
c      SUBROUTINE PlotLineShapes(istart,iend,MAXPIXEL,npixel,pixel)
c      USE mod_out985
c      USE mod_out985_variables
c      IMPLICIT none

c      INTEGER istart,iend,MAXPIXEL,npixel
c      TYPE(type_view) :: pixel(MAXPIXEL)

c      RETURN
c 9    STOP
c      END
c
c ======================================================================
c
c subroutine: GetSchematics
c
c      SUBROUTINE GetSchematics(xin,yin,zin,
c     .                         mode,MAXSURFACE,MAXPOINTS,icolour,
c     .                         opt,nobj,obj,
c     .                         nsur,npts,hsur,vsur,len1,len2)
c      IMPLICIT none

c...  Input:
c      INTEGER mode,MAXSURFACE,MAXPOINTS,nobj,icolour,
c     .        nsur,npts(MAXSURFACE),hsur(MAXSURFACE)
c      REAL    xin,yin,zin,rdum1,rdum2,rdum3,rdum4
c      REAL*8  vsur(3,MAXPOINTS,0:MAXSURFACE),len1,len2
c      CHARACTER buffer*1024,file*1024
c      TYPE(type_options985) :: opt
c      TYPE(type_3D_object)  :: obj(nobj)
c
c      RETURN
c 98   WRITE(0,*) ' FILE = ',file(1:LEN_TRIM(file))
c      CALL ER('GetSchematics','Problem with file',*99)
c 99   STOP
c      END
c
c ======================================================================
c
c subroutine: SelectTetrahedrons
c
c
      SUBROUTINE SelectTetrahedrons(nsur,npts2,vsur,
     .                              MAXSURFACE,MAXPOINTS,status)    
      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER, INTENT(IN) :: nsur,MAXSURFACE,MAXPOINTS,status,
     .                       npts2(0:MAXSURFACE)
      REAL*8 , INTENT(IN) :: vsur(3,MAXPOINTS,0:MAXSURFACE)

      REAL*8 CalcPerp

      INTEGER iobj,npts,isur,old_nobj
      REAL*8  r(2),phi(2),v(3,1000),p(3),t,pdist


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: TestTetrahedrons
c
c
      SUBROUTINE TestTetrahedrons(option,fname,nsur,npts,vsur,hsur,
     .                            MAXSURFACE,MAXPOINTS,status)
c      USE mod_eirene06_locals
      IMPLICIT none

      INTEGER option,nsur,MAXSURFACE,MAXPOINTS,status,
     .        npts(0:MAXSURFACE),hsur(0:MAXSURFACE),max_srf
      LOGICAL check
      REAL*8  vsur(3,MAXPOINTS,0:MAXSURFACE)
      CHARACTER fname*(*)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: SolidPlot985
c
c
      SUBROUTINE SolidPlot985(opt,nobj,obj,iplot)
      USE MOD_OUT985
      USE mod_out985_ribbon
      USE mod_filament
      USE mod_options
      USE mod_interface
      IMPLICIT none

c...  Input:
      TYPE(type_options985) :: opt
      INTEGER, INTENT(IN) :: nobj, iplot
      TYPE(type_3D_object) :: obj(nobj)

      RETURN
c 98   WRITE(0,*) 'ERROR DrawSolidPlot: File not found'
c      WRITE(0,*) '   FNAME= "'//fname(1:LEN_TRIM(fname))//'"'
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE WireframePlot985(iopt,MAXPIXEL,npixel,pixel,image,
     .                            iplot)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      INTEGER, INTENT(IN) :: iplot
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)

c      INCLUDE 'params'      
c      INCLUDE 'comgra'
c      INCLUDE 'colours'
c      INCLUDE 'slout'

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      CHARACTER glabel*512,caption*1024
      CHARACTER TITLE*80,TITLE2*80,JOB*72,GRAPH*80,GRAPH1*256
      CHARACTER*36 REF,NVIEW,PLANE,ANLY,TABLE,ZLABS(-2:93+1)
      CHARACTER*36 XLAB
      CHARACTER*72 YLAB
      CHARACTER*72 SMOOTH
      CHARACTER*256 cdum1

      INTEGER axis1(3)
      REAL    angle1(3),PI

      INTEGER CH1

      INTEGER nxbin,nybin,ipixel,nbin,isrf,isrf1,isrf2,isid,
     .        idet,nplot,iobj,i1,i2,i3,i4,ix,iy,iplot1,cnt
      REAL    deltar,deltaz,qmin,qmax,qval,ang1,ang2,ang,
     .        frac1,frac2,xcen,ycen,xnear,ynear,count
      REAL    XXMIN,XXMAX,YYMIN,YYMAX,dangle
      REAL*8  angle,p1(3,6),p2(3,6),x1,x2,z1,mat(3,3)

      REAL    CXMIN,CXMAX,CYMIN,CYMAX
      REAL    MAP1X,MAP2X,MAP1Y,MAP2Y
      
      INTEGER dat1,dat2
      REAL, ALLOCATABLE :: xdat(:),ydat(:)
      CHARACTER xlabel*256,ylabel*256,tag_x*2,tag_y*2,file*512,tag*7


c *TEMP*
      INTEGER, ALLOCATABLE :: nv(:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)

      PI = 3.141592

c...  PLOT A: 3D surface geometry plot for a given X,Y,Z along a particular view
c     -Translate and rotate as necessary:
c     -Project onto the R,Z plane:

      glabel = job  (CH1(job  ):LEN_TRIM(job  ))//'     '//
     .         graph(CH1(graph):LEN_TRIM(graph))

      ! jdemod - changed the : to ) to balance the parentheses in the format specifiers
      !          NOTE: these write statements will produce lines of 512 characters each
      !                right-padded with blanks 
      WRITE(nview ,'(512(A))') (' ',i1=1,LEN(nview )) 
      WRITE(plane ,'(512(A))') (' ',i1=1,LEN(plane )) 
      WRITE(anly  ,'(512(A))') (' ',i1=1,LEN(anly  )) 
      WRITE(table ,'(512(A))') (' ',i1=1,LEN(table )) 
      WRITE(xlab  ,'(512(A))') (' ',i1=1,LEN(xlab  )) 
      WRITE(ylab  ,'(512(A))') (' ',i1=1,LEN(ylab  )) 
      WRITE(smooth,'(512(A))') (' ',i1=1,LEN(smooth)) 
      WRITE(graph ,'(512(A))') (' ',i1=1,LEN(graph )) 
      WRITE(job   ,'(512(A))') (' ',i1=1,LEN(job   )) 

      XLAB = 'none'
      YLAB = 'none'
c      XLAB = '   R  (M)'
c      YLAB = '   Z  (M)'

      CALL THICK (1)
      CALL THICK2(1)

c      slopt = 1

c...  Set the "zoom" for the plot:
      xxmin =  0.0  ! Not sure if these are any good...
      xxmax =  6.0
      yymin = -3.0
      yymax =  3.0
      DO iplot1 = iplot+1, opt%nplots
        READ(opt%plots(iplot1),*) cdum1
        IF (cdum1(1:4).EQ.'zoom') THEN
          READ(opt%plots(iplot1),*) cdum1,xcen,ycen,xnear,ynear
          xxmin = xcen - xnear
          xxmax = xcen + xnear
          yymin = ycen - ynear
          yymax = ycen + ynear
        ELSE
          EXIT
        ENDIF   
      ENDDO      

      axis1 (1) = 1
      axis1 (2) = 2
      axis1 (3) = 3
      angle1(1) = 0.0
      angle1(2) = 0.0
      angle1(3) = 0.0

      nplot = 0
      iplot1 = iplot1 - 1
      cnt = 0
      DO WHILE (.TRUE.) 
c...    Wireframe:
        cnt = cnt + 1
        iplot1 = iplot1 + 1
        READ(opt%plots(iplot1),*) cdum1
c        WRITE(0,*) 'AXIS:',cdum1(1:LEN_TRIM(cdum1))
        IF (cdum1(1:4).EQ.'axis') THEN
          READ(opt%plots(iplot1),*) cdum1,(axis1(i1),angle1(i1),i1=1,3)
          nplot = nplot + 1
          IF (nplot.EQ.5) THEN
           nplot = 1
           CALL FRAME
          ENDIF
        ELSE
          EXIT
        ENDIF   

c...    Setup transformation matrix:
        CALL Calc_Transform2(mat,0.0D0,1,0)
        DO i1 = 1, 3
          angle = DBLE(angle1(i1)*3.141592/180.0)
          CALL Calc_Transform2(mat,angle,axis1(i1),1)
        ENDDO

c...    Use GRTSET_TRIM:
c        slopt4 = 1
c...    Stopping resizing of scale font in ghost1.o6a:
        iopt_ghost = 1

        SELECTCASE (nplot)
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

c        CALL GRTSET_TRIM(TITLE,REF,nVIEW,PLANE,glabel,
c     >                   xXMIN,xXMAX,yYMIN,yYMAX,
c     .                   TABLE,XLAB,YLAB,
c     .                   0,smooth,0,ANLY,1)
 
c...    Draw polygons:
        CXMIN = XXMIN
        CXMAX = XXMAX
        CYMIN = YYMIN
        CYMAX = YYMAX
        
        
        write(0,*) ' plot' ,map1x,map2x,map1y,map2y
          write(0,*) '     ' ,cxmin,cxmax,cymin,cymax
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)

        count = 0.0
        DO iobj = 1, nobj
c          IF (iobj.NE.1035) CYCLE  ! *TEMP*
c          IF (iobj.NE.471) CYCLE  ! *TEMP*
c          IF (iobj.NE.471.AND.iobj.NE.472.AND.iobj.NE.458) CYCLE  ! *TEMP*

          CALL LINCOL(1+obj(iobj)%colour)
c          CALL LINCOL(ncols+obj(iobj)%colour)            

          DO isid = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)

c            IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE

c            IF (isid.NE.2) CYCLE

c            IF (obj(iobj)%flag(isid).NE.-1.AND.
c     .          obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE  ! *TEMP*

c            IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE
c            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE
c            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL.AND.
c     .          obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE

c              WRITE(6,*) 'BOUNDARY?',obj(iobj)%ik,obj(iobj)%ir

c            IF (obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL) CYCLE   ! Recent filter...

c            IF (obj(iobj)%imap(1,isid).NE.0) CYCLE
c            IF (obj(iobj)%ik.NE.16) CYCLE


            IF (obj(iobj)%type.EQ.OP_INTEGRATION_VOLUME) THEN
              count = count + 
     .                1.0 / REAL(MAX(obj(iobj)%nsur,obj(iobj)%nside))

c              IF (count.GT.100.0) CYCLE
c              IF (count.GT.15000.0) CYCLE
            ENDIF

            IF (obj(iobj)%nside.NE.0) THEN
              isrf1 = obj(iobj)%iside(isid,1)
              isrf2 = obj(iobj)%iside(isid,2)
            ELSE
              isrf1 = 1
              isrf2 = 1
            ENDIF

            IF (obj(iobj)%gsur(isid).EQ.GT_TC) THEN
              dangle =  15.0 * PI / 180.0
              IF (.NOT..FALSE..AND.(.TRUE..OR.nplot.EQ.1)) THEN
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
c...            Filter:
                count = 0.0
c                IF (obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY) CYCLE
c                IF (obj(iobj)%imap(1,isid).NE.iobj) CYCLE
c                WRITE(0,*) 'GO MAN',iobj,isid
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
c...    Surface highlighting:
        DO iobj = 1, 0 ! nobj
          CALL LINCOL(1+obj(iobj)%colour+2)
c          CALL LINCOL(ncols+obj(iobj)%colour+2) 
          DO isid = 1, MAX(obj(iobj)%nsur,obj(iobj)%nside)

            IF ((obj(iobj)%tsur(isid).NE.SP_VESSEL_WALL.AND.
     .           obj(iobj)%tsur(isid).NE.SP_GRID_BOUNDARY).OR.
     .          obj(iobj)%esurf(isid).NE.21) CYCLE

            IF (obj(iobj)%nside.NE.0) THEN
              isrf1 = obj(iobj)%iside(isid,1)
              isrf2 = obj(iobj)%iside(isid,2)
            ELSE
              isrf1 = 1
              isrf2 = 1
            ENDIF

            IF     (obj(iobj)%gsur(isid).EQ.GT_TC) THEN
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
                STOP 'OUT OF DATE'
              ENDIF
            ELSE
              CALL ER('Plot985','Unknown surface geometry type',*98)
            ENDIF
          ENDDO
        ENDDO 
c...    Draw pixel views:
        IF (.NOT..TRUE.) THEN
          CALL LINCOL(1+55) 
          DO i1 = 1, MIN(1000,npixel)
            p1(1,1) = pixel(i1)%global_v1(1)
            p1(2,1) = pixel(i1)%global_v1(2)
            p1(3,1) = pixel(i1)%global_v1(3)
            p2(1,1) = pixel(i1)%global_v2(1)
            p2(2,1) = pixel(i1)%global_v2(2)
            p2(3,1) = pixel(i1)%global_v2(3)
            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))
            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ELSEIF (.TRUE.) THEN
c        ELSEIF (.TRUE..AND.nplot.GT.1) THEN
          CALL LINCOL(1+55) 
c          DO i1 = 1, MIN(1,nchord)
          DO i1 = 1, MIN(100000,nchord)
c          DO i1 = 1, MIN(1000,npixel)
            p1(1,1) = s_chord(i1)%v1(1)
            p1(2,1) = s_chord(i1)%v1(2)
            p1(3,1) = s_chord(i1)%v1(3)
            p2(1,1) = s_chord(i1)%v2(1)
            p2(2,1) = s_chord(i1)%v2(2)
            p2(3,1) = s_chord(i1)%v2(3)
c            IF (cnt.EQ.1) THEN
c              WRITE(0,*) '   v1=',s_chord(i1)%v1(1:3),i1
c              WRITE(0,*) '   v2=',s_chord(i1)%v2(1:3)
c            ENDIF
            call transform_vect(mat,p1(1,1))
            call transform_vect(mat,p2(1,1))
            CALL POSITN (SNGL(p1(1,1)),SNGL(p1(2,1)))
            CALL JOIN   (SNGL(p2(1,1)),SNGL(p2(2,1))) 
          ENDDO
        ENDIF
c...    Frame:
        CALL LINCOL(1)
c        CALL DrawFrame
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


      ENDDO


      RETURN
98    WRITE(0,*) 'OBJECT,SIDE,TYPE=',i1,i2,obj(i1)%gsur(i2)
99    STOP
      END
c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE ImagePlot985(iopt,MAXPIXEL,npixel,pixel,image,
     .                        iplot)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      INTEGER, INTENT(IN) :: iplot
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)


      RETURN
99    STOP
      END

c
c ======================================================================
c
c subroutine: DumpSurfaceCount
c
      SUBROUTINE DumpSurfaceCount(iplot)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER, INTENT(INOUT) :: iplot
      
      INTEGER       fp,iobj,isrf,isrf1,isrf2
      REAL, ALLOCATABLE :: obj_count(:)
      CHARACTER*256 cdum1,fname
      
      iplot = iplot + 1
      READ(opt%plots(iplot),*) cdum1
      IF (cdum1(1:4).EQ.'name') THEN
        READ(opt%plots(iplot),*) cdum1,fname
      ELSE
         WRITE(0,*) 'ERROR DumpSurfaceCount: Expecting file name' 
        iplot = iplot - 1
        RETURN
      ENDIF
     
      write(0,*) 'fname '//TRIM(fname)

      fp = 99
      OPEN(fp,FILE=fname(1:LEN_TRIM(fname)),
     .     FORM='FORMATTED',STATUS='REPLACE',ERR=98)     

      WRITE(fp,'(2A7,A12)')
     .  '*  iobj','isrf','count'           

      ALLOCATE(obj_count(nobj))
      obj_count = 0.0
      DO iobj = 1, nobj
        IF (obj(iobj)%reflec(1).NE.0) CYCLE
        isrf1 = obj(iobj)%iside(1,1)
        isrf2 = obj(iobj)%iside(1,2)
        DO isrf = isrf1, isrf2
          WRITE(fp,'(2I7,F12.2)') iobj,isrf,srf(isrf)%count
          obj_count(iobj) = obj_count(iobj) + srf(isrf)%count
       ENDDO
      ENDDO


      
      CLOSE(fp)
      
      RETURN
 98   WRITE(0,*) 'ERROR DumpSurfaceCount: Unable to open file' 
 99   STOP
      END
c
c ======================================================================
c
c subroutine: Output985
c
      SUBROUTINE Output985(iopt,MAXPIXEL,npixel,pixel,image)
      USE mod_out985
      USE mod_out985_variables
      USE mod_interface
      IMPLICIT none

      INTEGER iopt,MAXPIXEL,npixel         ! Put this into _variables... 
      TYPE(type_view) :: pixel(MAXPIXEL)
      REAL*8 image(1100,1100)

      INTEGER   iplot,option
      LOGICAL   debug
      CHARACTER buffer*1024

      debug = .FALSE.

      IF (debug) THEN
        WRITE(0,*) 'PLOT LIST:',opt%nplots
        DO iplot = 1, opt%nplots
          WRITE(0,*) TRIM(opt%plots(iplot))
        ENDDO
      ENDIF

      DO iplot = 1, opt%nplots
        WRITE(buffer,'(1024X)')
c        WRITE(0,*) 'PLOTS:',opt%plots(iplot)
        READ(opt%plots(iplot),*) buffer
        IF (buffer(1:4).NE.'plot') CYCLE
        READ(opt%plots(iplot),*) opt%plots(iplot),option
        SELECTCASE (option)
          CASE (000)
          CASE (001)
c....       Wireframe (all plots currently...):
            CALL WireframePlot985
     .             (iopt,MAXPIXEL,npixel,pixel,image,iplot)
            CALL FRAME
          CASE (002)
c....       Solid:
c            WRITE(0,*) 'DRAWING SOLID 3D PLOT'
            CALL SolidPlot985(opt,nobj,obj,iplot)
            CALL FRAME
          CASE (003)
c....       Image plot, and line shapes for the moment:
            CALL ImagePlot985
     .             (iopt,MAXPIXEL,npixel,pixel,image,iplot)
            CALL FRAME
          CASE (004)
c....       Image plot, and line shapes for the moment:
            CALL DumpSurfaceCount(iplot)
          CASE DEFAULT
            CALL WN('Output985','Unrecognised plot option')
            WRITE(0,*) '  BUFFER= "',buffer(1:LEN_TRIM(buffer)),'"'
        ENDSELECT
      ENDDO      

c      CALL FRAME

      RETURN
 99   STOP
      END
