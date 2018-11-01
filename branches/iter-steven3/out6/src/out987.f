c     -*-Fortran-*-
c
c
c
c ======================================================================
c
c subroutine: ExtractQuantity
c
      SUBROUTINE ExtractQuantity(ngauge,gauge,gauge_tag,tdata,iopt)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER, INTENT(IN) :: ngauge,iopt
      REAL   , INTENT(IN) :: gauge(5,ngauge),tdata(*)
      CHARACTER*256, INTENT(IN) :: gauge_tag(ngauge)

      INTEGER iobj,igauge,count,iside,isrf,ivtx
      REAL    obj_centroid(3,nobj),p(3),val,dist

      WRITE(0,*) 'nobj=',nobj
c...  Calculate the centres of each tetrahedron:
      DO iobj = 1, nobj
        IF (obj(iobj)%nside.NE.4) CYCLE
        p = 0.0
        count = 0
        DO iside = 1, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside,1))
          DO ivtx = 1, srf(isrf)%nvtx
            count = count + 1
            p(1:3) = p(1:3) + SNGL(vtx(1:3,srf(isrf)%ivtx(ivtx)))
          ENDDO
        ENDDO
        p(1:3) = p(1:3) / REAL(count)
        obj_centroid(1:3,iobj) = p(1:3)
      ENDDO

c...  Loop over gauges:
      DO igauge = 1, ngauge

c...    Convert gauge r,z,phi position to x,y,z:
        IF (ABS(gauge(4,igauge)+90.0).LT.1.0E-6.OR.
     .      ABS(gauge(4,igauge)-90.0).LT.1.0E-6) THEN
          p(1) = 0.0
          p(3) = gauge(2,igauge) * SIGN(1.0,gauge(4,igauge))
        ELSE
          p(1) = gauge(2,igauge) * COS(gauge(4,igauge)* 3.14159 / 180.0)
          p(3) = gauge(2,igauge) * SIN(gauge(4,igauge)* 3.14159 / 180.0)
        ENDIF
        p(2) = gauge(3,igauge)

        WRITE(0,*) 'iguage=',gauge(5,igauge),igauge
        WRITE(0,*) '       ',p(1:3)

c...    Check which tetrahedrons are within the gauge volume:
        count = 0
        val   = 0.0
        DO iobj = 1, nobj
          IF (obj(iobj)%nside.NE.4) CYCLE
          dist = SQRT((p(1)-SNGL(obj_centroid(1,iobj)))**2 + 
     .                (p(2)-SNGL(obj_centroid(2,iobj)))**2 + 
     .                (p(3)-SNGL(obj_centroid(3,iobj)))**2)
          IF (dist.LT.gauge(5,igauge)) THEN
            count = count + 1
c            WRITE(0,*) '    ---> go!',count,iobj,nobj
            val  = val + tdata(iobj)
          ENDIF   
        ENDDO
        IF (count.GT.0) val = val / REAL(count)

        WRITE(0,10)
     .    'GAUGE:',iopt,gauge(2:4,igauge),p(1:3),
     .             gauge(5,igauge),count,val,TRIM(gauge_tag(igauge))
        WRITE(6,10)
     .    'GAUGE:',iopt,gauge(2:4,igauge),p(1:3),
     .             gauge(5,igauge),count,val,TRIM(gauge_tag(igauge))
10      FORMAT(A,I4,2(3F8.2),2X,F8.2,I6,1P,E10.2,0P,2X,A)

      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SortPolygonPoints
c
      SUBROUTINE SortPolygonPoints(nv,rv,zv,problem_detected)
      IMPLICIT none

      INTEGER, INTENT(INOUT) :: nv
      LOGICAL, INTENT(INOUT) :: problem_detected
      REAL   , INTENT(INOUT) :: rv(nv),zv(nv)

      REAL*8, PARAMETER :: DTOL = 0.00001D0

      INTEGER i1,i2,i3,i4
      REAL*8  drv(nv),dzv(nv),t12,t34

      drv = DBLE(rv)
      dzv = DBLE(zv)

      i4 = -1

      DO i1 = 1, 3              

        i2 = i1 + 1
        IF (i2.EQ.nv) i2 = 1
        i3 = i2 + 1
        IF (i3.EQ.nv) i3 = 1

        CALL CalcInter(drv(i1),dzv(i1),drv(i2),dzv(i2),
     .                 drv(i3),dzv(i3),drv(nv),dzv(nv),t12,t34)

c        WRITE(0,*) 'og->',i1,i2,i3,nv,i4
c        WRITE(0,*) 'og->',drv(i1),dzv(i1)
c        WRITE(0,*) 'og->',drv(i2),dzv(i2)
c        WRITE(0,*) 'og->',drv(i3),dzv(i3)
c        WRITE(0,*) 'og->',drv(nv),dzv(nv)
c        WRITE(0,*) 'og->',t12,t34

        IF (t12.GE.0.0D0-DTOL.AND.t12.LE.1.0D0+DTOL.AND.
     .      t34.GE.0.0D0-DTOL.AND.t34.LE.1.0D0+DTOL) THEN
          IF (i4.EQ.-1) THEN
            i4 = i1
          ELSE
            i4 = -1
            EXIT
c            CALL ER('SortPolygonPoints','Something wrong',*99)
          ENDIF
        ENDIF 

      ENDDO
      IF (i4.EQ.-1) THEN
c...    No intersection found, which is strange, so reduce the 
c.      polygon to a triangle:
        nv = 3
        problem_detected = .TRUE.
      ELSE
c     .   CALL ER('SortPolygonPoints','Intersection not found',*99)

c      WRITE(0,*) 'i4=',i4

        DO i1 = i4+2, nv
          rv(i1) = SNGL(drv(i1-1))
          zv(i1) = SNGL(dzv(i1-1))
        ENDDO
        rv(i4+1) = SNGL(drv(nv))
        zv(i4+1) = SNGL(dzv(nv))      
      ENDIF

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: SpliceTetrahedrons
c
      SUBROUTINE SpliceTetrahedrons(ntet,tet_axis,tet_value,
     .                              itet,nv,rv,zv)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER, INTENT(INOUT) :: ntet
      INTEGER, INTENT(IN   ) :: tet_axis
      REAL   , INTENT(IN   ) :: tet_value
      INTEGER, INTENT(OUT  ) :: itet(*),nv(ntet)
      REAL   , INTENT(OUT  ) :: rv(4,ntet),zv(4,ntet)

      INTEGER AddVertex,AddSurface
      REAL    ATAN2C

      INTEGER iobj,iside,isrf,ivtx,ivtx2,i1,i2
      LOGICAL check_lt,check_gt,debug,status
      REAL    p1(3),val

c     For surface intersection code:
      INTEGER n,npts,point_kill
      REAL*8  v(3,MAXINTER),d(MAXINTER)
      REAL*8  DTOL, v1(3),v2(3),pts(3,10)

      TYPE(type_3D_object) :: newobj
      TYPE(type_surface  ) :: newsrf
      REAL*8                  newvtx(3,4)

      real za02as,t1
      external za02as

      point_kill = 0

      debug = .FALSE.

      status = .FALSE.

      DTOL = 1.0D-10

c..   Count the number of tetrahedrons that bound the surface of interest:
      IF (ntet.EQ.1) THEN 
        t1 = ZA02AS (1)
        write(0,*) 'counting bounding tetrahedrons'
        ntet = 0
        DO iobj = 1, nobj
          check_lt = .FALSE.
          check_gt = .FALSE.
          DO iside = 1, obj(iobj)%nside
            isrf = ABS(obj(iobj)%iside(iside,1))
            DO ivtx = 1, srf(isrf)%nvtx
              p1 = SNGL(vtx(1:3,srf(isrf)%ivtx(ivtx)))
              SELECTCASE (tet_axis)
                CASE(1)  ! phi
                  IF (ABS(p1(1)).LT.1.0E-6) THEN
                    val = 90.0 * SIGN(1.0,p1(3))
                  ELSE
                    val = ATAN2C(p1(3),p1(1)) * 180.0 / 3.141592
                  ENDIF
                CASE(2)  ! Y (height)
                  val = p1(2)
                CASE DEFAULT
                  CALL ER('SpliceTetrahedrons','Unknown option',*99)
              ENDSELECT

              IF (val.LT.tet_value) check_lt = .TRUE.
              IF (val.GT.tet_value) check_gt = .TRUE.
              IF (debug)
     .          WRITE(0,'(A,4I8,2X,3F14.7,2X,3F8.3,2L4)') 
     .            'val:',iobj,iside,isrf,ivtx,val,tet_value,
     .            obj(iobj)%phi,p1,
     .            check_lt,check_gt
            ENDDO
          ENDDO
          IF (check_lt.AND.check_gt) THEN
c...        The tetrahedron crosses the surface so add it to the list:
            ntet = ntet + 1
            itet(ntet) = iobj
          ENDIF 
c          IF (ntet.EQ.200) RETURN
        ENDDO
        write(0,*) 'done',ZA02AS(1)-t1,ntet
        RETURN
      ENDIF

c     STOP 'whit'

c...  Add a surface to the list of objects that corresponds to the plane,
c     so that the stardard surface intersection code in RAY can be called:
      SELECTCASE (tet_axis)
c       ----------------------------------------------------------------
        CASE(1)
          IF (ABS(tet_value+90.0).LT.1.0E-6.OR.
     .        ABS(tet_value-90.0).LT.1.0E-6) THEN
            p1(1) =  0.0
            p1(3) = 20.0 * SIGN(1.0,tet_value)
          ELSE
            p1(1) = 20.0 * COS(tet_value * 3.14159 / 180.0)
            p1(3) = 20.0 * SIN(tet_value * 3.14159 / 180.0)
          ENDIF
          newvtx(1:3,1) = (/       0.0D0 ,  20.0D0,       0.0D0 /)
c          newvtx(1:3,2) = (/-DBLE(p1(1)) ,  20.0D0,-DBLE(p1(3)) /)
          newvtx(1:3,2) = (/ DBLE(p1(1)) ,  20.0D0, DBLE(p1(3)) /)
          newvtx(1:3,3) = (/ DBLE(p1(1)) , -20.0D0, DBLE(p1(3)) /)
c          newvtx(1:3,3) = (/-DBLE(p1(1)) , -20.0D0,-DBLE(p1(3)) /)
          newvtx(1:3,4) = (/       0.0D0 , -20.0D0,       0.0D0 /)
c       ----------------------------------------------------------------
        CASE(2)
          newvtx(1:3,1) = (/  50.0D0 , DBLE(tet_value),  50.0D0 /)
          newvtx(1:3,2) = (/  50.0D0 , DBLE(tet_value), -50.0D0 /)
          newvtx(1:3,3) = (/ -50.0D0 , DBLE(tet_value), -50.0D0 /)
          newvtx(1:3,4) = (/ -50.0D0 , DBLE(tet_value),  50.0D0 /)
c       ----------------------------------------------------------------
        CASE DEFAULT
          CALL ER('SpliceTetrahedrons','Unknown option',*99)
c       ----------------------------------------------------------------
      ENDSELECT
      IF (debug) THEN
        WRITE(0,*) newvtx(1:3,1)
        WRITE(0,*) newvtx(1:3,2)
        WRITE(0,*) newvtx(1:3,3)
        WRITE(0,*) newvtx(1:3,4)
      ENDIF
      newsrf%type = SP_PLANAR_POLYGON
      newsrf%nvtx = 4
      DO i1 = 1, 4
        newsrf%ivtx(i1) = AddVertex(newvtx(1,i1))
      ENDDO
      newobj%nside        = 1
      newobj%iside(1,1:2) = AddSurface(newsrf)
      newobj%gsur(1:4)    = GT_TD
      newobj%nsur         = 0
      newobj%ipts(2,1)    = 0
      nobj = nobj + 1
      obj(nobj) = newobj

      pts = 0.0D0

      t1 = ZA02AS (1)
      write(0,*) 'surface intersection calculation'

      DO i1 = 1, ntet

        IF (MOD(i1,ntet/10).EQ.0) write(0,*) 'counting',i1,ntet

        iobj = itet(i1)  ! Set index to the next tetrahedron in the list
        npts = 0
     
        DO iside = 1, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside,1))
          DO ivtx = 1, srf(isrf)%nvtx
            ivtx2 = ivtx + 1
            IF (ivtx2.GT.srf(isrf)%nvtx) ivtx2 = 1
            v1 = vtx(1:3,srf(isrf)%ivtx(ivtx ))
            v2 = vtx(1:3,srf(isrf)%ivtx(ivtx2))

            n = 0
            CALL LineThroughSurface(v1,v2,nobj,1,nsrf,n,v,d,0,DTOL)
            IF (debug)            
     .         WRITE(0,'(A,4I6,2(2X,3F8.3),I6)') 
     .           '   :',iobj,iside,isrf,ivtx,v1,v2,n
            IF (n.EQ.1) THEN
              DO i2 = 1, npts
                IF (DABS(v(1,1)-pts(1,i2)).LT.DTOL.AND.
     .              DABS(v(2,1)-pts(2,i2)).LT.DTOL.AND.
     .              DABS(v(3,1)-pts(3,i2)).LT.DTOL) EXIT
              ENDDO
              IF (i2.EQ.npts+1) THEN
                npts = npts + 1
                pts(1:3,npts) = v(1:3,1)
                IF (debug) THEN
                  DO i2 = 1, npts
                    WRITE(0,'(A,I4,3F12.7)') '  -- ',i2,pts(1:3,i2)
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
 
          ENDDO
        ENDDO

        IF     (npts.EQ.3) THEN
          SELECTCASE (tet_axis)
c           ------------------------------------------------------------
            CASE(1)
              nv(    i1) = 3
              rv(1:3,i1) = SNGL(DSQRT(pts(1,1:3)**2 + pts(3,1:3)**2))
              zv(1:3,i1) = SNGL(pts(2,1:3))
c           ------------------------------------------------------------
            CASE(2)
              nv(    i1) = 3
              rv(1:3,i1) =  SNGL(pts(1,1:3))
              zv(1:3,i1) = -SNGL(pts(3,1:3))
c           ------------------------------------------------------------
            CASE DEFAULT
              CALL ER('SpliceTetrahedrons','Unknown option',*99)           
c           ------------------------------------------------------------
          ENDSELECT
        ELSEIF (npts.EQ.4) THEN 
          SELECTCASE (tet_axis)
c           ------------------------------------------------------------
            CASE(1)
              nv(    i1) = 4
              rv(1:4,i1) = SNGL(DSQRT(pts(1,1:4)**2 + pts(3,1:4)**2))
              zv(1:4,i1) = SNGL(pts(2,1:4))
c           ------------------------------------------------------------
            CASE(2)
              nv(    i1) = 4
              rv(1:4,i1) =  SNGL(pts(1,1:4))
              zv(1:4,i1) = -SNGL(pts(3,1:4))
c           ------------------------------------------------------------
            CASE DEFAULT
              CALL ER('SpliceTetrahedrons','Unknown option',*99)           
c           ------------------------------------------------------------
          ENDSELECT
c...      Now, need to order the points:

c...     
          CALL SortPolygonPoints(nv(i1),rv(1:4,i1),zv(1:4,i1),status)
        ELSE
c...      An error, kill point:
          nv(i1)     = 3
          rv(1:4,i1) = 0.0
          zv(1:4,i1) = 0.0
          point_kill = point_kill + 1
c        ELSE
c          CALL ER('SpliceTetrahedrons','Incorrect number of points',*99)
        ENDIF

      ENDDO

      write(0,*) 'done',ZA02AS(1)-t1,ntet

      WRITE(0,*) 'ntet=',ntet

      nvtx = nvtx - 4
      nsrf = nsrf - 1
      nobj = nobj - 1

      IF (status) THEN
        WRITE(0,*) '  --------------------------------------------'
        WRITE(0,*) '  WARNING SpliceTetrahedrons: Sorting problems'
        WRITE(0,*) '  --------------------------------------------'        
      ENDIF
      IF (point_kill.GT.0) THEN
        WRITE(0,*) '  --------------------------------------------'
        WRITE(0,*) '  WARNING SpliceTetrahedrons: kill=',point_kill
        WRITE(0,*) '  --------------------------------------------'        
      ENDIF

      RETURN
 99   WRITE(0,*) 'NPTS=',npts
      STOP
      END
c
c ======================================================================
c
c subroutine: GetTetrahedrons
c
      SUBROUTINE GetTetrahedrons(ntri)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER, INTENT(OUT) :: ntri

      INTEGER GetNumberOfObjects

      IF (nobj.GT.0) THEN
        WRITE(0,*) 'MESSAGE GetTetrahedrons: Not reloading since '//
     .             'objects exist'
        ntri = nobj
        RETURN
      ENDIF

      MAX3D = GetNumberOfObjects('eirene.transfer')
      ALLOCATE(obj(MAX3D+1))  ! Need to the space for an extra object for SpliceTetrahedrons
      nobj = 0
      CALL ALLOC_SURFACE(-1,MP_INITIALIZE)

      opt%obj_type  (1) = 6
      opt%obj_option(1) = 2
      opt%obj_fname (1) = 'tetrahedrons.raw'
      CALL ProcessTetrahedronGrid(1)

      WRITE(0,*) 'NOBJECTS=',nobj,MAX3D

c     STOP 'wtf?'

      ntri = nobj

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DrawGrid
c
c
      SUBROUTINE DrawGrid(iopt1)
      USE mod_eirene06_parameters ! 04
      USE mod_eirene06
      IMPLICIT none
      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      INTEGER iopt1

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      REAL       DTOL
      PARAMETER (DTOL=1.0D-06)

      INTEGER iopt,i1,i2,nline,lastcolour,v1,v2,ik,ir,id,idum1
      LOGICAL drawseparatrix
      INTEGER, ALLOCATABLE :: lines(:,:),lcolour(:)

      drawseparatrix = .FALSE.

      iopt = iopt1                      ! FIX REQUIRED FOR COMPILER HANGUP ON CHANGING SIGN OF IOPT, SEE BELOW...

      IF (iopt.LT.0) THEN
        iopt = -iopt
        drawseparatrix = .TRUE.
      ENDIF

      WRITE(0,*) 'IOPT:',iopt

c      CALL THICK2(10)

      IF (iopt.LT.500) THEN

        ALLOCATE(lines(3*ntri,2))
        ALLOCATE(lcolour(3*ntri))

        nline = 0
        DO i1 = 1, ntri
          DO v1 = 1, 3
            v2 = v1 + 1
            IF (v1.EQ.3) v2 = 1
c *TEMP*
            IF     (iopt.EQ.95) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.     ! Walls
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 2
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

c              drawseparatrix = .TRUE.

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.   ! Targets
     .            tri(i1)%sur(v1).NE.0.AND.
     .            tri(i1)%sideindex(2,v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.96) THEN

c             All triangles:
              nline = nline + 1
              lcolour(nline) = ncols + 21
              lines(nline,1)=tri(i1)%ver(v1)
              lines(nline,2)=tri(i1)%ver(v2)

            ELSEIF (iopt.EQ.97) THEN
              IF (i1.EQ.6509.OR.i1.EQ.6513) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.98) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.99) THEN
c...          Magnetic grid triangles + wall surfaces:
              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN


                nline = nline + 1
                lcolour(nline) = 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            v1.NE.1) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSE

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%index(2).EQ.irsep.AND.
     .            tri(i1)%sideindex(1,v1).EQ.14) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines(nline,1)=tri(i1)%ver(v1)
                lines(nline,2)=tri(i1)%ver(v2)
              ENDIF

           ENDIF

c            nline = nline + 1
c            lines(nline,1)=tri(i1)%ver(v1)
c            lines(nline,2)=tri(i1)%ver(v2)

          ENDDO
        ENDDO
c...    Remove duplicates:
        DO i1 = 1, nline-1
          DO i2 = i1+1, nline
            IF (lines(i1,1).NE.-999.0.AND.lines(i2,1).NE.-999.0) THEN
c * IMPROVE THIS CHECK! *
              IF ((DABS(ver(lines(i1,1),1)-
     .                  ver(lines(i2,2),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,1),2)-
     .                  ver(lines(i2,2),2)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),1)-
     .                  ver(lines(i2,1),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),2)-
     .                  ver(lines(i2,1),2)).LT.DTOL).OR.
     .            (DABS(ver(lines(i1,1),1)-
     .                  ver(lines(i2,1),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,1),2)-
     .                  ver(lines(i2,1),2)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),1)-
     .                  ver(lines(i2,2),1)).LT.DTOL.AND.
     .             DABS(ver(lines(i1,2),2)-
     .                  ver(lines(i2,2),2)).LT.DTOL)) THEN
                IF (lcolour(i1).EQ.1) THEN
                  lines(i1,1) = -999.0
                ELSE
                  lines(i2,1) = -999.0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i1 = nline, 1, -1
          IF (lines(i1,1).EQ.-999.0) THEN
c            WRITE(0,*) 'DELETING:',i1
            DO i2 = i1, nline-1
              lines(i2,1) = lines(i2+1,1)
              lines(i2,2) = lines(i2+1,2)
              lcolour(i2) = lcolour(i2+1)
            ENDDO
            nline = nline - 1
          ENDIF
        ENDDO
      ENDIF




      IF (.TRUE.) THEN
c...    Plot polygons:
        CALL PSPACE (map1x,map2x,map1y,map2y)      
        CALL MAP    (cxmin,cxmax,cymin,cymax)
        lastcolour = -1
        DO i1 = 1, nline
          IF (lastcolour.NE.lcolour(i1)) THEN
            CALL LINCOL(lcolour(i1)) 
            lastcolour = lcolour(i1)
          ENDIF
          CALL POSITN(SNGL(ver(lines(i1,1),1)),SNGL(ver(lines(i1,1),2)))
          CALL JOIN  (SNGL(ver(lines(i1,2),1)),SNGL(ver(lines(i1,2),2)))
        ENDDO
        DEALLOCATE(lines)
        DEALLOCATE(lcolour)
      ENDIF



      IF (drawseparatrix) THEN

        CALL THICK2(1)
        CALL THICK (1)
        CALL BROKEN(5,5,5,5)

c        CALL LINCOL(ncols+2) 
        CALL LINCOL(ncols+4) 
  
        ir = irsep
        DO ik = 1, nks(ir)
          id = korpg(ik,ir)
          CALL POSITN(rvertp(1,id),zvertp(1,id))
          CALL JOIN  (rvertp(4,id),zvertp(4,id))
        ENDDO

        IF (nrs.EQ.65) THEN
          ir = 38
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            CALL POSITN(rvertp(2,id),zvertp(2,id))
            CALL JOIN  (rvertp(3,id),zvertp(3,id))
          ENDDO
        ENDIF

        IF (nrs.EQ.72) THEN
          ir = 47
          DO ik = 1, nks(ir)
            id = korpg(ik,ir)
            CALL POSITN(rvertp(2,id),zvertp(2,id))
            CALL JOIN  (rvertp(3,id),zvertp(3,id))
          ENDDO
        ENDIF

c        IF (irwall.GT.23) THEN
c          ir = 31
c          DO ik = 1, nks(ir)
c            id = korpg(ik,ir)
c            CALL POSITN(rvertp(1,id),zvertp(1,id))
c            CALL JOIN  (rvertp(4,id),zvertp(4,id))
c          ENDDO
c          ir = 40
c          DO ik = 1, nks(ir)
c            id = korpg(ik,ir)
c            CALL POSITN(rvertp(1,id),zvertp(1,id))
c            CALL JOIN  (rvertp(4,id),zvertp(4,id))
c          ENDDO
c        ENDIF

        CALL FULL
      ENDIF




c...  Frame:
      CALL DrawFrame


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DrawColourScale
c
c
      SUBROUTINE DrawColourScale(mode,colmode,qmin,qmax,label)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'slout'

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

      IF     (mode.EQ.0) THEN  ! Don't draw a scale
      ELSEIF (mode.EQ.1) THEN  ! Vertical scale drawn to the left of the plot
        dspot = 0.016         
        minx = map2x + 0.02
        maxx = map2x + 0.04
        miny =  HI
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
          IF     (qval.GT.-1.0.AND.qval.LT.10.0) THEN
            WRITE(nums,'(F4.1)') qval
          ELSEIF (ABS(qval).LT.100.0) THEN
            WRITE(nums,'(I3)') NINT(qval)
          ELSEIF (ABS(qval).LT.1000.0) THEN
            WRITE(nums,'(I4)') NINT(qval)
          ELSE
            WRITE(nums,'(1P,E9.1)') qval
          ENDIF
c          IF (qmax.GT.1.0.AND.qmax.LT.100.0) THEN
c            WRITE(nums,'(F4.1)') 
c     .        qval
cc            WRITE(nums,'(F4.1,A,I3,A)') 
cc     .        qval,' (',NINT(qval/qmax*100.0),'%)'
c          ELSE
c            WRITE(nums,'(1P,E10.2,0P,A,I3,A)') 
c     .        qval,' (',NINT(qval/qmax*100.0),'%)'
c          ENDIF
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
          minx = HI                ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.030     ! map1y - 0.04      ! HI
          maxy = map1y - 0.010     ! map1y - 0.02      ! map2y - dspot/2.0
        ELSE
          minx = HI                ! map2x + 0.02
          maxx = map2x - dspot/2.0 ! map2x + 0.04
          miny = map1y - 0.075     ! map1y - 0.04      ! HI
          maxy = map1y - 0.055     ! map1y - 0.02      ! map2y - dspot/2.0
        ENDIF
        dscale = 2.0
        DO rscale = 100.0-dscale, 0.0, -dscale
          qval = (rscale + 0.5 * dscale) / 100.0 * (qmax - qmin) + qmin
          IF (qmin.NE.qmax) CALL SetCol255_04(colmode,qval,qmin,qmax)
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
          IF     (qval.GT.-0.999.AND.qval.LT.1.0) THEN
c          IF     (qmax.GT.0.1.AND.qmax.LE.1.0) THEN
            WRITE(nums,'(F5.2)') qval
          ELSEIF (qval.GT.-0.999.AND.qval.LT.10.0) THEN
c          IF     (qmax.GT.0.1.AND.qmax.LE.1.0) THEN
            WRITE(nums,'(F4.1)') qval
          ELSEIF (ABS(qval).LT.999.0) THEN
c          ELSEIF ((qmax.GT. 1.0.AND.qmax.LE. 999.9).OR.
c     .            (qmax.LE.-1.0.AND.qmax.GT.-999.9)) THEN
            WRITE(nums,'(I4)') NINT(qval)
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

        frac = (qval - qmin) / (qmax - qmin + 1.0E-10)
        frac5 = 100.0*frac
        fmod5 = AMOD(frac5,2.0)
        frac = MIN(0.98,(frac5-fmod5)/100.0)

        bright = 1.0-(0.98-frac)**20

        IF (mode.NE.last_mode.OR.frac.NE.last_frac.OR.
     .      bright.NE.last_bright) 
     .    CALL ColSet(1.0-0.75*frac,1.0,bright,255)
c     .    CALL ColSet(0.75*frac+0.25,1.0,bright,255)

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
c
      SUBROUTINE NextLine(fp,ntally,icount,rdum,binary)
      IMPLICIT none

      INTEGER   fp,ntally,icount,i1
      LOGICAL   binary
      REAL      rdum(*)   
      CHARACTER buffer*256

      REAL*8    ddum(ntally)

      DO WHILE (.TRUE.) 
        IF (binary) THEN
        ELSE
          READ(fp,'(A256)',END=98) buffer               
        ENDIF
        IF (.NOT.binary.AND.buffer(1:1).EQ.'*') CYCLE 
        IF (binary) THEN
          READ(fp      ,ERR=97) icount,(rdum(i1),i1=1,ntally)          
        ELSE
          READ(buffer,*,ERR=97) icount,(ddum(i1),i1=1,ntally)          
        ENDIF
        DO i1 = 1, ntally
          IF (binary) THEN
            IF (rdum(i1).GT.1.0D+30) THEN
              STOP 'NOT SURE WHAT TO DO HERE'
            ENDIF
          ELSE
            IF (ddum(i1).GT.1.0D+30) THEN
              WRITE(0,*) 'WARNING NextLine: EIRENE data beyond '//
     .                   'single precision size limit, setting to zero'
              rdum(i1) = 0.0
            ELSE
              rdum(i1) = SNGL(ddum(i1))
            ENDIF
          ENDIF
        ENDDO
        RETURN
      ENDDO
      
 97   WRITE(0,*) 'buffer >'//TRIM(buffer)//'<'
      CALL ER('NextLine','Data format error',*99)
 98   CALL ER('NextLine','Unexpected end-of-file',*99)
 99   STOP
      END

c      SUBROUTINE NextLine(fp1,ntally,icount,rdum)
c      IMPLICIT none
c
c      INTEGER   fp,fp1,ntally,icount,i1
c      REAL      rdum(*)   
c      CHARACTER buffer*512
c      LOGICAL output
c
cc      output = .FALSE.
cc      IF (fp1.EQ.44) THEN
cc        fp = 99
cc        output = .TRUE.
cc      ELSE
c        fp = fp1
cc      ENDIF
c
c      DO WHILE (.TRUE.) 
c        READ(fp,'(A512)',END=98) buffer               
cc        IF (output) WRITE(0,*) 'BUFFER:',icount,buffer(1:50)
c        IF (buffer(1:1).EQ.'*') CYCLE 
c        READ(buffer,*,ERR=97) icount,(rdum(i1),i1=1,ntally)          
cc        IF (output) WRITE(0,*) 'BUFFER:',icount,buffer(1:50)
c        RETURN
c      ENDDO
c      
c 97   CALL ER('NextLine','Data format error',*99)
c 98   CALL ER('NextLine','Unexpected end-of-file',*99)
c 99   STOP
c      END
c
c ======================================================================
c
      SUBROUTINE Plot987(job,graph,ref,title,iopt,
     .                   xxmin,xxmax,yymin,yymax,ft,fp,zadj,
     .                   ismoth,ignors,itec,avs,navs,nizs)
      USE mod_interface 
      USE mod_eirene06_parameters
      USE mod_eirene06 
      USE mod_out985_clean 
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'cedge2d'
      INCLUDE 'pindata'
      INCLUDE 'dynam2'
      INCLUDE 'dynam3'
      INCLUDE 'comgra'
      INCLUDE 'colours'
      INCLUDE 'slcom'
      INCLUDE 'slout'

      REAL       TOL
      PARAMETER (TOL=1.0E-06)

      REAL*8     DTOL
      PARAMETER (DTOL=1.0D-06)

      INTEGER CH1
      REAL    GetMach

      INTEGER   ismoth,IGNORS(MAXNGS),ITEC,NAVS,iopt,nizs
      REAL      XXMIN,XXMAX,YYMIN,YYMAX,ft,fp,zadj,AVS(0:100)
      CHARACTER TITLE*(*),JOB*(*),GRAPH*(*),REF*(*)

      COMMON /GHOSTCOM/ iopt_ghost
      INTEGER           iopt_ghost

      COMMON /PLOT982COM/ hardscale,colscale
      LOGICAL             hardscale,colscale

      CHARACTER table*36,nview*36,anly*36,plane*36,fname*256,
     .          YLABEL*256,XLABEL*36,smooth*64,glabel*512,
     .          file_id*256

      INTEGER i1,i2,v1,v2,nc,ik,ir,id,iz,tet_axis,ntet,
     .        nline,lastcolour,scaleopt,colouropt,posopt,ngauge
      LOGICAL setqmin,setqmax,inside,scale_set,double_null,
     .        tetrahedrons_dump,file_exist
      REAL    qmin,qmax,frac,frac5,fmod5,scalefact,fact,
     .        posx,posy,poswidth,posheight,rdum(4),taus,tet_value,
     .        gauge(5,100),val,p(3),r
      CHARACTER label*512,cdum1*512,cdum2*512,tet_id*256
      CHARACTER*256 gauge_tag(100)

      REAL, POINTER :: gdata(:,:)

      INTEGER, ALLOCATABLE :: lines2(:,:),nv(:),lcolour(:),itet(:)
      REAL   , ALLOCATABLE :: tdata(:),tdata1(:)
      REAL, TARGET, ALLOCATABLE :: gdata1(:,:)
      REAL   , ALLOCATABLE :: rv(:,:),zv(:,:),cq(:)

c...  ADAS:
      CHARACTER ADASID*80,adasex*3,graph3*80
      integer   adasyr,ISELE,ISELR,ISELX,iseld,ierr,ircode
      REAL      wlngth

      INTEGER numplots
      REAL    dx,dy,dist,save_map2x
      DATA    numplots /0/
 
      LOGICAL :: reset_origin = .TRUE.

      real za02as,t1
      external za02as

      SAVE

c      WRITE(0,*) 'DATA:',job
c      WRITE(0,*) 'DATA:',graph
c      WRITE(0,*) 'DATA:',title
c      WRITE(0,*) 'DATA:',ref

      double_null = .FALSE.   ! *** THIS DOESN'T WORK BECAUSE RINGTYPE IS NOT PASSED TO OUT! ***
      DO ir = 2, irwall-1
        IF (ringtype(ir).EQ.PFZ) double_null = .TRUE.
      ENDDO

      WRITE(glabel,'(512(A:))') (' ',i1=1,LEN(glabel)) 

      IF (numplots.EQ.0) 
     .  glabel = job  (CH1(job  ):LEN_TRIM(job  ))//'     '//
     .           graph(CH1(graph):LEN_TRIM(graph))

      WRITE(nview ,'(512(A:))') (' ',i1=1,LEN(nview )) 
      WRITE(plane ,'(512(A:))') (' ',i1=1,LEN(plane )) 
      WRITE(anly  ,'(512(A:))') (' ',i1=1,LEN(anly  )) 
      WRITE(table ,'(512(A:))') (' ',i1=1,LEN(table )) 
      WRITE(xlabel,'(512(A:))') (' ',i1=1,LEN(xlabel)) 
      WRITE(ylabel,'(512(A:))') (' ',i1=1,LEN(ylabel)) 
      WRITE(smooth,'(512(A:))') (' ',i1=1,LEN(smooth)) 
c      WRITE(graph ,'(512(A:))') (' ',i1=1,LEN(graph )) 
c      WRITE(job   ,'(512(A:))') (' ',i1=1,LEN(job   )) 


      xlabel = 'z (m)'
      ylabel = 'r (m)'
c      xlabel = 'none'
c      ylabel = 'none'

c...  Use GRTSET_TRIM:
      slopt4 = 1
c...  Stopping resizing of scale font in ghost1.o6a:
      iopt_ghost = 1

c      cxmin=xxmin     
c      cxmax=xxmax 
c      cymin=yymin 
c      cymax=yymax 

c      CALL SLSET (0.05,0.95,0.05,0.95,
c     .            xxmin,xxmax,yymin,yymax)
c      CALL PSPACE (map1x,map2x,map1y,map2y)

      scale_set = .FALSE.
      scaleopt   = 2
      colouropt  = 1
      scalefact  = 1.0
      label = 'default'

      qmin =  HI
      qmax = -HI
c...  Read scale information:
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:12).EQ.'Scale'.OR.cdum1(8:12).EQ.'scale'.OR.
     .      cdum1(8:12).EQ.'SCALE') THEN
        scale_set = .TRUE.
c        WRITE(0,*) 'CDUM1:>'//TRIM(cdum1)//'<'
        READ(cdum1,*) cdum2,scaleopt,colouropt,scalefact,qmin,qmax,
     .                label
        IF (qmin.EQ.-99.0) qmin =  HI
        IF (qmax.EQ.-99.0) qmax = -HI
c        WRITE(0,*) 'SCALE:',scaleopt,colouropt,scalefact,
c     .              label(1:LEN_TRIM(label))
        IF (label.EQ.'default') 
     .    label = graph(CH1(graph):LEN_TRIM(graph))
      ELSE
        BACKSPACE 5
      ENDIF

c...  Plot location:      
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:15).EQ.'Position'.OR.cdum1(8:15).EQ.'position'.OR.
     .      cdum1(8:15).EQ.'POSITION') THEN
        READ(cdum1,*) cdum2,posopt,posx,posy,poswidth,posheight

c        WRITE(0,*) 'POSITION',posopt,posx,posy,poswidth,posheight

        IF     (posopt.EQ.1) THEN
          map1x = 0.05
          map2x = map1x + 0.80 * (xxmax-xxmin) / (yymax-yymin)
c          map2x = map1x + 0.75 * (xxmax-xxmin) / (yymax-yymin)
          map1y = 0.15
          map2y = map1y + 0.80
c          map2y = map1y + 0.75
c          map1x = 0.05
c          map2x = map1x + 0.40 * (xxmax-xxmin) / (yymax-yymin)
c          map1y = 0.15
c          map2y = map1y + 0.40
        ELSEIF (posopt.EQ.2) THEN
          map1x = posx
          map2x = map1x + poswidth * (xxmax-xxmin) / (yymax-yymin)
          map1y = posy
          map2y = map1y + poswidth
        ELSEIF (posopt.EQ.3) THEN
          map1x = posx
          map2x = map1x + poswidth 
          map1y = posy
          map2y = map1y + posheight
        ELSE
          STOP 'BAD POISTION OPTION'
        ENDIF
      ELSE
        BACKSPACE 5
        dist = 0.70 ! 0.80
        dx = dist * (xxmax-xxmin) / (yymax-yymin)
        IF (reset_origin) THEN
          map1x = 0.05 
        ELSE
          map1x = save_map2x + 0.04
        ENDIF
c        map1x = 0.05 + REAL(numplots) * dx
        map2x = map1x + dx
        map1y = 0.20
        map2y = map1y + dist
        save_map2x = map2x
      ENDIF

      CALL FULL
c      CALL THICK2(6)
      CALL THICK2(1)
      CALL THICK(1)

c      write(0,*) 'map',map1x,map2x,map1y,map2y

      IF (numplots.NE.0) ylabel = 'none'




c...  Decide if tetrahedrons are being plotted:      
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:12).EQ.'Tetra'.OR.cdum1(8:12).EQ.'tetra'.OR.
     .      cdum1(8:12).EQ.'TETRA') THEN
        READ(cdum1,*) cdum2,tet_axis,tet_value  ! 1 = phi, 2 = y, 3 = r
        tetrahedrons = .TRUE.
      ELSE
        BACKSPACE 5
        tetrahedrons = .FALSE.
      ENDIF
c...  Decide if tetrahedrons are being plotted:      
      READ(5,'(A512)') cdum1
      IF   (cdum1(8:11).EQ.'Dump'.OR.cdum1(8:11).EQ.'dump'.OR.
     .      cdum1(8:11).EQ.'DUMP') THEN
        READ(cdum1,*) cdum2,tet_id
        tetrahedrons_dump = .TRUE.
      ELSE
        BACKSPACE 5
        tetrahedrons_dump = .FALSE.
      ENDIF

c...  Set the file name of the triangle/tetrahedron data file:
      fname = 'default'
      file_id = ' '

      READ(5,'(A512)') cdum1
      IF   (cdum1(8:11).EQ.'File'.OR.cdum1(8:11).EQ.'file'.OR.
     .      cdum1(8:11).EQ.'FILE') THEN
        READ(cdum1,*) cdum2,file_id

        fname = 'eirene.'//TRIM(file_id)//'.transfer'

        INQUIRE(FILE=fname,EXIST=file_exist)                                

        IF (.NOT.file_exist) THEN
        
          CALL CollectTransferFile(TRIM(file_id),TRIM(fname))

        ENDIF

      ELSE
        BACKSPACE 5
      ENDIF



      CALL GRTSET_TRIM (TITLE//' '//TRIM(file_id),' ',' ',' ',glabel,
     >                  xXMIN,xXMAX,yYMIN,yYMAX,
     .                  ' ',xlabel,ylabel,
     .                  0,' ',0,' ',1)

c      CALL GRTSET_TRIM (TITLE,REF,nVIEW,PLANE,glabel,
c     >                  xXMIN,xXMAX,yYMIN,yYMAX,
c     .                  TABLE,XLABEL,YLABEL,
c     .                  0,smooth,0,ANLY,1)



c     ==================================================================
      IF (iopt.GE.500) THEN
c     ==================================================================

        ALLOCATE(gdata1(MAXNKS,MAXNRS))
        gdata1 = 0.0

c        WRITE(0,*) 'DEBUG: IOPT',iopt

        SELECTCASE (iopt)
          CASE (501)
            gdata => e2dion
          CASE (502)
            gdata => e2drec
          CASE (520)  
            gdata => e2dnbs
          CASE (521) 
            gdata => e2dtebs
          CASE (522) 
            gdata => e2dtibs
          CASE (530)  ! Cross-field metric THETAG 
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              gdata1(1:nks(ir),ir) = thetag(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (531)  ! Plot of BRATIO
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              DO ik = 1, nks(ir)
                gdata1(ik,ir) = 1.0 / bratio(ik,ir) ! ksb(ik,ir) - ksb(ik-1,ir)
              ENDDO
            ENDDO
            gdata => gdata1

          CASE (540) 
            gdata1 = 0.0
            DO ir = 2, nrs ! irsep-1
              gdata1(1:nks(ir),ir) = pinion(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (542) 
            gdata1 = 0.0
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = pinalpha(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (599) 
            gdata => e2dcxrec  ! Testing
          CASE (600)  ! Atom density
            gdata => pinatom
          CASE (601)  ! Molecule density
            gdata => pinmol
          CASE (700)  ! Mach no. (absolute)
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              DO ik = 1, nks(ir)
                gdata1(ik,ir) = GetMach
     .                    (kvhs(ik,ir)/qtim,ktebs(ik,ir),ktibs(ik,ir)) * 
     .                    SIGN(1.0,kvhs(ik,ir))
              ENDDO
            ENDDO
            gdata => gdata1
          CASE (701)  ! Velocity (absolute)
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(:,ir) = ABS(kvhs(:,ir)) / qtim
            ENDDO
            gdata => gdata1
          CASE (708) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = ktebs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1	
          CASE (710) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = knbs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (720) 
            gdata1 = 0.0
c            DO ir = irsep, nrs
            DO ir = 2, nrs
              IF (idring(ir).EQ.BOUNDARY) CYCLE
              gdata1(1:nks(ir),ir) = ktibs(1:nks(ir),ir)
            ENDDO
            gdata => gdata1
          CASE (799:875) 
              iz = iopt - 800
              READ(5,'(A512)') cdum1
c             ----------------------------------------------------------
              IF (cdum1(8:11).EQ.'Adas'.OR.cdum1(8:11).EQ.'ADAS'.OR.
     .            cdum1(8:11).EQ.'adas') THEN
c...            Load PLRP data from ADAS:
                BACKSPACE 5
                CALL RDG1 (GRAPH3,ADASID,adasyr,adasex,
     .                     ISELE,ISELR,ISELX,ISELD,IERR)
                WRITE(0,*) 'ADAS SETTINGS:',adasid,adasyr,adasex
                WRITE(0,*) '             :',isele,iselr,iselx,iseld
                IF (IERR.NE.0) THEN
                  WRITE(6,*) '987: ERROR READING ADAS DETAILS, '//
     .                       'IERR = ',IERR
                  IERR = 0
                  GOTO 99
                ENDIF
                WRITE(0,*) '             :',cion,iz
                CALL LDADAS(cion,IZ,ADASID,ADASYR,ADASEX,ISELE,ISELR,
     .                      ISELX,gdata1,Wlngth,IRCODE)
                WRITE(0,*) 'ADAS DATA:',iz,wlngth,ircode


                WRITE(6,*) 'data dump de James'
                WRITE(6,*) 'wavelength',wlngth/10.0

                write(6,'(2A6,4A11,A)') 
     .            'cell','ring','s (m)','smax-s (m)','p (m)',
     .            'pmax-p (m)','  emission (arb.)'
                DO ir = irsep, nrs
                  write(6,*) nks(ir)
                  DO ik = 1, nks(ir)
                    write(6,'(2I6,4F11.3,1P,E10.2,0P)') 
     .                ik,ir,
     .                kss(ik,ir),ksmaxs(ir)-kss(ik,ir),
     .                kps(ik,ir),kpmaxs(ir)-kps(ik,ir),
     .                gdata1(ik,ir)
                  ENDDO
                ENDDO

c             ----------------------------------------------------------
              ELSEIF (cdum1(8:10).EQ.'Ion'.OR.cdum1(8:10).EQ.'ION'.OR.
     .                cdum1(8:10).EQ.'ion') THEN
c                WRITE(0,*) 'Loading ionisation data'
                DO ir = 2, nrs
                  IF (idring(ir).EQ.BOUNDARY) CYCLE
                  gdata1(1:nks(ir),ir) = tizs(1:nks(ir),ir,iz)
                ENDDO
c             ----------------------------------------------------------
              ELSEIF (cdum1(8:12).EQ.'Power'.OR.cdum1(8:12).EQ.'POWER'
     .                .OR.cdum1(8:12).EQ.'power') THEN
                WRITE(0,*) 'Loading total radiated power data',iz,
     .                     MIN(cion,nizs)
                IF (iz.EQ.-1) THEN
                  DO iz = 0, MIN(cion,nizs)
                    DO ir = 2, nrs
                      IF (idring(ir).EQ.BOUNDARY) CYCLE
                      gdata1(1:nks(ir),ir) = gdata1(1:nks(ir),ir   ) + 
     .                                       powls (1:nks(ir),ir,iz)
                    ENDDO
                  ENDDO
                ELSE
                  DO ir = 2, nrs
                    IF (idring(ir).EQ.BOUNDARY) CYCLE
                    gdata1(1:nks(ir),ir) = powls(1:nks(ir),ir,iz) 
                  ENDDO
                ENDIF
                gdata1 = gdata1 * absfac
c             ----------------------------------------------------------
              ELSEIF (cdum1(8:15).EQ.'Legrange'.OR.    ! Net force on impurities, from OUT 
     .                cdum1(8:15).EQ.'LEGRANGE'.OR.    ! plot 669/670
     .                cdum1(8:15).EQ.'legrange') THEN
                FACT = QTIM**2 * EMI / CRMI
                DO ir = 2, nrs ! irsep-1
                  DO ik = 1, nks(ir)
                    TAUS = CRMI * KTIBS(IK,IR)**1.5 * SQRT(1.0/CRMB) /
     +                     (6.8E-14 * (1 + CRMB / CRMI) * KNBS(IK,IR) *
     +                     REAL(IZ)**2.0 * RIZB**2 * 15.0)
                    RDUM(1) = AMU * CRMI * KVHS(IK,IR) / QTIM / TAUS
                    RDUM(2) = KFIGS(IK,IR) * KBETAS(IZ) * ECH / FACT
                    RDUM(3) = KFEGS(IK,IR) * KALPHS(IZ) * ECH / FACT
                    RDUM(4) = REAL(IZ) * KES(IK,IR) * ECH / FACT
                    WRITE(6,'(A,2I6,5E10.2)') 
     .                'FORCES:',ik,ir,rdum(1:4),
     .                ABS(sum(rdum(1:4)))*scalefact
                    gdata1(ik,ir) = ABS(SUM(rdum(1:4)))
                  ENDDO
                ENDDO
                gdata => gdata1
c             ----------------------------------------------------------
              ELSE
c...            Load impurity density data:
c                WRITE(0,*) 'DEBUG: here 2'
c                   IF (iz.EQ.0) 
c     .                WRITE(0,*) sdlims(1:nks(109),109,iz)
                BACKSPACE 5
                DO ir = 2, nrs
                  IF (idring(ir).EQ.BOUNDARY) CYCLE
                  gdata1(1:nks(ir),ir) = sdlims(1:nks(ir),ir,iz)
c                  write(6,*) 'iz madness',iz
c                  do ik = 1, nks(ir)
c                    write(6,*) 'sdlims',sdlims(ik,ir,iz)
c                  enddo
                ENDDO
              ENDIF
              gdata => gdata1
c             ----------------------------------------------------------
          CASE DEFAULT 
            CALL ER('Plot987','Unrecognized option',*99)
        ENDSELECT

        ALLOCATE(nv(MAXNKS*MAXNRS))
        ALLOCATE(rv(4,MAXNKS*MAXNKS))
        ALLOCATE(zv(4,MAXNKS*MAXNKS))
        ALLOCATE(cq(MAXNKS*MAXNRS))

c...    Load up data:
        nc = 0
        nv = 0
c       DO ir = irsep, nrs
        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
c          IF (ir.LT.irsep.OR.ir.GE.irtrap) CYCLE
          IF ((iopt.EQ.520.OR.iopt.EQ.521.OR.iopt.EQ.522.OR.
     .         iopt.EQ.599).AND.
     .        ir.LT.irsep) CYCLE

          DO ik = 1, nks(ir)
            IF (ir.LT.irsep.AND.ik.EQ.nks(ir)) CYCLE
            id = korpg(ik,ir)
            nc = nc + 1
            nv(nc) = nvertp(id)
            DO i1 = 1, nvertp(id)
              rv(i1,nc) = rvertp(i1,id)
              zv(i1,nc) = zvertp(i1,id)
            ENDDO
            cq(nc) = gdata(ik,ir)
          ENDDO
        ENDDO
c     ==================================================================
      ELSEIF (iopt.LT.500) THEN
c     ==================================================================

        IF (.NOT.tetrahedrons) THEN
          CALL LoadTriangles_06  
        ELSE
          CALL GetTetrahedrons(ntri)
          write(0,*) 'ntri=',ntri
        ENDIF

        t1 = ZA02AS (1)
        write(0,*) 'loading data' 

        ALLOCATE(tdata(ntri+1))
        tdata = 0.0
        IF (iopt.EQ.1 ) CALL LoadTriangleData(2,1,1,1,tdata,fname)  ! D  density
        IF (.FALSE..AND.iopt.EQ.1) THEN
          DO i1 = 1, ntri 
            IF (ver(tri(i1)%ver(1),1).GT.1.80D0) THEN
c     .          ver(tri(i1)%ver(1),1).LT.1.98D0) THEN
              IF (ver(tri(i1)%ver(1),2).LT.-1.80D0) 
     .          WRITE(0,*) 'LOWER DIVERTOR n_D :',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT.-0.3D0.AND.
     .            ver(tri(i1)%ver(1),2).LT. 0.3D0) 
     .          WRITE(0,*) 'MIDPLANE       n_D :',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT. 1.80D0) 
     .          WRITE(0,*) 'UPPER DIVERTOR n_D :',tdata(i1)
            ENDIF
          ENDDO
        ENDIF
        IF (iopt.EQ.2 ) CALL LoadTriangleData(3,1,1,1,tdata,fname)  ! D2 density
        IF (.FALSE..AND.iopt.EQ.2) THEN
          DO i1 = 1, ntri 
            IF (ver(tri(i1)%ver(1),1).GT. 1.80D0) THEN
              IF (ver(tri(i1)%ver(1),2).LT.-1.80D0) 
     .          WRITE(0,*) 'LOWER DIVERTOR n_D2:',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT.-0.3D0.AND.
     .            ver(tri(i1)%ver(1),2).LT. 0.3D0) 
     .          WRITE(0,*) 'MIDPLANE       n_D2:',tdata(i1)
              IF (ver(tri(i1)%ver(1),2).GT. 1.80D0) 
     .          WRITE(0,*) 'UPPER DIVERTOR n_D2:',tdata(i1)
            ENDIF
          ENDDO
        ENDIF

        IF (iopt.EQ.3 ) CALL LoadTriangleData(6,1,7,1,tdata,fname)  ! Dalpha
        IF (iopt.EQ.4 ) CALL LoadTriangleData(1,1,5,1,tdata,fname)  ! Ionisation
        IF (.FALSE..AND.iopt.EQ.4) THEN
          DO i1 = 1, ntri
            IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .          tri(i1)%index(2).GE.irsep)
     .        tdata(i1) = 0.0
          ENDDO
        ENDIF

        IF (iopt.EQ.5 ) CALL LoadTriangleData(1,1,10,0,tdata,fname)  ! D+ density
        IF (iopt.EQ.6 ) CALL LoadTriangleData(1,1,15,0,tdata,fname)  ! D+ temperature
        IF (iopt.EQ.7 ) THEN                                             ! D  pressure from D2 energy density
          CALL LoadTriangleData(2,1,6,1,tdata,fname) 
          fact = ECH * 0.667 ! * 7.502
          DO i1 = 1, ntri 
            tdata(i1) = tdata(i1) * fact
          ENDDO
        ENDIF
        IF (iopt.EQ.8 ) THEN                                             ! D2 pressure from D2 energy density
          CALL LoadTriangleData(3,1,6,1,tdata,fname) 
          fact = ECH * 0.667 ! * 7.502
          DO i1 = 1, ntri 
            tdata(i1) = tdata(i1) * fact
          ENDDO
        ENDIF
        IF (iopt.EQ.9 ) THEN
          CALL LoadTriangleData(3,1,7,0,tdata,fname)  ! D2 average energy
          DO i1 = 1, ntri 
            tdata(i1) = tdata(i1) * 0.667
          ENDDO
        ENDIF

        IF (iopt.EQ.10) THEN 
          ALLOCATE(tdata1(ntri))
          CALL LoadTriangleData(2,1,1,1,tdata ,fname)  ! D  density
          CALL LoadTriangleData(3,1,1,1,tdata1,fname)  ! D2 density
          tdata(1:ntri) = tdata(1:ntri) + 2.0 * tdata1(1:ntri)
          DEALLOCATE(tdata1)
        ENDIF

        IF (iopt.EQ.21) CALL LoadTriangleData(2,1,1,0,tdata,fname)  ! D  density, no volume scaling
        IF (iopt.EQ.22) CALL LoadTriangleData(3,1,1,0,tdata,fname)  ! D2 density, no volume scaling
        IF (iopt.EQ.23) CALL LoadTriangleData(6,2,6,1,tdata,fname)  ! Dgamma (total)

        IF (iopt.EQ.60) CALL LoadTriangleData(5,1,1,1,tdata,fname)  ! Balmer alpha
        IF (iopt.EQ.61) CALL LoadTriangleData(5,2,1,1,tdata,fname)  ! Lyman alpha
        IF (iopt.EQ.62) CALL LoadTriangleData(5,3,1,1,tdata,fname)  ! Lyman beta
        IF (iopt.EQ.63) CALL LoadTriangleData(5,4,1,1,tdata,fname)  ! Lyman gamma
        IF (iopt.EQ.64) CALL LoadTriangleData(5,5,1,1,tdata,fname)  ! Lyman delta
        IF (iopt.EQ.65) CALL LoadTriangleData(5,6,1,1,tdata,fname)  ! Lyman epsilon

        IF (iopt.EQ.80) CALL LoadTriangleData(2,2,1,1,tdata,fname)  ! He(1|1) density
        IF (iopt.EQ.81) CALL LoadTriangleData(2,3,1,1,tdata,fname)  ! He(2|1)
        IF (iopt.EQ.82) CALL LoadTriangleData(2,4,1,1,tdata,fname)  ! He(2|3)

        write(0,*) 'done',ZA02AS(1)-t1

c...    Load up data:
        IF (tetrahedrons) THEN 

c...      Decide if tetrahedrons are being plotted:      
          ngauge = 0
20        READ(5,'(A512)') cdum1
          IF   (cdum1(8:12).EQ.'Gauge'.OR.cdum1(8:12).EQ.'gauge'.OR.
     .          cdum1(8:12).EQ.'GAUGE') THEN
            ngauge = ngauge + 1
            READ(cdum1,*) cdum2,gauge(1:5,ngauge),gauge_tag(ngauge)  ! 1=opt, 2-4=r,z,phi, 5=radius of gauge (sphere) volume (opt=1)
            GOTO 20
          ELSE
            BACKSPACE 5
          ENDIF
          IF (ngauge.GT.0) THEN
            CALL ExtractQuantity(ngauge,gauge(1:5,1:ngauge),
     .                           gauge_tag(1:ngauge),tdata,iopt)
          ENDIF

          ntet = 1
          ALLOCATE(itet(  ntri))
          ALLOCATE(nv  (  ntet))
          ALLOCATE(rv  (4,ntet))
          ALLOCATE(zv  (4,ntet))
          CALL SpliceTetrahedrons(ntet,tet_axis,tet_value,
     .                            itet,nv,rv,zv)
          DEALLOCATE(nv  )
          DEALLOCATE(rv  )
          DEALLOCATE(zv  )
          ALLOCATE(nv  (  ntet))
          ALLOCATE(rv  (4,ntet))
          ALLOCATE(zv  (4,ntet))
          ALLOCATE(cq  (  ntet))
          CALL SpliceTetrahedrons(ntet,tet_axis,tet_value,
     .                            itet,nv,rv,zv)

          DO i1 = 1, ntet
            cq(i1) = tdata(itet(i1))
          ENDDO
          nc = ntet

          IF (tetrahedrons_dump) THEN
           CALL inOpenInterface('idl.dump_tet_'//TRIM(tet_id),ITF_WRITE)
           CALL inPutData(rv(1,1:nc),'R1','m')
           CALL inPutData(zv(1,1:nc),'Z1','m')
           CALL inPutData(rv(2,1:nc),'R2','m')
           CALL inPutData(zv(2,1:nc),'Z2','m')
           CALL inPutData(rv(3,1:nc),'R3','m')
           CALL inPutData(zv(3,1:nc),'Z3','m')
           CALL inPutData(rv(4,1:nc),'R4','m')
           CALL inPutData(zv(4,1:nc),'Z4','m')
           CALL inPutData(nv(  1:nc),'NV','N/A')
           CALL inPutData(cq(  1:nc),'VAL','N/A')
           CALL inPutData(iopt      ,'OPT','N/A')        
           CALL inPutData(tet_axis  ,'TET_AXIS','N/A')        
           CALL inPutData(tet_value ,'TET_VALUE','N/A')        
           CALL inCloseInterface
          ENDIF

        ELSE
          ALLOCATE(nv(  ntri))
          ALLOCATE(rv(3,ntri))
          ALLOCATE(zv(3,ntri))
          ALLOCATE(cq(  ntri))
          nc = 0
          nv = 0
          DO i1 = 1, ntri
            IF ((iopt.EQ.5.OR.iopt.EQ.6).AND.
     .          (tri(i1)%type.NE.MAGNETIC_GRID.OR.
     .           tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .           tri(i1)%index(2).LT.irsep)) CYCLE
            nc = nc + 1
            DO v1 = 1, 3
              nv(nc) = nv(nc) + 1
              rv(v1,nc) = SNGL(ver(tri(i1)%ver(v1),1))
              zv(v1,nc) = SNGL(ver(tri(i1)%ver(v1),2))
            ENDDO
            cq(nc) = tdata(i1)
          ENDDO
        ENDIF
c     ==================================================================
      ENDIF
c     ==================================================================




c...  Scale the data:
      IF (ALLOCATED(cq)) THEN
        IF (scalefact.EQ.-1.0) THEN
         cq = LOG10(MAX(1.0E-10,cq))
        ELSE
         cq = cq * scalefact
        ENDIF
      ENDIF

      IF (iopt.LT.90.OR.iopt.GE.500) THEN
c...    Find plotting quantity minimum and maximum:
        setqmin = .TRUE.
        setqmax = .TRUE.
        IF (qmin.NE. HI) setqmin = .FALSE.
        IF (qmax.NE.-HI) setqmax = .FALSE.
        IF (setqmin.OR.setqmax) THEN
          DO i1 = 1, nc
c...        Decide if the cell is within the viewing range:
            inside = .FALSE.
            DO i2 = 1, nv(i1)
              IF (rv(i2,i1).GE.xxmin.AND.rv(i2,i1).LE.xxmax.AND.
     .            zv(i2,i1).GE.yymin.AND.zv(i2,i1).LE.yymax) 
     .          inside=.TRUE.
            ENDDO
            IF (inside) THEN
              IF (setqmin.AND.
     .            cq(i1).NE.0.0.AND.qmin.GT.cq(i1)) qmin=cq(i1)
              IF (setqmax.AND.
     .            cq(i1).NE.0.0.AND.qmax.LT.cq(i1)) qmax=cq(i1)
            ENDIF
          ENDDO
          IF (setqmin.AND.qmin.GT.-1.0E-10) qmin = MAX(qmin,0.01*qmax)
        ENDIF
c        WRITE(0,*) 'QMIN,QMAX:',qmin,qmax

c...    Draw polygons:
        CALL PSPACE(MAP1X,MAP2X,MAP1Y,MAP2Y)
        CALL MAP   (CXMIN,CXMAX,CYMIN,CYMAX)
        CALL HSV
c         hardscale = .TRUE.
        IF (qmin.NE.qmax) THEN  ! This check is required to avoid a strange seg fault - SL, 22/01/2010
          DO i1 = 1, nc
            IF (cq(i1).LT.qmin) cq(i1) = qmin
c            IF (cq(i1).GE.qmin) THEN
c              WRITE(0,*) colouropt,i1,cq(i1),qmin,qmax
              CALL SetCol255_04(colouropt,cq(i1),qmin,qmax)
c              CALL SetCol255_04(colouropt,cq(i1),qmin,qmax)
              CALL FILCOL(255)
              CALL LINCOL(255) 
c              WRITE(6,*) 'PLOT:',rv(1,i1),zv(1,i1),cq(i1),nv(i1)
c              WRITE(6,*) '    :',rv(2,i1),zv(2,i1)
c              WRITE(6,*) '    :',rv(3,i1),zv(3,i1)
c              WRITE(6,*) '    :',rv(4,i1),zv(4,i1)
              CALL PTPLOT(rv(1,i1),zv(1,i1),1,nv(i1),1)
c            ENDIF
          ENDDO
        ENDIF

        IF (scale_set) THEN
c...      Process tags:
          IF (label(1:14).EQ.'<charge state>') 
     .     WRITE(label,'(A,I2,A)') '+',iz,' '//label(15:LEN_TRIM(label))
        ENDIF
        CALL DrawColourScale(scaleopt,colouropt,qmin,qmax,label)

        IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
        DEALLOCATE(nv)
        DEALLOCATE(rv)
        DEALLOCATE(zv)
        DEALLOCATE(cq)
      ENDIF



c...  Draw Vessel and grid outline:

c     ==================================================================
      IF (iopt.GE.500) THEN
c     ==================================================================
c...    Magnetic grid:
        nline = 0
        ALLOCATE(lines2(4*MAXNKS*MAXNRS,2))
        ALLOCATE(lcolour(4*MAXNKS*MAXNRS))

        nver = 0
        CALL ALLOC_VERTEX(20*MAXNKS) ! Borrowed from mod_triangle

        DO ir = 1, nrs
          IF (idring(ir).EQ.BOUNDARY) CYCLE
          DO ik = 1, nks(ir)

            id = korpg(ik,ir)

            IF     (idring(irouts(ik,ir)).EQ.BOUNDARY.OR.
     .              irouts(ik,ir).EQ.ir) THEN  ! The connection map for CDN is not working...
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(2,id))
              ver(nver,2) = DBLE(zvertp(2,id))
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(3,id))
              ver(nver,2) = DBLE(zvertp(3,id))
              nline = nline + 1
              lines2(nline,1) = nver - 1
              lines2(nline,2) = nver 
              lcolour(nline) = ncols + 1
            ELSEIF (idring(irins(ik,ir)).EQ.BOUNDARY.OR.ir.EQ.irsep.OR.
     .            (double_null.AND.
     .             ir.EQ.irouts(1                 ,MAX(1,irsep2)).OR.
     .             ir.EQ.irouts(nks(MAX(1,irsep2)),MAX(1,irsep2)))) THEN
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(1,id))
              ver(nver,2) = DBLE(zvertp(1,id))
              nver = nver + 1
              ver(nver,1) = DBLE(rvertp(4,id))
              ver(nver,2) = DBLE(zvertp(4,id))
              nline = nline + 1
              lines2(nline,1) = nver - 1
              lines2(nline,2) = nver 
              IF (idring(irins(ik,ir)).EQ.BOUNDARY) THEN
                lcolour(nline) = ncols + 1
              ELSE
                lcolour(nline) = ncols + 3
              ENDIF
            ENDIF

          ENDDO

c...      Targets:
          IF (ir.GT.irsep) THEN
            id = korpg(1,ir)
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(1,id))
            ver(nver,2) = DBLE(zvertp(1,id))
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(2,id))
            ver(nver,2) = DBLE(zvertp(2,id))
            nline = nline + 1
            lines2(nline,1) = nver - 1
            lines2(nline,2) = nver 
            lcolour(nline) = ncols + 1

            id = korpg(nks(ir),ir)
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(3,id))
            ver(nver,2) = DBLE(zvertp(3,id))
            nver = nver + 1
            ver(nver,1) = DBLE(rvertp(4,id))
            ver(nver,2) = DBLE(zvertp(4,id))
            nline = nline + 1
            lines2(nline,1) = nver - 1
            lines2(nline,2) = nver 
            lcolour(nline) = ncols + 1
          ENDIF

        ENDDO

c...    Wall:
        DO i1 = 1, wallpts
          IF (wallpt(i1,18).NE.0.0) CYCLE

          nver = nver + 1
          ver(nver,1) = DBLE(wallpt(i1,20))
          ver(nver,2) = DBLE(wallpt(i1,21))
          nver = nver + 1
          ver(nver,1) = DBLE(wallpt(i1,22))
          ver(nver,2) = DBLE(wallpt(i1,23))
          nline = nline + 1
          lines2(nline,1) = nver - 1
          lines2(nline,2) = nver 
          lcolour(nline) = 1
        ENDDO
c     ==================================================================
      ELSEIF (iopt.LT.500.AND..NOT.tetrahedrons) THEN
c     ==================================================================
        ALLOCATE(lines2 (3*ntri,2))
        ALLOCATE(lcolour(3*ntri))
        nline = 0
        DO i1 = 1, ntri
          DO v1 = 1, 3
            v2 = v1 + 1
            IF (v1.EQ.3) v2 = 1
c *TEMP*
            IF     (iopt.EQ.95) THEN

              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.     ! Walls
     .             tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.   ! Targets
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.96) THEN
c             All triangles:
              nline = nline + 1
              lcolour(nline) = ncols + 1
              lines2(nline,1)=tri(i1)%ver(v1)
              lines2(nline,2)=tri(i1)%ver(v2)
            ELSEIF (iopt.EQ.97) THEN
              IF (i1.EQ.6268.OR.i1.EQ.6416) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSEIF (iopt.EQ.98) THEN
c              IF (i1.NE.1.AND.i1.NE.380.AND.
c     .            i1.NE.423.AND.i1.NE.448.AND.i1.NE.450) CYCLE
c              IF (tri(i1)%type.EQ.VACUUM_GRID) WRITE(0,*) ' *** VAC!'
              IF (tri(i1)%type.EQ.VACUUM_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF
              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF
            ELSEIF (iopt.EQ.99) THEN
c...          Magnetic grid triangles + wall surfaces:
              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            v1.NE.1) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

            ELSE
              IF (tri(i1)%type.EQ.VACUUM_GRID.AND.
c     .            (tri(i1)%sur(v1).EQ.4.OR.
c     .             tri(i1)%sur(v1).GT.10)) THEN
     .             tri(i1)%sur(v1).NE.0) THEN

                nline = nline + 1
                lcolour(nline) = 1  ! ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%sur(v1).NE.0) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 3
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

              IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .            tri(i1)%index(2).EQ.irsep.AND.
     .            tri(i1)%sideindex(1,v1).EQ.14) THEN
                nline = nline + 1
                lcolour(nline) = ncols + 1
                lines2(nline,1)=tri(i1)%ver(v1)
                lines2(nline,2)=tri(i1)%ver(v2)
              ENDIF

           ENDIF
c            nline = nline + 1
c            lines2(nline,1)=tri(i1)%ver(v1)
c            lines2(nline,2)=tri(i1)%ver(v2)
          ENDDO
        ENDDO
c...    Remove duplicates:
        DO i1 = 1, nline-1
          DO i2 = i1+1, nline
            IF (lines2(i1,1).NE.-999.0.AND.lines2(i2,1).NE.-999.0) THEN
c              IF (lines2(i1,1).LE.0.OR.lines2(i1,2).LE.0.OR.
c     .            lines2(i2,1).LE.0.OR.lines2(i2,2).LE.0) CYCLE
c              IF (lines2(i1,1).GT.288.OR.lines2(i1,2).GT.288.OR.  
c     .            lines2(i2,1).GT.288.OR.lines2(i2,2).GT.288) CYCLE

              IF
c * IMPROVE THIS CHECK! *
     .          ((DABS(ver(lines2(i1,1),1)-
     .                 ver(lines2(i2,2),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,1),2)-
     .                 ver(lines2(i2,2),2)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),1)-
     .                 ver(lines2(i2,1),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),2)-
     .                 ver(lines2(i2,1),2)).LT.DTOL).OR.
     .           (DABS(ver(lines2(i1,1),1)-
     .                 ver(lines2(i2,1),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,1),2)-
     .                 ver(lines2(i2,1),2)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),1)-
     .                 ver(lines2(i2,2),1)).LT.DTOL.AND.
     .            DABS(ver(lines2(i1,2),2)-
     .                 ver(lines2(i2,2),2)).LT.DTOL)) THEN
                IF (lcolour(i1).EQ.1) THEN
                  lines2(i1,1) = -999.0
                ELSE
                  lines2(i2,1) = -999.0
                ENDIF
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        DO i1 = nline, 1, -1
          IF (lines2(i1,1).EQ.-999.0) THEN
c            WRITE(0,*) 'DELETING:',i1
            DO i2 = i1, nline-1
              lines2(i2,1) = lines2(i2+1,1)
              lines2(i2,2) = lines2(i2+1,2)
              lcolour(i2) = lcolour(i2+1)
            ENDDO
            nline = nline - 1
          ENDIF
        ENDDO
c     ==================================================================
      ENDIF
c     ==================================================================



      IF (iopt.EQ.95) THEN
        CALL DrawGrid(iopt)        
      ELSE
c...    Plot polygons:
        CALL PSPACE (map1x,map2x,map1y,map2y)      
        CALL MAP    (cxmin,cxmax,cymin,cymax)
        IF (.NOT.tetrahedrons) THEN
          lastcolour = -1
          DO i1 = 1, nline
            IF (lastcolour.NE.lcolour(i1)) THEN
              CALL LINCOL(lcolour(i1)) 
              lastcolour = lcolour(i1)
            ENDIF

c            IF (lines2(i1,1).LE.0.OR.lines2(i1,2).LE.0.OR.
c     .          lines2(i2,1).LE.0.OR.lines2(i2,2).LE.0) CYCLE
c            IF (lines2(i1,1).GT.288.OR.lines2(i1,2).GT.288.OR. 
c     .          lines2(i2,1).GT.288.OR.lines2(i2,2).GT.288) CYCLE

            CALL POSITN(SNGL(ver(lines2(i1,1),1)),
     .                  SNGL(ver(lines2(i1,1),2)))
            CALL JOIN  (SNGL(ver(lines2(i1,2),1)),
     .                  SNGL(ver(lines2(i1,2),2)))
          ENDDO
        ENDIF
c...    Frame:
        CALL DrawFrame
      ENDIF

c...     Draw the views on the LOS plot:
c         CALL PSPACE (map1x,map2x,map1y,map2y)
c         CALL MAP    (cxmin,cxmax,0.0,1.0)
c         CALL BROKEN(6,6,6,6)
c         CALL LinCol(1)
c         DO i2 = 1, nline
c           CALL POSITN (lines2(i2),0.0)
c           CALL JOIN   (lines2(i2),1.0)        
c         ENDDO

c...  Add a caption to the plot:
      READ(5,'(A256)') cdum1
      IF   (cdum1(8:14).EQ.'Noframe'.OR.cdum1(8:12).EQ.'noframe'.OR.
     .      cdum1(8:14).EQ.'NOFRAME') THEN
        numplots = numplots + 1
        reset_origin = .FALSE.
      ELSE
        numplots = 0
        BACKSPACE 5
        CALL FRAME
        reset_origin = .TRUE.
      ENDIF

c...  Clear arrays:
c      CALL DEALLOC_ALL  
      CALL DEALLOC_VERTEX  
      CALL DEALLOC_SURFACE   
      CALL DEALLOC_TRIANGLE

      IF (ALLOCATED(itet   )) DEALLOCATE(itet   )
      IF (ALLOCATED(gdata1 )) DEALLOCATE(gdata1 )
      IF (ALLOCATED(lines2 )) DEALLOCATE(lines2 )
      IF (ALLOCATED(lcolour)) DEALLOCATE(lcolour)

      RETURN
99    WRITE(0,*) 'IOPT =',iopt
      STOP
      END
