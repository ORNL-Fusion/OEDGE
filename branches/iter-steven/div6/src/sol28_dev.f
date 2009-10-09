c     -*-Fortran-*-
c
c ====================================================================
c
      SUBROUTINE GetVertex(iobj,i,x,y)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN ) :: iobj,i
      REAL*8 , INTENT(OUT) :: x,y

      INTEGER isrf,ivtx

c...  The first vertex of each side corresponds to i:      
      
      isrf = obj(iobj)%iside(i)

      IF (srf(ABS(isrf))%type.NE.SPR_LINE_SEGMENT) 
     .  CALL ER('GetVertex','Routine only handles line segments',*99)

      IF (isrf.LT.0) THEN
        isrf = -isrf
        ivtx = srf(isrf)%ivtx(srf(isrf)%nvtx)
        x = vtx(1,ivtx)
        y = vtx(2,ivtx)
      ELSE
        ivtx = srf(isrf)%ivtx(1)
        x = vtx(1,ivtx)
        y = vtx(2,ivtx)
      ENDIF

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ClipWallToGrid(nwall,rwall,zwall,MAXNSEG)
      USE mod_geometry
      USE mod_sol28_global
      IMPLICIT none

      INTEGER, INTENT(IN)    :: MAXNSEG
      INTEGER, INTENT(INOUT) :: nwall
      REAL   , INTENT(INOUT) :: rwall(MAXNSEG,2),zwall(MAXNSEG,2)

      INTEGER GetTube       
 
      INTEGER, PARAMETER :: MAXNLIST = 1000
      REAL*8 , PARAMETER :: DTOL = 1.0D-07

      INTEGER fp,iobj,itube,nlist,ilist(MAXNLIST,2),clist(MAXNLIST,2),
     .        tube_set,i1,i2,i3,swall(nwall),iwall,mlist(MAXNLIST)
      LOGICAL debug
      REAL*8  x1,x2,y1,y2,x3,x4,y3,y4,s12,s34,s12max,
     .        xlist(MAXNLIST,2),ylist(MAXNLIST,2),store_x2,store_y2

      fp = 88
      debug = .TRUE.

      CALL DumpData_OSM('output.clipping','Trying to clip grid')

c...  Assume wall is clockwise-specified.  Cut the grid by extending 
c     the poloidal surfaces of end cells that lie on the radial 
c     fluid grid boundary surfaces.  Delete everything that's between
c     cuts, as appropriate:

      tube_set = -1
      nlist = 0

      DO iobj = 1, nobj
        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID.OR.  ! Only for fluid grid objects
     .      grp(obj(iobj)%group)%type  .NE.GRP_QUADRANGLE.OR.
     .      GetTube(iobj,IND_OBJECT).LT.grid%isep) CYCLE

        IF (obj(iobj)%omap(2).EQ.-1.OR.obj(iobj)%omap(4).EQ.-1) THEN
          itube = GetTube(iobj,IND_OBJECT)
          IF     (tube_set.NE.itube) THEN
            tube_set = itube         ! There's an assumption here that 
            nlist = nlist + 1        ! the cells on a ring are sequential
            ilist(nlist,1:2) = iobj  ! which shoudl be fine...
          ELSEIF (tube_set.EQ.itube) THEN
            ilist(nlist,2  ) = iobj
          ENDIF
        ENDIF
      ENDDO

      IF (debug) THEN
        DO i1 = 1, nlist
          WRITE(fp,*) 'CLIP LIST:',ilist(i1,1:2)
        ENDDO
      ENDIF

c...  Collect the cuts:
      clist = 0
      mlist = 0
      xlist = 0.0D0
      ylist = 0.0D0
      i1 = 0
      DO WHILE(i1.LT.nlist)
        i1 = i1 + 1
        DO i2 = 1, 2
          iobj = ilist(i1,i2)
          IF (obj(iobj)%omap(2).EQ.-1) THEN
            mlist(i1) = 1
            IF (i2.EQ.1) THEN
              CALL GetVertex(iobj,3,x1,y1)
              CALL GetVertex(iobj,2,x2,y2)
            ELSE
              CALL GetVertex(iobj,2,x1,y1)
              CALL GetVertex(iobj,3,x2,y2)
            ENDIF
          ELSE
            mlist(i1) = 2
            IF (i2.EQ.1) THEN
              CALL GetVertex(iobj,4,x1,y1)
              CALL GetVertex(iobj,1,x2,y2)
            ELSE
              CALL GetVertex(iobj,1,x1,y1)
              CALL GetVertex(iobj,4,x2,y2)
            ENDIF
          ENDIF
c         Increase the length of the line segment in case there's
c         a nominal mismatch between the wall and target 
c         specifications or if the line segment is very short:
          store_x2 = x2
          store_y2 = y2
          x2 = x1 + 100.0D0 * (x2 - x1)
          y2 = y1 + 100.0D0 * (y2 - y1)

          IF (debug) THEN
            WRITE(fp,*) ' --------------------'
            WRITE(fp,*) ' X,Y1=',x1,y1
            WRITE(fp,*) ' X,Y2=',x2,y2
          ENDIF
c         Search the wall for intersections:
          s12max = 1.0D+10
          DO iwall = 1, nwall
            x3 = DBLE(rwall(iwall,1))
            y3 = DBLE(zwall(iwall,1))
            x4 = DBLE(rwall(iwall,2))
            y4 = DBLE(zwall(iwall,2))
            CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
            IF (debug) THEN
              WRITE(fp,*) '  CALCINTER:',i1,i2,iwall
              WRITE(fp,*) '  CALCINTER:',s12,s34
              WRITE(fp,*) '  CALCINTER:',x3,y3
              WRITE(fp,*) '  CALCINTER:',x4,y4
            ENDIF
            IF (s12.GT.0.0D0.AND.s12.LT.1.0D0.AND.
     .          s34.GT.0.0D0.AND.s34.LT.1.0D0.AND.
     .          s12.LT.s12max) THEN
              s12max = s12
              clist(i1,i2) = iwall
              xlist(i1,i2) = store_x2
              ylist(i1,i2) = store_y2
              IF (debug) WRITE(fp,*) '  *** CUT ***',i1,i2,s12,iwall
            ENDIF
          ENDDO
          IF (clist(i1,i2).EQ.0) THEN
c...        Problem with this cut pair, so delete them from the 
c           list (have to complete the wall by hand at the moment):            
            IF (debug) WRITE(fp,*) ' CUT NOT FOUND, DELETING CUT',i1
            DO i3 = i1, nlist-1
              ilist(i3,:) = ilist(i3+1,:)
            ENDDO
            i1 = i1 - 1
            nlist = nlist - 1
            EXIT
          ENDIF
        ENDDO
      ENDDO

      IF (debug) THEN
        DO i1 = 1, nlist
          WRITE(fp,'(A,4I6,2X,4F14.7)') 
     .      'CLIP LIST:',ilist(i1,:),clist(i1,:),xlist(i1,:),ylist(i1,:)
        ENDDO
      ENDIF

c...  Check if a line segment is cut more than once at either end:
      DO i2 = 1, 2       
        swall = 0
        DO i1 = 1, nlist
          IF (swall(clist(i1,i2)).EQ.0) THEN
            swall(clist(i1,i2)) = 1
          ELSE
            CALL ER('ClipWallToGrid','Wall segment cut at the same '//
     .              'end more than once',*99)
          ENDIF
        ENDDO
      ENDDO

      swall = 0
      DO i1 = 1, nlist
c       Overwrite the approriate end vertices in the wall segments:
        IF (mlist(i1).EQ.1) THEN
          rwall(clist(i1,1),1) = SNGL(xlist(i1,1))
          zwall(clist(i1,1),1) = SNGL(ylist(i1,1))
          rwall(clist(i1,2),2) = SNGL(xlist(i1,2))
          zwall(clist(i1,2),2) = SNGL(ylist(i1,2))
        ELSE
          rwall(clist(i1,1),2) = SNGL(xlist(i1,1))
          zwall(clist(i1,1),2) = SNGL(ylist(i1,1))
          rwall(clist(i1,2),1) = SNGL(xlist(i1,2))
          zwall(clist(i1,2),1) = SNGL(ylist(i1,2))
        ENDIF
c       Select segments that aren't "behind" the targets, and so won't
c       be deleted:
        IF (mlist(i1).EQ.1) THEN
          iwall = clist(i1,1)-1           ! Unlikely that clist(,2) = clist(,1) - 1, 
        ELSE                              ! but if it happens there will be trouble...
          iwall = clist(i1,1)+1           
        ENDIF
        DO WHILE(iwall.NE.clist(i1,2))  
          IF (mlist(i1).EQ.1) THEN
            iwall = iwall + 1 
            IF (iwall.EQ.nwall+1) iwall = 1
          ELSE
            iwall = iwall - 1 
            IF (iwall.EQ.0) iwall = nwall
          ENDIF
          swall(iwall) = 1
          IF (debug) WRITE(fp,*) 'SAVING:',i1,iwall
        ENDDO
      ENDDO

c...  Delete wall segments:
      DO iwall = nwall, 1, -1
        IF (swall(iwall).EQ.1) CYCLE
        IF (debug) WRITE(fp,*) 'DELETING:',iwall 
        DO i1 = iwall, nwall-1
          rwall(i1,:) = rwall(i1+1,:)
          zwall(i1,:) = zwall(i1+1,:)
          swall(i1  ) = swall(i1+1  )
        ENDDO
        nwall = nwall - 1
      ENDDO

      IF (debug) THEN
        DO iwall = 1, nwall
          WRITE(fp,'(A,2F14.7,2X,2F14.7)') 
     .      'WALL:',rwall(iwall,1),zwall(iwall,1),
     .              rwall(iwall,2),zwall(iwall,2)
        ENDDO
      ENDIF


      RETURN
 99   WRITE(fp,*) 'I1   =',i1
      WRITE(fp,*) 'I2   =',i2
      WRITE(fp,*) 'S12  =',s12
      WRITE(fp,*) 'IWALL=',iwall
      STOP
      END
c
c ====================================================================
c
      SUBROUTINE OutputVoidSpecification(fp)
      USE mod_sol28_global
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: fp

      INTEGER ivoid

      WRITE(fp,*)
      WRITE(fp,'(A)') 'EIRENE TRIANGLE GRID SETTINGS:'
      WRITE(fp,'(A4,3(2X,A8),A10,2X,A18,2X,A4,2X,3A10)') 
     .  'zone','grid1,2','wall1,2','add1,2','res','hole','code',
     .  'ne','Te','Ti'
      DO ivoid = 1, opt_eir%nvoid
        WRITE(fp,'(I4,3(2X,2I4),F10.4,2X,2F9.4,2X,I4,2X,1P,3E10.2,0P)') 
     .    opt_eir%void_zone(    ivoid),
     .    opt_eir%void_grid(1:2,ivoid),
     .    opt_eir%void_wall(1:2,ivoid),
     .    opt_eir%void_add (1:2,ivoid),
     .    opt_eir%void_res (    ivoid),
     .    opt_eir%void_hole(1:2,ivoid),
     .    opt_eir%void_code(    ivoid),
     .    opt_eir%void_ne  (    ivoid),
     .    opt_eir%void_te  (    ivoid),
     .    opt_eir%void_ti  (    ivoid)
      ENDDO

      RETURN
 99   STOP
      END


c
c ======================================================================
c
      LOGICAL FUNCTION PointInVoid(x1,y1,nseg,seg,MAXNSEG,
     .                             npts,pts,MAXNPTS,debug,fp)
      IMPLICIT none

      INTEGER, INTENT(IN) :: MAXNPTS,MAXNSEG,nseg,npts,seg(0:MAXNSEG,4),
     .                       fp
      LOGICAL, INTENT(IN) :: debug
      REAL*8 , INTENT(IN) :: x1,y1,pts(MAXNPTS,2)

      REAL*8 , PARAMETER :: DTOL=1.0D-07

      INTEGER iseg,ninter
      REAL*8  x2,x3,x4,y2,y3,y4,s12,s34

      PointInVoid = .FALSE.

      x2 = x1 + 20.0D0
      y2 = y1
      ninter = 0
      DO iseg = 1, nseg
        x3 = pts(seg(iseg,1),1)
        y3 = pts(seg(iseg,1),2)
        x4 = pts(seg(iseg,2),1)
        y4 = pts(seg(iseg,2),2)
        CALL CalcInter(x1,y1,x2,y2,x3,y3,x4,y4,s12,s34) 
        IF (s12.GT.DTOL.AND.s34.GT.0.0D0.AND.s34.LT.1.0D0) 
     .    ninter = ninter + 1
        IF (debug) WRITE(fp,'(4X,A,2F14.7,I4,2F12.5)')
     .      '           :',s12,s34,ninter,x1,y1
      ENDDO  

      IF (ninter.GT.0.AND.MOD(ninter+1,2).EQ.0) PointInVoid = .TRUE.

      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE SetupVoidProcessing
      USE mod_sol28_global
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,ivoid,isrf,i1,i2,idefault
      LOGICAL void_specified

      TYPE(type_options_eirene) :: opt_tmp

c...  Output:
      CALL OutputVoidSpecification(eirfp)

      CALL ProcessVoid(-999)  ! Initialisation

c      default_index = 0     
c      DO ivoid = 1, opt_eir%nvoid
c        IF (opt_eir%void_zone.EQ.-1) default_index = ivoid
c      ENDDO

c...  Assign defaults, if necessary:
      idefault  = 0
      void_specified = .FALSE.
      DO i1 = 1, opt_eir%nvoid
        IF (opt_eir%void_zone(i1).EQ.-1) idefault = i1
        IF (opt_eir%void_zone(i1).GE. 1) void_specified = .TRUE.
      ENDDO
      IF (void_specified.AND.idefault.NE.0) 
     .  CALL ER('SetupVoidProcessing','Cannot specify default void '//
     .          'setup and also set specific void parameters',*99)
      IF (idefault.NE.0) THEN
        i1 = idefault
        opt_tmp = opt_eir
        DO isrf = 1, nsurface
          IF (surface(isrf)%type    .NE.NON_DEFAULT_STANDARD  .OR.
     .        surface(isrf)%subtype .NE.MAGNETIC_GRID_BOUNDARY.OR.
     .        surface(isrf)%index(6).LT.opt_tmp%void_grid(1,i1)) CYCLE
          opt_eir%nvoid = opt_eir%nvoid + 1
          opt_eir%void_zone(  opt_eir%nvoid) = opt_eir%nvoid - 
     .                                         opt_tmp%nvoid 
          opt_eir%void_grid(:,opt_eir%nvoid) = surface(isrf)%index(6)
          opt_eir%void_wall(:,opt_eir%nvoid) = opt_tmp%void_wall(:,i1)
          opt_eir%void_add (:,opt_eir%nvoid) = opt_tmp%void_add (:,i1)
          opt_eir%void_res (  opt_eir%nvoid) = opt_tmp%void_res (  i1)
          opt_eir%void_hole(:,opt_eir%nvoid) = opt_tmp%void_hole(:,i1)
          opt_eir%void_code(  opt_eir%nvoid) = opt_tmp%void_code(  i1)
          opt_eir%void_ne  (  opt_eir%nvoid) = opt_tmp%void_ne  (  i1)
          opt_eir%void_te  (  opt_eir%nvoid) = opt_tmp%void_te  (  i1)
          opt_eir%void_ti  (  opt_eir%nvoid) = opt_tmp%void_ti  (  i1)
        ENDDO
      ENDIF

c...  Output:
      CALL OutputVoidSpecification(eirfp)

c...  Check void blacklist:
      DO i1 = 1, opt_eir%nvoid
        IF (opt_eir%void_zone(i1).NE.-2) CYCLE
        DO i2 = 1, opt_eir%nvoid
          IF (i1.EQ.i2) CYCLE
          IF (opt_eir%void_grid(1,i1).LE.opt_eir%void_zone(i2).AND.
     .        opt_eir%void_grid(2,i1).GE.opt_eir%void_zone(i2)) 
     .      opt_eir%void_zone(i2) = -999  ! Delete
        ENDDO
      ENDDO

c...  Output:
      CALL OutputVoidSpecification(eirfp)

c...  Delete tagged void regions:
      DO i1 = opt_eir%nvoid, 1, -1
        IF (opt_eir%void_zone(i1).NE.-999) CYCLE
        DO i2 = i1, opt_eir%nvoid-1
          opt_eir%void_zone(  i2) = opt_eir%void_zone(  i2+1)-1
          opt_eir%void_grid(:,i2) = opt_eir%void_grid(:,i2+1)
          opt_eir%void_wall(:,i2) = opt_eir%void_wall(:,i2+1)
          opt_eir%void_add (:,i2) = opt_eir%void_add (:,i2+1)
          opt_eir%void_res (  i2) = opt_eir%void_res (  i2+1)
          opt_eir%void_hole(:,i2) = opt_eir%void_hole(:,i2+1)
          opt_eir%void_code(  i2) = opt_eir%void_code(  i2+1)
          opt_eir%void_ne  (  i2) = opt_eir%void_ne  (  i2+1)
          opt_eir%void_te  (  i2) = opt_eir%void_te  (  i2+1)
          opt_eir%void_ti  (  i2) = opt_eir%void_ti  (  i2+1)
        ENDDO
        opt_eir%nvoid = opt_eir%nvoid - 1
      ENDDO

c...  Output:
      CALL OutputVoidSpecification(eirfp)


      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE ProcessVoid(izone)
      USE mod_sol28_global
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER, INTENT(IN) :: izone

      LOGICAL PointInVoid

      INTEGER, PARAMETER :: MAXNSEG = 10000, MAXNPTS = 20000
      REAL*8 , PARAMETER :: DTOL=1.0D-07

      INTEGER   fp,ivoid,isrf,isrf1,isrf2,i1,i2,itri,v1,v2,code,
     .          nseg,seg(0:MAXNSEG,4),icnt,nhole,npts,pass,
     .          i3,i4,i5,iseg1,iseg2,ilink,tmp_nseg,ninter
      LOGICAL   debug,cont,link
      CHARACTER command*512
      REAL      area,ne,te,ti,res(nsurface)
      REAL*8    x1,x2,y1,y2,len,t,tstep,xhole(50),yhole(50),
     .          pts(MAXNPTS,2)

      SAVE


      IF (izone.LT.0) THEN
        IF (izone.EQ.-999) res = 0.0  ! Initialisation
        RETURN
      ENDIF

      debug = .TRUE.
      fp = 88

      nseg = 0
      seg  = 0
      npts = 0       
      nhole = 0
      area = 0.0
      ne = 0.0
      te = 0.0
      ti = 0.0

      IF (debug) WRITE(fp,*) 'HERE IN ASSEMBLE VOID',izone

      DO ivoid = 1, opt_eir%nvoid
        IF (opt_eir%void_zone(ivoid).NE.izone) CYCLE

        IF (debug) WRITE(fp,*) '  PROCESSING VOID SETUP',ivoid
c       ----------------------------------------------------------------
c...    Examine the outer radial boundary surfaces of the fluid grid and 
c       collect the associated line segments:
c
c       Map the surface indices provided in the input file to the surface
c       indices defined in the EIRENE setup code:
        isrf1 = -1
        isrf2 = -1
        i1 = opt_eir%void_grid(1,ivoid)
        i2 = opt_eir%void_grid(2,ivoid)
        DO isrf = 1, nsurface
c          IF (debug) WRITE(fp,*) 'GRID SURFACE:',isrf
          IF (surface(isrf)%type   .NE.NON_DEFAULT_STANDARD  .OR.
     .        surface(isrf)%subtype.NE.MAGNETIC_GRID_BOUNDARY) CYCLE
          IF (debug) WRITE(fp,*) 'GRID SURFACE: OK',isrf,
     .                            surface(isrf)%index(6)
          IF (isrf1.EQ.-1.AND.surface(isrf)%index(6).GE.i1) isrf1 = isrf
          IF (                surface(isrf)%index(6).LE.i2) isrf2 = isrf
        ENDDO

        IF (debug) WRITE(fp,*) 'GRID:',ivoid,i1,i2,isrf1,isrf2

c       The boundary surface indices have been identified, so search the
c       list of triangles for sides that match up with these surfaces:
        IF (isrf1.NE.-1.AND.isrf2.NE.-1) THEN
          DO itri = 1, ntri
            DO v1 = 1, 3
              isrf = tri(itri)%sur(v1)
              IF (isrf.EQ.0) CYCLE
              IF (surface(isrf)%type    .EQ.NON_DEFAULT_STANDARD  .AND.
     .            surface(isrf)%subtype .EQ.MAGNETIC_GRID_BOUNDARY.AND.
     .            surface(isrf)%index(6).GT.0.AND. 
     .            isrf.GE.isrf1.AND.isrf.LE.isrf2) THEN
                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                v2 = v1 + 1
                IF (v1.EQ.3) v2 = 1
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v1),1) 
                pts(npts,2) = ver(tri(itri)%ver(v1),2)
                npts = npts + 1
                pts(npts,1) = ver(tri(itri)%ver(v2),1)
                pts(npts,2) = ver(tri(itri)%ver(v2),2)
                IF (debug) THEN
                  WRITE(fp,'(A,5I6)') 'GRID:',nseg,itri,isrf,isrf1,isrf2
                  WRITE(fp,*) '    PTS1=',pts(npts-1,1:2)
                  WRITE(fp,*) '    PTS2=',pts(npts  ,1:2)
                ENDIF
              ENDIF
            ENDDO              
          ENDDO      
        ENDIF

c       ------------------------------------------------------------------
c...    Search through the list of standard wall line segments and build
c       a list for each zone that completes the individual voids:
        i1 = opt_eir%void_wall(1,ivoid)
        i2 = opt_eir%void_wall(2,ivoid)

        IF     (i1.EQ.-1.AND.i2.EQ.-1) THEN
c         Find the wall segments for this zone automatically -- just keep
c         mindlessly filing through the wall segments until the path is 
c         closed:
          cont = .TRUE.
          pass = 0
          DO WHILE(cont)
            pass = pass + 1
            IF (pass.EQ.100) STOP 'NOT PASSING...'
            cont = .FALSE.
c           Check if there is already a link to this segment:
            tmp_nseg = nseg
            DO iseg1 = 1, tmp_nseg
              IF (seg(iseg1,3).EQ.1) CYCLE
c             Check both ends of the current focus segment:
              DO ilink = 1, 2
                link = .FALSE.
                DO iseg2 = 1, nseg
                  IF (iseg1.EQ.iseg2) CYCLE
                  i3 = seg(iseg1,ilink)
c                 Check both ends of the test segment:
                  i4 = seg(iseg2,1)
                  i5 = seg(iseg2,2)
                  IF ((DABS(pts(i3,1)-pts(i4,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i4,2)).LT.DTOL).OR.
     .                (DABS(pts(i3,1)-pts(i5,1)).LT.DTOL.AND.
     .                 DABS(pts(i3,2)-pts(i5,2)).LT.DTOL)) THEN
                    link = .TRUE.
                    IF (ilink.EQ.2) seg(iseg1,3) = 1
                    EXIT
                  ENDIF
                ENDDO
                IF (debug) WRITE(fp,*) '  PASS:',pass,iseg1,link
                IF (.NOT.link) EXIT
              ENDDO

c             Both ends of the segment are attached so start
c             looking at the next segment:
              IF (link) CYCLE

              cont = .TRUE.
              DO isrf = 1, nsurface
                IF (surface(isrf)%type.NE.VESSEL_WALL.OR.    
     .              seg(iseg1,4).EQ.isrf) CYCLE
                WRITE(fp,*) 'WTF:',res(74)
                x1 = surface(isrf)%v(1,1)
                y1 = surface(isrf)%v(2,1)
                x2 = surface(isrf)%v(1,2)
                y2 = surface(isrf)%v(2,2)
                IF ((DABS(pts(i3,1)-x1).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y1).LT.DTOL).OR.
     .              (DABS(pts(i3,1)-x2).LT.DTOL.AND.
     .               DABS(pts(i3,2)-y2).LT.DTOL)) THEN

                  IF (res(isrf).EQ.0.0) THEN
                    res(isrf) = opt_eir%void_res(ivoid)
                  ELSE
                    WRITE(fp,*) 'ISRF,RES:',isrf,nsurface,res(isrf)
                    WRITE(fp,*) 'X1,Y1:',x1,y1
                    WRITE(fp,*) 'X2,Y2:',x2,y2
                    STOP 'NOT READY FOR THE RES...'
                  ENDIF

                  len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
                  IF (len.GT.DBLE(res(isrf))) THEN
                    tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
                  ELSE
                    tstep = 1.0D0
                  ENDIF

                  IF (debug) THEN
                    WRITE(fp,'(A,2I6,2X,2I6,2F10.4)') 
     .                ' NEW WALL SEG:',iseg1,ilink, 
     .                surface(isrf)%index(1:2),res(isrf),tstep

                    WRITE(fp,*) '    I3    =',i3
                    WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
                    WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
                    WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
                    WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
                    WRITE(fp,*) '    X1,Y1:',x1,y1
                    WRITE(fp,*) '    X2,Y2:',x2,y2

c                    STOP 'dfsdfsd'
                  ENDIF
              
c                  WRITE(eirfp,*) x1,y1
c                  WRITE(eirfp,*) x2,y2
                  DO t = 0.0D0, 0.9999999D0, tstep 
                    nseg = nseg + 1
                    seg(nseg,1) = npts + 1
                    seg(nseg,2) = npts + 2
                    seg(nseg,4) = isrf
                    npts = npts + 1
                    pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                    pts(npts,2) = y1 + t * (y2 - y1) 
                    npts = npts + 1
                    pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                    pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
                  ENDDO
                  EXIT
                ENDIF
              ENDDO
              IF (isrf.EQ.nsurface+1) 
     .          CALL ER('ProcessVoid','No link to wall surface',*99)
            ENDDO

          ENDDO
        ELSEIF (i1.GT.0.AND.i2.GT.0) THEN
c         Select the wall segments based on the index values in I1 and I2:
c            DO i3 = 1, nsurface
c              IF (surface(i3)%type.EQ.VESSEL_WALL.AND.
c     .            ((eirtri(i2,2).EQ.2.0.AND.
c     .              surface(i3)%index(1).GE.NINT(eirtri(i2,3)).AND.
c     .              surface(i3)%index(1).LE.NINT(eirtri(i2,4))).OR.
c     .             (eirtri(i2,2).EQ.3.0.AND.
c     .              surface(i3)%index(2).GE.NINT(eirtri(i2,3)).AND.
c     .              surface(i3)%index(2).LE.NINT(eirtri(i2,4))))) THEN
          STOP 'WORKING ON IT'
        ENDIF
      ENDDO


c...  Sort segments:
      DO i2 = 1, nseg-1
        DO i3 = i2+1, nseg
          IF     (DABS(pts(seg(i2,2),1)-pts(seg(i3,1),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,1),2)).LT.DTOL) THEN
            IF (i3.EQ.i2+1) THEN
c...          Do nothing, all okay:
              EXIT
            ELSE
              seg(0   ,1:4) = seg(i2+1,1:4)
              seg(i2+1,1:4) = seg(i3  ,1:4)
              seg(i3  ,1:4) = seg(0   ,1:4)
              EXIT
            ENDIF
          ELSEIF (DABS(pts(seg(i2,2),1)-pts(seg(i3,2),1)).LT.DTOL.AND.
     .            DABS(pts(seg(i2,2),2)-pts(seg(i3,2),2)).LT.DTOL) THEN
            IF (i3.EQ.i2+1) THEN
              seg(0 ,1) = seg(i3,1)
              seg(i3,1) = seg(i3,2) ! Swap the order of the points
              seg(i3,2) = seg(0 ,1)
              EXIT
            ELSE
              seg(0   ,1:4) = seg(i2+1,1:4)
              seg(i2+1,1  ) = seg(i3  ,2  )  ! Swap the order of the points
              seg(i2+1,2  ) = seg(i3  ,1  )
              seg(i2+1,3:4) = seg(i3  ,3:4)
              seg(i3  ,1:4) = seg(0   ,1:4)
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF (i3.EQ.nseg+1) 
     .    CALL ER('ProcessVoid','Zone perimeter gap detected',*99)
      ENDDO
    
      IF (debug) THEN
        DO i1 = 1, nseg
          WRITE(fp,*) 'PERIMETER:',i1,pts(seg(i1,1),1:2)
          WRITE(fp,*) '         :',i1,pts(seg(i1,2),1:2)
        ENDDO
      ENDIF


c     ------------------------------------------------------------------
c...  The perimeter of the void / zone is now defined, so look to see if 
c     any additional line segments or holes should be included:     
      tmp_nseg = nseg
      DO ivoid = 1, opt_eir%nvoid
        IF (opt_eir%void_zone(ivoid).NE.izone) CYCLE

c...    Set requested triangle area:
        IF (area.EQ.0.0.AND.opt_eir%void_res(ivoid).NE.0.0) 
     .    area = 0.5 * opt_eir%void_res(ivoid)**2

c...    Holes:
        IF     (opt_eir%void_hole(1,ivoid).EQ.-1.0.AND.
     .          opt_eir%void_hole(2,ivoid).EQ.-1.0) THEN
          DO isrf = 1, nsurface
            IF (surface(isrf)%type.NE.HOLE_IN_GRID) CYCLE
            x1 = surface(isrf)%v(1,1)  
            y1 = surface(isrf)%v(2,1)  
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'HOLE TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
              nhole = nhole + 1
              xhole(nhole) = x1
              yhole(nhole) = y1
              IF (debug) WRITE(fp,*) '   ADDING HOLE:',nhole,x1,y1
            ENDIF
          ENDDO

        ELSEIF (opt_eir%void_hole(1,ivoid).NE.0.0.AND.
     .          opt_eir%void_hole(2,ivoid).NE.0.0) THEN
          nhole = nhole + 1
          xhole(nhole) = opt_eir%void_hole(1,ivoid)
          yhole(nhole) = opt_eir%void_hole(2,ivoid)
        ENDIF

c...    Plasma conditions for all triangles in zone -- but note that
c       this does not imply any associated recycling in EIRENE:
        IF (opt_eir%void_ne(ivoid).GT.0.0) THEN
          ne = opt_eir%void_ne(ivoid)
          te = opt_eir%void_te(ivoid)
          ti = opt_eir%void_ti(ivoid)
        ENDIF

c...    Look for wall surfaces that are inside each zone:
        IF (opt_eir%void_add(1,ivoid).EQ.-1.AND.
     .      opt_eir%void_add(2,ivoid).EQ.-1) THEN
          DO isrf = 1, nsurface
            IF (surface(isrf)%type    .NE.VESSEL_WALL.OR.
     .          surface(isrf)%index(2).EQ.0          .OR.
     .          res(isrf).NE.0.0) CYCLE

            x1 = surface(isrf)%v(1,1)  ! Only check the first point of each surface,
            y1 = surface(isrf)%v(2,1)  ! assuming the wall isn't ill posed...
            IF (debug) WRITE(fp,'(4X,A,3I4)')
     .        'ADD  TRYING:',isrf,surface(isrf)%index(1:2)
            IF (PointInVoid(x1,y1,tmp_nseg,seg,MAXNSEG,
     .                      npts,pts,MAXNPTS,debug,fp)) THEN
c...          The point is inside the void perimeter, so add the segment to the list:
              res(isrf) = opt_eir%void_res(ivoid)
c...          Determine if the segment needs to be sub-divided:
              x1 = surface(isrf)%v(1,1)
              y1 = surface(isrf)%v(2,1)
              x2 = surface(isrf)%v(1,2)
              y2 = surface(isrf)%v(2,2)
              len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)
              IF (len.GT.DBLE(res(isrf))) THEN
                tstep = 1.0D0 / DBLE(INT(len/DBLE(res(isrf))) + 1)
              ELSE
                tstep = 1.0D0
              ENDIF
              IF (debug) THEN
                WRITE(fp,*) ' NEW ADD SEG:',isrf
c                WRITE(fp,*) '    I3    =',i3
c                WRITE(fp,*) '    PTS1,2=',pts(i3,1:2)
c                WRITE(fp,*) '    ISEG11=',pts(seg(iseg1,1),1:2)
c                WRITE(fp,*) '    ISEG12=',pts(seg(iseg1,2),1:2)
c                WRITE(fp,*) '    ISRF,RES:',isrf,res(isrf)
c                WRITE(fp,*) '    X1,Y1:',x1,y1
c                WRITE(fp,*) '    X2,Y2:',x2,y2
c                STOP 'dfsdfsd'
              ENDIF
              DO t = 0.0D0, 0.9999999D0, tstep 
                nseg = nseg + 1
                seg(nseg,1) = npts + 1
                seg(nseg,2) = npts + 2
                seg(nseg,4) = isrf
                npts = npts + 1
                pts(npts,1) = x1 + t * (x2 - x1)             ! Orientation is correct
                pts(npts,2) = y1 + t * (y2 - y1) 
                npts = npts + 1
                pts(npts,1) = x1 + (t + tstep) * (x2 - x1) 
                pts(npts,2) = y1 + (t + tstep) * (y2 - y1)
              ENDDO
            ENDIF
          ENDDO
        ENDIF    
 
      ENDDO


c...  Eliminate duplicate verticies:
      DO i2 = 1, npts
        DO i3 = i2+1, npts
          IF (pts(i2,1).NE.-999.0.AND.
     .        DABS(pts(i2,1)-pts(i3,1)).LT.DTOL.AND.
     .        DABS(pts(i2,2)-pts(i3,2)).LT.DTOL) THEN
            pts(i3,1) = -999.0
            pts(i3,2) = -999.0
            DO i4 = 1, nseg
              IF (seg(i4,1).EQ.i3) seg(i4,1) = i2
              IF (seg(i4,2).EQ.i3) seg(i4,2) = i2
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c...  Delete points:
      DO i2 = npts, 1, -1
        IF (pts(i2,1).EQ.-999.0) THEN
          DO i3 = i2, npts-1
            pts(i3,1) = pts(i3+1,1)
            pts(i3,2) = pts(i3+1,2)
          ENDDO
          DO i4 = 1, nseg
            IF (seg(i4,1).GE.i2) seg(i4,1) = seg(i4,1) - 1
            IF (seg(i4,2).GE.i2) seg(i4,2) = seg(i4,2) - 1
          ENDDO
          npts = npts - 1
        ENDIF
      ENDDO


      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=99)      
      WRITE(fp,*) npts,2,1,0
      DO i2 = 1, npts
c        WRITE(fp,'(I6,2F12.7)') i2,pts(i2,1),pts(i2,2)
        WRITE(fp,'(I6,2F19.14)') i2,pts(i2,1),pts(i2,2)
      ENDDO
      WRITE(fp,*) nseg,0
      DO i2 = 1, nseg
        IF (i2.EQ.nseg) THEN
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),0
        ELSE
          WRITE(fp,'(4I6)') i2,seg(i2,1),seg(i2,2),1
        ENDIF
      ENDDO
      WRITE(fp,*) nhole
      DO i2 = 1, nhole
        WRITE(fp,'(I6,2F12.7)') i2,xhole(i2),yhole(i2)
      ENDDO
      CLOSE (fp)

c...  Call triangle:
      IF (area.EQ.0.0) area = 0.01
      WRITE(command,10) 'triangle -p -q -a',area,' -Y triangle.poly'
 10   FORMAT(A,F10.8,A)
      WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      WRITE(0    ,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
      CALL CIssue(command(1:LEN_TRIM(command)),code)
      WRITE(eirfp,*) 'RETURN_CODE:',code

      icnt = icnt + 1

      CALL ReadPolyFile_06(izone,ne,te,ti)


      RETURN
 99   WRITE(0,*) 'I3    =',i3
      WRITE(0,*) 'PTS1,2=',pts(i3,1:2)
      STOP
      END
c
c ====================================================================
c
      SUBROUTINE RoutineInDevelopment
      IMPLICIT none
      RETURN
 99   STOP
      END
c
c ====================================================================
c
      SUBROUTINE LocalGridRefinement
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj




      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE CalculateTeProfile(inode1,inode2,s,target)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: inode1,inode2,target
      REAL*8 , INTENT(IN) :: s(0:icmax+1)

      INTEGER ion,ic,ic1,ic2,i1,i2,count
      INTEGER icstep
      LOGICAL cont
      REAL*8  te1,te2,Psol,L,k,stotal,vsign,adjust

      ion = 1

      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
c      te1 = DBLE(node(inode1)%te)
      IF (opt%bc(target).EQ.3) THEN
        IF (target.EQ.LO) te1 = te(ictarg(LO))
        IF (target.EQ.HI) te1 = te(ictarg(HI))
      ELSE
        te1 = DBLE(node(inode1)%te)
      ENDIF
      te2 = DBLE(node(inode2)%te)
c      IF (opt%bc(target).EQ.3) THEN
c        IF (target.EQ.LO) te2 = te(ictarg(LO))
c        IF (target.EQ.HI) te2 = te(ictarg(HI))
c      ELSE
c        te2 = DBLE(node(inode2)%te)
c      ENDIF

      k = DBLE(opt%te_kappa(target))
      L = 0.5D0 * DBLE(tube%smax)

c...  Set total power into flux tube:
      SELECTCASE (opt%te_ano_psol(target)) 
        CASE (1000)           ! ie ENEANO array from ref_plasma...? 
          Psol = 0.0D0  ! ???
c          STOP 'OKAY, READY'
        CASE (0)
c...      Conduction model estimate, all in upstream (initial guess only):
          Psol = 2.0D0/7.0D0 * (te2**3.5D0 - te1**3.5D0) * k / L   ! Stangeby, p 190
c          WRITE(0,*) 'PSOL:',psol,ic1,ic2
        CASE (1)
c...      From reference solution:
          Psol = ref_tube%eneano(target)
        CASE DEFAULT
          STOP 'USER OPTION NOT READY, YOU NINNY'
      ENDSELECT


      adjust = 0.0D0
      count = 0

              WRITE(logfp,'(A,I6)') 
     .          'TARGET:',target

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        count = count + 1

        CALL EvolveTeProfile(inode1,inode2,s,target,Psol,te1,te2)


c...    Analyse electron profiles and modify parameters, if necessary:

        IF (DABS(te(ic1)-te1).GT.MIN(0.001D0,0.05D0*te1)) THEN
c          vsign = SIGN(1.0D0,te(ic1)-te1) 
c          IF (target.EQ.LO) vsign = vsign * SIGN(1.0D0,te2-te1)
          vsign = SIGN(1.0D0,te(ic1)-te1) * SIGN(1.0D0,te2-te1)
          IF (adjust.EQ.0.0D0) THEN
            adjust = 0.3D0 * vsign
          ELSEIF ((adjust.GT.0.0D0.AND.vsign.EQ.-1.0D0).OR.
     .            (adjust.LT.0.0D0.AND.vsign.EQ. 1.0D0)) THEN
            adjust = -0.3D0 * adjust
          ENDIF
          cont = .TRUE.
        ENDIF

c...    Analyse electron profiles and modify parameters, if necessary:

        SELECTCASE (opt%bc(target))  ! (opt%te_ano_psol(target))  ! *** CHANGE VARIABLE *** 
          CASE (1)
c...        Adjust the anomalous power term on the half ring:
c...        Scale Psol to improve match to specified target temperature:
            SELECTCASE (1)
              CASE (1)  ! Power into SOL:
                Psol = Psol * (1.0D0 + adjust)
              CASE DEFAULT
                STOP 'NOT READY, SORRY'
            ENDSELECT
            IF (log.GE.2) THEN
              WRITE(logfp,'(A,I3,3F11.4,1P,2E14.6,0P)') 
     .          'ADJUST:',count,te(ic1),te1,vsign,adjust,Psol
            ENDIF

          CASE (3)

            node(inode2)%te = node(inode2)%te - adjust
            te2 = DBLE(node(inode2)%te)

            IF (log.GE.2) THEN
              WRITE(logfp,'(A,3F11.4,1P,E14.6,0P,F11.4)') 
     .          'ADJUST:',te(ic1),te1,vsign,adjust,te2
            ENDIF

          CASE DEFAULT
            CALL ER('EvolveTeProfile','Bad BC',*99)
        ENDSELECT

c...    Analyse ion profiles:

c...    If convection being used, need to run for extra iterations
c       to make sure the solution is converged:

        IF (count.EQ.50) THEN
          te(ic1:ic2) = te1
          WRITE(0,*) '*** ENERGY MODEL FAILED ***'
          WRITE(logfp,*) '*** ENERGY MODEL FAILED ***'
          EXIT
c          STOP 'EXCESS ITERATIONS'
        ENDIF

      ENDDO





c...  Store Psol information:
      tube%eneano(target) = Psol   ! Keep this up...?

c...  Linear distribution of anomalous power for now...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(1)
        CASE(2)
c...      Everything in at the midplane:
          eneano(ic2+icstep) = Psol / sdelta(ic2+icstep)  ! Approximate..!
c          WRITE(0,*) 'ENEA:',psol,ic2+icstep
        CASE(3)
c...      Distributed uniformly along the half ring:
          IF (ic1.GT.ic2) THEN
            ic  = ic1
            ic1 = ic2
            ic2 = ic
          ENDIF
          ic1 = MAX(ic1,1)
          ic2 = MIN(ic2,icmax)
          IF (target.EQ.LO) ic2 = ic2 - 1
          IF (target.EQ.HI) ic1 = ic1 + 1
          stotal = SUM(sdelta(ic1:ic2))
          eneano(ic1:ic2) = Psol / stotal
c          WRITE(0,*) 'ENEA:',psol,ic1,ic2,icstep,stotal
        CASE DEFAULT
          CALL ER('EvolveTeProfile','Bad TE_ANO',*99)
      ENDSELECT


      IF (log.GE.2) THEN
        WRITE(logfp,'(A)') 'DONE TE ITERATIONS'
      ENDIF

      IF (target.EQ.HI.AND.opt%bc(target).EQ.1) THEN
        WRITE(logfp,*) 'INTEGRATING ENERGY SOURCES AGAIN:'
        CALL IntegrateSources(3)

        ic1 = ictarg(LO)
        ic2 = ictarg(HI)

c        qe(ic1) = -eneint(icmid  ,1)
c        qe(ic2) = -(eneint(TOTAL,1)-eneint(icmid,1))

        WRITE(logfp,*) 'QE   :',qe(0),qe(icmax+1)
        WRITE(logfp,*) 'ISAT :',isat(ictarg(LO:HI),ion)
        WRITE(logfp,*) 'TE   :',te(ictarg(LO:HI))

        gamma(LO) = qe(0)       / (isat(ic1,ion)*ECH) / te(ic1)
        gamma(HI) = qe(icmax+1) / (isat(ic2,ion)*ECH) / te(ic2)

c        gamma = 3.0D0
c        qe(ic1) = 3.0D0 * (ECH*te(ic1)) * isat(ic1,ion)
c        qe(ic2) = 3.0D0 * (ECH*te(ic2)) * isat(ic1,ion)

        WRITE(logfp,*) 'GAMMA:',gamma(LO:HI)
        WRITE(logfp,*) '...and again...'

        CALL IntegrateSources(3)
c        gamma(LO:HI) = 3.0D0
c        STOP 'sdgsdg'
      ENDIF


      RETURN
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE EvolveTeProfile(inode1,inode2,s,target,Psol,te1,te2)
      USE mod_sol28_params
      USE mod_sol28_solver
      IMPLICIT none
 
      INTEGER, INTENT(IN) :: inode1,inode2,target
      REAL*8 , INTENT(IN) :: s(0:icmax+1),Psol,te1,te2

      INTEGER :: ic,ic1,ic2,icstep,count,ion
      REAL*8  :: tgrad,x,k,frac,tarsign,
     .           flux,qano(0:S28_MAXNKS+1),
     .           qe_src(0:S28_MAXNKS+1),
     .           qconv_e(0:S28_MAXNKS+1),qconv_i(0:S28_MAXNKS+1)
      LOGICAL :: cont

c      te1 = DBLE(node(inode1)%te)
c      te2 = DBLE(node(inode2)%te)
      ic1 = node(inode1)%icell
      ic2 = node(inode2)%icell
      icstep = SIGN(1,ic1-ic2)
      tarsign = DBLE(icstep)

c      IF (opt%bc(target).GE.2.AND.
c     .    (ic1.EQ.ictarg(target).AND.te(ic1).EQ.0.0D0)) THEN
c        WRITE(0,*) 'CYCLING -------------------------'
c        RETURN
c      ENDIF


c      WRITE(logfp,*) 'TE: ',ic1,ic2,icstep,target
c      WRITE(logfp,*) '    ',te1,te2

      ion = 1


c...  Integrate e/i cf/vol sources/sinks:
      qe_src = 0.0D0
      IF (icstep.EQ.-1) THEN
        qe_src(ic1:ic2) = eneint(ic2,1) - eneint(ic1:ic2,1)  ! Should I set the INT arrays from 0:icmax+1, and have 
        IF (ic1.EQ.0) qe_src(ic1) = eneint(ic2,1)            ! TOTAL = icmax+1 to make things like this simpler?
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ELSE
        qe_src(ic2:ic1) = eneint(ic2:ic1,1) - eneint(ic2,1)  
        IF (ic1.EQ.icmax+1) qe_src(ic1) = eneint(TOTAL,1)-eneint(ic2,1)            
c        WRITE(0,*) 'QE:',qe(0:icmax+1)
c        WRITE(0,*) 'QE:',eneint(ic2,1),ic1,icmax+1
c        WRITE(0,*) 'QE:',enesrc(0:icmax+1,1)
      ENDIF

c...  Initializations:
      k = DBLE(opt%te_kappa(target))
      x = 5.0D0 / 2.0D0
      qano  = 0.0D0
      te(ic2) = te2
 
c...  Setup anomalous energy input:    NEED TO ADD EFFICIENCY, ie DON'T ALWAYS HAVE TO CALCULATE EVERYTHING...
      SELECTCASE (opt%te_ano(target))
        CASE(1000)
        CASE(0)
c...      None:
        CASE(1)
c...      From reference solution - ANO source incorporated into QE_SRC via ENEINT:
        CASE(2)
c...      All power in at the symmetry point:
          qano = Psol 
        CASE(3)
c...      Power distributed uniformly about the symmetry point(!):
          DO ic = ic2+icstep, ic1, icstep  ! No power assigned to the symmetry point cell...
            frac = (s(ic ) - (s(ic2) - 0.5D0 * sdelta(ic2))) / 
     .             (s(ic1) - (s(ic2) - 0.5D0 * sdelta(ic2)))
            qano(ic) = Psol * frac                      
          ENDDO
c        CASE(3)
c...      Between x-points:
c        CASE(4)
c...      B-field depenence... use CalcProfile?   What if heat is not centred at symmetry point? 
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

c...  Electron convected heat flux:
      SELECTCASE (opt%te_conv(target))
        CASE(0)
c...      None:
          DO ic = ic2, ic1, icstep
            qconv_e(ic) = 0.0D0
          ENDDO
        CASE(1)
          te(ic1) = te1
          DO ic = ic2, ic1, icstep
c            WRITE(logfp,*) 'CONV:',ic,te(ic),ne(ic),vi(ic,ion)
            IF (te(ic).NE.0.0D0) THEN
              flux = ne(ic) * vi(ic,ion) * tarsign
              qconv_e(ic) = 5.0D0 / 2.0D0 * ECH * te(ic) * flux
            ENDIF
          ENDDO
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

c...  Ion convected heat flux:
      SELECTCASE (0)  ! (opt%ti_conv(target))  ! Don't need this when solving the electron channel...
        CASE(0)
c...      None:
          DO ic = ic2, ic1, icstep
            qconv_i(ic) = 0.0D0
          ENDDO
        CASE(1)
          te(ic1) = te1                        ! *** ASSUMES Ti = Te, on this line and 4 below ***
          DO ic = ic2, ic1, icstep
            IF (te(ic).NE.0.0D0) THEN
              flux = ne(ic) * vi(ic,ion) * tarsign
              qconv_i(ic) = (5.0D0 / 2.0D0 * ECH * te(ic) +  
     .                       0.5D0 * mi(ion) * vi(ic,ion)**2) * flux
            ENDIF
          ENDDO
        CASE DEFAULT
          STOP 'USER TE OPTION NOT READY'
      ENDSELECT

      qconv = qconv_e + qconv_i

c...  Electron-ion equilibration:
      IF (.FALSE.) THEN
      ENDIF

c...  Calculate Te profile:
      te(ic1) = 0.0D0
      DO ic = ic2+icstep, ic1, icstep
        qcond(ic) = qano(ic) + qe_src(ic) - qconv(ic) 
        tgrad = -qcond(ic) / (k * te(ic-icstep)**x)                ! Reference...
c...    Note: This is not strictly correct since the convection is a cell centered quantity
c       but the temperature evolution occurs between cell centres, i.e. the convecion
c       term should be a wieghted average, or some such, of the convection
c       across the 2 cells -- similar inconcistencies can likely be found elsewhere
c       and will likely need to be resolved before a full cross-field transport model
c       can work... might even run into trouble before then.  For a similar slop, the
c       treatment of the last half-cell on each ring also needs careful review.
        te(ic) = te(ic-icstep) + (s(ic-icstep) - s(ic)) * tgrad

        qcond(ic) = qcond(ic) !* tarsign
        qconv(ic) = qconv(ic) !* tarsign

c        qe(ic) = qcond(ic) + qconv(ic)


        IF (log.GE.2) THEN
          IF (ic.EQ.ic2+icstep) THEN
            WRITE(logfp,'(A4,A8,2X,3A12,2X,2A12)')
     .        'IND','Te','qcond','qconv_e','qconv_i','qano','qe_src'   
          ENDIF
          WRITE(logfp,'(I4,F8.2,2X,1P,3D12.4,2X,2D12.4,1P,
     .                  2E10.2)')
     .      ic,
     .      te(ic),
     .      qcond(ic),qconv_e(ic),qconv_i(ic),
     .      qano(ic),qe_src(ic),
     .      tgrad,(s(ic-icstep) - s(ic))
        ENDIF
c          IF (target.EQ.LO)
c     .      WRITE(0,'(A,2I6,2F10.2,1P,6E10.2,2X,2E10.2,0P)') 
c     .        '  Te->',ic,ic-icstep,s(ic),te(ic),vi(ic,ion),flux,
c     .        qano(ic),qe(ic),qconv(ic),qcond(ic),vi(ic,ion),ne(ic)
c        IF (te(ic).LE.0.1D0*te1) EXIT  
c        IF (te(ic).LE.0.5D0*te1.OR.te(ic).GT.1.5D0*te2) THEN
c          te(ic) = MIN(MAX(te(ic),te1),te2)
c          EXIT  ! There is no way to avoid the near target dip in 
c        ENDIF
        IF (te(ic).LE.0.5D0*MIN(te1,te2)) EXIT  ! There is no way to avoid the near target dip in 
c        IF (te(ic).LE.0.5D0*te1) EXIT  ! There is no way to avoid the near target dip in 
c        IF (te(ic).LT.te1) EXIT        ! Te if Qe is ill-posed?
      ENDDO


c      WRITE(0,*) 'QCOND:',qcond(ic1),qcond(ic2)
c      WRITE(0,*) 'TE   :',te(ic1),te(ic2)
c      WRITE(0,*) 'PSOL :',psol,ic1,ic2,opt%te_ano_psol(target)



      RETURN
 99   STOP
      END





