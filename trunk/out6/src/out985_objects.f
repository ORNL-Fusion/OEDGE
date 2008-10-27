c     -*-Fortran-*-
c
c ======================================================================
c
      SUBROUTINE LoadImageReconstruction(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER, INTENT(IN) :: ielement

      INTEGER AddVertex,AddSurface

      INTEGER   fp,count,istart,idum1,i1,i2,ivol,fobj,
     .          iobj,iobj2,iside,iside2,isrf,isrf2
      LOGICAL   map_found
      REAL*8    val,newvtx(3,4),v1(3),v2(3),v3(3),v4(3)
      CHARACTER file*1024,buffer*2048
      TYPE(type_surface) newsrf

      REAL*8, PARAMETER :: DTOL = 1.0D-10

      fp = 99


      WRITE(0,*) 'LOADING IMAGE RECONSTRUCTIONS'//
     .  opt%obj_fname(ielement)


c...  Load up surface and vertex arrays:

      istart = nsrf + 1
      fobj   = nobj + 1
      ivol   = 0

      file = opt%obj_fname(ielement)

      OPEN(fp,FILE=TRIM(file),FORM='FORMATTED',STATUS='OLD',ERR=98)     
 
      count = 0

      READ(fp,*)  ! Clear index line

      DO WHILE (.TRUE.)       

        READ(fp,'(A2048)',END=10) buffer         

        newvtx = 0.0D0
        READ(buffer,*,ERR=10) 
     .    idum1,val,idum1,(newvtx(1:2,i1),i1=1,idum1)

c        newvtx(1:2,1) = newvtx(1:2,1) + 0.0000005
c        newvtx(1:2,2) = newvtx(1:2,2) + 0.0000005
c        newvtx(1:2,3) = newvtx(1:2,3) + 0.0000005
c        newvtx(1:2,4) = newvtx(1:2,4) + 0.0000005

        IF (nobj+1.GT.MAX3D) 
     .    CALL ER('LoadImageRecon','Insufficient array bounds '//
     .            'for all objects',*99)    

        nobj = nobj + 1
        ivol = ivol + 1

c...    Extract vertices:
c        WRITE(0,*) 'IMAGE RECON. NOBJ:',nobj
c        WRITE(0,*) '   :',idum1,val

        obj(nobj)%index       = ielement  ! nobj
        obj(nobj)%type        = OP_INTEGRATION_VOLUME
        obj(nobj)%subtype     = OP_INVERSION_GRID  ! OK for now...
        obj(nobj)%mode        = 0      
        obj(nobj)%surface     = 1      ! SOLID ???
        obj(nobj)%wedge1      = 0
        obj(nobj)%wedge2      = 0
        obj(nobj)%colour      = 3
        obj(nobj)%orientation = 1      ! CW
        obj(nobj)%ik          = 0
        obj(nobj)%ir          = 0
        obj(nobj)%in          = 0
        obj(nobj)%ivolume     = ivol
        obj(nobj)%nside       = idum1
        obj(nobj)%reflec      = 0
        obj(nobj)%quantity(1) = 1.0D0 ! val ! / 1.0E+18 ! 1.0 ! val


c        WRITE(0,*) 'quantit:',nobj,obj(nobj)%quantity(1)
c..     Defunct:
        obj(nobj)%nsur        = 0
        obj(nobj)%ipts(2,1)   = 0
        obj(nobj)%nmap(1)     = 0

        DO i1 = 1, obj(nobj)%nside
          i2 = i1 + 1
          IF (i2.GT.obj(nobj)%nside) i2 = 1  
          newsrf%type = -1
          newsrf%nvtx =  2
          newsrf%ivtx(1) = AddVertex(newvtx(1,i1))
          newsrf%ivtx(2) = AddVertex(newvtx(1,i2))
          obj(nobj)%iside(i1,1) = AddSurface(newsrf)
          obj(nobj)%iside(i1,2) = obj(nobj)%iside(i1,1)
          obj(nobj)%gsur(i1)    = GT_TC
          obj(nobj)%tsur(i1)    = -1
c          WRITE(0,*) '   -',i1,obj(nobj)%iside(i1,1)
        ENDDO

      ENDDO
 10   CONTINUE
      CLOSE(fp)


      IF (.TRUE.) THEN
c...    Force rebuild of connection map:

        DO iobj = fobj, nobj
          IF (MOD(iobj-fobj+1,1000).EQ.0) 
     .      WRITE(0,*) '  BUILDING CONNECTION MAP',iobj,nobj-fobj+1
c          IF (iobj.EQ.5) STOP ' debugging...'
          DO iside = 1, obj(iobj)%nside
            isrf = obj(iobj)%iside(iside,1)     
            v1(1:3) = vtx(1:3,srf(isrf)%ivtx(1))
            v2(1:3) = vtx(1:3,srf(isrf)%ivtx(2))

            obj(iobj)%tsur(iside) = SP_GRID_BOUNDARY
            obj(iobj)%imap(1,iside) = iobj
            obj(iobj)%isur(1,iside) = iside

c            WRITE(0,*) '  IOBJ :',iobj,iside,isrf
c            WRITE(0,*) '       :',v1(1:2)
c            WRITE(0,*) '       :',v2(1:2)

            map_found = .FALSE. 
            DO iobj2 = MAX(fobj,iobj-1050),MIN(nobj,iobj+1050) ! Crude...
              IF (iobj.EQ.iobj2) CYCLE
              DO iside2 = 1, obj(iobj2)%nside
                isrf2 = obj(iobj2)%iside(iside2,1)     
                v3(1:3) = vtx(1:3,srf(isrf2)%ivtx(1))
                v4(1:3) = vtx(1:3,srf(isrf2)%ivtx(2))
c                WRITE(0,*) '  IOBJ2:',iobj2,iside2,isrf2
c                WRITE(0,*) '       :',v1(1:2)
c                WRITE(0,*) '       :',v2(1:2)
c                WRITE(0,*) '       :',v3(1:2)
c                WRITE(0,*) '       :',v4(1:2)
                IF ((DABS(v1(1)-v3(1)).LT.DTOL.AND.
     .               DABS(v1(2)-v3(2)).LT.DTOL.AND.
     .               DABS(v2(1)-v4(1)).LT.DTOL.AND.
     .               DABS(v2(2)-v4(2)).LT.DTOL).OR.
     .              (DABS(v1(1)-v4(1)).LT.DTOL.AND.
     .               DABS(v1(2)-v4(2)).LT.DTOL.AND.
     .               DABS(v2(1)-v3(1)).LT.DTOL.AND.
     .               DABS(v2(2)-v3(2)).LT.DTOL)) THEN
                  obj(iobj)%tsur(iside) = SP_GRID_SURFACE
                  obj(iobj)%imap(1,iside) = iobj2
                  obj(iobj)%isur(1,iside) = iside2
                  map_found = .TRUE.
c                  WRITE(0,*) '  MAP:',iobj,iside,iobj2,iside2
                  EXIT
                ENDIF
              ENDDO
              IF (map_found) EXIT
            ENDDO 
c            WRITE(0,*) '  :',iobj,iside,
c     .        obj(iobj)%imap(1,iside),obj(iobj)%isur(1,iside)
c            IF (iobj.EQ.5) STOP 'sdfsdf'
          ENDDO
        ENDDO
c          WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c          WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c          WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)
      ENDIF

      WRITE(0,*) 'DONE'

      WRITE(6,*) 'DESPERATE'
      DO iobj = 1, nobj
        DO iside = 1, 4
          isrf = obj(iobj)%iside(iside,1)     
          WRITE(6,*) iobj,iside,vtx(1:2,srf(isrf)%ivtx(1))
          WRITE(6,*) iobj,iside,vtx(1:2,srf(isrf)%ivtx(2))
c          obj(iobj)%v(1,iside) = vtx(1,srf(isrf)%ivtx(1))
c          obj(iobj)%v(2,iside) = vtx(2,srf(isrf)%ivtx(1))
c          obj(iobj)%v(3,iside) = 0.0D0
c          obj(nobj)%npts(iside) = 2
c          obj(nobj)%ipts(1,iside) = 1
c          obj(nobj)%ipts(2,iside) = 2
c          obj(nobj)%nmap(iside)   = 1
        ENDDO
c        obj(nobj)%nside = 0
c        obj(nobj)%nsur  = 4
      ENDDO
c      nvtx = 0
c      nsrf = 0


      RETURN
 98   WRITE(0,*) 'ERROR LoadImageReconstruction: File not found'
      WRITE(0,*) '   '//TRIM(file)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE LoadLineSegmentFile(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

c...  Input:
      INTEGER ielement

      INTEGER AddVertex,AddSurface

      INTEGER   fp,count,istart,n1,n2,idum1
      REAL*8    newvtx(3,2)
      CHARACTER file*1024,buffer*2048
      TYPE(type_surface) newsrf


      fp = 99


      WRITE(0,*) 'LOADING LINE SEGMENTS'//
     .  opt%obj_fname(ielement)


c...  Load up surface and vertex arrays:

      istart = nsrf + 1

      file = opt%obj_fname(ielement)

      OPEN(fp,FILE=file(1:LEN_TRIM(file)),
     .     FORM='FORMATTED',STATUS='OLD',ERR=98)     
 
      n1 = opt%obj_n(ielement,1)
      n2 = opt%obj_n(ielement,2)

      count = 0

      DO WHILE (.TRUE.)       

        READ(fp,'(A2048)',END=10) buffer         

        IF (buffer(1:1).EQ.'*'.OR.LEN_TRIM(buffer).EQ.0) THEN
          count = 0
          CYCLE
        ENDIF

c...    Extract vertices:
        newvtx(1:3,1) = 0.0D0

        SELECTCASE (1)
          CASE(1)
            READ(buffer,*,ERR=10) newvtx(1,1),newvtx(2,1)
          CASE(2)
            READ(buffer,*,ERR=10) newvtx(2,1),newvtx(1,1)
          CASE DEFAULT
            CALL ER('LoadLineSegmentFile','Unknown orientation '//
     .              'option',*99)
        ENDSELECT


        IF (count.GT.0.AND.
     .      ((count.GE.n1).OR.(-1.EQ.n1)).AND.
     .      ((count.LE.n2).OR.(-1.EQ.n2))) THEN

          newsrf%type = SP_LINE_SEGMENT   ! LEFT OFF -- NEED TO ADD THIS IN?
          newsrf%nvtx = 2
          newsrf%ivtx(1) = AddVertex(newvtx(1,1))
          newsrf%ivtx(2) = AddVertex(newvtx(1,2))

          idum1 = AddSurface(newsrf)
        ENDIF

        newvtx(1:3,2) = newvtx(1:3,1)

        count = count + 1

      ENDDO
 10   CONTINUE
      CLOSE(fp)

      WRITE(0,*) 'VERTICES? SURFACES?',nvtx,nsrf
      WRITE(0,*) '                   ',istart

c...  Assign object(s):    

      WRITE(0,*) 'NOBJ:',nobj,MAX3D

      IF (nobj+1.GT.MAX3D) 
     .  CALL ER('LoadVesselStructures','Insufficient array bounds '//
     .          'for all objects',*99)    

      IF (istart.GT.nsrf) THEN
        WRITE(0,*) 'LoadVesselStructures: Strange, no objects loaded'
        RETURN
      ENDIF

      nobj = nobj + 1
      WRITE(0,*) 'VESSEL STRUCTURE IOBJ:',nobj

      obj(nobj)%index       = ielement  ! nobj
      obj(nobj)%type        = OP_EMPTY
      obj(nobj)%mode        = 0      
      obj(nobj)%surface     = 1      ! SOLID
      obj(nobj)%wedge1      = 0
      obj(nobj)%wedge2      = 0
      obj(nobj)%colour      = 1
      obj(nobj)%orientation = 1      ! CW
      obj(nobj)%ik          = 0
      obj(nobj)%ir          = 0
      obj(nobj)%in          = -1  ! What should this be?
      obj(nobj)%ivolume     = 0
      obj(nobj)%nside       = 1
      obj(nobj)%iside(1,1)  = istart ! Start index of range of surfaces in surface array, from loading code above
      obj(nobj)%iside(1,2)  = nsrf   ! End index of range of surfaces in surface array
      obj(nobj)%gsur(1)     = GT_TC
      obj(nobj)%tsur(1)     = SP_VESSEL_WALL
      obj(nobj)%reflec(1)   = opt%obj_reflec(ielement)
c..   Defunct:
      obj(nobj)%nsur        = 0
      obj(nobj)%ipts(2,1)   = 0
      obj(nobj)%nmap(1)     = 0

      WRITE(0,*) 'DONE'


      RETURN
 98   WRITE(0,*) 'ERROR LoadLineSegmentFile: File not found'
      WRITE(0,*) '   '//file(1:LEN_TRIM(file))
 99   STOP
      END
c
c ====================================================================== 
c
c
c             side 3
c         3------------4            1  2  3  4  5 ....   nxbin
c         |            |            nxbin+1 ......     2*nxbin
c         |            |            .... 
c         |            |            ....           nybin*nxbin
c  side 2 |            | side 4
c         |            |
c         |            |
c         |            |
c         2------------1
c             side 1
c
c
c 
c
c
      SUBROUTINE BuildInversionMesh(ielement)
      USE mod_out985
      USE mod_out985_variables
      IMPLICIT none

      INTEGER ielement



      LOGICAL CheckInversionCell

      INTEGER i1,i2,i3,ix,iy,ncells,isector,nsector,nxbin,nybin,ivol,
     .        istart,iend,fobj,iobj,isid,iobj2,ik,ir,iside
      LOGICAL outofgrid
      REAL    angle,dangle,ang
      REAL*8  p1(3,8),p2(3,8),xcen,ycen,xwidth,ywidth,xdelta,ydelta,
     .        x1,z1,xorigin,yorigin,distsep,distwal,
     .        frac,maxdist,maxdiststd,maxdistxpt,
     .        x(4),y(4),x2(4),y2(4)

      REAL, PARAMETER :: RADDEG = 57.29577952

      REAL*8, PARAMETER :: DTOL = 1.0D-10

c...    Number of objects per toroidal segment:

c      WRITE(0,*) 'BUILDING RADIAL MAP'


      nsector = opt%obj_nsector
      IF (nsector.EQ.-1) nsector = grd_ntorseg ! eirntorseg

      IF (nsector.EQ.0) THEN
        WRITE(0,*) 'OVER-RIDING nsector IN BUILDINVERSIONMESH'
        nsector = 48
      ENDIF

      dangle = 360.0 / REAL(nsector) / RADDEG

      fobj = nobj + 1

      IF (ielement.NE.0) THEN
        nxbin = opt%obj_n(ielement,1)
        nybin = opt%obj_n(ielement,2)
        ncells = nxbin * nybin

        xcen = 0.5D0*DBLE(opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
        ycen = 0.5D0*DBLE(opt%obj_z(ielement,1) + opt%obj_z(ielement,2))
c        xcen = 0.5D0 * (opt%obj_r(ielement,1) + opt%obj_r(ielement,2))
c        ycen = 0.5D0 * (opt%obj_z(ielement,1) + opt%obj_z(ielement,2))

        xwidth = DBLE(opt%obj_r(ielement,2) - opt%obj_r(ielement,1))
        ywidth = DBLE(opt%obj_z(ielement,2) - opt%obj_z(ielement,1))
c        xwidth = opt%obj_r(ielement,2) - opt%obj_r(ielement,1)
c        ywidth = opt%obj_z(ielement,2) - opt%obj_z(ielement,1)

        xdelta = xwidth / DBLE(nxbin)
        ydelta = ywidth / DBLE(nybin) 

        xorigin = xcen - 0.5D0 * DBLE(nxbin) * xdelta
        yorigin = ycen - 0.5D0 * DBLE(nybin) * ydelta

        WRITE(6,*) 'XCEN,YCEN',xcen,ycen
        WRITE(6,*) 'delta',xdelta,ydelta
        WRITE(6,*) 'wid',xwidth,ywidth
        WRITE(6,*) 'origin',xorigin,yorigin
         

        IF (opt%obj_option(ielement).EQ.6) nybin = 2 * nybin - 1
      ELSE
c...    Original code:
        STOP 'CODE OBSOLETE'
      ENDIF


      IF ((ielement.NE.0.AND.
     .     (opt%obj_option(ielement).EQ.2.OR.
     .      opt%obj_option(ielement).EQ.4.OR.
     .      opt%obj_option(ielement).EQ.6))) THEN
c...    (Perfect) toroidal symmetry:

        ivol = 0  ! Integration volume index 

        maxdiststd = 0.15D0
        maxdistxpt = 0.25D0

        DO iy = nybin, 1, -1
          DO ix = 1, nxbin

c...        Check if the inversion cell is inside the OSM fluid grid:
            IF (opt%obj_option(ielement).EQ.4) THEN
              xcen = xorigin + xdelta * (DBLE(ix) - 0.5D0)
              ycen = yorigin + ydelta * (DBLE(iy) - 0.5D0)           
              IF (.NOT.CheckInversionCell(1,xcen,ycen)) CYCLE
            ENDIF

            IF (nobj+1.GT.MAX3D) 
     .        CALL ER('BuildObjects','Insufficient array bounds '//
     .                'for all objects',*98)     

            ivol = ivol + 1

            nobj = nobj + 1
            obj(nobj)%index       = ielement  ! nobj
            obj(nobj)%type        = OP_INTEGRATION_VOLUME
            obj(nobj)%subtype     = OP_INVERSION_GRID
            obj(nobj)%mode        = 0      
            obj(nobj)%surface     = 1      ! SOLID
            obj(nobj)%wedge1      = 0
            obj(nobj)%wedge2      = 0
            obj(nobj)%colour      = 3
            obj(nobj)%orientation = 1      ! CW
            obj(nobj)%ik          = ix
            obj(nobj)%ir          = iy
            obj(nobj)%in          = 0
            obj(nobj)%ivolume     = ivol
            obj(nobj)%nsur        = 4
            obj(nobj)%gsur(1:4)   = GT_TC
            obj(nobj)%nver        = 4
            obj(nobj)%reflec(1:4) = 0
 
            obj(nobj)%quantity = 1.0
c            obj(nobj)%quantity(1) = 1.0

c...        Checkerboard:
            IF (.FALSE.) THEN
              IF (MOD(ix,2).EQ.0) THEN
                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%quantity(1) = 1.0
                ELSE
                  obj(nobj)%quantity(1) = 0.0
                ENDIF
              ELSE
                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%quantity(1) = 0.0
                ELSE
                  obj(nobj)%quantity(1) = 1.0
                ENDIF
              ENDIF
            ENDIF
c            IF (ix.GT.nxbin/2) obj(nobj)%quantity(1) = 0.0
c            obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)

            SELECTCASE (opt%obj_option(ielement)) 
              CASE(2,4)  ! Rectangles
c...            Vertices:
                obj(nobj)%v(1,1) = xorigin + xdelta * DBLE(ix)   
                obj(nobj)%v(2,1) = yorigin + ydelta * DBLE(iy-1)
                obj(nobj)%v(3,1) = 0.0D0
                obj(nobj)%v(1,2) = xorigin + xdelta * DBLE(ix-1)
                obj(nobj)%v(2,2) = yorigin + ydelta * DBLE(iy-1)
                obj(nobj)%v(3,2) = 0.0D0
                obj(nobj)%v(1,3) = xorigin + xdelta * DBLE(ix-1)
                obj(nobj)%v(2,3) = yorigin + ydelta * DBLE(iy)
                obj(nobj)%v(3,3) = 0.0D0
                obj(nobj)%v(1,4) = xorigin + xdelta * DBLE(ix)
                obj(nobj)%v(2,4) = yorigin + ydelta * DBLE(iy)
                obj(nobj)%v(3,4) = 0.0D0

c...            Surface 1:
                IF (iy.EQ.1) THEN    
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj
                  obj(nobj)%isur(1,1) = 1
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + nxbin
                  obj(nobj)%isur(1,1) = 3
                ENDIF
          
                obj(nobj)%npts(1) = 2
                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
c...            Surface 2:
                IF (ix.EQ.1) THEN
                  obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(2) = 1
                  obj(nobj)%imap(1,2) = nobj
                  obj(nobj)%isur(1,2) = 2
                ELSE
                  obj(nobj)%tsur(2) = SP_GRID_SURFACE
                  obj(nobj)%nmap(2) = 1
                  obj(nobj)%imap(1,2) = nobj - 1
                  obj(nobj)%isur(1,2) = 4
                ENDIF
                obj(nobj)%npts(2) = 2
                obj(nobj)%ipts(1,2) = 2
                obj(nobj)%ipts(2,2) = 3
c...            Surface 4:
                IF (iy.EQ.nybin) THEN   
                  obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(3) = 1
                  obj(nobj)%imap(1,3) = nobj
                  obj(nobj)%isur(1,3) = 3
                ELSE
                  obj(nobj)%tsur(3) = SP_GRID_SURFACE
                  obj(nobj)%nmap(3) = 1
                  obj(nobj)%imap(1,3) = nobj - nxbin
                  obj(nobj)%isur(1,3) = 1
                ENDIF
                obj(nobj)%npts(3) = 2
                obj(nobj)%ipts(1,3) = 3
                obj(nobj)%ipts(2,3) = 4
c...            Surface 4:
                IF (ix.EQ.nxbin) THEN
                  obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(4) = 1
                  obj(nobj)%imap(1,4) = nobj
                  obj(nobj)%isur(1,4) = 4
                ELSE
                  obj(nobj)%tsur(4) = SP_GRID_SURFACE
                  obj(nobj)%nmap(4) = 1
                  obj(nobj)%imap(1,4) = nobj + 1
                  obj(nobj)%isur(1,4) = 2
                ENDIF
                obj(nobj)%npts(4) = 2
                obj(nobj)%ipts(1,4) = 4
                obj(nobj)%ipts(2,4) = 1

              CASE(6) ! Diamonds
c...            Vertices:
                obj(nobj)%v(1,1) = xorigin + xdelta * (DBLE(ix  )-0.5D0)
                obj(nobj)%v(2,1) = yorigin + ydelta *  DBLE(iy-1)
                obj(nobj)%v(3,1) = 0.0D0
                obj(nobj)%v(1,2) = xorigin + xdelta *  DBLE(ix-1)
                obj(nobj)%v(2,2) = yorigin + ydelta * (DBLE(iy-1)+0.5D0)
                obj(nobj)%v(3,2) = 0.0D0
                obj(nobj)%v(1,3) = xorigin + xdelta * (DBLE(ix-1)+0.5D0)
                obj(nobj)%v(2,3) = yorigin + ydelta *  DBLE(iy  )
                obj(nobj)%v(3,3) = 0.0D0
                obj(nobj)%v(1,4) = xorigin + xdelta *  DBLE(ix  )
                obj(nobj)%v(2,4) = yorigin + ydelta * (DBLE(iy  )-0.5D0)
                obj(nobj)%v(3,4) = 0.0D0

                IF (MOD(iy,2).EQ.0) THEN 
                  obj(nobj)%v(1,1:4) = obj(nobj)%v(1,1:4) + 0.5D0*xdelta   
                ELSE
c                  obj(nobj)%v(2,1:4) = obj(nobj)%v(2,1:4) + 0.5D0*ydelta   
                ENDIF

                obj(nobj)%v(2,1:4) = obj(nobj)%v(2,1:4) - 
     .                               0.5D0*DBLE(iy-1)*ydelta   

                obj(nobj)%npts(1:4) = 2
                obj(nobj)%nmap(1:4) = 1

                obj(nobj)%ipts(1,1) = 1
                obj(nobj)%ipts(2,1) = 2
                obj(nobj)%ipts(1,2) = 2
                obj(nobj)%ipts(2,2) = 3
                obj(nobj)%ipts(1,3) = 3
                obj(nobj)%ipts(2,3) = 4
                obj(nobj)%ipts(1,4) = 4
                obj(nobj)%ipts(2,4) = 1

              CASE DEFAULT
                CALL ER('BuildInversionMesh','Unknown option',*99)
            ENDSELECT  
        
          ENDDO
        ENDDO
        

        IF (opt%obj_option(ielement).EQ.4.OR.
     .      opt%obj_option(ielement).EQ.6) THEN
c...      Force rebuild of connection map:

          DO iobj = fobj, nobj
            DO isid = 1, obj(iobj)%nsur

              x(1:4) = obj(iobj)%v(1,1:4)
              y(1:4) = obj(iobj)%v(2,1:4)

              obj(iobj)%tsur(isid) = SP_GRID_BOUNDARY
              obj(iobj)%imap(1,isid) = iobj
              obj(iobj)%isur(1,isid) = isid
 
              SELECTCASE (isid)
                CASE(1)
                  DO iobj2 = MIN(nobj,iobj+1), MIN(nobj,iobj+nxbin+1) 
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(1)-x2(4)).LT.DTOL.AND.
     .                  DABS(y(1)-y2(4)).LT.DTOL.AND.
     .                  DABS(x(2)-x2(3)).LT.DTOL.AND.
     .                  DABS(y(2)-y2(3)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 3                   
                      EXIT
                    ENDIF
                  ENDDO

                CASE(2)
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(2)-x2(1)).LT.DTOL.AND.
     .                  DABS(y(2)-y2(1)).LT.DTOL.AND.
     .                  DABS(x(3)-x2(4)).LT.DTOL.AND.
     .                  DABS(y(3)-y2(4)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 4                    
                      EXIT
                    ENDIF
                  ENDDO

                CASE(3)
                  DO iobj2 = MAX(fobj,iobj-1),MAX(fobj,iobj-nxbin-1),-1
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(3)-x2(2)).LT.DTOL.AND.
     .                  DABS(y(3)-y2(2)).LT.DTOL.AND.
     .                  DABS(x(4)-x2(1)).LT.DTOL.AND.
     .                  DABS(y(4)-y2(1)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 1                  
                      EXIT
                    ENDIF
                  ENDDO               

                CASE(4)
                  DO iobj2 = MIN(nobj,iobj+1),MIN(nobj,iobj+nxbin+1)
                    x2(1:4) = obj(iobj2)%v(1,1:4)
                    y2(1:4) = obj(iobj2)%v(2,1:4)
                    IF (DABS(x(1)-x2(2)).LT.DTOL.AND.
     .                  DABS(y(1)-y2(2)).LT.DTOL.AND.
     .                  DABS(x(4)-x2(3)).LT.DTOL.AND.
     .                  DABS(y(4)-y2(3)).LT.DTOL) THEN
                      obj(iobj)%tsur(isid) = SP_GRID_SURFACE
                      obj(iobj)%imap(1,isid) = iobj2
                      obj(iobj)%isur(1,isid) = 2                   
                      EXIT
                    ENDIF
                  ENDDO
                CASEDEFAULT
                  CALL ER('BuildInversionMesh','Too many sides',*99)
              ENDSELECT       

            ENDDO
 
c            WRITE(0,*) 'CEL:',iobj,obj(iobj)%tsur(1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%imap(1,1:4)
c            WRITE(0,*) '   :',iobj,obj(iobj)%isur(1,1:4)

          ENDDO

        ENDIF


      ELSEIF ((ielement.NE.0.AND.
     .         opt%obj_option(ielement).EQ.1)) THEN
c     .        opt%ob_invgrd.EQ.1) THEN

c...    Toroidal approximation via discretization:

        istart = 1
        iend   = nsector 

c...    *HACK* (more hacks below)
        WRITE(0,*)
        WRITE(0,*) '--------------------------------------------------'
        WRITE(0,*) ' TOROIDAL INVERSION MESH HACK!'
        WRITE(0,*) '--------------------------------------------------'
        istart = 0
        iend   = 2
        ivol   = 0 

        IF (iend-istart+1.GT.nsector) 
     .    CALL ER('BuildInversionMesh','An excess of toroidal '//
     .            'sectors requested',*99)

        DO isector = istart, iend

          ang = REAL(isector - 1) * dangle

c...      *HACK* 
c          ivol = 0  ! Integration volume index 

          DO iy = nybin, 1, -1
            DO ix = 1, nxbin

              IF (nobj+1.GT.MAX3D) 
     .          CALL ER('BuildObjects','Insufficient array bounds '//
     .                  'for all objects',*98)     

              ivol = ivol + 1

              nobj = nobj + 1
              obj(nobj)%index       = ielement  ! nobj
              obj(nobj)%type        = OP_INTEGRATION_VOLUME
              obj(nobj)%mode        = 0      
              obj(nobj)%surface     = 1      ! SOLID
              obj(nobj)%phi         = ang
              obj(nobj)%wedge1      = 0
              obj(nobj)%wedge2      = 0
              obj(nobj)%colour      = 3
              obj(nobj)%orientation = 1      ! CW
              obj(nobj)%ik          = ix
              obj(nobj)%ir          = iy
              obj(nobj)%in          = 0
              obj(nobj)%ivolume     = ivol
              obj(nobj)%nsur        = 6
              obj(nobj)%gsur(1:6)   = GT_TD
              obj(nobj)%nver        = 8

c             obj(nobj)%quantity = 1.0

c...          *HACK* 
              IF (iy.EQ.nybin/2.AND.
     .            ix.GT.nxbin/3.AND.ix.LT.nxbin-nxbin/3+1.AND.
     .            isector.EQ.1) THEN
                obj(nobj)%quantity(1) = 1.0
              ELSE
                obj(nobj)%quantity(1) = 0.0
              ENDIF

c              IF (MOD(ix,2).EQ.0) THEN
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 1.0
c                ELSE
c                  obj(nobj)%quantity(1) = 0.0
c                ENDIF
c              ELSE
c                IF (MOD(iy,2).EQ.0) THEN 
c                  obj(nobj)%quantity(1) = 0.0
c                ELSE
c                  obj(nobj)%quantity(1) = 1.0
c                ENDIF
c              ENDIF
c              obj(3)%quantity(1) = 10.0
c            obj(nobj)%quantity(1) = REAL(ivol)
c            obj(nobj)%quantity(2) = 1.0

              obj(nobj)%nmap = 0

c...          ...
              p1(1,1)   = xorigin + xdelta * DBLE(ix)   ! Clean this up...
              p1(2,1)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,1)   = DBLE(p1(1,1))*DTAN(DBLE(-0.5*dangle))
              p1(1,1+4) = p1(1,1)
              p1(2,1+4) = p1(2,1)
              p1(3,1+4) = DBLE(p1(1,1))*DTAN(DBLE(+0.5*dangle))

              p1(1,2)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,2)   = yorigin + ydelta * DBLE(iy-1)
              p1(3,2)   = DBLE(p1(1,2))*DTAN(DBLE(-0.5*dangle))
              p1(1,2+4) = p1(1,2)
              p1(2,2+4) = p1(2,2)
              p1(3,2+4) = DBLE(p1(1,2))*DTAN(DBLE(+0.5*dangle))

              p1(1,3)   = xorigin + xdelta * DBLE(ix-1)
              p1(2,3)   = yorigin + ydelta * DBLE(iy)
              p1(3,3)   = DBLE(p1(1,3))*DTAN(DBLE(-0.5*dangle))
              p1(1,3+4) = p1(1,3)
              p1(2,3+4) = p1(2,3)
              p1(3,3+4) = DBLE(p1(1,3))*DTAN(DBLE(+0.5*dangle))

              p1(1,4)   = xorigin + xdelta * DBLE(ix)
              p1(2,4)   = yorigin + ydelta * DBLE(iy)
              p1(3,4)   = DBLE(p1(1,4))*DTAN(DBLE(-0.5*dangle))
              p1(1,4+4) = p1(1,4)
              p1(2,4+4) = p1(2,4)
              p1(3,4+4) = DBLE(p1(1,4))*DTAN(DBLE(+0.5*dangle))
c... 
              DO i1 = 1, 8
                obj(nobj)%v(1:3,i1) = p1(1:3,i1)
              ENDDO
c...          Rotate vertices:
              DO i1 = 1, 8
                x1 = obj(nobj)%v(1,i1)
                z1 = obj(nobj)%v(3,i1)
                obj(nobj)%v(1,i1) = DCOS(DBLE(ang)) * x1 -
     .                              DSIN(DBLE(ang)) * z1
                obj(nobj)%v(3,i1) = DSIN(DBLE(ang)) * x1 +
     .                              DCOS(DBLE(ang)) * z1
              ENDDO
c...          Surface 1:
              IF (isector.EQ.istart) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(1) = SP_GRID_SURFACE
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj + ncells * (nsector - 1)
                  obj(nobj)%isur(1,1) = 6
                ELSE
                  obj(nobj)%tsur(1) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(1) = 1
                  obj(nobj)%imap(1,1) = nobj
                  obj(nobj)%isur(1,1) = 1
                ENDIF
              ELSE
                obj(nobj)%tsur(1) = SP_GRID_SURFACE
                obj(nobj)%nmap(1) = 1
                obj(nobj)%imap(1,1) = nobj - ncells
                obj(nobj)%isur(1,1) = 6
              ENDIF
              obj(nobj)%npts(1) = 4
              obj(nobj)%ipts(1,1) = 1
              obj(nobj)%ipts(2,1) = 2
              obj(nobj)%ipts(3,1) = 3
              obj(nobj)%ipts(4,1) = 4
c...          Surface 2:
              IF (iy.EQ.1) THEN     ! *** THIS IS NOT TRUE IF TARGET TRANSPARENT!
                obj(nobj)%tsur(2) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj
                obj(nobj)%isur(1,2) = 2
              ELSE
                obj(nobj)%tsur(2) = SP_GRID_SURFACE
                obj(nobj)%nmap(2) = 1
                obj(nobj)%imap(1,2) = nobj + nxbin
                obj(nobj)%isur(1,2) = 4
              ENDIF
              obj(nobj)%npts(2) = 4
              obj(nobj)%ipts(1,2) = 2   ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,2) = 1
              obj(nobj)%ipts(3,2) = 5
              obj(nobj)%ipts(4,2) = 6
c...          Surface 3:
              IF (ix.EQ.1) THEN
                obj(nobj)%tsur(3) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj
                obj(nobj)%isur(1,3) = 3
              ELSE
                obj(nobj)%tsur(3) = SP_GRID_SURFACE
                obj(nobj)%nmap(3) = 1
                obj(nobj)%imap(1,3) = nobj - 1
                obj(nobj)%isur(1,3) = 5
              ENDIF
              obj(nobj)%npts(3) = 4
              obj(nobj)%ipts(1,3) = 3
              obj(nobj)%ipts(2,3) = 2
              obj(nobj)%ipts(3,3) = 6
              obj(nobj)%ipts(4,3) = 7
c...          Surface 4:
              IF (iy.EQ.nybin) THEN   
                obj(nobj)%tsur(4) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj
                obj(nobj)%isur(1,4) = 4
              ELSE
                obj(nobj)%tsur(4) = SP_GRID_SURFACE
                obj(nobj)%nmap(4) = 1
                obj(nobj)%imap(1,4) = nobj - nxbin
                obj(nobj)%isur(1,4) = 2
              ENDIF
              obj(nobj)%npts(4) = 4
              obj(nobj)%ipts(1,4) = 4    ! If this is right, then the rotation is backwards! 
              obj(nobj)%ipts(2,4) = 3
              obj(nobj)%ipts(3,4) = 7
              obj(nobj)%ipts(4,4) = 8
c...          Surface 5:
              IF (ix.EQ.nxbin) THEN
                obj(nobj)%tsur(5) = SP_GRID_BOUNDARY
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj
                obj(nobj)%isur(1,5) = 5
              ELSE
                obj(nobj)%tsur(5) = SP_GRID_SURFACE
                obj(nobj)%nmap(5) = 1
                obj(nobj)%imap(1,5) = nobj + 1
                obj(nobj)%isur(1,5) = 3
              ENDIF
              obj(nobj)%npts(5) = 4
              obj(nobj)%ipts(1,5) = 1
              obj(nobj)%ipts(2,5) = 4
              obj(nobj)%ipts(3,5) = 8
              obj(nobj)%ipts(4,5) = 5
c...          Surface 6:
              IF (isector.EQ.iend) THEN
                IF (iend-istart+1.EQ.nsector) THEN   
c...              Full torus:
                  obj(nobj)%tsur(6) = SP_GRID_SURFACE
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj - ncells * (nsector - 1)
                  obj(nobj)%isur(1,6) = 1
                ELSE
                  obj(nobj)%tsur(6) = SP_GRID_BOUNDARY
                  obj(nobj)%nmap(6) = 1
                  obj(nobj)%imap(1,6) = nobj
                  obj(nobj)%isur(1,6) = 6
                ENDIF
              ELSE
                obj(nobj)%tsur(6) = SP_GRID_SURFACE
                obj(nobj)%nmap(6) = 1
                obj(nobj)%imap(1,6) = nobj + ncells
                obj(nobj)%isur(1,6) = 1
              ENDIF
              obj(nobj)%npts(6) = 4
              obj(nobj)%ipts(1,6) = 8
              obj(nobj)%ipts(2,6) = 7
              obj(nobj)%ipts(3,6) = 6
              obj(nobj)%ipts(4,6) = 5
            ENDDO
          ENDDO
 
        ENDDO

      ELSE
        CALL ER('BuildInversionMesh','Unknown option',*99)
      ENDIF

c...  Local mesh refinement adjustment based on mask, from a previous 
c     iteration:


c...  Generalized cell mapping method:

c      WRITE(0,*) '???',obj(8)%gsur(1)
      WRITE(6,*) 'DESPERATE'
      DO iobj = 1, nobj
        DO i1 = 1, 4
          i2 = i1 + 1
          IF (i2.GT.4) i2 = 1
          WRITE(6,*) iobj,i1,obj(iobj)%v(1:2,i1)
          WRITE(6,*) iobj,i1,obj(iobj)%v(1:2,i2)
        ENDDO
      ENDDO

 98   RETURN
 99   STOP
      END

