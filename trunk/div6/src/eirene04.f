c     -*Fortran*-
c
c ======================================================================
c
c OSM-EIRENE04 interface:
c
c Write input file
c Write triangle
c     -WriteEIRENE_04
c
c
c
c Read triangle
c Move data to regular grid
c
c
c ======================================================================
c
c subroutine: ProcessTriangles_04
c
      SUBROUTINE ProcessTriangles_04(mode)
      USE mod_eirene04
      IMPLICIT none

      INTEGER mode

      LOGICAL PointOnLine

      REAL       TOL        ,DTOL
c      PARAMETER (TOL=1.0E-05,DTOL=1.0D-07)
      PARAMETER (TOL=1.0E-06,DTOL=1.0D-07)

      INTEGER i1,i2,i3,v1,v2,v3,v4,knot,ring,side,target,
     .        xupdate(10),yupdate(10),ix,iy,iscan
      LOGICAL test,output
      REAL    xmin,xmax,ymin,ymax,xval,yval
      REAL*8  x(0:2),y(0:2),s,t

      INTEGER, ALLOCATABLE :: xregion(:),yregion(:),nregion(:,:),
     .                        iregion(:,:,:)

      DATA (xupdate(i1),i1=1,10) /0, -1,  0,  1, -1, 1, -1, 0, 1, 0/ , 
     .     (yupdate(i1),i1=1,10) /0, -1, -1, -1,  0, 0,  1, 1, 1, 0/
  
      WRITE(0,*) 'PROCESSING TRIANGLES'  

      output = .FALSE.

      IF (mode.EQ.-1) GOTO 10

      WRITE(0,*) '  REMOVING DUPLICATE VERTICIES'

c...  Eliminate duplicate verticies:    ! SPEED:? SORT VERTICIES INTO REGIONS AND ONLY SCAN OVER NEIGHBOUR REGIONS? 
      DO i1 = 1, nver
        DO i2 = i1+1, nver
          IF (ver(i1,1).NE.-999.0.AND.
     .        ABS(ver(i1,1)-ver(i2,1)).LT.TOL.AND.
     .        ABS(ver(i1,2)-ver(i2,2)).LT.TOL) THEN
            ver(i2,1) = -999.0
            ver(i2,2) = -999.0
c...        Search all triangles for the vertex to be removed and update:
            DO i3 = 1, ntri
              DO v1 = 1, 3
                IF (tri(i3)%ver(v1).EQ.i2) tri(i3)%ver(v1) = i1
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO
c...  Remove vertices tagged for deletion:
      DO i1 = nver, 1, -1
        IF (ver(i1,1).EQ.-999.0) THEN
          DO i2 = i1, nver-1
            ver(i2,1) = ver(i2+1,1)
            ver(i2,2) = ver(i2+1,2)
          ENDDO
          DO i2 = 1, ntri
            DO v1 = 1, 3
              IF (tri(i2)%ver(v1).GE.i1) 
     .          tri(i2)%ver(v1) = tri(i2)%ver(v1) - 1
            ENDDO
          ENDDO
          nver = nver - 1
        ENDIF
      ENDDO

c...  Check if 2 close lying points got merged by accident, and if 
c     yes, then remove the triangle:
      DO i1 = ntri, 1, -1
        DO v1 = 1, 3   
          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) THEN
            WRITE(0,*) 'KILLING TRIANLGE-FIX!:',i1
            DO i2 = i1, ntri-1
              tri(i2) = tri(i2+1)
            ENDDO
            ntri = ntri - 1
            EXIT
          ENDIF
        ENDDO
      ENDDO


 10   CONTINUE

      WRITE(0,*) '  ASSIGNING REGIONS' 
      ALLOCATE(xregion(ntri)) ! Move to triangle construct 
      ALLOCATE(yregion(ntri))
      ALLOCATE(nregion(10,10)) 
      ALLOCATE(iregion(10,10,ntri))  ! I used to divide NTRI by 2 or 3 to reduce the array size, but this causes
                                     ! problems if the grid is very dense in a particular region of interest 

c...  Find horizontal and vertical extents of the grid:
      xmin =  1.0E+31
      xmax = -1.0E+31
      ymin =  1.0E+31
      ymax = -1.0E+31
      DO i1 = 1, ntri      
        DO v1 = 1, 3
          xmin = MIN(xmin,ver(tri(i1)%ver(v1),1))
          xmax = MAX(xmax,ver(tri(i1)%ver(v1),1))
          ymin = MIN(ymin,ver(tri(i1)%ver(v1),2))
          ymax = MAX(ymax,ver(tri(i1)%ver(v1),2))
        ENDDO
      ENDDO
      xmin = xmin - 0.01
      xmax = xmax + 0.01
      ymin = ymin - 0.01
      ymax = ymax + 0.01
c...  Bin me baby:
      v1 = 1
      DO i1 = 1, ntri      
        xval = ver(tri(i1)%ver(v1),1)  ! ***BETTER TO USE TRIANGLE CENTERS, BUT NOT STORED YET... 
        yval = ver(tri(i1)%ver(v1),2)
        xregion(i1) = INT((xval - xmin) / (xmax - xmin) * 10.0) + 1
        yregion(i1) = INT((yval - ymin) / (ymax - ymin) * 10.0) + 1
      ENDDO
c...  Build list of regions:
      WRITE(0,*) '  BUILDING REGION LISTS' 
      nregion = 0
      iregion = 0
      DO i1 =  1, ntri
        ix = xregion(i1)
        iy = yregion(i1)
        nregion(ix,iy)                = nregion(ix,iy) + 1
        iregion(ix,iy,nregion(ix,iy)) = i1
      ENDDO
      DO i1 = 10, 1, -1
        WRITE(0,'(2X,10I5)') (nregion(i2,i1),i2=1,10)
      ENDDO



      WRITE(0,*) '  BUILDING CONNECTION MAP' 

c...  Build connection map:
      DO i1 = 1, ntri

        DO v1 = 1, 3

          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1

c...      Find matching cell based on vertex coordinates:
          ix = xregion(i1)
          iy = yregion(i1)
          iscan = 0
          i3 = 0
          DO WHILE (iscan.LE.10)

c...        Advance the loop index 'i2' (sorry that this is convoluted, but it needs 
c           to be fast 'cause some grids are big):
            i3 = i3 + 1
            IF     (iscan.EQ.0) THEN
c...          Check if a map is already assigned, and if yes, look there 
c             first (iscan=0):
              IF (tri(i1)%map(v1).GT.0) THEN
                i2 = tri(i1)%map(v1)
              ELSE
                i3 = 0
                iscan = 1
                CYCLE
              ENDIF
            ELSEIF (iscan.EQ.10) THEN
c...          Map not found so search the whole triangle mesh to be sure that 
c             there isn't one -- a large triangle can cause problems for 0<iscan<10 (iscan=10):
              IF (i2.LT.ntri) THEN
                i2 = i3
              ELSE
c...            All triangles have been searched but no mapping found, trigger exit condition:
                iscan = 999
                CYCLE
              ENDIF
            ELSE
c...          First, search the region that the triangle is in (iscan=1), and 
c             then the neighbouring regions:
              IF (i3.LE.nregion(ix,iy)) THEN
                i2 = iregion(ix,iy,i3)
              ELSE 
c...            Change search region, making sure the new region is valid:
                i3 = 0
                ix = 0
                DO WHILE (ix.LT.1.OR.ix.GT.10.OR.
     .                    iy.LT.1.OR.iy.GT.10)
                  iscan = iscan + 1
                  ix = xregion(i1) + xupdate(iscan)
                  iy = yregion(i1) + yupdate(iscan)
                ENDDO        
                CYCLE
              ENDIF
            ENDIF

            IF (i1.EQ.i2) CYCLE

c            IF (i1.EQ.2) THEN
c              WRITE(0,'(A,8I6)') ' SEARCH:',
c     .          iscan,i1,v1,i2,xregion(i1),xregion(i2),
c     .          yregion(i1),yregion(i2)
c            ENDIF

            DO v3 = 1, 3
              v4 = v3 + 1        
              IF (v3.EQ.3) v4 = 1
              IF ((tri(i1)%ver(v1).EQ.tri(i2)%ver(v3).AND.
     .             tri(i1)%ver(v2).EQ.tri(i2)%ver(v4)).OR.
     .            (tri(i1)%ver(v1).EQ.tri(i2)%ver(v4).AND.
     .             tri(i1)%ver(v2).EQ.tri(i2)%ver(v3))) THEN

c...            Adding mapping:
                tri(i1)%map(v1) = i2
                tri(i1)%sid(v1) = v3
                tri(i2)%map(v3) = i1
                tri(i2)%sid(v3) = v1

c...            Neighbour has been found, trigger exit condtion:
                iscan = 999
                EXIT

              ENDIF
            ENDDO

          ENDDO

        ENDDO
      ENDDO


      WRITE(0,*) '  MAPPING SIDES TO SURFACES' 

c...  Map triangles to surfaces:     
      DO i1 = 1, ntri
        DO v1 = 1, 3
          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1          

          DO i2 = 1, nadd            
            IF     (add(i2)%type.EQ.NON_DEFAULT_STANDARD) THEN

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN

                knot   = tri(i1)%index(1)  ! Parameter ...?
                ring   = tri(i1)%index(2)  ! Parameter ...?
                side   = tri(i1)%sideindex(1,v1)
                target = tri(i1)%sideindex(2,v1)

                IF     (add(i2)%subtype .EQ.STRATUM.AND.
     .                  add(i2)%index(1).LE.ring   .AND.
     .                  add(i2)%index(2).GE.ring   .AND.
     .                  add(i2)%index(3).EQ.target) THEN
                  tri(i1)%map(v1) = 0 ! This should only be set to 0 if the target is opaque...
                  tri(i1)%sid(v1) = 0 ! ditto
                  tri(i1)%sur(v1) = add(i2)%num

                ELSEIF (add(i2)%subtype .EQ.MAGNETIC_GRID_BOUNDARY) THEN
c                  i3 = 1
c                  DO WHILE (add(i2)%index(i3).NE.0)
                    IF (add(i2)%index(1).LE.knot.AND.
     .                  add(i2)%index(2).GE.knot.AND.
     .                  add(i2)%index(3).EQ.ring.AND.
     .                  add(i2)%index(4).EQ.side) THEN
                      tri(i1)%sur(v1) = add(i2)%num 
                      EXIT
                    ENDIF
c                    i3 = i3 + 2
c                  ENDDO

                ELSEIF (add(i2)%subtype .EQ.ADDITIONAL) THEN ! *** WRONG PLACE ***?
                ENDIF

              ENDIF

            ELSEIF (add(i2)%type.EQ.VESSEL_WALL) THEN
              test = .TRUE.  ! ***REMOVE***
c...          Assign surface end points:
              x(0) = DBLE(add(i2)%v(1,1))
              y(0) = DBLE(add(i2)%v(2,1))
              x(1) = DBLE(add(i2)%v(1,2))
              y(1) = DBLE(add(i2)%v(2,2))
c...          Side vertex 1:
              x(2) = DBLE(ver(tri(i1)%ver(v1),1))
              y(2) = DBLE(ver(tri(i1)%ver(v1),2))
              test = test.AND.PointOnLine(x,y,s,t,1,output)
c...          Side vertex 2:
              IF (test) THEN
                x(2) = DBLE(ver(tri(i1)%ver(v2),1))
                y(2) = DBLE(ver(tri(i1)%ver(v2),2))
                test = test.AND.PointOnLine(x,y,s,t,1,output)
c...            Assign surface index, as appropriate:
                IF (test) THEN
                  IF (tri(i1)%map(v1).EQ.-1) THEN
                    tri(i1)%map(v1) = 0
                    tri(i1)%sid(v1) = 0
                  ENDIF
c...              Need to identify which non-default standard surface
c                 should be defined: ... 
                  IF (add(i2)%index(3).NE.0) THEN
                    tri(i1)%sur(v1) = add(i2)%index(3)
                    tri(i1)%sideindex(3,v1)=add(i2)%index(1)  ! Store xVESM index of surface
                    tri(i1)%sideindex(4,v1)=add(i2)%index(2)  ! Store additional surface index 
c                    tri(i1)%sur(v1) = add(add(i2)%index(3))%num
                  ELSE
                    CALL ER('...','Surface index not assigned?1',*99)
c                    tri(i1)%sur(v1) = 4 ! Temp
                  ENDIF
                ENDIF
              ENDIF

            ELSE
              CALL ER('ProcessTriangles_04','Invalid surface type',*99)
            ENDIF
          ENDDO
        ENDDO
      ENDDO


      WRITE(0,*) '  CHECKING ASSIGNMENTS' 

      DO i1 = 1, ntri
        DO v1 = 1, 3                             ! Perhaps eliminate maps through solid surfaces? 
          IF (tri(i1)%map(v1).EQ.0) THEN
            IF ( tri(i1)%sur(v1).EQ.0.OR.
     .          (tri(i1)%sur(v1).NE.0.AND.
     .           add(MAX(1,tri(i1)%sur(v1)))%iliin.EQ.-1)) THEN 
c            IF (tri(i1)%map(v1).EQ.-1) THEN 
c...          If this shows up again, it may be related to DTOL in PointOnLine:
              WRITE(0,*) 'PROBLEMS WITH MAP',i1
              CALL WriteTriangleFiles
              CALL DumpGrid('PROBLEM WITH TRIANGLE MAP')
            ENDIF
          ENDIF
c...      Check if 2 close lying points got merged by accident:
          v2 = v2 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) 
     .      CALL ER('ProcessTriangles_04','2 sided triangle',*98)

        ENDDO
      ENDDO

      WRITE(0,*) 'DONE'

      DEALLOCATE(xregion) ! Move to triangle construct since this could be useful elsewhere...?
      DEALLOCATE(yregion)
      DEALLOCATE(nregion)
      DEALLOCATE(iregion)

      RETURN
 98   WRITE(0,*) 'TRI:',i1,v1
 99   STOP
      END
c
c ======================================================================
c
c subroutine: ReadPolyFile
c
      SUBROUTINE ReadPolyFile(zone)
      USE mod_eirene04
      IMPLICIT none

      INTEGER zone,fp,i1,i2,idum1,idum2

c
c     Need to work from a list of additional surfaces that are loaded from the fluid
c     code and from the existing triangle data set (the pre-defined sides), rather than
c     from the fluid code grid as I have done here.  Perhaps just lost a selection 
c     of line segments and order them, then build the domain.  The additional surfaces
c     will have to be assigned a volume index in the input file -- that is how they will be
c     selected... 
c
c     Only additional surfaces that are not part of a particular zone are passed to EIRENE
c     as additional surfaces, the rest go as non-default standard surfaces (with reflection
c     models specified...  
c
c     Additional surface has: OSM index, EIRENE zone index, non-default surface index (?) so
c     that an additional surface that is part of a zone is not simply assigned as some 
c     anonymouns part of the zone wall (?), ... should every wall segment be assigned its own
c     non-default standard index, and then this used to map data back and forth? 
c
c     Need to assign a non-default standard surface for each region of IRWALL...
c
c     Need to generalize the reading of target data on the EIRENE side so that I have control 
c     over strata definitions... 
c
c
c

      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.1.ele',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)      
      READ(fp,*) idum1
      DO i1 = 1, idum1
        READ(fp,*,ERR=98,END=98) idum2,(tri(i1+ntri)%ver(i2),i2=1,3)
        DO i2 = 1, 3
          tri(i1+ntri)%ver(i2) = tri(i1+ntri)%ver(i2) + nver
          tri(i1+ntri)%map(i2) = -1 
          tri(i1+ntri)%sid(i2) = -1
          tri(i1+ntri)%sur(i2) =  0
        ENDDO
        tri(i1+ntri)%type = VACUUM_GRID
        tri(i1+ntri)%zone = zone
      ENDDO
      CLOSE (fp) 
      ntri = ntri + idum1

      fp = 99      
      OPEN(UNIT=fp,FILE='triangle.1.node',ACCESS='SEQUENTIAL',
     .     STATUS='OLD',ERR=98)      
      READ(fp,*) idum1
      DO i1 = 1, idum1
        READ(fp,*,ERR=98,END=98) idum2,ver(i1+nver,1),ver(i1+nver,2)
      ENDDO
      CLOSE (fp) 
      nver = nver + idum1


c...  Build connection map:

c...  Associate with non-default standard surface and apply reflection models
c     from former additional surfaces: 

c...  Eliminate redundant verticies:


      RETURN
 98   CALL ER('ReadPolyFile','Problems with file access',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: DefineEireneSurfaces
c
      INTEGER FUNCTION NewEireneSurface(type)
      USE mod_eirene04
      IMPLICIT none

      INTEGER type,index,
     .        defiliin,defilside,defilswch,defiltor,defilcol,defilcell
      REAL    defrecyct,defrecycf

      REAL       ECH
      PARAMETER (ECH=1.602E-19)

c...  Assign defaults to surface properties:
      defiliin  = 1
      defilside = 0
      defilswch = 0
      defiltor  = 0
      defilcell = 0
      defilcol  = 1
      defrecyct = 1.0
      defrecycf = 1.0

      IF     (type.EQ.VESSEL_WALL) THEN

        nadd = nadd + 1
        add(nadd)%type   = type
        add(nadd)%subtype = -1
        add(nadd)%surtxt   = '* vessel wall (default)'
        add(nadd)%num    = 0
        add(nadd)%iliin  = defiliin
        add(nadd)%ilside = defilside
        add(nadd)%ilswch = defilswch
        add(nadd)%iltor  = defiltor
        add(nadd)%ilcell = defilcell
        add(nadd)%ilcol  = defilcol
        add(nadd)%reflect = LOCAL
        add(nadd)%ewall = -wtemp * 1.38E-23 / ECH
        add(nadd)%material = wmater
        add(nadd)%recyct = defrecyct  ! 0.0  ! defrecyct  
        add(nadd)%recycf = defrecycf
c...    Assume a 2-point line segment:
        add(nadd)%nsur = 1
        add(nadd)%npts(1) = 2
        add(nadd)%ipts(1,1) = 1
        add(nadd)%ipts(2,1) = 2
        add(nadd)%nver = 2

      ELSEIF (type.EQ.NON_DEFAULT_STANDARD) THEN

        nadd = nadd + 1
        add(nadd)%type   = type
        add(nadd)%subtype = -1
        add(nadd)%surtxt  = '* non-default standard (default)'
        add(nadd)%num    = 0
        add(nadd)%iliin  = defiliin
        add(nadd)%ilside = defilside
        add(nadd)%ilswch = defilswch
        add(nadd)%iltor  = defiltor
        add(nadd)%ilcell = defilcell
        add(nadd)%ilcol  = defilcol
        add(nadd)%reflect = GLOBAL
        add(nadd)%ewall    = 0.0
        add(nadd)%material = 0.0
        add(nadd)%recyct = defrecyct
        add(nadd)%recycf = defrecycf
        add(nadd)%nsur = 0
        add(nadd)%nver = 0

      ELSE
        CALL ER('NewEireneSurface','Invalid type',*99)
      ENDIF

      NewEireneSurface = nadd

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: RefineTriangles
c
c     *** NOT CURRENTLY IN USE ***
c
      SUBROUTINE RefineTriangles
      USE mod_eirene04
      IMPLICIT none

      REAL*8 CalcTriangleArea

      INTEGER i1,i2,i3,iside,zone,iside2,itri2
      LOGICAL cont
      REAL    xcen,ycen
      REAL*8  x(3),y(3),area,sidelen,maxsidelen,limit

      RETURN

      limit = 0.0001
      zone = 3

c...  Loop over triangles:

      cont = .TRUE.
      DO WHILE (cont)
        cont = .FALSE.

        DO i1 = ntri, 1, -1

c...      Estimate triangle center (fix this up!):
          xcen = 0.0
          ycen = 0.0
          DO i2 = 1, 3
            xcen = xcen + 0.3333 * ver(tri(i1)%ver(i2),1)
            ycen = ycen + 0.3333 * ver(tri(i1)%ver(i2),2)
          ENDDO

c...      Selection rules:
          IF (tri(i1)%type.NE.VACUUM_GRID.OR.
     .        tri(i1)%zone.NE.zone.OR.
     .        xcen.LT. 0.75.OR.
     .        ycen.GT.-0.50) CYCLE         ! *HARD*

          x(1) = DBLE(ver(tri(i1)%ver(1),1))
          y(1) = DBLE(ver(tri(i1)%ver(1),2))
          x(2) = DBLE(ver(tri(i1)%ver(2),1))
          y(2) = DBLE(ver(tri(i1)%ver(2),2))
          x(3) = DBLE(ver(tri(i1)%ver(3),1))
          y(3) = DBLE(ver(tri(i1)%ver(3),2))

          area = CalcTriangleArea(x(1),y(1),x(2),y(2),x(3),y(3))

          WRITE(0,*) 'I1,AREA=',i1,ntri,area

c          STOP 'sdfsd'

          IF (area.GT.limit) THEN
c...        Split triangle in half:

            cont = .TRUE.

c...        Find best side -- the longest and where it is not neighbouring
c           a triangle of a different type or in a different zone:        
            iside = 0
            maxsidelen = 0.0
            DO i2 = 1, 3
              i3 = i2 + 1
              IF (i2.EQ.3) i3 = 1
              sidelen = DSQRT( (x(i2) - x(i3))**2 + (y(i2) - y(i3))**2)
              WRITE(0,*) '.......=',i1,i2,i3


              IF (tri(i1)%map(i2).LE.0.OR. ! *TEMP* is this okay
     .            tri(tri(i1)%map(i2))%type.EQ.VACUUM_GRID.AND.
     .            tri(tri(i1)%map(i2))%zone.EQ.zone.AND.            ! *HARD*
     .            sidelen.GT.maxsidelen) THEN            
                iside = i2
                maxsidelen = sidelen
              ENDIF
            ENDDO
            IF (iside.EQ.0) 
     .        CALL ER('RefineTriangles','No can find good side',*99)

            iside2 = tri(i1)%sid(iside)
            itri2  = tri(i1)%map(iside)

            WRITE(0,*) 'SPLITTING:',i1,iside
            WRITE(0,*) 'SPLITTING:',itri2,iside2

c...        Split original triangle:
            i2 = iside
            i3 = i2 + 1
            IF (i2.EQ.3) i3 = 1

c...        New vertex:
            nver = nver + 1
            ver(nver,1) = SNGL(0.5 * (x(i2) + x(i3)))
            ver(nver,2) = SNGL(0.5 * (y(i2) + y(i3)))
            ver(nver,3) = 0.0
c...        New triangle:
            ntri = ntri + 1
            tri(ntri)%type = VACUUM_GRID 
            tri(ntri)%zone = tri(i1)%zone
            DO i2 = 1, 3
              tri(ntri)%map(i2) = -1   ! * NEED TO DO BETTER *
              tri(ntri)%sid(i2) = -1
              tri(ntri)%sur(i2) =  0
            ENDDO
 
            IF     (iside.EQ.1) THEN
              tri(ntri)%ver(1) = nver
              tri(ntri)%ver(2) = tri(i1)%ver(2)
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(2) = nver              
            ELSEIF (iside.EQ.2) THEN
              tri(ntri)%ver(1) = tri(i1)%ver(1)
              tri(ntri)%ver(2) = nver
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(3) = nver              
            ELSEIF (iside.EQ.3) THEN
              tri(ntri)%ver(1) = nver
              tri(ntri)%ver(2) = tri(i1)%ver(2)
              tri(ntri)%ver(3) = tri(i1)%ver(3)
c...          Modify old triangle:
              tri(i1)%ver(3) = nver              
            ENDIF

c...        Split neigbouring triangle:

            IF (itri2.GT.0) THEN
c...          New triangle:
              ntri = ntri + 1
              tri(ntri)%type = VACUUM_GRID
              tri(ntri)%zone = tri(i1)%zone
              DO i2 = 1, 3
                tri(ntri)%map(i2) = -1 
                tri(ntri)%sid(i2) = -1
                tri(ntri)%sur(i2) =  0
              ENDDO
 
              IF     (iside2.EQ.1) THEN
                tri(ntri)%ver(1) = nver
                tri(ntri)%ver(2) = tri(itri2)%ver(2)
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(2) = nver              
              ELSEIF (iside2.EQ.2) THEN
                tri(ntri)%ver(1) = tri(itri2)%ver(1)
                tri(ntri)%ver(2) = nver
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(3) = nver              
              ELSEIF (iside2.EQ.3) THEN
                tri(ntri)%ver(1) = nver
                tri(ntri)%ver(2) = tri(itri2)%ver(2)
                tri(ntri)%ver(3) = tri(itri2)%ver(3)
c...            Modify old triangle:
                tri(itri2)%ver(3) = nver              
              ENDIF
            ENDIF

c            EXIT  ! * TEMP *

          ENDIF

        ENDDO

c... *TEMP*
        CALL ProcessTriangles_04(-1)

      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WritePolyFile
c
      SUBROUTINE WritePolyFile(eirntri,MAXNRS,eirtri)
      USE mod_eirene04
      IMPLICIT none

      INTEGER eirntri,MAXNRS
      REAL    eirtri(MAXNRS,20)
      
      INTEGER code

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      INTEGER fp,i1,i2,i3,i4,i5,v1,v2,id,ik,ir,npts,
     .        nseg,seg(0:1000,2),
     .        icnt,nhole
      LOGICAL zone,hole
      REAL    pts(1000,2),area
      REAL*8  x1,x2,y1,y2,len,t,tstep,xhole(50),yhole(50)
      character*256 command

      WRITE(0,*) 'HERE IN WRITEPOLYFILE'

      icnt = 0
c...CHECK BOUNDS ON SEG AND PTS AS THEY ARE ASSIGNED

      DO i1 = 1, eirntri
        IF (eirtri(i1,1).EQ.0.0) CYCLE

        nseg = 0
        npts = 0
        nhole = 0

c...    Build list of line segments for passing to TRIANGLE:
        DO i2 = i1+1, eirntri 
          IF (eirtri(i2,1).NE.0.0) EXIT          

c          WRITE(0,*) 'DATA:',i1,i2,eirtri(i2,2)

          IF     (eirtri(i2,2).EQ.1.0) THEN  ! PARAMETER
c...        Pull line segments from existing list of triangles, which initially 
c           will just be those from the triangularization of the standard fluid
c           code magnetic grid:
            DO i3 = 1, ntri
              DO v1 = 1, 3
                IF (tri(i3)%sur(v1).GE.NINT(eirtri(i2,3)).AND.
     .              tri(i3)%sur(v1).LE.NINT(eirtri(i2,4))) THEN 
                  nseg = nseg + 1
                  seg(nseg,1) = npts + 1
                  seg(nseg,2) = npts + 2
                  v2 = v1 + 1
                  IF (v1.EQ.3) v2 = 1
                  IF (.FALSE..AND.npts.EQ.0) THEN
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v2),1)  ! Switch orientation??? -- for some! 
                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v1),1)
                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
                  ELSE
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v1),1)  ! Switch orientation??? -- for some! 
                    pts(npts,2) = ver(tri(i3)%ver(v1),2)
                    npts = npts + 1
                    pts(npts,1) = ver(tri(i3)%ver(v2),1)
                    pts(npts,2) = ver(tri(i3)%ver(v2),2)
                  ENDIF
                ENDIF
              ENDDO              
            ENDDO

          ELSEIF (eirtri(i2,2).EQ.2.0.OR.eirtri(i2,2).EQ.3.0) THEN ! PARMETER
c...        Pull line segments from the list of additional surfaces: 

c            WRITE(0,*) 'I!:',i1

            DO i3 = 1, nadd
              IF (add(i3)%type.EQ.VESSEL_WALL.AND.
     .            ((eirtri(i2,2).EQ.2.0.AND.
     .              add(i3)%index(1).GE.NINT(eirtri(i2,3)).AND.
     .              add(i3)%index(1).LE.NINT(eirtri(i2,4))).OR.
     .             (eirtri(i2,2).EQ.3.0.AND.
     .              add(i3)%index(2).GE.NINT(eirtri(i2,3)).AND.
     .              add(i3)%index(2).LE.NINT(eirtri(i2,4))))) THEN

c                IF (eirtri(i2,2).EQ.2.0) 
c     . WRITE(0,*) 'TAKING:',add(i3)%index(1),i3

                x1 = DBLE(add(i3)%v(1,1))
                y1 = DBLE(add(i3)%v(2,1))
                x2 = DBLE(add(i3)%v(1,2))
                y2 = DBLE(add(i3)%v(2,2))

                len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)

                IF (eirtri(i1,3).GT.0.0.AND.
     .              len.GT.DBLE(eirtri(i1,3))) THEN
                  tstep = 1.0D0 / DBLE(INT(len/DBLE(eirtri(i1,3))) + 1)
                ELSE
                  tstep = 1.0D0
                ENDIF

c                WRITE(0,*) x1,y1
c                WRITE(0,*) x2,y2
                DO t = 0.0D0, 0.9999999D0, tstep 
                  nseg = nseg + 1
                  seg(nseg,1) = npts + 1
                  seg(nseg,2) = npts + 2
                  npts = npts + 1
                  pts(npts,1) = REAL(x1 + t * (x2 - x1))             ! Orientation is correct
                  pts(npts,2) = REAL(y1 + t * (y2 - y1)) 
                  npts = npts + 1
                  pts(npts,1) = REAL(x1 + (t + tstep) * (x2 - x1)) 
                  pts(npts,2) = REAL(y1 + (t + tstep) * (y2 - y1)) 
                ENDDO
              ENDIF
            ENDDO
 
          ELSEIF (eirtri(i2,2).EQ.4.0) THEN
c...        Holes:
            nhole = nhole + 1
            xhole(nhole) = eirtri(i2,3)
            yhole(nhole) = eirtri(i2,4)

          ELSE
            CALL ER('WritePolyFile','Invalid triangle grid segment',*99)
          ENDIF

        ENDDO


c...    Eliminate duplicate verticies:
        IF (.TRUE.) THEN
          DO i2 = 1, npts
            DO i3 = i2+1, npts
              IF (pts(i2,1).NE.-999.0.AND.
     .            ABS(pts(i2,1)-pts(i3,1)).LT.TOL.AND.
     .            ABS(pts(i2,2)-pts(i3,2)).LT.TOL) THEN
                pts(i3,1) = -999.0
                pts(i3,2) = -999.0
                DO i4 = 1, nseg
                  IF (seg(i4,1).EQ.i3) seg(i4,1) = i2
                  IF (seg(i4,2).EQ.i3) seg(i4,2) = i2
                ENDDO
              ENDIF
            ENDDO
          ENDDO
c...      Delete points:
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
c...      Sort segments:
          IF (.FALSE.) THEN
            DO i2 = 1, nseg-1
              DO i3 = i2+1, nseg
                IF     (ABS(pts(seg(i2,2),1)-
     .                      pts(seg(i3,1),1)).LT.TOL.AND.
     .                  ABS(pts(seg(i2,2),2)-
     .                      pts(seg(i3,1),2)).LT.TOL) THEN
                  IF (i3.EQ.i2+1) THEN
c...                Do nothing, all okay:
                    EXIT
                  ELSE
                    seg(0   ,1) = seg(i2+1,1)
                    seg(0   ,2) = seg(i2+1,2)
                    seg(i2+1,1) = seg(i3  ,1)
                    seg(i2+1,2) = seg(i3  ,2)
                    seg(i3  ,1) = seg(0   ,1)
                    seg(i3  ,2) = seg(0   ,2)
                    EXIT
                  ENDIF
                ELSEIF (ABS(pts(seg(i2,2),1)-
     .                      pts(seg(i3,2),1)).LT.TOL.AND.
     .                  ABS(pts(seg(i2,2),2)-
     .                      pts(seg(i3,2),2)).LT.TOL) THEN
                  IF (i3.EQ.i2+1) THEN
                    seg(0 ,1) = seg(i3,1)
                    seg(i3,1) = seg(i3,2)
                    seg(i3,2) = seg(0 ,1)
                    EXIT
                  ELSE
                    seg(0   ,1) = seg(i2+1,1)
                    seg(0   ,2) = seg(i2+1,2)
                    seg(i2+1,1) = seg(i3  ,2)
                    seg(i2+1,2) = seg(i3  ,1)
                    seg(i3  ,1) = seg(0   ,1)
                    seg(i3  ,2) = seg(0   ,2)
                    EXIT
                  ENDIF
                ENDIF
              ENDDO
              IF (i3.EQ.nseg+1) THEN
                WRITE(0,*) 'WARNING: SEGMENT GAP DETECTED'
              ENDIF
            ENDDO
          ENDIF

        ENDIF

        fp = 99      
        OPEN(UNIT=fp,FILE='triangle.poly',ACCESS='SEQUENTIAL',
     .       STATUS='REPLACE',ERR=99)      
        WRITE(fp,*) npts,2,1,0
        DO i2 = 1, npts
          WRITE(fp,'(I6,2F12.7)') i2,pts(i2,1),pts(i2,2)
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

c...    Call triangle:
        area = 0.01
        IF (eirtri(i1,3).NE.0.0) area = 0.5 * eirtri(i1,3)**2
        WRITE(command,10) 'triangle -p -q -a',area,' -Y triangle.poly'
 10     FORMAT(A,F10.8,A)
        WRITE(0,*) 'COMMAND: "'//command(1:LEN_TRIM(command))//'"'
        CALL CIssue(command(1:LEN_TRIM(command)),code)
        WRITE(0,*) 'RETURN_CODE:',code

        icnt = icnt + 1
c        IF (icnt.EQ.1) STOP 'STOP: CHECK POLY FILES'


        CALL ReadPolyFile(NINT(eirtri(i1,2)))

      ENDDO



      WRITE(0,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteTriangleFiles
c
      SUBROUTINE WriteTriangleFiles
      USE mod_eirene04
      IMPLICIT none

c *TEMP* -- for target data
c      INCLUDE 'params'
c      INCLUDE 'comtor'
c      INCLUDE 'cgeom'
c      INCLUDE 'slcom'

c      REAL GetMach,GetJsat,GetFlux 

      INTEGER fp,i1,i2,v1,ik1,ir1,it
c      INTEGER fp,i1,i2,v1,ik1,ir1,in,it,region,tarside
      LOGICAL found
      REAL    version
c      REAL    tarte,tarti,tarne,tarv,tarflux,tarisat,tarM

      REAL sumflux1,sumflux2

      REAL, ALLOCATABLE :: tdata(:)      

      WRITE(0,*) 'WRITING TRIANGLE FILES'

      version = 1.00

      fp = 99

      ALLOCATE(tdata(ntri))

      IF (photons.EQ.-1) THEN
c...    Load ionisation data from previous EIRENE call:
        STOP 'NEED TO MAKE COMPATIBLE WITH POSSIBLE BINARY FILE A'
        CALL LoadTriangleData(7,0,13,0,tdata,'default')  
      ELSE
        tdata = -999.0
      ENDIF

c...  Dump triangles (for OUT, not EIRENE):
      OPEN(UNIT=fp,FILE='triangles.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,3(2F10.6,2X))')
     .    i1,(ver(tri(i1)%ver(i2),1),ver(tri(i1)%ver(i2),2),i2=1,3)
      ENDDO
      CLOSE(fp)      


c...  Dump vertices:
      OPEN(UNIT=fp,FILE='triangles.npco_char',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triangles.points',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) nver
      DO i1 = 1, nver
        WRITE(fp,'(I6,3F12.6)') i1,ver(i1,1)*100.0,ver(i1,2)*100.0,0.0
      ENDDO
      CLOSE(fp)      

c...  Dump sides:
      OPEN(UNIT=fp,FILE='triangles.elemente',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triangles.sides',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,3(4X,3I6),4X,2I6)') 
     .    i1,(tri(i1)%ver(v1),v1=1,3),
     .    (tri(i1)%sideindex(1,v1),v1=1,3),
     .    (tri(i1)%sideindex(2,v1),v1=1,3),
     .    tri(i1)%index(1),tri(i1)%index(2)
      ENDDO
      CLOSE(fp)      

c...  Dump connection map:
      OPEN(UNIT=fp,FILE='triangles.neighbors',ACCESS='SEQUENTIAL',
c      OPEN(UNIT=fp,FILE='triangles.map',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,4X,3(3I6,4X),2I6,2X,2I4)') i1,
     .    (tri(i1)%map(v1),tri(i1)%sid(v1),tri(i1)%sur(v1),v1=1,3),
     .    tri(i1)%index(1),tri(i1)%index(2),
     .    tri(i1)%type,tri(i1)%zone
      ENDDO
      CLOSE(fp)      

c...  Dump plasma data:
      OPEN(UNIT=fp ,FILE='triangles.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A,F4.2,A)') 
     .  '* VERSION ',version,' OF '//
     .  fluid_code(1:LEN_TRIM(fluid_code))//
     .  ' PLASMA FILE FOR TRIANGULAR GRID'
c      WRITE(fp,'(A)') '* VERSION 1.0 OSM PLASMA FILE FOR TRIANGULAR '//
c     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '* Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A7,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*      ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I7,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             4E12.4,0P,6X,2I4)') i1,
     .    tri(i1)%plasma(1),           ! Te (eV)
     .    tri(i1)%plasma(2),           ! Ti (eV)
     .    tri(i1)%plasma(3)*1.0E-06,   ! ne (cm-3)
     .    tri(i1)%plasma(4)*100.0,     ! vx (cm s-1)
     .    tri(i1)%plasma(5)*100.0,     ! vy 
     .    tri(i1)%plasma(6)*100.0,     ! vz
     .    tri(i1)%bfield(1),           ! Bx (Tesla)
     .    tri(i1)%bfield(2),           ! By
     .    tri(i1)%bfield(3),           ! Bz
     .    tdata(i1),                   ! Ionisation rate from previous run (w or w/o photons)...
     .    tri(i1)%index(1),tri(i1)%index(2)
      ENDDO

      sumflux1= 0.0
      sumflux2= 0.0

c...  Target data:
      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,*) ntardat
      DO i1 = 1, ntri
        IF (tri(i1)%type.NE.MAGNETIC_GRID) CYCLE
        ik1 = tri(i1)%index(1)
        ir1 = tri(i1)%index(2) 
        DO v1 = 1, 3
          IF (v1.EQ.2.OR.tri(i1)%sur(v1).EQ.0) CYCLE

c...      Find corresponding target data, as set in ProcessFluidCode, based
c         on the fluid code cell/ring indeces:
          found = .FALSE.
          DO it = 1, ntardat
            IF (tardat(it,2).NE.REAL(ik1).OR.
     .          tardat(it,3).NE.REAL(ir1)) CYCLE

            IF (.NOT.found) THEN
              found = .TRUE.         

c              IF (tardat(it,7).LT.0.0) sumflux1 = sumflux1 + tardat(it,7)  ! ...debugging...
c              IF (tardat(it,7).GT.0.0) sumflux2 = sumflux2 + tardat(it,7)

              WRITE(fp,'(I7,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
     .                   F6.2,1P,E10.2,0P,6X,2I4)') 
c...            Target quantities:
     .          i1,v1,                  ! Triangle index and side index
     .          tardat(it,7 ),          ! ion flux to surface for species 1 (Amps)
     .          tardat(it,6 ),          ! Te (eV)
     .          tardat(it,8 ),          ! Ti (ev)
     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
     .          tardat(it,11),          ! Mach no.        (not read)
     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
     .          ik1,ir1                 ! Fluid grid indeces, for debugging only

            ELSE
              CALL ER('WriteTriangleFiles','Target data appears to  '//
     .                'be over-specified',*99)
            ENDIF
          ENDDO

        ENDDO
      ENDDO
      CLOSE(fp)      

c      WRITE(0,*) 'SUMFLUX:',sumflux1,sumflux2
c      STOP 'sdfsd'

      OPEN(UNIT=fp,FILE='triangles.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM E-FIELD FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,3A12)') 
     .  '* Index','Ex','Ey','Ez'
      WRITE(fp,'(A7,3A12)')
     .  '*      ','V m-1','V m-1','V m-1'
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
c...    Dump efield data:
        WRITE(fp,'(I7,1P,3E12.4,10X,0P,2I4)') i1,
     .    tri(i1)%efield(1),
     .    tri(i1)%efield(2),
     .    tri(i1)%efield(3),
     .    tri(i1)%index(1),tri(i1)%index(2) 
      ENDDO
      CLOSE(fp)      

      WRITE(0,*) 'DONE'

      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)

      RETURN
96    WRITE(0,*) 'WRITETRIANGEFILES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END

c
c ======================================================================
c
c subroutine: NewTriangle
c
      SUBROUTINE AssignPlasmaQuantities
      USE mod_eirene04
      IMPLICIT none

      INTEGER i1,i2

      DO i1 = 1, ntri
        DO i2 = 1, ncell
          IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .        tri(i1)%index(1).EQ.cell(i2)%index(1).AND.
     .        tri(i1)%index(2).EQ.cell(i2)%index(2)) THEN     

            tri(i1)%plasma(1:6) = cell(i2)%plasma(1:6) 
            tri(i1)%bfield(1:3) = cell(i2)%bfield(1:3) 
            tri(i1)%efield(1:3) = cell(i2)%efield(1:3) 

          ENDIF
        ENDDO
      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: NewTriangle
c
      SUBROUTINE NewTriangle(icell)
      USE mod_eirene04
      IMPLICIT none

      INTEGER icell

      ntri = ntri + 1

      tri(ntri)%type = MAGNETIC_GRID 
      tri(ntri)%index(1) = cell(icell)%index(1)
      tri(ntri)%index(2) = cell(icell)%index(2)

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: AssignVertex
c
      SUBROUTINE AssignVertex(ivert,ipoint,iside,ilist,xlist,ylist)
      USE mod_eirene04
      IMPLICIT none

      INTEGER ivert,iside,ipoint,ilist(50,2)
      REAL*8  xlist(0:50,2),ylist(0:50,2)

      IF (ilist(ipoint,iside).EQ.0) THEN
        nver = nver + 1
        ver(nver,1) = xlist(ipoint,iside)
        ver(nver,2) = ylist(ipoint,iside)
        tri(ntri)%ver(ivert) = nver
      ELSE
        tri(ntri)%ver(ivert) = ilist(ipoint,iside)
      ENDIF

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildFluidGridTriangles
c
      SUBROUTINE BuildFluidGridTriangles
      USE mod_eirene04
      IMPLICIT none

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

c      REAL    GetMach,GetJsat,GetFlux 
      LOGICAL PointOnLine

      INTEGER i1,i2,i3,i4,nlist(2),ilist(50,2)
      LOGICAL output,status
      REAL*8  x(3),y(3),s,t,xlist(0:50,2),ylist(0:50,2),slist(0:50,2)

c      REAL*8, ALLOCATABLE :: xvertex(:),yvertex(:)

      WRITE(0,*) 'BUILDING SUPER TRIANGLES'

      CALL OutputData(85,'Building super triangles')



c...  Process cells and build list of triangles and verticies:
      ntri = 0
      nver = 0

      DO i1 = 1, ncell

c        IF (.NOT.(cell(i1)%index(1).EQ.111.AND.
c     .            cell(i1)%index(2).EQ.18 .OR.
c     .            cell(i1)%index(1).EQ.112.AND.
c     .            cell(i1)%index(2).EQ.19 .OR.
c     .            cell(i1)%index(2).EQ.14 .OR.
c     .            cell(i1)%index(1).EQ.111.AND.
c     .            cell(i1)%index(2).EQ.19 )) CYCLE

        nlist = 0
        ilist = 0

        DO i2 = 1, 2
c...      Build list of points on the 14 and 23 surfaces:
          IF (i2.EQ.1) THEN
            x(1) = cell(i1)%r(1)
            y(1) = cell(i1)%z(1)
            x(2) = cell(i1)%r(4) 
            y(2) = cell(i1)%z(4)
          ELSE
            x(1) = cell(i1)%r(2)
            y(1) = cell(i1)%z(2)
            x(2) = cell(i1)%r(3) 
            y(2) = cell(i1)%z(3)
          ENDIF
c... 
          nlist(i2) = nlist(i2) + 1
          xlist(nlist(i2),i2) = x(1)
          ylist(nlist(i2),i2) = y(1)
          slist(nlist(i2),i2) = 0.0D0
c...      Search all other cells for points that are on the 14 and 23 sides of the 
c         current focus cell:
          DO i3 = 1, ncell
            IF (i1.EQ.i3) CYCLE   ! Speed this up, perhaps with a quick distance check, or apply a zone above...

            output = .FALSE.
c            IF (i2.EQ.2.AND.
c     .          cell(i1)%index(1).EQ.1 .AND.
c     .          cell(i1)%index(2).EQ.62.AND.
c     .          (cell(i3)%index(2).EQ.63.OR.
c     .           .FALSE.)) output = .TRUE.
c            IF (i2.EQ.2.AND.cell(i3)%index(2).EQ.19) output = .TRUE.


            IF (cell(i3)%index(1).EQ.1) THEN  ! *THIS MAY NOT WORK FOR EDGE2D*
c...          Only check points 1 and 2 on the 1st cells of a given ring.  In theory, 
c             only points on side 23 of neighbouring cells need to be checked:
              IF (i2.EQ.1) THEN
                x(3) = cell(i3)%r(2)
                y(3) = cell(i3)%z(2)
              ELSE
                x(3) = cell(i3)%r(1)
                y(3) = cell(i3)%z(1)
              ENDIF
              IF (PointOnLine(x,y,s,t,2,output)) THEN
c                WRITE(0,*) 'A:',i2,i3
                nlist(i2) = nlist(i2) + 1
                xlist(nlist(i2),i2) = x(3)
                ylist(nlist(i2),i2) = y(3)
                slist(nlist(i2),i2) = s
              ENDIF
            ENDIF
c...        Check ...:
            IF (i2.EQ.1) THEN
              x(3) = cell(i3)%r(3)
              y(3) = cell(i3)%z(3)
            ELSE
              x(3) = cell(i3)%r(4)
              y(3) = cell(i3)%z(4)
            ENDIF
            IF (PointOnLine(x,y,s,t,2,output)) THEN
c              WRITE(0,*) 'B:',i2,i3
              nlist(i2) = nlist(i2) + 1
              xlist(nlist(i2),i2) = x(3)
              ylist(nlist(i2),i2) = y(3)
              slist(nlist(i2),i2) = s
            ENDIF
          ENDDO
c...
          nlist(i2) = nlist(i2) + 1
          xlist(nlist(i2),i2) = x(2)
          ylist(nlist(i2),i2) = y(2)
          slist(nlist(i2),i2) = 1.0D0

c          DO i3 = 1, nlist(i2)
c            WRITE(0,'(A,I3,3F10.4)') 
c     .        'LIST A:',i2,xlist(i3,i2),ylist(i3,i2),slist(i3,i2)
c          ENDDO

        ENDDO


c...    Process list of points:
        DO i2 = 1, 2
c...      Eliminate duplicates:
          DO i3 = 1, nlist(i2)-1
            DO i4 = i3+1, nlist(i2)
              IF (slist(i3,i2).EQ.slist(i4,i2)) slist(i4,i2) = -999.0D0
            ENDDO
          ENDDO
          DO i3 = nlist(i2), 1, -1
            IF (slist(i3,i2).EQ.-999.0D0) THEN 
c              WRITE(0,*) 'ELIMINATING TRIANGLE POINT!', ! This should be unecessary...
c     .          cell(i1)%index(1),cell(i1)%index(2)
              DO i4 = i3, nlist(i2)-1
                xlist(i4,i2) = xlist(i4+1,i2)
                ylist(i4,i2) = ylist(i4+1,i2)
                slist(i4,i2) = slist(i4+1,i2)
              ENDDO
              nlist(i2) = nlist(i2) - 1
            ENDIF
          ENDDO
c...      Sort:
          i3 = 1
          DO WHILE (i3.LT.nlist(i2))
            status = .FALSE.
            DO i4 = i3+1, nlist(i2)
              IF (slist(i4,i2).LT.slist(i3,i2)) THEN
                WRITE(0,*) 'SORTING TRIANGLE SIDE!',    ! This sorting should be unnecessary...
     .            cell(i1)%index(1),cell(i1)%index(2)
                status = .TRUE.
                xlist(0 ,i2) = xlist(i3,i2)
                ylist(0 ,i2) = ylist(i3,i2)
                slist(0 ,i2) = slist(i3,i2)
                xlist(i3,i2) = xlist(i4,i2)
                ylist(i3,i2) = ylist(i4,i2)
                slist(i3,i2) = slist(i4,i2)
                xlist(i4,i2) = xlist(0 ,i2)
                ylist(i4,i2) = ylist(0 ,i2)
                slist(i4,i2) = slist(0 ,i2)
              ENDIF
            ENDDO
            IF (.NOT.status) i3 = i3 + 1
          ENDDO

c          DO i3 = 1, nlist(i2)
c            WRITE(0,'(A,I3,3F10.4)') 
c     .        'LIST B:',i2,xlist(i3,i2),ylist(i3,i2),slist(i3,i2)
c          ENDDO

        ENDDO

c...    Check if vertices in list are already in the vertex list:
        DO i2 = 1, 2
          DO i3 = 1, nver
            DO i4 = 1, nlist(i2)
              IF (ABS(ver(i3,1)-xlist(i4,i2)).LT.TOL.AND.
     .            ABS(ver(i3,2)-ylist(i4,i2)).LT.TOL) ilist(i4,i2) = i3
            ENDDO
          ENDDO
        ENDDO           


c...    Build triangles:

        DO i2 = 1, MAX(nlist(1),nlist(2))-1

          IF (nlist(1).GT.i2.AND.nlist(2).GT.i2) THEN

            CALL NewTriangle(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,4)
            IF (i2.EQ.1) 
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,1)

            CALL AssignVertex(1,i2  ,2,ilist,xlist,ylist)
            CALL AssignVertex(2,i2  ,1,ilist,xlist,ylist)
            CALL AssignVertex(3,i2+1,1,ilist,xlist,ylist)

            CALL NewTriangle(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            IF (i2.EQ.nlist(1)-1.AND.nlist(1).EQ.nlist(2)) 
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)

            CALL AssignVertex(1,i2+1,1,ilist,xlist,ylist)
            CALL AssignVertex(2,i2+1,2,ilist,xlist,ylist)
            CALL AssignVertex(3,i2  ,2,ilist,xlist,ylist)

          ELSEIF (nlist(1).GT.i2) THEN

            CALL NewTriangle(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,4)
            IF (i2.EQ.nlist(1)-1) 
     .        tri(ntri)%sideindex(2,3) = cell(i1)%sideindex(2,3)

            CALL AssignVertex(1,nlist(2),2,ilist,xlist,ylist)
            CALL AssignVertex(2,i2      ,1,ilist,xlist,ylist)
            CALL AssignVertex(3,i2+1    ,1,ilist,xlist,ylist)

          ELSEIF (nlist(2).GT.i2) THEN

c            WRITE(0,*) 'COOL B:',i2  

            CALL NewTriangle(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            IF (i2.EQ.nlist(2)-1)
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)

            CALL AssignVertex(1,nlist(1),1,ilist,xlist,ylist)
            CALL AssignVertex(2,i2+1    ,2,ilist,xlist,ylist)
            CALL AssignVertex(3,i2      ,2,ilist,xlist,ylist)

          ELSE
            CALL ER('BuildFluidGridTriangles','Unknown situation',*99)
          ENDIF

        ENDDO

      ENDDO


c...  Remove vertex duplicates:




c      CALL WriteTriangleFiles
c      CALL DumpGrid('PROBLEM WITH TRIANGLE MAP')
c      STOP 'CRAPPO!'





c      STOP 'DUMPING'
      WRITE(0,*) 'DONE'

      RETURN
96    WRITE(0,*) 'BUILDSUPERTRIANGES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END
c
c
c ======================================================================
c
c  subroutine: SaveTriangles
c
      subroutine SaveTriangles
      USE mod_eirene04
      IMPLICIT none

      INTEGER fp,i1,i2

      WRITE(0,*) 'ERROR: Old triangle .raw file not supported'
      STOP 

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.0,ntri,nver,nadd
      WRITE(fp,ERR=98) (tri(i1),i1=1,ntri)
      WRITE(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      WRITE(fp,ERR=98) (add(i1),i1=1,nadd)
      CLOSE (fp)
      
      RETURN
 98   CALL ER('SaveTriangles','Problems writing data file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteEireneInputFile_04
c
      SUBROUTINE WriteEireneInputFile_04
      USE mod_eirene04
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'slcom'

      INTEGER    MAXSTRAT,MAXSDATA
      PARAMETER (MAXSTRAT=10,MAXSDATA=10)
      COMMON /EIRCOM/
     .        eir_07nsrc,eir_07opt ,eir_07stra,
     .        eir_07ind1,
     .        eir_07ind2,eir_07ind3,
     .        eir_07wght           ,eir_07data
      INTEGER eir_07nsrc,eir_07opt ,eir_07stra(MAXSTRAT),
     .        eir_07ind1(MAXSTRAT),
     .        eir_07ind2(MAXSTRAT),eir_07ind3(MAXSTRAT)
      REAL    eir_07wght           ,eir_07data(MAXSTRAT,MAXSDATA)

      INTEGER   ik,ik1,ik2,ir,i1,i2,i3,fp1,fp2,in,icnt,
     .          add1,ilst(1024)
      LOGICAL   output,firstcall
      REAL      x0,y0,r,zaa,roa,fact
c      REAL      x0,y0,r,vcel(MAXASCDAT),zaa,roa,fact
      CHARACTER buffer*200,geostr*4

      DATA firstcall /.TRUE./
      SAVE
c
c     Check whether DIVIMP input option requests EIRENE data file:
c
c      IF (eirdata.NE.1) RETURN

      output = .FALSE.

      IF (output) WRITE(0,*) 'WRITING EIRENE INPUT FILE'
c
c     Initialization:
c
      fp1   = 97
      fp2   = 98
c      fp1   = 80
c      fp2   = 81
c      fp1   = EIRIN
c      fp2   = EIROUT

c...  eirene.input is assumed to be a template input file ...:
      OPEN(UNIT=fp1,FILE='eirene.template',FORM='FORMATTED',
     .     ERR=95,STATUS='OLD')
      OPEN(UNIT=fp2,FILE='eirene.input',FORM='FORMATTED',
     .     ERR=96,STATUS='REPLACE')
c      OPEN(UNIT=fp1,FORM='FORMATTED',ERR=95,STATUS='OLD')
c      OPEN(UNIT=fp2,FORM='FORMATTED',ERR=95,STATUS='REPLACE')

c      fp2 = 0
c      REWIND(fp1)

      CALL MS('WriteInputFile','Using xVESM to store wall data')

cc      IF (iflexopt(6).EQ.11) THEN
c        eirtemp1 = -ctargt * 1.38E-23 / ECH
c        eirtemp2 = -cwallt * 1.38E-23 / ECH
cc      ELSE
cc        eirtemp1 = ctargt * 1.38E-23 / ECH
cc        eirtemp2 = cwallt * 1.38E-23 / ECH
cc      ENDIF


c...  Correct simulated pressure gauge volumes for cyclindrical approximation:

c
c     This is a somewhat convoluted loop at the moment, which reads
c     through an existing EIRENE input data file that serves as a
c     template for the new data file being written.  Grid and neutral
c     wall specific data are substituted into the template, as well
c     as any EIRENE options/settings that are specifiable
c     from DIVIMP (such as EIRENE execution time):
c




10    CONTINUE

      CALL ReadLine(fp1,buffer,1,*50,*98)

20    CONTINUE

      IF (output) WRITE(0,*) '  BUFFER:'//buffer(1:6)

      IF (buffer(1:6).EQ.'*** 0.') THEN
c
c Need to remove the requirement that the template file have an intitial
c seciton labelled *** 0...
c
c       This section has been added to EIRENE and contains options
c       for the new EIRENE code that are related to the
c       generalization of the grid:
c
        WRITE(0,*) 'DEFUNCT - HALTING CODE'
        STOP

c        WRITE(fp2,'(A)') '*** 0. DIVIMP RELATED SETUP DATA (DIVIMP)'
c        WRITE(fp2,'(2A,I6)') '''Geometry option  (GEOMOPT)  ',
c     .    '0-standard   1-from DIVIMP    ''',eirgeom
c        WRITE(fp2,'(2A,I6)') '''Grid option      (GRIDOPT)  ',
c     .    '0-structured 1-generalized    ''',eirgrid
cc        WRITE(fp2,'(2A,I6)') '''AddUsr option    (ADDOPT)   ',
cc     .    '0-execute    1-do not execute ''',eiradd
c        WRITE(fp2,'(2A,I6)') '''Wall data option (NEUTOPT)  ',
c     .    '0-standard   1-accurate       ''',eirneut
c        WRITE(fp2,'(2A,I6)') '''Debug option     (DEBUGOPT) ',
c     .    '0-off                         ''',eirdebug
c        WRITE(fp2,'(2A,I6)') '''CX D2+ production(CXD2OPT)  ',
c     .    '0-off 1-Dalpha only 2-full    ''',eircxd2
c        WRITE(fp2,'(2A,I6)') '''User ID          (OPTUSER)  ',
c     .    '                              ''',optuser
c        WRITE(fp2,'( A,I6)') '''Transparent toroidal spans of non-'//
c     .                      'standard surfaces       ''',eirntrans
c        IF (eirntorseg.NE.0) THEN
c          fact = -360.0 / eirzaa
c        ELSE
c          fact = 100.0
c        ENDIF
c        DO i1 = 1, eirntrans
c          WRITE(fp2,'(F6.1,2F12.4)') 
c     .      eirtrans(i1,1),(eirtrans(i1,i2)*fact,i2=2,3)
c        ENDDO

c        icnt = 0
c        DO i1 = 1, eirnaout
c          IF (eiraout(i1,1).EQ.1.0) THEN
c            add1 = MAX(0,(eiraout(i1,2) - 1) * asc_ncell) 
c            DO i2 = 3,6
c              IF (eiraout(i1,i2).GT.0.0) THEN
c                icnt = icnt + 1
c                ilst(icnt) = eiraout(i1,i2) + add1 + 1 + eirnpgdat
c              ENDIF
c            ENDDO
c          ENDIF
c        ENDDO
c        WRITE(fp2,'( A,I6)') '''Additional cells for recodring '//
c     .                       'stratum data''',icnt
c        WRITE(fp2,'(6I6)') (ilst(i1),i1=1,icnt)        


c        IF (eirzaa.LT.0.0) THEN
c          WRITE(fp2,'(A,F9.3)') '''Global volume scaling''',eirtorfrac
c        ELSE
c          WRITE(fp2,'(A,F9.3)') '''Global volume scaling''',1.0
c        ENDIF

c        IF (grdnmod.NE.0) THEN
c          WRITE(fp2,'(A,I9)') '''Radial connection map''',1
c        ELSE
c          WRITE(fp2,'(A,I9)') '''Radial connection map''',0
c        ENDIF


c
c       Advance the template file until the next section tag
c       is found:
c
c21      CALL ReadLine(fp1,buffer,1,*97,*98)
c        IF (buffer(1:3).NE.'***') GOTO 21
c        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 1.') THEN

        fp04 = fp2

c        IF (eirphoton.EQ.2) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,3)
c        ELSE
          CALL WriteBlock01_04(fp1,fp2)
22        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 22
          BACKSPACE fp1
c        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 2.') THEN

        CALL WriteBlock02_04

24      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 24
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 3a') THEN

        CALL WriteBlock03a_04

26      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 26
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 3b'.OR.
     .        buffer(1:6).EQ.'*** 3B') THEN

        CALL WriteBlock03b_04

25      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 25
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 4.') THEN

c        IF (eirphoton.GT.0) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,98)
c        ELSE
          CALL WriteBlock04_04

41        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 41
          BACKSPACE fp1
c        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 5.') THEN

c        IF     (eirphoton.EQ.1) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,116)
c        ELSEIF (eirphoton.EQ.2) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,118)
c        ELSE
          CALL WriteBlock05_04

39        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 39
          BACKSPACE fp1
c        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 6.') THEN

        IF (photons.GT.0) THEN
c        IF (eirphoton.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,14)
        ELSE
          CALL WriteBlock06_04(fp1,fp2)
        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 7.') THEN

        IF (photons.GT.0) THEN
c        IF (eirphoton.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,186)
        ELSE
          CALL WriteLine(fp2,buffer)
          DO WHILE (.TRUE.) 
            CALL ReadLine(fp1,buffer,1,*97,*98)
            IF (buffer(1:3).NE.'***') THEN
              CALL WriteLine(fp2,buffer)
            ELSE
              BACKSPACE fp1
              EXIT
            ENDIF
          ENDDO
        ENDIF

      ELSEIF (buffer(1:7).EQ.'*** 10.') THEN

        IF (photons.GT.0) THEN
c        IF (eirphoton.GT.0) THEN
          WRITE(0,*) 'HERE IN PHOTON CODE!'
          CALL WriteLine(fp2,buffer)
          DO WHILE (.TRUE.) 
            CALL ReadLine(fp1,buffer,1,*97,*98)
            IF (buffer(1:3).NE.'***') THEN
              CALL WriteLine(fp2,buffer)
            ELSE
              BACKSPACE fp1
              EXIT
            ENDIF
          ENDDO
c          CALL Transferline2(fp1,fp2,buffer,40)
        ELSE
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,32)
        ENDIF

      ELSEIF (.FALSE..AND.buffer(1:7).EQ.'*** 10.') THEN
c OFF
c        WRITE(fp2,35) '*** 10. DATA FOR ADDITIONAL TALLIES - DIVIMP'
c        WRITE(fp2,36) 0,0,eirntally,0,0,0
c        WRITE(fp2,35) '*** 10A.'
c        WRITE(fp2,35) '*** 10B.'
c        WRITE(fp2,35) '*** 10C.'
c        DO i1 = 1, eirntally
c          WRITE(fp2,37) eirtally(i1,1)(1:LEN_TRIM(eirtally(i1,1)))
c          WRITE(fp2,37) eirtally(i1,2)(1:69)//eirtally(i1,5)(1:3)
c          WRITE(fp2,37) eirtally(i1,3)(1:24),
c     .                  eirtally(i1,4)(1:24)
c        ENDDO
c        WRITE(fp2,35) '*** 10D.'
c        WRITE(fp2,35) '*** 10E.'        

c35      FORMAT(A)
c36      FORMAT(6(I6:))
c37      FORMAT(6(A:))

c...    Read to end of input block in template file:
c40      CALL ReadLine(fp1,buffer,1,*50,*98)
c        IF (buffer(1:6).NE.'*** 11') GOTO 40
c        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 11') THEN

        CALL WriteBlock11_04

29      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 29
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 12') THEN

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,1)

      ELSEIF (buffer(1:6).EQ.'*** 13') THEN

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,1)

c        IF (eirdtimv.NE.0.0) THEN
c          WRITE(fp2,'(A)') '*** 13. DATA FOR NONLI. AND/OR TIME DEP.'
c          WRITE(fp2,'( I6)') 199999
c          WRITE(fp2,'(2I6)') 0,1
c          WRITE(fp2,'(1P,2E12.4)') eirdtimv,0.0
c          WRITE(fp2,'(A)') '** 13A. DATA FOR SNAPSHOT TALLIES'
c          WRITE(fp2,'( I6)') 0
c          CALL ReadLine(fp1,buffer,1,*50,*98)
c        ELSE
c          CALL WriteLine(fp2,buffer)
c          CALL TransferLine(fp1,fp2,buffer,1)
c        ENDIF

      ELSEIF (buffer(1:6).EQ.'*** 14') THEN

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,16)
        GOTO 50

      ELSE
c
c       Input block identifier is not recognized, so
c       copy the entire section as is:
c
        CALL WriteLine(fp2,buffer)

        GOTO 10
      ENDIF

      CALL ReadLine(fp1,buffer,1,*97,*98)

      IF (buffer(1:3).NE.'***') THEN
        WRITE(0,*) 'BAD BUFFER:'//buffer(1:3)
        CALL ER('WriteInputFile','Invalid template format',*99)
      ENDIF

      GOTO 20

50    CONTINUE

      CLOSE (fp1)
      CLOSE (fp2)

c      IF (.FALSE..AND.eirneut.EQ.0) THEN
c        nvesm = 0
c        write(0,*)
c        write(0,*) ' TEMPORARY BLANKING OF NEUTRAL WALL DATA! '
c        write(0,*)
c      ENDIF

c      STOP 'WRITE INPUT FILE'
      RETURN
c
c     Error code:
c 
95    WRITE(0,*) 'FILE ERROR A'
96    WRITE(0,*) 'FILE ERROR B'
      STOP
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(50,*) '  Last line read: '
      WRITE(50,*) '  "',buffer,'"'
c99    WRITE(EROUT,*) '  Last line read: '
c      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock01
c
c
c
c
      SUBROUTINE WriteBlock01_04
      USE mod_eirene04
      IMPLICIT none

      INTEGER   ntime

      WRITE(fp04,90) '*** 1. DATA FOR OPERATING MODE (DIVIMP)'

      ntime = 0
      IF (dtimv.NE.0.0) ntime = 1

      nfile = 111
      IF (photons.EQ.2) nfile = 311

      IF     (niter.GE.1) THEN
c...    BGK or photons:
        WRITE(fp04,91) 2,0,time,nfile,0,niter,0,ntime
        WRITE(fp04,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp04,90) 'FFFFF FFFF'
      ELSEIF (.TRUE.) THEN
c...    Standard (no BGK or photons):
        WRITE(fp04,91) 2,0,time,nfile,0,0,0,ntime
c        WRITE(fp04,91) 2,0,time,nfile,0,1,0,ntime
c        WRITE(fp04,91) 1,1,0,0,1,9,1,0,5  ! NGSTAL=1
        WRITE(fp04,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp04,90) 'FFFFF FFFF'
      ELSE
        CALL ER('WriteBlock01_04','Trouble',*99)
      ENDIF


      RETURN
90    FORMAT(A)
91    FORMAT(3I6,I6.5,20(I6:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock02_04
c
c
c
c
      SUBROUTINE WriteBlock02_04
      USE mod_eirene04
      IMPLICIT none

      WRITE(fp04,90) '*** 2. DATA FOR STANDARD MESH (DIVIMP)'

      IF (.TRUE.) THEN
        WRITE(fp04,91) 1,1,1,1
        WRITE(fp04,90) 'T'
        WRITE(fp04,90) 'FFFFF TFF'
c        WRITE(fp04,91) ntri,0,0,0,nver
        WRITE(fp04,91) ntri+1,0,0,0,nver
        WRITE(fp04,90) 'CASE cmod'
        WRITE(fp04,90) 'F'
        WRITE(fp04,90) 'TFFF'
        WRITE(fp04,91) 1,0
        WRITE(fp04,92) 0.0,0.0,0.0
        WRITE(fp04,90) 'F'
        WRITE(fp04,90) 'FTFF'
        WRITE(fp04,91) 1,1,100
        WRITE(fp04,92) 0.0,0.0,360.0
        WRITE(fp04,90) 'F'
        WRITE(fp04,91) 0
c        WRITE(fp04,90) 'T'
c        WRITE(fp04,91) 1
c        WRITE(fp04,93) 1.0E+00
        WRITE(fp04,90) 'F'
        WRITE(fp04,91) 0
      ELSE
        CALL ER('WriteBlock02_04','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock04_04
c
c
c
c
      SUBROUTINE WriteBlock04_04
      USE mod_eirene04
      IMPLICIT none

      INTEGER ncra,ncrm,iscd1,iscd2

      WRITE(fp04,90) '*** 4. DATA FOR SPECIES SPECIFICATION AND '//
     .               'ATOMIC PHYSICS MODULE (OSM)'

      IF     (.TRUE.) THEN

        WRITE(fp04,90) '* ATOMIC REACTION CARDS  NREACI='
        WRITE(fp04,91) 33
        WRITE(fp04,95) '  1 AMJUEL H.4 2.1.5    EI ',0,1
        WRITE(fp04,95) '  2 CONST  H.4          EI ',0,1     ! Dummy ionisation rate
        WRITE(fp04,92) -200.0,0.0,0.0,0.0,0.0,0.0
        WRITE(fp04,92)    0.0,0.0,0.0					  
        WRITE(fp04,95) '  3 AMJUEL H.102.1.5    EI ',0, 1
        WRITE(fp04,95) '  4 HYDHEL H.1 3.1.8    CX ',1, 1
        WRITE(fp04,95) '  5 HYDHEL H.3 3.1.8    CX ',1, 1
        WRITE(fp04,95) '  6 AMJUEL H.2 2.26B0   EI ',0,56
        WRITE(fp04,95) '  6 METHAN H.2 2.23      EI',0,12
        WRITE(fp04,95) '  7 METHAN H.1 3.2       CX',1,12
        WRITE(fp04,95) '  7 METHAN H.3 3.2       CX',1,12
        IF (opacity.EQ.5.AND.photons.EQ.0) THEN
          WRITE(fp04,95) '  8 H2VIBR H.4 2.1.8a   RC ',0,1
        ELSE
          WRITE(fp04,95) '  8 AMJUEL H.4 2.1.8    RC ',0,1				  
        ENDIF
        WRITE(fp04,95) '  9 HYDHEL H.2 2.2.9    EI ',0, 2,0.0,0.0,0.0
        WRITE(fp04,95) ' 10 HYDHEL H.2 2.2.5    DS ',0, 2,0.0,0.0,0.0
        WRITE(fp04,95) ' 11 HYDHEL H.2 2.2.10   DS ',0, 2,0.0,0.0,0.0
        WRITE(fp04,95) ' 13 AMJUEL H.0 0.3T     EL ',1, 2				  
        WRITE(fp04,95) ' 13 AMJUEL H.1 0.3T     EL ',1, 2				  
        WRITE(fp04,95) ' 13 AMJUEL H.3 0.3T     EL ',1, 2,0.0,0.0,0.0
        WRITE(fp04,95) ' 14 AMJUEL H.4 2.2.12   EI ',0, 2
        WRITE(fp04,95) ' 15 AMJUEL H.4 2.2.11   EI ',0, 2
        WRITE(fp04,95) ' 16 AMJUEL H.4 2.2.14   EI ',0, 2
        WRITE(fp04,95) ' 17 AMJUEL H.8 2.2.14   EI ',0, 2
        WRITE(fp04,95) ' 18 AMJUEL H.3 3.2.3    CX ',1, 2
        WRITE(fp04,96) ' 19 PHOTON H.2 HLya121.5669a RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 20 PHOTON H.2 HLyb102.5722a RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 21 PHOTON H.2 HBaa656.4667a RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 22 PHOTON P.1 HLya121.5669a OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 23 PHOTON P.1 HLyb102.5722a OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 24 PHOTON P.1 HBaa656.4667a OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 25 PHOTON H.2 HLyg97.2536a  RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 26 PHOTON P.1 HLyg97.2536a  OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 27 PHOTON H.2 HLyd94.9742a  RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 28 PHOTON P.1 HLyd94.9742a  OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 29 PHOTON H.2 HLye93.7803a  RC ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0
        WRITE(fp04,96) ' 30 PHOTON P.1 HLye93.7803a  OT ',0,1,0.,0.,0.
        WRITE(fp04,91) 3,2,1,0

        WRITE(fp04,90) ' 31 CONST  H.2           EL  2  2'				! 17 -> 31
        WRITE(fp04,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0 
        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0 
        WRITE(fp04,90) ' 32 CONST  H.2           EL  2  2'	                        ! 18 -> 32
        WRITE(fp04,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
        WRITE(fp04,90) ' 33 CONST  H.2           EL  4  4'                              ! 19 -> 33
        WRITE(fp04,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0

        WRITE(fp04,90) '** 4a NEUTRAL ATOMS SPECIES CARDS: NATMI='
        WRITE(fp04,91) 1							    
        ncra = 3
        ncrm = 5
        IF (bgk.EQ.3) THEN
          ncra = ncra + 2
          ncrm = ncrm + 2
        ENDIF

c        WRITE(fp04,94) 1,'D(n=1)  ',2,1,1,0,0,-4,0,ncra ! Causes atoms to be pumped, but not officially...
        WRITE(fp04,94) 1,'D(n=1)  ',2,1,1,0,1,-4,0,ncra			    
        WRITE(fp04,91) 1,115,114,0  ,30000
        WRITE(fp04,93) 2.0,0.0,0.0,0.0,1.0
        IF (photons.EQ.-1) THEN
          WRITE(fp04,91) 2,115,114,0  ,30000
          WRITE(fp04,93) 999.0,0.0,0.0,0.0,1.0
        ELSE
          WRITE(fp04,91) 2,115,114,0  ,30000
          WRITE(fp04,93) 2.0,0.0,0.0,0.0,1.0
        ENDIF
        WRITE(fp04,91) 4,114,111,114,01001				    
        WRITE(fp04,93) 0.0,0.0,0.0,0.0,1.0
        IF (bgk.EQ.3) THEN
          WRITE(fp04,91) 31,214,0,0,01001,0,111
          WRITE(fp04,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp04,91) 33,414,0,0,01001,0,112				    
          WRITE(fp04,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF

        WRITE(fp04,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
        WRITE(fp04,91) 1							    
        WRITE(fp04,94) 1,'D2      ',4,2,2,0,1,1,0,ncrm,0,0
        WRITE(fp04,91)  9,115,113,0,0			    
        WRITE(fp04,93) -1.5400E+01,0.0,0.0,0.0,0.0
        WRITE(fp04,91) 10,115,121,0,0
        WRITE(fp04,93) -1.0500E+01,0.0,3.0,3.0,0.0
        WRITE(fp04,91) 11,115,111,114,0
        WRITE(fp04,93) -2.5000E+01,0.0,5.0,5.0,0.0
        WRITE(fp04,91) 13,114,  0,  0,01001				    
c        WRITE(fp04,92) 0.0,0.0,0.0,0.0,1.0E-10
        WRITE(fp04,92) 0.0,0.0,0.0,0.0
        WRITE(fp04,91) 18,114,111,113,01001				    
        WRITE(fp04,92) 0.0,0.0,0.0,0.0
        IF (bgk.EQ.3) THEN
          WRITE(fp04,91) 32,314,0,0,01001,0,112
          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp04,91) 33,514,0,0,01001,0,111
          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF

        WRITE(fp04,90) '**4c TEST ION SPECIES CARDS:  NIONI ION '//
     .                 'SPECIES ARE CONSIDERED, NIONI='
        WRITE(fp04,91) 1						
        WRITE(fp04,94) 1,'D2+     ',4,2,2,1,0,-1,0,3,-1		
        WRITE(fp04,91) 14,115,111,114,0
        WRITE(fp04,92) -1.0400E+01,0.0,4.3,4.3,0.0
        WRITE(fp04,91) 15,115,124,0,0
        WRITE(fp04,92) -1.5500E+01,0.0,0.25,0.25,0.0
        WRITE(fp04,91) 16,115,121,000,30000
        WRITE(fp04,92) 16.0,0.0,10.0,0.0,0.0

        IF (photons.EQ.1.OR.photons.EQ.2) THEN

          iscd1 = 214
          iscd2 = 0
          IF (bgk.EQ.3) THEN
            iscd1 = 614
            iscd2 = 4
          ENDIF

          WRITE(fp04,90) '** 4d photons'
          WRITE(fp04,91) 6			
          WRITE(fp04,94) 1,'Ba-alpha',2,1,0,0,1,+1,0,0
          WRITE(fp04,94) 2,'Ly-alpha',2,1,0,0,1,+1,0,1
          WRITE(fp04,91) 22,iscd1,( 3+iscd2)*100+14,0,0
c          WRITE(fp04,91) 22,iscd1,314,0,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0
          WRITE(fp04,94) 3,'Ly-beta ',2,1,0,0,1,+1,0,1
          WRITE(fp04,91) 23,iscd1,( 5+iscd2)*100+14,0,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0
          WRITE(fp04,94) 4,'Ly-gamma',2,1,0,0,1,+1,0,1
          WRITE(fp04,91) 26,iscd1,( 7+iscd2)*100+14,0,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0
          WRITE(fp04,94) 5,'Ly-delta',2,1,0,0,1,+1,0,1
          WRITE(fp04,91) 28,iscd1,( 9+iscd2)*100+14,0,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0
          WRITE(fp04,94) 6,'Ly-epsi ',2,1,0,0,1,+1,0,1
          WRITE(fp04,91) 30,iscd1,(11+iscd2)*100+14,0,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0
        ELSE
          WRITE(fp04,90) '** 4d photons'
          WRITE(fp04,91) 0
        ENDIF

      ELSEIF (.FALSE.) THEN

c        WRITE(fp04,90) '* ATOMIC REACTION CARDS  NREACI='
c        WRITE(fp04,91) 22
c        WRITE(fp04,90) '  1 AMJUEL H.4 2.1.5     EI  0  1'
c        WRITE(fp04,90) '  2 AMJUEL H.102.1.5     EI  0  1'
c        WRITE(fp04,90) '  3 HYDHEL H.1 3.1.8     CX  1  1'
c        WRITE(fp04,90) '  3 HYDHEL H.3 3.1.8     CX  1  1'
c        WRITE(fp04,90) '  4 AMJUEL H.4 2.2.9     EI  0  2'
c        WRITE(fp04,90) '  5 AMJUEL H.4 2.2.5     DS  0  2'
c        WRITE(fp04,90) '  6 AMJUEL H.4 2.2.10    DS  0  2'
c        WRITE(fp04,90) '  7 AMJUEL H.4 2.2.12    DS  0  2'
c        WRITE(fp04,90) '  8 AMJUEL H.4 2.2.11    DS  0  2'
c        WRITE(fp04,90) '  9 AMJUEL H.4 2.2.14    DS  0  2'
c        WRITE(fp04,90) ' 10 AMJUEL H.8 2.2.14    DS  0  2'
c        WRITE(fp04,90) ' 12 AMJUEL H.0 0.3T     EL   1  2'
c        WRITE(fp04,90) ' 12 AMJUEL H.1 0.3T     EL   1  2'
c        WRITE(fp04,90) ' 12 AMJUEL H.3 0.3T     EL   1  2 0.00000E+00'//
c     .                 ' 0.00000E+00 0.00000E+00'
c        WRITE(fp04,90) ' 13 HYDHEL H.2 2.3.9     EI  0  4'
c        WRITE(fp04,90) ' 14 METHAN H.2 2.23      EI  0 12'
c        IF (.TRUE.) THEN
c          WRITE(fp04,90) ' 15 H2VIBR H.4 2.1.8a    RC  0  1'
c          WRITE(fp04,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ELSE
c          WRITE(fp04,90) ' 15 AMJUEL H.4 2.1.8     RC  0  1'
c          WRITE(fp04,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ENDIF
c        WRITE(fp04,90) ' 17 CONST  H.2           EL  2  2'				
c        WRITE(fp04,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp04,90) ' 18 CONST  H.2           EL  2  2'				
c        WRITE(fp04,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp04,90) ' 19 CONST  H.2           EL  4  4'
c        WRITE(fp04,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp04,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp04,90) ' 20 AMJUEL H.3 3.2.3     CX  1  2'
c        WRITE(fp04,90) ' 21 AMJUEL H.9 3.1.8     CX  1  1'
c        WRITE(fp04,90) ' 22 AMJUEL H.2 3.1.8FJ   CX  1  1'

c        WRITE(fp04,90) '*NEUTRAL ATOMS SPECIES CARDS: NATMI='
c        WRITE(fp04,91) 1							    
c        ncra = 2
c        ncrm = 5
c        IF (niter.GE.1) THEN
c          ncra = ncra + 2
c          ncrm = ncrm + 2
c        ENDIF

c        WRITE(fp04,94) 1,'D       ',2,1,1,0,1,-4,0,ncra			    
c        WRITE(fp04,91) 1,115,114,0  ,30000,000		    
c        WRITE(fp04,92) 2.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp04,91) 3,114,111,114,01001				    
c        WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        IF (niter.GE.1) THEN
c          WRITE(fp04,91) 17,214,0,0,01001,0,111
c          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp04,91) 19,414,0,0,01001,0,112				    
c          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF

c        WRITE(fp04,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
c        WRITE(fp04,91) 1							    
c        WRITE(fp04,94) 1,'D2      ',4,2,2,0,0,2,0,ncrm
c        WRITE(fp04,91) 4,115,113,0				    
c        WRITE(fp04,92) -1.5400E+01,0.0000E+00				    
c        WRITE(fp04,91) 5,115,121,000				    
c        WRITE(fp04,92) -1.0500E+01,0.0000E+00,3.0000E+00,3.0000E+00	    
c        WRITE(fp04,91) 6,115,111,114				    
c        WRITE(fp04,92) -2.5000E+01,0.0000E+00,5.0000E+00,5.0000E+00	    
c        IF (niter.GE.1) THEN
c          WRITE(fp04,91) 18,314,0,0,01001,0,112
c          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp04,91) 19,514,0,0,01001,0,111
c          WRITE(fp04,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF
c        WRITE(fp04,91) 20,114,111,113,01001				    
c        WRITE(fp04,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp04,91) 12,114,  0,  0,01001				    
c        WRITE(fp04,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c
c        WRITE(fp04,90) '** 4c TEST ION SPECIES CARDS:  NIONI ION '//
c     .                 'SPECIES ARE CONSIDERED, NIONI='
c        WRITE(fp04,91) 1						
c        WRITE(fp04,94) 1,'D2+     ',4,2,2,1,0,-4,0,3,-1		
c        WRITE(fp04,91) 7,115,111,114                        
c        WRITE(fp04,92) -1.0400E+01,0.0000E+00,4.3000E+00,4.3000E+00
c        WRITE(fp04,91) 8,115,124,000                        
c        WRITE(fp04,92) -1.5500E+01,0.0000E+00,0.2500E+00,0.2500E+00
c        WRITE(fp04,91) 9,115,121,000,30002                  
c        WRITE(fp04,92)  1.0000E+01,0.0000E+00,0.5000E+00,0.5000E+00

c        WRITE(fp04,90) '** 4d photons'
c        WRITE(fp04,91) 0

      ELSE
        CALL ER('WriteBlock04_04','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
94    FORMAT(I2,1X,A8,12(I3:))
95    FORMAT(A,2I3:,1P,3E12.5)
96    FORMAT(A,2I3:,1P,3E12.4)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock11_04
c
c
c
c
      SUBROUTINE WriteBlock11_04
      USE mod_eirene04
      IMPLICIT none

      WRITE(fp04,90) '*** 11. DATA FOR NUMERICAL/GRAPHICAL OUTPUT(OSM)'

      IF (.TRUE.) THEN
        WRITE(fp04,90) 'FTFFF TTFFF FFTFF TTTTT T'
        WRITE(fp04,90) 'Ttttt ttttt tttt'
        WRITE(fp04,91) 4
        WRITE(fp04,91) 14,0
        WRITE(fp04,91) -2,0
        WRITE(fp04,91) -3,0
        WRITE(fp04,91) -4,0
        WRITE(fp04,91) 0
        WRITE(fp04,90) 'TTFTT FFTFT ftFFF F'
        WRITE(fp04,91) 1,ntri,1,1,1,1
        WRITE(fp04,90) 'F PEI                      1 001002'
        WRITE(fp04,90) 'F LPT                      1 003008'
        WRITE(fp04,90) 'F ENTRANCE AND COVER       1 009010'	  
        WRITE(fp04,90) 'F SOUFFLET                 2 033033 038038'
        WRITE(fp04,90) 'F VERTICAL PORT            2 039061 068068'
        WRITE(fp04,90) 'F'
        WRITE(fp04,90) 'F'
        WRITE(fp04,90) 'F'
        WRITE(fp04,92) 230.0,230.0, 80.0,0.0,-750.0
        WRITE(fp04,92)  95.0, 95.0,800.0,0.0,   0.0,750.0
        WRITE(fp04,92)  45.0, 20.0
        WRITE(fp04,91) 1,100,1,2,3,4,5,6,9,0,1
        WRITE(fp04,91) 0
      ELSE
        CALL ER('WriteBlock11_04','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock03b
c
c
c
c
      SUBROUTINE WriteBlock03b_04
      USE mod_eirene04
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'slcom'

      INTEGER fp1,fp2,mode

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

      INTEGER    i1,i2,i3,i4,i5,ir,iv,id,cvesm,ndivadsur,nvp,
     .           nvoidside,iliin,iacell,defiliin,
c     .           nvoidside,voidside(MAXPTS),iliin,iacell,defiliin,
     .           defilside,ilside,defilswch,ilswch,fp3,s1,v1,cell2,
     .           ilcell,in,nsdtor,nlimi1,nlimi2,nlimi3,nlimi4,nlimi5,
     .           sbgki,cbgki,haddi,hstdi,nlimi,id0,id1,id2,
     .           addtorseg1,addtorseg2,numseg,iltor,addsection,
     .           numsec,nsec,defilcell,defiltor
      REAL       x0,y0,r,material(4),recyct,recycf,rtmp,ztmp,
     .           x(4),y(4),z(4),
     .           defrecyct,defrecycf,rlbnd,zval,ewall
      LOGICAL    status,done


      INTEGER isec,ishift
      REAL    zshift


      CHARACTER  buffer*200

c      INTEGER    MAXSUR
c      PARAMETER (MAXSUR=MAXASC)
c      INTEGER nsur,csur(MAXSUR),nfunny,isur(MAXSUR)
c      REAL    xsur(4,MAXSUR),ysur(4,MAXSUR),zsur(4,MAXSUR)

c      COMMON  /EIRWALCOM/ walln,wallr,wallz,wallw,ebgki
c      INTEGER walln,wallw(MAXPTS),ebgki
c      REAL    wallr(MAXPTS,2),wallz(MAXPTS,2),wallt(MAXPTS)

      DATA material / 9642., 1206., 18474., 904./

      DATA addsection,addtorseg1,addtorseg2 /0,0,0/


      COMMON /SURFACEMAPCOM/ surfacemap,surfacesrc
      INTEGER surfacemap(5000),surfacesrc(5000)            

      SAVE












      WRITE(fp04,80) '*** 3b. DATA FOR ADDITIONAL SURFACES (DIVIMP)'
      WRITE(fp04,81) 0
c      WRITE(fp04,81) nadd

      DO i1 = 1, 0
c      DO i1 = 1, nadd
        IF (add(i1)%type.EQ.VESSEL_WALL) THEN

c...  NSUR = 1 and NPTS = 2:...
          rlbnd = 2.0

          IF (.TRUE.) THEN
            WRITE(fp04,80) '*'
            WRITE(fp04,83) rlbnd,1.0,1.0E-05,1.0E+05
            WRITE(fp04,81) add(i1)%iliin ,add(i1)%ilside,
     .                     add(i1)%ilswch,0             ,
     .                     add(i1)%iltor ,2             ,
     .                     0             ,add(i1)%ilcell,0

            IF     (rlbnd.EQ.2.0) THEN

              WRITE(fp04,82) (add(i1)%v(1,add(i1)%ipts(i2,1))*100.0,
     .                        add(i1)%v(2,add(i1)%ipts(i2,1))*100.0,  
     .                        add(i1)%v(3,add(i1)%ipts(i2,1))*100.0,  
     .                        i2=1,2)

c            ELSEIF (rlbnd.EQ.4.0) THEN
c              WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=1,2)
c              WRITE(fp2,95) (x(i2),y(i2),z(i2),i2=4,3,-1)
            ELSE
              CALL ER('Block3b','Unknown RLBND value',*99)
            ENDIF

            WRITE(fp04,81) 1,0,0,0
            WRITE(fp04,82) add(i1)%material,add(i1)%ewall
            WRITE(fp04,82) add(i1)%recycf  ,add(i1)%recyct,
     .                    0.0,1.0,0.5,1.0

          ENDIF


        ENDIF
      ENDDO



 80   FORMAT(A)
 81   FORMAT(20(I6:))
 82   FORMAT(1P,20(E12.4:))
 83   FORMAT(1P,20(E12.5:))       






      RETURN
96    CALL ER('WriteInputFile','Cannot create dump file',*99)
97    CALL ER('WriteInputFile','Unexpected end of file',*99)
98    CALL ER('WriteInputFile','Problems reading template file',*99)
99    WRITE(50,*) '  Last line read: '
      WRITE(50,*) '  "',buffer,'"'
c99    WRITE(EROUT,*) '  Last line read: '
c      WRITE(EROUT,*) '  "',buffer,'"'
      STOP
      END




c
c ======================================================================
c
c subroutine: WriteBlock03a_04
c
c
c
c
      SUBROUTINE WriteBlock03a_04
      USE mod_eirene04
      IMPLICIT none

      INTEGER i1,nstsi,instsi

c      WRITE(0,*) 'WRITING BLOCK3a'

c...  Count the number of non-default standard surfaces that have been
c     defined:
      nstsi = 0
      DO i1 = 1, nadd
c        WRITE(0,*) 'NSTSI=',nstsi 
        IF (add(i1)%type.EQ.NON_DEFAULT_STANDARD) nstsi = nstsi + 1
      ENDDO

      WRITE(fp04,90) '*** 3a. DATA FOR NON-DEFAULT SURFACES (OSM)'
      WRITE(fp04,91) nstsi
      instsi = 0
      DO i1 = 1, nadd
        IF (add(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
        instsi = instsi + 1

c        fp04 = 0

        WRITE(fp04,90) add(i1)%surtxt(1:LEN_TRIM(add(i1)%surtxt))
        WRITE(fp04,91) instsi,1,1
        WRITE(fp04,91) add(i1)%iliin ,add(i1)%ilside,
     .                 add(i1)%ilswch,0             ,
     .                 add(i1)%iltor ,add(i1)%ilcol ,
     .                 0,0,0,0

        IF (add(i1)%reflect.EQ.LOCAL) THEN
          WRITE(fp04,91) 1,0
          WRITE(fp04,92) add(i1)%material,add(i1)%ewall 
          WRITE(fp04,92) add(i1)%recycf,add(i1)%recyct,0.0,1.0,0.5,1.0
        ENDIF
      ENDDO

c      WRITE(0,*) 'DONE'

c      STOP 'sdfsd'

90    FORMAT(A)
91    FORMAT(20(I6:))
92    FORMAT(1P,20(E12.4:))
93    FORMAT(1P,20(E12.5:))

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock05
c
c
c
c
      SUBROUTINE WriteBlock05_04
      USE mod_eirene04
      IMPLICIT none

      INTEGER   fp1,fp2,nmass,i1,ibgk,nplsi,iscd2
      CHARACTER buffer*200
      CHARACTER*8 psym(2,5)

c     nmass = NINT(crmb)
      nmass = NINT(2.0) ! Replace with assignment of mod_eirene04 variable with CRMB...

      IF (nmass.NE.1.AND.nmass.NE.2) 
     .  CALL ER('WriteBlock05','EIRENE can only run with H and D',*99)

      psym(1,1) = 'H+'
      psym(2,1) = 'D+'
      psym(1,2) = 'H(B)'
      psym(2,2) = 'D(B)'
      psym(1,3) = 'H2(B)'
      psym(2,3) = 'D2(B)'
      psym(1,4) = 'HH2(B)'
      psym(2,4) = 'DD2(B)'
      psym(1,5) = 'H2H(B)'
      psym(2,5) = 'D2D(B)'

      iscd2 = 214

      nplsi = 1
      IF (photons.EQ.1.OR.photons.EQ.2) nplsi = nplsi + 13
      IF (bgk.EQ.3) nplsi = nplsi + 4

      WRITE(fp04,90) '*** 5. DATA FOR PLASMA-BACKGROUND (OSM)'
      WRITE(fp04,90) '*BULK ION SPECIES CARDS:  NPLSI ION SPECIES '//
     .               'ARE CONSIDERED, NPLSI='
      WRITE(fp04,91) nplsi
      WRITE(fp04,94) 1,psym(nmass,1),nmass,1,1,1,1,-4,0,1      ! D+
      WRITE(fp04,91) 8,115,111
!      WRITE(fp04,92) 0.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
      WRITE(fp04,92) 16.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
c      WRITE(fp04,92) 16.0,0.0,0.0,0.0,1.0E-15                      ! No volume recombination

      ibgk = 0

      IF (bgk.EQ.3) THEN
c...    BGK:
        ibgk = 4
        iscd2 = 614
        WRITE(fp04,94) 2,psym(nmass,2),  nmass,1,1,0,1,-1,0,0
        WRITE(fp04,94) 3,psym(nmass,3),2*nmass,2,2,0,1,-1,0,0
        WRITE(fp04,94) 4,psym(nmass,4),  nmass,1,1,0,1,-1,0,0
        WRITE(fp04,94) 5,psym(nmass,5),2*nmass,2,2,0,1,-1,0,0
      ENDIF

      IF (photons.EQ.1.OR.photons.EQ.2) THEN
c...    Photons:

        IF (photons.EQ.1) THEN
          WRITE(fp04,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0
          WRITE(fp04,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1
          WRITE(fp04,91) 19,0,210,iscd2,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        ELSE
          WRITE(fp04,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0,0,0,0,0,0,
     .                   'FORT.13   ',1
          WRITE(fp04,91) 2

          WRITE(fp04,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                   'COLRAD    ',1
          WRITE(fp04,91) 19,0,210,iscd2,0
          WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
          WRITE(fp04,95) 2,'AMJUEL H.122.1.5b   OT'
        ENDIF

        WRITE(fp04,94) 4+ibgk,'D_2c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp04,91) 19,0,210,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp04,95) 1,'AMJUEL H.122.1.8b   OT'

        WRITE(fp04,94) 5+ibgk,'D_3g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp04,91) 20,0,310,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp04,94) 6+ibgk,'D_3c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp04,91) 20,0,310,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp04,95) 1,'AMJUEL H.122.1.8a   OT'

        WRITE(fp04,94) 7+ibgk,'D_4g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp04,91) 25,0,410,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp04,94) 8+ibgk,'D_4c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp04,91) 25,0,410,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp04,95) 1,'AMJUEL H.122.1.8c   OT'

        WRITE(fp04,94) 9+ibgk,'D_5g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp04,91) 27,0,510,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp04,94) 10+ibgk,'D_5c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp04,91) 27,0,510,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp04,95) 1,'AMJUEL H.122.1.8d   OT'

        WRITE(fp04,94) 11+ibgk,'D_6g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp04,91) 29,0,610,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp04,94) 12+ibgk,'D_6c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp04,91) 29,0,610,iscd2,0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp04,95) 1,'AMJUEL H.122.1.8e   OT'

        WRITE(fp04,94) 13+ibgk,'Lyman_a ',nmass,1,1,0,1,-4,0,0
        WRITE(fp04,94) 14+ibgk,'Lyman_b ',nmass,1,1,0,1,-4,0,0
      ENDIF

c...
      WRITE(fp04,91) 5,-5,5,5,5
      WRITE(fp04,92) 1.0, 90.0,1.5,2.5,4.3,72.0             ! Te (bogus, not actually used in eirene)
      DO i1 = 1, nplsi
        WRITE(fp04,92) 1.0,200.0,3.0,1.0,4.3,72.0,REAL(i1)  ! Ti
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp04,92) 0.0,0.0,3.0,1.0,4.6,72.0,REAL(i1)    ! ni
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0, 0.0,REAL(i1)    ! vi
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0, 0.0
        WRITE(fp04,92) 0.0,0.0,0.0,0.0,0.0, 0.0
      ENDDO
      WRITE(fp04,92) 0.1,  0.1,0.0,0.0,0.0,72.0             ! B-field


c      IF (.FALSE..AND.bgk.EQ.0.AND.photons.EQ.0) THEN
c...    Standard (no BGK or photons):
c        WRITE(fp04,91) 1
c        WRITE(fp04,91) 5,-5,5,5,5
c        WRITE(fp04,92) 1.0, 90.0,1.5,2.5,4.3,72.0
c        WRITE(fp04,92) 1.0,200.0,3.0,1.0,4.3,72.0
c        WRITE(fp04,92) 0.0,  0.0,3.0,1.0,4.6,72.0
c        WRITE(fp04,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp04,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp04,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp04,92) 0.1,  0.1,0.0,0.0,0.0,72.0
c      ELSE
c        CALL ER('WriteBlock05','Invalid EIR_07OPT',*99)
c      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6))
92    FORMAT(1P,20(E12.4))
94    FORMAT(I2,1X,A8,12(I3),1X,A10,1X,I2)
95    FORMAT(I6,1X,A)

      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock06
c
c
c
c
      SUBROUTINE WriteBlock06_04(fp1,fp2)
      USE mod_eirene04
      IMPLICIT none

c      INCLUDE 'params'
c      INCLUDE 'cgeom'
c      INCLUDE 'comtor'
c      INCLUDE 'pindata'
c      INCLUDE 'slcom'

      INTEGER   fp1,fp2
      CHARACTER buffer*200

      WRITE(fp2,90) '*** 6. DATA FOR GENERAL REFLECTION MODEL (DIVIMP)'

      WRITE(fp2,90) 'TF'

      IF (trim.EQ.1) THEN
c. *HARDCODED* Need to read the DIVIMP execution directory from an enviroment
c              variable:
        WRITE(fp2,90) 'PATH  ./TRIM/'
        WRITE(fp2,90) 'D_on_Mo'               
        WRITE(fp2,90) 'D_on_Fe'               
        WRITE(fp2,90) 'D_on_C'               
        WRITE(fp2,90) 'D_on_Be'               
      ENDIF

      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0
      WRITE(fp2,91) 1.0,50.0,0.1

c      WRITE(fp2,91) eirermin,50.0,0.1,eirrinteg,eireinteg

c...  Advance input stream to the start of the next input block:
40    CALL ReadLine(fp1,buffer,1,*97,*98)
      IF (buffer(1:3).NE.'***') GOTO 40
      BACKSPACE fp1

      RETURN
90    FORMAT(A)
91    FORMAT(1P,10(E12.4:),0P)
97    CALL ER('WriteBlock06','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock06','Problems reading template file',*99)
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock07
c
c
c
c
c
c ======================================================================
c
c Write EIRENE geometry file:
c
c ======================================================================
c
c subroutine: WriteGeometryFile
c
c
c
c
c ======================================================================
c
c
c ======================================================================
c
c
c ======================================================================
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
              WRITE(0,*) 'WARNING NextLine: EIRENE data beyond '//
     .                   'single precision size limit, setting to zero'
              rdum(i1) = 0.0
c              STOP 'NOT SURE WHAT TO DO HERE'
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


