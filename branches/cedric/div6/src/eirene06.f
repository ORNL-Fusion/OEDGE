c     -*-Fortran-*-
c
c ======================================================================
c
c OSM-EIRENE06 interface:
c
c Write input file
c Write triangle
c     -WriteEIRENE_06
c
c ======================================================================
c
      SUBROUTINE LoadTriangles_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2
      REAL    version

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ(fp,ERR=98) version,ntri,nver,nsurface

      IF (version.NE.1.0)
     .  CALL ER('LoadTriangles','Unsupporting version',*99)

      CALL ALLOC_VERTEX(nver)
      CALL ALLOC_SURFACE(nsurface)
      CALL ALLOC_TRIANGLE(ntri)
      READ(fp,ERR=98) (tri(i1),i1=1,ntri)
      READ(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      READ(fp,ERR=98) (surface(i1),i1=1,nsurface)

c      READ(fp,ERR=98) tri,ver,add
      CLOSE (fp)
      
      RETURN
 98   CALL ER('LoadTriangles','Problems reading data file',*99)
 99   STOP
      END
c
c ======================================================================
c
      SUBROUTINE GetObjCentre(iobj,cen)
      USE mod_eirene06
c      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,ivtx
c      INTEGER iobj,iside,isrf,ivtx
      REAL*8  cen(3),count


      cen = 0.0D0
      DO ivtx = 1, 3
        cen(1:3) = cen(1:3) + ver(tri(iobj)%ver(ivtx),1:3)
      ENDDO
      cen = cen / 3.0D0

c      count = 0.0D0
c      DO iside = 1, obj(iobj)%nside
c        isrf = ABS(obj(iobj)%iside(iside))
c        DO ivtx = 1, srf(isrf)%nvtx
c          count = count + 1.0D0
c          cen(1:3) = cen(1:3) + vtx(1:3,srf(isrf)%ivtx(ivtx))
c        ENDDO
c      ENDDO
c      cen = cen / count

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteEireneObjects
c
      SUBROUTINE WriteEireneObjects
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      IMPLICIT none

      REAL*8 gmCalcSurfaceArea  ! This will go in mod_geometry...

      INTEGER fp,ivtx,iobj
      LOGICAL found
      REAL    version

      INTEGER   ipts(4),i1,v1,istart,iend,iside,isrf,ipla,
     .          iobj1(4),iside1(4),isrf1(4),
     .          ik,ir,it,ctardat,max_ik,max_ir,max_is,max_plasma,
     .          tar_maxik,tar_maxir,itarget
      REAL      sumflux1,sumflux2,frac
      INTEGER, ALLOCATABLE :: tar_objects(:)
      REAL   , ALLOCATABLE :: tdata(:)      
      REAL*8 , ALLOCATABLE :: tar_area(:), tar_totarea(:,:)

      WRITE(eirfp,*) 'WRITING EIRENE OBJECT FILES'

      version = 1.00

      fp = 99

      ALLOCATE(tdata(nobj))

      IF (photons.EQ.-1) THEN
c...    Load ionisation data from previous EIRENE call:
        CALL LoadTriangleData(7,0,13,0,tdata)  
      ELSE
        tdata = -999.0
      ENDIF

c...  Dump vertices:
      OPEN(UNIT=fp,FILE='objects.npco_char',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      istart = 1
      iend   = nvtx
      IF (nvtxmap.NE.0) THEN
        WRITE(fp,*) nvtxmap
        DO ivtx = istart, iend
          IF (vtxmap(ivtx).NE.0) THEN
            WRITE(fp,'(I10,3F14.6)') vtxmap(ivtx),vtx(1:3,ivtx)*100.0D0 ! *** Increase accuracy...? ***
          ENDIF
        ENDDO
      ELSE 
        WRITE(fp,*) iend-istart+1
        DO ivtx = istart, iend
          WRITE(fp,'(I6,3F14.6)') ivtx-istart+1,vtx(1:3,ivtx)*100.0D0  ! *** Increase accuracy...? ***
        ENDDO
      ENDIF
      CLOSE(fp)      

c...  Dump sides:
      OPEN(UNIT=fp,FILE='objects.elemente',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      

      istart = 1
      iend   = nobj

      max_ik     = -1
      max_ir     = -1
      max_is     = -1
      max_plasma = -1
      DO iobj = istart, iend
        max_ik     = MAX(max_ik    ,obj(iobj)%index(IND_IK))
        max_ir     = MAX(max_ir    ,obj(iobj)%index(IND_IR))
        max_is     = MAX(max_is    ,obj(iobj)%index(IND_IS))
        max_plasma = MAX(max_plasma,obj(iobj)%index(IND_PLASMA))
      ENDDO

      WRITE(fp,'(5I9)') iend-istart+1,max_ik,max_ir,max_is,max_plasma
      DO iobj = istart, iend
c...    Get vertices associated with the base triangle (always side #1 if 
c       grid prepared properly):
        isrf = obj(iobj)%iside(1)
        IF (isrf.LT.0) THEN
c...      Need to reverse order of points:
          ipts(1) = srf(-isrf)%ivtx(3)  ! This must match the convention in CheckTetrahedronStructure
          ipts(2) = srf(-isrf)%ivtx(2)
          ipts(3) = srf(-isrf)%ivtx(1)
        ELSE
          ipts(1:3) = srf(isrf)%ivtx(1:3)
        ENDIF
c...    Forth point, or apex point of sorts, from any of the other sides:
        isrf = ABS(obj(iobj)%iside(2))
        DO ivtx = 1, 3
          IF (srf(isrf)%ivtx(ivtx).NE.ipts(1).AND.
     .        srf(isrf)%ivtx(ivtx).NE.ipts(2).AND.
     .        srf(isrf)%ivtx(ivtx).NE.ipts(3)) EXIT
        ENDDO
        IF (ivtx.EQ.4) 
     .    CALL ER('WriteEireneObjects','4th tetrahedron point '//
     .            'not identified',*99)
        ipts(4) = srf(isrf)%ivtx(ivtx)
        IF (nvtxmap.NE.0) THEN  
          DO ivtx = 1, 4
            IF (vtxmap(ipts(ivtx)).EQ.0) THEN
              STOP 'STOP: PROBLEM WITH VTXMAP'
            ELSE
              ipts(ivtx) = vtxmap(ipts(ivtx))
            ENDIF
          ENDDO          
        ENDIF
        WRITE(fp,'(I9,4X,4I8,4X,4I3,4X,4I6)') 
     .    iobj-istart+1,
     .    (ipts(v1),v1=1,4),  ! (tri(i1)%sideindex(1,v1),v1=1,3),
     .    (0       ,v1=1,4),  ! (tri(i1)%sideindex(2,v1),v1=1,3),
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    obj(iobj)%index(IND_IS),
     .    obj(iobj)%index(IND_PLASMA)
      ENDDO
      CLOSE(fp)      

c...  Connection map:
      OPEN(UNIT=fp,FILE='objects.neighbors',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) nobj
      DO iobj = istart, iend
c...    Collect connection map information:
        DO iside = 1, obj(iobj)%nside
          iobj1 (iside) = obj(iobj)%omap(iside)                    
          iside1(iside) = obj(iobj)%smap(iside)
          isrf = ABS(obj(iobj)%iside(iside))
          isrf1 (iside) = srf(isrf)%index(IND_SURFACE)             ! Surface (block 2A)
        ENDDO
        IF (tetrahedrons.AND.
     .      iobj1(1).EQ.0.AND.iside1(1).EQ.0.AND.isrf1(1).EQ.0) THEN   ! *HACK* temporary for toroidal end surfaces...
          isrf1(1) = 8  
        ENDIF
c        IF (iobj1(2).EQ.0.AND.iside1(2).EQ.0.AND.isrf1(2).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(2) = 8 
c        ENDIF
c        IF (iobj1(3).EQ.0.AND.iside1(3).EQ.0.AND.isrf1(3).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(3) = 8 
c        ENDIF
c        IF (iobj1(4).EQ.0.AND.iside1(4).EQ.0.AND.isrf1(4).EQ.0) THEN   ! *HACK* temporary for toroidal surfaces...
c          isrf1(4) = 8 
c        ENDIF
c        WRITE(fp,'(I9,4X,4(I9,F14.3,I4,4X),2I6,2X,2I4)') iobj-istart+1,
        WRITE(fp,'(I9,4X,4(I9,I4,I4,4X),2I6,4X,2I4)') iobj-istart+1,
     .    (iobj1(v1),iside1(v1),isrf1(v1),v1=1,4),  
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    0,0
c        WRITE(fp,'(13I10)') iobj-istart+1,
c     .    (iobj1(v1),iside1(v1),isrf1(v1),v1=1,4)
      ENDDO
      CLOSE(fp)      


      WRITE(eirfp,*) ' *** VFIELD AND BFIELD NEED CORRECTING ***'

c...  Plasma data:
      OPEN(UNIT=fp ,FILE='objects.plasma',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A,F4.2,A)') 
     .  '* VERSION ',version,' OF '//
     .  fluid_code(1:LEN_TRIM(fluid_code))//
     .  ' PLASMA FILE FOR TETRAHEDRON GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A9,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*   Index','Te','Ti','ne','vx',
     .  'vy','vz','Bx','By','Bz'
      WRITE(fp,'(A9,2A8,2X,A10,2X,3A10,2X,3A12)') 
     .  '*        ','(eV)','(eV)','(cm-3)','(cm s-1)',
     .  '(cm s-1)','(cm s-1)','(Tesla)','(Tesla)','(Tesla)'
      WRITE(fp,*) nobj
      DO iobj = istart, iend
        ipla = obj(iobj)%index(IND_PLASMA)
        WRITE(fp,'(I9,2F8.2,2X,1P,E10.2,2X,3E10.2,2X,
     .             3E12.4,2X,E12.4,0P,6X,3I4,I6)') iobj-istart+1,
     .    plasma(1,ipla),           ! Te (eV)
     .    plasma(2,ipla),           ! Ti (eV)
     .    plasma(3,ipla)*1.0E-06,   ! ne (cm-3)
     .    plasma(4,ipla)*100.0,     ! vx (cm s-1)
     .    plasma(5,ipla)*100.0,     ! vy 
     .    plasma(6,ipla)*100.0,     ! vz
     .    bfield(1,ipla),           ! Bx (Tesla)
     .    bfield(2,ipla),           ! By
     .    bfield(3,ipla),           ! Bz
     .    tdata(iobj),  ! Ionisation rate from previous run (w or w/o photons)...
     .    obj(iobj)%index(IND_IK),
     .    obj(iobj)%index(IND_IR),
     .    obj(iobj)%index(IND_IS),
     .    obj(iobj)%index(IND_PLASMA)
      ENDDO

      sumflux1= 0.0
      sumflux2= 0.0

c...  Target data:
      ctardat = 0
      tar_maxik = 0
      tar_maxir = 0
c...  First count up how many surfaces are targets:
      DO iobj = istart, iend
        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
        DO iside = 1, obj(iobj)%nside  ! This will always be iside=1 for now...
          isrf = ABS(obj(iobj)%iside(iside))
          IF (srf(isrf)%index(IND_TARGET).NE.0) THEN 
            ctardat = ctardat + 1
            tar_maxik = MAX(tar_maxik,obj(iobj)%index(IND_IK))
            tar_maxir = MAX(tar_maxir,obj(iobj)%index(IND_IR))
          ENDIF
        ENDDO
      ENDDO
      IF (ctardat.EQ.0.OR.tar_maxik.EQ.0.OR.tar_maxir.EQ.0) 
     .  CALL ER('WriteEireneObjects','Problem with target mapping',*99)
c...  Dynamic allocation because this number could be large for tetrahedrons, and
c     then loop again over all objects, recording the relevant information:
      ALLOCATE(tar_objects(ctardat))
      ALLOCATE(tar_area   (ctardat))
      ALLOCATE(tar_totarea(tar_maxik,tar_maxir))
      ctardat = 0
      tar_totarea = 0.0D0
      DO iobj = istart, iend
        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
        DO iside = 1, obj(iobj)%nside 
          isrf = ABS(obj(iobj)%iside(iside))
          IF (srf(isrf)%index(IND_TARGET).NE.0) THEN
            ctardat = ctardat + 1
            tar_objects(ctardat) = iobj
            tar_area   (ctardat) = gmCalcSurfaceArea(isrf)
            ik = obj(iobj)%index(IND_IK)
            ir = obj(iobj)%index(IND_IR) 
            tar_totarea(ik,ir) = tar_totarea(ik,ir) + tar_area(ctardat)
          ENDIF
        ENDDO
      ENDDO
      WRITE(fp,'(A)') '* TARGET DATA'
      WRITE(fp,*) ctardat
c helium
c      WRITE(fp,*) ctardat*2
c      DO iobj = istart, iend
c        IF (grp(obj(iobj)%group)%origin.NE.GRP_MAGNETIC_GRID) CYCLE
      DO itarget = 1, ctardat
        iobj = tar_objects(itarget)
        ik = obj(iobj)%index(IND_IK)
        ir = obj(iobj)%index(IND_IR) 
        DO iside = 1, obj(iobj)%nside
          isrf = ABS(obj(iobj)%iside(iside))
          IF (srf(isrf)%index(IND_TARGET).EQ.0) CYCLE
c...      Find corresponding target data, as set in ProcessFluidCode, based
c         on the fluid code cell/ring indices:
          found = .FALSE.
          DO it = 1, ntardat
            IF (tardat(it,2).NE.REAL(ik).OR.
     .          tardat(it,3).NE.REAL(ir)) CYCLE
            IF (.NOT.found) THEN
              found = .TRUE.         
              frac = SNGL(tar_area(itarget) / tar_totarea(ik,ir))
c              WRITE(0,*) 'FRAC:',frac,ik,ir
              WRITE(fp,'(I9,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
     .                   F6.2,1P,E10.2,0P,6X,3I4)') 
c...            Target quantities:
     .          iobj-istart+1,iside,    ! Triangle index and side index
     .          tardat(it,7 )*frac,     ! ion flux to surface for species 1 (Amps)
     .          tardat(it,6 ),          ! Te (eV)                                 
     .          tardat(it,8 ),          ! Ti (ev)
     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
     .          tardat(it,11),          ! Mach no.        (not read)
     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
     .          ik,ir,iside             ! Fluid grid indices, for debugging only

c helium
c              WRITE(fp,'(I9,I6,1P,E10.2,0P,2F8.2,1P,2E10.2,0P,
c     .                   F6.2,1P,E10.2,0P,6X,2I4)') 
c...            Target quantities:
c     .          iobj-istart+1,iside,    ! Triangle index and side index
c     .          tardat(it,7 ),          ! ion flux to surface for species 1 (Amps)
c     .          tardat(it,6 ),          ! Te (eV)
c     .          tardat(it,8 ),          ! Ti (ev)
c     .          tardat(it,9 )*1.0E-06,  ! ni (cm-3)
c     .          tardat(it,10)*100.0,    ! v_para (cm s-1) (not read by EIRENE as yet)  
c     .          tardat(it,11),          ! Mach no.        (not read)
c     .          tardat(it,12)*1.0E-04,  ! jsat (A cm-2)   (not read)
c     .          ik,ir                   ! Fluid grid indices, for debugging only
            ELSE
              CALL ER('WriteEireneObjects','Target data appears to  '//
     .                'be over-specified',*99)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      DEALLOCATE(tar_objects)
      DEALLOCATE(tar_area)
      DEALLOCATE(tar_totarea)
      CLOSE(fp)      

c      OPEN(UNIT=fp,FILE='objects.efield',ACCESS='SEQUENTIAL',
c     .     STATUS='REPLACE',ERR=96)      
c...  Header:
c      WRITE(fp,'(A)') '* VERSION 1.0 OSM E-FIELD FILE FOR TRIANGULAR '//
c     .                'GRID'
c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A)') '* BULK PLASMA DATA (PART DEUX)'
c      WRITE(fp,'(A)') '*'
c      WRITE(fp,'(A7,3A12)') 
c     .  '* Index','Ex','Ey','Ez'
c      WRITE(fp,'(A7,3A12)')
c     .  '*      ','V m-1','V m-1','V m-1'
c      WRITE(fp,*) nobj
c      DO i1 = 1, nobj
cc...    Dump efield data:
c        WRITE(fp,'(I7,1P,3E12.4,10X,0P,2I4)') i1,
c     .    obj(i1)%efield(1),
c     .    obj(i1)%efield(2),
c     .    obj(i1)%efield(3),
c     .    obj(i1)%index(1),obj(i1)%index(2) 
c      ENDDO
c      CLOSE(fp)      

      WRITE(eirfp,*) 'DONE'

      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)
      IF (ALLOCATED(vtxmap)) DEALLOCATE(vtxmap)
      IF (ALLOCATED(plasma)) DEALLOCATE(plasma)
      IF (ALLOCATED(bfield)) DEALLOCATE(bfield)

      RETURN
96    WRITE(0,*) 'WRITETRIANGEFILES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END
c
c
c ======================================================================
c
c  subroutine: BinItems
c
c
c ...don't do this on the fly because...
c
      SUBROUTINE BinItems(n,cen,nx,ny,nz,nlist,ilist,glist)
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN)  :: n,nx,ny,nz
      REAL                 :: cen(3,n)
      INTEGER, INTENT(OUT) :: nlist(nx,ny,nz),ilist(nx,ny,nz),glist(n)


      INTEGER i,ix,iy,iz,count
      REAL    xmin,xmax,ymin,ymax,zmin,zmax,x,z,dx,dy,dz

      REAL   , PARAMETER :: PI = 3.1415926536, TOL = 1.0E-07

c...  Convert x coordinate to r and z coordinate to phi:
      IF (.TRUE.) THEN
c     
c       z <---------|
c                 / |
c                /TH|
c               /   |
c              /   .|.
c             /     .
c                   x
c
         DO i = 1, n
           x = cen(1,i)
           z = cen(3,i)
           cen(1,i) = SQRT(x**2 + z**2)
           IF (ABS(x).LT.TOL) THEN
             IF (z.GT.0.0) cen(3,i) =       PI / 2.0
             IF (z.LT.0.0) cen(3,i) = 3.0 * PI / 2.0
           ELSE
             cen(3,i) = ATAN(z / x)
             IF (x       .LT.0.0) cen(3,i) =cen(3,i)+PI      ! Hopefully this is not
             IF (cen(3,i).LT.0.0) cen(3,i) =cen(3,i)+PI*2.0  ! compiler dependant...
           ENDIF
          cen(3,i) = cen(3,i) * 180.0 / PI
          IF (cen(3,i).GT.360.0+TOL) THEN
            CALL ER('BinItems','PHI > 360.0',*99)
          ENDIF
          IF ((ABS(cen(3,i)      ).LT.TOL).OR.
     .        (ABS(cen(3,i)-360.0).LT.TOL)) cen(3,i) = 0.0 
        ENDDO

      ENDIF


      nlist = 0
      ilist = 0
      glist = 0


      WRITE(geofp,*) '    BINNING'

c.... Find extents:
      xmin =  1.0E+20
      xmax = -1.0E+20
      ymin =  1.0E+20
      ymax = -1.0E+20
      zmin =  1.0E+20
      zmax = -1.0E+20
      DO i = 1, n
        xmin = MIN(xmin,cen(1,i))
        xmax = MAX(xmax,cen(1,i))
        ymin = MIN(ymin,cen(2,i))
        ymax = MAX(ymax,cen(2,i))
        zmin = MIN(zmin,cen(3,i))
        zmax = MAX(zmax,cen(3,i))
      ENDDO
      xmin = xmin * (1.0 - 0.001 * SIGN(1.0,xmin))
      xmax = xmax * (1.0 + 0.001 * SIGN(1.0,xmax))
      ymin = ymin * (1.0 - 0.001 * SIGN(1.0,ymin))
      ymax = ymax * (1.0 + 0.001 * SIGN(1.0,ymax))
      zmin = zmin * (1.0 - 0.001 * SIGN(1.0,zmin))
      zmax = zmax * (1.0 + 0.001 * SIGN(1.0,zmax))
      dx = (xmax - xmin) / REAL(nx)      
      dy = (ymax - ymin) / REAL(ny)      
      dz = (zmax - zmin) / REAL(nz)      
      IF (ABS(dz).LT.1.0E-6) dz = 1.0

c...  Assign zones (a bit painful doing it this way, but can't think
c     of how else to produce GLIST as a linear array, minimizing 
c     memory requirements):
      DO i = 1, n
        ix = INT((cen(1,i) - xmin) / dx) + 1
        iy = INT((cen(2,i) - ymin) / dy) + 1
        iz = INT((cen(3,i) - zmin) / dz) + 1
        nlist(ix,iy,iz) = nlist(ix,iy,iz) + 1         
      ENDDO

      count = 1
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx
            ilist(ix,iy,iz) = count
            count = count + nlist(ix,iy,iz)
          ENDDO
        ENDDO
      ENDDO

      WRITE(geofp,*) '    CHECK:',ilist(nx,ny,nz)+nlist(nx,ny,nz)-1,n

      nlist = 0
      DO i = 1, n
        ix = INT((cen(1,i) - xmin) / dx) + 1
        iy = INT((cen(2,i) - ymin) / dy) + 1
        iz = INT((cen(3,i) - zmin) / dz) + 1
        nlist(ix,iy,iz) = nlist(ix,iy,iz) + 1         
        glist(ilist(ix,iy,iz)+nlist(ix,iy,iz)-1) = i
      ENDDO      

      WRITE(geofp,*) '    DONE'


      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: RemoveDuplicateSurfaces
c
c
c ...don't do this on the fly because...
c
      SUBROUTINE RemoveDuplicateSurfaces
      USE mod_geometry
      IMPLICIT none

      INTEGER iobj,isrf,isrf1,ivtx,i1,i2,nx,ny,nz,ix,iy,iz,count,iside
      INTEGER, ALLOCATABLE :: nlist(:,:,:),ilist(:,:,:),glist(:),
     .                        mlist(:)

      REAL, ALLOCATABLE :: cen(:,:)  ! Does this really help memory management..?


      LOGICAL, PARAMETER :: aggressive = .FALSE.
      REAL   , PARAMETER :: PI = 3.1415926536, TOL = 1.0E-07


      WRITE(geofp,*) '  REMOVING DUPLICATE SURFACES'


! NEED SURFACE MAPPING INFO FIRST...  NO YOU DON'T SINCE NEEDED INFO
! IS RECORDED WHEN SURFACE IS CREATED -- OBJ ASSOCIATION! 

      SELECTCASE (0)
        CASE(0)

          WRITE(geofp,*) '    FINDING SURFACE CENTRES'

!...      Find geometric center of surface vertices:
          ALLOCATE(cen(3,nsrf))
          cen = 0.0D0
          DO isrf = 1, nsrf
            DO ivtx = 1, srf(isrf)%nvtx
              cen(1:3,isrf) =cen(1:3,isrf)+vtx(1:3,srf(isrf)%ivtx(ivtx))
            ENDDO
            cen(1:3,isrf) = cen(1:3,isrf) / REAL(srf(isrf)%nvtx)
          ENDDO


          nx = 30
          ny = 30
          nz = 23  ! 10

          ALLOCATE(nlist(nx,ny,nz))
          ALLOCATE(ilist(nx,ny,nz))
          ALLOCATE(glist(nsrf))
          ALLOCATE(mlist(nsrf))
          mlist = 0

          CALL BinItems(nsrf,cen,nx,ny,nz,nlist,ilist,glist)


          WRITE(geofp,*) '    SEARCHING FOR DUPLICATES'

c...      Search for matching surfaces:
          DO iz = 1, nz
            DO iy = 1, ny
              DO ix = 1, nx

                DO i1 = 0, nlist(ix,iy,iz)-2
                  isrf = glist(ilist(ix,iy,iz) + i1)

c                  IF (srf(isrf)%nvtx.NE.3) CYCLE
            
                  IF (mlist(isrf).NE.0) CYCLE
                  DO i2 = i1+1, nlist(ix,iy,iz)-1
                    isrf1 = glist(ilist(ix,iy,iz) + i2)
                    IF (MatchSurface(isrf,isrf1)) mlist(isrf1) = isrf
                  ENDDO
                ENDDO 

              ENDDO ! IX
            ENDDO   ! IY
          ENDDO     ! IZ

          ix = 0
          DO isrf = 1, nsrf
            IF (mlist(isrf).NE.0) ix = ix + 1
          ENDDO
          WRITE(geofp,*) '    ',ix,' DUPLICATES OF',nsrf,' FOUND'

          WRITE(geofp,*) '    UPDATING OBJECTS'

c...      Update object surface pointers:
          DO iobj = 1, nobj        
            DO iside = 1, obj(iobj)%nside            
              isrf = ABS(obj(iobj)%iside(iside))
              IF (mlist(isrf).NE.0) obj(iobj)%iside(iside)= -mlist(isrf)
            ENDDO
          ENDDO

c...      Update surface links:
          DO isrf = 1, nsrf
            IF (srf(isrf)%link.NE.0) THEN            
              WRITE(0,*) 'SOME WORK NEEDED HERE TO SORT OUT '//
     .                   'LINK UPDATES, STOPPING'
              STOP
            ENDIF
          ENDDO

c...      Delete extraneous surfaces:
          IF (aggressive) THEN
            WRITE(geofp,*) '    DELETING SURFACES'
          ENDIF

          IF (.FALSE.) THEN
            WRITE(geofp,*) '    CHECKING'
            DO isrf = 150000, nsrf-1
              DO isrf1 = isrf+1, nsrf
                IF (MatchSurface(isrf,isrf1).AND.mlist(isrf1).EQ.0) THEN
                  WRITE(geofp,*) 'CURSES: DUPLICATE FOUND'
                  WRITE(geofp,*) '  ISRF,1:',isrf,isrf1

                  DO iz = 1, nz
                    DO iy = 1, ny
                      DO ix = 1, nx
                        DO i1 = 0, nlist(ix,iy,iz)-1
                          IF (isrf .EQ.glist(ilist(ix,iy,iz)+i1)) 
     .                      WRITE(geofp,*) ' 0:',ix,iy,iz,cen(1:3,isrf)
                          IF (isrf1.EQ.glist(ilist(ix,iy,iz)+i1)) 
     .                      WRITE(geofp,*) ' 1:',ix,iy,iz,cen(1:3,isrf1)
                        ENDDO
                      ENDDO
                    ENDDO
                  ENDDO

                ENDIF
              ENDDO
            ENDDO
          ENDIF

        CASE DEFAULT
          STOP 'UNRECOGNIZED CASE IN DELETESURFACE'
      ENDSELECT

c...  Clear memory (move above once debugging is done?):
      IF (ALLOCATED(cen  )) DEALLOCATE(cen)     
      IF (ALLOCATED(nlist)) DEALLOCATE(nlist)
      IF (ALLOCATED(ilist)) DEALLOCATE(ilist)
      IF (ALLOCATED(glist)) DEALLOCATE(glist)
      IF (ALLOCATED(mlist)) DEALLOCATE(mlist)

      WRITE(geofp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: RemoveDuplicateVertices
c
      SUBROUTINE RemoveDuplicateVertices
      USE mod_geometry
      IMPLICIT none

      INTEGER ivtx,ivtx1,isrf,i1,count,nx,ny,nz,ix,iy,iz,i2,nremove
      INTEGER, ALLOCATABLE :: nlist(:,:,:),ilist(:,:,:),glist(:),
     .                        mlist(:)

      REAL, ALLOCATABLE :: cen(:,:)  ! Does this really help memory management..?

      LOGICAL, PARAMETER :: aggressive = .FALSE.

      WRITE(geofp,*) '  REMOVING DUPLICATE VERTICES'

!...  Find geometric center of surface vertices:
      ALLOCATE(cen(3,nvtx))
      cen = 0.0D0
      DO ivtx = 1, nvtx
        cen(1:3,ivtx) = vtx(1:3,ivtx)
      ENDDO


      nx = 30
      ny = 30
      nz = 23  ! 10

      ALLOCATE(nlist(nx,ny,nz))
      ALLOCATE(ilist(nx,ny,nz))
      ALLOCATE(glist(nvtx))
      ALLOCATE(mlist(nvtx))
      nlist = 0
      ilist = 0
      glist = 0
      mlist = 0

      CALL BinItems(nvtx,cen,nx,ny,nz,nlist,ilist,glist)

c...  Search for matching vertices:
      DO iz = 1, nz
        DO iy = 1, ny
          DO ix = 1, nx

            DO i1 = 0, nlist(ix,iy,iz)-2
              ivtx = glist(ilist(ix,iy,iz) + i1)
              IF (mlist(ivtx).NE.0) CYCLE
              DO i2 = i1+1, nlist(ix,iy,iz)-1
                ivtx1 = glist(ilist(ix,iy,iz) + i2)

                IF (mlist(ivtx1).EQ.0.AND.
     .              DABS(vtx(1,ivtx)-vtx(1,ivtx1)).LT.1.0D-07.AND.
     .              DABS(vtx(2,ivtx)-vtx(2,ivtx1)).LT.1.0D-07.AND.
     .              DABS(vtx(3,ivtx)-vtx(3,ivtx1)).LT.1.0D-07) 
     .            mlist(ivtx1) = ivtx  ! Tag for reassignment

              ENDDO
            ENDDO 

          ENDDO ! IX
        ENDDO   ! IY
      ENDDO     ! IZ


      nremove = 0
      DO ivtx = 1, nvtx
        IF (mlist(ivtx).NE.0) nremove = nremove + 1
      ENDDO
      WRITE(geofp,*) '    REMOVING ',nremove,' OF',nvtx


c...  Replace reference to tagged vertices in surface lists:
      DO isrf = 1, nsrf
        DO ivtx = 1, srf(isrf)%nvtx            
          i1 = srf(isrf)%ivtx(ivtx)
          IF (mlist(i1).NE.0) srf(isrf)%ivtx(ivtx) = mlist(i1)
        ENDDO
        srf(isrf)%svtx = SUM(srf(isrf)%ivtx(1:srf(isrf)%nvtx))
      ENDDO

      nvtxmap = 0
      IF (aggressive) THEN
        WRITE(geofp,*) '    DELETING DUPLICATE VERTICES...'

c       Leaving this for now...
c        1-build list of vertecies to delete and pass to DeleteVertex
c        2-in DeleteVerted, build index map for non-deleted verteces
c        3-move data to get rid of delted vertices
c        4-do a single scan through surface definitions to get
c          rid of references to deleted verteces
c        the above bit of code can then be deleted I think...
c        adapt for deleteing surfaces quickly/efficiently... 
        DO ivtx = 1, nvtx
 
        ENDDO

c        count = 0
c        ivtx = nvtx + 1
c        DO WHILE (ivtx.GT.1)
c          ivtx = ivtx - 1
c          IF (mlist(ivtx).NE.0) THEN
c            count = count + 1
c            IF (MOD(count,1000).EQ.0) 
c     .        WRITE(geofp,*) '      ',count,' OF',nremove
c            CALL DeleteVertex(ivtx)
c            DO i1 = ivtx, nvtx-1
c              mlist(i1) = mlist(i1+1)
c            ENDDO
c          ENDIF
c        ENDDO
      ELSE
        IF (ALLOCATED(vtxmap)) DEALLOCATE(vtxmap)
        ALLOCATE(vtxmap(nvtx))
        vtxmap = 0
        DO ivtx = 1, nvtx
          IF (mlist(ivtx).EQ.0) THEN
            nvtxmap = nvtxmap + 1
            vtxmap(ivtx) = nvtxmap
          ENDIF 
        ENDDO
      ENDIF


c...  Clear memory (move above once debugging is done?):
      IF (ALLOCATED(cen)) DEALLOCATE(cen)     

      IF (ALLOCATED(nlist)) DEALLOCATE(nlist)
      IF (ALLOCATED(ilist)) DEALLOCATE(ilist)
      IF (ALLOCATED(glist)) DEALLOCATE(glist)
      IF (ALLOCATED(mlist)) DEALLOCATE(mlist)

      WRITE(geofp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: BuildConnectionMap
c
      SUBROUTINE BuildConnectionMap(istart,iend)
      USE mod_eirene06_locals
      USE mod_eirene06
      USE mod_geometry
      IMPLICIT none

      INTEGER, INTENT(IN) :: istart,iend

      INTEGER iobj,iside,isrf,ivtx,iobj1,iside1

      INTEGER, PARAMETER :: ns = 26
C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER s(3,ns)
C     INTEGER s(3,ns) / 0, 0,-1,  0,-1, 0, -1, 0, 0,  
      DATA    s       / 0, 0,-1,  0,-1, 0, -1, 0, 0,  
     .                  0, 0,+1,  0,+1, 0, +1, 0, 0, 
     .                 -1,-1, 0, +1,-1, 0, -1,+1, 0, +1,+1, 0, 
     .                 -1,-1,-1, +1,-1,-1, -1,+1,-1, +1,+1,-1, 
     .                 -1,-1,+1, +1,-1,+1, -1,+1,+1, +1,+1,+1, 
     .                 -1, 0,-1,  0,-1,-1, 
     .                 +1, 0,-1,  0,+1,-1,
     .                 -1, 0,+1,  0,-1,+1, 
     .                 +1, 0,+1,  0,+1,+1 /


      WRITE(geofp,*) '  BUILDING CONNECTION MAP'


c...  Clean up duplicate verticies, necessary for connection
c     map search:
      CALL RemoveDuplicateVertices
      CALL RemoveDuplicateSurfaces

c...  Removing duplicate surfaces essentially builds the
c     connection map since the entire grid is searched
c     to find matching surfaces, in order to eliminate
c     redundant data.  Note however that this will need to
c     be supplimented for local mesh refinement plasma
c     grids since neighbouring cells may no share common
c     surfaces (but all still fine for triangle/tetrahedron
c     grids):


      CALL BuildConnectionMap_New
      RETURN
 

c...  OLD CODE:
      WRITE(eirfp,*) '    MAPPING'
      DO iobj = istart, iend
        obj(iobj)%omap = 0
        obj(iobj)%smap = 0
      ENDDO
      DO iobj = istart, iend
        DO iside = 1, obj(iobj)%nside         
          isrf = obj(iobj)%iside(iside)
          IF (isrf.LT.0) THEN
            iobj1  = srf(-isrf)%obj
            iside1 = srf(-isrf)%side
            obj(iobj )%omap(iside ) = iobj1
            obj(iobj )%smap(iside ) = iside1
            obj(iobj1)%omap(iside1) = iobj 
            obj(iobj1)%smap(iside1) = iside
          ENDIF
        ENDDO
      ENDDO

      WRITE(eirfp,*) '  DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: BuildNewTriangleObjects
c
      SUBROUTINE BuildNewTriangleObjects
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf) newsrf

      INTEGER iver,itri,ivtx,isrf,iside
      REAL*8 a(3)

c...  Copy vertices:
      DO iver = 1, nver
        a(1:2) = DBLE(ver(iver,1:2))
        a(3)   = 0.0D0
        ntryvtx = ntryvtx + 1
        tryvtx(1:3,ntryvtx) = a(1:3)
c        ivtx = AddVertex(a)
      ENDDO 
c...  Build new triangle objects and surfaces:        

      
      ngrp = 2
      grp(1)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(1)%type   = GRP_TRIANGLE
      grp(2)%origin = GRP_VACUUM_GRID
      grp(2)%type   = GRP_TRIANGLE

      DO itri = 1, ntri
c...     
        try(itri)%group             = tri(itri)%type      ! Just works by luck...
        try(itri)%index(IND_IK    ) = tri(itri)%index(1)
        try(itri)%index(IND_IR    ) = tri(itri)%index(2)
        try(itri)%index(IND_IS    ) = 0
        try(itri)%index(IND_ZONE  ) = tri(itri)%zone
c...    Plasma data:
c        IF (try(itri)%index(IND_IK).NE.0) THEN
          try(itri)%index(IND_PLASMA) = itri
c        ELSE
c          try(itri)%index(IND_PLASMA) = 0
c        ENDIF
        plasma(1:20,itri) = tri(itri)%plasma(1:20)
        bfield(1:3 ,itri) = tri(itri)%bfield(1:3)
c...    Create surfaces:
        try(itri)%nside = 3

        DO iside = 1, try(itri)%nside
          try(itri)%omap(iside) = tri(itri)%map(iside)
          try(itri)%smap(iside) = tri(itri)%sid(iside)
c          try(itri)%map(iside) = REAL(tri(itri)%map(iside))+
c     .                           REAL(tri(itri)%sid(iside))/100.0
        ENDDO     

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,1)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,1)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(1)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(1)
        newsrf%ivtx(2) = tri(itri)%ver(2)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(1) = ntrysrf

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,2)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,2)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(2)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(2)
        newsrf%ivtx(2) = tri(itri)%ver(3)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(2) = ntrysrf

        newsrf%index(IND_STDGRD)  = tri(itri)%sideindex(1,3)
        newsrf%index(IND_TARGET)  = tri(itri)%sideindex(2,3)
        newsrf%index(IND_SURFACE) = tri(itri)%sur(3)
        newsrf%type    = SPR_LINE_SEGMENT
        newsrf%link    = 0
        newsrf%nvtx    = 2
        newsrf%ivtx(1) = tri(itri)%ver(3)
        newsrf%ivtx(2) = tri(itri)%ver(1)
        ntrysrf = ntrysrf + 1
        trysrf(ntrysrf) = newsrf
        try(itri)%iside(3) = ntrysrf
      ENDDO


      RETURN
 99   STOP
      END
c
c
c Read triangle
c Move data to regular grid
c
c ======================================================================
c
c subroutine: ProcessTetrahedrons_06
c
      SUBROUTINE ProcessTetrahedrons_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      USE mod_eirene06_locals
      USE mod_geometry
      IMPLICIT none

      TYPE(type_srf   ) newsrf
      TYPE(type_object) newobj


c...  Bricks:

      INTEGER itry,nseg,i1,i2,isid,itet,ivtx,save_nobj,save_isrf,
     .        iobj,isrf,istart,iend,isector,idum1
      LOGICAL regular_grid
      LOGICAL :: output = .FALSE.
      LOGICAL :: hack   = .FALSE.

      REAL*8 a(3,3),b(3,7),c(3,4),ang,ang1,ang2,dang(500,2),
     .       theta,frac

      REAL*8     DTOL
      PARAMETER (DTOL=1.0D-07)

C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER s(4,14)
C     INTEGER s(4,14) /3, 2, 1, 7,   
      DATA    s       /3, 2, 1, 7,   
     .                 1, 2, 5, 7,   5, 4, 1, 7,   
     .                 2, 3, 6, 7,   6, 5, 2, 7,   
     .                 3, 1, 4, 7,   4, 6, 3, 7,   
     .                 4, 5, 6, 7,
!...  Reversed:
     .                 1, 2, 4, 7,   4, 2, 5, 7,
     .                 2, 3, 5, 7,   5, 3, 6, 7,
     .                 3, 1, 6, 7,   6, 1, 4, 7 /

      INTEGER ishift,iside,itry1
      LOGICAL, ALLOCATABLE :: trycheck(:)

C     Krieger IPP/07 - SUN compiler chokes on this syntax, wants
C     variable declaration and initialization separately
      INTEGER t(3,4)
C     INTEGER t(3,4) /1, 2, 3,  1, 4, 2,  2, 4, 3,  3, 4, 1 /
      DATA    t      /1, 2, 3,  1, 4, 2,  2, 4, 3,  3, 4, 1 /

      WRITE(eirfp,*) 'BUILDING TETRAHEDRONS',output,eirfp

      IF (.TRUE.) THEN
c...    Convert legacy triangle objects to generalized geometry objects:
        ntry = ntri
        ntrysrf = 0
        ntryvtx = 0
        ALLOCATE(try(ntry))
        ALLOCATE(trysrf(3*ntry))
        ALLOCATE(tryvtx(3,6*ntry))
        ALLOCATE(plasma(20,ntry))
        ALLOCATE(bfield(3 ,ntry))
        CALL BuildNewTriangleObjects
      ENDIF

      WRITE(eirfp,*) '  NTRY:',ntry
      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Start of tetrahedron objects:
      istart = nobj + 1
      
c...  Build bricks:
      ngrp = 12
      grp(3)%origin = GRP_MAGNETIC_GRID  ! *** NEED AddGroup and CopyGroup functions ***
      grp(3)%type   = GRP_TETRAHEDRON
      grp(4)%origin = GRP_VACUUM_GRID
      grp(4)%type   = GRP_TETRAHEDRON


      regular_grid = .TRUE.


      IF (regular_grid) THEN   
c...    Toroidal distribution:
        ang = 360.0D0 / DBLE(ntorseg) * torfrac
        nseg = ntorseg
        dang(1,1) = 0.0D0
        dang(1,2) = dang(1,1) + ang
        DO isector = 2, nseg
          dang(isector,1) = dang(isector-1,2)
          dang(isector,2) = dang(isector  ,1) + ang
        ENDDO
      ELSE
        WRITE(0,*)
        WRITE(0,*) '=================================================='
        WRITE(0,*) '   WARNING: PARTICLE BALANCE ONLY SET FOR'
        WRITE(0,*) '            FOR REGULAR TOROIDAL GRID'
        WRITE(0,*) '=================================================='
        WRITE(0,*)
c        nseg = 12
c        dang(1 ,1) =  0.0D0    !  0.0D0       
c        dang(2 ,1) = 15.0D0    !  30.0D0      
c        dang(3 ,1) = 20.0D0    !  60.0D0      
c        dang(4 ,1) = 45.0D0    !  65.0D0      
c        dang(5 ,1) = 60.0D0    !  80.0D0      
c        dang(6 ,1) = 85.0D0    !  85.0D0      
c        dang(7 ,1) = 90.0D0    !  90.0D0      
c        dang(8 ,1) = 95.0D0    !  95.0D0      
c        dang(9 ,1) = 120.0D0   !  100.0D0     
c        dang(10,1) = 135.0D0   !  105.0D0     
c        dang(11,1) = 150.0D0   !  120.0D0     
c        dang(12,1) = 175.0D0   !  150.0D0     
c        dang(12,2) = 180.0D0   !  180.0D0     
        nseg = 12
        dang(1 ,1) = -30.0D0   !  0.0D0       
        dang(2 ,1) = -25.0D0   !  30.0D0      
        dang(3 ,1) = -20.0D0   !  60.0D0      
        dang(4 ,1) = -15.0D0   !  65.0D0      
        dang(5 ,1) = -10.0D0   !  80.0D0      
        dang(6 ,1) = -5.0D0    !  85.0D0      
        dang(7 ,1) =  0.0D0    !  90.0D0      
        dang(8 ,1) =  5.0D0    !  95.0D0      
        dang(9 ,1) =  10.0D0   !  100.0D0     
        dang(10,1) =  15.0D0   !  105.0D0     
        dang(11,1) =  20.0D0   !  120.0D0     
        dang(12,1) =  25.0D0   !  150.0D0     
        dang(12,2) =  30.0D0   !  180.0D0     
        DO isector = 1, nseg-1
          dang(isector,2) = dang(isector+1,1)
        ENDDO
      ENDIF

      DO isector = 1, nseg
        WRITE(0,'(A,I6,2F12.6)') 
     .    'ANGLES:',isector,dang(isector,1:2)
      ENDDO

      dang = dang * D_DEGRAD

c      ang1 = 360.0D0 / DBLE(ntorseg) * D_DEGRAD

      ALLOCATE(trycheck(ntry))
      trycheck = .FALSE.

      DO itry = 1, ntry  ! ntry   

c        IF (try(itry)%index(IND_IK).NE.1.OR.
c     .      try(itry)%index(IND_IR).NE.5) CYCLE  ! separatrix on sonnet_13018_250.sm
c        IF (hack) CYCLE
c        hack = .TRUE.

c...    Assemble brick vertices:
        DO i1 = 1, 3
          isrf = try(itry)%iside(i1)
          i2 = 1
          IF (isrf.LT.0) i2 = 2 ! Side orientation is switched, so use other end point
          a(1,i1) = tryvtx(1,trysrf(ABS(isrf))%ivtx(i2))
          a(2,i1) = tryvtx(2,trysrf(ABS(isrf))%ivtx(i2))
          a(3,i1) = 0.0D0
        ENDDO

c...    Check orientation...?
        IF (output) THEN
          WRITE(eirfp,*) 'A:',a(1:2,1)
          WRITE(eirfp,*) 'A:',a(1:2,2)
          WRITE(eirfp,*) 'A:',a(1:2,3)
        ENDIF

c...    Expand toroidally:
c        b(1,1:3) = a(1,1:3)
c        b(2,1:3) = a(2,1:3)
c        b(3,1:3) = a(1,1:3) * DTAN(-0.5D0*ang1)
c        b(1,4:6) = a(1,1:3)
c        b(2,4:6) = a(2,1:3)
c        b(3,4:6) = a(1,1:3) * DTAN(+0.5D0*ang1)

        isector = 1
        IF (output) THEN
          WRITE(eirfp,*) 'DANG:',dang(isector,1:2)/D_DEGRAD
        ENDIF
        b(1,1:3) = a(1,1:3) * DCOS(dang(isector,1))
        b(2,1:3) = a(2,1:3)
        b(3,1:3) = a(1,1:3) * DSIN(dang(isector,1))
        b(1,4:6) = a(1,1:3) * DCOS(dang(isector,2))
        b(2,4:6) = a(2,1:3)
        b(3,4:6) = a(1,1:3) * DSIN(dang(isector,2))

c        b(1,1:3) = a(1,1:3) * DCOS(-dang(isector,1))
c        b(2,1:3) = a(2,1:3)
c        b(3,1:3) = a(1,1:3) * DSIN(-dang(isector,1))
c        b(1,4:6) = a(1,1:3) * DCOS(+dang(isector,2))
c        b(2,4:6) = a(2,1:3)
c        b(3,4:6) = a(1,1:3) * DSIN(+dang(isector,2))

c       Center of brick:
        b(1:3,7) = 0.0D0
        DO i1 = 1, 6
          b(1:3,7) = b(1:3,7) + b(1:3,i1) / 6.0D0
c        DO i1 = 1, 3                               ! Bug, sort of...
c          b(1:2,7) = b(1:2,7) + b(1:2,i1) / 3.0D0  
        ENDDO

c...    Check orientation...?
        IF (output) THEN
          WRITE(eirfp,*) 'B:',b(1:3,1)
          WRITE(eirfp,*) 'B:',b(1:3,2)
          WRITE(eirfp,*) 'B:',b(1:3,3)
          WRITE(eirfp,*) 'B:',b(1:3,4)
          WRITE(eirfp,*) 'B:',b(1:3,5)
          WRITE(eirfp,*) 'B:',b(1:3,6)
          WRITE(eirfp,*) 'B:',b(1:3,7)
        ENDIF


c...    Determine group assignment (ad hoc at the moment):
c...    Make tetrahedrons:
        DO itet = 1, 8  ! 8 tetrahedrons for each triangle
          newobj = try(itry)
          newobj%segment(1) = itet
          newobj%group = try(itry)%group + 2 
          newobj%phi = SNGL(0.5D0*(dang(1,1) + dang(1,2)) / D_DEGRAD)
          newobj%nside = 4
          newobj%index(IND_IS) = 1
          newobj%index(IND_IK) = try(itry)%index(IND_IK)
          newobj%index(IND_IR) = try(itry)%index(IND_IR)

          ishift = 0
          IF (itet.GE.2.AND.itet.LE.7) THEN
            iside = INT(REAL(itet) / 2.0)
            itry1 = try(itry)%omap(iside)
c            itry1 = INT(try(itry)%map(iside))
            IF (itry1.NE.0) THEN
              IF (trycheck(itry1)) ishift = 7
            ENDIF
c            WRITE(eirfp,'(A,4I6,2I6)') 'TRYCHECK:',
c     .        itry,itet,iside,ishift,
c     .        try(itry)%smap(iside),itry1
          ENDIF

          c(1:3,1) = b(1:3,s(1,itet+ishift))
          c(1:3,2) = b(1:3,s(2,itet+ishift))
          c(1:3,3) = b(1:3,s(3,itet+ishift))
          c(1:3,4) = b(1:3,s(4,itet+ishift))

          IF (output) THEN
            WRITE(eirfp,*)
            WRITE(eirfp,*) itet,ishift
            WRITE(eirfp,*) s(1:4,itet+ishift)
            WRITE(eirfp,*) 'C:',c(1:3,1)
            WRITE(eirfp,*) 'C:',c(1:3,2)
            WRITE(eirfp,*) 'C:',c(1:3,3)
            WRITE(eirfp,*) 'C:',c(1:3,4)
          ENDIF

c...      Assign index mapping:
          DO isid = 1, newobj%nside
c           Vertices:
            IF (isid.EQ.1.AND.itet.GE.2.AND.itet.LE.7) THEN
              i1 = INT(REAL(itet-1)/2.0+0.51)                
              IF (output) WRITE(eirfp,*) 'IIII:',itet,i1
              newsrf%index = trysrf(ABS(try(itry)%iside(i1)))%index
          
              IF (output) THEN
                WRITE(eirfp,*) ' :',
     .            trysrf(ABS(try(itry)%iside(i1)))%index(1:3)
              ENDIF
            ELSE
              newsrf%index = 0
            ENDIF
            newsrf%type = SPR_PLANAR_POLYGON
            newsrf%obj  = nobj + 1
            newsrf%side = isid
            newsrf%nvtx = 3
            DO ivtx = 1, 3
c              IF (DABS(c(1,t(ivtx,isid))).LT.1.0D-07) THEN
c                WRITE(0,*) 'TETRAHEDRONS: SUSPICIOUS X-VAL A',nobj+1
c                STOP 
c              ENDIF
              newsrf%ivtx(ivtx) = AddVertex(c(1,t(ivtx,isid)))
            ENDDO
            
            newobj%iside(isid) = AddSurface(newsrf)

            IF (output) THEN
              WRITE(eirfp,*) 'I.:',isid
              WRITE(eirfp,*) 'IN:',newsrf%index(1:3)
              WRITE(eirfp,*) 'I0:',t(1:3,isid)
              WRITE(eirfp,*) 'I1:',newsrf%ivtx(1:3)
              WRITE(eirfp,*) 'I1:',srf(ABS(newobj%iside(isid)))%
     .                                     ivtx(1:3)

              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(1))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(2))
              WRITE(eirfp,*) 'V:',vtx(1:2,newsrf%ivtx(3))
            ENDIF

          ENDDO
c...      Center of tetrahedron -- needs to be done properly:
          newobj%x = SNGL(0.25D0 * SUM(c(1,1:4)))
          newobj%y = SNGL(0.25D0 * SUM(c(2,1:4)))
          newobj%z = SNGL(0.25D0 * SUM(c(3,1:4)))

c          IF (itry.EQ.1) THEN
c            WRITE(0,*) 'Cx:',c(1,1:4)
c            WRITE(0,*) 'Cy:',c(2,1:4)
c            WRITE(0,*) 'Cz:',c(3,1:4)
c            WRITE(0,*) 'Cen:',newobj%x,newobj%y,newobj%z
c          ENDIF

c...      Add object:
          idum1 = AddObject(newobj)

        ENDDO  ! Tetrahedron loop

        trycheck(itry) = .TRUE.
         
      ENDDO  ! Triangle loop

c...  Done with these, clear some memory:
      DEALLOCATE(trycheck)
      DEALLOCATE(try)
      DEALLOCATE(trysrf)
      DEALLOCATE(tryvtx)

      WRITE(eirfp,*) '  NTRY:',ntry
      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx
      WRITE(eirfp,*) '  TOROIDAL REPLICATION'

c...  Toroidal replication:
      save_nobj = nobj
c      DO ang = dang,  0.0D0*D_DEGRAD, dang
c      DO ang = 179.9D0*D_DEGRAD, 359.9D0*D_DEGRAD, dang
c      DO ang = dang, 179.9D0*D_DEGRAD, dang
c      DO ang = dang, 359.9D0*D_DEGRAD, dang
c      DO ang = DBLE(torus1)+dang, DBLE(torus2*0.99999)*D_DEGRAD, ang1

c      isector = 1   
c      DO ang = DBLE(torus1)+ang1, DBLE(torus2*0.99999)*D_DEGRAD, ang1
c        isector = isector + 1

      DO isector = 2, nseg ! nseg  
        ang1 = dang(isector,1) - dang(1,1)
        ang2 = dang(isector,2) - dang(1,2)
c        ang  = 0.5D0 * (ang1+ ang2) 
c        ang1 = ang - (dang(isector,1) - dang(1,1)) 
c        ang2 = ang + (dang(isector,2) - dang(1,2)) 

        WRITE(eirfp,'(A,I6,3F12.6)') '   TOROIDALIZING:',isector,
     .    SNGL(0.5D0 * (dang(isector,1)+dang(isector,2)) / D_DEGRAD),
     .    ang1/D_DEGRAD,ang2/D_DEGRAD

        DO iobj = 1, save_nobj ! save_nobj
          newobj = obj(iobj)
          newobj%index(IND_IS) = isector
          newobj%phi = 
     .      SNGL(0.5D0*(dang(isector,1) + dang(isector,2)) / D_DEGRAD)
c          newobj%phi = SNGL(ang) / D_DEGRAD 

          DO isid = 1, newobj%nside
            isrf = ABS(newobj%iside(isid))  ! I don't need to do any trickery here to make sure
            save_isrf = 0                   ! that new sides are CCW from the outside of the object
                                            ! because the properly oriented side will be duplicated
            DO WHILE (.TRUE.)               ! first and the subsequent CW sides will then be mapped 
                                            ! to those sides...
              newsrf = srf(isrf)
              newsrf%obj  = nobj + 1
              newsrf%side = isid
              DO i1 = 1, newsrf%nvtx
c...            Find toroidal angle of vertex relative to the minimum
c               toroidal angle in sector 1, since these tetrahedrons
c               may be stretched

                 

                a(1:3,1) = vtx(1:3,newsrf%ivtx(i1))

                IF (DABS(a(1,1)).LT.DTOL) THEN
c                IF (a(1,1).EQ.0.0D0) THEN
                  ang = 0.0D0
                ELSE
                  ang = DATAN(a(3,1) / a(1,1))
                ENDIF
          
c                WRITE(0,'(A,I9,I4,12X,4F12.4)') 
c     .            'ANG:',iobj,i1,ang/D_DEGRAD,a(1:3,1)

                frac = (ang - dang(1,1)) / (dang(1,2) - dang(1,1))

                ang = (1.0D0 - frac) * (dang(isector,1) - dang(1,1)) + 
     .                         frac  * (dang(isector,2) - dang(1,2))

c                WRITE(0,'(A,I9,I4,2F12.4)') 
c     .            '   :',iobj,i1,frac,ang/D_DEGRAD

                a(1,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
                a(2,2) = a(2,1)
                a(3,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
c                IF (DABS(a(1,2)).LT.1.0D-07) THEN
c                  WRITE(0,*) 'TETRAHEDRONS: SUSPICIOUS X-VAL B',nobj+1
c                  STOP 
c                ENDIF
                newsrf%ivtx(i1) = AddVertex(a(1,2))

c                WRITE(0,'(A,I9,I4,24X,3F12.4)') 
c     .            '   :',iobj,i1,a(1:3,1) 

c                a(1:3,1) = vtx(1:3,newsrf%ivtx(i1))
c                a(1  ,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
c                a(2  ,2) = a(2,1)
c                a(3  ,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
c                newsrf%ivtx(i1) = AddVertex(a(1,2))
              ENDDO


              IF (save_isrf.EQ.0) THEN 
                newobj%iside(isid) = AddSurface(newsrf)
              ELSE
                isrf = AddSurface(newsrf)
                srf(save_isrf)%link = isrf
              ENDIF

              IF (srf(isrf)%link.EQ.0) THEN
                EXIT
              ELSE
                save_isrf = isrf
                isrf = srf(isrf)%link 
                WRITE(eirfp,*) 'ACTUALLY USING THIS CODE?'
                STOP
              ENDIF

            ENDDO

          ENDDO

c...      Barf, need to tidy this up since it's a repeat of the above code:
          a(1,1) = DBLE(obj(iobj)%x)
          a(2,1) = DBLE(obj(iobj)%y)
          a(3,1) = DBLE(obj(iobj)%z)
          IF (DABS(a(1,1)).LT.DTOL) THEN
            ang = 0.0D0
          ELSE
            ang = DATAN(a(3,1)/a(1,1))
          ENDIF
          frac = (ang - dang(1,1)) / (dang(1,2) - dang(1,1))
          ang = (1.0D0 - frac) * (dang(isector,1) - dang(1,1)) + 
     .                   frac  * (dang(isector,2) - dang(1,2))
          a(1,2) = DCOS(ang) * a(1,1) - DSIN(ang) * a(3,1)
          a(2,2) = a(2,1)
          a(3,2) = DSIN(ang) * a(1,1) + DCOS(ang) * a(3,1)
          newobj%x = SNGL(a(1,2))
          newobj%y = SNGL(a(2,2))
          newobj%z = SNGL(a(3,2))
          IF (iobj.EQ.1) THEN
            WRITE(0,*) 'Cen:',newobj%x,newobj%y,newobj%z
          ENDIF

          idum1 = AddObject(newobj)

c          if (iobj.Eq.10)    STOP 'sdfsdfd'

        ENDDO
      ENDDO

 

C     IN RAY/OUT, I DON'T YET HAVE A SITUATION WHERE ORIENTATION IS A BIG ISSUE, SO FOR EIRENE
c     I'LL JUST DO THE ORDERING ON THE FLY, CHECKING IF THE POINTS ARE CO- OR COUNTER CLOCKWISE
c     WITH RESPECT TO THE CENTER OF THE TETRAHEDRON, BUT NEED TO COME UP WITH A TEST FIRST, 
c     PERHAPS THE NORMAL WITH RESPECT TO THE CENTER VECTOR? 

      iend = nobj

      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx

c...  Build connection map:
      CALL BuildConnectionMap(istart,iend)

c.... For Eirene, need to make sure the tetrahedrons are arranged according to
c     convention:
      CALL FixTetrahedrons(istart,iend)

c...  Plasma association:

c...  Impose filament structures:
      IF (.TRUE.) THEN
        WRITE(eirfp,*) '  NOBJ:',nobj
        WRITE(eirfp,*) '  NSRF:',nsrf
        WRITE(eirfp,*) '  NVTX:',nvtx
        CALL ResolveFilament   
        iend = nobj
        CALL BuildConnectionMap(istart,iend)
        CALL FixTetrahedrons(istart,iend)
      ENDIF

      CALL CheckTetrahedronStructure

      WRITE(eirfp,*) '  NOBJ:',nobj
      WRITE(eirfp,*) '  NSRF:',nsrf
      WRITE(eirfp,*) '  NVTX:',nvtx
      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c
c ======================================================================
c
c subroutine: ProcessTriangles_06
c
      SUBROUTINE ProcessTriangles_06(mode)
      USE mod_eirene06_parameters
      USE mod_eirene06
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
  
      WRITE(eirfp,*) 'PROCESSING TRIANGLES'  

      output = .FALSE.

      IF (mode.EQ.-1) GOTO 10

      WRITE(eirfp,*) '  REMOVING DUPLICATE VERTICIES'

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
            WRITE(eirfp,*) 'KILLING TRIANLGE-FIX!:',i1
            DO i2 = i1, ntri-1
              tri(i2) = tri(i2+1)
            ENDDO
            ntri = ntri - 1
            EXIT
          ENDIF
        ENDDO
      ENDDO


 10   CONTINUE

      WRITE(eirfp,*) '  ASSIGNING REGIONS' 
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
      WRITE(eirfp,*) '  BUILDING REGION LISTS' 
      nregion = 0
      iregion = 0
      DO i1 =  1, ntri
        ix = xregion(i1)
        iy = yregion(i1)
        nregion(ix,iy)                = nregion(ix,iy) + 1
        iregion(ix,iy,nregion(ix,iy)) = i1
      ENDDO
      DO i1 = 10, 1, -1
        WRITE(eirfp,'(2X,10I5)') (nregion(i2,i1),i2=1,10)
      ENDDO



      WRITE(eirfp,*) '  BUILDING CONNECTION MAP',ntri

c...  Build connection map:
      DO i1 = 1, ntri
        tri(i1)%map(1:3) = 0
        tri(i1)%sid(1:3) = 0
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
c              WRITE(eirfp,'(A,8I6)') ' SEARCH:',
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


      WRITE(eirfp,*) '  MAPPING SIDES TO SURFACES' 

c...  Map triangles to surfaces:     
      DO i1 = 1, ntri
        tri(i1)%sur(1:3) = 0
        DO v1 = 1, 3
          v2 = v1 + 1
          IF (v1.EQ.3) v2 = 1          

          tri(i1)%sur(v1) = 0

          DO i2 = 1, nsurface            
            IF     (surface(i2)%type.EQ.NON_DEFAULT_STANDARD) THEN

              IF (tri(i1)%type.EQ.MAGNETIC_GRID) THEN

                knot   = tri(i1)%index(1)  ! Parameter ...?
                ring   = tri(i1)%index(2)  ! Parameter ...?
                side   = tri(i1)%sideindex(1,v1)
                target = tri(i1)%sideindex(2,v1)

                IF     (surface(i2)%subtype .EQ.STRATUM.AND.
     .                  surface(i2)%index(1).LE.ring   .AND.
     .                  surface(i2)%index(2).GE.ring   .AND.
     .                  surface(i2)%index(3).EQ.target) THEN
                  tri(i1)%map(v1) = 0 ! This should only be set to 0 if the target is opaque...
                  tri(i1)%sid(v1) = 0 ! ditto
                  tri(i1)%sur(v1) = surface(i2)%num

                ELSEIF (surface(i2)%subtype.EQ.
     .                  MAGNETIC_GRID_BOUNDARY) THEN
c                  i3 = 1
c                  DO WHILE (surface(i2)%index(i3).NE.0)
                    IF (surface(i2)%index(1).LE.knot.AND.
     .                  surface(i2)%index(2).GE.knot.AND.
     .                  surface(i2)%index(3).EQ.ring.AND.
     .                  surface(i2)%index(4).EQ.side) THEN
                      tri(i1)%sur(v1) = surface(i2)%num 
                      EXIT
                    ENDIF
c                    i3 = i3 + 2
c                  ENDDO

                ELSEIF (surface(i2)%subtype .EQ.ADDITIONAL) THEN ! *** WRONG PLACE ***?
                ENDIF

              ENDIF

            ELSEIF (surface(i2)%type.EQ.VESSEL_WALL) THEN
              test = .TRUE.  ! ***REMOVE***
c...          Assign surface end points:
              x(0) = DBLE(surface(i2)%v(1,1))
              y(0) = DBLE(surface(i2)%v(2,1))
              x(1) = DBLE(surface(i2)%v(1,2))
              y(1) = DBLE(surface(i2)%v(2,2))
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
                  IF (surface(i2)%index(3).NE.0) THEN
                    tri(i1)%sur(v1) = surface(i2)%index(3)
                    tri(i1)%sideindex(3,v1)=surface(i2)%index(1)  ! Store xVESM index of surface
                    tri(i1)%sideindex(4,v1)=surface(i2)%index(2)  ! Store additional surface index 
c                    tri(i1)%sur(v1) = surface(surface(i2)%index(3))%num
                  ELSE
                    CALL ER('...','Surface index not assigned?1',*99)
c                    tri(i1)%sur(v1) = 4 ! Temp
                  ENDIF
                ENDIF
              ENDIF

            ELSE
              CALL ER('ProcessTriangles_06','Invalid surface type',*99)
            ENDIF
          ENDDO
        ENDDO
      ENDDO


      WRITE(eirfp,*) '  CHECKING ASSIGNMENTS' 

      DO i1 = 1, ntri
        DO v1 = 1, 3                             ! Perhaps eliminate maps through solid surfaces? 
          IF (tri(i1)%map(v1).EQ.0) THEN
            IF ( tri(i1)%sur(v1).EQ.0.OR.
     .          (tri(i1)%sur(v1).NE.0.AND.
     .           surface(MAX(1,tri(i1)%sur(v1)))%iliin.EQ.-1)) THEN 
c            IF (tri(i1)%map(v1).EQ.-1) THEN 
c...          If this shows up again, it may be related to DTOL in PointOnLine:
              WRITE(eirfp,*) 'PROBLEMS WITH MAP',i1
              CALL WriteEireneTriangles
              CALL SaveTriangles_06
              CALL DumpGrid('PROBLEM #1 WITH TRIANGLE MAP')
            ENDIF
          ENDIF
          IF (tri(i1)%map(v1).EQ.-1) THEN 
c...        If this shows up again, it may be related to DTOL in PointOnLine:
            WRITE(eirfp,*) 'PROBLEMS WITH MAP',i1
            CALL WriteEireneTriangles
            CALL SaveTriangles_06
            CALL DumpGrid('PROBLEM #2 WITH TRIANGLE MAP')
          ENDIF
c...      Check if 2 close lying points got merged by accident:
          v2 = v2 + 1
          IF (v1.EQ.3) v2 = 1
          IF (tri(i1)%ver(v1).EQ.tri(i1)%ver(v2)) 
     .      CALL ER('ProcessTriangles_06','2 sided triangle',*98)

        ENDDO
      ENDDO

      WRITE(eirfp,*) 'DONE'

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
      SUBROUTINE ReadPolyFile_06(zone)
      USE mod_eirene06_parameters
      USE mod_eirene06
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
        tri(i1+ntri)%index = 0
        tri(i1+ntri)%sideindex = 0
        tri(i1+ntri)%plasma = 0.0
        tri(i1+ntri)%bfield = 0.0
        tri(i1+ntri)%efield = 0.0
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
      INTEGER FUNCTION NewEireneSurface_06(type)
      USE mod_eirene06_parameters
      USE mod_eirene06
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

      nsurface = nsurface + 1
      surface(nsurface)%type   = type
      surface(nsurface)%index  = 0
      surface(nsurface)%subtype = -1
      surface(nsurface)%num    = 0
      surface(nsurface)%iliin  = defiliin
      surface(nsurface)%ilside = defilside
      surface(nsurface)%ilswch = defilswch
      surface(nsurface)%iltor  = defiltor
      surface(nsurface)%ilcell = defilcell
      surface(nsurface)%ilcol  = defilcol
      surface(nsurface)%recyct = defrecyct
      surface(nsurface)%recycf = defrecycf

      IF     (type.EQ.VESSEL_WALL) THEN

        surface(nsurface)%surtxt   = '* vessel wall (default)'
        surface(nsurface)%reflect = LOCAL
        surface(nsurface)%ewall = -wtemp * 1.38E-23 / ECH
        surface(nsurface)%material = wmater
c...    Assume a 2-point line segment:
        surface(nsurface)%nsur = 1
        surface(nsurface)%npts(1) = 2
        surface(nsurface)%ipts(1,1) = 1
        surface(nsurface)%ipts(2,1) = 2
        surface(nsurface)%nver = 2

      ELSEIF (type.EQ.NON_DEFAULT_STANDARD) THEN

        surface(nsurface)%surtxt  = '* non-default standard (default)'
        surface(nsurface)%reflect = GLOBAL
        surface(nsurface)%ewall    = 0.0
        surface(nsurface)%material = 0.0
        surface(nsurface)%nsur = 0
        surface(nsurface)%nver = 0

      ELSE
        CALL ER('NewEireneSurface_06','Invalid type',*99)
      ENDIF

      NewEireneSurface_06 = nsurface

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
      SUBROUTINE RefineTriangles_06
      USE mod_eirene06_parameters
      USE mod_eirene06
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

          WRITE(eirfp,*) 'I1,AREA=',i1,ntri,area

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
              WRITE(eirfp,*) '.......=',i1,i2,i3


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

            WRITE(eirfp,*) 'SPLITTING:',i1,iside
            WRITE(eirfp,*) 'SPLITTING:',itri2,iside2

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
        CALL ProcessTriangles_06(-1)

      ENDDO

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WritePolyFile
c
      SUBROUTINE WritePolyFile_06(eirntri,MAXNRS,eirtri)
      USE mod_eirene06_parameters
      USE mod_eirene06
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

      WRITE(eirfp,*) 'HERE IN WRITEPOLYFILE'

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

c          WRITE(eirfp,*) 'DATA:',i1,i2,eirtri(i2,2)

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

c            WRITE(eirfp,*) 'I!:',i1

            DO i3 = 1, nsurface
              IF (surface(i3)%type.EQ.VESSEL_WALL.AND.
     .            ((eirtri(i2,2).EQ.2.0.AND.
     .              surface(i3)%index(1).GE.NINT(eirtri(i2,3)).AND.
     .              surface(i3)%index(1).LE.NINT(eirtri(i2,4))).OR.
     .             (eirtri(i2,2).EQ.3.0.AND.
     .              surface(i3)%index(2).GE.NINT(eirtri(i2,3)).AND.
     .              surface(i3)%index(2).LE.NINT(eirtri(i2,4))))) THEN

                x1 = DBLE(surface(i3)%v(1,1))
                y1 = DBLE(surface(i3)%v(2,1))
                x2 = DBLE(surface(i3)%v(1,2))
                y2 = DBLE(surface(i3)%v(2,2))

                len = DSQRT((x1 - x2)**2 + (y1 - y2)**2)

                IF (eirtri(i1,3).GT.0.0.AND.
     .              len.GT.DBLE(eirtri(i1,3))) THEN
                  tstep = 1.0D0 / DBLE(INT(len/DBLE(eirtri(i1,3))) + 1)
                ELSE
                  tstep = 1.0D0
                ENDIF

c                WRITE(eirfp,*) x1,y1
c                WRITE(eirfp,*) x2,y2
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
        WRITE(eirfp,*) 'COMMAND: >'//command(1:LEN_TRIM(command))//'<'
        CALL CIssue(command(1:LEN_TRIM(command)),code)
        WRITE(eirfp,*) 'RETURN_CODE:',code

        icnt = icnt + 1
c        IF (icnt.EQ.1) STOP 'STOP: CHECK POLY FILES'

        CALL ReadPolyFile_06(NINT(eirtri(i1,2)))


      ENDDO

c      STOP 'TEST'

      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteTriangleFiles
c
      SUBROUTINE WriteEireneTriangles
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2,v1,ik1,ir1,it,idum1,itri
      LOGICAL found
      REAL    version,rdum1

      REAL sumflux1,sumflux2

      REAL, ALLOCATABLE :: tdata(:)      

      WRITE(eirfp,*) 'WRITING TRIANGLE FILES'

      version = 1.00

      fp = 99

      ALLOCATE(tdata(ntri))

      IF (photons.EQ.-1) THEN
c...    Load ionisation data from previous EIRENE call:
        CALL LoadTriangleData(7,0,13,0,tdata)  
      ELSE
        tdata = -999.0
      ENDIF

c...  Dump triangles (for OUT, not EIRENE):
      OPEN(UNIT=fp,FILE='objects.dat',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
        WRITE(fp,'(I6,3(2F10.6,2X))')
     .    i1,(ver(tri(i1)%ver(i2),1),ver(tri(i1)%ver(i2),2),i2=1,3)
      ENDDO
      CLOSE(fp)      


c...  Dump vertices:
      OPEN(UNIT=fp,FILE='objects.npco_char',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
      WRITE(fp,*) nver
      DO i1 = 1, nver
        WRITE(fp,'(I6,3F12.6)') i1,ver(i1,1)*100.0,ver(i1,2)*100.0,0.0
      ENDDO
      CLOSE(fp)      

c...  Dump sides:
      OPEN(UNIT=fp,FILE='objects.elemente',ACCESS='SEQUENTIAL',
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
      OPEN(UNIT=fp,FILE='objects.neighbors',ACCESS='SEQUENTIAL',
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

c...  Loading magnetic field data from idl/magnetics/b.pro for "field everywhere"
c     grids for Detlev:
      IF (.FALSE..AND..NOT.tetrahedrons) THEN
        OPEN(UNIT=fp ,FILE='objects.bfield',ACCESS='SEQUENTIAL',
     .       STATUS='OLD',ERR=95)      
        READ(fp,*)
        READ(fp,*)
        DO itri = 1, ntri
          READ(fp,*) (rdum1,i1=1,4),idum1,idum1,(rdum1,i1=1,3),
     .               tri(itri)%bfield(1),tri(itri)%bfield(3),
     .               tri(itri)%bfield(2)
        ENDDO
        CLOSE(fp)
      ENDIF

      OPEN(UNIT=fp ,FILE='objects.plasma',ACCESS='SEQUENTIAL',
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
c         on the fluid code cell/ring indices:
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
     .          ik1,ir1                 ! Fluid grid indices, for debugging only

            ELSE
              CALL ER('WriteTriangleFiles','Target data appears to  '//
     .                'be over-specified',*99)
            ENDIF
          ENDDO

        ENDDO
      ENDDO
      CLOSE(fp)      

c      WRITE(eirfp,*) 'SUMFLUX:',sumflux1,sumflux2

      OPEN(UNIT=fp,FILE='objects.efield',ACCESS='SEQUENTIAL',
     .     STATUS='REPLACE',ERR=96)      
c...  Header:
      WRITE(fp,'(A)') '* VERSION 1.0 OSM E-FIELD FILE FOR TRIANGULAR '//
     .                'GRID'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A)') '* BULK PLASMA DATA (PART DEUX)'
      WRITE(fp,'(A)') '*'
      WRITE(fp,'(A7,4A12)') 
     .  '* Index','Ex','Ey','Ez','E_pot'
      WRITE(fp,'(A7,4A12)')
     .  '*      ','V m-1','V m-1','V m-1','V'
      WRITE(fp,*) ntri
      DO i1 = 1, ntri
c...    Dump efield data:
        WRITE(fp,'(I7,1P,4E12.4,10X,0P,2I4)') i1,
     .    tri(i1)%efield(1),
     .    tri(i1)%efield(2),
     .    tri(i1)%efield(3),
     .    tri(i1)%plasma(20),  ! temporary
     .    tri(i1)%index(1),tri(i1)%index(2) 
      ENDDO
      CLOSE(fp)      

      WRITE(eirfp,*) 'DONE'

      IF (ALLOCATED(tdata)) DEALLOCATE(tdata)

      RETURN
95    WRITE(0,*) 'WRITETRIANGEFILES: B-FIELD FILE NOT FOUND'
      STOP
96    WRITE(0,*) 'WRITETRIANGEFILES: PROBLEMS WITH FILE ACCESS'
      STOP
99    STOP
      END

c
c ======================================================================
c
c subroutine: NewTriangle
c
      SUBROUTINE AssignPlasmaQuantities_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER i1,i2

      WRITE(eirfp,*) 'ASSIGNING PLASMA QUANTITIES'

c      DO i1 = 1, ntri
c        tri(i1)%plasma = 0.0
c        tri(i1)%bfield = 0.0
c        tri(i1)%efield = 0.0
c      ENDDO

      DO i1 = 1, ntri
        DO i2 = 1, ncell
          IF (tri(i1)%type.EQ.MAGNETIC_GRID.AND.
     .        tri(i1)%index(1).EQ.cell(i2)%index(1).AND.
     .        tri(i1)%index(2).EQ.cell(i2)%index(2)) THEN     
            tri(i1)%plasma(1:6) = cell(i2)%plasma(1:6) 
            tri(i1)%bfield(1:3) = cell(i2)%bfield(1:3) 
            tri(i1)%efield(1:3) = cell(i2)%efield(1:3) 
            tri(i1)%plasma(20)  = cell(i2)%e_pot  
          ENDIF
        ENDDO
      ENDDO

      WRITE(eirfp,*) 'DONE'

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: NewTriangle
c
      SUBROUTINE NewTriangle_06(icell)
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER icell

      ntri = ntri + 1

      tri(ntri)%type = MAGNETIC_GRID 
      tri(ntri)%index(1) = cell(icell)%index(1)
      tri(ntri)%index(2) = cell(icell)%index(2)
      tri(ntri)%sideindex(1,1:3) = 0
      tri(ntri)%sideindex(2,1:3) = 0
      tri(ntri)%sideindex(3,1:3) = 0

      RETURN
 99   STOP
      END
c
c ======================================================================
c
c subroutine: AssignVertex
c
      SUBROUTINE AssignVertex_06(ivert,ipoint,iside,ilist,xlist,ylist)
      USE mod_eirene06_parameters
      USE mod_eirene06
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
      SUBROUTINE BuildFluidGridTriangles_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      REAL       TOL
      PARAMETER (TOL=1.0E-07)

c      REAL    GetMach,GetJsat,GetFlux 
      LOGICAL PointOnLine

      INTEGER i1,i2,i3,i4,nlist(2),ilist(50,2)
      LOGICAL output,status
      REAL*8  x(3),y(3),s,t,xlist(0:50,2),ylist(0:50,2),slist(0:50,2)

c      REAL*8, ALLOCATABLE :: xvertex(:),yvertex(:)

      WRITE(eirfp,*) 'BUILDING SUPER TRIANGLES'

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
c                WRITE(eirfp,*) 'A:',i2,i3
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
c              WRITE(eirfp,*) 'B:',i2,i3
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
c            WRITE(eirfp,'(A,I3,3F10.4)') 
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
c              WRITE(eirfp,*) 'ELIMINATING TRIANGLE POINT!', ! This should be unnecessary...
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
                WRITE(eirfp,*) 'SORTING TRIANGLE SIDE!',    ! This sorting should be unnecessary...
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
c            WRITE(eirfp,'(A,I3,3F10.4)') 
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

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,4)
            IF (i2.EQ.1) 
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,1)

            CALL AssignVertex_06(1,i2  ,2,ilist,xlist,ylist)
            CALL AssignVertex_06(2,i2  ,1,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2+1,1,ilist,xlist,ylist)

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            IF (i2.EQ.nlist(1)-1.AND.nlist(1).EQ.nlist(2)) 
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)

            CALL AssignVertex_06(1,i2+1,1,ilist,xlist,ylist)
            CALL AssignVertex_06(2,i2+1,2,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2  ,2,ilist,xlist,ylist)

          ELSEIF (nlist(1).GT.i2) THEN

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,4)
            IF (i2.EQ.nlist(1)-1) 
     .        tri(ntri)%sideindex(2,3) = cell(i1)%sideindex(2,3)

            CALL AssignVertex_06(1,nlist(2),2,ilist,xlist,ylist)
            CALL AssignVertex_06(2,i2      ,1,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2+1    ,1,ilist,xlist,ylist)

          ELSEIF (nlist(2).GT.i2) THEN

c            WRITE(eirfp,*) 'COOL B:',i2  

            CALL NewTriangle_06(i1)
            tri(ntri)%sideindex(1,2) = cell(i1)%sideindex(1,2)
            IF (i2.EQ.nlist(2)-1)
     .        tri(ntri)%sideindex(2,1) = cell(i1)%sideindex(2,3)

            CALL AssignVertex_06(1,nlist(1),1,ilist,xlist,ylist)
            CALL AssignVertex_06(2,i2+1    ,2,ilist,xlist,ylist)
            CALL AssignVertex_06(3,i2      ,2,ilist,xlist,ylist)

          ELSE
            CALL ER('BuildFluidGridTriangles','Unknown situation',*99)
          ENDIF

        ENDDO

      ENDDO


c...  Remove vertex duplicates:




c      CALL WriteEireneTriangles
c      CALL DumpGrid('PROBLEM WITH TRIANGLE MAP')
c      STOP 'CRAPPO!'





c      STOP 'DUMPING'
      WRITE(eirfp,*) 'DONE'

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
      subroutine SaveTriangles_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp,i1,i2
      REAL*8  cen(3)

      fp = 99
      OPEN(UNIT=fp,FILE='triangles.raw',ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.0,ntri,nver,nsurface
      WRITE(fp,ERR=98) (tri(i1),i1=1,ntri)
      WRITE(fp,ERR=98) ((ver(i1,i2),i2=1,3),i1=1,nver)
      WRITE(fp,ERR=98) (surface(i1),i1=1,nsurface)
      CLOSE (fp)


      IF (.TRUE.) THEN
c...    Dump grid data to an external file:
        OPEN (UNIT=fp,FILE='objects.centre',ACCESS='SEQUENTIAL',
     .        STATUS='REPLACE')      
        WRITE(fp,*) 'SHOT: 000000   TIME: 0000'
        WRITE(fp,'(2A6,2A10)') 'IK','IR','R (m)','Z (m)'
        DO i1 = 1, ntri
          CALL GetObjCentre(i1,cen)
          WRITE(fp,'(2I6,2F10.6)') 
     .      tri(i1)%index(1),tri(i1)%index(2),
     .      SNGL(cen(1)),SNGL(cen(2))
        ENDDO
        CLOSE(fp)
      ENDIF
      
      RETURN
 98   CALL ER('SaveTriangles','Problems writing data file',*99)
 99   STOP
      END
c
c ======================================================================
c
c  subroutine: CopyBlock
c
      SUBROUTINE CopyBlock(fp1,fp2)
      IMPLICIT none
      INTEGER fp1,fp2
      CHARACTER buffer*200

      DO WHILE (.TRUE.) 
        CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') THEN
          CALL WriteLine(fp2,buffer)
        ELSE
          BACKSPACE fp1
          EXIT
        ENDIF
      ENDDO

      RETURN
 97   CALL ER('CopyBlock','Unexpected end of file',*99)
 98   CALL ER('CopyBlock','Problems reading template file',*99)
 99   STOP
      END
c
c ======================================================================
c
c subroutine: WriteEireneInputFile_06
c
      SUBROUTINE WriteEireneInputFile_06
      USE mod_eirene06
      IMPLICIT none

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

      IF (output) WRITE(eirfp,*) 'WRITING EIRENE INPUT FILE 06'
c
c     Initialization:
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

c      CALL MS('WriteInputFile','Using xVESM to store wall data')

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


      ELSEIF (buffer(1:6).EQ.'*** 1.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        fp06 = fp2

c        IF (eirphoton.EQ.2) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,3)
c        ELSE
          CALL WriteBlock01_06(fp1,fp2)
22        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 22
          BACKSPACE fp1
c        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 2.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteBlock02_06

24      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 24
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 3a') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteBlock03a_06

26      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 26
        BACKSPACE fp1

      ELSEIF (buffer(1:6).EQ.'*** 3b'.OR.
     .        buffer(1:6).EQ.'*** 3B') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

c        CALL WriteLine(fp2,buffer)
c        CALL CopyBlock(fp1,fp2)
        CALL WriteBlock03b_06

25      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 25
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 4.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

c        IF (eirphoton.GT.0) THEN
c          CALL WriteLine(fp2,buffer)
c          CALL Transferline2(fp1,fp2,buffer,98)
c        ELSE

        IF (.FALSE..AND.tetrahedrons) THEN  ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock04_06
41        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 41
          BACKSPACE fp1
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 5.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (.FALSE..AND.tetrahedrons) THEN    ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock05_06
39        CALL ReadLine(fp1,buffer,1,*97,*98)
          IF (buffer(1:3).NE.'***') GOTO 39
          BACKSPACE fp1
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 6.') THEN 
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,14)
        ELSEIF (.FALSE..AND.tetrahedrons) THEN     ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock06_06(fp1,fp2)
        ENDIF

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 7.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          CALL WriteLine(fp2,buffer)
          CALL Transferline2(fp1,fp2,buffer,186)
c        ELSEIF (beam.NE.0) THEN
        ELSE
          CALL WriteBlock07_06
c        ELSE
c          CALL WriteLine(fp2,buffer)
c          CALL CopyBlock(fp1,fp2)
        ENDIF

c...  Advance input stream to the start of the next input block:
40      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 40
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:7).EQ.'*** 10.') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (photons.GT.0) THEN
          WRITE(eirfp,*) 'HERE IN PHOTON CODE!'
        ELSE
        ENDIF

        CALL WriteLine(fp2,buffer)
        CALL CopyBlock(fp1,fp2)

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 11') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        IF (.FALSE..AND.tetrahedrons) THEN   ! TETRAHEDRONS OFF
          CALL WriteLine(fp2,buffer)
          CALL CopyBlock(fp1,fp2)
        ELSE
          CALL WriteBlock11_06
        ENDIF

29      CALL ReadLine(fp1,buffer,1,*97,*98)
        IF (buffer(1:3).NE.'***') GOTO 29
        BACKSPACE fp1

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 12') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,1)

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 13') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

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

        IF (output) WRITE(eirfp,*) 'DONE'
      ELSEIF (buffer(1:6).EQ.'*** 14') THEN
        IF (output) WRITE(eirfp,*) buffer(1:6)

        CALL WriteLine(fp2,buffer)
        CALL Transferline2(fp1,fp2,buffer,16)
        IF (output) WRITE(eirfp,*) 'DONE'

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
      SUBROUTINE WriteBlock01_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER   ntime,iteration
      DATA iteration /0/
      SAVE

      iteration = iteration + 1

      WRITE(fp06,'(A,I6)') 
     .  '*** 1. DATA FOR OPERATING MODE (OSM), CALL ',iteration

      ntime = 0
      IF (dtimv.NE.0.0) ntime = 1

      nfile = 111
      IF (photons.EQ.2) nfile = 311
      IF (tetrahedrons) nfile = 0

      IF     (niter.GE.1) THEN
c...    BGK or photons:
        WRITE(fp06,91) 2,0,time,nfile,0,niter,0,ntime
        WRITE(fp06,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp06,90) 'FFFFF FFFF'
      ELSEIF (.TRUE.) THEN
c...    Standard (no BGK or photons):
        WRITE(fp06,91) 2,0,time,nfile,0,0,0,ntime
c        WRITE(fp06,91) 2,0,time,nfile,0,1,0,ntime
c        WRITE(fp06,91) 1,1,0,0,1,9,1,0,5  ! NGSTAL=1
        WRITE(fp06,91) 1,1,0,0,1,9,0,0,5  
        WRITE(fp06,90) 'FFFFF FFFF'
      ELSE
        CALL ER('WriteBlock01_06','Trouble',*99)
      ENDIF


      RETURN
90    FORMAT(A)
91    FORMAT(3I6,I6.5,20(I6:))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock02_06
c
c
c
c
      SUBROUTINE WriteBlock02_06
      USE mod_geometry
      USE mod_eirene06
      IMPLICIT none
 
      WRITE(fp06,90) '*** 2. DATA FOR STANDARD MESH (DIVIMP)'

      IF (tetrahedrons) THEN
        WRITE(fp06,91) 1,1,1,1
        WRITE(fp06,90) 'T'
        WRITE(fp06,90) 'FFFFF FTF'
        IF (nvtxmap.NE.0) THEN
          WRITE(fp06,94) nobj+1,0,0,0,0,nvtxmap
        ELSE
          WRITE(fp06,94) nobj+1,0,0,0,0,nvtx
        ENDIF
        WRITE(fp06,90) 'CASE cmod'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,0
        WRITE(fp06,92) 0.0,0.0,0.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,1
        WRITE(fp06,92) 0.0,0.0,1.0,1.0,1.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
      ELSEIF (.TRUE.) THEN
        WRITE(fp06,91) 1,1,1,1
        WRITE(fp06,90) 'T'
        WRITE(fp06,90) 'FFFFF TFF'
c        WRITE(fp06,91) ntri,0,0,0,nver
        WRITE(fp06,91) ntri+1,0,0,0,nver
        WRITE(fp06,90) 'CASE cmod'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'TFFF'
        WRITE(fp06,91) 1,0
        WRITE(fp06,92) 0.0,0.0,0.0
        WRITE(fp06,90) 'F'
c        WRITE(fp06,90) 'TFFF'  ! cylindrical (MAST)
c        WRITE(fp06,91) 0,0,0
c        WRITE(fp06,92) 0.0,0.0,3.768*100.0
        WRITE(fp06,90) 'FTFF'  ! toroidal 
        WRITE(fp06,91) 1,1,100
        WRITE(fp06,92) 0.0,0.0,360.0
        WRITE(fp06,90) 'F'
        WRITE(fp06,91) 0
        IF (.TRUE..OR.beam.EQ.1) THEN
          WRITE(fp06,90) 'T'  ! Additional cell for beam particles to launch into
          WRITE(fp06,91) 1    !  (PB set it up this way...)
          WRITE(fp06,93) 1.0E+00
        ELSE
          WRITE(fp06,90) 'F'
          WRITE(fp06,91) 0
        ENDIF
      ELSE
        CALL ER('WriteBlock02_06','Trouble',*99)
      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6))
94    FORMAT(20(I8))
92    FORMAT(1P,20(E12.4))
93    FORMAT(1P,20(E12.5))
99    STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock04_06
c
c
c
c
      SUBROUTINE WriteBlock04_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER ncra,ncrm,iscd1,iscd2,natmi

      WRITE(fp06,90) '*** 4. DATA FOR SPECIES SPECIFICATION AND '//
     .               'ATOMIC PHYSICS MODULE (OSM)'

      IF     (.TRUE.) THEN

        WRITE(fp06,90) '* ATOMIC REACTION CARDS  NREACI='
        WRITE(fp06,91) 35
        WRITE(fp06,95) '  1 AMJUEL H.4 2.1.5    EI ',0,1
        WRITE(fp06,95) '  2 CONST  H.4          EI ',0,1     ! Dummy ionisation rate
        WRITE(fp06,92) -200.0,0.0,0.0,0.0,0.0,0.0
        WRITE(fp06,92)    0.0,0.0,0.0					  
        WRITE(fp06,95) '  3 AMJUEL H.102.1.5    EI ',0, 1
        WRITE(fp06,95) '  4 HYDHEL H.1 3.1.8    CX ',1, 1
        WRITE(fp06,95) '  5 HYDHEL H.3 3.1.8    CX ',1, 1
c       WRITE(fp06,95) '  6 AMJUEL H.2 2.26B0   EI ',0,56  ! For iron...
        WRITE(fp06,95) '  6 METHAN H.2 2.23      EI',0,12
        WRITE(fp06,95) '  7 METHAN H.1 3.2       CX',1,12
        WRITE(fp06,95) '  7 METHAN H.3 3.2       CX',1,12
        IF (opacity.EQ.5.AND.photons.EQ.0) THEN
          WRITE(fp06,95) '  8 H2VIBR H.4 2.1.8a   RC ',0,1
        ELSE
          WRITE(fp06,95) '  8 AMJUEL H.4 2.1.8    RC ',0,1				  
        ENDIF
        WRITE(fp06,95) '  9 HYDHEL H.2 2.2.9    EI ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 10 HYDHEL H.2 2.2.5    DS ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 11 HYDHEL H.2 2.2.10   DS ',0, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 13 AMJUEL H.0 0.3T     EL ',1, 2				  
        WRITE(fp06,95) ' 13 AMJUEL H.1 0.3T     EL ',1, 2				  
        WRITE(fp06,95) ' 13 AMJUEL H.3 0.3T     EL ',1, 2,0.0,0.0,0.0
        WRITE(fp06,95) ' 14 AMJUEL H.4 2.2.12   EI ',0, 2
        WRITE(fp06,95) ' 15 AMJUEL H.4 2.2.11   EI ',0, 2
        WRITE(fp06,95) ' 16 AMJUEL H.4 2.2.14   EI ',0, 2
        WRITE(fp06,95) ' 17 AMJUEL H.8 2.2.14   EI ',0, 2
        WRITE(fp06,95) ' 18 AMJUEL H.3 3.2.3    CX ',1, 2
        IF (.TRUE..OR.photons.NE.0) THEN   
c          CALL ER('WriteBlock04_06','Need to check PHOTON setup',*99)
          WRITE(fp06,96) ' 19 PHOTON H.2 HLya121.5669a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 20 PHOTON H.2 HLyb102.5722a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 21 PHOTON H.2 HBaa656.4667a RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 22 PHOTON P.1 HLya121.5669a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 23 PHOTON P.1 HLyb102.5722a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 24 PHOTON P.1 HBaa656.4667a OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 25 PHOTON H.2 HLyg97.2536a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 26 PHOTON P.1 HLyg97.2536a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 27 PHOTON H.2 HLyd94.9742a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 28 PHOTON P.1 HLyd94.9742a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 29 PHOTON H.2 HLye93.7803a  RC ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
          WRITE(fp06,96) ' 30 PHOTON P.1 HLye93.7803a  OT ',0,1,0.,0.,0.
          WRITE(fp06,91) 3,2,1,0
        ENDIF
        IF (.TRUE..OR.bgk.EQ.3) THEN
c          CALL ER('WriteBlock04_06','Need to check BGK setup',*99)
          WRITE(fp06,90) ' 31 CONST  H.2           EL  2  2'				! 17 -> 31
          WRITE(fp06,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0 
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0 
          WRITE(fp06,90) ' 32 CONST  H.2           EL  2  2'	                        ! 18 -> 32
          WRITE(fp06,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
          WRITE(fp06,90) ' 33 CONST  H.2           EL  4  4'                              ! 19 -> 33
          WRITE(fp06,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
          WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
        ENDIF
        IF (.TRUE..OR.beam.EQ.1) THEN 
c          CALL ER('WriteBlock04_06','Need to check BEAM setup',*99)
          WRITE(fp06,96) ' 34 AMJUEL H.1 3.1.6FJ  PI ',1,1,0.0,0.0,0.0  ! For beams...
          WRITE(fp06,95) ' 35 AMJUEL H.102.1.8     RC',0,1,1.36E+01     ! Rec e cooling...
        ENDIF

        natmi = 1
        IF (beam.EQ.1) natmi = natmi + 1

        WRITE(fp06,90) '** 4a NEUTRAL ATOMS SPECIES CARDS: NATMI='
        WRITE(fp06,91) natmi
        ncra = 3
        ncrm = 5
        IF (bgk.EQ.3) THEN
          ncra = ncra + 2
          ncrm = ncrm + 2
        ENDIF

        WRITE(fp06,94) 1,'D(n=1)  ',2,1,1,0,1,-4,0,ncra			    
        WRITE(fp06,91) 1,115,114,0  ,30000
        WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        IF (photons.EQ.-1) THEN
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 999.0,0.0,0.0,0.0,1.0
        ELSE
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
        ENDIF
        WRITE(fp06,91) 4,114,111,114,01001				    
        WRITE(fp06,93) 0.0,0.0,0.0,0.0,1.0
        IF (bgk.EQ.3) THEN
          WRITE(fp06,91) 31,214,0,0,01001,0,111
          WRITE(fp06,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp06,91) 33,414,0,0,01001,0,112				    
          WRITE(fp06,93) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF
        IF (beam.EQ.1) THEN
          WRITE(fp06,94) 2,'D_B     ',2,1,1,0,1,-4,0,4			    
          WRITE(fp06,91) 1,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 2,115,114,0  ,30000
          WRITE(fp06,93) 2.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 4,114,111,114,01001				    
          WRITE(fp06,93) 0.0,0.0,0.0,0.0,1.0
          WRITE(fp06,91) 34,114,114,  0,01001 
          WRITE(fp06,93) 0.0,4.0,0.0
        ENDIF

        WRITE(fp06,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
        WRITE(fp06,91) 1							    
        WRITE(fp06,94) 1,'D2      ',4,2,2,0,1,1,0,ncrm,0,0
        WRITE(fp06,91)  9,115,113,0,0			    
        WRITE(fp06,93) -1.5400E+01,0.0,0.0,0.0,0.0
        WRITE(fp06,91) 10,115,121,0,0
        WRITE(fp06,93) -1.0500E+01,0.0,3.0,3.0,0.0
        WRITE(fp06,91) 11,115,111,114,0
        WRITE(fp06,93) -2.5000E+01,0.0,5.0,5.0,0.0
        WRITE(fp06,91) 13,114,  0,  0,01001				    
c        WRITE(fp06,92) 0.0,0.0,0.0,0.0,1.0E-10
        WRITE(fp06,92) 0.0,0.0,0.0,0.0
        WRITE(fp06,91) 18,114,111,113,01001				    
        WRITE(fp06,92) 0.0,0.0,0.0,0.0
        IF (bgk.EQ.3) THEN
          WRITE(fp06,91) 32,314,0,0,01001,0,112
          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
          WRITE(fp06,91) 33,514,0,0,01001,0,111
          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
        ENDIF

        WRITE(fp06,90) '**4c TEST ION SPECIES CARDS:  NIONI ION '//
     .                 'SPECIES ARE CONSIDERED, NIONI='
        WRITE(fp06,91) 1						
        WRITE(fp06,94) 1,'D2+     ',4,2,2,1,0,-1,0,3,-1		
        WRITE(fp06,91) 14,115,111,114,0
        WRITE(fp06,92) -1.0500E+01,0.0,4.3,4.3,0.0
        WRITE(fp06,91) 15,115,124,0,0
        WRITE(fp06,92) -1.5500E+01,0.0,0.25,0.25,0.0
        WRITE(fp06,91) 16,115,121,000,30000
        WRITE(fp06,92) 16.0,0.0,10.0,0.0,0.0

        IF (photons.EQ.1.OR.photons.EQ.2) THEN

          iscd1 = 214
          iscd2 = 0
          IF (bgk.EQ.3) THEN
            iscd1 = 614
            iscd2 = 4
          ENDIF

          WRITE(fp06,90) '** 4d photons'
          WRITE(fp06,91) 6			
          WRITE(fp06,94) 1,'Ba-alpha',2,1,0,0,1,+1,0,0
          WRITE(fp06,94) 2,'Ly-alpha',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 22,iscd1,( 3+iscd2)*100+14,0,0
c          WRITE(fp06,91) 22,iscd1,314,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 3,'Ly-beta ',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 23,iscd1,( 5+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 4,'Ly-gamma',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 26,iscd1,( 7+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 5,'Ly-delta',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 28,iscd1,( 9+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
          WRITE(fp06,94) 6,'Ly-epsi ',2,1,0,0,1,+1,0,1
          WRITE(fp06,91) 30,iscd1,(11+iscd2)*100+14,0,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0
        ELSE
          WRITE(fp06,90) '** 4d photons'
          WRITE(fp06,91) 0
        ENDIF

      ELSEIF (.FALSE.) THEN

c        WRITE(fp06,90) '* ATOMIC REACTION CARDS  NREACI='
c        WRITE(fp06,91) 22
c        WRITE(fp06,90) '  1 AMJUEL H.4 2.1.5     EI  0  1'
c        WRITE(fp06,90) '  2 AMJUEL H.102.1.5     EI  0  1'
c        WRITE(fp06,90) '  3 HYDHEL H.1 3.1.8     CX  1  1'
c        WRITE(fp06,90) '  3 HYDHEL H.3 3.1.8     CX  1  1'
c        WRITE(fp06,90) '  4 AMJUEL H.4 2.2.9     EI  0  2'
c        WRITE(fp06,90) '  5 AMJUEL H.4 2.2.5     DS  0  2'
c        WRITE(fp06,90) '  6 AMJUEL H.4 2.2.10    DS  0  2'
c        WRITE(fp06,90) '  7 AMJUEL H.4 2.2.12    DS  0  2'
c        WRITE(fp06,90) '  8 AMJUEL H.4 2.2.11    DS  0  2'
c        WRITE(fp06,90) '  9 AMJUEL H.4 2.2.14    DS  0  2'
c        WRITE(fp06,90) ' 10 AMJUEL H.8 2.2.14    DS  0  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.0 0.3T     EL   1  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.1 0.3T     EL   1  2'
c        WRITE(fp06,90) ' 12 AMJUEL H.3 0.3T     EL   1  2 0.00000E+00'//
c     .                 ' 0.00000E+00 0.00000E+00'
c        WRITE(fp06,90) ' 13 HYDHEL H.2 2.3.9     EI  0  4'
c        WRITE(fp06,90) ' 14 METHAN H.2 2.23      EI  0 12'
c        IF (.TRUE.) THEN
c          WRITE(fp06,90) ' 15 H2VIBR H.4 2.1.8a    RC  0  1'
c          WRITE(fp06,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ELSE
c          WRITE(fp06,90) ' 15 AMJUEL H.4 2.1.8     RC  0  1'
c          WRITE(fp06,90) ' 16 AMJUEL H.102.1.8     RC  0  1  1.3600E 01'
c        ENDIF
c        WRITE(fp06,90) ' 17 CONST  H.2           EL  2  2'				
c        WRITE(fp06,92) -2.1091E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 18 CONST  H.2           EL  2  2'				
c        WRITE(fp06,92) -2.0589E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 19 CONST  H.2           EL  4  4'
c        WRITE(fp06,92) -2.0357E+01,0.2500E+00,0.0,0.0,0.0,0.0
c        WRITE(fp06,92)  0.0       ,0.0       ,0.0,0.0,0.0,0.0
c        WRITE(fp06,90) ' 20 AMJUEL H.3 3.2.3     CX  1  2'
c        WRITE(fp06,90) ' 21 AMJUEL H.9 3.1.8     CX  1  1'
c        WRITE(fp06,90) ' 22 AMJUEL H.2 3.1.8FJ   CX  1  1'

c        WRITE(fp06,90) '*NEUTRAL ATOMS SPECIES CARDS: NATMI='
c        WRITE(fp06,91) 1							    
c        ncra = 2
c        ncrm = 5
c        IF (niter.GE.1) THEN
c          ncra = ncra + 2
c          ncrm = ncrm + 2
c        ENDIF

c        WRITE(fp06,94) 1,'D       ',2,1,1,0,1,-4,0,ncra			    
c        WRITE(fp06,91) 1,115,114,0  ,30000,000		    
c        WRITE(fp06,92) 2.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp06,91) 3,114,111,114,01001				    
c        WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        IF (niter.GE.1) THEN
c          WRITE(fp06,91) 17,214,0,0,01001,0,111
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp06,91) 19,414,0,0,01001,0,112				    
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF

c        WRITE(fp06,90) '** 4b NEUTRAL MOLECULES SPECIES CARDS: NMOLI='
c        WRITE(fp06,91) 1							    
c        WRITE(fp06,94) 1,'D2      ',4,2,2,0,0,2,0,ncrm
c        WRITE(fp06,91) 4,115,113,0				    
c        WRITE(fp06,92) -1.5400E+01,0.0000E+00				    
c        WRITE(fp06,91) 5,115,121,000				    
c        WRITE(fp06,92) -1.0500E+01,0.0000E+00,3.0000E+00,3.0000E+00	    
c        WRITE(fp06,91) 6,115,111,114				    
c        WRITE(fp06,92) -2.5000E+01,0.0000E+00,5.0000E+00,5.0000E+00	    
c        IF (niter.GE.1) THEN
c          WRITE(fp06,91) 18,314,0,0,01001,0,112
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c          WRITE(fp06,91) 19,514,0,0,01001,0,111
c          WRITE(fp06,92) 0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        ENDIF
c        WRITE(fp06,91) 20,114,111,113,01001				    
c        WRITE(fp06,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c        WRITE(fp06,91) 12,114,  0,  0,01001				    
c        WRITE(fp06,92)  0.0000E+00,0.0000E+00,0.0000E+00,0.0000E+00,1.0
c
c        WRITE(fp06,90) '** 4c TEST ION SPECIES CARDS:  NIONI ION '//
c     .                 'SPECIES ARE CONSIDERED, NIONI='
c        WRITE(fp06,91) 1						
c        WRITE(fp06,94) 1,'D2+     ',4,2,2,1,0,-4,0,3,-1		
c        WRITE(fp06,91) 7,115,111,114                        
c        WRITE(fp06,92) -1.0400E+01,0.0000E+00,4.3000E+00,4.3000E+00
c        WRITE(fp06,91) 8,115,124,000                        
c        WRITE(fp06,92) -1.5500E+01,0.0000E+00,0.2500E+00,0.2500E+00
c        WRITE(fp06,91) 9,115,121,000,30002                  
c        WRITE(fp06,92)  1.0000E+01,0.0000E+00,0.5000E+00,0.5000E+00

c        WRITE(fp06,90) '** 4d photons'
c        WRITE(fp06,91) 0

      ELSE
        CALL ER('WriteBlock04_06','Trouble',*99)
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
c subroutine: WriteBlock11_06
c
c
c
c
      SUBROUTINE WriteBlock11_06
      USE mod_eirene06
      IMPLICIT none

      WRITE(fp06,90) '*** 11. DATA FOR NUMERICAL/GRAPHICAL OUTPUT (OSM)'

      IF (.TRUE.) THEN
        WRITE(fp06,90) 'FTFFF FTFFF FFTFF TTTTT T'
c        WRITE(fp06,90) 'FTFFF TTFFF FFTFF TTTTT T'
c        WRITE(fp06,90) 'Ftttt ttttt tttt'
        WRITE(fp06,90) 'Ttttt ttttt tttt'
        WRITE(fp06,91) 4
        WRITE(fp06,91) 14,0
        WRITE(fp06,91) -2,0
        WRITE(fp06,91) -3,0
        WRITE(fp06,91) -4,0
        WRITE(fp06,91) 0
        WRITE(fp06,90) 'TTFTT FFTFT ftFFF F'
        WRITE(fp06,91) 1,ntri,1,1,1,1
        WRITE(fp06,90) 'F PEI                      1 001002'
        WRITE(fp06,90) 'F LPT                      1 003008'
        WRITE(fp06,90) 'F ENTRANCE AND COVER       1 009010'	  
        WRITE(fp06,90) 'F SOUFFLET                 2 033033 038038'
        WRITE(fp06,90) 'F VERTICAL PORT            2 039061 068068'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'F'
        WRITE(fp06,90) 'F'
        WRITE(fp06,92) 230.0,230.0, 80.0,0.0,-750.0
        WRITE(fp06,92)  95.0, 95.0,800.0,0.0,   0.0,750.0
        WRITE(fp06,92)  45.0, 20.0
        WRITE(fp06,91) 0,0,1,2,3,4,5,6,9,0,1
c        WRITE(fp06,91) 1,10,1,2,3,4,5,6,9,0,1
        WRITE(fp06,91) 0
      ELSE
        CALL ER('WriteBlock11_06','Trouble',*99)
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
c subroutine: WriteBlock03a_06
c
c
c
c
      SUBROUTINE WriteBlock03a_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none

      INTEGER i1,nstsi,instsi,def_ilcell,ilcell,n

c      WRITE(eirfp,*) 'WRITING BLOCK3a'

      def_ilcell = 0

c...  Count the number of non-default standard surfaces that have been
c     defined:
      nstsi = 0
      DO i1 = 1, nsurface
c        WRITE(eirfp,*) 'NSTSI=',nstsi 
        IF (surface(i1)%type.EQ.NON_DEFAULT_STANDARD) nstsi = nstsi + 1
      ENDDO

      WRITE(fp06,90) '*** 3a. DATA FOR NON-DEFAULT SURFACES (OSM)'
      WRITE(fp06,91) nstsi
      instsi = 0
      DO i1 = 1, nsurface
        IF (surface(i1)%type.NE.NON_DEFAULT_STANDARD) CYCLE
        instsi = instsi + 1

        ilcell = def_ilcell          
        IF (surface(i1)%ilswch.EQ.10000) ilcell = 1000  ! Special for neutral beams I believe

        n = LEN_TRIM(surface(i1)%surtxt)
        WRITE(fp06,90) surface(i1)%surtxt(1:n)  ! TRIM() doesn't work...
        WRITE(fp06,91) instsi,1,1
        WRITE(fp06,91) surface(i1)%iliin ,surface(i1)%ilside,
     .                 surface(i1)%ilswch,0             ,
     .                 surface(i1)%iltor ,surface(i1)%ilcol ,
     .                 0,ilcell,0,0

        IF (surface(i1)%reflect.EQ.LOCAL) THEN
          WRITE(fp06,91) 1,0
          WRITE(fp06,92) surface(i1)%material,surface(i1)%ewall 
          WRITE(fp06,92) surface(i1)%recycf,surface(i1)%recyct,
     .                   0.0,1.0,0.5,1.0
        ENDIF
      ENDDO

c      WRITE(eirfp,*) 'DONE'

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
c subroutine: WriteBlock03b
c
c
c
c
      SUBROUTINE WriteBlock03b_06
      USE mod_eirene06_parameters
      USE mod_eirene06
      IMPLICIT none


c      DATA material / 9642., 1206., 18474., 904./


      WRITE(fp06,80) '*** 3b. DATA FOR ADDITIONAL SURFACES (OSM)'
      IF (beam.EQ.1) THEN
        WRITE(fp06,81) 2

        WRITE(fp06,80) '* testing'
        WRITE(fp06,80) '  2.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     1     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) '  1.9400E 01  1.0000E 02 -1.0000E-05'//
     .                 '  1.9400E 01 -1.0000E 02  1.0000E-05'

        WRITE(fp06,80) '* universe, absorbing'
        WRITE(fp06,80) '  0.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     0     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) ' -5.2500E 04 -2.0000E 02  0.0000E 00'//
     .                 '  0.0000E 00  1.0000E 00  1.0000E 00'
        WRITE(fp06,80) '  0.0000E 00  0.0000E 00  0.0000E 00'//
     .                 '  0.0000E 00'

      ELSE
        WRITE(fp06,81) 1
        WRITE(fp06,80) '* testing'
        WRITE(fp06,80) '  2.0000E 00  1.0000E 00'
        WRITE(fp06,80) '     2     2     0     0     1     2     0'//
     .                 '     0     0     0'
        WRITE(fp06,80) '  1.9400E 01  1.0000E 02 -1.0000E-05'//
     .                 '  1.9400E 01 -1.0000E 02  1.0000E-05'
c        WRITE(fp06,81) 0
      ENDIF


 80   FORMAT(A)
 81   FORMAT(20(I6))
 82   FORMAT(1P,20(E12.4))
 83   FORMAT(1P,20(E12.5))       


      RETURN
      STOP
      END
c
c ======================================================================
c
c subroutine: WriteBlock05
c
c
c
c
      SUBROUTINE WriteBlock05_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER   fp1,fp2,nmass,i1,ibgk,nplsi,iscd2
      CHARACTER buffer*200
      CHARACTER*8 psym(2,5)

c     nmass = NINT(crmb)
      nmass = NINT(2.0) ! Replace with assignment of mod_eirene06 variable with CRMB...


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

      WRITE(fp06,90) '*** 5. DATA FOR PLASMA-BACKGROUND (OSM) 2006'
      WRITE(fp06,90) '*BULK ION SPECIES CARDS:  NPLSI ION SPECIES '//
     .               'ARE CONSIDERED, NPLSI='
      WRITE(fp06,91) nplsi
      WRITE(fp06,94) 1,psym(nmass,1),nmass,1,1,1,1,-4,0,1      ! D+
      WRITE(fp06,91) 8,115,111,0,30000
!      WRITE(fp06,92) 0.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
      WRITE(fp06,92) 35.0,0.0,0.0,0.0,1.0                      ! eirsrcmul*eirscale(11)
c      WRITE(fp06,92) 16.0,0.0,0.0,0.0,1.0E-15                      ! No volume recombination

      ibgk = 0

      IF (bgk.EQ.3) THEN
c...    BGK:
        ibgk = 4
        iscd2 = 614
        WRITE(fp06,94) 2,psym(nmass,2),  nmass,1,1,0,1,-1,0,0
        WRITE(fp06,94) 3,psym(nmass,3),2*nmass,2,2,0,1,-1,0,0
        WRITE(fp06,94) 4,psym(nmass,4),  nmass,1,1,0,1,-1,0,0
        WRITE(fp06,94) 5,psym(nmass,5),2*nmass,2,2,0,1,-1,0,0
      ENDIF

      IF (photons.EQ.1.OR.photons.EQ.2) THEN
c...    Photons:

        IF (photons.EQ.1) THEN
          WRITE(fp06,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0
          WRITE(fp06,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1
          WRITE(fp06,91) 19,0,210,iscd2,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        ELSE
          WRITE(fp06,94) 2+ibgk,'D_1     ',nmass,1,1,0,1,-4,0,0,0,0,0,0,
     .                   'FORT.13   ',1
          WRITE(fp06,91) 2

          WRITE(fp06,94) 3+ibgk,'D_2g    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                   'COLRAD    ',1
          WRITE(fp06,91) 19,0,210,iscd2,0
          WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
          WRITE(fp06,95) 2,4,0,'AMJUEL H.122.1.5b   OT'
        ENDIF

        WRITE(fp06,94) 4+ibgk,'D_2c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 19,0,210,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8b   OT'

        WRITE(fp06,94) 5+ibgk,'D_3g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 20,0,310,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 6+ibgk,'D_3c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 20,0,310,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8a   OT'

        WRITE(fp06,94) 7+ibgk,'D_4g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 25,0,410,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 8+ibgk,'D_4c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 25,0,410,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8c   OT'

        WRITE(fp06,94) 9+ibgk,'D_5g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 27,0,510,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 10+ibgk,'D_5c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 27,0,510,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8d   OT'

        WRITE(fp06,94) 11+ibgk,'D_6g    ',nmass,1,1,0,1,-4,0,1
        WRITE(fp06,91) 29,0,610,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  

        WRITE(fp06,94) 12+ibgk,'D_6c    ',nmass,1,1,0,1,-4,0,1,0,0,0,0,
     .                 'COLRAD    ',1
        WRITE(fp06,91) 29,0,610,iscd2,0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0                  
        WRITE(fp06,95) 1,4,0,'AMJUEL H.122.1.8e   OT'

        WRITE(fp06,94) 13+ibgk,'Lyman_a ',nmass,1,1,0,1,-4,0,0
        WRITE(fp06,94) 14+ibgk,'Lyman_b ',nmass,1,1,0,1,-4,0,0
      ENDIF

c...
      WRITE(fp06,91) 5,-5,5,5,5
      WRITE(fp06,92) 1.0, 90.0,1.5,2.5,4.3,72.0             ! Te (bogus, not actually used in eirene)
      DO i1 = 1, nplsi
        WRITE(fp06,92) 1.0,200.0,3.0,1.0,4.3,72.0,REAL(i1)  ! Ti
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp06,92) 0.0,0.0,3.0,1.0,4.6,72.0,REAL(i1)    ! ni
      ENDDO
      DO i1 = 1, nplsi
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0,REAL(i1)    ! vi
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0
        WRITE(fp06,92) 0.0,0.0,0.0,0.0,0.0, 0.0
      ENDDO
      WRITE(fp06,92) 0.1,  0.1,0.0,0.0,0.0,72.0             ! B-field


c      IF (.FALSE..AND.bgk.EQ.0.AND.photons.EQ.0) THEN
c...    Standard (no BGK or photons):
c        WRITE(fp06,91) 1
c        WRITE(fp06,91) 5,-5,5,5,5
c        WRITE(fp06,92) 1.0, 90.0,1.5,2.5,4.3,72.0
c        WRITE(fp06,92) 1.0,200.0,3.0,1.0,4.3,72.0
c        WRITE(fp06,92) 0.0,  0.0,3.0,1.0,4.6,72.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06,92) 0.0,  0.0,0.0,0.0,0.0, 0.0
c        WRITE(fp06.1,0.0,0.0,0.0,72.0
c      ELSE
c        CALL ER('WriteBlock05','Invalid EIR_07OPT',*99)
c      ENDIF

      RETURN
90    FORMAT(A)
91    FORMAT(20(I6))
92    FORMAT(1P,20(E12.4))
94    FORMAT(I2,1X,A8,12(I3),1X,A10,1X,I2)
95    FORMAT(3I6,1X,A)

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
      SUBROUTINE WriteBlock06_06(fp1,fp2)
      USE mod_eirene06
      IMPLICIT none

      INTEGER   fp1,fp2
      CHARACTER buffer*200

      WRITE(fp2,90) '*** 6. DATA FOR GENERAL REFLECTION MODEL (DIVIMP)'

      WRITE(fp2,90) 'TF'

      IF (trim_data.EQ.1) THEN
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
c subroutine: WriteBlock06
c
c
c
c
      SUBROUTINE WriteBlock07_06
      USE mod_eirene06
      IMPLICIT none

      INTEGER fp1

      INTEGER   i1
      CHARACTER buffer*200

      WRITE(fp06,90) '*** 7. DATA FOR PRIMARY SOURCES OF NEUTRALS (OSM)'
      WRITE(fp06,91) nstrata
      WRITE(fp06,91) (strata(i1)%indsrc,i1=1,nstrata)
      WRITE(fp06,92) alloc
      DO i1 = 1, nstrata
        WRITE(fp06,90) strata(i1)%txtsou(1:LEN_TRIM(strata(i1)%txtsou)) 
        WRITE(fp06,90) 'FFFFF'  
        WRITE(fp06,91) strata(i1)%npts,strata(i1)%ninitl,
     .                 strata(i1)%nemods,1
        WRITE(fp06,92) strata(i1)%flux
        WRITE(fp06,90) strata(i1)%species_tag
        WRITE(fp06,91) strata(i1)%nspez
        WRITE(fp06,90) strata(i1)%distrib
        WRITE(fp06,91) 1
        WRITE(fp06,91) strata(i1)%inum,strata(i1)%indim,strata(i1)%insor
        WRITE(fp06,92) strata(i1)%sorwgt,strata(i1)%sorlim,
     .                 strata(i1)%sorind,0.0,1000.0
        WRITE(fp06,91) strata(i1)%nrsor,0,0,0,strata(i1)%nasor
        WRITE(fp06,92) strata(i1)%sorad(1:6)
        WRITE(fp06,92) strata(i1)%sorene,strata(i1)%soreni,
     .                 0.0,0.0,0.0,0.0
        WRITE(fp06,92) strata(i1)%sorcos,strata(i1)%sormax
      ENDDO

c      WRITE(eirfp,*) 'STRATA:',strata(3)%indim,strata(3)%insor


      RETURN
90    FORMAT(A)
91    FORMAT(20I6)
92    FORMAT(1P,20E12.4)
97    CALL ER('WriteBlock07','Unexpected end of file'        ,*99)
98    CALL ER('WriteBlock07','Problems reading template file',*99)
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
