!     -*-Mode:f90-*-
      MODULE mod_geometry
      IMPLICIT none
      PRIVATE

      LOGICAL, PUBLIC :: obj_modified, srf_modified, vtx_modified  ! *** make PRIVATE when the index mapping routines are moved to this module

      REAL*8, PUBLIC, PARAMETER :: D_DEGRAD=1.74532925199D-02
   
      PUBLIC :: Alloc_Vtx ,  
     .          Alloc_Srf ,  
     .          Alloc_Obj , DeleteVertex
   
      PUBLIC :: AddVertex ,   
     .          AddSurface,   
     .          AddObject,
     .          MatchSurface,
     .          LoadObjects,
     .          SaveGeometryData,
     .          SetupCell2Obj,
     .          geoClean

   
      INTEGER, PUBLIC, PARAMETER :: MAX3DVER=4 ,     ! ...
     .                              MAX3DSIDE=6      ! ...this must be an even number, see type_object definition...
   
      INTEGER, PUBLIC, PARAMETER :: MP_INITIALIZE = 1, 
     .                              MP_INCREASE_SIZE = 2
   
      INTEGER, PUBLIC, PARAMETER :: SPR_PLANAR_POLYGON = 4,  ! *** CHANCE TO SRF_ from SPR_ ...
     .                              SPR_LINE_SEGMENT   = 5
   
      INTEGER, PUBLIC, PARAMETER :: GRP_TRIANGLE      = 1, 
     .                              GRP_TETRAHEDRON   = 2, 
     .                              GRP_QUADRANGLE    = 3, 
     .                              GRP_MAGNETIC_GRID = 20,
     .                              GRP_VACUUM_GRID   = 21,
     .                              OBJ_MAXNINDEX     = 10

!..   Vertex data:
      INTEGER, PARAMETER :: VTX_SIZE = 10000
      REAL   , PARAMETER :: VTX_STEP = 1.5
      INTEGER,PUBLIC :: nvtx
      INTEGER :: maxvtx
      REAL*8, PUBLIC, ALLOCATABLE :: vtx(:,:)

!..   Surface data:
      TYPE, PUBLIC :: type_srf  ! 44 bytes
        INTEGER :: index(5) ! ???       ! Generaly index array
        INTEGER :: type  ! INT*1        ! Category: SP_PLANAR_POLYGON=4, SP_LINE_SEGMENT=5
        INTEGER :: mode  ! DROP?        ! vertex indexing scheme (sequential, or not, or use IVTX(1)=-1 instead?)
        INTEGER :: link  ! DROP?        ! Connection to another surface
        INTEGER :: obj   ! FLEX         ! Object association (for connection map building, waste of space, replace...)
        INTEGER :: side  ! FLEX         !   ...same...
        INTEGER :: nvtx                 ! Number of vertices
        INTEGER :: ivtx(MAX3DVER)       ! Index of vertices in VTX array
        INTEGER :: svtx  ! FLEX         ! Sum of vertex indeces
      ENDTYPE type_srf

      REAL    , PARAMETER :: SRF_VERSION = 1.0
      INTEGER , PARAMETER :: SRF_SIZE = 10000
      REAL    , PARAMETER :: SRF_STEP = 1.5
      INTEGER             :: maxsrf
      INTEGER , PUBLIC    :: nsrf    ! make private ???
      TYPE(type_srf), PUBLIC, ALLOCATABLE :: srf(:)  ! make private ???


      TYPE, PUBLIC :: type_group
!       Categorization:
        INTEGER :: origin   ! Magnetic/standard grid or vacuum grid
        INTEGER :: type     ! Triangle or tetrahedron (thus far)
        INTEGER :: colour   ! Not in use, yet 
      ENDTYPE type_group

      INTEGER         , PUBLIC :: ngrp
      TYPE(type_group), PUBLIC :: grp(100)


! Index parameter:
!
!    1 IND_IK
!    2 IND_IR
!    3 IND_IS
!    4 IND_PLASMA
!    4 IND_IMPURITY
!    4 IND_FIELD
!    5 IND_GROUP
!    6 IND_TRIANGLE
!    7 IND_INTEGRAL
!    7 IND_INTEGRAND ??? 
!    7 IND_ZONE

! 1 IND_CELL 
! 2 IND_FLUID
! 3 IND_KINETIC
! 4 IND_FIELD
! 
  
      TYPE, PUBLIC :: type_object
!       Association:
        INTEGER   :: group             ! Group association
        INTEGER   :: index(OBJ_MAXNINDEX)         ! General index array
!        INTEGER   :: object           ! Absolute index in full geometry model
!        INTEGER*2 :: group            ! Group association
!        INTEGER*2 :: zone             ! Zone (for domain decomposition)
!        INTEGER*2 :: ring             ! Ring on poloidal magnetic 
!        INTEGER*2 :: cell             ! Cell on ring (just for debugging?)
!        INTEGER   :: parent           ! Primal seed...
!        INTEGER*2 :: segment          ! Toroidal segment 
!        INTEGER   :: index(5)         ! General index array -- 5 enough?
!       Location:
        REAL      :: x                 ! Object geometry centre x-axis coordinate (=R)  ! KEEP?
        REAL      :: y                 !                        y-axis            (=Z)  ! KEEP?
        REAL      :: z                 !                        z-axis                  ! KEEP?
        REAL      :: phi               ! Toroidal angle                                 ! KEEP?
        INTEGER*2 :: segment(2)        ! Toroidal segment span of object  ! KEEP?
!       Connection map:
        INTEGER   :: omap(MAX3DSIDE)   ! Object map for side              ! BETTER WAY, VIA SRF?
        INTEGER*2 :: smap(MAX3DSIDE)   ! Side map for side                                
!       Geometry:
        INTEGER   :: nside             ! Number of sides 
        INTEGER   :: iside(MAX3DSIDE)  ! Index of (primary) surface(s) in SRF array
        REAL      :: volume            ! Object volume   ! KEEP?
      ENDTYPE type_object

      REAL   , PARAMETER :: OBJ_VERSION  = 1.0
      INTEGER, PARAMETER :: OBJ_SIZE = 10000
      REAL   , PARAMETER :: OBJ_STEP = 1.5
      INTEGER            :: maxobj
      INTEGER, PUBLIC    :: nobj    ! make private ???
      TYPE(type_object), PUBLIC, ALLOCATABLE :: obj(:)  ! make private ???

      INTEGER, PUBLIC :: nvtxmap
      INTEGER, PUBLIC, ALLOCATABLE :: vtxmap(:)    

      INTEGER, PUBLIC :: geofp

      INTEGER, PUBLIC, ALLOCATABLE :: 
     .  obj_tube    (:)

      REAL*8, PUBLIC, ALLOCATABLE :: 
     .  obj_centroid(:,:),
     .  obj_phi     (:),
     .  obj_volume  (:),
     .  obj_distance(:,:)

      INTEGER, PUBLIC, PARAMETER :: 
     .  MODE_OBJ_CENTRE   = 1,  ! *** SHOULD CHANGE TO _CENTROID? ***
     .  MODE_OBJ_PHI      = 2, 
     .  MODE_OBJ_VOLUME   = 3,
     .  MODE_OBJ_TUBE     = 4,
     .  MODE_OBJ_DISTANCE = 5

      INTEGER, PUBLIC, ALLOCATABLE :: 
     .  srf_obj     (:), 
     .  srf_side    (:),
     .  srf_centroid(:),
     .  srf_area    (:),
     .  srf_sum     (:)  ! Wait, actually need this all the time..? 

      INTEGER, PUBLIC, PARAMETER :: 
     .  MODE_SRF_OBJ    = 101,
     .  MODE_SRF_SIDE   = 102, 
     .  MODE_SRF_AREA   = 103,
     .  MODE_SRF_SUM    = 104

      INTEGER, PUBLIC, ALLOCATABLE :: 
     .  cell2obj(:)

      INTEGER, PUBLIC :: geo_output
     
      CONTAINS
!
! ======================================================================
!
      SUBROUTINE GEO_ER(message)
      IMPLICIT none
      CHARACTER message*(*)
      WRITE(0,*) 'GEO_ER: '//message
      WRITE(0,*) 'HALTING CODE'
      STOP
      END SUBROUTINE GEO_ER
! ======================================================================
!
      SUBROUTINE Alloc_Vtx(memsize,mode)
      INTEGER :: memsize,mode
      REAL*8, ALLOCATABLE :: tmpvtx(:,:)
      INTEGER :: istat,i1

      vtx_modified = .TRUE.

      SELECTCASE(mode)

        CASE (MP_INITIALIZE) 
          IF (ALLOCATED(vtx)) THEN
            WRITE(0,*) 'ERROR ALLOC_VTX: VTX ARRAY ALREADY ALLOCATED'
            STOP
          ELSE
            nvtx = 0
            maxvtx = memsize       
            IF (memsize.EQ.-1) maxvtx = VTX_SIZE
            ALLOCATE(vtx(3,maxvtx),STAT=istat)
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_VTX: PROBLEM ALLOCATING VTX '//
     .                   'ARRAY'
              STOP      
            ENDIF
          ENDIF

        CASE (MP_INCREASE_SIZE)
          IF (.NOT.ALLOCATED(vtx)) THEN
            WRITE(0,*) 'ERROR ALLOC_VTX: VTX ARRAY NOT ALLOCATED'
            STOP
          ELSE
            IF (geo_output.GE.2) 
     .        WRITE(0,*) 'ALLOC_VTX: INCREASING SIZE',nvtx
            ALLOCATE(tmpvtx(3,nvtx),STAT=istat)     
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_VTX: BAD TMPVTX'
              STOP      
            ENDIF
            DO i1 = 1, 3
              tmpvtx(i1,1:nvtx) = vtx(i1,1:nvtx)
            ENDDO
            DEALLOCATE(vtx)
            IF (memsize.EQ.-1) THEN
              maxvtx = INT(maxvtx * VTX_STEP)
            ELSE
              maxvtx = maxvtx + memsize
            ENDIF
            ALLOCATE(vtx(3,maxvtx))
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_VTX: CANNOT REALLOCATE VTX'
              STOP      
            ENDIF
            DO i1 = 1, 3
              vtx(i1,1:nvtx) = tmpvtx(i1,1:nvtx)
            ENDDO
            DEALLOCATE(tmpvtx)
          ENDIF

        CASE DEFAULT
          WRITE(0,*) 'ERROR ALLOC_VTX: UNRECOGNIZED MODE VALUE'
          STOP      

      ENDSELECT
      RETURN
      END SUBROUTINE ALLOC_VTX
!
! ======================================================================
!
      SUBROUTINE Alloc_srf(memsize,mode)
      INTEGER :: memsize,mode
      TYPE(type_srf), ALLOCATABLE :: tmpsrf(:)
      INTEGER :: istat

      srf_modified = .TRUE.

      SELECTCASE (mode)

        CASE (MP_INITIALIZE)
          IF (ALLOCATED(srf)) THEN
            WRITE(0,*) 'ERROR ALLOC_SRF: SRF ARRAY ALREADY '//
     .                 'ALLOCATED'
            STOP
          ELSE
            nsrf = 0
            maxsrf = memsize       
            IF (memsize.EQ.-1) maxsrf = SRF_SIZE
            ALLOCATE(    srf(maxsrf),STAT=istat)
            IF (istat.NE.0) CALL GEO_ER('ERROR ALLOC_SRF: '//
     .                                  'PROBLEM ALLOCATING SRF ARRAY')
          ENDIF

        CASE (MP_INCREASE_SIZE)
          IF (.NOT.ALLOCATED(srf)) THEN
            WRITE(0,*) 'ERROR ALLOC_SRF: SRF ARRAY NOT ALLOCATED'
            STOP
          ELSE
            IF (geo_output.GE.2)
     .        WRITE(0,*) 'ALLOC_SRF: INCREASING SIZE',nsrf
            ALLOCATE(tmpsrf(nsrf),STAT=istat)     
            IF (istat.NE.0) CALL GEO_ER('ERROR ALLOC_SRF: BAD '//
     .                                  'TMPSRF')
            tmpsrf    (1:nsrf) = srf    (1:nsrf)
            DEALLOCATE(srf)
            IF (memsize.EQ.-1) THEN
              maxsrf = INT(maxsrf * SRF_STEP)
            ELSE
              maxsrf = maxsrf + memsize
            ENDIF
            ALLOCATE(srf(maxsrf))
            IF (istat.NE.0) CALL GEO_ER('ERROR ALLOC_SRF: CANNOT '//
     .                                  'REALLOCATE SRF')
            srf(1:nsrf) = tmpsrf(1:nsrf)
            DEALLOCATE(tmpsrf)
          ENDIF

        CASE DEFAULT
          WRITE(0,*) 'ERROR ALLOC_SRF: UNRECOGNIZED MODE VALUE'
          STOP      
      ENDSELECT
      RETURN
99    STOP
      END SUBROUTINE ALLOC_SRF
!
! ======================================================================
!
      SUBROUTINE Alloc_obj(memsize,mode)
      INTEGER :: memsize,mode
      TYPE(type_object), ALLOCATABLE :: tmpobj(:)
      INTEGER :: istat,iobj

      obj_modified = .TRUE.

      SELECTCASE (mode)

        CASE (MP_INITIALIZE)
          IF (ALLOCATED(obj)) THEN
            WRITE(0,*) 'ERROR ALLOC_OBJ: OBJ ARRAY ALREADY ALLOCATED'
            STOP
          ELSE
            nobj = 0
            maxobj = memsize       
            IF (memsize.EQ.-1) maxobj = OBJ_SIZE
            ALLOCATE(obj(maxobj),STAT=istat)
            DO iobj = 1, maxobj
              obj(iobj)%index = 0
              obj(iobj)%omap  = 0
              obj(iobj)%smap  = 0
              obj(iobj)%nside = 0
              obj(iobj)%iside = 0
            ENDDO  
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_OBJ: PROBLEM ALLOCATING OBJ '//
     .                   'ARRAY'
              STOP      
            ENDIF
          ENDIF

        CASE (MP_INCREASE_SIZE)
          IF (.NOT.ALLOCATED(obj)) THEN
            WRITE(0,*) 'ERROR ALLOC_OBJ: OBJ ARRAY NOT ALLOCATED'
            STOP
          ELSE
            IF (geo_output.GE.2) 
     .        WRITE(0,*) 'ALLOC_OBJ: INCREASING SIZE',nobj
            ALLOCATE(tmpobj(nobj),STAT=istat)     
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_OBJ: BAD TMPOBJ'
              STOP      
            ENDIF
            tmpobj(1:nobj) = obj(1:nobj)
            DEALLOCATE(obj)
            IF (memsize.EQ.-1) THEN
              maxobj = INT(maxobj * OBJ_STEP)
            ELSE
              maxobj = maxobj + memsize
            ENDIF
            ALLOCATE(obj(maxobj))
            DO iobj = 1, maxobj
              obj(iobj)%index = 0
              obj(iobj)%omap  = 0
              obj(iobj)%smap  = 0
              obj(iobj)%nside = 0
              obj(iobj)%iside = 0
            ENDDO  
            IF (istat.NE.0) THEN
              WRITE(0,*) 'ERROR ALLOC_OBJ: CANNOT REALLOCATE OBJ'
              STOP      
            ENDIF
            obj(1:nobj) = tmpobj(1:nobj)
            DEALLOCATE(tmpobj)
          ENDIF

        CASE DEFAULT
          WRITE(0,*) 'ERROR ALLOC_OBJ: UNRECOGNIZED MODE VALUE'
          STOP      
      ENDSELECT
      RETURN
      END SUBROUTINE ALLOC_OBJ
!
! ======================================================================
!
      INTEGER FUNCTION AddVertex(newvtx)
      IMPLICIT none

      REAL*8 newvtx(3)

      REAL*8 , PARAMETER :: DTOL = 1.0D-07  ! Tolerance

      INTEGER ivtx

      IF (.NOT.ALLOCATED(vtx)) THEN
!...    Initial allocation of memory for VTX (vertex list) array:
        nvtx = 0
        CALL Alloc_vtx(-1,MP_INITIALIZE)
      ELSEIF (nvtx+1.GE.maxvtx) THEN
!...    The size of VTX is not sufficient, so make it bigger:
        CALL Alloc_vtx(-1,MP_INCREASE_SIZE)
      ENDIF

!...  Scan over existing vertices to see if it already exists:
      DO ivtx = MAX(1,nvtx-128), nvtx   
        IF (DABS(vtx(1,ivtx)-newvtx(1)).LT.DTOL.AND.   
     .      DABS(vtx(2,ivtx)-newvtx(2)).LT.DTOL.AND.   
     .      DABS(vtx(3,ivtx)-newvtx(3)).LT.DTOL) THEN
!...      Found one:
!         WRITE(0,*) 'ADDVERTEX: FOUND DUPLICATE!'
          AddVertex = ivtx
          RETURN
        ENDIF
      ENDDO

!...  Add the vertex to the list:
      nvtx = nvtx + 1
      vtx(1:3,nvtx) = newvtx(1:3)
      AddVertex = nvtx

      vtx_modified = .TRUE.

      RETURN
 99   STOP

      END FUNCTION AddVertex
!
! ======================================================================
!
      SUBROUTINE DeleteVertex(ivtx_delete)
      IMPLICIT none

      INTEGER, INTENT(IN) :: ivtx_delete

      INTEGER i1,ivtx,isrf

!      DO i1 = ivtx_delete, nvtx-1
!        vtx(1:3,i1) = vtx(1:3,i1+1)
!      ENDDO
      vtx(1,ivtx_delete:nvtx-1) = vtx(1,ivtx_delete+1:nvtx)
      vtx(2,ivtx_delete:nvtx-1) = vtx(2,ivtx_delete+1:nvtx)
      vtx(3,ivtx_delete:nvtx-1) = vtx(3,ivtx_delete+1:nvtx)
      nvtx = nvtx - 1

!     Reduce memory use if possible!

      DO isrf = 1, nsrf
        DO ivtx = 1, srf(isrf)%nvtx
          IF (srf(isrf)%ivtx(ivtx).GE.ivtx_delete)   
     .      srf(isrf)%ivtx(ivtx) = srf(isrf)%ivtx(ivtx) - 1      
        ENDDO
        srf(isrf)%svtx = SUM(srf(isrf)%ivtx(1:srf(isrf)%nvtx))
      ENDDO

      END SUBROUTINE DeleteVertex
!
! ======================================================================
!
! function: AddSurface
!
!
      INTEGER FUNCTION AddSurface(newsrf)
      IMPLICIT none

      TYPE(type_srf) newsrf
      INTEGER imatch,i1,i2,isrf,sum_newsrf
      LOGICAL match


      IF (.NOT.ALLOCATED(srf)) THEN
!...    Check that array allocated.  This must be allocated/initialized
!       before the call to this array to be sure that NSRF=0:
        nsrf = 0
        CALL ALLOC_SRF(-1,MP_INITIALIZE)
      ELSEIF (nsrf+1.GE.maxsrf) THEN
!...    The size of SRF is not sufficient, so make it bigger:
        CALL ALLOC_SRF(-1,MP_INCREASE_SIZE)
      ENDIF

!...  Scan over existing surfaces to see if it already exists:
      match = .FALSE.
      sum_newsrf = SUM(newsrf%ivtx(1:newsrf%nvtx))
!     IF (newsrf%nvtx.EQ.3) WRITE(0,*) 'NSRF, CHECK:',nsrf,sum_newsrf
!      DO isrf = 1, 0
      DO isrf = MAX(1,nsrf-128), nsrf
!      DO isrf = 1, nsrf

!        IF (srf(isrf)%nvtx.LT.3) CYCLE

!        WRITE(0,*) 'ISRF CHECK:',isrf,srf(isrf)%svtx

        IF (sum_newsrf.NE.srf(isrf)%svtx) CYCLE  ! Check more things?  Make region lists to speed up search?

!       WRITE(0,*) 'CHECKING!'
!       WRITE(0,*) newsrf%ivtx(1:3)
!       DO i1 = 1, 3
!         WRITE(0,*) '   ',vtx(1:3,newsrf%ivtx(i1))
!       ENDDO 
!       WRITE(0,*) srf(isrf)%ivtx(1:3)
!       DO i1 = 1, 3
!         WRITE(0,*) '   ',vtx(1:3,srf(isrf)%ivtx(i1))
!       ENDDO 

        imatch = 0
        DO i1 = 1, newsrf%nvtx
          DO i2 = 1, srf(isrf)%nvtx
            IF (newsrf%ivtx(i1).EQ.srf(isrf)%ivtx(i2)) imatch=imatch + 1
          ENDDO
        ENDDO

!        WRITE(0,*) 'IMATCH:',imatch

        IF (imatch.EQ.newsrf%nvtx) THEN
          match = .TRUE.
          EXIT
        ENDIF
      ENDDO

      IF (match) THEN
!...    Refer to existing surface with matching verticies:
        AddSurface = -isrf
!        STOP 'HOLY CRAP!'
      ELSE
!...    Add surface to the list:
        nsrf = nsrf + 1
        srf(nsrf) = newsrf
        srf(nsrf)%svtx = sum_newsrf  ! Do I still need this...?
        srf(nsrf)%link = 0
        AddSurface = nsrf
      ENDIF

      srf_modified = .TRUE.

      RETURN
 99   STOP
      END FUNCTION AddSurface

!
! ======================================================================
!
! function: MatchSurface
!
!
      LOGICAL FUNCTION MatchSurface(isrf1,isrf2)
      IMPLICIT none

      INTEGER, INTENT(IN) :: isrf1,isrf2
      INTEGER mode,c,i1,i2
  
      mode = 0

      MatchSurface = .FALSE.     

      SELECTCASE (mode)
!...    Strict:
        CASE (0)  

          IF     (isrf1.EQ.isrf2) THEN 
            MatchSurface = .TRUE.

!            WRITE(0,*) 'HERE! A'

          ELSEIF (srf(isrf1)%svtx.EQ.srf(isrf2)%svtx.AND.     ! More checks...?
     .            srf(isrf1)%nvtx.EQ.srf(isrf2)%nvtx) THEN

            c = 0
            DO i1 = 1, srf(isrf1)%nvtx
              DO i2 = 1, srf(isrf2)%nvtx
                IF (srf(isrf1)%ivtx(i1).EQ.srf(isrf2)%ivtx(i2)) c=c+1
              ENDDO
            ENDDO

            IF (c.EQ.srf(isrf1)%nvtx) MatchSurface = .TRUE.

!            IF (c.EQ.srf(isrf1)%nvtx) WRITE(0,*) 'HERE! B'

          ENDIF

        CASE DEFAULT
          STOP 'UNRECOGNIZED MODE'
      ENDSELECT

!      WRITE(0,*) '    MATCH:',isrf1,isrf2,MatchSurface
!      WRITE(0,*) '         :',srf(isrf1)%svtx,srf(isrf2)%svtx



      RETURN
99    STOP
      END FUNCTION MatchSurface
!
! ======================================================================
!
! function: AddObject
!
!
      INTEGER FUNCTION AddObject(newobj)
      IMPLICIT none
      
      TYPE(type_object) newobj

      AddObject = -1

      IF (.NOT.ALLOCATED(obj)) THEN
!...    Check that array allocated.  This must be allocated/initialized
!       before the call to this array to be sure that NSRF=0:
        nobj = 0
        CALL Alloc_obj(-1,MP_INITIALIZE)
      ELSEIF (nobj+1.GE.maxobj) THEN
!...    The size of SRF is not sufficient, so make it bigger:
        CALL Alloc_obj(-1,MP_INCREASE_SIZE)
      ENDIF

!...  Add surface to the list:
      nobj = nobj + 1
      obj(nobj) = newobj

      AddObject = nobj

      obj_modified = .TRUE.

      RETURN
 99   STOP
      END FUNCTION AddObject
!
! ----------------------------------------------------------------------
!  subroutine: SaveGeometryData
!
      SUBROUTINE SaveGeometryData(fname)
      IMPLICIT none

      CHARACTER*(*) fname
      INTEGER fp,i1

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='REPLACE',ERR=98)            
      WRITE(fp,ERR=98) 1.0,ngrp,nobj,nsrf,nvtx
      WRITE(fp,ERR=98) (obj(    i1),i1=1,nobj)
      WRITE(fp,ERR=98) (grp(    i1),i1=1,ngrp)
      WRITE(fp,ERR=98) (srf(    i1),i1=1,nsrf)
      WRITE(fp,ERR=98) (vtx(1:3,i1),i1=1,nvtx)
      CLOSE (fp)
      
      RETURN
 98   CALL GEO_ER('SaveGeometryData: Problems writing data file')
 99   STOP
      END SUBROUTINE SaveGeometryData
!
! ----------------------------------------------------------------------
!
      SUBROUTINE LoadObjects(fname,status)
      IMPLICIT none

      INTEGER   status
      CHARACTER fname*(*)

      INTEGER fp,i1,ngrp1,nobj1,nsrf1,nvtx1
      REAL    version

      status = 0

      fp = 99
      OPEN(UNIT=fp,FILE=fname(1:LEN_TRIM(fname)),ACCESS='SEQUENTIAL',
     .     FORM='UNFORMATTED',STATUS='OLD',ERR=98)            
      READ(fp,ERR=98) version,ngrp1,nobj1,nsrf1,nvtx1

      IF (version.NE.1.0)
     .  CALL GEO_ER('LoadObjects: Unsupported version')

      CALL Alloc_obj(nobj1,MP_INITIALIZE)
      CALL Alloc_srf(nsrf1,MP_INITIALIZE)
      CALL Alloc_vtx(nvtx1,MP_INITIALIZE)

      READ(fp,ERR=98) (obj(    i1),i1=1,nobj1)
      READ(fp,ERR=98) (grp(    i1),i1=1,ngrp1)
      READ(fp,ERR=98) (srf(    i1),i1=1,nsrf1)
      READ(fp,ERR=98) (vtx(1:3,i1),i1=1,nvtx1)

      ngrp = ngrp1
      nobj = nobj1
      nsrf = nsrf1
      nvtx = nvtx1

      CLOSE (fp)
      
      RETURN
 98   CALL GEO_ER('LoadObjects: Problems reading data file')
 99   status = -1
      WRITE(0,*) '    FILE NAME: >'//fname(1:LEN_TRIM(fname))//'<'
      RETURN
      END SUBROUTINE LoadObjects
!
! ----------------------------------------------------------------------
!
      SUBROUTINE SetupCell2Obj(ncell,index)   ! *** REPLACED *** by GetObject function, need to modify the code...
      IMPLICIT none

      INTEGER, INTENT(IN) :: ncell, index
      INTEGER iobj

      IF (ALLOCATED(cell2obj)) DEALLOCATE(cell2obj)
      ALLOCATE(cell2obj(ncell))
      DO iobj = 1, nobj
        IF (grp(obj(iobj)%group)%origin.EQ.GRP_MAGNETIC_GRID.AND.
     .      grp(obj(iobj)%group)%type  .EQ.GRP_QUADRANGLE) 
     .    cell2obj(obj(iobj)%index(index)) = iobj
      ENDDO

      RETURN
99    STOP
      END SUBROUTINE SetupCell2Obj
!
! ----------------------------------------------------------------------
!
      SUBROUTINE geoClean
      IMPLICIT none
      nobj = 0
      nsrf = 0
      nvtx = 0
      IF (ALLOCATED(obj)) DEALLOCATE(obj) 
      IF (ALLOCATED(srf)) DEALLOCATE(srf) 
      IF (ALLOCATED(vtx)) DEALLOCATE(vtx)
      RETURN
99    STOP
      END SUBROUTINE geoClean
!
! ======================================================================
!

      END MODULE MOD_GEOMETRY

