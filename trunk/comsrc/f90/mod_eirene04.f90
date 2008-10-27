!     -*-Fortran-*-
  MODULE MOD_EIRENE04

    IMPLICIT NONE 

    PRIVATE

    PUBLIC :: ALLOC_SURFACE ,DEALLOC_SURFACE
    PUBLIC :: ALLOC_TRIANGLE,DEALLOC_TRIANGLE
    PUBLIC :: ALLOC_VERTEX  ,DEALLOC_VERTEX
    PUBLIC :: ALLOC_CELL    ,DEALLOC_CELL

    INTEGER    MAXSUR  ,MAXPTS  ,MAXLNK   ,MAXVER
    PARAMETER (MAXSUR=4,MAXPTS=4,MAXLNK=10,MAXVER=8)

    TYPE, PUBLIC :: type_surface  ! Change to surface 
      INTEGER :: type,subtype,num,index(10),orientation,zone
      INTEGER :: iliin,ilside,ilswch,ilcol,ilcell,iltor,reflect
      REAL    :: ewall,material,recyct,recycf
      CHARACTER*256 :: surtxt
!       Geometry:
      INTEGER :: nsur,nver
      INTEGER :: npts(MAXSUR),ipts(MAXPTS,MAXSUR)
      INTEGER :: nlnk(MAXSUR),ilnk(MAXLNK,MAXSUR)
      REAL*8  :: v(3,MAXVER)
    ENDTYPE type_surface

    TYPE, PUBLIC :: type_triangle
      INTEGER :: type,index(10),sideindex(10,3),zone
      INTEGER :: ver(3),map(3),sid(3),sur(3)
      REAL    :: bfield(3),efield(3),plasma(20)
    ENDTYPE type_triangle

    TYPE, PUBLIC :: type_cell
      INTEGER :: type,index(10),sideindex(10,4),zone
      INTEGER :: surface(4)
      REAL    :: bfield(3),efield(3),plasma(20)
      REAL*8  :: r(4),z(4)
    ENDTYPE type_cell

!...Code identifier:
    CHARACTER, PUBLIC :: fluid_code*256

!...Surface types:
    INTEGER, PUBLIC :: VESSEL_WALL    , NON_DEFAULT_STANDARD
    PARAMETER         (VESSEL_WALL = 1, NON_DEFAULT_STANDARD = 2)

!...Non-default standard surface sub-types:
    INTEGER, PUBLIC :: STRATUM    , MAGNETIC_GRID_BOUNDARY    , ADDITIONAL
    PARAMETER         (STRATUM = 1, MAGNETIC_GRID_BOUNDARY = 2, ADDITIONAL = 3)

!...Reflection model:
    INTEGER, PUBLIC :: GLOBAL    , LOCAL
    PARAMETER         (GLOBAL = 1, LOCAL = 2)

!...Surface types:
    INTEGER, PUBLIC :: MAGNETIC_GRID    , VACUUM_GRID
    PARAMETER         (MAGNETIC_GRID = 1, VACUUM_GRID = 2)

    TYPE(type_surface), PUBLIC, ALLOCATABLE, SAVE :: add(:)

!...Fluid code magnetic grid cells:
    INTEGER, PUBLIC, SAVE :: ncell
    TYPE(type_cell), PUBLIC, ALLOCATABLE, SAVE :: cell(:)
    INTEGER, PUBLIC, SAVE :: ntardat
    real, PUBLIC, ALLOCATABLE, SAVE :: tardat(:,:)

!...  Fluid code defined EIRENE geometry surfaces:
    INTEGER, PUBLIC, SAVE :: nadd 

!...  
    INTEGER, PUBLIC, SAVE :: fp04

!...Triangles:
    INTEGER, PUBLIC, SAVE :: ntri,nver
    REAL, PUBLIC, ALLOCATABLE, SAVE    :: ver(:,:)
    TYPE(type_triangle), PUBLIC, ALLOCATABLE, SAVE :: tri(:)

!...Block  1 variables:
    INTEGER, PUBLIC, SAVE :: time,niter,nfile

!...Block  3 variables:
    REAL, PUBLIC, SAVE :: wtemp,ttemp,wmater,tmater

!...Block  4 variables:
    INTEGER, PUBLIC, SAVE :: opacity, photons, bgk

!...Block  6 variables:
    INTEGER, PUBLIC, SAVE :: trim

!...Block 13 variables:
    INTEGER, PUBLIC, SAVE :: dtimv

    CONTAINS


    SUBROUTINE ALLOC_SURFACE(nadd)
    INTEGER, INTENT(IN) :: nadd
    IF (ALLOCATED(add)) THEN  
      DEALLOCATE (add)
    ENDIF
    ALLOCATE (add(nadd))
    RETURN
    END SUBROUTINE ALLOC_SURFACE

    SUBROUTINE DEALLOC_SURFACE
    IF (.NOT.ALLOCATED(add)) RETURN
    DEALLOCATE (add)
    RETURN
    END SUBROUTINE DEALLOC_SURFACE

    SUBROUTINE ALLOC_TRIANGLE(ntri)
    INTEGER, INTENT(IN) :: ntri
    IF (ALLOCATED(tri)) THEN  
      DEALLOCATE (tri)
    ENDIF
    ALLOCATE (tri(ntri))
    RETURN
    END SUBROUTINE ALLOC_TRIANGLE

    SUBROUTINE DEALLOC_TRIANGLE
    IF (.NOT.ALLOCATED(tri)) RETURN
    DEALLOCATE (tri)
    RETURN
    END SUBROUTINE DEALLOC_TRIANGLE

    SUBROUTINE ALLOC_VERTEX(nver)
    INTEGER, INTENT(IN) :: nver
    IF (ALLOCATED(ver)) THEN  
      DEALLOCATE (ver)
    ENDIF
    ALLOCATE (ver(nver,3))
    RETURN
    END SUBROUTINE ALLOC_VERTEX

    SUBROUTINE DEALLOC_VERTEX
    IF (.NOT.ALLOCATED(ver)) RETURN
    DEALLOCATE (ver)
    RETURN
    END SUBROUTINE DEALLOC_VERTEX

    SUBROUTINE ALLOC_CELL(ncell,ntar)
    INTEGER, INTENT(IN) :: ncell, ntar
    IF (ALLOCATED(cell))   DEALLOCATE (cell)
    IF (ALLOCATED(tardat)) DEALLOCATE (tardat)
    ALLOCATE (cell(ncell))
    ntardat = 0
    ALLOCATE (tardat(ntar,20))
    RETURN
    END SUBROUTINE ALLOC_CELL

    SUBROUTINE DEALLOC_CELL
    IF (.NOT.ALLOCATED(cell)) RETURN
    DEALLOCATE (cell)
    RETURN
    END SUBROUTINE DEALLOC_CELL


  END MODULE MOD_EIRENE04
