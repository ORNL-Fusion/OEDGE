!     -*-Fortran-*-
!
! ======================================================================
!
      MODULE MOD_EIRENE06_PARAMETERS

!...  Surface types:
      INTEGER, PUBLIC, PARAMETER :: VESSEL_WALL          = 1, 
     .                              NON_DEFAULT_STANDARD = 2,
     .                              HOLE_IN_GRID         = 3  ! OK, not really a "surface"...
    
!...  Non-default standard surface sub-types:
      INTEGER, PUBLIC, PARAMETER :: STRATUM                = 1, 
     .                              MAGNETIC_GRID_BOUNDARY = 2, 
     .                              ADDITIONAL             = 3
    
!...  Reflection model:
      INTEGER, PUBLIC, PARAMETER :: GLOBAL = 1, 
     .                              LOCAL  = 2
    
!...  Surface types:
      INTEGER, PUBLIC, PARAMETER :: MAGNETIC_GRID = 1, 
     .                              VACUUM_GRID   = 2

!...  Array sizes:
      INTEGER, PUBLIC, PARAMETER :: VOID_MAXNSEG  = 1000, 
     .                              VOID_MAXNHOLE = 100 ,
     .                              MAXNHISTORY   = 100 ,
     .                              MAXNGAUGE     = 20  ,
     .                              MAXNSTRATA    = 50

      END MODULE MOD_EIRENE06_PARAMETERS
!
! ======================================================================
!
      MODULE MOD_EIRENE06
      USE mod_eirene06_parameters
      IMPLICIT NONE 
      
      PRIVATE


      INTEGER    MAXSUR  ,MAXPTS  ,MAXLNK   ,MAXVER
      PARAMETER (MAXSUR=4,MAXPTS=4,MAXLNK=10,MAXVER=8)

      
      PUBLIC :: ALLOC_SURFACE ,DEALLOC_SURFACE
      PUBLIC :: ALLOC_TRIANGLE,DEALLOC_TRIANGLE
      PUBLIC :: ALLOC_VERTEX  ,DEALLOC_VERTEX
      PUBLIC :: ALLOC_CELL    ,DEALLOC_CELL
      
      
      TYPE, PUBLIC :: type_surface  ! Change to surface 
        INTEGER :: type,subtype,num,index(10),orientation,zone
        INTEGER :: iliin,ilside,ilswch,ilcol,ilcell,iltor,reflect
        REAL    :: ewall,material,recyct,recycf
        INTEGER :: ilspt,isrs
        REAL    :: recycs,recycc
        INTEGER :: hard           ! Force this surface to be separate from other surfaces (i.e. not automatically grouped with other surfaces with the same properties), for flux reporting from EIRENE
        CHARACTER*256 :: surtxt
        CHARACTER*256 :: sector   ! Toroidal sector list, don't change array size as initialisation in the code assumes *256
!       Geometry:  ! defunct?
        INTEGER :: nsur,nver
        INTEGER :: npts(MAXSUR),ipts(MAXPTS,MAXSUR)
        INTEGER :: nlnk(MAXSUR),ilnk(MAXLNK,MAXSUR)
        REAL*8  :: v(3,MAXVER)
      ENDTYPE type_surface
      
      TYPE, PUBLIC :: type_triangle
        INTEGER :: type
        INTEGER :: index(10)
        INTEGER :: sideindex(10,3)
        INTEGER :: zone
        INTEGER :: ver(3)
        INTEGER :: map(3)
        INTEGER :: sid(3)
        INTEGER :: sur(3)
        REAL    :: bfield(4)
        REAL    :: efield(3)
        REAL    :: plasma(20)  ! 20 = e_pot temp!
      ENDTYPE type_triangle

!     Void regions between the fluid grid and the first wall:
      TYPE, PUBLIC :: type_void
        INTEGER :: index
        INTEGER :: index_boundary
        REAL    :: resolution
        REAL*8  :: seg_resolution(VOID_MAXNSEG)
        INTEGER :: nseg_boundary
        REAL*8  :: seg_boundary1 (2,VOID_MAXNSEG)
        REAL*8  :: seg_boundary2 (2,VOID_MAXNSEG)
        INTEGER :: nseg_interior
        REAL*8  :: seg_interior1(2,VOID_MAXNSEG)
        REAL*8  :: seg_interior2(2,VOID_MAXNSEG)
        INTEGER :: nhole
        REAL*8  :: hole(VOID_MAXNHOLE)         
      ENDTYPE type_void

      TYPE, PUBLIC :: type_eirene_cell
        INTEGER :: type,index(10),sideindex(10,4),zone
        INTEGER :: surface(4)
        REAL    :: bfield(4),efield(3),plasma(20),e_pot
        REAL*8  :: r(4),z(4)
      ENDTYPE type_eirene_cell

      TYPE, PUBLIC :: type_strata
c...    Quantities set in OSM input file:
        REAL      :: type
        INTEGER   :: Z              ! atomic number
        INTEGER   :: A              ! atomic mass
        INTEGER   :: npts
        REAL      :: flux
        REAL      :: flux_fraction
        INTEGER   :: species
        INTEGER   :: species_index
        REAL      :: energy
        INTEGER   :: target
        INTEGER   :: range_cell(2)
        INTEGER   :: range_tube(2)
        REAL      :: cos
        REAL      :: cos_max
        CHARACTER :: note*512
c...    Quantities set in EIRENE interface routines:
        INTEGER   :: indsrc  ! ...
        CHARACTER :: txtsou*512
        INTEGER   :: ninitl
        INTEGER   :: nemods
        CHARACTER :: species_tag*4
        INTEGER   :: nspez
        CHARACTER :: distrib*5
        INTEGER   :: inum
        INTEGER   :: indim
        INTEGER   :: insor
        REAL      :: sorwgt
        REAL      :: sorlim
        REAL      :: sorind
        INTEGER   :: nrsor
        INTEGER   :: nasor
        REAL      :: sorad(6)
        REAL      :: sorene
        REAL      :: soreni
        REAL      :: sorcos
        REAL      :: sormax
c...    Quantities returned by EIRENE:
        INTEGER   :: ipanu
        REAL      :: fluxt
        REAL      :: ptrash
        REAL      :: etrash
      ENDTYPE type_strata

!...  Code identifier:
      CHARACTER, PUBLIC :: fluid_code*256
      
      
      TYPE(type_surface), PUBLIC, ALLOCATABLE, SAVE :: surface(:)
      
!...  Fluid code magnetic grid cells:
      INTEGER, PUBLIC, SAVE :: ncell
      TYPE(type_eirene_cell), PUBLIC, ALLOCATABLE, SAVE :: cell(:)

      INTEGER, PUBLIC, SAVE :: ntardat
      REAL, PUBLIC, ALLOCATABLE, SAVE :: tardat(:,:)

      INTEGER, PUBLIC, SAVE :: eir_pass

      INTEGER, PUBLIC, SAVE :: nvoid
      TYPE(type_void), PUBLIC, ALLOCATABLE, SAVE :: void(:)

      INTEGER, PUBLIC, SAVE :: nstrata ! ,osm_nstrata
      REAL   , PUBLIC, SAVE :: alloc
      TYPE(type_strata), PUBLIC, SAVE :: strata(100) ! ,osm_strata(100)
      
!...  Fluid code defined EIRENE geometry surfaces:
      INTEGER, PUBLIC, SAVE :: nsurface,default_surface,core_boundary
      
      INTEGER, PUBLIC, SAVE :: fp06
      
!...  Triangles:
      INTEGER, PUBLIC, SAVE :: ntri,nver
      REAL*8, PUBLIC, ALLOCATABLE, SAVE :: ver(:,:)  ! FIX ...need more precision...
c      REAL, PUBLIC, ALLOCATABLE, SAVE :: ver(:,:)
      TYPE(type_triangle), PUBLIC, ALLOCATABLE, SAVE :: tri(:)

!...  Block  1 variables:
      INTEGER, PUBLIC, SAVE :: time,niter,nfile,ncall,time_iteration
    
!...  Block  3 variables:
      REAL, PUBLIC, SAVE :: wtemp,ttemp,wmater,tmater,torus1,torus2
    
!...  Block  4 variables:
      INTEGER, PUBLIC, SAVE :: opacity, photons, bgk, ntorseg, beam  ! Do I need these "SAVE's"?
      REAL   , PUBLIC, SAVE :: torfrac
    
!...  Block  6 variables:
      INTEGER, PUBLIC, SAVE :: trim_data
    
!...  Block 13 variables:
      REAL, PUBLIC, SAVE :: dtimv,time0
    
      LOGICAL, PUBLIC, SAVE :: tetrahedrons, helium, time_dependent

!...  i/o:
      INTEGER, PUBLIC, SAVE :: eirfp 


      CONTAINS
c
c ----------------------------------------------------------------------    
c    
      SUBROUTINE ALLOC_SURFACE(nsurface)
      INTEGER, INTENT(IN) :: nsurface
      IF (ALLOCATED(surface)) THEN  
        DEALLOCATE(surface)
      ENDIF
      ALLOCATE (surface(nsurface))
      RETURN
      END SUBROUTINE ALLOC_SURFACE
c
c ----------------------------------------------------------------------    
c        
      SUBROUTINE DEALLOC_SURFACE
      IF (.NOT.ALLOCATED(surface)) RETURN
      DEALLOCATE(surface)
      RETURN
      END SUBROUTINE DEALLOC_SURFACE
c
c ----------------------------------------------------------------------    
c        
      SUBROUTINE ALLOC_TRIANGLE(ntri)
      INTEGER, INTENT(IN) :: ntri
      IF (ALLOCATED(tri)) THEN  
        DEALLOCATE (tri)
      ENDIF
      ALLOCATE (tri(ntri))
      RETURN
      END SUBROUTINE ALLOC_TRIANGLE
c
c ----------------------------------------------------------------------    
c      
      SUBROUTINE DEALLOC_TRIANGLE
      IF (.NOT.ALLOCATED(tri)) RETURN
      DEALLOCATE (tri)
      RETURN
      END SUBROUTINE DEALLOC_TRIANGLE
c
c ----------------------------------------------------------------------    
c      
      SUBROUTINE ALLOC_VERTEX(nver)
      INTEGER, INTENT(IN) :: nver
      IF (ALLOCATED(ver)) THEN  
        DEALLOCATE (ver)
      ENDIF
      ALLOCATE (ver(nver,3))
      RETURN
      END SUBROUTINE ALLOC_VERTEX
c
c ----------------------------------------------------------------------    
c      
      SUBROUTINE DEALLOC_VERTEX
      IF (.NOT.ALLOCATED(ver)) RETURN
      DEALLOCATE (ver)
      RETURN
      END SUBROUTINE DEALLOC_VERTEX
c
c ----------------------------------------------------------------------    
c      
      SUBROUTINE ALLOC_CELL(ncell,ntar)
      INTEGER, INTENT(IN) :: ncell, ntar
      IF (ALLOCATED(cell))   DEALLOCATE (cell)
      IF (ALLOCATED(tardat)) DEALLOCATE (tardat)
      ALLOCATE (cell(ncell))
      ntardat = 0
c      nstrata = 0
      ALLOCATE (tardat(ntar,20))
      RETURN
      END SUBROUTINE ALLOC_CELL
c
c ----------------------------------------------------------------------    
c      
      SUBROUTINE DEALLOC_CELL
      IF (ALLOCATED(cell)) DEALLOCATE (cell)
      IF (ALLOCATED(tardat)) DEALLOCATE (tardat)
      RETURN
      END SUBROUTINE DEALLOC_CELL

      END MODULE MOD_EIRENE06
!
! ========================================================================
!
      MODULE mod_eirene06_locals
      USE mod_geometry
      IMPLICIT none
       
      PUBLIC

      INTEGER, PARAMETER :: IND_IK     = 1,  ! Object index parameters
     .                      IND_IR     = 2,  
     .                      IND_IS     = 3,  
     .                      IND_ZONE   = 4,  
     .                      IND_PLASMA = 5,  
     .                      IND_BFIELD = 6   ! Vaccum zone outside standard grid, from external call to TRIANGLE

      INTEGER, PARAMETER :: IND_STDGRD  = 1,  ! Magnetic fluid grid side index, i.e. 12, 23, 34, 41
     .                      IND_TARGET  = 2,  ! Target (block 7 stratum in Eirene input file)
     .                      IND_SURFACE = 3   ! Surface (block 2A non-default surface in Eirene input file)

!     Triangle objects:
      INTEGER, SAVE :: ntry,ntrysrf,ntryvtx
      TYPE(type_object) , ALLOCATABLE, SAVE :: try(:)
      TYPE(type_srf)    , ALLOCATABLE, SAVE :: trysrf(:)
      REAL*8            , ALLOCATABLE, SAVE :: tryvtx(:,:)


!     Surface quantities:
      REAL, ALLOCATABLE, SAVE :: flux(:,:)


!     Volume quantities:      
      REAL, ALLOCATABLE, SAVE :: bfield(:,:)
      REAL, ALLOCATABLE, SAVE :: efield(:,:)
      REAL, ALLOCATABLE, SAVE :: plasma(:,:)


      END MODULE mod_eirene06_locals
!
! ======================================================================
!
      MODULE mod_eirene_history
      USE mod_eirene06_parameters
      IMPLICIT none
       
      PUBLIC

      TYPE, PUBLIC :: type_eirene_history
        INTEGER :: iiter
        INTEGER :: ngauge
        INTEGER :: gauge_nstrata    
        REAL    :: gauge_p_atm      (MAXNSTRATA,MAXNGAUGE)  ! [mTorr]
        REAL    :: gauge_p_mol      (MAXNSTRATA,MAXNGAUGE)
        REAL    :: gauge_parden_atm (MAXNSTRATA,MAXNGAUGE)
        REAL    :: gauge_parden_mol (MAXNSTRATA,MAXNGAUGE)
        REAL    :: gauge_egyden_atm (MAXNSTRATA,MAXNGAUGE)
        REAL    :: gauge_egyden_mol (MAXNSTRATA,MAXNGAUGE)
      ENDTYPE type_eirene_history


      INTEGER :: nhistory
      TYPE(type_eirene_history), SAVE :: history(MAXNHISTORY)

      END MODULE mod_eirene_history
!
! ======================================================================
!
