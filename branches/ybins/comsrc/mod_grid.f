!     -*-Fortran-*-
      MODULE mod_grid_divimp
      IMPLICIT none
      PUBLIC

!...  Variables to store the knot indices in the sonnet grid file, for use when
!     loading B2 data from Rhozansky:
      INTEGER, SAVE :: divimp_maxnks,divimp_maxnrs
      INTEGER, ALLOCATABLE, SAVE :: divimp_ik(:,:),divimp_ir(:,:)

      REAL*8 , ALLOCATABLE, SAVE :: d_rvertp(:,:),d_zvertp(:,:)

      CONTAINS

      SUBROUTINE ALLOC_GRID(MAXNKS,MAXNRS)
      INTEGER, INTENT(IN) :: MAXNKS,MAXNRS
      IF (ALLOCATED(divimp_ik)) THEN
        WRITE(0,*) 'ERROR MOD_GRID: Grid arrays already allocated'
        WRITE(0,*) 'HALTING CODE'
        STOP
      ELSE
        ALLOCATE(divimp_ik(MAXNKS,MAXNRS))
        ALLOCATE(divimp_ir(MAXNKS,MAXNRS))      
        divimp_ik = 0
        divimp_ir = 0
      ENDIF
      divimp_maxnks = MAXNKS
      divimp_maxnrs = MAXNRS
      RETURN
      END SUBROUTINE ALLOC_GRID
      

      SUBROUTINE DEALLOC_GRID
      IF (ALLOCATED(divimp_ik)) DEALLOCATE(divimp_ik)       
      IF (ALLOCATED(divimp_ir)) DEALLOCATE(divimp_ir)       
      RETURN
      END SUBROUTINE DEALLOC_GRID

      END MODULE mod_grid_divimp
!
! ======================================================================
!
      MODULE mod_grid
      IMPLICIT none
      PUBLIC


!      PUBLIC :: ALLOC_GRID, DEALLOC_GRID


      TYPE type_grid_cell
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        REAL*8  :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE type_grid_cell

      TYPE type_grid_body
        INTEGER :: irsep
        INTEGER :: irsep2
        INTEGER :: irwall
        INTEGER :: irtrap
        INTEGER :: nrs
        INTEGER :: ikti
        INTEGER :: ikto
        INTEGER :: nks(1000)
        REAL    :: psin(1000)
      ENDTYPE type_grid_body

      TYPE type_grid_wall
        INTEGER :: ptt
        INTEGER :: ptc
        REAL*8  :: pt1(2)
        REAL*8  :: pt2(2)
      ENDTYPE type_grid_wall

      INTEGER, PARAMETER :: GRD_FORMAT_SONNET = 1,  
     .                      GRD_FORMAT_GRID   = 2, 
     .                      XPT_SEARCH = 1,         
     .                      R_INWARD   = 2,         
     .                      P_FORWARD  = 3,         
     .                      R_OUTWARD  = 4,         
     .                      P_BACKWARD = 5,
     .                      LIMITER_SEARCH = 6         

      INTEGER, ALLOCATABLE :: imap(:,:)

      INTEGER :: nknot
      TYPE(type_grid_cell),ALLOCATABLE :: knot(:)
      TYPE(type_grid_body) :: grid_load

      INTEGER :: n_grid_wall
      TYPE(type_grid_wall) :: grid_wall(1000) 

      INTEGER   :: grd_format
      CHARACTER :: grd_filename*1024

!...  Variables to store the knot indices in the sonnet grid file, for use when
!     loading B2 data from Rhozansky:
!      INTEGER, SAVE :: divimp_maxnks,divimp_maxnrs
!      INTEGER, ALLOCATABLE, SAVE :: divimp_ik(:,:),divimp_ir(:,:)


      END MODULE mod_grid
