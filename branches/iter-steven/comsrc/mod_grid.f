!     -*-Fortran-*-
      MODULE mod_grid
      IMPLICIT none
      PUBLIC


!      PUBLIC :: ALLOC_GRID, DEALLOC_GRID


      TYPE type_grid_cell
        INTEGER :: index,ik,ir,nv,rzone,zzone,xpt,map
        REAL*8  :: rcen,zcen,bratio,rv(4),zv(4)
      ENDTYPE type_grid_cell

      INTEGER, PARAMETER :: GRD_FORMAT_SONNET = 1


      INTEGER   :: grd_format
      CHARACTER :: grd_filename*1024

!...  Variables to store the knot indices in the sonnet grid file, for use when
!     loading B2 data from Rhozansky:

      INTEGER, SAVE :: divimp_maxnks,divimp_maxnrs
      INTEGER, ALLOCATABLE, SAVE :: divimp_ik(:,:),divimp_ir(:,:)

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


      END MODULE mod_grid
