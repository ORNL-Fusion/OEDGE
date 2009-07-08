!     -*-F90-*-
MODULE mod_grid
  IMPLICIT none
  PRIVATE



  PUBLIC :: ALLOC_SONNET_INDEX_GRID, DEALLOC_SONNET_INDEX_GRID



  !...  Variables to store the knot indices in the sonnet grid file, for use when
  !     loading B2 data from Rhozansky:
  INTEGER, PUBLIC, ALLOCATABLE, SAVE :: sonnetik(:,:),sonnetir(:,:)



CONTAINS



  SUBROUTINE ALLOC_SONNET_INDEX_GRID(MAXNKS,MAXNRS)
    INTEGER, INTENT(IN) :: MAXNKS,MAXNRS
    !      INTEGER, INTENT(IN) :: nrs,nks(nrs)
    !      INTEGER :: maxnks,ir
    IF (ALLOCATED(sonnetik)) THEN
       WRITE(0,*) 'ERROR MOD_GRID: Grid arrays already allocated'
       WRITE(0,*) 'HALTING CODE'
       STOP
    ELSE
       !        maxnks = 0                   ! Find shorthand way of doing this... 
       !        DO ir = 1, nrs
       !          maxnks = MAX(maxnks,nks(ir))
       !        ENDDO
       !        ALLOCATE(sonnetik(maxnks,nrs+3))       
       !        ALLOCATE(sonnetir(maxnks,nrs+3))      
       ALLOCATE(sonnetik(MAXNKS,MAXNRS))
       ALLOCATE(sonnetir(MAXNKS,MAXNRS))      
    ENDIF
    RETURN
  END SUBROUTINE ALLOC_SONNET_INDEX_GRID


  SUBROUTINE DEALLOC_SONNET_INDEX_GRID
    IF (ALLOCATED(sonnetik)) THEN
       DEALLOCATE(sonnetik)       
       DEALLOCATE(sonnetir)       
    ENDIF
    RETURN
  END SUBROUTINE DEALLOC_SONNET_INDEX_GRID


END MODULE mod_grid
