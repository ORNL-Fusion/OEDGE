      MODULE COMNNL

      USE PRECISION
      USE PARMMOD
!pb      USE COMPRT , ONLY : NPARTT, MPARTT

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_COMNNL, DEALLOC_COMNNL, INIT_COMNNL

      REAL(DP), PUBLIC, SAVE ::
     R  DTIMV,  DTIMVI, DTIMVN, TIME0

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  RPART(:,:), RPARTC(:,:), RPARTW(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I  IPART(:,:), IPARTC(:,:),
     I  NPRNLS(:)

      INTEGER, PUBLIC, SAVE ::
     I  NPRNLI, IPRNLI, IPRNLS, IPRNL,
     I  NPTST,  NTMSTP, ITMSTP

      INTEGER, PUBLIC, SAVE ::
     I  NCMNNL, MCMNNL


      CONTAINS


      SUBROUTINE ALLOC_COMNNL

      IF (ALLOCATED(RPART)) RETURN

      ALLOCATE (RPART(NPRNL,NPARTT))
      ALLOCATE (RPARTC(NPRNL,NPARTT))
      ALLOCATE (RPARTW(0:NPRNL))

      ALLOCATE (IPART(NPRNL,MPARTT))
      ALLOCATE (IPARTC(NPRNL,MPARTT))
      ALLOCATE (NPRNLS(NSTRA))

      WRITE (55,*) ' COMNNL ',(2*NPRNL*NPARTT+NPRNL+1)*8 +
     .                        (2*NPRNL*MPARTT+NSTRA)*4

      CALL INIT_COMNNL

      RETURN
      END SUBROUTINE ALLOC_COMNNL


      SUBROUTINE DEALLOC_COMNNL

      IF (.NOT.ALLOCATED(RPART)) RETURN

      DEALLOCATE (RPART)
      DEALLOCATE (RPARTC)
      DEALLOCATE (RPARTW)

      DEALLOCATE (IPART)
      DEALLOCATE (IPARTC)
      DEALLOCATE (NPRNLS)

      RETURN
      END SUBROUTINE DEALLOC_COMNNL


      SUBROUTINE INIT_COMNNL

      RPART  = 0._DP
      RPARTC = 0._DP
      RPARTW = 0._DP

      IPART  = 0
      IPARTC = 0
      NPRNLS = 0

      RETURN
      END SUBROUTINE INIT_COMNNL

      END MODULE COMNNL