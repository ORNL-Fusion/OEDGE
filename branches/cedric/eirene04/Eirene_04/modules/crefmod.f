      MODULE CREFMOD

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CREFMOD, DEALLOC_CREFMOD, INIT_CREFMOD

      REAL(DP), ALLOCATABLE, PUBLIC :: HFTR0(:,:,:),
     .                               HFTR1(:,:,:,:),
     .                               HFTR2(:,:,:,:,:),
     .                               HFTR3(:,:,:,:,:,:)


      CONTAINS


      SUBROUTINE ALLOC_CREFMOD

      IF (ALLOCATED(HFTR0)) RETURN

      ALLOCATE (HFTR0(NHD1,NHD2,NHD6))
      ALLOCATE (HFTR1(NHD1,NHD2,NHD3,NHD6))
      ALLOCATE (HFTR2(NHD1,NHD2,NHD3,NHD4,NHD6))
      ALLOCATE (HFTR3(NHD1,NHD2,NHD3,NHD4,NHD5,NHD6))

      WRITE (55,'(A,T25,I15)')
     .      ' CREFMOD ',(NHD1*NHD2*NHD6*(1+NHD3*(1+NHD4*
     .                   (1+NHD5))))*8

      CALL INIT_CREFMOD

      RETURN
      END SUBROUTINE ALLOC_CREFMOD


      SUBROUTINE DEALLOC_CREFMOD

      IF (.NOT.ALLOCATED(HFTR0)) RETURN

      DEALLOCATE (HFTR0)
      DEALLOCATE (HFTR1)
      DEALLOCATE (HFTR2)
      DEALLOCATE (HFTR3)

      RETURN
      END SUBROUTINE DEALLOC_CREFMOD


      SUBROUTINE INIT_CREFMOD

      HFTR0 = 0._DP
      HFTR1 = 0._DP
      HFTR2 = 0._DP
      HFTR3 = 0._DP

      RETURN
      END SUBROUTINE INIT_CREFMOD

      END MODULE CREFMOD

