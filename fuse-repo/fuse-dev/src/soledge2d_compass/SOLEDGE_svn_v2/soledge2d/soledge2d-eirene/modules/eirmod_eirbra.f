      MODULE EIRMOD_EIRBRA
 
C  NEUTRAL SOURCE TERMS: SNI,SMO,SEE,SEI (EIRENE ---> BRAAMS)
 
      USE EIRMOD_PRECISION
!pb      USE PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_EIRBRA, EIRENE_DEALLOC_EIRBRA, 
     P          EIRENE_INIT_EIRBRA
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SNI(:,:,:,:), SMO(:,:,:,:),
     R SEE(:,:,:),   SEI(:,:,:)
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R VOLSUMN(:), VOLSUMM(:), VOLSUMEI(:), VOLSUMEE(:)
 
      INTEGER, SAVE :: NDXD, NDYD, NFLD, NSTRAD
 
C AK
      REAL(DP),PUBLIC,ALLOCATABLE,SAVE ::
     .                srcstrn(:),flxspci(:,:),
     .                srccrfc(:,:)
c*** srcstrn: intensities of different neutral sources (fluxt)
c*** flxspci: fluxes of plasma ions to the recycling surfaces
c*** srccrfc: source correction factors (from infcop)
 
C AK END
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_EIRBRA(NDXP, NDYP, NFL, NSTRA)
 
      INTEGER, INTENT(IN) :: NDXP, NDYP, NFL, NSTRA
 
      IF (ALLOCATED(SNI)) RETURN
 
      NDXD = NDXP+1
      NDYD = NDYP+1
      NFLD = NFL
      NSTRAD = NSTRA
 
      ALLOCATE (SNI(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SMO(0:NDXD,0:NDYD,NFLD,NSTRAD))
      ALLOCATE (SEE(0:NDXD,0:NDYD,NSTRAD))
      ALLOCATE (SEI(0:NDXD,0:NDYD,NSTRAD))
 
      ALLOCATE (VOLSUMN(NSTRA))
      ALLOCATE (VOLSUMM(NSTRA))
      ALLOCATE (VOLSUMEI(NSTRA))
      ALLOCATE (VOLSUMEE(NSTRA))
 
      WRITE (55,'(A,T25,I15)')
     .             ' EIRBRA ',((NDXD+1)*(NDYD+1)*NSTRAD*2*NFLD+
     .                        4*NSTRA)*8
 
      CALL EIRENE_INIT_EIRBRA
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_EIRBRA
 
 
      SUBROUTINE EIRENE_DEALLOC_EIRBRA
 
      IF (.NOT.ALLOCATED(SNI)) RETURN
 
      DEALLOCATE (SNI)
      DEALLOCATE (SMO)
      DEALLOCATE (SEE)
      DEALLOCATE (SEI)
 
      DEALLOCATE (VOLSUMN)
      DEALLOCATE (VOLSUMM)
      DEALLOCATE (VOLSUMEI)
      DEALLOCATE (VOLSUMEE)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_EIRBRA
 
 
      SUBROUTINE EIRENE_INIT_EIRBRA
 
      SNI = 0._DP
      SMO = 0._DP
      SEE = 0._DP
      SEI = 0._DP
 
      VOLSUMN  = 0._DP
      VOLSUMM  = 0._DP
      VOLSUMEI = 0._DP
      VOLSUMEE = 0._DP
 
      RETURN
      END SUBROUTINE EIRENE_INIT_EIRBRA
 
      END MODULE EIRMOD_EIRBRA
 
 
