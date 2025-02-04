      MODULE EIRMOD_BRAEIR
 
C  PLASMA DATA: NI,TE,TI,VV,UU,PR,UP,RR,FNIX,FNIY.. (BRAAMS ---> EIRENE)
 
      USE EIRMOD_PRECISION
!pb      USE PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_BRAEIR, EIRENE_DEALLOC_BRAEIR, 
     P          EIRENE_INIT_BRAEIR
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  DNIB(:,:,:),    TEB(:,:),       TIB(:,:),   VVB(:,:,:),
     R  UUB(:,:,:),     PRB(:,:),       UPB(:,:,:), RRB(:,:),
     R  FNIXB(:,:,:),   FNIYB(:,:,:),   FEIXB(:,:), FEIYB(:,:),
     R  FEEXB(:,:),     FEEYB(:,:),     VOLB(:,:),  BFELDB(:,:),
     R  VPARXB(:,:,:),  VPARYB(:,:,:),  VRADXB(:,:,:), VRADYB(:,:,:),
     R  DELTAE_PARXB(:,:), DELTAE_PARYB(:,:), 
     R  DELTAE_RADXB(:,:), DELTAE_RADYB(:,:),
     R  DELTAI_PARXB(:,:), DELTAI_PARYB(:,:), 
     R  DELTAI_RADXB(:,:), DELTAI_RADYB(:,:),
     R  DELTA_SHEATHXB(:,:), DELTA_SHEATHYB(:,:)

      INTEGER, SAVE :: NDXP, NDYP, NFL
 
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_BRAEIR(NDXD, NDYD, NFLD)
 
      INTEGER, INTENT(IN) :: NDXD, NDYD, NFLD
 
      IF (ALLOCATED(DNIB)) RETURN
 
      NDXP = NDXD+1
      NDYP = NDYD+1
      NFL = NFLD
 
      ALLOCATE (DNIB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (TEB(0:NDXP,0:NDYP))
      ALLOCATE (TIB(0:NDXP,0:NDYP))
      ALLOCATE (VVB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (UUB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (PRB(0:NDXP,0:NDYP))
      ALLOCATE (UPB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (RRB(0:NDXP,0:NDYP))
      ALLOCATE (FNIXB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (FNIYB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (FEIXB(0:NDXP,0:NDYP))
      ALLOCATE (FEIYB(0:NDXP,0:NDYP))
      ALLOCATE (FEEXB(0:NDXP,0:NDYP))
      ALLOCATE (FEEYB(0:NDXP,0:NDYP))
      ALLOCATE (VOLB(0:NDXP,0:NDYP))
      ALLOCATE (BFELDB(0:NDXP,0:NDYP))

      ALLOCATE (VPARXB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (VPARYB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (VRADXB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (VRADYB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (DELTAE_PARXB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAE_PARYB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAE_RADXB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAE_RADYB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAI_PARXB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAI_PARYB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAI_RADXB(0:NDXP,0:NDYP))
      ALLOCATE (DELTAI_RADYB(0:NDXP,0:NDYP))
      ALLOCATE (DELTA_SHEATHXB(0:NDXP,0:NDYP))
      ALLOCATE (DELTA_SHEATHYB(0:NDXP,0:NDYP))
 
      WRITE (55,'(A,T25,I15)')
     .      ' BRAEIR ',(NDXP+1)*(NDYP+1)*(10*NFL+20)*8
 
      CALL EIRENE_INIT_BRAEIR
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_BRAEIR
 
 
      SUBROUTINE EIRENE_DEALLOC_BRAEIR
 
      IF (.NOT.ALLOCATED(DNIB)) RETURN
 
      DEALLOCATE (DNIB)
      DEALLOCATE (TEB)
      DEALLOCATE (TIB)
      DEALLOCATE (VVB)
      DEALLOCATE (UUB)
      DEALLOCATE (PRB)
      DEALLOCATE (UPB)
      DEALLOCATE (RRB)
      DEALLOCATE (FNIXB)
      DEALLOCATE (FNIYB)
      DEALLOCATE (FEIXB)
      DEALLOCATE (FEIYB)
      DEALLOCATE (FEEXB)
      DEALLOCATE (FEEYB)
      DEALLOCATE (VOLB)
      DEALLOCATE (BFELDB)

      DEALLOCATE (VPARXB)
      DEALLOCATE (VPARYB)
      DEALLOCATE (VRADXB)
      DEALLOCATE (VRADYB)
      DEALLOCATE (DELTAE_PARXB)
      DEALLOCATE (DELTAE_PARYB)
      DEALLOCATE (DELTAE_RADXB)
      DEALLOCATE (DELTAE_RADYB)
      DEALLOCATE (DELTAI_PARXB)
      DEALLOCATE (DELTAI_PARYB)
      DEALLOCATE (DELTAI_RADXB)
      DEALLOCATE (DELTAI_RADYB)
      DEALLOCATE (DELTA_SHEATHXB)
      DEALLOCATE (DELTA_SHEATHYB)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_BRAEIR
 
 
      SUBROUTINE EIRENE_INIT_BRAEIR
 
      DNIB    = 0.D0
      TEB     = 0.D0
      TIB     = 0.D0
      VVB     = 0.D0
      UUB     = 0.D0
      PRB     = 0.D0
      UPB     = 0.D0
      RRB     = 0.D0
      FNIXB   = 0.D0
      FNIYB   = 0.D0
      FEIXB   = 0.D0
      FEIYB   = 0.D0
      FEEXB   = 0.D0
      FEEYB   = 0.D0
      VOLB    = 0.D0
      BFELDB  = 0.D0

      VPARXB  = 0.D0
      VPARYB  = 0.D0
      VRADXB  = 0.D0
      VRADYB  = 0.D0
      DELTAE_PARXB = 0.D0
      DELTAE_PARYB = 0.D0
      DELTAE_RADXB = 0.D0
      DELTAE_RADYB = 0.D0
      DELTAI_PARXB = 0.D0
      DELTAI_PARYB = 0.D0
      DELTAI_RADXB = 0.D0
      DELTAI_RADYB = 0.D0
      DELTA_SHEATHXB = 0.D0
      DELTA_SHEATHYB = 0.D0
 
      RETURN
      END SUBROUTINE EIRENE_INIT_BRAEIR
 
      END MODULE EIRMOD_BRAEIR
