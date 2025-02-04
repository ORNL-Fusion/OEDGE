      MODULE BRAEIR

C  PLASMA DATA: NI,TE,TI,VV,UU,PR,UP,RR,FNIX,FNIY.. (BRAAMS ---> EIRENE)

      USE PRECISION
!pb      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_BRAEIR, DEALLOC_BRAEIR, INIT_BRAEIR

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R  DNIB(:,:,:),    TEB(:,:),       TIB(:,:),   VVB(:,:,:),
     R  UUB(:,:,:),     PRB(:,:),       UPB(:,:,:), RRB(:,:),
     R  FNIXB(:,:,:),   FNIYB(:,:,:),   FEIXB(:,:), FEIYB(:,:),
     R  FEEXB(:,:),     FEEYB(:,:),     VOLB(:,:),  BFELDB(:,:),
     R  FNIX_YB(:,:,:), FNIY_XB(:,:,:), 
     R  UUDIAB(:,:,:),  VVDIAB(:,:,:) 

      INTEGER, SAVE :: NDXP, NDYP, NFL


      CONTAINS


      SUBROUTINE ALLOC_BRAEIR(NDXD, NDYD, NFLD)

      INTEGER, INTENT(IN) :: NDXD, NDYD, NFLD

      IF (ALLOCATED(DNIB)) RETURN
c slmod begin
c...  I think the +1 is unnecessary since NDXD is from NDXP, which is
c     set equal to NDX+1 already.  Same for NDYP:
      NDXP = NDXD
      NDYP = NDYD
c
c      NDXP = NDXD+1
c      NDYP = NDYD+1
c slmod end
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
      ALLOCATE (FNIX_YB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (FNIY_XB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (UUDIAB(0:NDXP,0:NDYP,NFL))
      ALLOCATE (VVDIAB(0:NDXP,0:NDYP,NFL))

      WRITE (55,*) ' BRAEIR ',(NDXP+1)*(NDYP+1)*(10*NFL+10)*8

      CALL INIT_BRAEIR

      RETURN
      END SUBROUTINE ALLOC_BRAEIR


      SUBROUTINE DEALLOC_BRAEIR

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
      DEALLOCATE (FNIX_YB)
      DEALLOCATE (FNIY_XB)
      DEALLOCATE (UUDIAB)
      DEALLOCATE (VVDIAB)

      RETURN
      END SUBROUTINE DEALLOC_BRAEIR


      SUBROUTINE INIT_BRAEIR

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
      FNIX_YB = 0.D0
      FNIY_XB = 0.D0
      UUDIAB  = 0.D0
      VVDIAB  = 0.D0

      RETURN
      END SUBROUTINE INIT_BRAEIR

      END MODULE BRAEIR
