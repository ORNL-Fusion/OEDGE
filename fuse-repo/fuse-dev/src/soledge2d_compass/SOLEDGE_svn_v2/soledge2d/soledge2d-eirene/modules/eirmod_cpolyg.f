      MODULE EIRMOD_CPOLYG
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CPOLYG, EIRENE_DEALLOC_CPOLYG, 
     P          EIRENE_INIT_CPOLYG
 
      INTEGER, PUBLIC, SAVE :: NCPLYG, NCPLY2, MCPLYG
 
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RCPLYG(:,:), RCPLY2(:)
 
      REAL(DP), PUBLIC, POINTER, SAVE ::
C NCPLYG, REAL
     R VPLX(:,:),  VPLY(:,:),
     R VVTX(:,:),  VVTY(:,:),
     R PLNX(:,:),  PLNY(:,:),
     R BGL(:,:),   BGLP(:,:),
     R PPLNX(:,:), PPLNY(:,:),
     R XPCOR,  YPCOR,  ZPCOR,  PLREFL
 
      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICPLYG(:)
 
      INTEGER, PUBLIC, POINTER, SAVE ::
     I NRPLG,  NPPLG
 
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_CPOLYG
 
      IF (ALLOCATED(RCPLYG)) RETURN
 
      NCPLYG = N1STS*10*N2NDPLG
      NCPLY2 = 4
      MCPLYG = 2
 
      ALLOCATE (RCPLYG(10*N1STS,N2NDPLG))
      ALLOCATE (RCPLY2(NCPLY2))
      ALLOCATE (ICPLYG(MCPLYG))
 
      WRITE (55,'(A,T25,I15)')
     .       ' CPOLYG ',(NCPLYG+NCPLY2)*8 + MCPLYG*4
 
C NCPLYG, REAL
      VPLX  => RCPLYG(1+0*N1STS :  1*N1STS,:)
      VPLY  => RCPLYG(1+1*N1STS :  2*N1STS,:)
      VVTX  => RCPLYG(1+2*N1STS :  3*N1STS,:)
      VVTY  => RCPLYG(1+3*N1STS :  4*N1STS,:)
      PLNX  => RCPLYG(1+4*N1STS :  5*N1STS,:)
      PLNY  => RCPLYG(1+5*N1STS :  6*N1STS,:)
      BGL   => RCPLYG(1+6*N1STS :  7*N1STS,:)
      BGLP  => RCPLYG(1+7*N1STS :  8*N1STS,:)
      PPLNX => RCPLYG(1+8*N1STS :  9*N1STS,:)
      PPLNY => RCPLYG(1+9*N1STS : 10*N1STS,:)
 
      XPCOR  => RCPLY2(1)
      YPCOR  => RCPLY2(2)
      ZPCOR  => RCPLY2(3)
      PLREFL => RCPLY2(4)
 
C MCPLYG, INTEGER
      NRPLG => ICPLYG(1)
      NPPLG => ICPLYG(2)
 
      CALL EIRENE_INIT_CPOLYG
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_CPOLYG
 
 
      SUBROUTINE EIRENE_DEALLOC_CPOLYG
 
      IF (.NOT.ALLOCATED(RCPLYG)) RETURN
 
      DEALLOCATE (RCPLYG)
      DEALLOCATE (RCPLY2)
      DEALLOCATE (ICPLYG)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CPOLYG
 
 
      SUBROUTINE EIRENE_INIT_CPOLYG
 
      RCPLYG = 0._DP
      RCPLY2 = 0._DP
      ICPLYG = 0
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CPOLYG
 
      END MODULE EIRMOD_CPOLYG
