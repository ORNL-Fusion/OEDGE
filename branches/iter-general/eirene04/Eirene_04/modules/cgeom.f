      MODULE CGEOM

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CGEOM, DEALLOC_CGEOM, INIT_CGEOM,
     P          CELL_ELEM, CELL_LIST

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: 
     R        RCGM1(:), RCGM2(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R VOLADD(:), VOLCOR(:), VOLG(:), VOLTAL(:), VOLTOT,
     R AREA(:),   AREAG(:),  CELDIA(:), XCOM(:), YCOM(:),
     R XPOINT(:), YPOINT(:) 

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R XPOL(:,:), YPOL(:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NPOINT(:,:),   NSTGRD(:), NGHPLS(:,:,:),
     I NGHPOL(:,:,:), NCLTAL(:), INDPOINT(:,:)

      INTEGER, PUBLIC, SAVE :: NCGM1, NCGM2, NNODES

      TYPE :: CELL_ELEM
        INTEGER :: NOCELL
        TYPE(CELL_ELEM), POINTER :: NEXT_CELL
      END TYPE CELL_ELEM

      TYPE :: CELL_LIST
        TYPE(CELL_ELEM), POINTER :: PCELL
      END TYPE CELL_LIST
        
      TYPE(CELL_LIST), ALLOCATABLE, SAVE, PUBLIC :: COORCELL(:)


      CONTAINS


      SUBROUTINE ALLOC_CGEOM

      IF (ALLOCATED(RCGM1)) RETURN

      NCGM1 = NADD+NBMAX+7*NRAD+NRTAL+1+NLMPGS
      NCGM2 = 2*N1STS*N2NDPLG

      ALLOCATE (RCGM1(NCGM1))
      ALLOCATE (RCGM2(2*N1STS,N2NDPLG))

      ALLOCATE (NPOINT(2,NPPART))
      ALLOCATE (NSTGRD(NRAD))
      ALLOCATE (NGHPLS(4,N1STS,N2NDPLG))
      ALLOCATE (NGHPOL(4,N1STS,N2NDPLG))
      ALLOCATE (NCLTAL(NRAD))
      ALLOCATE (INDPOINT(N1STS,N2NDPLG))
      
      ALLOCATE (COORCELL(NRAD))

      WRITE (55,'(A,T25,I15)')
     .      ' CGEOM ',(NCGM1+2*N1STS*N2NDPLG)*8 +
     .                (2*NPPART+2*NRAD+8*N1STS*N2NDPLG)*4

      VOLADD => RCGM1(1 : NADD)
      VOLTAL => RCGM1(1+NADD : NADD+NRTAL)
      VOLCOR => RCGM1(1+NADD+NRTAL : NADD+NRTAL+NBMAX)
      VOLG   => RCGM1(1+NADD+NRTAL+NBMAX : NADD+NRTAL+NBMAX+NRAD)
      AREA   => RCGM1(1+NADD+NRTAL+NBMAX+NRAD : NADD+NRTAL+NBMAX+2*NRAD)
      CELDIA => RCGM1(1+NADD+NRTAL+NBMAX+2*NRAD :
     .                  NADD+NRTAL+NBMAX+3*NRAD)
      XCOM   => RCGM1(1+NADD+NRTAL+NBMAX+3*NRAD :
     .                  NADD+NRTAL+NBMAX+4*NRAD)
      YCOM   => RCGM1(1+NADD+NRTAL+NBMAX+4*NRAD :
     .                  NADD+NRTAL+NBMAX+5*NRAD)
      XPOINT => RCGM1(1+NADD+NRTAL+NBMAX+5*NRAD :
     .                  NADD+NRTAL+NBMAX+6*NRAD)
      YPOINT => RCGM1(1+NADD+NRTAL+NBMAX+6*NRAD :
     .                  NADD+NRTAL+NBMAX+7*NRAD)
      AREAG  => RCGM1(1+NADD+NRTAL+NBMAX+7*NRAD :
     .                  NADD+NRTAL+NBMAX+7*NRAD+NLMPGS)
      VOLTOT => RCGM1(1+NADD+NRTAL+NBMAX+7*NRAD+NLMPGS)

      XPOL => RCGM2(1+0*N1STS:1*N1STS,:)
      YPOL => RCGM2(1+1*N1STS:2*N1STS,:)

      CALL INIT_CGEOM

      RETURN
      END SUBROUTINE ALLOC_CGEOM


      SUBROUTINE DEALLOC_CGEOM

      IF (.NOT.ALLOCATED(RCGM1)) RETURN

      DEALLOCATE (RCGM1)
      DEALLOCATE (RCGM2)

      DEALLOCATE (NPOINT)
      DEALLOCATE (NSTGRD)
      DEALLOCATE (NGHPLS)
      DEALLOCATE (NGHPOL)
      DEALLOCATE (NCLTAL)
      DEALLOCATE (INDPOINT)

      DEALLOCATE (COORCELL)

      RETURN
      END SUBROUTINE DEALLOC_CGEOM


      SUBROUTINE INIT_CGEOM

      INTEGER :: I

      RCGM1    = 0._DP
      RCGM2    = 0._DP

      NPOINT   = 0
      NSTGRD   = 0
      NGHPLS   = 0
      NGHPOL   = 0
      NCLTAL   = 0
      INDPOINT = 0

      DO I=1,NRAD
        NULLIFY (COORCELL(I)%PCELL)
      END DO

      RETURN
      END SUBROUTINE INIT_CGEOM

      END MODULE CGEOM





