      MODULE CGEOM

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CGEOM, DEALLOC_CGEOM, INIT_CGEOM

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: 
     R        RCGM1(:), RCGM2(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R VOLADD(:), VOLCOR(:), VOLG(:), VOLTAL(:), VOLTOT,
     R AREA(:),   AREAG(:),  CELDIA(:), XCOM(:), YCOM(:) 

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R XPOL(:,:), YPOL(:,:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I NPOINT(:,:),   NSTGRD(:), NGHPLS(:,:,:),
     I NGHPOL(:,:,:), NCLTAL(:)

      INTEGER, PUBLIC, SAVE :: NCGM1, NCGM2

c slmod begin
      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I RVRTAG(:,:), PVRTAG(:,:)
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R XVERT(:,:,:), YVERT(:,:,:)


      INTEGER, PUBLIC, SAVE ::
     I ASCNCELL, ASC3DMODE, ASCCODE, ASCNCUT

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I ASCCELL(:), ASCNVERTEX(:), ASCREGION(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R ASCZMIN3D(:), ASCZMAX3D(:), ASCVERTEX(:,:)


      INTEGER, PUBLIC, SAVE ::
     I EIRNSDTOR, EIRNTRANS, MAXNTOR

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EIRTRANS(:,:)


      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R EIRSDTOR(:)



c slmod end

      CONTAINS


      SUBROUTINE ALLOC_CGEOM

      IF (ALLOCATED(RCGM1)) RETURN

      NCGM1 = NADD+NBMAX+5*NRAD+NRTAL+1+NLMPGS
      NCGM2 = 2*N1STS*N2NDPLG

      ALLOCATE (RCGM1(NCGM1))
      ALLOCATE (RCGM2(2*N1STS,N2NDPLG))

      ALLOCATE (NPOINT(2,NPPART))
      ALLOCATE (NSTGRD(NRAD))
      ALLOCATE (NGHPLS(4,N1STS,N2NDPLG))
      ALLOCATE (NGHPOL(4,N1STS,N2NDPLG))
      ALLOCATE (NCLTAL(NRAD))

      WRITE (55,*) ' CGEOM ',(NCGM1+2*N1STS*N2NDPLG)*8 +
     .                       (2*NPPART+2*NRAD+8*N1STS*N2NDPLG)*4

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
      AREAG  => RCGM1(1+NADD+NRTAL+NBMAX+5*NRAD :
     .                  NADD+NRTAL+NBMAX+5*NRAD+NLMPGS)
      VOLTOT => RCGM1(1+NADD+NRTAL+NBMAX+5*NRAD+NLMPGS)

      XPOL => RCGM2(1+0*N1STS:1*N1STS,:)
      YPOL => RCGM2(1+1*N1STS:2*N1STS,:)
c slmod begin
c
c
c


      ALLOCATE (RVRTAG(N1ST,N2ND))
      ALLOCATE (PVRTAG(N1ST,N2ND))

      ALLOCATE (XVERT(N1ST,N2ND,2))
      ALLOCATE (YVERT(N1ST,N2ND,2))


      ALLOCATE (ASCCELL(NLIM))
      ALLOCATE (ASCNVERTEX(NLIM))
      ALLOCATE (ASCREGION(NLIM))

      ALLOCATE (ASCZMIN3D(NLIM))
      ALLOCATE (ASCZMAX3D(NLIM))
      ALLOCATE (ASCVERTEX(40,NLIM))

      ALLOCATE (EIRTRANS(50,3))

      MAXNTOR=100
      ALLOCATE (EIRSDTOR(MAXNTOR))
c slmod end
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

      RETURN
      END SUBROUTINE DEALLOC_CGEOM


      SUBROUTINE INIT_CGEOM

      RCGM1  = 0._DP
      RCGM2  = 0._DP

      NPOINT = 0
      NSTGRD = 0
      NGHPLS = 0
      NGHPOL = 0
      NCLTAL = 0

      RETURN
      END SUBROUTINE INIT_CGEOM

      END MODULE CGEOM