      MODULE CSDVI

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CSDVI, DEALLOC_CSDVI, INIT_CSDVI

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        SDVI1(:,:), SDVI2(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R SIGMA(:,:),    SGMS(:),
     R SIGMAW(:,:),   SGMWS(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SIGMAC(:,:,:), SGMCS(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ISDVI(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I IIH(:),    IGH(:),
     I IIHW(:),   IGHW(:),
     I ICOV(:),
     I NSIGI,     NSIGVI,   NSIGSI, NSIGCI, NSIGI_BGK, NSIGI_COP, 
     I NSIGI_SPC,
C SPEED UP OF SUBROUTINE STATIS
     I IMETCL(:), ICLMT(:), NCLMT, NCLMTS

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I IIHC(:,:), IGHC(:,:)

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L LMETSP(:)

      INTEGER, PUBLIC, SAVE ::
     I NSDVI1, NSDVI2, NSDVC1, NSDVC2, NSDVI, MSDVI


      CONTAINS


      SUBROUTINE ALLOC_CSDVI (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(SDVI1)) RETURN

        NSDVI1 = NSD*(NRTAL+1)
        NSDVI2 = NSDW*(NLIMPS+1)
        NSDVC1 = 3*NCV*NRTAL
        NSDVC2 = 3*NCV
        NSDVI  = NSDVI1+NSDVI2+NSDVC1+NSDVC2
        MSDVI  = NSD*2+NSDW*2+NCV+7+2*NRTAL+2
  
  
        ALLOCATE (SDVI1(NSD,NRTAL+1))
        ALLOCATE (SDVI2(NSDW,NLIMPS+1))
        ALLOCATE (SIGMAC(0:2,NCV,NRTAL))
        ALLOCATE (SGMCS(0:2,NCV))
  
        ALLOCATE (ISDVI(MSDVI))
        ALLOCATE (IIHC(2,NCV))
        ALLOCATE (IGHC(2,NCV))
  
        WRITE (55,'(A,T25,I15)')
     .        ' CSDVI ',(NSD*(NRTAL+1) + NSDW*(NLIMPS+1) +
     .                   (NRTAL+1)*3*NCV)*8 +
     .                  (MSDVI+4*NCV)*4

        SIGMA  => SDVI1(:,1:NRTAL)
        SGMS   => SDVI1(:,NRTAL+1)
  
        SIGMAW => SDVI2(:,1:NLIMPS)
        SGMWS  => SDVI2(:,NLIMPS+1)
  
        NSIGI     => ISDVI(1)
        NSIGVI    => ISDVI(2)
        NSIGSI    => ISDVI(3)
        NSIGCI    => ISDVI(4)
        NSIGI_BGK => ISDVI(5)
        NSIGI_COP => ISDVI(6)
        NSIGI_SPC => ISDVI(7)
        NCLMT     => ISDVI(8)
        NCLMTS    => ISDVI(9)
        IIH       => ISDVI(10                  : 9+  NSD)
        IGH       => ISDVI(10+  NSD            : 9+2*NSD)
        IIHW      => ISDVI(10+2*NSD            : 9+2*NSD+  NSDW)
        IGHW      => ISDVI(10+2*NSD+  NSDW     : 9+2*NSD+2*NSDW)
        ICOV      => ISDVI(10+2*NSD+2*NSDW     : 9+2*NSD+2*NSDW+NCV)
        IMETCL    => ISDVI(10+2*NSD+2*NSDW+NCV : 
     .                     9+2*NSD+2*NSDW+NCV+NRTAL)
        ICLMT     => ISDVI(10+2*NSD+2*NSDW+NCV+NRTAL : MSDVI)

      ELSE IF (ICAL == 2) THEN
  
        IF (ALLOCATED(LMETSP)) RETURN

!pb        ALLOCATE (LMETSP(N1MX+NSNV))
        ALLOCATE (LMETSP(NSPZTOT+NSNV))

      END IF

      CALL INIT_CSDVI (ICAL)

      RETURN
      END SUBROUTINE ALLOC_CSDVI


      SUBROUTINE DEALLOC_CSDVI

      IF (.NOT.ALLOCATED(SDVI1)) RETURN

      DEALLOCATE (SDVI1)
      DEALLOCATE (SDVI2)
      DEALLOCATE (SIGMAC)
      DEALLOCATE (SGMCS)

      DEALLOCATE (ISDVI)
      DEALLOCATE (IIHC)
      DEALLOCATE (IGHC)

      DEALLOCATE (LMETSP)

      RETURN
      END SUBROUTINE DEALLOC_CSDVI


      SUBROUTINE INIT_CSDVI (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        SDVI1  = 0._DP
        SDVI2  = 0._DP
        SIGMAC = 0._DP
        SGMCS  = 0._DP
  
        ISDVI  = 0
        IIHC   = 0
        IGHC   = 0

      ELSE IF (ICAL == 2) THEN

        LMETSP = .FALSE.

      END IF

      RETURN
      END SUBROUTINE INIT_CSDVI

      END MODULE CSDVI
