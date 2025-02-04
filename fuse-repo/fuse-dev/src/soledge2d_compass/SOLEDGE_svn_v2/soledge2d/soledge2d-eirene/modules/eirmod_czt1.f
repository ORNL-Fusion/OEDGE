      MODULE EIRMOD_CZT1
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
 
      IMPLICIT NONE
 
      PRIVATE
 
      PUBLIC :: EIRENE_ALLOC_CZT1, EIRENE_DEALLOC_CZT1, EIRENE_INIT_CZT1
 
      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCZT1(:), RCZT2(:)
 
      REAL(DP), PUBLIC, POINTER, SAVE ::
     R RSQDVI(:), CVRSSI(:), ALMASI(:),
     R RSQDVP(:), CVRSSP(:), ALMASP(:),
     R RSQDVA(:), CVRSSA(:), ALMASA(:),
     R RSQDVM(:), CVRSSM(:), ALMASM(:),
     R DEFCX(:),  EEFCX(:),
     R DEFEL(:),  EEFEL(:),
     R DEFPI(:),  EEFPI(:),
     R DENE,      DENI
 
      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R ZT1(:,:),  ZRG(:,:)
 
      INTEGER, PUBLIC, SAVE :: NZT1, NZT2
 
 
      CONTAINS
 
 
      SUBROUTINE EIRENE_ALLOC_CZT1(ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
      INTEGER IND
 
      IF (ICAL == 1) THEN
 
        IF (ALLOCATED(RCZT1)) RETURN
 
        NZT1=2+NION*3+NPLS*3+NATM*3+NMOL*3
 
        ALLOCATE (RCZT1(NZT1))
        ALLOCATE (ZT1(NPLS,NRAD))
        ALLOCATE (ZRG(NPLS,NRAD))
 
        WRITE (55,'(A,T25,I15)')
     .        ' CZT1(1) ',(NZT1+2*NPLS*NRAD)*8
 
        IND = 0
        RSQDVI => RCZT1(1+IND+0*NION : IND+1*NION)
        CVRSSI => RCZT1(1+IND+1*NION : IND+2*NION)
        ALMASI => RCZT1(1+IND+2*NION : IND+3*NION)
 
        IND = IND+3*NION
        RSQDVP => RCZT1(1+IND+0*NPLS : IND+1*NPLS)
        CVRSSP => RCZT1(1+IND+1*NPLS : IND+2*NPLS)
        ALMASP => RCZT1(1+IND+2*NPLS : IND+3*NPLS)
 
        IND = IND+3*NPLS
        RSQDVA => RCZT1(1+IND+0*NATM : IND+1*NATM)
        CVRSSA => RCZT1(1+IND+1*NATM : IND+2*NATM)
        ALMASA => RCZT1(1+IND+2*NATM : IND+3*NATM)
 
        IND = IND+3*NATM
        RSQDVM => RCZT1(1+IND+0*NMOL : IND+1*NMOL)
        CVRSSM => RCZT1(1+IND+1*NMOL : IND+2*NMOL)
        ALMASM => RCZT1(1+IND+2*NMOL : IND+3*NMOL)
 
        IND = IND+3*NMOL
        DENE   => RCZT1(1+IND)
        DENI   => RCZT1(2+IND)
 
      ELSE IF (ICAL == 2) THEN
 
        IF (ALLOCATED(RCZT2)) RETURN
 
        NZT2=2*NREL+2*NRCX+2*NRPI
 
        ALLOCATE (RCZT2(NZT2))
 
        WRITE (55,'(A,T25,I15)')
     .         ' CZT1(2) ',NZT2*8
 
        IND = 0
        DEFCX  => RCZT2(1+IND+0*NRCX : IND+1*NRCX)
        EEFCX  => RCZT2(1+IND+1*NRCX : IND+2*NRCX)
 
        IND = IND+2*NRCX
        DEFEL  => RCZT2(1+IND+0*NREL : IND+1*NREL)
        EEFEL  => RCZT2(1+IND+1*NREL : IND+2*NREL)
 
        IND = IND+2*NREL
        DEFPI  => RCZT2(1+IND+0*NRPI : IND+1*NRPI)
        EEFPI  => RCZT2(1+IND+1*NRPI : IND+2*NRPI)
 
      END IF
 
      CALL EIRENE_INIT_CZT1(ICAL)
 
      RETURN
      END SUBROUTINE EIRENE_ALLOC_CZT1
 
 
      SUBROUTINE EIRENE_DEALLOC_CZT1
 
      IF (.NOT.ALLOCATED(RCZT1)) RETURN
 
      DEALLOCATE (RCZT1)
      DEALLOCATE (RCZT2)
      DEALLOCATE (ZT1)
      DEALLOCATE (ZRG)
 
      RETURN
      END SUBROUTINE EIRENE_DEALLOC_CZT1
 
 
      SUBROUTINE EIRENE_INIT_CZT1(ICAL)
 
      INTEGER, INTENT(IN) :: ICAL
 
      IF (ICAL == 1) THEN
 
        RCZT1 = 0._DP
        ZT1   = 0._DP
        ZRG   = 0._DP
 
      ELSE IF (ICAL == 2) THEN
 
        RCZT2 = 0._DP
 
      END IF
 
      RETURN
      END SUBROUTINE EIRENE_INIT_CZT1
 
      END MODULE EIRMOD_CZT1
 
 
 
 
 
