      MODULE CREF

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_CREF, DEALLOC_CREF, INIT_CREF

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCREF(:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICREF(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R RPROB0,    ERMIN,     ERCUT,
     R ENAR(:),   DENAR(:),  WIAR(:),   DWIAR(:),  RAAR(:),
     R DRAAR(:),
     R TM(:),     TC(:),     WM(:),     WC(:),    ERDC(:),  HFTR3F(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R RINTEG(:), EINTEG(:), AINTEG(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I INE, INEM, INW, INWM, INR, INRM, NFLR

      INTEGER, PUBLIC, SAVE :: NCREF, MCREF

      LOGICAL, PUBLIC, SAVE :: LTRMOL

      CHARACTER(72), PUBLIC, ALLOCATABLE, SAVE :: REFFIL(:)


      CONTAINS


      SUBROUTINE ALLOC_CREF

      IF (ALLOCATED(RCREF)) RETURN
      
      NCREF = 3+12+11+7+6+5+4+5*NHD6+NHD5
      MCREF  = 7

      ALLOCATE (RCREF(NCREF))
      ALLOCATE (ICREF(MCREF))
      ALLOCATE (REFFIL(NHD6))

      ALLOCATE (RINTEG(0:NLIMPS))
      ALLOCATE (EINTEG(0:NLIMPS))
      ALLOCATE (AINTEG(0:NLIMPS))

      WRITE(55,*) ' CREF ',(NCREF+3*(NLIMPS+1))*8 + MCREF*4 + NHD6*72

      RPROB0    => RCREF(1)
      ERMIN     => RCREF(2)
      ERCUT     => RCREF(3)
      ENAR      => RCREF(4:15)
      DENAR     => RCREF(16:26)
      WIAR      => RCREF(27:33)
      DWIAR     => RCREF(34:39)
      RAAR      => RCREF(40:44)
      DRAAR     => RCREF(45:48)
      TM        => RCREF(49+0*NHD6 : 48+1*NHD6)
      TC        => RCREF(49+1*NHD6 : 48+2*NHD6)
      WM        => RCREF(49+2*NHD6 : 48+3*NHD6)
      WC        => RCREF(49+3*NHD6 : 48+4*NHD6)
      ERDC      => RCREF(49+4*NHD6 : 48+5*NHD6)
      HFTR3F    => RCREF(49+5*NHD6 : 48+5*NHD6+NHD5)

      INE  => ICREF(1)
      INEM => ICREF(2)
      INW  => ICREF(3)
      INWM => ICREF(4)
      INR  => ICREF(5)
      INRM => ICREF(6)
      NFLR => ICREF(7)

      CALL INIT_CREF

      RETURN
      END SUBROUTINE ALLOC_CREF


      SUBROUTINE DEALLOC_CREF

      IF (.NOT.ALLOCATED(RCREF)) RETURN

      DEALLOCATE (RCREF)
      DEALLOCATE (ICREF)
      DEALLOCATE (REFFIL)

      DEALLOCATE (RINTEG)
      DEALLOCATE (EINTEG)
      DEALLOCATE (AINTEG)

      RETURN
      END SUBROUTINE DEALLOC_CREF


      SUBROUTINE INIT_CREF

      RCREF  = 0._DP
      ICREF  = 0
      REFFIL = ' '

      RINTEG = 0._DP
      EINTEG = 0._DP
      AINTEG = 0._DP

      RETURN
      END SUBROUTINE INIT_CREF


      END MODULE CREF