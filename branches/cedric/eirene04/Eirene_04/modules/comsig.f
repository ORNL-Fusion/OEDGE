      MODULE COMSIG

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_COMSIG, DEALLOC_COMSIG, INIT_COMSIG

      INTEGER, PUBLIC, SAVE ::
     P         NCMSIG, MCMSIG

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE :: RCMSIG(:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R        FUFFER(:,:), ENERGY(:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R        XCHORD(:), YCHORD(:), ZCHORD(:),
     R        XPIVOT(:), YPIVOT(:), ZPIVOT(:),
     R        TILINE(:), TINP(:),   EMIN1(:),  EMAX1(:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE :: ICMSIG(:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I         IPIVOT(:), ICHORD(:),
     I         NSPSTR(:), NSPSPZ(:),
     I         NSPBLC(:), NSPADD(:),
     I         NCHTAL(:), NSPINI(:), NSPEND(:),
     I         NSPSCL(:), NSPNEW(:)

      INTEGER, PUBLIC, SAVE ::
     I         NCHORI,    NCHORD,    NCHENI


      CONTAINS

      SUBROUTINE ALLOC_COMSIG

      IF (ALLOCATED(RCMSIG)) RETURN

      NCMSIG=10*NCHOR
      MCMSIG=11*NCHOR

      ALLOCATE (RCMSIG(NCMSIG))
      ALLOCATE (FUFFER(NCHOR,NCHEN))
      ALLOCATE (ENERGY(NCHEN))
      ALLOCATE (ICMSIG(MCMSIG))

      WRITE (55,'(A,T25,I15)')
     .      ' COMSIG ',(NCMSIG+(NCHOR+1)*NCHEN)+8 + MCMSIG*4

      XCHORD => RCMSIG( 0*NCHOR+1 :  1*NCHOR)
      YCHORD => RCMSIG( 1*NCHOR+1 :  2*NCHOR)
      ZCHORD => RCMSIG( 2*NCHOR+1 :  3*NCHOR)
      XPIVOT => RCMSIG( 3*NCHOR+1 :  4*NCHOR)
      YPIVOT => RCMSIG( 4*NCHOR+1 :  5*NCHOR)
      ZPIVOT => RCMSIG( 5*NCHOR+1 :  6*NCHOR)
      TILINE => RCMSIG( 6*NCHOR+1 :  7*NCHOR)
      TINP   => RCMSIG( 7*NCHOR+1 :  8*NCHOR)
      EMIN1  => RCMSIG( 8*NCHOR+1 :  9*NCHOR)
      EMAX1  => RCMSIG( 9*NCHOR+1 : 10*NCHOR)

      IPIVOT => ICMSIG( 0*NCHOR+1 :  1*NCHOR)
      ICHORD => ICMSIG( 1*NCHOR+1 :  2*NCHOR)
      NSPSTR => ICMSIG( 2*NCHOR+1 :  3*NCHOR)
      NSPSPZ => ICMSIG( 3*NCHOR+1 :  4*NCHOR)
      NSPBLC => ICMSIG( 4*NCHOR+1 :  5*NCHOR)
      NSPADD => ICMSIG( 5*NCHOR+1 :  6*NCHOR)
      NCHTAL => ICMSIG( 6*NCHOR+1 :  7*NCHOR)
      NSPINI => ICMSIG( 7*NCHOR+1 :  8*NCHOR)
      NSPEND => ICMSIG( 8*NCHOR+1 :  9*NCHOR)
      NSPSCL => ICMSIG( 9*NCHOR+1 : 10*NCHOR)
      NSPNEW => ICMSIG(10*NCHOR+1 : 11*NCHOR)

      RETURN
      END SUBROUTINE ALLOC_COMSIG


      SUBROUTINE DEALLOC_COMSIG

      IF (.NOT.ALLOCATED(RCMSIG)) RETURN

      DEALLOCATE (RCMSIG)
      DEALLOCATE (FUFFER)
      DEALLOCATE (ENERGY)
      DEALLOCATE (ICMSIG)

      RETURN
      END SUBROUTINE DEALLOC_COMSIG


      SUBROUTINE INIT_COMSIG

      RCMSIG = 0._DP
      FUFFER = 0._DP
      ENERGY = 0._DP
      ICMSIG = 0
C
C  SET DEFAULTS FOR LINE OF SIGHT INTEGRATION (BLOCK 12)
C
      NSPBLC = 1
      NSPADD = 0
      NSPNEW = 0

      RETURN
      END SUBROUTINE INIT_COMSIG

      END MODULE COMSIG