!PB 02.03.06: NLRAY ADDED, 
!             NLRAY=.TRUE. : store trajectories and use ray-tracing for this stratum
      MODULE COMSOU

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: ALLOC_COMSOU, DEALLOC_COMSOU, INIT_COMSOU

      REAL(DP), PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     R        RCMSOU(:,:)

      REAL(DP), PUBLIC, POINTER, SAVE ::
     R FLUX(:),     SCALV(:),    RAYFRAC(:),
     R SORENI(:),   SORENE(:),
     R SORVDX(:),   SORVDY(:),   SORVDZ(:),
     R SORCOS(:),   SORMAX(:),
     R SORCTX(:),   SORCTY(:),   SORCTZ(:),
     R SORWGT(:,:), SOREXP(:,:),
     R SORLIM(:,:), SORIND(:,:), SORIFL(:,:),
     R SORAD1(:,:), SORAD2(:,:), SORAD3(:,:),
     R SORAD4(:,:), SORAD5(:,:), SORAD6(:,:)

      REAL(DP), PUBLIC, ALLOCATABLE, SAVE ::
     R SREC(:,:),   EIO(:,:),    EEL(:,:),   MOM(:,:)

      INTEGER, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     I         ICMSOU(:,:)

      INTEGER, PUBLIC, POINTER, SAVE ::
     I IVLSF(:),   ISCLS(:),   ISCLT(:),   ISCL1(:),  ISCL2(:),
     I ISCL3(:),   ISCLB(:),   ISCLA(:),
     I NRSOR(:,:), NPSOR(:,:), NTSOR(:,:),
     I NBSOR(:,:), NASOR(:,:), NISOR(:,:),
     I INDIM(:,:), INSOR(:,:), ISTOR(:,:),
     I NSPEZ(:),   NPTS(:),    NINITL(:),  NEMODS(:),
     I NAMODS(:),  NSRFSI(:),  NRAYEN(:)

      INTEGER, PUBLIC, ALLOCATABLE, SAVE ::
     I INGRDA(:,:,:), INGRDE(:,:,:)

      INTEGER, PUBLIC, SAVE ::
     I NSTRAI, NOMSOU, MOMSOU, LOMSOU

      LOGICAL, PUBLIC, TARGET, ALLOCATABLE, SAVE ::
     L         LCMSOU(:,:)

      LOGICAL, PUBLIC, POINTER, SAVE ::
     L NLPNT(:),  NLLNE(:),  NLSRF(:),  NLVOL(:),  NLCNS(:),
     L NLMOL(:),  NLATM(:),  NLION(:),  NLPLS(:),  NLPHOT(:),
     L NLAVRP(:), NLAVRT(:), NLSRON(:), NLRAY(:)

      LOGICAL, PUBLIC, ALLOCATABLE, SAVE ::
     L NLSYMP(:), NLSYMT(:)


      CONTAINS


      SUBROUTINE ALLOC_COMSOU (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        IF (ALLOCATED(RCMSOU)) RETURN

        NOMSOU=11*NSTRA*NSRFS+13*NSTRA
        MOMSOU=9*NSTRA*NSRFS+15*NSTRA
        LOMSOU=14*NSTRA

        ALLOCATE (RCMSOU(13+11*NSRFS,NSTRA))
        ALLOCATE (ICMSOU(15+9*NSRFS,NSTRA))
        ALLOCATE (LCMSOU(14,NSTRA))

        ALLOCATE (INGRDA(NSRFS,NSTRA,3))
        ALLOCATE (INGRDE(NSRFS,NSTRA,3))

        ALLOCATE (NLSYMP(0:NSTRA))
        ALLOCATE (NLSYMT(0:NSTRA))

        WRITE (55,'(A,T25,I15)')
     .        ' COMSOU ',NOMSOU*8 + (MOMSOU+NSRFS*NSTRA*6)*4 +
     .                  (LOMSOU+2*(NSTRA+1))*4

        FLUX   => RCMSOU( 1,:)
        SCALV  => RCMSOU( 2,:)
        SORENI => RCMSOU( 3,:)
        SORENE => RCMSOU( 4,:)
        SORVDX => RCMSOU( 5,:)
        SORVDY => RCMSOU( 6,:)
        SORVDZ => RCMSOU( 7,:)
        SORCOS => RCMSOU( 8,:)
        SORMAX => RCMSOU( 9,:)
        SORCTX => RCMSOU(10,:)
        SORCTY => RCMSOU(11,:)
        SORCTZ => RCMSOU(12,:)
        RAYFRAC => RCMSOU(13,:)
        SORWGT => RCMSOU(14+ 0*NSRFS : 13+ 1*NSRFS,:)
        SOREXP => RCMSOU(14+ 1*NSRFS : 13+ 2*NSRFS,:)
        SORLIM => RCMSOU(14+ 2*NSRFS : 13+ 3*NSRFS,:)
        SORIND => RCMSOU(14+ 3*NSRFS : 13+ 4*NSRFS,:)
        SORIFL => RCMSOU(14+ 4*NSRFS : 13+ 5*NSRFS,:)
        SORAD1 => RCMSOU(14+ 5*NSRFS : 13+ 6*NSRFS,:)
        SORAD2 => RCMSOU(14+ 6*NSRFS : 13+ 7*NSRFS,:)
        SORAD3 => RCMSOU(14+ 7*NSRFS : 13+ 8*NSRFS,:)
        SORAD4 => RCMSOU(14+ 8*NSRFS : 13+ 9*NSRFS,:)
        SORAD5 => RCMSOU(14+ 9*NSRFS : 13+10*NSRFS,:)
        SORAD6 => RCMSOU(14+10*NSRFS : 13+11*NSRFS,:)

        IVLSF  => ICMSOU( 1,:)
        ISCLS  => ICMSOU( 2,:)
        ISCLT  => ICMSOU( 3,:)
        ISCL1  => ICMSOU( 4,:)
        ISCL2  => ICMSOU( 5,:)
        ISCL3  => ICMSOU( 6,:)
        ISCLB  => ICMSOU( 7,:)
        ISCLA  => ICMSOU( 8,:)
        NSPEZ  => ICMSOU( 9,:)
        NPTS   => ICMSOU(10,:)
        NINITL => ICMSOU(11,:)
        NEMODS => ICMSOU(12,:)
        NAMODS => ICMSOU(13,:)
        NSRFSI => ICMSOU(14,:)
        NRAYEN => ICMSOU(15,:)
        NRSOR  => ICMSOU(16+ 0*NSRFS : 15+ 1*NSRFS,:)
        NPSOR  => ICMSOU(16+ 1*NSRFS : 15+ 2*NSRFS,:)
        NTSOR  => ICMSOU(16+ 2*NSRFS : 15+ 3*NSRFS,:)
        NBSOR  => ICMSOU(16+ 3*NSRFS : 15+ 4*NSRFS,:)
        NASOR  => ICMSOU(16+ 4*NSRFS : 15+ 5*NSRFS,:)
        NISOR  => ICMSOU(16+ 5*NSRFS : 15+ 6*NSRFS,:)
        INDIM  => ICMSOU(16+ 6*NSRFS : 15+ 7*NSRFS,:)
        INSOR  => ICMSOU(16+ 7*NSRFS : 15+ 8*NSRFS,:)
        ISTOR  => ICMSOU(16+ 8*NSRFS : 15+ 9*NSRFS,:)

        NLPNT  => LCMSOU( 1,:)
        NLLNE  => LCMSOU( 2,:)
        NLSRF  => LCMSOU( 3,:)
        NLVOL  => LCMSOU( 4,:)
        NLCNS  => LCMSOU( 5,:)
        NLMOL  => LCMSOU( 6,:)
        NLATM  => LCMSOU( 7,:)
        NLION  => LCMSOU( 8,:)
        NLPLS  => LCMSOU( 9,:)
        NLPHOT => LCMSOU(10,:)
        NLAVRP => LCMSOU(11,:)
        NLAVRT => LCMSOU(12,:)
        NLSRON => LCMSOU(13,:)
        NLRAY  => LCMSOU(14,:)

      ELSE IF (ICAL == 2) THEN

        IF (ALLOCATED(SREC)) RETURN

        ALLOCATE (SREC(0:NPLS,0:NREC))
        ALLOCATE (EIO(0:NPLS,0:NREC))
        ALLOCATE (EEL(0:NPLS,0:NREC))
        ALLOCATE (MOM(0:NPLS,0:NREC))

      END IF

      CALL INIT_COMSOU(ICAL)

      RETURN
      END SUBROUTINE ALLOC_COMSOU


      SUBROUTINE DEALLOC_COMSOU

      IF (ALLOCATED(RCMSOU)) THEN

        DEALLOCATE (RCMSOU)
        DEALLOCATE (ICMSOU)
        DEALLOCATE (LCMSOU)

        DEALLOCATE (INGRDA)
        DEALLOCATE (INGRDE)

        DEALLOCATE (NLSYMP)
        DEALLOCATE (NLSYMT)
      END IF

      IF (ALLOCATED(SREC)) THEN
        DEALLOCATE (SREC)
        DEALLOCATE (EIO)
        DEALLOCATE (EEL)
        DEALLOCATE (MOM)
      END IF

      RETURN
      END SUBROUTINE DEALLOC_COMSOU


      SUBROUTINE INIT_COMSOU (ICAL)

      INTEGER, INTENT(IN) :: ICAL

      IF (ICAL == 1) THEN

        RCMSOU = 0._DP
        ICMSOU = 0
        LCMSOU = .FALSE.
        NLSRON = .TRUE.

        INGRDA = 0
        INGRDE = 0

        NLSYMP = .FALSE.
        NLSYMT = .FALSE.

      ELSE IF (ICAL == 2) THEN

        SREC   = 0._DP
        EIO    = 0._DP
        EEL    = 0._DP
        MOM    = 0._DP

      END IF

      RETURN
      END SUBROUTINE INIT_COMSOU

      END MODULE COMSOU
