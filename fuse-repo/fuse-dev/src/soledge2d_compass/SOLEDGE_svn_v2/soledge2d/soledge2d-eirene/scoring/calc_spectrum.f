!pb  24.04.07:  allow for logarithmic equidistant energy bins
 
      SUBROUTINE EIRENE_CALC_SPECTRUM (WT,IND,ISC)
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CESTIM
      USE EIRMOD_COMPRT
      USE EIRMOD_CUPD
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_CCONA
      USE EIRMOD_COMUSR
 
      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: IND, ISC
      REAL(DP), INTENT(IN) :: WT
      INTEGER :: ISPC, I, IS, IC, IRDO, IRD, IAT, IML, IIO, IPL
      REAL(DP) :: ADD, WV, DIST, WTR, SPCVX, SPCVY, SPCVZ, CDYN, EB
      REAL(DP), ALLOCATABLE, SAVE :: CNDYNA(:), CNDYNM(:), CNDYNI(:),
     .                               CNDYNP(:)
      TYPE(EIRENE_SPECTRUM), POINTER :: P
 
      IF ((ISC == 0) .AND. (IND .NE. 1)) RETURN
 
      SELECT CASE (ITYP)
      CASE (0)
        IS = IPHOT
        CDYN = 1._DP
      CASE (1)
        IS = IATM
        IF (.NOT.ALLOCATED(CNDYNA)) THEN
          ALLOCATE (CNDYNA(NATM))
          DO IAT=1,NATMI
            CNDYNA(IAT)=AMUA*RMASSA(IAT)
          END DO
        END IF
        CDYN = CNDYNA(IATM)
      CASE (2)
        IS = IMOL
        IF (.NOT.ALLOCATED(CNDYNM)) THEN
          ALLOCATE (CNDYNM(NMOL))
          DO IML=1,NMOLI
            CNDYNM(IML)=AMUA*RMASSM(IML)
          END DO
        END IF
        CDYN = CNDYNM(IMOL)
      CASE (3)
        IS = IION
        IF (.NOT.ALLOCATED(CNDYNI)) THEN
          ALLOCATE (CNDYNI(NION))
          DO IIO=1,NIONI
            CNDYNI(IIO)=AMUA*RMASSI(IIO)
          END DO
        END IF
        CDYN = CNDYNI(IION)
      CASE (4)
        IS = IPLS
        IF (.NOT.ALLOCATED(CNDYNP)) THEN
          ALLOCATE (CNDYNP(NPLS))
          DO IPL=1,NPLSI
            CNDYNP(IPL)=AMUA*RMASSP(IPL)
          END DO
        END IF
        CDYN = CNDYNP(IPLS)
      END SELECT
 
      IF (ISC == 0) THEN    ! SURFACE
 
        DO ISPC=1,NADSPC
          P => ESTIML(ISPC)%PSPC
          IF ((P%ISRFCLL == ISC) .AND.
     .        (P%ISPCSRF == MSURF) .AND.
     .        (P%IPRTYP == ITYP) .AND.
     .        ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN
 
            SELECT CASE(ESTIML(ISPC)%PSPC%ISPCTYP)
            CASE (1)
              ADD = WT
            CASE (2)
              ADD = WT*E0
            CASE DEFAULT
              ADD = 0._DP
            END SELECT
 
            EB = E0
 
            IF (ESTIML(ISPC)%PSPC%LOG) EB=LOG10(EB)
 
            IF (EB < ESTIML(ISPC)%PSPC%SPCMIN) THEN
              I = 0
            ELSEIF (EB >= ESTIML(ISPC)%PSPC%SPCMAX) THEN
              I = ESTIML(ISPC)%PSPC%NSPC + 1
            ELSE
              I = (EB - ESTIML(ISPC)%PSPC%SPCMIN) *
     .             ESTIML(ISPC)%PSPC%SPCDELI + 1
            END IF
            ESTIML(ISPC)%PSPC%SPC(I) = ESTIML(ISPC)%PSPC%SPC(I) + ADD
            ESTIML(ISPC)%PSPC%ESP_MIN= MIN(ESTIML(ISPC)%PSPC%ESP_MIN,EB)
            ESTIML(ISPC)%PSPC%ESP_MAX= MAX(ESTIML(ISPC)%PSPC%ESP_MAX,EB)
            ESTIML(ISPC)%PSPC%IMETSP = 1
          END IF
        END DO
 
      ELSE     ! CELL
 
        WV=WEIGHT/VEL
        DO IC=1,NCOU
          DIST=CLPD(IC)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(IC)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)
 
 
          DO ISPC=1,NADSPC
            P => ESTIML(ISPC)%PSPC
            IF ((P%ISRFCLL > 0) .AND.
     .          (((P%ISRFCLL == 1).AND.(P%ISPCSRF == IRD)) .OR.      ! scoring cell
     .           ((P%ISRFCLL == 2).AND.(P%ISPCSRF == IRDO))) .AND.   ! geometry cell
     .          (P%IPRTYP == ITYP) .AND.
     .          ((P%IPRSP == IS) .OR. (P%IPRSP == 0))) THEN
 
              SELECT CASE(ESTIML(ISPC)%PSPC%ISPCTYP)
              CASE (1)
                ADD = WTR
              CASE (2)
                ADD = WTR*E0
              CASE (3)
                ADD = WTR*VEL*CDYN
              CASE DEFAULT
                ADD = 0._DP
              END SELECT
 
              EB = E0
              IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
                SPCVX = ESTIML(ISPC)%PSPC%SPCVX
                SPCVY = ESTIML(ISPC)%PSPC%SPCVY
                SPCVZ = ESTIML(ISPC)%PSPC%SPCVZ
                EB = EB * (SPCVX*VELX+SPCVY*VELY+SPCVZ*VELZ)
              END IF
 
              IF (ESTIML(ISPC)%PSPC%LOG) EB=LOG10(EB)
 
              IF (EB < ESTIML(ISPC)%PSPC%SPCMIN) THEN
                I = 0
              ELSEIF (EB >= ESTIML(ISPC)%PSPC%SPCMAX) THEN
                I = ESTIML(ISPC)%PSPC%NSPC + 1
              ELSE
                I = (EB - ESTIML(ISPC)%PSPC%SPCMIN) *
     .               ESTIML(ISPC)%PSPC%SPCDELI + 1
              END IF
              ESTIML(ISPC)%PSPC%SPC(I) = ESTIML(ISPC)%PSPC%SPC(I) + ADD
              ESTIML(ISPC)%PSPC%ESP_MIN =
     .               MIN(ESTIML(ISPC)%PSPC%ESP_MIN,EB)
              ESTIML(ISPC)%PSPC%ESP_MAX =
     .               MAX(ESTIML(ISPC)%PSPC%ESP_MAX,EB)
              ESTIML(ISPC)%PSPC%IMETSP = 1
            END IF
          END DO
 
        END DO
 
      END IF
 
      RETURN
      END SUBROUTINE EIRENE_CALC_SPECTRUM
