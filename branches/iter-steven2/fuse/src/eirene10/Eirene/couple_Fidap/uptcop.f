C
      SUBROUTINE UPTCOP(XSTOR2,XSTORV2,WV,IFLAG)
C
C  UPDATE TALLIES COPV FOR COUPLING TO OTHER CODES
C
C  THIS VERSION: COUPLING TO FIDAP, FOR RADIATION TRANSFER
C
C  COPV(1): TOTAL SAMPLED EMISSION, SAME AS EPPHT IN LOCATE.F
C                                   COLLISION ESTIMATOR
C  COPV(2): TOTAL ABSORPTION, SAME AS EPHPHT IN UPDATE.F
C                             COLLISION ESTIMATOR IN CELL OF EMISSION
C                             TRACKLENGTH ESTIMATOR ELSE
c  copv(1), copv(2) are redundant, in principle, but they (and their
c                   variances), are used in infcop, entry if3cop.
C
C  COPV(3): SAMPLED EMISSION, THICK PART, ESTIMATOR AS COPV(1)
C  COPV(4): ABSORPTION OF THICK PART, ESTIMATORS AS WITH COPV(2)
C  COPV(5): SAMPLED EMISSION, THIN PART, ESTIMATOR AS COPV(2)
C  COPV(6): ABSORPTION OF THIN PART, ESTIMATORS AS WITH COPV(2)
C
C  COPV(7), COPV(8), COPV(9):  VEC-I
C           NOW: SUM COPV(3)+....COPV(6), FOR VARIANCE (STATIS_COP) ON ADDV(5)
C                IN STATIS_COP: THIS SUM IS ON SIGMA_COP(NCPVI+1)
C
C  USER SUPPLIED TRACKLENGTH ESTIMATOR, VOLUME AVERAGED
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE CPOLYG
      USE CGRID
      USE CSPEZ
      USE CZT1
      USE CGEOM
      USE COMPRT
      USE CSDVI
      USE COMXS
      USE PHOTON

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XSTOR2(MSTOR1,MSTOR2,N2ND+N3RD),
     .                        XSTORV2(NSTORV,N2ND+N3RD), WV
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: P, WTRSIG, EION, V0_PARB, PARMOM_0, DIST, WTR, XLCRIT
      INTEGER :: IFIRST, ISP, ICOU, IRDO, IRD, IAOT, IROT, UPDF
      INTEGER, SAVE :: NMTSP
      LOGICAL, SAVE :: LTHICK
CDR
      DATA IFIRST/0/
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
C
        NMTSP=NPHOTI+NATMI+NMOLI+NIONI+NPLSI+NADVI+NALVI
C
      ENDIF

C
C  WV=WEIGHT/VEL
C
C  PHOTONS
      IF (ITYP.EQ.0) THEN

!  EMISSION

        IF (NSTCLL > 0) THEN
          IRD=NCLTAL(NSTCLL)
!  TOTAL SAMPLED EMISSION
          COPV(1,IRD) = COPV(1,IRD) + STEMIS
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+1)=.TRUE.
C
C  COPV(3): THICK FRACTION OF SAMPLED EMISSION
C  COPV(5): THIN  FRACTION OF SAMPLED EMISSION
C

!  ZMFP < ZMFPTHI*XLCRIT  => THICK FRACTION
!  XLCRIT = TEMPERATURE GRADIENT LENGTH (CM) IN FLIGHT DIRECTION

          XLCRIT = TEDTEDX(IRD)*VELX+TEDTEDY(IRD)*VELY+TEDTEDZ(IRD)*VELZ

          IF (ABS(XLCRIT)*ZMFPI > 1._dp/ZMFPTHI) THEN
! XLCRIT/zmfp > ZMFPHT  => THICK FRACTION 
!  EMISSION PART 1
            COPV(3,IRD) = COPV(3,IRD) + STEMIS  
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+3)=.TRUE.
            ISP = 4
            LTHICK = .TRUE.
          ELSE
!  XLCRIT/zmfp < ZMFPHT  => THIN FRACTION 
!  EMISSION PART 2
            COPV(5,IRD) = COPV(5,IRD) + STEMIS  
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+5)=.TRUE.
            ISP = 6
            LTHICK = .FALSE.
          END IF
          IF (IMETCL(IRD) == 0) THEN
            NCLMT = NCLMT+1
            ICLMT(NCLMT) = IRD
            IMETCL(IRD) = NCLMT
          END IF
          NSTCLL = -1
        END IF

!  EMISSION DONE, NOW: ABSORPTION

C  COPV(2): TOTAL ABSORPTION
C  COPV(4): THICK FRACTION OF ABSORPTION
C  COPV(6): THIN  FRACTION OF ABSORPTION
C
        DO ICOU=1,NCOU
          DIST=CLPD(ICOU)
          WTR=WV*DIST
          IRDO=NRCELL+NUPC(ICOU)*NR1P2+NBLCKA
          IRD=NCLTAL(IRDO)

!  xlcrit = temperature gradient length (cm) in flight direction
          XLCRIT = TEDTEDX(IRD)*VELX+TEDTEDY(IRD)*VELY+TEDTEDZ(IRD)*VELZ
          IF (ABS(XLCRIT)*ZMFPI > 1._dp/ZMFPTHI) THEN
            ISP = 4
            LTHICK = .TRUE.
          ELSE
            ISP = 6
            LTHICK = .FALSE.
          END IF

          IF (IMETCL(IRD) == 0) THEN
            NCLMT = NCLMT+1
            ICLMT(NCLMT) = IRD
            IMETCL(IRD) = NCLMT
          END IF
C
          IF (LGVAC(IRDO,0)) CYCLE
C
          if (ncou.gt.1) then
          XSTOR(1:mstor1,1:mstor2) = XSTOR2(1:mstor1,1:mstor2,ICOU)
          XSTORV(1:nstorv) = XSTORV2(1:nstorv,ICOU)
          endif

          WTRSIG=WTR*(SIGTOT-SIGBGK)

C  COLLISION ESTIMATOR FOR COPV(2), COPV(4), COPV(6)

          IF ((LAST_EVENT%IFLAG == 1) .AND.
     .        (LAST_EVENT%NCELL == IRD)) THEN

!  use collision estimator for first cell ("brick") along track
!  in case of a collision sample 1 (the full weight)
!  in case of no collision sample 0 
            IF ((IFLAG == 4).OR.(IFLAG == 5)) THEN
!  TOTAL
              COPV(2,IRD) = COPV(2,IRD) - WEIGHT*E0
              COPV(12,IRD) = COPV(12,IRD) - WEIGHT*E0
!  FRACTION ACCORDING TO ISP
              COPV(ISP,IRD) = COPV(ISP,IRD) - WEIGHT*E0
!pb              if (nltrc) write(iunout,*) ' start cell '
!pb              if (nltrc) write(iunout,*) ' ird ',ird,LAST_EVENT%NCELL
!pb              if (nltrc) write(iunout,*) ' weight*e0 ',weight*e0
            ELSE
!             add 0 ==> nothing to be done
!pb              if (nltrc) write(iunout,*) ' 0 added, ird ', ird
            END IF

          ELSE

C TRACKLENGTH ESTIMATOR FOR COPV(2), COPV(4), COPV(6)

!  TOTAL, TRACKLENGTH
            COPV(2,IRD) = COPV(2,IRD) - WTRSIG*E0
            COPV(13,IRD) = COPV(13,IRD) - WTRSIG*E0
!  FRACTION ACCORDING TO ISP, TRACKLENGTH
            COPV(ISP,IRD) = COPV(ISP,IRD) - WTRSIG*E0
C
!pb            if (nltrc) write (iunout,*) ' normal case '
!pb            if (nltrc) write (iunout,*) ' ird, wtrsig*e0 ',ird,WTRSIG*E0

          END IF
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+2)=.TRUE.
          LMETSP(NSPTOT+NADVI+NALVI+NCLVI+ISP)=.TRUE.

!  I_x, I_y, I_z , radiative  vector flux, to be related to  kappa* grad(Te)
          IF (LTHICK) THEN
            COPV(7,IRD) = COPV(7,IRD) + WTR*E0*VEL*VELX
            COPV(8,IRD) = COPV(8,IRD) + WTR*E0*VEL*VELY
            COPV(9,IRD) = COPV(9,IRD) + WTR*E0*VEL*VELZ
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+7)=.TRUE.
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+8)=.TRUE.
            LMETSP(NSPTOT+NADVI+NALVI+NCLVI+9)=.TRUE.
          END IF
C
C
        END DO
C
      ENDIF
C
      RETURN
      END
