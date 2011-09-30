C
C
      SUBROUTINE CDEFN(AL,PL,K,COU,NTE,CF,LEXP,LTEST,LSUM)
C
C  EIRENE DEFAULT ATOMIC DATA FOR INTERACTION WITH HYDROGEN
C  SAME AS CDEF, BUT FOR 2 PARAMETER FITTING EXPRESSIONS
C
C  INPUT:
C    AL(J),J=1,NTE, PL(I),I=1,NTE
C  OUTPUT:
C    MAXWELLIAN RATES, AL=LN(KT), KT IN (EV)
C                      PL=LN(NE), NE IN (1/CM**3)
C  K
C       K>0: A&M DATA FROM FILES HYDHEL, METHANE OR AMJUEL
C           K:  NUMBER OF REACTION IN EIRENE "CREAC"-ARRAY
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL(*), PL(*)
      REAL(DP), INTENT(OUT) :: COU(0:9,*), CF(9,0:9)
      INTEGER, INTENT(IN) :: K, NTE
      LOGICAL, INTENT(IN) :: LEXP, LTEST, LSUM
      REAL(DP) :: DUMMP(9)
      REAL(DP) :: CCXM1, CCXM2, EXPO1, EXPO2, FPAR1, FPAR2, FPAR3, S01,
     .          S02, DS12, CTEST, EXTRAP
      INTEGER :: I, JJ, KK, IFEX, J, II, ICELL
C
C  K>0:
C  DATA FROM ARRAY CREAC(9,0:9,K)
C
      IF (K.GT.0) THEN
C  DATA FROM A&M DATA FILES
        IF (.NOT.REACDAT(K)%LRTC) THEN
          WRITE (IUNOUT,*) ' NO DATA AVAILABLE FOR RATE COEFFICIENT',K
          CALL EXIT_OWN(1)
        END IF
        IF (REACDAT(K)%RTC%IFIT == 1) THEN
          WRITE (IUNOUT,*) ' ONLY 1D FIT AVAILABLE FOR RATE COEFF.',K
          CALL EXIT_OWN(1)
        END IF
        DO 11 J=1,9
          DO 12 II=1,9
            CF(II,J)=REACDAT(K)%RTC%POLY%DBLPOL(II,J)
12        CONTINUE
11      CONTINUE
        IF (LTEST) THEN
          DO 13 J=1,9
            CTEST=0.
            DO 14 II=1,9
              CTEST=CTEST+ABS(CF(II,J))
14          CONTINUE
            IF (CTEST.LE.EPS60) GOTO 990
13        CONTINUE
        ENDIF
      ELSE
        GOTO 990
      ENDIF
C
      IF (LSUM) THEN
        DO 20 ICELL=1,NTE
          COU(1,ICELL)=0.
20      CONTINUE
C
        DO 25 ICELL=1,NTE
          IF (AL(ICELL).LT.REACDAT(K)%RTC%POLY%RCMN) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(<S*V>)
            S01=REACDAT(K)%RTC%POLY%RCMN
            S02=LOG(2._DP)+REACDAT(K)%RTC%POLY%RCMN
            DS12=S02-S01
            EXPO1=0.
            EXPO2=0.
            DO 1 J=1,9
              JJ=J-1
              DO 1 I=1,9
                II=I-1
                EXPO1=EXPO1+S01**II*PL(ICELL)**JJ*CF(I,J)
                EXPO2=EXPO2+S02**II*PL(ICELL)**JJ*CF(I,J)
1           CONTINUE
            CCXM1=EXPO1
            CCXM2=EXPO2
            FPAR1=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
            FPAR2=      (CCXM2-CCXM1)/DS12
            FPAR3=0.D0
C
            IFEX=5
            COU(1,ICELL)=EXTRAP(AL(ICELL),IFEX,FPAR1,FPAR2,FPAR3)
            if (.not.lexp) cou(1,icell)=log(cou(1,icell))
C
          ELSE
C
            DO 22 JJ=9,1,-1
              DUMMP(JJ)=CF(9,JJ)
              DO 23 KK=8,1,-1
                DUMMP(JJ)=DUMMP(JJ)*AL(ICELL)+CF(kk,jj)
23            CONTINUE
22          CONTINUE
            cou(1,icell)=dummp(9)
            DO 24 JJ=8,1,-1
              cou(1,icell)=cou(1,icell)*PL(icell)+DUMMP(JJ)
24          CONTINUE
C
C           DO 22 J=1,9
C             JJ=J-1
C             DO 22 I=1,9
C               II=I-1
C               COU(1,ICELL)=COU(1,ICELL)+AL(ICELL)**II*
C    .                                    PL(ICELL)**JJ*CF(I,J)
C22          CONTINUE
C
            if (lexp) COU(1,ICELL)=EXP(MAX(-100._DP,COU(1,ICELL)))
          ENDIF
25      CONTINUE
      ENDIF
C
C
      RETURN
C
C
990   CONTINUE
      WRITE (iunout,*) 
     .  'ERROR IN SUBROUTINE CDEFN: ZERO FIT COEFFICIENTS'
      WRITE (iunout,*) 'J,K = ',J,K,'  EXIT CALLED!'
      CALL EXIT_OWN(1)
      END
