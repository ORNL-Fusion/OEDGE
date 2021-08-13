C
C
      FUNCTION CROSS(AL,K,IR,TEXT)
C
C  CROSS SECTION
C    AL=LN(ELAB), ELAB IN (EV)
C    RETURN CROSS SECTION IN CM**2
C
C  K>0 :  DATA FROM ARRAY CREAC(9,0:9,K)
C
C  K<0 :  DEFAULT MODEL
C
C  K=-1:  H + H+ --> H+ + H   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 3.1.8
C
      USE PRECISION
      USE PARMMOD
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL
      INTEGER, INTENT(IN) :: K, IR
      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      REAL(DP) :: CF1(9),B(8)
      REAL(DP) :: RMIN, FP1, FP2, S01, S02, DS12, EXPO1, EXPO2, CROSS,
     .            CCXM1, CCXM2, EXPO, EXTRAP, E, XI
      INTEGER :: IF8, II, I

      DATA CF1
     ./-3.274124E+01,-8.916457E-02,-3.016991E-02,
     .  9.205482E-03, 2.400267E-03,-1.927122E-03,
     .  3.654750E-04,-2.788866E-05, 7.422296E-07/
      DATA RMIN
     ./-2.3025851E+00/
      DATA FP1,FP2
     ./-3.2945896E+01,-1.713112E-01/
C
      IF (K.EQ.-1) THEN
C  ELAB BELOW MIMINUM ENERGY FOR FIT:
        IF (AL.LT.RMIN) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFMN=5
          CROSS=EXTRAP(AL,5,FP1,FP2,0._DP)
C  ELAB ABOVE MAXIMUM ENERGY FOR FIT: NOT IN USE FOR K=-1
C       ELSEIF (ELAB.GT.RMAX) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFMX=5
C         CROSS=EXTRAP(AL,5,FP4,FP5,0.D0)
C
        ELSE
          EXPO=CF1(9)
          DO 10 II=1,8
            IF8=9-II
            EXPO=EXPO*AL+CF1(IF8)
10        CONTINUE
          CROSS=EXP(MAX(-100._DP,EXPO))
        ENDIF
C
      ELSEIF (K.GT.0) THEN

        IF (IFTFLG(K,1) == 0) THEN

C  ELAB BELOW MINIMUM ENERGY FOR FIT:
          IF (AL.LT.RCMN(K,1)) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
            IF (IFEXMN(K,1).LT.0) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(SIGMA)
              S01=RCMN(K,1)
              S02=LOG(2._DP)+RCMN(K,1)
              DS12=S02-S01
              EXPO1=CREAC(9,0,K)
              EXPO2=CREAC(9,0,K)
              DO 1 II=1,8
                IF8=9-II
                EXPO1=EXPO1*S01+CREAC(IF8,0,K)
                EXPO2=EXPO2*S02+CREAC(IF8,0,K)
 1            CONTINUE
              CCXM1=EXPO1
              CCXM2=EXPO2
              FPARM(K,1,1)=CCXM1+(CCXM2-CCXM1)/DS12*(-S01)
              FPARM(K,2,1)=      (CCXM2-CCXM1)/DS12
              FPARM(K,3,1)=0.D0
C     
              IFEXMN(K,1)=5
            ENDIF
            CROSS=EXTRAP(AL,IFEXMN(K,1),FPARM(K,1,1),FPARM(K,2,1),
     .                                               FPARM(K,3,1))
C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
          ELSEIF (AL.GT.RCMX(K,1)) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
            CROSS=EXTRAP(AL,IFEXMX(K,1),FPARM(K,4,1),FPARM(K,5,1),
     .                                               FPARM(K,6,1))
          ELSE
            EXPO=CREAC(9,0,K)
            DO 100 II=1,8
              IF8=9-II
              EXPO=EXPO*AL+CREAC(IF8,0,K)
 100        CONTINUE
            CROSS=EXP(MAX(-100._DP,EXPO))
          ENDIF

        ELSE IF (IFTFLG(K,1) == 3) THEN
C  default extrapolation ifexmn=-1 not yet available
C  ELAB BELOW MINIMUM ENERGY FOR FIT:
          IF (AL.LT.RCMN(K,1)) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
            CROSS=EXTRAP(AL,IFEXMN(K,1),FPARM(K,1,1),FPARM(K,2,1),
     .                                               FPARM(K,3,1))
C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
          ELSEIF (AL.GT.RCMX(K,1)) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
            CROSS=EXTRAP(AL,IFEXMX(K,1),FPARM(K,4,1),FPARM(K,5,1),
     .                                               FPARM(K,6,1))
          ELSE
            E = EXP(AL)
            XI = CREAC(1,0,K)
            B(1:8) = CREAC(2:9,0,K)
            CROSS = B(1)*LOG(E/XI)
            DO I=1,7
              CROSS = CROSS + B(I+1)*(1.D0-XI/E)**I 
            END DO
            CROSS = CROSS * 1.D-13/(XI*E)
          ENDIF
        ELSE
          WRITE (6,*) ' WRONG FITTING FLAG IN CROSS '
          WRITE (6,*) ' K = ',K,' IFTFLG = ',IFTFLG(K,1)
          WRITE (6,*) 'REACTION NO. ',IR
          CALL EXIT_OWN(1)
        END IF
      ELSE
        WRITE (6,*) 'ERROR IN CROSS: K= ',K
        WRITE (6,*) 'CALLED FROM ',TEXT
        WRITE (6,*) 'REACTION NO. ',IR
        CALL EXIT_OWN(1)
      ENDIF
      RETURN
      END
