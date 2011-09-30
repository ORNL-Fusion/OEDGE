C 0406: default resonant cx for He in He+/He++ plasma added:
C       Janev (HYDHEL) ,1987, reactions 5.3.1 and 6.3.1
C
      FUNCTION CROSS(AL,K,IR,TEXT)
C
C  CROSS SECTION
C    AL=LN(ELAB), ELAB IN (EV)
C    RETURN CROSS SECTION IN CM**2
C
C  K>0 :  DATA FROM ARRAY REACDAT, I.E. FROM EXTERNAL DATABASE
C
C  K<0 :  DEFAULT MODEL DEFINED IN SETUP_DEFAULT_REACTIONS
C
C  K=-1:  H + H+ --> H+ + H   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 3.1.8
C
C  K=-2:  He + He+ --> He+ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 5.3.1
C
C  K=-3:  He + He++ --> H++ + He   CROSS SECTION, JANEV, 3.1.8
C         LINEAR EXTRAPOLATION AT LOW ENERGY END FOR LN(SIGMA)
C         IDENTICAL TO hydhel.tex, H.1, 6.3.1
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMPRT, ONLY: IUNOUT

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL
      INTEGER, INTENT(IN) :: K, IR
      CHARACTER(LEN=*), INTENT(IN) :: TEXT
      REAL(DP) :: B(8), FP(6)
      REAL(DP) :: S01, S02, DS12, EXPO1, EXPO2, CROSS,
     .            CCXM1, CCXM2, EXPO, EXTRAP, E, XI, SNGL_POLY
      INTEGER :: IF8, II, I
      type(poly_data), pointer :: rp

C
      IF ((K >= -10) .AND. (K <= NREAC)) THEN
        IF (IFTFLG(K,1) == 0) THEN

          RP => REACDAT(K)%CRS%POLY
          EXPO = SNGL_POLY(RP%DBLPOL,AL,RP%RCMN,RP%RCMX,RP%FPARM,
     .                     RP%IFEXMN,RP%IFEXMX)
          CROSS = EXP(MAX(-100._DP,EXPO))

          CROSS = CROSS*FACREA(K,1)

        ELSE IF (IFTFLG(K,1) == 3) THEN
C  default extrapolation ifexmn=-1 not yet available
C  ELAB BELOW MINIMUM ENERGY FOR FIT:
          IF (AL.LT.REACDAT(K)%CRS%POLY%RCMN) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMN(K)
            FP = REACDAT(K)%CRS%POLY%FPARM
            CROSS=EXTRAP(AL,REACDAT(K)%CRS%POLY%IFEXMN,
     .                   FP(1),FP(2),FP(3))
            CROSS = CROSS*FACREA(K,1)
C  ELAB ABOVE MAXIMUM ENERGY FOR FIT:
          ELSEIF (AL.GT.REACDAT(K)%CRS%POLY%RCMX) THEN
C  USE ASYMPTOTIC EXPRESSION NO. IFEXMX(K,1)
            FP = REACDAT(K)%CRS%POLY%FPARM
            CROSS=EXTRAP(AL,REACDAT(K)%CRS%POLY%IFEXMX,
     .                   FP(4),FP(5),FP(6))
            CROSS = CROSS*FACREA(K,1)
          ELSE
            E = EXP(AL)
            XI = REACDAT(K)%CRS%POLY%DBLPOL(1,1)
            B(1:8) = REACDAT(K)%CRS%POLY%DBLPOL(2:9,1)
            CROSS = B(1)*LOG(E/XI)
            DO I=1,7
              CROSS = CROSS + B(I+1)*(1.D0-XI/E)**I
            END DO
            CROSS = CROSS * 1.D-13/(XI*E)
            CROSS = CROSS*FACREA(K,1)
          ENDIF
        ELSE
          WRITE (iunout,*) ' WRONG FITTING FLAG IN CROSS '
          WRITE (iunout,*) ' K = ',K,' IFTFLG = ',IFTFLG(K,1)
          WRITE (iunout,*) 'REACTION NO. ',IR
          CALL EXIT_OWN(1)
        END IF
      ELSE
        WRITE (iunout,*) 'ERROR IN CROSS: K= ',K
        WRITE (iunout,*) 'CALLED FROM ',TEXT
        WRITE (iunout,*) 'REACTION NO. ',IR
        CALL EXIT_OWN(1)
      ENDIF
      RETURN
      END
