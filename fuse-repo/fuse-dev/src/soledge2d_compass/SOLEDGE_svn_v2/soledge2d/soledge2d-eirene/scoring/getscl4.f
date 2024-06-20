C
C
      SUBROUTINE EIRENE_GETSCL4 (ISTRA,FA,FM,FI,FPH)
C
C  FIND SCALING FACTORS TO ENFORCE PARTICLE BALANCE
C  SIMPLE VERSION: NOT SPLIT BY SPECIES, ONLY BY TYPE
C  MODIFIED JAN/95: INCLUDE SURFACE TALLIES IN MATRIX, NOT IN
C  INHOMOGENITY
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_CCONA
      USE EIRMOD_COUTAU
      USE EIRMOD_COMPRT, ONLY: IUNOUT
 
      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), INTENT(OUT) :: FA, FM, FI, FPH
      REAL(DP) :: FC(4), P(4,5), B(4), PP(4,4), FFC(3), BB(3)
      REAL(DP) :: DTB1, DTB2, DTA, EIRENE_DETER, DTB3, FNEN, DTB4, 
     .            EIRENE_DETER4X4
      REAL(DP) :: P11, P12, P13, P21, P22, P23, P31, P32, P33,
     .            B1, B2, B3
      INTEGER :: ICOL, IROW, I, J, I1, I2, J1, J2, IC, JOUT, IOUT
      LOGICAL :: LCOLM(4), LROW(4)
C
C
      FC(1)=1.
      FC(2)=1.
      FC(3)=1.
      FC(4)=1.
C
C P(..,1)*FC(1)
C P(..,2)*FC(2)
C P(..,3)*FC(3)
C P(..,4)*FC(4)
C
      P(1,1)=PAATI(0,ISTRA)+POTATI(0,ISTRA)+PRFAAI(0,ISTRA)+
     .       PGENAI(0,ISTRA)
      P(1,2)=PMATI(0,ISTRA)+PRFMAI(0,ISTRA)
      P(1,3)=PIATI(0,ISTRA)+PRFIAI(0,ISTRA)
      P(1,4)=0._DP
      P(2,1)=PAMLI(0,ISTRA)+PRFAMI(0,ISTRA)
      P(2,2)=PMMLI(0,ISTRA)+POTMLI(0,ISTRA)+PRFMMI(0,ISTRA)+
     .       PGENMI(0,ISTRA)
      P(2,3)=PIMLI(0,ISTRA)+PRFIMI(0,ISTRA)
      P(2,4)=0._DP
      P(3,1)=PAIOI(0,ISTRA)+PRFAII(0,ISTRA)
      P(3,2)=PMIOI(0,ISTRA)+PRFMII(0,ISTRA)
      P(3,3)=PIIOI(0,ISTRA)+POTIOI(0,ISTRA)+PRFIII(0,ISTRA)+
     .       PGENII(0,ISTRA)
      P(3,4)=0._DP
      P(4,1)=0._DP
      P(4,2)=0._DP
      P(4,3)=0._DP
      P(4,4)=1._DP
C
      B(1)=-(PPATI(0,ISTRA)+WTOTA(0,ISTRA))
      B(2)=-(PPMLI(0,ISTRA)+WTOTM(0,ISTRA))
      B(3)=-(PPIOI(0,ISTRA)+WTOTI(0,ISTRA))
      B(4)=1._DP
C
      ICOL=0
      IROW=0
      DO 1 I=1,4
        LROW(I)=P(I,1)**2+P(I,2)**2+P(I,3)**2+P(I,4)**2.GT.EPS30
        LCOLM(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2+P(4,I)**2.GT.EPS30
        IF (LROW(I)) IROW=IROW+1
        IF (LCOLM(I)) ICOL=ICOL+1
1     CONTINUE
C
      IF (IROW.EQ.0) THEN
C  NO ROW IS NON ZERO, I.E. NO PARTICLES FOLLOWED
        GOTO 1000
C
      ELSEIF (IROW.EQ.1) THEN
C  ONLY ONE ROW (NO. I) IS NON ZERO, I.E., ONLY ATOMS, ONLY MOLECULES
C                                    OR  ONLY TEST IONS ARE FOLLOWED
         DO 10 I=1,4
           IF (LROW(I)) THEN
             IF (LCOLM(1)) THEN
               FC(1)=(B(I)-P(I,2)-P(I,3)-P(I,4))/P(I,1)
             ELSEIF (LCOLM(2)) THEN
               FC(2)=(B(I)-P(I,3)-P(I,4))/P(I,2)
             ELSEIF (LCOLM(3)) THEN
               FC(3)=(B(I)-P(I,4))/P(I,3)
             ELSE
               FC(4)=B(I)/P(I,4)
             ENDIF
           ENDIF
10       CONTINUE
C
C
      ELSEIF (IROW.EQ.2) THEN
C  TWO ROWS ARE NON ZERO
        J1=0
        J2=0
C  DETERMINE THE INDICES FOR THE NON ZERO ROWS
        DO 20 J=1,4
          IF (LROW(J)) THEN
            IF (J1.EQ.0) THEN
              J1=J
            ELSE
              J2=J
            ENDIF
          ENDIF
20      CONTINUE
C
        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NON ZERO
          DO 30 I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
30        CONTINUE
C
        ELSE
C  MORE THAN ONE COLUMN IS NON ZERO
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE FIRST TWO NON ZERO COLUMNS
          DO 40 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
40        CONTINUE
C
! nur spalten rechts von spalte i2 koennen noch werte enthalten
          DO I=I2+1,4
            B(J1)=B(J1)-P(J1,I)
            B(J2)=B(J2)-P(J2,I)
          END DO
C
          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF
        ENDIF
C
C
      ELSEIF (IROW.EQ.3) THEN
C  THREE ROWS ARE NON ZERO
C
C  DETERMINE THE ROW WHICH IS COMPLETELY ZERO
        JOUT = 0
        DO J=1,4
          IF (.NOT.LROW(J)) JOUT=J
        END DO
 
        IF (ICOL.EQ.1) THEN
C  ONLY ONE COLUMN IS NON ZERO
          J1=1
          IF (JOUT.EQ.1) J1=2
          DO  I=1,4
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
          END DO
 
        ELSEIF (ICOL.EQ.2) THEN
C  TWO COLUMNS ARE NON ZERO
          J1=1
          IF (JOUT.EQ.J1) J1=J1+1
          J2=J1+1
          IF (JOUT.EQ.J2) J2=J2+1
 
C  DETERMINE THE INDICES FOR THE FIRST TWO NON ZERO COLUMNS
          I1=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
                EXIT
              ENDIF
            ENDIF
          END DO
 
          FNEN=P(J1,I1)*P(J2,I2)-P(J2,I1)*P(J1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(J1)*P(J2,I2)-B(J2)*P(J1,I2))/FNEN
            IF (ABS(P(J1,I2)).GT.EPS12) THEN
              FC(I2)=(B(J1)-P(J1,I1)*FC(I1))/P(J1,I2)
            ELSEIF (ABS(P(J2,I2)).GT.EPS12) THEN
              FC(I2)=(B(J2)-P(J2,I1)*FC(I1))/P(J2,I2)
            ENDIF
          ENDIF
 
        ELSE
 
C  AT LEAST THREE COLUMNS ARE NON ZERO
          IF (ICOL.EQ.4) THEN
            DO J=1,4
              B(J)=B(J)-P(J,4)
            END DO
            IOUT=4
          ELSE
            IOUT=0
            DO I=1,4
              IF (.NOT.LCOLM(I)) IOUT=I
            END DO
          END IF
 
          J1=0
          DO J=1,4
            IF (J.EQ.JOUT) CYCLE
            J1=J1+1
            I1=0
            BB(J1)=B(J)
            DO I=1,4
              IF (I.EQ.IOUT) CYCLE
              I1=I1+1
              PP(J1,I1)=P(J,I)
            END DO
          END DO
 
          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=BB(1)
          B2=BB(2)
          B3=BB(3)
          dta=EIRENE_deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=EIRENE_deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=EIRENE_deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=EIRENE_deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          Ffc(1)=dtb1/(dta+1.d-30)
          Ffc(2)=dtb2/(dta+1.d-30)
          Ffc(3)=dtb3/(dta+1.d-30)
 
          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
              IF (IC == 3) EXIT
            END IF
          END DO
        END IF
C
C
      ELSEIF (IROW.EQ.4) THEN
C  ALL ROWS ARE NON ZERO
C
        IF (ICOL.EQ.1) THEN
          DO 50 I=1,4
            IF (LCOLM(I)) FC(I)=B(I)/P(1,I)
50        CONTINUE
C
        ELSEIF (ICOL.EQ.2) THEN
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
          DO 60 I=1,4
            IF (LCOLM(I)) THEN
              IF (I1.EQ.0) THEN
                I1=I
              ELSE
                I2=I
              ENDIF
            ENDIF
60        CONTINUE
C
          FNEN=P(1,I1)*P(2,I2)-P(2,I1)*P(1,I2)
          IF (ABS(FNEN).GT.EPS12) THEN
            FC(I1)=(B(1)*P(2,I2)-B(2)*P(1,I2))/FNEN
            IF (ABS(P(1,I2)).GT.EPS12) THEN
              FC(I2)=(B(1)-P(1,I1)*FC(I1))/P(1,I2)
            ELSEIF (ABS(P(2,I2)).GT.EPS12) THEN
              FC(I2)=(B(2)-P(2,I1)*FC(I1))/P(2,I2)
            ENDIF
          ENDIF
C
        ELSEIF (ICOL.EQ.3) THEN
 
          IC=0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              PP(1:3,IC) = P(1:3,I)
            END IF
          END DO
 
          P11=PP(1,1)
          P21=PP(2,1)
          P31=PP(3,1)
          P12=PP(1,2)
          P22=PP(2,2)
          P32=PP(3,2)
          P13=PP(1,3)
          P23=PP(2,3)
          P33=PP(3,3)
          B1=B(1)
          B2=B(2)
          B3=B(3)
          dta=EIRENE_deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=EIRENE_deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=EIRENE_deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=EIRENE_deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          Ffc(1)=dtb1/(dta+1.d-30)
          Ffc(2)=dtb2/(dta+1.d-30)
          Ffc(3)=dtb3/(dta+1.d-30)
 
          IC = 0
          DO I=1,4
            IF (LCOLM(I)) THEN
              IC = IC + 1
              FC(I) = FFC(IC)
            END IF
          END DO
 
        ELSE
C  THE WHOLE MATRIX IS TO BE USED
          pp(1:4,1:4) = p(1:4,1:4)
          dta=EIRENE_deter4x4(pp)
 
          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,1) = b(1:4)
          dtb1=EIRENE_deter4x4(pp)
 
          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,2) = b(1:4)
          dtb2=EIRENE_deter4x4(pp)
 
          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,3) = b(1:4)
          dtb3=EIRENE_deter4x4(pp)
 
          pp(1:4,1:4) = p(1:4,1:4)
          pp(1:4,4) = b(1:4)
          dtb4=EIRENE_deter4x4(pp)
 
          fc(1)=dtb1/(dta+1.d-30)
          fc(2)=dtb2/(dta+1.d-30)
          fc(3)=dtb3/(dta+1.d-30)
          fc(4)=dtb4/(dta+1.d-30)
        ENDIF
      ENDIF
C
1000  CONTINUE
 
!  FOR THE TIME BEING
 
      CALL EIRENE_LEER(1)
      WRITE (iunout,*)
     .  'EIRENE RECOMMENDED RESCALING OF VOLUME AVERAGED '
      WRITE (iunout,*) 'TALLIES DUE TO STATISTICAL ERRORS IN BALANCE '
      CALL EIRENE_MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .             FC(1),FC(2),FC(3),FC(4))
      CALL EIRENE_LEER(2)
C
      FA=FC(1)
      FM=FC(2)
      FI=FC(3)
      FPH=FC(4)
 
 
C
      RETURN
      END
