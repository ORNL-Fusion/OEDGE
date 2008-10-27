C
C
      SUBROUTINE GETSCL (ISTRA,FA,FM,FI,FPH)
C
C  FIND SCALING FACTORS TO ENFORCE PARTICLE BALANCE
C  SIMPLE VERSION: NOT SPLIT BY SPECIES, ONLY BY TYPE
C  MODIFIED JAN/95: INCLUDE SURFACE TALLIES IN MATRIX, NOT IN
C  INHOMOGENITY
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COUTAU

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: ISTRA
      REAL(DP), INTENT(OUT) :: FA, FM, FI, FPH
      REAL(DP) :: FC(3), P(3,4), B(3)
      REAL(DP) :: DTB1, DTB2, DTA, DETER, DTB3, FNEN
      REAL(DP) :: P11, P12, P13, P21, P22, P23, P31, P32, P33,
     .          B1, B2, B3
      INTEGER :: ICOL, IROW, I, J, I1, I2, J1, J2
      LOGICAL :: LCOLM(3), LROW(3)
C
C
      FC(1)=1.
      FC(2)=1.
      FC(3)=1.
C
C P(..,1)*FC(1)
C P(..,2)*FC(2)
C P(..,3)*FC(3)
C
      P(1,1)=PAATI(0,ISTRA)+POTATI(0,ISTRA)+PRFAAI(0,ISTRA)+
     .       PGENAI(0,ISTRA)
      P(1,2)=PMATI(0,ISTRA)+PRFMAI(0,ISTRA)
      P(1,3)=PIATI(0,ISTRA)+PRFIAI(0,ISTRA)
      P(2,1)=PAMLI(0,ISTRA)+PRFAMI(0,ISTRA)
      P(2,2)=PMMLI(0,ISTRA)+POTMLI(0,ISTRA)+PRFMMI(0,ISTRA)+
     .       PGENMI(0,ISTRA)
      P(2,3)=PIMLI(0,ISTRA)+PRFIMI(0,ISTRA)
      P(3,1)=PAIOI(0,ISTRA)+PRFAII(0,ISTRA)
      P(3,2)=PMIOI(0,ISTRA)+PRFMII(0,ISTRA)
      P(3,3)=PIIOI(0,ISTRA)+POTIOI(0,ISTRA)+PRFIII(0,ISTRA)+
     .       PGENII(0,ISTRA)
C
      B(1)=-(PPATI(0,ISTRA)+WTOTA(0,ISTRA))
      B(2)=-(PPMLI(0,ISTRA)+WTOTM(0,ISTRA))
      B(3)=-(PPIOI(0,ISTRA)+WTOTI(0,ISTRA))
C
      ICOL=0
      IROW=0
      DO 1 I=1,3
        LROW(I)=P(I,1)**2+P(I,2)**2+P(I,3)**2.GT.EPS30
        LCOLM(I)=P(1,I)**2+P(2,I)**2+P(3,I)**2.GT.EPS30
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
         DO 10 I=1,3
           IF (LROW(I)) THEN
             IF (LCOLM(1)) THEN
               FC(1)=(B(1)-P(I,2)-P(I,3))/P(I,1)
             ELSEIF (LCOLM(2)) THEN
               FC(2)=(B(2)-P(I,3))/P(I,2)
             ELSE
C  NOTHING TO BE DONE, ALL FC'S ARE 1.
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
        DO 20 J=1,3
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
          DO 30 I=1,3
            IF (LCOLM(I)) FC(I)=B(J1)/P(J1,I)
30        CONTINUE
C
        ELSE
C  MORE THAN ONE COLUMN IS NON ZERO
          I1=0
          I2=0
          IF (ICOL.EQ.2) THEN
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
            DO 40 I=1,3
              IF (LCOLM(I)) THEN
                IF (I1.EQ.0) THEN
                  I1=I
                ELSE
                  I2=I
                ENDIF
              ENDIF
40          CONTINUE
C
          ELSE
            I1=1
            I2=2
            B(J1)=B(J1)-P(J1,3)
            B(J2)=B(J2)-P(J2,3)
          ENDIF
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
C
        IF (ICOL.EQ.1) THEN
          DO 50 I=1,3
            IF (LCOLM(I)) FC(I)=B(1)/P(1,I)
50        CONTINUE
C
        ELSEIF (ICOL.EQ.2) THEN
          I1=0
          I2=0
C  DETERMINE THE INDICES FOR THE NON ZERO COLUMNS
          DO 60 I=1,3
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
        ELSE
C  THE WHOLE MATRIX IS TO BE USED
          P11=P(1,1)
          P21=P(2,1)
          P31=P(3,1)
          P12=P(1,2)
          P22=P(2,2)
          P32=P(3,2)
          P13=P(1,3)
          P23=P(2,3)
          P33=P(3,3)
          B1=B(1)
          B2=B(2)
          B3=B(3)
          dta=deter(p11,p21,p31,
     .              p12,p22,p32,
     .              p13,p23,p33)
          dtb1=deter(b1,b2,b3,
     .               p12,p22,p32,
     .               p13,p23,p33)
          dtb2=deter(p11,p21,p31,
     .               b1,b2,b3,
     .               p13,p23,p33)
          dtb3=deter(p11,p21,p31,
     .               p12,p22,p32,
     .               b1,b2,b3)
          fc(1)=dtb1/(dta+1.d-30)
          fc(2)=dtb2/(dta+1.d-30)
          fc(3)=dtb3/(dta+1.d-30)
        ENDIF
      ENDIF
C
1000  CONTINUE

!  FOR THE TIME BEING
      FPH=1._DP

      CALL LEER(1)
      WRITE (6,*) 'EIRENE RECOMMENDED RESCALING OF VOLUME AVERAGED '
      WRITE (6,*) 'TALLIES DUE TO STATISTICAL ERRORS IN BALANCE '
      CALL MASR4 ('FATM,FMOL,FION,FPHOT            ',
     .             FC(1),FC(2),FC(3),fph)
      CALL LEER(2)
C
      FA=FC(1)
      FM=FC(2)
      FI=FC(3)

      
C
      RETURN
      END
