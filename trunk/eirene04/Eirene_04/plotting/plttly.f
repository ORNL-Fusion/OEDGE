C
C
      SUBROUTINE PLTTLY (X,Y,VBAR,YMN,YMX,IR1,IR2,IRS,NKURV,TXTTAL,
     .                   TXTSPC,TXTUNT,TXTRUN,TXHEAD,
     .                   LBAR,XMI,XMA,YMNLG,YMXLG,LPLOT,LHIST,IERR,
     .                   N1BAR,N1DIM)

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE CPLMSK

      IMPLICIT NONE
 
      INTEGER, INTENT(IN) :: N1BAR, N1DIM
      REAL(DP), INTENT(IN) :: X(*), VBAR(N1BAR,*),
     .                      YMN(*), YMX(*), YMNLG(*), YMXLG(*)
      REAL(DP), INTENT(INOUT) :: Y(N1DIM,*)
      INTEGER, INTENT(IN) :: IR1(*), IR2(*), IRS(*), NKURV
      LOGICAL, INTENT(IN) :: LBAR(*), LPLOT(*)
      CHARACTER(LEN=*), INTENT(IN) :: TXTTAL(*),TXTSPC(*),TXTUNT(*),
     .                                TXTRUN, TXHEAD
      REAL(DP) :: YA, YMINY, YMY, AA, FM, ST1, ST2, FP, DMINY, DMAXY,
     .          XMI, XMA
      INTEGER :: I1, I2, IS, ICURV, IT, ISY, NP, J, NPS, IKURV, 
     .           I, IERR
      LOGICAL :: LHIST
      CHARACTER(10) :: CHR
C
      IKURV=0
      DO 1 I=1,NKURV
        IF (LPLOT(I)) IKURV=IKURV+1
1     CONTINUE
      IF (IKURV.EQ.0) RETURN
C
      MINX=XMI
      MAXX=XMA
      MINY=1.D30
      MAXY=-1.D30
      DO 3 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 3
        MINY=MIN(MINY,REAL(YMN(I),SP))
        MAXY=MAX(MAXY,REAL(YMX(I),SP))
3     CONTINUE
C
      IF (LOGY) THEN
        DMAXY=MAXY
        DMINY=MINY
        DMAXY=MAX(1.E-36_DP,DMAXY)
        MAXY=DMAXY
        DMINY=MAX(DMINY,DMAXY/1.E12_DP)
        MINY=DMINY
      ENDIF
C
      IF (MINY.GE.MAXY) THEN
        IERR=1
        RETURN
      ENDIF
C
      CALL GRNXTB(1)
      CALL GRSCLC (0.,0.,39.5,28.7)
      CALL GRSCLV (0.,0.,39.5,28.7)
      IT=LEN(TXTRUN)
      CALL GRTXT (1.,27.5,IT,TXTRUN)
      IT=LEN(TXHEAD)
      CALL GRTXT (1.,26.75,IT,TXHEAD)
      CALL GRTXT (1.,26.,15,'MIN. ABSZISSA =')
      WRITE (CHR,'(1P,E10.3)') XMI
      CALL GRTXTC (10,CHR)
      CALL GRTXT (1.,25.5,15,'MAX. ABSZISSA =')
      WRITE (CHR,'(1P,E10.3)') XMA
      CALL GRTXTC (10,CHR)
      YA=24.75
      DO 2 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 2
        CALL GRNWPN(I)
        ISY=I+1
        CALL GRJMPS (0.5,REAL(YA,KIND(1.E0)),ISY)
        IT=LEN(TXTTAL(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTTAL(I))
        YA=YA-0.5
        IT=LEN(TXTSPC(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTSPC(I))
        YA=YA-0.5
        IT=LEN(TXTUNT(I))
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),IT,TXTUNT(I))
        YA=YA-0.5
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),11,'MAX. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMXLG(I)
        CALL GRTXTC (10,CHR)
        YA=YA-0.5
        CALL GRTXT (1.5,REAL(YA,KIND(1.E0)),11,'MIN. VALUE =')
        WRITE (CHR,'(1P,E10.3)') YMNLG(I)
        CALL GRTXTC (10,CHR)
        YA=YA-1.0
2     CONTINUE
      CALL GRNWPN(1)
C
C
C  PLOT X AND Y AXIS
C
      CALL PLTMSK(IERR)
      IF (IERR.GT.0) RETURN
C
C  PREPARE ARRAYS FOR PLOTTING ON LOG. SCALE
C
      IF (LOGY) THEN
        YMINY=MAX(10.D0**MINLY,10.D0**(MAXLY-12))
        DO 20 ICURV=1,NKURV
          IF (.NOT.LPLOT(ICURV)) GOTO 20
          DO 21 I=IR1(ICURV),IR2(ICURV)-1,IRS(ICURV)
21          Y(I,ICURV)=LOG10(MAX(YMINY,Y(I,ICURV)))
20      CONTINUE
      ENDIF
C
C  PLOT !
C
      DO 50 I=1,NKURV
        IF (.NOT.LPLOT(I)) GOTO 50
        CALL GRNWPN(I)
        I1=IR1(I)
        IS=IRS(I)
        I2=IR2(I)-IS
C
        IF (LHIST) THEN
C  PLOT LINES
          IF (LOGY) THEN
            YMY=MINLY
          ELSE
            YMY=MINY
          ENDIF
          CALL GRJMP(REAL(X(I1),KIND(1.E0)),REAL(YMY,KIND(1.E0)))
          CALL GRDRW(REAL(X(I1),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          CALL GRDRW(REAL(X(I1+IS),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          DO J=I1+IS,I2,IS
            CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)))
            CALL GRDRW (REAL(X(J+IS),KIND(1.E0)),
     .                  REAL(Y(J,I),KIND(1.E0)))
          END DO
          CALL GRDRW(REAL(X(I2+IS),KIND(1.E0)),REAL(YMY,KIND(1.E0)))
C  PLOT SYMBOLS
          ISY=I+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO J=I1,I2,IS*NPS
            CALL GRJMPS (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                   REAL(Y(J,I),KIND(1.E0)),ISY)
          END DO
C
        ELSEIF (.NOT.LHIST) THEN
C  PLOT LINES
          CALL GRJMP(REAL(X(I1),KIND(1.E0)),REAL(Y(I1,I),KIND(1.E0)))
          DO 33 J=I1+IS,I2,IS
33          CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)))
C  PLOT SYMBOLS
          ISY=I+1
          NP=(I2-I1+1)/IS
          NPS=MAX0(NP/7,1)
          DO J=I1,I2,IS*NPS
            CALL GRJMPS (REAL(X(J),KIND(1.E0)),REAL(Y(J,I),KIND(1.E0)),
     .                   ISY)
          END DO
        ENDIF
C
C  PLOT ERROR BARS
          IF (LBAR(I)) THEN
            DO 40 J=I1,I2,IS
              FP=1.+VBAR(J,I)*0.01
              FM=1.-VBAR(J,I)*0.01
              IF (LOGY) THEN
                AA=10.**Y(J,I)
                ST1=LOG10(MAX(1.E-30_DP,AA*FM))
                ST2=LOG10(MAX(1.E-30_DP,AA*FP))
                ST1=MAX(REAL(MINLY,DP),MIN(REAL(MAXLY,DP),ST1))
                ST2=MAX(REAL(MINLY,DP),MIN(REAL(MAXLY,DP),ST2))
              ELSE
                ST1=Y(J,I)*FM
                ST2=Y(J,I)*FP
                ST1=MAX(REAL(MINY,DP),MIN(REAL(MAXY,DP),ST1))
                ST2=MAX(REAL(MINY,DP),MIN(REAL(MAXY,DP),ST2))
              ENDIF
            IF (LHIST) THEN
              CALL GRJMP (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                    REAL(ST1,KIND(1.E0)))
              CALL GRDRW (REAL(0.5*(X(J)+X(J+IS)),KIND(1.E0)),
     .                    REAL(ST2,KIND(1.E0)))
            ELSE
              CALL GRJMP (REAL(X(J),KIND(1.E0)),REAL(ST1,KIND(1.E0)))
              CALL GRDRW (REAL(X(J),KIND(1.E0)),REAL(ST2,KIND(1.E0)))
            ENDIF
40        CONTINUE
        ENDIF
50    CONTINUE
      CALL GRNWPN(1)
      CALL GRCHRC (0.3,0.,16)
C
      RETURN
      END
