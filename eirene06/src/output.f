C EIRENE06 COMPILATION
C ===== SOURCE: algtal.f
C
C
      SUBROUTINE ALGTAL

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CGEOM
      USE COMPRT
      USE CTEXT
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE
C
      REAL(DP), ALLOCATABLE :: VEC1(:), VEC2(:), RESULT(:,:)
      REAL(DP), ALLOCATABLE :: OP(:), WEI(:), SUMWEI(:)
      REAL(DP) :: CONST(20)
      INTEGER :: IIND(20), IZIF(4,20)
      INTEGER :: I, ITL, IALV, NOP, IOP, K, ILIMPS,
     .           IINDEX, II, J, IALS, IN
      LOGICAL :: LFREE1, LFREE2
      LOGICAL, ALLOCATABLE :: LLIMPS(:)
      LOGICAL :: LLMPS
      CHARACTER(72) :: HCHR
      CHARACTER(1) :: OPER(20)
C
C
C     CALCULATE ALGEBRAIC VOLUME TALLIES
C
C
      IF (NALVI+NALSI <= 0) RETURN

      IF (.NOT.LALGV.AND.NALVI.GT.0) THEN
        WRITE (iunout,*) ' ALGV IS SWITCHED OFF '
        WRITE (iunout,*) 
     .    ' NO ALGEBRAIC VOLUME TALLIES CAN BE CALCULATED '
        GOTO 300
      END IF

      IF (NALVI > 0) THEN
        ALLOCATE (VEC1(MAX(NSBOX_TAL,NLIMPS)))
        ALLOCATE (VEC2(MAX(NSBOX_TAL,NLIMPS)))
        ALLOCATE (RESULT(2,MAX(NSBOX_TAL,NLIMPS)))
      END IF

      DO 200 IALV=1,NALVI
C
        HCHR=CHRTAL(IALV)
C
        CALL ALGEBR (HCHR,OPER,IZIF,CONST,NOP)
C
        DO 1 IOP=1,NOP
          IIND(IOP)=0
C         WRITE (iunout,*) IOP,OPER(IOP),(IZIF(J,IOP),J=1,4)
1       CONTINUE
        LFREE1=.TRUE.
        LFREE2=.TRUE.

        IF (ANY(IZIF(2,1:NOP)<0).OR.ANY(IZIF(4,1:NOP)<0)) THEN
          ALLOCATE(OP(MAX(NSBOX,NLIMPS)))
          ALLOCATE(WEI(MAX(NSBOX,NLIMPS)))
          ALLOCATE(SUMWEI(MAX(NSBOX_TAL,NLIMPS)))
        END IF
C
        DO 100 IOP=1,NOP
C
C  1. OPERAND
C
C  TALLY HOLEN
          IF (IZIF(2,IOP)*IZIF(4,IOP) < 0.D0) GOTO 94
          IF (IZIF(2,IOP).GT.0) THEN
            IF (.NOT.LIVTALV(IZIF(2,IOP))) GOTO 95
            IF (IZIF(2,IOP).GT.NTALV) GOTO 90
            IF (IZIF(1,IOP).GT.NFSTVI(IZIF(2,IOP))) GOTO 91
            DO 10 I=1,NSBOX_TAL
              VEC1(I)=ESTIMV(NADDV(IZIF(2,IOP))+IZIF(1,IOP),I)
10          CONTINUE
C
          ELSEIF (IZIF(2,IOP).LT.0) THEN
            ITL=IABS(IZIF(2,IOP))
            IF (ITL.GT.NTALI) GOTO 90
            IF (IZIF(1,IOP).GT.NFRSTP(ITL)) GOTO 91
            K=IZIF(1,IOP)
            SELECT CASE (ITL)
            CASE (1)
              OP(1:NSBOX)  = TEIN(1:NSBOX)
              WEI(1:NSBOX) = DEIN(1:NSBOX)*VOL(1:NSBOX)
            CASE (2)
              OP(1:NSBOX)  = TIIN(MPLSTI(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (3)
              OP(1:NSBOX) = DEIN(1:NSBOX)
              WEI(1:NSBOX) = VOL(1:NSBOX)
            CASE (4)
              OP(1:NSBOX) = DIIN(K,1:NSBOX)
              WEI(1:NSBOX) = VOL(1:NSBOX)
            CASE (5)
              OP(1:NSBOX) = VXIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (6)
              OP(1:NSBOX) = VYIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (7)
              OP(1:NSBOX) = VZIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (8)
              OP(1:NSBOX) = BXIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (9)
              OP(1:NSBOX) = BYIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (10)
              OP(1:NSBOX) = BZIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (11)
              OP(1:NSBOX) = BFIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (12)
              OP(1:NSBOX) = ADIN(K,1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (13)
              OP(1:NSBOX) = EDRIFT(K,1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (14)
              OP(1:NSBOX) = VOL(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (15)
              OP(1:NSBOX) = WGHT(K,1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (16)
              OP(1:NSBOX) = BXPERP(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (17)
              OP(1:NSBOX) = BYPERP(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE DEFAULT
              WRITE (iunout,*) ' WRONG TALLY NUMBER IN ALGTAL IALV = ',
     .                           IALV
              WRITE (iunout,*) ' NO ALGBRAIC TALLY CALCULATED '
              CALL LEER(1)
              EXIT
            END SELECT

            SUMWEI = EPS60
            VEC1 = 0._DP
            DO I=1,NSBOX
              IN=NCLTAL(I)
              VEC1(IN) = VEC1(IN) + WEI(I)*OP(I)
              SUMWEI(IN) = SUMWEI(IN) + WEI(I)
            END DO
            VEC1(1:NSBOX_TAL) = VEC1(1:NSBOX_TAL)/SUMWEI(1:NSBOX_TAL)
C
C  KONSTANTE WURDE EINGELESEN
          ELSEIF (IZIF(1,IOP).LT.0) THEN
            DO 15 I=1,NSBOX_TAL
              VEC1(I)=CONST(IOP)
15          CONTINUE
C
          ELSE
C  ZWISCHENERGEBNIS HOLEN
            IF (IIND(IZIF(1,IOP)).EQ.1) THEN
              DO 20 I=1,NSBOX_TAL
                VEC1(I)=RESULT(1,I)
                IIND(IZIF(1,IOP))=0
20            CONTINUE
              LFREE1=.TRUE.
            ELSEIF (IIND(IZIF(1,IOP)).EQ.2) THEN
              DO 21 I=1,NSBOX_TAL
                VEC1(I)=RESULT(2,I)
                IIND(IZIF(1,IOP))=0
21            CONTINUE
              LFREE2=.TRUE.
            ELSE
              GOTO 92
            ENDIF
          ENDIF
C
C  2. OPERAND
C
C  TALLY HOLEN
          IF (IZIF(4,IOP).GT.0) THEN
            IF (IZIF(4,IOP).GT.NTALV) GOTO 90
            IF (.NOT.LIVTALV(IZIF(4,IOP))) GOTO 95
            IF (IZIF(3,IOP).GT.NFSTVI(IZIF(4,IOP))) GOTO 91
            DO 30 I=1,NSBOX_TAL
              VEC2(I)=ESTIMV(NADDV(IZIF(4,IOP))+IZIF(3,IOP),I)
30          CONTINUE
C
          ELSEIF (IZIF(4,IOP).LT.0) THEN
            ITL=IABS(IZIF(4,IOP))
            IF (ITL.GT.NTALI) GOTO 90
            IF (IZIF(3,IOP).GT.NFRSTP(ITL)) GOTO 91
            K=IZIF(3,IOP)
            SELECT CASE (ITL)
            CASE (1)
              OP(1:NSBOX)  = TEIN(1:NSBOX)
              WEI(1:NSBOX) = DEIN(1:NSBOX)*VOL(1:NSBOX)
            CASE (2)
              OP(1:NSBOX)  = TIIN(MPLSTI(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (3)
              OP(1:NSBOX) = DEIN(1:NSBOX)
              WEI(1:NSBOX) = VOL(1:NSBOX)
            CASE (4)
              OP(1:NSBOX) = DIIN(K,1:NSBOX)
              WEI(1:NSBOX) = VOL(1:NSBOX)
            CASE (5)
              OP(1:NSBOX) = VXIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (6)
              OP(1:NSBOX) = VYIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (7)
              OP(1:NSBOX) = VZIN(MPLSV(K),1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (8)
              OP(1:NSBOX) = BXIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (9)
              OP(1:NSBOX) = BYIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (10)
              OP(1:NSBOX) = BZIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (11)
              OP(1:NSBOX) = BFIN(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (12)
              OP(1:NSBOX) = ADIN(K,1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (13)
              OP(1:NSBOX) = EDRIFT(K,1:NSBOX)
              WEI(1:NSBOX) = DIIN(K,1:NSBOX)*VOL(1:NSBOX)
            CASE (14)
              OP(1:NSBOX) = VOL(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (15)
              OP(1:NSBOX) = WGHT(K,1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (16)
              OP(1:NSBOX) = BXPERP(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE (17)
              OP(1:NSBOX) = BYPERP(1:NSBOX)
              WEI(1:NSBOX) = 1._DP
            CASE DEFAULT
              WRITE (iunout,*) ' WRONG TALLY NUMBER IN ALGTAL IALV = ',
     .                           IALV
              WRITE (iunout,*) ' NO ALGBRAIC TALLY CALCULATED '
              CALL LEER(1)
              EXIT
            END SELECT

            SUMWEI = EPS60
            VEC2 = 0._DP
            DO I=1,NSBOX
              IN=NCLTAL(I)
              VEC2(IN) = VEC2(IN) + WEI(I)*OP(I)
              SUMWEI(IN) = SUMWEI(IN) + WEI(I)
            END DO
            VEC2(1:NSBOX_TAL) = VEC2(1:NSBOX_TAL)/SUMWEI(1:NSBOX_TAL)
C
          ELSEIF (IZIF(3,IOP).LT.0) THEN
            DO 35 I=1,NSBOX_TAL
              VEC2(I)=CONST(IOP)
35          CONTINUE
C
          ELSE
C  ZWISCHENERGEBNIS HOLEN
            IF (IIND(IZIF(3,IOP)).EQ.1) THEN
              DO 40 I=1,NSBOX_TAL
                VEC2(I)=RESULT(1,I)
                IIND(IZIF(3,IOP))=0
40            CONTINUE
              LFREE1=.TRUE.
            ELSEIF (IIND(IZIF(3,IOP)).EQ.2) THEN
              DO 41 I=1,NSBOX_TAL
                VEC2(I)=RESULT(2,I)
                IIND(IZIF(3,IOP))=0
41            CONTINUE
              LFREE2=.TRUE.
            ELSE
              GOTO 92
            ENDIF
          ENDIF
C
C
C  BERECHNE ZWISCHENERGEBNIS UND SPEICHERE AUF RESULT(II,....)
C
          IF (LFREE1) THEN
            II=1
            IIND(IOP)=1
            LFREE1=.FALSE.
          ELSEIF (LFREE2) THEN
            II=2
            IIND(IOP)=2
            LFREE2=.FALSE.
          ELSE
            GOTO 92
          ENDIF
C
          IF (OPER(IOP).EQ.'+') THEN
            DO 50 I=1,NSBOX_TAL
              RESULT(II,I)=VEC1(I)+VEC2(I)
50          CONTINUE
          ELSEIF (OPER(IOP).EQ.'-') THEN
            DO 60 I=1,NSBOX_TAL
              RESULT(II,I)=VEC1(I)-VEC2(I)
60          CONTINUE
          ELSEIF (OPER(IOP).EQ.'*') THEN
            DO 70 I=1,NSBOX_TAL
              RESULT(II,I)=VEC1(I)*VEC2(I)
70          CONTINUE
          ELSEIF (OPER(IOP).EQ.'/') THEN
            DO 81 I=1,NSBOX_TAL
              IF (VEC2(I).NE.0.D0) GOTO 82
81          CONTINUE
C  DIVISION BY ZERO TALLY. ALGEBR. TALLY IRRELEVANT. RETURN ZERO TALLY
            DO 83 I=1,NSBOX_TAL
              RESULT(II,I)=0.
83          CONTINUE
            GOTO 120
82          DO 80 I=1,NSBOX_TAL
              RESULT(II,I)=VEC1(I)/(VEC2(I)+EPS30)
80          CONTINUE
          ELSEIF (OPER(IOP).EQ.'^') THEN
            DO 85 I=1,NSBOX_TAL
              RESULT(II,I)=VEC1(I)**VEC2(I)
85          CONTINUE
          ELSE
            GOTO 93
          ENDIF
C
          GOTO 100
C
90        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) ' TALLY NUMBER OUT OF RANGE '
          WRITE (iunout,*) 
     .      ' CHECK INPUT FOR ADDITIONAL VOLUME TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          GOTO 200
C
91        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) ' SPECIES INDEX OUT OF RANGE '
          WRITE (iunout,*) 
     .      ' CHECK INPUT FOR ADDITIONAL VOLUME TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          GOTO 200
C
92        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) ' WRONG NUMBER OF INTERMEDIATE RESULT FOUND '
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
93        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) ' OPERATOR NOT FORESEEN '
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
94        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) 
     .      ' ARGUMENTS OF OPERATION HAVE DIFFERENT SPACING '
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
95        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) 
     .      ' OPERAND OF ALGEBRAIC EXPRESSION IS SWITCHED OFF'
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
C
100     CONTINUE
C
C  STORE RESULT IN ALGV
120     DO 150 J=1,NSBOX_TAL
          ALGV(IALV,J)=RESULT(II,J)
150     CONTINUE

        IF (ALLOCATED(OP)) THEN
          DEALLOCATE(OP)
          DEALLOCATE(WEI)
          DEALLOCATE(SUMWEI)
        END IF
C
200   CONTINUE
C
C
C     CALCULATE ALGEBRAIC SURFACE TALLIES
C
 300  CONTINUE

      IF (.NOT.LALGS.AND.NALSI.GT.0) THEN
        WRITE (iunout,*) ' ALGS IS SWITCHED OFF '
        WRITE (iunout,*) 
     .    ' NO ALGEBRAIC SURFACE TALLIES CAN BE CALCULATED '
        RETURN
      END IF

!pb      IF (NLIMPS.GT.NRAD) GOTO 999

      IF (NALSI > 0) THEN
        IF (.NOT.ALLOCATED(VEC1)) THEN
          ALLOCATE (VEC1(NLIMPS))
          ALLOCATE (VEC2(NLIMPS))
          ALLOCATE (RESULT(2,NLIMPS))
        END IF
        ALLOCATE (LLIMPS(NLIMPS))
      END IF
C
      DO 500 IALS=1,NALSI
C
        HCHR=CHRTLS(IALS)
C
        CALL ALGEBR (HCHR,OPER,IZIF,CONST,NOP)
C
        LLMPS=.FALSE.
        DO J=1,NLIMPS
          VEC1(J)=0.D0
          VEC2(J)=0.D0
          LLIMPS(J)=.FALSE.
        ENDDO
        DO 301 IOP=1,NOP
          IIND(IOP)=0
C         WRITE (iunout,*) IOP,OPER(IOP),(IZIF(J,IOP),J=1,4)
301    CONTINUE
       LFREE1=.TRUE.
       LFREE2=.TRUE.
C
        DO 400 IOP=1,NOP
C
C  1. OPERAND
C
C  TALLY HOLEN
          IF (IZIF(2,IOP).GT.0) THEN
            IF (IZIF(2,IOP).GT.NTALS) GOTO 390
            IF (.NOT.LIVTALS(IZIF(2,IOP))) GOTO 395
            IF (IZIF(1,IOP).GT.NFRSTW(IZIF(2,IOP))) THEN
              IZIF(1,IOP)=IZIF(1,IOP)-NFRSTW(IZIF(2,IOP))
              IF (IZIF(1,IOP).GT.NFRSTW(IZIF(2,IOP))*NLIMPS) GOTO 391
              ILIMPS=(IZIF(1,IOP)-1)/NFRSTW(IZIF(2,IOP))+1
              ISPZ=IZIF(1,IOP)-(ILIMPS-1)*NFRSTW(IZIF(2,IOP))
              IINDEX=NADDW(IZIF(2,IOP))*NLMPGS+IZIF(1,IOP)+NESTM1
              IZIF(1,IOP)=ISPZ
              VEC1(1)=ESTIMS(NADDW(IZIF(2,IOP))+IZIF(1,IOP),ILIMPS)
              LLIMPS(ILIMPS)=.TRUE.
              LLMPS=.TRUE.
            ELSE
            DO 310 I=1,NLIMPS
              IINDEX=NADDW(IZIF(2,IOP))*NLMPGS+(I-1)*NFRSTW(IZIF(2,IOP))
     .              +IZIF(1,IOP)+NESTM1
              VEC1(I)=ESTIMS(NADDW(IZIF(2,IOP))+IZIF(1,IOP),I)
310         CONTINUE
            ENDIF
C
C
C  KONSTANTE WURDE EINGELESEN
          ELSEIF (IZIF(1,IOP).LT.0) THEN
            DO 315 I=1,NLIMPS
              VEC1(I)=CONST(IOP)
315          CONTINUE
C
          ELSE
C  ZWISCHENERGEBNIS HOLEN
            IF (IIND(IZIF(1,IOP)).EQ.1) THEN
              DO 320 I=1,NLIMPS
                VEC1(I)=RESULT(1,I)
                IIND(IZIF(1,IOP))=0
320           CONTINUE
              LFREE1=.TRUE.
            ELSEIF (IIND(IZIF(1,IOP)).EQ.2) THEN
              DO 321 I=1,NLIMPS
                VEC1(I)=RESULT(2,I)
                IIND(IZIF(1,IOP))=0
321           CONTINUE
              LFREE2=.TRUE.
            ELSE
              GOTO 392
            ENDIF
          ENDIF
C
C  2. OPERAND
C
C  TALLY HOLEN
          IF (IZIF(4,IOP).GT.0) THEN
            IF (IZIF(4,IOP).GT.NTALS) GOTO 390
            IF (.NOT.LIVTALS(IZIF(4,IOP))) GOTO 395
            IF (IZIF(3,IOP).GT.NFRSTW(IZIF(4,IOP))) THEN
              IZIF(3,IOP)=IZIF(3,IOP)-NFRSTW(IZIF(4,IOP))
              IF (IZIF(3,IOP).GT.NFRSTW(IZIF(4,IOP))*NLIMPS) GOTO 391
              ILIMPS=(IZIF(3,IOP)-1)/NFRSTW(IZIF(4,IOP))+1
              ISPZ=IZIF(3,IOP)-(ILIMPS-1)*NFRSTW(IZIF(4,IOP))
              IINDEX=NADDW(IZIF(4,IOP))*NLMPGS+IZIF(3,IOP)+NESTM1
              IZIF(3,IOP)=ISPZ
              VEC2(1)=ESTIMS(NADDW(IZIF(4,IOP))+IZIF(3,IOP),ILIMPS)
              LLIMPS(ILIMPS)=.TRUE.
              LLMPS=.TRUE.
            ELSE
            DO 330 I=1,NLIMPS
              IINDEX=NADDW(IZIF(4,IOP))*NLMPGS+(I-1)*NFRSTW(IZIF(4,IOP))
     .              +IZIF(3,IOP)+NESTM1
              VEC2(I)=ESTIMS(NADDW(IZIF(4,IOP))+IZIF(3,IOP),I)
330         CONTINUE
            ENDIF
C
C
C  KONSTANTE WURDE EINGELESEN
          ELSEIF (IZIF(3,IOP).LT.0) THEN
            DO 335 I=1,NLIMPS
              VEC2(I)=CONST(IOP)
335         CONTINUE
C
          ELSE
C  ZWISCHENERGEBNIS HOLEN
            IF (IIND(IZIF(3,IOP)).EQ.1) THEN
              DO 340 I=1,NLIMPS
                VEC2(I)=RESULT(1,I)
                IIND(IZIF(3,IOP))=0
340            CONTINUE
              LFREE1=.TRUE.
            ELSEIF (IIND(IZIF(3,IOP)).EQ.2) THEN
              DO 341 I=1,NLIMPS
                VEC2(I)=RESULT(2,I)
                IIND(IZIF(3,IOP))=0
341           CONTINUE
              LFREE2=.TRUE.
            ELSE
              GOTO 392
            ENDIF
          ENDIF
C
C  BERECHNE ZWISCHENERGEBNIS UND SPEICHERE AUF RESULT(II,....)
C
          IF (LFREE1) THEN
            II=1
            IIND(IOP)=1
            LFREE1=.FALSE.
          ELSEIF (LFREE2) THEN
            II=2
            IIND(IOP)=2
            LFREE2=.FALSE.
          ELSE
            GOTO 392
          ENDIF
C
          IF (OPER(IOP).EQ.'+') THEN
            DO 350 I=1,NLIMPS
              RESULT(II,I)=VEC1(I)+VEC2(I)
350         CONTINUE
          ELSEIF (OPER(IOP).EQ.'-') THEN
            DO 360 I=1,NLIMPS
              RESULT(II,I)=VEC1(I)-VEC2(I)
360         CONTINUE
          ELSEIF (OPER(IOP).EQ.'*') THEN
            DO 370 I=1,NLIMPS
              RESULT(II,I)=VEC1(I)*VEC2(I)
370         CONTINUE
          ELSEIF (OPER(IOP).EQ.'/') THEN
            DO 381 I=1,NLIMPS
              IF (VEC2(I).NE.0.D0) GOTO 382
381         CONTINUE
C  DIVISION BY ZERO TALLY. ALGEBR. TALLY IRRELEVANT. RETURN ZERO TALLY
            DO 383 I=1,NLIMPS
              RESULT(II,I)=0.
383         CONTINUE
            GOTO 420
382         DO 380 I=1,NLIMPS
              RESULT(II,I)=VEC1(I)/(VEC2(I)+EPS30)
380         CONTINUE
          ELSEIF (OPER(IOP).EQ.'^') THEN
            DO 385 I=1,NLIMPS
              RESULT(II,I)=VEC1(I)**VEC2(I)
385         CONTINUE
          ELSE
            GOTO 393
          ENDIF
C
          GOTO 400
C
390       CONTINUE
          WRITE (iunout,*) ' TALLY NUMBER OUT OF RANGE '
          WRITE (iunout,*) 
     .      ' CHECK INPUT FOR ADDITIONAL SURFACE TALLY NO. ',IALS
          WRITE (iunout,*) CHRTLS(IALS)
          GOTO 500
C
391       CONTINUE
          WRITE (iunout,*) ' SPECIES INDEX OUT OF RANGE '
          WRITE (iunout,*) 
     .      ' CHECK INPUT FOR ADDITIONAL SURFACE TALLY NO. ',IALS
          WRITE (iunout,*) CHRTLS(IALS)
          GOTO 500
C
392       CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGEBR '
          WRITE (iunout,*) ' WRONG NUMBER OF INTERMEDIATE RESULT FOUND '
          WRITE (iunout,*) CHRTLS(IALS)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 500
393       CONTINUE
C
          WRITE (iunout,*) ' OPERATOR NOT FORESEEN '
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALS
          WRITE (iunout,*) CHRTLS(IALS)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 500
C
395       CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE ALGTAL '
          WRITE (iunout,*) 
     .      ' OPERAND OF ALGEBRAIC EXPRESSION IS SWITCHED OFF'
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALS
          WRITE (iunout,*) CHRTLS(IALS)
          WRITE (iunout,'(1X,A,4I4)') 
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 500
C
C
400     CONTINUE
C
C  STORE RESULT IN ALGS
420     CONTINUE
        IF (LLMPS) THEN
          DO 440 J=1,NLIMPS
            IF (LLIMPS(J)) ALGS(IALS,J)=RESULT(II,1)
440       CONTINUE
        ELSE
          DO 450 J=1,NLIMPS
            ALGS(IALS,J)=RESULT(II,J)
450       CONTINUE
        ENDIF
C
500   CONTINUE

      DEALLOCATE (VEC1)
      DEALLOCATE (VEC2)
      DEALLOCATE (RESULT)
      IF (NALSI > 0) DEALLOCATE (LLIMPS)
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'STORAGE CONFLICT IN ALGTAL, BECAUSE NRAD<NLIMPS'
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: halpha.f
*DK HALPHA
      SUBROUTINE HALPHA (IST,IAD1,IAD2,IAD3,IAD4,IAD5,IADS,IADF)
c
C  SUBROUTINE FOR H-ALPHA EMISSIVITY. (BALMER SERIES)
C  CALLED FROM EIRENE, SECTION DIAGNO, SUBR. SIGAL
C  THE H-ALPHA EMISSIVITY PROFILE (PHOTONS/S/CM**3) IS COMPUTED
C  AND WRITTEN ONTO TALLIES ADDV(IAD1,...),... FOR STRATUM NO. IST
C  IAD1: CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  IAD2: CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  IAD3: CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  IAD4: CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  IAD5: CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  IADF: FULCHER
C  IADS: SUM OVER ALL CONTRINUTIONS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE COMSIG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMPRT
      USE COMSOU
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI
c slmod begin
      USE CTETRA
c slmod end

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IAD1, IAD2, IAD3, IAD4, IAD5, IADS, IST,
     .                       IADF
      REAL(DP) :: DA31(0:8,0:8),
     .            DP31(0:8,0:8),
     .            DM31(0:8,0:8),
     .            DI31(0:8,0:8),
     .            DN31(0:8,0:8),
     .            DNFUL(0:8,0:8)
      REAL(DP) :: RHMH2(0:8), RH2PH2(0:8,0:8), FP(6)
      REAL(DP) :: DATM3, DION3, DMOL3, DPLS3, DNML3, RATIO7, TEI, DEJ,
     .          SIGADD1, SIGADD2, SIGADD3, SIGADD4, SIGADD5, TEF, DEF,
     .          RHMP2, DA, DM,
     .          DNFU3, POWALFF,
     .          DPP, DI, DN, RATIO2, SIGADD,
     .          FAC21, FAC31, FAC41, FAC51, FAC61,
     .          FAC32, FAC42, FAC52, FAC62,
     .          FAC43, FAC53, FAC63,
     .          POWALF, POWALF1, POWALF2, POWALF3, POWALF4,
     .          POWALF5, DE, TE, FACFUL, SIGADDF, RCMIN, RCMAX
      INTEGER :: IRC, IFIRST, NCELC, IERROR, IR, I, J, JFEXMN, JFEXMX
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(6) :: CISTRA
C
      SAVE
C
      DATA IFIRST/0/
C  RADIATIVE TRANSITION RATES (1/S)
C  BALMER ALPHA
      FAC32=4.410E7
C  BALMER BETA
      FAC42=8.419E6
C  BALMER GAMMA
      FAC52=2.530E6
C  BALMER DELTA
      FAC62=9.732E5
C
C  LYMAN ALPHA
      FAC21=4.699E8
C  LYMAN BETA
      FAC31=5.575E7
C  LYMAN GAMMA
      FAC41=1.278E7
C  LYMAN DELTA
      FAC51=4.125E6
C  LYMAN EPSILON
      FAC61=1.644E6
C
C  PASCHEN ALPHA
      FAC43=8.986E6
C  PASCHEN BETA
      FAC53=2.201E6
C  PASCHEN GAMMA
      FAC63=7.783E5
C
C  FULCHER
      FACFUL=2.53E7
C
C  INITIALIZE ATOMIC DATA ARRAYS
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        IERROR=0
C
C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DPP,DM,DI, AND DN ARRAY
C
        IR=NREACI+1
        IF (IR.GT.NREAC) THEN
          WRITE (iunout,*) 'FROM SUBROUTINE HALPHA: '
          CALL MASPRM('NREAC',5,NREAC,'IR',2,IR,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
C
        FILNAM='AMJUEL  '
        H123='H.12'
        CRC='OT '
        FP = 0._DP
        RCMIN = -HUGE(1._DP)
        RCMAX =  HUGE(1._DP)
        JFEXMN = 0
        JFEXMX = 0
C
C  H(n=3)/H(n=1)
        REAC='2.1.5a   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DA31(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=3)/H+
        REAC='2.1.8a   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DP31(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=3)/H2(g)
        REAC='2.2.5a   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DM31(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=3)/H2+(g)
        REAC='2.2.14a  '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DI31(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=3)/H-
        REAC='7.2a     '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DN31(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H2(n=3,D)/nH2
        REAC='2.2.5fu  '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DNFUL(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C
C  NOW READ RATIO OF DENSITIES:
C
C  FIRST: H-/H2
        FILNAM='AMJUEL  '
        H123='H.11'
        REAC='7.0a     '
        CRC='OT '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO I=1,9
          RHMH2(I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(I,1)
        ENDDO

C  NEXT : H2+/H2
        FILNAM='AMJUEL  '
        H123='H.12'
        REAC='2.0c     '
        CRC='OT '
C  2.0c INCLUDES ION CONVERION (CX) ON H2(V)
C  OLD VERSION (WITHOUT THIS CX) SHOULD BE RECOVERED BY
C  READING 2.0b INSTEAD, AND OMITTING THE H- CHANNEL 5.
C
CDR: IF CX ON H2(V) IS NOT INCLUDED IN A SPECIFIC NEUTRAL TRANSPORT EQUATION
CDR  (SEE INPUT BLOCK 4), THEN IT SHOULD NOT BE INCLUDED HERE EITHER.
C       REAC='2.0b    '
CDR
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO I=1,9
          DO J=1,9
            RH2PH2(I-1,J-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(I,J)
          ENDDO
        ENDDO
      ENDIF
C
C  END OF INITIALIZATION
C
      IF (IESTR.EQ.IST) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (iunout,*) 'ERROR IN HALPHA: DATA FOR STRATUM ISTRA= ',IST
        WRITE (iunout,*) 'ARE NOT AVAILABLE. H-ALPHA ABANDONNED'
        RETURN
      ENDIF
C
C  LOOP OVER 2D MESH
C
      POWALF=0.
      POWALF1=0.
      POWALF2=0.
      POWALF3=0.
      POWALF4=0.
      POWALF5=0.
      POWALFF=0.

      IF (MAX(IAD1,IAD2,IAD3,IAD4,IAD5,IADS,IADF) > NADV) GOTO 999
      ADDV(IAD1,1:NRTAL) = 0.D0
      ADDV(IAD2,1:NRTAL) = 0.D0
      ADDV(IAD3,1:NRTAL) = 0.D0
      ADDV(IAD4,1:NRTAL) = 0.D0
      ADDV(IAD5,1:NRTAL) = 0.D0
      ADDV(IADF,1:NRTAL) = 0.D0
      ADDV(IADS,1:NRTAL) = 0.D0
C
      DO 1000 NCELL=1,NSBOX
C
C  LOCAL PLASMA DATA
C
        IF (NSTGRD(NCELL) > 0) CYCLE

        IF (LGVAC(NCELL,NPLS+1)) THEN
          TE=TVAC
          DE=DVAC
          NCELC=NCLTAL(NCELL)
        ELSE
          TE=TEIN(NCELL)
          DE=DEIN(NCELL)
          NCELC=NCLTAL(NCELL)
        ENDIF
!pb        IF (NCELC <= 0) CYCLE
       
C
        SIGADD1=0.
        SIGADD2=0.
        SIGADD3=0.
        SIGADD4=0.
        SIGADD5=0.
        SIGADDF=0.
        IF (LGVAC(NCELL,NPLS+1)) GOTO 500
        DEF=LOG(DE*1.D-8)
        TEF=LOG(TE)
        DATM3=0.
        DPLS3=0.
        DMOL3=0.
        DION3=0.
        DNML3=0.
        DNFU3=0.
        DO 150 J=0,8
          DEJ=DEF**J
          DO 150 I=0,8
            TEI=TEF**I
            DATM3=DATM3+DA31(I,J)*TEI*DEJ
            DPLS3=DPLS3+DP31(I,J)*TEI*DEJ
            DMOL3=DMOL3+DM31(I,J)*TEI*DEJ
            DION3=DION3+DI31(I,J)*TEI*DEJ
            DNML3=DNML3+DN31(I,J)*TEI*DEJ
            DNFU3=DNFU3+DNFUL(I,J)*TEI*DEJ
150     CONTINUE
        DATM3=EXP(DATM3)
        DPLS3=EXP(DPLS3)
        DMOL3=EXP(DMOL3)
        DION3=EXP(DION3)
        DNML3=EXP(DNML3)
        DNFU3=EXP(DNFU3)

C  RATIO OF DENSITIES: H- TO H2, COLL. EQUIL. IN VIBRATION
C  (ONLY TE-DEPENDENT)

        RATIO7=0
        DO 160 I=0,8
          TEI=TEF**I
          RATIO7=RATIO7+RHMH2(I)*TEI
160     CONTINUE
        RATIO7=EXP(RATIO7)

C  RATIO OF DENSITIES: H2+ TO H2, INCL. ION CONVERSION
C  RATIO OF DENSITIES: H2+ TO H2, EXCL. ION CONVERSION

        RATIO2=0
        DO 170 J=0,8
          DEJ=DEF**J
          DO 170 I=0,8
            TEI=TEF**I
            RATIO2=RATIO2+RH2PH2(I,J)*TEI*DEJ
170     CONTINUE
        RATIO2=EXP(RATIO2)
C
C  CHANNEL 1
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENA (IONIZATION)
C
        DO 200 IATM=1,NATMI
          IF (NCHARA(IATM).NE.1) GOTO 200
          DA=DATM3*PDENA(IATM,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD1=SIGADD1+DA*FAC32
200     CONTINUE
C
C  CHANNEL 2
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN DIIN (RECOMBINATION)
C
        DO 205 IPLS=1,NPLSI
          IF (NCHARP(IPLS).NE.1.OR.NCHRGP(IPLS).NE.1) GOTO 205
          DPP=DPLS3*DIIN(IPLS,NCELL)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD2=SIGADD2+DPP*FAC32
205     CONTINUE
C
C  CHANNEL 3
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENM (DISSOCIATION OF H2)
C
        DO 210 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 210
          DM=DMOL3*PDENM(IMOL,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD3=SIGADD3+DM*FAC32
210     CONTINUE
C
C  CHANNEL 4
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENI: (DISSOCIATION OF H2+)
C
C  REVISED: USE (PDENM * DENSITY RATIO H2+/H2) NOW, INSTEAD OF PDENI
C
C       DO 215 IION=1,NIONI
C         IF (NCHARI(IION).NE.2) GOTO 215
C         DI=DION3*PDENI(IION,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
C         SIGADD4=SIGADD4+DI*FAC32
C215     CONTINUE
        DO 215 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 215
          DI=DION3*PDENM(IMOL,NCELC)*RATIO2
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD4=SIGADD4+DI*FAC32
215     CONTINUE
C
C  CHANNEL 5
C  H ALPHA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN H- DENSITY (CHARGE EXCHANGE RECOMBINATION)
C
C  REVISED: USE (PDENM * DENSITY RATIO H-/H2) NOW, INSTEAD OF PDENI
C
        DO 220 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 220
          DN=DNML3*PDENM(IMOL,NCELC)*RATIO7
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD5=SIGADD5+DN*FAC32
220     CONTINUE
C
C
C  CHANNEL 6
C  FULCHER SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENM (FULCHER)
C
        DO IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) CYCLE
          DM=DNFU3*PDENM(IMOL,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 3-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADDF=SIGADDF+DM*FACFUL
        END DO
C
500     CONTINUE
C
        SIGADD=SIGADD1+SIGADD2+SIGADD3+SIGADD4+SIGADD5
C
        ADDV(IAD1,NCELC)=ADDV(IAD1,NCELC)+SIGADD1*VOL(NCELL)
        ADDV(IAD2,NCELC)=ADDV(IAD2,NCELC)+SIGADD2*VOL(NCELL)
        ADDV(IAD3,NCELC)=ADDV(IAD3,NCELC)+SIGADD3*VOL(NCELL)
        ADDV(IAD4,NCELC)=ADDV(IAD4,NCELC)+SIGADD4*VOL(NCELL)
        ADDV(IAD5,NCELC)=ADDV(IAD5,NCELC)+SIGADD5*VOL(NCELL)
        ADDV(IADF,NCELC)=ADDV(IADF,NCELC)+SIGADDF*VOL(NCELL)
        ADDV(IADS,NCELC)=ADDV(IADS,NCELC)+SIGADD*VOL(NCELL)
C
C
C       FACT=E0032*ELCHA = 3.028E-19
        POWALF1=POWALF1+SIGADD1*3.028E-19*VOL(NCELL)
        POWALF2=POWALF2+SIGADD2*3.028E-19*VOL(NCELL)
        POWALF3=POWALF3+SIGADD3*3.028E-19*VOL(NCELL)
        POWALF4=POWALF4+SIGADD4*3.028E-19*VOL(NCELL)
        POWALF5=POWALF5+SIGADD5*3.028E-19*VOL(NCELL)
        POWALFF=POWALFF+SIGADDF*3.028E-19*VOL(NCELL)
C
        POWALF =POWALF +SIGADD *3.028E-19*VOL(NCELL)
C
1000  CONTINUE

      ADDV(IAD1,1:NSBOX_TAL)=ADDV(IAD1,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD2,1:NSBOX_TAL)=ADDV(IAD2,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD3,1:NSBOX_TAL)=ADDV(IAD3,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD4,1:NSBOX_TAL)=ADDV(IAD4,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IAD5,1:NSBOX_TAL)=ADDV(IAD5,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IADF,1:NSBOX_TAL)=ADDV(IADF,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)
      ADDV(IADS,1:NSBOX_TAL)=ADDV(IADS,1:NSBOX_TAL)/VOLTAL(1:NSBOX_TAL)

      CALL LEER(2)
      CALL FTCRI(IST,CISTRA)
      IF (IST.GT.0)
     .CALL MASBOX('SUBR. HALPHA CALLED, FOR STRATUM NO. '//CISTRA)
      IF (IST.EQ.0)
     .CALL MASBOX('SUBR. HALPHA CALLED, FOR SUM OVER STRATA')
      CALL LEER(1)
      WRITE (iunout,*) ' AFTER INTEGRATION OVER COMPUTATIONAL DOMAIN'
      WRITE (iunout,*) ' TOTAL RADIATED POWER (WATT) BY HALPHA:',POWALF
      WRITE (iunout,*) ' COUPL. TO GROUNDSTATE                :',POWALF1
      WRITE (iunout,*) ' COUPLING TO CONTINUUM                :',POWALF2
      WRITE (iunout,*) ' COUPLING TO MOLECULES                :',POWALF3
      WRITE (iunout,*) ' COUPLING TO MOL.IONS                 :',POWALF4
      WRITE (iunout,*) ' COUPLING TO NEG.IONS                 :',POWALF5
      WRITE (iunout,*) ' COUPLING TO FULCHER                  :',POWALFF
      CALL LEER(2)

c slmod begin
c ... Calling INTTAL causes some cases with tetrahedrons to crash (i-tet-000a), sans error
c     message, but others are fine (i-lfs-0002b).  There's no response from WRITE(0 messages
c     in INTTAL at all, so the code is dying somewhere in the netherverse between the 
c     CALL and SUBROUTINE.  I checked the variable being passed and they all seem to be 
c     properly declared.  I had this problem on the Culham computes as well, but that as also
c     the Intel compiler.  I didn't compile with any additional debugging flags (just bounds
c     checking, as standard) and run the case again. The integrals aren't required but OSM/
c     DIVIMP, so just avoiding this code (in the HGAMMA routine as well). - SL, 12/03/2010
c      WRITE(0,*) 'DIM 1',SIZE(addv,1)
c      WRITE(0,*) 'DIM 2',SIZE(addv,2)
c      WRITE(0,*) 'NADV ',nadv
c      WRITE(0,*) 'SIZE VOL ',SIZE(vol)
c      WRITE(0,*) 'IAD1,NADV,NSBOX_TAL',IAD1,NADV,NSBOX_TAL
c      WRITE(0,*) 'NSBOX,NR1TAL,NP2TAL,NT3TAL,NBMLT',
c     .           NR1TAL,NP2TAL,NT3TAL,NBMLT
      IF (NLTET) THEN
        WRITE(0,*) 'WARNING HALPHA: NOT CALLING INTTAL DUE TO '//
     .             'COMPILER BUG (?)'
      ELSE 
        CALL INTTAL (ADDV,VOLTAL,IAD1,NADV,NSBOX_TAL,ADDVI(IAD1,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IAD2,NADV,NSBOX_TAL,ADDVI(IAD2,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IAD3,NADV,NSBOX_TAL,ADDVI(IAD3,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IAD4,NADV,NSBOX_TAL,ADDVI(IAD4,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IAD5,NADV,NSBOX_TAL,ADDVI(IAD5,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IADF,NADV,NSBOX_TAL,ADDVI(IADF,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
        CALL INTTAL (ADDV,VOLTAL,IADS,NADV,NSBOX_TAL,ADDVI(IADS,IST),
     .               NR1TAL,NP2TAL,NT3TAL,NBMLT)
       ENDIF
c
c      CALL INTTAL (ADDV,VOLTAL,IAD1,NADV,NSBOX_TAL,ADDVI(IAD1,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IAD2,NADV,NSBOX_TAL,ADDVI(IAD2,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IAD3,NADV,NSBOX_TAL,ADDVI(IAD3,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IAD4,NADV,NSBOX_TAL,ADDVI(IAD4,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IAD5,NADV,NSBOX_TAL,ADDVI(IAD5,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IADF,NADV,NSBOX_TAL,ADDVI(IADF,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c      CALL INTTAL (ADDV,VOLTAL,IADS,NADV,NSBOX_TAL,ADDVI(IADS,IST),
c     .             NR1TAL,NP2TAL,NT3TAL,NBMLT)
c slmod end
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ENDIF
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. HALPHA '
      WRITE (iunout,*) 'NO STORAGE AVAILBALE ON ADDITIONAL TALLY ADDV '
      WRITE (iunout,*) 'STORAGE REQUESTED FOR IADV= ',
     .             IAD1,IAD2,IAD3,IAD4,
     .             IAD5,IADS,IADF
      WRITE (iunout,*) 'CHECK INPUT BLOCK 10A '
      CALL EXIT_OWN(1)
      END
C ===== SOURCE: hgamma.f
!pb  24.11.06: mistyped SGSM_BGK changed to SGMS_BGK
 
*DK HGAMMA
      SUBROUTINE HGAMMA (IST,IAD1,IAD2,IAD3,IAD4,IAD5,IADS)
c
C  SUBROUTINE FOR H-GAMMA EMISSIVITY. (BALMER SERIES)
C  CALLED FROM EIRENE, SECTION DIAGNO, SUBR. SIGAL
C  THE H-GAMMA EMISSIVITY PROFILE (PHOTONS/S/CM**3) IS COMPUTED
C  AND WRITTEN ONTO TALLIES ADDV(IAD1,...),... FOR STRATUM NO. IST
C  IAD1: CONTRIBUTION LINEAR IN H   -ATOM      DENSITY
C  IAD2: CONTRIBUTION LINEAR IN H+  -ION       DENSITY
C  IAD3: CONTRIBUTION LINEAR IN H2  -MOLEC.    DENSITY
C  IAD4: CONTRIBUTION LINEAR IN H2+ -MOLEC.ION DENSITY
C  IAD5: CONTRIBUTION LINEAR IN H-  -NEG. ION  DENSITY
C  IADF: ....
C  IADS: SUM OVER ALL CONTRINUTIONS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CUPD
      USE COMSIG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMPRT
      USE COMSOU
      USE CLGIN
      USE COUTAU
      USE COMXS
      USE CSPEI
c slmod begin
      USE CTETRA
c slmod end

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: IST, IAD1, IAD2, IAD3, IAD4, IAD5, IADS

      REAL(DP) :: DA51(0:8,0:8)
      REAL(DP) :: DP51(0:8,0:8)
      REAL(DP) :: DM51(0:8,0:8)
      REAL(DP) :: DI51(0:8,0:8)
      REAL(DP) :: DN51(0:8,0:8)

      REAL(DP) :: RHMH2(0:8), RH2PH2(0:8,0:8), FP(6)
      REAL(DP) :: DNML5, DION5, DMOL5, DEJ, RATIO2, RATIO7, TEI, DPLS5,
     .          SIGADD1, SIGADD2, SIGADD3, SIGADD4, SIGADD5, DATM5, TEF,
     .          DEF, DA, DPP, DM, DI, SGSM_BGK, SIGADD, DN, FAC21,
     .          FAC41, FAC51, FAC32, FAC42, FAC52, POWALF1, POWALF2,
     .          POWALF3, POWALF4, POWALF5, DE, TE, FAC43, FAC53, POWALF,
     .          FAC31, RCMIN, RCMAX
      INTEGER :: IRC, IFIRST, NCELC, IERROR, IR, I, J, JFEXMN, JFEXMX
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      CHARACTER(8) :: FILNAM
      CHARACTER(4) :: H123
      CHARACTER(9) :: REAC
      CHARACTER(3) :: CRC
      CHARACTER(6) :: CISTRA
C
      SAVE
C
      DATA IFIRST/0/
C  RADIATIVE TRANSITION RATES (1/S)
C  BALMER ALPHA
      FAC32=4.410E7
C  BALMER BETA
      FAC42=8.419E6
C  BALMER GAMMA
      FAC52=2.530E6
C
C  LYMAN ALPHA
      FAC21=4.699E8
C  LYMAN BETA
      FAC31=5.575E7
C  LYMAN GAMMA
      FAC41=1.278E7
C  LYMAN DELTA
      FAC51=4.125E6
C
C  PASCHEN ALPHA
      FAC43=8.986E6
C  PASCHEN BETA
      FAC53=2.201E6
C
C  INITIALIZE ATOMIC DATA ARRAYS
C
      IF (IFIRST.EQ.0) THEN
        IFIRST=1
        IERROR=0
C
C  READ REDUCED POPULATION COEFFICIENT FOR HYDR. ATOMS FROM FILE AMJUEL
C  AND PUT THEM FROM CREAC(..,..,IR) ONTO DA,DP,DM,DI, AND DN ARRAY
C
        IR=NREACI+1
        IF (IR.GT.NREAC) THEN
          WRITE (iunout,*) 'FROM SUBROUTINE HGAMMA: '
          CALL MASPRM('NREAC',5,NREAC,'IR',2,IR,IERROR)
          CALL EXIT_OWN(1)
        ENDIF
C
        FILNAM='AMJUEL  '
        H123='H.12'
        CRC='OT '
C
C  H(n=5)/H(n=1)
        REAC='2.1.5d   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DA51(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=5)/H+
        REAC='2.1.8d   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DP51(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=5)/H2(g)
        REAC='2.2.5d   '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DM51(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=5)/H2+(g)
        REAC='2.2.14d  '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DI51(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C  H(n=5)/H-
        REAC='7.2d     '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO J=1,9
          DO I=1,9
            DN51(J-1,I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(J,I)
          ENDDO
        ENDDO
C
C  NOW READ RATIO OF DENSITIES:
C
C  FIRST: H-/H2
        FILNAM='AMJUEL  '
        H123='H.11'
        REAC='7.0a     '
        CRC='OT '
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO I=1,9
          RHMH2(I-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(I,1)
        ENDDO

C  NEXT : H2+/H2
        FILNAM='AMJUEL  '
        H123='H.12'
        REAC='2.0c     '
        CRC='OT '
C  2.0c INCLUDES ION CONVERION (CX) ON H2(V)
C  OLD VERSION (WITHOUT THIS CX) SHOULD BE RECOVERED BY
C  READING 2.0b INSTEAD, AND OMITTING THE H- CHANNEL 5.
C
CDR: IF CX ON H2(V) IS NOT INCLUDED IN A SPECIFIC NEUTRAL TRANSPORT EQUATION
CDR  (SEE INPUT BLOCK 4), THEN IT SHOULD NOT BE INCLUDED HERE EITHER.
c slmod begin
c        WRITE(0,*) 'NOTE: REAC=2.0b IS OFF, ASK DETLEV'
c        REAC='2.0b    '
c slmod end
CDR
        REACDAT(NREACI+1)%LOTH = .FALSE.
        CALL SLREAC(NREACI+1,FILNAM,H123,REAC,CRC,
     .              RCMIN, RCMAX, FP, JFEXMN, JFEXMX,'  ',0)
        DO I=1,9
          DO J=1,9
            RH2PH2(I-1,J-1)=REACDAT(NREACI+1)%OTH%POLY%DBLPOL(I,J)
          ENDDO
        ENDDO
      ENDIF
C
C  END OF INITIALIZATION
C
      IF (IESTR.EQ.IST) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL RSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
      ELSE
        WRITE (iunout,*) 'ERROR IN HGAMMA: DATA FOR STRATUM ISTRA= ',IST
        WRITE (iunout,*) 'ARE NOT AVAILABLE. H-GAMMA ABANDONNED'
        RETURN
      ENDIF
C
C  LOOP OVER 2D MESH
C
      POWALF=0.
      POWALF1=0.
      POWALF2=0.
      POWALF3=0.
      POWALF4=0.
      POWALF5=0.

      IF (MAX(IAD1,IAD2,IAD3,IAD4,IAD5,IADS     ) > NADV) GOTO 999
      ADDV(IAD1,1:NRTAL) = 0.D0
      ADDV(IAD2,1:NRTAL) = 0.D0
      ADDV(IAD3,1:NRTAL) = 0.D0
      ADDV(IAD4,1:NRTAL) = 0.D0
      ADDV(IAD5,1:NRTAL) = 0.D0

      ADDV(IADS,1:NRTAL) = 0.D0
C
      DO 1000 NCELL=1,NSBOX
C
C  LOCAL PLASMA DATA
C
        IF (LGVAC(NCELL,NPLS+1)) THEN
          TE=TVAC
          DE=DVAC
          NCELC=NCLTAL(NCELL)
        ELSE
          TE=TEIN(NCELL)
          DE=DEIN(NCELL)
          NCELC=NCLTAL(NCELL)
        ENDIF
C
        SIGADD1=0.
        SIGADD2=0.
        SIGADD3=0.
        SIGADD4=0.
        SIGADD5=0.

        IF (LGVAC(NCELL,NPLS+1)) GOTO 500
        DEF=LOG(DE*1.D-8)
        TEF=LOG(TE)
        DATM5=0.
        DPLS5=0.
        DMOL5=0.
        DION5=0.
        DNML5=0.

        DO 150 J=0,8
          DEJ=DEF**J
          DO 150 I=0,8
            TEI=TEF**I
            DATM5=DATM5+DA51(I,J)*TEI*DEJ
            DPLS5=DPLS5+DP51(I,J)*TEI*DEJ
            DMOL5=DMOL5+DM51(I,J)*TEI*DEJ
            DION5=DION5+DI51(I,J)*TEI*DEJ
            DNML5=DNML5+DN51(I,J)*TEI*DEJ

150     CONTINUE
        DATM5=EXP(DATM5)
        DPLS5=EXP(DPLS5)
        DMOL5=EXP(DMOL5)
        DION5=EXP(DION5)
        DNML5=EXP(DNML5)


C  RATIO OF DENSITIES: H- TO H2, COLL. EQUIL. IN VIBRATION
C  (ONLY TE-DEPENDENT)

        RATIO7=0
        DO 160 I=0,8
          TEI=TEF**I
          RATIO7=RATIO7+RHMH2(I)*TEI
160     CONTINUE
        RATIO7=EXP(RATIO7)

C  RATIO OF DENSITIES: H2+ TO H2, INCL. ION CONVERSION
C  RATIO OF DENSITIES: H2+ TO H2, EXCL. ION CONVERSION

        RATIO2=0
        DO 170 J=0,8
          DEJ=DEF**J
          DO 170 I=0,8
            TEI=TEF**I
            RATIO2=RATIO2+RH2PH2(I,J)*TEI*DEJ
170     CONTINUE
        RATIO2=EXP(RATIO2)
C
C  CHANNEL 1
C  H GAMMA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENA (IONIZATION)
C
        DO 200 IATM=1,NATMI
          IF (NCHARA(IATM).NE.1) GOTO 200
          DA=DATM5*PDENA(IATM,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD1=SIGADD1+DA*FAC52
200     CONTINUE
C
C  CHANNEL 2
C  H GAMMA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN DIIN (RECOMBINATION)
C
        DO 205 IPLS=1,NPLSI
          IF (NCHARP(IPLS).NE.1.OR.NCHRGP(IPLS).NE.1) GOTO 205
          DPP=DPLS5*DIIN(IPLS,NCELL)
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD2=SIGADD2+DPP*FAC52
205     CONTINUE
C
C  CHANNEL 3
C  H GAMMA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENM (DISSOCIATION OF H2)
C
        DO 210 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 210
          DM=DMOL5*PDENM(IMOL,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD3=SIGADD3+DM*FAC52
210     CONTINUE
C
C  CHANNEL 4
C  H GAMMA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN PDENI: (DISSOCIATION OF H2+)
C
C  REVISED: USE (PDENM * DENSITY RATIO H2+/H2) NOW, INSTEAD OF PDENI
C
C       DO 215 IION=1,NIONI
C         IF (NCHARI(IION).NE.2) GOTO 215
C         DI=DION5*PDENI(IION,NCELC)
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
C         SIGADD4=SIGADD4+DI*FAC52
C215     CONTINUE
        DO 215 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 215
          DI=DION5*PDENM(IMOL,NCELC)*RATIO2
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD4=SIGADD4+DI*FAC52
215     CONTINUE
C
C  CHANNEL 5
C  H GAMMA SOURCE RATE:  PHOTONS/SEC/CM**3
C  LINEAR IN H- DENSITY (CHARGE EXCHANGE RECOMBINATION)
C
C  REVISED: USE (PDENM * DENSITY RATIO H-/H2) NOW, INSTEAD OF PDENI
C
        DO 220 IMOL=1,NMOLI
          IF (NCHARM(IMOL).NE.2) GOTO 220
          DN=DNML5*PDENM(IMOL,NCELL)*RATIO7
C  RADIATIVE TRANSITION PROB. LEVEL 5-->2 (1/SEC)
C  SIGADD: PHOTONS/SEC/CM**3
          SIGADD5=SIGADD5+DN*FAC52
220     CONTINUE
C
500     CONTINUE
C
C
        SIGADD=SIGADD1+SIGADD2+SIGADD3+SIGADD4+SIGADD5
C
        ADDV(IAD1,NCELC)=ADDV(IAD1,NCELC)+SIGADD1
        ADDV(IAD2,NCELC)=ADDV(IAD2,NCELC)+SIGADD2
        ADDV(IAD3,NCELC)=ADDV(IAD3,NCELC)+SIGADD3
        ADDV(IAD4,NCELC)=ADDV(IAD4,NCELC)+SIGADD4
        ADDV(IAD5,NCELC)=ADDV(IAD5,NCELC)+SIGADD5

        ADDV(IADS,NCELC)=ADDV(IADS,NCELC)+SIGADD
C
C
        POWALF1=POWALF1+SIGADD1*4.577E-19*VOL(NCELL)
        POWALF2=POWALF2+SIGADD2*4.577E-19*VOL(NCELL)
        POWALF3=POWALF3+SIGADD3*4.577E-19*VOL(NCELL)
        POWALF4=POWALF4+SIGADD4*4.577E-19*VOL(NCELL)
        POWALF5=POWALF5+SIGADD5*4.577E-19*VOL(NCELL)

C
        POWALF =POWALF +SIGADD *4.577E-19*VOL(NCELL)
C
1000  CONTINUE

      CALL LEER(2)
      CALL FTCRI(IST,CISTRA)
      IF (IST.GT.0)
     .CALL MASBOX('SUBR. HGAMMA CALLED, FOR STRATUM NO. '//CISTRA)
      IF (IST.EQ.0)
     .CALL MASBOX('SUBR. HGAMMA CALLED, FOR SUM OVER STRATA')
      CALL LEER(1)
      WRITE (iunout,*) ' RADIATED POWER BY HGAMMA:',POWALF
      WRITE (iunout,*) ' COUPL. TO GROUNDSTATE   :',POWALF1
      WRITE (iunout,*) ' COUPLING TO CONTINUUM   :',POWALF2
      WRITE (iunout,*) ' COUPLING TO MOLECULES   :',POWALF3
      WRITE (iunout,*) ' COUPLING TO MOL.IONS    :',POWALF4
      WRITE (iunout,*) ' COUPLING TO NEG.IONS    :',POWALF5
      CALL LEER(2)

c slmod begin
c ... See comments for HALPHA:
      IF (NLTET) THEN
      ELSE
        CALL INTTAL (ADDV,VOL,IAD1,NADV,NSBOX,ADDVI(IAD1,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (ADDV,VOL,IAD2,NADV,NSBOX,ADDVI(IAD2,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (ADDV,VOL,IAD3,NADV,NSBOX,ADDVI(IAD3,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (ADDV,VOL,IAD4,NADV,NSBOX,ADDVI(IAD4,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (ADDV,VOL,IAD5,NADV,NSBOX,ADDVI(IAD5,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)
        CALL INTTAL (ADDV,VOL,IADS,NADV,NSBOX,ADDVI(IADS,IST),
     .               NR1ST,NP2ND,NT3RD,NBMLT)      
      ENDIF
c
c      CALL INTTAL (ADDV,VOL,IAD1,NADV,NSBOX,ADDVI(IAD1,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c      CALL INTTAL (ADDV,VOL,IAD2,NADV,NSBOX,ADDVI(IAD2,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c      CALL INTTAL (ADDV,VOL,IAD3,NADV,NSBOX,ADDVI(IAD3,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c      CALL INTTAL (ADDV,VOL,IAD4,NADV,NSBOX,ADDVI(IAD4,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c      CALL INTTAL (ADDV,VOL,IAD5,NADV,NSBOX,ADDVI(IAD5,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c      CALL INTTAL (ADDV,VOL,IADS,NADV,NSBOX,ADDVI(IADS,IST),
c     .             NR1ST,NP2ND,NT3RD,NBMLT)
c slmod end
C
      IF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.IST.EQ.0) THEN
        IESTR=IST
        CALL WRSTRT(IST,NSTRAI,NESTM1,NESTM2,NADSPC,
     .              ESTIMV,ESTIMS,ESTIML,
     .              NSDVI1,SDVI1,NSDVI2,SDVI2,
     .              NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .              NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .              NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .              NSIGI_SPC,TRCFLE)
C
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU, IUNOUT)
        WRITE (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE)   WRITE (iunout,*) 'WRITE 11  IRC= ',IRC
      ENDIF
C
      RETURN
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. HGAMMA '
      WRITE (iunout,*) 'NO STORAGE AVAILBALE ON ADDITIONAL TALLY ADDV '
      WRITE (iunout,*) 'STORAGE REQUESTED FOR IADV= ',
     .             IAD1,IAD2,IAD3,IAD4,
     .             IAD5,IADS
      WRITE (iunout,*) 'CHECK INPUT BLOCK 10A '
      CALL EXIT_OWN(1)
      END







C ===== SOURCE: movie.f
C
      SUBROUTINE MOVIE
C
C  MODIFICATIONS TO INPUT, IN ORDER TO PERFORM MANY TIME STEPS
C  AND TO PRODUCE A MOVIE WITH TEST PARTICLE HISTORIES
C
      USE PRECISION
      USE PARMMOD
      USE CPLOT
      USE CTRCEI
      USE COMPRT
      USE COMNNL
      USE COMSOU

      IMPLICIT NONE

      INTEGER :: ISTR

C  TRY TO MAINTAIN A SWARM OF NPRNLI SOURCE PARTICLES AT EACH TIMESTEP,
C  COLD START FROM CENSUS

      PLHST=.TRUE.
      I1TRC=1
      I2TRC=NPRNLI
C
      NVOLPR=0
      NVOLPL=0
C
      DO ISTR=1,NSTRAI
        NINITL(ISTR)=0
      ENDDO
C
      RETURN
      END
C ===== SOURCE: outeir.f
C  11.01.05:   text re. switching off generation limit corrected
C              electron particle balance wtote instead of wtotp
C
      SUBROUTINE OUTEIR(INDOUT)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CADGEO
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CSPEZ
      USE CTRCEI
      USE CGEOM
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE COMPRT
      USE COMNNL
      USE COMSOU
      USE CTEXT
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE
C
      INTEGER, INTENT(IN) :: INDOUT
      REAL(DP) :: VECTOR(NRAD)
      REAL(DP) :: TOTA(0:NATM),DIFA(0:NATM,0:NSTRA),
     .            DIFRA(0:NATM,0:NSTRA)
      REAL(DP) :: TOTM(0:NMOL),DIFM(0:NMOL,0:NSTRA),
     .            DIFRM(0:NMOL,0:NSTRA)
      REAL(DP) :: TOTI(0:NION),DIFI(0:NION,0:NSTRA),
     .            DIFRI(0:NION,0:NSTRA)
      REAL(DP) :: TOTPH(0:NPHOT),DIFPH(0:NPHOT,0:NSTRA),
     .            DIFRPH(0:NPHOT,0:NSTRA)
      REAL(DP) :: DIFR, TOTT, PGAINP, PLOSSP, DIFT, EGAINE, SPA, ELOSSE,
     .            EGAINP, ELOSSP, PGAINE, PLOSSE, SMEAN, OUTAUI, TALAV,
     .            TALTOT, DIF
      INTEGER :: IALS, IALV, J, JJ, IS, NFTI, NFTE, I0, K, NF, N, IALG,
     .           ITAL, IPRV, KMAX, KK, ILAST, ICOUNT, I, IINDEX, ISPC, 
     .           IT
      INTEGER :: IADTYP(0:4)
      LOGICAL :: LCOVN(NCV)
C
      CHARACTER(72) :: TXTTL
      CHARACTER(24) :: TXTSP, TXTUN
      CHARACTER(6) :: CISTRA
      CHARACTER(10) :: TEXTYP(0:4)
C
C
C ISTRA IS THE STRATUM NUMBER. ISTRA=0 STANDS FOR: SUM OVER STRATA
      ISTRA=INDOUT
C
C MAKE SURE TO PRINT COVARIANCE TALLIES ONLY ONCE
      DO N=1,NSIGCI
        LCOVN(N)=.FALSE.
      ENDDO
C
      IF (XMCP(ISTRA).LT.1.) THEN
        WRITE (IUNOUT,*) 'OUTEIR CALLED, ISTRA, XMCP: ',
     .                    ISTRA,XMCP(ISTRA)
        WRITE (IUNOUT,*) 'PRINTOUT ABANDONNED FOR THIS STRATUM '
        CALL LEER(1)
        GOTO 1000
      ENDIF
C
C
      IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
      ELSEIF (NFILEN.EQ.1.OR.NFILEN.EQ.2) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSEIF ((NFILEN.EQ.6.OR.NFILEN.EQ.7).AND.ISTRA.EQ.0) THEN
        IESTR=ISTRA
        CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .             ESTIMV,ESTIMS,ESTIML,
     .             NSDVI1,SDVI1,NSDVI2,SDVI2,
     .             NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .             NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .             NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .             NSIGI_SPC,TRCFLE)
        IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
          CALL SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .               NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
        ENDIF
      ELSE
        WRITE (iunout,*) 'ERROR IN OUTEIR: DATA FOR STRATUM ISTRA= ',
     .                    ISTRA
        WRITE (iunout,*) 'ARE NOT AVAILABLE. PRINTOUT ABANDONNED'
        RETURN
      ENDIF
C
C
      CALL PAGE
      IF (ISTRA.NE.0) THEN
        CALL FTCRI(ISTRA,CISTRA)
        CALL MASBOX ('THIS IS STRATUM NUMBER ISTRA='//CISTRA)
        WRITE (iunout,*) TXTSOU(ISTRA)
        CALL LEER(1)
        CALL MASAGE ('TYPE OF PRIMARY SOURCE:                      ')
        CALL MASL4('NLPNT,NLLNE,NLSRF,NLVOL         ',
     .              NLPNT(ISTRA),NLLNE(ISTRA),NLSRF(ISTRA),NLVOL(ISTRA))
        CALL MASL5('NLPHOT,NLATM,NLMOL,NLION,NLPLS          ',
     .              NLPHOT(ISTRA),NLATM(ISTRA),NLMOL(ISTRA),
     .              NLION (ISTRA),NLPLS(ISTRA))
        CALL MASAGE ('SPECIES INDEX OF SOURCE PARTICLES             ')
        CALL MASJ1 ('NSPEZ=  ',NSPEZ(ISTRA))
      ELSE
        CALL MASBOX ('THIS IS THE SUM OVER THE STRATA')
        CALL LEER(1)
      ENDIF
      CALL MASAGE ('NUMBER OF MONTE-CARLO HISTORIES               ')
      CALL MASR1('NHIST=  ',XMCP(ISTRA))
      CALL LEER(3)
C
C  PRINT THOSE VOLUME AVERAGED TALLIES, WHICH HAVE BEEN SELECTED
C  AND WHICH ARE DIFFERENT FROM ZERO
      IALG=0
      DO 100 IPRV=1,NVOLPR
        ITAL=NPRTLV(IPRV)
        IF ((ITAL > 0) .AND. (ITAL <= NTALV)) THEN
          IF (.NOT.LIVTALV(ITAL)) THEN
            CALL LEER(2)
            WRITE (iunout,*) TXTTAL(1,ITAL)
            CALL MASAGE
     .        ('TALLY SWITCHED OFF => NO OUTPUT AVAILABLE      ')
            CALL LEER(2)
            CYCLE
          END IF
        END IF
C
C  USER SUPPLIED POST PROCESSED TALLY, ITAL=0
C
        IF (ITAL.EQ.0) THEN
C   CALL TO TALUSR: A POST PROCESSED USER SUPPLIED TALLY
          ICOUNT=1
120       CALL TALUSR(ICOUNT,VECTOR,TALTOT,TALAV,
     .                TXTTL,TXTSP,TXTUN,ILAST,*121)
          WRITE (iunout,*) 'USER SUPPLIED POST PROCESSED TALLY NO. ',
     .                ICOUNT
          CALL PRTTAL(TXTTL,TXTSP,TXTUN,
     .                VECTOR,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                NFLAGV(IPRV),NTLVFL(IPRV))
          CALL LEER(2)
          CALL MASAGE
     .       ('TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") ')
          CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)       
          CALL LEER(3)
121       ICOUNT=ICOUNT+1
          IF (ICOUNT.LE.ILAST) GOTO 120
        ELSEIF (ITAL.GT.0.AND.ITAL.NE.NTALR) THEN
          I0=0
          IF (NFRSTI(ITAL).GT.1) I0=1
          NFTI=1
          NFTE=NFSTVI(ITAL)
          IF (NSPEZV(IPRV,1).GT.0) THEN
            NFTI=NSPEZV(IPRV,1)
            NFTE=MAX(NFTI,NSPEZV(IPRV,2))
          ENDIF
          NF=NFIRST(ITAL)
          DO 119 K=NFTI,NFTE
            IF (K.GT.NFSTVI(ITAL)) GOTO 106
            CALL FETCH_OUTAU (OUTAUI,ITAL,K,ISTRA,IUNOUT)
            IF (OUTAUI.EQ.0.D0) GOTO 118
C
            DO 110 I=1,NSBOX_TAL
              VECTOR(I)=ESTIMV(NADDV(ITAL)+K,I)
110         CONTINUE
            TALTOT=OUTAUI
            TALAV=TALTOT/VOLTOT
C
            CALL PRTTAL(TXTTAL(K,ITAL),TXTSPC(K,ITAL),TXTUNT(K,ITAL),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                  NFLAGV(IPRV),NTLVFL(IPRV))
            CALL LEER(2)
            CALL MASAGE
     .      ('TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") ')
            CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
            CALL LEER(3)
C
C  CHECK IF STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY "ITAL"
            DO 101 N=1,NSIGVI
              IF (IIH(N).NE.ITAL) GOTO 101
              IF (IGH(N).NE.K.AND.IGH(N).NE.0) GOTO 101
              GOTO 111
101         CONTINUE
            GOTO 106
111         CONTINUE
            DO 112 I=1,NSBOX_TAL
              VECTOR(I)=SIGMA(N,I)
112         CONTINUE
            SMEAN=SGMS(N)
C
            CALL MASAGE
     .      ('RELATIVE STANDARD DEVIATION                   ')
            IF (IGH(N).GT.0) TXTSP=TXTSPC(K,ITAL)
            IF (IGH(N).EQ.0) TXTSP='TOTAL                   '
            CALL PRTTAL
     .           (TXTTAL(K,ITAL),TXTSP,'%                       ',
     .            VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .            NFLAGV(IPRV),NTLVFL(IPRV))
            CALL LEER(2)
            CALL MASAGE
     .      ('STANDARD DEVIATION OF MEAN VALUE (%)          ')
            CALL MASR1 ('MEAN    ',SMEAN)
            CALL LEER(3)
C
106         CONTINUE
C  CHECK IF BGK-STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY "ITAL"
            IF (NSIGI_BGK.GT.0) THEN
              IF (ITAL.EQ.NTALB) THEN
                KMAX=NBGVI_STAT
                KK=K
              ELSEIF (ITAL.EQ.1) THEN
                KMAX=NATMI
                KK=NBGVI+K
              ELSEIF (ITAL.EQ.5) THEN
                KMAX=NATMI
                KK=NBGVI+NATMI+K
              ELSEIF (ITAL.EQ.2) THEN
                KMAX=NMOLI
                KK=NBGVI+2*NATMI+K
              ELSEIF (ITAL.EQ.6) THEN
                KMAX=NMOLI
                KK=NBGVI+2*NATMI+NMOLI+K
              ELSE
                KMAX=0
              ENDIF
              IF (K.GT.KMAX) GOTO 119
              DO I=1,NSBOX_TAL
                VECTOR(I)=SIGMA_BGK(KK,I)
              ENDDO
              SMEAN=SGMS_BGK(KK)
C
              CALL MASAGE
     .        ('RELATIVE STANDARD DEVIATION (BGK)              ')
              TXTSP=TXTSPC(K,ITAL)
              CALL PRTTAL(TXTTAL(K,ITAL),TXTSP,'%                     ',
     .                    VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
              CALL LEER(2)
              CALL MASAGE
     .        ('STANDARD DEVIATION OF MEAN VALUE (%)          ')
              CALL MASR1 ('MEAN    ',SMEAN)
              CALL LEER(3)
              GOTO 119
            ENDIF
C
C  CHECK IF COP-STANDARD DEVIATION IS AVAILABLE FOR THIS TALLY "ITAL"
            IF (NSIGI_COP.GT.0) THEN
              IF (ITAL.EQ.NTALM) THEN
                KMAX=NCPVI_STAT
              ELSE
                KMAX=0
              ENDIF
              IF (K.GT.KMAX) GOTO 119
              DO I=1,NSBOX_TAL
                VECTOR(I)=SIGMA_COP(K,I)
              ENDDO
              SMEAN=SGMS_COP(K)
C
              CALL MASAGE
     .        ('RELATIVE STANDARD DEVIATION (COPV)             ')
              TXTSP=TXTSPC(K,ITAL)
              CALL PRTTAL(TXTTAL(K,ITAL),TXTSP,'%                     ',
     .                    VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
              CALL LEER(2)
              CALL MASAGE
     .        ('STANDARD DEVIATION OF MEAN VALUE (%)          ')
              CALL MASR1 ('MEAN    ',SMEAN)
              CALL LEER(3)
              GOTO 119
            ENDIF
C
C  CHECK IF COVARIANCE IS AVAILABLE
            DO 103 N=1,NSIGCI
              IF (LCOVN(N)) GOTO 103
              IF  (IIHC(1,N).NE.ITAL.AND.IIHC(2,N).NE.ITAL) GOTO 103
              IF  (IGHC(1,N).NE.K.AND.IGHC(1,N).NE.0.AND.
     .             IGHC(2,N).NE.K.AND.IGHC(2,N).NE.0) GOTO 103
              LCOVN(N)=.TRUE.
              GOTO 113
103         CONTINUE
            GOTO 119
113         CONTINUE
C
            DO 114 I=1,NSBOX_TAL
              VECTOR(I)=SIGMAC(0,N,I)
114         CONTINUE
            SMEAN=SGMCS(0,N)
C
            CALL MASAGE
     .      ('COVARIANCE                                    ')
            IF (IGHC(1,N).GT.0) TXTSP=TXTSPC(IGHC(1,N),IIHC(1,N))
            IF (IGHC(1,N).EQ.0) TXTSP='TOTAL                   '
            WRITE (iunout,*) 'BETWEEN ESTIMATOR: '
            ISPZ=MAX(1,IGHC(1,N))
            CALL PRTTAL(TXTTAL(ISPZ,IIHC(1,N)),TXTSP,
     .                  TXTUNT(ISPZ,IIHC(1,N)),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                 -1,0)
            WRITE (iunout,*) 'AND ESTIMATOR: '
            ISPZ=MAX(1,IGHC(2,N))
            CALL PRTTAL(TXTTAL(ISPZ,IIHC(2,N)),TXTSP,
     .                  TXTUNT(ISPZ,IIHC(2,N)),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                  NFLAGV(IPRV),0)
            CALL LEER(2)
            CALL MASAGE
     .      ('COVARIANCE OF MEAN VALUES                      ')
            CALL MASR1 ('MEAN    ',SMEAN)
            CALL LEER(3)
C
            DO 115 I=1,NSBOX_TAL
              VECTOR(I)=SIGMAC(1,N,I)
115         CONTINUE
            SMEAN=SGMCS(1,N)
C
            CALL MASAGE
     .      ('STANDARD DEVIATION FOR 1ST TALLY               ')
            IF (IGHC(1,N).GT.0) TXTSP=TXTSPC(IGHC(1,N),IIHC(1,N))
            IF (IGHC(1,N).EQ.0) TXTSP='TOTAL                   '
            ISPZ=MAX(1,IGHC(1,N))
            CALL PRTTAL(TXTTAL(ISPZ,IIHC(1,N)),TXTSP,
     .                  TXTUNT(ISPZ,IIHC(1,N)),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                  NFLAGV(IPRV),0)
            CALL LEER(2)
            CALL MASAGE
     .      ('STANDARD DEVIATION OF MEAN VALUE              ')
            CALL MASR1 ('MEAN    ',SMEAN)
            CALL LEER(3)
C
            DO 116 I=1,NSBOX_TAL
              VECTOR(I)=SIGMAC(2,N,I)
116         CONTINUE
            SMEAN=SGMCS(2,N)
C
            CALL MASAGE
     .      ('STANDARD DEVIATION FOR 2ND TALLY              ')
            IF (IGHC(2,N).GT.0) TXTSP=TXTSPC(IGHC(2,N),IIHC(2,N))
            IF (IGHC(2,N).EQ.0) TXTSP='TOTAL                   '
            ISPZ=MAX(1,IGHC(2,N))
            CALL PRTTAL(TXTTAL(ISPZ,IIHC(2,N)),TXTSP,
     .                  TXTUNT(ISPZ,IIHC(2,N)),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                  NFLAGV(IPRV),0)
            CALL LEER(2)
            CALL MASAGE
     .      ('STANDARD DEVIATION OF MEAN VALUE              ')
            CALL MASR1 ('MEAN    ',SMEAN)
            CALL LEER(3)
            GOTO 119
C
118         CONTINUE
            CALL PRTTAL(TXTTAL(K,ITAL),TXTSPC(K,ITAL),TXTUNT(K,ITAL),
     .                VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,-1,0)
            CALL MASAGE('IDENTICAL ZERO, NOT PRINTED                  ')
            CALL LEER(2)
119       CONTINUE
C
C  ALGEBRAIC TALLY, ITAL=NTALR
C
        ELSEIF (ITAL.EQ.NTALR) THEN
C
C   REDO ALGEBRAIC EXPRESSION IN TALLIES, IN CASE NFILEN=2
          IF ((NFILEN.EQ.2..OR.NFILEN.EQ.7).AND.IALG.EQ.0) THEN
            IALG=1
            CALL ALGTAL
C
            DO 130 IALV=1,NALVI
              CALL INTTAL (ALGV,VOLTAL,IALV,NALV,NSBOX_TAL,
     .                     ALGVI(IALV,ISTRA),
     .                     NR1TAL,NP2TAL,NT3TAL,NBMLT)
130         CONTINUE
C
            DO 132 IALS=1,NALSI
              ALGSI(IALS,ISTRA)=0.
              DO 131 J=1,NLIMPS
                ALGSI(IALS,ISTRA)=ALGSI(IALS,ISTRA)+ALGS(IALS,J)
131           CONTINUE
132         CONTINUE
          ENDIF
C
          I0=0
          IF (NFRSTI(ITAL).GT.1) I0=1
          NFTI=1
          NFTE=NFSTVI(ITAL)
          IF (NSPEZV(IPRV,1).GT.0) THEN
            NFTI=NSPEZV(IPRV,1)
            NFTE=MAX(NFTI,NSPEZV(IPRV,2))
          ENDIF
          NF=NFIRST(ITAL)
          DO 159 K=NFTI,NFTE
            CALL FETCH_OUTAU (OUTAUI,ITAL,K,ISTRA,IUNOUT)
            IF (OUTAUI.EQ.0.D0) GOTO 158
C
            DO 155 I=1,NSBOX_TAL
              IINDEX=NADDV(ITAL)*NRAD+(I-1)*NF+K
              VECTOR(I)=ESTIMV(NADDV(ITAL)+K,I)
155         CONTINUE
            TALTOT=OUTAUI
            TALAV=TALTOT/VOLTOT
            CALL PRTTAL(TXTTAL(K,ITAL),TXTSPC(K,ITAL),TXTUNT(K,ITAL),
     .                  VECTOR,NR1TAL,NP2TAL,NT3TAL,NBMLT,NSBOX_TAL,
     .                  NFLAGV(IPRV),NTLVFL(IPRV))
            CALL LEER(2)
            CALL MASAGE
     .         ('TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") ')
                CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
            CALL LEER(3)
            GOTO 159
C
158         CONTINUE
            CALL PRTTAL(TXTTAL(K,ITAL),TXTSPC(K,ITAL),TXTUNT(K,ITAL),
     .                  VECTOR,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,-1,0)
            CALL MASAGE('IDENTICAL ZERO, NOT PRINTED                  ')
            CALL LEER(2)
159       CONTINUE
C
        ENDIF
C
100   CONTINUE

C  SPECTRA IN SELECTED CELLS

      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO ISPC=1,NADSPC
        IF (ESTIML(ISPC)%PSPC%ISRFCLL /= 0) THEN
          CALL LEER (1)
          WRITE (iunout,'(A,I6)') ' SPECTRUM CALCULATED FOR CELL ',
     .      ESTIML(ISPC)%PSPC%ISPCSRF
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))') 
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IT = ESTIML(ISPC)%PSPC%ISPCTYP
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
          WRITE (iunout,'(A20,A8)') ' TYPE OF PARTICLE : ',
     .                       TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
          IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
            WRITE (iunout,'(A10,10X,A16)') ' SPECIES :',
     .                                     'SUM OVER SPECIES'
          ELSE
            WRITE (iunout,'(A10,10X,A8)') ' SPECIES :',
     .             TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .                   ESTIML(ISPC)%PSPC%IPRSP)
          END IF
          WRITE (iunout,'(A22,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .           ESTIML(ISPC)%PSPC%SPCINT
          IF (NSIGI_SPC > 0)
     .      WRITE (iunout,'(A22,ES12.4)') ' STANDARD DEVIATION   ',
     .           ESTIML(ISPC)%PSPC%SGMS
        END IF
      END DO

c slmod begin      
      IF (NSPCPR > 0) CALL OUTSPEC(ISTRA)
c
c      IF (NSPCPR > 0) CALL OUTSPEC
c slmod end
C
C   OUTPUT OF VOLUME AVERAGED TALLIES FINISHED
C
      IF (TRCBLPH.OR.TRCBLA.OR.TRCBLM.OR.TRCBLI.OR.
     .    TRCBLP .OR.TRCBLE) THEN
        CALL LEER (3)
        CALL HEADNG('                               ',31)
        CALL HEADNG('= GLOBAL BALANCES (AMP/WATT) =',31)
        CALL LEER (3)
      ENDIF
C
C
C  TEST PARTICLE FLUX BALANCE - AND ENERGY FLUX BALANCE
C
      CALL MASAGE ('TEST PARTICLE INFLUX FROM  SOURCE (AMP)       ')
      CALL MASAGE ('INFLUX/(1.602*E-19) IS THE "ATOMIC" FLUX, PART./S')
      CALL MASR1 ('INFLUX= ',FLUXT(ISTRA))
      CALL LEER(2)
C
C  NEUTRAL ATOMS PARTICLE BALANCE
C
      DIFA(0,ISTRA)=0.
      TOTA(0)=0.
      IF (.NOT.TRCBLA) GOTO 410
      IF (.NOT.LOGATM(0,ISTRA)) GOTO 410
      CALL HEADNG('PARTICLE FLUX BALANCE (AMP), NEUTRAL ATOMS',43)
      CALL LEER(1)
C   ATOMS FROM PRIMARY SOURCE
      IF (ANY(WTOTA(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM PRIMARY SOURCE              ')
        CALL MASYR1 ('WTOTA  = ',
     .                WTOTA,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',WTOTA(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPPAT) THEN
        IF (LMSPPAT) THEN
        CALL MASAGE ('ATOMS FROM RECOMBINING BULK IONS               ')
        CALL MASAGE ('PPAT   SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PPATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS FROM RECOMBINING BULK IONS               ')
        CALL MASYR1 ('PPATI  = ',
     .                PPATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPATI(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPAAT) THEN
        IF (LMSPAAT) THEN
        CALL MASAGE ('ATOMS BORN BY ATOM - PLASMA INTERACTIONS       ')
        CALL MASAGE ('PAAT   SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PAATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS BORN BY ATOM - PLASMA INTERACTIONS       ')
        CALL MASYR1 ('PAATI  = ',
     .                PAATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PAATI(0,ISTRA))
      ENDIF
      IF (.NOT.LPMAT) THEN
        IF (LMSPMAT) THEN
        CALL MASAGE ('ATOMS BORN BY MOLECULE - PLASMA INTERACTIONS   ')
        CALL MASAGE ('PMAT   SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PMATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS BORN BY MOLECULE - PLASMA INTERACTIONS   ')
        CALL MASYR1 ('PMATI  = ',
     .                PMATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMATI(0,ISTRA))
      ENDIF
      IF (.NOT.LPIAT) THEN
        IF (LMSPIAT) THEN
        CALL MASAGE ('ATOMS BORN BY TEST ION - PLASMA INTERACTIONS   ')
        CALL MASAGE ('PIAT   SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PIATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS BORN BY TEST ION - PLASMA INTERACTIONS   ')
        CALL MASYR1 ('PIATI  = ',
     .                PIATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PIATI(0,ISTRA))
      ENDIF
      IF (.NOT.LPPHAT) THEN
        IF (LMSPPHAT) THEN
        CALL MASAGE ('ATOMS BORN BY PHOTON - PLASMA INTERACTIONS     ')
        CALL MASAGE ('PPHAT  SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PPHATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS BORN BY PHOTON - PLASMA INTERACTIONS     ')
        CALL MASYR1 ('PPHATI = ',
     .                PPHATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPHATI(0,ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LPGENA) THEN
        IF (LMSPGENA) THEN
        CALL MASAGE ('ATOMS LOST DUE TO GENERATION-OR FLUID LIMIT'    )
        CALL MASAGE ('PGENA  SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PGENAI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMS LOST DUE TO GENERATION-OR FLUID LIMIT')
        CALL MASYR1 ('PGENAI = ',
     .                PGENAI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PGENAI(0,ISTRA))
      ENDIF
C   ESCAPING FLUX
      IF (.NOT.LPOTAT) THEN
        IF (LMSPOTAT) THEN
        CALL MASAGE ('ATOMIC EFFLUX ONTO THE SURFACES                ')
        CALL MASAGE ('POTAT  SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(POTATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC EFFLUX ONTO THE SURFACES                ')
        CALL MASYR1 ('POTATI = ',
     .                POTATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',POTATI(0,ISTRA))
      ENDIF
C   REFLECTED FLUX
      IF (.NOT.LPRFAAT) THEN
        IF (LMSPRFAAT) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: ATOMS       ')
        CALL MASAGE ('PRFAAT SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PRFAAI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: ATOMS       ')
        CALL MASYR1 ('PRFAAI = ',
     .                PRFAAI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFAAI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFMAT) THEN
        IF (LMSPRFMAT) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: MOLECULES   ')
        CALL MASAGE ('PRFMAT SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PRFMAI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: MOLECULES   ')
        CALL MASYR1 ('PRFMAI = ',
     .                PRFMAI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFMAI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFIAT) THEN
        IF (LMSPRFIAT) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: TEST IONS   ')
        CALL MASAGE ('PRFIAT SWITCHED OFF => MISSING IN BALANCE      ')
        ENDIF
      ELSE IF (ANY(PRFIAI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: TEST IONS   ')
        CALL MASYR1 ('PRFIAI = ',
     .                PRFIAI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFIAI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFPHAT) THEN
        IF (LMSPRFPHAT) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: PHOTONS')
        CALL MASAGE ('PRFPHAT SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFPHAI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('ATOMIC INFLUX FROM SURFACES, ORIG: PHOTONS')
        CALL MASYR1 ('PRFPHAI= ',
     .                PRFPHAI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFPHAI(0,ISTRA))
      ENDIF
      IF (.NOT.LSPTAT) THEN
        IF (LMSSPTAT) THEN
        CALL MASAGE ('FLUX SPUTTERED BY ATOMS:                       ')
        CALL MASAGE ('SPTAT   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(SPTATI(1:NATM,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('FLUX SPUTTERED BY ATOMS:                       ')
        CALL MASYR1 ('SPTATI = ',
     .                SPTATI,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',SPTATI(0,ISTRA))
      ENDIF
      DO 405 IATM=0,NATMI
        DIFA(IATM,ISTRA)=
     +       POTATI(IATM,ISTRA)+PPATI(IATM,ISTRA)+
     +       PGENAI(IATM,ISTRA)+
     +       WTOTA(IATM,ISTRA)+
     +       PAATI (IATM,ISTRA)+PRFAAI (IATM,ISTRA)+
     +       PMATI (IATM,ISTRA)+PRFMAI (IATM,ISTRA)+
     +       PIATI (IATM,ISTRA)+PRFIAI (IATM,ISTRA)+
     +       PPHATI(IATM,ISTRA)+PRFPHAI(IATM,ISTRA)
        TOTA(IATM)=
     +  ABS(POTATI(IATM,ISTRA))+ABS(PAATI(IATM,ISTRA))+
     +       ABS(PPATI(IATM,ISTRA))+ABS(PMATI(IATM,ISTRA))+
     +       ABS(PIATI(IATM,ISTRA))+ABS(WTOTA(IATM,ISTRA))+
     +       ABS(PRFAAI(IATM,ISTRA))+ABS(PRFMAI(IATM,ISTRA))+
     +       ABS(PRFIAI(IATM,ISTRA))+ABS(PGENAI(IATM,ISTRA))+
     +       ABS(PPHATI(IATM,ISTRA))+ABS(PRFPHAI(IATM,ISTRA))
        TOTA(IATM)=TOTA(IATM)+EPS60
        DIF=SIGN (1._DP,DIFA(IATM,ISTRA))*
     .      MAX(0._DP,ABS(DIFA(IATM,ISTRA))/TOTA(IATM)-EPS10)
        DIFRA(IATM,ISTRA)=DIF*100.
        DIFA(IATM,ISTRA)=DIF*TOTA(IATM)
405   CONTINUE
      CALL MASAGE ('ABSOLUTE ERRORS IN PARTICLE BALANCE            ')
      CALL MASYR1 ('DIFA =   ',
     .              DIFA,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFA(0,ISTRA))
      CALL MASAGE ('RELATIVE ERRORS IN PARTICLE BALANCE (%)        ')
      CALL MASYR1 ('DIFRA =  ',
     .              DIFRA,LOGATM,ISTRA,0,NATM,0,NSTRA,TEXTS(NSPH+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFRA(0,ISTRA))
      CALL LEER(2)
C
C  NEUTRAL MOLECULES PARTICLE BALANCE
C
410   CONTINUE
      DIFM(0,ISTRA)=0.
      TOTM(0)=0.
      IF (.NOT.TRCBLM) GOTO 420
      IF (.NOT.LOGMOL(0,ISTRA)) GOTO 420
      CALL LEER (2)
      CALL HEADNG('PARTICLE FLUX BALANCE (AMP), NEUTRAL MOLECULES',47)
      CALL LEER(1)
C  MOLECULES FROM PRIMARY SOURCE
      IF (ANY(WTOTM(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM PRIMARY SOURCE          ')
        CALL MASYR1 ('WTOTM  = ',
     .                WTOTM,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',WTOTM(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPPML) THEN
        IF (LMSPPML) THEN
        CALL MASAGE ('MOLECULES FROM RECOMBINING BULK IONS          ')
        CALL MASAGE ('PPML    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES FROM RECOMBINING BULK IONS          ')
        CALL MASYR1 ('PPMLI  = ',
     .                PPMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPMLI(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPAML) THEN
        IF (LMSPAML) THEN
        CALL MASAGE ('MOLECULES BORN BY ATOM - PLASMA INTERACTIONS    ')
        CALL MASAGE ('PAML    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PAMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES BORN BY ATOM - PLASMA INTERACTIONS    ')
        CALL MASYR1 ('PAMLI  = ',
     .                PAMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMMLI(0,ISTRA))
      ENDIF
      IF (.NOT.LPMML) THEN
        IF (LMSPMML) THEN
        CALL MASAGE ('MOLECULES BORN BY MOLECULE - PLASMA INTERACTIONS')
        CALL MASAGE ('PMML    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PMMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES BORN BY MOLECULE - PLASMA INTERACTIONS')
        CALL MASYR1 ('PMMLI  = ',
     .                PMMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMMLI(0,ISTRA))
      ENDIF
      IF (.NOT.LPIML) THEN
        IF (LMSPIML) THEN
        CALL MASAGE ('MOLECULES BORN BY TEST ION - PLASMA INTERACTIONS')
        CALL MASAGE ('PIML    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PIMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES BORN BY TEST ION - PLASMA INTERACTIONS')
        CALL MASYR1 ('PIMLI  = ',
     .                PIMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PIMLI(0,ISTRA))
      ENDIF
      IF (.NOT.LPPHML) THEN
        IF (LMSPPHML) THEN
        CALL MASAGE ('MOLECULES BORN BY PHOTON - PLASMA INTERACTIONS')
        CALL MASAGE ('PPHML   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPHMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES BORN BY PHOTON - PLASMA INTERACTIONS')
        CALL MASYR1 ('PPHMLI = ',
     .                PPHMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPHMLI(0,ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LPGENM) THEN
        IF (LMSPGENM) THEN
        CALL MASAGE ('MOLECULES LOST DUE TO GENERATION-OR FLUID LIMIT')
        CALL MASAGE ('PGENM   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PGENMI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULES LOST DUE TO GENERATION-OR FLUID LIMIT')
        CALL MASYR1 ('PGENMI = ',
     .                PGENMI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PGENMI(0,ISTRA))
      ENDIF
C  ESCAPING MOLECULES
      IF (.NOT.LPOTML) THEN
        IF (LMSPOTML) THEN
        CALL MASAGE ('MOLECULAR EFFLUX ONTO THE SURFACES             ')
        CALL MASAGE ('POTML   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(POTMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR EFFLUX ONTO THE SURFACES             ')
        CALL MASYR1 ('POTMLI = ',
     .                POTMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',POTMLI(0,ISTRA))
      ENDIF
C  REFLECTED MOLECULES
      IF (.NOT.LPRFAML) THEN
        IF (LMSPRFAML) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: ATOMS    ')
        CALL MASAGE ('PRFAML  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFAMI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: ATOMS    ')
        CALL MASYR1 ('PRFAMI = ',
     .                PRFAMI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFAMI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFMML) THEN
        IF (LMSPRFMML) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: MOLECULES')
        CALL MASAGE ('PRFMML  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFMMI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: MOLECULES')
        CALL MASYR1 ('PRFMMI = ',
     .                PRFMMI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFMMI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFIML) THEN
        IF (LMSPRFIML) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: TEST IONS')
        CALL MASAGE ('PRFIML  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFIMI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: TEST IONS')
        CALL MASYR1 ('PRFIMI = ',
     .                PRFIMI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFIMI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFPHML) THEN
        IF (LMSPRFPHML) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: PHOTONS')
        CALL MASAGE ('PRFPHML SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFPHMI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('MOLECULAR INFLUX FROM SURFACES, ORIG: PHOTONS')
        CALL MASYR1 ('PRFPHMI= ',
     .                PRFPHMI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFPHMI(0,ISTRA))
      ENDIF
      IF (.NOT.LSPTML) THEN
        IF (LMSSPTML) THEN
        CALL MASAGE ('FLUX SPUTTERED BY MOLECULES:                   ')
        CALL MASAGE ('SPTML   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(SPTMLI(1:NMOL,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('FLUX SPUTTERED BY MOLECULES:                   ')
        CALL MASYR1 ('SPTMLI = ',
     .                SPTMLI,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',SPTMLI(0,ISTRA))
      ENDIF
      DO 415 IMOL=0,NMOLI
        DIFM(IMOL,ISTRA)=
     +     POTMLI(IMOL,ISTRA)+PMMLI(IMOL,ISTRA)+PPMLI(IMOL,ISTRA)+
     +     PIMLI(IMOL,ISTRA)+PGENMI(IMOL,ISTRA)+
     +     WTOTM(IMOL,ISTRA)+PRFAMI(IMOL,ISTRA)+PRFMMI(IMOL,ISTRA)+
     +     PRFIMI(IMOL,ISTRA)+PAMLI(IMOL,ISTRA)+
     +     PPHMLI(IMOL,ISTRA)+PRFPHMI(IMOL,ISTRA)
        TOTM(IMOL)=
     +     ABS(POTMLI(IMOL,ISTRA))+ABS(PMMLI(IMOL,ISTRA))+
     +     ABS(PPMLI(IMOL,ISTRA))+ABS(WTOTM(IMOL,ISTRA))+
     +     ABS(PIMLI(IMOL,ISTRA))+ABS(PGENMI(IMOL,ISTRA))+
     +     ABS(PRFAMI(IMOL,ISTRA))+ABS(PRFMMI(IMOL,ISTRA))+
     +     ABS(PAMLI(IMOL,ISTRA))+
     +     ABS(PRFIMI(IMOL,ISTRA))+
     +     ABS(PPHMLI(IMOL,ISTRA))+ABS(PRFPHMI(IMOL,ISTRA))
        TOTM(IMOL)=TOTM(IMOL)+EPS60
        DIF=SIGN (1._DP,DIFM(IMOL,ISTRA))*
     .      MAX(0._DP,ABS(DIFM(IMOL,ISTRA))/TOTM(IMOL)-EPS10)
        DIFRM(IMOL,ISTRA)=DIF*100.
        DIFM(IMOL,ISTRA)=DIF*TOTM(IMOL)
415   CONTINUE
      CALL MASAGE ('ABSOLUTE ERRORS IN PARTICLE BALANCE            ')
      CALL MASYR1 ('DIFM =   ',
     .              DIFM,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFM(0,ISTRA))
      CALL MASAGE ('RELATIVE ERRORS IN PARTICLE BALANCE (%)        ')
      CALL MASYR1 ('DIFRM =  ',
     .              DIFRM,LOGMOL,ISTRA,0,NMOL,0,NSTRA,TEXTS(NSPA+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFRM(0,ISTRA))
      CALL LEER(2)
C
C  TEST IONS PARTICLE BALANCE
C
420   CONTINUE
      DIFI(0,ISTRA)=0.
      TOTI(0)=0.
      IF (.NOT.TRCBLI) GOTO 430
      IF (.NOT.LOGION(0,ISTRA)) GOTO 430
      CALL LEER (2)
      CALL HEADNG('PARTICLE FLUX BALANCE (AMP), TEST IONS',39)
      CALL LEER(1)
C  TEST IONS FROM PRIMARY SOURCE
      IF (ANY(WTOTI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS INFLUX FROM PRIMARY SOURCE          ')
        CALL MASYR1 ('IOFLUX = ',
     .                WTOTI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',WTOTI(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPPIO) THEN
        IF (LMSPPIO) THEN
        CALL MASAGE ('TEST IONS FROM RECOMBINING BULK IONS          ')
        CALL MASAGE ('PPIO    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS FROM RECOMBINING BULK IONS          ')
        CALL MASYR1 ('PPIOI  = ',
     .                PPIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPIOI(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPIIO) THEN
        IF (LMSPIIO) THEN
        CALL MASAGE ('TEST IONS BORN BY TEST ION - PLASMA INTERACTIONS')
        CALL MASAGE ('PIIO    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PIIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS BORN BY TEST ION - PLASMA INTERACTIONS')
        CALL MASYR1 ('PIIOI  = ',
     .                PIIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PIIOI(0,ISTRA))
      ENDIF
      IF (.NOT.LPAIO) THEN
        IF (LMSPAIO) THEN
        CALL MASAGE ('TEST IONS BORN BY ATOM - PLASMA INTERACTIONS')
        CALL MASAGE ('PAIO    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PAIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS BORN BY ATOM - PLASMA INTERACTIONS')
        CALL MASYR1 ('PAIOI  = ',
     .                PAIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PAIOI(0,ISTRA))
      ENDIF
      IF (.NOT.LPMIO) THEN
        IF (LMSPMIO) THEN
        CALL MASAGE ('TEST IONS BORN BY MOLECULE - PLASMA INTERACTIONS')
        CALL MASAGE ('PMIO    SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PMIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS BORN BY MOLECULE - PLASMA INTERACTIONS')
        CALL MASYR1 ('PMIOI  = ',
     .                PMIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMIOI(0,ISTRA))
      ENDIF
      IF (.NOT.LPPHIO) THEN
        IF (LMSPPHIO) THEN
        CALL MASAGE ('TEST IONS BORN BY PHOTON - PLASMA INTERACTIONS')
        CALL MASAGE ('PPHIO   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPHIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS BORN BY PHOTON - PLASMA INTERACTIONS')
        CALL MASYR1 ('PPHIOI = ',
     .                PPHIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPHIOI(0,ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LPGENI) THEN
        IF (LMSPGENI) THEN
        CALL MASAGE ('TEST ION LOST DUE TO GENERATION-OR FLUID LIMIT')
        CALL MASAGE ('PGENI   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PGENII(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS LOST DUE TO GENERATION-OR FLUID LIMIT')
        CALL MASYR1 ('PGENII = ',
     .                PGENII,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PGENII(0,ISTRA))
      ENDIF
C  ESCAPING TEST IONS
      IF (.NOT.LPOTIO) THEN
        IF (LMSPOTIO) THEN
        CALL MASAGE ('TEST IONS EFFLUX ONTO THE SURFACES             ')
        CALL MASAGE ('POTIO   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(POTIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST IONS EFFLUX ONTO THE SURFACES             ')
        CALL MASYR1 ('POTIOI = ',
     .                POTIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',POTIOI(0,ISTRA))
      ENDIF
C  REFLECTED TEST IONS
      IF (.NOT.LPRFAIO) THEN
        IF (LMSPRFAIO) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: ATOMS     ')
        CALL MASAGE ('PRFAIO  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFAII(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: ATOMS     ')
        CALL MASYR1 ('PRFAII = ',
     .                PRFAII,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFAII(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFMIO) THEN
        IF (LMSPRFMIO) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: MOLECULES ')
        CALL MASAGE ('PRFMIO  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFMII(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: MOLECULES ')
        CALL MASYR1 ('PRFMII = ',
     .                PRFMII,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFMII(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFIIO) THEN
        IF (LMSPRFIIO) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: TEST IONS')
        CALL MASAGE ('PRFIIO  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFIII(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: TEST IONS')
        CALL MASYR1 ('PRFIII = ',
     .                PRFIII,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFIII(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFPHIO) THEN
        IF (LMSPRFPHIO) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: PHOTONS ')
        CALL MASAGE ('PRFPHIO SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFPHII(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('TEST ION INFLUX FROM SURFACES, ORIG: PHOTONS ')
        CALL MASYR1 ('PRFPHII= ',
     .              PRFPHII,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFPHII(0,ISTRA))
      ENDIF
      IF (.NOT.LSPTIO) THEN
        IF (LMSSPTIO) THEN
        CALL MASAGE ('FLUX SPUTTERED BY TEST IONS:                   ')
        CALL MASAGE ('SPTIO   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(SPTIOI(1:NION,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('FLUX SPUTTERED BY TEST IONS:                   ')
        CALL MASYR1 ('SPTIOI = ',
     .                SPTIOI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',SPTIOI(0,ISTRA))
      ENDIF
      DIFI(0,ISTRA)=0.
      TOTI(0)=0.
      DO 425 IION=0,NIONI
        DIFI(IION,ISTRA)=
     +     PIIOI(IION,ISTRA)+POTIOI(IION,ISTRA)+PPIOI(IION,ISTRA)+
     +     WTOTI(IION,ISTRA)+PMIOI(IION,ISTRA)+PAIOI(IION,ISTRA)+
     +     PRFAII(IION,ISTRA)+PRFMII(IION,ISTRA)+
     +     PRFIII(IION,ISTRA)+
     +     PPHIOI(IION,ISTRA)+PRFPHII(IION,ISTRA)
        TOTI(IION)=
     +     ABS(PIIOI(IION,ISTRA))+ABS(POTIOI(IION,ISTRA))+
     +     ABS(PPIOI(IION,ISTRA))+ABS(WTOTI(IION,ISTRA))+
     +     ABS(PMIOI(IION,ISTRA))+ABS(PAIOI(IION,ISTRA))+
     +     ABS(PRFAII(IION,ISTRA))+ABS(PRFMII(IION,ISTRA))+
     +     ABS(PRFIII(IION,ISTRA))+
     +     ABS(PPHIOI(IION,ISTRA))+ABS(PRFPHII(IION,ISTRA))
        TOTI(IION)=TOTI(IION)+EPS60
        DIF=SIGN (1._DP,DIFI(IION,ISTRA))*
     .      MAX(0._DP,ABS(DIFI(IION,ISTRA))/TOTI(IION)-EPS10)
        DIFRI(IION,ISTRA)=DIF*100.
        DIFI(IION,ISTRA)=DIF*TOTI(IION)
425   CONTINUE
      CALL MASAGE ('ABSOLUTE ERRORS IN PARTICLE BALANCE            ')
      CALL MASYR1 ('DIFI =   ',
     .              DIFI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFI(0,ISTRA))
      CALL MASAGE ('RELATIVE ERRORS IN PARTICLE BALANCE (%)        ')
      CALL MASYR1 ('DIFRI =  ',
     .              DIFRI,LOGION,ISTRA,0,NION,0,NSTRA,TEXTS(NSPAM+1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFRI(0,ISTRA))
      CALL LEER(2)
C
C  PHOTONS PARTICLE BALANCE
C
430   CONTINUE
      DIFPH(0,ISTRA)=0.
      TOTPH(0)=0.
      IF (.NOT.TRCBLPH) GOTO 439
      IF (.NOT.LOGPHOT(0,ISTRA)) GOTO 439
      CALL LEER(2)
      CALL HEADNG('PARTICLE FLUX BALANCE (AMP), PHOTONS',36)
      CALL LEER(1)
C   PHOTONS FROM PRIMARY SOURCE
      IF (ANY(WTOTPH(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM PRIMARY SOURCE            ')
        CALL MASYR1 ('WTOTPH = ',
     .                WTOTPH,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',WTOTPH(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPPPHT) THEN
        IF (LMSPPPHT) THEN
        CALL MASAGE ('PHOTONS FROM RECOMBINING BULK IONS             ')
        CALL MASAGE ('PPPHT   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS FROM RECOMBINING BULK IONS             ')
        CALL MASYR1 ('PPPHTI = ',
     .                PPPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPPHTI(0,ISTRA))
      ENDIF
C
      IF (.NOT.LPPHPHT) THEN
        IF (LMSPPHPHT) THEN
        CALL MASAGE ('PHOTONS BORN BY PHOTON - PLASMA INTERACTIONS   ')
        CALL MASAGE ('PPHPHT  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PPHPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS BORN BY PHOTON - PLASMA INTERACTIONS   ')
        CALL MASYR1 ('PPHPHTI= ',
     .                PPHPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPHPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPAPHT) THEN
        IF (LMSPAPHT) THEN
        CALL MASAGE ('PHOTONS BORN BY ATOM - PLASMA INTERACTIONS     ')
        CALL MASAGE ('PAPHT   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PAPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS BORN BY ATOM - PLASMA INTERACTIONS     ')
        CALL MASYR1 ('PAPHTI = ',
     .                PAPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PAPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPMPHT) THEN
        IF (LMSPMPHT) THEN
        CALL MASAGE ('PHOTONS BORN BY MOLECULE - PLASMA INTERACTIONS ')
        CALL MASAGE ('PMPHT   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PMPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS BORN BY MOLECULE - PLASMA INTERACTIONS ')
        CALL MASYR1 ('PMPHTI = ',
     .                PMPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPIPHT) THEN
        IF (LMSPIPHT) THEN
        CALL MASAGE ('PHOTONS BORN BY TEST ION - PLASMA INTERACTIONS ')
        CALL MASAGE ('PIPHT   SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PIPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS BORN BY TEST ION - PLASMA INTERACTIONS ')
        CALL MASYR1 ('PIPHTI = ',
     .                PIPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PIPHTI(0,ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LPGENPH) THEN
        IF (LPGENPH) THEN
        CALL MASAGE ('PHOTONS ABSORBED DUE TO GENERAT.-OR FLUID LIMIT')
        CALL MASAGE ('PGENPH  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PGENPHI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONS ABSORBED DUE TO GENERAT.-OR FLUID LIMIT')
        CALL MASYR1 ('PGENPHI= ',
     .                PGENPHI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PGENPHI(0,ISTRA))
      ENDIF
C   ESCAPING FLUX
      IF (.NOT.LPOTPHT) THEN
        IF (LMSPOTPHT) THEN
        CALL MASAGE ('PHOTONIC EFFLUX ONTO THE SURFACES              ')
        CALL MASAGE ('POTPHT  SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(POTPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC EFFLUX ONTO THE SURFACES              ')
        CALL MASYR1 ('POTPHTI= ',
     .                POTPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',POTPHTI(0,ISTRA))
      ENDIF
C   REFLECTED FLUX
      IF (.NOT.LPRFAPHT) THEN
        IF (LMSPRFAPHT) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: ATOMS     ')
        CALL MASAGE ('PRFAPHT SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFAPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: ATOMS     ')
        CALL MASYR1 ('PRFAPHTI=',
     .                PRFAPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFAPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFMPHT) THEN
        IF (LMSPRFMPHT) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: MOLECULES ')
        CALL MASAGE ('PRFMPHT SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFMPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: MOLECULES ')
        CALL MASYR1 ('PRFMPHTI=',
     .                PRFMPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFMPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFIPHT) THEN
        IF (LMSPRFIPHT) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: TEST IONS ')
        CALL MASAGE ('PRFIPHT SWITCHED OFF => MISSING IN BALANCE     ')
        ENDIF
      ELSE IF (ANY(PRFIPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: TEST IONS ')
        CALL MASYR1 ('PRFIPHTI = ',
     .                PRFIPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFIPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LPRFPHPHT) THEN
        IF (LMSPRFPHPHT) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: PHOTONS   ')
        CALL MASAGE ('PRFPHPHT SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(PRFPHPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('PHOTONIC INFLUX FROM SURFACES, ORIG: PHOTONS   ')
        CALL MASYR1 ('PRFPHPHTI',
     .                PRFPHPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PRFPHPHTI(0,ISTRA))
      ENDIF
      IF (.NOT.LSPTPHT) THEN
        IF (LMSSPTPHT) THEN
        CALL MASAGE ('FLUX SPUTTERED BY PHOTONS:                     ')
        CALL MASAGE ('SPTPHT   SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(SPTPHTI(1:NPHOT,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('FLUX SPUTTERED BY PHOTONS:                     ')
        CALL MASYR1 ('SPTPHTI = ',
     .                SPTPHTI,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS( +1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',SPTPHTI(0,ISTRA))
      ENDIF
      DO 435 IPHOT=0,NPHOTI
        DIFPH(IPHOT,ISTRA)=
     +  POTPHTI(IPHOT,ISTRA)+PAPHTI(IPHOT,ISTRA)+PPPHTI(IPHOT,ISTRA)+
     +       PMPHTI(IPHOT,ISTRA)+PIPHTI(IPHOT,ISTRA)+
     +       PGENPHI(IPHOT,ISTRA)+
     +       WTOTPH(IPHOT,ISTRA)+
     +       PPHPHTI(IPHOT,ISTRA)+
     +       PRFAPHTI(IPHOT,ISTRA)+
     +       PRFMPHTI(IPHOT,ISTRA)+
     +       PRFIPHTI(IPHOT,ISTRA)+
     +       PRFPHPHTI(IPHOT,ISTRA)
        TOTPH(IPHOT)=
     +  ABS(POTPHTI(IPHOT,ISTRA))+ABS(PAPHTI(IPHOT,ISTRA))+
     +       ABS(PPPHTI(IPHOT,ISTRA))+ABS(PMPHTI(IPHOT,ISTRA))+
     +       ABS(PIPHTI(IPHOT,ISTRA))+ABS(WTOTPH(IPHOT,ISTRA))+
     +       ABS(PRFAPHTI(IPHOT,ISTRA))+ABS(PRFMPHTI(IPHOT,ISTRA))+
     +       ABS(PRFIPHTI(IPHOT,ISTRA))+
     +       ABS(PGENPHI(IPHOT,ISTRA))+
     +       ABS(PPHPHTI(IPHOT,ISTRA))+
     +       ABS(PRFPHPHTI(IPHOT,ISTRA))
        TOTPH(IPHOT)=TOTPH(IPHOT)+EPS60
        DIF=SIGN (1._DP,DIFPH(IPHOT,ISTRA))*
     .      MAX(0._DP,ABS(DIFPH(IPHOT,ISTRA))/TOTPH(IPHOT)-EPS10)
        DIFRPH(IPHOT,ISTRA)=DIF*100.
        DIFPH(IPHOT,ISTRA)=DIF*TOTPH(IPHOT)
435   CONTINUE
      CALL MASAGE ('ABSOLUTE ERRORS IN PARTICLE BALANCE            ')
      CALL MASYR1 ('DIFPH=   ',
     .              DIFPH,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(   +1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFPH(0,ISTRA))
      CALL MASAGE ('RELATIVE ERRORS IN PARTICLE BALANCE (%)        ')
      CALL MASYR1 ('DIFRPH=  ',
     .              DIFRPH,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,TEXTS(  +1))
      CALL MASAGE ('SUM OVER SPECIES                               ')
      CALL MASR1 ('TOTAL=  ',DIFRPH(0,ISTRA))
      CALL LEER(2)
C
C  TOTAL TEST PARTICLE BALANCE
C
439   CONTINUE
      CALL MASAGE ('TRASH: TEST PARTICLES KILLED DUE TO ERRORS     ')
      CALL MASR1 ('PTRASH= ',PTRASH(ISTRA))
      CALL LEER(3)
      IF (TRCBLA.OR.TRCBLM.OR.TRCBLI) THEN
        CALL MASAGE ('TOTAL ERROR IN PARTICLE FLUXBALANCE            ')
        DIFT=DIFA(0,ISTRA)+DIFM(0,ISTRA)+DIFI(0,ISTRA)+DIFPH(0,ISTRA)+
     .       PTRASH(ISTRA)
        TOTT=TOTA(0)+TOTM(0)+TOTI(0)+TOTPH(0)+ABS(PTRASH(ISTRA))+EPS60
        DIFR=DIFT/TOTT*100.
        CALL MASR2 ('DIF: ABS, REL(%)',DIFT,DIFR)
        CALL LEER(2)
      ENDIF
C
C
C   ENERGY FLUX BALANCE,  ATOMS
C
      DIFA(0,ISTRA)=0.
      TOTA(0)=0.
      IF (.NOT.TRCBLA) GOTO 440
      IF (.NOT.LOGATM(0,ISTRA)) GOTO 440
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), NEUTRAL ATOMS',41)
      CALL LEER(1)
      CALL MASAGE ('ENERGY FLUX FROM PRIMARY SOURCE                ')
      CALL MASR1 ('ETOTA=  ',ETOTA(ISTRA))
      CALL MASAGE ('ENERGY FLUX FROM RECOMBINING BULK IONS         ')
      CALL MASR1 ('EPATI=  ',EPATI(ISTRA))
      IF (LOGATM(0,ISTRA)) THEN
        IF (LEAAT) THEN
          CALL MASAGE('ENERGY FLUX FROM ATOM PLASMA INTERACTION       ')
          CALL MASR1 ('EAATI=  ',EAATI(ISTRA))
        ELSE IF (LMSEAAT) THEN
          CALL MASAGE('ENERGY FLUX FROM ATOM PLASMA INTERACTION       ')
          CALL MASAGE
     .      ('EAAT     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ENDIF
      IF (LOGMOL(0,ISTRA)) THEN
        IF (LEMAT) THEN
          CALL MASAGE('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
          CALL MASR1 ('EMATI=  ',EMATI(ISTRA))
        ELSE IF (LMSEMAT) THEN
          CALL MASAGE('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
          CALL MASAGE
     .      ('EMAT     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ENDIF
      IF (LOGION(0,ISTRA)) THEN
        IF (LEIAT) THEN
          CALL MASAGE('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS  ')
          CALL MASR1 ('EIATI=  ',EIATI(ISTRA))
        ELSE IF (LMSEIAT) THEN
          CALL MASAGE('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS  ')
          CALL MASAGE
     .      ('EIAT     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ENDIF
      IF (LOGPHOT(0,ISTRA)) THEN
        IF (LEPHAT) THEN
          CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTIONS')
          CALL MASR1 ('EPHATI= ',EPHATI(ISTRA))
        ELSE IF (LMSEPHAT) THEN
          CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTIONS')
          CALL MASAGE
     .      ('EPHAT    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LEGENA) THEN
        IF (LMSEGENA) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASAGE ('EGENA    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EGENAI(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASR1 ('EGENAI= ',EGENAI(0,ISTRA))
      ENDIF
C   EFFLUX TO SURFACES
      IF (LEOTAT) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASR1 ('EOTATI= ',EOTATI(0,ISTRA))
      ELSE IF (LMSEOTAT) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASAGE ('EOTAT    SWITCHED OFF => MISSING IN BALANCE    ')
      ENDIF
C   REFLECTED FROM SURFACES
      CALL MASAGE ('REFLECTED FROM NON TRANSPARENT SURFACES         ')
      CALL MASR1 ('ERFATI= ',ERFAAI(0,ISTRA)+ERFMAI(0,ISTRA)+
     +                       ERFIAI(0,ISTRA)+ERFPHAI(0,ISTRA))
      IF (LMSERFAAT) CALL MASAGE
     .      ('ERFAAT SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFMAT) CALL MASAGE
     .      ('ERFMAT SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFIAT) CALL MASAGE
     .      ('ERFIAT SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFPHAT) CALL MASAGE
     .      ('ERFPHAT SWITCHED OFF => MISSING IN BALANCE     ')
      CALL MASAGE ('ABSOLUTE AND RELATIVE ERROR IN BALANCE          ')
      DIFA(0,ISTRA)=EOTATI(0,ISTRA)+EPATI(ISTRA)+EAATI(ISTRA)+
     +     ETOTA(ISTRA)+EMATI(ISTRA)+EIATI(ISTRA)+EGENAI(0,ISTRA)+
     +     ERFAAI(0,ISTRA)+ERFMAI(0,ISTRA)+ERFIAI(0,ISTRA)+
     +     EPHATI(ISTRA)+ERFPHAI(0,ISTRA)

      TOTA(0)=ABS(EOTATI(0,ISTRA))+ABS(EPATI(ISTRA))+ABS(EAATI(ISTRA))+
     +     ABS(ETOTA(ISTRA))+ABS(EMATI(ISTRA))+ABS(EIATI(ISTRA))+
     +     ABS(EGENAI(0,ISTRA))+
     +     ABS(ERFAAI(0,ISTRA)+ERFMAI(0,ISTRA)+ERFIAI(0,ISTRA))+
     +     ABS(EPHATI(ISTRA))+ABS(ERFPHAI(0,ISTRA))

      TOTA(0)=TOTA(0)+EPS60
      DIFRA(0,ISTRA)=SIGN(1._DP,DIFA(0,ISTRA))*
     *        MAX(0._DP,ABS(DIFA(0,ISTRA))/TOTA(0)*100._DP-1.E-5_DP)
      DIFA(0,ISTRA)=DIFRA(0,ISTRA)/100.*TOTA(0)
      CALL MASR2 ('DIF: ABS, REL(%)',DIFA(0,ISTRA),DIFRA(0,ISTRA))
      CALL LEER(2)
C
440   CONTINUE
C   ENERGY FLUX BALANCE,  MOLECULES
      DIFM(0,ISTRA)=0.
      TOTM(0)=0.
      IF (.NOT.TRCBLM) GOTO 450
      IF (.NOT.LOGMOL(0,ISTRA)) GOTO 450
C
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), NEUTRAL MOLECULES',45)
      CALL LEER(1)
      CALL MASAGE ('FROM PRIMARY SOURCE                            ')
      CALL MASR1 ('ETOTM = ',ETOTM(ISTRA))
      IF (.NOT.LEPML) THEN
        IF (LMSEPML) THEN
        CALL MASAGE ('FROM RECOMBINING BULK IONS                     ')
        CALL MASAGE ('EPML SWITCHED OFF => MISSING IN BALANCE        ')
        ENDIF
      ELSEIF (EPMLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('FROM RECOMBINING BULK IONS                     ')
        CALL MASR1 ('EPMLI=  ',EPMLI(ISTRA))
      END IF
      IF (.NOT.LEAML) THEN
        IF (LMSEAML) THEN
        CALL MASAGE ('ENERGY FLUX FROM ATOM PLASMA INTERACTIONS     ')
        CALL MASAGE ('EAML     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EAMLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM ATOM PLASMA INTERACTIONS     ')
        CALL MASR1 ('EAMLI = ',EAMLI(ISTRA))
      ENDIF
      IF (LEMML) THEN
        CALL MASAGE ('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASR1 ('EMMLI = ',EMMLI(ISTRA))
      ELSE IF (LMSEMML) THEN
        CALL MASAGE ('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASAGE ('EMML SWITCHED OFF => MISSING IN BALANCE        ')
      END IF
      IF (.NOT.LEIML) THEN
        IF (LMSEIML) THEN
        CALL MASAGE ('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS ')
        CALL MASAGE ('EIML     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EIMLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS ')
        CALL MASR1 ('EIMLI = ',EIMLI(ISTRA))
      ENDIF
      IF (.NOT.LEPHML) THEN
        IF (LMSEPHML) THEN
        CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTIONS')
        CALL MASAGE ('EPHML    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPHMLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTIONS')
        CALL MASR1 ('EPHMLI= ',EPHMLI(ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LEGENM) THEN
        IF (LMSEGENM) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASAGE ('EGENM    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EGENMI(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASR1 ('EGENMI= ',EGENMI(0,ISTRA))
      ENDIF
C   EFFLUX ONTO SURFACES
      IF (LEOTML) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASR1 ('EOTMLI= ',EOTMLI(0,ISTRA))
      ELSE IF (LMSEOTML) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASAGE ('EOTML SWITCHED OFF => MISSING IN BALANCE       ')
      END IF
      CALL MASAGE ('REFLECTED FROM NON TRANSPARENT SURFACES        ')
      CALL MASR1 ('ERFMLI= ',ERFAMI(0,ISTRA)+ERFMMI(0,ISTRA)+
     +                       ERFIMI(0,ISTRA)+ERFPHMI(0,ISTRA))
      IF (LMSERFAML) CALL MASAGE
     .      ('ERFAML SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFMML) CALL MASAGE
     .      ('ERFMML SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFIML) CALL MASAGE
     .      ('ERFIML SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFPHML) CALL MASAGE
     .      ('ERFPHML SWITCHED OFF => MISSING IN BALANCE     ')
      CALL MASAGE ('ABSOLUTE AND RELATIVE ERROR IN BALANCE          ')
      DIFM(0,ISTRA)=
     +     EAMLI(ISTRA)+EOTMLI(0,ISTRA)+EIMLI(ISTRA)+EMMLI(ISTRA)+
     +     ETOTM(ISTRA)+EPMLI(ISTRA)+
     +     ERFAMI(0,ISTRA)+ERFMMI(0,ISTRA)+ERFIMI(0,ISTRA)+
     +     EPHMLI(ISTRA)+ERFPHMI(0,ISTRA)

      TOTM(0)=ABS(EAMLI(ISTRA))+ABS(EOTMLI(0,ISTRA))+ABS(EIMLI(ISTRA))+
     +        ABS(EMMLI(ISTRA))+ABS(ETOTM(ISTRA))+ABS(EPMLI(ISTRA))+
     +        ABS(ERFAMI(0,ISTRA)+ERFMMI(0,ISTRA)+ERFIMI(0,ISTRA))+
     +        ABS(EPHMLI(ISTRA))+ABS(ERFPHMI(0,ISTRA))

      TOTM(0)=TOTM(0)+EPS60
      DIFRM(0,ISTRA)=SIGN(1._DP,DIFM(0,ISTRA))*
     *        MAX(0._DP,ABS(DIFM(0,ISTRA))/TOTM(0)*100._DP-1.E-5_DP)
      DIFM(0,ISTRA)=DIFRM(0,ISTRA)/100.*TOTM(0)
      CALL MASR2 ('DIF: ABS, REL(%)',DIFM(0,ISTRA),DIFRM(0,ISTRA))
      CALL LEER(2)
C
450   CONTINUE
C   ENERGY FLUX BALANCE,  TEST IONS
      DIFI(0,ISTRA)=0.
      TOTI(0)=0.
      IF (.NOT.TRCBLI) GOTO 460
      IF (.NOT.LOGION(0,ISTRA)) GOTO 460
C
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), TEST IONS',37)
      CALL LEER(1)
      CALL MASAGE ('FROM PRIMARY SOURCE                            ')
      CALL MASR1 ('ETOTI = ',ETOTI(ISTRA))
      IF (.NOT.LEPIO) THEN
        IF (LMSEPIO) THEN
        CALL MASAGE ('FROM RECOMBINING BULK IONS                     ')
        CALL MASAGE ('EPIO SWITCHED OFF => MISSING IN BALANCE        ')
        ENDIF
      ELSEIF (EPIOI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('FROM RECOMBINING BULK IONS                     ')
        CALL MASR1 ('EPIOI=  ',EPIOI(ISTRA))
      END IF
      IF (.NOT.LEAIO) THEN
        IF (LMSEAIO) THEN
        CALL MASAGE ('ENERGY GAINED FROM ATOM PLASMA INTERACTIONS   ')
        CALL MASAGE ('EAIO     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EAIOI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAINED FROM ATOM PLASMA INTERACTIONS   ')
        CALL MASR1 ('EAIOI = ',EAIOI(ISTRA))
      ENDIF
      IF (LEMIO) THEN
        CALL MASAGE ('ENERGY GAINED FROM MOLECULE PLASMA INTERACTIONS')
        CALL MASR1 ('EMIOI = ',EMIOI(ISTRA))
      ELSE IF (LMSEMIO) THEN
        CALL MASAGE ('ENERGY GAINED FROM MOLECULE PLASMA INTERACTIONS')
        CALL MASAGE ('EMIO SWITCHED OFF => MISSING IN BALANCE        ')
      END IF
      IF (.NOT.LEIIO) THEN
        IF (LMSEIIO) THEN
        CALL MASAGE ('ENERGY GAINED FROM TEST ION PLASMA INTERACTIONS')
        CALL MASAGE ('EIIO     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EIIOI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAINED FROM TEST ION PLASMA INTERACTIONS')
        CALL MASR1 ('EIIOI = ',EIIOI(ISTRA))
      ENDIF
      IF (.NOT.LEPHIO) THEN
        IF (LMSEPHIO) THEN
        CALL MASAGE ('ENERGY GAINED FROM PHOTON PLASMA INTERACTIONS')
        CALL MASAGE ('EPHIO    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPHIOI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAINED FROM PHOTON PLASMA INTERACTIONS')
        CALL MASR1 ('EPHIOI= ',EPHIOI(ISTRA))
      ENDIF
      IF (EELFI(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAINED FROM ELECTRIC FIELDS')
        CALL MASR1 ('EELFI = ',EELFI(0,ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LEGENI) THEN
        IF (LMSEGENI) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASAGE ('EGENI    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EGENII(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASR1 ('EGENII= ',EGENII(0,ISTRA))
      ENDIF
C   EFFLUX ONTO SURFACES
      IF (.NOT.LEOTIO) THEN
        IF (LMSEOTIO) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASAGE ('EOTIOI   SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EOTIOI(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES      ')
        CALL MASR1 ('EOTIOI= ',EOTIOI(0,ISTRA))
      ENDIF
C   REFLECTED FROM SURFACES
      CALL MASAGE ('REFLECTED FROM NON TRANSPARENT SURFACES        ')
      CALL MASR1 ('ERFIOI= ',ERFAII(0,ISTRA)+ERFMII(0,ISTRA)+
     +                       ERFIII(0,ISTRA)+ERFPHII(0,ISTRA))
      IF (LMSERFAIO) CALL MASAGE
     .      ('ERFAIO SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFMIO) CALL MASAGE
     .      ('ERFMIO SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFIIO) CALL MASAGE
     .      ('ERFIIO SWITCHED OFF => MISSING IN BALANCE      ')
      IF (LMSERFPHIO) CALL MASAGE
     .      ('ERFPHIO SWITCHED OFF => MISSING IN BALANCE     ')
C
      CALL MASAGE ('ABSOLUTE AND RELATIVE ERROR IN BALANCE          ')
      DIFI(0,ISTRA)=
     +     EAIOI(ISTRA)+EOTIOI(0,ISTRA)+EIIOI(ISTRA)+EMIOI(ISTRA)+
     +     ETOTI(ISTRA)+EPIOI(ISTRA)+EELFI(0,ISTRA)+
     +     ERFAII(0,ISTRA)+ERFMII(0,ISTRA)+ERFIII(0,ISTRA)+
     +     EPHIOI(ISTRA)+ERFPHII(0,ISTRA)

      TOTI(0)=ABS(EAIOI(ISTRA))+ABS(EOTIOI(0,ISTRA))+ABS(EIIOI(ISTRA))+
     +        ABS(EMIOI(ISTRA))+ABS(ETOTI(ISTRA))+ABS(EPIOI(ISTRA))+
     +        ABS(ERFAII(0,ISTRA)+ERFMII(0,ISTRA)+ERFIII(0,ISTRA))+
     +        ABS(EPHIOI(ISTRA))+ABS(ERFPHII(0,ISTRA))

      TOTI(0)=TOTI(0)+EPS60
      DIFRI(0,ISTRA)=SIGN(1._DP,DIFI(0,ISTRA))*
     *        MAX(0._DP,ABS(DIFI(0,ISTRA))/TOTI(0)*100._DP-1.E-5_DP)
      DIFI(0,ISTRA)=DIFRI(0,ISTRA)/100.*TOTI(0)
      CALL MASR2 ('DIF: ABS, REL(%)',DIFI(0,ISTRA),DIFRI(0,ISTRA))
      CALL LEER(2)

460   CONTINUE
C
C   ENERGY FLUX BALANCE,  PHOTONS
C
      DIFPH(0,ISTRA)=0.
      TOTPH(0)=0.
      IF (.NOT.TRCBLPH) GOTO 469
      IF (.NOT.LOGPHOT(0,ISTRA)) GOTO 469
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), PHOTONS',36)
      CALL LEER(1)
      CALL MASAGE ('ENERGY FLUX FROM PRIMARY SOURCE                ')
      CALL MASR1 ('ETOTPH= ',ETOTPH(ISTRA))
      IF (LEPPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM RECOMBINING BULK IONS         ')
        CALL MASR1 ('EPPHTI= ',EPPHTI(ISTRA))
      ELSE IF (LMSEPPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM RECOMBINING BULK IONS         ')
        CALL MASAGE ('EPPHT SWITCHED OFF => MISSING IN BALANCE       ')
      END IF
      IF (.NOT.LEPHPHT) THEN
        IF (LMSEPHPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTION    ')
        CALL MASAGE ('EPHPHT   SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPHPHTI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM PHOTON PLASMA INTERACTION    ')
        CALL MASR1 ('EPHPHTI=',EPHPHTI(ISTRA))
      ENDIF
      IF (.NOT.LEAPHT) THEN
        IF (LMSEAPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM ATOM PLASMA INTERACTION       ')
        CALL MASAGE ('EAPHT    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EAPHTI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM ATOM PLASMA INTERACTION       ')
        CALL MASR1 ('EAPHTI= ',EAPHTI(ISTRA))
      ENDIF
      IF (.NOT.LEMPHT) THEN
        IF (LMSEMPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASAGE ('EMPHT    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EMPHTI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASR1 ('EMPHTI= ',EMPHTI(ISTRA))
      ENDIF
      IF (.NOT.LEIPHT) THEN
        IF (LMSEIPHT) THEN
        CALL MASAGE ('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS  ')
        CALL MASAGE ('EIPHT    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EIPHTI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY FLUX FROM TEST ION PLASMA INTERACTIONS  ')
        CALL MASR1 ('EIPHTI= ',EIPHTI(ISTRA))
      ENDIF
C   GENERATION LIMIT
      IF (.NOT.LEGENPH) THEN
        IF (LMSEGENPH) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASAGE ('EGENPH   SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EGENPHI(0,ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY ABSORBED DUE TO GENERATION LIMIT        ')
        CALL MASR1 ('EGENPHI=',EGENPHI(0,ISTRA))
      ENDIF
      IF (LEOTPHT) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES       ')
        CALL MASR1 ('EOTPHTI=',EOTPHTI(0,ISTRA))
      ELSE IF (LMSEOTPHT) THEN
        CALL MASAGE ('ENERGY FLUX ONTO NON TRANSPARENT SURFACES       ')
        CALL MASAGE ('EOTPHT SWITCHED OFF => MISSING IN BALANCE      ')
      END IF
      CALL MASAGE ('REFLECTED FROM NON TRANSPARENT SURFACES         ')
      CALL MASR1 ('ERFPHTI=',ERFAPHTI(0,ISTRA)+ERFMPHTI(0,ISTRA)+
     +                       ERFIPHTI(0,ISTRA)+ERFPHPHTI(0,ISTRA))
      IF (LMSERFAPHT) CALL MASAGE
     .      ('ERFAPHT SWITCHED OFF => MISSING IN BALANCE     ')
      IF (LMSERFMPHT) CALL MASAGE
     .      ('ERFMPHT SWITCHED OFF => MISSING IN BALANCE     ')
      IF (LMSERFIPHT) CALL MASAGE
     .      ('ERFIPHT SWITCHED OFF => MISSING IN BALANCE     ')
      IF (LMSERFPHPHT) CALL MASAGE
     .      ('ERFPHPHT SWITCHED OFF => MISSING IN BALANCE    ')
      CALL MASAGE ('ABSOLUTE AND RELATIVE ERROR IN BALANCE          ')
      DIFPH(0,ISTRA)=EOTPHTI(0,ISTRA)+EPPHTI(ISTRA)+EAPHTI(ISTRA)+
     +     ETOTPH(ISTRA)+EMPHTI(ISTRA)+EIPHTI(ISTRA)+EGENPHI(0,ISTRA)+
     +     ERFAPHTI(0,ISTRA)+ERFMPHTI(0,ISTRA)+ERFIPHTI(0,ISTRA)+
     +     EPHPHTI(ISTRA)+ERFPHPHTI(0,ISTRA)

      TOTPH(0)=ABS(EOTPHTI(0,ISTRA))+ABS(EPPHTI(ISTRA))+
     +     ABS(EAPHTI(ISTRA))+
     +     ABS(ETOTPH(ISTRA))+ABS(EMPHTI(ISTRA))+ABS(EIPHTI(ISTRA))+
     +     ABS(EGENPHI(0,ISTRA))+
     +     ABS(ERFAPHTI(0,ISTRA)+ERFMPHTI(0,ISTRA)+ERFIPHTI(0,ISTRA))+
     +     ABS(EPHPHTI(ISTRA))+ABS(ERFPHPHTI(0,ISTRA))

      TOTPH(0)=TOTPH(0)+EPS60
      DIFRPH(0,ISTRA)=SIGN(1._DP,DIFPH(0,ISTRA))*
     *        MAX(0._DP,ABS(DIFPH(0,ISTRA))/TOTPH(0)*100._DP-1.E-5_DP)
      DIFPH(0,ISTRA)=DIFRPH(0,ISTRA)/100.*TOTPH(0)
      CALL MASR2 ('DIF: ABS, REL(%)',DIFPH(0,ISTRA),DIFRPH(0,ISTRA))
      CALL LEER(2)
C
C  TOTAL TEST PARTICLE ENERGY FLUX BALANCE
C
469   CONTINUE
      CALL MASAGE ('ETRASH: ENERGY ABSORBED DUE TO ERRORS          ')
      CALL MASR1 ('ETRASH= ',ETRASH(ISTRA))
      CALL LEER(3)
      IF (TRCBLA.OR.TRCBLM.OR.TRCBLI) THEN
        CALL MASAGE ('TOTAL ERROR IN ENERGY FLUXBALANCE              ')
        DIFT=DIFA(0,ISTRA)+DIFM(0,ISTRA)+DIFI(0,ISTRA)+DIFPH(0,ISTRA)+
     +       ETRASH(ISTRA)

        TOTT=TOTA(0)+TOTM(0)+TOTI(0)+TOTPH(0)+ABS(ETRASH(ISTRA))+EPS60
        DIFR=DIFT/TOTT*100.
        CALL MASR2 ('DIF: ABS, REL(%)',DIFT,DIFR)
        CALL LEER(2)
      ENDIF
C
600   CONTINUE
C
C  INTEGRATED BULK ION SOURCE TERMS
C
      IF (.NOT.TRCBLP) GOTO 700
      IF (.NOT.LOGPLS(0,ISTRA)) GOTO 700
      CALL HEADNG ('PARTICLE FLUX BALANCE (AMP), BULK IONS',39)
      CALL LEER(1)
C  PRIMARY SOURCE ORIGINATING FROM BULK PARTICLES
      IF (ANY(WTOTP(1:NPLS,ISTRA).NE.0.D0)) THEN
        IF (.NOT.LPPPL) THEN
        IF (LMSPPPL) THEN
        CALL MASAGE ('BULK PARTICLE FLUX                             ')
        CALL MASAGE ('(DEFINING THE TEST PARTICLE SOURCE             ')
        CALL MASAGE ('STRENGTH (AMP))                                ')
        CALL MASAGE ('PPPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
        ELSE IF (ANY(PPPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('BULK PARTICLE FLUX                             ')
        CALL MASAGE ('(DEFINING THE TEST PARTICLE SOURCE             ')
        CALL MASAGE ('STRENGTH (AMP))                                ')
        CALL MASYR1 ('PPPLI =  ',PPPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPPLI(0,ISTRA))
        ENDIF
      ENDIF
C  ATOMS PLASMA INTERACTION
      IF (.NOT.LPAPL) THEN
        IF (LMSPAPL) THEN
        CALL MASAGE ('BULK PARTICLES BORN BY ATOM PLASMA INTERACTIONS')
        CALL MASAGE ('PAPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(PAPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('BULK PARTICLES BORN BY ATOM PLASMA INTERACTIONS')
        CALL MASYR1 ('PAPLI  = ',PAPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PAPLI(0,ISTRA))
      ENDIF
C  MOLECULE PLASMA INTERACTION
      IF (.NOT.LPMPL) THEN
        IF (LMSPMPL) THEN
        CALL MASAGE ('BULK IONS BORN BY MOLECULE PLASMA INTERACTIONS')
        CALL MASAGE ('PMPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(PMPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('BULK IONS BORN BY MOLECULE PLASMA INTERACTIONS')
        CALL MASYR1 ('PMPLI  = ',PMPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PMPLI(0,ISTRA))
      ENDIF
C  TEST ION PLASMA INTERACTION
      IF (.NOT.LPIPL) THEN
        IF (LMSPIPL) THEN
        CALL MASAGE ('BULK IONS BORN BY TEST ION PLASMA INTERACTIONS')
        CALL MASAGE ('PIPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(PIPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('BULK IONS BORN BY TEST ION PLASMA INTERACTIONS')
        CALL MASYR1 ('PIPLI  = ',PIPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PIPLI(0,ISTRA))
      ENDIF
C  PHOTON PLASMA INTERACTION
      IF (.NOT.LPPHPL) THEN
        IF (LMSPPHPL) THEN
        CALL MASAGE ('BULK IONS BORN BY PHOTON PLASMA INTERACTIONS')
        CALL MASAGE ('PPHPL    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(PPHPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('BULK IONS BORN BY PHOTON PLASMA INTERACTIONS')
        CALL MASYR1 ('PPHPLI = ',PPHPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',PPHPLI(0,ISTRA))
      ENDIF
      IF (.NOT.LSPTPL) THEN
        IF (LMSSPTPL) THEN
        CALL MASAGE ('FLUX SPUTTERED BY BULK IONS:                   ')
        CALL MASAGE ('SPTPL    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (ANY(SPTPLI(1:NPLS,ISTRA).NE.0.D0)) THEN
        CALL MASAGE ('FLUX SPUTTERED BY BULK IONS:                   ')
        CALL MASYR1 ('SPTPLI = ',SPTPLI,
     .                LOGPLS,ISTRA,0,NPLS,0,NSTRA,TEXTS(NSPAMI+1))
        CALL MASAGE ('SUM OVER SPECIES                               ')
        CALL MASR1 ('TOTAL=  ',SPTPLI(0,ISTRA))
      ENDIF
      CALL LEER(2)
      PLOSSP=PPPLI(0,ISTRA)
      PGAINP=PAPLI(0,ISTRA)+PMPLI(0,ISTRA)+PIPLI(0,ISTRA)+
     +       PPHPLI(0,ISTRA)

      CALL MASAGE ('TOTAL LOSS  ---  TOTAL GAIN                    ')
      CALL MASR2 ('LOSS--GAIN      ',PLOSSP,PGAINP)
      CALL LEER(2)
C
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), BULK IONS',37)
      CALL LEER(1)
      IF (ETOTP(ISTRA).NE.0.D0) THEN
      IF (.NOT.LEPPL) THEN
        IF (LMSEPPL) THEN
        CALL MASAGE ('BULK ION ENERGY FLUX  BEING                    ')
        CALL MASAGE ('NEUTRALIZED  (DEFINING THE TEST PARTICLES      ')
        CALL MASAGE ('SOURCE STRENGTH (WATT))                        ')
        CALL MASAGE ('EPPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPPLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('BULK ION ENERGY FLUX  BEING                    ')
        CALL MASAGE ('NEUTRALIZED  (DEFINING THE TEST PARTICLES      ')
        CALL MASAGE ('SOURCE STRENGTH (WATT))                        ')
        CALL MASR1 ('EPPLI=  ',EPPLI(ISTRA))
      ENDIF
      ENDIF
      IF (.NOT.LEAPL) THEN
        IF (LMSEAPL) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO ATOM PLASMA INTERACTIONS  ')
        CALL MASAGE ('EAPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EAPLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO ATOM PLASMA INTERACTIONS  ')
        CALL MASR1 ('EAPLI=  ',EAPLI(ISTRA))
      ENDIF
      IF (.NOT.LEMPL) THEN
        IF (LMSEMPL) THEN
        CALL MASAGE ('ENERGY GAIN FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASAGE ('EMPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EMPLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN FROM MOLECULE PLASMA INTERACTIONS  ')
        CALL MASR1 ('EMPLI=  ',EMPLI(ISTRA))
      ENDIF
      IF (.NOT.LEIPL) THEN
        IF (LMSEIPL) THEN
        CALL MASAGE ('ENERGY GAIN FROM TEST ION PLASMA INTERACTIONS  ')
        CALL MASAGE ('EIPL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EIPLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN FROM TEST ION PLASMA INTERACTIONS  ')
        CALL MASR1 ('EIPLI=  ',EIPLI(ISTRA))
      ENDIF
      IF (.NOT.LEPHPL) THEN
        IF (LMSEPHPL) THEN
        CALL MASAGE ('ENERGY GAIN FROM PHOTON PLASMA INTERACTIONS  ')
        CALL MASAGE ('EPHPL    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPHPLI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN FROM PHOTON PLASMA INTERACTIONS  ')
        CALL MASR1 ('EPHPLI= ',EPHPLI(ISTRA))
      ENDIF
      CALL LEER(2)
      ELOSSP=EPPLI(ISTRA)
      EGAINP=EAPLI(ISTRA)+EMPLI(ISTRA)+EIPLI(ISTRA)+EPHPLI(ISTRA)
      CALL MASAGE ('TOTAL LOSS  ---   TOTAL GAIN                   ')
      CALL MASR2 ('LOSS--GAIN      ',ELOSSP,EGAINP)
      CALL LEER(2)
C
700   CONTINUE
C
C  INTEGRATED ELECTRON SOURCE TERMS
C
      IF (.NOT.TRCBLE) GOTO 800
C
      CALL HEADNG ('PARTICLE FLUX BALANCE (AMP), ELECTRONS',39)
      CALL LEER(1)
      IF (WTOTE(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('BULK ELECTRON PARTICLE FLUX BEING NEUTRALIZED  ')
        CALL MASR1 ('PTOTE=  ',WTOTE(ISTRA))
      ENDIF
C  ATOMS PLASMA INTERACTION
      IF (.NOT.LPAEL) THEN
        IF (LMSPAEL) THEN
        CALL MASAGE ('ELECTRONS BORN BY ATOM PLASMA INTERACTIONS')
        CALL MASAGE ('PAEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (PAELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ELECTRONS BORN BY ATOM PLASMA INTERACTIONS')
        CALL MASR1 ('PAELI=  ',PAELI(ISTRA))
      ENDIF
C  MOLECULE PLASMA INTERACTION
      IF (.NOT.LPMEL) THEN
        IF (LMSPMEL) THEN
        CALL MASAGE ('ELECTRONS BORN BY MOLECULE PLASMA INTERACTIONS')
        CALL MASAGE ('PMEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (PMELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ELECTRONS BORN BY MOLECULE PLASMA INTERACTIONS')
        CALL MASR1 ('PMELI=  ',PMELI(ISTRA))
      ENDIF
C  TEST ION PLASMA INTERACTION
      IF (.NOT.LPIEL) THEN
        IF (LMSPIEL) THEN
        CALL MASAGE ('ELECTRONS BORN BY TEST ION PLASMA INTERACTIONS')
        CALL MASAGE ('PIEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (PIELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ELECTRONS BORN BY TEST ION PLASMA INTERACTIONS')
        CALL MASR1 ('PIELI=  ',PIELI(ISTRA))
      ENDIF
C  PHOTON PLASMA INTERACTION
      IF (.NOT.LPPHEL) THEN
        IF (LMSPPHEL) THEN
        CALL MASAGE ('ELECTRONS BORN BY PHOTON PLASMA INTERACTIONS')
        CALL MASAGE ('PPHEL    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (PPHELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ELECTRONS BORN BY PHOTON PLASMA INTERACTIONS')
        CALL MASR1 ('PPHELI=  ',PPHELI(ISTRA))
      ENDIF
      CALL LEER(2)
      PLOSSE=WTOTE(ISTRA)
      PGAINE=PAELI(ISTRA)+PMELI(ISTRA)+PIELI(ISTRA)+PPHELI(ISTRA)
      CALL MASAGE ('TOTAL LOSS  ---  TOTAL GAIN                    ')
      CALL MASR2 ('LOSS--GAIN      ',PLOSSE,PGAINE)
      CALL LEER(2)
C
      CALL HEADNG ('ENERGY FLUX BALANCE (WATT), ELECTRONS',37)
      CALL LEER(1)
C  TO BE WRITTEN: ETOTE FOR RECOMBINATION
C     IF (ETOTE(ISTRA).NE.0.D0) THEN
C       CALL MASAGE ('BULK ELECTRON ENERGY FLUX BEING LOST           ')
C       CALL MASAGE ('UPON RECOMBINATION                             ')
C       CALL MASR1 ('ETOTE=  ',ETOTE(ISTRA))
C     ENDIF
      IF (.NOT.LEAEL) THEN
        IF (LMSEAEL) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO ATOM PLASMA INTERACTIONS')
        CALL MASAGE ('EAEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EAELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO ATOM PLASMA INTERACTIONS')
        CALL MASR1 ('EAELI=  ',EAELI(ISTRA))
      ENDIF
      IF (.NOT.LEMEL) THEN
        IF (LMSEMEL) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO MOLECULE PLASMA INTERACTIONS')
        CALL MASAGE ('EMEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EMELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO MOLECULE PLASMA INTERACTIONS')
        CALL MASR1 ('EMELI=  ',EMELI(ISTRA))
      ENDIF
      IF (.NOT.LEIEL) THEN
        IF (LMSEIEL) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO TEST ION PLASMA INTERACTIONS')
        CALL MASAGE ('EIEL     SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EIELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO TEST ION PLASMA INTERACTIONS')
        CALL MASR1 ('EIELI=  ',EIELI(ISTRA))
      ENDIF
      IF (.NOT.LEPHEL) THEN
        IF (LMSEPHEL) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO PHOTON PLASMA INTERACTIONS')
        CALL MASAGE ('EPHEL    SWITCHED OFF => MISSING IN BALANCE    ')
        ENDIF
      ELSE IF (EPHELI(ISTRA).NE.0.D0) THEN
        CALL MASAGE ('ENERGY GAIN DUE TO PHOTON PLASMA INTERACTIONS')
        CALL MASR1 ('EPHELI=  ',EPHELI(ISTRA))
      ENDIF
      ELOSSE=EMELI(ISTRA)+EAELI(ISTRA)+EIELI(ISTRA)+EPHELI(ISTRA)
      EGAINE=0.
      CALL MASAGE ('TOTAL LOSS  ---   TOTAL GAIN                   ')
      CALL MASR2 ('LOSS--GAIN      ',ELOSSE,EGAINE)
      CALL LEER(2)
800   CONTINUE
C
      WRITE (iunout,*) 
     .  'SURFACES, AT WHICH TEST PARTICLES FLUXES ARE REDUCED'
      IF (LSPUMP) THEN
        WRITE (iunout,'(1X,A3,1X,A8,A12)') 
     .        'NO.','SPECIES ',' PUMPED FLUX'
        DO J=1,NLIMPS
          JJ=J
          IF (J.GT.NLIM) JJ=-(J-NLIM)
          SPA=0.D0
          DO IS=1,NSPTOT
            IF (SPUMP(IS,J).GT.0.) THEN
              WRITE (iunout,'(1X,I3,1X,A8,1PE12.4)') 
     .               JJ,TEXTS(IS),SPUMP(IS,J)
              SPA=SPA+SPUMP(IS,J)*NPRT(IS)
            ENDIF
          ENDDO
          IF (SPA.GT.0.D0) CALL LEER(1)
        ENDDO
      ELSE IF (LMSSPUMP) THEN
        CALL MASAGE ('SPUMP SWITCHED OFF => MISSING IN BALANCE       ')
      END IF

      IF (SUM(ABS(CEMETERYV)) > 0._DP) THEN
        CALL LEER(2)
        WRITE (iunout,*) ' CEMETERYV != 0 '
        WRITE (iunout,*) ' CHECK SUBROUTINE UPDATE FOR BUGS! '
        CALL LEER(2)
      END IF
C
C   DETAILED OUTPUT OF FLUXES ONTO AND FROM SURFACES
      IF (NSURPR.GT.0) THEN
        CALL OUTFLX('FLUXES AT SURFACES              ',ISTRA)
      ENDIF

      IF (SUM(ABS(CEMETERYS)) > 0._DP) THEN
        CALL LEER(2)
        WRITE (iunout,*) ' CEMETERYS != 0 '
        WRITE (iunout,*) ' CHECK SUBROUTINE UPDATE FOR BUGS! '
        CALL LEER(2)
      END IF
C
      CALL PAGE
C
1000  CONTINUE
6666  FORMAT (3X,1A8,8X,12(A4,2X,A8,3X))
7777  FORMAT (1X,3A8)
      RETURN
      END
C ===== SOURCE: outflx.f
C
C
C 22.10.03; ein falsches write statement rausgenommen (bei reemitted
C           from incident test ions, totals.
C 27.03.04; iliin=-3 option (only net fluxes on transp. surfaces) re-enforced
C           simultaneous changes in escape.f.
C           not active for bulk particle fluxes updated in subr. locate
C           not yet in repository
C 05.05.04:  printout of spectrum: text improved
C
c 12.12.05:  some more lines shifted 2 columns to the right 
c            (still not complete, done until line 467),
c            This is for surface do loop "do 10000"
c 16.12.05:  variance for adds tally:  printout activated
c
c 16.01.06:  bug fix: suma1, suma2, etc... initialized (=0)
c            otherwise problems due to new options for deactivation of tallies
C
      SUBROUTINE OUTFLX(A,ISTRA)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CSPEZ
      USE CTRCEI
      USE CTEXT
      USE CSDVI
      USE CLGIN
      USE COUTAU
      USE CTRIG
c slmod begin - tet res
      USE CTETRA
c slmod end

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTRA
      CHARACTER(32), INTENT(IN) :: A
      REAL(DP) :: SUMA1(0:NATM,0:NSTRA), VARA1(0:NATM,0:NSTRA),
     .          SUMM1(0:NMOL,0:NSTRA), VARM1(0:NMOL,0:NSTRA),
     .          SUMI1(0:NION,0:NSTRA), VARI1(0:NION,0:NSTRA),
     .          SUMPH1(0:NPHOT,0:NSTRA), VARPH1(0:NPHOT,0:NSTRA),
     .          SUMP1(0:NPLS,0:NSTRA), VARP1(0:NPLS,0:NSTRA),
     .          SUMA2(0:NATM,0:NSTRA), VARA2(0:NATM,0:NSTRA),
     .          SUMM2(0:NMOL,0:NSTRA), VARM2(0:NMOL,0:NSTRA),
     .          SUMI2(0:NION,0:NSTRA), VARI2(0:NION,0:NSTRA),
     .          SUMPH2(0:NPHOT,0:NSTRA), VARPH2(0:NPHOT,0:NSTRA),
     .          SUMP2(0:NPLS,0:NSTRA), VARP2(0:NPLS,0:NSTRA),
     .          SUMS(0:NADS,0:NSTRA),  VARS(0:NADS,0:NSTRA),
     .          SUML(0:NALS,0:NSTRA),  VARL(0:NALS,0:NSTRA)
      REAL(DP) :: HELP(NRAD), HELPP(NLMPGS)
      REAL(DP) :: PGAINE, SUM1, DUMMY, SUMMT, SUMMEM, SUMMTI, SUMMEI,
     .          SUMMTA, SUMMTM, SUMMA, SUMMS, SUMMM, SUMMI, SUMMP,
     .          SUMA, SUMM, SUMI, SUMME, SUMMTP, SUMMEP, SUMP, TTTT,
     .          SUMMEA, SUMMEPH, SUMMTPH, SUMMPH, SUMPH, EN
      INTEGER :: NR, NP, NT, MSURFG, J, NCELL, NTOTAL, N1, N2, N3, I,
     .           ITRII, IPLGN, ISPR, ITEXT, ISTS, NFTI, NFTE, NTCO,
     .           NF, K, ITALS, I0, IADS, IALS, IION, IPLS, IATM, N,
     .           IMOL, IPHOT, ISPC, IOUT, ISF, IE, IT
      INTEGER :: IADTYP(0:4)
      LOGICAL :: LGVRA1(0:NATM,0:NSTRA),
     .           LGVRM1(0:NMOL,0:NSTRA),
     .           LGVRI1(0:NION,0:NSTRA),
     .           LGVRPH1(0:NPHOT,0:NSTRA),
     .           LGVRP1(0:NPLS,0:NSTRA)
      LOGICAL :: LGVRA2(0:NATM,0:NSTRA),
     .           LGVRM2(0:NMOL,0:NSTRA),
     .           LGVRI2(0:NION,0:NSTRA),
     .           LGVRPH2(0:NPHOT,0:NSTRA),
     .           LGVRP2(0:NPLS,0:NSTRA)
      LOGICAL :: LGVARS(0:NADS,0:NSTRA),
     .           LGVARL(0:NALS,0:NSTRA)
      LOGICAL :: LOGADS(0:NADS,0:NSTRA),
     .           LOGALS(0:NALS,0:NSTRA)
      LOGICAL :: PRINTED(NLIMPS)
      CHARACTER(8) :: TEXTA(NADS), TEXTL(NALS)
      CHARACTER(24) :: TXT24
      CHARACTER(72) :: TXT72
      CHARACTER(10) :: TEXTYP(0:4)

      PRINTED = .FALSE.
      CALL LEER(1)
      WRITE (iunout,9999) A
C
C  SURFACE LOOP
C
      DO 10000 ISPR=1,NSURPR
C
C
        ITEXT=0
        I=NPRSRF(ISPR)
        CALL LEER(2)
        IF (IGJUM0(I).NE.0) THEN
C  CHECK, IF SURFACE "I" IS STILL THERE AS ONE SIDE OF A TRIANGLE
          IF (LEVGEO.EQ.4) THEN
            DO ITRII=1,NTRII
            DO IPLGN=1,3
              ISTS=ABS(INMTI(IPLGN,ITRII))
              IF (ISTS.EQ.I) GOTO 5
            ENDDO
            ENDDO
          ENDIF
          WRITE (iunout,*) ' SURFACE NO. ',I,' : OUT '
          PRINTED(I) = .TRUE.
          GOTO 10000
        ENDIF
5       CONTINUE
C  PRINT SURFACE AREA (NOT FOR "TIME SURFACE")
        CALL MASBOX(TXTSFL(I))
        IF (NTIME.GE.1.AND.I.EQ.NLIM+NSTSI) GOTO 1
        IF (SAREA(I).NE.666.) THEN
          WRITE (iunout,'(A22,1P,1E12.4)') ' SURFACE AREA (CM**2) ',
     .           SAREA(I)
        ELSE
          WRITE (iunout,*) 'SURFACE AREA (CM**2) ','?'
        ENDIF
1       CONTINUE
C
C
c slmod begin
        IF (I.GT.NLIM.AND.LEVGEO.LE.5.AND.NLMPGS.NE.NLIMPS) THEN
c
c        IF (I.GT.NLIM.AND.LEVGEO.LE.4.AND.NLMPGS.NE.NLIMPS) THEN
c slmod end
C
C  SPATIAL RESOLUTION ON NON DEFAULT STANDARD SURFACE?
          HELP=0.D0
          ISTS=I-NLIM
C
          ITALS=NPRTLS(ISPR)
          IF (ITALS.LE.0.OR.ITALS.GT.NTALS) GOTO 11
          I0=0
          IF (NFRTWI(ITALS).GT.1) I0=1
          NFTI=1
          NFTE=NFSTWI(ITALS)
          IF (NSPEZS(ISPR,1).GT.0) THEN
            NFTI=NSPEZS(ISPR,1)
            NFTE=MAX(NFTI,NSPEZS(ISPR,2))
          ENDIF
          NF=NFRSTW(ITALS)
          DO 10 K=NFTI,NFTE
            IF (K.GT.NFSTWI(ITALS)) THEN
              CALL LEER(1)
              WRITE (iunout,*) 'SPECIES INDEX OUT OF RANGE IN OUTFLX'
              WRITE (iunout,*) 'ISTS,ITALS, K, ',ISTS,ITALS,K
              CALL LEER(1)
              GOTO 10
            ENDIF
C
            DO J=1,NLMPGS
              HELPP(J)=ESTIMS(NADDW(ITALS)+K,J)
            END DO
C
            IF (LEVGEO.LE.3) THEN
              IF (INUMP(ISTS,2).NE.0) then
C  POLOIDAL SURFACE
                sum1=0
                np=1
                do nr=1,nr1st
                do nt=1,nt3rd
                  MSURFG=NR+(NT-1)*NR1P2
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=NR1ST
                N2=1
                N3=NT3RD
              ELSEIF (INUMP(ISTS,1).NE.0) then
C  RADIAL SURFACE
                sum1=0
                NR=1
                do np=1,np2nd
                do nt=1,nt3rd
                  MSURFG=NP+(NT-1)*NP2T3
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=1
                N2=NP2ND
                N3=NT3RD
              ELSEIF (INUMP(ISTS,3).NE.0) then
C  TOROIDAL SURFACE
                sum1=0
                nt=1
                do nr=1,nr1st
                do np=1,np2nd
                  MSURFG=Nr+(Np-1)*Nr1p2
                  MSURFG=NLIM+NSTS+MSURFG+(ISTS-1)*NGITT
                  NCELL=NR+((NP-1)+(NT-1)*NP2T3)*NR1P2
                  HELP(ncell)=HELPP(MSURFG)
                  sum1=sum1+HELPP(msurfg)
                ENDDO
                ENDDO
                N1=NR1ST
                N2=NP2ND
                N3=1
              ENDIF
            ELSE IF (LEVGEO.EQ.4) THEN
              sum1=0.D0
              ntco=0
              DO NP=1,3
                DO NR=1,NTRII
                  IF (INMTI(NP,NR) == NLIM+ISTS) THEN
                    MSURFG=NLIM+NSTS+INSPAT(NP,NR)
                    NTCO=NTCO+1
                    HELP(NTCO)=HELPP(MSURFG)
                    SUM1=SUM1+HELPP(MSURFG)
                  END IF
                END DO
              END DO
              N1=NTCO+1
              N2=1
              N3=1
c slmod begin - tet res
            ELSE IF (LEVGEO.EQ.5) THEN
              sum1=0.D0
              ntco=0
              DO NP=1,4
                DO NR=1,NTET
                  IF (INMTIT(NP,NR) == NLIM+ISTS) THEN
                    MSURFG=NLIM+NSTS+INSPATT(NP,NR)
                    NTCO=NTCO+1
                    HELP(NTCO)=HELPP(MSURFG)
                    SUM1=SUM1+HELPP(MSURFG)
                  END IF
                END DO
              END DO
              N1=NTCO+1
              N2=1
              N3=1
c slmod end
            END IF
            write (iunout,*) 'test ',sum1
            NTOTAL=N1*N2*N3
            CALL INTVOL (HELP,1,1,NTOTAL,DUMMY,N1,N2,N3,1)
            IF (ABS(DUMMY) > EPS60) THEN
              CALL PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,NFLAGS(ISPR),
     .                  NTLSFL(ISPR),
     .                  IRPTA(ISTS,1),IRPTE(ISTS,1),IRPTA(ISTS,2),
     .                  IRPTE(ISTS,2),1,1)
            ELSE
              CALL PRTTLS(TXTTLW(K,ITALS),TXTSPW(K,ITALS),
     .                  TXTUNW(K,ITALS),
     .                  HELP,N1,N2,N3,1,NTOTAL,-1,NTLSFL(ISPR),
     .                  IRPTA(ISTS,1),IRPTE(ISTS,1),IRPTA(ISTS,2),
     .                  IRPTE(ISTS,2),1,1)
              CALL MASAGE
     .            ('IDENTICAL ZERO, NOT PRINTED                  ')
              CALL LEER(2)
            END IF
10        CONTINUE
11        CONTINUE
        ENDIF


C  SPECTRA
        IF (NTLSFL(ISPR) > 0) THEN
          ISF=I
          IF (I < 0) ISF=ABS(I)+NLIM
          IOUT = NTLSFL(ISPR)

          TEXTYP(0) = 'PHOTONS   '
          TEXTYP(1) = 'ATOMS     '
          TEXTYP(2) = 'MOLECULES '
          TEXTYP(3) = 'TEST IONS '
          TEXTYP(4) = 'BULK IONS '
          IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

          DO ISPC=1,NADSPC
            IF ((ESTIML(ISPC)%PSPC%ISRFCLL == 0) .AND.
     .          (ESTIML(ISPC)%PSPC%ISPCSRF == ISF)) THEN
              WRITE (IOUT,*)
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' SPECTRUM CALCULATED FOR SURFACE ',I
              IT = ESTIML(ISPC)%PSPC%ISPCTYP
              IF (IT == 1) THEN
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'PARTICLE FLUX IN AMP'
              ELSE
                WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                            'ENERGY FLUX IN WATT '
              END IF
              WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .              TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
              IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
                WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
              ELSE
                WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .                 TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .                       ESTIML(ISPC)%PSPC%IPRSP)
              END IF
              WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .               ESTIML(ISPC)%PSPC%SPCMIN
              WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .               ESTIML(ISPC)%PSPC%SPCMAX
              WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ',
     .               ESTIML(ISPC)%PSPC%NSPC
              WRITE (IOUT,*)
              IF (ESTIML(ISPC)%PSPC%SPCINT > EPS60) THEN
                IF (NSIGI_SPC == 0) THEN
                  DO IE=1, ESTIML(ISPC)%PSPC%NSPC
                    EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .                   (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
                    WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .                 ESTIML(ISPC)%PSPC%SPC(IE)
                  END DO
                ELSE
                  DO IE=1, ESTIML(ISPC)%PSPC%NSPC
                    EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .                   (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
                    WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
     .                   ESTIML(ISPC)%PSPC%SPC(IE),
     .                   ESTIML(ISPC)%PSPC%SDV(IE)
                  END DO
                END IF
              ELSE
                WRITE (IOUT,*) ' SPECTRUM IDENTICAL 0 '
              END IF
              WRITE (IOUT,*)
              WRITE (IOUT,*) ' INTEGRAL OF SPECTRUM ',
     .               ESTIML(ISPC)%PSPC%SPCINT
              IF (NSIGI_SPC > 0)
     .          WRITE (IOUT,*) ' STANDARD DEVIATION  ',
     .               ESTIML(ISPC)%PSPC%SGMS
            END IF
          END DO
        END IF

        IF (PRINTED(I)) CYCLE
C
C  *****************************************************
C   INCIDENT FLUXES, POSITIVE PARTIAL FLUXES, NET FLUXES
C  *****************************************************
C
C
C   SURFACE AVERAGED TALLY NO.1 AND NO.26
C
        SUMMT=0.
        SUMME=0.
        SUMMTP=0.
        SUMMEP=0.
        SUMA=0.
        SUMM=0.
        SUMI=0.
        SUMP=0.
        SUMPH=0.
        SUMA1(:,ISTRA)=0._DP
        SUMA2(:,ISTRA)=0._DP
        DO 20 IATM=1,NATMI
          LGVRA1(IATM,ISTRA)=.FALSE.
          LGVRA2(IATM,ISTRA)=.FALSE.
          IF (LPOTAT) SUMA1(IATM,ISTRA)=POTAT(IATM,I)
          SUMMT=SUMMT+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
          SUMA=SUMA+SUMA1(IATM,ISTRA)
          IF (LEOTAT) SUMA2(IATM,ISTRA)=EOTAT(IATM,I)
          SUMME=SUMME+SUMA2(IATM,ISTRA)
20      CONTINUE
C
        DO 21 N=1,NSIGSI
          IF (IIHW(N).EQ.1) THEN
            DO 22 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 22
              VARA1(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
22          CONTINUE
          ELSEIF (IIHW(N).EQ.26) THEN
            DO 522 IATM=1,NATMI
              IF (IGHW(N).NE.IATM) GOTO 522
              VARA2(IATM,ISTRA)=SIGMAW(N,I)
              LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
522         CONTINUE
          ENDIF
21      CONTINUE
C
C   SURFACE AVERAGED TALLY NO.7 AND NO.32
C
        SUMM1(:,ISTRA)=0._DP
        SUMM2(:,ISTRA)=0._DP
        DO 30 IMOL=1,NMOLI
          LGVRM1(IMOL,ISTRA)=.FALSE.
          LGVRM2(IMOL,ISTRA)=.FALSE.
          IF (LPOTML) SUMM1(IMOL,ISTRA)=POTML(IMOL,I)
          SUMMT=SUMMT+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
          SUMM=SUMM+SUMM1(IMOL,ISTRA)
          IF (LEOTML) SUMM2(IMOL,ISTRA)=EOTML(IMOL,I)
          SUMME=SUMME+SUMM2(IMOL,ISTRA)
30      CONTINUE
C
        DO 23 N=1,NSIGSI
          IF (IIHW(N).EQ.7) THEN
            DO 24 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 24
              VARM1(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
24          CONTINUE
          ELSEIF (IIHW(N).EQ.32) THEN
            DO 524 IMOL=1,NMOLI
              IF (IGHW(N).NE.IMOL) GOTO 524
              VARM2(IMOL,ISTRA)=SIGMAW(N,I)
              LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
524         CONTINUE
          ENDIF
23      CONTINUE
C
C   SURFACE AVERAGED TALLY NO.13 AND NO.38
C
        SUMI1(:,ISTRA)=0._DP
        SUMI2(:,ISTRA)=0._DP
        DO 40 IION=1,NIONI
          LGVRI1(IION,ISTRA)=.FALSE.
          LGVRI2(IION,ISTRA)=.FALSE.
          IF (LPOTIO) SUMI1(IION,ISTRA)=POTIO(IION,I)
          SUMMT=SUMMT+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
          SUMI=SUMI+SUMI1(IION,ISTRA)
          IF (LEOTIO) SUMI2(IION,ISTRA)=EOTIO(IION,I)
          SUMME=SUMME+SUMI2(IION,ISTRA)
40      CONTINUE
C
        DO 25 N=1,NSIGSI
          IF (IIHW(N).EQ.13) THEN
            DO 26 IION=1,NIONI
              IF (IGHW(N).NE.IION) GOTO 26
              VARI1(IION,ISTRA)=SIGMAW(N,I)
              LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
26          CONTINUE
          ELSEIF (IIHW(N).EQ.38) THEN
            DO 526 IION=1,NIONI
              IF (IGHW(N).NE.IION) GOTO 526
              VARI2(IION,ISTRA)=SIGMAW(N,I)
              LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
526         CONTINUE
          ENDIF
25      CONTINUE
C
C   SURFACE AVERAGED TALLY NO 19 and NO 44
C
        SUMPH1(:,ISTRA)=0._DP
        SUMPH2(:,ISTRA)=0._DP
        DO IPHOT=1,NPHOTI
          LGVRPH1(IPHOT,ISTRA)=.FALSE.
          LGVRPH2(IPHOT,ISTRA)=.FALSE.
          IF (LPOTPHT) SUMPH1(IPHOT,ISTRA)=POTPHT(IPHOT,I)
          SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
          SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
          IF (LEOTPHT) SUMPH2(IPHOT,ISTRA)=EOTPHT(IPHOT,I)
          SUMME=SUMME+SUMPH2(IPHOT,ISTRA)
        ENDDO
  
        DO N=1,NSIGSI
          IF (IIHW(N).EQ.19) THEN
            DO IPHOT=1,NPHOTI
              IF (IGHW(N).NE.IPHOT) CYCLE
              VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
              LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
            ENDDO
          ELSEIF (IIHW(N).EQ.44) THEN
            DO IPHOT=1,NPHOTI
              IF (IGHW(N).NE.IPHOT) CYCLE
              VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
              LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
            ENDDO
          ENDIF
        ENDDO
C
C   SURFACE AVERAGED TALLY NO.25 AND NO.50
C
        SUMP1(:,ISTRA)=0._DP
        SUMP2(:,ISTRA)=0._DP
        DO 50 IPLS=1,NPLSI
         LGVRP1(IPLS,ISTRA)=.FALSE.
         LGVRP2(IPLS,ISTRA)=.FALSE.
         IF (LPOTPL) SUMP1(IPLS,ISTRA)=POTPL(IPLS,I)
         SUMMTP=SUMMTP+SUMP1(IPLS,ISTRA)*NPRT(NSPAMI+IPLS)
         SUMP=SUMP+SUMP1(IPLS,ISTRA)
         IF (LEOTPL) SUMP2(IPLS,ISTRA)=EOTPL(IPLS,I)
         SUMMEP=SUMMEP+SUMP2(IPLS,ISTRA)
50      CONTINUE
C
        DO 27 N=1,NSIGSI
          IF (IIHW(N).EQ.25) THEN
            DO 28 IPLS=1,NPLSI
              IF (IGHW(N).NE.IPLS) GOTO 28
              VARP1(IPLS,ISTRA)=SIGMAW(N,I)
              LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
28          CONTINUE
          ELSEIF (IIHW(N).EQ.50) THEN
            DO 528 IPLS=1,NPLSI
              IF (IGHW(N).NE.IPLS) GOTO 528
              VARP2(IPLS,ISTRA)=SIGMAW(N,I)
              LGVRP2(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
528         CONTINUE
          ENDIF
27      CONTINUE
C
C
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMP)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO FLUXES INCIDENT ON THIS SURFACE             ')
        CALL LEER(1)
      ELSE
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (iunout,*) 'FLUX INCIDENT ON SURFACE:'
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (iunout,*) 'INCIDENT: ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (iunout,*) 'INCIDENT: MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (iunout,*) 'INCIDENT: TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (iunout,*) 'INCIDENT: PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (iunout,*) 'INCIDENT: BULK IONS (RECYCLING SOURCE)'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (iunout,*) 
     .      'TOTAL INCIDENT "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (iunout,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('TOT.PFLX',SUMMT)
          CALL MASR1 ('TOT.EFLX',SUMME)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          IF (SUMA.NE.0.D0.OR.SUMM.NE.0.D0.OR.SUMI.NE.0.D0
     .                    .OR.SUMPH.NE.0.D0)
     .    WRITE (iunout,*) 
     .      'PARTIAL PARTICLE AND ENERGY CURRENTS, POSITIVE '
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (iunout,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (iunout,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (iunout,*) 'TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (iunout,*) 'PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (iunout,*) 
     .        'FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (iunout,*) 'BULK IONS'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (iunout,*) 
     .      'TOTAL POSITIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (iunout,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('POS.PFLX',SUMMT)
          CALL MASR1 ('POS.EFLX',SUMME)
C
        ELSEIF (ILIIN(I).EQ.-3) THEN
          IF (SUMA.NE.0.D0.OR.SUMM.NE.0.D0.OR.SUMI.NE.0.D0
     .                    .OR.SUMPH.NE.0.D0)
     .    WRITE (iunout,*) 'NET PARTICLE AND ENERGY CURRENTS '
C  SURFACE AVERAGED TALLY NO. 1
          IF (SUMA.NE.0.D0) THEN
            WRITE (iunout,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 26
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 7
          IF (SUMM.NE.0.D0) THEN
            WRITE (iunout,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 32
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 13
          IF (SUMI.NE.0.D0) THEN
            WRITE (iunout,*) 'TEST IONS'
            CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 38
            CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
            CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .                   TEXTS(NSPAM+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 19
          IF (SUMPH.NE.0.D0) THEN
            WRITE (iunout,*) 'PHOTONS'
            CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 44
            CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
            CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,
     .                   NSTRA,TEXTS(0+1))
          ENDIF
C  SURFACE AVERAGED TALLY NO. 25
          IF (SUMP.NE.0.D0) THEN
            WRITE (iunout,*) 
     .        'FLUX INCIDENT ON SURFACE (RECYCLING SOURCE):'
            WRITE (iunout,*) 'BULK IONS'
            CALL MASYR1('P-FLUX:  ',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
C  SURFACE AVERAGED TALLY NO. 50
            CALL MASYR1('E-FLUX:  ',SUMP2,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
            CALL MASYR1('ST.DEV.% ',VARP2,LGVRP2,ISTRA,0,NPLS,0,NSTRA,
     .                   TEXTS(NSPAMI+1))
          ENDIF
          CALL LEER (1)
          WRITE (iunout,*) 'TOTAL NET "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0)
     .    WRITE (iunout,*) '(EXCLUDING BULK IONS (RECYCLING SOURCE) '
          CALL MASR1 ('NET PFLX',SUMMT)
          CALL MASR1 ('NET EFLX',SUMME)


        ENDIF
C
C  INDEPENDENT OF VALUE AND SIGN OF ILIIN:
C
        IF (SUMMTP.GT.0.D0) THEN
          CALL LEER (1)
          WRITE (iunout,*) 
     .      'TOTAL INCIDENT RECYCLING SOURCE "ATOMIC" FLUXES'
          WRITE (iunout,*) 'BULK IONS'
          CALL MASR1 ('SRC.PFLX',SUMMTP)
          CALL MASR1 ('SRC.EFLX',SUMMEP)
        ENDIF
      ENDIF
C
C  ******************************************
C   REEMITTED FLUXES, NEGATIVE PARTIAL FLUXES
C  ******************************************
C
C   FIRST: FROM INCIDENT ATOMS
C
      SUMMTA=0.
      SUMMEA=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.2 AND NO.27
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFAAT) SUMA1(IATM,ISTRA)=PRFAAT(IATM,I)
        SUMMTA=SUMMTA+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFAAT) SUMA2(IATM,ISTRA)=ERFAAT(IATM,I)
        SUMMEA=SUMMEA+SUMA2(IATM,ISTRA)
102   CONTINUE
C
      DO 121 N=1,NSIGSI
        IF (IIHW(N).EQ.2) THEN
          DO 122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
122       CONTINUE
        ELSEIF (IIHW(N).EQ.27) THEN
          DO 622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
622       CONTINUE
        ENDIF
121   CONTINUE
C
C
C   SURFACE AVERAGED TALLY NO.8 AND NO.33
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFAML) SUMM1(IMOL,ISTRA)=PRFAML(IMOL,I)
        SUMMTA=SUMMTA+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFAML) SUMM2(IMOL,ISTRA)=ERFAML(IMOL,I)
        SUMMEA=SUMMEA+SUMM2(IMOL,ISTRA)
103   CONTINUE
C
      DO 123 N=1,NSIGSI
        IF (IIHW(N).EQ.8) THEN
          DO 124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
124       CONTINUE
        ELSEIF (IIHW(N).EQ.33) THEN
          DO 624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
624       CONTINUE
        ENDIF
123   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.14 AND NO.39
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFAIO) SUMI1(IION,ISTRA)=PRFAIO(IION,I)
        SUMMTA=SUMMTA+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFAIO) SUMI2(IION,ISTRA)=ERFAIO(IION,I)
        SUMMEA=SUMMEA+SUMI2(IION,ISTRA)
104   CONTINUE
C
      DO 125 N=1,NSIGSI
        IF (IIHW(N).EQ.14) THEN
          DO 126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
126       CONTINUE
        ELSEIF (IIHW(N).EQ.39) THEN
          DO 626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
626       CONTINUE
        ENDIF
125   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.20 AND NO.45
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFAPHT) SUMPH1(IPHOT,ISTRA)=PRFAPHT(IPHOT,I)
        SUMMTA=SUMMTA+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFAPHT) SUMPH2(IPHOT,ISTRA)=ERFAPHT(IPHOT,I)
        SUMMEA=SUMMEA+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.20) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.45) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT ATOMS '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT ATOMS:'
C  SURFACE AVERAGED TALLY NO. 2
        IF (SUMA.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 27
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 8
        IF (SUMM.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 33
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 14
        IF (SUMI.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 39
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 20
        IF (SUMPH.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 45
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTA)
        CALL MASR1 ('TOT.EFLX',SUMMEA)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 2
          IF (SUMA.NE.0.D0) THEN
            WRITE (iunout,*) 
     .        'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (iunout,*) 'ATOMS'
            CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 27
            CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .                   TEXTS(NSPH+1))
            CALL LEER (1)
            CALL MASR1 ('NEG.PFLX',SUMMTA)
            CALL MASR1 ('NEG.EFLX',SUMMEA)
          ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT MOLECULES
      SUMMTM=0.
      SUMMEM=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.3 AND NO.28
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 1102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFMAT) SUMA1(IATM,ISTRA)=PRFMAT(IATM,I)
        SUMMTM=SUMMTM+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFMAT) SUMA2(IATM,ISTRA)=ERFMAT(IATM,I)
        SUMMEM=SUMMEM+SUMA2(IATM,ISTRA)
1102  CONTINUE
C
      DO 1121 N=1,NSIGSI
        IF (IIHW(N).EQ.3) THEN
          DO 1122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
1122      CONTINUE
        ELSEIF (IIHW(N).EQ.28) THEN
          DO 1622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 1622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
1622      CONTINUE
        ENDIF
1121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.9 AND NO.34
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 1103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFMML) SUMM1(IMOL,ISTRA)=PRFMML(IMOL,I)
        SUMMTM=SUMMTM+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFMML) SUMM2(IMOL,ISTRA)=ERFMML(IMOL,I)
        SUMMEM=SUMMEM+SUMM2(IMOL,ISTRA)
1103  CONTINUE
C
      DO 1123 N=1,NSIGSI
        IF (IIHW(N).EQ.9) THEN
          DO 1124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
1124      CONTINUE
        ELSEIF (IIHW(N).EQ.34) THEN
          DO 1624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 1624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
1624      CONTINUE
        ENDIF
1123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.15 AND NO.40
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 1104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFMIO) SUMI1(IION,ISTRA)=PRFMIO(IION,I)
        SUMMTM=SUMMTM+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFMIO) SUMI2(IION,ISTRA)=ERFMIO(IION,I)
        SUMMEM=SUMMEM+SUMI2(IION,ISTRA)
1104  CONTINUE
C
      DO 1125 N=1,NSIGSI
        IF (IIHW(N).EQ.15) THEN
          DO 1126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
1126      CONTINUE
        ELSEIF (IIHW(N).EQ.40) THEN
          DO 1626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 1626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
1626      CONTINUE
        ENDIF
1125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.21 AND NO.46
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFMPHT) SUMPH1(IPHOT,ISTRA)=PRFMPHT(IPHOT,I)
        SUMMTM=SUMMTM+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFMPHT) SUMPH2(IPHOT,ISTRA)=ERFMPHT(IPHOT,I)
        SUMMEM=SUMMEM+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.21) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.46) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT MOLECULES '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT MOLECULES:'
C  SURFACE AVERAGED TALLY NO. 3
        IF (SUMA.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 28
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 9
        IF (SUMM.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 34
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 15
        IF (SUMI.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 40
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 21
        IF (SUMPH.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 46
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTM)
        CALL MASR1 ('TOT.EFLX',SUMMEM)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 9
          IF (SUMM.NE.0.D0) THEN
            IF (ITEXT.EQ.0) WRITE (iunout,*) 
     .         'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE'
            ITEXT=1
            WRITE (iunout,*) 'MOLECULES'
            CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 34
            CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .                   TEXTS(NSPA+1))
            CALL LEER (1)
            CALL MASR1 ('NEG.PFLX',SUMMTM)
            CALL MASR1 ('NEG.EFLX',SUMMEM)
          ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT TEST-IONS
      SUMMTI=0.
      SUMMEI=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.4 AND NO.29
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 2102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFIAT) SUMA1(IATM,ISTRA)=PRFIAT(IATM,I)
        SUMMTI=SUMMTI+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFIAT) SUMA2(IATM,ISTRA)=ERFIAT(IATM,I)
        SUMMEI=SUMMEI+SUMA2(IATM,ISTRA)
2102  CONTINUE
C
      DO 2121 N=1,NSIGSI
        IF (IIHW(N).EQ.4) THEN
          DO 2122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
2122      CONTINUE
        ELSEIF (IIHW(N).EQ.29) THEN
          DO 2622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 2622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
2622      CONTINUE
        ENDIF
2121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.10 AND NO.35
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 2103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFIML) SUMM1(IMOL,ISTRA)=PRFIML(IMOL,I)
        SUMMTI=SUMMTI+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFIML) SUMM2(IMOL,ISTRA)=ERFIML(IMOL,I)
        SUMMEI=SUMMEI+SUMM2(IMOL,ISTRA)
2103  CONTINUE
C
      DO 2123 N=1,NSIGSI
        IF (IIHW(N).EQ.10) THEN
          DO 2124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
2124      CONTINUE
        ELSEIF (IIHW(N).EQ.35) THEN
          DO 2624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 2624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
2624      CONTINUE
        ENDIF
2123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.16 AND NO.41
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 2104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFIIO) SUMI1(IION,ISTRA)=PRFIIO(IION,I)
        SUMMTI=SUMMTI+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFIIO) SUMI2(IION,ISTRA)=ERFIIO(IION,I)
        SUMMEI=SUMMEI+SUMI2(IION,ISTRA)
2104  CONTINUE
C
      DO 2125 N=1,NSIGSI
        IF (IIHW(N).EQ.16) THEN
          DO 2126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
2126      CONTINUE
        ELSEIF (IIHW(N).EQ.41) THEN
          DO 2626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 2626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
2626      CONTINUE
        ENDIF
2125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO 22 AND NO.47
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFIPHT) SUMPH1(IPHOT,ISTRA)=PRFIPHT(IPHOT,I)
        SUMMTI=SUMMTI+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFIPHT) SUMPH2(IPHOT,ISTRA)=ERFIPHT(IPHOT,I)
        SUMMEI=SUMMEI+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.22) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.47) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT TEST IONS '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT TEST IONS:'
C  SURFACE AVERAGED TALLY NO. 4
        IF (SUMA.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 29
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 10
        IF (SUMM.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 35
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 16
        IF (SUMI.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 41
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 22
        IF (SUMPH.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 47
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTI)
        CALL MASR1 ('TOT.EFLX',SUMMEI)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 16
        IF (SUMI.NE.0.D0) THEN
        IF (ITEXT.EQ.0) WRITE (iunout,*) 
     .    'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE '
        ITEXT=1
        WRITE (iunout,*) 'TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 41
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL LEER (1)
        CALL MASR1 ('NEG.PFLX',SUMMTI)
        CALL MASR1 ('NEG.EFLX',SUMMEI)
        ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT PHOTONS
      SUMMTPH=0.
      SUMMEPH=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.5 AND NO.30
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 4102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFPHAT) SUMA1(IATM,ISTRA)=PRFPHAT(IATM,I)
        SUMMTPH=SUMMTPH+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFPHAT) SUMA2(IATM,ISTRA)=ERFPHAT(IATM,I)
        SUMMEPH=SUMMEPH+SUMA2(IATM,ISTRA)
4102  CONTINUE
C
      DO 4121 N=1,NSIGSI
        IF (IIHW(N).EQ.5) THEN
          DO 4122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
4122      CONTINUE
        ELSEIF (IIHW(N).EQ.30) THEN
          DO 4622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 4622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
4622      CONTINUE
        ENDIF
4121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.11 AND NO.36
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 4103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFPHML) SUMM1(IMOL,ISTRA)=PRFPHML(IMOL,I)
        SUMMTPH=SUMMTPH+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFPHML) SUMM2(IMOL,ISTRA)=ERFPHML(IMOL,I)
        SUMMEPH=SUMMEPH+SUMM2(IMOL,ISTRA)
4103  CONTINUE
C
      DO 4123 N=1,NSIGSI
        IF (IIHW(N).EQ.11) THEN
          DO 4124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
4124      CONTINUE
        ELSEIF (IIHW(N).EQ.36) THEN
          DO 4624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 4624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
4624      CONTINUE
        ENDIF
4123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.17 AND NO.42
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 4104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFPHIO) SUMI1(IION,ISTRA)=PRFPHIO(IION,I)
        SUMMTPH=SUMMTPH+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFPHIO) SUMI2(IION,ISTRA)=ERFPHIO(IION,I)
        SUMMEPH=SUMMEPH+SUMI2(IION,ISTRA)
4104   CONTINUE
C
      DO 4125 N=1,NSIGSI
        IF (IIHW(N).EQ.17) THEN
          DO 4126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
4126      CONTINUE
        ELSEIF (IIHW(N).EQ.42) THEN
          DO 4626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 4626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
4626      CONTINUE
        ENDIF
4125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO 23 AND NO.48
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFPHPHT) SUMPH1(IPHOT,ISTRA)=PRFPHPHT(IPHOT,I)
        SUMMTI=SUMMTI+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFPHPHT) SUMPH2(IPHOT,ISTRA)=ERFPHPHT(IPHOT,I)
        SUMMEI=SUMMEI+SUMPH2(IPHOT,ISTRA)
      ENDDO

      DO N=1,NSIGSI
        IF (IIHW(N).EQ.23) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.48) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0.AND.ILIIN(I).GT.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'NO FLUXES REEMITTED FROM INCIDENT PHOTONS   '
        CALL LEER(1)
      ELSEIF (ILIIN(I).NE.-3) THEN
        CALL LEER(1)
C
        IF (ILIIN(I).GT.0) THEN
        WRITE (iunout,*) 'FLUX REEMITTED FROM INCIDENT PHOTONS:'
C  SURFACE AVERAGED TALLY NO. 5
        IF (SUMA.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 30
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 11
        IF (SUMM.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 36
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 17
        IF (SUMI.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 42
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 48
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        CALL MASR1 ('TOT.PFLX',SUMMTPH)
        CALL MASR1 ('TOT.EFLX',SUMMEPH)
C
        ELSEIF (ILIIN(I).LT.0) THEN
C  SURFACE AVERAGED TALLY NO. 23
        IF (SUMPH.NE.0.D0) THEN
        IF (ITEXT.EQ.0) WRITE (iunout,*) 
     .    'PARTIAL PARTICLE AND ENERGY CURRENTS, NEGATIVE '
        ITEXT=1
        WRITE (iunout,*) 'PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 48
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL LEER (1)
        CALL MASR1 ('NEG.PFLX',SUMMTPH)
        CALL MASR1 ('NEG.EFLX',SUMMEPH)
        ENDIF
        ENDIF
      ENDIF
C
C
C   REEMITTED FLUXES, NEXT: FROM INCIDENT BULK-IONS
      SUMMTP=0.
      SUMMEP=0.
      SUMA=0.
      SUMM=0.
      SUMI=0.
      SUMPH=0.
C
C   SURFACE AVERAGED TALLY NO.6 AND NO.31
C
      SUMA1(:,ISTRA)=0._DP
      SUMA2(:,ISTRA)=0._DP
      DO 3102 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        LGVRA2(IATM,ISTRA)=.FALSE.
        IF (LPRFPAT) SUMA1(IATM,ISTRA)=PRFPAT(IATM,I)
        SUMMTP=SUMMTP+SUMA1(IATM,ISTRA)*NPRT(NSPH+IATM)
        SUMA=SUMA+SUMA1(IATM,ISTRA)
        IF (LERFPAT) SUMA2(IATM,ISTRA)=ERFPAT(IATM,I)
        SUMMEP=SUMMEP+SUMA2(IATM,ISTRA)
3102  CONTINUE
C
      DO 3121 N=1,NSIGSI
        IF (IIHW(N).EQ.6) THEN
          DO 3122 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3122
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
3122      CONTINUE
        ELSEIF (IIHW(N).EQ.31) THEN
          DO 3622 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 3622
            VARA2(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA2(IATM,ISTRA)=LOGATM(IATM,ISTRA)
3622      CONTINUE
        ENDIF
3121  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.12 AND NO.37
C
      SUMM1(:,ISTRA)=0._DP
      SUMM2(:,ISTRA)=0._DP
      DO 3103 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        LGVRM2(IMOL,ISTRA)=.FALSE.
        IF (LPRFPML) SUMM1(IMOL,ISTRA)=PRFPML(IMOL,I)
        SUMMTP=SUMMTP+SUMM1(IMOL,ISTRA)*NPRT(NSPA+IMOL)
        SUMM=SUMM+SUMM1(IMOL,ISTRA)
        IF (LERFPML) SUMM2(IMOL,ISTRA)=ERFPML(IMOL,I)
        SUMMEP=SUMMEP+SUMM2(IMOL,ISTRA)
3103  CONTINUE
C
      DO 3123 N=1,NSIGSI
        IF (IIHW(N).EQ.12) THEN
          DO 3124 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3124
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
3124      CONTINUE
        ELSEIF (IIHW(N).EQ.37) THEN
          DO 3624 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 3624
            VARM2(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM2(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
3624      CONTINUE
        ENDIF
3123  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.18 AND NO.43
C
      SUMI1(:,ISTRA)=0._DP
      SUMI2(:,ISTRA)=0._DP
      DO 3104 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        LGVRI2(IION,ISTRA)=.FALSE.
        IF (LPRFPIO) SUMI1(IION,ISTRA)=PRFPIO(IION,I)
        SUMMTP=SUMMTP+SUMI1(IION,ISTRA)*NPRT(NSPAM+IION)
        SUMI=SUMI+SUMI1(IION,ISTRA)
        IF (LERFPIO) SUMI2(IION,ISTRA)=ERFPIO(IION,I)
        SUMMEP=SUMMEP+SUMI2(IION,ISTRA)
3104  CONTINUE
C
      DO 3125 N=1,NSIGSI
        IF (IIHW(N).EQ.18) THEN
          DO 3126 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3126
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
3126      CONTINUE
        ELSEIF (IIHW(N).EQ.43) THEN
          DO 3626 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 3626
            VARI2(IION,ISTRA)=SIGMAW(N,I)
            LGVRI2(IION,ISTRA)=LOGION(IION,ISTRA)
3626      CONTINUE
        ENDIF
3125  CONTINUE
C
C   SURFACE AVERAGED TALLY NO.24 AND NO.49
C
      SUMPH1(:,ISTRA)=0._DP
      SUMPH2(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        LGVRPH2(IPHOT,ISTRA)=.FALSE.
        IF (LPRFPPHT) SUMPH1(IPHOT,ISTRA)=PRFPPHT(IPHOT,I)
        SUMMTP=SUMMTP+SUMPH1(IPHOT,ISTRA)*NPRT(0+IPHOT)
        SUMPH=SUMPH+SUMPH1(IPHOT,ISTRA)
        IF (LERFPPHT) SUMPH2(IPHOT,ISTRA)=ERFPPHT(IPHOT,I)
        SUMMEP=SUMMEP+SUMPH2(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.24) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ELSEIF (IIHW(N).EQ.49) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH2(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH2(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
      TTTT=ABS(SUMA)+ABS(SUMM)+ABS(SUMI)+ABS(SUMPH)
      IF (TTTT.EQ.0.D0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 'NO FLUXES RE-EMITTED FROM INCIDENT BULK IONS '
        CALL LEER(1)
      ELSE
        CALL LEER(1)
C
        WRITE (iunout,*) 'FLUX RE-EMITTED FROM INCIDENT BULK IONS:'
        WRITE (iunout,*) '(RECYCLING SOURCE) '
C  SURFACE AVERAGED TALLY NO. 6
        IF (SUMA.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: ATOMS'
        CALL MASYR1('P-FLUX:  ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
C  SURFACE AVERAGED TALLY NO. 31
        CALL MASYR1('E-FLUX:  ',SUMA2,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA2,LGVRA2,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 12
        IF (SUMM.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: MOLECULES'
        CALL MASYR1('P-FLUX:  ',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
C  SURFACE AVERAGED TALLY NO. 37
        CALL MASYR1('E-FLUX:  ',SUMM2,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM2,LGVRM2,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 18
        IF (SUMI.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: TEST IONS'
        CALL MASYR1('P-FLUX:  ',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
C  SURFACE AVERAGED TALLY NO. 43
        CALL MASYR1('E-FLUX:  ',SUMI2,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI2,LGVRI2,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        ENDIF
C  SURFACE AVERAGED TALLY NO. 24
        IF (SUMPH.NE.0.D0) THEN
        WRITE (iunout,*) 'REEMITTED: PHOTONS'
        CALL MASYR1('P-FLUX:  ',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
C  SURFACE AVERAGED TALLY NO. 49
        CALL MASYR1('E-FLUX:  ',SUMPH2,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH2,LGVRPH2,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        ENDIF
        CALL LEER (1)
        WRITE (iunout,*) 
     .    'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
        WRITE (iunout,*) 
     .    'CONTRIB. FROM INCIDENT BULK IONS (REC. SOURCE)'
        CALL MASR1 ('TOT.PFLX',SUMMTP)
        CALL MASR1 ('TOT.EFLX',SUMMEP)
C
      ENDIF
C
      IF (ABS(SUMMTA)+ABS(SUMMTM)+ABS(SUMMTI)+
     .    ABS(SUMMTPH).NE.0.D0) THEN
        CALL LEER (1)
C
        IF (ILIIN(I).GT.0) THEN
          WRITE (iunout,*) 
     .      'TOTAL REEMITTED "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0) WRITE (iunout,*) 
     .      '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
          CALL LEER (1)
C
        ELSEIF (ILIIN(I).LT.0.AND.ILIIN(I).NE.-3) THEN
          WRITE (iunout,*) 
     .      'TOTAL NEGATIVE "ATOMIC" FLUXES, AMPERE AND WATT'
          IF (SUMMTP.GT.0.D0) WRITE (iunout,*) 
     .      '(EXCLUDING CONTRIB. FROM INCIDENT BULK IONS)'
          CALL MASR1 ('TOT.PFLX',SUMMTA+SUMMTM+SUMMTI+SUMMTPH)
          CALL MASR1 ('TOT.EFLX',SUMMEA+SUMMEM+SUMMEI+SUMMEPH)
        ENDIF
      ENDIF
C
C
C   SPUTTERED FLUXES
C
C   SURFACE AVERAGED TALLY NO.51
C
      SUMMT=0.
      SUMMA=0.
      SUMA1(:,ISTRA)=0._DP
      DO 202 IATM=1,NATMI
        LGVRA1(IATM,ISTRA)=.FALSE.
        IF (LSPTAT) SUMA1(IATM,ISTRA)=SPTAT(IATM,I)
        SUMMT=SUMMT+SUMA1(IATM,ISTRA)
        SUMMA=SUMMA+SUMA1(IATM,ISTRA)
202     CONTINUE
C
      DO 221 N=1,NSIGSI
        IF (IIHW(N).EQ.51) THEN
          DO 222 IATM=1,NATMI
            IF (IGHW(N).NE.IATM) GOTO 222
            VARA1(IATM,ISTRA)=SIGMAW(N,I)
            LGVRA1(IATM,ISTRA)=LOGATM(IATM,ISTRA)
222       CONTINUE
        ENDIF
221   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.52
C
      SUMMM=0.
      SUMM1(:,ISTRA)=0._DP
      DO 203 IMOL=1,NMOLI
        LGVRM1(IMOL,ISTRA)=.FALSE.
        IF (LSPTML) SUMM1(IMOL,ISTRA)=SPTML(IMOL,I)
        SUMMT=SUMMT+SUMM1(IMOL,ISTRA)
        SUMMM=SUMMM+SUMM1(IMOL,ISTRA)
203   CONTINUE
C
      DO 223 N=1,NSIGSI
        IF (IIHW(N).EQ.52) THEN
          DO 224 IMOL=1,NMOLI
            IF (IGHW(N).NE.IMOL) GOTO 224
            VARM1(IMOL,ISTRA)=SIGMAW(N,I)
            LGVRM1(IMOL,ISTRA)=LOGMOL(IMOL,ISTRA)
224       CONTINUE
        ENDIF
223   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.53
C
      SUMMI=0.
      SUMI1(:,ISTRA)=0._DP
      DO 204 IION=1,NIONI
        LGVRI1(IION,ISTRA)=.FALSE.
        IF (LSPTIO) SUMI1(IION,ISTRA)=SPTIO(IION,I)
        SUMMT=SUMMT+SUMI1(IION,ISTRA)
        SUMMI=SUMMI+SUMI1(IION,ISTRA)
204   CONTINUE
C
      DO 225 N=1,NSIGSI
        IF (IIHW(N).EQ.53) THEN
          DO 226 IION=1,NIONI
            IF (IGHW(N).NE.IION) GOTO 226
            VARI1(IION,ISTRA)=SIGMAW(N,I)
            LGVRI1(IION,ISTRA)=LOGION(IION,ISTRA)
226       CONTINUE
        ENDIF
225   CONTINUE
C
C   SURFACE AVERAGED TALLY NO.54
C
      SUMMPH=0.
      SUMPH1(:,ISTRA)=0._DP
      DO IPHOT=1,NPHOTI
        LGVRPH1(IPHOT,ISTRA)=.FALSE.
        IF (LSPTPHT) SUMPH1(IPHOT,ISTRA)=SPTPHT(IPHOT,I)
        SUMMT=SUMMT+SUMPH1(IPHOT,ISTRA)
        SUMMPH=SUMMPH+SUMPH1(IPHOT,ISTRA)
      ENDDO
C
      DO N=1,NSIGSI
        IF (IIHW(N).EQ.54) THEN
          DO IPHOT=1,NPHOTI
            IF (IGHW(N).NE.IPHOT) CYCLE
            VARPH1(IPHOT,ISTRA)=SIGMAW(N,I)
            LGVRPH1(IPHOT,ISTRA)=LOGPHOT(IPHOT,ISTRA)
          ENDDO
        ENDIF
      ENDDO
C
C   SURFACE AVERAGED TALLY NO.55
C
      SUMMP=0.
      SUMP1(:,ISTRA)=0._DP
      DO 205 IPLS=1,NPLSI
        LGVRP1(IPLS,ISTRA)=.FALSE.
        IF (LSPTPL) SUMP1(IPLS,ISTRA)=SPTPL(IPLS,I)
        SUMMT=SUMMT+SUMP1(IPLS,ISTRA)
        SUMMP=SUMMP+SUMP1(IPLS,ISTRA)
205   CONTINUE
C
      DO 227 N=1,NSIGSI
        IF (IIHW(N).EQ.55) THEN
          DO 228 IPLS=1,NPLSI
            IF (IGHW(N).NE.IPLS) GOTO 228
            VARP1(IPLS,ISTRA)=SIGMAW(N,I)
            LGVRP1(IPLS,ISTRA)=LOGPLS(IPLS,ISTRA)
228       CONTINUE
        ENDIF
227   CONTINUE
C
      IF (SUMMT.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO FLUXES SPUTTERED FROM THIS SURFACE         ')
        CALL LEER(1)
        GOTO 206
      ENDIF
      CALL LEER(1)
      WRITE (iunout,*) 'FLUX SPUTTERED FROM SURFACE:'
      IF (SUMMA.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 51
        CALL MASYR1('ATOMS    ',SUMA1,LOGATM,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
        CALL MASYR1('ST.DEV.% ',VARA1,LGVRA1,ISTRA,0,NATM,0,NSTRA,
     .               TEXTS(NSPH+1))
      ENDIF
      IF (SUMMM.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 52
        CALL MASYR1('MOLECULES',SUMM1,LOGMOL,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
        CALL MASYR1('ST.DEV.% ',VARM1,LGVRM1,ISTRA,0,NMOL,0,NSTRA,
     .               TEXTS(NSPA+1))
      ENDIF
      IF (SUMMI.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 53
        CALL MASYR1('TEST IONS',SUMI1,LOGION,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
        CALL MASYR1('ST.DEV.% ',VARI1,LGVRI1,ISTRA,0,NION,0,NSTRA,
     .               TEXTS(NSPAM+1))
      ENDIF
      IF (SUMMPH.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 54
        CALL MASYR1('PHOTONS',SUMPH1,LOGPHOT,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
        CALL MASYR1('ST.DEV.% ',VARPH1,LGVRPH1,ISTRA,0,NPHOT,0,NSTRA,
     .               TEXTS(0+1))
      ENDIF
      IF (SUMMP.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 55
        CALL MASYR1('BULK IONS',SUMP1,LOGPLS,ISTRA,0,NPLS,0,NSTRA,
     .               TEXTS(NSPAMI+1))
        CALL MASYR1('ST.DEV.% ',VARP1,LGVRP1,ISTRA,0,NPLS,0,NSTRA,
     .               TEXTS(NSPAMI+1))
      ENDIF
      CALL LEER (1)
      CALL MASR1 ('TOT. FLX',SUMMT)

206   CONTINUE
C
C   ADDITIONAL SURFACE AVERAGED TALLIES
C
C   SURFACE AVERAGED TALLY NO. 57
C
      SUMMS=0.
      SUMS(:,ISTRA)=0._DP
      DO 302 IADS=1,NADSI
        TEXTA(IADS)=TXTSPW(IADS,NTLSA)
        LGVARS(IADS,ISTRA)=.FALSE.
        IF (LADDS) SUMS(IADS,ISTRA)=ADDS(IADS,I)
        LOGADS(IADS,ISTRA)=LADDS .AND. (ADDS(IADS,I).NE.0.)
        SUMMS=SUMMS+SUMS(IADS,ISTRA)
302   CONTINUE
C
      DO 321 N=1,NSIGSI
        IF (IIHW(N).NE.57) GOTO 321
        IF (IGHW(N).EQ.0) GOTO 321
        DO 322 IADS=1,NADSI
          IF (IGHW(N).NE.IADS) GOTO 322
          VARS(IADS,ISTRA)=SIGMAW(N,I)
          LGVARS(IADS,ISTRA)=LOGADS(IADS,ISTRA)
322     CONTINUE
321   CONTINUE
C
      IF (SUMMS.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO ADDITIONAL SURFACE TALLIES AT THIS SURFACE ')
        CALL LEER(1)
        GOTO 305
      ENDIF
      CALL LEER(1)
      WRITE (iunout,*) 'ADDITIONAL SURFACE TALLIES:'
      IF (SUMMS.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 57
        CALL MASYR1('ADD.TALLY',SUMS,LOGADS,ISTRA,0,NADS,0,NSTRA,TEXTA)
        CALL MASYR1('ST.DEV.% ',VARS,LGVARS,ISTRA,0,NADS,0,NSTRA,TEXTA)
      ENDIF
      CALL LEER (1)
C
305   CONTINUE
C
C
C   ALGEBRAIC SURFACE AVERAGED TALLIES
C
C   SURFACE AVERAGED TALLY NO. 58
C
      SUMMS=0.
      SUML(:,ISTRA)=0._DP
      DO 402 IALS=1,NALSI
        TEXTL(IALS)=TXTSPW(IALS,NTLSR)
        LGVARL(IALS,ISTRA)=.FALSE.
        IF (LALGS) SUML(IALS,ISTRA)=ALGS(IALS,I)
        LOGALS(IALS,ISTRA)=LALGS .AND. (ALGS(IALS,I).NE.0.)
        SUMMS=SUMMS+SUML(IALS,ISTRA)
402   CONTINUE
C
      DO 421 N=1,NSIGSI
        IF (IIHW(N).NE.58) GOTO 421
        IF (IGHW(N).EQ.0) GOTO 421
        DO 422 IALS=1,NALSI
          IF (IGHW(N).NE.IALS) GOTO 422
          VARL(IALS,ISTRA)=SIGMAW(N,I)
          LGVARL(IALS,ISTRA)=LOGALS(IALS,ISTRA)
422     CONTINUE
421   CONTINUE
C
      IF (SUMMS.EQ.0.D0) THEN
        CALL LEER(1)
        CALL MASAGE ('NO ALGEBRAIC SURFACE TALLIES AT THIS SURFACE ')
        CALL LEER(1)
        GOTO 405
      ENDIF
      CALL LEER(1)
      WRITE (iunout,*) 'ALGEBRAIC SURFACE TALLIES:'
      DO IALS=1,NALSI
        WRITE (iunout,*) 'NO ',IALS,'  ',TXTTLW(IALS,NTLSR)
      ENDDO
      IF (SUMMS.NE.0.D0) THEN
C  SURFACE AVERAGED TALLY NO. 58
        CALL MASYR1('         ',SUML,LOGALS,ISTRA,0,NALS,0,NSTRA,TEXTL)
      ENDIF
      CALL LEER (1)
C
405   CONTINUE
C

C  SPECTRA
      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO ISPC=1,NADSPC
        IF ((ESTIML(ISPC)%PSPC%ISRFCLL == 0) .AND.
     .      (ESTIML(ISPC)%PSPC%ISPCSRF == I)) THEN
          CALL LEER (1)
          WRITE (iunout,'(A33)') ' SPECTRUM CALCULATED FOR SURFACE '
          IT = ESTIML(ISPC)%PSPC%ISPCTYP
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A40)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
          END IF
          WRITE (iunout,'(A20,A8)') ' TYPE OF PARTICLE : ',
     .                       TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
          IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
            WRITE (iunout,'(A10,10X,A16)') ' SPECIES :',
     .                                     'SUM OVER SPECIES'
          ELSE
            WRITE (iunout,'(A10,10X,A8)') ' SPECIES :',
     .             TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .                   ESTIML(ISPC)%PSPC%IPRSP)
          END IF
          WRITE (iunout,'(A22,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .           ESTIML(ISPC)%PSPC%SPCINT
          IF (NSIGI_SPC > 0)
     .      WRITE (iunout,'(A22,ES12.4)') ' STANDARD DEVIATION   ',
     .           ESTIML(ISPC)%PSPC%SGMS
        END IF
      END DO

        PRINTED(I) = .TRUE.
10000 CONTINUE  ! END OF LOOP OVER SURFACES, FOR WHICH PRINTOUT WAS REQUESTED

      RETURN
9999  FORMAT (1X,A32)
      END
C ===== SOURCE: outlst.f
      SUBROUTINE OUTLST

      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COMPRT, ONLY: IUNOUT
      USE CLAST
      USE COMXS

      IMPLICIT NONE

      REAL(DP) :: SMMEAN
      INTEGER :: IWR, IRCX, IREL

      call leer(1)
      iwr=0
      do ircx=1,nrcxi
        if (iflrcx(ircx).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOCX '
            write (iunout,*) 
     .        'IRCX, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLS'
            iwr=1
          endif
          SMMEAN=xcmean(ircx)/(ncmean(ircx)+eps60)
          write (iunout,*) ircx,SMMEAN,NCMEAN(IRCX)
        endif
      enddo
      call leer(1)
      iwr=0
      do irel=1,nreli
        if (iflrel(irel).gt.0) then
          if (iwr.eq.0) then
            WRITE (iunout,*) 'REJECTION SAMPLING EFFICIENCY IN VELOEL '
            write (iunout,*) 
     .        'IREL, MEAN NO. OF SAMPLING, TOTAL NO. OF CALLS'
            iwr=1
          endif
          SMMEAN=xemean(irel)/(nemean(irel)+eps60)
          write (iunout,*) irel,SMMEAN,NEMEAN(IREL)
        endif
      enddo
      call leer(1)
      RETURN
      END
C ===== SOURCE: outpla.f
C
C
      SUBROUTINE OUTPLA(ICAL)
C  This routine prints background tallies as requested in input block 11.
C  ical=0: called after initialization phase
C  ical=1: called in post processing phase from iteration loop
C          Some background parameters
C          may have been modified due to iteration loop.
C          Print only the modified tallies.
C
C  PRINT INPUT TALLIES ONTO OUTPUT FILE IUNOUT 
C
      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CSPEZ
      USE CTRCEI
      USE CGEOM
      USE CTEXT
      USE COUTAU
      USE CSPEI
      USE CINIT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICAL
      REAL(DP), ALLOCATABLE :: HELPP(:),HELPW(:)
      REAL(DP) :: TALTYP(NTALI)
      REAL(DP) :: TALAV, HELPI, TALTOT, TOTAL
      INTEGER :: IR, IP, IT, I, NBLCKA, IB, IPRV, ITAL, NXM, NYM, NZM,
     .           ITALI, K, NF, NFTI, NFTE

C  INDICATOR FOR THE TALLIES THAT MAY HAVE BEEN MODIFIED IN POST PROCESSING
C   CURRENTLY: DENSITY, BULK ION TEMP, AND BULK ION DRIFT VELOCITY
      INTEGER :: JPRTAL(5) = (/-2,-4,-5,-6,-7/)
C
C  TYPE OF TALLY: TALTYP=0: #              (#-UNITS)
C                 TALTYP=1: #-DENSITY      (#-UNITS/CM**3)
C                 TALTYP=2: VOLUME         (CM**3)
C                 TALTYP=3: DIMENSIONLESS  (1)
C                 TALTYP=4: UNKNOWN        (?)
      TALTYP(1)=0
      TALTYP(2)=0
      TALTYP(3)=1
      TALTYP(4)=1
      TALTYP(5)=0
      TALTYP(6)=0
      TALTYP(7)=0
      TALTYP(8)=3
      TALTYP(9)=3
      TALTYP(10)=3
      TALTYP(11)=3
      TALTYP(12)=4
      TALTYP(13)=0
      TALTYP(14)=2
      TALTYP(15)=3
      TALTYP(16)=3
      TALTYP(17)=3

      IF (ICAL == 1) THEN
!  IS ANY DENSITY MODEL DEFINED ?
        IF (ALL(CDENMODEL == REPEAT(' ',LEN(CDENMODEL)))) RETURN
!  IS OUTPUT OF AN INPUT TALLY ASKED FOR?
        DO IPRV=1,NVOLPR
          ITAL=NPRTLV(IPRV)
          IF (ANY(JPRTAL == ITAL)) THEN
            CALL LEER (2)
            CALL HEADNG 
     .        ('BACKGROUND TALLIES CHANGED IN POST PROCESSING',45)
            EXIT
          END IF
        END DO
!  NO OUTPUT OF INPUT TALLIES REQUIRED
        IF (IPRV > NVOLPR) RETURN
      END IF

      IF (NVOLPR > 0) THEN
        ALLOCATE (HELPP(NRAD))
        ALLOCATE (HELPW(NRAD))
      END IF
C
C  PRINT INPUT VOLUME AVERAGED TALLIES, WHICH HAVE BEEN SELECTED
C
      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)
      DO 100 IPRV=1,NVOLPR
        ITAL=NPRTLV(IPRV)
        IF ((ICAL == 1) .AND. (ALL(JPRTAL .NE. ITAL))) CYCLE
        IF (ITAL.LT.0) THEN
          ITALI=-ITAL
          NF=NFRSTP(ITALI)
          NFTI=1
          NFTE=NFSTPI(ITALI)
          IF (NSPEZV(IPRV,1).GT.0) THEN
            NFTI=NSPEZV(IPRV,1)
            NFTE=MAX(NFTI,NSPEZV(IPRV,2))
          ENDIF
          DO 119 K=NFTI,NFTE
            IF (K.GT.NFSTPI(ITALI)) THEN
              CALL LEER(1)
              WRITE (iunout,*) 
     .          'SPECIES INDEX OUT OF RANGE IN SUBR. OUTPLA'
              WRITE (iunout,*) 'ITALI, K, ',ITALI,K
              CALL LEER(1)
              GOTO 119
            ENDIF
            IF ((ICAL == 1) .AND. (VERIFY(CDENMODEL(K),' ') == 0)) CYCLE
            SELECT CASE (ITALI)
            CASE (1)
              HELPP(1:NSBOX) = TEIN(1:NSBOX)
            CASE (2)
              HELPP(1:NSBOX) = TIIN(MPLSTI(K),1:NSBOX)
            CASE (3)
              HELPP(1:NSBOX) = DEIN(1:NSBOX)
            CASE (4)
              HELPP(1:NSBOX) = DIIN(K,1:NSBOX)
            CASE (5)
              HELPP(1:NSBOX) = VXIN(MPLSV(K),1:NSBOX)
            CASE (6)
              HELPP(1:NSBOX) = VYIN(MPLSV(K),1:NSBOX)
            CASE (7)
              HELPP(1:NSBOX) = VZIN(MPLSV(K),1:NSBOX)
            CASE (8)
              HELPP(1:NSBOX) = BXIN(1:NSBOX)
            CASE (9)
              HELPP(1:NSBOX) = BYIN(1:NSBOX)
            CASE (10)
              HELPP(1:NSBOX) = BZIN(1:NSBOX)
            CASE (11)
              HELPP(1:NSBOX) = BFIN(1:NSBOX)
            CASE (12)
              HELPP(1:NSBOX) = ADIN(K,1:NSBOX)
            CASE (13)
              HELPP(1:NSBOX) = EDRIFT(K,1:NSBOX)
            CASE (14)
              HELPP(1:NSBOX) = VOL(1:NSBOX)
            CASE (15)
              HELPP(1:NSBOX) = WGHT(K,1:NSBOX)
            CASE (16)
              HELPP(1:NSBOX) = BXPERP(1:NSBOX)
            CASE (17)
              HELPP(1:NSBOX) = BYPERP(1:NSBOX)
            CASE DEFAULT
              WRITE (iunout,*) ' WRONG TALLY NUMBER, ITAL = ',ITAL
              WRITE (iunout,*) ' NO OUTPUT PERFORMED '
              CALL LEER(1)
              GOTO 100
            END SELECT
C
            TOTAL=0.D0
            DO 121 IB=1,NBMLT
            NBLCKA=NSTRD*(IB-1)
            DO 121 IR=1,NXM
            DO 121 IP=1,NYM
            DO 121 IT=1,NZM
              I=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
              IF (ITALI.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DEIN(I)*VOL(I)
              ELSEIF (ITALI.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.EQ.3.OR.ITALI.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
                HELPW(I)=VOL(I)
              ELSEIF (ITALI.EQ.5.OR.ITALI.EQ.6.OR.ITALI.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.GE.8.AND.ITALI.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.16.OR.ITALI.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.12.OR.ITALI.EQ.14.OR.ITALI.EQ.15) THEN
C  ADDITIONAL TALLY, CELL VOLUME ,WEIGHT FUNCTION " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.13) THEN
C  ION DRIFT ENERGY
                HELPW(I)=DIIN(K,I)*VOL(I)
              ENDIF
              TOTAL=TOTAL+HELPW(I)
121         CONTINUE
C
C  SAME LOOP AGAIN, OVER ADDITIONAL CELL REGION
            DO 122 I=NSURF+1,NSURF+NRADD
              IF (ITALI.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DEIN(I)*VOL(I)
              ELSEIF (ITALI.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.EQ.3.OR.ITALI.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
                HELPW(I)=VOL(I)
              ELSEIF (ITALI.EQ.5.OR.ITALI.EQ.6.OR.ITALI.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.GE.8.AND.ITALI.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.16.OR.ITALI.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.12.OR.ITALI.EQ.14.OR.ITALI.EQ.15) THEN
C  ADDITIONAL TALLY (NO.12)
C  CELL VOLUME      (NO.14)
C  WEIGHT FUNCTION  (NO.15)
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.13) THEN
C  ION DRIFT ENERGY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ENDIF
              TOTAL=TOTAL+HELPW(I)
122         CONTINUE
C
            IF (ITALI.NE.NTALO) THEN
              CALL INTTAL (HELPP,HELPW,1,1,NSBOX,HELPI,
     .                     NR1ST,NP2ND,NT3RD,NBMLT)
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL INTVOL (HELPP,      1,1,NSBOX,HELPI,
     .                     NR1ST,NP2ND,NT3RD,NBMLT)
            ENDIF
            TALTOT=HELPI
            TALAV=HELPI/(TOTAL+EPS60)
C
            IF (TALTOT.EQ.0.D0) GOTO 118
C
            IF (ITALI.NE.NTALO) THEN
              CALL PRTTAL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL PRTVOL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
            ENDIF
            CALL LEER(2)
            IF (TALTYP(ITALI).EQ.0) THEN
              CALL MASAGE
     .        ('WEIGHTED MEAN VALUE ("UNITS")               ')
              CALL MASR1 ('MEAN:   ',TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.1) THEN
              CALL MASAGE
     .        ('WEIGHTED TOTAL ("UNITS*CM**3"), MEAN ("UNITS")')
              CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.2) THEN
              CALL MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS")                         ')
              CALL MASR1 ('TOTAL:  ',TALTOT)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.3) THEN
              CALL MASAGE
     .        ('ARITHMETIC MEAN ("UNITS")                          ')
              CALL MASR1 ('MEAN:   ',TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.4) THEN
              CALL MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS"), MEAN ("UNITS")')
              CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL LEER(3)
            ENDIF
            GOTO 119
C
118         CONTINUE
            CALL PRTTAL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),TXTPUN(K,ITALI),
     .                  HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,-1,0)
            CALL MASAGE('IDENTICAL ZERO, NOT PRINTED                  ')
            CALL LEER(2)
119       CONTINUE
        ENDIF
100   CONTINUE
      CALL LEER(2)

      IF (NVOLPR > 0) THEN
        DEALLOCATE (HELPP)
        DEALLOCATE (HELPW)
      END IF
C
      RETURN
      END
c slmod begin


C ===== SOURCE: outspec.f
!pb  25.10.06:  format specifications corrected


      SUBROUTINE OUTSPEC(ISTRA)
c
c      SUBROUTINE OUTSPEC
c slmod end

      USE PRECISION
      USE PARMMOD
      USE COMPRT, ONLY: IUNOUT
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CTRCEI
      USE CTEXT
      USE CSDVI
c slmod begin
      USE MOD_INTERFACE
c slmod end
      IMPLICIT NONE
      INTEGER :: IADTYP(0:4)
      INTEGER :: IOUT, ISPC, I, IT, IE
      REAL(DP) :: EN
      CHARACTER(10) :: TEXTYP(0:4)
c slmod begin
      INTEGER,INTENT(IN) :: ISTRA
      CHARACTER FILE*128,TAG*3,UNITS*32
c slmod end

C  SPECTRA

      IOUT = 20
      OPEN (UNIT=IOUT,FILE='spectra.out')

      TEXTYP(0) = 'PHOTONS   '
      TEXTYP(1) = 'ATOMS     '
      TEXTYP(2) = 'MOLECULES '
      TEXTYP(3) = 'TEST IONS '
      TEXTYP(4) = 'BULK IONS '
      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)

      DO ISPC=1,NADSPC
        I = ESTIML(ISPC)%PSPC%ISPCSRF
        IT = ESTIML(ISPC)%PSPC%ISPCTYP

        WRITE (IOUT,*)
        WRITE (IOUT,*)
        WRITE (IOUT,*)

        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0)  THEN
          IF (I > NLIM) THEN
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' NONDEFAULT STANDARD SURFACE ',I-NLIM
          ELSE
            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                     ' ADDITIONAL SURFACE ',I
          END IF        
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))') 
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
          ELSE IF (IT == 2) THEN
            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
          END IF

        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 1)  THEN
          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
     .                   ' SCORING CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))') 
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF

        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 2)  THEN
          WRITE (IOUT,'(A,A)') ' SPECTRUM CALCULATED FOR',
     .                   ' GEOMETRICAL CELL ',I
          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
            WRITE (iunout,'(A,3(ES12.4,A1))') 
     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
          END IF
          IF (IT == 1) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
          ELSEIF (IT == 2) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
          ELSEIF (IT == 3) THEN
            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
          END IF
        END IF

        WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
     .         TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
        IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
          WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
     .                'SUM OVER SPECIES'
        ELSE
          WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
     .          TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
     .          ESTIML(ISPC)%PSPC%IPRSP)
        END IF

        WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
     .         ESTIML(ISPC)%PSPC%SPCMIN
        WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
     .         ESTIML(ISPC)%PSPC%SPCMAX
        WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ',
     .         ESTIML(ISPC)%PSPC%NSPC
        WRITE (IOUT,*)
        IF (ESTIML(ISPC)%PSPC%SPCINT > EPS60) THEN
          IF (NSIGI_SPC == 0) THEN
            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE)
            END DO
          ELSE
            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
              EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
              WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
     .               ESTIML(ISPC)%PSPC%SPC(IE),
     .               ESTIML(ISPC)%PSPC%SDV(IE)
            END DO
          END IF
        ELSE
          WRITE (IOUT,'(A)') ' SPECTRUM IDENTICAL 0 '
        END IF
        WRITE (IOUT,*)
        WRITE (IOUT,'(A,ES12.4)') ' INTEGRAL OF SPECTRUM ',
     .         ESTIML(ISPC)%PSPC%SPCINT
        IF (NSIGI_SPC > 0)
     .    WRITE (IOUT,'(A,ES12.4)') ' STANDARD DEVIATION  ',
     .                   ESTIML(ISPC)%PSPC%SGMS
c slmod begin
c         WRITE(0,*) 'ISTRA=',istra
        IF (istra.EQ.0) THEN
          tag = 'sum'
        ELSE
          WRITE(tag,'(I3.3)') istra 
        ENDIF
          WRITE(file,'(A,I4.4,A,I7.7,A)') 
     .      'idl.eirene_spectrum_',ispc,'_',i,'_'//tag
c        WRITE(0,*) TRIM(file)
        UNITS='?'
        IF (IT == 1) UNITS='Amps/BIN(eV)'
        IF (IT == 2) UNITS='Watt/BIN(eV)'
        CALL inOpenInterface(TRIM(file),ITF_WRITE)
        IF (NSIGI_SPC == 0) THEN
          DO IE=1, ESTIML(ISPC)%PSPC%NSPC
            EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .           (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
            CALL inPutData(EN                       ,'BIN' ,'eV')
c            IF (ispc.eq.5) 
c     .        WRITE(0,*) 'debug: spec',ie,ispc,ESTIML(ISPC)%PSPC%SPC(IE)    ! *** LEFT OFF ***
            CALL inPutData(ESTIML(ISPC)%PSPC%SPC(IE),'FLUX',TRIM(UNITS))  
            CALL inPutData(-1.0D0                   ,'STDE','N/A')
          END DO
        ELSE
          DO IE=1, ESTIML(ISPC)%PSPC%NSPC
            EN = ESTIML(ISPC)%PSPC%SPCMIN +
     .           (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
            CALL inPutData(EN                       ,'BIN' ,'eV')
            CALL inPutData(ESTIML(ISPC)%PSPC%SPC(IE),'FLUX',TRIM(UNITS))
            CALL inPutData(ESTIML(ISPC)%PSPC%SDV(IE),'STDE','N/A')
          END DO
        END IF
        CALL inPutData(ESTIML(ISPC)%PSPC%SPCINT,'INTEGRAL','?')      
        CALL inPutData(ESTIML(ISPC)%PSPC%IPRTYP,'SPECIES_TYPE','N/A')
        CALL inPutData(it,'SPECTRUM_TYPE','N/A')
c        write(0,*) 'index check',i,nlim
        CALL inPutData(I,'INDEX','N/A')
c        IF (I >  NLIM) CALL inPutData(-(I-NLIM),'INDEX','N/A')
c        IF (I <= NLIM) CALL inPutData(I        ,'INDEX','N/A')
        CALL inPutData(ESTIML(ISPC)%PSPC%SPCMIN,'MIN_VALUE','eV')
        CALL inPutData(ESTIML(ISPC)%PSPC%SPCMAX,'MAX_VALUE','eV')
        CALL inPutData(ISTRA                   ,'STRATUM'   ,'N/A')
        CALL inCloseInterface
c slmod end
      END DO

      RETURN
      END SUBROUTINE OUTSPEC

c
c      SUBROUTINE OUTSPEC
c
c      USE PRECISION
c      USE PARMMOD
c      USE COMPRT, ONLY: IUNOUT
c      USE COMUSR
c      USE CESTIM
c      USE CCONA
c      USE CTRCEI
c      USE CTEXT
c      USE CSDVI
c
c      IMPLICIT NONE
c      INTEGER :: IADTYP(0:4)
c      INTEGER :: IOUT, ISPC, I, IT, IE
c      REAL(DP) :: EN
c      CHARACTER(10) :: TEXTYP(0:4)
c
cC  SPECTRA
c
c      IOUT = 20
c      OPEN (UNIT=IOUT,FILE='spectra.out')
c
c      TEXTYP(0) = 'PHOTONS   '
c      TEXTYP(1) = 'ATOMS     '
c      TEXTYP(2) = 'MOLECULES '
c      TEXTYP(3) = 'TEST IONS '
c      TEXTYP(4) = 'BULK IONS '
c      IADTYP(0:4) = (/ 0, NSPH, NSPA, NSPAM, NSPAMI /)
c
c      DO ISPC=1,NADSPC
c        I = ESTIML(ISPC)%PSPC%ISPCSRF
c        IT = ESTIML(ISPC)%PSPC%ISPCTYP
c
c        WRITE (IOUT,*)
c        WRITE (IOUT,*)
c        WRITE (IOUT,*)
c
c        IF (ESTIML(ISPC)%PSPC%ISRFCLL == 0)  THEN
c          IF (I > NLIM) THEN
c            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
c     .                     ' NONDEFAULT STANDARD SURFACE ',I-NLIM
c          ELSE
c            WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
c     .                     ' ADDITIONAL SURFACE ',I
c          END IF        
c          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
c            WRITE (iunout,'(A,3(ES12.4,A1))') 
c     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
c     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
c          END IF
c          IF (IT == 1) THEN
c            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
c     .                'INCIDENT PARTICLE FLUX IN AMP/BIN(EV)   '
c          ELSE IF (IT == 2) THEN
c            WRITE (IOUT,'(A,A)') ' TYPE OF SPECTRUM : ',
c     .                'INCIDENT ENERGY FLUX IN WATT/BIN(EV)    '
c          END IF
c
c        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 1)  THEN
c          WRITE (IOUT,'(A,A,I6)') ' SPECTRUM CALCULATED FOR',
c     .                   ' SCORING CELL ',I
c          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
c            WRITE (iunout,'(A,3(ES12.4,A1))') 
c     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
c     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
c          END IF
c          IF (IT == 1) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
c          ELSEIF (IT == 2) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
c          ELSEIF (IT == 3) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
c          END IF
c
c        ELSE IF (ESTIML(ISPC)%PSPC%ISRFCLL == 2)  THEN
c          WRITE (IOUT,'(A,A)') ' SPECTRUM CALCULATED FOR',
c     .                   ' GEOMETRICAL CELL ',I
c          IF (ESTIML(ISPC)%PSPC%IDIREC > 0) THEN
c            WRITE (iunout,'(A,3(ES12.4,A1))') 
c     .      ' IN DIRECTION (',ESTIML(ISPC)%PSPC%SPCVX,',',
c     .      ESTIML(ISPC)%PSPC%SPCVY,',',ESTIML(ISPC)%PSPC%SPCVZ,')'
c          END IF
c          IF (IT == 1) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL PARTICLE DENSITY IN #/CM**3/BIN(EV)   '
c          ELSEIF (IT == 2) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL ENERGY DENSITY IN EV/CM**3/BIN(EV)    '
c          ELSEIF (IT == 3) THEN
c            WRITE (iunout,'(A20,A)') ' TYPE OF SPECTRUM : ',
c     .        'SPECTRAL MOMENTUM DENSITY IN (G*CM/S)/CM**3/BIN(EV)    '
c          END IF
c        END IF
c
c        WRITE (IOUT,'(A20,A9)') ' TYPE OF PARTICLE : ',
c     .         TEXTYP(ESTIML(ISPC)%PSPC%IPRTYP)
c        IF (ESTIML(ISPC)%PSPC%IPRSP == 0) THEN
c          WRITE (IOUT,'(A10,10X,A16)') ' SPECIES :',
c     .                'SUM OVER SPECIES'
c        ELSE
c          WRITE (IOUT,'(A10,10X,A8)') ' SPECIES :',
c     .          TEXTS(IADTYP(ESTIML(ISPC)%PSPC%IPRTYP)+
c     .          ESTIML(ISPC)%PSPC%IPRSP)
c        END IF
c
c        WRITE (IOUT,'(A15,5X,ES12.4)') ' MINIMAL ENERGY ',
c     .         ESTIML(ISPC)%PSPC%SPCMIN
c        WRITE (IOUT,'(A15,5X,ES12.4)') ' MAXIMAL ENERGY ',
c     .         ESTIML(ISPC)%PSPC%SPCMAX
c        WRITE (IOUT,'(A16,4x,I6)') ' NUMBER OF BINS ',
c     .         ESTIML(ISPC)%PSPC%NSPC
c        WRITE (IOUT,*)
c        IF (ESTIML(ISPC)%PSPC%SPCINT > EPS60) THEN
c          IF (NSIGI_SPC == 0) THEN
c            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
c              EN = ESTIML(ISPC)%PSPC%SPCMIN +
c     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
c              WRITE (IOUT,'(I6,2ES12.4)') IE,EN,
c     .               ESTIML(ISPC)%PSPC%SPC(IE)
c            END DO
c          ELSE
c            DO IE=1, ESTIML(ISPC)%PSPC%NSPC
c              EN = ESTIML(ISPC)%PSPC%SPCMIN +
c     .             (IE-0.5)*ESTIML(ISPC)%PSPC%SPCDEL
c              WRITE (IOUT,'(I6,3ES12.4)') IE,EN,
c     .               ESTIML(ISPC)%PSPC%SPC(IE),
c     .               ESTIML(ISPC)%PSPC%SDV(IE)
c            END DO
c          END IF
c        ELSE
c          WRITE (IOUT,'(A)') ' SPECTRUM IDENTICAL 0 '
c        END IF
c        WRITE (IOUT,*)
c        WRITE (IOUT,'(A,ES12.4)') ' INTEGRAL OF SPECTRUM ',
c     .         ESTIML(ISPC)%PSPC%SPCINT
c        IF (NSIGI_SPC > 0)
c     .    WRITE (IOUT,'(A,ES12.4)') ' STANDARD DEVIATION  ',
c     .                   ESTIML(ISPC)%PSPC%SGMS
c      END DO
c
c      RETURN
c      END SUBROUTINE OUTSPEC
c slmod end
