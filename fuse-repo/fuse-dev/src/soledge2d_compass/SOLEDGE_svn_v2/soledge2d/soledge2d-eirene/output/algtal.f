C
C
      SUBROUTINE EIRENE_ALGTAL
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_CESTIM
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CGEOM
      USE EIRMOD_COMPRT
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_CSPEI
 
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
        CALL EIRENE_ALGEBR (HCHR,OPER,IZIF,CONST,NOP)
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
!pb          IF (IZIF(2,IOP)*IZIF(4,IOP) < 0.D0) GOTO 94
          IF ((IZIF(2,IOP)*IZIF(4,IOP) < 0.D0) .AND.
     .        (NRAD /= NRTAL)) GOTO 94
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
              CALL EIRENE_LEER(1)
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
              CALL EIRENE_LEER(1)
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
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
          WRITE (iunout,*) ' TALLY NUMBER OUT OF RANGE '
          WRITE (iunout,*)
     .      ' CHECK INPUT FOR ADDITIONAL VOLUME TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          GOTO 200
C
91        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
          WRITE (iunout,*) ' SPECIES INDEX OUT OF RANGE '
          WRITE (iunout,*)
     .      ' CHECK INPUT FOR ADDITIONAL VOLUME TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          GOTO 200
C
92        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
          WRITE (iunout,*) ' WRONG NUMBER OF INTERMEDIATE RESULT FOUND '
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)')
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
93        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
          WRITE (iunout,*) ' OPERATOR NOT FORESEEN '
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)')
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
94        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
          WRITE (iunout,*)
     .      ' ARGUMENTS OF OPERATION HAVE DIFFERENT SPACING '
          WRITE (iunout,*) ' NO CALCULATION IS DONE FOR TALLY NO. ',IALV
          WRITE (iunout,*) CHRTAL(IALV)
          WRITE (iunout,'(1X,A,4I4)')
     .          (OPER(J),(IZIF(K,J),K=1,4),J=1,NOP)
          GOTO 200
C
95        CONTINUE
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
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
        CALL EIRENE_ALGEBR (HCHR,OPER,IZIF,CONST,NOP)
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
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGEBR '
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
          WRITE (iunout,*) ' ERROR IN SUBROUTINE EIRENE_ALGTAL '
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
      CALL EIRENE_EXIT_OWN(1)
      END
