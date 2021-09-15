C
C
      SUBROUTINE CDEF(AL,JI,JE,K,COU,NTE,CF,LEXP,LTEST,LSUM)
C
C  EIRENE ATOMIC DATA , DEFAULT OR FROM FILE: POLYNOM FIT FORMAT
C
C  INPUT:
C    AL(J),J=1,NTE
C  OUTPUT:
C    MAXWELLIAN RATES, AL=LN(KT), KT IN (EV)
C  K
C       K>0: A&M DATA FROM FILES HYDHEL, METHANE OR AMJUEL
C           K:  NUMBER OF REACTION IN EIRENE "CREAC"-ARRAY
C
C       K<0: HARD WIRED EIRENE DEFAULT ATOMIC AND MOLECULAR DATA PACKAGE
C       EACH NUMBER ABS(K) CORRESPONDS TO ONE SPECIFIC REACTION DATA FIT
C       (SEE COMMENTS BELOW)
C  JI,JE
C  FOR JI<=J<=JE RETURN:
C       1<=J<=9:  JTH ENERGY COEFFICIENT OF TWO PARAM. FITS
C                 AT TEMPERATUR KT (EV)
C                 (E.G.: RATES FOR HEAVY PARTICLE INTERACTIONS)
C                 IF ONLY ONE FIT AVAILABLE, ITS INDEX IS J=1
C                 E.Q. FOR ELECT. IMP. RATES AS FUNCTION OF TE
C
C  FIT FROM JANEV ET AL, PPPL-TM-368, 1985  (PREPRINT) OR:
C           SPRINGER SERIES ON ATOMS AND PLASMAS, VOL 4, 1987
C
      USE PRECISION
      USE PARMMOD
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: AL(*)
      REAL(DP), INTENT(OUT) :: COU(0:9,*),CF(9,0:9)
      INTEGER, INTENT(IN) :: JI, JE, K, NTE
      LOGICAL, INTENT(IN) :: LEXP, LTEST, LSUM
      REAL(DP) :: CFP(9,11)
      REAL(DP) :: CF1(9),CF2(9),CF3(9),
     .          CF4(9),CF5(9),CF6(9),CF7(9),
     .          CF8(9),CF9(9),CF10(9),CF11(9)
      REAL(DP) :: CTEST
      INTEGER :: ICELL, J, II
C
      EQUIVALENCE (CFP(1,1),CF1(1)),(CFP(1,2),CF2(1)),
     .            (CFP(1,3),CF3(1)),(CFP(1,4),CF4(1)),
     .            (CFP(1,5),CF5(1)),(CFP(1,6),CF6(1)),
     .            (CFP(1,7),CF7(1)),(CFP(1,8),CF8(1)),
     .            (CFP(1,9),CF9(1)),(CFP(1,10),CF10(1)),
     .            (CFP(1,11),CF11(1))
C
C  K>0:
C  DATA FROM ARRAY CREAC(9,-1:9,K)
C
C -K=1:   E + HE --> 2E + HE+
C  RATE COEFFICIENT, JANEV, 2.3.9
      DATA CF1
     ./-4.409865E+01,2.391597E+01,-1.075323E+01,3.058039,
     . -5.685119E-01,6.795391E-02,-5.009056E-03,2.067236E-04,
     . -3.649161E-06/
C -K=2:   FREE
C -K=3:   FREE
C -K=4:   E + H --> H+ + 2E
C  RATE COEFFICIENT, JANEV, 2.1.5
      DATA CF4
     ./-3.271397E+01,1.353656E+01,-5.739329E 00,1.563155E 00,
     . -2.877056E-01,3.482560E-02,-2.631976E-03,1.119544E-04,
     . -2.039150E-06/
C -K=5:  E + H2 --> H + H + E
C  RATE COEFFICIENT, JANEV, 2.2.5, PREPRINT (CORRECT), NOT "BOOK"
      DATA CF5
     ./-2.7872175E+01,1.0522527E+01,-4.9732123E+00,
     .  1.4511982E+00,-3.0627906E-01,4.4333795E-02,
     . -4.0963442E-03, 2.1596703E-04,-4.9285453E-06/
C -K=6:  E + H2 --> H+ + H + 2E
C  RATE COEFFICIENT, JANEV, 2.2.10
      DATA CF6
     ./-3.834597E+01,1.426322E+01,-5.826467E+00,
     .  1.727941E+00,-3.598121E-01,4.822199E-02,
     . -3.909403E-03,1.738777E-04,-3.252845E-06/
C -K=7: E + H2 --> H2+(VIB) + 2E
C  RATE COEFFICIENT, JANEV, 2.2.9
      DATA CF7
     ./-3.568640E+01,1.733469E+01,-7.767469E+00,
     .  2.211579E+00,-4.169840E-01,5.088290E-02,
     . -3.832738E-03,1.612863E-04,-2.893392E-06/
C -K=8: E + H2+(VIB) --> H + H+ + E
C  RATE COEFFICIENT, JANEV, 2.2.12
      DATA CF8
     ./-1.781416E+01,2.277799E+00,-1.266868E+00,
     .  4.296170E-01,-9.609908E-02,1.387958E-02,
     . -1.231349E-03,6.042383E-05,-1.247521E-06/
C -K=9: E + H2+(VIB) --> H+ + H+ + 2E
C  RATE COEFFICIENT, JANEV, 2.2.11
      DATA CF9
     ./-3.746192E+01,1.559355E+01,-6.693238E+00,
     .  1.981700E+00,-4.044820E-01,5.352392E-02,
     . -4.317452E-03,1.918499E-04,-3.591779E-06/
C -K=10: E + H2+(VIB) --> H + H(N)
C  RATE COEFFICIENT, JANEV, 2.2.14
      DATA CF10
     ./-1.670436E+01,-6.035645E-01,-1.942746E-08,
     . -2.005952E-07,2.962996E-08,2.134293E-08,
     . -6.353973E-09,6.152557E-10,-2.025362E-11/
C -K=11:  FREE
C
      IF (K.LT.0) THEN
C  EIRENE DEFAULT DATA, JI=JE, LEXP=TRUE,LSUM=FALSE,LTEST=FALSE
        DO 10 J=JI,JE
          DO 10 II=1,9
            CF(II,J)=CFP(II,-K+J-JI)
10      CONTINUE
C
      ELSEIF (K.GT.0) THEN
C  DATA FROM A&M DATA FILES
!pb     IF (AL.LT.RCMN(K,2).OR.AL.GT.RCMX(K,2)) THEN
!pb       WRITE (6,*) 'WARNING FROM SUBR. CDEF '
!pb       WRITE (6,*) 'PLASMA TEMPERATURE OUT OF SAVE RANGE FOR'
!pb       WRITE (6,*) 'RATE COEFFICIENT NR. ',K, exp(al)
!pb       WRITE (6,*) 'CONSULT D.REITER '
!pb     ENDIF
        DO 11 J=JI,JE
          DO 12 II=1,9
            CF(II,J)=CREAC(II,J,K)
12        CONTINUE
11      CONTINUE
        IF (LTEST) THEN
          DO 13 J=JI,JE
            CTEST=0.
            DO 14 II=1,9
              CTEST=CTEST+ABS(CF(II,J))
14          CONTINUE
            IF (CTEST.LE.EPS60) GOTO 990
13        CONTINUE
        ENDIF
      ENDIF
C
      IF (LSUM) THEN
        DO 20 J=JI,JE
          DO 20 ICELL=1,NTE
            COU(J,ICELL)=CF(9,J)
20      CONTINUE
C
        DO 21 J=JI,JE
          DO 22 II=8,1,-1
            DO 23 ICELL=1,NTE
              COU(J,ICELL)=COU(J,ICELL)*AL(ICELL)+CF(II,J)
23          CONTINUE
22        CONTINUE
21      CONTINUE
      ENDIF
C
      IF (LEXP) THEN
        DO 24 ICELL=1,NTE
          COU(JE,ICELL)=EXP(MAX(-100._DP,COU(JE,ICELL)))
24      CONTINUE
      ENDIF
C
      RETURN
C
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN SUBROUTINE CDEF: ZERO FIT COEFFICIENTS'
      WRITE (6,*) 'J,K = ',J,K,'  EXIT CALLED!'
      CALL EXIT
      END
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
        DO 11 J=1,9
          DO 12 II=1,9
            CF(II,J)=CREAC(II,J,K)
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
          IF (AL(ICELL).LT.RCMN(K,2)) THEN
C  DETERMINE EXTRAPOLATION COEFFICIENTS FOR LINEAR EXTRAP. IN LN(<S*V>)
            S01=RCMN(K,2)
            S02=ALOG(2.)+RCMN(K,2)
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
      WRITE (6,*) 'ERROR IN SUBROUTINE CDEFN: ZERO FIT COEFFICIENTS'
      WRITE (6,*) 'J,K = ',J,K,'  EXIT CALLED!'
      CALL EXIT
      END
C
C
      SUBROUTINE CONDENSE
C
C  CONDENSE COLLISION KERNEL, IF SOME SECONDARIES ARE NOT FOLLOWED
C  BY EIRENE, I.E., IF NFOLA(IATM), NFOLM(IMOL), NFOLI(IION) LT 0
C  FOR SOME TEST PARTICLE SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ISP, IION, ICOL, IATM, IMOL, IRDS, IMDS

      DO 10 IATM=1,NATMI
C  NRCA=0 ?
        DO 100 ICOL=1,NRCA(IATM)
100     CONTINUE
10    CONTINUE
C
C
      DO 20 IMOL=1,NMOLI
C  currently: only electron impact collisions
        DO 200 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 220 IION=1,NIONI
            ISP=NSPAM+IION
            IF (PIODS(IRDS,IION).GT.0) THEN
              IF (NFOLI(IION).LT.0) THEN
                WRITE (6,*) 'TEST ION ',TEXTS(ISP),' CAN BE CONDENSED'
              ENDIF
            ENDIF
220       CONTINUE
200     CONTINUE
20    CONTINUE
C
      RETURN
      END
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
     .          CCXM1, CCXM2, EXPO, EXTRAP, E, XI
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
              S02=ALOG(2.)+RCMN(K,1)
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
          CALL EXIT
        END IF
      ELSE
        WRITE (6,*) 'ERROR IN CROSS: K= ',K
        WRITE (6,*) 'CALLED FROM ',TEXT
        WRITE (6,*) 'REACTION NO. ',IR
c slmod begin (sl)
        WRITE (0,*) 'ERROR IN CROSS: K= ',K
c slmod end
      ENDIF
      RETURN
      END



c     ***********************************************************
c     --Subroutinen--


      subroutine energyloc(xx,a,x,b)
c     ***********************************************************
c     * Ermittlung des Feldindex b mit vogegebener Energie x,   *
c     * so dass x zwischen Energie(b) und Energie(b+1)          *
c     * Modifikation: Ergebnis: Energien zw. 0.1 eV und 100 eV  *
c     * *********************************************************
      USE PRECISION
      implicit none
      integer a,b
      REAL(DP) :: x,xx(a)
      integer bl,bm,bu

      bl=0
      bu=a+1
      
80    if(bu-bl.gt.1)then
        bm=(bu+bl)*0.5
        if((xx(a).ge.xx(1)).eqv.(x.ge.xx(bm)))then
          bl=bm
        else
          bu=bm
        endif
      goto 80
      endif

      if(x.le.xx(1))then
        b=1
      else if(x.eq.xx(a))then
        b=a
      else
        b=bl
      endif

      return
      end
C
C
      FUNCTION EXTRAP(ELAB,IFLAG,FP1,FP2,FP3)
C
C  NOTE:
C  INPUT:  ELAB IS LOG OF RELATIVE ENERGY, OR LOG OF TEMP
C  OUTPUT: EXTRAP IS NOT LOG, BUT THE TRUE VALUE
C
C  FUNCTION FOR EXTRAPOLATING SINGLE PARAMETER FITS BEYOND THEIR
C  RANGE OF VALIDITY
C  TYPE  IFLAG=1--4: JANEV ET AL. , SPRINGER, 1987, P13
C  TYPE  IFLAG=5  BACHMANN ET AL., IPP-REPORT, .....ELASTIC
C
      USE PRECISION
      
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ELAB, FP1, FP2, FP3
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: X, EL, EXTRAP

      IF (IFLAG.EQ.1) THEN
C  NON ZERO THRESHOLD
        EXTRAP=0.
      ELSEIF (IFLAG.EQ.2) THEN
C  EXTRAPOLATION AT HIGH ENERGY END FOR REACTIONS WITH NON ZERO THRESHOLD
C  FP1 SHOULD BE = E_THRESHOLD (EV)
        EL=EXP(ELAB)
        X=EL/FP1
        EXTRAP=FP2*X**FP3*LOG(X)
      ELSEIF (IFLAG.EQ.3) THEN
        EXTRAP=EXP(FP1+FP2*ELAB)
      ELSEIF (IFLAG.EQ.4) THEN
C
C  OUT
C       EXTRAP=EXP((FP1+FP2*ELAB)**2)
C
      ELSEIF (IFLAG.EQ.5) THEN
C  LINEAR OR QUADRATIC EXTRAPOLATION IN LN(SIGMA)
        EXTRAP=EXP(FP1+FP2*ELAB+FP3*ELAB**2)
      ELSE
        GOTO 999
      ENDIF
      RETURN
999   CONTINUE
      WRITE (6,*) 'ERROR IN EXTRAP. EXIT CALLED '
      CALL EXIT
      END


      FUNCTION FEELEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) ELEIC(9), FEELEI1, PLS, DEIMIN, EE, FTABEI1, DSUB, ELEI
      INTEGER :: J, I, KK, II

      FEELEI1=0.D0
      KK=NELREI(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEELEI1=-EIONHE
        CASE (-2)
            FEELEI1=EELDS1(IREI,1)
        CASE (-3)
            FEELEI1=-1.5*TEIN(K)
        CASE (-4)
            FEELEI1=-EIONH
        CASE (-5)
            FEELEI1=-10.5
        CASE (-6)
            FEELEI1=-25.0
        CASE (-7)
            FEELEI1=EELDS1(IREI,1)
        CASE (-8)
            FEELEI1=-10.5
        CASE (-9)
            FEELEI1=-15.5
        CASE (-10)
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            FEELEI1=-0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        IF (JELREI(IREI) == 1) THEN
          ELEI = CREAC(9,1,KK)
          DO II=8,1,-1
            ELEI = ELEI*TEINL(K) + CREAC(II,1,KK)
          END DO
          ELEI=EXP(MAX(-100._DP,ELEI))
          FEELEI1=-ELEI*DEIN(K)*FACREA(KK)/(FTABEI1(IREI,K)+EPS60)
        ELSE
          ELEIC(1:JELREI(IREI)) = CREAC(9,1:JELREI(IREI),KK)
          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          DO J=1,JELREI(IREI)
            DO II=8,1,-1
              ELEIC(J)=ELEIC(J)*TEINL(K)+CREAC(II,J,KK)
            END DO
          END DO
          ELEI = ELEIC(9)
          DO I=8,1,-1
            ELEI=ELEI*PLS+ELEIC(I)
          END DO
          EE=MAX(-100._DP,ELEI+FACREA(KK)+DEINL(K))
          FEELEI1=-EXP(EE)/(FTABEI1(IREI,K)+EPS60)
        END IF
      END IF

      RETURN
      END


      FUNCTION FEELRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: ELRC1(9), PLS, DELE, EE, FEELRC1, FTABRC1, DSUB, 
     .            DEIMIN, ELRC
      INTEGER :: J, I, KK, II

      FEELRC1=0.D0
      KK=NELRRC(IRRC)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        CASE (-2)
            FEELRC1=EELRC1(IRRC,1)*FTABRC1(IRRC,K)
        CASE (-3)
            FEELRC1=-1.5*TEIN(K)*FTABRC1(IRRC,K)
        END SELECT
      ELSE IF (KK > 0) THEN
        IF (JELRRC(IRRC) == 1) THEN
          ELRC = CREAC(9,1,KK)
          DO II=8,1,-1
            ELRC = ELRC*TEINL(K) + CREAC(II,1,KK)
          END DO
          ELRC=ELRC+FACREA(KK)
          ELRC=EXP(MAX(-100._DP,ELRC))
          FEELRC1=-ELRC*DEIN(K)
        ELSE
          ELRC1(1:JELRRC(IRRC)) = CREAC(9,1:JELRRC(IRRC),KK)
          DSUB=LOG(1.D8)
          DEIMIN=LOG(1.D8)
          PLS=MAX(DEIMIN,DEINL(K))-DSUB
          DO J=1,JELRRC(IRRC)
            DO II=8,1,-1
              ELRC1(J)=ELRC1(J)*TEINL(K)+CREAC(II,J,KK)
            END DO
          END DO
          ELRC = ELRC1(9)
          DO I=8,1,-1
            ELRC=ELRC*PLS+ELRC1(I)
          END DO
          EE=MAX(-100._DP,ELRC+DEINL(K)+FACREA(KK))
          FEELRC1=-EXP(EE)
        END IF
        IF (DELPOT(KK).NE.0.D0) THEN
          DELE=DELPOT(KK)
          FEELRC1=FEELRC1+DELE*FTABRC1(IRRC,K)
        END IF
      END IF

      RETURN
      END


      FUNCTION FEHVDS1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: FEHVDS1, EHVDS, FTABEI1
      INTEGER :: II, KK

      FEHVDS1=0.D0
      KK=NREAHV(IREI)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
            FEHVDS1=EHVDS1(IREI,1)
        CASE (-2)
            FEHVDS1=6.
        CASE (-3)
            FEHVDS1=10.0
        CASE (-4)
            FEHVDS1=8.6
        CASE (-5)
            FEHVDS1=0.5
        CASE (-6)
            FEHVDS1=0.88*TEIN(K)
        END SELECT
      ELSE IF (KK > 0) THEN
        EHVDS = CREAC(9,1,KK)
        DO II=8,1,-1
          EHVDS = EHVDS*TEINL(K) + CREAC(II,1,KK)
        END DO
        EHVDS=EXP(MAX(-100._DP,EHVDS+FACREA(KK)))
        FEHVDS1=EHVDS*DEIN(K)/(FTABEI1(IREI,K)+EPS60)
      END IF

      RETURN
      END


      FUNCTION FEPLCX3 (IRCX,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRCX, K
      REAL(DP) :: PLS, ADD, ELCX, FEPLCX3
      INTEGER :: II, KK

      FEPLCX3=0.D0
      KK=NELRCX(IRCX)
      IF (KK < 0) THEN
        SELECT CASE (KK)
C  DEFAULT CX MODEL
        CASE (-1)
            FEPLCX3=1.5*TIIN(IPLS,K)+EDRIFT(IPLS,K)
        CASE (-2)
C  MEAN ENERGY FROM DRIFTING MONOENERGETIC
            FEPLCX3=EPLCX3(IRCX,1,1)+EDRIFT(IPLS,K)
        CASE (-3)
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
            FEPLCX3=1.5*TIIN(IPLS,K)+EDRIFT(IPLS,K)
        END SELECT
      ELSE
C  MEAN ENERGY FROM SINGLE PARAMETER FIT KK
        ELCX = CREAC(9,1,KK)
        PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
        DO II=8,1,-1
          ELCX = ELCX*PLS + CREAC(II,1,KK)
        END DO
        ADD=EPLCX3(IRCX,1,1)
        FEPLCX3=ELCX*DIIN(IPLS,K)*ADD
      END IF

      RETURN
      END


      FUNCTION FEPLEL3 (IREL,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREL, K
      REAL(DP) :: PLS, ADD, EPEL, FEPLEL3
      INTEGER :: II, KK

      FEPLEL3=0.D0
      KK=NELREL(IREL)
      IF (KK < 0) THEN
        SELECT CASE (KK)
        CASE (-1)
C  DEFAULT EL MODEL
C  OUT
        CASE (-2)
C  MEAN ENERGY FROM DRIFTING MONOENERGETIC
            FEPLEL3=EPLEL3(IREL,1,1)+EDRIFT(IPLS,K)
        CASE (-3)
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
            FEPLEL3=1.5*TIIN(IPLS,K)+EDRIFT(IPLS,K)
        END SELECT
      ELSE
C  MEAN ENERGY FROM SINGLE PARAMETER FIT KK
        EPEL = CREAC(9,1,KK)
        PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
        DO II=8,1,-1
          EPEL = EPEL*PLS + CREAC(II,1,KK)
        END DO
        ADD=EPLEL3(IREL,1,1)
        FEPLEL3=EPEL*DIIN(IPLS,K)*ADD
      END IF

      RETURN
      END


      FUNCTION FEPLPI3 (IRPI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRPI, K
      REAL(DP) :: FEPLPI3
      INTEGER :: KK

      FEPLPI3=0.D0
      KK=NELRPI(IRPI)
      SELECT CASE (KK)
      CASE (-1)
          FEPLPI3=EPLPI3(IRPI,1,1)
      CASE (-2)
          FEPLPI3=1.5*TIIN(IPLS,K)+EDRIFT(IPLS,K)
      END SELECT

      RETURN
      END
C
C
      FUNCTION FI(R,ER,B,IFLAG,P,DFI)

C  EVALUATE EFFECTIVE POTENTIAL FUNCTION FI AT R
C  EVALUATE DFI(R)/DR AT R
C  RETURN FI=FI(R), DFI=DFI(R)/DR
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES,  H+ + H2,  HE+ + HE
C

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: R, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: DFI
      INTEGER, INTENT(IN) :: IFLAG

      REAL(DP) :: DSS, U, DU, SS, RI, RIQ, G1, REFF, RLIM, RR1, RR, 
     .            EMRR, G2, EX2, B2, V, DV, EX, RLOW, R2, R3, FI, G, 
     .            RMI, EPS

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        R2=R*R
        R3=R2*R
C  FIND V=V(R) and DV=DV(R)/DR
        IF (R.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSEIF (R.LT.RLOW) THEN
          R2=RLOW*RLOW
          R3=R2*RLOW
          EX=EXP(-RLOW)
          EX2=EX*EX
          SS=(1.+RLOW+R2/3.)*EX-1.
          RI=1./RLOW
          V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          DV=0.D0
        ELSE
          EX=EXP(-R)
          EX2=EX*EX
          SS=(1.+R+R2/3.)*EX-1.
          RI=1./R
          V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          RIQ=RI*RI
          DSS=-R/3.*(1.+R)*EX
          U=RI-(1.+RI)*EX2-(1.+R)*EX
          DU=-RIQ+(RIQ+2.*RI+2.)*EX2+R*EX
          DV=27.211*((DU*SS-U*DSS)/(SS*SS)-RIQ)
        ENDIF
C  FIND FI=FI(R) AND DFI=DFI(R)/DR
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSEIF (IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        R2=R*R
        R3=R2*R
C
        RR=R*RMI
        EMRR=1.-RR
C
C       G2=1.00+(1.0-G2)*MAX(0.D0,-EMRR)/EMRR
C       REFF=-G1*G2*EMRR
        IF (RR.LT.1.0) THEN
          G=G1
        ELSE
          G=G1*G2
        ENDIF
        REFF=-G*EMRR
C
        IF (REFF.GT.160.D0) THEN
          V=0.D0
          DV=0.D0
        ELSE
          EX=EXP(-REFF)
          EX2=EX*EX
          V=EPS*(EX2-(EX+EX))
          DV=-2.*EPS*RMI*G*(EX2-EX)
        ENDIF
        FI=1.-V/ER-B2/R2
        DFI=-DV/ER+2.*B2/R3
C
      ELSE
        WRITE (6,*) 'ERROR IN FUNCTION FI. IFLAG INVALID.'
        WRITE (6,*) 'IFLAG = ',IFLAG
        CALL EXIT
      ENDIF
C
      RETURN
      END
C
C
      SUBROUTINE FIVEC(ER,B,IFLAG,P)

C  VECTORIZED VERSION OF FUNCTION FI FOR GAUSS MEHLER QUADRATURE
C  EVALUATE EFFECTIVE POTENTIAL FUNCTION AT AR(I),I=1,NFI
C  NOTE: NFI.LE.128 IS NOT CHECKED, BUT USED
C  RETURN FI(AR(I)) IN THE ARRAY AFI(I),I=1,NFI
C     --------------
C  IFLAG=1:  H+ + H
C  IFLAG=2:  H+ + NOBLE GASES , H+ + H2, HE+ + HE
C

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      REAL(DP) :: RI, G1, SS, EX, EX2, RR, RLIM, REFF, EMRR, G2, V, 
     .          R, R2, B2, RLOW, RR1, G, RMI, EPS
      INTEGER :: IFI

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI
      COMMON /CFI/ AR,AFI,NFI

      B2=B*B
C
      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.D-7
        DO 1 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
C  FIND V=V(R)
          IF (R.GT.160.D0) THEN
            V=0.D0
          ELSEIF (R.LT.RLOW) THEN
            R2=RLOW*RLOW
            EX=EXP(-RLOW)
            EX2=EX*EX
            SS=(1.+RLOW+R2/3.)*EX-1.
            RI=1./RLOW
            V=27.211*((RI-(1.+RI)*EX2-(1.+RLOW)*EX)/SS+RI)
          ELSE
            EX=EXP(-R)
            EX2=EX*EX
            SS=(1.+R+R2/3.)*EX-1.
            RI=1./R
            V=27.211*((RI-(1.+RI)*EX2-(1.+R)*EX)/SS+RI)
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
1       CONTINUE
C
      ELSEIF(IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES, (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1./P(4)
        DO 2 IFI=1,NFI
          R=AR(IFI)
          R2=R*R
          RR=R*RMI
          EMRR=1.-RR
          IF (RR.LT.1) THEN
            G=G1
          ELSE
            G=G1*G2
          ENDIF
          REFF=-G*EMRR
C
          IF (REFF.GT.160.D0) THEN
            V=0.D0
          ELSE
            EX=EXP(-REFF)
            EX2=EX*EX
            V=EPS*(EX2-(EX+EX))
          ENDIF
          AFI(IFI)=1.-V/ER-B2/R2
2       CONTINUE
C
      ELSE
        WRITE (6,*) 'ERROR IN FUNCTION FIVEC. IFLAG INVALID.'
        WRITE (6,*) 'IFLAG = ',IFLAG
        CALL EXIT
      ENDIF
C
      RETURN
      END
C
C
C       FUNCTION FPATHA(K)
C       FUNCTION FPATHM(K)
C       FUNCTION FPATHI(K)
C       FUNCTION EXTRAP(....)
C       SUBROUTINE CDEF(.....)
C       SUBROUTINE SETAMD
C       SUBROUTINE XSECTA
C       SUBROUTINE XSECTM
C       SUBROUTINE XSECTI
C       SUBROUTINE XSECTP
C
C
      FUNCTION FPATHA (K,CFLAG)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR NEUTRAL
C   "BEAM ATOMS" OF VELOCITY VEL IN DRIFTING MAXWELLIAN PLASMA-BACKGROUND
C   IN CELL K
C   OUTPUT: COMMON COMLCA
C           CFLAG: FLAG FOR SAMPLING OF POST COLLISION STATES
C           CFLAG(1,...): EI
C           CFLAG(3,...): CX
C           CFLAG(4,...): II
C           CFLAG(5,...): EL
C           CFLAG(6,...): RC
C
C  CFLAG(...,1):
C      =0:   VI: DELTA COLLISION IN VELOCITY SPACE (BUT DIFFERENT
C                                                   SPECIES ALLOWED)
C      =1:   VI: MONOENERGETIC AND ISOTROPIC IN CENTER OF MASS SYSTEM
C      =2:   VI: MAXWELL PLUS DRIFT
C      =3:   VI: SIGMA-V-WEIGHTED MAXWELLIAN PLUS DRIFT
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(6,3)
      INTEGER, INTENT(IN) :: K

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLS)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)
      REAL(DP) :: CEL, RMN, RLMS, ER, RMI, RMSI, CXS, VEFFQ, TBCX, VEFF,
     .          FEPLCX3, AU_TO_CM2 , SIGHABER, FEPLEL3, TBEL,
     .          ELTHDUM, CTCHDUM, SIGMAX, FTABEI1, EHEAVY, FEELEI1,
     .          FEHVDS1, DENEL, FPATHA, VX, VY, VZ, PVELQ0, ELAB,
     .          VRELQ, VREL, FEPLPI3, CII, CROSS, ELB, PLS, TBPI, EXPO
      INTEGER :: IBGK, IRCX, IAEL, IREL, IPL, IAT, IAEI, IRDS, IAPI,
     .           IRPI, KK, IREAC, IACX, II, IF8, JAN, J
C
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTOR=0.D0
      XSTORV=0.D0
      FPATHA=1.D10
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
C
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSI
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT COLLISION - RATE - COEFFICIENT
C  NO MASS SCALING NEEDED FOR BULK ELECTRONS
C
20    IF (LGAEI(IATM,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 30
      DO 10 IAEI=1,NAEII(IATM)
        IRDS=LGAEI(IATM,IAEI)
        IF (MODCOL(1,2,NSPH+IATM,1).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IRDS)=TABDS1(IRDS,K)
          ELSE
            SIGVEI(IRDS)=FTABEI1(IRDS,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IRDS,5)=EELDS1(IRDS,K)
          EHEAVY=EHVDS1(IRDS,K)
        ELSE
          ESIGEI(IRDS,5)=FEELEI1(IRDS,K)
          EHEAVY=FEHVDS1(IRDS,K)
        END IF
C
        ESIGEI(IRDS,1)=EATDS(IRDS,0,1)*E0+EATDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,2)=EMLDS(IRDS,0,1)*E0+EMLDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,3)=EIODS(IRDS,0,1)*E0+EIODS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,4)=EPLDS(IRDS,  1)*E0+EPLDS(IRDS,  2)*EHEAVY

        SIGEIT=SIGEIT+SIGVEI(IRDS)
10    CONTINUE
C
C  ION IMPACT ON ATOM IATM, ION SPEZIES IPLS=1,NPLSI
C  30--->40
C
30    IF (LGAPI(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 40
      DO 36 IAPI=1,NAPII(IATM)
        IRPI=LGAPI(IATM,IAPI,0)
        IPLS=LGAPI(IATM,IAPI,1)
        IF (LGVAC(K,IPLS)) GOTO 36
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(4,2,NSPH+IATM,IPLS).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVPI(IRPI)=TABPI3(IRPI,K,1)
          ELSE
            KK=NREAPI(IRPI)
            PLS=TIINL(IPLS,K)+ADDPI(IRPI,IPLS)
            TBPI = CREAC(9,1,KK)
            DO II=8,1,-1
              TBPI = TBPI*PLS + CREAC(II,1,KK)
            END DO
            TBPI=EXP(MAX(-100._DP,TBPI))*DIIN(IPLS,K)
            SIGVPI(IRPI)=TBPI
          END IF
        ELSEIF (MODCOL(4,2,NSPH+IATM,IPLS).EQ.2) THEN
C  BEAM - MAXWELL
C
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFPI(IRPI))
          IF (NSTORDR >= NRAD) THEN
            TBPI3(1:NSTORDT) = TABPI3(IRPI,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NREAPI(IRPI)
            PLS=TIINL(IPLS,K)+ADDPI(IRPI,IPLS)
            DO J=JAN+1,9
              TBPI3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                TBPI3(J)=TBPI3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) TBPI3(1)=TBPI3(1)+DIINL(IPLS,K)
          EXPO=TBPI3(9)
          DO 33 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+TBPI3(IF8)
33        CONTINUE
          SIGVPI(IRPI)=EXP(EXPO)
        ELSEIF (MODCOL(4,2,NSPH+IATM,IPLS).EQ.3) THEN
C  BEAM - BEAM
          VRELQ=ZTI(IPLS)+PVELQ(IPLS)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFPI(IRPI)
          IREAC=MODCOL(4,1,NSPH+IATM,IPLS)
          CII=CROSS(ELAB,IREAC,IRPI,'FPATHA II')
          SIGVPI(IRPI)=CII*VREL*DENIO(IPLS)
        ELSE
          GOTO 991
        ENDIF
C
C  2.A ELECTRON ENERGY LOSS / COLLISION (EV)
C
        ESIGPI(IRPI,2)=0.D0
C
        SIGPIT=SIGPIT+SIGVPI(IRPI)
C
C  2.B BULK ION ENERGY LOSS / COLLISION (EV)
C
        IF (NSTORDR >= NRAD) THEN
          ESIGPI(IRPI,1)=EPLPI3(IRPI,K,1)
        ELSE
          ESIGPI(IRPI,1)=FEPLPI3(IRPI,K)
        END IF
        CFLAG(4,1)=2
36    CONTINUE
C
C
C  CHARGE EXCHANGE OF ATOMS IATM  WITH IONS OF SPEZIES IPLS=1,NPLSI
C  40--->50
C
40    CONTINUE
      IF (LGACX(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
      DO 41 IACX=1,NACXI(IATM)
        IRCX=LGACX(IATM,IACX,0)
        IPLS=LGACX(IATM,IACX,1)
        IF (LGVAC(K,IPLS)) GOTO 41
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(3,2,NSPH+IATM,IPLS).EQ.1) THEN
C  MODEL 1:
C  MAXWELLIAN RATE. IGNORE NEUTRAL VELOCITY
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
            KK=NREACX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            TBCX = CREAC(9,1,KK)
            DO II=8,1,-1
              TBCX = TBCX*PLS + CREAC(II,1,KK)
            END DO
            TBCX=EXP(MAX(-100._DP,TBCX))*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,NSPH+IATM,IPLS).EQ.2) THEN
C  MODEL 2:
C  BEAM - MAXWELLIAN RATE
          IF (TIIN(IPLS,K).LT.TVAC) THEN
C     HERE: T_I IS SO LOW, THAT ALL ION ENERGY IS IN DRIFT MOTION.
C           HENCE: USE BEAM-BEAM RATE INSTEAD.
            VRELQ=PVELQ(IPLS)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,NSPH+IATM,IPLS)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX1 ')
            SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              JAN=NSTORDT
            ELSE
              JAN=0
            END IF
            IF (JAN < 9) THEN
              KK=NREACX(IRCX)
              PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
              DO J=JAN+1,9
                TBCX3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                  TBCX3(J)=TBCX3(J)*PLS+CREAC(II,J,KK)
                END DO
              END DO
            END IF
            IF (JAN < 1) TBCX3(1)=TBCX3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
            EXPO=TBCX3(9)
            DO 43 II=1,8
              IF8=9-II
              EXPO=EXPO*ELB+TBCX3(IF8)
43          CONTINUE
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,NSPH+IATM,IPLS).EQ.3) THEN
C  MODEL 3
C  BEAM - BEAM RATE,  BUT WITH EFFECTIVE INTERACTION ENERGY
          VEFFQ=ZTI(IPLS)+PVELQ(IPLS)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,NSPH+IATM,IPLS)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX2')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSEIF (MODCOL(3,2,NSPH+IATM,IPLS).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLS)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,NSPH+IATM,IPLS)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHA CX3')
          SIGVCX(IRCX)=CXS*VREL*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,NSPH+IATM,IPLS).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=2
        ELSEIF (MODCOL(3,4,NSPH+IATM,IPLS).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NELRCX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            DO J=JAN+1,9
              EPCX3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                EPCX3(J)=EPCX3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) EPCX3(1)=EPCX3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
          EXPO=EPCX3(9)
          DO 45 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+EPCX3(IF8)
45        CONTINUE
          ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
          ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          CFLAG(3,1)=3
        ELSEIF (MODCOL(3,4,NSPH+IATM,IPLS).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C
C  ELASTIC COLLISIONS OF ATOMS IATM  WITH IONS OF SPEZIES IPLS=1,NPLSI
C  50--->60
C
50    CONTINUE
      IF (LGAEL(IATM,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 60
      DO 51 IAEL=1,NAELI(IATM)
        IREL=LGAEL(IATM,IAEL,0)
        IPLS=LGAEL(IATM,IAEL,1)
        IBGK=NPBGKP(IPLS,1)
        IF (LGVAC(K,IPLS)) GOTO 51
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(5,2,NSPH+IATM,IPLS).EQ.1) THEN
C  MAXWELLIAN DISTRIBUTION OF RELATIVE SPEEDS
          IF (NSTORDR >= NRAD) THEN
            SIGVEL(IREL)=TABEL3(IREL,K,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
            TBEL = CREAC(9,1,KK)
            DO II=8,1,-1
              TBEL = TBEL*PLS + CREAC(II,1,KK)
            END DO
            TBEL=EXP(MAX(-100._DP,TBEL))*DIIN(IPLS,K)
            SIGVEL(IREL)=TBEL
          END IF
        ELSEIF (MODCOL(5,2,IATM,IPLS).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLS,K).LT.TVAC) THEN
C  TEMPERATURE TOO LOW, USE: BEAM_ATOM - BEAM_DRIFT RATECOEFF.
            VRELQ=PVELQ(IPLS)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,NSPH+IATM,IPLS)
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL1')
            SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBEL3(1:NSTORDT) = TABEL3(IREL,K,1:NSTORDT)
              JAN=NSTORDT
            ELSE
              JAN=0
            END IF
            IF (JAN < 9) THEN
              KK=NREAEL(IREL)
              PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
              DO J=JAN+1,9
                TBEL3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                  TBEL3(J)=TBEL3(J)*PLS+CREAC(II,J,KK)
                END DO
              END DO
            END IF
            IF (JAN < 1) TBEL3(1)=TBEL3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFEL(IREL))
            EXPO=TBEL3(9)
            DO 53 II=1,8
              IF8=9-II
              EXPO=EXPO*ELB+TBEL3(IF8)
53          CONTINUE
            SIGVEL(IREL)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(5,2,NSPH+IATM,IPLS).EQ.3) THEN
C  BEAM - BEAM APPROXIMATION: USE THERMAL SPEED OF IPLS AS BEAM-SPEED
          VRELQ=ZTI(IPLS)+PVELQ(IPLS)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,NSPH+IATM,IPLS)
CH ANALOG HABERSCHEID KRAM FUER SIGMA
          IF (LHABER) THEN
            RMN=RMASSA(IATM)
            RMI=RMASSP(IPLS)
            RMSI=1./(RMN+RMI)
            RLMS=RMN*RMI*RMSI
            ER=RLMS*VRELQ*CVELI2
            AU_TO_CM2=5.29177E-9**2
            CALL SCATANG (ER,0.5_DP,ELTHDUM,CTCHDUM,SIGHABER)
            CEL= SIGHABER*AU_TO_CM2
          ELSE
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL2')
          END IF
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSEIF (MODCOL(5,2,NSPH+IATM,IPLS).EQ.4) THEN
C  MODEL 4
C  BEAM - BEAM RATE, IGNORE THERMAL ION ENERGY
          VRELQ=PVELQ(IPLS)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,NSPH+IATM,IPLS)
          CEL=CROSS(ELAB,IREAC,IREL,'FPATHA EL3')
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSE
          GOTO 995
        ENDIF
        SIGELT=SIGELT+SIGVEL(IREL)
        IF (IBGK.NE.0) SIGBGK=SIGBGK+SIGVEL(IREL)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(5,4,NSPH+IATM,IPLS).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLEL3(IREL,K)
          END IF
          CFLAG(5,1)=2
        ELSEIF (MODCOL(5,4,NSPH+IATM,IPLS).EQ.2) THEN
          IF (NSTORDR >= NRAD) THEN
            EPEL3(1:NSTORDT) = EPLEL3(IREL,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NELREL(IREL)
            PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
            DO J=JAN+1,9
              EPEL3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                EPEL3(J)=EPEL3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) EPEL3(1)=EPEL3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFEL(IREL))
          EXPO=EPEL3(9)
          DO 55 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+EPEL3(IF8)
55        CONTINUE
          ESIGEL(IREL,1)=EXP(EXPO)/SIGVEL(IREL)
          ESIGEL(IREL,1)=ESIGEL(IREL,1)+EDRIFT(IPLS,K)
          CFLAG(5,1)=3
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ISOTROPIC ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLCX3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLCX3(IREL,K)
          END IF
          CFLAG(5,1)=1
        ELSE
          GOTO 995
        ENDIF
51    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
100   CONTINUE
C
      SIGMAX=MAXVAL(XSTOR(:,1:4))
      WHERE (XSTOR(:,1:4) .LE. SIGMAX*1.D-10 )
        XSTOR(:,1:4) = 0.D0
      END WHERE
C
      SIGTOT=SIGEIT+SIGPIT+SIGCXT+SIGELT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHA=VEL/SIGTOT
        ZMFPI=1./FPATHA
      ENDIF
C
      RETURN
990   CONTINUE
      IAT=NSPH+IATM
      WRITE (6,*) 'ERROR IN FPATHA: INCONSISTENT ELEC. IMP. IONIZ. DATA'
      WRITE (6,*) 'IAT,MODCOL(1,J,IAT,1) '
      WRITE (6,*) IATM,(MODCOL(1,J,IAT,1),J=1,4)
      CALL EXIT
991   CONTINUE
      IAT=IATM
      WRITE (6,*) 'ERROR IN FPATHA: INCONSISTENT ION IMP. IONIZ. DATA'
      WRITE (6,*) 'IAT,IPL,MODCOL(4,J,IAT,IPL) '
      DO 993 IPL=1,NPLSI
        WRITE (6,*) IATM,IPL,(MODCOL(4,J,NSPH+IATM,IPL),J=1,4)
993   CONTINUE
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN FPATHA: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (6,*) 'IAT,IPL,MODCOL(3,J,IAT,IPL) '
      DO 994 IPL=1,NPLSI
        WRITE (6,*) IATM,IPL,(MODCOL(3,J,NSPH+IATM,IPL),J=1,4)
994   CONTINUE
      CALL EXIT
995   CONTINUE
      WRITE (6,*) 'ERROR IN FPATHA: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (6,*) 'IAT,IPL,MODCOL(5,J,IAT,IPL) '
      DO 996 IPL=1,NPLSI
        WRITE (6,*) IATM,IPL,(MODCOL(5,J,NSPH+IATM,IPL),J=1,4)
996   CONTINUE
      CALL EXIT
      END
C
      FUNCTION FPATHI (K,CFLAG)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR
C   "BEAM TEST IONS" OF VELOCITY VEL IN MAXWELLIAN PLASMA-BACKGROUND
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(6,3)
      INTEGER, INTENT(IN) :: K

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLS)
      REAL(DP) :: TBCX3(9), TBPI3(9), TBEL3(9)
      REAL(DP) :: EPCX3(9), EPPI3(9), EPEL3(9)
      REAL(DP) :: ELB, EXPO, FEPLCX3, ELAB, VEFF, CROSS, CXS, SIGMAX,
     .          VX, VY, VZ, PVELQ0, FPATHI, TEEI, DENEL, FTABEI1, PLS,
     .          TBCX, VEFFQ, FEELEI1, EHEAVY, FEHVDS1
      INTEGER :: J, IF8, JAN, IREAC, IPL, IIO, IRDS, IIDS, KK, II,
     .           IRCX, IICX, IIPI
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTOR=0.D0
      XSTORV=0.D0
      FPATHI=1.D10
C
      IF (LGVAC(K,0)) RETURN
C
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
      TEEI=TEIN(K)
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSI
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C
C  TEST ION ELECTR. IMP. RATE COEFFICIENT
C
      IF (LGIEI(IION,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 25
      DO 10 IIDS=1,NIDSI(IION)
        IRDS=LGIEI(IION,IIDS)
        IF (MODCOL(1,2,NSPAM+IION,1).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IRDS)=TABDS1(IRDS,K)
          ELSE
            SIGVEI(IRDS)=FTABEI1(IRDS,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IRDS,5)=EELDS1(IRDS,K)
          EHEAVY=EHVDS1(IRDS,K)
        ELSE
          ESIGEI(IRDS,5)=FEELEI1(IRDS,K)
          EHEAVY=FEHVDS1(IRDS,K)
        END IF
C
        ESIGEI(IRDS,1)=EATDS(IRDS,0,1)*E0+EATDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,2)=EMLDS(IRDS,0,1)*E0+EMLDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,3)=EIODS(IRDS,0,1)*E0+EIODS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,4)=EPLDS(IRDS,  1)*E0+EPLDS(IRDS,  2)*EHEAVY
C
        SIGEIT=SIGEIT+SIGVEI(IRDS)
10    CONTINUE
C
C  IONIZATION OF TEST ION BY BULK ION IMPACT, ION SPEZIES IPLS=1,NPLSI
C
25    CONTINUE
      DO 30 IIPI=1,NIPII(IION)
        SIGVPI(IIPI)=0.
        SIGPIT=SIGPIT+SIGVPI(IIPI)
30    CONTINUE
C
C   CHARGE EXCHANGE RATE COEFFICIENT FOR TEST ION IION
C   WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C
      IF (LGICX(IION,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
40    CONTINUE
      DO 41 IICX=1,NICXI(IION)
        IRCX=LGICX(IION,IICX,0)
        IPLS=LGICX(IION,IICX,1)
        IF (LGVAC(K,IPLS)) GOTO 41
C
        IF (MODCOL(3,2,NSPAM+IION,IPLS).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
            KK=NREACX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            TBCX = CREAC(9,1,KK)
            DO II=8,1,-1
              TBCX = TBCX*PLS + CREAC(II,1,KK)
            END DO
            TBCX=EXP(MAX(-100._DP,TBCX))*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,NSPAM+IION,IPLS).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLS,K).LT.TVAC) THEN
            VEFFQ=PVELQ(IPLS)
            VEFF=SQRT(VEFFQ)
            ELAB=LOG(VEFFQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,NSPAM+IION,IPLS)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHI CX')
            SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              JAN=NSTORDT
            ELSE
              JAN=0
            END IF
            IF (JAN < 9) THEN
              KK=NREACX(IRCX)
              PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
              DO J=JAN+1,9
                TBCX3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                  TBCX3(J)=TBCX3(J)*PLS+CREAC(II,J,KK)
                END DO
              END DO
            END IF
            IF (JAN < 1) TBCX3(1)=TBCX3(1)+DIINL(IPLS,K)
C  MINIMUM ENERGY: 0.1EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
            EXPO=TBCX3(9)
            DO 43 II=1,8
              IF8=9-II
              EXPO=EXPO*ELB+TBCX3(IF8)
43          CONTINUE
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,NSPAM+IION,IPLS).EQ.3) THEN
C  BEAM - BEAM
          VEFFQ=ZTI(IPLS)+PVELQ(IPLS)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,NSPAM+IION,IPLS)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHI CX')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
        IF (MODCOL(3,4,NSPAM+IION,IPLS).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=2
        ELSEIF (MODCOL(3,4,NSPAM+IION,IPLS).EQ.2) THEN
          IF (NSTORDR >= NRAD) THEN
            EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NELRCX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            DO J=JAN+1,9
              EPCX3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                EPCX3(J)=EPCX3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) EPCX3(1)=EPCX3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
          EXPO=EPCX3(9)
          DO 45 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+EPCX3(IF8)
45        CONTINUE
          ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
          ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          CFLAG(3,1)=3
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C
50    CONTINUE
C
C     TOTAL
C
      SIGMAX=MAXVAL(XSTOR(:,1:4))
      WHERE (XSTOR(:,1:4) .LE. SIGMAX*1.D-10 )
        XSTOR(:,1:4) = 0.D0
      END WHERE
C
100   CONTINUE
      SIGTOT=SIGPIT+SIGCXT+SIGEIT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHI=VEL/SIGTOT
        ZMFPI=1./FPATHI
      ENDIF
C
      RETURN
990   CONTINUE
      IIO=NSPAM+IION
      WRITE (6,*) 'ERROR IN FPATHI: INCONSISTENT ELEC. IMP. DATA'
      WRITE (6,*) 'IIO,MODCOL(1,J,IIO,1) '
      WRITE (6,*) IION,(MODCOL(1,J,IIO,1),J=1,4)
      CALL EXIT
992   CONTINUE
      IIO=NSPAM+IION
      WRITE (6,*) 'ERROR IN FPATHI: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (6,*) 'IIO,IPL,MODCOL(3,J,IIO,IPL) '
      DO 994 IPL=1,NPLSI
        WRITE (6,*) IION,IPL,(MODCOL(3,J,IIO,IPL),J=1,4)
994   CONTINUE
      CALL EXIT
      END
C
C
      FUNCTION FPATHM (K,CFLAG)
C
C   CALCULATE MEAN FREE PATH AND REACTION RATES FOR NEUTRAL
C   "BEAM MOLECULES" OF VELOCITY VEL IN MAXWELLIAN PLASMA-BACKGROUND
C   IN CELL K
C   CFLAG:  AS IN FUNCTION FPATHA
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CZT1
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(6,3)
      INTEGER, INTENT(IN) :: K

      REAL(DP) :: DENIO(NPLS), ZTI(NPLS)
      REAL(DP) :: PVELQ(NPLS)
      REAL(DP) :: TBPI3(9), TBCX3(9), TBEL3(9)
      REAL(DP) :: EPPI3(9), EPCX3(9), EPEL3(9)

      REAL(DP) :: VREL, VRELQ, TBEL, FEPLCX3, CROSS, EXPO, ELB, CEL,
     .          SIGMAX, FEPLEL3, VX, VY, VZ, EHEAVY, FTABEI1, FPATHM,
     .          PVELQ0, TEEI, DENEL, FEELEI1, VEFF, VEFFQ, CXS, ELAB,
     .          TBCX, FEHVDS1, PLS
      INTEGER :: IBGK, IREL, IMEL, J, JAN, IF8, IML, IPL, IMDS, IRDS,
     .           II, IREAC, IMCX, IMPI, KK, IRCX
C
C
C  SET DEFAULTS: NO REACTIONS
C
      XSTOR=0.D0
      XSTORV=0.D0
      FPATHM=1.D10
C
      IF (LGVAC(K,0)) RETURN
C
C   LOCAL PLASMA PARAMETERS
C
      DENEL=DEIN(K)
      DO 2 IPLS=1,NPLSI
        ZTI(IPLS)=ZT1(IPLS,K)
2       DENIO(IPLS)=DIIN(IPLS,K)
      TEEI=TEIN(K)
      PVELQ0=VEL*VEL
      DO 3 IPLS=1,NPLSI
        IF (NLDRFT) THEN
          IF (INDPRO(4) == 8) THEN
            CALL VECUSR (2,VX,VY,VZ,IPLS)
          ELSE
            VX=VXIN(IPLS,K)
            VY=VYIN(IPLS,K)
            VZ=VZIN(IPLS,K)
          END IF
          PVELQ(IPLS)=(VELX*VEL-VX)**2+
     .                (VELY*VEL-VY)**2+
     .                (VELZ*VEL-VZ)**2
        ELSE
          PVELQ(IPLS)=PVELQ0
        ENDIF
3     CONTINUE
C
C
C  ELECTRON IMPACT RATE COEFFICIENT
C
      IF (LGMEI(IMOL,0).EQ.0.OR.LGVAC(K,NPLS+1)) GOTO 25
      DO 10 IMDS=1,NMDSI(IMOL)
        IRDS=LGMEI(IMOL,IMDS)
        IF (MODCOL(1,2,NSPA+IMOL,1).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            SIGVEI(IRDS)=TABDS1(IRDS,K)
          ELSE
            SIGVEI(IRDS)=FTABEI1(IRDS,K)
          END IF
        ELSE
          GOTO 990
        ENDIF
C
        IF (NSTORDR >= NRAD) THEN
          ESIGEI(IRDS,5)=EELDS1(IRDS,K)
          EHEAVY=EHVDS1(IRDS,K)
        ELSE
          ESIGEI(IRDS,5)=FEELEI1(IRDS,K)
          EHEAVY=FEHVDS1(IRDS,K)
        ENDIF
C
        ESIGEI(IRDS,1)=EATDS(IRDS,0,1)*E0+EATDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,2)=EMLDS(IRDS,0,1)*E0+EMLDS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,3)=EIODS(IRDS,0,1)*E0+EIODS(IRDS,0,2)*EHEAVY
        ESIGEI(IRDS,4)=EPLDS(IRDS,  1)*E0+EPLDS(IRDS,  2)*EHEAVY
C
        SIGEIT=SIGEIT+SIGVEI(IRDS)
10    CONTINUE
C
C
C  IONIZATION OF MOLECULE BY ION IMPACT, ION SPEZIES IPLS=1,NPLSI
C
25    CONTINUE
      DO 30 IMPI=1,NMPII(IMOL)
        SIGVPI(IMPI)=0.
        SIGPIT=SIGPIT+SIGVPI(IMPI)
30    CONTINUE
C
C   CHARGE EXCHANGE RATE COEFFICIENT FOR MOLECULE IMOL
C   WITH BULK IONS OF SPEZIES IPLS=1,NPLSI
C
      IF (LGMCX(IMOL,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 50
40    CONTINUE
      DO 41 IMCX=1,NMCXI(IMOL)
        IRCX=LGMCX(IMOL,IMCX,0)
        IPLS=LGMCX(IMOL,IMCX,1)
        IF (LGVAC(K,IPLS)) GOTO 41
C
        IF (MODCOL(3,2,NSPA+IMOL,IPLS).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVCX(IRCX)=TABCX3(IRCX,K,1)
          ELSE
            KK=NREACX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            TBCX = CREAC(9,1,KK)
            DO II=8,1,-1
              TBCX = TBCX*PLS + CREAC(II,1,KK)
            END DO
            TBCX=EXP(MAX(-100._DP,TBCX))*DIIN(IPLS,K)
            SIGVCX(IRCX)=TBCX
          END IF
        ELSEIF (MODCOL(3,2,NSPA+IMOL,IPLS).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLS,K).LT.TVAC) THEN
            VEFFQ=PVELQ(IPLS)
            VEFF=SQRT(VEFFQ)
            ELAB=LOG(VEFFQ)+DEFCX(IRCX)
            IREAC=MODCOL(3,1,NSPA+IMOL,IPLS)
            CXS=CROSS(ELAB,IREAC,IRCX,'FPATHM CX')
            SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBCX3(1:NSTORDT) = TABCX3(IRCX,K,1:NSTORDT)
              JAN=NSTORDT
            ELSE
              JAN=0
            END IF
            IF (JAN < 9) THEN
              KK=NREACX(IRCX)
              PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
              DO J=JAN+1,9
                TBCX3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                  TBCX3(J)=TBCX3(J)*PLS+CREAC(II,J,KK)
                END DO
              END DO
            END IF
            IF (JAN < 1) TBCX3(1)=TBCX3(1)+DIINL(IPLS,K)
C  MINIMUM ENERGY: 0.1EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
            EXPO=TBCX3(9)
            DO 43 II=1,8
              IF8=9-II
              EXPO=EXPO*ELB+TBCX3(IF8)
43          CONTINUE
            SIGVCX(IRCX)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(3,2,NSPA+IMOL,IPLS).EQ.3) THEN
C  BEAM  - BEAM
          VEFFQ=ZTI(IPLS)+PVELQ(IPLS)
          VEFF=SQRT(VEFFQ)
          ELAB=LOG(VEFFQ)+DEFCX(IRCX)
          IREAC=MODCOL(3,1,NSPA+IMOL,IPLS)
          CXS=CROSS(ELAB,IREAC,IRCX,'FPATHM CX')
          SIGVCX(IRCX)=CXS*VEFF*DENIO(IPLS)
        ELSE
          GOTO 992
        ENDIF
C
        SIGCXT=SIGCXT+SIGVCX(IRCX)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(3,4,NSPA+IMOL,IPLS).EQ.1) THEN
C  MODEL 1:
C  MEAN ENERGY FROM DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=2
        ELSEIF (MODCOL(3,4,NSPA+IMOL,IPLS).EQ.2) THEN
C  MODEL 2:
C  MEAN ENERGY FROM CROSS SECTION WEIGHTED DRIFTING MAXWELLIAN
          IF (NSTORDR >= NRAD) THEN
            EPCX3(1:NSTORDT) = EPLCX3(IRCX,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NELRCX(IRCX)
            PLS=TIINL(IPLS,K)+ADDCX(IRCX,IPLS)
            DO J=JAN+1,9
              EPCX3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                EPCX3(J)=EPCX3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) EPCX3(1)=EPCX3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFCX(IRCX))
          EXPO=EPCX3(9)
          DO 45 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+EPCX3(IF8)
45        CONTINUE
          ESIGCX(IRCX,1)=EXP(EXPO)/SIGVCX(IRCX)
          ESIGCX(IRCX,1)=ESIGCX(IRCX,1)+EDRIFT(IPLS,K)
          CFLAG(3,1)=3
        ELSEIF (MODCOL(3,4,NSPA+IMOL,IPLS).EQ.3) THEN
C  MODEL 3:
C  MEAN ENERGY FROM DRIFTING ONE SPEED DISTRIBUTION
          IF (NSTORDR >= NRAD) THEN
            ESIGCX(IRCX,1)=EPLCX3(IRCX,K,1)
          ELSE
            ESIGCX(IRCX,1)=FEPLCX3(IRCX,K)
          END IF
          CFLAG(3,1)=1
        ELSE
          GOTO 992
        ENDIF
41    CONTINUE
C
C  ELASTIC COLLISIONS OF MOLECULE IMOL WITH IONS OF SPEZIES IPLS=1,NPLS
C  50--->60
C
50    CONTINUE
      IF (LGMEL(IMOL,0,0).EQ.0.OR.LGVAC(K,0)) GOTO 60
      DO 51 IMEL=1,NMELI(IMOL)
        IREL=LGMEL(IMOL,IMEL,0)
        IPLS=LGMEL(IMOL,IMEL,1)
        IBGK=NPBGKP(IPLS,1)
        IF (LGVAC(K,IPLS)) GOTO 51
C
C  1.) RATE COEFFICIENT
C
        IF (MODCOL(5,2,NSPA+IMOL,IPLS).EQ.1) THEN
C  MAXWELL
          IF (NSTORDR >= NRAD) THEN
            SIGVEL(IREL)=TABEL3(IREL,K,1)
          ELSE
            KK=NREAEL(IREL)
            PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
            TBEL = CREAC(9,1,KK)
            DO II=8,1,-1
              TBEL = TBEL*PLS + CREAC(II,1,KK)
            END DO
            TBEL=EXP(MAX(-100._DP,TBEL))*DIIN(IPLS,K)
            SIGVEL(IREL)=TBEL
          END IF
        ELSEIF (MODCOL(5,2,NSPA+IMOL,IPLS).EQ.2) THEN
C  BEAM - MAXWELL
          IF (TIIN(IPLS,K).LT.TVAC) THEN
            VRELQ=PVELQ(IPLS)
            VREL=SQRT(VRELQ)
            ELAB=LOG(VRELQ)+DEFEL(IREL)
            IREAC=MODCOL(5,1,NSPA+IMOL,IPLS)
            CEL=CROSS(ELAB,IREAC,IREL,'FPATHM EL')
            SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
          ELSE
            IF (NSTORDR >= NRAD) THEN
              TBEL3(1:NSTORDT) = TABEL3(IREL,K,1:NSTORDT)
              JAN=NSTORDT
            ELSE
              JAN=0
            END IF
            IF (JAN < 9) THEN
              KK=NREAEL(IREL)
              PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
              DO J=JAN+1,9
                TBEL3(J)=CREAC(9,J,KK)
                DO II=8,1,-1
                  TBEL3(J)=TBEL3(J)*PLS+CREAC(II,J,KK)
                END DO
              END DO
            END IF
            IF (JAN < 1) TBEL3(1)=TBEL3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
            ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFEL(IREL))
            EXPO=TBEL3(9)
            DO 53 II=1,8
              IF8=9-II
              EXPO=EXPO*ELB+TBEL3(IF8)
53          CONTINUE
            SIGVEL(IREL)=EXP(EXPO)
          ENDIF
        ELSEIF (MODCOL(5,2,NSPA+IMOL,IPLS).EQ.3) THEN
C  BEAM - BEAM
          VRELQ=ZTI(IPLS)+PVELQ(IPLS)
          VREL=SQRT(VRELQ)
          ELAB=LOG(VRELQ)+DEFEL(IREL)
          IREAC=MODCOL(5,1,NSPA+IMOL,IPLS)
          CEL=CROSS(ELAB,IREAC,IREL,'FPATHM EL')
          SIGVEL(IREL)=CEL*VREL*DENIO(IPLS)
        ELSE
          GOTO 995
        ENDIF
        SIGELT=SIGELT+SIGVEL(IREL)
        IF (IBGK.NE.0) SIGBGK=SIGBGK+SIGVEL(IREL)
C
C  2.) BULK ION ENERGY LOSS RATE:
C
        IF (MODCOL(5,4,NSPA+IMOL,IPLS).EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            ESIGEL(IREL,1)=EPLEL3(IREL,K,1)
          ELSE
            ESIGEL(IREL,1)=FEPLEL3(IREL,K)
          END IF
          CFLAG(5,1)=2
        ELSEIF (MODCOL(5,4,NSPA+IMOL,IPLS).EQ.2) THEN
          IF (NSTORDR >= NRAD) THEN
            EPEL3(1:NSTORDT) = EPLEL3(IREL,K,1:NSTORDT)
            JAN=NSTORDT
          ELSE
            JAN=0
          END IF
          IF (JAN < 9) THEN
            KK=NELREL(IREL)
            PLS=TIINL(IPLS,K)+ADDEL(IREL,IPLS)
            DO J=JAN+1,9
              EPEL3(J)=CREAC(9,J,KK)
              DO II=8,1,-1
                EPEL3(J)=EPEL3(J)*PLS+CREAC(II,J,KK)
              END DO
            END DO
          END IF
          IF (JAN < 1) EPEL3(1)=EPEL3(1)+DIINL(IPLS,K)
C  MINIMUM PROJECTILE ENERGY: 0.1 EV
          ELB=MAX(-2.3_DP,LOG(PVELQ(IPLS))+EEFEL(IREL))
          EXPO=EPEL3(9)
          DO 55 II=1,8
            IF8=9-II
            EXPO=EXPO*ELB+EPEL3(IF8)
55        CONTINUE
          ESIGEL(IREL,1)=EXP(EXPO)/SIGVEL(IREL)
          ESIGEL(IREL,1)=ESIGEL(IREL,1)+EDRIFT(IPLS,K)
          CFLAG(5,1)=3
        ELSE
          GOTO 995
        ENDIF
51    CONTINUE
C
60    CONTINUE
C
C     TOTAL
C
      SIGMAX=MAXVAL(XSTOR(:,1:4))
      WHERE (XSTOR(:,1:4) .LE. SIGMAX*1.D-10 )
        XSTOR(:,1:4) = 0.D0
      END WHERE
C
100   CONTINUE
      SIGTOT=SIGPIT+SIGCXT+SIGEIT+SIGELT
      IF (SIGTOT.GT.1.D-20) THEN
        FPATHM=VEL/SIGTOT
        ZMFPI=1./FPATHM
      ENDIF
C
      RETURN
990   CONTINUE
      IML=NSPA+IMOL
      WRITE (6,*) 'ERROR IN FPATHM: INCONSISTENT ELEC. IMP. DATA'
      WRITE (6,*) 'IML,MODCOL(1,J,IML,1) '
      WRITE (6,*) IMOL,(MODCOL(1,J,IML,1),J=1,4)
      CALL EXIT
991   CONTINUE
      IML=NSPA+IMOL
      WRITE (6,*) 'ERROR IN FPATHM: INCONSISTENT ION IMP. IONIZ. DATA'
      WRITE (6,*) 'IML,IPL,MODCOL(4,J,IML,IPL) '
      DO 993 IPL=1,NPLSI
        WRITE (6,*) IMOL,IPL,(MODCOL(4,J,IML,IPL),J=1,4)
993   CONTINUE
      CALL EXIT
992   CONTINUE
      IML=NSPA+IMOL
      WRITE (6,*) 'ERROR IN FPATHM: INCONSISTENT CHARGE EXCHANGE DATA'
      WRITE (6,*) 'IML,IPL,(MODCOL(3,J,IML,IPL),J=1,4) '
      DO 994 IPL=1,NPLSI
        WRITE (6,*) IMOL,IPL,(MODCOL(3,J,IML,IPL),J=1,4)
994   CONTINUE
      CALL EXIT
995   CONTINUE
      IML=NSPA+IMOL
      WRITE (6,*) 'ERROR IN FPATHM: INCONSISTENT ELASTIC COLL. DATA'
      WRITE (6,*) 'IML,IPL,(MODCOL(5,J,IML,IPL),J=1,4) '
      DO 996 IPL=1,NPLSI
        WRITE (6,*) IMOL,IPL,(MODCOL(5,J,IML,IPL),J=1,4)
996   CONTINUE
      CALL EXIT
      END
C
C
      FUNCTION FPATHPH (K,CFLAG)
C
C
      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: CFLAG(6,3)
      INTEGER, INTENT(IN) :: K
      REAL(DP) :: FPATHPH

      FPATHPH = HUGE(1._DP)
      
      RETURN
      END FUNCTION FPATHPH


      FUNCTION FTABEI1 (IREI,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IREI, K
      REAL(DP) :: TBEIC(9), FTABEI1, DEIMIN, DSUB, PLS, TBEI
      INTEGER :: J, I, II, KK

      TBEI=0.D0
      KK = NREAEI(IREI)
      IF (JEREAEI(IREI) == 1) THEN
        TBEI = CREAC(9,1,KK)
        DO II=8,1,-1
          TBEI = TBEI*TEINL(K) + CREAC(II,1,KK)
        END DO
        TBEI = TBEI + FACREA(KK)
        TBEI=EXP(MAX(-100._DP,TBEI))*DEIN(K)
      ELSE
        TBEIC(1:JEREAEI(IREI)) = CREAC(9,1:JEREAEI(IREI),KK)
        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
        PLS=MAX(DEIMIN,DEINL(K))-DSUB
        DO J=1,JEREAEI(IREI)
          DO II=8,1,-1
            TBEIC(J)=TBEIC(J)*TEINL(K)+CREAC(II,J,KK)
          END DO
        END DO
        TBEI = TBEIC(9)
        DO I=8,1,-1
          TBEI=TBEI*PLS+TBEIC(I)
        END DO
        TBEI=MAX(-100._DP,TBEI+FACREA(KK)+DEINL(K))
        TBEI=EXP(TBEI)
      ENDIF

      FTABEI1 = TBEI

      RETURN
      END


      FUNCTION FTABRC1 (IRRC,K)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IRRC, K
      REAL(DP) :: TBRCC(9), PLS(1), COUN(0:9,1), CF(9,0:9)
      REAL(DP) :: DEIMIN, DSUB, FTABRC1, ZX, TBRC
      INTEGER :: II, KK

      TBRC=0.D0
      KK = NREARC(IRRC)

      IF (KK == 0) THEN
        ZX=EIONH/MAX(1.E-5_DP,TEIN(K))
        TBRC=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(K)
      ELSEIF (JEREARC(IRRC) == 1) THEN
        TBRC = CREAC(9,1,KK)
        DO II=8,1,-1
          TBRC = TBRC*TEINL(K) + CREAC(II,1,KK)
        END DO
        TBRC=TBRC+FACREA(KK)
        TBRC=EXP(MAX(-100._DP,TBRC))*DEIN(K)
      ELSE
        DSUB=LOG(1.D8)
        DEIMIN=LOG(1.D8)
        PLS(1)=MAX(DEIMIN,DEINL(K))-DSUB
        CALL CDEFN(TEINL(K),PLS,KK,COUN,1,CF,.TRUE.,.FALSE.,
     ,             .TRUE.)
        TBRC=COUN(1,1)*DEIN(K)*FACREA(KK)
      ENDIF

      FTABRC1 = TBRC

      RETURN
      END
C
C
      SUBROUTINE GAUMEH(RS,ER,B,IFLAG,P,N,IGAUS,RESULT)
C
C  IGAUS=1
C  GAUSS MEHLER QUADRATURE, N=5,10,20, W(X)=1./SQRT(1+X)/SQRT(1-X)
C                           INTEGRATION FROM A=-1 TO B=1
C                           TRANSFORMED TO A=0,B=1, AND
C                           W(X)=1./SQRT(X)/SQRT(1-X)
C

      USE PRECISION

      IMPLICIT NONE 
      REAL(DP), INTENT(IN) :: RS, ER, B, P(*)
      REAL(DP), INTENT(OUT) :: RESULT
      INTEGER, INTENT(IN) :: IFLAG, N, IGAUS

      REAL(DP) :: AR(128), AFI(128)
      INTEGER :: NFI
      COMMON /CFI/ AR,AFI,NFI
      REAL(DP) :: X5A(5),X5B(5),X10A(10),X10B(10),X20A(20),X20B(20)
      REAL(DP) :: X5(5),X10(10),X20(20)
C
      REAL(DP) :: XG10(10),WG10(10),XG10A(10),XG10B(10)
      REAL(DP) :: SUM, F, X, W5, W10, W20
      INTEGER :: I, IFIRST, IFI
C
      SAVE
      DATA X5 /0.95105654E+00,0.58778542E+00,0.31391647E-06,
     .       - 0.58778363E+00,-0.95105600E+00/,
     .W5 /0.62831837E+00/
      DATA X10/9.876884E-01, 8.910065E-01,7.071068E-01,4.539905E-01,
     .         1.564344E-01,-1.564345E-01,-4.539905E-01,-7.071068E-01,
     .        -8.910065E-01,-9.876884E-01/,
     .W10/3.141593E-01/
      DATA X20/9.969173E-01,9.723699E-01,9.238795E-01,8.526402E-01,
     .         7.604060E-01,6.494480E-01,5.224985E-01,3.826834E-01,
     .         2.334453E-01,
     .         7.845905E-02,-7.845914E-02,-2.334454E-01,-3.826835E-01,
     .        -5.224986E-01,-6.494481E-01,-7.604060E-01,-8.526402E-01,
     .        -9.238796E-01,-9.723700E-01,-9.969174E-01/,
     .W20/1.570796E-01/
C
C  IGAUS=2
C  GAUSS MEHLER QUADRATURE, N=10, W(X)=1./SQRT(1-X)
C                           INTEGRATION FROM A=0 TO B=1
      DATA XG10/
     .    .013695585480651072296,
     .    .070758123420104485088,
     .    .16782834791297652828,
     .    .29588270759990962383,
     .    .44298868539955669668,
     .    .59543571523425259033,
     .    .73901490631777356542,
     .    .86034375925702257262,
     .    .94811360601967719799,
     .    .99414369156320387287/
      DATA WG10/
     .    .035228014278304236584,
     .    .081202859600773882686,
     .    .12534409666821812627,
     .    .16655348315340949082,
     .    .20386023963448093937,
     .    .23638906392303662412,
     .    .26337727689835336360,
     .    .28419221863676482994,
     .    .29834597294520589473,
     .    .30550677426145261144/
      DATA IFIRST /0/
      IF (IFIRST.EQ.0) THEN
C SET ROOTS AND WEIGHTS FOR GAUSS QUADRATURE RULES
        IFIRST=1
        DO 5 I=1,5
          X=0.5*(1.+X5(I))
          X5A(I)=SQRT(X*(1.-X))
          X5B(I)=1./X
5       CONTINUE
        DO 10 I=1,10
          X=0.5*(1.+X10(I))
          X10A(I)=SQRT(X*(1.-X))
          X10B(I)=1./X
10      CONTINUE
        DO 11 I=1,10
          X=XG10(I)
          XG10A(I)=SQRT(1.-X)
          XG10B(I)=1./X
11      CONTINUE
        DO 20 I=1,20
          X=0.5*(1.+X20(I))
          X20A(I)=SQRT(X*(1.-X))
          X20B(I)=1./X
20      CONTINUE
      ENDIF
C
      NFI=N
      IF (N.EQ.5) THEN
        DO 50 IFI=1,5
          AR(IFI)=RS*X5B(IFI)
50      CONTINUE
C
        CALL FIVEC(ER,B,IFLAG,P)
        SUM=0.D0
        DO 51 IFI=1,5
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X5A(IFI)/SQRT(F)
51      CONTINUE
        RESULT=W5*SUM
C
      ELSEIF (N.EQ.10) THEN
        IF (IGAUS.EQ.1) THEN
          DO 100 IFI=1,10
            AR(IFI)=RS*X10B(IFI)
100       CONTINUE
C
          CALL FIVEC(ER,B,IFLAG,P)
          SUM=0.D0
          DO 101 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+X10A(IFI)/SQRT(F)
101       CONTINUE
          RESULT=W10*SUM
        ELSEIF (IGAUS.EQ.2) THEN
          DO 110 IFI=1,10
            AR(IFI)=RS*XG10B(IFI)
110       CONTINUE
C
          CALL FIVEC(ER,B,IFLAG,P)
          SUM=0.D0
          DO 111 IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+XG10A(IFI)/SQRT(F)*WG10(IFI)
111       CONTINUE
          RESULT=SUM
        ENDIF
C
      ELSEIF (N.EQ.20) THEN
        DO 200 IFI=1,20
          AR(IFI)=RS*X20B(IFI)
200     CONTINUE
C
        CALL FIVEC(ER,B,IFLAG,P)
        SUM=0.D0
        DO 201 IFI=1,20
          F=AFI(IFI)
          IF (F.GT.0.D0)
     .    SUM=SUM+X20A(IFI)/SQRT(F)
201     CONTINUE
        RESULT=W20*SUM
C
      ELSE
        WRITE (6,*) 'ERROR IN GAUMEH, WRONG PARAMETER N'
        CALL EXIT
      ENDIF
      RETURN
      END
C
C
      FUNCTION RSTERN (ER,B,IFLAG,P)
C
C     IFLAG=1:  H+  + H, PURELY REPULSIVE POTENTIAL
C     IFLAG=2:  MORSE POTENTIAL, He+ + He, H+ + Noble Gases, H+ + H2
C
C               PARAMETERS: P(1):
C                           P(2):
C                           P(3):
C
C  RSTERN(ER,B) IS THE LARGEST ROOT OF THE EQUATION:
C
C     FI(R):=1.-V(R)/ER-(B/R)**2=0.
C
C   HERE:  ER COLLISION ENERGY (EV)
C          B  IMPACT PARAMETER
C          V  INTERACTION POTENTIAL (EV)
C
C     DIMENSION RV0(3),RVM(3),RVW(3),VM(3),VW(3),VSW(3)

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG

      REAL(DP) :: FIW, RV, RTSAF, BQ, TOL, RSTERN, FI, DFI, FITEST,
     .          rup, rlw, vw, rw, r0
      SAVE
C     DATA IFIRST /0/
C     DATA RV0/0.,0.99699,2.18039/
C     DATA RVM/0.,1.4556 ,2.835539/
C     DATA RVW/0.,1.99515,3.490687/
C     DATA VM /0.,-2.     ,-2.7/
C     DATA VW /0.,-1.5    ,-2.025/
C     DATA VSW/0.,1.284688,1.4283/
C
C
      IF (IFLAG.EQ.1) THEN
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAF
C  FOR IFLAG=1, RLW=B IS A LOWER BOUND AND RUP=INF IS AN UPPER BOUND.
C  FIRST TRY TO FIND BETTER BOUNDS:
        RUP=B
C
10      RLW=RUP
        RUP=RLW*2.
        FITEST=FI(RUP,ER,B,IFLAG,P,DFI)
        IF (FITEST.LT.0.D0) GOTO 10
C
        TOL=1.D-7*(RUP-RLW)
        RSTERN=RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
C
      ELSEIF (IFLAG.EQ.2) THEN
C  GENERALISED MORSE POTENTIAL
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAFE
C  FOR IFLAG=2 AND IFLAG=3, THE FOLLOWING INTERVAL CAN BE USED:
C  R0 IS THE ROOT OF THE INTERACTION POTENTIAL V(R): V(R0)=0.
C  RM IS THE RADIUS OF THE MINIMUM OF V(R)           V(RM)= VM = -EPS
C  RW IS THE RADIUS OF THE POINT OF INFLECTION       V(RM)= VW 
C        RM=P(4)
         R0=P(5)
         RW=P(6)
         VW=P(8)
C
        BQ=B*B
        RV=R0
C
C  CASE 1: B < R0,
C
        IF (RV.GT.B) THEN
          RLW=B
          RUP=RV
        ENDIF
C
C  CASE 2: R0 < B,
C
        IF (RV.LT.B) THEN
          RLW=RV
          RUP=B
C  TRY TO IMPROVE LOWER BOUND RLW
          FIW=1.-VW/ER-BQ/(RW**2)
          IF (FIW.LT.0.D0) RLW=RW
        ENDIF
C
        TOL=1.D-7*(RUP-RLW)
        RSTERN=RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
      ELSE
        GOTO 990
      ENDIF
C
      RETURN
990   CONTINUE
      WRITE (6,*) 'ERROR IN FUNCTION RSTERN '
      CALL EXIT
      END
C
C
      FUNCTION RTSAF(X1,X2,XACC,ER,B,IFLAG,P)

C  TAKEN FROM: NUMERICAL RECIPES, W.H.PRESS ET AL.,
C              CAMBRIDGE UNIV. PRESS, 1989, P258
C  MODIFIED, TO FIND LARGEST ROOT, IN CASE MORE THAN ONE ROOTS
C  AND TO SPEED UP

      USE PRECISION

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: X1, X2, XACC, ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      INTEGER, PARAMETER :: MAXIT=100
      REAL(DP) :: DF, XDIST, F, FI, TEMP, DX, DXOLD, XH, RTSAF, XL
      INTEGER :: J

      XL=X1
      XH=X2
C
      F=FI(XH,ER,B,IFLAG,P,DF)
      XDIST=(XH-XL)
      DX=MIN(XDIST,F/DF)
      DXOLD=DX
      RTSAF=XH-DX
C
      F=FI(RTSAF,ER,B,IFLAG,P,DF)

      IF(F.LT.0.D0) THEN
        XL=RTSAF
      ELSE
        XH=RTSAF
      ENDIF
C
      DO 11 J=1,MAXIT
        IF(((RTSAF-XH)*DF-F)*((RTSAF-XL)*DF-F).GT.0.
     *      .OR. ABS(F+F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5*(XH-XL)
          RTSAF=XL+DX
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=RTSAF
          RTSAF=RTSAF-DX
        ENDIF
C
        IF(ABS(DX).LT.XACC) RETURN
C
        F=FI(RTSAF,ER,B,IFLAG,P,DF)
        IF(F.LT.0.D0) THEN
          XL=RTSAF
        ELSE
          XH=RTSAF
        ENDIF
11    CONTINUE
      WRITE (6,*) 'RTSAF EXCEEDING MAXIMUM ITERATIONS'
      RETURN
      END
C
C
      SUBROUTINE SCATANG (INENERGY,INRAN,ELTHETA,CTTHETA,SIGM)
c     ***********************************************************
c     *                       --- TH I ---                      *
c     *                                                         *
c     *Programm zur Berechnung der Streuwinkel bei vorgegebener *
c     *Energie und Zufallszahl (Monte-Carlo-Simulation) für eine*
c     *atomare Reaktion.                                        *
c     *Die verschiedenen Reaktionen werden mittels der Daten von*
c     *P.S. Krstic und D.R. Schultz (Atomic and Plasma-Material *
c     *Interaction Data for Fusion), Differntielle Querschnitte *
c     *mit zugehörigen Winkeln, berechnet.                      *
c     *Das Programm berechnet Reaktionen zwischen Energien von  *
c     *0.1 eV und 100 eV (CM). Dabei sind 31 Energiewerte aus   *
c     *den o.g.Daten vorgegeben.                                *
c     *Input:    Dateien der jeweiligen Reaktion                *
c     *          Reaktionsenergy: inenergy                      *
c     *          Zufallszahl: inran                             *
c     *Output:   Verwendete Energie: energy                     *
c     *          Streuwinkel: eltheta (elastisch)               *
c     *                       cttheta (ladungsaustausch)        *
c     *        (weiteres s.u.)                                *
c     *                                                         *
c     *(19.9.2000)                          Torsten Haberscheidt*
c     ***********************************************************

c     --Parameter--
      

      USE PRECISION
      IMPLICIT NONE
      integer i,j,k,l,m,n,dat,IFIRST,IUN
      REAL(DP) :: pi,z,w
!      parameter (pi=3.1415926535D+000)
      REAL(DP) :: energy(31),inenergy,hasard
      REAL(DP) :: theta(768,31),el(770),ct(770),dtheta(770)
      REAL(DP) :: sig,sigma(31),gam,gamma(31),SIGM
      REAL(DP) :: mom,moment(31),vis,viscos(31)
      REAL(DP) :: normel(770), normct(770)
      REAL(DP) :: elangle(768,31), ctangle(768,31)
      REAL(DP) :: inran, eltheta, cttheta
      CHARACTER(10) :: FILNAM
      DATA IFIRST /0/
      SAVE
      
c     ***********************************************************
c     -- Einlesen der Daten --

      IF (IFIRST == 0) THEN
        IFIRST = 1

        PI=4.D0*ATAN(1.D0)
        iun=23

        do 10 i=1,31
          IF (I <= 10) THEN
            WRITE (FILNAM,'(a3,i1,a4)') 'el-',i-1,'.dat'
          ELSE
            WRITE (FILNAM,'(a3,i2,a4)') 'el-',i-1,'.dat'
          END IF
        
          open (iun,file=filnam)

          read(iUN,*)
          read(iUN,100) energy(i)
100       format(T19,E9.4)

          do 20 j=1,13
            read(iUN,*)
20        continue
      
          dat=0

          do 30 k=1,768
            read(iun,200,end=5) theta(k,i), el(k), ct(k)
            dat= dat+1
30        continue

200       format(T5, E11.6, T20, E11.6, T35, E11.6)
5         continue
          close(iUN)

c       *******************************************************
c       --Verarbeitung der Daten--

c       Winkelintervalle
      
          dtheta(1) = 0.5*(theta(1,i) + theta(2,i))

          do 40 k=2,dat-1
            dtheta(k) = 0.5*(theta(k+1,i) - theta(k-1,i))
40        continue

          dtheta(dat) = pi - 0.5*(theta(dat-1,i) + theta(dat,i))

c       Streuquerschnitte

          sig=0
          gam=0
          mom=0
          vis=0

          do 50 k=1,dat
            sig = sig + (dtheta(k) * el(k))
            gam = gam + (dtheta(k) * ct(k))
            mom = mom + (dtheta(k) * el(k))*(1-cos(theta(k,i)))
            vis = vis + (dtheta(k) * el(k))*(sin(theta(k,i)))**2
50        continue

          sigma(i) = sig
          gamma(i) = gam
          moment(i)= mom
          viscos(i)= vis

c       Normierung

          do 60 k=1,dat
            normel(k) = (dtheta(k) * el(k)) / sig
            normct(k) = (dtheta(k) * ct(k)) / gam
60        continue

c       Generieren von Zahlen [0;1]

          elangle(1,i) = normel(1)
          ctangle(1,i) = normct(1)

          do 70 k=2,dat
            elangle(k,i) = elangle(k-1,i) + normel(k)
            ctangle(k,i) = ctangle(k-1,i) + normct(k)
70        continue



10      continue

      END IF   ! END OF IFIRST-BLOCK

c     ***********************************************************

c     Berechnung des Monte-Carlo-Winkels mit linearer Interpolation

      call energyloc(energy,31,inenergy,m)

c     elastic
      if(inran.lt.elangle(1,m))then
        eltheta = inran*theta(1,m)/elangle(1,m)
      else if(inran.eq.elangle(dat,m))then
        eltheta = pi
      else
        call thetaloc(elangle(1:dat,m),dat,inran,n)
        z=(inran - elangle(n,m))/(elangle(n+1,m) - elangle(n,m))+n
      
        if(z.gt.(dat-1))then
          eltheta = (z-n)*(pi - theta(n,m)) + theta(n,m)
        else
          eltheta = (z-n)*(theta(n+1,m) - theta(n,m))+theta(n,m)
        endif
      endif

c     chargetransfer
      if(inran.lt.ctangle(1,m))then
        cttheta = inran*theta(1,m)/ctangle(1,m)
      else if(inran.eq.ctangle(dat,m))then
        cttheta = pi
      else
        call thetaloc(ctangle(1:dat,m),dat,inran,l)
        w=(inran - ctangle(l,m))/(ctangle(l+1,m) - ctangle(l,m))+l
      
        if(w.gt.(dat-1))then
          cttheta = (w-l)*(pi - theta(l,m)) + theta(l,m)
        else
          cttheta = (w-l)*(theta(l+1,m) - theta(l,m))+theta(l,m)
        endif
      endif

      SIGM=SIGMA(M)

c     ***********************************************************
c     -- Output des Programms--

c     Verwendete Energie in eV: energy(m)
c     Berechnete Streuwinkel: eltheta, cttheta

c     Totale Streuquerschnitte + Momente für alle Energien:
c       sigma(i), gamma(i), moment(i), viscos(i)
c     mit dem Feldindex i der Energie (zugeh. Vorgabe i=m).

c     Mögliche Energien: Feld energy(i)
c     Winkel der diff.Querschnitte abh.von Energie: theta(k,i)
c     mit k=1..768

c     ***********************************************************

      RETURN
      end
C
C
      SUBROUTINE SETAMD(ICAL)
C
C  SET ATOMIC AND MOLECULAR DATA: DRIVER
C
      USE PRECISION
      USE PARMMOD
      USE COMXS
      USE COMSOU
      USE CZT1

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ICAL
      INTEGER :: I, IRPI, IRDS

      IF (ICAL == 0) THEN
        NRCX=0
        NREL=0
        NRPI=0
        NRDS=0
        NREC=0
        NBGV=0
        CALL XSECTA_PARAM
        CALL XSECTM_PARAM
        CALL XSECTI_PARAM
        CALL XSECTP_PARAM
        CALL XSECTPH_PARAM
        
        NRCX=MAX(1,NRCX)
        NREL=MAX(1,NREL)
        NRPI=MAX(1,NRPI)
        NRDS=MAX(1,NRDS)
        NREC=MAX(1,NREC)
        NBGV=MAX(1,NBGV)
        
        CALL SET_PARMMOD(2)
        CALL ALLOC_COMXS(2)
        CALL ALLOC_COMSOU(2)
        CALL ALLOC_CZT1(2)
        RETURN
      ELSE
        CALL INIT_CMDTA(2)
      END IF

      NRCXI=0
      NRELI=0
      NRPII=0
      NREII=0
      NRRCI=0
      NRBGI=0
      CALL XSECTA
      CALL XSECTM
      CALL XSECTI
      CALL XSECTP
      CALL XSECTPH
      CALL CONDENSE

      IPATDS = 0
      IPMLDS = 0
      IPIODS = 0 
      IPPLDS = 0
      DO IRDS=1,NRDS
        ipatds(IRDS,0)=COUNT(PATDS(IRDS,1:) > 0)
        IPATDS(IRDS,1:ipatds(IRDS,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATDS(IRDS,1:) > 0)
        ipmlds(IRDS,0)=COUNT(PMLDS(IRDS,1:) > 0)
        IPMLDS(IRDS,1:ipmlds(IRDS,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLDS(IRDS,1:) > 0)
        ipiods(IRDS,0)=COUNT(PIODS(IRDS,1:) > 0)
        IPIODS(IRDS,1:ipiods(IRDS,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIODS(IRDS,1:) > 0)
        ipplds(IRDS,0)=COUNT(PPLDS(IRDS,1:) > 0)
        IPPLDS(IRDS,1:ipplds(IRDS,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLDS(IRDS,1:) > 0)
      END DO

      IPATPI = 0
      IPMLPI = 0
      IPIOPI = 0
      IPPLPI = 0
      DO IRPI=1,NRPI
        ipatpi(IRPI,0)=COUNT(PATPI(IRPI,1:) > 0)
        IPATPI(IRPI,1:ipatpi(IRPI,0))=PACK( (/ (i,i=1,natm) /),
     .                                     PATPI(IRPI,1:) > 0)
        ipmlpi(IRPI,0)=COUNT(PMLPI(IRPI,1:) > 0)
        IPMLPI(IRPI,1:ipmlpi(IRPI,0))=PACK( (/ (i,i=1,nmol) /),
     .                                     PMLPI(IRPI,1:) > 0)
        ipiopi(IRPI,0)=COUNT(PIOPI(IRPI,1:) > 0)
        IPIOPI(IRPI,1:ipiopi(IRPI,0))=PACK( (/ (i,i=1,nion) /),
     .                                     PIOPI(IRPI,1:) > 0)
        ipplpi(IRPI,0)=COUNT(PPLPI(IRPI,1:) > 0)
        IPPLPI(IRPI,1:ipplpi(IRPI,0))=PACK( (/ (i,i=1,npls) /),
     .                                     PPLPI(IRPI,1:) > 0)
      END DO

      RETURN
      END

      

      subroutine thetaloc(xx,a,x,b)
c     ***********************************************************
c     * Ermittlung des Feldindex b mit vorgegebener RAN-Zahl x, *
c     * so dass x zwischen RAN(b) und RAN(b+1)                  *
c     * Modifikation: für x=RAN(bmax) erhält man bmax           *
c     * *********************************************************
      USE PRECISION
      implicit none
      integer a,b
      REAL(DP) :: x,xx(a)
      integer bl,bm,bu

      bl=0
      bu=a+1
      
90    if(bu-bl.gt.1)then
        bm=(bu+bl)*0.5
        if((xx(a).ge.xx(1)).eqv.(x.ge.xx(bm)))then
          bl=bm
        else
          bu=bm
        endif
      goto 90
      endif

      if(x.eq.xx(1))then
        b=1
      else if(x.eq.xx(a))then
        b=a
      else
        b=bl
      endif

      return
      end
C
C
      SUBROUTINE XSECTA
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMSOU
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: FACTKK, CHRDF0, EELEC, RMASS, DSUB, DEIMIN, EHEAVY,
     .          EBULK, TMASS, PMASS
      INTEGER :: II, IML, IM, IIO, NTE, ISTORE, ISCND, ISCDE, IFRST,
     .           IAT, IREI, IATM, IDSC1, J, IPLS1, IPLS, IION1, NRC,
     .           KK, ISPZB, IAEL, ITYPB, IREL, IBGK, IA, ISP, IP,
     .           IAPI, IRPI, IACX, IDSC, IPL, IAEI, IESTM, IDEZ, IRCX
      CHARACTER(8) :: TEXTS1, TEXTS2

      SAVE
C
C
C   ELECTRON IMPACT COLLISIONS:
C
C  FIND SPECIES INDEX OF ION AFTER IONIZATION EVENT FOR THE DEFAULT
C  ELECTRON IMPACT IONIZATION MODELS FROM INPUT MASS AND
C  AND CHARGE NUMBER
C
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
70      CONTINUE
      END IF
C
      DO 100 IATM=1,NATMI
        IDSC1=0
        LGAEI(IATM,0)=0
C
        DO NRC=1,NRCA(IATM)
          KK=IREACA(IATM,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
        IF (NRCA(IATM).EQ.0.AND.NCHARA(IATM).LE.2) THEN
C
C  DEFAULT H,D,T OR HE ELEC. IMP. IONIZATION MODEL
C
          IION1=0
          IPLS1=0
          DO 52 IPLS=1,NPLSI
            IF (NCHARP(IPLS).EQ.NCHARA(IATM).AND.
     .          NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .          NCHRGP(IPLS).EQ.1) THEN
              IPLS1=IPLS
C
              IDSC1=IDSC1+1
              NREII=NREII+1
              IF (NREII.GT.NRDS) GOTO 993
              IREI=NREII
              LGAEI(IATM,IDSC1)=IREI
C
              PELDS(IREI)=1.
              PPLDS(IREI,IPLS1)=1.
              EPLDS(IREI,1)=1.D0
              EPLDS(IREI,2)=0.D0
              GOTO 50
            ENDIF
52        CONTINUE
          GOTO 100
C
50        CONTINUE
          NTE=NSBOX
          IF (NSTORDR < NRAD) NTE=1
          IF (NCHARA(IATM).EQ.1) THEN
            ISTORE=-4
            EELEC=-EIONH
            CALL CDEF (TEINL,1,1,ISTORE,COUN,NTE,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
          ELSEIF (NCHARA(IATM).EQ.2) THEN
            ISTORE=-1
            EELEC=-EIONHE
            CALL CDEF (TEINL,1,1,ISTORE,COUN,NTE,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
          ENDIF
C
          IF (NSTORDR >= NRAD) THEN
            DO 80 J=1,NSBOX
              TABDS1(IREI,J)=COUN(1,J)*DEIN(J)
80          CONTINUE
C  NO RADIATION LOSS INCLUDED
            EELDS1(IREI,1:NSBOX)=EELEC
          ELSE
            CREAC(1:9,1,ISTORE)=CF(1:9,1)
            FACREA(ISTORE) = 0._DP
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
          ENDIF
          MODCOL(1,2,NSPH+IATM,1)=1
          MODCOL(1,4,NSPH+IATM,1)=1
C
C  TRACKLENGTH ESTIMATOR FOR ALL COLLISION RATE CONTRIBUTIONS
C
          IESTEI(IREI,1)=0
          IESTEI(IREI,2)=0
          IESTEI(IREI,3)=0
C
          NAEII(IATM)=IDSC1
C
C  NON DEFAULT ELEC. IMP. COLLISION MODEL,
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 90 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
            IAT=NSPH+IATM
            RMASS=RMASSA(IATM)
            IFRST=ISCD1A(IATM,NRC)
            ISCND=ISCD2A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EHEAVY=ESCD1A(IATM,NRC)+ESCD2A(IATM,NRC)
            EELEC=EELECA(IATM,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IF (NREII.GT.NRDS) GOTO 993
            IREI=NREII
            LGAEI(IATM,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IAT,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
90        CONTINUE
          NAEII(IATM)=IDSC1
       ENDIF
C
        NAEIIM(IATM)=NAEII(IATM)-1
        LGAEI(IATM,0)=NAEII(IATM)
C
        DO 95 IAEI=1,NAEII(IATM)
          IREI=LGAEI(IATM,IAEI)
          CALL XSTEI_1(IREI)
95      CONTINUE
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:
C
      DO 200 IATM=1,NATMI
        IDSC=0
        LGACX(IATM,0,0)=0
        LGACX(IATM,0,1)=0
C
C   HYDROGENIC AND HELIUM DEFAULT MODEL 100 --- 140
C
        IF (NRCA(IATM).EQ.0) THEN
          DO 122 IPLS=1,NPLSI
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND IPLS
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 122
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS).AND.NCHRGP(IPLS).EQ.1)
     .          GOTO 125
124           CONTINUE
              GOTO 122
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
              IDSC=IDSC+1
              NRCXI=NRCXI+1
              IRCX=NRCXI
              LGACX(IATM,IDSC,0)=IRCX
              LGACX(IATM,IDSC,1)=IPLS
              N1STX(IRCX,1)=1
              N1STX(IRCX,2)=IAT
              N1STX(IRCX,3)=1
              N2NDX(IRCX,1)=4
              N2NDX(IRCX,2)=IPL
              N2NDX(IRCX,3)=1
C  PROJECTILE MASS IS 1.
C  TARGET     MASS IS 1.
              PMASS=1.*PMASSA
              TMASS=1.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-1
              MODCOL(3,1,NSPH+IATM,IPLS)=-1
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
              DEFCX(IRCX)=LOG(CVELI2*PMASS)
              EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  TRACKLENGTH ESTIMATOR FOR ALL COLLISION RATE CONTRIBUTIONS
C
              IESTCX(IRCX,1)=0
              IESTCX(IRCX,2)=0
              IESTCX(IRCX,3)=0
C
C  DEFAULT BULK ION ENERGY LOSS RATE = 1.5*TI+EDRIFT PER COLLISION
C
              IF (NSTORDR >= NRAD) THEN
                DO 127 J=1,NSBOX
                  EPLCX3(IRCX,J,1)=1.5*TIIN(IPLS,J)+EDRIFT(IPLS,J)
127             CONTINUE
              ELSE
                NELRCX(IRCX) = -1
              END IF
C
              MODCOL(3,2,NSPH+IATM,IPLS)=3
              MODCOL(3,4,NSPH+IATM,IPLS)=3
C
            ENDIF
122       CONTINUE
C
          NACXI(IATM)=IDSC
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IDSC=IDSC+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGACX(IATM,IDSC,0)=IRCX
            LGACX(IATM,IDSC,1)=IPLS
            FDLMCX(IRCX)=FLDLMA(IATM,NRC)
            IAT=NSPH+IATM
            IPL=IPLS
            RMASS=RMASSA(IATM)
            IFRST=ISCD1A(IATM,NRC)
            ISCND=ISCD2A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EBULK=EBULKA(IATM,NRC)
            CALL XSTCX(RMASS,IRCX,IAT,IPL,
     .                 IFRST,ISCND,EBULK,ISCDE,IESTM,KK,FACTKK)
C
130       CONTINUE
C
          NACXI(IATM)=IDSC
C  NO CX MODEL DEFINED
        ELSE
          NACXI(IATM)=0
        ENDIF
C
        NACXIM(IATM)=NACXI(IATM)-1
C
        LGACX(IATM,0,0)=0.
        DO 180 IACX=1,NACXI(IATM)
          LGACX(IATM,0,0)=LGACX(IATM,0,0)+LGACX(IATM,IACX,0)
180     CONTINUE
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IATM=1,NATMI
        IDSC=0
        LGAEL(IATM,0,0)=0
        LGAEL(IATM,0,1)=0
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCA(IATM).EQ.0) THEN
          NAELI(IATM)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 230 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGAEL(IATM,IDSC,0)=IREL
            LGAEL(IATM,IDSC,1)=IPLS
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TESTPARTICLES
            IF (IBGKA(IATM,NRC).NE.0) THEN
              IF (NPBGKA(IATM).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKA(IATM)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKA(IATM)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKA(IATM,NRC),1,3)
              ISPZB=IDEZ(IBGKA(IATM,NRC),3,3)
              IF (ITYPB.NE.1.OR.ISPZB.NE.IATM) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKA(IATM,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C  BGK-COLLISION PARAMETERS DONE
C
            IAT=NSPH+IATM
            IPL=IPLS
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EBULK=EBULKA(IATM,NRC)
            CALL XSTEL(IREL,IAT,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE

          NAELI(IATM)=IDSC
        ENDIF
C
        NAELIM(IATM)=NAELI(IATM)-1
C
        LGAEL(IATM,0,0)=0.
        DO 280 IAEL=1,NAELI(IATM)
          LGAEL(IATM,0,0)=LGAEL(IATM,0,0)+LGAEL(IATM,IAEL,0)
280     CONTINUE
C
300   CONTINUE
C
C   GENERAL HEAVY PARTICLE IMPACT COLLISIONS
C
      CALL XSTAPI(COUN,PLS)
C
C
      DO 1000 IATM=1,NATMI
C
        IF (TRCAMD) THEN
          CALL MASBOX ('ATOMIC SPECIES IATM = '//TEXTS(IATM))
          CALL LEER(1)
C
          IF (LGAEI(IATM,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 870 IAEI=1,NAEII(IATM)
              IREI=LGAEI(IATM,IAEI)
              CALL XSTEI_2(IREI)
870         CONTINUE
          ENDIF
C
C
          CALL LEER(2)
          IF (LGACX(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 890 IACX=1,NACXI(IATM)
              IRCX=LGACX(IATM,IACX,0)
              IPL =LGACX(IATM,IACX,1)
              CALL XSTCX_2(IRCX,IPL)
890         CONTINUE
          ENDIF
          IF (LGAEL(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 895 IAEL=1,NAELI(IATM)
              IREL=LGAEL(IATM,IAEL,0)
              IPL =LGAEL(IATM,IAEL,1)
              CALL XSTEL_2(IREL,IPL)
895         CONTINUE
          ENDIF
          CALL LEER(2)
          IF (LGAPI(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO GENERAL ION IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 885 IAPI=1,NAPII(IATM)
              IRPI=LGAPI(IATM,IAPI,0)
              IPLS=LGAPI(IATM,IAPI,1)
              CALL LEER(1)
              WRITE (6,*) 'GENERAL ION IMPACT REACTION NO. IRPI= ',IRPI
              CALL LEER(1)
              WRITE (6,*) 'INCIDENT BULK ION: IPLS:'
              WRITE (6,*) 'IPLS= ',TEXTS(NSPAMI+IPLS)
              CALL LEER(1)
              WRITE (6,*) 'ELECTRONS: PELPI,EELPI'
              WRITE (6,*) 'EL      ', PELPI(IRPI),0.D0
              CALL LEER(1)
              WRITE (6,*) 'BULK ION SECONDARIES:'
              IF (IPPLPI(IRPI,0).GT.0) THEN
                WRITE (6,*) 'BULK IONS: PPLPI,EPLPI '
                DO IP=1,IPPLPI(IRPI,0)
                  IPL=IPPLPI(IRPI,IP)
                  WRITE (6,*) TEXTS(NSPAMI+IPL),PPLPI(IRPI,IPL),
     .                                          EPLPI(IRPI,IPL)
                ENDDO
              ELSE
                WRITE (6,*) 'NONE'
              ENDIF
              CALL LEER(1)
              WRITE (6,*) 'TEST PARTICLE SECONDARIES:'
              IF (P2NPI(IRPI).EQ.0.D0) THEN
                WRITE (6,*) 'NONE'
              ENDIF
              DO IA=1,IPATPI(IRPI,0)
                IAT=IPATPI(IRPI,IA)
                ISP=IAT
                WRITE (6,*) 'ATOM     IATM= ',
     .                      TEXTS(ISP),PATPI(IRPI,IAT)
              ENDDO
              DO IM=1,IPMLPI(IRPI,0)
                IML=IPMLPI(IRPI,IM)
                ISP=NSPA+IML
                WRITE (6,*) 'MOLECULE IMOL= ',
     .                      TEXTS(ISP),PMLPI(IRPI,IML)
              ENDDO
              DO II=1,IPIOPI(IRPI,0)
                IIO=IPIOPI(IRPI,II)
                ISP=NSPAM+IIO
                WRITE (6,*) 'TEST ION IION= ',
     .                      TEXTS(ISP),PIOPI(IRPI,IIO)
              ENDDO
            CALL LEER(1)
885         CONTINUE
          ENDIF
        ENDIF
1000  CONTINUE
C
      RETURN
C
991   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
993   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTA.'
      WRITE (6,*) 'PARAMETER NRDS IS TOO SMALL.'
      WRITE (6,*) 'I.E. TOO MANY EI. REACTIONS REQUESTED. '
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTA.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IATM = ',IATM
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
999   CONTINUE
      WRITE (6,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IATM,IREL '
      WRITE (6,*) IATM,IREL,IPLS
      CALL EXIT
      END
C
C
      SUBROUTINE XSECTA_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMSOU
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: IPLS, NRC, IDSC, IATM, IAT, KK, IPL
C
      DO 100 IATM=1,NATMI
C
        IF (NRCA(IATM).EQ.0.AND.NCHARA(IATM).LE.2) THEN
C
C  DEFAULT H,D,T OR HE ELEC. IMP. IONIZATION MODEL
C
          DO 52 IPLS=1,NPLSI
            IF (NCHARP(IPLS).EQ.NCHARA(IATM).AND.
     .          NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .          NCHRGP(IPLS).EQ.1) THEN
              NRDS=NRDS+1
              GOTO 50
            ENDIF
52        CONTINUE
          GOTO 100
C
50        CONTINUE
C
C  NON DEFAULT ELEC. IMP. COLLISION MODEL,
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 90 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
            NRDS=NRDS+1
90        CONTINUE
        ENDIF
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:
C
      DO 200 IATM=1,NATMI
C
C   HYDROGENIC AND HELIUM DEFAULT MODEL 100 --- 140
C
        IF (NRCA(IATM).EQ.0) THEN
          DO 122 IPLS=1,NPLSI
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND IPLS
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 122
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS).AND.NCHRGP(IPLS).EQ.1)
     .          GOTO 125
124           CONTINUE
              GOTO 122
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
              NRCX=NRCX+1
C
            ENDIF
122       CONTINUE
C
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
C  NO CX MODEL DEFINED
        ELSE
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IATM=1,NATMI
C
C   AT PRESENT NO DEFAULT MODEL
C
        IF (NRCA(IATM).EQ.0) THEN
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 230 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
            NREL=NREL+1
C
C  SPECIAL TREATMENT: BGK COLLISIONS AMONGST TESTPARTICLES
            IF (IBGKA(IATM,NRC).NE.0) THEN
              IF (NPBGKA(IATM).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE

          NAELI(IATM)=IDSC
        ENDIF
C
300   CONTINUE
C
      CALL XSTAPI_PARAM
C
      RETURN
C
      END
C
C
      SUBROUTINE XSECTI
C
C  TABLE FOR REACTION RATES FOR TEST IONS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: COUN(0:9,NSTORDR), PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: FACTKK, CHRDF0, RMASS, ACCMAS, ACCINV, DSUB, DEIMIN,
     .          EHEAVY, EELEC, EBULK
      INTEGER :: ICOUNT, IA1, IP2, IPLS, ITEST, IIO, IRDS, IION, IDSC1,
     .           NRC, J, IPLS1, IPLS2, IATM, KK, IATM1, IATM2, ITYPB,
     .           ISPZB, III, IDSC, IREL, IBGK, IIDS, IERR, IMOL, IIEL,
     .           IIEI, IREI, IESTM, IFRST, ISCND, ISCDE, IPL, IICX,
     .           IDSC2, IDEZ, IRCX
      CHARACTER(8) :: TEXTS1, TEXTS2
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
10      CONTINUE
      END IF
C
C
C  SET TEST IONIC SPECIES ATOMIC AND MOLECULAR DATA;
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL.
C  FOR HYDROGENIC MOLECULE IONS ONLY THIS MEANS: ZERO MFP,
C  INSTANTANOUS DECAY INTO ATOMS OR BULK IONS
C  FOR ALL OTHER SPECIES: INFINITE MFP, I.E. NO COLLISIONS
C
C
      DO 100 IION=1,NIONI
        IDSC1=0
        LGIEI(IION,0)=0
C
        DO NRC=1,NRCI(IION)
          KK=IREACI(IION,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
C  H2+:
          IF (NMASSI(IION).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
C  HD+:
          ELSEIF (NMASSI(IION).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
C  D2+:
          ELSEIF (NMASSI(IION).EQ.4) THEN
C  TEST: D2+ OR HT+, USE TEXTS(IION)
            IF (INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
C  D2+ TEST ION IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0) THEN
C  HT+ TEST ION IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (6,*) 'TEST ION NO ',IION,' COULD NOT BE IDENTIFIED'
              WRITE (6,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT+:
          ELSEIF (NMASSI(IION).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
C  T2+:
          ELSEIF (NMASSI(IION).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS (MAY BE SPLITTED INTO 1A AND 1B)
          IF (IATM1.NE.IATM2) THEN
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
7000      ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGIEI(IION,IDSC1)=IRDS
          PATDS(IRDS,IA1)=PATDS(IRDS,IA1)+1.
          PPLDS(IRDS,IP2)=PPLDS(IRDS,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IRDS,NSPH+IA1)=P2ND(IRDS,NSPH+IA1)+1.
          EATDS(IRDS,IA1,1)=RMASSA(IA1)/ACCMAS
          EATDS(IRDS,IA1,2)=1./RMASSA(IA1)/ACCINV
          EPLDS(IRDS,      1)=RMASSP(IP2)/ACCMAS
          EPLDS(IRDS,      2)=1./RMASSP(IP2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IA1,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IA1,2)
          PELDS(IRDS)=0.
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-8,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 73 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)*FACTKK
73          CONTINUE
            EELDS1(IRDS,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 8.6 EV
            EHVDS1(IRDS,1:NSBOX)=8.6
          ELSE
            CALL CDEF (TEINL,1,1,-8,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-8) = CF(1:9,1)
            FACREA(-8) = LOG(FACTKK)
            NREAEI(IRDS) = -8
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -8
            NREAHV(IRDS) = -4
          END IF
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 7000
          ENDIF
C  SECOND PROCESS
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGIEI(IION,IDSC1)=IRDS
          PPLDS(IRDS,IPLS1)=PPLDS(IRDS,IPLS1)+1.
          PPLDS(IRDS,IPLS2)=PPLDS(IRDS,IPLS2)+1.
          EPLDS(IRDS,1)=1.0
          EPLDS(IRDS,2)=1.0
          PELDS(IRDS)=1.
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-9,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 71 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
71          CONTINUE
            EELDS1(IRDS,1:NSBOX)=-15.5
C  TRANSFERRED KINETIC ENERGY: 0.5 EV
            EHVDS1(IRDS,1:NSBOX)=0.5
          ELSE
            CALL CDEF (TEINL,1,1,-9,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-9) = CF(1:9,1)
            FACREA(-9) = 0._DP
            NREAEI(IRDS) = -9
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -9
            NREAHV(IRDS) = -5
          END IF
C  THIRD PROCESS
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGIEI(IION,IDSC1)=IRDS
          PATDS(IRDS,IATM1)=PATDS(IRDS,IATM1)+1.
          PATDS(IRDS,IATM2)=PATDS(IRDS,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IRDS,NSPH+IATM1)=P2ND(IRDS,NSPH+IATM1)+1.
          P2ND(IRDS,NSPH+IATM2)=P2ND(IRDS,NSPH+IATM2)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IRDS,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)+EATDS(IRDS,IATM2,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)+EATDS(IRDS,IATM2,2)
          PELDS(IRDS)=-1.
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-10,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 72 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
72          CONTINUE
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            EELDS1(IRDS,1:NSBOX)=-0.88*TEIN(1:NSBOX)
C  TRANSFERRED KINETIC ENERGY: INGOING ELECTRON ENERGY
            EHVDS1(IRDS,1:NSBOX)=0.88*TEIN(1:NSBOX)
          ELSE
            CALL CDEF (TEINL,1,1,-10,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-10) = CF(1:9,1)
            FACREA(-10) = 0._DP
            NREAEI(IRDS) = -10
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -10
            NREAHV(IRDS) = -6
          END IF
C
76        CONTINUE
C
          NIDSI(IION)=IDSC1
C
          MODCOL(1,2,NSPAM+IION,1)=1
          MODCOL(1,4,NSPAM+IION,1)=1
C
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=-NCHRGI(IION)
            IIO=NSPAM+IION
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EHEAVY=ESCD1I(IION,NRC)+ESCD2I(IION,NRC)
            EELEC=EELECI(IION,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IF (NREII.GT.NRDS) GOTO 993
            IRDS=NREII
            LGIEI(IION,IDSC1)=IRDS
            CALL XSTEI(RMASS,IRDS,IIO,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
90        CONTINUE
C
          NIDSI(IION)=IDSC1
C
        ENDIF
C
        NIDSIM(IION)=NIDSI(IION)-1
        LGIEI(IION,0)=NIDSI(IION)
C
100   CONTINUE

      DO 200 IION=1,NIONI
        IDSC2=0
        LGICX(IION,0,0)=0
        LGICX(IION,0,1)=0
C
C  NO DEFAULT CX RATES
C
        IF (NRCI(IION).EQ.0) THEN
          NICXI(IION)=0
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IPLS=IDEZ(IBULKI(IION,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGICX(IION,IDSC2,0)=IRCX
            LGICX(IION,IDSC2,1)=IPLS
            IIO=NSPAM+IION
            IPL=IPLS
            RMASS=RMASSI(IION)
            IFRST=ISCD1I(IION,NRC)
            ISCND=ISCD2I(IION,NRC)
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL XSTCX(RMASS,IRCX,IIO,IPL,
     .                 IFRST,ISCND,EBULK,ISCDE,IESTM,KK,FACTKK)
C
130       CONTINUE
C
          NICXI(IION)=IDSC2
        ENDIF
C
        NICXIM(IION)=NICXI(IION)-1
        LGICX(IION,0,0)=0
        DO IICX=1,NICXI(IION)
          LGICX(IION,0,0)=LGICX(IION,0,0)+LGICX(IION,IICX,0)
        ENDDO
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IION=1,NIONI
        IDSC=0
        LGIEL(IION,0,0)=0
        LGIEL(IION,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
          NIELI(IION)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKI(IION,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGIEL(IION,IDSC,0)=IREL
            LGIEL(IION,IDSC,1)=IPLS

C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKI(IION)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKI(IION)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKI(IION,NRC),1,3)
              ISPZB=IDEZ(IBGKI(IION,NRC),3,3)
              IF (ITYPB.NE.3.OR.ISPZB.NE.IION) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKI(IION,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C
C
            III=NSPAM+IION
            IPL=IPLS
            ISCDE=ISCDEI(IION,NRC)
            IESTM=IESTMI(IION,NRC)
            EBULK=EBULKI(IION,NRC)
            CALL XSTEL(IREL,III,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE
C
          NIELI(IION)=IDSC
        ENDIF
C
        NIELIM(IION)=NIELI(IION)-1
C
        LGIEL(IION,0,0)=0.
        DO 280 IIEL=1,NIELI(IION)
          LGIEL(IION,0,0)=LGIEL(IION,0,0)+LGIEL(IION,IIEL,0)
280     CONTINUE
C
300   CONTINUE
C
      DO 1000 IION=1,NIONI
C
        DO 500 IIEI=1,NIDSI(IION)
          IREI=LGIEI(IION,IIEI)
          CALL XSTEI_1(IREI)
500     CONTINUE
C
        IF (TRCAMD) THEN
          CALL MASBOX ('TEST ION SPECIES IION = '//TEXTS(NSPAM+IION))
          CALL LEER(1)
          IF (LGICX(IION,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 215 IICX=1,NICXI(IION)
              IRCX=LGICX(IION,IICX,0)
              IPL =LGICX(IION,IICX,1)
              CALL XSTCX_2(IRCX,IPL)
215         CONTINUE
          ENDIF
C
          IF (LGIEI(IION,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELECTRON IMPACT COLLISION'
            CALL LEER(1)
          ELSE
            DO 210 IIDS=1,NIDSI(IION)
              IREI=LGIEI(IION,IIDS)
              CALL XSTEI_2(IREI)
210         CONTINUE
          ENDIF
          CALL LEER(1)
          IF (LGIEL(IION,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 815 IIEL=1,NIELI(IION)
              IREL=LGIEL(IION,IIEL,0)
              IPL =LGIEL(IION,IIEL,1)
              CALL XSTEL_2(IREL,IPL)
815         CONTINUE
          ENDIF
        ENDIF
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE'
      CALL EXIT
991   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IION,IPLS ',KK,IION,IPLS
993   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (6,*) 'PARAMETER NRDS IS TOO SMALL.'
      WRITE (6,*) 'I.E. TOO MANY DISS. REACTIONS REQUESTED. '
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IION = ',IION
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (6,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (6,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      WRITE (6,*) 'ERROR FLAG = ',IERR
      CALL EXIT
999   CONTINUE
      WRITE (6,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IMOL,IREL '
      WRITE (6,*) IMOL,IREL,IPLS
      CALL EXIT
      RETURN
      END
C
C
      SUBROUTINE XSECTI_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ICOUNT, NRC, KK, IION, IPLS, ITEST, IATM, IATM1,
     .           IATM2, IDSC2, IPLS1, IPLS2
C
      DO 100 IION=1,NIONI
C
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
C  H2+:
          IF (NMASSI(IION).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHARP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
C  HD+:
          ELSEIF (NMASSI(IION).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
C  D2+:
          ELSEIF (NMASSI(IION).EQ.4) THEN
C  TEST: D2+ OR HT+, USE TEXTS(IION)
            IF (INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
C  D2+ TEST ION IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0) THEN
C  HT+ TEST ION IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (6,*) 'TEST ION NO ',IION,' COULD NOT BE IDENTIFIED'
              WRITE (6,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT+:
          ELSEIF (NMASSI(IION).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
C  T2+:
          ELSEIF (NMASSI(IION).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS (MAY BE SPLITTED INTO 1A AND 1B)
          IF (IATM1.NE.IATM2) THEN
            ICOUNT=1
          ELSE
            ICOUNT=2
          ENDIF
C
7000      CONTINUE
          NRDS=NRDS+1
          IF (ICOUNT.EQ.1) THEN
            ICOUNT=2
            GOTO 7000
          ENDIF
C  SECOND PROCESS
          NRDS=NRDS+1
C  THIRD PROCESS
          NRDS=NRDS+1
C
76        CONTINUE
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            NRDS=NRDS+1
90        CONTINUE
C
        ENDIF
C
100   CONTINUE

      DO 200 IION=1,NIONI
C
C  NO DEFAULT CX RATES
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
          NICXI(IION)=IDSC2
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IION=1,NIONI
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            NREL=NREL+1

C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE
C
        ENDIF
C
300   CONTINUE
C
      RETURN
C
      END
C
C
      SUBROUTINE XSECTM
C
C  TABLE FOR REACTION RATES FOR MOLECULES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: ACCMAS, ACCINV, FACTKK, DEIMIN, EELEC, CHRDF0, DSUB,
     .          RMASS, EBULK, EHEAVY
      INTEGER :: ITEST, IATM, IPLS, IION, IA1, IP2, ION, IRDS, ICOUNT,
     .           IION3, IDSC1, NRC, KK, J, IMOL, IPLS1, IPLS2, IPLS3,
     .           IATM1, IATM2, ITYPB, ISPZB, IMEL, IDSC, IREL, IBGK,
     .           IMDS, IERR, IMCX, ISCND, ISCDE, IESTM, IML, IFRST,
     .           IRCX, IPL, IDSC2, IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
10      CONTINUE
      END IF
C
C
C  SET MOLECULAR DATA
C
C  STORE "DEFAULT DISSOCIATION MODEL" DATA
C  IN EACH CELL. THIS DEFAULT MODEL (JANEV/LANGER, PPPL) MAY BE USED
C  FOR HYDROGENIC MOLECULES ONLY.
C
C
      DO 100 IMOL=1,NMOLI
        IDSC1=0
        LGMEI(IMOL,0)=0
C
        DO NRC=1,NRCM(IMOL)
          KK=IREACM(IMOL,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
        IF (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
C  TO ALL HYDROGENIC MOLECULES
C
          EELEC=-EIONH2
C
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
          IPLS3=0
          IION3=0
C  H2:
          IF (NMASSM(IMOL).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
            DO 25 IION=1,NIONI
              IF (NMASSI(IION).EQ.2.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
25          CONTINUE
C  HD:
          ELSEIF (NMASSM(IMOL).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
            DO 35 IION=1,NIONI
              IF (NMASSI(IION).EQ.3.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
35          CONTINUE
C  D2:  (OR HT ?)
          ELSEIF (NMASSM(IMOL).EQ.4) THEN
C  TEST: D2 OR HT, USE TEXTS(IMOL)
            IF (INDEX(TEXTS(NSPA+IMOL),'D').NE.0) THEN
C  D2 MOLECULE IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
              DO 45 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .              INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
                  IION3=IION
                ENDIF
45            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPA+IMOL),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPA+IMOL),'T').NE.0) THEN
C  HT MOLECULE IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
              DO 48 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .             (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0)) THEN
                  IION3=IION
                ENDIF
48            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (6,*) 'MOLECULE NO ',IMOL,' COULD NOT BE IDENTIFIED'
              WRITE (6,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT:
          ELSEIF (NMASSM(IMOL).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
            DO 55 IION=1,NIONI
              IF (NMASSI(IION).EQ.5.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
55          CONTINUE
C  T2:
          ELSEIF (NMASSM(IMOL).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
            DO 65 IION=1,NIONI
              IF (NMASSI(IION).EQ.6.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
65          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2*IION3
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGMEI(IMOL,IDSC1)=IRDS
          PATDS(IRDS,IATM1)=PATDS(IRDS,IATM1)+1.
          PATDS(IRDS,IATM2)=PATDS(IRDS,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IRDS,NSPH+IATM1)=P2ND(IRDS,NSPH+IATM1)+1.
          P2ND(IRDS,NSPH+IATM2)=P2ND(IRDS,NSPH+IATM2)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IRDS,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)+EATDS(IRDS,IATM2,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)+EATDS(IRDS,IATM2,2)
          PELDS(IRDS)=0.
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-5,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 70 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
70          CONTINUE
            EELDS1(IRDS,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 6 EV
            EHVDS1(IRDS,1:NSBOX)=6.
          ELSE
            CALL CDEF (TEINL,1,1,-5,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-5) = CF(1:9,1)
            FACREA(-5) = 0._DP
            NREAEI(IRDS)=-5
            JEREAEI(IRDS)=1
            NELREI(IRDS)=-5
            NREAHV(IRDS)=-2
          END IF
C  SECOND PROCESS (MAY BE SPLITTED INTO 2A AND 2B)
          IF (IATM1.NE.IATM2) THEN
            FACTKK=0.5
            ICOUNT=1
          ELSE
            FACTKK=1.D0
            ICOUNT=2
          ENDIF
C
          IA1=IATM1
          IP2=IPLS2
73        ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGMEI(IMOL,IDSC1)=IRDS
          PATDS(IRDS,IA1)=PATDS(IRDS,IA1)+1.
          PPLDS(IRDS,IP2)=PPLDS(IRDS,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IRDS,NSPH+IA1)=P2ND(IRDS,NSPH+IA1)+1.
          EATDS(IRDS,IA1,1)=RMASSA(IA1)/ACCMAS
          EATDS(IRDS,IA1,2)=1./RMASSA(IA1)/ACCINV
          EPLDS(IRDS,      1)=RMASSP(IP2)/ACCMAS
          EPLDS(IRDS,      2)=1./RMASSP(IP2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IA1,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IA1,2)
          PELDS(IRDS)=1.0
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-6,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 71 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)*FACTKK
71          CONTINUE
            EELDS1(IRDS,1:NSBOX)=-25.0
C  TRANSFERRED KINETIC ENERGY: 10 EV
            EHVDS1(IRDS,1:NSBOX)=10.0
          ELSE
            CALL CDEF (TEINL,1,1,-6,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-6) = CF(1:9,1)
            FACREA(-6) = LOG(FACTKK)
            NREAEI(IRDS) = -6
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -6
            NREAHV(IRDS) = -3
          END IF
          IF (ICOUNT.EQ.1) THEN
            IA1=IATM2
            IP2=IPLS1
            ICOUNT=2
            GOTO 73
          ENDIF
C
C  THIRD PROCESS
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGMEI(IMOL,IDSC1)=IRDS
          ION=NSPAM+IION3
          PIODS(IRDS,IION3)=PIODS(IRDS,IION3)+1.
          P2ND(IRDS,ION)=P2ND(IRDS,ION)+1.
          EIODS(IRDS,IION3,1)=1.
          EIODS(IRDS,IION3,2)=0.
          EIODS(IRDS,0,1)=1.
          EIODS(IRDS,0,2)=0.
          PELDS(IRDS)=1.0
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,-7,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            DO 72 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
72          CONTINUE
            EELDS1(IRDS,1:NSBOX)=EELEC
          ELSE
            CALL CDEF (TEINL,1,1,-7,COUN,1,CF,.TRUE.,.FALSE.,
     .                 .TRUE.)
            CREAC(1:9,1,-7) = CF(1:9,1)
            FACREA(-7) = 0._DP
            NREAEI(IRDS) = -7
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -7
            EELDS1(IRDS,1)=EELEC
          END IF
C
76        CONTINUE
          NMDSI(IMOL)=IDSC1
C
          MODCOL(1,2,NSPA+IMOL,1)=1
          MODCOL(1,4,NSPA+IMOL,1)=1
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.
            IML=NSPA+IMOL
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EHEAVY=ESCD1M(IMOL,NRC)+ESCD2M(IMOL,NRC)
            EELEC=EELECM(IMOL,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IF (NREII.GT.NRDS) GOTO 993
            IRDS=NREII
            LGMEI(IMOL,IDSC1)=IRDS
            CALL XSTEI(RMASS,IRDS,IML,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,
     .                 ISCDE,EELEC,IESTM,KK,FACTKK,PLS)
90        CONTINUE
          NMDSI(IMOL)=IDSC1
        ENDIF
C
        NMDSIM(IMOL)=NMDSI(IMOL)-1
        LGMEI(IMOL,0)=NMDSI(IMOL)
100   CONTINUE
C
C
C   CHARGE EXCHANGE:

      DO 200 IMOL=1,NMOLI
        IDSC2=0
        LGMCX(IMOL,0,0)=0
        LGMCX(IMOL,0,1)=0
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMCXI(IMOL)=0
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 130 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IRCX=NRCXI
            LGMCX(IMOL,IDSC2,0)=IRCX
            LGMCX(IMOL,IDSC2,1)=IPLS
            IML=NSPA+IMOL
            IPL=IPLS
            RMASS=RMASSM(IMOL)
            IFRST=ISCD1M(IMOL,NRC)
            ISCND=ISCD2M(IMOL,NRC)
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL XSTCX(RMASS,IRCX,IML,IPL,
     .                 IFRST,ISCND,EBULK,ISCDE,IESTM,KK,FACTKK)
C
130       CONTINUE
C
          NMCXI(IMOL)=IDSC2
        ENDIF
C
        NMCXIM(IMOL)=NMCXI(IMOL)-1
C
        LGMCX(IMOL,0,0)=0
        DO 161 IMCX=1,NMCXI(IMOL)
          LGMCX(IMOL,0,0)=LGMCX(IMOL,0,0)+LGMCX(IMOL,IMCX,0)
161     CONTINUE
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IMOL=1,NMOLI
        IDSC=0
        LGMEL(IMOL,0,0)=0
        LGMEL(IMOL,0,1)=0
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMELI(IMOL)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 230 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC=IDSC+1
            NRELI=NRELI+1
            IREL=NRELI
            LGMEL(IMOL,IDSC,0)=IREL
            LGMEL(IMOL,IDSC,1)=IPLS
C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
                NRBGI=NRBGI+3
                IBGK=NRBGI/3
                NPBGKM(IMOL)=IBGK
              ENDIF
              IF (NPBGKP(IPLS,1).EQ.0) THEN
                NPBGKP(IPLS,1)=NPBGKM(IMOL)
              ELSE
                GOTO 999
              ENDIF
C  SELF OR CROSS COLLISION?
              ITYPB=IDEZ(IBGKM(IMOL,NRC),1,3)
              ISPZB=IDEZ(IBGKM(IMOL,NRC),3,3)
              IF (ITYPB.NE.2.OR.ISPZB.NE.IMOL) THEN
C  CROSS COLLISION !
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKM(IMOL,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
C
            IML=NSPA+IMOL
            IPL=IPLS
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL XSTEL(IREL,IML,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE
C
          NMELI(IMOL)=IDSC
        ENDIF
C
        NMELIM(IMOL)=NMELI(IMOL)-1
C
        LGMEL(IMOL,0,0)=0.
        DO 280 IMEL=1,NMELI(IMOL)
          LGMEL(IMOL,0,0)=LGMEL(IMOL,0,0)+LGMEL(IMOL,IMEL,0)
280     CONTINUE
C
C
300   CONTINUE
C
      DO 1000 IMOL=1,NMOLI
C
        DO 500 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          CALL XSTEI_1(IRDS)
500     CONTINUE
C
        IF (TRCAMD) THEN
          CALL MASBOX ('MOLECULAR SPECIES IMOL = '//TEXTS(NSPA+IMOL))
          CALL LEER(1)
          IF (LGMCX(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 880 IMCX=1,NMCXI(IMOL)
              IRCX=LGMCX(IMOL,IMCX,0)
              IPL =LGMCX(IMOL,IMCX,1)
              CALL XSTCX_2(IRCX,IPL)
880         CONTINUE
          ENDIF
C
          IF (LGMEI(IMOL,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 870 IMDS=1,NMDSI(IMOL)
              IRDS=LGMEI(IMOL,IMDS)
              CALL XSTEI_2(IRDS)
870         CONTINUE
          ENDIF
C
          IF (LGMEL(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (6,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 895 IMEL=1,NMELI(IMOL)
              IREL=LGMEL(IMOL,IMEL,0)
              IPL =LGMEL(IMOL,IMEL,1)
              CALL XSTEL_2(IREL,IPL)
895         CONTINUE
          ENDIF
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE '
      CALL EXIT
991   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IMOL,IPLS ',KK,IMOL,IPLS
993   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTM.'
      WRITE (6,*) 'PARAMETER NRDS IS TOO SMALL.'
      WRITE (6,*) 'I.E. TOO MANY DISS. REACTIONS REQUESTED. '
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTM.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IMOL = ',IMOL
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (6,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (6,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (6,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      WRITE (6,*) 'ERROR FLAG = ',IERR
      CALL EXIT
999   CONTINUE
      WRITE (6,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IMOL,IREL '
      WRITE (6,*) IMOL,IREL,IPLS
      CALL EXIT
      RETURN
C
      END
C
C
      SUBROUTINE XSECTM_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: ITEST, ICOUNT, NRC, IION, IMOL, IATM, IPLS, IATM1,
     .           IATM2, IPLS1, IPLS2, IPLS3, IION3, KK
C
      DO 100 IMOL=1,NMOLI
C
        IF (NRCM(IMOL).EQ.0.AND.NCHARM(IMOL).EQ.2) THEN
C  APPLY THE DEFAULT MODEL FOR H2 DISSOCIATION AND IONIZATION
C  TO ALL HYDROGENIC MOLECULES
C
C  FIRST: FIND SECONDARY SPECIES INDICES:
          IATM1=0
          IATM2=0
          IPLS1=0
          IPLS2=0
          IPLS3=0
          IION3=0
C  H2:
          IF (NMASSM(IMOL).EQ.2) THEN
            DO 21 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
21          CONTINUE
            DO 23 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
23          CONTINUE
            DO 25 IION=1,NIONI
              IF (NMASSI(IION).EQ.2.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
25          CONTINUE
C  HD:
          ELSEIF (NMASSM(IMOL).EQ.3) THEN
            DO 31 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.1) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.2) THEN
                IATM2=IATM
              ENDIF
31          CONTINUE
            DO 33 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
33          CONTINUE
            DO 35 IION=1,NIONI
              IF (NMASSI(IION).EQ.3.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
35          CONTINUE
C  D2:  (OR HT ?)
          ELSEIF (NMASSM(IMOL).EQ.4) THEN
C  TEST: D2 OR HT, USE TEXTS(IMOL)
            IF (INDEX(TEXTS(NSPA+IMOL),'D').NE.0) THEN
C  D2 MOLECULE IDENTIFIED
              DO 41 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.2) THEN
                  IATM1=IATM
                  IATM2=IATM
                ENDIF
41            CONTINUE
              DO 43 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                  IPLS2=IPLS
                ENDIF
43            CONTINUE
              DO 45 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .              INDEX(TEXTS(NSPAM+IION),'D').NE.0) THEN
                  IION3=IION
                ENDIF
45            CONTINUE
            ELSEIF (INDEX(TEXTS(NSPA+IMOL),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPA+IMOL),'T').NE.0) THEN
C  HT MOLECULE IDENTIFIED
              DO 46 IATM=1,NATMI
                IF (NMASSA(IATM).EQ.1) THEN
                  IATM1=IATM
                ELSEIF (NMASSA(IATM).EQ.3) THEN
                  IATM2=IATM
                ENDIF
46            CONTINUE
              DO 47 IPLS=1,NPLSI
                IF (NMASSP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS1=IPLS
                ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                  IPLS2=IPLS
                ENDIF
47            CONTINUE
              DO 48 IION=1,NIONI
                IF (NMASSI(IION).EQ.4.AND.NCHARI(IION).EQ.2.AND.
     .             (INDEX(TEXTS(NSPAM+IION),'H').NE.0.OR.
     .              INDEX(TEXTS(NSPAM+IION),'T').NE.0)) THEN
                  IION3=IION
                ENDIF
48            CONTINUE
            ELSE
              CALL LEER(2)
              WRITE (6,*) 'MOLECULE NO ',IMOL,' COULD NOT BE IDENTIFIED'
              WRITE (6,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL LEER(2)
            ENDIF
C  DT:
          ELSEIF (NMASSM(IMOL).EQ.5) THEN
            DO 51 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.2) THEN
                IATM1=IATM
              ELSEIF (NMASSA(IATM).EQ.3) THEN
                IATM2=IATM
              ENDIF
51          CONTINUE
            DO 53 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.2.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
              ELSEIF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS2=IPLS
              ENDIF
53          CONTINUE
            DO 55 IION=1,NIONI
              IF (NMASSI(IION).EQ.5.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
55          CONTINUE
C  T2:
          ELSEIF (NMASSM(IMOL).EQ.6) THEN
            DO 61 IATM=1,NATMI
              IF (NMASSA(IATM).EQ.3) THEN
                IATM1=IATM
                IATM2=IATM
              ENDIF
61          CONTINUE
            DO 63 IPLS=1,NPLSI
              IF (NMASSP(IPLS).EQ.3.AND.NCHRGP(IPLS).EQ.1) THEN
                IPLS1=IPLS
                IPLS2=IPLS
              ENDIF
63          CONTINUE
            DO 65 IION=1,NIONI
              IF (NMASSI(IION).EQ.6.AND.NCHARI(IION).EQ.2) THEN
                IION3=IION
              ENDIF
65          CONTINUE
          ENDIF
          ITEST=IATM1*IATM2*IPLS1*IPLS2*IION3
          IF (ITEST.EQ.0) GOTO 76
C
C  SET DEFAULT MODEL: 3 ELECTRON IMPACT PROCESSES
C
C  FIRST PROCESS
          NRDS=NRDS+1
C  SECOND PROCESS (MAY BE SPLITTED INTO 2A AND 2B)
          IF (IATM1.NE.IATM2) THEN
            ICOUNT=1
          ELSE
            ICOUNT=2
          ENDIF
C
73        CONTINUE
          NRDS=NRDS+1
          IF (ICOUNT.EQ.1) THEN
            ICOUNT=2
            GOTO 73
          ENDIF
C
C  THIRD PROCESS
          NRDS=NRDS+1
C
76        CONTINUE
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            NRDS=NRDS+1
90        CONTINUE
        ENDIF
C
100   CONTINUE
C
C
C   CHARGE EXCHANGE:

      DO 200 IMOL=1,NMOLI
C
        IF (NRCM(IMOL).EQ.0) THEN
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 130 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.3) GOTO 130
C
            NRCX=NRCX+1
C
130       CONTINUE
C
        ENDIF
C
200   CONTINUE
C
C   ELASTIC COLLISIONS
C
      DO 300 IMOL=1,NMOLI
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCM(IMOL).EQ.0) THEN
          NMELI(IMOL)=0
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 230 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            NREL=NREL+1
C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE
C
        ENDIF
C
C
C
300   CONTINUE
C
      RETURN
C
      END
C
C
      SUBROUTINE XSECTP
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI
      IMPLICIT NONE
C
      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: DELE, FCTKKL, EEMAX, ZX, DSUB, DEIMIN, RMASS2, FACTKK
      INTEGER :: IIRC, IION3, IPLS3, IATM3, IMOL3, KK, NRC, IATM,
     .           IRRC, J, IDSC, IPLS, NSERC5, KREAD, I, MODC, IATM1,
     .           ITYP, ISPZ, IDEZ
      SAVE
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
70      CONTINUE
      END IF
C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
C
        IDSC=0
        LGPRC(IPLS,0)=0
C
        DO NRC=1,NRCP(IPLS)
          KK=IREACP(IPLS,NRC)
          IF (ISWR(KK).LE.0.OR.ISWR(KK).GT.6) GOTO 994
        ENDDO
C
        IF (NRCP(IPLS).EQ.0) THEN
C
          IF (NCHARP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL
C  HYDR. RECOMBINATION RATE-COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            DO 52 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .                            NCHRGP(IPLS).EQ.1) THEN
C
                IDSC=IDSC+1
                NRRCI=NRRCI+1
                IF (NRRCI.GT.NREC) GOTO 992
                IRRC=NRRCI
                LGPRC(IPLS,IDSC)=IRRC
                IF (NSTORDR >= NRAD) THEN
                DO 51 J=1,NSBOX
                  ZX=EIONH/MAX(1.E-5_DP,TEIN(J))
                  TABRC1(IRRC,J)=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(J)
c slmod begin (sl)
                  IF (IOPT1.EQ.1) TABRC1(IRRC,J)=TABRC1(IRRC,J)*
     .                                           TABRCM(IRRC,J)
c slmod end
                  EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
51              CONTINUE
                ELSE
                  NREARC(IRRC) = 0
                  JEREARC(IRRC) = 0
                  NELRRC(IRRC) = -1
                END IF
                IATM1=IATM
                NATPRC(IRRC)=IATM1
                NIOPRC(IRRC)=0
                NPLPRC(IRRC)=0
                NMLPRC(IRRC)=0
              ENDIF
52          CONTINUE
C
            MODCOL(6,2,NSPAMI+IPLS,1)=1
            MODCOL(6,4,NSPAMI+IPLS,1)=1
C
            NPRCI(IPLS)=IDSC
          ENDIF
C
C  NON DEFAULT MODEL:  240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
            IF (ISWR(KK).NE.6) GOTO 82
C
            FACTKK=FREACP(IPLS,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  RECOMBINATION MODEL FOR BULK IONS
            IDSC=IDSC+1
            NRRCI=NRRCI+1
            IF (NRRCI.GT.NREC) GOTO 992
            IRRC=NRRCI
            LGPRC(IPLS,IDSC)=IRRC
C
            ITYP=IDEZ(ISCD1P(IPLS,NRC),1,3)
            ISPZ=IDEZ(ISCD1P(IPLS,NRC),3,3)
            IF (ITYP.EQ.3) THEN
              NIOPRC(IRRC)=ISPZ
              RMASS2=RMASSI(ISPZ)
            ELSEIF (ITYP.EQ.4) THEN
              NPLPRC(IRRC)=ISPZ
              RMASS2=RMASSP(ISPZ)
            ELSEIF (ITYP.EQ.1) THEN
              NATPRC(IRRC)=ISPZ
              RMASS2=RMASSA(ISPZ)
            ELSEIF (ITYP.EQ.2) THEN
              NMLPRC(IRRC)=ISPZ
              RMASS2=RMASSM(ISPZ)
            ENDIF
C  CHECK MASS CONSERVATION
            IF (RMASSP(IPLS).NE.RMASS2) GOTO 993
C
C  1.) CROSS SECTION(TE)
C           NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3)
C
C  2.A) RATE COEFFICIENT = CONST.
C           TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
            IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                     .TRUE.)
                DO 92 J=1,NSBOX
                  TABRC1(IRRC,J)=TABRC1(IRRC,J)+
     +                                COUN(1,J)*DEIN(J)*FACTKK
c slmod begin (sl)
c...              FIX: Do this multiplication below somewhere?  
                  IF (IOPT1.EQ.1) TABRC1(IRRC,J)=TABRC1(IRRC,J)*
     .                                           TABRCM(IRRC,J)
c slmod end
92              CONTINUE
              ELSE
                NREARC(IRRC) = KK
                JEREARC(IRRC) = 1
                FACREA(KK) = LOG(FACTKK)
              END IF
              MODCOL(6,2,NSPAMI+IPLS,1)=1
C           ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
            ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
              IF (NSTORDR >= NRAD) THEN
                CALL CDEFN(TEINL,PLS,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     ,                     .TRUE.)
                DO 93 J=1,NSBOX
                  TABRC1(IRRC,J)=COUN(1,J)*DEIN(J)*FACTKK
c slmod begin (sl)
                  IF (IOPT1.EQ.1) TABRC1(IRRC,J)=TABRC1(IRRC,J)*
     .                                           TABRCM(IRRC,J)
c slmod end
93              CONTINUE
              ELSE
                FACREA(KK) = FACTKK
                NREARC(IRRC) = KK
                JEREARC(IRRC) = 2
              END IF
              MODCOL(6,2,NSPAMI+IPLS,1)=1
            ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ELECTRON ENERGY LOSS RATE
C
            NSERC5=IDEZ(ISCDEP(IPLS,NRC),5,5)
            IF (NSERC5.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 101 J=1,NSBOX
                  EELRC1(IRRC,J)=EELECP(IPLS,NRC)*TABRC1(IRRC,J)
101             CONTINUE
              ELSE
                NELRRC(IRRC) = -2
                EELRC1(IRRC,1)=EELECP(IPLS,NRC)
              END IF
              MODCOL(6,4,NSPAMI+IPLS,1)=1
            ELSEIF (NSERC5.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ELECTRON = -1.5*TE*RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 102 J=1,NSBOX
                  EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
102             CONTINUE
              ELSE
                NELRRC(IRRC) = -3
              END IF
              MODCOL(6,4,NSPAMI+IPLS,1)=1
            ELSEIF (NSERC5.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
              KREAD=EELECP(IPLS,NRC)
              MODC=IDEZ(MODCLF(KREAD),5,5)
              IF (MODC.EQ.1) THEN
                IF (NSTORDR >= NRAD) THEN
                  CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                       .FALSE.,.TRUE.)
                  DO 104 J=1,NSBOX
                    EELRC1(IRRC,J)=-COUN(1,J)*DEIN(J)*FACTKK
104               CONTINUE
                ELSE
                  NELRRC(IRRC)=KREAD
                  JELRRC(IRRC)=1
                  FACREA(KREAD) = LOG(FACTKK)
                END IF
                MODCOL(6,4,NSPAMI+IPLS,1)=1
C  4.D)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C             ELSEIF (MODC.EQ.2) THEN
C        IRRELEVANT
C               MODCOL(6,4,NSPAMI+IPLS,1)=2
C  4.E)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
              ELSEIF (MODC.EQ.3) THEN
                IF (NSTORDR >= NRAD) THEN
                  CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                       .FALSE.,.TRUE.)
                  DO 106 J=1,NSBOX
                    EELRC1(IRRC,J)=COUN(9,J)
106               CONTINUE
                  DO 108 I=8,1,-1
                    DO 107 J=1,NSBOX
                      EELRC1(IRRC,J)=EELRC1(IRRC,J)*PLS(J)+
     .                                    COUN(I,J)
107                 CONTINUE
108               CONTINUE
                  FCTKKL=LOG(FACTKK)
                  DO 109 J=1,NSBOX
                    EEMAX=MAX(-100._DP,EELRC1(IRRC,J)+FCTKKL+DEINL(J))
                    EELRC1(IRRC,J)=-EXP(EEMAX)
c slmod begin (sl)
                    IF (IOPT1.EQ.1) EELRC1(IRRC,J)=EELRC1(IRRC,J)*
     .                                             TABREM(IRRC,J)
c slmod end
109               CONTINUE
                ELSE
                  NELRRC(IRRC)=KREAD
                  JELRRC(IRRC)=9
                  FACREA(KREAD) = LOG(FACTKK)
                END IF
                MODCOL(6,4,NSPAMI+IPLS,1)=1
              ENDIF
              IF (DELPOT(KREAD).NE.0.D0) THEN
                DELE=DELPOT(KREAD)
                IF (NSTORDR >= NRAD) THEN
                  DO 110 J=1,NSBOX
                    EELRC1(IRRC,J)=EELRC1(IRRC,J)+
     .                                  DELE*TABRC1(IRRC,J)
110               CONTINUE
                END IF
              ENDIF
            ENDIF
C
82        CONTINUE
          NPRCI(IPLS)=IDSC
C
C  NO MODEL DEFINED
        ELSE
          NPRCI(IPLS)=0
        ENDIF
C
        NPRCIM(IPLS)=NPRCI(IPLS)-1
        LGPRC(IPLS,0)=NPRCI(IPLS)
C
        IF (TRCAMD) THEN
          CALL MASBOX ('BULK ION SPECIES IPLS = '//TEXTS(NSPAMI+IPLS))
          CALL LEER(1)
          IF (LGPRC(IPLS,0).EQ.0) THEN
            WRITE (6,*) 'NO RECOMBINATION '
          ELSE
            DO 220 IIRC=1,NPRCI(IPLS)
              IRRC=LGPRC(IPLS,IIRC)
              WRITE (6,*) 'RECOMBINATION NO. IRRC= ',IRRC
              WRITE (6,*) 'RECOMBINATION INTO SPECIES:'
              IION3=NIOPRC(IRRC)
              IF (IION3.NE.0) WRITE (6,*) 'TEST ION IION= ',
     .                                     TEXTS(NSPAM+IION3)
              IPLS3=NPLPRC(IRRC)
              IF (IPLS3.NE.0) WRITE (6,*) 'BULK ION IPLS= ',
     .                                     TEXTS(NSPAMI+IPLS3)
              IATM3=NATPRC(IRRC)
              IF (IATM3.NE.0) WRITE (6,*) 'ATOM     IATM= ',
     .                                     TEXTS(NSPH+IATM3)
              IMOL3=NMLPRC(IRRC)
              IF (IMOL3.NE.0) WRITE (6,*) 'MOLECULE IMOL= ',
     .                                     TEXTS(NSPA+IMOL3)
C             WRITE (6,*) 'ELECTRONS: PELPRC,EELRC1'
C             IF (NSTORDR >= NRAD) THEN
C               WRITE (6,*) 'EL      ',1.,EELRC1(IRRC,1)
C             ELSE
C               WRITE (6,*) 'EL      ',1.,FEELRC1(IRRC,1)
C             END IF
220         CONTINUE
          ENDIF
          CALL LEER(1)
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR RECOMBINATION'
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'NREC TOO SMALL, CHECK PARAMETER STATEMENTS'
      CALL EXIT
993   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTP: EXIT CALLED  '
      WRITE (6,*) 'MASS CONSERVATION VIOLATED, IPLS,IRRC ',IPLS,IRRC
      CALL EXIT
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTP.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IPLS = ',IPLS
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
C
      END
C
C
      SUBROUTINE XSECTPH
C
C  TABLE FOR REACTION RATES FOR PHOTONS
C
      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE

      RETURN
      END SUBROUTINE XSECTPH
C
C
      SUBROUTINE XSECTPH_PARAM

      USE PRECISION
      USE PARMMOD

      IMPLICIT NONE


      RETURN

      END SUBROUTINE XSECTPH_PARAM
C
C
      SUBROUTINE XSECTP_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: NRC, IATM, IPLS, KK
C
C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
C
        IF (NRCP(IPLS).EQ.0) THEN
C
          IF (NCHARP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL
C  HYDR. RECOMBINATION RATE-COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            DO 52 IATM=1,NATMI
              IF (NMASSP(IPLS).EQ.NMASSA(IATM).AND.
     .                            NCHRGP(IPLS).EQ.1) THEN
C
                NREC=NREC+1
              ENDIF
52          CONTINUE
C
          ENDIF
C
C  NON DEFAULT MODEL:  240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
            IF (ISWR(KK).NE.6) GOTO 82
C
C  RECOMBINATION MODEL FOR BULK IONS
            NREC=NREC+1
C
82        CONTINUE
C
C  NO MODEL DEFINED
        ELSE
        ENDIF
C
C
1000  CONTINUE
C
      RETURN
C
      END
C
C
      SUBROUTINE XSTAPI(COUN,PLS)
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), INTENT(OUT) :: PLS(*), COUN(0:9,*)
      REAL(DP) :: CF(9,0:9)
      REAL(DP) :: FCTKKL, P2N, FACTKK, TMASS, ADDT, ADDTL, PMASS, CHRDIF
      INTEGER :: NSEPI4, NSEPI5, IAPI, I, NEND, J, IAT, IO, ION, IA,
     .           IML, IM, MODC, KK, IDEZ, NRC, IDSC, IRPI, IIO, IPL
      SAVE
C
C   ION IMPACT COLLISIONS
C
      DO 1000 IATM=1,NATMI
        IDSC=0
        LGAPI(IATM,0,0)=0
        LGAPI(IATM,0,1)=0
C
C  NO DEFAULT MODEL
C

        IF (NRCA(IATM).EQ.0) THEN
          NAPII(IATM)=0
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) GOTO 130
C
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
C  INCIDENT BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 990
            IDSC=IDSC+1
            NRPII=NRPII+1
            IF (NRPII.GT.NRPI) GOTO 999
            IRPI=NRPII
            NREAPI(IRPI) = KK
            LGAPI(IATM,IDSC,0)=IRPI
            LGAPI(IATM,IDSC,1)=IPLS
            PPLPI(IRPI,IPLS)=PPLPI(IRPI,IPLS)-1.D0
C  SECONDARY INDEX, FIRST SECONDARY
            ITYP=IDEZ(ISCD1A(IATM,NRC),1,3)
            ISPZ=IDEZ(ISCD1A(IATM,NRC),3,3)
            IF (ITYP.EQ.3) PIOPI(IRPI,ISPZ)=PIOPI(IRPI,ISPZ)+1.D0
            IF (ITYP.EQ.4) PPLPI(IRPI,ISPZ)=PPLPI(IRPI,ISPZ)+1.D0
C  SECONDARY INDEX, SECOND SECONDARY
            ITYP=IDEZ(ISCD2A(IATM,NRC),1,3)
            ISPZ=IDEZ(ISCD2A(IATM,NRC),3,3)
            IF (ITYP.EQ.3) PIOPI(IRPI,ISPZ)=PIOPI(IRPI,ISPZ)+1.D0
            IF (ITYP.EQ.4) PPLPI(IRPI,ISPZ)=PPLPI(IRPI,ISPZ)+1.D0
C
            CHRDIF=0.
            DO 133 IIO=1,NIONI
              CHRDIF=CHRDIF+PIOPI(IRPI,IIO)*NCHRGI(IIO)
133         CONTINUE
            DO 134 IPL=1,NPLSI
              CHRDIF=CHRDIF+PPLPI(IRPI,IPL)*NCHRGP(IPL)
134         CONTINUE
            PELPI(IRPI)=PELPI(IRPI)+CHRDIF
C
            PPLPI(IRPI,IPLS)=PPLPI(IRPI,IPLS)+1.D0
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
            PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT: TARGET AT REST)
            TMASS=MASST(KK)*PMASSA
C
            ADDT=PMASS/RMASSP(IPLS)
            ADDTL=LOG(ADDT)
            ADDPI(IRPI,IPLS) = ADDTL
C
C CROSS SECTION (E-LAB)
            IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
              MODCOL(4,1,NSPH+IATM,IPLS)=KK
              MODCOL(4,2,NSPH+IATM,IPLS)=3
              IF (FACTKK.NE.1.D0)
     .        WRITE (6,*) 'FREACA NOT READY FOR CROSS SECTION IN XSTAPI'
            ENDIF
C
C RATE COEFFICIENT
            MODC=IDEZ(MODCLF(KK),3,5)
            IF (MODC.GE.1.AND.MODC.LE.2) THEN
              MODCOL(4,2,NSPH+IATM,IPLS)=MODC
              IF (MODC.EQ.1) NEND=1
              IF (MODC.EQ.2) NEND=NSTORDT
              IF (NSTORDR >= NRAD) THEN
                DO 142 J=1,NSBOX
                  PLS(J)=TIINL(IPLS,J)+ADDTL
142             CONTINUE
                IF (MODC.EQ.1) THEN
                  CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.TRUE.,
     .                      .FALSE.,.TRUE.)
                  DO 145 J=1,NSBOX
                    TABPI3(IRPI,J,1)=COUN(1,J)*DIIN(IPLS,J)*FACTKK
145               CONTINUE
                ELSEIF (MODC.EQ.2) THEN
                  CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.FALSE.,
     .                      .FALSE.,.TRUE.)
                  FCTKKL=LOG(FACTKK)
                  DO 146 J=1,NSBOX
                    TABPI3(IRPI,J,1)=COUN(1,J)+DIINL(IPLS,J)+FCTKKL
146               CONTINUE
                ENDIF
                DO 144 I=2,NEND
                  DO 143 J=1,NSBOX
                    TABPI3(IRPI,J,I)=COUN(I,J)
143               CONTINUE
144             CONTINUE
              ELSE
                FACREA(KK) = LOG(FACTKK)
              END IF
            ENDIF
C
            DEFPI(IRPI)=LOG(CVELI2*PMASS)
            EEFPI(IRPI)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4A. BULK ION ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION
            NSEPI4=IDEZ(ISCDEA(IATM,NRC),4,5)
            IF (NSEPI4.EQ.0) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 151 J=1,NSBOX
                  EPLPI3(IRPI,J,1)=EBULKA(IATM,NRC)
151             CONTINUE
              ELSE
                NELRPI(IRPI)=-1
                EPLPI3(IRPI,1,1)=EBULKA(IATM,NRC)
              END IF
              MODCOL(4,4,NSPH+IATM,IPLS)=1
            ELSEIF (NSEPI4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 252 J=1,NSBOX
                  EPLPI3(IRPI,J,1)=1.5*TIIN(IPLS,J)+EDRIFT(IPLS,J)
252             CONTINUE
              ELSE
                NELRPI(IRPI) = -2
              END IF
              MODCOL(4,4,NSPH+IATM,IPLS)=1
            ELSE
              WRITE (6,*) 'NSEPI4 ILL DEFINED IN XSTAPI '
              CALL EXIT
            ENDIF
C
C  4B. BULK ELECTRON ENERGY LOSS RATE
C
C  SET NET ENERGY LOSS RATE OF ELECTRON (IF ANY INVOLVED)
            NSEPI5=IDEZ(ISCDEA(IATM,NRC),5,5)
            IF (NSEPI5.EQ.0) THEN
C             MODCOL(4,4,NSPH+IATM,IPLS)=1
            ELSE
              WRITE (6,*) 'NSEPI5 ILL DEFINED IN XSTAPI '
              CALL EXIT
            ENDIF
C
130       CONTINUE
C
          NAPII(IATM)=IDSC
C  NO MODEL DEFINED
        ELSE
          NAPII(IATM)=0
        ENDIF
C
        NAPIIM(IATM)=NAPII(IATM)-1
C
        LGAPI(IATM,0,0)=0
        DO 180 IAPI=1,NAPII(IATM)
          LGAPI(IATM,0,0)=LGAPI(IATM,0,0)+LGAPI(IATM,IAPI,0)
180     CONTINUE
C
        DO 500 IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          DO 510 IAT=1,NATMI
            IA=IAT
            PATPI(IRPI,0)=PATPI(IRPI,0)+
     +                          PATPI(IRPI,IAT)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NATMI+IML
            PMLPI(IRPI,0)=PMLPI(IRPI,0)+
     +                          PMLPI(IRPI,IML)
520       CONTINUE
          DO 530 ION=1,NIONI
            IO=NSPAM+ION
            PIOPI(IRPI,0)=PIOPI(IRPI,0)+
     +                          PIOPI(IRPI,ION)
530       CONTINUE
          DO 540 IPL=1,NPLSI
            PPLPI(IRPI,0)=PPLPI(IRPI,0)+
     +                          PPLPI(IRPI,IPL)
540       CONTINUE
C
C
          P2NPI(IRPI)=PATPI(IRPI,0)+PMLPI(IRPI,0)+
     .                   PIOPI(IRPI,0)
          P2N=P2NP(IRPI,NSPAMI)
          DO 550 ISPZ=1,NSPAMI
            IF (P2N.GT.0.D0)
     .      P2NP(IRPI,ISPZ)=P2NP(IRPI,ISPZ)/P2N
550       CONTINUE
500     CONTINUE
1000  CONTINUE
C
      RETURN
C
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSTAPI: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSTAPI: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      CALL EXIT
999   CONTINUE
      WRITE (6,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EXIT
      RETURN
C
      END
C
C
      SUBROUTINE XSTAPI_PARAM

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      INTEGER :: KK, NRC
C
C   ION IMPACT COLLISIONS
C
      DO 1000 IATM=1,NATMI
C
C  NO DEFAULT MODEL
C
       IF (NRCA(IATM).EQ.0) THEN
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 130 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) GOTO 130
            NRPI=NRPI+1
C
130       CONTINUE
C
C  NO MODEL DEFINED
        ELSE
        ENDIF
C
1000  CONTINUE
C
      RETURN
C
      END
C
C
      SUBROUTINE XSTCX(RMASS,IRCX,ISP,IPL,ISCD1,ISCD2,
     .                 EBULK,ISCDE,IESTM,KK,FACTKK)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EBULK, FACTKK
      INTEGER, INTENT(IN) :: IRCX, ISP, IPL, ISCD1, ISCD2, ISCDE,
     .                       IESTM, KK
      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: ADD, ADDL, RMTEST, RMBULK, FCTKKL, ADDTL, CHRDIF,
     .            ADDT, TMASS, PMASS
      INTEGER :: ITYP1, ITYP2, ISPZ1, IERR, ISPZ2, IATM, IPLS, KREAD,
     .           IDEZ, J, NEND, MODC, NSECX4, I, IPL2, IIO2
      CHARACTER(8) :: TEXTS1, TEXTS2

      SAVE
C
C  SET NON DEFAULT CHARGE EXCHANGE COLLISION PROCESS NO. IRCX
C
      IF (IPL.LE.0.OR.IPL.GT.NPLSI) GOTO 990
      IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
      RMBULK=RMASSP(IPL)
      RMTEST=RMASS
C
C  1ST SECONDARY INDEX
      N1STX(IRCX,1)=IDEZ(ISCD1,1,3)
      N1STX(IRCX,2)=IDEZ(ISCD1,3,3)
      N1STX(IRCX,3)=0
      IF (N1STX(IRCX,1).LT.4) N1STX(IRCX,3)=1
C
      IF (N1STX(IRCX,1).EQ.1) THEN
        IF (RMBULK.NE.RMASSA(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.2) THEN
        IF (RMBULK.NE.RMASSM(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.3) THEN
        IF (RMBULK.NE.RMASSI(N1STX(IRCX,2))) GOTO 992
      ELSEIF (N1STX(IRCX,1).EQ.4) THEN
        IF (RMBULK.NE.RMASSP(N1STX(IRCX,2))) GOTO 992
      ENDIF
C  2ND SECONDARY INDEX
      N2NDX(IRCX,1)=IDEZ(ISCD2,1,3)
      N2NDX(IRCX,2)=IDEZ(ISCD2,3,3)
      N2NDX(IRCX,3)=N1STX(IRCX,3)
      IF (N2NDX(IRCX,1).LT.4) N2NDX(IRCX,3)=N2NDX(IRCX,3)+1
C
      IF (N2NDX(IRCX,1).EQ.1) THEN
        IF (RMTEST.NE.RMASSA(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.2) THEN
        IF (RMTEST.NE.RMASSM(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.3) THEN
        IF (RMTEST.NE.RMASSI(N2NDX(IRCX,2))) GOTO 992
      ELSEIF (N2NDX(IRCX,1).EQ.4) THEN
        IF (RMTEST.NE.RMASSP(N2NDX(IRCX,2))) GOTO 992
      ENDIF
C
      CHRDIF=-NCHRGP(IPL)
      IF (N1STX(IRCX,1).EQ.3) THEN
        IIO2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N1STX(IRCX,1).EQ.4) THEN
        IPL2=N1STX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.3) THEN
        IIO2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGI(IIO2)
      ENDIF
      IF (N2NDX(IRCX,1).EQ.4) THEN
        IPL2=N2NDX(IRCX,2)
        CHRDIF=CHRDIF+NCHRGP(IPL2)
      ENDIF
      IF (CHRDIF.NE.0) GOTO 990
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT; TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
      ADDT=PMASS/RMASSP(IPL)
      ADDTL=LOG(ADDT)
      NREACX(IRCX) = KK
      ADDCX(IRCX,IPL) = ADDTL
C
C CROSS SECTION (E-LAB)
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(3,1,ISP,IPL)=KK
        MODCOL(3,2,ISP,IPL)=3
        IF (FACTKK.NE.1.D0)
     .  WRITE (6,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTCX'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
        MODCOL(3,2,ISP,IPL)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN
          DO 242 J=1,NSBOX
            PLS(J)=TIINL(IPL,J)+ADDTL
242       CONTINUE
          IF (MODC.EQ.1) THEN
            CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.TRUE.,
     .                .FALSE.,.TRUE.)
            DO 245 J=1,NSBOX
              TABCX3(IRCX,J,1)=COUN(1,J)*DIIN(IPL,J)*FACTKK
245         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
            CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.FALSE.,
     .                .FALSE.,.TRUE.)
            FCTKKL=LOG(FACTKK)
            DO 246 J=1,NSBOX
              TABCX3(IRCX,J,1)=COUN(1,J)+DIINL(IPL,J)+FCTKKL
246         CONTINUE
          ENDIF
          DO 244 I=2,NEND
            DO 243 J=1,NSBOX
              TABCX3(IRCX,J,I)=COUN(I,J)
243         CONTINUE
244       CONTINUE
        ELSE
          FACREA(KK) = LOG(FACTKK)
        END IF
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
        IF (MODCOL(3,2,ISP,IPL).NE.3) GOTO 996
      ENDIF
      DEFCX(IRCX)=LOG(CVELI2*PMASS)
      EEFCX(IRCX)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
      NSECX4=IDEZ(ISCDE,4,5)
      IF (NSECX4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLCX3(IRCX,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
            ENDDO
          ELSE
            NELRCX(IRCX) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
251         CONTINUE
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(3,4,ISP,IPL)=3
      ELSEIF (NSECX4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
          DO 252 J=1,NSBOX
            EPLCX3(IRCX,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
  252     CONTINUE
          ELSE
            NELRCX(IRCX) = -3
          END IF
        ELSE
          WRITE (6,*) 'WARNING FROM SUBR. XSTCX '
          WRITE (6,*) 'MODIFIED TREATMENT OF CHARGE EXCHANGE '
          WRITE (6,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (6,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLCX3(IRCX,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
          ELSE
            NELRCX(IRCX) = -2
            EPLCX3(IRCX,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(3,4,ISP,IPL)=1
C     ELSEIF (NSECX4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        NELRCX(IRCX) = KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(3,4,ISP,IPL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPL,J)+ADDTL
253         CONTINUE
            CALL CDEF (PLS,1,NEND,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                 .FALSE.,.TRUE.)
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
                DO 254 J=1,NSBOX
                  EPLCX3(IRCX,J,1)=COUN(1,J)*DIIN(IPL,J)*ADD
254             CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
                DO 257 J=1,NSBOX
                  EPLCX3(IRCX,J,1)=COUN(1,J)+DIINL(IPL,J)+ADDL
257             CONTINUE
            ENDIF
            DO 256 I=2,NEND
              DO 255 J=1,NSBOX
                EPLCX3(IRCX,J,I)=COUN(I,J)
255           CONTINUE
256         CONTINUE
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLCX3(IRCX,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD) = ADDL
            END IF
          END IF
        ENDIF
      ELSE
        IERR=5
        GOTO 996
      ENDIF

C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTCX(IRCX,1)=IDEZ(IESTM,1,3)
      IESTCX(IRCX,2)=IDEZ(IESTM,2,3)
      IESTCX(IRCX,3)=IDEZ(IESTM,3,3)
C
      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)
      IF (IESTCX(IRCX,1).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (6,*) 'IRCX = ',IRCX
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,1)=0
      ENDIF
      IF (IESTCX(IRCX,2).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (6,*) 'IRCX = ',IRCX
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,2)=0
      ENDIF
      IF (IESTCX(IRCX,3).NE.0.AND.(ITYP1.NE.1.OR.ITYP2.NE.4)) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR EN.-BALANCE '
        WRITE (6,*) 'IRCX = ',IRCX
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTCX(IRCX,3)=0
      ENDIF
      RETURN
C
      ENTRY XSTCX_2(IRCX,IPL)
C
      CALL LEER(1)
      WRITE (6,*) 'CHARGE EXCHANGE REACTION NO. IRCX= ',IRCX
      CALL LEER(1)
      WRITE (6,*) 'CHARGE EXCHANGE WITH BULK IONS IPLS:'
      WRITE (6,*) '1ST AND 2ND NEXT GEN. SPECIES I2ND1, I2ND2:'
      ITYP1=N1STX(IRCX,1)
      ITYP2=N2NDX(IRCX,1)
      ISPZ1=N1STX(IRCX,2)
      ISPZ2=N2NDX(IRCX,2)
      IF (ITYP1.EQ.1) TEXTS1=TEXTS(NSPH+ISPZ1)
      IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
      IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
      IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
      IF (ITYP2.EQ.1) TEXTS2=TEXTS(NSPH+ISPZ2)
      IF (ITYP2.EQ.2) TEXTS2=TEXTS(NSPA+ISPZ2)
      IF (ITYP2.EQ.3) TEXTS2=TEXTS(NSPAM+ISPZ2)
      IF (ITYP2.EQ.4) TEXTS2=TEXTS(NSPAMI+ISPZ2)
      WRITE (6,*) 'IPLS= ',TEXTS(NSPAMI+IPL),'I2ND1= ',TEXTS1,
     .                    'I2ND2= ',TEXTS2
      CALL LEER(1)
      RETURN
C
990   CONTINUE
      WRITE (6,*) 'ERROR IN XSTACX: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE '
      CALL EXIT
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSTACX: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
      CALL EXIT
993   CONTINUE
      WRITE (6,*) 'ERROR IN XSTACX: EXIT CALLED  '
      WRITE (6,*) 'EBULK_ION .LE.0, BUT MONOENERGETIC DISTRIBUTION?'
      WRITE (6,*) 'CHECK ENERGY FLAG ISCDEA'
      WRITE (6,*) 'KK,IATM,IPLS,ISCDEA ',KK,IATM,IPLS,ISCDEA
      CALL EXIT
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSTACX: EXIT CALLED  '
      WRITE (6,*) 'NO CROSS SECTION AVAILABLE FOR NON DEFAULT CX'
      WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
      WRITE (6,*) 'EITHER PROVIDE CROSS SECTION OR USE DIFFERENT '
      WRITE (6,*) 'POST COLLISION SAMPLING FLAG ISCDEA'
      CALL EXIT
      END
C
C
      SUBROUTINE XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, CHRDF0, EELEC, FACTKK
      REAL(DP), INTENT(OUT) :: PLS(*)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ISCDE, IESTM, KK
      REAL(DP) :: COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL, FEHVDS1, EE, TB, FEELEI1, EN,
     .          P2N, EA, EI, ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV
      INTEGER :: MODC, KREAD, IM, IA, IERR, J, IPL, I, IP, IRAD, IO,
     .           ION, ISPZ, III, IDEZ, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML

      ITYP=IDEZ(IFRST,1,3)
      INUM=IDEZ(IFRST,2,3)
      ISPE=IDEZ(IFRST,3,3)
C
C ACCUMULATED MASS OF SECONDARIES: ACCMAS (AMU)
      ACCMAS=0.D0
      ACCMSA=0.D0
      ACCMSM=0.D0
      ACCMSI=0.D0
      ACCMSP=0.D0
      ACCINV=0.D0
      ACCINA=0.D0
      ACCINM=0.D0
      ACCINI=0.D0
      ACCINP=0.D0
C
      ICOUNT=1
 85   CONTINUE
      IF (ITYP.EQ.1) THEN
        IAT=ISPE
        IAA=NSPH+IAT
        PATDS(IREI,IAT)=PATDS(IREI,IAT)+INUM
        P2ND(IREI,IAA)=P2ND(IREI,IAA)+INUM
        ACCMAS=ACCMAS+INUM*RMASSA(IAT)
        ACCMSA=ACCMSA+INUM*RMASSA(IAT)
        ACCINV=ACCINV+INUM/RMASSA(IAT)
        ACCINA=ACCINA+INUM/RMASSA(IAT)
        EATDS(IREI,IAT,1)=RMASSA(IAT)
        EATDS(IREI,IAT,2)=1./RMASSA(IAT)
      ELSEIF (ITYP.EQ.2) THEN
        IML=ISPE
        IMM=NSPA+IML
        PMLDS(IREI,IML)=PMLDS(IREI,IML)+INUM
        P2ND(IREI,IMM)=P2ND(IREI,IMM)+INUM
        ACCMAS=ACCMAS+INUM*RMASSM(IML)
        ACCMSM=ACCMSM+INUM*RMASSM(IML)
        ACCINV=ACCINV+INUM/RMASSM(IML)
        ACCINM=ACCINM+INUM/RMASSM(IML)
        EMLDS(IREI,IML,1)=RMASSM(IML)
        EMLDS(IREI,IML,2)=1./RMASSM(IML)
      ELSEIF (ITYP.EQ.3) THEN
        IIO=ISPE
        III=NSPAM+IIO
        PIODS(IREI,IIO)=PIODS(IREI,IIO)+INUM
        P2ND(IREI,III)=P2ND(IREI,III)+INUM
        ACCMAS=ACCMAS+INUM*RMASSI(IIO)
        ACCMSI=ACCMSI+INUM*RMASSI(IIO)
        ACCINV=ACCINV+INUM/RMASSI(IIO)
        ACCINI=ACCINI+INUM/RMASSI(IIO)
        EIODS(IREI,IIO,1)=RMASSI(IIO)
        EIODS(IREI,IIO,2)=1./RMASSI(IIO)
      ELSEIF (ITYP.EQ.4) THEN
        IPL=ISPE
        PPLDS(IREI,IPL)=PPLDS(IREI,IPL)+INUM
        ACCMAS=ACCMAS+INUM*RMASSP(IPL)
        ACCMSP=ACCMSP+INUM*RMASSP(IPL)
        ACCINV=ACCINV+INUM/RMASSP(IPL)
        ACCINP=ACCINP+INUM/RMASSP(IPL)
      ENDIF
C
      IF (ISCND.NE.0.AND.ICOUNT.EQ.1) THEN
        ITYP=IDEZ(ISCND,1,3)
        INUM=IDEZ(ISCND,2,3)
        ISPE=IDEZ(ISCND,3,3)
        ICOUNT=2
        GOTO 85
      ENDIF
C
      IF (ABS(ACCMAS-RMASS).GT.1.D-20) THEN
        WRITE (6,*) 'MASS CONSERVATION VIOLATED FOR REACT. KK'
        WRITE (6,*) KK,ISP
        CALL EXIT
      ENDIF
C
      DO IAT=1,NATMI
        EATDS(IREI,IAT,1)=EATDS(IREI,IAT,1)/ACCMAS
        EATDS(IREI,IAT,2)=EATDS(IREI,IAT,2)/ACCINV
      ENDDO
      EATDS(IREI,0,1)=ACCMSA/ACCMAS
      EATDS(IREI,0,2)=ACCINA/ACCINV
      DO IML=1,NMOLI
        EMLDS(IREI,IML,1)=EMLDS(IREI,IML,1)/ACCMAS
        EMLDS(IREI,IML,2)=EMLDS(IREI,IML,2)/ACCINV
      ENDDO
      EMLDS(IREI,0,1)=ACCMSM/ACCMAS
      EMLDS(IREI,0,2)=ACCINM/ACCINV
      DO IIO=1,NIONI
        EIODS(IREI,IIO,1)=EIODS(IREI,IIO,1)/ACCMAS
        EIODS(IREI,IIO,2)=EIODS(IREI,IIO,2)/ACCINV
      ENDDO
      EIODS(IREI,0,1)=ACCMSI/ACCMAS
      EIODS(IREI,0,2)=ACCINI/ACCINV
      EPLDS(IREI,1)=ACCMSP/ACCMAS
      EPLDS(IREI,2)=ACCINP/ACCINV
C
      CHRDIF=CHRDF0
      DO 83 IIO=1,NIONI
        CHRDIF=CHRDIF+PIODS(IREI,IIO)*NCHRGI(IIO)
83    CONTINUE
      DO 84 IPL=1,NPLSI
        CHRDIF=CHRDIF+PPLDS(IREI,IPL)*NCHRGP(IPL)
84    CONTINUE
      PELDS(IREI)=PELDS(IREI)+CHRDIF
C
C
C  1.) CROSS SECTION(TE) : NOT NEEDED
C
C
C  2.) RATE COEFFICIENT (CM**3/S) * ELECTRON DENSITY (CM**-3)
C
C
C  2.A) RATE COEFFICIENT = CONST.
C     TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
      IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
        IF (NSTORDR >= NRAD) THEN
          CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .              .TRUE.)
          DO 92 J=1,NSBOX
            TABDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK
c slmod begin (sl)
            IF (IOPT1.EQ.1) TABDS1(IREI,J)=TABDS1(IREI,J)*
     .                                     TABDSM(IREI,J)
c slmod end
92        CONTINUE
        ELSE
          FACREA(KK) = LOG(FACTKK)
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,ISP,1)=1
C     ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,ISP,1)=2
      ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          CALL CDEF (TEINL,1,9,KK,COUN,NSBOX,CF,.FALSE.,.FALSE.,
     .              .TRUE.)
          DO 93 J=1,NSBOX
            TABDS1(IREI,J)=COUN(9,J)
93        CONTINUE
          DO 95 I=8,1,-1
            DO 94 J=1,NSBOX
              TABDS1(IREI,J)=TABDS1(IREI,J)*PLS(J)+COUN(I,J)
94          CONTINUE
95        CONTINUE
          FCTKKL=LOG(FACTKK)
          DO 97 J=1,NSBOX
            TB=MAX(-100._DP,TABDS1(IREI,J)+FCTKKL+DEINL(J))
            TABDS1(IREI,J)=EXP(TB)
c slmod begin (sl)
            IF (IOPT1.EQ.1) TABDS1(IREI,J)=TABDS1(IREI,J)*
     .                                     TABDSM(IREI,J)
c slmod end
97        CONTINUE
        ELSE
          FACREA(KK) = LOG(FACTKK)
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ENDIF
        MODCOL(1,2,ISP,1)=1
      ELSE
        IERR=1
        GOTO 996
      ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ENERGY RATES
C
C  4.A: ELECTRON ENERGY LOSS RATES
C
      EFLAG=IDEZ(ISCDE,5,5)
      IF (EFLAG.EQ.0) THEN
C  4.A1) ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 101 J=1,NSBOX
                  EELDS1(IREI,J)=EELEC
101             CONTINUE
              ELSE
                NELREI(IREI)=-2
                EELDS1(IREI,1)=EELEC
              END IF
              MODCOL(1,4,ISP,1)=1
      ELSEIF (EFLAG.EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 103 J=1,NSBOX
                  EELDS1(IREI,J)=-1.5*TEIN(J)
103             CONTINUE
              ELSE
                NELREI(IREI)=-3
              END IF
              MODCOL(1,4,ISP,1)=1
      ELSEIF (EFLAG.EQ.3) THEN
C  4.A2) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELEC
                MODC=IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                    DO 102 J=1,NSBOX
                      EELDS1(IREI,J)=-COUN(1,J)*DEIN(J)*FACTKK/
     .                               (TABDS1(IREI,J)+EPS60)
102                 CONTINUE
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                    FACREA(KREAD)=LOG(FACTKK)
                  ENDIF
                  MODCOL(1,4,ISP,1)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                         .FALSE.,.TRUE.)
                    DO 106 J=1,NSBOX
                      EELDS1(IREI,J)=COUN(9,J)
106                 CONTINUE
                    DO 108 I=8,1,-1
                      DO 107 J=1,NSBOX
                        EELDS1(IREI,J)=EELDS1(IREI,J)*PLS(J)+COUN(I,J)
107                   CONTINUE
108                 CONTINUE
                    FCTKKL=LOG(FACTKK)
                    DO 110 J=1,NSBOX
                      EE=MAX(-100._DP,EELDS1(IREI,J)+FCTKKL+DEINL(J))
                      EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
c slmod begin (sl)
                      IF (IOPT1.EQ.1) EELDS1(IREI,J)=EELDS1(IREI,J)*
     .                                               TABDEM(IREI,J)
c slmod end
110                 CONTINUE
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                    FACREA(KREAD)=LOG(FACTKK)
                  END IF
                  MODCOL(1,4,ISP,1)=1
                ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  4.B: HEAVY PARTICLE ENERGY GAIN RATE
C
      EFLAG=IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.B1)  RATE = CONST.*RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVDS1(IREI,J)=EHEAVY
201       CONTINUE
        ELSE
          NREAHV(IREI)=-1
          EHVDS1(IREI,1)=EHEAVY
        END IF
C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION
      ELSEIF (EFLAG.EQ.3) THEN
C  4.B3)  ENERGY RATE = EN.WEIGHTED RATE(TE)
        KREAD=EHEAVY
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                .FALSE.,.TRUE.)
            DO 202 J=1,NSBOX
              EHVDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK/
     .                       (TABDS1(IREI,J)+EPS60)
202         CONTINUE
          ELSE
            NREAHV(IREI)=KREAD
            FACREA(KREAD)=LOG(FACTKK)
          END IF
        ELSE
          WRITE (6,*) 'INVALID OPTION IN XSTEI '
          CALL EXIT
        ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEI(IREI,1)=IDEZ(IESTM,1,3)
      IESTEI(IREI,2)=IDEZ(IESTM,2,3)
      IESTEI(IREI,3)=IDEZ(IESTM,3,3)
C
      IF (IESTEI(IREI,1).NE.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (6,*) 'IREI = ',IREI
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,1)=0
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (6,*) 'IREI = ',IREI
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,2)=0
      ENDIF
      CALL LEER(1)
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_1(IREI)
C
          DO 510 IAT=1,NATMI
            IA=NSPH+IAT
            PATDS(IREI,0)=PATDS(IREI,0)+
     +                          PATDS(IREI,IAT)
            P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                          P2ND(IREI,IA)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NSPA+IML
            PMLDS(IREI,0)=PMLDS(IREI,0)+
     +                          PMLDS(IREI,IML)
            P2ND(IREI,IM)=P2ND(IREI,IM-1)+
     +                          P2ND(IREI,IM)
520       CONTINUE
          DO 530 ION=1,NIONI
            IO=NSPAM+ION
            PIODS(IREI,0)=PIODS(IREI,0)+
     +                          PIODS(IREI,ION)
            P2ND(IREI,IO)=P2ND(IREI,IO-1)+
     +                          P2ND(IREI,IO)
530       CONTINUE
          DO 540 IPL=1,NPLSI
            PPLDS(IREI,0)=PPLDS(IREI,0)+
     +                          PPLDS(IREI,IPL)
540       CONTINUE
C
          P2NDS(IREI)=PATDS(IREI,0)+PMLDS(IREI,0)+
     .                   PIODS(IREI,0)
          P2N=P2ND(IREI,NSPAMI)
          DO 550 ISPZ=NSPH+1,NSPAMI
            IF (P2N.GT.0.D0)
     .      P2ND(IREI,ISPZ)=P2ND(IREI,ISPZ)/P2N
550       CONTINUE
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_2(IREI)
C
      CALL LEER(1)
      WRITE (6,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL LEER(1)
      EI=1.D30
      EA=-1.D30
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELDS1(IREI,IRAD)
        ELSE
          EN=FEELEI1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE
      WRITE (6,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (6,*) 'ELECTRONS: PELDS, CONSTANT ENERGY: EEL'
        WRITE (6,'(1X,A8,2(1PE12.4))') 'EL      ',PELDS(IREI),EI
      ELSE
        WRITE (6,*) 'ELECTRONS: PELDS, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (6,'(1X,A8,3(1PE12.4))') 'EL      ',PELDS(IREI),EI,EA
      ENDIF
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVDS1(IREI,IRAD)
        ELSE
          EN=FEHVDS1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
876   CONTINUE
      IF (PPLDS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'BULK IONS: PPLDS '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (PPLDS(IREI,IPL).NE.0.D0)
     .      WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IP),PPLDS(IREI,IPL)
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EPLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EPLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2),
     .                                 ' * EHEAVY '
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (6,*) 'BULK IONS: NONE'
      ENDIF
      CALL LEER(1)
C
      WRITE (6,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NDS(IREI).EQ.0.D0) THEN
        WRITE (6,*) 'NONE'
        CALL LEER(1)
        RETURN
      ENDIF
C
      IF (PATDS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'ATOMS    : PATDS '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IA),PATDS(IREI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EATDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EATDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLDS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'MOLECULES: PMLDS '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLDS(IREI,IML).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IM),PMLDS(IREI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EMLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EMLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIODS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'TEST IONS: PIODS '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIODS(IREI,IIO).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IO),PIODS(IREI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EIODS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EIODS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      RETURN
C
C
C-----------------------------------------------------------------------
C
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (6,*) IREI,KK
      CALL EXIT
997   CONTINUE
      WRITE (6,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (6,*) IREI
      CALL EXIT
      END
C
C
      SUBROUTINE XSTEL(IREL,ISP,IPL,
     .                 EBULK,ISCDE,IESTM,
     .                 KK,FACTKK)
C
C  RETURNS:
C    MODCOL(5,...)
C    TABEL3(IREL,NCELL,...)
C    EPLEL3(IREL,NCELL,...)
C    DEFEL(IREL)
C    EEFEL(IREL)
C    IESTEL(IREL,...)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: EBULK, FACTKK
      INTEGER, INTENT(IN) :: IREL, ISP, IPL, ISCDE, IESTM, KK
      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: FCTKKL, ADD, ADDL, ADDT, ADDTL, PMASS, TMASS
      INTEGER :: I, NSEEL4, NEND, J, KREAD, MODC, IDEZ

      SAVE
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT; TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
      ADDT=PMASS/RMASSP(IPL)
      ADDTL=LOG(ADDT)
      ADDEL(IREL,IPL) = ADDTL
      NREAEL(IREL) = KK
C
C POTENTIAL
      IF (IDEZ(MODCLF(KK),1,5).EQ.1) THEN
        MODCOL(5,0,ISP,IPL)=KK
      ENDIF
C
C CROSS SECTION (E-LAB), IN FUNCTION CROSS, K=KK
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(5,1,ISP,IPL)=KK
        MODCOL(5,2,ISP,IPL)=3
        IF (FACTKK.NE.1.D0)
     .    WRITE (6,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTEL'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
        MODCOL(5,2,ISP,IPL)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN
          DO 242 J=1,NSBOX
            PLS(J)=TIINL(IPL,J)+ADDTL
242       CONTINUE
          IF (MODC.EQ.1) THEN
            CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.TRUE.,
     .                .FALSE.,.TRUE.)
            DO 245 J=1,NSBOX
              TABEL3(IREL,J,1)=COUN(1,J)*DIIN(IPL,J)*FACTKK
245         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
            CALL CDEF (PLS,1,NEND,KK,COUN,NSBOX,CF,.FALSE.,
     .                .FALSE.,.TRUE.)
            FCTKKL=LOG(FACTKK)
            DO 246 J=1,NSBOX
              TABEL3(IREL,J,1)=COUN(1,J)+DIINL(IPL,J)+FCTKKL
246         CONTINUE
          ENDIF
          DO 244 I=2,NEND
            DO 243 J=1,NSBOX
              TABEL3(IREL,J,I)=COUN(I,J)
243         CONTINUE
244       CONTINUE
        ELSE
          FACREA(KK) = LOG(FACTKK)
        END IF
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
        IF (MODCOL(5,2,ISP,IPL).NE.3) GOTO 993
      ENDIF
C
      DEFEL(IREL)=LOG(CVELI2*PMASS)
      EEFEL(IREL)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4. BULK ION ENERGY LOSS RATE
C
      NSEEL4=IDEZ(ISCDE,4,5)
      IF (NSEEL4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK ION = CONST.*RATECOEFF.
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
            ENDDO
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK
251         CONTINUE
          ELSE
            NELREL(IREL) = -1
            EPLEL3(IREL,1,1)=EBULK
          END IF
        END IF
        MODCOL(5,4,ISP,IPL)=1
      ELSEIF (NSEEL4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
252         CONTINUE
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          WRITE (6,*) 'WARNING FROM SUBR. XSTEL '
          WRITE (6,*) 'MODIFIED TREATMENT OF ELASTIC COLLISIONS '
          WRITE (6,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (6,*) 'RATHER THEN WITH T = TIIN '
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
          ELSE
            NELREL(IREL) = -2
            EPLEL3(IREL,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(5,4,ISP,IPL)=1
C     ELSEIF (NSEEL4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEEL4.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        NELREL(IREL)=KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(5,4,ISP,IPL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPL,J)+ADDTL
253         CONTINUE
            CALL CDEF (PLS,1,NEND,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                .FALSE.,.TRUE.)
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              DO 257 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF
            DO 256 I=2,NEND
              DO 255 J=1,NSBOX
                EPLEL3(IREL,J,I)=COUN(I,J)
255           CONTINUE
256         CONTINUE
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLEL3(IREL,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD) = ADDL
            END IF
          END IF
        ENDIF
      ELSE
        GOTO 993
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEL(IREL,1)=IDEZ(IESTM,1,3)
      IESTEL(IREL,2)=IDEZ(IESTM,2,3)
      IESTEL(IREL,3)=IDEZ(IESTM,3,3)
C
      IF (IESTEL(IREL,2).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: TR.L.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (6,*) 'IREL = ',IREL
        WRITE (6,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,2)=1
      ENDIF
      IF (IESTEL(IREL,3).EQ.0.AND.NPBGKP(IPL,1).EQ.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: TR.L.EST NOT AVAILABLE FOR EN.-BALANCE '
        WRITE (6,*) 'IREL = ',IREL
        WRITE (6,*) 'AUTOMATICALLY RESET TO COLLISION ESTIMATOR '
        IESTEL(IREL,3)=1
      ENDIF
      RETURN
C
      ENTRY XSTEL_2(IREL,IPL)
C
      CALL LEER(1)
      WRITE (6,*) 'ELASTIC COLLISION NO. IREL= ',IREL
      CALL LEER(1)
      WRITE (6,*) 'ELASTIC COLLISION WITH BULK IONS IPLS:'
      WRITE (6,*) 'IPLS= ',TEXTS(NSPAMI+IPL)
C
      CALL LEER(1)
      RETURN
C
993   CONTINUE
      WRITE ( 6,*) 'ERROR IN XSTEL, SPECIES ISP: '
      WRITE (6,*) ISP,IREL
      CALL EXIT
      END
