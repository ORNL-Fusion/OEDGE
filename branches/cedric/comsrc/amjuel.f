c     -*-Fortran-*-
c
c ======================================================================
c
c subroutine: DumpRates
c
      SUBROUTINE DumpRates(iopt)
      IMPLICIT none

      INTEGER i1,iopt,nrates(3)
      REAL    dat(20),ne,nel,te,GetEAD

      DATA nrates /4, 5, 3/

      CALL LoadEIRENEAtomicData

      WRITE(6,*) 'AMJUEL RATES:',iopt

      DO nel = 17.0, 21.0, 0.25
 
        ne = 10.0**nel
 
        DO te = 2.0, 5.0, 0.2

          IF     (iopt.EQ.1) THEN
c...        H rates:
            dat(1)=GetEAD(te,ne, 1,'H.4 ') * 1.0E-06          
            dat(2)=GetEAD(te,ne, 2,'H.4 ') * 1.0E-06          
            dat(3)=GetEAD(te,ne, 9,'H.1 ')         
            dat(4)=GetEAD(te,ne,10,'H.3 ') * 1.0E-06          
          ELSEIF (iopt.EQ.2) THEN
c...        H2 rates:
            dat(1)=GetEAD(te,ne,11,'H.1 ')         
            dat(2)=GetEAD(te,ne,12,'H.3 ') * 1.0E-06          
            dat(3)=GetEAD(te,ne,13,'H.4 ') * 1.0E-06          
            dat(4)=GetEAD(te,ne,14,'H.4 ') * 1.0E-06          
            dat(5)=GetEAD(te,ne,18,'H.3 ') * 1.0E-06          
          ELSEIF (iopt.EQ.3) THEN
c...        H2+ rates:
            dat(1)=GetEAD(te,ne,15,'H.4 ') * 1.0E-06          
            dat(2)=GetEAD(te,ne,16,'H.4 ') * 1.0E-06          
            dat(3)=GetEAD(te,ne,17,'H.3 ') * 1.0E-06          
          ELSE
c...        Error:
            WRITE(6,*) 'DumpRates: Unknown option',iopt
            RETURN
          ENDIF

          WRITE(6,'(F8.4,1P,E12.4,10(E10.2:))') 
     .      te,ne,(dat(i1),i1=1,nrates(iopt))
        ENDDO
      ENDDO


      RETURN
99    STOP
      END
c
c ======================================================================
c
c subroutine: CalcRadiatedPower
c
c input:
c   mode  - set = 1
c
c output:
c   array - hydrogenic radiation power loss (so positive values indicate
c           power loss)
c
c
      SUBROUTINE CalcRadiatedPower(array,mode)
      IMPLICIT none

      INCLUDE 'params'
      INCLUDE 'cgeom'
      INCLUDE 'comtor'
      INCLUDE 'pindata'
      INCLUDE 'slcom'

c...  Input:
      INTEGER mode
c...  Output:
      REAL    array(MAXNKS,MAXNRS)

      REAL GetEAD

      INTEGER ik,ir
      REAL    rion,rpow,prad

c...  Only needs to be called once -- move to start of OUT:
      CALL LoadEIRENEAtomicData

      DO ir = 1, nrs
        DO ik = 1, nks(ir)

        IF (mode.EQ.1) THEN
          rion = GetEAD(ktebs(ik,ir),knbs(ik,ir),1,'H.4 ') * 1.0E-06
          rpow = GetEAD(ktebs(ik,ir),knbs(ik,ir),5,'H.10') * 1.0E-06*ECH
        ELSE
          CALL ER('CalcRadiatedPower','Invalid MODE',*99)
        ENDIF

        prad = (rpow - 13.6 * ECH * rion) * knbs(ik,ir) * pinatom(ik,ir)

        array(ik,ir) = prad

c        WRITE(6,90) 'PRAD : ',ik,ir,ktebs(ik,ir),knbs(ik,ir),rion,rpow,
c     .              prad,-rpow*knbs(ik,ir)*pinatom(ik,ir),
c     .              -rpow*knbs(ik,ir)*pinatom(ik,ir)*kvols(ik,ir),
c     .              pinqe(ik,ir)*kvols(ik,ir),
c     .              kvols(ik,ir)*1.00E+6
c90      FORMAT(A,2I4,F6.1,1P,E10.2,2X,3E10.2,2X,2E10.2,2X,2E12.4,0P)
        ENDDO
      ENDDO

      RETURN
99    STOP
      END
c
c ======================================================================
c
c function: GetEAD
c
c
c
c
      REAL FUNCTION GetEAD(te,ne,ir,h123)
      IMPLICIT none

      INCLUDE 'PARMMOD'
      INCLUDE 'params'
      INCLUDE 'slcom'

      INTEGER   ir
      REAL      te,ne
      CHARACTER h123*4

      INTEGER i,j
      REAL*8  cf(9,0:9),lnte(1),lnne(1),coun(0:9,1),rdum1,
     .        dpls3,def,tef,dej,tei,DP31(0:8,0:8),
     .        ratio7,RHMH2(0:8)

      LOGICAL loadEAD
      DATA loadEAD /.TRUE./

      SAVE

c      IF (.NOT.loadEAD) CALL ER('GetEAD','EIRENE atomic data not'//
c     .                                   'loaded',*99)

      IF (loadEAD) THEN
        CALL LoadEIRENEAtomicData
        loadEAD = .FALSE.
      ENDIF

c...check to make sure ir is legit...

      lnte(1) = LOG(DBLE(te))
      lnne(1) = LOG(DBLE(ne * 1.0E-06)) - LOG(1.0D+08)

      IF     (h123.EQ.'H.1 ') THEN
c...    Fudgesicle - crude interpolation:
        IF     (ir.EQ.9 ) THEN
c...      p + D -> D + p:
          lnte(1) = DLOG10(DBLE(1.5D0*te))
          coun(1,1) = (lnte(1) - 0.0D0) / 2.0D0 * (3.0D-15 - 6.0D-15) +
     .                6.0D-15
        ELSEIF (ir.EQ.11) THEN
c...      p + D2 -> p + D2
          lnte(1) = DLOG10(DBLE(1.5D0*te))
          coun(1,1) = (lnte(1) - 0.0D0) / 2.0D0 * (1.5D-15 - 6.0D-15) +
     .                6.0D-15
        ELSE
          CALL ER('GetEAD','Sorry, no guess for reaction',*99)
        ENDIF 
      ELSEIF (h123.EQ.'H.3 ') THEN
        CALL CDEF (lnte,1,9,ir,COUN,1,CF,.FALSE.,.FALSE.,.TRUE.)
        coun(1,1) = EXP(coun(1,1))
c        coun(1,1) = coun(1,1) * DBLE(ne) * 1.0D-06
      ELSEIF (h123.EQ.'H.4 ') THEN
        CALL CDEFN(lnte,lnne,ir,coun,1,cf,.TRUE.,.FALSE.,.TRUE.)
      ELSEIF (h123.EQ.'H.10') THEN
        CALL CDEF (lnte,1,9,ir,COUN,1,CF,.FALSE.,.FALSE.,.TRUE.)
        rdum1=COUN(9,1)
        DO I=8,1,-1
          rdum1=rdum1*lnne(1)+COUN(I,1)
        ENDDO
        rdum1=MAX(-100.D0,rdum1)
        rdum1=EXP(rdum1)

c          WRITE(6,'(A,1P,E12.4,0P,2F8.4)')
c     .      'MARK: EELDS1 A= ',
c     .      rdum1,lnte(1),lnne(1)

        coun(1,1) = rdum1

      ELSEIF (h123.EQ.'H.11') THEN
        DO I=1,9
          RHMH2(I-1)=CREAC(I,1,IR)
        ENDDO
        RATIO7=0.0
        TEF=LOG(DBLE(TE))
        DO 160 I=0,8
          TEI=TEF**I
          RATIO7=RATIO7+RHMH2(I)*TEI
160     CONTINUE
        RATIO7=EXP(RATIO7)
        coun(1,1) = RATIO7

      ELSEIF (h123.EQ.'H.12') THEN

        IF (ir.EQ.20.OR.ir.EQ.21.OR.ir.EQ.25.OR.ir.EQ.26.OR.
     .      ir.EQ.22) THEN
c...      Recombination emission data, also MAR H2+/H2 ratio:
          DO J=1,9
            DO I=1,9
              DP31(J-1,I-1)=CREAC(J,I,IR)
            ENDDO
          ENDDO
          DEF=LOG(DBLE(NE*1.0E-6)*1.D-8)
          TEF=LOG(DBLE(TE))
          DPLS3=0.
          DO J=0,8
            DEJ=DEF**J
            DO I=0,8
              TEI=TEF**I
              DPLS3=DPLS3+DP31(I,J)*TEI*DEJ
            ENDDO
          ENDDO
          DPLS3=EXP(DPLS3)
          coun(1,1) = dpls3
        ELSE
          CALL ER('GetEAD','Invalid reaction',*99)
        ENDIF

      ELSE
        CALL ER('GetEAD','Unsupported option',*99)
      ENDIF

      GetEAD = SNGL(coun(1,1))

      RETURN
99    STOP
      END


c
c ======================================================================
c
c subroutine: LoadEIRENEAtomicData
c
c note: MASST,MASSP are not used here, althought they are specified in
c       the EIRENE in put file
c
c       DELPOT is also not used here, but is required for
c       AMJUEL H.102.1.8    RC  (see EIRENE input file)
c
c       SUBROUTINE SLREAC (IR,FILNAM,H123,REAC,CRC)
c       FILNAM: read a&m data from file filnam, e.g. AMJUEL, HYDHEL, METHAN, CONST
c       IR    : store data on eirene array CREAC(...,...,IR)
c       H123  : identifyer for reaction in filnam, e.g. H.1, H.2, H.3, ....
c       REAC  : number of reaction in filnam, e.g. 2.2.5
c       CRC   : type of data process e.g. EI, CX, OT, etc
c
c       1 AMJUEL H.4 2.1.5    EI   0  1
c
c       CALL SLREAC(ir,'<filnam>','<H123>','<reac>','<crc>')
c
c
      SUBROUTINE LoadEIRENEAtomicData
      IMPLICIT none

      INCLUDE 'PARMMOD'

      INTEGER i1

c      loadEAD = .TRUE.

      CALL DZero(creac(1,0,-10),9*(1+9)*(11+NREAC))
      CALL DZero(fparm,NREAC*6*2)
      CALL IZero(ifexmx,NREAC*2)
      CALL IZero(ifexmn,NREAC*2)

c...  Check that the database files are around:
      OPEN  (UNIT=29,FILE='AMJUEL', ERR=95)
      CLOSE (UNIT=29)
      OPEN  (UNIT=29,FILE='METHANE',ERR=96)
      CLOSE (UNIT=29)
      OPEN  (UNIT=29,FILE='HYDHEL' ,ERR=97)
      CLOSE (UNIT=29)
      OPEN  (UNIT=29,FILE='H2VIBR' ,ERR=98)
      CLOSE (UNIT=29)

c...Set bounds on data (when to activate extrapolation):
      DO i1 = 1, NREAC
        RCMN(i1,1)=-20.
        RCMX(i1,1)= 20.
        RCMN(i1,2)=-20.
        RCMX(i1,2)= 20.
      ENDDO

c...  H + e -> H+ + 2e
      CALL SLREAC( 1,'AMJUEL  ','H.4 ','2.1.5    ','EI ')
      CALL SLREAC(28,'H2VIBR  ','H.4 ','2.1.5a   ','EI ')

      CALL SLREAC( 2,'AMJUEL  ','H.4 ','2.1.5o   ','EI ')
      CALL SLREAC( 3,'AMJUEL  ','H.4 ','2.1.8    ','RC ')
      CALL SLREAC( 4,'AMJUEL  ','H.4 ','2.1.8o   ','RC ')
      CALL SLREAC(27,'H2VIBR  ','H.4 ','2.1.8a   ','RC ')



      CALL SLREAC( 5,'AMJUEL  ','H.10','2.1.5    ','EI ')
      CALL SLREAC( 6,'AMJUEL  ','H.10','2.1.5o   ','EI ')
      CALL SLREAC( 7,'AMJUEL  ','H.10','2.1.8    ','RC ')
      CALL SLREAC( 8,'AMJUEL  ','H.10','2.1.8o   ','RC ')

c...  p + H -> H + p
      CALL SLREAC( 9,'AMJUEL  ','H.1 ','3.1.8    ','CX ')
      CALL SLREAC(10,'AMJUEL  ','H.3 ','3.1.8    ','CX ')
c...  p + D2 -> p + D2
      CALL SLREAC(11,'AMJUEL  ','H.1 ','0.3      ','EL ')
      CALL SLREAC(12,'AMJUEL  ','H.3 ','0.3      ','EL ')
c...  e + H2 + e -> e + H2+
      CALL SLREAC(19,'AMJUEL  ','H.4 ','2.2.9    ','DS ')
c...  e + H2 -> e + H + H+
      CALL SLREAC(13,'AMJUEL  ','H.4 ','2.2.5    ','DS ')
c...  e + H2 -> e + H + H+
      CALL SLREAC(14,'AMJUEL  ','H.4 ','2.2.10   ','DS ')
c...  e + H2+ -> e + H + H+
      CALL SLREAC(15,'AMJUEL  ','H.4 ','2.2.12   ','DS ')
c...  e + H2+ -> e + H+ + H+
      CALL SLREAC(16,'AMJUEL  ','H.4 ','2.2.11   ','DS ')
c...  e + H2+ -> e + H + H
      CALL SLREAC(17,'AMJUEL  ','H.4 ','2.2.14   ','DS ')
c...  p + H2(v) -> H + H2+
      CALL SLREAC(18,'AMJUEL  ','H.3 ','3.2.3    ','CX ')



c...  REC EMISSION RATES - p(5)/n+:
      CALL SLREAC(20,'AMJUEL  ','H.12','2.1.8d   ','OT ')
c...  REC EMISSION RATES - p(3)/n+:
      CALL SLREAC(21,'AMJUEL  ','H.12','2.1.8a   ','OT ')
c...  REC EMISSION RATES - p(2)/n+:
      CALL SLREAC(25,'AMJUEL  ','H.12','2.1.8b   ','OT ')
c...  REC EMISSION RATES - p(4)/n+:
      CALL SLREAC(26,'AMJUEL  ','H.12','2.1.8c   ','OT ')


c...  D2+/D2 poulation ratio: 
      CALL SLREAC(22,'AMJUEL  ','H.11','2.0a     ','OT ')
c      CALL SLREAC(22,'AMJUEL  ','H.12','2.0c     ','OT ')

c...  Radiative and 3-body recombination rates:
      CALL SLREAC(23,'AMJUEL  ','H.4 ','2.1.8a   ','RC ')
      CALL SLREAC(24,'AMJUEL  ','H.4 ','2.1.8b   ','RC ')


      RETURN
95    CALL ER('LoadEIRENEAtomicData','AMJUEL not found' ,*99)
96    CALL ER('LoadEIRENEAtomicData','METHANE not found',*99)
97    CALL ER('LoadEIRENEAtomicData','HYDHEL not found' ,*99)
98    CALL ER('LoadEIRENEAtomicData','H2VIBR not ound' ,*99)
99    STOP
      END
c
c ======================================================================
c     
c slmod begin
      DOUBLE PRECISION FUNCTION EXTRAP(ELAB,IFLAG,FP1,FP2,FP3)
      IMPLICIT none
      INTEGER          IFLAG
      DOUBLE PRECISION ELAB,FP1,FP2,FP3

      DOUBLE PRECISION EL,X
c
c      FUNCTION EXTRAP(ELAB,IFLAG,FP1,FP2,FP3)
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c slmod end
C  NOTE:
C  INPUT:  ELAB IS LOG OF RELATIVE ENERGY, OR LOG OF TEMP
C  OUTPUT: EXTRAP IS NOT LOG, BUT THE TRUE VALUE
C
C  FUNCTION FOR EXTRAPOLATING SINGLE PARAMETER FITS BEYOND THEIR
C  RANGE OF VALIDITY
C  TYPE  IFLAG=1--4: JANEV ET AL. , SPRINGER, 1987, P13
C  TYPE  IFLAG=5  BACHMANN ET AL., IPP-REPORT, .....ELASTIC
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
C
      SUBROUTINE CDEF(AL,JI,JE,K,COU,NTE,CF,LEXP,LTEST,LSUM)
c slmod begin
      IMPLICIT none
      INTEGER          JI,JE,K,J,II,ICELL,NTE
      DOUBLE PRECISION AL,COU,CFP,CTEST,CF,
     .                 CF1,CF2,CF3,CF4,CF5,CF6,CF7,CF8,CF9,CF10,CF11
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c slmod end
      INCLUDE 'PARMMOD'
      INCLUDE 'COMXS'
      INCLUDE 'CCONA'
      LOGICAL LEXP,LTEST,LSUM
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
      DIMENSION AL(*),COU(0:9,*),CF(9,0:9)
      DIMENSION CFP(9,11)
      DIMENSION CF1(9),CF2(9),CF3(9),
     .          CF4(9),CF5(9),CF6(9),CF7(9),
     .          CF8(9),CF9(9),CF10(9),CF11(9)
C
      EQUIVALENCE (CFP(1,1),CF1(1)),(CFP(1,2),CF2(1)),
     .            (CFP(1,3),CF3(1)),(CFP(1,4),CF4(1)),
     .            (CFP(1,5),CF5(1)),(CFP(1,6),CF6(1)),
     .            (CFP(1,7),CF7(1)),(CFP(1,8),CF8(1)),
     .            (CFP(1,9),CF9(1)),(CFP(1,10),CF10(1)),
     .            (CFP(1,11),CF11(1))
C
C  K>0:
C  DATA FROM ARRAY CREAC(9,0:9,K)
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
          COU(JE,ICELL)=EXP(MAX(-100.D0,COU(JE,ICELL)))
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
      SUBROUTINE CDEFN(AL,PL,K,COU,NTE,CF,LEXP,LTEST,LSUM)
c slmod begin
      IMPLICIT none
      INTEGER          K,NTE,J,II,CTEST,ICELL,JJ,I,IFEX,KK
      DOUBLE PRECISION AL,PL,COU,CF,S01,S02,DS12,EXPO1,EXPO2,CCXM1,
     .                 CCXM2,FPAR1,FPAR2,FPAR3,DUMMP,EXTRAP
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c slmod end
      INCLUDE 'PARMMOD'
      INCLUDE 'COMXS'
      INCLUDE 'CCONA'
      LOGICAL LEXP,LTEST,LSUM
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
      DIMENSION AL(*),PL(*),COU(0:9,*),CF(9,0:9), DUMMP(9)
C
C  K>0:
C  DATA FROM ARRAY CREAC(9,0:9,K)
C
      IF (K.GT.0) THEN
C  DATA FROM A&M DATA FILES
        DO 11 J=1,9
          DO 12 II=1,9
            CF(II,J)=CREAC(II,J,K)
c            IF (k.EQ.3)
c     .      WRITE(6,*) 'MARK: CF= ',j,ii,creac(ii,j,k)
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
            if (.not.lexp) cou(1,icell)=dlog(cou(1,icell))
c            IF (k.EQ.3)
c     .      WRITE(6,*) 'MARK: EXT=',al(icell),pl(icell),cou(1,icell),
c     .        exp(al(icell)),exp(pl(icell)),rcmn(k,2)

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
c            IF (k.EQ.3)
c     .      WRITE(6,*) 'MARK: CFV=',al(icell),pl(icell),cou(1,icell),
c     .        exp(al(icell)),exp(pl(icell))
C
C           DO 22 J=1,9
C             JJ=J-1
C             DO 22 I=1,9
C               II=I-1
C               COU(1,ICELL)=COU(1,ICELL)+AL(ICELL)**II*
C    .                                    PL(ICELL)**JJ*CF(I,J)
C22          CONTINUE
C
            if (lexp) COU(1,ICELL)=EXP(MAX(-100.D0,COU(1,ICELL)))
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
c
      SUBROUTINE SLREAC (IR,FILNAM,H123,REAC,CRC)
c
C  input
C    FILNAM: read a&m data from file filnam, e.g. AMJUEL, HYDHEL, METHAN, CONST
c    IR    : store data on eirene array CREAC(...,...,IR)
c    H123  : identifyer for reaction in filnam, e.g. H.1, H.2, H.3, ....
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    CRC   : type of data process e.g. EI, CX, OT, etc
C  output
c    ISWR  : eirene flag for type of process  (1,2,...7)
c    CREAC : eirene storage array for a&m data CREAC(9,0:9,IR)
c    MODCOL: see below
c
c slmod begin
      IMPLICIT none
      INTEGER          IUNIN,IR,I0,ISW,IREAC,IC,J,IND,I,I0P1,I1,I2,IH,K
      DOUBLE PRECISION CONST

      DATA IUNIN /0/
c
c      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
c slmod end
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C  INPUT (PARAMETERLIST):
C    IR: REACTION NUMBER IN CREAC
C    FILNAM: HYDHEL, METHANE, AMJUEL OR CONST.
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (4 DIGITS NMLK)
C    FIRST DEZIMAL  K           =1  CROSS SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    SECOND DEZIMAL L           =1  <SIGMA V> FOR ONE
C                                   PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  <SIGMA V> FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  <SIGMA V> FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   L           =0  ELSE
C    THIRD  DEZIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FOURTH DEZIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
C                                   PROJECTILE ENERGY OR ELECTRON
C                                   DENSITY) AVAILABLE
C                                   (ON CREAC(..,1,IR))
C                               =2  DELTA E FOR
C                                   9 PROJECTILE ENERGIES AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                               =3  DELTA E FOR
C                                   9 ELECTRON DENSITIES  AVAILABLE
C                                   (ON CREAC(..,J,IR),J=1,9)
C                   N           =0  ELSE
C
      INCLUDE 'PARMMOD'
      INCLUDE 'COMPRT'
      INCLUDE 'COMXS'
      CHARACTER ZEILE*80
      CHARACTER FILNAM*8,H123*4,REAC*9,CRC*3
      CHARACTER AMJUEL*6,HYDHEL*6,METHANE*7,H2VIBR*6
      CHARACTER*2 CHR
      CHARACTER*3 CHRL,CHRR
      LOGICAL LCONST,LGEMIN,LGEMAX
C
      LGEMIN=.FALSE.
      LGEMAX=.FALSE.
      ISWR(IR)=0
      CONST=0.
      CHR='l0'
      I0=0
C
      AMJUEL='AMJUEL'
      HYDHEL='HYDHEL'
      METHANE='METHANE'
      H2VIBR='H2VIBR'
C
      IF (INDEX(CRC,'EI').NE.0.OR.
     .    INDEX(CRC,'DS').NE.0) ISWR(IR)=1
      IF (INDEX(CRC,'CX').NE.0) ISWR(IR)=3
      IF (INDEX(CRC,'II').NE.0.OR.
     .    INDEX(CRC,'PI').NE.0) ISWR(IR)=4
      IF (INDEX(CRC,'EL').NE.0) ISWR(IR)=5
      IF (INDEX(CRC,'RC').NE.0) ISWR(IR)=6
      IF (INDEX(CRC,'OT').NE.0) ISWR(IR)=7
C
      IF (INDEX(FILNAM,'AMJUEL').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='AMJUEL')
      ELSEIF (INDEX(FILNAM,'METHAN').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='METHANE')
      ELSEIF (INDEX(FILNAM,'HYDHEL').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='HYDHEL')
      ELSEIF (INDEX(FILNAM,'H2VIBR').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='H2VIBR')
      ELSE
        LCONST=.TRUE.
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF
C
      IREAC=INDEX(REAC,' ')-1
      IF (IREAC.LT.0) IREAC=LEN(REAC)

C  H.1
      IF (ISW.EQ.1) THEN
        CHR='a0'
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+1
C  H.2
      ELSEIF (ISW.EQ.2) THEN
        CHR='b0'
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10
C  H.3
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+20
        I0=1
C  H.4
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+30
        I0=1
C  H.5
      ELSEIF (ISW.EQ.5) THEN
        CHR='e0'
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
C  H.6
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+200
        I0=1
C  H.7
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+300
        I0=1
C  H.8
      ELSEIF (ISW.EQ.8) THEN
        CHR='h0'
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
C  H.9
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+2000
        I0=1
C  H.10
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+3000
        I0=1
C  H.11
      ELSEIF (ISW.EQ.11) THEN
        CHR='k0'
        CHRL='kl0'
        CHRR='kr0'
        I0=1
C  H.12
      ELSEIF (ISW.EQ.12) THEN
        I0=1
      ENDIF
C
      IF (LCONST) THEN
C
C  READ 9 FIT COEFFICIENTS FROM INPUT FILE
        READ (IUNIN,6664) (CREAC(IC,I0,IR),IC=1,9)
        RETURN
C
C  READ FROM DATA FILE
C
      ELSEIF (.NOT.LCONST) THEN
100     READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'##BEGIN DATA HERE##').EQ.0) GOTO 100

1       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,H123).EQ.0) GOTO 1
C
2       READ (29,'(A80)',END=990) ZEILE
        IF (INDEX(ZEILE,'H.').NE.0) GOTO 990
        IF (INDEX(ZEILE,'Reaction ').EQ.0.or.
     .      INDEX(ZEILE,REAC(1:ireac)).EQ.0) GOTO 2
      ENDIF
C
C  SINGLE PARAM. FIT, ISW=1,2,5,8,11
      IF (ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.ISW.EQ.11)
     .THEN
        IF (.NOT.LCONST) THEN
3         READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,CHR).EQ.0) GOTO 3
          DO 9 J=0,2
            IND=0
            DO 4 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
4           CONTINUE
            READ (29,'(A80)',END=990) ZEILE
9         CONTINUE
C
C  READ ASYMPTOTICS, IF AVAILABLE
C  I0P1=1 FOR CROSS SECTION
C  I0P1=2 FOR (WEIGHTED) RATE
          I0P1=I0+1
          IF (INDEX(ZEILE,CHRL).NE.0.AND.IFEXMN(IR,I0P1).EQ.0) THEN
c            WRITE(0,*) 'READING FPARM',ir
            IND=0
            DO 5 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
5           CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0.AND.IFEXMX(IR,I0P1).EQ.0) THEN
c            WRITE(0,*) 'READING FPARM',ir
            IND=0
            DO 7 I=4,6
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
7           CONTINUE
            LGEMAX=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
c
          if (lgemin.and.ifexmn(ir,I0P1).eq.0) then
c slmod begin - tmp
c            WRITE(0,*) 'ASSIGNING RCMN?',ir
c slmod end
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmn(IR,I0P1)
            rcmn(ir,I0P1)=log(rcmn(ir,I0P1))
            ifexmn(ir,I0P1)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
          if (lgemax.and.ifexmx(ir,I0P1).eq.0) then
c slmod begin - tmp
c            WRITE(0,*) 'ASSIGNING RCMX?',ir
c slmod end
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmx(IR,I0P1)
            rcmx(ir,I0P1)=log(rcmx(ir,I0P1))
            ifexmx(ir,I0P1)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
C
        ENDIF
C
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
        DO 11 J=0,2
16        READ (29,'(A80)',END=990) ZEILE
          IF (INDEX(ZEILE,'Index').EQ.0) GOTO 16
          READ (29,'(1X)')
          DO 17 I=1,9
c slmod begin - pgi
c Not sure what the problem was here, but unless the following change was
c made, the module could not be compiled using -O (-g worked fine) with the
c PGI HPF compiler.
            I1=J*3+1
            I2=J*3+3
            READ (29,*) IH,(CREAC(I,K,IR),K=I1,I2)
c
c            READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
c slmod end
17        CONTINUE
11      CONTINUE
C   NO ASYMPTOTICS AVAILABLE YET
C
      ENDIF
C
      CLOSE (UNIT=29)
C
      RETURN
C
990   WRITE (6,*) ' NO DATA FOUND FOR REACTION ',H123,' ',REAC,
     .            ' IN DATA SET ',FILNAM
      WRITE (6,*) ' IR,MODCLF(IR) ',IR,MODCLF(IR)
      CLOSE (UNIT=29)
      CALL EXIT
991   WRITE (6,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (6,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CALL EXIT
6664  FORMAT (6E12.4)
      END











