C
C
      SUBROUTINE RDTRIM
C
C  THIS SUBROUTINE READS SELECTIVELY SOME
C  REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C
      USE PRECISION
      USE PARMMOD
      USE CREFMOD
      USE CREF
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, IUN, IFILE, I, IWWW
C
C
      INE=12
      INW=7
      INR=5
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (6,*) 'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT
      ENDIF
C
      DO 7 IFILE=1,NFLR
        IUN=20
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED')
C
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        READ (IUN,*)
        DO 2 I1=1,INE
          DO 3 I2=1,INW
            READ (IUN,*)
            READ (IUN,*)
            READ (IUN,*) TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
     .                   enar(i1),wiar(i2),HFTR0(I1,I2,IFILE)
C  FIND NEAREST INTEGER FOR NUCLEAR MASS NUMBER
            IWWW=NINT(WM(IFILE))
            WM(IFILE)=IWWW
            READ (IUN,*)
            READ (IUN,*) (HFTR1(I1,I2,I3,IFILE),I3=1,INR)
            READ (IUN,*)
            DO 5 I3=1,INR
              READ (IUN,*) (HFTR2(I1,I2,I3,I4,IFILE),I4=1,INR)
5           CONTINUE
            READ (IUN,*)
            DO 6 I3=1,INR
            DO 6 I4=1,INR
              READ (IUN,*) (HFTR3(I1,I2,I3,I4,I5,IFILE),I5=1,INR)
6           CONTINUE
3         CONTINUE
2       CONTINUE
        CLOSE (UNIT=IUN)
7     CONTINUE
C
      INEM=INE-1
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      PID180=ATAN(1.)/45.
      DO 12 I=1,INW
12      WIAR(I)=COS(WIAR(I)*PID180)
      INWM=INW-1
      DO 13 I=1,INWM
13      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      RETURN
      END
C
C
      SUBROUTINE REFDAT(TMM,TCC,WMM,WCC)
C
C  THIS SUBROUTINE READS REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C    IFILE=1  H ON FE
C    IFILE=2  D ON FE
C    IFILE=3  H ON C
C    IFILE=4  D ON C
C    IFILE=5  HE ON FE
C    IFILE=6  HE ON C
C    IFILE=7  T ON FE
C    IFILE=8  T ON C
C    IFILE=9  D ON W
C    IFILE=10 HE ON W
C    IFILE=11 H ON W
C    IFILE=12 T ON W
C
      USE PRECISION
      USE PARMMOD
      USE CREFMOD
      USE CREF
      USE CSPEI

      IMPLICIT NONE
C
      REAL(DP), INTENT(OUT) :: TMM(*), TCC(*), WMM(*), WCC(*)
      REAL(DP) :: TML(12), TCL(12), WML(12), WCL(12),
     .          FELD(1092)
      REAL(DP) :: PID180
      INTEGER :: I1, I2, I3, I4, I5, NRECL, IUN, IFILE, I, J
C
      NFLR=NHD6
      IF (NHD6.GT.12) THEN
        WRITE (6,*) 'STORAGE ERROR. NHD6 MUST BE LESS OR EQUAL 12 '
        WRITE (6,*) 'NHD6= ',NHD6
        CALL EXIT
      ENDIF
C
      NRECL=1092
      INE=12
      INEM=INE-1
      TML(1)=1.
      TCL(1)=1.
      WML(1)=56.
      WCL(1)=26.
      TML(2)=2.
      TCL(2)=1.
      WML(2)=56.
      WCL(2)=26.
      TML(3)=1.
      TCL(3)=1.
      WML(3)=12.
      WCL(3)=6.
      TML(4)=2.
      TCL(4)=1.
      WML(4)=12.
      WCL(4)=6.
      TML(5)=4.
      TCL(5)=2.
      WML(5)=56.
      WCL(5)=26.
      TML(6)=4.
      TCL(6)=2.
      WML(6)=12.
      WCL(6)=6.
      TML(7)=3.
      TCL(7)=1.
      WML(7)=56.
      WCL(7)=26.
      TML(8)=3.
      TCL(8)=1.
      WML(8)=12.
      WCL(8)=6.
      TML(9)=2.
      TCL(9)=1.
      WML(9)=184.
      WCL(9)=74.
      TML(10)=4.
      TCL(10)=2.
      WML(10)=184.
      WCL(10)=74.
      TML(11)=1.
      TCL(11)=1.
      WML(11)=184.
      WCL(11)=74.
      TML(12)=3.
      TCL(12)=1.
      WML(12)=184.
      WCL(12)=74.
      DO 10 J=1,NHD6
        TMM(J)=TML(J)
        TCC(J)=TCL(J)
        WMM(J)=WML(J)
        WCC(J)=WCL(J)
10    CONTINUE
      ENAR(1)=1.
      ENAR(2)=2.
      ENAR(3)=5.
      ENAR(4)=10.
      ENAR(5)=20.
      ENAR(6)=50.
      ENAR(7)=100.
      ENAR(8)=200.
      ENAR(9)=500.
      ENAR(10)=1000.
      ENAR(11)=2000.
      ENAR(12)=5000.
      DO 11 I=1,INEM
11      DENAR(I)=1./(ENAR(I+1)-ENAR(I))
      INW=7
      INWM=INW-1
      PID180=ATAN(1.)/45.
      WIAR(1)=COS(0.D0)
      WIAR(2)=COS(30.*PID180)
      WIAR(3)=COS(45.*PID180)
      WIAR(4)=COS(60.*PID180)
      WIAR(5)=COS(70.*PID180)
      WIAR(6)=COS(80.*PID180)
      WIAR(7)=COS(85.*PID180)
      DO 12 I=1,INWM
12      DWIAR(I)=1./(WIAR(I+1)-WIAR(I))
C
      INR=5
      INRM=INR-1
      RAAR(1)=0.1
      RAAR(2)=0.3
      RAAR(3)=0.5
      RAAR(4)=0.7
      RAAR(5)=0.9
      DO 15 I=1,INRM
15      DRAAR(I)=1./(RAAR(I+1)-RAAR(I))
C
      IF (INE*INW*NFLR.GT.NH0 .OR.
     .    INE*INW*INR*NFLR.GT.NH1  .OR.
     .    INE*INW*INR*INR*NFLR.GT.NH2  .OR.
     .    INE*INW*INR*INR*INR*NFLR.GT.NH3) THEN
        WRITE (6,*) 'ERROR IN PARAMETER STATEMENT FOR REFLECTION DATA'
        CALL EXIT
      ENDIF
C
      IUN=21
      REWIND IUN
C
660   FORMAT (1X,1P,10E12.4)
661   FORMAT (4E20.12)
      DO 1 IFILE=1,NFLR
        DO 2 I1=1,INE
          I=0
          READ (IUN,661) (FELD(J),J=1,NRECL)
          DO 3 I2=1,INW
            I=I+1
            HFTR0(I1,I2,IFILE)=FELD(I)
3         CONTINUE
          DO 4 I2=1,INW
            DO 4 I3=1,INR
              I=I+1
              HFTR1(I1,I2,I3,IFILE)=FELD(I)
4         CONTINUE
          DO 5 I2=1,INW
            DO 5 I3=1,INR
              DO 5 I4=1,INR
                I=I+1
                HFTR2(I1,I2,I3,I4,IFILE)=FELD(I)
5         CONTINUE
          DO 6 I2=1,INW
            DO 6 I3=1,INR
              DO 6 I4=1,INR
                DO 6 I5=1,INR
                  I=I+1
                  HFTR3(I1,I2,I3,I4,I5,IFILE)=FELD(I)
6         CONTINUE
2       CONTINUE
1     CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE SLREAC (IR,FILNAM,H123,REAC,CRC)
c
C  input
C    FILNAM: read a&m data from file filnam, e.g. AMJUEL, HYDHEL, METHAN, CONST
c    IR    : store data on eirene array CREAC(...,...,IR)
c    H123  : identifyer for data type in filnam, e.g. H.1, H.2, H.3, ...
c    REAC  : number of reaction in filnam, e.g. 2.2.5
c    CRC   : type of process, e.g. EI, CX, OT, etc
C  internal 
C    ISW   <-- H123
C    IO    derived from ISW, initial value of 2nd index in CREAC
C  output
c    ISWR  : eirene flag for type of process  (1,2,...7)
c    CREAC : eirene storage array for a&m data CREAC(9,0:9,IR)
c    MODCLF: see below
c    DELPOT: ionisation potential (for H.10 data),
c            currently handeled in input.f. not nice!
c    IFTFLG: eirene flag for type of fitting  ("fit-flag=...")
C            DEFAULTS: =2 FOR POTENTIAL (GEN. MORSE)
C                      =0 FOR ALL OTHERS (POLYNOM, DOUBLE POLYNOM)
c
C  READ A&M DATA FROM THE FILES INTO EIRENE ARRAY CREAC
C
C
C  OUTPUT (IN COMMON COMXS):
C    READ DATA FROM "FILNAM" INTO ARRAY "CREAC"
C    DEFINE PARAMETER MODCLF(IR) (5 DIGITS NMLKJ)
C    FIRST DEZIMAL  J           =1  POTENTIAL AVAILABLE
C                                   (ON CREAC(..,-1,IR))
C                   J           =0  ELSE
C    SECOND DEZIMAL K           =1  CROSS SECTION AVAILABLE
C                                   (ON CREAC(..,0,IR))
C                   K           =0  ELSE
C    THIRD  DEZIMAL L           =1  <SIGMA V> FOR ONE
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
C    FOURTH DEZIMAL M               DATA FOR MOMENTUM EXCHANGE
C                                   TO BE WRITTEN
C    FIFTH  DEZIMAL N           =1  DELTA E FOR ONE PARAMETER E (E.G.
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
      USE PRECISION
      USE PARMMOD
      USE COMPRT
      USE COMXS

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: IR
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(9), INTENT(IN) :: REAC
      CHARACTER(3), INTENT(IN) :: CRC
      REAL(DP) :: CONST
      INTEGER :: I, IND, J, K, IH, I0P1, I0, IC, IREAC, ISW, INDFF, 
     .           IFLG, INC
      CHARACTER(80) :: ZEILE
      CHARACTER(6) :: AMJUEL, HYDHEL, H2VIBR
      CHARACTER(7) :: METHANE
      CHARACTER(2) :: CHR
      CHARACTER(3) :: CHRL, CHRR
      LOGICAL :: LCONST,LGEMIN,LGEMAX
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
C  ADD ONE MORE BLANK, IF POSSIBLE
      IF (IREAC.LT.9) IREAC=IREAC+1
C  H.0
      IF (ISW.EQ.0) THEN
        CHR='p0'
        CHRL='pl0'
        CHRR='pr0'
        I0=-1
        MODCLF(IR)=MODCLF(IR)+1
        IFLG=0
C  DEFAULT POTENTIAL: GENERALISED MORSE
        IFTFLG(IR,IFLG)=2
C  H.1
      ELSEIF (ISW.EQ.1) THEN
        CHR='a0'
        CHRL='al0'
        CHRR='ar0'
        I0=0
        MODCLF(IR)=MODCLF(IR)+10
        IFLG=1
C  DEFAULT CROSS SECTION: 8TH ORDER POLYNOM OF LN(SIGMA)
        IFTFLG(IR,IFLG)=0
C  H.2
      ELSEIF (ISW.EQ.2) THEN
        CHR='b0'
        CHRL='bl0'
        CHRR='br0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+100
        IFLG=2
C  DEFAULT RATE COEFFICIENT: 8TH ORDER POLYNOM OF LN(<SIGMA V>) FOR E0=0.
        IFTFLG(IR,IFLG)=0
C  H.3
      ELSEIF (ISW.EQ.3) THEN
        MODCLF(IR)=MODCLF(IR)+200
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.4
      ELSEIF (ISW.EQ.4) THEN
        MODCLF(IR)=MODCLF(IR)+300
        I0=1
        IFLG=2
        IFTFLG(IR,IFLG)=0
C  H.5
      ELSEIF (ISW.EQ.5) THEN
        CHR='e0'
        CHRL='el0'
        CHRR='er0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+1000
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.6
      ELSEIF (ISW.EQ.6) THEN
        MODCLF(IR)=MODCLF(IR)+2000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.7
      ELSEIF (ISW.EQ.7) THEN
        MODCLF(IR)=MODCLF(IR)+3000
        I0=1
        IFLG=3
        IFTFLG(IR,IFLG)=0
C  H.8
      ELSEIF (ISW.EQ.8) THEN
        CHR='h0'
        CHRL='hl0'
        CHRR='hr0'
        I0=1
        MODCLF(IR)=MODCLF(IR)+10000
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.9
      ELSEIF (ISW.EQ.9) THEN
        MODCLF(IR)=MODCLF(IR)+20000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.10
      ELSEIF (ISW.EQ.10) THEN
        MODCLF(IR)=MODCLF(IR)+30000
        I0=1
        IFLG=4
        IFTFLG(IR,IFLG)=0
C  H.11
      ELSEIF (ISW.EQ.11) THEN
        CHR='k0'
        CHRL='kl0'
        CHRR='kr0'
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
C  H.12
      ELSEIF (ISW.EQ.12) THEN
        I0=1
        IFLG=5
        IFTFLG(IR,IFLG)=0
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
C  SINGLE PARAM. FIT, ISW=0,1,2,5,8,11
      IF (ISW.EQ.0.OR.ISW.EQ.1.OR.ISW.EQ.2.OR.ISW.EQ.5.OR.ISW.EQ.8.OR.
     .    ISW.EQ.11) THEN
        IF (.NOT.LCONST) THEN
3         READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,CHR)+INDFF.EQ.0) GOTO 3
          IF (INDFF > 0) THEN
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 3
          ENDIF
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
          IF (ISW.EQ.0) GOTO 12 ! NO ASYMPTOTICS FOR POTENTIALS
          IF (INDEX(ZEILE,CHRL).NE.0.AND.IFEXMN(IR,I0P1).EQ.0) THEN
            IND=0
            DO 5 I=1,3
              INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
              IF (INC.GT.0) THEN
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
              ENDIF
5           CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0.AND.IFEXMX(IR,I0P1).EQ.0) THEN
            IND=0
            DO 7 I=4,6
              INC=INDEX(ZEILE((IND+1):80),CHR(1:1))
              IF (INC.GT.0) THEN
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
              ENDIF
7           CONTINUE
            LGEMAX=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
c
          if (lgemin.and.ifexmn(ir,I0P1).eq.0) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmn(IR,I0P1)
            rcmn(ir,I0P1)=log(rcmn(ir,I0P1))
            ifexmn(ir,I0P1)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
          if (lgemax.and.ifexmx(ir,I0P1).eq.0) then
            IND=INDEX(ZEILE,'=')
            READ (ZEILE((IND+2):80),'(E12.5)') rcmx(IR,I0P1)
            rcmx(ir,I0P1)=log(rcmx(ir,I0P1))
            ifexmx(ir,I0P1)=5
            READ (29,'(A80)',END=990) ZEILE
          endif
C
C  ANY OTHER ASYMPTOTICS INFO ON FILE?  SEARCH FOR Tmin, or Emin
          IF ((INDEX(ZEILE,'Tmin').NE.0.and.I0P1==2).or.
     .        (INDEX(ZEILE,'Emin').NE.0.and.I0P1==1)) then
            IND=INDEX(ZEILE,'n')
            READ (ZEILE((IND+2):80),'(E9.2)') rcmn(IR,I0P1)
            rcmn(ir,I0P1)=log(rcmn(ir,I0P1))
C  extrapolation from subr. CROSS
            if (I0P1.eq.1.and.iswr(ir).eq.1) ifexmn(ir,1)=1
            if (I0P1.eq.1.and.iswr(ir).eq.3) ifexmn(ir,1)=-1
            if (I0P1.eq.1.and.iswr(ir).eq.5) ifexmn(ir,1)=-1
C  extrapolation from subr. CDEF
C   ??      if (I0PT.eq.2) ifexmn(ir,1)=-1
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
12        CONTINUE
        ENDIF
C
C  TWO PARAM. FIT, ISW=3,4,6,7,9,10,12
      ELSEIF (ISW.EQ.3.OR.ISW.EQ.4.OR.ISW.EQ.6.OR.ISW.EQ.7.OR.
     .        ISW.EQ.9.OR.ISW.EQ.10.OR.ISW.EQ.12) THEN
        DO 11 J=0,2
16        READ (29,'(A80)',END=990) ZEILE
          INDFF=INDEX(ZEILE,'fit-flag')
          IF (INDEX(ZEILE,'Index')+INDFF.EQ.0) GOTO 16
          IF (INDFF > 0) THEN
            READ (ZEILE((INDFF+8):80),*) IFTFLG(IR,IFLG)
            GOTO 16
          ENDIF
          READ (29,'(1X)')
          DO 17 I=1,9
            READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
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
      CLOSE (UNIT=29)
      CALL EXIT
6664  FORMAT (6E12.4)
      END
c === ROUTINE: wrgeom
C
      SUBROUTINE WRGEOM(TRCFLE)
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CPOLYG
      USE CGRID
      USE CGEOM
      USE CLGIN
      USE CTRIG
      IMPLICIT NONE
      LOGICAL TRCFLE
C
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      WRITE (12) RCGM1,RCGM2,NPOINT,NSTGRD,NGHPLS,NGHPOL,NCLTAL
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCGM,ICGM '
      WRITE (12) RCGRID,ICGRID
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCGRID,ICGRID'
      WRITE (12) RCPLYG,RCPLY2,ICPLYG
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCPLYG,ICPLYG'
      WRITE (12) XTRIAN,YTRIAN,VTRIX,VTRIY,PTRIX,PTRIY,
     .           NECKE,NCHBAR,NSEITE,INMTI,NRKNOT,NTRII
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCTRIG,ICTRIG'
      WRITE (12) 
     R RLWMN,RLWMX,EWALL,EWBIN,TRANSP,FSHEAT,
     R ZNML,ZNCL,
     R RECYCF,RECYCT,RECPRM,EXPPL,EXPEL,EXPIL,
     R RECYCS,RECYCC,SPTPRM,
     I ILSWCH,ILEQUI,ILTOR,ILSIDE,ILIIN,ILREF,
     I ILSPT,ILCOL,ILFIT,ILCELL,ILBOX,ILPLG,ISPUT,
     I NLIMII,NLIMIE,ISWICH,ILBLCK,ILACLL,JUMLIM,
     I NSTSI,INUMP,IRPTA,IRPTE,ISRF,ISRT,ISRS,ISRC,
     I INMP1I,INMP2I,INMP3I,
     I IGFIL,IGJUM0,IGJUM1,IGJUM2,IGJUM3
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCLGN,ICLGN,LCLGN'
      WRITE (12) RADGEO,IADGEO,NLIMI
      IF (TRCFLE) WRITE (6,*) 'WRITE 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
C
      ENTRY RGEOM(TRCFLE)
      OPEN (UNIT=12,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 12
      READ (12) RCGM1,RCGM2,NPOINT,NSTGRD,NGHPLS,NGHPOL,NCLTAL
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCGM,ICGM '
      READ (12) RCGRID,ICGRID
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCGRID,ICGRID'
      READ (12) RCPLYG,RCPLY2,ICPLYG
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCPLYG,ICPLYG'
      READ (12) XTRIAN,YTRIAN,VTRIX,VTRIY,PTRIX,PTRIY,
     .           NECKE,NCHBAR,NSEITE,INMTI,NRKNOT,NTRII
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCTRIG,ICTRIG'
      READ (12) 
     R RLWMN,RLWMX,EWALL,EWBIN,TRANSP,FSHEAT,
     R ZNML,ZNCL,
     R RECYCF,RECYCT,RECPRM,EXPPL,EXPEL,EXPIL,
     R RECYCS,RECYCC,SPTPRM,
     I ILSWCH,ILEQUI,ILTOR,ILSIDE,ILIIN,ILREF,
     I ILSPT,ILCOL,ILFIT,ILCELL,ILBOX,ILPLG,ISPUT,
     I NLIMII,NLIMIE,ISWICH,ILBLCK,ILACLL,JUMLIM,
     I NSTSI,INUMP,IRPTA,IRPTE,ISRF,ISRT,ISRS,ISRC,
     I INMP1I,INMP2I,INMP3I,
     I IGFIL,IGJUM0,IGJUM1,IGJUM2,IGJUM3
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCLGN,ICLGN,LCLGN'
      READ (12) RADGEO,IADGEO,NLIMI
      IF (TRCFLE) WRITE (6,*) 'READ 12: RCADG,ICADG'
      CLOSE (UNIT=12)
      RETURN
      END
c === ROUTINE: wrmesh
      SUBROUTINE WRMESH
      USE PRECISION
      USE PARMMOD
      USE CADGEO
      USE CPLOT
      USE CPOLYG
      USE CGEOM
      USE CLGIN
      IMPLICIT NONE


      INTEGER, PARAMETER :: MAXPOIN=300
      REAL(DP) :: CONT(MAXPOIN,2), partcont(maxpoin,2,2), maxlen
      REAL(DP) :: XPE, YPE, HELP, XT, YT, PHI1, X1, X2, Y1, Y2, PHI2
      INTEGER  :: MAXCONT, ICONT, IPOIN, IWST, IWEN, IWL, IWP,
     .            IWAN, IMN, I, NCONT, J, IUHR, ISTORE, IP, IH, IFOUND
      INTEGER  :: IDIAG(MAXPOIN)
      REAL(SP) :: xmin,xmax,ymin,ymax,deltax,deltay,delta,xcm,ycm
      REAL(SP) :: XP,YP

C INITIALISIERUNG DER PLOTDATEN
      xmin = CH2X0-CH2MX
      ymin = CH2Y0-CH2MY
      xmax = CH2X0+CH2MX
      ymax = CH2Y0+CH2MY
      deltax = abs(xmax-xmin)
      deltay = abs(ymax-ymin)
      delta = max(deltax,deltay)
      xcm = 24. * deltax/delta
      ycm = 24. * deltay/delta

C ANZAHL DER KONTOUREN BESTIMMEN
C ILPLG WIRD IM INPUT-BLOCK 3 EINGELESEN
      CALL LEER(2)
      WRITE (6,*) 'SUBROUTINE WRMESH CALLED '
      CALL LEER(1)
      NCONT = 0
      DO I=1,NLIMI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO
      DO I=NLIM+1,NLIM+NSTSI
        NCONT = MAX(NCONT,ABS(ILPLG(I)))
      ENDDO

      WRITE (6,*) 'NUMBER OF CONTOURS FOR FEM MESH: ',NCONT
      CALL LEER(1)
      if (ncont == 0) return

      call grnxtf
      call grsclc(3.,3.,3.+real(xcm,kind(1.e0)),3.+real(ycm,kind(1.e0)))
      call grsclv(real(xmin,kind(1.e0)),real(ymin,kind(1.e0)),
     .            real(xmax,kind(1.e0)),real(ymax,kind(1.e0)))

      DO ICONT = 1,NCONT
        IPOIN = 0
        MAXLEN = 0.
C AKTUELLE KONTOUR BESTIMMEN, STUECKE MIT ILPLG=ICONT GEHOEREN ZUR
C AKTUELLEN CONTOUR, ANFANGS UND ENDPUNKT DIESES STUECKES WERDEN AUF
C PARTCONT GESPEICHERT
        DO I=1,NLIMI
          IF (ABS(ILPLG(I)) .EQ. ICONT) THEN
            IUHR=ILPLG(I)
C 0 < RLB(I) < 2
C 2-PUNKT OPTION WIRD IM TIMEA0 AUF RLB=1 ZURUECKGEFUEHRT
            IF ((RLB(I) .GT. 0.) .AND. (RLB(I) .LT. 2.) .AND.
     >          (P3(1,I) .EQ. 1.D55 .OR. P3(2,I) .EQ. 1.D55
     >          .OR. P3(3,I) .EQ. 1.D55)) THEN
              IPOIN = IPOIN + 1
              IF (A3LM(I) .EQ. 0.) THEN
C               X,Y-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(2,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(2,I)
                idiag(ipoin)=i
              ELSEIF (A2LM(I) .EQ. 0.) THEN
C               X,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(1,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(1,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ELSEIF (A1LM(I) .EQ. 0.) THEN
C               Y,Z-KOORDINATEN
                PARTCONT(IPOIN,1,1) = P1(2,I)
                PARTCONT(IPOIN,1,2) = P1(3,I)
                PARTCONT(IPOIN,2,1) = P2(2,I)
                PARTCONT(IPOIN,2,2) = P2(3,I)
                idiag(ipoin)=i
              ENDIF
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
            ELSE
C  ERROR
              WRITE(6,*) 'FALSCHE ANGABE FUER RLB, RLB = ',RLB(I),
     >                    ILPLG(I),I
            ENDIF
          ENDIF
        ENDDO

        DO I=1,NSTSI
          IF (ABS(ILPLG(NLIM+I)) .EQ. ICONT) THEN
            IUHR=ILPLG(NLIM+I)
            IF (INUMP(I,2) .NE. 0) THEN
C             POLOIDAL
              DO J=IRPTA(I,1),IRPTE(I,1)-1
c slmod begin
c...            Some development required:
                IF (GRIDOPT.EQ.1) THEN
                  WRITE(0,*) 'GRIDOPT: OPTION NOT SUPPORTED 001'
                  STOP
                ENDIF
c slmod end
                IF ((XPOL(J,INUMP(I,2)) .NE. XPOL(J+1,INUMP(I,2))) .OR.
     >              (YPOL(J,INUMP(I,2)) .NE. YPOL(J+1,INUMP(I,2)))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,1,2) = YPOL(J,INUMP(I,2))
                  PARTCONT(IPOIN,2,1) = XPOL(J+1,INUMP(I,2))
                  PARTCONT(IPOIN,2,2) = YPOL(J+1,INUMP(I,2))
                idiag(ipoin)=-i
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSEIF (INUMP(I,1) .NE. 0) THEN
C             RADIAL
              DO J=IRPTA(I,2),IRPTE(I,2)-1
c slmod begin
c...            Some (more) development required:
                IF (GRIDOPT.EQ.1) THEN
                  WRITE(0,*) 'GRIDOPT: OPTION NOT SUPPORTED 002'
                  STOP
                ENDIF
c slmod end
                IF ((XPOL(INUMP(I,1),J) .NE. XPOL(INUMP(I,1),J+1)) .OR.
     >              (YPOL(INUMP(I,1),J) .NE. YPOL(INUMP(I,1),J+1))) THEN
                  IPOIN = IPOIN + 1
                  PARTCONT(IPOIN,1,1) = XPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,1,2) = YPOL(INUMP(I,1),J)
                  PARTCONT(IPOIN,2,1) = XPOL(INUMP(I,1),J+1)
                  PARTCONT(IPOIN,2,2) = YPOL(INUMP(I,1),J+1)
                idiag(ipoin)=-i
              maxlen = maxlen +
     >               sqrt((partcont(ipoin,1,1)-partcont(ipoin,2,1))**2
     >                   +(partcont(ipoin,1,2)-partcont(ipoin,2,2))**2)
                ENDIF
              ENDDO
            ELSE
C  ERROR
              WRITE(6,*) 'CASE NOT FORESEEN: INUMP: ',
     >                     (INUMP(I,J),J=1,3)
            ENDIF
          ENDIF
        ENDDO
        IF (IPOIN.LE.0) THEN
          WRITE(6,*) 'CONTOUR ',ICONT,' NOT FOUND'
          GOTO 1000
        ENDIF

        call grnwpn(icont)
C STUECKE DER AKTUELLEN KONTOUR WERDEN SORTIERT
        DO I=1,IPOIN-1
          XPE = PARTCONT(I,2,1)
          YPE = PARTCONT(I,2,2)
          IFOUND=0
          DO J=I+1,IPOIN
            IF ((XPE .EQ. PARTCONT(J,1,1)) .AND.
     >          (YPE .EQ. PARTCONT(J,1,2))) THEN
              IFOUND=1
              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih
            ELSEIF ((XPE .EQ. PARTCONT(J,2,1)) .AND.
     >              (YPE .EQ. PARTCONT(J,2,2))) THEN
              IFOUND=1
              HELP = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              HELP = PARTCONT(I+1,1,1)
              PARTCONT(I+1,1,1) = PARTCONT(J,1,1)
              PARTCONT(J,1,1) = HELP
              HELP = PARTCONT(I+1,1,2)
              PARTCONT(I+1,1,2) = PARTCONT(J,1,2)
              PARTCONT(J,1,2) = HELP

              HELP = PARTCONT(I+1,2,1)
              PARTCONT(I+1,2,1) = PARTCONT(J,2,1)
              PARTCONT(J,2,1) = HELP
              HELP = PARTCONT(I+1,2,2)
              PARTCONT(I+1,2,2) = PARTCONT(J,2,2)
              PARTCONT(J,2,2) = HELP

              ih=idiag(i+1)
              idiag(i+1)=idiag(j)
              idiag(j)=ih
            ENDIF
          ENDDO
          IF (IFOUND.EQ.0) THEN
            WRITE (6,*) 'NO MATCHING POINT FOUND FOR CONTOUR ',ICONT
            write(6,*) i,idiag(i),partcont(i,1,1),partcont(i,1,2),
     >                            partcont(i,2,1),partcont(i,2,2)
            WRITE (6,*) 'USE NEXT POINT '
            IP=I+1
            write(6,*) iP,idiag(iP),partcont(iP,1,1),partcont(iP,1,2),
     >                              partcont(iP,2,1),partcont(iP,2,2)
          ENDIF
        ENDDO
        IF ((PARTCONT(1,1,1) .NE. PARTCONT(IPOIN,2,1)) .OR.
     >      (PARTCONT(1,1,2) .NE. PARTCONT(IPOIN,2,2))) THEN
          WRITE(6,*) 'CONTOUR ',ICONT,' IS NOT CLOSED'
        ELSE
          WRITE(6,*) 'CLOSED CONTOUR ',ICONT
        ENDIF
        do i=1,ipoin
          write(6,*) i,idiag(i),partcont(i,1,1),partcont(i,1,2),
     >                          partcont(i,2,1),partcont(i,2,2)
        enddo
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grjmp(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        DO I=2,IPOIN
          XP = PARTCONT(I,1,1)
          YP = PARTCONT(I,1,2)
          call grdrw(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))
        ENDDO
        XP = PARTCONT(1,1,1)
        YP = PARTCONT(1,1,2)
        call grDRW(REAL(XP,KIND(1.E0)),REAL(YP,KIND(1.E0)))

C  BERECHNUNG VON DELTA ALS MITTLERE LAENGE DER TEILSTUECKE
C  DELTA IST MASS FUER DIE GROESSE DER DREIECKE
        IMN=0
        IF (ICONT .EQ. 1) THEN
          maxlen = maxlen / REAL(IPOIN,KIND(1.D0))
          WRITE(78,*) maxlen
          WRITE(78,*)
        endif

C BESTIMMUNG DES UHRZEIGERSINNS DER KONTOUR
C KONTOUR MUSS FUER DIE TRIANGULIERUNG FOLGENDERMASSEN AUSGEGEBEN
C WERDEN:
C  - IM UHRZEIGERSINN FUER INNERE BEGRENZUNGEN DES GEBIETES (POSITIV)
C  - GEGEN UHRZEIGERSINN FUER AEUSSERE BEGRENZUNGEN DES GEBIETES (NEGATIV)
        YMIN=PARTCONT(1,1,2)
        DO I=1,IPOIN
          IF (PARTCONT(I,2,2) .LT. YMIN) THEN
            YMIN = PARTCONT(I,2,2)
            IMN=I
          ENDIF
        ENDDO
 
        IF (IMN .EQ. 0) THEN
          XT = PARTCONT(1,1,1)
          YT = PARTCONT(1,1,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IPOIN,1,1)
          Y1 = PARTCONT(IPOIN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(1,2,1)
          Y2 = PARTCONT(1,2,2)
        ELSE
C  SONDERFALL IMN=IPOIN ENTFAELLT, DA ERSTER PUNKT GLEICH LETZTER
C  PUNKT GILT
          XT = PARTCONT(IMN,2,1)
          YT = PARTCONT(IMN,2,2)
C  PUNKT, DER IM UMLAUF DER VORHERGEHENDE IST
          X1 = PARTCONT(IMN,1,1)
          Y1 = PARTCONT(IMN,1,2)
C  PUNKT, DER IM UMLAUF DER NAECHSTE IST
          X2 = PARTCONT(IMN+1,2,1)
          Y2 = PARTCONT(IMN+1,2,2)

        ENDIF
 
C  BESTIMME POLARWINKEL VON (X1,Y1) UND (X2,Y2) MIT (XT,YT) ALS URSPRUNG
        PHI1 = ATAN2 (Y1-YT,X1-XT)
        PHI2 = ATAN2 (Y2-YT,X2-XT)

        IF (PHI2 .GT. PHI1) THEN
C  ABSPEICHERUNG ERFOLGTE IM UHRZEIGERSINN
          ISTORE = 1
        ELSE
          ISTORE = -1
        ENDIF
C  IUHR=ILPLG > 0 ==> IM UHRZEIGERSINN AUSGEBEN
C  IUHR=ILPLG < 0 ==> ENTGEGEN DEM UHRZEIGERSINN AUSGEBEN
        IWAN=1
        IWEN=IPOIN
        IWST=1
        IWP=1
        IWL=2
        IF (ISTORE*IUHR .LT. 0.) THEN
          IWAN=IPOIN
          IWEN=1
          IWST=-1
          IWP=2
          IWL=1
        ENDIF
        WRITE(78,*) IPOIN+1
        DO I=iwan,iwen,iwst
          WRITE(78,'(1P,2(2X,E21.14))')
     >          PARTCONT(I,IWP,1),PARTCONT(I,IWP,2)
        ENDDO
        WRITE(78,'(1P,2(2X,E21.14))') PARTCONT(IWEN,IWL,1),
     >                               PARTCONT(IWEN,IWL,2)
1000    CONTINUE
      ENDDO
C  END OF NCONT LOOP
      call leer(1)
      write (6,*) 'input file fort.78 for FEM mesh generator written '
      call leer(2)
      call grnwpn(1)
      call grnxtf
      END
c === ROUTINE: algebr
C
C
      SUBROUTINE WRPLAM(TRCFLE)
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CZT1
      USE COMSOU
      USE CSTEP
      USE COMXS
      IMPLICIT NONE
      LOGICAL TRCFLE
C
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      WRITE (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,
     R           DION,DATM,DMOL,DPLS,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMUSR,ICMUSR,LCMUSR '
      CALL WRITE_CMDTA
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMDTA,ICMDTA'
      CALL WRITE_CMAMF
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMAMF,ICMAMF'
      WRITE (13) RCZT1,ZT1,ZRG
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCZT1,ZT1,ZRG'
      WRITE (13) RCMSOU,SREC,EIO,EEL,
     .           ICMSOU,INGRDA,INGRDE,NSTRAI,
     .           LCMSOU,NLSYMP,NLSYMT
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCMSOU,ICMSOU,LCMSOU'
      WRITE (13) FLSTEP,ELSTEP,FLTOT,VF,QUOT,ADD,QUOTI,ADDIV,
     .           TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .           IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .           ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (6,*) 'WRITE 13: RCSTEP,ICSTEP'
      CLOSE (UNIT=13)
      RETURN
C
      ENTRY RPLAM(TRCFLE)
      OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 13
      READ (13)
C  REAL
     R           TEIN,TIIN,DEIN,DIIN,VXIN,VYIN,VZIN,
     R           BXIN,BYIN,BZIN,BFIN,ADIN,EDRIFT,
     R           VOL,WGHT,
     R           FLXOUT,SAREA,
     R           TEINL,TIINL,DEINL,DIINL,BVIN,PARMOM,
     R           RMASSI,RMASSA,RMASSM,RMASSP,
     R           DIOD,DATD,DMLD,DPLD,
     R           DION,DATM,DMOL,DPLS,
     R           TVAC,DVAC,VVAC,ALLOC,
     T           TEXTS,
C  MUSR, INTEGER
     I           NSPH  ,NPHOTI,NPHOTIM,NFOLPH,NGENPH,
     I           NSPA  ,NATMI,NATMIM,NMASSA,NCHARA,NFOLA,NGENA,
     I           NSPAM ,NMOLI,NMOLIM,NMASSM,NCHARM,NFOLM,NGENM,
     I           NSPAMI,NIONI,NIONIM,NMASSI,NCHARI,NCHRGI,NFOLI,NGENI,
     I           NSPTOT,NPLSI,NPLSIM,NMASSP,NCHARP,NCHRGP,NBITS,
     I           NSNVI,NCPVI,NADVI,NBGVI,NALVI,NCLVI,NADSI,NALSI,NAINI,
     I           NPRT,ISPEZ,ISPEZI,
C  LUSR, LOGICAL
     L           LGVAC,LGDFT
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMUSR,ICMUSR,LCMUSR '
      CALL READ_CMDTA
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMDTA,ICMDTA'
      CALL READ_CMAMF
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMAMF,ICMAMF'
      READ (13) RCZT1,ZT1,ZRG
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCZT1,ZT1,ZRG'
      READ (13) RCMSOU,SREC,EIO,EEL,
     .          ICMSOU,INGRDA,INGRDE,NSTRAI,
     .          LCMSOU,NLSYMP,NLSYMT
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCMSOU,ICMSOU,LCMSOU'
      READ (13) FLSTEP,ELSTEP,FLTOT,VF,QUOT,ADD,QUOTI,ADDIV,
     .           TESTEP,TISTEP,RRSTEP,VXSTEP,VYSTEP,VZSTEP,DISTEP,
     .           IRSTEP,IPSTEP,ITSTEP,IASTEP,IBSTEP,IGSTEP,
     .           ISTUF,NSMAX,NSPSTI,NSPSTE
      IF (TRCFLE) WRITE (6,*) 'READ 13: RCSTEP,ICSTEP'
      CLOSE (UNIT=13)
      RETURN
      END
C
      SUBROUTINE WRREC
C
C  EVALUATE EIRENE RECOMMENDATIONS FOR A NEXT RUN OF THE SAME MODEL
C
      USE PRECISION
      USE PARMMOD
      USE CAI
      USE CCONA
      USE CTRCEI
      USE COMSOU
      USE COUTAU

      IMPLICIT NONE

      REAL(DP) :: WSUM, FTOT, XNSUM
      INTEGER :: NREQ, ISTRA
      REAL(DP) :: WTOTT(NSTRA),WMEAN(NSTRA),WREC(NSTRA),XNEXP(NSTRA),
     .          CPUFAC(NSTRA)

      OPEN (UNIT=14,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 14
C
C  FIRSTLY: STRATIFIED SOURCE SAMPLING
C
      IF (NSTRAI.EQ.1) THEN
        RATIO(1)=1.
        GOTO 350
      ENDIF
C
      WSUM=0.
      NREQ=0
      FTOT=0.
      DO 100 ISTRA=1,NSTRAI
        IF (XMCP(ISTRA).LE.0.D0) GOTO 100
        WTOTT(ISTRA)=-WTOTP(0,ISTRA)+WTOTA(0,ISTRA)+
     .                WTOTM(0,ISTRA)+WTOTI(0,ISTRA)
        WTOTT(ISTRA)=WTOTT(ISTRA)/(FLXFAC(ISTRA)+EPS60)
        WMEAN(ISTRA)=WTOTT(ISTRA)/(XMCP(ISTRA)+EPS60)
        WSUM=WSUM+WTOTT(ISTRA)
        NREQ=NREQ+NPTS(ISTRA)
        FTOT=FTOT+FLUXT(ISTRA)
100   CONTINUE
C
C  PROPORTIONAL ALLOCATION: RECOMMENDED REL. WEIGHT PER STRATUM: WREC
C                           EXPECTED REL. NO OF PARTICLES NEEDED: XNEXP
C                           RECOMMENDED NO. OF PARTICLES: NRECOM
C  NOT THE NO. OF PARTICLES BUT THE SUM OF BIRTH WEIGHTS PER STRATUM
C  IS ALLOCATED  PROPORTIONAL TO THE RELATIVE STRATUM POPULATION
C
C
      XNSUM=0.
      DO 200 ISTRA=1,NSTRAI
        IF (XMCP(ISTRA).LE.0.D0) GOTO 200
        WREC(ISTRA)=FLUXT(ISTRA)/(FTOT+EPS60)
        XNEXP(ISTRA)=WREC(ISTRA)/(WMEAN(ISTRA)+EPS60)
C  ACCOUNT FOR DIFFERENT COMPUTING SPEED AT DIFFERENT STRATA
C  ASSUME THEREFORE THAT ALLOCATED CPU TIME WAS PROPORTIONAL NPTS(ISTRA)
        CPUFAC(ISTRA)=XMCP(ISTRA)/(DBLE(NPTS(ISTRA))+EPS60)
        XNEXP(ISTRA)=XNEXP(ISTRA)/(CPUFAC(ISTRA)+EPS60)
        XNSUM=XNSUM+XNEXP(ISTRA)
200   CONTINUE
      DO 300 ISTRA=1,NSTRAI
        RATIO(ISTRA)=0.
        IF (XMCP(ISTRA).LE.0.D0) GOTO 300
C  SCALE XNEXP TO CONSERVE TOTAL NUMBER OF REQUESTED TRACKS
        XNEXP(ISTRA)=XNEXP(ISTRA)*DBLE(NREQ)/(XNSUM+EPS60)
C  CONVERT XNEXP TO AN INTERGER
        NRECOM(ISTRA)=IDINT(REAL(XNEXP(ISTRA),KIND(1.D0)))
        IF (XNEXP(ISTRA)-DBLE(NRECOM(ISTRA)).GT.0.5)
     .      NRECOM(ISTRA)=NRECOM(ISTRA)+1
        RATIO(ISTRA)=WREC(ISTRA)/(WTOTT(ISTRA)/(WSUM+EPS60)+EPS60)
300   CONTINUE
C
350   CONTINUE
      IF (TRCFLE) WRITE (6,*) 'WRITE 14: RATIO,NRECOM '
      WRITE (14) RATIO,NRECOM
C
      IF (.NOT.TRCREC.OR.NSTRAI.EQ.1) GOTO 1000
      CALL PAGE
      WRITE (6,*) '=================================================='
      WRITE (6,*) '= RECOMMENDED INPUT MODIFICATIONS FOR STRATIFIED ='
      WRITE (6,*) '= SOURCE SAMPLING                                ='
      WRITE (6,*) '=================================================='
      CALL LEER(3)
      WRITE (6,*) 'NPTS OLD = OLD INPUT VALUE FOR NPTS '
      WRITE (6,*) 'NPTS REC = EIRENE RECOMMENDATION FOR NEXT RUN '
      WRITE (6,*) 'RATIO    = RECOM. REL.WEIGHT / ACTUAL REL. WEIGHT'
      WRITE (6,*) 'VALUES OF RATIO CLOSE TO ONE MEAN '
      WRITE (6,*) '"ALMOST PROPORTIONAL ALLOCATION OF WEIGHTS"'
      CALL LEER(2)
      WRITE (6,*) 'STRAT. NO, NPTS OLD, NPTS REC, RATIO '
      DO 400 ISTRA=1,NSTRAI
        WRITE (6,'(1X,I2,8X,I6,4X,I6,5X,1P,E12.4)')
     .                ISTRA,NPTS(ISTRA),NRECOM(ISTRA),RATIO(ISTRA)
400   CONTINUE
      CALL LEER(2)
1000  CONTINUE
C
C  STRATIFIED SOURCE SAMPLING ACCESSMENT FINISHED
C
C  SECONDLY: WEIGHT WINDOWS
C
C  TO BE WRITTEN
C
      RETURN
C
      ENTRY RREC
C
      OPEN (UNIT=14,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 14
      IF (TRCFLE) WRITE (6,*) 'READ 14: RATIO,NRECOM '
      READ (14) RATIO,NRECOM
C
      RETURN
      END
C
C
      SUBROUTINE WRSNAP
C
C  SAVE SNAPSHOT POPULATION AT END OF TIMESTEP
C
      USE PRECISION
      USE PARMMOD
      USE CTRCEI
      USE COMPRT
      USE COMNNL
      USE COMSOU

      IMPLICIT NONE

      INTEGER :: I, J
C
      OPEN (UNIT=15,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15
C
      IF (TRCFLE) WRITE (6,*) 'WRITE 15: IPRNL,FLUX,DTIMV '
      WRITE (15) IPRNL,FLUX(NSTRAI),DTIMV
      WRITE (15) ((RPARTC(I,J),J=1,NPARTT),I=1,IPRNL)
      WRITE (15)  (RPARTW(I)              ,I=0,IPRNL)
      WRITE (15) ((IPARTC(I,J),J=1,MPARTT),I=1,IPRNL)
C
      RETURN
C
      ENTRY RSNAP
C
      OPEN (UNIT=15,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
      REWIND 15
      READ (15) IPRNL,FLUX(NSTRAI),DTIMV
      IF (TRCFLE) WRITE (6,*) 'READ 15: IPRNL,FLUX,DTIMV '
      READ (15) ((RPARTC(I,J),J=1,NPARTT),I=1,IPRNL)
      READ (15)  (RPARTW(I)              ,I=0,IPRNL)
      READ (15) ((IPARTC(I,J),J=1,MPARTT),I=1,IPRNL)
C
      RETURN
      END
