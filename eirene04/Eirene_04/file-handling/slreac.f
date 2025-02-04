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
c    IFTFLG: eirene flag for type of fitting expression ("fit-flag=...")
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
      USE PHOTON

      IMPLICIT NONE

      INTEGER,      INTENT(IN) :: IR
      CHARACTER(8), INTENT(IN) :: FILNAM
      CHARACTER(4), INTENT(IN) :: H123
      CHARACTER(LEN=*), INTENT(IN) :: REAC
      CHARACTER(3), INTENT(IN) :: CRC
      CHARACTER(11) :: REACSTR
      REAL(DP) :: CONST
      INTEGER :: I, IND, J, K, IH, I0P1, I0, IC, IREAC, ISW, INDFF, 
     .           IFLG, INC, IANF
      CHARACTER(80) :: ZEILE
      CHARACTER(6) :: AMJUEL, HYDHEL, H2VIBR, SPECTR
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
      SPECTR='SPECTR  '
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
      ELSEIF (INDEX(FILNAM,'SPECTR').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='SPECTR')
      ELSEIF (INDEX(FILNAM,'PHOTON').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='PHOTON')
      ELSEIF (INDEX(FILNAM,'PHTNEW').NE.0) THEN
        LCONST=.FALSE.
        OPEN (UNIT=29,FILE='PHTNEW')
      ELSEIF (INDEX(FILNAM,'CONST').NE.0) THEN
        LCONST=.TRUE.
      ELSE
        WRITE (6,*) ' NO SPECIFICATION FOR FILENAME IN REACTION CARD'
        WRITE (6,*) ' CHOOSE EITHER '
        WRITE (6,*) ' AMJUEL, METHAN, HYDHEL, H2VIBR, SPECTR '
        WRITE (6,*) ' OR '
        WRITE (6,*) ' CONST FOR ENTERING REACTION COEFFICIENTS VIA '
        WRITE (6,*) ' EIRENE INPUT-FILE '
        CALL EXIT_OWN(1)
      ENDIF
C
      IF (H123(4:4).EQ.' ') THEN
        READ (H123(3:3),'(I1)') ISW
      ELSE
        READ (H123(3:4),'(I2)') ISW
      ENDIF

C
      IF (INDEX(FILNAM,'PHOTON').NE.0) THEN
        CALL READ_PHOTDBK (IR,REAC,ISW) 
        RETURN
!      ELSEIF (INDEX(FILNAM,'PHTNEW').NE.0) THEN
!        CALL READ_PHOTDBK_NEU (IR,REAC,ISW) 
!        RETURN
      END IF
      REACSTR=REPEAT(' ',11)
      IANF=VERIFY(REAC,' ')
      IF (IANF > 0) THEN
        IREAC=INDEX(REAC(IANF:),' ')-1
        IF (IREAC.LT.0) IREAC=LEN(REAC(IANF:))
        REACSTR(2:IREAC+1)=REAC(IANF:IREAC+IANF-1)
C  ADD ONE MORE BLANK, IF POSSIBLE
        IREAC=IREAC+2
      ELSE
        IF (.NOT.LCONST) THEN
           WRITE (6,*) ' NO REACTION SPECIFIED IN REACTION CARD ',IR
           CALL EXIT_OWN (1)
        END IF
      END IF

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
        IND=INDEX(REACSTR,'FT')
        IF (IND /= 0) THEN
          READ (REACSTR(IND+2:),*) IFTFLG(IR,IFLG)
        END IF

        IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
C
C  READ ONLY ONE FIT COEFFICIENT FROM INPUT FILE
          READ (IUNIN,6664) CREAC(1,I0,IR)
        ELSE
C
C  READ 9 FIT COEFFICIENTS FROM INPUT FILE
          READ (IUNIN,6664) (CREAC(IC,I0,IR),IC=1,9)
        END IF
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
     .      INDEX(ZEILE,REACSTR(1:ireac)).EQ.0) GOTO 2
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
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
            IND=INDEX(ZEILE,CHR(1:1))
            CREAC(:,I0,IR)=0._DP
            READ (ZEILE((IND+2):80),'(E20.12)') CREAC(1,I0,IR)
          ELSE
            DO 9 J=0,2
              IND=0
              DO 4 I=1,3
                IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
                READ (ZEILE((IND+2):80),'(E20.12)') CREAC(J*3+I,I0,IR)
4             CONTINUE
              READ (29,'(A80)',END=990) ZEILE
9           CONTINUE
          END IF
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
C       ELSEIF (LCONST) THEN
C  NOTHING TO BE DONE
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
          IF (MOD(IFTFLG(IR,IFLG),100) == 10) THEN
            CREAC(:,:,IR)=0._DP 
            READ (29,*) IH,CREAC(1,1,IR)
            EXIT
          ELSE
            DO 17 I=1,9
              READ (29,*) IH,(CREAC(I,K,IR),K=J*3+1,J*3+3)
17          CONTINUE
          END IF
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
      CALL EXIT_OWN(1)
991   WRITE (6,*) ' INVALID CONSTANT IN SLREAC. CONST= ',CONST
      WRITE (6,*) ' CHECK "REACTION CARDS" FOR REACTION NO. ',IR
      CLOSE (UNIT=29)
      CALL EXIT_OWN(1)
6664  FORMAT (6E12.4)
      END
