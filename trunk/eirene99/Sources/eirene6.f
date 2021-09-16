      SUBROUTINE XSTAPI(COUN,PLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CLOGAU'
      INCLUDE 'CTRCEI'
      INCLUDE 'CZT1'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'COMPRT'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      INCLUDE 'CSPEI'
      DIMENSION PLS(*),COUN(0:9,*),CF(9,0:9)
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
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.4) GOTO 130
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
            ADDT=PMASS/RMASSP(IPLS)
            ADDTL=LOG(ADDT)
            ADDPI(IRPI,IPLS) = ADDTL
C
C CROSS SECTION (E-LAB)
            IF (IDEZ(MODCLF(KK),1,4).EQ.1) THEN
              MODCOL(4,1,IATM,IPLS)=KK
              MODCOL(4,2,IATM,IPLS)=3
              IF (FACTKK.NE.1.D0)
     .        WRITE (6,*) 'FREACA NOT READY FOR CROSS SECTION IN XSTAPI'
            ENDIF
C
C RATE COEFFICIENT
            MODC=IDEZ(MODCLF(KK),2,4)
            IF (MODC.GE.1.AND.MODC.LE.2) THEN
              MODCOL(4,2,IATM,IPLS)=MODC
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
                CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
              END IF
            ENDIF
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
              MODCOL(4,4,IATM,IPLS)=1
            ELSEIF (NSEPI4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 252 J=1,NSBOX
                  EPLPI3(IRPI,J,1)=1.5*TIIN(IPLS,J)+EDRIFT(IPLS,J)
252             CONTINUE
              ELSE
                NELRPI(IRPI) = -2
              END IF
              MODCOL(4,4,IATM,IPLS)=1
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
C             MODCOL(4,4,IATM,IPLS)=1
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
      SUBROUTINE XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CGRID'
      INCLUDE 'CCONA'
      DIMENSION COUN(0:9,NRAD),CF(9,0:9),PLS(NRAD)
      CHARACTER*8 TEXT
c slmod begin - not tr (yet?)
      COMMON /MULCOM/ IOPT1
c slmod end
C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
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
85    IF (ITYP.EQ.1) THEN
        IAT=ISPE
        IAA=IAT
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
        IMM=NATMI+IML
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
      IF (IDEZ(MODCLF(KK),2,4).EQ.1) THEN
        CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .            .TRUE.)
        IF (NSTORDR >= NRAD) THEN
          IF (output)
     .    WRITE(0,*) 'MARK: IREI 01= ',irei
          DO 92 J=1,NSBOX
            TABDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK
c slmod begin - tr
            IF (IOPT1.EQ.1)
     .      TABDS1(IREI,J)=TABDS1(IREI,J)*TABDSM(IREI,J)
c slmod end
92        CONTINUE
        ELSE
          CREAC(1:9,1,KK) = CF(1:9,1)
          CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,ISP,1)=1
C     ELSEIF (IDEZ(MODCLF(KK),2,4).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,ISP,1)=2
      ELSEIF (IDEZ(MODCLF(KK),2,4).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        CALL CDEF (TEINL,1,9,KK,COUN,NSBOX,CF,.FALSE.,.FALSE.,
     .            .TRUE.)
        IF (NSTORDR >= NRAD) THEN
          DO 93 J=1,NSBOX
            TABDS1(IREI,J)=COUN(9,J)
93        CONTINUE
          DO 95 I=8,1,-1
            DO 94 J=1,NSBOX
              TABDS1(IREI,J)=TABDS1(IREI,J)*PLS(J)+COUN(I,J)
94          CONTINUE
95        CONTINUE
c slmod begin - tr
          IF (output)
     .    WRITE(0,*) 'MARK: IREI 02= ',irei
c          DO J=1,NSBOX
c            IF (IOPT1.EQ.1)
c     .      TABDS1(IREI,J)=TABDS1(IREI,J)*TABDSM(IREI,J)
c          ENDDO
c slmod end

          FCTKKL=LOG(FACTKK)
          DO 97 J=1,NSBOX
            TB=MAX(-100.D0,TABDS1(IREI,J)+FCTKKL+DEINL(J))
            TABDS1(IREI,J)=EXP(TB)
c slmod begin - tr
            IF (IOPT1.EQ.1)
     .      TABDS1(IREI,J)=TABDS1(IREI,J)*TABDSM(IREI,J)
c slmod end
97        CONTINUE
        ELSE
          CREAC(1:9,1:9,KK) = CF(1:9,1:9)
          CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
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
                IF (output)
     .          WRITE(6,*) 'MARK: EELDS1 A ',irei
                DO 101 J=1,NSBOX
                  EELDS1(IREI,J)=EELEC
101             CONTINUE
              ELSE
                NELREI(IREI)=-2
                IF (output)
     .          WRITE(6,*) 'MARK: EELDS1 B'
                EELDS1(IREI,1)=EELEC
              END IF
              MODCOL(1,4,ISP,1)=1
      ELSEIF (EFLAG.EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                IF (output)
     .          WRITE(6,*) 'MARK: EELDS1 C'
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
                MODC=IDEZ(MODCLF(KREAD),4,4)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                    IF (output)
     .              WRITE(6,*) 'MARK: EELDS1 D'
                    DO 102 J=1,NSBOX
                      EELDS1(IREI,J)=-COUN(1,J)*DEIN(J)*FACTKK/
     .                               (TABDS1(IREI,J)+EPS60)
102                 CONTINUE
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                    CREAC(1,1,KREAD)=CREAC(1,1,KREAD)+LOG(FACTKK)
                  ENDIF
                  MODCOL(1,4,ISP,1)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                         .FALSE.,.TRUE.)
                    IF (output)
     .              WRITE(6,*) 'MARK: EELDS1 E ',irei,kread
                    DO 106 J=1,NSBOX
                      EELDS1(IREI,J)=COUN(9,J)
106                 CONTINUE
                    DO 108 I=8,1,-1
                      DO 107 J=1,NSBOX
                        EELDS1(IREI,J)=EELDS1(IREI,J)*PLS(J)+COUN(I,J)
107                   CONTINUE
108                 CONTINUE
                    FCTKKL=LOG(FACTKK)
c slmod begin - not tr
c                    DO J = 1, NSBOX
c                      WRITE(6,'(A,2I5,1P,E12.4,0P,F8.4,1P,5E12.4,0P)')
c     .                  'MARK: EELDS1 B= ',
c     .                  irei,j,eelds1(irei,j),tein(j),dein(j),
c     .                  exp(eelds1(irei,j))
c                    ENDDO
c slmod end
                    DO 110 J=1,NSBOX
                      EE=MAX(-100.D0,EELDS1(IREI,J)+FCTKKL+DEINL(J))
                      EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
c slmod begin - tr
                      IF (IOPT1.EQ.1)
     .                EELDS1(IREI,J)=EELDS1(IREI,J)*TABDEM(IREI,J)

c                      IF (output)
c     .                WRITE(6,'(A,2I5,1P,E12.4,0P,F8.4,1P,5E12.4,0P)')
c     .                  'MARK: EELDS1 A= ',
c     .                  irei,j,eelds1(irei,j),tein(j),dein(j),
c     .                  eelds1(irei,j)/dein(j),exp(ee)*TABDS1(IREI,J),
c     .                  exp(ee),
c     .                  exp(ee-deinl(j))
c slmod end
110                 CONTINUE
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                    CREAC(1,1,KREAD)=CREAC(1,1,KREAD)+LOG(FACTKK)
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
        MODC=IDEZ(MODCLF(KREAD),4,4)
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
            CREAC(1,1,KREAD)=CREAC(1,1,KREAD)+LOG(FACTKK)
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
            IA=IAT
            PATDS(IREI,0)=PATDS(IREI,0)+
     +                          PATDS(IREI,IAT)
            P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                          P2ND(IREI,IA)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NATMI+IML
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
          DO 550 ISPZ=1,NSPAMI
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
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IAT),PATDS(IREI,IAT)
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
c slmod begin - f90 - not tr
c TEXT and NRC are not defined in this subroutine.
      WRITE (6,*) KK,IDEZ(MODCLF(KK),2,4)
c
c      WRITE (6,*) TEXT,NRC,KK
c slmod end
      CALL EXIT
997   CONTINUE
      WRITE (6,*) 'ERROR IN XSTEI: ISCDE FLAG'
c slmod begin - f90 - not tr
c TEXT and NRC are not defined in this subroutine.
c
c      WRITE (6,*) TEXT,NRC
c slmod end
      CALL EXIT
      END
C
      SUBROUTINE XSTCX(RMASS,IRCX,ISP,IPL,ISCD1,ISCD2,
     .                 EBULK,ISCDE,IESTM,KK,FACTKK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CGRID'
      INCLUDE 'CZT1'
      INCLUDE 'CCONA'
      DIMENSION PLS(NRAD),COUN(0:9,NRAD),CF(9,0:9)
      CHARACTER*8 TEXTS1,TEXTS2
      SAVE
c slmod begin - f90 - not tr
      IPLS = -1
      IATM = -1
c slmod end
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
      IF (IDEZ(MODCLF(KK),1,4).EQ.1) THEN
        MODCOL(3,1,ISP,IPL)=KK
        MODCOL(3,2,ISP,IPL)=3
        IF (FACTKK.NE.1.D0)
     .  WRITE (6,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTCX'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),2,4)
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
          CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
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
        IF (NSTORDR >= NRAD) THEN
          DO 251 J=1,NSBOX
            EPLCX3(IRCX,J,1)=EBULK
251       CONTINUE
        ELSE
          NELRCX(IRCX) = -2
          EPLCX3(IRCX,1,1)=EBULK
        END IF
        MODCOL(3,4,ISP,IPL)=1
      ELSEIF (NSECX4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
        DO 252 J=1,NSBOX
          EPLCX3(IRCX,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
252     CONTINUE
        ELSE
          NELRCX(IRCX) = -3
        END IF
        MODCOL(3,4,ISP,IPL)=1
      ELSEIF (NSECX4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        NELRCX(IRCX) = KREAD
        MODC=IDEZ(MODCLF(KREAD),4,4)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(3,4,ISP,IPL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          DO 253 J=1,NSBOX
            PLS(J)=TIINL(IPL,J)+ADDTL
253       CONTINUE
          CALL CDEF (PLS,1,NEND,KREAD,COUN,NSBOX,CF,.FALSE.,
     .              .FALSE.,.TRUE.)
          IF (MODC.EQ.1) THEN
            ADD=FACTKK/ADDT
            IF (NSTORDR >= NRAD) THEN
              DO 254 J=1,NSBOX
                EPLCX3(IRCX,J,1)=COUN(1,J)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSE
              EPLCX3(IRCX,1,1)=ADD
            ENDIF
          ELSEIF (MODC.EQ.2) THEN
            ADDL=LOG(FACTKK)-ADDTL
            IF (NSTORDR >= NRAD) THEN
              DO 257 J=1,NSBOX
                EPLCX3(IRCX,J,1)=COUN(1,J)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ELSE
              CREAC(1,1,KREAD) = CREAC(1,1,KREAD) + ADDL
            END IF
          ENDIF
          IF (NSTORDR >= NRAD) THEN
            DO 256 I=2,NEND
              DO 255 J=1,NSBOX
                EPLCX3(IRCX,J,I)=COUN(I,J)
255           CONTINUE
256         CONTINUE
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
      IF (ITYP1.EQ.1) TEXTS1=TEXTS(ISPZ1)
      IF (ITYP1.EQ.2) TEXTS1=TEXTS(NSPA+ISPZ1)
      IF (ITYP1.EQ.3) TEXTS1=TEXTS(NSPAM+ISPZ1)
      IF (ITYP1.EQ.4) TEXTS1=TEXTS(NSPAMI+ISPZ1)
      IF (ITYP2.EQ.1) TEXTS2=TEXTS(ISPZ2)
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
      SUBROUTINE XSTEL(IREL,ISP,IPL,
     .                 EBULK,ISCDE,IESTM,
     .                 KK,FACTKK)
C  RETURNS:
C    MODCOL(5,...)
C    TABEL3(IREL,NCELL,...)
C    EPLEL3(IREL,NCELL,...)
C    DEFEL(IREL)
C    EEFEL(IREL)
C    IESTEL(IREL,...)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'PARMMOD'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CCONA'
      INCLUDE 'CGRID'
      INCLUDE 'CZT1'
      DIMENSION PLS(NRAD),COUN(0:9,NRAD),CF(9,0:9)
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
C CROSS SECTION (E-LAB), IN FUNCTION CROSS, K=KK
      IF (IDEZ(MODCLF(KK),1,4).EQ.1) THEN
        MODCOL(5,1,ISP,IPL)=KK
        MODCOL(5,2,ISP,IPL)=3
        IF (FACTKK.NE.1.D0) THEN
          WRITE (6,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTEL'
c          GOTO 993
        ENDIF
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),2,4)
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
          CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
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
        IF (NSTORDR >= NRAD) THEN
          DO 251 J=1,NSBOX
            EPLEL3(IREL,J,1)=EBULK
251       CONTINUE
        ELSE
          NELREL(IREL) = -1
          EPLEL3(IREL,1,1)=EBULK
        END IF
        MODCOL(5,4,ISP,IPL)=1
      ELSEIF (NSEEL4.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 252 J=1,NSBOX
            EPLEL3(IREL,J,1)=1.5*TIIN(IPL,J)+EDRIFT(IPL,J)
252       CONTINUE
        ELSE
          NELREL(IREL) = -2
        END IF
        MODCOL(5,4,ISP,IPL)=1
C     ELSEIF (NSEEL4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEEL4.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        NELREL(IREL)=KREAD
        MODC=IDEZ(MODCLF(KREAD),4,4)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(5,4,ISP,IPL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          DO 253 J=1,NSBOX
            PLS(J)=TIINL(IPL,J)+ADDTL
253       CONTINUE
          CALL CDEF (PLS,1,NEND,KREAD,COUN,NSBOX,CF,.FALSE.,
     .              .FALSE.,.TRUE.)
          IF (MODC.EQ.1) THEN
            ADD=FACTKK/ADDT
            IF (NSTORDR >= NRAD) THEN
              DO 254 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSE
              EPLEL3(IREL,1,1)=ADD
            END IF
          ELSEIF (MODC.EQ.2) THEN
            ADDL=LOG(FACTKK)-ADDTL
            IF (NSTORDR >= NRAD) THEN
              DO 257 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ELSE
              CREAC(1,1,KREAD) = CREAC(1,1,KREAD) + ADDL
            END IF
          ENDIF
          IF (NSTORDR >= NRAD) THEN
            DO 256 I=2,NEND
              DO 255 J=1,NSBOX
                EPLEL3(IREL,J,I)=COUN(I,J)
255           CONTINUE
256         CONTINUE
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
      WRITE (6,*) ISP,IREL,IDEZ(MODCLF(KK),2,4)
      CALL EXIT
      END
C
      SUBROUTINE XSECTM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  TABLE FOR REACTION RATES FOR MOLECULES
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'COMUSR'
      INCLUDE 'CTRCEI'
      INCLUDE 'COMXS'
      INCLUDE 'CTEXT'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      INCLUDE 'COUTAU'
      INCLUDE 'CSPEI'
      DIMENSION PLS(NRAD),COUN(0:9,NRAD),CF(9,0:9)
      CHARACTER*8 TEXTS1,TEXTS2
C
c slmod begin - f90 - not tr
      IERR = -1
c slmdo end
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      DO 10 J=1,NSBOX
        PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
10    CONTINUE
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
          P2ND(IRDS,IATM1)=P2ND(IRDS,IATM1)+1.
          P2ND(IRDS,IATM2)=P2ND(IRDS,IATM2)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IRDS,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)+EATDS(IRDS,IATM2,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)+EATDS(IRDS,IATM2,2)
          PELDS(IRDS)=0.
          CALL CDEF (TEINL,1,1,-5,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 70 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
70          CONTINUE
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 F'
            EELDS1(IRDS,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 6 EV
            EHVDS1(IRDS,1:NSBOX)=6.
          ELSE
            CREAC(1:9,1,-5) = CF(1:9,1)
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
73        ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGMEI(IMOL,IDSC1)=IRDS
          PATDS(IRDS,IATM1)=PATDS(IRDS,IATM1)+1.
          PPLDS(IRDS,IPLS2)=PPLDS(IRDS,IPLS2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSP(IPLS2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSP(IPLS2)
          P2ND(IRDS,IATM1)=P2ND(IRDS,IATM1)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EPLDS(IRDS,      1)=RMASSP(IPLS2)/ACCMAS
          EPLDS(IRDS,      2)=1./RMASSP(IPLS2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)
          PELDS(IRDS)=1.0
          CALL CDEF (TEINL,1,1,-6,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)*FACTKK
71          CONTINUE
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 G'
            EELDS1(IRDS,1:NSBOX)=-25.0
C  TRANSFERRED KINETIC ENERGY: 10 EV
            EHVDS1(IRDS,1:NSBOX)=10.0
          ELSE
            CREAC(1:9,1,-6) = CF(1:9,1)
            CREAC(1,1,-6) = CREAC(1,1,-6) + LOG(FACTKK)
            NREAEI(IRDS) = -6
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -6
            NREAHV(IRDS) = -3
          END IF
          IF (ICOUNT.EQ.1) THEN
            IATM1=IATM2
            IPLS2=IPLS1
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
          CALL CDEF (TEINL,1,1,-7,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
72          CONTINUE
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 H'
            EELDS1(IRDS,1:NSBOX)=EELEC
          ELSE
            CREAC(1:9,1,-7) = CF(1:9,1)
            NREAEI(IRDS) = -7
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -7
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 I'
            EELDS1(IRDS,1)=EELEC
          END IF
C
76        CONTINUE
          NMDSI(IMOL)=IDSC1
C
          MODCOL(1,2,NATMI+IMOL,1)=1
          MODCOL(1,4,NATMI+IMOL,1)=1
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 90 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.1) GOTO 90
            CHRDF0=0.
            IML=NATMI+IMOL
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
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.3) GOTO 130
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IF (NRCXI.GT.NRCX) GOTO 997
            IRCX=NRCXI
            LGMCX(IMOL,IDSC2,0)=IRCX
            LGMCX(IMOL,IDSC2,1)=IPLS
            IML=NATMI+IMOL
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
        IDSC3=0
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
            FACTKK=FREACM(IMOL,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.5) GOTO 230
C  BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKM(IMOL,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 991
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
            IDSC3=IDSC3+1
            NRELI=NRELI+1
c            IF (output) THEN
c            WRITE(0,*) 'MARK: BULK PARTICLE INDEX 6'
c            WRITE(6,*) 'MARK: BULK PARTICLE INDEX 6'
c            ENDIF
            IF (NRELI.GT.NREL) GOTO 998
            IREL=NRELI
            LGMEL(IMOL,IDSC3,0)=IREL
            LGMEL(IMOL,IDSC3,1)=IPLS
C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKM(IMOL,NRC).NE.0) THEN
              IF (NPBGKM(IMOL).EQ.0) THEN
                NRBGI=NRBGI+3
c                IF (output) THEN
c                WRITE(0,*) 'MARK: BGK 6'
c                WRITE(6,*) 'MARK: BGK 6'
c                ENDIF
                IF (NRBGI.GT.NBGV) GOTO 998
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
c                WRITE(0,*) 'MARK: CROSS!'
                IF (NPBGKP(IPLS,2).EQ.0) THEN
                  NPBGKP(IPLS,2)=IBGKM(IMOL,NRC)
                ELSE
                  GOTO 999
                ENDIF
              ENDIF
            ENDIF
            RMASS=RMASSM(IMOL)
            IML=NATMI+IMOL
            IPL=IPLS
            ISCDE=ISCDEM(IMOL,NRC)
            IESTM=IESTMM(IMOL,NRC)
            EBULK=EBULKM(IMOL,NRC)
            CALL XSTEL(IREL,IML,IPL,EBULK,
     .                 ISCDE,IESTM,KK,FACTKK)
C
230       CONTINUE
C
          NMELI(IMOL)=IDSC3
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
            WRITE (6,*) 'NO ELECTRON IMPACT DISSOCIATION '
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
      WRITE (6,*) 'I.E. TOO MANY DISS. REACTIONS REQUIRED. '
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (6,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (6,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (6,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      WRITE (6,*) 'ERROR FLAG = ',IERR
      CALL EXIT
997   CONTINUE
      WRITE (6,*) 'INSUFFICIENT STORAGE FOR CX: NRCX=',NRCX
      CALL EXIT
998   CONTINUE
      WRITE (6,*) 'INSUFFICIENT STORAGE FOR EL: NREL,NRBGI=',NREL,NRBGI
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
      SUBROUTINE XSECTI
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  TABLE FOR REACTION RATES FOR TEST IONS
      INCLUDE 'PARMMOD'
      INCLUDE 'CGRID'
      INCLUDE 'COMUSR'
      INCLUDE 'CTRCEI'
      INCLUDE 'COMXS'
      INCLUDE 'CTEXT'
      INCLUDE 'CSPEI'
      INCLUDE 'CCONA'
      INCLUDE 'CZT1'
      DIMENSION COUN(0:9,NRAD),PLS(NRAD),CF(9,0:9)
      CHARACTER*8 TEXTS1,TEXTS2
c slmod begin - f90 - not tr
      IERR = -1
c slmdo end
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      DO 10 J=1,NSBOX
        PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
10    CONTINUE
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
7000      ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IF (NREII.GT.NRDS) GOTO 993
          IRDS=NREII
          LGIEI(IION,IDSC1)=IRDS
          PATDS(IRDS,IATM1)=PATDS(IRDS,IATM1)+1.
          PPLDS(IRDS,IPLS2)=PPLDS(IRDS,IPLS2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSP(IPLS2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSP(IPLS2)
          P2ND(IRDS,IATM1)=P2ND(IRDS,IATM1)+1.
          P2ND(IRDS,IATM2)=P2ND(IRDS,IATM2)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EPLDS(IRDS,      1)=RMASSP(IPLS2)/ACCMAS
          EPLDS(IRDS,      2)=1./RMASSP(IPLS2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)
          PELDS(IRDS)=0.
          CALL CDEF (TEINL,1,1,-8,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 73 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)*FACTKK
73          CONTINUE
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 J'
            EELDS1(IRDS,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 8.6 EV
            EHVDS1(IRDS,1:NSBOX)=8.6
          ELSE
            CREAC(1:9,1,-8) = CF(1:9,1)
            CREAC(1,1,-8) = CREAC(1,1,-8) + LOG(FACTKK)
            NREAEI(IRDS) = -8
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -8
            NREAHV(IRDS) = -4
          END IF
          IF (ICOUNT.EQ.1) THEN
            IATM1=IATM2
            IPLS2=IPLS1
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
          CALL CDEF (TEINL,1,1,-9,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
71          CONTINUE
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 K'
            EELDS1(IRDS,1:NSBOX)=-15.5
C  TRANSFERRED KINETIC ENERGY: 0.5 EV
            EHVDS1(IRDS,1:NSBOX)=0.5
          ELSE
            CREAC(1:9,1,-9) = CF(1:9,1)
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
          P2ND(IRDS,IATM1)=P2ND(IRDS,IATM1)+1.
          P2ND(IRDS,IATM2)=P2ND(IRDS,IATM2)+1.
          EATDS(IRDS,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IRDS,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IRDS,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IRDS,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IRDS,0,    1)=EATDS(IRDS,IATM1,1)+EATDS(IRDS,IATM2,1)
          EATDS(IRDS,0,    2)=EATDS(IRDS,IATM1,2)+EATDS(IRDS,IATM2,2)
          PELDS(IRDS)=-1.
          CALL CDEF (TEINL,1,1,-10,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .               .TRUE.)
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              TABDS1(IRDS,J)=COUN(1,J)*DEIN(J)
72          CONTINUE
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            IF (output)
     .      WRITE(6,*) 'MARK: EELDS1 L'
            EELDS1(IRDS,1:NSBOX)=-0.88*TEIN(1:NSBOX)
C  TRANSFERRED KINETIC ENERGY: INGOING ELECTRON ENERGY
            EHVDS1(IRDS,1:NSBOX)=0.88*TEIN(1:NSBOX)
          ELSE
            CREAC(1:9,1,-10) = CF(1:9,1)
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
          MODCOL(1,2,NATMI+NMOLI+IION,1)=1
          MODCOL(1,4,NATMI+NMOLI+IION,1)=1
C
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.1) GOTO 90
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
            IF (output)
     .      WRITE(0,*) 'MARK: CALLING XSTEI 02'
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
            FACTKK=FREACI(IION,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (ISWR(KK).NE.3) GOTO 130
            IPLS=IDEZ(IBULKI(IION,NRC),3,3)
            IDSC2=IDSC2+1
            NRCXI=NRCXI+1
            IF (NRCXI.GT.NRCX) GOTO 997
            IRCX=NRCXI
            LGICX(IION,IDSC2,0)=IRCX
            LGICX(IION,IDSC2,1)=IPLS
            IIO=NATMI+NMOLI+IION
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
      DO 1000 IION=1,NIONI
C
        DO 500 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          CALL XSTEI_1(IRDS)
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
              IRDS=LGIEI(IION,IIDS)
              CALL XSTEI_2(IRDS)
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
      WRITE (6,*) 'I.E. TOO MANY DISS. REACTIONS REQUIRED. '
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (6,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (6,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      WRITE (6,*) 'ERROR FLAG = ',IERR
      CALL EXIT
997   CONTINUE
      WRITE (6,*) 'INSUFFICIENT STORAGE FOR CX: NRCX=',NRCX
      CALL EXIT
998   CONTINUE
      WRITE (6,*) 'INSUFFICIENT STORAGE FOR EL: NREL=',NREL
      CALL EXIT
      RETURN
      END
C
      SUBROUTINE XSECTP
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR BULK ION SPECIES
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CTRCEI'
      INCLUDE 'CZT1'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CCONA'
      INCLUDE 'CGRID'
      INCLUDE 'CTEXT'
      INCLUDE 'CSPEI'
C
      DIMENSION PLS(NRAD),COUN(0:9,NRAD),CF(9,0:9)
      SAVE
c slmod begin - not tr (yet)
      COMMON /MULCOM/ IOPT1
c slmod end
C
      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      DO 70 J=1,NSBOX
        PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
70    CONTINUE
C
C   RECOMBINATION
C
      DO 1000 IPLS=1,NPLSI
        IF (output)
     .  WRITE(0,*) 'MARK: XSECTP: RECOMBINATION, IPLS= ',IPLS,NRCP(IPLS)
C
        IDSC=0
        LGPRC(IPLS,0)=0
C
        IF (NRCP(IPLS).EQ.0) THEN
c          IF (output)
c     .    WRITE(0,*) 'MARK: XSECTP: NRCP = 0'
C
          IF (NCHARP(IPLS).EQ.1.AND.NCHRGP(IPLS).EQ.1) THEN
C
C  DEFAULT HYDROGENIC RECOMBINATION MODEL
C  HYDR. RECOMBINATION RATE-COEFFICIENT (1/S/CCM) E + H+ --> H + RAD.
C  GORDEEV ET. AL., PIS'MA ZH. EHKSP. TEOR. FIZ. 25 (1977) 223.
C
            IF (output)
     .      WRITE(0,*) 'MARK: XSECTP: DEFAULT MODEL'
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
                IF (output) THEN
                WRITE(0,*) 'MARK: EELRC1 A'
                WRITE(6,*) 'MARK: EELRC1 A'
                ENDIF
                DO 51 J=1,NSBOX
                  ZX=EIONH/MAX(1.D-5,TEIN(J))
                  TABRC1(IRRC,J)=1.27E-13*ZX**1.5/(ZX+0.59)*DEIN(J)
c slmod begin - tr
                  IF (IOPT1.EQ.1)
     .            TABRC1(IRRC,J)=TABRC1(IRRC,J)*TABRCM(IRRC,J)
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
            MODCOL(6,2,NATMI+NMOLI+NIONI+IPLS,1)=1
            MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=1
C
            NPRCI(IPLS)=IDSC
          ENDIF
C
C  NON DEFAULT MODEL:  240--
C
        ELSEIF (NRCP(IPLS).GT.0) THEN
          IF (output)
     .    WRITE(0,*) 'MARK: XSECTP: NRCP > 0'
          DO 82 NRC=1,NRCP(IPLS)
            KK=IREACP(IPLS,NRC)
            FACTKK=FREACP(IPLS,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
C  RECOMBINATION MODEL FOR BULK IONS
            IF (ISWR(KK).EQ.6) THEN
              IF (output)
     .        WRITE(0,*) 'MARK: XSECTP: RECOMMENDED MODEL'
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
              ELSEIF (ITYP.EQ.4) THEN
                NPLPRC(IRRC)=ISPZ
              ELSEIF (ITYP.EQ.1) THEN
                NATPRC(IRRC)=ISPZ
              ELSEIF (ITYP.EQ.2) THEN
                NMLPRC(IRRC)=ISPZ
              ENDIF
C
C  1.) CROSS SECTION(TE)
C             NOT NEEDED
C  2.  RATE COEFFICIENT (CM**3/S) * DENSITY (CM**-3)
C
C  2.A) RATE COEFFICIENT = CONST.
C             TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
              IF (output)
     .        WRITE(0,*) 'MARK: XSECTP: RATE COEFFICIENT, IDEZ =',
     .                   IDEZ(MODCLF(KK),2,4)
              IF (IDEZ(MODCLF(KK),2,4).EQ.1) THEN
                CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                     .TRUE.)
                IF (output) THEN
                WRITE(0,*) 'MARK: EELRC1 B'
                WRITE(6,*) 'MARK: EELRC1 B'
                ENDIF
                IF (NSTORDR >= NRAD) THEN
                  DO 92 J=1,NSBOX
                    TABRC1(IRRC,J)=TABRC1(IRRC,J)+
     +                                  COUN(1,J)*DEIN(J)*FACTKK
c slmod begin - tr
                    IF (IOPT1.EQ.1)
     .              TABRC1(IRRC,J)=TABRC1(IRRC,J)*TABRCM(IRRC,J)
c slmod end
92                CONTINUE
                ELSE
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 1
                  CREAC(1:9,1,KK) = CF(1:9,1)
                  CREAC(1,1,KK) = CREAC(1,1,KK) + LOG(FACTKK)
                END IF
                MODCOL(6,2,NATMI+NMOLI+NIONI+IPLS,1)=1
C             ELSEIF (IDEZ(MODCLF(KK),2,4).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM): IRRELEVANT
              ELSEIF (IDEZ(MODCLF(KK),2,4).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
                IF (output)
     .          WRITE(0,*) 'MARK: XSECTP: CALLING CDEFN'
                CALL CDEFN(TEINL,PLS,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                      .TRUE.)
                IF (NSTORDR >= NRAD) THEN
                  IF (output)
     .            WRITE(0,*) 'MARK: XSECTP: NSTORDR >= NRAD'
                  DO 93 J=1,NSBOX
                    TABRC1(IRRC,J)=COUN(1,J)*DEIN(J)*FACTKK
c slmod begin - tr
                    IF (IOPT1.EQ.1)
     .              TABRC1(IRRC,J)=TABRC1(IRRC,J)*TABRCM(IRRC,J)
c slmod end
93                CONTINUE
                ELSE
                  IF (output)
     .            WRITE(0,*) 'MARK: XSECTP: ELSE'
                  CREAC(1:9,1:9,KK) = CREAC(1:9,1:9,KK) * FACTKK
                  NREARC(IRRC) = KK
                  JEREARC(IRRC) = 2
                END IF
                IF (output)
     .          WRITE(0,*) 'MARK: XSECTP: MODCOL'
                MODCOL(6,2,NATMI+NMOLI+NIONI+IPLS,1)=1
              ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ELECTRON ENERGY LOSS RATE
C
              IF (output)
     .        WRITE(0,*) 'MARK: XSECTP: ELECTRON ENERGY LOSS RATE'
              NSERC5=IDEZ(ISCDEP(IPLS,NRC),5,5)
              IF (NSERC5.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
                IF (output) THEN
                WRITE(0,*) 'MARK: EELRC1 C'
                WRITE(6,*) 'MARK: EELRC1 C'
                ENDIF
                IF (NSTORDR >= NRAD) THEN
                  DO 101 J=1,NSBOX
                    EELRC1(IRRC,J)=EELECP(IPLS,NRC)*TABRC1(IRRC,J)
101               CONTINUE
                ELSE
                  NELRRC(IRRC) = -2
                  EELRC1(IRRC,1)=EELECP(IPLS,NRC)
                END IF
                MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=1
              ELSEIF (NSERC5.EQ.1) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ELECTRON = -1.5*TE*RATECOEFF.
                IF (output) THEN
                WRITE(0,*) 'MARK: EELRC1 D'
                WRITE(6,*) 'MARK: EELRC1 D'
                ENDIF
                IF (NSTORDR >= NRAD) THEN
                  DO 102 J=1,NSBOX
                    EELRC1(IRRC,J)=-1.5*TEIN(J)*TABRC1(IRRC,J)
102               CONTINUE
                ELSE
                  NELRRC(IRRC) = -3
                END IF
                MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=1
              ELSEIF (NSERC5.EQ.3) THEN
C  4.C)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELECP(IPLS,NRC)
                MODC=IDEZ(MODCLF(KREAD),4,4)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >= NRAD) THEN
                    IF (output) THEN
                    WRITE(0,*) 'MARK: EELRC1 E'
                    WRITE(6,*) 'MARK: EELRC1 E'
                    ENDIF
                    CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                    DO 104 J=1,NSBOX
                      EELRC1(IRRC,J)=-COUN(1,J)*DEIN(J)*FACTKK
104                 CONTINUE
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=1
                    CREAC(1,1,KREAD) = CREAC(1,1,KREAD) + LOG(FACTKK)
                  END IF
                  MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=1
C  4.D)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C               ELSEIF (MODC.EQ.2) THEN
C        IRRELEVANT
C                 MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=2
C  4.E)  ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                         .FALSE.,.TRUE.)
                    IF (output) THEN
                    WRITE(0,*) 'MARK: EELRC1 F ',kread
                    WRITE(6,*) 'MARK: EELRC1 F ',kread
                    ENDIF
                    DO 106 J=1,NSBOX
                      EELRC1(IRRC,J)=COUN(9,J)
106                 CONTINUE
                    DO 108 I=8,1,-1
                      DO 107 J=1,NSBOX
                        EELRC1(IRRC,J)=EELRC1(IRRC,J)*PLS(J)+
     .                                      COUN(I,J)
107                   CONTINUE
108                 CONTINUE
c slmod begin - tr
c                    IF (output) THEN
c                    DO J=1,NSBOX
c                      WRITE(6,'(A,2I5,5G12.4)') 'MARK: EELRC1 A= ',
c     .                  irrc,j,eelrc1(irrc,j),tein(j),dein(j),coun(1,j),
c     .                  eelrc1(irrc,j)/dein(j)
c                    ENDDO
c                    ENDIF
c slmod end
                    FCTKKL=LOG(FACTKK)
                    DO 109 J=1,NSBOX
                      EE=MAX(-100.D0,EELRC1(IRRC,J)+FCTKKL+DEINL(J))
                      EELRC1(IRRC,J)=-EXP(EE)
c slmod begin - tr
                      IF (IOPT1.EQ.1)
     .                EELRC1(IRRC,J)=EELRC1(IRRC,J)*TABREM(IRRC,J)

c                      IF (output)
c     .                WRITE(6,'(A,2I5,1P,E12.4,0P,F8.4,1P,2E12.4,0P)')
c     .                  'MARK: EELRC1 B= ',
c     .                  irrc,j,eelrc1(irrc,j),tein(j),dein(j),
c     .                  -eelrc1(irrc,j)/dein(j)
c slmod end
109                 CONTINUE
                  ELSE
                    NELRRC(IRRC)=KREAD
                    JELRRC(IRRC)=9
                    CREAC(1,1,KREAD) = CREAC(1,1,KREAD) + LOG(FACTKK)
                  END IF
                  MODCOL(6,4,NATMI+NMOLI+NIONI+IPLS,1)=1
                ENDIF
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    IF (output) THEN
                    WRITE(0,*) 'MARK: EELRC1 G',kread,delpot(kread)
                    WRITE(6,*) 'MARK: EELRC1 G'
                    ENDIF
                    DO 110 J=1,NSBOX
                      EELRC1(IRRC,J)=EELRC1(IRRC,J)+
     .                                    DELE*TABRC1(IRRC,J)
110                 CONTINUE
                  END IF
                ENDIF
              ENDIF
            ENDIF
            IF (output)
     .      WRITE(0,*) 'MARK: XSECTP: END OF LOOP'
C
82        CONTINUE
          NPRCI(IPLS)=IDSC
C
C  NO MODEL DEFINED
        ELSE
          NPRCI(IPLS)=0
        ENDIF
        IF (output)
     .  WRITE(0,*) 'MARK: XSECTP: BEYOND MODEL CODE'
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
     .                                     TEXTS(IATM3)
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
C
      END
C
      SUBROUTINE CONDENSE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  CONDENSE COLLISION KERNEL, IF SOME SECONDARIES ARE NOT FOLLOWED
C  BY EIRENE, I.E., IF NFOLA(IATM), NFOLM(IMOL), NFOLI(IION) LT 0
C  FOR SOME TEST PARTICLE SPECIES
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CTRCEI'
      INCLUDE 'CZT1'
      INCLUDE 'COMXS'
      INCLUDE 'COMUSR'
      INCLUDE 'CCONA'
      INCLUDE 'CGRID'
      INCLUDE 'CTEXT'
      INCLUDE 'CSPEI'
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
                WRITE (6,*) 'TEST ION ',TEXTS(ISP),' MUST BE CONDENSED'
              ENDIF
            ENDIF
220       CONTINUE
200     CONTINUE
20    CONTINUE
C
      RETURN
      END
C
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
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
            IND=0
            DO 5 I=1,3
              IND=IND+INDEX(ZEILE((IND+1):80),CHR(1:1))
              READ (ZEILE((IND+3):80),'(E20.12)') FPARM(IR,I,I0P1)
5           CONTINUE
            LGEMIN=.true.
            READ (29,'(A80)',END=990) ZEILE
          ENDIF
          IF (INDEX(ZEILE,CHRR).NE.0.AND.IFEXMX(IR,I0P1).EQ.0) THEN
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
c slmod begin - pgi - not tr
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
C
      SUBROUTINE REFDAT(TMM,TCC,WMM,WCC)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
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
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CREF'
      INCLUDE 'CSPEI'
      DIMENSION TMM(*),TCC(*),WMM(*),WCC(*)
      DIMENSION TML(12),TCL(12),WML(12),WCL(12)
      DIMENSION FELD(1092)
      DIMENSION HFTR0(NHD1,NHD2,NHD6),
     .          HFTR1(NHD1,NHD2,NHD3,NHD6),
     .          HFTR2(NHD1,NHD2,NHD3,NHD4,NHD6),
     .          HFTR3(NHD1,NHD2,NHD3,NHD4,NHD5,NHD6)
      EQUIVALENCE (RWK(NID2+1),HFTR0(1,1,1)),
     .            (RWK(NID2+1+NH0),HFTR1(1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1),HFTR2(1,1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1+NH2),HFTR3(1,1,1,1,1,1))
C
      IF (output)
     .WRITE(6,*) 'MARK: REFDAT: ASSIGNING NFLR 01= ',NHD6
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
c...  UPDATED MAR 1, 2004
c      TML(1)=2.
c      TCL(1)=1.
c      WML(1)=96.
c      WCL(1)=42.
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
      SUBROUTINE RDTRIM
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C  THIS SUBROUTINE READS SELECTIVELY SOME
C  REFLECTION DATA PRODUCED BY MONTE CARLO CODES
C
      INCLUDE 'PARMMOD'
      INCLUDE 'CREF'
      INCLUDE 'CSPEI'
      DIMENSION HFTR0(NHD1,NHD2,NHD6),
     .          HFTR1(NHD1,NHD2,NHD3,NHD6),
     .          HFTR2(NHD1,NHD2,NHD3,NHD4,NHD6),
     .          HFTR3(NHD1,NHD2,NHD3,NHD4,NHD5,NHD6)
      EQUIVALENCE (RWK(NID2+1),HFTR0(1,1,1)),
     .            (RWK(NID2+1+NH0),HFTR1(1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1),HFTR2(1,1,1,1,1)),
     .            (RWK(NID2+1+NH0+NH1+NH2),HFTR3(1,1,1,1,1,1))
C  LAST INDEX: NID2=NH0+NH1+NH2+NH3
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
      IF (output)
     .WRITE(6,*) 'MARK: RDTRIM: NFLR= ',NFLR
      DO 7 IFILE=1,NFLR
        IUN=20
c slmod begin - new
        OPEN (UNIT=IUN,FILE=REFFIL(IFILE),ACCESS='SEQUENTIAL',
     .        FORM='FORMATTED',ERR=99)
c slmod end
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
C     DO IFILE=1,NFLR
C       WRITE (9,'(1X,A72,///1X)') REFFIL(IFILE)
C       DO 52 I1=1,INE
C         DO 53 I2=1,INW
C           WRITE (9,'(//1X,F5.0,1X,F6.2,1X,F5.0,1X,F6.2,1P,3E10.2)')
C    .             TC(IFILE),TM(IFILE),WC(IFILE),WM(IFILE),
C    .             ENAR(I1),WIAR(I2),HFTR0(I1,I2,IFILE)
C           WRITE (9,'(1X)')
C           WRITE (9,'(1P,5E13.5)') (HFTR1(I1,I2,I3,IFILE),I3=1,INR)
C           WRITE (9,'(1X)')
C           WRITE (9,'(1P,5E13.5)') ((HFTR2(I1,I2,I3,I4,IFILE),
C    .                                I4=1,INR),I3=1,INR)
C           WRITE (9,'(1X)')
C           WRITE (9,'(1P,5E13.5)') (((HFTR3(I1,I2,I3,I4,I5,IFILE),
C    .                              I5=1,INR),I4=1,INR),I3=1,INR)
53        CONTINUE
52      CONTINUE
C       WRITE (9,'(///1X)')
C     ENDDO  ! IFILE
C
      RETURN
c slmod begin - new
 99   WRITE(0,*) 'ERROR RDTRIM: DATA FILE NOT FOUND.  HALT.'
      WRITE(0,*) IFILE,'FILE='//REFFIL(IFILE)(1:LEN_TRIM(REFFIL(IFILE)))
      DO I1 = 1, NFLR
        WRITE(0,*) 'I1=',I1,REFFIL(I1)(1:LEN_TRIM(REFFIL(I1)))
      ENDDO
      STOP 
c slmod end
      END
