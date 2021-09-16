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
     .            EBULK, TMASS, PMASS
      INTEGER :: II, IML, IM, IIO, NTE, ISTORE, ISCND, ISCDE, IFRST,
     .           IAT, IREI, IATM, IDSC1, J, IPLS1, IPLS, IION1, NRC,
     .           KK, ISPZB, IAEL, ITYPB, IREL, IBGK, IA, ISP, IP,
     .           IAPI, IRPI, IACX, IDSC, IPL, IAEI, IESTM, IRCX, IPLSTI
      INTEGER, EXTERNAL :: IDEZ
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
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
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
                IPLSTI=MPLSTI(IPLS)
                DO 127 J=1,NSBOX
                  EPLCX3(IRCX,J,1)=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
127             CONTINUE
                NELRCX(IRCX) = -1
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
          CALL MASBOX ('ATOMIC SPECIES IATM = '//TEXTS(NSPH+IATM))
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
                ISP=NSPH+IAT
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
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTA.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IATM = ',IATM
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (6,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IATM,IREL '
      WRITE (6,*) IATM,IREL,IPLS
      CALL EXIT_OWN(1)
      END
