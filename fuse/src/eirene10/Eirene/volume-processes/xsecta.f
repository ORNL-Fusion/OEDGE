c 24.11.05: chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
C  6.12.05: comments changed: default cx only for H on p. No He default cx
C  2.5.06:  default resonant cx added for He on He+ and He on He++
C           also modified: cross.f, xsecta_param.f
! 30.08.06: array PLS and COUN changed to allocatable arrays
! 30.08.06: data structure for reaction data redefined          
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
! 02.03.07: remove escd2* arrays
! 22.03.07: PI reactions revised
! 25.03.07: 3rd and 4th secondary introduced
C
      SUBROUTINE XSECTA
C
C       SET UP TABLES (E.G. OF REACTION RATE ) FOR ATOMIC SPECIES
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
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

      REAL(DP) :: CF(9,0:9)
      REAL(DP), ALLOCATABLE :: PLS(:)
      REAL(DP) :: FACTKK, CHRDF0, EELEC, RMASS, DSUB, DEIMIN, EHEAVY,
     .            EBULK, TMASS, PMASS, COU, RATE_COEFF, ERATE
      INTEGER :: II, IML, IM, IIO, NTE, ISTORE, ISCND, ISCDE, IFRST,
     .           IAT, IREI, IATM, IDSC1, J, IPLS1, IPLS, IION1, NRC,
     .           KK, ISPZB, IAEL, ITYPB, IREL, IBGK, IA, ISP, IP,
     .           IAPI, IRPI, IACX, IDSC, IPL, IAEI, IESTM, IRCX, IPLSTI,
     .           ISTORE_MDCL, ITHRD, IFRTH
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2

C
C
C   ELECTRON IMPACT COLLISIONS:
C
C  FIND SPECIES INDEX OF ION AFTER IONIZATION EVENT FOR THE DEFAULT
C  ELECTRON IMPACT IONIZATION MODELS FROM INPUT MASS AND
C  AND CHARGE NUMBER
C
C
      ALLOCATE (PLS(NSTORDR))

      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 70 J=1,NSBOX
          PLS(J)=MAX(DEIMIN,DEINL(J))
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
c  hydrogenic atoms
            ISTORE=-4
            EELEC=-EIONH
          ELSEIF (NCHARA(IATM).EQ.2) THEN
c  helium atoms
            ISTORE=-1
            EELEC=-EIONHE
          ENDIF
C
          IF (NSTORDR >= NRAD) THEN
            DO 80 J=1,NSBOX
              COU = RATE_COEFF(ISTORE,TEINL(J),0._DP,.TRUE.,0,ERATE)
              TABDS1(IREI,J)=COU*DEIN(J)
80          CONTINUE
C  NO RADIATION LOSS INCLUDED
            EELDS1(IREI,1:NSBOX)=EELEC
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
          ELSE
            NREAEI(IREI) = ISTORE
            JEREAEI(IREI) = 1
            NELREI(IREI) = ISTORE
          ENDIF
          FACREA(ISTORE,1) = 1._DP
          FACREA(ISTORE,2) = 0._DP

          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
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
            ITHRD=ISCD3A(IATM,NRC)
            IFRTH=ISCD4A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EHEAVY=ESCD1A(IATM,NRC)
            EELEC=EELECA(IATM,NRC)
            IDSC1=IDSC1+1
            NREII=NREII+1
            IREI=NREII
            LGAEI(IATM,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IAT,
     .                 IFRST,ISCND,ITHRD,IFRTH,EHEAVY,CHRDF0,
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
C  TENTATIVELY ASSUME: NO CHARGE EXCHANGE BETWEEN IATM AND ANY IPLS
      DO 200 IATM=1,NATMI
        IDSC=0
        LGACX(IATM,0,0)=0
        LGACX(IATM,0,1)=0
C
C   DEFAULT MODEL 100 --- 129: RESONANT CX FOR H  + P,
C                 130 --- 139: RESONANT CX FOR HE + HE+,
C                 140 --- 149: RESONANT CX FOR HE + HE++,
C
        IF (NRCA(IATM).EQ.0) THEN
          DO 155 IPLS=1,NPLSI
C CHECK: "ATOMIC" COLLISION PARTNERS ONLY
            IF (NPRT(NSPAMI+IPLS).NE.1.OR.NPRT(NSPH+IATM).NE.1) GOTO 155
C
            IF (NCHARA(IATM).EQ.1.AND.NCHARP(IPLS).EQ.1.AND.
     .          NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HYDROGENIC PARTICLE WITH HYDROGENIC ION
C
C  FIND BULK SECONDARIES
              DO 121 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 123
121           CONTINUE
              GOTO 155
123           DO 124 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 125
124           CONTINUE
              GOTO 155
C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
125           CONTINUE
C  PROJECTILE MASS IS 1.
C  TARGET     MASS IS 1.
              PMASS=1.*PMASSA
              TMASS=1.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-1
              ISTORE_MDCL = -1
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.1) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE+ ION
C
C  FIND BULK SECONDARIES
              DO 131 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.1)
     .          GOTO 133
131           CONTINUE
              GOTO 155
133           DO 134 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 135
134           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
135           CONTINUE
C  PROJECTILE MASS IS 4.
C  TARGET     MASS IS 4.
              PMASS=4.*PMASSA
              TMASS=4.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-2
              ISTORE_MDCL = -2
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSEIF (NCHARA(IATM).EQ.2.AND.NCHARP(IPLS).EQ.2.AND.
     .              NCHRGP(IPLS).EQ.2) THEN
C  NEUTRAL HELIUM PARTICLE WITH HE++ ION
C
C  FIND BULK SECONDARIES
              DO 141 IPL=1,NPLSI
                IF (NMASSA(IATM).EQ.NMASSP(IPL).AND.NCHRGP(IPL).EQ.2)
     .          GOTO 143
141           CONTINUE
              GOTO 155
143           DO 144 IAT=1,NATMI
                IF (NMASSA(IAT).EQ.NMASSP(IPLS))
     .          GOTO 145
144           CONTINUE
              GOTO 155

C  CHARGE EXCHANGE BETWEEN IATM AND IPLS RESULTS IN IPL AND IAT
145           CONTINUE
C  PROJECTILE MASS IS 4.
C  TARGET     MASS IS 4.
              PMASS=4.*PMASSA
              TMASS=4.*PMASSA
C
C  CROSS SECTION (E-LAB): IN FUNCTION CROSS, K=-3
              ISTORE_MDCL = -3
C
C             TABCX3(IRCX,...)= NOT AVAILABLE FOR DEFAULT MODEL
C
            ELSE
              GOTO 155
            ENDIF


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
            MODCOL(3,1,IRCX)=ISTORE_MDCL

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
              DO 150 J=1,NSBOX
                EPLCX3(IRCX,J,1)=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
150           CONTINUE
              NELRCX(IRCX) = -1
            ELSE
              NELRCX(IRCX) = -1
            END IF
C
            MODCOL(3,2,IRCX)=3
            MODCOL(3,4,IRCX)=3
C
155       CONTINUE
C
          NACXI(IATM)=IDSC
C
C  NON DEFAULT CX MODEL:
C
        ELSEIF (NRCA(IATM).GT.0) THEN
          DO 160 NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.3) GOTO 160
C
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            CHRDF0=0.D0
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
     .                 IFRST,ISCND,EBULK,CHRDF0,
     .                 ISCDE,IESTM,KK,FACTKK)
C
160       CONTINUE
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

      DO IATM=1,NATMI
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
          DO NRC=1,NRCA(IATM)
            KK=IREACA(IATM,NRC)
            IF (ISWR(KK).NE.4) CYCLE
            FACTKK=FREACA(IATM,NRC)
            IF (FACTKK.EQ.0.D0) FACTKK=1.
            IF (MASSP(KK).LE.0.OR.MASST(KK).LE.0) GOTO 992
C  INCIDENT BULK PARTICLE INDEX
            IPLS=IDEZ(IBULKA(IATM,NRC),3,3)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) GOTO 990
            IDSC=IDSC+1
            NRPII=NRPII+1
            IRPI=NRPII
            NREAPI(IRPI) = KK
            LGAPI(IATM,IDSC,0)=IRPI
            LGAPI(IATM,IDSC,1)=IPLS

            IPL=IPLS
            RMASS=RMASSA(IATM)
            IFRST=ISCD1A(IATM,NRC)
            ISCND=ISCD2A(IATM,NRC)
            ITHRD=ISCD3A(IATM,NRC)
            IFRTH=ISCD4A(IATM,NRC)
            ISCDE=ISCDEA(IATM,NRC)
            IESTM=IESTMA(IATM,NRC)
            EBULK=EBULKA(IATM,NRC)
            EHEAVY=ESCD1A(IATM,NRC)
            CALL XSTPI (RMASS,IRPI,IPL,EBULK,EHEAVY,
     .                  IFRST,ISCND,ITHRD,IFRTH,
     .                  ISCDE,IESTM,KK,FACTKK)
          END DO
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
        DO IAPI=1,NAPII(IATM)
          LGAPI(IATM,0,0)=LGAPI(IATM,0,0)+LGAPI(IATM,IAPI,0)
        END DO

        DO IAPI=1,NAPII(IATM)
          IRPI=LGAPI(IATM,IAPI,0)
          CALL XSTPI_1(IRPI)
        END DO
      END DO
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
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS '
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
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
            CALL LEER(1)
          ELSE
            DO 890 IACX=1,NACXI(IATM)
              IRCX=LGACX(IATM,IACX,0)
              IPL =LGACX(IATM,IACX,1)
              CALL XSTCX_2(IRCX,IPL)
890         CONTINUE
          ENDIF

          CALL LEER(2)
          IF (LGAEL(IATM,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
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
            WRITE (iunout,*) 'NO GENERAL ION IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 885 IAPI=1,NAPII(IATM)
              IRPI=LGAPI(IATM,IAPI,0)
              IPL =LGAPI(IATM,IAPI,1)
              CALL XSTPI_2(IRPI,IPL)
885         CONTINUE
          ENDIF

        ENDIF
1000  CONTINUE

      DEALLOCATE (PLS)
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTA: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IATM,IPLS ',KK,IATM,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTA.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IATM = ',IATM
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IATM,IREL '
      WRITE (iunout,*) IATM,IREL,IPLS
      CALL EXIT_OWN(1)
      END
