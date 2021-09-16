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
     .           IDSC2, IRCX
      INTEGER, EXTERNAL :: IDEZ
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
            NREAEI(IRDS) = -8
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -8
            NREAHV(IRDS) = -4
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
            NREAEI(IRDS) = -9
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -9
            NREAHV(IRDS) = -5
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
            NREAEI(IRDS) = -10
            JEREAEI(IRDS) = 1
            NELREI(IRDS) = -10
            NREAHV(IRDS) = -6
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
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (6,*) 'KK,IION,IPLS ',KK,IION,IPLS
994   CONTINUE
      WRITE (6,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (6,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (6,*) 'IION = ',IION
      WRITE (6,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (6,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (6,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (6,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (6,*) 'OR OTHER COLLISION DATA INCONSISTENY '
!pb      WRITE (6,*) 'ERROR FLAG = ',IERR
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (6,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IION,IREL '
      WRITE (6,*) IION,IREL,IPLS
      CALL EXIT_OWN(1)
      RETURN
      END
