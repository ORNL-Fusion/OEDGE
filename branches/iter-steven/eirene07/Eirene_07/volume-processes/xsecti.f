C  27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if iion is a molecular ion
c          otherwise he+ atomic ions could be confused with d2+ molecular
c          ions and then get assigned the wrong default collision model
c 24.11.05 chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
C
      SUBROUTINE XSECTI
C
C  TABLE FOR REACTION RATES FOR TEST IONS
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CTEXT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: FACTKK, CHRDF0, RMASS, ACCMAS, ACCINV, DSUB, DEIMIN,
     .          EHEAVY, EELEC, EBULK, COU, RATE_COEFF
      INTEGER :: ICOUNT, IA1, IP2, IPLS, ITEST, IIO, IION, IDSC1,
     .           NRC, J, IPLS1, IPLS2, IATM, KK, IATM1, IATM2, ITYPB,
     .           ISPZB, III, IDSC, IREL, IBGK, IIDS, IERR, IMOL, IIEL,
     .           IIEI, IREI, IESTM, IFRST, ISCND, ISCDE, IPL, IICX,
     .           IDSC2, IRCX
      INTEGER, EXTERNAL :: IDEZ
      CHARACTER(8) :: TEXTS1, TEXTS2
C
!pb      DSUB=LOG(1.D8)
      DEIMIN=LOG(1.D8)
      IF (NSTORDR >= NRAD) THEN
        DO 10 J=1,NSBOX
!pb          PLS(J)=MAX(DEIMIN,DEINL(J))-DSUB
          PLS(J)=MAX(DEIMIN,DEINL(J))
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
C   FIRST: DEAL WITH EI (ELECTRON IMPACT) COLLISIONS
C 
        IF (NRCI(IION).EQ.0.AND.NCHARI(IION).EQ.2.AND.
     .      NPRT(NSPAM+IION).GT.1) THEN
C  APPLY THE DEFAULT MODEL FOR H2+ DISSOCIATION
C  USE NPRT.GT.1 TO DISTINGUISH ATOMIC FROM MOLECULAR IONS

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
              WRITE (iunout,*) 'TEST ION NO ',IION,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
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
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATDS(IREI,IA1)=PATDS(IREI,IA1)+1.
          PPLDS(IREI,IP2)=PPLDS(IREI,IP2)+1.
          ACCMAS=ACCMAS+RMASSA(IA1)
          ACCMAS=ACCMAS+RMASSP(IP2)
          ACCINV=ACCINV+1./RMASSA(IA1)
          ACCINV=ACCINV+1./RMASSP(IP2)
          P2ND(IREI,NSPH+IA1)=P2ND(IREI,NSPH+IA1)+1.
          EATDS(IREI,IA1,1)=RMASSA(IA1)/ACCMAS
          EATDS(IREI,IA1,2)=1./RMASSA(IA1)/ACCINV
          EPLDS(IREI,      1)=RMASSP(IP2)/ACCMAS
          EPLDS(IREI,      2)=1./RMASSP(IP2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IA1,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IA1,2)
          PELDS(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 73 J=1,NSBOX
              COU = RATE_COEFF(-8,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)*FACTKK
73          CONTINUE
            EELDS1(IREI,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 8.6 EV
            EHVDS1(IREI,1:NSBOX)=8.6
            NREAEI(IREI) = -8
            JEREAEI(IREI) = 1
            NELREI(IREI) = -8
            NREAHV(IREI) = -4
          ELSE
            FACREA(-8) = FACTKK
            NREAEI(IREI) = -8
            JEREAEI(IREI) = 1
            NELREI(IREI) = -8
            NREAHV(IREI) = -4
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
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PPLDS(IREI,IPLS1)=PPLDS(IREI,IPLS1)+1.
          PPLDS(IREI,IPLS2)=PPLDS(IREI,IPLS2)+1.
          EPLDS(IREI,1)=1.0
          EPLDS(IREI,2)=1.0
          PELDS(IREI)=1.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = RATE_COEFF(-9,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
71          CONTINUE
            EELDS1(IREI,1:NSBOX)=-15.5
C  TRANSFERRED KINETIC ENERGY: 0.5 EV
            EHVDS1(IREI,1:NSBOX)=0.5
            NREAEI(IREI) = -9
            JEREAEI(IREI) = 1
            NELREI(IREI) = -9
            NREAHV(IREI) = -5
          ELSE
            FACREA(-9) = 1._DP
            NREAEI(IREI) = -9
            JEREAEI(IREI) = 1
            NELREI(IREI) = -9
            NREAHV(IREI) = -5
          END IF
C  THIRD PROCESS
          ACCMAS=0.D0
          ACCINV=0.D0
          IDSC1=IDSC1+1
          NREII=NREII+1
          IREI=NREII
          LGIEI(IION,IDSC1)=IREI
          PATDS(IREI,IATM1)=PATDS(IREI,IATM1)+1.
          PATDS(IREI,IATM2)=PATDS(IREI,IATM2)+1.
          ACCMAS=ACCMAS+RMASSA(IATM1)
          ACCMAS=ACCMAS+RMASSA(IATM2)
          ACCINV=ACCINV+1./RMASSA(IATM1)
          ACCINV=ACCINV+1./RMASSA(IATM2)
          P2ND(IREI,NSPH+IATM1)=P2ND(IREI,NSPH+IATM1)+1.
          P2ND(IREI,NSPH+IATM2)=P2ND(IREI,NSPH+IATM2)+1.
          EATDS(IREI,IATM1,1)=RMASSA(IATM1)/ACCMAS
          EATDS(IREI,IATM2,1)=RMASSA(IATM2)/ACCMAS
          EATDS(IREI,IATM1,2)=1./RMASSA(IATM1)/ACCINV
          EATDS(IREI,IATM2,2)=1./RMASSA(IATM2)/ACCINV
          EATDS(IREI,0,    1)=EATDS(IREI,IATM1,1)+EATDS(IREI,IATM2,1)
          EATDS(IREI,0,    2)=EATDS(IREI,IATM1,2)+EATDS(IREI,IATM2,2)
          PELDS(IREI)=-1.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = RATE_COEFF(-10,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
72          CONTINUE
C  FOR THE FACTOR -0.88 SEE: EIRENE MANUAL, INPUT BLOCK 4, EXAMPLES
            EELDS1(IREI,1:NSBOX)=-0.88*TEIN(1:NSBOX)
C  TRANSFERRED KINETIC ENERGY: INGOING ELECTRON ENERGY
            EHVDS1(IREI,1:NSBOX)=0.88*TEIN(1:NSBOX)
            NREAEI(IREI) = -10
            JEREAEI(IREI) = 1
            NELREI(IREI) = -10
            NREAHV(IREI) = -6
          ELSE
            FACREA(-10) = 1._DP
            NREAEI(IREI) = -10
            JEREAEI(IREI) = 1
            NELREI(IREI) = -10
            NREAHV(IREI) = -6
          END IF
C
76        CONTINUE
C
          NIDSI(IION)=IDSC1
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
            IREI=NREII
            LGIEI(IION,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IIO,
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

C   SECONDLY: DEAL WITH CX (CHARGE EXCHANGE) COLLISIONS

      DO 200 IION=1,NIONI
        IDSC2=0
        LGICX(IION,0,0)=0
        LGICX(IION,0,1)=0
C
C  THERE ARE CURRENTLY NO DEFAULT CX RATES FOR TEST IONS
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
            CHRDF0=-NCHRGI(IION)
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
            CALL XSTCX(RMASS,IRCX,IIO,IPL,IFRST,ISCND,EBULK,
     .                       CHRDF0,ISCDE,IESTM,KK,FACTKK)
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
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
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
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISION'
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
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
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
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IION,IPLS ',KK,IION,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTI.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IION = ',IION
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTI: EXIT CALLED  '
      WRITE (iunout,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (iunout,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (iunout,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IION,IREL '
      WRITE (iunout,*) IION,IREL,IPLS
      CALL EXIT_OWN(1)
      RETURN
      END
