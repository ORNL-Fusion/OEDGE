C 27.6.05  irds --> irei
c 24.11.05 use nprt(ispz) to check if imol is really a molecule.
c          otherwise He atoms could be confused with d2 molecules,
c          if they accidentally are specified in the molecule block 4b
c 24.11.05 chrdf0 in parameterlist for call to xstcx
c          (was ok already for call to xstei)
! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
!           setting of modcol corrected
C
      SUBROUTINE XSECTM
C
C  TABLE FOR REACTION RATES FOR MOLECULES
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
      USE COUTAU
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP) :: PLS(NSTORDR), CF(9,0:9)
      REAL(DP) :: ACCMAS, ACCINV, FACTKK, DEIMIN, EELEC, CHRDF0, DSUB,
     .          RMASS, EBULK, EHEAVY, COU, RATE_COEFF
      INTEGER :: ITEST, IATM, IPLS, IION, IA1, IP2, ION, ICOUNT,
     .           IION3, IDSC1, NRC, KK, J, IMOL, IPLS1, IPLS2, IPLS3,
     .           IATM1, IATM2, ITYPB, ISPZB, IMEL, IDSC, IREL, IBGK,
     .           IMDS, IERR, IMCX, ISCND, ISCDE, IESTM, IML, IFRST,
     .           IRCX, IREI, IPL, IDSC2
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
C  CHECK IF THIS REALLY IS A  MOLECULE: USE NPRT(ISPZ).GT.1?

        IF (NPRT(NSPA+IMOL).LE.1) THEN
          WRITE (IUNOUT,*) 'SEVERE INPUT ERROR DETECTED IN XSECTM: '
          WRITE (IUNOUT,*) 'IMOL= ',IMOL,' CARRIES ONLY ONE FLUX UNIT'
          WRITE (IUNOUT,*) 'EXIT CALLED FROM XSECTM '
        ENDIF

C  YES, "IMOL" IS A MOLECULE !

      
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
              WRITE (iunout,*) 'MOLECULE NO ',IMOL,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
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
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
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
          PELDS(IREI)=0.
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 70 J=1,NSBOX
              COU = RATE_COEFF(-5,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
70          CONTINUE
            EELDS1(IREI,1:NSBOX)=-10.5
C  TRANSFERRED KINETIC ENERGY: 6 EV
            EHVDS1(IREI,1:NSBOX)=6.
            NREAEI(IREI)=-5
            JEREAEI(IREI)=1
            NELREI(IREI)=-5
            NREAHV(IREI)=-2
          ELSE
            FACREA(-5) = 0._DP
            NREAEI(IREI)=-5
            JEREAEI(IREI)=1
            NELREI(IREI)=-5
            NREAHV(IREI)=-2
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
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
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
          PELDS(IREI)=1.0
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 71 J=1,NSBOX
              COU = RATE_COEFF(-6,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)*FACTKK
71          CONTINUE
            EELDS1(IREI,1:NSBOX)=-25.0
C  TRANSFERRED KINETIC ENERGY: 10 EV
            EHVDS1(IREI,1:NSBOX)=10.0
            NREAEI(IREI) = -6
            JEREAEI(IREI) = 1
            NELREI(IREI) = -6
            NREAHV(IREI) = -3
          ELSE
            FACREA(-6) = LOG(FACTKK)
            NREAEI(IREI) = -6
            JEREAEI(IREI) = 1
            NELREI(IREI) = -6
            NREAHV(IREI) = -3
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
          IREI=NREII
          LGMEI(IMOL,IDSC1)=IREI
          ION=NSPAM+IION3
          PIODS(IREI,IION3)=PIODS(IREI,IION3)+1.
          P2ND(IREI,ION)=P2ND(IREI,ION)+1.
          EIODS(IREI,IION3,1)=1.
          EIODS(IREI,IION3,2)=0.
          EIODS(IREI,0,1)=1.
          EIODS(IREI,0,2)=0.
          PELDS(IREI)=1.0
          MODCOL(1,2,IREI)=1
          MODCOL(1,4,IREI)=1
          IF (NSTORDR >= NRAD) THEN
            DO 72 J=1,NSBOX
              COU = RATE_COEFF(-7,TEINL(J),0._DP,.TRUE.,0)
              TABDS1(IREI,J)=COU*DEIN(J)
72          CONTINUE
            EELDS1(IREI,1:NSBOX)=EELEC
            NREAEI(IREI) = -7
            JEREAEI(IREI) = 1
            NELREI(IREI) = -7
          ELSE
            FACREA(-7) = 0._DP
            NREAEI(IREI) = -7
            JEREAEI(IREI) = 1
            NELREI(IREI) = -7
            EELDS1(IREI,1)=EELEC
          END IF
C
76        CONTINUE
          NMDSI(IMOL)=IDSC1
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
            IREI=NREII
            LGMEI(IMOL,IDSC1)=IREI
            CALL XSTEI(RMASS,IREI,IML,
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
            CHRDF0=0.
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
     .                 IFRST,ISCND,EBULK,CHRDF0,
     .                 ISCDE,IESTM,KK,FACTKK)
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
          IREI=LGMEI(IMOL,IMDS)
          CALL XSTEI_1(IREI)
500     CONTINUE
C
        IF (TRCAMD) THEN
          CALL MASBOX ('MOLECULAR SPECIES IMOL = '//TEXTS(NSPA+IMOL))
          CALL LEER(1)
          IF (LGMCX(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO CHARGE EXCHANGE WITH BULK IONS'
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
            WRITE (iunout,*) 'NO ELECTRON IMPACT COLLISIONS '
            CALL LEER(1)
          ELSE
            DO 870 IMDS=1,NMDSI(IMOL)
              IREI=LGMEI(IMOL,IMDS)
              CALL XSTEI_2(IREI)
870         CONTINUE
          ENDIF
C
          IF (LGMEL(IMOL,0,0).EQ.0) THEN
            CALL LEER(1)
            WRITE (iunout,*) 'NO ELASTIC COLLISIONS WITH BULK IONS'
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
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR CHARGE EXCHANGE '
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ELASTIC COLLISION '
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      WRITE (iunout,*) 'KK,IMOL,IPLS ',KK,IMOL,IPLS
994   CONTINUE
      WRITE (iunout,*) 'ERROR DETECTED IN XSECTM.'
      WRITE (iunout,*) 'REACTION NO. KK= ',KK, 'NOT READ FROM FILE '
      WRITE (iunout,*) 'IMOL = ',IMOL
      WRITE (iunout,*) 'ISWR(KK) = ',ISWR(KK)
      WRITE (iunout,*) 'EXIT CALLED      '
      CALL EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSECTM: EXIT CALLED  '
      WRITE (iunout,*) 'NO COLLISION DATA AVAILABLE FOR THE CHOICE  '
      WRITE (iunout,*) 'OF POST COLLISION SAMPLING FLAG ISCDEA'
      WRITE (iunout,*) 'OR OTHER COLLISION DATA INCONSISTENY '
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'SPECIES CONFLICT FOR BGK COLLISIONS. IMOL,IREL '
      WRITE (iunout,*) IMOL,IREL,IPLS
      CALL EXIT_OWN(1)
      RETURN
C
      END
