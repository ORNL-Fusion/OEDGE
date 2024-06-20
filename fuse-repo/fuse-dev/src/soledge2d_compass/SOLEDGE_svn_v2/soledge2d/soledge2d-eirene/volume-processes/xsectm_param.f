C  april 07:  setting of NRPI  introduced (ion impact collisions)
C
C
      SUBROUTINE EIRENE_XSECTM_PARAM
 
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_CZT1
      USE EIRMOD_CTRCEI
      USE EIRMOD_CTEXT
      USE EIRMOD_COUTAU
      USE EIRMOD_COMXS
      USE EIRMOD_CSPEI
 
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
              CALL EIRENE_LEER(2)
              WRITE (iunout,*) 'MOLECULE NO ',IMOL,
     .                         ' COULD NOT BE IDENTIFIED'
              WRITE (iunout,*) 'NO DEFAULT A&M DATA ASSIGNED'
              CALL EIRENE_LEER(2)
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
C   ION IMPACT COLLISIONS
C
      DO 1000 IMOL=1,NMOLI
C
C  NO DEFAULT MODEL
C
       IF (NRCM(IMOL).EQ.0) THEN
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCM(IMOL).GT.0) THEN
          DO 150 NRC=1,NRCM(IMOL)
            KK=IREACM(IMOL,NRC)
            IF (ISWR(KK).NE.4) GOTO 150
            NRPI=NRPI+1
C
150       CONTINUE
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
