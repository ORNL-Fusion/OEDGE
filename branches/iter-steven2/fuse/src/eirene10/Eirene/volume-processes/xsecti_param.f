C
C
      SUBROUTINE XSECTI_PARAM

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

      INTEGER :: ICOUNT, NRC, KK, IION, IPLS, ITEST, IATM, IATM1,
     .           IATM2, IPLS1, IPLS2
C
      DO 100 IION=1,NIONI
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
            ICOUNT=1
          ELSE
            ICOUNT=2
          ENDIF
C
7000      CONTINUE
          NRDS=NRDS+1
          IF (ICOUNT.EQ.1) THEN
            ICOUNT=2
            GOTO 7000
          ENDIF
C  SECOND PROCESS
          NRDS=NRDS+1
C  THIRD PROCESS
          NRDS=NRDS+1
C
76        CONTINUE
C
C  NON DEFAULT MODEL SPECIFIED IN INPUT BLOCK 4
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 90 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.1) GOTO 90
C
            NRDS=NRDS+1
90        CONTINUE
C
        ENDIF
C
100   CONTINUE

      DO 200 IION=1,NIONI
C
C  NO DEFAULT CX RATES
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT CX MODEL:
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 130 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
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
      DO 300 IION=1,NIONI
C
C  DEFAULT EL MODEL: NOT AVAILABLE
C
        IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT EL MODEL:  240--
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 230 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
            IF (ISWR(KK).NE.5) GOTO 230
C
            NREL=NREL+1

C  BGK SELF AND CROSS COLLISIONS?
            IF (IBGKI(IION,NRC).NE.0) THEN
              IF (NPBGKI(IION).EQ.0) THEN
                NBGK=NBGK+3
              ENDIF
            ENDIF
C
230       CONTINUE
C
        ENDIF
C
300   CONTINUE
C
C
C   ION IMPACT COLLISIONS
C
      DO 1000 IION=1,NIONI
C
C  NO DEFAULT MODEL
C
       IF (NRCI(IION).EQ.0) THEN
C
C  NON DEFAULT ION IMPACT MODEL:  130--190
C
        ELSEIF (NRCI(IION).GT.0) THEN
          DO 150 NRC=1,NRCI(IION)
            KK=IREACI(IION,NRC)
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
