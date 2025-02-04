C
C
      SUBROUTINE OUTPLA
C
C  PRINT INPUT TALLIES ONTO OUTPUT FILE IUNOUT=6
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CGRID
      USE CSPEZ
      USE CTRCEI
      USE CGEOM
      USE CTEXT
      USE COUTAU
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), ALLOCATABLE :: HELPP(:),HELPW(:)
      REAL(DP) :: TALTYP(NTALI)
      REAL(DP) :: TALAV, HELPI, TALTOT, TOTAL
      INTEGER :: IR, IP, IT, I, NBLCKA, IB, IPRV, ITAL, NXM, NYM, NZM,
     .           ITALI, K, NF, NFTI, NFTE
C
C  TYPE OF TALLY: TALTYP=0: #              (#-UNITS)
C                 TALTYP=1: #-DENSITY      (#-UNITS/CM**3)
C                 TALTYP=2: VOLUME         (CM**3)
C                 TALTYP=3: DIMENSIONLESS  (1)
C                 TALTYP=4: UNKNOWN        (?)
      TALTYP(1)=0
      TALTYP(2)=0
      TALTYP(3)=1
      TALTYP(4)=1
      TALTYP(5)=0
      TALTYP(6)=0
      TALTYP(7)=0
      TALTYP(8)=3
      TALTYP(9)=3
      TALTYP(10)=3
      TALTYP(11)=3
      TALTYP(12)=4
      TALTYP(13)=0
      TALTYP(14)=2
      TALTYP(15)=3
      TALTYP(16)=3
      TALTYP(17)=3

      IF (NVOLPR > 0) THEN
        ALLOCATE (HELPP(NRAD))
        ALLOCATE (HELPW(NRAD))
      END IF
C
C  PRINT INPUT VOLUME AVERAGED TALLIES, WHICH HAVE BEEN SELECTED
C
      NXM=MAX(1,NR1STM)
      NYM=MAX(1,NP2NDM)
      NZM=MAX(1,NT3RDM)
      DO 100 IPRV=1,NVOLPR
        ITAL=NPRTLV(IPRV)
        IF (ITAL.LT.0) THEN
          ITALI=-ITAL
          NF=NFRSTP(ITALI)
          NFTI=1
          NFTE=NFSTPI(ITALI)
          IF (NSPEZV(IPRV,1).GT.0) THEN
            NFTI=NSPEZV(IPRV,1)
            NFTE=MAX(NFTI,NSPEZV(IPRV,2))
          ENDIF
          DO 119 K=NFTI,NFTE
            IF (K.GT.NFSTPI(ITALI)) THEN
              CALL LEER(1)
              WRITE (6,*) 'SPECIES INDEX OUT OF RANGE IN SUBR. OUTPLA'
              WRITE (6,*) 'ITALI, K, ',ITALI,K
              CALL LEER(1)
              GOTO 119
            ENDIF
            SELECT CASE (ITALI)
            CASE (1)
              HELPP(1:NSBOX) = TEIN(1:NSBOX)
            CASE (2)
              HELPP(1:NSBOX) = TIIN(MPLSTI(K),1:NSBOX)
            CASE (3)
              HELPP(1:NSBOX) = DEIN(1:NSBOX)
            CASE (4)
              HELPP(1:NSBOX) = DIIN(K,1:NSBOX)
            CASE (5)
              HELPP(1:NSBOX) = VXIN(MPLSV(K),1:NSBOX)
            CASE (6)
              HELPP(1:NSBOX) = VYIN(MPLSV(K),1:NSBOX)
            CASE (7)
              HELPP(1:NSBOX) = VZIN(MPLSV(K),1:NSBOX)
            CASE (8)
              HELPP(1:NSBOX) = BXIN(1:NSBOX)
            CASE (9)
              HELPP(1:NSBOX) = BYIN(1:NSBOX)
            CASE (10)
              HELPP(1:NSBOX) = BZIN(1:NSBOX)
            CASE (11)
              HELPP(1:NSBOX) = BFIN(1:NSBOX)
            CASE (12)
              HELPP(1:NSBOX) = ADIN(K,1:NSBOX)
            CASE (13)
              HELPP(1:NSBOX) = EDRIFT(K,1:NSBOX)
            CASE (14)
              HELPP(1:NSBOX) = VOL(1:NSBOX)
            CASE (15)
              HELPP(1:NSBOX) = WGHT(K,1:NSBOX)
            CASE (16)
              HELPP(1:NSBOX) = BXPERP(1:NSBOX)
            CASE (17)
              HELPP(1:NSBOX) = BYPERP(1:NSBOX)
            CASE DEFAULT
              WRITE (6,*) ' WRONG TALLY NUMBER, ITAL = ',ITAL
              WRITE (6,*) ' NO OUTPUT PERFORMED '
              CALL LEER(1)
              GOTO 100
            END SELECT
C
            TOTAL=0.D0
            DO 121 IB=1,NBMLT
            NBLCKA=NSTRD*(IB-1)
            DO 121 IR=1,NXM
            DO 121 IP=1,NYM
            DO 121 IT=1,NZM
              I=IR + ((IP-1)+(IT-1)*NP2T3)*NR1P2 + NBLCKA
              IF (ITALI.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DEIN(I)*VOL(I)
              ELSEIF (ITALI.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.EQ.3.OR.ITALI.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
                HELPW(I)=VOL(I)
              ELSEIF (ITALI.EQ.5.OR.ITALI.EQ.6.OR.ITALI.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.GE.8.AND.ITALI.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.16.OR.ITALI.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.12.OR.ITALI.EQ.14.OR.ITALI.EQ.15) THEN
C  ADDITIONAL TALLY, CELL VOLUME ,WEIGHT FUNCTION " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.13) THEN
C  ION DRIFT ENERGY
                HELPW(I)=DIIN(K,I)*VOL(I)
              ENDIF
              TOTAL=TOTAL+HELPW(I)
121         CONTINUE
C
C  SAME LOOP AGAIN, OVER ADDITIONAL CELL REGION
            DO 122 I=NSURF+1,NSURF+NRADD
              IF (ITALI.EQ.1) THEN
C  ELECTR. TEMPERATURE: NE*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DEIN(I)*VOL(I)
              ELSEIF (ITALI.EQ.2) THEN
C  ION TEMPERTURE: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.EQ.3.OR.ITALI.EQ.4) THEN
C  PARTICLE DENSITY PROFILES: VOLUME WEIGHTED AVERAGES
                HELPW(I)=VOL(I)
              ELSEIF (ITALI.EQ.5.OR.ITALI.EQ.6.OR.ITALI.EQ.7) THEN
C  ION DRIFT VELOCITY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ELSEIF (ITALI.GE.8.AND.ITALI.LE.11) THEN
C  B-FIELD UNIT VECTOR, B-FIELD STRENGTH   " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.16.OR.ITALI.LE.17) THEN
C  B_PERP-FIELD: " 1 - WEIGHTED" AVERAGES
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.12.OR.ITALI.EQ.14.OR.ITALI.EQ.15) THEN
C  ADDITIONAL TALLY (NO.12)
C  CELL VOLUME      (NO.14)
C  WEIGHT FUNCTION  (NO.15)
                HELPW(I)=1.D0
                IF (NSTGRD(I).GT.0) HELPW(I)=0.D0
              ELSEIF (ITALI.EQ.13) THEN
C  ION DRIFT ENERGY: NI(K)*VOLUME WEIGHTED AVERAGES
                HELPW(I)=DIIN(K,I)*VOL(I)
              ENDIF
              TOTAL=TOTAL+HELPW(I)
122         CONTINUE
C
            IF (ITALI.NE.NTALO) THEN
              CALL INTTAL (HELPP,HELPW,1,1,NSBOX,HELPI,
     .                     NR1ST,NP2ND,NT3RD,NBMLT)
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL INTVOL (HELPP,      1,1,NSBOX,HELPI,
     .                     NR1ST,NP2ND,NT3RD,NBMLT)
            ENDIF
            TALTOT=HELPI
            TALAV=HELPI/(TOTAL+EPS60)
C
            IF (TALTOT.EQ.0.D0) GOTO 118
C
            IF (ITALI.NE.NTALO) THEN
              CALL PRTTAL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
            ELSEIF (ITALI.EQ.NTALO) THEN
              CALL PRTVOL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),
     .                    TXTPUN(K,ITALI),
     .                    HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,
     .                    NFLAGV(IPRV),NTLVFL(IPRV))
            ENDIF
            CALL LEER(2)
            IF (TALTYP(ITALI).EQ.0) THEN
              CALL MASAGE
     .        ('WEIGHTED MEAN VALUE ("UNITS")               ')
              CALL MASR1 ('MEAN:   ',TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.1) THEN
              CALL MASAGE
     .        ('WEIGHTED TOTAL ("UNITS*CM**3"), MEAN ("UNITS")')
              CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.2) THEN
              CALL MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS")                         ')
              CALL MASR1 ('TOTAL:  ',TALTOT)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.3) THEN
              CALL MASAGE
     .        ('ARITHMETIC MEAN ("UNITS")                          ')
              CALL MASR1 ('MEAN:   ',TALAV)
              CALL LEER(3)
            ELSEIF (TALTYP(ITALI).EQ.4) THEN
              CALL MASAGE
     .        ('ARITHMETIC TOTAL ("UNITS"), MEAN ("UNITS")')
              CALL MASR2 ('TOTAL, MEAN:    ',TALTOT,TALAV)
              CALL LEER(3)
            ENDIF
            GOTO 119
C
118         CONTINUE
            CALL PRTTAL(TXTPLS(K,ITALI),TXTPSP(K,ITALI),TXTPUN(K,ITALI),
     .                  HELPP,NR1ST,NP2ND,NT3RD,NBMLT,NSBOX,-1,0)
            CALL MASAGE('IDENTICAL ZERO, NOT PRINTED                  ')
            CALL LEER(2)
119       CONTINUE
        ENDIF
100   CONTINUE
      CALL LEER(2)

      IF (NVOLPR > 0) THEN
        DEALLOCATE (HELPP)
        DEALLOCATE (HELPW)
      END IF
C
      RETURN
      END
