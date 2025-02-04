C
C
      SUBROUTINE XSTEL(IREL,ISP,IPL,
     .                 EBULK,ISCDE,IESTM,
     .                 KK,FACTKK)
C
C  RETURNS:
C    MODCOL(5,...)
C    TABEL3(IREL,NCELL,...)
C    EPLEL3(IREL,NCELL,...)
C    DEFEL(IREL)
C    EEFEL(IREL)
C    IESTEL(IREL,...)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE CZT1
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: EBULK, FACTKK
      INTEGER, INTENT(IN) :: IREL, ISP, IPL, ISCDE, IESTM, KK
      REAL(DP) :: PLS(NSTORDR), COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: FCTKKL, ADD, ADDL, ADDT, ADDTL, PMASS, TMASS
      INTEGER :: I, NSEEL4, NEND, J, KREAD, MODC, IDEZ, IERR, IPLTI

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
      IPLTI = MPLSTI(IPL)
C
C POTENTIAL
      IF (IDEZ(MODCLF(KK),1,5).EQ.1) THEN
        MODCOL(5,0,ISP,IPL)=KK
      ENDIF
C
C CROSS SECTION (E-LAB), IN FUNCTION CROSS, K=KK
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(5,1,ISP,IPL)=KK
        MODCOL(5,2,ISP,IPL)=3
        IF (FACTKK.NE.1.D0)
     .    WRITE (6,*) 'FREAC NOT READY FOR CROSS SECTION IN XSTEL'
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
        MODCOL(5,2,ISP,IPL)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN
          DO 242 J=1,NSBOX
            PLS(J)=TIINL(IPLTI,J)+ADDTL
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
          FACREA(KK) = LOG(FACTKK)
        END IF
      ELSE
C  NO RATE COEFFICIENT. IS THERE A CROSS SECTION AT LEAST?
        IF (MODCOL(5,2,ISP,IPL).NE.3) GOTO 993
      ENDIF
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
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 251 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK
251         CONTINUE
            NELREL(IREL) = -1
          ELSE
            NELREL(IREL) = -1
            EPLEL3(IREL,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(5,4,ISP,IPL)=3
      ELSEIF (NSEEL4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLEL3(IREL,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELREL(IREL) = -3
          ELSE
            NELREL(IREL) = -3
          END IF
        ELSE
          WRITE (6,*) 'WARNING FROM SUBR. XSTEL '
          WRITE (6,*) 'MODIFIED TREATMENT OF ELASTIC COLLISIONS '
          WRITE (6,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (6,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLEL3(IREL,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELREL(IREL) = -2
          ELSE
            NELREL(IREL) = -2
            EPLEL3(IREL,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(5,4,ISP,IPL)=1
C     ELSEIF (NSEEL4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEEL4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (6,*) 'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (6,*) 'IN ELASTIC COLLISION IREL= ',IREL
            WRITE (6,*) 'BECAUSE NO ENERGY WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTEL(IREL,3)=1
          MODCOL(5,4,ISP,IPL)=2
        ELSE  ! ION ENERGY AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELREL(IREL) = KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(5,4,ISP,IPL)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPLTI,J)+ADDTL
253         CONTINUE
            CALL CDEF (PLS,1,NEND,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                 .FALSE.,.TRUE.)
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              DO 257 J=1,NSBOX
                EPLEL3(IREL,J,1)=COUN(1,J)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF
            DO 256 I=2,NEND
              DO 255 J=1,NSBOX
                EPLEL3(IREL,J,I)=COUN(I,J)
255           CONTINUE
256         CONTINUE
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLEL3(IREL,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        IERR=5
        GOTO 993
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEL(IREL,1)=IDEZ(IESTM,1,3)
      IESTEL(IREL,2)=IDEZ(IESTM,2,3)
      IF (IESTEL(IREL,3).EQ.0) IESTEL(IREL,3)=IDEZ(IESTM,3,3)
C
C
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
      WRITE (6,*) ISP,IREL
      CALL EXIT_OWN(1)
      END
