! 30.08.06: data structure for reaction data redefined
! 12.10.06: modcol revised
! 22.11.06: flag for shift of first parameter to rate_coeff introduced
C
C
      SUBROUTINE XSTPI (RMASS,IRPI,ISP,EBULK,EHEAVY,
     .                  IFRST,ISCND,ITHRD,IFRTH,
     .                  ISCDE,IESTM,KK,FACTKK)
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
      USE COMPRT
      USE COMXS
      USE CSPEI

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: EBULK, EHEAVY, FACTKK, RMASS
      INTEGER, INTENT(IN) :: IRPI, ISP, IFRST, ISCND, ITHRD,IFRTH, 
     .                       ISCDE, KK, IESTM 
      REAL(DP) :: CF(9),CFF(9)
      REAL(DP), ALLOCATABLE :: PLS(:)
      REAL(DP) :: ADD,ADDL,FCTKKL, P2N, TMASS, ADDT, ADDTL, PMASS,
     .            CHRDIF, COU, RATE_COEFF, ACCMAS, XLFTMAS,
     .            ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMSA, ACCINA, 
     .            ACCINM, ACCMSP, ACCINV, EFLAG, FEHVPI3, FEELPI1,
     .            ENERGY_RATE_COEFF, EI, EA, EN, ERATE
      INTEGER :: NSEPI4, NSEPI5, IAPI, I, NEND, J, IO, ION, IA, INUM,
     .           IML, IM, MODC, NRC, IIO, IPLTI, IP, IAT, IPL, II, IS,
     .           ICOUNT, IAA, IMM, III, KREAD, IERR, IMIN, IMAX, IRAD
      INTEGER, EXTERNAL :: IDEZ
      SAVE
C
C   ION IMPACT COLLISIONS
C
      ALLOCATE (PLS(NSTORDR))

      IPLTI=MPLSTI(ISP)
      XLFTMAS=RMASS+RMASSP(ISP)

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

      DO ICOUNT=1,4
        IF (ICOUNT == 1) THEN
C  SECONDARY INDEX, FIRST SECONDARY
          ITYP=IDEZ(IFRST,1,3)
          INUM=IDEZ(IFRST,2,3)
          ISPZ=IDEZ(IFRST,3,3)
        ELSE IF (ICOUNT == 2) THEN
C  SECONDARY INDEX, SECOND SECONDARY
          IF (ISCND == 0) EXIT
          ITYP=IDEZ(ISCND,1,3)
          INUM=IDEZ(ISCND,2,3)
          ISPZ=IDEZ(ISCND,3,3)
        ELSE IF (ICOUNT == 3) THEN
C  SECONDARY INDEX, SECOND SECONDARY
          IF (ITHRD == 0) EXIT
          ITYP=IDEZ(ITHRD,1,3)
          INUM=IDEZ(ITHRD,2,3)
          ISPZ=IDEZ(ITHRD,3,3)
        ELSE IF (ICOUNT == 4) THEN
C  SECONDARY INDEX, SECOND SECONDARY
          IF (IFRTH == 0) EXIT
          ITYP=IDEZ(IFRTH,1,3)
          INUM=IDEZ(IFRTH,2,3)
          ISPZ=IDEZ(IFRTH,3,3)
        END IF
        IF (ITYP.EQ.1) THEN 
          IAT=ISPZ
          IAA=NSPH+IAT
          PATPI(IRPI,IAT)=PATPI(IRPI,IAT)+INUM
          P2NP(IRPI,IAA)=P2NP(IRPI,IAA)+INUM
          ACCMAS=ACCMAS+INUM*RMASSA(IAT)
          ACCMSA=ACCMSA+INUM*RMASSA(IAT)
          ACCINV=ACCINV+INUM/RMASSA(IAT)
          ACCINA=ACCINA+INUM/RMASSA(IAT)
          EATPI(IRPI,IAT,1)=RMASSA(IAT)
          EATPI(IRPI,IAT,2)=1./RMASSA(IAT)
        ELSE IF (ITYP.EQ.2) THEN 
          IML=ISPZ
          IMM=NSPA+IML
          PMLPI(IRPI,IML)=PMLPI(IRPI,IML)+INUM
          P2NP(IRPI,IMM)=P2NP(IRPI,IMM)+INUM
          ACCMAS=ACCMAS+INUM*RMASSM(IML)
          ACCMSM=ACCMSM+INUM*RMASSM(IML)
          ACCINV=ACCINV+INUM/RMASSM(IML)
          ACCINM=ACCINM+INUM/RMASSM(IML)
          EMLPI(IRPI,IML,1)=RMASSM(IML)
          EMLPI(IRPI,IML,2)=1./RMASSM(IML)
        ELSE IF (ITYP.EQ.3) THEN 
          IIO=ISPZ
          III=NSPAM+IIO
          PIOPI(IRPI,IIO)=PIOPI(IRPI,IIO)+INUM
          P2NP(IRPI,III)=P2NP(IRPI,III)+INUM
          ACCMAS=ACCMAS+INUM*RMASSI(IIO)
          ACCMSI=ACCMSI+INUM*RMASSI(IIO)
          ACCINV=ACCINV+INUM/RMASSI(IIO)
          ACCINI=ACCINI+INUM/RMASSI(IIO)
          EIOPI(IRPI,IIO,1)=RMASSI(IIO)
          EIOPI(IRPI,IIO,2)=1./RMASSI(IIO)
        ELSE IF (ITYP.EQ.4) THEN 
          IPL=ISPZ
          PPLPI(IRPI,IPL)=PPLPI(IRPI,IPL)+INUM
          ACCMAS=ACCMAS+INUM*RMASSP(IPL)
          ACCMSP=ACCMSP+INUM*RMASSP(IPL)
          ACCINV=ACCINV+INUM/RMASSP(IPL)
          ACCINP=ACCINP+INUM/RMASSP(IPL)
        END IF
      END DO

      IF (ABS(ACCMAS-XLFTMAS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTPI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT SPECIES ',TEXTS(ISP)
        WRITE (iunout,*) 
     .    'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
        CALL EXIT_OWN(1)
      ENDIF
C
      DO IAT=1,NATMI
        EATPI(IRPI,IAT,1)=EATPI(IRPI,IAT,1)/ACCMAS
        EATPI(IRPI,IAT,2)=EATPI(IRPI,IAT,2)/ACCINV
      ENDDO
      EATPI(IRPI,0,1)=ACCMSA/ACCMAS
      EATPI(IRPI,0,2)=ACCINA/ACCINV
      DO IML=1,NMOLI
        EMLPI(IRPI,IML,1)=EMLPI(IRPI,IML,1)/ACCMAS
        EMLPI(IRPI,IML,2)=EMLPI(IRPI,IML,2)/ACCINV
      ENDDO
      EMLPI(IRPI,0,1)=ACCMSM/ACCMAS
      EMLPI(IRPI,0,2)=ACCINM/ACCINV
      DO IIO=1,NIONI
        EIOPI(IRPI,IIO,1)=EIOPI(IRPI,IIO,1)/ACCMAS
        EIOPI(IRPI,IIO,2)=EIOPI(IRPI,IIO,2)/ACCINV
      ENDDO
      EIOPI(IRPI,0,1)=ACCMSI/ACCMAS
      EIOPI(IRPI,0,2)=ACCINI/ACCINV

      EPLPI(IRPI,1)=ACCMSP/ACCMAS
      EPLPI(IRPI,2)=ACCINP/ACCINV
C     
      CHRDIF=0.

      CHRDIF = CHRDIF-PPLPI(IRPI,ISP)*NCHRGP(ISP)
      DO 133 IIO=1,NIONI
        CHRDIF=CHRDIF+PIOPI(IRPI,IIO)*NCHRGI(IIO)
133   CONTINUE
      DO 134 IP=1,NPLSI
        CHRDIF=CHRDIF+PPLPI(IRPI,IP)*NCHRGP(IP)
134   CONTINUE
      PELPI(IRPI)=PELPI(IRPI)+CHRDIF
C
C
C  TARGET MASS IN <SIGMA*V> FORMULA: MAXW. BULK PARTICLE
C  (= PROJECTILE MASS IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      PMASS=MASSP(KK)*PMASSA
C  PROJECTILE MASS IN <SIGMA*V> FORMULA: MONOENERG. TEST PARTICLE
C  (= TARGET PARTICLE IN CROSS SECTION MEASUREMENT: TARGET AT REST)
      TMASS=MASST(KK)*PMASSA
C
      ADDT=PMASS/RMASSP(ISP)
      ADDTL=LOG(ADDT)
      ADDPI(IRPI,ISP) = ADDTL
C
C CROSS SECTION (E-LAB)
      IF (IDEZ(MODCLF(KK),2,5).EQ.1) THEN
        MODCOL(4,1,IRPI)=KK
C  TENTATIVLEY ASSUME: SIGMA * V_EFF MODEL FOR RATE COEFFICIENT
        MODCOL(4,2,IRPI)=3
      ENDIF
C
C RATE COEFFICIENT
      MODC=IDEZ(MODCLF(KK),3,5)
      IF (MODC.GE.1.AND.MODC.LE.2) THEN
C  RATE COEFFICIENT IS AVAILABLE FROM DATABASE
C  MODC=1:  VS. T_IPLS, AND FOR V_TEST=0
C  MODC=2:  VS. T_IPLS,E_TEST (DOUBLE PARAMETER DATASET OR FIT)
        MODCOL(4,2,IRPI)=MODC
        IF (MODC.EQ.1) NEND=1
        IF (MODC.EQ.2) NEND=NSTORDT
        IF (NSTORDR >= NRAD) THEN  !  USE PRECOMPUTED TABLES: TABPI3
          DO 142 J=1,NSBOX
            PLS(J)=TIINL(IPLTI,J)+ADDTL
142       CONTINUE
          IF (MODC.EQ.1) THEN
            DO 145 J=1,NSBOX
              IF (LGVAC(J,ISP)) CYCLE
              COU = RATE_COEFF(KK,PLS(J),0._DP,.TRUE.,0,ERATE)
              TABPI3(IRPI,J,1)=COU*DIIN(ISP,J)*FACTKK
145         CONTINUE
          ELSEIF (MODC.EQ.2) THEN
            FCTKKL=LOG(FACTKK)
            DO J=1,NSBOX
              IF (LGVAC(J,ISP)) CYCLE
              CALL PREP_RTCS (KK,3,1,NEND,PLS(J),CF)
              TABPI3(IRPI,J,1:NEND)=CF(1:NEND)
              TABPI3(IRPI,J,1)=TABPI3(IRPI,J,1)+
     .                         DIINL(ISP,J)+FCTKKL
            END DO
          END IF
        ELSE  ! DO NOT USE PRECOMPUTED TABLES. COMPUTE DURING TRACKING
        END IF
      ENDIF

      FACREA(KK,1) = FACTKK
      FACREA(KK,2) = LOG(FACTKK)
C     
      DEFPI(IRPI)=LOG(CVELI2*PMASS)
      EEFPI(IRPI)=LOG(CVELI2*TMASS)
C
C  3. BULK ION MOMENTUM LOSS RATE
C
C
C  4A. HEAVY BULK PARTICLE ENERGY LOSS RATE
C
C  SET ENERGY LOSS RATE OF IMPACTING ION

      NSEPI4=IDEZ(ISCDE,4,5)
      IF (NSEPI4.EQ.0) THEN
C  4.A)  ENERGY LOSS RATE OF IMP. BULK PARTICLE = CONST.*RATECOEFF.
C        SAMPLE COLLIDING ION FROM DRIFTING MONOENERGETIC ISOTROPIC DISTRIBUTION
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO J=1,NSBOX
              EPLPI3(IRPI,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
            ENDDO
            NELRPI(IRPI) = -3
          ELSE
            NELRPI(IRPI) = -3
          END IF
        ELSE
          IF (NSTORDR >= NRAD) THEN
            DO 151 J=1,NSBOX
              EPLPI3(IRPI,J,1)=EBULK+EDRIFT(IPL,J)
151         CONTINUE
            NELRPI(IRPI) = -2
          ELSE
            NELRPI(IRPI) = -2
            EPLPI3(IRPI,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(4,4,IRPI)=3
      ELSEIF (NSEPI4.EQ.1) THEN
C  4.B) ENERGY LOSS RATE OF IMP. ION = 1.5*TI* RATECOEFF.
C       SAMPLE COLLIDING ION FROM DRIFTING MAXWELLIAN
        IF (EBULK.LE.0.D0) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 252 J=1,NSBOX
              EPLPI3(IRPI,J,1)=1.5*TIIN(IPLTI,J)+EDRIFT(IPL,J)
252         CONTINUE
            NELRPI(IRPI) = -3
          ELSE
            NELRPI(IRPI) = -3
          END IF
        ELSE
          WRITE (iunout,*) 'WARNING FROM SUBR. XSTPI '
          WRITE (iunout,*) 'MODIFIED TREATMENT OF BULK ION IMPACT '
          WRITE (iunout,*) 'SAMPLE FROM MAXWELLIAN WITH T = ',EBULK/1.5
          WRITE (iunout,*) 'RATHER THEN WITH T = TIIN '
          CALL LEER(1)
          IF (NSTORDR >= NRAD) THEN
            DO 2511 J=1,NSBOX
              EPLPI3(IRPI,J,1)=EBULK+EDRIFT(IPL,J)
2511        CONTINUE
            NELRPI(IRPI) = -2
          ELSE
            NELRPI(IRPI) = -2
            EPLPI3(IRPI,1,1)=EBULK
          END IF
        ENDIF
        MODCOL(4,4,IRPI)=1
C     ELSEIF (NSEPI4.EQ.2) THEN
C  use i-integral expressions. to be written
      ELSEIF (NSEPI4.EQ.3) THEN
C  4.B)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
C  4.C)  ENERGY LOSS RATE OF IMP. ION = EN.WEIGHTED RATE
        KREAD=EBULK
        IF (KREAD.EQ.0) THEN
c  data for mean ion energy loss are not available
c  use collision estimator for energy balance
          IF (IDEZ(IESTM,3,3).NE.1) THEN
            WRITE (iunout,*) 
     .        'COLLISION ESTIMATOR ENFORCED FOR ION ENERGY '
            WRITE (iunout,*) 'IN PI COLLISION IRPI= ',IRPI
            WRITE (iunout,*) 'BECAUSE NO ENERGY WEIGHTED RATE AVAILABLE'
          ENDIF
          IESTPI(IRPI,3)=1
          MODCOL(4,4,IRPI)=2
        ELSE
C  ION ENERGY AVERAGED RATE AVAILABLE AS REACTION NO. "KREAD"
        NELRPI(IRPI) = KREAD
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.GE.1.AND.MODC.LE.2) THEN
          MODCOL(4,4,IRPI)=MODC
          IF (MODC.EQ.1) NEND=1
          IF (MODC.EQ.2) NEND=NSTORDT
          IF (NSTORDR >= NRAD) THEN
            PLS=0._DP
            DO 253 J=1,NSBOX
              PLS(J)=TIINL(IPLTI,J)+ADDTL
253         CONTINUE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              DO 254 J=1,NSBOX
                IF (LGVAC(J,IPL)) CYCLE
                EPLPI3(IRPI,J,1)=ENERGY_RATE_COEFF(KREAD,PLS(J),0._DP,
     .                           .FALSE.,0)*DIIN(IPL,J)*ADD
254           CONTINUE
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              DO 257 J=1,NSBOX
                IF (LGVAC(J,IPL)) CYCLE
                CALL PREP_RTCS (KREAD,5,1,NEND,PLS(J),CFF)
                EPLPI3(IRPI,J,1:NEND) = CFF(1:NEND)
                EPLPI3(IRPI,J,1) = EPLPI3(IRPI,J,1)+DIINL(IPL,J)+ADDL
257           CONTINUE
            ENDIF
          ELSE
            IF (MODC.EQ.1) THEN
              ADD=FACTKK/ADDT
              EPLPI3(IRPI,1,1)=ADD
            ELSEIF (MODC.EQ.2) THEN
              ADDL=LOG(FACTKK)-ADDTL
              FACREA(KREAD,1) = EXP(ADDL)
              FACREA(KREAD,2) = ADDL
            END IF
          END IF
        ENDIF
        ENDIF
      ELSE
        WRITE (iunout,*) 'NSEPI4 ILL DEFINED IN XSTPI '
        CALL EXIT_OWN(1)
      ENDIF
C
C  4B. BULK ELECTRON ENERGY LOSS RATE
C
C  SET NET ENERGY LOSS RATE OF ELECTRON (IF ANY INVOLVED)
      NSEPI5=IDEZ(ISCDE,5,5)
      IF (NSEPI5.EQ.0) THEN
C       MODCOL(4,4,IRPI)=1
      ELSE
        WRITE (iunout,*) 'NSEPI5 ILL DEFINED IN XSTPI '
        CALL EXIT_OWN(1)
      ENDIF
C
C  4.C: HEAVY PARTICLE ENERGY GAIN RATE
C
      EFLAG=IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.C1)  RATE = CONST.*RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVPI3(IRPI,J,1)=EHEAVY
201       CONTINUE
          NRHVPI(IRPI)=-1
        ELSE
          NRHVPI(IRPI)=-1
          EHVPI3(IRPI,1,1)=EHEAVY
        END IF
C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION
      ELSEIF (EFLAG.EQ.3) THEN
C  4.C3)  ENERGY RATE = EN.WEIGHTED RATE(TE)
        KREAD=EHEAVY
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
!pb              IF (LGVAC(J,ISP)) CYCLE
              IF (LGVAC(J,NPLS+1)) CYCLE
              EHVPI3(IRPI,J,1)=ENERGY_RATE_COEFF(KREAD,TEINL(J),0._DP,
     .                .TRUE.,0)*DEIN(J)*FACTKK/(TABPI3(IRPI,J,1)+EPS60)
202         CONTINUE
            NRHVPI(IRPI)=KREAD
          ELSE
            NRHVPI(IRPI)=KREAD
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTPI '
          CALL EXIT_OWN(1)
        ENDIF
        FACREA(KREAD,1)=FACTKK
        FACREA(KREAD,2)=LOG(FACTKK)
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTPI(IRPI,1)=IDEZ(IESTM,1,3)
      IESTPI(IRPI,2)=IDEZ(IESTM,2,3)
      IF (IESTPI(IRPI,3).EQ.0) IESTPI(IRPI,3)=IDEZ(IESTM,3,3)
C
      IF (IESTPI(IRPI,1).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTPI(IRPI,1)=0
      ENDIF
      IF (IESTPI(IRPI,2).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IRPI,2)=0
      ENDIF
      IF (IESTPI(IRPI,3).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR EN.-BALANCE '
        WRITE (iunout,*) 'IRPI = ',IRPI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IRPI,3)=0
      ENDIF

      DEALLOCATE (PLS)

      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTPI_1(IRPI)
C
C  SET TOTAL NUMBER OF SECONDARIES BY TYPE OF SECONDARY: P..PI(IRPI,0)
C  AND
C  CONVERT SECONDARY SPECIES DISTRIBUTION P2NP(IRPI)  INTO
C  CUMMULATIVE DISTRIBUTION (NOT YET NORMALIZED)

      DO 510 IAT=1,NATMI
        IA=NSPH+IAT
        PATPI(IRPI,0)=PATPI(IRPI,0)+
     +                      PATPI(IRPI,IAT)
        P2NP(IRPI,IA)=P2NP(IRPI,IA-1)+
     +                      P2NP(IRPI,IA)
510   CONTINUE
      DO 520 IML=1,NMOLI
        IM=NSPA+IML
        PMLPI(IRPI,0)=PMLPI(IRPI,0)+
     +                      PMLPI(IRPI,IML)
        P2NP(IRPI,IM)=P2NP(IRPI,IM-1)+
     +                      P2NP(IRPI,IM)
520   CONTINUE
      DO 530 ION=1,NIONI
        IO=NSPAM+ION
        PIOPI(IRPI,0)=PIOPI(IRPI,0)+
     +                      PIOPI(IRPI,ION)
        P2NP(IRPI,IO)=P2NP(IRPI,IO-1)+
     +                      P2NP(IRPI,IO)
530   CONTINUE
      DO 540 IPL=1,NPLSI
        PPLPI(IRPI,0)=PPLPI(IRPI,0)+
     +                      PPLPI(IRPI,IPL)
540   CONTINUE
C
C
      P2NPI(IRPI)=PATPI(IRPI,0)+PMLPI(IRPI,0)+
     .            PIOPI(IRPI,0)


      P2N=P2NP(IRPI,NSPAMI)
      DO 550 ISPZ=NSPH+1,NSPAMI
        IF (P2N.GT.0.D0)
     .  P2NP(IRPI,ISPZ)=P2NP(IRPI,ISPZ)/P2N
550   CONTINUE
C
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTPI_2(IRPI,ISP)

      CALL LEER(2)
      WRITE (iunout,*) 'GENERAL ION IMPACT REACTION NO. IRPI= ', IRPI
      CALL LEER(1)
      EI=1.D30
      EA=-1.D30
      imin=0
      imax=0
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,ISP)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELPI1(IRPI,IRAD)
        ELSE
          EN=FEELPI1(IRPI,IRAD)
        END IF
        if (en < ei) imin=irad
        if (en > ea) imax=irad
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE

      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (iunout,*) 'ELECTRONS: PELPI, CONSTANT ENERGY: EEL'
        WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELPI(IRPI),EI
      ELSE
        WRITE (iunout,*) 
     .    'ELECTRONS: PELPI, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (iunout,'(1X,A8,3(1PE12.4))') 'EL      ',PELPI(IRPI),EI,EA
      ENDIF
      write (iunout,*) ' imin = ', imin, ' imax = ',imax
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,ISP)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVPI3(IRPI,IRAD,1)
        ELSE
          EN=FEHVPI3(IRPI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
876   CONTINUE
      IF (PPLPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLPI '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (IPL.EQ.ISP) THEN
            WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLPI(IRPI,IPL)-1
          ELSEIF (PPLPI(IRPI,IPL).NE.0.D0) THEN
            WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLPI(IRPI,IPL)
          ENDIF
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EPLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLPI(IRPI,1),
     .                                 ' * E0 + ',EPLPI(IRPI,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EPLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLPI(IRPI,1),
     .                                 ' * E0 + ',EPLPI(IRPI,2),
     .                                 ' * EHEAVY '
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (iunout,*) 'BULK IONS: NONE'
      ENDIF
      CALL LEER(1)
C
      WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NPI(IRPI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL LEER(1)
        RETURN
      ENDIF
C
      IF (PATPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATPI '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATPI(IRPI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATPI(IRPI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATPI(IRPI,0,1),
     .                                 ' * E0 + ',EATPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATPI(IRPI,0,1),
     .                                 ' * E0 + ',EATPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLPI '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLPI(IRPI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLPI(IRPI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLPI(IRPI,0,1),
     .                                 ' * E0 + ',EMLPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLPI(IRPI,0,1),
     .                                 ' * E0 + ',EMLPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIOPI(IRPI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIOPI '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIOPI(IRPI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIOPI(IRPI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIOPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIOPI(IRPI,0,1),
     .                                 ' * E0 + ',EIOPI(IRPI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIOPI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIOPI(IRPI,0,1),
     .                                 ' * E0 + ',EIOPI(IRPI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF


      CALL LEER(1)
      IF (IESTPI(IRPI,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART.-BALANCE '
      IF (IESTPI(IRPI,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM.-BALANCE '
      IF (IESTPI(IRPI,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN.-BALANCE '
      CALL LEER(1)

      RETURN
C
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED  '
      WRITE (iunout,*) 'INVALID SPECIES INDEX FOR ION IMPACT COLLISION'
      CALL EXIT_OWN(1)
992   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: EXIT CALLED  '
      WRITE (iunout,*) 
     .  'MASS NUMBERS OF INTERACTING PARTICLES INCONSISTENT'
      CALL EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTPI: ISCDE FLAG'
      WRITE (iunout,*) IRPI
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'INSUFFICIENT STORAGE FOR PI: NRPI=',NRPI
      CALL EXIT_OWN(1)
      RETURN
C
      END
