!pb  28.06.06: bug fix for NSTORAM=0 and MODC=1
!pb            FACREA=FACTKK instead of FACREA=log(FACTKK) as single 
!pb            polynomial fit is linear in it parameter
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  22.11.06: DELPOT introduced
C
C
      SUBROUTINE XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
c   isp:  incident test particle species identifier



C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE COMPRT, ONLY: IUNOUT
      USE CCONA
      USE CGRID
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, CHRDF0, EELEC, FACTKK
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ISCDE, IESTM, KK
      REAL(DP) :: CF(9,0:9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL, FEHVDS1, EE, TB, FEELEI1, EN,
     .          P2N, EA, EI, ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV, COU, RATE_COEFF, 
     .          ENERGY_RATE_COEFF, DELE
      INTEGER :: MODC, KREAD, IM, IA, IERR, J, IPL, I, IP, IRAD, IO,
     .           ION, ISPZ, III, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML, IMIN, IMAX
      INTEGER, EXTERNAL :: IDEZ

      ITYP=IDEZ(IFRST,1,3)
      INUM=IDEZ(IFRST,2,3)
      ISPE=IDEZ(IFRST,3,3)
C
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
C
      ICOUNT=1
 85   CONTINUE
      IF (ITYP.EQ.1) THEN
        IAT=ISPE
        IAA=NSPH+IAT
        PATDS(IREI,IAT)=PATDS(IREI,IAT)+INUM
        P2ND(IREI,IAA)=P2ND(IREI,IAA)+INUM
        ACCMAS=ACCMAS+INUM*RMASSA(IAT)
        ACCMSA=ACCMSA+INUM*RMASSA(IAT)
        ACCINV=ACCINV+INUM/RMASSA(IAT)
        ACCINA=ACCINA+INUM/RMASSA(IAT)
        EATDS(IREI,IAT,1)=RMASSA(IAT)
        EATDS(IREI,IAT,2)=1./RMASSA(IAT)
      ELSEIF (ITYP.EQ.2) THEN
        IML=ISPE
        IMM=NSPA+IML
        PMLDS(IREI,IML)=PMLDS(IREI,IML)+INUM
        P2ND(IREI,IMM)=P2ND(IREI,IMM)+INUM
        ACCMAS=ACCMAS+INUM*RMASSM(IML)
        ACCMSM=ACCMSM+INUM*RMASSM(IML)
        ACCINV=ACCINV+INUM/RMASSM(IML)
        ACCINM=ACCINM+INUM/RMASSM(IML)
        EMLDS(IREI,IML,1)=RMASSM(IML)
        EMLDS(IREI,IML,2)=1./RMASSM(IML)
      ELSEIF (ITYP.EQ.3) THEN
        IIO=ISPE
        III=NSPAM+IIO
        PIODS(IREI,IIO)=PIODS(IREI,IIO)+INUM
        P2ND(IREI,III)=P2ND(IREI,III)+INUM
        ACCMAS=ACCMAS+INUM*RMASSI(IIO)
        ACCMSI=ACCMSI+INUM*RMASSI(IIO)
        ACCINV=ACCINV+INUM/RMASSI(IIO)
        ACCINI=ACCINI+INUM/RMASSI(IIO)
        EIODS(IREI,IIO,1)=RMASSI(IIO)
        EIODS(IREI,IIO,2)=1./RMASSI(IIO)
      ELSEIF (ITYP.EQ.4) THEN
        IPL=ISPE
        PPLDS(IREI,IPL)=PPLDS(IREI,IPL)+INUM
        ACCMAS=ACCMAS+INUM*RMASSP(IPL)
        ACCMSP=ACCMSP+INUM*RMASSP(IPL)
        ACCINV=ACCINV+INUM/RMASSP(IPL)
        ACCINP=ACCINP+INUM/RMASSP(IPL)
      ENDIF
C
      IF (ISCND.NE.0.AND.ICOUNT.EQ.1) THEN
        ITYP=IDEZ(ISCND,1,3)
        INUM=IDEZ(ISCND,2,3)
        ISPE=IDEZ(ISCND,3,3)
        ICOUNT=2
        GOTO 85
      ENDIF
C
      IF (ABS(ACCMAS-RMASS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTEI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT SPECIES ',TEXTS(ISP)
        WRITE (iunout,*) 'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
        CALL EXIT_OWN(1)
      ENDIF
C
      DO IAT=1,NATMI
        EATDS(IREI,IAT,1)=EATDS(IREI,IAT,1)/ACCMAS
        EATDS(IREI,IAT,2)=EATDS(IREI,IAT,2)/ACCINV
      ENDDO
      EATDS(IREI,0,1)=ACCMSA/ACCMAS
      EATDS(IREI,0,2)=ACCINA/ACCINV
      DO IML=1,NMOLI
        EMLDS(IREI,IML,1)=EMLDS(IREI,IML,1)/ACCMAS
        EMLDS(IREI,IML,2)=EMLDS(IREI,IML,2)/ACCINV
      ENDDO
      EMLDS(IREI,0,1)=ACCMSM/ACCMAS
      EMLDS(IREI,0,2)=ACCINM/ACCINV
      DO IIO=1,NIONI
        EIODS(IREI,IIO,1)=EIODS(IREI,IIO,1)/ACCMAS
        EIODS(IREI,IIO,2)=EIODS(IREI,IIO,2)/ACCINV
      ENDDO
      EIODS(IREI,0,1)=ACCMSI/ACCMAS
      EIODS(IREI,0,2)=ACCINI/ACCINV
      EPLDS(IREI,1)=ACCMSP/ACCMAS
      EPLDS(IREI,2)=ACCINP/ACCINV
C
      CHRDIF=CHRDF0
      DO 83 IIO=1,NIONI
        CHRDIF=CHRDIF+PIODS(IREI,IIO)*NCHRGI(IIO)
83    CONTINUE
      DO 84 IPL=1,NPLSI
        CHRDIF=CHRDIF+PPLDS(IREI,IPL)*NCHRGP(IPL)
84    CONTINUE
      PELDS(IREI)=PELDS(IREI)+CHRDIF
C
C
C  1.) CROSS SECTION(TE) : NOT NEEDED
C
C
C  2.) RATE COEFFICIENT (CM**3/S) * ELECTRON DENSITY (CM**-3)
C
C
C  2.A) RATE COEFFICIENT = CONST.
C     TO BE WRITTEN
C  2.B) RATE COEFFICIENT(TE)
      IF (IDEZ(MODCLF(KK),3,5).EQ.1) THEN
        IF (NSTORDR >= NRAD) THEN
C  RATE:  (1/S)
C  RATE COEFFICIENT: (CM^3/S)
          DO J=1,NSBOX
            COU = RATE_COEFF(KK,TEINL(J),0._DP,.TRUE.,0)
            TABDS1(IREI,J)=COU*FACTKK
            IF (IFTFLG(KK,2) < 100) 
     .        TABDS1(IREI,J)=TABDS1(IREI,J)*DEIN(J)
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ELSE
          FACREA(KK) = FACTKK
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,IREI)=1
C     ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,IREI)=2
      ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)
          DO J=1,NSBOX
            COU = RATE_COEFF(KK,TEINL(J),PLS(J),.FALSE.,1)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DEINL(J)
            TB=MAX(-100._DP,TB)
            TABDS1(IREI,J)=EXP(TB)
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ELSE
          FACREA(KK) = FACTKK
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ENDIF
        MODCOL(1,2,IREI)=1
      ELSE
        IERR=1
        GOTO 996
      ENDIF
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ENERGY RATES
C
C  4.A: ELECTRON ENERGY LOSS RATES
C
      EFLAG=IDEZ(ISCDE,5,5)
      IF (EFLAG.EQ.0) THEN
C  4.A1) ENERGY LOSS RATE OF IMP. ELECTRON = CONST.*RATECOEFF.
              IF (NSTORDR >= NRAD) THEN
                DO 101 J=1,NSBOX
                  EELDS1(IREI,J)=EELEC
101             CONTINUE
                NELREI(IREI)=-2
              ELSE
                NELREI(IREI)=-2
                EELDS1(IREI,1)=EELEC
              END IF
              MODCOL(1,4,IREI)=1
      ELSEIF (EFLAG.EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 103 J=1,NSBOX
                  EELDS1(IREI,J)=-1.5*TEIN(J)
103             CONTINUE
                NELREI(IREI)=-3
              ELSE
                NELREI(IREI)=-3
              END IF
              MODCOL(1,4,IREI)=1
      ELSEIF (EFLAG.EQ.3) THEN
C  4.A2) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELEC
                MODC=IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    DO 102 J=1,NSBOX
                      EELDS1(IREI,J)=-ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .                                0._DP,.TRUE.,0)*DEIN(J)*FACTKK/
     .                                (TABDS1(IREI,J)+EPS60) 
102                 CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                    FACREA(KREAD)=FACTKK
                  ENDIF
                  MODCOL(1,4,IREI)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    DO J = 1, NSBOX
                      EE = ENERGY_RATE_COEFF(KREAD,TEINL(J),PLS(J),
     .                                       .FALSE.,1)
                      EE = MAX(-100._DP,EE+FCTKKL+DEINL(J))
                      EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
                    END DO
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                    FACREA(KREAD)=LOG(FACTKK)
                  END IF
                  MODCOL(1,4,IREI)=1
                ENDIF
                IF (DELPOT(KREAD).NE.0.D0) THEN
                  DELE=DELPOT(KREAD)
                  IF (NSTORDR >= NRAD) THEN
                    DO J=1,NSBOX
                      EELDS1(IREI,J)=EELDS1(IREI,J)+DELE
                    END DO  
                  END IF
                ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  4.B: HEAVY PARTICLE ENERGY GAIN RATE
C
      EFLAG=IDEZ(ISCDE,3,5)
      IF (EFLAG.EQ.0) THEN
C  4.B1)  RATE = CONST.*RATECOEFF.
        IF (NSTORDR >= NRAD) THEN
          DO 201 J=1,NSBOX
            EHVDS1(IREI,J)=EHEAVY
201       CONTINUE
          NREAHV(IREI)=-1
        ELSE
          NREAHV(IREI)=-1
          EHVDS1(IREI,1)=EHEAVY
        END IF
C     ELSEIF (EFLAG.EQ.1) THEN
C        NOT A VALID OPTION
      ELSEIF (EFLAG.EQ.3) THEN
C  4.B3)  ENERGY RATE = EN.WEIGHTED RATE(TE)
        KREAD=EHEAVY
        MODC=IDEZ(MODCLF(KREAD),5,5)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
              EHVDS1(IREI,J)=ENERGY_RATE_COEFF(KREAD,TEINL(J),0._DP,
     .                .TRUE.,0)*DEIN(J)*FACTKK/(TABDS1(IREI,J)+EPS60)
202         CONTINUE
            NREAHV(IREI)=KREAD
          ELSE
            NREAHV(IREI)=KREAD
            FACREA(KREAD)=LOG(FACTKK)
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTEI '
          CALL EXIT_OWN(1)
        ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEI(IREI,1)=IDEZ(IESTM,1,3)
      IESTEI(IREI,2)=IDEZ(IESTM,2,3)
      IESTEI(IREI,3)=IDEZ(IESTM,3,3)
C
      IF (IESTEI(IREI,1).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,1)=0
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        CALL LEER(1)
        WRITE (iunout,*) 
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,2)=0
      ENDIF
      CALL LEER(1)
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_1(IREI)
C
          DO 510 IAT=1,NATMI
            IA=NSPH+IAT
            PATDS(IREI,0)=PATDS(IREI,0)+
     +                          PATDS(IREI,IAT)
            P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                          P2ND(IREI,IA)
510       CONTINUE
          DO 520 IML=1,NMOLI
            IM=NSPA+IML
            PMLDS(IREI,0)=PMLDS(IREI,0)+
     +                          PMLDS(IREI,IML)
            P2ND(IREI,IM)=P2ND(IREI,IM-1)+
     +                          P2ND(IREI,IM)
520       CONTINUE
          DO 530 ION=1,NIONI
            IO=NSPAM+ION
            PIODS(IREI,0)=PIODS(IREI,0)+
     +                          PIODS(IREI,ION)
            P2ND(IREI,IO)=P2ND(IREI,IO-1)+
     +                          P2ND(IREI,IO)
530       CONTINUE
          DO 540 IPL=1,NPLSI
            PPLDS(IREI,0)=PPLDS(IREI,0)+
     +                          PPLDS(IREI,IPL)
540       CONTINUE
C
          P2NDS(IREI)=PATDS(IREI,0)+PMLDS(IREI,0)+
     .                   PIODS(IREI,0)
          P2N=P2ND(IREI,NSPAMI)
          DO 550 ISPZ=NSPH+1,NSPAMI
            IF (P2N.GT.0.D0)
     .      P2ND(IREI,ISPZ)=P2ND(IREI,ISPZ)/P2N
550       CONTINUE
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY XSTEI_2(IREI)
C
      CALL LEER(1)
      WRITE (iunout,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL LEER(1)
      EI=1.D30
      EA=-1.D30
      imin=0
      imax=0
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELDS1(IREI,IRAD)
        ELSE
          EN=FEELEI1(IREI,IRAD)
        END IF
        if (en < ei) imin=irad
        if (en > ea) imax=irad
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE
      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (iunout,*) 'ELECTRONS: PELDS, CONSTANT ENERGY: EEL'
        WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELDS(IREI),EI
      ELSE
        WRITE (iunout,*) 
     .    'ELECTRONS: PELDS, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (iunout,'(1X,A8,3(1PE12.4))') 'EL      ',PELDS(IREI),EI,EA
      ENDIF
      write (iunout,*) ' imin = ', imin, ' imax = ',imax
C
      EI=1.D30
      EA=-1.D30
      DO 876 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 876
        IF (NSTORDR >= NRAD) THEN
          EN=EHVDS1(IREI,IRAD)
        ELSE
          EN=FEHVDS1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
876   CONTINUE
      IF (PPLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLDS '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (PPLDS(IREI,IPL).NE.0.D0)
     .      WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLDS(IREI,IPL)
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EPLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EPLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2),
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
      IF (P2NDS(IREI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL LEER(1)
        RETURN
      ENDIF
C
      IF (PATDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATDS '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATDS(IREI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLDS '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLDS(IREI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLDS(IREI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLDS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIODS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIODS '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIODS(IREI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIODS(IREI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIODS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIODS '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      RETURN
C
C
C-----------------------------------------------------------------------
C
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (iunout,*) IREI,KK
      CALL EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (iunout,*) IREI
      CALL EXIT_OWN(1)
      END
