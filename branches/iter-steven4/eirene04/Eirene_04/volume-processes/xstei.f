C
C
      SUBROUTINE XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CGRID
      USE COMXS

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, CHRDF0, EELEC, FACTKK
      REAL(DP), INTENT(IN) :: PLS(NRAD)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ISCDE, IESTM, KK
      REAL(DP) :: COUN(0:9,NSTORDR), CF(9,0:9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL, FEHVDS1, EE, TB, FEELEI1, EN,
     .          P2N, EA, EI, ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV
      INTEGER :: MODC, KREAD, IM, IA, IERR, J, IPL, I, IP, IRAD, IO,
     .           ION, ISPZ, III, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML
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
        WRITE (6,*) 'MASS CONSERVATION VIOLATED FOR REACT. KK'
        WRITE (6,*) KK,ISP
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
          IF (MOD(IFTFLG(KK,2),100) == 10) THEN
C  RATE:  (1/S)
            COUN(1,1:NSBOX)=CREAC(1,1,KK)
          ELSE
C  RATE COEFFICIENT: (CM^3/S)
            CALL CDEF (TEINL,1,1,KK,COUN,NSBOX,CF,.TRUE.,.FALSE.,
     .                .TRUE.)
          END IF

          IF (IFTFLG(KK,2) < 100) THEN
            DO J=1,NSBOX
              TABDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK
            END DO
          ELSE
            DO J=1,NSBOX
              TABDS1(IREI,J)=COUN(1,J)*FACTKK
            END DO
          END IF
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ELSE
          FACREA(KK) = LOG(FACTKK)
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,ISP,1)=1
C     ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,ISP,1)=2
      ELSEIF (IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          CALL CDEF (TEINL,1,9,KK,COUN,NSBOX,CF,.FALSE.,.FALSE.,
     .              .TRUE.)
          DO 93 J=1,NSBOX
            TABDS1(IREI,J)=COUN(9,J)
93        CONTINUE
          DO 95 I=8,1,-1
            DO 94 J=1,NSBOX
              TABDS1(IREI,J)=TABDS1(IREI,J)*PLS(J)+COUN(I,J)
94          CONTINUE
95        CONTINUE
          FCTKKL=LOG(FACTKK)
          IF (IFTFLG(KK,2) < 100) THEN
            DO 97 J=1,NSBOX
              TB=MAX(-100._DP,TABDS1(IREI,J)+FCTKKL+DEINL(J))
              TABDS1(IREI,J)=EXP(TB)
97          CONTINUE
          ELSE
            DO J=1,NSBOX
              TB=MAX(-100._DP,TABDS1(IREI,J)+FCTKKL)
              TABDS1(IREI,J)=EXP(TB)
            END DO
          END IF
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ELSE
          FACREA(KK) = LOG(FACTKK)
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ENDIF
        MODCOL(1,2,ISP,1)=1
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
              MODCOL(1,4,ISP,1)=1
      ELSEIF (EFLAG.EQ.1) THEN
              IF (NSTORDR >= NRAD) THEN
                DO 103 J=1,NSBOX
                  EELDS1(IREI,J)=-1.5*TEIN(J)
103             CONTINUE
                NELREI(IREI)=-3
              ELSE
                NELREI(IREI)=-3
              END IF
              MODCOL(1,4,ISP,1)=1
      ELSEIF (EFLAG.EQ.3) THEN
C  4.A2) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE)
                KREAD=EELEC
                MODC=IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                         .FALSE.,.TRUE.)
                    DO 102 J=1,NSBOX
                      EELDS1(IREI,J)=-COUN(1,J)*DEIN(J)*FACTKK/
     .                               (TABDS1(IREI,J)+EPS60)
102                 CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                    FACREA(KREAD)=LOG(FACTKK)
                  ENDIF
                  MODCOL(1,4,ISP,1)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    CALL CDEF (TEINL,1,9,KREAD,COUN,NSBOX,CF,.FALSE.,
     .                         .FALSE.,.TRUE.)
                    DO 106 J=1,NSBOX
                      EELDS1(IREI,J)=COUN(9,J)
106                 CONTINUE
                    DO 108 I=8,1,-1
                      DO 107 J=1,NSBOX
                        EELDS1(IREI,J)=EELDS1(IREI,J)*PLS(J)+COUN(I,J)
107                   CONTINUE
108                 CONTINUE
                    FCTKKL=LOG(FACTKK)
                    DO 110 J=1,NSBOX
                      EE=MAX(-100._DP,EELDS1(IREI,J)+FCTKKL+DEINL(J))
                      EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
110                 CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                    FACREA(KREAD)=LOG(FACTKK)
                  END IF
                  MODCOL(1,4,ISP,1)=1
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
            CALL CDEF (TEINL,1,1,KREAD,COUN,NSBOX,CF,.TRUE.,
     .                .FALSE.,.TRUE.)
            DO 202 J=1,NSBOX
              EHVDS1(IREI,J)=COUN(1,J)*DEIN(J)*FACTKK/
     .                       (TABDS1(IREI,J)+EPS60)
202         CONTINUE
            NREAHV(IREI)=KREAD
          ELSE
            NREAHV(IREI)=KREAD
            FACREA(KREAD)=LOG(FACTKK)
          END IF
        ELSE
          WRITE (6,*) 'INVALID OPTION IN XSTEI '
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
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (6,*) 'IREI = ',IREI
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,1)=0
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        CALL LEER(1)
        WRITE (6,*) 'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (6,*) 'IREI = ',IREI
        WRITE (6,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
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
      WRITE (6,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL LEER(1)
      EI=1.D30
      EA=-1.D30
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELDS1(IREI,IRAD)
        ELSE
          EN=FEELEI1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE
      WRITE (6,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (6,*) 'ELECTRONS: PELDS, CONSTANT ENERGY: EEL'
        WRITE (6,'(1X,A8,2(1PE12.4))') 'EL      ',PELDS(IREI),EI
      ELSE
        WRITE (6,*) 'ELECTRONS: PELDS, ENERGY RANGE: EEL_MIN,EEL_MAX'
        WRITE (6,'(1X,A8,3(1PE12.4))') 'EL      ',PELDS(IREI),EI,EA
      ENDIF
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
        WRITE (6,*) 'BULK IONS: PPLDS '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (PPLDS(IREI,IPL).NE.0.D0)
     .      WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IP),PPLDS(IREI,IPL)
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EPLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EPLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2),
     .                                 ' * EHEAVY '
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (6,*) 'BULK IONS: NONE'
      ENDIF
      CALL LEER(1)
C
      WRITE (6,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NDS(IREI).EQ.0.D0) THEN
        WRITE (6,*) 'NONE'
        CALL LEER(1)
        RETURN
      ENDIF
C
      IF (PATDS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'ATOMS    : PATDS '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IA),PATDS(IREI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EATDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EATDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLDS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'MOLECULES: PMLDS '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLDS(IREI,IML).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IM),PMLDS(IREI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EMLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EMLDS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIODS(IREI,0).GT.0.D0) THEN
        WRITE (6,*) 'TEST IONS: PIODS '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIODS(IREI,IIO).NE.0.D0)
     .    WRITE (6,'(1X,A8,1PE12.4)') TEXTS(IO),PIODS(IREI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (6,*) 'ENERGY: EIODS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2)*EI
        ELSE
          WRITE (6,*) 'ENERGY: EIODS '
          WRITE (6,'(1X,1PE12.4,A8,1PE12.4,A10)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (6,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (6,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      RETURN
C
C
C-----------------------------------------------------------------------
C
996   CONTINUE
      WRITE (6,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (6,*) IREI,KK
      CALL EXIT_OWN(1)
997   CONTINUE
      WRITE (6,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (6,*) IREI
      CALL EXIT_OWN(1)
      END
