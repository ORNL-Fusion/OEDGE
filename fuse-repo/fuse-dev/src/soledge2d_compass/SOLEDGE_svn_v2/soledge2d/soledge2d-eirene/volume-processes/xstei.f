!pb  28.06.06: bug fix for NSTORAM=0 and MODC=1
!pb            FACREA=FACTKK instead of FACREA=log(FACTKK) as single
!pb            polynomial fit is linear in it parameter
!pb  30.08.06: data structure for reaction data redefined
!pb  12.10.06: modcol revised
!pb  22.11.06: flag for shift of first parameter to rate_coeff introduced
!pb  22.11.06: DELPOT introduced
!dr  30.01.07: if lgvac(..,npls+1)  cycle (do not evaluate rates in vacuum)
!pb  20.04.07: allow for third and fourth secondary
!pb  June  07: LHCOL: indicate: direct coupling to col-rad code, rather
!pb            than reading fits or data tables.
!dr  to be done: also for recombination, and generalize to other species (He,...)
!dr             currently: lable H.4 2.1.5 or H.10 2.1.5 not used.
C
      SUBROUTINE EIRENE_XSTEI(RMASS,IREI,ISP,
     .                 IFRST,ISCND,ITHRD,IFRTH,
     .                 EHEAVY,CHRDF0,ISCDE,EELEC,IESTM,
     .                 KK,FACTKK,PLS)
!
!  Set rate coefficients and reaction kinetics for e + ISP collisions.
!
!
c   rmass: mass of incident test particle
c   irei: counting index for this particular electron impact collision
c   isp:  incident test particle species identifier
 
C
C  SET NON DEFAULT ELECTRON IMPACT COLLISION PROCESS NO. IREI
C
      USE EIRMOD_PRECISION
      USE EIRMOD_PARMMOD
      USE EIRMOD_COMUSR
      USE EIRMOD_COMPRT, ONLY: IUNOUT
      USE EIRMOD_CCONA
      USE EIRMOD_CGRID
      USE EIRMOD_COMXS
 
      IMPLICIT NONE
 
      REAL(DP), INTENT(IN) :: RMASS, EHEAVY, CHRDF0, EELEC, FACTKK
      REAL(DP), INTENT(IN) :: PLS(NSTORDR)
      INTEGER, INTENT(IN) :: IREI, ISP, IFRST, ISCND, ITHRD, IFRTH,
     .                       ISCDE, IESTM, KK
      REAL(DP) :: CF(9,0:9)
      REAL(DP) :: EFLAG, CHRDIF, FCTKKL, EIRENE_FEHVDS1, EE, TB, 
     .          EIRENE_FEELEI1, EN,
     .          P2N, EA, EI, ACCINI, ACCINP, ACCMSM, ACCMSI, ACCMAS,
     .          ACCMSA, ACCINA, ACCINM, ACCMSP, ACCINV, COU, 
     .          EIRENE_RATE_COEFF,
     .          EIRENE_ENERGY_RATE_COEFF, DELE, ERATE
      INTEGER :: MODC, KREAD, IM, IA, IERR, J, IPL, I, IP, IRAD, IO,
     .           ION, ISPZ, III, INUM, ITYP, ISPE, ICOUNT, IAT,
     .           IMM, IIO, IAA, IML, IMIN, IMAX
      INTEGER, EXTERNAL :: EIRENE_IDEZ
      LOGICAL :: LHCOL
 
      ITYP=EIRENE_IDEZ(IFRST,1,3)
      INUM=EIRENE_IDEZ(IFRST,2,3)
      ISPE=EIRENE_IDEZ(IFRST,3,3)
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

      IF ((ISPE < 1) .OR. (ISPE > MAXSPC(ITYP))) GOTO 994

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
        ITYP=EIRENE_IDEZ(ISCND,1,3)
        INUM=EIRENE_IDEZ(ISCND,2,3)
        ISPE=EIRENE_IDEZ(ISCND,3,3)
        ICOUNT=2
        GOTO 85
      ELSEIF (ITHRD.NE.0.AND.ICOUNT.EQ.2) THEN
        ITYP=EIRENE_IDEZ(ITHRD,1,3)
        INUM=EIRENE_IDEZ(ITHRD,2,3)
        ISPE=EIRENE_IDEZ(ITHRD,3,3)
        ICOUNT=3
        GOTO 85
      ELSEIF (IFRTH.NE.0.AND.ICOUNT.EQ.3) THEN
        ITYP=EIRENE_IDEZ(IFRTH,1,3)
        INUM=EIRENE_IDEZ(IFRTH,2,3)
        ISPE=EIRENE_IDEZ(IFRTH,3,3)
        ICOUNT=4
        GOTO 85
      ENDIF
C
      IF (ABS(ACCMAS-RMASS).GT.1.D-10) THEN
        WRITE (IUNOUT,*) 'MESSAGE FROM XSTEI.F: '
        WRITE (IUNOUT,*) 'FOR INCIDENT SPECIES ',TEXTS(ISP)
        WRITE (iunout,*)
     .    'MASS CONSERVATION VIOLATED FOR REACT. KK= ',KK
        CALL EIRENE_EXIT_OWN(1)
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
      LHCOL = .FALSE.
C
C  2.A) RATE COEFFICIENT = CONST.
C     TO BE WRITTEN
      IF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.1) THEN
C  2.B) RATE COEFFICIENT(TE)
        IF (NSTORDR >= NRAD) THEN
C  RATE:  (1/S)
C  RATE COEFFICIENT: (CM^3/S)
          DO J=1,NSBOX
            IF (LGVAC(J,NPLS+1)) CYCLE
            COU = EIRENE_RATE_COEFF(KK,TEINL(J),0._DP,.TRUE.,0,ERATE)
            TABDS1(IREI,J)=COU*FACTKK
            IF (IFTFLG(KK,2) < 100)
     .        TABDS1(IREI,J)=TABDS1(IREI,J)*DEIN(J)
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ELSE
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 1
        ENDIF
        MODCOL(1,2,IREI)=1
C     ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.2) THEN
C  2.C) RATE COEFFICIENT(TE,EBEAM)
C  TO BE WRITTEN
C       MODCOL(1,2,IREI)=2
      ELSEIF (EIRENE_IDEZ(MODCLF(KK),3,5).EQ.3) THEN
C  2.D) RATE COEFFICIENT(TE,NE)
        IF (NSTORDR >= NRAD) THEN
          FCTKKL=LOG(FACTKK)
          KREAD=EELEC
          IF ((REACDAT(KK)%RTC%IFIT == 5) .AND.
     .        (EIRENE_IDEZ(ISCDE,5,5) == 3)) THEN
            IF (REACDAT(KREAD)%LRTCEW) THEN
              IF (REACDAT(KREAD)%RTCEW%IFIT == 5) LHCOL=.TRUE.
            END IF
          END IF
          DO J=1,NSBOX
            IF (LGVAC(J,NPLS+1)) CYCLE
            COU = EIRENE_RATE_COEFF(KK,TEINL(J),PLS(J),.FALSE.,1,ERATE)
            TB = COU + FCTKKL
            IF (IFTFLG(KK,2) < 100) TB = TB + DEINL(J)
            TB=MAX(-100._DP,TB)
            TABDS1(IREI,J)=EXP(TB)
            IF (LHCOL) THEN
              EE = MAX(-100._DP,ERATE+FCTKKL+DEINL(J))
              EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
            END IF
          END DO
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ELSE
          NREAEI(IREI) = KK
          JEREAEI(IREI) = 9
        ENDIF
        MODCOL(1,2,IREI)=1
      ELSE
        IERR=1
        GOTO 996
      ENDIF
 
      FACREA(KK,1) = FACTKK
      FACREA(KK,2) = LOG(FACTKK)
C
C  3. ELECTRON MOMENTUM LOSS RATE
C
C
C  4. ENERGY RATES
C
C  4.A: ELECTRON ENERGY LOSS RATES
C
      EFLAG=EIRENE_IDEZ(ISCDE,5,5)
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
                  IF (LGVAC(J,NPLS+1)) CYCLE
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
                IF ((KREAD < 1) .OR. (KREAD > NREACI)) GOTO 998
                MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
                IF (MODC.EQ.1) THEN
                  IF (NSTORDR >=NRAD) THEN
                    DO 102 J=1,NSBOX
                      IF (LGVAC(J,NPLS+1)) CYCLE
                      EELDS1(IREI,J)=-EIRENE_ENERGY_RATE_COEFF(KREAD,
     .                                TEINL(J),
     .                                0._DP,.TRUE.,0)*DEIN(J)*FACTKK/
     .                                (TABDS1(IREI,J)+EPS60)
102                 CONTINUE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=1
                  ENDIF
                  MODCOL(1,4,IREI)=1
C  4.A3) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,EBEAM)
C        TO BE WRITTEN
C  4.A4) ENERGY LOSS RATE OF IMP. ELECTRON = EN.WEIGHTED RATE(TE,NE)
                ELSEIF (MODC.EQ.3) THEN
                  IF (NSTORDR >= NRAD) THEN
                    FCTKKL=LOG(FACTKK)
                    IF (.NOT.LHCOL) THEN
                      DO J = 1, NSBOX
                        IF (LGVAC(J,NPLS+1)) CYCLE
                        EE = EIRENE_ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .                                                PLS(J),.FALSE.,1)
                        EE = MAX(-100._DP,EE+FCTKKL+DEINL(J))
                        EELDS1(IREI,J)=-EXP(EE)/(TABDS1(IREI,J)+EPS60)
                      END DO
                    END IF
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  ELSE
                    NELREI(IREI)=KREAD
                    JELREI(IREI)=9
                  END IF
                  MODCOL(1,4,IREI)=1
                ENDIF
                FACREA(KREAD,1)=FACTKK
                FACREA(KREAD,2)=LOG(FACTKK)
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
      EFLAG=EIRENE_IDEZ(ISCDE,3,5)
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
        MODC=EIRENE_IDEZ(MODCLF(KREAD),5,5)
        FACREA(KREAD,1)=FACTKK
        FACREA(KREAD,2)=LOG(FACTKK)
        IF (MODC.EQ.1) THEN
          IF (NSTORDR >= NRAD) THEN
            DO 202 J=1,NSBOX
              IF (LGVAC(J,NPLS+1)) CYCLE
              EHVDS1(IREI,J)=EIRENE_ENERGY_RATE_COEFF(KREAD,TEINL(J),
     .             0._DP,.TRUE.,0)*DEIN(J)*FACTKK/(TABDS1(IREI,J)+EPS60)
202         CONTINUE
            NREAHV(IREI)=KREAD
          ELSE
            NREAHV(IREI)=KREAD
          END IF
        ELSE
          WRITE (iunout,*) 'INVALID OPTION IN XSTEI '
          CALL EIRENE_EXIT_OWN(1)
        ENDIF
      ELSE
        IERR=2
        GOTO 997
      ENDIF
C
C  ESTIMATOR FOR CONTRIBUTION TO COLLISION RATES FROM THIS REACTION
      IESTEI(IREI,1)=EIRENE_IDEZ(IESTM,1,3)
      IESTEI(IREI,2)=EIRENE_IDEZ(IESTM,2,3)
      IESTEI(IREI,3)=EIRENE_IDEZ(IESTM,3,3)
C
      IF (IESTEI(IREI,1).NE.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR PART.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,1)=0
        CALL EIRENE_LEER(1)
      ENDIF
      IF (IESTEI(IREI,2).NE.0) THEN
        CALL EIRENE_LEER(1)
        WRITE (iunout,*)
     .    'WARNING: COLL.EST NOT AVAILABLE FOR MOM.-BALANCE '
        WRITE (iunout,*) 'IREI = ',IREI
        WRITE (iunout,*) 'AUTOMATICALLY RESET TO TRACKLENGTH ESTIMATOR '
        IESTEI(IREI,2)=0
        CALL EIRENE_LEER(1)
      ENDIF
 
 
 
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY EIRENE_XSTEI_1(IREI)
C
C  SET TOTAL NUMBER OF SECONDARIES BY TYPE OF SECONDARY: P..DS(IREI,0)
C  AND
C  CONVERT SECONDARY SPECIES DISTRIBUTION P2ND(IREI)  INTO
C  CUMMULATIVE DISTRIBUTION (NOT YET NORMALIZED)
 
      DO 510 IAT=1,NATMI
        IA=NSPH+IAT
        PATDS(IREI,0)=PATDS(IREI,0)+
     +                      PATDS(IREI,IAT)
        P2ND(IREI,IA)=P2ND(IREI,IA-1)+
     +                      P2ND(IREI,IA)
510   CONTINUE
      DO 520 IML=1,NMOLI
        IM=NSPA+IML
        PMLDS(IREI,0)=PMLDS(IREI,0)+
     +                      PMLDS(IREI,IML)
        P2ND(IREI,IM)=P2ND(IREI,IM-1)+
     +                      P2ND(IREI,IM)
520   CONTINUE
      DO 530 ION=1,NIONI
        IO=NSPAM+ION
        PIODS(IREI,0)=PIODS(IREI,0)+
     +                      PIODS(IREI,ION)
        P2ND(IREI,IO)=P2ND(IREI,IO-1)+
     +                      P2ND(IREI,IO)
530   CONTINUE
      DO 540 IPL=1,NPLSI
        PPLDS(IREI,0)=PPLDS(IREI,0)+
     +                      PPLDS(IREI,IPL)
540   CONTINUE
C
C  TOTAL NUMBER OF SECONDARIES
      P2NDS(IREI)=PATDS(IREI,0)+PMLDS(IREI,0)+
     .            PIODS(IREI,0)
 
C  NORMALIZE SECONDARY SPECIES DISTRIBUTION P2ND
      P2N=P2ND(IREI,NSPAMI)
      DO 550 ISPZ=NSPH+1,NSPAMI
        IF (P2N.GT.0.D0)
     .  P2ND(IREI,ISPZ)=P2ND(IREI,ISPZ)/P2N
550   CONTINUE
C
      RETURN
C
C-----------------------------------------------------------------------
C
      ENTRY EIRENE_XSTEI_2(IREI)
C
      CALL EIRENE_LEER(2)
      WRITE (iunout,*) 'ELEC. IMPACT REACTION NO. IREI= ',IREI
      CALL EIRENE_LEER(1)
      EI=1.D30
      EA=-1.D30
      imin=0
      imax=0
      DO 875 IRAD=1,NSBOX
        IF (LGVAC(IRAD,NPLS+1)) GOTO 875
        IF (NSTORDR >= NRAD) THEN
          EN=EELDS1(IREI,IRAD)
        ELSE
          EN=EIRENE_FEELEI1(IREI,IRAD)
        END IF
        if (en < ei) imin=irad
        if (en > ea) imax=irad
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
875   CONTINUE
 
      WRITE (iunout,*) 'BACKGROUND SECONDARIES:'
      IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
        WRITE (iunout,*) 'ELECTRONS: PELEI, CONSTANT ENERGY: EEL'
        WRITE (iunout,'(1X,A8,2(1PE12.4))') 'EL      ',PELDS(IREI),EI
      ELSE
        WRITE (iunout,*)
     .    'ELECTRONS: PELEI, ENERGY RANGE: EEL_MIN,EEL_MAX'
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
          EN=EIRENE_FEHVDS1(IREI,IRAD)
        END IF
        EI=MIN(EI,EN)
        EA=MAX(EA,EN)
876   CONTINUE
      IF (PPLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'BULK IONS: PPLEI '
        DO 874 IPL=1,NPLSI
          IP=NSPAMI+IPL
          IF (PPLDS(IREI,IPL).NE.0.D0)
     .      WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IP),PPLDS(IREI,IPL)
874     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EPLEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EPLEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EPLDS(IREI,1),
     .                                 ' * E0 + ',EPLDS(IREI,2),
     .                                 ' * EHEAVY '
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ELSE
        WRITE (iunout,*) 'BULK IONS: NONE'
      ENDIF
      CALL EIRENE_LEER(1)
C
      WRITE (iunout,*) 'TEST PARTICLE SECONDARIES:'
      IF (P2NDS(IREI).EQ.0.D0) THEN
        WRITE (iunout,*) 'NONE'
        CALL EIRENE_LEER(1)
        RETURN
      ENDIF
C
      IF (PATDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'ATOMS    : PATEI '
        DO 871 IAT=1,NATMI
          IA=NSPH+IAT
          IF (PATDS(IREI,IAT).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IA),PATDS(IREI,IAT)
871     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EATEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EATEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EATDS(IREI,0,1),
     .                                 ' * E0 + ',EATDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PMLDS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'MOLECULES: PMLEI '
        DO 872 IML=1,NMOLI
          IM=NSPA+IML
          IF (PMLDS(IREI,IML).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IM),PMLDS(IREI,IML)
872     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EMLEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EMLEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EMLDS(IREI,0,1),
     .                                 ' * E0 + ',EMLDS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
      IF (PIODS(IREI,0).GT.0.D0) THEN
        WRITE (iunout,*) 'TEST IONS: PIOEI '
        DO 873 IIO=1,NIONI
          IO=NSPAM+IIO
          IF (PIODS(IREI,IIO).NE.0.D0)
     .    WRITE (iunout,'(1X,A8,1PE12.4)') TEXTS(IO),PIODS(IREI,IIO)
873     CONTINUE
        IF (ABS((EI-EA)/(EA+EPS60)).LE.EPS10) THEN
          WRITE (iunout,*) 'ENERGY: EIOEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2)*EI
        ELSE
          WRITE (iunout,*) 'ENERGY: EIOEI '
          WRITE (iunout,'(1X,1PE12.4,A8,1PE12.4,A10)') EIODS(IREI,0,1),
     .                                 ' * E0 + ',EIODS(IREI,0,2),
     .                                 ' * EHEAVY'
          WRITE (iunout,*) 'ENERGY RANGE: EHEAVY_MIN, EHEAVY_MAX'
          WRITE (iunout,'(1X,2(1PE12.4))') EI,EA
        ENDIF
      ENDIF
 
 
      CALL EIRENE_LEER(1)
      IF (IESTEI(IREI,1).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR PART.-BALANCE '
      IF (IESTEI(IREI,2).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR MOM.-BALANCE '
      IF (IESTEI(IREI,3).NE.0)
     .   WRITE (IUNOUT,*) 'COLLISION ESTIMATOR FOR EN.-BALANCE '
      CALL EIRENE_LEER(1)
      RETURN
C
C
C-----------------------------------------------------------------------
C
994   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: EXIT CALLED '
      WRITE (iunout,*)
     .  'SPECIES INDEX OF SECONDARY PARTICLE OUT OF RANGE'
      WRITE (iunout,*) 'KK ',KK
      CALL EIRENE_EXIT_OWN(1)
996   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI, MODCLF(KK) ',MODCLF(KK)
      WRITE (iunout,*) IREI,KK
      CALL EIRENE_EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: ISCDE FLAG'
      WRITE (iunout,*) IREI
      CALL EIRENE_EXIT_OWN(1)
998   CONTINUE
      WRITE (iunout,*) 'ERROR IN XSTEI: INVALID KREAD'
      WRITE (iunout,*) IREI,KREAD
      CALL EIRENE_EXIT_OWN(1)
      END
