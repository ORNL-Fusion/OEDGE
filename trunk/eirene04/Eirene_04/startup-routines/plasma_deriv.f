C
      SUBROUTINE PLASMA_DERIV

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CTRCEI
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CZT1
      USE CGEOM
      USE CSPEI
      USE COMXS

      IMPLICIT NONE

      REAL(DP) :: ZTII, ZTNI, FCT2, FCRG, FCT1, VDION, ZTEI, ZTNE, 
     .            EMPLS, FCT0, TMPLS, DNPLS, DI, AM1, TEF, DEF,
     .            TEI, DEJ, TVACL, DVACL
      REAL(DP) :: COEF(0:8), COEF2D(0:8,0:8)
      REAL(DP), ALLOCATABLE :: DEINTF(:), SUMNI(:), SUMMNI(:)
      INTEGER :: I, IR, IN, IP, IPM, J, IPLS, IOLD, IO, ISW, IREAC, I1,
     .           I0, IPLSTI, IPLSV, IOLDTI, IOLDV


      DO IPLS=1,NPLSI
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF (INDEX(CDENMODEL(IPLS),'FORT.13') > 0) THEN
          CALL ALLOC_BCKGRND
          ALLOCATE(DEINTF(NRAD))
          OPEN (UNIT=13,ACCESS='SEQUENTIAL',FORM='UNFORMATTED')
          REWIND 13
          READ (13,IOSTAT=IO) TEINTF,TIINTF,DEINTF,DIINTF,
     .                        VXINTF,VYINTF,VZINTF
          IF (TRCFLE) WRITE (6,*) 'READ 13: RCMUSR, IO= ',IO
          CLOSE (UNIT=13)
          IF (IO.EQ.0) THEN
            TIIN(IPLSTI,:)=TIINTF(IPLSTI,:)
            DIIN(IPLS,:)=DIINTF(IPLS,:)
            VXIN(IPLSV,:)=VXINTF(IPLSV,:)
            VYIN(IPLSV,:)=VYINTF(IPLSV,:)
            VZIN(IPLSV,:)=VZINTF(IPLSV,:)
          ENDIF
          DEALLOCATE(DEINTF)
        ELSEIF (INDEX(CDENMODEL(IPLS),'BOLTZMANN') > 0) THEN
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
          TIIN(IPLSTI,:)=TIIN(IOLDTI,:)
          VXIN(IPLSV,:)=VXIN(IOLDV,:)
          VYIN(IPLSV,:)=VYIN(IOLDV,:)
          VZIN(IPLSV,:)=VZIN(IOLDV,:)
          DIIN(IPLS,1:NSURF)=MAX(DVAC,DIIN(IOLD,1:NSURF)*
     .        TDMPAR(IPLS)%TDM%G_BOLTZ*
     .        EXP(-TDMPAR(IPLS)%TDM%DELTAE/TIIN(IOLDTI,1:NSURF)))
        END IF
      END DO
C
C  COMPUTE SOME 'DERIVED' PLASMA DATA PROFILES FROM THE INPUT PROFILES
C
C  SET ELECTRON-DENSITY FROM QUASI-NEUTRALITY, AND DRIFT ENERGY (EV)
      TVACL=LOG(TVAC)
      DVACL=LOG(DVAC*1.E-8_DP)
      LGVAC=.TRUE.
      DO 5102 J=1,NSBOX
        DEIN(J)=0.
        DO IPLS=1,NPLSI
          DEIN(J)=DEIN(J)+DBLE(NCHRGP(IPLS))*DIIN(IPLS,J)
        END DO
C  SET 'LOG OF TEMPERATURE AND DENSITY' ARRAYS
        ZTEI=MAX(TVAC,MIN(TEIN(J),1.E10_DP))
        TEINL(J)=LOG(ZTEI)
        ZTNE=MAX(DVAC,MIN(DEIN(J),1.E20_DP))
        DEINL(J)=LOG(ZTNE)
        TMPLS=TEIN(J)
        DNPLS=DEIN(J)
        LGVAC(J,NPLS+1)=TMPLS.LE.TVAC.OR.DNPLS.LE.DVAC
        LGVAC(J,0)     =LGVAC(J,0).AND.LGVAC(J,NPLS+1)
5102  CONTINUE
C
C   SPECIAL DENSITY MODELS
C
      DO IPLS=1,NPLSI
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        IF ((LEN_TRIM(CDENMODEL(IPLS)) > 0) .AND. 
     .      (NCHRGP(IPLS) .NE. 0)) THEN
          WRITE (6,*) ' CHARGE NUMBER MUST BE 0 WITH CORONA DENSITY',
     .                ' MODEL '
          WRITE (6,*) ' IPLS = ',IPLS
          CALL EXIT_OWN(1)
        END IF

        SELECT CASE (CDENMODEL(IPLS))

        CASE ('SAHA      ')
!PB   TO BE WRITTEN
          WRITE (6,*) ' DENSITY PROFILE ACCORDING TO SAHA IS NOT ',
     .                'AVAILABLE '
          WRITE (6,*) ' PLEASE CHOOSE DIFFERENT OPTION FOR ION ',
     .                'DENSITY ',IPLS
          CALL EXIT_OWN(1)

        CASE ('CORONA    ')
          IOLD=TDMPAR(IPLS)%TDM%ISP(1)
          IOLDTI=MPLSTI(IOLD)
          IOLDV=MPLSV(IOLD)
          TIIN(IPLSTI,:)=TIIN(IOLDTI,:)
          VXIN(IPLSV,:)=VXIN(IOLDV,:)
          VYIN(IPLSV,:)=VYIN(IOLDV,:)
          VZIN(IPLSV,:)=VZIN(IOLDV,:)
          CALL SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(1),
     .                 TDMPAR(IPLS)%TDM%H2(1),
     .                 TDMPAR(IPLS)%TDM%REACTION(1),
     .                 TDMPAR(IPLS)%TDM%CR(1))
          COEF(0:8)=CREAC(1:9,1,NREACI+1)
          AM1=1._DP/TDMPAR(IPLS)%TDM%A_CORONA
          DO IR=1,NSURF
            IF (LGVAC(IR,NPLS+1)) THEN
              TEF=TVACL
            ELSE
              TEF=TEINL(IR)
            END IF
            DI=0._DP
            DO I=0,8
              TEI=TEF**I
              DI=DI+COEF(I)*TEI
            END DO
            DI=EXP(DI)
            DIIN(IPLS,IR)=MAX(DIIN(IOLD,IR)*DI*DEIN(IR)*AM1,DVAC)
          END DO

        CASE ('COLRAD    ')
          IF (.NOT.ALLOCATED(SUMNI)) THEN
            ALLOCATE (SUMNI(NRAD))
            ALLOCATE (SUMMNI(NRAD))
          END IF
          SUMNI = 0._DP
          SUMMNI = 0._DP
          DIIN(IPLS,:)=0._DP
          TIIN(IPLSTI,:)=0._DP
          VXIN(IPLSV,:)=0._DP
          VYIN(IPLSV,:)=0._DP
          VZIN(IPLSV,:)=0._DP
          DO IREAC=1,TDMPAR(IPLS)%TDM%NRE
            IOLD=TDMPAR(IPLS)%TDM%ISP(IREAC)
            IOLDTI=MPLSTI(IOLD)
            IOLDV=MPLSV(IOLD)
            CALL SLREAC (NREACI+1,TDMPAR(IPLS)%TDM%FNAME(IREAC),
     .                   TDMPAR(IPLS)%TDM%H2(IREAC),
     .                   TDMPAR(IPLS)%TDM%REACTION(IREAC),
     .                   TDMPAR(IPLS)%TDM%CR(IREAC))
            I1=INDEX(TDMPAR(IPLS)%TDM%H2(IREAC),'.')
            READ (TDMPAR(IPLS)%TDM%H2(IREAC)(I1+1:),*) ISW
            I0=1
            IF (ISW==0) I0=-1
            IF (ISW==1) I0=0
            SELECT CASE (ISW)
            CASE (0,1,2,5,8,11)
              COEF(0:8)=CREAC(1:9,I0,NREACI+1)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                DI=0._DP
                DO I=0,8
                  TEI=TEF**I
                  DI=DI+COEF(I)*TEI
                END DO
                DI=EXP(DI)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+DIIN(IOLD,IR)*DI
                TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+ 
     .                          DIIN(IOLD,IR)*TIIN(IOLDTI,IR)
                VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VXIN(IOLDV,IR)
                VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VYIN(IOLDV,IR)
                VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VZIN(IOLDV,IR)
                SUMNI(IR) = SUMNI(IR) + DIIN(IOLD,IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*DIIN(IOLD,IR)
              END DO
            CASE (3,4,6,7,9,10,12)
              COEF2D(0:8,0:8)=CREAC(1:9,1:9,NREACI+1)
              DO IR=1,NSURF
                IF (LGVAC(IR,NPLS+1)) CYCLE
                TEF=TEINL(IR)
                DEF=LOG(DEIN(IR)*1.E-8_DP)
                DI=0._DP
                DO J=0,8
                  DEJ=DEF**J
                  DO I=0,8
                    TEI=TEF**I
                    DI=DI+COEF2D(I,J)*TEI*DEJ
                  END DO
                END DO
                DI=EXP(DI)
                DIIN(IPLS,IR)=DIIN(IPLS,IR)+DIIN(IOLD,IR)*DI
                TIIN(IPLSTI,IR)=TIIN(IPLSTI,IR)+
     .                          DIIN(IOLD,IR)*TIIN(IOLDTI,IR)
                VXIN(IPLSV,IR)=VXIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VXIN(IOLDV,IR)
                VYIN(IPLSV,IR)=VYIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VYIN(IOLDV,IR)
                VZIN(IPLSV,IR)=VZIN(IPLSV,IR)+
     .                        NMASSP(IOLD)*DIIN(IOLD,IR)*VZIN(IOLDV,IR)
                SUMNI(IR) = SUMNI(IR) + DIIN(IOLD,IR)
                SUMMNI(IR) = SUMMNI(IR) + NMASSP(IOLD)*DIIN(IOLD,IR)
              END DO  ! IR
            CASE DEFAULT
              WRITE (6,*) ' H.',ISW,' NOT FORESEEN IN COLRAD DENSITY',
     .                    ' MODEL '
            END SELECT
          END DO   ! IREAC
          DIIN(IPLS,:)=MAX(DVAC,DIIN(IPLS,:))
!          SUMNI(:)=MAX(DVAC,SUMNI(:))
!          SUMMNI(:)=MAX(DVAC,SUMMNI(:))
          TIIN(IPLSTI,:)=MAX(TVAC,TIIN(IPLSTI,:)/(SUMNI(:)+eps60))
          VXIN(IPLSV,:)=MAX(VVAC,VXIN(IPLSV,:)/(SUMMNI(:)+eps60))
          VYIN(IPLSV,:)=MAX(VVAC,VYIN(IPLSV,:)/(SUMMNI(:)+eps60))
          VZIN(IPLSV,:)=MAX(VVAC,VZIN(IPLSV,:)/(SUMMNI(:)+eps60))
        CASE DEFAULT
!  NOTHING TO BE DONE HERE, ALREADY COMPLETED
        END SELECT
      END DO
      IF (ALLOCATED(SUMNI)) THEN
        DEALLOCATE (SUMNI)
        DEALLOCATE (SUMMNI)
      END IF
C
C  SET  DRIFT ENERGY (EV)
      DO J=1,NSBOX
        DO IPLS=1,NPLSI
          IPLSV=MPLSV(IPLS)
          IF (NLDRFT) THEN
            IF (INDPRO(4) == 8) THEN
              IF(IPLS.EQ.1) THEN
                EDRIFT(IPLS,J)=CVRSSP(IPLS)*VDION(J)**2
              ELSE
                WRITE(6,*)'WARNING! IPLS>1 NO DRIFT!'
                EDRIFT(IPLS,J)=0.D0
              END IF
            ELSE
              EDRIFT(IPLS,J)=CVRSSP(IPLS)*
     .              (VXIN(IPLSV,J)**2+VYIN(IPLSV,J)**2+VZIN(IPLSV,J)**2)
            END IF
          ELSE
            EDRIFT(IPLS,J)=0.D0
          ENDIF
        END DO
      END DO
C
C  SET B_PERP
C
      DO J=1,NSBOX
        IF (ABS(BXIN(J)) > EPS10) THEN
           BYPERP(J) = 1._DP
           BXPERP(J) = -BYIN(J)/BXIN(J)
        ELSEIF (ABS(BYIN(J)) > EPS10) THEN
           BXPERP(J) = 0._DP
           BYPERP(J) = -BXIN(J)/BYIN(J)
        ELSE
           BXPERP(J) = 1._DP
           BYPERP(J) = 0._DP
        END IF
C  CHECK ORIENTATION
        IF (BXIN(J)*BYPERP(J)-BXPERP(J)*BYIN(J) < 0._DP) THEN
           BXPERP(J) = -BXPERP(J)
           BYPERP(J) = -BYPERP(J)
        END IF
      END DO

      DO 5103 J=1,NSBOX
C  SET 'VACUUM REGION FLAGS'
C  LGVAC(...,0)  VACUUM, NO REACTION RATES AT ALL
C  LGVAC(...,IPLS)       NO REACTION RATES FOR BACKGROUND SPECIES IPLS
C  LGVAC(...,NPLS+1)     NO REACTION RATES FOR BACKGROUND ELECTRONS
C                        BUT PERHAPS FOR NEUTRAL BACKGROUND
        DO 5106 IPLS=1,NPLSI
          IPLSTI=MPLSTI(IPLS)
          EMPLS=1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J)
          DNPLS=DIIN(IPLS,J)
          LGVAC(J,IPLS)=EMPLS.LE.TVAC.OR.DNPLS.LE.DVAC
          LGVAC(J,0)   =LGVAC(J,0).AND.LGVAC(J,IPLS)
5106    CONTINUE
5103  CONTINUE
C
      DO 5105 IPLS=1,NPLSI
C  FACTOR FOR MOST PROBABLE SPEED
        FCT0=1./RMASSP(IPLS)*2.*CVEL2A*CVEL2A
C  FACTOR FOR MEAN SPEED
        FCT1=1./RMASSP(IPLS)*8./PIA*CVEL2A*CVEL2A
C  FACTOR FOR ROOT MEAN SQUARE SPEED
        FCT2=1./RMASSP(IPLS)*3.*CVEL2A*CVEL2A
        FCRG=CVEL2A/SQRT(RMASSP(IPLS))
        IPLSTI=MPLSTI(IPLS)
        IPLSV=MPLSV(IPLS)
        DO 5105 J=1,NSBOX
          ZTII=MAX(TVAC,MIN(TIIN(IPLSTI,J),1.E10_DP))
          TIINL(IPLSTI,J)=LOG(ZTII)
          BVIN(IPLSV,J)=BXIN(J)*VXIN(IPLSV,J)+
     .                  BYIN(J)*VYIN(IPLSV,J)+
     .                  BZIN(J)*VZIN(IPLSV,J)
          PARMOM(IPLS,J)=BVIN(IPLSV,J)*SIGN(1._DP,BVIN(IPLSV,J))*
     .                   AMUA*RMASSP(IPLS)
C
C  ZT1: FOR "EFFECTIVE" PLASMA PARTICLE VELOCITY IN CROSS SECTIONS
C       FOR HEAVY PARTICLE INTERACTIONS
C       SQRT(ZT1) IS THE MEAN VELOCITY AT TI=ZTII, TAKEN AS
C       ROOT MEAN SQUARE SPEED
C
          ZT1(IPLS,J)=FCT2*ZTII
C
C  ZRGQ: VARIANCE FOR SAMPLING FROM MAXWELLIAN VELOCITY DISTRIBUTION
C  ZRG=SQRT(ZRGQ) = STANDARD DEVIATION
C
          ZRG(IPLS,J)=FCRG*SQRT(ZTII)
C
          ZTNI=MAX(DVAC,MIN(DIIN(IPLS,J),1.E20_DP))
          DIINL(IPLS,J)=LOG(ZTNI)
5105  CONTINUE
C
      IF (LEVGEO.EQ.3) THEN
        DO 5161 I=1,NPPLG-1
          DO 5162 IP=NPOINT(2,I),NPOINT(1,I+1)-1
            IPM=IP-1
            DO 5163 IPLS=0,NPLS+1
              DO 5163 IR=1,NR1STM
                IN=IR+IPM*NR1ST
                LGVAC(IN,IPLS)=.TRUE.
5163        CONTINUE
5162      CONTINUE
5161    CONTINUE
      ENDIF
C
      RETURN
      END

