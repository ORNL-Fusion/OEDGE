C
C
      SUBROUTINE STATIS_BGK
C
C  STANDARD DEVIATION FOR TALLIES NEEDED FOR BGK ITERATION
C  CURRENTLY: BGKV ,  ON SIGMA_BGK(I,...), I=1,NBGVI
C             PDENA,                       I=NBGVI+1,NBGVI+NATMI
C             EDENA,                       I=NBGVI+NATMI+1,NBGVI+2*NATMI
C             PDENM,                       I=NBGVI+2*NATMI+1, .....+NMOLI
C             EDENM                        I=NBGVI+2*NATMI+NMOLI+1, ....+2*NMOLI
C(SEE ALSO: SUBR. OUTEIR)
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE CCONA
      USE CGRID
      USE CSDVI
      USE CSDVI_BGK
      USE COUTAU

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: XN, FSIG, ZFLUX
      INTEGER, INTENT(IN) :: NBIN, NRIN, NPIN, NTIN, NSIN
      LOGICAL, INTENT(IN) :: LP,LT

      INTEGER, ALLOCATABLE, SAVE :: IND(:,:),   IIND(:),    INDSS(:,:)
      REAL(DP), ALLOCATABLE, SAVE :: SD(:), SDD(:)
      REAL(DP) :: XNM, ZFLUXQ, SD2, SD2S, DS, SG, D2S, DSA, SG2, D, DD,
     .          DA, SD1, SD1S
      INTEGER :: IMO, IBGV1, IBGV2, IR, ICO, IAT, NR1, NP2, NT3,
     .           IBGV, J, IIN, IRU, I
      INTEGER, SAVE :: NSB
C
!pb      SAVE
C
      ENTRY STATS0_BGK
C
      IF (.NOT.ALLOCATED(IND)) THEN
        AllOCATE (IND(NRTAL,8))
        AllOCATE (IIND(NRTAL))
        AllOCATE (INDSS(NRTAL,8))
        AllOCATE (SD(0:NRTAL))
        AllOCATE (SDD(0:NRTAL))
        SD=0._DP
        SDD=0._DP
      END IF

      CALL INDTAL(IND,NRTAL,NR1TAL,NP2TAL,NT3TAL,NBMLT)
      DO IR=1,NSBOX_TAL
        IIND(IR)=0
        IIN=0
        DO J=1,8
          IF (IND(IR,J).NE.0) THEN
            IIND(IR)=IIND(IR)+1
            IIN=IIN+1
            INDSS(IR,IIN)=J
          ENDIF
        ENDDO
      ENDDO

      RETURN
C
      ENTRY STATS1_BGK(NBIN,NRIN,NPIN,NTIN,NSIN,LP,LT)
C
      IF (NBGVI.EQ.0) RETURN
      NSB=NBIN
      NR1=NRIN
      NP2=NPIN
      NT3=NTIN
C
      IF (NCLMTS < NCLMT) NCLMTS = NCLMT
      DO I=1,NCLMT
        IR = ICLMT(I)
        DO IIN=2,IIND(IR)
          J=INDSS(IR,IIN)
          IRU=IND(IR,J)
          IF (IMETCL(IRU) == 0) THEN
            NCLMTS = NCLMTS+1
            IMETCL(IRU) = NCLMTS
            ICLMT(NCLMTS) = IRU
          END IF
        END DO
      END DO
C
C  STATISTICS FOR BGKV
C
      DO IBGV=1,NBGVI
        IF (LMETSP(NSPAN(NTALB)+IBGV-1)) THEN
          SD1S=0.
!          SD = 0.D0
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=BGKV(IBGV,IR)-SDVIA_BGK(IBGV,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV,IR)=BGKV(IBGV,IR)
            SD(IR) = SD1
            DO IIN=2,IIND(IR)
              J=INDSS(IR,IIN)
              IRU=IND(IR,J)
              SD(IRU)=SD(IRU)+SD1
            END DO
          END DO

          DO ICO = 1,NCLMTS
            IR = ICLMT(ICO)
            SD1=SD(IR)
            SIGMA_BGK(IBGV,IR)=SIGMA_BGK(IBGV,IR)+SD1*SD1
            SD(IR)=0._DP
          END DO
          SGMS_BGK(IBGV)=SGMS_BGK(IBGV)+SD1S*SD1S
        END IF
      END DO

C  STATISTICS FOR PDENA AND EDENA
C
      DO IAT=1,NATMI
        IF (LMETSP(NSPAN(1)+IAT-1)) THEN
          IBGV1=NBGVI+IAT
          IBGV2=NBGVI+NATMI+IAT
          SD1S=0.
          SD2S=0.
!pb          SD = 0.D0
!pb          SDD = 0.D0
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=PDENA(IAT,IR)-SDVIA_BGK(IBGV1,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV1,IR)=PDENA(IAT,IR)

            SD2=EDENA(IAT,IR)-SDVIA_BGK(IBGV2,IR)
            SD2S=SD2S+SD2
            SDVIA_BGK(IBGV2,IR)=EDENA(IAT,IR)

            SD(IR) = SD1
            SDD(IR) = SD2
            DO IIN=2,IIND(IR)
              J=INDSS(IR,IIN)
              IRU=IND(IR,J)
              SD(IRU)=SD(IRU)+SD1
              SDD(IRU)=SDD(IRU)+SD2
            END DO
          END DO

          DO ICO = 1,NCLMTS
            IR = ICLMT(ICO)
            SD1=SD(IR)
            SD2=SDD(IR)
            SIGMA_BGK(IBGV1,IR)=SIGMA_BGK(IBGV1,IR)+SD1*SD1
            SIGMA_BGK(IBGV2,IR)=SIGMA_BGK(IBGV2,IR)+SD2*SD2
            SD(IR)=0._DP
            SDD(IR)=0._DP
          END DO
          SGMS_BGK(IBGV1)=SGMS_BGK(IBGV1)+SD1S*SD1S
          SGMS_BGK(IBGV2)=SGMS_BGK(IBGV2)+SD2S*SD2S
        END IF
      END DO

C  STATISTICS FOR PDENM AND EDENM
C
      DO IMO=1,NMOLI
        IF (LMETSP(NSPAN(2)+IMO-1)) THEN
          IBGV1=NBGVI+2*NATMI+IMO
          IBGV2=NBGVI+2*NATMI+NMOLI+IMO
          SD1S=0.
          SD2S=0.
!          SD = 0.D0
!          SDD = 0.D0
          DO ICO = 1,NCLMT
            IR = ICLMT(ICO)
            SD1=PDENM(IMO,IR)-SDVIA_BGK(IBGV1,IR)
            SD1S=SD1S+SD1
            SDVIA_BGK(IBGV1,IR)=PDENM(IMO,IR)

            SD2=EDENM(IMO,IR)-SDVIA_BGK(IBGV2,IR)
            SD2S=SD2S+SD2
            SDVIA_BGK(IBGV2,IR)=EDENM(IMO,IR)

            SD(IR) = SD1
            SDD(IR) = SD2
            DO IIN=2,IIND(IR)
              J=INDSS(IR,IIN)
              IRU=IND(IR,J)
              SD(IRU)=SD(IRU)+SD1
              SDD(IRU)=SDD(IRU)+SD2
            END DO
          END DO

          DO ICO = 1,NCLMTS
            IR = ICLMT(ICO)
            SD1=SD(IR)
            SD2=SDD(IR)
            SIGMA_BGK(IBGV1,IR)=SIGMA_BGK(IBGV1,IR)+SD1*SD1
            SIGMA_BGK(IBGV2,IR)=SIGMA_BGK(IBGV2,IR)+SD2*SD2
            SD(IR)=0._DP
            SDD(IR)=0._DP
          END DO
          SGMS_BGK(IBGV1)=SGMS_BGK(IBGV1)+SD1S*SD1S
          SGMS_BGK(IBGV2)=SGMS_BGK(IBGV2)+SD2S*SD2S
        END IF
      END DO
C
C
1020  CONTINUE
      RETURN
C
      ENTRY STATS2_BGK(XN,FSIG,ZFLUX)
C
C  1. FALL  ALLE BEITRAEGE GLEICHES VORZEICHEN: SIG ZWISCHEN 0 UND 1
C  2. FALL  NEGATIVE UND POSITIVE BEITRAGE KOMMEN VOR:
C           LT. FORMEL SIND AUCH WERTE GROESSER 1  MOEGLICH.
C
      XNM=XN-1.
      IF (XNM.LE.0.) RETURN
      ZFLUXQ=ZFLUX*ZFLUX
C
      IF (NBGVI.EQ.0) GOTO 2200
C
C   STATISTICS FOR BGKV
      DO 2112 IBGV=1,NBGVI
!        SD=0.
        DS=0.
        DO IR=1,NSB
          SD1=BGKV(IBGV,IR)
          DS=DS+SD1
          DO IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO 2111 IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
2111    CONTINUE
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
2112  CONTINUE
C
C   STATISTICS FOR PDENA
      DO IAT=1,NATMI
        IBGV=NBGVI+IAT
!pb        SD=0.
        DS=0.
        DO IR=1,NSB
          SD1=PDENA(IAT,IR)
          DS=DS+SD1
          DO IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR EDENA
      DO IAT=1,NATMI
        IBGV=NBGVI+NATMI+IAT
!pb        SD=0.
        DS=0.
        DO IR=1,NSB
          SD1=EDENA(IAT,IR)
          DS=DS+SD1
          DO IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR PDENM
      DO IMO=1,NMOLI
        IBGV=NBGVI+2*NATMI+IMO
!pb        SD=0.
        DS=0.
        DO IR=1,NSB
          SD1=PDENM(IMO,IR)
          DS=DS+SD1
          DO IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
C   STATISTICS FOR EDENM
      DO IMO=1,NMOLI
        IBGV=NBGVI+2*NATMI+NMOLI+IMO
!pb        SD=0.
        DS=0.
        DO IR=1,NSB
          SD1=EDENM(IMO,IR)
          DS=DS+SD1
          DO IIN=1,IIND(IR)
            J=INDSS(IR,IIN)
            IRU=IND(IR,J)
            SD(IRU)=SD(IRU)+SD1
          END DO
        END DO

        DO IR=1,NSB
          D=SD(IR)
          DD=D*D
          DA=ABS(D)
          SG2=MAX(0._DP,SIGMA_BGK(IBGV,IR)-DD/XN)
C RELATIV STANDARD DEVIATION
          SG=SQRT(SG2)/(DA+EPS60)
          SIGMA_BGK(IBGV,IR)=SG*FSIG
C CUMULATED VARIANCE FOR SUM OVER STRATA
          STV_BGK(IBGV,IR)=STV_BGK(IBGV,IR)+SG2*ZFLUXQ/XNM/XN
          EE_BGK(IBGV,IR)=EE_BGK(IBGV,IR)+D*ZFLUX/XN
          SD(IR)=0._DP
        END DO
        D2S=DS*DS
        DSA=ABS(DS)
        SG2=MAX(0._DP,SGMS_BGK(IBGV)-D2S/XN)
        SG=SQRT(SG2)/(DSA+EPS60)
        SGMS_BGK(IBGV)=SG*FSIG
C
        STVS_BGK(IBGV)=STVS_BGK(IBGV)+SG2*ZFLUXQ/XNM/XN
        EES_BGK(IBGV)=EES_BGK(IBGV)+DS*ZFLUX/XN
      END DO
C
2200  CONTINUE
      RETURN
      END
