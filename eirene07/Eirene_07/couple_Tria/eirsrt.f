C

      SUBROUTINE EIRSRT(LSTOP,LTIME,DELTAT,FLUXES,
     .                  B2BRM,B2RD,B2Q,B2VP)

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CESTIM
      USE BRASPOI
      USE COMSOU
      USE CLOGAU
      USE CTRCEI
      USE CGRID
      USE CPOLYG
      USE CGEOM
      USE CSPEI
      USE COUTAU
      USE CSPEZ
      USE CINIT
      USE CCONA
      USE COMXS
      USE CSDVI
      USE CSDVI_BGK
      USE CSDVI_COP
      USE CZT1
      USE CCOUPL
      USE COMPRT
      USE COMNNL
      USE BRASCL

      IMPLICIT NONE
C
      REAL(DP), INTENT(IN) :: FLUXES(NSTRA)
      REAL(DP), INTENT(IN) :: DELTAT, B2BRM, B2RD, B2Q, B2VP
      LOGICAL, INTENT(IN) :: LSTOP, LTIME

      REAL(DP) :: FLUXS(NSTRA)
      REAL(DP) :: FTABEI1, FEELEI1, FLXI, ESIG, RESET_SECOND, DUMMY,
     .          SECOND_OWN, DTIMVO
      INTEGER :: IN, IAEI, IRDS, IIDS, ICPV, IMDS, IFIRST, K, JC, NDXY,
     .           J, IRC, NREC10, NREC11, IPLSTI
      REAL(DP), ALLOCATABLE :: OUTAU(:)
      INTEGER, ALLOCATABLE :: IHELP(:)
C
      TYPE(CELLSIM), POINTER :: CPSIM
      TYPE(CELLMUL), POINTER :: CPMUL

C
      SAVE
      DATA IFIRST/0/
C
      IF (LTIME) THEN
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        DUMMY=RESET_SECOND()
        IF(IFIRST.EQ.0) THEN
C
          CALL GRSTRT(35,8)
C
C  READ FORMATTED INPUT FILE IUNIN
C  AND RUN EIRENE FOR ONE TIME-CYCLE: ITIMV=1
C  WITH OR WITHOUT INITIAL DISTRIBUTION ON FILE FT15 (NFILE-J FLAG)
C  AS FINAL STRATUM
C  EXPECT PLASMA DATA ON FORT.31 (NLPLAS=.FALSE.)
C
          CALL EIRENE(DELTAT,.FALSE.,.FALSE.,1,.TRUE.)
C
C  EIRENE RUN DONE. CENSUS ARRAY WRITTEN
C  NOW ITIMV=ITIMV+1, NLPLAS=.TRUE.
C
          DO 3 ISTRA=1,NSTRAI
            FLUXS(ISTRA)=FLUX(ISTRA)
3         CONTINUE
          IFIRST=1
        ELSE
C
C  NOW: NLPLAS=.TRUE., I.E., PLASMA DATA EXPECTED ON BRAEIR
C  NOW: ITIMV=ITIMV+1
C  BUT: COMMON BRAEIR REDONE IN EXTERNAL CODE.
C  REACTIVATE INDEX MAPPING, EVEN WITHOUT READING INPUT BLOCK 14 AGAIN
          NCUTB_SAVE=NCUTB
C
          DTIMVO=DTIMV
          DTIMVN=DELTAT
C
C-----------------------------------------------------------------------
C
C  STRATA 1 TO NTARGI ARE SCALED IN PLASMA CODE  (RECYCLING STRATA)
C
C     RETURN TO PLASMA CODE THE PROFILES PER UNIT SOURCE STRENGTH
C     IE. THE PROFILES ARE SCALED BY 1./FLUX(ISTRA) BEFORE RETURN
C
C  STRATA NTARGI+1 TO NSTRAI-1  ARE SCALED BY EIRENE
C
C     (EG. GAS PUFF, VOLUME RECOMBINATION, ETC.)
C     THEY MAY BE RESCALED BY PLASMA CODE FACTORS: FLUXES(ISTRA)
C     RETURN TO PLASMA CODE THE PROFILES SCALED WITH
C     SOURCE STRENGTH: FLUX(ISTRA) (AMP)
C
C  STRATUM NSTRAI IS RESCALED WITH RATIO OF OLD TO NEW TIMESTEP
C
C     RETURN TO PLASMA CODE THE PROFILES WITH FLUX(ISTRA) (AMP)
C
          DO ISTRA=NTARGI+1,NSTRAI-1
            IF (FLUXES(ISTRA).NE.0.) THEN
              FLUX(ISTRA)=FLUXS(ISTRA)*FLUXES(ISTRA)*ELCHA
            ELSE
              FLUX(ISTRA)=FLUXS(ISTRA)
            ENDIF
          ENDDO
C
          IF (DTIMVN.NE.DTIMVO) THEN
            FLUX(NSTRAI)=FLUX(NSTRAI)*DTIMVO/DTIMVN
C
            WRITE (iunout,*) 'FLUX IS RESCALED BY DTIMV_OLD/DTIMV_NEW '
            CALL MASR1('FLUX    ',FLUX(NSTRAI))
            CALL LEER(1)
          ENDIF
C
C-----------------------------------------------------------------------
C
          DTIMV=DTIMVN
C
C  RUN EIRENE ON TIMESTEP DTIMV
C  THEN CALL INTERFACING ROUTINE AT ENTRY IF3COP (FROM EIRENE MAIN)
C
          IITER=1
          IPRNLI=0
          CALL EIRENE_COUPLE (LSTOP,1)
          IF (LSTOP) THEN
            CALL GREND
          ENDIF
        ENDIF
        CALL LEER(2)
        WRITE(*,*) 'EIRENE USED ',SECOND_OWN(),' CPU SECONDS'
        CALL LEER(2)
C
        RETURN
C
      ELSEIF (.NOT.LTIME) THEN
C
        IF (IFIRST.GE.1) GOTO 10000
C
        CALL ALLOC_COMUSR(1)
        CALL SETPRM
C
        NREC10=1500
        OPEN (UNIT=10,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC10)
        NREC11=NOUTAU
        OPEN (UNIT=11,ACCESS='DIRECT',FORM='UNFORMATTED',RECL=8*NREC11)
C
        IRC=3
        READ (11,REC=IRC) RCCPL
        ALLOCATE (IHELP(NOUTAU))
        JC=0
        IRC=IRC+1
        READ (11,REC=IRC) IHELP
        DO K=1,NPTRGT
          DO J=1,10*NSTEP
            JC=JC+1
            ICCPL1(J,K)=IHELP(JC)
            IF (JC == NOUTAU) THEN
              IRC=IRC+1
              READ (11,REC=IRC) IHELP
              JC=0
            END IF
          END DO
        END DO
        DEALLOCATE (IHELP)
        IRC=IRC+1
        READ (11,REC=IRC) ICCPL2
        IRC=IRC+1
        READ (11,REC=IRC) LCCPL
        IF (LTRCFL) WRITE (iunout,*) 
     .    'READ DATA FOR SHORT CYCLE, SAVED ON  FT11 ',IRC
C
        NLPLG=LNLPLG
        NLDRFT=LNLDRF
        TRCFLE=LTRCFL
        NSTRAI=NSTRI
        DO 1 ISTRA=1,NSTRAI
          NLVOL(ISTRA)=LNLVOL(ISTRA)
1       CONTINUE
        NMODE=NMODEI
        NFILEN=NFILNN
C
C
        INDPRO(1)=6
        INDPRO(2)=6
        INDPRO(3)=6
        INDPRO(4)=6
        INDPRO(5)=6
C
        CALL SETCON
        CALL RGEOM(TRCFLE)
        CALL RPLAM(TRCFLE,0)
C
        IF (TRCFLE) 
     .    WRITE (iunout,*) 'READ DATA FOR EIRENE RECALL OPTION'
        IRC=1
        READ (11,REC=IRC) LOGATM,LOGION,LOGMOL,LOGPLS,LOGPHOT
        IF (TRCFLE) WRITE (iunout,*) 
     .    'DATA FOR RECALL OPTION READ FROM  FT11, IRC=1 '
        IRC=2
        ALLOCATE (OUTAU(NOUTAU))
        CALL WRITE_COUTAU (OUTAU, IUNOUT)
        READ (11,REC=IRC) OUTAU
        DEALLOCATE (OUTAU)
        IF (TRCFLE) WRITE (iunout,*) 
     .    'DATA FOR RECALL OPTION READ FROM  FT11, IRC=2 '
C
        CALL INTER0
C
        NDXY=(NDXA-1)*NR1ST+NDYA
C
        IF (.NOT.ALLOCATED(SPLODA)) CALL ALLOC_BRASCL
        CALL INIT_BRASCL1
C
C  INITIAL: ATOMS, EI-PROCESSES
C
        DO 21 IATM=1,NATMI
        DO 21 IPLS=1,NPLSI
        DO 21 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 21
          DO 22 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLODA(IN,IATM,IPLS)=SPLODA(IN,IATM,IPLS)+
     .                        TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLODA(IN,IATM,IPLS)=SPLODA(IN,IATM,IPLS)+
     .                        FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
22        CONTINUE
21      CONTINUE
        DO 23 IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          DO 24 IN=1,NDXY
            SEIODA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
24        CONTINUE
23      CONTINUE
C
        DO 25 IATM=1,NATMI
        DO 25 IAEI=1,NAEII(IATM)
          IRDS=LGAEI(IATM,IAEI)
          DO 25 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEEODA(IN,IATM)=SEEODA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEEODA(IN,IATM)=SEEODA(IN,IATM)+FEELEI1(IRDS,IN)*
     .                                        FTABEI1(IRDS,IN)
            END IF
25      CONTINUE
C
C  INITIAL: TEST IONS, EI-PROCESSES
C
        DO 26 IION=1,NIONI
        DO 26 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          DO 26 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEEODI(IN,IION)=SEEODI(IN,IION)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEEODI(IN,IION)=SEEODI(IN,IION)+FEELEI1(IRDS,IN)*
     .                                        FTABEI1(IRDS,IN)
            END IF
26      CONTINUE
C
        DO 27 IION=1,NIONI
        DO 27 IPLS=1,NPLSI
        DO 27 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 27
          DO 28 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLODI(IN,IION,IPLS)=SPLODI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLODI(IN,IION,IPLS)=SPLODI(IN,IION,IPLS)+
     .                             FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ENDIF
28        CONTINUE
27      CONTINUE
C
        DO 29 IION=1,NIONI
        DO 29 IPLS=1,NPLSI
        DO 29 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 30 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEIODI(IN,IION)=SEIODI(IN,IION)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEIODI(IN,IION)=SEIODI(IN,IION)+FTABEI1(IRDS,IN)*ESIG
            END IF
30        CONTINUE
29      CONTINUE
C
C
C  INITIAL: MOLECULES, EI-PROCESSES
C
        DO 35 IMOL=1,NMOLI
        DO 35 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 35 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEEODM(IN,IMOL)=SEEODM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEEODM(IN,IMOL)=SEEODM(IN,IMOL)+FEELEI1(IRDS,IN)*
     .                                        FTABEI1(IRDS,IN)
            ENDIF
35      CONTINUE
C
        DO 47 IMOL=1,NMOLI
        DO 47 IPLS=1,NPLSI
        DO 47 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 47
          DO 48 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLODM(IN,IMOL,IPLS)=SPLODM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLODM(IN,IMOL,IPLS)=SPLODM(IN,IMOL,IPLS)+
     .                             FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
48        CONTINUE
47      CONTINUE
C
        DO 49 IMOL=1,NMOLI
        DO 49 IPLS=1,NPLSI
        DO 49 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 50 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEIODM(IN,IMOL)=SEIODM(IN,IMOL)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEIODM(IN,IMOL)=SEIODM(IN,IMOL)+FTABEI1(IRDS,IN)*ESIG
            END IF
50        CONTINUE
49      CONTINUE
C
C
C
        DO 60 ISTRA=1,NSTRAI
          IF (ISTRA.EQ.IESTR) THEN
C  NOTHING TO BE DONE
          ELSEIF ((NFILEN.EQ.1.OR.NFILEN.EQ.2).AND.ISTRA.NE.IESTR) THEN
            IESTR=ISTRA
            CALL RSTRT(ISTRA,NSTRAI,NESTM1,NESTM2,NADSPC,
     .                 ESTIMV,ESTIMS,ESTIML,
     .                 NSDVI1,SDVI1,NSDVI2,SDVI2,
     .                 NSDVC1,SIGMAC,NSDVC2,SGMCS,
     .                 NSBGK,SIGMA_BGK,NBGV_STAT,SGMS_BGK,
     .                 NSCOP,SIGMA_COP,NCPV_STAT,SGMS_COP,
     .                 NSIGI_SPC,TRCFLE)
            IF (NLSYMP(ISTRA).OR.NLSYMT(ISTRA)) THEN
              CALL SYMET(ESTIMV,NTALV,NRTAL,NR1TAL,NP2TAL,NT3TAL,
     .                   NADDV,NFIRST,NLSYMP(ISTRA),NLSYMT(ISTRA))
            ENDIF
          ELSE
            WRITE (iunout,*) 'ERROR IN EIRSRT: STRATUM ISTRA= ',ISTRA
            WRITE (iunout,*) 'IS NOT AVAILABLE. EXIT CALLED'
            CALL EXIT_OWN(1)
          ENDIF
C
C  SAVE EIRENE TALLIES, SCALE PER UNIT FLUX (AMP), ON COMMON BRASCL
C  WTOTP IS NEGATIVE IN EIRENE (SINK FOR IONS)
C  ALL STRATA WHICH ARE NOT SPECIFIED BY INPUT BLOCK 14 (FROM
C  PLASMA CODE DATA) ARE NOT RESCALED HERE
C
          IF (ISTRA.LE.NTARGI.AND.WTOTP(0,ISTRA).NE.0.) THEN
            FLXI=-1./WTOTP(0,ISTRA)
          ELSEIF (ISTRA.LE.NTARGI.AND.WTOTP(0,ISTRA).EQ.0.) THEN
            GOTO 60
          ELSEIF (ISTRA.GT.NTARGI) THEN
            FLXI=1.
          ENDIF

          NULLIFY(PAPLS(ISTRA)%PMUL)
          NULLIFY(PMPLS(ISTRA)%PMUL)
          NULLIFY(PIPLS(ISTRA)%PMUL)

          NULLIFY(EAELS(ISTRA)%PSIM)
          NULLIFY(EMELS(ISTRA)%PSIM)
          NULLIFY(EIELS(ISTRA)%PSIM)
          NULLIFY(EAPLS(ISTRA)%PSIM)
          NULLIFY(EMPLS(ISTRA)%PSIM)
          NULLIFY(EIPLS(ISTRA)%PSIM)

          NULLIFY(PDENAS(ISTRA)%PMUL)
          NULLIFY(PDENMS(ISTRA)%PMUL)
          NULLIFY(PDENIS(ISTRA)%PMUL)
          NULLIFY(EDENAS(ISTRA)%PMUL)

          NULLIFY(COPVS(ISTRA)%PMUL)

          NULLIFY(MAPLS(ISTRA)%PMUL)
          NULLIFY(MMPLS(ISTRA)%PMUL)
          NULLIFY(MIPLS(ISTRA)%PMUL)
          NULLIFY(MPHPLS(ISTRA)%PMUL)

          DO IPLS=1,NPLSI
            DO IN=1,NDXY
              IF (PAPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = PAPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => PAPLS(ISTRA)%PMUL
                PAPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (PMPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = PMPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => PMPLS(ISTRA)%PMUL
                PMPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (PIPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = PIPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => PIPLS(ISTRA)%PMUL
                PIPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (MAPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = MAPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => MAPLS(ISTRA)%PMUL
                MAPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (MMPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = MMPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => MMPLS(ISTRA)%PMUL
                MMPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (MIPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = MIPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => MIPLS(ISTRA)%PMUL
                MIPLS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (MPHPL(IPLS,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IPLS
                CPMUL%ICM = IN
                CPMUL%VALUEM = MPHPL(IPLS,IN)*FLXI
                CPMUL%NXTMUL => MPHPLS(ISTRA)%PMUL
                MPHPLS(ISTRA)%PMUL => CPMUL
              ENDIF
            ENDDO
          ENDDO
!PB       DO IN=1,NDXY
          DO IN=1,NSBOX_TAL
            IF (EAEL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EAEL(IN)*FLXI
              CPSIM%NXTSIM => EAELS(ISTRA)%PSIM
              EAELS(ISTRA)%PSIM => CPSIM
            ENDIF
            IF (EMEL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EMEL(IN)*FLXI
              CPSIM%NXTSIM => EMELS(ISTRA)%PSIM
              EMELS(ISTRA)%PSIM => CPSIM
            ENDIF
            IF (EIEL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EIEL(IN)*FLXI
              CPSIM%NXTSIM => EIELS(ISTRA)%PSIM
              EIELS(ISTRA)%PSIM => CPSIM
            ENDIF
            IF (EAPL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EAPL(IN)*FLXI
              CPSIM%NXTSIM => EAPLS(ISTRA)%PSIM
              EAPLS(ISTRA)%PSIM => CPSIM
            ENDIF
            IF (EMPL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EMPL(IN)*FLXI
              CPSIM%NXTSIM => EMPLS(ISTRA)%PSIM
              EMPLS(ISTRA)%PSIM => CPSIM
            ENDIF
            IF (EIPL(IN) .NE. 0.D0) THEN
              ALLOCATE(CPSIM)
              CPSIM%ICS = IN
              CPSIM%VALUES = EIPL(IN)*FLXI
              CPSIM%NXTSIM => EIPLS(ISTRA)%PSIM
              EIPLS(ISTRA)%PSIM => CPSIM
            ENDIF
          ENDDO
          DO IATM=1,NATMI
!PB         DO IN=1,NDXY
            DO IN=1,NSBOX_TAL
              IF (PDENA(IATM,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IATM
                CPMUL%ICM = IN
                CPMUL%VALUEM = PDENA(IATM,IN)*FLXI
                CPMUL%NXTMUL => PDENAS(ISTRA)%PMUL
                PDENAS(ISTRA)%PMUL => CPMUL
              ENDIF
              IF (EDENA(IATM,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IATM
                CPMUL%ICM = IN
                CPMUL%VALUEM = EDENA(IATM,IN)*FLXI
                CPMUL%NXTMUL => EDENAS(ISTRA)%PMUL
                EDENAS(ISTRA)%PMUL => CPMUL
              ENDIF
            ENDDO
          ENDDO
          DO IMOL=1,NMOLI
!PB         DO IN=1,NDXY
            DO IN=1,NSBOX_TAL
              IF (PDENM(IMOL,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IMOL
                CPMUL%ICM = IN
                CPMUL%VALUEM = PDENM(IMOL,IN)*FLXI
                CPMUL%NXTMUL => PDENMS(ISTRA)%PMUL
                PDENMS(ISTRA)%PMUL => CPMUL
              ENDIF
            ENDDO
          ENDDO
          DO IION=1,NIONI
!PB         DO IN=1,NDXY
            DO IN=1,NSBOX_TAL
              IF (PDENI(IION,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = IION
                CPMUL%ICM = IN
                CPMUL%VALUEM = PDENI(IION,IN)*FLXI
                CPMUL%NXTMUL => PDENIS(ISTRA)%PMUL
                PDENIS(ISTRA)%PMUL => CPMUL
              ENDIF
            ENDDO
          ENDDO
          DO ICPV=1,NCPVI
!PB         DO IN=1,NDXY
            DO IN=1,NSBOX_TAL
              IF (COPV(ICPV,IN) .NE. 0.D0) THEN
                ALLOCATE(CPMUL)
                CPMUL%IART = ICPV
                CPMUL%ICM = IN
                CPMUL%VALUEM = COPV(ICPV,IN)*FLXI
                CPMUL%NXTMUL => COPVS(ISTRA)%PMUL
                COPVS(ISTRA)%PMUL => CPMUL
              ENDIF
            ENDDO
          ENDDO
C
C
60      CONTINUE
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        CALL INTER3(LSTOP,IFIRST,1,NSTRAI,0)
C
        IFIRST=IFIRST+1
        RETURN
C
10000   CONTINUE
C
        CALL INTER1
C
        CALL PLASMA
C
        CALL PLASMA_DERIV(0)
C
        CALL SETAMD
C
        IF (.NOT.ALLOCATED(SPLNWA)) CALL ALLOC_BRASCL
        CALL INIT_BRASCL2
C
C  NEW: ATOMS, EI PROCESSES
C
        DO 101 IATM=1,NATMI
        DO 101 IPLS=1,NPLSI
          DO 102 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 101
            DO 102 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
              ELSE
                SPLNWA(IN,IATM,IPLS)=SPLNWA(IN,IATM,IPLS)+
     .                          FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
              END IF
102       CONTINUE
101     CONTINUE
C
        DO 103 IPLS=1,NPLSI
          IPLSTI = MPLSTI(IPLS)
          DO 104 IN=1,NDXY
            SEINWA(IN,IPLS)=DIIN(IPLS,IN)*
     .                      (1.5*TIIN(IPLSTI,IN)+EDRIFT(IPLS,IN))
104       CONTINUE
103     CONTINUE
C
        DO 105 IATM=1,NATMI
          DO 105 IAEI=1,NAEII(IATM)
            IRDS=LGAEI(IATM,IAEI)
            DO 105 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWA(IN,IATM)=SEENWA(IN,IATM)+FEELEI1(IRDS,IN)*
     .                                          FTABEI1(IRDS,IN)
              END IF
105     CONTINUE
C
C  NEW: TEST IONS, EI PROCESSES
C
        DO 106 IION=1,NIONI
          DO 106 IIDS=1,NIDSI(IION)
            IRDS=LGIEI(IION,IIDS)
            DO 106 IN=1,NDXY
              IF (NSTORDR >= NRAD) THEN
                SEENWI(IN,IION)=SEENWI(IN,IION)+EELDS1(IRDS,IN)*
     .                                          TABDS1(IRDS,IN)
              ELSE
                SEENWI(IN,IION)=SEENWI(IN,IION)+FEELEI1(IRDS,IN)*
     .                                          FTABEI1(IRDS,IN)
              END IF
106     CONTINUE
C
        DO 107 IION=1,NIONI
        DO 107 IPLS=1,NPLSI
        DO 107 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 107
          DO 108 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWI(IN,IION,IPLS)=SPLNWI(IN,IION,IPLS)+
     .                             FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
108       CONTINUE
107     CONTINUE
C
        DO 109 IION=1,NIONI
        DO 109 IPLS=1,NPLSI
        DO 109 IIDS=1,NIDSI(IION)
          IRDS=LGIEI(IION,IIDS)
          ESIG=EPLDS(IRDS,2)
          DO 110 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWI(IN,IION)=SEINWI(IN,IION)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWI(IN,IION)=SEINWI(IN,IION)+FTABEI1(IRDS,IN)*ESIG
            END IF
110       CONTINUE
109     CONTINUE
C
C  NEW: MOLECULES, EI PROCESSES
C
        DO 115 IMOL=1,NMOLI
        DO 115 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          DO 116 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+EELDS1(IRDS,IN)*
     .                                        TABDS1(IRDS,IN)
            ELSE
              SEENWM(IN,IMOL)=SEENWM(IN,IMOL)+FEELEI1(IRDS,IN)*
     .                                        FTABEI1(IRDS,IN)
            END IF
116       CONTINUE
115     CONTINUE
C
        DO 117 IMOL=1,NMOLI
        DO 117 IPLS=1,NPLSI
        DO 117 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          IF (PPLDS(IRDS,IPLS).EQ.0.) GOTO 117
          DO 118 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                             TABDS1(IRDS,IN)*PPLDS(IRDS,IPLS)
            ELSE
              SPLNWM(IN,IMOL,IPLS)=SPLNWM(IN,IMOL,IPLS)+
     .                             FTABEI1(IRDS,IN)*PPLDS(IRDS,IPLS)
            END IF
118       CONTINUE
117     CONTINUE
C
        DO 119 IMOL=1,NMOLI
        DO 119 IPLS=1,NPLSI
        DO 119 IMDS=1,NMDSI(IMOL)
          IRDS=LGMEI(IMOL,IMDS)
          ESIG=EPLDS(IRDS,2)
          DO 120 IN=1,NDXY
            IF (NSTORDR >= NRAD) THEN
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+TABDS1(IRDS,IN)*ESIG
            ELSE
              SEINWM(IN,IMOL)=SEINWM(IN,IMOL)+FTABEI1(IRDS,IN)*ESIG
            END IF
120       CONTINUE
119     CONTINUE
C
        B2BREM=B2BRM
        B2RAD=B2RD
        B2QIE=B2Q
        B2VDP=B2VP
        CALL INTER3(LSTOP,IFIRST,1,NSTRAI,0)
C
        IFIRST=IFIRST+1

        IF (LSTOP) THEN
          CALL DEALLOC_COMUSR
          CALL DEALLOC_CESTIM
          CALL DEALLOC_BRASCL
          CALL DEALLOC_BRASPOI
        END IF

        RETURN
C
      ENDIF

      END
