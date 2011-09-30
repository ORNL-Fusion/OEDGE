cdr  111107: "istep out of range" error message removed once again.
!pb  220307: LEVGEO=6 --> LEVGEO=10
!pb  271006: use flux set by user defined sampling routine 
!pb  100106: ENTRY SAMVOL_REINIT added for reinitialsation of Eirene
!pb  181206: calculate bremsstrahlung
!pb  240806: set output values for DIWL and SHWL
cdr  200806: tiwl(*), ... instead of tiwl(npls),... to unify code.
cdr  060406: check "istep out of range" moved to correct place
c  031105
c  iplsti moved after check of validity of ipls, to produce legal exit
c         rather than code crash
C  JET 2005, PATCH 1: NEW ARGUMENTS EFWL AND SHWL IN PARAMETER LIST
c                     AT ENTRY SMVOL1 AND SMUSR1
C
      SUBROUTINE SAMVOL

      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CCONA
      USE CLOGAU
      USE CINIT
      USE CPOLYG
      USE CGRID
      USE CZT1
      USE CTRCEI
      USE CGEOM
      USE COMPRT
      USE COMSOU
      USE COMXS
      USE CSPEI
      USE CTRIG
      USE CTETRA
      USE PHOTON
      IMPLICIT NONE
C

      REAL(DP), INTENT(OUT) :: TEWL, SHWL, TIWL(*), DIWL(*),
     .                         VXWL(*), VYWL(*), VZWL(*),
     .                         EFWL(*), WEISPZ(*)
      INTEGER, INTENT(IN) :: NVLM
      REAL(DP), ALLOCATABLE, SAVE :: FREC(:,:,:), VSOURC(:,:), VSMXI(:)
      REAL(DP), ALLOCATABLE, SAVE :: RQ21(:), PS21(:)
      REAL(DP), ALLOCATABLE, SAVE :: ASIMP(:,:)
      REAL(DP) :: D(3)
      INTEGER, ALLOCATABLE, SAVE  :: ISOURC(:,:), ICMX(:),
     .                               IFREC(:)
      REAL(DP) :: ZEP1, X1, Y1, X2, Y2, X3, Y3, RR, RRI, WINK,
     .            ZRM1, CNORM, EPR, ELR, RRD, RRN, ZZ, X01, Z1, Z2, Z3,
     .            REC, BX, BY, BZ, ADD, FTABRC1, CNDYNP, VX, VY, VZ,
     .            VPARA, EELRC, FEELRC1, SUMM, EISUMM, EISUM, SUM,
     .            X4, Y4, Z4, MOMPARA, BREMS, TOT_BREMS(NPLS), Z, ngffmh
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER :: IC1, IC2, ICOUNT, ICELL, IAUSR, IBUSR, IRUSR, IPUSR,
     .           ITUSR, IN, IIRC, IRC, IRRC, J, IT1, IT2, ISTEP, IFRC,
     .           IR2, IP1, IP2, ICTOT, IND, IR, IP, IT, IC, ISRFSI, I,
     .           ICC, IR1, ISR, ISTR, IL, IU, IM, MXREC, MXPLS, IFPLS,
     .           IPLSTI, IPLSV, KK, ICCT
      INTEGER, SAVE :: ISTROLD=-1
      LOGICAL, ALLOCATABLE, SAVE :: LPLSSR(:)
C
C  AT ENTRY SAMVL0:
C    DEFINE THE CUMULATIVE DISTRIBUTION FUNCTION
C    FREC(IPLS,IRRC,ICELL) FOR EACH VOLUME SOURCE DISTRIBUTION, FOR SAMPLING
C    THE CELL INDEX ICELL OF THE VOLUME SOURCE PARTICLE.
C
C    A FEW GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING
C    OF THE STARTING POINT IN EACH CELL ARE COMPUTED
C
C    THE SOURCE STRENGTH FLUX(ISTRA) IS MODIFIED FOR THE
C    STRATA WITH NLVOL(ISTRA)=.TRUE.
C
C  AT ENTRY SAMVL1:
C    THE INITIAL CO-ORDINATES OF A TEST FLIGHT ARE SAMPLED,
C    AND THE CELL NUMBERS ARE COMPUTED
C
      ENTRY SAMVL0
C

      IF (.NOT.ALLOCATED(FREC)) THEN

C  LPLSSR(IPLS):
C  IDENTIFY THOSE IPLS WHICH NEED A VOLUME SOURCE DISTRIBUTION

        ALLOCATE (LPLSSR(NPLSI))
        LPLSSR = .FALSE.
        DO ISTR=1,NSTRAI
          IF (NLVOL(ISTR) .AND. NLPLS(ISTR) .AND. (NPTS(ISTR) > 0)
     .        .AND. (FLUX(ISTR) > 0._DP)) THEN
            IPLS = NSPEZ(ISTR)
            IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
              LPLSSR = .TRUE.
            ELSE
              LPLSSR(IPLS) = .TRUE.
            END IF
          END IF
        END DO
        MXREC=MAXVAL(NPRCI(1:NPLSI))
        MXPLS=COUNT(LPLSSR(1:NPLSI))
        ALLOCATE (FREC(0:MXPLS,0:MXREC,0:NRAD))
        ALLOCATE (VSOURC(NSRFS,0:NRAD))
        ALLOCATE (VSMXI(NSRFS))
        ALLOCATE (RQ21(N1ST))
        ALLOCATE (PS21(N2ND))
        IF (LEVGEO == 3) ALLOCATE (ASIMP(2,NRAD))
        ALLOCATE (ISOURC(NSRFS,0:NRAD))
        ALLOCATE (ICMX(NSRFS))
        ALLOCATE (IFREC(NPLS))
      END IF
C
      FREC=0.
      SREC=0.
      IFREC=0
C
      IFPLS=0
      DO 2 IPLS=1,NPLSI
        IF (LGPRC(IPLS,0).EQ.0) GOTO 2
        IF (.NOT.LPLSSR(IPLS)) GOTO 2
        IFPLS=IFPLS+1
        IFREC(IPLS)=IFPLS
        DO 3 IIRC=1,NPRCI(IPLS)
          IRRC=LGPRC(IPLS,IIRC)
          KK=NREARC(IRRC)
          ICCT=NREACT(KK)
          DO 3 J=1,NSBOX
            ADD=0.
C  EXCLUDE DEAD CELLS (GRID CUTS, ISOLATED CELLS FROM COUPLE_.., ETC)
C  EXCLUDE IPLS-VACUUM CELLS
            IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
              IF (NSTORDR >= NRAD) THEN
                ADD=TABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              ELSE
                ADD=FTABRC1(IRRC,J)*DIIN(IPLS,J)*VOL(J)*ELCHA
              END IF
            END IF
            IF (ICCT > 0) 
     .        ADD = ADD*(XINTLEFT(ICCT,J) + 
     .                   XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))
            FREC(IFPLS,IIRC,J)  =FREC(IFPLS,IIRC,J-1)+ADD
            SREC(IPLS,IRRC)     =SREC(IPLS,IRRC)+ADD
3       CONTINUE
2     CONTINUE

C  SUM OVER SPECIES AND RECOMBINATION TYP INDICES
      DO 4 IPLS=1,NPLSI
        IF (LGPRC(IPLS,0).EQ.0) GOTO 4
        IF (.NOT.LPLSSR(IPLS)) GOTO 4
        IFPLS=IFREC(IPLS)
        DO 5 IIRC=1,NPRCI(IPLS)
          IRRC=LGPRC(IPLS,IIRC)
          SREC(IPLS,0)=SREC(IPLS,0)+SREC(IPLS,IRRC)
          SREC(0,IRRC)=SREC(0,IRRC)+SREC(IPLS,IRRC)
          SREC(0,0)   =SREC(0,0)   +SREC(IPLS,IRRC)
          DO 5 J=1,NSBOX
            FREC(IFPLS,0,J)=FREC(IFPLS,0,J)+FREC(IFPLS,IIRC,J)
5         CONTINUE
4     CONTINUE
C
C
      IF (TRCSOU.AND.IFPLS.GT.0) THEN
        EIO=0.
        EEL=0.
        MOM=0.
C
        DO 7 IPLS=1,NPLSI
          CNDYNP=AMUA*RMASSP(IPLS)
          IF (LGPRC(IPLS,0).EQ.0) GOTO 7
          IF (.NOT.LPLSSR(IPLS)) GOTO 7
          IFPLS=IFREC(IPLS)
          IPLSTI = MPLSTI(IPLS)
          IPLSV = MPLSV(IPLS)
          DO 6 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)
            ICCT=NREACT(KK)
            DO 6 J=1,NSBOX
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
                REC=FREC(IFPLS,IIRC,J)-FREC(IFPLS,IIRC,J-1)
                IF (REC.LE.0.D0) GOTO 6
                ADD=(1.5*TIIN(IPLSTI,J)+EDRIFT(IPLS,J))*REC
                IF (ICCT > 0) 
     .            ADD = ADD*(XINTLEFT(ICCT,J) + 
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))
                EIO(IPLS,IRRC)=EIO(IPLS,IRRC)-ADD
                EIO(IPLS,0)   =EIO(IPLS,0   )-ADD
                IF (INDPRO(5) == 8) THEN
                  CALL VECUSR(1,BX,BY,BZ,1)
                ELSE
                  BX=BXIN(J)
                  BY=BYIN(J)
                  BZ=BZIN(J)
                END IF
                IF (INDPRO(4) == 8) THEN
                  CALL VECUSR(2,VX,VY,VZ,IPLSV)
                  VPARA=VX*BX+VY*BY+VZ*BZ
                  MOMPARA=VPARA*CNDYNP*SIGN(1._DP,VPARA)
                ELSE IF (INDPRO(5) == 8) THEN
                  VX = VXIN(IPLSV,J)
                  VY = VYIN(IPLSV,J)
                  VZ = VZIN(IPLSV,J)
                  VPARA=VX*BX+VY*BY+VZ*BZ
                  MOMPARA=VPARA*CNDYNP*SIGN(1._DP,VPARA)
                ELSE
                  MOMPARA=PARMOM(IPLS,J)
                ENDIF
                ADD=MOMPARA*REC
                IF (ICCT > 0) 
     .            ADD = ADD*(XINTLEFT(ICCT,J) + 
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))
                MOM(IPLS,IRRC)=MOM(IPLS,IRRC)-ADD
                MOM(IPLS,0)   =MOM(IPLS,0   )-ADD
              ENDIF
6         CONTINUE
C
          DO 8 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            KK=NREARC(IRRC)
            ICCT=NREACT(KK)
            DO 8 J=1,NSBOX
              ADD=0.D0
              IF (NSTGRD(J).EQ.0.AND..NOT.LGVAC(J,IPLS)) THEN
                IF (NSTORDR >= NRAD) THEN
                  EELRC = EELRC1(IRRC,J)
                ELSE
                  EELRC = FEELRC1(IRRC,J)
                END IF
                ADD=EELRC*DIIN(IPLS,J)*VOL(J)*ELCHA
                IF (ICCT > 0) 
     .            ADD = ADD*(XINTLEFT(ICCT,J) + 
     .                       XINT_INF(ICCT,J) - XINTRIGHT(ICCT,J))
              ENDIF
              EEL(IPLS,IRRC)=EEL(IPLS,IRRC)+ADD
              EEL(IPLS,0   )=EEL(IPLS,0   )+ADD
8         CONTINUE
7       CONTINUE

        TOT_BREMS = 0._DP
        DO IPLS=1,NPLSI
          IF (NCHRGP(IPLS) == 0) CYCLE
          Z = NCHRGP(IPLS)
          DO J = 1, NSBOX
            IF (LGVAC(J,NPLS+1).OR.LGVAC(J,IPLS)) CYCLE
            BREMS = 1.54E-32_DP * TEIN(J)**0.5 * Z**2 *
     .              ngffmh(Z**2 * 13.6_DP/TEIN(J))
            TOT_BREMS(IPLS) = TOT_BREMS(IPLS) + 
     .                        BREMS*DEIN(J)*DIIN(IPLS,J)*VOL(J)
          END DO
        END DO
C
        CALL LEER(1)
        WRITE (iunout,*) 'DIAGNOSTICS FROM SUBR. SAMVOL: '
        CALL LEER(1)
        WRITE (iunout,*) 'VOLUME RECOMBINATION RATES INTEGRATED OVER'
        WRITE (iunout,*) 'ENTIRE COMPUTATIONAL GRID '
        CALL LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION PARTICLE LOSS (AMP): '
        ITYP=4
        DO 10 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 11 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL MASAJR('IPLS,IRRC, SREC         ',
     .                   TEXTS(ISPZ),IRRC,-SREC(IPLS,IRRC))
11        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL MASAJR('IPLS,TOT., SREC(IPLS,0) ',
     .                   TEXTS(ISPZ),0   ,-SREC(IPLS,0))
          ENDIF
10      CONTINUE
        CALL LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ION ENERGY LOSS (WATT): '
        DO 12 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 13 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL MASAJR('IPLS,IRRC,EIO           ',
     .                   TEXTS(ISPZ),IRRC,EIO(IPLS,IRRC))
13        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL MASAJR('IPLS,TOT.,EIO(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EIO(IPLS,0))
          ENDIF
12      CONTINUE
        CALL LEER(1)
        WRITE (iunout,*) 'RECOMBINATION ELECTRON ENERGY LOSS (WATT): '
        DO 14 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 15 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL MASAJR('IPLS,IRRC,EEL           ',
     .                   TEXTS(ISPZ),IRRC,EEL(IPLS,IRRC))
15        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL MASAJR('IPLS,TOT.,EEL(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,EEL(IPLS,0))
          ENDIF
14      CONTINUE
        CALL LEER(1)
        WRITE (iunout,*) 'RECOMBINATION PARALLEL MOMENTUM LOSS : '
        DO 16 IPLS=1,NPLSI
          IF (.NOT.LPLSSR(IPLS)) CYCLE
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          DO 17 IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            CALL MASAJR('IPLS,IRRC,MOM           ',
     .                   TEXTS(ISPZ),IRRC,MOM(IPLS,IRRC))
17        CONTINUE
          IF (NPRCI(IPLS).GT.1) THEN
            CALL MASAJR('IPLS,TOT.,MOM(IPLS,0)   ',
     .                   TEXTS(ISPZ),0   ,MOM(IPLS,0))
          ENDIF
16      CONTINUE
        CALL LEER(1)

        WRITE (iunout,*) 'BREMSSTRAHLUNG: '
        DO IPLS=1,NPLSI
          ISPZ=ISPEZ(ITYP,IPHOT,IATM,IMOL,IION,IPLS)
          CALL MASAJR('IPLS,TOT.BREMSSTRAHLUNG ',
     .                 TEXTS(ISPZ),0   ,TOT_BREMS(IPLS))          
        END DO

      ENDIF    !trcsou
C
C  SET TOTAL SOURCE STRENGTH FOR STRATA WITH NLVOL(ISTRA)=.TRUE.,
C
      DO 50 ISTRA=1,NSTRAI
        IF (NLVOL(ISTRA).AND.NLPLS(ISTRA)) THEN
          IPLS=NSPEZ(ISTRA)
          IF (IPLS.LE.0.OR.IPLS.GT.NPLSI) THEN
            WRITE (iunout,*) 'SOURCE SPECIES INDEX NSPEZ OUT OF RANGE'
            WRITE (iunout,*) 'ISTRA, NSPEZ(ISTRA) ',ISTRA,NSPEZ(ISTRA)
            CALL EXIT_OWN(1)
          ENDIF
          IPLSTI = MPLSTI(IPLS)
          SUMM=0.D0
          EISUMM=0.D0
          DO 53 ISRFSI=1,NSRFSI(ISTRA)
            ISR=ISRFSI
            SUM=0.D0
            EISUM=0.D0
            IF (SORLIM(ISR,ISTRA).LT.0) THEN
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR USER SPECIFIED VOLUME SOURCE
              CALL SM0USR(ISR,ISTRA,
     .                    SORAD1(ISR,ISTRA),SORAD2(ISR,ISTRA),
     .                    SORAD3(ISR,ISTRA),SORAD4(ISR,ISTRA),
     .                    SORAD5(ISR,ISTRA),SORAD6(ISR,ISTRA))
!pb assume flux is set in samusr
              SUMM=FLUX(ISTRA)
            ELSE
C  INITIALIZE SAMPLING DISTRIBUTIONS FOR DEFAULT VOLUME RECOMBINATION SOURCES
C  ACCOUNT FOR INGRDA(ISRFSI,ISTRA,...), INGRDE(ISRFSI,ISTRA,...)
              I=ISTRA
              ICC=0
              IRC=-1
              IF (NR1ST.GT.1) THEN
              IF (INGRDA(ISR,I,1).LE.0..OR.INGRDE(ISR,I,1).LE.0.D0) THEN
                CALL LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*) 
     .            'NEW INPUT FOR INGRDA(.,.,1),INGRDE(.,.,1)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(ISR,I,1)=1
                INGRDE(ISR,I,1)=MAX0(1,NR1ST)
                CALL LEER(1)
              ENDIF
              ENDIF
              IF (NP2ND.GT.1) THEN
              IF (INGRDA(ISR,I,2).LE.0..OR.INGRDE(ISR,I,2).LE.0.D0) THEN
                CALL LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*) 
     .            'NEW INPUT FOR INGRDA(.,.,2),INGRDE(.,.,2)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(ISR,I,2)=1
                INGRDE(ISR,I,2)=MAX0(1,NP2ND)
                CALL LEER(1)
              ENDIF
              ENDIF
              IF (NT3RD.GT.1) THEN
              IF (INGRDA(ISR,I,3).LE.0..OR.INGRDE(ISR,I,3).LE.0.D0) THEN
                CALL LEER(1)
                WRITE (iunout,*) 'WARNING FROM SAMVL0, ISTRA= ',ISTRA
                WRITE (iunout,*) 
     .            'NEW INPUT FOR INGRDA(.,.,3),INGRDE(.,.,3)'
                WRITE (iunout,*) 'AUTOMATIC CORRECTION CARRIED OUT '
                INGRDA(ISR,I,3)=1
                INGRDE(ISR,I,3)=MAX0(1,NT3RD)
                CALL LEER(1)
              ENDIF
              ENDIF
              IF (NPRCI(IPLS).EQ.0) THEN
                WRITE (iunout,*) 'NO DEFAULT VOLUME SOURCE DISTRIBUTION'
                WRITE (iunout,*) 'DEFINED. SUBSTRATUM TURNED OFF'
                WRITE (iunout,*) 'IPLS,ISRFSI,ISTRA ',IPLS,ISRFSI,ISTRA
                SORWGT(ISR,ISTRA)=0.D0
                GOTO 53
              ENDIF
              IF (NLRAD) THEN
                IR1=MAX0(1,INGRDA(ISR,ISTRA,1))
                IR2=MIN0(NR1ST,INGRDE(ISR,ISTRA,1))
              ELSE
                IR1=1
                IR2=2
              ENDIF
              IF (NLPOL) THEN
                IP1=MAX0(1,INGRDA(ISR,ISTRA,2))
                IP2=MIN0(NP2ND,INGRDE(ISR,ISTRA,2))
              ELSE
                IP1=1
                IP2=2
              ENDIF
              IF (NLTOR) THEN
                IT1=MAX0(1,INGRDA(ISR,ISTRA,3))
                IT2=MIN0(NT3RD,INGRDE(ISR,ISTRA,3))
              ELSE
                IT1=1
                IT2=2
              ENDIF

              ISTEP=SORIND(ISR,ISTRA)
              IFPLS=IFREC(IPLS)
              DO 52 IIRC=1,NPRCI(IPLS)
                IRRC=LGPRC(IPLS,IIRC)
                IF (ISTEP.EQ.IRRC) THEN
C  RECOMBINATION PROCESS IRRC IDENTIFIED
                  IRC=IRRC
                  IFRC=IIRC
                ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL RECOMBINATION PROCESSES FOR SPECIES IPLS
                  IRC=0
                  IFRC=0
                  ISTEP=-1
                ELSE
C  TRY OTHER RECOMBINATION PROCESS ASSIGNED TO IPLS
                  GOTO 52
                ENDIF
                DO 51 IR=IR1,IR2-1
                  DO 51 IP=IP1,IP2-1
                    DO 51 IT=IT1,IT2-1
                      NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                      ADD=FREC(IFPLS,IFRC,NCELL)-
     .                    FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                      IF (ADD.GT.0.D0) THEN
                        ICC=ICC+1
                        SUM=SUM+ADD
                        EISUM=EISUM-
     .                   (1.5*TIIN(IPLSTI,NCELL)+EDRIFT(IPLS,NCELL))*ADD
                      ENDIF
51              CONTINUE
52            CONTINUE
              IF (SUM.EQ.0.D0) THEN
                WRITE (IUNOUT,*) 'NO VOL. RECOMBINATION SOURCE FOR: '
                WRITE (IUNOUT,*) 'ISTRA, ISRFSI, IPLS, ISTEP ',
     .                            ISTRA, ISR   , IPLS, ISTEP
                WRITE (IUNOUT,*) 'EITHER: ISTEP OUT OF RANGE IN SAMVOL'
                WRITE (IUNOUT,*) 'OR:  DENSITY OF RECOMBINING IPLS = 0 '
                SORWGT(ISR,ISTRA)=0.D0
                GOTO 53
              ENDIF
              SORWGT(ISR,ISTRA)=SUM
              CALL LEER(1)
              WRITE (iunout,*) 'SUB-STRATUM WEIGHT REDEFINED '
              CALL MASJ2R ('ISRFSI,ISTRA,SORWGT     ',ISRFSI,ISTRA,SUM)
              IF (TRCSOU) THEN
                CALL MASJ3 ('IRRC,IPLS,ICMX          ',
     .                       IRC ,IPLS,ICC)
                CALL LEER(1)
              ENDIF
              SUMM=SUMM+SUM
              EISUMM=EISUMM+EISUM
            ENDIF
53        CONTINUE
C
          IF (SUMM.GT.0.D0) THEN
            FLUX(ISTRA)=SUMM
            WRITE (iunout,*) 'SOURCE STRENGTH REDEFINED '
            CALL MASJR2('ISTRA, FLUX, EIFLUX     ',
     .                   ISTRA,FLUX(ISTRA),EISUMM)
            CALL LEER(1)
          ELSE 
            FLUX(ISTRA)=0.D0
            WRITE (iunout,*) 'SOURCE ISTRA= ',ISTRA,' TURNED OFF '
            CALL LEER(1)
          ENDIF
        ENDIF
50    CONTINUE
C
C  PREPARE SOME GEOMETRICAL CONSTANTS FOR RANDOM SAMPLING IN STANDARD MESH CELLS
      IF (LEVGEO.EQ.2) THEN
        IF (NLPOL) THEN
          DO 54 IP=1,NP2NDM
            PS21(IP)=PSURF(IP+1)-PSURF(IP)
54        CONTINUE
        ENDIF
        DO 55 IR=1,NR1STM
          RQ21(IR)=RQ(IR+1)-RQ(IR)
55      CONTINUE
      ENDIF
C
      IF (LEVGEO.EQ.3) THEN
        IT=1
        DO 56 IR=1,NR1ST-1
        DO 56 IP=1,NP2ND-1
          IND=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
          X1=XPOL(IR,IP)
          X2=XPOL(IR,IP+1)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR,IP)
          Y2=YPOL(IR,IP+1)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(1,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
          X1=XPOL(IR+1,IP)
          X2=XPOL(IR,IP)
          X3=XPOL(IR+1,IP+1)
          Y1=YPOL(IR+1,IP)
          Y2=YPOL(IR,IP)
          Y3=YPOL(IR+1,IP+1)
          ASIMP(2,IND)=0.5*(X1*(Y2-Y3)+X2*(Y3-Y1)+X3*(Y1-Y2))
56      CONTINUE
      ENDIF
C
      RETURN
C
C  AT THIS POINT: CALLED FROM PARTICLE LOOP TO INITIALIZE TEST FLIGHT
C
      ENTRY SAMVL1(NVLM,TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
C  USER SUPPLIED SOURCE
C
      IF (SORLIM(NVLM,ISTRA).LT.0) THEN
        CALL SM1USR(NVLM,X0,Y0,Z0,
     .              SORAD1(NVLM,ISTRA),SORAD2(NVLM,ISTRA),
     .              SORAD3(NVLM,ISTRA),SORAD4(NVLM,ISTRA),
     .              SORAD5(NVLM,ISTRA),SORAD6(NVLM,ISTRA),
     .              IRUSR,IPUSR,ITUSR,IAUSR,IBUSR,
     .              TIWL,TEWL,DIWL,VXWL,VYWL,VZWL,EFWL,SHWL,WEISPZ)
        NRCELL=IRUSR
        NPCELL=IPUSR
        NTCELL=ITUSR
        NACELL=IAUSR
        NBLOCK=IBUSR
        NBLCKA=NSTRD*(NBLOCK-1)+NACELL
        NCELL=NRCELL+((NPCELL-1)+(NTCELL-1)*NP2T3)*NR1P2+NBLCKA
C
        MTSURF=0
        NLSRFZ=.FALSE.
        MPSURF=0
        NLSRFY=.FALSE.
        MRSURF=0
        NLSRFX=.FALSE.
        EFWL(1:NPLSI)=0._DP
        RETURN
      ENDIF
C
C  VOLUME RECOMBINATION SOURCE
C
C  TENTATIVELY ASSUME: A BULK ION WILL BE GENERATED
      LGPART=.TRUE.
      ITYP=4
C
      IF (.NOT.NLPLS(ISTRA)) GOTO 999

      IF (ISTROLD /= ISTRA) THEN
        ISTROLD=ISTRA
        IPLS=NSPEZ(ISTRA)
        DO ISRFSI=1,NSRFSI(ISTRA)
          ISR=ISRFSI
          ICC=0
          VSOURC(ISR,0)=0.D0
          IF (NLRAD) THEN
            IR1=MAX0(1,INGRDA(ISR,ISTRA,1))
            IR2=MIN0(NR1ST,INGRDE(ISR,ISTRA,1))
          ELSE
            IR1=1
            IR2=2
          ENDIF
          IF (NLPOL) THEN
            IP1=MAX0(1,INGRDA(ISR,ISTRA,2))
            IP2=MIN0(NP2ND,INGRDE(ISR,ISTRA,2))
          ELSE
            IP1=1
            IP2=2
          ENDIF
          IF (NLTOR) THEN
            IT1=MAX0(1,INGRDA(ISR,ISTRA,3))
            IT2=MIN0(NT3RD,INGRDE(ISR,ISTRA,3))
          ELSE
            IT1=1
            IT2=2
          ENDIF

          ISTEP=SORIND(ISR,ISTRA)
          IFPLS=IFREC(IPLS)
          DO IIRC=1,NPRCI(IPLS)
            IRRC=LGPRC(IPLS,IIRC)
            IF (ISTEP.EQ.IRRC) THEN
              IRC=IRRC
              IFRC=IIRC
            ELSEIF (ISTEP.EQ.0) THEN
C  SUM OVER ALL RECOMBINATION PROCESSES FOR SPECIES IPLS
              IRC=0
              IFRC=0
              ISTEP=-1
            ELSE
              CYCLE
            ENDIF
            DO IR=IR1,IR2-1
              DO IP=IP1,IP2-1
                DO IT=IT1,IT2-1
                  NCELL=IR+((IP-1)+(IT-1)*NP2T3)*NR1P2
                  ADD=FREC(IFPLS,IFRC,NCELL)-
     .                FREC(IFPLS,IFRC,NCELL-1)
C  INDIRECT ADDRESSING
                  IF (ADD.GT.0.D0) THEN
                    ICC=ICC+1
                    ISOURC(ISR,ICC)=NCELL
                    VSOURC(ISR,ICC)=VSOURC(ISR,ICC-1)+ADD
                  ENDIF
                END DO
              END DO
            END DO
          END DO ! IIRC
          ICMX(ISR)=ICC
          VSMXI(ISR) = 1._DP / VSOURC(ISR,ICC)
        END DO ! ISRFSI
C
      END IF
C
C  FIND CELL NUMBER: NCELL
C
      IF (INDIM(NVLM,ISTRA) .GE. 0) THEN
!PB  choose cell according to cell contribution to total source strength 
        IC1=0
        IC2=ICMX(NVLM)
        ZEP1=RANF_EIRENE()*VSOURC(NVLM,IC2)

        IL=0
        IU=IC2

c  binary search
        DO WHILE (IU-IL.gt.1)
          IM=(IU+IL)*0.5
          IF(ZEP1.GE.VSOURC(NVLM,IM)) THEN
            IL=IM
          ELSE
            IU=IM
          ENDIF
        END DO
c
        ICELL=IU

        NCELL=ISOURC(NVLM,ICELL)
      ELSE

!pb  uniform distribution 
        IC1=0
        IC2=ICMX(NVLM)

        ICELL = MIN(INT(1+RANF_EIRENE()*(IC2-1)),IC2)
        NCELL = ISOURC(NVLM,ICELL)

        WEIGHT=(VSOURC(NVLM,ICELL)-VSOURC(NVLM,ICELL-1))*VSMXI(NVLM)*IC2

      END IF
C     
      IF (NCELL.GT.NSURF) GOTO 991
C
C  A CELL NUMBER NCELL HAS NOW BEEN SAMPLED
C
      CALL NCELLN(NCELL,NRCELL,NPCELL,NTCELL,NACELL,NBLOCK,
     .            NR1ST,NP2ND,NT3RD,NBMLT,NLRAD,NLPOL,NLTOR)
C
C  FIND TOROIDAL CO-ORDINATE IN NTCELL
C
      IF (.NOT.NLTOR) THEN
C       NTCELL=1
        IPERID=1
        IF (NLTRZ) THEN
          Z0=0.
        ELSEIF (NLTRA) THEN
C  TACTICALLY ASSUME: PARTICLE STARTS IN LOCAL TOR. BASIS CELL NO.1
          ZRM1=ZSURF(1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
          IPERID=1
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
C         IPERID=LEARCA(PHI,ZSURF,1,NTTRA,1,'SAMVOL      ')
        ELSEIF (NLTRT) THEN
          GOTO 999
        ENDIF
      ELSEIF (NLTOR) THEN
        IPERID=NTCELL
C  SAMPLE IN CELL NTCELL
        IF (NLTRZ) THEN
          Z0=ZSURF(NTCELL)+RANF_EIRENE()*(ZSURF(NTCELL+1)-ZSURF(NTCELL))
        ELSEIF (NLTRT) THEN
          PHI=ZSURF(NTCELL)+RANF_EIRENE()*
     .        (ZSURF(NTCELL+1)-ZSURF(NTCELL))
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ELSEIF (NLTRA) THEN
          ZRM1=ZFULL*(NTCELL-1)
          PHI=ZRM1+RANF_EIRENE()*ZFULL
C         Z0=??, TO BE FOUND FROM X01,PHI LATER
        ENDIF
      ENDIF
C
C  FIND RADIAL AND POLOIDAL CO-ORDINATE
C
      IF (LEVGEO.EQ.1) THEN
        X0=RSURF(NRCELL)+RANF_EIRENE()*(RSURF(NRCELL+1)-RSURF(NRCELL))
        IF (NLPOL) THEN
          Y0=PSURF(NPCELL)+RANF_EIRENE()*(PSURF(NPCELL+1)-PSURF(NPCELL))
        ELSE
          Y0=YIA+RANF_EIRENE()*(YAA-YIA)
        END IF
      ELSEIF (LEVGEO.EQ.2) THEN
        IF (NLCRC) THEN
C  POLOIDAL CO-ORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL CO-ORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          X0=RR*COS(WINK)
          Y0=RR*SIN(WINK)
        ELSEIF (NLELL) THEN
CDR NOT READY. STRICKLY, THETA AND R ARE CORRELATED. USE
CDR            MARGINAL AND CONDITIONAL DISTRIBUTION F1(R) AND
CDR            F2(PHI, GIVEN R)
C  POLOIDAL CO-ORDINATE
          IF (NLPOL) THEN
            WINK=PSURF(NPCELL)+RANF_EIRENE( )*PS21(NPCELL)
          ELSEIF (.NOT.NLPOL) THEN
            WINK=RANF_EIRENE( )*PI2A
          ENDIF
C  RADIAL CO-ORDINATE
          RR=SQRT(RQ(NRCELL)+RANF_EIRENE( )*RQ21(NRCELL))
C
          RRI=RSURF(NRCELL)
          RRD=RSURF(NRCELL+1)-RRI
          RRN=(RR-RRI)/RRD
C
          ELR=ELL(NRCELL)+RRN*(ELL(NRCELL+1)-ELL(NRCELL))
          EPR=EP1(NRCELL)+RRN*(EP1(NRCELL+1)-EP1(NRCELL))
          X0=RR*COS(WINK)+EPR
          Y0=RR*SIN(WINK)*ELR
        ELSEIF (NLTRI) THEN
          GOTO 999
        ENDIF
      ELSEIF (LEVGEO.EQ.3.OR.LEVGEO.EQ.4) THEN
        IF (LEVGEO.EQ.3) THEN
          IF (.NOT.NLPOL) THEN
            GOTO 999
          ENDIF
          IN = NRCELL + (NPCELL-1)*NR1ST
          ZEP1=AREA(IN)*RANF_EIRENE()
!pb          IF (ZEP1.LE.ASIMP(1,NCELL)) THEN
          IF (ZEP1.LE.ASIMP(1,IN)) THEN
C   PUNKT IN DREIECK 1
            X1=XPOL(NRCELL,NPCELL)
            X2=XPOL(NRCELL,NPCELL+1)
            X3=XPOL(NRCELL+1,NPCELL+1)
            Y1=YPOL(NRCELL,NPCELL)
            Y2=YPOL(NRCELL,NPCELL+1)
            Y3=YPOL(NRCELL+1,NPCELL+1)
          ELSE
C   PUNKT IN DREIECK 2
            X1=XPOL(NRCELL+1,NPCELL)
            X2=XPOL(NRCELL,NPCELL)
            X3=XPOL(NRCELL+1,NPCELL+1)
            Y1=YPOL(NRCELL+1,NPCELL)
            Y2=YPOL(NRCELL,NPCELL)
            Y3=YPOL(NRCELL+1,NPCELL+1)
          ENDIF
          IPOLG=NPCELL
        ELSEIF (LEVGEO.EQ.4) THEN
          X1=XTRIAN(NECKE(1,NCELL))
          X2=XTRIAN(NECKE(2,NCELL))
          X3=XTRIAN(NECKE(3,NCELL))
          Y1=YTRIAN(NECKE(1,NCELL))
          Y2=YTRIAN(NECKE(2,NCELL))
          Y3=YTRIAN(NECKE(3,NCELL))
        ENDIF
        Z1=0.
        Z2=0.
        Z3=0.
        CALL FPOLYT_3(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X0,Y0,ZZ)
      ELSEIF (LEVGEO.EQ.5) THEN
        X1=XTETRA(NTECK(1,NCELL))
        Y1=YTETRA(NTECK(1,NCELL))
        Z1=ZTETRA(NTECK(1,NCELL))
        X2=XTETRA(NTECK(2,NCELL))
        Y2=YTETRA(NTECK(2,NCELL))
        Z2=ZTETRA(NTECK(2,NCELL))
        X3=XTETRA(NTECK(3,NCELL))
        Y3=YTETRA(NTECK(3,NCELL))
        Z3=ZTETRA(NTECK(3,NCELL))
        X4=XTETRA(NTECK(4,NCELL))
        Y4=YTETRA(NTECK(4,NCELL))
        Z4=ZTETRA(NTECK(4,NCELL))
        CALL FPOLYT_4(X1,Y1,Z1,X2,Y2,Z2,X3,Y3,Z3,X4,Y4,Z4,X0,Y0,Z0)
      ELSEIF (LEVGEO.EQ.10) THEN
        WRITE (iunout,*) 'ERROR EXIT FROM SAMVOL. LEVGEO ',LEVGEO
        CALL EXIT_OWN(1)
      ENDIF
C
      IF (NLTRA) THEN
C  FIND Z0 FROM X01,PHI IN LOCAL TOROIDAL CELL NTCELL
        X01=X0+RMTOR
        CALL FZRTRI(X0,Z0,NTCELL,X01,PHI,NTCELL)
      ENDIF
C
      MTSURF=0
      NLSRFZ=.FALSE.
      MPSURF=0
      NLSRFY=.FALSE.
      MRSURF=0
      NLSRFX=.FALSE.
C
C  NEXT: ANALOG SPECIES INDEX DISTRIBUTION: WEISPZ(IPL)
C
      DO 630 ISPZ=1,NSPZ
        WEISPZ(ISPZ)=-1.
630   CONTINUE
C
C  NOT IN USE ANYMORE
C  CURRENTLY: ONLY SINGLE SPECIES VOLUME SOURCES POSSIBLE
C  MULTI SPECIES VOL-SOURCES HAVE TO BE TREATED BY STRATIFIED SAMPLING
C     IF (NSPEZ(ISTRA).LE.0) THEN
C       IF (NCELL.EQ.1) THEN
C         DO 640 IPL=1,NPLSI
C           IREC=0
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,1))/
C    .                  (FREC(0,  0,1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
640       CONTINUE
C       ELSE
C         DO 645 IPL=1,NPLSI
C           IFPLS=IFREC(IPLS)
C           WEISPZ(IPL)=(FREC(IFPLS,0,NCELL)-FREC(IFPLS,0,NCELL-1))/
C     .                 (FREC(0,    0,NCELL)-FREC(0,    0,NCELL-1))
C           IF (WEISPZ(IPL).LT.0) GOTO 991
645       CONTINUE
C       ENDIF
C     ENDIF
C
      CRTX=SORAD4(NVLM,ISTRA)
      CRTY=SORAD5(NVLM,ISTRA)
      CRTZ=SORAD6(NVLM,ISTRA)
      CNORM=SQRT(CRTX**2+CRTY**2+CRTZ**2)+EPS60
      CRTX=CRTX/CNORM
      CRTY=CRTY/CNORM
      CRTZ=CRTZ/CNORM
!PB
      TEWL=TEIN(NCELL)
      TIWL(1:NPLSI)=TIIN(MPLSTI(1:NPLSI),NCELL)
      DIWL(1:NPLSI)=DIIN(1:NPLSI,NCELL)
      VXWL(1:NPLSI)=VXIN(MPLSV(1:NPLSI),NCELL)
      VYWL(1:NPLSI)=VYIN(MPLSV(1:NPLSI),NCELL)
      VZWL(1:NPLSI)=VZIN(MPLSV(1:NPLSI),NCELL)
      EFWL(1:NPLSI)=0._DP
      SHWL=0._DP
C
      RETURN
C
990   CONTINUE
      WRITE (iunout,*) 'ERROR IN SAMVOL'
      CALL EXIT_OWN(1)
991   CONTINUE
      WRITE (iunout,*) 'SAMPLING ERROR IN SAMVOL'
      WRITE (iunout,*) 'NCELL,NSURF,NSBOX ',NCELL,NSURF,NSBOX
      CALL EXIT_OWN(1)
997   CONTINUE
      WRITE (iunout,*) 'SORIND (=IRRC) OUT OF RANGE IN SAMVOL'
      WRITE (iunout,*) 'IRRC,NREC ',IRRC,NREC
      CALL EXIT_OWN(1)
999   CONTINUE
      WRITE (iunout,*) 'UNWRITTEN OPTION IN SAMVOL'
      CALL EXIT_OWN(1)

C     the following ENTRY is for reinitialization of EIRENE (DMH)

      ENTRY SAMVOL_REINIT

      ISTROLD = -1

      DEALLOCATE (LPLSSR)
      DEALLOCATE (FREC)
      DEALLOCATE (VSOURC)
      DEALLOCATE (VSMXI)
      DEALLOCATE (RQ21)
      DEALLOCATE (PS21)
      IF (LEVGEO == 3) DEALLOCATE (ASIMP)
      DEALLOCATE (ISOURC)
      DEALLOCATE (ICMX)
      DEALLOCATE (IFREC)

      return

      END
