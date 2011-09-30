c jan. 08:  lower ceiling for FLX (flux dependence of Y_chem)
c           increased from 1e4 to 1e19. Otherwise at high twall (e.g. 700K)
c           the Y_chem scaled unphysically with 1/FLX.
c           The old ceiling was ok only for low TWALL, but completely wrong
c           for high TWALL. So now the new ceiling is used throughout.
c
c june 05:  new: modchm=6: Haasz/Davis 1998 formula (no flx. dep)
c           new: modchm=7: Haasz/Davis 1998 formula with flux. dep
c
c     new: user defined sputter model: modpys/modchm=9 (was: =3)
c
c           merging of flx. dep A6/A7 in roth formula for chem. sput. removed again
c           now modchm=2 is back to "flux dep option A6" (Roth) (as it already
c           was the case in all eirene version up to 2004)
c           modchm=3 is the "flux dep option A7" (Roth, 1999)
c           modchm=4 is the "flux dep option A8" (Roth, Nucl.Fus 44 (2004) L21-L)
c
c
C cvs repository jan.05  bug fix, ispz--> ispz+nsph in several places,
C                        due to photon species offset nsph
C cvs repository sept.04
C
C phys sputtering:
C still to be done:
C new formula and fitting parameters (file SPUTER_2001):
C Eckstein, 2001, ATOMIC AND PLASMA-MATERIAL INTERACTION DATA FOR FUSION,
C VOL 7B, IAEA 2001, P18 ff
C Current phys. sputter model file is SPUTER_1993, ref. see below
C
C chem sputtering:
C roth formula, revisions in aug.04 (no change in results):
C
C test integration over gaussian Etherm, as in original roth-routine 1998.
C test q=0.1 for D on C, as in original roth paper psi98 and in "warrier code"
C (this appears in 2 places !!)
C      Comp. Phys. Communic. 160 (2004) 46
C      Data file "SPUTER", stream 33, has q=0.08 for D on C
C      (seems to be more recent than phys. sputter data quoted by Roth 1998,
C       so we use this q=0.08 here)
C
C jan 05: (E.Tsitrone)
C test flux dependence, two expressions for low and high flux, A6 and A7, loc.cit.
C and connect both to provide smooth dependence as ftc. of flux,
C now same as in Warrier code and in E.Tsitrone's code
C
C june 05: connecting the flux dep. A6 and A7 removed again.
C          keep again as independent options for flx dep. sputtering yield
C          as in all previous versions (and as in Roth's original code).
C
C sept 05: use database name for opening SPUTER database
c may 06:  modifications for: photons do not sputter !
c march 07: some species flags for chemical sputtering: 
c           programming cleaned up (no change in model)

      SUBROUTINE SPUTER
C
C  GIVEN A PARTICLE: E0,VELX,VELY,VELZ,   HITS A SURFACE
C  MSURF, CRTX,CRTY,CRTZ (NORMAL UNIT VECTOR AT POINT OF
C  INCIDENCE) THEN THIS ROUTINE RETURNS DATA FOR THE
C  A) PHYSICAL SPUTTERED PARTICLE  (YIELD=YIELD1)
C  B) CHEMICALLY ERODED  PARTICLE  (YIELD=YIELD2)
C
C
C  ENTRY SPUTR0:  INITIALIZE ARRAYS
C  ENTRY SPUTR1:  CARRY OUT SPUTTERING (CALLED FROM LOCATE AND ESCAPE)
C
C  INPUT:
C    STREAM FORT.33: FILE SPUTER
C
C    WMIN:  NOT IN USE
C    FMASS: MASS NUMBER OF PROJECTILE
C    FCHAR: NUCLEAR CHARGE NUMBER OF PROJECTILE
C    FLXSP: INCIDENT FLUX (FOR FLUX DEP. CHEM.YIELD) #/CM^2/S
C    ISPZ : SPECIES INDEX OF PROJECTILE
C    MSURF : SURFACE NUMBER, needed for surface flags:
C    COSIN : CRTX*VELX+CRTY*VELY+CRTZ*VELZ
C            recys(ispz,msurf)
C            recyc(ispz,msurf)
C            itsput(1,msurf) = MODPYS
C                 MODPYS = 0:  NO PHYS. SPUTTERING (DEFAULT)
C                 MODPYS = 1:  CONSTANT YIELD: RECYCS
C                 MODPYS = 2:  ROTH/BOHDANSKY/ECKSTEIN MODEL:
C                              "REVISED BOHDANSKY FORMULAR",
C                              ECKSTEIN, W., et.al., IPP 9/117 (Garching, 1993)
C                              FITTING PARAMETERS READ FROM FILE "sputer",
C                              STREAM 33, LOC.CIT.
C                 MODPYS = 9:  USER SUPPLIED SPUTTER MODEL, CALL SPTUSR
C
C            itsput(2,msurf) = MODCHM
C                 MODCHM = 0:  NO CHEM. SPUTTERING (DEFAULT)
C                 MODCHM = 1:  CONSTANT YIELD: RECYCC
C                 MODCHM = 2:  ROTH FORMULA, PSI 1998, SAN DIEGO
C                              J.ROTH, J.NUCL.MAT 266-269 (1999) 51-57
C                              FLUX DEP. OPTION A6
C                 MODCHM = 6:  HAASZ/DAVIS FORMULA w/o flux dep.
C                 MODCHM = 7:  HAASZ/DAVIS FORMULA with flux dep.
C                 MODCHM = 9:  USER SUPPLIED SPUTTER MODEL, CALL SPTUSR
C            EWALL(MSURF)= ENWALL
C    IGASP: SPECIES INDEX FLAG FOR PHYS. SPUTTERED PARTICLE
C    IGASC: SPECIES INDEX FLAG FOR CHEM. SPUTTERED PARTICLE
C  OUTPUT:
C      YIELD1: THE PHYSICAL SPUTTERING YIELD PER INCIDENT PARTICLE
C      YIELD2: THE CHEMICAL SPUTTERING YIELD PER INCIDENT PARTICLE
C      ISPZP: SPECIES INDEX ISPZ OF PHYSICALLY SPUTTERED PARTICLE
C      ISPZC: SPECIES INDEX ISPZ OF CHEMICALLY SPUTTERED PARTICLE
C    ESPTP: ENERGY OF PHYS. SPUTTERED PARTICLE
C    VSPTP: VELOCITY OF PHYS. SPUTTERED PARTICLE
C    VXSPTP:
C    VYSPTP:   UNIT SPEED VECTOR OF PHYS. SPUTTERED PARTICLE
C    VZSPTP:
C    ESPTC: ENERGY OF CHEM. SPUTTERED PARTICLE
C    VSPTC: VELOCITY OF CHEM. SPUTTERED PARTICLE
C    VXSPTC:
C    VYSPTC:   UNIT SPEED VECTOR OF CHEM. SPUTTERED PARTICLE
C    VZSPTC:
C
      USE PRECISION
      USE PARMMOD
      USE COMUSR
      USE CADGEO
      USE CCONA
      USE CRAND
      USE CZT1
      USE CTRCEI
      USE COMPRT
      USE CLGIN
      USE CINIT
      USE CPES

      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: WMIN, FMASS, FCHAR, FLXSP
      REAL(DP), INTENT(OUT) :: YIELD1, YIELD2, ESPTC, VSPTC, ESPTP,
     .                         VSPTP, VXSPTP, VYSPTP, VZSPTP,
     .                         VXSPTC, VYSPTC, VZSPTC
      INTEGER, INTENT(IN) :: IGASC, IGASP
      INTEGER, INTENT(OUT) :: ISPZC, ISPZP
C
C  DATA FOR PHYSICAL SPUTTERING: IDENTIFY TARGET-PROJECTIL
      REAL(DP), SAVE :: ETH(28,11),Q(28,11),M2M1(28,11),ES(28)
      REAL(DP), SAVE :: D(3),EDAM(3),EDES(3)
      REAL(DP) :: RSQDV, CVRSS, RT, FTHOMP, VX, UB, GAMMA, EMAX,
     .          F1, F2, F3, SQE, QQP, ANGFAC, CAOPT, F, VY, VZ, YDES,
     .          EDESE0, EDAME0, QSE, YDAM, YSURF, ETHERM, ETHEKT,
     .          PRFCC, FLX, G2, G3, YTHERM, ERELKT, C, EREL,
     .          ENWALL, TWALL, PRFCS, ETHE0, E0ETF, SE, COSIN, QQS
      REAL(DP), EXTERNAL :: RANF_EIRENE
      INTEGER, SAVE :: ETF(28,11)
      INTEGER, SAVE :: NPROJ(11),NTARG(28),NTAMU(28)
      INTEGER, ALLOCATABLE, SAVE :: IPROJ(:),IPROJS(:),ITARG(:),
     .                              ISPZSP_DEF(:)
      INTEGER :: MS, ITYPP, IATMP, IATMC, MSS, IMOLC, ITYPC,
     .           I28, I11, IT, MODCHM, IPS, IPR,
     .           MODPYS, IIO, IPL, ILIM, ISP, IAT, IP, ISTSI, ISURF,
     .           ITA, NT, IA, NA, IFILE
      integer :: isam
      real(dp) :: wg(5),pm(5),final,flxlim
      real(dp) :: YHAASZ97M, YHAASZ97
      CHARACTER(20), SAVE :: TTARG(28)
      CHARACTER(8) :: TEXT

      DATA NPROJ/1,2,3,4,12,16,0,20,40,84,131/
      DATA NTARG/703,904,1105,1206,2713,2814,4822,5123,
     .           5224,5525,5626,5927,5928,6429,7031,7332,
     .           9140,9341,9642,10646,10847,11549,18173,
     .           18474,19578,19779,20782,23892/
      DATA TTARG/'LITHIUM             ', 'BERYLLIUM           ',
     .           'BOR                 ', 'GRAPHITE            ',
     .           'ALUMINIUM           ', 'SILIZIUM            ',
     .           'TITANIUM            ', 'VANADIUM            ',
     .           'CHROMIUM            ', 'MANGANESE           ',
     .           'IRON                ', 'COBALT              ',
     .           'NICKEL              ', 'COPPER              ',
     .           'GALLIUM             ', 'GERMANIUM           ',
     .           'ZIRCONIUM           ', 'NIOBIUM             ',
     .           'MOLYBDENUM          ', 'PALADIUM            ',
     .           'SILVER              ', 'INDIUM              ',
     .           'TANTALUM            ', 'TUNGSTEN            ',
     .           'PLATINUM            ', 'GOLD                ',
     .           'LEAD                ', 'URANIUM             '/
C  CHEMICAL EROSION DATA
      DATA D /250.,125.,83./
      DATA EDAM /15.,15.,15./
      DATA EDES /2.,2.,2./
C  PURE CARBON
      DATA EREL /1.8/
C  SI,TI,W DOPED CARBON
C     DATA EREL /1.5/
C  B DOPED CARBON
C     DATA EREL /1.2/
C
      ENTRY SPUTR0
C
C  INITIALIZE SPUTER OPTION MODPYS=2
C
      IF (MY_PE == 0) THEN
        DO IFILE=1, NDBNAMES
          IF (INDEX(DBHANDLE(IFILE),'SPUTER') /= 0) EXIT
        END DO

        IF (IFILE > NDBNAMES) THEN
          WRITE (IUNOUT,*) ' NO DATABASENAME FOR SPUTTERING DEFINED '
          WRITE (IUNOUT,*) ' CALCULATION ABANDONNED '
          CALL EXIT_OWN(1)
        END IF

        OPEN (UNIT=33,FILE=DBFNAME(IFILE))
        READ(33,*)
        READ(33,*)
        READ(33,*)
        READ(33,*)
        READ(33,*)
        READ(33,*)
        READ(33,*)
        DO I28=1,28
          READ(33,*)
          READ(33,*)
          READ(33,'(4X,E5.2)') ES(I28)
          READ(33,'(12X,11(E8.2,1X))') (M2M1(I28,I11),I11=1,11)
          READ(33,'(12X,11(I8,  1X))') (ETF(I28,I11),I11=1,11)
          READ(33,'(12X,11(E8.2,1X))') (ETH(I28,I11),I11=1,11)
          READ(33,'(12X,11(E8.2,1X))') (Q(I28,I11),I11=1,11)
        ENDDO
        CLOSE (UNIT=33)
      END IF

      if (nprs > 1) call broadsput(es,m2m1,etf,eth,q,28,11)

C  ASSIGN SPUTER DATA TO EIRENE PROJECTILE-TARGET COMBINATIONS
C  FOR THIS PARTICULAR RUN
C
C  SELF SPUTTERING: IDENTIFY TARGET ATOMIC MASS NUMBER
      DO IT=1,28
        NTAMU(IT)=NTARG(IT)/100
      ENDDO
C
      IF (.NOT.ALLOCATED(IPROJ)) THEN
        ALLOCATE (IPROJ(NSPZ))
        ALLOCATE (IPROJS(NSPZ))
        ALLOCATE (ITARG(NLIMPS))
        ALLOCATE (ISPZSP_DEF(NLIMPS))
      END IF

C  SPUTTERING BY PHOTONS
      DO 1 ISP=1,NSPH
        IPROJ(ISP)=0
        IPROJS(ISP)=0
1     CONTINUE
C
C  SPUTTERING BY ATOMS
      DO 10 ISP=NSPH+1,NSPA
        IAT=ISP-NSPH
        IPROJ(ISP)=0
        IPROJS(ISP)=0
        DO IP=1,11
          IF (NMASSA(IAT).EQ.NPROJ(IP)) IPROJ(ISP)=IP
        ENDDO
C  ANY TARGET DATA FOR SELF SPUTTERING WITH IAT?
        DO IT=1,28
          IF (NMASSA(IAT).EQ.NTAMU(IT)) IPROJS(ISP)=IT
        ENDDO
10    CONTINUE
C
C  SPUTTERING BY MOLECULES
      DO 20 ISP=NSPA+1,NSPAM
        IPROJ(ISP)=0
        IPROJS(ISP)=0
20    CONTINUE
C
C  SPUTTERING BY TEST IONS
      DO 30 ISP=NSPAM+1,NSPAMI
        IIO=ISP-NSPAM
        IPROJ(ISP)=0
        IPROJS(ISP)=0
C  TEST: MOLECULAR OR ATOMIC ION?
        IF (NPRT(ISP).NE.1) GOTO 30
        DO IP=1,11
          IF (NMASSI(IIO).EQ.NPROJ(IP)) IPROJ(ISP)=IP
        ENDDO
C  ANY TARGET DATA FOR SELF SPUTTERING WITH IIO?
        DO IT=1,28
          IF (NMASSI(IIO).EQ.NTAMU(IT)) IPROJS(ISP)=IT
        ENDDO
30    CONTINUE
C
C  SPUTTERING BY BULK IONS
      DO 40 ISP=NSPAMI+1,NSPTOT
        IPL=ISP-NSPAMI
        IPROJ(ISP)=0
        IPROJS(ISP)=0
C  TEST: MOLECULAR OR ATOMIC ION?
        IF (NPRT(ISP).NE.1) GOTO 40
        DO IP=1,11
          IF (NMASSP(IPL).EQ.NPROJ(IP))
     .    IPROJ(ISP)=IP
        ENDDO
C  ANY TARGET DATA FOR SELF SPUTTERING WITH IPL?
        DO IT=1,28
          IF (NMASSP(IPL).EQ.NTAMU(IT))
     .    IPROJS(ISP)=IT
        ENDDO
40    CONTINUE

      ITARG=0
      ISPZSP_DEF=0
      DO ILIM=1,NLIMI
        NT=100.*ZNML(ILIM)+ZNCL(ILIM)+EPS10
        ITARG(ILIM)=0
        DO IT=1,28
          IF (NT.EQ.NTARG(IT)) ITARG(ILIM)=IT
        ENDDO
        DO IA=1,NATMI
          NA=100*NMASSA(IA)+NCHARA(IA)
          IF (NT.EQ.NA) ISPZSP_DEF(ILIM)=IA
        ENDDO
      ENDDO
      DO ILIM=NLIM+1,NLIM+NSTSI
        NT=100.*ZNML(ILIM)+ZNCL(ILIM)+EPS10
        ITARG(ILIM)=0
        DO IT=1,28
          IF (NT.EQ.NTARG(IT)) ITARG(ILIM)=IT
        ENDDO
        DO IA=1,NATMI
          NA=100*NMASSA(IA)+NCHARA(IA)
          IF (NT.EQ.NA) ISPZSP_DEF(ILIM)=IA
        ENDDO
      ENDDO
C
      IF (TRCREF) THEN
        WRITE (iunout,*) 
     .    'PRINTOUT FROM SUBR. SPUTER, AFTER INITIALISATION'
C
        WRITE (iunout,*)
        WRITE (iunout,*) 'ISPZ,IPROJ(ISPZ),IPROJS(ISPZ)'
        DO ISP=1,NSPTOT
          WRITE (iunout,*) TEXTS(ISP),IPROJ(ISP),IPROJS(ISP)
        ENDDO
        WRITE (iunout,*)
        WRITE (iunout,*) 
     .    'SURFACE NUMBER, TARGET MATERIAL, SPUTTERED ATOM'
        DO ILIM=1,NLIMI
          IF (ILIIN(ILIM).LE.0) THEN
            WRITE (iunout,*) ILIM, ' TRANSPARENT SURFACE '
          ELSEIF (ILIIN(ILIM).GE.3) THEN
            WRITE (iunout,*) ILIM, ' PERIODICITY- OR MIRROR SURFACE '
          ELSEIF (IGJUM0(ILIM).EQ.1) THEN
            WRITE(iunout,*) ILIM, ' SURFACE OUT '
          ELSEIF (ITARG(ILIM).GT.0.AND.ITARG(ILIM).LE.28) THEN
            WRITE(iunout,*) ILIM,' ',TTARG(ITARG(ILIM)),ISPZSP_DEF(ILIM)
          ELSE
            WRITE(iunout,*) ILIM, ' TARGET MATERIAL NOT IDENTIFIED '
          ENDIF
        ENDDO
        DO ISTSI=1,NSTSI
          ISURF=NLIM+ISTSI
          IF (ILIIN(ISURF).LE.0) THEN
            WRITE(iunout,*) -ISTSI,' TRANSPARENT SURFACE '
          ELSEIF (ILIIN(ISURF).GE.3) THEN
            WRITE(iunout,*) -ISTSI,' PERIODICITY- OR MIRROR SURFACE '
          ELSEIF (ITARG(ISURF).GT.0.AND.ITARG(ISURF).LE.28) THEN
            WRITE(iunout,*) -ISTSI,' ',TTARG(ITARG(ISURF)),
     .                      ISPZSP_DEF(ISURF)
          ELSE
            WRITE(iunout,*) -ISTSI,' TARGET MATERIAL NOT IDENTIFIED '
          ENDIF
        ENDDO
        WRITE (iunout,*)
C       DO I28=1,28
C         DO I11=1,11
C           WRITE (iunout,*) NTARG(I28),NPROJ(I11)
C           WRITE (iunout,*) 'ES   ',ES(I28)
C           WRITE (iunout,*) 'M2M1 ',M2M1(I28,I11)
C           WRITE (iunout,*) 'ETF  ',ETF(I28,I11)
C           WRITE (iunout,*) 'ETH  ',ETH(I28,I11)
C           WRITE (iunout,*) 'Q    ',Q(I28,I11)
C           WRITE (iunout,*)
C         ENDDO
C       ENDDO
      ENDIF
C
      CALL SP0USR
C
      RETURN
C
      ENTRY SPUTR1(WMIN,FMASS,FCHAR,FLXSP,
     .             IGASP,
     .             YIELD1,
     .             ISPZP,ESPTP,VSPTP,VXSPTP,VYSPTP,VZSPTP,
     .             IGASC,
     .             YIELD2,
     .             ISPZC,ESPTC,VSPTC,VXSPTC,VYSPTC,VZSPTC)

C  TENTATIVELY ASSUME: NO SPUTTERED PARTICLES
      YIELD1=0.D0
      YIELD2=0.D0
      ISPZP=0
      ISPZC=0

C  CURRENTLY: PHOTONS DO NOT SPUTTER
      IF (ISPZ.LE.NSPH) RETURN
C
C
C  SURFACE NUMBER :  MSURF  (MSURF=0: DEFAULT MODEL)
C  SPECIES INDEX  :  ISPZ   (INCIDENT SPECIES)
C  IDENTIFY PROJECTILE AND TARGET
C
      ITA=ITARG(MSURF)
      IPS=IPROJS(ISPZ)
      IF (IPS.GT.0.AND.IPS.EQ.ITA) THEN
        IPR=7
      ELSE
        IPR=IPROJ(ISPZ)
      ENDIF
C
      MODPYS=ISPUT(1,MSURF)
      MODCHM=ISPUT(2,MSURF)
      ENWALL=EWALL(MSURF)
      IF (ENWALL.GT.0.D0) THEN
        TWALL=ENWALL*0.66667
      ELSE
        TWALL=-ENWALL
      ENDIF
      SE=0.D0
C
      COSIN=CRTX*VELX+CRTY*VELY+CRTZ*VELZ
      IF (COSIN.LE.0.D0) GOTO 999
C
C  NO PHYSICAL SPUTTERING?
C
      IF (MODPYS.LE.0) GOTO 5000
C
C  FIRST: PHYSICAL SPUTTERING AND SUBLIMATION
C
      IF (MODPYS.EQ.1) THEN
C
C  CONSTANT SPUTTER YIELD: RECYCS
        YIELD1=RECYCS(ISPZ,MSURF)
C  NO SUBLIMATION YIELD IN THIS MODEL:
        QQS=0.
C
      ELSEIF (MODPYS.EQ.2) THEN
C
C   ECKSTEIN/ROTH/BOHDANSKY/MODEL: IPP 9/82, FEB. 1993
C
        IF (IPR.GT.0.AND.ITA.GT.0) THEN
          PRFCS=RECYCS(ISPZ,MSURF)
C  NO SPUTTERING BELOW THRESHOLD
          IF (E0.LE.ETH(ITA,IPR).OR.PRFCS.LE.0.D0) GOTO 5000
          ETHE0=ETH(ITA,IPR)/E0
          E0ETF=E0/ETF(ITA,IPR)
          SQE=SQRT(E0ETF)
C  YIELD FACTOR FOR PHYS. SPUTTERING
          QQP=Q(ITA,IPR)
C  YIELD FACTOR FOR SUBLIMATION
          QQS=0.
c         IF (IFLAG.EQ.2) THEN
c           QQS=54.*FMASS**1.2*EXP(-0.78/TWALL)
c         ENDIF
          F1=(QQP+QQS)*(1.-ETHE0**0.666667)*(1.-ETHE0)*(1.-ETHE0)
c  replace Thomas-Fermi potential by Kr-C potential
c  Thomas-Fermi potential
c         F2=3.441*SQE*LOG(E0ETF+2.718)
c         F3=1.+6.355*SQE+E0ETF*(6.882*SQE-1.708)
c  Kr-C potential
          F2=0.5*LOG(1.+1.2288*E0ETF)
          F3=E0ETF+0.1728*SQE+0.008*E0ETF**0.1504
C
          SE=F2/F3
          YIELD1=F1*SE
          YIELD1=MAX(0._DP,YIELD1)*PRFCS
C INCIDENT ANGULAR DEPENDENCE OF YIELD: YAMAMURA FIT: LOC.CIT.,P 10
C MAXIMUM (ABOUT: 3.36) AT COSIN=0.26.
C ANGFAC -->0. FOR COSIN -->0.
C ANGFAC -->1. FOR COSIN -->1.
C         AOPT=75.
C         CAOPT=COS(AOPT*PIA/180.D0)
          CAOPT=0.26
          F=2.
          ANGFAC=COSIN**(-F)*EXP(F*(1.-1./COSIN)*CAOPT)
          YIELD1=YIELD1*ANGFAC
        ENDIF
C
      ELSEIF (MODPYS.EQ.9) THEN
C  USER SUPPLIED SPUTER MODEL
        CALL SP1USR
C
      ENDIF
C
C  FIND TYPE AND SPECIES INDEX OF SPUTTERED PARTICLE
C  ONLY ATOMS
C
      IF (YIELD1.GT.0.D0) THEN
        IF (IGASP.GT.0.AND.IGASP.LE.NATMI) THEN
          ITYPP=1
          IATMP=IGASP
          ISPZP=IATMP+NSPH
        ELSEIF (IGASP.GT.NATMI) THEN
          WRITE (iunout,*) 'ERROR: ISPZP OUT OF RANGE IN SPUTER ?'
          ISPZP=0
          GOTO 5000
        ELSEIF (IGASP.LT.0) THEN
C  DEFAULT: TRY TO FIND SPECIES INDEX OF SPUTTERED PARTICLE AUTOMATICALLY
          ISPZP=ISPZSP_DEF(MSURF)
          IF (ISPZP.LE.0) THEN
            MS=MSURF
            IF (MSURF.GT.NLIM) MS=-(MSURF-NLIM)
            WRITE (iunout,*) 'WARNING FROM SPUTER '
            WRITE (iunout,*) 'SPECIES INDEX, PHYS. SPUTER ? MSURF ',MS
            ISPZP=0
            GOTO 5000
          ENDIF
          ITYPP=1
          IATMP=ISPZP
          ISPZP=IATMP+NSPH
        ELSE
C  IGASP=0 :
C  RETURN ONLY SPUTTER YIELD, NOT THE COORDINATES OF SPUTTERED PARTICLE
          ISPZP=0
          GOTO 5000
        ENDIF
      ELSE
        GOTO 5000
      ENDIF
C
C  SAMPLE ENERGY FROM THOMPSON DISTRIBUTION IN CASE OF PHYS. SPUTTERING
C  USE ENERGY AND ANGULAR DISTRIBUTION SAME AS FOR CHEMICAL SPUTTERING
C  IN CASE OF RADIATION ENHANCED SUBLIMATION
C
      IF (QQS.GT.0.D0) THEN
        RT=QQS/(QQS+QQP)
        IF (RANF_EIRENE( ).LE.RT) THEN
          RSQDV=RSQDVA(IATMP)
          CVRSS=CVRSSA(IATMP)
          ESPTP=ENWALL
          GOTO 1000
        ENDIF
      ENDIF
C
C  COMPUTE SURFACE BINDING ENERGY U0(IFLAG) FROM ETH(IFLAG) FOR THOMPSON
C  DISTRIBUTION
C
      GAMMA=4.*FMASS*ZNML(MSURF)/(FMASS+ZNML(MSURF))**2
      EMAX=GAMMA*E0
      UB=ES(ITA)
      ESPTP=FTHOMP(UB,EMAX)
      RSQDV=RSQDVA(IATMP)
      GOTO 1000
C
C  ANGULAR DISTRIBUTIONS
C
1000  CONTINUE
C
      IF (ESPTP.LE.0.D0) GOTO 1100
C
      VSPTP=RSQDV*SQRT(ESPTP)
C
C  SAMPLE SPEED VECTOR FROM COSINE
C
      IF (INIV4.LE.0) CALL FCOSIN
      VX=FC1(INIV4)
      VY=FC2(INIV4)
      VZ=FC3(INIV4)
      INIV4=INIV4-1
      CALL ROTATF (VXSPTP,VYSPTP,VZSPTP,VX,VY,VZ,CRTX,CRTY,CRTZ)
      GOTO 5000
C
1100  CONTINUE
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
      CALL VELOCS (TWALL,0._DP,0._DP,0._DP,0._DP,0._DP,RSQDV,CVRSS,
     .             -CRTX,-CRTY,-CRTZ,
     .             ESPTP,VXSPTP,VYSPTP,VZSPTP,VSPTP)
C
2000  CONTINUE
C
C   PHYSICAL SPUTTERING DONE
C
C   CHEMICAL SPUTTERING, REEMITTED PARTICLES ARE COSINE DISTRIBUTED AND
C   THERMAL
C
5000  CONTINUE
C
C  NO CHEMICAL SPUTTERING?
C
      IF (MODCHM.LE.0) RETURN
C
      SELECT CASE (MODCHM)

      CASE(1)
C  CONSTANT SPUTTER YIELD: RECYCS
C  IS INCIDENT PARTICLE "HYDROGENIC" AND "ATOMIC"?
C  IS TARGET SURFACE CARBON?
C
        IF (IPR.GT.0.AND.IPR.LE.3.AND.ITA.EQ.4) THEN
          YIELD2=RECYCC(ISPZ,MSURF)
        ELSE
C  NO CHEM. SPUTTERING DATA FOR THIS TARGET-PROJECTILE COMBINATION
          GOTO 20000
        ENDIF
C
      CASE(2,3,4)
C
C   ROTH/PACHER MODEL: PSI 1998, SAN DIEGO (J.NUCL.MAT)
C
        IF (IPR.GT.0.AND.IPR.LE.3.AND.ITA.EQ.4) THEN
C  FLXSP IS IN #/CM^2/S. CONVERT TO #/M^2/S
          FLX=FLXSP*1.D4
C  CEILING OF FLX: 1E19 #/S/M**2. FOR LOWER FLX AND AT HIGH TWALL
C                                 THE FORMULA BECOMES UNPHYSICAL (PROTO 1/FLX)
          FLX=MAX(1.E19_DP,FLX)
C
          PRFCC=RECYCC(ISPZ,MSURF)
C  TO BE WRITTEN
C    ROTH, J.NUCL.MAT 99: SAMPLE ETHERM FROM GAUSSIAN, MEAN: 1.7, SIGMA 0.3
C  THIS PRESENT IMPLEMENTATION:
C    ROTH, NUCL.FUS. 96 (P1647): USE ONLY THE MEAN OF ETHERM
cdr june 2004:
cdr test integration expression from original roth-code
cdr       pm(1)=1.865
cdr       pm(2)=1.7
cdr       pm(3)=1.535
cdr       pm(4)=1.38
cdr       pm(5)=1.26
cdr       wg(1)=1./4.5
cdr       wg(2)=1./4.5
cdr       wg(3)=1./4.5
cdr       wg(4)=1./9.
cdr       wg(5)=1./9.
cdr       final=0.
cdr       do 6000 isam=1,5
cdr
cdr no, take only mean value pm=1.7, same as e.g. in
cdr all older eirene versions and warrier-code
          pm(2)=1.7
          wg(2)=1.
          isam=2
          final=0.
cdr
          ETHERM=pm(isam)
          ETHEKT=EXP(-ETHERM/TWALL)
C
          ERELKT=EXP(-EREL/TWALL)

cdr  continuous merging of option A6 and A7, as in Warrier code. Out!
cdr        FLXLIM=1.D30*EXP(-1.4/TWALL)
cdr        IF (FLX.LE.FLXLIM) THEN ! A6,  else: A7

C  EXPRESSION A.6 FOR C, WEAK-FLUX DEPENDENCE
            IF (MODCHM.EQ.2) C=1._DP/(1._DP+3.E7_DP*EXP(-1.4_DP/TWALL))
C  EXPRESSION A.7 FOR C, STRONG-FLUX DEPENDENCE
            IF (MODCHM.EQ.3) C=1._DP/(1._DP+3.E-23_DP*FLX)
C  EXPRESSION A.8 FOR C, new FLUX DEPENDENCE, roth, itpa 2003 st. petersburg
            IF (MODCHM.EQ.4) C=1._DP/(1._DP+(1.67E-22_DP*FLX)**0.54)
cdr       ENDIF
C
          G2=2.D-32*FLX+ETHEKT
          G3=2.D-32*FLX+(1._DP+2.D+29/FLX*ERELKT)*ETHEKT
          YTHERM=C/G3*0.033_DP*ETHEKT
          IF (SE.GT.0.) THEN
cdr in warrier code: qqp=0.1, for D on C
cdr         qqp=0.1
cdr
            QSE=QQP*SE
          ELSE
c  Kr-C potential
            E0ETF=E0/ETF(ITA,IPR)
            SQE=SQRT(E0ETF)
            QQP=Q(ITA,IPR)
            F2=0.5_DP*LOG(1._DP+1.2288_DP*E0ETF)
            F3=E0ETF+0.1728_DP*SQE+0.008_DP*E0ETF**0.1504_DP
            SE=F2/F3
cdr in warrier code: qqp=0.1, for D on C
cdr         qqp=0.1
cdr
            QSE=QQP*SE
          ENDIF
c
          YDAM=0.
          IF (E0.GE.EDAM(IPR)) THEN
            EDAME0=EDAM(IPR)/E0
            YDAM=QSE*(1._DP-EDAME0**0.666667_DP)*
     .               (1._DP-EDAME0)*(1._DP-EDAME0)
           ENDIF
C
          YDES=0.
          IF (E0.GE.EDES(IPR)) THEN
            EDESE0=EDES(IPR)/E0
            YDES=QSE*(1._DP-EDESE0**0.666667_DP)*
     .               (1._DP-EDESE0)*(1._DP-EDESE0)
          ENDIF
C
          YSURF=0._DP
          IF (E0 < 1000._DP) YSURF=C*G2/G3*YDES/
     .                             (1._DP+EXP((E0-65._DP)/40._DP))
          YIELD2=YTHERM*(1._DP+D(IPR)*YDAM)+YSURF
          YIELD2=MAX(0._DP,YIELD2)*PRFCC
cdr
          final=final+yield2*wg(isam)
cdr
6000      continue
          yield2=final
cdr
        ELSE
C  NO CHEM. SPUTTERING DATA FOR THIS TARGET-PROJECTILE COMBINATION
          GOTO 20000
        ENDIF
C
      CASE(6)
C  Haasz-Davis formula, 1998
         yield2=yhaasz97m(e0,twall)
      CASE(7)
C  Haasz-Davis formula, 1998, with flx. dep from Roth, Nucl.Fus 2004
         C=1._DP/(1._DP+(1.67E-22_DP*FLX)**0.54)
         yield2=C * yhaasz97m(e0,twall)
      CASE(9)
C  USER SUPPLIED SPUTER MODEL
        CALL SP1USR
      CASE DEFAULT
        write (iunout,*) 'error in sputer.f. modchm ?? ',modchm
        call exit_own(1)
      END SELECT
C
C  FIND TYPE AND SPECIES OF CHEM. SPUTTERED MOLECULE
C  ATOMS OR MOLECULES
C
      IF (YIELD2.GT.0.D0) THEN
        IF (IGASC.GT.0.AND.IGASC.LE.NATMI+NMOLI) THEN
          IF (ISPZC.GT.NATMI) THEN
            ITYPC=2
            IMOLC=IGASC-NATMI
            ISPZC=IMOLC+NSPA
          ELSE
            ITYPC=1
            IATMC=IGASC
            ISPZC=IATMC+NSPH
          ENDIF
        ELSEIF (IGASC.GT.NATMI+NMOLI) THEN
          WRITE (iunout,*) 'ERROR: ISPZC OUT OF RANGE IN SPUTER ?'
          ISPZC=0
          GOTO 20000
        ELSEIF (IGASC.LT.0) THEN
C  DEFAULT: TRY TO FIND SPECIES INDEX OF SPUTTERED PARTICLE AUTOMATICALLY
          ISPZC=ISPZSP_DEF(MSURF)
          IF (ISPZC.LE.0) THEN
            MS=MSURF
            IF (MSURF.GT.NLIM) MS=-(MSURF-NLIM)
            WRITE (iunout,*) 'WARNING FROM SPUTER '
            WRITE (iunout,*) 'SPECIES INDEX, CHEM. SPUTER ? MSURF ',MS
            ISPZC=0
            GOTO 20000
          ENDIF
          ITYPC=1
          IATMC=ISPZC
          ISPZC=IATMC+NSPH
        ELSE
C  IGASC=0 :
C  RETURN ONLY SPUTTER YIELD, NOT THE COORDINATES OF SPUTTERED PARTICLE
          ISPZC=0
          GOTO 20000
        ENDIF
      ELSE
        GOTO 20000
      ENDIF
C
      ESPTC=ENWALL
      IF (ITYPC.EQ.1) THEN
        RSQDV=RSQDVA(IATMC)
        CVRSS=CVRSSA(IATMC)
      ELSE
        RSQDV=RSQDVM(IMOLC)
        CVRSS=CVRSSM(IMOLC)
      ENDIF
C
C  ANGULAR DISTRIBUTIONS
C
10000 CONTINUE
C
      IF (ESPTC.LE.0.D0) GOTO 11000
C
      VSPTC=RSQDV*SQRT(ESPTC)
C
C  SAMPLE SPEED VECTOR FROM COSINE
C
      IF (INIV4.LE.0) CALL FCOSIN
      VX=FC1(INIV4)
      VY=FC2(INIV4)
      VZ=FC3(INIV4)
      INIV4=INIV4-1
      CALL ROTATF (VXSPTC,VYSPTC,VZSPTC,VX,VY,VZ,CRTX,CRTY,CRTZ)
      RETURN
C
11000 CONTINUE
C  SAMPLE FROM MAXWELLIAN FLUX AROUND INNER (!) NORMAL AT TEMP. TW (EV)
      CALL VELOCS (TWALL,0._DP,0._DP,0._DP,0._DP,0._DP,RSQDV,CVRSS,
     .             -CRTX,-CRTY,-CRTZ,
     .             ESPTC,VXSPTC,VYSPTC,VZSPTC,VSPTC)
C
20000 RETURN
C
999   CONTINUE
      WRITE (iunout,*) 'ERROR IN SUBR. SPUTER '
      MSS=MSURF
      IF (MSS.GT.NLIM) MSS=-(MSURF-NLIM)
      WRITE (iunout,*) 'MSURF = ',MSS
      WRITE (iunout,*) 'COSIN.LT.0. ', COSIN
      WRITE (iunout,*) 'DO NOT SPUTTER FOR PARTICLE NO. NPANU= ',NPANU
      RETURN
      END



      FUNCTION YHAASZ97(E0,TEMP_eV)
c  this is the haasz/davis chemical sputtering yield (1997).
c  program provided via Toronto group (p.c. stangeby/d. elder)
c  d.reiter: single--> double prec.
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  E0  (eV)       -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *  D.REITER: MODIFIED: TEMP_EV (EV)                                 *
C  *                                                                   *
C  *********************************************************************
C
      use precision
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: E0, TEMP_EV
      REAL(DP) ::TEMP
      REAL(DP) ::FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >           FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >           FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4)
      REAL(DP) :: POLY_C(4),YFIT,FITE0
      REAL(DP) :: YHAASZ97
      INTEGER I
C
C     Poly. fit c. /       a0,      a1,      a2,      a3
C
      DATA FITC300 / -0.03882, 0.07432,-0.03470, 0.00486/
      DATA FITC350 / -0.05185, 0.10126,-0.05065, 0.00797/
      DATA FITC400 / -0.06089, 0.12186,-0.06240, 0.01017/
      DATA FITC450 / -0.08065, 0.16884,-0.09224, 0.01625/
      DATA FITC500 / -0.08872, 0.19424,-0.10858, 0.01988/
      DATA FITC550 / -0.08728, 0.20002,-0.11420, 0.02230/
      DATA FITC600 / -0.05106, 0.13146,-0.07514, 0.01706/
      DATA FITC650 /  0.07373,-0.13263, 0.09571,-0.01672/
      DATA FITC700 /  0.02722,-0.03599, 0.02064, 0.00282/
      DATA FITC750 /  0.09052,-0.18253, 0.12362,-0.02109/
      DATA FITC800 /  0.02604,-0.05480, 0.04025,-0.00484/
      DATA FITC850 /  0.03478,-0.08537, 0.06883,-0.01404/
      DATA FITC900 /  0.02173,-0.06399, 0.05862,-0.01380/
      DATA FITC950 / -0.00086,-0.01858, 0.02897,-0.00829/
      DATA FITC1000/ -0.01551, 0.01359, 0.00600,-0.00353/
C
C in calling program (eirene), temp_EV is in eV
c convert to K
      TEMP=TEMP_EV*11604.
c
C Find right polynomial fit coefficients for a given temperature
C
c
      IF      (TEMP.LE.300.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
         ENDDO
      ELSE IF (TEMP.GT.300.0 .AND. TEMP.LE.350.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
         ENDDO
      ELSE IF (TEMP.GT.350.0 .AND. TEMP.LE.400.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
         ENDDO
      ELSE IF (TEMP.GT.400.0 .AND. TEMP.LE.450.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
         ENDDO
      ELSE IF (TEMP.GT.450.0 .AND. TEMP.LE.500.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
         ENDDO
      ELSE IF (TEMP.GT.500.0 .AND. TEMP.LE.550.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
         ENDDO
      ELSE IF (TEMP.GT.550.0 .AND. TEMP.LE.600.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
         ENDDO
      ELSE IF (TEMP.GT.600.0 .AND. TEMP.LE.650.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
         ENDDO
      ELSE IF (TEMP.GT.650.0 .AND. TEMP.LE.700.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
         ENDDO
      ELSE IF (TEMP.GT.700.0 .AND. TEMP.LE.750.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
         ENDDO
      ELSE IF (TEMP.GT.750.0 .AND. TEMP.LE.800.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
         ENDDO
      ELSE IF (TEMP.GT.800.0 .AND. TEMP.LE.850.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
         ENDDO
      ELSE IF (TEMP.GT.850.0 .AND. TEMP.LE.900.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
         ENDDO
      ELSE IF (TEMP.GT.900.0 .AND. TEMP.LE.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
         ENDDO
      ELSE IF (TEMP.GT.950.0) THEN
         DO I = 1,4
           POLY_C(I) = FITC1000(I)
         ENDDO
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0)  THEN
         FITE0 = 10.
      ELSE IF (E0.GT.200.0) THEN
         FITE0 = 200.
      ELSE
         FITE0 = E0
      ENDIF
C
      YFIT = 0.0
      DO I = 1,4
        YFIT = YFIT + POLY_C(I)*LOG10(FITE0)**(I-1)
      ENDDO

      YHAASZ97 = YFIT

CW    WRITE(iunout,*) 'YHAASZ97 = ',YHAASZ97

      RETURN
      END
c
c
c
      FUNCTION YHAASZ97M(E0,TEMP_eV)
c  this is the haasz/davis chemical sputtering yield (1998).
c  program provided via Toronto group (p.c. stangeby/d. elder)
c  d.reiter: single--> double prec.
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J.Davis and parameterized by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV   *
C  *  E0  (eV)       -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *  D.REITER: MODIFIED: TEMP_EV (EV)                                 *
C  *                                                                   *
C  *********************************************************************
C
      use precision
      IMPLICIT NONE

      REAL(DP), INTENT(IN) :: E0, TEMP_EV
      real(dp) :: TEMP
      real(dp) :: YHAASZ97M, YDAVIS98, YHAASZ97
      real(dp) :: m1,m2,m3,reducf,FRAC

      DATA m1/602.39/, m2/202.24/, m3/43.561/, reducf/0.2/
C
C in calling program (eirene), TEMP_EV is in eV
c convert to K
      TEMP=TEMP_EV*11604.

      IF (E0 .GE. 10) THEN
         YHAASZ97M = YHAASZ97(E0,TEMP_eV)
      ELSEIF (E0 .LT. 10. .AND. E0 .GE. 5.) THEN
         FRAC = (E0-5.)/5.
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3)
         YHAASZ97M = FRAC*YHAASZ97(E0,TEMP_eV)+ (1.-FRAC)*YDAVIS98
      ELSEIF (E0 .LT. 5.) THEN
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1)**2 + m3)
         YHAASZ97M = YDAVIS98
      ENDIF

      RETURN
      END
